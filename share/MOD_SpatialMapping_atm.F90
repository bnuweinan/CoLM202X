#include <define.h>
#ifdef CCPL

MODULE MOD_SpatialMapping_atm

!--------------------------------------------------------------------------------
! !DESCRIPTION:
!
!    Spatial Mapping module for atmospheric forcing data with irregular grids or atmosphere model grids.
!
!  Created by Nan Wei, Feb 2025, template from MOD_SpatialMapping.F90
!--------------------------------------------------------------------------------

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   ! ------
   type :: atm_grid_list_type
      integer :: ng
      integer, allocatable :: elmid(:)
   END type atm_grid_list_type

   type :: patm_grid_list_type
      integer :: ng
      integer, allocatable :: ilat(:)
      integer, allocatable :: ilon(:)
      integer, allocatable :: elmid(:)
   END type patm_grid_list_type

   type :: spatial_mapping_type

      integer, allocatable :: elmid(:)     ! atmospheric grid element index on each IO process

      type(atm_grid_list_type), allocatable :: glist (:)

      integer :: npset
      integer, allocatable :: npart(:)
      type(pointer_int32_2d), allocatable :: address (:)

      logical  :: has_missing_value = .false.
      real(r8) :: missing_value     = spval

      type(pointer_real8_1d), allocatable :: areapart(:) ! intersection area
      real(r8), allocatable               :: areapset(:)
      real(r8), allocatable               :: areagrid(:)

   CONTAINS

      procedure, PUBLIC :: build_arealweighted => spatial_mapping_build_arealweighted

      ! 1) from pixelset to grid
      procedure, PUBLIC  :: pset2grid       => spatial_mapping_pset2grid
      procedure, PUBLIC  :: pset2grid_max   => spatial_mapping_pset2grid_max

      procedure, PUBLIC  :: get_sumarea  => spatial_mapping_get_sumarea

      ! 2) from grid to pixelset
      procedure, PUBLIC  :: grid2pset    => spatial_mapping_grid2pset

      ! 3) between grid and intersections
      procedure, PUBLIC  :: grid2part => spatial_mapping_grid2part
      procedure, PUBLIC  :: part2grid => spatial_mapping_part2grid
      procedure, PUBLIC  :: normalize => spatial_mapping_normalize

      ! 4) intersections to pixelset
      procedure, PUBLIC  :: part2pset => spatial_mapping_part2pset

      procedure, PUBLIC  :: allocate_part   => spatial_mapping_allocate_part
      procedure, PUBLIC  :: deallocate_part => spatial_mapping_deallocate_part
      procedure, PUBLIC  :: forc_free_mem   => forc_free_mem_spatial_mapping

      final :: spatial_mapping_free_mem

   END type spatial_mapping_type

!-----------------------
CONTAINS

   !------------------------------------------
   SUBROUTINE spatial_mapping_build_arealweighted (this, fgrid, pixelset, elmid_atm, mesh_atm_pio, numelm_atm, mesh_atm, filter)

   USE MOD_Precision
   USE MOD_Namelist
   USE MOD_Block
   USE MOD_Pixel
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_Mesh
   USE MOD_Mesh_atm
   USE MOD_Utils
   USE MOD_UserDefFun
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(grid_type),           intent(in) :: fgrid
   type(pixelset_type),       intent(in) :: pixelset
   type(block_data_int32_2d), intent(in) :: elmid_atm
   integer, intent(in) :: mesh_atm_pio(:)
   integer, intent(in) :: numelm_atm
   type (irregular_elm_type_atm), intent(inout) :: mesh_atm (:)
   logical, intent(in) :: filter(:)

   ! Local variables
   type(pointer_real8_1d), allocatable :: afrac(:)
   type(atm_grid_list_type),   allocatable :: gfrom(:)
   type(pointer_int32_1d), allocatable :: list_lat(:)
   integer,  allocatable :: ng_lat(:)
   type(pointer_real8_1d), allocatable :: pafrac(:)
   type(patm_grid_list_type), allocatable :: gplist(:)
   type(patm_grid_list_type), allocatable :: gpfrom(:)
   type(pointer_int32_2d), allocatable :: paddress(:)
   integer,  allocatable :: ys(:), yn(:), xw(:), xe(:)
   integer,  allocatable :: xlist(:), ylist(:), elist(:)
   integer,  allocatable :: ipt(:)
   logical,  allocatable :: msk(:)

   integer  :: ie, iset, iio
   integer  :: ng, ig, ng_all, iloc
   integer  :: npxl, ipxl, ilat, ilon
   integer  :: iworker, iproc, idest, isrc, nrecv
   integer  :: rmesg(2), smesg(2)
   integer  :: iy, ix, xblk, yblk, xloc, yloc
   integer  :: ipxstt, ipxend
   real(r8) :: lat_s, lat_n, lon_w, lon_e, area
   logical  :: skip, is_new


#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (p_is_io .and. numelm_atm > 0) THEN
         allocate (this%elmid (numelm_atm))
         DO ie = 1, numelm_atm
            this%elmid(ie) = mesh_atm(ie)%indx
         ENDDO
      ENDIF

IF(trim(DEF_file_mesh_atm) /= trim(DEF_file_mesh)) then

      IF (p_is_worker) THEN

         this%npset = pixelset%nset

         allocate (pafrac (pixelset%nset))
         allocate (gpfrom (pixelset%nset))

         allocate (ys (pixel%nlat))
         allocate (yn (pixel%nlat))
         allocate (xw (pixel%nlon))
         allocate (xe (pixel%nlon))

         DO ilat = 1, pixel%nlat
            ys(ilat) = find_nearest_south (pixel%lat_s(ilat), fgrid%nlat, fgrid%lat_s)
            yn(ilat) = find_nearest_north (pixel%lat_n(ilat), fgrid%nlat, fgrid%lat_n)
         ENDDO

         DO ilon = 1, pixel%nlon
            xw(ilon) = find_nearest_west (pixel%lon_w(ilon), fgrid%nlon, fgrid%lon_w)
            xe(ilon) = find_nearest_east (pixel%lon_e(ilon), fgrid%nlon, fgrid%lon_e)
         ENDDO

         allocate (list_lat (fgrid%nlat))
         DO iy = 1, fgrid%nlat
            allocate (list_lat(iy)%val (1000))
         ENDDO

         allocate (ng_lat (fgrid%nlat)); ng_lat(:) = 0

         DO iset = 1, pixelset%nset

            ie = pixelset%ielm(iset)
            npxl = pixelset%ipxend(iset) - pixelset%ipxstt(iset) + 1

            ipxstt = pixelset%ipxstt(iset)
            ipxend = pixelset%ipxend(iset)

            ! deal with 2m WMO patch
            IF (ipxstt==-1 .and. ipxend==-1) THEN
               ipxstt = 1
               ipxend = mesh(ie)%npxl
               npxl   = mesh(ie)%npxl
            ENDIF

            allocate (pafrac(iset)%val (npxl))
            allocate (gpfrom(iset)%ilat(npxl))
            allocate (gpfrom(iset)%ilon(npxl))

            gpfrom(iset)%ng = 0

            DO ipxl = ipxstt, ipxend

               ilat = mesh(ie)%ilat(ipxl)
               ilon = mesh(ie)%ilon(ipxl)

               DO iy = ys(ilat), yn(ilat), fgrid%yinc

                  lat_s = max(fgrid%lat_s(iy), pixel%lat_s(ilat))
                  lat_n = min(fgrid%lat_n(iy), pixel%lat_n(ilat))

                  IF ((lat_n-lat_s) < 1.0e-6_r8) THEN
                     CYCLE
                  ENDIF

                  ix = xw(ilon)
                  DO WHILE (.true.)

                     IF (ix == xw(ilon)) THEN
                        lon_w = pixel%lon_w(ilon)
                     ELSE
                        lon_w = fgrid%lon_w(ix)
                     ENDIF

                     IF (ix == xe(ilon)) THEN
                        lon_e = pixel%lon_e(ilon)
                     ELSE
                        lon_e = fgrid%lon_e(ix)
                     ENDIF

                     skip = .false.
                     IF (.not. (lon_between_floor (lon_w, pixel%lon_w(ilon), lon_e) &
                        .and. lon_between_ceil (lon_e, lon_w, pixel%lon_e(ilon)))) THEN
                        skip = .true.
                     ELSE
                        IF (lon_e > lon_w) THEN
                           IF ((lon_e-lon_w) < 1.0e-6_r8) THEN
                              skip = .true.
                           ENDIF
                        ELSE
                           IF ((lon_e+360.0_r8-lon_w) < 1.0e-6_r8) THEN
                              skip = .true.
                           ENDIF
                        ENDIF
                     ENDIF

                     IF (.not. skip) THEN

                        area = areaquad (lat_s, lat_n, lon_w, lon_e)

                        CALL insert_into_sorted_list2 ( ix, iy, &
                           gpfrom(iset)%ng, gpfrom(iset)%ilon, gpfrom(iset)%ilat, &
                           iloc, is_new)

                        IF (is_new) THEN
                           IF (iloc < gpfrom(iset)%ng) THEN
                              pafrac(iset)%val(iloc+1:gpfrom(iset)%ng) &
                                 = pafrac(iset)%val(iloc:gpfrom(iset)%ng-1)
                           ENDIF

                           pafrac(iset)%val(iloc) = area
                        ELSE
                           pafrac(iset)%val(iloc) = pafrac(iset)%val(iloc) + area
                        ENDIF

                        IF (gpfrom(iset)%ng == size(gpfrom(iset)%ilat)) THEN
                           CALL expand_list (gpfrom(iset)%ilat, 0.2_r8)
                           CALL expand_list (gpfrom(iset)%ilon, 0.2_r8)
                           CALL expand_list (pafrac(iset)%val,  0.2_r8)
                        ENDIF

                        CALL insert_into_sorted_list1 ( &
                           ix, ng_lat(iy), list_lat(iy)%val, iloc)

                        IF (ng_lat(iy) == size(list_lat(iy)%val)) THEN
                           CALL expand_list (list_lat(iy)%val, 0.2_r8)
                        ENDIF

                     ENDIF

                     IF (ix == xe(ilon))  EXIT
                     ix = mod(ix,fgrid%nlon) + 1
                  ENDDO
               ENDDO

            ENDDO
         ENDDO

         deallocate (ys)
         deallocate (yn)
         deallocate (xw)
         deallocate (xe)

         ng_all = sum(ng_lat)
         allocate (xlist(ng_all))
         allocate (ylist(ng_all))

         ig = 0
         DO iy = 1, fgrid%nlat
            DO ix = 1, ng_lat(iy)
               ig = ig + 1
               xlist(ig) = list_lat(iy)%val(ix)
               ylist(ig) = iy
            ENDDO
         ENDDO

         deallocate (ng_lat)
         DO iy = 1, fgrid%nlat
            deallocate (list_lat(iy)%val)
         ENDDO
         deallocate (list_lat)

#ifdef USEMPI
         allocate (ipt (ng_all))
         allocate (msk (ng_all))
         DO ig = 1, ng_all
            xblk = fgrid%xblk(xlist(ig))
            yblk = fgrid%yblk(ylist(ig))
            ipt(ig) = gblock%pio(xblk,yblk)
         ENDDO
#endif

         allocate (gplist (0:p_np_io-1))
         DO iproc = 0, p_np_io-1
#ifdef USEMPI
            msk = (ipt == p_address_io(iproc))
            ng  = count(msk)
#else
            ng  = ng_all
#endif

            gplist(iproc)%ng = ng

            IF (ng > 0) THEN
               allocate (gplist(iproc)%ilat (ng))
               allocate (gplist(iproc)%ilon (ng))
               allocate (gplist(iproc)%elmid(ng))

#ifdef USEMPI
               gplist(iproc)%ilon = pack(xlist, msk)
               gplist(iproc)%ilat = pack(ylist, msk)
#else
               gplist(iproc)%ilon = xlist
               gplist(iproc)%ilat = ylist
#endif
            ENDIF
         ENDDO

#ifdef USEMPI
         deallocate (ipt)
         deallocate (msk)
#endif
         deallocate (xlist)
         deallocate (ylist)         

         allocate (paddress  (pixelset%nset))

         DO iset = 1, pixelset%nset
    
            ng = gpfrom(iset)%ng
            allocate (paddress(iset)%val (2,ng))

            DO ig = 1, gpfrom(iset)%ng
               ilon = gpfrom(iset)%ilon(ig)
               ilat = gpfrom(iset)%ilat(ig)
               xblk = fgrid%xblk(ilon)
               yblk = fgrid%yblk(ilat)

#ifdef USEMPI
               iproc = p_itis_io(gblock%pio(xblk,yblk))
#else
               iproc = 0
#endif

               paddress(iset)%val(1,ig) = iproc
               paddress(iset)%val(2,ig) = find_in_sorted_list2 ( &
                  ilon, ilat, gplist(iproc)%ng, gplist(iproc)%ilon, gplist(iproc)%ilat)
            ENDDO

            deallocate (gpfrom(iset)%ilon)
            deallocate (gpfrom(iset)%ilat)
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iproc = 0, p_np_io-1
            idest = p_address_io(iproc)
            smesg = (/p_iam_glb, gplist(iproc)%ng/)

            CALL mpi_send (smesg, 2, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err)

            IF (gplist(iproc)%ng > 0) THEN
               CALL mpi_send (gplist(iproc)%ilon, gplist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
               if(p_err == 0) deallocate(gplist(iproc)%ilon)
               CALL mpi_send (gplist(iproc)%ilat, gplist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
               if(p_err == 0) deallocate(gplist(iproc)%ilat)
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_io) THEN

         allocate (gplist (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)
            iproc = p_itis_worker(isrc)

            gplist(iproc)%ng = nrecv

            IF (nrecv > 0) THEN
               allocate (gplist(iproc)%ilon (nrecv))
               allocate (gplist(iproc)%ilat (nrecv))
               allocate (gplist(iproc)%elmid(nrecv))

               CALL mpi_recv (gplist(iproc)%ilon, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (gplist(iproc)%ilat, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO ig = 1, nrecv
                  ilon = gplist(iproc)%ilon(ig)
                  ilat = gplist(iproc)%ilat(ig)
                  xblk = fgrid%xblk (ilon)
                  yblk = fgrid%yblk (ilat)
                  xloc = fgrid%xloc (ilon)
                  yloc = fgrid%yloc (ilat)

                  gplist(iproc)%elmid(ig) = elmid_atm%blk(xblk,yblk)%val(xloc,yloc)
               END DO
               
               deallocate (gplist(iproc)%ilon)
               deallocate (gplist(iproc)%ilat)
            ENDIF
         ENDDO

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            idest = p_address_worker(iproc)
            smesg = (/p_iam_glb, gplist(iproc)%ng/)

            CALL mpi_send (smesg, 2, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err)

            IF (gplist(iproc)%ng > 0) THEN
               CALL mpi_send (gplist(iproc)%elmid, gplist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
               if(p_err == 0) deallocate (gplist(iproc)%elmid)
            ENDIF
         END DO
         deallocate (gplist)

      ENDIF

      IF (p_is_worker) THEN

         DO iio = 0, p_np_io-1

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)
            iproc = p_itis_io(isrc)

            IF (nrecv > 0) THEN
               CALL mpi_recv (gplist(iproc)%elmid, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
            ENDIF
         ENDDO

      ENDIF
#endif

      IF (p_is_worker) THEN

         allocate (afrac (pixelset%nset))
         allocate (gfrom (pixelset%nset))

         DO iset = 1, pixelset%nset

            ng = gpfrom(iset)%ng
            allocate (gpfrom(iset)%elmid(ng))

            DO ig = 1, ng
               iproc = paddress(iset)%val(1,ig)
               iloc  = paddress(iset)%val(2,ig)
               gpfrom(iset)%elmid(ig) = gplist(iproc)%elmid(iloc)
            ENDDO

            deallocate (paddress(iset)%val)
            allocate (afrac(iset)%val  (ng))
            allocate (gfrom(iset)%elmid(ng))
            allocate (msk(ng))

            gfrom(iset)%ng = 0

            DO ig = 1, ng

               IF (gpfrom(iset)%elmid(ig) > 0) then

                  msk = (gpfrom(iset)%elmid == gpfrom(iset)%elmid(ig))
                  CALL insert_into_sorted_list1 ( &
                     gpfrom(iset)%elmid(ig), gfrom(iset)%ng, gfrom(iset)%elmid, iloc, is_new)

                  IF (is_new) THEN
                     IF (iloc < gfrom(iset)%ng) THEN
                        afrac(iset)%val(iloc+1:gfrom(iset)%ng) &
                           = afrac(iset)%val(iloc:gfrom(iset)%ng-1)
                     ENDIF

                     afrac(iset)%val(iloc) = sum(pafrac(iset)%val(1:ng),msk)
                  ELSE
                     write(*,*)"gpfrom elmid double count,error",gpfrom(iset)%elmid(ig);stop
                  ENDIF

                  WHERE(msk) gpfrom(iset)%elmid = -1
               ENDIF

            ENDDO

            deallocate (msk)
            deallocate (pafrac(iset)%val)
            deallocate (gpfrom(iset)%elmid)

         ENDDO
         deallocate (paddress)
         deallocate (pafrac)

         DO iproc = 0, p_np_io-1
            IF (gplist(iproc)%ng > 0) deallocate (gplist(iproc)%elmid)
         ENDDO

         deallocate (gplist)
         deallocate (gpfrom)

      ENDIF

ELSE

      IF (p_is_worker) THEN

         this%npset = pixelset%nset
         allocate (afrac (pixelset%nset))
         allocate (gfrom (pixelset%nset))

         DO iset = 1, pixelset%nset

            ng = 1
            allocate (afrac(iset)%val  (ng))
            allocate (gfrom(iset)%elmid(ng))

            gfrom(iset)%ng = ng
            gfrom(iset)%elmid = pixelset%eindex(iset)
            afrac(iset)%val   = 0.0

            ie = pixelset%ielm(iset)

            ipxstt = pixelset%ipxstt(iset)
            ipxend = pixelset%ipxend(iset)

            ! deal with 2m WMO patch
            IF (ipxstt==-1 .and. ipxend==-1) THEN
               ipxstt = 1
               ipxend = mesh(ie)%npxl
            ENDIF

            DO ipxl = ipxstt, ipxend
               afrac(iset)%val = afrac(iset)%val + areaquad (&
                  pixel%lat_s(mesh(ie)%ilat(ipxl)), &
                  pixel%lat_n(mesh(ie)%ilat(ipxl)), &
                  pixel%lon_w(mesh(ie)%ilon(ipxl)), &
                  pixel%lon_e(mesh(ie)%ilon(ipxl)) )
            ENDDO

         ENDDO

      ENDIF

ENDIF

      IF (p_is_worker) THEN 

         allocate (elist(1000))
         ng_all = 0

         DO iset = 1, pixelset%nset
            DO ig = 1, gfrom(iset)%ng

               CALL insert_into_sorted_list1 ( &
                  gfrom(iset)%elmid(ig), ng_all, elist, iloc)
               IF (ng_all == size(elist)) CALL expand_list (elist, 0.2_r8)

            ENDDO
         ENDDO

#ifdef USEMPI
         allocate (ipt (ng_all))
         allocate (msk (ng_all))
         DO ig = 1, ng_all
            ipt(ig) = mesh_atm_pio(elist(ig))
         ENDDO
#endif

         allocate (this%glist (0:p_np_io-1))
         DO iproc = 0, p_np_io-1
#ifdef USEMPI
            msk = (ipt == p_address_io(iproc))
            ng  = count(msk)
#else
            ng  = ng_all
#endif

            this%glist(iproc)%ng = ng

            IF (ng > 0) THEN
               allocate (this%glist(iproc)%elmid(ng))

#ifdef USEMPI
               this%glist(iproc)%elmid = pack(elist(1:ng_all), msk)
#else
               this%glist(iproc)%elmid = elist(1:ng_all)
#endif
            ENDIF
         ENDDO 

#ifdef USEMPI
         deallocate (ipt)
         deallocate (msk)
#endif
         deallocate (elist)

         allocate (this%address  (pixelset%nset))
         allocate (this%areapart (pixelset%nset))

         allocate (this%npart (pixelset%nset))

         DO iset = 1, pixelset%nset

            ng = gfrom(iset)%ng

            this%npart(iset) = ng

            allocate (this%address(iset)%val (2,ng))
            allocate (this%areapart(iset)%val (ng))

            this%areapart(iset)%val = afrac(iset)%val(1:ng)

            IF (pixelset%has_shared) THEN
               this%areapart(iset)%val = this%areapart(iset)%val * pixelset%pctshared(iset)
            ENDIF

            DO ig = 1, gfrom(iset)%ng
#ifdef USEMPI
               iproc = p_itis_io(mesh_atm_pio(gfrom(iset)%elmid(ig)))
#else
               iproc = 0
#endif

               this%address(iset)%val(1,ig) = iproc
               this%address(iset)%val(2,ig) = find_in_sorted_list1 ( &
                  gfrom(iset)%elmid(ig), this%glist(iproc)%ng, this%glist(iproc)%elmid)
            ENDDO
         ENDDO

         DO iset = 1, pixelset%nset
            deallocate (afrac(iset)%val  )
            deallocate (gfrom(iset)%elmid)
         ENDDO

         deallocate (afrac)
         deallocate (gfrom)

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN

         DO iproc = 0, p_np_io-1
            idest = p_address_io(iproc)
            smesg = (/p_iam_glb, this%glist(iproc)%ng/)

            CALL mpi_send (smesg, 2, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err)

            IF (this%glist(iproc)%ng > 0) THEN
               CALL mpi_send (this%glist(iproc)%elmid, this%glist(iproc)%ng, MPI_INTEGER, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_io) THEN

         DO ie = 1, numelm_atm
            mesh_atm(ie)%landmask = 0
         ENDDO

         allocate (this%glist (0:p_np_worker-1))

         DO iworker = 0, p_np_worker-1

            CALL mpi_recv (rmesg, 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)
            iproc = p_itis_worker(isrc)

            this%glist(iproc)%ng = nrecv

            IF (nrecv > 0) THEN
               allocate (this%glist(iproc)%elmid (nrecv))

               CALL mpi_recv (this%glist(iproc)%elmid, nrecv, MPI_INTEGER, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)

               DO ig = 1, nrecv
                  iloc = findloc_ud (this%elmid == this%glist(iproc)%elmid(ig))
                  mesh_atm(iloc)%landmask = 1
               ENDDO
            ENDIF

         ENDDO

      ENDIF
#endif

      IF (p_is_worker) THEN
         IF (this%npset > 0) THEN
            allocate (this%areapset (this%npset))
            this%areapset(:) = 0.
         ENDIF
         DO iset = 1, this%npset
            IF (this%npart(iset) > 0) THEN
               this%areapset(iset) = sum(this%areapart(iset)%val)
            ENDIF
         ENDDO
      ENDIF

      CALL this%get_sumarea (this%areagrid, numelm_atm, filter)

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE spatial_mapping_build_arealweighted

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_pset2grid (this, pdata, gdata, spv, msk)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_UserDefFun
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8), intent(in) :: pdata(:)
   real(r8), allocatable, intent(inout) :: gdata(:)

   real(r8), intent(in), optional :: spv
   logical,  intent(in), optional :: msk(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

               IF (present(spv)) THEN
                  pbuff(iproc)%val(:) = spv
               ELSE
                  pbuff(iproc)%val(:) = 0.0
               ENDIF
            ENDIF
         ENDDO

         DO iset = 1, this%npset

            IF (present(spv)) THEN
               IF (pdata(iset) == spv) CYCLE
            ENDIF

            IF (present(msk)) THEN
               IF (.not. msk(iset)) CYCLE
            ENDIF

            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               IF (present(spv)) THEN
                  IF (pbuff(iproc)%val(iloc) /= spv) THEN
                     pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) &
                        + pdata(iset) * this%areapart(iset)%val(ipart)
                  ELSE
                     pbuff(iproc)%val(iloc) = &
                        pdata(iset) * this%areapart(iset)%val(ipart)
                  ENDIF
               ELSE
                  pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) &
                     + pdata(iset) * this%areapart(iset)%val(ipart)
               ENDIF
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         IF (present(spv)) THEN
            gdata = spv
         ELSE
            gdata = 0.0_r8
         ENDIF

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  iloc = findloc_ud (this%elmid == this%glist(iproc)%elmid(ig))

                  IF (present(spv)) THEN
                     IF (gbuff(ig) /= spv) THEN
                        IF (gdata(iloc) /= spv) THEN
                           gdata(iloc) = gdata(iloc) + gbuff(ig)
                        ELSE
                           gdata(iloc) = gbuff(ig)
                        ENDIF
                     ENDIF
                  ELSE
                     gdata(iloc) = gdata(iloc) + gbuff(ig)
                  ENDIF
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO

         IF (.not. present(spv)) THEN

            WHERE (this%areagrid > 0.)
               gdata = gdata / this%areagrid
            ELSEWHERE
               gdata = spval
            ENDWHERE

         ELSE

            WHERE ((this%areagrid > 0.) .and. (gdata /= spv))
               gdata = gdata / this%areagrid
            ELSEWHERE
               gdata = spval
            ENDWHERE

         ENDIF

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE spatial_mapping_pset2grid

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_pset2grid_max (this, pdata, gdata, spv, msk)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_UserDefFun
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8), intent(in) :: pdata(:)
   real(r8), allocatable, intent(inout) :: gdata(:)

   real(r8), intent(in), optional :: spv
   logical,  intent(in), optional :: msk(:)

   ! Local variables
   integer  :: iproc, idest, isrc
   integer  :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))
               pbuff(iproc)%val(:) = spval
            ENDIF
         ENDDO

         DO iset = 1, this%npset

            IF (present(spv)) THEN
               IF (pdata(iset) == spv) CYCLE
            ENDIF

            IF (present(msk)) THEN
               IF (.not. msk(iset)) CYCLE
            ENDIF

            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               IF (pbuff(iproc)%val(iloc) /= spval) THEN
                  pbuff(iproc)%val(iloc) = max(pdata(iset), pbuff(iproc)%val(iloc))
               ELSE
                  pbuff(iproc)%val(iloc) = pdata(iset)
               ENDIF
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         gdata = spval

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  IF (gbuff(ig) /= spval) THEN
                     iloc = findloc_ud (this%elmid == this%glist(iproc)%elmid(ig))

                     IF (gdata(iloc) /= spval) THEN
                        gdata(iloc) = max(gdata(iloc), gbuff(ig))
                     ELSE
                        gdata(iloc) = gbuff(ig)
                     ENDIF
                  ENDIF
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE spatial_mapping_pset2grid_max

   ! ------------------------------
   SUBROUTINE spatial_mapping_get_sumarea (this, sumarea, numelm_atm, filter)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_UserDefFun
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8), allocatable, intent(out) :: sumarea(:)
   integer, intent(in) :: numelm_atm
   logical, intent(in), optional :: filter(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))
               pbuff(iproc)%val(:) = 0.0
            ENDIF
         ENDDO

         DO iset = 1, this%npset

            IF (present(filter)) THEN
               IF (.not. filter(iset)) CYCLE
            ENDIF

            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)
               pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) + this%areapart(iset)%val(ipart)
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io .and. numelm_atm > 0) THEN

         allocate (sumarea(numelm_atm))
         sumarea = 0.0_r8

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  iloc = findloc_ud (this%elmid == this%glist(iproc)%elmid(ig))

                  sumarea(iloc) = sumarea(iloc) + gbuff(ig)
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE spatial_mapping_get_sumarea

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_grid2pset (this, gdata, pdata)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8), intent(in)  :: gdata(:)
   real(r8), intent(out) :: pdata(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  iloc = findloc_ud (this%elmid == this%glist(iproc)%elmid(ig))
                  gbuff(ig) = gdata(iloc)
               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

         DO iset = 1, this%npset

            IF (this%areapset(iset) > 0.) THEN

               pdata(iset) = 0.

               DO ipart = 1, this%npart(iset)
                  iproc = this%address(iset)%val(1,ipart)
                  iloc  = this%address(iset)%val(2,ipart)

                  IF (this%areapart(iset)%val(ipart) > 0.) THEN
                     pdata(iset) = pdata(iset) &
                        + pbuff(iproc)%val(iloc) * this%areapart(iset)%val(ipart)
                  ENDIF
               ENDDO

               pdata(iset) = pdata(iset) / this%areapset(iset)

            ELSE
               pdata(iset) = spval
            ENDIF

         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)

      ENDIF

   END SUBROUTINE spatial_mapping_grid2pset

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_grid2part (this, gdata, sdata)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8),                  intent(in)    :: gdata(:)
   type(pointer_real8_1d),    intent(inout) :: sdata(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_io) THEN

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

               DO ig = 1, this%glist(iproc)%ng
                  iloc = findloc_ud (this%elmid == this%glist(iproc)%elmid(ig))
                  gbuff(ig) = gdata(iloc)
               ENDDO

#ifdef USEMPI
               idest = p_address_worker(iproc)
               CALL mpi_send (gbuff, this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)

               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

      ENDIF

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_io(iproc)
               CALL mpi_recv (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               pbuff(0)%val = gbuff
               deallocate (gbuff)
#endif
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               sdata(iset)%val(ipart) = pbuff(iproc)%val(iloc)
            ENDDO
         ENDDO

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)

      ENDIF

   END SUBROUTINE spatial_mapping_grid2part

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_part2grid (this, sdata, gdata)

   USE MOD_Precision
   USE MOD_Block
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(pointer_real8_1d), intent(in)    :: sdata(:)
   real(r8), allocatable,  intent(inout) :: gdata(:)

   ! Local variables
   integer :: iproc, idest, isrc
   integer :: ig, ilon, ilat, xblk, yblk, xloc, yloc, iloc, iset, ipart
   integer :: iblkme

   real(r8), allocatable :: gbuff(:)
   type(pointer_real8_1d), allocatable :: pbuff(:)

      IF (p_is_worker) THEN

         allocate (pbuff (0:p_np_io-1))

         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               allocate (pbuff(iproc)%val (this%glist(iproc)%ng))
               pbuff(iproc)%val(:) = 0.0
            ENDIF
         ENDDO

         DO iset = 1, this%npset
            DO ipart = 1, this%npart(iset)
               iproc = this%address(iset)%val(1,ipart)
               iloc  = this%address(iset)%val(2,ipart)

               pbuff(iproc)%val(iloc) = pbuff(iproc)%val(iloc) &
                  + sdata(iset)%val(ipart) * this%areapart(iset)%val(ipart)
            ENDDO
         ENDDO

#ifdef USEMPI
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               idest = p_address_io(iproc)
               CALL mpi_send (pbuff(iproc)%val, this%glist(iproc)%ng, MPI_DOUBLE, &
                  idest, mpi_tag_data, p_comm_glb, p_err)
            ENDIF
         ENDDO
#endif

      ENDIF

      IF (p_is_io) THEN

         gdata = 0.0_r8

         DO iproc = 0, p_np_worker-1
            IF (this%glist(iproc)%ng > 0) THEN

               allocate (gbuff (this%glist(iproc)%ng))

#ifdef USEMPI
               isrc = p_address_worker(iproc)
               CALL mpi_recv (gbuff, this%glist(iproc)%ng, MPI_DOUBLE, &
                  isrc, mpi_tag_data, p_comm_glb, p_stat, p_err)
#else
               gbuff = pbuff(0)%val
#endif

               DO ig = 1, this%glist(iproc)%ng
                  iloc = findloc_ud (this%elmid == this%glist(iproc)%elmid(ig))
                  gdata(iloc) = gdata(iloc) + gbuff(ig)
               ENDDO

               deallocate (gbuff)
            ENDIF
         ENDDO

         WHERE (this%areagrid > 0)
            gdata = gdata / this%areagrid
         ELSEWHERE
            gdata = this%missing_value
         ENDWHERE

      ENDIF

      IF (p_is_worker) THEN
         DO iproc = 0, p_np_io-1
            IF (this%glist(iproc)%ng > 0) THEN
               deallocate (pbuff(iproc)%val)
            ENDIF
         ENDDO
         deallocate (pbuff)
      ENDIF

   END SUBROUTINE spatial_mapping_part2grid

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_normalize (this, gdata, sdata)

   USE MOD_Precision
   USE MOD_Block
   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   real(r8), allocatable,  intent(in)    :: gdata(:)
   type(pointer_real8_1d), intent(inout) :: sdata(:)

   ! Local variables
   integer :: iblkme, xblk, yblk, iset, ipart

   real(r8), allocatable :: sumdata(:)
   type(pointer_real8_1d), allocatable :: scaldata(:)


      IF (p_is_io)     allocate (sumdata(size(this%elmid)))
      IF (p_is_worker) CALL this%allocate_part  (scaldata)

      CALL this%part2grid (sdata, sumdata)

      IF (p_is_io) THEN

         WHERE (sumdata /= this%missing_value)
            sumdata = gdata / sumdata
         ENDWHERE

      ENDIF

      CALL this%grid2part (sumdata, scaldata)

      IF (p_is_worker) THEN

         DO iset = 1, this%npset
            DO ipart = 1, this%npart(iset)
               IF (this%areapart(iset)%val(ipart) > 0.) THEN
                  sdata(iset)%val(ipart) = sdata(iset)%val(ipart) * scaldata(iset)%val(ipart)
               ELSE
                  sdata(iset)%val(ipart) = this%missing_value
               ENDIF
            ENDDO
         ENDDO

      ENDIF

      IF (p_is_io)     deallocate(sumdata)
      IF (p_is_worker) CALL this%deallocate_part(scaldata)

   END SUBROUTINE spatial_mapping_normalize

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_part2pset (this, sdata, pdata)

   USE MOD_Precision
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_SPMD_Task
   USE MOD_Vars_Global, only: spval
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(pointer_real8_1d), intent(in)  :: sdata(:)
   real(r8),               intent(out) :: pdata(:)

   ! Local variables
   integer :: iset

      IF (p_is_worker) THEN

         pdata(:) = spval

         DO iset = 1, this%npset
            IF (this%areapset(iset) > 0) THEN
               pdata(iset) = sum(sdata(iset)%val * this%areapart(iset)%val) / this%areapset(iset)
            ENDIF
         ENDDO

      ENDIF

   END SUBROUTINE spatial_mapping_part2pset

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_allocate_part (this, datapart)

   USE MOD_SPMD_Task
   USE MOD_DataType
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(pointer_real8_1d), allocatable :: datapart (:)

   ! Local variables
   integer :: iset

      IF (p_is_worker) THEN

         IF (this%npset > 0) THEN
            allocate (datapart (this%npset))
         ENDIF

         DO iset = 1, this%npset
            IF (this%npart(iset) > 0) THEN
               allocate (datapart(iset)%val (this%npart(iset)))
            ENDIF
         ENDDO

      ENDIF

   END SUBROUTINE spatial_mapping_allocate_part

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_deallocate_part (this, datapart)

   USE MOD_SPMD_Task
   USE MOD_DataType
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   type(pointer_real8_1d), allocatable :: datapart (:)

   ! Local variables
   integer :: iset

      IF (p_is_worker) THEN

         DO iset = 1, this%npset
            IF (this%npart(iset) > 0) THEN
               deallocate (datapart(iset)%val)
            ENDIF
         ENDDO

         IF (this%npset > 0) THEN
            deallocate (datapart)
         ENDIF

      ENDIF

   END SUBROUTINE spatial_mapping_deallocate_part

   !-----------------------------------------------------
   SUBROUTINE spatial_mapping_free_mem (this)

   USE MOD_SPMD_Task
   IMPLICIT NONE

   type (spatial_mapping_type) :: this

   ! Local variables
   integer :: iproc, iset

      IF (allocated (this%elmid))   deallocate (this%elmid)

      IF (allocated(this%glist)) THEN
         DO iproc = lbound(this%glist,1), ubound(this%glist,1)
            IF (allocated(this%glist(iproc)%elmid)) deallocate (this%glist(iproc)%elmid)
         ENDDO

         deallocate (this%glist)
      ENDIF

      IF (p_is_worker) THEN

         IF (allocated(this%npart)) deallocate(this%npart)

         IF (allocated(this%address)) THEN
            DO iset = lbound(this%address,1), ubound(this%address,1)
               IF (allocated(this%address(iset)%val)) THEN
                  deallocate (this%address(iset)%val)
               ENDIF
            ENDDO

            deallocate (this%address)
         ENDIF

         IF (allocated(this%areapart)) THEN
            DO iset = lbound(this%areapart,1), ubound(this%areapart,1)
               IF (allocated(this%areapart(iset)%val)) THEN
                  deallocate (this%areapart(iset)%val)
               ENDIF
            ENDDO

            deallocate (this%areapart)
         ENDIF

         IF (allocated(this%areapset)) deallocate(this%areapset)

      ENDIF

      IF (allocated (this%areagrid)) deallocate(this%areagrid)

   END SUBROUTINE spatial_mapping_free_mem

   SUBROUTINE forc_free_mem_spatial_mapping(this)

   USE MOD_SPMD_Task
   IMPLICIT NONE

   class (spatial_mapping_type) :: this

   ! Local variables
   integer :: iproc, iset

      IF (allocated (this%elmid))   deallocate (this%elmid)

      IF (allocated(this%glist)) THEN
         DO iproc = lbound(this%glist,1), ubound(this%glist,1)
            IF (allocated(this%glist(iproc)%elmid)) deallocate (this%glist(iproc)%elmid)
         ENDDO

         deallocate (this%glist)
      ENDIF

      IF (p_is_worker) THEN

         IF (allocated(this%npart)) deallocate(this%npart)

         IF (allocated(this%address)) THEN
            DO iset = lbound(this%address,1), ubound(this%address,1)
               IF (allocated(this%address(iset)%val)) THEN
                  deallocate (this%address(iset)%val)
               ENDIF
            ENDDO

            deallocate (this%address)
         ENDIF

         IF (allocated(this%areapart)) THEN
            DO iset = lbound(this%areapart,1), ubound(this%areapart,1)
               IF (allocated(this%areapart(iset)%val)) THEN
                  deallocate (this%areapart(iset)%val)
               ENDIF
            ENDDO

            deallocate (this%areapart)
         ENDIF

         IF (allocated(this%areapset)) deallocate(this%areapset)

      ENDIF

      IF (allocated (this%areagrid)) deallocate(this%areagrid)

   END SUBROUTINE forc_free_mem_spatial_mapping

END MODULE MOD_SpatialMapping_atm

#endif
