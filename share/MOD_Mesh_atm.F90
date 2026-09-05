#include <define.h>
#ifdef CCPL

MODULE MOD_Mesh_atm

!------------------------------------------------------------------------------------
! !DESCRIPTION:
!
!    MESH_atm refers to the atmospheric grids in atmospheric models or forcing data.
!    This is particularly needed if 
!    (1) in offline modes, CoLM uses atmospheric forcing data with irregular grids.
!    (2) in online modes, CoLM imports and exports variables as atmospheric grids.
!
!  Created by Nan Wei, Feb 2025, template from MOD_Mesh.F90
!------------------------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE

   ! ---- data types ----
   type :: irregular_elm_type_atm

      integer*8 :: indx
      integer   :: xblk, yblk

      integer :: npxl
      integer, allocatable :: ilon(:)
      integer, allocatable :: ilat(:)

      integer :: landmask

   END type irregular_elm_type_atm

CONTAINS

   SUBROUTINE copy_elm_atm (elm_from, elm_to)

   IMPLICIT NONE
   type (irregular_elm_type_atm), intent(in)  :: elm_from
   type (irregular_elm_type_atm), intent(out) :: elm_to

      elm_to%indx = elm_from%indx
      elm_to%npxl = elm_from%npxl
      elm_to%xblk = elm_from%xblk
      elm_to%yblk = elm_from%yblk

      IF (allocated(elm_to%ilat)) deallocate(elm_to%ilat)
      IF (allocated(elm_to%ilon)) deallocate(elm_to%ilon)

      allocate (elm_to%ilat (elm_to%npxl))
      allocate (elm_to%ilon (elm_to%npxl))
      elm_to%ilon = elm_from%ilon
      elm_to%ilat = elm_from%ilat

      elm_to%landmask = elm_from%landmask

   END SUBROUTINE copy_elm_atm

   ! --------------------------------
   SUBROUTINE mesh_atm_build (gforc, DEF_file_mesh_atm, numelm_atm, mesh_atm, nelm_blk_atm, mesh_atm_pio, elmid_atm)

   USE MOD_SPMD_Task
   USE MOD_NetCDFBlock
   USE MOD_Block
   USE MOD_Grid
   USE MOD_Utils
   USE MOD_DataType

   IMPLICIT NONE

   type (grid_type),  intent(in) :: gforc
   character(len=256),intent(in) :: DEF_file_mesh_atm

   integer,intent(out) :: numelm_atm                                        ! only available on IO processes
   type (irregular_elm_type_atm), intent(out), allocatable :: mesh_atm (:)  ! only available on IO processes

   integer, intent(out), allocatable :: nelm_blk_atm(:,:)
   integer, intent(out), allocatable :: mesh_atm_pio(:)
   type(block_data_int32_2d), intent(out) :: elmid_atm

   ! Local Variables
   integer  :: nelm, ie, je, dsp
   integer  :: iblkme, iblk, jblk, xloc, yloc
   integer  :: xp, yp, xblk, yblk, npxl, ipxl, ix, iy
   integer  :: iworker, iproc, idest, isrc, iloc, iloc_max(2)
   integer  :: ysp, ynp, nyp
   integer  :: xwp, xep, nxp
   logical  :: is_new
   integer  :: smesg(6), rmesg(6), blktag, nsend, nrecv, irecv
   integer  :: iblk_p, jblk_p, nelm_glb

   integer, allocatable :: nelm_worker(:)

   integer*8 :: elmid
   integer*8, allocatable :: elist(:), elist2(:,:), sbuf64(:), elist_recv(:)

   integer, allocatable :: iaddr(:), elmindx(:), order(:)
   integer, allocatable :: xlist(:), ylist(:), npxl_(:), xlist_recv(:), ylist_recv(:), sbuf(:)
   logical, allocatable :: msk(:), work_done(:)

   integer, allocatable :: xlist2(:,:), ylist2(:,:), ipt2(:,:), npxl_blk(:,:), blkdsp(:,:), blkcnt(:,:)
   logical, allocatable :: msk2(:,:)

   type(irregular_elm_type_atm), allocatable :: meshtmp (:)

      IF (p_is_io) THEN
         CALL allocate_block_data (gforc, elmid_atm)
      ENDIF

      CALL ncio_read_block (DEF_file_mesh_atm, 'elmindex', gforc, elmid_atm)

      ! Step 1: How many atm elms in each block?
      IF (p_is_io) THEN

         nelm = 0

         allocate (nelm_worker (0:p_np_worker-1))
         nelm_worker(:) = 0

         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)

            allocate (msk2(gforc%xcnt(iblk),gforc%ycnt(jblk)))
            allocate (elist2(gforc%xcnt(iblk),gforc%ycnt(jblk)))
            elist2 = elmid_atm%blk(iblk,jblk)%val

            DO yloc = 1, gforc%ycnt(jblk)
               DO xloc = 1, gforc%xcnt(iblk)

                  elmid = elist2(xloc,yloc)

                  IF (elmid > 0) THEN

                     iworker = mod(elmid, p_np_worker)
                     nelm_worker(iworker) = nelm_worker(iworker) + 1
 
                     msk2 = (elist2 == elmid)
                     WHERE(msk2) elist2 = -1

                  ENDIF

               ENDDO
            ENDDO

            deallocate (msk2)
            deallocate (elist2)

#ifdef USEMPI
            DO iworker = 0, p_np_worker-1
               IF (nelm_worker(iworker) > 0) THEN
                  idest = p_address_worker(iworker)
                  smesg(1:2) = (/p_iam_glb, nelm_worker(iworker)/)
                  ! send(01)
                  CALL mpi_send (smesg(1:2), 2, MPI_INTEGER, &
                     idest, mpi_tag_size, p_comm_glb, p_err)
               ENDIF
            ENDDO
#endif

            nelm = nelm + sum(nelm_worker)
            nelm_worker(:) = 0
         ENDDO

#ifdef USEMPI
         DO iworker = 0, p_np_worker-1
            idest = p_address_worker(iworker)
            ! send(02)
            smesg(1:2) = (/p_iam_glb, 0/)
            CALL mpi_send (smesg(1:2), 2, MPI_INTEGER, &
               idest, mpi_tag_size, p_comm_glb, p_err)
         ENDDO
#endif

         deallocate (nelm_worker)

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         nelm = 0
         allocate(work_done(0:p_np_io-1))
         work_done(:) = .false.
         DO WHILE (.not. all(work_done))
            ! recv(01,02)
            CALL mpi_recv (rmesg(1:2), 2, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_size, p_comm_glb, p_stat, p_err)

            isrc  = rmesg(1)
            nrecv = rmesg(2)

            IF (nrecv > 0) THEN
               nelm = nelm + nrecv
            ELSE
               work_done(p_itis_io(isrc)) = .true.
            ENDIF
         ENDDO

         deallocate(work_done)
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! Step 2: Build atm pixel list for each atm elm.
      IF (p_is_worker) THEN
         IF (nelm > 0) THEN
            allocate (meshtmp (nelm))
            allocate (elist (nelm))
            allocate (iaddr (nelm))
         ENDIF
         nelm = 0
      ENDIF

      IF (p_is_io) THEN

         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            IF (gforc%xcnt(iblk) <= 0) CYCLE
            IF (gforc%ycnt(jblk) <= 0) CYCLE

            ysp = gforc%ydsp(jblk) + 1
            ynp = gforc%ydsp(jblk) + gforc%ycnt(jblk)
            nyp = ynp - ysp + 1

            xwp = gforc%xdsp(iblk) + 1
            xep = gforc%xdsp(iblk) + gforc%xcnt(iblk)
            IF (xep > gforc%nlon) xep = xep - gforc%nlon
            nxp = xep - xwp + 1
            IF (nxp <= 0) nxp = nxp + gforc%nlon

            allocate (elist2 (nxp,nyp))
            allocate (xlist2 (nxp,nyp))
            allocate (ylist2 (nxp,nyp))
            allocate (msk2   (nxp,nyp))

            DO iy = ysp, ynp
               yloc = gforc%yloc(iy)

               ix = xwp
               DO WHILE (.true.)
                  xloc = gforc%xloc(ix)

                  elmid = elmid_atm%blk(iblk,jblk)%val(xloc,yloc)

                  xlist2(xloc,yloc) = ix
                  ylist2(xloc,yloc) = iy
                  elist2(xloc,yloc) = elmid

                  IF (ix == xep) EXIT
                  ix = mod(ix,gforc%nlon) + 1
               ENDDO
            ENDDO

#ifdef USEMPI
            allocate (sbuf (nxp*nyp))
            allocate (ipt2 (nxp,nyp))
            allocate (sbuf64 (nxp*nyp))

            blktag = iblkme
            ipt2 = mod(elist2, p_np_worker)
            DO iproc = 0, p_np_worker-1
               msk2  = (ipt2 == iproc) .and. (elist2 > 0)
               nsend = count(msk2)
               IF (nsend > 0) THEN

                  idest = p_address_worker(iproc)

                  smesg(1:3) = (/p_iam_glb, nsend, blktag/)
                  ! send(03)
                  CALL mpi_send (smesg(1:3), 3, MPI_INTEGER, &
                     idest, mpi_tag_mesg, p_comm_glb, p_err)

                  sbuf64(1:nsend) = pack(elist2, msk2)
                  ! send(04)
                  CALL mpi_send (sbuf64(1:nsend), nsend, MPI_INTEGER8, &
                     idest, blktag, p_comm_glb, p_err)

                  sbuf(1:nsend) = pack(xlist2, msk2)
                  ! send(05)
                  CALL mpi_send (sbuf(1:nsend), nsend, MPI_INTEGER, &
                     idest, blktag, p_comm_glb, p_err)

                  sbuf(1:nsend) = pack(ylist2, msk2)
                  ! send(06)
                  CALL mpi_send (sbuf(1:nsend), nsend, MPI_INTEGER, &
                     idest, blktag, p_comm_glb, p_err)

               ENDIF
            ENDDO

            deallocate (sbuf  )
            deallocate (ipt2  )
            deallocate (sbuf64)
#endif

            deallocate (elist2)
            deallocate (xlist2)
            deallocate (ylist2)
            deallocate (msk2  )

         ENDDO

#ifdef USEMPI
         DO iworker = 0, p_np_worker-1
            idest = p_address_worker(iworker)
            ! send(07)
            smesg(1:3) = (/p_iam_glb, 0, 0/)
            CALL mpi_send (smesg(1:3), 3, MPI_INTEGER, &
               idest, mpi_tag_mesg, p_comm_glb, p_err)
         ENDDO
#endif

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN

         allocate(work_done(0:p_np_io-1))
         work_done(:) = .false.
         DO WHILE (.not. all(work_done))
            ! recv(03,07)
            CALL mpi_recv (rmesg(1:3), 3, MPI_INTEGER, &
               MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)

            isrc   = rmesg(1)
            nrecv  = rmesg(2)
            blktag = rmesg(3)
            IF (nrecv > 0) THEN

               allocate (elist_recv (nrecv))
               ! recv(04)
               CALL mpi_recv (elist_recv, nrecv, MPI_INTEGER8, &
                  isrc, blktag, p_comm_glb, p_stat, p_err)

               allocate (xlist_recv (nrecv))
               ! recv(05)
               CALL mpi_recv (xlist_recv, nrecv, MPI_INTEGER, &
                  isrc, blktag, p_comm_glb, p_stat, p_err)

               allocate (ylist_recv (nrecv))
               ! recv(06)
               CALL mpi_recv (ylist_recv, nrecv, MPI_INTEGER, &
                  isrc, blktag, p_comm_glb, p_stat, p_err)

               allocate (msk(nrecv))

               DO irecv = 1, nrecv

                  elmid = elist_recv(irecv)

                  IF (elmid > 0) THEN

                     CALL insert_into_sorted_list1 (elmid, nelm, elist, iloc, is_new)

                     msk  = (elist_recv == elmid)
                     npxl = count(msk)

                     IF (is_new) THEN
                        IF (iloc < nelm) THEN
                           iaddr(iloc+1:nelm) = iaddr(iloc:nelm-1)
                        ENDIF
                        iaddr(iloc) = nelm

                        meshtmp(iaddr(iloc))%indx = elmid
                        meshtmp(iaddr(iloc))%npxl = npxl
                     ELSE
                        meshtmp(iaddr(iloc))%npxl = meshtmp(iaddr(iloc))%npxl + npxl
                     ENDIF

                     allocate (xlist(npxl))
                     allocate (ylist(npxl))
                     xlist = pack(xlist_recv, msk)
                     ylist = pack(ylist_recv, msk)

                     CALL append_to_list (meshtmp(iaddr(iloc))%ilon, xlist)
                     CALL append_to_list (meshtmp(iaddr(iloc))%ilat, ylist)

                     WHERE(msk) elist_recv = -1
                     deallocate (xlist)
                     deallocate (ylist)
                  ENDIF

               ENDDO

               deallocate (msk)
               deallocate (elist_recv)
               deallocate (xlist_recv)
               deallocate (ylist_recv)
            ELSE
               work_done(p_itis_io(isrc)) = .true.
            ENDIF
         ENDDO

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      IF (allocated(elist)) deallocate (elist)
      IF (allocated(iaddr)) deallocate (iaddr)

      ! Step 3: Which block each atm elm locates at.
      IF (p_is_worker) THEN

         allocate (npxl_blk (gblock%nxblk,gblock%nyblk))
         allocate (nelm_blk_atm (gblock%nxblk,gblock%nyblk))

         nelm_blk_atm(:,:) = 0

         DO ie = 1, nelm

            npxl_blk (:,:) = 0

            DO ipxl = 1, meshtmp(ie)%npxl
               xp = meshtmp(ie)%ilon(ipxl)
               yp = meshtmp(ie)%ilat(ipxl)

               xblk = gforc%xblk(xp)
               yblk = gforc%yblk(yp)

               npxl_blk(xblk,yblk) = npxl_blk(xblk,yblk) + 1
            ENDDO

            iloc_max = maxloc(npxl_blk)
            meshtmp(ie)%xblk = iloc_max(1)
            meshtmp(ie)%yblk = iloc_max(2)

            nelm_blk_atm(iloc_max(1), iloc_max(2)) = &
               nelm_blk_atm(iloc_max(1), iloc_max(2)) + 1

         ENDDO

         deallocate (npxl_blk)
      ENDIF

#ifdef USEMPI
      IF (.not. p_is_worker) THEN
         allocate (nelm_blk_atm (gblock%nxblk,gblock%nyblk))
         nelm_blk_atm(:,:) = 0
      ENDIF

      CALL mpi_allreduce (MPI_IN_PLACE, nelm_blk_atm, gblock%nxblk*gblock%nyblk, &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)

      allocate(mesh_atm_pio(sum(nelm_blk_atm)))
      mesh_atm_pio = 0

      IF (p_is_worker) THEN
         DO ie = 1, nelm
            mesh_atm_pio(meshtmp(ie)%indx) = gblock%pio(meshtmp(ie)%xblk,meshtmp(ie)%yblk)
         ENDDO
      ENDIF

      CALL mpi_allreduce (MPI_IN_PLACE, mesh_atm_pio, size(mesh_atm_pio), &
         MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)
#endif

      ! Step 4: IF MPI is used, sending atm elms from worker to their IO processes.

      IF (p_is_io) THEN

         allocate (blkdsp (gblock%nxblk, gblock%nyblk))
         blkdsp(1,1) = 0
         DO iblk = 1, gblock%nxblk
            DO jblk = 1, gblock%nyblk
               IF ((iblk /= 1) .or. (jblk /= 1)) THEN
                  IF (jblk == 1) THEN
                     iblk_p = iblk - 1
                     jblk_p = gblock%nyblk
                  ELSE
                     iblk_p = iblk
                     jblk_p = jblk - 1
                  ENDIF

                  IF (gblock%pio(iblk_p,jblk_p) == p_iam_glb) THEN
                     blkdsp(iblk,jblk) = blkdsp(iblk_p,jblk_p) + nelm_blk_atm(iblk_p,jblk_p)
                  ELSE
                     blkdsp(iblk,jblk) = blkdsp(iblk_p,jblk_p)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO

      ENDIF

#ifdef USEMPI
      IF (p_is_worker) THEN
         DO iblk = 1, gblock%nxblk
            DO jblk = 1, gblock%nyblk

               idest = gblock%pio(iblk,jblk)

               nsend = 0
               npxl  = 0
               DO ie = 1, nelm
                  IF ((meshtmp(ie)%xblk == iblk) .and. (meshtmp(ie)%yblk == jblk)) THEN
                     nsend = nsend + 1
                     npxl  = npxl  + meshtmp(ie)%npxl
                  ENDIF
               ENDDO

               IF (nsend > 0) THEN

                  allocate (elist (nsend))
                  allocate (npxl_ (nsend))
                  allocate (xlist (npxl ))
                  allocate (ylist (npxl ))

                  nsend = 0
                  npxl  = 0
                  DO ie = 1, nelm
                     IF ((meshtmp(ie)%xblk == iblk) .and. (meshtmp(ie)%yblk == jblk)) THEN

                        nsend = nsend + 1

                        elist(nsend) = meshtmp(ie)%indx
                        npxl_(nsend) = meshtmp(ie)%npxl

                        xlist(npxl+1:npxl+meshtmp(ie)%npxl) = meshtmp(ie)%ilon
                        ylist(npxl+1:npxl+meshtmp(ie)%npxl) = meshtmp(ie)%ilat

                        npxl = npxl + meshtmp(ie)%npxl
                     ENDIF
                  ENDDO

                  blktag = p_iam_glb + 1000

                  ! send(09)
                  smesg(1:6) = (/p_iam_glb, blktag, iblk, jblk, nsend, npxl/)
                  CALL mpi_send (smesg(1:6), 6, MPI_INTEGER, idest, mpi_tag_mesg, p_comm_glb, p_err)

                  ! send(10)
                  CALL mpi_send (elist, nsend, MPI_INTEGER8, idest, blktag, p_comm_glb, p_err)
                  CALL mpi_send (npxl_, nsend, MPI_INTEGER,  idest, blktag, p_comm_glb, p_err)

                  ! send(11)
                  CALL mpi_send (xlist, npxl, MPI_INTEGER, idest, blktag, p_comm_glb, p_err)
                  CALL mpi_send (ylist, npxl, MPI_INTEGER, idest, blktag, p_comm_glb, p_err)

                  deallocate (elist)
                  deallocate (npxl_)
                  deallocate (xlist)
                  deallocate (ylist)

               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (p_is_io) THEN

         numelm_atm = sum(nelm_blk_atm, mask = gblock%pio == p_iam_glb)

         IF (numelm_atm > 0) THEN

            allocate (mesh_atm (numelm_atm))

            allocate (blkcnt (gblock%nxblk, gblock%nyblk))
            blkcnt(:,:) = 0

            DO WHILE (sum(blkcnt) < numelm_atm)

               ! recv(09)
               CALL mpi_recv (rmesg(1:6), 6, MPI_INTEGER, MPI_ANY_SOURCE, mpi_tag_mesg, p_comm_glb, p_stat, p_err)
               isrc   = rmesg(1)
               blktag = rmesg(2)
               xblk   = rmesg(3)
               yblk   = rmesg(4)
               nrecv  = rmesg(5)
               npxl   = rmesg(6)

               allocate (elist (nrecv))
               allocate (npxl_ (nrecv))
               allocate (xlist (npxl ))
               allocate (ylist (npxl ))

               ! recv(10)
               CALL mpi_recv (elist, nrecv, MPI_INTEGER8, isrc, blktag, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (npxl_, nrecv, MPI_INTEGER,  isrc, blktag, p_comm_glb, p_stat, p_err)

               ! recv(11)
               CALL mpi_recv (xlist, npxl, MPI_INTEGER, isrc, blktag, p_comm_glb, p_stat, p_err)
               CALL mpi_recv (ylist, npxl, MPI_INTEGER, isrc, blktag, p_comm_glb, p_stat, p_err)

               dsp = 0
               DO ie = 1, nrecv

                  je = blkdsp(xblk,yblk) + blkcnt(xblk,yblk) + ie

                  mesh_atm(je)%indx = elist(ie)
                  mesh_atm(je)%xblk = xblk
                  mesh_atm(je)%yblk = yblk
                  mesh_atm(je)%npxl = npxl_(ie)

                  allocate (mesh_atm(je)%ilon (npxl_(ie)))
                  allocate (mesh_atm(je)%ilat (npxl_(ie)))

                  mesh_atm(je)%ilon = xlist(dsp+1:dsp+npxl_(ie))
                  mesh_atm(je)%ilat = ylist(dsp+1:dsp+npxl_(ie))

                  dsp = dsp + npxl_(ie)

               ENDDO

               blkcnt(xblk,yblk) = blkcnt(xblk,yblk) + nrecv

               deallocate (elist)
               deallocate (npxl_)
               deallocate (xlist)
               deallocate (ylist)

            ENDDO

         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)

#endif

      ! Step 4-2: sort elms.
      IF (p_is_io) THEN
         IF (allocated (meshtmp)) THEN
            DO ie = 1, size(meshtmp)
               IF (allocated(meshtmp(ie)%ilon))  deallocate (meshtmp(ie)%ilon)
               IF (allocated(meshtmp(ie)%ilat))  deallocate (meshtmp(ie)%ilat)
            ENDDO
            deallocate (meshtmp)
         ENDIF

         IF (numelm_atm > 0) THEN
            allocate (meshtmp (numelm_atm))
            DO ie = 1, numelm_atm
               CALL copy_elm_atm(mesh_atm(ie), meshtmp(ie))
            ENDDO

            DO iblkme = 1, gblock%nblkme
               iblk = gblock%xblkme(iblkme)
               jblk = gblock%yblkme(iblkme)

               IF (blkcnt(iblk,jblk) > 0) THEN
                  allocate (elmindx (blkcnt(iblk,jblk)))
                  allocate (order   (blkcnt(iblk,jblk)))

                  DO ie = blkdsp(iblk,jblk)+1, blkdsp(iblk,jblk)+blkcnt(iblk,jblk)
                     elmindx(ie-blkdsp(iblk,jblk)) = mesh_atm(ie)%indx
                  ENDDO

                  order = (/ (ie, ie = 1, blkcnt(iblk,jblk)) /)
                  CALL quicksort (blkcnt(iblk,jblk), elmindx, order)

                  DO ie = 1, blkcnt(iblk,jblk)
                     CALL copy_elm_atm (meshtmp(blkdsp(iblk,jblk)+order(ie)), &
                        mesh_atm(blkdsp(iblk,jblk)+ie))
                  ENDDO

                  deallocate (elmindx)
                  deallocate (order  )
               ENDIF

            ENDDO
         ENDIF
      ENDIF

      IF (allocated(blkdsp)) deallocate(blkdsp)
      IF (allocated(blkcnt)) deallocate(blkcnt)

      IF (allocated (meshtmp)) THEN
         DO ie = 1, size(meshtmp)
            IF (allocated(meshtmp(ie)%ilon))  deallocate (meshtmp(ie)%ilon)
            IF (allocated(meshtmp(ie)%ilon))  deallocate (meshtmp(ie)%ilat)
         ENDDO

         deallocate (meshtmp )
      ENDIF

      IF (p_is_master) THEN
         write(*,'(A)') 'Making ATM mesh elements from :'//trim(DEF_file_mesh_atm)
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)

      IF (p_is_io) THEN

         CALL mpi_reduce (numelm_atm, nelm_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_io, p_err)
         IF (p_iam_io == p_root) THEN
            write(*,'(A,I12,A)') 'Total   : ', nelm_glb, ' ATM elements.'
         ENDIF

      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#endif

   END SUBROUTINE mesh_atm_build

   ! --------------------------------
   SUBROUTINE mesh_atm_free_mem (numelm_atm, mesh_atm, nelm_blk_atm, mesh_atm_pio, elmid_atm)

   USE MOD_Block
   USE MOD_DataType

   IMPLICIT NONE
   integer, intent(in) :: numelm_atm

   type (irregular_elm_type_atm),allocatable :: mesh_atm (:)  ! only available on IO processes

   integer, allocatable :: nelm_blk_atm(:,:)
   integer, allocatable :: mesh_atm_pio(:)

   type(block_data_int32_2d) :: elmid_atm

   ! Local variables
   integer :: ie
   integer :: iblkme, iblk, jblk

      IF (allocated(mesh_atm)) THEN
         DO ie = 1, numelm_atm
            deallocate (mesh_atm(ie)%ilon)
            deallocate (mesh_atm(ie)%ilat)
         ENDDO

         deallocate (mesh_atm)
      ENDIF
     
      IF (allocated(nelm_blk_atm)) deallocate(nelm_blk_atm)
      IF (allocated(mesh_atm_pio)) deallocate(mesh_atm_pio)

      IF (allocated (elmid_atm%blk)) then
         DO iblkme = 1, gblock%nblkme
            iblk = gblock%xblkme(iblkme)
            jblk = gblock%yblkme(iblkme)
            IF (allocated (elmid_atm%blk(iblk,jblk)%val)) deallocate (elmid_atm%blk(iblk,jblk)%val)
         ENDDO
         deallocate(elmid_atm%blk)
      ENDIF

   END SUBROUTINE mesh_atm_free_mem

END MODULE MOD_Mesh_atm

#endif
