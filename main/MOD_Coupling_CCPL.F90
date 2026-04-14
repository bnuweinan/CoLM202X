#include <define.h>

#ifdef CCPL
module MOD_Coupling_CCPL

   use MOD_Precision
   use MOD_Forcing
   use CCPL_interface_mod

   implicit none

   integer, public                 :: CoLM_comp_id
   integer, private                :: decomp_id, grid_h2d_id
   real(r8), allocatable, private  :: dlon(:), dlat(:)   ! lon & lat of each atm element in degrees on each IO proc
   integer,  allocatable, private  :: mask(:)            ! landmask of each atm element on each IO proc based on CoLM landtype data (1=land, 0=ocean)
   integer,  allocatable, private  :: eindex_proc(:)

contains

   subroutine register_component_coupling_configuration

   use MOD_SPMD_Task
   use MOD_Namelist
   use MOD_Vars_2DForcing
   use MOD_Vars_lnd2atm
   use MOD_NetCDFSerial

   implicit none

   integer          :: export_interface_id, import_interface_id
   integer          :: fields_id_import(18), fields_id_export(23),timer_id
   integer          :: ie,ipxl,numelm_atm_glob,coupling_freq
   real(r8),allocatable :: lat_atm(:,:,:),lon_atm(:,:,:),lat_glob(:),lon_glob(:)
   integer :: cpl_year,cpl_mon,cpl_day,cpl_sec
   integer :: s_year, s_month, s_day, s_seconds

      !valid datetime between CoLM namelist and coupler config
      s_year    = DEF_simulation_time%start_year
      s_month   = DEF_simulation_time%start_month
      s_day     = DEF_simulation_time%start_day
      s_seconds = DEF_simulation_time%start_sec

      call CCPL_get_start_time(CoLM_comp_id,cpl_year,cpl_mon,cpl_day,cpl_sec)
      if (cpl_year /= s_year .or. cpl_mon /= s_month .or. &
          cpl_day /= s_day .or. cpl_sec /= s_seconds)then
         if (p_is_master) then
            write(0,*) 'Error: datetime from coupler env_run.nml not match CoLM namelist'
            write(0,*) 'CoLM NAMELIST DATETIME:',s_year,s_month,s_day,s_seconds
            write(0,*) 'COUPLER DATETIME:',cpl_year,cpl_mon,cpl_day,cpl_sec
         endif

         call CCPL_abort("error in validing CoLM datetime from env_run.nml")

      endif

      if (p_is_io) then
         allocate (dlon(numelm_atm))
         allocate (dlat(numelm_atm))
         allocate (mask(numelm_atm))
         allocate (eindex_proc(numelm_atm))

         call ncio_read_serial(DEF_file_mesh_atm, 'lat_atm', lat_atm)
         call ncio_read_serial(DEF_file_mesh_atm, 'lon_atm', lon_atm)

         allocate(lat_glob(size(lat_atm)))
         allocate(lon_glob(size(lon_atm)))
         lat_glob = pack(lat_atm,.true.)
         lon_glob = pack(lon_atm,.true.)

         do ie = 1, numelm_atm
            dlat(ie)        = lat_glob(mesh_atm(ie)%indx)
            dlon(ie)        = lon_glob(mesh_atm(ie)%indx)
            mask(ie)        = mesh_atm(ie)%landmask
            eindex_proc(ie) = mesh_atm(ie)%indx
         end do

         deallocate (lat_atm)
         deallocate (lon_atm)
         deallocate (lat_glob)
         deallocate (lon_glob)

      else
         numelm_atm = 0
         allocate (dlon(numelm_atm))
         allocate (dlat(numelm_atm))
         allocate (mask(numelm_atm))
         allocate (eindex_proc(numelm_atm))
      end if

      call mpi_allreduce (numelm_atm, numelm_atm_glob, 1, MPI_INTEGER, MPI_SUM, p_comm_glb, p_err)

      call CCPL_set_normal_time_step(CoLM_comp_id, nint(DEF_simulation_time%timestep), annotation="setting the time step for CoLM")

      grid_h2d_id = CCPL_register_H2D_grid_via_local_data(CoLM_comp_id, "CoLM_Forcing_H2D_grid", "LON_LAT", "degrees",  &
                    "cyclic", numelm_atm_glob, numelm_atm, eindex_proc, 0.0, 360., -90.0, 90.0, dlon, dlat, &!mask, &
                    annotation="register CoLM Forcing H2D grid")

      decomp_id = CCPL_register_normal_parallel_decomp("decomp_CoLM_Forcing_grid", grid_h2d_id, numelm_atm, eindex_proc,&
                  annotation="register decomp for CoLM Forcing grid")

      !------------register CoLM field instances to C-Coupler3--------------
      call allocate_atm_Forcing (numelm_atm)
      call allocate_vars_lnd2atm (numelm_atm)

      fields_id_import(1) = CCPL_register_field_instance(forc_xy_us, "a2l_us", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="m/s", annotation="register u speed")
      fields_id_import(2) = CCPL_register_field_instance(forc_xy_vs, "a2l_vs", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="m/s", annotation="register v speed")
      fields_id_import(3) = CCPL_register_field_instance(forc_xy_t,  "a2l_t",  decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="K", annotation="register temperature at reference height")
      fields_id_import(4) = CCPL_register_field_instance(forc_xy_q,  "a2l_q",  decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="kg/kg", annotation="register specific humidity at reference height")
      fields_id_import(5) = CCPL_register_field_instance(forc_xy_prl,"a2l_prl",decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="mm/s", annotation="register large scale precipitation")
      fields_id_import(6) = CCPL_register_field_instance(forc_xy_prc,"a2l_prc",decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="mm/s", annotation="register convective precipitation")
      fields_id_import(7) = CCPL_register_field_instance(forc_xy_psrf,"a2l_psrf", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="Pa", annotation="register atmospheric pressure at the surface")
      fields_id_import(8) = CCPL_register_field_instance(forc_xy_pbot,"a2l_pbot", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="Pa", annotation="register atm model bottom layer pressure")
      fields_id_import(9) = CCPL_register_field_instance(forc_xy_solarin, "a2l_solarin", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="W/m2", annotation="register atm total solar rad onto srf")
      fields_id_import(10) = CCPL_register_field_instance(forc_xy_sols, "a2l_sols", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="W/m2", annotation="register atm vis direct beam solar rad onto srf")
      fields_id_import(11) = CCPL_register_field_instance(forc_xy_soll, "a2l_soll", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="W/m2", annotation="register atm nir direct beam solar rad onto srf")
      fields_id_import(12) = CCPL_register_field_instance(forc_xy_solsd, "a2l_solsd", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="W/m2", annotation="register atm vis diffuse solar rad onto srf")
      fields_id_import(13) = CCPL_register_field_instance(forc_xy_solld, "a2l_solld", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="W/m2", annotation="register atm nir diffuse solar rad onto srf")
      fields_id_import(14) = CCPL_register_field_instance(forc_xy_frl, "a2l_frl", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="W/m2", annotation="register atmospheric infrared (longwave) radiation")
      fields_id_import(15) = CCPL_register_field_instance(forc_xy_hgt_u, "a2l_hgt", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="m", annotation="register observational height")
      fields_id_import(16) = CCPL_register_field_instance(forc_xy_hpbl, "a2l_hpbl", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="m", annotation="register atmospheric boundary layer height")
      fields_id_import(17) = CCPL_register_field_instance(forc_xy_srflag, "a2l_srflag", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="0-1", annotation="register snow fraction in precipitation")
      fields_id_import(18) = CCPL_register_field_instance(topo_grid, "a2l_topography", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="m", annotation="register topography from atmospheric model")

      fields_id_export(1) = CCPL_register_field_instance(l2a_z0m, "l2a_z0m", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="m", annotation="register effective roughness")
      fields_id_export(2) = CCPL_register_field_instance(l2a_ustar, "l2a_ustar", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="m/s", annotation="register u* in similarity theory")
      fields_id_export(3) = CCPL_register_field_instance(l2a_fm, "l2a_fm", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="-", annotation="register integral of profile FUNCTION for momentum")
      fields_id_export(4) = CCPL_register_field_instance(l2a_fh, "l2a_fh", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="-", annotation="register integral of profile FUNCTION for heat")
      fields_id_export(5) = CCPL_register_field_instance(l2a_fq, "l2a_fq", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="-", annotation="register integral of profile FUNCTION for moisture")
      fields_id_export(6) = CCPL_register_field_instance(l2a_rib, "l2a_rib", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="-", annotation="register bulk Richardson number in surface layer")
      fields_id_export(7) = CCPL_register_field_instance(l2a_tref, "l2a_tref", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="K", annotation="register 2 m height air temperature")
      fields_id_export(8) = CCPL_register_field_instance(l2a_qref, "l2a_qref", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="kg/kg", annotation="register 2 m height air specific humidity")
      fields_id_export(9) = CCPL_register_field_instance(l2a_us10m, "l2a_us10m", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="m/s", annotation="register 10 m u wind speed")
      fields_id_export(10) = CCPL_register_field_instance(l2a_vs10m, "l2a_vs10m", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="m/s", annotation="register 10 m v wind speed")
      fields_id_export(11) = CCPL_register_field_instance(l2a_scv, "l2a_scv", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="mm", annotation="register snow water equivalent")
      fields_id_export(12) = CCPL_register_field_instance(l2a_snowdp, "l2a_snowdp", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="m", annotation="snow depth")
      fields_id_export(13) = CCPL_register_field_instance(l2a_trad, "l2a_trad", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="K", annotation="register radiative temperature of surface")
      fields_id_export(14) = CCPL_register_field_instance(l2a_fsno, "l2a_fsno", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="0-1", annotation="register fraction of snow cover on ground")
      fields_id_export(15) = CCPL_register_field_instance(l2a_qg, "l2a_qg", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="kg/kg", annotation="register ground specific humidity")
      fields_id_export(16) = CCPL_register_field_instance(l2a_lfevpa, "l2a_lfevpa", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="adj. MCV", annotation="register latent heat flux from canopy height to atmosphere")
      fields_id_export(17) = CCPL_register_field_instance(l2a_fsena, "l2a_fsena", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="adj. MCV", annotation="sensible heat from canopy height to atmosphere")
      fields_id_export(18) = CCPL_register_field_instance(l2a_alb_sols, "l2a_alb_sols", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="0-1", annotation="register vis direct beam solar averaged albedo")
      fields_id_export(19) = CCPL_register_field_instance(l2a_alb_soll, "l2a_alb_soll", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="0-1", annotation="register nir direct beam solar averaged albedo")
      fields_id_export(20) = CCPL_register_field_instance(l2a_alb_solsd, "l2a_alb_solsd", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="0-1", annotation="register vis diffuse solar averaged albedo")
      fields_id_export(21) = CCPL_register_field_instance(l2a_alb_solld, "l2a_alb_solld", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="0-1", annotation="register nir diffuse solar averaged albedo")
      fields_id_export(22) = CCPL_register_field_instance(l2a_olrg, "l2a_olrg", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="W/m2", annotation="register outgoing long-wave radiation from ground+canopy")
      fields_id_export(23) = CCPL_register_field_instance(l2a_landmask, "l2a_landmask", decomp_id, grid_h2d_id, 0, usage_tag=CCPL_TAG_CPL_REST, &
                            field_unit="0-1", annotation="register CoLM landmask")

      !--------register coupling frequency to C-Coupler3-------------
      coupling_freq = nint(DEF_simulation_time%timestep)
      timer_id = CCPL_define_single_timer(CoLM_comp_id, "seconds", coupling_freq, 0, 0, annotation="define a single timer for CoLM")

      !--------register export & import interface to C-Coupler3------
      export_interface_id = CCPL_register_export_interface("send_data_to_atm", size(fields_id_export), fields_id_export, timer_id, &
                            annotation="register interface for sending data to atm")
      import_interface_id = CCPL_register_import_interface("receive_data_from_atm", size(fields_id_import), fields_id_import, timer_id, 0, &
                            annotation="register interface for receiving data from atm")

      call CCPL_end_coupling_configuration(CoLM_comp_id, annotation = "component CoLM ends configuration")

   end subroutine register_component_coupling_configuration

end module MOD_Coupling_CCPL
#endif
