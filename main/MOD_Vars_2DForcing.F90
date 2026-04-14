#include <define.h>

MODULE MOD_Vars_2DForcing
!-----------------------------------------------------------------------
!  Meteorogical Forcing
!
!  Created by Yongjiu Dai, 03/2014
!-----------------------------------------------------------------------
#ifndef CCPL 

   USE MOD_DataType
   IMPLICIT NONE
   SAVE

!-----------------------------------------------------------------------
   type(block_data_real8_2d) :: forc_xy_pco2m  ! CO2 concentration in atmos. (pascals)
   type(block_data_real8_2d) :: forc_xy_po2m   ! O2 concentration in atmos. (pascals)
   type(block_data_real8_2d) :: forc_xy_us     ! wind in eastward direction [m/s]
   type(block_data_real8_2d) :: forc_xy_vs     ! wind in northward direction [m/s]
   type(block_data_real8_2d) :: forc_xy_t      ! temperature at reference height [kelvin]
   type(block_data_real8_2d) :: forc_xy_q      ! specific humidity at reference height [kg/kg]
   type(block_data_real8_2d) :: forc_xy_prc    ! convective precipitation [mm/s]
   type(block_data_real8_2d) :: forc_xy_prl    ! large scale precipitation [mm/s]
   type(block_data_real8_2d) :: forc_xy_psrf   ! atmospheric pressure at the surface [pa]
   type(block_data_real8_2d) :: forc_xy_pbot   ! atm bottom level pressure (or reference height) (pa)
   type(block_data_real8_2d) :: forc_xy_sols   ! atm vis direct beam solar rad onto srf [W/m2]
   type(block_data_real8_2d) :: forc_xy_soll   ! atm nir direct beam solar rad onto srf [W/m2]
   type(block_data_real8_2d) :: forc_xy_solsd  ! atm vis diffuse solar rad onto srf [W/m2]
   type(block_data_real8_2d) :: forc_xy_solld  ! atm nir diffuse solar rad onto srf [W/m2]
   type(block_data_real8_2d) :: forc_xy_frl    ! atmospheric infrared (longwave) radiation [W/m2]
   type(block_data_real8_2d) :: forc_xy_hgt_u  ! observational height of wind [m]
   type(block_data_real8_2d) :: forc_xy_hgt_t  ! observational height of temperature [m]
   type(block_data_real8_2d) :: forc_xy_hgt_q  ! observational height of humidity [m]
   type(block_data_real8_2d) :: forc_xy_rhoair ! air density [kg/m3]
   type(block_data_real8_2d) :: forc_xy_hpbl   ! atmospheric boundary layer height [m]

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_2D_Forcing

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_2D_Forcing (grid)
   ! -------------------------------------------------------------------
   ! Allocates memory for CoLM 2d [lon_points,lat_points] variables
   ! -------------------------------------------------------------------
   USE MOD_SPMD_Task
   USE MOD_Grid
   USE MOD_DataType
   IMPLICIT NONE

   type(grid_type), intent(in) :: grid

      IF (p_is_io) THEN

         CALL allocate_block_data (grid, forc_xy_pco2m ) ! CO2 concentration in atmos. (pascals)
         CALL allocate_block_data (grid, forc_xy_po2m  ) ! O2 concentration in atmos. (pascals)
         CALL allocate_block_data (grid, forc_xy_us    ) ! wind in eastward direction [m/s]
         CALL allocate_block_data (grid, forc_xy_vs    ) ! wind in northward direction [m/s]
         CALL allocate_block_data (grid, forc_xy_t     ) ! temperature at reference height [kelvin]
         CALL allocate_block_data (grid, forc_xy_q     ) ! specific humidity at reference height [kg/kg]
         CALL allocate_block_data (grid, forc_xy_prc   ) ! convective precipitation [mm/s]
         CALL allocate_block_data (grid, forc_xy_prl   ) ! large scale precipitation [mm/s]
         CALL allocate_block_data (grid, forc_xy_psrf  ) ! atmospheric pressure at the surface [pa]
         CALL allocate_block_data (grid, forc_xy_pbot  ) ! atm bottom level pressure (or reference height) (pa)
         CALL allocate_block_data (grid, forc_xy_sols  ) ! atm vis direct beam solar rad onto srf [W/m2]
         CALL allocate_block_data (grid, forc_xy_soll  ) ! atm nir direct beam solar rad onto srf [W/m2]
         CALL allocate_block_data (grid, forc_xy_solsd ) ! atm vis diffuse solar rad onto srf [W/m2]
         CALL allocate_block_data (grid, forc_xy_solld ) ! atm nir diffuse solar rad onto srf [W/m2]
         CALL allocate_block_data (grid, forc_xy_frl   ) ! atmospheric infrared (longwave) radiation [W/m2]
         CALL allocate_block_data (grid, forc_xy_hgt_u ) ! observational height of wind [m]
         CALL allocate_block_data (grid, forc_xy_hgt_t ) ! observational height of temperature [m]
         CALL allocate_block_data (grid, forc_xy_hgt_q ) ! observational height of humidity [m]
         CALL allocate_block_data (grid, forc_xy_rhoair) ! air density [kg/m3]
         CALL allocate_block_data (grid, forc_xy_hpbl  ) ! atmospheric boundary layer height [m]
      ENDIF

   END SUBROUTINE allocate_2D_Forcing

#else

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

!-----------------------------------------------------------------------
   real(r8), allocatable :: forc_xy_pco2m(:)
   real(r8), allocatable :: forc_xy_po2m(:)
   real(r8), allocatable :: forc_xy_us(:)
   real(r8), allocatable :: forc_xy_vs(:)
   real(r8), allocatable :: forc_xy_t(:)
   real(r8), allocatable :: forc_xy_q(:)
   real(r8), allocatable :: forc_xy_prc(:)
   real(r8), allocatable :: forc_xy_prl(:)
   real(r8), allocatable :: forc_xy_psrf(:)
   real(r8), allocatable :: forc_xy_pbot(:)
   real(r8), allocatable :: forc_xy_sols(:)
   real(r8), allocatable :: forc_xy_soll(:)
   real(r8), allocatable :: forc_xy_solsd(:)
   real(r8), allocatable :: forc_xy_solld(:)
   real(r8), allocatable :: forc_xy_solarin(:)
   real(r8), allocatable :: forc_xy_frl(:)
   real(r8), allocatable :: forc_xy_hgt_u(:)
   real(r8), allocatable :: forc_xy_hgt_t(:)
   real(r8), allocatable :: forc_xy_hgt_q(:)
   real(r8), allocatable :: forc_xy_hpbl(:)
   real(r8), allocatable :: forc_xy_srflag(:)
   real(r8), allocatable :: topo_grid(:)

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_atm_Forcing

CONTAINS

!----------------------------------------------------------------------

   SUBROUTINE allocate_atm_Forcing (numelm_atm)
   ! -------------------------------------------------------------------
   ! Allocates memory for atmospheric forcing variables from atm models.
   ! Note that numelm_atm = 0 for non-IO proc, but still need to execute "allocate" statement due to CCPL use.
   ! -------------------------------------------------------------------

   use MOD_SPMD_Task
   use MOD_Vars_Global
   IMPLICIT NONE

   integer, intent(in) :: numelm_atm

         allocate (forc_xy_pco2m  (numelm_atm))
         allocate (forc_xy_po2m   (numelm_atm))
         allocate (forc_xy_us     (numelm_atm))
         allocate (forc_xy_vs     (numelm_atm))
         allocate (forc_xy_t      (numelm_atm))
         allocate (forc_xy_q      (numelm_atm))
         allocate (forc_xy_prc    (numelm_atm))
         allocate (forc_xy_prl    (numelm_atm))
         allocate (forc_xy_psrf   (numelm_atm))
         allocate (forc_xy_pbot   (numelm_atm))
         allocate (forc_xy_sols   (numelm_atm))
         allocate (forc_xy_soll   (numelm_atm))
         allocate (forc_xy_solsd  (numelm_atm))
         allocate (forc_xy_solld  (numelm_atm))
         allocate (forc_xy_solarin(numelm_atm))
         allocate (forc_xy_frl    (numelm_atm))
         allocate (forc_xy_hgt_u  (numelm_atm))
         allocate (forc_xy_hgt_t  (numelm_atm))
         allocate (forc_xy_hgt_q  (numelm_atm))
         allocate (forc_xy_hpbl   (numelm_atm))
         allocate (forc_xy_srflag (numelm_atm))
         allocate (topo_grid      (numelm_atm))

      IF (p_is_io) THEN 

         forc_xy_pco2m = spval
         forc_xy_po2m = spval
         forc_xy_us = spval
         forc_xy_vs = spval
         forc_xy_t = spval
         forc_xy_q = spval
         forc_xy_prc = spval
         forc_xy_prl = spval
         forc_xy_psrf = spval
         forc_xy_pbot = spval
         forc_xy_sols = spval
         forc_xy_soll = spval
         forc_xy_solsd = spval
         forc_xy_solld = spval
         forc_xy_solarin = spval
         forc_xy_frl = spval
         forc_xy_hgt_u = spval
         forc_xy_hgt_t = spval
         forc_xy_hgt_q = spval
         forc_xy_hpbl = spval
         forc_xy_srflag = spval
         topo_grid = spval

      END IF 

   END SUBROUTINE allocate_atm_Forcing

#endif

END MODULE MOD_Vars_2DForcing
! ---------- EOP ------------
