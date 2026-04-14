#include <define.h>

#ifdef CCPL

MODULE MOD_Vars_lnd2atm
!-----------------------------------------------------------------------
!  Variables exported from CoLM to atmospheric models.
!
!  Created by Nan Wei, 2025.02.
!-----------------------------------------------------------------------

   USE MOD_Precision

   IMPLICIT NONE
   SAVE

!-----------------------------------------------------------------------
   real(r8), allocatable :: l2a_z0m(:)      ! effective roughness [m]
   real(r8), allocatable :: l2a_ustar(:)    ! u* in similarity theory [m/s]
   real(r8), allocatable :: l2a_fm(:)       ! integral of profile FUNCTION for momentum
   real(r8), allocatable :: l2a_fh(:)       ! integral of profile FUNCTION for heat
   real(r8), allocatable :: l2a_fq(:)       ! integral of profile FUNCTION for moisture
   real(r8), allocatable :: l2a_rib(:)      ! bulk Richardson number in surface layer
   real(r8), allocatable :: l2a_tref(:)     ! 2 m height air temperature [kelvin]
   real(r8), allocatable :: l2a_qref(:)     ! 2 m height air specific humidity [kg/kg]
   real(r8), allocatable :: l2a_us10m(:)    ! 10 m u wind speed [m/s]
   real(r8), allocatable :: l2a_vs10m(:)    ! 10 m v wind speed [m/s]
   real(r8), allocatable :: l2a_scv(:)      ! snow cover, water equivalent [mm]
   real(r8), allocatable :: l2a_snowdp(:)   ! snow depth [meter]
   real(r8), allocatable :: l2a_trad(:)     ! radiative temperature of surface [K]
   real(r8), allocatable :: l2a_fsno(:)     ! fraction of snow cover on ground
   real(r8), allocatable :: l2a_qg(:)       ! ground specific humidity [kg/kg]
   real(r8), allocatable :: l2a_lfevpa(:)   ! latent heat flux from canopy height to atmosphere [adj. MCV]
   real(r8), allocatable :: l2a_fsena(:)    ! sensible heat from canopy height to atmosphere [adj. MCV]
   real(r8), allocatable :: l2a_alb_sols(:) ! vis direct beam solar averaged albedo [-]
   real(r8), allocatable :: l2a_alb_soll(:) ! nir direct beam solar averaged albedo [-]
   real(r8), allocatable :: l2a_alb_solsd(:)! vis diffuse solar averaged albedo [-]
   real(r8), allocatable :: l2a_alb_solld(:)! nir diffuse solar averaged albedo [-]
   real(r8), allocatable :: l2a_olrg(:)     ! outgoing long-wave radiation from ground+canopy [W/m2]

   real(r8), allocatable :: l2a_landmask(:) ! landmask in CoLM

!-------------not used in MCV-----------------
   real(r8), allocatable :: l2a_taux(:)
   real(r8), allocatable :: l2a_tauy(:)
   real(r8), allocatable :: l2a_fevpa(:)

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_vars_lnd2atm
   PUBLIC :: deallocate_vars_lnd2atm
   PUBLIC :: collect_data_from_patch_to_atmelm

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_vars_lnd2atm (numelm_atm)
   ! -------------------------------------------------------------------
   ! Allocates memory for land variables exported to atm models.
   ! Note that numelm_atm = 0 for non-IO proc, but still need to execute "allocate" statement due to CCPL use.
   ! -------------------------------------------------------------------

   IMPLICIT NONE

   integer, intent(in) :: numelm_atm

   !   IF (p_is_io) THEN

         allocate (l2a_z0m      (numelm_atm))
         allocate (l2a_ustar    (numelm_atm))
         allocate (l2a_fm       (numelm_atm))
         allocate (l2a_fh       (numelm_atm))
         allocate (l2a_fq       (numelm_atm))
         allocate (l2a_rib      (numelm_atm))
         allocate (l2a_tref     (numelm_atm))
         allocate (l2a_qref     (numelm_atm))
         allocate (l2a_us10m    (numelm_atm))
         allocate (l2a_vs10m    (numelm_atm))
         allocate (l2a_scv      (numelm_atm))
         allocate (l2a_snowdp   (numelm_atm))
         allocate (l2a_trad     (numelm_atm))
         allocate (l2a_fsno     (numelm_atm))
         allocate (l2a_qg       (numelm_atm))
         allocate (l2a_lfevpa   (numelm_atm))
         allocate (l2a_fsena    (numelm_atm))
         allocate (l2a_alb_sols (numelm_atm))
         allocate (l2a_alb_soll (numelm_atm))
         allocate (l2a_alb_solsd(numelm_atm))
         allocate (l2a_alb_solld(numelm_atm))
         allocate (l2a_olrg     (numelm_atm))

         allocate (l2a_taux     (numelm_atm))
         allocate (l2a_tauy     (numelm_atm))
         allocate (l2a_fevpa    (numelm_atm))

         allocate (l2a_landmask (numelm_atm))

   !   ENDIF

   END SUBROUTINE allocate_vars_lnd2atm

   SUBROUTINE deallocate_vars_lnd2atm ()
   ! -------------------------------------------------------------------
   ! Allocates memory for land variables exported to atm models.
   ! Note that numelm_atm = 0 for non-IO proc, but still need to execute "allocate" statement due to CCPL use.
   ! -------------------------------------------------------------------

   IMPLICIT NONE

   !   IF (p_is_io) THEN

         if (allocated(l2a_z0m      )) deallocate(l2a_z0m      )
         if (allocated(l2a_ustar    )) deallocate(l2a_ustar    )
         if (allocated(l2a_fm       )) deallocate(l2a_fm       )
         if (allocated(l2a_fh       )) deallocate(l2a_fh       )
         if (allocated(l2a_fq       )) deallocate(l2a_fq       )
         if (allocated(l2a_rib      )) deallocate(l2a_rib      )
         if (allocated(l2a_tref     )) deallocate(l2a_tref     )
         if (allocated(l2a_qref     )) deallocate(l2a_qref     )
         if (allocated(l2a_us10m    )) deallocate(l2a_us10m    )
         if (allocated(l2a_vs10m    )) deallocate(l2a_vs10m    )
         if (allocated(l2a_scv      )) deallocate(l2a_scv      )
         if (allocated(l2a_snowdp   )) deallocate(l2a_snowdp   )
         if (allocated(l2a_trad     )) deallocate(l2a_trad     )
         if (allocated(l2a_fsno     )) deallocate(l2a_fsno     )
         if (allocated(l2a_qg       )) deallocate(l2a_qg       )
         if (allocated(l2a_lfevpa   )) deallocate(l2a_lfevpa   )
         if (allocated(l2a_fsena    )) deallocate(l2a_fsena    )
         if (allocated(l2a_alb_sols )) deallocate(l2a_alb_sols )
         if (allocated(l2a_alb_soll )) deallocate(l2a_alb_soll )
         if (allocated(l2a_alb_solsd)) deallocate(l2a_alb_solsd)
         if (allocated(l2a_alb_solld)) deallocate(l2a_alb_solld)
         if (allocated(l2a_olrg     )) deallocate(l2a_olrg     )

         if (allocated(l2a_taux     )) deallocate(l2a_taux     )
         if (allocated(l2a_tauy     )) deallocate(l2a_tauy     )
         if (allocated(l2a_fevpa    )) deallocate(l2a_fevpa    )

         if (allocated(l2a_landmask )) deallocate(l2a_landmask )
   !   ENDIF

   END SUBROUTINE deallocate_vars_lnd2atm

   SUBROUTINE collect_data_from_patch_to_atmelm (mg2p_forc, numelm_atm, mesh_atm)
   
   use MOD_Mesh_atm
   use MOD_SpatialMapping_atm
   use MOD_Vars_TimeInvariants
   use MOD_Vars_TimeVariables
   use MOD_Vars_1DFluxes
   use MOD_Vars_1DForcing
   use MOD_Vars_2DForcing
   use MOD_SPMD_Task
   use MOD_Namelist
   use MOD_Const_Physical
   use MOD_Vars_Global
   use MOD_FrictionVelocity
   use MOD_TurbulenceLEddy

   IMPLICIT NONE

   type (spatial_mapping_type), intent(in) :: mg2p_forc
   integer, intent(in) :: numelm_atm 
   type (irregular_elm_type_atm), allocatable, intent(in) :: mesh_atm (:)

   ! Local Variables
   real(r8) rhoair,thm,th,thv,ur,displa_av,zldis,hgt_u,hgt_t,hgt_q
   real(r8) hpbl ! atmospheric boundary layer height [m]
   real(r8) z0m_av,z0h_av,z0q_av,us,vs,tm,qm,psrf,taux_e,tauy_e,fsena_e,fevpa_e
   real(r8) r_ustar_e, r_tstar_e, r_qstar_e, r_zol_e, r_ustar2_e, r_fm10m_e
   real(r8) r_fm_e, r_fh_e, r_fq_e, r_rib_e, r_us10m_e, r_vs10m_e
   real(r8) obu,fh2m,fq2m
   real(r8) um,thvstar,beta,zii,wc,wc2
   integer  i

      CALL mg2p_forc%pset2grid (z0m,         l2a_z0m      ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (tref,        l2a_tref     ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (qref,        l2a_qref     ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (scv,         l2a_scv      ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (snowdp,      l2a_snowdp   ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (trad,        l2a_trad     ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (fsno,        l2a_fsno     ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (fsena,       l2a_fsena    ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (lfevpa,      l2a_lfevpa   ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (alb(1,1,:),  l2a_alb_sols ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (alb(2,1,:),  l2a_alb_soll ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (alb(1,2,:),  l2a_alb_solsd,  msk = patchmask)
      CALL mg2p_forc%pset2grid (alb(2,2,:),  l2a_alb_solld,  msk = patchmask)
      CALL mg2p_forc%pset2grid (olrg,        l2a_olrg     ,  msk = patchmask)
      
      CALL mg2p_forc%pset2grid (taux,        l2a_taux     ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (tauy,        l2a_tauy     ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (fevpa,       l2a_fevpa    ,  msk = patchmask)

      if (p_is_io) then
         do i = 1, numelm_atm

            l2a_landmask(i) = mesh_atm(i)%landmask
            if (mesh_atm(i)%landmask == 1) then
               z0m_av  = l2a_z0m       (i)
               hgt_u   = forc_xy_hgt_u (i)
               hgt_t   = forc_xy_hgt_t (i)
               hgt_q   = forc_xy_hgt_q (i)
               us      = forc_xy_us    (i)
               vs      = forc_xy_vs    (i)
               tm      = forc_xy_t     (i)
               qm      = forc_xy_q     (i)
               psrf    = forc_xy_pbot  (i)
               taux_e  = l2a_taux      (i)
               tauy_e  = l2a_tauy      (i)
               fsena_e = l2a_fsena     (i)
               fevpa_e = l2a_fevpa     (i)
               IF (DEF_USE_CBL_HEIGHT) THEN !//TODO: Shaofeng, 2023.05.18
                  hpbl = forc_xy_hpbl  (i)
               ENDIF

               z0h_av = z0m_av
               z0q_av = z0m_av

               displa_av = 2./3.*z0m_av/0.07

               hgt_u = max(hgt_u, 5.+displa_av)
               hgt_t = max(hgt_t, 5.+displa_av)
               hgt_q = max(hgt_q, 5.+displa_av)

               zldis = hgt_u-displa_av

               rhoair = (psrf - 0.378*qm*psrf/(0.622+0.378*qm)) / (rgas*tm)

               r_ustar_e = sqrt(max(1.e-6,sqrt(taux_e**2+tauy_e**2))/rhoair)
               r_tstar_e = -fsena_e/(rhoair*r_ustar_e)/cpair
               r_qstar_e = -fevpa_e/(rhoair*r_ustar_e)

               thm = tm + 0.0098*hgt_t
               th  = tm*(100000./psrf)**(rgas/cpair)
               thv = th*(1.+0.61*qm)

               r_zol_e = zldis*vonkar*grav * (r_tstar_e*(1.+0.61*qm)+0.61*th*r_qstar_e) &
                  / (r_ustar_e**2*thv)

               IF(r_zol_e >= 0.)THEN   !stable
                  r_zol_e = min(2.,max(r_zol_e,1.e-6))
               ELSE                       !unstable
                  r_zol_e = max(-100.,min(r_zol_e,-1.e-6))
               ENDIF

               beta = 1.
               zii = 1000.

               thvstar=r_tstar_e*(1.+0.61*qm)+0.61*th*r_qstar_e
               ur = sqrt(us*us+vs*vs)
               IF(r_zol_e >= 0.)THEN
                  um = max(ur,0.1)
               ELSE
                  IF (DEF_USE_CBL_HEIGHT) THEN !//TODO: Shaofeng, 2023.05.18
                     zii = max(5.*hgt_u,hpbl)
                  ENDIF !//TODO: Shaofeng, 2023.05.18
                  wc = (-grav*r_ustar_e*thvstar*zii/thv)**(1./3.)
                  wc2 = beta*beta*(wc*wc)
                  um = max(0.1,sqrt(ur*ur+wc2))
               ENDIF

               obu = zldis/r_zol_e
               IF (DEF_USE_CBL_HEIGHT) THEN
                  CALL moninobuk_leddy(hgt_u,hgt_t,hgt_q,displa_av,z0m_av,z0h_av,z0q_av,&
                     obu,um, hpbl, r_ustar2_e,fh2m,fq2m,r_fm10m_e,r_fm_e,r_fh_e,r_fq_e) !Shaofeng, 2023.05.20
               ELSE
                  CALL moninobuk(hgt_u,hgt_t,hgt_q,displa_av,z0m_av,z0h_av,z0q_av,&
                    obu,um,r_ustar2_e,fh2m,fq2m,r_fm10m_e,r_fm_e,r_fh_e,r_fq_e) !Shaofeng, 2023.05.20
               ENDIF

               ! bug found by chen qiying 2013/07/01
               r_rib_e = r_zol_e /vonkar * r_ustar2_e**2 / (vonkar/r_fh_e*um**2)
               r_rib_e = min(5.,r_rib_e)

               r_us10m_e = us/um * r_ustar2_e /vonkar * r_fm10m_e
               r_vs10m_e = vs/um * r_ustar2_e /vonkar * r_fm10m_e

               ! Assign values to atm element.
               l2a_ustar (i) = r_ustar_e
               l2a_rib   (i) = r_rib_e
               l2a_fm    (i) = r_fm_e
               l2a_fh    (i) = r_fh_e
               l2a_fq    (i) = r_fq_e
               l2a_us10m (i) = r_us10m_e
               l2a_vs10m (i) = r_vs10m_e
               l2a_qg    (i) = fevpa_e/(rhoair*vonkar/r_fq_e*r_ustar2_e) + qm

            else

               l2a_ustar (i) = spval
               l2a_rib   (i) = spval
               l2a_fm    (i) = spval
               l2a_fh    (i) = spval
               l2a_fq    (i) = spval
               l2a_us10m (i) = spval
               l2a_vs10m (i) = spval
               l2a_qg    (i) = spval

            endif
         end do
      endif

      CALL mg2p_forc%pset2grid (fsena/forc_rhoair/cpair,  l2a_fsena    ,  msk = patchmask)
      CALL mg2p_forc%pset2grid (lfevpa/forc_rhoair/hvap,  l2a_lfevpa   ,  msk = patchmask)

   END SUBROUTINE collect_data_from_patch_to_atmelm

END MODULE MOD_Vars_lnd2atm
#endif
! ---------- EOP ------------
