#include "cppdefs.h"
      MODULE mod_ocean
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  2D Primitive Variables.                                             !
!                                                                      !
!  rubar      Right-hand-side of 2D U-momentum equation (m4/s2).       !
!  rvbar      Right-hand-side of 2D V-momentum equation (m4/s2).       !
!  rzeta      Right-hand-side of free surface equation (m3/s).         !
!  ubar       Vertically integrated U-momentum component (m/s).        !
!  vbar       Vertically integrated V-momentum component (m/s).        !
!  zeta       Free surface (m).                                        !
!                                                                      !
!  3D Primitive Variables.                                             !
!                                                                      !
!  pden       Potential Density anomaly (kg/m3).                       !
!  rho        Density anomaly (kg/m3).                                 !
!  ru         Right-hand-side of 3D U-momentum equation (m4/s2).       !
!  rv         Right hand side of 3D V-momentum equation (m4/s2).       !
!  t          Tracer type variables (active and passive).              !
!  u          3D U-momentum component (m/s).                           !
!  v          3D V-momentum component (m/s).                           !
!  W          S-coordinate (omega*Hz/mn) vertical velocity (m3/s).     !
!                                                                      !
!  Biology Variables.                                                  !
!                                                                      !
!  pH         Surface concentration of hydrogen ions.                  !
!                                                                      !
!  Sediment Variables.                                                 !
!                                                                      !
!  bed        Sediment properties in each bed layer:                   !
!               bed(:,:,:,ithck) => layer thickness                    !
!               bed(:,:,:,iaged) => layer age                          !
!               bed(:,:,:,iporo) => layer porosity                     !
!               bed(:,:,:,idiff) => layer bio-diffusivity              !
!  bedldu     Bed load transport in xi  direction at u points (kg/m/s) !
!  bedldv     Bed load transport in eta direction at v points (kg/m/s) !
!  bed_frac   Sediment fraction of each size class in each bed layer   !
!               (nondimensional: 0-1.0).  Sum of bed_frac = 1.0        !
!  bed_mass   Sediment mass of each size class in each bed layer       !
!               (kg/m2).                                               !
!  bottom     Exposed sediment layer properties:                       !
!               bottom(:,:,isd50) => mean grain diameter               !
!               bottom(:,:,idens) => mean grain density                !
!               bottom(:,:,iwsed) => mean settling velocity            !
!               bottom(:,:,itauc) => mean critical erosion stress      !
!               bottom(:,:,irlen) => ripple length                     !
!               bottom(:,:,irhgt) => ripple height                     !
!               bottom(:,:,ibwav) => bed wave excursion amplitude      !
!               bottom(:,:,izNik) => Nikuradse bottom roughness        !
!               bottom(:,:,izbio) => biological bottom roughness       !
!               bottom(:,:,izbfm) => bed form bottom roughness         !
!               bottom(:,:,izbld) => bed load bottom roughness         !
!               bottom(:,:,izapp) => apparent bottom roughness         !
!               bottom(:,:,izwbl) => wave bottom roughness             !
!               bottom(:,:,izdef) => default bottom roughness          !
!               bottom(:,:,iactv) => active layer thickness            !
!               bottom(:,:,ishgt) => saltation height                  !
!
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_OCEAN
!
!  Nonlinear model state.
!
          real(r8), pointer :: rubar(:,:,:)
          real(r8), pointer :: rvbar(:,:,:)
          real(r8), pointer :: rzeta(:,:,:)
          real(r8), pointer :: ubar(:,:,:)
          real(r8), pointer :: vbar(:,:,:)
          real(r8), pointer :: zeta(:,:,:)
#ifdef SOLVE3D
          real(r8), pointer :: pden(:,:,:)
          real(r8), pointer :: rho(:,:,:)
          real(r8), pointer :: ru(:,:,:,:)
          real(r8), pointer :: rv(:,:,:,:)
          real(r8), pointer :: t(:,:,:,:,:)
          real(r8), pointer :: u(:,:,:,:)
          real(r8), pointer :: v(:,:,:,:)
          real(r8), pointer :: W(:,:,:)
          real(r8), pointer :: wvel(:,:,:)
# if defined BIO_FASHAM && defined CARBON
          real(r8), pointer :: pH(:,:)
# endif
# if defined SEDIMENT
          real(r8), pointer :: bed(:,:,:,:)
          real(r8), pointer :: bed_frac(:,:,:,:)
          real(r8), pointer :: bed_mass(:,:,:,:,:)
# endif
# if defined SEDIMENT || defined BBL_MODEL
          real(r8), pointer :: bottom(:,:,:)
#  ifdef BEDLOAD
          real(r8), pointer :: bedldu(:,:,:)
          real(r8), pointer :: bedldv(:,:,:)
#  endif
# endif
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
          real(r8), pointer :: tl_rubar(:,:,:)
          real(r8), pointer :: tl_rvbar(:,:,:)
          real(r8), pointer :: tl_rzeta(:,:,:)
          real(r8), pointer :: tl_ubar(:,:,:)
          real(r8), pointer :: tl_vbar(:,:,:)
          real(r8), pointer :: tl_zeta(:,:,:)
# ifdef SOLVE3D
          real(r8), pointer :: tl_pden(:,:,:)
          real(r8), pointer :: tl_rho(:,:,:)
          real(r8), pointer :: tl_ru(:,:,:,:)
          real(r8), pointer :: tl_rv(:,:,:,:)
          real(r8), pointer :: tl_t(:,:,:,:,:)
          real(r8), pointer :: tl_u(:,:,:,:)
          real(r8), pointer :: tl_v(:,:,:,:)
          real(r8), pointer :: tl_W(:,:,:)
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
          real(r8), pointer :: ad_rubar(:,:,:)
          real(r8), pointer :: ad_rvbar(:,:,:)
          real(r8), pointer :: ad_rzeta(:,:,:)
          real(r8), pointer :: ad_ubar(:,:,:)
          real(r8), pointer :: ad_vbar(:,:,:)
          real(r8), pointer :: ad_zeta(:,:,:)
# ifdef SOLVE3D
          real(r8), pointer :: ad_pden(:,:,:)
          real(r8), pointer :: ad_rho(:,:,:)
          real(r8), pointer :: ad_ru(:,:,:,:)
          real(r8), pointer :: ad_rv(:,:,:,:)
          real(r8), pointer :: ad_t(:,:,:,:,:)
          real(r8), pointer :: ad_u(:,:,:,:)
          real(r8), pointer :: ad_v(:,:,:,:)
          real(r8), pointer :: ad_W(:,:,:)
# endif
#endif

#ifdef FOUR_DVAR
!
!  Variational data assimilation background error covariance and
!  descent conjugate vectors (directions).
!
          real(r8), pointer :: b_ubar(:,:)
          real(r8), pointer :: b_vbar(:,:)
          real(r8), pointer :: b_zeta(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: b_t(:,:,:,:)
          real(r8), pointer :: b_u(:,:,:)
          real(r8), pointer :: b_v(:,:,:)
# endif
          real(r8), pointer :: d_ubar(:,:)
          real(r8), pointer :: d_vbar(:,:)
          real(r8), pointer :: d_zeta(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: d_t(:,:,:,:)
          real(r8), pointer :: d_u(:,:,:)
          real(r8), pointer :: d_v(:,:,:)
# endif
#endif

#if defined FORWARD_READ && \
   (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
# ifdef FORWARD_RHS
          real(r8), pointer :: rubarG(:,:,:)
          real(r8), pointer :: rvbarG(:,:,:)
          real(r8), pointer :: rzetaG(:,:,:)
# endif
          real(r8), pointer :: ubarG(:,:,:)
          real(r8), pointer :: vbarG(:,:,:)
          real(r8), pointer :: zetaG(:,:,:)
# ifdef SOLVE3D
#  ifdef FORWARD_RHS
          real(r8), pointer :: ruG(:,:,:,:)
          real(r8), pointer :: rvG(:,:,:,:)
#  endif
          real(r8), pointer :: tG(:,:,:,:,:)
          real(r8), pointer :: uG(:,:,:,:)
          real(r8), pointer :: vG(:,:,:,:)
# endif
#endif

        END TYPE T_OCEAN

        TYPE (T_OCEAN), allocatable :: OCEAN(:)

      CONTAINS

      SUBROUTINE allocate_ocean (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine allocates all variables in the module for all nested   !
!  grids.                                                              !
!                                                                      !
!=======================================================================
!
      USE mod_param
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( OCEAN(Ngrids) )
!
!  Nonlinear model state.
!
      allocate ( OCEAN(ng) % rubar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % rvbar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % rzeta(LBi:UBi,LBj:UBj,2) )

      allocate ( OCEAN(ng) % ubar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % vbar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % zeta(LBi:UBi,LBj:UBj,3) )

#ifdef SOLVE3D
      allocate ( OCEAN(ng) % pden(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % rho(LBi:UBi,LBj:UBj,N(ng)) )

      allocate ( OCEAN(ng) % ru(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( OCEAN(ng) % rv(LBi:UBi,LBj:UBj,0:N(ng),2) )

      allocate ( OCEAN(ng) % t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng)) )
      allocate ( OCEAN(ng) % u(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % v(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % W(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( OCEAN(ng) % wvel(LBi:UBi,LBj:UBj,0:N(ng)) )

# if defined BIO_FASHAM && defined CARBON
      allocate ( OCEAN(ng) % pH(LBi:UBi,LBj:UBj) )
# endif

# if defined SEDIMENT
      allocate ( OCEAN(ng) % bed(LBi:UBi,LBj:UBj,Nbed,MBEDP) )
      allocate ( OCEAN(ng) % bed_frac(LBi:UBi,LBj:UBj,Nbed,NST) )
      allocate ( OCEAN(ng) % bed_mass(LBi:UBi,LBj:UBj,Nbed,2,NST) )
# endif
# if defined SEDIMENT || defined BBL_MODEL
      allocate ( OCEAN(ng) % bottom(LBi:UBi,LBj:UBj,MBOTP) )
#  ifdef BEDLOAD
      allocate ( OCEAN(ng) % bedldu(LBi:UBi,LBj:UBj,NST) )
      allocate ( OCEAN(ng) % bedldv(LBi:UBi,LBj:UBj,NST) )
#  endif
# endif
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      allocate ( OCEAN(ng) % tl_rubar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % tl_rvbar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % tl_rzeta(LBi:UBi,LBj:UBj,2) )

      allocate ( OCEAN(ng) % tl_ubar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % tl_vbar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % tl_zeta(LBi:UBi,LBj:UBj,3) )

# ifdef SOLVE3D
      allocate ( OCEAN(ng) % tl_pden(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % tl_rho(LBi:UBi,LBj:UBj,N(ng)) )

      allocate ( OCEAN(ng) % tl_ru(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( OCEAN(ng) % tl_rv(LBi:UBi,LBj:UBj,0:N(ng),2) )

      allocate ( OCEAN(ng) % tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng)) )
      allocate ( OCEAN(ng) % tl_u(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % tl_v(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % tl_W(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      allocate ( OCEAN(ng) % ad_rubar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % ad_rvbar(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % ad_rzeta(LBi:UBi,LBj:UBj,2) )

      allocate ( OCEAN(ng) % ad_ubar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % ad_vbar(LBi:UBi,LBj:UBj,3) )
      allocate ( OCEAN(ng) % ad_zeta(LBi:UBi,LBj:UBj,3) )

# ifdef SOLVE3D
      allocate ( OCEAN(ng) % ad_pden(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % ad_rho(LBi:UBi,LBj:UBj,N(ng)) )

      allocate ( OCEAN(ng) % ad_ru(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( OCEAN(ng) % ad_rv(LBi:UBi,LBj:UBj,0:N(ng),2) )

      allocate ( OCEAN(ng) % ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng)) )
      allocate ( OCEAN(ng) % ad_u(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % ad_v(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % ad_W(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif
#endif

#ifdef FOUR_DVAR
!
!  Variational data assimilation background error covariance and
!  descent conjugate vectors (directions).
!
      allocate ( OCEAN(ng) % b_ubar(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % b_vbar(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % b_zeta(LBi:UBi,LBj:UBj) )

# ifdef SOLVE3D
      allocate ( OCEAN(ng) % b_t(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      allocate ( OCEAN(ng) % b_u(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % b_v(LBi:UBi,LBj:UBj,N(ng)) )
# endif

      allocate ( OCEAN(ng) % d_ubar(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % d_vbar(LBi:UBi,LBj:UBj) )
      allocate ( OCEAN(ng) % d_zeta(LBi:UBi,LBj:UBj) )

# ifdef SOLVE3D
      allocate ( OCEAN(ng) % d_t(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
      allocate ( OCEAN(ng) % d_u(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( OCEAN(ng) % d_v(LBi:UBi,LBj:UBj,N(ng)) )
# endif
#endif

#if defined FORWARD_READ && \
   (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
# ifdef FORWARD_RHS
      allocate ( OCEAN(ng) % rubarG(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % rvbarG(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % rzetaG(LBi:UBi,LBj:UBj,2) )
# endif
      allocate ( OCEAN(ng) % ubarG(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % vbarG(LBi:UBi,LBj:UBj,2) )
      allocate ( OCEAN(ng) % zetaG(LBi:UBi,LBj:UBj,2) )

# ifdef SOLVE3D
#  ifdef FORWARD_RHS
      allocate ( OCEAN(ng) % ruG(LBi:UBi,LBj:UBj,0:N(ng),2) )
      allocate ( OCEAN(ng) % rvG(LBi:UBi,LBj:UBj,0:N(ng),2) )
#  endif
      allocate ( OCEAN(ng) % tG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng)) )
      allocate ( OCEAN(ng) % uG(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( OCEAN(ng) % vG(LBi:UBi,LBj:UBj,N(ng),2) )
# endif
#endif

      RETURN
      END SUBROUTINE allocate_ocean

      SUBROUTINE initialize_ocean (ng, tile, model)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax

      real(r8), parameter :: IniVal = 0.0_r8
!
#include "tile.h"
!
!  Set array initialization range.
!
#ifdef _OPENMP
      IF (WESTERN_EDGE) THEN
        Imin=LBi
      ELSE
        Imin=Istr
      END IF
      IF (EASTERN_EDGE) THEN
        Imax=UBi
      ELSE
        Imax=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
        Jmin=LBj
      ELSE
        Jmin=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
        Jmax=UBj
      ELSE
        Jmax=Jend
      END IF
#else
      Imin=LBi
      Imax=UBi
      Jmin=LBj
      Jmax=UBj
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear model state.
!
      IF ((model.eq.0).or.(model.eq.iNLM)) THEN
 
        OCEAN(ng) % rubar(Imin:Imax,Jmin:Jmax,1:2) = IniVal
        OCEAN(ng) % rvbar(Imin:Imax,Jmin:Jmax,1:2) = IniVal
        OCEAN(ng) % rzeta(Imin:Imax,Jmin:Jmax,1:2) = IniVal

        OCEAN(ng) % ubar(Imin:Imax,Jmin:Jmax,1:3) = IniVal
        OCEAN(ng) % vbar(Imin:Imax,Jmin:Jmax,1:3) = IniVal
        OCEAN(ng) % zeta(Imin:Imax,Jmin:Jmax,1:3) = IniVal

#ifdef SOLVE3D
        OCEAN(ng) % pden(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
        OCEAN(ng) % rho(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal

        OCEAN(ng) % ru(Imin:Imax,Jmin:Jmax,0:N(ng),1:2) = IniVal
        OCEAN(ng) % rv(Imin:Imax,Jmin:Jmax,0:N(ng),1:2) = IniVal

        OCEAN(ng) % t(Imin:Imax,Jmin:Jmax,1:N(ng),1:3,1:NT(ng)) = IniVal
        OCEAN(ng) % u(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
        OCEAN(ng) % v(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
 
        OCEAN(ng) % W(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
        OCEAN(ng) % wvel(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal

# if defined BIO_FASHAM && defined CARBON
        OCEAN(ng) % pH(Imin:Imax,Jmin:Jmax) = 8.0_r8
# endif

# ifdef SEDIMENT
        OCEAN(ng) % bed(Imin:Imax,Jmin:Jmax,1:Nbed,1:MBEDP) = IniVal
        OCEAN(ng) % bed_frac(Imin:Imax,Jmin:Jmax,1:Nbed,1:NST) = IniVal
        OCEAN(ng) % bed_mass(Imin:Imax,Jmin:Jmax,1:Nbed,1:2,1:NST) =    &
     &              IniVal
# endif

# if defined SEDIMENT || defined BBL_MODEL
        OCEAN(ng) % bottom(Imin:Imax,Jmin:Jmax,1:MBOTP) = IniVal
#  ifdef BEDLOAD
        OCEAN(ng) % bedldu(Imin:Imax,Jmin:Jmax,1:NST) = IniVal
        OCEAN(ng) % bedldv(Imin:Imax,Jmin:Jmax,1:NST) = IniVal
#  endif
# endif

#endif
      END IF

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      IF ((model.eq.0).or.(model.eq.iTLM).or.(model.eq.iRPM)) THEN

        OCEAN(ng) % tl_rubar(Imin:Imax,Jmin:Jmax,1:2) = IniVal
        OCEAN(ng) % tl_rvbar(Imin:Imax,Jmin:Jmax,1:2) = IniVal
        OCEAN(ng) % tl_rzeta(Imin:Imax,Jmin:Jmax,1:2) = IniVal

        OCEAN(ng) % tl_ubar(Imin:Imax,Jmin:Jmax,1:3) = IniVal
        OCEAN(ng) % tl_vbar(Imin:Imax,Jmin:Jmax,1:3) = IniVal
        OCEAN(ng) % tl_zeta(Imin:Imax,Jmin:Jmax,1:3) = IniVal

# ifdef SOLVE3D
        OCEAN(ng) % tl_pden(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
        OCEAN(ng) % tl_rho(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal

        OCEAN(ng) % tl_ru(Imin:Imax,Jmin:Jmax,0:N(ng),1:2) = IniVal
        OCEAN(ng) % tl_rv(Imin:Imax,Jmin:Jmax,0:N(ng),1:2) = IniVal

        OCEAN(ng) % tl_t(Imin:Imax,Jmin:Jmax,1:N(ng),1:3,1:NT(ng)) =    &
     &              IniVal

        OCEAN(ng) % tl_u(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
        OCEAN(ng) % tl_v(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
 
        OCEAN(ng) % tl_W(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
# endif
      END IF
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      IF ((model.eq.0).or.(model.eq.iADM)) THEN

        OCEAN(ng) % ad_rubar(Imin:Imax,Jmin:Jmax,1:2) = IniVal
        OCEAN(ng) % ad_rvbar(Imin:Imax,Jmin:Jmax,1:2) = IniVal
        OCEAN(ng) % ad_rzeta(Imin:Imax,Jmin:Jmax,1:2) = IniVal

        OCEAN(ng) % ad_ubar(Imin:Imax,Jmin:Jmax,1:3) = IniVal
        OCEAN(ng) % ad_vbar(Imin:Imax,Jmin:Jmax,1:3) = IniVal
        OCEAN(ng) % ad_zeta(Imin:Imax,Jmin:Jmax,1:3) = IniVal

# ifdef SOLVE3D
        OCEAN(ng) % ad_pden(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
        OCEAN(ng) % ad_rho(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal

        OCEAN(ng) % ad_ru(Imin:Imax,Jmin:Jmax,0:N(ng),1:2) = IniVal
        OCEAN(ng) % ad_rv(Imin:Imax,Jmin:Jmax,0:N(ng),1:2) = IniVal

        OCEAN(ng) % ad_t(Imin:Imax,Jmin:Jmax,1:N(ng),1:3,1:NT(ng)) =    &
     &              IniVal

        OCEAN(ng) % ad_u(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
        OCEAN(ng) % ad_v(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
 
        OCEAN(ng) % ad_W(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
# endif
      END IF
#endif

#ifdef FOUR_DVAR
!
!  Variational data assimilation background error covariance and
!  descent conjugate vectors (directions).
!
      OCEAN(ng) % b_ubar(Imin:Imax,Jmin:Jmax) = IniVal
      OCEAN(ng) % b_vbar(Imin:Imax,Jmin:Jmax) = IniVal
      OCEAN(ng) % b_zeta(Imin:Imax,Jmin:Jmax) = IniVal

# ifdef SOLVE3D
      OCEAN(ng) % b_t(Imin:Imax,Jmin:Jmax,1:N(ng),1:NT(ng)) = IniVal
      OCEAN(ng) % b_u(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      OCEAN(ng) % b_v(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
# endif

      OCEAN(ng) % d_ubar(Imin:Imax,Jmin:Jmax) = IniVal
      OCEAN(ng) % d_vbar(Imin:Imax,Jmin:Jmax) = IniVal
      OCEAN(ng) % d_zeta(Imin:Imax,Jmin:Jmax) = IniVal

# ifdef SOLVE3D
      OCEAN(ng) % d_t(Imin:Imax,Jmin:Jmax,1:N(ng),1:NT(ng)) = IniVal
      OCEAN(ng) % d_u(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      OCEAN(ng) % d_v(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
# endif
#endif

#if defined FORWARD_READ && \
   (defined TANGENT || defined TL_IOMS || defined ADJOINT)
!
!  Latest two records of the nonlinear trajectory used to interpolate
!  the background state in the tangent linear and adjoint models.
!
# ifdef FORWARD_RHS
      OCEAN(ng) % rubarG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      OCEAN(ng) % rvbarG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      OCEAN(ng) % rzetaG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
# endif
      OCEAN(ng) % ubarG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      OCEAN(ng) % vbarG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      OCEAN(ng) % zetaG(Imin:Imax,Jmin:Jmax,1:2) = IniVal

# ifdef SOLVE3D
#  ifdef FORWARD_RHS
      OCEAN(ng) % ruG(Imin:Imax,Jmin:Jmax,0:N(ng),1:2) = IniVal
      OCEAN(ng) % rvG(Imin:Imax,Jmin:Jmax,0:N(ng),1:2) = IniVal
#  endif
      OCEAN(ng) % tG(Imin:Imax,Jmin:Jmax,1:N(ng),1:2,1:NT(ng)) = IniVal
      OCEAN(ng) % uG(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
      OCEAN(ng) % vG(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
# endif
#endif

      RETURN
      END SUBROUTINE initialize_ocean

      END MODULE mod_ocean
