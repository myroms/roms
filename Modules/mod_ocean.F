#include "cppdefs.h"
      MODULE mod_ocean
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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
!  t          Tracer type variables (usually, potential temperature    !
!               and salinity).                                         !
!  u          3D U-momentum component (m/s).                           !
!  v          3D V-momentum component (m/s).                           !
!  W          S-coordinate (omega*Hz/mn) vertical velocity (m3/s).     !
!                                                                      !
!  Sediment Variables.                                                 !
!                                                                      !
!  bed        Sediment properties in each bed layer:                   !
!               bed(:,:,:,ithck) => layer thickness                    !
!               bed(:,:,:,iaged) => layer age                          !
!               bed(:,:,:,iporo) => layer porosity                     !
!               bed(:,:,:,idiff) => layer bio-diffusivity              !
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
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_OCEAN
          
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
# if defined SEDIMENT
          real(r8), pointer :: bed(:,:,:,:)
          real(r8), pointer :: bed_frac(:,:,:,:)
          real(r8), pointer :: bed_mass(:,:,:,:)
# endif
# if defined SEDIMENT || defined BBL_MODEL
          real(r8), pointer :: bottom(:,:,:)
# endif
#endif

        END TYPE T_OCEAN

        TYPE (T_OCEAN), allocatable :: OCEAN(:)

      CONTAINS

      SUBROUTINE allocate_ocean (ng, LBi, UBi, LBj, UBj)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( OCEAN(Ngrids) )
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

# if defined SEDIMENT
      allocate ( OCEAN(ng) % bed(LBi:UBi,LBj:UBj,Nbed,MBEDP) )
      allocate ( OCEAN(ng) % bed_frac(LBi:UBi,LBj:UBj,Nbed,NST) )
      allocate ( OCEAN(ng) % bed_mass(LBi:UBi,LBj:UBj,Nbed,NST) )
# endif
# if defined SEDIMENT || defined BBL_MODEL
      allocate ( OCEAN(ng) % bottom(LBi:UBi,LBj:UBj,MBOTP) )
# endif
#endif

      RETURN
      END SUBROUTINE allocate_ocean

      SUBROUTINE initialize_ocean (ng, tile)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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
      integer, intent(in) :: ng, tile
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

      OCEAN(ng) % u(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
      OCEAN(ng) % v(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
 
      OCEAN(ng) % W(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
      OCEAN(ng) % wvel(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal

      OCEAN(ng) % t(Imin:Imax,Jmin:Jmax,1:N(ng),1:2,1:NT(ng)) = IniVal

# if defined SEDIMENT || defined BBL_MODEL
      OCEAN(ng) % bed(Imin:Imax,Jmin:Jmax,1:Nbed,1:MBEDP) = IniVal
      OCEAN(ng) % bed_frac(Imin:Imax,Jmin:Jmax,1:Nbed,1:NST) = IniVal
      OCEAN(ng) % bed_mass(Imin:Imax,Jmin:Jmax,1:Nbed,1:NST) = IniVal
# endif

# if defined SEDIMENT || defined BBL_MODEL
      OCEAN(ng) % bottom(Imin:Imax,Jmin:Jmax,1:MBOTP) = IniVal
# endif

#endif
      RETURN
      END SUBROUTINE initialize_ocean

      END MODULE mod_ocean
