#include "cppdefs.h"
      MODULE mod_grid
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Hz         Thicknesses (m) of vertical RHO-points.                  !
!  Huon       Total U-momentum flux term, Hz*u/pn.                     !
!  Hvom       Total V-momentum flux term, Hz*v/pm.                     !
!  Tcline     Width (m) of surface or bottom boundary layer where      !
!               higher vertical resolution is required during          !
!               stretching.                                            !
!  angler     Angle (radians) between XI-axis and true EAST at         !
!               RHO-points.                                            !
!  dmde       ETA-derivative of inverse metric factor pm,              !
!               d(1/pm)/d(ETA).                                        !
!  dndx       XI-derivative  of inverse metric factor pn,              !
!               d(1/pn)/d(XI).                                         !
!  f          Coriolis parameter (1/s).                                !
!  fomn       Compound term, f/(pm*pn) at RHO points.                  !
!  grdscl     Grid scale used to adjust horizontal mixing according    !
!               to grid area.                                          !
!  h          Bottom depth (m) at RHO-points.                          !
!  hinv       Inverse of Bottom depth (1/m) at RHO-points.             !
!  latr       Latitude (degrees_north) at RHO-points.                  !
!  lonr       Longitude (degrees_east) at RHO-points.                  !
!  omm        RHO-grid area (meters2).                                 !
!  om_p       PSI-grid spacing (meters) in the XI-direction.           !
!  om_r       RHO-grid spacing (meters) in the XI-direction.           !
!  om_u       U-grid spacing (meters) in the XI-direction.             !
!  om_v       V-grid spacing (meters) in the XI-direction.             !
!  on_p       PSI-grid spacing (meters) in the ETA-direction.          !
!  on_r       RHO-grid spacing (meters) in the ETA-direction.          !
!  on_u       U-grid spacing (meters) in the ETA-direction.            !
!  on_v       V-grid spacing (meters) in the ETA-direction.            !
!  pm         Coordinate transformation metric "m" (1/meters)          !
!               associated with the differential distances in XI.      !
!  pmask      Slipperiness mask at PSI-points:                         !
!               (0=Land, 1=Sea, 1-gamma2=boundary).                    !
!  pmon_p     Compound term, pm/pn at PSI-points.                      !
!  pmon_r     Compound term, pm/pn at RHO-points.                      !
!  pmon_u     Compound term, pm/pn at U-points.                        !
!  pmon_v     Compound term, pm/pn at V-points.                        !
!  pn         Coordinate transformation metric "n" (1/meters)          !
!               associated with the differential distances in ETA.     !
!  pnom_p     Compound term, pn/pm at PSI-points.                      !
!  pnom_r     Compound term, pn/pm at RHO-points.                      !
!  pnom_u     Compound term, pn/pm at U-points.                        !
!  pnom_v     Compound term, pn/pm at V-points.                        !
!  rmask      Mask at RHO-points (0=Land, 1=Sea).                      !
!  umask      Mask at U-points (0=Land, 1=Sea).                        !
!  vmask      Mask at V-points (0=Land, 1=Sea).                        !
!  xp         XI-coordinates (m) at PSI-points.                        !
!  xr         XI-coordinates (m) at RHO-points.                        !
!  yp         ETA-coordinates (m) at PSI-points.                       !
!  yr         ETA-coordinates (m) at RHO-points.                       !
!  zice       Depth of ice shelf cavity (m, negative) at               !
!               RHO-points.                                            !
!  z_r        Actual depths (m) at horizontal RHO-points and           !
!               vertical RHO-points.                                   !
!  z_w        Actual depths (m) at horizontal RHO-points and           !
!               vertical W-points.                                     !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_GRID

          real(r8), pointer :: angler(:,:)
#if defined CURVGRID && defined UV_ADV
          real(r8), pointer :: dmde(:,:)
          real(r8), pointer :: dndx(:,:)
#endif
          real(r8), pointer :: f(:,:)
          real(r8), pointer :: fomn(:,:)
          real(r8), pointer :: grdscl(:,:)
          real(r8), pointer :: h(:,:)
          real(r8), pointer :: hinv(:,:)
          real(r8), pointer :: latr(:,:)
          real(r8), pointer :: lonr(:,:)
          real(r8), pointer :: omn(:,:)
          real(r8), pointer :: om_p(:,:)
          real(r8), pointer :: om_r(:,:)
          real(r8), pointer :: om_u(:,:)
          real(r8), pointer :: om_v(:,:)
          real(r8), pointer :: on_p(:,:)
          real(r8), pointer :: on_r(:,:)
          real(r8), pointer :: on_u(:,:)
          real(r8), pointer :: on_v(:,:)
          real(r8), pointer :: pm(:,:)
          real(r8), pointer :: pn(:,:)
          real(r8), pointer :: pmon_p(:,:)
          real(r8), pointer :: pmon_r(:,:)
          real(r8), pointer :: pmon_u(:,:)
          real(r8), pointer :: pmon_v(:,:)
          real(r8), pointer :: pnom_p(:,:)
          real(r8), pointer :: pnom_r(:,:)
          real(r8), pointer :: pnom_u(:,:)
          real(r8), pointer :: pnom_v(:,:)
          real(r8), pointer :: xp(:,:)
          real(r8), pointer :: xr(:,:)
          real(r8), pointer :: yp(:,:)
          real(r8), pointer :: yr(:,:)
#ifdef SOLVE3D
          real(r8), pointer :: Hz(:,:,:)
          real(r8), pointer :: Huon(:,:,:)
          real(r8), pointer :: Hvom(:,:,:)
          real(r8), pointer :: z_r(:,:,:)
          real(r8), pointer :: z_w(:,:,:)
# ifdef ICESHELF
          real(r8), pointer :: zice(:,:)
# endif
#endif
#ifdef MASKING
          real(r8), pointer :: pmask(:,:)
          real(r8), pointer :: rmask(:,:)
          real(r8), pointer :: umask(:,:)
          real(r8), pointer :: vmask(:,:)
#endif

        END TYPE T_GRID

        TYPE (T_GRID), allocatable :: GRID(:)

      CONTAINS

      SUBROUTINE allocate_grid (ng, LBi, UBi, LBj, UBj)
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
!
!  Local variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( GRID(Ngrids) )
!
      allocate ( GRID(ng) % angler(LBi:UBi,LBj:UBj) )

#if defined CURVGRID && defined UV_ADV
      allocate ( GRID(ng) % dmde(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % dndx(LBi:UBi,LBj:UBj) )
#endif

      allocate ( GRID(ng) % f(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % fomn(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % grdscl(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % h(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % hinv(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % latr(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % lonr(LBi:UBi,LBj:UBj))
      allocate ( GRID(ng) % omn(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % om_p(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % om_r(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % om_u(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % om_v(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % on_p(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % on_r(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % on_u(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % on_v(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pm(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pn(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pmon_p(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pmon_r(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pmon_u(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pmon_v(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pnom_p(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pnom_r(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pnom_u(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % pnom_v(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % xp(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % xr(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % yp(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % yr(LBi:UBi,LBj:UBj) )

#ifdef SOLVE3D

      allocate ( GRID(ng) % Hz(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % Huon(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % Hvom(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % z_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % z_w(LBi:UBi,LBj:UBj,0:N(ng)) )

# ifdef ICESHELF
      allocate ( GRID(ng) % zice(LBi:UBi,LBj:UBj) )
# endif

#endif

#ifdef MASKING
      allocate ( GRID(ng) % pmask(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % rmask(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % umask(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % vmask(LBi:UBi,LBj:UBj) )
#endif

      RETURN
      END SUBROUTINE allocate_grid

      SUBROUTINE initialize_grid (ng, tile)
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
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: IstrR, IendR, JstrR, JendR, IstrU, JstrV

      real(r8), parameter :: IniVal = 0.0_r8
!
#include "tile.h"
#ifdef DISTRIBUTE
      IstrR=LBi
      IendR=UBi
      JstrR=LBj
      JendR=UBj
#else
# include "set_bounds.h"
#endif
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
!
      GRID(ng) % angler(IstrR:IendR,JstrR:JendR) = IniVal
#if defined CURVGRID && defined UV_ADV
      GRID(ng) % dmde(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % dndx(IstrR:IendR,JstrR:JendR) = IniVal
#endif
      GRID(ng) % f(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % fomn(IstrR:IendR,JstrR:JendR) = IniVal

      GRID(ng) % grdscl(IstrR:IendR,JstrR:JendR) = IniVal

      GRID(ng) % h(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % hinv(IstrR:IendR,JstrR:JendR) = IniVal

      GRID(ng) % latr(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % lonr(IstrR:IendR,JstrR:JendR) = IniVal

      GRID(ng) % omn(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % om_p(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % om_r(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % om_u(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % om_v(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % on_p(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % on_r(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % on_u(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % on_v(IstrR:IendR,JstrR:JendR) = IniVal

      GRID(ng) % pm(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % pn(IstrR:IendR,JstrR:JendR) = IniVal

      GRID(ng) % pmon_p(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % pmon_r(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % pmon_u(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % pmon_v(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % pnom_p(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % pnom_r(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % pnom_u(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % pnom_v(IstrR:IendR,JstrR:JendR) = IniVal

      GRID(ng) % xp(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % xr(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % yp(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % yr(IstrR:IendR,JstrR:JendR) = IniVal

#ifdef SOLVE3D
      GRID(ng) % Hz(IstrR:IendR,JstrR:JendR,1:N(ng)) = IniVal

      GRID(ng) % Huon(IstrR:IendR,JstrR:JendR,1:N(ng)) = IniVal
      GRID(ng) % Hvom(IstrR:IendR,JstrR:JendR,1:N(ng)) = IniVal

      GRID(ng) % z_r(IstrR:IendR,JstrR:JendR,1:N(ng)) = IniVal
      GRID(ng) % z_w(IstrR:IendR,JstrR:JendR,0:N(ng)) = IniVal

# ifdef ICESHELF
      GRID(ng) % zice(IstrR:IendR,JstrR:JendR) = IniVal
# endif

#endif

#ifdef MASKING
      GRID(ng) % pmask(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % rmask(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % umask(IstrR:IendR,JstrR:JendR) = IniVal
      GRID(ng) % vmask(IstrR:IendR,JstrR:JendR) = IniVal
#endif

      RETURN
      END SUBROUTINE initialize_grid

      END MODULE mod_grid
