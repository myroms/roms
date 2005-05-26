#include "cppdefs.h"
      MODULE mod_grid
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  IJwaterR   Water points IJ couter for RHO-points masked variables.  !
!  IJwaterU   Water points IJ couter for   U-points masked variables.  !
!  IJwaterV   Water points IJ couter for   V-points masked variables.  !
!  Hz         Thicknesses (m) of vertical RHO-points.                  !
!  Huon       Total U-momentum flux term, Hz*u/pn.                     !
!  Hvom       Total V-momentum flux term, Hz*v/pm.                     !
!  IcePress   Pressure under the ice shelf at RHO-points.              !
!  Tcline     Width (m) of surface or bottom boundary layer where      !
!               higher vertical resolution is required during          !
!               stretching.                                            !
!  angler     Angle (radians) between XI-axis and true EAST at         !
!               RHO-points.                                            !
!  bed_thick0 Initial sum all bed layer thicknesses (m)                !
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
!  latp       Latitude (degrees_north) at PSI-points.                  !
!  latr       Latitude (degrees_north) at RHO-points.                  !
!  latu       Latitude (degrees_north) at U-points.                    !
!  latv       Latitude (degrees_north) at V-points.                    !
!  lonp       Longitude (degrees_east) at PSI-points.                  !
!  lonr       Longitude (degrees_east) at RHO-points.                  !
!  lonu       Longitude (degrees_east) at U-points.                    !
!  lonv       Longitude (degrees_east) at V-points.                    !
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
!  xu         XI-coordinates (m) at U-points.                          !
!  xv         XI-coordinates (m) at V-points.                          !
!  yp         ETA-coordinates (m) at PSI-points.                       !
!  yr         ETA-coordinates (m) at RHO-points.                       !
!  yu         ETA-coordinates (m) at U-points.                         !
!  yv         ETA-coordinates (m) at V-points.                         !
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
!
!  Nonlinear model state.
!
# if defined MASKING && defined PROPAGATOR
          integer, pointer :: IJwaterR(:,:)
          integer, pointer :: IJwaterU(:,:)
          integer, pointer :: IJwaterV(:,:)
# endif

          real(r8), pointer :: angler(:,:)
#ifdef SEDIMENT
          real(r8), pointer :: bed_thick0(:,:)
#endif
#if defined CURVGRID && defined UV_ADV
          real(r8), pointer :: dmde(:,:)
          real(r8), pointer :: dndx(:,:)
#endif
          real(r8), pointer :: f(:,:)
          real(r8), pointer :: fomn(:,:)
          real(r8), pointer :: grdscl(:,:)
          real(r8), pointer :: h(:,:)
          real(r8), pointer :: hinv(:,:)
          real(r8), pointer :: latp(:,:)
          real(r8), pointer :: latr(:,:)
          real(r8), pointer :: latu(:,:)
          real(r8), pointer :: latv(:,:)
          real(r8), pointer :: lonp(:,:)
          real(r8), pointer :: lonr(:,:)
          real(r8), pointer :: lonu(:,:)
          real(r8), pointer :: lonv(:,:)
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
          real(r8), pointer :: xu(:,:)
          real(r8), pointer :: xv(:,:)
          real(r8), pointer :: yp(:,:)
          real(r8), pointer :: yr(:,:)
          real(r8), pointer :: yu(:,:)
          real(r8), pointer :: yv(:,:)
#ifdef SOLVE3D
          real(r8), pointer :: Hz(:,:,:)
          real(r8), pointer :: Huon(:,:,:)
          real(r8), pointer :: Hvom(:,:,:)
          real(r8), pointer :: z_r(:,:,:)
          real(r8), pointer :: z_w(:,:,:)
# ifdef ICESHELF
          real(r8), pointer :: IcePress(:,:)
          real(r8), pointer :: zice(:,:)
# endif
#endif
#ifdef MASKING
          real(r8), pointer :: pmask(:,:)
          real(r8), pointer :: rmask(:,:)
          real(r8), pointer :: umask(:,:)
          real(r8), pointer :: vmask(:,:)
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
          real(r8), pointer :: tl_h(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: tl_Hz(:,:,:)
          real(r8), pointer :: tl_Huon(:,:,:)
          real(r8), pointer :: tl_Hvom(:,:,:)
          real(r8), pointer :: tl_z_r(:,:,:)
          real(r8), pointer :: tl_z_w(:,:,:)
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
          real(r8), pointer :: ad_h(:,:)
# ifdef SOLVE3D
          real(r8), pointer :: ad_Hz(:,:,:)
          real(r8), pointer :: ad_Huon(:,:,:)
          real(r8), pointer :: ad_Hvom(:,:,:)
          real(r8), pointer :: ad_z_r(:,:,:)
          real(r8), pointer :: ad_z_w(:,:,:)
# endif
#endif
        END TYPE T_GRID

        TYPE (T_GRID), allocatable :: GRID(:)

      CONTAINS

      SUBROUTINE allocate_grid (ng, LBi, UBi, LBj, UBj)
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
!  Nonlinear model state.
!
#if defined MASKING && defined PROPAGATOR
      allocate ( GRID(ng) % IJwaterR(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % IJwaterU(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % IJwaterV(LBi:UBi,LBj:UBj) )
#endif

      allocate ( GRID(ng) % angler(LBi:UBi,LBj:UBj) )
#ifdef SEDIMENT
      allocate ( GRID(ng) % bed_thick0(LBi:UBi,LBj:UBj) )
#endif
#if defined CURVGRID && defined UV_ADV
      allocate ( GRID(ng) % dmde(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % dndx(LBi:UBi,LBj:UBj) )
#endif

      allocate ( GRID(ng) % f(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % fomn(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % grdscl(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % h(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % hinv(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % latp(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % latr(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % latu(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % latv(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % lonp(LBi:UBi,LBj:UBj))
      allocate ( GRID(ng) % lonr(LBi:UBi,LBj:UBj))
      allocate ( GRID(ng) % lonu(LBi:UBi,LBj:UBj))
      allocate ( GRID(ng) % lonv(LBi:UBi,LBj:UBj))
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
      allocate ( GRID(ng) % xu(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % xv(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % yp(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % yr(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % yu(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % yv(LBi:UBi,LBj:UBj) )

#ifdef SOLVE3D

      allocate ( GRID(ng) % Hz(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % Huon(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % Hvom(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % z_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % z_w(LBi:UBi,LBj:UBj,0:N(ng)) )

# ifdef ICESHELF
      allocate ( GRID(ng) % IcePress(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % zice(LBi:UBi,LBj:UBj) )
# endif

#endif

#ifdef MASKING
      allocate ( GRID(ng) % pmask(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % rmask(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % umask(LBi:UBi,LBj:UBj) )
      allocate ( GRID(ng) % vmask(LBi:UBi,LBj:UBj) )
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      allocate ( GRID(ng) % tl_h(LBi:UBi,LBj:UBj) )
# ifdef SOLVE3D
      allocate ( GRID(ng) % tl_Hz(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % tl_Huon(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % tl_Hvom(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % tl_z_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % tl_z_w(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      allocate ( GRID(ng) % ad_h(LBi:UBi,LBj:UBj) )
# ifdef SOLVE3D
      allocate ( GRID(ng) % ad_Hz(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % ad_Huon(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % ad_Hvom(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % ad_z_r(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( GRID(ng) % ad_z_w(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif
#endif

      RETURN
      END SUBROUTINE allocate_grid

      SUBROUTINE initialize_grid (ng, tile)
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
!  Nonlinear model state.
!
#if defined MASKING && defined PROPAGATOR
      GRID(ng) % IJwaterR(Imin:Imax,Jmin:Jmax) = 0
      GRID(ng) % IJwaterU(Imin:Imax,Jmin:Jmax) = 0
      GRID(ng) % IJwaterV(Imin:Imax,Jmin:Jmax) = 0
#endif

      GRID(ng) % angler(Imin:Imax,Jmin:Jmax) = IniVal
#ifdef SEDIMENT
      GRID(ng) % bed_thick0(Imin:Imax,Jmin:Jmax) = IniVal
#endif
#if defined CURVGRID && defined UV_ADV
      GRID(ng) % dmde(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % dndx(Imin:Imax,Jmin:Jmax) = IniVal
#endif
      GRID(ng) % f(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % fomn(Imin:Imax,Jmin:Jmax) = IniVal

      GRID(ng) % grdscl(Imin:Imax,Jmin:Jmax) = IniVal

      GRID(ng) % h(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % hinv(Imin:Imax,Jmin:Jmax) = IniVal

      GRID(ng) % latp(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % latr(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % latu(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % latv(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % lonp(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % lonr(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % lonu(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % lonv(Imin:Imax,Jmin:Jmax) = IniVal

      GRID(ng) % omn(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % om_p(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % om_r(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % om_u(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % om_v(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % on_p(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % on_r(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % on_u(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % on_v(Imin:Imax,Jmin:Jmax) = IniVal

      GRID(ng) % pm(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % pn(Imin:Imax,Jmin:Jmax) = IniVal

      GRID(ng) % pmon_p(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % pmon_r(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % pmon_u(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % pmon_v(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % pnom_p(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % pnom_r(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % pnom_u(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % pnom_v(Imin:Imax,Jmin:Jmax) = IniVal

      GRID(ng) % xp(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % xr(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % xu(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % xv(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % yp(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % yu(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % yv(Imin:Imax,Jmin:Jmax) = IniVal

#ifdef SOLVE3D
      GRID(ng) % Hz(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal

      GRID(ng) % Huon(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      GRID(ng) % Hvom(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal

      GRID(ng) % z_r(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      GRID(ng) % z_w(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal

# ifdef ICESHELF
      GRID(ng) % IcePress(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % zice(Imin:Imax,Jmin:Jmax) = IniVal
# endif

#endif

#ifdef MASKING
      GRID(ng) % pmask(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % rmask(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % umask(Imin:Imax,Jmin:Jmax) = IniVal
      GRID(ng) % vmask(Imin:Imax,Jmin:Jmax) = IniVal
#endif

#if defined TANGENT || defined TL_IOMS
!
!  Tangent linear model state.
!
      GRID(ng) % tl_h(Imin:Imax,Jmin:Jmax) = IniVal
# ifdef SOLVE3D
      GRID(ng) % tl_Hz(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal

      GRID(ng) % tl_Huon(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      GRID(ng) % tl_Hvom(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal

      GRID(ng) % tl_z_r(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      GRID(ng) % tl_z_w(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
# endif
#endif

#ifdef ADJOINT
!
!  Adjoint model state.
!
      GRID(ng) % ad_h(Imin:Imax,Jmin:Jmax) = IniVal
# ifdef SOLVE3D
      GRID(ng) % ad_Hz(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal

      GRID(ng) % ad_Huon(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      GRID(ng) % ad_Hvom(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal

      GRID(ng) % ad_z_r(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      GRID(ng) % ad_z_w(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
# endif
#endif

      RETURN
      END SUBROUTINE initialize_grid

      END MODULE mod_grid
