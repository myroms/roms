#include "cppdefs.h"
      MODULE mod_mixing
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Horizontal and vertical mixing coefficients:                        !
!                                                                      !
!  Akt          Vertical mixing coefficient (m2/s) for tracers.        !
!  dAktdz       Vertical gradient in mixing coefficient (m/s) for      !
!                 tracer 1, used in float random walk calculations     !
!  Akv          Vertical mixing coefficient (m2/s) for momentum.       !
!                                                                      !
!  Variables associated with the equation of state:                    !
!                                                                      !
!  alpha        Surface thermal expansion coefficient (1/Celsius).     !
!  beta         Surface saline contraction coefficient (1/PSU).        !
!  bvf          Brunt-Vaisala frequency squared (1/s2).                !
!  neutral      Coefficient to convert "in situ" density to neutral    !
!                 surface.                                             !
!                                                                      !
!  tke          Turbulent energy squared (m2/s2) at horizontal         !
!                 at W-points.                                         !
!  gls          Turbulent energy squared times turbulent length        !
!                 scale (m3/s2) at W-points.                           !
!                                                                      !
!  Large/McWilliams/Doney interior vertical mixing variables:          !
!                                                                      !
!  alfaobeta    Ratio of thermal expansion and saline contraction      !
!                 coefficients (Celsius/PSU) used in double            !
!                 diffusion.                                           !
!                                                                      !
!  Water clarity parameters:                                           !
!                                                                      !
!  Jwtype       Water clarity (Jerlov water type classification).      !
!                                                                      !
!  Large/McWilliams/Doney oceanic boundary layer variables:            !
!                                                                      !
!  ghats        Boundary layer nonlocal transport (T units/m).         !
!  hbbl         Depth of bottom oceanic boundary layer (m).            !
!  hsbl         Depth of surface oceanic boundary layer (m).           !
!  kbbl         Index of grid level above bottom  boundary layer.      !
!  ksbl         Index of grid level below surface boundary layer.      !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_MIXING

#ifdef SOLVE3D
# if defined LMD_SKPP || defined SOLAR_SOURCE
          integer,  pointer :: Jwtype(:,:)
# endif
# if defined LMD_SKPP || defined LMD_BKPP
          integer,  pointer :: ksbl(:,:)
#  ifdef LMD_BKPP
          integer,  pointer :: kbbl(:,:)
#  endif
# endif
#endif

#if defined UV_VIS2 || !defined SOLVE3D
          real(r8), pointer :: visc2_p(:,:)
          real(r8), pointer :: visc2_r(:,:)
#endif
#ifdef UV_VIS4
          real(r8), pointer :: visc4_p(:,:)
          real(r8), pointer :: visc4_r(:,:)
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
          real(r8), pointer :: diff2(:,:,:)
# endif
# ifdef TS_DIF4
          real(r8), pointer :: diff4(:,:,:)
# endif
          real(r8), pointer :: Akv(:,:,:)
          real(r8), pointer :: Akt(:,:,:,:)
# ifdef FLOAT_VWALK
          real(r8), pointer :: dAktdz(:,:,:)
# endif
# if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
          real(r8), pointer :: alpha(:,:)
          real(r8), pointer :: beta(:,:)
# endif
# ifdef BV_FREQUENCY
          real(r8), pointer :: bvf(:,:,:)
# endif
# if defined MIX_ISO_TS || defined GENT_McWILLIAMS
          real(r8), pointer :: neutral(:,:,:)
# endif
# if defined GLS_MIXING || defined MY25_MIXING
          real(r8), pointer :: tke(:,:,:,:)
          real(r8), pointer :: gls(:,:,:,:)
          real(r8), pointer :: Lscale(:,:,:)
          real(r8), pointer :: Akk(:,:,:)
#  ifdef GLS_MIXING
          real(r8), pointer :: Akp(:,:,:)
#  endif
# endif
# if defined LMD_MIXING && defined LMD_DDMIX
          real(r8), pointer :: alfaobeta(:,:,:)
# endif
# if defined LMD_SKPP || defined LMD_BKPP
          real(r8), pointer :: hsbl(:,:)
#  ifdef LMD_BKPP
          real(r8), pointer :: hbbl(:,:)
#  endif
#  ifdef LMD_NONLOCAL
          real(r8), pointer :: ghats(:,:,:,:)
#  endif
# endif
#endif

        END TYPE T_MIXING

        TYPE (T_MIXING), allocatable :: MIXING(:)

      CONTAINS
        
      SUBROUTINE allocate_mixing (ng, LBi, UBi, LBj, UBj)
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
!  Allocate module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) allocate ( MIXING(Ngrids) )
! 
#if defined UV_VIS2 || !defined SOLVE3D
      allocate ( MIXING(ng) % visc2_p(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % visc2_r(LBi:UBi,LBj:UBj) )
#endif

#ifdef UV_VIS4
      allocate ( MIXING(ng) % visc4_p(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % visc4_r(LBi:UBi,LBj:UBj) )
#endif

#ifdef SOLVE3D
# ifdef TS_DIF2
      allocate ( MIXING(ng) % diff2(LBi:UBi,LBj:UBj,NT(ng)) )
# endif

# ifdef TS_DIF4
      allocate ( MIXING(ng) % diff4(LBi:UBi,LBj:UBj,NT(ng)) )
# endif

      allocate ( MIXING(ng) % Akv(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % Akt(LBi:UBi,LBj:UBj,0:N(ng),NAT) )

# ifdef FLOAT_VWALK
      allocate ( MIXING(ng) % dAktdz(LBi:UBi,LBj:UBj,N(ng)) )
# endif

# if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
      allocate ( MIXING(ng) % alpha(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % beta(LBi:UBi,LBj:UBj) )
# endif

# ifdef BV_FREQUENCY
      allocate ( MIXING(ng) % bvf(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif

# if defined MIX_ISO_TS || defined GENT_McWILLIAMS
      allocate ( MIXING(ng) % neutral(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif

# if defined GLS_MIXING || defined MY25_MIXING
      allocate ( MIXING(ng) % tke(LBi:UBi,LBj:UBj,0:N(ng),3) )
      allocate ( MIXING(ng) % gls(LBi:UBi,LBj:UBj,0:N(ng),3) )
      allocate ( MIXING(ng) % Lscale(LBi:UBi,LBj:UBj,0:N(ng)) )
      allocate ( MIXING(ng) % Akk(LBi:UBi,LBj:UBj,0:N(ng)) )

#  ifdef GLS_MIXING
      allocate ( MIXING(ng) % Akp(LBi:UBi,LBj:UBj,0:N(ng)) )
#  endif
# endif

# if defined LMD_MIXING && defined LMD_DDMIX
      allocate ( MIXING(ng) % alfaobeta(LBi:UBi,LBj:UBj,0:N(ng)) )
# endif

# if defined LMD_SKPP || defined SOLAR_SOURCE
      allocate ( MIXING(ng) % Jwtype(LBi:UBi,LBj:UBj) )
# endif

# if defined LMD_SKPP || defined LMD_BKPP
      allocate ( MIXING(ng) % ksbl(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % hsbl(LBi:UBi,LBj:UBj) )

#  ifdef LMD_BKPP
      allocate ( MIXING(ng) % kbbl(LBi:UBi,LBj:UBj) )
      allocate ( MIXING(ng) % hbbl(LBi:UBi,LBj:UBj) )
#  endif

#  ifdef LMD_NONLOCAL
      allocate ( MIXING(ng) % ghats(LBi:UBi,LBj:UBj,0:N(ng),NAT) )
#  endif

# endif
#endif

      RETURN
      END SUBROUTINE allocate_mixing

      SUBROUTINE initialize_mixing (ng, tile)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine allocates and initializes all variables in module      !
!  "mod_mixing" for all nested grids.                                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: itrc

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
#if defined UV_VIS2 || !defined SOLVE3D
      MIXING(ng) % visc2_p(Imin:Imax,Jmin:Jmax) = visc2(ng)
      MIXING(ng) % visc2_r(Imin:Imax,Jmin:Jmax) = visc2(ng)
#endif

#ifdef UV_VIS4
      MIXING(ng) % visc4_p(Imin:Imax,Jmin:Jmax) = visc4(ng)
      MIXING(ng) % visc4_r(Imin:Imax,Jmin:Jmax) = visc4(ng)
#endif

#ifdef SOLVE3D

# ifdef TS_DIF2
      DO itrc=1,NT(ng)
        MIXING(ng) % diff2(Imin:Imax,Jmin:Jmax,itrc) = tnu2(itrc,ng)
      END DO
# endif

# ifdef TS_DIF4
      DO itrc=1,NT(ng)
        MIXING(ng) % diff4(Imin:Imax,Jmin:Jmax,itrc) = tnu4(itrc,ng)
      END DO
# endif

      MIXING(ng) % Akv(Imin:Imax,Jmin:Jmax,0) = IniVal
      MIXING(ng) % Akv(Imin:Imax,Jmin:Jmax,1:N(ng)-1) = Akv_bak(ng)
      MIXING(ng) % Akv(Imin:Imax,Jmin:Jmax,N(ng)) = IniVal

      MIXING(ng) % Akt(Imin:Imax,Jmin:Jmax,0,1:NAT) = IniVal
      DO itrc=1,NAT
        MIXING(ng) % Akt(Imin:Imax,Jmin:Jmax,1:N(ng)-1,itrc) =          &
     &                  Akt_bak(itrc,ng)
      END DO
      MIXING(ng) % Akt(Imin:Imax,Jmin:Jmax,N(ng),1:NAT) = IniVal

# ifdef FLOAT_VWALK
      MIXING(ng) % dAktdz(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
# endif

# if defined LMD_SKPP || defined LMD_BKPP || defined BULK_FLUXES
      MIXING(ng) % alpha(Imin:Imax,Jmin:Jmax) = IniVal
      MIXING(ng) % beta(Imin:Imax,Jmin:Jmax) = IniVal
# endif

# ifdef BV_FREQUENCY
      MIXING(ng) % bvf(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
# endif

# if defined MIX_ISO_TS || defined GENT_McWILLIAMS
      MIXING(ng) % neutral(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
# endif

# if defined GLS_MIXING || defined MY25_MIXING
      MIXING(ng) % tke(Imin:Imax,Jmin:Jmax,0:N(ng),1:3) = gls_Kmin(ng)
      MIXING(ng) % gls(Imin:Imax,Jmin:Jmax,0:N(ng),1:3) = gls_Pmin(ng)

      MIXING(ng) % Lscale(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal

      MIXING(ng) % Akk(Imin:Imax,Jmin:Jmax,0) = IniVal
      MIXING(ng) % Akk(Imin:Imax,Jmin:Jmax,1:N(ng)-1) = Akk_bak(ng)
      MIXING(ng) % Akk(Imin:Imax,Jmin:Jmax,N(ng)) = IniVal

#  ifdef GLS_MIXING
      MIXING(ng) % Akp(Imin:Imax,Jmin:Jmax,0) = IniVal
      MIXING(ng) % Akp(Imin:Imax,Jmin:Jmax,1:N(ng)-1) = Akk_bak(ng)
      MIXING(ng) % Akp(Imin:Imax,Jmin:Jmax,N(ng)) = IniVal
#  endif
# endif

# if defined LMD_MIXING && defined LMD_DDMIX
      MIXING(ng) % alfaobeta(Imin:Imax,Jmin:Jmax,0:N(ng)) = IniVal
# endif

# if defined LMD_SKPP || defined SOLAR_SOURCE
      MIXING(ng) % Jwtype(Imin:Imax,Jmin:Jmax) = lmd_Jwt(ng)
# endif

# if defined LMD_SKPP || defined LMD_BKPP
      MIXING(ng) % ksbl(Imin:Imax,Jmin:Jmax) = 0
      MIXING(ng) % hsbl(Imin:Imax,Jmin:Jmax) = IniVal
#  ifdef LMD_BKPP
      MIXING(ng) % kbbl(Imin:Imax,Jmin:Jmax) = 0
      MIXING(ng) % hbbl(Imin:Imax,Jmin:Jmax) = IniVal
#  endif

#  ifdef LMD_NONLOCAL
      MIXING(ng) % ghats(Imin:Imax,Jmin:Jmax,0:N(ng),1:NAT) = IniVal
#  endif
# endif

#endif

      RETURN
      END SUBROUTINE initialize_mixing

      END MODULE mod_mixing
