#include "cppdefs.h"
      MODULE mod_clima
#ifdef CLIMATOLOGY
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Sea surface height climatology.                                     !
!                                                                      !
!   ssh         Climatology for sea surface height (m).                !
!   sshG        Latest two-time snapshots of input "ssh" grided        !
!                 data used for interpolation.                         !
!                                                                      !
!  2D momentum climatology.                                            !
!                                                                      !
!   ubarclm     Vertically integrated U-momentum climatology (m/s).    !
!   ubarclmG    Latest two-time snapshots of input "ubarclm" grided    !
!                 data used for interpolation.                         !
!   vbarclm     Vertically integrated V-momentum climatology (m/s).    !
!   vbarclmG    Latest two-time snapshots of input "vbarclm" grided    !
!                 data used for interpolation.                         !
!                                                                      !
!  Tracer climatology.                                                 !
!                                                                      !
!   tclm        Climatology for tracer type variables (usually,        !
!                 temperature: degC; salinity: PSU).                   !
!   tclmG       Latest two-time snapshots of input "tclm" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  3D momentum climatology.                                            !
!                                                                      !
!   uclm        3D U-momentum climatology (m/s).                       !
!   uclmG       Latest two-time snapshots of input "uclm" grided       !
!                 data used for interpolation.                         !
!   vclm        3D V-momentum climatology (m/s).                       !
!   vclmG       Latest two-time snapshots of input "vclm" grided       !
!                 data used for interpolation.                         !
!                                                                      !
!  Nudging variables.                                                  !
!                                                                      !
!   M2nudgcof   Time-scale (1/sec) coefficients for nudging towards    !
!                 2D momentum data.                                    !
!   M3nudgcof   Time-scale (1/sec) coefficients for nudging towards    !
!                 3D momentum data.                                    !
!   Tnudgcof    Time-scale (1/sec) coefficients for nudging towards    !
!                 tracer data.                                         !
!   Znudgcof    Time-scale (1/sec) coefficients for nudging towards    !
!                 sea surface height data.                             !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_CLIMA

# ifdef ZCLIMATOLOGY
          real(r8), pointer :: ssh(:,:)
#  ifndef ANA_SSH
          real(r8), pointer :: sshG(:,:,:)
#  endif
# endif
# ifdef ZCLM_NUDGING
          real(r8), pointer :: Znudgcof(:,:)
# endif
# ifdef M2CLIMATOLOGY
          real(r8), pointer :: ubarclm(:,:)
          real(r8), pointer :: vbarclm(:,:)
#  ifndef ANA_M2CLIMA
          real(r8), pointer :: ubarclmG(:,:,:)
          real(r8), pointer :: vbarclmG(:,:,:)
#  endif
# endif
# ifdef M2CLM_NUDGING
          real(r8), pointer :: M2nudgcof(:,:)
# endif
# ifdef SOLVE3D
#  ifdef TCLIMATOLOGY
          real(r8), pointer :: tclm(:,:,:,:)
#   ifndef ANA_TCLIMA
          real(r8), pointer :: tclmG(:,:,:,:,:)
#   endif
#  endif
#  ifdef TCLM_NUDGING
          real(r8), pointer :: Tnudgcof(:,:,:)
#  endif
#  ifdef M3CLIMATOLOGY
          real(r8), pointer :: uclm(:,:,:)
          real(r8), pointer :: vclm(:,:,:)
#   ifndef ANA_M3CLIMA
          real(r8), pointer :: uclmG(:,:,:,:)
          real(r8), pointer :: vclmG(:,:,:,:)
#   endif
#  endif
#  ifdef M3CLM_NUDGING
          real(r8), pointer :: M3nudgcof(:,:)
#  endif
# endif

        END TYPE T_CLIMA

        TYPE (T_CLIMA), allocatable :: CLIMA(:)

      CONTAINS

      SUBROUTINE allocate_clima (ng, LBi, UBi, LBj, UBj)
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
      IF (ng.eq.1) allocate ( CLIMA(Ngrids) )
!
# ifdef ZCLIMATOLOGY
      allocate ( CLIMA(ng) % ssh(LBi:UBi,LBj:UBj) )
#  ifndef ANA_SSH
      allocate ( CLIMA(ng) % sshG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef ZCLM_NUDGING
      allocate ( CLIMA(ng) % Znudgcof(LBi:UBi,LBj:UBj) )
# endif

# ifdef M2CLIMATOLOGY
      allocate ( CLIMA(ng) % ubarclm(LBi:UBi,LBj:UBj) )
      allocate ( CLIMA(ng) % vbarclm(LBi:UBi,LBj:UBj) )
#  ifndef ANA_M2CLIMA
      allocate ( CLIMA(ng) % ubarclmG(LBi:UBi,LBj:UBj,2) )
      allocate ( CLIMA(ng) % vbarclmG(LBi:UBi,LBj:UBj,2) )
#  endif
# endif

# ifdef M2CLM_NUDGING
      allocate ( CLIMA(ng) % M2nudgcof(LBi:UBi,LBj:UBj) )
# endif

# ifdef SOLVE3D
#  ifdef TCLIMATOLOGY
      allocate ( CLIMA(ng) % tclm(LBi:UBi,LBj:UBj,N(ng),NT(ng)) )
#   ifndef ANA_TCLIMA
      allocate ( CLIMA(ng) % tclmG(LBi:UBi,LBj:UBj,N(ng),2,NT(ng)) )
#   endif
#  endif

#  ifdef TCLM_NUDGING
      allocate ( CLIMA(ng) % Tnudgcof(LBi:UBi,LBj:UBj,NT(ng)) )
#  endif

#  ifdef M3CLIMATOLOGY
      allocate ( CLIMA(ng) % uclm(LBi:UBi,LBj:UBj,N(ng)) )
      allocate ( CLIMA(ng) % vclm(LBi:UBi,LBj:UBj,N(ng)) )
#   ifndef ANA_M3CLIMA
      allocate ( CLIMA(ng) % uclmG(LBi:UBi,LBj:UBj,N(ng),2) )
      allocate ( CLIMA(ng) % vclmG(LBi:UBi,LBj:UBj,N(ng),2) )
#   endif
#  endif

#  ifdef M3CLM_NUDGING
      allocate ( CLIMA(ng) % M3nudgcof(LBi:UBi,LBj:UBj) )
#  endif
# endif

      RETURN
      END SUBROUTINE allocate_clima

      SUBROUTINE initialize_clima (ng, tile)
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
      integer :: Imin, Imax, Jmin, Jmax

      real(r8), parameter :: IniVal = 0.0_r8
!
# include "tile.h"
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
# ifdef ZCLIMATOLOGY
      CLIMA(ng) % ssh(Imin:Imax,Jmin:Jmax) = IniVal
#  ifndef ANA_SSH
      CLIMA(ng) % sshG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif
# endif
# ifdef ZCLM_NUDGING
      CLIMA(ng) % Znudgcof(Imin:Imax,Jmin:Jmax) = IniVal
# endif

# ifdef M2CLIMATOLOGY
      CLIMA(ng) % ubarclm(Imin:Imax,Jmin:Jmax) = IniVal
      CLIMA(ng) % vbarclm(Imin:Imax,Jmin:Jmax) = IniVal
#  ifndef ANA_M2CLIMA
      CLIMA(ng) % ubarclmG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
      CLIMA(ng) % vbarclmG(Imin:Imax,Jmin:Jmax,1:2) = IniVal
#  endif
# endif
# ifdef M2CLM_NUDGING
      CLIMA(ng) % M2nudgcof(Imin:Imax,Jmin:Jmax) = IniVal
# endif

# ifdef SOLVE3D

#  ifdef M3CLIMATOLOGY
      CLIMA(ng) % uclm(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
      CLIMA(ng) % vclm(Imin:Imax,Jmin:Jmax,1:N(ng)) = IniVal
#   ifndef ANA_M3CLIMA
      CLIMA(ng) % uclmG(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
      CLIMA(ng) % vclmG(Imin:Imax,Jmin:Jmax,1:N(ng),1:2) = IniVal
#   endif
#  endif
#  ifdef M3CLM_NUDGING
      CLIMA(ng) % M3nudgcof(Imin:Imax,Jmin:Jmax) = IniVal
#  endif

#  ifdef TCLIMATOLOGY
      CLIMA(ng) % tclm(Imin:Imax,Jmin:Jmax,1:N(ng),1:NT(ng)) = IniVal
#   ifndef ANA_TCLIMA
      CLIMA(ng) % tclmG(Imin:Imax,Jmin:Jmax,1:N(ng),1:2,1:NT(ng)) = IniVal
#   endif
#  endif
#  ifdef TCLM_NUDGING
      CLIMA(ng) % Tnudgcof(Imin:Imax,Jmin:Jmax,1:NT(ng)) = IniVal
#  endif

# endif

      RETURN
      END SUBROUTINE initialize_clima
#endif
      END MODULE mod_clima
