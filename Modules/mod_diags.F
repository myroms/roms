#include "cppdefs.h"
      MODULE mod_diags
#ifdef DIAGNOSTICS
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Diagnostics fields for output.                                      !
!                                                                      !
!  DiaBio2d  Diagnostics for 2D biological terms.                      !
!  DiaBio3d  Diagnostics for 3D biological terms.                      !
!  DiaTrc    Diagnostics for tracer terms.                             !
!  DiaU2d    Diagnostics for 2D U-momentum terms.                      !
!  DiaV2d    Diagnostics for 2D V-momentum terms.                      !
!  DiaU3d    Diagnostics for 3D U-momentum terms.                      !
!  DiaV3d    Diagnostics for 3D V-momentum terms.                      !
!                                                                      !
!  Diagnostics fields work arrays.                                     !
!                                                                      !
!  DiaTwrk    Diagnostics work array for tracer terms.                 !
!  DiaU2wrk   Diagnostics work array for 2D U-momentum terms.          !
!  DiaV2wrk   Diagnostics work array for 2D V-momentum terms.          !
!  DiaRUbar   Diagnostics RHS array for 2D U-momentum terms.           !
!  DiaRVbar   Diagnostics RHS array for 2D V-momentum terms.           !
!  DiaU2int   Diagnostics array for 2D U-momentum terms                !
!               integrated over short barotropic timesteps.            !
!  DiaV2int   Diagnostics array for 2D U-momentum terms                !
!               integrated over short barotropic timesteps.            !
!  DiaRUfrc   Diagnostics forcing array for 2D U-momentum terms.       !
!  DiaRVfrc   Diagnostics forcing array for 2D V-momentum terms.       !
!  DiaU3wrk   Diagnostics work array for 3D U-momentum terms.          !
!  DiaV3wrk   Diagnostics work array for 3D V-momentum terms.          !
!  DiaRU      Diagnostics RHS array for 3D U-momentum terms.           !
!  DiaRV      Diagnostics RHS array for 3D V-momentum terms.           !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        TYPE T_DIAGS

# ifdef BIO_FASHAM
          real(r8), pointer :: DiaBio2d(:,:,:)
          real(r8), pointer :: DiaBio3d(:,:,:,:)
# endif

# ifdef DIAGNOSTICS_TS
          real(r8), pointer :: DiaTrc(:,:,:,:,:)
          real(r8), pointer :: DiaTwrk(:,:,:,:,:)
# endif

# ifdef DIAGNOSTICS_UV
          real(r8), pointer :: DiaU2d(:,:,:)
          real(r8), pointer :: DiaV2d(:,:,:)
          real(r8), pointer :: DiaU2wrk(:,:,:)
          real(r8), pointer :: DiaV2wrk(:,:,:)
          real(r8), pointer :: DiaRUbar(:,:,:,:)
          real(r8), pointer :: DiaRVbar(:,:,:,:)

#  ifdef SOLVE3D
          real(r8), pointer :: DiaU2int(:,:,:)
          real(r8), pointer :: DiaV2int(:,:,:)
          real(r8), pointer :: DiaRUfrc(:,:,:,:)
          real(r8), pointer :: DiaRVfrc(:,:,:,:)
          real(r8), pointer :: DiaU3d(:,:,:,:)
          real(r8), pointer :: DiaV3d(:,:,:,:)
          real(r8), pointer :: DiaU3wrk(:,:,:,:)
          real(r8), pointer :: DiaV3wrk(:,:,:,:)
          real(r8), pointer :: DiaRU(:,:,:,:,:)
          real(r8), pointer :: DiaRV(:,:,:,:,:)
#  endif

# endif

        END TYPE T_DIAGS

        TYPE (T_DIAGS), allocatable :: DIAGS(:)

      CONTAINS

      SUBROUTINE allocate_diags (ng, LBi, UBi, LBj, UBj)
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
      IF (ng.eq.1 ) allocate ( DIAGS(Ngrids) )
!
# ifdef BIO_FASHAM
      allocate ( DIAGS(ng) % DiaBio2d(LBi:UBi,LBj:UBj,NDbio2d) )
      allocate ( DIAGS(ng) % DiaBio3d(LBi:UBi,LBj:UBj,N(ng),NDbio3d) )
# endif

# ifdef DIAGNOSTICS_TS
      allocate ( DIAGS(ng) % DiaTrc(LBi:UBi,LBj:UBj,N(ng),NT(ng),NDT) )
      allocate ( DIAGS(ng) % DiaTwrk(LBi:UBi,LBj:UBj,N(ng),NT(ng),NDT) )
# endif

# ifdef DIAGNOSTICS_UV
      allocate ( DIAGS(ng) % DiaU2d(LBi:UBi,LBj:UBj,NDM2d) )
      allocate ( DIAGS(ng) % DiaV2d(LBi:UBi,LBj:UBj,NDM2d) )
      allocate ( DIAGS(ng) % DiaU2wrk(LBi:UBi,LBj:UBj,NDM2d) )
      allocate ( DIAGS(ng) % DiaV2wrk(LBi:UBi,LBj:UBj,NDM2d) )
      allocate ( DIAGS(ng) % DiaRUbar(LBi:UBi,LBj:UBj,2,NDM2d-1) )
      allocate ( DIAGS(ng) % DiaRVbar(LBi:UBi,LBj:UBj,2,NDM2d-1) )

#  ifdef SOLVE3D
      allocate ( DIAGS(ng) % DiaU2int(LBi:UBi,LBj:UBj,NDM2d) )
      allocate ( DIAGS(ng) % DiaV2int(LBi:UBi,LBj:UBj,NDM2d) )
      allocate ( DIAGS(ng) % DiaRUfrc(LBi:UBi,LBj:UBj,3,NDM2d-1) )
      allocate ( DIAGS(ng) % DiaRVfrc(LBi:UBi,LBj:UBj,3,NDM2d-1) )
      allocate ( DIAGS(ng) % DiaU3d(LBi:UBi,LBj:UBj,N(ng),NDM3d) )
      allocate ( DIAGS(ng) % DiaV3d(LBi:UBi,LBj:UBj,N(ng),NDM3d) )
      allocate ( DIAGS(ng) % DiaU3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d) )
      allocate ( DIAGS(ng) % DiaV3wrk(LBi:UBi,LBj:UBj,N(ng),NDM3d) )
      allocate ( DIAGS(ng) % DiaRU(LBi:UBi,LBj:UBj,N(ng),2,NDrhs) )
      allocate ( DIAGS(ng) % DiaRV(LBi:UBi,LBj:UBj,N(ng),2,NDrhs) )
#  endif

# endif

      RETURN
      END SUBROUTINE allocate_diags

      SUBROUTINE initialize_diags (ng, tile)
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
# ifdef BIO_FASHAM
      DIAGS(ng) % DiaBio2d(Imin:Imax,Jmin:Jmax,1:NDbio2d) = IniVal
      DIAGS(ng) % DiaBio3d(Imin:Imax,Jmin:Jmax,1:N(ng),1:NDbio3d) =     &
     &                    IniVal
# endif

# ifdef DIAGNOSTICS_TS
      DIAGS(ng) % DiaTrc(Imin:Imax,Jmin:Jmax,1:N(ng),1:NT(ng),1:NDT) =  &
     &                   IniVal
      DIAGS(ng) % DiaTwrk(Imin:Imax,Jmin:Jmax,1:N(ng),1:NT(ng),1:NDT) = &
     &                    IniVal
# endif

# ifdef DIAGNOSTICS_UV
      DIAGS(ng) % DiaU2d(Imin:Imax,Jmin:Jmax,1:NDM2d) = IniVal
      DIAGS(ng) % DiaV2d(Imin:Imax,Jmin:Jmax,1:NDM2d) = IniVal
      DIAGS(ng) % DiaU2wrk(Imin:Imax,Jmin:Jmax,1:NDM2d) = IniVal
      DIAGS(ng) % DiaV2wrk(Imin:Imax,Jmin:Jmax,1:NDM2d) = IniVal
      DIAGS(ng) % DiaRUbar(Imin:Imax,Jmin:Jmax,2,1:NDM2d-1) = IniVal
      DIAGS(ng) % DiaRVbar(Imin:Imax,Jmin:Jmax,2,1:NDM2d-1) = IniVal

#  ifdef SOLVE3D
      DIAGS(ng) % DiaU2int(Imin:Imax,Jmin:Jmax,1:NDM2d) = IniVal
      DIAGS(ng) % DiaV2int(Imin:Imax,Jmin:Jmax,1:NDM2d) = IniVal
      DIAGS(ng) % DiaRUfrc(Imin:Imax,Jmin:Jmax,1:3,1:NDM2d-1) = IniVal
      DIAGS(ng) % DiaRVfrc(Imin:Imax,Jmin:Jmax,1:3,1:NDM2d-1) = IniVal
      DIAGS(ng) % DiaU3d(Imin:Imax,Jmin:Jmax,1:N(ng),1:NDM3d) = IniVal
      DIAGS(ng) % DiaV3d(Imin:Imax,Jmin:Jmax,1:N(ng),1:NDM3d) = IniVal
      DIAGS(ng) % DiaU3wrk(Imin:Imax,Jmin:Jmax,1:N(ng),1:NDM3d) = IniVal
      DIAGS(ng) % DiaV3wrk(Imin:Imax,Jmin:Jmax,1:N(ng),1:NDM3d) = IniVal
      DIAGS(ng) % DiaRU(Imin:Imax,Jmin:Jmax,N(ng),2,1:NDrhs) = IniVal
      DIAGS(ng) % DiaRV(Imin:Imax,Jmin:Jmax,N(ng),2,1:NDrhs) = IniVal
#  endif

# endif

      RETURN
      END SUBROUTINE initialize_diags
#endif
      END MODULE mod_diags
