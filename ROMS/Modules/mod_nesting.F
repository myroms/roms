      MODULE mod_nesting
#ifdef NESTING
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module defines generic structures used for nesting, composed,  !
!  and mosaic grids.                                                   !
!                                                                      !
!=======================================================================
!                                                 
        USE mod_kinds

        implicit none
!
! Generic structure for 2D state fields.
!
        TYPE T_NEST2D
          integer, pointer :: id             ! state field id
          integer, pointer :: gtype          ! staggered grid type
          integer, pointer :: tile(:)        ! tile partition

          integer, pointer :: LBi(:)         ! lower bound I-dimension
          integer, pointer :: UBi(:)         ! upper bound I-dimension
          integer, pointer :: LBj(:)         ! lower bound J-dimension
          integer, pointer :: UBj(:)         ! upper bound J-dimension

          integer, pointer :: Istr(:)        ! starting tile I-direction
          integer, pointer :: Iend(:)        ! ending   tile I-direction
          integer, pointer :: Jstr(:)        ! starting tile J-direction
          integer, pointer :: Jend(:)        ! ending   tile J-direction

          integer, pointer :: IstrR(:)       ! starting RHO  I-direction
          integer, pointer :: IstrU(:)       ! starting U    I-direction
          integer, pointer :: IendR(:)       ! ending   RHO  I-direction

          integer, pointer :: JstrR(:)       ! starting RHO  J-direction
          integer, pointer :: JstrV(:)       ! starting V    J-direction
          integer, pointer :: JendR(:)       ! ending   RHO  J-direction

          real (r8), pointer :: x(:,:)       ! X-positions
          real (r8), pointer :: y(:,:)       ! Y-positions
# ifdef MASKING
          real (r8), pointer :: mask(:,:)    ! land-sea masking
# endif
          real (r8), pointer :: s(:,:)       ! state array(i,j)
        END TYPE T_NEST2D

# ifdef SOLVE3D
!
! Generic structure for 3D state fields.
!
        TYPE T_NEST3D
          integer, pointer :: id             ! state field id
          integer, pointer :: gtype          ! staggered grid type
          integer, pointer :: tile(:)        ! tile partition

          integer, pointer :: LBi(:)         ! lower bound I-dimension
          integer, pointer :: UBi(:)         ! upper bound I-dimension
          integer, pointer :: LBj(:)         ! lower bound J-dimension
          integer, pointer :: UBj(:)         ! upper bound J-dimension
          integer, pointer :: LBk            ! lower bound K-dimension
          integer, pointer :: UBk            ! upper bound K-dimension

          integer, pointer :: Istr(:)        ! starting tile I-direction
          integer, pointer :: Iend(:)        ! ending   tile I-direction
          integer, pointer :: Jstr(:)        ! starting tile J-direction
          integer, pointer :: Jend(:)        ! ending   tile J-direction

          integer, pointer :: IstrR(:)       ! starting RHO  I-direction
          integer, pointer :: IstrU(:)       ! starting U    I-direction
          integer, pointer :: IendV(:)       ! ending   RHO  I-direction

          integer, pointer :: JstrR(:)       ! starting RHO  J-direction
          integer, pointer :: JstrV(:)       ! starting V    J-direction
          integer, pointer :: JendR(:)       ! ending   RHO  J-direction

          integer, pointer :: Kstr           ! starting K-index
          integer, pointer :: Kend           ! ending   K-index

          real (r8), pointer :: x(:,:)       ! X-positions
          real (r8), pointer :: y(:,:)       ! Y-positions
#  ifdef MASKING
          real (r8), pointer :: mask(:,:)    ! land-sea masking
#  endif
          real (r8), pointer :: s(:,:,:)     ! state array(i,j,k)
        END TYPE T_NEST3D
!
! Generic structure for 4D state fields.
!
        TYPE T_NEST4D
          integer, pointer :: id             ! state field id
          integer, pointer :: gtype          ! staggered grid type
          integer, pointer :: tile(:)        ! tile partition

          integer, pointer :: LBi(:)         ! lower bound I-dimension
          integer, pointer :: UBi(:)         ! upper bound I-dimension
          integer, pointer :: LBj(:)         ! lower bound J-dimension
          integer, pointer :: UBj(:)         ! upper bound J-dimension
          integer, pointer :: LBk            ! lower bound K-dimension
          integer, pointer :: UBk            ! upper bound K-dimension
          integer, pointer :: LBl            ! lower bound L-dimension
          integer, pointer :: UBl            ! upper bound L-dimension

          integer, pointer :: Istr(:)        ! starting tile I-direction
          integer, pointer :: Iend(:)        ! ending   tile I-direction
          integer, pointer :: Jstr(:)        ! starting tile J-direction
          integer, pointer :: Jend(:)        ! ending   tile J-direction

          integer, pointer :: IstrR(:)       ! starting RHO  I-direction
          integer, pointer :: IstrU(:)       ! starting U    I-direction
          integer, pointer :: IendV(:)       ! ending   RHO  I-direction

          integer, pointer :: JstrR(:)       ! starting RHO  J-direction
          integer, pointer :: JstrV(:)       ! starting V    J-direction
          integer, pointer :: JendR(:)       ! ending   RHO  J-direction

          integer, pointer :: Kstr           ! starting K-index
          integer, pointer :: Kend           ! ending   K-index
          integer, pointer :: Lstr           ! starting L-index
          integer, pointer :: Lend           ! ending   L-index

          real (r8), pointer :: x(:,:)       ! X-positions
          real (r8), pointer :: y(:,:)       ! Y-positions
#  ifdef MASKING
          real (r8), pointer :: mask(:,:)    ! land-sea masking
#  endif
          real (r8), pointer :: s(:,:,:,:)   ! state array(i,j,k,l)
        END TYPE T_NEST4D
# endif

      CONTAINS

      SUBROUTINE allocate_nesting2d (field2d, id, gtype, Tindex)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 2D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
# ifdef SOLVE3D
      USE mod_coupling
# endif
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, target :: id, gtype

      integer, intent(in), optional :: Tindex

      TYPE (T_NEST2D), pointer :: field2d(:)
!
!  Local variable declarations.
!
      integer :: ng
!
!-----------------------------------------------------------------------
!  Generic 2D structure.
!-----------------------------------------------------------------------
!
!  Allocate.
!
      IF (ASSOCIATED(field2d)) THEN
        deallocate ( field2d )
      END IF
      allocate ( field2d(Ngrids) )
!
!  Load field descriptors.
!
      DO ng=1,Ngrids
!
!  Load id, grid type, tile, array bounds, and starting/ending
!  computational indices.
!
        field2d(ng) % id    => id
        field2d(ng) % gtype => gtype
        field2d(ng) % tile  => BOUNDS(ng) % tile
        field2d(ng) % LBi   => BOUNDS(ng) % LBi
        field2d(ng) % UBi   => BOUNDS(ng) % UBi
        field2d(ng) % LBj   => BOUNDS(ng) % LBj
        field2d(ng) % UBj   => BOUNDS(ng) % UBj
        field2d(ng) % Istr  => BOUNDS(ng) % Istr
        field2d(ng) % Iend  => BOUNDS(ng) % Iend
        field2d(ng) % Jstr  => BOUNDS(ng) % Jstr
        field2d(ng) % Jend  => BOUNDS(ng) % Jend
        field2d(ng) % IstrR => BOUNDS(ng) % IstrR
        field2d(ng) % IstrU => BOUNDS(ng) % IstrU
        field2d(ng) % IendR => BOUNDS(ng) % IendR
        field2d(ng) % JstrR => BOUNDS(ng) % JstrR
        field2d(ng) % JstrV => BOUNDS(ng) % JstrV
        field2d(ng) % JendR => BOUNDS(ng) % JendR
!
!  Associate the appropriate grid arrays.
!
        IF (spherical) THEN
          IF (gtype.eq.u2dvar) THEN
            field2d(ng) % x => GRID(ng) % lonu
            field2d(ng) % y => GRID(ng) % latu
          ELSE IF (gtype.eq.v2dvar) THEN
            field2d(ng) % x => GRID(ng) % lonv
            field2d(ng) % y => GRID(ng) % latv
          ELSE
            field2d(ng) % x => GRID(ng) % lonr
            field2d(ng) % y => GRID(ng) % latr
          END IF
        ELSE
          IF (gtype.eq.u2dvar) THEN
            field2d(ng) % x => GRID(ng) % xu
            field2d(ng) % y => GRID(ng) % yu
          ELSE IF (gtype.eq.v2dvar) THEN
            field2d(ng) % x => GRID(ng) % xv
            field2d(ng) % y => GRID(ng) % yv
          ELSE
            field2d(ng) % x => GRID(ng) % xr
            field2d(ng) % y => GRID(ng) % yr
          END IF
        END IF

# ifdef MASKING
!
!  Associate the appropriate Land/Sea mask.
!
        IF (gtype.eq.u2dvar) THEN
          field2d(ng) % mask  => GRID(ng) % umask
        ELSE IF (gtype.eq.v2dvar) THEN
          field2d(ng) % mask  => GRID(ng) % vmask
        ELSE
          field2d(ng) % mask  => GRID(ng) % rmask
        END IF
# endif
!
!  Associate the appropriate state 2D array.
!
        IF (id.eq.idFsur) THEN
          field2d(ng) % s => OCEAN(ng) % zeta(:,:,Tindex)
        ELSE IF (id.eq.idUbar) THEN
          field2d(ng) % s => OCEAN(ng) % ubar(:,:,Tindex)
        ELSE IF (id.eq.idVbar) THEN
          field2d(ng) % s => OCEAN(ng) % vbar(:,:,Tindex)
        ELSE IF (id.eq.idRzet) THEN
          field2d(ng) % s => OCEAN(ng) % rzeta(:,:,Tindex)
# ifdef SOLVE3D
        ELSE IF (id.eq.idZavg) THEN
          field2d(ng) % s => COUPLING(ng) % Zt_avg1
        ELSE IF (id.eq.idUfx1) THEN
          field2d(ng) % s => COUPLING(ng) % DU_avg1
        ELSE IF (id.eq.idVfx1) THEN
          field2d(ng) % s => COUPLING(ng) % DV_avg1
# endif
        END IF

      END DO

      RETURN
      END SUBROUTINE allocate_nesting2d

# ifdef SOLVE3D

      SUBROUTINE allocate_nesting3d (field3d, id, gtype, Tindex)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 3D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, target :: id, gtype

      integer, intent(in), optional :: Tindex

      TYPE (T_NEST3D), pointer :: field3d(:)
!
!  Local variable declarations.
!
      integer :: ng
      integer, target :: LBk, UBk, Kstr, Kend
!
!-----------------------------------------------------------------------
!  Generic 3D structure.
!-----------------------------------------------------------------------
!
!  Allocate.
!
      IF (ASSOCIATED(field3d)) THEN
        deallocate ( field3d )
      END IF
      allocate ( field3d(Ngrids) )
!
!  Load field descriptors.
!
      DO ng=1,Ngrids
!
!  Load id, grid type, tile, array bounds, and starting/ending
!  computational indices.
!
        field3d(ng) % id    => id
        field3d(ng) % gtype => gtype
        field3d(ng) % tile  => BOUNDS(ng) % tile
        field3d(ng) % LBi   => BOUNDS(ng) % LBi
        field3d(ng) % UBi   => BOUNDS(ng) % UBi
        field3d(ng) % LBj   => BOUNDS(ng) % LBj
        field3d(ng) % UBj   => BOUNDS(ng) % UBj
        field3d(ng) % Istr  => BOUNDS(ng) % Istr
        field3d(ng) % Iend  => BOUNDS(ng) % Iend
        field3d(ng) % Jstr  => BOUNDS(ng) % Jstr
        field3d(ng) % Jend  => BOUNDS(ng) % Jend
        field3d(ng) % IstrR => BOUNDS(ng) % IstrR
        field3d(ng) % IstrU => BOUNDS(ng) % IstrU
        field3d(ng) % IendR => BOUNDS(ng) % IendR
        field3d(ng) % JstrR => BOUNDS(ng) % JstrR
        field3d(ng) % JstrV => BOUNDS(ng) % JstrV
        field3d(ng) % JendR => BOUNDS(ng) % JendR
!
!  Load third-dimension bounds.
!
        IF (gtype.eq.w3dvar) THEN
          LBk=0
          UBk=N(ng)
          Kstr=0
          Kend=N(ng)
          field3d(ng) % LBk   => LBk
          field3d(ng) % UBk   => UBk
          field3d(ng) % Kstr  => Kstr
          field3d(ng) % Kend  => Kend
        ELSE
          LBk=1
          UBk=N(ng)
          Kstr=1
          Kend=N(ng)
          field3d(ng) % LBk   => LBk
          field3d(ng) % UBk   => UBk
          field3d(ng) % Kstr  => Kstr
          field3d(ng) % Kend  => Kend
        END IF
!
!  Associate the appropriate grid arrays.
!
        IF (spherical) THEN
          IF (gtype.eq.u3dvar) THEN
            field3d(ng) % x => GRID(ng) % lonu
            field3d(ng) % y => GRID(ng) % latu
          ELSE IF (gtype.eq.v3dvar) THEN
            field3d(ng) % x => GRID(ng) % lonv
            field3d(ng) % y => GRID(ng) % latv
          ELSE
            field3d(ng) % x => GRID(ng) % lonr
            field3d(ng) % y => GRID(ng) % latr
          END IF
        ELSE
          IF (gtype.eq.u3dvar) THEN
            field3d(ng) % x => GRID(ng) % xu
            field3d(ng) % y => GRID(ng) % yu
          ELSE IF (gtype.eq.v3dvar) THEN
            field3d(ng) % x => GRID(ng) % xv
            field3d(ng) % y => GRID(ng) % yv
          ELSE
            field3d(ng) % x => GRID(ng) % xr
            field3d(ng) % y => GRID(ng) % yr
          END IF
        END IF

#  ifdef MASKING
!
!  Associate the appropriate Land/Sea mask.
!
        IF (gtype.eq.u3dvar) THEN
          field3d(ng) % mask  => GRID(ng) % umask
        ELSE IF (gtype.eq.v3dvar) THEN
          field3d(ng) % mask  => GRID(ng) % vmask
        ELSE
          field3d(ng) % mask  => GRID(ng) % rmask
        END IF
#  endif
!
!  Associate the appropriate state 3D array.
!
        IF (id.eq.idUvel) THEN
          field3d(ng) % s => OCEAN(ng) % u(:,:,:,Tindex)
        ELSE IF (id.eq.idVvel) THEN
          field3d(ng) % s => OCEAN(ng) % v(:,:,:,Tindex)
        ELSE IF (id.eq.idRu3d) THEN
          field3d(ng) % s => OCEAN(ng) % ru(:,:,:,Tindex)
        ELSE IF (id.eq.idRv3d) THEN
          field3d(ng) % s => OCEAN(ng) % rv(:,:,:,Tindex)
        ELSE IF (id.eq.idOvel) THEN
          field3d(ng) % s => OCEAN(ng) % W
        ELSE IF (id.eq.idWvel) THEN
          field3d(ng) % s => OCEAN(ng) % wvel
        ELSE IF (id.eq.idDano) THEN
          field3d(ng) % s => OCEAN(ng) % rho
        ELSE IF (id.eq.idVvis) THEN
          field3d(ng) % s => MIXING(ng) % Akv
#  if defined GLS_MIXING || defined MY25_MIXING
        ELSE IF (id.eq.idMtke) THEN
          field3d(ng) % s => MIXING(ng) % tke(:,:,:,Tindex)
        ELSE IF (id.eq.idMtls) THEN
          field3d(ng) % s => MIXING(ng) % gls(:,:,:,Tindex)
        ELSE IF (id.eq.idVmKK) THEN
          field3d(ng) % s => MIXING(ng) % Akk
#   ifdef GLS_MIXING
        ELSE IF (id.eq.idVvel) THEN
          field3d(ng) % s => MIXING(ng) % Akp
#   endif
#  endif
        END IF

      END DO

      RETURN
      END SUBROUTINE allocate_nesting3d

      SUBROUTINE allocate_nesting4d (field4d, id, gtype, Tindex)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes nesting structure for 4D     !
!  state variables.                                                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_mixing
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, target :: id, gtype

      integer, intent(in), optional :: Tindex

      TYPE (T_NEST4D), pointer :: field4d(:)
!
!  Local variable declarations.
!
      integer :: ng
      integer, target :: LBk, UBk, LBl, UBl
      integer, target :: Kstr, Kend, Lstr, Lend
!
!-----------------------------------------------------------------------
!  Generic 4D structure.
!-----------------------------------------------------------------------
!
!  Allocate.
!
      IF (associated(field4d)) THEN
        deallocate ( field4d )
      END IF
      allocate ( field4d(Ngrids) )
!
!  Load field descriptors.
!
      DO ng=1,Ngrids
!
!  Load id, grid type, tile, array bounds, and starting/ending
!  computational indices.
!
        field4d(ng) % id    => id
        field4d(ng) % gtype => gtype
        field4d(ng) % tile  => BOUNDS(ng) % tile
        field4d(ng) % LBi   => BOUNDS(ng) % LBi
        field4d(ng) % UBi   => BOUNDS(ng) % UBi
        field4d(ng) % LBj   => BOUNDS(ng) % LBj
        field4d(ng) % UBj   => BOUNDS(ng) % UBj
        field4d(ng) % Istr  => BOUNDS(ng) % Istr
        field4d(ng) % Iend  => BOUNDS(ng) % Iend
        field4d(ng) % Jstr  => BOUNDS(ng) % Jstr
        field4d(ng) % Jend  => BOUNDS(ng) % Jend
        field4d(ng) % IstrR => BOUNDS(ng) % IstrR
        field4d(ng) % IstrU => BOUNDS(ng) % IstrU
        field4d(ng) % IendR => BOUNDS(ng) % IendR
        field4d(ng) % JstrR => BOUNDS(ng) % JstrR
        field4d(ng) % JstrV => BOUNDS(ng) % JstrV
        field4d(ng) % JendR => BOUNDS(ng) % JendR
!
!  Load third-dimension bounds.
!
        IF (gtype.eq.r3dvar) THEN
          LBk=1
          UBk=N(ng)
          Kstr=1
          Kend=N(ng)
          field4d(ng) % LBk   => LBk
          field4d(ng) % UBk   => UBk
          field4d(ng) % Kstr  => Kstr
          field4d(ng) % Kend  => Kend
        ELSE IF
        END IF
!
!  Associate the appropriate grid arrays.
!
        IF (spherical) THEN
          IF (gtype.eq.u3dvar) THEN
            field4d(ng) % x => GRID(ng) % lonu
            field4d(ng) % y => GRID(ng) % latu
          ELSE IF (gtype.eq.v3dvar) THEN
            field4d(ng) % x => GRID(ng) % lonv
            field4d(ng) % y => GRID(ng) % latv
          ELSE
            field4d(ng) % x => GRID(ng) % lonr
            field4d(ng) % y => GRID(ng) % latr
          END IF
        ELSE
          IF (gtype.eq.u3dvar) THEN
            field4d(ng) % x => GRID(ng) % xu
            field4d(ng) % y => GRID(ng) % yu
          ELSE IF (gtype.eq.v3dvar) THEN
            field4d(ng) % x => GRID(ng) % xv
            field4d(ng) % y => GRID(ng) % yv
          ELSE
            field4d(ng) % x => GRID(ng) % xr
            field4d(ng) % y => GRID(ng) % yr
          END IF
        END IF

#  ifdef MASKING
!
!  Associate the appropriate Land/Sea mask.
!
        IF (gtype.eq.u3dvar) THEN
          field4d(ng) % mask  => GRID(ng) % umask
        ELSE IF (gtype.eq.v3dvar) THEN
          field4d(ng) % mask  => GRID(ng) % vmask
        ELSE
          field4d(ng) % mask  => GRID(ng) % rmask
        END IF
#  endif
!
!  Associate the appropriate state 4D array.  Notice that temperature
!  triggers associating all tracers.
!
        IF (id.eq.idTvar(itemp)) THEN
          LBl=1
          UBl=NT(ng)
          Lstr=1
          Lend=NT(ng)
          field4d(ng) % LBl   => LBl
          field4d(ng) % UBl   => UBl
          field4d(ng) % Lstr  => Lstr
          field4d(ng) % Lend  => Lend
          field4d(ng) % s     => OCEAN(ng) % t(:,:,:,Tindex,:)
        END IF

      END DO

      RETURN
      END SUBROUTINE allocate_nesting4d
# endif

#endif
      END MODULE mod_nesting
