#include "cppdefs.h"
      MODULE roms_vchange_mod

#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  It has several routines to process nonlinear and linear variable    !
!  changes needed for the H(x) operator, which computes the model      !
!  state variables at the observation locations.                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE roms_obs_mod
!
# ifdef DISTRIBUTE
#  ifdef ADJOINT
      USE mp_exchange_mod,  ONLY : ad_mp_exchange3d
#  endif
      USE mp_exchange_mod,  ONLY : mp_exchange3d
# endif
      USE tracer_index_mod, ONLY : tracer_index
!
!-----------------------------------------------------------------------
!  Observation Variable Change Operator: CLASS(obs_VarCha)
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: obs_VarCha
!
!  Geometry pointers.
!
        real (r8), pointer :: mask(:,:)      => NULL() ! land-sea mask
# ifdef SOLVE3D
        real (r8), pointer :: depth(:,:,:)   => NULL() ! depths (m)
# endif
!
!  Generic field arrays: Allocatable arrays are safer and preferred
!  because pointer association to array sections has a unity lower
!  bound for each dimension. This is because the right-hand side is
!  treated as a new array expression that destroys the original
!  parallel decomposition
!
        real (r8), pointer :: field2d(:,:)   => NULL()
# ifdef SOLVE3D
        real (r8), pointer :: field3d(:,:,:) => NULL()
# endif
# ifdef SOLVE3D
!
!  Generic A-grid and C-grid vector components.
!
        real (r8), pointer :: Ua(:,:,:)      => NULL()
        real (r8), pointer :: Va(:,:,:)      => NULL()
!
        real (r8), pointer :: Uc(:,:,:)      => NULL()
        real (r8), pointer :: Vc(:,:,:)      => NULL()
# endif
!
        CONTAINS
!
!  Execution.
!
        PROCEDURE :: apply2d        => var_change_apply2d
# ifdef SOLVE3D
        PROCEDURE :: apply3d        => var_change_apply3d
# endif
# ifdef ADJOINT
        PROCEDURE :: apply2d_ad     => var_change_apply2d_ad
#  ifdef SOLVE3D
        PROCEDURE :: apply3d_ad     => var_change_apply3d_ad
#  endif
# endif
# ifdef SOLVE3D
!
!  Vertical interpolation: k-level to constant depth.
!
        PROCEDURE :: interp_k2z     => var_change_interp_k2z
#  ifdef ADJOINT
        PROCEDURE :: interp_k2z_ad  => var_change_interp_k2z_ad
#  endif
!
!  Horizontal velocity component from C-grid to A-grid and vice versa.
!
        PROCEDURE :: uv_a2c         => var_change_uv_a2c
        PROCEDURE :: uv_c2a         => var_change_uv_c2a
#  ifdef ADJOINT
        PROCEDURE :: uv_a2c_ad      => var_change_uv_a2c
        PROCEDURE :: uv_c2a_ad      => var_change_uv_c2a
#  endif
# endif

      END TYPE obs_VarCha
!
!-----------------------------------------------------------------------
!  Module variables.
!-----------------------------------------------------------------------
!
!  Nonlinear and linear variable changes object.
!
      TYPE (obs_VarCha), allocatable :: VarCha(:)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      PUBLIC :: vchange_initialize
      PUBLIC :: vchange_finalize
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      CONTAINS
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!  It allocates the variable change object.
!
      SUBROUTINE vchange_initialize
!
!  Allocate variable changes object.
!
      IF (.not.allocated(VarCha)) THEN
        allocate ( VarCha(Ngrids) )
      END IF
!
      RETURN
      END SUBROUTINE vchange_initialize
!
!-----------------------------------------------------------------------
!  It destoys variable change object.
!
      SUBROUTINE vchange_finalize (ng, model)
!
      integer,          intent(in   ) :: ng         ! nested grid
      integer,          intent(in   ) :: model      ! kernel identifier
!
!  Deallocate generic arrays.
!
      IF (associated(VarCha(ng)%field2d)) nullify (VarCha(ng)%field2d)
# ifdef SOLVE3D
      IF (associated(VarCha(ng)%field3d)) nullify (VarCha(ng)%field3d)
      IF (associated(VarCha(ng)%Ua))      nullify (VarCha(ng)%Ua)
      IF (associated(VarCha(ng)%Va))      nullify (VarCha(ng)%Va)
      IF (associated(VarCha(ng)%Uc))      nullify (VarCha(ng)%Uc)
      IF (associated(VarCha(ng)%Vc))      nullify (VarCha(ng)%Vc)
# endif
!
!  Nullify pointers.
!
      IF (associated(VarCha(ng)%mask))   nullify (VarCha(ng)%mask)
# ifdef SOLVE3D
      IF (associated(VarCha(ng)%depth))  nullify (VarCha(ng)%depth)
# endif
!
      RETURN
      END SUBROUTINE vchange_finalize
!
! <><><><><><><><><><><><><><><><><><><><><><><> CLASS OBS_VARCHA <><><>
!
!  Given an observation type, it retuns a pointer for the 2D state
!  variable associated with the observation, which is used to compute
!  the model at the observtion locations, H(x). If appropriate, it
!  applies the necessary nonlinear or linear variable changes. If no
!  variable changes, the transformation is an identity operator.
!
      SUBROUTINE var_change_apply2d (self, ng, tile, model, varname)
!
      CLASS (obs_VarCha), intent(inout) :: self     ! VarChan object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: tile     ! domain partition
      integer,            intent(in   ) :: model    ! kernel identifier
      character (len=*),  intent(in   ) :: varname  ! state variable
!
      integer                           :: is, ie, js, je, ks, ke
      integer                           :: itrc, nz
!
      real (r8),    target, allocatable :: field(:,:,:)
      real (r8),            allocatable :: Zeven(:)
!
!  Initialize.
!
      is=BOUNDS(ng)%LBi(tile)
      ie=BOUNDS(ng)%UBi(tile)
      js=BOUNDS(ng)%LBj(tile)
      je=BOUNDS(ng)%UBj(tile)
      ks=1
      ke=N(ng)
!
!  Set generic variable change array, VarCha%field2d, for the state
!  variable associated with the observations. Then, compute the 
!  appropriate variable change(s), if any. Otherwise, the default is
!  the identity operator.
!
# ifdef SOLVE3D
      IF (TRIM(varname).eq.'Ucodar') THEN
        self%Ua => NULL()
        self%Va => NULL()
        self%Uc => NULL()
        self%Vc => NULL()
      END IF
# endif
!
      SELECT CASE (TRIM(varname))
!
        CASE ('adt', 'ssh', 'ADT', 'SSH')
!
          self%field2d => NULL()
!
          IF (model.eq.iNLM) THEN
            self%field2d(is:,js:) => OCEAN(ng)%zeta(is:ie,js:je,        &
     &                                              KOUT)
# ifdef ADJOINT
          ELSE IF (model.eq.iADM) THEN
            self%field2d => OCEAN(ng)%f_zeta
# endif
# ifdef TLM_OBS
          ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
            self%field2d(is:,js:) => OCEAN(ng)%tl_zeta(is:ie,js:je,     &
     &                                                 KOUT)
# endif
# ifndef VERIFICATION
          ELSE IF (model.eq.iSTD) THEN
            self%field2d(is:,js:) => OCEAN(ng)%e_zeta(is:ie,js:je,1)
# endif
          END IF

# ifdef SOLVE3D
!
        CASE ('sss', 'sst', 'SSS', 'SST')
!
          self%field2d => NULL()
!
          itrc=tracer_index(varname)
          IF (model.eq.iNLM) THEN
            self%field2d(is:,js:) => OCEAN(ng)%t(is:ie,js:je,           &
     &                                           N(ng),NOUT,itrc)
#  ifdef ADJOINT
          ELSE IF (model.eq.iADM) THEN
            self%field2d(is:,js:) => OCEAN(ng)%f_t(is:ie,js:je,         &
     &                                             N(ng),itrc)
#  endif
#  ifdef TLM_OBS
          ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
            self%field2d(is:,js:) => OCEAN(ng)%tl_t(is:ie,js:je,        &
     &                                              N(ng),NOUT,itrc)
#  endif
#  ifndef VERIFICATION
          ELSE IF (model.eq.iSTD) THEN
            self%field2d(is:,js:) => OCEAN(ng)%e_t(is:ie,js:je,         &
     &                                             N(ng),1,itrc)
#  endif
          END IF
!
        CASE ('Ucodar', 'Vcodar')
!
          self%field2d => NULL()
!
!  Transform u- and v-momentum component from C-grid to A-grid. The
!  A-grid velocities are computed once during the "Ucodar" case.
!
          IF (TRIM(varname).eq.'Ucodar') THEN
!
            self%Ua => NULL()
            self%Va => NULL()
            allocate ( self%Ua(is:ie, js:je, ks:ke) )
            allocate ( self%Va(is:ie, js:je, ks:ke) )
!
            self%Uc => NULL()
            self%Vc => NULL()
!
            IF (model.eq.iNLM) THEN
              self%Uc(is:ie,js:je,ks:ke)=OCEAN(ng)%u(is:ie,js:je,       &
     &                                               ks:ke,NOUT)
              self%Vc(is:ie,js:je,ks:ke)=OCEAN(ng)%v(is:ie,js:je,       &
     &                                               ks:ke,NOUT)
#  ifdef TLM_OBS
            ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
              self%Uc(is:ie,js:je,ks:ke)=OCEAN(ng)%tl_u(is:ie,js:je,    &
     &                                                  ks:ke,NOUT)
              self%Vc(is:ie,js:je,ks:ke)=OCEAN(ng)%tl_v(is:ie,js:je,    &
     &                                                  ks:ke,NOUT)
#  endif
            END IF
            CALL self%uv_c2a (ng, tile, model)
          END IF
!
!  Interpolate A-grid u-component to constant depth, Zeven = -2m.
!
          nz=1
!
          IF (.not.allocated(field)) THEN
            allocate ( field(is:ie, js:je, nz) )
            field=0.0_r8
          END IF
          IF (.not.allocated(Zeven)) THEN
            allocate ( Zeven(nz) )
            Zeven=-2.0_r8
          END IF
!
          IF (TRIM(varname).eq.'Ucodar') THEN
            self%mask  => GRID(ng)%rmask
            self%depth => GRID(ng)%z_r
            CALL self%interp_k2z (ng, tile, model,                      &
     &                            is, ie, js, je,                       &
     &                            varname, Zeven,                       &
     &                            self%Ua, field)
          ELSE IF (TRIM(varname).eq.'Vcodar') THEN
            CALL self%interp_k2z (ng, tile, model,                      &
     &                            is, ie, js, je,                       &
     &                            varname, Zeven,                       &
     &                            self%Va, field)
            nullify (self%mask, self%depth)
          END IF
!
          self%field2d(is:ie,js:je)=field(is:ie,js:je,1)
          deallocate (field)
# endif
      END SELECT
!
      RETURN
      END SUBROUTINE var_change_apply2d

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Given an observation type, it retuns a pointer for the 3D state
!  variable associated with the observation, which is used to compute
!  the model at the observtion locations, H(x). If appropriate, it
!  applies the necessary nonlinear or linear variable changes. If no
!  variable changes, the transformation is an identity operator.
!
      SUBROUTINE var_change_apply3d (self, ng, tile, model, varname)
!
      CLASS (obs_VarCha), intent(inout) :: self     ! VarChan object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: tile     ! domain partition
      integer,            intent(in   ) :: model    ! kernel identifier
      character (len=*),  intent(in   ) :: varname
!
      integer                           :: is, ie, js, je, ks, ke
      integer                           :: itrc
!
!  Initialize.
!
      is=BOUNDS(ng)%LBi(tile)
      ie=BOUNDS(ng)%UBi(tile)
      js=BOUNDS(ng)%LBj(tile)
      je=BOUNDS(ng)%UBj(tile)
      ks=1
      ke=N(ng)
!
!  Set pointer (VarCha%field3d) to state variable associated with the
!  observations. Compute the appropriate variable change(s), if any.
!  Otherwise, the default is the identity operator.
!
      SELECT CASE (TRIM(varname))
!
        CASE ('salt', 'temp', 'tracer')
!
          self%field3d => NULL()
!
          itrc=tracer_index(varname)
          IF (model.eq.iNLM) THEN
            self%field3d(is:,js:,ks:) => OCEAN(ng)%t(is:ie,js:je,       &
     &                                               ks:ke,NOUT,itrc)
#  ifdef ADJOINT
          ELSE IF (model.eq.iADM) THEN
            self%field3d(is:,js:,ks:) => OCEAN(ng)%f_t(is:ie,js:je,     &
     &                                                 ks:ke,itrc)
#  endif
#  ifdef TLM_OBS
          ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
            self%field3d(is:,js:,ks:) => OCEAN(ng)%tl_t(is:ie,js:je,    &
     &                                                  ks:ke,NOUT,itrc)
#  endif
#  ifndef VERIFICATION
          ELSE IF (model.eq.iSTD) THEN
            self%field3d(is:,js:,ks:) => OCEAN(ng)%e_t(is:ie,js:je,     &
     &                                                 ks:ke,1,itrc)
#  endif
          END IF
      END SELECT
!
      RETURN
      END SUBROUTINE var_change_apply3d
# endif

# ifdef ADJOINT
!
!-----------------------------------------------------------------------
!  Adjoint of applying the necessary linear 2D variable changes.
!
      SUBROUTINE var_change_apply2d_ad (ad_self, ng, tile, model,       &
     &                                  varname)
!
      CLASS (obs_VarCha), intent(inout) :: ad_self  ! VarChan object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: tile     ! domain partition
      integer,            intent(in   ) :: model    ! kernel identifier
      character (len=*),  intent(in   ) :: varname
!
      integer                           :: LBi, UBi, LBj, UBj
      integer                           :: i, j, nz
!
      real (r8),    target, allocatable :: ad_field(:,:,:)
      real (r8),            allocatable :: Zeven(:)
!
!  Initialize.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set pointer (ad_self%field2d) to state variable associated with the
!  observations. Compute the appropriate variable  change(s), if any.
!  Otherwise, the default is the identity operator.
!
      SELECT CASE (TRIM(varname))

#  ifdef SOLVE3D
        CASE ('Ucodar', 'Vcodar')
!
!  Adjoint of pointer association.
!
          nz=1
!
          IF (.not.allocated(ad_field)) THEN
            allocate ( ad_field(LBi:UBi, LBj:UBj, nz) )
            ad_field=0.0_r8
          END IF
          IF (.not.allocated(Zeven)) THEN
            allocate ( Zeven(nz) )
            Zeven=-2.0_r8
          END IF
!
          DO j=LBj,UBj
            DO i=LBi,UBi
!>            self%field2d(i,j)=field(i,j,1)
!>
              ad_field(i,j,1)=ad_field(i,j,1)+ad_self%field2d(i,j)
            END DO
          END DO
!
!  Adjoint of interpolating A-grid components to constant depth,
!  Zeven = -2m.
!
          IF (TRIM(varname).eq.'Ucodar') THEN
!>          CALL self%interp_k2z (ng, tile, model,                      &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            varname, Zeven,                       &
!>   &                            self%Ua, field)
!>
            CALL ad_self%interp_k2z_ad (ng, tile, model,                &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  varname, Zeven,                 &
     &                                  ad_self%Ua, ad_field)
          ELSE IF (TRIM(varname).eq.'Vcodar') THEN
!>          CALL self%interp_k2z (ng, tile, model,                      &
!>   &                            LBi, UBi, LBj, UBj,                   &
!>   &                            varname, Zeven,                       &
!>   &                            self%Va, field)
!>
            CALL ad_self%interp_k2z_ad (ng, tile, model,                &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  varname, Zeven,                 &
     &                                  ad_self%Va, ad_field)
          END IF
          deallocate (ad_field)
!
!  Adjoint of transforming u- and v-momentum component from C-grid to
!  A-grid.
!
!>        IF (TRIM(varname).eq.'Ucodar') THEN
!>
          IF (TRIM(varname).eq.'Vcodar') THEN
!>          CALL self%uv_c2a (ng, tile, model)
!>
            CALL ad_self%uv_c2a_ad (ng, tile, model)
!
            DO k=1,N(ng)
              DO j=LBj,UBj
                DO i=LBi,UBi
                  ad_self%Ua(i,j,k)=ad_self%Ua(i,j,k)+ad_self%Uc(i,j,k)
                  ad_self%Va(i,j,k)=ad_self%Va(i,j,k)+ad_self%Vc(i,j,k)
                END DO
              END DO
            END DO
          END IF
#  endif
      END SELECT
!
      RETURN
      END SUBROUTINE var_change_apply2d_ad

#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint of applying the necessary linear 3D variable changes.
!
      SUBROUTINE var_change_apply3d_ad (ad_self, ng, tile, model,       &
     &                                  varname)
!
      CLASS (obs_VarCha), intent(inout) :: ad_self  ! VarChan object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: tile     ! domain partition
      integer,            intent(in   ) :: model    ! kernel identifier
      character (len=*),  intent(in   ) :: varname
!
      integer                           :: itrc
!
!  Set pointer (ad_self%field2d) to state variable associated with the
!  observations. Compute the appropriate variable  change(s), if any.
!  Otherwise, the default is the identity operator.
!
      SELECT CASE (TRIM(varname))
        CASE ('salt', 'temp', 'tracer')
          itrc=tracer_index(varname)
!>
!> Currently, no variable changes are needed for 3D variables.
!>
      END SELECT
!
      RETURN
      END SUBROUTINE var_change_apply3d_ad
#  endif
# endif

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  It vertically interpolates a 3D terrain-following field to specied
!  constant depth or depths.
!
      SUBROUTINE var_change_interp_k2z (self, ng, tile, model,          &
     &                                  LBi, UBi, LBj, UBj,             &
     &                                  name, Zeven, Fk, Fz)
!
      CLASS (obs_VarCha), intent(inout) :: self     ! VarChan object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: tile     ! domain partition
      integer,            intent(in   ) :: model    ! kernel identifier
      integer,            intent(in   ) :: LBi, UBi ! array I-bounds
      integer,            intent(in   ) :: LBj, UBj ! array J-bounds
      character (len=*),  intent(in   ) :: name     ! field name
      real (r8),          intent(in   ) :: Zeven(:) ! requested depth(s)
      real (r8),          intent(in   ) :: Fk(LBi:,LBj:,:)
      real (r8),          intent(out  ) :: Fz(LBi:,LBj:,:)
!
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: IstrR, IendR, JstrR, JendR
      integer                           :: i, iz, j, k, k1, k2
      real (r8)                         :: Zbot, Ztop, dz, r1, r2
!
!  Initialize.
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      IstrR=BOUNDS(ng)%IstrR(tile)
      IendR=BOUNDS(ng)%IendR(tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JendR=BOUNDS(ng)%JendR(tile)
!
!  Interpolate 3D field to requested constant depth(s).
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
#  ifdef MASKING
          IF (self%mask(i,j).lt.0.5_r8) CYCLE
#  endif
          DO iz=1,SIZE(Zeven)
            Ztop=self%depth(i,j,N(ng))
            Zbot=self%depth(i,j,1    )
            IF (Zeven(iz).ge.Ztop) THEN        ! If shallower, assign
              k1=N(ng)                         ! to top grid cell. The
              k2=N(ng)                         ! field is located on the
              r1=1.0_r8                        ! upper cell half above
              r2=0.0_r8                        ! its middle depth.
              Fz(i,j,iz)=Fk(i,j,k1)
            ELSE IF (Zbot.ge.Zeven(iz)) THEN
              r1=0.0_r8                        ! If deeper, ignore.
              r2=0.0_r8                        ! return zero value.
              Fz(i,j,iz)=0.0_r8
            ELSE
              DO k=N(ng),2,-1                  ! Otherwise, interpolate
                Ztop=self%depth(i,j,k  )       ! to requested depth
                Zbot=self%depth(i,j,k-1)
                IF ((Ztop.gt.Zeven(iz)).and.                            &
     &              (Zeven(iz).ge.Zbot)) THEN
                  k1=k-1
                  k2=k
                END IF
              END DO
              dz=self%depth(i,j,k2)-self%depth(i,j,k1)
              r2=(Zeven(iz)-self%depth(i,j,k1))/dz
              r1=1.0_r8-r2
              IF ((r1+r2).gt.0.0_r8) THEN
                Fz(i,j,iz)=r1*Fk(i,j,k1)+r2*Fk(i,j,k2)
              END IF
            END IF
          END DO
        END DO
      END DO
!
      RETURN
      END SUBROUTINE var_change_interp_k2z
!
!-----------------------------------------------------------------------
!  It transforms vector components from A-grid (cell center) to
!  staggered Arakawa C-grid.
!
      SUBROUTINE var_change_uv_a2c (self, ng, tile, model)
!
      CLASS (obs_VarCha), intent(inout) :: self     ! VarChan object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: tile     ! domain partition
      integer,            intent(in   ) :: model    ! kernel identifier
!
      integer                           :: LBi, UBi, LBj, UBj
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: IstrR, IendR, JstrR, JendR
      integer                           :: i, j, k
!
      real (r8),            allocatable :: Urho(:,:), Vrho(:,:)
!
!  Initialize.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      IstrR=BOUNDS(ng)%IstrR(tile)
      IendR=BOUNDS(ng)%IendR(tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JendR=BOUNDS(ng)%JendR(tile)
!
#  ifdef NESTING
      IminS=Istr-4
      ImaxS=Iend+3
      JminS=Jstr-4
      JmaxS=Jend+3
#  else
      IminS=Istr-3
      ImaxS=Iend+3
      JminS=Jstr-3
      JmaxS=Jend+3
#  endif
!
      allocate ( Urho(IminS:ImaxS, JminS:JmaxS) )
      Urho=0.0_r8
      allocate ( Vrho(IminS:ImaxS, JminS:JmaxS) )
      Vrho=0.0_r8
!
!  A-grid to C-grid vector transformation.
!
      K_LOOP : DO k=1,N(ng)
!
!  Rotate from geographical Eastward and Northward directions to
!  computational (XI,ETA) directions.
!
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
            Urho(i,j)=self%Ua(i,j,k)*GRID(ng)%CosAngler(i,j)+           &
     &                self%Va(i,j,k)*GRID(ng)%SinAngler(i,j)
            Vrho(i,j)=self%Va(i,j,k)*GRID(ng)%CosAngler(i,j)-           &
     &                self%Ua(i,j,k)*GRID(ng)%SinAngler(i,j)
          END DO
        END DO
!
!  Compute staggered C-grid components.
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
            self%Uc(i,j,k)=0.5_r8*(Urho(i-1,j)+Urho(i,j))
#  ifdef MASKING
            self%Uc(i,j,k)=self%Uc(i,j,k)*GRID(ng)%umask(i,j)
#  endif
          END DO
        END DO
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
            self%Vc(i,j,k)=0.5_r8*(Vrho(i,j-1)+Vrho(i,j))
#  ifdef MASKING
            self%Vc(i,j,k)=self%Vc(i,j,k)*GRID(ng)%vmask(i,j)
#  endif
          END DO
        END DO
      END DO K_LOOP
!
!  Deallocate local arrays.
!
      IF (allocated(Urho)) deallocate (Urho)
      IF (allocated(Vrho)) deallocate (Vrho)

#  ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(g),                  &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    self%Uc, self%Vc)
#  endif
!
      RETURN
      END SUBROUTINE var_change_uv_a2c
!
!-----------------------------------------------------------------------
!  It transforms vector components from staggered Arakawa C-grid to
!  A-grid (cell center).
!
      SUBROUTINE var_change_uv_c2a (self, ng, tile, model)
!
      CLASS (obs_VarCha), intent(inout) :: self   ! VarChan object
      integer,            intent(in   ) :: ng     ! nested grid
      integer,            intent(in   ) :: tile   ! domain partition
      integer,            intent(in   ) :: model  ! kernel identifier
!
      integer                           :: LBi, UBi, LBj, UBj
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: IstrR, IendR, JstrE, JendR
      integer                           :: i, j, k
!
      real (r8)                         :: UaCosA, UaSinA, VaCosA, VaSinA
!
!  Initialize.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      IstrR=BOUNDS(ng)%IstrR(tile)
      IendR=BOUNDS(ng)%IendR(tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JendR=BOUNDS(ng)%JendR(tile)
!
!  C-grid to A-grid vector transformation.
!
      K_LOOP : DO k=1,N(ng)
!
!  Compute A-grid (cell center) vector components and apply lateral
!  boundary conditions.
!
        DO j=JstrR,JendR
          DO i=Istr,Iend
            self%Ua(i,j,k)=0.5_r8*(self%Uc(i,j,k)+self%Uc(i+1,j,k))
            IF (.not.EWperiodic(ng)) THEN
              IF (DOMAIN(ng)%Western_Edge(tile)) THEN
                self%Ua(Istr-1,j,k)=self%Ua(Istr,j,k)
              END IF
              IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
                self%Ua(Iend+1,j,k)=self%Ua(Iend,j,k)
              END IF
            END IF
          END DO
        END DO
!
        DO j=Jstr,Jend
          DO i=IstrR,IendR
            self%Va(i,j,k)=0.5_r8*(self%Vc(i,j,k)+self%Vc(i,j+1,k))
            IF (.not.NSperiodic(ng)) THEN
              IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
                self%Va(i,Jstr-1,k)=self%Va(i,Jstr,k)
              END IF
              IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
                self%Va(i,Jend+1,k)=self%Va(i,Jend,k)
              END IF
            END IF
          END DO
        END DO
!
!  Rotate from computational to gegraphical Eastward and Northward
!  directions.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            UaCosA=self%Ua(i,j,k)*GRID(ng)%CosAngler(i,j)
            UaSinA=self%Ua(i,j,k)*GRID(ng)%SinAngler(i,j)
            VaCosA=self%Va(i,j,k)*GRID(ng)%CosAngler(i,j)
            VaSinA=self%Va(i,j,k)*GRID(ng)%SinAngler(i,j)
            self%Ua(i,j,k)=UaCosA-VaSinA
            self%Va(i,j,k)=VaCosA+UaSinA
#  ifdef MASKING
            self%Ua(i,j,k)=self%Ua(i,j,k)*GRID(ng)%rmask(i,j)
            self%Va(i,j,k)=self%Va(i,j,k)*GRID(ng)%rmask(i,j)
#  endif
          END DO
        END DO
      END DO K_LOOP

#  ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
      CALL mp_exchange3d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    self%Ua, self%Va)
#  endif
!
      RETURN
      END SUBROUTINE var_change_uv_c2a

#  ifdef ADJOINT
!
!-----------------------------------------------------------------------
!  Adjoint of vertically interpolating a 3D terrain-following field to
!  specified constant depth or depths.
!
      SUBROUTINE var_change_interp_k2z_ad (ad_self, ng, tile, model,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     name, Zeven, ad_Fk, ad_Fz)
!
      CLASS (obs_VarCha), intent(inout) :: ad_self  ! VarChan object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: tile     ! domain partition
      integer,            intent(in   ) :: model    ! kernel identifier
      integer,            intent(in   ) :: LBi, UBi ! array I-bounds
      integer,            intent(in   ) :: LBj, UBj ! array J-bounds
      character (len=*),  intent(in   ) :: name     ! field name
      real (r8),          intent(in   ) :: Zeven(:) ! requested depth(s)
      real (r8),          intent(inout) :: ad_Fk(LBi:,LBj:,:)
      real (r8),          intent(inout) :: ad_Fz(LBi:,LBj:,:)
!
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: IstrR, IendR, JstrR, JendR
      integer                           :: i, iz, j, k, k1, k2
      real (r8)                         :: Zbot, Ztop, dz, r1, r2
!
!  Initialize.
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      IstrR=BOUNDS(ng)%IstrR(tile)
      IendR=BOUNDS(ng)%IendR(tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JendR=BOUNDS(ng)%JendR(tile)
!
!  Interpolate 3D field to requested constant depth(s).
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
#  ifdef MASKING
          IF (ad_self%mask(i,j).lt.0.5_r8) CYCLE
#  endif
          DO iz=1,SIZE(Zeven)
            Ztop=ad_self%depth(i,j,N(ng))
            Zbot=ad_self%depth(i,j,1    )
            IF (Zeven(iz).ge.Ztop) THEN        ! If shallower, assign
              k1=N(ng)                         ! to top grid cell. The
              k2=N(ng)                         ! field is located on the
              r1=1.0_r8                        ! upper cell half above
              r2=0.0_r8                        ! its middle depth.
!>            Fz(i,j,iz)=Fk(i,j,k1)
!>
              ad_Fk(i,j,k1)=ad_Fk(i,j,k1)+ad_Fz(i,j,iz)
              ad_Fz(i,j,iz)=0.0_r8
            ELSE IF (Zbot.ge.Zeven(iz)) THEN
              r1=0.0_r8                        ! If deeper, ignore.
              r2=0.0_r8                        ! return zero value.
!>            Fz(i,j,iz)=0.0_r8
!>
              ad_Fz(i,j,iz)=0.0_r8
            ELSE
              DO k=N(ng),2,-1                  ! Otherwise, interpolate
                Ztop=ad_self%depth(i,j,k    )  ! to requested depth
                Zbot=ad_self%depth(i1,j1,k-1)
                IF ((Ztop.gt.Zeven(iz)).and.                            &
     &              (Zeven(i).ge.Zbot)) THEN
                  k1=k-1
                  k2=k
                END IF
              END DO
              dz=ad_self%depth(i,j,k2)-ad_self%depth(i,j,k1)
              r2=(Zeven(iz)-ad_self%depth(i1,j1,k1))/dz
              r1=1.0_r8-r2
              IF ((r1+r2).gt.0.0_r8) THEN
!>              Fz(i,j,iz)=r1*Fk(i,j,k1)+r2*Fk(i,j,k2)
!>
                ad_Fk(i,j,k1)=ad_Fk(i,j,k1)+r1*ad_Fz(i,j,iz)
                ad_Fk(i,j,k2)=ad_Fk(i,j,k2)+r2*ad_Fz(i,j,iz)
                ad_Fz(i,j,iz)=0.0_r8
              END IF
            END IF
          END DO
        END DO
      END DO
!
      RETURN
      END SUBROUTINE var_change_interp_k2z_ad
!
!-----------------------------------------------------------------------
!  Adjoint of transforming vector components from A-grid (cell center)
!  to staggered Arakawa C-grid.
!
      SUBROUTINE var_change_uv_a2c_ad (ad_self, ng, tile, model)
!
      CLASS (obs_VarCha), intent(inout) :: ad_self ! VarChan object
      integer,            intent(in   ) :: ng      ! nested grid
      integer,            intent(in   ) :: tile    ! domain partition
      integer,            intent(in   ) :: model   ! kernel identifier
!
      integer                           :: LBi, UBi, LBj, UBj
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: IstrR, IendR, JstrR, JendR
      integer                           :: i, j, k
!
      real (r8),            allocatable :: ad_Urho(:,:), ad_Vrho(:,:)
      real (r8)                         :: adfac, adfac1, adfac2
!
!  Initialize.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      IstrR=BOUNDS(ng)%IstrR(tile)
      IendR=BOUNDS(ng)%IendR(tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JendR=BOUNDS(ng)%JendR(tile)
!
#   ifdef NESTING
      IminS=Istr-4
      ImaxS=Iend+3
      JminS=Jstr-4
      JmaxS=Jend+3
#   else
      IminS=Istr-3
      ImaxS=Iend+3
      JminS=Jstr-3
      JmaxS=Jend+3
#   endif
!
      allocate ( ad_Urho(IminS:ImaxS, JminS:JmaxS) )
      ad_Urho=0.0_r8
      allocate ( ad_Vrho(IminS:ImaxS, JminS:JmaxS) )
      ad_Vrho=0.0_r8
!
      adfac=0.0_r8
      adfac1=0.0_r8
      adfac2=0.0_r8

#   ifdef DISTRIBUTE
!
!  Adjoint of exchange boundary data.
!
!>    CALL mp_exchange3d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    self%Uc, self%Vc)
!>
      CALL ad_mp_exchange3d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_self%Uc, ad_self%Vc)
#   endif
!
!  Adjoint of variable change from A-grid to C-grid.
!
      K_LOOP :DO k=1,N(ng)
!
!  Adjoint of compute staggered C-grid components.
!
        DO j=Jstr,JendR
          DO i=IstrR,IendR
#   ifdef MASKING
!>          self%Vc(i,j,k)=self%Vc(i,j,k)*GRID(ng)%vmask(i,j)
!>
            ad_self%Vc(i,j,k)=ad_self%Vc(i,j,k)*GRID(ng)%vmask(i,j)
#   endif
!>          self%Vc(i,j,k)=0.5_r8*(Vrho(i,j-1)+Vrho(i,j))
!>
            adfac=0.5_r8*ad_self%Vc(i,j,k)
            ad_Vrho(i,j-1)=ad_Vrho(i,j-1)+adfac
            ad_Vrho(i,j  )=ad_Vrho(i,j  )+adfac
            ad_self%Vc(i,j,k)=0.0_r8
          END DO
        END DO
!
        DO j=JstrR,JendR
          DO i=Istr,IendR
#   ifdef MASKING
!>          self%Uc(i,j,k)=self%Uc(i,j,k)*GRID(ng)%umask(i,j)
!>
            ad_self%Uc(i,j,k)=ad_self%Uc(i,j,k)*GRID(ng)%umask(i,j)
#   endif
!>          self%Uc(i,j,k)=0.5_r8*(Urho(i-1,j)+Urho(i,j))
!>
            adfac=0.5_r8*ad_self%Uc(i,j,k)
            ad_Urho(i-1,j)=ad_Urho(i-1,j)+adfac
            ad_Urho(i  ,j)=ad_Urho(i  ,j)+adfac
            ad_self%Uc(i,j,k)=0.0_r8
          END DO
        END DO
!
!  Adjoint of rotate vector components to computations (XI,ETA)
!  directions.
!
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
!>          Vrho(i,j)=self%Va(i,j,k)*GRID(ng)%CosAngler(i,j)-           &
!>   &                self%Ua(i,j,k)*GRID(ng)%SinAngler(i,j)
!>
            adfac1=GRID(ng)%CosAngler(i,j)*ad_Vrho(i,j)
            adfac2=GRID(ng)%SinAngler(i,j)*ad_Vrho(i,j)
            ad_self%Va(i,j,k)=ad_self%Va(i,j,k)+adfac1
            ad_self%Ua(i,j,k)=ad_self%Ua(i,j,k)-adfac2
            ad_Vrho(i,j)=0.0_r8
!>          Urho(i,j)=self%Ua(i,j,k)*GRID(ng)%CosAngler(i,j)+           &
!>   &                self%Va(i,j,k)*GRID(ng)%SinAngler(i,j)
!>
            adfac1=GRID(ng)%CosAngler(i,j)*ad_Urho(i,j)
            adfac2=GRID(ng)%SinAngler(i,j)*ad_Urho(i,j)
            ad_self%Ua(i,j,k)=ad_self%Ua(i,j,k)+adfac1
            ad_self%Va(i,j,k)=ad_self%Va(i,j,k)+adfac2
            ad_Urho(i,j)=0.0_r8
          END DO
        END DO
      END DO K_LOOP
!
!  Deallocate local arrays.
!
      IF (allocated(ad_Urho)) deallocate (ad_Urho)
      IF (allocated(ad_Vrho)) deallocate (ad_Vrho)
!
      RETURN
      END SUBROUTINE var_change_uv_a2c_ad
!
!-----------------------------------------------------------------------
!  Adjoint of transforming vector components from staggered Arakawa
!  C-grid to A-grid (cell center).
!
      SUBROUTINE var_change_uv_c2a_ad (ad_self, ng, tile, model)
!
      CLASS (obs_VarCha), intent(inout) :: ad_self ! VarChan object
      integer,            intent(in   ) :: ng      ! nested grid
      integer,            intent(in   ) :: tile    ! domain partition
      integer,            intent(in   ) :: model   ! kernel identifier
!
      integer                           :: LBi, UBi, LBj, UBj
      integer                           :: Istr, Iend, Jstr, Jend
      integer                           :: IstrR, IendR, JstrE, JendR
      integer                           :: i, j, k
!
      real (r8)                         :: adfac
      real (r8)                         :: ad_UaCosA, ad_UaSinA
      real (r8)                         :: ad_VaCosA, ad_VaSinA
!
!  Initialize.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      Istr=BOUNDS(ng)%Istr(tile)
      Iend=BOUNDS(ng)%Iend(tile)
      Jstr=BOUNDS(ng)%Jstr(tile)
      Jend=BOUNDS(ng)%Jend(tile)
!
      IstrR=BOUNDS(ng)%IstrR(tile)
      IendR=BOUNDS(ng)%IendR(tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JendR=BOUNDS(ng)%JendR(tile)
!
      adfac=0.0_r8
      ad_UaCosA=0.0_r8
      ad_UaSinA=0.0_r8
      ad_VaCosA=0.0_r8
      ad_VaSinA=0.0_r8

#   ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
!>    CALL mp_exchange3d (ng, tile, model, 2,                           &
!>   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!>   &                    NghostPoints,                                 &
!>   &                    EWperiodic(ng), NSperiodic(ng),               &
!>   &                    self%Ua, self%Va)
!>
      CALL ad_mp_exchange3d (ng, tile, model, 2,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_self%Ua, ad_self%Va)
#   endif
!
!  Adjoint of transform vector components from staggered Arakawa C-grid
!  to A-grid (cell center).
!
      K_LOOP : DO k=1,N(ng)
!
!  Adjoint of rotate from computational to gegraphical Eastward and
!  Northward directions.
!
        DO j=JstrR,JendR
          DO i=IstrR,IendR
#   ifdef MASKING
!>          self%Va(i,j,k)=self%Va(i,j,k)*GRID(ng)%rmask(i,j)
!>
            ad_self%Va(i,j,k)=ad_self%Va(i,j,k)*GRID(ng)%rmask(i,j)
!>          self%Ua(i,j,k)=self%Ua(i,j,k)*GRID(ng)%rmask(i,j)
!>
            ad_self%Ua(i,j,k)=ad_self%Ua(i,j,k)*GRID(ng)%rmask(i,j)
#   endif
!>          seld%Va(i,j,k)=VaCosA+UaSinA
!>
            ad_UaSinA=ad_UaSinA+ad_self%Va(i,j,k)
            ad_VaCosA=ad_VaCosA+ad_self%Va(i,j,k)
            ad_self%Va(i,j,k)=0.0_r8
!>          self%Ua(i,j,k)=UaCosA-VaSinA
!>
            ad_VaSinA=ad_VaSinA-ad_self%Ua(i,j,k)
            ad_UaCosA=ad_UaCosA+ad_self%Ua(i,j,k)
            ad_self%Ua(i,j,k)=0.0_r8
!>          VaSinA=self%Va(i,j,k)*GRID(ng)%SinAngler(i,j)
!>
            ad_self%Va(i,j,k)=ad_self%Va(i,j,k)+                        &
     &                        GRID(ng)%SinAngler(i,j)*ad_VaSinA
            ad_VaSinA=0.0_r8

!>          VaCosA=self%Va(i,j,k)*GRID(ng)%CosAngler(i,j)
!>
            ad_self%Va(i,j,k)=ad_self%Va(i,j,k)+                        &
     &                        GRID(ng)%CosAngler(i,j)*ad_VaCosA
            ad_VaCosA=0.0_r8
!>          UaSinA=self%Ua(i,j,k)*GRID(ng)%SinAngler(i,j)
!>
            ad_self%Ua(i,j,k)=ad_self%Ua(i,j,k)+                        &
     &                        GRID(ng)%SinAngler(i,j)*ad_UaSinA
            ad_UaSinA=0.0_r8
!>          UaCosA=self%Ua(i,j,k)*GRID(ng)%CosAngler(i,j)
!>
            ad_self%Ua(i,j,k)=ad_self%Ua(i,j,k)+                        &
     &                        GRID(ng)%CosAngler(i,j)*ad_UaCosA
            ad_UaCosA=0.0_r8
          END DO
       END DO
!
!  Adjoint of compute A-grid (cell center) vector components.
!
       DO j=Jstr,Jend
          DO i=IstrR,IendR
            IF (.not.NSperiodic(ng)) THEN
              IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
!>              self%Va(i,Jend+1,k)=self%Va(i,Jend,k)
!>
                ad_self%Va(i,Jend,k)=ad_self%Va(i,Jend  ,k)+            &
     &                               ad_self%Va(i,Jend+1,k)
                ad_self%Va(i,Jend+1,k)=0.0_r8
              END IF
              IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
!>              self%Va(i,Jstr-1,k)=self%Va(i,Jstr,k)
!>
                ad_self%Va(i,Jstr,k)=ad_self%Va(i,Jstr  ,k)+            &
     &                               ad_self%Va(i,Jstr-1,k)
                ad_self%Va(i,Jstr-1,k)=0.0_r8
              END IF
            END IF
!>          self%Va(i,j,k)=0.5_r8*(self%Vc(i,j,k)+self%Vc(i,j+1,k))
!>
            adfac=0.5_r8*ad_self%Va(i,j,k)
            ad_self%Vc(i,j  ,k)=ad_self%Vc(i,j  ,k)+adfac
            ad_self%Vc(i,j+1,k)=ad_self%Vc(i,j+1,k)+adfac
            ad_self%Va(i,j,k)=0.0_r8
          END DO
        END DO
!
        DO j=JstrR,JendR
          DO i=Istr,Iend
            IF (.not.EWperiodic(ng)) THEN
              IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
!>              self%Ua(Iend+1,j,k)=self%Ua(Iend,j,k)
!>
                ad_self%Ua(Iend,j,k)=ad_self%Ua(Iend  ,j,k)+            &
     &                               ad_self%Ua(Iend+1,j,k)
                ad_self%Ua(Iend+1,j,k)=0.0_r8
              END IF
              IF (DOMAIN(ng)%Western_Edge(tile)) THEN
!>              self%Ua(Istr-1,j,k)=self%Ua(Istr,j,k)
!>
                ad_self%Ua(Istr,j,k)=ad_self%Ua(Istr  ,j,k)+            &
     &                               ad_self%Ua(Istr-1,j,k)
                ad_self%Ua(Istr-1,j,k)=0.0_r8
              END IF
            END IF
!>          self%Ua(i,j,k)=0.5_r8*(self%Uc(i,j,k)+self%Uc(i+1,j,k))
!>
            adfac=0.5_r8*ad_self%Ua(i,j,k)
            ad_self%Uc(i  ,j,k)=ad_self%Uc(i  ,j,k)+adfac
            ad_self%Uc(i+1,j,k)=ad_self%Uc(i+1,j,k)+adfac
            ad_self%Ua(i,j,k)=0.0_r8
          END DO
        END DO
      END DO K_LOOP
!
      RETURN
      END SUBROUTINE var_change_uv_c2a_ad
#  endif
# endif
!
!-----------------------------------------------------------------------
!
#endif
      END MODULE roms_vchange_mod
