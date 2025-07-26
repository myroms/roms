#include "cppdefs.h"
      MODULE roms_hofx_mod

#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  It contains several routines to compute the H(x) operator that      !
!  interpolates the observartion to the nonlinear (NLM) and linear     !
!  (TLM, RPM, ADM) model locations.                                    !
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
      USE distribute_mod,   ONLY : mp_collect
# endif
# if defined ADJOINT && defined DISTRIBUTE
      USE mp_exchange_mod,  ONLY : ad_mp_exchange2d, ad_mp_exchange3d
# endif
      USE strings_mod,      ONLY : FoundError
      USE tracer_index_mod, ONLY : tracer_index
!
      implicit none
!
!-----------------------------------------------------------------------
!  Observation Operator HofX Object: CLASS(obs_hofx).
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: obs_hofx
!
!  Total number of aggregated observations.

        integer :: ndatum
!
!  Kernel variables lower and upper bound array dimensions, and parallel
!  tile partition indices.
!
        integer :: tile                             ! domain partition
        integer :: LBi,   UBi,   LBj,   UBj         ! array bounds
        integer :: IstrR, IendR, JstrR, JendR       ! partition bounds
!
!  Contol vector to store state variables at the aggregated observation
!  locations.
!
        real(r8), allocatable :: state(:)           ! state vector
        integer,  allocatable :: ObsVetting(:)      ! reject/accept flag
!
!  Pointers for geometry variables.
!
        real(r8), pointer :: angle(:,:)   => NULL() ! curvilinear angle
        real(r8), pointer :: mask(:,:)    => NULL() ! land-sea mask
# ifdef SOLVE3D
        real(r8), pointer :: depth(:,:,:) => NULL() ! depths (m)
# endif
!
        CONTAINS
!
!  Constructors and destructors.
!
        PROCEDURE :: create           => obs_hofx_create
        PROCEDURE :: destroy          => obs_hofx_destroy
!
!  Operations.
!
        PROCEDURE :: add              => obs_hofx_add
        PROCEDURE :: extract          => obs_hofx_extract
        PROCEDURE :: interp2d         => obs_hofx_interp2d
        PROCEDURE :: interp2d_ad      => obs_hofx_interp2d_ad
        PROCEDURE :: interp3d         => obs_hofx_interp3d
        PROCEDURE :: interp3d_ad      => obs_hofx_interp3d_ad
        PROCEDURE :: nc_write         => obs_hofx_nc_write
        PROCEDURE :: rms              => obs_hofx_rms
        PROCEDURE :: update           => obs_hofx_update
        PROCEDURE :: zeros            => obs_hofx_zeros
        
      END TYPE obs_hofx
!
!-----------------------------------------------------------------------
!  Module variables.
!-----------------------------------------------------------------------
!
!  Background error (standard deviation) identifier.

      integer, parameter :: iSTD = 10
!
!  H(x) operators per nested grid, [1:Ngrids].
!
      TYPE (obs_hofx), allocatable :: nlm_hofx(:)     ! nonlinear
# ifdef ADJOINT
      TYPE (obs_hofx), allocatable :: adm_hofx(:)     ! adjoint
# endif
# ifndef VERIFICATION
      TYPE (obs_hofx), allocatable :: berr_hofx(:)    ! background error
# endif
# if defined TANGENT || defined TL_IOMS
      TYPE (obs_hofx), allocatable :: tlm_hofx(:)     ! tangent linear
# endif
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      PRIVATE
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      CONTAINS
!
! <><><><><><><><><><><><><><><><><><><><><><><><> CLASS OBS_HOFX <><><>
!
      SUBROUTINE obs_hofx_create (self, ng, tile, model, ntotal)
!
      CLASS (obs_hofx), intent(inout) :: self     ! H(x) object
      integer,          intent(in   ) :: ng       ! nested grid number
      integer,          intent(in   ) :: tile     ! domain partition
      integer,          intent(in   ) :: model    ! kernel identifier
      integer,          intent(in   ) :: ntotal   ! total observations
!
!  Initialize scalar variables.
!
      self%ndatum=ntotal
      self%tile=tile
!
      self%LBi=BOUNDS(ng)%LBi(tile)
      self%UBi=BOUNDS(ng)%UBi(tile)
      self%LBj=BOUNDS(ng)%LBj(tile)
      self%UBj=BOUNDS(ng)%UBj(tile)
!
      self%IstrR=BOUNDS(ng)%IstrR(tile)
      self%IendR=BOUNDS(ng)%IendR(tile)
      self%JstrR=BOUNDS(ng)%JstrR(tile)
      self%JendR=BOUNDS(ng)%JendR(tile)
!
!  Allocate contol state vector.
!
      IF (.not.allocated(self%state)) THEN
        allocate ( self%state(ntotal) )
      END IF
!
!  Allocate processing flag to reject (zero) or accept (unity) 
!  observations.
!
      IF (.not.allocated(self%ObsVetting)) THEN
        allocate ( self%ObsVetting(ntotal) )
      END IF
!
!  Set geometry arrays pointers.
!
      self%angle => GRID(ng)%angler
      self%mask  => GRID(ng)%rmask
# ifdef SOLVE3D
      self%depth => GRID(ng)%z_r
# endif
!
!
!  Initialize H(x) vector to zero to facilitate collective operatations
!  in distributed-memory.
!
      CALL self%zeros ()
!
      RETURN
      END SUBROUTINE obs_hofx_create
!
!-----------------------------------------------------------------------
!
!  It destroys H(x) operator object.
!
      SUBROUTINE obs_hofx_destroy (self)
!
      CLASS (obs_hofx), intent(inout) :: self     ! H(x) object
!
!  Deallocate H(x) object variables.
!
      IF (allocated(self%state))      deallocate (self%state)
      IF (allocated(self%ObsVetting)) deallocate (self%ObsVetting)
!
!  Nullify pointers.
!
      IF (associated(self%angle)) nullify (self%angle)
      IF (associated(self%mask))  nullify (self%mask)
# ifdef SOLVE3D
      IF (associated(self%depth)) nullify (self%depth)
# endif
!
      RETURN
      END SUBROUTINE obs_hofx_destroy
!
!-----------------------------------------------------------------------
!
!  It adds two state vectors: self = self + c * rhs.
!
      SUBROUTINE obs_hofx_add (self, c, rhs)
!
      CLASS (obs_hofx), intent(inout) :: self     ! LHS H(x) object
      CLASS (obs_hofx), intent(in   ) :: rhs      ! RHS H(x) object
      real (r8),        intent(in   ) :: c        ! scale factor
!
      integer                         :: i
!
!  Add states vectors.
!
      DO i=1,SIZE(self%state)
        self%state(i)=self%state(i)+rhs%state(i)
      END DO
!
      RETURN
      END SUBROUTINE obs_hofx_add
!
!-----------------------------------------------------------------------
!
!  It computes the H(x) operator, which interpolates model variables to
!  observation locations in space and time. If adjoint model kernel, it
!  computes transpose(H)x.
!
      SUBROUTINE obs_hofx_extract (self, obs, ng, model)
!
      CLASS (obs_hofx), intent(inout) :: self     ! H(x) object
      TYPE (obs_pool),  intent(inout) :: obs      ! observation object
      integer,          intent(in   ) :: ng       ! nested grid number
      integer,          intent(in   ) :: model    ! kernel identifier
!
      TYPE (obs_data),        pointer :: depot
      integer                         :: TrcIndex, iobs, ivar
      character (len=40)              :: pool_name

      character (len=*),    parameter :: MyFile =                       &
     &    __FILE__//", roms_hofx_extract"
!
      SourceFile=MyFile
!
!  If observations are available for the (time-0.5*dt) to (time+0.5*dt)
!  window, interpolate model variables to the observation locations.
!  The "ivar" is a counter for "obs%pool(iobs)%nvars", and the logic
!  below must be generalized.
!
      IF (obs%has(ng, model)) THEN
        DO iobs=1,SIZE(obs%pool)
          IF ((obs%pool(iobs)%NstrObs.gt.0).and.                        &
     &        (obs%pool(iobs)%NendObs.gt.0)) THEN
            pool_name=obs%pool(iobs)%name
            SELECT CASE (TRIM(pool_name))
              CASE ('SSH')
                ivar=1
                CALL obs_pool_get (obs, TRIM(pool_name), depot) 
                IF (model.ne.iADM) THEN
                  CALL self%interp2d (depot, 'SSH', ng, model, ivar,    &
     &                                KOUT)
# ifdef ADJOINT
                ELSE
                  CALL self%interp2d_ad (depot, 'SSH', ng, model, ivar)
# endif
                END IF
# ifdef SOLVE3D
              CASE ('SSS', 'salt', 'SST', 'temp')
                ivar=1
                CALL obs_pool_get (obs, TRIM(pool_name), depot) 
                TrcIndex=tracer_index(obs%pool(iobs)%name)
                IF (model.ne.iADM) THEN
                  CALL self%interp3d (depot, 'tracer', ng, model, ivar, &
     &                                NOUT,                             &
     &                                itrc = TrcIndex)
#  ifdef ADJOINT
                ELSE
                  CALL self%interp3d_ad (depot, 'tracer', ng, model,    &
     &                                   ivar,                          &
     &                                   itrc = TrcIndex)
#  endif
                END IF
              CASE ('uv')
                CALL obs_pool_get (obs, TRIM(pool_name), depot) 
                ivar=1
                IF (model.ne.iADM) THEN
                  CALL self%interp3d (depot, 'u', ng, model, ivar,      &
     &                                NOUT)
#  ifdef ADJOINT
                ELSE
                  CALL self%interp3d_ad (depot, 'u', ng, model, ivar)
#  endif
                END IF
!
                ivar=2
                CALL obs_pool_get (obs, TRIM(pool_name), depot) 
                IF (model.ne.iADM) THEN
                  CALL self%interp3d (depot, 'v', ng, model, ivar,      &
     &                                NOUT)
#  ifdef ADJOINT
                ELSE
                  CALL self%interp3d_ad (depot, 'v', ng, model, ivar)
#  endif
                END IF
# endif
              CASE DEFAULT
                IF (Master) WRITE (stdout,10) TRIM(pool_name), iobs
                exit_flag=5
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
            END SELECT
          END IF
        END DO
      END IF
!
  10  FORMAT (/,' ROMS_HOFX::extract: Cannot find name = ',a,2x,        &
     &        'in obs%pool(',i0,').')
!
      RETURN
      END SUBROUTINE obs_hofx_extract
!
!-----------------------------------------------------------------------
!
!  It interpolates 2D state variables at the observation locations.
!  Observations are assumed to be at RHO-points (A-grid). The staggering
!  of state variables is solely done for numerical accuracy and to
!  suppress spurious modes.
!
!  The interpolation weights matrix, Hmat(1:4), is as follows:
!
!                             4____________3  (i2,j2)
!                             |            |
!                             | * (p,q)    |
!                             |            |
!                             |____________|
!                    (i1,j1)  1            2
!
      SUBROUTINE obs_hofx_interp2d (self, obs, name, ng, model, ivar,   &
     &                              Tindex)
!
      CLASS (obs_hofx),  intent(inout) :: self     ! H(x) object
      TYPE  (obs_data),  intent(in   ) :: obs      ! observation store
      character (len=*), intent(in   ) :: name     ! state variable name
      integer,           intent(in   ) :: ng       ! nested grid number
      integer,           intent(in   ) :: model    ! kernel identifier
      integer,           intent(in   ) :: ivar     ! variable counter
      integer,           intent(in   ) :: Tindex   ! time level index
!
      integer                         :: i, ic, iobs, i1, i2, j1, j2
      real (r8)                       :: p1, p2, q1, q2, wsum
      real (r8),         dimension(4) :: Hmat
      real (r8),              pointer :: field(:,:)

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", roms_hofx_interp2d"
!
      SourceFile=MyFile
!
!  Set pointer to kernel variable.
!
      SELECT CASE (TRIM(name))
        CASE ('ssh', 'SSH')
          IF (model.eq.iNLM) THEN
            field => OCEAN(ng)%zeta(:,:,Tindex)
# ifdef TLM_OBS
          ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
            field => OCEAN(ng)%tl_zeta(:,:,Tindex)
# endif
# ifndef VERIFICATION
          ELSE IF (model.eq.iSTD) THEN
            field => OCEAN(ng)%e_zeta(:,:,1)          ! NSA=1 array
# endif
          END IF
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) TRIM(obs%name)
          exit_flag=5
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END SELECT
!
!  Interpolate 2D model state variable at the observation locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
            (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng))) THEN
!
          iobs=obs%pack_mapping(i,ivar)               ! index mapping
!
          i1=INT(obs%Xgrid(i))
          j1=INT(obs%Ygrid(i))
          i2=i1+1
          j2=j1+1
          IF (i2.gt.Lm(ng)+1) THEN
            i2=i1                 ! Observation at the eastern boundary
          END IF
          IF (j2.gt.Mm(ng)+1) THEN
            j2=j1                 ! Observation at the northern boundary
          END IF
          p2=REAL(i2-i1,r8)*(obs%Xgrid(i)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(obs%Ygrid(i)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          Hmat(1)=p1*q1
          Hmat(2)=p2*q1
          Hmat(3)=p2*q2
          Hmat(4)=p1*q2
# ifdef MASKING
          Hmat(1)=Hmat(1)*self%mask(i1,j1)
          Hmat(2)=Hmat(2)*self%mask(i2,j1)
          Hmat(3)=Hmat(3)*self%mask(i2,j2)
          Hmat(4)=Hmat(4)*self%mask(i1,j2)
          wsum=0.0_r8
          DO ic=1,4
            wsum=wsum+Hmat(ic)
          END DO
          IF (wsum.gt.0.0_r8) THEN
            wsum=1.0_r8/wsum
            DO ic=1,4
              Hmat(ic)=Hmat(ic)*wsum
            END DO
          END IF
# endif
          self%state(iobs)=Hmat(1)*field(i1,j1)+                        &
     &                     Hmat(2)*field(i2,j1)+                        &
     &                     Hmat(3)*field(i2,j2)+                        &
     &                     Hmat(4)*field(i1,j2)
# ifdef MASKING
          IF (wsum.gt.0.0_r8) self%obsVetting(iobs)=1
# else
          self%ObsVetting(iobs)=1
# endif
        END IF
      END DO
!
!  Nullify pointer.
!
      IF (associated(field)) nullify (field)
!
  10  FORMAT (/,' ROMS_HOFX::interp2d: Cannot find name = ',a,2x,        &
     &        'in obs%data observation depot object.')
!
      RETURN
      END SUBROUTINE obs_hofx_interp2d

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!
!  It interpolates 3D state variables at the observation locations.
!  Observations are assumed to be at RHO-points (A-grid). The staggering
!  of state variables is solely done for numerical accuracy and to
!  suppress spurious modes.
!
!  The interpolation weights matrix, Hmat(1:8), is as follows:
!
!                               8______________7
!                               /.            /| (i2,j2,k2)
!                              / .           / |
!                            5/_____________/6 |
!                             |  .          |  |
!                             |  . *(p,q,r) |  |         Grid Cell
!                             | 4...........|..|3
!                             | .           |  /
!                             |.            | /
!                  (i1,j1,k1) |_____________|/
!                             1           2
!
      SUBROUTINE obs_hofx_interp3d (self, obs, name, ng, model, ivar,   &
     &                              Tindex, itrc)
!
      CLASS (obs_hofx),  intent(inout) :: self     ! H(x) object
      TYPE  (obs_data),  intent(inout) :: obs      ! observation store
      character (len=*), intent(in   ) :: name     ! state variable name
      integer,           intent(in   ) :: ng       ! nested grid number
      integer,           intent(in   ) :: model    ! kernel identifier
      integer,           intent(in   ) :: ivar     ! variable counter
      integer,           intent(in   ) :: Tindex   ! time level index
      integer, optional, intent(in   ) :: itrc     ! tracer index
!
      logical                          :: IsVector
      integer                          :: i, ic, iobs, k
      integer                          :: i1, i2, j1, j2, k1, k2
      real (r8)                        :: p1, p2, q1, q2, r1, r2
      real (r8)                        :: w11, w12, w21, w22, wsum
      real (r8)                        :: Zbot, Ztop, dz
      real (r8),          dimension(8) :: Hmat
      real (r8),               pointer :: field(:,:,:) => NULL()

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", roms_hofx_interp3d"
!
      SourceFile=MyFile
!
!  Set pointer to kernel variable.
!
      IsVector=.FALSE.
!
      SELECT CASE (TRIM(name))
        CASE ('tracer', 'TRACER')
          IF (model.eq.iNLM) THEN
            field => OCEAN(ng)%t(:,:,:,itrc,Tindex)
#  ifdef TLM_OBS
          ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
            field => OCEAN(ng)%tl_t(:,:,:,itrc,Tindex)
#  endif
#  ifndef VERIFICATION
          ELSE IF (model.eq.iSTD) THEN
            field => OCEAN(ng)%e_t(:,:,:,itrc,1)      ! NSA=1 array
#  endif
          END IF
        CASE ('ua', 'Ua')
          IF (model.eq.iNLM) THEN
            field => OCEAN(ng)%ua
#  ifdef TLM_OBS
          ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
            field => OCEAN(ng)%tl_ua
#  endif
#  ifndef VERIFICATION
          ELSE IF (model.eq.iSTD) THEN
            field => OCEAN(ng)%e_u(:,:,:,1)           ! NSA=1 array
#  endif
          END IF
        CASE ('va', 'Va')
          IF (model.eq.iNLM) THEN
            field => OCEAN(ng)%va
#  ifdef TLM_OBS
          ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
            field => OCEAN(ng)%tl_va
#  endif
#  ifndef VERIFICATION
          ELSE IF (model.eq.iSTD) THEN
            field => OCEAN(ng)%e_v(:,:,:,1)           ! NSA=1 array
#  endif
          END IF
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) TRIM(obs%name)
          exit_flag=5
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END SELECT
!
!  Interpolate 2D model state variable at the observation locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
            (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng))) THEN
!
          iobs=obs%pack_mapping(i,ivar)               ! index mapping
!
          i1=INT(obs%Xgrid(i))
          j1=INT(obs%Ygrid(i))
          i2=i1+1
          j2=j1+1
          IF (i2.gt.Lm(ng)+1) THEN
            i2=i1                 ! Observation at the eastern boundary
          END IF
          IF (j2.gt.Mm(ng)+1) THEN
            j2=j1                 ! Observation at the northern boundary
          END IF
          p2=REAL(i2-i1,r8)*(obs%Xgrid(i)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(obs%Ygrid(i)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          w11=p1*q1
          w21=p2*q1
          w22=p2*q2
          w12=p1*q2
          IF (obs%Zgrid(i).gt.0.0_r8) THEN
            k1=MAX(1,INT(obs%Zgrid(i)))        ! Positions in fractional
            k2=MIN(INT(obs%Zgrid(i))+1,N(ng))  ! levels
            r2=REAL(k2-k1,r8)*(obs%Zgrid(i)-REAL(k1,r8))
            r1=1.0_r8-r2
          ELSE
            Ztop=self%depth(i1,j1,N(ng))
            Zbot=self%depth(i1,j1,1    )
            IF (obs%Zgrid(i).ge.Ztop) THEN
              k1=N(ng)                         ! If shallower, assign to
              k2=N(ng)                         ! top grid cell. The
              r1=1.0_r8                        ! observation is located
              r2=0.0_r8                        ! on the upper cell half
              obs%Zgrid(i)=REAL(N(ng),r8)      ! above its middle depth.
            ELSE IF (Zbot.ge.obs%Zgrid(i)) THEN
              r1=0.0_r8                        ! If deeper, ignore.
              r2=0.0_r8
              self%ObsVetting(iobs)=0
            ELSE
              DO k=N(ng),2,-1                  ! Otherwise, interpolate
                Ztop=self%depth(i1,j1,k  )     ! to fractional level
                Zbot=self%depth(i1,j1,k-1)
                IF ((Ztop.gt.obs%Zgrid(i)).and.                         &
     &              (obs%Zgrid(i).ge.Zbot)) THEN
                  k1=k-1
                  k2=k
                END IF
              END DO
              dz=self%depth(i1,j1,k2)-self%depth(i1,j1,k1)
              r2=(obs%Zgrid(i)-self%depth(i1,j1,k1))/dz
              r1=1.0_r8-r2
              obs%Zgrid(i)=REAL(k1,r8)+r2      ! overwrite
            END IF
          END IF
          IF ((r1+r2).gt.0.0_r8) THEN
            Hmat(1)=w11*r1
            Hmat(2)=w21*r1
            Hmat(3)=w22*r1
            Hmat(4)=w12*r1
            Hmat(5)=w11*r2
            Hmat(6)=w21*r2
            Hmat(7)=w22*r2
            Hmat(8)=w12*r2
#  ifdef MASKING
            Hmat(1)=Hmat(1)*self%mask(i1,j1)
            Hmat(2)=Hmat(2)*self%mask(i2,j1)
            Hmat(3)=Hmat(3)*self%mask(i2,j2)
            Hmat(4)=Hmat(4)*self%mask(i1,j2)
            Hmat(5)=Hmat(5)*self%mask(i1,j1)
            Hmat(6)=Hmat(6)*self%mask(i2,j1)
            Hmat(7)=Hmat(7)*self%mask(i2,j2)
            Hmat(8)=Hmat(8)*self%mask(i1,j2)
            wsum=0.0_r8
            DO ic=1,8
              wsum=wsum+Hmat(ic)
            END DO
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              DO ic=1,8
                Hmat(ic)=Hmat(ic)*wsum
              END DO
            END IF
#  endif
            self%state(iobs)=Hmat(1)*field(i1,j1,k1)+                  &
     &                       Hmat(2)*field(i2,j1,k1)+                  &
     &                       Hmat(3)*field(i2,j2,k1)+                  &
     &                       Hmat(4)*field(i1,j2,k1)+                  &
     &                       Hmat(5)*field(i1,j1,k2)+                  &
     &                       Hmat(6)*field(i2,j1,k2)+                  &
     &                       Hmat(7)*field(i2,j2,k2)+                  &
     &                       Hmat(8)*field(i1,j2,k2)
#  ifdef MASKING
            IF (wsum.gt.0.0_r8) self%ObsVetting(iobs)=1
#  else
            self%ObsVetting(iobs)=1
#  endif
#  ifndef ALLOW_BOTTOM_OBS
!
!  Reject observations that lie in the lower bottom grid cell (k=1) to
!  avoid clustering due shallowing of bathymetry during smoothing and
!  coarse level half-thickness (-h < Zobs < self%depth(:,:,1)) in deep
!  water.
!
            IF ((obs%Zgrid(i).gt.0.0_r8).and.                           &
     &          (obs%Zgrid(i).le.1.0_r8)) THEN
              self%ObsVetting(iobs)=0
            END IF
#  endif
          END IF
        END IF
      END DO
!
!  Nullify pointer.
!
      IF (associated(field)) nullify (field)
!
  10  FORMAT (/,' ROMS_HOFX::interp3d: Cannot find name = ',a,2x,        &
     &        'in obs%data observation depot object.')
!
      RETURN
      END SUBROUTINE obs_hofx_interp3d
# endif

# ifdef ADJOINT
!
!-----------------------------------------------------------------------
!
!  Adjoint od interpolating 2D state variables at the observation
!  locations. Observations are assumed to be at RHO-points (A-grid).
!  The staggering of state variables is solely done for numerical
!  accuracy and to suppress spurious modes.
!
!  The interpolation weights matrix, Hmat(1:4), is as follows:
!
!                             4____________3  (i2,j2)
!                             |            |
!                             | * (p,q)    |
!                             |            |
!                             |____________|
!                    (i1,j1)  1            2
!
      SUBROUTINE obs_hofx_interp2d_ad (self, obs, name, ng, model, ivar)
!
      CLASS (obs_hofx),  intent(inout) :: self     ! H(x) object
      TYPE  (obs_data),  intent(in   ) :: obs      ! observation store
      character (len=*), intent(in   ) :: name     ! state variable name
      integer,           intent(in   ) :: ng       ! nested grid number
      integer,           intent(in   ) :: model    ! kernel identifier
      integer,           intent(in   ) :: ivar     ! variable counter
!
      integer                         :: i, ic, iobs, i1, i2, j1, j2
      real (r8)                       :: p1, p2, q1, q2, wsum
      real (r8),         dimension(4) :: Hmat
      real (r8),              pointer :: ad_field(:,:)

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", roms_hofx_interp2d_ad"
!
      SourceFile=MyFile
!
!  Set pointer to kernel variable.
!
      SELECT CASE (TRIM(name))
        CASE ('ssh', 'SSH')
          ad_field => OCEAN(ng)%f_zeta
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) TRIM(obs%name)
          exit_flag=5
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END SELECT
!
!  Interpolate 2D model state variable at the observation locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
            (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng))) THEN
!
          iobs=obs%pack_mapping(i,ivar)               ! index mapping
!
          i1=INT(obs%Xgrid(i))
          j1=INT(obs%Ygrid(i))
          i2=i1+1
          j2=j1+1
          IF (i2.gt.Lm(ng)+1) THEN
            i2=i1                  ! Observation at the eastern boundary
          END IF
          IF (j2.gt.Mm(ng)+1) THEN
            j2=j1                  ! Observation at the northern boundary
          END IF
          p2=REAL(i2-i1,r8)*(obs%Xgrid(i)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(obs%Ygrid(i)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          Hmat(1)=p1*q1
          Hmat(2)=p2*q1
          Hmat(3)=p2*q2
          Hmat(4)=p1*q2
#  ifdef MASKING
          Hmat(1)=Hmat(1)*self%mask(i1,j1)
          Hmat(2)=Hmat(2)*self%mask(i2,j1)
          Hmat(3)=Hmat(3)*self%mask(i2,j2)
          Hmat(4)=Hmat(4)*self%mask(i1,j2)
          wsum=0.0_r8
          DO ic=1,4
            wsum=wsum+Hmat(ic)
          END DO
          IF (wsum.gt.0.0_r8) THEN
            wsum=1.0_r8/wsum
            DO ic=1,4
              Hmat(ic)=Hmat(ic)*wsum
            END DO
          END IF
#  endif
!^        self%state(iobs)=Hmat(1)*field(i1,j1)+                        &
!^   &                     Hmat(2)*field(i2,j1)+                        &
!^   &                     Hmat(3)*field(i2,j2)+                        &
!^   &                     Hmat(4)*field(i1,j2)
!^
          ad_field(i1,j1)=ad_field(i1,j1)+                              &
     &                    Hmat(1)*self%state(iobs)
          ad_field(i2,j1)=ad_field(i2,j1)+                              &
     &                    Hmat(2)*self%state(iobs)
          ad_field(i2,j2)=ad_field(i2,j2)+                              &
     &                    Hmat(3)*self%state(iobs)
          ad_field(i1,j2)=ad_field(i1,j2)+                              &
     &                    Hmat(4)*self%state(iobs)
          self%state(iobs)=0.0_r8
#  ifdef MASKING
          IF (wsum.gt.0.0_r8) self%ObsVetting(iobs)=1
#  else
          self%ObsVetting(iobs)=1
#  endif
        END IF
      END DO

#  ifdef DISTRIBUTE
!
!  Exchange tile data.
!
      CALL ad_mp_exchange2d (ng, self%tile, model, 1,                   &
     &                       self%LBi, self%UBi, self%LBj, self%UBj,    &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_field)
#  endif
!
!  Nullify pointer.
!
      IF (associated(ad_field)) nullify (ad_field)
!
  10  FORMAT (/,' ROMS_HOFX::interp2d_ad: Cannot find name = ',a,2x,    &
     &        'in obs%data observation depot object.')
!
      RETURN
      END SUBROUTINE obs_hofx_interp2d_ad

#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!
!  It interpolates 3D state variables at the observation locations.
!  Observations are assumed to be at RHO-points (A-grid). The staggering
!  of state variables is solely done for numerical accuracy and to
!  suppress spurious modes.
!
!  The interpolation weights matrix, Hmat(1:8), is as follows:
!
!                               8______________7
!                               /.            /| (i2,j2,k2)
!                              / .           / |
!                            5/_____________/6 |
!                             |  .          |  |
!                             |  . *(p,q,r) |  |         Grid Cell
!                             | 4...........|..|3
!                             | .           |  /
!                             |.            | /
!                  (i1,j1,k1) |_____________|/
!                             1           2
!
      SUBROUTINE obs_hofx_interp3d_ad (self, obs, name, ng, model,      &
     &                                 ivar, itrc)
!
      CLASS (obs_hofx),  intent(inout) :: self     ! H(x) object
      TYPE  (obs_data),  intent(inout) :: obs      ! observation store
      character (len=*), intent(in   ) :: name     ! state variable name
      integer,           intent(in   ) :: ng       ! nested grid number
      integer,           intent(in   ) :: model    ! kernel identifier
      integer,           intent(in   ) :: ivar     ! variable counter
      integer, optional, intent(in   ) :: itrc     ! tracer index
!
      integer                          :: i, ic, iobs, k
      integer                          :: i1, i2, j1, j2, k1, k2
      real (r8)                        :: p1, p2, q1, q2, r1, r2
      real (r8)                        :: w11, w12, w21, w22, wsum
      real (r8)                        :: Zbot, Ztop, dz
      real (r8),          dimension(8) :: Hmat
      real (r8),               pointer :: ad_field(:,:,:) => NULL()

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", roms_hofx_interp3d_ad"
!
      SourceFile=MyFile
!
!  Set pointer to kernel variable.
!
      SELECT CASE (TRIM(name))
        CASE ('tracer', 'TRACER')
          ad_field => OCEAN(ng)%f_t(:,:,:,itrc)
        CASE ('ua', 'Ua')
          ad_field => OCEAN(ng)%f_u
        CASE ('va', 'Va')
          ad_field => OCEAN(ng)%f_v
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) TRIM(obs%name)
          exit_flag=5
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END SELECT
!
!  Interpolate 2D model state variable at the observation locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
            (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng))) THEN
!
          iobs=obs%pack_mapping(i,ivar)               ! index mapping
!
          i1=INT(obs%Xgrid(i))
          j1=INT(obs%Ygrid(i))
          i2=i1+1
          j2=j1+1
          IF (i2.gt.Lm(ng)+1) THEN
            i2=i1                 ! Observation at the eastern boundary
          END IF
          IF (j2.gt.Mm(ng)+1) THEN
            j2=j1                 ! Observation at the northern boundary
          END IF
          p2=REAL(i2-i1,r8)*(obs%Xgrid(i)-REAL(i1,r8))
          q2=REAL(j2-j1,r8)*(obs%Ygrid(i)-REAL(j1,r8))
          p1=1.0_r8-p2
          q1=1.0_r8-q2
          w11=p1*q1
          w21=p2*q1
          w22=p2*q2
          w12=p1*q2
          IF (obs%Zgrid(i).gt.0.0_r8) THEN
            k1=MAX(1,INT(obs%Zgrid(i)))        ! Positions in fractional
            k2=MIN(k1+1,N(ng))                 ! levels
            r2=REAL(k2-k1,r8)*(obs%Zgrid(i)-REAL(k1,r8))
            r1=1.0_r8-r2
          ELSE
            Ztop=self%depth(i1,j1,N(ng))
            Zbot=self%depth(i1,j1,1    )
            IF (obs%Zgrid(i).ge.Ztop) THEN
              k1=N(ng)                         ! If shallower, assign to
              k2=N(ng)                         ! top grid cell. The
              r1=1.0_r8                        ! observation is located
              r2=0.0_r8                        ! on the upper cell half
              obs%Zgrid(i)=REAL(N(ng),r8)      ! above its middle depth.
            ELSE IF (Zbot.ge.obs%Zgrid(iobs)) THEN
              r1=0.0_r8                        ! If deeper, ignore.
              r2=0.0_r8
              self%ObsVetting(iobs)=0
            ELSE
              DO k=N(ng),2,-1                  ! Otherwise, interpolate
                Ztop=self%depth(i1,j1,k  )     ! to fractional level
                Zbot=self%depth(i1,j1,k-1)
                IF ((Ztop.gt.obs%Zgrid(i)).and.                         &
     &              (obs%Zgrid(i).ge.Zbot)) THEN
                  k1=k-1
                  k2=k
                END IF
              END DO
              dz=self%depth(i1,j1,k2)-self%depth(i1,j1,k1)
              r2=(obs%Zgrid(i)-self%depth(i1,j1,k1))/dz
              r1=1.0_r8-r2
              obs%Zgrid(i)=REAL(k1,r8)+r2      ! overwrite
            END IF
          END IF
          IF ((r1+r2).gt.0.0_r8) THEN
            Hmat(1)=w11*r1
            Hmat(2)=w21*r1
            Hmat(3)=w22*r1
            Hmat(4)=w12*r1
            Hmat(5)=w11*r2
            Hmat(6)=w21*r2
            Hmat(7)=w22*r2
            Hmat(8)=w12*r2
#   ifdef MASKING
            Hmat(1)=Hmat(1)*self%mask(i1,j1)
            Hmat(2)=Hmat(2)*self%mask(i2,j1)
            Hmat(3)=Hmat(3)*self%mask(i2,j2)
            Hmat(4)=Hmat(4)*self%mask(i1,j2)
            Hmat(5)=Hmat(5)*self%mask(i1,j1)
            Hmat(6)=Hmat(6)*self%mask(i2,j1)
            Hmat(7)=Hmat(7)*self%mask(i2,j2)
            Hmat(8)=Hmat(8)*self%mask(i1,j2)
            wsum=0.0_r8
            DO ic=1,8
              wsum=wsum+Hmat(ic)
            END DO
            IF (wsum.gt.0.0_r8) THEN
              wsum=1.0_r8/wsum
              DO ic=1,8
                Hmat(ic)=Hmat(ic)*wsum
              END DO
            END IF
#   endif
!^          self%state(iobs)=Hmat(1)*field(i1,j1,k1)+                  &
!^   &                       Hmat(2)*field(i2,j1,k1)+                  &
!^   &                       Hmat(3)*field(i2,j2,k1)+                  &
!^   &                       Hmat(4)*field(i1,j2,k1)+                  &
!^   &                       Hmat(5)*field(i1,j1,k2)+                  &
!^   &                       Hmat(6)*field(i2,j1,k2)+                  &
!^   &                       Hmat(7)*field(i2,j2,k2)+                  &
!^   &                       Hmat(8)*field(i1,j2,k2)
!^
            ad_field(i1,j1,k1)=ad_field(i1,j1,k1)+                     &
     &                         Hmat(1)*self%state(iobs)
            ad_field(i2,j1,k1)=ad_field(i2,j1,k1)+                     &
     &                         Hmat(2)*self%state(iobs)
            ad_field(i2,j2,k1)=ad_field(i2,j2,k1)+                     &
     &                         Hmat(3)*self%state(iobs)
            ad_field(i1,j2,k1)=ad_field(i1,j2,k1)+                     &
     &                         Hmat(4)*self%state(iobs)
            ad_field(i1,j1,k2)=ad_field(i1,j1,k2)+                     &
     &                         Hmat(5)*self%state(iobs)
            ad_field(i2,j1,k2)=ad_field(i2,j1,k2)+                     &
     &                         Hmat(6)*self%state(iobs)
            ad_field(i2,j2,k2)=ad_field(i2,j2,k2)+                     &
     &                         Hmat(7)*self%state(iobs)
            ad_field(i1,j2,k2)=ad_field(i1,j2,k2)+                     &
     &                         Hmat(8)*self%state(iobs)
            self%state(iobs)=0.0_r8
#   ifdef MASKING
            IF (wsum.gt.0.0_r8) self%ObsVetting(iobs)=1
#   else
            self%ObsVetting(iobs)=1
#   endif
          END IF
        END IF
      END DO

#   ifdef DISTRIBUTE
!
!  Exchange tile data.
!
      CALL ad_mp_exchange3d (ng, self%tile, model, 1,                   &
     &                       self%LBi, self%UBi, self%LBj, self%UBj,    &
     &                       1, N(ng), NghostPoints,                    &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_field)
#   endif
!
!  Nullify pointer.
!
      IF (associated(ad_field)) nullify (ad_field)
!
  10  FORMAT (/,' ROMS_HOFX::interp3d_ad: Cannot find name = ',a,2x,    &
     &        'in obs%data observation depot object.')
!
      RETURN
      END SUBROUTINE obs_hofx_interp3d_ad
#  endif
# endif
!
!-----------------------------------------------------------------------
!
!  It writes H(x) variables into single or multiple output NetCDF
!  file(s).
!
      SUBROUTINE obs_hofx_nc_write (self, ng, model, ncname)
!
      CLASS (obs_hofx),   intent(in   ) :: self   ! H(x) object
      integer,            intent(in   ) :: ng     ! nested grid number
      integer,            intent(in   ) :: model  ! kernel identifier
      character (len=*) , intent(in   ) :: ncname ! NetCDF filename
!
      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", roms_hofx_nc_write"
!
      SourceFile=MyFile
!
      RETURN
      END SUBROUTINE obs_hofx_nc_write
!
!-----------------------------------------------------------------------
!
!  It computes the root-mean-square value of the H(x) vector.
!
      SUBROUTINE obs_hofx_rms (self, rms)
!
      CLASS (obs_hofx), intent(in   ) :: self     ! H(x) object
      real (r8),        intent(out  ) :: rms      ! root-mean-square
!
      integer                         :: i
      real (r8)                       :: my_sum
!
!  Compute RMS value.
!
      my_sum=0.0_r8
      DO i=1,self%ndatum
        my_sum=my_sum+self%state(i)*self%state(i)
      END DO
!
!  Normalize by the number observations and take the square root.
!
      rms=SQRT(my_sum/REAL(self%ndatum,r8))
!
      RETURN
      END SUBROUTINE obs_hofx_rms
!
!-----------------------------------------------------------------------
!
!  It collects H(x) values from all tasks in the distributed-memory
!  group.
!
      SUBROUTINE obs_hofx_update (self, ng, model)
!
      CLASS (obs_hofx), intent(inout) :: self     ! H(x) object
      integer,          intent(in   ) :: ng       ! nested grid number
      integer,          intent(in   ) :: model    ! kernel identifier
!
      integer,              parameter :: Ivalue = 0
      real (r8),            parameter :: Rvalue = 0.0_r8

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", roms_hofx_update"
!
      SourceFile=MyFile

# ifdef DISTRIBUTE
!
!  Collect H(x) values and ObsVetting flag.
!
      CALL mp_collect (ng, model, self%ndatum, Rvalue, self%state)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      CALL mp_collect (ng, model, self%ndatum, Ivalue, self%ObsVetting)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
      RETURN
      END SUBROUTINE obs_hofx_update
!
!-----------------------------------------------------------------------
!
!  It zeroth out H(x) vector variables.
!
      SUBROUTINE obs_hofx_zeros (self)
!
      CLASS (obs_hofx),  intent(inout) :: self    ! H(x) object
!
      integer,               parameter :: Ivalue = 0
      real (r8),             parameter :: Rvalue = 0.0_r8
!
!  Initialize object vector variables to zero. In distributed-memory
!  applications, the collection of data in all tasks is donne by
!  summation.
!
      self%state = Rvalue
      self%ObsVetting = Ivalue
!
      RETURN
      END SUBROUTINE obs_hofx_zeros
!
!------------------------------------------------------------------------
!
#endif
      END MODULE roms_hofx_mod
