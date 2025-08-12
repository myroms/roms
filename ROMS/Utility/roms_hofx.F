#include "cppdefs.h"
      MODULE roms_hofx_mod

#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group            Andrew M. Moore   !
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
      USE distribute_mod,   ONLY : mp_aggregate2d,                      &
     &                             mp_collect,                          &
     &                             mp_reduce
# endif
# if defined ADJOINT && defined DISTRIBUTE
      USE mp_exchange_mod,  ONLY : ad_mp_exchange2d
#  ifdef SOLVE3D
      USE mp_exchange_mod,  ONLY : ad_mp_exchange3d
#  endif
# endif
      USE strings_mod,      ONLY : FoundError
      USE tracer_index_mod, ONLY : tracer_index
!
      implicit none
!
!-----------------------------------------------------------------------
!  Observation Filter Operator Object: CLASS(obs_filter).
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: obs_filter
!
!  Kernel variables lower and upper bound array dimensions, and parallel
!  tile partition indices.
!
        integer :: tile                             ! domain partition
        integer :: LBi,  UBi,  LBj,  UBj            ! array bounds
        integer :: Imin, Imax, Jmin, Jmax           ! global bounds
!
!  Application global fields required for spatial averaging operators.
!
        real (r8), pointer :: area(:,:) => NULL()   ! grid area
        real (r8), pointer :: dx(:,:)   => NULL()   ! X-grid spacing
        real (r8), pointer :: dy(:,:)   => NULL()   ! Y-grid spacing
# ifdef MASKING
        real (r8), pointer :: mask(:,:) => NULL()   ! land-sea mask
# endif
!
!  Generic global state variables to facilitate area averaging. The
!  area averaging is only supported for 2D fields (SSH, SSS, SST, Usur).
!
        real (r8), allocatable :: global2d(:,:)     ! 2D state variable
!
        CONTAINS
!
        PROCEDURE :: create   => obs_filter_create
        PROCEDURE :: destroy  => obs_filter_destroy

      END TYPE obs_filter
!
!-----------------------------------------------------------------------
!  Observation Operator H(x) Object: CLASS(obs_hofx).
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: obs_hofx
!
!  Total number of aggregated observations.
!
        integer :: ndatum
!
!  Kernel variables lower and upper bound array dimensions, and parallel
!  tile partition indices.
!
        integer :: tile                              ! domain partition
        integer :: LBi,   UBi,   LBj,   UBj          ! array bounds
        integer :: IstrR, IendR, JstrR, JendR        ! partition bounds
!
!  Tally of processed observations.
!
        integer (i8b) :: ObsCount                    ! approved
        integer (i8b) :: ObsReject                   ! rejected
!
!  H(x) operator name.
!
        character (len=:), allocatable :: name
!
!  Contol vector to store state variables at the aggregated observation
!  locations.
!
        real (r8), allocatable :: state(:)           ! state vector
        integer,   allocatable :: ObsVetting(:)      ! reject/accept flag
!
!  Pointers for geometry variables.
!
        real (r8), pointer :: angle(:,:)   => NULL() ! curvilinear angle
        real (r8), pointer :: mask(:,:)    => NULL() ! land-sea mask
# ifdef SOLVE3D
        real (r8), pointer :: depth(:,:,:) => NULL() ! depths (m)
# endif
!
        CONTAINS
!
!  Constructors and destructors.
!
        PROCEDURE :: create              => obs_hofx_create
        PROCEDURE :: destroy             => obs_hofx_destroy
!
!  H(x) Operators: Extract model at observation locations drivers.
!
        PROCEDURE :: extract             => obs_hofx_extract
        PROCEDURE :: extract2d           => obs_hofx_extract2d
# ifdef SOLVE3D
        PROCEDURE :: extract3d           => obs_hofx_extract3d
# endif
!
!  H(x) operators: Interpolation without filters.
!
        PROCEDURE :: interp2d            => obs_hofx_interp2d
# ifdef ADJOINT
        PROCEDURE :: interp2d_ad         => obs_hofx_interp2d_ad
# endif
# ifdef SOLVE3D
        PROCEDURE :: interp3d            => obs_hofx_interp3d
#  ifdef ADJOINT
        PROCEDURE :: interp3d_ad         => obs_hofx_interp3d_ad
#  endif
# endif
!
!  H(x) Operators: Spatial and/or temporal H(X) average filters.
!
        PROCEDURE :: area_avg2d          => obs_hofx_area_avg2d
# ifdef ADJOINT
        PROCEDURE :: area_avg2d_ad       => obs_hofx_area_avg2d_ad
# endif
        PROCEDURE :: area_time_avg2d     => obs_hofx_area_time_avg2d
# ifdef ADJOINT
        PROCEDURE :: area_time_avg2d_ad  => obs_hofx_area_time_avg2d_ad
# endif
!
        PROCEDURE :: time_avg2d          => obs_hofx_time_avg2d
# ifdef ADJOINT
        PROCEDURE :: time_avg2d_ad       => obs_hofx_time_avg2d_ad
# endif
# ifdef SOLVE3D
        PROCEDURE :: time_avg3d          => obs_hofx_time_avg3d
#  ifdef ADJOINT
        PROCEDURE :: time_avg3d_ad       => obs_hofx_time_avg3d_ad
#  endif
# endif
!
!  H(x) Utility: Operations, exchanges, and  writing.
!
        PROCEDURE :: add                 => obs_hofx_add
        PROCEDURE :: rms                 => obs_hofx_rms
        PROCEDURE :: nc_write            => obs_hofx_nc_write
        PROCEDURE :: update              => obs_hofx_update
        PROCEDURE :: zeros               => obs_hofx_zeros

      END TYPE obs_hofx
!
!-----------------------------------------------------------------------
!  Module variables.
!-----------------------------------------------------------------------
!
!  Background error (standard deviation) identifier.
!
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
!  Geometry for area-averaged observations.
!
      TYPE (obs_filter), allocatable :: geom_hofx(:)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      PUBLIC  :: hofx_initialize
      PUBLIC  :: hofx_run
      PUBLIC  :: hofx_finalize
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      CONTAINS
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!  It allocates and initializes the H(x) objects.
!
      SUBROUTINE hofx_initialize (model)
!
      integer,          intent(in   ) :: model      ! kernel identifier
!
      logical                         :: IsAreaAveraged
      integer                         :: i, ng, ntotal, tile
!
!  Get total observations count and area-averaged switch.
!
      ntotal=obs_store%ndatum
      IsAreaAveraged=.FALSE.
      DO i=1,SIZE(obs_store%pool)
        IF (obs_store%pool(i)%IsAreaAveraged) THEN
          IsAreaAveraged=.TRUE.
          EXIT
        END IF
      END DO
!
!  If applicable, allocate and initialize global grid geometry needed
!  for area-averaged Hofx operators.
!
      IF (IsAreaAveraged) THEN
        IF (.not.allocated(geom_hofx)) THEN
          allocate ( geom_hofx(Ngrids) )
        END IF
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL geom_hofx(ng)% create (ng, tile, model)
          END DO
        END DO
      END IF
!
!  Allocate and create nonlinear H(x) object.
!
      IF (model.eq.iNLM) THEN
        IF (.not.allocated(nlm_hofx)) THEN
          allocate ( nlm_hofx(Ngrids) )
        END IF
# ifndef VERIFICATION
        IF (.not.allocated(berr_hofx)) THEN
          allocate ( berr_hofx(Ngrids) )
        END IF
# endif
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL nlm_hofx(ng)%create  (ng, tile, model, ntotal)
# ifndef VERIFICATION
            CALL berr_hofx(ng)%create (ng, tile, model, ntotal)
# endif
          END DO
        END DO
      END IF

# if defined FOUR_DVAR && defined TLM_OBS
!
!  Allocate and create tangent linear H(x) object.
!
      IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
        IF (.not.allocated(nlm_hofx)) THEN
          allocate ( tlm_hofx(Ngrids) )
        END IF
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL tlm_hofx(ng)%create  (ng, tile, model, ntotal)
          END DO
        END DO
      END IF
# endif

# if defined FOUR_DVAR && defined ADJOINT
!
!  Allocate and create adjoint H^(x) object.
!
      IF (model.eq.iADM) THEN
        IF (.not.allocated(nlm_hofx)) THEN
          allocate ( adm_hofx(Ngrids) )
        END IF
!
        DO ng=1,Ngrids
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL adm_hofx(ng)%create  (ng, tile, model, ntotal)
          END DO
        END DO
      END IF
# endif
!
      RETURN
      END SUBROUTINE hofx_initialize
!
!-----------------------------------------------------------------------
!  It computes the H(x) for the specified model kernel.
!
      SUBROUTINE hofx_run (ng, tile, model)
!
      integer,            intent(in   ) :: ng     ! nested grid
      integer,            intent(in   ) :: tile   ! domain partition
      integer,            intent(in   ) :: model  ! kernel identifier
!
      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_initialize"
!
      SourceFile=MyFile
!
!  Compute model at the observation locations.
!
      SELECT CASE (model)
        CASE (iNLM)                                  ! nonlinear
          CALL obs_hofx_extract (nlm_hofx(ng),                          &
     &                           obs_store, geom_hofx(ng), ng, model)

# ifndef VERIFICATION
        CASE (iSTD)                                  ! background error
          CALL obs_hofx_extract (berr_hofx(ng),                         &
     &                           obs_store, geom_hofx(ng), ng, model)
# endif
# if defined FOUR_DVAR && defined TLM_OBS
        CASE (iTLM, iRPM)                            ! tangent linear
          CALL obs_hofx_extract (tlm_hofx(ng),                          &
     &                           obs_store, geom_hofx(ng), ng, model)
# endif
# if defined FOUR_DVAR && defined ADJOINT
        CASE (iADM)                                  ! adjoint
          CALL obs_hofx_extract (adm_hofx(ng),                          &
     &                           obs_store, geom_hofx(ng), ng, model)
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) model
          exit_flag=5
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!

  10  FORMAT (/,' ROMS_HOFX::HOFX_RUN: illegal identifier, model = ',i0)
!
      RETURN
      END SUBROUTINE hofx_run
!
!-----------------------------------------------------------------------
!  It finalizes the H(x) computation.
!
      SUBROUTINE hofx_finalize (model)
!
      integer,           intent(in   ) :: model  ! kernel identifier
!
      integer                          :: ng

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_finalize"
!
      SourceFile=MyFile
!
!  If appropiate, collect H(x) data from all tasks.  Then, write out to
!  output NetCDF files.
!
      SELECT CASE (model)
        CASE (iNLM)                                  ! nonlinear
          DO ng=1,Ngrids
# ifdef DISTRIBUTE
            CALL obs_hofx_update (nlm_hofx(ng), ng, model)
# endif
            CALL obs_hofx_nc_write (nlm_hofx(ng), obs_store, ng, model)
          END DO

# ifndef VERIFICATION
        CASE (iSTD)                                  ! background error
          DO ng=1,Ngrids
#  ifdef DISTRIBUTE
            CALL obs_hofx_update (berr_hofx(ng), ng, model)
#  endif
            CALL obs_hofx_nc_write (berr_hofx(ng), obs_store, ng, model)
          END DO
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) model
          exit_flag=5
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!

  10  FORMAT (/,' ROMS_HOFX::HOFX_FINALIZE: illegal identifier, '       &
     &        'model = ',i0)
!
      RETURN
      END SUBROUTINE hofx_finalize
!
! <><><><><><><><><><><><><><><><><><><><><><><> CLASS OBS_FILTER <><><>
!
!  It allocates observation filter object.
!
      SUBROUTINE obs_filter_create (self, ng, tile, model)
!
      CLASS (obs_filter), intent(inout) :: self   ! H(x) object
      integer,            intent(in   ) :: ng     ! nested grid
      integer,            intent(in   ) :: tile   ! domain partition
      integer,            intent(in   ) :: model  ! kernel identifier
!
      real (r8)                         :: size2d
!
!  Set tiled and global array bounds.
!
      self%tile=tile
!
      self%LBi=BOUNDS(ng)%LBi(tile)
      self%UBi=BOUNDS(ng)%UBi(tile)
      self%LBj=BOUNDS(ng)%LBj(tile)
      self%UBj=BOUNDS(ng)%UBj(tile)
!
      self%Imin=0
      self%Imax=Lm(ng)+1
      self%Jmin=0
      self%Jmax=Mm(ng)+1
      size2d=REAL((Lm(ng)+2)*(Mm(ng)+2),r8)
!
!  Set global area, grid spacing mask required in the spatial averaging
!  filter.
!
# ifdef DISTRIBUTE
      IF (.not.associated(self%area)) THEN
        allocate ( self%area(self%Imin:self%Imax, self%Jmin:self%Jmax) )
        self%area = 0.0_r8
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (.not.associated(self%dx)) THEN
        allocate ( self%dx(self%Imin:self%Imax, self%Jmin:self%Jmax) )
        self%dx = 0.0_r8
        Dmem(ng)=Dmem(ng)+size2d
      END IF
      IF (.not.associated(self%dy)) THEN
        allocate ( self%dy(self%Imin:self%Imax, self%Jmin:self%Jmax) )
        self%dy = 0.0_r8
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  ifdef MASKING
      IF (.not.associated(self%mask)) THEN
        allocate ( self%mask(self%Imin:self%Imax, self%Jmin:self%Jmax) )
        self%mask = 0.0_r8
        Dmem(ng)=Dmem(ng)+size2d
      END IF
#  endif
!
      CALL mp_aggregate2d (ng, model, r2dvar,                           &
     &                     self%LBi, self%UBi, self%LBj, self%UBj,      &
     &                     self%Imin, self%Imax, self%Jmin, self%Jmax,  &
     &                     GRID(ng)%omn, self%area)
      CALL mp_aggregate2d (ng, model, r2dvar,                           &
     &                     self%LBi, self%UBi, self%LBj, self%UBj,      &
     &                     self%Imin, self%Imax, self%Jmin, self%Jmax,  &
     &                     GRID(ng)%om_r, self%dx)
      CALL mp_aggregate2d (ng, model, r2dvar,                           &
     &                     self%LBi, self%UBi, self%LBj, self%UBj,      &
     &                     self%Imin, self%Imax, self%Jmin, self%Jmax,  &
     &                     GRID(ng)%on_r, self%dy)
#  ifdef MASKING
      CALL mp_aggregate2d (ng, model, r2dvar,                           &
     &                     self%LBi, self%UBi, self%LBj, self%UBj,      &
     &                     self%Imin, self%Imax, self%Jmin, self%Jmax,  &
     &                     GRID(ng)%rmask, self%mask)
#  endif
# else
      self%area => GRID(ng)%omn
      self%dx   => GRID(ng)%om_r
      self%dy   => GRID(ng)%on_r
#  ifdef MASKING
      self%mask => GRID(ng)%rmask
#  endif
# endif
!
!  Allocate generic global variables.
!
      IF (.not.allocated(self%global2d)) THEN
        allocate ( self%global2d(self%Imin:self%Imax,                   &
     &                           self%Jmin:self%Jmax) )
        self%global2d = 0.0_r8
        Dmem(ng)=Dmem(ng)+size2d
      END IF
!
      RETURN
      END SUBROUTINE obs_filter_create
!
!-----------------------------------------------------------------------
!  It destroys observation filter object.
!
      SUBROUTINE obs_filter_destroy (self)
!
      CLASS (obs_filter), intent(inout) :: self   ! H(x) object
!
!  It nullifies observation filter object arrays.
!
      IF (associated(self%area)) nullify (self%area)
      IF (associated(self%dx))   nullify (self%dx)
      IF (associated(self%dy))   nullify (self%dy)
#  ifdef MASKING
      IF (associated(self%mask)) nullify (self%mask)
#  endif
!
      IF (allocated(self%global2d)) deallocate (self%global2d)
!
      RETURN
      END SUBROUTINE obs_filter_destroy
!
! <><><><><><><><><><><><><><><><><><><><><><><><> CLASS OBS_HOFX <><><>
!
!  It allocates H(x) object.
!
      SUBROUTINE obs_hofx_create (self, ng, tile, model, ntotal)
!
      CLASS (obs_hofx), intent(inout) :: self     ! H(x) object
      integer,          intent(in   ) :: ng       ! nested grid
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
        Dmem(ng)=Dmem(ng)+REAL(ntotal,r8)
      END IF
!
!  Allocate processing flag to reject (zero) or accept (unity)
!  observations.
!
      IF (.not.allocated(self%ObsVetting)) THEN
        allocate ( self%ObsVetting(ntotal) )
        Dmem(ng)=Dmem(ng)+REAL(ntotal,r8)
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
! <><><><><><><><><><><><><><><><><><><><><><><><> HOFX OPERATORS <><><>
!
!  It computes the H(x) operator, which interpolates model variables to
!  observation locations in space and time. If adjoint model kernel, it
!  computes its transpose, H^(x). Area-averaged and/or time-averaged
!  operators are available.
!
      SUBROUTINE obs_hofx_extract (self, obs, geom, ng, model)
!
      CLASS (obs_hofx),  intent(inout) :: self     ! H(x) object
      TYPE (obs_pool),   intent(inout) :: obs      ! observation object
      TYPE (obs_filter), intent(inout) :: geom     ! geometry object
      integer,           intent(in   ) :: ng       ! nested grid
      integer,           intent(in   ) :: model    ! kernel identifier
!
      TYPE (obs_data),         pointer :: depot
      integer (i8b),       allocatable :: reject(:)
      integer                          :: i, iobs, mvars, nvars
      integer                          :: is, ie, ObsSum, ObsReject

      character (len=3),   allocatable :: op_handle(:)
      character (len=40)               :: varname
      character (len=*),     parameter :: MyFile =                       &
     &    __FILE__//", roms_hofx_extract"
!
      SourceFile=MyFile
!
!  If observations are available for the (time-0.5*dt) to (time+0.5*dt)
!  interval or time-averaged window, it extracts model control variables
!  to the observation locations.
!
      IF (obs%has(ng, model)) THEN
        mvars=MAXVAL(obs%pool(:)%nvars)
        IF (.not.allocated(reject)) THEN
          allocate ( reject(mvars) )
        END IF
        IF (.not.allocated(op_handle)) THEN
          allocate ( op_handle(mvars) )
          DO i=1,mvars
            op_handle(i)='SUM'
          END DO
        END IF
        ObsSum=0
        ObsReject=0
!
        IF (Master) WRITE(stdout,10) tdays(ng), time_code(ng)
!
        DO iobs=1,SIZE(obs%pool)
          IF ((obs%pool(iobs)%NstrObs.gt.0).and.                        &
     &        (obs%pool(iobs)%NendObs.gt.0)) THEN
            nvars=obs%pool(iobs)%nvars
            reject=0
!
            SELECT CASE (obs%pool(iobs)%name)
              CASE ('SSH', 'SSS', 'SST')
                CALL obs_pool_get (obs, obs%pool(iobs)%name, depot)
                CALL self%extract2d (depot, geom, ng, model, reject,    &
     &                               obs%pool(iobs)%name)
# ifdef SOLVE3D
              CASE ('salt', 'temp', 'uv')
                CALL obs_pool_get (obs, obs%pool(iobs)%name, depot)
                CALL self%extract3d (depot, ng, model, reject,          &
     &                               obs%pool(iobs)%name)
# endif
              CASE DEFAULT
                IF (Master) WRITE (stdout,20) obs%pool(iobs)%name, iobs
                exit_flag=5
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
            END SELECT
!
!  Store processing tally.
!
# ifdef DISTRIBUTE
            CALL mp_reduce (ng, model, mvars, reject, op_handle)
# endif
            is=obs%pool(iobs)%NstrObs
            ie=obs%pool(iobs)%NendObs
            DO i=1,nvars
              IF (obs%pool(iobs)%IsTimeAveraged.and.                    &
     &            (obs%pool(iobs)%TimeStepStr.eq.iic(ng))) THEN
                self%ObsCount=self%ObsCount+(ie-is+1)
                self%ObsReject=self%ObsReject+reject(i)
              ELSE IF (.not.obs%pool(iobs)%IsTimeAveraged) THEN
                self%ObsCount=self%ObsCount+(ie-is+1)
                self%ObsReject=self%ObsReject+reject(i)
              END IF
            END DO
!
!  Report.
!
            IF (Master) THEN
              DO i=1,nvars
                IF (obs%pool(iobs)%IsAreaAveraged.and.                  &
     &              obs%pool(iobs)%IsTimeAveraged) THEN
                  varname=TRIM(obs%pool(iobs)%vname(i)) //              &
     &                    ' (AreaTimeAvg)'
                ELSE IF (obs%pool(iobs)%IsAreaAveraged) THEN
                  varname=TRIM(obs%pool(iobs)%vname(i)) //              &
     &                    ' (AreaAvg)'
                ELSE IF (obs%pool(iobs)%IsTimeAveraged) THEN
                  varname=TRIM(obs%pool(iobs)%vname(i)) //              &
     &                    ' (TimeAvg)'
                ELSE
                  varname=TRIM(obs%pool(iobs)%vname(i))
                END IF
                WRITE (stdout,30) TRIM(varname),                        &
     &                            obs%pool(iobs)%TimeStepStr,           &
     &                            obs%pool(iobs)%TimeStepEnd,           &
     &                            is, ie, ie-is+1, reject(i)
                ObsSum=ObsSum+ie-is+1
                ObsReject=ObsReject+reject(i)
              END DO
            END IF
          END IF
        END DO
        IF (Master) WRITE (stdout,40) ObsSum, ObsReject,                &
     &                                self%ObsCount, self%ObsReject
!
!  Deallocate local arrays.
!
        IF (allocated(reject))    deallocate (reject)
        IF (allocated(op_handle)) deallocate (op_handle)
      END IF
!
  10  FORMAT (/,' Number of State Observations Processed:',2x,          &
     &        'ObsTime = ',f12.4,',',t68,a,/,/,                         &
     &        3x,'Variable',14x,'IstrStep',3x,'IendStep',4x,'IstrObs',  &
     &        4x,'IendObs',6x,'Count',3x,'Rejected',/)
  20  FORMAT (/,' ROMS_HOFX::extract: Cannot find name = ',a,2x,        &
     &        'in obs%pool(',i0,').')
  30  FORMAT (3x,a,t23,6(1x,i10))
  40  FORMAT (/,3x,'Total',t67,2(1x,i10),                               &
     &        /,3x,'Obs Tally',t67,2(1x,i10),/)
!
      RETURN
      END SUBROUTINE obs_hofx_extract
!
!-----------------------------------------------------------------------
!  It computes 2D control variables H(x) operator. If adjoint kernel,
!  model=iADM, it computes the transport operator, H^(x). Area-averaged
!  and/or time-averaged operators are available.
!
      SUBROUTINE obs_hofx_extract2d (self, obs, geom, ng, model,        &
     &                               reject, name)
!
      CLASS (obs_hofx),  intent(inout) :: self      ! H(x) object
      TYPE (obs_data),   intent(inout) :: obs       ! observation depot
      TYPE (obs_filter), intent(inout) :: geom      ! geometry object
      integer,           intent(in   ) :: ng        ! nested grid
      integer,           intent(in   ) :: model     ! kernel identifier
      integer (i8b),     intent(out  ) :: reject(:) ! obs reject count
      character (len=*), intent(in   ) :: name      ! obs pool name
!
      logical                          :: IsAreaTimeAveraged
      logical                          :: IsAreaAveraged
      logical                          :: IsTimeAveraged
      integer                          :: LBi, UBi, LBj, UBj
      integer                          :: iTstr, iTend, ivar
# ifdef SOLVE3D
      integer                          :: itrc
# endif
!
      real (r8),               pointer :: field(:,:) => NULL()

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", roms_hofx_extract2d"
!
      SourceFile=MyFile
!
!  Set pointer to kernel variable.
!
      VAR_LOOP : DO ivar=1,obs%nvars
        SELECT CASE (TRIM(obs%vname(ivar)))
          CASE ('ssh', 'SSH')
            IF (model.eq.iNLM) THEN
              field => OCEAN(ng)%zeta(:,:,KOUT)
# ifdef ADJOINT
            ELSE IF (model.eq.iADM) THEN
              field => OCEAN(ng)%f_zeta
# endif
# ifdef TLM_OBS
            ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
              field => OCEAN(ng)%tl_zeta(:,:,KOUT)
# endif
# ifndef VERIFICATION
            ELSE IF (model.eq.iSTD) THEN
              field => OCEAN(ng)%e_zeta(:,:,1)          ! NSA=1 array
# endif
            END IF
# ifdef SOLVE3D
          CASE ('sss', 'sst', 'SSS', 'SST')
            itrc=tracer_index(obs%name)
            IF (model.eq.iNLM) THEN
              field => OCEAN(ng)%t(:,:,N(ng),itrc,NOUT)
#  ifdef ADJOINT
            ELSE IF (model.eq.iADM) THEN
              field => OCEAN(ng)%f_t(:,:,N(ng),itrc)
#  endif
#  ifdef TLM_OBS
            ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
              field => OCEAN(ng)%tl_t(:,:,N(ng),itrc, NOUT)
#  endif
#  ifndef VERIFICATION
            ELSE IF (model.eq.iSTD) THEN
              field => OCEAN(ng)%e_t(:,:,N(ng),itrc,1) ! NSA=1 array
#  endif
            END IF
          CASE ('u', 'Ua', 'Ucodar')
            IF (model.eq.iNLM) THEN
              field => OCEAN(ng)%ua(:,:,N(ng))
#  ifdef ADJOINT
            ELSE IF (model.eq.iADM) THEN
              field => OCEAN(ng)%f_u(:,:,N(ng))
#  endif
#  ifdef TLM_OBS
            ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
              field => OCEAN(ng)%tl_ua(:,:,N(ng))
#  endif
#  ifndef VERIFICATION
            ELSE IF (model.eq.iSTD) THEN
              field => OCEAN(ng)%e_u(:,:,N(ng), 1)      ! NSA=1 array
#  endif
            END IF
          CASE ('v', 'Va', 'Vcodar')
            IF (model.eq.iNLM) THEN
              field => OCEAN(ng)%va(:,:,N(ng))
#  ifdef ADJOINT
            ELSE IF (model.eq.iADM) THEN
              field => OCEAN(ng)%f_v(:,:,N(ng))
#  endif
#  ifdef TLM_OBS
            ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
              field => OCEAN(ng)%tl_va(:,:,N(ng))
#  endif
#  ifndef VERIFICATION
            ELSE IF (model.eq.iSTD) THEN
              field => OCEAN(ng)%e_v(:,:,N(ng), 1)      ! NSA=1 array
#  endif
            END IF
# endif
          CASE DEFAULT
            IF (Master) WRITE (stdout,10) TRIM(obs%name)
            exit_flag=5
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END SELECT
!
!  Set switch for filters.
!
        IsAreaTimeAveraged = obs%IsAreaAveraged.and.                    &
     &                       obs%IsTimeAveraged
        IsAreaAveraged     = .not.obs%IsTimeAveraged.and.               &
     &                            obs%IsAreaAveraged
        IsTimeAveraged     = .not.obs%IsAreaAveraged.and.               &
     &                            obs%IsTimeAveraged
        IF (obs%IsTimeAveraged) THEN
          CALL obs_data_averaging_timesteps (obs, ng, model,            &
     &                                       iTstr, iTend)
          obs%TimeStepStr=iTstr
          obs%TimeStepEnd=iTend
        ELSE
          obs%TimeStepStr=iic(ng)
          obs%TimeStepEnd=iic(ng)
        END IF
!
!  Compute 2D control variable H(x) operator.
!
        LBi=self%LBi
        UBi=self%UBi
        LBj=self%LBj
        UBj=self%UBj
!
        IF (model.ne.iADM) THEN
          IF (IsAreaTimeAveraged) THEN
            CALL self%area_time_avg2d (obs, geom, ng, model,            &
     &                                 ivar, reject(ivar),              &
     &                                 iTstr, iTend,                    &
     &                                 LBi, UBi, LBj, UBj, field)
          ELSE IF (IsAreaAveraged) THEN
            CALL self%area_avg2d (obs, geom, ng, model,                 &
     &                            ivar, reject(ivar),                   &
     &                            LBi, UBi, LBj, UBj, field)
          ELSE IF (IsTimeAveraged) THEN
            CALL self%time_avg2d (obs, ng, model,                       &
     &                            ivar, reject(ivar),                   &
     &                            iTstr, iTend,                         &
     &                            LBi, UBi, LBj, UBj, field)
          ELSE
            CALL self%interp2d (obs, ng, model,                         &
     &                          ivar, reject(ivar),                     &
     &                          LBi, UBi, LBj, UBj, field)
          END IF
# ifdef ADJOINT
        ELSE IF (model.eq.iADM) THEN
          IF (IsAreaTimeAveraged) THEN
            CALL self%area_time_avg2d_ad (obs, geom, ng, model,         &
     &                                    ivar, reject(ivar),           &
     &                                    iTstr, iTend,                 &
     &                                    LBi, UBi, LBj, UBj, field)
          ELSE IF (IsAreaAveraged) THEN
            CALL self%area_avg2d_ad (obs, geom, ng, model,              &
     &                               ivar, reject(ivar),                &
     &                               LBi, UBi, LBj, UBj, field)
          ELSE IF (IsTimeAveraged) THEN
            CALL self%time_avg2d_ad (obs, ng, model,                    &
     &                               ivar, reject(ivar),                &
     &                               iTstr, iTend,                      &
     &                               LBi, UBi, LBj, UBj, field)
          ELSE
            CALL self%interp2d_ad (obs, ng, model,                      &
     &                             ivar, reject(ivar),                  &
     &                             LBi, UBi, LBj, UBj, field)
          END IF
# endif
        END IF
      END DO VAR_LOOP
!
!  Nullify pointer.
!
      IF (associated(field)) nullify (field)
!
  10  FORMAT (/,' ROMS_HOFX::extract2d: Cannot find name = ',a,2x,      &
     &        'in observation depot.')
!
      RETURN
      END SUBROUTINE obs_hofx_extract2d

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  It computes 3D control variables H(x) operator. If adjoint kernel,
!  model=iADM, it computes the transport operator, H^(x). Area-averaged
!  and/or time-averaged operators are available.
!
      SUBROUTINE obs_hofx_extract3d (self, obs, ng, model,              &
     &                               reject, name)
!
      CLASS (obs_hofx),  intent(inout) :: self      ! H(x) object
      TYPE (obs_data),   intent(inout) :: obs       ! observation depot
      integer,           intent(in   ) :: ng        ! nested grid
      integer,           intent(in   ) :: model     ! kernel identifier
      integer (i8b),     intent(out  ) :: reject(:) ! obs reject count
      character (len=*), intent(in   ) :: name      ! obs pool name
!
      logical                          :: IsTimeAveraged
      integer                          :: LBi, UBi, LBj, UBj
      integer                          :: iTstr, iTend
      integer                          :: itrc, ivar
!
      real (r8),               pointer :: field(:,:,:) => NULL()

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", roms_hofx_extract3d"
!
      SourceFile=MyFile
!
!  Set pointer to kernel variable.
!
      VAR_LOOP : DO ivar=1,obs%nvars
        SELECT CASE (TRIM(obs%vname(ivar)))
          CASE ('salt', 'temp', 'tracer')
            itrc=tracer_index(obs%name)
            IF (model.eq.iNLM) THEN
              field => OCEAN(ng)%t(:,:,:,itrc,NOUT)
#  ifdef ADJOINT
            ELSE IF (model.eq.iADM) THEN
              field => OCEAN(ng)%f_t(:,:,:,itrc)
#  endif
#  ifdef TLM_OBS
            ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
              field => OCEAN(ng)%tl_t(:,:,:,itrc, NOUT)
#  endif
#  ifndef VERIFICATION
            ELSE IF (model.eq.iSTD) THEN
              field => OCEAN(ng)%e_t(:,:,:,itrc,1)      ! NSA=1 array
#  endif
            END IF
          CASE ('u', 'Ua', 'Ucodar')
            IF (model.eq.iNLM) THEN
              field => OCEAN(ng)%ua
#  ifdef ADJOINT
            ELSE IF (model.eq.iADM) THEN
              field => OCEAN(ng)%f_u
#  endif
#  ifdef TLM_OBS
            ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
              field => OCEAN(ng)%tl_ua
#  endif
#  ifndef VERIFICATION
            ELSE IF (model.eq.iSTD) THEN
              field => OCEAN(ng)%e_u(:,:,:, 1)          ! NSA=1 array
#  endif
            END IF
          CASE ('v', 'Va', 'Vcodar')
            IF (model.eq.iNLM) THEN
              field => OCEAN(ng)%va
#  ifdef ADJOINT
            ELSE IF (model.eq.iADM) THEN
              field => OCEAN(ng)%f_v
#  endif
#  ifdef TLM_OBS
            ELSE IF ((model.eq.iTLM).or.(model.eq.iRPM)) THEN
              field => OCEAN(ng)%tl_va
#  endif
#  ifndef VERIFICATION
            ELSE IF (model.eq.iSTD) THEN
              field => OCEAN(ng)%e_v(:,:,:, 1)          ! NSA=1 array
#  endif
            END IF
          CASE DEFAULT
            IF (Master) WRITE (stdout,10) TRIM(obs%vname(ivar))
            exit_flag=5
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END SELECT
!
!  Set switch for filters.
!
        IsTimeAveraged = .not.obs%IsAreaAveraged.and.                   &
     &                        obs%IsTimeAveraged
        IF (obs%IsTimeAveraged) THEN
          CALL obs_data_averaging_timesteps (obs, ng, model,            &
     &                                       iTstr, iTend)
          obs%TimeStepStr=iTstr
          obs%TimeStepEnd=iTend
        ELSE
          obs%TimeStepStr=iic(ng)
          obs%TimeStepEnd=iic(ng)
        END IF
!
!  Compute 3D control variable H(x) operator.
!
        LBi=self%LBi
        UBi=self%UBi
        LBj=self%LBj
        UBj=self%UBj
!
        IF (model.ne.iADM) THEN
          IF (IsTimeAveraged) THEN
            CALL self%time_avg3d (obs, ng, model,                       &
     &                            ivar, reject(ivar),                   &
     &                            iTstr, iTend,                         &
     &                            LBi, UBi, LBj, UBj, field)
          ELSE
            CALL self%interp3d (obs, ng, model,                         &
     &                          ivar, reject(ivar),                     &
     &                          LBi, UBi, LBj, UBj, field)
          END IF
#  ifdef ADJOINT
        ELSE IF (model.eq.iADM) THEN
          IF (IsTimeAveraged) THEN
            CALL self%time_avg3d_ad (obs, ng, model,                    &
     &                               ivar, reject(ivar),                &
     &                               iTstr, iTend,                      &
     &                               LBi, UBi, LBj, UBj, field)
          ELSE
            CALL self%interp3d_ad (obs, ng, model,                      &
     &                             ivar, reject(ivar),                  &
     &                             LBi, UBi, LBj, UBj, field)
          END IF
#  endif
        END IF
      END DO VAR_LOOP
!
!  Nullify pointer.
!
      IF (associated(field)) nullify (field)
!
  10  FORMAT (/,' ROMS_HOFX::extract3d: Cannot find name = ',a,2x,      &
     &        'in observation depot.')
!
      RETURN
      END SUBROUTINE obs_hofx_extract3d
# endif
!
!-----------------------------------------------------------------------
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
      SUBROUTINE obs_hofx_interp2d (self, obs, ng, model,               &
     &                              ivar, reject,                       &
     &                              LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),  intent(inout) :: self       ! H(x) object
      TYPE  (obs_data),  intent(in   ) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(out  ) :: reject     ! obs reject count
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(in   ) :: field(LBi:,LBj:)
!
      integer                          :: i, ic, iobs, i1, i2, j1, j2
      real (r8)                        :: p1, p2, q1, q2, wsum
      real (r8),          dimension(4) :: Hmat
!
!  Interpolate 2D model state variable at the observation locations.
!
      reject=0
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
          IF (wsum.gt.0.0_r8) THEN
            self%obsVetting(iobs)=1
          ELSE
            self%obsVetting(iobs)=0
            reject=reject+1
          END IF
# else
          self%ObsVetting(iobs)=1
# endif
        END IF
      END DO
!
      RETURN
      END SUBROUTINE obs_hofx_interp2d

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
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
      SUBROUTINE obs_hofx_interp3d (self, obs, ng, model,               &
     &                              ivar, reject,                       &
     &                              LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),  intent(inout) :: self       ! H(x) object
      TYPE  (obs_data),  intent(inout) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(out  ) :: reject     ! obs reject count
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(in   ) :: field(LBi:,LBj:,:)
!
      integer                          :: i, ic, iobs, k
      integer                          :: i1, i2, j1, j2, k1, k2
      real (r8)                        :: p1, p2, q1, q2, r1, r2
      real (r8)                        :: w11, w12, w21, w22, wsum
      real (r8)                        :: Zbot, Ztop, dz
      real (r8),          dimension(8) :: Hmat
!
!  Interpolate 3D model state variable at the observation locations.
!
      reject=0
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
              reject=reject+1
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
            self%state(iobs)=Hmat(1)*field(i1,j1,k1)+                   &
     &                       Hmat(2)*field(i2,j1,k1)+                   &
     &                       Hmat(3)*field(i2,j2,k1)+                   &
     &                       Hmat(4)*field(i1,j2,k1)+                   &
     &                       Hmat(5)*field(i1,j1,k2)+                   &
     &                       Hmat(6)*field(i2,j1,k2)+                   &
     &                       Hmat(7)*field(i2,j2,k2)+                   &
     &                       Hmat(8)*field(i1,j2,k2)
#  ifdef MASKING
            IF (wsum.gt.0.0_r8) THEN
              self%ObsVetting(iobs)=1
            ELSE
              self%ObsVetting(iobs)=0
              reject=reject+1
            END IF
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
              reject=reject+1
            END IF
#  endif
          END IF
        END IF
      END DO
!
      RETURN
      END SUBROUTINE obs_hofx_interp3d
# endif

# ifdef ADJOINT
!
!-----------------------------------------------------------------------
!  Adjoint of interpolating 2D state variables at the observation
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
      SUBROUTINE obs_hofx_interp2d_ad (ad_self, obs, ng, model,         &
     &                                 ivar, reject,                    &
     &                                 LBi, UBi, LBj, UBj, ad_field)
!
      CLASS (obs_hofx),  intent(inout) :: ad_self    ! H^(x) object
      TYPE  (obs_data),  intent(in   ) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(out  ) :: reject     ! obs reject count
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(inout) :: ad_field(LBi:,LBj:)
!
      integer                          :: i, ic, iobs, i1, i2, j1, j2
      real (r8)                        :: p1, p2, q1, q2, wsum
      real (r8),          dimension(4) :: Hmat
!
!  Interpolate 2D model state variable at the observation locations.
!
      reject=0
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
          Hmat(1)=Hmat(1)*ad_self%mask(i1,j1)
          Hmat(2)=Hmat(2)*ad_self%mask(i2,j1)
          Hmat(3)=Hmat(3)*ad_self%mask(i2,j2)
          Hmat(4)=Hmat(4)*ad_self%mask(i1,j2)
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
     &                    Hmat(1)*ad_self%state(iobs)
          ad_field(i2,j1)=ad_field(i2,j1)+                              &
     &                    Hmat(2)*ad_self%state(iobs)
          ad_field(i2,j2)=ad_field(i2,j2)+                              &
     &                    Hmat(3)*ad_self%state(iobs)
          ad_field(i1,j2)=ad_field(i1,j2)+                              &
     &                    Hmat(4)*ad_self%state(iobs)
          ad_self%state(iobs)=0.0_r8
#  ifdef MASKING
          IF (wsum.gt.0.0_r8) THEN
            ad_self%ObsVetting(iobs)=1
          ELSE
            ad_self%ObsVetting(iobs)=0
            reject=reject+1
          END IF
#  else
          ad_self%ObsVetting(iobs)=1
#  endif
        END IF
      END DO

#  ifdef DISTRIBUTE
!
!  Exchange tile data.
!
      CALL ad_mp_exchange2d (ng, ad_self%tile, model, 1,                &
     &                       ad_self%LBi, ad_self%UBi,                  &
     &                       ad_self%LBj, ad_self%UBj,                  &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_field)
#  endif
!
  10  FORMAT (/,' ROMS_HOFX::interp2d_ad: Cannot find name = ',a,2x,    &
     &        'in obs%data observation depot object.')
!
      RETURN
      END SUBROUTINE obs_hofx_interp2d_ad

#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
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
!                             1             2
!
      SUBROUTINE obs_hofx_interp3d_ad (ad_self, obs, ng, model,         &
     &                                 ivar, reject,                    &
     &                                 LBi, UBi, LBj, UBj, ad_field)
!
      CLASS (obs_hofx),  intent(inout) :: ad_self    ! H^(x) object
      TYPE  (obs_data),  intent(inout) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(out  ) :: reject     ! obs reject count
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(inout) :: ad_field(LBi:,LBj:,:)
!
      integer                          :: i, ic, iobs, k
      integer                          :: i1, i2, j1, j2, k1, k2
      real (r8)                        :: p1, p2, q1, q2, r1, r2
      real (r8)                        :: w11, w12, w21, w22, wsum
      real (r8)                        :: Zbot, Ztop, dz
      real (r8),          dimension(8) :: Hmat
!
!  Interpolate 2D model state variable at the observation locations.
!
      reject=0
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
            Ztop=ad_self%depth(i1,j1,N(ng))
            Zbot=ad_self%depth(i1,j1,1    )
            IF (obs%Zgrid(i).ge.Ztop) THEN
              k1=N(ng)                         ! If shallower, assign to
              k2=N(ng)                         ! top grid cell. The
              r1=1.0_r8                        ! observation is located
              r2=0.0_r8                        ! on the upper cell half
              obs%Zgrid(i)=REAL(N(ng),r8)      ! above its middle depth.
            ELSE IF (Zbot.ge.obs%Zgrid(iobs)) THEN
              r1=0.0_r8                        ! If deeper, ignore.
              r2=0.0_r8
              ad_self%ObsVetting(iobs)=0
              reject=reject+1
            ELSE
              DO k=N(ng),2,-1                  ! Otherwise, interpolate
                Ztop=ad_self%depth(i1,j1,k  )  ! to fractional level
                Zbot=ad_self%depth(i1,j1,k-1)
                IF ((Ztop.gt.obs%Zgrid(i)).and.                         &
     &              (obs%Zgrid(i).ge.Zbot)) THEN
                  k1=k-1
                  k2=k
                END IF
              END DO
              dz=ad_self%depth(i1,j1,k2)-ad_self%depth(i1,j1,k1)
              r2=(obs%Zgrid(i)-ad_self%depth(i1,j1,k1))/dz
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
            Hmat(1)=Hmat(1)*ad_self%mask(i1,j1)
            Hmat(2)=Hmat(2)*ad_self%mask(i2,j1)
            Hmat(3)=Hmat(3)*ad_self%mask(i2,j2)
            Hmat(4)=Hmat(4)*ad_self%mask(i1,j2)
            Hmat(5)=Hmat(5)*ad_self%mask(i1,j1)
            Hmat(6)=Hmat(6)*ad_self%mask(i2,j1)
            Hmat(7)=Hmat(7)*ad_self%mask(i2,j2)
            Hmat(8)=Hmat(8)*ad_self%mask(i1,j2)
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
!^          self%state(iobs)=Hmat(1)*field(i1,j1,k1)+                   &
!^   &                       Hmat(2)*field(i2,j1,k1)+                   &
!^   &                       Hmat(3)*field(i2,j2,k1)+                   &
!^   &                       Hmat(4)*field(i1,j2,k1)+                   &
!^   &                       Hmat(5)*field(i1,j1,k2)+                   &
!^   &                       Hmat(6)*field(i2,j1,k2)+                   &
!^   &                       Hmat(7)*field(i2,j2,k2)+                   &
!^   &                       Hmat(8)*field(i1,j2,k2)
!^
            ad_field(i1,j1,k1)=ad_field(i1,j1,k1)+                      &
     &                         Hmat(1)*ad_self%state(iobs)
            ad_field(i2,j1,k1)=ad_field(i2,j1,k1)+                      &
     &                         Hmat(2)*ad_self%state(iobs)
            ad_field(i2,j2,k1)=ad_field(i2,j2,k1)+                      &
     &                         Hmat(3)*ad_self%state(iobs)
            ad_field(i1,j2,k1)=ad_field(i1,j2,k1)+                      &
     &                         Hmat(4)*ad_self%state(iobs)
            ad_field(i1,j1,k2)=ad_field(i1,j1,k2)+                      &
     &                         Hmat(5)*ad_self%state(iobs)
            ad_field(i2,j1,k2)=ad_field(i2,j1,k2)+                      &
     &                         Hmat(6)*ad_self%state(iobs)
            ad_field(i2,j2,k2)=ad_field(i2,j2,k2)+                      &
     &                         Hmat(7)*ad_self%state(iobs)
            ad_field(i1,j2,k2)=ad_field(i1,j2,k2)+                      &
     &                         Hmat(8)*ad_self%state(iobs)
            ad_self%state(iobs)=0.0_r8
#   ifdef MASKING
            IF (wsum.gt.0.0_r8) THEN
              ad_self%ObsVetting(iobs)=1
            ELSE
              ad_self%ObsVetting(iobs)=0
              reject=reject+1
            END IF
#   else
            ad_self%ObsVetting(iobs)=1
#   endif
          END IF
        END IF
      END DO

#   ifdef DISTRIBUTE
!
!  Exchange tile data.
!
      CALL ad_mp_exchange3d (ng, ad_self%tile, model, 1,                &
     &                       ad_self%LBi, ad_self%UBi,                  &
     &                       ad_self%LBj, ad_self%UBj,                  &
     &                       1, N(ng), NghostPoints,                    &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_field)
#   endif
!
      RETURN
      END SUBROUTINE obs_hofx_interp3d_ad
#  endif
# endif
!
!<><><><><><><><><><><><><><><><><><><><><><><> HOFX WITH FILTERS <><><>
!
!  It computes 2D state variables at the area-averaged observation
!  locations. The 2D state variable is area-averaged according to the
!  "spatialAverage" half-length scale, LS, specified in the IODA-type
!  observation file. Area-averaged observations are assumed to be at
!  RHO-points (A-grid).
!
!  The state variable area-average computed below is over a square
!  region of 2*LS by 2*LS.
!
!                             4____________3  (i2,j2)
!                             |            |
!                             | * (p,q)    |
!                             |            |
!                             |____________|
!                    (i1,j1)  1            2
!
      SUBROUTINE obs_hofx_area_avg2d (self, obs, geom, ng, model,       &
     &                                ivar, reject,                     &
     &                                LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),   intent(inout) :: self      ! H(x) object
      TYPE  (obs_data),   intent(in   ) :: obs       ! observation depot
      TYPE  (obs_filter), intent(inout) :: geom      ! geometry
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: model     ! kernel identifier
      integer,            intent(in   ) :: ivar      ! variable counter
      integer (i8b),      intent(out  ) :: reject    ! obs reject count
      integer,            intent(in   ) :: LBi, UBi  ! field I-bounds
      integer,            intent(in   ) :: LBj, UBj  ! field J-bounds
      real (r8),          intent(in   ) :: field(LBi:,LBj:)
!
      integer                         :: iobs, i1, j1, m, nscope
      integer                         :: i, is, ie, j, js, je
      real (r8)                       :: dfr, dfx, dfy
      real (r8)                       :: area_sum, resol

# ifdef DISTRIBUTE
!
!  Gather state variable from all the tasks in the group to build a
!  global array.
!
      CALL mp_aggregate2d (ng, model, r2dvar,                           &
     &                     geom%LBi, geom%UBi, geom%LBj, geom%UBj,      &
     &                     geom%Imin, geom%Imax, geom%Jmin, geom%Jmax,  &
     &                     field, geom%global2d)
# endif
!
!  Compute area-averaged 2D model state variable at the area_averaged
!  observation locations.
!
      reject=0
      DO m=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(m)).and.                            &
            (obs%Xgrid(m).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(m)).and.                            &
     &      (obs%Ygrid(m).lt.rYmax(ng))) THEN
!
          iobs=obs%pack_mapping(m,ivar)               ! index mapping
!
          i1=INT(obs%Xgrid(m))
          j1=INT(obs%Ygrid(m))
!
!  Set the area of influence in grid units using the local average grid
!  resolution.
!
          resol=0.5_r8*(geom%dx(i1,j1)+geom%dy(i1,j1))
          nscope=INT(obs%spatialAverage/resol)
          is=MAX(geom%Imin, i1-nscope)
          ie=MIN(i1+nscope, geom%Imax)
          js=MAX(geom%Jmin, j1-nscope)
          je=MIN(j1+nscope, geom%Jmax)
!
!  Average state variable over the area of influence.
!
          area_sum=0.0_r8
          self%state(iobs)=0.0_r8
          DO j=js,je
            dfy=REAL(j-j1, r8)
            DO i=is,ie
# ifdef MASKING
              IF (geom%mask(i,j).lt.1.0_r8) CYCLE
# endif
              dfx=REAL(i-i1, r8)
              dfr=REAL(nscope, r8)
              IF (SQRT(dfx*dfx+dfy*dfy).le.dfr) THEN
                area_sum=area_sum+geom%area(i,j)
                self%state(iobs)=self%state(iobs)+                      &
# ifdef DISTRIBUTE
     &                           geom%area(i,j)*geom%global2d(i,j)
# else
     &                           geom%area(i,j)*field(i,j)
# endif
              END IF
            END DO
          END DO
!
!  Compute area-weighted average.
!
          IF (area_sum.gt.0.0_r8) THEN
            self%state(iobs)=self%state(iobs)/area_sum
            self%ObsVetting(iobs)=1
          ELSE
            self%ObsVetting(iobs)=0
            reject=reject+1
          END IF
        END IF
      END DO
!
  10  FORMAT (/,' ROMS_HOFX::AREA_AVG2d: Cannot find name = ',a,2x,     &
     &        'in obs%data observation depot object.')
!
      RETURN
      END SUBROUTINE obs_hofx_area_avg2d
!
!-----------------------------------------------------------------------
!  It computes 2D state variables at the area-averaged and time-averaged
!  observation locations. The 2D state variable is area-averaged using
!  the "spatialAverage" half-length scale, LS, and then time-averaged
!  over the specified period, as defined in the IODA observation file.
!  Area- and time-averaged observations are assumed to be at RHO points
!  (A-grid).
!
!  The state variable area-average computed below is over a square
!  region of 2*LS by 2*LS.
!
!                             4____________3  (i2,j2)
!                             |            |
!                             | * (p,q)    |
!                             |            |
!                             |____________|
!                    (i1,j1)  1            2
!
      SUBROUTINE obs_hofx_area_time_avg2d (self, obs, geom, ng, model,  &
     &                                     ivar, reject, iTstr, iTend,  &
     &                                     LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),   intent(inout) :: self      ! H(x) object
      TYPE  (obs_data),   intent(in   ) :: obs       ! observation depot
      TYPE  (obs_filter), intent(inout) :: geom      ! geometry object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: model     ! kernel identifier
      integer,            intent(in   ) :: ivar      ! variable counter
      integer (i8b),      intent(out  ) :: reject    ! obs reject count
      integer,            intent(in   ) :: iTstr     ! starting timestep
      integer,            intent(in   ) :: iTend     ! ending   timestep
      integer,            intent(in   ) :: LBi, UBi  ! field I-bounds
      integer,            intent(in   ) :: LBj, UBj  ! field J-bounds
      real (r8),          intent(in   ) :: field(LBi:,LBj:)
!
      integer                           :: iobs, i1, j1, m, nscope
      integer                           :: i, is, ie, j, js, je
      real (r8)                         :: dfr, dfx, dfy
      real (r8)                         :: Tsteps, area_sum, resol

      character (len=*),      parameter :: MyFile =                     &
     &    __FILE__//", roms_hofx_area_time_avg2d"
!
      SourceFile=MyFile

# ifdef DISTRIBUTE
!
!  Gather state variable from all the tasks in the group to build a
!  global array.
!
      CALL mp_aggregate2d (ng, model, r2dvar,                           &
     &                     geom%LBi, geom%UBi, geom%LBj, geom%UBj,      &
     &                     geom%Imin, geom%Imax, geom%Jmin, geom%Jmax,  &
     &                     field, geom%global2d)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
!
!  Compute area-averaged 2D model state variable at the area_averaged
!  observation locations.
!
      reject=0
      DO m=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(m)).and.                            &
            (obs%Xgrid(m).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(m)).and.                            &
     &      (obs%Ygrid(m).lt.rYmax(ng))) THEN
!
          iobs=obs%pack_mapping(m,ivar)               ! index mapping
!
          i1=INT(obs%Xgrid(m))
          j1=INT(obs%Ygrid(m))
!
!  Set the area of influence in grid units using the local average grid
!  resolution.
!
          resol=0.5_r8*(geom%dx(i1,j1)+geom%dy(i1,j1))
          nscope=INT(obs%spatialAverage/resol)
          is=MAX(geom%Imin, i1-nscope)
          ie=MIN(i1+nscope, geom%Imax)
          js=MAX(geom%Jmin, j1-nscope)
          je=MIN(j1+nscope, geom%Jmax)
!
          Tsteps=REAL(iTend-iTstr, r8)
!
!  Average state variable over the area of influence and time window.
!
          area_sum=0.0_r8
          IF ((iTstr.le.iic(ng)).and.(iic(ng).le.iTend)) THEN
            DO j=js,je
              dfy=REAL(j-j1, r8)
              DO i=is,ie
# ifdef MASKING
                IF (geom%mask(i,j).lt.1.0_r8) CYCLE
# endif
                dfx=REAL(i-i1, r8)
                dfr=REAL(nscope, r8)
                IF (SQRT(dfx*dfx+dfy*dfy).le.dfr) THEN
                  area_sum=area_sum+geom%area(i,j)
                  self%state(iobs)=self%state(iobs)+                    &
# ifdef DISTRIBUTE
     &                             geom%area(i,j)*geom%global2d(i,j)
# else
     &                             geom%area(i,j)*field(i,j)
# endif
                END IF
              END DO
            END DO
          END IF
!
!  Compute area-weighted average.
!
          IF ((area_sum.gt.0.0_r8).and.(iic(ng).eq.iTend)) THEN
            self%state(iobs)=self%state(iobs)/(area_sum*Tsteps)
          END IF
        END IF
        self%ObsVetting(iobs)=1
      END DO
!
!  Clear/reset global array.
!
      geom%global2d=0.0_r8
!
  10  FORMAT (/,' ROMS_HOFX::AREA_TIME_AVG2D: Cannot find name = ',a,   &
     &        2x,'in obs%data observation depot object.')
!
      RETURN
      END SUBROUTINE obs_hofx_area_time_avg2d

# ifdef ADJOINT
!
!-----------------------------------------------------------------------
!  Adjoint of computing 2D state variables at the area-averaged
!  observation locations. The 2D state variable is area-averaged
!  according to the "spatialAverage" half-length scale, LS, specified
!  in the IODA-type observation file. Area-averaged observations are
!  assumed to be at RHO-points (A-grid).
!
!  The state variable area-average computed below is over a square
!  region of 2*LS by 2*LS.
!
!                             4____________3  (i2,j2)
!                             |            |
!                             | * (p,q)    |
!                             |            |
!                             |____________|
!                    (i1,j1)  1            2
!
      SUBROUTINE obs_hofx_area_avg2d_ad (ad_self, obs, geom, ng,        &
     &                                   model, ivar, reject,           &
     &                                   LBi, UBi, LBj, UBj, ad_field)
!
      CLASS (obs_hofx),   intent(inout) :: ad_self   ! H'(x) object
      TYPE  (obs_data),   intent(in   ) :: obs       ! observation depot
      TYPE  (obs_filter), intent(inout) :: geom      ! geometry object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: model     ! kernel identifier
      integer,            intent(in   ) :: ivar      ! variable counter
      integer (i8b),      intent(out  ) :: reject    ! obs reject count
      integer,            intent(in   ) :: LBi, UBi  ! field I-bounds
      integer,            intent(in   ) :: LBj, UBj  ! field J-bounds
      real (r8),          intent(inout) :: ad_field(LBi:,LBj:)
!
      integer                         :: iobs, i1, j1, m, nscope
      integer                         :: i, is, ie, j, js, je
      real (r8)                       :: dfr, dfx, dfy
      real (r8)                       :: area_sum, resol
!
!  Clear work global 2D array to hold the adjoint of the area of
!  influence data.
!
      geom%global2d=0.0_r8
!
!  Compute area-averaged 2D model state variable at the area_averaged
!  observation locations.
!
      reject=0
      DO m=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(m)).and.                            &
            (obs%Xgrid(m).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(m)).and.                            &
     &      (obs%Ygrid(m).lt.rYmax(ng))) THEN
!
          iobs=obs%pack_mapping(m,ivar)               ! index mapping
!
          i1=INT(obs%Xgrid(m))
          j1=INT(obs%Ygrid(m))
!
!  Set the area of influence in grid units using the local average grid
!  resolution.
!
          resol=0.5_r8*(geom%dx(i1,j1)+geom%dy(i1,j1))
          nscope=INT(obs%spatialAverage/resol)
          is=MAX(geom%Imin, i1-nscope)
          ie=MIN(i1+nscope, geom%Imax)
          js=MAX(geom%Jmin, j1-nscope)
          je=MIN(j1+nscope, geom%Jmax)
!
!  Adjoint of average state variable over the area of influence.
!
          area_sum=0.0_r8               ! compute area weights
          DO j=js,je
            dfy=REAL(j-j1, r8)
            DO i=is,ie
#  ifdef MASKING
              IF (geom%mask(i,j).lt.1.0_r8) CYCLE
#  endif
              dfx=REAL(i-i1, r8)
              dfr=REAL(nscope, r8)
              IF (SQRT(dfx*dfx+dfy*dfy).le.dfr) THEN
                area_sum=area_sum+geom%area(i,j)
              END IF
            END DO
          END DO
!
          IF (area_sum.gt.0.0_r8) THEN
!>          self%state(iobs)=ad_self%state(iobs)/area_sum
!>          self%ObsVetting(iobs)=1
!>
            ad_self%state(iobs)=ad_self%state(iobs)/area_sum
            ad_self%ObsVetting(iobs)=1
          ELSE
!>          ad_self%ObsVetting(iobs)=0
!>
            ad_self%ObsVetting(iobs)=0
            reject=reject+1
          END IF
!
          DO j=js,je
            dfy=REAL(j-j1, r8)
            DO i=is,ie
#  ifdef MASKING
              IF (geom%mask(i,j).lt.1.0_r8) CYCLE
#  endif
              dfx=REAL(i-i1, r8)
              dfr=REAL(nscope, r8)
              IF (SQRT(dfx*dfx+dfy*dfy).le.dfr) THEN
!>              self%state(iobs)=self%state(iobs)+                      &
!>   &                           geom%area(i,j)*geom%global2d(i,j)
!>
                geom%global2d(i,j)=geom%global2d(i,j)+                  &
     &                             geom%area(i,j)*ad_self%state(iobs)
              END IF
            END DO
          END DO
          ad_self%state(iobs)=0.0_r8
        END IF
      END DO
!
!  Add contribution to 2D state adjoint variable.
!
      DO j=ad_self%JstrR,ad_self%JendR
        DO i=ad_self%JstrR,ad_self%JendR
          ad_field(i,j)=geom%global2d(i,j)
          geom%global2d(i,j)=0.0_r8
        END DO
      END DO
!
      RETURN
      END SUBROUTINE obs_hofx_area_avg2d_ad
!
!-----------------------------------------------------------------------
!  Adjoint of computing 2D state variables at the area-averaged and
!  time-averaged observation locations. The 2D state variable is
!  area-averaged using the "spatialAverage" half-length scale, LS, and
!  then time-averaged over the specified period, as defined in the IODA
!  observation file. Area- and time-averaged observations are assumed
!  to be at RHO points (A-grid).
!
!  The state variable area-average computed below is over a square
!  region of 2*LS by 2*LS.
!
!                             4____________3  (i2,j2)
!                             |            |
!                             | * (p,q)    |
!                             |            |
!                             |____________|
!                    (i1,j1)  1            2
!
      SUBROUTINE obs_hofx_area_time_avg2d_ad (ad_self, obs, geom,       &
     &                                        ng, model, ivar, reject,  &
     &                                        iTstr, iTend,             &
     &                                        LBi, UBi, LBj, UBj,       &
     &                                        ad_field)
!
      CLASS (obs_hofx),   intent(inout) :: ad_self   ! H^(x) object
      TYPE  (obs_data),   intent(in   ) :: obs       ! observation depot
      TYPE  (obs_filter), intent(inout) :: geom      ! geometry object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: model     ! kernel identifier
      integer,            intent(in   ) :: ivar      ! variable counter
      integer (i8b),      intent(out  ) :: reject    ! obs reject count
      integer,            intent(in   ) :: iTstr     ! starting timestep
      integer,            intent(in   ) :: iTend     ! ending   timestep
      integer,            intent(in   ) :: LBi, UBi  ! field I-bounds
      integer,            intent(in   ) :: LBj, UBj  ! field J-bounds
      real (r8),          intent(inout) :: ad_field(LBi:,LBj:)
!
      integer                           :: iobs, i1, j1, m, nscope
      integer                           :: i, is, ie, j, js, je
      real (r8)                         :: dfr, dfx, dfy
      real (r8)                         :: Tsteps, area_sum, resol
!
!  Clear work global 2D array to hold the adjoint of the area of
!  influence data.
!
      geom%global2d=0.0_r8
!
!  Compute area-averaged 2D model state variable at the area_averaged
!  observation locations.
!
      reject=0
      DO m=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(m)).and.                            &
            (obs%Xgrid(m).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(m)).and.                            &
     &      (obs%Ygrid(m).lt.rYmax(ng))) THEN
!
          iobs=obs%pack_mapping(m,ivar)               ! index mapping
!
          i1=INT(obs%Xgrid(m))
          j1=INT(obs%Ygrid(m))
!
!  Set the area of influence in grid units using the local average grid
!  resolution.
!
          resol=0.5_r8*(geom%dx(i1,j1)+geom%dy(i1,j1))
          nscope=INT(obs%spatialAverage/resol)
          is=MAX(geom%Imin, i1-nscope)
          ie=MIN(i1+nscope, geom%Imax)
          js=MAX(geom%Jmin, j1-nscope)
          je=MIN(j1+nscope, geom%Jmax)
!
          Tsteps=REAL(iTend-iTstr, r8)
!
!  Adjoint of average state variable over the area of influence and time
!  window.
!
          area_sum=0.0_r8               ! compute area weights
          DO j=js,je
            dfy=REAL(j-j1, r8)
            DO i=is,ie
#  ifdef MASKING
              IF (geom%mask(i,j).lt.1.0_r8) CYCLE
#  endif
              dfx=REAL(i-i1, r8)
              dfr=REAL(nscope, r8)
              IF (SQRT(dfx*dfx+dfy*dfy).le.dfr) THEN
                area_sum=area_sum+geom%area(i,j)
              END IF
            END DO
          END DO
!
          IF ((area_sum.gt.0.0_r8).and.(iic(ng).eq.iTend)) THEN
!>          self%state(iobs)=self%state(iobs)/(area_sum*Tsteps)
!>
            ad_self%state(iobs)=ad_self%state(iobs)/(area_sum*Tsteps)
          END IF
          ad_self%ObsVetting(iobs)=1
!
          IF ((iTstr.le.iic(ng)).and.(iic(ng).le.iTend)) THEN
            DO j=js,je
              dfy=REAL(j-j1, r8)
              DO i=is,ie
#  ifdef MASKING
                IF (geom%mask(i,j).lt.1.0_r8) CYCLE
#  endif
                dfx=REAL(i-i1, r8)
                dfr=REAL(nscope, r8)
                IF (SQRT(dfx*dfx+dfy*dfy).le.dfr) THEN
!>                self%state(iobs)=self%state(iobs)+                    &
!>   &                             geom%area(i,j)*geom%global2d(i,j)
!>
                  geom%global2d(i,j)=geom%global2d(i,j)+                &
     &                               geom%area(i,j)*ad_self%state(iobs)
                END IF
              END DO
            END DO
            ad_self%state(iobs)=0.0_r8
          END IF
        END IF
      END DO
!
!  Add contribution to 2D state adjoint variable.
!
      DO j=ad_self%JstrR,ad_self%JendR
        DO i=ad_self%JstrR,ad_self%JendR
          ad_field(i,j)=geom%global2d(i,j)
          geom%global2d(i,j)=0.0_r8
        END DO
      END DO
!
  10  FORMAT (/,' ROMS_HOFX::AREA_TIME_AVG2D_AD: Cannot find name = ',  &
     &        a,2x,'in obs%data observation depot object.')
!
      RETURN
      END SUBROUTINE obs_hofx_area_time_avg2d_ad
# endif
!
!-----------------------------------------------------------------------
!  It interpolates 2D state variables at the time-averaged observation
!  locations. The 2D state variable is time-averaged according to the
!  starting and ending time window values specified in the IODA-type
!  observation file. Time-averaged observations are assumed to be at
!  RHO-points (A-grid).
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
      SUBROUTINE obs_hofx_time_avg2d (self, obs, ng, model,             &
     &                                ivar, reject, iTstr, iTend,       &
     &                                LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),  intent(inout) :: self       ! H(x) object
      TYPE  (obs_data),  intent(in   ) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(out  ) :: reject     ! obs reject count
      integer,           intent(in   ) :: iTstr      ! starting timestep
      integer,           intent(in   ) :: iTend      ! ending   timestep
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(in   ) :: field(LBi:,LBj:)
!
      integer                          :: i, ic, iobs, i1, i2, j1, j2
      real (r8)                        :: p1, p2, q1, q2, wsum
      real (r8)                        :: Tsteps
      real (r8),          dimension(4) :: Hmat
!
!  Interpolate 2D model state variable at the time-averaged observation
!  locations.
!
      Tsteps=REAL(iTend-iTstr, r8)
!
      reject=0
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
     &      (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
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
!
!  Average resulting H(x) over the specified time window.
!
          IF ((iTstr.le.iic(ng)).and.(iic(ng).le.iTend)) THEN
            self%state(iobs)=self%state(iobs)+                          &
     &                       Hmat(1)*field(i1,j1)+                      &
     &                       Hmat(2)*field(i2,j1)+                      &
     &                       Hmat(3)*field(i2,j2)+                      &
     &                       Hmat(4)*field(i1,j2)
          END IF
          IF (iic(ng).eq.iTend) THEN
            self%state(iobs)=self%state(iobs)/Tsteps
          END IF
# ifdef MASKING
          IF (wsum.gt.0.0_r8) THEN
            self%obsVetting(iobs)=1
          ELSE
            self%ObsVetting(iobs)=0
            reject=reject+1
          END IF
# else
          self%ObsVetting(iobs)=1
# endif
        END IF
      END DO
!
      RETURN
      END SUBROUTINE obs_hofx_time_avg2d

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  It interpolates 3D state variables at the time-averaged observation
!  locations. The 3D state variable is time-averaged according to the
!  starting and ending time window values specified in the IODA-type
!  observation file. Time-averaged observations are assumed to be at
!  RHO-points (A-grid).
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
      SUBROUTINE obs_hofx_time_avg3d (self, obs, ng, model,             &
     &                                ivar, reject, iTstr, iTend,       &
     &                                LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),  intent(inout) :: self       ! H(x) object
      TYPE  (obs_data),  intent(inout) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(out  ) :: reject     ! obs reject count
      integer,           intent(in   ) :: iTstr      ! starting timestep
      integer,           intent(in   ) :: iTend      ! ending   timestep
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(in   ) :: field(LBi:,LBj:,:)
!
      integer                          :: i, ic, iobs, k
      integer                          :: i1, i2, j1, j2, k1, k2
      real (r8)                        :: p1, p2, q1, q2, r1, r2
      real (r8)                        :: w11, w12, w21, w22, wsum
      real (r8)                        :: Zbot, Ztop, dz
      real (r8)                        :: Tsteps
      real (r8),          dimension(8) :: Hmat
!
!  Interpolate 3D model state variable at the time-averaged observation
!  locations.
!
      Tsteps=REAL(iTend-iTstr, r8)
!
      reject=0
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
              reject=reject+1
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
!
!  Average resulting H(x) over the specified time window.
!
            IF ((iTstr.le.iic(ng)).and.(iic(ng).le.iTend)) THEN
              self%state(iobs)=self%state(iobs)+                        &
     &                         Hmat(1)*field(i1,j1,k1)+                 &
     &                         Hmat(2)*field(i2,j1,k1)+                 &
     &                         Hmat(3)*field(i2,j2,k1)+                 &
     &                         Hmat(4)*field(i1,j2,k1)+                 &
     &                         Hmat(5)*field(i1,j1,k2)+                 &
     &                         Hmat(6)*field(i2,j1,k2)+                 &
     &                         Hmat(7)*field(i2,j2,k2)+                 &
     &                         Hmat(8)*field(i1,j2,k2)
            END IF
            IF (iic(ng).eq.iTend) THEN
              self%state(iobs)=self%state(iobs)/Tsteps
            END IF
#  ifdef MASKING
            IF (wsum.gt.0.0_r8) THEN
              self%ObsVetting(iobs)=1
            ELSE
              self%ObsVetting(iobs)=0
              reject=reject+1
            END IF
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
              reject=reject+1
            END IF
#  endif
          END IF
        END IF
      END DO
!
      RETURN
      END SUBROUTINE obs_hofx_time_avg3d
# endif

# ifdef ADJOINT
!
!-----------------------------------------------------------------------
!  Adjoint of interpolating 2D state variables at the time-averaged
!  observation locations. The 2D state variable is time-averaged
!  according to the starting and ending time window values specified
!  in the IODA-type observation file. Time-averaged observations are
!  assumed to be at RHO-points (A-grid).
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
      SUBROUTINE obs_hofx_time_avg2d_ad (ad_self, obs, ng, model,       &
     &                                   ivar, reject, iTstr, iTend,    &
     &                                   LBi, UBi, LBj, UBj, ad_field)
!
      CLASS (obs_hofx),  intent(inout) :: ad_self    ! H^(x) object
      TYPE  (obs_data),  intent(in   ) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(out  ) :: reject     ! obs reject count
      integer,           intent(in   ) :: iTstr      ! starting timestep
      integer,           intent(in   ) :: iTend      ! ending   timestep
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(inout) :: ad_field(LBi:,LBj:)
!
      integer                          :: i, ic, iobs, i1, i2, j1, j2
      real (r8)                        :: p1, p2, q1, q2, wsum
      real (r8)                        :: Tsteps
      real (r8),          dimension(4) :: Hmat
!
!  Adjoint of interpolate 2D model state variable at the time-averaged
!  observation locations.
!
      Tsteps=REAL(iTend-iTstr, r8)
!
      reject=0
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
     &      (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
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
#  ifdef MASKING
          Hmat(1)=Hmat(1)*ad_self%mask(i1,j1)
          Hmat(2)=Hmat(2)*ad_self%mask(i2,j1)
          Hmat(3)=Hmat(3)*ad_self%mask(i2,j2)
          Hmat(4)=Hmat(4)*ad_self%mask(i1,j2)
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
!
!  Adjoint of average resulting H^(x) over the specified time window.
!
          IF (iic(ng).eq.iTend) THEN
!>          self%state(iobs)=self%state(iobs)/Tsteps
!>
            ad_self%state(iobs)=ad_self%state(iobs)/Tsteps
          END IF
!
          IF ((iTstr.le.iic(ng)).and.(iic(ng).le.iTend)) THEN
!>          self%state(iobs)=self%state(iobs)+                          &
!>   &                       Hmat(1)*field(i1,j1)+                      &
!>   &                       Hmat(2)*field(i2,j1)+                      &
!>   &                       Hmat(3)*field(i2,j2)+                      &
!>   &                       Hmat(4)*field(i1,j2)
!>
             ad_field(i1,j1)=ad_field(i1,j1)+                           &
     &                       Hmat(1)*ad_self%state(iobs)
             ad_field(i2,j1)=ad_field(i2,j1)+                           &
     &                       Hmat(2)*ad_self%state(iobs)
             ad_field(i2,j2)=ad_field(i2,j2)+                           &
     &                       Hmat(3)*ad_self%state(iobs)
             ad_field(i1,j2)=ad_field(i1,j2)+                           &
     &                       Hmat(4)*ad_self%state(iobs)
          END IF
          IF (iic(ng).eq.iTend) THEN
            ad_self%state(iobs)=0.0_r8
          END IF
#  ifdef MASKING
          IF (wsum.gt.0.0_r8) THEN
            ad_self%ObsVetting(iobs)=1
          ELSE
            ad_self%ObsVetting(iobs)=0
            reject=reject+1
          END IF
#  else
          ad_self%ObsVetting(iobs)=1
#  endif
        END IF
      END DO
!
      RETURN
      END SUBROUTINE obs_hofx_time_avg2d_ad

#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  Adjoint of interpolating 3D state variables at the time-averaged
!  observation locations. The 3D state variable is time-averaged
!  according to the starting and ending time window values specified
!  in the IODA-type observation file. Time-averaged observations are
!  assumed to be at RHO-points (A-grid).
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
      SUBROUTINE obs_hofx_time_avg3d_ad (ad_self, obs, ng, model,       &
     &                                   ivar, iTstr, iTend, reject,    &
     &                                   LBi, UBi, LBj, UBj, ad_field)
!
      CLASS (obs_hofx),  intent(inout) :: ad_self    ! H^(x) object
      TYPE  (obs_data),  intent(inout) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(out  ) :: reject     ! obs reject count
      integer,           intent(in   ) :: iTstr      ! starting timestep
      integer,           intent(in   ) :: iTend      ! ending   timestep
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(inout) :: ad_field(LBi:,LBj:,:)
!
      integer                          :: i, ic, iobs, k
      integer                          :: i1, i2, j1, j2, k1, k2
      real (r8)                        :: p1, p2, q1, q2, r1, r2
      real (r8)                        :: w11, w12, w21, w22, wsum
      real (r8)                        :: Zbot, Ztop, dz
      real (r8)                        :: Tsteps
      real (r8),          dimension(8) :: Hmat
!
!  Adjoint of interpolate 3D model state variable at the time-averaged
!  observation locations.
!
      Tsteps=REAL(iTend-iTstr, r8)
!
      reject=0
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
            Ztop=ad_self%depth(i1,j1,N(ng))
            Zbot=ad_self%depth(i1,j1,1    )
            IF (obs%Zgrid(i).ge.Ztop) THEN
              k1=N(ng)                         ! If shallower, assign to
              k2=N(ng)                         ! top grid cell. The
              r1=1.0_r8                        ! observation is located
              r2=0.0_r8                        ! on the upper cell half
              obs%Zgrid(i)=REAL(N(ng),r8)      ! above its middle depth.
            ELSE IF (Zbot.ge.obs%Zgrid(i)) THEN
              r1=0.0_r8                        ! If deeper, ignore.
              r2=0.0_r8
              ad_self%ObsVetting(iobs)=0
              reject=reject+1
            ELSE
              DO k=N(ng),2,-1                  ! Otherwise, interpolate
                Ztop=ad_self%depth(i1,j1,k  )  ! to fractional level
                Zbot=ad_self%depth(i1,j1,k-1)
                IF ((Ztop.gt.obs%Zgrid(i)).and.                         &
     &              (obs%Zgrid(i).ge.Zbot)) THEN
                  k1=k-1
                  k2=k
                END IF
              END DO
              dz=ad_self%depth(i1,j1,k2)-ad_self%depth(i1,j1,k1)
              r2=(obs%Zgrid(i)-ad_self%depth(i1,j1,k1))/dz
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
            Hmat(1)=Hmat(1)*ad_self%mask(i1,j1)
            Hmat(2)=Hmat(2)*ad_self%mask(i2,j1)
            Hmat(3)=Hmat(3)*ad_self%mask(i2,j2)
            Hmat(4)=Hmat(4)*ad_self%mask(i1,j2)
            Hmat(5)=Hmat(5)*ad_self%mask(i1,j1)
            Hmat(6)=Hmat(6)*ad_self%mask(i2,j1)
            Hmat(7)=Hmat(7)*ad_self%mask(i2,j2)
            Hmat(8)=Hmat(8)*ad_self%mask(i1,j2)
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
!
!  Adjoint of average resulting H^(x) over the specified time window.
!
            IF (iic(ng).eq.iTend) THEN
!>            self%state(iobs)=self%state(iobs)/Tsteps
!>
              ad_self%state(iobs)=ad_self%state(iobs)/Tsteps
            END IF
!
            IF ((iTstr.le.iic(ng)).and.(iic(ng).le.iTend)) THEN
!>            self%state(iobs)=self%state(iobs)+                        &
!>   &                         Hmat(1)*field(i1,j1,k1)+                 &
!>   &                         Hmat(2)*field(i2,j1,k1)+                 &
!>   &                         Hmat(3)*field(i2,j2,k1)+                 &
!>   &                         Hmat(4)*field(i1,j2,k1)+                 &
!>   &                         Hmat(5)*field(i1,j1,k2)+                 &
!>   &                         Hmat(6)*field(i2,j1,k2)+                 &
!>   &                         Hmat(7)*field(i2,j2,k2)+                 &
!>   &                         Hmat(8)*field(i1,j2,k2)
!>
              ad_field(i1,j1,k1)=ad_field(i1,j1,k1)+                    &
     &                           Hmat(1)*ad_self%state(iobs)
              ad_field(i2,j1,k1)=ad_field(i2,j1,k1)+                    &
     &                           Hmat(2)*ad_self%state(iobs)
              ad_field(i2,j2,k1)=ad_field(i2,j2,k1)+                    &
     &                           Hmat(3)*ad_self%state(iobs)
              ad_field(i1,j2,k1)=ad_field(i1,j2,k1)+                    &
     &                           Hmat(4)*ad_self%state(iobs)
              ad_field(i1,j1,k2)=ad_field(i1,j1,k2)+                    &
     &                           Hmat(5)*ad_self%state(iobs)
              ad_field(i2,j1,k2)=ad_field(i2,j1,k2)+                    &
     &                           Hmat(6)*ad_self%state(iobs)
              ad_field(i2,j2,k2)=ad_field(i2,j2,k2)+                    &
     &                           Hmat(7)*ad_self%state(iobs)
              ad_field(i1,j2,k2)=ad_field(i1,j2,k2)+                    &
     &                           Hmat(8)*ad_self%state(iobs)
            END IF
            IF (iic(ng).eq.iTend) THEN
              ad_self%state(iobs)=0.0_r8
            END IF
#   ifdef MASKING
            IF (wsum.gt.0.0_r8) THEN
              ad_self%ObsVetting(iobs)=1
            ELSE
              ad_self%ObsVetting(iobs)=0
              reject=reject+1
            END IF
#   else
            ad_self%ObsVetting(iobs)=1
#   endif
#   ifndef ALLOW_BOTTOM_OBS
!
!  Reject observations that lie in the lower bottom grid cell (k=1) to
!  avoid clustering due shallowing of bathymetry during smoothing and
!  coarse level half-thickness (-h < Zobs < self%depth(:,:,1)) in deep
!  water.
!
            IF ((obs%Zgrid(i).gt.0.0_r8).and.                           &
     &          (obs%Zgrid(i).le.1.0_r8)) THEN
              ad_self%ObsVetting(iobs)=0
              reject=reject+1
            END IF
#   endif
          END IF
        END IF
      END DO
!
      RETURN
      END SUBROUTINE obs_hofx_time_avg3d_ad
#  endif
# endif
!
! <><><><><><><><><><><><><><><><><><><><><><><><><> HOFX UTILITY <><><>
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
!  It writes H(x) variables into single or multiple output NetCDF
!  file(s).
!
      SUBROUTINE obs_hofx_nc_write (self, obs, ng, model)
!
      CLASS (obs_hofx),   intent(in   ) :: self   ! H(x) object
      TYPE (obs_pool),    intent(in   ) :: obs    ! observation depot
      integer,            intent(in   ) :: ng     ! nested grid
      integer,            intent(in   ) :: model  ! kernel identifier
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
!  It collects H(x) values from all tasks in the distributed-memory
!  group.
!
      SUBROUTINE obs_hofx_update (self, ng, model)
!
      CLASS (obs_hofx), intent(inout) :: self     ! H(x) object
      integer,          intent(in   ) :: ng       ! nested grid
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
