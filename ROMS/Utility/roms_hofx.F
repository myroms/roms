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
      USE mod_fourdvar
      USE mod_grid
      USE mod_iounits
      USE mod_netcdf
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
!
      USE roms_obs_mod
      USE roms_vchange_mod
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
!  Cost Function values per observations in the pool depot.
!
        TYPE (obs_cost), allocatable :: ObsCost(:)
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

# ifdef FOUR_DVAR
!
!  Data assimilation cost functtion or data penalty functional.
!
        PROCEDURE :: cost_function       => obs_hofx_cost_function
# endif
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
# ifdef BGQC
        PROCEDURE :: background_QC       => obs_hofx_background_QC
# endif
        PROCEDURE :: increment           => obs_hofx_increment
        PROCEDURE :: innovation          => obs_hofx_innovation
        PROCEDURE :: nc_write            => obs_hofx_nc_write
        PROCEDURE :: residual            => obs_hofx_residual
        PROCEDURE :: rms                 => obs_hofx_rms
        PROCEDURE :: stats               => obs_hofx_stats
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
      PUBLIC :: hofx_initialize
      PUBLIC :: hofx_run
      PUBLIC :: hofx_finalize
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
          CALL nlm_hofx(ng)%extract (obs_store, VarCha(ng),             &
     &                               geom_hofx(ng), ng, model)

# ifndef VERIFICATION
        CASE (iSTD)                                  ! background error
          CALL berr_hofx(ng)%extract (obs_store, VarCha(ng),            &
     &                                geom_hofx(ng), ng, model)
# endif
# if defined FOUR_DVAR && defined TLM_OBS
        CASE (iTLM, iRPM)                            ! tangent linear
          CALL tlm_hofx(ng)%extract (obs_store, VarCha(ng),             &
     &                               geom_hofx(ng), ng, model)
# endif
# if defined FOUR_DVAR && defined ADJOINT
        CASE (iADM)                                  ! adjoint
          CALL adm_hofx(ng)%extract (obs_store, VarCha(ng),             &
     &                               geom_hofx(ng), ng, model)
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
      IF (Master) WRITE (stdout,10) KernelString(model)
!
      SELECT CASE (model)
        CASE (iNLM)                                  ! nonlinear
          DO ng=1,Ngrids
# ifdef DISTRIBUTE
            CALL nlm_hofx(ng)%update (ng, model)
# endif
            IF (wrtNLmod(ng)) THEN
              IF (Nrun.eq.1) THEN
                CALL nlm_hofx(ng)%nc_write (obs_store, ng, model,       &
# ifdef VERIFICATION
     &                                      'hofx')
# else
     &                                      'hofxInitial')
# endif
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
!
                CALL nlm_hofx(ng)%nc_write (obs_store, ng, model,       &
     &                                      'ObsVetting')
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
# ifdef FOUR_DVAR
!
                CALL nlm_hofx(ng)%innovation (obs_store, ng, model)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
# endif
              END IF
# ifdef FOUR_DVAR
!
              IF (outer.eq.Nouter) THEN
                CALL nlm_hofx(ng)%nc_write (obs_store, ng, model,       &
     &                                      'hofxFinal')
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
!
                CALL nlm_hofx(ng)%increment (obs_store, ng, model)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
!
                CALL nlm_hofx(ng)%residual (obs_store, ng, model)
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
              END IF
# endif
            END IF
          END DO

# ifndef VERIFICATION
        CASE (iSTD)                                  ! background error
          DO ng=1,Ngrids
#  ifdef DISTRIBUTE
            CALL berr_hofx(ng)%update (ng, model)
#  endif
          END DO
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,20) model
          exit_flag=5
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Report observation-model statistics.
!
      IF (model.eq.iNLM) THEN
        DO ng=1,Ngrids
          CALL nlm_hofx(ng)%stats (obs_store, ng, model)
        END DO
      END IF
!
  10  FORMAT (/,2x,'ROMS_HOFX        - ',a,': ',                        &
     &        'Finalizing and writing H(x)')
  20  FORMAT (/,' ROMS_HOFX::HOFX_FINALIZE: illegal identifier, '       &
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
      integer                         :: i, mvars, npool
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
!  Allocate cost function object.
!
      npool=SIZE(obs_store%pool)
      IF (.not.allocated(self%ObsCost)) THEN
        allocate ( self%ObsCost(npool) )
      END IF
      DO i=1,npool
        mvars=obs_store%pool(i)%nvars
        IF (.not.allocated(self%ObsCost(i)%value)) THEN
          allocate ( self%ObsCost(i)%value(mvars) )
          self%ObsCost(i)%value=0.0_r8
        END IF
      END DO
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

# ifndef VERIFICATION
!
! <><><><><><><><><><><><><><><><><><><><><><><><> COST FUNCTIONS <><><>
!
!  It computes the data assimilation cost function according to the
!  selected 4D-Var formulation methodology in ROMS: Indirect Representer
!  Approach (R4D-Var), Restricted B-preconditioning Lanczos (RBL4D-Var),
!  Saddle-Point (SP4D-Var), or Incremental Approach (I4D-Var).
!
!  Currently, we use the dual formulation RBL4D-Var because of its
!  advanced capabilities. The R4D-Var and I4D-Var algorithms are now
!  deprecated and no futher developed, but are keep for historical
!  purposes.

#  if defined R4DVAR    || defined R4DVAR_ANA_SENSITIVITY || \
      defined TL_R4DVAR
!
!  Data penalty functional:
!
!         Jdata = transpose[H(X) - y] * O^(-1) * [H(X) - y]
!
!         y  : observations vector
!       H(X) : representer model at observation locations
!         O  : observations error covariance
!
#  elif defined RBL4DVAR || defined RBL4DVAR_ANA_SENSITIVITY || \
        defined SP4DVAR  || defined TL_RBL4DVAR
!
!  Observation cost function (Jo) as the misfit (square difference)
!  between model and observations:
!
!         Jo = 1/2 transpose[H(X) - y] * O^(-1) * [H(X) - y]
!
!         y  : observations vector
!       H(X) : nolinear model at observation locations
!         O  : observations error covariance
!
#  else
!
!  Incremental, strong constraint observation cost function (Jo):
!
!         Jo = 1/2 transpose[H(dX) - d] * O^(-1) * [H(dX) - d]
!
!         d  : innovation vector, [y - H(Xb)]
!      H(Xb) : background at observation locations
!          y : observations vector
!      H(dX) : increment at observation locations
!         O  : observations error covariance
!
#  endif

      SUBROUTINE obs_hofx_cost_function (self, obs, ng, model)
!
      CLASS (obs_hofx),   intent(inout) :: self     ! H(x) object
      CLASS (obs_pool),   intent(inout) :: obs      ! observation object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: model    ! kernel identifier
!
      integer                           :: nlocs, npool, mvars
      integer                           :: iloc, iobs, ipool, ivar
      real (r8)                         :: ObsError, ObsCost, misfit
#  ifdef DATALESS_LOOPS
      real (r8)                         :: nlmObsCost
#  endif

#  if defined R4DVAR    || defined R4DVAR_ANA_SENSITIVITY || \
      defined TL_R4DVAR
!
!  Compute data penalty functional, Jdata.
!
      IF (model.eq.iRPM) THEN
        npool=SIZE(obs%pool)
        DO ipool=1,npool
          mvars=obs%pool(ipool)%nvars
          nlocs=obs%pool(ipool)%nlocs
          ObsCost=0.0_r8
          DO ivar=1,mvars
            DO iloc=1,nlocs
              iobs=obs%pool(ipool)%pack_mapping(iloc,ivar)
              ObsError=obs%pool(ipool)%error(iloc,ivar)
              IF ((self%ObsVetting(iobs).gt.0).and.                     &
     &            (ObsError.ne.0.0_r8)) THEN
                misfit=self%state(iobs)-                                &
     &                 obs%pool(ipool)%value(iloc,ivar)
                ObsCost=ObsCost+                                        &
     &                  misfit*misfit/(ObsError*ObsError)
#   ifdef DATALESS_LOOPS
                misfit=nlm_hofx(ng)%state(iobs)-                        &
                nlmObsCost=nlmObsCost+                                  &
     &                     misfit*misfit/(ObsError*ObsError)
#   endif
              END IF
            END DO
            self%ObsCost(ipool)%value(ivar)=ObsCost
#   ifdef DATALESS_LOOPS
            nlm_hofx(ng)%ObsCost(ipool)%value(ivar)=nlmObsCost
#   endif
          END DO
        END DO
      END IF

#  elif defined RBL4DVAR || defined RBL4DVAR_ANA_SENSITIVITY || \
        defined SP4DVAR  || defined TL_RBL4DVAR
!
!  Compute observation cost function, Jo.
!
      IF (model.eq.iNLM) THEN
        npool=SIZE(obs%pool)
        DO ipool=1,npool
          mvars=obs%pool(ipool)%nvars
          nlocs=obs%pool(ipool)%nlocs
          ObsCost=0.0_r8
          DO ivar=1,mvars
            DO iloc=1,nlocs
              iobs=obs%pool(ipool)%pack_mapping(iloc,ivar)
              ObsError=obs%pool(ipool)%error(iloc,ivar)
              IF ((self%ObsVetting(iobs).gt.0).and.                     &
     &            (ObsError.ne.0.0_r8)) THEN
                misfit=self%state(iobs)-                                &
     &                 obs%pool(ipool)%value(iloc,ivar)
                ObsCost=ObsCost+                                        &
     &                  misfit*misfit/(ObsError*ObsError)
              END IF
            END DO
            self%ObsCost(ipool)%value(ivar)=ObsCost
          END DO
        END DO
      END IF
#  else
!
!  Incremental, strong constraint observation cost function, Jo.
!
      IF (model.eq.iNLM) THEN
        npool=SIZE(obs%pool)
        DO ipool=1,npool
          mvars=obs%pool(ipool)%nvars
          nlocs=obs%pool(ipool)%nlocs
          ObsCost=0.0_r8
          DO ivar=1,mvars
            DO iloc=1,nlocs
              iobs=obs%pool(ipool)%pack_mapping(iloc,ivar)
              ObsError=obs%pool(ipool)%error(iloc,ivar)
              IF ((self%ObsVetting(iobs).gt.0).and.                     &
     &            (ObsError.ne.0.0_r8)) THEN
                misfit=self%state(iobs)-                                &
     &                 obs%pool(ipool)%value(iloc,ivar)
                ObsCost=ObsCost+                                        &
     &                  misfit*misfit/(ObsError*ObsError)
              END IF
            END DO
            self%ObsCost(ipool)%value(ivar)=ObsCost
          END DO
        END DO
      ELSE IF (model.eq.iTLM) THEN
        npool=SIZE(obs%pool)
        DO ipool=1,npool
          mvars=obs%pool(ipool)%nvars
          nlocs=obs%pool(ipool)%nlocs
          ObsCost=0.0_r8
          DO ivar=1,mvars
            DO iloc=1,nlocs
              iobs=obs%pool(ipool)%pack_mapping(iloc,ivar)
              ObsError=obs%pool(ipool)%error(iloc,ivar)
              IF ((self%ObsVetting(iobs).gt.0).and.                     &
     &            (ObsError.ne.0.0_r8)) THEN
                misfit=nlm_hofx(ng)%state(iobs)+self%state(iobs)-       &
     &                 obs%pool(ipool)%value(iloc,ivar)
                ObsCost=ObsCost+                                        &
                        misfit*misfit/(ObsError*ObsError)
              END IF
            END DO
            self%ObsCost(ipool)%value(ivar)=ObsCost
          END DO
        END DO
      END IF
#  endif
!
      RETURN
      END SUBROUTINE obs_hofx_cost_function
# endif

!
! <><><><><><><><><><><><><><><><><><><><><><><><> HOFX OPERATORS <><><>
!
!  It computes the H(x) operator, which interpolates model variables to
!  observation locations in space and time. If adjoint model kernel, it
!  computes its transpose, H^(x). Area-averaged and/or time-averaged
!  operators are available.
!
      SUBROUTINE obs_hofx_extract (self, obs, vchange, geom, ng, model)
!
      CLASS (obs_hofx),   intent(inout) :: self     ! H(x) object
      CLASS (obs_pool),   intent(inout) :: obs      ! observation object
      CLASS (obs_VarCha), intent(inout) :: vchange  ! variable change
      TYPE (obs_filter),  intent(inout) :: geom     ! geometry object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: model    ! kernel identifier
!
      TYPE (obs_data),         pointer :: depot
      integer (i8b),       allocatable :: tally(:)
      integer                          :: i, iobs, mvars, nvars
      integer                          :: is, ie, mc, ObsSum, ObsReject

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
        IF (.not.allocated(tally)) THEN
          allocate ( tally(2*mvars) )
          tally=0         
        END IF
        IF (.not.allocated(op_handle)) THEN
          allocate ( op_handle(2*mvars) )
          DO i=1,2*mvars
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
!
            SELECT CASE (obs%pool(iobs)%name)
              CASE ('SSH')
                CALL obs_pool_get (obs, obs%pool(iobs)%name, depot)
                CALL self%extract2d (depot, vchange, geom, ng, model,   &
     &                               tally, obs%pool(iobs)%name)
# ifdef SOLVE3D
              CASE ('SSS', 'SST', 'uv_codar')
                CALL obs_pool_get (obs, obs%pool(iobs)%name, depot)
                CALL self%extract2d (depot, vchange, geom, ng, model,   &
     &                               tally, obs%pool(iobs)%name)
              CASE ('salt', 'temp')
                CALL obs_pool_get (obs, obs%pool(iobs)%name, depot)
                CALL self%extract3d (depot, vchange, ng, model,         &
     &                               tally, obs%pool(iobs)%name)
# endif
              CASE DEFAULT
                IF (Master) WRITE (stdout,20) obs%pool(iobs)%name, iobs
                exit_flag=5
                IF (FoundError(exit_flag, NoError,                      &
     &                         __LINE__, MyFile)) RETURN
            END SELECT
!
!  Store processing: tally(1,nvars)=Count and tally(2,nvars)=Reject.
!  To facilitate procesing the tally is a remapped vector.
!
# ifdef DISTRIBUTE
            CALL mp_reduce (ng, model, 2*mvars, tally, op_handle)
# endif
            is=obs%pool(iobs)%NstrObs
            ie=obs%pool(iobs)%NendObs
            DO i=1,nvars
              mc=1+(i-1)*2
              self%ObsCount=self%ObsCount+tally(mc)
              self%ObsReject=self%ObsReject+tally(mc+1)
            END DO
!
!  Report.
!
            IF (Master) THEN
              DO i=1,nvars
                mc=1+(i-1)*2
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
     &                            is, ie, tally(mc), tally(mc+1)
                ObsSum=ObsSum+tally(mc)
                ObsReject=ObsReject+tally(mc+1)
              END DO
            END IF
          END IF
        END DO
        IF (Master) WRITE (stdout,40) ObsSum, ObsReject,                &
     &                                self%ObsCount, self%ObsReject
!
!  Deallocate local arrays.
!
        IF (allocated(tally))     deallocate (tally)
        IF (allocated(op_handle)) deallocate (op_handle)
      END IF
!
  10  FORMAT (/,' Number of State Observations Processed:',2x,          &
     &        'ObsTime = ',f12.4,',',t68,a,/,/,                         &
     &        3x,'Variable',15x,'IstrObs',4x,'IendObs',6x,'Count',      &
     &        3x,'Rejected',/)
  20  FORMAT (/,' ROMS_HOFX::extract: Cannot find name = ',a,2x,        &
     &        'in obs%pool(',i0,').')
  30  FORMAT (3x,a,t23,4(1x,i10))
  40  FORMAT (/,3x,'Total',t45,2(1x,i10),                               &
     &        /,3x,'Obs Tally',t45,2(1x,i10),/)
!
      RETURN
      END SUBROUTINE obs_hofx_extract
!
!-----------------------------------------------------------------------
!  It computes 2D control variables H(x) operator. If adjoint kernel,
!  model=iADM, it computes the transport operator, H^(x). Area-averaged
!  and/or time-averaged operators are available.
!
      SUBROUTINE obs_hofx_extract2d (self, obs, vchange, geom,          &
     &                               ng, model, tally, name)
!
      CLASS (obs_hofx),   intent(inout) :: self      ! H(x) object
      CLASS (obs_data),   intent(inout) :: obs       ! observation depot
      CLASS (obs_VarCha), intent(inout) :: vchange   ! variable change
      TYPE (obs_filter),  intent(inout) :: geom      ! geometry object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: model     ! kernel identifier
      integer (i8b),      intent(inout) :: tally(:)  ! obs count/reject
      character (len=*),  intent(in   ) :: name      ! obs pool name
!
      logical                           :: IsAreaTimeAveraged
      logical                           :: IsAreaAveraged
      logical                           :: IsTimeAveraged
      integer                           :: LBi, UBi, LBj, UBj
      integer                           :: ivar, tile
# ifdef SOLVE3D
      integer                           :: itrc
# endif
!
      real (r8),                pointer :: field(:,:) => NULL()

      character (len=*),      parameter :: MyFile =                     &
     &    __FILE__//", roms_hofx_extract2d"
!
      SourceFile=MyFile
!
!  Initialize.
!
      tile=self%tile
      LBi=self%LBi
      UBi=self%UBi
      LBj=self%LBj
      UBj=self%UBj
!
!  Compute 2D H(x) operator.
!
      VAR_LOOP : DO ivar=1,obs%nvars
!
!  Get pointer (vchange%field2d) of the state variable associated with
!  the observations. If necessary apply the required variable changes.
!
        CALL vchange%apply2d (ng, tile, model, obs%vname(ivar))
!
!  Set switch for filters.
!
        IsAreaTimeAveraged = obs%IsAreaAveraged.and.                    &
     &                       obs%IsTimeAveraged
        IsAreaAveraged     = .not.obs%IsTimeAveraged.and.               &
     &                            obs%IsAreaAveraged
        IsTimeAveraged     = .not.obs%IsAreaAveraged.and.               &
     &                            obs%IsTimeAveraged
!
!  Compute 2D control variable H(x) operator.
!
        IF (model.ne.iADM) THEN
          IF (IsAreaTimeAveraged) THEN
            CALL self%area_time_avg2d (obs, geom, ng, model,            &
     &                                 ivar, tally,                     &
     &                                 LBi, UBi, LBj, UBj,              &
     &                                 vchange%field2d)
          ELSE IF (IsAreaAveraged) THEN
            CALL self%area_avg2d (obs, geom, ng, model,                 &
     &                            ivar, tally,                          &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vchange%field2d)
          ELSE IF (IsTimeAveraged) THEN
            CALL self%time_avg2d (obs, ng, model,                       &
     &                            ivar, tally,                          &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vchange%field2d)
          ELSE
            CALL self%interp2d (obs, ng, model,                         &
     &                          ivar, tally,                            &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vchange%field2d)
          END IF
# ifdef ADJOINT
        ELSE IF (model.eq.iADM) THEN                ! transpose, H^(x)
          IF (IsAreaTimeAveraged) THEN
            CALL self%area_time_avg2d_ad (obs, geom, ng, model,         &
     &                                    ivar, tally,                  &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    vchange%field2d)

          ELSE IF (IsAreaAveraged) THEN
            CALL self%area_avg2d_ad (obs, geom, ng, model,              &
     &                               ivar, tally,                       &
     &                               LBi, UBi, LBj, UBj,                &
     &                               vchange%field2d)
          ELSE IF (IsTimeAveraged) THEN
            CALL self%time_avg2d_ad (obs, ng, model,                    &
     &                               ivar, tally,                       &
     &                               LBi, UBi, LBj, UBj,                &
     &                               vchange%field2d)
          ELSE
            CALL self%interp2d_ad (obs, ng, model,                      &
     &                             ivar, tally,                         &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             vchange%field2d)
          END IF
!
!  Adjoint of variables changes, if any.
!
          CALL vchange%apply2d_ad (ng, tile, model, obs%vname(ivar))
# endif
        END IF
      END DO VAR_LOOP
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
      SUBROUTINE obs_hofx_extract3d (self, obs, vchange,                &
     &                               ng, model, tally, name)
!
      CLASS (obs_hofx),   intent(inout) :: self      ! H(x) object
      CLASS (obs_data),   intent(inout) :: obs       ! observation depot
      CLASS (obs_VarCha), intent(inout) :: vchange   ! variable change
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: model     ! kernel identifier
      integer (i8b),      intent(inout) :: tally(:)  ! obs count/reject
      character (len=*),  intent(in   ) :: name      ! obs pool name
!
      logical                           :: IsTimeAveraged
      integer                           :: LBi, UBi, LBj, UBj, tile
      integer                           :: iTstr, iTend
      integer                           :: itrc, ivar
!
      character (len=*),      parameter :: MyFile =                     &
     &    __FILE__//", roms_hofx_extract3d"
!
      SourceFile=MyFile
!
!  Initialize.
!
      tile=self%tile
      LBi=self%LBi
      UBi=self%UBi
      LBj=self%LBj
      UBj=self%UBj
!
!  Compute the 3D H(x) operator.
!
      VAR_LOOP : DO ivar=1,obs%nvars
!
!  Get pointer (vchange%field3d) of the state variable associated with
!  the observations. If necessary apply the required variable changes.
!
        CALL vchange%apply3d (ng, tile, model, obs%vname(ivar))
!
!  Set switch for filters.
!
        IsTimeAveraged = .not.obs%IsAreaAveraged.and.                   &
     &                        obs%IsTimeAveraged
!
!  Compute 3D control variable H(x) operator.
!
        IF (model.ne.iADM) THEN
          IF (IsTimeAveraged) THEN
            CALL self%time_avg3d (obs, ng, model,                       &
     &                            ivar, tally,                          &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            vchange%field3d)
          ELSE
            CALL self%interp3d (obs, ng, model,                         &
     &                          ivar, tally,                             &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          vchange%field3d)

          END IF
#  ifdef ADJOINT
        ELSE IF (model.eq.iADM) THEN                ! transpose, H^(x)
          IF (IsTimeAveraged) THEN
            CALL self%time_avg3d_ad (obs, ng, model,                    &
     &                               ivar, tally,                       &
     &                               LBi, UBi, LBj, UBj,                &
     &                               vchange%field3d)

          ELSE
            CALL self%interp3d_ad (obs, ng, model,                      &
     &                             ivar, tally,                         &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             vchange%field3d)

          END IF
!
!  Adjoint of variables changes, if any.
!
          CALL vchange%apply3d_ad (ng, tile, model, obs%vname(ivar))
#  endif
        END IF
      END DO VAR_LOOP
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
     &                              ivar, tally,                        &
     &                              LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),  intent(inout) :: self       ! H(x) object
      CLASS (obs_data),  intent(in   ) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(inout) :: tally(:)   ! obs count/reject
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(in   ) :: field(LBi:,LBj:)
!
      integer                          :: i, ic, iobs, mc
      integer                          :: i1, i2, j1, j2
      real (r8)                        :: p1, p2, q1, q2, wsum
      real (r8),          dimension(4) :: Hmat
!
!  Interpolate 2D model state variable at the observation locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
            (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(i,ivar)            ! index mapping
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
            tally(mc+1)=tally(mc+1)+1              ! obs rejected
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
     &                              ivar, tally,                        &
     &                              LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),  intent(inout) :: self       ! H(x) object
      CLASS (obs_data),  intent(inout) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(inout) :: tally(:)   ! obs count/reject
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(in   ) :: field(LBi:,LBj:,:)
!
      integer                          :: i, ic, iobs, k, mc
      integer                          :: i1, i2, j1, j2, k1, k2
      real (r8)                        :: p1, p2, q1, q2, r1, r2
      real (r8)                        :: w11, w12, w21, w22, wsum
      real (r8)                        :: Zbot, Ztop, dz
      real (r8),          dimension(8) :: Hmat
!
!  Interpolate 3D model state variable at the observation locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
            (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(i,ivar)            ! index mapping
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
              tally(mc+1)=tally(mc+1)+1        ! obs rejected
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
              tally(mc+1)=tally(mc+1)+1            ! obs rejected
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
              tally(mc+1)=tally(mc+1)+1            ! obs rejected
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
     &                                 ivar, tally,                     &
     &                                 LBi, UBi, LBj, UBj, ad_field)
!
      CLASS (obs_hofx),  intent(inout) :: ad_self    ! H^(x) object
      CLASS (obs_data),  intent(in   ) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(inout) :: tally(:)   ! obs count/reject
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(inout) :: ad_field(LBi:,LBj:)
!
      integer                          :: i, ic, iobs, mc
      integer                          :: i1, i2, j1, j2
      real (r8)                        :: p1, p2, q1, q2, wsum
      real (r8),          dimension(4) :: Hmat
!
!  Interpolate 2D model state variable at the observation locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
            (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(i,ivar)            ! index mapping
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
            tally(mc+1)=tally(mc+1)+1              ! obs rejected
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
     &                                 ivar, tally,                     &
     &                                 LBi, UBi, LBj, UBj, ad_field)
!
      CLASS (obs_hofx),  intent(inout) :: ad_self    ! H^(x) object
      CLASS (obs_data),  intent(inout) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(inout) :: tally(:)   ! obs count/reject
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(inout) :: ad_field(LBi:,LBj:,:)
!
      integer                          :: i, ic, iobs, k, mc
      integer                          :: i1, i2, j1, j2, k1, k2
      real (r8)                        :: p1, p2, q1, q2, r1, r2
      real (r8)                        :: w11, w12, w21, w22, wsum
      real (r8)                        :: Zbot, Ztop, dz
      real (r8),          dimension(8) :: Hmat
!
!  Interpolate 2D model state variable at the observation locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
            (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(i,ivar)            ! index mapping
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
              tally(mc+1)=tally(mc+1)+1        ! obs rejected
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
              tally(mc+1)=tally(mc+1)+1            ! obs rejected
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
     &                                ivar, tally,                      &
     &                                LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),   intent(inout) :: self      ! H(x) object
      CLASS (obs_data),   intent(in   ) :: obs       ! observation depot
      TYPE  (obs_filter), intent(inout) :: geom      ! geometry
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: model     ! kernel identifier
      integer,            intent(in   ) :: ivar      ! variable counter
      integer (i8b),      intent(inout) :: tally(:)  ! obs count/reject
      integer,            intent(in   ) :: LBi, UBi  ! field I-bounds
      integer,            intent(in   ) :: LBj, UBj  ! field J-bounds
      real (r8),          intent(in   ) :: field(LBi:,LBj:)
!
      integer                           :: iobs, i1, j1, m, mc, nscope
      integer                           :: i, is, ie, j, js, je
      real (r8)                         :: dfr, dfx, dfy
      real (r8)                         :: area_sum, resol
      real (r8)                         :: my_field
!
      character (len=*),      parameter :: MyFile =                     &
     &    __FILE__//", roms_hofx_area_avg2d"
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
      DO m=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(m)).and.                            &
            (obs%Xgrid(m).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(m)).and.                            &
     &      (obs%Ygrid(m).lt.rYmax(ng))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(m,ivar)            ! index mapping
!
          i1=INT(obs%Xgrid(m))
          j1=INT(obs%Ygrid(m))
!
!  Set the area of influence in grid units using the local average grid
!  resolution.
!
          resol=0.5_r8*(geom%dx(i1,j1)+geom%dy(i1,j1))
          nscope=INT(obs%spatialAverage(ivar)/resol)
          is=MAX(geom%Imin, i1-nscope)
          ie=MIN(i1+nscope, geom%Imax)
          js=MAX(geom%Jmin, j1-nscope)
          je=MIN(j1+nscope, geom%Jmax)
!
!  Average state variable over the area of influence. If model standard
!  deviations, square their values to obtain variances. Afterwards, take
!  the square root to convert back to standard deviations.
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
                IF (model.eq.iSTD) THEN     ! sum over the square values
# ifdef DISTRIBUTE
                  my_field=geom%global2d(i,j)*geom%global2d(i,j)
# else
                  my_field=field(i,j)*field(i,j)
# endif
                ELSE
# ifdef DISTRIBUTE
                  my_field=geom%global2d(i,j)
# else
                  my_field=field(i,j)
# endif
                END IF
                self%state(iobs)=self%state(iobs)+                      &
     &                           geom%area(i,j)*my_field
              END IF
            END DO
          END DO
!
!  Compute area-weighted average.
!
          IF (area_sum.gt.0.0_r8) THEN
            self%state(iobs)=self%state(iobs)/area_sum
            IF (model.eq.iSTD) THEN         ! take the square root
              self%state(iobs)=SQRT(self%state(iobs))
            END IF
            self%ObsVetting(iobs)=1
          ELSE
            self%ObsVetting(iobs)=0
            tally(mc+1)=tally(mc+1)+1              ! obs rejected
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
     &                                     ivar, tally,                 &
     &                                     LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),   intent(inout) :: self      ! H(x) object
      CLASS (obs_data),   intent(in   ) :: obs       ! observation depot
      TYPE  (obs_filter), intent(inout) :: geom      ! geometry object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: model     ! kernel identifier
      integer,            intent(in   ) :: ivar      ! variable counter
      integer (i8b),      intent(inout) :: tally(:)  ! obs count/reject
      integer,            intent(in   ) :: LBi, UBi  ! field I-bounds
      integer,            intent(in   ) :: LBj, UBj  ! field J-bounds
      real (r8),          intent(in   ) :: field(LBi:,LBj:)
!
      integer                           :: iobs, i1, j1, m, mc, nscope
      integer                           :: i, is, ie, j, js, je
      integer                           :: iTstr, iTend
      real (r8)                         :: dfr, dfx, dfy
      real (r8)                         :: Tsteps, area_sum, resol
      real (r8)                         :: my_field
!
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
      DO m=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(m)).and.                            &
            (obs%Xgrid(m).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(m)).and.                            &
     &      (obs%Ygrid(m).lt.rYmax(ng)).and.                            &
     &      (obs%TimeStr(m).le.time(ng)).and.                           &
     &      (time(ng).le.obs%TimeEnd(m))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(m,ivar)            ! index mapping
!
          i1=INT(obs%Xgrid(m))
          j1=INT(obs%Ygrid(m))
!
!  Set the area of influence in grid units using the local average grid
!  resolution.
!
          resol=0.5_r8*(geom%dx(i1,j1)+geom%dy(i1,j1))
          nscope=INT(obs%spatialAverage(ivar)/resol)
          is=MAX(geom%Imin, i1-nscope)
          ie=MIN(i1+nscope, geom%Imax)
          js=MAX(geom%Jmin, j1-nscope)
          je=MIN(j1+nscope, geom%Jmax)
!
          iTstr=MAX(ntstart(ng),                                        &
     &              INT((obs%TimeStr(m)-dstart*day2sec)/dt(ng))+1)
          iTend=MIN(ntend(ng),                                          &
     &              INT((obs%TimeEnd(m)-dstart*day2sec)/dt(ng)))
          Tsteps=REAL(iTend-iTstr, r8)
!
!  Average state variable over the area of influence and time window.
!  If model standard deviations, square their values to obtain
!  variances. Afterwards, take the square root to convert back to
!  standard deviations.
!
          area_sum=0.0_r8
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
                IF (model.eq.iSTD) THEN   ! sum over the square values
# ifdef DISTRIBUTE
                  my_field=geom%global2d(i,j)*geom%global2d(i,j)
# else
                  my_field=field(i,j)*field(i,j)
# endif
                ELSE
# ifdef DISTRIBUTE
                  my_field=geom%global2d(i,j)
# else
                  my_field=field(i,j)
# endif
                END IF
                self%state(iobs)=self%state(iobs)+                      &
     &                           geom%area(i,j)*my_field
              END IF
            END DO
          END DO
!
!  Compute area-weighted average.
!
          IF ((area_sum.gt.0.0_r8).and.(iic(ng).eq.iTend)) THEN
            self%state(iobs)=self%state(iobs)/(area_sum*Tsteps)
            IF (model.eq.iSTD) THEN         ! take the square root
              self%state(iobs)=SQRT(self%state(iobs))
            END IF
          END IF
          self%ObsVetting(iobs)=1
        END IF
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
     &                                   model, ivar, tally,            &
     &                                   LBi, UBi, LBj, UBj, ad_field)
!
      CLASS (obs_hofx),   intent(inout) :: ad_self   ! H'(x) object
      CLASS (obs_data),   intent(in   ) :: obs       ! observation depot
      TYPE  (obs_filter), intent(inout) :: geom      ! geometry object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: model     ! kernel identifier
      integer,            intent(in   ) :: ivar      ! variable counter
      integer (i8b),      intent(inout) :: tally(:)  ! obs count/reject
      integer,            intent(in   ) :: LBi, UBi  ! field I-bounds
      integer,            intent(in   ) :: LBj, UBj  ! field J-bounds
      real (r8),          intent(inout) :: ad_field(LBi:,LBj:)
!
      integer                         :: iobs, i1, j1, m, mc, nscope
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
      DO m=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(m)).and.                            &
            (obs%Xgrid(m).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(m)).and.                            &
     &      (obs%Ygrid(m).lt.rYmax(ng))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(m,ivar)            ! index mapping
!
          i1=INT(obs%Xgrid(m))
          j1=INT(obs%Ygrid(m))
!
!  Set the area of influence in grid units using the local average grid
!  resolution.
!
          resol=0.5_r8*(geom%dx(i1,j1)+geom%dy(i1,j1))
          nscope=INT(obs%spatialAverage(ivar)/resol)
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
            tally(mc+1)=tally(mc+1)+1              ! obs rejected
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
     &                                        ng, model, ivar, tally,   &
     &                                        LBi, UBi, LBj, UBj,       &
     &                                        ad_field)
!
      CLASS (obs_hofx),   intent(inout) :: ad_self   ! H^(x) object
      CLASS (obs_data),   intent(in   ) :: obs       ! observation depot
      TYPE  (obs_filter), intent(inout) :: geom      ! geometry object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: model     ! kernel identifier
      integer,            intent(in   ) :: ivar      ! variable counter
      integer (i8b),      intent(inout) :: tally(:)  ! obs count/reject
      integer,            intent(in   ) :: LBi, UBi  ! field I-bounds
      integer,            intent(in   ) :: LBj, UBj  ! field J-bounds
      real (r8),          intent(inout) :: ad_field(LBi:,LBj:)
!
      integer                           :: iobs, i1, j1, m, mc, nscope
      integer                           :: i, is, ie, j, js, je
      integer                           :: iTstr, iTend
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
      DO m=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(m)).and.                            &
            (obs%Xgrid(m).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(m)).and.                            &
     &      (obs%Ygrid(m).lt.rYmax(ng)).and.                            &
     &      (obs%TimeStr(m).le.time(ng)).and.                           &
     &      (time(ng).le.obs%TimeEnd(m))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(m,ivar)            ! index mapping
!
          i1=INT(obs%Xgrid(m))
          j1=INT(obs%Ygrid(m))
!
!  Set the area of influence in grid units using the local average grid
!  resolution.
!
          resol=0.5_r8*(geom%dx(i1,j1)+geom%dy(i1,j1))
          nscope=INT(obs%spatialAverage(ivar)/resol)
          is=MAX(geom%Imin, i1-nscope)
          ie=MIN(i1+nscope, geom%Imax)
          js=MAX(geom%Jmin, j1-nscope)
          je=MIN(j1+nscope, geom%Jmax)
!
          iTstr=MAX(ntstart(ng),                                        &
     &              INT((obs%TimeStr(m)-dstart*day2sec)/dt(ng))+1)
          iTend=MIN(ntend(ng),                                          &
     &              INT((obs%TimeEnd(m)-dstart*day2sec)/dt(ng)))
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
     &                                ivar, tally,                      &
     &                                LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),  intent(inout) :: self       ! H(x) object
      CLASS (obs_data),  intent(in   ) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(inout) :: tally(:)   ! obs count/reject
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(in   ) :: field(LBi:,LBj:)
!
      integer                          :: i, ic, iobs, mc
      integer                          :: i1, i2, j1, j2
      integer                          :: iTstr, iTend
      real (r8)                        :: p1, p2, q1, q2, wsum
      real (r8)                        :: Tsteps, my_field
      real (r8),          dimension(4) :: Hmat
!
!  Interpolate 2D model state variable at the time-averaged observation
!  locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
     &      (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng)).and.                            &
     &      (obs%TimeStr(i).le.time(ng)).and.                           &
     &      (time(ng).le.obs%TimeEnd(i))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(i,ivar)            ! index mapping
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
!  Average resulting H(x) over the specified time window. If model
!  standard deviations, square their values to obtain variances.
!  Afterwards, take the square root to convert back to standard
!  deviations. Probably, unnecessary since it is time independent.
!
          iTstr=MAX(ntstart(ng),                                        &
     &              INT((obs%TimeStr(i)-dstart*day2sec)/dt(ng))+1)
          iTend=MIN(ntend(ng),                                          &
     &              INT((obs%TimeEnd(i)-dstart*day2sec)/dt(ng)))
          Tsteps=REAL(iTend-iTstr, r8)
!
          IF (model.eq.iSTD) THEN           ! use squared values
            self%state(iobs)=self%state(iobs)+                          &
     &                       Hmat(1)*field(i1,j1)*field(i1,j1)+         &
     &                       Hmat(2)*field(i2,j1)*field(i2,j1)+         &
     &                       Hmat(3)*field(i2,j2)*field(i2,j2)+         &
     &                       Hmat(4)*field(i1,j2)*field(i1,j2)
          ELSE
            self%state(iobs)=self%state(iobs)+                          &
     &                       Hmat(1)*field(i1,j1)+                      &
     &                       Hmat(2)*field(i2,j1)+                      &
     &                       Hmat(3)*field(i2,j2)+                      &
     &                       Hmat(4)*field(i1,j2)
          END IF
          IF (iic(ng).eq.iTend) THEN
            self%state(iobs)=self%state(iobs)/Tsteps
            IF (model.eq.iSTD) THEN         ! take the square root
              self%state(iobs)=SQRT(self%state(iobs))
            END IF
          END IF
# ifdef MASKING
          IF (wsum.gt.0.0_r8) THEN
            self%obsVetting(iobs)=1
          ELSE
            self%ObsVetting(iobs)=0
            tally(mc+1)=tally(mc+1)+1              ! obs rejected
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
     &                                ivar, tally,                      &
     &                                LBi, UBi, LBj, UBj, field)
!
      CLASS (obs_hofx),  intent(inout) :: self       ! H(x) object
      CLASS (obs_data),  intent(inout) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(inout) :: tally(:)   ! obs count/reject
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(in   ) :: field(LBi:,LBj:,:)
!
      integer                          :: i, ic, iobs, k, mc
      integer                          :: i1, i2, j1, j2, k1, k2
      integer                          :: iTstr, iTend
      real (r8)                        :: p1, p2, q1, q2, r1, r2
      real (r8)                        :: w11, w12, w21, w22, wsum
      real (r8)                        :: Zbot, Ztop, dz
      real (r8)                        :: Tsteps
      real (r8),          dimension(8) :: Hmat
!
!  Interpolate 3D model state variable at the time-averaged observation
!  locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
            (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng)).and.                            &
     &      (obs%TimeStr(i).le.time(ng)).and.                           &
     &      (time(ng).le.obs%TimeEnd(i))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(i,ivar)            ! index mapping
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
              tally(mc+1)=tally(mc+1)+1        ! obs rejected
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
!  Average resulting H(x) over the specified time window.  If model
!  standard deviations, square their values to obtain variances.
!  Afterwards, take the square root to convert back to standard
!  deviations. Probably, unnecessary since it is time independent.
!
            iTstr=MAX(ntstart(ng),                                      &
     &                INT((obs%TimeStr(i)-dstart*day2sec)/dt(ng))+1)
            iTend=MIN(ntend(ng),                                        &
     &                INT((obs%TimeEnd(i)-dstart*day2sec)/dt(ng)))
            Tsteps=REAL(iTend-iTstr, r8)
!
            IF (model.eq.iSTD) THEN         ! use squared values
              self%state(iobs)=self%state(iobs)+                        &
     &                         Hmat(1)*field(i1,j1,k1)*field(i1,j1,k1)+ &
     &                         Hmat(2)*field(i2,j1,k1)*field(i2,j1,k1)+ &
     &                         Hmat(3)*field(i2,j2,k1)*field(i2,j2,k1)+ &
     &                         Hmat(4)*field(i1,j2,k1)*field(i1,j2,k1)+ &
     &                         Hmat(5)*field(i1,j1,k2)*field(i1,j1,k2)+ &
     &                         Hmat(6)*field(i2,j1,k2)*field(i2,j1,k2)+ &
     &                         Hmat(7)*field(i2,j2,k2)*field(i2,j2,k2)+ &
     &                         Hmat(8)*field(i1,j2,k2)*field(i1,j2,k2)
            ELSE
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
              IF (model.eq.iSTD) THEN        ! take the square root
                self%state(iobs)=SQRT(self%state(iobs))
              END IF
            ELSE
              self%ObsVetting(iobs)=0
              tally(mc+1)=tally(mc+1)+1            ! obs rejected
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
              tally(mc+1)=tally(mc+1)+1            ! obs rejected
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
     &                                   ivar, tally,                   &
     &                                   LBi, UBi, LBj, UBj, ad_field)
!
      CLASS (obs_hofx),  intent(inout) :: ad_self    ! H^(x) object
      CLASS (obs_data),  intent(in   ) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(inout) :: tally(:)   ! obs count/reject
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(inout) :: ad_field(LBi:,LBj:)
!
      integer                          :: i, ic, iobs, mc
      integer                          :: i1, i2, j1, j2
      integer                          :: iTstr, iTend
      real (r8)                        :: p1, p2, q1, q2, wsum
      real (r8)                        :: Tsteps
      real (r8),          dimension(4) :: Hmat
!
!  Adjoint of interpolate 2D model state variable at the time-averaged
!  observation locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
     &      (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng)).and.                            &
     &      (obs%TimeStr(i).le.time(ng)).and.                           &
     &      (time(ng).le.obs%TimeEnd(i))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(i,ivar)            ! index mapping
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
          iTstr=MAX(ntstart(ng),                                        &
     &              INT((obs%TimeStr(i)-dstart*day2sec)/dt(ng))+1)
          iTend=MIN(ntend(ng),                                          &
     &              INT((obs%TimeEnd(i)-dstart*day2sec)/dt(ng)))
          Tsteps=REAL(iTend-iTstr, r8)
!
          IF (iic(ng).eq.iTend) THEN
!>          self%state(iobs)=self%state(iobs)/Tsteps
!>
            ad_self%state(iobs)=ad_self%state(iobs)/Tsteps
          END IF
!
!>        self%state(iobs)=self%state(iobs)+                            &
!>   &                     Hmat(1)*field(i1,j1)+                        &
!>   &                     Hmat(2)*field(i2,j1)+                        &
!>   &                     Hmat(3)*field(i2,j2)+                        &
!>   &                     Hmat(4)*field(i1,j2)
!>
           ad_field(i1,j1)=ad_field(i1,j1)+                             &
     &                     Hmat(1)*ad_self%state(iobs)
           ad_field(i2,j1)=ad_field(i2,j1)+                             &
     &                     Hmat(2)*ad_self%state(iobs)
           ad_field(i2,j2)=ad_field(i2,j2)+                             &
     &                     Hmat(3)*ad_self%state(iobs)
           ad_field(i1,j2)=ad_field(i1,j2)+                             &
     &                     Hmat(4)*ad_self%state(iobs)
          IF (iic(ng).eq.iTend) THEN
            ad_self%state(iobs)=0.0_r8
          END IF
#  ifdef MASKING
          IF (wsum.gt.0.0_r8) THEN
            ad_self%ObsVetting(iobs)=1
          ELSE
            ad_self%ObsVetting(iobs)=0
            tally(mc+1)=tally(mc+1)+1              ! obs rejected
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
     &                                   ivar, tally,                   &
     &                                   LBi, UBi, LBj, UBj, ad_field)
!
      CLASS (obs_hofx),  intent(inout) :: ad_self    ! H^(x) object
      CLASS (obs_data),  intent(inout) :: obs        ! observation depot
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: ivar       ! variable counter
      integer (i8b),     intent(inout) :: tally(:)   ! obs count/reject
      integer,           intent(in   ) :: LBi, UBi   ! field I-bounds
      integer,           intent(in   ) :: LBj, UBj   ! field J-bounds
      real (r8),         intent(inout) :: ad_field(LBi:,LBj:,:)
!
      integer                          :: i, ic, iobs, k, mc
      integer                          :: i1, i2, j1, j2, k1, k2
      integer                          :: iTstr, iTend
      real (r8)                        :: p1, p2, q1, q2, r1, r2
      real (r8)                        :: w11, w12, w21, w22, wsum
      real (r8)                        :: Zbot, Ztop, dz
      real (r8)                        :: Tsteps
      real (r8),          dimension(8) :: Hmat
!
!  Adjoint of interpolate 3D model state variable at the time-averaged
!  observation locations.
!
      DO i=obs%NstrObs,obs%NendObs
        IF ((rXmin(ng).le.obs%Xgrid(i)).and.                            &
            (obs%Xgrid(i).lt.rXmax(ng)).and.                            &
     &      (rYmin(ng).le.obs%Ygrid(i)).and.                            &
     &      (obs%Ygrid(i).lt.rYmax(ng)).and.                            &
     &      (obs%TimeStr(i).le.time(ng)).and.                           &
     &      (time(ng).le.obs%TimeStepEnd(i))) THEN
!
          mc=1+(ivar-1)*2                          ! count/reject index
          tally(mc)=tally(mc)+1                    ! obs counted
          iobs=obs%pack_mapping(i,ivar)            ! index mapping
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
              tally(mc+1)=tally(mc+1)+1        ! obs reject
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
            iTstr=MAX(ntstart(ng),                                      &
     &                INT((obs%TimeStr(i)-dstart*day2sec)/dt(ng))+1)
            iTend=MIN(ntend(ng),                                        &
     &                INT((obs%TimeEnd(i)-dstart*day2sec)/dt(ng)))
            Tsteps=REAL(iTend-iTstr, r8)
!
            IF (iic(ng).eq.iTend) THEN
!>            self%state(iobs)=self%state(iobs)/Tsteps
!>
              ad_self%state(iobs)=ad_self%state(iobs)/Tsteps
            END IF
!
!>          self%state(iobs)=self%state(iobs)+                          &
!>   &                       Hmat(1)*field(i1,j1,k1)+                   &
!>   &                       Hmat(2)*field(i2,j1,k1)+                   &
!>   &                       Hmat(3)*field(i2,j2,k1)+                   &
!>   &                       Hmat(4)*field(i1,j2,k1)+                   &
!>   &                       Hmat(5)*field(i1,j1,k2)+                   &
!>   &                       Hmat(6)*field(i2,j1,k2)+                   &
!>   &                       Hmat(7)*field(i2,j2,k2)+                   &
!>   &                       Hmat(8)*field(i1,j2,k2)
!>
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
            IF (iic(ng).eq.iTend) THEN
              ad_self%state(iobs)=0.0_r8
            END IF
#   ifdef MASKING
            IF (wsum.gt.0.0_r8) THEN
              ad_self%ObsVetting(iobs)=1
            ELSE
              ad_self%ObsVetting(iobs)=0
              tally(mc+1)=tally(mc+1)+1            ! obs rejected
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
              tally(mc+1)=tally(mc+1)+1            ! obs rejected
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

# ifdef BGQC
!
!-----------------------------------------------------------------------
!  It rejects observations that fail the background quality contol by
!  setting ObsVetting to zero.
!
      SUBROUTINE obs_hofx_background_QC (self, obs, berr, ng, model)
!
      CLASS (obs_hofx), intent(inout) :: self   ! H(x) object
      CLASS (obs_pool), intent(in   ) :: obs    ! observation depot
      CLASS (obs_hofx), intent(in   ) :: Berr   ! H(x) background error
      integer,          intent(in   ) :: ng     ! nested grid
      integer,          intent(in   ) :: model  ! kernel identifier
!
      integer                         :: nlocs, npool, mvars
      integer                         :: i, iloc, iobs, ivar
      real (r8)                       :: BerrInv, obs_err, obs_val
      real (r8)                       :: df1, df2, thresh
      real (r8),          allocatable :: Bthresh(:)
!
!  Perform background quality control of the observations. Here, "self"
!  is the nonlinear H(x), except in the inverse representers alogorithm
!  where H(x) is for the finite amplitude tangent linear model.
!
      npool=SIZE(obs%pool)
      DO i=1,npool
        mvars=obs%pool(i)%nvars
        nlocs=obs%pool(i)%nlocs
!
        IF (.not.allocated(Bthresh)) THEN
          allocate ( Bthresh(nlocs) )
          Bthresh=bgqc_large
        END IF
!
        DO ivar=1,mvars
          CALL obs%pool(i)%qc_threshold (ng, model, ivar, Bthresh)
          DO iloc=1,nlocs
            iobs=obs%pool(i)%pack_mapping(iloc,ivar)
            IF (self%ObsVetting(iobs).ne.0) THEN
              obs_val=obs%pool(i)%value(iloc,ivar)
              obs_err=obs%pool(i)%error(iloc,ivar)
              BerrInv=1.0_r8/Berr%state(iobs)
!
              df1=(obs_val-self%state(iobs))*BerrInv
              df2=(1.0_r8/obs_err)/(BerrInv*BerrInv)
!
              thresh=Bthresh(iloc)*(1.0_r8+df2)
              IF (df1*df1.gt.thresh) THEN
                self%ObsVetting(iobs)=0
              END IF
            END IF
          END DO
        END DO
        IF (allocated(Bthresh)) deallocate (Bthresh)
      END DO
!
      RETURN
      END SUBROUTINE obs_hofx_background_QC
# endif
!
!-----------------------------------------------------------------------
!  It computes and writes the increment (analysis minus background) at
!  the observation locations into multiple output NetCDF-4 files.
!
      SUBROUTINE obs_hofx_increment (self, obs, ng, model)
!
      CLASS (obs_hofx), intent(in   ) :: self   ! H(x) analysis object
      CLASS (obs_pool), intent(in   ) :: obs    ! observation depot
      integer,          intent(in   ) :: ng     ! nested grid
      integer,          intent(in   ) :: model  ! kernel identifier
!
      integer                         :: nlocs, npool, mvars
      integer                         :: i, iloc, iobs, ivar
      real (r8)                       :: Wmin, Wmax
      real (r8),          allocatable :: work(:)
!
      character (len=*),    parameter :: MyFile =                       &
     &    __FILE__//", roms_hofx_increment"
!
      SourceFile=MyFile
!
!  Process increment at the observation locations.
!
      npool=SIZE(obs%pool)
      DO i=1,npool
        mvars=obs%pool(i)%nvars
        nlocs=obs%pool(i)%nlocs
        allocate ( work(nlocs) )
        DO ivar=1,mvars
          work=0.0_r8
!
!  Read in background.
!
          CALL netcdf_get_fvar (ng, model, obs%pool(i)%ncname_out,      &
     &                          TRIM(obs%pool(i)%vars(ivar)), work,     &
     &                          GrpName = 'hofxInitial')
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Compute increment.
!
          Wmin=spval
          Wmax=-spval
          DO iloc=1,nlocs
            iobs=obs%pool(i)%pack_mapping(iloc,ivar)
            work(iloc)=REAL(self%ObsVetting(iobs), r8)*                 &
     &                 (self%state(iobs)-work(iloc))
            IF (self%ObsVetting(iobs).gt.0) THEN
              Wmin=MIN(Wmin, work(iloc))
              Wmax=MAX(Wmax, work(iloc))
            END IF
          END DO
!
!  Write out increment.
!
          CALL netcdf_put_fvar (ng, model, obs%pool(i)%ncname_out,      &
     &                          TRIM(obs%pool(i)%vars(ivar)), work,     &
     &                          (/1/), (/nlocs/),                       &
     &                          GrpName = 'Increment')
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          IF (Master) THEN
             WRITE (stdout,10) 'Increment',                             &
     &                         TRIM(obs%pool(i)%vars(ivar)),            &
     &                         Wmin, Wmax, TRIM(obs%pool(i)%ncname_out)
          END IF
          deallocate (work)
        END DO
      END DO
!
 10   FORMAT (t22,'Group:',a,'::',a,/,t23,'(Min = ',1p,e12.5,           &
     &        ' Max = ', 1p,e12.5,')   File: ',a)
!
      RETURN
      END SUBROUTINE obs_hofx_increment
!
!-----------------------------------------------------------------------
!  It computes and writes the innovation (observations minus analysis)
!  at the observation locations into multiple output NetCDF-4 files.
!
      SUBROUTINE obs_hofx_innovation (self, obs, ng, model)
!
      CLASS (obs_hofx), intent(in   ) :: self   ! H(x) backgrouns object
      CLASS (obs_pool), intent(in   ) :: obs    ! observation depot
      integer,          intent(in   ) :: ng     ! nested grid
      integer,          intent(in   ) :: model  ! kernel identifier
!
      integer                         :: nlocs, npool, mvars
      integer                         :: i, iloc, iobs, ivar
      real (r8)                       :: Wmin, Wmax
      real (r8),          allocatable :: work(:)
!
      character (len=*),    parameter :: MyFile =                       &
     &    __FILE__//", roms_hofx_innovation"
!
      SourceFile=MyFile
!
!  Process innovation at the observation locations.
!
      npool=SIZE(obs%pool)
      DO i=1,npool
        mvars=obs%pool(i)%nvars
        nlocs=obs%pool(i)%nlocs
        allocate ( work(nlocs) )
        DO ivar=1,mvars
          work=0.0_r8
!
!  Compute innovation.
!
          Wmin=spval
          Wmax=-spval
          DO iloc=1,nlocs
            iobs=obs%pool(i)%pack_mapping(iloc,ivar)
            work(iloc)=REAL(self%ObsVetting(iobs), r8)*                 &
     &                 (obs%pool(i)%value(iloc,ivar)-self%state(iobs))
            IF (self%ObsVetting(iobs).gt.0) THEN
              Wmin=MIN(Wmin, work(iloc))
              Wmax=MAX(Wmax, work(iloc))
            END IF
          END DO
!
!  Write out innovation.
!
          CALL netcdf_put_fvar (ng, model, obs%pool(i)%ncname_out,      &
     &                          TRIM(obs%pool(i)%vars(ivar)), work,     &
     &                          (/1/), (/nlocs/),                       &
     &                          GrpName = 'Innovation')
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          IF (Master) THEN
             WRITE (stdout,10) 'Innovation',                            &
     &                         TRIM(obs%pool(i)%vars(ivar)),            &
     &                         Wmin, Wmax, TRIM(obs%pool(i)%ncname_out)
          END IF
          deallocate (work)
        END DO
      END DO
!
 10   FORMAT (t22,'Group:',a,'::',a,/,t23,'(Min = ',1p,e12.5,           &
     &        ' Max = ', 1p,e12.5,')   File: ',a)
!
      RETURN
      END SUBROUTINE obs_hofx_innovation
!
!-----------------------------------------------------------------------
!  It writes generic H(x) variables into multiple output NetCDF-4 files
!  file(s).
!
      SUBROUTINE obs_hofx_nc_write (self, obs, ng, model, grpname)
!
      CLASS (obs_hofx),  intent(in   ) :: self     ! H(x) object
      CLASS (obs_pool),  intent(in   ) :: obs      ! observation depot
      integer,           intent(in   ) :: ng       ! nested grid
      integer,           intent(in   ) :: model    ! kernel identifier
      character (len=*), intent(in   ) :: grpname  ! NetCDF group name
!
      integer                          :: nlocs, npool, mvars
      integer                          :: i, iloc, iobs, ivar
      integer (i8b),       allocatable :: iwork(:)
      real (r8)                        :: Wmin, Wmax
      real (r8),           allocatable :: work(:)
!
      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", roms_hofx_nc_write"
!
      SourceFile=MyFile
!
!  Write out generic field.
!
      SELECT CASE (TRIM(GrpName))
        CASE ('ObsVetting')
          npool=SIZE(obs%pool)
          DO i=1,npool
            mvars=obs%pool(i)%nvars
            nlocs=obs%pool(i)%nlocs
            allocate ( iwork(nlocs) )
!
            DO ivar=1,mvars
              iwork=0_i8b
              DO iloc=1,nlocs
                iobs=obs%pool(i)%pack_mapping(iloc,ivar)
                iwork(iloc)=self%ObsVetting(iobs)
              END DO
!
              CALL netcdf_put_ivar (ng, model, obs%pool(i)%ncname_out,  &
     &                              obs%pool(i)%vars(ivar), iwork,      &
     &                              (/1/), (/nlocs/),                   &
     &                              GrpName = TRIM(grpname))
              IF (FoundError(exit_flag, NoError,                        &
     &            __LINE__, MyFile)) RETURN
!
              IF (Master) THEN
                WRITE (stdout,10) TRIM(grpname),                        &
     &                            TRIM(obs%pool(i)%vars(ivar)),         &
     &                            MINVAL(iwork), MAXVAL(iwork),         &
     &                            TRIM(obs%pool(i)%ncname_out)
              END IF
              deallocate (iwork)
            END DO
          END DO
!
        CASE DEFAULT
          npool=SIZE(obs%pool)
          DO i=1,npool
            mvars=obs%pool(i)%nvars
            nlocs=obs%pool(i)%nlocs
            allocate ( work(nlocs) )
!
            Wmin=spval
            Wmax=-spval
            DO ivar=1,mvars
              work=0.0_r8
              DO iloc=1,nlocs
                iobs=obs%pool(i)%pack_mapping(iloc,ivar)
                work(iloc)=self%state(iobs)
                IF (self%ObsVetting(iobs).gt.0) THEN
                  Wmin=MIN(Wmin, work(iloc))
                  Wmax=MAX(Wmax, work(iloc))
                END IF
              END DO
!
              CALL netcdf_put_fvar (ng, model, obs%pool(i)%ncname_out,  &
     &                              obs%pool(i)%vars(ivar), work,       &
     &                              (/1/), (/nlocs/),                   &
     &                              GrpName = TRIM(grpname))
              IF (FoundError(exit_flag, NoError,                        &
     &            __LINE__, MyFile)) RETURN
!
              IF (Master) THEN
                WRITE (stdout,20) TRIM(grpname),                        &
     &                            TRIM(obs%pool(i)%vars(ivar)),         &
     &                            Wmin, Wmax,                           &
     &                            TRIM(obs%pool(i)%ncname_out)
              END IF
              deallocate (work)
            END DO
          END DO
      END SELECT
!
 10   FORMAT (t22,'Group:',a,'::',a,/,t23,'(Min = ',i12,              &
     &        ' Max = ',i12,')   File: ',a)
 20   FORMAT (t22,'Group:',a,'::',a,/,t23,'(Min = ',1p,e12.5,         &
     &        ' Max = ', 1p,e12.5,')   File: ',a)
!
      RETURN
      END SUBROUTINE obs_hofx_nc_write
!
!-----------------------------------------------------------------------
!  It computes and writes the residual (observations minus analysis)
!  at the observation locations into multiple output NetCDF-4 files.
!
      SUBROUTINE obs_hofx_residual (self, obs, ng, model)
!
      CLASS (obs_hofx), intent(in   ) :: self   ! H(x) analysis object
      CLASS (obs_pool), intent(in   ) :: obs    ! observation depot
      integer,          intent(in   ) :: ng     ! nested grid
      integer,          intent(in   ) :: model  ! kernel identifier
!
      integer                         :: nlocs, npool, mvars
      integer                         :: i, iloc, iobs, ivar
      real (r8)                       :: Wmin, Wmax
      real (r8),          allocatable :: work(:)
!
      character (len=*),    parameter :: MyFile =                       &
     &    __FILE__//", roms_hofx_residual"
!
      SourceFile=MyFile
!
!  Process residual at the observation locations.
!
      npool=SIZE(obs%pool)
      DO i=1,npool
        mvars=obs%pool(i)%nvars
        nlocs=obs%pool(i)%nlocs
        allocate ( work(nlocs) )
        DO ivar=1,mvars
          work=0.0_r8
!
!  Compute residual.
!
          Wmin=spval
          Wmax=-spval
          DO iloc=1,nlocs
            iobs=obs%pool(i)%pack_mapping(iloc,ivar)
            work(iloc)=REAL(self%ObsVetting(iobs), r8)*                 &
     &                 (obs%pool(i)%value(iloc,ivar)-self%state(iobs))
            IF (self%ObsVetting(iobs).gt.0) THEN
              Wmin=MIN(Wmin, work(iloc))
              Wmax=MAX(Wmax, work(iloc))
            END IF
          END DO
!
!  Write out residual.
!
          CALL netcdf_put_fvar (ng, model, obs%pool(i)%ncname_out,      &
     &                          TRIM(obs%pool(i)%vars(ivar)), work,     &
     &                          (/1/), (/nlocs/),                       &
     &                          GrpName = 'Residual')
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          IF (Master) THEN
             WRITE (stdout,10) 'Residual',                              &
     &                         TRIM(obs%pool(i)%vars(ivar)),            &
     &                         Wmin, Wmax, TRIM(obs%pool(i)%ncname_out)
          END IF
          deallocate (work)
        END DO
      END DO
!
 10   FORMAT (t22,'Group:',a,'::',a,/,t23,'(Min = ',1p,e12.5,           &
     &        ' Max = ', 1p,e12.5,')   File: ',a)
!
      RETURN
      END SUBROUTINE obs_hofx_residual
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
!  It computes several statistics between model and observarions
!  (Oke et al., 2002).
!
      SUBROUTINE obs_hofx_stats (self, obs, ng, model)
!
      CLASS (obs_hofx),   intent(in   ) :: self    ! H(x) object
      CLASS (obs_pool),   intent(in   ) :: obs     ! observation depot
      integer,            intent(in   ) :: ng      ! nested grid
      integer,            intent(in   ) :: model   ! kernel identifier
!
      integer                           :: i, ic, iloc, iobs, ivar
      integer                           :: mvars, nlocs, npool
      integer,   dimension(NobsVar(ng)) :: Ncount
!
      real (r8)                         :: cff1, cff2, mod_val, obs_val
      real (r8), dimension(NobsVar(ng)) :: CC, MB, MSE, SDE
      real (r8), dimension(NobsVar(ng)) :: mod_min, mod_max
      real (r8), dimension(NobsVar(ng)) :: mod_mean, mod_std
      real (r8), dimension(NobsVar(ng)) :: obs_min, obs_max
      real (r8), dimension(NobsVar(ng)) :: obs_mean, obs_std
!
      character (len=11)                :: svar_name(NobsVar(ng))
      character (len=11)                :: text(NobsVar(ng))
!
!  Initialize.
!
      CC=0.0_r8
      MB=0.0_r8
      MSE=0.0_r8
      SDE=0.0_r8
      mod_min=Large
      mod_max=-Large
      mod_mean=0.0_r8
      obs_min=Large
      obs_max=-Large
      obs_mean=0.0_r8
      mod_std=0.0_r8
      obs_std=0.0_r8
      Ncount=0
!
!  Compute statistics.
!
      ic=0
      npool=SIZE(obs%pool)
      DO i=1,npool
        mvars=obs%pool(i)%nvars
        nlocs=obs%pool(i)%nlocs
        DO ivar=1,mvars
          ic=ic+1
          svar_name(ic)=obs%pool(i)%vname(ivar)
          text(ic)='-----------'
!
!  Compute min, max, and mean values.
!
          DO iloc=1,nlocs
            iobs=obs%pool(i)%pack_mapping(iloc,ivar)
            IF (self%ObsVetting(iobs).gt.0) THEN
              mod_val=self%state(iobs)
              obs_val=obs%pool(i)%value(iloc,ivar)
              Ncount(ic)=Ncount(ic)+1
              mod_min(ic)=MIN(mod_min(ic), mod_val)
              obs_min(ic)=MIN(obs_min(ic), obs_val)
              mod_max(ic)=MAX(mod_max(ic), mod_val)
              obs_max(ic)=MAX(obs_max(ic), obs_val)
              mod_mean(ic)=mod_mean(ic)+mod_val
              obs_mean(ic)=obs_mean(ic)+obs_val
            END IF
          END DO
          mod_mean(ic)=mod_mean(ic)/REAL(Ncount(ic),r8)
          obs_mean(ic)=obs_mean(ic)/REAL(Ncount(ic),r8)
!
!  Compute standard deviation and cross-correlation between model and
!  observations (CC).
!
          DO iloc=1,nlocs
            iobs=obs%pool(i)%pack_mapping(iloc,ivar)
            IF (self%ObsVetting(iobs).gt.0) THEN
              mod_val=self%state(iobs)
              obs_val=obs%pool(i)%value(iloc,ivar)
              cff1=mod_val-mod_mean(ic)
              cff2=obs_val-obs_mean(ic)
              mod_std(ic)=mod_std(ic)+cff1*cff1
              obs_std(ic)=obs_std(ic)+cff2*cff2
              CC(ic)=CC(ic)+cff1*cff2
            END IF
          END DO
          mod_std(ic)=SQRT(mod_std(ic)/REAL(Ncount(ic)-1,r8))
          obs_std(ic)=SQRT(obs_std(ic)/REAL(Ncount(ic)-1,r8))
          CC(ic)=(CC(ic)/REAL(Ncount(ic),r8))/(mod_std(ic)*obs_std(ic))
!
!  Compute model bias (MB), standard deviation error (SDE), and mean
!  squared error (MSE).
!
          MB(ic)=mod_mean(ic)-obs_mean(ic)
          SDE(ic)=mod_std(ic)-obs_std(ic)
          MSE(ic)=MB(ic)*MB(ic)+                                        &
     &            SDE(ic)*SDE(ic)+                                      &
     &            2.0_r8*mod_std(ic)*obs_std(ic)*(1.0_r8-CC(ic))
        END DO
      END DO
!
!  Report.
!
      IF (Master) THEN
        WRITE (stdout,10)
        WRITE (stdout,20) (svar_name(i),i=1,ic)
        WRITE (stdout,30) (text(i),i=1,ic)
        WRITE (stdout,40) 'Observation Min   ', (obs_min (i),i=1,ic)
        WRITE (stdout,40) 'Observation Max   ', (obs_max (i),i=1,ic)
        WRITE (stdout,40) 'Observation Mean  ', (obs_mean(i),i=1,ic)
        WRITE (stdout,40) 'Observation STD   ', (obs_std (i),i=1,ic)
        WRITE (stdout,40) 'Model Min         ', (mod_min (i),i=1,ic)
        WRITE (stdout,40) 'Model Max         ', (mod_max (i),i=1,ic)
        WRITE (stdout,40) 'Model Mean        ', (mod_mean(i),i=1,ic)
        WRITE (stdout,40) 'Model STD         ', (mod_std (i),i=1,ic)
        WRITE (stdout,40) 'Model Bias        ', (MB(i),i=1,ic)
        WRITE (stdout,40) 'STD Error         ', (SDE(i),i=1,ic)
        WRITE (stdout,40) 'Cross-Correlation ', (CC(i),i=1,ic)
        WRITE (stdout,40) 'Mean Squared Error', (MSE(i),i=1,ic)
        WRITE (stdout,50) 'Observation Count ', (Ncount(i),i=1,ic)
      END IF
!
 10   FORMAT (/,' Model-Observations Comparison Statistics:',/)
 20   FORMAT (t22,*(a11,1x))
 30   FORMAT (t22,*(a11,1x),/)
 40   FORMAT (a,3x,*(1p,e11.4,0p,1x))
 50   FORMAT (a,3x,*(i11,1x))
!
      RETURN
      END SUBROUTINE obs_hofx_stats
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
