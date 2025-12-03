#include "cppdefs.h"
      MODULE roms_obs_mod

#if (defined FOUR_DVAR || defined VERIFICATION) && defined OBSERVATIONS
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  It contains several routines to process input 4D-Var observation    !
!  file(s).                                                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_netcdf
      USE mod_scalars
!
      USE dateclock_mod,        ONLY : datestr
#ifdef DISTRIBUTE
      USE distribute_mod,       ONLY : mp_assemble, mp_bcasti
#endif
      USE get_hash_mod,         ONLY : get_hash
      USE mod_ncparam,          ONLY : date_str
      USE roms_interpolate_mod, ONLY : hindices
      USE strings_mod,          ONLY : FoundError, nc_error
      USE unique_mod,           ONLY : unique
!
      implicit none
!
!-----------------------------------------------------------------------
!  Enhanced NetCDF-4 group object.
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: nc_group

        integer :: gid                             ! group ID
        integer :: nvars                           ! number of variables
!
        integer, allocatable :: vid(:)             ! variables ID

        character (len=:), allocatable :: name     ! group name
        character (len=:), allocatable :: vname(:) ! variables name

      END TYPE nc_group
!
!-----------------------------------------------------------------------
!  Observation Cost Function Object.
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: obs_cost

        integer :: nvars                       ! number of variables
!
        real (r8), allocatable :: value(:)     ! cost function value

      END TYPE obs_cost
!
!-----------------------------------------------------------------------
!  Observation Data Type Structures/Objects: CLASS(obs_data).
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: obs_data

        integer :: nlocs                       ! number of observations
        integer :: nsurvey                     ! number of times surveys
        integer :: nvars                       ! number of variables
        integer :: poolID                      ! observation pool ID
        integer :: timeWindow                  ! time-averaged windows
!
!  Switches indicating area-averaged and time-averaged observations.
!
        logical :: IsAreaAveraged
        logical :: IsTimeAveraged
!
!  Observation origin identifier. Repetitive observations have their
!  provenance multiplied by -1 so they can be quickly identified during
!  processing.
!
        integer, allocatable :: provenance(:)
!
!  Control vector state variable identifier.
!
        integer, allocatable :: stateID(:)
!
!  Observation survey time indices as they appear in the "time"
!  variable.
!
        integer, allocatable :: surveyIndex(:)
!
!  Starting and ending indices for time-averaging are determined by the
!  specified "TimeStr" and "TimeEnd" range values. A value of zero
!  indicates that there are no time-averaged observations available for
!  that period.
!
        integer, allocatable :: beginIndex(:)    ! time-averaging start
        integer, allocatable :: finishIndex(:)   ! time-averagind end
!
!  Starting and ending timestep and times for time-averaged operator.
!
        integer :: TimeStepStr                   ! starting timestep
        integer :: TimeStepEnd                   ! ending timestep
!
        real (r8), allocatable :: TimeStr(:)     ! starting time
        real (r8), allocatable :: TimeEnd(:)     ! ending time
!
!  Starting and ending observation vector indices available for
!  requested model time interval, time-0.5*dt to time+0.5*dt. They
!  are zero, if no observations are available in the time window.
!
        integer :: NstrObs
        integer :: NendObs
!
!  Mapping indices for the packed observations vector. It is used
!  when computing H(x), the cost function, and minimization where
!  and the observations are clustered in a 1D array.
!
        integer, allocatable :: pack_mapping(:,:)
!
!  Half-length spatial scale (m) for area-averaged operator.
!
        real (r8) :: spatialAverage
!
!  Observation spatial and time coordinates. The fractional values are
!  used for efficient interpolation in the H(x) operator.
!
        real (r8), allocatable :: lon(:)        ! observation longitude
        real (r8), allocatable :: lat(:)        ! observation latitude
        real (r8), allocatable :: depth(:)      ! depth of observation
        real (r8), allocatable :: time(:)       ! time  of observation
        real (r8), allocatable :: surveyTime(:) ! survey time
!
        real (r8), allocatable :: Xgrid(:)      ! fractional X-grid
        real (r8), allocatable :: Ygrid(:)      ! fractional Y-grid
        real (r8), allocatable :: Zgrid(:)      ! fractional Z-grid
!
!  If the observations are scalar, the second dimension is unity.
!  Otherwise, it is used for vectors or multiple observations at the
!  exact spatial and time coordinates.
!
        real (r8), allocatable :: error(:,:)    ! observation error(s)
        real (r8), allocatable :: value(:,:)    ! observation value(s)
!
        character (len=:), allocatable :: vars(:) ! NetCDF variable name
!
!  Observation statistics.
!
        integer (i8b), allocatable :: error_hash(:)
        integer (i8b), allocatable :: value_hash(:)
!
        real (r8),     allocatable :: error_avg(:)
        real (r8),     allocatable :: error_min(:)
        real (r8),     allocatable :: error_max(:)

        real (r8),     allocatable :: value_avg(:)
        real (r8),     allocatable :: value_min(:)
        real (r8),     allocatable :: value_max(:)
!
        real (r8) :: time_min
        real (r8) :: time_max
!
        character (len=22) :: date_min
        character (len=22) :: date_max
!
!  Internal observation data object name.
!
        character (len=:), allocatable :: name
        character (len=:), allocatable :: vname(:)
!
!  Input observation filename. It supports native ROMS 4D-Var NetCDF
!  files and IODA enhanced NetCDF-4 files with Groups.
!
        character (len=:), allocatable :: ncname_inp
!
!  Output observation plus H(x) operator enhanced NetCDF-4 file
!  information.
!
        integer :: ncid                                ! NetCDF file ID
        integer :: ngroups                             ! groups size
!
        TYPE (nc_group), allocatable :: group(:)       ! group object
!
        character (len=:),  allocatable :: ncname_out  ! output filename
!
        CONTAINS
!
        PROCEDURE :: destroy             => obs_data_destroy
        PROCEDURE :: frac_coords         => obs_data_frac_coords
        PROCEDURE :: averaging_timesteps => obs_data_averaging_timesteps
        PROCEDURE :: populate_ioda       => obs_data_populate_ioda
        PROCEDURE :: qc_threshold        => obs_data_qc_threshold
        PROCEDURE :: stats               => obs_data_stats
!
      END TYPE obs_data
!
!-----------------------------------------------------------------------
!  Collection Observation Data Types Object: CLASS(obs_pool).
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: obs_pool
!
!  Total number of observations in pool depot.
!
        integer :: ndatum
!
!  Observation ID associated with control state variable, [1:ntype].
!  (zeta=1, ubar=2, vbar=3, u=4, v=5, temperature=6, salinity=7)
!
        integer, allocatable :: stateID(:)
!
!  Observation set object. Each element corresponds to an observation
!  associated with a state variable, [1:ntype].
!
        TYPE (obs_data), allocatable :: pool(:)
!
        CONTAINS
!
        PROCEDURE :: create           => obs_pool_create
        PROCEDURE :: destroy          => obs_pool_destroy
        PROCEDURE :: get              => obs_pool_get
        PROCEDURE :: has              => obs_pool_has
        PROCEDURE :: load_ioda        => obs_pool_load_ioda
        PROCEDURE :: load_native      => obs_pool_load_native
        PROCEDURE :: nc_create        => obs_pool_nc_create
        PROCEDURE :: nc_read          => obs_pool_nc_read
        PROCEDURE :: nc_write         => obs_pool_nc_write
        PROCEDURE :: populate_native  => obs_pool_populate_native
        PROCEDURE :: size             => obs_pool_size
        PROCEDURE :: tally            => obs_pool_tally
        PROCEDURE :: unpack_i         => obs_pool_unpack_i
        PROCEDURE :: unpack_r         => obs_pool_unpack_r
!
      END TYPE obs_pool
!
!-----------------------------------------------------------------------
!  Module variables.
!-----------------------------------------------------------------------
!
!  Define observation storage object. Currently, there is unecessary to
!  have different in nested grid applications since fractional grid
!  coordinates can be recomputed from they geographical locations.
!
      TYPE (obs_pool) :: obs_store
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      PUBLIC  :: obs_initialize
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      CONTAINS
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!  It allocates the observation set object. Then, reads observation
!  file(s), and populates data into the depot object. It should be
!  called once during ROMS initialization phase.
!
      SUBROUTINE obs_initialize (ng, model, nfiles)
!
      integer,          intent(in   ) :: ng       ! nested grid
      integer,          intent(in   ) :: model    ! kernel identifier
      integer,          intent(in   ) :: nfiles   ! number of files
!
      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_initialize"
!
      SourceFile=MyFile
!
!  Creates and populates observation pool object.
!
      CALL obs_pool_create (obs_store, ng, model, nfiles)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      RETURN
      END SUBROUTINE obs_initialize
!
!<><><><><><><><><><><><><><><><><><><><><><><><>  CLASS OBS_DATA <><><>
!
!  It destroys a observation data object.
!
      SUBROUTINE obs_data_destroy (self)
!
      CLASS (obs_data), intent(inout) :: self     ! observation object

!  Deallocate observation data object variables.
!
      IF (allocated(self%pack_mapping)) deallocate (self%pack_mapping)
      IF (allocated(self%provenance))   deallocate (self%provenance)
      IF (allocated(self%lon))          deallocate (self%lon)
      IF (allocated(self%lat))          deallocate (self%lat)
      IF (allocated(self%depth))        deallocate (self%depth)
      IF (allocated(self%time))         deallocate (self%time)
      IF (allocated(self%stateID))      deallocate (self%stateID)
      IF (allocated(self%surveyIndex))  deallocate (self%surveyIndex)
      IF (allocated(self%surveyTime))   deallocate (self%surveyTime)
      IF (allocated(self%vars))         deallocate (self%vars)
      IF (allocated(self%Xgrid))        deallocate (self%Xgrid)
      IF (allocated(self%Ygrid))        deallocate (self%Ygrid)
      IF (allocated(self%Zgrid))        deallocate (self%Zgrid)
!
      IF (allocated(self%error))        deallocate (self%error)
      IF (allocated(self%value))        deallocate (self%value)
!
      IF (allocated(self%error_avg))    deallocate (self%error_avg)
      IF (allocated(self%error_hash))   deallocate (self%error_hash)
      IF (allocated(self%error_min))    deallocate (self%error_min)
      IF (allocated(self%error_max))    deallocate (self%error_max)

      IF (allocated(self%value_avg))    deallocate (self%value_avg)
      IF (allocated(self%value_hash))   deallocate (self%value_hash)
      IF (allocated(self%value_min))    deallocate (self%value_min)
      IF (allocated(self%value_max))    deallocate (self%value_max)
!
      IF (allocated(self%name))         deallocate (self%name)
      IF (allocated(self%vname))        deallocate (self%vname)
      IF (allocated(self%ncname_inp))   deallocate (self%ncname_inp)
      IF (allocated(self%ncname_out))   deallocate (self%ncname_out)
!
      RETURN
      END SUBROUTINE obs_data_destroy
!
!-----------------------------------------------------------------------
!  If not available, it computes observation object fractional
!  coordinates.
!
      SUBROUTINE obs_data_frac_coords (self, ng, model, tile)
!
      CLASS (obs_data), intent(inout) :: self     ! observation object
      integer,          intent(in   ) :: ng       ! nested grid
      integer,          intent(in   ) :: model    ! kernel identifier
      integer,          intent(in   ) :: tile     ! domain partition
!
      logical                         :: rectangular
      integer                         :: LBi, UBi, LBj, UBj
      integer                         :: IstrR, IendR, JstrR, JendR
!
!  Compute observation fractional coordinates.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
      IstrR=BOUNDS(ng)%IstrR(tile)
      IendR=BOUNDS(ng)%IendR(tile)
      JstrR=BOUNDS(ng)%JstrR(tile)
      JendR=BOUNDS(ng)%JendR(tile)
!
      rectangular=.FALSE.
!
      self%Xgrid=0.0_r8
      self%Ygrid=0.0_r8
!
      CALL hindices (ng, LBi, UBi, LBj, UBj,                            &
     &               IstrR, IendR, JstrR, JendR,                        &
     &               GRID(ng)%angler, GRID(ng)%lonr, GRID(ng)%latr,     &
     &               1, self%nlocs, 1, 1,                               &
     &               1, self%nlocs, 1, 1,                               &
     &               self%lon, self%lat,                                &
     &               self%Xgrid, self%Ygrid,                            &
     &               spval, rectangular)

#ifdef DISTRIBUTE
!
!  Collect fractional coordinates.
!
      CALL mp_assemble (ng, model, self%nlocs, spval, self%Xgrid)
      CALL mp_assemble (ng, model, self%nlocs, spval, self%Ygrid)
#endif
!
      RETURN
      END SUBROUTINE obs_data_frac_coords
!
!-----------------------------------------------------------------------
!  If time-averaged observations, return current averaging window
!  timesteps.
!
      SUBROUTINE obs_data_averaging_timesteps (self, ng, model,         &
     &                                         iTstr, iTend)
!
      CLASS (obs_data), intent(inout) :: self     ! observation object
      integer,          intent(in   ) :: ng       ! nested grid
      integer,          intent(in   ) :: model    ! kernel identifier
      integer,          intent(out  ) :: iTstr    ! starting timestep
      integer,          intent(out  ) :: iTend    ! ending   timestep
!
      integer                         :: i
!
!  Given current model time, get averaging window starting and ending
!  times.
!
      iTstr=-1
      iTend=-1
!
      IF (self%IsTimeAveraged) THEN
        DO i=1,self%timeWindow
          IF (model.eq.iADM) THEN
            IF ((self%TimeStr(i).le.time(ng)).and.                      &
     &          (time(ng).lt.self%TimeEnd(i))) THEN
              iTstr=INT((self%TimeStr(i)-dstart*day2sec)/dt(ng))+1
              iTend=INT((self%TimeEnd(i)-dstart*day2sec)/dt(ng))
            END IF
          ELSE
            IF ((self%TimeStr(i).le.time(ng)).and.                      &
     &          (time(ng).lt.self%TimeEnd(i))) THEN
              iTstr=INT((self%TimeStr(i)-dstart*day2sec)/dt(ng))+1
              iTend=INT((self%TimeEnd(i)-dstart*day2sec)/dt(ng))
              EXIT
            END IF
          END IF
        END DO
      END IF
!
      RETURN
      END SUBROUTINE obs_data_averaging_timesteps
!
!-----------------------------------------------------------------------
!  If not available, it computes observation object fractional
!  coordinates.
!
      SUBROUTINE obs_data_populate_ioda (self, ng, model, mvar, IOtype)
!
      CLASS (obs_data),  intent(inout) :: self    ! observation object
      integer,           intent(in   ) :: ng      ! nested grid
      integer,           intent(in   ) :: model   ! kernel identifier
      integer,           intent(in   ) :: mvar    ! number of variables
      integer,           intent(in   ) :: IOtype  ! NetCDF IO type
!
      logical                          :: Loaded, isDayUnits
      integer,             allocatable :: wcount(:)
      integer                          :: i, ic, is, j, lstr, ncount
      real (r8),          dimension(4) :: Tmin, Tmax
      character (len=22), dimension(4) :: Dmin, Dmax

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_data_populate_ioda"
!
      SourceFile=MyFile
!
!  Read in observation variables in Group MetaData.
!
      IF (Master) WRITE (stdout,10) TRIM(self%ncname_inp)
!
      IF (IOtype.eq.io_nf90) THEN
        CALL netcdf_inq_var (ng, model, self%ncname_inp,                &
                             GrpName = 'MetaData')
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        DO i=1,n_var
          Loaded=.FALSE.
          SELECT CASE (TRIM(var_name(i)))
            CASE ('dateTime')
              CALL netcdf_get_time (ng, model, self%ncname_inp,         &
     &                              var_name(i), Rclock%DateNumber,     &
     &                              self%time,                          &
     &                              GrpName = 'MetaData',               &
     &                              min_val = Tmin(1),                  &
     &                              max_val = Tmax(1))
              Loaded=.TRUE.
              ncount=SIZE(self%time)
              CALL datestr (Tmin(1), .FALSE., Dmin(1))
              CALL datestr (Tmax(1), .FALSE., Dmax(1))
            CASE ('depth')
              CALL netcdf_get_fvar (ng, model, self%ncname_inp,         &
     &                              var_name(i), self%depth,            &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%depth)
            CASE ('latitude')
              CALL netcdf_get_fvar (ng, model, self%ncname_inp,         &
     &                              var_name(i), self%lat,              &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%lat)
            CASE ('longitude')
              CALL netcdf_get_fvar (ng, model, self%ncname_inp,         &
     &                              var_name(i), self%lon,              &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%lon)
            CASE ('provenance')
              CALL netcdf_get_ivar (ng, model, self%ncname_inp,         &
     &                              var_name(i), self%provenance,       &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%provenance)
            CASE ('stateID')
              CALL netcdf_get_ivar (ng, model, self%ncname_inp,         &
     &                              var_name(i), self%stateID,          &
     &                              GrpName = 'MetaData')
            CASE ('surveyIndex')
              CALL netcdf_get_ivar (ng, model, self%ncname_inp,         &
     &                              var_name(i), self%surveyIndex,      &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%surveyIndex)
            CASE ('surveyTime')
              CALL netcdf_get_time (ng, model, self%ncname_inp,         &
     &                              var_name(i), Rclock%DateNumber,     &
     &                              self%surveyTime,                    &
     &                              GrpName = 'MetaData',               &
     &                              min_val = Tmin(2),                  &
     &                              max_val = Tmax(2))
              Loaded=.TRUE.
              ncount=SIZE(self%surveyTime)
              CALL datestr (Tmin(2), .FALSE., Dmin(2))
              CALL datestr (Tmax(2), .FALSE., Dmax(2))
            CASE ('dateTimeAverageBegin')
              CALL netcdf_get_time (ng, model, self%ncname_inp,         &
     &                              var_name(i), Rclock%DateNumber,     &
     &                              self%TimeStr,                       &
     &                              GrpName = 'MetaData',               &
     &                              min_val = Tmin(3),                  &
     &                              max_val = Tmax(3))
              Loaded=.TRUE.
              ncount=SIZE(self%TimeStr)
              CALL datestr (Tmin(3), .FALSE., Dmin(3))
              CALL datestr (Tmax(3), .FALSE., Dmax(3))
            CASE ('dateTimeAverageEnd')
              CALL netcdf_get_time (ng, model, self%ncname_inp,         &
     &                              var_name(i), Rclock%DateNumber,     &
     &                              self%TimeEnd,                       &
     &                              GrpName = 'MetaData',               &
     &                              min_val = Tmin(4),                  &
     &                              max_val = Tmax(4))
              Loaded=.TRUE.
              ncount=SIZE(self%TimeEnd)
              CALL datestr (Tmin(4), .FALSE., Dmin(4))
              CALL datestr (Tmax(4), .FALSE., Dmax(4))
            CASE ('spatialAverage')
              CALL netcdf_get_fvar (ng, model, self%ncname_inp,         &
     &                              var_name(i), self%spatialAverage,   &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              self%IsAreaAveraged=.TRUE.
              ncount=1
            CASE ('x_grid')
              CALL netcdf_get_fvar (ng, model, self%ncname_inp,         &
     &                              var_name(i), self%Xgrid,            &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%Xgrid)
            CASE ('y_grid')
              CALL netcdf_get_fvar (ng, model, self%ncname_inp,         &
     &                              var_name(i), self%Ygrid,            &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%Ygrid)
            CASE ('z_grid')
              CALL netcdf_get_fvar (ng, model, self%ncname_inp,         &
     &                              var_name(i), self%Zgrid,            &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%Zgrid)
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IF (Master.and.Loaded) WRITE (stdout,20) 'MetaData',          &
     &                                 TRIM(var_name(i)), ncount
        END DO
# if defined PIO_LIB && defined DISTRIBUTE
      ELSE IF (IOtype.eq.io_pio) THEN
!
!  HGA: Is possible to read Group variables with PIO?
!
# endif
      END IF
!
!  Read the observation variables in Group ObsError. If the IODA file is
!  created with the "roms2ioda.m" script, the observation error is in
!  terms of standard deviations and has the same units as its value.
!
      IF (IOtype.eq.io_nf90) THEN
        CALL netcdf_inq_var (ng, model, self%ncname_inp,                &
                             GrpName = 'ObsError')
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        DO i=1,n_var
          CALL netcdf_get_fvar (ng, model, self%ncname_inp,             &
     &                          var_name(i), self%error(:,i),           &
     &                          GrpName = 'ObsError')
          self%error(:,i)=self%error(:,i)*self%error(:,i)
          ncount=SIZE(self%error,1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IF (Master) WRITE (stdout,30) 'ObsError', TRIM(var_name(i)),  &
     &                                  ncount, self%stateID(i)
        END DO
# if defined PIO_LIB && defined DISTRIBUTE
      ELSE IF (IOtype.eq.io_pio) THEN
!
!  HGA: Is possible to read Group variables with PIO?
!
# endif
      END IF
!
!  Read in observation variables in Group ObsValue.
!
      IF (IOtype.eq.io_nf90) THEN
        CALL netcdf_inq_var (ng, model, self%ncname_inp,                &
                             GrpName = 'ObsValue')
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        lstr=0
        DO i=1,n_var
          CALL netcdf_get_fvar (ng, model, self%ncname_inp,             &
     &                          var_name(i), self%value(:,i),           &
     &                          GrpName = 'ObsValue')
          lstr=MAX(lstr, LEN_TRIM(var_name(i)))
          ncount=SIZE(self%value,1)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          IF (Master) WRITE (stdout,30) 'ObsValue', TRIM(var_name(i)),  &
     &                                  ncount, self%stateID(i)
        END DO
# if defined PIO_LIB && defined DISTRIBUTE
      ELSE IF (IOtype.eq.io_pio) THEN
!
!  HGA: Is possible to read Group variables with PIO?
!
# endif
      END IF
!
!  Load observation NetCDF variable name(s).
!
      IF (.not.allocated(self%vars)) THEN
        allocate ( character(LEN=lstr) :: self%vars(mvar) )
      END IF
      DO i=1,n_var
        self%vars(i)=TRIM(var_name(i))
      END DO
!
!  If time-averaged observations count the number of observations for
!  each averaging time window and the starting indices as they appear
!  in the observation time coordinate.
!
      IF (self%IsTimeAveraged) THEN
        IF (.not.allocated(wcount)) THEN
          allocate ( wcount(self%timeWindow) )
          wcount=0
        END IF
!
        ic=0
        is=1
        DO j=1,self%timeWindow
          wcount(j)=ic
          DO i=is,self%nlocs
            IF ((self%TimeStr(j).le.self%time(i)).and.                  &
     &          (self%time(i).lt.self%TimeEnd(j))) THEN
              wcount(j)=wcount(j)+1
              IF (i.eq.is) THEN
                self%beginIndex(j)=i
              END IF
              self%finishIndex(j)=wcount(j)
            ELSE
              is=i
              ic=is-1
              EXIT
            END IF
          END DO
        END DO
      END IF
!
!  Compute statistics.
!
      isDayUnits=.FALSE.
      CALL self%stats (isDayUnits)
!
!  Report information.
!
      IF (Master) THEN
        WRITE (stdout,'(1x)')
        DO i=1,mvar
          WRITE (stdout,40) TRIM(self%vars(i)),                         &
     &                      self%value_min(i), self%value_max(i),       &
     &                      self%value_avg(i), self%value_hash(i),      &
     &                      self%error_min(i), self%error_max(i),       &
     &                      self%error_avg(i)
          WRITE (stdout,50) 'DateTime Range:   ',                       &
     &                      TRIM(Dmin(1)), TRIM(Dmax(1))
          WRITE (stdout,50) 'SurveyTime Range: ',                       &
     &                      TRIM(Dmin(2)), TRIM(Dmax(2))
        END DO
      END IF
!
!  Deallocate local variable.
!
      IF (allocated(wcount)) deallocate (wcount)
!
  10  FORMAT (/,2x,'ROMS_OBS         - Processing: ',a)
  20  FORMAT (t23,'Group: ',a,'::',a,t75,'Count = ',i0)
  30  FORMAT (t23,'Group: ',a,'::',a,t75,'Count = ',i0,t90,'ID = ',i0)
  40  FORMAT (t21,'- ',a,/,t23,'(ValMin = ',1p,e12.5,' ValMax = ',      &
     &        1p,e12.5,' ValAvg = ',1p,e12.5,' CheckSum = ',i0,')',     &
     &        /,t23,'(ErrMin = ',1p,e12.5,' ErrMax = ',1p,e12.5,        &
     &        ' ErrAvg = ',1p,e12.5,')')
  50  FORMAT (t23,a,a,' to ',a)
!
      RETURN
      END SUBROUTINE obs_data_populate_ioda
!
!-----------------------------------------------------------------------
!  It sets the background quality control threshold in terms of the
!  squared number of standard deviations to use for accepting or
!  rejecting observations according specified input parameters.
!
      SUBROUTINE obs_data_qc_threshold (self, ng, model, ivar, Bthresh)
!
      CLASS (obs_data), intent(in   ) :: self       ! observation object
      integer,          intent(in   ) :: ng         ! nested grid
      integer,          intent(in   ) :: model      ! kernel identifier
      integer,          intent(in   ) :: ivar       ! variable counter
      real (r8),        intent(out  ) :: Bthresh(:) ! threshold
!
      integer                         :: i, iloc
!
!  Background quality control in terms of state variable indices.
!
      IF (bgqc_type(ng).eq.1) THEN
        DO i=1,MstateVar
          IF (self%stateID(ivar).eq.i) THEN
            Bthresh=S_bgqc(i,ng)
            EXIT
          END IF
        END DO
!
!  Background quality control in term of observation provenance.
!
      ELSE IF (bgqc_type(ng).eq.2) THEN
        DO iloc=1,self%nlocs
          Bthresh(iloc)=bgqc_large            ! do not reject
          DO i=1,Nprovenance(ng)
            IF (self%provenance(iloc).eq.Iprovenance(i,ng)) THEN
              Bthresh(iloc)=P_bgqc(i,ng)
              EXIT
            END IF
          END DO
        END DO
      END IF
!
      RETURN
      END SUBROUTINE obs_data_qc_threshold
!
!-----------------------------------------------------------------------
!  It computes and reports observation object statistics.
!
      SUBROUTINE obs_data_stats (self, isDayUnits)
!
      CLASS (obs_data), intent(inout) :: self       ! observation object
      logical,          intent(in   ) :: isDayUnits ! days or seconds
!
      integer (i8b)                   :: checksum
      integer                         :: i
!
      character (len=*),    parameter :: MyFile =                       &
     &    __FILE__//", obs_data_stats"
!
      SourceFile=MyFile
!
!  Compute observation statistics.
!
      self%time_min=MINVAL(self%time)
      self%time_max=MAXVAL(self%time)
!
      IF (isDayUnits) THEN
        CALL datestr (self%time_min, .TRUE., self%date_min)
        CALL datestr (self%time_max, .TRUE., self%date_max)
      ELSE
        CALL datestr (self%time_min, .FALSE., self%date_min)
        CALL datestr (self%time_max, .FALSE., self%date_max)
      END IF
!
      DO i=1,self%nvars
        self%error_avg(i)=SUM(self%error(:,i))/REAL(self%nlocs, r8)
        self%error_min(i)=MINVAL(self%error(:,i))
        self%error_max(i)=MAXVAL(self%error(:,i))
        CALL get_hash (self%error(:,i), self%nlocs, checksum, .FALSE.)
        self%error_hash(i)=checksum
!
        self%value_avg(i)=SUM(self%value(:,i))/REAL(self%nlocs, r8)
        self%value_min(i)=MINVAL(self%value(:,i))
        self%value_max(i)=MAXVAL(self%value(:,i))
        CALL get_hash (self%value(:,i), self%nlocs, checksum, .FALSE.)
        self%value_hash(i)=checksum
      END DO
!
      RETURN
      END SUBROUTINE obs_data_stats
!
!<><><><><><><><><><><><><><><><><><><><><><><><>  CLASS OBS_POOL <><><>
!
!  It creates a set observation-type objects. A single or multiple
!  data assimilation observation files are allowed. In IODA, we have
!  a file for each observation type.
!
      SUBROUTINE obs_pool_create (self, ng, model, nfiles)
!
      CLASS (obs_pool), intent(inout) :: self     ! observation set
      integer,          intent(in   ) :: ng       ! nested grid
      integer,          intent(in   ) :: model    ! kernel identifier
      integer,          intent(in   ) :: nfiles   ! number of files
!
# ifdef DISTRIBUTE
      integer,           dimension(2) :: ibuffer
# endif
      character (len=*),    parameter :: MyFile =                       &
     &    __FILE__//", obs_pool_create"
!
      SourceFile=MyFile
!
!  Inquire observation file(s) to determine the size of the observation
!  pool vector and allocate.
!
      CALL self%size (ng, model, nfiles)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Read in observation files, allocate arrays in pool structure, and
!  populate data into object.
!
      CALL self%nc_read (ng, model, nfiles)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Count total number of observations.
!
      CALL self%tally ()
!
!  If input multiple IODA observation files, create output IODA files
!  containing processed observations plus H(x) operators. The files are
!  created serially by the output thread.
!
      IF (nfiles.gt.1) THEN
        CALL self%nc_create (ng, model)
# ifdef DISTRIBUTE
!
        ibuffer(1)=exit_flag
        ibuffer(2)=ioerror
        CALL mp_bcasti (ng, model, ibuffer)
        exit_flag=ibuffer(1)
        ioerror=ibuffer(2)
# endif
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
! HGA: For now test writing of output IODA files here.
!
      IF (nfiles.gt.1) THEN
        CALL self%nc_write (ng, model)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END IF
!
      RETURN
      END SUBROUTINE obs_pool_create
!
!-----------------------------------------------------------------------
!  It destroys a set observation-type objects.
!
      SUBROUTINE obs_pool_destroy (self)
!
      CLASS (obs_pool), intent(inout) :: self     ! observation set
!
      integer :: i
!
!  Deallocate observation collection object.
!
      DO i=1,SIZE(self%pool)
        CALL self%pool(i)%destroy ()
      END DO
      deallocate (self%pool)
!
      RETURN
      END SUBROUTINE obs_pool_destroy
!
!-----------------------------------------------------------------------
!  Given internal observation name, it returns an observation data
!  object pointer.
!
      SUBROUTINE obs_pool_get (self, name, field)
!
      CLASS (obs_pool), target, intent(in ) :: self   ! observation set
      character (len=*),        intent(in ) :: name   ! observation name
      TYPE (obs_data), pointer, intent(out) :: field  ! data pointer
!
      integer                               :: iobs

      character (len=*),          parameter :: MyFile =                 &
     &    __FILE__//", obs_pool_get"
!
      SourceFile=MyFile
!
!  Get requested field pointer.
!
      DO iobs=1,SIZE(self%pool)
        IF (TRIM(name).eq.self%pool(iobs)%name) THEN
          field => self%pool(iobs)
          RETURN
        END IF
      END DO
!
!  Error observation data was not found.
!
      IF (Master) THEN
        WRITE (stdout,10) TRIM(name)
  10    FORMAT (/, 'ROMS_OBS - obs_pool::get: cannot find field ',a)
      END IF
      exit_flag=5
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      END SUBROUTINE obs_pool_get
!
!-----------------------------------------------------------------------
!  It inquires about available observation for current model time, and
!  their starting and ending indices on pool set.
!
      FUNCTION obs_pool_has (self, ng, model) RESULT (foundit)
!
      CLASS (obs_pool), intent(inout) :: self     ! observation set
      integer,          intent(in   ) :: ng       ! nested grid
      integer,          intent(in   ) :: model    ! observation name
!
      logical                         :: foundit
      integer                         :: Nstr, Nend, iobs, i
      integer                         :: nsurvey, nwindow
      real (r8)                       :: Tstr, Tend
!
!  Check if requested observation name exists in the pool object.
!
      foundit=.FALSE.
      Tstr=time(ng)-0.5_r8*dt(ng)
      Tend=time(ng)+0.5_r8*dt(ng)
!
      DO iobs=1,SIZE(self%pool)
        nsurvey=self%pool(iobs)%nsurvey
        nwindow=self%pool(iobs)%timeWindow
        self%pool(iobs)%NstrObs=0
        self%pool(iobs)%NendObs=0
!
        IF (model.eq.iADM) THEN                   ! backward search
          IF (self%pool(iobs)%IsTimeAveraged) THEN
             DO i=nwindow,1,-1
              IF ((self%pool(iobs)%TimeStr(i).le.time(ng)).and.         &
     &            (time(ng).lt.self%pool(iobs)%TimeEnd(i))) THEN
                IF ((self%pool(iobs)%beginIndex(i).gt.0).and.           &
     &              (self%pool(iobs)%finishIndex(i).gt.0)) THEN
                  Nstr=self%pool(iobs)%beginIndex(i)
                  Nend=self%pool(iobs)%finishIndex(i)
                  self%pool(iobs)%NstrObs=Nstr
                  self%pool(iobs)%NendObs=Nend
                  foundit=.TRUE.
                  CYCLE
                END IF
              END IF
            END DO
          ELSE
            DO i=nsurvey,1,-1
              IF ((Tstr.le.self%pool(iobs)%surveyTime(i)).and.          &
     &            (self%pool(iobs)%surveyTime(i).lt.Tend)) THEN
                Nstr=self%pool(iobs)%surveyIndex(i)
                IF (i.ne.1) THEN
                  Nend=self%pool(iobs)%surveyIndex(i-1)-1
                ELSE
                  Nend=self%pool(iobs)%surveyIndex(i)
                END IF
                self%pool(iobs)%NstrObs=Nstr
                self%pool(iobs)%NendObs=Nend
                foundit=.TRUE.
                CYCLE
              END IF
            END DO
          END IF
        ELSE                                      ! forward search
          IF (self%pool(iobs)%IsTimeAveraged) THEN
            DO i=1,nwindow
              IF ((self%pool(iobs)%TimeStr(i).le.time(ng)).and.         &
     &            (time(ng).lt.self%pool(iobs)%TimeEnd(i))) THEN
                IF ((self%pool(iobs)%beginIndex(i).gt.0).and.           &
     &              (self%pool(iobs)%finishIndex(i).gt.0)) THEN
                  Nstr=self%pool(iobs)%beginIndex(i)
                  Nend=self%pool(iobs)%finishIndex(i)
                  self%pool(iobs)%NstrObs=Nstr
                  self%pool(iobs)%NendObs=Nend
                  foundit=.TRUE.
                  CYCLE
                END IF
              END IF
            END DO
          ELSE
            DO i=1,nsurvey
              IF ((Tstr.le.self%pool(iobs)%surveyTime(i)).and.          &
     &            (self%pool(iobs)%surveyTime(i).lt.Tend)) THEN
                Nstr=self%pool(iobs)%surveyIndex(i)
                IF (i.ne.nsurvey) THEN
                  Nend=self%pool(iobs)%surveyIndex(i+1)-1
                ELSE
                  Nend=self%pool(iobs)%surveyIndex(i)
                END IF
                self%pool(iobs)%NstrObs=Nstr
                self%pool(iobs)%NendObs=Nend
                foundit=.TRUE.
                CYCLE
              END IF
            END DO
          END IF
        END IF
      END DO
!
      RETURN
      END FUNCTION obs_pool_has
!
!-----------------------------------------------------------------------
!  It allocates the observation arrays in the appropiate pool vector
!  element, and then reads the specified IODA observation file to
!  populate the data.
!
      SUBROUTINE obs_pool_load_ioda (self, ng, model, mobs, msurvey,    &
     &                               mvar, IOtype, ncname)
!
      CLASS (obs_pool),  intent(inout) :: self     ! observation set
      integer,           intent(in   ) :: ng       ! nested grid
      integer,           intent(in   ) :: model    ! kernel identifier
      integer,           intent(in   ) :: mobs     ! observations number
      integer,           intent(in   ) :: msurvey  ! time surveys
      integer,           intent(in   ) :: mvar     ! number of variables
      integer,           intent(in   ) :: IOtype   ! NetCDF IO type
      character (len=*), intent(in   ) :: ncname   ! NetCDF filename
!
      logical                          :: haveSST, haveSSS, haveVector
      integer,         dimension(mvar) :: myStateID
      integer                          :: Idot, Islash
      integer                          :: i, iobs, lstr, myID
      integer                          :: mtimes, mvatt
!
      character (len=40)               :: poolname, varname(mvar)
      character (len=80)               :: string

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_pool_load_ioda"
!
      SourceFile=MyFile
!
!  Get the associated state variable ID.
!
      haveSST=.FALSE.
      haveSSS=.FALSE.
!
      SELECT CASE (IOtype)
        CASE (io_nf90)
          CALL netcdf_get_ivar (ng, model, TRIM(ncname),                &
     &                          'stateID', myStateID,                   &
     &                          GrpName = 'MetaData')
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL netcdf_inq_var (ng, model, TRIM(ncname),                 &
     &                         myVarName = 'stateID',                   &
     &                         GrpName = 'MetaData',                    &
     &                         nVarAtt = mvatt)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!

          DO i=1,mvatt
            IF (TRIM(var_Aname(i)).eq.'surface_level') THEN
              IF (ANY(myStateID.eq.6)) haveSST=.TRUE.
              IF (ANY(myStateID.eq.7)) haveSSS=.TRUE.
            END IF
          END DO

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL pio_netcdf_get_ivar (ng, model, TRIM(ncname),            &
     &                              'stateID', myStateID,               &
     &                              GrpName = 'MetaData')
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) IOtype
          exit_flag=3
      END SELECT
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Determine pool vector element to process, U- and V-velocity.
!
      haveVector=.FALSE.
      IF (mvar.gt.1) THEN
        IF (ANY(myStateID.eq.4).and.ANY(myStateID.eq.5)) THEN
          haveVector=.TRUE.
          DO i=1,SIZE(self%stateID)
            IF (self%stateID(i).eq.45) THEN
              iobs=i
              myID=45
              EXIT
            END IF
          END DO
        END IF
      ELSE
        DO i=1,SIZE(self%stateID)
          IF ((haveSST.and.(self%stateID(i).eq.60)).or.                 &
     &        (haveSSS.and.(self%stateID(i).eq.70))) THEN
            iobs=i
            myID=self%stateID(i)
            EXIT
          ELSE IF (myStateID(1).eq.self%stateID(i)) THEN
            iobs=i
            myID=self%stateID(i)
            EXIT
          END IF
        END DO
      END IF
!
!  Initialize scalar variables in observation pool object.
!
      self%pool(iobs)%nlocs   = mobs
      self%pool(iobs)%nsurvey = msurvey
      self%pool(iobs)%nvars   = mvar
      self%pool(iobs)%poolID  = myID
!
      IF (.not.allocated(self%pool(iobs)%ncname_inp)) THEN
        lstr=LEN_TRIM(ncname)
        allocate ( character(LEN=lstr) :: self%pool(iobs)%ncname_inp )
        self%pool(iobs)%ncname_inp = TRIM(ncname)
      END IF
!
      IF (.not.allocated(self%pool(iobs)%ncname_out)) THEN
        Idot  =INDEX(ncname,CHAR(46),BACK=.TRUE.)-1
        Islash=INDEX(ncname,CHAR(47),BACK=.TRUE.)+1
        string=ncname(Islash:Idot)//'_hofx.nc4'
        lstr=LEN_TRIM(string)
        allocate ( character(LEN=lstr) :: self%pool(iobs)%ncname_out )
        self%pool(iobs)%ncname_out = TRIM(string)
      END IF
!
!  Allocate pool vector element structure, Group MetaData.
!
      IF (.not.allocated(self%pool(iobs)%provenance)) THEN
        allocate ( self%pool(iobs)%provenance(mobs) )
        self%pool(iobs)%provenance = 0
        Dmem(ng)=Dmem(ng)+REAL(mobs,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%lon)) THEN
        allocate ( self%pool(iobs)%lon(mobs) )
        self%pool(iobs)%lon = 0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mobs,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%lat)) THEN
        allocate ( self%pool(iobs)%lat(mobs) )
        self%pool(iobs)%lat = 0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mobs,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%depth)) THEN
        allocate ( self%pool(iobs)%depth(mobs) )
        self%pool(iobs)%depth = 0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mobs,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%time)) THEN
        allocate ( self%pool(iobs)%time(mobs) )
        self%pool(iobs)%time = 0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mobs,r8)
      END IF

      IF (self%pool(iobs)%timeWindow.gt.0) THEN
        mtimes=self%pool(iobs)%timeWindow
        IF (.not.allocated(self%pool(iobs)%TimeStr)) THEN
          allocate ( self%pool(iobs)%TimeStr(mtimes) )
          self%pool(iobs)%TimeStr=0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mtimes,r8)
        END IF
        IF (.not.allocated(self%pool(iobs)%TimeEnd)) THEN
          allocate ( self%pool(iobs)%TimeEnd(mtimes) )
          self%pool(iobs)%TimeEnd=0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mtimes,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%beginIndex)) THEN
          allocate ( self%pool(iobs)%beginIndex(mtimes) )
          self%pool(iobs)%beginIndex=0
          Dmem(ng)=Dmem(ng)+REAL(mtimes,r8)
        END IF
        IF (.not.allocated(self%pool(iobs)%finishIndex)) THEN
          allocate ( self%pool(iobs)%finishIndex(mtimes) )
          self%pool(iobs)%finishIndex=0
          Dmem(ng)=Dmem(ng)+REAL(mtimes,r8)
        END IF
      END IF

      IF (.not.allocated(self%pool(iobs)%stateID)) THEN
        allocate ( self%pool(iobs)%stateID(mvar) )
        self%pool(iobs)%stateID = 0
        Dmem(ng)=Dmem(ng)+REAL(mvar,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%surveyIndex)) THEN
        allocate ( self%pool(iobs)%surveyIndex(msurvey) )
        self%pool(iobs)%surveyIndex = 0
        Dmem(ng)=Dmem(ng)+REAL(msurvey,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%surveyTime)) THEN
        allocate ( self%pool(iobs)%surveyTime(msurvey) )
        self%pool(iobs)%surveyTime = 0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(msurvey,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%Xgrid)) THEN
        allocate ( self%pool(iobs)%Xgrid(mobs) )
        self%pool(iobs)%Xgrid = 0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mobs,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%Ygrid)) THEN
        allocate ( self%pool(iobs)%Ygrid(mobs) )
        self%pool(iobs)%Ygrid = 0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mobs,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%Zgrid)) THEN
        allocate ( self%pool(iobs)%Zgrid(mobs) )
        self%pool(iobs)%Zgrid = 0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mobs,r8)
      END IF
!
!  Allocate pool vector element structure, Group ObsError.
!
      IF (.not.allocated(self%pool(iobs)%error)) THEN
        allocate ( self%pool(iobs)%error(mobs,mvar) )
        self%pool(iobs)%error = 0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mobs*mvar,r8)
      END IF
!
!  Allocate pool vector element structure, Group ObsValue.
!
      IF (.not.allocated(self%pool(iobs)%value)) THEN
        allocate ( self%pool(iobs)%value(mobs,mvar) )
        self%pool(iobs)%value = 0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mobs*mvar,r8)
      END IF
!
!  Allocate statistics variables.
!
      IF (.not.allocated(self%pool(iobs)%error_hash)) THEN
        allocate ( self%pool(iobs)%error_hash(mvar) )
        self%pool(iobs)%error_hash=0_i8b
        Dmem(ng)=Dmem(ng)+REAL(mvar,r8)
      END IF
      IF (.not.allocated(self%pool(iobs)%value_hash)) THEN
        allocate ( self%pool(iobs)%value_hash(mvar) )
        self%pool(iobs)%value_hash=0_i8b
        Dmem(ng)=Dmem(ng)+REAL(mvar,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%error_avg)) THEN
        allocate ( self%pool(iobs)%error_avg(mvar) )
        self%pool(iobs)%error_avg=0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mvar,r8)
      END IF
      IF (.not.allocated(self%pool(iobs)%value_avg)) THEN
        allocate ( self%pool(iobs)%value_avg(mvar) )
        self%pool(iobs)%value_avg=0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mvar,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%error_min)) THEN
        allocate ( self%pool(iobs)%error_min(mvar) )
        self%pool(iobs)%error_min=0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mvar,r8)
      END IF
      IF (.not.allocated(self%pool(iobs)%value_min)) THEN
        allocate ( self%pool(iobs)%value_min(mvar) )
        self%pool(iobs)%value_min=0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mvar,r8)
      END IF

      IF (.not.allocated(self%pool(iobs)%error_max)) THEN
        allocate ( self%pool(iobs)%error_max(mvar) )
        self%pool(iobs)%error_max=0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mvar,r8)
      END IF
      IF (.not.allocated(self%pool(iobs)%value_max)) THEN
        allocate ( self%pool(iobs)%value_max(mvar) )
        self%pool(iobs)%value_max=0.0_r8
        Dmem(ng)=Dmem(ng)+REAL(mvar,r8)
      END IF
!
!  Allocate internal variable name(s).
!
      SELECT CASE (self%stateID(iobs))
        CASE (1)
          poolname='SSH'
          varname(1)='SSH'
        CASE (6)
          poolname='temp'
          varname(1)='temp'
        CASE (7)
          poolname='salt'
          varname(1)='salt'
        CASE (45)
          poolname='uv'
          varname(1)='Ucodar'
          Varname(2)='Vcodar'
        CASE (60)
          poolname='SST'
          varname(1)='SST'
        CASE (67)
          poolname='TS'
          varname(1)='temp'
          varname(2)='salt'
        CASE (70)
          poolname='SSS'
          varname(1)='SSS'
      END SELECT
!
      lstr=LEN_TRIM(poolname)
      allocate ( character(LEN=lstr) :: self%pool(iobs)%name )
      self%pool(iobs)%name=TRIM(poolname)
!
      lstr=0
      DO i=1,mvar
        lstr=MAX(lstr, LEN_TRIM(varname(i)))
      END DO
      allocate ( character(LEN=lstr) :: self%pool(iobs)%vname(mvar) )
      DO i=1,mvar
        self%pool(iobs)%vname(i)=TRIM(varname(i))
      END DO
!
!  Read in IODA NetCDF file and populate observation data.
!
      CALL self%pool(iobs)%populate_ioda (ng, model, mvar, IOtype)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' OBS_POOL_LOAD_IODA - Illegal output file type,',        &
     &        ' io_type = ',i0,/,22x,                                   &
     &        'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE obs_pool_load_ioda
!
!-----------------------------------------------------------------------
!  It allocates the observation arrays in the appropiate pool vector
!  element, and then reads the specified native observation file to
!  populate the data.
!
      SUBROUTINE obs_pool_load_native (self, ng, model, mobs, msurvey,  &
     &                                 IOtype, ncname)
!
      CLASS (obs_pool),  intent(inout) :: self     ! observation set
      integer,           intent(in   ) :: ng       ! nested grid
      integer,           intent(in   ) :: model    ! kernel identifier
      integer,           intent(in   ) :: mobs     ! observations number
      integer,           intent(in   ) :: msurvey  ! time surveys
      integer,           intent(in   ) :: IOtype   ! NetCDF IO type
      character (len=*), intent(in   ) :: ncname   ! NetCDF filename
!
      integer,     dimension(5+NT(ng)) :: mvar
      integer,         dimension(mobs) :: obs_type
      integer                          :: i, id, iobs, lstr, mSSS, mSST
      integer                          :: Lsur, mymvar, mysize
      real (r8),       dimension(mobs) :: obs_Zgrid

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_pool_load_native"
!
      SourceFile=MyFile
!
!  Read in "obs_type" variable to determine the observed variables and
!  its IDs.
!
      SELECT CASE (IOtype)
        CASE (io_nf90)
          CALL netcdf_get_ivar (ng, model, ncname,                      &
     &                          'obs_type', obs_type)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL netcdf_get_fvar (ng, model, ncname,                      &
     &                          'obs_Zgrid', obs_Zgrid)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

# if defined PIO_LIB && defined DISTRIBUTE
        CASE (io_pio)
          CALL pio_netcdf_get_ivar (ng, model, ncname,                  &
     &                              'obs_type', obs_type)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
          CALL pio_netcdf_get_fvar (ng, model, ncname,                  &
     &                              'obs_Zgrid', obs_type)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
# endif
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) IOtype
          exit_flag=3
      END SELECT
!
!  Determine the number of observations per control vector variable.
!
!     Lsur=N(ng)
      Lsur=50                    ! HGA: temporary for debugging
      mSSS=0
      mSST=0
      mvar=0
!
      DO i=1,mobs
        id=obs_type(i)
        IF (id.eq.6) THEN
          IF (INT(obs_Zgrid(i)).eq.Lsur) THEN
            mSST=mSST+1
          ELSE
            mvar(id)=mvar(id)+1
          END IF
        ELSE
          mvar(id)=mvar(id)+1
        END IF
      END DO
!
!  Allocate observation pool vector object arrays and initialize
!  scalars.
!
      DO iobs=1,SIZE(self%stateID)
        id=self%stateID(iobs)
        IF (id.eq.45) THEN              ! U-, V-Velocity
          mysize=mvar(4)
          mymvar=2
        ELSE IF (id.eq.60) THEN         ! SST
          mysize=mSST
          mymvar=1
        ELSE IF (id.eq.70) THEN         ! SSS
          mysize=mSSS
          mymvar=1
        ELSE
          mysize=mvar(id)               ! SSH, T, S
          mymvar=1
        END IF
!
        self%pool(iobs)%nlocs   = mysize
        self%pool(iobs)%nvars   = mymvar
        self%pool(iobs)%poolID  = self%stateID(iobs)
!
        IF (.not.allocated(self%pool(iobs)%ncname_inp)) THEN
          lstr=LEN_TRIM(ncname)
          allocate ( character(LEN=lstr) :: self%pool(iobs)%ncname_inp )
          self%pool(iobs)%ncname_inp = TRIM(ncname)
        END IF
!
        IF (.not.allocated(self%pool(iobs)%provenance)) THEN
          allocate ( self%pool(iobs)%provenance(mysize) )
          self%pool(iobs)%provenance = 0
          Dmem(ng)=Dmem(ng)+REAL(mysize,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%lon)) THEN
          allocate ( self%pool(iobs)%lon(mysize) )
          self%pool(iobs)%lon = 0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mysize,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%lat)) THEN
          allocate ( self%pool(iobs)%lat(mysize) )
          self%pool(iobs)%lat = 0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mysize,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%depth)) THEN
          allocate ( self%pool(iobs)%depth(mysize) )
          self%pool(iobs)%depth = 0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mysize,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%time)) THEN
          allocate ( self%pool(iobs)%time(mysize) )
          self%pool(iobs)%time = 0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mysize,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%stateID)) THEN
          allocate ( self%pool(iobs)%stateID(mymvar) )
          self%pool(iobs)%stateID = 0
          Dmem(ng)=Dmem(ng)+REAL(mymvar,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%Xgrid)) THEN
          allocate ( self%pool(iobs)%Xgrid(mysize) )
          self%pool(iobs)%Xgrid = 0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mysize,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%Ygrid)) THEN
          allocate ( self%pool(iobs)%Ygrid(mysize) )
          self%pool(iobs)%Ygrid = 0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mysize,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%Zgrid)) THEN
          allocate ( self%pool(iobs)%Zgrid(mysize) )
          self%pool(iobs)%Zgrid = 0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mysize,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%error)) THEN
          allocate ( self%pool(iobs)%error(mysize,mymvar) )
          self%pool(iobs)%error = 0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mysize*mymvar,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%value)) THEN
          allocate ( self%pool(iobs)%value(mysize,mymvar) )
          self%pool(iobs)%value = 0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mysize*mymvar,r8)
        END IF
!
!  Allocate statistics variables.
!
        IF (.not.allocated(self%pool(iobs)%error_hash)) THEN
          allocate ( self%pool(iobs)%error_hash(mymvar) )
          self%pool(iobs)%error_hash=0_i8b
          Dmem(ng)=Dmem(ng)+REAL(mymvar,r8)
        END IF
        IF (.not.allocated(self%pool(iobs)%value_hash)) THEN
          allocate ( self%pool(iobs)%value_hash(mymvar) )
          self%pool(iobs)%value_hash=0_i8b
          Dmem(ng)=Dmem(ng)+REAL(mymvar,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%error_avg)) THEN
          allocate ( self%pool(iobs)%error_avg(mymvar) )
          self%pool(iobs)%error_avg=0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mymvar,r8)
        END IF
        IF (.not.allocated(self%pool(iobs)%value_avg)) THEN
          allocate ( self%pool(iobs)%value_avg(mymvar) )
          self%pool(iobs)%value_avg=0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mymvar,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%error_min)) THEN
          allocate ( self%pool(iobs)%error_min(mymvar) )
          self%pool(iobs)%error_min=0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mymvar,r8)
        END IF
        IF (.not.allocated(self%pool(iobs)%value_min)) THEN
          allocate ( self%pool(iobs)%value_min(mymvar) )
          self%pool(iobs)%value_min=0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mymvar,r8)
        END IF

        IF (.not.allocated(self%pool(iobs)%error_max)) THEN
          allocate ( self%pool(iobs)%error_max(mymvar) )
          self%pool(iobs)%error_max=0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mymvar,r8)
        END IF
        IF (.not.allocated(self%pool(iobs)%value_max)) THEN
          allocate ( self%pool(iobs)%value_max(mymvar) )
          self%pool(iobs)%value_max=0.0_r8
          Dmem(ng)=Dmem(ng)+REAL(mymvar,r8)
        END IF

      END DO
!
!  Read in native NetCDF file and populate obervation data.
!
      CALL self%populate_native (ng, model, mobs, IOTYPE, ncname,       &
     &                           obs_type, obs_Zgrid)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
  10  FORMAT (' OBS_POOL_LOAD_NATIVE - Illegal output file type,',      &
     &        ' io_type = ',i0,/,24x,                                   &
     &        'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE obs_pool_load_native
!
!-----------------------------------------------------------------------
!  It creates ouput IODA-type observations plus H(x) variables NetCDF4
!  enhanced file(s). It inquires input IODA-type observation files to
!  copy their schema, and adds Groups associated with the H(x) operator.
!
      SUBROUTINE obs_pool_nc_create (self, ng, model)
!
      CLASS (obs_pool),  intent(inout) :: self     ! observation set
      integer,           intent(in   ) :: ng       ! nested grid
      integer,           intent(in   ) :: model    ! kernel identifier
!
      integer                          :: ncid_inp, ncid_out, npool
      integer                          :: dimid, grpid, varid, vtype
      integer                          :: indx, gindex, qcindex
      integer                          :: i, iobs, j, k, kc
      integer                          :: lstr, vlstr
      integer                          :: ngroups, nvars, nvdim
      integer,            dimension(5) :: vdims
!
      character (len=80)               :: string

# ifdef VERIFICATION
      character (len=15), dimension(2) :: new_group =                   &
     &                                    (/ 'ObsVetting     ',         &
     &                                       'hofx           ' /)
# else
      character (len=15), dimension(7) :: new_group =                   &
     &                                    (/ 'ObsVetting     ',         &
     &                                       'BackgroundError',         &
     &                                       'hofxInitial    ',         &
     &                                       'hofxFinal      ',         &
     &                                       'Innovation     ',         &
     &                                       'Increment      ',         &
     &                                       'Residual       ' /)
#endif
      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_pool_nc_create"
!
      SourceFile=MyFile
!
!  Loop over all available files.
!
      npool=SIZE(self%pool)
      IF (Master) WRITE (stdout,'(1x)')
!
      FILE_LOOP : DO iobs=1,npool
!
!  Open input observation file.
!
        CALL netcdf_open (ng, model, self%pool(iobs)%ncname_inp, 0,     &
     &                    ncid_inp)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about input file dimensions and format.
!
        CALL netcdf_get_dim (ng, model, self%pool(iobs)%ncname_inp,     &
     &                       ncid = ncid_inp)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Inquire about input file groups.
!
        IF (ncformat.eq.nf90_format_netcdf4) THEN
          CALL netcdf_get_grp (ng, model, self%pool(iobs)%ncname_inp,   &
     &                         ncid = ncid_inp)
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          IF (Master) WRITE (stdout,10) ncformat,                       &
     &                                  self%pool(iobs)%ncname_inp
          exit_flag=5
          RETURN
        END IF
!
!  Inquire about the variables.
!
        CALL netcdf_inq_var (ng, model, self%pool(iobs)%ncname_inp,     &
     &                       ncid = ncid_inp)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!-----------------------------------------------------------------------
!  Copy NetCDF4 file schema (serial, only output thread).
!-----------------------------------------------------------------------
!
        DEFINE : IF ((ncformat.eq.nf90_format_netcdf4).and.             &
     &               OutThread) THEN
!
!  Create NetCDF4 file.
!
          IF (nc_error(nf90_create(self%pool(iobs)%ncname_out,          &
     &                             nf90_netcdf4, ncid_out),             &
     &                 nf90_noerr, __LINE__, MyFile)) RETURN
          self%pool(iobs)%ncID=ncid_out
          WRITE (stdout,20) ng, self%pool(iobs)%ncname_out
!
!  Define dimensions.
!
          DO i=1,n_dim
            IF (nc_error(nf90_def_dim(ncid_out, TRIM(dim_name(i)),      &
     &                                dim_size(i), dimid),              &
     &                   nf90_noerr, __LINE__, MyFile)) RETURN
          END DO
!
!  Allocate ouput NetCDF group structure in pool object.
!
          ngroups=n_grp+SIZE(new_group)
          self%pool(iobs)%ngroups=ngroups
          IF (.not.allocated(self%pool(iobs)%group)) THEN
            allocate ( self%pool(iobs)%group(ngroups) )
          END IF
!
!  Define variables.
!
          IF (n_grp.gt.0) THEN
            DO k=1,n_grp
              nvars=n_gvar(k)
              IF (nc_error(nf90_def_grp(ncid_out, TRIM(grp_name(k)),    &
     &                                  grpid),                         &
     &                     nf90_noerr, __LINE__, MyFile)) RETURN
!
              self%pool(iobs)%group(k)%gid=grpid
              IF (.not.allocated(self%pool(iobs)%group(k)%name)) THEN
                lstr=LEN_TRIM(grp_name(k))
                allocate ( character(LEN=lstr) ::                       &
     &                     self%pool(iobs)%group(k)%name )
                self%pool(iobs)%group(k)%name=TRIM(grp_name(k))
              END IF
              self%pool(iobs)%group(k)%nvars=nvars
!
              IF (TRIM(grp_name(k)).eq.'ObsValue') THEN
                gindex=k     ! used to duplicate variables in new groups
              ELSE IF (TRIM(grp_name(k)).eq.'PreQC') THEN
                qcindex=k    ! used to duplicate interger flags
              END IF
!
              vlstr=0
              DO j=1,nvars
                vlstr=MAX(vlstr, LEN_TRIM(grp_vname(j,k)))
              END DO
!
              vdims=0
              DO j=1,nvars
                nvdim=grp_nvdim(j,k)
                vdims(1:nvdim)=grp_vdim(1:nvdim,j,k)
                vtype=grp_vtype(j,k)
                IF (nc_error(nf90_def_var(grpid, TRIM(grp_vname(j,k)),  &
     &                                    vtype, vdims(1:nvdim),        &
     &                                    varid),                       &
     &                       nf90_noerr, __LINE__, MyFile)) RETURN
!
                IF (.not.allocated(self%pool(iobs)%group(k)%vid)) THEN
                  allocate ( self%pool(iobs)%group(k)%vid(nvars) )
                  self%pool(iobs)%group(k)%vid=0
                END IF
                self%pool(iobs)%group(k)%vid(j)=varid
!
                IF (.not.allocated(self%pool(iobs)%group(k)%vname)) THEN
                  allocate ( character(LEN=vlstr) ::                    &
     &                       self%pool(iobs)%group(k)%vname(nvars) )
                END IF
                self%pool(iobs)%group(k)%vname(j)=TRIM(grp_vname(j,k))
!
                DO i=1,grp_nvatt(j,k)
                  IF (nc_error(nf90_copy_att(grp_id(k), grp_vid(j,k),   &
     &                                       TRIM(grp_aname(i,j,k)),    &
     &                                       grpid, varid),             &
     &                         nf90_noerr, __LINE__, MyFile)) RETURN
                END DO
              END DO
            END DO
          ELSE
            WRITE (stdout,30) n_grp, self%pool(iobs)%ncname_out
            exit_flag=5
            RETURN
          END IF
!
!  Add new groups: hofx_initial, hofx_final, innovation, increment, and
!  residual.
!
          DO k=1,SIZE(new_group)
            IF (nc_error(nf90_def_grp(ncid_out, TRIM(new_group(k)),     &
     &                                grpid),                           &
     &                   nf90_noerr, __LINE__, MyFile)) RETURN
!
            IF (k.eq.1) THEN
              indx=qcindex         ! integer variables like PreQC
            ELSE
              indx=gindex          ! float variables like ObsValue
            END IF
            nvars=n_gvar(indx)
!
            kc=n_grp+k
            self%pool(iobs)%group(kc)%gid=grpid
            IF (.not.allocated(self%pool(iobs)%group(kc)%name)) THEN
              lstr=LEN_TRIM(new_group(k))
              allocate ( character(LEN=lstr) ::                         &
     &                   self%pool(iobs)%group(kc)%name )
              self%pool(iobs)%group(kc)%name=TRIM(new_group(k))
            END IF
            self%pool(iobs)%group(kc)%nvars=n_gvar(indx)
!
            vlstr=0
            DO j=1,nvars
              vlstr=MAX(vlstr, LEN_TRIM(grp_vname(j,indx)))
            END DO
!
            IF (.not.allocated(self%pool(iobs)%group(kc)%vid)) THEN
              allocate ( self%pool(iobs)%group(kc)%vid(nvars) )
              self%pool(iobs)%group(kc)%vid=0
            END IF
            IF (.not.allocated(self%pool(iobs)%group(kc)%vname)) THEN
              allocate ( character(LEN=vlstr) ::                        &
     &                   self%pool(iobs)%group(kc)%vname(nvars) )
            END IF
!
            vdims=0
            DO j=1,nvars
              nvdim=grp_nvdim(j,indx)
              vdims(1:nvdim)=grp_vdim(1:nvdim,j,indx)
              vtype=grp_vtype(j,indx)
              IF (nc_error(nf90_def_var(grpid, TRIM(grp_vname(j,indx)), &
     &                                  vtype, vdims(1:nvdim), varid),  &
     &                     nf90_noerr, __LINE__, MyFile)) RETURN
!
              self%pool(iobs)%group(kc)%vid(j)=varid
              self%pool(iobs)%group(kc)%vname(j)=TRIM(grp_vname(j,indx))
!
              DO i=1,grp_nvatt(j,indx)
                IF (TRIM(grp_aname(i,j,indx)).eq.'long_name') THEN
                  SELECT CASE (k)
                    CASE (1)
                      string='Observation screening flag, '//           &
     &                       '0:reject and 1:accept'
                    CASE (2)
                      string='background error at observation locations'
                    CASE (3)
                      string='model background at observation locations'
                    CASE (4)
                      string='model analysis at observation locations'
                    CASE (5)
                      string='observation minus background'
                    CASE (6)
                      string='analysis minus background'
                    CASE (7)
                      string='observation minus analysis'
                  END SELECT
                  IF (nc_error(nf90_put_att(grpid, varid, 'long_name',  &
     &                                      TRIM(string)),              &
     &                         nf90_noerr, __LINE__, MyFile)) RETURN
                ELSE
                  IF (nc_error(nf90_copy_att(grp_id(indx),              &
     &                                       grp_vid(j,indx),           &
     &                                       TRIM(grp_aname(i,j,indx)), &
     &                                       grpid, varid),             &
     &                         nf90_noerr, __LINE__, MyFile)) RETURN
                END IF
              END DO
            END DO
          END DO
!
!  Copy global attributes.
!
          DO i=1,n_gatt
            IF (nc_error(nf90_inq_attname(ncid_inp, nf90_global, i,     &
     &                                    string),                      &
     &                   nf90_noerr, __LINE__, MyFile)) RETURN
            IF (TRIM(string).eq.'history') THEN
              WRITE (history,'(a,1x,a,", ",a)') 'ROMS, Version',        &
     &                                          TRIM(version),          &
     &                                          TRIM(date_str)
              IF (nc_error(nf90_put_att(ncid_out, nf90_global,          &
     &                                  TRIM(string), TRIM(history)),   &
     &                     nf90_noerr, __LINE__, MyFile)) RETURN
            ELSE
              IF (nc_error(nf90_copy_att(ncid_inp, nf90_global,         &
     &                                   TRIM(string),                  &
     &                                   ncid_out, nf90_global),        &
     &                     nf90_noerr, __LINE__, MyFile)) RETURN
            END IF
          END DO
!
!  End Definition.
!
          IF (nc_error(nf90_enddef(ncid_out),                           &
     &                 nf90_noerr, __LINE__, MyFile)) RETURN
!
!  Close output to syncronize to disk and allow parallel access.
!
          IF (nc_error(nf90_close(ncid_out),                             &
     &                 nf90_noerr, __LINE__, MyFile)) RETURN
        END IF DEFINE
!
!  Close input file.
!
        CALL netcdf_close (ng, model, ncid_inp,                         &
     &                     self%pool(iobs)%ncname_inp, .FALSE.)

# ifdef DISTRIBUTE
!
!  In distribute-memory applications, only the master process have the
!  ouput file information and allocated the "nc_group" structure. Thus,
!  inquire created output file and allocate values for other processes.
!
        CALL netcdf_get_grp (ng, model, self%pool(iobs)%ncname_out)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        IF (MyRank.ne.MyMaster) THEN
          IF (.not.allocated(self%pool(iobs)%group)) THEN
            allocate ( self%pool(iobs)%group(n_grp) )
          END IF
          DO k=1,n_grp
            nvars=n_gvar(k)
            IF (.not.allocated(self%pool(iobs)%group(k)%name)) THEN
              lstr=LEN_TRIM(grp_name(k))
              allocate ( character(LEN=lstr) ::                         &
     &                   self%pool(iobs)%group(k)%name )
              self%pool(iobs)%group(k)%name=TRIM(grp_name(k))
            END IF
            self%pool(iobs)%group(k)%nvars=nvars
!
            vlstr=0
            DO j=1,nvars
              vlstr=MAX(vlstr, LEN_TRIM(grp_vname(j,k)))
            END DO
!
            IF (.not.allocated(self%pool(iobs)%group(k)%vid)) THEN
              allocate ( self%pool(iobs)%group(k)%vid(nvars) )
              self%pool(iobs)%group(k)%vid=0
            END IF
!
            IF (.not.allocated(self%pool(iobs)%group(k)%vname)) THEN
              allocate ( character(LEN=vlstr) ::                        &
     &                   self%pool(iobs)%group(k)%vname(nvars) )
            END IF
!
            DO j=1,nvars
              self%pool(iobs)%group(k)%vid(j)=grp_vid(j,k)
              self%pool(iobs)%group(k)%vname(j)=TRIM(grp_vname(j,k))
            END DO
          END DO
        END IF
# endif
      END DO FILE_LOOP
!
  10  FORMAT (/,'OBS_POOL_NC_CREATE - non-IODA-type observation file,', &
     &        ' Format = ',i0,/,t21,'File: ',a)
  20  FORMAT (2x,'OBS_POOL         - NC_CREATE: creating H(x) file,',   &
     &        t56,'Grid ',i2.2,': ',a)
  30  FORMAT (/,'OBS_POOL_NC_CREATE - Cannot find Groups in ',          &
     &        'observation NetCDF file, n_grp = ',i0,/,t21,a)
!
      RETURN
      END SUBROUTINE obs_pool_nc_create
!
!-----------------------------------------------------------------------
!  It reads the input observation file(s) to populate the data into
!  the pool object. Each element of the pool vector stores a scalar
!  or vector variable(s) associated with the data assimilation control
!  vector.
!
      SUBROUTINE obs_pool_nc_read (self, ng, model, nfiles)
!
      CLASS (obs_pool), intent(inout) :: self     ! observation set
      integer,          intent(in   ) :: ng       ! nested grid
      integer,          intent(in   ) :: model    ! kernel identifier
      integer,          intent(in   ) :: nfiles   ! number of files
!
      logical                         :: isNative
      integer                         :: mobs, msurvey, mvar
      integer                         :: IOtype, i, nf
      character (len=256)             :: ncname

      character (len=*),    parameter :: MyFile =                       &
     &    __FILE__//", obs_pool_nc_read"
!
      SourceFile=MyFile
!
!  Initialize.
!
      isNative=.FALSE.
      mobs=0
      msurvey=0
      mvar=0
!
!  Inquire input observation NetCDF file(s).
!
      FILES_LOOP : DO nf=1,nfiles
        ncname=OBS2(nf,ng)%name
        IOtype=OBS2(nf,ng)%IOtype
!
        self%pool(nf)%timeWindow=0
        self%pool(nf)%spatialAverage=0.0_r8
        self%pool(nf)%IsAreaAveraged=.FALSE.
        self%pool(nf)%IsTimeAveraged=.FALSE.
!
        SELECT CASE (IOtype)
          CASE (io_nf90)
            CALL netcdf_get_dim (ng, model, ncname)
# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_get_dim (ng, model, ncname)
# endif
          CASE DEFAULT
            IF (Master) WRITE (stdout,10) IOtype
            exit_flag=3
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Get number of observations from its dimensions, and determine if
!  native or IODA file.
!
        DO i=1,n_dim
          IF (TRIM(dim_name(i)).eq.'datum') then
            mobs=dim_size(i)
            isNative=.TRUE.                             ! only in native
          ELSE IF (TRIM(dim_name(i)).eq.'survey') THEN
            msurvey=dim_size(i)                         ! in both
          ELSE IF (TRIM(dim_name(i)).eq.'Location') THEN
            mobs=dim_size(i)
            isNative=.FALSE.                            ! only in IODA
          ELSE IF (TRIM(dim_name(i)).eq.'nvars') THEN
            mvar=dim_size(i)                            ! only in IODA
          ELSE IF (TRIM(dim_name(i)).eq.'timeWindow') THEN
            self%pool(nf)%timeWindow=dim_size(i)        ! only in IODA
            self%pool(nf)%IsTimeAveraged=.TRUE.
          END IF
        END DO
!
!  Populate data into observation pool object.
!
        IF (isNative) THEN
          CALL self%load_native (ng, model, mobs, msurvey,              &
     &                           IOtype, TRIM(ncname))
        ELSE
          CALL self%load_ioda (ng, model, mobs, msurvey, mvar,          &
     &                         IOtype, TRIM(ncname))
        END IF
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END DO FILES_LOOP
!
  10  FORMAT (' OBS_POOL_READ - Illegal output file type, io_type = ',  &
     &        i0,/,17x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE obs_pool_nc_read
!
!-----------------------------------------------------------------------
!  It writes observation data ouput IODA-type files. The H(x) operator
!  variables will be written elsewhere.
!
      SUBROUTINE obs_pool_nc_write (self, ng, model)
!
      CLASS (obs_pool),  intent(inout) :: self     ! observation set
      integer,           intent(in   ) :: ng       ! nested grid
      integer,           intent(in   ) :: model    ! kernel identifier
!
      integer,             allocatable :: Iwork(:)
      integer                          :: i, j, k, ncid, npool
      integer                          :: nlocs, nsurvey, nvars, nwindow
      character (len=40)               :: grpname, varname
      character (len=256)              :: ncname

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_pool_nc_write"
!
      SourceFile=MyFile
!
!  Loop over all available files.
!
      npool=SIZE(self%pool)
!
      FILE_LOOP : DO i=1,npool
!
!  Open output observation plus H(x) file in the pool depot.
!
        ncname=TRIM(self%pool(i)%ncname_out)
        CALL netcdf_open (ng, model, TRIM(ncname), 1, ncid)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        self%pool(i)%ncid=ncid
!
!  Write out observation variables.
!
        DO k=1,self%pool(i)%ngroups
          grpname=TRIM(self%pool(i)%group(k)%name)
          nlocs=self%pool(i)%nlocs
          nsurvey=self%pool(i)%nsurvey
          nvars=self%pool(i)%nvars
          nwindow=self%pool(i)%timeWindow
          DO j=1,self%pool(i)%group(k)%nvars
            varname=TRIM(self%pool(i)%group(k)%vname(j))
            IF (TRIM(grpname).eq.'MetaData') THEN
              SELECT CASE (TRIM(varname))
                CASE ('dateTime')
                  IF (.not.allocated(Iwork)) THEN
                    allocate ( Iwork(nlocs) )
                  END IF
                  CALL netcdf_get_ivar (ng, model,                      &
     &                                  self%pool(i)%ncname_inp,        &
     &                                  TRIM(varname), Iwork,           &
     &                                  GrpName = TRIM(grpname))
                  IF (FoundError(exit_flag, NoError,                    &
     &                           __LINE__, MyFile)) RETURN
!
                  CALL netcdf_put_ivar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  Iwork,                          &
     &                                  (/1/), (/nlocs/),               &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                  deallocate (Iwork)
                CASE ('dateTimeAverageBegin',                           &
     &                'dateTimeAverageEnd')
                  IF (.not.allocated(Iwork)) THEN
                    allocate ( Iwork(nwindow) )
                  END IF
                  CALL netcdf_get_ivar (ng, model,                      &
     &                                  self%pool(i)%ncname_inp,        &
     &                                  TRIM(varname), Iwork,           &
     &                                  GrpName = TRIM(grpname))
                  IF (FoundError(exit_flag, NoError,                    &
     &                           __LINE__, MyFile)) RETURN
!
                  CALL netcdf_put_ivar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  Iwork,                          &
     &                                  (/1/), (/nwindow/),             &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                  deallocate (Iwork)
                CASE ('depth')
                  CALL netcdf_put_fvar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  self%pool(i)%depth,             &
     &                                  (/1/), (/nlocs/),               &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                CASE ('latitude')
                  CALL netcdf_put_fvar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  self%pool(i)%lat,               &
     &                                  (/1/), (/nlocs/),               &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                CASE ('longitude')
                  CALL netcdf_put_fvar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  self%pool(i)%lon,               &
     &                                  (/1/), (/nlocs/),               &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                CASE ('provenance')
                  CALL netcdf_put_ivar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  self%pool(i)%provenance,        &
     &                                  (/1/), (/nlocs/),               &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                CASE ('sequenceNumber')
                  IF (.not.allocated(Iwork)) THEN
                    allocate ( Iwork(nlocs) )
                  END IF
                  CALL netcdf_get_ivar (ng, model,                      &
     &                                  self%pool(i)%ncname_inp,        &
     &                                  TRIM(varname), Iwork,           &
     &                                  GrpName = TRIM(grpname))
                  IF (FoundError(exit_flag, NoError,                    &
     &                           __LINE__, MyFile)) RETURN
!
                  CALL netcdf_put_ivar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  Iwork,                          &
     &                                  (/1/), (/nlocs/),               &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                  deallocate (Iwork)
                CASE ('stateID')
                  CALL netcdf_put_ivar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  self%pool(i)%stateID,           &
     &                                  (/1/), (/nvars/),               &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                CASE ('spatialAverage')
                  CALL netcdf_put_fvar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  self%pool(i)%spatialAverage,    &
     &                                  (/0/), (/0/),                   &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                CASE ('surveyIndex')
                  CALL netcdf_put_ivar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  self%pool(i)%surveyIndex,       &
     &                                  (/1/), (/nsurvey/),             &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                CASE ('surveyTime')
                  IF (.not.allocated(Iwork)) THEN
                    allocate ( Iwork(nsurvey) )
                  END IF
                  CALL netcdf_get_ivar (ng, model,                      &
     &                                  self%pool(i)%ncname_inp,        &
     &                                  TRIM(varname), Iwork,           &
     &                                  GrpName = TRIM(grpname))
                  IF (FoundError(exit_flag, NoError,                    &
     &                           __LINE__, MyFile)) RETURN

                  CALL netcdf_put_ivar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  Iwork,                          &
     &                                  (/1/), (/nsurvey/),             &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                  deallocate (Iwork)
                CASE ('x_grid')
                  CALL netcdf_put_fvar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  self%pool(i)%Xgrid,             &
     &                                  (/1/), (/nlocs/),               &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                CASE ('y_grid')
                  CALL netcdf_put_fvar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  self%pool(i)%Ygrid,             &
     &                                  (/1/), (/nlocs/),               &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
                CASE ('z_grid')
                  CALL netcdf_put_fvar (ng, model, TRIM(ncname),        &
     &                                  TRIM(varname),                  &
     &                                  self%pool(i)%Zgrid,             &
     &                                  (/1/), (/nlocs/),               &
     &                                  ncid = ncid,                    &
     &                                  GrpName = TRIM(grpname))
              END SELECT
            ELSE IF (TRIM(grpname).eq.'ObsError') THEN
              CALL netcdf_put_fvar (ng, model, TRIM(ncname),            &
     &                              TRIM(varname),                      &
     &                              self%pool(i)%error(:,j),            &
     &                              (/1/), (/nlocs/),                   &
     &                              ncid = ncid,                        &
     &                              GrpName = TRIM(grpname))
            ELSE IF (TRIM(grpname).eq.'ObsValue') THEN
              CALL netcdf_put_fvar (ng, model, TRIM(ncname),            &
     &                              TRIM(varname),                      &
     &                              self%pool(i)%value(:,j),            &
     &                              (/1/), (/nlocs/),                   &
     &                              ncid = ncid,                        &
     &                              GrpName = TRIM(grpname))
            ELSE IF (TRIM(grpname).eq.'PreQC') THEN
              IF (.not.allocated(Iwork)) THEN
                allocate ( Iwork(nlocs) )
               END IF
               CALL netcdf_get_ivar (ng, model,                         &
     &                               self%pool(i)%ncname_inp,           &
     &                               TRIM(varname), Iwork,              &
     &                               GrpName = TRIM(grpname))
               IF (FoundError(exit_flag, NoError,                       &
     &                        __LINE__, MyFile)) RETURN
!
               CALL netcdf_put_ivar (ng, model, TRIM(ncname),           &
     &                               TRIM(varname),                     &
     &                               Iwork,                             &
     &                               (/1/), (/nlocs/),                  &
     &                               ncid = ncid,                       &
     &                               GrpName = TRIM(grpname))
               deallocate (Iwork)
            END IF
!
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) THEN
              IF (Master) WRITE (stdout,10) TRIM(grpname),              &
     &                                      TRIM(varname), TRIM(ncname)
              exit_flag=3
              RETURN
            END IF
          END DO
        END DO
!
!  Close output NetCDF.  HGA: Need to check why I need to close to
!  use ncdump while ROMS is still running.
!
        CALL netcdf_close (ng, model, ncid, ncname, .FALSE.)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END DO FILE_LOOP
!
 10   FORMAT (/,'OBS_POOL - NC_WRITE: unable to write, Group = ',a,     &
     &        ', Variable = ',a,/,11x,'File: ',a)
!
      RETURN
      END SUBROUTINE obs_pool_nc_write
!
!-----------------------------------------------------------------------
!  It allocates the observation arrays in the appropiate pool vector
!  element, and then reads the specified native observation file to
!  populate the data.
!
      SUBROUTINE obs_pool_populate_native (self, ng, model, mobs,       &
     &                                     IOtype, ncname,              &
     &                                     typeID, Zgrid)
!
      CLASS (obs_pool),  intent(inout) :: self       ! observation set
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: model      ! kernel identifier
      integer,           intent(in   ) :: mobs       ! datum dimension
      integer,           intent(in   ) :: IOtype     ! NetCDF IO type
      character (len=*), intent(in   ) :: ncname     ! NetCDF filename
      integer,           intent(in   ) :: typeID(:)  ! observation type
      real (r8),         intent(in   ) :: Zgrid(:)   ! fractional Zgrid
!
      logical                          :: Loaded, isDayUnits
      integer,         dimension(mobs) :: obs_idata
      integer,             allocatable :: surveyIndex(:)
      integer                          :: Ninp, Nout
      integer                          :: i, iobs, lstr, mvar, npool
      real (r8),       dimension(mobs) :: obs_rdata
      real (r8),           allocatable :: surveyTime(:)
      real (r8),          dimension(2) :: Tmin, Tmax
      character (len=22), dimension(2) :: Dmin, Dmax
      character (len=40)               :: name, varname(2), vars(2)

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_pool_populate_native"
!
      SourceFile=MyFile
!
!  Read in observation data and populate.
!
      IF (IOtype.eq.io_nf90) THEN
        CALL netcdf_inq_var (ng, model, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        DO i=1,n_var
          Loaded=.FALSE.
          SELECT CASE (TRIM(var_name(i)))
            CASE ('obs_provenance')
              CALL netcdf_get_ivar (ng, model, ncname,                  &
     &                              var_name(i), obs_idata)
!
              CALL self%unpack_i (ng, mobs, typeID, Zgrid,              &
     &                            obs_idata, TRIM(var_name(i)))
            CASE ('obs_time')
              CALL netcdf_get_time (ng, model, ncname,                  &
     &                              var_name(i), Rclock%DateNumber,     &
     &                              obs_rdata)
!
              CALL self%unpack_r (ng, mobs, typeID, Zgrid,              &
     &                            obs_rdata, TRIM(var_name(i)))
            CASE ('obs_depth',                                          &
     &            'obs_lon',                                            &
     &            'obs_lat',                                            &
     &            'obs_Xgrid',                                          &
     &            'obs_Ygrid',                                          &
     &            'obs_Zgrid',                                          &
     &            'obs_error',                                          &
     &            'obs_value')
              CALL netcdf_get_fvar (ng, model, ncname,                  &
     &                              TRIM(var_name(i)), obs_rdata)
!
              CALL self%unpack_r (ng, mobs, typeID, Zgrid,              &
     &                            obs_rdata, TRIM(var_name(i)))
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
# if defined PIO_LIB && defined DISTRIBUTE
      ELSE IF (IOtype.eq.io_pio) THEN
        CALL pio_netcdf_inq_var (ng, model, ncname)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        DO i=1,n_var
          Loaded=.FALSE.
          SELECT CASE (TRIM(var_name(i)))
            CASE ('obs_provenance')
              CALL pio_netcdf_get_ivar (ng, model, ncname,              &
     &                                  var_name(i), obs_idata)
!
              CALL self%unpack_i (ng, mobs, typeID, Zgrid,              &
     &                            obs_idata, TRIM(var_name(i)))
            CASE ('obs_time')
              CALL pio_netcdf_get_time (ng, model, ncname,              &
     &                                  var_name(i), Rclock%DateNumber, &
     &                                  obs_rdata)
!
              CALL pio_self%unpack_r (ng, mobs, typeID, Zgrid,          &
     &                                obs_rdata, TRIM(var_name(i)))
            CASE ('obs_depth',                                          &
     &            'obs_lon',                                            &
     &            'obs_lat',                                            &
     &            'obs_Xgrid',                                          &
     &            'obs_Ygrid',                                          &
     &            'obs_Zgrid',                                          &
     &            'obs_error',                                          &
     &            'obs_value')
              CALL pio_netcdf_get_fvar (ng, model, ncname,              &
     &                                  TRIM(var_name(i)), obs_rdata)
!
              CALL self%unpack_r (ng, mobs, typeID, Zgrid,              &
     &                            obs_rdata, TRIM(var_name(i)))
          END SELECT
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END DO
# endif
      END IF
!
!  Determine surveyTime and surveyIndex.
!
      npool=SIZE(self%pool)

      DO iobs=1,npool
        Ninp=SIZE(self%pool(iobs)%time)
        CALL unique (Ninp, self%pool(iobs)%time,                        &
     &               Nout, surveyTime, surveyIndex)
!
        self%pool(iobs)%nsurvey=Nout
!
        IF (.not.allocated(self%pool(iobs)%surveyTime)) THEN
          allocate ( self%pool(iobs)%surveyTime(Nout) )
        END IF
        self%pool(iobs)%surveyTime=surveyTime
!
        IF (.not.allocated(self%pool(iobs)%surveyIndex)) THEN
          allocate ( self%pool(iobs)%surveyIndex(Nout) )
        END IF
        self%pool(iobs)%surveyTime=surveyTime
!
        deallocate (SurveyTime, SurveyIndex)
!
        Tmin(1)=MINVAL(self%pool(iobs)%time)
        Tmax(1)=MAXVAL(self%pool(iobs)%time)
        CALL datestr (Tmin(1), .TRUE., Dmin(1))
        CALL datestr (Tmax(1), .TRUE., Dmax(1))
!
        Tmin(2)=MINVAL(self%pool(iobs)%surveyTime)
        Tmax(2)=MAXVAL(self%pool(iobs)%surveyTime)
        CALL datestr (Tmin(2), .TRUE., Dmin(2))
        CALL datestr (Tmax(2), .TRUE., Dmax(2))
!
        IF (.not.allocated(self%pool(iobs)%ncname_inp)) THEN
          lstr=LEN_TRIM(ncname)
          allocate ( character(LEN=lstr) :: self%pool(iobs)%ncname_inp )
          self%pool(iobs)%ncname_inp=TRIM(ncname)
        END IF
!
        SELECT CASE (self%stateID(iobs))
          CASE (1)
            mvar=1
            name='SSH'
            varname(mvar)='SSH'
            vars(mvar)='absoluteDynamicTopography'
            self%pool(iobs)%stateID(mvar)=1
          CASE (6)
            mvar=1
            name='temp'
            varname(mvar)='temp'
            vars(mvar)='waterTemperature'
            self%pool(iobs)%stateID(mvar)=6
          CASE (7)
            mvar=1
            name='salt'
            varname(mvar)='salt'
            vars(mvar)='Salinity'
            self%pool(iobs)%stateID(mvar)=7
          CASE (45)
            mvar=2
            name='uv'
            varname(1)='u'
            varname(2)='v'
            vars(1)='waterSurfaceZonalVelocity'
            vars(2)='waterSurfaceMeridionalVelocity'
            self%pool(iobs)%stateID(1)=4
            self%pool(iobs)%stateID(2)=5
          CASE (60)
            mvar=1
            name='SST'
            varname(mvar)='SST'
            vars(mvar)='seaSurfaceTemperature'
            self%pool(iobs)%stateID(1)=6
          CASE (70)
            mvar=1
            name='SSS'
            varname(mvar)='SSS'
            vars(mvar)='seaSurfaceSalinity'
            self%pool(iobs)%stateID(1)=7
        END SELECT
!
        IF (.not.allocated(self%pool(iobs)%name)) THEN
          lstr=LEN_TRIM(name)
          allocate ( character(LEN=lstr) :: self%pool(iobs)%name )
          self%pool(iobs)%name=TRIM(name)
        END IF
!
        IF (.not.allocated(self%pool(iobs)%vname)) THEN
          lstr=0
          DO i=1,mvar
            lstr=MAX(lstr, LEN_TRIM(varname(i)))
          END DO
          allocate (character(LEN=lstr) :: self%pool(iobs)%vname(mvar) )
          DO i=1,mvar
            self%pool(iobs)%vname(i)=TRIM(varname(i))
          END DO
        END IF
!
        IF (.not.allocated(self%pool(iobs)%vars)) THEN
          lstr=0
          DO i=1,mvar
            lstr=MAX(lstr, LEN_TRIM(vars(i)))
          END DO
          allocate ( character(LEN=lstr) :: self%pool(iobs)%vars(mvar) )
          DO i=1,mvar
            self%pool(iobs)%vars(i)=TRIM(vars(i))
          END DO
        END IF
!
!  Compute statistics.
!
        isDayUnits=.TRUE.
        CALL self%pool(iobs)%stats (isDayUnits)
!
!  Report processing information.
!
        IF (Master) THEN
          WRITE (stdout,10) TRIM(ncname)
          WRITE (stdout,20) 'MetaData', 'obs_time',      Ninp
          WRITE (stdout,20) 'MetaData', 'obs_depth',     Ninp
          WRITE (stdout,20) 'MetaData', 'obs_lat',       Ninp
          WRITE (stdout,20) 'MetaData', 'obs_lon',       Ninp
          WRITE (stdout,20) 'MetaData', 'obs_provenace', Ninp
          WRITE (stdout,20) 'MetaData', 'surveyIndex',   Nout
          WRITE (stdout,20) 'MetaData', 'obs_survey',    Nout
          WRITE (stdout,20) 'MetaData', 'obs_Xgrid',     Ninp
          WRITE (stdout,20) 'MetaData', 'obs_Ygrid',     Ninp
          WRITE (stdout,20) 'MetaData', 'obs_Zgrid',     Ninp
          DO i=1,mvar
            WRITE (stdout,30) 'ObsError', TRIM(vars(i)), Ninp,          &
     &                        self%pool(iobs)%stateID(i)
          END DO
          DO i=1,mvar
            WRITE (stdout,30) 'ObsValue', TRIM(vars(i)), Ninp,          &
     &                        self%pool(iobs)%stateID(i)
          END DO
!
          WRITE (stdout,'(1x)')
          DO i=1,mvar
            WRITE (stdout,40) TRIM(self%pool(iobs)%vars(i)),            &
     &                        self%pool(iobs)%value_min(i),             &
     &                        self%pool(iobs)%value_max(i),             &
     &                        self%pool(iobs)%value_avg(i),             &
     &                        self%pool(iobs)%value_hash(i),            &
     &                        self%pool(iobs)%error_min(i),             &
     &                        self%pool(iobs)%error_max(i),             &
     &                        self%pool(iobs)%error_avg(i)
            WRITE (stdout,50) 'DateTime Range:   ',                     &
     &                        TRIM(Dmin(1)), TRIM(Dmax(1))
            WRITE (stdout,50) 'SurveyTime Range: ',                     &
     &                        TRIM(Dmin(2)), TRIM(Dmax(2))
          END DO
        END IF
      END DO
!
  10  FORMAT (/,2x,'ROMS_OBS         - Processing: ',a)
  20  FORMAT (t23,a,'::',a,t75,'Count = ',i0)
  30  FORMAT (t23,a,'::',a,t75,'Count = ',i0,t90,'ID = ',i0)
  40  FORMAT (t21,'- ',a,/,t23,'(ValMin = ',1p,e12.5,' ValMax = ',      &
     &        1p,e12.5,' ValAvg = ',1p,e12.5,' CheckSum = ',i0,')',     &
     &        /,t23,'(ErrMin = ',1p,e12.5,' ErrMax = ',1p,e12.5,        &
     &        ' ErrAvg = ',1p,e12.5,')')
  50  FORMAT (t23,a,a,' to ',a)
!
      RETURN
      END SUBROUTINE obs_pool_populate_native
!
!-----------------------------------------------------------------------
!  Define the size of the observation pool vector as the number of
!  variables in the state control vector that are directly or
!  indirectly linked to the measurements. Retrieve the state variable
!  stateID (1:nvars) from the input observations file(s), nfiles.
!  Then, allocate observation pool vector.
!
      SUBROUTINE obs_pool_size (self, ng, model, nfiles)
!
      CLASS (obs_pool), intent(inout) :: self     ! observation set
      integer,          intent(in   ) :: ng       ! nested grid
      integer,          intent(in   ) :: model    ! kernel identifier
      integer,          intent(in   ) :: nfiles   ! number of files
!
      logical                         :: haveID(5+NT(ng), nfiles)
      logical,      dimension(nfiles) :: haveSSS, haveSST, haveVector
      logical,      dimension(nfiles) :: add_Sid, add_Tid, isNative
      integer,      dimension(nfiles) :: mobs, msurvey, mvar
      integer                         :: stateID(5+NT(ng), nfiles)
      integer                         :: uniqueVal(5+NT(ng))
      integer,            allocatable :: myStateID(:)
      integer,            allocatable :: obs_type(:)
      integer                         :: i, ic, mdatum, mvatt, nf, nvars
      integer                         :: id, mSalt, mTemp, mSSS, mSST
      integer                         :: Lsur
      real (r8),          allocatable :: obs_Zgrid(:)

      character (len=*),    parameter :: MyFile =                       &
     &    __FILE__//", obs_pool_size"
!
      SourceFile=MyFile
!
!  Initialize.
!
      add_Sid=.FALSE.
      add_Tid=.FALSE.
      isNative=.FALSE.
      haveID=.FALSE.
      haveSSS=.FALSE.
      haveSST=.FALSE.
      haveVector=.FALSE.
!
      mobs=0
      mSSS=0
      mSST=0
      msurvey=0
      mvar=0
      mvatt=0
      stateID=0
!
!     Lsur=N(ng)
      Lsur=50                    ! HGA: teporary for debugging
!
!  Inquire input observation NetCDF file(s).
!
      FILES_LOOP : DO nf=1,nfiles
        SELECT CASE (OBS2(nf,ng)%IOtype)
          CASE (io_nf90)
            CALL netcdf_get_dim (ng, model, OBS2(nf,ng)%name)
# if defined PIO_LIB && defined DISTRIBUTE
          CASE (io_pio)
            CALL pio_netcdf_get_dim (ng, model, OBS2(nf,ng)%name)
# endif
          CASE DEFAULT
            IF (Master) WRITE (stdout,10) OBS2(i,ng)%IOtype
            exit_flag=3
        END SELECT
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Get number of observation from its dimensions, and determine if
!  native or IODA file.
!
        DO i=1,n_dim
          IF (TRIM(dim_name(i)).eq.'datum') then
            mdatum=dim_size(i)
            mobs(nf)=mobs(nf)+mdatum
            isNative(nf)=.TRUE.                         ! only in native
          ELSE IF (TRIM(dim_name(i)).eq.'survey') THEN
            msurvey(nf)=msurvey(nf)+dim_size(i)         ! only in native
          ELSE IF (TRIM(dim_name(i)).eq.'Location') THEN
            mobs(nf)=mobs(nf)+dim_size(i)
            isNative(nf)=.FALSE.                        ! only in IODA
          ELSE IF (TRIM(dim_name(i)).eq.'nvars') THEN
            mvar(nf)=mvar(nf)+dim_size(i)               ! only in IODA
          END IF
        END DO
!
!  If native observation file, read in "obs_type" and "obs_Zgrid" to
!  determine the observed variables and its IDs.
!
        IF (isNative(nf)) THEN
          IF (.not.allocated(obs_type)) THEN
            allocate ( obs_type(mdatum) )
            obs_type=0
          END IF
          IF (.not.allocated(obs_Zgrid)) THEN
            allocate ( obs_Zgrid(mdatum) )
            obs_Zgrid=0.0_r8
          END IF
!
          SELECT CASE (OBS2(nf,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_get_ivar (ng, model, OBS2(nf,ng)%name,        &
     &                              'obs_type', obs_type)
              IF (FoundError(exit_flag,NoError,__LINE__,MyFile)) RETURN
!
              CALL netcdf_get_fvar (ng, model, OBS2(nf,ng)%name,        &
     &                              'obs_Zgrid', obs_Zgrid)
              IF (FoundError(exit_flag,NoError,__LINE__,MyFile)) RETURN

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_get_ivar (ng, model, OBS2(nf,ng)%name,    &
     &                                  'obs_type', obs_type)
              IF (FoundError(exit_flag,NoError,__LINE__,MyFile)) RETURN
!
              CALL pio_netcdf_get_fvar (ng, model, OBS2(nf,ng)%name,    &
     &                                  'obs_Zgrid', obs_Zgrid)
              IF (FoundError(exit_flag,NoError,__LINE__,MyFile)) RETURN
# endif
          END SELECT
!
          uniqueVal=0
          DO i=1,mdatum
            id=obs_type(i)
            uniqueVal(id)=uniqueVal(id)+1
            haveID(id,nf)=.TRUE.
            IF ((id.eq.6).and.(INT(obs_Zgrid(i)).eq.Lsur)) THEN
              haveSST(nf)=.TRUE.
              mSST=mSST+1
            ELSE IF ((id.eq.6).and.(INT(obs_Zgrid(i)).eq.Lsur)) THEN
              haveSSS(nf)=.TRUE.
              mSSS=mSSS+1
            END IF
          END DO
          mvar(nf)=mvar(nf)+COUNT(uniqueVal.gt.0)
          mTemp=uniqueVal(6)
          mSalt=uniqueVal(7)
          IF (haveSST(nf).and.(mSST.gt.0).and.(mTemp.gt.mSST)) THEN
            add_Tid(nf)=.TRUE.
          END IF
          IF (haveSSS(nf).and.(mSSS.gt.0).and.(mSalt.gt.mSSS)) THEN
            add_Sid(nf)=.TRUE.
          END IF
          deallocate (obs_type)
          deallocate (obs_Zgrid)
!
!  Otherwise, get variable ID for IODA file.
!
        ELSE
          allocate ( myStateID(mvar(nf)) )
          myStateID=0
!
          SELECT CASE (OBS2(nf,ng)%IOtype)
            CASE (io_nf90)
              CALL netcdf_get_ivar (ng, model, OBS2(nf,ng)%name,        &
     &                              'stateID', myStateID,               &
     &                              GrpName = 'MetaData')
              IF (FoundError(exit_flag,NoError,__LINE__,MyFile)) RETURN
!
              CALL netcdf_inq_var (ng, model, OBS2(nf,ng)%name,         &
     &                             myVarName = 'stateID',               &
     &                             GrpName = 'MetaData',                &
     &                             nVarAtt = mvatt)
              IF (FoundError(exit_flag,NoError,__LINE__,MyFile)) RETURN
!
              DO i=1,mvatt
                IF (TRIM(var_Aname(i)).eq.'surface_level') THEN
                  IF (ANY(myStateID.eq.6)) THEN
                    haveSST(nf)=.TRUE.
                  END IF
                  IF (ANY(myStateID.eq.7)) THEN
                    haveSSS(nf)=.TRUE.
                  END IF
                END IF
              END DO

# if defined PIO_LIB && defined DISTRIBUTE
            CASE (io_pio)
              CALL pio_netcdf_get_ivar (ng, model, OBS2(nf,ng)%name,    &
     &                                  'stateID', myStateID,           &
     &                                  GrpName = 'MetaData')
              IF (FoundError(exit_flag,NoError,__LINE__,MyFile)) RETURN
# endif
          END SELECT
!
          DO i=1,mvar(nf)
            id=myStateID(i)
            haveID(id,nf)=.TRUE.
          END DO
          deallocate (myStateID)
        END IF
!
!  Set number of state variables and their IDs.
!
        ic=0
        DO i=1,5+NT(ng)
          IF (haveID(i,nf)) THEN
            ic=ic+1
            StateID(ic,nf)=i
          END IF
        END DO
!
!  Set number of control vector variables. Since we have multiple
!  observation files, the logic in StateID already accounts for
!  SST and SSS, and there is not need for augmentation of "nvars",
!  as in the native case.
!
        nvars=COUNT(StateID.gt.0)
        IF (isNative(nf)) THEN
          IF (haveSST(nf)) nvars=nvars+1
          IF (haveSSS(nf)) nvars=nvars+1
        END IF

      END DO FILES_LOOP
!
!  If velocity observations, myStateID=4 and myStateID=5, adjust the
!  "nvars" values because vector components are stored in the same pool
!  vector element since they share the same location and time coordinates.
!  Note that stateID=45 indicates that the element now contains the
!  velocity vector components.
!
      DO nf=1,nfiles
        IF (ANY(StateID(:,nf).eq.4).and.ANY(StateID(:,nf).eq.5)) THEN
          haveVector(nf)=.TRUE.                          ! has ID = 4, 5
          nvars=nvars-1
        END IF
      END DO
!
!  Alocate and initialize the stateID vector.
!
      IF (.not.allocated(self%stateID)) THEN
        allocate ( self%stateID(nvars) )
      END IF
!
      ic=0
      DO nf=1,nfiles
        DO i=1,5+NT(ng)
          IF (StateID(i,nf).gt.0) THEN
            IF (haveVector(nf).and.                                     &
     &          (StateID(i,nf).eq.4).or.(StateID(i,nf).eq.5)) THEN
              IF (StateID(i,nf).eq.5) CYCLE            ! skip assigment
              ic=ic+1
              self%stateID(ic)=45                      ! has ID = 4, 5
            ELSE IF (haveSST(nf).and.(stateID(i,nf).eq.6)) THEN
              ic=ic+1
              self%stateID(ic)=60                      ! has ID = 6, z=0
              IF (isNative(nf).and.add_Tid(nf)) THEN
                ic=ic+1
                self%stateID(ic)=6                     ! has ID = 6
              END IF
            ELSE IF (haveSSS(nf).and.(stateID(i,nf).eq.7)) THEN
              ic=ic+1
              self%stateID(ic)=70                      ! has ID = 7, z=0
              IF (isNative(nf).and.add_Sid(nf)) THEN
                ic=ic+1
                self%stateID(ic)=7                     ! has ID = 7
              END IF
            ELSE
              ic=ic+1
              self%stateID(ic)=StateID(i,nf)
            END IF
          END IF
        END DO
      END DO
!
!  Allocate observation pool vector.
!
      IF (.not.allocated(self%pool)) THEN
        allocate ( self%pool(nvars) )
      END IF
!
  10  FORMAT (' OBS_POOL_SIZE - Illegal output file type, io_type = ',  &
     &        i0,/,17x,'Check KeyWord ''OUT_LIB'' in ''roms.in''.')
!
      RETURN
      END SUBROUTINE obs_pool_size
!
!-----------------------------------------------------------------------
!  It determines the total number of observatos in the pool depot.
!
      SUBROUTINE obs_pool_tally (self)
!
      CLASS (obs_pool),  intent(inout) :: self       ! observation set
!
      integer,             allocatable :: ObsCount(:)
      integer                          :: mvars, ndatum, npool
      integer                          :: i, ic, j, m, nc
!
!  Count the total number of observations.
!
      ndatum=0
      npool=SIZE(self%pool)
      IF (.not.allocated(ObsCount)) THEN
        allocate ( ObsCount(npool) )
        ObsCount=0
      END IF
!
      DO i=1,npool
        ObsCount(i)=self%pool(i)%nlocs
        ndatum=ndatum+ObsCount(i)*self%pool(i)%nvars
      END DO
      self%ndatum=ndatum
!
!  Set mapping indices to packed observations vector. It is used in the
!  H(x) object.
!
      ic=0
      DO i=1,npool
        mvars=self%pool(i)%nvars
        IF (.not.allocated(self%pool(i)%pack_mapping)) THEN
          allocate ( self%pool(i)%pack_mapping(ObsCount(i),mvars) )
          self%pool(i)%pack_mapping=0
        END IF
        DO m=1,mvars
          nc=0
          DO j=1,self%pool(i)%nlocs
            ic=ic+1
            nc=nc+1
            self%pool(i)%pack_mapping(nc,m)=ic
          END DO
        END DO
      END DO
      IF (allocated(ObsCount)) deallocate (ObsCount)
!
      RETURN
      END SUBROUTINE obs_pool_tally
!
!-----------------------------------------------------------------------
!  It unpacks native observation integer data into pool object.
!
      SUBROUTINE obs_pool_unpack_i (self, ng, ND, typeID, Zgrid, idata, &
     &                              VarName)
!
      CLASS (obs_pool),  intent(inout) :: self       ! observation set
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: ND         ! datum dimension
      integer,           intent(in   ) :: typeID(ND) ! observation type
      real (r8),         intent(in   ) :: Zgrid(ND)  ! fractional Zgrid
      integer,           intent(in   ) :: idata(ND)  ! obsevation data
      character (len=*), intent(in   ) :: VarName    ! variable name
!
      logical                          :: mask(ND)
      integer                          :: Lsur, iobs, mvar, npool
      integer                          :: poolID

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_pool_unpack_i"
!
      SourceFile=MyFile
!
!  Unpack native observation integer data.
!
!     Lsur=N(ng)
      Lsur=50                    ! HGA: temporary for debugging
!
      mask=.FALSE.
      mvar=0
      npool=SIZE(self%pool)
!
      DO iobs=1,npool
        poolID=self%pool(iobs)%poolID
!
        SELECT CASE (poolID)
          CASE (1)                                      ! SSH
            mask=(typeID.eq.1)
            mvar=1
          CASE (6)                                      ! Temperature
            mask=(typeID.eq.6).and.(INT(Zgrid).ne.Lsur)
            mvar=1
          CASE (7)                                      ! Salinity
            mask=(typeID.eq.7)
            mvar=1
          CASE (45)                                     ! U-,V-velocity
            mask=(typeID.eq.4)                          ! U-Velocity
            mvar=2
          CASE (60)                                     ! SST, z=0
            mask=(typeID.eq.6).and.(INT(Zgrid).eq.Lsur)
            mvar=1
          CASE (70)                                     ! SSS, z=0
            mask=(typeID.eq.7).and.(INT(Zgrid).eq.Lsur)
            mvar=1
          CASE DEFAULT
            IF (Master) WRITE (stdout,10) poolID
            exit_flag=5
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END SELECT
!
!  Load observation data to the appropiate array in the pool object.
!
        SELECT CASE (TRIM(VarName))
          CASE ('obs_provenance')
            self%pool(iobs)%provenance = PACK(idata, mask)
          CASE DEFAULT
            IF (Master) WRITE (stdout,20) TRIM(VarName)
            exit_flag=5
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END SELECT
      END DO
!
  10  FORMAT (/,'ROMS_OBS - obs_pool::unpack_i: cannot find ',          &
     &        'poolID = ',i0)
  20  FORMAT (/,'ROMS_OBS - obs_pool::unpack_i: cannot find ',          &
     &        'VarName = ',a)
!
      RETURN
      END SUBROUTINE obs_pool_unpack_i
!
!-----------------------------------------------------------------------
!  It unpacks native observation floatin-point data into pool object.
!
      SUBROUTINE obs_pool_unpack_r (self, ng, ND, typeID, Zgrid, rdata, &
     &                              VarName)
!
      CLASS (obs_pool),  intent(inout) :: self       ! observation set
      integer,           intent(in   ) :: ng         ! nested grid
      integer,           intent(in   ) :: ND         ! datum dimension
      integer,           intent(in   ) :: typeID(ND) ! observation type
      real (r8),         intent(in   ) :: Zgrid(ND)  ! fractional Zgrid
      real (r8),         intent(in   ) :: rdata(ND)  ! obsevation data
      character (len=*), intent(in   ) :: VarName    ! variable name
!
      logical                          :: mask(ND), msk2(ND)
      integer                          :: Lsur, i, iobs, mvar, npool
      integer                          :: poolID

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_pool_unpack_r"
!
      SourceFile=MyFile
!
!  Unpack native observation integer data.
!
!     Lsur=N(ng)
      Lsur=50                    ! HGA: temporary for debugging
!
      mask=.FALSE.
      msk2=.FALSE.
      mvar=0
      npool=SIZE(self%pool)
!
      DO iobs=1,npool
        poolID=self%pool(iobs)%poolID
!
        SELECT CASE (poolID)
          CASE (1)                                      ! SSH
            mask=(typeID.eq.1)
            mvar=1
          CASE (6)                                      ! Temperature
            mask=(typeID.eq.6).and.(INT(Zgrid).ne.Lsur)
            mvar=1
          CASE (7)                                      ! Salinity
            mask=(typeID.eq.7)
            mvar=1
          CASE (45)                                     ! U-,V-velocity
            mask=(typeID.eq.4)                          ! U-velocity
            msk2=(typeID.eq.5)                          ! V-velocity
            mvar=2
          CASE (60)                                     ! SST. z=0
            mask=(typeID.eq.6).and.(INT(Zgrid).eq.Lsur)
            mvar=1
          CASE (70)                                     ! SSS, z=0
            mask=(typeID.eq.7).and.(INT(Zgrid).eq.Lsur)
            mvar=1
          CASE DEFAULT
            IF (Master) WRITE (stdout,10) poolID
            exit_flag=5
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END SELECT
!
!  Load observation data to the appropiate array in the pool object.
!
        SELECT CASE (TRIM(VarName))
          CASE ('obs_time')
            self%pool(iobs)%time  = PACK(rdata, mask)
          CASE ('obs_depth')
            self%pool(iobs)%depth = PACK(rdata, mask)
          CASE ('obs_lon')
            self%pool(iobs)%lon   = PACK(rdata, mask)
          CASE ('obs_lat')
            self%pool(iobs)%lat   = PACK(rdata, mask)
          CASE ('obs_Xgrid')
            self%pool(iobs)%Xgrid = PACK(rdata, mask)
          CASE ('obs_Ygrid')
            self%pool(iobs)%Ygrid = PACK(rdata, mask)
          CASE ('obs_Zgrid')
            self%pool(iobs)%Zgrid = PACK(rdata, mask)
          CASE ('obs_error')
            DO i=1,mvar
              self%pool(iobs)%error(:,i)   = PACK(rdata, mask)
              IF (i.eq.2) THEN
                self%pool(iobs)%error(:,i) = PACK(rdata, msk2)
              END IF
            END DO
          CASE ('obs_value')
            DO i=1,mvar
              self%pool(iobs)%value(:,i)   = PACK(rdata, mask)
              IF (i.eq.2) THEN
                self%pool(iobs)%value(:,i) = PACK(rdata, msk2)
              END IF
            END DO
          CASE DEFAULT
            IF (Master) WRITE (stdout,20) TRIM(VarName)
            exit_flag=5
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END SELECT
      END DO
!
  10  FORMAT (/,'ROMS_OBS - obs_pool::unpack_r: cannot find ',          &
     &        'poolID = ',i0)
  20  FORMAT (/,'ROMS_OBS - obs_pool::unpack_r: cannot find ',          &
     &        'VarName = ',a)
!
      RETURN
      END SUBROUTINE obs_pool_unpack_r
#endif
      END MODULE roms_obs_mod
