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
      USE distribute_mod,       ONLY : mp_assemble
#endif
      USE get_hash_mod,         ONLY : get_hash
      USE roms_interpolate_mod, ONLY : hindices
      USE unique_mod,           ONLY : unique
!     
      implicit none
!
!-----------------------------------------------------------------------
!  Observation Data Type Structures/Objects: CLASS(obs_data).
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: obs_data

        integer :: nlocs                        ! number of observations
        integer :: nsurvey                      ! number of time surveys
        integer :: nvars                        ! number of variables
        integer :: poolID                       ! observation pool ID
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
!  Observation spatial and time coordinates. The fractional values are
!  used for efficient interpolation in the H(x) operator.
!
        real(r8), allocatable :: lon(:)         ! observation longitude
        real(r8), allocatable :: lat(:)         ! observation latitude
        real(r8), allocatable :: depth(:)       ! depth of observation
        real(r8), allocatable :: time(:)        ! time  of observation
        real(r8), allocatable :: surveyTime(:)  ! survey time
!
        real(r8), allocatable :: Xgrid(:)       ! fractional X-grid
        real(r8), allocatable :: Ygrid(:)       ! fractional Y-grid
        real(r8), allocatable :: Zgrid(:)       ! fractional Z-grid
!
!  If the observations are scalar, the second dimension is unity.
!  Otherwise, it is used for vectors or multiple observations at the
!  exact spatial and time coordinates.
!
        real(r8), allocatable :: error(:,:)     ! observation error(s)
        real(r8), allocatable :: value(:,:)     ! observation value(s)
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
!
!  Input observation filename. It supportss native ROMS 4D-Var NetCDF
!  files and IODA enhanced NetCDF-4 files with Groups.
!
        character (len=:), allocatable :: ncname
!
        CONTAINS
!
        PROCEDURE :: destroy          => obs_data_destroy
        PROCEDURE :: frac_coords      => obs_data_frac_coords
        PROCEDURE :: populate_ioda    => obs_data_populate_ioda
        PROCEDURE :: stats            => obs_data_stats
!
      END TYPE obs_data
!
!-----------------------------------------------------------------------
!  Collection Observation Data Types Object: CLASS(obs_pool).
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: obs_pool
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
        PROCEDURE :: populate_native  => obs_pool_populate_native
        PROCEDURE :: read             => obs_pool_read
        PROCEDURE :: size             => obs_pool_size
        PROCEDURE :: unpack_i         => obs_pool_unpack_i
        PROCEDURE :: unpack_r         => obs_pool_unpack_r
!
      END TYPE obs_pool
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      PUBLIC  :: roms_obs_initialize
      PRIVATE
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      CONTAINS
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!  It initialzes the observation set object.
!
      SUBROUTINE roms_obs_initialize (self, ng, model, nfiles)
!
      CLASS (obs_pool), intent(inout) :: self     ! observation set
      integer,          intent(in   ) :: ng       ! nested grid number
      integer,          intent(in   ) :: model    ! kernel identifier
      integer,          intent(in   ) :: nfiles   ! number of files
!
      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", roms_obs_initialize"
!
      SourceFile=MyFile
!
!  Creates observation pool object.
!
      CALL self%create (ng, model, nfiles)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      RETURN
      END SUBROUTINE roms_obs_initialize
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
      IF (allocated(self%provenance))  deallocate (self%provenance)
      IF (allocated(self%lon))         deallocate (self%lon)
      IF (allocated(self%lat))         deallocate (self%lat)
      IF (allocated(self%depth))       deallocate (self%depth)
      IF (allocated(self%time))        deallocate (self%time)
      IF (allocated(self%surveyIndex))  deallocate (self%surveyIndex)
      IF (allocated(self%surveyTime))  deallocate (self%surveyTime)
      IF (allocated(self%Xgrid))       deallocate (self%Xgrid)
      IF (allocated(self%Ygrid))       deallocate (self%Ygrid) 
      IF (allocated(self%Zgrid))       deallocate (self%Zgrid)
!
      IF (allocated(self%error))       deallocate (self%error)
      IF (allocated(self%value))       deallocate (self%value)
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
      integer,          intent(in   ) :: ng       ! nested grid number
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
!  If not available, it computes observation object fractional
!  coordinates.
!
      SUBROUTINE obs_data_populate_ioda (self, ng, model, mvar, IOtype)
!
      CLASS (obs_data),  intent(inout) :: self    ! observation object
      integer,           intent(in   ) :: ng      ! nested grid number
      integer,           intent(in   ) :: model   ! kernel identifier
      integer,           intent(in   ) :: mvar    ! number of variables
      integer,           intent(in   ) :: IOtype  ! NetCDF IO type
!
      logical                          :: Loaded, isDayUnits
      integer                          :: i, lstr, ncount
      real(r8),           dimension(2) :: Tmin, Tmax
      character (len=22), dimension(2) :: Dmin, Dmax

      character (len=*),     parameter :: MyFile =                      &
     &    __FILE__//", obs_data_populate_ioda"
!
      SourceFile=MyFile
!
!  Read in observation variables in Group MetaData.
!
      IF (Master) WRITE (stdout,10) TRIM(self%ncname)
!      
      IF (IOtype.eq.io_nf90) THEN
        CALL netcdf_inq_var (ng, model, self%ncname,                    &
                             GrpName = 'MetaData')
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        DO i=1,n_var
          Loaded=.FALSE.
          SELECT CASE (TRIM(var_name(i)))
            CASE ('dateTime')
              CALL netcdf_get_time (ng, model, self%ncname,             &
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
              CALL netcdf_get_fvar (ng, model, self%ncname,             &
     &                              var_name(i), self%depth,            &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%depth)
            CASE ('latitude')
              CALL netcdf_get_fvar (ng, model, self%ncname,             &
     &                              var_name(i), self%lat,              &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%lat)
            CASE ('longitude')
              CALL netcdf_get_fvar (ng, model, self%ncname,             &
     &                              var_name(i), self%lon,              &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%lon)
            CASE ('provenance')
              CALL netcdf_get_ivar (ng, model, self%ncname,             &
     &                              var_name(i), self%provenance,       &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%provenance)
            CASE ('stateID')
              CALL netcdf_get_ivar (ng, model, self%ncname,             &
     &                              var_name(i), self%stateID,          &
     &                              GrpName = 'MetaData')
            CASE ('surveyIndex')
              CALL netcdf_get_ivar (ng, model, self%ncname,             &
     &                              var_name(i), self%surveyIndex,      &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%surveyIndex)
            CASE ('surveyTime')
              CALL netcdf_get_time (ng, model, self%ncname,             &
     &                              var_name(i), Rclock%DateNumber,     &
     &                              self%surveyTime,                    &
     &                              GrpName = 'MetaData',               &
     &                              min_val = Tmin(2),                  &
     &                              max_val = Tmax(2))
              Loaded=.TRUE.
              ncount=SIZE(self%surveyTime)
              CALL datestr (Tmin(2), .FALSE., Dmin(2))
              CALL datestr (Tmax(2), .FALSE., Dmax(2))
            CASE ('x_grid')
              CALL netcdf_get_fvar (ng, model, self%ncname,             &
     &                              var_name(i), self%Xgrid,            &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%Xgrid)
            CASE ('y_grid')
              CALL netcdf_get_fvar (ng, model, self%ncname,             &
     &                              var_name(i), self%Ygrid,            &
     &                              GrpName = 'MetaData')
              Loaded=.TRUE.
              ncount=SIZE(self%Ygrid)
            CASE ('z_grid')
              CALL netcdf_get_fvar (ng, model, self%ncname,             &
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
!  Read in observation variables in Group ObsError.
!
      IF (IOtype.eq.io_nf90) THEN
        CALL netcdf_inq_var (ng, model, self%ncname,                    &
                             GrpName = 'ObsError')
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        DO i=1,n_var
          CALL netcdf_get_fvar (ng, model, self%ncname,                 &
     &                          var_name(i), self%error(:,i),           &
     &                          GrpName = 'ObsError')
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
        CALL netcdf_inq_var (ng, model, self%ncname,                    &
                             GrpName = 'ObsValue')
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
        lstr=0
        DO i=1,n_var
          CALL netcdf_get_fvar (ng, model, self%ncname,                 &
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
        CALL get_hash (self%error(:,i), self%nvars, checksum, .FALSE.)
        self%error_hash(i)=checksum
!
        self%value_avg(i)=SUM(self%value(:,i))/REAL(self%nlocs, r8)
        self%value_min(i)=MINVAL(self%value(:,i))
        self%value_max(i)=MAXVAL(self%value(:,i))
        CALL get_hash (self%value(:,i), self%nvars, checksum, .FALSE.)
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
      integer,          intent(in   ) :: ng       ! nested grid number
      integer,          intent(in   ) :: model    ! kernel identifier
      integer,          intent(in   ) :: nfiles   ! number of files
!
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
      CALL self%read (ng, model, nfiles)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
      RETURN
      END SUBROUTINE obs_pool_create
!
!-----------------------------------------------------------------------
!
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
!
!  Given internal observation name, it returns an observation data
!  object pointer.
!
      SUBROUTINE obs_pool_get (self, name, field)
!
      CLASS (obs_pool), target, intent(in ) :: self   ! observation set
      character (len=*),        intent(in ) :: name   ! observation name
      TYPE (obs_data), pointer, intent(out) :: field  ! data pointer
!
      integer                               :: i

      character (len=*),          parameter :: MyFile =                 &
     &    __FILE__//", obs_pool_get"
!
      SourceFile=MyFile
!
!  Get requested field pointer.
!
      DO i=1,SIZE(self%pool)
        IF (TRIM(name).eq.self%pool(i)%name) THEN
          field => self%pool(i)
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
!
!  It inquires if provided observation name exists in the pool object.
!
      FUNCTION obs_pool_has (self, name, findex)  RESULT (foundit)
!
      CLASS (obs_pool),  intent(in   ) :: self   ! observation set
      character (len=*), intent(in   ) :: name   ! observation name
      integer, optional, intent(out  ) :: findex ! field index in set
!
      logical                          :: foundit
      integer                          :: i
!
!  Check if requested observation name exists in the pool object.
!
      foundit=.FALSE.
      DO i=1,SIZE(self%pool)
        IF (TRIM(name).eq.self%pool(i)%name) THEN
          foundit=.TRUE.
          IF (PRESENT(findex)) findex=1
          RETURN
        END IF
      END DO
!
      END FUNCTION obs_pool_has
!
!-----------------------------------------------------------------------
!
!  It allocates the observation arrays in the appropiate pool vector
!  element, and then reads the specified IODA observation file to
!  populate the data.
!
      SUBROUTINE obs_pool_load_ioda (self, ng, model, mobs, msurvey,    &
     &                               mvar, IOtype, ncname)
!
      CLASS (obs_pool),  intent(inout) :: self     ! observation set
      integer,           intent(in   ) :: ng       ! nested grid number
      integer,           intent(in   ) :: model    ! kernel identifier
      integer,           intent(in   ) :: mobs     ! observations number
      integer,           intent(in   ) :: msurvey  ! time surveys
      integer,           intent(in   ) :: mvar     ! number of variables
      integer,           intent(in   ) :: IOtype   ! NetCDF IO type
      character (len=*), intent(in   ) :: ncname   ! NetCDF filename
!
      logical                          :: haveSST, haveSSS, haveVector
      integer,         dimension(mvar) :: myStateID
      integer                          :: i, iobs, lstr, myID, mvatt

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
!  Determine pool vector element to process.
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
      IF (.not.allocated(self%pool(iobs)%ncname)) THEN
        lstr=LEN_TRIM(ncname)
        allocate ( character(LEN=lstr) :: self%pool(iobs)%ncname )
        self%pool(iobs)%ncname = TRIM(ncname)
      END IF
!
!  Allocate pool vector element structure, Group MetaData.
!
      IF (.not.allocated(self%pool(iobs)%provenance)) THEN
        allocate ( self%pool(iobs)%provenance(mobs) )
        self%pool(iobs)%provenance = 0
      END IF

      IF (.not.allocated(self%pool(iobs)%lon)) THEN
        allocate ( self%pool(iobs)%lon(mobs) )
        self%pool(iobs)%lon = 0.0_r8
      END IF

      IF (.not.allocated(self%pool(iobs)%lat)) THEN
        allocate ( self%pool(iobs)%lat(mobs) )
        self%pool(iobs)%lat = 0.0_r8
      END IF

      IF (.not.allocated(self%pool(iobs)%depth)) THEN
        allocate ( self%pool(iobs)%depth(mobs) )
        self%pool(iobs)%depth = 0.0_r8
      END IF

      IF (.not.allocated(self%pool(iobs)%time)) THEN
        allocate ( self%pool(iobs)%time(mobs) )
        self%pool(iobs)%time = 0.0_r8
      END IF

      IF (.not.allocated(self%pool(iobs)%stateID)) THEN
        allocate ( self%pool(iobs)%stateID(mvar) )
        self%pool(iobs)%stateID = 0
      END IF

      IF (.not.allocated(self%pool(iobs)%surveyIndex)) THEN
        allocate ( self%pool(iobs)%surveyIndex(msurvey) )
        self%pool(iobs)%surveyIndex = 0
      END IF

      IF (.not.allocated(self%pool(iobs)%surveyTime)) THEN
        allocate ( self%pool(iobs)%surveyTime(msurvey) )
        self%pool(iobs)%surveyTime = 0.0_r8
      END IF

      IF (.not.allocated(self%pool(iobs)%Xgrid)) THEN
        allocate ( self%pool(iobs)%Xgrid(mobs) )
        self%pool(iobs)%Xgrid = 0.0_r8
      END IF

      IF (.not.allocated(self%pool(iobs)%Ygrid)) THEN
        allocate ( self%pool(iobs)%Ygrid(mobs) )
        self%pool(iobs)%Ygrid = 0.0_r8
      END IF

      IF (.not.allocated(self%pool(iobs)%Zgrid)) THEN
        allocate ( self%pool(iobs)%Zgrid(mobs) )
        self%pool(iobs)%Zgrid = 0.0_r8
      END IF
!
!  Allocate pool vector element structure, Group ObsError.
!
      IF (.not.allocated(self%pool(iobs)%error)) THEN
        allocate ( self%pool(iobs)%error(mobs,mvar) )
        self%pool(iobs)%error = 0.0_r8
      END IF
!
!  Allocate pool vector element structure, Group ObsValue.
!
      IF (.not.allocated(self%pool(iobs)%value)) THEN
        allocate ( self%pool(iobs)%value(mobs,mvar) )
        self%pool(iobs)%value = 0.0_r8
      END IF
!
!  Allocate statistics variables.
!
      IF (.not.allocated(self%pool(iobs)%error_hash)) THEN
        allocate ( self%pool(iobs)%error_hash(mvar) )
        self%pool(iobs)%error_hash=0_i8b
      END IF
      IF (.not.allocated(self%pool(iobs)%value_hash)) THEN
        allocate ( self%pool(iobs)%value_hash(mvar) )
        self%pool(iobs)%value_hash=0_i8b
      END IF

      IF (.not.allocated(self%pool(iobs)%error_avg)) THEN
        allocate ( self%pool(iobs)%error_avg(mvar) )
        self%pool(iobs)%error_avg=0.0_r8
      END IF
      IF (.not.allocated(self%pool(iobs)%value_avg)) THEN
        allocate ( self%pool(iobs)%value_avg(mvar) )
        self%pool(iobs)%value_avg=0.0_r8
      END IF

      IF (.not.allocated(self%pool(iobs)%error_min)) THEN
        allocate ( self%pool(iobs)%error_min(mvar) )
        self%pool(iobs)%error_min=0.0_r8
      END IF
      IF (.not.allocated(self%pool(iobs)%value_min)) THEN
        allocate ( self%pool(iobs)%value_min(mvar) )
        self%pool(iobs)%value_min=0.0_r8
      END IF

      IF (.not.allocated(self%pool(iobs)%error_max)) THEN
        allocate ( self%pool(iobs)%error_max(mvar) )
        self%pool(iobs)%error_max=0.0_r8
      END IF
      IF (.not.allocated(self%pool(iobs)%value_max)) THEN
        allocate ( self%pool(iobs)%value_max(mvar) )
        self%pool(iobs)%value_max=0.0_r8
      END IF
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
!
!  It allocates the observation arrays in the appropiate pool vector
!  element, and then reads the specified native observation file to
!  populate the data.
!
      SUBROUTINE obs_pool_load_native (self, ng, model, mobs, msurvey,  &
     &                                 IOtype, ncname)
!
      CLASS (obs_pool),  intent(inout) :: self     ! observation set
      integer,           intent(in   ) :: ng       ! nested grid number
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
        IF (id.eq.45) THEN
          mysize=mvar(4)
          mymvar=2
        ELSE IF (id.eq.60) THEN
          mysize=mSST
          mymvar=1
        ELSE IF (id.eq.70) THEN
          mysize=mSSS
          mymvar=1
        ELSE
          mysize=mvar(id)
          mymvar=1
        END IF
!
        self%pool(iobs)%nlocs   = mysize
        self%pool(iobs)%nvars   = mymvar
        self%pool(iobs)%poolID  = self%stateID(iobs)
!
        IF (.not.allocated(self%pool(iobs)%ncname)) THEN
          lstr=LEN_TRIM(ncname)
          allocate ( character(LEN=lstr) :: self%pool(iobs)%ncname )
          self%pool(iobs)%ncname = TRIM(ncname)
        END IF
!
        IF (.not.allocated(self%pool(iobs)%provenance)) THEN
          allocate ( self%pool(iobs)%provenance(mysize) )
          self%pool(iobs)%provenance = 0
        END IF

        IF (.not.allocated(self%pool(iobs)%lon)) THEN
          allocate ( self%pool(iobs)%lon(mysize) )
          self%pool(iobs)%lon = 0.0_r8
        END IF

        IF (.not.allocated(self%pool(iobs)%lat)) THEN
          allocate ( self%pool(iobs)%lat(mysize) )
          self%pool(iobs)%lat = 0.0_r8
        END IF

        IF (.not.allocated(self%pool(iobs)%depth)) THEN
          allocate ( self%pool(iobs)%depth(mysize) )
          self%pool(iobs)%depth = 0.0_r8
        END IF

        IF (.not.allocated(self%pool(iobs)%time)) THEN
          allocate ( self%pool(iobs)%time(mysize) )
          self%pool(iobs)%time = 0.0_r8
        END IF

        IF (.not.allocated(self%pool(iobs)%stateID)) THEN
          allocate ( self%pool(iobs)%stateID(mymvar) )
          self%pool(iobs)%stateID = 0
        END IF

        IF (.not.allocated(self%pool(iobs)%Xgrid)) THEN
          allocate ( self%pool(iobs)%Xgrid(mysize) )
          self%pool(iobs)%Xgrid = 0.0_r8
        END IF

        IF (.not.allocated(self%pool(iobs)%Ygrid)) THEN
          allocate ( self%pool(iobs)%Ygrid(mysize) )
          self%pool(iobs)%Ygrid = 0.0_r8
        END IF

        IF (.not.allocated(self%pool(iobs)%Zgrid)) THEN
          allocate ( self%pool(iobs)%Zgrid(mysize) )
          self%pool(iobs)%Zgrid = 0.0_r8
        END IF

        IF (.not.allocated(self%pool(iobs)%error)) THEN
          allocate ( self%pool(iobs)%error(mysize,mymvar) )
          self%pool(iobs)%error = 0.0_r8
        END IF

        IF (.not.allocated(self%pool(iobs)%value)) THEN
          allocate ( self%pool(iobs)%value(mysize,mymvar) )
          self%pool(iobs)%value = 0.0_r8
        END IF
!
!  Allocate statistics variables.
!
        IF (.not.allocated(self%pool(iobs)%error_hash)) THEN
          allocate ( self%pool(iobs)%error_hash(mymvar) )
          self%pool(iobs)%error_hash=0_i8b
        END IF
        IF (.not.allocated(self%pool(iobs)%value_hash)) THEN
          allocate ( self%pool(iobs)%value_hash(mymvar) )
          self%pool(iobs)%value_hash=0_i8b
        END IF

        IF (.not.allocated(self%pool(iobs)%error_avg)) THEN
          allocate ( self%pool(iobs)%error_avg(mymvar) )
          self%pool(iobs)%error_avg=0.0_r8
        END IF
        IF (.not.allocated(self%pool(iobs)%value_avg)) THEN
          allocate ( self%pool(iobs)%value_avg(mymvar) )
          self%pool(iobs)%value_avg=0.0_r8
        END IF

        IF (.not.allocated(self%pool(iobs)%error_min)) THEN
          allocate ( self%pool(iobs)%error_min(mymvar) )
          self%pool(iobs)%error_min=0.0_r8
        END IF
        IF (.not.allocated(self%pool(iobs)%value_min)) THEN
          allocate ( self%pool(iobs)%value_min(mymvar) )
          self%pool(iobs)%value_min=0.0_r8
        END IF

        IF (.not.allocated(self%pool(iobs)%error_max)) THEN
          allocate ( self%pool(iobs)%error_max(mymvar) )
          self%pool(iobs)%error_max=0.0_r8
        END IF
        IF (.not.allocated(self%pool(iobs)%value_max)) THEN
          allocate ( self%pool(iobs)%value_max(mymvar) )
          self%pool(iobs)%value_max=0.0_r8
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
!
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
      character (len=40)               :: name, vars(2)

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
        IF (.not.allocated(self%pool(iobs)%ncname)) THEN
          lstr=LEN_TRIM(ncname)
          allocate ( character(LEN=lstr) :: self%pool(iobs)%ncname )
          self%pool(iobs)%ncname=TRIM(ncname)
        END IF
!
        SELECT CASE (self%stateID(iobs))
          CASE (1)
            mvar=1
            name='SSH'
            vars(mvar)='absoluteDynamicTopography'
            self%pool(iobs)%stateID(mvar)=1
          CASE (6)
            mvar=1
            name='temp'
            vars(mvar)='waterTemperature'
            self%pool(iobs)%stateID(mvar)=6
          CASE (7)
            mvar=1
            name='salt'
            vars(mvar)='Salinity'
            self%pool(iobs)%stateID(mvar)=7
          CASE (45)
            mvar=2
            name='uv'
            vars(1)='waterSurfaceZonalVelocity'
            vars(2)='waterSurfaceMeridionalVelocity'
            self%pool(iobs)%stateID(1)=4
            self%pool(iobs)%stateID(2)=5
          CASE (60)
            mvar=1
            name='SST'
            vars(mvar)='seaSurfaceTemperature'
            self%pool(iobs)%stateID(1)=6
          CASE (70)
            mvar=1
            name='SSS'
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
!
!  It reads the input observation file(s) to populate the data into
!  the pool object. Each element of the pool vector stores a scalar
!  or vector variable(s) associated with the data assimilation control
!  vector.
!
      SUBROUTINE obs_pool_read (self, ng, model, nfiles)
!
      CLASS (obs_pool), intent(inout) :: self     ! observation set
      integer,          intent(in   ) :: ng       ! nested grid number
      integer,          intent(in   ) :: model    ! kernel identifier
      integer,          intent(in   ) :: nfiles   ! number of files
!
      logical                         :: isNative
      integer                         :: mobs, msurvey, mvar
      integer                         :: IOtype, i, nf
      character (len=256)             :: ncname

      character (len=*),    parameter :: MyFile =                       &
     &    __FILE__//", obs_pool_read"
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
      END SUBROUTINE obs_pool_read
!
!-----------------------------------------------------------------------
!
!  Define the size of the observation pool vector as the number of
!  variables in the state control vector that are directly or
!  indirectly linked to the measurements. Retrieve the state variable
!  stateID (1:nvars) from the input observations file(s), nfiles.
!  Then, allocate observation pool vector.
!
      SUBROUTINE obs_pool_size (self, ng, model, nfiles)
!
      CLASS (obs_pool), intent(inout) :: self     ! observation set
      integer,          intent(in   ) :: ng       ! nested grid number
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
!  If native observation file, read in "obs_type" variable to determine
!  the observed variables and its IDs.
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
     &                                  'obs_Zgrid', obs_Zgris)
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
!
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
          CASE (1)
            mask=(typeID.eq.1)
            mvar=1
          CASE (6)
            mask=(typeID.eq.6).and.(INT(Zgrid).ne.Lsur)
            mvar=1
          CASE (7)
            mask=(typeID.eq.7)
            mvar=1
          CASE (45)                   
            mask=(typeID.eq.4)
            mvar=2
          CASE (60)
            mask=(typeID.eq.6).and.(INT(Zgrid).eq.Lsur)
            mvar=1
          CASE (70)
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
!
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
          CASE (1)
            mask=(typeID.eq.1)
            mvar=1
          CASE (6)
            mask=(typeID.eq.6).and.(INT(Zgrid).ne.Lsur)
            mvar=1
          CASE (7)
            mask=(typeID.eq.7)
            mvar=1
          CASE (45)                   
            mask=(typeID.eq.4)
            msk2=(typeID.eq.5)
            mvar=2
          CASE (60)
            mask=(typeID.eq.6).and.(INT(Zgrid).eq.Lsur)
            mvar=1
          CASE (70)
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
