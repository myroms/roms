#include "cppdefs.h"
      MODULE mod_ncparam
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This MODULE contains all the variables associated with input and    !
!  output  NetCDF  files.  The IO model is very generic and easy to    !
!  change or expand.  The NetCDF files can be in any language.  All    !
!  the IO information is managed using the following variables:        !
!                                                                      !
!  Vname      Input/output variables names and attributes:             !
!               Vname(1,*)  => field variable name.                    !
!               Vname(2,*)  => long-name attribute.                    !
!               Vname(3,*)  => units attribute.                        !
!               Vname(4,*)  => field type attribute.                   !
!               Vname(5,*)  => associated time variable name.          !
!  Tname      Input/output associated time variables names.            !
!                                                                      !
!  Linfo      Input/output fields logical information:                 !
!               Linfo(1,*)  => switch indicating grided data.          !
!               Linfo(2,*)  => switch indicating time cycling.         !
!               Linfo(3,*)  => switch indicating only one-time         !
!                              record available.                       !
!                                                                      !
!  Iinfo      Input/output fields integer information:                 !
!               Iinfo(1,*)  => variable grid type.                     !
!               Iinfo(2,*)  => field variable NetCDF ID.               !
!               Iinfo(3,*)  => associated time variable NetCDF ID.     !
!               Iinfo(4,*)  => number of time records.                 !
!               Iinfo(5,*)  => rolling two-time levels index.          !
!               Iinfo(6,*)  => latest processed time record .          !
!                                                                      !
!  Finfo      Input/output field floating-point information:           !
!               Finfo(1,*)  => starting time (days) of data.           !
!               Finfo(2,*)  => Ending time (days) of data.             !
!               Finfo(3,*)  => length (days) of time cycling.          !
!               Finfo(4,*)  => Scale to convert time to day units.     !
!               Finfo(5,*)  => latest monotonic time (sec).            !
!               Finfo(6,*)  => minimum value for current data.         !
!               Finfo(7,*)  => maximum value for current data.         !
!  Fscale     Scale to convert input data to model units.              !
!  Fpoint     Latest two-time records of input point data.             !
!  Tintrp     Time (sec) of latest field snapshots used for            !
!               interpolation.                                         !
!  Vtime      Latest two-time values of processed input data.          !
!                                                                      !
!=======================================================================
!
        USE mod_param        

        implicit none
!
!  Maximum number of variables in a generic NetCDF file (MV) and
!  maximum number of variables in information arrays (NV).
!
        integer, parameter :: MV = 250
        integer, parameter :: NV = 250
!
!  Input/output grid-type variables.
!
        integer, parameter :: p2dvar = 1         ! 2D PSI-variable
        integer, parameter :: r2dvar = 2         ! 2D RHO-variable
        integer, parameter :: u2dvar = 3         ! 2D U-variable
        integer, parameter :: v2dvar = 4         ! 2D V-variable
        integer, parameter :: p3dvar = 5         ! 3D PSI-variable
        integer, parameter :: r3dvar = 6         ! 3D RHO-variable
        integer, parameter :: u3dvar = 7         ! 3D U-variable
        integer, parameter :: v3dvar = 8         ! 3D V-variable
        integer, parameter :: w3dvar = 9         ! 3D W-variable
        integer, parameter :: b3dvar = 10        ! 3D BED-sediment
!
!  Number of horizontal water points.
!
        integer  :: Nxyp(Ngrids)   ! PSI water points
        integer  :: Nxyr(Ngrids)   ! RHO water points
        integer  :: Nxyu(Ngrids)   ! U water points
        integer  :: Nxyv(Ngrids)   ! V water points
!
!  Switches indicating which variables are written to output files.
!
        logical  :: Hout(NV,Ngrids) ! history file switches
        logical  :: Sout(NV,Ngrids) ! station file switches
!
!  Input/output identification indices.
!
        integer  :: idXgrd = -1    ! XI-grid position
        integer  :: idYgrd = -2    ! ETA-grid position
        integer  :: idZgrd = -3    ! S-grid position
        integer  :: iddpth = -4    ! depth
        integer  :: idglon = -5    ! longitude
        integer  :: idglat = -6    ! latitude

        integer  :: idAbed         ! bed wave excursion amplitude
        integer  :: idAice         ! fraction of cell covered by ice
        integer  :: idCfra         ! cloud fraction
        integer  :: idDano         ! density anomaly
        integer  :: iddQdT         ! heat flux sensitivity to SST
        integer  :: idevap         ! evaporation rate
        integer  :: idFsur         ! free-surface
        integer  :: idHbbl         ! depth of bottom boundary layer
        integer  :: idHice         ! depth of ice cover
        integer  :: idHrip         ! bed ripple height
        integer  :: idHsbl         ! depth of surface boundary layer
        integer  :: idHsno         ! depth of snow cover
        integer  :: idLhea         ! net latent heat flux
        integer  :: idLrad         ! net longwave radiation flux
        integer  :: idLrip         ! bed ripple length
        integer  :: idMtke         ! turbulent kinetic energy
        integer  :: idMtls         ! turbulent length scale
        integer  :: idOvel         ! omega vertical velocity
        integer  :: idQair         ! surface air humidity
        integer  :: idPair         ! surface air pressure
        integer  :: idPbar         ! streamfunction
        integer  :: idRdir         ! river runoff direction
        integer  :: idRepo         ! river runoff ETA-positions
        integer  :: idRflg         ! river runoff flag
        integer  :: idRtra         ! river runoff mass transport
        integer  :: idRu2d         ! RHS of 2D U-momentum
        integer  :: idRv2d         ! RHS of 2D V-momentum
        integer  :: idRxpo         ! river runoff XI-positions
        integer  :: idRvsh         ! river runoff transport profile
        integer  :: idRzet         ! RHS of free-surface
        integer  :: idrain         ! rainfall rate
        integer  :: idSden         ! bottom sediment grain density
        integer  :: idSdif         ! vertical S-diffusion coefficient
        integer  :: idSrad         ! net shortwave radiation flux
        integer  :: idSSHc         ! SSH climatology
        integer  :: idSSHe         ! SSH error variance
        integer  :: idSSHo         ! SSH observations
        integer  :: idSSSc         ! SSS climatology
        integer  :: idSSTc         ! SST climatology
        integer  :: idSSTe         ! SST error variance
        integer  :: idSSTo         ! SST observations
        integer  :: idShea         ! net sensible heat flux
        integer  :: idSsiz         ! bottom sediment grain size
        integer  :: idTair         ! surface air temperature
        integer  :: idTdif         ! vertical T-diffusion coefficient
        integer  :: idTice         ! temperature of ice surface
        integer  :: idtime         ! ocean time
        integer  :: idTper         ! tidal period
        integer  :: idTvan         ! tidal current angle
        integer  :: idTvma         ! maximum tidal current
        integer  :: idTvmi         ! minimum tidal current
        integer  :: idTvph         ! tidal current phase
        integer  :: idTzam         ! tidal elevation amplitude
        integer  :: idTzph         ! tidal elevation phase
        integer  :: idUair         ! surface U-wind
        integer  :: idUbar         ! 2D U-velocity
        integer  :: idUbcl         ! 2D U-velocity climatology
        integer  :: idUbed         ! bed wave orbital U-velocity
        integer  :: idUbms         ! bottom U-momentum stress
        integer  :: idUbot         ! bottom U-momentum above bed
        integer  :: idUbws         ! bottom U-wave stress
        integer  :: idUclm         ! 3D U-velocity climatology
        integer  :: idUice         ! ice U-velocity
        integer  :: idUobs         ! 3D U-velocity observations
        integer  :: idUsms         ! surface U-momentum stress
        integer  :: idUsur         ! surface U-velocity observations
        integer  :: idUVer         ! 3D velocity error variance
        integer  :: idUVse         ! surface velocity error variance
        integer  :: idUvel         ! 3D U-velocity
        integer  :: idVair         ! surface V-wind
        integer  :: idVbar         ! 2D V-velocity
        integer  :: idVbcl         ! 2D V-velocity climatology
        integer  :: idVbed         ! bed wave orbital V-velocity
        integer  :: idVbms         ! bottom V-momentum stress
        integer  :: idVbot         ! bottom V-momentum above bed
        integer  :: idVbws         ! bottom V-wave stress
        integer  :: idVclm         ! 3D V-velocity climatology
        integer  :: idVice         ! ice V-velocity
        integer  :: idVobs         ! 3D V-velocity observations
        integer  :: idVsms         ! surface V-momentum stress
        integer  :: idVsur         ! surface V-velocity observations
        integer  :: idVvel         ! 3D V-velocity
        integer  :: idVvis         ! vertical viscosity coefficient
        integer  :: idWamp         ! wind-induced wave amplitude
        integer  :: idWdir         ! wind-induced wave direction
        integer  :: idWper         ! wind-induced wave period
        integer  :: idWvel         ! true vertical velocity
        integer  :: idZapp         ! apparent hydraulic bottom roughness
        integer  :: idZnot         ! hydraulic bottom roughness

        integer, allocatable :: idRtrc(:)     ! river runoff for tracers
        integer, allocatable :: idTbot(:)     ! bottom flux for tracers
        integer, allocatable :: idTbry(:,:)   ! tracers boundary conditions
        integer, allocatable :: idTclm(:)     ! tracers climatology
        integer, allocatable :: idTerr(:)     ! tracers error variance
        integer, allocatable :: idTobs(:)     ! tracers observations
        integer, allocatable :: idTsur(:)     ! surface flux for tracers

        integer  :: idU2bc(4)      ! 2D U-velocity boundary conditions
        integer  :: idU3bc(4)      ! 3D U-velocity boundary conditions
        integer  :: idV2bc(4)      ! 2D V-velocity boundary conditions
        integer  :: idV3bc(4)      ! 3D V-velocity boundary conditions
        integer  :: idZbry(4)      ! free-surface boundary conditions
!
!  Assimilation state variables indices (order is important).
!
        integer  :: isFsur = 1                ! free-surface
        integer  :: isUbar = 2                ! 2D U-velocity
        integer  :: isVbar = 3                ! 2D V-velocity
        integer  :: isUvel = 4                ! 3D U-velocity
        integer  :: isVvel = 5                ! 3D V-velocity
        integer, allocatable :: isTvar(:)     ! tracers
        integer, allocatable :: idSvar(:)     ! full state vector indices
!
!  Input/Output NetCDF files IDs.
!
        integer  :: ncSSHid(Ngrids)     ! SSH observations
        integer  :: ncSSTid(Ngrids)     ! SST observations
        integer  :: ncTOBSid(Ngrids)    ! tracer observations
        integer  :: ncVOBSid(Ngrids)    ! currents observations
        integer  :: ncVSURid(Ngrids)    ! surface currents observations
        integer  :: ncavgid(Ngrids)     ! output averages
        integer  :: ncbryid(Ngrids)     ! input boundary conditions
        integer  :: ncclmid(Ngrids)     ! input climatology
        integer  :: ncesid(Ngrids)      ! error subspace
        integer  :: ncfltid(Ngrids)     ! output floats
        integer  :: ncfrcid(NV,Ngrids)  ! input forcing
        integer  :: ncfwdid(Ngrids)     ! forward solution
        integer  :: ncgrdid(Ngrids)     ! input grid
        integer  :: nchisid(Ngrids)     ! output history
        integer  :: ncrstid(Ngrids)     ! input/output restart
        integer  :: ncstaid(Ngrids)     ! output stations

        integer  :: idefavg(Ngrids)     ! averages file creation flag
        integer  :: idefhis(Ngrids)     ! history file creation flag
!
!  Output NetCDF variables IDs.
!
        integer, allocatable :: idTvar(:)     ! tracers variables
        integer, allocatable :: avgTid(:,:)   ! averages tracers IDs
        integer, allocatable :: fltTid(:,:)   ! floats tracers IDs
        integer, allocatable :: hisTid(:,:)   ! history tracers IDs
        integer, allocatable :: rstTid(:,:)   ! restart tracers IDs
        integer, allocatable :: staTid(:,:)   ! stations tracers IDs

        integer  :: avgVid(NV,Ngrids)    ! averages variables IDs
        integer  :: fltVid(-6:NV,Ngrids) ! floats variables IDs
        integer  :: hisVid(NV,Ngrids)    ! history variables IDs
        integer  :: rstVid(NV,Ngrids)    ! restart variables IDs
        integer  :: staVid(NV,Ngrids)    ! stations variables IDs

        integer  :: tavgindx(Ngrids)     ! averages time record index
        integer  :: tfltindx(Ngrids)     ! floats time record index
        integer  :: thisindx(Ngrids)     ! history time record index
        integer  :: trstindx(Ngrids)     ! restart time record index
        integer  :: tstaindx(Ngrids)     ! stations time record index
!
!  Current NetCDF file information.
!
        integer  :: nvars          ! number of variables
        integer  :: nvdims(MV)     ! number of variable dimensions
        integer  :: tsize          ! size of unlimited dimension
        integer  :: vdims(5,MV)    ! variables dimensions IDs
        integer  :: vflag(MV)      ! variable flag: 1: full field
                                   !               -1: water points only
        character (len=80 ) :: type
!
!  Input/Output information variables.
!
        logical  :: Linfo(3,NV,Ngrids)

        integer  :: Iinfo(6,NV,Ngrids)

        real(r8) :: Finfo(7,NV,Ngrids)
        real(r8) :: Fpoint(2,NV,Ngrids)
        real(r8) :: Fscale(NV,Ngrids)
        real(r8) :: Tintrp(2,NV,Ngrids)
        real(r8) :: Vtime(2,NV,Ngrids)

        character (len=5  ) :: version = '2.0  '
        character (len=15 ) :: varnam(MV)
        character (len=44 ) :: date_str
        character (len=46 ) :: Vname(5,0:NV)
        character (len=46 ) :: Tname(0:NV)
        character (len=120) :: history

      CONTAINS

      SUBROUTINE initialize_ncparam
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine allocates and initializes all variables in module      !
!  "mod_ncparam" for all nested grids.                                 !
!                                                                      !
!=======================================================================
!
        USE mod_parallel
        USE mod_iounits
        USE mod_scalars
#ifdef SEDIMENT
        USE mod_sediment
#endif
!
!  Local variable declarations.
!
        logical :: load

        integer, parameter :: inp = 10
        integer :: Lvar, i, gtype, ng, varid

        real(r8) :: scale

        character (len=80), dimension(7) :: Vinfo
!
!-----------------------------------------------------------------------
!  Allocate tracer type variables.
!-----------------------------------------------------------------------
!
        allocate ( idRtrc(MT) )
        allocate ( idTbot(MT) )
        allocate ( idTbry(4,MT) )
        allocate ( idTclm(MT) )
        allocate ( idTerr(MT) )
        allocate ( idTobs(MT) )
        allocate ( idTsur(MT) )
        allocate ( idTvar(MT) )
        allocate ( isTvar(MT) )
      
        allocate ( avgTid(MT,Ngrids) )
        allocate ( fltTid(MT,Ngrids) )
        allocate ( hisTid(MT,Ngrids) )
        allocate ( rstTid(MT,Ngrids) )
        allocate ( staTid(MT,Ngrids) )

        allocate ( idSvar(MAXVAL(NSV)) )
!
!-----------------------------------------------------------------------
!  Initialize several variables.
!-----------------------------------------------------------------------
!
!  Initialize NetCDF files IDs to closed state.
!
        DO ng=1,Ngrids
          DO i=1,NV
            ncfrcid(i,ng)=-1
          END DO
          ncSSHid(ng)=-1
          ncSSTid(ng)=-1
          ncTOBSid(ng)=-1
          ncVOBSid(ng)=-1
          ncVSURid(ng)=-1
          ncavgid(ng)=-1
          ncbryid(ng)=-1
          ncclmid(ng)=-1
          ncesid(ng)=-1
          ncfltid(ng)=-1
          ncfwdid(ng)=-1
          ncgrdid(ng)=-1
          nchisid(ng)=-1
          ncrstid(ng)=-1
          ncstaid(ng)=-1

          tavgindx(ng)=0
          tfltindx(ng)=0
          thisindx(ng)=0
          trstindx(ng)=0
          tstaindx(ng)=0

          idefavg(ng)=-1
          idefhis(ng)=-1
        END DO
!
!  Set IDs for tracers state variables.
!
        DO i=1,MT
          isTvar(i)=5+i
        END DO
!
!  Initialize IO information variables.
!
        DO ng=1,Ngrids
          DO i=1,NV
            Linfo(1,i,ng)=.false.
            Linfo(2,i,ng)=.false.
            Linfo(3,i,ng)=.false.
            Hout(i,ng)=.false.
            Sout(i,ng)=.false.
            Iinfo(1,i,ng)=0
            Iinfo(2,i,ng)=-1
            Iinfo(3,i,ng)=-1
            Iinfo(4,i,ng)=0
            Iinfo(5,i,ng)=0
            Iinfo(6,i,ng)=0
            Finfo(1,i,ng)=0.0_r8
            Finfo(2,i,ng)=0.0_r8
            Finfo(3,i,ng)=0.0_r8
            Finfo(5,i,ng)=0.0_r8
            Finfo(6,i,ng)=0.0_r8
            Finfo(7,i,ng)=0.0_r8
            Fscale(i,ng)=1.0_r8
            Fpoint(1,i,ng)=0.0_r8
            Fpoint(2,i,ng)=0.0_r8
            Tintrp(1,i,ng)=0.0_r8
            Tintrp(2,i,ng)=0.0_r8
            Vtime(1,i,ng)=0.0_r8
            Vtime(2,i,ng)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Define names of variables for Input/Output NetCDF files.
!-----------------------------------------------------------------------
!
!  Open input variable information file.
!
        OPEN (inp, FILE=TRIM(varname), FORM='formatted', STATUS='old',  &
     &        ERR=10)
        GOTO 20
  10    IF (Master) WRITE(stdout,50) TRIM(varname)
        STOP
  20    CONTINUE
!
!  Read in variable information.  Ignore blank and comment [char(33)=!]
!  input lines.
!
        varid=0
        DO WHILE (.true.)
          READ (inp,*,ERR=30,END=40) Vinfo(1)
          Lvar=LEN_TRIM(Vinfo(1))
          IF ((Lvar.gt.0).and.(Vinfo(1)(1:1).ne.CHAR(33))) THEN
            READ (inp,*,ERR=30) Vinfo(2)
            READ (inp,*,ERR=30) Vinfo(3)
            READ (inp,*,ERR=30) Vinfo(4)
            READ (inp,*,ERR=30) Vinfo(5)
            READ (inp,*,ERR=30) Vinfo(6)
            READ (inp,*,ERR=30) Vinfo(7)
            READ (inp,*,ERR=30) scale
!
!  Determine staggered C-grid variable.
!
            SELECT CASE (TRIM(ADJUSTL(Vinfo(7))))
              CASE ('p2dvar')
                gtype=p2dvar
              CASE ('r2dvar')
                gtype=r2dvar
              CASE ('u2dvar')
                gtype=u2dvar
              CASE ('v2dvar')
                gtype=v2dvar
              CASE ('p3dvar')
                gtype=p3dvar
              CASE ('r3dvar')
                gtype=r3dvar
              CASE ('u3dvar')
                gtype=u3dvar
              CASE ('v3dvar')
                gtype=v3dvar
              CASE ('w3dvar')
                gtype=w3dvar
              CASE ('b3dvar')
                gtype=b3dvar
              CASE DEFAULT
                gtype=0
            END SELECT
!
!  Assign identification indices.
!
            load=.true.
            varid=varid+1
            SELECT CASE (TRIM(ADJUSTL(Vinfo(6))))
              CASE ('idtime')
                idtime=varid
              CASE ('idFsur')
                idFsur=varid
              CASE ('idUbar')
                idUbar=varid
              CASE ('idVbar')
                idVbar=varid
              CASE ('idUvel')
                idUvel=varid
              CASE ('idVvel')
                idVvel=varid
              CASE ('idWvel')
                idWvel=varid
              CASE ('idOvel')
                idOvel=varid
              CASE ('idDano')
                idDano=varid
              CASE ('idTvar(itemp)')
                idTvar(itemp)=varid
              CASE ('idTvar(isalt)')
                idTvar(isalt)=varid
              CASE ('idUsms')
                idUsms=varid
              CASE ('idVsms')
                idVsms=varid
              CASE ('idUbms')
                idUbms=varid
              CASE ('idVbms')
                idVbms=varid
              CASE ('idUbws')
                idUbws=varid
              CASE ('idVbws')
                idVbws=varid
              CASE ('idTsur(itemp)')
                idTsur(itemp)=varid
              CASE ('iddQdT')
                iddQdT=varid
              CASE ('idTsur(isalt)')
                idTsur(isalt)=varid
              CASE ('idTbot(itemp)')
                idTbot(itemp)=varid
              CASE ('idTbot(isalt)')
                idTbot(isalt)=varid
              CASE ('idMtke')
                idMtke=varid
              CASE ('idMtls')
                idMtls=varid
              CASE ('idVvis')
                idVvis=varid
              CASE ('idTdif')
                idTdif=varid
              CASE ('idSdif')
                idSdif=varid
              CASE ('idZbry(iwest)')
                idZbry(iwest)=varid
              CASE ('idZbry(ieast)')
                idZbry(ieast)=varid
              CASE ('idZbry(isouth)')
                idZbry(isouth)=varid
              CASE ('idZbry(inorth)')
                idZbry(inorth)=varid
              CASE ('idU2bc(iwest)')
                idU2bc(iwest)=varid
              CASE ('idU2bc(ieast)')
                idU2bc(ieast)=varid
              CASE ('idU2bc(isouth)')
                idU2bc(isouth)=varid
              CASE ('idU2bc(inorth)')
                idU2bc(inorth)=varid
              CASE ('idV2bc(iwest)')
                idV2bc(iwest)=varid
              CASE ('idV2bc(ieast)')
                idV2bc(ieast)=varid
              CASE ('idV2bc(isouth)')
                idV2bc(isouth)=varid
              CASE ('idV2bc(inorth)')
                idV2bc(inorth)=varid
              CASE ('idU3bc(iwest)')
                idU3bc(iwest)=varid
              CASE ('idU3bc(ieast)')
                idU3bc(ieast)=varid
              CASE ('idU3bc(isouth)')
                idU3bc(isouth)=varid
              CASE ('idU3bc(inorth)')
                idU3bc(inorth)=varid
              CASE ('idV3bc(iwest)')
                idV3bc(iwest)=varid
              CASE ('idV3bc(ieast)')
                idV3bc(ieast)=varid
              CASE ('idV3bc(isouth)')
                idV3bc(isouth)=varid
              CASE ('idV3bc(inorth)')
                idV3bc(inorth)=varid
              CASE ('idTbry(iwest,itemp)')
                idTbry(iwest,itemp)=varid
              CASE ('idTbry(ieast,itemp)')
                idTbry(ieast,itemp)=varid
              CASE ('idTbry(isouth,itemp)')
                idTbry(isouth,itemp)=varid
              CASE ('idTbry(inorth,itemp)')
                idTbry(inorth,itemp)=varid
              CASE ('idTbry(iwest,isalt)')
                idTbry(iwest,isalt)=varid
              CASE ('idTbry(ieast,isalt)')
                idTbry(ieast,isalt)=varid
              CASE ('idTbry(isouth,isalt)')
                idTbry(isouth,isalt)=varid
              CASE ('idTbry(inorth,isalt)')
                idTbry(inorth,isalt)=varid
              CASE ('idPair')
                idPair=varid
              CASE ('idTair')
                idTair=varid
              CASE ('idQair')
                idQair=varid
              CASE ('idCfra')
                idCfra=varid
              CASE ('idSrad')
                idSrad=varid
              CASE ('idLrad')
                idLrad=varid
              CASE ('idLhea')
                idLhea=varid
              CASE ('idShea')
                idShea=varid
              CASE ('idrain')
                idrain=varid
              CASE ('idevap')
                idevap=varid
              CASE ('idUair')
                idUair=varid
              CASE ('idVair')
                idVair=varid
              CASE ('idWamp')
                idWamp=varid
              CASE ('idWdir')
                idWdir=varid
              CASE ('idWper')
                idWper=varid
              CASE ('idTper')
                idTper=varid
              CASE ('idTzam')
                idTzam=varid
              CASE ('idTzph')
                idTzph=varid
              CASE ('idTvph')
                idTvph=varid
              CASE ('idTvan')
                idTvan=varid
              CASE ('idTvma')
                idTvma=varid
              CASE ('idTvmi')
                idTvmi=varid
              CASE ('idRxpo')
                idRxpo=varid
              CASE ('idRepo')
                idRepo=varid
              CASE ('idRdir')
                idRdir=varid
              CASE ('idRvsh')
                idRvsh=varid
              CASE ('idRtra')
                idRtra=varid
              CASE ('idRflg')
                idRflg=varid
              CASE ('idRtrc(itemp)')
                idRtrc(itemp)=varid
              CASE ('idRtrc(isalt)')
                idRtrc(isalt)=varid
              CASE ('idHsbl')
                idHsbl=varid
              CASE ('idHbbl')
                idHbbl=varid
              CASE ('idUbot')
                idUbot=varid
              CASE ('idVbot')
                idVbot=varid
              CASE ('idUbed')
                idUbed=varid
              CASE ('idVbed')
                idVbed=varid
              CASE ('idAbed')
                idAbed=varid
              CASE ('idZnot')
                idZnot=varid
              CASE ('idZapp')
                idZapp=varid
              CASE ('idHrip')
                idHrip=varid
              CASE ('idLrip')
                idLrip=varid
              CASE ('idSSHc')
                idSSHc=varid
              CASE ('idUbcl')
                idUbcl=varid
              CASE ('idVbcl')
                idVbcl=varid
              CASE ('idUclm')
                idUclm=varid
              CASE ('idVclm')
                idVclm=varid
              CASE ('idSSSc')
                idSSSc=varid
              CASE ('idSSTc')
                idSSTc=varid
              CASE ('idTclm(itemp)')
                idTclm(itemp)=varid
              CASE ('idTclm(isalt)')
                idTclm(isalt)=varid
              CASE ('idSSHo')
                idSSHo=varid
              CASE ('idSSHe')
                idSSHe=varid
              CASE ('idUobs')
                idUobs=varid
              CASE ('idVobs')
                idVobs=varid
              CASE ('idUVer')
                idUVer=varid
              CASE ('idUsur')
                idUsur=varid
              CASE ('idVsur')
                idVsur=varid
              CASE ('idUVse')
                idUVse=varid
              CASE ('idSSTo')
                idSSTo=varid
              CASE ('idSSTe')
                idSSTe=varid
              CASE ('idTobs(itemp)')
                idTobs(itemp)=varid
              CASE ('idTerr(itemp)')
                idTerr(itemp)=varid
              CASE ('idTobs(isalt)')
                idTobs(isalt)=varid
              CASE ('idTerr(isalt)')
                idTerr(isalt)=varid
              CASE ('idSvar(isFsur)')
                idSvar(isFsur)=varid
              CASE ('idSvar(isUbar)')
                idSvar(isUbar)=varid
              CASE ('idSvar(isVbar)')
                idSvar(isVbar)=varid
#ifdef SOLVE3D
              CASE ('idSvar(isUvel)')
                idSvar(isUvel)=varid
              CASE ('idSvar(isVvel)')
                idSvar(isVvel)=varid
              CASE ('idSvar(isTvar(itemp))')
                idSvar(isTvar(itemp))=varid
# ifdef SALINITY
              CASE ('idSvar(isTvar(isalt))')
                idSvar(isTvar(isalt))=varid
# endif
#endif
              CASE ('idSsiz')
                idSsiz=varid
              CASE ('idSden')
                idSden=varid
#ifdef SEDIMENT
              CASE ('idTvar(idmud(i))')
                load=.true.
              CASE ('idTvar(isand(i))')
                load=.true.                
              CASE ('idfrac')
                load=.true.
              CASE ('idSbed(ithck)')
                idSbed(ithck)=varid
              CASE ('idSbed(iaged)')
                idSbed(iaged)=varid
              CASE ('idSbed(iporo)')
                idSbed(iporo)=varid
              CASE ('idSbed(idiff)')
                idSbed(idiff)=varid
              CASE ('idBott(irlen)')
                idBott(irlen)=varid
              CASE ('idBott(irhgt)')
                idBott(irhgt)=varid
              CASE ('idBott(ishgt)')
                idBott(ishgt)=varid
              CASE ('idBott(izNik)')
                idBott(izNik)=varid
              CASE ('idBott(izbio)')
                idBott(izbio)=varid
              CASE ('idBott(izbld)')
                idBott(izbld)=varid
              CASE ('idBott(izbfm)')
                idBott(izbfm)=varid
              CASE ('idBott(izapp)')
                idBott(izapp)=varid
#endif
#ifdef BIO_FASHAM
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iPhyt)')
                idTvar(iPhyt)=varid
              CASE ('idTvar(iZoop)')
                idTvar(iZoop)=varid
              CASE ('idTvar(iLDeN)')
                idTvar(iLDeN)=varid
              CASE ('idTvar(iSDeN)')
                idTvar(iSDeN)=varid
              CASE ('idTvar(iChlo)')
                idTvar(iChlo)=varid
              CASE ('idTbry(iwest,iNO3_)')
                idTbry(iwest,iNO3_)=varid
              CASE ('idTbry(ieast,iNO3_)')
                idTbry(ieast,iNO3_)=varid
              CASE ('idTbry(isouth,iNO3_)')
                idTbry(isouth,iNO3_)=varid
              CASE ('idTbry(inorth,iNO3_)')
                idTbry(inorth,iNO3_)=varid
              CASE ('idTbry(iwest,iNH4_)')
                idTbry(iwest,iNH4_)=varid
              CASE ('idTbry(ieast,iNH4_)')
                idTbry(ieast,iNH4_)=varid
              CASE ('idTbry(isouth,iNH4_)')
                idTbry(isouth,iNH4_)=varid
              CASE ('idTbry(inorth,iNH4_)')
                idTbry(inorth,iNH4_)=varid
              CASE ('idTbry(iwest,iPhyt)')
                idTbry(iwest,iPhyt)=varid
              CASE ('idTbry(ieast,iPhyt)')
                idTbry(ieast,iPhyt)=varid
              CASE ('idTbry(isouth,iPhyt)')
                idTbry(isouth,iPhyt)=varid
              CASE ('idTbry(inorth,iPhyt)')
                idTbry(inorth,iPhyt)=varid
              CASE ('idTbry(iwest,iZoop)')
                idTbry(iwest,iZoop)=varid
              CASE ('idTbry(ieast,iZoop)')
                idTbry(ieast,iZoop)=varid
              CASE ('idTbry(isouth,iZoop)')
                idTbry(isouth,iZoop)=varid
              CASE ('idTbry(inorth,iZoop)')
                idTbry(inorth,iZoop)=varid
              CASE ('idTbry(iwest,iSDeN)')
                idTbry(iwest,iSDeN)=varid
              CASE ('idTbry(ieast,iSDeN)')
                idTbry(ieast,iSDeN)=varid
              CASE ('idTbry(isouth,iSDeN)')
                idTbry(isouth,iSDeN)=varid
              CASE ('idTbry(inorth,iSDeN)')
                idTbry(inorth,iSDeN)=varid
              CASE ('idTbry(iwest,iLDeN)')
                idTbry(iwest,iLDeN)=varid
              CASE ('idTbry(ieast,iLDeN)')
                idTbry(ieast,iLDeN)=varid
              CASE ('idTbry(isouth,iLDeN)')
                idTbry(isouth,iLDeN)=varid
              CASE ('idTbry(inorth,iLDeN)')
                idTbry(inorth,iLDeN)=varid
              CASE ('idTbry(iwest,iChlo)')
                idTbry(iwest,iChlo)=varid
              CASE ('idTbry(ieast,iChlo)')
                idTbry(ieast,iChlo)=varid
              CASE ('idTbry(isouth,iChlo)')
                idTbry(isouth,iChlo)=varid
              CASE ('idTbry(inorth,iChlo)')
                idTbry(inorth,iChlo)=varid
#endif
#ifdef ECOSIM
              CASE ('idTvar(iNO3_)')
                idTvar(iNO3_)=varid
              CASE ('idTvar(iNH4_)')
                idTvar(iNH4_)=varid
              CASE ('idTvar(iPO4_)')
                idTvar(iPO4_)=varid
              CASE ('idTvar(iFeO_)')
                idTvar(iFeO_)=varid
              CASE ('idTvar(iSiO_)')
                idTvar(iSiO_)=varid
              CASE ('idTvar(iDIC_)')
                idTvar(iDIC_)=varid
              CASE ('idTvar(iBacC(1))')
                idTvar(iBacC(1))=varid
              CASE ('idTvar(iBacF(1))')
                idTvar(iBacF(1))=varid
              CASE ('idTvar(iBacN(1))')
                idTvar(iBacN(1))=varid
              CASE ('idTvar(iBacP(1))')
                idTvar(iBacP(1))=varid
              CASE ('idTvar(iCDMC(1))')
                idTvar(iCDMC(1))=varid
              CASE ('idTvar(iCDMC(2))')
                idTvar(iCDMC(2))=varid
              CASE ('idTvar(iDOMC(1))')
                idTvar(iDOMC(1))=varid
              CASE ('idTvar(iDOMC(2))')
                idTvar(iDOMC(2))=varid
              CASE ('idTvar(iDOMN(1))')
                idTvar(iDOMN(1))=varid
              CASE ('idTvar(iDOMN(2))')
                idTvar(iDOMN(2))=varid
              CASE ('idTvar(iDOMP(1))')
                idTvar(iDOMP(1))=varid
              CASE ('idTvar(iDOMP(2))')
                idTvar(iDOMP(2))=varid
              CASE ('idTvar(iFecC(1))')
                idTvar(iFecC(1))=varid
              CASE ('idTvar(iFecC(2))')
                idTvar(iFecC(2))=varid
              CASE ('idTvar(iFecF(1))')
                idTvar(iFecF(1))=varid
              CASE ('idTvar(iFecF(2))')
                idTvar(iFecF(2))=varid
              CASE ('idTvar(iFecN(1))')
                idTvar(iFecN(1))=varid
              CASE ('idTvar(iFecN(2))')
                idTvar(iFecN(2))=varid
              CASE ('idTvar(iFecP(1))')
                idTvar(iFecP(1))=varid
              CASE ('idTvar(iFecP(2))')
                idTvar(iFecP(2))=varid
              CASE ('idTvar(iFecS(1))')
                idTvar(iFecS(1))=varid
              CASE ('idTvar(iFecS(2))')
                idTvar(iFecS(2))=varid
              CASE ('idTvar(iPhyC(1))')
                idTvar(iPhyC(1))=varid
              CASE ('idTvar(iPhyC(2))')
                idTvar(iPhyC(2))=varid
              CASE ('idTvar(iPhyC(3))')
                idTvar(iPhyC(3))=varid
              CASE ('idTvar(iPhyC(4))')
                idTvar(iPhyC(4))=varid
              CASE ('idTvar(iPhyF(1))')
                idTvar(iPhyF(1))=varid
              CASE ('idTvar(iPhyF(2))')
                idTvar(iPhyF(2))=varid
              CASE ('idTvar(iPhyF(3))')
                idTvar(iPhyF(3))=varid
              CASE ('idTvar(iPhyF(4))')
                idTvar(iPhyF(4))=varid
              CASE ('idTvar(iPhyN(1))')
                idTvar(iPhyN(1))=varid
              CASE ('idTvar(iPhyN(2))')
                idTvar(iPhyN(2))=varid
              CASE ('idTvar(iPhyN(3))')
                idTvar(iPhyN(3))=varid
              CASE ('idTvar(iPhyN(4))')
                idTvar(iPhyN(4))=varid
              CASE ('idTvar(iPhyP(1))')
                idTvar(iPhyP(1))=varid
              CASE ('idTvar(iPhyP(2))')
                idTvar(iPhyP(2))=varid
              CASE ('idTvar(iPhyP(3))')
                idTvar(iPhyP(3))=varid
              CASE ('idTvar(iPhyP(4))')
                idTvar(iPhyP(4))=varid
              CASE ('idTvar(iPhyS(1))')
                idTvar(iPhyS(1))=varid
              CASE ('idTvar(iPhyS(2))')
                idTvar(iPhyS(2))=varid
              CASE ('idTvar(iPigs(1,1))')
                idTvar(iPigs(1,1))=varid
              CASE ('idTvar(iPigs(2,1))')
                idTvar(iPigs(2,1))=varid
              CASE ('idTvar(iPigs(3,1))')
                idTvar(iPigs(3,1))=varid
              CASE ('idTvar(iPigs(4,1))')
                idTvar(iPigs(4,1))=varid
              CASE ('idTvar(iPigs(1,3))')
                idTvar(iPigs(1,3))=varid
              CASE ('idTvar(iPigs(2,3))')
                idTvar(iPigs(2,3))=varid
              CASE ('idTvar(iPigs(3,3))')
                idTvar(iPigs(3,3))=varid
              CASE ('idTvar(iPigs(1,4))')
                idTvar(iPigs(1,4))=varid
              CASE ('idTvar(iPigs(2,4))')
                idTvar(iPigs(2,4))=varid
              CASE ('idTvar(iPigs(3,4))')
                idTvar(iPigs(3,4))=varid
              CASE ('idTvar(iPigs(1,5))')
                idTvar(iPigs(1,5))=varid
              CASE ('idTvar(iPigs(2,5))')
                idTvar(iPigs(2,5))=varid
              CASE ('idTvar(iPigs(3,5))')
                idTvar(iPigs(3,5))=varid
              CASE ('idTvar(iPigs(4,5))')
                idTvar(iPigs(4,5))=varid
              CASE ('idTvar(iPigs(4,7))')
                idTvar(iPigs(4,7))=varid
#endif
              CASE DEFAULT
                load=.false.
            END SELECT
!
!  Load variable data into information arrays.
!
            IF (load) THEN
              IF (varid.gt.MV) THEN
                WRITE (stdout,60) MV, varid
                STOP
              END IF
              DO i=1,5
                Vname(i,varid)=TRIM(ADJUSTL(Vinfo(i)))
              END DO            
              DO ng=1,Ngrids
                Iinfo(1,varid,ng)=gtype
                Fscale(varid,ng)=scale
              END DO
#ifdef SEDIMENT
!
!  Adjust information for all sediment tracers.
!

              SELECT CASE (Vinfo(1))
                CASE ('mud_')
                  IF (NCS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NCS
                      varid=varid+1
                      idTvar(idmud(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),90)                         &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),100)                        &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),110)                        &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),120)                        &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),110)                        &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('sand_')
                  IF (NNS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NNS
                      varid=varid+1
                      idTvar(isand(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),90)                         &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),100)                        &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),110)                        &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),120)                        &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),110)                        &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('mudfrac_')
                  IF (NCS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NCS
                      varid=varid+1
                      idfrac(i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),90)                         &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),100)                        &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),110)                        &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),120)                        &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),110)                        &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('sandfrac_')
                  IF (NNS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NNS
                      varid=varid+1
                      idfrac(NCS+i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),90)                         &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),100)                        &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),110)                        &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),120)                        &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),110)                        &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
              END SELECT
#endif
            ELSE
              varid=varid-1
!!            WRITE (stdout,70) TRIM(ADJUSTL(Vinfo(1))),                &
!!   &                          TRIM(ADJUSTL(Vinfo(6)))
            END IF
          END IF
        END DO
        GOTO 40
  30    WRITE (stdout,80) TRIM(ADJUSTL(Vinfo(1)))
        STOP
  40    CLOSE (inp)
!
  50    FORMAT (/,' MOD_NCPARAM - Unable to open variable information', &
     &          ' file: ',/,15x,a,/,15x,'Default file is located in',   &
     &          ' source directory.')
  60    FORMAT (/,' MOD_NCPARAM - too small dimension ',                &
     &          'parameter, MV = ',2i5,/,11x,                           &
     &          'change file  mod_ncparam.F  and recompile.')
  70    FORMAT (/,' MOD_NCPARM - Cannot load information for ',         &
     &          'variable: ',a,/,15x,'Need CASE construct for: ',a)
  80    FORMAT (/,' MOD_NCPARM - Error while reading information ',     &
     &          'for variable: ',a)
  90    FORMAT (a,i2.2)
 100    FORMAT (a,a,i2.2)
 110    FORMAT (a)
 120    FORMAT (a,a)

        RETURN
        END SUBROUTINE initialize_ncparam

      END MODULE mod_ncparam
