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
        integer, parameter :: MV = 350
        integer, parameter :: NV = 350
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

        integer  :: idAice         ! fraction of cell covered by ice
        integer  :: idCfra         ! cloud fraction
        integer  :: idDano         ! density anomaly
        integer  :: iddQdT         ! heat flux sensitivity to SST
        integer  :: idevap         ! evaporation rate
        integer  :: idFsur         ! free-surface
        integer  :: idHbbl         ! depth of bottom boundary layer
        integer  :: idHice         ! depth of ice cover
        integer  :: idHsbl         ! depth of surface boundary layer
        integer  :: idHsno         ! depth of snow cover
        integer  :: idLdwn         ! downwelling longwave radiation flux
        integer  :: idLhea         ! net latent heat flux
        integer  :: idLrad         ! net longwave radiation flux
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
        integer  :: idUbcs         ! bottom max U-momentum-wave stress
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
        integer  :: idVbcs         ! bottom max V-momentum-wave stress
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
        integer  :: idWorb         ! wind-induced wave orbital velocity
        integer  :: idWvel         ! true vertical velocity

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
!  Time-averaged quadratic terms IDs.
!
        integer  :: idUUav                    ! <uu> for 3D momentum 
        integer  :: idUVav                    ! <uv> for 3D momentum 
        integer  :: idVVav                    ! <vv> for 3D momentum 

        integer  :: idUTav(NAT)               ! <uT> for active tracers
        integer  :: idVTav(NAT)               ! <vT> for active tracers
!
!  Biological 2D diagnostic variable IDs.
!
        integer, allocatable :: iDbio2(:)     ! 2D biological terms

        integer  :: iCOfx = 1                 ! air-sea CO2 flux
        integer  :: iDNIT = 2                 ! denitrification flux
!
!  Biological 3D diagnostic variable IDs.
!
        integer, allocatable :: iDbio3(:)     ! 3D biological terms

        integer  :: iPPro = 1                 ! primary productivity
        integer  :: iNO3u = 2                 ! NO3 uptake
!
!  Tracer/Momentum Diagnostic variable IDs.
!
        integer, allocatable :: idDtrc(:,:)   ! tracers terms
        integer, allocatable :: idDu2d(:)     ! 2D u-momentum terms
        integer, allocatable :: idDv2d(:)     ! 2D v-momentum terms
        integer, allocatable :: idDu3d(:)     ! 3D u-momentum terms
        integer, allocatable :: idDv3d(:)     ! 3D v-momentum terms
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
        integer  :: ncdiaid(Ngrids)     ! output diagnostics
        integer  :: ncesid(Ngrids)      ! error subspace
        integer  :: ncfltid(Ngrids)     ! output floats
        integer  :: ncfrcid(NV,Ngrids)  ! input forcing
        integer  :: ncfwdid(Ngrids)     ! forward solution
        integer  :: ncgrdid(Ngrids)     ! input grid
        integer  :: nchisid(Ngrids)     ! output history
        integer  :: ncrstid(Ngrids)     ! input/output restart
        integer  :: ncstaid(Ngrids)     ! output stations

        integer  :: idefavg(Ngrids)     ! averages file creation flag
        integer  :: idefdia(Ngrids)     ! diagnostics file creation flag
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
        integer  :: diaVid(NV,Ngrids)    ! diagnostics variables IDs
        integer  :: fltVid(-6:NV,Ngrids) ! floats variables IDs
        integer  :: hisVid(NV,Ngrids)    ! history variables IDs
        integer  :: rstVid(NV,Ngrids)    ! restart variables IDs
        integer  :: staVid(NV,Ngrids)    ! stations variables IDs

        integer  :: tavgindx(Ngrids)     ! averages time record index
        integer  :: tdiaindx(Ngrids)     ! diagnostics time record index
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

        real(r8), allocatable :: Vtime(:,:,:)

        character (len=5  ) :: version = '2.1  '
        character (len=15 ) :: varnam(MV)
        character (len=44 ) :: date_str
        character (len=46 ) :: Tname(0:NV)
        character (len=80 ) :: Vname(5,0:NV)
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
#if defined SEDIMENT || defined BBL_MODEL
        USE mod_sediment
#endif
!
!  Local variable declarations.
!
        logical :: load

        integer, parameter :: inp = 10
        integer :: Lvar, gtype, i, j, ng, varid

        real(r8) :: scale

        character (len=80), dimension(7) :: Vinfo
!
!-----------------------------------------------------------------------
!  Allocate several indices variables.
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

        allocate ( iDbio2(NDbio2d) )
        allocate ( iDbio3(NDbio3d) )

        allocate ( idDtrc(MT,NDT) )
        allocate ( idDu2d(NDM2d) )
        allocate ( idDv2d(NDM2d) )
        allocate ( idDu3d(NDM3d) )
        allocate ( idDv3d(NDM3d) )

        allocate ( idSvar(MAXVAL(NSV)) )

#if defined SSH_TIDES || defined UV_TIDES
        allocate ( Vtime(MAX(2,MTC),NV,Ngrids) )
#else
        allocate ( Vtime(2,NV,Ngrids) )
#endif
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
          ncdiaid(ng)=-1
          ncesid(ng)=-1
          ncfltid(ng)=-1
          ncfwdid(ng)=-1
          ncgrdid(ng)=-1
          nchisid(ng)=-1
          ncrstid(ng)=-1
          ncstaid(ng)=-1

          tavgindx(ng)=0
          tdiaindx(ng)=0
          tfltindx(ng)=0
          thisindx(ng)=0
          trstindx(ng)=0
          tstaindx(ng)=0

          idefavg(ng)=-1
          idefdia(ng)=-1
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
              CASE ('idUbcs')
                idUbcs=varid
              CASE ('idVbws')
                idVbws=varid
              CASE ('idVbcs')
                idVbcs=varid
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
              CASE ('idLdwn')
                idLdwn=varid
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
              CASE ('idWorb')
                idWorb=varid
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
# ifdef AVERAGES_QUADRATIC
              CASE ('idUTav(itemp)')
                idUTav(itemp)=varid
              CASE ('idVTav(itemp)')
                idVTav(itemp)=varid
              CASE ('idUUav')
                idUUav=varid
              CASE ('idUVav')
                idUVav=varid
              CASE ('idVVav')
                idVVav=varid
# endif
# ifdef SALINITY
              CASE ('idSvar(isTvar(isalt))')
                idSvar(isTvar(isalt))=varid
#  ifdef AVERAGES_QUADRATIC
              CASE ('idUTav(isalt)')
                idUTav(isalt)=varid
              CASE ('idVTav(isalt)')
                idVTav(isalt)=varid
#  endif
# endif
#endif
#if defined SEDIMENT || defined BBL_MODEL
              CASE ('idTvar(idmud(i))')
                load=.true.
              CASE ('idTvar(isand(i))')
                load=.true.                
              CASE ('idfrac')
                load=.true.
              CASE ('idBmas')
                load=.true.
              CASE ('idSbed(ithck)')
                idSbed(ithck)=varid
              CASE ('idSbed(iaged)')
                idSbed(iaged)=varid
              CASE ('idSbed(iporo)')
                idSbed(iporo)=varid
              CASE ('idSbed(idiff)')
                idSbed(idiff)=varid
              CASE ('idBott(isd50)')
                idBott(isd50)=varid
              CASE ('idBott(idens)')
                idBott(idens)=varid
              CASE ('idBott(iwsed)')
                idBott(iwsed)=varid
              CASE ('idBott(itauc)')
                idBott(itauc)=varid
              CASE ('idBott(irlen)')
                idBott(irlen)=varid
              CASE ('idBott(irhgt)')
                idBott(irhgt)=varid
              CASE ('idBott(ibwav)')
                idBott(ibwav)=varid
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
              CASE ('idBott(izwbl)')
                idBott(izwbl)=varid
              CASE ('idBott(izdef)')
                idBott(izdef)=varid
              CASE ('idBott(iactv)')
                idBott(iactv)=varid
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
              CASE ('idTvar(iTIC_)')
                idTvar(iTIC_)=varid
              CASE ('idTvar(iTAlk)')
                idTvar(iTAlk)=varid
              CASE ('idTvar(iLDeC)')
                idTvar(iLDeC)=varid
              CASE ('idTvar(iSDeC)')
                idTvar(iSDeC)=varid
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
              CASE ('idTbry(iwest,iSDeC)')
                idTbry(iwest,iSDeC)=varid
              CASE ('idTbry(ieast,iSDeC)')
                idTbry(ieast,iSDeC)=varid
              CASE ('idTbry(isouth,iSDeC)')
                idTbry(isouth,iSDeC)=varid
              CASE ('idTbry(inorth,iSDeC)')
                idTbry(inorth,iSDeC)=varid
              CASE ('idTbry(iwest,iLDeC)')
                idTbry(iwest,iLDeC)=varid
              CASE ('idTbry(ieast,iLDeC)')
                idTbry(ieast,iLDeC)=varid
              CASE ('idTbry(isouth,iLDeC)')
                idTbry(isouth,iLDeC)=varid
              CASE ('idTbry(inorth,iLDeC)')
                idTbry(inorth,iLDeC)=varid
              CASE ('idTbry(iwest,iTIC_)')
                idTbry(iwest,iTIC_)=varid
              CASE ('idTbry(ieast,iTIC_)')
                idTbry(ieast,iTIC_)=varid
              CASE ('idTbry(isouth,iTIC_)')
                idTbry(isouth,iTIC_)=varid
              CASE ('idTbry(inorth,iTIC_)')
                idTbry(inorth,iTIC_)=varid
              CASE ('idTbry(iwest,iTAlk)')
                idTbry(iwest,iTAlk)=varid
              CASE ('idTbry(ieast,iTAlk)')
                idTbry(ieast,iTAlk)=varid
              CASE ('idTbry(isouth,iTAlk)')
                idTbry(isouth,iTAlk)=varid
              CASE ('idTbry(inorth,iTAlk)')
                idTbry(inorth,iTAlk)=varid
# ifdef DIAGNOSTICS_BIO
              CASE ('iDbio2(iCOfx)')
                iDbio2(iCOfx)=varid
              CASE ('iDbio2(iDNIT)')
                iDbio2(iDNIT)=varid
              CASE ('iDbio3(iPPro)')
                iDbio3(iPPro)=varid
              CASE ('iDbio3(iNO3u)')
                iDbio3(iNO3u)=varid
# endif
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
              CASE ('idTvar(iBacC)')
                load=.true.
              CASE ('idTvar(iBacF)')
                load=.true.
              CASE ('idTvar(iBacN)')
                load=.true.
              CASE ('idTvar(iBacP)')
                load=.true.
              CASE ('idTvar(iCDMC)')
                load=.true.
              CASE ('idTvar(iDOMC)')
                load=.true.
              CASE ('idTvar(iDOMN)')
                load=.true.
              CASE ('idTvar(iDOMP)')
                load=.true.
              CASE ('idTvar(iFecC)')
                load=.true.
              CASE ('idTvar(iFecF)')
                load=.true.
              CASE ('idTvar(iFecN)')
                load=.true.
              CASE ('idTvar(iFecP)')
                load=.true.
              CASE ('idTvar(iFecS)')
                load=.true.
              CASE ('idTvar(iPhyC)')
                load=.true.
              CASE ('idTvar(iPhyF)')
                load=.true.
              CASE ('idTvar(iPhyN)')
                load=.true.
              CASE ('idTvar(iPhyP)')
                load=.true.
              CASE ('idTvar(iPhyS)')
                load=.true.
              CASE ('idTvar(iPigs)')
#endif
#ifdef DIAGNOSTICS_UV
              CASE ('idDu2d(M2pgrd)')
                idDu2d(M2pgrd)=varid
              CASE ('idDv2d(M2pgrd)')
                idDv2d(M2pgrd)=varid
              CASE ('idDu2d(M2strs)')
                idDu2d(M2strs)=varid
              CASE ('idDv2d(M2strs)')
                idDv2d(M2strs)=varid
              CASE ('idDu2d(M2rate)')
                idDu2d(M2rate)=varid
              CASE ('idDv2d(M2rate)')
                idDv2d(M2rate)=varid
# ifdef UV_ADV
              CASE ('idDu2d(M2hadv)')
                idDu2d(M2hadv)=varid
              CASE ('idDv2d(M2hadv)')
                idDv2d(M2hadv)=varid
# endif
# ifdef UV_COR
              CASE ('idDu2d(M2fcor)')
                idDu2d(M2fcor)=varid
              CASE ('idDv2d(M2fcor)')
                idDv2d(M2fcor)=varid
# endif
# if defined UV_VIS2 || defined UV_VIS4
              CASE ('idDu2d(M2hvis)')
                idDu2d(M2hvis)=varid
              CASE ('idDv2d(M2hvis)')
                idDv2d(M2hvis)=varid
# endif
# ifdef SOLVE3D
              CASE ('idDu3d(M3pgrd)')
                idDu3d(M3pgrd)=varid
              CASE ('idDv3d(M3pgrd)')
                idDv3d(M3pgrd)=varid
              CASE ('idDu3d(M3vvis)')
                idDu3d(M3vvis)=varid
              CASE ('idDv3d(M3vvis)')
                idDv3d(M3vvis)=varid
              CASE ('idDu3d(M3rate)')
                idDu3d(M3rate)=varid
              CASE ('idDv3d(M3rate)')
                idDv3d(M3rate)=varid
#  ifdef UV_ADV
              CASE ('idDu3d(M3hadv)')
                idDu3d(M3hadv)=varid
              CASE ('idDv3d(M3hadv)')
                idDv3d(M3hadv)=varid
              CASE ('idDu3d(M3vadv)')
                idDu3d(M3vadv)=varid
              CASE ('idDv3d(M3vadv)')
                idDv3d(M3vadv)=varid
#  endif
#  ifdef UV_COR
              CASE ('idDu3d(M3fcor)')
                idDu3d(M3fcor)=varid
              CASE ('idDv3d(M3fcor)')
                idDv3d(M3fcor)=varid
#  endif
#  if defined UV_VIS2 || defined UV_VIS4
              CASE ('idDu3d(M3hvis)')
                idDu3d(M3hvis)=varid
              CASE ('idDv3d(M3hvis)')
                idDv3d(M3hvis)=varid
#  endif
# endif
#endif
#ifdef DIAGNOSTICS_TS
              CASE ('idDtrc(iTrate)')
                load=.true.
              CASE ('idDtrc(iThadv)')
                load=.true.
              CASE ('idDtrc(iTvadv)')
                load=.true.
# if defined TS_DIF2 || defined TS_DIF4
              CASE ('idDtrc(iThdif)')
                load=.true.
# endif
              CASE ('idDtrc(iTvdif)')
                load=.true.
#endif
              CASE DEFAULT
                load=.false.
            END SELECT
!
!  Load variable data into information arrays.
!
            IF (load) THEN
              load=.false.
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
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
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
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
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
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
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
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('mudmass_')
                  IF (NCS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NCS
                      varid=varid+1
                      idBmas(i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
                CASE ('sandmass_')
                  IF (NNS.gt.0) THEN
                    varid=varid-1
                    DO i=1,NNS
                      varid=varid+1
                      idBmas(i)=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i2.2)')                 &
     &                      TRIM(ADJUSTL(Vinfo(1))), i
                      WRITE (Vname(2,varid),'(a,a,i2.2)')               &
     &                      TRIM(ADJUSTL(Vinfo(2))), ', size class ', i
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END DO
                  END IF
              END SELECT
#endif
#ifdef ECOSIM
!
!  Adjust information for all EcoSim tracers.
!
              SELECT CASE (Vinfo(1))
                CASE ('Bac_C')
                  varid=varid-1
                  DO i=1,Nbac
                    varid=varid+1
                    idTvar(iBacC(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Bac_C', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'bacteria, carbon group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Bac_F')
                  varid=varid-1
                  DO i=1,Nbac
                    varid=varid+1
                    idTvar(iBacF(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Bac_F', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'bacteria, iron group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Bac_N')
                  varid=varid-1
                  DO i=1,Nbac
                    varid=varid+1
                    idTvar(iBacN(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Bac_N', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'bacteria, nitrogen group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Bac_P')
                  varid=varid-1
                  DO i=1,Nbac
                    varid=varid+1
                    idTvar(iBacP(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Bac_P', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'bacteria, phosphorus group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('CDM_C')
                  varid=varid-1
                  DO i=1,Ndom
                    varid=varid+1
                    idTvar(iCDMC(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'CDM_C', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'color degradational matter, carbon group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('DOM_C')
                  varid=varid-1
                  DO i=1,Ndom
                    varid=varid+1
                    idTvar(iDOMC(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'DOM_C', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'dissolved organic matter, carbon group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('DOM_N')
                  varid=varid-1
                  DO i=1,Ndom
                    varid=varid+1
                    idTvar(iDOMN(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'DOM_N', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'dissolved organic matter, nitrogen group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('DOM_P')
                  varid=varid-1
                  DO i=1,Ndom
                    varid=varid+1
                    idTvar(iDOMP(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'DOM_P', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                  'dissolved organic matter, phosphorus group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Fec_C')
                  varid=varid-1
                  DO i=1,Nfec
                    varid=varid+1
                    idTvar(iFecC(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Fec_C', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'fecal matter, carbon group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Fec_F')
                  varid=varid-1
                  DO i=1,Nfec
                    varid=varid+1
                    idTvar(iFecF(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Fec_F', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'fecal matter, iron group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Fec_N')
                  varid=varid-1
                  DO i=1,Nfec
                    varid=varid+1
                    idTvar(iFecN(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Fec_N', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'fecal matter, nitrogen group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Fec_P')
                  varid=varid-1
                  DO i=1,Nfec
                    varid=varid+1
                    idTvar(iFecP(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Fec_P', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'fecal matter, phosphorus group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Fec_S')
                  varid=varid-1
                  DO i=1,Nfec
                    varid=varid+1
                    idTvar(iFecS(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Fec_S', i
                    WRITE (Vname(2,varid),'(a,i1)')                     &
     &                    'fecal matter, silica group ', i
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Phy_C')
                  varid=varid-1
                  DO i=1,Nphy
                    varid=varid+1
                    idTvar(iPhyC(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Phy_C', i
                    WRITE (Vname(2,varid),'(a,a)')                      &
     &                    TRIM(PhyName(i)), ', carbon group'
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Phy_F')
                  varid=varid-1
                  DO i=1,Nphy
                    varid=varid+1
                    idTvar(iPhyF(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Phy_F', i
                    WRITE (Vname(2,varid),'(a,a)')                      &
     &                    TRIM(PhyName(i)), ', iron group'
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Phy_N')
                  varid=varid-1
                  DO i=1,Nphy
                    varid=varid+1
                    idTvar(iPhyN(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Phy_N', i
                    WRITE (Vname(2,varid),'(a,a)')                      &
     &                    TRIM(PhyName(i)), ', nitrogen group'
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Phy_P')
                  varid=varid-1
                  DO i=1,Nphy
                    varid=varid+1
                    idTvar(iPhyP(i))=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,i1)') 'Phy_P', i
                    WRITE (Vname(2,varid),'(a,a)')                      &
     &                    TRIM(PhyName(i)), ', phosphorus group'
                    WRITE (Vname(3,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('Phy_S')
                  j=0
                  DO i=1,Nphy
                    IF (PHY(i).le.2) j=j+1
                  END DO                
                  IF (j.gt.0) varid=varid-1
                  DO i=1,Nphy
                    IF (PHY(i).le.2) THEN
                      varid=varid+1
                      idTvar(iPhyS(i))=varid
                      DO ng=1,Ngrids
                        Fscale(varid,ng)=scale
                        Iinfo(1,varid,ng)=gtype
                      END DO
                      WRITE (Vname(1,varid),'(a,i1)') 'Phy_S', i
                      WRITE (Vname(2,varid),'(a,a)')                    &
     &                      TRIM(PhyName(i)), ', silica group'
                      WRITE (Vname(3,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(3)))
                      WRITE (Vname(4,varid),'(a,a)')                    &
     &                      TRIM(Vname(1,varid)), ', scalar, series'
                      WRITE (Vname(5,varid),'(a)')                      &
     &                      TRIM(ADJUSTL(Vinfo(5)))
                    END IF
                  END DO
                CASE ('Pigs_XY')
                  varid=varid-1
                  DO j=1,Npig 
                    DO i=1,Nphy
                      IF (PIG(PHY(i),j).gt.0) THEN
                        varid=varid+1
                        idTvar(iPigs(i,j))=varid
                        DO ng=1,Ngrids
                          Fscale(varid,ng)=scale
                          Iinfo(1,varid,ng)=gtype
                        END DO
                        WRITE (Vname(1,varid),'(a,i1,i1)') 'Pigs_', i, j
                        WRITE (Vname(2,varid),'(a,a,a)')                &
     &                        TRIM(PhyName(i)), ', ', TRIM(PigName(j))
                        WRITE (Vname(3,varid),'(a)')                    &
     &                        TRIM(ADJUSTL(Vinfo(3)))
                        WRITE (Vname(4,varid),'(a,a)')                  &
     &                        TRIM(Vname(1,varid)), ', scalar, series'
                        WRITE (Vname(5,varid),'(a)')                    &
     &                        TRIM(ADJUSTL(Vinfo(5)))
                      END IF
                    END DO
                  END DO
              END SELECT
#endif
#ifdef DIAGNOSTICS_TS
!
!  Adjust information for tracer diagnostic variables.  This needs to be 
!  done last because it needs all the tracers variable names.
!

              SELECT CASE (Vinfo(1))
                CASE ('_rate')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTrate)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('_hadv')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iThadv)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
                CASE ('_vadv')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTvadv)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
# if defined TS_DIF2 || defined TS_DIF4
                CASE ('_hdiff')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iThdif)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                      &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
# endif
                CASE ('_vdiff')
                  varid=varid-1
                  DO i=1,MT
                    varid=varid+1
                    idDtrc(i,iTvdif)=varid
                    DO ng=1,Ngrids
                      Fscale(varid,ng)=scale
                      Iinfo(1,varid,ng)=gtype
                    END DO
                    WRITE (Vname(1,varid),'(a,a)')                      &
     &                    TRIM(ADJUSTL(Vname(1,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(1)))
                    WRITE (Vname(2,varid),'(a,", ",a)')                 &
     &                    TRIM(ADJUSTL(Vname(2,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(2)))
                    WRITE (Vname(3,varid),'(a,1x,a)')                   &
     &                    TRIM(ADJUSTL(Vname(3,idTvar(i)))),            &
     &                    TRIM(ADJUSTL(Vinfo(3)))
                    WRITE (Vname(4,varid),'(a,a)')                      &
     &                    TRIM(Vname(1,varid)), ', scalar, series'
                    WRITE (Vname(5,varid),'(a)')                        &
     &                    TRIM(ADJUSTL(Vinfo(5)))
                  END DO
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
     &          'parameter, MV = ',2i5,/,15x,                           &
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
