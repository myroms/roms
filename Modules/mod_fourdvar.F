#include "cppdefs.h"
      MODULE mod_fourdvar

#if defined FOUR_DVAR || defined GRADIENT_CHECK || defined TLM_CHECK
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  Variational data assimilation variables:                            !
!                                                                      !
!  ADmodVal     Adjoint model values at observation locations.         !
!  BackFun      Cost function misfit (mean squared difference) between !
!                 model and background state.                          !
!  DTdif        Time-step (s) for diffusion equation space convolution.!
!  CostFun      Cost function misfit (mean squared difference) between !
!                 model and observations.                              !
!  CostFunOld   Cost function misfit (mean squared difference) between !
!                 model and observations from previous iteration.      !
!  CovNorm      Isotropic covariance normalization factor used in the  !
!                 generalized diffusion equation convolution.          !
!  Hmat         Interpolation weights matrix, Hmat(8,:),               !
!                                                                      !
!                              8____________7                          !
!                              /.          /| (i2,j2,k2)               !
!                             / .         / |                          !
!                           5/___________/6 |                          !
!                            |  .        |  |                          !
!                            |  .        |  |         Grid Cell        !
!                            | 4.........|..|3                         !
!                            | .         |  /                          !
!                            |.          | /                           !
!                 (i1,j1,k1) |___________|/                            !
!                            1           2                             !
!                                                                      !
!  Ldecay       Isotropic covariance horizontal decorrelation scale    !
!                 (m).                                                 !
!  NLmodVal     Nonlinear model values at observation locations.       !
!  ObsCount     Current observation counts per state variable.         !
!  ObsErr       Observation error.                                     !
!  ObsNorm      Observation norm used to scale cost function and       !
!                 adjoint misfit forcing terms.                        !
!  ObsType      Observation type identifier.                           !
!  ObsVal       Observation values.                                    !
!  ObsVar       Global observation variance for each state variable.   !
!  TLmodVal     Tangent linear model values at observation locations.  !
!  Xobs         Observations X-locations (grid coordinates).           !
!  Yobs         Observations Y-locations (grid coordinates).           !
!  Zobs         Observations Z-locations (grid coordinates).           !
!                                                                      !
!=======================================================================
!
        USE mod_param
!
        implicit none

        integer,  allocatable :: ObsType(:)

        real(r8), allocatable :: ObsErr(:)
        real(r8), allocatable :: ObsNorm(:)
        real(r8), allocatable :: ObsVal(:)
        real(r8), allocatable :: Xobs(:)
        real(r8), allocatable :: Yobs(:)
# ifdef SOLVE3D
        real(r8), allocatable :: Zobs(:)
# endif          
        real(r8), allocatable :: Hmat(:,:)
        real(r8), allocatable :: ADmodVal(:)
        real(r8), allocatable :: NLmodVal(:,:)
# ifdef TLM_OBS
        real(r8), allocatable :: TLmodVal(:,:)
# endif
!
        TYPE T_FOURDVAR

          integer , pointer :: NobsSurvey(:)
          integer , pointer :: ObsCount(:)

          real(r8), pointer :: BackFun(:)
          real(r8), pointer :: CostFun(:)
          real(r8), pointer :: CostFunOld(:)
          real(r8), pointer :: ObsVar(:)
          real(r8), pointer :: SurveyTime(:)
  
        END TYPE T_FOURDVAR

        TYPE (T_FOURDVAR), allocatable :: FOURDVAR(:)
!
!-----------------------------------------------------------------------
!  Observations parameters.
!-----------------------------------------------------------------------
!
!  Maximum number of observations to process.
!
        integer :: Mobs
!
!  Number of model state variables to process.
!
        integer, dimension(Ngrids) :: NstateVar
!
!  Number of interpolation weights and (I,J,K) indices offsets.
!
# ifdef SOLVE3D
        integer, parameter :: Nweights = 8

        integer, parameter, dimension(Nweights) :: Ioffset =            &
     &                      (/ 0, 1, 0, 1, 0, 1, 0, 1 /)
        integer, parameter, dimension(Nweights) :: Joffset =            &
     &                      (/ 0, 0, 1, 1, 0, 0, 1, 1 /)
        integer, parameter, dimension(Nweights) :: Koffset =            &
     &                      (/ 0, 0, 0, 0, 1, 1, 1, 1 /)
# else
        integer, parameter :: Nweights = 4

        integer, parameter, dimension(Nweights) :: Ioffset =            &
     &                      (/ 0, 1, 0, 1 /)
        integer, parameter, dimension(Nweights) :: Joffset =            &
     &                      (/ 0, 0, 1, 1 /)
# endif
!        
!  (I,J) lower and upper bound ranges for processing of RHO-type
!  observations.
!
        integer, dimension(Ngrids) :: rILB
        integer, dimension(Ngrids) :: rIUB
        integer, dimension(Ngrids) :: rJLB
        integer, dimension(Ngrids) :: rJUB
!
!  (I,J) lower and upper bound ranges for processing of U-type
!  observations.
!
        integer, dimension(Ngrids) :: uILB
        integer, dimension(Ngrids) :: uIUB
        integer, dimension(Ngrids) :: uJLB
        integer, dimension(Ngrids) :: uJUB
!
!  (I,J) lower and upper bound ranges for processing of U-type
!  observations.
!
        integer, dimension(Ngrids) :: vILB
        integer, dimension(Ngrids) :: vIUB
        integer, dimension(Ngrids) :: vJLB
        integer, dimension(Ngrids) :: vJUB
!
!  Size of observation NetCDF file unlimited dimension.
!
        integer, dimension(Ngrids) :: Ndatum         
!
!  Number of observations surveys available.
!
        integer, dimension(Ngrids) :: Nsurvey
!
!  Observation surveys counter.
!
        integer, dimension(Ngrids) :: ObsSurvey
!
!  Current number of observations processed.
!
        integer, dimension(Ngrids) :: Nobs
!
!  Current starting and ending observation file index processed.
!
        integer, dimension(Ngrids) :: NstrObs
        integer, dimension(Ngrids) :: NendObs
!
!  Switch to activate the processing of model state at the observation
!  locations.
!
        logical, dimension(Ngrids) :: ProcessObs
!
!   Switch to activate writting of nonlinear model values at
!   observations locations into observations NetCDF file.
!
        logical, dimension(Ngrids) :: wrtNLmod
!
!   Switch to activate writting of tangent linear model values at
!   observations locations into observations NetCDF file.
!
        logical, dimension(Ngrids) :: wrtTLmod
!
!-----------------------------------------------------------------------
!  Descent algorithm parameters.
!-----------------------------------------------------------------------
!
!  Number of iteration between each steepest descent step.
!
        integer :: NiterSD
!
!  Iteration of the last steepest descent.
!
        integer :: IterSD
!
!  Conjugate gradient algorithm: [0] Fletcher-Reeves [1] Polak-Riviere
!
        integer :: ICG
!
!  Pass number into the descent algorithm within the same iteration.
!
        integer :: Ipass
!
!  Conjugate gradient alogorithm step size.
!
        real(r8) :: CGstepI            ! initial step size
        real(r8) :: CGstepF            ! first guess step size
        real(r8) :: CGstepR            ! refined step size
!
!  Step size factors.
!
        real(r8) :: StepSizeB          ! denominator
        real(r8) :: StepSizeT          ! numerator
!
!  Conjugate gradient orthogonality tolerance: dot product orthogonality
!  between current and previous search.
!
        real(r8) :: CGtol
!
!  Cost function percentage change between succesive iterations for
!  convergence.
!
        real(r8) :: CostFunFac
!
!  Brunt-Vaisala expansion polynomial used for cost function
!  normalization.
!
        integer :: npN2
        real(r8), allocatable :: pcoef_N2(:)
!
!------------------------------------------------------------------------
!  Covariance parameters
!------------------------------------------------------------------------
!
!  Isotropic covariance normalization factor used in the generalized
!  diffusion equation space convolution.
!
        real(r8), dimension(Ngrids) :: CovNorm
!
!  Time-step (seconds) for diffusion equation space convolution.
!
        real(r8), dimension(Ngrids) :: DTdif
!
!  Isotropic covariance horizontal decorrelation scale (m).
!
        real(r8), dimension(Ngrids) :: Ldecay
!
!------------------------------------------------------------------------
!  Dot product parameters.
!------------------------------------------------------------------------
!
!  Dot product between tangent linear and adjoint vectors.
!
        real(r8) :: DotProduct
        real(r8) :: adDotProduct
!
!  Tangent linear model linearization check dot products.
!
        integer :: ig1count            ! counter for g1 dot product
        integer :: ig2count            ! counter for g2 dot product

        real(r8), dimension(1000) :: g1
        real(r8), dimension(1000) :: g2

      CONTAINS

      SUBROUTINE initialize_fourdvar
!
!=======================================================================
!  Copyright (c) 2004 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine initializes several variables in module "mod_fourdvar" !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
!
!  Local variable declarations.
!
      integer :: i, ndims, ng, ngatts, recdim, status, varid
      integer :: dimid, dimsiz

      real(r8), parameter :: IniVal = 0.0_r8

      character (len=20) :: dimnam
      character (len=80) :: fname
!
!-----------------------------------------------------------------------
!  Inquire observations NetCDF and determine the maxumum dimension of
!  several observations arrays.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
!
!  Open observation NetCDF.
!
        status=nf_open (TRIM(OBSname(ng)),nf_nowrite,ncOBSid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(OBSname(ng))
          exit_flag=4
          ioerror=status
          RETURN
        END IF
!
!  Inquire about the size of the "datum" unlimitted dimension and
!  "survey" dimension.
!
        Ndatum(ng)=0
        Nsurvey(ng)=0
        status=nf_inq(ncOBSid(ng),ndims,nvars,ngatts,recdim)
        IF (status.eq.nf_noerr) THEN
          DO i=1,ndims
            dimid=i
            status=nf_inq_dim(ncOBSid(ng),dimid,dimnam,dimsiz)
            IF (TRIM(dimnam).eq.'datum') then
              Ndatum(ng)=dimsiz
            ELSE IF (TRIM(dimnam).eq.'survey') THEN
              Nsurvey(ng)=dimsiz
            END IF
          END DO
          IF (Ndatum(ng).eq.0) THEN
            WRITE (stdout,20) 'datum', TRIM(OBSname(ng))
            exit_flag=4
          END IF
          IF (Nsurvey(ng).eq.0) THEN
            WRITE (stdout,20) 'survey', TRIM(OBSname(ng))
            exit_flag=4
          END IF
        ELSE
          WRITE(stdout,30) OBSname(ng)
        END IF
      END DO
!
!  Allocate module structure.
!
      allocate ( FOURDVAR(Ngrids) )
!
!  Allocate vectors to store the number of observations per survey
!  and their times.
!
      DO ng=1,Ngrids

        allocate ( FOURDVAR(ng) % NobsSurvey(Nsurvey(ng)) )
        FOURDVAR(ng) % NobsSurvey = 0

        allocate ( FOURDVAR(ng) % SurveyTime(Nsurvey(ng)) )
        FOURDVAR(ng) % NobsSurvey = IniVal
!
!  Read in number of observations available per survey.
!
        status=nf_inq_varid(ncOBSid(ng),TRIM(Vname(1,idNobs)),varid)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,40) TRIM(Vname(1,idNobs)), TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
        status=nf_get_vara_int(ncOBSid(ng), varid, 1, Nsurvey(ng),      &
     &                         FOURDVAR(ng)%NobsSurvey)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,50) TRIM(Vname(1,idNobs)), TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
!
!  Determine maximum size of observation arrays.
!
        Mobs=0
        DO i=1,Nsurvey(ng)
          Mobs=MAX(Mobs, FOURDVAR(ng)%NobsSurvey(i))
        END DO
!
!  Close observation NetCDF file.
!
        status=nf_close (ncOBSid(ng))
        ncOBSid(ng)=-1
      END DO
!
!  Allocate and initialize model and observation variables. 
!
      allocate ( ObsType(Mobs) )
      ObsType = 0

      allocate ( ObsErr(Mobs) )
      ObsErr = IniVal

      allocate ( ObsNorm(Mobs) )
      ObsNorm = IniVal

      allocate ( ObsVal(Mobs) )
      ObsVal = IniVal

      allocate ( Xobs(Mobs) )
      Xobs = IniVal
 
      allocate ( Yobs(Mobs) )
      Yobs = IniVal

# ifdef SOLVE3D
      allocate ( Zobs(Mobs) )
      Zobs = IniVal
# endif
 
      allocate ( Hmat(Nweights,Mobs) )
      Hmat = IniVal

      allocate ( ADmodVal(Mobs) )
      ADmodVal = IniVal

      allocate ( NLmodVal(Mobs,2) )
      ADmodVal = IniVal

# ifdef TLM_OBS
      allocate ( TLmodVal(Mobs,2) )
      ADmodVal = IniVal
# endif
!
!-----------------------------------------------------------------------
!  Allocate remaining structure variables.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids

# ifdef SOLVE3D
        NstateVar(ng)=5+NT(ng)
# else
        NstateVar(ng)=3
# endif

        allocate ( FOURDVAR(ng) % BackFun(0:NstateVar(ng)) )
        FOURDVAR(ng) % BackFun = IniVal

        allocate ( FOURDVAR(ng) % CostFun(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostFun = IniVal

        allocate ( FOURDVAR(ng) % CostFunOld(0:NstateVar(ng)) )
        FOURDVAR(ng) % CostFunOld = IniVal

        allocate ( FOURDVAR(ng) % ObsVar(NstateVar(ng)) )
        FOURDVAR(ng) % ObsVar = IniVal

        allocate ( FOURDVAR(ng) % ObsCount(NstateVar(ng)) )
        FOURDVAR(ng) % ObsCount = 0

      END DO
!
!-----------------------------------------------------------------------
!  Set lower and upper observation processing (I,J) ranges: Either the
!  full grid or only interior points will be considered.
!-----------------------------------------------------------------------
!   
      DO ng=1,Ngrids
# ifdef FULL_GRID
        rILB(ng)=0
        rIUB(ng)=Lm(ng)+1
        rJLB(ng)=0
        rJUB(ng)=Mm(ng)+1

        uILB(ng)=1
        uIUB(ng)=Lm(ng)+1
        uJLB(ng)=0
        uJUB(ng)=Mm(ng)+1

        vILB(ng)=0
        vIUB(ng)=Lm(ng)+1
        vJLB(ng)=1
        vJUB(ng)=Mm(ng)+1
# else
        rILB(ng)=1
        rIUB(ng)=Lm(ng)
        rJLB(ng)=1
        rJUB(ng)=Mm(ng)

        uILB(ng)=2
        uIUB(ng)=Lm(ng)
        uJLB(ng)=1
        uJUB(ng)=Mm(ng)

        vILB(ng)=1
        vIUB(ng)=Lm(ng)
        vJLB(ng)=2
        vJUB(ng)=Mm(ng)
# endif
      END DO
!
!-----------------------------------------------------------------------
!  Initialize various variables.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        ProcessObs(ng)=.FALSE.
        wrtNLmod(ng)=.FALSE.
        wrtTLmod(ng)=.FALSE.
      END DO
!
 10   FORMAT (/,' MOD_FOURDVAR - unable to open input NetCDF file: ',a)
 20   FORMAT (/,' MOD_FOURDVAR - error inquiring dimension: ',a,2x,      &
     &          ' in input NetCDF file: ',a)
 30   FORMAT (/,' MOD_FOURDVAR - unable to inquire size of dimension',   &
     &          ' in variable: ',a)
 40   FORMAT (/,' MOD_FOURDVAR - error inquiring variable ID for: ',a,   &
     &        2x,' in input NetCDF file: ',a)
 50   FORMAT (/,' MOD_FOURDVAR - error while reading variable: ',a,2x,   &
     &          ' in input NetCDF file: ',a)

      RETURN
      END SUBROUTINE initialize_fourdvar
#endif
      END MODULE mod_fourdvar
