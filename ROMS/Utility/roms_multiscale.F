#include "cppdefs.h"
      MODULE roms_multiscale_mod
#ifdef MULTI_SCALE_B
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group            Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
! The module implements routines and objects for modeling multiscale   !
! background-error covariance (B) through implicit diffusion-based     !
! operators, as outlined by Weaver et al. (2016, 2018) and Mirouze     !
! et al. (2016).                                                       !
!                                                                      !
! The implicit two-dimensional diffusion solver, applied level by      !
! level, utilizes Chebyshev iterations (CI), which require estimates   !
! of the maximum and minimum eigenvalues of the K-Laplacian operator.  !
! These eigenvalues are precomputed for a fixed grid and specified     !
! k-values (diffusion coefficient with units of length squared) using  !
! the Lanczos variant of the conjugate gradients (CG) iterative        !
! algorithm, initialized with random vectors, to obtain the necessary  !
! eigenvalue estimates.                                                !
!                                                                      !
! This combined CG and CI methodology is employed for modeling the     !
! implicit multiscale backgrond-error covariance.                      !
!                                                                      !
! Horizontal spatially varying isotropic or anisotropic correlation    !
! scales (HscaleX, HscaleY) can be specified for each variable in the  !
! control vector (zeta, u, v, T, S). When provided, these correlation  !
! scales are read from an input NetCDF file. Error correlations are    !
! assumed to be separable in the horizontal and vertical directions.   !
!                                                                      !
! References:                                                          !
!                                                                      !
! Mirouze, I., E. Blockley, D.J. Lea, M.J. Martin, and M.J. Bell,      !
!   2016: A multiple length scale correlation operator for ocean data  !
!   assimilation, Tellus A, 68, 29744, doi:10.3402/tellusa.v68.29744.  !
!                                                                      !
! Weaver, A.T., J. Tshimanga, and A. Piacentini, 2016: Correlation     !
!   operators based on an implicitly formulated diffusion equation     !
!   solved with the Chebyshev iteration, Q. J. R. Meteorol. Soc., 142, !
!   455-471, doi:10.1002/qj.2664.                                      !
!                                                                      !
! Weaver, A.T., S. Gurol, J. Tshimanga, M. Chrust, and A. Piacentini,  !
!   2018: "Time"-parallel diffusion-based correlation operators. Q. J. !
!   R. Meteorol. Soc., 144, 2067-2088, doi:10.1002/qj.3302.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
# ifdef DISTRIBUTE
      USE distribute_mod,       ONLY : mp_reduce
# endif
      USE strings_mod,          ONLY : FoundError
!
      implicit none
!
!-----------------------------------------------------------------------
!  ROMS multiscale CLASS object.
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: multiscale
!
!  Kernel variables lower and upper bound array dimensions, and parallel
!  tile partition indices.
!
        integer :: tile                              ! domain partition
        integer :: LBi,  UBi,  LBj,  UBj             ! array bounds
# ifdef ADJUST_BOUNDARY
        integer :: LBij,  UBij                       ! boundary bounds
# endif
!
!  Multiscale background covariance matrix parameters:
!
!    Nscale   number of correlation scales to be combined
!    Mlap     nunber of k-Laplacian smoother applications
!    NiterCG  number of CG iterations to obtain each scale eigenvalue
!    NiterCI  number of implicit Chebyshev solver iterations
!    Bwgt     weights assigned to each correlation function
!
        integer :: Nscale
        integer, allocatable :: Mlap(:,:)         ! [MstateVar, Nscale]
        integer, allocatable :: NiterCG(:)        ! [Nscale]
        integer, allocatable :: NiterCI(:)        ! [Nscale]
!
        real(r8), allocatable :: Bwgt(:,:)        ! [MstateVar, Nscale]
!
!  The maximum and minimum eigenvalues of the K-Laplacian for each 2D
!  variable in the control vector as functions of both the number of
!  iterations and the specified correlation functions, [Mlap,Nscales].
!
        real(r8), pointer :: eigMaxR2D(:,:)       ! free surface
        real(r8), pointer :: eigMinR2D(:,:)
!
        real(r8), pointer :: eigMaxU2D(:,:)       ! 2D U-momentum
        real(r8), pointer :: eigMinU2D(:,:)
!
        real(r8), pointer :: eigMaxV2D(:,:)       ! 2D V-momentum
        real(r8), pointer :: eigMinV2D(:,:)

# ifdef ADJUST_WSTRESS
!
        real(r8), pointer :: eigMaxSUS(:,:)       ! surface U-stress
        real(r8), pointer :: eigMinSUS(:,:)
!
        real(r8), pointer :: eigMaxSVS(:,:)       ! surface V-stress
        real(r8), pointer :: eigMinSVS(:,:)
# endif

# ifdef SOLVE3D
!
!  The maximum and minimum eigenvalues of the K-Laplacian for each 3D
!  variable in the control vector as functions of both the number of
!  iterations, vertical levels, and the specified correlation functions,
!  [Mlap,N,Nscales] and [Mlap,N,NT,Nscale].
!
        real(r8), pointer :: eigMaxR3D(:,:,:,:)   ! tracer variables
        real(r8), pointer :: eigMinR3D(:,:,:,:)
!
        real(r8), pointer :: eigMaxU3D(:,:,:)     ! 3D U-momentum
        real(r8), pointer :: eigMinU3D(:,:,:)
!
        real(r8), pointer :: eigMaxV3D(:,:,:)     ! 3D V-momentum
        real(r8), pointer :: eigMinV3D(:,:,:)

#  ifdef ADJUST_STFLUX
!
        real(r8), pointer :: eigMaxSTF(:,:,:)     ! surf. Tracer flux
        real(r8), pointer :: eigMinSTF(:,:,:)
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
!  The maximum and minimum eigenvalues of the K-Laplacian for each 2D
!  variable lateral boundary adjustment in the control vector as
!  functions of both the number of iterations and the specified
!  correlation functions, [Mlap,4,Nscales].
!
        real(r8), pointer :: eigMaxR1D(:,:,:)     ! free surface
        real(r8), pointer :: eigMinR1D(:,:,:)
!
        real(r8), pointer :: eigMaxU1D(:,:,:)     ! 2D U-momentum
        real(r8), pointer :: eigMinU1D(:,:,:)
!
        real(r8), pointer :: eigMaxV1D(:,:,:)     ! 3D V-momentum
        real(r8), pointer :: eigMinV1D(:,:,:)

#  ifdef SOLVE3D
!
!  The maximum and minimum eigenvalues of the K-Laplacian for each 3D
!  variable lateral boundary adjustment in the control vector as
!  functions of both the number of iterations and the specified
!  correlation functions, [Mlap,N,4,Nscales] and [Mlap,N,NT,4,Nscales].
!
        real(r8), pointer :: eigMaxR1DZ(:,:,:,:,:)   ! tracers
        real(r8), pointer :: eigMinR1DZ(:,:,:,:,:)
!
        real(r8), pointer :: eigMaxU1DZ(:,:,:,:)     ! 3D U-momentum
        real(r8), pointer :: eigMinU1DZ(:,:,:,:)
!
        real(r8), pointer :: eigMaxV1DZ(:,:,:,:)     ! 3D V-momentum
        real(r8), pointer :: eigMinV1DZ(:,:,:,:)
#  endif
# endif
!
!  K-Laplacian Conjugate Gradient (CG) iterative solver horizontal
!  vectors, 2D variables in control vector [LBi:UBi,LBj:UBj].
!
        real(r8), allocatable :: cg2d_p(:,:)         ! direction, p(k)
        real(r8), allocatable :: cg2d_q(:,:)         ! auxiliary, A*p(k)
        real(r8), allocatable :: cg2d_r(:,:)         ! residual,  r(k)
        real(r8), allocatable :: cg2d_x(:,:)         ! estimate,  x(k)

# ifdef ADJUST_BOUNDARY
!
!  K-Laplacian Conjugate Gradient (CG) iterative solver horizontal
!  vectors, 1D lateral boundary adjustment [LBij:UBij].
!
        real(r8), allocatable :: cgB1d_p(:)          ! direction, p(k)
        real(r8), allocatable :: cgB1d_q(:)          ! auxiliary, A*p(k)
        real(r8), allocatable :: cgB1d_r(:)          ! residual,  r(k)
        real(r8), allocatable :: cgB1d_x(:)          ! estimate,  x(k)

#  ifdef SOLVE3D
!
!  K-Laplacian Conjugate Gradient (CG) iterative solver horizontal
!  vectors, 2D lateral boundary adjustment [LBij:UBij,N(ng)].
!
        real(r8), allocatable :: cgB2d_p(:,:)        ! direction, p(k)
        real(r8), allocatable :: cgB2d_q(:,:)        ! auxiliary, A*p(k)
        real(r8), allocatable :: cgB2d_r(:,:)        ! residual,  r(k)
        real(r8), allocatable :: cgB2d_x(:,:)        ! estimate,  x(k)
#  endif
# endif

# ifdef SOLVE3D
!
!  K-Laplacian Conjugate Gradient (CG) iterative solver horizontal
!  vectors, 3D variables in control vector [LBi:UBi,LBj:UBj,N(ng)].
!
        real(r8), allocatable :: cg3d_p(:,:,:)       ! direction, p(k)
        real(r8), allocatable :: cg3d_q(:,:,:)       ! auxiliary, A*p(k)
        real(r8), allocatable :: cg3d_r(:,:,:)       ! residual,  r(k)
        real(r8), allocatable :: cg3d_x(:,:,:)       ! estimate,  x(k)
# endif
!
!  Chebyshev iterations (CI) polynomial-based solver horizontal vectors,
!  2D variables in control vector, [LBi:UBi,LBj:UBj].
!
        real(r8), allocatable :: ci2d_p(:,:)         ! direction, p(k)
        real(r8), allocatable :: ci2d_q(:,:)         ! auxiliary, A*p(k)
        real(r8), allocatable :: ci2d_r(:,:)         ! residual,  r(k)
        real(r8), allocatable :: ci2d_x(:,:)         ! estimate,  x(k)

# ifdef ADJUST_BOUNDARY
!
!  Chebyshev iterations (CI) polynomial-based solver horizontal vectors,
!  lateral 1D boundary adjustment in control vector, [LBij:UBij].
!
        real(r8), allocatable :: ciB1d_p(:)          ! direction, p(k)
        real(r8), allocatable :: ciB1d_q(:)          ! auxiliary, A*p(k)
        real(r8), allocatable :: ciB1d_r(:)          ! residual,  r(k)
        real(r8), allocatable :: ciB1d_x(:)          ! estimate,  x(k)

#  ifdef SOLVE3D
!
!  Chebyshev iterations (CI) polynomial-based solver horizontal vectors,
!  lateral 2D boundary adjustment in control vector, [LBij:UBij].
!
        real(r8), allocatable :: ciB2d_p(:,:)        ! direction, p(k)
        real(r8), allocatable :: ciB2d_q(:,:)        ! auxiliary, A*p(k)
        real(r8), allocatable :: ciB2d_r(:,:)        ! residual,  r(k)
        real(r8), allocatable :: ciB2d_x(:,:)        ! estimate,  x(k)
#  endif
# endif

# ifdef SOLVE3D
!
!  Chebyshev iterations (CI) polynomial-based solver horizontal vectors,
!  3D variable in control vector, [LBi:UBi,LBj:UBj,N(ng)].
!
        real(r8), allocatable :: ci3d_p(:,:,:)       ! direction, p(k)
        real(r8), allocatable :: ci3d_q(:,:,:)       ! auxiliary, A*p(k)
        real(r8), allocatable :: ci3d_r(:,:,:)       ! residual,  r(k)
        real(r8), allocatable :: ci3d_x(:,:,:)       ! estimate,  x(k)
# endif

# ifdef READ_SCALES
!
!  Spatially varying isotropic or anisotropic horizontal Daley length
!  scales (in meters) are specified for each correlation function and
!  combined when modeling multiscale background-error covariance matrix
!  effects. These length scales are read from a specified input NetCDF
!  file. The algorithm requires (scaleX, scaleY) values for each control
!  vector variable. Users are responsible for determining appropriate
!  scale values, such as the first Rossby radius of deformation or
!  alternative methodologies.
!
        real(r8), pointer :: bZscaleX(:,:,:)     ! free surface
        real(r8), pointer :: bZscaleY(:,:,:)
!
        real(r8), pointer :: bU2DscaleX(:,:,:)   ! 2D U-momentum
        real(r8), pointer :: bU2DscaleY(:,:,:)
!
        real(r8), pointer :: bV2DscaleX(:,:,:)   ! 2D V-momentum
        real(r8), pointer :: bV2DscaleY(:,:,:)

#  ifdef ADJUST_WSTRESS
!
        real(r8), pointer :: bSUSscaleX(:,:,:)   ! surface U-stress
        real(r8), pointer :: bSUSscaleY(:,:,:)
!
        real(r8), pointer :: bSVSscaleX(:,:,:)   ! surface V-stress
        real(r8), pointer :: bSVSscaleY(:,:,:)
#  endif

#  ifdef SOLVE3D
!
        real(r8), pointer :: bTscaleX(:,:,:,:)   ! tracer variables
        real(r8), pointer :: bTscaleY(:,:,:,:)
!
        real(r8), pointer :: bU3DscaleX(:,:,:)   ! 3D U-momentum
        real(r8), pointer :: bU3DscaleY(:,:,:)
!
        real(r8), pointer :: bV3DscaleX(:,:,:)   ! 3D V-momentum
        real(r8), pointer :: bV3DscaleY(:,:,:)

#   ifdef ADJUST_STFLUX
!
        real(r8), pointer :: bSTFscaleX(:,:,:,:) ! surface tracer flux
        real(r8), pointer :: bSTFscaleY(:,:,:,:)
#   endif
#  endif
# endif
!
        CONTAINS
!
        PROCEDURE :: create        => multiscale_create
        PROCEDURE :: tl_CG_2d      => multiscale_CG_2d_tl
        PROCEDURE :: tl_CI_2d      => multiscale_CI_2d_tl
        PROCEDURE :: ad_CI_2d      => multiscale_CI_2d_ad
# ifdef SOLVE3D
!
        PROCEDURE :: tl_CG_3d      => multiscale_CG_3d_tl
        PROCEDURE :: tl_CI_3d      => multiscale_CI_3d_tl
        PROCEDURE :: ad_CI_3d      => multiscale_CI_3d_ad
!
        PROCEDURE :: tl_Vdiff      => multiscale_Vdiff_tl
        PROCEDURE :: ad_Vdiff      => multiscale_Vdiff_ad
# endif

# ifdef ADJUST_BOUNDARY
!
        PROCEDURE :: tl_CG_b1d     => multiscale_CG_b1d_tl
        PROCEDURE :: tl_CI_b1d     => multiscale_CI_b1d_tl
        PROCEDURE :: ad_CI_b1d     => multiscale_CI_b1d_ad

#  ifdef SOLVE3D
!
        PROCEDURE :: tl_CG_b2d     => multiscale_CG_b2d_tl
        PROCEDURE :: tl_CI_b2d     => multiscale_CI_b2d_tl
        PROCEDURE :: ad_CI_b2d     => multiscale_CI_b2d_ad
!
        PROCEDURE :: tl_bry_Vdiff  => multiscale_bry_Vdiff_tl
        PROCEDURE :: ad_bry_Vdiff  => multiscale_bry_Vdiff_ad
#  endif
# endif
!
      END TYPE multiscale
!
!-----------------------------------------------------------------------
!  Module variables.
!-----------------------------------------------------------------------
!
      TYPE (multiscale), allocatable :: B_ms(:)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      PUBLIC  :: bms_initialize
      PUBLIC  :: bms_finalize
!
# ifdef MULTI_SCALE_DEBUG
      PRIVATE :: multiscale_print
# endif
# ifdef ADJUST_BOUNDARY
      PRIVATE :: dot_prod_1d
# endif
      PRIVATE :: dot_prod_2d
      PRIVATE
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      CONTAINS
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
!  It allocates and initializes the background-error covariance (B)
!  multiscale object.
!
      SUBROUTINE bms_initialize (model)
!
      integer,          intent(in   ) :: model      ! kernel identifier
!
      integer                         :: ng, tile
!
!  Allocate and create B-multiscale object.
!
      IF (.not.allocated(B_ms)) THEN
        allocate ( B_ms(Ngrids) )
      END IF
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL multiscale_create (B_ms(ng), ng, tile, model)
        END DO
      END DO 
!
      RETURN
      END SUBROUTINE bms_initialize
!
!-----------------------------------------------------------------------
!  It deallocates background-error covariance multiscale object.
!
      SUBROUTINE bms_finalize (model)
!
      integer,          intent(in   ) :: model      ! kernel identifier
!
      integer                         :: ng, tile

      RETURN
      END SUBROUTINE bms_finalize
!
! <><><><><><><><><><><><><><><><><><><><><><><> CLASS MULTISCALE <><><>
!
      SUBROUTINE multiscale_create (self, ng, tile, model)
!
      CLASS (multiscale), intent(inout) :: self   ! B-multiscale object
      integer,            intent(in   ) :: ng     ! nested grid
      integer,            intent(in   ) :: tile   ! domain partition
      integer,            intent(in   ) :: model  ! kernel identifier
!
      integer                           :: Nlap, is, ms
      integer                           :: LBi, UBi, LBj, UBj
# ifdef ADJUST_BOUNDARY
      integer                           :: LBij, UBij
# endif
!
      real (r8), parameter              :: IniVal = 0.0_r8
      real (r8)                         :: size2d
!
!  Set tiled arrays bounds.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
# ifdef ADJUST_BOUNDARY
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
# endif
!
      size2d=REAL((UBi-LBi+1)*(UBj-LBj+1),r8)
!
!  Initiziale CLASS parameters.
!
      self%tile=tile
!
      self%LBi=LBi
      self%UBi=UBi
      self%LBj=LBj
      self%UBj=LBj
# ifdef ADJUST_BOUNDARY
      self%LBij=LBij
      self%UBij=UBij
# endif
!
      self%Nscale=Nscale(ng)
!
!  Allocate CLASS parameters.
!
      IF (.not.allocated(self % Mlap)) THEN
        allocate ( self % Mlap(MstateVar,Nscale(ng)) )
        DO ms=1,Nscale(ng)
          DO is=1,MstateVar
            self % Mlap(is,ms) = Mlap(is,ms,ng)
          END DO
        END DO
      END IF
!
      IF (.not.allocated(self % NiterCG)) THEN
        allocate ( self % NiterCG(Nscale(ng)) )
        DO ms=1,Nscale(ng)
          self % NiterCG(ms) = NiterCG(ms,ng)
        END DO
      END IF
!
      IF (.not.allocated(self % NiterCI)) THEN
        allocate ( self % NiterCI(Nscale(ng)) )
        DO ms=1,Nscale(ng)
          self % NiterCI(ms) = NiterCI(ms,ng)
        END DO
      END IF
!
      IF (.not.allocated(self % Bwgt)) THEN
        allocate ( self % Bwgt(MstateVar,Nscale(ng)) )
        DO ms=1,Nscale(ng)
          DO is=1,MstateVar
            self % Bwgt(is,ms) = Bwgt(is,ms,ng)
          END DO
        END DO
      END IF
!
      Nlap=MAXVAL(self % Mlap)
!
!  Allocate the maximum and minimum eigenvalues of the K-Laplacian.
!
      IF (.not.associated(self % eigMaxR2D)) THEN
        allocate ( self % eigMaxR2D(Nlap,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*Nscale(ng),r8)
        self % eigMaxR2d = IniVal
      END IF
      IF (.not.associated(self % eigMinR2D)) THEN
        allocate ( self % eigMinR2D(Nlap,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*Nscale(ng),r8)
        self % eigMinR2D = IniVal
      END IF
!
      IF (.not.associated(self % eigMaxU2D)) THEN
        allocate ( self % eigMaxU2D(Nlap,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*Nscale(ng),r8)
        self % eigMaxU2D = IniVal
      END IF
      IF (.not.associated(self % eigMinU2D)) THEN
        allocate ( self % eigMinU2D(Nlap,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*Nscale(ng),r8)
        self % eigMinU2D = IniVal
      END IF
!
      IF (.not.associated(self % eigMaxV2D)) THEN
        allocate ( self % eigMaxV2D(Nlap,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*Nscale(ng),r8)
        self % eigMaxV2D = IniVal
      END IF
      IF (.not.associated(self % eigMinV2D)) THEN
        allocate ( self % eigMinV2D(Nlap,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*Nscale(ng),r8)
        self % eigMinV2D = IniVal
      END IF

# ifdef ADJUST_WSTRESS
!
      IF (.not.associated(self % eigMaxSUS)) THEN
        allocate ( self % eigMaxSUS(Nlap,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*Nscale(ng),r8)
        self % eigMaxSUS = IniVal
      END IF
      IF (.not.associated(self % eigMinSUS)) THEN
        allocate ( self % eigMinSUS(Nlap,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*Nscale(ng),r8)
        self % eigMinSUS = IniVal
      END IF
!
      IF (.not.associated(self % eigMaxSVS)) THEN
        allocate ( self % eigMaxSVS(Nlap,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*Nscale(ng),r8)
        self % eigMaxSVS = IniVal
      END IF
      IF (.not.associated(self % eigMinSVS)) THEN
        allocate ( self % eigMinSVS(Nlap,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*Nscale(ng),r8)
        self % eigMinSVS = IniVal
      END IF
# endif

# ifdef SOLVE3D
!
      IF (.not.associated(self % eigMaxR3D)) THEN
        allocate ( self % eigMaxR3D(Nlap,N(ng),NT(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*NT(ng)*Nscale(ng),r8)
        self % eigMaxR3D = IniVal
      END IF
      IF (.not.associated(self % eigMinR3D)) THEN
        allocate ( self % eigMinR3D(Nlap,N(ng),NT(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*NT(ng)*Nscale(ng),r8)
        self % eigMinR3D = IniVal
      END IF
!
      IF (.not.associated(self % eigMaxU3D)) THEN
        allocate ( self % eigMaxU3D(Nlap,N(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*Nscale(ng),r8)
        self % eigMaxU3D = IniVal
      END IF
      IF (.not.associated(self % eigMinU3D)) THEN
        allocate ( self % eigMinU3D(Nlap,N(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*Nscale(ng),r8)
        self % eigMinU3D = IniVal
      END IF
!
      IF (.not.associated(self % eigMaxV3D)) THEN
        allocate ( self % eigMaxV3D(Nlap,N(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*Nscale(ng),r8)
        self % eigMaxV3D = IniVal
      END IF
      IF (.not.associated(self % eigMinV3D)) THEN
        allocate ( self % eigMinV3D(Nlap,N(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*Nscale(ng),r8)
        self % eigMinV3D = IniVal
      END IF

#  ifdef ADJUST_STFLUX
!
      IF (.not.associated(self % eigMaxSTF)) THEN
        allocate ( self % eigMaxSTF(Nlap,NT(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*NT(ng)*Nscale(ng),r8)
        self % eigMaxSTF = IniVal
      END IF
      IF (.not.associated(self % eigMinSTF)) THEN
        allocate ( self % eigMinSTF(Nlap,NT(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*NT(ng)*Nscale(ng),r8)
        self % eigMinSTF = IniVal
      END IF
#  endif
# endif

# ifdef ADJUST_BOUNDARY
!
      IF (.not.associated(self % eigMaxR1D)) THEN
        allocate ( self % eigMaxR1D(Nlap,4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*4*Nscale(ng),r8)
        self % eigMaxR1D = IniVal
      END IF
      IF (.not.associated(self % eigMinR1D)) THEN
        allocate ( self % eigMinR1D(Nlap,4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*4*Nscale(ng),r8)
        self % eigMinR1D = IniVal
      END IF
!
      IF (.not.associated(self % eigMaxU1D)) THEN
        allocate ( self % eigMaxU1D(Nlap,4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*4*Nscale(ng),r8)
        self % eigMaxU1D = IniVal
      END IF
      IF (.not.associated(self % eigMinU1D)) THEN
        allocate ( self % eigMinU1D(Nlap,4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*4*Nscale(ng),r8)
        self % eigMinU1D = IniVal
      END IF
!
      IF (.not.associated(self % eigMaxV1D)) THEN
        allocate ( self % eigMaxV1D(Nlap,4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*4*Nscale(ng),r8)
        self % eigMaxV1D = IniVal
      END IF
      IF (.not.associated(self % eigMinV1D)) THEN
        allocate ( self % eigMinV1D(Nlap,4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*4*Nscale(ng),r8)
        self % eigMinV1D = IniVal
      END IF

#  ifdef SOLVE3D
!
      IF (.not.associated(self % eigMaxR1DZ)) THEN
        allocate ( self % eigMaxR1DZ(Nlap,N(ng),NT(ng),4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*NT(ng)*4*Nscale(ng),r8)
        self % eigMaxR1DZ = IniVal
      END IF
      IF (.not.associated(self % eigMinR1DZ)) THEN
        allocate ( self % eigMinR1DZ(Nlap,N(ng),NT(ng),4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*NT(ng)*4*Nscale(ng),r8)
        self % eigMinR1DZ = IniVal
      END IF
!
      IF (.not.associated(self % eigMaxU1DZ)) THEN
        allocate ( self % eigMaxU1DZ(Nlap,N(ng),4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*4*Nscale(ng),r8)
        self % eigMaxU1DZ = IniVal
      END IF
      IF (.not.associated(self % eigMinU1DZ)) THEN
        allocate ( self % eigMinU1DZ(Nlap,N(ng),4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*4*Nscale(ng),r8)
        self % eigMinU1DZ = IniVal
      END IF
!
      IF (.not.associated(self % eigMaxV1DZ)) THEN
        allocate ( self % eigMaxV1DZ(Nlap,N(ng),4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*4*Nscale(ng),r8)
        self % eigMaxV1DZ = IniVal
      END IF
      IF (.not.associated(self % eigMinV1DZ)) THEN
        allocate ( self % eigMinV1DZ(Nlap,N(ng),4,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+REAL(Nlap*N(ng)*4*Nscale(ng),r8)
        self % eigMinV1DZ = IniVal
      END IF
#  endif
# endif
!
!  K-Laplacian Conjugate Gradient (CG) iterative solver horizontal
!  vectors.
!
      IF (.not.allocated(self % cg2d_p)) THEN
        allocate ( self % cg2d_p(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % cg2d_p = IniVal
      END IF
      IF (.not.allocated(self % cg2d_q)) THEN
        allocate ( self % cg2d_q(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % cg2d_q = IniVal
      END IF
      IF (.not.allocated(self % cg2d_r)) THEN
        allocate ( self % cg2d_r(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % cg2d_r = IniVal
      END IF
      IF (.not.allocated(self % cg2d_x)) THEN
        allocate ( self % cg2d_x(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % cg2d_x = IniVal
      END IF

# ifdef ADJUST_BOUNDARY
!
      IF (.not.allocated(self % cgB1d_p)) THEN
        allocate ( self % cgB1d_p(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % cgB1d_p = IniVal
      END IF
      IF (.not.allocated(self % cgB1d_q)) THEN
        allocate ( self % cgB1d_q(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % cgB1d_q = IniVal
      END IF
      IF (.not.allocated(self % cgB1d_r)) THEN
        allocate ( self % cgB1d_r(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % cgB1d_r = IniVal
      END IF
      IF (.not.allocated(self % cgB1d_x)) THEN
        allocate ( self % cgB1d_x(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % cgB1d_x = IniVal
      END IF

#  ifdef SOLVE3D
!
      IF (.not.allocated(self % cgB2d_p)) THEN
        allocate ( self % cgB2d_p(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % cgB2d_p = IniVal
      END IF
      IF (.not.allocated(self % cgB2d_q)) THEN
        allocate ( self % cgB2d_q(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % cgB2d_q = IniVal
      END IF
      IF (.not.allocated(self % cgB2d_r)) THEN
        allocate ( self % cgB2d_r(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % cgB2d_r = IniVal
      END IF
      IF (.not.allocated(self % cgB2d_x)) THEN
        allocate ( self % cgB2d_x(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % cgB2d_x = IniVal
      END IF
#  endif
# endif

# ifdef SOLVE3D
!
      IF (.not.allocated(self % cg3d_p)) THEN
        allocate ( self % cg3d_p(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % cg3d_p = IniVal
      END IF
      IF (.not.allocated(self % cg3d_q)) THEN
        allocate ( self % cg3d_q(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % cg3d_q = IniVal
      END IF
      IF (.not.allocated(self % cg3d_r)) THEN
        allocate ( self % cg3d_r(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % cg3d_r = IniVal
      END IF
      IF (.not.allocated(self % cg3d_x)) THEN
        allocate ( self % cg3d_x(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % cg3d_x = IniVal
      END IF
# endif
!
!  Chebyshev iterations (CI) polynomial-based solver horizontal vectors.
!
      IF (.not.allocated(self % ci2d_p)) THEN
        allocate ( self % ci2d_p(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % ci2d_p = IniVal
      END IF
      IF (.not.allocated(self % ci2d_q)) THEN
        allocate ( self % ci2d_q(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % ci2d_q = IniVal
      END IF
      IF (.not.allocated(self % ci2d_r)) THEN
        allocate ( self % ci2d_r(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % ci2d_r = IniVal
      END IF
      IF (.not.allocated(self % ci2d_x)) THEN
        allocate ( self % ci2d_x(LBi:UBi,LBj:UBj) )
        Dmem(ng)=Dmem(ng)+size2d
        self % ci2d_x = IniVal
      END IF

# ifdef ADJUST_BOUNDARY
!
      IF (.not.allocated(self % ciB1d_p)) THEN
        allocate ( self % ciB1d_p(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % ciB1d_p = IniVal
      END IF
      IF (.not.allocated(self % ciB1d_q)) THEN
        allocate ( self % ciB1d_q(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % ciB1d_q = IniVal
      END IF
      IF (.not.allocated(self % ciB1d_r)) THEN
        allocate ( self % ciB1d_r(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % ciB1d_r = IniVal
      END IF
      IF (.not.allocated(self % ciB1d_x)) THEN
        allocate ( self % ciB1d_x(LBij:UBij) )
        Dmem(ng)=Dmem(ng)+REAL(UBij-LBij+1,r8)
        self % ciB1d_x = IniVal
      END IF

#  ifdef SOLVE3D
!
      IF (.not.allocated(self % ciB2d_p)) THEN
        allocate ( self % ciB2d_p(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % ciB2d_p = IniVal
      END IF
      IF (.not.allocated(self % ciB2d_q)) THEN
        allocate ( self % ciB2d_q(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % ciB2d_q = IniVal
      END IF
      IF (.not.allocated(self % ciB2d_r)) THEN
        allocate ( self % ciB2d_r(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % ciB2d_r = IniVal
      END IF
      IF (.not.allocated(self % ciB2d_x)) THEN
        allocate ( self % ciB2d_x(LBij:UBij,N(ng)) )
        Dmem(ng)=Dmem(ng)+REAL((UBij-LBij+1)*N(ng),r8)
        self % ciB2d_x = IniVal
      END IF
#  endif
# endif

# ifdef SOLVE3D
!
      IF (.not.allocated(self % ci3d_p)) THEN
        allocate ( self % ci3d_p(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % ci3d_p = IniVal
      END IF
      IF (.not.allocated(self % ci3d_q)) THEN
        allocate ( self % ci3d_q(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % ci3d_q = IniVal
      END IF
      IF (.not.allocated(self % ci3d_r)) THEN
        allocate ( self % ci3d_r(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % ci3d_r = IniVal
      END IF
      IF (.not.allocated(self % ci3d_x)) THEN
        allocate ( self % ci3d_x(LBi:UBi,LBj:UBj,N(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng),r8)
        self % ci3d_x = IniVal
      END IF
# endif

# ifdef READ_SCALES
!
!  Spatially varying isotropic or anisotropic horizontal Daley length
!  scales.
!
      IF (.not.associated(self % bZscaleX)) THEN
        allocate ( self % bZscaleX(LBi:UBi,LBj:UBj,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(Nscale(ng),r8)
        self % bZscaleX = IniVal
      END IF
      IF (.not.associated(self % bZscaleY)) THEN
        allocate ( self % bZscaleY(LBi:UBi,LBj:UBj,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(Nscale(ng),r8)
        self % bZscaleY = IniVal
      END IF
!
      IF (.not.associated(self % bU2DscaleX)) THEN
        allocate ( self % bU2DscaleX(LBi:UBi,LBj:UBj,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(Nscale(ng),r8)
        self % bU2DscaleX = IniVal
      END IF
      IF (.not.associated(self % bU2DscaleY)) THEN
        allocate ( self % bU2DscaleY(LBi:UBi,LBj:UBj,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(Nscale(ng),r8)
        self % bU2DscaleY = IniVal
      END IF
!
      IF (.not.associated(self % bV2DscaleX)) THEN
        allocate ( self % bV2DscaleX(LBi:UBi,LBj:UBj,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(Nscale(ng),r8)
        self % bV2DscaleX = IniVal
      END IF
      IF (.not.associated(self % bV2DscaleY)) THEN
        allocate ( self % bV2DscaleY(LBi:UBi,LBj:UBj,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(Nscale(ng),r8)
        self % bV2DscaleY = IniVal
      END IF

#  ifdef ADJUST_WSTRESS
!
      IF (.not.associated(self % bSUSscaleX)) THEN
        allocate ( self % bSUSscaleX(LBi:UBi,LBj:UBj,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(Nscale(ng),r8)
        self % bSUSscaleX = IniVal
      END IF
      IF (.not.associated(self % bSUSscaleY)) THEN
        allocate ( self % bSUSscaleY(LBi:UBi,LBj:UBj,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(Nscale(ng),r8)
        self % bSUSscaleY = IniVal
      END IF
!
      IF (.not.associated(self % bSVSscaleX)) THEN
        allocate ( self % bSVSscaleX(LBi:UBi,LBj:UBj,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(Nscale(ng),r8)
        self % bSVSscaleX = IniVal
      END IF
      IF (.not.associated(self % bSVSscaleY)) THEN
        allocate ( self % bSVSscaleY(LBi:UBi,LBj:UBj,Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(Nscale(ng),r8)
        self % bSVSscaleY = IniVal
      END IF
#  endif

#  ifdef SOLVE3D
!
      IF (.not.associated(self % bTscaleX)) THEN
        allocate ( self % bTscaleX(LBi:UBi,LBj:UBj,                     &
     &                             N(ng),NT(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng)*NT(ng)*Nscale(ng),r8)
        self % bTscaleX = IniVal
      END IF
      IF (.not.associated(self % bTscaleY)) THEN
        allocate ( self % bTscaleY(LBi:UBi,LBj:UBj,                     &
     &                             N(ng),NT(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng)*NT(ng)*Nscale(ng),r8)
        self % bTscaleY = IniVal
      END IF
!
      IF (.not.associated(self % bU3DscaleX)) THEN
        allocate ( self % bU3DscaleX(LBi:UBi,LBj:UBj,N(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng)*Nscale(ng),r8)
        self % bU3DscaleX = IniVal
      END IF
      IF (.not.associated(self % bU3DscaleY)) THEN
        allocate ( self % bU3DscaleY(LBi:UBi,LBj:UBj,N(ng)*Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng)*Nscale(ng),r8)
        self % bU3DscaleY = IniVal
      END IF
!
      IF (.not.associated(self % bV3DscaleX)) THEN
        allocate ( self % bV3DscaleX(LBi:UBi,LBj:UBj,N(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng)*Nscale(ng),r8)
        self % bV3DscaleX = IniVal
      END IF
      IF (.not.associated(self % bV3DscaleY)) THEN
        allocate ( self % bV3DscaleY(LBi:UBi,LBj:UBj,N(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(N(ng)*Nscale(ng),r8)
        self % bV3DscaleY = IniVal
      END IF

#   ifdef ADJUST_STFLUX
!
      IF (.not.associated(self % bSTFscaleX)) THEN
        allocate ( self % bSTFscaleX(LBi:UBi,LBj:UBj,                   &
     &                               NT(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(NT(ng)*Nscale(ng),r8)
        self % bSTFscaleX = IniVal
      END IF
      IF (.not.associated(self % bSTFscaleY)) THEN
        allocate ( self % bSTFscaleY(LBi:UBi,LBj:UBj,                   &
     &                               NT(ng),Nscale(ng)) )
        Dmem(ng)=Dmem(ng)+size2d*REAL(NT(ng)*Nscale(ng),r8)
        self % bSTFscaleY = IniVal
      END IF
#   endif
#  endif
# endif
!
      RETURN
      END SUBROUTINE multiscale_create
!
!-----------------------------------------------------------------------
!  This routine computes the extreme eigenvalues required by the
!  Chebyshev Iterations (CI) solver, which is applied to implicit
!  diffusion operators in the modeling of the multiscale background-
!  error covariance for 2D variables in the control vector.
!
!  The minimum and maximum eigenvalues of the K-Laplacian operator
!  are determined using a Conjugate Gradient (CG) solution to the
!  corresponding linear system, A x = b.
!
      SUBROUTINE multiscale_CG_2d_tl (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCG, Lweak,        &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCG   ! CG iterations
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCG, iterDiff, info
      integer                           :: i, j
# ifdef MULTI_SCALE_DEBUG
      integer                           :: status
# endif
!
      real (r8)                         :: dotn, dotpq, dotr, dotrd
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_scale
      real (r8), dimension(0:NiterCG)               :: cg_a
      real (r8), dimension(0:NiterCG+1)             :: cg_b
      real (r8), dimension(NiterCG+1)               :: cg_Rv
      real (r8), dimension(2*(NiterCG+1)-2)         :: work
      real (r8), dimension(NiterCG)                 :: zwork
      real (r8), dimension(NiterCG+1,NiterCG+1)     :: zgv
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CG_2d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      self%cg2d_r=0.0_r8
      self%cg2d_p=0.0_r8
      self%cg2d_q=0.0_r8
      self%cg2d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')                          ! free surface
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR2D
          eigMax => self%eigMaxR2D
        CASE ('ubar', 'ubar_eastward')         ! 2D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU2D
          eigMax => self%eigMaxU2D
        CASE ('vbar', 'vbar_northward')        ! 2D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV2D
          eigMax => self%eigMaxV2D
        CASE ('sustr')                         ! surface U-stress
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSUS
          eigMax => self%eigMaxSUS
        CASE ('svstr')                         ! surface V-stress
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSVS
          eigMax => self%eigMaxSVS
# ifdef SOLVE3D
        CASE ('shflux')                        ! surface net heat flux
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSTF(:,:,itemp)
          eigMax => self%eigMaxSTF(:,:,itemp)
        CASE ('ssflux')                        ! surface net salt flux
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSTF(:,:,isalt)
          eigMax => self%eigMaxSTF(:,:,isalt)
# endif
      END SELECT
!
!  Set control variable squared root area scale (2D).
!
      SELECT CASE (ctype)
        CASE (r2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
# ifdef MASKING
              tl_A(i,j)=GRID(ng)%rmask(i,j)*tl_A(i,j)
# endif
            END DO
          END DO
        CASE (u2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j))
# ifdef MASKING
              tl_A(i,j)=GRID(ng)%umask(i,j)*tl_A(i,j)
# endif
            END DO
          END DO
        CASE (v2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j))
# ifdef MASKING
              tl_A(i,j)=GRID(ng)%vmask(i,j)*tl_A(i,j)
# endif
            END DO
          END DO
      END SELECT
!
!  Advance in K-space implicit pseudo-diffusion 2D equation using 
!  Weaver et al. (2016) algorithm 1.
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            self%cg2d_x(i,j)=0.0_r8                       ! step 1
            self%cg2d_r(i,j)=-tl_A(i,j)                   ! step 2
            self%cg2d_p(i,j)=tl_A(i,j)                    ! step 3                     
          END DO
        END DO
        cg_b(0)=0.0_r8                                    ! step 4
!
        dotn = dot_prod_2d (ng, tile, model, ctype,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      self%cg2d_r,                                &
     &                      self%cg2d_r)
!
!  Invert the implicit diffusion equation using CG iterations.
!
!     A x = b  linear system
!
        KLAP_ITER : DO iterCG=0,NiterCG                   ! steps 6 - 16
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j)=tl_A(i,j)/tl_scale(i,j)
            END DO
          END DO
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
          SELECT CASE (ctype)
            CASE (r2dvar)
              CALL tl_multiscale_Klap_r2d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
            CASE (u2dvar)
              CALL tl_multiscale_Klap_u2d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
            CASE (v2dvar)
              CALL tl_multiscale_Klap_v2d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
          END SELECT
!
!  Conjugate Gradient (CG) algorithm, (Weaver et al., 2016).
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j)=tl_A(i,j)*tl_scale(i,j)
              self%cg2d_q(i,j)=tl_A(i,j)                  ! step 7
            END DO
          END DO
!
          dotr = dot_prod_2d (ng, tile, model, ctype,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        self%cg2d_r,                              &
     &                        self%cg2d_r)
          dotpq = dot_prod_2d (ng, tile, model, ctype,                  &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         self%cg2d_p,                             &
     &                         self%cg2d_q)
!
          cg_a(iterCG)=dotr/dotpq                         ! step 8

# ifdef MULTI_SCALE_DEBUG
!
          deps=SQRT(dotr/dotn)
          status=multiscale_print(ng, ifield, ms, iterDiff, iterCG,     & 
     &                            NiterCG, 'CG_2d', deps)
# endif
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%cg2d_x(i,j)=self%cg2d_x(i,j)+                        &
     &                         cg_a(iterCG)*                            &
     &                         self%cg2d_p(i,j)           ! step 9
              self%cg2d_r(i,j)=self%cg2d_r(i,j)+                        &
     &                         cg_a(iterCG)*                            &
     &                         self%cg2d_q(i,j)           ! step 10
            END DO
          END DO
!
!  Reorthogonalize r(k+1) with respect to r(0:k).           step 11
!
          dotrd=dotr
          dotr = dot_prod_2d (ng, tile, model, ctype,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        self%cg2d_r,                              &
     &                        self%cg2d_r)
!
          cg_b(iterCG+1)=dotr/dotrd                       ! step 12
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%cg2d_p(i,j)=-self%cg2d_r(i,j)+                       &
     &                         cg_b(iterCG+1)*                          &
     &                         self%cg2d_p(i,j)           ! step 13
              tl_A(i,j)=self%cg2d_p(i,j)
            END DO
          END DO

        END DO KLAP_ITER
!
!  Update control vector 2D variable with the implicit diffusion
!  solution iterate.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            tl_A(i,j)=self%cg2d_x(i,j)
          END DO
        END DO
!
!  Use the LAPACK routine DSTEQR to compute the Ritz vectors and 
!  Ritz values of the tridiagonal matrix.
!
        DO i=1,NiterCG+1                                  ! step 14
          cg_Rv(i)=1.0_r8/cg_a(i-1)
          IF (i.gt.1) THEN
            cg_Rv(i)=cg_Rv(i)+cg_b(i-1)/cg_a(i-2)
          END IF
        END DO
        DO i=1,NiterCG
          zwork(i)=-SQRT(cg_b(i))/cg_a(i-1)
        END DO
!
        CALL dsteqr ('I', NiterCG+1, cg_Rv, zwork,                      &
     &               zgv, NiterCG+1, work, info)
        IF (info.ne.0) THEN
          WRITE (stdout,*) ' multiscale_CG_2d_tl - Error in DSTEQR:',   &
     &                     ' info = ', info
          exit_flag=8
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        eigMax(iterDiff,ms)=MAXVAL(cg_Rv)  ! maximum Ritz eigenvalue
        eigMin(iterDiff,ms)=MINVAL(cg_Rv)  ! minimum Ritz eigenvalue
!
!  Reset iterative conjugate gradient arrays.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            self%cg2d_p(i,j)=0.0_r8
            self%cg2d_r(i,j)=0.0_r8
            self%cg2d_q(i,j)=0.0_r8
          END DO
        END DO
!
      END DO KDIFF_ITER
!
!  Convert the solution back to physical space.
!  (cf Eqn 18 of Weaver et al 2016, QJRMS, 143,455-471).
!  Not needed in the CI algorithm as the tl_scale factors cancel.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          tl_A(i,j)=self%cg2d_x(i,j)/tl_scale(i,j)
        END DO
      END DO
!
      RETURN
      END SUBROUTINE multiscale_CG_2d_tl
!
!-----------------------------------------------------------------------
!  This routine inverts the implicit diffusion linear operator, Ax = b,
!  using Chebyshev Iterations (CI) to model the multiscale propagation
!  of background-error covariance for 2D variables in the control
!  vector. It operates in conjunction with its adjoint to enforce
!  symmetry-preserving correlation functions. The algorithm requires
!  the extreme eigenvalues of matrix A, which are obtained from the
!  Conjugate Gradient (CG) solver. Additionally, initial values for
!  the estimate x(0) and its residual r(0) are required.
!
      SUBROUTINE multiscale_CI_2d_tl (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCI, Lweak,        &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CI iterations
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j
# ifdef MULTI_SCALE_DEBUG
      integer                           :: status
# endif
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_scale
      real (r8), dimension(0:NiterCI)               :: ci_alpha
      real (r8), dimension(0:NiterCI+1)             :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_2d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      self%ci2d_r=0.0_r8
      self%ci2d_p=0.0_r8
      self%ci2d_q=0.0_r8
      self%ci2d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms). Also, select pre-computed extreme
!  eigenvalues of matrix A.
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')                          ! free surface
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR2D(:,ms)
          eigMax => self%eigMaxR2D(:,ms)
        CASE ('ubar', 'ubar_eastward')         ! 2D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU2D(:,ms)
          eigMax => self%eigMaxU2D(:,ms)
        CASE ('vbar', 'vbar_northward')        ! 2D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV2D(:,ms)
          eigMax => self%eigMaxV2D(:,ms)
        CASE ('sustr')                         ! surface U-stress
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSUS(:,ms)
          eigMax => self%eigMaxSUS(:,ms)
        CASE ('svstr')                         ! surface V-stress
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSVS(:,ms)
          eigMax => self%eigMaxSVS(:,ms)
# ifdef SOLVE3D
        CASE ('shflux')                        ! surface net heat flux
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSTF(:,itemp,ms)
          eigMax => self%eigMaxSTF(:,itemp,ms)
        CASE ('ssflux')                        ! surface net salt flux
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSTF(:,isalt,ms)
          eigMax => self%eigMaxSTF(:,isalt,ms)
# endif
      END SELECT
!
!  Set control variable squared root area scale (2D).
!
      SELECT CASE (ctype)
        CASE (r2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
            END DO
          END DO
        CASE (u2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j))
            END DO
          END DO
        CASE (v2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j))
            END DO
          END DO
      END SELECT
!
!  Advance in K-space implicit pseudo-diffusion equation using 
!  Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        ci_sigma=0.5_r8*(eigMax(iterDiff)+                              &
     &                   eigMin(iterDiff))              ! step 1
        ci_delta=0.5_r8*(eigMax(iterDiff)-                              &
     &                   eigMin(iterDiff))              ! step 2
        ci_alpham1=1.0_r8                               ! step 3
        ci_beta(0)=0.0_r8                               ! step 4
!
!  Invert the implicit diffusion equation using Chebyshev Iterations.
!
!     A x = b  linear system
!
        DO iterCI=0,NiterCI                             ! steps 5 to 8
          IF (iterCI.eq.0) THEN
            ci_alpha(0)=1.0_r8/ci_sigma
            cff=ci_delta*ci_alpha(0)
            ci_beta(1)=0.5_r8*cff*cff
          ELSE
            ci_alpha(iterCI)=1.0_r8/(ci_sigma-                          &
     &                               ci_beta(iterCI)/                   &
     &                               ci_alpha(iterCI-1))
            cff=0.5_r8*ci_delta*ci_alpha(iterCI)
            ci_beta(iterCI+1)=cff*cff
          END IF
        END DO
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            tl_A(i,j)=tl_scale(i,j)*tl_A(i,j)
            self%ci2d_x(i,j)=0.0_r8                     ! step 9
            self%ci2d_r(i,j)=-tl_A(i,j)                 ! step 10
            self%ci2d_p(i,j)=tl_A(i,j)                  ! step 11
          END DO
        END DO

# ifdef MULTI_SCALE_DEBUG
!
        dotn = dot_prod_2d (ng, tile, model, ctype,                     &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      self%ci2d_r, self%ci2d_r)
# endif
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=0,NiterCI                 ! steps 12 to 17
!
          DO j=Jstr,Jend
            DO i=Istr,Iend
              tl_A(i,j)=tl_A(i,j)/tl_scale(i,j)
            END DO
          END DO
!
          SELECT CASE (ctype)                           ! step 12
            CASE (r2dvar)
              CALL tl_multiscale_Klap_r2d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
            CASE (u2dvar)
              CALL tl_multiscale_Klap_u2d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
            CASE (v2dvar)
              CALL tl_multiscale_Klap_v2d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
          END SELECT
!                                                        step 13
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j)=tl_A(i,j)*tl_scale(i,j)
              self%ci2d_q(i,j)=tl_A(i,j)
            END DO
          END DO

# ifdef MULTI_SCALE_DEBUG
!
          dotr = dot_prod_2d (ng, tile, model, ctype,                   &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        self%ci2d_r, self%ci2d_r)
          deps=SQRT(dotr/dotn)
!
          status=multiscale_print(ng, ifield, ms, iterDiff, iterCI,     & 
     &                            NiterCI, 'CI_2d', deps)
# endif
!                                                        steps 14, 15
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%ci2d_x(i,j)=self%ci2d_x(i,j)+                        &
     &                         ci_alpha(iterCI)*self%ci2d_p(i,j)
              self%ci2d_r(i,j)=self%ci2d_r(i,j)+                        &
     &                         ci_alpha(iterCI)*self%ci2d_q(i,j)
            END DO
          END DO
!                                                        step 16
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%ci2d_p(i,j)=-self%ci2d_r(i,j)+                       &
     &                         ci_beta(iterCI+1)*self%ci2d_p(i,j)
              tl_A(i,j)=self%ci2d_p(i,j)
            END DO
          END DO

        END DO KLAP_ITER
!
!  Reset iterative conjugate gradient arrays.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            self%ci2d_p(i,j)=0.0_r8
            self%ci2d_r(i,j)=0.0_r8
            self%ci2d_q(i,j)=0.0_r8
          END DO
        END DO
!
!  Update control vector 2D variable with the implicit diffusion
!  solution iterate.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            tl_A(i,j)=self%ci2d_x(i,j)/tl_scale(i,j)
          END DO
        END DO
!
      END DO KDIFF_ITER
!
      RETURN
      END SUBROUTINE multiscale_CI_2d_tl
!
!-----------------------------------------------------------------------
!  This routine serves as the adjoint to the inversion of the implicit
!  diffusion linear operator, Ax = b, utilizing Chebyshev Iterations
!  (CI) to model the multiscale propagation of background-error
!  covariance for 2D variables in the control vector. It is needed to
!  enforce symmetry-preserving correlation functions. The algorithm
!  requires the extreme eigenvalues of matrix A, which are obtained
!  from the Conjugate Gradient (CG) solver.
!
      SUBROUTINE multiscale_CI_2d_ad (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCI, Lweak,        &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CI iterations
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_scale
      real (r8), dimension(0:NiterCI)               :: ci_alpha
      real (r8), dimension(0:NiterCI+1)             :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_2d_ad"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      self%ci2d_r=0.0_r8
      self%ci2d_p=0.0_r8
      self%ci2d_q=0.0_r8
      self%ci2d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms). Also, select pre-computed extreme
!  eigenvalues of matrix A.
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')                          ! free surface
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR2D(:,ms)
          eigMax => self%eigMaxR2D(:,ms)
        CASE ('ubar', 'ubar_eastward')         ! 2D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU2D(:,ms)
          eigMax => self%eigMaxU2D(:,ms)
        CASE ('vbar', 'vbar_northward')        ! 2D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV2D(:,ms)
          eigMax => self%eigMaxV2D(:,ms)
        CASE ('sustr')                         ! surface U-stress
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSUS(:,ms)
          eigMax => self%eigMaxSUS(:,ms)
        CASE ('svstr')                         ! surface V-stress
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSVS(:,ms)
          eigMax => self%eigMaxSVS(:,ms)
# ifdef SOLVE3D
        CASE ('shflux')                        ! surface net heat flux
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSTF(:,itemp,ms)
          eigMax => self%eigMaxSTF(:,itemp,ms)
        CASE ('ssflux')                        ! surface net salt flux
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinSTF(:,isalt,ms)
          eigMax => self%eigMaxSTF(:,isalt,ms)
# endif
      END SELECT
!
!  Set control variable squared root area scale (2D).
!
      SELECT CASE (ctype)
        CASE (r2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
            END DO
          END DO
        CASE (u2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j))
            END DO
          END DO
        CASE (v2dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j))
            END DO
          END DO
      END SELECT
!
!  Adjoint of advance in K-space implicit pseudo-diffusion equation
!  using  Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=Mlap,1,-1
!
        ci_sigma=0.5_r8*(eigMax(iterDiff)+                              &
     &                   eigMin(iterDiff))              ! step 1
        ci_delta=0.5_r8*(eigMax(iterDiff)-                              &
     &                   eigMin(iterDiff))              ! step 2
        ci_alpham1=1.0_r8                               ! step 3
        ci_beta(0)=0.0_r8                               ! step 4
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            self%ci2d_p(i,j)=0.0_r8
            self%ci2d_r(i,j)=0.0_r8
            self%ci2d_q(i,j)=0.0_r8
            self%ci2d_x(i,j)=0.0_r8
          END DO
        END DO
!
!  Adjoint of invert the implicit diffusion equation using Chebyshev
!  Iterations.
!
!     A x = b  linear system
!
        DO iterCI=0,NiterCI                             ! steps 5 to 8
          IF (iterCI.eq.0) THEN
            ci_alpha(0)=1.0_r8/ci_sigma
            cff=ci_delta*ci_alpha(0)
            ci_beta(1)=0.5_r8*cff*cff
          ELSE
            ci_alpha(iterCI)=1.0_r8/(ci_sigma-                          &
     &                               ci_beta(iterCI)/                   &
     &                               ci_alpha(iterCI-1))
            cff=0.5_r8*ci_delta*ci_alpha(iterCI)
            ci_beta(iterCI+1)=cff*cff
          END IF
        END DO
!
!  Adjoint of update control vector 2D variable with the implicit
!  diffusion solution iterate.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
!^          tl_A(i,j)=self%ci2d_x(i,j)/tl_scale(i,j)
!^
            self%ci2d_x(i,j)=ad_A(i,j)/ad_scale(i,j)
            ad_A(i,j)=0.0_r8
          END DO
        END DO
!
!  Adjoint of reset iterative conjugate gradient arrays.
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
!^          self%ci2d_q(i,j)=0.0_r8
!^
            self%ci2d_q(i,j)=0.0_r8
!^          self%ci2d_r(i,j)=0.0_r8
!^
            self%ci2d_r(i,j)=0.0_r8
!^          self%ci2d_p(i,j)=0.0_r8
!^
            self%ci2d_p(i,j)=0.0_r8
          END DO
        END DO
!
!  Adjoint of compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=NiterCI,0,-1              ! steps 17 to 12
!
          DO j=Jmin,Jmax                                ! step 16
            DO i=Imin,Imax
!^            tl_A(i,j)=self%ci2d_p(i,j)
!^
              self%ci2d_p(i,j)=self%ci2d_p(i,j)+ad_A(i,j)
              ad_A(i,j)=0.0_r8
!^            self%ci2d_p(i,j)=-self%ci2d_r(i,j)+                       &
!^   &                         ci_beta(iterCI+1)*self%ci2d_p(i,j)
!^
              self%ci2d_r(i,j)=self%ci2d_r(i,j)-self%ci2d_p(i,j)
              self%ci2d_p(i,j)=ci_beta(iterCI+1)*self%ci2d_p(i,j)
            END DO
          END DO
!                                                        steps 15, 14
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            self%ci2d_r(i,j)=self%ci2d_r(i,j)+                        &
!^   &                         ci_alpha(iterCI)*self%ci2d_q(i,j)
!^
              self%ci2d_q(i,j)=self%ci2d_q(i,j)+                        &
     &                         ci_alpha(iterCI)*self%ci2d_r(i,j)
!^            self%ci2d_x(i,j)=self%ci2d_x(i,j)+                        &
!^   &                         ci_alpha(iterCI)*self%ci2d_p(i,j)
!^
              self%ci2d_p(i,j)=self%ci2d_p(i,j)+                        &
     &                         ci_alpha(iterCI)*self%ci2d_x(i,j)
            END DO
          END DO
!                                                        step 13
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            self%ci2d_q(i,j)=tl_A(i,j)
!^
              ad_A(i,j)=ad_A(i,j)+self%ci2d_q(i,j)
              self%ci2d_q(i,j)=0.0_r8
!^            tl_A(i,j)=tl_A(i,j)*tl_scale(i,j)
!^
              ad_A(i,j)=ad_A(i,j)*ad_scale(i,j)
            END DO
          END DO
!
          SELECT CASE (ctype)                           ! step 12
            CASE (r2dvar)
              CALL ad_multiscale_Klap_r2d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     ad_A)
            CASE (u2dvar)
              CALL ad_multiscale_Klap_u2d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     ad_A)
            CASE (v2dvar)
              CALL ad_multiscale_Klap_v2d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     ad_A)
          END SELECT
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            tl_A(i,j)=tl_A(i,j)/tl_scale(i,j)
!^
              ad_A(i,j)=ad_A(i,j)/ad_scale(i,j)
            END DO
          END DO

        END DO KLAP_ITER
!
        DO j=Jmin,Jmax
          DO i=Imin,Imax
!^          self%cg2d_p(i,j)=tl_A(i,j)                    ! step 3
!^
            ad_A(i,j)=ad_A(i,j)+self%ci2d_p(i,j)
            self%ci2d_p(i,j)=0.0_r8
!^          self%cg2d_r(i,j)=-tl_A(i,j)                   ! step 2
!^
            ad_A(i,j)=ad_A(i,j)-self%ci2d_r(i,j)!
!^          self%cg2d_x(i,j)=0.0_r8                       ! step 1
!^
            self%ci2d_r(i,j)=0.0_r8
          END DO
        END DO

      END DO KDIFF_ITER
!
      RETURN
      END SUBROUTINE multiscale_CI_2d_ad

# ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  This routine computes the extreme eigenvalues required by the
!  Chebyshev Iterations (CI) solver, which is applied to implicit
!  diffusion operators in the modeling of the multiscale background-
!  error covariance for 3D variables in the control vector.
!
!  The minimum and maximum eigenvalues of the K-Laplacian operator
!  are determined using a Conjugate Gradient (CG) solution to the
!  corresponding linear system.
!
      SUBROUTINE multiscale_CG_3d_tl (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCG, Lweak,        &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCG   ! CG iterations
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCG, iterDiff, info
      integer                           :: i, j, k
#  ifdef MULTI_SCALE_DEBUG
      integer                           :: status
#  endif
!
      real (r8)                         :: dotpq, dotr, dotrd
      real (r8), dimension(N(ng))       :: dotn
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_scale
      real (r8), dimension(0:NiterCG,N(ng))         :: cg_a
      real (r8), dimension(0:NiterCG+1,N(ng))       :: cg_b
      real (r8), dimension(NiterCG+1)               :: cg_Rv
      real (r8), dimension(2*(NiterCG+1)-2)         :: work
      real (r8), dimension(NiterCG)                 :: zwork
      real (r8), dimension(NiterCG+1,NiterCG+1)     :: zgv
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CG_3d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u3dvar)
          Imin=IstrU
        CASE (v3dvar)
          Jmin=JstrV
      END SELECT
!
      self%cg3d_r=0.0_r8
      self%cg3d_p=0.0_r8
      self%cg3d_q=0.0_r8
      self%cg3d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward')               ! 3D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU3D(:,:,ms)
          eigMax => self%eigMaxU3D(:,:,ms)
        CASE ('v', 'v_northward')              ! 3D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV3D(:,:,ms)
          eigMax => self%eigMaxV3D(:,:,ms)
        CASE ('temp')                          ! temperature
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR3D(:,:,itemp,ms)
          eigMax => self%eigMaxR3D(:,:,itemp,ms)
        CASE ('salt')                          ! salinity
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR3D(:,:,isalt,ms)
          eigMax => self%eigMaxR3D(:,:,isalt,ms)
      END SELECT
!
!  Set control variable squared root area scale (2D).
!
      SELECT CASE (ctype)
        CASE (r3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*                    &
     &                           GRID(ng)%on_r(i,j))
#  ifdef MASKING
              DO k=1,N(ng)
                tl_A(i,j,k)=tl_A(i,j,k)*GRID(ng)%rmask(i,j)
              END DO
#  endif
            END DO
          END DO
        CASE (u3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*                    &
     &                           GRID(ng)%on_u(i,j))
#  ifdef MASKING
              DO k=1,N(ng)
                tl_A(i,j,k)=tl_A(i,j,k)*GRID(ng)%umask(i,j)
              END DO
#  endif
            END DO
          END DO
        CASE (v3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*                    &
     &                           GRID(ng)%on_v(i,j))
#  ifdef MASKING
              DO k=1,N(ng)
                tl_A(i,j,k)=tl_A(i,j,k)*GRID(ng)%vmask(i,j)
              END DO
#  endif
            END DO
          END DO
      END SELECT
!
!  Advance in K-space implicit pseudo-diffusion 3D equation using 
!  Weaver et al. (2016) algorithm 1.
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        LEVEL_LOOP1 : DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%cg3d_x(i,j,k)=0.0_r8                    ! step 1
              self%cg3d_r(i,j,k)=-tl_A(i,j,k)              ! step 2
              self%cg3d_p(i,j,k)=tl_A(i,j,k)               ! step 3
            END DO
          END DO
          cg_b(0,k)=0.0_r8                                 ! step 4
!
          dotn(k) = dot_prod_2d (ng, tile, model, ctype,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           self%cg3d_r(:,:,k),                    &
     &                           self%cg3d_r(:,:,k))
        END DO LEVEL_LOOP1
!
!  Invert the implicit diffusion equation using CG iterations.
!
!     A x = b  linear system
!
        KLAP_ITER : DO iterCG=0,NiterCG                   ! steps 6 - 16
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                tl_A(i,j,k)=tl_A(i,j,k)/tl_scale(i,j)
              END DO
            END DO
          END DO LEVEL_LOOP2
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
          SELECT CASE (ctype)
            CASE (r3dvar)
              CALL tl_multiscale_Klap_r3d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
            CASE (u3dvar)
              CALL tl_multiscale_Klap_u3d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
            CASE (v3dvar)
              CALL tl_multiscale_Klap_v3d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
          END SELECT
!
!  Conjugate Gradient (CG) algorithm, (Weaver et al., 2016).
!
          LEVEL_LOOP3 : DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                tl_A(i,j,k)=tl_A(i,j,k)*tl_scale(i,j)
                self%cg3d_q(i,j,k)=tl_A(i,j,k)            ! step 7
              END DO
            END DO
!
            dotr = dot_prod_2d (ng, tile, model, ctype,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          self%cg3d_r(:,:,k),                     &
     &                          self%cg3d_r(:,:,k))
            dotpq = dot_prod_2d (ng, tile, model, ctype,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           self%cg3d_p(:,:,k),                    &
     &                           self%cg3d_q(:,:,k))
!
            cg_a(iterCG,k)=dotr/dotpq                     ! step 8

#  ifdef MULTI_SCALE_DEBUG
!
            deps=SQRT(dotr/dotn(k))
            status=multiscale_print(ng, ifield, ms, iterDiff, iterCG,   & 
     &                              NiterCG, 'CG_3d', deps, k)
#  endif
!
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                self%cg3d_x(i,j,k)=self%cg3d_x(i,j,k)+                  &
     &                             cg_a(iterCG,k)*                      &
     &                             self%cg3d_p(i,j,k)     ! step 9
                self%cg3d_r(i,j,k)=self%cg3d_r(i,j,k)+                  &
     &                             cg_a(iterCG,k)*                      &
     &                             self%cg3d_q(i,j,k)     ! step 10
              END DO
            END DO
!
!  Reorthogonalize r(k+1) with respect to r(0:k).
!
            dotrd=dotr
            dotr = dot_prod_2d (ng, tile, model, ctype,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          self%cg3d_r(:,:,k),                     &
     &                          self%cg3d_r(:,:,k))
!
            cg_b(iterCG+1,k)=dotr/dotrd                     ! step 12
!
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                self%cg3d_p(i,j,k)=-self%cg3d_r(i,j,k)+                 &
     &                             cg_b(iterCG+1,k)*                    &
     &                             self%cg3d_p(i,j,k)       ! step 13
                tl_A(i,j,k)=self%cg3d_p(i,j,k)
              END DO
            END DO
          END DO LEVEL_LOOP3

        END DO KLAP_ITER
!
!  Update control vector 3D variable with the implicit diffusion
!  solution iterate.
!
        LEVEL_LOOP4 : DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j,k)=self%cg3d_x(i,j,k)
            END DO
          END DO
!
!  Use the LAPACK routine DSTEQR to compute the Ritz vectors and
!  Ritz values of the tridiagonal matrix.
!
          DO i=1,NiterCG+1
            cg_Rv(i)=1.0_r8/cg_a(i-1,k)
            IF (i.gt.1) THEN
              cg_Rv(i)=cg_Rv(i)+cg_b(i-1,k)/cg_a(i-2,k)
            END IF
          END DO
          DO i=1,NiterCG
            zwork(i)=-SQRT(cg_b(i,k))/cg_a(i-1,k)
          END DO
!
          CALL dsteqr ('I', NiterCG+1, cg_Rv, zwork,                    &
     &                 zgv, NiterCG+1, work, info)
          IF (info.ne.0) THEN
            WRITE (stdout,*) ' multiscale_CG_2d_tl - Error in DSTEQR:', &
     &                       ' info = ', info
            exit_flag=8
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          eigMax(iterDiff,k)=MAXVAL(cg_Rv) ! maximum Ritz eigenvalue
          eigMin(iterDiff,k)=MINVAL(cg_Rv) ! minimum Ritz eigenvalue
!
!  Reset iterative conjugate gradient arrays.
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%cg3d_p(i,j,k)=0.0_r8
              self%cg3d_r(i,j,k)=0.0_r8
              self%cg3d_q(i,j,k)=0.0_r8
            END DO
          END DO
        END DO LEVEL_LOOP4

      END DO KDIFF_ITER
!
!  Convert the solution back to physical space.
!  (cf Eqn 18 of Weaver et al 2016, QJRMS, 143,455-471).
!  Not needed in the CI algorithm as the tl_scale factors cancel.
!
      DO k=1,N(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
            tl_A(i,j,k)=self%cg3d_x(i,j,k)/tl_scale(i,j)
          END DO
        END DO
      END DO
!
      RETURN
      END SUBROUTINE multiscale_CG_3d_tl
!
!-----------------------------------------------------------------------
!  This routine inverts the implicit diffusion linear operator, Ax = b,
!  using Chebyshev Iterations (CI) to model the multiscale propagation
!  of background-error covariance for 3D variables in the control
!  vector. It operates in conjunction with its adjoint to enforce
!  symmetry-preserving correlation functions. The algorithm requires
!  the extreme eigenvalues of matrix A, which are obtained from the
!  Conjugate Gradient (CG) solver. Additionally, initial values for
!  the estimate x(0) and its residual r(0) are required.
!
      SUBROUTINE multiscale_CI_3d_tl (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCI, Lweak,        &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CI iterations
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j, k
#  ifdef MULTI_SCALE_DEBUG
      integer                           :: status
#  endif
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: tl_scale
      real (r8), dimension(0:NiterCI,N(ng))         :: ci_alpha
      real (r8), dimension(0:NiterCI+1,N(ng))       :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_3d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u3dvar)
          Imin=IstrU
        CASE (v3dvar)
          Jmin=JstrV
      END SELECT
!
      self%ci3d_r=0.0_r8
      self%ci3d_p=0.0_r8
      self%ci3d_q=0.0_r8
      self%ci3d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms). Also, select pre-computed extreme
!  eigenvalues of matrix A.
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
         CASE ('u', 'u_eastward')              ! 3D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU3D(:,:,ms)
          eigMax => self%eigMaxU3D(:,:,ms)
        CASE ('v', 'v_northward')              ! 3D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV3D(:,:,ms)
          eigMax => self%eigMaxV3D(:,:,ms)
        CASE ('temp')                          ! temperature
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR3D(:,:,itemp,ms)
          eigMax => self%eigMaxR3D(:,:,itemp,ms)
        CASE ('salt')                          ! salinity
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR3D(:,:,isalt,ms)
          eigMax => self%eigMaxR3D(:,:,isalt,ms)
      END SELECT
!
!  Set control variable squared root area scale (2D).
!
      SELECT CASE (ctype)
        CASE (r3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
            END DO
          END DO
        CASE (u3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j))
            END DO
          END DO
        CASE (v3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j))
            END DO
          END DO
      END SELECT
!
!  Advance in K-space implicit pseudo-diffusion equation using 
!  Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        LEVEL_LOOP1 : DO k=1,N(ng)
          ci_sigma=0.5_r8*(eigMax(iterDiff,k)+                          &
     &                     eigMin(iterDiff,k))          ! step 1
          ci_delta=0.5_r8*(eigMax(iterDiff,k)-                          &
     &                     eigMin(iterDiff,k))          ! step 2
          ci_alpham1=1.0_r8                             ! step 3
          ci_beta(0,k)=0.0_r8                           ! step 4
!
!  Invert the implicit diffusion equation using Chebyshev Iterations.
!
!     A x = b  linear system
!
          DO iterCI=0,NiterCI                           ! steps 5 to 8
            IF (iterCI.eq.0) THEN
              ci_alpha(0,k)=1.0_r8/ci_sigma
              cff=ci_delta*ci_alpha(0,k)
              ci_beta(1,k)=0.5_r8*cff*cff
            ELSE
              ci_alpha(iterCI,k)=1.0_r8/(ci_sigma-                      &
     &                                   ci_beta(iterCI,k)/             &
     &                                   ci_alpha(iterCI-1,k))
              cff=0.5_r8*ci_delta*ci_alpha(iterCI,k)
              ci_beta(iterCI+1,k)=cff*cff
            END IF
          END DO
!                                            why r = -tl_A & p = +tl_A?
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j,k)=tl_scale(i,j)*tl_A(i,j,k)
              self%ci3d_x(i,j,k)=0.0_r8                 ! step 9
              self%ci3d_r(i,j,k)=-tl_A(i,j,k)           ! step 10
              self%ci3d_p(i,j,k)=tl_A(i,j,k)            ! step 11
            END DO
          END DO

#  ifdef MULTI_SCALE_DEBUG
!
          dotn(k) = dot_prod_2d (ng, tile, model, ctype,                &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           self%ci3d_r(:,:,k),                    &
     &                           self%ci3d_r(:,:,k))
#  endif
        END DO LEVEL_LOOP1
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=0,NiterCI                 ! steps 12 to 17
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
                tl_A(i,j,k)=tl_A(i,j,k)/tl_scale(i,j)
              END DO
            END DO
          END DO LEVEL_LOOP2
!
          SELECT CASE (ctype)                           ! step 12
            CASE (r3dvar)
              CALL tl_multiscale_Klap_r3d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
            CASE (u3dvar)
              CALL tl_multiscale_Klap_u3d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
            CASE (v3dvar)
              CALL tl_multiscale_Klap_v3d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     tl_A)
          END SELECT
!                                                        step 13
          LEVEL_LOOP3 : DO k=1,N(ng)
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                tl_A(i,j,k)=tl_A(i,j,k)*tl_scale(i,j)
                self%ci3d_q(i,j,k)=tl_A(i,j,k)
              END DO
            END DO

#  ifdef MULTI_SCALE_DEBUG
!
            dotr = dot_prod_2d (ng, tile, model, ctype,                 &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          self%ci3d_r(:,:,k),                     &
     &                          self%ci3d_r(:,:,k))
            deps=sqrt(dotr/dotn(k))
!
            status=multiscale_print(ng, ifield, ms, iterDiff, iterCI,   & 
     &                              NiterCI, 'CI_3d', deps, k)
#  endif
!                                                        steps 14, 15
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                self%ci3d_x(i,j,k)=self%ci3d_x(i,j,k)+                  &
     &                             ci_alpha(iterCI,k)*                  &
     &                             self%ci3d_p(i,j,k)
                self%ci3d_r(i,j,k)=self%ci3d_r(i,j,k)+                  &
     &                             ci_alpha(iterCI,k)*                  &
     &                             self%ci3d_q(i,j,k)
              END DO
            END DO
!                                                        step 16
            DO j=Jmin,Jmax
              DO i=Imin,Imax
                self%ci3d_p(i,j,k)=-self%ci3d_r(i,j,k)+                &
     &                             ci_beta(iterCI+1,k)*                &
     &                             self%ci3d_p(i,j,k)
                tl_A(i,j,k)=self%ci3d_p(i,j,k)
              END DO
            END DO
          END DO LEVEL_LOOP3

        END DO KLAP_ITER
!
!  Reset iterative conjugate gradient arrays.
!
        LEVEL_LOOP4 : DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              self%ci3d_p(i,j,k)=0.0_r8
              self%ci3d_r(i,j,k)=0.0_r8
              self%ci3d_q(i,j,k)=0.0_r8
            END DO
          END DO
!
!  Update control vector 3D variable with the implicit diffusion
!  solution iterate.
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              tl_A(i,j,k)=self%ci3d_x(i,j,k)/tl_scale(i,j)
            END DO
          END DO
        END DO LEVEL_LOOP4
!
      END DO KDIFF_ITER
!
      RETURN
      END SUBROUTINE multiscale_CI_3d_tl
!
!-----------------------------------------------------------------------
!  This routine serves as the adjoint to the inversion of the implicit
!  diffusion linear operator, Ax = b, utilizing Chebyshev Iterations
!  (CI) to model the multiscale propagation of background-error
!  covariance for 3D variables in the control vector. It is needed to
!  enforce symmetry-preserving correlation functions. The algorithm
!  requires the extreme eigenvalues of matrix A, which are obtained
!  from the Conjugate Gradient (CG) solver.
!
      SUBROUTINE multiscale_CI_3d_ad (self, ng, tile, model, ifield,    &
     &                                ctype, ms, NiterCI, Lweak,        &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CI iterations
      logical,            intent(in   ) :: Lweak     ! weak constraint
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:,:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j, k
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(IminS:ImaxS,JminS:JmaxS) :: ad_scale
      real (r8), dimension(0:NiterCI,N(ng))         :: ci_alpha
      real (r8), dimension(0:NiterCI+1,N(ng))       :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_3d_ad"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u3dvar)
          Imin=IstrU
        CASE (v3dvar)
          Jmin=JstrV
      END SELECT
!
      self%ci3d_r=0.0_r8
      self%ci3d_p=0.0_r8
      self%ci3d_q=0.0_r8
      self%ci3d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms). Also, select pre-computed extreme
!  eigenvalues of matrix A.
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
         CASE ('u', 'u_eastward')              ! 3D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU3D(:,:,ms)
          eigMax => self%eigMaxU3D(:,:,ms)
        CASE ('v', 'v_northward')              ! 3D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV3D(:,:,ms)
          eigMax => self%eigMaxV3D(:,:,ms)
        CASE ('temp')                          ! temperature
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR3D(:,:,itemp,ms)
          eigMax => self%eigMaxR3D(:,:,itemp,ms)
        CASE ('salt')                          ! salinity
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR3D(:,:,isalt,ms)
          eigMax => self%eigMaxR3D(:,:,isalt,ms)
      END SELECT
!
!  Set control variable squared root area scale (2D).
!
      SELECT CASE (ctype)
        CASE (r3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_r(i,j)*GRID(ng)%on_r(i,j))
            END DO
          END DO
        CASE (u3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_u(i,j)*GRID(ng)%on_u(i,j))
            END DO
          END DO
        CASE (v3dvar)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
              ad_scale(i,j)=SQRT(GRID(ng)%om_v(i,j)*GRID(ng)%on_v(i,j))
            END DO
          END DO
      END SELECT
!
!  Adjoint of advance in K-space implicit pseudo-diffusion equation
!  using Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=Mlap,1,-1
!
!  Initialize implicit diffusion operator for the Chebyshev Iterations
!  solver. Non-adjointable steps.
!
        LEVEL_LOOP4 : DO k=1,N(ng)
          ci_sigma=0.5_r8*(eigMax(iterDiff,k)+                          &
     &                     eigMin(iterDiff,k))          ! step 1
          ci_delta=0.5_r8*(eigMax(iterDiff,k)-                          &
     &                     eigMin(iterDiff,k))          ! step 2
          ci_alpham1=1.0_r8                             ! step 3
          ci_beta(0,k)=0.0_r8                           ! step 4
!
          DO iterCI=0,NiterCI                           ! steps 5 to 8
            IF (iterCI.eq.0) THEN
              ci_alpha(0,k)=1.0_r8/ci_sigma
              cff=ci_delta*ci_alpha(0,k)
              ci_beta(1,k)=0.5_r8*cff*cff
            ELSE
              ci_alpha(iterCI,k)=1.0_r8/(ci_sigma-                      &
     &                                   ci_beta(iterCI,k)/             &
     &                                   ci_alpha(iterCI-1,k))
              cff=0.5_r8*ci_delta*ci_alpha(iterCI,k)
              ci_beta(iterCI+1,k)=cff*cff
            END IF
          END DO
!
!  Adjoint of update control vector 3D variable with the implicit
!  diffusion solution iterate.
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            tl_A(i,j,k)=self%ci3d_x(i,j,k)/tl_scale(i,j)
!^
              self%ci3d_x(i,j,k)=ad_A(i,j,k)/ad_scale(i,j)
              ad_A(i,j,k)=0.0_r8
            END DO
          END DO
!
!  Reset iterative conjugate gradient arrays.
!
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            self%ci3d_q(i,j,k)=0.0_r8
!^
              self%ci3d_q(i,j,k)=0.0_r8
!^            self%ci3d_r(i,j,k)=0.0_r8
!^
              self%ci3d_r(i,j,k)=0.0_r8
!^            self%ci3d_p(i,j,k)=0.0_r8
!^
              self%ci3d_p(i,j,k)=0.0_r8
            END DO
          END DO
        END DO LEVEL_LOOP4
!
!  Adjoint of compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=NiterCI,0,-1              ! steps 17 to 12
!
          LEVEL_LOOP3 : DO k=1,N(ng)
!                                                        step 16
            DO j=Jmin,Jmax
              DO i=Imin,Imax
!^              tl_A(i,j,k)=self%ci3d_p(i,j,k)
!^
                self%ci3d_p(i,j,k)=self%ci3d_p(i,j,k)+ad_A(i,j,k)
                ad_A(i,j,k)=0.0_r8
!^              self%ci3d_p(i,j,k)=-self%ci3d_r(i,j,k)+                 &
!^   &                             ci_beta(iterCI+1,k)*                 &
!^   &                             self%ci3d_p(i,j,k)
!^
                self%ci3d_r(i,j,k)=self%ci3d_r(i,j,k)-                  &
     &                             self%ci3d_p(i,j,k)
                self%ci3d_p(i,j,k)=ci_beta(iterCI+1,k)*                 &
     &                             self%ci3d_p(i,j,k)
              END DO
            END DO
!                                                        steps 15, 14
            DO j=Jmin,Jmax
              DO i=Imin,Imax
!^              self%ci3d_r(i,j,k)=self%ci3d_r(i,j,k)+                  &
!^   &                             ci_alpha(iterCI,k)*                  &
!^   &                             self%ci3d_q(i,j,k)
!^
                self%ci3d_q(i,j,k)=self%ci3d_q(i,j,k)+                  &
     &                             ci_alpha(iterCI,k)*                  &
     &                             self%ci3d_r(i,j,k)
!^              self%ci3d_x(i,j,k)=self%ci3d_x(i,j,k)+                  &
!^   &                             ci_alpha(iterCI,k)*                  &
!^   &                             self%ci3d_p(i,j,k)
!^
                self%ci3d_p(i,j,k)=self%ci3d_p(i,j,k)+                  &
     &                             ci_alpha(iterCI,k)*                  &
     &                             self%ci3d_x(i,j,k)
              END DO
            END DO
!                                                        step 13
            DO j=Jmin,Jmax
              DO i=Imin,Imax
!^              self%ci3d_q(i,j,k)=tl_A(i,j,k)
!^
                ad_A(i,j,k)=ad_A(i,j,k)+self%ci3d_q(i,j,k)
                self%ci3d_q(i,j,k)=0.0_r8
!^              tl_A(i,j,k)=tl_A(i,j,k)*tl_scale(i,j)
!^
                ad_A(i,j,k)=ad_A(i,j,k)*ad_scale(i,j)
              END DO
            END DO

          END DO LEVEL_LOOP3
!
          SELECT CASE (ctype)                           ! step 12
            CASE (r3dvar)
              CALL ad_multiscale_Klap_r3d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     ad_A)
            CASE (u3dvar)
              CALL ad_multiscale_Klap_u3d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     ad_A)
            CASE (v3dvar)
              CALL ad_multiscale_Klap_v3d (self, ng, tile, model,       &
     &                                     ifield, ctype, ms, Lweak,    &
     &                                     LBi, UBi, LBj, UBj,          &
     &                                     IminS, ImaxS, JminS, JmaxS,  &
     &                                     ad_A)
          END SELECT
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            DO j=Jstr,Jend
              DO i=Istr,Iend
!^              tl_A(i,j,k)=tl_A(i,j,k)/tl_scale(i,j)
!^
                ad_A(i,j,k)=ad_A(i,j,k)/ad_scale(i,j)
              END DO
            END DO
          END DO LEVEL_LOOP2

        END DO KLAP_ITER
!
        LEVEL_LOOP1 : DO k=1,N(ng)
          DO j=Jmin,Jmax
            DO i=Imin,Imax
!^            self%ci3d_p(i,j,k)=tl_A(i,j,k)            ! step 11
!^
              ad_A(i,j,k)=ad_A(i,j,k)+self%ci3d_p(i,j,k)
              self%ci3d_p(i,j,k)=0.0_r8
!^            self%ci3d_r(i,j,k)=-tl_A(i,j,k)           ! step 10
!^
              ad_A(i,j,k)=ad_A(i,j,k)-self%ci3d_r(i,j,k)
              self%ci3d_r(i,j,k)=0.0_r8
!^            self%ci3d_x(i,j,k)=0.0_r8                 ! step 9
!^
              self%ci3d_x(i,j,k)=0.0_r8
!^            tl_A(i,j,k)=tl_scale(i,j)*tl_A(i,j,k)
!^
              ad_A(i,j,k)=ad_scale(i,j)*ad_A(i,j,k)
            END DO
          END DO
        END DO LEVEL_LOOP1
!
      END DO KDIFF_ITER
!
      RETURN
      END SUBROUTINE multiscale_CI_3d_ad
!
!-----------------------------------------------------------------------
!  It advances the vertical pseudo-difusion operator using an implicit
!  solver using parabolic splines or inverting a tridiagonal matrix.
!  Note that error correlations are considered separable in the
!  horizontal and vertical directions.
!
      SUBROUTINE multiscale_Vdiff_tl (self, ng, tile, model, ifield,    &
     &                                ctype, NVsteps,                   &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                DTsizeV, Kv, tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: tl_A(LBi:,LBj:,:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
!
      real (r8)                         :: cff
!
      real (r8), pointer                :: mask(:,:) => NULL()
!
      real(r8), dimension(LBi:UBi,LBj:UBj,N(ng),2)         :: tl_Awrk
#  if defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng))   :: oHz
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: FC
#  else
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#  endif
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: tl_DC
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
#  ifdef MASKING
        CASE (r3dvar)
          mask => GRID(ng)%rmask
#  endif
        CASE (u3dvar)
          Imin=IstrU
#  ifdef MASKING
          mask => GRID(ng)%umask
#  endif
        CASE (v3dvar)
          Jmin=JstrV
#  ifdef MASKING
          mask => GRID(ng)%vmask
#  endif
      END SELECT
!
      Nold=1
      Nnew=2
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
#  ifdef SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
          END DO
#  else
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
            FC(i,j,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-        &
     &                                    GRID(ng)%z_r(i,j,k  ))
          END DO
          FC(i,j,0)=0.0_r8
#  endif
        END DO
      END DO
!
!  Set initial conditions.
!
      DO k=1,N(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
!^          Awrk(i,j,k,Nold)=A(i,j,k)
!^
            tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
          END DO
        END DO
      END DO

#  ifdef SPLINES_VCONV
!
!  Integrate vertical diffusion equation implicitly using parabolic
!  splines.
!
      DO step=1,NVsteps
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
        DO j=Jmin,Jmax
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Imin,Imax
              FC(i,k)=cff1*GRID(ng)%Hz(i,j,k  )-                        &
     &                DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                        &
     &                DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Imin,Imax
            CF(i,0)=0.0_r8
!^          DC(i,0)=0.0_r8
!^
            tl_DC(i,0)=0.0_r8
          END DO
!
!  LU decomposition and forward substitution.
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Imin,Imax
              BC(i,k)=cff1*(GRID(ng)%Hz(i,j,k)+GRID(ng)%Hz(i,j,k+1))+   &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
!^            DC(i,k)=cff*(Awrk(i,j,k+1,Nold)-Awrk(i,j,k,Nold)-         &
!^   &                     FC(i,k)*DC(i,k-1))
!^
              tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
     &                        tl_Awrk(i,j,k  ,Nold)-                    &
     &                        FC(i,k)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Backward substitution.
!
          DO i=Imin,Imax
!^          DC(i,N(ng))=0.0_r8
!^
            tl_DC(i,N(ng))=0.0_r8
          END DO
          DO k=N(ng)-1,1,-1
            DO i=Imin,Imax
!^            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
            END DO
          END DO
!
          DO k=1,N(ng)
            DO i=Imin,Imax
!^            DC(i,k)=DC(i,k)*Kv(i,j,k)
!^
              tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nold)+                        &
!^   &                         DTsizeV*oHz(i,j,k)*                      &
!^   &                         (DC(i,k)-DC(i,k-1))
!^
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
     &                            DTsizeV*oHz(i,j,k)*                   &
     &                            (tl_DC(i,k)-tl_DC(i,k-1))
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO

#  else
!
!  Integrate vertical diffusion equation implicitly.
!
      DO step=1,NVsteps
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the tangent diffusion equation into tl_DC.
!
        DO j=Jmin,Jmax
          DO k=1,N(ng)
            DO i=Imin,Imax
              BC(i,k)=GRID(ng)%Hz(i,j,k)-FC(i,j,k)-FC(i,j,k-1)
!^            DC(i,k)=Awrk(i,j,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
              tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*GRID(ng)%Hz(i,j,k)
            END DO
          END DO
!
!  Solve the tangent linear tridiagonal system.
!
          DO i=Imin,Imax
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
!^          DC(i,1)=cff*DC(i,1)
!^
            tl_DC(i,1)=cff*tl_DC(i,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Imin,Imax
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
!^            DC(i,k)=cff*(DC(i,k)-FC(i,j,k-1)*DC(i,k-1))
!^
              tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Imin,Imax
!^          DC(i,N(ng))=(DC(i,N(ng))-                                   &
!^   &                   FC(i,j,N(ng)-1)*DC(i,N(ng)-1))/                &
!^   &                  (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^
            tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
     &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
     &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^          Awrk(i,j,N(ng),Nnew)=DC(i,N(ng))
!^
            tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
#   ifdef MASKING
!^          Awrk(i,j,N(ng),Nnew)=Awrk(i,j,N(ng),Nnew)*mask(i,j)
!^
            tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*mask(i,j)
#   endif
          END DO
!
          DO k=N(ng)-1,1,-1
            DO i=Imin,Imax
!^            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
              tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^            Awrk(i,j,k,Nnew)=DC(i,k)
!^
              tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
#   ifdef MASKING
!^            Awrk(i,j,k,Nnew)=Awrk(i,j,k,Nnew)*mask(i,j)
!^
              tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*mask(i,j)
#   endif
            END DO
          END DO
        END DO
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#  endif
!
!  Load convolved solution.
!
      DO k=1,N(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
!^          A(i,j,k)=Awrk(i,j,k,Nold)
!^
            tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
          END DO
        END DO
      END DO

#  ifdef DISTRIBUTE
!
!  Exchange boundary data.
!
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    A)
!^
      CALL mp_exchange3d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    tl_A)
#  endif
!
      RETURN
      END SUBROUTINE multiscale_Vdiff_tl
!
!-----------------------------------------------------------------------
!  It advances the vertical pseudo-difusion adjoint operator using an
!  implicit solver using parabolic splines or inverting a tridiagonal
!  matrix. Note that error correlations are considered separable in the
!  horizontal and vertical directions.
!
      SUBROUTINE multiscale_Vdiff_ad (self, ng, tile, model, ifield,    &
     &                                ctype, NVsteps,                   &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                DTsizeV, Kv, ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: ad_A(LBi:,LBj:,:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
!
      real (r8)                         :: adfac, cff
!
      real (r8), pointer                :: mask(:,:) => NULL()
!
      real(r8), dimension(LBi:UBi,LBj:UBj,N(ng),2)         :: ad_Awrk
#  if defined SPLINES_VCONV
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,N(ng))   :: oHz
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: FC
#  else
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS,0:N(ng)) :: FC
#  endif
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: BC
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: CF
      real(r8), dimension(IminS:ImaxS,0:N(ng))             :: ad_DC
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
#  ifdef MASKING
        CASE (r3dvar)
          mask => GRID(ng)%rmask
#  endif
        CASE (u3dvar)
          Imin=IstrU
#  ifdef MASKING
          mask => GRID(ng)%umask
#  endif
        CASE (v3dvar)
          Jmin=JstrV
#  ifdef MASKING
          mask => GRID(ng)%vmask
#  endif
      END SELECT
!
      Nold=1
      Nnew=2
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      DO j=Jmin,Jmax
        DO i=Imin,Imax
#  ifdef SPLINES_VCONV
          DO k=1,N(ng)
            oHz(i,j,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
          END DO
#  else
          FC(i,j,N(ng))=0.0_r8
          DO k=1,N(ng)-1
            FC(i,j,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-        &
     &                                    GRID(ng)%z_r(i,j,k  ))
          END DO
          FC(i,j,0)=0.0_r8
#  endif
        END DO
      END DO
!
!  Adjoint of load convolved solution.
!
#  ifdef DISTRIBUTE
!^    CALL mp_exchange3d (ng, tile, model, 1,                           &
!^   &                    LBi, UBi, LBj, UBj, 1, N(ng),                 &
!^   &                    NghostPoints,                                 &
!^   &                    EWperiodic(ng), NSperiodic(ng),               &
!^   &                    tl_A)
!^
      CALL ad_mp_exchange3d (ng, tile, model, 1,                        &
     &                       LBi, UBi, LBj, UBj, 1, N(ng),              &
     &                       NghostPoints,                              &
     &                       EWperiodic(ng), NSperiodic(ng),            &
     &                       ad_A)
!
#  endif

      DO k=1,N(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
!^          tl_A(i,j,k)=tl_Awrk(i,j,k,Nold)
!^
            ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+ad_A(i,j,k)
            ad_A(i,j,k)=0.0_r8
          END DO
        END DO
      END DO

#  ifdef SPLINES_VCONV
!
!  Integrate adjoint vertical diffusion equation implicitly using
!  parabolic splines.
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
        DO j=Jmin,Jmax
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
!  Compute basic state time-invariant coefficients.
!
          cff1=1.0_r8/6.0_r8
          DO k=1,N(ng)-1
            DO i=Imin,Imax
              FC(i,k)=cff1*GRID(ng)%Hz(i,j,k  )-                        &
     &                DTsizeV*Kv(i,j,k-1)*oHz(i,j,k  )
              CF(i,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                        &
     &                DTsizeV*Kv(i,j,k+1)*oHz(i,j,k+1)
            END DO
          END DO
          DO i=Imin,Imax
            CF(i,0)=0.0_r8
          END DO
!
          cff1=1.0_r8/3.0_r8
          DO k=1,N(ng)-1
            DO i=Imin,Imax
              BC(i,k)=cff1*(GRID(ng)%Hz(i,j,k)+GRID(ng)%Hz(i,j,k+1))+   &
     &                DTsizeV*Kv(i,j,k)*(oHz(i,j,k)+oHz(i,j,k+1))
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
              CF(i,k)=cff*CF(i,k)
            END DO
          END DO
!
!  Adjoint backward substitution.
!
          DO k=1,N(ng)
            DO i=Imin,Imax
!^            tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nold)+                  &
!^   &                            DTsizeV*oHz(i,j,k)*                   &
!^   &                            (tl_DC(i,k)-tl_DC(i,k-1))
!^
              adfac=DTsizeV*oHz(i,j,k)*ad_Awrk(i,j,k,Nnew)
              ad_DC(i,k-1)=ad_DC(i,k-1)-adfac
              ad_DC(i,k  )=ad_DC(i,k  )+adfac
              ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+                  &
     &                            ad_Awrk(i,j,k,Nnew)
              ad_Awrk(i,j,k,Nnew)=0.0_r8
!^            tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^
              ad_DC(i,k)=ad_DC(i,k)*Kv(i,j,k)
            END DO
          END DO
          DO k=1,N(ng)-1
            DO i=Imin,Imax
!^            tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
              ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
            END DO
          END DO
          DO i=Imin,Imax
!^          tl_DC(i,N(ng))=0.0_r8
!^
            ad_DC(i,N(ng))=0.0_r8
          END DO
!
!  Adjoint LU decomposition and forward substitution.
!
          DO k=N(ng)-1,1,-1
            DO i=Imin,Imax
              cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!^            tl_DC(i,k)=cff*(tl_Awrk(i,j,k+1,Nold)-                    &
!^   &                        tl_Awrk(i,j,k  ,Nold)-                    &
!^   &                        FC(i,k)*tl_DC(i,k-1))
!^
              adfac=cff*ad_DC(i,k)
              ad_Awrk(i,j,k  ,Nold)=ad_Awrk(i,j,k  ,Nold)-adfac
              ad_Awrk(i,j,k+1,Nold)=ad_Awrk(i,j,k+1,Nold)+adfac
              ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
!
          DO i=Imin,Imax
!^          tl_DC(i,0)=0.0_r8
!^
            ad_DC(i,0)=0.0_r8
          END DO
        END DO
      END DO

#  else
!
!  Integrate adjoint vertical diffusion equation implicitly.
!
      DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
        DO j=Jmin,Jmax
!
!  Compute diagonal matrix coefficients BC.
!
          DO k=1,N(ng)
            DO i=Imin,Imax
              BC(i,k)=GRID(ng)%Hz(i,j,k)-FC(i,j,k)-FC(i,j,k-1)
            END DO
          END DO
!
!  Compute new solution by back substitution.
!
          DO i=Imin,Imax
            cff=1.0_r8/BC(i,1)
            CF(i,1)=cff*FC(i,j,1)
          END DO
          DO k=2,N(ng)-1
            DO i=Imin,Imax
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
              CF(i,k)=cff*FC(i,j,k)
            END DO
          END DO
!^        DO k=N(ng)-1,1,-1
!^
          DO k=1,N(ng)-1
            DO i=Imin,Imax
#   ifdef MASKING
!^            tl_Awrk(i,j,k,Nnew)=tl_Awrk(i,j,k,Nnew)*mask(i,j)
!^
              ad_Awrk(i,j,k,Nnew)=ad_Awrk(i,j,k,Nnew)*mask(i,j)
#   endif
!^            tl_Awrk(i,j,k,Nnew)=tl_DC(i,k)
!^
              ad_DC(i,k)=ad_DC(i,k)+                                    &
     &                   ad_Awrk(i,j,k,Nnew)
              ad_Awrk(i,j,k,Nnew)=0.0_r8
!^            tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
              ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
            END DO
          END DO
!
          DO i=Imin,Imax
#   ifdef MASKING
!^          tl_Awrk(i,j,N(ng),Nnew)=tl_Awrk(i,j,N(ng),Nnew)*mask(i,j)
!^
            ad_Awrk(i,j,N(ng),Nnew)=ad_Awrk(i,j,N(ng),Nnew)*mask(i,j)
#   endif
!^          tl_Awrk(i,j,N(ng),Nnew)=tl_DC(i,N(ng))
!^
            ad_DC(i,N(ng))=ad_DC(i,N(ng))+                              &
     &                     ad_Awrk(i,j,N(ng),Nnew)
            ad_Awrk(i,j,N(ng),Nnew)=0.0_r8
!^          tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                             &
!^   &                      FC(i,j,N(ng)-1)*tl_DC(i,N(ng)-1))/          &
!^   &                     (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
!^
            adfac=ad_DC(i,N(ng))/                                       &
     &            (BC(i,N(ng))-FC(i,j,N(ng)-1)*CF(i,N(ng)-1))
            ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-FC(i,j,N(ng)-1)*adfac
            ad_DC(i,N(ng)  )=adfac
          END DO
!
!  Solve the adjoint tridiagonal system.
!
!^        DO k=2,N(ng)-1
!^
          DO k=N(ng)-1,2,-1
            DO i=Imin,Imax
              cff=1.0_r8/(BC(i,k)-FC(i,j,k-1)*CF(i,k-1))
!^            tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,j,k-1)*tl_DC(i,k-1))
!^
              adfac=cff*ad_DC(i,k)
              ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,j,k-1)*adfac
              ad_DC(i,k  )=adfac
            END DO
          END DO
          DO i=Imin,Imax
            cff=1.0_r8/BC(i,1)
!^          tl_DC(i,1)=cff*tl_DC(i,1)
!^
            ad_DC(i,1)=cff*ad_DC(i,1)
          END DO
!
!  Adjoint of right-hand-side terms for the diffusion equation.
!
          DO k=1,N(ng)
            DO i=Imin,Imax
!^            tl_DC(i,k)=tl_Awrk(i,j,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
              ad_Awrk(i,j,k,Nold)=ad_Awrk(i,j,k,Nold)+                  &
     &                            GRID(ng)%Hz(i,j,k)*ad_DC(i,k)
              ad_DC(i,k)=0.0_r8
            END DO
          END DO
        END DO
      END DO
#  endif
!
!  Set adjoint initial conditions.
!
      DO k=1,N(ng)
        DO j=Jmin,Jmax
          DO i=Imin,Imax
!^          tl_Awrk(i,j,k,Nold)=tl_A(i,j,k)
!^
            ad_A(i,j,k)=ad_A(i,j,k)+ad_Awrk(i,j,k,Nold)
            ad_Awrk(i,j,k,Nold)=0.0_r8
          END DO
        END DO
      END DO
!
      RETURN
      END SUBROUTINE multiscale_Vdiff_ad
# endif /* SOLVE3D */

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  This routine computes the extreme eigenvalues required by the
!  Chebyshev Iterations (CI) solver, which is applied to implicit
!  diffusion operators in the modeling of the multiscale background-
!  error covariance for the boundary adjustments of 2D variables in
!  the control vector.
!
!  The minimum and maximum eigenvalues of the K-Laplacian operator
!  are determined using a Conjugate Gradient (CG) solution to the
!  corresponding linear system, A x = b.
!
      SUBROUTINE multiscale_CG_b1d_tl (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCG,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCG   ! CG iterations
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBij:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCG, iterDiff, info
      integer                           :: i, j, k
#  ifdef MULTI_SCALE_DEBUG
      integer                           :: status
#  endif
!
      real (r8)                         :: dotn, dotpq, dotr, dotrd
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(LBij:UBij)           :: tl_scale
      real (r8), dimension(0:NiterCG)           :: cg_a
      real (r8), dimension(0:NiterCG+1)         :: cg_b
      real (r8), dimension(NiterCG+1)           :: cg_Rv
      real (r8), dimension(2*(NiterCG+1)-2)     :: work
      real (r8), dimension(NiterCG)             :: zwork
      real (r8), dimension(NiterCG+1,NiterCG+1) :: zgv
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CG_b1d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      self%cgB1d_r=0.0_r8
      self%cgB1d_p=0.0_r8
      self%cgB1d_q=0.0_r8
      self%cgB1d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')                          ! free surface
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU1D(:,ibry,ms)
          eigMax => self%eigMaxU1D(:,ibry,ms)
        CASE ('ubar','ubar_eastward')          ! 2D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU1D(:,ibry,ms)
          eigMax => self%eigMaxU1D(:,ibry,ms)
        CASE ('vbar', 'vbar_northward')        ! 2D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV1D(:,ibry,ms)
          eigMax => self%eigMaxV1D(:,ibry,ms)
      END SELECT
!
!  Set control variable squared root area scale.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_r(i,j))
#  ifdef MASKING
                 tl_A(j)=GRID(ng)%rmask(i,j)*tl_A(j)
#  endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_r(i,j))
#  ifdef MASKING
                 tl_A(i)=GRID(ng)%rmask(i,j)*tl_A(i)
#  endif
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_u(i,j))
#  ifdef MASKING
                 tl_A(j)=GRID(ng)%umask(i,j)*tl_A(j)
#  endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_u(i,j))
#  ifdef MASKING
                 tl_A(i)=GRID(ng)%umask(i,j)*tl_A(i)
#  endif
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_v(i,j))
#  ifdef MASKING
                 tl_A(j)=GRID(ng)%vmask(i,j)*tl_A(j)
#  endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_v(i,j))
#  ifdef MASKING
                 tl_A(i)=GRID(ng)%vmask(i,j)*tl_A(i)
#  endif
               END DO
             END IF
        END SELECT
      END IF
!
!  Advance in K-space implicit pseudo-diffusion 2D equation using 
!  Weaver et al. (2016) algorithm 1.
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmin
              self%cgB1d_x(j)=0.0_r8                      ! step 1
              self%cgB1d_r(j)=-tl_A(j)                    ! step 2
              self%cgB1d_p(j)=tl_A(j)                     ! step 3
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imin
              self%cgB1d_x(i)=0.0_r8                      ! step 1
              self%cgB1d_r(i)=-tl_A(i)                    ! step 2
              self%cgB1d_p(i)=tl_A(i)                     ! step 3
            END DO
        END SELECT
        cg_b(0)=0.0_r8                                    ! step 4
!
        dotn = dot_prod_1d (ng, tile, model, ctype, ibry,               &
     &                      LBij, UBij,                                 &
     &                      self%cgB1d_r,                               &
     &                      self%cgB1d_r)
!
!  Invert the implicit diffusion equation using CG iterations.
!
!     A x = b  linear system
!
        KLAP_ITER : DO iterCG=0,NiterCG                   ! steps 6 - 16
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                tl_A(j)=tl_A(j)/tl_scale(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                tl_A(i)=tl_A(i)/tl_scale(i)
              END DO
          END SELECT
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
          CALL tl_multiscale_Klap_b1d (self, ng, tile, model,           &
     &                                 ifield, ibry, ctype, ms,         &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 tl_A)
!
!  Conjugate Gradient (CG) algorithm, (Weaver et al., 2016).
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                tl_A(j)=tl_A(j)*tl_scale(j)
                self%cgB1d_q(j)=tl_A(j)                   ! step 7
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                tl_A(i)=tl_A(i)*tl_scale(i)
                self%cgB1d_q(i)=tl_A(i)                   ! step 7
              END DO
          END SELECT
!
          dotr = dot_prod_1d (ng, tile, model, ctype, ibry,             &
     &                        LBij, UBij,                               &
     &                        self%cgB1d_r,                             &
     &                        self%cgB1d_r)
          dotpq = dot_prod_1d (ng, tile, model, ctype, ibry,            &
     &                         LBij, UBij,                              &
     &                         self%cgB1d_p,                            &
     &                         self%cgB1d_q)
!
          cg_a(iterCG)=dotr/dotpq                         ! step 8

#  ifdef MULTI_SCALE_DEBUG
!
          deps=SQRT(dotr/dotn)
          status=multiscale_print(ng, ifield, ms, iterDiff, iterCG,     & 
     &                            NiterCG, 'CG_b1d', deps)
#  endif
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                self%cgB1d_x(j)=self%cgB1d_x(j)+                        &
     &                          cg_a(iterCG)*                           &
     &                          self%cgB1d_p(j)           ! step 9
                self%cgB1d_r(j)=self%cgB1d_r(j)+                        &
     &                          cg_a(iterCG)*                           &
     &                          self%cgB1d_q(j)           ! step 10
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                self%cgB1d_x(i)=self%cgB1d_x(i)+                        &
     &                          cg_a(iterCG)*                           &
     &                          self%cgB1d_p(i)           ! step 9
                self%cgB1d_r(i)=self%cgB1d_r(i)+                        &
     &                          cg_a(iterCG)*                           &
     &                          self%cgB1d_q(i)           ! step 10
              END DO
          END SELECT
!
!  Reorthogonalize r(k+1) with respect to r(0:k).           step 11
!
          dotrd=dotr
          dotr = dot_prod_1d (ng, tile, model, ctype, ibry,             &
     &                        LBij, UBij,                               &
     &                        self%cgB1d_r,                             &
     &                        self%cgB1d_r)
!
          cg_b(iterCG+1)=dotr/dotrd                       ! step 12
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                self%cgB1d_p(j)=-self%cgB1d_r(j)+                       &
     &                          cg_b(iterCG+1)*                         &
     &                          self%cgB1d_p(j)           ! step 13
                tl_A(j)=self%cgB1d_p(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                self%cgB1d_p(i)=-self%cgB1d_r(i)+                       &
     &                          cg_b(iterCG+1)*                         &
     &                          self%cgB1d_p(i)           ! step 13
                tl_A(i)=self%cgB1d_p(i)
              END DO
          END SELECT
!
        END DO KLAP_ITER
!
!  Update control vector 2D variable lateral boundary consition with the
!  implicit diffusion solution iterate.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmin
              tl_A(j)=self%cgB1d_x(j)
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imin
              tl_A(i)=self%cgB1d_x(i)
            END DO
        END SELECT
!
!  Use the LAPACK routine DSTEQR to compute the Ritz vectors and 
!  Ritz values of the tridiagonal matrix.
!
        DO i=1,NiterCG+1                                  ! step 14
          cg_Rv(i)=1.0_r8/cg_a(i-1)
          IF (i.gt.1) THEN
            cg_Rv(i)=cg_Rv(i)+cg_b(i-1)/cg_a(i-2)
          END IF
        END DO
        DO i=1,NiterCG
          zwork(i)=-SQRT(cg_b(i))/cg_a(i-1)
        END DO
!
        CALL dsteqr ('I', NiterCG+1, cg_Rv, zwork,                      &
     &               zgv, NiterCG+1, work, info)
        IF (info.ne.0) THEN
          WRITE (stdout,*) ' multiscale_CG_b1d_tl - Error in DSTEQR:',  &
     &                     ' info = ', info
          exit_flag=8
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        END IF
!
        eigMax(iterDiff)=MAXVAL(cg_Rv)  ! maximum Ritz eigenvalue
        eigMin(iterDiff)=MINVAL(cg_Rv)  ! minimum Ritz eigenvalue
!
!  Reset iterative conjugate gradient arrays.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmin
              self%cgB1d_p(j)=0.0_r8
              self%cgB1d_r(j)=0.0_r8
              self%cgB1d_q(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imin
              self%cgB1d_p(i)=0.0_r8
              self%cgB1d_r(i)=0.0_r8
              self%cgB1d_q(i)=0.0_r8
            END DO
        END SELECT
!
      END DO KDIFF_ITER
!
!  Convert the solution back to physical space.
!  (cf Eqn 18 of Weaver et al 2016, QJRMS, 143,455-471).
!  Not needed in the CI algorithm as the tl_scale factors cancel.
!
      SELECT CASE (ibry)
        CASE (iwest, ieast)
          DO j=Jmin,Jmin
            tl_A(j)=self%cgB1d_x(j)/tl_scale(j)
          END DO
        CASE (isouth, inorth)
          DO i=Imin,Imin
            tl_A(i)=self%cgB1d_x(i)/tl_scale(i)
          END DO
      END SELECT
!
      RETURN
      END SUBROUTINE multiscale_CG_b1d_tl
!
!-----------------------------------------------------------------------
!  This routine inverts the implicit diffusion linear operator, Ax = b,
!  using Chebyshev Iterations (CI) to model the multiscale propagation
!  of background-error covariance for the boundary adjustment of 2D
!  variables in the control vector. It operates in conjunction with
!  its adjoint to enforce symmetry-preserving correlation functions.
!  The algorithm requires the extreme eigenvalues of matrix A, which
!  are obtained from the Conjugate Gradient (CG) solver. Additionally,
!  initial values for the estimate x(0) and its residual r(0) are
!  required.
!
      SUBROUTINE multiscale_CI_b1d_tl (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCI,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CG iterations
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBij:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j
#  ifdef MULTI_SCALE_DEBUG
      integer                           :: status
#  endif
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(LBij:UBij)   :: tl_scale
      real (r8), dimension(0:NiterCI)   :: ci_alpha
      real (r8), dimension(0:NiterCI+1) :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_b1d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      self%ciB1d_r=0.0_r8
      self%ciB1d_p=0.0_r8
      self%ciB1d_q=0.0_r8
      self%ciB1d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')                          ! free surface
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR1D(:,ibry,ms)
          eigMax => self%eigMaxR1D(:,ibry,ms)
        CASE ('ubar', 'ubar_eastward')         ! 2D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU1D(:,ibry,ms)
          eigMax => self%eigMaxU1D(:,ibry,ms)
        CASE ('vbar', 'vbar_northward')        ! 2D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV1D(:,ibry,ms)
          eigMax => self%eigMaxV1D(:,ibry,ms)
      END SELECT
!
!  Set control variable squared root area scale.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Advance in K-space implicit pseudo-diffusion equation using 
!  Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        ci_sigma=0.5_r8*(eigMax(iterDiff)+                              &
     &                   eigMin(iterDiff))              ! step 1
        ci_delta=0.5_r8*(eigMax(iterDiff)-                              &
     &                   eigMin(iterDiff))              ! step 2
        ci_alpham1=1.0_r8                               ! step 3
        ci_beta(0)=0.0_r8                               ! step 4
!
!  Invert the implicit diffusion equation using Chebyshev Iterations.
!
!     A x = b  linear system
!
        DO iterCI=0,NiterCI                             ! steps 5 to 8
          IF (iterCI.eq.0) THEN
            ci_alpha(0)=1.0_r8/ci_sigma
            cff=ci_delta*ci_alpha(0)
            ci_beta(1)=0.5_r8*cff*cff
          ELSE
            ci_alpha(iterCI)=1.0_r8/(ci_sigma-                          &
     &                               ci_beta(iterCI)/                   &
     &                               ci_alpha(iterCI-1))
            cff=0.5_r8*ci_delta*ci_alpha(iterCI)
            ci_beta(iterCI+1)=cff*cff
          END IF
        END DO
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmin
              tl_A(j)=tl_scale(j)*tl_A(j)
              self%ciB1d_x(j)=0.0_r8                    ! step 9
              self%ciB1d_r(j)=-tl_A(j)                  ! step 10
              self%ciB1d_p(j)=tl_A(j)                   ! step 11
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imin
              tl_A(i)=tl_scale(i)*tl_A(i)
              self%ciB1d_x(i)=0.0_r8                    ! step 9
              self%ciB1d_r(i)=-tl_A(i)                  ! step 10
              self%ciB1d_p(i)=tl_A(i)                   ! step 11
            END DO
        END SELECT

#  ifdef MULTI_SCALE_DEBUG
!
        dotn = dot_prod_1d (ng, tile, model, ctype, ibry,               &
     &                      LBij, UBij,                                 &
     &                      self%ciB1d_r,                               &
     &                      self%ciB1d_r)
#  endif
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=0,NiterCI                 ! steps 12 to 17
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                tl_A(j)=tl_A(j)/tl_scale(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                tl_A(i)=tl_A(i)/tl_scale(i)
              END DO
          END SELECT
!
          CALL tl_multiscale_Klap_b1d (self, ng, tile, model,           &
     &                                 ifield, ibry, ctype, ms,         &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 tl_A)            ! step 12
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                tl_A(j)=tl_A(j)*tl_scale(j)
                self%ciB1d_q(j)=tl_A(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                tl_A(i)=tl_A(i)*tl_scale(i)
                self%ciB1d_q(i)=tl_A(i)
              END DO
          END SELECT

#  ifdef MULTI_SCALE_DEBUG
!
          dotr = dot_prod_1d (ng, tile, model, ctype, ibry,             &
     &                        LBij, UBij,                               &
     &                        self%ciB1d_r,                             &
     &                        self%ciB1d_r)
          deps=SQRT(dotr/dotn)
!
          status=multiscale_print(ng, ifield, ms, iterDiff, iterCI,     & 
     &                            NiterCI, 'CI_B1d', deps)
#  endif
!                                                         steps 14, 15
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                self%ciB1d_x(j)=self%ciB1d_x(j)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_p(j)
                self%ciB1d_r(j)=self%ciB1d_r(j)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_q(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                self%ciB1d_x(i)=self%ciB1d_x(i)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_p(i)
                self%ciB1d_r(i)=self%ciB1d_r(i)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_q(i)
              END DO
          END SELECT
!                                                         step 16
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                self%ciB1d_p(j)=-self%ciB1d_r(j)+                       &
     &                          ci_beta(iterCI+1)*self%ciB1d_p(j)
                tl_A(j)=self%ciB1d_p(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                self%ciB1d_p(i)=-self%ciB1d_r(i)+                       &
     &                          ci_beta(iterCI+1)*self%ciB1d_p(i)
                tl_A(i)=self%ciB1d_p(i)
              END DO
          END SELECT
!
        END DO KLAP_ITER
!
!  Reset iterative conjugate gradient arrays.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmin
              self%ciB1d_p(j)=0.0_r8
              self%ciB1d_r(j)=0.0_r8
              self%ciB1d_q(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imin
              self%ciB1d_p(i)=0.0_r8
              self%ciB1d_r(i)=0.0_r8
              self%ciB1d_q(i)=0.0_r8
            END DO
        END SELECT
!
!  Update control vector lateral boundary variable with the implicit
!  diffusion solution iterate.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmin
              tl_A(j)=self%ciB1d_x(j)/tl_scale(j)
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imin
              tl_A(i)=self%ciB1d_x(i)/tl_scale(i)
            END DO
        END SELECT
!
      END DO KDIFF_ITER
!
      RETURN
      END SUBROUTINE multiscale_CI_b1d_tl
!
!-----------------------------------------------------------------------
!  This routine serves as the adjoint to the inversion of the implicit
!  diffusion linear operator, Ax = b, using Chebyshev Iterations (CI)
!  to model the multiscale propagation of background-error covariance
!  for the boundary adjustment of 2D variables in the control vector.
!  The algorithm requires the extreme eigenvalues of matrix A, which
!  are obtained from the Conjugate Gradient (CG) solver. Additionally,
!  initial values for the estimate x(0) and its residual r(0) are
!  required.
!
      SUBROUTINE multiscale_CI_b1d_ad (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCI,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CG iterations
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBij:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: adfac, cff
      real (r8)                         :: ci_alpham1, ci_delta, ci_sigma
!
      real (r8), pointer                :: eigMin(:) => NULL()
      real (r8), pointer                :: eigMax(:) => NULL()
!
      real (r8), dimension(LBij:UBij)   :: ad_scale
      real (r8), dimension(0:NiterCI)   :: ci_alpha
      real (r8), dimension(0:NiterCI+1) :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_b1d_ad"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      self%ciB1d_r=0.0_r8
      self%ciB1d_p=0.0_r8
      self%ciB1d_q=0.0_r8
      self%ciB1d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('zeta')                          ! free surface
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR1D(:,ibry,ms)
          eigMax => self%eigMaxR1D(:,ibry,ms)
        CASE ('ubar', 'ubar_eastward')         ! 2D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU1D(:,ibry,ms)
          eigMax => self%eigMaxU1D(:,ibry,ms)
        CASE ('vbar', 'vbar_northward')        ! 2D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV1D(:,ibry,ms)
          eigMax => self%eigMaxV1D(:,ibry,ms)
      END SELECT
!
!  Set control variable squared root area scale.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Adjoint of advance in K-space implicit pseudo-diffusion equation
!  using Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=Mlap,1,-1
!
        ci_sigma=0.5_r8*(eigMax(iterDiff)+                              &
     &                   eigMin(iterDiff))              ! step 1
        ci_delta=0.5_r8*(eigMax(iterDiff)-                              &
     &                   eigMin(iterDiff))              ! step 2
        ci_alpham1=1.0_r8                               ! step 3
        ci_beta(0)=0.0_r8                               ! step 4
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmin
              self%ciB1d_p(j)=0.0_r8
              self%ciB1d_r(j)=0.0_r8
              self%ciB1d_q(j)=0.0_r8
              self%ciB1d_x(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imin
              self%ciB1d_p(i)=0.0_r8
              self%ciB1d_r(i)=0.0_r8
              self%ciB1d_q(i)=0.0_r8
              self%ciB1d_x(i)=0.0_r8
            END DO
        END SELECT
!
!  Adjoint of invert the implicit diffusion equation using Chebyshev
!  Iterations.
!
!     A x = b  linear system
!
        DO iterCI=0,NiterCI                             ! steps 5 to 8
          IF (iterCI.eq.0) THEN
            ci_alpha(0)=1.0_r8/ci_sigma
            cff=ci_delta*ci_alpha(0)
            ci_beta(1)=0.5_r8*cff*cff
          ELSE
            ci_alpha(iterCI)=1.0_r8/(ci_sigma-                          &
     &                               ci_beta(iterCI)/                   &
     &                               ci_alpha(iterCI-1))
            cff=0.5_r8*ci_delta*ci_alpha(iterCI)
            ci_beta(iterCI+1)=cff*cff
          END IF
        END DO
!
!  Adjoint of update control vector lateral boundary variable with the
!  implicit diffusion solution iterate.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmin
!^            tl_A(j)=self%ciB1d_x(j)/tl_scale(j)
!^
              self%ciB1d_x(j)=ad_A(j)/ad_scale(j)
              ad_A(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imin
!^            tl_A(i)=self%ciB1d_x(i)/tl_scale(i)
!^
              self%ciB1d_x(i)=ad_A(i)/ad_scale(i)
              ad_A(i)=0.0_r8
            END DO
        END SELECT
!
!  Adjoint of reset iterative conjugate gradient arrays.
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmin
!^            self%ciB1d_q(j)=0.0_r8
!^
              self%ciB1d_q(j)=0.0_r8
!^            self%ciB1d_r(j)=0.0_r8
!^
              self%ciB1d_r(j)=0.0_r8
!^            self%ciB1d_p(j)=0.0_r8
!^
              self%ciB1d_p(j)=0.0_r8
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imin
!^            self%ciB1d_q(i)=0.0_r8
!^
              self%ciB1d_q(i)=0.0_r8
!^            self%ciB1d_r(i)=0.0_r8
!^
              self%ciB1d_r(i)=0.0_r8
!^            self%ciB1d_p(i)=0.0_r8
!^
              self%ciB1d_p(i)=0.0_r8
            END DO
        END SELECT
!
!  Adjoint of compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=NiterCI,0,-1              ! steps 17 to 12
!
          SELECT CASE (ibry)                            ! step 16
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
!^              tl_A(j)=self%ciB1d_p(j)
!^
                self%ciB1d_p(j)=self%ciB1d_p(j)+ad_A(j)
                ad_A(j)=0.0_r8
!^              self%ciB1d_p(j)=-self%ciB1d_r(j)+                       &
!^   &                          ci_beta(iterCI+1)*self%ciB1d_p(j)
!^
                self%ciB1d_r(j)=self%ciB1d_r(j)-self%ciB1d_p(j)
                self%ciB1d_p(j)=ci_beta(iterCI+1)*self%ciB1d_p(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
!^              tl_A(i)=self%ciB1d_p(i)
!^
                self%ciB1d_p(i)=self%ciB1d_p(i)+ad_A(i)
                ad_A(i)=0.0_r8
!^              self%ciB1d_p(i)=-self%ciB1d_r(i)+                       &
!^   &                          ci_beta(iterCI+1)*self%ciB1d_p(i)
!^
                self%ciB1d_r(i)=self%ciB1d_r(i)-self%ciB1d_p(i)
                self%ciB1d_p(i)=ci_beta(iterCI+1)*self%ciB1d_p(i)
              END DO
          END SELECT
!                                                         steps 15, 14
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
!^              self%ciB1d_r(j)=self%ciB1d_r(j)+                        &
!^   &                          ci_alpha(iterCI)*self%ciB1d_q(j)
!^
                self%ciB1d_q(j)=self%ciB1d_q(j)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_r(j)
!^              self%ciB1d_x(j)=self%ciB1d_x(j)+                        &
!^   &                          ci_alpha(iterCI)*self%ciB1d_p(j)
!^
                self%ciB1d_p(j)=self%ciB1d_p(j)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_x(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
!^              self%ciB1d_r(i)=self%ciB1d_r(i)+                        &
!^   &                          ci_alpha(iterCI)*self%ciB1d_q(i)
!^
                self%ciB1d_q(i)=self%ciB1d_q(i)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_r(i)
!^              self%ciB1d_x(i)=self%ciB1d_x(i)+                        &
!^   &                          ci_alpha(iterCI)*self%ciB1d_p(i)
!^
                self%ciB1d_p(i)=self%ciB1d_p(i)+                        &
     &                          ci_alpha(iterCI)*self%ciB1d_x(i)
              END DO
          END SELECT
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
!^              self%ciB1d_q(j)=tl_A(j)
!^
                ad_A(j)=ad_A(j)+self%ciB1d_q(j)
                self%ciB1d_q(j)=0.0_r8
!^              tl_A(j)=tl_A(j)*tl_scale(j)
!^
                ad_A(j)=ad_A(j)*ad_scale(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
!^              self%ciB1d_q(i)=tl_A(i)
!^
                ad_A(i)=ad_A(i)+self%ciB1d_q(i)
                self%ciB1d_q(i)=0.0_r8
!^              tl_A(i)=tl_A(i)*tl_scale(i)
!^
                ad_A(i)=ad_A(i)*ad_scale(i)
              END DO
          END SELECT
!
          CALL ad_multiscale_Klap_b1d (self, ng, tile, model,           &
     &                                 ifield, ibry, ctype, ms,         &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 ad_A)            ! step 12
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
!^              tl_A(j)=tl_A(j)/tl_scale(j)
!^
                ad_A(j)=ad_A(j)/ad_scale(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
!^              tl_A(i)=tl_A(i)/tl_scale(i)
!^
                ad_A(i)=ad_A(i)/ad_scale(i)
              END DO
          END SELECT
!
        END DO KLAP_ITER
!
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmin
!^            self%ciB1d_p(j)=tl_A(j)                   ! step 11
!^
              ad_A(j)=ad_A(j)+self%ciB1d_p(j)
              self%ciB1d_p(j)=0.0_r8
!^            self%ciB1d_r(j)=-tl_A(j)                  ! step 10
!^
              ad_A(j)=ad_A(j)-self%ciB1d_r(j)
!^            self%ciB1d_x(j)=0.0_r8                    ! step 9
!^
              self%ciB1d_x(j)=0.0_r8
!^            tl_A(j)=tl_scale(j)*tl_A(j)
!^
              ad_A(j)=ad_scale(j)*ad_A(j)
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imin
!^            self%ciB1d_p(i)=tl_A(i)                   ! step 11
!^
              ad_A(i)=ad_A(i)+self%ciB1d_p(i)
              self%ciB1d_p(i)=0.0_r8
!^            self%ciB1d_r(i)=-tl_A(i)                  ! step 10
!^
              ad_A(i)=ad_A(i)-self%ciB1d_r(i)
!^            self%ciB1d_x(i)=0.0_r8                    ! step 9
!^
              self%ciB1d_x(i)=0.0_r8
!^            tl_A(i)=tl_scale(i)*tl_A(i)
!^
              ad_A(i)=ad_scale(i)*ad_A(i)
            END DO
        END SELECT
!
      END DO KDIFF_ITER
!
      RETURN
      END SUBROUTINE multiscale_CI_b1d_ad

#  ifdef SOLVE3D
!
!-----------------------------------------------------------------------
!  This routine computes the extreme eigenvalues required by the
!  Chebyshev Iterations (CI) solver, which is applied to implicit
!  diffusion operators in the modeling of the multiscale background-
!  error covariance for the boundary adjustments of 3D variables in
!  the control vector.
!
!  The minimum and maximum eigenvalues of the K-Laplacian operator
!  are determined using a Conjugate Gradient (CG) solution to the
!  corresponding linear system, A x = b.
!
      SUBROUTINE multiscale_CG_b2d_tl (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCG,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCG   ! CG iterations
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCG, iterDiff, info
      integer                           :: i, j, k
#   ifdef MULTI_SCALE_DEBUG
      integer                           :: status
#   endif
!
      real (r8)                         :: dotpq, dotr, dotrd
      real (r8), dimension(N(ng))       :: dotn
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(LBij:UBij)           :: tl_scale
      real (r8), dimension(0:NiterCG,N(ng))     :: cg_a
      real (r8), dimension(0:NiterCG+1,N(ng))   :: cg_b
      real (r8), dimension(NiterCG+1)           :: cg_Rv
      real (r8), dimension(2*(NiterCG+1)-2)     :: work
      real (r8), dimension(NiterCG)             :: zwork
      real (r8), dimension(NiterCG+1,NiterCG+1) :: zgv
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CG_b2d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      self%cgB2d_r=0.0_r8
      self%cgB2d_p=0.0_r8
      self%cgB2d_q=0.0_r8
      self%cgB2d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward')               ! 3D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU1DZ(:,:,ibry,ms)
          eigMax => self%eigMaxU1DZ(:,:,ibry,ms)
        CASE ('v', 'v_northward')              ! 3D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV1DZ(:,:,ibry,ms)
          eigMax => self%eigMaxV1DZ(:,:,ibry,ms)
        CASE ('temp')                          ! temperature
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR1DZ(:,:,itemp,ibry,ms)
          eigMax => self%eigMaxR1DZ(:,:,itemp,ibry,ms)
        CASE ('salt')                          ! salinity
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR1DZ(:,:,isalt,ibry,ms)
          eigMax => self%eigMaxR1DZ(:,:,isalt,ibry,ms)
      END SELECT
!
!  Set control variable squared root area scale.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_r(i,j))
#   ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(j,k)=GRID(ng)%rmask(i,j)*tl_A(j,k)
                 END DO
#   endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_r(i,j))
#   ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(i,k)=GRID(ng)%rmask(i,j)*tl_A(i,k)
                 END DO
#   endif
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_u(i,j))
#   ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(j,k)=GRID(ng)%umask(i,j)*tl_A(j,k)
                 END DO
#   endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_u(i,j))
#   ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(i,k)=GRID(ng)%umask(i,j)*tl_A(i,k)
                 END DO
#   endif
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_v(i,j))
#   ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(j,k)=GRID(ng)%vmask(i,j)*tl_A(j,k)
                 END DO
#   endif
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_v(i,j))
#   ifdef MASKING
                 DO k=1,N(ng)
                   tl_A(i,k)=GRID(ng)%vmask(i,j)*tl_A(i,k)
                 END DO
#   endif
               END DO
             END IF
        END SELECT
      END IF
!
!  Advance in K-space implicit pseudo-diffusion 2D equation using 
!  Weaver et al. (2016) algorithm 1.
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        LEVEL_LOOP1 : DO k=1,N(ng)
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                self%cgB2d_x(j,k)=0.0_r8                  ! step 1
                self%cgB2d_r(j,k)=-tl_A(j,k)              ! step 2
                self%cgB2d_p(j,k)=tl_A(j,k)               ! step 3
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                self%cgB2d_x(i,k)=0.0_r8                  ! step 1
                self%cgB2d_r(i,k)=-tl_A(i,k)              ! step 2
                self%cgB2d_p(i,k)=tl_A(i,k)               ! step 3
              END DO
          END SELECT
          cg_b(0,k)=0.0_r8                                ! step 4
!
          dotn(k) = dot_prod_1d (ng, tile, model, ctype, ibry,          &
     &                           LBij, UBij,                            &
     &                           self%cgB2d_r(:,k),                     &
     &                           self%cgB2d_r(:,k))
        END DO LEVEL_LOOP1
!
!  Invert the implicit diffusion equation using CG iterations.
!
!     A x = b  linear system
!
        KLAP_ITER : DO iterCG=0,NiterCG                   ! steps 6 - 16
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
                  tl_A(j,k)=tl_A(j,k)/tl_scale(j)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
                  tl_A(i,k)=tl_A(i,k)/tl_scale(i)
                END DO
            END SELECT
          END DO LEVEL_LOOP2
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
          CALL tl_multiscale_Klap_b2d (self, ng, tile, model,           &
     &                                 ifield, ibry, ctype, ms,         &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 tl_A)
!
!  Conjugate Gradient (CG) algorithm, (Weaver et al., 2016).
!
          LEVEL_LOOP3 : DO k=1,N(ng)
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
                  tl_A(j,k)=tl_A(j,k)*tl_scale(j)
                  self%cgB2d_q(j,k)=tl_A(j,k)             ! step 7
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
                  tl_A(i,k)=tl_A(i,k)*tl_scale(i)
                  self%cgB2d_q(i,k)=tl_A(i,k)             ! step 7
                END DO
            END SELECT
!
            dotr = dot_prod_1d (ng, tile, model, ctype, ibry,           &
     &                          LBij, UBij,                             &
     &                          self%cgB2d_r(:,k),                      &
     &                          self%cgB2d_r(:,k))
            dotpq = dot_prod_1d (ng, tile, model, ctype, ibry,          &
     &                           LBij, UBij,                            &
     &                           self%cgB2d_p(:,k),                     &
     &                           self%cgB2d_q(:,k))
!
            cg_a(iterCG,k)=dotr/dotpq                     ! step 8

#   ifdef MULTI_SCALE_DEBUG
!
            deps=SQRT(dotr/dotn)
            status=multiscale_print(ng, ifield, ms, iterDiff, iterCG,   & 
     &                              NiterCG, 'CG_b2d', deps,k)
#   endif
!
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
                  self%cgB2d_x(j,k)=self%cgB2d_x(j,k)+                  &
     &                              cg_a(iterCG,k)*                     &
     &                              self%cgB2d_p(j,k)     ! step 9
                  self%cgB2d_r(j,k)=self%cgB2d_r(j,k)+                  &
     &                              cg_a(iterCG,k)*                     &
     &                              self%cgB2d_q(j,k)     ! step 10
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
                  self%cgB2d_x(i,k)=self%cgB2d_x(i,k)+                  &
     &                              cg_a(iterCG,k)*                     &
     &                              self%cgB2d_p(i,k)     ! step 9
                  self%cgB2d_r(i,k)=self%cgB2d_r(i,k)+                  &
     &                              cg_a(iterCG,k)*                     &
     &                              self%cgB2d_q(i,k)     ! step 10
                END DO
            END SELECT
!
!  Reorthogonalize r(k+1) with respect to r(0:k).           step 11
!
            dotrd=dotr
            dotr = dot_prod_1d (ng, tile, model, ctype, ibry,           &
     &                          LBij, UBij,                             &
     &                          self%cgB2d_r(:,k),                      &
     &                          self%cgB2d_r(:,k))
!
            cg_b(iterCG+1,k)=dotr/dotrd                   ! step 12
!
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
                  self%cgB2d_p(j,k)=-self%cgB2d_r(j,k)+                 &
     &                              cg_b(iterCG+1,k)*                   &
     &                              self%cgB2d_p(j,k)     ! step 13
                  tl_A(j,k)=self%cgB2d_p(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
                  self%cgB2d_p(i,k)=-self%cgB2d_r(i,k)+                 &
     &                              cg_b(iterCG+1,k)*                   &
     &                              self%cgB2d_p(i,k)     ! step 13
                  tl_A(i,k)=self%cgB2d_p(i,k)
                END DO
            END SELECT
          END DO LEVEL_LOOP3
!
        END DO KLAP_ITER
!
!  Update control vector 2D variable with the implicit diffusion
!  solution iterate.
!
        LEVEL_LOOP4 : DO k=1,N(ng)
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                tl_A(j,k)=self%cgB2d_x(j,k)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                tl_A(i,k)=self%cgB2d_x(i,k)
              END DO
          END SELECT
!
!  Use the LAPACK routine DSTEQR to compute the Ritz vectors and 
!  Ritz values of the tridiagonal matrix.
!
          DO i=1,NiterCG+1                                ! step 14
            cg_Rv(i)=1.0_r8/cg_a(i-1,k)
            IF (i.gt.1) THEN
              cg_Rv(i)=cg_Rv(i)+cg_b(i-1,k)/cg_a(i-2,k)
            END IF
          END DO
          DO i=1,NiterCG
            zwork(i)=-SQRT(cg_b(i,k))/cg_a(i-1,k)
          END DO
!
          CALL dsteqr ('I', NiterCG+1, cg_Rv, zwork,                    &
     &                 zgv, NiterCG+1, work, info)
          IF (info.ne.0) THEN
            WRITE (stdout,*) ' multiscale_CG_b2d_tl - Error in DSTEQR:',&
     &                       ' info = ', info
            exit_flag=8
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
!
          eigMax(iterDiff,k)=MAXVAL(cg_Rv)  ! maximum Ritz eigenvalue
          eigMin(iterDiff,k)=MINVAL(cg_Rv)  ! minimum Ritz eigenvalue
!
!  Reset iterative conjugate gradient arrays.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                self%cgB2d_p(j,k)=0.0_r8
                self%cgB2d_r(j,k)=0.0_r8
                self%cgB2d_q(j,k)=0.0_r8
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                self%cgB2d_p(i,k)=0.0_r8
                self%cgB2d_r(i,k)=0.0_r8
                self%cgB2d_q(i,k)=0.0_r8
              END DO
          END SELECT
        END DO LEVEL_LOOP4
!
      END DO KDIFF_ITER
!
!  Convert the solution back to physical space.
!  (cf Eqn 18 of Weaver et al 2016, QJRMS, 143,455-471).
!  Not needed in the CI algorithm as the tl_scale factors cancel.
!
      SELECT CASE (ibry)
        CASE (iwest, ieast)
          DO k=1,N(ng)
            DO j=Jmin,Jmin
              tl_A(j,k)=self%cgB2d_x(j,k)/tl_scale(j)
            END DO
          END DO
        CASE (isouth, inorth)
          DO k=1,N(ng)
            DO i=Imin,Imin
              tl_A(i,k)=self%cgB2d_x(i,k)/tl_scale(i)
            END DO
          END DO
      END SELECT
!
      RETURN
      END SUBROUTINE multiscale_CG_b2d_tl
!
!-----------------------------------------------------------------------
!  This routine inverts the implicit diffusion linear operator, Ax = b,
!  using Chebyshev Iterations (CI) to model the multiscale propagation
!  of background-error covariance for the boundary adjustment of 3D
!  variables in the control vector. It operates in conjunction with
!  its adjoint to enforce symmetry-preserving correlation functions.
!  The algorithm requires the extreme eigenvalues of matrix A, which
!  are obtained from the Conjugate Gradient (CG) solver. Additionally,
!  initial values for the estimate x(0) and its residual r(0) are
!  required.
!
      SUBROUTINE multiscale_CI_b2d_tl (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCI,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CI iterations
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: tl_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j, k
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: cff, ci_alpham1, ci_delta
      real (r8)                         :: ci_sigma
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(LBij:UBij)         :: tl_scale
      real (r8), dimension(0:NiterCI,N(ng))   :: ci_alpha
      real (r8), dimension(0:NiterCI+1,N(ng)) :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_b2d_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      self%ciB2d_r=0.0_r8
      self%ciB2d_p=0.0_r8
      self%ciB2d_q=0.0_r8
      self%ciB2d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward')               ! 3D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU1DZ(:,:,ibry,ms)
          eigMax => self%eigMaxU1DZ(:,:,ibry,ms)
        CASE ('v', 'v_northward')              ! 3D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV1DZ(:,:,ibry,ms)
          eigMax => self%eigMaxV1DZ(:,:,ibry,ms)
        CASE ('temp')                          ! temperature
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR1DZ(:,:,itemp,ibry,ms)
          eigMax => self%eigMaxR1DZ(:,:,itemp,ibry,ms)
        CASE ('salt')                          ! salinity
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR1DZ(:,:,isalt,ibry,ms)
          eigMax => self%eigMaxR1DZ(:,:,isalt,ibry,ms)
      END SELECT
!
!  Set control variable squared root area scale.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 tl_scale(j)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 tl_scale(i)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Advance in K-space implicit pseudo-diffusion equation using 
!  Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=1,Mlap
!
        LEVEL_LOOP1 :  DO k=1,N(ng)
          ci_sigma=0.5_r8*(eigMax(iterDiff,k)+                          &
     &                     eigMin(iterDiff,k))          ! step 1
          ci_delta=0.5_r8*(eigMax(iterDiff,k)-                          &
     &                     eigMin(iterDiff,k))          ! step 2
          ci_alpham1=1.0_r8                             ! step 3
          ci_beta(0,k)=0.0_r8                           ! step 4
!
!  Invert the implicit diffusion equation using Chebyshev Iterations.
!
!     A x = b  linear system
!
          DO iterCI=0,NiterCI                           ! steps 5 to 8
            IF (iterCI.eq.0) THEN
              ci_alpha(0,k)=1.0_r8/ci_sigma
              cff=ci_delta*ci_alpha(0,k)
              ci_beta(1,k)=0.5_r8*cff*cff
            ELSE
              ci_alpha(iterCI,k)=1.0_r8/(ci_sigma-                      &
     &                                   ci_beta(iterCI,k)/             &
     &                                   ci_alpha(iterCI-1,k))
              cff=0.5_r8*ci_delta*ci_alpha(iterCI,k)
              ci_beta(iterCI+1,k)=cff*cff
            END IF
          END DO
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                tl_A(j,k)=tl_scale(j)*tl_A(j,k)
                self%ciB2d_x(j,k)=0.0_r8               ! step 9
                self%ciB2d_r(j,k)=-tl_A(j,k)           ! step 10
                self%ciB2d_p(j,k)=tl_A(j,k)            ! step 11
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                tl_A(i,k)=tl_scale(i)*tl_A(i,k)
                self%ciB2d_x(i,k)=0.0_r8               ! step 9
                self%ciB2d_r(i,k)=-tl_A(i,k)           ! step 10
                self%ciB2d_p(i,k)=tl_A(i,k)            ! step 11
              END DO
          END SELECT

#   ifdef MULTI_SCALE_DEBUG
!
          dotn(k) = dot_prod_1d (ng, tile, model, ctype, ibry,          &
     &                           LBij, UBij,                            &
     &                           self%ciB2d_r(:,k),                     &
     &                           self%ciB2d_r(:,k))
#   endif
        END DO LEVEL_LOOP1
!
!  Compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=0,NiterCI                 ! steps 12 to 17
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
                  tl_A(j,k)=tl_A(j,k)/tl_scale(j)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
                  tl_A(i,k)=tl_A(i,k)/tl_scale(i)
                END DO
            END SELECT
          END DO LEVEL_LOOP2
!
          CALL tl_multiscale_Klap_b2d (self, ng, tile, model,           &
     &                                 ifield, ibry, ctype, ms,         &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 tl_A)            ! step 12
!
          LEVEL_LOOP3 : DO k=1,N(ng)
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
                  tl_A(j,k)=tl_A(j,k)*tl_scale(j)
                  self%ciB2d_q(j,k)=tl_A(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
                  tl_A(i,k)=tl_A(i,k)*tl_scale(i)
                  self%ciB2d_q(i,k)=tl_A(i,k)
                END DO
            END SELECT

#   ifdef MULTI_SCALE_DEBUG
!
            dotr = dot_prod_1d (ng, tile, model, ctype, ibry,           &
     &                          LBij, UBij,                             &
     &                          self%ciB2d_r(:,k),                      &
     &                          self%ciB2d_r(:,k))
            deps=SQRT(dotr/dotn)
!
            status=multiscale_print(ng, ifield, ms, iterDiff, iterCI,   & 
     &                              NiterCI, 'CI_B2d', deps, k)
#   endif
!                                                         steps 14, 15
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
                  self%ciB2d_x(j,k)=self%ciB2d_x(j,k)+                    &
     &                              ci_alpha(iterCI,k)*                   &
     &                              self%ciB2d_p(j,k)
                  self%ciB2d_r(j,k)=self%ciB2d_r(j,k)+                    &
     &                              ci_alpha(iterCI,k)*                   &
     &                              self%ciB2d_q(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
                  self%ciB2d_x(i,k)=self%ciB2d_x(i,k)+                    &
     &                              ci_alpha(iterCI,k)*                   &
     &                              self%ciB2d_p(i,k)
                  self%ciB2d_r(i,k)=self%ciB2d_r(i,k)+                    &
     &                              ci_alpha(iterCI,k)*                   &
     &                              self%ciB2d_q(i,k)
                END DO
            END SELECT
!                                                         step 16
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
                  self%ciB2d_p(j,k)=-self%ciB2d_r(j,k)+                  &
     &                              ci_beta(iterCI+1,k)*                 &
     &                              self%ciB2d_p(j,k)
                  tl_A(j,k)=self%ciB2d_p(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
                  self%ciB2d_p(i,k)=-self%ciB2d_r(i,k)+                  &
     &                              ci_beta(iterCI+1,k)*                 &
     &                              self%ciB2d_p(i,k)
                    tl_A(i,k)=self%ciB2d_p(i,k)
                END DO
            END SELECT
          END DO LEVEL_LOOP3
!
        END DO KLAP_ITER
!
!  Reset iterative conjugate gradient arrays.
!
        LEVEL_LOOP4 : DO k=1,N(ng)
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                self%ciB2d_p(j,k)=0.0_r8
                self%ciB2d_r(j,k)=0.0_r8
                self%ciB2d_q(j,k)=0.0_r8
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                self%ciB2d_p(i,k)=0.0_r8
                self%ciB2d_r(i,k)=0.0_r8
                self%ciB2d_q(i,k)=0.0_r8
              END DO
          END SELECT
!
!  Update control vector lateral boundary variable with the implicit
!  diffusion solution iterate.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                tl_A(j,k)=self%ciB2d_x(j,k)/tl_scale(j)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                tl_A(i,k)=self%ciB2d_x(i,k)/tl_scale(i)
              END DO
          END SELECT
        END DO LEVEL_LOOP4
!
      END DO KDIFF_ITER
!
      RETURN
      END SUBROUTINE multiscale_CI_b2d_tl
!
!-----------------------------------------------------------------------
!  This routine serves as the adjoint to the inversion of the implicit
!  diffusion linear operator, Ax = b, using Chebyshev Iterations (CI)
!  to model the multiscale propagation of background-error covariance
!  for the boundary adjustment of 3D variables in the control vector.
!  The algorithm requires the extreme eigenvalues of matrix A, which
!  are obtained from the Conjugate Gradient (CG) solver. Additionally,
!  initial values for the estimate x(0) and its residual r(0) are
!  required.
!
      SUBROUTINE multiscale_CI_b2d_ad (self, ng, tile, model, ifield,   &
     &                                 ibry, ctype, ms, NiterCI,        &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: NiterCI   ! CG iterations
      integer,            intent(in   ) :: LBij, UBij
      integer,            intent(in   ) :: IminS, ImaxS, JminS, JmaxS
      real (r8),          intent(inout) :: ad_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Mlap, iterCI, iterDiff
      integer                           :: i, j, k
!
      real (r8)                         :: dotn, dotr, deps
      real (r8)                         :: adfac, cff
      real (r8)                         :: ci_alpham1, ci_delta, ci_sigma
!
      real (r8), pointer                :: eigMin(:,:) => NULL()
      real (r8), pointer                :: eigMax(:,:) => NULL()
!
      real (r8), dimension(LBij:UBij)         :: ad_scale
      real (r8), dimension(0:NiterCI,N(ng))   :: ci_alpha
      real (r8), dimension(0:NiterCI+1,N(ng)) :: ci_beta
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_CI_b2d_ad"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar)
          Imin=IstrU
        CASE (v2dvar)
          Jmin=JstrV
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      self%ciB2d_r=0.0_r8
      self%ciB2d_p=0.0_r8
      self%ciB2d_q=0.0_r8
      self%ciB2d_x=0.0_r8
!
!  Select number of K-Laplacian inverse operator applications (Mlap)
!  for requested variable in the 2D state/control vector and 
!  multiscale index/counter (ms).
!
!  Mlap MUST be greater than 2 and EVEN. Choose Mlap > or O(10) to
!  approximate a Gaussian correlation functions.
!
      SELECT CASE (TRIM(StateVarName(ifield)))
        CASE ('u', 'u_eastward')               ! 3D u-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinU1DZ(:,:,ibry,ms)
          eigMax => self%eigMaxU1DZ(:,:,ibry,ms)
        CASE ('v', 'v_northward')              ! 3D v-momentum
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinV1DZ(:,:,ibry,ms)
          eigMax => self%eigMaxV1DZ(:,:,ibry,ms)
        CASE ('temp')                          ! temperature
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR1DZ(:,:,itemp,ibry,ms)
          eigMax => self%eigMaxR1DZ(:,:,itemp,ibry,ms)
        CASE ('salt')                          ! salinity
          Mlap=self%Mlap(ifield,ms)
          eigMin => self%eigMinR1DZ(:,:,isalt,ibry,ms)
          eigMax => self%eigMaxR1DZ(:,:,isalt,ibry,ms)
      END SELECT
!
!  Set control variable squared root area scale.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ctype)
           CASE (r2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,r2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,r2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_r(i,j))
               END DO
             END IF
           CASE (u2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,u2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,u2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_u(i,j))
               END DO
             END IF
           CASE (v2dvar)
             IF ((ibry.eq.iwest).or.(ibry.eq.ieast)) THEN
               i=BOUNDS(ng)%edge(ibry,v2dvar)
               DO j=Jmin,Jmax
                 ad_scale(j)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             ELSE IF ((ibry.eq.isouth).or.(ibry.eq.inorth)) THEN
               j=BOUNDS(ng)%edge(ibry,v2dvar)
               DO i=Imin,Imax
                 ad_scale(i)=SQRT(GRID(ng)%on_v(i,j))
               END DO
             END IF
        END SELECT
      END IF
!
!  Adjoint of advance in K-space implicit pseudo-diffusion equation
!  using Weaver et al. (2016) algorithm 3 (fixed number of iterations
!  to a pre-determined value K).
!
      KDIFF_ITER : DO iterDiff=Mlap,1,-1
!
        LEVEL_LOOP4 :  DO k=1,N(ng)
          ci_sigma=0.5_r8*(eigMax(iterDiff,k)+                          &
     &                     eigMin(iterDiff,k))          ! step 1
          ci_delta=0.5_r8*(eigMax(iterDiff,k)-                          &
     &                     eigMin(iterDiff,k))          ! step 2
          ci_alpham1=1.0_r8                             ! step 3
          ci_beta(0,k)=0.0_r8                           ! step 4
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
                self%ciB2d_p(j,k)=0.0_r8
                self%ciB2d_r(j,k)=0.0_r8
                self%ciB2d_q(j,k)=0.0_r8
                self%ciB2d_x(j,k)=0.0_r8
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
                self%ciB2d_p(i,k)=0.0_r8
                self%ciB2d_r(i,k)=0.0_r8
                self%ciB2d_q(i,k)=0.0_r8
                self%ciB2d_x(i,k)=0.0_r8
              END DO
          END SELECT
!
!  Adjoint of invert the implicit diffusion equation using Chebyshev
!  Iterations.
!
!     A x = b  linear system
!
          DO iterCI=0,NiterCI                           ! steps 5 to 8
            IF (iterCI.eq.0) THEN
              ci_alpha(0,k)=1.0_r8/ci_sigma
              cff=ci_delta*ci_alpha(0,k)
              ci_beta(1,k)=0.5_r8*cff*cff
            ELSE
              ci_alpha(iterCI,k)=1.0_r8/(ci_sigma-                      &
     &                                   ci_beta(iterCI,k)/             &
     &                                   ci_alpha(iterCI-1,k))
              cff=0.5_r8*ci_delta*ci_alpha(iterCI,k)
              ci_beta(iterCI+1,k)=cff*cff
            END IF
          END DO
!
!  Adjoint of update control vector lateral boundary variable with the
!  implicit diffusion solution iterate.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
!^              tl_A(j,k)=self%ciB2d_x(j,k)/tl_scale(j)
!^
                self%ciB2d_x(j,k)=ad_A(j,k)/ad_scale(j)
                ad_A(j,k)=0.0_r8
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
!^              tl_A(i,k)=self%ciB2d_x(i,k)/tl_scale(i)
!^
                self%ciB2d_x(i,k)=ad_A(i,k)/ad_scale(i)
                ad_A(i,k)=0.0_r8
              END DO
          END SELECT
!
!  Adjoint of reset iterative conjugate gradient arrays.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
!^              self%ciB2d_q(j,k)=0.0_r8
!^
                self%ciB2d_q(j,k)=0.0_r8
!^              self%ciB2d_r(j,k)=0.0_r8
!^
                self%ciB2d_r(j,k)=0.0_r8
!^              self%ciB2d_p(j,k)=0.0_r8
!^
                self%ciB2d_p(j,k)=0.0_r8
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
!^              self%ciB2d_q(i,k)=0.0_r8
!^
                self%ciB2d_q(i,k)=0.0_r8
!^              self%ciB2d_r(i,k)=0.0_r8
!^
                self%ciB2d_r(i,k)=0.0_r8
!^              self%ciB2d_p(i,k)=0.0_r8
!^
                self%ciB2d_p(i,k)=0.0_r8
              END DO
          END SELECT
        END DO LEVEL_LOOP4
!
!  Adjoint of compute K-Laplacian, [1 + Del(K*Del)], operator.
!
        KLAP_ITER : DO iterCI=NiterCI,0,-1              ! steps 17 to 12
!
          LEVEL_LOOP3 : DO k=1,N(ng)
!                                                         step 16
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
!^                tl_A(j,k)=self%ciB2d_p(j,k)
!^
                  self%ciB2d_p(j,k)=self%ciB2d_p(j,k)+ad_A(j,k)
                  ad_A(j,k)=0.0_r8
!^                self%ciB2d_p(j,k)=-self%ciB2d_r(j,k)+                 &
!^   &                              ci_beta(iterCI+1,k)*                &
!&   &                              self%ciB2d_p(j,k)
!^
                  self%ciB2d_r(j,k)=self%ciB2d_r(j,k)-                  &
     &                              self%ciB2d_p(j,k)
                  self%ciB2d_p(j,k)=ci_beta(iterCI+1,k)*                &
     &                              self%ciB2d_p(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
!^                tl_A(i,k)=self%ciB2d_p(i,k)
!^
                  self%ciB2d_p(i,k)=self%ciB2d_p(i,k)+ad_A(i,k)
                  ad_A(i,k)=0.0_r8
!^                self%ciB2d_p(i,k)=-self%ciB2d_r(i,k)+                 &
!^   &                              ci_beta(iterCI+1,k)*                &
!^   &                              self%ciB2d_p(i,k)
!^
                  self%ciB2d_r(i,k)=self%ciB2d_r(i,k)-                  &
     &                              self%ciB2d_p(i,k)
                  self%ciB2d_p(i,k)=ci_beta(iterCI+1,k)*                &
     &                              self%ciB2d_p(i,k)
                END DO
            END SELECT
!                                                         steps 15, 14
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
!^                self%ciB2d_r(j,k)=self%ciB2d_r(j,k)+                  &
!^   &                              ci_alpha(iterCI,k)*                 &
!^   &                              self%ciB2d_q(j,k)
!^
                  self%ciB2d_q(j,k)=self%ciB2d_q(j,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_r(j,k)
!^                self%ciB2d_x(j,k)=self%ciB2d_x(j,k)+                  &
!^   &                              ci_alpha(iterCI,k)*                 &
!^   &                              self%ciB2d_p(j,k)
!^
                  self%ciB2d_p(j,k)=self%ciB2d_p(j,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_x(j,k)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
!^                self%ciB2d_r(i,k)=self%ciB2d_r(i,k)+                  &
!^   &                              ci_alpha(iterCI,k)*                 &
!^   &                              self%ciB2d_q(i,k)
!^
                  self%ciB2d_q(i,k)=self%ciB2d_q(i,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_r(i,k)
!^                self%ciB2d_x(i,k)=self%ciB2d_x(i,k)+                  &
!^   &                              ci_alpha(iterCI,k)*                 &
!^   &                              self%ciB2d_p(i,k)
!^
                  self%ciB2d_p(i,k)=self%ciB2d_p(i,k)+                  &
     &                              ci_alpha(iterCI,k)*                 &
     &                              self%ciB2d_x(i,k)
                END DO
            END SELECT
!
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
!^                self%ciB2d_q(j,k)=tl_A(j,k)
!^
                  ad_A(j,k)=ad_A(j,k)+self%ciB2d_q(j,k)
                  self%ciB2d_q(j,k)=0.0_r8
!^                tl_A(j,k)=tl_A(j,k)*tl_scale(j)
!^
                  ad_A(j,k)=ad_A(j,k)*ad_scale(j)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
!^                self%ciB2d_q(i,k)=tl_A(i,k)
!^
                  ad_A(i,k)=ad_A(i,k)+self%ciB2d_q(i,k)
                  self%ciB2d_q(i,k)=0.0_r8
!^                tl_A(i,k)=tl_A(i,k)*tl_scale(i)
!^
                  ad_A(i,k)=ad_A(i,k)*ad_scale(i)
                END DO
            END SELECT
!
          END DO LEVEL_LOOP3
!
          CALL ad_multiscale_Klap_b2d (self, ng, tile, model,           &
     &                                 ifield, ibry, ctype, ms,         &
     &                                 LBij, UBij,                      &
     &                                 IminS, ImaxS, JminS, JmaxS,      &
     &                                 ad_A)            ! step 12
!
          LEVEL_LOOP2 : DO k=1,N(ng)
            SELECT CASE (ibry)
              CASE (iwest, ieast)
                DO j=Jmin,Jmin
!^                tl_A(j,k)=tl_A(j,k)/tl_scale(j)
!^
                  ad_A(j,k)=ad_A(j,k)/ad_scale(j)
                END DO
              CASE (isouth, inorth)
                DO i=Imin,Imin
!^                tl_A(i,k)=tl_A(i,k)/tl_scale(i)
!^
                  ad_A(i,k)=ad_A(i,k)/ad_scale(i)
                END DO
            END SELECT
          END DO LEVEL_LOOP2
!
        END DO KLAP_ITER
!
        LEVEL_LOOP1 : DO k=1,N(ng)
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              DO j=Jmin,Jmin
!^              self%ciB2d_p(j,k)=tl_A(j,k)             ! step 11
!^
                ad_A(j,k)=ad_A(j,k)+self%ciB2d_p(j,k)
                self%ciB2d_p(j,k)=0.0_r8
!^              self%ciB2d_r(j,k)=-tl_A(j,k)            ! step 10
!^
                ad_A(j,k)=ad_A(j,k)-self%ciB2d_r(j,k)
!^              self%ciB2d_x(j,k)=0.0_r8                ! step 9
!^
                self%ciB2d_x(j,k)=0.0_r8
!^              tl_A(j,k)=tl_scale(j)*tl_A(j,k)
!^
                ad_A(j,k)=ad_scale(j)*ad_A(j,k)
              END DO
            CASE (isouth, inorth)
              DO i=Imin,Imin
!^              self%ciB2d_p(i,k)=tl_A(i,k)             ! step 11
!^
                ad_A(i,k)=ad_A(i,k)+self%ciB2d_p(i,k)
                self%ciB2d_p(i,k)=0.0_r8
!^              self%ciB2d_r(i,k)=-tl_A(i,k)            ! step 10
!^
                ad_A(i,k)=ad_A(i,k)-self%ciB2d_r(i,k)
!^              self%ciB2d_x(i,k)=0.0_r8                ! step 9
!^
                self%ciB2d_x(i,k)=0.0_r8
!^              tl_A(i,k)=tl_scale(i)*tl_A(i,k)
!^
                ad_A(i,k)=ad_scale(i)*ad_A(i,k)
              END DO
          END SELECT
        END DO LEVEL_LOOP1
!
      END DO KDIFF_ITER
!
      RETURN
      END SUBROUTINE multiscale_CI_b2d_ad
!
!-----------------------------------------------------------------------
!  It advances the vertical pseudo-difusion operator using an implicit
!  solver using parabolic splines or inverting a tridiagonal matrix for
!  adjusting lateral boundaries in the control vector.
!
      SUBROUTINE multiscale_bry_Vdiff_tl (self, ng, tile, model,        &
     &                                    ifield, ibry, ctype, NVsteps, &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    LBij, UBij,                   &
     &                                    DTsizeV, Kv, tl_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: LBij, UBij
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: tl_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
!
      real (r8)                         :: cff, cff1
!
      real (r8), pointer                :: mask(:,:) => NULL()
!
      real(r8), dimension(LBij:UBij,N(ng),2) :: tl_Awrk
#   if defined SPLINES_VCONV
      real(r8), dimension(LBij:UBij,N(ng))   :: oHz
#   endif
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: BC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: CF
      real(r8), dimension(LBij:UBij,0:N(ng)) :: tl_DC
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_bry_Vdiff_tl"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
#   ifdef MASKING
        CASE (r3dvar)
          mask => GRID(ng)%rmask
#   endif
        CASE (u3dvar)
          Imin=IstrU
#   ifdef MASKING
          mask => GRID(ng)%umask
#   endif
        CASE (v3dvar)
          Jmin=JstrV
#   ifdef MASKING
          mask => GRID(ng)%vmask
#   endif
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      Nold=1
      Nnew=2
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            i=BOUNDS(ng)%edge(ibry,ctype)
            DO j=Jmin,Jmax
#   ifdef SPLINES_VCONV
              DO k=1,N(ng)
                oHz(j,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
              END DO
#   else
              FC(j,N(ng))=0.0_r8
              DO k=1,N(ng)-1
                FC(j,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-      &
     &                                      GRID(ng)%z_r(i,j,k  ))
              END DO
              FC(j,0)=0.0_r8
#   endif
            END DO
          CASE (isouth, inorth)
            j=BOUNDS(ng)%edge(ibry,ctype)
            DO i=Imin,Imax
#   ifdef SPLINES_VCONV
              DO k=1,N(ng)
                oHz(i,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
              END DO
#   else
              FC(i,N(ng))=0.0_r8
              DO k=1,N(ng)-1
                FC(i,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-      &
     &                                      GRID(ng)%z_r(i,j,k  ))
              END DO
              FC(i,0)=0.0_r8
#   endif
            END DO
        END SELECT
      END IF
!
!  Set initial conditions.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO k=1,N(ng)
              DO j=Jmin,Jmax
!^              Awrk(j,k,Nold)=A(j,k)
!^
                tl_Awrk(j,k,Nold)=tl_A(j,k)
              END DO
            END DO
          CASE (isouth, inorth)
            DO k=1,N(ng)
              DO i=Imin,Imax
!^              Awrk(i,k,Nold)=A(i,k)
!^
                tl_Awrk(i,k,Nold)=tl_A(i,k)
              END DO
            END DO
        END SELECT
      END IF

#   ifdef SPLINES_VCONV
!
!  Integrate vertical diffusion equation implicitly using parabolic
!  splines.
!
      DO step=1,NVsteps
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
        IF (Lboundary(ibry)) THEN
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              cff1=1.0_r8/6.0_r8
              DO j=Jmin,Jmax
                DO k=1,N(ng)-1
                  FC(j,k)=cff1*GRID(ng)%Hz(i,j,k  )-                    &
     &                    DTsizeV*Kv(i,j,k-1)*oHz(j,k  )
                  CF(j,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                    &
     &                    DTsizeV*Kv(i,j,k+1)*oHz(j,k+1)
                END DO
                CF(j,0)=0.0_r8
!^              DC(j,0)=0.0_r8
!^
                tl_DC(j,0)=0.0_r8
              END DO
            CASE (isouth, inorth)
              i=BOUNDS(ng)%edge(ibry,ctype)
              cff1=1.0_r8/6.0_r8
              DO i=Imin,Imax
                DO k=1,N(ng)-1
                  FC(i,k)=cff1*GRID(ng)%Hz(i,j,k  )-                    &
     &                    DTsizeV*Kv(i,j,k-1)*oHz(i,k  )
                  CF(i,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                    &
     &                    DTsizeV*Kv(i,j,k+1)*oHz(i,k+1)
                END DO
                CF(i,0)=0.0_r8
!^              DC(i,0)=0.0_r8
!^
                tl_DC(i,0)=0.0_r8
              END DO
          END SELECT
!
!  LU decomposition and forward substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              cff1=1.0_r8/3.0_r8
              DO k=1,N(ng)-1
                DO j=Jmin,Jmax
                  BC(j,k)=cff1*(GRID(ng)%Hz(i,j,k  )+                   &
     &                          GRID(ng)%Hz(i,j,k+1))+                  &
     &                    DTsizeV*Kv(i,j,k)*                            &
     &                    (oHz(j,k)+oHz(j,k+1))
                  cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
                  CF(j,k)=cff*CF(j,k)
!^                DC(j,k)=cff*(Awrk(j,k+1,Nold)-                        &
!^   &                         Awrk(j,k  ,Nold)-                        &
!^   &                         FC(j,k)*DC(j,k-1))
!^
                  tl_DC(j,k)=cff*(tl_Awrk(j,k+1,Nold)-                  &
     &                            tl_Awrk(j,k  ,Nold)-                  &
     &                            FC(j,k)*tl_DC(j,k-1))
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              cff1=1.0_r8/3.0_r8
              DO k=1,N(ng)-1
                DO i=Imin,Imax
                  BC(i,k)=cff1*(GRID(ng)%Hz(i,j,k  )+                   &
     &                          GRID(ng)%Hz(i,j,k+1))+                  &
     &                    DTsizeV*Kv(i,j,k)*                            &
     &                    (oHz(i,k)+oHz(i,k+1))
                  cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
                  CF(i,k)=cff*CF(i,k)
!^                DC(i,k)=cff*(Awrk(i,k+1,Nold)-                        &
!^   &                         Awrk(i,k  ,Nold)-                        &
!^   &                         FC(i,k)*DC(i,k-1))
!^
                  tl_DC(i,k)=cff*(tl_Awrk(i,k+1,Nold)-                  &
     &                            tl_Awrk(i,k  ,Nold)-                  &
     &                            FC(i,k)*tl_DC(i,k-1))
                END DO
              END DO
          END SELECT
!
!  Backward substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              DO j=Jmin,Jmax
!^              DC(j,N(ng))=0.0_r8
!^
                tl_DC(j,N(ng))=0.0_r8
              END DO
              DO k=N(ng)-1,1,-1
                DO j=Jmin,Jmax
!^                DC(j,k)=DC(j,k)-CF(j,k)*DC(j,k+1)
!^
                  tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
                END DO
              END DO
              DO k=1,N(ng)
                DO j=Jmin,Jmax
!^                DC(j,k)=DC(j,k)*Kv(i,j,k)
!^
                  tl_DC(j,k)=tl_DC(j,k)*Kv(i,j,k)
!^                Awrk(j,k,Nnew)=Awrk(j,k,Nold)+                        &
!^   &                           DTsizeV*oHz(j,k)*                      &
!^   &                           (DC(j,k)-DC(j,k-1))
!^
                  tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                  &
     &                              DTsizeV*oHz(j,k)*                   &
     &                              (tl_DC(j,k)-tl_DC(j,k-1))
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              DO i=Imin,Imax
!^              DC(i,N(ng))=0.0_r8
!^
                tl_DC(i,N(ng))=0.0_r8
              END DO
              DO k=N(ng)-1,1,-1
                DO i=Imin,Imax
!^                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
                  tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
                END DO
              END DO
              DO k=1,N(ng)
                DO i=Imin,Imax
!^                DC(i,k)=DC(i,k)*Kv(i,j,k)
!^
                  tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^                Awrk(i,k,Nnew)=Awrk(i,k,Nold)+                        &
!^   &                           DTsizeV*oHz(i,k)*                      &
!^   &                           (DC(i,k)-DC(i,k-1))
!^
                  tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                  &
     &                              DTsizeV*oHz(i,k)*                   &
     &                              (tl_DC(i,k)-tl_DC(i,k-1))
                END DO
              END DO
          END SELECT
        END IF
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO

#   else
!
!  Integrate vertical diffusion equation implicitly.
!
      DO step=1,NVsteps
!
!  Compute diagonal matrix coefficients BC and load right-hand-side
!  terms for the diffusion equation into DC.
!
        IF (Lboundary(ibry)) THEN
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              DO k=1,N(ng)
                DO j=Jmin,Jmax
                  BC(j,k)=GRID(ng)%Hz(i,j,k)-FC(j,k)-FC(j,k-1)
!^                DC(j,k)=Awrk(j,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
                  tl_DC(j,k)=tl_Awrk(j,k,Nold)*GRID(ng)%Hz(i,j,k)
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              DO k=1,N(ng)
                DO i=Imin,Imax
                  BC(i,k)=GRID(ng)%Hz(i,j,k)-FC(i,k)-FC(i,k-1)
!^                DC(i,k)=Awrk(i,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
                  tl_DC(i,k)=tl_Awrk(i,k,Nold)*GRID(ng)%Hz(i,j,k)
                END DO
              END DO
          END SELECT
!
!  Solve the tridiagonal system.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              DO j=Jmin,Jmax
                cff=1.0_r8/BC(j,1)
                CF(j,1)=cff*FC(j,1)
!^              DC(j,1)=cff*DC(j,1)
!^
                tl_DC(j,1)=cff*tl_DC(j,1)
              END DO
              DO k=2,N(ng)-1
                DO j=Jmin,Jmax
                  cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
                  CF(j,k)=cff*FC(j,k)
!^                DC(j,k)=cff*(DC(j,k)-FC(j,k-1)*DC(j,k-1))
!^
                  tl_DC(j,k)=cff*(tl_DC(j,k)-FC(j,k-1)*tl_DC(j,k-1))
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              DO i=Imin,Imax
                cff=1.0_r8/BC(i,1)
                CF(i,1)=cff*FC(i,1)
!^              DC(i,1)=cff*DC(i,1)
!^
                tl_DC(i,1)=cff*tl_DC(i,1)
              END DO
              DO k=2,N(ng)-1
                DO i=Imin,Imax
                  cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
                  CF(i,k)=cff*FC(i,k)
!^                DC(i,k)=cff*(DC(i,k)-FC(i,k-1)*DC(i,k-1))
!^
                  tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,k-1)*tl_DC(i,k-1))
                END DO
              END DO
          END SELECT
!
!  Compute new solution by back substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              DO j=Jmin,Jmax
!^              DC(j,N(ng))=(DC(j,N(ng))-                               &
!^   &                       FC(j,N(ng)-1)*DC(j,N(ng)-1))/              &
!^   &                      (BC(j,N(ng))-                               &
!^   &                       FC(j,N(ng)-1)*CF(j,N(ng)-1))
!^
                tl_DC(j,N(ng))=(tl_DC(j,N(ng))-                         &
     &                          FC(j,N(ng)-1)*tl_DC(j,N(ng)-1))/        &
     &                         (BC(j,N(ng))-                            &
     &                          FC(j,N(ng)-1)*CF(j,N(ng)-1))
!^              Awrk(j,N(ng),Nnew)=DC(j,N(ng))
!^
                tl_Awrk(j,N(ng),Nnew)=tl_DC(j,N(ng))
#    ifdef MASKING
!^              Awrk(j,N(ng),Nnew)=Awrk(j,N(ng),Nnew)*mask(i,j)
!^
                tl_Awrk(j,N(ng),Nnew)=tl_Awrk(j,N(ng),Nnew)*mask(i,j)
#    endif
              END DO
              DO k=N(ng)-1,1,-1
                DO j=Jmin,Jmax
!^                DC(j,k)=DC(j,k)-CF(j,k)*DC(j,k+1)
!^
                  tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!^                Awrk(j,k,Nnew)=DC(j,k)
!^
                  tl_Awrk(j,k,Nnew)=tl_DC(j,k)
#    ifdef MASKING
!^                Awrk(j,k,Nnew)=Awrk(j,k,Nnew)*mask(i,j)
!^
                  tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nnew)*mask(i,j)
#    endif
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              DO i=Imin,Imax
!^              DC(i,N(ng))=(DC(i,N(ng))-                               &
!^   &                       FC(i,N(ng)-1)*DC(i,N(ng)-1))/              &
!^   &                      (BC(i,N(ng))-                               &
!^   &                       FC(i,N(ng)-1)*CF(i,N(ng)-1))
!^
                tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                         &
     &                          FC(i,N(ng)-1)*tl_DC(i,N(ng)-1))/        &
     &                         (BC(i,N(ng))-                            &
     &                          FC(i,N(ng)-1)*CF(i,N(ng)-1))
!^              Awrk(i,N(ng),Nnew)=DC(i,N(ng))
!^
                tl_Awrk(i,N(ng),Nnew)=tl_DC(i,N(ng))
#    ifdef MASKING
!^              Awrk(i,N(ng),Nnew)=Awrk(i,N(ng),Nnew)*mask(i,j)
!^
                tl_Awrk(i,N(ng),Nnew)=tl_Awrk(i,N(ng),Nnew)*mask(i,j)
#    endif
              END DO
              DO k=N(ng)-1,1,-1
                DO i=Imin,Imax
!^                DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
!^
                  tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^                Awrk(i,k,Nnew)=DC(i,k)
!^
                  tl_Awrk(i,k,Nnew)=tl_DC(i,k)
#    ifdef MASKING
!^                Awrk(i,k,Nnew)=Awrk(i,k,Nnew)*mask(i,j)
!^
                  tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nnew)*mask(i,j)
#    endif
                END DO
              END DO
          END SELECT
        END IF
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
      END DO
#   endif
!
!  Load convolved solution.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            i=BOUNDS(ng)%edge(ibry,ctype)
            DO k=1,N(ng)
              DO j=Jmin,Jmax
!^              A(j,k)=Awrk(j,k,Nold)
!^
                tl_A(j,k)=tl_Awrk(j,k,Nold)
              END DO
            END DO
          CASE (isouth, inorth)
            j=BOUNDS(ng)%edge(ibry,ctype)
            DO k=1,N(ng)
              DO i=Imin,Imax
!^              A(i,k)=Awrk(i,k,Nold)
!^
                tl_A(i,k)=tl_Awrk(i,k,Nold)
              END DO
            END DO
        END SELECT
      END IF

#   ifdef DISTRIBUTE
!
!^    CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij, 1, N(ng),                     &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        A)
!^
      CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
     &                        LBij, UBij, 1, N(ng),                     &
     &                        NghostPoints,                             &
     &                        EWperiodic(ng), NSperiodic(ng),           &
     &                        tl_A)
#   endif
!
      RETURN
      END SUBROUTINE multiscale_bry_Vdiff_tl
!
!-----------------------------------------------------------------------
!  It advances the adjoint vertical pseudo-difusion operator using an
!  implicit solver using parabolic splines or inverting a tridiagonal
!  matrix for adjusting lateral boundaries in the control vector.
!
      SUBROUTINE multiscale_bry_Vdiff_ad (self, ng, tile, model,        &
     &                                    ifield, ibry, ctype, NVsteps, &
     &                                    LBi, UBi, LBj, UBj,           &
     &                                    LBij, UBij,                   &
     &                                    DTsizeV, Kv, ad_A)
!
      CLASS (multiscale), intent(inout) :: self      ! multiscale object
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ifield    ! state field ID 
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: NVsteps   ! integration steps
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      integer,            intent(in   ) :: LBij, UBij
      real (r8),          intent(in   ) :: DTsizeV
      real (r8),          intent(in   ) :: Kv(LBi:,LBj:,0:)
      real (r8),          intent(inout) :: ad_A(LBij:,:)
!
      logical, dimension(4)             :: Lboundary
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: Nnew, Nold, Nsav
      integer                           :: i, j, k, step
!
      real (r8)                         :: adfac, cff, cff1
!
      real (r8), pointer                :: mask(:,:) => NULL()
!
      real(r8), dimension(LBij:UBij,N(ng),2) :: ad_Awrk
#   if defined SPLINES_VCONV
      real(r8), dimension(LBij:UBij,N(ng))   :: oHz
#   endif
      real(r8), dimension(LBij:UBij,0:N(ng)) :: FC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: BC
      real(r8), dimension(LBij:UBij,0:N(ng)) :: CF
      real(r8), dimension(LBij:UBij,0:N(ng)) :: ad_DC
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", multiscale_bry_Vdiff_ad"
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
#   ifdef MASKING
        CASE (r3dvar)
          mask => GRID(ng)%rmask
#   endif
        CASE (u3dvar)
          Imin=IstrU
#   ifdef MASKING
          mask => GRID(ng)%umask
#   endif
        CASE (v3dvar)
          Jmin=JstrV
#   ifdef MASKING
          mask => GRID(ng)%vmask
#   endif
      END SELECT
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      Nold=1
      Nnew=2
!
!  Compute vertical metric factor.  Notice that "z_r" and "Hz" are
!  assumed to be time invariant in the vertical diffusion operator.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            i=BOUNDS(ng)%edge(ibry,ctype)
            DO j=Jmin,Jmax
#   ifdef SPLINES_VCONV
              DO k=1,N(ng)
                oHz(j,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
              END DO
#   else
              FC(j,N(ng))=0.0_r8
              DO k=1,N(ng)-1
                FC(j,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-      &
     &                                      GRID(ng)%z_r(i,j,k  ))
              END DO
              FC(j,0)=0.0_r8
#   endif
            END DO
          CASE (isouth, inorth)
            j=BOUNDS(ng)%edge(ibry,ctype)
            DO i=Imin,Imax
#   ifdef SPLINES_VCONV
              DO k=1,N(ng)
                oHz(i,k)=1.0_r8/GRID(ng)%Hz(i,j,k)
              END DO
#   else
              FC(i,N(ng))=0.0_r8
              DO k=1,N(ng)-1
                FC(i,k)=-DTsizeV*Kv(i,j,k)/(GRID(ng)%z_r(i,j,k+1)-      &
     &                                      GRID(ng)%z_r(i,j,k  ))
              END DO
              FC(i,0)=0.0_r8
#   endif
            END DO
        END SELECT
      END IF
!
!  Adjoint of load convolved solution.
!
#   ifdef DISTRIBUTE
!^    CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij, 1, N(ng),                     &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        tl_A
!^
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ibry,              &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
#   endif

      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            i=BOUNDS(ng)%edge(ibry,ctype)
            DO k=1,N(ng)
              DO j=Jmin,Jmax
!^              tl_A(j,k)=tl_Awrk(j,k,Nold)
!^
                ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                    &
     &                            ad_A(j,k)
                ad_A(j,k)=0.0_r8
              END DO
            END DO
          CASE (isouth, inorth)
            j=BOUNDS(ng)%edge(ibry,ctype)
            DO k=1,N(ng)
              DO i=Imin,Imax
!^              tl_A(i,k)=tl_Awrk(i,k,Nold)
!^
                ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                    &
     &                            ad_A(i,k)
                ad_A(i,k)=0.0_r8
              END DO
            END DO
        END SELECT
      END IF

#   ifdef SPLINES_VCONV
!
!  Integrate adjoint vertical diffusion equation implicitly using
!  parabolic splines.
!
      STEP_LOOP : DO step=1,NVsteps
!
!  Update integration indices.
!
        Nsav=Nold
        Nold=Nnew
        Nnew=Nsav
!
!  Use conservative, parabolic spline reconstruction of vertical
!  diffusion derivatives.  Then, time step vertical diffusion term
!  implicitly.
!
!  Compute basic state time-invariant coefficients.
!
        IF (Lboundary(ibry)) THEN
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              cff1=1.0_r8/6.0_r8
              DO j=Jmin,Jmax
                DO k=1,N(ng)-1
                  FC(j,k)=cff1*GRID(ng)%Hz(i,j,k  )-                    &
     &                    DTsizeV*Kv(i,j,k-1)*oHz(j,k  )
                  CF(j,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                    &
     &                    DTsizeV*Kv(i,j,k+1)*oHz(j,k+1)
                END DO
                CF(j,0)=0.0_r8
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              cff1=1.0_r8/6.0_r8
              DO i=Imin,Imax
                DO k=1,N(ng)-1
                  FC(i,k)=cff1*GRID(ng)%Hz(i,j,k  )-                    &
     &                    DTsizeV*Kv(i,j,k-1)*oHz(i,k  )
                  CF(i,k)=cff1*GRID(ng)%Hz(i,j,k+1)-                    &
     &                    DTsizeV*Kv(i,j,k+1)*oHz(i,k+1)
                END DO
                CF(i,0)=0.0_r8
              END DO
          END SELECT
!
!  Adjoint of LU decomposition and forward substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              cff1=1.0_r8/3.0_r8
              DO k=1,N(ng)-1
                DO j=Jmin,Jmax
                  BC(j,k)=cff1*(GRID(ng)%Hz(i,j,k  )+                   &
     &                          GRID(ng)%Hz(i,j,k+1))+                  &
     &                    DTsizeV*Kv(i,j,k)*                            &
     &                    (oHz(j,k)+oHz(j,k+1))
                  cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
                  CF(j,k)=cff*CF(j,k)
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              cff1=1.0_r8/3.0_r8
              DO k=1,N(ng)-1
                DO i=Imin,Imax
                  BC(i,k)=cff1*(GRID(ng)%Hz(i,j,k  )+                   &
     &                          GRID(ng)%Hz(i,j,k+1))+                  &
     &                    DTsizeV*Kv(i,j,k)*                            &
     &                    (oHz(i,k)+oHz(i,k+1))
                  cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
                  CF(i,k)=cff*CF(i,k)
                END DO
              END DO
          END SELECT
!
!  Adjoint backward substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              DO k=1,N(ng)
                DO j=Jmin,Jmax
!^                tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nold)+                  &
!^   &                              DTsizeV*oHz(j,k)*                   &
!^   &                              (tl_DC(j,k)-tl_DC(j,k-1))
!^
                  adfac=DTsizeV*oHz(j,k)*ad_Awrk(j,k,Nnew)
                  ad_DC(j,k-1)=ad_DC(j,k-1)-adfac
                  ad_DC(j,k  )=ad_DC(j,k  )+adfac
                  ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                  &
     &                              ad_Awrk(j,k,Nnew)
                  ad_Awrk(j,k,Nnew)=0.0_r8
!^                tl_DC(j,k)=tl_DC(j,k)*Kv(i,j,k)
!^
                  ad_DC(j,k)=ad_DC(j,k)*Kv(i,j,k)
                END DO
              END DO
              DO k=1,N(ng)-1
                DO j=Jmin,Jmax
!^                tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!^
                  ad_DC(j,k+1)=ad_DC(j,k+1)-CF(j,k)*ad_DC(j,k)
                END DO
              END DO
              DO j=Jmin,Jmax
!^              tl_DC(j,N(ng))=0.0_r8
!^
                ad_DC(j,N(ng))=0.0_r8
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              DO k=1,N(ng)
                DO i=Imin,Imax
!^                tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nold)+                  &
!^   &                              DTsizeV*oHz(i,k)*                   &
!^   &                              (tl_DC(i,k)-tl_DC(i,k-1))
!^
                  adfac=DTsizeV*oHz(i,k)*ad_Awrk(i,k,Nnew)
                  ad_DC(i,k-1)=ad_DC(i,k-1)-adfac
                  ad_DC(i,k  )=ad_DC(i,k  )+adfac
                  ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                  &
     &                              ad_Awrk(i,k,Nnew)
                  ad_Awrk(i,k,Nnew)=0.0_r8
!^                tl_DC(i,k)=tl_DC(i,k)*Kv(i,j,k)
!^
                  ad_DC(i,k)=ad_DC(i,k)*Kv(i,j,k)
                END DO
              END DO
              DO k=1,N(ng)-1
                DO i=Imin,Imax
!^                tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
                  ad_DC(i,k+1)=ad_DC(i,k+1)-CF(i,k)*ad_DC(i,k)
                END DO
              END DO
              DO i=Imin,Imax
!^              tl_DC(i,N(ng))=0.0_r8
!^
                ad_DC(i,N(ng))=0.0_r8
              END DO
          END SELECT
!
!  Adjoint of LU decomposition and forward substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              DO k=N(ng)-1,1,-1
                DO j=Jmin,Jmax
                  cff=1.0_r8/(BC(j,k)-FC(j,k)*CF(j,k-1))
!^                tl_DC(j,k)=cff*(tl_Awrk(j,k+1,Nold)-                  &
!^   &                            tl_Awrk(j,k  ,Nold)-                  &
!^   &                            FC(j,k)*tl_DC(j,k-1))
!^
                  adfac=cff*ad_DC(j,k)
                  ad_Awrk(j,k  ,Nold)=ad_Awrk(j,k  ,Nold)-adfac
                  ad_Awrk(j,k+1,Nold)=ad_Awrk(j,k+1,Nold)+adfac
                  ad_DC(j,k-1)=ad_DC(j,k-1)-FC(j,k)*adfac
                  ad_DC(j,k)=0.0_r8
                END DO
              END DO
              DO j=Jmin,Jmax
!^              tl_DC(j,0)=0.0_r8
!^
                ad_DC(j,0)=0.0_r8
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              DO k=N(ng)-1,1,-1
                DO i=Imin,Imax
                  cff=1.0_r8/(BC(i,k)-FC(i,k)*CF(i,k-1))
!^                tl_DC(i,k)=cff*(tl_Awrk(i,k+1,Nold)-                  &
!^   &                            tl_Awrk(i,k  ,Nold)-                  &
!^   &                            FC(i,k)*tl_DC(i,k-1))
!^
                  adfac=cff*ad_DC(i,k)
                  ad_Awrk(i,k  ,Nold)=ad_Awrk(i,k  ,Nold)-adfac
                  ad_Awrk(i,k+1,Nold)=ad_Awrk(i,k+1,Nold)+adfac
                  ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k)*adfac
                  ad_DC(i,k)=0.0_r8
                END DO
              END DO
              DO i=Imin,Imax
!^              tl_DC(i,0)=0.0_r8
!^
                ad_DC(i,0)=0.0_r8
              END DO
          END SELECT
        END IF
      END DO STEP_LOOP

#   else
!
!  Integrate adjoint vertical diffusion equation implicitly.
!
      STEP_LOOP : DO step=1,NVsteps
!
!  Adjoint of update integration indices.
!
        Nsav=Nnew
        Nnew=Nold
        Nold=Nsav
!
!  Compute diagonal matrix coefficients BC and CF.
!
        IF (Lboundary(ibry)) THEN
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              DO k=1,N(ng)
                DO j=Jmin,Jmax
                  BC(j,k)=GRID(ng)%Hz(i,j,k)-FC(j,k)-FC(j,k-1)
                END DO
              END DO
              DO j=Jmin,Jmax
                cff=1.0_r8/BC(j,1)
                CF(j,1)=cff*FC(j,1)
              END DO
              DO k=2,N(ng)-1
                DO j=Jmin,Jmax
                  cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
                  CF(j,k)=cff*FC(j,k)
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              DO k=1,N(ng)
                DO i=Imin,Imax
                  BC(i,k)=GRID(ng)%Hz(i,j,k)-FC(i,k)-FC(i,k-1)
                END DO
              END DO
              DO i=Imin,Imax
                cff=1.0_r8/BC(i,1)
                CF(i,1)=cff*FC(i,1)
              END DO
              DO k=2,N(ng)-1
                DO i=Imin,Imax
                  cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
                  CF(i,k)=cff*FC(i,k)
                END DO
              END DO
          END SELECT
!
!  Adjoint of compute new solution by back substitution.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
!^            DO k=N(ng)-1,1,-1
!^
              DO k=1,N(ng)-1
                DO j=Jmin,Jmax
#    ifdef MASKING
!^                tl_Awrk(j,k,Nnew)=tl_Awrk(j,k,Nnew)*mask(i,j)
!^
                  ad_Awrk(j,k,Nnew)=ad_Awrk(j,k,Nnew)*mask(i,j)
#    endif
!^                tl_Awrk(j,k,Nnew)=tl_DC(j,k)
!^
                  ad_DC(j,k)=ad_DC(j,k)+                                &
     &                       ad_Awrk(j,k,Nnew)
                  ad_Awrk(j,k,Nnew)=0.0_r8
!^                tl_DC(j,k)=tl_DC(j,k)-CF(j,k)*tl_DC(j,k+1)
!^
                  ad_DC(j,k+1)=-CF(j,k)*ad_DC(j,k)
                END DO
              END DO
              DO j=Jmin,Jmax
#    ifdef MASKING
!^              tl_Awrk(j,N(ng),Nnew)=tl_Awrk(j,N(ng),Nnew)*mask(i,j)
!^
                ad_Awrk(j,N(ng),Nnew)=ad_Awrk(j,N(ng),Nnew)*mask(i,j)

#    endif
!^              tl_Awrk(j,N(ng),Nnew)=tl_DC(j,N(ng))
!^
                ad_DC(j,N(ng))=ad_DC(j,N(ng))+                          &
     &                         ad_Awrk(j,N(ng),Nnew)
                ad_Awrk(j,N(ng),Nnew)=0.0_r8
!^              tl_DC(j,N(ng))=(tl_DC(j,N(ng))-                         &
!^   &                          FC(j,N(ng)-1)*tl_DC(j,N(ng)-1))/        &
!^   &                         (BC(j,N(ng))-                            &
!^   &                          FC(j,N(ng)-1)*CF(j,N(ng)-1))
!^
                adfac=ad_DC(j,N(ng))/                                   &
     &                (BC(j,N(ng))-                                     &
     &                 FC(j,N(ng)-1)*CF(j,N(ng)-1))
                ad_DC(j,N(ng)-1)=ad_DC(j,N(ng)-1)-                      &
     &                           FC(j,N(ng)-1)*adfac
                ad_DC(j,N(ng)  )=adfac
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
!^            DO k=N(ng)-1,1,-1
!^
              DO k=1,N(ng)-1
                DO i=Imin,Imax
#    ifdef MASKING
!^                tl_Awrk(i,k,Nnew)=tl_Awrk(i,k,Nnew)*mask(i,j)
!^
                  ad_Awrk(i,k,Nnew)=ad_Awrk(i,k,Nnew)*mask(i,j)
#    endif
!^                tl_Awrk(i,k,Nnew)=tl_DC(i,k)
!^
                  ad_DC(i,k)=ad_DC(i,k)+                                &
     &                       ad_Awrk(i,k,Nnew)
                  ad_Awrk(i,k,Nnew)=0.0_r8
!^                tl_DC(i,k)=tl_DC(i,k)-CF(i,k)*tl_DC(i,k+1)
!^
                  ad_DC(i,k+1)=-CF(i,k)*ad_DC(i,k)
                END DO
              END DO
              DO i=Imin,Imax
#    ifdef MASKING
!^              tl_Awrk(i,N(ng),Nnew)=tl_Awrk(i,N(ng),Nnew)*mask(i,j)
!^
                ad_Awrk(i,N(ng),Nnew)=ad_Awrk(i,N(ng),Nnew)*mask(i,j)
#    endif
!^              tl_Awrk(i,N(ng),Nnew)=tl_DC(i,N(ng))
!^
                ad_DC(i,N(ng))=ad_DC(i,N(ng))+                          &
     &                         ad_Awrk(i,N(ng),Nnew)
                ad_Awrk(i,N(ng),Nnew)=0.0_r8
!^              tl_DC(i,N(ng))=(tl_DC(i,N(ng))-                         &
!^   &                          FC(i,N(ng)-1)*tl_DC(i,N(ng)-1))/        &
!^   &                         (BC(i,N(ng))-                            &
!^   &                          FC(i,N(ng)-1)*CF(i,N(ng)-1))
!^
                adfac=ad_DC(i,N(ng))/                                   &
     &                (BC(i,N(ng))-                                     &
     &                 FC(i,N(ng)-1)*CF(i,N(ng)-1))
                ad_DC(i,N(ng)-1)=ad_DC(i,N(ng)-1)-                      &
     &                           FC(i,N(ng)-1)*adfac
                ad_DC(i,N(ng)  )=adfac
              END DO
            END SELECT
!
!  Solve the adjoint tridiagonal system.
!
           SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
!^            DO k=2,N(ng)-1
!^
              DO k=N(ng)-1,2,-1
                DO j=Jmin,Jmax
                  cff=1.0_r8/(BC(j,k)-FC(j,k-1)*CF(j,k-1))
!^                tl_DC(j,k)=cff*(tl_DC(j,k)-FC(j,k-1)*tl_DC(j,k-1))
!^
                  adfac=cff*ad_DC(j,k)
                  ad_DC(j,k-1)=ad_DC(j,k-1)-FC(j,k-1)*adfac
                  ad_DC(j,k  )=adfac
                END DO
              END DO
              DO j=Jmin,Jmax
                cff=1.0_r8/BC(j,1)
!^              tl_DC(j,1)=cff*tl_DC(j,1)
!^
                ad_DC(j,1)=cff*ad_DC(j,1)
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
!^            DO k=2,N(ng)-1
!^
              DO k=N(ng)-1,2,-1
                DO i=Imin,Imax
                  cff=1.0_r8/(BC(i,k)-FC(i,k-1)*CF(i,k-1))
!^                tl_DC(i,k)=cff*(tl_DC(i,k)-FC(i,k-1)*tl_DC(i,k-1))
!^
                  adfac=cff*ad_DC(i,k)
                  ad_DC(i,k-1)=ad_DC(i,k-1)-FC(i,k-1)*adfac
                  ad_DC(i,k  )=adfac
                END DO
              END DO
              DO i=Imin,Imax
                cff=1.0_r8/BC(i,1)
!^              tl_DC(i,1)=cff*tl_DC(i,1)
!^
                ad_DC(i,1)=cff*ad_DC(i,1)
              END DO
          END SELECT
!
!  Adjoint of right-hand-side terms for the diffusion equation.
!
          SELECT CASE (ibry)
            CASE (iwest, ieast)
              i=BOUNDS(ng)%edge(ibry,ctype)
              DO k=1,N(ng)
                DO j=Jmin,Jmax
!^                tl_DC(j,k)=tl_Awrk(j,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
                  ad_Awrk(j,k,Nold)=ad_Awrk(j,k,Nold)+                  &
     &                              GRID(ng)%Hz(i,j,k)*ad_DC(j,k)
                  ad_DC(j,k)=0.0_r8
                END DO
              END DO
            CASE (isouth, inorth)
              j=BOUNDS(ng)%edge(ibry,ctype)
              DO k=1,N(ng)
                DO i=Imin,Imax
!^                tl_DC(i,k)=tl_Awrk(i,k,Nold)*GRID(ng)%Hz(i,j,k)
!^
                  ad_Awrk(i,k,Nold)=ad_Awrk(i,k,Nold)+                  &
     &                              GRID(ng)%Hz(i,j,k)*ad_DC(i,k)
                  ad_DC(i,k)=0.0_r8
                END DO
              END DO
          END SELECT
        END IF
      END DO STEP_LOOP
#   endif
!
!  Set adjoint initial conditions.
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            i=BOUNDS(ng)%edge(ibry,ctype)
            DO k=1,N(ng)
              DO j=Jmin,Jmax
!^              tl_Awrk(j,k,Nold)=tl_A(j,k)
!^
                ad_A(j,k)=ad_A(j,k)+ad_Awrk(j,k,Nold)
                ad_Awrk(j,k,Nold)=0.0_r8
              END DO
            END DO
          CASE (isouth, inorth)
            j=BOUNDS(ng)%edge(ibry,ctype)
            DO k=1,N(ng)
              DO i=Imin,Imax
!^              tl_Awrk(i,k,Nold)=tl_A(i,k)
!^
                ad_A(i,k)=ad_A(i,k)+ad_Awrk(i,k,Nold)
                ad_Awrk(i,k,Nold)=0.0_r8
              END DO
            END DO
        END SELECT
      END IF

#   ifdef DISTRIBUTE
!
!^    CALL mp_exchange3d_bry (ng, tile, model, 1, ibry,                 &
!^   &                        LBij, UBij, 1, N(ng),                     &
!^   &                        NghostPoints,                             &
!^   &                        EWperiodic(ng), NSperiodic(ng),           &
!^   &                        tl_A)
!^
      CALL ad_mp_exchange3d_bry (ng, tile, model, 1, ibry,              &
     &                           LBij, UBij, 1, N(ng),                  &
     &                           NghostPoints,                          &
     &                           EWperiodic(ng), NSperiodic(ng),        &
     &                           ad_A)
#   endif
!
      RETURN
      END SUBROUTINE multiscale_bry_Vdiff_ad
#  endif /* SOLVE3D */
# endif /* ADJUST_BOUNDARY */

# ifdef MULTI_SCALE_DEBUG
!
!-----------------------------------------------------------------------
!  It prints valuable diagnostic information about the convergence of
!  the CG and CI algorithms for tuning the implicit multiscale diffusion
!  operator.
!  
      FUNCTION multiscale_print (ng, ifield, ms, iterDiff, Ndiff,       &
     &                           iter, Niter, label, deps,              &
     &                           level) RESULT (status)
!
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: ifield    ! state field ID
      integer,            intent(in   ) :: ms        ! multiscale index
      integer,            intent(in   ) :: iterDiff  ! diffusion step
      integer,            intent(in   ) :: Ndiff     ! diff applications
      integer,            intent(in   ) :: iter      ! solver iteration
      integer,            intent(in   ) :: Niter     ! CG/CI iterations
      integer, optional,  intent(in   ) :: level     ! vertical level
      real (r8),          intent(in   ) :: deps      ! convergence norm
      character (len=*),  intent(in   ) :: label     ! solver label
!
      integer, parameter                :: iunit = 45
      integer                           :: status
!
!  Report CG or CI algorithm convergence, normalized ratio of initial
!  and iterative inner products.
!
      status=0
!
      IF (Master) THEN
        IF ((iterDiff.eq.1).and.(iter.eq.0)) THEN
          WRITE (iunit,10) label, TRIM(Vname(1,idSvar(ifield))), ms
        END IF
        IF (iter.eq.Niter) THEN
          IF (PRESENT(level)) THEN
            WRITE (iunit,20) iterDiff, iter, deps, level
          ELSE
            WRITE (iunit,30) iterDiff, iter, deps
          END IF
        END IF
        IF (iterDiff.eq.Ndiff) THEN
          IF (PRESENT(level)) THEN
            WRITE (iunit,40) iterDiff, iter, deps, level
          ELSE
            WRITE (iunit,50) iterDiff, iter, deps
          END IF
          FLUSH (iunit)
        END IF
      END IF
!
  10  FORMAT (a,': Convergence for ''',a,'', ', B-MultiScale = ',i0,/)
  20  FORMAT (7x,'K-LAP  => iterDiff = ',i3.3,2x,'iter = ',i3.3,2x,     &
     &        'deps = ',1p,e12.5,2x,'level = ',i0) 
  30  FORMAT (7x,'K-LAP  => iterDiff = ',i3.3,2x,'iter = ',i3.3,2x,     &
     &        'deps = ',1p,e12.5) 
  40  FORMAT (7x,'K-DIFF => iterDiff = ',i3.3,2x,'iter = ',i3.3,2x,     &
     &        'deps = ',1p,e12.5,2x,'level = ',i0) 
  50  FORMAT (7x,'K-DIFF => iterDiff = ',i3.3,2x,'iter = ',i3.3,2x,     &
     &        'deps = ',1p,e12.5) 
!
      RETURN
      END FUNCTION multiscale_print
# endif

# ifdef ADJUST_BOUNDARY
!
!-----------------------------------------------------------------------
!  It computes the dot product between two 1D tiled boundary vectors.
!  
      FUNCTION dot_prod_1d (ng, tile, model, ctype, ibry,               &
     &                      LBij, UBij, B1, B2) RESULT (DotProd)
!
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: ibry      ! boundary ID
      integer,            intent(in   ) :: LBij, UBij
      real (r8),          intent(in   ) :: B1(LBij:)
      real (r8),          intent(in   ) :: B2(LBij:)
!
      logical, dimension(4)             :: Lboundary
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: NSUB, i, j
      real (r8)                         :: DotProd
      real (r8)                         :: cff, my_DotProd
#  ifdef DISTRIBUTE
      character (len=3)                 :: op_handle
#  endif
!
!  Initialize.
!
      Lboundary(iwest )=DOMAIN(ng)%Western_Edge (tile)
      Lboundary(ieast )=DOMAIN(ng)%Eastern_Edge (tile)
      Lboundary(isouth)=DOMAIN(ng)%Southern_Edge(tile)
      Lboundary(inorth)=DOMAIN(ng)%Northern_Edge(tile)
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar, u3dvar)
          Imin=IstrU
        CASE (v2dvar, v3dvar)
          Jmin=JstrV
      END SELECT
!
!  Compute dot product between A1 and A2 arrays.
!
      my_DotProd=0.0_r8
!
      IF (Lboundary(ibry)) THEN
        SELECT CASE (ibry)
          CASE (iwest, ieast)
            DO j=Jmin,Jmax
              cff=B1(j)*B2(j)
              my_DotProd=my_DotProd+cff
            END DO
          CASE (isouth, inorth)
            DO i=Imin,Imax
              cff=B1(i)*B2(i)
              my_DotProd=my_DotProd+cff
            END DO
        END SELECT
      END IF
!
!  Perform parallel global reduction operation.
!
#  ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
#  else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
#  endif
!$OMP CRITICAL (DOT_PROD)
      IF (tile_count.eq.0) THEN
        DotProd=0.0_r8
      END IF
      DotProd=DotProd+my_DotProd
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#  ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, DotProd, op_handle)
#  endif
      END IF
!$OMP END CRITICAL (DOT_PROD)
!
      RETURN
      END FUNCTION dot_prod_1d
# endif
!
!-----------------------------------------------------------------------
!  It compute the dot product between two 2D tiled arrays.
!  
      FUNCTION dot_prod_2d (ng, tile, model, ctype,                     &
     &                      LBi, UBi, LBj, UBj, A1, A2) RESULT (DotProd)
!
      integer,            intent(in   ) :: ng        ! nested grid
      integer,            intent(in   ) :: tile      ! domain partition
      integer,            intent(in   ) :: model     ! kernel ID
      integer,            intent(in   ) :: ctype     ! C-grid type
      integer,            intent(in   ) :: LBi, UBi, LBj, UBj
      real (r8),          intent(in   ) :: A1(LBi:,LBj:)
      real (r8),          intent(in   ) :: A2(LBi:,LBj:)
!
      integer                           :: Istr, IstrU, Iend, Imin, Imax
      integer                           :: Jstr, JstrV, Jend, Jmin, Jmax
      integer                           :: NSUB, i, j
      real (r8)                         :: DotProd
      real (r8)                         :: cff, my_DotProd
# ifdef DISTRIBUTE
      character (len=3)                 :: op_handle
# endif
!
!  Initialize.
!
      Istr =BOUNDS(ng)%Istr (tile)   ! tile computational range
      IstrU=BOUNDS(ng)%IstrU(tile)
      Iend =BOUNDS(ng)%Iend (tile)
      JstrV=BOUNDS(ng)%JstrV(tile)
      Jend =BOUNDS(ng)%Jend (tile)
!
      Imin=Istr
      Imax=Iend
      Jmin=Jstr
      Jmax=Jend
      SELECT CASE (ctype)
        CASE (u2dvar, u3dvar)
          Imin=IstrU
        CASE (v2dvar, v3dvar)
          Jmin=JstrV
      END SELECT
!
!  Compute dot product between A1 and A2 arrays.
!
      my_DotProd=0.0_r8
      DO j=Jmin,Jmax
        DO i=Imin,Imax
          cff=A1(i,j)*A2(i,j)
          my_DotProd=my_DotProd+cff
        END DO
      END DO
!
!  Perform parallel global reduction operation.
!
# ifdef DISTRIBUTE
      NSUB=1                             ! distributed-memory
# else
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
# endif
!$OMP CRITICAL (DOT_PROD)
      IF (tile_count.eq.0) THEN
        DotProd=0.0_r8
      END IF
      DotProd=DotProd+my_DotProd
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        op_handle='SUM'
        CALL mp_reduce (ng, model, 1, DotProd, op_handle)
# endif
      END IF
!$OMP END CRITICAL (DOT_PROD)
!
      RETURN
      END FUNCTION dot_prod_2d
#endif /* MULTI_SCALE_B */
      END MODULE roms_multiscale_mod
