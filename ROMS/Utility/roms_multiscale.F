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
      USE ad_bc_2d_mod,     ONLY : ad_dabc_r2d_tile,                    &
     &                             ad_dabc_u2d_tile,                    &
     &                             ad_dabc_v2d_tile
# ifdef ADJUST_BOUNDARY
      USE ad_bc_bry2d_mod,  ONLY : ad_bc_r2d_bry_tile,                  &
     &                             ad_bc_u2d_bry_tile,                  &
     &                             ad_bc_v2d_bry_tile
# endif
# ifdef SOLVE3D
      USE ad_bc_3d_mod,     ONLY : ad_dabc_r3d_tile,                    &
     &                             ad_dabc_u3d_tile,                    &
     &                             ad_dabc_v3d_tile
#  ifdef ADJUST_BOUNDARY
      USE ad_bc_bry3d_mod,  ONLY : ad_bc_r3d_bry_tile,                  &
     &                             ad_bc_u3d_bry_tile,                  &
     &                             ad_bc_v3d_bry_tile
#  endif
# endif
      USE bc_2d_mod,        ONLY : dabc_r2d_tile,                       &
     &                             dabc_u2d_tile,                       &
     &                             dabc_v2d_tile
# ifdef ADJUST_BOUNDARY
      USE bc_bry2d_mod,     ONLY : bc_r2d_bry_tile,                     &
     &                             bc_u2d_bry_tile,                     &
     &                             bc_v2d_bry_tile
# endif
# ifdef SOLVE3D
      USE bc_3d_mod,        ONLY : dabc_r3d_tile,                       &
     &                             dabc_u3d_tile,                       &
     &                             dabc_v3d_tile
#  ifdef ADJUST_BOUNDARY
      USE bc_bry3d_mod,     ONLY : bc_r3d_bry_tile,                     &
     &                             bc_u3d_bry_tile,                     &
     &                             bc_v3d_bry_tile
#  endif
# endif
      USE lapack_mod,       ONLY : DSTEQR
# ifdef DISTRIBUTE
      USE mp_exchange_mod,  ONLY : ad_mp_exchange2d
#  ifdef ADJUST_BOUNDARY
      USE mp_exchange_mod,  ONLY : ad_mp_exchange2d_bry
#  endif
#  ifdef SOLVE3D
      USE mp_exchange_mod,  ONLY : ad_mp_exchange3d
#   ifdef ADJUST_BOUNDARY
      USE mp_exchange_mod,  ONLY : ad_mp_exchange3d_bry
#   endif

#  endif
      USE mp_exchange_mod,  ONLY : mp_exchange2d
#  ifdef ADJUST_BOUNDARY
      USE mp_exchange_mod,  ONLY : mp_exchange2d_bry
#  endif
#  ifdef SOLVE3D
      USE mp_exchange_mod,  ONLY : mp_exchange3d
#   ifdef ADJUST_BOUNDARY
      USE mp_exchange_mod,  ONLY : mp_exchange3d_bry
#   endif
#  endif
      USE distribute_mod,   ONLY : mp_reduce
# endif
      USE strings_mod,      ONLY : FoundError
!
      implicit none
!
!-----------------------------------------------------------------------
!  ROMS multiscale background-error covariance CLASS object.
!-----------------------------------------------------------------------
!
      TYPE, PUBLIC :: multiscale
!
!  Kernel variables lower and upper bound array dimensions, and parallel
!  tile partition indices.
!
        integer :: tile                           ! domain partition
        integer :: LBi,  UBi,  LBj,  UBj          ! array bounds
# ifdef ADJUST_BOUNDARY
        integer :: LBij, UBij                     ! boundary bounds
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
        integer              :: Nscale
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
!  Implicit Chebyshev Iterations (CI) solver.
!
        PROCEDURE :: tl_CI_2d      => multiscale_CI_2d_tl
        PROCEDURE :: ad_CI_2d      => multiscale_CI_2d_ad
# ifdef SOLVE3D
        PROCEDURE :: tl_CI_3d      => multiscale_CI_3d_tl
        PROCEDURE :: ad_CI_3d      => multiscale_CI_3d_ad
# endif
# ifdef ADJUST_BOUNDARY
!
        PROCEDURE :: tl_CI_b1d     => multiscale_CI_b1d_tl
        PROCEDURE :: ad_CI_b1d     => multiscale_CI_b1d_ad
#  ifdef SOLVE3D
        PROCEDURE :: tl_CI_b2d     => multiscale_CI_b2d_tl
        PROCEDURE :: ad_CI_b2d     => multiscale_CI_b2d_ad
#  endif
# endif
!
!  Implicit Conjugate Gradient (CG) solver.
!
        PROCEDURE :: tl_CG_2d      => multiscale_CG_2d_tl
# ifdef SOLVE3D
        PROCEDURE :: tl_CG_3d      => multiscale_CG_3d_tl
# endif
# ifdef ADJUST_BOUNDARY
!
        PROCEDURE :: tl_CG_b1d     => multiscale_CG_b1d_tl
#  ifdef SOLVE3D
        PROCEDURE :: tl_CG_b2d     => multiscale_CG_b2d_tl
#  endif
# endif

# ifdef SOLVE3D
!
!  Implicit Vertical Diffusion solver.
!
        PROCEDURE :: tl_Vdiff      => multiscale_Vdiff_tl
        PROCEDURE :: ad_Vdiff      => multiscale_Vdiff_ad
#  ifdef ADJUST_BOUNDARY
!
        PROCEDURE :: tl_bry_Vdiff  => multiscale_bry_Vdiff_tl
        PROCEDURE :: ad_bry_Vdiff  => multiscale_bry_Vdiff_ad
#  endif
# endif
!
!  K-Laplacian Operator, [1 + Del(K*Del)].
!
        PROCEDURE :: tl_Klap_r2d   => multiscale_Klap_r2d_tl
        PROCEDURE :: ad_Klap_r2d   => multiscale_Klap_r2d_ad
        PROCEDURE :: tl_Klap_u2d   => multiscale_Klap_u2d_tl
        PROCEDURE :: ad_Klap_u2d   => multiscale_Klap_u2d_ad
        PROCEDURE :: tl_Klap_v2d   => multiscale_Klap_v2d_tl
        PROCEDURE :: ad_Klap_v2d   => multiscale_Klap_v2d_ad

# ifdef SOLVE3D
!
        PROCEDURE :: tl_Klap_r3d   => multiscale_Klap_r3d_tl
        PROCEDURE :: ad_Klap_r3d   => multiscale_Klap_r3d_ad
        PROCEDURE :: tl_Klap_u3d   => multiscale_Klap_u3d_tl
        PROCEDURE :: ad_Klap_u3d   => multiscale_Klap_u3d_ad
        PROCEDURE :: tl_Klap_v3d   => multiscale_Klap_v3d_tl
        PROCEDURE :: ad_Klap_v3d   => multiscale_Klap_v3d_ad
# endif
# ifdef ADJUST_BOUNDARY
!
        PROCEDURE :: tl_Klap_b1d   => multiscale_Klap_b1d_tl
        PROCEDURE :: ad_Klap_b1d   => multiscale_Klap_b1d_ad
#  ifdef SOLVE3D
!
        PROCEDURE :: tl_Klap_b2d   => multiscale_Klap_b2d_tl
        PROCEDURE :: ad_Klap_b2d   => multiscale_Klap_b2d_ad
#  endif
# endif
!
      END TYPE multiscale
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      TYPE (multiscale), allocatable :: B_ms(:)
!
      PUBLIC  :: multiscale_initialize
      PUBLIC  :: multiscale_finalize
!
# ifdef MULTI_SCALE_DEBUG
      PRIVATE :: multiscale_print
# endif
# ifdef ADJUST_BOUNDARY
      PRIVATE :: dot_prod_1d
# endif
      PRIVATE :: dot_prod_2d
!
      CONTAINS
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  It allocates and initializes the background-error covariance (B)
!  multiscale object.
!
     SUBROUTINE multiscale_initialize
!
      integer :: ng, tile
!
!  Allocate and create B-multiscale object.
!
      IF (.not.allocated(B_ms)) THEN
        allocate ( B_ms(Ngrids) )
      END IF
!
      DO ng=1,Ngrids
        DO tile=first_tile(ng),last_tile(ng),+1
          CALL multiscale_create (B_ms(ng), ng, tile)
        END DO
      END DO 
!
      RETURN
      END SUBROUTINE multiscale_initialize
!
!-----------------------------------------------------------------------
!  It deallocates the multiscale background-error covariance object.
!
      SUBROUTINE multiscale_finalize
!
      integer :: ng, tile
!
!  Deallocate background-error ocovariance derived-type structure.
!  Note that it also deallocates its internal arrays.
!
      IF (allocated(B_ms)) deallocate (B_ms)
!
      RETURN
      END SUBROUTINE multiscale_finalize
!
!-----------------------------------------------------------------------
!  It allocates and initializes the multiscale background-error
!  covariance object arrays used by the implicit diffusion operators.
!
      SUBROUTINE multiscale_create (self, ng, tile)
!
      CLASS (multiscale), intent(inout) :: self     ! multiscale object
      integer,            intent(in   ) :: ng       ! nested grid
      integer,            intent(in   ) :: tile     ! domain partition
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
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! <><><> CLASS MULTISCALE:  Chebyshev Iterations Solver <><><><><><><><>
!
# include "multiscale_CIsolver.h"

!
! <><><> CLASS MULTISCALE:  Conjugate Gradient Solver <><><><><><><><><>
!
# include "multiscale_CGsolver.h"

# ifdef SOLVE3D
!
! <><><> CLASS MULTISCALE:  Vertical Diffusion Solver <><><><><><><><><>
!
#  include "multiscale_Vdiff.h"
# endif

!
! <><><> CLASS MULTISCALE:  K-Laplacian Operators <><><><><><><><><><><>
!
# include "multiscale_Klaplacian.h"

!
!::::::: PRIVATE FUNCTIONS :::::::::::::::::::::::::::::::::::::::::::::

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
