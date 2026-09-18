#include "cppdefs.h"
      MODULE lapack_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2026 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  This module includes modernized versions of the selected routines   !
!  from thw Linear Algebra Package (LAPACK) library, which are used    !
!  in ROMS 4D-Var algorithms.                                          !
!                                                                      !
!  Adapted from LAPACK library version 2.0                             !
!                                                                      !
!  NOTES:                                                              !
!                                                                      !
!   - The LAPACK library was written originally in Fortran-77 and has  !
!     been extensively tested and used in various compilers and        !
!     applications. We are modernizing a few of the functions used by  !
!     ROMS at the request of NOAA and NCEP to remove the undesirable   !
!     GOTOs statements with the modern capabilities of the Fortran     !
!     standard (1995 and 2003).                                        !
!                                                                      !
!   - We need to declare the non-scalar arguments to routines as       !
!     assumed-size array, A(*), instead of assumed-shape array, A(:).  !
!     It is relatively common in linear algebra to pass sub-matrices   !
!     as an argument. However, this feature was declared obsolescent   !
!     in Fortran95. The colon we need an explicit interface. Currently,!
!     if fails to compile if colon are used intead of asterisk, even   !
!     if it is inside this module.                                     !
!                                                                      !
!   - July 2026: The LAPACK module was upated to include the DSTEVR    !
!     tridiagonal solver to replace the legacy/archaic DSTEQR. It uses !
!     the Multiple Relatively Robust Representations (MRRR) algorithm  !
!     to compute specific eigenvalues and eigenvectors. It evaluates   !
!     full or partial spectrums extremely fast, completely free of     !
!     obsolete control flows.                                          !
!     Source: https://github.com/Reference-LAPACK/lapack               !
!             98e5d7a23  (HEAD -> master, origin/master, origin/HEAD)  !
!             v3.12.1-351-g98e5d7a23 (Wed Jun 24 02:20:05 2026 -0600)  !
!                                                                      !
!======================================================================!
!
      USE mod_kinds, ONLY : dp
!
      implicit none
!
      PUBLIC  :: DSTEVR             ! modern
      PUBLIC  :: DSTEQR             ! legacy archaic routine
!
      PRIVATE :: DSTEBZ
      PRIVATE :: DSTEIN
      PRIVATE :: DSTEMR
      PRIVATE :: DSTERF
!
      PRIVATE :: DAXPY              ! BLAS
      PRIVATE :: DCOPY              ! BLAS
      PRIVATE :: DDOT               ! BLAS
      PRIVATE :: DNRM2              ! BLAS
      PRIVATE :: IDAMAX             ! BLAS
      PRIVATE :: LSAME              ! BLAS
!
      PRIVATE :: DISNAN             ! LAPACK
      PRIVATE :: DLAEBZ             ! LAPACK
      PRIVATE :: DLAE2              ! LAPACK
      PRIVATE :: DLAEV2             ! LAPACK
      PRIVATE :: DLAGTF             ! LAPACK
      PRIVATE :: DLAGTS             ! LAPACK
      PRIVATE :: DLAISNAN           ! LAPACK
      PRIVATE :: DLAMC1             ! LAPACK
      PRIVATE :: DLAMC2             ! LAPACK
      PRIVATE :: DLAMC3             ! LAPACK
      PRIVATE :: DLAMC4             ! LAPACK
      PRIVATE :: DLAMC5             ! LAPACK
      PRIVATE :: DLAMCH             ! LAPACK
      PRIVATE :: DLANEG             ! LAPACK
      PRIVATE :: DLANST             ! LAPACK
      PRIVATE :: DLAPY2             ! LAPACK
      PRIVATE :: DLAR1V             ! LAPACK
      PRIVATE :: DLARNV             ! LAPACK
      PRIVATE :: DLARRA             ! LAPACK
      PRIVATE :: DLARRB             ! LAPACK
      PRIVATE :: DLARRC             ! LAPACK
      PRIVATE :: DLARRD             ! LAPACK
      PRIVATE :: DLARRE             ! LAPACK
      PRIVATE :: DLARRF             ! LAPACK
      PRIVATE :: DLARRJ             ! LAPACK
      PRIVATE :: DLARRK             ! LAPACK
      PRIVATE :: DLARRR             ! LAPACK
      PRIVATE :: DLARRV             ! LAPACK
      PRIVATE :: DLARTG             ! LAPACK
      PRIVATE :: DLARUV             ! LAPACK
      PRIVATE :: DLASCL             ! LAPACK
      PRIVATE :: DLASET             ! LAPACK
      PRIVATE :: DLASQ2             ! LAPACK
      PRIVATE :: DLASQ3             ! LAPACK
      PRIVATE :: DLASQ4             ! LAPACK
      PRIVATE :: DLASQ5             ! LAPACK
      PRIVATE :: DLASQ6             ! LAPACK
      PRIVATE :: DLASR              ! LAPACK
      PRIVATE :: DLASRT             ! LAPACK
      PRIVATE :: DLASSQ             ! LAPACK
      PRIVATE :: DSCAL              ! LAPACK
      PRIVATE :: DSWAP              ! LAPACK
      PRIVATE :: IEEECK             ! LAPACK
      PRIVATE :: ILAENV             ! LAPACK
      PRIVATE :: IPARMQ             ! LAPACK
      PRIVATE :: XERBLA             ! LAPACK
      PRIVATE
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      CONTAINS
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!     Eigenvalues and Eigenvectors Drivers:
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
      SUBROUTINE DSTEVR (JOBZ, RANGE, N, D, E, VL, VU, IL, IU, ABSTOL,  &
     &                   M, W, Z, LDZ, ISUPPZ, WORK, LWORK, IWORK,      &
     &                   LIWORK, INFO)
!
!=======================================================================
!                                                                      !
!  DSTEVR computes selected eigenvalues and, optionally, eigenvectors  !
!  of a real symmetric tridiagonal matrix T.  Eigenvalues and          !
!  eigenvectors can be selected by specifying either a range of values !
!  or a range of indices for the desired eigenvalues.                  !
!                                                                      !
!  Whenever possible, DSTEVR calls DSTEMR to compute the               !
!  eigenspectrum using Relatively Robust Representations.              !
!  DSTEMR computes eigenvalues by the dqds algorithm, while            !
!  orthogonal eigenvectors are computed from various "good" L D L^T    !
!  representations (also known as Relatively Robust Representations).  !
!  Gram-Schmidt  orthogonalization is avoided as far as possible.      !
!                                                                      !
!  More specifically, the various steps of the algorithm are as        !
!  follows. For the i-th unreduced block of T,                         !
!                                                                      !
!     (a) Compute T - sigma_i = L_i D_i L_i^T, such that L_i D_i L_i^T !
!         is a relatively robust representation,                       !
!     (b) Compute the eigenvalues, lambda_j, of L_i D_i L_i^T to high  !
!         relative accuracy by the dqds algorithm,                     !
!     (c) If there is a cluster of close eigenvalues, "choose" sigma_i !
!         close to the cluster, and go to step (a),                    !
!     (d) Given the approximate eigenvalue lambda_j of L_i D_i L_i^T,  !
!         compute the corresponding eigenvector by forming a           !
!         rank-revealing twisted factorization.                        !
!                                                                      !
!  The desired accuracy of the output can be specified by the input    !
!  parameter ABSTOL.                                                   !
!                                                                      !
!  For more details, see "A new O(n^2) algorithm for the symmetric     !
!  tridiagonal eigenvalue/eigenvector problem", by Inderjit Dhillon,   !
!  Computer Science Division Technical Report No. UCB//CSD-97-971,     !
!  UC Berkeley, May 1997.                                              !
!                                                                      !
!                                                                      !
!  Note 1 : DSTEVR calls DSTEMR when the full spectrum is requested    !
!  on machines which conform to the ieee-754 floating point standard.  !
!  DSTEVR calls DSTEBZ and DSTEIN on non-ieee machines and             !
!  when partial spectrum requests are made.                            !
!                                                                      !
!  Normal execution of DSTEMR may create NaNs and infinities and       !
!  hence may abort due to a floating point exception in environments   !
!  which do not handle NaNs and infinities in the ieee standard        !
!  default manner.                                                     !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  JOBZ     [in]    CHARACTER*1                                        !
!           = 'N':  Compute eigenvalues only;                          !
!           = 'V':  Compute eigenvalues and eigenvectors.              !
!                                                                      !
!  RANGE    [in]    CHARACTER*1                                        !
!           = 'A': all eigenvalues will be found.                      !
!           = 'V': all eigenvalues in the half-open interval (VL,VU]   !
!                  will be found.                                      !
!           = 'I': the IL-th through IU-th eigenvalues will be found.  !
!           For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and     !
!           DSTEIN are called                                          !
!                                                                      !
!  N        [in]    INTEGER                                            !
!           The order of the matrix.  N >= 0.                          !
!                                                                      !
!  D        [in,out] REAL array, dimension (N)                         !
!           On entry, the N diagonal elements of the tridiagonal       !
!           matrix A.                                                  !
!           On exit, D may be multiplied by a constant factor chosen   !
!           to avoid over/underflow in computing the eigenvalues.      !
!                                                                      !
!  E        [in,out] REAL array, dimension (max(1,N-1))                !
!           On entry, the (n-1) subdiagonal elements of the            !
!           tridiagonal matrix A in elements 1 to N-1 of E.            !
!           On exit, E may be multiplied by a constant factor chosen   !
!           to avoid over/underflow in computing the eigenvalues.      !
!                                                                      !
!  VL       [in]    REAL                                               !
!           If RANGE='V', the lower bound of the interval to           !
!           be searched for eigenvalues. VL < VU.                      !
!           Not referenced if RANGE = 'A' or 'I'.                      !
!                                                                      !
!  VU       [in]    REAL                                               !
!           If RANGE='V', the upper bound of the interval to           !
!           be searched for eigenvalues. VL < VU.                      !
!           Not referenced if RANGE = 'A' or 'I'.                      !
!                                                                      !
!  IL       [in]    INTEGER                                            !
!           If RANGE='I', the index of the                             !
!           smallest eigenvalue to be returned.                        !
!           1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.  !
!           Not referenced if RANGE = 'A' or 'V'.                      !
!                                                                      !
!  IU       [in]    INTEGER                                            !
!           If RANGE='I', the index of the                             !
!           largest eigenvalue to be returned.                         !
!           1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.  !
!           Not referenced if RANGE = 'A' or 'V'.                      !
!                                                                      !
!  ABSTOL   [in]    REAL                                               !
!           The absolute error tolerance for the eigenvalues.          !
!           An approximate eigenvalue is accepted as converged         !
!           when it is determined to lie in an interval [a,b]          !
!           of width less than or equal to                             !
!                                                                      !
!                   ABSTOL + EPS *   max(|a|,|b|) ,                    !
!                                                                      !
!           where EPS is the machine precision.  If ABSTOL is less     !
!           than or equal to zero, then  EPS*|T|  will be used in its  !
!           place, where |T| is the 1-norm of the tridiagonal matrix   !
!           obtained by reducing A to tridiagonal form.                !
!                                                                      !
!           See "Computing Small Singular Values of Bidiagonal         !
!           Matrices with Guaranteed High Relative Accuracy," by       !
!           Demmel and Kahan, LAPACK Working Note #3.                  !
!                                                                      !
!           If high relative accuracy is important, set ABSTOL to      !
!           DLAMCH('Safe minimum').  Doing so will guarantee that      !
!           eigenvalues are computed to high relative accuracy when    !
!           possible in future releases.  The current code does not    !
!           make any guarantees about high relative accuracy, but      !
!           future releases will. See J. Barlow and J. Demmel,         !
!           "Computing Accurate Eigensystems of Scaled Diagonally      !
!           Dominant Matrices", LAPACK Working Note #7, for a          !
!           discussion of which matrices define their eigenvalues to   !
!           high relative accuracy.                                    !
!                                                                      !
!  M        [out]   INTEGER                                            !
!           The total number of eigenvalues found.  0 <= M <= N.       !
!           If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.    !
!                                                                      !
!  W        [out]   REAL array, dimension (N)                          !
!           The first M elements contain the selected eigenvalues in   !
!           ascending order.                                           !
!                                                                      !
!  Z        [out]   REAL array, dimension (LDZ, max(1,M))              !
!           If JOBZ = 'V', then if INFO = 0, the first M columns of Z  !
!           contain the orthonormal eigenvectors of the matrix A       !
!           corresponding to the selected eigenvalues, with the i-th   !
!           column of Z holding the eigenvector associated with W(i).  !
!                                                                      !
!           Note: the user must ensure that at least max(1,M) columns  !
!           are supplied in the array Z; if RANGE = 'V', the exact     !
!           value of M is not known in advance and an upper bound must !
!           be used.                                                   !
!                                                                      !
!  LDZ      [in]    INTEGER                                            !
!           The leading dimension of the array Z.  LDZ >= 1, and if    !
!           JOBZ = 'V', LDZ >= max(1,N).                               !
!                                                                      !
!  ISUPPZ   [out]   INTEGER array, dimension (2*max(1,M))              !
!           The support of the eigenvectors in Z, i.e., the indices    !
!           indicating the nonzero elements in Z. The i-th eigenvector !
!           is nonzero only in elements ISUPPZ(2*i-1) through          !
!           ISUPPZ(2*i).                                               !
!           Implemented only for RANGE = 'A' or 'I' and IU-IL = N-1    !
!                                                                      !
!  WORK     [out]   REAL array, dimension (MAX(1,LWORK))               !
!           On exit, if INFO = 0, WORK(1) returns the optimal (and     !
!           minimal) LWORK.                                            !
!                                                                      !
!  LWORK    [in]    INTEGER                                            !
!           The dimension of the array WORK.  LWORK >= max(1,20*N).    !
!                                                                      !
!           If LWORK = -1, then a workspace query is assumed; the      !
!           routine only calculates the optimal sizes of the WORK and  !
!           IWORK arrays, returns these values as the first entries of !
!           the WORK and IWORK arrays, and no error message related to !
!           LWORK or LIWORK is issued by XERBLA.                       !
!                                                                      !
!  IWORK    [out]   INTEGER array, dimension (MAX(1,LIWORK))           !
!           On exit, if INFO = 0, IWORK(1) returns the optimal (and    !
!           minimal) LIWORK.                                           !
!                                                                      !
!  LIWORK   [in]    INTEGER                                            !
!           The dimension of the array IWORK.  LIWORK >= max(1,10*N).  !
!                                                                      !
!           If LIWORK = -1, then a workspace query is assumed; the     !
!           routine only calculates the optimal sizes of the WORK and  !
!           IWORK arrays, returns these values as the first entries of !
!           the WORK and IWORK arrays, and no error message related to !
!           LWORK or LIWORK is issued by XERBLA.                       !
!                                                                      !
!  INFO     [out]   INTEGER                                            !
!           = 0:  successful exit                                      !
!           < 0:  if INFO = -i, the i-th argument had an illegal value !
!           > 0:  Internal error                                       !
!                                                                      !
!  Authors:  LAPACK driver routine                                     !
!  ========                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!  ============                                                        !
!                                                                      !
!   Inderjit Dhillon, IBM Almaden, USA                                 !
!   Osni Marques, LBNL/NERSC, USA                                      !
!   Ken Stanley, Computer Science Division, University of California   !
!                at Berkeley, USA                                      !
!                                                                      !
!  Online html documentation available at                              !
!              http://www.netlib.org/lapack/explore-html/              !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), intent(in) :: JOBZ, RANGE
!
      integer,   intent(in)    :: IL, IU, LDZ, LIWORK, LWORK, N
      integer,   intent(out)   :: INFO, M
      integer,   intent(out)   :: ISUPPZ(*), IWORK(*)
!
      real (dp), intent(in)    :: ABSTOL, VL, VU
      real (dp), intent(inout) :: D(*), E(*)
      real (dp), intent(out)   :: W(*), WORK(*), Z(LDZ,*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
!
      logical ::  ALLEIG, INDEIG, TEST, LQUERY, VALEIG, WANTZ, TRYRAC
!
      integer ::  I, IEEEOK, IMAX, INDIBL, INDIFL, INDISP
      integer ::  INDIWO, ISCALE, ITMP1, J, JJ, LIWMIN, LWMIN, NSPLIT
!
      real (dp) :: BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA, SMLNUM
      real (dp) :: TMP1, TNRM, VLL, VUU
!
      character (len=1) :: ORDER
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input parameters.
!
      IEEEOK = ILAENV(10, 'DSTEVR', 'N', 1, 2, 3, 4)
!
      WANTZ = LSAME(JOBZ, 'V')
      ALLEIG = LSAME(RANGE, 'A')
      VALEIG = LSAME(RANGE, 'V')
      INDEIG = LSAME(RANGE, 'I')
!
      LQUERY = ((LWORK.EQ.-1) .OR. (LIWORK.EQ.-1))
      LWMIN = MAX(1, 20*N)
      LIWMIN = MAX(1, 10*N)
!
      INFO = 0
      IF (.NOT.(WANTZ .OR. LSAME(JOBZ, 'N'))) THEN
         INFO = -1
      ELSE IF (.NOT.(ALLEIG .OR. VALEIG .OR. INDEIG)) THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      ELSE
         IF (VALEIG) THEN
            IF (N.GT.0 .AND. VU.LE.VL) INFO = -7
         ELSE IF (INDEIG) THEN
            IF (IL.LT.1 .OR. IL.GT.MAX(1, N)) THEN
               INFO = -8
            ELSE IF (IU.LT.MIN(N, IL) .OR. IU.GT.N) THEN
               INFO = -9
            END IF
         END IF
      END IF
      IF (INFO.EQ.0) THEN
         IF (LDZ.LT.1 .OR. (WANTZ .AND. LDZ.LT.N)) THEN
            INFO = -14
         END IF
      END IF
!
      IF (INFO.EQ.0) THEN
         WORK(1) = LWMIN
         IWORK(1) = LIWMIN
!
         IF (LWORK.LT.LWMIN .AND. .NOT.LQUERY) THEN
            INFO = -17
         ELSE IF (LIWORK.LT.LIWMIN .AND. .NOT.LQUERY) THEN
            INFO = -19
         END IF
      END IF
!
      IF (INFO.NE.0) THEN
         CALL XERBLA ('DSTEVR', -INFO)
         RETURN
      ELSE IF (LQUERY) THEN
         RETURN
      END IF
!
!  Quick return if possible
!
      M = 0
      IF (N.EQ.0) RETURN
!
      IF (N.EQ.1) THEN
         IF (ALLEIG .OR. INDEIG) THEN
            M = 1
            W(1) = D(1)
         ELSE
            IF (VL.LT.D(1) .AND. VU.GE.D(1)) THEN
               M = 1
               W(1) = D(1)
            END IF
         END IF
         IF (WANTZ) THEN
            Z(1, 1) = ONE
            RETURN
         END IF
      END IF
!
!  Get machine constants.
!
      SAFMIN = DLAMCH('Safe minimum')
      EPS = DLAMCH('Precision')
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT(SMLNUM)
      RMAX = MIN(SQRT(BIGNUM), ONE / SQRT(SQRT(SAFMIN)))
!
!  Scale matrix to allowable range, if necessary.
!
      ISCALE = 0
      IF (VALEIG) THEN
         VLL = VL
         VUU = VU
      END IF
!
      TNRM = DLANST('M', N, D, E)
      IF (TNRM.GT.ZERO .AND. TNRM.LT.RMIN) THEN
         ISCALE = 1
         SIGMA = RMIN / TNRM
      ELSE IF (TNRM.GT.RMAX) THEN
         ISCALE = 1
         SIGMA = RMAX / TNRM
      END IF
      IF (ISCALE.EQ.1) THEN
         CALL DSCAL (N, SIGMA, D, 1)
         CALL DSCAL (N-1, SIGMA, E(1), 1)
         IF (VALEIG) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF
!
!  Initialize indices into workspaces.  Note: These indices are used only
!  if DSTERF or DSTEMR fail.
!
!  IWORK(INDIBL:INDIBL+M-1) corresponds to IBLOCK in DSTEBZ and
!  stores the block indices of each of the M<=N eigenvalues.
!
      INDIBL = 1
!
!  IWORK(INDISP:INDISP+NSPLIT-1) corresponds to ISPLIT in DSTEBZ and
!  stores the starting and finishing indices of each block.
!
      INDISP = INDIBL + N
!
!  IWORK(INDIFL:INDIFL+N-1) stores the indices of eigenvectors
!  that corresponding to eigenvectors that fail to converge in
!  DSTEIN.  This information is discarded; if any fail, the driver
!  returns INFO > 0.
!
      INDIFL = INDISP + N
!
!  INDIWO is the offset of the remaining integer workspace.
!
      INDIWO = INDISP + N
!
!  If all eigenvalues are desired, then
!  call DSTERF or DSTEMR.  If this fails for some eigenvalue, then
!  try DSTEBZ.
!
      TEST = .FALSE.
      IF (INDEIG) THEN
         IF (IL.EQ.1 .AND. IU.EQ.N) THEN
            TEST = .TRUE.
         END IF
      END IF
      IF ((ALLEIG .OR. TEST) .AND. IEEEOK.EQ.1) THEN
         CALL DCOPY (N-1, E(1), 1, WORK(1), 1)
         IF (.NOT.WANTZ) THEN
            CALL DCOPY (N, D, 1, W, 1)
            CALL DSTERF (N, W, WORK, INFO)
         ELSE
            CALL DCOPY (N, D, 1, WORK(N+1), 1)
            IF (ABSTOL .LE. TWO*N*EPS) THEN
               TRYRAC = .TRUE.
            ELSE
               TRYRAC = .FALSE.
            END IF
            CALL DSTEMR (JOBZ, 'A', N, WORK(N+1), WORK, VL, VU, IL,     &
     &                   IU, M, W, Z, LDZ, N, ISUPPZ, TRYRAC,           &
     &                   WORK(2*N+1), LWORK-2*N, IWORK, LIWORK, INFO)
         END IF
         IF (INFO.EQ.0) THEN
            M = N
            GOTO 10
         END IF
         INFO = 0
      END IF
!
!  Otherwise, call DSTEBZ and, if eigenvectors are desired, DSTEIN.
!
      IF (WANTZ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
!
      CALL DSTEBZ (RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTOL, D, E,     &
     &             M, NSPLIT, W, IWORK(INDIBL), IWORK(INDISP), WORK,    &
     &             IWORK(INDIWO), INFO)
!
      IF (WANTZ) THEN
         CALL DSTEIN (N, D, E, M, W, IWORK(INDIBL), IWORK(INDISP),      &
     &                Z, LDZ, WORK, IWORK(INDIWO), IWORK(INDIFL),       &
     &                INFO)
      END IF
!
!  If matrix was scaled, then rescale eigenvalues appropriately.
!
   10 CONTINUE
      IF (ISCALE.EQ.1) THEN
         IF (INFO.EQ.0) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL (IMAX, ONE / SIGMA, W, 1)
      END IF
!
!  If eigenvalues are not in order, then sort them, along with
!  eigenvectors.
!
      IF (WANTZ) THEN
         DO J = 1, M - 1
            I = 0
            TMP1 = W(J)
            DO JJ = J + 1, M
               IF (W(JJ).LT.TMP1) THEN
                  I = JJ
                  TMP1 = W(JJ)
               END IF
            END DO
!
            IF (I.NE.0) THEN
               ITMP1 = IWORK(I)
               W(I) = W(J)
               IWORK(I) = IWORK(J)
               W(J) = TMP1
               IWORK(J) = ITMP1
               CALL DSWAP (N, Z(1,I), 1, Z(1,J), 1)
            END IF
         END DO
      END IF
!
!  Causes problems with tests 19 & 20:
!  IF (wantz .and. INDEIG) Z(1,1) = Z(1,1) / 1.002 + .002
!
      WORK(1) = LWMIN
      IWORK(1) = LIWMIN
!
      RETURN
      END SUBROUTINE DSTEVR
!
      SUBROUTINE DSTEQR (COMPZ, N, D, E, Z, LDZ, WORK, INFO)
!
!=======================================================================
!                                                                      !
!  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a  !
!  symmetric tridiagonal matrix using the implicit QL or QR method.    !
!  The eigenvectors of a full or band symmetric matrix can also be     !
!  found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this    !
!  matrix to tridiagonal form.                                         !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  COMPZ   (input) CHARACTER*1                                         !
!          = 'N':  Compute eigenvalues only.                           !
!          = 'V':  Compute eigenvalues and eigenvectors of the         !
!                  original symmetric matrix. On entry, Z must contain !
!                  the orthogonal matrix used to reduce the original   !
!                  matrix to tridiagonal form.                         !
!          = 'I':  Compute eigenvalues and eigenvectors of the         !
!                  tridiagonal matrix.  Z is initialized to the        !
!                  identity matrix.                                    !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The order of the matrix.  N >= 0.                           !
!                                                                      !
!  D       (input/output) DOUBLE PRECISION array, dimension (N)        !
!          On entry, the diagonal elements of the tridiagonal matrix.  !
!          On exit, if INFO = 0, the eigenvalues in ascending order.   !
!                                                                      !
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)      !
!          On entry, the (n-1) subdiagonal elements of the tridiagonal !
!          matrix.                                                     !
!          On exit, E has been destroyed.                              !
!                                                                      !
!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)   !
!          On entry, if  COMPZ = 'V', then Z contains the orthogonal   !
!          matrix used in the reduction to tridiagonal form.           !
!          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the  !
!          orthonormal eigenvectors of the original symmetric matrix,  !
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors !
!          of the symmetric tridiagonal matrix.                        !
!          If COMPZ = 'N', then Z is not referenced.                   !
!                                                                      !
!  LDZ     (input) INTEGER                                             !
!          The leading dimension of the array Z.  LDZ >= 1, and if     !
!          eigenvectors are desired, then  LDZ >= max(1,N).            !
!                                                                      !
!  WORK    (workspace) DOUBLE PRECISION array, dimension               !
!          (max(1,2*N-2)).                                             !
!          If COMPZ = 'N', then WORK is not referenced.                !
!                                                                      !
!  INFO    (output) INTEGER                                            !
!          = 0:  successful exit                                       !
!          < 0:  if INFO = -i, the i-th argument had an illegal value  !
!          > 0:  the algorithm has failed to find all the eigenvalues  !
!                in a total of 30*N iterations; if INFO = i, then i    !
!                elements of E have not converged to zero; on exit, D  !
!                and E contain the elements of a symmetric tridiagonal !
!                matrix which is orthogonally similar to the original  !
!                matrix.                                               !
!                                                                      !
!  -- LAPACK routine (version 2.0) --                                  !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     September 30, 1994                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: LDZ, N
      integer,   intent(out)   :: INFO
!
      real (dp), intent(inout) :: WORK(*)
      real (dp), intent(inout) :: D(*), E(*), Z(LDZ,*)
!
      character (len=1), intent(in) :: COMPZ
!
!  Local variable declarations.
!
      integer, parameter :: MAXIT = 30
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: THREE = 3.0_dp
!
      integer :: I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND, LENDM1
      integer :: LENDP1, LENDSV, LM1, LSV, M, MM, MM1, NM1, NMAXIT
!
      real (dp) :: ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2
      real (dp) :: S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input parameters.
!
      INFO = 0
!
      IF (LSAME(COMPZ, 'N')) THEN
        ICOMPZ = 0
      ELSE IF (LSAME(COMPZ, 'V')) THEN
        ICOMPZ = 1
      ELSE IF (LSAME(COMPZ, 'I')) THEN
        ICOMPZ = 2
      ELSE
        ICOMPZ = -1
      END IF
!
      IF (ICOMPZ.LT.0) THEN
        INFO = -1
      ELSE IF (N.LT.0) THEN
        INFO = -2
      ELSE IF ((LDZ.LT.1).OR.(ICOMPZ.GT.0 .AND. LDZ.LT.MAX(1,N))) THEN
        INFO = -6
      END IF
!
      IF (INFO.NE.0) THEN
        CALL XERBLA ('DSTEQR', -INFO)
        RETURN
      END IF
!
!  Quick return if possible
!
      IF (N.EQ.0) RETURN
!
      IF (N.EQ.1) THEN
        IF (ICOMPZ.EQ.2) Z(1,1) = ONE
        RETURN
      END IF
!
!  Determine the unit roundoff and over/underflow thresholds.
!
      EPS = DLAMCH('E')
      EPS2 = EPS**2
      SAFMIN = DLAMCH('S')
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT(SAFMAX) / THREE
      SSFMIN = SQRT(SAFMIN) / EPS2
!
!  Compute the eigenvalues and eigenvectors of the tridiagonal
!  matrix.
!
      IF (ICOMPZ .EQ. 2) THEN
        CALL DLASET ('Full', N, N, ZERO, ONE, Z, LDZ)
      END IF
!
      NMAXIT = N*MAXIT
      JTOT = 0
!
!  Determine where the matrix splits and choose QL or QR iteration
!  for each block, according to whether top or bottom diagonal
!  element is smaller.
!
      L1 = 1
      NM1 = N - 1
!
   10 CONTINUE
      IF (L1 .GT. N) GOTO 160
      IF (L1 .GT. 1) E(L1-1) = ZERO
      IF (L1 .LE. NM1) THEN
        DO M = L1, NM1
          TST = ABS(E(M))
          IF (TST .EQ. ZERO) GOTO 30
          IF (TST .LE. (SQRT(ABS(D(M)))*SQRT(ABS(D(M+1))))*EPS) THEN
            E(M) = ZERO
            GOTO 30
          END IF
        END DO
      END IF
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF (LEND .EQ. L) GOTO 10
!
!  Scale submatrix in rows and columns L to LEND.
!
      ANORM = DLANST('I', LEND-L+1, D(L), E(L))
      ISCALE = 0
      IF (ANORM .EQ. ZERO) GOTO 10
      IF (ANORM .GT. SSFMAX) THEN
        ISCALE = 1
        CALL DLASCL ('G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D(L), N,    &
     &               INFO)
        CALL DLASCL ('G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E(L), N,      &
     &               INFO)
      ELSE IF (ANORM .LT. SSFMIN) THEN
        ISCALE = 2
        CALL DLASCL ('G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D(L), N,    &
     &               INFO)
        CALL DLASCL ('G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E(L), N,      &
     &               INFO)
      END IF
!
!  Choose between QL and QR iteration
!
      IF (ABS(D(LEND)) .LT. ABS(D(L))) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF (LEND .GT. L) THEN
!
!  QL Iteration.
!
!  Look for small subdiagonal element.
!
   40   CONTINUE
        IF (L .NE. LEND) THEN
          LENDM1 = LEND - 1
          DO M = L, LENDM1
            TST = ABS(E(M))**2
            IF (TST .LE. (EPS2*ABS(D(M)))*ABS(D(M+1))+SAFMIN) GOTO 60
          END DO
        END IF
!
        M = LEND
!
   60   CONTINUE
        IF (M .LT. LEND) E(M) = ZERO
        P = D(L)
        IF (M .EQ. L) GOTO 80
!
!  If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!  to compute its eigensystem.
!
        IF (M .EQ. L+1) THEN
          IF (ICOMPZ .GT. 0) THEN
            CALL DLAEV2 (D(L), E(L), D(L+1), RT1, RT2, C, S)
            WORK(L) = C
            WORK(N-1+L) = S
            CALL DLASR ('R', 'V', 'B', N, 2, WORK(L),                   &
     &                  WORK(N-1+L), Z(1, L), LDZ)
          ELSE
            CALL DLAE2 (D(L), E(L), D(L+1), RT1, RT2)
          END IF
          D(L) = RT1
          D(L+1) = RT2
          E(L) = ZERO
          L = L + 2
          IF (L .LE. LEND) GOTO 40
          GOTO 140
        END IF
!
        IF (JTOT .EQ. NMAXIT) GOTO 140
        JTOT = JTOT + 1
!
!  Form shift.
!
        G = (D(L+1)-P) / (TWO*E(L))
        R = DLAPY2(G, ONE)
        G = D(M) - P + (E(L) / (G + SIGN(R, G)))
!
        S = ONE
        C = ONE
        P = ZERO
!
!  Inner loop
!
        MM1 = M - 1
        DO I = MM1, L, -1
          F = S*E(I)
          B = C*E(I)
          CALL DLARTG (G, F, C, S, R)
          IF (I .NE. M-1) E(I+1) = R
          G = D(I+1) - P
          R = (D(I)-G)*S + TWO*C*B
          P = S*R
          D(I+1) = G+P
          G = C*R - B
!
!  If eigenvectors are desired, then save rotations.
!
          IF (ICOMPZ.GT.0) THEN
            WORK(I) = C
            WORK(N-1+I) = -S
          END IF
        END DO
!
!  If eigenvectors are desired, then apply saved rotations.
!
        IF (ICOMPZ.GT.0) THEN
          MM = M - L + 1
          CALL DLASR ('R', 'V', 'B', N, MM, WORK(L), WORK(N-1+L),       &
     &                Z(1,L), LDZ)
        END IF
!
        D(L) = D(L) - P
        E(L) = G
        GOTO 40
!
!  Eigenvalue found.
!
   80   CONTINUE
        D(L) = P
 !
        L = L + 1
        IF (L .LE. LEND) GOTO 40
        GOTO 140
!
      ELSE
!
!  QR Iteration
!
!  Look for small superdiagonal element.
!
   90   CONTINUE
        IF (L .NE. LEND) THEN
          LENDP1 = LEND + 1
          DO M = L, LENDP1, -1
            TST = ABS(E(M-1))**2
            IF (TST .LE. (EPS2*ABS(D(M)))*ABS(D(M-1))+SAFMIN) GOTO 110
          END DO
        END IF
!
        M = LEND
!
  110   CONTINUE
        IF (M .GT. LEND) E(M-1) = ZERO
        P = D(L)
        IF (M .EQ. L) GOTO 130
!
!  If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!  to compute its eigensystem.
!
        IF (M .EQ. L-1) THEN
          IF (ICOMPZ .GT. 0) THEN
            CALL DLAEV2 (D(L-1), E(L-1), D(L), RT1, RT2, C, S)
            WORK(M) = C
            WORK(N-1+M) = S
            CALL DLASR ('R', 'V', 'F', N, 2, WORK(M),                   &
     &                  WORK(N-1+M), Z(1,L-1), LDZ)
          ELSE
            CALL DLAE2 (D(L-1), E(L-1), D(L), RT1, RT2)
          END IF
          D(L-1) = RT1
          D(L) = RT2
          E(L-1) = ZERO
          L = L - 2
          IF (L .GE. LEND) GOTO 90
          GOTO 140
        END IF
!
        IF (JTOT .EQ. NMAXIT) GOTO 140
        JTOT = JTOT + 1
!
!  Form shift.
!
        G = (D(L-1)-P) / (TWO*E(L-1))
        R = DLAPY2(G, ONE)
        G = D(M) - P + (E(L-1) / (G+SIGN(R, G)))
!
        S = ONE
        C = ONE
        P = ZERO
!
!  Inner loop
!
        LM1 = L - 1
        DO I = M, LM1
          F = S*E(I)
          B = C*E(I)
          CALL DLARTG (G, F, C, S, R)
          IF (I .NE. M) E(I-1) = R
          G = D(I) - P
          R = (D(I+1)-G)*S + TWO*C*B
          P = S*R
          D(I) = G + P
          G = C*R - B
!
!  If eigenvectors are desired, then save rotations.
!
          IF (ICOMPZ.GT.0) THEN
            WORK(I) = C
            WORK(N-1+I) = S
          END IF
        END DO
!
!  If eigenvectors are desired, then apply saved rotations.
!
        IF (ICOMPZ.GT.0) THEN
          MM = L - M + 1
          CALL DLASR ('R', 'V', 'F', N, MM, WORK(M), WORK(N-1+M),       &
     &                Z(1,M), LDZ)
        END IF
!
        D(L) = D(L) - P
        E(LM1) = G
        GOTO 90
!
!  Eigenvalue found.
!
  130   CONTINUE
        D(L) = P
!
        L = L - 1
        IF (L .GE. LEND) GOTO 90
        GOTO 140
!
      END IF
!
!  Undo scaling if necessary.
!
  140 CONTINUE
      IF (ISCALE .EQ. 1) THEN
        CALL DLASCL ('G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,         &
     &               D(LSV), N, INFO)
        CALL DLASCL ('G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E(LSV),   &
     &               N, INFO)
      ELSE IF (ISCALE .EQ. 2) THEN
        CALL DLASCL ('G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,         &
     &               D(LSV), N, INFO)
        CALL DLASCL ('G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E(LSV),   &
     &               N, INFO)
      END IF
!
!  Check for no convergence to an eigenvalue after a total
!  of N*MAXIT iterations.
!
      IF (JTOT .LT. NMAXIT) GOTO 10
      DO I = 1, N - 1
        IF (E(I) .NE. ZERO) INFO = INFO + 1
      END DO
      GOTO 190
!
!  Order eigenvalues and eigenvectors.
!
  160 CONTINUE
      IF (ICOMPZ .EQ. 0) THEN
!
!    Use Quick Sort.
!
        CALL DLASRT ('I', N, D, INFO)
!
      ELSE
!
!    Use Selection Sort to minimize swaps of eigenvectors.
!
        DO II = 2, N
          I = II - 1
          K = I
          P = D(I)
          DO  J = II, N
            IF (D(J) .LT. P) THEN
               K = J
               P = D(J)
            END IF
          END DO
          IF (K .NE. I) THEN
            D(K) = D(I)
            D(I) = P
            CALL DSWAP (N, Z(1,I), 1, Z(1,K), 1)
          END IF
        END DO
      END IF
!
  190 CONTINUE
      RETURN
!
      END SUBROUTINE DSTEQR
!
      SUBROUTINE DSTEBZ (RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E, &
     &                   M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,     &
     &                   INFO)
!
!=======================================================================
!                                                                      !
!  DSTEBZ computes the eigenvalues of a symmetric tridiagonal          !
!  matrix T.  The user may ask for all eigenvalues, all eigenvalues    !
!  in the half-open interval (VL, VU], or the IL-th through IU-th      !
!  eigenvalues.                                                        !
!                                                                      !
!  To avoid overflow, the matrix must be scaled so that its largest    !
!  element is no greater than overflow**(1/2) * underflow**(1/4) in    !
!  absolute value, and for greatest accuracy, it should not be much    !
!  smaller than that.                                                  !
!                                                                      !
!  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal       !
!  Matrix", Report CS41, Computer Science Dept., Stanford              !
!  University, July 21, 1966.                                          !
!                                                                      !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  RANGE    [in]      CHARACTER*1                                      !
!           = 'A': ("All")   all eigenvalues will be found.            !
!           = 'V': ("Value") all eigenvalues in the half-open interval !
!                            (VL, VU] will be found.                   !
!           = 'I': ("Index") the IL-th through IU-th eigenvalues (of   !
!                            the entire matrix) will be found.         !
!                                                                      !
!  ORDER    [in]      CHARACTER*1                                      !
!           = 'B': ("By Block") the eigenvalues will be grouped by     !
!                               split-off block (see IBLOCK, ISPLIT)   !
!                               and ordered from smallest to largest   !
!                               within the block.                      !
!           = 'E': ("Entire matrix")vthe eigenvalues for the entire    !
!                               matrix will be ordered from smallest   !
!                               to largest.                            !
!  N        [in]      INTEGER                                          !
!           The order of the tridiagonal matrix T.  N >= 0.            !
!                                                                      !
!  VL       [in]      REAL                                             !
!           If RANGE='V', the lower bound of the interval to be        !
!           searched for eigenvalues.  Eigenvalues less than or equal  !
!           to VL, or greater than VU, will not be returned.  VL < VU. !
!           Not referenced if RANGE = 'A' or 'I'.                      !
!                                                                      !
!  VU       [in]      REAL                                             !
!           If RANGE='V', the upper bound of the interval to be        !
!           searched for eigenvalues.  Eigenvalues less than or equal  !
!           to VL, or greater than VU, will not be returned.  VL < VU. !
!           Not referenced if RANGE = 'A' or 'I'.                      !
!                                                                      !
!  IL       [in]      INTEGER                                          !
!           If RANGE='I', the index of the                             !
!           smallest eigenvalue to be returned.                        !
!           1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.  !
!           Not referenced if RANGE = 'A' or 'V'.                      !
!                                                                      !
!  IU       [in]      INTEGER                                          !
!           If RANGE='I', the index of the                             !
!           largest eigenvalue to be returned.                         !
!           1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.  !
!           Not referenced if RANGE = 'A' or 'V'.                      !
!                                                                      !
!  ABSTOL   [in]      REAL                                             !
!           The absolute tolerance for the eigenvalues. An eigenvalue  !
!           (or cluster) is considered to be located if it has been    !
!           determined to lie in an interval whose width is ABSTOL or  !
!           less.  If ABSTOL is less than or equal to zero, then       !
!           ULP*|T| will be used, where |T| means the 1-norm of T.     !
!                                                                      !
!           Eigenvalues will be computed most accurately when ABSTOL   !
!           is set to twice the underflow threshold 2*DLAMCH('S'),     !
!           not zero.                                                  !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           The n diagonal elements of the tridiagonal matrix T.       !
!                                                                      !
!  E        [in]      REAL array, dimension (N-1)                      !
!           The (n-1) off-diagonal elements of the tridiagonal         !
!           matrix T.                                                  !
!                                                                      !
!  M        [out]     INTEGER                                          !
!           The actual number of eigenvalues found. 0 <= M <= N.       !
!           (See also the description of INFO=2,3.)                    !
!                                                                      !
!  NSPLIT   [out]     INTEGER                                          !
!           The number of diagonal blocks in the matrix T.             !
!           1 <= NSPLIT <= N.                                          !
!                                                                      !
!  W        [out]     REAL array, dimension (N)                        !
!           On exit, the first M elements of W will contain the        !
!           eigenvalues.  (DSTEBZ may use the remaining N-M elements   !
!           as workspace.)                                             !
!                                                                      !
!  IBLOCK   [out]     INTEGER array, dimension (N)                     !
!           At each row/column j where E(j) is zero or small, the      !
!           matrix T is considered to split into a block diagonal      !
!           matrix.                                                    !
!           On exit, if INFO = 0, IBLOCK(i) specifies to which         !
!           block (from 1 to the number of blocks) the eigenvalue W(i) !
!           belongs.  (DSTEBZ may use the remaining N-M elements as    !
!           workspace.)                                                !
!                                                                      !
!  ISPLIT   [out]     INTEGER array, dimension (N)                     !
!           The splitting points, at which T breaks up into            !
!           submatrices.  The first submatrix consists of rows/columns !
!           1 to ISPLIT(1), the second of rows/columns ISPLIT(1)+1     !
!           through ISPLIT(2), etc., and the NSPLIT-th consists of     !
!           rows/columns ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.  !
!           (Only the first NSPLIT elements will actually be used, but !
!           since the user cannot know a priori what value NSPLIT will !
!           have, N words must be reserved for ISPLIT.)                !
!                                                                      !
!  WORK     [out]     REAL array, dimension (4*N)                      !
!                                                                      !
!  IWORK    [out]     INTEGER array, dimension (3*N)                   !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           = 0:  successful exit                                      !
!           < 0:  if INFO = -i, the i-th argument had an illegal value !
!           > 0:  some or all of the eigenvalues failed to converge or !
!                 were not computed:                                   !
!                 =1 or 3: Bisection failed to converge for some       !
!                          eigenvalues; these eigenvalues are flagged  !
!                          by a negative block number.  The effect is  !
!                          that the eigenvalues may not be as accurate !
!                          as the absolute and relative tolerances.    !
!                          This is generally caused by unexpectedly    !
!                          inaccurate arithmetic.                      !
!                 =2 or 3: RANGE='I' only: Not all of the eigenvalues  !
!                          IL:IU were found.                           !
!                          Effect: M < IU+1-IL                         !
!                          Cause:  non-monotonic arithmetic, causing   !
!                                  the Sturm sequence to be            !
!                                  non-monotonic.                      !
!                          Cure:   recalculate, using RANGE='A', and   !
!                                  pick out eigenvalues IL:IU. In some !
!                                  cases, increasing the PARAMETER     !
!                                  "FUDGE" may make things work.       !
!                 = 4:     RANGE='I', and the Gershgorin interval      !
!                          initially used was too small. No            !
!                          eigenvalues were computed.                  !
!                          Probable cause: your machine has sloppy     !
!                                          floating-point arithmetic.  !
!                          Cure: Increase the PARAMETER "FUDGE",       !
!                                recompile, and try again.             !
!                                                                      !
!  Internal Parameters:                                                !
!                                                                      !
!   RELFAC  REAL, default = 2.0                                        !
!           The relative tolerance.  An interval (a,b] lies within     !
!           "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),    !
!           where "ulp" is the machine precision (distance from 1 to   !
!           the next larger floating point number.)                    !
!                                                                      !
!   FUDGE   REAL, default = 2.0                                        !
!           A "fudge factor" to widen the Gershgorin intervals.        !
!           Ideally, a value of 1 should work, but on machines with    !
!           sloppy arithmetic, this needs to be larger.  The default   !
!           for publicly released versions should be large enough to   !
!           handle the worst machine around.  Note that this has no    !
!           effect on accuracy of the solution.                        !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), intent(in) :: ORDER, RANGE
!
      integer,   intent(in)  :: IL, IU, N
      integer,   intent(out) :: INFO, M, NSPLIT
      integer,   intent(out) :: IBLOCK(*), ISPLIT(*), IWORK(*)
!
      real (dp), intent(in)  :: ABSTOL, VL, VU
      real (dp), intent(in)  :: D(*), E(*)
      real (dp), intent(out) :: W(*), WORK(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: HALF = 0.5_dp
!
      real (dp), parameter :: FUDGE = 2.1_dp
      real (dp), parameter :: RELFAC = 2.0_dp
!
      logical :: NCNVRG, TOOFEW
!
      integer :: IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO
      integer :: IM, IN, IOFF, IORDER, IOUT, IRANGE, ITMAX
      integer :: ITMP1, IW, IWOFF, J, JB, JDISC, JE, NB, NWL, NWU
!
      integer, dimension(1) :: IDUMMA
!
      real (dp) :: ATOLI, BNORM, GL, GU, PIVMIN, RTOLI, SAFEMN
      real (dp) :: TMP1, TMP2, TNORM, ULP, WKILL, WL, WLU, WU, WUL
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
!
!  Decode RANGE.
!
      IF (LSAME(RANGE, 'A')) THEN
         IRANGE = 1
      ELSE IF (LSAME(RANGE, 'V')) THEN
         IRANGE = 2
      ELSE IF (LSAME(RANGE, 'I')) THEN
         IRANGE = 3
      ELSE
         IRANGE = 0
      END IF
!
!  Decode ORDER
!
      IF (LSAME(ORDER, 'B')) THEN
         IORDER = 2
      ELSE IF (LSAME(ORDER, 'E')) THEN
         IORDER = 1
      ELSE
         IORDER = 0
      END IF
!
!  Check for Errors.
!
      IF (IRANGE.LE.0) THEN
         INFO = -1
      ELSE IF (IORDER.LE.0) THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      ELSE IF (IRANGE.EQ.2) THEN
         IF (VL.GE.VU)  INFO = -5
      ELSE IF (IRANGE.EQ.3 .AND. (IL.LT.1 .OR. IL.GT.MAX(1,N))) THEN
         INFO = -6
      ELSE IF (IRANGE.EQ.3 .AND. (IU.LT.MIN(N,IL) .OR. IU.GT.N)) THEN
         INFO = -7
      END IF
!
      IF (INFO.NE.0) THEN
         CALL XERBLA ('DSTEBZ', -INFO)
         RETURN
      END IF
!
!  Initialize error flags.
!
      INFO = 0
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.
!
!  Quick return if possible.
!
      M = 0
      IF (N.EQ.0) RETURN
!
!  Simplifications.
!
      IF (IRANGE.EQ.3 .AND. IL.EQ.1 .AND. IU.EQ.N) IRANGE = 1
!
!  Get machine constants
!  NB is the minimum vector length for vector bisection, or 0
!  if only scalar is to be done.
!
      SAFEMN = DLAMCH('S')
      ULP = DLAMCH('P')
      RTOLI = ULP*RELFAC
      NB = ILAENV(1, 'DSTEBZ', ' ', N, -1, -1, -1)

      IF (NB.LE.1) NB = 0
!
!  Special Case when N=1.
!
      IF (N.EQ.1) THEN
         NSPLIT = 1
         ISPLIT(1) = 1
         IF (IRANGE.EQ.2 .AND. (VL.GE.D(1) .OR. VU.LT.D(1))) THEN
            M = 0
         ELSE
            W(1) = D(1)
            IBLOCK(1) = 1
            M = 1
         END IF
         RETURN
      END IF
!
!  Compute Splitting Points.
!
      NSPLIT = 1
      WORK(N) = ZERO
      PIVMIN = ONE
!
      DO 10 J = 2, N
         TMP1 = E(J-1)**2
         IF (ABS(D(J)*D(J-1))*ULP**2+SAFEMN.GT.TMP1) THEN
            ISPLIT(NSPLIT) = J - 1
            NSPLIT = NSPLIT + 1
            WORK(J-1) = ZERO
         ELSE
            WORK(J-1) = TMP1
            PIVMIN = MAX(PIVMIN, TMP1)
         END IF
   10 CONTINUE
      ISPLIT(NSPLIT) = N
      PIVMIN = PIVMIN*SAFEMN
!
!  Compute Interval and ATOLI.
!
      IF (IRANGE.EQ.3) THEN
!
!  RANGE='I': Compute the interval containing eigenvalues
!             IL through IU.
!
!  Compute Gershgorin interval for entire (split) matrix
!  and use it as the initial interval.
!
         GU = D(1)
         GL = D(1)
         TMP1 = ZERO
!
         DO 20 J = 1, N - 1
            TMP2 = SQRT(WORK(J))
            GU = MAX(GU, D(J)+TMP1+TMP2)
            GL = MIN(GL, D(J)-TMP1-TMP2)
            TMP1 = TMP2
   20    CONTINUE
!
         GU = MAX(GU, D(N)+TMP1)
         GL = MIN(GL, D(N)-TMP1)
         TNORM = MAX(ABS(GL), ABS(GU))
         GL = GL - FUDGE*TNORM*ULP*N - FUDGE*TWO*PIVMIN
         GU = GU + FUDGE*TNORM*ULP*N + FUDGE*PIVMIN
!
!  Compute Iteration parameters.
!
         ITMAX = INT((LOG(TNORM+PIVMIN)-LOG(PIVMIN)) /                  &
     &           LOG(TWO)) + 2
         IF (ABSTOL.LE.ZERO) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
!
         WORK(N+1) = GL
         WORK(N+2) = GL
         WORK(N+3) = GU
         WORK(N+4) = GU
         WORK(N+5) = GL
         WORK(N+6) = GU
         IWORK(1) = -1
         IWORK(2) = -1
         IWORK(3) = N + 1
         IWORK(4) = N + 1
         IWORK(5) = IL - 1
         IWORK(6) = IU
!
         CALL DLAEBZ (3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D,   &
     &                E, WORK, IWORK(5), WORK(N+1), WORK(N+5), IOUT,    &
     &                IWORK, W, IBLOCK, IINFO)
!
         IF (IWORK(6).EQ.IU) THEN
            WL = WORK(N+1)
            WLU = WORK(N+3)
            NWL = IWORK(1)
            WU = WORK(N+4)
            WUL = WORK(N+2)
            NWU = IWORK(4)
         ELSE
            WL = WORK(N+2)
            WLU = WORK(N+4)
            NWL = IWORK(2)
            WU = WORK(N+3)
            WUL = WORK(N+1)
            NWU = IWORK(3)
         END IF
!
         IF (NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N) THEN
            INFO = 4
            RETURN
         END IF
      ELSE
!
!  RANGE='A' or 'V' -- Set ATOLI
!
         TNORM = MAX(ABS(D(1))+ABS(E(1)), ABS(D(N))+ABS(E(N-1)))
!
         DO 30 J = 2, N - 1
            TNORM = MAX(TNORM, ABS(D(J))+ABS(E(J-1))+ABS(E(J)))
   30    CONTINUE
!
         IF (ABSTOL.LE.ZERO) THEN
            ATOLI = ULP*TNORM
         ELSE
            ATOLI = ABSTOL
         END IF
!
         IF (IRANGE.EQ.2) THEN
            WL = VL
            WU = VU
         ELSE
            WL = ZERO
            WU = ZERO
         END IF
      END IF
!
!  Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
!  NWL accumulates the number of eigenvalues .le. WL,
!  NWU accumulates the number of eigenvalues .le. WU
!
      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0
!
      DO 70 JB = 1, NSPLIT
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT(JB)
         IN = IEND - IOFF
!
         IF (IN.EQ.1) THEN
!
!  Special Case -- IN=1.
!
            IF (IRANGE.EQ.1 .OR. WL.GE.D(IBEGIN)-PIVMIN) NWL = NWL + 1
            IF (IRANGE.EQ.1 .OR. WU.GE.D(IBEGIN)-PIVMIN) NWU = NWU + 1

            IF (IRANGE.EQ.1 .OR. (WL.LT.D(IBEGIN)-PIVMIN .AND.          &
     &                            WU.GE.D(IBEGIN)-PIVMIN)) THEN
               M = M + 1
               W(M) = D(IBEGIN)
               IBLOCK(M) = JB
            END IF
         ELSE
!
!  General Case -- IN > 1
!
!  Compute Gershgorin Interval
!  and use it as the initial interval
!
            GU = D(IBEGIN)
            GL = D(IBEGIN)
            TMP1 = ZERO
!
            DO 40 J = IBEGIN, IEND - 1
               TMP2 = ABS(E(J))
               GU = MAX(GU, D(J)+TMP1+TMP2)
               GL = MIN(GL, D(J)-TMP1-TMP2)
               TMP1 = TMP2
   40       CONTINUE
!
            GU = MAX(GU, D(IEND)+TMP1)
            GL = MIN(GL, D(IEND)-TMP1)
            BNORM = MAX(ABS(GL), ABS(GU))
            GL = GL - FUDGE*BNORM*ULP*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*BNORM*ULP*IN + FUDGE*PIVMIN
!
!  Compute ATOLI for the current submatrix.
!
            IF (ABSTOL.LE.ZERO) THEN
               ATOLI = ULP*MAX(ABS(GL), ABS(GU))
            ELSE
               ATOLI = ABSTOL
            END IF
!
            IF (IRANGE.GT.1) THEN
               IF (GU.LE.WL) THEN
                  NWL = NWL + IN
                  NWU = NWU + IN
                  GOTO 70
               END IF
               GL = MAX(GL, WL)
               GU = MIN(GU, WU)
               IF (GL.GE.GU) GOTO 70
            END IF
!
!  Set Up Initial Interval.
!
            WORK(N+1) = GL
            WORK(N+IN+1) = GU
            CALL DLAEBZ (1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,     &
     &                   D(IBEGIN), E(IBEGIN), WORK(IBEGIN),            &
     &                   IDUMMA, WORK(N+1), WORK(N+2*IN+1), IM,         &
     &                   IWORK, W(M+1), IBLOCK(M+1), IINFO)
!
            NWL = NWL + IWORK(1)
            NWU = NWU + IWORK(IN+1)
            IWOFF = M - IWORK(1)
!
!  Compute Eigenvalues.
!
            ITMAX = INT((LOG(GU-GL+PIVMIN)-LOG(PIVMIN)) /               &
     &              LOG(TWO)) + 2
            CALL DLAEBZ (2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN, &
     &                   D(IBEGIN), E(IBEGIN), WORK(IBEGIN),            &
     &                   IDUMMA, WORK(N+1), WORK(N+2*IN+1), IOUT,       &
     &                   IWORK, W(M+1), IBLOCK(M+1), IINFO)
!
!  Copy Eigenvalues Into W and IBLOCK
!  Use -JB for block number for unconverged eigenvalues.
!
            DO 60 J = 1, IOUT
               TMP1 = HALF*(WORK(J+N)+WORK(J+IN+N))
!
!  Flag non-convergence.
!
               IF (J.GT.IOUT-IINFO) THEN
                  NCNVRG = .TRUE.
                  IB = -JB
               ELSE
                  IB = JB
               END IF
               DO 50 JE = IWORK(J)+1+IWOFF, IWORK(J+IN)+IWOFF
                  W(JE) = TMP1
                  IBLOCK(JE) = IB
   50          CONTINUE
   60       CONTINUE
!
            M = M + IM
         END IF
   70 CONTINUE
!
!  If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
!  If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
!
      IF (IRANGE.EQ.3) THEN
         IM = 0
         IDISCL = IL - 1 - NWL
         IDISCU = NWU - IU
!
         IF (IDISCL.GT.0 .OR. IDISCU.GT.0) THEN
            DO 80 JE = 1, M
               IF (W(JE).LE.WLU .AND. IDISCL.GT.0) THEN
                  IDISCL = IDISCL - 1
               ELSE IF (W(JE).GE.WUL .AND. IDISCU.GT.0) THEN
                  IDISCU = IDISCU - 1
               ELSE
                  IM = IM + 1
                  W(IM) = W(JE)
                  IBLOCK(IM) = IBLOCK(JE)
               END IF
   80       CONTINUE
            M = IM
         END IF
         IF (IDISCL.GT.0 .OR. IDISCU.GT.0) THEN
!
!  Code to deal with effects of bad arithmetic:
!  Some low eigenvalues to be discarded are not in (WL,WLU],
!  or high eigenvalues to be discarded are not in (WUL,WU]
!  so just kill off the smallest IDISCL/largest IDISCU
!  eigenvalues, by simply finding the smallest/largest
!  eigenvalue(s).
!
!  (If N(w) is monotone non-decreasing, this should never happen.)
!
            IF (IDISCL.GT.0) THEN
               WKILL = WU
               DO 100 JDISC = 1, IDISCL
                  IW = 0
                  DO 90 JE = 1, M
                     IF (IBLOCK(JE).NE.0 .AND.                          &
     &                   (W(JE).LT.WKILL .OR. IW.EQ.0)) THEN
                        IW = JE
                        WKILL = W(JE)
                     END IF
   90             CONTINUE
                  IF (IW .NE. 0) IBLOCK(IW) = 0
  100          CONTINUE
            END IF
!
            IF (IDISCU.GT.0) THEN
               WKILL = WL
               DO 120 JDISC = 1, IDISCU
                  IW = 0
                  DO 110 JE = 1, M
                     IF (IBLOCK(JE).NE.0 .AND.                          &
     &                   (W(JE).GT.WKILL .OR. IW.EQ.0)) THEN
                        IW = JE
                        WKILL = W(JE)
                     END IF
  110             CONTINUE
                  IF (IW .NE. 0) IBLOCK(IW) = 0
  120          CONTINUE
            END IF
!
            IM = 0
            DO 130 JE = 1, M
               IF (IBLOCK(JE).NE.0) THEN
                  IM = IM + 1
                  W(IM) = W(JE)
                  IBLOCK(IM) = IBLOCK(JE)
               END IF
  130       CONTINUE
            M = IM
         END IF

         IF (IDISCL.LT.0 .OR. IDISCU.LT.0) THEN
            TOOFEW = .TRUE.
         END IF
      END IF
!
!  If ORDER='B', do nothing -- the eigenvalues are already sorted
!  by block.
!  If ORDER='E', sort the eigenvalues from smallest to largest
!
      IF (IORDER.EQ.1 .AND. NSPLIT.GT.1) THEN
         DO 150 JE = 1, M - 1
            IE = 0
            TMP1 = W(JE)
            DO 140 J = JE + 1, M
               IF (W(J).LT.TMP1) THEN
                  IE = J
                  TMP1 = W(J)
               END IF
  140       CONTINUE
!
            IF (IE.NE.0) THEN
               ITMP1 = IBLOCK(IE)
               W(IE) = W(JE)
               IBLOCK(IE) = IBLOCK(JE)
               W(JE) = TMP1
               IBLOCK(JE) = ITMP1
            END IF
  150    CONTINUE
      END IF
!
      INFO = 0
      IF (NCNVRG) INFO = INFO + 1
      IF (TOOFEW) INFO = INFO + 2
!
      RETURN
      END SUBROUTINE DSTEBZ
!
      SUBROUTINE DSTEIN (N, D, E, M, W, IBLOCK, ISPLIT, Z, LDZ, WORK,   &
     &                   IWORK, IFAIL, INFO)
!
!=======================================================================
!                                                                      !
!  DSTEIN computes the eigenvectors of a real symmetric tridiagonal    !
!  matrix T corresponding to specified eigenvalues, using inverse      !
!  iteration.                                                          !
!                                                                      !
!  The maximum number of iterations allowed for each eigenvector is    !
!  specified by an internal parameter MAXITS (currently set to 5).     !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]    INTEGER                                            !
!           The order of the matrix.  N >= 0.                          !
!                                                                      !
!  D        [in]    REAL array, dimension (N)                          !
!           The n diagonal elements of the tridiagonal matrix T.       !
!                                                                      !
!  E        [in]    REAL array, dimension (N-1)                        !
!           The (n-1) subdiagonal elements of the tridiagonal matrix   !
!           T, in elements 1 to N-1.                                   !
!                                                                      !
!  M        [in]    INTEGER                                            !
!           The number of eigenvectors to be found.  0 <= M <= N.      !
!                                                                      !
!  W        [in]    REAL array, dimension (N)                          !
!           The first M elements of W contain the eigenvalues for      !
!           which eigenvectors are to be computed.  The eigenvalues    !
!           should be grouped by split-off block and ordered from      !
!           smallest to largest within the block.  (The output array   !
!           W from DSTEBZ with ORDER = 'B' is expected here.)          !
!                                                                      !
!  IBLOCK   [in]    INTEGER array, dimension (N)                       !
!           The submatrix indices associated with the corresponding    !
!           eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs   !
!           to the first submatrix from the top, =2 if W(i) belongs    !
!           to the second submatrix, etc.  (The output array IBLOCK    !
!           from DSTEBZ is expected here.)                             !
!                                                                      !
!  ISPLIT   [in]    INTEGER array, dimension (N)                       !
!           The splitting points, at which T breaks up into            !
!           submatrices. The first submatrix consists of rows/columns  !
!           1 to ISPLIT(1), the second of rows/columns ISPLIT(1)+1     !
!           through ISPLIT(2), etc.                                    !
!           (The output array ISPLIT from DSTEBZ is expected here.)    !
!                                                                      !
!  Z        [out]   REAL array, dimension (LDZ, M)                     !
!           The computed eigenvectors.  The eigenvector associated     !
!           with the eigenvalue W(i) is stored in the i-th column of   !
!           Z.  Any vector which fails to converge is set to its       !
!           current iterate after MAXITS iterations.                   !
!                                                                      !
!  LDZ      [in]    INTEGER                                            !
!           The leading dimension of the array Z.  LDZ >= max(1,N).    !
!                                                                      !
!  WORK     [out]   REAL array, dimension (5*N)                        !
!                                                                      !
!  IWORK    [out]   INTEGER array, dimension (N)                       !
!                                                                      !
!  IFAIL    [out]   INTEGER array, dimension (M)                       !
!           On normal exit, all elements of IFAIL are zero.            !
!           If one or more eigenvectors fail to converge after         !
!           MAXITS iterations, then their indices are stored in        !
!           array IFAIL.                                               !
!                                                                      !
!  INFO     [out]   INTEGER                                            !
!           = 0: successful exit.                                      !
!           < 0: if INFO = -i, the i-th argument had an illegal value  !
!           > 0: if INFO = i, then i eigenvectors failed to converge   !
!                in MAXITS iterations.  Their indices are stored in    !
!                array IFAIL.                                          !
!                                                                      !
!  Internal Parameters:                                                !
!                                                                      !
!  MAXITS   INTEGER, default = 5                                       !
!           The maximum number of iterations performed.                !
!                                                                      !
!  EXTRA    INTEGER, default = 2                                       !
!           The number of iterations performed after norm growth       !
!           criterion is satisfied, should be at least 1.              !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)  :: LDZ, M, N
      integer,   intent(in)  :: IBLOCK(*), ISPLIT(*)
      integer,   intent(out) :: INFO
      integer,   intent(out) :: IFAIL(*), IWORK(*)
!
      real (dp), intent(in)  ::  D(*), E(*), W(*)
      real (dp), intent(out) :: WORK(*), Z(LDZ,*)
!
!  Imported variable declarations.
!
      integer, parameter :: MAXITS = 5
      integer, parameter :: EXTRA = 2
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TEN = 10.0_dp
      real (dp), parameter :: ODM3 = 1.0E-3_dp
      real (dp), parameter :: ODM1 = 1.0E-1_dp
!
      integer :: B1, BLKSIZ, BN, GPIND, I, IINFO, INDRV1
      integer :: INDRV2, INDRV3, INDRV4, INDRV5, ITS, J, J1
      integer :: JBLK, JMAX, NBLK, NRMCHK
!
      integer, dimension(4) :: ISEED
!
      real (dp) :: DTPCRT, EPS, EPS1, NRM, ONENRM, ORTOL, PERTOL
      real (dp) :: SCL, SEP, TOL, XJ, XJM, ZTR
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input parameters.
!
      INFO = 0
      DO I = 1, M
         IFAIL(I) = 0
      END DO
!
      IF (N.LT.0) THEN
         INFO = -1
      ELSE IF (M.LT.0 .OR. M.GT.N) THEN
         INFO = -4
      ELSE IF (LDZ.LT.MAX(1, N)) THEN
         INFO = -9
      ELSE
         DO J = 2, M
            IF (IBLOCK(J).LT.IBLOCK(J-1)) THEN
               INFO = -6
               GOTO 30
            END IF
            IF (IBLOCK(J).EQ.IBLOCK(J-1) .AND. W(J).LT.W(J-1)) THEN
               INFO = -5
               GOTO 30
            END IF
         END DO
   30    CONTINUE
      END IF
!
      IF (INFO.NE.0) THEN
         CALL XERBLA ('DSTEIN', -INFO)
         RETURN
      END IF
!
!  Quick return if possible.
!
      IF (N.EQ.0 .OR. M.EQ.0) THEN
         RETURN
      ELSE IF (N.EQ.1) THEN
         Z(1, 1) = ONE
         RETURN
      END IF
!
!  Get machine constants.
!
      EPS = DLAMCH('Precision')
!
!  Initialize seed for random number generator DLARNV.
!
      DO I = 1, 4
         ISEED(I) = 1
      END DO
!
!  Initialize pointers.
!
      INDRV1 = 0
      INDRV2 = INDRV1 + N
      INDRV3 = INDRV2 + N
      INDRV4 = INDRV3 + N
      INDRV5 = INDRV4 + N
!
!  Compute eigenvectors of matrix blocks.
!
      J1 = 1
      DO 160 NBLK = 1, IBLOCK(M)
!
!  Find starting and ending indices of block nblk.
!
         IF (NBLK.EQ.1) THEN
            B1 = 1
         ELSE
            B1 = ISPLIT(NBLK-1) + 1
         END IF
         BN = ISPLIT(NBLK)
         BLKSIZ = BN - B1 + 1
         IF (BLKSIZ.EQ.1) GOTO 60
         GPIND = J1
!
!  Compute reorthogonalization criterion and stopping criterion.
!
         ONENRM = ABS(D(B1)) + ABS(E(B1))
         ONENRM = MAX(ONENRM, ABS(D(BN))+ABS(E(BN-1)))
         DO 50 I = B1 + 1, BN - 1
            ONENRM = MAX(ONENRM, ABS(D(I))+ABS(E(I-1))+ABS(E(I)))
   50    CONTINUE
         ORTOL = ODM3*ONENRM
!
         DTPCRT = SQRT(ODM1 / BLKSIZ)
!
!  Loop through eigenvalues of block nblk.
!
   60    CONTINUE
         JBLK = 0
         DO 150 J = J1, M
            IF (IBLOCK(J).NE.NBLK) THEN
               J1 = J
               GOTO 160
            END IF
            JBLK = JBLK + 1
            XJ = W(J)
!
!  Skip all the work if the block size is one.
!
            IF (BLKSIZ.EQ.1) THEN
               WORK(INDRV1+1) = ONE
               GOTO 120
            END IF
!
!  If eigenvalues j and j-1 are too close, add a relatively
!  small perturbation.
!
            IF (JBLK.GT.1) THEN
               EPS1 = ABS(EPS*XJ)
               PERTOL = TEN*EPS1
               SEP = XJ - XJM
               IF (SEP.LT.PERTOL) THEN
                 XJ = XJM + PERTOL
               END IF
            END IF
!
            ITS = 0
            NRMCHK = 0
!
!  Get random starting vector.
!
            CALL DLARNV (2, ISEED, BLKSIZ, WORK(INDRV1+1))
!
!  Copy the matrix T so it will not be destroyed in factorization.
!
            CALL DCOPY (BLKSIZ, D(B1), 1, WORK(INDRV4+1), 1)
            CALL DCOPY (BLKSIZ-1, E(B1), 1, WORK(INDRV2+2), 1)
            CALL DCOPY (BLKSIZ-1, E(B1), 1, WORK(INDRV3+1), 1)
!
!  Compute LU factors with partial pivoting  (PT = LU).
!
            TOL = ZERO
            CALL DLAGTF (BLKSIZ, WORK(INDRV4+1), XJ,                    &
     &                   WORK(INDRV2+2),                                &
     &                   WORK(INDRV3+1), TOL, WORK(INDRV5+1), IWORK,    &
     &                   IINFO)
!
!  Update iteration count.
!
   70       CONTINUE
            ITS = ITS + 1
            IF (ITS.GT.MAXITS) GOTO 100
!
!  Normalize and scale the righthand side vector Pb.
!
            JMAX = IDAMAX(BLKSIZ, WORK(INDRV1+1), 1)
            SCL = BLKSIZ*ONENRM*                                        &
     &            MAX(EPS, ABS(WORK(INDRV4+BLKSIZ))) /                  &
     &            ABS(WORK(INDRV1+JMAX))
            CALL DSCAL (BLKSIZ, SCL, WORK(INDRV1+1), 1)
!
!  Solve the system LU = Pb.
!
            CALL DLAGTS (-1, BLKSIZ, WORK(INDRV4+1),                    &
     &                   WORK(INDRV2+2),                                &
     &                   WORK(INDRV3+1), WORK(INDRV5+1), IWORK,         &
     &                   WORK(INDRV1+1), TOL, IINFO)
!
!  Reorthogonalize by modified Gram-Schmidt if eigenvalues are
!  close enough.
!
            IF (JBLK.EQ.1) GOTO 90

            IF (ABS(XJ-XJM).GT.ORTOL) THEN
               GPIND = J
            END IF

            IF (GPIND.NE.J) THEN
               DO 80 I = GPIND, J - 1
                  ZTR = -DDOT(BLKSIZ, WORK(INDRV1+1), 1, Z(B1,I), 1)
                  CALL DAXPY (BLKSIZ, ZTR, Z(B1,I), 1,                  &
     &                        WORK(INDRV1+1), 1)
   80          CONTINUE
            END IF
!
!  Check the infinity norm of the iterate.
!
   90       CONTINUE
            JMAX = IDAMAX(BLKSIZ, WORK(INDRV1+1), 1)
            NRM = ABS(WORK(INDRV1+JMAX))
!
!  Continue for additional iterations after norm reaches
!  stopping criterion.
!
            IF (NRM.LT.DTPCRT) GOTO 70
            NRMCHK = NRMCHK + 1
            IF (NRMCHK.LT.EXTRA+1) GOTO 70
!
            GOTO 110
!
!  If stopping criterion was not satisfied, update info and
!  store eigenvector number in array ifail.
!
  100       CONTINUE
            INFO = INFO + 1
            IFAIL(INFO) = J
!
!  Accept iterate as jth eigenvector.
!
  110       CONTINUE
            SCL = ONE / DNRM2(BLKSIZ, WORK(INDRV1+1), 1)
            JMAX = IDAMAX(BLKSIZ, WORK(INDRV1+1), 1)
            IF (WORK(INDRV1+JMAX).LT.ZERO) SCL = -SCL
            CALL DSCAL (BLKSIZ, SCL, WORK(INDRV1+1), 1)
  120       CONTINUE
            DO 130 I = 1, N
               Z(I, J) = ZERO
  130       CONTINUE
            DO 140 I = 1, BLKSIZ
               Z(B1+I-1, J) = WORK(INDRV1+I)
  140       CONTINUE
!
!  Save the shift to check eigenvalue spacing at next
!  iteration.
!
            XJM = XJ
!
  150    CONTINUE
  160 CONTINUE
!
      RETURN
      END SUBROUTINE DSTEIN
!
      SUBROUTINE DSTEMR (JOBZ, RANGE, N, D, E, VL, VU, IL, IU,          &
     &                   M, W, Z, LDZ, NZC, ISUPPZ, TRYRAC,             &
     &                   WORK, LWORK, IWORK, LIWORK, INFO)
!
!=======================================================================
!                                                                      !
!  DSTEMR computes selected eigenvalues and, optionally, eigenvectors  !
!  of a real symmetric tridiagonal matrix T. Any such unreduced matrix !
!  has a well defined set of pairwise different real eigenvalues, the  !
!  corresponding real eigenvectors are pairwise orthogonal.            !
!                                                                      !
!  The spectrum may be computed either completely or partially by      !
!  specifying either an interval (VL,VU] or a range of indices IL:IU   !
!  for the desired eigenvalues.                                        !
!                                                                      !
!  Depending on the number of desired eigenvalues, these are computed  !
!  either by bisection or the dqds algorithm. Numerically orthogonal   !
!  eigenvectors are computed by the use of various suitable L D L^T    !
!  factorizations near clusters of close eigenvalues (referred to as   !
!  RRRs, Relatively Robust Representations). An informal sketch of the !
!  algorithm follows.                                                  !
!                                                                      !
!  For each unreduced block (submatrix) of T,                          !
!     (a) Compute T - sigma I  = L D L^T, so that L and D              !
!         define all the wanted eigenvalues to high relative accuracy. !
!         This means that small relative changes in the entries of D   !
!         and L cause only small relative changes in the eigenvalues   !
!         and eigenvectors. The standard (unfactored) representation   !
!         of the tridiagonal matrix T does not have this property in   !
!         general.                                                     !
!     (b) Compute the eigenvalues to suitable accuracy.                !
!         If the eigenvectors are desired, the algorithm attains full  !
!         accuracy of the computed eigenvalues only right before       !
!         the corresponding vectors have to be computed, see steps c)  !
!         and d).                                                      !
!     (c) For each cluster of close eigenvalues, select a new shift    !
!         close to the cluster, find a new factorization, and refine   !
!         the shifted eigenvalues to suitable accuracy.                !
!     (d) For each eigenvalue with a large enough relative separation  !
!         compute the corresponding eigenvector by forming a rank      !
!         revealing twisted factorization. Go back to (c) for any      !
!         clusters that remain.                                        !
!                                                                      !
!  For more details, see:                                              !
!                                                                      !
!  - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple           !
!    representations to compute orthogonal eigenvectors of symmetric   !
!    tridiagonal matrices," Linear Algebra and its Applications,       !
!    387(1), pp. 1-28, August 2004.                                    !
!                                                                      !
!  - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors  !
!    and Relative Gaps," SIAM Journal on Matrix Analysis and           !
!    Applications, Vol. 25, 2004.  Also LAPACK Working Note 154.       !
!                                                                      !
!  - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric       !
!    tridiagonal eigenvalue/eigenvector problem", Computer Science     !
!    Division Technical Report No. UCB/CSD-97-971, UC Berkeley,        !
!    May 1997.                                                         !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  1.DSTEMR works only on machines which follow IEEE-754               !
!  floating-point standard in their handling of infinities and NaNs.   !
!  This permits the use of efficient inner loops avoiding a check for  !
!  zero divisors.                                                      !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  JOBZ     [in]    CHARACTER*1                                        !
!           = 'N':  Compute eigenvalues only;                          !
!           = 'V':  Compute eigenvalues and eigenvectors.              !
!                                                                      !
!  RANGE    [in]    CHARACTER*1                                        !
!           = 'A': all eigenvalues will be found.                      !
!           = 'V': all eigenvalues in the half-open interval (VL,VU]   !
!                  will be found.                                      !
!           = 'I': the IL-th through IU-th eigenvalues will be found.  !
!                                                                      !
!  N        [in]    INTEGER                                            !
!           The order of the matrix.  N >= 0.                          !
!                                                                      !
!  D        [in,out] REAL array, dimension (N)                         !
!           On entry, the N diagonal elements of the tridiagonal       !
!           matrix T. On exit, D is overwritten.                       !
!                                                                      !
!  E        [in,out] REAL array, dimension (N)                         !
!           On entry, the (N-1) subdiagonal elements of the            !
!           tridiagonal matrix T in elements 1 to N-1 of E. E(N) need  !
!           not be set on input, but is used internally as workspace.  !
!           On exit, E is overwritten.                                 !
!                                                                      !
!  VL       [in]    REAL                                               !
!           If RANGE='V', the lower bound of the interval to           !
!           be searched for eigenvalues. VL < VU.                      !
!           Not referenced if RANGE = 'A' or 'I'.                      !
!                                                                      !
!  VU       [in]    REAL                                               !
!           If RANGE='V', the upper bound of the interval to           !
!           be searched for eigenvalues. VL < VU.                      !
!           Not referenced if RANGE = 'A' or 'I'.                      !
!                                                                      !
!  IL       [in]    INTEGER                                            !
!           If RANGE='I', the index of the                             !
!           smallest eigenvalue to be returned.                        !
!           1 <= IL <= IU <= N, if N > 0.                              !
!           Not referenced if RANGE = 'A' or 'V'.                      !
!                                                                      !
!  IU       [in]    INTEGER                                            !
!           If RANGE='I', the index of the                             !
!           largest eigenvalue to be returned.                         !
!           1 <= IL <= IU <= N, if N > 0.                              !
!           Not referenced if RANGE = 'A' or 'V'.                      !
!                                                                      !
!  M        [out]   INTEGER                                            !
!           The total number of eigenvalues found.  0 <= M <= N.       !
!           If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.    !
!                                                                      !
!  W        [out]   REAL array, dimension (N)                          !
!           The first M elements contain the selected eigenvalues in   !
!           ascending order.                                           !
!                                                                      !
!  Z        [out]   REAL array, dimension (LDZ, max(1,M))              !
!           If JOBZ = 'V', and if INFO = 0, then the first M columns   !
!           of Z contain the orthonormal eigenvectors of the matrix T  !
!           corresponding to the selected eigenvalues, with the i-th   !
!           column of Z holding the eigenvector associated with W(i).  !
!           If JOBZ = 'N', then Z is not referenced.                   !
!                                                                      !
!           Note: the user must ensure that at least max(1,M) columns  !
!           are supplied in the array Z; if RANGE = 'V', the exact     !
!           value of M is not known in advance and can be computed     !
!           with a workspace query by setting NZC = -1, see below.     !
!                                                                      !
!  LDZ      [in]    INTEGER                                            !
!           The leading dimension of the array Z.  LDZ >= 1, and if    !
!           JOBZ = 'V', then LDZ >= max(1,N).                          !
!                                                                      !
!  NZC      [in]    INTEGER                                            !
!           The number of eigenvectors to be held in the array Z.      !
!           If RANGE = 'A', then NZC >= max(1,N).                      !
!           If RANGE = 'V', then NZC >= the number of eigenvalues in   !
!                                       (VL,VU].                       !
!           If RANGE = 'I', then NZC >= IU-IL+1.                       !
!                                                                      !
!           If NZC = -1, then a workspace query is assumed; the        !
!           routine calculates the number of columns of the array Z    !
!           that are needed to hold the eigenvectors.                  !
!           This value is returned as the first entry of the Z array,  !
!           and no error message related to NZC is issued by XERBLA.   !
!                                                                      !
!  ISUPPZ   [out]   INTEGER array, dimension (2*max(1,M))              !
!           The support of the eigenvectors in Z, i.e., the indices    !
!           indicating the nonzero elements in Z. The i-th computed    !
!           eigenvector is nonzero only in elements ISUPPZ(2*i-1)      !
!           through ISUPPZ(2*i). This is relevant in the case when     !
!           the matrix is split. ISUPPZ is only accessed when JOBZ is  !
!           'V' and N > 0.                                             !
!                                                                      !
!  TRYRAC   [in,out]  LOGICAL                                          !
!           If TRYRAC = .TRUE., indicates that the code should check   !
!           whether the tridiagonal matrix defines its eigenvalues to  !
!           high relative accuracy.  If so, the code uses relative-    !
!           accuracy preserving algorithms that might be (a bit)       !
!           slower depending on the matrix. If the matrix does not     !
!            define its eigenvalues to high relative accuracy, the     !
!           code can uses possibly faster algorithms.                  !
!           If TRYRAC = .FALSE., the code is not required to guarantee !
!           relatively accurate eigenvalues and can use the fastest    !
!           possible techniques.                                       !
!           On exit, a .TRUE. TRYRAC will be set to .FALSE. if the     !
!           matrix does not define its eigenvalues to high relative    !
!           accuracy.                                                  !
!                                                                      !
!  WORK     [out]   REAL array, dimension (LWORK)                      !
!           On exit, if INFO = 0, WORK(1) returns the optimal          !
!           (and minimal) LWORK.                                       !
!                                                                      !
!  LWORK    [in]    INTEGER                                            !
!           The dimension of the array WORK. LWORK >= max(1,18*N)      !
!           if JOBZ = 'V', and LWORK >= max(1,12*N) if JOBZ = 'N'.     !
!           If LWORK = -1, then a workspace query is assumed; the      !
!           routine only calculates the optimal size of the WORK       !
!           array, returns this value as the first entry of the WORK   !
!           array, and no error message related to LWORK is issued by  !
!           XERBLA.                                                    !
!                                                                      !
!  IWORK    [out]   INTEGER array, dimension (LIWORK)                  !
!           On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK. !
!                                                                      !
!  LIWORK   [in]    INTEGER                                            !
!           The dimension of the array IWORK.  LIWORK >= max(1,10*N)   !
!           if the eigenvectors are desired, and LIWORK >= max(1,8*N)  !
!           if only the eigenvalues are to be computed.                !
!           If LIWORK = -1, then a workspace query is assumed; the     !
!           routine only calculates the optimal size of the IWORK      !
!           array, returns this value as the first entry of the IWORK  !
!           array, and no error message related to LIWORK is issued by !
!           XERBLA.                                                    !
!                                                                      !
!  INFO     [out]   INTEGER                                            !
!           On exit, INFO                                              !
!           = 0:  successful exit                                      !
!           < 0:  if INFO = -i, the i-th argument had an illegal value !
!           > 0:  if INFO = 1X, internal error in DLARRE,              !
!                 if INFO = 2X, internal error in DLARRV.              !
!                 Here, the digit X = ABS(IINFO) < 10, where IINFO     !
!                 is the nonzero error code returned by DLARRE or      !
!                 DLARRV, respectively.                                !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!  Aravindh Krishnamoorthy, FAU, Erlangen, Germany                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), intent(in) :: JOBZ, RANGE
!
      logical,   intent(inout) :: TRYRAC
!
      integer,   intent(in)    :: IL, IU, LDZ, LIWORK, LWORK, N, NZC
      integer,   intent(out)   :: INFO, M
      integer,   intent(out)   :: ISUPPZ(*), IWORK(*)
!
      real (dp), intent(in)    :: VL, VU
      real (dp), intent(inout) :: D(*), E(*)
      real (dp), intent(out)   :: W(*), WORK(*), Z(LDZ,*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: FOUR = 4.0_dp
      real (dp), parameter :: MINRGP = 1.0E-3_dp
!
      logical :: ALLEIG, INDEIG, LQUERY, VALEIG, WANTZ, ZQUERY, LAESWAP
!
      integer :: I, IBEGIN, IEND, IFIRST, IIL, IINDBL, IINDW
      integer :: IINDWK, IINFO, IINSPL, IIU, ILAST, IN, INDD
      integer :: INDE2, INDERR, INDGP, INDGRS, INDWRK, ITMP
      integer :: ITMP2, J, JBLK, JJ, LIWMIN, LWMIN, NSPLIT
      integer :: NZCMIN, OFFSET, WBEGIN, WEND
!
      real (dp) :: BIGNUM, CS, EPS, PIVMIN, R1, R2, RMAX, RMIN
      real (dp) :: RTOL1, RTOL2, SAFMIN, SCALE, SMLNUM, SN
      real (dp) :: THRESH, TMP, TNRM, WL, WU
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input parameters.
!
      WANTZ = LSAME(JOBZ, 'V')
      ALLEIG = LSAME(RANGE, 'A')
      VALEIG = LSAME(RANGE, 'V')
      INDEIG = LSAME(RANGE, 'I')
!
      LQUERY = ((LWORK.EQ.-1).OR.(LIWORK.EQ.-1))
      ZQUERY = (NZC.EQ.-1)
      LAESWAP = .FALSE.
!
!  DSTEMR needs WORK of size 6*N, IWORK of size 3*N.
!  In addition, DLARRE needs WORK of size 6*N, IWORK of size 5*N.
!  Furthermore, DLARRV needs WORK of size 12*N, IWORK of size 7*N.
!
      IF (WANTZ) THEN
         LWMIN = 18*N
         LIWMIN = 10*N
      ELSE
!
!  Need less workspace if only the eigenvalues are wanted
!
         LWMIN = 12*N
         LIWMIN = 8*N
      END IF

      WL = ZERO
      WU = ZERO
      IIL = 0
      IIU = 0
      NSPLIT = 0

      IF (VALEIG) THEN
!
!  We do not reference VL, VU in the cases RANGE = 'I','A'
!  The interval (WL, WU] contains all the wanted eigenvalues.
!  It is either given by the user or computed in DLARRE.
!
         WL = VL
         WU = VU
      ELSE IF (INDEIG) THEN
!
!  We do not reference IL, IU in the cases RANGE = 'V','A'
!
         IIL = IL
         IIU = IU
      END IF
!
      INFO = 0
      IF (.NOT.(WANTZ .OR. LSAME(JOBZ, 'N'))) THEN
         INFO = -1
      ELSE IF (.NOT.(ALLEIG .OR. VALEIG .OR. INDEIG)) THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      ELSE IF (VALEIG .AND. N.GT.0 .AND. WU.LE.WL) THEN
         INFO = -7
      ELSE IF (INDEIG .AND. (IIL.LT.1 .OR. IIL.GT.N)) THEN
         INFO = -8
      ELSE IF (INDEIG .AND. (IIU.LT.IIL .OR. IIU.GT.N)) THEN
         INFO = -9
      ELSE IF (LDZ.LT.1 .OR. (WANTZ .AND. LDZ.LT.N)) THEN
         INFO = -13
      ELSE IF (LWORK.LT.LWMIN .AND. .NOT.LQUERY) THEN
         INFO = -17
      ELSE IF (LIWORK.LT.LIWMIN .AND. .NOT.LQUERY) THEN
         INFO = -19
      END IF
!
!  Get machine constants.
!
      SAFMIN = DLAMCH('Safe minimum')
      EPS = DLAMCH('Precision')
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT(SMLNUM)
      RMAX = MIN(SQRT(BIGNUM), ONE / SQRT(SQRT(SAFMIN)))
!
      IF (INFO.EQ.0) THEN
         WORK(1) = LWMIN
         IWORK(1) = LIWMIN
!
         IF (WANTZ .AND. ALLEIG) THEN
            NZCMIN = N
         ELSE IF (WANTZ .AND. VALEIG) THEN
            CALL DLARRC ('T', N, VL, VU, D, E, SAFMIN,                  &
     &                   NZCMIN, ITMP, ITMP2, INFO)
         ELSE IF (WANTZ .AND. INDEIG) THEN
            NZCMIN = IIU-IIL+1
         ELSE
            NZCMIN = 0                ! WANTZ .EQ. FALSE.
         END IF
         IF (ZQUERY .AND. INFO.EQ.0) THEN
            Z(1,1) = NZCMIN
         ELSE IF (NZC.LT.NZCMIN .AND. .NOT.ZQUERY) THEN
            INFO = -14
         END IF
      END IF

      IF (INFO.NE.0) THEN
!
         CALL XERBLA ('DSTEMR', -INFO)
!
         RETURN
      ELSE IF (LQUERY .OR. ZQUERY) THEN
         RETURN
      END IF
!
!  Handle N = 0, 1, and 2 cases immediately.
!
      M = 0
      IF (N.EQ.0) RETURN
!
      IF (N.EQ.1) THEN
         IF (ALLEIG .OR. INDEIG) THEN
            M = 1
            W(1) = D(1)
         ELSE
            IF (WL.LT.D(1) .AND. WU.GE.D(1)) THEN
               M = 1
               W(1) = D(1)
            END IF
         END IF
         IF (WANTZ.AND.(.NOT.ZQUERY)) THEN
            Z(1, 1) = ONE
            ISUPPZ(1) = 1
            ISUPPZ(2) = 1
         END IF
         RETURN
      END IF
!
      IF (N.EQ.2) THEN
         IF (.NOT.WANTZ) THEN
            CALL DLAE2 (D(1), E(1), D(2), R1, R2)
         ELSE IF (WANTZ.AND.(.NOT.ZQUERY)) THEN
            CALL DLAEV2 (D(1), E(1), D(2), R1, R2, CS, SN)
         END IF
!
!  D/S/LAE2 and D/S/LAEV2 outputs satisfy |R1| >= |R2|. However,
!  the following code requires R1 >= R2. Hence, we correct
!  the order of R1, R2, CS, SN if R1 < R2 before further processing.
!
         IF (R1.LT.R2) THEN
            E(2) = R1
            R1 = R2
            R2 = E(2)
            LAESWAP = .TRUE.
         END IF
         IF (ALLEIG.OR.                                                 &
     &       (VALEIG.AND.(R2.GT.WL).AND.                                &
     &                   (R2.LE.WU)).OR.                                &
     &       (INDEIG.AND.(IIL.EQ.1))) THEN
            M = M+1
            W(M) = R2
            IF (WANTZ.AND.(.NOT.ZQUERY)) THEN
               IF (LAESWAP) THEN
                  Z(1, M) = CS
                  Z(2, M) = SN
               ELSE
                  Z(1, M) = -SN
                  Z(2, M) = CS
               END IF
!
!  Note: At most one of SN and CS can be zero.
!
               IF (SN.NE.ZERO) THEN
                  IF (CS.NE.ZERO) THEN
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 2
                  ELSE
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 1
                  END IF
               ELSE
                  ISUPPZ(2*M-1) = 2
                  ISUPPZ(2*M) = 2
               END IF
            END IF
         END IF
         IF (ALLEIG.OR.                                                 &
     &      (VALEIG.AND.(R1.GT.WL).AND.                                 &
     &                  (R1.LE.WU)).OR.                                 &
     &      (INDEIG.AND.(IIU.EQ.2))) THEN
            M = M+1
            W(M) = R1
            IF (WANTZ.AND.(.NOT.ZQUERY)) THEN
               IF (LAESWAP) THEN
                  Z(1, M) = -SN
                  Z(2, M) = CS
               ELSE
                  Z(1, M) = CS
                  Z(2, M) = SN
               END IF
!
!  Note: At most one of SN and CS can be zero.
!
               IF (SN.NE.ZERO) THEN
                  IF (CS.NE.ZERO) THEN
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 2
                  ELSE
                     ISUPPZ(2*M-1) = 1
                     ISUPPZ(2*M) = 1
                  END IF
               ELSE
                  ISUPPZ(2*M-1) = 2
                  ISUPPZ(2*M) = 2
               END IF
            END IF
         END IF

      ELSE
!
!  Continue with general N.
!
         INDGRS = 1
         INDERR = 2*N + 1
         INDGP = 3*N + 1
         INDD = 4*N + 1
         INDE2 = 5*N + 1
         INDWRK = 6*N + 1
!
         IINSPL = 1
         IINDBL = N + 1
         IINDW = 2*N + 1
         IINDWK = 3*N + 1
!
!  Scale matrix to allowable range, if necessary.
!  The allowable range is related to the PIVMIN parameter; see the
!  comments in DLARRD.  The preference for scaling small values
!  up is heuristic; we expect users matrices not to be close to the
!  RMAX threshold.
!
         SCALE = ONE
         TNRM = DLANST('M', N, D, E)
         IF (TNRM.GT.ZERO .AND. TNRM.LT.RMIN) THEN
            SCALE = RMIN / TNRM
         ELSE IF (TNRM.GT.RMAX) THEN
            SCALE = RMAX / TNRM
         END IF
         IF (SCALE.NE.ONE) THEN
            CALL DSCAL (N, SCALE, D, 1)
            CALL DSCAL (N-1, SCALE, E, 1)
            TNRM = TNRM*SCALE
            IF (VALEIG) THEN
!
!  If eigenvalues in interval have to be found,
!  scale (WL, WU] accordingly.
!
               WL = WL*SCALE
               WU = WU*SCALE
            END IF
         END IF
!
!  Compute the desired eigenvalues of the tridiagonal after splitting
!  into smaller subblocks if the corresponding off-diagonal elements
!  are small
!  THRESH is the splitting parameter for DLARRE
!  A negative THRESH forces the old splitting criterion based on the
!  size of the off-diagonal. A positive THRESH switches to splitting
!  which preserves relative accuracy.
!
         IF (TRYRAC) THEN
!
!  Test whether the matrix warrants the more expensive relative
!  approach.
!
            CALL DLARRR (N, D, E, IINFO)
         ELSE
!
!  The user does not care about relative accurately eigenvalues.
!
            IINFO = -1
         END IF
!
!  Set the splitting criterion.
!
         IF (IINFO.EQ.0) THEN
            THRESH = EPS
         ELSE
            THRESH = -EPS
!
!  Relative accuracy is desired but T does not guarantee it.
!
            TRYRAC = .FALSE.
         END IF
!
         IF (TRYRAC) THEN
!
!  Copy original diagonal, needed to guarantee relative accuracy.
!
            CALL DCOPY (N,D,1,WORK(INDD),1)
         END IF
!
!  Store the squares of the offdiagonal values of T.
!
         DO J = 1, N-1
            WORK(INDE2+J-1) = E(J)**2
         END DO
!
!  Set the tolerance parameters for bisection.
!
         IF (.NOT.WANTZ) THEN

!           DLARRE computes the eigenvalues to full precision.

            RTOL1 = FOUR * EPS
            RTOL2 = FOUR * EPS
         ELSE
!
!  DLARRE computes the eigenvalues to less than full precision.
!  DLARRV will refine the eigenvalue approximations, and we can
!  need less accurate initial bisection in DLARRE.
!  Note: these settings do only affect the subset case and DLARRE
!
            RTOL1 = SQRT(EPS)
            RTOL2 = MAX(SQRT(EPS)*5.0E-3_dp, FOUR * EPS)
         END IF

         CALL DLARRE (RANGE, N, WL, WU, IIL, IIU, D, E,                 &
     &                WORK(INDE2), RTOL1, RTOL2, THRESH, NSPLIT,        &
     &                IWORK(IINSPL), M, W, WORK(INDERR),                &
     &                WORK(INDGP), IWORK(IINDBL),                       &
     &                IWORK(IINDW), WORK(INDGRS), PIVMIN,               &
     &                WORK(INDWRK), IWORK(IINDWK), IINFO)
         IF (IINFO.NE.0) THEN
            INFO = 10 + ABS(IINFO)
            RETURN
         END IF
!
!  Note that if RANGE .NE. 'V', DLARRE computes bounds on the desired
!  part of the spectrum. All desired eigenvalues are contained in
!  (WL,WU]
!
         IF (WANTZ) THEN
!
!  Compute the desired eigenvectors corresponding to the computed
!  eigenvalues
!
            CALL DLARRV (N, WL, WU, D, E,                               &
     &                   PIVMIN, IWORK(IINSPL), M,                      &
     &                   1, M, MINRGP, RTOL1, RTOL2,                    &
     &                   W, WORK(INDERR), WORK(INDGP), IWORK(IINDBL),   &
     &                   IWORK(IINDW), WORK(INDGRS), Z, LDZ,            &
     &                   ISUPPZ, WORK(INDWRK), IWORK(IINDWK), IINFO)
            IF (IINFO.NE.0) THEN
               INFO = 20 + ABS(IINFO)
               RETURN
            END IF
         ELSE
!
!  DLARRE computes eigenvalues of the (shifted) root representation
!  DLARRV returns the eigenvalues of the unshifted matrix.
!  However, if the eigenvectors are not desired by the user, we need
!  to apply the corresponding shifts from DLARRE to obtain the
!  eigenvalues of the original matrix.
!
            DO 20 J = 1, M
               ITMP = IWORK(IINDBL+J-1)
               W(J) = W(J) + E(IWORK(IINSPL+ITMP-1))
 20         CONTINUE
         END IF
!
         IF (TRYRAC) THEN
!
!  Refine computed eigenvalues so that they are relatively accurate
!  with respect to the original matrix T.
!
            IBEGIN = 1
            WBEGIN = 1
            DO 39  JBLK = 1, IWORK(IINDBL+M-1)
               IEND = IWORK(IINSPL+JBLK-1)
               IN = IEND - IBEGIN + 1
               WEND = WBEGIN - 1
!
!  Check if any eigenvalues have to be refined in this block
!
 36            CONTINUE
               IF (WEND.LT.M) THEN
                  IF (IWORK(IINDBL+WEND).EQ.JBLK) THEN
                     WEND = WEND + 1
                     GOTO 36
                  END IF
               END IF
               IF (WEND.LT.WBEGIN) THEN
                  IBEGIN = IEND + 1
                  GOTO 39
               END IF

               OFFSET = IWORK(IINDW+WBEGIN-1)-1
               IFIRST = IWORK(IINDW+WBEGIN-1)
               ILAST = IWORK(IINDW+WEND-1)
               RTOL2 = FOUR * EPS

               CALL DLARRJ (IN,                                         &
     &                      WORK(INDD+IBEGIN-1), WORK(INDE2+IBEGIN-1),  &
     &                      IFIRST, ILAST, RTOL2, OFFSET, W(WBEGIN),    &
     &                      WORK(INDERR+WBEGIN-1),                      &
     &                      WORK(INDWRK), IWORK(IINDWK), PIVMIN,        &
     &                      TNRM, IINFO)
               IBEGIN = IEND + 1
               WBEGIN = WEND + 1
 39         CONTINUE
         END IF
!
!  If matrix was scaled, then rescale eigenvalues appropriately.
!
         IF (SCALE.NE.ONE) THEN
            CALL DSCAL (M, ONE / SCALE, W, 1)
         END IF

      END IF

!
!  If eigenvalues are not in increasing order, then sort them,
!  possibly along with eigenvectors.
!
      IF (NSPLIT.GT.1 .OR. N.EQ.2) THEN
         IF (.NOT. WANTZ) THEN
            CALL DLASRT ('I', M, W, IINFO)
            IF (IINFO.NE.0) THEN
               INFO = 3
               RETURN
            END IF
         ELSE
            DO J = 1, M - 1
               I = 0
               TMP = W(J)
               DO JJ = J + 1, M
                  IF (W(JJ).LT.TMP) THEN
                     I = JJ
                     TMP = W(JJ)
                  END IF
               END DO
               IF (I.NE.0) THEN
                  W(I) = W(J)
                  W(J) = TMP
                  IF (WANTZ) THEN
                     CALL DSWAP (N, Z(1, I), 1, Z(1, J), 1)
                     ITMP = ISUPPZ(2*I-1)
                     ISUPPZ(2*I-1) = ISUPPZ(2*J-1)
                     ISUPPZ(2*J-1) = ITMP
                     ITMP = ISUPPZ(2*I)
                     ISUPPZ(2*I) = ISUPPZ(2*J)
                     ISUPPZ(2*J) = ITMP
                  END IF
               END IF
            END DO
         END IF
      END IF
!
      WORK(1) = LWMIN
      IWORK(1) = LIWMIN
!
      RETURN
      END SUBROUTINE DSTEMR
!
      SUBROUTINE DSTERF (N, D, E, INFO)
!
!=======================================================================
!                                                                      !
!  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix   !
!  using the Pal-Walker-Kahan variant of the QL or QR algorithm.       !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix.  N >= 0.                          !
!                                                                      !
!  D        [in,out]  REAL array, dimension (N)                        !
!           On entry, the N diagonal elements of the tridiagonal       !
!           matrix.                                                    !
!           On exit, if INFO = 0, the eigenvalues in ascending order.  !
!                                                                      !
!  E        [in,out]  REAL array, dimension (N-1)                      !
!           On entry, the (n-1) subdiagonal elements of the            !
!           tridiagonal matrix.                                        !
!           On exit, E has been destroyed.                             !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           = 0:  successful exit                                      !
!           < 0:  if INFO = -i, the i-th argument had an illegal value !
!           > 0:  the algorithm failed to find all of the eigenvalues  !
!                 in a total of 30*N iterations; if INFO = i, then i   !
!                 elements of E have not converged to zero.            !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: N
      integer,   intent(out)   :: INFO
!
      real (dp), intent(inout) :: D(*), E(*)
!
!  Local variable declarations.
!
      integer, parameter :: MAXIT = 30
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: THREE = 3.0_dp
!
      integer :: I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M, NMAXIT
!
      real (dp) :: ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC
      real (dp) :: OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN
      real (dp) :: SIGMA, SSFMAX, SSFMIN, RMAX
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input parameters.
!
      INFO = 0
!
!  Quick return if possible
!
      IF (N.LT.0) THEN
         INFO = -1
         CALL XERBLA ('DSTERF', -INFO)
         RETURN
      END IF
      IF (N.LE.1) RETURN
!
!  Determine the unit roundoff for this environment.
!
      EPS = DLAMCH('E')
      EPS2 = EPS**2
      SAFMIN = DLAMCH('S')
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT(SAFMAX) / THREE
      SSFMIN = SQRT(SAFMIN) / EPS2
      RMAX = DLAMCH('O')
!
!  Compute the eigenvalues of the tridiagonal matrix.
!
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
!
!  Determine where the matrix splits and choose QL or QR iteration
!  for each block, according to whether top or bottom diagonal
!  element is smaller.
!
      L1 = 1
!
   10 CONTINUE
      IF (L1.GT.N) GOTO 170

      IF (L1.GT.1) E(L1-1) = ZERO

      DO 20 M = L1, N - 1
         IF (ABS(E(M)).LE.(SQRT(ABS(D(M)))*SQRT(ABS(D(M+1))))*EPS) THEN
            E(M) = ZERO
            GOTO 30
         END IF
   20 CONTINUE
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF (LEND.EQ.L) GOTO 10
!
!  Scale submatrix in rows and columns L to LEND
!
      ANORM = DLANST('M', LEND-L+1, D(L), E(L))
      ISCALE = 0

      IF (ANORM.EQ.ZERO) GOTO 10

      IF ((ANORM.GT.SSFMAX)) THEN
         ISCALE = 1
         CALL DLASCL ('G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D(L), N,   &
     &                INFO)
         CALL DLASCL ('G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E(L), N,     &
     &                INFO)
      ELSE IF (ANORM.LT.SSFMIN) THEN
         ISCALE = 2
         CALL DLASCL ('G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D(L), N,   &
     &                INFO)
         CALL DLASCL ('G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E(L), N,     &
     &                INFO)
      END IF
!
      DO 40 I = L, LEND - 1
         E(I) = E(I)**2
   40 CONTINUE
!
!  Choose between QL and QR iteration.
!
      IF (ABS(D(LEND)).LT.ABS(D(L))) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF (LEND.GE.L) THEN
!
!  QL Iteration
!
!  Look for small subdiagonal element.
!
   50    CONTINUE
         IF (L.NE.LEND) THEN
            DO 60 M = L, LEND - 1
               IF (ABS(E(M)).LE.EPS2*ABS(D(M)*D(M+1))) GOTO 70
   60       CONTINUE
         END IF
         M = LEND
!
   70    CONTINUE
         IF (M.LT.LEND) E(M) = ZERO
         P = D(L)

         IF (M.EQ.L) GOTO 90
!
!  If remaining matrix is 2 by 2, use DLAE2 to compute its
!  eigenvalues.
!
         IF (M.EQ.L+1) THEN
            RTE = SQRT(E(L))
            CALL DLAE2 (D(L), RTE, D(L+1), RT1, RT2)
            D(L) = RT1
            D(L+1) = RT2
            E(L) = ZERO
            L = L + 2
            IF (L.LE.LEND) GOTO 50
            GOTO 150
         END IF
!
         IF (JTOT.EQ.NMAXIT) GOTO 150
         JTOT = JTOT + 1
!
!  Form shift.
!
         RTE = SQRT(E(L))
         SIGMA = (D(L+1)-P) / (TWO*RTE)
         R = DLAPY2(SIGMA, ONE)
         SIGMA = P - (RTE / (SIGMA+SIGN(R, SIGMA)))
!
         C = ONE
         S = ZERO
         GAMMA = D(M) - SIGMA
         P = GAMMA*GAMMA
!
!  Inner loop.
!
         DO 80 I = M - 1, L, -1
            BB = E(I)
            R = P + BB
            IF (I.NE.M-1) E(I+1) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D(I)
            GAMMA = C*(ALPHA-SIGMA) - S*OLDGAM
            D(I+1) = OLDGAM + (ALPHA-GAMMA)
            IF (C.NE.ZERO) THEN
               P = (GAMMA*GAMMA) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
!
         E(L) = S*P
         D(L) = SIGMA + GAMMA
         GOTO 50
!
!  Eigenvalue found.
!
   90    CONTINUE
         D(L) = P
!
         L = L + 1
         IF (L.LE.LEND) GOTO 50
         GOTO 150
!
      ELSE
!
!  QR Iteration
!
!  Look for small superdiagonal element.
!
  100    CONTINUE
         DO 110 M = L, LEND + 1, -1
            IF (ABS(E(M-1)).LE.EPS2*ABS(D(M)*D(M-1))) GOTO 120
  110    CONTINUE
         M = LEND
!
  120    CONTINUE
         IF (M.GT.LEND) E(M-1) = ZERO
         P = D(L)
         IF (M.EQ.L) GOTO 140
!
!  If remaining matrix is 2 by 2, use DLAE2 to compute its
!  eigenvalues.
!
         IF (M.EQ.L-1) THEN
            RTE = SQRT(E(L-1))
            CALL DLAE2 (D(L), RTE, D(L-1), RT1, RT2)
            D(L) = RT1
            D(L-1) = RT2
            E(L-1) = ZERO
            L = L - 2
            IF (L.GE.LEND) GOTO 100
            GOTO 150
         END IF
!
         IF (JTOT.EQ.NMAXIT) GOTO 150
         JTOT = JTOT + 1
!
!  Form shift.
!
         RTE = SQRT(E(L-1))
         SIGMA = (D(L-1)-P) / (TWO*RTE)
         R = DLAPY2(SIGMA, ONE)
         SIGMA = P - (RTE / (SIGMA+SIGN(R, SIGMA)))
!
         C = ONE
         S = ZERO
         GAMMA = D(M) - SIGMA
         P = GAMMA*GAMMA
!
!  Inner loop.
!
         DO 130 I = M, L - 1
            BB = E(I)
            R = P + BB
            IF (I.NE.M) E(I-1) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D(I+1)
            GAMMA = C*(ALPHA-SIGMA) - S*OLDGAM
            D(I) = OLDGAM + (ALPHA-GAMMA)
            IF (C.NE.ZERO) THEN
               P = (GAMMA*GAMMA) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
!
         E(L-1) = S*P
         D(L) = SIGMA + GAMMA
         GOTO 100
!
!  Eigenvalue found.
!
  140    CONTINUE
         D(L) = P
!
         L = L - 1
         IF (L.GE.LEND) GOTO 100
         GOTO 150
!
      END IF
!
!  Undo scaling if necessary.
!
  150 CONTINUE
      IF (ISCALE.EQ.1) THEN
         CALL DLASCL ('G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,        &
     &                D(LSV), N, INFO)
      END IF

      IF (ISCALE.EQ.2) THEN
         CALL DLASCL ('G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,        &
     &                D(LSV), N, INFO)
      END IF
!
!  Check for no convergence to an eigenvalue after a total
!  of N*MAXIT iterations.
!
      IF (JTOT.LT.NMAXIT) GOTO 10
      DO 160 I = 1, N - 1
         IF (E(I).NE.ZERO) INFO = INFO + 1
  160 CONTINUE
      GOTO 180
!
!  Sort eigenvalues in increasing order.
!
  170 CONTINUE
      CALL DLASRT ('I', N, D, INFO)
!
  180 CONTINUE
!
      RETURN
      END SUBROUTINE DSTERF
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!     BLAS Library Support Routines:
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
      SUBROUTINE DAXPY (N, DA, DX, INCX, DY, INCY)
!
!=======================================================================
!                                                                      !
!  DAXPY constant times a vector plus a vector. Uses unrolled loops    !
!  for increments equal to one.                                        !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           number of elements in input vector(s)                      !
!                                                                      !
!  DA       [in]      REAL                                             !
!           On entry, DA specifies the scalar alpha.                   !
!                                                                      !
!  DX       [in]      REAL array, dimension (1 + (N-1) * ABS(INCX))    !
!                                                                      !
!  INCX     [in]      INTEGER                                          !
!           storage spacing between elements of DX                     !
!                                                                      !
!  DY       [in,out]  REAL array, dimension (1 + (N-1) * ABS(INCY))    !
!                                                                      !
!  INCY     [in]      INTEGER                                          !
!           storage spacing between elements of DY                     !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!      Jack Dongarra, linpack, 3/11/78.                                !
!      modified 12/3/93, array(1) declarations changed to array(*)     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)  :: INCX, INCY, N
!
      real (dp), intent(in)  :: DA
      real (dp), intent(in)  :: DX(*)
      real (dp), intent(out) :: DY(*)
!
!  Local variable declarations.
!
      integer :: I, IX, IY, M, MP1
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0_dp) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!  Code for both increments equal to 1.
!
!  clean-up loop
!
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DY(I) + DA*DX(I)
            END DO
         END IF
         IF (N.LT.4) RETURN
         MP1 = M + 1
         DO I = MP1,N,4
            DY(I) = DY(I) + DA*DX(I)
            DY(I+1) = DY(I+1) + DA*DX(I+1)
            DY(I+2) = DY(I+2) + DA*DX(I+2)
            DY(I+3) = DY(I+3) + DA*DX(I+3)
         END DO
      ELSE
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
         END DO
      END IF
!
      RETURN
      END SUBROUTINE DAXPY
!
      SUBROUTINE DCOPY (N, DX, INCX, DY, INCY)
!
!=======================================================================
!                                                                      !
!  DCOPY copies a vector, x, to a vector, y. It uses unrolled loops    !
!  for increments equal to 1.                                          !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N       [in]     INTEGER                                            !
!          number of elements in input vector(s)                       !
!                                                                      !
!  DX      [in]     REAL array, dimension (1 + (N - 1)*abs(INCX))      !
!                                                                      !
!  INCX    [in]     INTEGER                                            !
!          storage spacing between elements of DX                      !
!                                                                      !
!  DY      [out]    REAL array, dimension (1 + (N - 1)*abs(INCY))      !
!                                                                      !
!  INCY    [in]     INTEGER                                            !
!          storage spacing between elements of DX                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)  :: INCX, INCY, N
!
      real (dp), intent(in)  :: DX(*)
      real (dp), intent(out) :: DY(*)
!
!  Local variable declarations.
!
      integer :: I, IX, IY, M, MP1
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (N.LE.0) RETURN
!
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!  Code for both increments equal to 1.
!
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,7
            DY(I) = DX(I)
            DY(I+1) = DX(I+1)
            DY(I+2) = DX(I+2)
            DY(I+3) = DX(I+3)
            DY(I+4) = DX(I+4)
            DY(I+5) = DX(I+5)
            DY(I+6) = DX(I+6)
         END DO
      ELSE
!
!  Code for unequal increments or equal increments not equal to 1.
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
!
      RETURN
      END SUBROUTINE DCOPY
!
      FUNCTION DDOT (N, DX, INCX, DY, INCY) RESULT (DOT)
!
!=======================================================================
!                                                                      !
!  DDOT forms the dot product of two vectors. Uses unrolled loops for  !
!  increments equal to one.                                            !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           number of elements in input vector(s)                      !
!                                                                      !
!  DX       [in]      REAL array, dimension (1 + (N-1) * ABS(INCX))    !
!                                                                      !
!  INCX     [in]      INTEGER                                          !
!           storage spacing between elements of DX                     !
!                                                                      !
!  DY       [in,out]  REAL array, dimension (1 + (N-1) * ABS(INCY))    !
!                                                                      !
!  INCY     [in]      INTEGER                                          !
!           storage spacing between elements of DY                     !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!      jack dongarra, linpack, 3/11/78.                                !
!      modified 12/3/93, array(1) declarations changed to array(*)     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in) :: INCX, INCY, N
!
      real (dp), intent(in) :: DX(*), DY(*)
!
!  Local variable declarations.
!
      integer   :: I, IX, IY, M, MP1
!
      real (dp) :: DOT, DTEMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      DOT = 0.0_dp
      DTEMP = 0.0_dp
!
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!  Code for both increments equal to 1.
!
!  Clean-up loop.
!
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF (N.LT.5) THEN
               DOT=DTEMP
            RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
           DTEMP = DTEMP + DX(I)*DY(I) +                                &
     &                     DX(I+1)*DY(I+1) +                            &
     &                     DX(I+2)*DY(I+2) +                            &
     &                     DX(I+3)*DY(I+3) +                            &
     &                     DX(I+4)*DY(I+4)
         END DO
      ELSE
!
!  Code for unequal increments or equal increments
!  not equal to 1.
!
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DOT = DTEMP
!
      RETURN
      END FUNCTION DDOT
!
      FUNCTION DNRM2 (N, X, INCX) RESULT (NRM2)
!
!=======================================================================
!                                                                      !
!  DNRM2 returns the euclidean norm of a vector via the function name, !
!  so that                                                             !
!                                                                      !
!     DNRM2 := SQRT(X'*X)                                              !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           number of elements in input vector(s)                      !
!                                                                      !
!  X        [in]      REAL array, dimension (1 + (N-1)*ABS(INCX))      !
!                                                                      !
!  INCX     [in]      INTEGER, storage spacing between elements of X   !
!           If INCX > 0, X(1+(i-1)*INCX) = x(i) for 1 <= i <= N        !
!           If INCX < 0, X(1-(N-i)*INCX) = x(i) for 1 <= i <= N        !
!           If INCX = 0, X is not a vector so there is no need to call !
!           this subroutine. If you call it anyway, it will count X(1) !
!           in the vector norm N times.                                !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Edward Anderson, Lockheed Martin                                    !
!  August 2016                                                         !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Weslley Pereira, University of Colorado Denver, USA                 !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  Anderson E. (2017)                                                  !
!  Algorithm 978: Safe Scaling in the Level 1 BLAS                     !
!  ACM Trans Math Softw 44:1--28                                       !
!  https://doi.org/10.1145/3061665                                     !
!                                                                      !
!  Blue, James L. (1978)                                               !
!  A Portable Fortran Program to Find the Euclidean Norm of a Vector   !
!  ACM Trans Math Softw 4:15--23                                       !
!  https://doi.org/10.1145/355769.355771                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in) :: INCX, N
!
      real (dp), intent(in) :: X(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE  = 1.0_dp
      real (dp), parameter :: MAXN = HUGE(0.0_dp)
!
!                             BLUE scaling constants
!
      real (dp), parameter :: TSML = real(RADIX(0._dp),dp)**            &
     &           CEILING((MINEXPONENT(0.0_dp)-1)*0.5_dp)
      real (dp), parameter :: TBIG = real(RADIX(0._dp),dp)**            &
     &           FLOOR((MAXEXPONENT(0.0_dp)-DIGITS(0._dp)+1)*0.5_dp)
      real (dp), parameter :: SSML = real(RADIX(0._dp),dp)**            &
     &           (-FLOOR((MINEXPONENT(0._dp)-DIGITS(0._dp))*0.5_dp))
      real (dp), parameter :: SBIG = real(RADIX(0._dp),dp)**            &
     &           (-CEILING((MAXEXPONENT(0._dp)+DIGITS(0._dp)-1)*0.5_dp))
!
      logical :: NOTBIG
      integer :: I, IX
!
      real (dp) :: ABIG, AMED, ASML, AX, SCL, SUMSQ, YMAX, YMIN
      real (dp) :: NRM2
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      NRM2 = ZERO
      IF (N <= 0) RETURN
!
      SCL = ONE
      SUMSQ = ZERO
!
!  Compute the sum of squares in 3 accumulators:
!     ABIG -- sums of squares scaled down to avoid overflow
!     ASML -- sums of squares scaled up to avoid underflow
!     AMED -- sums of squares that do not require scaling
!  The thresholds and multipliers are
!     TBIG -- values bigger than this are scaled down by SBIG
!     TSML -- values smaller than this are scaled up by SSML
!
      NOTBIG = .TRUE.
      ASML = ZERO
      AMED = ZERO
      ABIG = ZERO
      IX = 1
      IF (INCX < 0) IX = 1 - (N-1)*INCX

      DO I = 1, N
        AX = ABS(X(IX))
        IF (AX > TBIG) THEN
          ABIG = ABIG + (AX*SBIG)**2
          NOTBIG = .FALSE.
        ELSE IF (AX < TSML) THEN
          IF (NOTBIG) ASML = ASML + (AX*SSML)**2
        ELSE
          AMED = AMED + AX**2
        END IF
        IX = IX + INCX
      END DO
!
!  Combine ABIG and AMED or AMED and ASML if more than one
!  accumulator was used.
!
      IF (ABIG > ZERO) THEN
!
!  Combine ABIG and AMED if ABIG > 0.
!
        IF ((AMED > ZERO) .OR. (AMED > MAXN) .OR. (AMED /= AMED)) THEN
          ABIG = ABIG + (AMED*SBIG)*SBIG
        END IF
        SCL = ONE / SBIG
        SUMSQ = ABIG

      ELSE IF (ASML > ZERO) THEN
!
!  Combine AMED and ASML if ASML > 0.
!
        IF ((AMED > ZERO) .OR. (AMED > MAXN) .OR. (AMED /= AMED)) THEN
          AMED = SQRT(AMED)
          ASML = SQRT(ASML) / SSML
          IF (ASML > AMED) THEN
            YMIN = AMED
            YMAX = ASML
          ELSE
            YMIN = ASML
            YMAX = AMED
          END IF
          SCL = ONE
          SUMSQ = YMAX**2*(ONE + (YMIN/YMAX)**2)
        ELSE
          SCL = ONE / SSML
          SUMSQ = ASML
        END IF

      ELSE
!
!  Otherwise, all values are mid-range.
!
        SCL = ONE
        SUMSQ = AMED
      END IF
      NRM2 = SCL*SQRT(SUMSQ)
!
      RETURN
      END FUNCTION DNRM2
!
      FUNCTION IDAMAX (N, DX, INCX) RESULT (IAMAX)
!
!=======================================================================
!                                                                      !
!  IDAMAX finds the index of the first element having maximum absolute !
!  value.                                                              !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           number of elements in input vector(s)                      !
!                                                                      !
!  DX       [in]      REAL array, dimension (1 + (N-1) * ABS(INCX))    !
!                                                                      !
!  INCX     [in]      INTEGER                                          !
!           storage spacing between elements of DX                     !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!      jack dongarra, linpack, 3/11/78.                                !
!      modified 3/93 to return if incx .le. 0.                         !
!      modified 12/3/93, array(1) declarations changed to array(*)     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in) :: INCX, N
!
      real (dp), intent(in) :: DX(*)
!
!  Local variable declarations.
!
      integer   :: I, IX
!
      real (dp) :: DMAX, IAMAX
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IAMAX = 1
      IF (N.EQ.1) RETURN

      IF (INCX.EQ.1) THEN
!
!  Code for increment equal to 1.
!
        DMAX = ABS(DX(1))
        DO I = 2,N
          IF (ABS(DX(I)).GT.DMAX) THEN
            IAMAX = I
            DMAX = ABS(DX(I))
          END IF
        END DO
      ELSE
!
!  Code for increment not equal to 1.
!
        IX = 1
        DMAX = ABS(DX(1))
        IX = IX + INCX
        DO I = 2,N
          IF (ABS(DX(IX)).GT.DMAX) THEN
            IAMAX = I
            DMAX = ABS(DX(IX))
          END IF
          IX = IX + INCX
        END DO
      END IF
!
      RETURN
      END FUNCTION IDAMAX
!
      FUNCTION LSAME (CA, CB) RESULT (IsSAME)
!
!=======================================================================
!                                                                      !
!  LSAME returns .TRUE. if CA is the same letter as CB regardless of   !
!  case.                                                               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  CA      (input) CHARACTER*1                                         !
!  CB      (input) CHARACTER*1                                         !
!          CA and CB specify the single characters to be compared.     !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     September 30, 1994                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), intent(in) :: CA, CB
!
!  Local variable declarations.
!
      logical :: IsSAME
!
      integer :: INTA, INTB, ZCODE
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test if the characters are equal.
!
      IsSAME = CA .EQ. CB
      IF (IsSAME) RETURN
!
!  Now test for equivalence if both characters are alphabetic.
!
      ZCODE = ICHAR('Z')
!
!  Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!  machines, on which ICHAR returns a value with bit 8 set.
!  ICHAR('A') on Prime machines returns 193 which is the same as
!  ICHAR('A') on an EBCDIC machine.
!
      INTA = ICHAR(CA)
      INTB = ICHAR(CB)
!
      IF (ZCODE.EQ.90 .OR. ZCODE.EQ.122) THEN
!
!  ASCII is assumed - ZCODE is the ASCII code of either lower or
!  upper case 'Z'.
!
        IF (INTA.GE.97 .AND. INTA.LE.122) INTA = INTA - 32
        IF (INTB.GE.97 .AND. INTB.LE.122) INTB = INTB - 32
!
      ELSE IF (ZCODE.EQ.233 .OR. ZCODE.EQ.169) THEN
!
!  EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!  upper case 'Z'.
!
        IF (INTA.GE.129 .AND. INTA.LE.137 .OR.                          &
     &      INTA.GE.145 .AND. INTA.LE.153 .OR.                          &
     &      INTA.GE.162 .AND. INTA.LE.169) INTA = INTA + 64
        IF (INTB.GE.129 .AND. INTB.LE.137 .OR.                          &
     &      INTB.GE.145 .AND. INTB.LE.153 .OR.                          &
     &      INTB.GE.162 .AND. INTB.LE.169) INTB = INTB + 64
!
      ELSE IF (ZCODE.EQ.218 .OR. ZCODE.EQ.250) THEN
!
!  ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!  plus 128 of either lower or upper case 'Z'.
!
        IF (INTA.GE.225 .AND. INTA.LE.250) INTA = INTA - 32
        IF (INTB.GE.225 .AND. INTB.LE.250) INTB = INTB - 32
      END IF
!
      IsSAME = INTA .EQ. INTB
!
      RETURN
      END FUNCTION LSAME
!
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!     LAPACK Library Support Routines:
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!
      FUNCTION DISNAN (DIN) RESULT (LISNAN)
!
!=======================================================================
!                                                                      !
!  DISNAN returns .TRUE. if its argument is NaN, and .FALSE.           !
!  otherwise.  To be replaced by the Fortran 2003 intrinsic in the     !
!  future.                                                             !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  DIN      [in]      REAL                                             !
!           Input to test for NaN.                                     !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in) :: DIN
!
!  Local variable declarations.
!
      logical :: LISNAN
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      LISNAN = DLAISNAN(DIN,DIN)
!
      RETURN
      END FUNCTION DISNAN
!
      SUBROUTINE DLAEBZ (IJOB, NITMAX, N, MMAX, MINP, NBMIN, ABSTOL,    &
     &                   RELTOL, PIVMIN, D, E, E2, NVAL, AB, C, MOUT,   &
     &                   NAB, WORK, IWORK, INFO)
!
!=======================================================================
!                                                                      !
!  DLAEBZ contains the iteration loops which compute and use the       !
!  function N(w), which is the count of eigenvalues of a symmetric     !
!  tridiagonal matrix T less than or equal to its argument  w.  It     !
!  performs a choice of two types of loops:                            !
!                                                                      !
!  IJOB=1, followed by                                                 !
!  IJOB=2: It takes as input a list of intervals and returns a list of !
!          sufficiently small intervals whose union contains the same  !
!          eigenvalues as the union of the original intervals.         !
!          The input intervals are (AB(j,1),AB(j,2)], j=1,...,MINP.    !
!          The output interval (AB(j,1),AB(j,2)] will contain          !
!          eigenvalues NAB(j,1)+1,...,NAB(j,2), where 1 <= j <= MOUT.  !
!                                                                      !
!  IJOB=3: It performs a binary search in each input interval          !
!          (AB(j,1),AB(j,2)] for a point  w(j)  such that              !
!          N(w(j))=NVAL(j), and uses  C(j)  as the starting point of   !
!          the search.  If such a w(j) is found, then on output        !
!          AB(j,1)=AB(j,2)=w.  If no such w(j) is found, then on       !
!          output (AB(j,1),AB(j,2)] will be a small interval           !
!          containing the point where N(w) jumps through NVAL(j),      !
!          unless that point lies outside the initial interval.        !
!                                                                      !
!  Note that the intervals are in all cases half-open intervals,       !
!  i.e., of the form  (a,b] , which includes  b  but not  a .          !
!                                                                      !
!  To avoid underflow, the matrix should be scaled so that its largest !
!  element is no greater than  overflow**(1/2) * underflow**(1/4)      !
!  in absolute value.  To assure the most accurate computation         !
!  of small eigenvalues, the matrix should be scaled to be             !
!  not much smaller than that, either.                                 !
!                                                                      !
!  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal       !
!  Matrix", Report CS41, Computer Science Dept., Stanford              !
!  University, July 21, 1966                                           !
!                                                                      !
!  Note: the arguments are, in general, *not* checked for unreasonable !
!  values.                                                             !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  IJOB     [in]      INTEGER                                          !
!           Specifies what is to be done:                              !
!           = 1:  Compute NAB for the initial intervals.               !
!           = 2:  Perform bisection iteration to find eigenvalues of T.!
!           = 3:  Perform bisection iteration to invert N(w), i.e.,    !
!                 to find a point which has a specified number of      !
!                 eigenvalues of T to its left.                        !
!           Other values will cause DLAEBZ to return with INFO=-1.     !
!                                                                      !
!  NITMAX   [in]      INTEGER                                          !
!           The maximum number of "levels" of bisection to be          !
!           performed, i.e., an interval of width W will not be made   !
!           smaller than 2^(-NITMAX) * W.  If not all intervals        !
!           have converged after NITMAX iterations, then INFO is set   !
!           to the number of non-converged intervals.                  !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The dimension n of the tridiagonal matrix T. It must be at !
!           least 1.                                                   !
!                                                                      !
!  MMAX     [in]      INTEGER                                          !
!           The maximum number of intervals.  If more than MMAX        !
!           intervals are generated, then DLAEBZ will quit with        !
!           INFO=MMAX+1.                                               !
!                                                                      !
!  MINP     [in]      INTEGER                                          !
!           The initial number of intervals.  It may not be greater    !
!           than MMAX.                                                 !
!                                                                      !
!  NBMIN    [in]      INTEGER                                          !
!           The smallest number of intervals that should be processed  !
!           using a vector loop.  If zero, then only the scalar loop   !
!           will be used.                                              !
!                                                                      !
!  ABSTOL   [in]      REAL                                             !
!           The minimum (absolute) width of an interval.  When an      !
!           interval is narrower than ABSTOL, or than RELTOL times the !
!           larger (in magnitude) endpoint, then it is considered to   !
!           be sufficiently small, i.e., converged.  This must be at   !
!           least zero.                                                !
!                                                                      !
!  RELTOL   [in]      REAL                                             !
!           The minimum relative width of an interval.  When an        !
!           interval is narrower than ABSTOL, or than RELTOL times the !
!           larger (in magnitude) endpoint, then it is considered to   !
!           be sufficiently small, i.e., converged.  Note: this should !
!           always be at least radix*machine epsilon.                  !
!                                                                      !
!  PIVMIN   [in]      REAL                                             !
!           The minimum absolute value of a "pivot" in the Sturm       !
!           sequence loop.                                             !
!           This must be at least  max |e(j)**2|*safe_min  and at      !
!           least safe_min, where safe_min is at least                 !
!           the smallest number that can divide one without overflow.  !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           The diagonal elements of the tridiagonal matrix T.         !
!                                                                      !
!  E        [in]      REAL array, dimension (N)                        !
!           The offdiagonal elements of the tridiagonal matrix T in    !
!           positions 1 through N-1.  E(N) is arbitrary.               !
!                                                                      !
!  E2       [in]      REAL array, dimension (N)                        !
!           The squares of the offdiagonal elements of the tridiagonal !
!           matrix T.  E2(N) is ignored.                               !
!                                                                      !
!  NVAL     [in,out]  INTEGER array, dimension (MINP)                  !
!           If IJOB=1 or 2, not referenced.                            !
!           If IJOB=3, the desired values of N(w).  The elements of    !
!           NVAL will be reordered to correspond with the intervals    !
!           in AB. Thus, NVAL(j) on output will not, in general be the !
!           same as NVAL(j) on input, but it will correspond with the  !
!           interval (AB(j,1),AB(j,2)] on output.                      !
!                                                                      !
!  AB       [in,out]  REAL array, dimension (MMAX,2)                   !
!           The endpoints of the intervals.  AB(j,1) is  a(j), the     !
!           left endpoint of the j-th interval, and AB(j,2) is b(j),   !
!           the right endpoint of the j-th interval.  The input        !
!           intervals will, in general, be modified, split, and        !
!           reordered by the calculation.                              !
!                                                                      !
!  C        [in,out]  REAL array, dimension (MMAX)                     !
!           If IJOB=1, ignored.                                        !
!           If IJOB=2, workspace.                                      !
!           If IJOB=3, then on input C(j) should be initialized to the !
!           first search point in the binary search.                   !
!                                                                      !
!  MOUT     [out]     INTEGER                                          !
!           If IJOB=1, the number of eigenvalues in the intervals.     !
!           If IJOB=2 or 3, the number of intervals output.            !
!           If IJOB=3, MOUT will equal MINP.                           !
!                                                                      !
!  NAB      [in,out]  INTEGER array, dimension (MMAX,2)                !
!           If IJOB=1, then on output NAB(i,j) will be set to          !
!                      N(AB(i,j)).                                     !
!           If IJOB=2, then on input, NAB(i,j) should be set. It must  !
!                      satisfy the condition:                          !
!              N(AB(i,1)) <= NAB(i,1) <= NAB(i,2) <= N(AB(i,2)),       !
!              which means that in interval i only eigenvalues         !
!              NAB(i,1)+1,...,NAB(i,2) will be considered.  Usually,   !
!              NAB(i,j)=N(AB(i,j)), from a previous call to DLAEBZ     !
!              withIJOB=1.                                             !
!              On output, NAB(i,j) will contain                        !
!              MAX(na(k),MIN(nb(k),N(AB(i,j)))), where k is the index  !
!              of the input interval that the output interval          !
!              (AB(j,1),AB(j,2)] came from, and na(k) and nb(k) are    !
!              the input values of NAB(k,1) and NAB(k,2).              !
!           If IJOB=3, then on output, NAB(i,j) contains N(AB(i,j)),   !
!              unless N(w) > NVAL(i) for all search points  w , in     !
!              which case NAB(i,1) will not be modified, i.e., the     !
!              output value will be the same as the input value        !
!              (modulo reorderings -- see NVAL and AB), or unless      !
!              N(w) < NVAL(i) for all search points  w , in which      !
!              case NAB(i,2) will not be modified.  Normally, NAB      !
!              should be set to some distinctive value(s) before       !
!              DLAEBZ is called.                                       !
!                                                                      !
!  WORK     [out]     REAL array, dimension (MMAX)                     !
!           Workspace.                                                 !
!                                                                      !
!  IWORK    [out]     INTEGER array, dimension (MMAX)                  !
!           Workspace.                                                 !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           = 0:       All intervals converged.                        !
!           = 1--MMAX: The last INFO intervals did not converge.       !
!           = MMAX+1:  More than MMAX intervals were generated.        !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!   This routine is intended to be called only by other LAPACK         !
!   routines, thus the interface is less user-friendly.  It is         !
!   intended for two purposes:                                         !
!                                                                      !
!   (a) finding eigenvalues.  In this case, DLAEBZ should have one or  !
!       more initial intervals set up in AB, and DLAEBZ should be      !
!       called with IJOB=1.  This sets up NAB, and also counts the     !
!       eigenvalues. Intervals with no eigenvalues would usually be    !
!       thrown out at this point.  Also, if not all the eigenvalues    !
!       in an interval i are desired, NAB(i,1) can be increased or     !
!       NAB(i,2) decreased.                                            !
!       For example, set NAB(i,1)=NAB(i,2)-1 to get the largest        !
!       eigenvalue.  DLAEBZ is then called with IJOB=2 and MMAX        !
!       no smaller than the value of MOUT returned by the call with    !
!       IJOB=1.  After this (IJOB=2) call, eigenvalues NAB(i,1)+1      !
!       through NAB(i,2) are approximately AB(i,1) (or AB(i,2)) to the !
!       tolerance specified by ABSTOL and RELTOL.                      !
!                                                                      !
!   (b) finding an interval (a',b'] containing eigenvalues             !
!       w(f),...,w(l). In this case, start with a Gershgorin interval  !
!       (a,b). Set up AB to contain 2 search intervals, both initially !
!       (a,b). One NVAL element should contain  f-1  and the other     !
!       should contain  l, while C should contain a and b, resp.       !
!       NAB(i,1) should be -1 and NAB(i,2) should be N+1, to flag an   !
!       error if the desired interval does not lie in (a,b).  DLAEBZ   !
!       is then called with IJOB=3.  On exit, if w(f-1) < w(f), then   !
!       one of the intervals --j -- will have AB(j,1)=AB(j,2) and      !
!       NAB(j,1)=NAB(j,2)=f-1, while if, to the specified tolerance,   !
!       w(f-k)=...=w(f+r), k > 0 and r >= 0, then the interval will    !
!       have  N(AB(j,1))=NAB(j,1)=f-k and N(AB(j,2))=NAB(j,2)=f+r.     !
!       The cases w(l) < w(l+1) and w(l-r)=...=w(l+k) are handled      !
!       similarly.                                                     !
!                                                                      !
!  -- LAPACK auxiliary routine                                         !
!  -- LAPACK is a software package provided by Univ. of Tennessee,     !
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and       !
!  -- NAG Ltd.                                                         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: IJOB, MINP, MMAX, N, NBMIN, NITMAX
      integer,   intent(inout) :: NAB(MMAX,*), NVAL(*)
      integer,   intent(out)   :: INFO, MOUT
      integer,   intent(out)   :: IWORK(*)
!
      real (dp), intent(in)    :: ABSTOL, PIVMIN, RELTOL
      real (dp), intent(in)    :: D(*), E(*), E2(*)
      real (dp), intent(inout) :: AB(MMAX,*), C(*)
      real (dp), intent(out)   :: WORK(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: HALF = 0.5_dp
!
      integer   :: ITMP1, ITMP2, J, JI, JIT, JP, KF, KFNEW, KL, KLNEW
!
      real (dp) :: TMP1, TMP2
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Check for Errors.
!
      INFO = 0
      IF (IJOB.LT.1 .OR. IJOB.GT.3) THEN
         INFO = -1
         RETURN
      END IF
!
!  Initialize NAB
!
      IF (IJOB.EQ.1) THEN
!
!  Compute the number of eigenvalues in the initial intervals.
!
         MOUT = 0
         DO JI = 1, MINP
            DO JP = 1, 2
               TMP1 = D(1) - AB(JI, JP)
               IF (ABS(TMP1).LT.PIVMIN) TMP1 = -PIVMIN
               NAB(JI, JP) = 0
               IF (TMP1.LE.ZERO) NAB(JI, JP) = 1
!
               DO J = 2, N
                  TMP1 = D(J) - E2(J-1) / TMP1 - AB(JI,JP)
                  IF (ABS(TMP1).LT.PIVMIN) TMP1 = -PIVMIN
                  IF (TMP1.LE.ZERO) NAB(JI,JP) = NAB(JI,JP) + 1
               END DO
            END DO
            MOUT = MOUT + NAB(JI, 2) - NAB(JI, 1)
         END DO
         RETURN
      END IF
!
!  Initialize for loop.
!
!     KF and KL have the following meaning:
!        Intervals 1,...,KF-1 have converged.
!        Intervals KF,...,KL  still need to be refined.
!
      KF = 1
      KL = MINP
!
!  If IJOB=2, initialize C.
!  If IJOB=3, use the user-supplied starting point.
!
      IF (IJOB.EQ.2) THEN
         DO JI = 1, MINP
            C(JI) = HALF*(AB(JI, 1)+AB(JI, 2))
         END DO
      END IF
!
!  Iteration loop.
!
      DO 130 JIT = 1, NITMAX
!
!  Loop over intervals.
!
         IF (KL-KF+1.GE.NBMIN .AND. NBMIN.GT.0) THEN
!
!  Begin of Parallel Version of the loop.
!
            DO JI = KF, KL
!
!  Compute N(c), the number of eigenvalues less than c.
!
               WORK(JI) = D(1) - C(JI)
               IWORK(JI) = 0
               IF (WORK(JI).LE.PIVMIN) THEN
                  IWORK(JI) = 1
                  WORK(JI) = MIN(WORK(JI), -PIVMIN)
               END IF
!
               DO J = 2, N
                  WORK(JI) = D(J) - E2(J-1) / WORK(JI) - C(JI)
                  IF (WORK(JI).LE.PIVMIN) THEN
                     IWORK(JI) = IWORK(JI) + 1
                     WORK(JI) = MIN(WORK(JI), -PIVMIN)
                  END IF
               END DO
            END DO
!
            IF (IJOB.LE.2) THEN
!
!  IJOB=2: Choose all intervals containing eigenvalues.
!
               KLNEW = KL
               DO JI = KF, KL
!
!  Insure that N(w) is monotone.
!
                  IWORK(JI) = MIN(NAB(JI,2),                          &
     &                              MAX(NAB(JI,1), IWORK(JI)))
!
!  Update the Queue -- add intervals if both halves
!  contain eigenvalues.
!
                  IF (IWORK(JI).EQ.NAB(JI, 2)) THEN
!
!  No eigenvalue in the upper interval: just use the lower interval.
!
                     AB(JI, 2) = C(JI)
!
                  ELSE IF (IWORK(JI).EQ.NAB(JI, 1)) THEN
!
!  No eigenvalue in the lower interval: just use the upper interval.
!
                     AB(JI, 1) = C(JI)
                  ELSE
                     KLNEW = KLNEW + 1
                     IF (KLNEW.LE.MMAX) THEN
!
!  Eigenvalue in both intervals -- add upper to queue.
!
                        AB(KLNEW, 2) = AB(JI, 2)
                        NAB(KLNEW, 2) = NAB(JI, 2)
                        AB(KLNEW, 1) = C(JI)
                        NAB(KLNEW, 1) = IWORK(JI)
                        AB(JI, 2) = C(JI)
                        NAB(JI, 2) = IWORK(JI)
                     ELSE
                        INFO = MMAX + 1
                     END IF
                  END IF
               END DO
               IF (INFO.NE.0) RETURN
               KL = KLNEW
            ELSE
!
!  IJOB=3: Binary search.  Keep only the interval containing
!                          w   s.t. N(w) = NVAL
!
               DO JI = KF, KL
                  IF (IWORK(JI).LE.NVAL(JI)) THEN
                     AB(JI, 1) = C(JI)
                     NAB(JI, 1) = IWORK(JI)
                  END IF
                  IF (IWORK(JI).GE.NVAL(JI)) THEN
                     AB(JI, 2) = C(JI)
                     NAB(JI, 2) = IWORK(JI)
                  END IF
               END DO
            END IF
!
         ELSE
!
!  End of Parallel Version of the loop.
!
!  Begin of Serial Version of the loop.
!
            KLNEW = KL
            DO JI = KF, KL
!
!  Compute N(w), the number of eigenvalues less than w.
!
               TMP1 = C(JI)
               TMP2 = D(1) - TMP1
               ITMP1 = 0
               IF (TMP2.LE.PIVMIN) THEN
                  ITMP1 = 1
                  TMP2 = MIN(TMP2, -PIVMIN)
               END IF
!
               DO J = 2, N
                  TMP2 = D(J) - E2(J-1) / TMP2 - TMP1
                  IF (TMP2.LE.PIVMIN) THEN
                     ITMP1 = ITMP1 + 1
                     TMP2 = MIN(TMP2, -PIVMIN)
                  END IF
               END DO
!
               IF (IJOB.LE.2) THEN
!
!  IJOB=2: Choose all intervals containing eigenvalues.
!
!  Insure that N(w) is monotone.
!
                  ITMP1 = MIN(NAB(JI,2),                                &
     &                        MAX(NAB(JI,1), ITMP1))
!
!  Update the Queue -- add intervals if both halves
!  contain eigenvalues.
!
                  IF (ITMP1.EQ.NAB(JI,2)) THEN
!
!  No eigenvalue in the upper interval: just use the lower interval.
!
                     AB(JI, 2) = TMP1
!
                  ELSE IF (ITMP1.EQ.NAB(JI,1)) THEN
!
!  No eigenvalue in the lower interval: just use the upper interval.
!
                     AB(JI, 1) = TMP1
                  ELSE IF (KLNEW.LT.MMAX) THEN
!
!  Eigenvalue in both intervals -- add upper to queue.
!
                     KLNEW = KLNEW + 1
                     AB(KLNEW, 2) = AB(JI, 2)
                     NAB(KLNEW, 2) = NAB(JI, 2)
                     AB(KLNEW, 1) = TMP1
                     NAB(KLNEW, 1) = ITMP1
                     AB(JI, 2) = TMP1
                     NAB(JI, 2) = ITMP1
                  ELSE
                     INFO = MMAX + 1
                     RETURN
                  END IF
               ELSE
!
!  IJOB=3: Binary search.  Keep only the interval
!                          containing  w  s.t. N(w) = NVAL
!
                  IF (ITMP1.LE.NVAL(JI)) THEN
                     AB(JI, 1) = TMP1
                     NAB(JI, 1) = ITMP1
                  END IF
                  IF (ITMP1.GE.NVAL(JI)) THEN
                     AB(JI, 2) = TMP1
                     NAB(JI, 2) = ITMP1
                  END IF
               END IF
            END DO
            KL = KLNEW
!
         END IF
!
!  Check for convergence.
!
         KFNEW = KF
         DO JI = KF, KL
            TMP1 = ABS(AB(JI, 2)-AB(JI, 1))
            TMP2 = MAX(ABS(AB(JI, 2)), ABS(AB(JI, 1)))
            IF (TMP1.LT.MAX(ABSTOL, PIVMIN, RELTOL*TMP2) .OR.           &
     &          NAB(JI, 1).GE.NAB(JI, 2)) THEN
!
!  Converged -- Swap with position KFNEW, then increment KFNEW.
!
               IF (JI.GT.KFNEW) THEN
                  TMP1 = AB(JI, 1)
                  TMP2 = AB(JI, 2)
                  ITMP1 = NAB(JI, 1)
                  ITMP2 = NAB(JI, 2)
                  AB(JI, 1) = AB(KFNEW, 1)
                  AB(JI, 2) = AB(KFNEW, 2)
                  NAB(JI, 1) = NAB(KFNEW, 1)
                  NAB(JI, 2) = NAB(KFNEW, 2)
                  AB(KFNEW, 1) = TMP1
                  AB(KFNEW, 2) = TMP2
                  NAB(KFNEW, 1) = ITMP1
                  NAB(KFNEW, 2) = ITMP2
                  IF (IJOB.EQ.3) THEN
                     ITMP1 = NVAL(JI)
                     NVAL(JI) = NVAL(KFNEW)
                     NVAL(KFNEW) = ITMP1
                  END IF
               END IF
               KFNEW = KFNEW + 1
            END IF
         END DO
         KF = KFNEW
!
!  Choose Midpoints.
!
         DO JI = KF, KL
            C(JI) = HALF*(AB(JI, 1)+AB(JI, 2))
         END DO
!
!  If no more intervals to refine, quit.
!
         IF (KF.GT.KL) GOTO 140
  130 CONTINUE
!
!  Converged.
!
  140 CONTINUE
      INFO = MAX(KL+1-KF, 0)
      MOUT = KL
!
      RETURN
      END SUBROUTINE DLAEBZ
!
      SUBROUTINE DLAE2 (A, B, C, RT1, RT2)
!
!=======================================================================
!                                                                      !
!  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix:       !
!                                                                      !
!     [  A   B  ]                                                      !
!     [  B   C  ].                                                     !
!                                                                      !
!  On return, RT1 is the eigenvalue of larger absolute value, and RT2  !
!  is the eigenvalue of smaller absolute value.                        !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  A       (input) DOUBLE PRECISION                                    !
!          The (1,1) element of the 2-by-2 matrix.                     !
!                                                                      !
!  B       (input) DOUBLE PRECISION                                    !
!          The (1,2) and (2,1) elements of the 2-by-2 matrix.          !
!                                                                      !
!  C       (input) DOUBLE PRECISION                                    !
!          The (2,2) element of the 2-by-2 matrix.                     !
!                                                                      !
!  RT1     (output) DOUBLE PRECISION                                   !
!          The eigenvalue of larger absolute value.                    !
!                                                                      !
!  RT2     (output) DOUBLE PRECISION                                   !
!          The eigenvalue of smaller absolute value.                   !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  RT1 is accurate to a few ulps barring over/underflow.               !
!                                                                      !
!  RT2 may be inaccurate if there is massive cancellation in the       !
!  determinant A*C-B*B; higher precision or correctly rounded or       !
!  correctly truncated arithmetic would be needed to compute RT2       !
!  accurately in all cases.                                            !
!                                                                      !
!  Overflow is possible only if RT1 is within a factor of 5 of         !
!  overflow.                                                           !
!                                                                      !
!  Underflow is harmless if the input data is 0 or exceeds             !
!  underflow_threshold / macheps.                                      !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in)  :: A, B, C
      real (dp), intent(out) :: RT1, RT2
!
!  Local variable declarations.
!
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: HALF = 0.5_dp
!
      real (dp) :: AB, ACMN, ACMX, ADF, DF, RT, SM, TB
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Compute the eigenvalues.
!
      SM = A + C
      DF = A - C
      ADF = ABS(DF)
      TB = B + B
      AB = ABS(TB)
!
      IF (ABS(A) .GT. ABS(C)) THEN
        ACMX = A
        ACMN = C
      ELSE
        ACMX = C
        ACMN = A
      END IF
!
      IF (ADF .GT. AB) THEN
        RT = ADF*SQRT(ONE+(AB/ADF)**2)
      ELSE IF (ADF .LT. AB) THEN
        RT = AB*SQRT(ONE+(ADF/AB)**2)
      ELSE
        RT = AB*SQRT(TWO)                 ! Includes case AB=ADF=0
      END IF
!
      IF (SM .LT. ZERO) THEN
        RT1 = HALF*(SM-RT)
!
!  Order of execution important.
!  To get fully accurate smaller eigenvalue,
!  next line needs to be executed in higher precision.
!
        RT2 = (ACMX / RT1)*ACMN - (B / RT1)*B
      ELSE IF (SM .GT. ZERO) THEN
        RT1 = HALF*(SM + RT)
!
!  Order of execution important.
!  To get fully accurate smaller eigenvalue,
!  next line needs to be executed in higher precision.
!
        RT2 = (ACMX / RT1)*ACMN - (B / RT1)*B
      ELSE
        RT1 = HALF*RT                     ! Includes case RT1 = RT2 = 0
        RT2 = -HALF*RT
      END IF
!
      RETURN
      END SUBROUTINE DLAE2
!
      SUBROUTINE DLAEV2 (A, B, C, RT1, RT2, CS1, SN1)
!
!=======================================================================
!                                                                      !
!  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric        !
!  matrix:                                                             !
!                                                                      !
!     [  A   B  ]                                                      !
!     [  B   C  ].                                                     !
!                                                                      !
!  On return, RT1 is the eigenvalue of larger absolute value, RT2 is   !
!  the eigenvalue of smaller absolute value, and (CS1,SN1) is the unit !
!  right eigenvector for RT1, giving the decomposition:                !
!                                                                      !
!     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]            !
!     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].           !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  A       (input) DOUBLE PRECISION                                    !
!          The (1,1) element of the 2-by-2 matrix.                     !
!                                                                      !
!  B       (input) DOUBLE PRECISION                                    !
!          The (1,2) element and the conjugate of the (2,1) element of !
!          the 2-by-2 matrix.                                          !
!                                                                      !
!  C       (input) DOUBLE PRECISION                                    !
!          The (2,2) element of the 2-by-2 matrix.                     !
!                                                                      !
!  RT1     (output) DOUBLE PRECISION                                   !
!          The eigenvalue of larger absolute value.                    !
!                                                                      !
!  RT2     (output) DOUBLE PRECISION                                   !
!          The eigenvalue of smaller absolute value.                   !
!                                                                      !
!  CS1     (output) DOUBLE PRECISION                                   !
!  SN1     (output) DOUBLE PRECISION                                   !
!          The vector (CS1, SN1) is a unit right eigenvector for RT1.  !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  RT1 is accurate to a few ulps barring over/underflow.               !
!                                                                      !
!  RT2 may be inaccurate if there is massive cancellation in the       !
!  determinant A*C-B*B; higher precision or correctly rounded or       !
!  correctly truncated arithmetic would be needed to compute RT2       !
!  accurately in all cases.                                            !
!                                                                      !
!  CS1 and SN1 are accurate to a few ulps barring over/underflow.      !
!                                                                      !
!  Overflow is possible only if RT1 is within a factor of 5 of         !
!  overflow.                                                           !
!                                                                      !
!  Underflow is harmless if the input data is 0 or exceeds             !
!  underflow_threshold / macheps.                                      !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in)  :: A, B, C
      real (dp), intent(out) :: CS1, RT1, RT2, SN1
!
!  Local variable declarations.
!
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: HALF = 0.5_dp
!
      integer :: SGN1, SGN2
      real (dp) :: AB, ACMN, ACMX, ACS, ADF, CS, CT, DF, RT, SM, TB, TN
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Compute the eigenvalues.
!
      SM = A + C
      DF = A - C
      ADF = ABS(DF)
      TB = B + B
      AB = ABS(TB)
!
      IF (ABS(A) .GT. ABS(C)) THEN
         ACMX = A
         ACMN = C
      ELSE
         ACMX = C
         ACMN = A
      END IF
!
      IF (ADF .GT. AB) THEN
        RT = ADF*SQRT(ONE+(AB/ADF)**2)
      ELSE IF (ADF .LT. AB) THEN
        RT = AB*SQRT(ONE+(ADF/AB)**2)
      ELSE
        RT = AB*SQRT(TWO)                 ! Includes case AB=ADF=0
      END IF
      IF (SM.LT.ZERO) THEN
        RT1 = HALF*(SM - RT)
        SGN1 = -1
!
!  Order of execution important.
!  To get fully accurate smaller eigenvalue,
!  next line needs to be executed in higher precision.
!
        RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE IF (SM .GT. ZERO) THEN
        RT1 = HALF*(SM + RT)
        SGN1 = 1
!
!  Order of execution important.
!  To get fully accurate smaller eigenvalue,
!  next line needs to be executed in higher precision.
!
        RT2 = (ACMX/RT1)*ACMN - (B/RT1)*B
      ELSE
        RT1 = HALF*RT                     ! Includes case RT1 = RT2 = 0
        RT2 = -HALF*RT
        SGN1 = 1
      END IF
!
!  Compute the eigenvector.
!
      IF (DF .GE. ZERO) THEN
        CS = DF + RT
        SGN2 = 1
      ELSE
        CS = DF - RT
        SGN2 = -1
      END IF
      ACS = ABS(CS)
      IF (ACS .GT. AB) THEN
        CT = -TB / CS
        SN1 = ONE / SQRT(ONE+CT*CT)
        CS1 = CT*SN1
      ELSE
        IF (AB .EQ. ZERO) THEN
          CS1 = ONE
          SN1 = ZERO
        ELSE
          TN = -CS / TB
          CS1 = ONE / SQRT(ONE+TN*TN)
          SN1 = TN*CS1
        END IF
      END IF
      IF (SGN1 .EQ. SGN2) THEN
        TN = CS1
        CS1 = -SN1
        SN1 = TN
      END IF
!
      RETURN
      END SUBROUTINE DLAEV2
!
      SUBROUTINE DLAGTF (N, A, LAMBDA, B, C, TOL, D, IN, INFO)
!
!=======================================================================
!                                                                      !
!  DLAGTF factorizes the matrix (T - lambda*I), where T is an N by H   !
!  tridiagonal matrix and lambda is a scalar, as                       !
!                                                                      !
!     T - lambda*I = PLU,                                              !
!                                                                      !
!  where P is a permutation matrix, L is a unit lower tridiagonal      !
!  matrix with at most one non-zero sub-diagonal elements per column   !
!  and U is an upper triangular matrix with at most two non-zero       !
!  super-diagonal elements per column.                                 !
!                                                                      !
!  The factorization is obtained by Gaussian elimination with partial  !
!  pivoting and implicit row scaling.                                  !
!                                                                      !
!  The parameter LAMBDA is included in the routine so that DLAGTF may  !
!  be used, in conjunction with DLAGTS, to obtain eigenvectors of T by !
!  inverse iteration.                                                  !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix T.                                 !
!                                                                      !
!  A        [in,out]  REAL array, dimension (N)                        !
!           On entry, A must contain the diagonal elements of T.       !
!           On exit, A is overwritten by the n diagonal elements of    !
!           the upper triangular matrix U of the factorization of T.   !
!                                                                      !
!  LAMBDA   [in]      REAL                                             !
!           On entry, the scalar lambda.                               !
!                                                                      !
!  B        [in,out]  REAL array, dimension (N-1)                      !
!           On entry, B must contain the (n-1) super-diagonal elements !
!           of T.                                                      !
!           On exit, B is overwritten by the (n-1) super-diagonal      !
!           elements of the matrix U of the factorization of T.        !
!                                                                      !
!  C        [in,out]  REAL array, dimension (N-1)                      !
!           On entry, C must contain the (n-1) sub-diagonal elements   !
!           of T.                                                      !
!           On exit, C is overwritten by the (n-1) sub-diagonal        !
!           elements of the matrix L of the factorization of T.        !
!                                                                      !
!  TOL      [in] TOL  REAL                                             !
!           On entry, a relative tolerance used to indicate whether    !
!           or not the matrix (T - lambda*I) is nearly singular. TOL   !
!           should normally be chose as approximately the largest      !
!           relative error in the elements of T. For example, if the   !
!           elements of T are correct to about 4 significant figures,  !
!           then TOL should be set to about 5*10**(-4). If TOL is      !
!           supplied as less than eps, where eps is the relative       !
!           machine precision, then the value eps is used in place     !
!           of TOL.                                                    !
!                                                                      !
!  D        [out]     REAL array, dimension (N-2)                      !
!           On exit, D is overwritten by the (n-2) second super-       !
!           diagonal elements of the matrix U of the factorization     !
!           of T.                                                      !
!                                                                      !
!  IN       [out]     INTEGER array, dimension (N)                     !
!           On exit, IN contains details of the permutation matrix P.  !
!           If an interchange occurred at the kth step of the          !
!           elimination, then IN(k) = 1, otherwise IN(k) = 0. The      !
!           element IN(n) returns the smallest positive integer j,     !
!           such that                                                  !
!                                                                      !
!              ABS(u(j,j)) <= norm((T - lambda*I)(j))*TOL,             !
!                                                                      !
!           where norm(A(j)) denotes the sum of the absolute values    !
!           of the jth row of the matrix A. If no such j exists then   !
!           IN(n) is returned as zero. If IN(n) is returned as         !
!           positive, then a diagonal element of U is small,           !
!           indicating that (T - lambda*I) is singular or nearly       !
!           singular.                                                  !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           = 0:  successful exit                                      !
!           < 0:  if INFO = -k, the kth argument had an illegal value  !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  -- LAPACK computational routine --                                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: N
      integer,   intent(out)   :: INFO
      integer,   intent(out)   :: IN(*)
!
      real (dp), intent(in)    :: LAMBDA, TOL
      real (dp), intent(inout) :: A(*), B(*), C(*)
      real (dp), intent(out)   :: D(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
!
      integer :: K
!
      real (dp) :: EPS, MULT, PIV1, PIV2, SCALE1, SCALE2, TEMP, TL
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
      IF (N.LT.0) THEN
         INFO = -1
         CALL XERBLA ('DLAGTF', -INFO)
         RETURN
      END IF
!
      IF (N.EQ.0) RETURN
!
      A(1) = A(1) - LAMBDA
      IN(N) = 0
      IF (N.EQ.1) THEN
         IF (A(1).EQ.ZERO) IN(1) = 1
         RETURN
      END IF
!
      EPS = DLAMCH('Epsilon')
!
      TL = MAX(TOL, EPS)
      SCALE1 = ABS(A(1)) + ABS(B(1))
      DO K = 1, N - 1
         A(K+1) = A(K+1) - LAMBDA
         SCALE2 = ABS(C(K)) + ABS(A(K+1))
         IF (K.LT.(N-1)) SCALE2 = SCALE2 + ABS(B(K+1))
         IF (A(K).EQ.ZERO) THEN
            PIV1 = ZERO
         ELSE
            PIV1 = ABS(A(K)) / SCALE1
         END IF
         IF (C(K).EQ.ZERO) THEN
            IN(K) = 0
            PIV2 = ZERO
            SCALE1 = SCALE2
            IF (K.LT.(N-1)) D(K) = ZERO
         ELSE
            PIV2 = ABS(C(K)) / SCALE2
            IF (PIV2.LE.PIV1) THEN
               IN(K) = 0
               SCALE1 = SCALE2
               C(K) = C(K) / A(K)
               A(K+1) = A(K+1) - C(K)*B(K)
               IF (K.LT.(N-1)) D(K) = ZERO
            ELSE
               IN(K) = 1
               MULT = A(K) / C(K)
               A(K) = C(K)
               TEMP = A(K+1)
               A(K+1) = B(K) - MULT*TEMP
               IF (K.LT.(N-1)) THEN
                  D(K) = B(K+1)
                  B(K+1) = -MULT*D(K)
               END IF
               B(K) = TEMP
               C(K) = MULT
            END IF
         END IF
         IF ((MAX(PIV1, PIV2).LE.TL) .AND. (IN(N).EQ.0)) IN(N) = K
      END DO
      IF ((ABS(A(N)).LE.SCALE1*TL) .AND. (IN(N).EQ.0)) IN(N) = N
!
      RETURN
      END SUBROUTINE DLAGTF
!
      SUBROUTINE DLAGTS (JOB, N, A, B, C, D, IN, Y, TOL, INFO)
!
!=======================================================================
!                                                                      !
!  DLAGTS may be used to solve one of the systems of equations         !
!                                                                      !
!     (T - lambda*I)*x = y   or   (T - lambda*I)**T*x = y,             !
!                                                                      !
!  where T is an N by N tridiagonal matrix, for x, following the       !
!  factorization of (T - lambda*I) as                                  !
!                                                                      !
!     (T - lambda*I) = P*L*U ,                                         !
!                                                                      !
!  by routine DLAGTF. The choice of equation to be solved is           !
!  controlled by the argument JOB, and in each case there is an option !
!  to perturb zero or very small diagonal elements of U, this option   !
!  being intended for use in applications such as inverse iteration.   !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  JOB      [in]      INTEGER                                          !
!           Specifies the job to be performed by DLAGTS as follows:    !
!           =  1: The equations  (T - lambda*I)x = y  are to be solved,!
!                 but diagonal elements of U are not to be perturbed.  !
!           = -1: The equations  (T - lambda*I)x = y  are to be solved !
!                 and, if overflow would otherwise occur, the diagonal !
!                 elements of U are to be perturbed. See argument TOL  !
!                 below.                                               !
!           =  2: The equations  (T - lambda*I)**Tx = y  are to be     !
!                 solved, but diagonal elements of U are not to be     !
!                 perturbed.                                           !
!           = -2: The equations  (T - lambda*I)**Tx = y  are to be     !
!                 solved and, if overflow would otherwise occur, the   !
!                 diagonal elements of U are to be perturbed. See      !
!                 argument TOL below.                                  !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix T.                                 !
!                                                                      !
!  A        [in]      REAL array, dimension (N)                        !
!           On entry, A must contain the diagonal elements of U as     !
!           returned from DLAGTF.                                      !
!                                                                      !
!  B        [in]      REAL array, dimension (N-1)                      !
!           On entry, B must contain the first super-diagonal elements !
!           of U as returned from DLAGTF.                              !
!                                                                      !
!  C        [in] C    REAL array, dimension (N-1)                      !
!           On entry, C must contain the sub-diagonal elements of L as !
!           returned from DLAGTF.                                      !
!                                                                      !
!  D        [in]      REAL array, dimension (N-2)                      !
!           On entry, D must contain the second super-diagonal         !
!           elements of U as returned from DLAGTF.                     !
!                                                                      !
!  IN       [in]      INTEGER array, dimension (N)                     !
!           On entry, IN must contain details of the matrix P as       !
!           returned from DLAGTF.                                      !
!                                                                      !
!  Y        [in,out]  REAL array, dimension (N)                        !
!           On entry, the right hand side vector y.                    !
!           On exit, Y is overwritten by the solution vector x.        !
!                                                                      !
!  TOL      [in,out]  REAL                                             !
!           On entry, with  JOB < 0, TOL should be the minimum         !
!           perturbation to be made to very small diagonal elements    !
!           of U. TOL should normally be chosen as about eps*norm(U),  !
!           where eps is the relative machine precision, but if TOL    !
!           is supplied as non-positive, then it is reset to           !
!           eps*MAX(ABS(u(i,j))). If JOB > 0 then TOL is not           !
!           referenced.                                                !
!           On exit, TOL is changed as described above, only if TOL    !
!           is non-positive on entry. Otherwise TOL is unchanged.      !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           = 0:  successful exit                                      !
!           < 0:  if INFO = -i, the i-th argument had an illegal value !
!           > 0:  overflow would occur when computing the INFO(th)     !
!                 element of the solution vector x. This can only      !
!                 occur when JOB is supplied as positive and either    !
!                 means that a diagonal element of U is very small,    !
!                 or that the elements of the right-hand side vector   !
!                 y are very large.                                    !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
       integer,   intent(in)    :: JOB, N
       integer,   intent(in)    :: IN(*)
       integer,   intent(out)   :: INFO
!
       real (dp), intent(in)    :: A(*), B(*), C(*), D(*)
       real (dp), intent(inout) :: TOL
       real (dp), intent(inout) :: Y(*)
!
!  Local variable declarations.
!
       real (dp), parameter :: ZERO = 0.0_dp
       real (dp), parameter :: ONE = 1.0_dp
!
       integer   :: K
!
       real (dp) :: ABSAK, AK, BIGNUM, EPS, PERT, SFMIN, TEMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
      IF ((ABS(JOB).GT.2) .OR. (JOB.EQ.0)) THEN
         INFO = -1
      ELSE IF (N.LT.0) THEN
         INFO = -2
      END IF
      IF (INFO.NE.0) THEN
         CALL XERBLA ('DLAGTS', -INFO)
         RETURN
      END IF
!
      IF (N.EQ.0) RETURN
!
      EPS = DLAMCH('Epsilon')
      SFMIN = DLAMCH('Safe minimum')
      BIGNUM = ONE / SFMIN
!
      IF (JOB.LT.0) THEN
         IF (TOL.LE.ZERO) THEN
            TOL = ABS(A(1))
            IF (N.GT.1) TOL = MAX(TOL,ABS(A(2)), ABS(B(1)))
            DO K = 3, N
               TOL = MAX(TOL, ABS(A(K)), ABS(B(K-1)), ABS(D(K-2)))
            END DO
            TOL = TOL*EPS
            IF (TOL.EQ.ZERO) TOL = EPS
         END IF
      END IF
!
      IF (ABS(JOB).EQ.1) THEN
         DO K = 2, N
            IF (IN(K-1).EQ.0) THEN
               Y(K) = Y(K) - C(K-1)*Y(K-1)
            ELSE
               TEMP = Y(K-1)
               Y(K-1) = Y(K)
               Y(K) = TEMP - C(K-1)*Y(K)
            END IF
         END DO

         IF (JOB.EQ.1) THEN
            DO K = N, 1, -1
               IF (K.LE.N-2) THEN
                  TEMP = Y(K) - B(K)*Y(K+1) - D(K)*Y(K+2)
               ELSE IF (K.EQ.N-1) THEN
                  TEMP = Y(K) - B(K)*Y(K+1)
               ELSE
                  TEMP = Y(K)
               END IF
               AK = A(K)
               ABSAK = ABS(AK)
               IF (ABSAK.LT.ONE) THEN
                  IF (ABSAK.LT.SFMIN) THEN
                     IF (ABSAK.EQ.ZERO .OR.                            &
     &                   ABS(TEMP)*SFMIN.GT.ABSAK) THEN
                        INFO = K
                        RETURN
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF (ABS(TEMP).GT.ABSAK*BIGNUM) THEN
                     INFO = K
                     RETURN
                  END IF
               END IF
               Y(K) = TEMP / AK
            END DO
         ELSE
            DO K = N, 1, -1
               IF (K.LE.N-2) THEN
                  TEMP = Y(K) - B(K)*Y(K+1) - D(K)*Y(K+2)
               ELSE IF (K.EQ.N-1) THEN
                  TEMP = Y(K) - B(K)*Y(K+1)
               ELSE
                  TEMP = Y(K)
               END IF
               AK = A(K)
               PERT = SIGN(TOL, AK)
   40          CONTINUE
               ABSAK = ABS(AK)
               IF (ABSAK.LT.ONE) THEN
                  IF (ABSAK.LT.SFMIN) THEN
                     IF (ABSAK.EQ.ZERO .OR.                             &
     &                   ABS(TEMP)*SFMIN.GT.ABSAK) THEN
                        AK = AK + PERT
                        PERT = 2*PERT
                        GOTO 40
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF (ABS(TEMP).GT.ABSAK*BIGNUM) THEN
                     AK = AK + PERT
                     PERT = 2*PERT
                     GOTO 40
                  END IF
               END IF
               Y(K) = TEMP / AK
            END DO
         END IF
      ELSE
!
!  Come to here if  JOB = 2 or -2.
!
         IF (JOB.EQ.2) THEN
            DO K = 1, N
               IF (K.GE.3) THEN
                  TEMP = Y(K) - B(K-1)*Y(K-1) - D(K-2)*Y(K-2)
               ELSE IF (K.EQ.2) THEN
                  TEMP = Y(K) - B(K-1)*Y(K-1)
               ELSE
                  TEMP = Y(K)
               END IF
               AK = A(K)
               ABSAK = ABS(AK)
               IF (ABSAK.LT.ONE) THEN
                  IF (ABSAK.LT.SFMIN) THEN
                     IF (ABSAK.EQ.ZERO .OR.                             &
     &                   ABS(TEMP)*SFMIN.GT.ABSAK) THEN
                        INFO = K
                        RETURN
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF (ABS(TEMP).GT.ABSAK*BIGNUM) THEN
                     INFO = K
                     RETURN
                  END IF
               END IF
               Y(K) = TEMP / AK
            END DO
         ELSE
            DO K = 1, N
               IF (K.GE.3) THEN
                  TEMP = Y(K) - B(K-1)*Y(K-1) - D(K-2)*Y(K-2)
               ELSE IF (K.EQ.2) THEN
                  TEMP = Y(K) - B(K-1)*Y(K-1)
               ELSE
                  TEMP = Y(K)
               END IF
               AK = A(K)
               PERT = SIGN(TOL, AK)
   70          CONTINUE
               ABSAK = ABS(AK)
               IF (ABSAK.LT.ONE) THEN
                  IF (ABSAK.LT.SFMIN) THEN
                     IF (ABSAK.EQ.ZERO .OR.                             &
     &                   ABS(TEMP)*SFMIN.GT.ABSAK) THEN
                        AK = AK + PERT
                        PERT = 2*PERT
                        GOTO 70
                     ELSE
                        TEMP = TEMP*BIGNUM
                        AK = AK*BIGNUM
                     END IF
                  ELSE IF (ABS(TEMP).GT.ABSAK*BIGNUM) THEN
                     AK = AK + PERT
                     PERT = 2*PERT
                     GOTO 70
                  END IF
               END IF
               Y(K) = TEMP / AK
            END DO
         END IF
!
         DO K = N, 2, -1
            IF (IN(K-1).EQ.0) THEN
               Y(K-1) = Y(K-1) - C(K-1)*Y(K)
            ELSE
               TEMP = Y(K-1)
               Y(K-1) = Y(K)
               Y(K) = TEMP - C(K-1)*Y(K)
            END IF
         END DO
      END IF
!
      RETURN
      END SUBROUTINE DLAGTS
!
      FUNCTION DLAISNAN (DIN1, DIN2) RESULT (LAISNAN)
!
!=======================================================================
!                                                                      !
!  This routine is not for general use.  It exists solely to avoid     !
!  over-optimization in DISNAN.                                        !
!                                                                      !
!  DLAISNAN checks for NaNs by comparing its two arguments for         !
!  inequality.  NaN is the only floating-point value where NaN != NaN  !
!  returns .TRUE.  To check for NaNs, pass the same variable as both   !
!  arguments.                                                          !
!                                                                      !
!  A compiler must assume that the two arguments are                   !
!  not the same variable, and the test will not be optimized away.     !
!  Interprocedural or whole-program optimization may delete this       !
!  test.  The ISNAN functions will be replaced by the correct          !
!  Fortran 03 intrinsic once the intrinsic is widely available.        !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  DIN1     [in]      REAL                                             !
!                                                                      !
!  DIN2     [in]      REAL                                             !
!           Two numbers to compare for inequality.                     !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(IN) :: DIN1, DIN2
!
!  Local variable declarations.
!
      logical :: LAISNAN
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      LAISNAN = (DIN1.NE.DIN2)
!
      RETURN
      END FUNCTION DLAISNAN
!
      SUBROUTINE DLAMC1 (BETA, T, RND, IEEE1)
!
!=======================================================================
!                                                                      !
!  DLAMC1 determines the machine parameters given by BETA, T, RND, and !
!  IEEE1.                                                              !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  BETA    (output) INTEGER                                            !
!          The base of the machine.                                    !
!                                                                      !
!  T       (output) INTEGER                                            !
!          The number of (BETA) digits in the mantissa.                !
!                                                                      !
!  RND     (output) LOGICAL                                            !
!          Specifies whether proper rounding  (RND = .TRUE.)  or       !
!          chopping  (RND = .FALSE.)  occurs in addition. This may     !
!          not be a reliable guide to the way in which the machine     !
!          performs its arithmetic.                                    !
!                                                                      !
!  IEEE1   (output) LOGICAL                                            !
!          Specifies whether rounding appears to be done in the IEEE   !
!          'round to nearest' style.                                   !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  The routine is based on the routine  ENVRON  by Malcolm and         !
!  incorporates suggestions by Gentleman and Marovich. See             !
!                                                                      !
!     Malcolm M. A. (1972) Algorithms to reveal properties of          !
!        floating-point arithmetic. Comms. of the ACM, 15, 949-951.    !
!                                                                      !
!     Gentleman W. M. and Marovich S. B. (1974) More on algorithms     !
!        that reveal properties of floating point arithmetic units.    !
!        Comms. of the ACM, 17, 276-277.                               !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(out) :: IEEE1, RND
!
      integer, intent(out) :: BETA, T
!
!  Local variable declarations.
!
      logical, save :: FIRST = .TRUE.
      logical, save :: LIEEE1, LRND
!
      integer, save :: LBETA, LT
!
      real (dp) :: A, B, C, F, ONE, QTR, SAVEC, T1, T2
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      FIRST_PASS : IF (FIRST) THEN
        FIRST = .FALSE.
        ONE = 1
!
!  LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!  IEEE1, T and RND.
!
!  Throughout this routine  we use the function  DLAMC3  to ensure
!  that relevant values are  stored and not held in registers,  or
!  are not affected by optimizers.
!
!  Compute  a = 2.0**m  with the  smallest positive integer m such that
!
!    fl(a + 1.0) = a.
!
        A = 1
        C = 1
!
        DO WHILE (C .EQ. ONE)
          A = 2*A
          C = DLAMC3(A, ONE)
          C = DLAMC3(C, -A)
        END DO
!
!  Now compute  b = 2.0**m  with the smallest positive integer m
!  such that
!
!    fl(a + b) .gt. a.
!
        B = 1
        C = DLAMC3(A, B)
!
        DO WHILE (C .EQ. A)
          B = 2*B
          C = DLAMC3(A, B)
        END DO
!
!  Now compute the base.  a and c  are neighbouring floating point
!  numbers  in the  interval  (beta**t, beta**(t + 1))  and so
!  their difference is beta. Adding 0.25 to c is to ensure that it
!  is truncated to beta and not (beta - 1).
!
        QTR = ONE / 4
        SAVEC = C
        C = DLAMC3(C, -A)
        LBETA = C + QTR
!
!  Now determine whether rounding or chopping occurs,  by adding a
!  bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
        B = LBETA
        F = DLAMC3(B/2, -B/100)
        C = DLAMC3(F, A)
        IF (C .EQ. A) THEN
          LRND = .TRUE.
        ELSE
          LRND = .FALSE.
        END IF
        F = DLAMC3(B/2, B/100)
        C = DLAMC3(F, A)
        IF ((LRND) .AND. (C.EQ.A)) LRND = .FALSE.
!
!  Try and decide whether rounding is done in the  IEEE  round to
!  nearest style. B/2 is half a unit in the last place of the two
!  numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!  zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!  A, but adding B/2 to SAVEC should change SAVEC.
!
        T1 = DLAMC3(B/2, A)
        T2 = DLAMC3(B/2, SAVEC)
        LIEEE1 = (T1.EQ.A) .AND. (T2.GT.SAVEC) .AND. LRND
!
!  Now find  the  mantissa, t.  It should  be the  integer part of
!  log to the base beta of a,  however it is safer to determine  t
!  by powering.  So we find t as the smallest positive integer for
!  which
!
!    fl(beta**t + 1.0) = 1.0.
!
        LT = 0
        A = 1
        C = 1
!
        DO WHILE (C .EQ. ONE)
          LT = LT + 1
          A = A*LBETA
          C = DLAMC3(A, ONE)
          C = DLAMC3(C, -A)
        END DO
!
      END IF FIRST_PASS
!
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
!
      RETURN
      END SUBROUTINE DLAMC1
!
      SUBROUTINE DLAMC2 (BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX)
!
!=======================================================================
!                                                                      !
!  DLAMC2 determines the machine parameters specified in its argument  !
!  list.                                                               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  BETA    (output) INTEGER                                            !
!          The base of the machine.                                    !
!                                                                      !
!  T       (output) INTEGER                                            !
!          The number of (BETA) digits in the mantissa.                !
!                                                                      !
!  RND     (output) LOGICAL                                            !
!          Specifies whether proper rounding  (RND = .TRUE.)  or       !
!          chopping  (RND = .FALSE.)  occurs in addition. This may     !
!          not be a reliable guide to the way in which the machine     !
!          performs its arithmetic.                                    !
!                                                                      !
!  EPS     (output) DOUBLE PRECISION                                   !
!          The smallest positive number such that                      !
!                                                                      !
!             fl(1.0 - EPS) .LT. 1.0,                                  !
!                                                                      !
!          where fl denotes the computed value.                        !
!                                                                      !
!  EMIN    (output) INTEGER                                            !
!          The minimum exponent before (gradual) underflow occurs.     !
!                                                                      !
!  RMIN    (output) DOUBLE PRECISION                                   !
!          The smallest normalized number for the machine, given by    !
!          BASE**(EMIN - 1), where  BASE  is the floating point        !
!          value of BETA.                                              !
!                                                                      !
!  EMAX    (output) INTEGER                                            !
!          The maximum exponent before overflow occurs.                !
!                                                                      !
!  RMAX    (output) DOUBLE PRECISION                                   !
!          The largest positive number for the machine, given by       !
!          BASE**EMAX * (1 - EPS), where  BASE  is the floating        !
!          point value of BETA.                                        !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  The computation of  EPS  is based on a routine PARANOIA by          !
!  W. Kahan of the University of California at Berkeley.               !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical,   intent(out) :: RND
!
      integer,   intent(out) :: BETA, EMAX, EMIN, T
!
      real (dp), intent(out) :: EPS, RMAX, RMIN
!
!  Local variable declarations.
!
      logical, save :: FIRST = .TRUE.
      logical, save :: IWARN = .FALSE.
      logical       :: IEEE, LIEEE1, LRND
!
      integer, save :: LBETA, LEMAX, LEMIN, LT
      integer       :: GNMIN, GPMIN, I, NGNMIN, NGPMIN
!
      real (dp), save :: LEPS, LRMAX, LRMIN
      real (dp)       :: A, B, C, HALF, ONE, RBASE,                     &
     &                   SIXTH, SMALL, THIRD, TWO, ZERO
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      FIRST_PASS : IF (FIRST) THEN
        FIRST = .FALSE.
        ZERO = 0
        ONE = 1
        TWO = 2
!
!  LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!  BETA, T, RND, EPS, EMIN and RMIN.
!
!  Throughout this routine  we use the function  DLAMC3  to ensure
!  that relevant values are stored  and not held in registers,  or
!  are not affected by optimizers.
!
!  DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
        CALL DLAMC1 (LBETA, LT, LRND, LIEEE1)
!
!  Start to find EPS.
!
        B = LBETA
        A = B**(-LT)
        LEPS = A
!
!  Try some tricks to see whether or not this is the correct EPS.
!
        B = TWO / 3
        HALF = ONE / 2
        SIXTH = DLAMC3(B, -HALF)
        THIRD = DLAMC3(SIXTH, SIXTH)
        B = DLAMC3(THIRD, -HALF)
        B = DLAMC3(B, SIXTH)
        B = ABS(B)
        IF (B .LT. LEPS) B = LEPS
!
        LEPS = 1

        DO WHILE ((LEPS.GT.B) .AND. (B.GT.ZERO))
          LEPS = B
          C = DLAMC3(HALF*LEPS, (TWO**5)*(LEPS**2))
          C = DLAMC3(HALF, -C)
          B = DLAMC3(HALF, C)
          C = DLAMC3(HALF, -B)
          B = DLAMC3(HALF, C)
        END DO
!
        IF (A .LT. LEPS) LEPS = A
!
!  Computation of EPS complete.
!
!  Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!  Keep dividing  A by BETA until (gradual) underflow occurs. This
!  is detected when we cannot recover the previous A.
!
        RBASE = ONE / LBETA
        SMALL = ONE
        DO I = 1, 3
          SMALL = DLAMC3(SMALL*RBASE, ZERO)
        END DO
        A = DLAMC3(ONE, SMALL)
        CALL DLAMC4 (NGPMIN, ONE, LBETA)
        CALL DLAMC4 (NGNMIN, -ONE, LBETA)
        CALL DLAMC4 (GPMIN, A, LBETA)
        CALL DLAMC4 (GNMIN, -A, LBETA)
        IEEE = .FALSE.
!
        IF ((NGPMIN.EQ.NGNMIN) .AND. (GPMIN.EQ.GNMIN)) THEN
          IF (NGPMIN .EQ. GPMIN) THEN
            LEMIN = NGPMIN
!            (Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX)
          ELSE IF ((GPMIN-NGPMIN) .EQ. 3) THEN
            LEMIN = NGPMIN - 1 + LT
            IEEE = .TRUE.
!            (Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers)
          ELSE
            LEMIN = MIN(NGPMIN, GPMIN)
!            (A guess; no known machine)
            IWARN = .TRUE.
          END IF
!
        ELSE IF ((NGPMIN.EQ.GPMIN) .AND. (NGNMIN.EQ.GNMIN)) THEN
          IF (ABS(NGPMIN-NGNMIN) .EQ. 1) THEN
            LEMIN = MAX(NGPMIN, NGNMIN)
!            (Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205)
          ELSE
            LEMIN = MIN(NGPMIN, NGNMIN)
!            (A guess; no known machine)
            IWARN = .TRUE.
          END IF
!
        ELSE IF ((ABS(NGPMIN-NGNMIN).EQ.1) .AND. (GPMIN.EQ.GNMIN)) THEN
          IF ((GPMIN-MIN(NGPMIN, NGNMIN)) .EQ. 3) THEN
            LEMIN = MAX(NGPMIN, NGNMIN) - 1 + LT
!            (Twos-complement machines with gradual underflow;
!              no known machine)
          ELSE
            LEMIN = MIN(NGPMIN, NGNMIN)
!            (A guess; no known machine)
            IWARN = .TRUE.
          END IF
!
        ELSE
          LEMIN = MIN(NGPMIN, NGNMIN, GPMIN, GNMIN)
!         (A guess; no known machine)
          IWARN = .TRUE.
        END IF
!
!  Comment out this if block if EMIN is okay.
!
        IF (IWARN) THEN
          FIRST = .TRUE.
          PRINT 10, LEMIN
        END IF
!
!  Assume IEEE arithmetic if we found denormalised  numbers above,
!  or if arithmetic seems to round in the  IEEE style,  determined
!  in routine DLAMC1. A true IEEE machine should have both  things
!  true; however, faulty machines may have one or the other.
!
        IEEE = IEEE .OR. LIEEE1
!
!  Compute  RMIN by successive division by  BETA. We could compute
!  RMIN as BASE**(EMIN - 1),  but some machines underflow during
!  this computation.
!
        LRMIN = 1
        DO I = 1, 1 - LEMIN
          LRMIN = DLAMC3(LRMIN*RBASE, ZERO)
        END DO
!
!  Finally, call DLAMC5 to compute EMAX and RMAX.
!
        CALL DLAMC5 (LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX)
      END IF FIRST_PASS
!
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
!
 10   FORMAT (/ / ' WARNING. The value EMIN may be incorrect:-',       &
     &        '  EMIN = ', i0, /,                                       &
     &        ' If, after inspection, the value EMIN looks',            &
     &        ' acceptable please comment out ',                        &
     &        /, ' the IF block as marked within the code of routine',  &
     &        ' DLAMC2,', / ' otherwise supply EMIN explicitly.', /)
!
      RETURN
      END SUBROUTINE DLAMC2
!
      FUNCTION DLAMC3 (A, B) RESULT (AplusB)
!
!=======================================================================
!                                                                      !
!  DLAMC3  is intended to force A and B to be stored prior to doing    !
!  the addition of A and B , for use in situations where optimizers    !
!  might hold one of these in a register.                              !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  A, B    (input) DOUBLE PRECISION                                    !
!          The values A and B.                                         !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in) ::  A, B
!
!  Local variable declarations.
!
      real (dp) :: AplusB
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      AplusB = A + B
!
      RETURN
      END FUNCTION DLAMC3
!
      SUBROUTINE DLAMC4 (EMIN, START, BASE)
!
!=======================================================================
!                                                                      !
!  DLAMC4 is a service routine for DLAMC2.                             !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  EMIN    (output) EMIN                                               !
!          The minimum exponent before (gradual) underflow, computed   !
!          by setting A = START and dividing by BASE until the         !
!          previous A can not be recovered.                            !
!                                                                      !
!  START   (input) DOUBLE PRECISION                                    !
!          The starting point for determining EMIN.                    !
!                                                                      !
!  BASE    (input) INTEGER                                             !
!          The base of the machine.                                    !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in)   :: BASE
      integer, intent(out)  :: EMIN
!
      real (dp), intent(in) :: START
!
!  Local variable declarations.
!
      integer   :: I
!
      real (dp) :: A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3 (A*RBASE, ZERO)
      C1 = A
      C2 = A
      D1 = A
      D2 = A
!
      DO WHILE ((C1.EQ.A) .AND. (C2.EQ.A) .AND.                         &
     &          (D1.EQ.A) .AND. (D2.EQ.A))
        EMIN = EMIN - 1
        A = B1
        B1 = DLAMC3 (A / BASE, ZERO)
        C1 = DLAMC3 (B1*BASE, ZERO)
        D1 = ZERO
        DO I = 1, BASE
          D1 = D1 + B1
        END DO
        B2 = DLAMC3 (A*RBASE, ZERO)
        C2 = DLAMC3 (B2 / RBASE, ZERO)
        D2 = ZERO
        DO I = 1, BASE
          D2 = D2 + B2
        END DO
      END DO
!
      RETURN
      END SUBROUTINE DLAMC4
!
      SUBROUTINE DLAMC5 (BETA, P, EMIN, IEEE, EMAX, RMAX)
!
!=======================================================================
!                                                                      !
!  DLAMC5 attempts to compute RMAX, the largest machine floating-point !
!  number, without overflow.  It assumes that EMAX + abs(EMIN) sum     !
!  approximately to a power of 2.  It will fail on machines where this !
!  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,!
!  EMAX = 28718).  It will also fail if the value supplied for EMIN is !
!  too large (i.e. too close to zero), probably with overflow.         !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  BETA    (input) INTEGER                                             !
!          The base of floating-point arithmetic.                      !
!                                                                      !
!  P       (input) INTEGER                                             !
!          The number of base BETA digits in the mantissa of a         !
!          floating-point value.                                       !
!                                                                      !
!  EMIN    (input) INTEGER                                             !
!          The minimum exponent before (gradual) underflow.            !
!                                                                      !
!  IEEE    (input) LOGICAL                                             !
!          A logical flag specifying whether or not the arithmetic     !
!          system is thought to comply with the IEEE standard.         !
!                                                                      !
!  EMAX    (output) INTEGER                                            !
!          The largest exponent before overflow                        !
!                                                                      !
!  RMAX    (output) DOUBLE PRECISION                                   !
!          The largest machine floating-point number.                  !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical,   intent(in)  :: IEEE
!
      integer,   intent(in)  :: BETA, EMIN, P
      integer,   intent(out) :: EMAX
!
      real (dp), intent(out) :: RMAX
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      integer :: EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
!
      real (dp) :: OLDY, RECBAS, Y, Z
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  First compute LEXP and UEXP, two powers of 2 that bound
!  ABS(EMIN). We then assume that EMAX + ABS(EMIN) will sum
!  approximately to the bound that is closest to ABS(EMIN).
!  (EMAX is the exponent of the required number RMAX).
!
      LEXP = 1
      EXBITS = 1
#ifdef REMOVE_LAPACK_GOTOS
      TRY_LOOP : DO                              ! forever DO-loop
        TRY = LEXP*2
        IF (TRY .LE. (-EMIN)) THEN
          LEXP = TRY
          EXBITS = EXBITS + 1
          CYCLE TRY_LOOP                         ! iterate DO-loop again
        ELSE
          EXIT TRY_LOOP                          ! terminate DO-loop
        END IF
      END DO TRY_LOOP
#else
   10 CONTINUE
      TRY = LEXP*2
      IF (TRY .LE. (-EMIN)) THEN
        LEXP = TRY
        EXBITS = EXBITS + 1
        GOTO 10
      END IF
#endif
      IF (LEXP .EQ. -EMIN) THEN
        UEXP = LEXP
      ELSE
        UEXP = TRY
        EXBITS = EXBITS + 1
      END IF
!
!  Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!  than or equal to EMIN. EXBITS is the number of bits needed to
!  store the exponent.
!
      IF ((UEXP+EMIN) .GT. (-LEXP-EMIN)) THEN
        EXPSUM = 2*LEXP
      ELSE
        EXPSUM = 2*UEXP
      END IF
!
!  EXPSUM is the exponent range, approximately equal to
!  EMAX - EMIN + 1 .
!
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
!
!  NBITS is the total number of bits needed to store a
!  floating-point number.
!
      IF ((MOD(NBITS, 2).EQ.1) .AND. (BETA.EQ.2)) THEN
!
!  Either there are an odd number of bits used to store a
!  floating-point number, which is unlikely, or some bits are
!  not used in the representation of numbers, which is possible,
!  (e.g. Cray machines) or the mantissa has an implicit bit,
!  (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!  most likely. We have to assume the last alternative.
!  If this is true, then we need to reduce EMAX by one because
!  there must be some way of representing zero in an implicit-bit
!  system. On machines like Cray, we are reducing EMAX by one
!  unnecessarily.
!
        EMAX = EMAX - 1
      END IF
!
      IF (IEEE) THEN
!
!  Assume we are on an IEEE machine which reserves one exponent
!  for infinity and NaN.
!
        EMAX = EMAX - 1
      END IF
!
!  Now create RMAX, the largest machine number, which should
!  be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!  First compute 1.0 - BETA**(-P), being careful that the
!  result is less than 1.0 .
!
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO I = 1, P
        Z = Z*RECBAS
        IF (Y .LT. ONE) OLDY = Y
        Y = DLAMC3(Y, Z)
      END DO
      IF (Y .GE. ONE) Y = OLDY
!
!  Now multiply by BETA**EMAX to get RMAX.
!
      DO I = 1, EMAX
        Y = DLAMC3(Y*BETA, ZERO)
      END DO
!
      RMAX = Y
!
      RETURN
      END SUBROUTINE DLAMC5
!
      FUNCTION DLAMCH (CMACH) RESULT (RMACH)
!
!=======================================================================
!                                                                      !
!  DLAMCH determines double precision machine parameters.              !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  CMACH   (input) CHARACTER*1                                         !
!          Specifies the value to be returned by DLAMCH:               !
!          = 'E' or 'e',   DLAMCH := eps                               !
!          = 'S' or 's ,   DLAMCH := sfmin                             !
!          = 'B' or 'b',   DLAMCH := base                              !
!          = 'P' or 'p',   DLAMCH := eps*base                          !
!          = 'N' or 'n',   DLAMCH := t                                 !
!          = 'R' or 'r',   DLAMCH := rnd                               !
!          = 'M' or 'm',   DLAMCH := emin                              !
!          = 'U' or 'u',   DLAMCH := rmin                              !
!          = 'L' or 'l',   DLAMCH := emax                              !
!          = 'O' or 'o',   DLAMCH := rmax                              !
!                                                                      !
!          where                                                       !
!                                                                      !
!          eps   = relative machine precision                          !
!          sfmin = safe minimum, such that 1/sfmin does not overflow   !
!          base  = base of the machine                                 !
!          prec  = eps*base                                            !
!          t     = number of (base) digits in the mantissa             !
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise !
!          emin  = minimum exponent before (gradual) underflow         !
!          rmin  = underflow threshold - base**(emin-1)                !
!          emax  = largest exponent before overflow                    !
!          rmax  = overflow threshold  - (base**emax)*(1-eps)          !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), intent(in) :: CMACH
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      logical, save   :: FIRST = .TRUE.
      logical         :: LRND
!
      integer         :: BETA, IMAX, IMIN, IT
!
      real (dp), save :: EPS, SFMIN, BASE, T, RND, EMIN, RMIN,           &
     &                   EMAX, RMAX, PREC
      real (dp)       :: RMACH, SMALL
!
!-----------------------------------------------------------------------
!  Executable Statements ..
!-----------------------------------------------------------------------
!
      IF (FIRST) THEN
        FIRST = .FALSE.
        CALL DLAMC2 (BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX)
        BASE = BETA
        T = IT
        IF (LRND) THEN
          RND = ONE
          EPS = (BASE**(1-IT)) / 2
        ELSE
          RND = ZERO
          EPS = BASE**(1-IT)
        END IF
        PREC = EPS*BASE
        EMIN = IMIN
        EMAX = IMAX
        SFMIN = RMIN
        SMALL = ONE / RMAX
        IF (SMALL.GE.SFMIN) THEN
!
!  Use SMALL plus a bit, to avoid the possibility of rounding
!  causing overflow when computing  1/sfmin.
!
          SFMIN = SMALL*(ONE+EPS)
        END IF
      END IF
!
      IF (LSAME(CMACH, 'E')) THEN
         RMACH = EPS
      ELSE IF (LSAME(CMACH, 'S')) THEN
         RMACH = SFMIN
      ELSE IF (LSAME(CMACH, 'B')) THEN
         RMACH = BASE
      ELSE IF (LSAME(CMACH, 'P')) THEN
         RMACH = PREC
      ELSE IF (LSAME(CMACH, 'N')) THEN
         RMACH = T
      ELSE IF (LSAME(CMACH, 'R')) THEN
         RMACH = RND
      ELSE IF (LSAME(CMACH, 'M')) THEN
         RMACH = EMIN
      ELSE IF (LSAME(CMACH, 'U')) THEN
         RMACH = RMIN
      ELSE IF (LSAME(CMACH, 'L')) THEN
         RMACH = EMAX
      ELSE IF (LSAME(CMACH, 'O')) THEN
         RMACH = RMAX
      END IF
!
      RETURN
      END FUNCTION DLAMCH
!
      FUNCTION DLANEG (N, D, LLD, SIGMA, PIVMIN, R) RESULT (LANEG)
!
!=======================================================================
!                                                                      !
!  DLANEG computes the Sturm count, the number of negative pivots      !
!  encountered while factoring tridiagonal T - sigma I = L D L^T.      !
!  This implementation works directly on the factors without forming   !
!  the tridiagonal matrix T.  The Sturm count is also the number of    !
!  eigenvalues of T less than sigma.                                   !
!                                                                      !
!  This routine is called from DLARRB.                                 !
!                                                                      !
!  The current routine does not use the PIVMIN parameter but rather    !
!  requires IEEE-754 propagation of Infinities and NaNs.  This         !
!  routine also has no input range restrictions but does require       !
!  default exception handling such that x/0 produces Inf when x is     !
!  non-zero, and Inf/Inf produces NaN.  For more information, see:     !
!                                                                      !
!    Marques, Riedy, and Voemel, "Benefits of IEEE-754 Features in     !
!    Modern Symmetric Tridiagonal Eigensolvers," SIAM Journal on       !
!    Scientific Computing, v28, n5, 2006.  DOI 10.1137/050641624       !
!    (Tech report version in LAWN 172 with the same title.)            !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix.                                   !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           The N diagonal elements of the diagonal matrix D.          !
!                                                                      !
!  LLD      [in]      REAL array, dimension (N-1)                      !
!           The (N-1) elements L(i)*L(i)*D(i).                         !
!                                                                      !
!  SIGMA    [in]      REAL                                             !
!           Shift amount in T - sigma I = L D L^T.                     !
!                                                                      !
!  PIVMIN   [in]      REAL                                             !
!           The minimum pivot in the Sturm sequence.  May be used      !
!           when zero pivots are encountered on non-IEEE-754           !
!           architectures.                                             !
!                                                                      !
!  R        [in]      INTEGER                                          !
!           The twist index for the twisted factorization that is used !
!           for the negcount.                                          !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!  Jason Riedy, University of California, Berkeley, USA                !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in) :: N, R
!
      real (dp), intent(in) :: PIVMIN, SIGMA
      real (dp), intent(in) :: D(*), LLD(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
!     Some architectures propagate Infinities and NaNs very slowly, so
!     the code computes counts in BLKLEN chunks.  Then a NaN can
!     propagate at most BLKLEN columns before being detected.  This is
!     not a general tuning parameter; it needs only to be just large
!     enough that the overhead is tiny in common cases.
!
      integer, parameter :: BLKLEN = 128
!
      logical   :: SAWNAN
!
      integer   :: BJ, J, NEG1, NEG2, NEGCNT
      integer   :: LANEG
!
      real (dp) :: BSAV, DMINUS, DPLUS, GAMMA, P, T, TMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      NEGCNT = 0
!
!  I) upper part: L D L^T - SIGMA I = L+ D+ L+^T
!
      T = -SIGMA
      DO BJ = 1, R-1, BLKLEN
         NEG1 = 0
         BSAV = T
         DO J = BJ, MIN(BJ+BLKLEN-1, R-1)
            DPLUS = D(J) + T
            IF (DPLUS.LT.ZERO) NEG1 = NEG1 + 1
            TMP = T / DPLUS
            T = TMP * LLD(J) - SIGMA
         END DO
         SAWNAN = DISNAN(T)
!
!  Run a slower version of the above loop if a NaN is detected.
!  A NaN should occur only with a zero pivot after an infinite
!  pivot.  In that case, substituting 1 for T/DPLUS is the
!  correct limit.
!
         IF (SAWNAN) THEN
            NEG1 = 0
            T = BSAV
            DO J = BJ, MIN(BJ+BLKLEN-1, R-1)
               DPLUS = D(J) + T
               IF (DPLUS.LT.ZERO) NEG1 = NEG1 + 1
               TMP = T / DPLUS
               IF (DISNAN(TMP)) TMP = ONE
               T = TMP * LLD(J) - SIGMA
            END DO
         END IF
         NEGCNT = NEGCNT + NEG1
      END DO
!
!  II) lower part: L D L^T - SIGMA I = U- D- U-^T
!
      P = D(N) - SIGMA
      DO BJ = N-1, R, -BLKLEN
         NEG2 = 0
         BSAV = P
         DO J = BJ, MAX(BJ-BLKLEN+1, R), -1
            DMINUS = LLD(J) + P
            IF (DMINUS.LT.ZERO) NEG2 = NEG2 + 1
            TMP = P / DMINUS
            P = TMP * D(J) - SIGMA
         END DO
         SAWNAN = DISNAN(P)
!
!  As above, run a slower version that substitutes 1 for Inf/Inf.
!
         IF (SAWNAN) THEN
            NEG2 = 0
            P = BSAV
            DO J = BJ, MAX(BJ-BLKLEN+1, R), -1
               DMINUS = LLD(J) + P
               IF (DMINUS.LT.ZERO) NEG2 = NEG2 + 1
               TMP = P / DMINUS
               IF (DISNAN(TMP)) TMP = ONE
               P = TMP * D(J) - SIGMA
            END DO
         END IF
         NEGCNT = NEGCNT + NEG2
      END DO
!
!  III) Twist index.  T was shifted by SIGMA initially.
!
      GAMMA = (T + SIGMA) + P
      IF (GAMMA.LT.ZERO) NEGCNT = NEGCNT+1
      LANEG = NEGCNT

      RETURN
      END FUNCTION DLANEG
!
      FUNCTION DLANST (NORM, N, D, E) RESULT (LANST)
!
!=======================================================================
!                                                                      !
!  DLANST returns the value of the one norm, or the Frobenius norm,    !
!  or the infinity norm, or the element of largest absolute value of   !
!  a real symmetric tridiagonal matrix A.                              !
!                                                                      !
!  Description:                                                        !
!                                                                      !
!  DLANST returns the value                                            !
!                                                                      !
!      ANORM = (MAX(ABS(A(i,j))), NORM = 'M' or 'm'                    !
!              (                                                       !
!              (norm1(A),         NORM = '1', 'O' or 'o'               !
!              (                                                       !
!              (normI(A),         NORM = 'I' or 'i'                    !
!              (                                                       !
!              (normF(A),         NORM = 'F', 'f', 'E' or 'e'          !
!  where                                                               !
!                                                                      !
!  norm1 denotes the one norm of a matrix (maximum column sum),        !
!  normI denotes the infinity norm of a matrix (maximum row sum), and  !
!  normF denotes the Frobenius norm of a matrix (square root of sum of !
!        squares).                                                     !
!  Note that MAX(ABS(A(i,j))) is not a matrix norm.                    !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  NORM    (input) CHARACTER*1                                         !
!          Specifies the value to be returned in DLANST as described   !
!          above.                                                      !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The order of the matrix A.  N >= 0.  When N = 0, DLANST is  !
!          set to zero.                                                !
!                                                                      !
!  D       (input) DOUBLE PRECISION array, dimension (N)               !
!          The diagonal elements of A.                                 !
!                                                                      !
!  E       (input) DOUBLE PRECISION array, dimension (N-1)             !
!          The (n-1) sub-diagonal or super-diagonal elements of A.     !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     February 29, 1992                                                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in) :: N
!
      real (dp), intent(in) :: D(*), E(*)
!
      character (len=1), intent(in) :: NORM
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      integer   :: I
      real (dp) :: ANORM, SCALE, SUM
      real (dp) :: LANST
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (N.LE.0) THEN
        ANORM = ZERO
      ELSE IF (LSAME(NORM, 'M')) THEN
!
!  Find MAX(ABS(A(i,j))).
!
        ANORM = ABS(D(N))
        DO I = 1, N - 1
          ANORM = MAX(ANORM, ABS(D(I)))
          ANORM = MAX(ANORM, ABS(E(I)))
        END DO
      ELSE IF (LSAME(NORM, 'O') .OR. NORM.EQ.'1' .OR.                   &
     &         LSAME(NORM, 'I')) THEN
!
!  Find norm1(A).
!
        IF (N.EQ.1) THEN
          ANORM = ABS(D(1))
        ELSE
          ANORM = MAX(ABS(D(1)) + ABS(E(1)),                            &
     &                ABS(E(N-1)) + ABS(D(N)))
          DO I = 2, N - 1
            ANORM = MAX(ANORM, ABS(D(I)) + ABS(E(I)) + ABS(E(I-1)))
          END DO
        END IF
      ELSE IF ((LSAME(NORM, 'F')) .OR. (LSAME(NORM, 'E'))) THEN
!
!  Find normF(A).
!
        SCALE = ZERO
        SUM = ONE
        IF (N .GT. 1) THEN
          CALL DLASSQ (N-1, E, 1, SCALE, SUM)
          SUM = 2*SUM
        END IF
        CALL DLASSQ (N, D, 1, SCALE, SUM)
        ANORM = SCALE*SQRT(SUM)
      END IF
!
      LANST = ANORM
!
      RETURN
      END FUNCTION DLANST
!
      FUNCTION DLAPY2 (X, Y) RESULT (LAPY2)
!
!=======================================================================
!                                                                      !
!  DLAPY2 returns SQRT(x**2+y**2), taking care not to cause            !
!  unnecessary overflow.                                               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  X       (input) DOUBLE PRECISION                                    !
!  Y       (input) DOUBLE PRECISION                                    !
!          X and Y specify the values x and y.                         !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in) :: X, Y
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      real (dp) :: W, XABS, YABS, Z
      real (dp) :: LAPY2
!
!-----------------------------------------------------------------------
!  Executable Statements ..
!-----------------------------------------------------------------------
!
      XABS = ABS(X)
      YABS = ABS(Y)
      W = MAX(XABS, YABS)
      Z = MIN(XABS, YABS)
      IF (Z .EQ. ZERO) THEN
        LAPY2 = W
      ELSE
        LAPY2 = W*SQRT(ONE+(Z / W)**2)
      END IF
!
      RETURN
      END FUNCTION DLAPY2
!
      SUBROUTINE DLAR1V (N, B1, BN, LAMBDA, D, L, LD, LLD,              &
     &                   PIVMIN, GAPTOL, Z, WANTNC, NEGCNT, ZTZ,        &
     &                   MINGMA, R, ISUPPZ, NRMINV, RESID, RQCORR, WORK)
!
!=======================================================================
!                                                                      !
!  DLAR1V computes the (scaled) r-th column of the inverse of          !
!  the sumbmatrix in rows B1 through BN of the tridiagonal matrix      !
!  L D L**T - sigma I. When sigma is close to an eigenvalue, the       !
!  computed vector is an accurate eigenvector. Usually, r corresponds  !
!  to the index where the eigenvector is largest in magnitude.         !
!                                                                      !
!  The following steps accomplish this computation :                   !
!                                                                      !
!  (a) Stationary qd transform, L D L^T - sigma I = L(+) D(+) L(+)^T,  !
!  (b) Progressive qd transform, L D L^T - sigma I = U(-) D(-) U(-)^T, !
!  (c) Computation of the diagonal elements of the inverse of          !
!      L D L^T - sigma I by combining the above transforms, and        !
!      choosing r as the index where the diagonal of the inverse is    !
!      (one of the) largest in magnitude.                              !
!  (d) Computation of the (scaled) r-th column of the inverse using    !
!      the twisted factorization obtained by combining the top part of !
!      the stationary and the bottom part of the progressive transform.!
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix L D L^T.                           !
!                                                                      !
!  B1       [in]      INTEGER                                          !
!           First index of the submatrix of L D L^T.                   !
!                                                                      !
!  BN       [in]      INTEGER                                          !
!           Last index of the submatrix of L D L^T.                    !
!                                                                      !
!  LAMBDA   [in]      REAL                                             !
!            The shift. In order to compute an accurate eigenvector,   !
!            LAMBDA should be a good approximation to an eigenvalue    !
!            of L D L^T.                                               !
!                                                                      !
!  L        [in]      REAL  array, dimension (N-1)                     !
!           The (N-1) subdiagonal elements of the unit bidiagonal      !
!           matrix L, in elements 1 to N-1.                            !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!            The N diagonal elements of the diagonal matrix D.         !
!                                                                      !
!  LD       [in]      REAL array, dimension (N-1)                      !
!            The N-1 elements L(i)*D(i).                               !
!                                                                      !
!  LLD      [in]      REAL array, dimension (N-1)                      !
!            The N-1 elements L(i)*L(i)*D(i).                          !
!                                                                      !
!  PIVMIN   [in]      REAL                                             !
!           The minimum pivot in the Sturm sequence.                   !
!                                                                      !
!  GAPTOL   [in]      REAL                                             !
!           Tolerance that indicates when eigenvector entries are      !
!           negligible w.r.t. their contribution to the residual.      !
!                                                                      !
!  Z        [in,out]  REAL array, dimension (N)                        !
!           On input, all entries of Z must be set to 0.               !
!           On output, Z contains the (scaled) r-th column of the      !
!           inverse. The scaling is such that Z(R) equals 1.           !
!                                                                      !
!  WANTNC   [in]      LOGICAL                                          !
!           Specifies whether NEGCNT has to be computed.               !
!                                                                      !
!  NEGCNT   [out]     INTEGER                                          !
!           If WANTNC is .TRUE. then NEGCNT = the number of            !
!           pivots < pivmin in the matrix factorization L D L^T,       !
!           and NEGCNT = -1 otherwise.                                 !
!                                                                      !
!  ZTZ      [out]     REAL                                             !
!           The square of the 2-norm of Z.                             !
!                                                                      !
!  MINGMA   [out]     REAL                                             !
!           The reciprocal of the largest (in magnitude) diagonal      !
!           element of the inverse of L D L^T - sigma I.               !
!                                                                      !
!  R        [in,out]  INTEGER                                          !
!           The twist index for the twisted factorization used to      !
!           compute Z.                                                 !
!           On input, 0 <= R <= N. If R is input as 0, R is set to     !
!           the index where (L D L^T - sigma I)^{-1} is largest        !
!           in magnitude. If 1 <= R <= N, R is unchanged.              !
!           On output, R contains the twist index used to compute Z.   !
!           Ideally, R designates the position of the maximum entry    !
!           in the eigenvector.                                        !
!                                                                      !
!  ISUPPZ   [out]     INTEGER array, dimension (2)                     !
!           The support of the vector in Z, i.e., the vector Z is      !
!           nonzero only in elements ISUPPZ(1) through ISUPPZ(2).      !
!                                                                      !
!  NRMINV   [out]     REAL                                             !
!           NRMINV = 1/SQRT(ZTZ)                                       !
!                                                                      !
!  RESID    [out]     REAL                                             !
!           The residual of the FP vector.                             !
!           RESID = ABS(MINGMA)/SQRT(ZTZ)                              !
!                                                                      !
!  RQCORR   [out]     REAL                                             !
!           The Rayleigh Quotient correction to LAMBDA.                !
!           RQCORR = MINGMA*TMP                                        !
!                                                                      !
!  WORK     [out]     REAL array, dimension (4*N)                      !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical,   intent(in)    :: WANTNC
!
      integer,   intent(in)    :: B1, BN, N
      integer,   intent(inout) :: R
      integer,   intent(out)   :: NEGCNT
      integer,   intent(out)   :: ISUPPZ(*)
!
      real (dp), intent(in)    :: GAPTOL, LAMBDA, PIVMIN
      real (dp), intent(in)    :: D(*), L(*), LD(*), LLD(*)
      real (dp), intent(inout) :: Z(*)
      real (dp), intent(out)   :: MINGMA, NRMINV, RESID, RQCORR, ZTZ
      real (dp), intent(out)   :: WORK(*)
!
!  Imported variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      logical   :: SAWNAN1, SAWNAN2
      integer   :: I, INDLPL, INDP, INDS, INDUMN, NEG1, NEG2, R1, R2
      real (dp) :: DMINUS, DPLUS, EPS, S, TMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      EPS = DLAMCH('Precision')
!
      IF (R.EQ.0) THEN
         R1 = B1
         R2 = BN
      ELSE
         R1 = R
         R2 = R
      END IF
!
      INDLPL = 0          ! storage for LPLUS
      INDUMN = N          ! Storage for UMINUS
      INDS = 2*N + 1
      INDP = 3*N + 1
!
      IF (B1.EQ.1) THEN
         WORK(INDS) = ZERO
      ELSE
         WORK(INDS+B1-1) = LLD(B1-1)
      END IF
!
!  Compute the stationary transform (using the differential form)
!  until the index R2.
!
      SAWNAN1 = .FALSE.
      NEG1 = 0
      S = WORK(INDS+B1-1) - LAMBDA
      DO I = B1, R1 - 1
         DPLUS = D(I) + S
         WORK(INDLPL+I) = LD(I) / DPLUS
         IF (DPLUS.LT.ZERO) NEG1 = NEG1 + 1
         WORK(INDS+I) = S*WORK(INDLPL+I)*L(I)
         S = WORK(INDS+I) - LAMBDA
      END DO
      SAWNAN1 = DISNAN(S)
      IF (SAWNAN1) GOTO 60
      DO I = R1, R2 - 1
         DPLUS = D(I) + S
         WORK(INDLPL+I) = LD(I) / DPLUS
         WORK(INDS+I) = S*WORK(INDLPL+I)*L(I)
         S = WORK(INDS+I) - LAMBDA
      END DO
      SAWNAN1 = DISNAN(S)
!
 60   CONTINUE
      IF (SAWNAN1) THEN
!
!  Runs a slower version of the above loop if a NaN is detected.
!
         NEG1 = 0
         S = WORK(INDS+B1-1) - LAMBDA
         DO I = B1, R1 - 1
            DPLUS = D(I) + S
            IF (ABS(DPLUS).LT.PIVMIN) DPLUS = -PIVMIN
            WORK(INDLPL+I) = LD(I) / DPLUS
            IF (DPLUS.LT.ZERO) NEG1 = NEG1 + 1
            WORK(INDS+I) = S*WORK(INDLPL+I)*L(I)
            IF (WORK(INDLPL+I).EQ.ZERO) WORK(INDS+I) = LLD(I)
            S = WORK(INDS+I) - LAMBDA
         END DO
         DO I = R1, R2 - 1
            DPLUS = D(I) + S
            IF (ABS(DPLUS).LT.PIVMIN) DPLUS = -PIVMIN
            WORK(INDLPL+I) = LD(I) / DPLUS
            WORK(INDS+I) = S*WORK(INDLPL+I)*L(I)
            IF (WORK(INDLPL+I).EQ.ZERO) WORK(INDS+I) = LLD(I)
            S = WORK(INDS+I) - LAMBDA
         END DO
      END IF
!
!  Compute the progressive transform (using the differential form)
!  until the index R1.
!
      SAWNAN2 = .FALSE.
      NEG2 = 0
      WORK(INDP+BN-1) = D(BN) - LAMBDA
      DO I = BN - 1, R1, -1
         DMINUS = LLD(I) + WORK(INDP+I)
         TMP = D(I) / DMINUS
         IF (DMINUS.LT.ZERO) NEG2 = NEG2 + 1
         WORK(INDUMN+I) = L(I)*TMP
         WORK(INDP+I-1) = WORK(INDP+I)*TMP - LAMBDA
      END DO
      TMP = WORK(INDP+R1-1)
      SAWNAN2 = DISNAN(TMP)

      IF (SAWNAN2) THEN
!
!  Runs a slower version of the above loop if a NaN is detected.
!
         NEG2 = 0
         DO I = BN-1, R1, -1
            DMINUS = LLD(I) + WORK(INDP+I)
            IF (ABS(DMINUS).LT.PIVMIN) DMINUS = -PIVMIN
            TMP = D(I) / DMINUS
            IF (DMINUS.LT.ZERO) NEG2 = NEG2 + 1
            WORK(INDUMN+I) = L(I)*TMP
            WORK(INDP+I-1) = WORK(INDP+I)*TMP - LAMBDA
            IF (TMP.EQ.ZERO) WORK(INDP+I-1) = D(I) - LAMBDA
         END DO
      END IF
!
!  Find the index (from R1 to R2) of the largest (in magnitude)
!  diagonal element of the inverse.
!
      MINGMA = WORK(INDS+R1-1) + WORK(INDP+R1-1)
      IF (MINGMA.LT.ZERO) NEG1 = NEG1 + 1
      IF (WANTNC) THEN
         NEGCNT = NEG1 + NEG2
      ELSE
         NEGCNT = -1
      END IF
      IF (ABS(MINGMA).EQ.ZERO) MINGMA = EPS*WORK(INDS+R1-1)
      R = R1
      DO I = R1, R2 - 1
         TMP = WORK(INDS+I) + WORK(INDP+I)
         IF (TMP.EQ.ZERO) TMP = EPS*WORK(INDS+I)
         IF (ABS(TMP).LE.ABS(MINGMA)) THEN
            MINGMA = TMP
            R = I + 1
         END IF
      END DO
!
!  Compute the FP vector: solve N^T v = e_r.
!
      ISUPPZ(1) = B1
      ISUPPZ(2) = BN
      Z(R) = ONE
      ZTZ = ONE
!
!  Compute the FP vector upwards from R.
!
      IF (.NOT.SAWNAN1 .AND. .NOT.SAWNAN2) THEN
         DO I = R-1, B1, -1
            Z(I) = -(WORK(INDLPL+I)*Z(I+1))
            IF ((ABS(Z(I))+ABS(Z(I+1)))*ABS(LD(I)).LT.GAPTOL) THEN
               Z(I) = ZERO
               ISUPPZ(1) = I + 1
               GOTO 220
            END IF
            ZTZ = ZTZ + Z(I)*Z(I)
         END DO
 220     CONTINUE
      ELSE
!
!  Run slower loop if NaN occurred.
!
         DO I = R - 1, B1, -1
            IF (Z(I+1).EQ.ZERO) THEN
               Z(I) = -(LD(I+1) / LD(I))*Z(I+2)
            ELSE
               Z(I) = -(WORK(INDLPL+I)*Z(I+1))
            END IF
            IF ((ABS(Z(I))+ABS(Z(I+1)))*ABS(LD(I)).LT.GAPTOL) THEN
               Z(I) = ZERO
               ISUPPZ(1) = I + 1
               GOTO 240
            END IF
            ZTZ = ZTZ + Z(I)*Z(I)
         END DO
 240     CONTINUE
      END IF
!
!  Compute the FP vector downwards from R in blocks of size BLKSIZ.
!
      IF (.NOT.SAWNAN1 .AND. .NOT.SAWNAN2) THEN
         DO I = R, BN-1
            Z(I+1) = -(WORK(INDUMN+I)*Z(I))
            IF ((ABS(Z(I))+ABS(Z(I+1)))*ABS(LD(I)).LT.GAPTOL) THEN
               Z(I+1) = ZERO
               ISUPPZ(2) = I
               GOTO 260
            END IF
            ZTZ = ZTZ + Z(I+1)*Z(I+1)
         END DO
 260     CONTINUE
      ELSE
!
!  Run slower loop if NaN occurred.
!
         DO I = R, BN - 1
            IF (Z(I).EQ.ZERO) THEN
               Z(I+1) = -(LD(I-1) / LD(I))*Z(I-1)
            ELSE
               Z(I+1) = -(WORK(INDUMN+I)*Z(I))
            END IF
            IF ((ABS(Z(I))+ABS(Z(I+1)))*ABS(LD(I)).LT.GAPTOL) THEN
               Z(I+1) = ZERO
               ISUPPZ(2) = I
               GOTO 280
            END IF
            ZTZ = ZTZ + Z(I+1)*Z(I+1)
         END DO
 280     CONTINUE
      END IF
!
!  Compute quantities for convergence test.
!
      TMP = ONE / ZTZ
      NRMINV = SQRT(TMP)
      RESID = ABS(MINGMA)*NRMINV
      RQCORR = MINGMA*TMP
!
      RETURN
      END SUBROUTINE DLAR1V
!
      SUBROUTINE DLARNV (IDIST, ISEED, N, X)
!
!=======================================================================
!                                                                      !
!  DLARNV returns a vector of n random real numbers from a uniform or  !
!  normal distribution.                                                !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  IDIST    [in]      INTEGER                                          !
!           Specifies the distribution of the random numbers:          !
!           = 1:  uniform (0,1)                                        !
!           = 2:  uniform (-1,1)                                       !
!           = 3:  normal (0,1)                                         !
!                                                                      !
!  ISEED    [in,out]  INTEGER array, dimension (4)                     !
!           On entry, the seed of the random number generator; the     !
!           array elements must be between 0 and 4095, and ISEED(4)    !
!           must be odd.                                               !
!           On exit, the seed is updated.                              !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The number of random numbers to be generated.              !
!                                                                      !
!  X        [out]     REAL array, dimension (N)                        !
!           The generated random numbers.                              !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  This routine calls the auxiliary routine DLARUV to generate random  !
!  real numbers from a uniform (0,1) distribution, in batches of up to !
!  128 using vectorisable code. The Box-Muller method is used to       !
!  transform numbers from a uniform to a normal distribution.          !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: IDIST, N
      integer,   intent(inout) :: ISEED(4)
!
      real (dp), intent(out)   :: X(*)
!
!  Local variable declarations.
!
      integer,   parameter :: LV = 128
!
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: TWOPI =                                  &
     &                     6.28318530717958647692528676655900576839_dp
!
      integer   :: I, IL, IL2, IV
!
      real (dp) :: U(LV)
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      DO IV = 1, N, LV / 2
         IL = MIN(LV / 2, N-IV+1)
         IF (IDIST.EQ.3) THEN
            IL2 = 2*IL
         ELSE
            IL2 = IL
         END IF
!
!  Use DLARUV to generate IL2 numbers from a uniform (0,1)
!  distribution (IL2 <= LV).
!
         CALL DLARUV (ISEED, IL2, U)
!
         IF (IDIST.EQ.1) THEN
!
!  Copy generated numbers.
!
            DO I = 1, IL
               X(IV+I-1) = U(I)
            END DO
         ELSE IF (IDIST.EQ.2) THEN
!
!  Convert generated numbers to uniform (-1,1) distribution.
!
            DO I = 1, IL
               X(IV+I-1) = TWO*U(I) - ONE
            END DO
         ELSE IF (IDIST.EQ.3) THEN
!
!  Convert generated numbers to normal (0,1) distribution.
!
            DO I = 1, IL
               X(IV+I-1) = SQRT(-TWO*LOG(U(2*I-1)))*                    &
     &                     COS(TWOPI*U(2*I))
            END DO
         END IF
      END DO
!
      RETURN
      END SUBROUTINE DLARNV
!
      SUBROUTINE DLARRA (N, D, E, E2, SPLTOL, TNRM,                     &
     &                   NSPLIT, ISPLIT, INFO)
!
!=======================================================================
!                                                                      !
!  Compute the splitting points with threshold SPLTOL.                 !
!  DLARRA sets any "small" off-diagonal elements to zero.              !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix. N > 0.                            !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           On entry, the N diagonal elements of the tridiagonal       !
!           matrix T.                                                  !
!                                                                      !
!  E        [in,out]  REAL array, dimension (N)                        !
!           On entry, the first (N-1) entries contain the subdiagonal  !
!           elements of the tridiagonal matrix T; E(N) need not be set.!
!           On exit, the entries E(ISPLIT(I)), 1 <= I <= NSPLIT,       !
!           are set to zero, the other entries of E are untouched.     !
!                                                                      !
!  E2       [in,out]  REAL array, dimension (N)                        !
!           On entry, the first (N-1) entries contain the SQUARES of   !
!           the subdiagonal elements of the tridiagonal matrix T;      !
!           E2(N) need not be set.                                     !
!           On exit, the entries E2(ISPLIT(I)), 1 <= I <= NSPLIT,      !
!           have been set to zero.                                     !
!                                                                      !
!  SPLTOL   [in]      REAL                                             !
!           The threshold for splitting. Two criteria can be used:     !
!           SPLTOL<0 : criterion based on absolute off-diagonal value  !
!           SPLTOL>0 : criterion that preserves relative accuracy      !
!                                                                      !
!  TNRM     [in]      REAL                                             !
!           The norm of the matrix.                                    !
!                                                                      !
!  NSPLIT   [out]     INTEGER                                          !
!           The number of blocks T splits into, 1 <= NSPLIT <= N.      !
!                                                                      !
!  ISPLIT   [out]     INTEGER array, dimension (N)                     !
!           The splitting points, at which T breaks up into blocks.    !
!           The first block consists of rows/columns 1 to ISPLIT(1),   !
!           the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),  !
!           etc., and the NSPLIT-th consists of rows/columns           !
!           ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.               !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           = 0:  successful exit.                                     !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: N
      integer,   intent(out)   :: INFO, NSPLIT
      integer,   intent(out)   :: ISPLIT(*)
!
      real (dp), intent(in)    :: SPLTOL, TNRM
      real (dp), intent(in)    :: D(*)
      real (dp), intent(inout) :: E(*), E2(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
!
      integer   :: I
      real (dp) :: EABS, TMP1
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
      NSPLIT = 1
!
!  Quick return if possible.
!
      IF (N.LE.0) THEN
         RETURN
      END IF
!
!  Compute splitting points.
!
      IF (SPLTOL.LT.ZERO) THEN
!
!  Criterion based on absolute off-diagonal value.
!
         TMP1 = ABS(SPLTOL)* TNRM
         DO I = 1, N-1
            EABS = ABS(E(I))
            IF (EABS .LE. TMP1) THEN
               E(I) = ZERO
               E2(I) = ZERO
               ISPLIT(NSPLIT) = I
               NSPLIT = NSPLIT + 1
            END IF
         END DO
      ELSE
!
!  Criterion that guarantees relative accuracy.
!
         DO I = 1, N-1
            EABS = ABS(E(I))
            IF (EABS .LE.                                               &
     &          SPLTOL * SQRT(ABS(D(I)))*SQRT(ABS(D(I+1)))) THEN
               E(I) = ZERO
               E2(I) = ZERO
               ISPLIT(NSPLIT) = I
               NSPLIT = NSPLIT + 1
            END IF
         END DO
      END IF
      ISPLIT(NSPLIT) = N
!
      RETURN
      END SUBROUTINE DLARRA
!
      SUBROUTINE DLARRB (N, D, LLD, IFIRST, ILAST, RTOL1,               &
     &                   RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK,     &
     &                   PIVMIN, SPDIAM, TWIST, INFO)
!
!=======================================================================
!                                                                      !
!  Given the relatively robust representation(RRR) L D L^T, DLARRB     !
!  does "limited" bisection to refine the eigenvalues of L D L^T,      !
!  W(IFIRST-OFFSET) through W(ILAST-OFFSET), to more accuracy.         !
!  Initial guesses for these eigenvalues are input in W, the           !
!  corresponding estimate of the error in these guesses and their      !
!  gaps are input in WERR and WGAP, respectively. During bisection,    !
!  intervals [left, right] are maintained by storing their mid-points  !
!  and semi-widths in the arrays W and WERR respectively.              !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix.                                   !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           The N diagonal elements of the diagonal matrix D.          !
!                                                                      !
!  LLD      [in]      REAL array, dimension (N-1)                      !
!           The (N-1) elements L(i)*L(i)*D(i).                         !
!                                                                      !
!  IFIRST   [in]      INTEGER                                          !
!           The index of the first eigenvalue to be computed.          !
!                                                                      !
!  ILAST    [in]      INTEGER                                          !
!           The index of the last eigenvalue to be computed.           !
!                                                                      !
!  RTOL1    [in]      REAL                                             !
!                                                                      !
!  RTOL2    [in]      REAL                                             !
!           Tolerance for the convergence of the bisection intervals.  !
!           An interval [LEFT,RIGHT] has converged if                  !
!           RIGHT-LEFT < MAX(RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|))     !
!           where GAP is the (estimated) distance to the nearest       !
!           eigenvalue.                                                !
!                                                                      !
!  OFFSET   [in]      INTEGER                                          !
!           Offset for the arrays W, WGAP and WERR, i.e., the          !
!           IFIRST-OFFSET through ILAST-OFFSET elements of these       !
!           arrays are to be used.                                     !
!                                                                      !
!  W        [in,out]  REAL array, dimension (N)                        !
!           On input, W(IFIRST-OFFSET) through W(ILAST-OFFSET) are     !
!           estimates of the eigenvalues of L D L^T indexed IFIRST     !
!           through ILAST.                                             !
!           On output, these estimates are refined.                    !
!                                                                      !
!  WGAP     [in,out]  REAL array, dimension (N-1)                      !
!           On input, the (estimated) gaps between consecutive         !
!           eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap    !
!           between eigenvalues I and I+1. Note that if IFIRST = ILAST !
!           then WGAP(IFIRST-OFFSET) must be set to ZERO.              !
!           On output, these gaps are refined.                         !
!                                                                      !
!  WERR     [in,out]  REAL array, dimension (N)                        !
!           On input, WERR(IFIRST-OFFSET) through WERR(ILAST-OFFSET)   !
!           are the errors in the estimates of the corresponding       !
!           elements in W.                                             !
!           On output, these errors are refined.                       !
!                                                                      !
!  WORK     [out]     REAL array, dimension (2*N)                      !
!           Workspace.                                                 !
!                                                                      !
!  IWORK    [out]     INTEGER array, dimension (2*N)                   !
!           Workspace.                                                 !
!                                                                      !
!  PIVMIN   [in]      REAL                                             !
!           The minimum pivot in the Sturm sequence.                   !
!                                                                      !
!  SPDIAM   [in]      REAL                                             !
!           The spectral diameter of the matrix.                       !
!                                                                      !
!  TWIST    [in]      INTEGER                                          !
!           The twist index for the twisted factorization that is used !
!           for the negcount.                                          !
!           TWIST = N: Compute negcount from                           !
!                      L D L^T - LAMBDA I = L+ D+ L+^T                 !
!           TWIST = 1: Compute negcount from                           !
!                      L D L^T - LAMBDA I = U- D- U-^T                 !
!           TWIST = R: Compute negcount from                           !
!                      L D L^T - LAMBDA I = N(r) D(r) N(r)             !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           Error flag.                                                !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)  :: IFIRST, ILAST, N, OFFSET, TWIST
      integer,   intent(out) :: INFO
      integer,   intent(out) :: IWORK(*)
!
      real (dp), intent(in)    :: PIVMIN, RTOL1, RTOL2, SPDIAM
      real (dp), intent(in)    :: D(*), LLD(*)
      real (dp), intent(inout) :: W(*), WERR(*), WGAP(*)
      real (dp), intent(out)   :: WORK(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: HALF = 0.5_dp
!
      integer :: MAXITR
      integer :: I, I1, II, IP, ITER, K, NEGCNT, NEXT, NINT
      integer :: OLNINT, PREV, R
!
      real (dp) :: BACK, CVRGD, GAP, LEFT, LGAP, MID, MNWDTH
      real (dp) :: RGAP, RIGHT, TMP, WIDTH
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
!
!  Quick return if possible.
!
      IF (N.LE.0) THEN
         RETURN
      END IF
!
      MAXITR = INT((LOG(SPDIAM+PIVMIN)-LOG(PIVMIN)) /  LOG(TWO)) + 2
      MNWDTH = TWO * PIVMIN
!
      R = TWIST
      IF ((R.LT.1).OR.(R.GT.N)) R = N
!
!  Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
!  The Sturm Count, Count(WORK(2*I-1)) is arranged to be I-1, while
!  Count(WORK(2*I)) is stored in IWORK(2*I).
!  The integer IWORK(2*I-1)  for an unconverged interval is set to
!  the index of the next unconverged interval, and is -1 or 0 for a
!  converged interval. Thus a linked list of unconverged intervals is
!  set up.
!
      I1 = IFIRST
!
!  The number of unconverged intervals.
!
      NINT = 0
!
!  The last unconverged interval found
!
      PREV = 0
!
      RGAP = WGAP(I1-OFFSET)
      DO 75 I = I1, ILAST
         K = 2*I
         II = I - OFFSET
         LEFT = W(II) - WERR(II)
         RIGHT = W(II) + WERR(II)
         LGAP = RGAP
         RGAP = WGAP(II)
         GAP = MIN(LGAP, RGAP)
!
!  Make sure that [LEFT,RIGHT] contains the desired eigenvalue
!  Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT
!
!        DO WHILE (NEGCNT(LEFT).GT.I-1)
!
         BACK = WERR(II)
 20      CONTINUE
         NEGCNT = DLANEG(N, D, LLD, LEFT, PIVMIN, R)
         IF (NEGCNT.GT.I-1) THEN
            LEFT = LEFT - BACK
            BACK = TWO*BACK
            GOTO 20
         END IF
!
!  Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT
!
!        DO WHILE (NEGCNT(RIGHT).LT.I)
!
         BACK = WERR(II)
 50      CONTINUE

         NEGCNT = DLANEG(N, D, LLD, RIGHT, PIVMIN, R)
         IF (NEGCNT.LT.I) THEN
            RIGHT = RIGHT + BACK
            BACK = TWO*BACK
            GOTO 50
         END IF
         WIDTH = HALF*ABS(LEFT - RIGHT)
         TMP = MAX(ABS(LEFT), ABS(RIGHT))
         CVRGD = MAX(RTOL1*GAP,RTOL2*TMP)
         IF (WIDTH.LE.CVRGD .OR. WIDTH.LE.MNWDTH) THEN
!
!  This interval has already converged and does not need refinement.
!  (Note that the gaps might change through refining the
!  eigenvalues, however, they can only get bigger.)
!  Remove it from the list.
!
            IWORK(K-1) = -1
!
!  Make sure that I1 always points to the first unconverged interval.
!
            IF ((I.EQ.I1).AND.(I.LT.ILAST)) I1 = I + 1
            IF ((PREV.GE.I1).AND.(I.LE.ILAST)) IWORK(2*PREV-1) = I + 1
         ELSE
!
!  Unconverged interval found.
!
            PREV = I
            NINT = NINT + 1
            IWORK(K-1) = I + 1
            IWORK(K) = NEGCNT
         END IF
         WORK(K-1) = LEFT
         WORK(K) = RIGHT
 75   CONTINUE
!
!  There are still unconverged intervals
!
!     DO WHILE  (NINT.GT.0)
!     and WHILE (ITER.LT.MAXITR)
!
      ITER = 0
 80   CONTINUE
      PREV = I1 - 1
      I = I1
      OLNINT = NINT

      DO 100 IP = 1, OLNINT
         K = 2*I
         II = I - OFFSET
         RGAP = WGAP(II)
         LGAP = RGAP
         IF (II.GT.1) LGAP = WGAP(II-1)
         GAP = MIN(LGAP, RGAP)
         NEXT = IWORK(K-1)
         LEFT = WORK(K-1)
         RIGHT = WORK(K)
         MID = HALF*(LEFT + RIGHT)
!
!  Semiwidth of interval
!
         WIDTH = RIGHT - MID
         TMP = MAX(ABS(LEFT), ABS(RIGHT))
         CVRGD = MAX(RTOL1*GAP,RTOL2*TMP)
         IF ((WIDTH.LE.CVRGD) .OR. (WIDTH.LE.MNWDTH) .OR.               &
     &       (ITER.EQ.MAXITR)) THEN
!
!  Reduce number of unconverged intervals
!
            NINT = NINT - 1
!
!  Mark interval as converged.
!
            IWORK(K-1) = 0
            IF (I1.EQ.I) THEN
               I1 = NEXT
            ELSE
!
!  Prev holds the last unconverged interval previously examined.
!
               IF (PREV.GE.I1) IWORK(2*PREV-1) = NEXT
            END IF
            I = NEXT
            GOTO 100
         END IF
         PREV = I
!
!  Perform one bisection step.
!
         NEGCNT = DLANEG(N, D, LLD, MID, PIVMIN, R)
         IF (NEGCNT.LE.I-1) THEN
            WORK(K-1) = MID
         ELSE
            WORK(K) = MID
         END IF
         I = NEXT
 100  CONTINUE
      ITER = ITER + 1
!
!  Do another loop if there are still unconverged intervals
!  However, in the last iteration, all intervals are accepted
!  since this is the best we can do.
!
      IF ((NINT.GT.0).AND.(ITER.LE.MAXITR)) GOTO 80
!
!  At this point, all the intervals have converged.
!
      DO I = IFIRST, ILAST
         K = 2*I
         II = I - OFFSET
!
!  All intervals marked by '0' have been refined.
!
         IF (IWORK(K-1).EQ.0) THEN
            W(II) = HALF*(WORK(K-1)+WORK(K))
            WERR(II) = WORK(K) - W(II)
         END IF
      END DO
!
      DO I = IFIRST+1, ILAST
         K = 2*I
         II = I - OFFSET
         WGAP(II-1) = MAX(ZERO,                                         &
     &                    W(II)-WERR(II)-W(II-1)-WERR(II-1))
      END DO
!
      RETURN
      END SUBROUTINE DLARRB
!
      SUBROUTINE DLARRC (JOBT, N, VL, VU, D, E, PIVMIN,                 &
     &                   EIGCNT, LCNT, RCNT, INFO)
!
!=======================================================================
!                                                                      !
!  Find the number of eigenvalues of the symmetric tridiagonal matrix  !
!  T that are in the interval (VL,VU] if JOBT = 'T', and of L D L^T    !
!  if JOBT = 'L'.                                                      !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  JOBT     [in]      CHARACTER*1                                      !
!           = 'T':  Compute Sturm count for matrix T.                  !
!           = 'L':  Compute Sturm count for matrix L D L^T.            !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix. N > 0.                            !
!                                                                      !
!  VL       [in]      REAL                                             !
!           The lower bound for the eigenvalues.                       !
!                                                                      !
!  VU       [in]      REAL                                             !
!           The upper bound for the eigenvalues.                       !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           JOBT = 'T': The N diagonal elements of the tridiagonal     !
!                       matrix T.                                      !
!           JOBT = 'L': The N diagonal elements of the diagonal        !
!                       matrix D.                                      !
!                                                                      !
!  E        [in]      REAL array, dimension (N)                        !
!           JOBT = 'T': The N-1 offdiagonal elements of the matrix T.  !
!           JOBT = 'L': The N-1 offdiagonal elements of the matrix L.  !
!                                                                      !
!  PIVMIN   [in]      REAL                                             !
!           The minimum pivot in the Sturm sequence for T.             !
!                                                                      !
!  EIGCNT   [out]     INTEGER                                          !
!           The number of eigenvalues of the symmetric tridiagonal     !
!           matrix T that are in the interval (VL,VU].                 !
!                                                                      !
!  LCNT     [out]     INTEGER                                          !
!                                                                      !
!  RCNT     [out]     INTEGER                                          !
!           The left and right negcounts of the interval.              !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), intent(in) :: JOBT
!
      integer,   intent(in)  :: N
      integer,   intent(out) :: EIGCNT, INFO, LCNT, RCNT
!
      real (dp), intent(in)  :: PIVMIN, VL, VU
      real (dp), intent(in)  :: D(*), E(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
!
      logical   :: MATT
!
      integer   :: I
!
      real (dp) :: LPIVOT, RPIVOT, SL, SU, TMP, TMP2
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
      LCNT = 0
      RCNT = 0
      EIGCNT = 0
!
!  Quick return if possible.
!
      IF (N.LE.0) THEN
         RETURN
      END IF
!
      MATT = LSAME(JOBT, 'T')
!
      IF (MATT) THEN
!
!  Sturm sequence count on T
!
         LPIVOT = D(1) - VL
         RPIVOT = D(1) - VU
         IF (LPIVOT.LE.ZERO) THEN
            LCNT = LCNT + 1
         END IF
         IF (RPIVOT.LE.ZERO) THEN
            RCNT = RCNT + 1
         END IF
         DO I = 1, N-1
            TMP = E(I)**2
            LPIVOT = (D(I+1)-VL) - TMP/LPIVOT
            RPIVOT = (D(I+1)-VU) - TMP/RPIVOT
            IF (LPIVOT.LE.ZERO) THEN
               LCNT = LCNT + 1
            END IF
            IF (RPIVOT.LE.ZERO) THEN
               RCNT = RCNT + 1
            END IF
         END DO
      ELSE
!
!  Sturm sequence count on L D L^T.
!
         SL = -VL
         SU = -VU
         DO I = 1, N - 1
            LPIVOT = D(I) + SL
            RPIVOT = D(I) + SU
            IF (LPIVOT.LE.ZERO) THEN
               LCNT = LCNT + 1
            END IF
            IF (RPIVOT.LE.ZERO) THEN
               RCNT = RCNT + 1
            END IF
            TMP = E(I) * D(I) * E(I)
!
            TMP2 = TMP / LPIVOT
            IF (TMP2.EQ.ZERO) THEN
               SL =  TMP - VL
            ELSE
               SL = SL*TMP2 - VL
            END IF
!
            TMP2 = TMP / RPIVOT
            IF (TMP2.EQ.ZERO) THEN
               SU =  TMP - VU
            ELSE
               SU = SU*TMP2 - VU
            END IF
         END DO
         LPIVOT = D(N) + SL
         RPIVOT = D(N) + SU
         IF (LPIVOT.LE.ZERO) THEN
            LCNT = LCNT + 1
         END IF
         IF (RPIVOT.LE.ZERO) THEN
            RCNT = RCNT + 1
         END IF
      END IF
      EIGCNT = RCNT - LCNT
!
      RETURN
      END SUBROUTINE DLARRC
!
      SUBROUTINE DLARRD (RANGE, ORDER, N, VL, VU, IL, IU, GERS,         &
     &                   RELTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT,      &
     &                   M, W, WERR, WL, WU, IBLOCK, INDEXW,            &
     &                   WORK, IWORK, INFO)
!
!=======================================================================
!                                                                      !
!  DLARRD computes the eigenvalues of a symmetric tridiagonal matrix T !
!  to suitable accuracy. This is an auxiliary code to be called from   !
!  DSTEMR. The user may ask for all eigenvalues, all eigenvalues in    !
!  the half-open interval (VL, VU], or the IL-th through IU-th         !
!  eigenvalues.                                                        !
!                                                                      !
!  To avoid overflow, the matrix must be scaled so that its largest    !
!  element is no greater than overflow**(1/2) * underflow**(1/4) in    !
!  absolute value, and for greatest accuracy, it should not be much    !
!  smaller than that.                                                  !
!                                                                      !
!  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal       !
!  Matrix", Report CS41, Computer Science Dept., Stanford              !
!  University, July 21, 1966.                                          !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  RANGE    [in]      CHARACTER*1                                      !
!           = 'A': ("All")   all eigenvalues will be found.            !
!           = 'V': ("Value") all eigenvalues in the half-open interval !
!                            (VL, VU] will be found.                   !
!           = 'I': ("Index") the IL-th through IU-th eigenvalues (of   !
!                            the entire matrix) will be found.         !
!                                                                      !
!  ORDER    [in]      CHARACTER*1                                      !
!           = 'B': ("By Block") the eigenvalues will be grouped by     !
!                               split-off block (see IBLOCK, ISPLIT)   !
!                               and ordered from smallest to largest   !
!                               within the block.                      !
!           = 'E': ("Entire matrix")                                   !
!                               the eigenvalues for the entire matrix  !
!                               will be ordered from smallest to       !
!                               largest.                               !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the tridiagonal matrix T.  N >= 0.            !
!                                                                      !
!  VL       [in]      REAL                                             !
!           If RANGE='V', the lower bound of the interval to be        !
!           searched for eigenvalues.  Eigenvalues less than or equal  !
!           to VL, or greater than VU, will not be returned.  VL < VU. !
!           Not referenced if RANGE = 'A' or 'I'.                      !
!                                                                      !
!  VU       [in]      REAL                                             !
!           If RANGE='V', the upper bound of the interval to be        !
!           searched for eigenvalues.  Eigenvalues less than or equal  !
!           to VL, or greater than VU, will not be returned.  VL < VU. !
!           Not referenced if RANGE = 'A' or 'I'.                      !
!                                                                      !
!  IL       [in]      INTEGER                                          !
!           If RANGE='I', the index of the smallest eigenvalue to be   !
!           returned.                                                  !
!           1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.  !
!           Not referenced if RANGE = 'A' or 'V'.                      !
!                                                                      !
!  IU       [in]      INTEGER                                          !
!           If RANGE='I', the index of the largest eigenvalue to be    !
!           returned.                                                  !
!           1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.  !
!           Not referenced if RANGE = 'A' or 'V'.                      !
!                                                                      !
!  GERS     [in]      REAL array, dimension (2*N)                      !
!           The N Gerschgorin intervals (the i-th Gerschgorin interval !
!           is (GERS(2*i-1), GERS(2*i)).                               !
!                                                                      !
!  RELTOL   [in]      REAL                                             !
!           The minimum relative width of an interval.  When an        !
!           interval is narrower than RELTOL times the larger (in      !
!           magnitude) endpoint, then it is considered to be           !
!           sufficiently small, i.e., converged.  Note: this should    !
!           always be at least radix*machine epsilon.                  !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           The n diagonal elements of the tridiagonal matrix T.       !
!                                                                      !
!  E        [in]      REAL array, dimension (N-1)                      !
!           The (N-1) off-diagonal elements of the tridiagonal         !
!           matrix T.                                                  !
!                                                                      !
!  E2       [in]      REAL array, dimension (N-1)                      !
!           The (N-1) squared off-diagonal elements of the tridiagonal !
!           matrix T.                                                  !
!                                                                      !
!  PIVMIN   [in]      REAL                                             !
!           The minimum pivot allowed in the Sturm sequence for T.     !
!                                                                      !
!  NSPLIT   [in]      INTEGER                                          !
!           The number of diagonal blocks in the matrix T.             !
!           1 <= NSPLIT <= N.                                          !
!                                                                      !
!  ISPLIT   [in]                                                       !
!           ISPLIT is INTEGER array, dimension (N)                     !
!           The splitting points, at which T breaks up into            !
!           submatrices. The first submatrix consists of rows/columns  !
!           1 to ISPLIT(1), the second of rows/columns ISPLIT(1)+1     !
!           through ISPLIT(2), etc., and the NSPLIT-th consists of     !
!           rows/columns ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.  !
!           (Only the first NSPLIT elements will actually be used, but !
!           since the user cannot know a priori what value NSPLIT will !
!           have, N words must be reserved for ISPLIT.)                !
!                                                                      !
!  M        [out]     INTEGER                                          !
!           The actual number of eigenvalues found. 0 <= M <= N.       !
!           (See also the description of INFO=2,3.)                    !
!                                                                      !
!  W        [out]     REAL array, dimension (N)                        !
!           On exit, the first M elements of W will contain the        !
!           eigenvalue approximations. DLARRD computes an interval     !
!           I_j = (a_j, b_j] that includes eigenvalue j. The           !
!           eigenvalue approximation is given as the interval midpoint !
!           W(j)= (a_j + b_j)/2. The corresponding error is bounded    !
!           by WERR(j) = ABS(a_j - b_j)/2.                             !
!                                                                      !
!  WERR     [out]     REAL array, dimension (N)                        !
!           The error bound on the corresponding eigenvalue            !
!           approximation in W.                                        !
!                                                                      !
!  WL       [out]     REAL                                             !
!                                                                      !
!  WU       [out]     REAL                                             !
!           The interval (WL, WU] contains all the wanted eigenvalues. !
!           If RANGE='V', then WL=VL and WU=VU.                        !
!           If RANGE='A', then WL and WU are the global Gerschgorin    !
!                         bounds on the spectrum.                      !
!           If RANGE='I', then WL and WU are computed by DLAEBZ from   !
!                         the index range specified.                   !
!                                                                      !
!  IBLOCK   [out]     INTEGER array, dimension (N)                     !
!           At each row/column j where E(j) is zero or small, the      !
!           matrix T is considered to split into a block diagonal      !
!           matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to      !
!           which block (from 1 to the number of blocks) the           !
!           eigenvalue W(i) belongs.  (DLARRD may use the remaining    !
!           N-M elements as workspace.)                                !
!                                                                      !
!  INDEXW   [out]     INTEGER array, dimension (N)                     !
!           The indices of the eigenvalues within each block           !
!           (submatrix); for example, INDEXW(i)= j and IBLOCK(i)=k     !
!           imply that the i-th eigenvalue W(i) is the j-th eigenvalue !
!           in block k.                                                !
!                                                                      !
!  WORK     [out]     REAL array, dimension (4*N)                      !
!                                                                      !
!  IWORK    [out]     INTEGER array, dimension (3*N)                   !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           = 0:  successful exit.                                     !
!           < 0:  if INFO = -i, the i-th argument had an illegal value !
!           > 0:  some or all of the eigenvalues failed to converge or !
!                 were not computed:                                   !
!                 =1 or 3: Bisection failed to converge for some       !
!                          eigenvalues; these eigenvalues are flagged  !
!                          by a negative block number.  The effect is  !
!                          that the eigenvalues may not be as accurate !
!                          as the absolute and relative tolerances.    !
!                          This is generally caused by unexpectedly    !
!                          inaccurate arithmetic.                      !
!                 =2 or 3: RANGE='I' only: Not all of the eigenvalues  !
!                          IL:IU were found.                           !
!                          Effect: M < IU+1-IL                         !
!                          Cause:  non-monotonic arithmetic, causing   !
!                                  Sturm sequence to be non-monotonic. !
!                          Cure:   recalculate, using RANGE='A', and   !
!                                  pick out eigenvalues IL:IU. In some !
!                                  cases, increasing the PARAMETER     !
!                                  "FUDGE" may make things work.       !
!                 = 4:    RANGE='I', and the Gershgorin interval       !
!                         initially used was too small. No eigenvalues !
!                         were computed.                               !
!                         Probable cause: your machine has sloppy      !
!                                         floating-point arithmetic.   !
!                         Cure: Increase the PARAMETER "FUDGE",        !
!                               recompile, and try again.              !
!                                                                      !
!  Internal Parameters:                                                !
!                                                                      !
!   FUDGE   DOUBLE PRECISION, default = 2                              !
!           A "fudge factor" to widen the Gershgorin intervals.        !
!           Ideally, a value of 1 should work, but on machines with    !
!           sloppy arithmetic, this needs to be larger.  The default   !
!           for publicly released versions should be large enough to   !
!           handle the worst machine around.  Note that this has no    !
!           effect on accuracy of the solution.                        !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  W. Kahan, University of California, Berkeley, USA                   !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), intent(in) :: ORDER, RANGE
!
      integer,   intent(in)  :: IL, IU, N, NSPLIT
      integer,   intent(in)  :: ISPLIT(*)
      integer,   intent(out) :: INFO, M
      integer,   intent(out) :: IBLOCK(*), INDEXW(*), IWORK(*)
!
      real (dp), intent(in)  :: PIVMIN, RELTOL, VL, VU
      real (dp), intent(in)  :: D(*), E(*), E2(*), GERS(*)
      real (dp), intent(out) :: WL, WU
      real (dp), intent(out) :: W(*), WERR(*), WORK(*)
!
!  Local variable declarations.
!
      integer,   parameter :: ALLRNG = 1
      integer,   parameter :: VALRNG = 2
      integer,   parameter :: INDRNG = 3
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: HALF = 0.5_dp
      real (dp), parameter :: FUDGE = TWO
!
      logical   :: NCNVRG, TOOFEW
!
      integer   :: I, IB, IBEGIN, IDISCL, IDISCU, IE, IEND, IINFO
      integer   :: IM, IN, IOFF, IOUT, IRANGE, ITMAX, ITMP1
      integer   :: ITMP2, IW, IWOFF, J, JBLK, JDISC, JE, JEE, NB
      integer   :: NWL, NWU
      INTEGER   :: IDUMMA(1)
!
      real (dp) :: ATOLI, EPS, GL, GU, RTOLI, TMP1, TMP2
      real (dp) :: TNORM, UFLOW, WKILL, WLU, WUL
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
      M = 0
!
!  Quick return if possible.
!
      IF (N.LE.0) THEN
         RETURN
      END IF
!
!  Decode RANGE.
!
      IF (LSAME(RANGE, 'A')) THEN
         IRANGE = ALLRNG
      ELSE IF (LSAME(RANGE, 'V')) THEN
         IRANGE = VALRNG
      ELSE IF (LSAME(RANGE, 'I')) THEN
         IRANGE = INDRNG
      ELSE
         IRANGE = 0
      END IF
!
!  Check for Errors.
!
      IF (IRANGE.LE.0) THEN
         INFO = -1
      ELSE IF (.NOT.(LSAME(ORDER,'B').OR.LSAME(ORDER,'E'))) THEN
         INFO = -2
      ELSE IF (N.LT.0) THEN
         INFO = -3
      ELSE IF (IRANGE.EQ.VALRNG) THEN
         IF (VL.GE.VU) INFO = -5
      ELSE IF (IRANGE.EQ.INDRNG .AND.                                   &
     &         (IL.LT.1 .OR. IL.GT.MAX(1,N))) THEN
         INFO = -6
      ELSE IF (IRANGE.EQ.INDRNG .AND.                                   &
     &         (IU.LT.MIN(N,IL) .OR. IU.GT.N)) THEN
         INFO = -7
      END IF
!
      IF (INFO.NE.0) THEN
         RETURN
      END IF
!
!  Initialize error flags.
!
      NCNVRG = .FALSE.
      TOOFEW = .FALSE.
!
!  Simplification:
!
      IF (IRANGE.EQ.INDRNG .AND. IL.EQ.1 .AND. IU.EQ.N) IRANGE = 1
!
!  Get machine constants
!
      EPS = DLAMCH('P')
      UFLOW = DLAMCH('U')
!
!  Special Case when N=1
!  Treat case of 1x1 matrix for quick return.
!
      IF (N.EQ.1) THEN
         IF ((IRANGE.EQ.ALLRNG).OR.                                     &
     &       ((IRANGE.EQ.VALRNG).AND.(D(1).GT.VL).AND.(D(1).LE.VU)).OR. &
     &       ((IRANGE.EQ.INDRNG).AND.(IL.EQ.1).AND.(IU.EQ.1))) THEN
            M = 1
            W(1) = D(1)
!
!  The computation error of the eigenvalue is zero.
!
            WERR(1) = ZERO
            IBLOCK(1) = 1
            INDEXW(1) = 1
         END IF
         RETURN
      END IF
!
!  NB is the minimum vector length for vector bisection, or 0
!  if only scalar is to be done.
!
      NB = ILAENV(1, 'DSTEBZ', ' ', N, -1, -1, -1)
      IF (NB.LE.1) NB = 0
!
!  Find global spectral radius
!
      GL = D(1)
      GU = D(1)
      DO I = 1,N
         GL = MIN(GL, GERS(2*I - 1))
         GU = MAX(GU, GERS(2*I))
      END DO
!
!  Compute global Gerschgorin bounds and spectral diameter.
!
      TNORM = MAX(ABS(GL), ABS(GU))
      GL = GL - FUDGE*TNORM*EPS*N - FUDGE*TWO*PIVMIN
      GU = GU + FUDGE*TNORM*EPS*N + FUDGE*TWO*PIVMIN
!
!  [JAN/28/2009] remove the line below since SPDIAM variable not use
!  SPDIAM = GU - GL
!
!  Input arguments for DLAEBZ:
!  The relative tolerance.  An interval (a,b] lies within
!  "relative tolerance" if  b-a < RELTOL*MAX(|a|,|b|),
!
      RTOLI = RELTOL
!
!  Set the absolute tolerance for interval convergence to zero to force
!  interval convergence based on relative size of the interval.
!  This is dangerous because intervals might not converge when RELTOL is
!  small. But at least a very small number should be selected so that
!  for strongly graded matrices, the code can get relatively accurate
!  eigenvalues.
!
      ATOLI = FUDGE*TWO*UFLOW + FUDGE*TWO*PIVMIN

      IF (IRANGE.EQ.INDRNG) THEN
!
!  RANGE='I': Compute an interval containing eigenvalues
!  IL through IU. The initial interval [GL,GU] from the global
!  Gerschgorin bounds GL and GU is refined by DLAEBZ.
!
         ITMAX = INT((LOG(TNORM+PIVMIN)-LOG(PIVMIN)) / LOG(TWO)) + 2
         WORK(N+1) = GL
         WORK(N+2) = GL
         WORK(N+3) = GU
         WORK(N+4) = GU
         WORK(N+5) = GL
         WORK(N+6) = GU
         IWORK(1) = -1
         IWORK(2) = -1
         IWORK(3) = N + 1
         IWORK(4) = N + 1
         IWORK(5) = IL - 1
         IWORK(6) = IU
!
         CALL DLAEBZ (3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN,      &
     &                D, E, E2, IWORK(5), WORK(N+1), WORK(N+5), IOUT,   &
     &                IWORK, W, IBLOCK, IINFO)
         IF (IINFO .NE. 0) THEN
            INFO = IINFO
            RETURN
         END IF
!
!  On exit, output intervals may not be ordered by ascending negcount.
!
         IF (IWORK(6).EQ.IU) THEN
            WL = WORK(N+1)
            WLU = WORK(N+3)
            NWL = IWORK(1)
            WU = WORK(N+4)
            WUL = WORK(N+2)
            NWU = IWORK(4)
         ELSE
            WL = WORK(N+2)
            WLU = WORK(N+4)
            NWL = IWORK(2)
            WU = WORK(N+3)
            WUL = WORK(N+1)
            NWU = IWORK(3)
         END IF
!
!  On exit, the interval [WL, WLU] contains a value with negcount NWL,
!  and [WUL, WU] contains a value with negcount NWU.
!
         IF (NWL.LT.0 .OR. NWL.GE.N .OR. NWU.LT.1 .OR. NWU.GT.N) THEN
            INFO = 4
            RETURN
         END IF

      ELSE IF (IRANGE.EQ.VALRNG) THEN
         WL = VL
         WU = VU

      ELSE IF (IRANGE.EQ.ALLRNG) THEN
         WL = GL
         WU = GU
      END IF
!
!  Find Eigenvalues -- Loop Over blocks and recompute NWL and NWU.
!  NWL accumulates the number of eigenvalues .le. WL,
!  NWU accumulates the number of eigenvalues .le. WU
!
      M = 0
      IEND = 0
      INFO = 0
      NWL = 0
      NWU = 0
!
      DO 70 JBLK = 1, NSPLIT
         IOFF = IEND
         IBEGIN = IOFF + 1
         IEND = ISPLIT(JBLK)
         IN = IEND - IOFF
!
         IF (IN.EQ.1) THEN
!
!  1x1 block.
!
            IF (WL.GE.D(IBEGIN)-PIVMIN) NWL = NWL + 1
            IF (WU.GE.D(IBEGIN)-PIVMIN) NWU = NWU + 1
            IF (IRANGE.EQ.ALLRNG .OR.                                   &
     &          (WL.LT.D(IBEGIN)-PIVMIN .AND.                           &
     &           WU.GE.D(IBEGIN)-PIVMIN)) THEN
               M = M + 1
               W(M) = D(IBEGIN)
               WERR(M) = ZERO
!
!  The gap for a single block does not matter for the later
!  algorithm and is assigned an arbitrary large value.
!
               IBLOCK(M) = JBLK
               INDEXW(M) = 1
            END IF
!
!  Disabled 2x2 case because of a failure on the following matrix:
!
!        RANGE = 'I', IL = IU = 4
!          Original Tridiagonal, d = [
!           -0.150102010615740E+00
!           -0.849897989384260E+00
!           -0.128208148052635E-15
!            0.128257718286320E-15
!          ];
!          e = [
!           -0.357171383266986E+00
!           -0.180411241501588E-15
!           -0.175152352710251E-15
!          ];
!

!         ELSE IF (IN.EQ.2) THEN
!
!  2x2 block.
!
!            DISC = SQRT((HALF*(D(IBEGIN)-D(IEND)))**2 + E(IBEGIN)**2)
!            TMP1 = HALF*(D(IBEGIN)+D(IEND))
!            L1 = TMP1 - DISC
!            IF (WL.GE. L1-PIVMIN) NWL = NWL + 1
!            IF (WU.GE. L1-PIVMIN) NWU = NWU + 1
!            IF (IRANGE.EQ.ALLRNG .OR.                                  &
!     &          (WL.LT.L1-PIVMIN .AND. WU.GE. L1-PIVMIN)) THEN
!               M = M + 1
!               W(M) = L1
!
!  The uncertainty of eigenvalues of a 2x2 matrix is very small.
!
!               WERR(M) = EPS * ABS(W(M)) * TWO
!               IBLOCK(M) = JBLK
!               INDEXW(M) = 1
!            END IF
!            L2 = TMP1 + DISC
!            IF (WL.GE. L2-PIVMIN) NWL = NWL + 1
!            IF (WU.GE. L2-PIVMIN) NWU = NWU + 1
!            IF (IRANGE.EQ.ALLRNG .OR.                                  &
!     &          (WL.LT.L2-PIVMIN .AND. WU.GE.L2-PIVMIN)) THEN
!               M = M + 1
!               W(M) = L2
!
!  The uncertainty of eigenvalues of a 2x2 matrix is very small.
!
!               WERR(M) = EPS * ABS(W(M)) * TWO
!               IBLOCK(M) = JBLK
!               INDEXW(M) = 2
!            END IF

         ELSE
!
!  General Case - block of size IN >= 2:
!  Compute local Gerschgorin interval and use it as the initial
!  interval for DLAEBZ.
!
            GU = D(IBEGIN)
            GL = D(IBEGIN)
            TMP1 = ZERO

            DO J = IBEGIN, IEND
               GL = MIN(GL, GERS(2*J - 1))
               GU = MAX(GU, GERS(2*J))
            END DO
!
!  [JAN/28/2009]
!  Change SPDIAM by TNORM in lines 2 and 3 thereafter
!  line 1: remove computation of SPDIAM (not useful anymore)
!
!           SPDIAM = GU - GL
!           GL = GL - FUDGE*SPDIAM*EPS*IN - FUDGE*PIVMIN
!           GU = GU + FUDGE*SPDIAM*EPS*IN + FUDGE*PIVMIN
!
            GL = GL - FUDGE*TNORM*EPS*IN - FUDGE*PIVMIN
            GU = GU + FUDGE*TNORM*EPS*IN + FUDGE*PIVMIN
!
            IF (IRANGE.GT.1) THEN
               IF (GU.LT.WL) THEN
!
!  The local block contains none of the wanted eigenvalues.
!
                  NWL = NWL + IN
                  NWU = NWU + IN
                  GOTO 70
               END IF
!
!  Refine search interval if possible, only range (WL,WU] matters.
!
               GL = MAX(GL, WL)
               GU = MIN(GU, WU)
               IF (GL.GE.GU) GOTO 70
            END IF
!
!  Find negcount of initial interval boundaries GL and GU
!
            WORK (N+1) = GL
            WORK (N+IN+1) = GU
            CALL DLAEBZ (1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,     &
     &                   D(IBEGIN), E(IBEGIN), E2(IBEGIN),              &
     &                   IDUMMA, WORK(N+1), WORK(N+2*IN+1), IM,         &
     &                   IWORK, W(M+1), IBLOCK(M+1), IINFO)
            IF (IINFO .NE. 0) THEN
               INFO = IINFO
               RETURN
            END IF
!
            NWL = NWL + IWORK(1)
            NWU = NWU + IWORK(IN+1)
            IWOFF = M - IWORK(1)
!
!  Compute Eigenvalues.
!
            ITMAX = INT((LOG(GU-GL+PIVMIN)-LOG(PIVMIN)) / LOG(TWO)) + 2
            CALL DLAEBZ (2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI,         &
     &                   PIVMIN,                                        &
     &                   D(IBEGIN), E(IBEGIN), E2(IBEGIN),              &
     &                   IDUMMA, WORK(N+1), WORK(N+2*IN+1), IOUT,       &
     &                   IWORK, W(M+1), IBLOCK(M+1), IINFO)
            IF (IINFO .NE. 0) THEN
               INFO = IINFO
               RETURN
            END IF
!
!  Copy eigenvalues into W and IBLOCK.
!  Use -JBLK for block number for unconverged eigenvalues.
!  Loop over the number of output intervals from DLAEBZ
!
            DO J = 1, IOUT
!
!  Eigenvalue approximation is middle point of interval
!
               TMP1 = HALF*(WORK(J+N)+WORK(J+IN+N))
!
!  Semi length of error interval.
!
               TMP2 = HALF*ABS(WORK(J+N)-WORK(J+IN+N))
               IF (J.GT.IOUT-IINFO) THEN
!
!  Flag non-convergence.
!
                  NCNVRG = .TRUE.
                  IB = -JBLK
               ELSE
                  IB = JBLK
               END IF
               DO JE = IWORK(J)+1+IWOFF, IWORK(J+IN)+IWOFF
                  W(JE) = TMP1
                  WERR(JE) = TMP2
                  INDEXW(JE) = JE - IWOFF
                  IBLOCK(JE) = IB
               END DO
            END DO
!
            M = M + IM
         END IF
   70 CONTINUE
!
!  If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
!  If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
!
      IF (IRANGE.EQ.INDRNG) THEN
         IDISCL = IL - 1 - NWL
         IDISCU = NWU - IU
!
         IF (IDISCL.GT.0) THEN
            IM = 0
            DO JE = 1, M
!
!  Remove some of the smallest eigenvalues from the left so that
!  at the end IDISCL =0. Move all eigenvalues up to the left.
!
               IF (W(JE).LE.WLU .AND. IDISCL.GT.0) THEN
                  IDISCL = IDISCL - 1
               ELSE
                  IM = IM + 1
                  W(IM) = W(JE)
                  WERR(IM) = WERR(JE)
                  INDEXW(IM) = INDEXW(JE)
                  IBLOCK(IM) = IBLOCK(JE)
               END IF
            END DO
            M = IM
         END IF
         IF (IDISCU.GT.0) THEN
!
!  Remove some of the largest eigenvalues from the right so that
!  at the end IDISCU =0. Move all eigenvalues up to the left.
!
            IM=M+1
            DO JE = M, 1, -1
               IF (W(JE).GE.WUL .AND. IDISCU.GT.0) THEN
                  IDISCU = IDISCU - 1
               ELSE
                  IM = IM - 1
                  W(IM) = W(JE)
                  WERR(IM) = WERR(JE)
                  INDEXW(IM) = INDEXW(JE)
                  IBLOCK(IM) = IBLOCK(JE)
               END IF
            END DO
            JEE = 0
            DO  JE = IM, M
               JEE = JEE + 1
               W(JEE) = W(JE)
               WERR(JEE) = WERR(JE)
               INDEXW(JEE) = INDEXW(JE)
               IBLOCK(JEE) = IBLOCK(JE)
            END DO
            M = M-IM+1
         END IF

         IF (IDISCL.GT.0 .OR. IDISCU.GT.0) THEN
!
!  Code to deal with effects of bad arithmetic. (If N(w) is
!  monotone non-decreasing, this should never happen.)
!  Some low eigenvalues to be discarded are not in (WL,WLU],
!  or high eigenvalues to be discarded are not in (WUL,WU]
!  so just kill off the smallest IDISCL/largest IDISCU
!  eigenvalues, by marking the corresponding IBLOCK = 0
!
            IF (IDISCL.GT.0) THEN
               WKILL = WU
               DO JDISC = 1, IDISCL
                  IW = 0
                  DO JE = 1, M
                     IF (IBLOCK(JE).NE.0 .AND.                          &
     &                   (W(JE).LT.WKILL .OR. IW.EQ.0)) THEN
                        IW = JE
                        WKILL = W(JE)
                     END IF
                  END DO
                  IBLOCK(IW) = 0
               END DO
            END IF
            IF (IDISCU.GT.0) THEN
               WKILL = WL
               DO JDISC = 1, IDISCU
                  IW = 0
                  DO JE = 1, M
                     IF (IBLOCK(JE).NE.0 .AND.                         &
     &                   (W(JE).GE.WKILL .OR. IW.EQ.0)) THEN
                        IW = JE
                        WKILL = W(JE)
                     END IF
                  END DO
                  IBLOCK(IW) = 0
               END DO
            END IF
!
!  Now erase all eigenvalues with IBLOCK set to zero.
!
            IM = 0
            DO JE = 1, M
               IF (IBLOCK(JE).NE.0) THEN
                  IM = IM + 1
                  W(IM) = W(JE)
                  WERR(IM) = WERR(JE)
                  INDEXW(IM) = INDEXW(JE)
                  IBLOCK(IM) = IBLOCK(JE)
               END IF
            END DO
            M = IM
         END IF
         IF (IDISCL.LT.0 .OR. IDISCU.LT.0) THEN
            TOOFEW = .TRUE.
         END IF
      END IF
!
      IF ((IRANGE.EQ.ALLRNG .AND. M.NE.N) .OR.                          &
     &    (IRANGE.EQ.INDRNG .AND. M.NE.IU-IL+1)) THEN
         TOOFEW = .TRUE.
      END IF
!
!  If ORDER='B', do nothing the eigenvalues are already sorted by
!  block.
!  If ORDER='E', sort the eigenvalues from smallest to largest.
!
      IF (LSAME(ORDER,'E') .AND. NSPLIT.GT.1) THEN
         DO JE = 1, M - 1
            IE = 0
            TMP1 = W(JE)
            DO J = JE + 1, M
               IF (W(J).LT.TMP1) THEN
                  IE = J
                  TMP1 = W(J)
               END IF
            END DO
            IF (IE.NE.0) THEN
               TMP2 = WERR(IE)
               ITMP1 = IBLOCK(IE)
               ITMP2 = INDEXW(IE)
               W(IE) = W(JE)
               WERR(IE) = WERR(JE)
               IBLOCK(IE) = IBLOCK(JE)
               INDEXW(IE) = INDEXW(JE)
               W(JE) = TMP1
               WERR(JE) = TMP2
               IBLOCK(JE) = ITMP1
               INDEXW(JE) = ITMP2
            END IF
         END DO
      END IF
!
      INFO = 0
      IF (NCNVRG) INFO = INFO + 1
      IF (TOOFEW) INFO = INFO + 2
!
      RETURN
      END SUBROUTINE DLARRD
!
      SUBROUTINE DLARRE (RANGE, N, VL, VU, IL, IU, D, E, E2,            &
     &                   RTOL1, RTOL2, SPLTOL, NSPLIT, ISPLIT, M,       &
     &                   W, WERR, WGAP, IBLOCK, INDEXW, GERS, PIVMIN,   &
     &                   WORK, IWORK, INFO)
!
!=======================================================================
!                                                                      !
!  To find the desired eigenvalues of a given real symmetric           !
!  tridiagonal matrix T, DLARRE sets any "small" off-diagonal          !
!  elements to zero, and for each unreduced block T_i, it finds        !
!  (a) a suitable shift at one end of the block's spectrum,            !
!  (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and   !
!  (c) eigenvalues of each L_i D_i L_i^T.                              !
!                                                                      !
!  The representations and eigenvalues found are then used by          !
!  DSTEMR to compute the eigenvectors of T.                            !
!  The accuracy varies depending on whether bisection is used to       !
!  find a few eigenvalues or the dqds algorithm (subroutine DLASQ2)    !
!  to compute all and then discard any unwanted one.                   !
!                                                                      !
!  As an added benefit, DLARRE also outputs the n Gerschgorin          !
!  intervals for the matrices L_i D_i L_i^T.                           !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  RANGE    [in]      CHARACTER*1                                      !
!           = 'A': ("All")   all eigenvalues will be found.            !
!           = 'V': ("Value") all eigenvalues in the half-open interval !
!                            (VL, VU] will be found.                   !
!           = 'I': ("Index") the IL-th through IU-th eigenvalues (of   !
!                            the entire matrix) will be found.         !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix, N > 0.                            !
!                                                                      !
!  VL       [in,out]  REAL                                             !
!           If RANGE='V', the lower bound for the eigenvalues.         !
!           Eigenvalues less than or equal to VL, or greater than VU,  !
!           will not be returned.  VL < VU.                            !
!           If RANGE='I' or ='A', DLARRE computes bounds on the        !
!           desired part of the spectrum.                              !
!                                                                      !
!  VU       [in,out]  REAL                                             !
!           If RANGE='V', the upper bound for the eigenvalues.         !
!           Eigenvalues less than or equal to VL, or greater than VU,  !
!           will not be returned.  VL < VU.                            !
!           If RANGE='I' or ='A', DLARRE computes bounds on the        !
!           desired part of the spectrum.                              !
!                                                                      !
!  IL       [in]      INTEGER                                          !
!           If RANGE='I', the index of the smallest eigenvalue to be   !
!           returned.                                                  !
!           1 <= IL <= IU <= N.                                        !
!                                                                      !
!  IU       [in]      INTEGER                                          !
!           If RANGE='I', the index of the largest eigenvalue to be    !
!           returned.                                                  !
!           1 <= IL <= IU <= N.                                        !
!                                                                      !
!  D        [in,out]  REAL  array, dimension (N)                       !
!           On entry, the N diagonal elements of the tridiagonal       !
!           matrix T.                                                  !
!           On exit, the N diagonal elements of the diagonal           !
!           matrices D_i.                                              !
!                                                                      !
!  E        [in,out]  REAL array, dimension (N)                        !
!           On entry, the first (N-1) entries contain the subdiagonal  !
!           elements of the tridiagonal matrix T; E(N) need not be set.!
!           On exit, E contains the subdiagonal elements of the unit   !
!           bidiagonal matrices L_i. The entries E(ISPLIT(I)),         !
!           1 <= I <= NSPLIT, contain the base points sigma_i on       !
!           output.                                                    !
!                                                                      !
!  E2       [in,out]  REAL array, dimension (N)                        !
!           On entry, the first (N-1) entries contain the SQUARES of   !
!           the subdiagonal elements of the tridiagonal matrix T;      !
!           E2(N) need not be set.                                     !
!           On exit, the entries E2(ISPLIT(I)), 1 <= I <= NSPLIT,      !
!           have been set to zero.                                     !
!                                                                      !
!  RTOL1    [in]      REAL                                             !
!                                                                      !
!  RTOL2    [in]      REAL                                             !
!           Parameters for bisection.                                  !
!           An interval [LEFT,RIGHT] has converged if                  !
!           RIGHT-LEFT < MAX(RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|))     !
!                                                                      !
!  SPLTOL   [in]      REAL                                             !
!           The threshold for splitting.                               !
!                                                                      !
!  NSPLIT   [out]     INTEGER                                          !
!           The number of blocks T splits into. 1 <= NSPLIT <= N.      !
!                                                                      !
!  ISPLIT   [out]     INTEGER array, dimension (N)                     !
!           The splitting points, at which T breaks up into blocks.    !
!           The first block consists of rows/columns 1 to ISPLIT(1),   !
!           the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),  !
!           etc., and the NSPLIT-th consists of rows/columns           !
!           ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.               !
!                                                                      !
!  M        [out]     INTEGER                                          !
!           The total number of eigenvalues (of all L_i D_i L_i^T)     !
!           found.                                                     !
!                                                                      !
!  W        [out]     REAL array, dimension (N)                        !
!           The first M elements contain the eigenvalues. The          !
!           eigenvalues of each of the blocks, L_i D_i L_i^T, are      !
!           sorted in ascending order (DLARRE may use the              !
!           remaining N-M elements as workspace).                      !
!                                                                      !
!  WERR     [out]     REAL array, dimension (N)                        !
!           The error bound on the corresponding eigenvalue in W.      !
!                                                                      !
!  WGAP     [out]     REAL array, dimension (N)                        !
!           The separation from the right neighbor eigenvalue in W.    !
!           The gap is only with respect to the eigenvalues of the     !
!           same block as each block has its own representation tree.  !
!           Exception: at the right end of a block we store the left   !
!           gap.                                                       !
!                                                                      !
!  IBLOCK   [out]     INTEGER array, dimension (N)                     !
!           The indices of the blocks (submatrices) associated with    !
!           the corresponding eigenvalues in W; IBLOCK(i)=1 if         !
!           eigenvalue W(i) belongs to the first block from the top,   !
!           =2 if W(i) belongs to the second block, etc.               !
!                                                                      !
!  INDEXW   [out]     INTEGER array, dimension (N)                     !
!           The indices of the eigenvalues within each block           !
!           (submatrix); for example, INDEXW(i)= 10 and IBLOCK(i)=2    !
!           imply that the i-th eigenvalue W(i) is the 10-th           !
!           eigenvalue in block 2.                                     !
!                                                                      !
!  GERS     [out]     REAL array, dimension (2*N)                      !
!           The N Gerschgorin intervals (the i-th Gerschgorin interval !
!           is (GERS(2*i-1), GERS(2*i)).                               !
!                                                                      !
!  PIVMIN   [out]     REAL                                             !
!           The minimum pivot in the Sturm sequence for T.             !
!                                                                      !
!  WORK     [out]     REAL array, dimension (6*N)                      !
!           Workspace.                                                 !
!                                                                      !
!  IWORK    [out]     INTEGER array, dimension (5*N)                   !
!           Workspace.                                                 !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           = 0:  successful exit.                                     !
!           > 0:  A problem occurred in DLARRE.                        !
!           < 0:  One of the called subroutines signaled an internal   !
!                 problem. Needs inspection of the corresponding       !
!                 parameter IINFO for further information.             !
!                                                                      !
!           =-1:  Problem in DLARRD.                                   !
!           = 2:  No base representation could be found in MAXTRY      !
!                 iterations. Increasing MAXTRY and recompilation      !
!                 might be a remedy.                                   !
!           =-3:  Problem in DLARRB when computing the refined root    !
!                 representation for DLASQ2.                           !
!           =-4:  Problem in DLARRB when preforming bisection on the   !
!                 desired part of the spectrum.                        !
!           =-5:  Problem in DLASQ2.                                   !
!           =-6:  Problem in DLASQ2.                                   !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  The base representations are required to suffer very little         !
!  element growth and consequently define all their eigenvalues to     !
!  high relative accuracy.                                             !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=1), intent(in) :: RANGE
!
      integer,   intent(in)    :: IL, IU, N
      integer,   intent(out)   :: INFO, M, NSPLIT
      integer,   intent(out)   :: IBLOCK(*), INDEXW(*), ISPLIT(*)
      integer,   intent(out)   :: IWORK(*)
!
      real (dp), intent(in)    :: RTOL1, RTOL2, SPLTOL
      real (dp), intent(inout) :: VL, VU
      real (dp), intent(inout) :: D(*), E(*), E2(*), GERS(*)
      real (dp), intent(out)   :: PIVMIN
      real (dp), intent(out)   :: W(*), WERR(*), WGAP(*), WORK(*)
!
!  Local variable declarations.
!
      integer  , parameter :: MAXTRY = 6
      integer  , parameter :: ALLRNG = 1
      integer  , parameter :: INDRNG = 2
      integer  , parameter :: VALRNG = 3
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: FOUR = 4.0_dp
      real (dp), parameter :: HNDRD = 100.0_dp
      real (dp), parameter :: PERT = 8.0_dp
      real (dp), parameter :: HALF = 0.5_dp
      real (dp), parameter :: FOURTH = 0.25_dp
      real (dp), parameter :: FAC= HALF
      real (dp), parameter :: MAXGROWTH = 64.0_dp
      real (dp), parameter :: FUDGE = 2.0_dp
!
      logical   :: FORCEB, NOREP, USEDQD
!
      integer   :: CNT, CNT1, CNT2, I, IBEGIN, IDUM, IEND, IINFO
      integer   :: IN, INDL, INDU, IRANGE, J, JBLK, MB, MM
      integer   :: WBEGIN, WEND
      INTEGER   :: ISEED(4)
!
      real (dp) :: AVGAP, BSRTOL, CLWDTH, DMAX, DPIVOT, EABS
      real (dp) :: EMAX, EOLD, EPS, GL, GU, ISLEFT, ISRGHT, RTL
      real (dp) :: RTOL, S1, S2, SAFMIN, SGNDEF, SIGMA, SPDIAM
      real (dp) :: TAU, TMP, TMP1
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
      NSPLIT = 0
      M = 0
!
!  Quick return if possible
!
      IF (N.LE.0) THEN
         RETURN
      END IF
!
!  Decode RANGE
!
      IF (LSAME(RANGE, 'A')) THEN
         IRANGE = ALLRNG
      ELSE IF (LSAME(RANGE, 'V')) THEN
         IRANGE = VALRNG
      ELSE IF (LSAME(RANGE, 'I')) THEN
         IRANGE = INDRNG
      END IF
!
!  Get machine constants
!
      SAFMIN = DLAMCH('S')
      EPS = DLAMCH('P')
!
!  Set parameters.
!
      RTL = SQRT(EPS)
      BSRTOL = SQRT(EPS)
!
!  Treat case of 1x1 matrix for quick return.
!
      IF (N.EQ.1) THEN
         IF ((IRANGE.EQ.ALLRNG).OR.                                     &
     &       ((IRANGE.EQ.VALRNG).AND.(D(1).GT.VL).AND.(D(1).LE.VU)).OR. &
     &       ((IRANGE.EQ.INDRNG).AND.(IL.EQ.1).AND.(IU.EQ.1))) THEN
            M = 1
            W(1) = D(1)
!
!  The computation error of the eigenvalue is zero.
!
            WERR(1) = ZERO
            WGAP(1) = ZERO
            IBLOCK(1) = 1
            INDEXW(1) = 1
            GERS(1) = D(1)
            GERS(2) = D(1)
         END IF
!
!  Store the shift for the initial RRR, which is zero in this case
!
         E(1) = ZERO
         RETURN
      END IF
!
!  General case: tridiagonal matrix of order > 1
!
!  Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter.
!  Compute maximum off-diagonal entry and pivmin.
!
      GL = D(1)
      GU = D(1)
      EOLD = ZERO
      EMAX = ZERO
      E(N) = ZERO
      DO I = 1,N
         WERR(I) = ZERO
         WGAP(I) = ZERO
         EABS = ABS(E(I))
         IF (EABS .GE. EMAX) THEN
            EMAX = EABS
         END IF
         TMP1 = EABS + EOLD
         GERS(2*I-1) = D(I) - TMP1
         GL =  MIN(GL, GERS(2*I - 1))
         GERS(2*I) = D(I) + TMP1
         GU = MAX(GU, GERS(2*I))
         EOLD  = EABS
      END DO
!
!  The minimum pivot allowed in the Sturm sequence for T.
!
      PIVMIN = SAFMIN * MAX(ONE, EMAX**2)
!
!  Compute spectral diameter. The Gerschgorin bounds give an
!  estimate that is wrong by at most a factor of SQRT(2).
!
      SPDIAM = GU - GL
!
!  Compute splitting points.
!
      CALL DLARRA (N, D, E, E2, SPLTOL, SPDIAM,                         &
     &             NSPLIT, ISPLIT, IINFO)
!
!  Can force use of bisection instead of faster DQDS.
!  Option left in the code for future multisection work.
!
      FORCEB = .FALSE.
!
!  Initialize USEDQD, DQDS should be used for ALLRNG unless someone
!  explicitly wants bisection.
!
      USEDQD = ((IRANGE.EQ.ALLRNG) .AND. (.NOT.FORCEB))

      IF ((IRANGE.EQ.ALLRNG) .AND. (.NOT. FORCEB)) THEN
!
!  Set interval [VL,VU] that contains all eigenvalues.
!
         VL = GL
         VU = GU
      ELSE
!
!  We call DLARRD to find crude approximations to the eigenvalues
!  in the desired range. In case IRANGE = INDRNG, we also obtain the
!  interval (VL,VU] that contains all the wanted eigenvalues.
!  An interval [LEFT,RIGHT] has converged if
!  RIGHT-LEFT.LT.RTOL*MAX(ABS(LEFT),ABS(RIGHT))
!  DLARRD needs a WORK of size 4*N, IWORK of size 3*N
!
         CALL DLARRD (RANGE, 'B', N, VL, VU, IL, IU, GERS,              &
     &                BSRTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT,         &
     &                MM, W, WERR, VL, VU, IBLOCK, INDEXW,              &
     &                WORK, IWORK, IINFO)
         IF (IINFO.NE.0) THEN
            INFO = -1
            RETURN
         END IF
!
!  Make sure that the entries M+1 to N in W, WERR, IBLOCK, INDEXW are 0.
!
         DO I = MM+1,N
            W(I) = ZERO
            WERR(I) = ZERO
            IBLOCK(I) = 0
            INDEXW(I) = 0
         END DO
      END IF
!
!  Loop over unreduced blocks.
!
      IBEGIN = 1
      WBEGIN = 1
      DO 170 JBLK = 1, NSPLIT
         IEND = ISPLIT(JBLK)
         IN = IEND - IBEGIN + 1
!
!  1 X 1.
!
         IF (IN.EQ.1) THEN
            IF ((IRANGE.EQ.ALLRNG).OR.((IRANGE.EQ.VALRNG).AND.          &
     &          (D(IBEGIN).GT.VL).AND.(D(IBEGIN).LE.VU)).OR.            &
     &          ((IRANGE.EQ.INDRNG).AND.(IBLOCK(WBEGIN).EQ.JBLK))) THEN
               M = M + 1
               W(M) = D(IBEGIN)
               WERR(M) = ZERO
!
!  The gap for a single block does not matter for the later
!  algorithm and is assigned an arbitrary large value
!
               WGAP(M) = ZERO
               IBLOCK(M) = JBLK
               INDEXW(M) = 1
               WBEGIN = WBEGIN + 1
            END IF
!
!  E(IEND) holds the shift for the initial RRR.
!
            E(IEND) = ZERO
            IBEGIN = IEND + 1
            GOTO 170
         END IF
!
!  Blocks of size larger than 1x1.
!
!  E(IEND) will hold the shift for the initial RRR, for now set it =0
!
         E(IEND) = ZERO
!
!  Find local outer bounds GL,GU for the block
!
         GL = D(IBEGIN)
         GU = D(IBEGIN)
         DO I = IBEGIN , IEND
            GL = MIN(GERS(2*I-1), GL)
            GU = MAX(GERS(2*I), GU)
         END DO
         SPDIAM = GU - GL

         IF (.NOT. ((IRANGE.EQ.ALLRNG).AND.(.NOT.FORCEB))) THEN
!
!  Count the number of eigenvalues in the current block.
!
            MB = 0
            DO I = WBEGIN,MM
               IF (IBLOCK(I).EQ.JBLK) THEN
                  MB = MB+1
               ELSE
                  GOTO 21
               END IF
            END DO
 21         CONTINUE

            IF (MB.EQ.0) THEN
!
!  No eigenvalue in the current block lies in the desired range
!  E(IEND) holds the shift for the initial RRR
!
               E(IEND) = ZERO
               IBEGIN = IEND + 1
               GOTO 170
            ELSE
!
!  Decide whether dqds or bisection is more efficient.
!
               USEDQD = ((MB .GT. FAC*IN) .AND. (.NOT.FORCEB))
               WEND = WBEGIN + MB - 1
!
!  Calculate gaps for the current block
!  In later stages, when representations for individual
!  eigenvalues are different, we use SIGMA = E(IEND).
!
               SIGMA = ZERO
               DO I = WBEGIN, WEND - 1
                  WGAP(I) = MAX(ZERO,                                   &
     &                          W(I+1)-WERR(I+1) - (W(I)+WERR(I)))
               END DO
               WGAP(WEND) = MAX(ZERO,                                   &
     &                          VU - SIGMA - (W(WEND)+WERR(WEND)))
!
!  Find local index of the first and last desired evalue.
!
               INDL = INDEXW(WBEGIN)
               INDU = INDEXW(WEND)
            END IF
         END IF
         IF (((IRANGE.EQ.ALLRNG) .AND. (.NOT. FORCEB)).OR.USEDQD) THEN
!
!  Case of DQDS.
!  Find approximations to the extremal eigenvalues of the block
!
            CALL DLARRK (IN, 1, GL, GU, D(IBEGIN),                      &
     &                   E2(IBEGIN), PIVMIN, RTL, TMP, TMP1, IINFO)
            IF (IINFO.NE.0) THEN
               INFO = -1
               RETURN
            END IF
            ISLEFT = MAX(GL, TMP - TMP1 - HNDRD * EPS* ABS(TMP - TMP1))

            CALL DLARRK (IN, IN, GL, GU, D(IBEGIN),                     &
     &                   E2(IBEGIN), PIVMIN, RTL, TMP, TMP1, IINFO)
            IF (IINFO.NE.0) THEN
               INFO = -1
               RETURN
            END IF
            ISRGHT = MIN(GU, TMP + TMP1 + HNDRD * EPS * ABS(TMP + TMP1))
!
!  Improve the estimate of the spectral diameter.
!
            SPDIAM = ISRGHT - ISLEFT
         ELSE
!
!  Case of bisection.
!  Find approximations to the wanted extremal eigenvalues.
!
            ISLEFT = MAX(GL, W(WBEGIN) - WERR(WBEGIN) -                 &
     &                   HNDRD * EPS*ABS(W(WBEGIN)- WERR(WBEGIN)))
            ISRGHT = MIN(GU,W(WEND) + WERR(WEND) +                      &
     &                   HNDRD * EPS * ABS(W(WEND)+ WERR(WEND)))
         END IF
!
!  Decide whether the base representation for the current block
!  L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I
!  should be on the left or the right end of the current block.
!  The strategy is to shift to the end which is "more populated"
!  Furthermore, decide whether to use DQDS for the computation of
!  the eigenvalue approximations at the end of DLARRE or bisection.
!  dqds is chosen if all eigenvalues are desired or the number of
!  eigenvalues to be computed is large compared to the blocksize.
!
         IF ((IRANGE.EQ.ALLRNG) .AND. (.NOT.FORCEB)) THEN
!
!  If all the eigenvalues have to be computed, we use dqd.
!
            USEDQD = .TRUE.
!
!  INDL is the local index of the first eigenvalue to compute
!
            INDL = 1
            INDU = IN
!
!  MB =  number of eigenvalues to compute.
!
            MB = IN
            WEND = WBEGIN + MB - 1
!
!  Define 1/4 and 3/4 points of the spectrum
!
            S1 = ISLEFT + FOURTH * SPDIAM
            S2 = ISRGHT - FOURTH * SPDIAM
         ELSE
!
!  DLARRD has computed IBLOCK and INDEXW for each eigenvalue
!  approximation.
!  Choose sigma.
!
            IF (USEDQD) THEN
               S1 = ISLEFT + FOURTH * SPDIAM
               S2 = ISRGHT - FOURTH * SPDIAM
            ELSE
               TMP = MIN(ISRGHT,VU) -  MAX(ISLEFT,VL)
               S1 =  MAX(ISLEFT,VL) + FOURTH * TMP
               S2 =  MIN(ISRGHT,VU) - FOURTH * TMP
            END IF
         END IF
!
!  Compute the negcount at the 1/4 and 3/4 points.
!
         IF (MB.GT.1) THEN
            CALL DLARRC ('T', IN, S1, S2, D(IBEGIN),                    &
     &                   E(IBEGIN), PIVMIN, CNT, CNT1, CNT2, IINFO)
         END IF

         IF (MB.EQ.1) THEN
            SIGMA = GL
            SGNDEF = ONE
         ELSE IF (CNT1 - INDL .GE. INDU - CNT2) THEN
            IF ((IRANGE.EQ.ALLRNG) .AND. (.NOT.FORCEB)) THEN
               SIGMA = MAX(ISLEFT,GL)
            ELSE IF (USEDQD) THEN
!
!  Use Gerschgorin bound as shift to get pos def matrix for dqds.
!
               SIGMA = ISLEFT
            ELSE
!
!  Use approximation of the first desired eigenvalue of the block as
!  shift.
!
               SIGMA = MAX(ISLEFT,VL)
            END IF
            SGNDEF = ONE
         ELSE
            IF ((IRANGE.EQ.ALLRNG) .AND. (.NOT.FORCEB)) THEN
               SIGMA = MIN(ISRGHT,GU)
            ELSE IF (USEDQD) THEN
!
!  Use Gerschgorin bound as shift to get neg def matrixe for dqds.
!
               SIGMA = ISRGHT
            ELSE
!
!  Use approximation of the first desired eigenvalue of the block as
!  shift.
!
               SIGMA = MIN(ISRGHT,VU)
            END IF
            SGNDEF = -ONE
         END IF
!
!  An initial SIGMA has been chosen that will be used for computing
!  T - SIGMA I = L D L^T
!  Define the increment TAU of the shift in case the initial shift
!  needs to be refined to obtain a factorization with not too much
!  element growth.
!
         IF (USEDQD) THEN
!
!  The initial SIGMA was to the outer end of the spectrum
!  the matrix is definite and we need not retreat.
!
            TAU = SPDIAM*EPS*N + TWO*PIVMIN
            TAU = MAX(TAU,TWO*EPS*ABS(SIGMA))
         ELSE
            IF (MB.GT.1) THEN
               CLWDTH = W(WEND) + WERR(WEND) - W(WBEGIN) - WERR(WBEGIN)
               AVGAP = ABS(CLWDTH / DBLE(WEND-WBEGIN))
               IF (SGNDEF.EQ.ONE) THEN
                  TAU = HALF*MAX(WGAP(WBEGIN),AVGAP)
                  TAU = MAX(TAU,WERR(WBEGIN))
               ELSE
                  TAU = HALF*MAX(WGAP(WEND-1),AVGAP)
                  TAU = MAX(TAU,WERR(WEND))
               END IF
            ELSE
               TAU = WERR(WBEGIN)
            END IF
         END IF
!
         DO 80 IDUM = 1, MAXTRY
!
!  Compute L D L^T factorization of tridiagonal matrix T - sigma I.
!  Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of
!  pivots in WORK(2*IN+1:3*IN).
!
            DPIVOT = D(IBEGIN) - SIGMA
            WORK(1) = DPIVOT
            DMAX = ABS(WORK(1))
            J = IBEGIN
            DO I = 1, IN - 1
               WORK(2*IN+I) = ONE / WORK(I)
               TMP = E(J)*WORK(2*IN+I)
               WORK(IN+I) = TMP
               DPIVOT = (D(J+1)-SIGMA) - TMP*E(J)
               WORK(I+1) = DPIVOT
               DMAX = MAX(DMAX, ABS(DPIVOT))
               J = J + 1
            END DO
!
!  Check for element growth.
!
            IF (DMAX .GT. MAXGROWTH*SPDIAM) THEN
               NOREP = .TRUE.
            ELSE
               NOREP = .FALSE.
            END IF
            IF (USEDQD .AND. .NOT.NOREP) THEN
!
!  Ensure the definiteness of the representation
!  All entries of D (of L D L^T) must have the same sign
!
               DO I = 1, IN
                  TMP = SGNDEF*WORK(I)
                  IF (TMP.LT.ZERO) NOREP = .TRUE.
               END DO
            END IF
            IF (NOREP) THEN
!
!  Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin
!  shift which makes the matrix definite. So we should end up
!  here really only in the case of IRANGE = VALRNG or INDRNG.
!
               IF (IDUM.EQ.MAXTRY-1) THEN
                  IF (SGNDEF.EQ.ONE) THEN
!
!  The fudged Gerschgorin shift should succeed
!
                     SIGMA = GL - FUDGE*SPDIAM*EPS*N - FUDGE*TWO*PIVMIN
                  ELSE
                     SIGMA = GU + FUDGE*SPDIAM*EPS*N + FUDGE*TWO*PIVMIN
                  END IF
               ELSE
                  SIGMA = SIGMA - SGNDEF * TAU
                  TAU = TWO * TAU
               END IF
            ELSE
!
!  An initial RRR is found.
!
               GOTO 83
            END IF
 80      CONTINUE
!
!  If the program reaches this point, no base representation could be
!  found in MAXTRY iterations.
!
         INFO = 2
         RETURN

 83      CONTINUE
!
!  At this point, we have found an initial base representation
!  T - SIGMA I = L D L^T with not too much element growth.
!  Store the shift.
!
         E(IEND) = SIGMA
!
!  Store D and L.
!
         CALL DCOPY (IN, WORK, 1, D(IBEGIN), 1)
         CALL DCOPY (IN-1, WORK(IN+1), 1, E(IBEGIN), 1)

         IF (MB.GT.1) THEN
!
!  Perturb each entry of the base representation by a small
!  (but random) relative amount to overcome difficulties with
!   glued matrices.
!
            DO I = 1, 4
               ISEED(I) = 1
            END DO

            CALL DLARNV (2, ISEED, 2*IN-1, WORK(1))
            DO I = 1,IN-1
               D(IBEGIN+I-1) = D(IBEGIN+I-1)*(ONE+EPS*PERT*WORK(I))
               E(IBEGIN+I-1) = E(IBEGIN+I-1)*(ONE+EPS*PERT*WORK(IN+I))
            END DO
            D(IEND) = D(IEND)*(ONE+EPS*FOUR*WORK(IN))
!
         END IF
!
!  Do not update the Gerschgorin intervals because keeping track
!  of the updates would be too much work in DLARRV.
!  We update W instead and use it to locate the proper Gerschgorin
!  intervals.
!
!  Compute the required eigenvalues of L D L' by bisection or dqds
!
         IF (.NOT.USEDQD) THEN
!
!  If DLARRD has been used, shift the eigenvalue approximations
!  according to their representation. This is necessary for
!  a uniform DLARRV since dqds computes eigenvalues of the
!  shifted representation. In DLARRV, W will always hold the
!  UNshifted eigenvalue approximation.
!
            DO J=WBEGIN,WEND
               W(J) = W(J) - SIGMA
               WERR(J) = WERR(J) + ABS(W(J)) * EPS
            END DO
!
!  Use DLARRB to reduce eigenvalue error of the approximations
!  from DLARRD.
!
            DO I = IBEGIN, IEND-1
               WORK(I) = D(I) * E(I)**2
            END DO
!
!  Use bisection to find EV from INDL to INDU.
!
            CALL DLARRB (IN, D(IBEGIN), WORK(IBEGIN),                   &
     &                   INDL, INDU, RTOL1, RTOL2, INDL-1,              &
     &                   W(WBEGIN), WGAP(WBEGIN), WERR(WBEGIN),         &
     &                   WORK(2*N+1), IWORK, PIVMIN, SPDIAM,            &
     &                   IN, IINFO)
            IF (IINFO .NE. 0) THEN
               INFO = -4
               RETURN
            END IF
!
!  DLARRB computes all gaps correctly except for the last one
!  Record distance to VU/GU.
!
            WGAP(WEND) = MAX(ZERO, (VU-SIGMA) - (W(WEND) + WERR(WEND)))
            DO I = INDL, INDU
               M = M + 1
               IBLOCK(M) = JBLK
               INDEXW(M) = I
            END DO
         ELSE
!
!  Use dqds to get all eigs (and then possibly delete unwanted
!  eigenvalues).
!  Note that dqds finds the eigenvalues of the L D L^T representation
!  of T to high relative accuracy. High relative accuracy
!  might be lost when the shift of the RRR is subtracted to obtain
!  the eigenvalues of T. However, T is not guaranteed to define its
!  eigenvalues to high relative accuracy anyway.
!  Set RTOL to the order of the tolerance used in DLASQ2
!  This is an ESTIMATED error, the worst case bound is 4*N*EPS
!  which is usually too large and requires unnecessary work to be
!  done by bisection when computing the eigenvectors
!
            RTOL = LOG(DBLE(IN)) * FOUR * EPS
            J = IBEGIN
            DO I = 1, IN - 1
               WORK(2*I-1) = ABS(D(J))
               WORK(2*I) = E(J)*E(J)*WORK(2*I-1)
               J = J + 1
            END DO
            WORK(2*IN-1) = ABS(D(IEND))
            WORK(2*IN) = ZERO
            CALL DLASQ2 (IN, WORK, IINFO)
            IF (IINFO .NE. 0) THEN
!
!  If IINFO = -5 then an index is part of a tight cluster
!  and should be changed. The index is in IWORK(1) and the
!  gap is in WORK(N+1)
!
               INFO = -5
               RETURN
            ELSE
!
!  Test that all eigenvalues are positive as expected.
!
               DO I = 1, IN
                  IF (WORK(I).LT.ZERO) THEN
                     INFO = -6
                     RETURN
                  END IF
               END DO
            END IF
            IF (SGNDEF.GT.ZERO) THEN
               DO I = INDL, INDU
                  M = M + 1
                  W(M) = WORK(IN-I+1)
                  IBLOCK(M) = JBLK
                  INDEXW(M) = I
               END DO
            ELSE
               DO I = INDL, INDU
                  M = M + 1
                  W(M) = -WORK(I)
                  IBLOCK(M) = JBLK
                  INDEXW(M) = I
               END DO
            END IF
!
!  The value of RTOL below should be the tolerance in DLASQ2.
!
            DO I = M - MB + 1, M
               WERR(I) = RTOL * ABS(W(I))
            END DO
!
!  Compute the right gap between the intervals.
!
            DO I = M - MB + 1, M - 1
               WGAP(I) = MAX(ZERO, W(I+1)-WERR(I+1) - (W(I)+WERR(I)))
            END DO
            WGAP(M) = MAX(ZERO, (VU-SIGMA) - (W(M) + WERR(M)))
         END IF
!
!  Oroceed with next block.
!
         IBEGIN = IEND + 1
         WBEGIN = WEND + 1
 170  CONTINUE
!
      RETURN
      END SUBROUTINE DLARRE
!
      SUBROUTINE DLARRF (N, D, L, LD, CLSTRT, CLEND,                    &
     &                   W, WGAP, WERR,                                 &
     &                   SPDIAM, CLGAPL, CLGAPR, PIVMIN, SIGMA,         &
     &                   DPLUS, LPLUS, WORK, INFO)
!
!=======================================================================
!                                                                      !
!  Given the initial representation L D L^T and its cluster of close   !
!  eigenvalues (in a relative measure), W(CLSTRT), W(CLSTRT+1), ...    !
!  W(CLEND), DLARRF finds a new relatively robust representation       !
!  L D L^T - SIGMA I = L(+) D(+) L(+)^T such that at least one of the  !
!  eigenvalues of L(+) D(+) L(+)^T is relatively isolated.             !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix (subblock, if the matrix split).   !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           The N diagonal elements of the diagonal matrix D.          !
!                                                                      !
!  L        [in]      REAL array, dimension (N-1)                      !
!           The (N-1) subdiagonal elements of the unit bidiagonal      !
!           matrix L.                                                  !
!                                                                      !
!  LD       [in]      REAL array, dimension (N-1)                      !
!           The (N-1) elements L(i)*D(i).                              !
!                                                                      !
!  CLSTRT   [in]      INTEGER                                          !
!           The index of the first eigenvalue in the cluster.          !
!                                                                      !
!  CLEND    [in]      INTEGER                                          !
!           The index of the last eigenvalue in the cluster.           !
!                                                                      !
!  W        [in]      REAL array, dimension  >=  (CLEND-CLSTRT+1)      !
!           The eigenvalue APPROXIMATIONS of L D L^T in ascending      !
!           order. W(CLSTRT) through W(CLEND) form the cluster of      !
!           relatively close eigenalues.                               !
!                                                                      !
!  WGAP     [in,out]  REAL array, dimension  >=  (CLEND-CLSTRT+1)      !
!           The separation from the right neighbor eigenvalue in W.    !
!                                                                      !
!  WERR     [in]      REAL array, dimension  >=  (CLEND-CLSTRT+1)      !
!           WERR contain the semiwidth of the uncertainty interval     !
!           of the corresponding eigenvalue APPROXIMATION in W.        !
!                                                                      !
!  SPDIAM   [in]      REAL                                             !
!           Estimate of the spectral diameter obtained from the        !
!           Gerschgorin intervals                                      !
!                                                                      !
!  CLGAPL   [in]      REAL                                             !
!                                                                      !
!  CLGAPR   [in]      REAL                                             !
!           Absolute gap on each end of the cluster.                   !
!           Set by the calling routine to protect against shifts too   !
!           close to eigenvalues outside the cluster.                  !
!                                                                      !
!  PIVMIN   [in]      REAL                                             !
!           The minimum pivot allowed in the Sturm sequence.           !
!                                                                      !
!  SIGMA    [out]     REAL                                             !
!           The shift used to form L(+) D(+) L(+)^T.                   !
!                                                                      !
!  DPLUS    [out]     REAL array, dimension (N)                        !
!           The N diagonal elements of the diagonal matrix D(+).       !
!                                                                      !
!  LPLUS    [out]     REAL array, dimension (N-1)                      !
!           The first (N-1) elements of LPLUS contain the subdiagonal  !
!           elements of the unit bidiagonal matrix L(+).               !
!                                                                      !
!  WORK     [out]     REAL array, dimension (2*N)                      !
!           Workspace.                                                 !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           Signals processing OK (=0) or failure (=1).                !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: CLSTRT, CLEND, N
      integer,   intent(out)   :: INFO
!
      real (dp), intent(in)    :: CLGAPL, CLGAPR, PIVMIN, SPDIAM
      real (dp), intent(in)    :: D(*), L(*), LD(*), W(*), WERR(*)
      real (dp), intent(inout) :: WGAP(*)
      real (dp), intent(out)   :: SIGMA
      real (dp), intent(out)   :: DPLUS(*), LPLUS(*), WORK(*)
!
!  Local variable declarations.
!
      integer,   parameter :: KTRYMAX = 1
      integer,   parameter :: SLEFT = 1
      integer,   parameter :: SRIGHT = 2
!
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: FOUR = 4.0_dp
      real (dp), parameter :: QUART = 0.25_dp
      real (dp), parameter :: MAXGROWTH1 = 8.0_dp
      real (dp), parameter :: MAXGROWTH2 = 8.0_dp
!
      logical   :: DORRR1, FORCER, NOFAIL, SAWNAN1, SAWNAN2, TRYRRR1
!
      integer   :: I, INDX, KTRY, SHIFT
!
      real (dp) :: AVGAP, BESTSHIFT, CLWDTH, EPS, FACT, FAIL
      real (dp) :: FAIL2, GROWTHBOUND, LDELTA, LDMAX, LSIGMA
      real (dp) :: MAX1, MAX2, MINGAP, OLDP, PROD, RDELTA, RDMAX
      real (dp) :: RRR1, RRR2, RSIGMA, S, SMLGROWTH, TMP, ZNM2
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
!
!  Quick return if possible.
!
      IF (N.LE.0) THEN
         RETURN
      END IF
!
      FACT = DBLE(2**KTRYMAX)
      EPS = DLAMCH('Precision')
      SHIFT = 0
      FORCER = .FALSE.
!
!  Note that we cannot guarantee that for any of the shifts tried,
!  the factorization has a small or even moderate element growth.
!  There could be Ritz values at both ends of the cluster and despite
!  backing off, there are examples where all factorizations tried
!  (in IEEE mode, allowing zero pivots & infinities) have INFINITE
!  element growth.
!  For this reason, we should use PIVMIN in this subroutine so that at
!  least the L D L^T factorization exists. It can be checked afterwards
!  whether the element growth caused bad residuals/orthogonality.
!
!  Decide whether the code should accept the best among all
!  representations despite large element growth or signal INFO=1
!  Setting NOFAIL to .FALSE. for quick fix for bug 113,
!
      NOFAIL = .FALSE.
!
!  Compute the average gap length of the cluster.
!
      CLWDTH = ABS(W(CLEND)-W(CLSTRT)) + WERR(CLEND) + WERR(CLSTRT)
      AVGAP = CLWDTH / DBLE(CLEND-CLSTRT)
      MINGAP = MIN(CLGAPL, CLGAPR)
!
!  Initial values for shifts to both ends of cluster.
!
      LSIGMA = MIN(W(CLSTRT),W(CLEND)) - WERR(CLSTRT)
      RSIGMA = MAX(W(CLSTRT),W(CLEND)) + WERR(CLEND)
!
!  Use a small fudge to make sure that we really shift to the outside.
!
      LSIGMA = LSIGMA - ABS(LSIGMA)* FOUR * EPS
      RSIGMA = RSIGMA + ABS(RSIGMA)* FOUR * EPS
!
!  Compute upper bounds for how much to back off the initial shifts.
!
      LDMAX = QUART * MINGAP + TWO * PIVMIN
      RDMAX = QUART * MINGAP + TWO * PIVMIN

      LDELTA = MAX(AVGAP,WGAP(CLSTRT))/FACT
      RDELTA = MAX(AVGAP,WGAP(CLEND-1))/FACT
!
!  Initialize the record of the best representation found
!
      S = DLAMCH('S')
      SMLGROWTH = ONE / S
      FAIL = DBLE(N-1)*MINGAP/(SPDIAM*EPS)
      FAIL2 = DBLE(N-1)*MINGAP/(SPDIAM*SQRT(EPS))
      BESTSHIFT = LSIGMA
!
!     WHILE (KTRY <= KTRYMAX)
!
      KTRY = 0
      GROWTHBOUND = MAXGROWTH1*SPDIAM

 5    CONTINUE
      SAWNAN1 = .FALSE.
      SAWNAN2 = .FALSE.
!
!  Ensure that we do not back off too much of the initial shifts
!
      LDELTA = MIN(LDMAX,LDELTA)
      RDELTA = MIN(RDMAX,RDELTA)
!
!  Compute the element growth when shifting to both ends of the cluster
!  accept the shift if there is no element growth at one of the two ends
!
!  Left end
!
      S = -LSIGMA
      DPLUS(1) = D(1) + S
      IF (ABS(DPLUS(1)).LT.PIVMIN) THEN
         DPLUS(1) = -PIVMIN
!
!  Need to set SAWNAN1 because refined RRR test should not be used
!  in this case.
!
         SAWNAN1 = .TRUE.
      END IF
      MAX1 = ABS(DPLUS(1))
      DO I = 1, N - 1
         LPLUS(I) = LD(I) / DPLUS(I)
         S = S*LPLUS(I)*L(I) - LSIGMA
         DPLUS(I+1) = D(I+1) + S
         IF (ABS(DPLUS(I+1)).LT.PIVMIN) THEN
            DPLUS(I+1) = -PIVMIN
!
!  Need to set SAWNAN1 because refined RRR test should not be used
!  in this case.
!
            SAWNAN1 = .TRUE.
         END IF
         MAX1 = MAX(MAX1,ABS(DPLUS(I+1)))
      END DO
      SAWNAN1 = SAWNAN1 .OR.  DISNAN(MAX1)

      IF (FORCER .OR. (MAX1.LE.GROWTHBOUND .AND. .NOT.SAWNAN1)) THEN
         SIGMA = LSIGMA
         SHIFT = SLEFT
         GOTO 100
      END IF
!
!  Right end.
!
      S = -RSIGMA
      WORK(1) = D(1) + S
      IF (ABS(WORK(1)).LT.PIVMIN) THEN
         WORK(1) = -PIVMIN
!
!  Need to set SAWNAN2 because refined RRR test should not be used
!  in this case.
!
         SAWNAN2 = .TRUE.
      END IF
      MAX2 = ABS(WORK(1))
      DO I = 1, N - 1
         WORK(N+I) = LD(I) / WORK(I)
         S = S*WORK(N+I)*L(I) - RSIGMA
         WORK(I+1) = D(I+1) + S
         IF (ABS(WORK(I+1)).LT.PIVMIN) THEN
            WORK(I+1) = -PIVMIN
!
!  Need to set SAWNAN2 because refined RRR test should not be used
!  in this case.
!
            SAWNAN2 = .TRUE.
         END IF
         MAX2 = MAX(MAX2,ABS(WORK(I+1)))
      END DO
      SAWNAN2 = SAWNAN2 .OR.  DISNAN(MAX2)

      IF (FORCER .OR. (MAX2.LE.GROWTHBOUND .AND. .NOT.SAWNAN2)) THEN
         SIGMA = RSIGMA
         SHIFT = SRIGHT
         GOTO 100
      END IF
!
!  If we are at this point, both shifts led to too much element growth.
!  Record the better of the two shifts (provided it didn't lead to NaN)
!
      IF (SAWNAN1.AND.SAWNAN2) THEN
         GOTO 50                    ! both MAX1 and MAX2 are NaN

      ELSE
         IF (.NOT.SAWNAN1) THEN
            INDX = 1
            IF (MAX1.LE.SMLGROWTH) THEN
               SMLGROWTH = MAX1
               BESTSHIFT = LSIGMA
            END IF
         END IF
         IF (.NOT.SAWNAN2) THEN
            IF (SAWNAN1 .OR. MAX2.LE.MAX1) INDX = 2
            IF (MAX2.LE.SMLGROWTH) THEN
               SMLGROWTH = MAX2
               BESTSHIFT = RSIGMA
            END IF
         END IF
      END IF
!
!  If we are here, both the left and the right shift led to
!  element growth. If the element growth is moderate, then
!  we may still accept the representation, if it passes a
!  refined test for RRR. This test supposes that no NaN occurred.
!  Moreover, we use the refined RRR test only for isolated clusters.
!
      IF ((CLWDTH.LT.MINGAP/REAL(128,dp)) .AND.                         &
     &    (MIN(MAX1,MAX2).LT.FAIL2) .AND.                               &
     &    (.NOT.SAWNAN1).AND.(.NOT.SAWNAN2)) THEN
         DORRR1 = .TRUE.
      ELSE
         DORRR1 = .FALSE.
      END IF
      TRYRRR1 = .TRUE.
      IF (TRYRRR1 .AND. DORRR1) THEN
        IF (INDX.EQ.1) THEN
         TMP = ABS(DPLUS(N))
         ZNM2 = ONE
         PROD = ONE
         OLDP = ONE
         DO I = N-1, 1, -1
            IF (PROD .LE. EPS) THEN
               PROD = ((DPLUS(I+1)*WORK(N+I+1))/                        &
     &                 (DPLUS(I )*WORK(N+I )))*OLDP
            ELSE
               PROD = PROD*ABS(WORK(N+I))
            END IF
            OLDP = PROD
            ZNM2 = ZNM2 + PROD**2
            TMP = MAX(TMP, ABS(DPLUS(I) * PROD))
         END DO
         RRR1 = TMP/(SPDIAM * SQRT(ZNM2))
         IF (RRR1.LE.MAXGROWTH2) THEN
            SIGMA = LSIGMA
            SHIFT = SLEFT
            GOTO 100
         END IF
      ELSE IF (INDX.EQ.2) THEN
         TMP = ABS(WORK(N))
         ZNM2 = ONE
         PROD = ONE
         OLDP = ONE
         DO I = N-1, 1, -1
            IF (PROD .LE. EPS) THEN
               PROD = ((WORK(I+1)*LPLUS(I+1))/(WORK(I)*LPLUS(I)))*OLDP
            ELSE
               PROD = PROD*ABS(LPLUS(I))
            END IF
            OLDP = PROD
            ZNM2 = ZNM2 + PROD**2
            TMP = MAX(TMP, ABS(WORK(I) * PROD))
         END DO
         RRR2 = TMP/(SPDIAM * SQRT(ZNM2))
         IF (RRR2.LE.MAXGROWTH2) THEN
            SIGMA = RSIGMA
            SHIFT = SRIGHT
            GOTO 100
         END IF
      END IF
      END IF

 50   CONTINUE

      IF (KTRY.LT.KTRYMAX) THEN
!
!  If we are here, both shifts failed also the RRR test.
!  Back off to the outside
!
         LSIGMA = MAX(LSIGMA - LDELTA, LSIGMA - LDMAX)
         RSIGMA = MIN(RSIGMA + RDELTA, RSIGMA + RDMAX)
         LDELTA = TWO * LDELTA
         RDELTA = TWO * RDELTA
         KTRY = KTRY + 1
         GOTO 5
      ELSE
!
!  None of the representations investigated satisfied our
!  criteria. Take the best one we found.
!
         IF ((SMLGROWTH.LT.FAIL*SPDIAM).OR.NOFAIL) THEN
            LSIGMA = BESTSHIFT
            RSIGMA = BESTSHIFT
            FORCER = .TRUE.
            GOTO 5
         ELSE
            INFO = 1
            RETURN
         END IF
      END IF

 100  CONTINUE

      IF (SHIFT.EQ.SLEFT) THEN
      ELSE IF (SHIFT.EQ.SRIGHT) THEN
!
!  Store new L and D back into DPLUS, LPLUS.
!
         CALL DCOPY (N, WORK, 1, DPLUS, 1)
         CALL DCOPY (N-1, WORK(N+1), 1, LPLUS, 1)
      END IF
!
      RETURN
      END SUBROUTINE DLARRF
!
      SUBROUTINE DLARRJ (N, D, E2, IFIRST, ILAST,                       &
     &                   RTOL, OFFSET, W, WERR, WORK, IWORK,            &
     &                   PIVMIN, SPDIAM, INFO)
!
!=======================================================================
!                                                                      !
!  Given the initial eigenvalue approximations of T, DLARRJ does       !
!  bisection to refine the eigenvalues of T, W(IFIRST-OFFSET)          !
!  through W(ILAST-OFFSET), to more accuracy. Initial guesses for      !
!  these eigenvalues are input in W, the corresponding estimate of     !
!  the error in these guesses in WERR. During bisection, intervals     !
!  [left, right] are maintained by storing their mid-points and        !
!  semi-widths in the arrays W and WERR respectively.                  !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix.                                   !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           The N diagonal elements of T.                              !
!                                                                      !
!  E2       [in]      REAL array, dimension (N-1)                      !
!           The Squares of the (N-1) subdiagonal elements of T.        !
!                                                                      !
!  IFIRST   [in]      INTEGER                                          !
!           The index of the first eigenvalue to be computed.          !
!                                                                      !
!  ILAST    [in]      INTEGER                                          !
!           The index of the last eigenvalue to be computed.           !
!                                                                      !
!  RTOL     [in]      REAL                                             !
!           Tolerance for the convergence of the bisection intervals.  !
!           An interval [LEFT,RIGHT] has converged if                  !
!           RIGHT-LEFT < RTOL*MAX(|LEFT|,|RIGHT|).                     !
!                                                                      !
!  OFFSET   [in]      INTEGER                                          !
!           Offset for the arrays W and WERR, i.e., the IFIRST-OFFSET  !
!           through ILAST-OFFSET elements of these arrays are to be    !
!           used.                                                      !
!                                                                      !
!  W        [in,out]  REAL array, dimension (N)                        !
!           On input, W(IFIRST-OFFSET) through W(ILAST-OFFSET) are     !
!           estimates of the eigenvalues of L D L^T indexed IFIRST     !
!           through ILAST.                                             !
!           On output, these estimates are refined.                    !
!                                                                      !
!  WERR     [in,out]  REAL array, dimension (N)                        !
!           On input, WERR(IFIRST-OFFSET) through WERR(ILAST-OFFSET)   !
!           are the errors in the estimates of the corresponding       !
!           elements in W.                                             !
!           On output, these errors are refined.                       !
!                                                                      !
!  WORK     [out]     REAL array, dimension (2*N)                      !
!           Workspace.                                                 !
!                                                                      !
!  IWORK    [out]     INTEGER array, dimension (2*N)                   !
!           Workspace.                                                 !
!                                                                      !
!  PIVMIN   [in]      REAL                                             !
!           The minimum pivot in the Sturm sequence for T.             !
!                                                                      !
!  SPDIAM   [in]      REAL                                             !
!           The spectral diameter of T.                                !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           Error flag.                                                !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: IFIRST, ILAST, N, OFFSET
      integer,   intent(out)   :: INFO
      integer,   intent(out)   :: IWORK(*)
!
      real (dp), intent(in)    :: PIVMIN, RTOL, SPDIAM
      real (dp), intent(in)    :: D(*), E2(*)
      real (dp), intent(inout) :: W(*), WERR(*)
      real (dp), intent(out)   :: WORK(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: HALF = 0.5_dp
!
      integer   :: MAXITR
      integer   :: CNT, I, I1, I2, II, ITER, J, K, NEXT, NINT
      integer   :: OLNINT, P, PREV, SAVI1
!
      real (dp) :: DPLUS, FAC, LEFT, MID, RIGHT, S, TMP, WIDTH
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
!
!  Quick return if possible
!
      IF (N.LE.0) THEN
         RETURN
      END IF
!
      MAXITR = INT((LOG(SPDIAM+PIVMIN)-LOG(PIVMIN)) /                   &
     &             LOG(TWO)) + 2
!
!  Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
!  The Sturm Count, Count(WORK(2*I-1)) is arranged to be I-1, while
!  COUNT(WORK(2*I)) is stored in IWORK(2*I).
!  The integer IWORK(2*I-1) for an unconverged interval is set to
!  the index of the next unconverged interval, and is -1 or 0 for a
!  converged interval. Thus, a linked list of unconverged intervals
! is set up.
!
      I1 = IFIRST
      I2 = ILAST
      NINT = 0         ! the number of unconverged intervals
      PREV = 0         ! the last unconverged interval found
!
      DO 75 I = I1, I2
         K = 2*I
         II = I - OFFSET
         LEFT = W(II) - WERR(II)
         MID = W(II)
         RIGHT = W(II) + WERR(II)
         WIDTH = RIGHT - MID
         TMP = MAX(ABS(LEFT), ABS(RIGHT))
!
!  The following test prevents the test of converged intervals
!
         IF (WIDTH.LT.RTOL*TMP) THEN
!
!  This interval has already converged and does not need refinement.
!  (Note that the gaps might change through refining the eigenvalues,
!  however, they can only get bigger.)
!  Remove it from the list.
!
            IWORK(K-1) = -1
!
!  Make sure that I1 always points to the first unconverged interval.
!
            IF ((I.EQ.I1).AND.(I.LT.I2)) I1 = I + 1
            IF ((PREV.GE.I1).AND.(I.LE.I2)) IWORK(2*PREV-1) = I + 1

         ELSE
!
!  Unconverged interval found.
!
            PREV = I
!
!  Make sure that [LEFT,RIGHT] contains the desired eigenvalue
!
!           DO WHILE (CNT(LEFT).GT.I-1)
!
            FAC = ONE
 20         CONTINUE
            CNT = 0
            S = LEFT
            DPLUS = D(1) - S
            IF (DPLUS.LT.ZERO) CNT = CNT + 1
            DO J = 2, N
               DPLUS = D(J) - S - E2(J-1)/DPLUS
               IF (DPLUS.LT.ZERO) CNT = CNT + 1
            END DO
            IF (CNT.GT.I-1) THEN
               LEFT = LEFT - WERR(II)*FAC
               FAC = TWO*FAC
               GOTO 20
            END IF
!
!           DO WHILE (CNT(RIGHT).LT.I)
!
            FAC = ONE
 50         CONTINUE
            CNT = 0
            S = RIGHT
            DPLUS = D(1) - S
            IF (DPLUS.LT.ZERO) CNT = CNT + 1
            DO J = 2, N
               DPLUS = D(J) - S - E2(J-1)/DPLUS
               IF (DPLUS.LT.ZERO) CNT = CNT + 1
            END DO
            IF (CNT.LT.I) THEN
               RIGHT = RIGHT + WERR(II)*FAC
               FAC = TWO*FAC
               GOTO 50
            END IF
            NINT = NINT + 1
            IWORK(K-1) = I + 1
            IWORK(K) = CNT
         END IF
         WORK(K-1) = LEFT
         WORK(K) = RIGHT
 75   CONTINUE

      SAVI1 = I1
!
!     DO WHILE (NINT.GT.0)
!     i.e. there are still unconverged intervals and
!     WHILE (ITER.LT.MAXITR)
!
      ITER = 0
 80   CONTINUE
      PREV = I1 - 1
      I = I1
      OLNINT = NINT

      DO 100 P = 1, OLNINT
         K = 2*I
         II = I - OFFSET
         NEXT = IWORK(K-1)
         LEFT = WORK(K-1)
         RIGHT = WORK(K)
         MID = HALF*(LEFT + RIGHT)
!
!  Semiwidth of interval.
!
         WIDTH = RIGHT - MID
         TMP = MAX(ABS(LEFT), ABS(RIGHT))

         IF ((WIDTH.LT.RTOL*TMP) .OR. (ITER.EQ.MAXITR)) THEN
!
!  Reduce number of unconverged intervals.
!
            NINT = NINT - 1
!
!  Mark interval as converged.
!
            IWORK(K-1) = 0
            IF (I1.EQ.I) THEN
               I1 = NEXT
            ELSE
!
!  PREV holds the last unconverged interval previously examined.
!
               IF (PREV.GE.I1) IWORK(2*PREV-1) = NEXT
            END IF
            I = NEXT
            GOTO 100
         END IF
         PREV = I
!
!  Perform one bisection step.
!
         CNT = 0
         S = MID
         DPLUS = D(1) - S
         IF (DPLUS.LT.ZERO) CNT = CNT + 1
         DO J = 2, N
            DPLUS = D(J) - S - E2(J-1)/DPLUS
            IF (DPLUS.LT.ZERO) CNT = CNT + 1
         END DO
         IF (CNT.LE.I-1) THEN
            WORK(K-1) = MID
         ELSE
            WORK(K) = MID
         END IF
         I = NEXT

 100  CONTINUE
      ITER = ITER + 1
!
!  Do another loop if there are still unconverged intervals
!  However, in the last iteration, all intervals are accepted
!  since this is the best we can do.
!
      IF ((NINT.GT.0).AND.(ITER.LE.MAXITR)) GOTO 80
!
!  At this point, all the intervals have converged.
!
      DO I = SAVI1, ILAST
         K = 2*I
         II = I - OFFSET
!
!  All intervals marked by '0' have been refined.
!
         IF (IWORK(K-1).EQ.0) THEN
            W(II) = HALF*(WORK(K-1) + WORK(K))
            WERR(II) = WORK(K) - W(II)
         END IF
      END DO
!
      RETURN
      END SUBROUTINE DLARRJ
!
      SUBROUTINE DLARRK (N, IW, GL, GU,                                 &
     &                   D, E2, PIVMIN, RELTOL, W, WERR, INFO)
!
!=======================================================================
!                                                                      !
!  DLARRK computes one eigenvalue of a symmetric tridiagonal matrix T  !
!  to suitable accuracy. This is an auxiliary code to be called from   !
!  DSTEMR.                                                             !
!                                                                      !
!  To avoid overflow, the matrix must be scaled so that its largest    !
!  element is no greater than overflow**(1/2) * underflow**(1/4) in    !
!  absolute value, and for greatest accuracy, it should not be much    !
!  smaller than that.                                                  !
!                                                                      !
!  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal       !
!  Matrix", Report CS41, Computer Science Dept., Stanford              !
!  University, July 21, 1966.                                          !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the tridiagonal matrix T,  N >= 0.            !
!                                                                      !
!  IW       [in]      INTEGER                                          !
!           The index of the eigenvalues to be returned.               !
!                                                                      !
!  GL       [in]      REAL                                             !
!                                                                      !
!  GU       [in]      REAL                                             !
!           An upper and a lower bound on the eigenvalue.              !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           The N diagonal elements of the tridiagonal matrix T.       !
!                                                                      !
!  E2       [in]      REAL array, dimension (N-1)                      !
!           The (N-1) squared off-diagonal elements of the tridiagonal !
!           matrix T.                                                  !
!                                                                      !
!  PIVMIN   [in]      REAL                                             !
!           The minimum pivot allowed in the Sturm sequence for T.     !
!                                                                      !
!  RELTOL   [in]      REAL                                             !
!           The minimum relative width of an interval.  When an        !
!           interval is narrower than RELTOL times the larger (in      !
!           magnitude) endpoint, then it is considered to be           !
!           sufficiently small, i.e., converged.  Note: this should    !
!           always be at least radix*machine epsilon.                  !
!                                                                      !
!  W        [out]     REAL                                             !
!                                                                      !
!  WERR     [out]     REAL                                             !
!           The error bound on the corresponding eigenvalue            !
!           approximation in W.                                        !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           = 0:       Eigenvalue converged                            !
!           = -1:      Eigenvalue did NOT converge                     !
!                                                                      !
!  Internal Parameters:                                                !
!                                                                      !
!   FUDGE   DOUBLE PRECISION, default = 2                              !
!           A "fudge factor" to widen the Gershgorin intervals.        !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)  :: IW, N
      integer,   intent(out) :: INFO
!
      real (dp), intent(in)  :: PIVMIN, RELTOL, GL, GU
      real (dp), intent(in)  :: D(*), E2(*)
      real (dp), intent(out) :: W, WERR
!
!  Local variable declarations.
!
      real (dp), parameter :: HALF = 0.5_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: FUDGE = TWO
      real (dp), parameter :: ZERO = 0.0_dp
!
      integer   :: I, IT, ITMAX, NEGCNT
!
      real (dp) :: ATOLI, EPS, LEFT, MID, RIGHT, RTOLI, TMP1
      real (dp) :: TMP2, TNORM
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Quick return if possible.
!
      IF (N.LE.0) THEN
         INFO = 0
         RETURN
      END IF
!
!  Get machine constants.
!
      EPS = DLAMCH('P')

      TNORM = MAX(ABS(GL), ABS(GU))
      RTOLI = RELTOL
      ATOLI = FUDGE*TWO*PIVMIN

      ITMAX = INT((LOG(TNORM+PIVMIN)-LOG(PIVMIN)) / LOG(TWO)) + 2

      INFO = -1

      LEFT = GL - FUDGE*TNORM*EPS*N - FUDGE*TWO*PIVMIN
      RIGHT = GU + FUDGE*TNORM*EPS*N + FUDGE*TWO*PIVMIN
      IT = 0

 10   CONTINUE
!
!  Check if interval converged or maximum number of iterations reached.
!
      TMP1 = ABS(RIGHT - LEFT)
      TMP2 = MAX(ABS(RIGHT), ABS(LEFT))
      IF (TMP1.LT.MAX(ATOLI, PIVMIN, RTOLI*TMP2)) THEN
         INFO = 0
         GOTO 30
      END IF
      IF (IT.GT.ITMAX) GOTO 30
!
!  Count number of negative pivots for mid-point.
!
      IT = IT + 1
      MID = HALF * (LEFT + RIGHT)
      NEGCNT = 0
      TMP1 = D(1) - MID
      IF (ABS(TMP1).LT.PIVMIN) TMP1 = -PIVMIN
      IF (TMP1.LE.ZERO) NEGCNT = NEGCNT + 1
!
      DO I = 2, N
         TMP1 = D(I) - E2(I-1) / TMP1 - MID
         IF (ABS(TMP1).LT.PIVMIN) TMP1 = -PIVMIN
         IF (TMP1.LE.ZERO) NEGCNT = NEGCNT + 1
      END DO

      IF (NEGCNT.GE.IW) THEN
         RIGHT = MID
      ELSE
         LEFT = MID
      END IF
      GOTO 10

 30   CONTINUE
!
!  Converged or maximum number of iterations reached.
!
      W = HALF * (LEFT + RIGHT)
      WERR = HALF * ABS(RIGHT - LEFT)
!
      RETURN
      END SUBROUTINE DLARRK
!
      SUBROUTINE DLARRR (N, D, E, INFO)
!
!=======================================================================
!                                                                      !
!  Perform tests to decide whether the symmetric tridiagonal matrix    !
!  T warrants expensive computations which guarantee high relative     !
!  accuracy in the eigenvalues.                                        !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix, N > 0.                            !
!                                                                      !
!  D        [in]      REAL array, dimension (N)                        !
!           The N diagonal elements of the tridiagonal matrix T.       !
!                                                                      !
!  E        [in,out]  REAL array, dimension (N)                        !
!           On entry, the first (N-1) entries contain the subdiagonal  !
!           elements of the tridiagonal matrix T; E(N) is set to ZERO. !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           INFO = 0 (default) : The matrix warrants computations      !
!                                preserving relative accuracy.         !
!           INFO = 1           : The matrix warrants computations      !
!                                guaranteeing only absolute accuracy.  !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: N
      integer,   intent(out)   :: INFO
!
      real (dp), intent(in)    :: D(*)
      real (dp), intent(inout) :: E(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: RELCOND = 0.999_dp
!
      logical   :: YESREL
!
      integer   :: I
!
      real (dp) :: EPS, SAFMIN, SMLNUM, RMIN, TMP, TMP2
      real (dp) :: OFFDIG, OFFDIG2
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Quick return if possible.
!
      IF (N.LE.0) THEN
         INFO = 0
         RETURN
      END IF
!
!  As a default, do NOT go for relative-accuracy preserving
!  computations.
!
      INFO = 1
!
      SAFMIN = DLAMCH('Safe minimum')
      EPS = DLAMCH('Precision')
      SMLNUM = SAFMIN / EPS
      RMIN = SQRT(SMLNUM)
!
!  Tests for relative accuracy
!
!  Test for Scaled Diagonal Dominance (SDD):
!  Scale the diagonal entries to one and check whether the sum of the
!  off-diagonals is less than one
!
!  The SDD relative error bounds have a 1/(1- 2*x) factor in them,
!  x = MAX(OFFDIG + OFFDIG2), so when x is close to 1/2, no relative
!  accuracy is promised.  In the notation of the code fragment below,
!  1/(1 - (OFFDIG + OFFDIG2)) is the condition number.
!  We do not think it is worth going into "SDD mode" unless the relative
!  condition number is reasonable, not 1/macheps.
!  The threshold should be compatible with other thresholds used in the
!  code. We set  OFFDIG + OFFDIG2 <= .999 =: RELCOND, it corresponds
!  to losing at most 3 decimal digits:
!        1 / (1 - (OFFDIG + OFFDIG2)) <= 1000
!  instead of the current OFFDIG + OFFDIG2 < 1
!
      YESREL = .TRUE.
      OFFDIG = ZERO
      TMP = SQRT(ABS(D(1)))
      IF (TMP.LT.RMIN) YESREL = .FALSE.
      IF (.NOT.YESREL) GOTO 11
!
      DO I = 2, N
         TMP2 = SQRT(ABS(D(I)))
         IF (TMP2.LT.RMIN) YESREL = .FALSE.
         IF (.NOT.YESREL) GOTO 11
         OFFDIG2 = ABS(E(I-1))/(TMP*TMP2)
         IF (OFFDIG+OFFDIG2.GE.RELCOND) YESREL = .FALSE.
         IF (.NOT.YESREL) GOTO 11
         TMP = TMP2
         OFFDIG = OFFDIG2
      END DO
 11   CONTINUE

      IF (YESREL) THEN
         INFO = 0
         RETURN
      ELSE
      END IF
!
!  MORE TO BE IMPLEMENTED
!

!
!  Test if the lower bidiagonal matrix L from T = L D L^T
!  (zero shift factor) is well conditioned
!

!
!  Test if the upper bidiagonal matrix U from T = U D U^T
!  (zero shift factor) is well conditioned.
!  In this case, the matrix needs to be flipped and, at the end
!  of the eigenvector computation, the flip needs to be applied
!  to the computed eigenvectors (and the support)
!

!
      RETURN
      END SUBROUTINE DLARRR
!
      SUBROUTINE DLARRV (N, VL, VU, D, L, PIVMIN,                       &
     &                   ISPLIT, M, DOL, DOU, MINRGP,                   &
     &                   RTOL1, RTOL2, W, WERR, WGAP,                   &
     &                   IBLOCK, INDEXW, GERS, Z, LDZ, ISUPPZ,          &
     &                   WORK, IWORK, INFO)
!
!=======================================================================
!                                                                      !
!  DLARRV computes the eigenvectors of the tridiagonal matrix          !
!  T = L D L^T given L, D and APPROXIMATIONS to the eigenvalues of     !
!  L D L^T.  The input eigenvalues should have been computed by        !
!  DLARRE.                                                             !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the matrix, N >= 0.                           !
!                                                                      !
!  VL       [in]      REAL                                             !
!           Lower bound of the interval that contains the desired      !
!           eigenvalues. VL < VU. Needed to compute gaps on the left   !
!           or right end of the extremal eigenvalues in the desired    !
!           RANGE.                                                     !
!                                                                      !
!  VU       [in]      REAL                                             !
!           Upper bound of the interval that contains the desired      !
!           eigenvalues. VL < VU.                                      !
!           Note: VU is currently not used by this implementation of   !
!           DLARRV, VU is passed to DLARRV because it could be used    !
!           to compute gaps on the right end of the extremal           !
!           eigenvalues. However, with not much initial accuracy in    !
!           LAMBDA and VU, the formula can lead to an overestimation   !
!           of the right gap and thus to inadequately early RQI        !
!           'convergence'. This is currently prevented this by forcing !
!           a small right gap. And so it turns out that VU is          !
!           currently not used by this implementation of DLARRV.       !
!                                                                      !
!  D        [in,out]  REAL  array, dimension (N)                       !
!           On entry, the N diagonal elements of the diagonal matrix D !
!           On exit, D may be overwritten.                             !
!                                                                      !
!  L        [in,out]  REAL array, dimension (N)                        !
!           On entry, the (N-1) subdiagonal elements of the unit       !
!           bidiagonal matrix L are in elements 1 to N-1 of L          !
!           (if the matrix is not split.) At the end of each block     !
!           is stored the corresponding shift as given by DLARRE.      !
!           On exit, L is overwritten.                                 !
!                                                                      !
!  PIVMIN   [in]      REAL                                             !
!           The minimum pivot allowed in the Sturm sequence.           !
!                                                                      !
!  ISPLIT   [in]      INTEGER array, dimension (N)                     !
!           The splitting points, at which T breaks up into blocks.    !
!           The first block consists of rows/columns 1 to              !
!           ISPLIT(1), the second of rows/columns ISPLIT(1)+1          !
!           through ISPLIT(2), etc.                                    !
!                                                                      !
!  M        [in]      INTEGER                                          !
!           The total number of input eigenvalues,  0 <= M <= N.       !
!                                                                      !
!  DOL      [in]      INTEGER                                          !
!                                                                      !
!  DOU      [in]      INTEGER                                          !
!           If the user wants to compute only selected eigenvectors    !
!           from all the eigenvalues supplied, he can specify an index !
!           range DOL:DOU.                                             !
!           Or else the setting DOL=1, DOU=M should be applied.        !
!           Note that DOL and DOU refer to the order in which the      !
!           eigenvalues are stored in W.                               !
!           If the user wants to compute only selected eigenpairs,     !
!           then the columns DOL-1 to DOU+1 of the eigenvector space   !
!           Z contain the computed eigenvectors. All other columns of  !
!           Z are set to zero.                                         !
!                                                                      !
!  MINRGP   [in]      REAL                                             !
!                                                                      !
!  RTOL1    [in,out]  REAL                                             !
!                                                                      !
!  RTOL2    [in,out]  REAL                                             !
!           Parameters for bisection.                                  !
!           An interval [LEFT,RIGHT] has converged if                  !
!           RIGHT-LEFT < MAX(RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|))     !
!                                                                      !
!  W        [in,out]  REAL array, dimension (N)                        !
!           The first M elements of W contain the APPROXIMATE          !
!           eigenvalues for which eigenvectors are to be computed.     !
!           The eigenvalues should be grouped by split-off block and   !
!           ordered from smallest to largest within the block (The     !
!           output array W from DLARRE is expected here). Furthermore, !
!           they are with respect to the shift of the corresponding    !
!           root representation for their block. On exit, W holds the  !
!           eigenvalues of the UNshifted matrix.                       !
!                                                                      !
!  WERR     [in,out]  REAL array, dimension (N)                        !
!           The first M elements contain the semiwidth of the          !
!           uncertainty interval of the corresponding eigenvalue in W. !
!                                                                      !
!  WGAP     [in,out]  REAL array, dimension (N)                        !
!           The separation from the right neighbor eigenvalue in W.    !
!                                                                      !
!  IBLOCK   [in]      INTEGER array, dimension (N)                     !
!           The indices of the blocks (submatrices) associated with    !
!           the corresponding eigenvalues in W.                        !
!           IBLOCK(i)=1 if eigenvalue W(i) belongs to the first block  !
!           from the top, =2 if W(i) belongs to the second block, etc. !
!                                                                      !
!  INDEXW   [in]      INTEGER array, dimension (N)                     !
!           The indices of the eigenvalues within each block           !
!           (submatrix); for example, INDEXW(i)= 10 and IBLOCK(i)=2    !
!           imply that the i-th eigenvalue W(i) is the 10-th           !
!           eigenvalue in the second block.                            !
!                                                                      !
!  GERS     [in]      REAL array, dimension (2*N)                      !
!           The N Gerschgorin intervals (the i-th Gerschgorin interval !
!           is (GERS(2*i-1), GERS(2*i)). The Gerschgorin intervals     !
!           should be computed from the original UNshifted matrix.     !
!                                                                      !
!  Z        [out] Z   REAL array, dimension (LDZ, MAX(1,M))            !
!           If INFO = 0, the first M columns of Z contain the          !
!           orthonormal eigenvectors of the matrix T                   !
!           corresponding to the input eigenvalues, with the i-th      !
!           column of Z holding the eigenvector associated with W(i).  !
!           Note: the user must ensure that at least max(1,M) columns  !
!           are supplied in the array Z.                               !
!                                                                      !
!  LDZ      [in]      INTEGER                                          !
!           The leading dimension of the array Z.  LDZ >= 1, and if    !
!           JOBZ = 'V', LDZ >= MAX(1,N).                               !
!                                                                      !
!  ISUPPZ   [out]     INTEGER array, dimension (2*max(1,M))            !
!           The support of the eigenvectors in Z, i.e., the indices    !
!           indicating the nonzero elements in Z. The I-th eigenvector !
!           is nonzero only in elements ISUPPZ(2*I-1) through          !
!           ISUPPZ(2*I).                                               !
!                                                                      !
!  WORK     [out]     REAL array, dimension (12*N)                     !
!                                                                      !
!  IWORK    [out]     INTEGER array, dimension (7*N)                   !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!           = 0:  successful exit                                      !
!                                                                      !
!           > 0:  A problem occurred in DLARRV.                        !
!           < 0:  One of the called subroutines signaled an internal   !
!                 problem. Needs inspection of the corresponding       !
!                 parameter IINFO for further information.             !
!                                                                      !
!           =-1:  Problem in DLARRB when refining a child eigenvalues. !
!           =-2:  Problem in DLARRF when computing the RRR of a child. !
!                 When a child is inside a tight cluster, it can be    !
!                 difficult to find an RRR. A partial remedy from the  !
!                 user point of view is to make the parameter MINRGP   !
!                 smaller and recompile. However, as the orthogonality !
!                 of the computed vectors is proportional to 1/MINRGP, !
!                 the user should be aware that he might be trading in !
!                 precision when he decreases MINRGP.                  !
!           =-3:  Problem in DLARRB when refining a single eigenvalue  !
!                 after the Rayleigh correction was rejected.          !
!           = 5:  The Rayleigh Quotient Iteration failed to converge   !
!                 to full accuracy in MAXITR steps.                    !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Contributors:                                                       !
!                                                                      !
!  Beresford Parlett, University of California, Berkeley, USA          !
!  Jim Demmel, University of California, Berkeley, USA                 !
!  Inderjit Dhillon, University of Texas, Austin, USA                  !
!  Osni Marques, LBNL/NERSC, USA                                       !
!  Christof Voemel, University of California, Berkeley, USA            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: DOL, DOU, LDZ, M, N
      integer,   intent(in)    :: IBLOCK(*), INDEXW(*), ISPLIT(*)
      integer,   intent(out)   :: INFO
      integer,   intent(out)   :: ISUPPZ(*), IWORK(*)
!
      real (dp), intent(in)    :: MINRGP, PIVMIN, VL, VU
      real (dp), intent(in)    :: GERS(*)
      real (dp), intent(inout) :: RTOL1, RTOL2
      real (dp), intent(inout) :: D(*), L(*), W(*), WERR(*), WGAP(*)
      real (dp), intent(out)   :: WORK(*), Z(LDZ,*)
!
!  Local variable declarations.
!
      integer, parameter :: MAXITR = 10
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: THREE = 3.0_dp
      real (dp), parameter :: FOUR = 4.0_dp
      real (dp), parameter :: HALF = 0.5_dp
!
      logical   :: ESKIP, NEEDBS, STP2II, TRYRQC, USEDBS, USEDRQ
!
      integer   :: DONE, I, IBEGIN, IDONE, IEND, II, IINDC1
      integer   :: IINDC2, IINDR, IINDWK, IINFO, IM, IN, INDEIG
      integer   :: INDLD, INDLLD, INDWRK, ISUPMN, ISUPMX, ITER
      integer   :: ITMP1, J, JBLK, K, MINIWSIZE, MINWSIZE, NCLUS
      integer   :: NDEPTH, NEGCNT, NEWCLS, NEWFST, NEWFTT, NEWLST
      integer   :: NEWSIZ, OFFSET, OLDCLS, OLDFST, OLDIEN, OLDLST
      integer   :: OLDNCL, P, PARITY, Q, WBEGIN, WEND, WINDEX
      integer   :: WINDMN, WINDPL, ZFROM, ZTO, ZUSEDL, ZUSEDU
      integer   :: ZUSEDW
!
      real (dp) :: BSTRES, BSTW, EPS, FUDGE, GAP, GAPTOL, GL, GU
      real (dp) :: LAMBDA, LEFT, LGAP, MINGMA, NRMINV, RESID
      real (dp) :: RGAP, RIGHT, RQCORR, RQTOL, SAVGAP, SGNDEF
      real (dp) :: SIGMA, SPDIAM, SSIGMA, TAU, TMP, TOL, ZTZ
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      INFO = 0
!
!  Quick return if possible
!
      IF ((N.LE.0).OR.(M.LE.0)) THEN
         RETURN
      END IF
!
!  The first N entries of WORK are reserved for the eigenvalues
!
      INDLD = N+1
      INDLLD= 2*N+1
      INDWRK= 3*N+1
      MINWSIZE = 12 * N
!
      DO I = 1, MINWSIZE
        WORK(I) = ZERO
      END DO
!
!  IWORK(IINDR+1:IINDR+N) hold the twist indices R for the
!  factorization used to compute the FP vector
!
      IINDR = 0
!
!  IWORK(IINDC1+1:IINC2+N) are used to store the clusters of the current
!  layer and the one above.
!
      IINDC1 = N
      IINDC2 = 2*N
      IINDWK = 3*N + 1

      MINIWSIZE = 7 * N
      DO I = 1, MINIWSIZE
        IWORK(I) = 0
      END DO

      ZUSEDL = 1
      IF (DOL.GT.1) THEN
        ZUSEDL = DOL-1           ! set lower bound for use of Z
      END IF
      ZUSEDU = M
      IF (DOU.LT.M) THEN
        ZUSEDU = DOU+1           ! set lower bound for use of Z
      END IF
!
!  The width of the part of Z that is used
!
      ZUSEDW = ZUSEDU - ZUSEDL + 1

      CALL DLASET ('Full', N, ZUSEDW, ZERO, ZERO, Z(1,ZUSEDL), LDZ)

      EPS = DLAMCH('Precision')
      RQTOL = TWO * EPS
!
!  Set expert flags for standard code.
!
      TRYRQC = .TRUE.

      IF ((DOL.EQ.1).AND.(DOU.EQ.M)) THEN
      ELSE
!
!  Only selected eigenpairs are computed. Since the other evalues
!  are not refined by RQ iteration, bisection has to compute to full
!   accuracy.
!
         RTOL1 = FOUR * EPS
         RTOL2 = FOUR * EPS
      END IF
!
!  The entries WBEGIN:WEND in W, WERR, WGAP correspond to the
!  desired eigenvalues. The support of the nonzero eigenvector
!  entries is contained in the interval IBEGIN:IEND.
!  Remark that if k eigenpairs are desired, then the eigenvectors
!  are stored in k contiguous columns of Z.
!
      DONE = 0          ! number of eigenvectors already computed
      IBEGIN = 1
      WBEGIN = 1

      DO 170 JBLK = 1, IBLOCK(M)
         IEND = ISPLIT(JBLK)
         SIGMA = L(IEND)
!
!  Find the eigenvectors of the submatrix indexed IBEGIN through IEND.
!
         WEND = WBEGIN - 1
 15      CONTINUE
         IF (WEND.LT.M) THEN
            IF (IBLOCK(WEND+1).EQ.JBLK) THEN
               WEND = WEND + 1
               GOTO 15
            END IF
         END IF
         IF (WEND.LT.WBEGIN) THEN
            IBEGIN = IEND + 1
            GOTO 170
         ELSE IF ((WEND.LT.DOL).OR.(WBEGIN.GT.DOU)) THEN
            IBEGIN = IEND + 1
            WBEGIN = WEND + 1
            GOTO 170
         END IF
!
!  Find local spectral diameter of the block.
!
         GL = GERS(2*IBEGIN-1)
         GU = GERS(2*IBEGIN)
         DO I = IBEGIN+1 , IEND
            GL = MIN(GERS(2*I-1), GL)
            GU = MAX(GERS(2*I), GU)
         END DO
         SPDIAM = GU - GL
!
!  OLDIEN is the last index of the previous block
!
         OLDIEN = IBEGIN - 1
!
!  Calculate the size of the current block
!
         IN = IEND - IBEGIN + 1
!
!  The number of eigenvalues in the current block
!
         IM = WEND - WBEGIN + 1
!
!  This is for a 1x1 block
!
         IF (IBEGIN.EQ.IEND) THEN
            DONE = DONE+1
            Z(IBEGIN, WBEGIN) = ONE
            ISUPPZ(2*WBEGIN-1) = IBEGIN
            ISUPPZ(2*WBEGIN) = IBEGIN
            W(WBEGIN) = W(WBEGIN) + SIGMA
            WORK(WBEGIN) = W(WBEGIN)
            IBEGIN = IEND + 1
            WBEGIN = WBEGIN + 1
            GOTO 170
         END IF
!
!  The desired (shifted) eigenvalues are stored in W(WBEGIN:WEND)
!  Note that these can be approximations, in this case, the corresp.
!  entries of WERR give the size of the uncertainty interval.
!  The eigenvalue approximations will be refined when necessary as
!  high relative accuracy is required for the computation of the
!  corresponding eigenvectors.
!
         CALL DCOPY (IM, W(WBEGIN), 1, WORK(WBEGIN), 1)
!
!   We store in W the eigenvalue approximations w.r.t. the original
!   matrix T.
!
         DO I=1,IM
            W(WBEGIN+I-1) = W(WBEGIN+I-1)+SIGMA
         END DO
!
!  NDEPTH is the current depth of the representation tree
!
         NDEPTH = 0
!
!  PARITY is either 1 or 0
!
         PARITY = 1
!
!  NCLUS is the number of clusters for the next level of the
!  representation tree, we start with NCLUS = 1 for the root
!
         NCLUS = 1
         IWORK(IINDC1+1) = 1
         IWORK(IINDC1+2) = IM
!
!  IDONE is the number of eigenvectors already computed in the current
!  block
!
         IDONE = 0
!
!  Generate the representation tree for the current block and
!  compute the eigenvectors.
!
!        loop : WHILE (IDONE.LT.IM)
!
   40    CONTINUE
         IF (IDONE.LT.IM) THEN
!
!  This is a crude protection against infinitely deep trees.
!
            IF (NDEPTH.GT.M) THEN
               INFO = -2
               RETURN
            END IF
!
!  Breadth first processing of the current level of the representation
!  tree: OLDNCL = number of clusters on current level
!
            OLDNCL = NCLUS
!
!  Reset NCLUS to count the number of child clusters.
!
            NCLUS = 0
!
            PARITY = 1 - PARITY
            IF (PARITY.EQ.0) THEN
               OLDCLS = IINDC1
               NEWCLS = IINDC2
            ELSE
               OLDCLS = IINDC2
               NEWCLS = IINDC1
            END IF
!
!  Process the clusters on the current level.
!
            DO 150 I = 1, OLDNCL
               J = OLDCLS + 2*I
!
!  OLDFST, OLDLST = first, last index of current cluster.
!                   cluster indices start with 1 and are relative
!                   to WBEGIN when accessing W, WGAP, WERR, Z
!
               OLDFST = IWORK(J-1)
               OLDLST = IWORK(J)
               IF (NDEPTH.GT.0) THEN
!
!  Retrieve relatively robust representation (RRR) of cluster
!  that has been computed at the previous level.
!  The RRR is stored in Z and overwritten once the eigenvectors
!  have been computed or when the cluster is refined
!
                  IF ((DOL.EQ.1).AND.(DOU.EQ.M)) THEN
!
!  Get representation from location of the leftmost evalue
!  of the cluster
!
                     J = WBEGIN + OLDFST - 1
                  ELSE
                     IF (WBEGIN+OLDFST-1.LT.DOL) THEN
!
!  Get representation from the left end of Z array.
!
                        J = DOL - 1
                     ELSE IF (WBEGIN+OLDFST-1.GT.DOU) THEN
!
!  Get representation from the right end of Z array
!
                        J = DOU
                     ELSE
                        J = WBEGIN + OLDFST - 1
                     END IF
                  END IF
                  CALL DCOPY (IN, Z(IBEGIN,J), 1, D(IBEGIN), 1)
                  CALL DCOPY (IN-1, Z(IBEGIN, J+1), 1, L(IBEGIN), 1)
                  SIGMA = Z(IEND, J+1)
!
!  Set the corresponding entries in Z to zero
!
                  CALL DLASET ('Full', IN, 2, ZERO, ZERO,               &
     &                         Z(IBEGIN, J), LDZ)
               END IF
!
!   Compute DL and DLL of current RRR.
!
               DO J = IBEGIN, IEND-1
                  TMP = D(J)*L(J)
                  WORK(INDLD-1+J) = TMP
                  WORK(INDLLD-1+J) = TMP*L(J)
               END DO

               IF (NDEPTH.GT.0) THEN
!
!  P and Q are index of the first and last eigenvalue to compute
!  within the current block
!
                  P = INDEXW(WBEGIN-1+OLDFST)
                  Q = INDEXW(WBEGIN-1+OLDLST)
!
!  Offset for the arrays WORK, WGAP and WERR, i.e., the P-OFFSET
!  through the Q-OFFSET elements of these arrays are to be used.
!  OFFSET = P-OLDFST
!
                  OFFSET = INDEXW(WBEGIN) - 1
!
!  Perform limited bisection (if necessary) to get approximate
!  eigenvalues to the precision needed.
!
                  CALL DLARRB (IN, D(IBEGIN),                           &
     &                         WORK(INDLLD+IBEGIN-1),                   &
     &                         P, Q, RTOL1, RTOL2, OFFSET,              &
     &                         WORK(WBEGIN),WGAP(WBEGIN),WERR(WBEGIN),  &
     &                         WORK(INDWRK), IWORK(IINDWK),             &
     &                         PIVMIN, SPDIAM, IN, IINFO)
                  IF (IINFO.NE.0) THEN
                     INFO = -1
                     RETURN
                  END IF
!
!  We also recompute the extremal gaps. W holds all eigenvalues
!  of the unshifted matrix and must be used for computation
!  of WGAP, the entries of WORK might stem from RRRs with
!  different shifts. The gaps from WBEGIN-1+OLDFST to
!  WBEGIN-1+OLDLST are correctly computed in DLARRB.
!  However, we only allow the gaps to become greater since
!  this is what should happen when we decrease WERR
!
                  IF (OLDFST.GT.1) THEN
                     WGAP(WBEGIN+OLDFST-2) =                            &
     &                   MAX(WGAP(WBEGIN+OLDFST-2),                     &
     &                       W(WBEGIN+OLDFST-1)-WERR(WBEGIN+OLDFST-1)-  &
     &                       W(WBEGIN+OLDFST-2)-WERR(WBEGIN+OLDFST-2))
                  END IF
                  IF (WBEGIN + OLDLST -1 .LT. WEND) THEN
                     WGAP(WBEGIN+OLDLST-1) =                            &
     &                   MAX(WGAP(WBEGIN+OLDLST-1),                     &
     &                       W(WBEGIN+OLDLST)-WERR(WBEGIN+OLDLST)-      &
     &                       W(WBEGIN+OLDLST-1)-WERR(WBEGIN+OLDLST-1))
                  END IF
!
!  Each time the eigenvalues in WORK get refined, we store
!  the newly found approximation with all shifts applied in W
!
                  DO J=OLDFST,OLDLST
                     W(WBEGIN+J-1) = WORK(WBEGIN+J-1)+SIGMA
                  END DO
               END IF
!
!  Process the current node.
!
               NEWFST = OLDFST
               DO 140 J = OLDFST, OLDLST
                  IF (J.EQ.OLDLST) THEN
!
!  We are at the right end of the cluster, this is also the
!  boundary of the child cluster.
!
                     NEWLST = J
                  ELSE IF (WGAP(WBEGIN + J -1).GE.                      &
     &                     MINRGP* ABS(WORK(WBEGIN + J -1))) THEN
!
!  The right relative gap is big enough, the child cluster
!  (NEWFST,..,NEWLST) is well separated from the following
!
                     NEWLST = J
                   ELSE
!
!  Inside a child cluster, the relative gap is not big enough.
!
                     GOTO 140
                  END IF
!
!  Compute size of child cluster found
!
                  NEWSIZ = NEWLST - NEWFST + 1
!
!  NEWFTT is the place in Z where the new RRR or the computed
!  eigenvector is to be stored.
!
                  IF ((DOL.EQ.1).AND.(DOU.EQ.M)) THEN
!
!  Store representation at location of the leftmost evalue
!  of the cluster.
!
                     NEWFTT = WBEGIN + NEWFST - 1
                  ELSE
                     IF (WBEGIN+NEWFST-1.LT.DOL) THEN
!
!  Store representation at the left end of Z array.
!
                        NEWFTT = DOL - 1
                     ELSE IF (WBEGIN+NEWFST-1.GT.DOU) THEN
!
!  Store representation at the right end of Z array.
!
                        NEWFTT = DOU
                     ELSE
                        NEWFTT = WBEGIN + NEWFST - 1
                     END IF
                  END IF

                  IF (NEWSIZ.GT.1) THEN
!
!  Current child is not a singleton but a cluster.
!  Compute and store new representation of child.
!
!  Compute left and right cluster gap.
!
!  LGAP and RGAP are not computed from WORK because
!  the eigenvalue approximations may stem from RRRs
!  different shifts. However, W hold all eigenvalues
!  of the unshifted matrix. Still, the entries in WGAP
!  have to be computed from WORK since the entries
!  in W might be of the same order so that gaps are not
!  exhibited correctly for very close eigenvalues.
!
                     IF (NEWFST.EQ.1) THEN
                        LGAP = MAX(ZERO,                                &
     &                             W(WBEGIN)-WERR(WBEGIN) - VL)
                     ELSE
                        LGAP = WGAP(WBEGIN+NEWFST-2)
                     END IF
                     RGAP = WGAP(WBEGIN+NEWLST-1)
!
!  Compute left- and rightmost eigenvalue of child
!  to high precision in order to shift as close
!  as possible and obtain as large relative gaps
!  as possible.
!
                     DO K =1,2
                        IF (K.EQ.1) THEN
                           P = INDEXW(WBEGIN-1+NEWFST)
                        ELSE
                           P = INDEXW(WBEGIN-1+NEWLST)
                        END IF
                        OFFSET = INDEXW(WBEGIN) - 1
                        CALL DLARRB (IN, D(IBEGIN),                     &
     &                               WORK(INDLLD+IBEGIN-1), P, P,       &
     &                               RQTOL, RQTOL, OFFSET,              &
     &                               WORK(WBEGIN), WGAP(WBEGIN),        &
     &                               WERR(WBEGIN), WORK(INDWRK),        &
     &                               IWORK(IINDWK), PIVMIN, SPDIAM,     &
     &                               IN, IINFO)
                     END DO
!
                     IF ((WBEGIN+NEWLST-1.LT.DOL).OR.                   &
     &                   (WBEGIN+NEWFST-1.GT.DOU)) THEN
!
!  If the cluster contains no desired eigenvalues skip the computation
!  of that branch of the representation tree.
!
!  We could skip before the refinement of the extremal eigenvalues of
!  the child, but then the representation tree could be different from
!  the one when nothing is skipped. For this reason we skip at this
!  place.
!
                        IDONE = IDONE + NEWLST - NEWFST + 1
                        GOTO 139
                     END IF
!
!  Compute RRR of child cluster.
!  Note that the new RRR is stored in Z DLARRF needs LWORK = 2*N
!
                     CALL DLARRF (IN, D(IBEGIN), L(IBEGIN),             &
     &                            WORK(INDLD+IBEGIN-1),                 &
     &                            NEWFST, NEWLST, WORK(WBEGIN),         &
     &                            WGAP(WBEGIN), WERR(WBEGIN),           &
     &                            SPDIAM, LGAP, RGAP, PIVMIN, TAU,      &
     &                            Z(IBEGIN,NEWFTT), Z(IBEGIN,NEWFTT+1), &
     &                            WORK(INDWRK), IINFO)
                     IF (IINFO.EQ.0) THEN
!
!  A new RRR for the cluster was found by DLARRF
!  update shift and store it.
!
                        SSIGMA = SIGMA + TAU
                        Z(IEND, NEWFTT+1) = SSIGMA
!
!  WORK() are the midpoints and WERR() the semi-width
!  Note that the entries in W are unchanged.
!
                        DO K = NEWFST, NEWLST
                           FUDGE = THREE*EPS*ABS(WORK(WBEGIN+K-1))
                           WORK(WBEGIN+K-1) = WORK(WBEGIN+K-1) - TAU
                           FUDGE = FUDGE +                              &
     &                             FOUR*EPS*ABS(WORK(WBEGIN+K-1))
!
!  Fudge errors
!
                           WERR(WBEGIN+K-1) = WERR(WBEGIN+K-1) +        &
     &                                        FUDGE
!
!  Gaps are not fudged. Provided that WERR is small when eigenvalues
!  are close, a zero gap indicates that a new representation is needed
!  for resolving the cluster. A fudge could lead to a wrong decision
!  of judging eigenvalues 'separated' which in reality are not. This
!  could have a negative impact on the orthogonality of the computed
!  eigenvectors.
!
                        END DO

                        NCLUS = NCLUS + 1
                        K = NEWCLS + 2*NCLUS
                        IWORK(K-1) = NEWFST
                        IWORK(K) = NEWLST
                     ELSE
                        INFO = -2
                        RETURN
                     END IF
                  ELSE
!
!  Compute eigenvector of singleton
!
                     ITER = 0
!
                     TOL = FOUR * LOG(DBLE(IN)) * EPS
!
                     K = NEWFST
                     WINDEX = WBEGIN + K - 1
                     WINDMN = MAX(WINDEX - 1,1)
                     WINDPL = MIN(WINDEX + 1,M)
                     LAMBDA = WORK(WINDEX)
                     DONE = DONE + 1
!
!  Check if eigenvector computation is to be skipped.
!
                     IF ((WINDEX.LT.DOL).OR.(WINDEX.GT.DOU)) THEN
                        ESKIP = .TRUE.
                        GOTO 125
                     ELSE
                        ESKIP = .FALSE.
                     END IF
                     LEFT = WORK(WINDEX) - WERR(WINDEX)
                     RIGHT = WORK(WINDEX) + WERR(WINDEX)
                     INDEIG = INDEXW(WINDEX)
!
!  Note that since we compute the eigenpairs for a child,
!  all eigenvalue approximations are w.r.t the same shift.
!  In this case, the entries in WORK should be used for
!  computing the gaps since they exhibit even very small
!  differences in the eigenvalues, as opposed to the
!  entries in W which might "look" the same.

                     IF (K .EQ. 1) THEN
!
!  In the case RANGE='I' and with not much initial
!  accuracy in LAMBDA and VL, the formula
!  LGAP = MAX(ZERO, (SIGMA - VL) + LAMBDA)
!  can lead to an overestimation of the left gap and
!  thus to inadequately early RQI 'convergence'.
!  Prevent this by forcing a small left gap.
!
                        LGAP = EPS*MAX(ABS(LEFT),ABS(RIGHT))
                     ELSE
                        LGAP = WGAP(WINDMN)
                     END IF
                     IF (K .EQ. IM) THEN
!
!  In the case RANGE='I' and with not much initial
!  accuracy in LAMBDA and VU, the formula
!  can lead to an overestimation of the right gap and
!  thus to inadequately early RQI 'convergence'.
!  Prevent this by forcing a small right gap.
!
                        RGAP = EPS*MAX(ABS(LEFT),ABS(RIGHT))
                     ELSE
                        RGAP = WGAP(WINDEX)
                     END IF
                     GAP = MIN(LGAP, RGAP)
                     IF ((K .EQ. 1).OR.(K .EQ. IM)) THEN
!
!  The eigenvector support can become wrong
!  because significant entries could be cut off due to a
!  large GAPTOL parameter in LAR1V. Prevent this.
!
                        GAPTOL = ZERO
                     ELSE
                        GAPTOL = GAP * EPS
                     END IF
                     ISUPMN = IN
                     ISUPMX = 1
!
!  Update WGAP so that it holds the minimum gap
!  to the left or the right. This is crucial in the
!  case where bisection is used to ensure that the
!  eigenvalue is refined up to the required precision.
!  The correct value is restored afterwards.
!
                     SAVGAP = WGAP(WINDEX)
                     WGAP(WINDEX) = GAP
!
!  We want to use the Rayleigh Quotient Correction
!  as often as possible since it converges quadratically
!  when we are close enough to the desired eigenvalue.
!  However, the Rayleigh Quotient can have the wrong sign
!  and lead us away from the desired eigenvalue. In this
!  case, the best we can do is to use bisection.
!
                     USEDBS = .FALSE.
                     USEDRQ = .FALSE.
!
!  Bisection is initially turned off unless it is forced.
!
                     NEEDBS =  .NOT.TRYRQC
 120                 CONTINUE
!
!  Check if bisection should be used to refine eigenvalue.
!
                     IF (NEEDBS) THEN
!
!  Take the bisection as new iterate
!
                        USEDBS = .TRUE.
                        ITMP1 = IWORK(IINDR+WINDEX)
                        OFFSET = INDEXW(WBEGIN) - 1
                        CALL DLARRB (IN, D(IBEGIN),                     &
     &                               WORK(INDLLD+IBEGIN-1), INDEIG,     &
     &                               INDEIG, ZERO, TWO*EPS, OFFSET,     &
     &                               WORK(WBEGIN), WGAP(WBEGIN),        &
     &                               WERR(WBEGIN), WORK(INDWRK),        &
     &                               IWORK(IINDWK), PIVMIN, SPDIAM,     &
     &                               ITMP1, IINFO)
                        IF (IINFO.NE.0) THEN
                           INFO = -3
                           RETURN
                        END IF
                        LAMBDA = WORK(WINDEX)
!
!  Reset twist index from inaccurate LAMBDA to
!  force computation of true MINGMA.

                        IWORK(IINDR+WINDEX) = 0
                     END IF
!
!  Given LAMBDA, compute the eigenvector.
!
                     CALL DLAR1V (IN, 1, IN, LAMBDA, D(IBEGIN),         &
     &                            L(IBEGIN), WORK(INDLD+IBEGIN-1),      &
     &                            WORK(INDLLD+IBEGIN-1),                &
     &                            PIVMIN, GAPTOL, Z(IBEGIN, WINDEX),    &
     &                            .NOT.USEDBS, NEGCNT, ZTZ, MINGMA,     &
     &                            IWORK(IINDR+WINDEX),                  &
     &                            ISUPPZ(2*WINDEX-1),                   &
     &                            NRMINV, RESID, RQCORR, WORK(INDWRK))
                     IF (ITER .EQ. 0) THEN
                        BSTRES = RESID
                        BSTW = LAMBDA
                     ELSE IF (RESID.LT.BSTRES) THEN
                        BSTRES = RESID
                        BSTW = LAMBDA
                     END IF
                     ISUPMN = MIN(ISUPMN,ISUPPZ(2*WINDEX-1))
                     ISUPMX = MAX(ISUPMX,ISUPPZ(2*WINDEX))
                     ITER = ITER + 1
!
!  sin alpha <= |resid|/gap
!  Note that both the residual and the gap are
!  proportional to the matrix, so ||T|| does not play
!  a role in the quotient
!
!  Convergence test for Rayleigh-Quotient iteration
!  (omitted when Bisection has been used).
!
                     IF (RESID.GT.TOL*GAP .AND. ABS(RQCORR).GT.         &
     &                   RQTOL*ABS(LAMBDA) .AND. .NOT. USEDBS) THEN
!
!  We need to check that the RQCORR update does not move the eigenvalue
!  away from the desired one and towards a neighbor. -> protection with
!  bisection.
!
                        IF (INDEIG.LE.NEGCNT) THEN
!
!  The wanted eigenvalue lies to the left.
!
                           SGNDEF = -ONE
                        ELSE
!
!  The wanted eigenvalue lies to the right.
!
                           SGNDEF = ONE
                        END IF
!
!  We only use the RQCORR if it improves the iterate reasonably.
!
                        IF ((RQCORR*SGNDEF.GE.ZERO) .AND.               &
     &                      (LAMBDA + RQCORR.LE. RIGHT) .AND.           &
     &                      (LAMBDA + RQCORR.GE. LEFT)) THEN
                           USEDRQ = .TRUE.
!
!  Store new midpoint of bisection interval in WORK.
!
                           IF (SGNDEF.EQ.ONE) THEN
!
!  The current LAMBDA is on the left of the true eigenvalue.
!
                              LEFT = LAMBDA
!
!  We prefer to assume that the error estimate is correct. We could
!  make the interval not as a bracket but to be modified if the RQCORR
!  chooses to. In this case, the RIGHT side should be modified as
!  follows:
!
!                             RIGHT = MAX(RIGHT, LAMBDA + RQCORR)
                           ELSE
!
!  The current LAMBDA is on the right of the true eigenvalue.
!
                              RIGHT = LAMBDA
!
!  See comment about assuming the error estimate is correct above.
!
!                             LEFT = MIN(LEFT, LAMBDA + RQCORR)
                           END IF
                           WORK(WINDEX) = HALF * (RIGHT + LEFT)
!
!  Take RQCORR since it has the correct sign and improves the iterate
!  reasonably.
                           LAMBDA = LAMBDA + RQCORR
!
!  Update width of error interval.
!
                           WERR(WINDEX) = HALF * (RIGHT-LEFT)
                        ELSE
                           NEEDBS = .TRUE.
                        END IF
                        IF (RIGHT-LEFT.LT.RQTOL*ABS(LAMBDA)) THEN
!
!  The eigenvalue is computed to bisection accuracy compute eigenvector
!  and stop.
!
                           USEDBS = .TRUE.
                           GOTO 120
                        ELSE IF (ITER.LT.MAXITR) THEN
                           GOTO 120
                        ELSE IF (ITER.EQ.MAXITR) THEN
                           NEEDBS = .TRUE.
                           GOTO 120
                        ELSE
                           INFO = 5
                           RETURN
                        END IF
                     ELSE
                        STP2II = .FALSE.
                        IF (USEDRQ .AND. USEDBS .AND.                   &
     &                      BSTRES.LE.RESID) THEN
                           LAMBDA = BSTW
                           STP2II = .TRUE.
                        END IF
                        IF (STP2II) THEN
!
!  Improve error angle by second step.
!
                           CALL DLAR1V (IN, 1, IN, LAMBDA,              &
     &                                  D(IBEGIN), L(IBEGIN),           &
     &                                  WORK(INDLD+IBEGIN-1),           &
     &                                  WORK(INDLLD+IBEGIN-1),          &
     &                                  PIVMIN, GAPTOL,                 &
     &                                  Z(IBEGIN,WINDEX),               &
     &                                  .NOT.USEDBS, NEGCNT, ZTZ,       &
     &                                  MINGMA, IWORK(IINDR+WINDEX),    &
     &                                  ISUPPZ(2*WINDEX-1),             &
     &                                  NRMINV, RESID, RQCORR,          &
     &                                  WORK(INDWRK))
                        END IF
                        WORK(WINDEX) = LAMBDA
                     END IF
!
!  Compute FP-vector support w.r.t. whole matrix.
!
                     ISUPPZ(2*WINDEX-1) = ISUPPZ(2*WINDEX-1)+OLDIEN
                     ISUPPZ(2*WINDEX) = ISUPPZ(2*WINDEX)+OLDIEN
                     ZFROM = ISUPPZ(2*WINDEX-1)
                     ZTO = ISUPPZ(2*WINDEX)
                     ISUPMN = ISUPMN + OLDIEN
                     ISUPMX = ISUPMX + OLDIEN
!
!  Ensure vector is ok if support in the RQI has changed.
!
                     IF (ISUPMN.LT.ZFROM) THEN
                        DO  II = ISUPMN,ZFROM-1
                           Z(II, WINDEX) = ZERO
                        END DO
                     END IF
                     IF (ISUPMX.GT.ZTO) THEN
                        DO II = ZTO+1,ISUPMX
                           Z(II, WINDEX) = ZERO
                        END DO
                     END IF
                     CALL DSCAL (ZTO-ZFROM+1, NRMINV,                   &
     &                           Z(ZFROM,WINDEX), 1)
 125                 CONTINUE
!
!  Update W.
!
                     W(WINDEX) = LAMBDA+SIGMA
!
!  Recompute the gaps on the left and right. But only allow them to
!  become larger and not smaller (which can only happen through "bad"
!  cancellation and does not reflect the theory where the initial gaps
!  are underestimated due to WERR being too crude.)
!
                     IF (.NOT.ESKIP) THEN
                        IF (K.GT.1) THEN
                           WGAP(WINDMN) = MAX(WGAP(WINDMN),             &
     &                                        W(WINDEX)-WERR(WINDEX)-   &
     &                                        W(WINDMN)-WERR(WINDMN))
                        END IF
                        IF (WINDEX.LT.WEND) THEN
                           WGAP(WINDEX) = MAX(SAVGAP,                   &
     &                                        W(WINDPL)-WERR(WINDPL)-   &
     &                                        W(WINDEX)-WERR(WINDEX))
                        END IF
                     END IF
                     IDONE = IDONE + 1
                  END IF
!
!  Here ends the code for the current child.
!
 139              CONTINUE
!
!  Proceed to any remaining child nodes.
!
                  NEWFST = J + 1
 140           CONTINUE
 150        CONTINUE
            NDEPTH = NDEPTH + 1
            GOTO 40
         END IF
         IBEGIN = IEND + 1
         WBEGIN = WEND + 1
 170  CONTINUE
!
      RETURN
      END SUBROUTINE DLARRV
!
      SUBROUTINE DLARTG (F, G, CS, SN, R)
!
!=======================================================================
!                                                                      !
!  DLARTG generate a plane rotation so that                            !
!                                                                      !
!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.      !
!     [ -SN  CS  ]     [ G ]     [ 0 ]                                 !
!                                                                      !
!  This is a slower, more accurate version of the BLAS1 routine DROTG, !
!  with the following other differences:                               !
!                                                                      !
!     F and G are unchanged on return.                                 !
!     If G=0, then CS=1 and SN=0.                                      !
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any      !
!        floating point operations (saves work in DBDSQR when          !
!        there are zeros on the diagonal).                             !
!                                                                      !
!  If F exceeds G in magnitude, CS will be positive.                   !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  F       (input) DOUBLE PRECISION                                    !
!          The first component of vector to be rotated.                !
!                                                                      !
!  G       (input) DOUBLE PRECISION                                    !
!          The second component of vector to be rotated.               !
!                                                                      !
!  CS      (output) DOUBLE PRECISION                                   !
!          The cosine of the rotation.                                 !
!                                                                      !
!  SN      (output) DOUBLE PRECISION                                   !
!          The sine of the rotation.                                   !
!                                                                      !
!  R       (output) DOUBLE PRECISION                                   !
!          The nonzero component of the rotated vector.                !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     September 30, 1994                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      real (dp), intent(in)  :: F, G
      real (dp), intent(out) :: CS, R, SN
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
!
      logical, save   :: FIRST = .TRUE.
!
      integer         :: COUNT, I
!
      real (dp), save :: SAFMX2, SAFMIN, SAFMN2
      real (dp)       :: EPS, F1, G1, SCALE
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (FIRST) THEN
        FIRST = .FALSE.
        SAFMIN = DLAMCH('S')
        EPS = DLAMCH('E')
        SAFMN2 = DLAMCH('B')**INT(LOG(SAFMIN / EPS) /                   &
     &           LOG(DLAMCH('B')) / TWO)
        SAFMX2 = ONE / SAFMN2
      END IF
      IF (G .EQ. ZERO) THEN
        CS = ONE
        SN = ZERO
        R = F
      ELSE IF (F .EQ. ZERO) THEN
        CS = ZERO
        SN = ONE
        R = G
      ELSE
        F1 = F
        G1 = G
        SCALE = MAX(ABS(F1), ABS(G1))
        IF (SCALE .GE. SAFMX2) THEN
          COUNT = 0
#ifdef REMOVE_LAPACK_GOTOS
          LOOP_1 : DO                            ! forever DO-loop
            COUNT = COUNT + 1
            F1 = F1*SAFMN2
            G1 = G1*SAFMN2
            SCALE = MAX(ABS(F1), ABS(G1))
            IF (SCALE .GE. SAFMX2) CYCLE LOOP_1  ! iterate DO-loop again
            R = SQRT(F1**2 + G1**2)
            CS = F1 / R
            SN = G1 / R
            DO I = 1, COUNT
              R = R*SAFMX2
            END DO
            EXIT LOOP_1                          ! terminate DO-loop
          END DO LOOP_1
#else
   10     CONTINUE
          COUNT = COUNT + 1
          F1 = F1*SAFMN2
          G1 = G1*SAFMN2
          SCALE = MAX(ABS(F1), ABS(G1))
          IF (SCALE .GE. SAFMX2) GOTO 10
          R = SQRT(F1**2 + G1**2)
          CS = F1 / R
          SN = G1 / R
          DO I = 1, COUNT
            R = R*SAFMX2
          END DO
#endif
        ELSE IF (SCALE .LE. SAFMN2) THEN
          COUNT = 0
#ifdef REMOVE_LAPACK_GOTOS
          LOOP_2 : DO                            ! forever DO-loop
            COUNT = COUNT + 1
            F1 = F1*SAFMX2
            G1 = G1*SAFMX2
            SCALE = MAX(ABS(F1), ABS(G1))
            IF (SCALE .LE. SAFMN2) CYCLE LOOP_2  ! iterate DO-loop again
            R = SQRT(F1**2 + G1**2)
            CS = F1 / R
            SN = G1 / R
            DO I = 1, COUNT
              R = R*SAFMN2
            END DO
            EXIT LOOP_2                          ! terminate DO-loop
          END DO LOOP_2
#else
   30     CONTINUE
          COUNT = COUNT + 1
          F1 = F1*SAFMX2
          G1 = G1*SAFMX2
          SCALE = MAX(ABS(F1), ABS(G1))
          IF (SCALE .LE. SAFMN2) GOTO 30
          R = SQRT(F1**2 + G1**2)
          CS = F1 / R
          SN = G1 / R
          DO I = 1, COUNT
            R = R*SAFMN2
          END DO
#endif
        ELSE
          R = SQRT(F1**2 + G1**2)
          CS = F1 / R
          SN = G1 / R
        END IF
        IF ((ABS(F).GT.ABS(G)) .AND. (CS.LT.ZERO)) THEN
          CS = -CS
          SN = -SN
          R = -R
        END IF
      END IF
!
      RETURN
      END SUBROUTINE DLARTG
!
      SUBROUTINE DLARUV (ISEED, N, X)
!
!=======================================================================
!                                                                      !
!  DLARUV returns a vector of n random real numbers from a uniform     !
!  (0,1) distribution (n <= 128).                                      !
!                                                                      !
!  This is an auxiliary routine called by DLARNV and ZLARNV.           !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  ISEED    [in,out]  INTEGER array, dimension (4)                     !
!           On entry, the seed of the random number generator; the     !
!           array elements must be between 0 and 4095, and ISEED(4)    !
!           must be odd.                                               !
!           On exit, the seed is updated.                              !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The number of random numbers to be generated. N <= 128.    !
!                                                                      !
!  X        [out]     REAL array, dimension (N)                        !
!           The generated random numbers.                              !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  This routine uses a multiplicative congruential method with modulus !
!  2**48 and multiplier 33952834046453 (see G.S.Fishman,               !
!  'Multiplicative congruential random number generators with modulus  !
!  2**b: an exhaustive analysis for b = 32 and a partial analysis for  !
!  b = 48', Math. Comp. 189, pp 331-344, 1990).                        !
!                                                                      !
!  48-bit integers are stored in 4 integer array elements with 12 bits !
!  per element. Hence the routine is portable across machines with     !
!  integers of 32 bits or more.                                        !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: N
      integer,   intent(inout) :: ISEED(4)
!
      real (dp), intent(out)   :: X(N)
!
!  Local variable declarations.
!
      integer,   parameter :: LV = 128
      integer,   parameter :: IPW2 = 4096
!
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: R = ONE / IPW2
!
      integer :: I, I1, I2, I3, I4, IT1, IT2, IT3, IT4, J
      integer :: MM(LV, 4)
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      MM(  1, 1:4) = [ 494,  322, 2508, 2549]
      MM(  2, 1:4) = [2637,  789, 3754, 1145]
      MM(  3, 1:4) = [ 255, 1440, 1766, 2253]
      MM(  4, 1:4) = [2008,  752, 3572,  305]
      MM(  5, 1:4) = [1253, 2859, 2893, 3301]
      MM(  6, 1:4) = [3344,  123,  307, 1065]
      MM(  7, 1:4) = [4084, 1848, 1297, 3133]
      MM(  8, 1:4) = [1739,  643, 3966, 2913]
      MM(  9, 1:4) = [3143, 2405,  758, 3285]
      MM( 10, 1:4) = [3468, 2638, 2598, 1241]
      MM( 11, 1:4) = [ 688, 2344, 3406, 1197]
      MM( 12, 1:4) = [1657,   46, 2922, 3729]
      MM( 13, 1:4) = [1238, 3814, 1038, 2501]
      MM( 14, 1:4) = [3166,  913, 2934, 1673]
      MM( 15, 1:4) = [1292, 3649, 2091,  541]
      MM( 16, 1:4) = [3422,  339, 2451, 2753]
      MM( 17, 1:4) = [1270, 3808, 1580,  949]
      MM( 18, 1:4) = [2016,  822, 1958, 2361]
      MM( 19, 1:4) = [ 154, 2832, 2055, 1165]
      MM( 20, 1:4) = [2862, 3078, 1507, 4081]
      MM( 21, 1:4) = [ 697, 3633, 1078, 2725]
      MM( 22, 1:4) = [1706, 2970, 3273, 3305]
      MM( 23, 1:4) = [ 491,  637,   17, 3069]
      MM( 24, 1:4) = [ 931, 2249,  854, 3617]
      MM( 25, 1:4) = [1444, 2081, 2916, 3733]
      MM( 26, 1:4) = [ 444, 4019, 3971,  409]
      MM( 27, 1:4) = [3577, 1478, 2889, 2157]
      MM( 28, 1:4) = [3944,  242, 3831, 1361]
      MM( 29, 1:4) = [2184,  481, 2621, 3973]
      MM( 30, 1:4) = [1661, 2075, 1541, 1865]
      MM( 31, 1:4) = [3482, 4058,  893, 2525]
      MM( 32, 1:4) = [ 657,  622,  736, 1409]
      MM( 33, 1:4) = [3023, 3376, 3992, 3445]
      MM( 34, 1:4) = [3618,  812,  787, 3577]
      MM( 35, 1:4) = [1267,  234, 2125,   77]
      MM( 36, 1:4) = [1828,  641, 2364, 3761]
      MM( 37, 1:4) = [ 164, 4005, 2460, 2149]
      MM( 38, 1:4) = [3798, 1122,  257, 1449]
      MM( 39, 1:4) = [3087, 3135, 1574, 3005]
      MM( 40, 1:4) = [2400, 2640, 3912,  225]
      MM( 41, 1:4) = [2870, 2302, 1216,   85]
      MM( 42, 1:4) = [3876,   40, 3248, 3673]
      MM( 43, 1:4) = [1905, 1832, 3401, 3117]
      MM( 44, 1:4) = [1593, 2247, 2124, 3089]
      MM( 45, 1:4) = [1797, 2034, 2762, 1349]
      MM( 46, 1:4) = [1234, 2637,  149, 2057]
      MM( 47, 1:4) = [3460, 1287, 2245,  413]
      MM( 48, 1:4) = [ 328, 1691,  166,   65]
      MM( 49, 1:4) = [2861,  496,  466, 1845]
      MM( 50, 1:4) = [1950, 1597, 4018,  697]
      MM( 51, 1:4) = [ 617, 2394, 1399, 3085]
      MM( 52, 1:4) = [2070, 2584,  190, 3441]
      MM( 53, 1:4) = [3331, 1843, 2879, 1573]
      MM( 54, 1:4) = [ 769,  336,  153, 3689]
      MM( 55, 1:4) = [1558, 1472, 2320, 2941]
      MM( 56, 1:4) = [2412, 2407,   18,  929]
      MM( 57, 1:4) = [2800,  433,  712,  533]
      MM( 58, 1:4) = [ 189, 2096, 2159, 2841]
      MM( 59, 1:4) = [ 287, 1761, 2318, 4077]
      MM( 60, 1:4) = [2045, 2810, 2091,  721]
      MM( 61, 1:4) = [1227,  566, 3443, 2821]
      MM( 62, 1:4) = [2838,  442, 1510, 2249]
      MM( 63, 1:4) = [ 209,   41,  449, 2397]
      MM( 64, 1:4) = [2770, 1238, 1956, 2817]
      MM( 65, 1:4) = [3654, 1086, 2201,  245]
      MM( 66, 1:4) = [3993,  603, 3137, 1913]
      MM( 67, 1:4) = [ 192,  840, 3399, 1997]
      MM( 68, 1:4) = [2253, 3168, 1321, 3121]
      MM( 69, 1:4) = [3491, 1499, 2271,  997]
      MM( 70, 1:4) = [2889, 1084, 3667, 1833]
      MM( 71, 1:4) = [2857, 3438, 2703, 2877]
      MM( 72, 1:4) = [2094, 2408,  629, 1633]
      MM( 73, 1:4) = [1818, 1589, 2365,  981]
      MM( 74, 1:4) = [ 688, 2391, 2431, 2009]
      MM( 75, 1:4) = [1407,  288, 1113,  941]
      MM( 76, 1:4) = [ 634,   26, 3922, 2449]
      MM( 77, 1:4) = [3231,  512, 2554,  197]
      MM( 78, 1:4) = [ 815, 1456,  184, 2441]
      MM( 79, 1:4) = [3524,  171, 2099,  285]
      MM( 80, 1:4) = [1914, 1677, 3228, 1473]
      MM( 81, 1:4) = [ 516, 2657, 4012, 2741]
      MM( 82, 1:4) = [ 164, 2270, 1921, 3129]
      MM( 83, 1:4) = [ 303, 2587, 3452,  909]
      MM( 84, 1:4) = [2144, 2961, 3901, 2801]
      MM( 85, 1:4) = [3480, 1970,  572,  421]
      MM( 86, 1:4) = [ 119, 1817, 3309, 4073]
      MM( 87, 1:4) = [3357,  676, 3171, 2813]
      MM( 88, 1:4) = [ 837, 1410,  817, 2337]
      MM( 89, 1:4) = [2826, 3723, 3039, 1429]
      MM( 90, 1:4) = [2332, 2803, 1696, 1177]
      MM( 91, 1:4) = [2089, 3185, 1256, 1901]
      MM( 92, 1:4) = [3780,  184, 3715,   81]
      MM( 93, 1:4) = [1700,  663, 2077, 1669]
      MM( 94, 1:4) = [3712,  499, 3019, 2633]
      MM( 95, 1:4) = [ 150, 3784, 1497, 2269]
      MM( 96, 1:4) = [2000, 1631, 1101,  129]
      MM( 97, 1:4) = [3375, 1925,  717, 1141]
      MM( 98, 1:4) = [1621, 3912,   51,  249]
      MM( 99, 1:4) = [3090, 1398,  981, 3917]
      MM(100, 1:4) = [3765, 1349, 1978, 2481]
      MM(101, 1:4) = [1149, 1441, 1813, 3941]
      MM(102, 1:4) = [3146, 2224, 3881, 2217]
      MM(103, 1:4) = [  33, 2411,   76, 2749]
      MM(104, 1:4) = [3082, 1907, 3846, 3041]
      MM(105, 1:4) = [2741, 3192, 3694, 1877]
      MM(106, 1:4) = [ 359, 2786, 1682,  345]
      MM(107, 1:4) = [3316,  382,  124, 2861]
      MM(108, 1:4) = [1749,   37, 1660, 1809]
      MM(109, 1:4) = [ 185,  759, 3997, 3141]
      MM(110, 1:4) = [2784, 2948,  479, 2825]
      MM(111, 1:4) = [2202, 1862, 1141,  157]
      MM(112, 1:4) = [2199, 3802,  886, 2881]
      MM(113, 1:4) = [1364, 2423, 3514, 3637]
      MM(114, 1:4) = [1244, 2051, 1301, 1465]
      MM(115, 1:4) = [2020, 2295, 3604, 2829]
      MM(116, 1:4) = [3160, 1332, 1888, 2161]
      MM(117, 1:4) = [2785, 1832, 1836, 3365]
      MM(118, 1:4) = [2772, 2405, 1990,  361]
      MM(119, 1:4) = [1217, 3638, 2058, 2685]
      MM(120, 1:4) = [1822, 3661,  692, 3745]
      MM(121, 1:4) = [1245,  327, 1194, 2325]
      MM(122, 1:4) = [2252, 3660,   20, 3609]
      MM(123, 1:4) = [3904,  716, 3285, 3821]
      MM(124, 1:4) = [2774, 1842, 2046, 3537]
      MM(125, 1:4) = [ 997, 3987, 2107,  517]
      MM(126, 1:4) = [2573, 1368, 3508, 3017]
      MM(127, 1:4) = [1148, 1848, 3525, 2141]
      MM(128, 1:4) = [ 545, 2366, 3801, 1537]
!
!  Quick return for N < 1
!
      IF (N < 1) THEN
         RETURN
      END IF
!
      I1 = ISEED(1)
      I2 = ISEED(2)
      I3 = ISEED(3)
      I4 = ISEED(4)
!
      DO 10 I = 1, MIN(N, LV)
!
  20     CONTINUE
!
!  Multiply the seed by i-th power of the multiplier modulo 2**48
!
         IT4 = I4*MM(I, 4)
         IT3 = IT4 / IPW2
         IT4 = IT4 - IPW2*IT3
         IT3 = IT3 + I3*MM(I,4) + I4*MM(I,3)
         IT2 = IT3 / IPW2
         IT3 = IT3 - IPW2*IT2
         IT2 = IT2 + I2*MM(I,4) + I3*MM(I,3) + I4*MM(I,2)
         IT1 = IT2 / IPW2
         IT2 = IT2 - IPW2*IT1
         IT1 = IT1 + I1*MM(I,4) + I2*MM(I,3) + I3*MM(I,2) + I4*MM(I,1)
         IT1 = MOD(IT1, IPW2)
!
!  Convert 48-bit integer to a real number in the interval (0,1).
!
         X(I) = R*(REAL(IT1,dp)+                                        &
     &             R*(REAL(IT2,dp)+R*(REAL(IT3,dp)+R*REAL(IT4,dp))))
!
         IF (X(I).EQ.ONE) THEN
!
!  If a real number has n bits of precision, and the first
!  n bits of the 48-bit integer above happen to be all 1 (which
!  will occur about once every 2**n calls), then X(I) will
!  be rounded to exactly 1.0.
!  Since X(I) is not supposed to return exactly 0.0 or 1.0,
!  the statistically correct thing to do in this situation is
!  simply to iterate again.
!  N.B. the case X(I) = 0.0 should not be possible.
!
            I1 = I1 + 2
            I2 = I2 + 2
            I3 = I3 + 2
            I4 = I4 + 2
            GOTO 20
         END IF
!
   10 CONTINUE
!
!  Return final value of seed
!
      ISEED(1) = IT1
      ISEED(2) = IT2
      ISEED(3) = IT3
      ISEED(4) = IT4
!
      RETURN
      END SUBROUTINE DLARUV
!
      SUBROUTINE DLASCL (TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO)
!
!=======================================================================
!                                                                      !
!  DLASCL multiplies the M by N real matrix A by the real scalar       !
!  CTO/CFROM. This is done without over/underflow as long as the       !
!  final result CTO*A(I,J)/CFROM does not over/underflow. TYPE         !
!  specifies that A may be full, upper triangular, lower triangular,   !
!  upper Hessenberg, or banded.                                        !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  TYPE    (input) CHARACTER*1                                         !
!          TYPE indices the storage type of the input matrix.          !
!          = 'G':  A is a full matrix.                                 !
!          = 'L':  A is a lower triangular matrix.                     !
!          = 'U':  A is an upper triangular matrix.                    !
!          = 'H':  A is an upper Hessenberg matrix.                    !
!          = 'B':  A is a symmetric band matrix with lower bandwidth   !
!                  KL and upper bandwidth KU and with the only the     !
!                  lower half stored.                                  !
!          = 'Q':  A is a symmetric band matrix with lower bandwidth   !
!                  KL and upper bandwidth KU and with the only the     !
!                  upper half stored.                                  !
!          = 'Z':  A is a band matrix with lower bandwidth KL and      !
!                  upper bandwidth KU.                                 !
!                                                                      !
!  KL      (input) INTEGER                                             !
!          The lower bandwidth of A.  Referenced only if TYPE = 'B',   !
!          'Q' or 'Z'.                                                 !
!                                                                      !
!  KU      (input) INTEGER                                             !
!          The upper bandwidth of A.  Referenced only if TYPE = 'B',   !
!          'Q' or 'Z'.                                                 !
!                                                                      !
!  CFROM   (input) DOUBLE PRECISION                                    !
!  CTO     (input) DOUBLE PRECISION                                    !
!          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed !
!          without over/underflow if the final result CTO*A(I,J)/CFROM !
!          can be represented without over/underflow.  CFROM must be   !
!          nonzero.                                                    !
!                                                                      !
!  M       (input) INTEGER                                             !
!          The number of rows of the matrix A.  M >= 0.                !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The number of columns of the matrix A.  N >= 0.             !
!                                                                      !
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,M)    !
!          The matrix to be multiplied by CTO/CFROM.  See TYPE for the !
!          storage type.                                               !
!                                                                      !
!  LDA     (input) INTEGER                                             !
!          The leading dimension of the array A.  LDA >= max(1,M).     !
!                                                                      !
!  INFO    (output) INTEGER                                            !
!          0  - successful exit                                        !
!          <0 - if INFO = -i, the i-th argument had an illegal value.  !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     February 29, 1992                                                !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
!
      integer,   intent (in)   :: KL, KU, LDA, M, N
      integer,   intent (out)  :: INFO
!
      real (dp), intent(in)    :: CFROM, CTO
      real (dp), intent(inout) :: A(LDA, *)
!
      character (len=1), intent(in) :: TYPE
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      logical   :: DONE
!
      integer   :: I, ITYPE, J, K1, K2, K3, K4
!
      real (dp) :: BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input arguments.
!
      INFO = 0
!
      IF (LSAME(TYPE, 'G')) THEN
        ITYPE = 0
      ELSE IF (LSAME(TYPE, 'L')) THEN
        ITYPE = 1
      ELSE IF (LSAME(TYPE, 'U')) THEN
        ITYPE = 2
      ELSE IF (LSAME(TYPE, 'H')) THEN
        ITYPE = 3
      ELSE IF (LSAME(TYPE, 'B')) THEN
        ITYPE = 4
      ELSE IF (LSAME(TYPE, 'Q')) THEN
        ITYPE = 5
      ELSE IF (LSAME(TYPE, 'Z')) THEN
        ITYPE = 6
      ELSE
        ITYPE = -1
      END IF
!
      IF (ITYPE .EQ. -1) THEN
        INFO = -1
      ELSE IF (CFROM .EQ. ZERO) THEN
        INFO = -4
      ELSE IF (M .LT. 0) THEN
        INFO = -6
      ELSE IF (N.LT.0 .OR. (ITYPE.EQ.4 .AND. N.NE.M) .OR.               &
     &                     (ITYPE.EQ.5 .AND. N.NE.M)) THEN
        INFO = -7
      ELSE IF (ITYPE.LE.3 .AND. LDA.LT.MAX(1, M)) THEN
        INFO = -9
      ELSE IF (ITYPE .GE. 4) THEN
        IF (KL.LT.0 .OR. KL.GT.MAX(M-1, 0)) THEN
          INFO = -2
        ELSE IF (KU.LT.0 .OR. KU.GT.MAX(N-1, 0) .OR.                    &
     &           ((ITYPE.EQ.4 .OR. ITYPE.EQ.5) .AND. KL.NE.KU)) THEN
          INFO = -3
        ELSE IF ((ITYPE.EQ.4 .AND. LDA.LT.KL+1) .OR.                    &
     &           (ITYPE.EQ.5 .AND. LDA.LT.KU+1) .OR.                    &
     &           (ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1)) THEN
          INFO = -9
        END IF
      END IF
!
      IF (INFO .NE. 0) THEN
         CALL XERBLA ('DLASCL', -INFO)
         RETURN
      END IF
!
!  Quick return if possible.
!
      IF ((N.EQ.0) .OR. (M.EQ.0))  RETURN
!
!  Get machine parameters
!
      SMLNUM = DLAMCH ('S')
      BIGNUM = ONE / SMLNUM
!
      CFROMC = CFROM
      CTOC = CTO
!
      DONE = .FALSE.
!
      ITERATE : DO WHILE (.NOT. DONE)
        CFROM1 = CFROMC*SMLNUM
        CTO1 = CTOC / BIGNUM
        IF (ABS(CFROM1).GT.ABS(CTOC) .AND. CTOC.NE.ZERO) THEN
          MUL = SMLNUM
          DONE = .FALSE.
          CFROMC = CFROM1
        ELSE IF (ABS(CTO1) .GT. ABS(CFROMC)) THEN
          MUL = BIGNUM
          DONE = .FALSE.
          CTOC = CTO1
        ELSE
          MUL = CTOC / CFROMC
          DONE = .TRUE.
        END IF
!
        IF (ITYPE .EQ. 0) THEN
!
!  Full matrix.
!
          DO J = 1, N
            DO I = 1, M
               A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        ELSE IF (ITYPE .EQ. 1) THEN
!
!  Lower triangular matrix.
!
          DO J = 1, N
            DO I = J, M
              A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        ELSE IF (ITYPE.EQ.2) THEN
!
!  Upper triangular matrix.
!
          DO J = 1, N
            DO I = 1, MIN(J, M)
              A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        ELSE IF (ITYPE .EQ. 3) THEN
!
!  Upper Hessenberg matrix.
!
           DO J = 1, N
             DO I = 1, MIN(J+1, M)
               A(I,J) = A(I,J)*MUL
             END DO
           END DO
!
        ELSE IF (ITYPE .EQ. 4) THEN
!
!  Lower half of a symmetric band matrix.
!
          K3 = KL + 1
          K4 = N + 1
          DO J = 1, N
            DO I = 1, MIN(K3, K4-J)
              A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        ELSE IF (ITYPE .EQ. 5) THEN
!
!  Upper half of a symmetric band matrix.
!
          K1 = KU + 2
          K3 = KU + 1
          DO J = 1, N
            DO I = MAX(K1-J, 1), K3
              A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        ELSE IF (ITYPE .EQ. 6) THEN
!
!  Band matrix.
!
          K1 = KL + KU + 2
          K2 = KL + 1
          K3 = 2*KL + KU + 1
          K4 = KL + KU + 1 + M
          DO J = 1, N
            DO I = MAX(K1-J, K2), MIN(K3, K4-J)
              A(I,J) = A(I,J)*MUL
            END DO
          END DO
!
        END IF
      END DO ITERATE
!
      RETURN
      END SUBROUTINE DLASCL
!
      SUBROUTINE DLASET (UPLO, M, N, ALPHA, BETA, A, LDA)
!
!=======================================================================
!                                                                      !
!  DLASET initializes an m-by-n matrix A to BETA on the diagonal and   !
!  ALPHA on the offdiagonals.                                          !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  UPLO    (input) CHARACTER*1                                         !
!          Specifies the part of the matrix A to be set.               !
!          = 'U':      Upper triangular part is set; the strictly      !
!                      lower triangular part of A is not changed.      !
!          = 'L':      Lower triangular part is set; the strictly      !
!                      upper triangular part of A is not changed.      !
!          Otherwise:  All of the matrix A is set.                     !
!                                                                      !
!  M       (input) INTEGER                                             !
!          The number of rows of the matrix A.  M >= 0.                !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The number of columns of the matrix A.  N >= 0.             !
!                                                                      !
!  ALPHA   (input) DOUBLE PRECISION                                    !
!          The constant to which the offdiagonal elements are to be    !
!          set.                                                        !
!                                                                      !
!  BETA    (input) DOUBLE PRECISION                                    !
!          The constant to which the diagonal elements are to be set.  !
!                                                                      !
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)    !
!          On exit, the leading m-by-n submatrix of A is set as        !
!          follows:                                                    !
!                                                                      !
!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,          !
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,          !
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,    !
!                                                                      !
!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).           !
!                                                                      !
!  LDA     (input) INTEGER                                             !
!          The leading dimension of the array A.  LDA >= max(1,M).     !
!                                                                      !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: LDA, M, N
!
      real (dp), intent(in)    :: ALPHA, BETA
      real (dp), intent(inout) :: A(LDA, *)
!
      character (len=1), intent(in) :: UPLO
!
!  Local variable declarations.
!
      integer :: I, J
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (LSAME(UPLO, 'U')) THEN
!
!  Set the strictly upper triangular or trapezoidal part of the
!  array to ALPHA.
!
        DO J = 2, N
          DO I = 1, MIN(J-1, M)
            A(I,J) = ALPHA
          END DO
        END DO
!
      ELSE IF (LSAME(UPLO, 'L')) THEN
!
!  Set the strictly lower triangular or trapezoidal part of the
!  array to ALPHA.
!
        DO J = 1, MIN(M, N)
          DO I = J + 1, M
            A(I,J) = ALPHA
          END DO
        END DO
!
      ELSE
!
!  Set the leading m-by-n submatrix to ALPHA.
!
        DO J = 1, N
          DO I = 1, M
            A(I,J) = ALPHA
          END DO
        END DO
      END IF
!
!  Set the first MIN(M,N) diagonal elements to BETA.
!
      DO I = 1, MIN(M, N)
        A(I,I) = BETA
      END DO
!
      RETURN
      END SUBROUTINE DLASET
!
      SUBROUTINE DLASQ2 (N, Z, INFO)
!
!=======================================================================
!                                                                      !
!  DLASQ2 computes all the eigenvalues of the symmetric positive       !
!  definite tridiagonal matrix associated with the qd array Z to high  !
!  relative accuracy are computed to high relative accuracy, in the    !
!  absence of denormalization, underflow and overflow.                 !
!                                                                      !
!  To see the relation of Z to the tridiagonal matrix, let L be a      !
!  unit lower bidiagonal matrix with subdiagonals Z(2,4,6,,..) and     !
!  let U be an upper bidiagonal matrix with 1's above and diagonal     !
!  Z(1,3,5,,..). The tridiagonal is L*U or, if you prefer, the         !
!  symmetric tridiagonal to which it is similar.                       !
!                                                                      !
!  Note : DLASQ2 defines a logical variable, IEEE, which is true       !
!  on machines which follow ieee-754 floating-point standard in their  !
!  handling of infinities and NaNs, and false otherwise. This variable !
!  is passed to DLASQ3.                                                !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The number of rows and columns in the matrix, N >= 0.      !
!                                                                      !
!  Z        [in,out]  REAL array, dimension (4*N)                      !
!           On entry Z holds the qd array.                             !
!           On exit, entries 1 to N hold the eigenvalues in decreasing !
!           order, Z(2*N+1) holds the trace, and Z(2*N+2) holds the    !
!           sum of the eigenvalues. If N > 2, then Z(2*N+3) holds the  !
!           iteration count, Z(2*N+4) holds NDIVS/NIN^2, and Z(2*N+5)  !
!           holds the percentage of shifts that failed.                !
!                                                                      !
!  INFO     [out]     INTEGER                                          !
!            = 0: successful exit                                      !
!            < 0: if the i-th argument is a scalar and had an illegal  !
!                 value, then INFO = -i, if the i-th argument is an    !
!                 array and the j-entry had an illegal value, then     !
!                 INFO = -(i*100+j).                                   !
!            > 0: the algorithm failed                                 !
!                 = 1, a split was marked by a positive value in E     !
!                 = 2, current block of Z not diagonalized after 100*N !
!                      iterations (in inner while loop). On exit Z     !
!                      holds a qd array with the same eigenvalues as   !
!                      the given Z.                                    !
!                 = 3, termination criterion of outer while loop not   !
!                      met (program created more than N unreduced      !
!                      blocks).                                        !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  Local Variables: I0:N0 defines a current unreduced segment of Z.    !
!  The shifts are accumulated in SIGMA. Iteration count is in ITER.    !
!  Ping-pong is controlled by PP (alternates between 0 and 1).         !
!                                                                      !
!  -- LAPACK computational routine --                                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: N
      integer,   intent(out)   :: INFO
!
      real (dp), intent(inout) :: Z(*)
!
!  Local variable declarations.
!
      real (dp), parameter :: CBIAS = 1.50_dp
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: HALF = 0.5_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: FOUR = 4.0_dp
      real (dp), parameter :: HUNDRD = 100.0_dp
!
      logical   :: IEEE
!
      integer   :: I0, I1, I4, IINFO, IPN4, ITER, IWHILA, IWHILB
      integer   :: K, KMIN, N0, N1, NBIG, NDIV, NFAIL, PP, SPLT
      integer   :: TTYPE
!
      real (dp) :: D, DEE, DEEMIN, DESIG, DMIN, DMIN1, DMIN2, DN
      real (dp) :: DN1, DN2, E, EMAX, EMIN, EPS, G, OLDEMN, QMAX
      real (dp) :: QMIN, S, SAFMIN, SIGMA, T, TAU, TEMP, TOL
      real (dp) :: TOL2, TRACE, ZMAX, TEMPE, TEMPQ
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input arguments.
!  (in case DLASQ2 is not called by DLASQ1)
!
      INFO = 0
      EPS = DLAMCH('Precision')
      SAFMIN = DLAMCH('Safe minimum')
      TOL = EPS*HUNDRD
      TOL2 = TOL**2
!
      IF (N.LT.0) THEN
         INFO = -1
         CALL XERBLA ('DLASQ2', 1)
         RETURN
      ELSE IF (N.EQ.0) THEN
         RETURN
      ELSE IF (N.EQ.1) THEN
!
!  1-by-1 case.
!
         IF (Z(1).LT.ZERO) THEN
            INFO = -201
            CALL XERBLA ('DLASQ2', 2)
         END IF
         RETURN
      ELSE IF (N.EQ.2) THEN
!
!  2-by-2 case.
!
         IF (Z(1).LT.ZERO) THEN
            INFO = -201
            CALL XERBLA ('DLASQ2', 2)
            RETURN
         ELSE IF (Z(2).LT.ZERO) THEN
            INFO = -202
            CALL XERBLA ('DLASQ2', 2)
            RETURN
         ELSE IF (Z(3).LT.ZERO) THEN
           INFO = -203
           CALL XERBLA ('DLASQ2', 2)
           RETURN
         ELSE IF (Z(3).GT.Z(1)) THEN
            D = Z(3)
            Z(3) = Z(1)
            Z(1) = D
         END IF
         Z(5) = Z(1) + Z(2) + Z(3)
         IF (Z(2).GT.Z(3)*TOL2) THEN
            T = HALF*((Z(1)-Z(3))+Z(2))
            S = Z(3)*(Z(2) / T)
            IF (S.LE.T) THEN
               S = Z(3)*(Z(2) / (T*(ONE+SQRT(ONE+S/T))))
            ELSE
               S = Z(3)*(Z(2) / (T+SQRT(T)*SQRT(T+S)))
            END IF
            T = Z(1) + (S+Z(2))
            Z(3) = Z(3)*(Z(1) / T)
            Z(1) = T
         END IF
         Z(2) = Z(3)
         Z(6) = Z(2) + Z(1)
         RETURN
      END IF
!
!  Check for negative data and compute sums of q and e.
!
      Z(2*N) = ZERO
      EMIN = Z(2)
      QMAX = ZERO
      ZMAX = ZERO
      D = ZERO
      E = ZERO
!
      DO K = 1, 2*(N-1), 2
         IF (Z(K).LT.ZERO) THEN
            INFO = -(200+K)
            CALL XERBLA ('DLASQ2', 2)
            RETURN
         ELSE IF (Z(K+1).LT.ZERO) THEN
            INFO = -(200+K+1)
            CALL XERBLA ('DLASQ2', 2)
            RETURN
         END IF
         D = D + Z(K)
         E = E + Z(K+1)
         QMAX = MAX(QMAX, Z(K))
         EMIN = MIN(EMIN, Z(K+1))
         ZMAX = MAX(QMAX, ZMAX, Z(K+1))
      END DO
      IF (Z(2*N-1).LT.ZERO) THEN
         INFO = -(200+2*N-1)
         CALL XERBLA ('DLASQ2', 2)
         RETURN
      END IF
      D = D + Z(2*N-1)
      QMAX = MAX(QMAX, Z(2*N-1))
      ZMAX = MAX(QMAX, ZMAX)
!
!  Check for diagonality.
!
      IF (E.EQ.ZERO) THEN
         DO K = 2, N
            Z(K) = Z(2*K-1)
         END DO
         CALL DLASRT ('D', N, Z, IINFO)
         Z(2*N-1) = D
         RETURN
      END IF
!
      TRACE = D + E
!
!  Check for zero data.
!
      IF (TRACE.EQ.ZERO) THEN
         Z(2*N-1) = ZERO
         RETURN
      END IF
!
!  Check whether the machine is IEEE conformable.
!
      IEEE = (ILAENV(10, 'DLASQ2', 'N', 1, 2, 3, 4).EQ.1)
!
!  Rearrange data for locality: Z=(q1,qq1,e1,ee1,q2,qq2,e2,ee2,...).
!
      DO K = 2*N, 2, -2
         Z(2*K) = ZERO
         Z(2*K-1) = Z(K)
         Z(2*K-2) = ZERO
         Z(2*K-3) = Z(K-1)
      END DO
!
      I0 = 1
      N0 = N
!
!  Reverse the qd-array, if warranted.
!
      IF (CBIAS*Z(4*I0-3).LT.Z(4*N0-3)) THEN
         IPN4 = 4*(I0+N0)
         DO I4 = 4*I0, 2*(I0+N0-1), 4
            TEMP = Z(I4-3)
            Z(I4-3) = Z(IPN4-I4-3)
            Z(IPN4-I4-3) = TEMP
            TEMP = Z(I4-1)
            Z(I4-1) = Z(IPN4-I4-5)
            Z(IPN4-I4-5) = TEMP
         END DO
      END IF
!
!  Initial split checking via dqd and Li test.
!
      PP = 0
!
      DO K = 1, 2

         D = Z(4*N0+PP-3)
         DO I4 = 4*(N0-1) + PP, 4*I0 + PP, -4
            IF (Z(I4-1).LE.TOL2*D) THEN
               Z(I4-1) = -ZERO
               D = Z(I4-3)
            ELSE
               D = Z(I4-3)*(D / (D+Z(I4-1)))
            END IF
         END DO
!
!  dqd maps Z to ZZ plus Li's test.
!
         EMIN = Z(4*I0+PP+1)
         D = Z(4*I0+PP-3)
         DO I4 = 4*I0 + PP, 4*(N0-1) + PP, 4
            Z(I4-2*PP-2) = D + Z(I4-1)
            IF (Z(I4-1).LE.TOL2*D) THEN
               Z(I4-1) = -ZERO
               Z(I4-2*PP-2) = D
               Z(I4-2*PP) = ZERO
               D = Z(I4+1)
            ELSE IF (SAFMIN*Z(I4+1).LT.Z(I4-2*PP-2) .AND.               &
     &               SAFMIN*Z(I4-2*PP-2).LT.Z(I4+1)) THEN
               TEMP = Z(I4+1) / Z(I4-2*PP-2)
               Z(I4-2*PP) = Z(I4-1)*TEMP
               D = D*TEMP
            ELSE
               Z(I4-2*PP) = Z(I4+1)*(Z(I4-1) / Z(I4-2*PP-2))
               D = Z(I4+1)*(D / Z(I4-2*PP-2))
            END IF
            EMIN = MIN(EMIN, Z(I4-2*PP))
         END DO
         Z(4*N0-PP-2) = D
!
!  Now find qmax.
!
         QMAX = Z(4*I0-PP-2)
         DO I4 = 4*I0 - PP + 2, 4*N0 - PP - 2, 4
            QMAX = MAX(QMAX, Z(I4))
         END DO
!
!  Prepare for the next iteration on K.
!
         PP = 1 - PP
      END DO
!
!  Initialise variables to pass to DLASQ3.
!
      TTYPE = 0
      DMIN1 = ZERO
      DMIN2 = ZERO
      DN    = ZERO
      DN1   = ZERO
      DN2   = ZERO
      G     = ZERO
      TAU   = ZERO
!
      ITER = 2
      NFAIL = 0
      NDIV = 2*(N0-I0)
!
      DO 160 IWHILA = 1, N + 1
         IF (N0.LT.1) GOTO 170
!
!  While array unfinished do
!
!  E(N0) holds the value of SIGMA when submatrix in I0:N0
!   splits from the rest of the array, but is negated.
!
         DESIG = ZERO
         IF (N0.EQ.N) THEN
            SIGMA = ZERO
         ELSE
            SIGMA = -Z(4*N0-1)
         END IF
         IF (SIGMA.LT.ZERO) THEN
            INFO = 1
            RETURN
         END IF
!
!  Find last unreduced submatrix's top index I0, find QMAX and
!  EMIN. Find Gershgorin-type bound if Q's much greater than E.
!
         EMAX = ZERO
         IF (N0.GT.I0) THEN
            EMIN = ABS(Z(4*N0-5))
         ELSE
            EMIN = ZERO
         END IF
         QMIN = Z(4*N0-3)
         QMAX = QMIN
         DO I4 = 4*N0, 8, -4
            IF (Z(I4-5).LE.ZERO)  GOTO 100
            IF (QMIN.GE.FOUR*EMAX) THEN
               QMIN = MIN(QMIN, Z(I4-3))
               EMAX = MAX(EMAX, Z(I4-5))
            END IF
            QMAX = MAX(QMAX, Z(I4-7)+Z(I4-5))
            EMIN = MIN(EMIN, Z(I4-5))
         END DO
         I4 = 4
!
  100    CONTINUE
         I0 = I4 / 4
         PP = 0
!
         IF (N0-I0.GT.1) THEN
            DEE = Z(4*I0-3)
            DEEMIN = DEE
            KMIN = I0
            DO I4 = 4*I0+1, 4*N0-3, 4
               DEE = Z(I4)*(DEE /(DEE+Z(I4-2)))
               IF (DEE.LE.DEEMIN) THEN
                  DEEMIN = DEE
                  KMIN = (I4+3)/4
               END IF
            END DO
            IF ((KMIN-I0)*2.LT.N0-KMIN .AND.                            &
     &          DEEMIN.LE.HALF*Z(4*N0-3)) THEN
               IPN4 = 4*(I0+N0)
               PP = 2
               DO I4 = 4*I0, 2*(I0+N0-1), 4
                  TEMP = Z(I4-3)
                  Z(I4-3) = Z(IPN4-I4-3)
                  Z(IPN4-I4-3) = TEMP
                  TEMP = Z(I4-2)
                  Z(I4-2) = Z(IPN4-I4-2)
                  Z(IPN4-I4-2) = TEMP
                  TEMP = Z(I4-1)
                  Z(I4-1) = Z(IPN4-I4-5)
                  Z(IPN4-I4-5) = TEMP
                  TEMP = Z(I4)
                  Z(I4) = Z(IPN4-I4-4)
                  Z(IPN4-I4-4) = TEMP
               END DO
            END IF
         END IF
!
!  Put -(initial shift) into DMIN.
!
         DMIN = -MAX(ZERO, QMIN-TWO*SQRT(QMIN)*SQRT(EMAX))
!
!  Now I0:N0 is unreduced.
!    PP = 0 for ping, PP = 1 for pong.
!    PP = 2 indicates that flipping was applied to the Z array and
!           and that the tests for deflation upon entry in DLASQ3
!           should not be performed.
!
         NBIG = 100*(N0-I0+1)
         DO IWHILB = 1, NBIG
            IF (I0.GT.N0) GOTO 150
!
!  While submatrix unfinished take a good dqds step.
!
            CALL DLASQ3 (I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX,       &
     &                   NFAIL,                                         &
     &                   ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2,         &
     &                   DN, DN1, DN2, G, TAU)
!
            PP = 1 - PP
!
!  When EMIN is very small check for splits.
!
            IF (PP.EQ.0 .AND. N0-I0.GE.3) THEN
               IF (Z(4*N0).LE.TOL2*QMAX .OR.                            &
     &             Z(4*N0-1).LE.TOL2*SIGMA) THEN
                  SPLT = I0 - 1
                  QMAX = Z(4*I0-3)
                  EMIN = Z(4*I0-1)
                  OLDEMN = Z(4*I0)
                  DO I4 = 4*I0, 4*(N0-3), 4
                     IF (Z(I4).LE.TOL2*Z(I4-3) .OR.                     &
     &                   Z(I4-1).LE.TOL2*SIGMA) THEN
                        Z(I4-1) = -SIGMA
                        SPLT = I4 / 4
                        QMAX = ZERO
                        EMIN = Z(I4+3)
                        OLDEMN = Z(I4+4)
                     ELSE
                        QMAX = MAX(QMAX, Z(I4+1))
                        EMIN = MIN(EMIN, Z(I4-1))
                        OLDEMN = MIN(OLDEMN, Z(I4))
                     END IF
                  END DO
                  Z(4*N0-1) = EMIN
                  Z(4*N0) = OLDEMN
                  I0 = SPLT + 1
               END IF
            END IF
!
         END DO
!
         INFO = 2
!
!  Maximum number of iterations exceeded, restore the shift
!  SIGMA and place the new d's and e's in a qd array.
!  This might need to be done for several blocks
!
         I1 = I0
         N1 = N0
 145     CONTINUE
         TEMPQ = Z(4*I0-3)
         Z(4*I0-3) = Z(4*I0-3) + SIGMA
         DO K = I0+1, N0
            TEMPE = Z(4*K-5)
            Z(4*K-5) = Z(4*K-5) * (TEMPQ / Z(4*K-7))
            TEMPQ = Z(4*K-3)
            Z(4*K-3) = Z(4*K-3) + SIGMA + TEMPE - Z(4*K-5)
         END DO
!
!  Prepare to do this on the previous block if there is one
!
         IF (I1.GT.1) THEN
            N1 = I1-1
            DO WHILE ((I1.GE.2) .AND. (Z(4*I1-5).GE.ZERO))
               I1 = I1 - 1
            END DO
            SIGMA = -Z(4*N1-1)
            GOTO 145
         END IF

         DO K = 1, N
            Z(2*K-1) = Z(4*K-3)
!
!  Only the block 1..N0 is unfinished.  The rest of the e
!  must be essentially zero, although sometimes other data
!  has been stored in them.
!
            IF (K.LT.N0) THEN
               Z(2*K) = Z(4*K-1)
            ELSE
               Z(2*K) = 0
            END IF
         END DO
         RETURN
!
!  End IWHILB
!
  150    CONTINUE
!
  160 CONTINUE
!
      INFO = 3
      RETURN
!
!  End IWHILA
!
  170 CONTINUE
!
!  Move q to the front.
!
      DO K = 2, N
         Z(K) = Z(4*K-3)
      END DO
!
!  Sort and compute sum of eigenvalues.
!
      CALL DLASRT ('D', N, Z, IINFO)
!
      E = ZERO
      DO K = N, 1, -1
         E = E + Z(K)
      END DO
!
!  Store trace, sum(eigenvalues) and information on performance.
!
      Z(2*N+1) = TRACE
      Z(2*N+2) = E
      Z(2*N+3) = DBLE(ITER)
      Z(2*N+4) = DBLE(NDIV) / DBLE(N**2)
      Z(2*N+5) = HUNDRD*NFAIL / DBLE(ITER)
!
      RETURN
      END SUBROUTINE DLASQ2
!
      SUBROUTINE DLASQ3 (I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX,       &
     &                   NFAIL,                                         &
     &                   ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2,         &
     &                   DN, DN1, DN2, G, TAU)
!
!=======================================================================
!                                                                      !
!  DLASQ3 checks for deflation, computes a shift (TAU) and calls dqds. !
!  In case of failure it changes shifts, and tries again until output  !
!  is positive.                                                        !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  I0       [in]      INTEGER                                          !
!           First index.                                               !
!                                                                      !
!  N0       [in,out]  INTEGER                                          !
!           Last index.                                                !
!                                                                      !
!  Z        [in,out]  REAL array, dimension (4*N0)                     !
!           Z holds the qd array.                                      !
!                                                                      !
!  PP       [in,out]  INTEGER                                          !
!           PP=0 for ping, PP=1 for pong.                              !
!           PP=2 indicates that flipping was applied to the Z array    !
!                and that the initial tests for deflation should not   !
!                beperformed.                                          !
!                                                                      !
!  DIM      [out]     REAL                                             !
!           Minimum value of d.                                        !
!                                                                      !
!  SIGMA    [out]     REAL                                             !
!           Sum of shifts used in current segment.                     !
!                                                                      !
!  DESIG    [in,out]  REAL                                             !
!           Lower order part of SIGMA.                                 !
!                                                                      !
!  QMAX     [in,out]  REAL                                             !
!           Maximum value of q.                                        !
!                                                                      !
!  NFAIL    [in,out]  INTEGER                                          !
!           Increment NFAIL by 1 each time the shift was too big.      !
!                                                                      !
!  ITER     [in,out]  INTEGER                                          !
!           Increment ITER by 1 for each iteration.                    !
!                                                                      !
!  NDIV     [in,out]  INTEGER                                          !
!           Increment NDIV by 1 for each division.                     !
!                                                                      !
!  IEEE     [in]      LOGICAL                                          !
!           Flag for IEEE or non IEEE arithmetic (passed to DLASQ5).   !
!                                                                      !
!  TTYPE    [in,out]  INTEGER                                          !
!           Shift type.                                                !
!                                                                      !
!  DMIN1    [in,out]  REAL                                             !
!                                                                      !
!  DMIN2    [in,out]  REAL                                             !
!                                                                      !
!  DN       [in,out]  REAL                                             !
!                                                                      !
!  DN1      [in,out]  REAL                                             !
!                                                                      !
!  DN2      [in,out]  REAL                                             !
!                                                                      !
!  G        [in,out]  REAL                                             !
!                                                                      !
!  TAU      [in,out]  REAL                                             !
!           These are passed as arguments in order to save their       !
!           values between calls to DLASQ3.                            !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  -- LAPACK computational routine --                                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical,   intent(in)    :: IEEE
!
      integer,   intent(in)    :: I0
      integer,   intent(inout) :: ITER, N0, NDIV, NFAIL, PP, TTYPE
!
      real (dp), intent(inout) :: DESIG, DMIN1, DMIN2, DN, DN1, DN2, G
      real (dp), intent(inout) :: QMAX, SIGMA, TAU
      real (dp), intent(inout) :: Z(*)
      real (dp), intent(out)   :: DMIN
!
!  Local variable declarations.
!
      real (dp), parameter :: CBIAS = 1.50_dp
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: QURTR = 0.25_dp
      real (dp), parameter :: HALF = 0.5_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: HUNDRD = 100.0_dp
!
      integer   :: IPN4, J4, N0IN, NN
      real (dp) :: EPS, S, T, TEMP, TOL, TOL2
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      N0IN = N0
      EPS = DLAMCH('Precision')
      TOL = EPS*HUNDRD
      TOL2 = TOL**2
!
!  Check for deflation.
!
   10 CONTINUE
!
      IF (N0.LT.I0) RETURN
      IF (N0.EQ.I0) GOTO 20
      NN = 4*N0 + PP
      IF (N0.EQ.(I0+1)) GOTO 40
!
!  Check whether E(N0-1) is negligible, 1 eigenvalue.
!
      IF (Z(NN-5).GT.TOL2*(SIGMA+Z(NN-3)) .AND.                         &
     &    Z(NN-2*PP-4).GT.TOL2*Z(NN-7)) GOTO 30
!
   20 CONTINUE
!
      Z(4*N0-3) = Z(4*N0+PP-3) + SIGMA
      N0 = N0 - 1
      GOTO 10
!
!  Check  whether E(N0-2) is negligible, 2 eigenvalues.
!
   30 CONTINUE
!
      IF (Z(NN-9).GT.TOL2*SIGMA .AND.                                   &
     &    Z(NN-2*PP-8).GT.TOL2*Z(NN-11)) GOTO 50
!
   40 CONTINUE
!
      IF (Z(NN-3).GT.Z(NN-7)) THEN
         S = Z(NN-3)
         Z(NN-3) = Z(NN-7)
         Z(NN-7) = S
      END IF
      T = HALF*((Z(NN-7)-Z(NN-3))+Z(NN-5))
      IF (Z(NN-5).GT.Z(NN-3)*TOL2 .AND. T.NE.ZERO) THEN
         S = Z(NN-3)*(Z(NN-5) / T)
         IF (S.LE.T) THEN
            S = Z(NN-3)*(Z(NN-5) /                                      &
     &                   (T*(ONE+SQRT(ONE+S / T))))
         ELSE
            S = Z(NN-3)*(Z(NN-5) / (T+SQRT(T)*SQRT(T+S)))
         END IF
         T = Z(NN-7) + (S+Z(NN-5))
         Z(NN-3) = Z(NN-3)*(Z(NN-7) / T)
         Z(NN-7) = T
      END IF
      Z(4*N0-7) = Z(NN-7) + SIGMA
      Z(4*N0-3) = Z(NN-3) + SIGMA
      N0 = N0 - 2
      GOTO 10
!
   50 CONTINUE
      IF (PP.EQ.2) PP = 0
!
!  Reverse the qd-array, if warranted.
!
      IF (DMIN.LE.ZERO .OR. N0.LT.N0IN) THEN
         IF (CBIAS*Z(4*I0+PP-3) .LT. Z(4*N0+PP-3)) THEN
            IPN4 = 4*(I0+N0)
            DO J4 = 4*I0, 2*(I0+N0-1), 4
               TEMP = Z(J4-3)
               Z(J4-3) = Z(IPN4-J4-3)
               Z(IPN4-J4-3) = TEMP
               TEMP = Z(J4-2)
               Z(J4-2) = Z(IPN4-J4-2)
               Z(IPN4-J4-2) = TEMP
               TEMP = Z(J4-1)
               Z(J4-1) = Z(IPN4-J4-5)
               Z(IPN4-J4-5) = TEMP
               TEMP = Z(J4)
               Z(J4) = Z(IPN4-J4-4)
               Z(IPN4-J4-4) = TEMP
            END DO
            IF (N0-I0.LE.4) THEN
               Z(4*N0+PP-1) = Z(4*I0+PP-1)
               Z(4*N0-PP) = Z(4*I0-PP)
            END IF
            DMIN2 = MIN(DMIN2, Z(4*N0+PP-1))
            Z(4*N0+PP-1) = MIN(Z(4*N0+PP-1), Z(4*I0+PP-1),              &
     &                         Z(4*I0+PP+3))
            Z(4*N0-PP) = MIN(Z(4*N0-PP), Z(4*I0-PP),                    &
     &                       Z(4*I0-PP+4))
            QMAX = MAX(QMAX, Z(4*I0+PP-3), Z(4*I0+PP+1))
            DMIN = -ZERO
         END IF
      END IF
!
!  Choose a shift.
!
      CALL DLASQ4 (I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN, DN1,    &
     &             DN2, TAU, TTYPE, G)
!
!  Call dqds until DMIN > 0.
!
   70 CONTINUE
!
      CALL DLASQ5 (I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1, DMIN2, DN,   &
     &             DN1, DN2, IEEE, EPS)
!
      NDIV = NDIV + (N0-I0+2)
      ITER = ITER + 1
!
!  Check status.
!
      IF (DMIN.GE.ZERO .AND. DMIN1.GE.ZERO) THEN
!
!  Success.
!
         GOTO 90
!
      ELSE IF (DMIN.LT.ZERO .AND. DMIN1.GT.ZERO .AND.                   &
     &         Z(4*(N0-1)-PP).LT.TOL*(SIGMA+DN1) .AND.                  &
     &         ABS(DN).LT.TOL*SIGMA) THEN
!
!  Convergence hidden by negative DN.
!
         Z(4*(N0-1)-PP+2) = ZERO
         DMIN = ZERO
         GOTO 90
      ELSE IF (DMIN.LT.ZERO) THEN
!
!  TAU too big. Select new TAU and try again.
!
         NFAIL = NFAIL + 1
         IF (TTYPE.LT.-22) THEN
!
!  Failed twice. Play it safe.
!
            TAU = ZERO
         ELSE IF (DMIN1.GT.ZERO) THEN
!
!  Late failure. Gives excellent shift.
!
            TAU = (TAU+DMIN)*(ONE-TWO*EPS)
            TTYPE = TTYPE - 11
         ELSE
!
!  Early failure. Divide by 4.
!
            TAU = QURTR*TAU
            TTYPE = TTYPE - 12
         END IF
         GOTO 70
      ELSE IF (DISNAN(DMIN)) THEN
!
!  NaN.
!
         IF (TAU.EQ.ZERO) THEN
            GOTO 80
         ELSE
            TAU = ZERO
            GOTO 70
         END IF
      ELSE
!
!  Possible underflow. Play it safe.
!
         GOTO 80
      END IF
!
!  Risk of underflow.
!
   80 CONTINUE
      CALL DLASQ6 (I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN, DN1, DN2)
      NDIV = NDIV + (N0-I0+2)
      ITER = ITER + 1
      TAU = ZERO
!
   90 CONTINUE
      IF (TAU.LT.SIGMA) THEN
         DESIG = DESIG + TAU
         T = SIGMA + DESIG
         DESIG = DESIG - (T-SIGMA)
      ELSE
         T = SIGMA + TAU
         DESIG = SIGMA - (T-TAU) + DESIG
      END IF
      SIGMA = T
!
      RETURN
      END SUBROUTINE DLASQ3
!
      SUBROUTINE DLASQ4 (I0, N0, Z, PP, N0IN, DMIN, DMIN1, DMIN2, DN,   &
     &                   DN1, DN2, TAU, TTYPE, G)
!
!=======================================================================
!                                                                      !
!  DLASQ4 computes an approximation TAU to the smallest eigenvalue     !
!  using values of d from the previous transform.                      !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  I0       [in]      INTEGER                                          !
!           First index.                                               !
!                                                                      !
!  N0       [in]      INTEGER                                          !
!           Last index.                                                !
!                                                                      !
!  Z        [in]      REAL array, dimension (4*N0)                     !
!           Z holds the qd array.                                      !
!                                                                      !
!  PP       [in]      INTEGER                                          !
!           PP=0 for ping, PP=1 for pong.                              !
!                                                                      !
!  N0IN     [in]      INTEGER                                          !
!           The value of N0 at start of EIGTEST.                       !
!                                                                      !
!  DMIN     [in]      REAL                                             !
!           Minimum value of d.                                        !
!                                                                      !
!  DMIN1    [in]      REAL                                             !
!           Minimum value of d, excluding D(N0).                       !
!                                                                      !
!  DMIN2    [in]      REAL                                             !
!           Minimum value of d, excluding D(N0) and D(N0-1).           !
!                                                                      !
!  DN       [in]      REAL                                             !
!           d(N)                                                       !
!                                                                      !
!  DN1      [in]      REAL                                             !
!           d(N-1)                                                     !
!                                                                      !
!  DN2      [in]      REAL                                             !
!           d(N-2)                                                     !
!                                                                      !
!  TAU      [out]     REAL                                             !
!           This is the shift.                                         !
!                                                                      !
!  TTYPE    [out]     INTEGER                                          !
!           Shift type.                                                !
!                                                                      !
!  G        [in,out]  REAL                                             !
!           G is passed as an argument in order to save its value      !
!           between calls to DLASQ4.                                   !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!   CNST1 = 9/16                                                       !
!                                                                      !
!  -- LAPACK computational routine --                                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: I0, N0, N0IN, PP
      integer,   intent(inout) :: TTYPE
!
      real (dp), intent(in)    :: DMIN, DMIN1, DMIN2, DN, DN1, DN2
      real (dp), intent(in)    :: Z(*)
      real (dp), intent(inout) :: G
      real (dp), intent(out)   :: TAU
!
!  Local variable declarations.
!
      real (dp), parameter :: CNST1 = 0.5630_dp
      real (dp), parameter :: CNST2 = 1.010_dp
      real (dp), parameter :: CNST3 = 1.050_dp
!
      real (dp), parameter :: QURTR = 0.25_dp
      real (dp), parameter :: THIRD = 0.3330_dp
      real (dp), parameter :: HALF = 0.5_dp
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
      real (dp), parameter :: TWO = 2.0_dp
      real (dp), parameter :: HUNDRD = 100.0_dp
!
      integer   :: I4, NN, NP
      real (dp) :: A2, B1, B2, GAM, GAP1, GAP2, S
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  A negative DMIN forces the shift to take that absolute value
!  TTYPE records the type of shift.
!
      IF (DMIN.LE.ZERO) THEN
         TAU = -DMIN
         TTYPE = -1
         RETURN
      END IF
!
      NN = 4*N0 + PP
      IF (N0IN.EQ.N0) THEN
!
!  No eigenvalues deflated.
!
         IF (DMIN.EQ.DN .OR. DMIN.EQ.DN1) THEN
!
            B1 = SQRT(Z(NN-3))*SQRT(Z(NN-5))
            B2 = SQRT(Z(NN-7))*SQRT(Z(NN-9))
            A2 = Z(NN-7) + Z(NN-5)
!
!  Cases 2 and 3.
!
            IF (DMIN.EQ.DN .AND. DMIN1.EQ.DN1) THEN
               GAP2 = DMIN2 - A2 - DMIN2*QURTR
               IF (GAP2.GT.ZERO .AND. GAP2.GT.B2) THEN
                  GAP1 = A2 - DN - (B2 / GAP2)*B2
               ELSE
                  GAP1 = A2 - DN - (B1+B2)
               END IF
               IF (GAP1.GT.ZERO .AND. GAP1.GT.B1) THEN
                  S = MAX(DN-(B1 / GAP1)*B1, HALF*DMIN)
                  TTYPE = -2
               ELSE
                  S = ZERO
                  IF (DN.GT.B1) S = DN - B1
                  IF (A2.GT.(B1+B2)) S = MIN(S, A2-(B1+B2))
                  S = MAX(S, THIRD*DMIN)
                  TTYPE = -3
               END IF
            ELSE
!
!  Case 4.
!
               TTYPE = -4
               S = QURTR*DMIN
               IF (DMIN.EQ.DN) THEN
                  GAM = DN
                  A2 = ZERO
                  IF (Z(NN-5) .GT. Z(NN-7)) RETURN
                  B2 = Z(NN-5) / Z(NN-7)
                  NP = NN - 9
               ELSE
                  NP = NN - 2*PP
                  GAM = DN1
                  IF (Z(NP-4) .GT. Z(NP-2)) RETURN
                  A2 = Z(NP-4) / Z(NP-2)
                  IF (Z(NN-9) .GT. Z(NN-11)) RETURN
                  B2 = Z(NN-9) / Z(NN-11)
                  NP = NN - 13
               END IF
!
!  Approximate contribution to norm squared from I < NN-1.
!
               A2 = A2 + B2
               DO I4 = NP, 4*I0 - 1 + PP, -4
                  IF (B2.EQ.ZERO) GOTO 20
                  B1 = B2
                  IF (Z(I4) .GT. Z(I4-2)) RETURN
                  B2 = B2*(Z(I4) / Z(I4-2))
                  A2 = A2 + B2
                  IF (HUNDRD*MAX(B2,B1).LT.A2 .OR. CNST1.LT.A2) GOTO 20
               END DO
   20          CONTINUE
               A2 = CNST3*A2
!
!  Rayleigh quotient residual bound.
!
               IF (A2.LT.CNST1) S = GAM*(ONE-SQRT(A2)) / (ONE+A2)
            END IF
         ELSE IF (DMIN.EQ.DN2) THEN
!
!  Case 5.
!
            TTYPE = -5
            S = QURTR*DMIN
!
!  Compute contribution to norm squared from I > NN-2.
!
            NP = NN - 2*PP
            B1 = Z(NP-2)
            B2 = Z(NP-6)
            GAM = DN2
            IF (Z(NP-8).GT.B2 .OR. Z(NP-4).GT.B1) RETURN
            A2 = (Z(NP-8) / B2)*(ONE+Z(NP-4) / B1)
!
!  Approximate contribution to norm squared from I < NN-2.
!
            IF (N0-I0.GT.2) THEN
               B2 = Z(NN-13) / Z(NN-15)
               A2 = A2 + B2
               DO I4 = NN - 17, 4*I0 - 1 + PP, -4
                  IF (B2.EQ.ZERO) GOTO 40
                  B1 = B2
                  IF (Z(I4) .GT. Z(I4-2))  RETURN
                  B2 = B2*(Z(I4) / Z(I4-2))
                  A2 = A2 + B2
                  IF (HUNDRD*MAX(B2,B1).LT.A2 .OR. CNST1.LT.A2) GOTO 40
               END DO
   40          CONTINUE
               A2 = CNST3*A2
            END IF
!
            IF (A2.LT.CNST1) S = GAM*(ONE-SQRT(A2)) / (ONE+A2)
         ELSE
!
!  Case 6, no information to guide us.
!
            IF (TTYPE.EQ.-6) THEN
               G = G + THIRD*(ONE-G)
            ELSE IF (TTYPE.EQ.-18) THEN
               G = QURTR*THIRD
            ELSE
               G = QURTR
            END IF
            S = G*DMIN
            TTYPE = -6
         END IF
!
      ELSE IF (N0IN.EQ.(N0+1)) THEN
!
!  One eigenvalue just deflated. Use DMIN1, DN1 for DMIN and DN.
!
         IF (DMIN1.EQ.DN1 .AND. DMIN2.EQ.DN2) THEN
!
!  Cases 7 and 8.
!
            TTYPE = -7
            S = THIRD*DMIN1
            IF (Z(NN-5).GT.Z(NN-7)) RETURN
            B1 = Z(NN-5) / Z(NN-7)
            B2 = B1
            IF (B2.EQ.ZERO) GOTO 60
            DO I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               A2 = B1
               IF (Z(I4).GT.Z(I4-2)) RETURN
               B1 = B1*(Z(I4) / Z(I4-2))
               B2 = B2 + B1
               IF (HUNDRD*MAX(B1, A2).LT.B2)  GOTO 60
            END DO
   60       CONTINUE
            B2 = SQRT(CNST3*B2)
            A2 = DMIN1 / (ONE+B2**2)
            GAP2 = HALF*DMIN2 - A2
            IF (GAP2.GT.ZERO .AND. GAP2.GT.B2*A2) THEN
               S = MAX(S, A2*(ONE-CNST2*A2*(B2 / GAP2)*B2))
            ELSE
               S = MAX(S, A2*(ONE-CNST2*B2))
               TTYPE = -8
            END IF
         ELSE
!
!  Case 9.
!
            S = QURTR*DMIN1
            IF (DMIN1.EQ.DN1) S = HALF*DMIN1
            TTYPE = -9
         END IF
!
      ELSE IF (N0IN.EQ.(N0+2)) THEN
!
!  Two eigenvalues deflated. Use DMIN2, DN2 for DMIN and DN.
!
!  Cases 10 and 11.
!
         IF (DMIN2.EQ.DN2 .AND. TWO*Z(NN-5).LT.Z(NN-7)) THEN
            TTYPE = -10
            S = THIRD*DMIN2
            IF (Z(NN-5).GT.Z(NN-7)) RETURN
            B1 = Z(NN-5) / Z(NN-7)
            B2 = B1
            IF (B2.EQ.ZERO) GOTO 80
            DO I4 = 4*N0 - 9 + PP, 4*I0 - 1 + PP, -4
               IF (Z(I4).GT.Z(I4-2)) RETURN
               B1 = B1*(Z(I4) / Z(I4-2))
               B2 = B2 + B1
               IF (HUNDRD*B1.LT.B2) GOTO 80
            END DO
   80       CONTINUE
            B2 = SQRT(CNST3*B2)
            A2 = DMIN2 / (ONE+B2**2)
            GAP2 = Z(NN-7) + Z(NN-9) -                              &
     &             SQRT(Z(NN-11))*SQRT(Z(NN-9)) - A2
            IF (GAP2.GT.ZERO .AND. GAP2.GT.B2*A2) THEN
               S = MAX(S, A2*(ONE-CNST2*A2*(B2 / GAP2)*B2))
            ELSE
               S = MAX(S, A2*(ONE-CNST2*B2))
            END IF
         ELSE
            S = QURTR*DMIN2
            TTYPE = -11
         END IF
      ELSE IF (N0IN.GT.(N0+2)) THEN
!
!  Case 12, more than two eigenvalues deflated. No information.
!
         S = ZERO
         TTYPE = -12
      END IF
!
      TAU = S
!
      RETURN
      END SUBROUTINE DLASQ4
!
      SUBROUTINE DLASQ5 (I0, N0, Z, PP, TAU, SIGMA, DMIN, DMIN1,        &
     &                   DMIN2, DN, DNM1, DNM2, IEEE, EPS)
!
!=======================================================================
!                                                                      !
!  DLASQ5 computes one dqds transform in ping-pong form, one version   !
!  for IEEE machines another for non IEEE machines.                    !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  I0       [in]      INTEGER                                          !
!           First index.                                               !
!                                                                      !
!  N0       [in]      INTEGER                                          !
!           Last index.                                                !
!                                                                      !
!  Z        [in,out]  REAL array, dimension (4*N)                      !
!           Z holds the qd array. EMIN is stored in Z(4*N0) to avoid   !
!           an extra argument.                                         !
!                                                                      !
!  PP       [in]      INTEGER                                          !
!           PP=0 for ping, PP=1 for pong.                              !
!                                                                      !
!  TAU      [in,out]  REAL                                             !
!           This is the shift.                                         !
!                                                                      !
!  SIGMA    [in]      REAL                                             !
!           This is the accumulated shift up to this step.             !
!                                                                      !
!  DMIN     [out]     REAL                                             !
!           Minimum value of d.                                        !
!                                                                      !
!  DMIN1    [out]     REAL                                             !
!           Minimum value of d, excluding D(N0).                       !
!                                                                      !
!  DMIN2    [out]     REAL                                             !
!           Minimum value of d, excluding D(N0) and D(N0-1).           !
!                                                                      !
!  DN       [out]     REAL                                             !
!           d(N0), the last value of d.                                !
!                                                                      !
!  DNM1     [out]     REAL                                             !
!           d(N0-1).                                                   !
!                                                                      !
!  DNM2     [out]     REAL                                             !
!           d(N0-2).                                                   !
!                                                                      !
!  IEEE     [in]      LOGICAL                                          !
!           Flag for IEEE or non IEEE arithmetic.                      !
!                                                                      !
!  EPS      [in]      REAL                                             !
!           This is the value of epsilon used.                         !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  -- LAPACK computational routine --                                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical,   intent(in)    :: IEEE
!
      integer,   intent(in)    :: I0, N0, PP
!
      real (dp), intent(in)    :: EPS, SIGMA
      real (dp), intent(inout) :: TAU
      real (dp), intent(inout) :: Z(*)
      real (dp), intent(out)   :: DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: HALF = 0.5_dp
!
      integer   :: J4, J4P2
      real (dp) ::  D, EMIN, TEMP, DTHRESH
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF ((N0-I0-1).LE.0) RETURN
!
      DTHRESH = EPS*(SIGMA+TAU)
      IF (TAU.LT.DTHRESH*HALF) TAU = ZERO
!
      IF (TAU.NE.ZERO) THEN
        J4 = 4*I0 + PP - 3
        EMIN = Z(J4+4)
        D = Z(J4) - TAU
        DMIN = D
        DMIN1 = -Z(J4)
!
        IF (IEEE) THEN
!
!  Code for IEEE arithmetic.
!
          IF (PP.EQ.0) THEN
            DO J4 = 4*I0, 4*(N0-3), 4
               Z(J4-2) = D + Z(J4-1)
               TEMP = Z(J4+1) / Z(J4-2)
               D = D*TEMP - TAU
               DMIN = MIN(DMIN, D)
               Z(J4) = Z(J4-1)*TEMP
               EMIN = MIN(Z(J4), EMIN)
            END DO
          ELSE
            DO J4 = 4*I0, 4*(N0-3), 4
               Z(J4-3) = D + Z(J4)
               TEMP = Z(J4+2) / Z(J4-3)
               D = D*TEMP - TAU
               DMIN = MIN(DMIN, D)
               Z(J4-1) = Z(J4)*TEMP
               EMIN = MIN(Z(J4-1), EMIN)
            END DO
          END IF
!
!  Unroll last two steps.
!
          DNM2 = D
          DMIN2 = DMIN
          J4 = 4*(N0-2) - PP
          J4P2 = J4 + 2*PP - 1
          Z(J4-2) = DNM2 + Z(J4P2)
          Z(J4) = Z(J4P2+2)*(Z(J4P2) / Z(J4-2))
          DNM1 = Z(J4P2+2)*(DNM2 / Z(J4-2)) - TAU
          DMIN = MIN(DMIN, DNM1)
!
          DMIN1 = DMIN
          J4 = J4 + 4
          J4P2 = J4 + 2*PP - 1
          Z(J4-2) = DNM1 + Z(J4P2)
          Z(J4) = Z(J4P2+2)*(Z(J4P2) / Z(J4-2))
          DN = Z(J4P2+2)*(DNM1 / Z(J4-2)) - TAU
          DMIN = MIN(DMIN, DN)
!
        ELSE
!
!  Code for non IEEE arithmetic.
!
          IF (PP.EQ.0) THEN
            DO J4 = 4*I0, 4*(N0-3), 4
              Z(J4-2) = D + Z(J4-1)
              IF (D.LT.ZERO) THEN
                RETURN
              ELSE
                Z(J4) = Z(J4+1)*(Z(J4-1) / Z(J4-2))
                D = Z(J4+1)*(D / Z(J4-2)) - TAU
              END IF
              DMIN = MIN(DMIN, D)
              EMIN = MIN(EMIN, Z(J4))
            END DO
          ELSE
            DO J4 = 4*I0, 4*(N0-3), 4
              Z(J4-3) = D + Z(J4)
              IF (D.LT.ZERO) THEN
                 RETURN
              ELSE
                 Z(J4-1) = Z(J4+2)*(Z(J4) / Z(J4-3))
                 D = Z(J4+2)*(D / Z(J4-3)) - TAU
              END IF
              DMIN = MIN(DMIN, D)
              EMIN = MIN(EMIN, Z(J4-1))
            END DO
          END IF
!
!  Unroll last two steps.
!
          DNM2 = D
          DMIN2 = DMIN
          J4 = 4*(N0-2) - PP
          J4P2 = J4 + 2*PP - 1
          Z(J4-2) = DNM2 + Z(J4P2)
          IF (DNM2.LT.ZERO) THEN
            RETURN
          ELSE
            Z(J4) = Z(J4P2+2)*(Z(J4P2) / Z(J4-2))
            DNM1 = Z(J4P2+2)*(DNM2 / Z(J4-2)) - TAU
          END IF
          DMIN = MIN(DMIN, DNM1)
!
          DMIN1 = DMIN
          J4 = J4 + 4
          J4P2 = J4 + 2*PP - 1
          Z(J4-2) = DNM1 + Z(J4P2)
          IF (DNM1.LT.ZERO) THEN
            RETURN
          ELSE
            Z(J4) = Z(J4P2+2)*(Z(J4P2) / Z(J4-2))
            DN = Z(J4P2+2)*(DNM1 / Z(J4-2)) - TAU
          END IF
          DMIN = MIN(DMIN, DN)
!
        END IF

      ELSE
!
!  This is the version that sets d's to zero if they are small enough.
!
        J4 = 4*I0 + PP - 3
        EMIN = Z(J4+4)
        D = Z(J4) - TAU
        DMIN = D
        DMIN1 = -Z(J4)

        IF (IEEE) THEN
!
!  Code for IEEE arithmetic.
!
          IF (PP.EQ.0) THEN
            DO J4 = 4*I0, 4*(N0-3), 4
              Z(J4-2) = D + Z(J4-1)
              TEMP = Z(J4+1) / Z(J4-2)
              D = D*TEMP - TAU
              IF (D.LT.DTHRESH) D = ZERO
              DMIN = MIN(DMIN, D)
              Z(J4) = Z(J4-1)*TEMP
              EMIN = MIN(Z(J4), EMIN)
            END DO
          ELSE
            DO J4 = 4*I0, 4*(N0-3), 4
              Z(J4-3) = D + Z(J4)
              TEMP = Z(J4+2) / Z(J4-3)
              D = D*TEMP - TAU
              IF (D.LT.DTHRESH) D = ZERO
              DMIN = MIN(DMIN, D)
              Z(J4-1) = Z(J4)*TEMP
              EMIN = MIN(Z(J4-1), EMIN)
            END DO
          END IF
!
!  Unroll last two steps.
!
          DNM2 = D
          DMIN2 = DMIN
          J4 = 4*(N0-2) - PP
          J4P2 = J4 + 2*PP - 1
          Z(J4-2) = DNM2 + Z(J4P2)
          Z(J4) = Z(J4P2+2)*(Z(J4P2) / Z(J4-2))
          DNM1 = Z(J4P2+2)*(DNM2 / Z(J4-2)) - TAU
          DMIN = MIN(DMIN, DNM1)
!
          DMIN1 = DMIN
          J4 = J4 + 4
          J4P2 = J4 + 2*PP - 1
          Z(J4-2) = DNM1 + Z(J4P2)
          Z(J4) = Z(J4P2+2)*(Z(J4P2) / Z(J4-2))
          DN = Z(J4P2+2)*(DNM1 / Z(J4-2)) - TAU
          DMIN = MIN(DMIN, DN)
!
        ELSE
!
!  Code for non IEEE arithmetic.
!
          IF (PP.EQ.0) THEN
            DO J4 = 4*I0, 4*(N0-3), 4
              Z(J4-2) = D + Z(J4-1)
              IF (D.LT.ZERO) THEN
                RETURN
              ELSE
                Z(J4) = Z(J4+1)*(Z(J4-1) / Z(J4-2))
                D = Z(J4+1)*(D / Z(J4-2)) - TAU
              END IF
              IF (D.LT.DTHRESH) D = ZERO
              DMIN = MIN(DMIN, D)
              EMIN = MIN(EMIN, Z(J4))
            END DO
          ELSE
            DO J4 = 4*I0, 4*(N0-3), 4
              Z(J4-3) = D + Z(J4)
              IF (D.LT.ZERO) THEN
                RETURN
              ELSE
                Z(J4-1) = Z(J4+2)*(Z(J4) / Z(J4-3))
                D = Z(J4+2)*(D / Z(J4-3)) - TAU
              END IF
              IF (D.LT.DTHRESH) D = ZERO
              DMIN = MIN(DMIN, D)
              EMIN = MIN(EMIN, Z(J4-1))
            END DO
          END IF
!
!  Unroll last two steps.
!
          DNM2 = D
          DMIN2 = DMIN
          J4 = 4*(N0-2) - PP
          J4P2 = J4 + 2*PP - 1
          Z(J4-2) = DNM2 + Z(J4P2)
          IF (DNM2.LT.ZERO) THEN
            RETURN
          ELSE
            Z(J4) = Z(J4P2+2)*(Z(J4P2) / Z(J4-2))
            DNM1 = Z(J4P2+2)*(DNM2 / Z(J4-2)) - TAU
          END IF
          DMIN = MIN(DMIN, DNM1)
!
          DMIN1 = DMIN
          J4 = J4 + 4
          J4P2 = J4 + 2*PP - 1
          Z(J4-2) = DNM1 + Z(J4P2)
          IF (DNM1.LT.ZERO) THEN
            RETURN
          ELSE
            Z(J4) = Z(J4P2+2)*(Z(J4P2) / Z(J4-2))
            DN = Z(J4P2+2)*(DNM1 / Z(J4-2)) - TAU
          END IF
          DMIN = MIN(DMIN, DN)
!
        END IF
      END IF
!
      Z(J4+2) = DN
      Z(4*N0-PP) = EMIN
!
      RETURN
      END SUBROUTINE DLASQ5
!
      SUBROUTINE DLASQ6 (I0, N0, Z, PP, DMIN, DMIN1, DMIN2, DN,         &
     &                   DNM1, DNM2)
!
!=======================================================================
!                                                                      !
!  DLASQ6 computes one dqd (shift equal to zero) transform in          !
!  ping-pong form, with protection against underflow and overflow.     !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  I0       [in]      INTEGER                                          !
!           First index.                                               !
!                                                                      !
!  N0       [in]      INTEGER                                          !
!           Last index.                                                !
!                                                                      !
!  Z        [in,out]  REAL array, dimension (4*N)                      !
!           Z holds the qd array. EMIN is stored in Z(4*N0) to avoid   !
!           an extra argument.                                         !
!                                                                      !
!  PP       [in]      INTEGER                                          !
!           PP=0 for ping, PP=1 for pong.                              !
!                                                                      !
!  DMIN     [out]     REAL                                             !
!           Minimum value of d.                                        !
!                                                                      !
!  DMIN1    [out]     REAL                                             !
!           Minimum value of d, excluding D(N0).                       !
!                                                                      !
!  DMIN2    [out]     REAL                                             !
!           Minimum value of d, excluding D(N0) and D(N0-1).           !
!                                                                      !
!  DN       [out]     REAL                                             !
!           d(N0), the last value of d.                                !
!                                                                      !
!  DNM1     [out]     REAL                                             !
!           d(N0-1).                                                   !
!                                                                      !
!  DNM2     [out]     REAL                                             !
!           d(N0-2).                                                   !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  -- LAPACK computational routine --                                  !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: I0, N0, PP
!
      real (dp), intent(inout) :: Z(*)
      real (dp), intent(out)   :: DMIN, DMIN1, DMIN2, DN, DNM1, DNM2
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
!
      integer   :: J4, J4P2
      real (dp) :: D, EMIN, SAFMIN, TEMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF ((N0-I0-1).LE.0) RETURN
!
      SAFMIN = DLAMCH('Safe minimum')
      J4 = 4*I0 + PP - 3
      EMIN = Z(J4+4)
      D = Z(J4)
      DMIN = D
!
      IF (PP.EQ.0) THEN
         DO J4 = 4*I0, 4*(N0-3), 4
            Z(J4-2) = D + Z(J4-1)
            IF (Z(J4-2).EQ.ZERO) THEN
               Z(J4) = ZERO
               D = Z(J4+1)
               DMIN = D
               EMIN = ZERO
            ELSE IF (SAFMIN*Z(J4+1).LT.Z(J4-2) .AND.                    &
     &               SAFMIN*Z(J4-2).LT.Z(J4+1)) THEN
               TEMP = Z(J4+1) / Z(J4-2)
               Z(J4) = Z(J4-1)*TEMP
               D = D*TEMP
            ELSE
               Z(J4) = Z(J4+1)*(Z(J4-1) / Z(J4-2))
               D = Z(J4+1)*(D / Z(J4-2))
            END IF
            DMIN = MIN(DMIN, D)
            EMIN = MIN(EMIN, Z(J4))
         END DO
      ELSE
         DO J4 = 4*I0, 4*(N0-3), 4
            Z(J4-3) = D + Z(J4)
            IF (Z(J4-3).EQ.ZERO) THEN
               Z(J4-1) = ZERO
               D = Z(J4+2)
               DMIN = D
               EMIN = ZERO
            ELSE IF (SAFMIN*Z(J4+2).LT.Z(J4-3) .AND.                    &
     &               SAFMIN*Z(J4-3).LT.Z(J4+2)) THEN
               TEMP = Z(J4+2) / Z(J4-3)
               Z(J4-1) = Z(J4)*TEMP
               D = D*TEMP
            ELSE
               Z(J4-1) = Z(J4+2)*(Z(J4) / Z(J4-3))
               D = Z(J4+2)*(D / Z(J4-3))
            END IF
            DMIN = MIN(DMIN, D)
            EMIN = MIN(EMIN, Z(J4-1))
         END DO
      END IF
!
!  Unroll last two steps.
!
      DNM2 = D
      DMIN2 = DMIN
      J4 = 4*(N0-2) - PP
      J4P2 = J4 + 2*PP - 1
      Z(J4-2) = DNM2 + Z(J4P2)
      IF (Z(J4-2).EQ.ZERO) THEN
         Z(J4) = ZERO
         DNM1 = Z(J4P2+2)
         DMIN = DNM1
         EMIN = ZERO
      ELSE IF (SAFMIN*Z(J4P2+2).LT.Z(J4-2) .AND.                        &
     &         SAFMIN*Z(J4-2).LT.Z(J4P2+2)) THEN
         TEMP = Z(J4P2+2) / Z(J4-2)
         Z(J4) = Z(J4P2)*TEMP
         DNM1 = DNM2*TEMP
      ELSE
         Z(J4) = Z(J4P2+2)*(Z(J4P2) / Z(J4-2))
         DNM1 = Z(J4P2+2)*(DNM2 / Z(J4-2))
      END IF
      DMIN = MIN(DMIN, DNM1)
!
      DMIN1 = DMIN
      J4 = J4 + 4
      J4P2 = J4 + 2*PP - 1
      Z(J4-2) = DNM1 + Z(J4P2)
      IF (Z(J4-2).EQ.ZERO) THEN
         Z(J4) = ZERO
         DN = Z(J4P2+2)
         DMIN = DN
         EMIN = ZERO
      ELSE IF (SAFMIN*Z(J4P2+2).LT.Z(J4-2) .AND.                        &
     &         SAFMIN*Z(J4-2).LT.Z(J4P2+2)) THEN
         TEMP = Z(J4P2+2) / Z(J4-2)
         Z(J4) = Z(J4P2)*TEMP
         DN = DNM1*TEMP
      ELSE
         Z(J4) = Z(J4P2+2)*(Z(J4P2) / Z(J4-2))
         DN = Z(J4P2+2)*(DNM1 / Z(J4-2))
      END IF
      DMIN = MIN(DMIN, DN)
!
      Z(J4+2) = DN
      Z(4*N0-PP) = EMIN
!
      RETURN
      END SUBROUTINE DLASQ6
!
      SUBROUTINE DLASR (SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA)
!
!=======================================================================
!                                                                      !
!  DLASR performs the transformation:                                  !
!                                                                      !
!     A := P*A,   when SIDE = 'L' or 'l'  ( Left-hand side)            !
!                                                                      !
!     A := A*P',  when SIDE = 'R' or 'r'  (Right-hand side)            !
!                                                                      !
!  where A is an m by n real matrix and P is an orthogonal matrix,     !
!  consisting of a sequence of plane rotations determined by the       !
!  parameters PIVOT and DIRECT as follows (z = m when SIDE = 'L'       !
!  or 'l' and z = n when SIDE = 'R' or 'r'):                           !
!                                                                      !
!  When  DIRECT = 'F' or 'f'  (Forward sequence) then                  !
!                                                                      !
!     P = P(z - 1)*...*P(2)*P(1),                                      !
!                                                                      !
!  and when DIRECT = 'B' or 'b'  (Backward sequence) then              !
!                                                                      !
!     P = P(1)*P(2)*...*P(z - 1),                                      !
!                                                                      !
!  where  P(k) is a plane rotation matrix for the following planes:    !
!                                                                      !
!     when  PIVOT = 'V' or 'v'  (Variable pivot),                      !
!        the plane (k, k + 1)                                          !
!                                                                      !
!     when  PIVOT = 'T' or 't'  (Top pivot),                           !
!        the plane (1, k + 1)                                          !
!                                                                      !
!     when  PIVOT = 'B' or 'b'  (Bottom pivot),                        !
!        the plane (k, z)                                              !
!                                                                      !
!  c(k) and s(k)  must contain the  cosine and sine that define        !
!  the matrix  P(k).  The two by two plane rotation part of the        !
!  matrix P(k), R(k), is assumed to be of the form                     !
!                                                                      !
!     R(k) = ( c(k)  s(k)).                                            !
!              (-s(k)  c(k))                                           !
!                                                                      !
!  This version vectorises across rows of the array A when             !
!  SIDE = 'L'.                                                         !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  SIDE    (input) CHARACTER*1                                         !
!          Specifies whether the plane rotation matrix P is applied    !
!          to A on the left or the right.                              !
!          = 'L':  Left, compute A := P*A                              !
!          = 'R':  Right, compute A:= A*P'                             !
!                                                                      !
!  DIRECT  (input) CHARACTER*1                                         !
!          Specifies whether P is a forward or backward sequence of    !
!          plane rotations.                                            !
!          = 'F':  Forward, P = P(z - 1)*...*P(2)*P(1)                 !
!          = 'B':  Backward, P = P(1)*P(2)*...*P(z - 1)                !
!                                                                      !
!  PIVOT   (input) CHARACTER*1                                         !
!          Specifies the plane for which P(k) is a plane rotation      !
!          matrix.                                                     !
!          = 'V':  Variable pivot, the plane (k,k+1)                   !
!          = 'T':  Top pivot, the plane (1,k+1)                        !
!          = 'B':  Bottom pivot, the plane (k,z)                       !
!                                                                      !
!  M       (input) INTEGER                                             !
!          The number of rows of the matrix A.  If m <= 1, an          !
!          immediate return is effected.                               !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The number of columns of the matrix A.  If n <= 1, an       !
!          immediate return is effected.                               !
!                                                                      !
!  C, S    (input) DOUBLE PRECISION arrays, dimension                  !
!                  (M-1) if SIDE = 'L'                                 !
!                  (N-1) if SIDE = 'R'                                 !
!          c(k) and s(k) contain the cosine and sine that define the   !
!          matrix P(k).  The two by two plane rotation part of the     !
!          matrix P(k), R(k), is assumed to be of the form             !
!          R(k) = ( c(k)  s(k)).                                       !
!                 (-s(k)  c(k))                                        !
!                                                                      !
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)    !
!          The m by n matrix A.  On exit, A is overwritten by P*A if   !
!          SIDE = 'R' or by A*P' if SIDE = 'L'.                        !
!                                                                      !
!  LDA     (input) INTEGER                                             !
!          The leading dimension of the array A.  LDA >= max(1,M).     !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: LDA, M, N
!
      real (dp), intent(in)    :: C(*), S(*)
      real (dp), intent(inout) :: A(LDA, *)
!
      character (len=1), intent(in) :: DIRECT, PIVOT, SIDE
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
      real (dp), parameter :: ONE = 1.0_dp
!
      integer :: I, INFO, J
!
      real (dp) :: CTEMP, STEMP, TEMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input parameters.
!
      INFO = 0
      IF (.NOT. (LSAME(SIDE, 'L') .OR. LSAME(SIDE, 'R'))) THEN
        INFO = 1
      ELSE IF (.NOT.(LSAME(PIVOT, 'V') .OR. LSAME(PIVOT, 'T') .OR.      &
     &         LSAME(PIVOT, 'B'))) THEN
        INFO = 2
      ELSE IF (.NOT.(LSAME(DIRECT, 'F') .OR. LSAME(DIRECT, 'B'))) THEN
        INFO = 3
      ELSE IF (M .LT. 0) THEN
        INFO = 4
      ELSE IF (N .LT. 0) THEN
        INFO = 5
      ELSE IF (LDA .LT. MAX(1, M)) THEN
        INFO = 9
      END IF
      IF (INFO .NE. 0) THEN
        CALL XERBLA ('DLASR ', INFO)
        RETURN
      END IF
!
!  Quick return if possible.
!
      IF ((M.EQ.0) .OR. (N.EQ.0)) RETURN
!
      IF (LSAME(SIDE, 'L')) THEN
!
!  Form P * A.
!
        IF (LSAME(PIVOT, 'V')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 1, M - 1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J+1,I)
                   A(J+1,I) = CTEMP*TEMP - STEMP*A(J,I)
                   A(J,  I) = STEMP*TEMP + CTEMP*A(J,I)
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = M - 1, 1, -1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J+1,I)
                  A(J+1,I) = CTEMP*TEMP - STEMP*A(J,I)
                  A(J,  I) = STEMP*TEMP + CTEMP*A(J,I)
                END DO
              END IF
            END DO
          END IF
        ELSE IF (LSAME(PIVOT, 'T')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 2, M
              CTEMP = C(J-1)
              STEMP = S(J-1)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J,I)
                  A(J,I) = CTEMP*TEMP - STEMP*A(1,I)
                  A(1,I) = STEMP*TEMP + CTEMP*A(1,I)
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = M, 2, -1
              CTEMP = C(J-1)
              STEMP = S(J-1)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J,I)
                  A(J,I) = CTEMP*TEMP - STEMP*A(1,I)
                  A(1,I) = STEMP*TEMP + CTEMP*A(1,I)
                END DO
              END IF
            END DO
          END IF
        ELSE IF (LSAME(PIVOT, 'B')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 1, M - 1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J,I)
                  A(J,I) = STEMP*A(M,I) + CTEMP*TEMP
                  A(M,I) = CTEMP*A(M,I) - STEMP*TEMP
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = M - 1, 1, -1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, N
                  TEMP = A(J,I)
                  A(J,I) = STEMP*A(M,I) + CTEMP*TEMP
                  A(M,I) = CTEMP*A(M,I) - STEMP*TEMP
                END DO
              END IF
            END DO
          END IF
        END IF
      ELSE IF (LSAME(SIDE, 'R')) THEN
!
!  Form A * P'.
!
        IF (LSAME(PIVOT, 'V')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 1, N - 1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J+1)
                  A(I,J+1) = CTEMP*TEMP - STEMP*A(I,J)
                  A(I,J ) = STEMP*TEMP + CTEMP*A(I,J)
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = N - 1, 1, -1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J+1)
                  A(I,J+1) = CTEMP*TEMP - STEMP*A(I,J)
                  A(I,J ) = STEMP*TEMP + CTEMP*A(I,J)
                END DO
              END IF
            END DO
          END IF
        ELSE IF (LSAME(PIVOT, 'T')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 2, N
              CTEMP = C(J-1)
              STEMP = S(J-1)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J)
                  A(I,J) = CTEMP*TEMP - STEMP*A(I,1)
                  A(I,1) = STEMP*TEMP + CTEMP*A(I,1)
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = N, 2, -1
              CTEMP = C(J-1)
              STEMP = S(J-1)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J)
                  A(I,J) = CTEMP*TEMP - STEMP*A(I,1)
                  A(I,1) = STEMP*TEMP + CTEMP*A(I,1)
                END DO
              END IF
            END DO
          END IF
        ELSE IF (LSAME(PIVOT, 'B')) THEN
          IF (LSAME(DIRECT, 'F')) THEN
            DO J = 1, N - 1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J)
                  A(I,J) = STEMP*A(I,N) + CTEMP*TEMP
                  A(I,N) = CTEMP*A(I,N) - STEMP*TEMP
                END DO
              END IF
            END DO
          ELSE IF (LSAME(DIRECT, 'B')) THEN
            DO J = N - 1, 1, -1
              CTEMP = C(J)
              STEMP = S(J)
              IF ((CTEMP.NE.ONE) .OR. (STEMP.NE.ZERO)) THEN
                DO I = 1, M
                  TEMP = A(I,J)
                  A(I,J) = STEMP*A(I,N) + CTEMP*TEMP
                  A(I,N) = CTEMP*A(I,N) - STEMP*TEMP
                END DO
              END IF
            END DO
          END IF
        END IF
      END IF
!
      RETURN
      END SUBROUTINE DLASR
!
      SUBROUTINE DLASRT (ID, N, D, INFO)
!
!=======================================================================
!                                                                      !
!  Sort the numbers in D in increasing order (if ID = 'I') or          !
!  in decreasing order (if ID = 'D').                                  !
!                                                                      !
!  Use Quick Sort, reverting to Insertion sort on arrays of            !
!  size <= 20. Dimension of STACK limits N to about 2**32.             !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  ID      (input) CHARACTER*1                                         !
!          = 'I': sort D in increasing order;                          !
!          = 'D': sort D in decreasing order.                          !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The length of the array D.                                  !
!                                                                      !
!  D       (input/output) DOUBLE PRECISION array, dimension (N)        !
!          On entry, the array to be sorted.                           !
!          On exit, D has been sorted into increasing order            !
!          (D(1) <= ... <= D(N)) or into decreasing order              !
!          (D(1) >= ... >= D(N)), depending on ID.                     !
!                                                                      !
!  INFO    (output) INTEGER                                            !
!          = 0:  successful exit                                       !
!          < 0:  if INFO = -i, the i-th argument had an illegal value  !
!                                                                      !
!  -- LAPACK routine (version 2.0) --                                  !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     September 30, 1994                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: N
      integer,   intent(out)   :: INFO
!
      real (dp), intent(inout) :: D(*)
!
      character (len=1), intent(in) :: ID
!
!  Local variable declarations.
!
      integer, parameter :: SELECT = 20
!
      integer :: DIR, ENDD, I, J, START, STKPNT
      integer :: STACK(2, 32)
!
      real (dp) :: D1, D2, D3, DMNMX, TMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
!  Test the input paramters.
!
      INFO = 0
      DIR = -1
      IF (LSAME(ID, 'D')) THEN
        DIR = 0
      ELSE IF (LSAME(ID, 'I')) THEN
        DIR = 1
      END IF
!
      IF (DIR .EQ. -1) THEN
        INFO = -1
      ELSE IF (N .LT. 0) THEN
        INFO = -2
      END IF
      IF (INFO .NE. 0) THEN
        CALL XERBLA ('DLASRT', -INFO)
        RETURN
      END IF
!
!  Quick return if possible
!
      IF (N .LE. 1) RETURN
!
      STKPNT = 1
      STACK(1,1) = 1
      STACK(2,1) = N
!
      SORT_LOOP : DO WHILE (STKPNT .GT. 0)
        START = STACK(1,STKPNT)
        ENDD  = STACK(2,STKPNT)
        STKPNT = STKPNT - 1
        IF ((ENDD-START.LE.SELECT) .AND. (ENDD-START.GT.0)) THEN
!
!  Do Insertion sort on D(START:ENDD).
!
          IF (DIR .EQ. 0) THEN
!
!    Sort into decreasing order.
!
            OUTER1 : DO I = START + 1, ENDD
              INNER1 : DO J = I, START + 1, -1
                IF (D(J) .GT. D(J-1)) THEN
                  DMNMX = D(J)
                  D(J ) = D(J-1)
                  D(J-1) = DMNMX
                ELSE
                  EXIT INNER1
                END IF
              END DO INNER1
            END DO OUTER1
!
          ELSE
!
!    Sort into increasing order.
!
            OUTER2 : DO I = START + 1, ENDD
              INNER2 : DO J = I, START + 1, -1
                IF (D(J) .LT. D(J-1)) THEN
                  DMNMX = D(J)
                  D(J ) = D(J-1)
                  D(J-1) = DMNMX
                ELSE
                  EXIT INNER2
                END IF
              END DO INNER2
            END DO OUTER2
!
          END IF
!
        ELSE IF (ENDD-START .GT. SELECT) THEN
!
!  Partition D(START:ENDD) and stack parts, largest one first.
!
!    Choose partition entry as median of 3
!
          D1 = D(START)
          D2 = D(ENDD)
          I = (START+ENDD) / 2
          D3 = D(I)
          IF (D1 .LT. D2) THEN
            IF (D3 .LT. D1) THEN
              DMNMX = D1
            ELSE IF (D3 .LT. D2) THEN
              DMNMX = D3
            ELSE
              DMNMX = D2
            END IF
          ELSE
            IF (D3 .LT. D2) THEN
              DMNMX = D2
            ELSE IF (D3 .LT. D1) THEN
              DMNMX = D3
            ELSE
              DMNMX = D1
            END IF
          END IF
!
          IF (DIR .EQ. 0) THEN
!
!    Sort into decreasing order.
!
            I = START - 1
            J = ENDD + 1
#ifdef REMOVE_LAPACK_GOTOS
            DecLOOP_3 : DO                           ! forever LOOP 3
              DecLOOP_1 : DO                         ! forever LOOP 1
                J = J - 1
                IF (D(J) .GE. DMNMX) EXIT DecLOOP_1  ! terminate LOOP 1
              END DO DecLOOP_1
!
              DecLOOP_2 : DO                         ! forever LOOP 2
                I = I + 1
                IF (D(I) .LE. DMNMX) EXIT DecLOOP_2  ! terminate LOOP 2
              END DO DecLOOP_2
!
              IF (I .LT. J) THEN
                TMP = D(I)
                D(I) = D(J)
                D(J) = TMP
                CYCLE DecLOOP_3                      ! iterate LOOP 3
              ELSE
                EXIT DecLOOP_3                       ! terminate LOOP 3
              END IF
            END DO DecLOOP_3
#else
   60       CONTINUE
   70       CONTINUE
            J = J - 1
            IF (D(J) .LT. DMNMX) GOTO 70
   80       CONTINUE
            I = I + 1
            IF (D(I) .GT. DMNMX) GOTO 80
            IF (I .LT. J) THEN
              TMP = D(I)
              D(I) = D(J)
              D(J) = TMP
              GOTO 60
            END IF
#endif
            IF (J-START .GT. ENDD-J-1) THEN
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = START
              STACK(2, STKPNT) = J
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = J + 1
              STACK(2, STKPNT) = ENDD
            ELSE
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = J + 1
              STACK(2, STKPNT) = ENDD
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = START
              STACK(2, STKPNT) = J
            END IF
          ELSE
!
!    Sort into increasing order.
!
            I = START - 1
            J = ENDD + 1
#ifdef REMOVE_LAPACK_GOTOS
            IncLOOP_3 : DO                           ! forever LOOP 3
              IncLOOP_1 : DO                         ! forever LOOP 1
                J = J - 1
                IF (D(J) .LE. DMNMX) EXIT IncLOOP_1  ! terminate LOOP 1
              END DO IncLOOP_1
!
              IncLOOP_2 : DO                         ! forever LOOP 2
                I = I + 1
                IF (D(I) .GE. DMNMX) EXIT IncLOOP_2  ! terminate LOOP 2
              END DO IncLOOP_2
!
              IF (I .LT. J) THEN
                TMP = D(I)
                D(I) = D(J)
                D(J) = TMP
                CYCLE IncLOOP_3                      ! iterate LOOP 3
              ELSE
                EXIT IncLOOP_3                       ! terminate LOOP 3
              END IF
            END DO IncLOOP_3
#else
   90       CONTINUE
  100       CONTINUE
            J = J - 1
            IF (D(J) .GT. DMNMX) GOTO 100
  110       CONTINUE
            I = I + 1
            IF (D(I) .LT. DMNMX) GOTO 110
            IF (I .LT. J) THEN
              TMP = D(I)
              D(I) = D(J)
              D(J) = TMP
              GOTO 90
            END IF
#endif
            IF (J-START .GT. ENDD-J-1) THEN
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = START
              STACK(2, STKPNT) = J
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = J + 1
              STACK(2, STKPNT) = ENDD
            ELSE
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = J + 1
              STACK(2, STKPNT) = ENDD
              STKPNT = STKPNT + 1
              STACK(1, STKPNT) = START
              STACK(2, STKPNT) = J
            END IF
          END IF
        END IF
      END DO SORT_LOOP
!
      RETURN
      END SUBROUTINE DLASRT
!
      SUBROUTINE DLASSQ (N, X, INCX, SCALE, SUMSQ)
!
!=======================================================================
!                                                                      !
!  DLASSQ  returns the values  scl  and  smsq  such that               !
!                                                                      !
!     (scl**2)*smsq = x(1)**2 +...+ x(n)**2 + (scale**2)*sumsq         !
!                                                                      !
!  where  x(i) = X(1 + (i - 1)*INCX). The value of  sumsq  is          !
!  assumed to be non-negative and  scl  returns the value              !
!                                                                      !
!     scl = max(scale, abs(x(i))).                                     !
!                                                                      !
!  scale and sumsq must be supplied in SCALE and SUMSQ and             !
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.       !
!                                                                      !
!  The routine makes only one pass through the vector x.               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N       (input) INTEGER                                             !
!          The number of elements to be used from the vector X.        !
!                                                                      !
!  X       (input) DOUBLE PRECISION                                    !
!          The vector for which a scaled sum of squares is computed.   !
!             x(i)  = X(1 + (i - 1)*INCX), 1 <= i <= n.                !
!                                                                      !
!  INCX    (input) INTEGER                                             !
!          The increment between successive values of the vector X.    !
!          INCX > 0.                                                   !
!                                                                      !
!  SCALE   (input/output) DOUBLE PRECISION                             !
!          On entry, the value  scale  in the equation above.          !
!          On exit, SCALE is overwritten with  scl, the scaling factor !
!          for the sum of squares.                                     !
!                                                                      !
!  SUMSQ   (input/output) DOUBLE PRECISION                             !
!          On entry, the value  sumsq  in the equation above.          !
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of !
!          squares from which  scl  has been factored out.             !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer,   intent(in)    :: INCX, N
!
      real (dp), intent(in)    :: X(*)
      real (dp), intent(inout) :: SCALE, SUMSQ
!
!  Local variable declarations.
!
      real (dp), parameter :: ZERO = 0.0_dp
!
      integer   :: IX
!
      real (dp) :: ABSXI
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (N .GT. 0) THEN
        DO IX = 1, 1 + (N-1)*INCX, INCX
          IF (X(IX) .NE. ZERO) THEN
            ABSXI = ABS(X(IX))
            IF (SCALE .LT. ABSXI) THEN
              SUMSQ = 1.0_dp + SUMSQ*(SCALE / ABSXI)**2
              SCALE = ABSXI
            ELSE
              SUMSQ = SUMSQ + (ABSXI / SCALE)**2
            END IF
          END IF
        END DO
      END IF
!
      RETURN
      END SUBROUTINE DLASSQ
!
      SUBROUTINE DSCAL (N, DA, DX, INCX)
!
!=======================================================================
!                                                                      !
!  DSCAL scales a vector by a constant uses unrolled loops for         !
!  increment equal to 1.                                               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N       [in]     INTEGER                                            !
!          number of elements in input vector(s)                       !
!                                                                      !
!  DA      [in]     REAL                                               !
!          On entry, DA specifies the scalar alpha.                    !
!                                                                      !
!  DX      [in,out] REAL array, dimension (1 + (N - 1)*abs(INCX))      !
!                                                                      !
!  INCX    [in]     INTEGER                                            !
!          storage spacing between elements of DX                      !
!                                                                      !
!=======================================================================
!
!  Imported variable definitions.
!
      integer,   intent(in)    :: INCX, N
!
      real (dp), intent(in)    :: DA
      real (dp), intent(inout) :: DX(*)
!
!  Local variable declarations.
!
      integer :: I, M, MP1, NINCX
!
      real (dp), parameter :: ONE = 1.0_dp
!
      IF (N.LE.0 .OR. INCX.LE.0 .OR. DA.EQ.ONE) RETURN
      IF (INCX.EQ.1) THEN
!
!  Code for increment equal to 1.
!
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
!
!  Code for increment not equal to 1.
!
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      END IF
!
      RETURN
      END SUBROUTINE DSCAL
!
      SUBROUTINE DSWAP (N, DX, INCX, DY, INCY)
!
!=======================================================================
!                                                                      !
!  DSWAP interchanges two vectors. It uses unrolled loops for          !
!  increments equal one.                                               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  N       (input) INTEGER                                             !
!          number of elements in input vector(s)                       !
!                                                                      !
!  DX      (input/output) DOUBLE PRECISION                             !
!          array, dimension (1 + (N - 1)*abs(INCX))                    !
!                                                                      !
!  INCX    (input) INTEGER                                             !
!          storage spacing between elements of DX                      !
!                                                                      !
!  DY      (input/output) DOUBLE PRECISION                             !
!          array, dimension (1 + (N - 1)*abs(INCY))                    !
!                                                                      !
!  INCY    (input) INTEGER                                             !
!          storage spacing between elements of DY                      !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     October 31, 1992                                                 !
!                                                                      !
!  Based on Jack Dongarra, LINPACK, 3/11/78.                           !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: INCX, INCY, N
!
      real (dp), intent(inout) :: DX(*), DY(*)
!
!  Local variable declarations.
!
      integer   :: I, IX, IY, M, MP1
!
      real (dp) :: DTEMP
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      IF (N .le. 0) RETURN
!
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
!
!  Code for both increments equal to 1.
!
        M = MOD(N,3)
        IF (M .NE. 0) THEN
          DO I = 1, M
            DTEMP = DX(I)
            DX(I) = DY(I)
            DY(I) = DTEMP
          END DO
          IF (N .LT. 3) RETURN
        END IF
        MP1 = M + 1
        DO I = MP1, N, 3
          DTEMP = DX(I)
          DX(I) = DY(I)
          DY(I) = DTEMP
          DTEMP = DX(I+1)
          DX(I+1) = DY(I+1)
          DY(I+1) = DTEMP
          DTEMP = DX(I+2)
          DX(I+2) = DY(I+2)
          DY(I+2) = DTEMP
        END DO
      ELSE
!
!  Code for unequal increments or equal increments not equal to 1.
!
        IX = 1
        IY = 1
        IF (INCX .LT. 0) IX = (-N+1)*INCX + 1
        IF (INCY .LT. 0) IY = (-N+1)*INCY + 1
        DO I = 1, N
          DTEMP = DX(IX)
          DX(IX) = DY(IY)
          DY(IY) = DTEMP
          IX = IX + INCX
          IY = IY + INCY
        END DO
      END IF
!
      RETURN
      END SUBROUTINE DSWAP
!
      FUNCTION IEEECK (ISPEC, ZERO, ONE) RESULT (ICHECK)
!
!=======================================================================
!                                                                      !
!  IEEECK is called from the ILAENV to verify that Infinity and        !
!  possibly NaN arithmetic is safe (i.e. will not trap).               !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  ISPEC    [in]      INTEGER                                          !
!           Specifies whether to test just for infinity arithmetic     !
!           or whether to test for infinity and NaN arithmetic.        !
!           = 0: Verify infinity arithmetic only.                      !
!           = 1: Verify infinity and NaN arithmetic                    !
!                                                                      !
!  ZERO     [in]      REAL                                             !
!           Must contain the value 0.0                                 !
!           This is passed to prevent the compiler from optimizing     !
!           away this code.                                            !
!                                                                      !
!  ONE      [in]      REAL                                             !
!           Must contain the value 1.0                                 !
!           This is passed to prevent the compiler from optimizing     !
!           away this code.                                            !
!                                                                      !
!  RETURN VALUE:      INTEGER                                          !
!           = 0:  Arithmetic failed to produce the correct answers     !
!           = 1:  Arithmetic produced the correct answers              !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_SUPPORT_INF,        &
     &                                         IEEE_SUPPORT_NAN
!
!  Imported variable declarations.
!
      integer,   intent(in) :: ISPEC
      REAL (dp), intent(in) :: ONE, ZERO
!
!  Local variable declarations.
!
      integer :: ICHECK
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      ICHECK = 1
!
      IF (.NOT.IEEE_SUPPORT_INF(ONE)) THEN
         ICHECK = 0
         RETURN
      END IF
!
!  Return if we were only asked to check infinity arithmetic
!
      IF (ISPEC.EQ.0) RETURN
!
      IF (.NOT.IEEE_SUPPORT_NAN(ONE)) THEN
         ICHECK = 0
         RETURN
      END IF
!
      RETURN
      END FUNCTION IEEECK
!
      FUNCTION ILAENV (ISPEC, NAME, OPTS, N1, N2, N3, N4) RESULT (LAENV)
!
!=======================================================================
!                                                                      !
!  ILAENV is called from the LAPACK routines to choose problem-        !
!  dependent parameters for the local environment. See ISPEC for       !
!  a description of the parameters.                                    !
!                                                                      !
!  ILAENV returns an INTEGER                                           !
!                                                                      !
!  If ILAENV >= 0: ILAENV returns the value of the parameter specified !
!                  by ISPEC.                                           !
!  If ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal    !
!                  value.                                              !
!                                                                      !
!  This version provides a set of parameters which should give good,   !
!  but not optimal, performance on many of the currently available     !
!  computers.  Users are encouraged to modify this subroutine to set   !
!  the tuning parameters for their particular machine using the option !
!  and problem size information in the arguments.                      !
!                                                                      !
!  This routine will not function correctly if it is converted to all  !
!  lower case.  Converting it to all upper case is allowed.            !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  ISPEC    [in]      INTEGER                                          !
!           Specifies the parameter to be returned as the value of     !
!           ILAENV.                                                    !
!           = 1: the optimal blocksize; if this value is 1, an         !
!                unblocked algorithm will give the best performance.   !
!           = 2: the minimum block size for which the block routine    !
!                should be used; if the usable block size is less than !
!                this value, an unblocked routine should be used.      !
!           = 3: the crossover point (in a block routine, for N less   !
!                than this value, an unblocked routine should be used) !
!           = 4: the number of shifts, used in the nonsymmetric        !
!                eigenvalue routines (DEPRECATED)                      !
!           = 5: the minimum column dimension for blocking to be used; !
!                rectangular blocks must have dimension at least k by  !
!                k by m, where k is given by ILAENV(2,...) and m by    !
!                ILAENV(5,...)                                         !
!           = 6: the crossover point for the SVD (when reducing an M   !
!                by N matrix to bidiagonal form, if MAX(M,N)/MIN(M,N)  !
!                exceeds this value, a QR factorization is used first  !
!                to reduce the matrix to a triangular form.)           !
!           = 7: the number of processors                              !
!           = 8: the crossover point for the multishift QR method      !
!                for nonsymmetric eigenvalue problems (DEPRECATED)     !
!           = 9: maximum size of the subproblems at the bottom of the  !
!                computation tree in the divide-and-conquer algorithm  !
!                (used by xGELSD and xGESDD)                           !
!           =10: IEEE infinity and NaN arithmetic can be trusted not   !
!                to trap,                                              !
!           =11: infinity arithmetic can be trusted not to trap        !
!           12 <= ISPEC <= 17:                                         !
!                xHSEQR or related subroutines,                        !
!                see IPARMQ for detailed explanation                   !
!                                                                      !
!  NAME     [in]      CHARACTER*(*)                                    !
!           The name of the calling subroutine, in either upper case   !
!           or lower case.                                             !
!                                                                      !
!  OPTS     [in]      CHARACTER*(*)                                    !
!           The character options to the subroutine NAME, concatenated !
!           into a single character string.  For example, UPLO = 'U',  !
!           TRANS = 'T', and DIAG = 'N' for a triangular routine would !
!           be specified as OPTS = 'UTN'.                              !
!                                                                      !
!  N        [in]      INTEGER                                          !
!                                                                      !
!  N2       [in]      INTEGER                                          !
!                                                                      !
!  N3       [in]      INTEGER                                          !
!                                                                      !
!  N4       [in]      INTEGER                                          !
!           Problem dimensions for the subroutine NAME; these may not  !
!           all be required.                                           !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!   The following conventions have been used when calling ILAENV from  !
!   the LAPACK routines:                                               !
!   1)  OPTS is a concatenation of all of the character options to     !
!       subroutine NAME, in the same order that they appear in the     !
!       argument list for NAME, even if they are not used in           !
!       determining the value of the parameter specified by ISPEC.     !
!   2)  The problem dimensions N1, N2, N3, N4 are specified in the     !
!       order that they appear in the argument list for NAME.  N1 is   !
!       used first, N2 second, and so on, and unused problem           !
!       dimensions are passed a value of -1.                           !
!   3)  The parameter value returned by ILAENV is checked for validity !
!       in the calling subroutine.  For example, ILAENV is used to     !
!       retrieve the optimal blocksize for STRTRI as follows:          !
!                                                                      !
!       NB = ILAENV(1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1)          !
!       IF (NB.LE.1) NB = MAX(1, N)                                    !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      CHARACTER (len=*), intent(in) :: NAME, OPTS
      integer,           intent(in) ::  ISPEC, N1, N2, N3, N4
!
!  Local variable declarations.
!
      logical :: CNAME, SNAME, TWOSTAGE
!
      integer :: I, IC, IZ, NB, NBMIN, NX
      integer :: LAENV
!
      character (len=1) :: C1
      character (len=2) :: C2, C4
      character (len=3) :: C3
      character (len=16) :: SUBNAM
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      GOTO ( 10,  10,  10,  80,  90, 100, 110, 120,                     &
     &      130, 140, 150, 160, 160, 160, 160, 160, 160) ISPEC
!
!  Invalid value for ISPEC.
!
      LAENV = -1
      RETURN
!
   10 CONTINUE
!
!  Convert NAME to upper case if the first character is lower case.
!
      LAENV = 1
      SUBNAM = NAME
      IC = ICHAR(SUBNAM(1:1))
      IZ = ICHAR('Z')
      IF (IZ.EQ.90 .OR. IZ.EQ.122) THEN
!
!  ASCII character set.
!
         IF (IC.GE.97 .AND. IC.LE.122) THEN
            SUBNAM(1:1) = CHAR(IC-32)
            DO I = 2, 6
              IC = ICHAR(SUBNAM(I:I))
              IF (IC.GE.97 .AND. IC.LE.122) SUBNAM(I:I) = CHAR(IC-32)
            END DO
         END IF
!
      ELSE IF (IZ.EQ.233 .OR. IZ.EQ.169) THEN
!
!  EBCDIC character set
!
         IF ((IC.GE.129 .AND. IC.LE.137) .OR.                           &
     &       (IC.GE.145 .AND. IC.LE.153) .OR.                           &
     &       (IC.GE.162 .AND. IC.LE.169)) THEN
            SUBNAM(1:1) = CHAR(IC+64)
            DO I = 2, 6
              IC = ICHAR(SUBNAM(I:I))
              IF ((IC.GE.129 .AND. IC.LE.137) .OR.                      &
     &            (IC.GE.145 .AND. IC.LE.153) .OR.                      &
     &            (IC.GE.162 .AND. IC.LE.169)) SUBNAM(I:I) = CHAR(IC+64)
            END DO
         END IF
!
      ELSE IF (IZ.EQ.218 .OR. IZ.EQ.250) THEN
!
!  Prime machines:  ASCII+128
!
         IF (IC.GE.225 .AND. IC.LE.250) THEN
            SUBNAM(1: 1) = CHAR(IC-32)
            DO I = 2, 6
              IC = ICHAR(SUBNAM(I: I))
              IF (IC.GE.225 .AND. IC.LE.250) SUBNAM(I:I) = CHAR(IC-32)
            END DO
         END IF
      END IF
!
      C1 = SUBNAM(1: 1)
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF (.NOT.(CNAME .OR. SNAME)) RETURN
      C2 = SUBNAM(2: 3)
      C3 = SUBNAM(4: 6)
      C4 = C3(2: 3)
      TWOSTAGE = LEN(SUBNAM).GE.11 .AND. SUBNAM(11:11).EQ.'2'
!
      GOTO (50, 60, 70) ISPEC
!
   50 CONTINUE
!
!  ISPEC = 1:  block size
!
!  In these examples, separate code is provided for setting NB for
!  real and complex.  We assume that NB will take the same value in
!  single or double precision.
!
      NB = 1
!
      IF (SUBNAM(2:6).EQ.'LAORH') THEN
!
!  This is for *LAORHR_GETRFNP routine
!
         IF (SNAME) THEN
             NB = 32
         ELSE
             NB = 32
         END IF
      ELSE IF (C2.EQ.'GE') THEN
         IF (C3.EQ.'TRF') THEN
            IF (SNAME) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF (C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.    &
     &            C3.EQ.'QLF') THEN
            IF (SNAME) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF (C3.EQ.'QR ') THEN
            IF (N3 .EQ. 1) THEN
               IF (SNAME) THEN
!
!  M*N
!
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               ELSE
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               END IF
            ELSE
               IF (SNAME) THEN
                  NB = 1
               ELSE
                  NB = 1
               END IF
            END IF
         ELSE IF (C3.EQ.'LQ ') THEN
            IF (N3 .EQ. 2) THEN
               IF (SNAME) THEN
!
!  M*N
!
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               ELSE
                  IF ((N1*N2.LE.131072).OR.(N1.LE.8192)) THEN
                     NB = N1
                  ELSE
                     NB = 32768/N2
                  END IF
               END IF
            ELSE
               IF (SNAME) THEN
                  NB = 1
               ELSE
                  NB = 1
               END IF
            END IF
         ELSE IF (C3.EQ.'HRD') THEN
            IF (SNAME) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF (C3.EQ.'BRD') THEN
            IF (SNAME) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF (C3.EQ.'TRI') THEN
            IF (SNAME) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF (SUBNAM(4:7).EQ.'QP3RK') THEN
            IF (SNAME) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      ELSE IF (C2.EQ.'PO') THEN
         IF (C3.EQ.'TRF') THEN
            IF (SNAME) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF (C2.EQ.'SY') THEN
         IF (C3.EQ.'TRF') THEN
            IF (SNAME) THEN
               IF (TWOSTAGE) THEN
                  NB = 192
               ELSE
                  NB = 64
               END IF
            ELSE
               IF (TWOSTAGE) THEN
                  NB = 192
               ELSE
                  NB = 64
               END IF
            END IF
         ELSE IF (SNAME .AND. C3.EQ.'TRD') THEN
            NB = 32
         ELSE IF (SNAME .AND. C3.EQ.'GST') THEN
            NB = 64
         END IF
      ELSE IF (CNAME .AND. C2.EQ.'HE') THEN
         IF (C3.EQ.'TRF') THEN
            IF (TWOSTAGE) THEN
               NB = 192
            ELSE
               NB = 64
            END IF
         ELSE IF (C3.EQ.'TRD') THEN
            NB = 32
         ELSE IF (C3.EQ.'GST') THEN
            NB = 64
         END IF
      ELSE IF (SNAME .AND. C2.EQ.'OR') THEN
         IF (C3(1:1).EQ.'G') THEN
            IF (C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR') THEN
               NB = 32
            END IF
         ELSE IF (C3(1:1).EQ.'M') THEN
            IF (C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR') THEN
               NB = 32
            END IF
         END IF
      ELSE IF (CNAME .AND. C2.EQ.'UN') THEN
         IF (C3(1:1).EQ.'G') THEN
            IF (C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR') THEN
               NB = 32
            END IF
         ELSE IF (C3(1:1).EQ.'M') THEN
            IF (C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR') THEN
               NB = 32
            END IF
         END IF
      ELSE IF (C2.EQ.'GB') THEN
         IF (C3.EQ.'TRF') THEN
            IF (SNAME) THEN
               IF (N4.LE.64) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF (N4.LE.64) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF (C2.EQ.'PB') THEN
         IF (C3.EQ.'TRF') THEN
            IF (SNAME) THEN
               IF (N2.LE.64) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF (N2.LE.64) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF (C2.EQ.'TR') THEN
         IF (C3.EQ.'TRI') THEN
            IF (SNAME) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF (C3.EQ.'EVC') THEN
            IF (SNAME) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF (C3.EQ.'SYL') THEN
!
!  The upper bound is to prevent overly aggressive scaling.
!
            IF (SNAME) THEN
               NB = MIN(MAX(48, INT((MIN(N1, N2) * 16) / 100)), 240)
            ELSE
               NB = MIN(MAX(24, INT((MIN(N1, N2) * 8) / 100)), 80)
            END IF
         END IF
      ELSE IF (C2.EQ.'LA') THEN
         IF (C3.EQ.'UUM') THEN
            IF (SNAME) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF (C3.EQ.'TRS') THEN
            IF (SNAME) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      ELSE IF (SNAME .AND. C2.EQ.'ST') THEN
         IF (C3.EQ.'EBZ') THEN
            NB = 1
         END IF
      ELSE IF (C2.EQ.'GG') THEN
         NB = 32
         IF (C3.EQ.'HD3') THEN
            IF (SNAME) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         END IF
      END IF
      LAENV = NB
      RETURN
!
   60 CONTINUE
!
!  ISPEC = 2:  minimum block size.
!
      NBMIN = 2
      IF (C2.EQ.'GE') THEN
         IF (C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.         &
     &       C3.EQ.'QLF') THEN
            IF (SNAME) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF (C3.EQ.'HRD') THEN
            IF (SNAME) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF (C3.EQ.'BRD') THEN
            IF (SNAME) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF (C3.EQ.'TRI') THEN
            IF (SNAME) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF (SUBNAM(4:7).EQ.'QP3RK') THEN
            IF (SNAME) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF

      ELSE IF (C2.EQ.'SY') THEN
         IF (C3.EQ.'TRF') THEN
            IF (SNAME) THEN
               NBMIN = 8
            ELSE
               NBMIN = 8
            END IF
         ELSE IF (SNAME .AND. C3.EQ.'TRD') THEN
            NBMIN = 2
         END IF
      ELSE IF (CNAME .AND. C2.EQ.'HE') THEN
         IF (C3.EQ.'TRD') THEN
            NBMIN = 2
         END IF
      ELSE IF (SNAME .AND. C2.EQ.'OR') THEN
         IF (C3(1:1).EQ.'G') THEN
            IF (C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR') THEN
               NBMIN = 2
            END IF
         ELSE IF (C3(1:1).EQ.'M') THEN
            IF (C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR') THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF (CNAME .AND. C2.EQ.'UN') THEN
         IF (C3(1:1).EQ.'G') THEN
            IF (C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR') THEN
               NBMIN = 2
            END IF
         ELSE IF (C3(1:1).EQ.'M') THEN
            IF (C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR') THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF (C2.EQ.'GG') THEN
         NBMIN = 2
         IF (C3.EQ.'HD3') THEN
            NBMIN = 2
         END IF
      END IF
      LAENV = NBMIN
      RETURN
!
   70 CONTINUE
!
!  ISPEC = 3:  crossover point
!
      NX = 0
      IF (C2.EQ.'GE') THEN
         IF (C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.         &
     &       C3.EQ.'QLF') THEN
            IF (SNAME) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF (C3.EQ.'HRD') THEN
            IF (SNAME) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF (C3.EQ.'BRD') THEN
            IF (SNAME) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF (SUBNAM(4:7).EQ.'QP3RK') THEN
            IF (SNAME) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF (C2.EQ.'SY') THEN
         IF (SNAME .AND. C3.EQ.'TRD') THEN
            NX = 32
         END IF
      ELSE IF (CNAME .AND. C2.EQ.'HE') THEN
         IF (C3.EQ.'TRD') THEN
            NX = 32
         END IF
      ELSE IF (SNAME .AND. C2.EQ.'OR') THEN
         IF (C3(1:1).EQ.'G') THEN
            IF (C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR') THEN
               NX = 128
            END IF
         END IF
      ELSE IF (CNAME .AND. C2.EQ.'UN') THEN
         IF (C3(1:1).EQ.'G') THEN
            IF (C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.         &
     &          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.         &
     &          C4.EQ.'BR') THEN
               NX = 128
            END IF
         END IF
      ELSE IF (C2.EQ.'GG') THEN
         NX = 128
         IF (C3.EQ.'HD3') THEN
            NX = 128
         END IF
      ELSE IF (C2.EQ.'LA') THEN
         IF (C3.EQ.'RFT') THEN
            NX = 64
         END IF
      END IF
      LAENV = NX
      RETURN
!
   80 CONTINUE
!
!  ISPEC = 4:  number of shifts (used by xHSEQR)
!
      LAENV = 6
      RETURN
!
   90 CONTINUE
!
!  ISPEC = 5:  minimum column dimension (not used)
!
      LAENV = 2
      RETURN
!
  100 CONTINUE
!
!  ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
!
      LAENV = INT(REAL(MIN(N1, N2))*1.6E0)
      RETURN
!
  110 CONTINUE
!
!  ISPEC = 7:  number of processors (not used)
!
      LAENV = 1
      RETURN
!
  120 CONTINUE
!
!  ISPEC = 8:  crossover point for multishift (used by xHSEQR)
!
      LAENV = 50
      RETURN
!
  130 CONTINUE
!
!  ISPEC = 9:  maximum size of the subproblems at the bottom of the
!              computation tree in the divide-and-conquer algorithm
!              (used by xGELSD and xGESDD)
!
      LAENV = 25
      RETURN
!
  140 CONTINUE
!
!  ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to
!              trap.
!
!     LAENV = 0
      LAENV = 1
      IF (LAENV.EQ.1) THEN
         LAENV = IEEECK(1, 0.0_dp, 1.0_dp)
      END IF
      RETURN
!
  150 CONTINUE
!
!  ISPEC = 11: ieee infinity arithmetic can be trusted not to trap
!
!     LAENV = 0
      LAENV = 1
      IF (LAENV.EQ.1) THEN
         LAENV = IEEECK(0, 0.0_dp, 1.0_dp)
      END IF
      RETURN
!
  160 CONTINUE
!
!  12 <= ISPEC <= 17: xHSEQR or related subroutines.
!
      LAENV = IPARMQ(ISPEC, NAME, OPTS, N1, N2, N3, N4)
!
      RETURN
      END FUNCTION ILAENV
!
      FUNCTION IPARMQ (ISPEC, NAME, OPTS, N, ILO, IHI, LWORK)           &
     &                RESULT (PARMQ)
!
!=======================================================================
!                                                                      !
!  IPARMQ sets problem and machine dependent parameters useful for     !
!  xHSEQR and related subroutines for eigenvalue problems. It is       !
!  called whenever IPARMQ is called with 12 <= ISPEC <= 16             !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  ISPEC    [in]      INTEGER                                          !
!           ISPEC specifies which tunable parameter IPARMQ should      !
!           return.                                                    !
!                                                                      !
!           ISPEC=12: (INMIN)  Matrices of order nmin or less          !
!                     are sent directly to xLAHQR, the implicit        !
!                     double shift QR algorithm.  NMIN must be         !
!                     at least 11.                                     !
!                                                                      !
!           ISPEC=13: (INWIN)  Size of the deflation window.           !
!                     This is best set greater than or equal to        !
!                     the number of simultaneous shifts NS.            !
!                     Larger matrices benefit from larger deflation    !
!                     windows.                                         !
!                                                                      !
!           ISPEC=14: (INIBL) Determines when to stop nibbling and     !
!                     invest in an (expensive) multi-shift QR sweep.   !
!                     If the aggressive early deflation subroutine     !
!                     finds LD converged eigenvalues from an order     !
!                     NW deflation window and LD > (NW*NIBBLE)/100,    !
!                     then the next QR sweep is skipped and early      !
!                     deflation is applied immediately to the          !
!                     remaining active diagonal block.  Setting        !
!                     IPARMQ(ISPEC=14) = 0 causes TTQRE to skip a      !
!                     multi-shift QR sweep whenever early deflation    !
!                     finds a converged eigenvalue.  Setting           !
!                     IPARMQ(ISPEC=14) greater than or equal to 100    !
!                     prevents TTQRE from skipping a multi-shift       !
!                     QR sweep.                                        !
!                                                                      !
!           ISPEC=15: (NSHFTS) The number of simultaneous shifts in    !
!                     a multi-shift QR iteration.                      !
!                                                                      !
!           ISPEC=16: (IACC22) IPARMQ is set to 0, 1 or 2 with the     !
!                     following meanings.                              !
!                     0:  During the multi-shift QR/QZ sweep,          !
!                         blocked eigenvalue reordering, blocked       !
!                         Hessenberg-triangular reduction,             !
!                         reflections and/or rotations are not         !
!                         accumulated when updating the                !
!                         far-from-diagonal matrix entries.            !
!                     1:  During the multi-shift QR/QZ sweep,          !
!                         blocked eigenvalue reordering, blocked       !
!                         Hessenberg-triangular reduction,             !
!                         reflections and/or rotations are             !
!                         accumulated, and matrix-matrix               !
!                         multiplication is used to update the         !
!                         far-from-diagonal matrix entries.            !
!                     2:  During the multi-shift QR/QZ sweep,          !
!                         blocked eigenvalue reordering, blocked       !
!                         Hessenberg-triangular reduction,             !
!                         reflections and/or rotations are             !
!                         accumulated, and 2-by-2 block structure      !
!                         is exploited during matrix-matrix            !
!                         multiplies.                                  !
!                     (If xTRMM is slower than xGEMM, then             !
!                      IPARMQ(ISPEC=16)=1 may be more efficient than   !
!                      IPARMQ(ISPEC=16)=2 despite the greater level of !
!                      arithmetic work implied by the latter choice.)  !
!                                                                      !
!           ISPEC=17: (ICOST) An estimate of the relative cost of      !
!                     flops within the near-the-diagonal shift chase   !
!                     compared to flops within the BLAS calls of a QZ  !
!                     sweep.                                           !
!                                                                      !
!  NAME     [in]      CHARACTER string                                 !
!           Name of the calling subroutine                             !
!                                                                      !
!  OPTS     [in]      CHARACTER string                                 !
!           This is a concatenation of the string arguments to         !
!           TTQRE.                                                     !
!                                                                      !
!  N        [in]      INTEGER                                          !
!           The order of the Hessenberg matrix H.                      !
!                                                                      !
!  ILO      [in]      INTEGER                                          !
!                                                                      !
!  IHI      [in]      INTEGER                                          !
!           It is assumed that H is already upper triangular           !
!           in rows and columns 1:ILO-1 and IHI+1:N.                   !
!                                                                      !
!  LWORK    [in]      INTEGER                                          !
!           The amount of workspace available.                         !
!                                                                      !
!  Authors:                                                            !
!                                                                      !
!  Univ. of Tennessee                                                  !
!  Univ. of California Berkeley                                        !
!  Univ. of Colorado Denver                                            !
!  NAG Ltd.                                                            !
!                                                                      !
!  Further Details:                                                    !
!                                                                      !
!  Little is known about how best to choose these parameters.          !
!  It is possible to use different values of the parameters            !
!  for each of CHSEQR, DHSEQR, SHSEQR and ZHSEQR.                      !
!                                                                      !
!  It is probably best to choose different parameters for              !
!  different matrices and different parameters at different            !
!  times during the iteration, but this has not been                   !
!  implemented --- yet.                                                !
!                                                                      !
!  The best choices of most of the parameters depend                   !
!  in an ill-understood way on the relative execution                  !
!  rate of xLAQR3 and xLAQR5 and on the nature of each                 !
!  particular eigenvalue problem.  Experiment may be the               !
!  only practical way to determine which choices are most              !
!  effective.                                                          !
!                                                                      !
!  Following is a list of default values supplied by IPARMQ.           !
!  These defaults may be adjusted in order to attain better            !
!  performance in any particular computational environment.            !
!                                                                      !
!    IPARMQ(ISPEC=12) The xLAHQR vs xLAQR0 crossover point.            !
!                     Default: 75. (Must be at least 11.)              !
!                                                                      !
!    IPARMQ(ISPEC=13) Recommended deflation window size.               !
!                     This depends on ILO, IHI and NS, the             !
!                     number of simultaneous shifts returned           !
!                     by IPARMQ(ISPEC=15).  The default for            !
!                     (IHI-ILO+1) <= 500 is NS.  The default           !
!                     for (IHI-ILO+1) > 500 is 3*NS/2.                 !
!                                                                      !
!    IPARMQ(ISPEC=14) Nibble crossover point.  Default: 14.            !
!                                                                      !
!    IPARMQ(ISPEC=15) Number of simultaneous shifts, NS.               !
!                     a multi-shift QR iteration.                      !
!                                                                      !
!                     If IHI-ILO+1 is ...                              !
!                                                                      !
!                     greater than      ...but less    ... the         !
!                     or equal to ...      than        default is      !
!                             0               30       NS =   2+       !
!                            30               60       NS =   4+       !
!                            60              150       NS =  10        !
!                           150              590       NS =  **        !
!                           590             3000       NS =  64        !
!                          3000             6000       NS = 128        !
!                          6000             infinity   NS = 256        !
!                                                                      !
!                 (+)  By default matrices of this order are           !
!                      passed to the implicit double shift routine     !
!                      xLAHQR.  See IPARMQ(ISPEC=12) above.   These    !
!                      values of NS are used only in case of a rare    !
!                      xLAHQR failure.                                 !
!                                                                      !
!                 (**) The asterisks (**) indicate an ad-hoc           !
!                      function increasing from 10 to 64.              !
!                                                                      !
!    IPARMQ(ISPEC=16) Select structured matrix multiply.               !
!                     (See ISPEC=16 above for details.)                !
!                     Default: 3.                                      !
!                                                                      !
!    IPARMQ(ISPEC=17) Relative cost heuristic for blocksize selection. !
!                     Expressed as a percentage.                       !
!                     Default: 10.                                     !
!                                                                      !
!  -- LAPACK auxiliary routine --                                      !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: IHI, ILO, ISPEC, LWORK, N
!
      character (len=*), intent(in) :: NAME, OPTS
!
!  Local variable declarations.
!
      integer,   parameter :: INMIN = 12
      integer,   parameter :: INWIN = 13
      integer,   parameter :: INIBL = 14
      integer,   parameter :: ISHFTS = 15
      integer,   parameter :: IACC22 = 16
      integer,   parameter :: ICOST = 17
!
      integer,   parameter :: NMIN = 75
      integer,   parameter :: K22MIN = 14
      integer,   parameter :: KACMIN = 14
      integer,   parameter :: NIBBLE = 14
      integer,   parameter :: KNWSWP = 500
      integer,   parameter :: RCOST = 10
!
      real (dp), parameter :: TWO = 2.0_dp
!
      integer           :: NH, NS
      integer           :: I, IC, IZ
      integer           :: PARMQ
!
      character (len=6) :: SUBNAM
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------

      IF ((ISPEC.EQ.ISHFTS) .OR. (ISPEC.EQ.INWIN) .OR.                  &
     &    (ISPEC.EQ.IACC22)) THEN
!
!  Set the number simultaneous shifts ====
!
         NH = IHI - ILO + 1
         NS = 2
         IF (NH.GE.30)   NS = 4
         IF (NH.GE.60)   NS = 10
         IF (NH.GE.150)  NS = MAX(10, NH/NINT(LOG(REAL(NH)) / LOG(TWO)))
         IF (NH.GE.590)  NS = 64
         IF (NH.GE.3000) NS = 128
         IF (NH.GE.6000) NS = 256
         NS = MAX(2, NS-MOD(NS,2))
      END IF
!
      IF (ISPEC.EQ.INMIN) THEN
!
!  Matrices of order smaller than NMIN get sent
!    .     to xLAHQR, the classic double shift algorithm.
!    .     This must be at least 11.
!
         PARMQ = NMIN
!
      ELSE IF (ISPEC.EQ.INIBL) THEN
!
!  INIBL: skip a multi-shift qr iteration and
!    .    whenever aggressive early deflation finds
!    .    at least (NIBBLE*(window size)/100) deflations.
!
         PARMQ = NIBBLE
!
      ELSE IF (ISPEC.EQ.ISHFTS) THEN
!
!  NSHFTS: The number of simultaneous shifts.
!
         PARMQ = NS
!
      ELSE IF (ISPEC.EQ.INWIN) THEN
!
!  NW: deflation window size.
!
         IF (NH.LE.KNWSWP) THEN
            PARMQ = NS
         ELSE
            PARMQ = 3*NS / 2
         END IF
!
      ELSE IF (ISPEC.EQ.IACC22) THEN
!
!  IACC22: Whether to accumulate reflections
!    .     before updating the far-from-diagonal elements
!    .     and whether to use 2-by-2 block structure while
!    .     doing it.  A small amount of work could be saved
!    .     by making this choice dependent also upon the
!    .     NH=IHI-ILO+1.
!
!  Convert NAME to upper case if the first character is lower case.
!
         PARMQ = 0
         SUBNAM = NAME
         IC = ICHAR(SUBNAM(1:1))
         IZ = ICHAR('Z')
         IF (IZ.EQ.90.OR.IZ.EQ.122) THEN
!
!  ASCII character set.
!
            IF (IC.GE.97 .AND. IC.LE.122) THEN
               SUBNAM(1:1) = CHAR(IC-32)
               DO I = 2, 6
                  IC = ICHAR(SUBNAM(I:I))
                  IF (IC.GE.97 .AND. IC.LE.122) THEN
                    SUBNAM(I:I) = CHAR(IC-32)
                  END IF
               END DO
            END IF
!
         ELSE IF (IZ.EQ.233 .OR. IZ.EQ.169) THEN
!
!  EBCDIC character set.
!
            IF ((IC.GE.129 .AND. IC.LE.137) .OR.                        &
     &          (IC.GE.145 .AND. IC.LE.153) .OR.                        &
     &          (IC.GE.162 .AND. IC.LE.169)) THEN
               SUBNAM(1: 1) = CHAR(IC+64)
               DO I = 2, 6
                  IC = ICHAR(SUBNAM(I:I))
                  IF ((IC.GE.129 .AND. IC.LE.137) .OR.                  &
     &                (IC.GE.145 .AND. IC.LE.153) .OR.                  &
     &                (IC.GE.162 .AND. IC.LE.169)) THEN
                     SUBNAM(I:I) = CHAR(IC+64)
                  END IF
               END DO
            END IF
!
         ELSE IF (IZ.EQ.218 .OR. IZ.EQ.250) THEN
!
!  Prime machines:  ASCII+128.
!
            IF (IC.GE.225 .AND. IC.LE.250) THEN
               SUBNAM(1:1) = CHAR(IC-32)
               DO I = 2, 6
                  IC = ICHAR(SUBNAM(I:I))
                  IF (IC.GE.225 .AND. IC.LE.250) THEN
                     SUBNAM(I:I) = CHAR(IC-32)
                  END IF
               END DO
            END IF
         END IF
!
         IF (SUBNAM(2:6).EQ.'GGHRD' .OR.                                &
     &       SUBNAM(2:6).EQ.'GGHD3') THEN
            PARMQ = 1
            IF (NH.GE.K22MIN) PARMQ = 2
         ELSE IF (SUBNAM(4:6).EQ.'EXC') THEN
            IF (NH.GE.KACMIN) PARMQ = 1
            IF (NH.GE.K22MIN) PARMQ = 2
         ELSE IF (SUBNAM(2:6).EQ.'HSEQR' .OR.                           &
     &            SUBNAM(2:5).EQ.'LAQR') THEN
            IF (NS.GE.KACMIN) PARMQ = 1
            IF (NS.GE.K22MIN) PARMQ = 2
         END IF
!
      ELSE IF (ISPEC.EQ.ICOST) THEN
!
!  Relative cost of near-the-diagonal chase vs BLAS updates.
!
         PARMQ = RCOST
      ELSE
!
!  Invalid value of ISPEC.
!
         PARMQ = -1
!
      END IF
!
      RETURN
      END FUNCTION IPARMQ
!
      SUBROUTINE XERBLA (SRNAME, INFO)
!
!=======================================================================
!                                                                      !
!  XERBLA  is an error handler for the LAPACK routines.                !
!  It is called by an LAPACK routine if an input parameter has an      !
!  invalid value.  A message is printed and execution stops.           !
!                                                                      !
!  Installers may consider modifying the STOP statement in order to    !
!  call system-specific exception-handling facilities.                 !
!                                                                      !
!  Arguments:                                                          !
!                                                                      !
!  SRNAME  (input) CHARACTER*6                                         !
!          The name of the routine which called XERBLA.                !
!                                                                      !
!  INFO    (input) INTEGER                                             !
!          The position of the invalid parameter in the parameter list !
!          of the calling routine.                                     !
!                                                                      !
!  -- LAPACK auxiliary routine (version 2.0) --                        !
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,      !
!     Courant Institute, Argonne National Lab, and Rice University     !
!     September 30, 1994                                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: INFO
!
      character (len=*), intent(in) :: SRNAME
!
!-----------------------------------------------------------------------
!  Executable Statements.
!-----------------------------------------------------------------------
!
      PRINT 10, SRNAME, INFO
  10  FORMAT (' LAPACK_MOD: ** On entry to ', a, ' parameter number ',  &
     &        i0,' had an illegal value.')
      STOP
!
      END SUBROUTINE XERBLA
!
      END MODULE lapack_mod
