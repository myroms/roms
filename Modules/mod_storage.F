#include "cppdefs.h"
      MODULE mod_storage

#if defined PROPAGATOR

!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  ROMS/TOMS Scratch Arrays:                                           !
!                                                                      !
!  Mstate       Maximum size of state vector.                          !
!  Nstate       Size of state vector (water points only) per grid.     ! 
!  NCV          Number of Lanczos vectors to compute.                  !
!  NEV          Number of eigenvalues to compute.                      !
!                                                                      !
!  Bvec         Lanczos/Arnoldi basis vectors.                         !
!  Rvalue       Ritz eigenvalues.                                      !
!  Rvector      Ritz eigenvectors.                                     !
!  SworkD       State work array for reverse communications.           !
!  SworkL       State work array.                                      !
!  norm         Euclidean norm.                                        !
!  resid        Initial/final residual vector.                         !
!  select       Logical switch of Ritz vectors to compute.             !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none

        integer, dimension(Ngrids) :: Nstate
        integer :: Mstate                           ! MAXVAL(Nstate)
        integer :: NCV                              ! 2*NEV
        integer :: NEV

        logical, allocatable :: select(:)           ! [1:NCV]

        real(r8), allocatable :: Bvec(:,:)          ! [1:Mstate,1:NCV]
        real(r8), allocatable :: Rvalue(:)          ! [1:NEV]
        real(r8), allocatable :: Rvector(:,:)       ! [1:Mstate,1:NEV]
        real(r8), allocatable :: SworkD(:)          ! [1:3*Mstate]
        real(r8), allocatable :: SworkL(:)          ! [1:NCV*(NCV+8)]
        real(r8), allocatable :: norm(:)            ! [1:NEV]
        real(r8), allocatable :: resid(:)           ! [1:Mstate]
   
      CONTAINS

      SUBROUTINE allocate_storage (ng)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine allocates and initialize module variables.             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng
!
!  Local variable declarations
!
      real(r8), parameter :: IniVal = 0.0_r8
!
!-----------------------------------------------------------------------
!  Allocate and initialize module variables.
!-----------------------------------------------------------------------
!
      IF (ng.eq.1) THEN

        allocate ( select(NCV) )
        select = .TRUE.

        allocate ( Bvec(Mstate,NCV) )
        Bvec = IniVal

        allocate ( Rvalue(NEV) )
        Rvalue = IniVal

        allocate ( Rvector(Mstate,NEV) )
        Rvector = IniVal

        allocate ( SworkD(3*Mstate) )
        SworkD = IniVal

        allocate ( SworkL(NCV*(NCV+8)) )
        SworkD = IniVal

        allocate ( norm(NEV) )
        norm = IniVal

        allocate ( resid(Mstate) )
        resid = IniVal

      END IF

      RETURN
      END SUBROUTINE allocate_storage
#endif
      END MODULE mod_storage
