#include "cppdefs.h"
      SUBROUTINE mod_arrays (allocate_vars)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2018 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine routine allocates and initializa model state arrays    !
!  for each nested and/or multiple connected grids.                    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_collect
#endif
      USE mod_netcdf, ONLY : Matts, Mdims, Mvars, NvarD, NvarA
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
      USE mod_average, ONLY : allocate_average, initialize_average
#endif
      USE mod_boundary, ONLY : allocate_boundary, initialize_boundary
      USE mod_clima, ONLY : allocate_clima, initialize_clima
#ifdef SOLVE3D
      USE mod_coupling, ONLY : allocate_coupling, initialize_coupling
#endif
#ifdef DIAGNOSTICS
      USE mod_diags, ONLY : allocate_diags, initialize_diags
#endif
      USE mod_forces, ONLY : allocate_forces, initialize_forces
      USE mod_grid, ONLY : allocate_grid, initialize_grid
      USE mod_mixing, ONLY : allocate_mixing, initialize_mixing
#ifdef NESTING
      USE mod_nesting, ONLY : allocate_nesting, initialize_nesting
#endif
      USE mod_ocean, ONLY : allocate_ocean, initialize_ocean
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sedbed, ONLY : allocate_sedbed, initialize_sedbed
#endif
      USE mod_sources, ONLY : allocate_sources
#if defined SSH_TIDES || defined UV_TIDES
      USE mod_tides, ONLY : allocate_tides, initialize_tides
#endif
#ifdef BBL_MODEL
      USE mod_bbl, ONLY : allocate_bbl, initialize_bbl
#endif
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: allocate_vars
!
!  Local variable declarations.
!
      logical :: LallocateClima
!
      integer :: ng, thread, tile
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
      integer :: Nlevels, Ntiles
!
      integer, parameter :: model = 0
!
      real(r8) :: bytefac, megabytefac, size1d, size2d
      real(r8) :: Avalue, sumAsize, sumDsize, totalAsize, totalDsize
!
      real(r8), parameter :: spv = 0.0_r8
!
      real(r8), allocatable ::  Asize(:,:)
      real(r8), allocatable ::  Dsize(:,:)
      real(r8), allocatable ::  IOsize(:,:)

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn on allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_on (ng, iNLM, 1, __LINE__, __FILE__)
        END DO
!$OMP BARRIER
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Allocate model structures.
!-----------------------------------------------------------------------
!
      IF (allocate_vars) then
#ifdef DISTRIBUTE
        tile=MyRank
#else
        tile=0
#endif

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI
        LallocateClima=.TRUE.
#else
        LallocateClima=.FALSE.
#endif

        DO ng=1,Ngrids
!$OMP MASTER
          LBi=BOUNDS(ng)%LBi(tile)
          UBi=BOUNDS(ng)%UBi(tile)
          LBj=BOUNDS(ng)%LBj(tile)
          UBj=BOUNDS(ng)%UBj(tile)
          LBij=BOUNDS(ng)%LBij
          UBij=BOUNDS(ng)%UBij
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
          CALL allocate_average (ng, LBi, UBi, LBj, UBj)
#endif
          CALL allocate_boundary (ng)
#ifdef BBL_MODEL
          CALL allocate_bbl (ng, LBi, UBi, LBj, UBj)
#endif
          IF (LallocateClima.or.Lclimatology(ng)) THEN
            CALL allocate_clima (ng, LBi, UBi, LBj, UBj)
          END IF
#ifdef SOLVE3D
          CALL allocate_coupling (ng, LBi, UBi, LBj, UBj)
#endif
#ifdef DIAGNOSTICS
          CALL allocate_diags (ng, LBi, UBi, LBj, UBj)
#endif
          CALL allocate_forces (ng, LBi, UBi, LBj, UBj)
          CALL allocate_grid (ng, LBi, UBi, LBj, UBj, LBij, UBij)
          CALL allocate_mixing (ng, LBi, UBi, LBj, UBj)
          CALL allocate_ocean (ng, LBi, UBi, LBj, UBj)
#if defined SEDIMENT || defined BBL_MODEL
          CALL allocate_sedbed (ng, LBi, UBi, LBj, UBj)
#endif
#if defined SSH_TIDES || defined UV_TIDES
          CALL allocate_tides (ng, LBi, UBi, LBj, UBj)
#endif
          IF (LuvSrc(ng).or.LwSrc(ng).or.ANY(LtracerSrc(:,ng))) THEN
            CALL allocate_sources (ng)
          END IF
!$OMP END MASTER
!$OMP BARRIER
        END DO
#ifdef NESTING
!
!  Allocate and initialized contact points boundaty structure. It
!  needs to be delayed to the end because we need "LBC_apply" to
!  allocated in "mod_boundary" for all nested grid.
!
!$OMP MASTER
        CALL allocate_nesting
!$OMP END MASTER
!$OMP BARRIER
#endif
      END IF
!
!-----------------------------------------------------------------------
!  Intialize variables within structures for each grid.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
#ifdef NESTING
        IF (ng.eq.1) THEN
          CALL initialize_nesting
        END IF
#endif
        DO tile=first_tile(ng),last_tile(ng),+1
#if defined AVERAGES    || \
   (defined AD_AVERAGES && defined ADJOINT) || \
   (defined RP_AVERAGES && defined TL_IOMS) || \
   (defined TL_AVERAGES && defined TANGENT)
          CALL initialize_average (ng, tile)
#endif
#ifdef BBL_MODEL
          CALL initialize_bbl (ng, tile)
#endif
          CALL initialize_boundary (ng, tile, model)

          IF (LallocateClima.or.Lclimatology(ng)) THEN
            CALL initialize_clima (ng, tile)
          END IF
#ifdef SOLVE3D
          CALL initialize_coupling (ng, tile, model)
#endif
#ifdef DIAGNOSTICS
          CALL initialize_diags (ng, tile)
#endif
          CALL initialize_forces (ng, tile, model)
          CALL initialize_grid (ng, tile, model)
          CALL initialize_mixing (ng, tile, model)
          CALL initialize_ocean (ng, tile, model)
#if defined SEDIMENT || defined BBL_MODEL
          CALL initialize_sedbed (ng, tile, model)
#endif
#if defined SSH_TIDES || defined UV_TIDES
          CALL initialize_tides (ng, tile)
#endif
        END DO
!$OMP BARRIER
      END DO
!
!-----------------------------------------------------------------------
!  Report estimate of dynamic memory and automatic memory requirements.
!-----------------------------------------------------------------------
!
!  Allocate and initialize.
!
!$OMP MASTER
      Ntiles=MAXVAL(NtileI)*MAXVAL(NtileJ)-1
      IF (.not.allocated(Asize)) THEN
        allocate ( Asize(0:Ntiles,Ngrids) )
        Asize=spv
      END IF
      IF (.not.allocated(Dsize)) THEN
        allocate ( Dsize(0:Ntiles,Ngrids) )
        Dsize=spv
      END IF
      IF (.not.allocated(IOsize)) THEN
        allocate ( IOsize(0:Ntiles,Ngrids) )
        IOsize=spv
      END IF
!
!  Determine size floating-point arrays in bytes.  We could use the
!  Fortran 2008 standard function STORAGE_SIZE.  However since ROMS
!  ROMS by default is double-precision, we just set it value to 64 bits
!  or 8 bytes (1 byte = 8 bits).  The number of array elements is
!  multiplied by the megabytes factor.
!
      bytefac=8.0_r8                    ! double precision (64 bits)
      megabytefac=bytefac*1.0E-6_r8     ! 1 Mb = 1.0E+6 bytes (SI units)
!
!  Add static memory requirements for processing NetCDF data.  The
!  variables are declared in "mod_netcdf".  Notice that a single
!  character has a size of eight bits (1 byte).
!
      Dmem(1)=Dmem(1)+REAL(Matts,r8)               ! att_kind
      Dmem(1)=Dmem(1)+2.0_r8*REAL(Mdims,r8)        ! dim_id,dim_size
      Dmem(1)=Dmem(1)+5.0_r8*REAL(Mvars,r8)        ! var_*
      Dmem(1)=Dmem(1)+REAL(NvarD*Mvars,r8)         ! var_dim
      Dmem(1)=Dmem(1)+2.0_r8*REAL(NvarD,r8)        ! var_Dids,var_Dsize
      Dmem(1)=Dmem(1)+2.0_r8*REAL(NvarA,r8)        ! var_Aint,var_Afloat
      Dmem(1)=Dmem(1)+0.125_r8*REAL(40*Matts,r8)   ! att_name
      Dmem(1)=Dmem(1)+0.125_r8*REAL(40*Mdims,r8)   ! dim_name
      Dmem(1)=Dmem(1)+0.125_r8*REAL(40*Mvars,r8)   ! dim_name
      Dmem(1)=Dmem(1)+0.125_r8*REAL(40*NvarA,r8)   ! var_Aname
      Dmem(1)=Dmem(1)+0.125_r8*REAL(40*NvarD,r8)   ! var_Dname
      Dmem(1)=Dmem(1)+0.125_r8*REAL(1024*NvarA,r8) ! var_Achar
!
!  Estimate automatic memory requirements (megabytes) by looking at the
!  routines that use it most, like step2d, step3d_t, or NetCDF I/O.
!  (see memory.txt for more information).
!
      DO ng=1,Ngrids
        DO tile=0,NtileI(ng)*NtileJ(ng)-1
#ifdef NESTING
          IminS=BOUNDS(ng)%Istr(tile)-4
          ImaxS=BOUNDS(ng)%Iend(tile)+3
          JminS=BOUNDS(ng)%Jstr(tile)-4
          JmaxS=BOUNDS(ng)%Jend(tile)+3
#else
          IminS=BOUNDS(ng)%Istr(tile)-3
          ImaxS=BOUNDS(ng)%Iend(tile)+3
          JminS=BOUNDS(ng)%Jstr(tile)-3
          JmaxS=BOUNDS(ng)%Jend(tile)+3
#endif
          size1d=REAL((ImaxS-IminS),r8)
          size2d=REAL((ImaxS-IminS)*(JmaxS-JminS),r8)
#ifdef SOLVE3D
          Asize(tile,ng)=megabytefac*                                   &
     &                   (4.0_r8*size1d*REAL(N(ng)+1,r8)+               &
     &                    7.0_r8*size2d+                                &
     &                    5.0_r8*size2d*REAL(N(ng),r8)+                 &
     &                    1.0_r8*size2d*REAL(N(ng)*NT(ng),r8))
# if !(defined PARALLEL_IO && defined DISTRIBUTE)
#  ifdef INLINE_2DIO
          Nlevels=1
#  else
          Nlevels=N(ng)+1
#  endif
          IOsize(tile,ng)=megabytefac*2.0_r8*                           &
     &                    REAL(2+(Lm(ng)+2)*(Mm(ng)+2)*(Nlevels),r8)
# else
          IOsize(tile,ng)=megabytefac*                                  &
     &                    REAL(2+(Lm(ng)+2)*(Mm(ng)+2)*(N(ng)+1),r8)
# endif
#else
          Asize(tile,ng)=megabytefac*                                   &
     &                   (38.0_r8*size2d)
          IOsize(tile,ng)=megabytefac*                                  &
     &                    REAL(2+(Lm(ng)+2)*(Mm(ng)+2),r8)
#endif
        END DO
      END DO
!
!  Determine total maximum value of dynamic-memory and automatic-memory
!  requirements, and convert number of array element to bytes.
!
#ifdef DISTRIBUTE
      Dsize(0:numthreads-1,1:Ngrids)=spv
      DO ng=1,Ngrids
        Dsize(MyRank,ng)=megabytefac*Dmem(ng)
        CALL mp_collect (ng, iNLM, numthreads, spv, Dsize(0:,ng))
      END DO
#else
      Dsize(0:numthreads-1,1:Ngrids)=spv
      DO ng=1,Ngrids
        Dsize(0,ng)=megabytefac*Dmem(ng)
      END DO
#endif
!
!  Report dynamic and automatic memory reuirements.
!
      IF (Master) THEN
        WRITE (stdout,"(/,80('>'))")
        totalAsize=0.0_r8
        totalDsize=0.0_r8
        DO ng=1,Ngrids
          sumAsize=0.0_r8
          sumDsize=0.0_r8
#ifdef SOLVE3D
          WRITE (stdout,10) ng, Lm(ng), Mm(ng), N(ng),                  &
     &                      NtileI(ng), NtileJ(ng)
#else
          WRITE (stdout,10) ng, Lm(ng), Mm(ng),                         &
     &                      NtileI(ng), NtileJ(ng)
#endif
          DO tile=0,NtileI(ng)*NtileJ(ng)-1
            Avalue=MAX(Asize(tile,ng), IOsize(tile,ng))
            sumAsize=sumAsize+Avalue
            sumDsize=sumDsize+Dsize(tile,ng)
            WRITE (stdout,20) tile, Dsize(tile,ng), Avalue,             &
     &                        Dsize(tile,ng)+Avalue
          END DO
          totalAsize=totalAsize+sumAsize
          totalDsize=totalDsize+sumDsize
          IF (Ngrids.gt.1) THEN
            WRITE (stdout,30) '  SUM', sumDsize, sumAsize,              &
     &                                 sumAsize+sumDsize
          ELSE
            WRITE (stdout,30) 'TOTAL', sumDsize, sumAsize,              &
     &                                 sumAsize+sumDsize
          END IF
        END DO
        IF (Ngrids.gt.1) THEN
            WRITE (stdout,30) 'TOTAL', totalDsize, totalAsize,          &
     &                                 totalAsize+totalDsize
        END IF
        WRITE (stdout,"(/,80('<'))")
      END IF
!
!  Deallocate dynamic and automatic memory local arrays.
!
      IF (allocated(Asize))  deallocate ( Asize )
      IF (allocated(Dsize))  deallocate ( Dsize )
      IF (allocated(IOsize)) deallocate ( IOsize )
!$OMP END MASTER
!$OMP BARRIER

#ifdef PROFILE
!
!-----------------------------------------------------------------------
!  Turn off allocation time wall clock.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iNLM, 1, __LINE__, __FILE__)
        END DO
!$OMP BARRIER
      END DO
#endif
!
#ifdef SOLVE3D
 10   FORMAT (/,' Dynamic and Automatic memory (MB) usage for Grid ',    &
     &          i2.2,':',2x,i0,'x',i0,'x',i0,2x,'tiling: ',i0,'x',i0,    &
     &          /,/,5x,'tile',10x,'Dynamic',8x,'Automatic',12x,          &
     &          'USAGE',/)
#else
 10   FORMAT (/,' Dynamic and Automatic memory (MB) usage for Grid ',    &
     &          i2.2,':',2x,i0,'x',i0,2x,'tiling: ',i0,'x',i0,           &
     &          /,/,5x,'tile',10x,'Dynamic',8x,'Automatic',12x,          &
     &          'USAGE',/)
#endif
 20   FORMAT (4x,i5,3(4x,f13.2))
 30   FORMAT (/,4x,a,3(4x,f13.2))
!
      RETURN
      END SUBROUTINE mod_arrays
