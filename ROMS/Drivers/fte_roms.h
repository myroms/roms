      MODULE roms_kernel_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group            Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  ROMS Finite Time Eigenmodes (FTE) Driver:                           !
!                                                                      !
!  This driver computes the finite time eigenmodes of the propagator   !
!  R(0,t)  linearized about a time  evolving  circulation.  They are   !
!  often referred to as the finite time normal modes and are used to   !
!  measure the asymptotic stability of the circulation.                !
!                                                                      !
!  These  routines  control the  initialization,  time-stepping, and   !
!  finalization of  ROMS  model following ESMF conventions:            !
!                                                                      !
!     ROMS_initialize                                                  !
!     ROMS_run                                                         !
!     ROMS_finalize                                                    !
!                                                                      !
!  WARNING: This driver cannot be run from the ESMF super-structure.   !
!  =======  Therefore, ESMF coupling is not possible.                  !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Moore, A.M. et al., 2004: A comprehensive ocean prediction and    !
!      analysis system based on the tangent linear and adjoint of a    !
!      regional ocean model, Ocean Modelling, 7, 227-258.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_arrays
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
#if defined PIO_LIB && defined DISTRIBUTE
      USE mod_pio_netcdf
#endif
      USE mod_scalars
      USE mod_stepping
      USE mod_storage
!
      USE propagator_mod
!
      USE close_io_mod,      ONLY : close_file, close_inp, close_out
#ifdef CHECKPOINTING
      USE def_gst_mod,       ONLY : def_gst
      USE get_gst_mod,       ONLY : get_gst
#endif
      USE inp_par_mod,       ONLY : inp_par
#ifdef MCT_LIB
# ifdef ATM_COUPLING
      USE mct_coupler_mod,   ONLY : initialize_ocn2atm_coupling
# endif
# ifdef WAV_COUPLING
      USE mct_coupler_mod,   ONLY : initialize_ocn2wav_coupling
# endif
#endif
      USE packing_mod,       ONLY : c_norm2
      USE packing_mod,       ONLY : r_norm2
      USE stdout_mod,        ONLY : Set_StdOutUnit, stdout_unit
      USE strings_mod,       ONLY : FoundError
#ifdef CHECKPOINTING
      USE wrt_gst_mod,       ONLY : wrt_gst
#endif
      USE wrt_rst_mod,       ONLY : wrt_rst
!
      implicit none
!
      PRIVATE :: IRAM_error
      PUBLIC  :: ROMS_initialize
      PUBLIC  :: ROMS_run
      PUBLIC  :: ROMS_finalize
!
      CONTAINS
!
      SUBROUTINE ROMS_initialize (first, mpiCOMM)
!
!=======================================================================
!                                                                      !
!  This routine allocates and initializes ROMS state variables         !
!  and internal and external parameters.                               !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      logical, intent(inout) :: first
!
      integer, intent(in), optional :: mpiCOMM
!
!  Local variable declarations.
!
      logical :: allocate_vars = .TRUE.
!
#ifdef DISTRIBUTE
      integer :: MyError, MySize
#endif
      integer :: chunk_size, ng, thread
#ifdef _OPENMP
      integer :: my_threadnum
#endif
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_initialize"

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Set distribute-memory (mpi) world communictor.
!-----------------------------------------------------------------------
!
      IF (PRESENT(mpiCOMM)) THEN
        OCN_COMM_WORLD=mpiCOMM
      ELSE
        OCN_COMM_WORLD=MPI_COMM_WORLD
      END IF
      CALL mpi_comm_rank (OCN_COMM_WORLD, MyRank, MyError)
      CALL mpi_comm_size (OCN_COMM_WORLD, MySize, MyError)
#endif
!
!-----------------------------------------------------------------------
!  On first pass, initialize model parameters a variables for all
!  nested/composed grids.  Notice that the logical switch "first"
!  is used to allow multiple calls to this routine during ensemble
!  configurations.
!-----------------------------------------------------------------------
!
      IF (first) THEN
        first=.FALSE.
!
!  Initialize parallel control switches. These scalars switches are
!  independent from standard input parameters.
!
        CALL initialize_parallel
!
!  Set the ROMS standard output unit to write verbose execution info.
!  Notice that the default standard out unit in Fortran is 6.
!
!  In some applications like coupling or disjointed mpi-communications,
!  it is advantageous to write standard output to a specific filename
!  instead of the default Fortran standard output unit 6. If that is
!  the case, it opens such formatted file for writing.
!
        IF (Set_StdOutUnit) THEN
          stdout=stdout_unit(Master)
          Set_StdOutUnit=.FALSE.
        END IF
!
!  Read in model tunable parameters from standard input. Allocate and
!  initialize variables in several modules after the number of nested
!  grids and dimension parameters are known.
!
        CALL inp_par (iTLM)
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
!
!  Set domain decomposition tile partition range.  This range is
!  computed only once since the "first_tile" and "last_tile" values
!  are private for each parallel thread/node.
!
!$OMP PARALLEL
#if defined _OPENMP
      MyThread=my_threadnum()
#elif defined DISTRIBUTE
      MyThread=MyRank
#else
      MyThread=0
#endif
      DO ng=1,Ngrids
        chunk_size=(NtileX(ng)*NtileE(ng)+numthreads-1)/numthreads
        first_tile(ng)=MyThread*chunk_size
        last_tile (ng)=first_tile(ng)+chunk_size-1
      END DO
!$OMP END PARALLEL
!
!  Initialize internal wall clocks. Notice that the timings does not
!  includes processing standard input because several parameters are
!  needed to allocate clock variables.
!
        IF (Master) THEN
          WRITE (stdout,10)
 10       FORMAT (/,' Process Information:',/)
        END IF
!
        DO ng=1,Ngrids
!$OMP PARALLEL
          DO thread=THREAD_RANGE
            CALL wclock_on (ng, iTLM, 0, __LINE__, MyFile)
          END DO
!$OMP END PARALLEL
        END DO
!
!  Allocate and initialize modules variables.
!
!$OMP PARALLEL
        CALL ROMS_allocate_arrays (allocate_vars)
        CALL ROMS_initialize_arrays
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      END IF

#if defined MCT_LIB && (defined ATM_COUPLING || defined WAV_COUPLING)
!
!-----------------------------------------------------------------------
!  Initialize coupling streams between model(s).
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
# ifdef ATM_COUPLING
        CALL initialize_ocn2atm_coupling (ng, MyRank)
# endif
# ifdef WAV_COUPLING
        CALL initialize_ocn2wav_coupling (ng, MyRank)
# endif
      END DO
#endif
!
!-----------------------------------------------------------------------
!  Initialize tangent linear for all grids first in order to compute
!  the size of the state vector, Nstate.  This size is computed in
!  routine "wpoints".
!-----------------------------------------------------------------------

#ifdef FORWARD_FLUXES
!
!  Set the BLK structure to contain the nonlinear model surface fluxes
!  needed by the tangent linear and adjoint models. Also, set switches
!  to process that structure in routine "check_multifile". Notice that
!  it is possible to split the solution into multiple NetCDF files to
!  reduce their size.
!
!  The switch LreadFRC is deactivated because all the atmospheric
!  forcing, including shortwave radiation, is read from the NLM
!  surface fluxes or is assigned during ESM coupling.  Such fluxes
!  are available from the QCK structure. There is no need for reading
!  and processing from the FRC structure input forcing-files.
!
      CALL edit_multifile ('QCK2BLK')
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      DO ng=1,Ngrids
        LreadBLK(ng)=.TRUE.
        LreadFRC(ng)=.FALSE.
      END DO
#endif
!
!  Initialize perturbation tangent linear model.
!
      DO ng=1,Ngrids
        LreadFWD(ng)=.TRUE.
!$OMP PARALLEL
        CALL tl_initial (ng)
!$OMP END PARALLEL
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
!
!
!  Allocate arrays associated with Generalized Stability Theory (GST)
!  analysis.
!
      CALL allocate_storage
!
!  Initialize various IO flags.
!
      Nrun=0
      DO ng=1,Ngrids
        LdefTLM(ng)=.TRUE.
        LwrtTLM(ng)=.TRUE.
        LwrtPER(ng)=.FALSE.
        LcycleTLM(ng)=.FALSE.
        nTLM(ng)=ntimes(ng)
      END DO
!
!  Initialize ARPACK parameters.
!
      Lrvec=.TRUE.                ! Compute Ritz vectors
      bmat='I'                    ! standard eigenvalue problem
      which='LM'                  ! compute NEV largest eigenvalues
      howmany='A'                 ! compute NEV Ritz vectors
      DO ng=1,Ngrids
        ido(ng)=0                 ! reverse communication flag
        info(ng)=0                ! random initial residual vector
        iparam(1,ng)=1            ! exact shifts
        iparam(3,ng)=MaxIterGST   ! maximum number of Arnoldi iterations
        iparam(4,ng)=1            ! block size in the recurrence
        iparam(7,ng)=1            ! type of eigenproblem being solved
      END DO
!
!  ARPACK debugging parameters.
!
      logfil=stdout               ! output logical unit
      ndigit=-3                   ! number of decimal digits
      msaupd=1                    ! iterations, timings, Ritz
      msaup2=1                    ! norms, Ritz values
      msaitr=0
      mseigt=0
      msapps=0
      msgets=0
      mseupd=0
!
!  Determine size of the eigenproblem (Nsize) and size of work space
!  array SworkL (LworkL).
!
      DO ng=1,Ngrids
        Nconv(ng)=0
        Nsize(ng)=Nend(ng)-Nstr(ng)+1
      END DO

#ifdef CHECKPOINTING
!
!  If restart, read in check pointing data GST restart NetCDF file.
!  Otherwise, create check pointing restart NetCDF file.
!
      DO ng=1,Ngrids
        IF (LrstGST) THEN
          CALL get_gst (ng, iTLM)
          ido(ng)=-2
        ELSE
          CALL def_gst (ng, iTLM)
        END IF
        IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END DO
#endif
!
      RETURN
      END SUBROUTINE ROMS_initialize
!
      SUBROUTINE ROMS_run (RunInterval)
!
!=======================================================================
!                                                                      !
!  This routine computes the eigenvectors of the propagator R(0,t)     !
!  for either autonomous  or  non-autonomous ocean circulations. A     !
!  single integration of an arbitrary perturbation state vector "u"    !
!  forward in time over the interval  [0,t]  by the tangent linear     !
!  model: R(0,t)*u.  The eigenspectrum of  R(0,t) is computed with     !
!  the Arnoldi algorithm using ARPACK library:                         !
!                                                                      !
!  Lehoucq, R.B., D.C. Sorensen, and C. Yang, 1997:  ARPACK user's     !
!    guide:  solution  of  large  scale  eigenvalue  problems with     !
!    implicit restarted Arnoldi Methods, Rice University, 140p.        !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations
!
      real(dp), intent(in) :: RunInterval            ! seconds
!
!  Local variable declarations.
!
      logical :: ITERATE, Lcomplex
#ifdef CHECKPOINTING
      logical :: LwrtGST
#endif
!
      integer :: Fcount, Is, Ie, i, icount, iter, ng, srec
      integer :: NconvRitz(Ngrids)
!
      real(r8) :: Enorm

      real(r8), dimension(2) :: my_norm, my_Ivalue, my_Rvalue
!
      TYPE (T_GST), allocatable :: state(:)
      TYPE (T_GST), allocatable :: tl_state(:)
!
      character (len=55) :: string

      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_run"
!
!-----------------------------------------------------------------------
!  Implicit Restarted Arnoldi Method (IRAM) for the computation of
!  optimal perturbation Ritz eigenfunctions.
!-----------------------------------------------------------------------
!
!  Allocate nested grid pointers for state vectors.
!
      IF (.not.allocated(state)) THEN
        allocate ( state(Ngrids) )
      END IF
      IF (.not.allocated(tl_state)) THEN
        allocate ( tl_state(Ngrids) )
      END IF
!
!  Iterate until either convergence or maximum iterations has been
!  exceeded.
!
      iter=0
      ITERATE=.TRUE.
#ifdef CHECKPOINTING
      LwrtGST=.TRUE.
#endif
!
      ITER_LOOP : DO WHILE (ITERATE)
        iter=iter+1
!
!  Reverse communication interface.
!
        DO ng=1,Ngrids
#ifdef PROFILE
          CALL wclock_on (ng, iTLM, 38, __LINE__, MyFile)
#endif
#ifdef DISTRIBUTE
          CALL pdnaupd (OCN_COMM_WORLD,                                 &
     &                  ido(ng), bmat, Nsize(ng), which, NEV,           &
     &                  Ritz_tol,                                       &
     &                  STORAGE(ng)%resid(Nstr(ng)), NCV,               &
     &                  STORAGE(ng)%Bvec(Nstr(ng),1), Nsize(ng),        &
     &                  iparam(1,ng), ipntr(1,ng),                      &
     &                  STORAGE(ng)%SworkD,                             &
     &                  SworkL(1,ng), LworkL, info(ng))
#else
          CALL dnaupd (ido(ng), bmat, Nsize(ng), which, NEV,            &
     &                 Ritz_tol,                                        &
     &                 STORAGE(ng)%resid, NCV,                          &
     &                 STORAGE(ng)%Bvec, Nsize(ng),                     &
     &                 iparam(1,ng), ipntr(1,ng),                       &
     &                 STORAGE(ng)%SworkD,                              &
     &                 SworkL(1,ng), LworkL, info(ng))
#endif
          Nconv(ng)=iaup2(4)
#ifdef PROFILE
          CALL wclock_off (ng, iTLM, 38, __LINE__, MyFile)
#endif
#ifdef CHECKPOINTING
!
!  If appropriate, write out check point data into GST restart NetCDF
!  file. Notice that the restart data is always saved if MaxIterGST
!  is reached without convergence. It is also saved when convergence
!  is achieved (ido=99).
!
          IF ((MOD(iter,nGST).eq.0).or.(iter.ge.MaxIterGST).or.         &
     &        (ANY(ido.eq.99))) THEN
            CALL wrt_gst (ng, iTLM)
            IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
          END IF
#endif
        END DO
!
!  Terminate computations if maximum number of iterations is reached.
!  This will faciliate splitting the analysis in several computational
!  cycles using the restart option.
!
        IF ((iter.ge.MaxIterGST).and.ANY(ido.ne.99)) THEN
          ITERATE=.FALSE.
          EXIT ITER_LOOP
        END IF
!
!  Perform matrix-vector operation:  R`(t,0)u
!
        IF (ANY(ABS(ido).eq.1)) THEN
          DO ng=1,Ngrids
            Fcount=TLM(ng)%load
            TLM(ng)%Nrec(Fcount)=0
            TLM(ng)%Rindex=0
          END DO
!
!  Set state vectors to process by the propagator via pointer
!  equivalence.
!
          DO ng=1,Ngrids
            IF (ASSOCIATED(state(ng)%vector)) THEN
              nullify (state(ng)%vector)
            END IF
            Is=ipntr(1,ng)
            Ie=Is+Nsize(ng)-1
            state(ng)%vector => STORAGE(ng)%SworkD(Is:Ie)

            IF (ASSOCIATED(tl_state(ng)%vector)) THEN
              nullify (tl_state(ng)%vector)
            END IF
            Is=ipntr(2,ng)
            Ie=Is+Nsize(ng)-1
            tl_state(ng)%vector => STORAGE(ng)%SworkD(Is:Ie)
          END DO

!$OMP PARALLEL
          CALL propagator_fte (RunInterval, state, tl_state)
!$OMP END PARALLEL
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
        ELSE
          IF (ANY(info.ne.0)) THEN
            DO ng=1,Ngrids
              IF (info(ng).ne.0) THEN
                IF (Master) THEN
                  CALL IRAM_error (info(ng), 1, string)
                  WRITE (stdout,10) 'DNAUPD', TRIM(string),             &
     &                              ', info = ', info(ng)
                END IF
                RETURN
              END IF
            END DO
          ELSE
!
!  Compute Ritz vectors (the only choice left is IDO=99).  They are
!  generated in ARPACK in decreasing magnitude of its eigenvalue.
!  The most significant is first.
!
            DO ng=1,Ngrids
              NconvRitz(ng)=iparam(5,ng)
              IF (Master) THEN
                WRITE (stdout,20) 'Number of converged Ritz values:',   &
     &                            iparam(5,ng)
                WRITE (stdout,20) 'Number of Arnoldi iterations:',      &
     &                            iparam(3,ng)
              END IF
#ifdef PROFILE
              CALL wclock_on (ng, iTLM, 38, __LINE__, MyFile)
#endif
#ifdef DISTRIBUTE
              CALL pdneupd (OCN_COMM_WORLD,                             &
     &                      Lrvec, howmany, pick(1,ng),                 &
     &                      RvalueR(1,ng), RvalueI(1,ng),               &
     &                      STORAGE(ng)%Rvector(Nstr(ng),1),            &
     &                      Nsize(ng), sigmaR, sigmaI,                  &
     &                      SworkEV(1,ng), bmat, Nsize(ng),             &
     &                      which, NEV, Ritz_tol,                       &
     &                      STORAGE(ng)%resid(Nstr(ng)), NCV,           &
     &                      STORAGE(ng)%Bvec(Nstr(ng),1), Nsize(ng),    &
     &                      iparam(1,ng), ipntr(1,ng),                  &
     &                      STORAGE(ng)%SworkD,                         &
     &                      SworkL(:,ng), LworkL, info(ng))
#else
              CALL dneupd (Lrvec, howmany, pick(1,ng),                  &
     &                     RvalueR(1,ng), RvalueI(1,ng),                &
     &                     STORAGE(ng)%Rvector, Nsize(ng),              &
     &                     sigmaR, sigmaI,                              &
     &                     SworkEV(1,ng), bmat, Nsize(ng),              &
     &                     which, NEV, Ritz_tol,                        &
     &                     STORAGE(ng)%resid, NCV,                      &
     &                     STORAGE(ng)%Bvec, Nsize(ng),                 &
     &                     iparam(1,ng), ipntr(1,ng),                   &
     &                     STORAGE(ng) % SworkD,                        &
     &                     SworkL(1,ng), LworkL, info(ng))
#endif
#ifdef PROFILE
              CALL wclock_off (ng, iTLM, 38, __LINE__, MyFile)
#endif
            END DO

            IF (ANY(info.ne.0)) THEN
              DO ng=1,Ngrids
                IF (info(ng).ne.0) THEN
                  IF (Master) THEN
                    CALL IRAM_error (info(ng), 2, string)
                    WRITE (stdout,10) 'DNEUPD', TRIM(string),           &
     &                                ', info = ', info(ng)
                  END IF
                  RETURN
                END IF
              END DO
            ELSE
!
!  Activate writing of each eigenvector into the tangent history NetCDF
!  file. The "ocean_time" is the eigenvector number. If writing one
!  eigenvector per file (LmultiGST=.TRUE.), both the real and imaginaryg
!  eigenvectors are stored in the same file. Then, this file will
!  contains four records for the initial and final perturbations of the
!  complex eigenvector.
!
              Nrun=0
              icount=0
              Lcomplex=.TRUE.

              DO i=1,MAXVAL(NconvRitz)
                DO ng=1,Ngrids
                  IF ((i.eq.1).or.LmultiGST) THEN
                    Fcount=TLM(ng)%load
                    TLM(ng)%Nrec(Fcount)=0
                    TLM(ng)%Rindex=0
                  END IF
                  IF (LmultiGST.and.Lcomplex) THEN
                    LdefTLM(ng)=.TRUE.
                    icount=icount+1
                    WRITE (TLM(ng)%name,30) TRIM(TLM(ng)%base), icount
                  END IF
                END DO
!
!  Compute and write Ritz eigenvectors.
!
                IF (ANY(RvalueI(i,:).eq.0.0_r8)) THEN
!
!  Ritz value is real.
!
                  DO ng=1,Ngrids
                    Is=Nstr(ng)
                    Ie=Nend(ng)
                    IF (ASSOCIATED(state(ng)%vector)) THEN
                      nullify (state(ng)%vector)
                    END IF

                    IF (ASSOCIATED(tl_state(ng)%vector)) THEN
                      nullify (tl_state(ng)%vector)
                    END IF
                    state(ng)%vector => STORAGE(ng)%Rvector(Is:Ie,i)
                    tl_state(ng)%vector => SworkR(Is:Ie)
                  END DO

!$OMP PARALLEL
                  CALL propagator_fte (RunInterval, state, tl_state)
!$OMP END PARALLEL
                  IF (FoundError(exit_flag, NoError,                    &
     &                           __LINE__, MyFile)) RETURN
!
                  DO ng=1,Ngrids
                    CALL r_norm2 (ng, iTLM, Nstr(ng), Nend(ng),         &
     &                            -RvalueR(i,ng),                       &
     &                            state(ng)%vector,                     &
     &                            tl_state(ng)%vector, Enorm)
                    norm(i,ng)=Enorm
                  END DO

                ELSE IF (Lcomplex) THEN
!
!  Ritz value is complex.
!
                  DO ng=1,Ngrids
                    Is=Nstr(ng)
                    Ie=Nend(ng)
                    IF (ASSOCIATED(state(ng)%vector)) THEN
                      nullify (state(ng)%vector)
                    END IF

                    IF (ASSOCIATED(tl_state(ng)%vector)) THEN
                      nullify (tl_state(ng)%vector)
                    END IF
                    state(ng)%vector => STORAGE(ng)%Rvector(Is:Ie,i)
                    tl_state(ng)%vector => SworkR(Is:Ie)
                  END DO

!$OMP PARALLEL
                  CALL propagator_fte (RunInterval, state, tl_state)
!$OMP END PARALLEL
                  IF (FoundError(exit_flag, NoError,                    &
     &                           __LINE__, MyFile)) RETURN
!
                  DO ng=1,Ngrids
                    CALL c_norm2 (ng, iTLM, Nstr(ng), Nend(ng),         &
     &                            -RvalueR(i,ng), RvalueI(i,ng),        &
     &                            STORAGE(ng)%Rvector(Nstr(ng):,i  ),   &
     &                            STORAGE(ng)%Rvector(Nstr(ng):,i+1),   &
     &                            tl_state(ng)%vector, Enorm)
                    norm(i,ng)=Enorm
                  END DO
!
                  DO ng=1,Ngrids
                    Is=Nstr(ng)
                    Ie=Nend(ng)
                    IF (ASSOCIATED(state(ng)%vector)) THEN
                      nullify (state(ng)%vector)
                    END IF

                    IF (ASSOCIATED(tl_state(ng)%vector)) THEN
                      nullify (tl_state(ng)%vector)
                    END IF
                    state(ng)%vector => STORAGE(ng)%Rvector(Is:Ie,i+1)
                    tl_state(ng)%vector => SworkR(Is:Ie)
                  END DO

!$OMP PARALLEL
                  CALL propagator_fte (RunInterval, state, tl_state)
!$OMP END PARALLEL
                  IF (FoundError(exit_flag, NoError,                    &
     &                           __LINE__, MyFile)) RETURN
!
                  DO ng=1,Ngrids
                    CALL c_norm2 (ng, iTLM, Nstr(ng), Nend(ng),         &
     &                            -RvalueR(i,ng), -RvalueI(i,ng),       &
     &                            STORAGE(ng)%Rvector(Nstr(ng):,i+1),   &
     &                            STORAGE(ng)%Rvector(Nstr(ng):,i  ),   &
     &                            tl_state(ng)%vector, Enorm)
                    norm(i  ,ng)=SQRT(norm(i,ng)*norm(i,ng)+            &
     &                                Enorm*Enorm)
                    norm(i+1,ng)=norm(i,ng)
                    Lcomplex=.FALSE.
                  END DO
                ELSE
                  Lcomplex=.TRUE.
                END IF
                IF (Master) THEN
                  DO ng=1,Ngrids
                    WRITE (stdout,40) i, norm(i,ng),                    &
     &                                RvalueR(i,ng), RvalueI(i,ng), i
                  END DO
                END IF
!
!  Write out Ritz eigenvalues and Ritz eigenvector Euclidean norm
!  (residual) to NetCDF file(s).  Notice that we write the same value
!  twice in the TLM file for the initial and final perturbation of
!  the eigenvector.
!
                SourceFile=MyFile
                DO ng=1,Ngrids
                  my_norm(1)=norm(i,ng)
                  my_norm(2)=my_norm(1)
                  my_Rvalue(1)=RvalueR(i,ng)
                  my_Rvalue(2)=my_Rvalue(1)
                  my_Ivalue(1)=RvalueI(i,ng)
                  my_Ivalue(2)=my_Ivalue(1)
                  IF (LmultiGST) THEN
                    IF (.not.Lcomplex.or.(RvalueI(i,ng).eq.0.0_r8)) THEN
                      srec=1
                    ELSE
                      srec=3
                    END IF
                  ELSE
                    srec=2*i-1
                  END IF
!
                  IF (LwrtTLM(ng)) THEN
                    SELECT CASE (TLM(ng)%IOtype)
                      CASE (io_nf90)
                        CALL netcdf_put_fvar (ng, iTLM,                 &
     &                                        TLM(ng)%name,             &
     &                                        'Ritz_rvalue',            &
     &                                        my_Rvalue,                &
     &                                        start = (/srec/),         &
     &                                        total = (/2/),            &
     &                                        ncid = TLM(ng)%ncid)

                        IF (FoundError(exit_flag, NoError,              &
     &                                 __LINE__, MyFile)) RETURN
!
                        CALL netcdf_put_fvar (ng, iTLM,                 &
     &                                        TLM(ng)%name,             &
     &                                        'Ritz_ivalue',            &
     &                                        my_Ivalue,                &
     &                                        start = (/srec/),         &
     &                                        total = (/2/),            &
     &                                        ncid = TLM(ng)%ncid)

                        IF (FoundError(exit_flag, NoError,              &
     &                                 __LINE__, MyFile)) RETURN
!
                        CALL netcdf_put_fvar (ng, iTLM,                 &
     &                                        TLM(ng)%name,             &
     &                                        'Ritz_norm',              &
     &                                        my_norm,                  &
     &                                        start = (/srec/),         &
     &                                        total = (/2/),            &
     &                                        ncid = TLM(ng)%ncid)

                        IF (FoundError(exit_flag, NoError,              &
     &                                 __LINE__, MyFile)) RETURN

#if defined PIO_LIB && defined DISTRIBUTE
                      CASE (io_pio)
                        CALL pio_netcdf_put_fvar (ng, iTLM,             &
     &                                            TLM(ng)%name,         &
     &                                            'Ritz_rvalue',        &
     &                                            my_Rvalue,            &
     &                                            start = (/srec/),     &
     &                                            total = (/2/),        &
     &                                        pioFile = TLM(ng)%pioFile)

                        IF (FoundError(exit_flag, NoError,              &
     &                                 __LINE__, MyFile)) RETURN
!
                        CALL pio_netcdf_put_fvar (ng, iTLM,             &
     &                                            TLM(ng)%name,         &
     &                                            'Ritz_ivalue',        &
     &                                            my_Ivalue,            &
     &                                            start = (/srec/),     &
     &                                            total = (/2/),        &
     &                                        pioFile = TLM(ng)%pioFile)

                        IF (FoundError(exit_flag, NoError,              &
     &                                 __LINE__, MyFile)) RETURN
!
                        CALL pio_netcdf_put_fvar (ng, iTLM,             &
     &                                            TLM(ng)%name,         &
     &                                            'Ritz_norm',          &
     &                                            my_norm,              &
     &                                            start = (/srec/),     &
     &                                            total = (/2/),        &
     &                                        pioFile = TLM(ng)%pioFile)

                        IF (FoundError(exit_flag, NoError,              &
     &                                 __LINE__, MyFile)) RETURN
#endif
                    END SELECT
!
                    IF (LmultiGST.and.Lcomplex) THEN
                      CALL close_file (ng, iTLM, TLM(ng), TLM(ng)%name)
                      IF (FoundError(exit_flag, NoError,                &
     &                               __LINE__, MyFile)) RETURN
                    END IF
                  END IF
                END DO
              END DO
            END IF
          END IF
          ITERATE=.FALSE.
        END IF

      END DO ITER_LOOP
!
 10   FORMAT (/,1x,'Error in ',a,1x,a,a,1x,i5,/)
 20   FORMAT (/,a,1x,i2,/)
 30   FORMAT (a,'_',i3.3,'.nc')
 40   FORMAT (1x,i4.4,'-th residual',1p,e14.6,0p,                       &
     &        '  Ritz values',1pe14.6,0p,1x,1pe14.6,2x,i4.4)
!
      RETURN
      END SUBROUTINE ROMS_run
!
      SUBROUTINE ROMS_finalize
!
!=======================================================================
!                                                                      !
!  This routine terminates ROMS nonlinear and adjoint models           !
!  execution.                                                          !
!                                                                      !
!=======================================================================
!
!  Local variable declarations.
!
      integer :: Fcount, ng, thread
!
      character (len=*), parameter :: MyFile =                          &
     &  __FILE__//", ROMS_finalize"
!
!-----------------------------------------------------------------------
!  If blowing-up, save latest model state into RESTART NetCDF file.
!-----------------------------------------------------------------------
!
!  If cycling restart records, write solution into the next record.
!
      IF (exit_flag.eq.1) THEN
        DO ng=1,Ngrids
          IF (LwrtRST(ng)) THEN
            IF (Master) WRITE (stdout,10)
 10         FORMAT (/,' Blowing-up: Saving latest model state into ',   &
     &                ' RESTART file',/)
            Fcount=RST(ng)%load
            IF (LcycleRST(ng).and.(RST(ng)%Nrec(Fcount).ge.2)) THEN
              RST(ng)%Rindex=2
              LcycleRST(ng)=.FALSE.
            END IF
            blowup=exit_flag
            exit_flag=NoError
#ifdef DISTRIBUTE
            CALL wrt_rst (ng, MyRank)
#else
            CALL wrt_rst (ng, -1)
#endif
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Stop model and time profiling clocks, report memory requirements, and
!  close output NetCDF files.
!-----------------------------------------------------------------------
!
!  Stop time clocks.
!
      IF (Master) THEN
        WRITE (stdout,20)
 20     FORMAT (/,'Elapsed wall CPU time for each process (seconds):',/)
      END IF
!
      DO ng=1,Ngrids
!$OMP PARALLEL
        DO thread=THREAD_RANGE
          CALL wclock_off (ng, iTLM, 0, __LINE__, MyFile)
        END DO
!$OMP END PARALLEL
      END DO
!
!  Report dynamic memory and automatic memory requirements.
!
!$OMP PARALLEL
      CALL memory
!$OMP END PARALLEL
!
!  Close IO files.
!
      DO ng=1,Ngrids
        CALL close_inp (ng, iTLM)
      END DO
      CALL close_out
!
      RETURN
      END SUBROUTINE ROMS_finalize
!
      SUBROUTINE IRAM_error (info, icall, string)
!
!=======================================================================
!                                                                      !
!  This routine decodes internal error messages from the Implicit      !
!  Restarted Arnoldi Method (IRAM) for the computation of optimal      !
!  perturbation Ritz eigenfunctions.                                   !
!                                                                      !
!=======================================================================
!
!  imported variable declarations.
!
      integer, intent(in) :: info, icall
!
      character (len=*), intent(out) :: string
!
!-----------------------------------------------------------------------
!  Decode error message from IRAM.
!-----------------------------------------------------------------------
!
      IF (info.eq.0)  THEN
        string='Normal exit                                            '
      ELSE IF (info.eq.1) THEN
        IF (icall.eq.1) THEN
          string='Maximum number of iterations taken                   '
        ELSE
          string='Could not reorder Schur vectors                      '
        END IF
      ELSE IF (info.eq.3) THEN
        string='No shifts could be applied during an IRAM cycle        '
      ELSE IF (info.eq.-1) THEN
        string='Nstate must be positive                                '
      ELSE IF (info.eq.-2) THEN
        string='NEV must be positive                                   '
      ELSE IF (info.eq.-3) THEN
        string='NCV must be greater NEV and less than or equal Nstate  '
      ELSE IF (info.eq.-4) THEN
        string='Maximum number of iterations must be greater than zero '
      ELSE IF (info.eq.-5) THEN
        string='WHICH must be one of LM, SM, LA, SA or BE              '
      ELSE IF (info.eq.-6) THEN
        string='BMAT must be one of I or G                             '
      ELSE IF (info.eq.-7) THEN
        string='Length of private work array SworkL is not sufficient  '
      ELSE IF (info.eq.-8) THEN
        IF (icall.eq.1) THEN
          string='Error return from LAPACK eigenvalue calculation      '
        ELSE
          string='Error in DLAHQR in the Shurn vectors calculation     '
        END IF
      ELSE IF (info.eq.-9) THEN
        IF (icall.eq.1) THEN
          string='Starting vector is zero'
        ELSE
          string='Error in DTREVC in the eigenvectors calculation      '
        END IF
      ELSE IF (info.eq.-10) THEN
        string='IPARAM(7) must be 1, 2, 3, 4, 5                        '
      ELSE IF (info.eq.-11) THEN
        string='IPARAM(7) = 1 and BMAT = G are incompatable            '
      ELSE IF (info.eq.-12) THEN
        IF (icall.eq.1) THEN
          string='IPARAM(1) must be equal to 0 or 1                    '
        ELSE
          string='HOWMANY = S not yet implemented                      '
        END IF
      ELSE IF (info.eq.-13) THEN
        string='HOWMANY must be one of A or P if Lrvec = .TRUE.        '
      ELSE IF (info.eq.-14) THEN
        string='Did not find any eigenvalues to sufficient accuaracy   '
      ELSE IF (info.eq.-15) THEN
        string='Different count of converge Ritz values in DNEUPD      '
      ELSE IF (info.eq.-9999) THEN
        string='Could not build and Arnoldi factorization              '
      END IF
!
      RETURN
      END SUBROUTINE IRAM_error

      END MODULE roms_kernel_mod
