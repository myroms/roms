#include "cppdefs.h"
      SUBROUTINE inp_par
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine reads in input model parameters from standard input.   !
!  It also writes out these parameters to standard output.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Local variable declarations.
!
      logical :: Lwrite

      integer :: inp, out, itrc, ng

      real(r8) :: cff

      character (len=80) :: Iname
!
!-----------------------------------------------------------------------
!  Read in and report input model parameters.
!-----------------------------------------------------------------------
!
!  Set input units.
!
#ifdef DISTRIBUTE
      Lwrite=Master
      inp=1
      out=stdout
#else
      Lwrite=Master
      inp=stdinp
      out=stdout
#endif
!
!  Get current data.
!
#ifdef DISTRIBUTE
      IF (Master) CALL get_date (date_str)
      CALL mp_bcasts (ng, date_str, 44)
#else
      CALL get_date (date_str)
#endif
!
!-----------------------------------------------------------------------
!  Read in physical model input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) WRITE (out,10) version, TRIM(date_str)
 10   FORMAT (/,' Model Input Parameters:  ROMS/TOMS version ',a,/,     &
     &       26x,a,/,1x,77('-'))
#ifdef DISTRIBUTE
!!    CALL my_getarg (1,Iname)
      IF (Master) CALL my_getarg (1,Iname)
      CALL mp_bcasts (ng, Iname, 80)
      OPEN (inp, FILE=TRIM(Iname), FORM='formatted', STATUS='old',      &
     &      ERR=20)
      GOTO 40
 20   WRITE (stdout,30)
      STOP
 30   FORMAT (/,' INP_PAR - Unable to open ROMS/TOMS input script ',    &
     &              'file.',                                            &
     &        /,11x,'In distributed-memory applications, the input',    &
     &        /,11x,'script file is processed in parallel. The Unix',   &
     &        /,11x,'routine GETARG is used to get script file name.',  &
     &        /,11x,'For example, in MPI applications make sure that',  &
     &        /,11x,'command line is something like:',/,                &
     &        /,11x,'mpirun -np 4 ocean ocean.in',/,/,11x,'and not',/,  &
     &        /,11x,'mpirun -np 4 ocean < ocean.in',/)
 40   CONTINUE
#endif
      CALL read_PhyPar (inp, out, Lwrite)

#ifdef BIOLOGY
!
!-----------------------------------------------------------------------
!  Read in biological model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (inp, FILE=TRIM(bparnam), FORM='formatted', STATUS='old')

      CALL read_BioPar (inp, out, Lwrite)
#endif
#ifdef SEDIMENT
!
!-----------------------------------------------------------------------
!  Read in sediment model input parameters.
!-----------------------------------------------------------------------
!
      OPEN (inp ,FILE=TRIM(sparnam), FORM='formatted', STATUS='old')

      CALL read_SedPar (inp, out, Lwrite)
#endif
#if defined ASSIMILATION || defined NUDGING
!
!-----------------------------------------------------------------------
!  Read in input assimilation parameters.
!-----------------------------------------------------------------------
!
      OPEN (inp, FILE=TRIM(aparnam), FORM='formatted', STATUS='old')

      CALL read_AssPar (inp, out, Lwrite)
#endif
#ifdef FLOATS
!
!-----------------------------------------------------------------------
!  Read in floats input parameters.
!-----------------------------------------------------------------------
!
      OPEN (inp, FILE=TRIM(fposnam), FORM='formatted', STATUS='old')

      CALL read_FloatsPar (inp, out, Lwrite)
#endif
#ifdef STATIONS
!
!-----------------------------------------------------------------------
!  Read in stations input parameters.
!-----------------------------------------------------------------------
!
      OPEN (inp, FILE=TRIM(sposnam), FORM='formatted', STATUS='old')

      CALL read_StaPar (inp, out, Lwrite)
#endif
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options and definitions.
!-----------------------------------------------------------------------
!
      IF (Master) CALL checkdefs
      IF (Master) CALL my_flush (out)
!
!-----------------------------------------------------------------------
!  Compute various constants.
!-----------------------------------------------------------------------
!
      gorho0=g/rho0
      DO ng=1,Ngrids
        dtfast(ng)=dt(ng)/REAL(ndtfast(ng),r8)
!
!  Take the square root of the biharmonic coefficients so it can
!  be applied to each harmonic operator.
!
        visc4(ng)=SQRT(ABS(visc4(ng)))
        tkenu4(ng)=SQRT(ABS(tkenu4(ng)))
        DO itrc=1,NAT
          tnu4(itrc,ng)=SQRT(ABS(tnu4(itrc,ng)))
        END DO
!
!  Compute inverse nudging coefficients (1/s) used in various tasks.
!
        IF (Znudg(ng).gt.0.0_r8) THEN
          Znudg(ng)=1.0_r8/(Znudg(ng)*86400.0_r8)
        ELSE
          Znudg(ng)=0.0_r8
        END IF
        IF (M2nudg(ng).gt.0.0_r8) THEN
          M2nudg(ng)=1.0_r8/(M2nudg(ng)*86400.0_r8)
        ELSE
          M2nudg(ng)=0.0_r8
        END IF
#ifdef SOLVE3D
        IF (M3nudg(ng).gt.0.0_r8) THEN
          M3nudg(ng)=1.0_r8/(M3nudg(ng)*86400.0_r8)
        ELSE
          M3nudg(ng)=0.0_r8
        END IF
        DO itrc=1,NAT
          IF (Tnudg(itrc,ng).gt.0.0_r8) THEN
            Tnudg(itrc,ng)=1.0_r8/(Tnudg(itrc,ng)*86400.0_r8)
          ELSE
            Tnudg(itrc,ng)=0.0_r8
          END IF
        END DO
#endif
!
!  Set scales for input forcing fields.  Notice that for consitency
!  with Boussinesq approximation, these scales are overwritten here
!  to use the provided value for "rho0".
!
        cff=1.0_r8/rho0
        Fscale(idUsms,ng)=cff
        Fscale(idVsms,ng)=cff
        Fscale(idUbms,ng)=cff
        Fscale(idVbms,ng)=cff
        Fscale(idUbws,ng)=cff
        Fscale(idVbws,ng)=cff
        cff=1.0_r8/(rho0*Cp)
        Fscale(idTsur(itemp),ng)=cff
        Fscale(idTbot(itemp),ng)=cff
        Fscale(idSrad,ng)=cff
        Fscale(idLrad,ng)=cff
        Fscale(idLhea,ng)=cff
        Fscale(idShea,ng)=cff
        Fscale(iddQdT,ng)=cff
      END DO
      RETURN
      END SUBROUTINE inp_par

      SUBROUTINE read_PhyPar (inp, out, Lwrite)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine reads in physical model input parameters.              !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
      USE mod_strings
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer, parameter :: Tunit = 2

      integer :: Lstr, Npts, Nval, i, itrc, k, ng, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(NAT,Ngrids) :: Ltracer

      real(r8), dimension(NAT,Ngrids) :: Rtracer

      real(r8), dimension(100) :: Rval

      character (len=1 ), parameter :: blank = ' '
      character (len=19) :: ref_att
      character (len=40) :: KeyWord
      character (len=80) :: line, fname
      character (len=80), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in physical model parameters. Then, load input data into module.
!  Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'VARNAME') THEN
            DO i=1,LEN(usrname)
              varname(i:i)=blank
            END DO
            varname=TRIM(ADJUSTL(Cval(Nval)))
            CALL initialize_ncparam
          ELSE IF (TRIM(KeyWord).eq.'NtileI') THEN
            Npts=load_i(Nval, Rval, Ngrids, NtileI)
#ifdef DISTRIBUTE
            NtileX(1:Ngrids)=1
#else
            NtileX(1:Ngrids)=NtileI(1:Ngrids)
#endif
          ELSE IF (TRIM(KeyWord).eq.'NtileJ') THEN
            Npts=load_i(Nval, Rval, Ngrids, NtileJ)
#ifdef DISTRIBUTE
            NtileE(1:Ngrids)=1
#else
            NtileE(1:Ngrids)=NtileJ(1:Ngrids)
#endif
          ELSE IF (TRIM(KeyWord).eq.'NTIMES') THEN
            Npts=load_i(Nval, Rval, 1, ntimes)
          ELSE IF (TRIM(KeyWord).eq.'DT') THEN
            Npts=load_r(Nval, Rval, Ngrids, dt)
          ELSE IF (TRIM(KeyWord).eq.'NDTFAST') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndtfast)
          ELSE IF (TRIM(KeyWord).eq.'ERstr') THEN
            Npts=load_i(Nval, Rval, 1, ERstr)
          ELSE IF (TRIM(KeyWord).eq.'ERend') THEN
            Npts=load_i(Nval, Rval, 1, ERend)
          ELSE IF (TRIM(KeyWord).eq.'NRREC') THEN
            Npts=load_i(Nval, Rval, 1, nrrec)
            LastRec=.false.
            IF (nrrec.lt.0) LastRec=.true.
          ELSE IF (TRIM(KeyWord).eq.'LCYCLE') THEN
            Npts=load_l(Nval, Cval, Ngrids, lcycle)
          ELSE IF (TRIM(KeyWord).eq.'NRST') THEN
            Npts=load_i(Nval, Rval, 1, nrst)
          ELSE IF (TRIM(KeyWord).eq.'NSTA') THEN
            Npts=load_i(Nval, Rval, Ngrids, nsta)
          ELSE IF (TRIM(KeyWord).eq.'NFLT') THEN
            Npts=load_i(Nval, Rval, Ngrids, nflt)
          ELSE IF (TRIM(KeyWord).eq.'NINFO') THEN
            Npts=load_i(Nval, Rval, Ngrids, ninfo)
          ELSE IF (TRIM(KeyWord).eq.'LDEFOUT') THEN
            Npts=load_l(Nval, Cval, Ngrids, ldefout)
          ELSE IF (TRIM(KeyWord).eq.'NHIS') THEN
            Npts=load_i(Nval, Rval, Ngrids, nhis)
          ELSE IF (TRIM(KeyWord).eq.'NDEFHIS') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefhis)
          ELSE IF (TRIM(KeyWord).eq.'NTSAVG') THEN
            Npts=load_i(Nval, Rval, Ngrids, ntsavg)
          ELSE IF (TRIM(KeyWord).eq.'NAVG') THEN
            Npts=load_i(Nval, Rval, Ngrids, navg)
          ELSE IF (TRIM(KeyWord).eq.'NDEFAVG') THEN
            Npts=load_i(Nval, Rval, Ngrids, ndefavg)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NAT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                tnu2(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NAT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                tnu4(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'VISC2') THEN
            Npts=load_r(Nval, Rval, Ngrids, visc2)
          ELSE IF (TRIM(KeyWord).eq.'VISC4') THEN
            Npts=load_r(Nval, Rval, Ngrids, visc4)
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NAT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                Akt_bak(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKV_BAK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Akv_bak)
          ELSE IF (TRIM(KeyWord).eq.'AKK_BAK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Akk_bak)
          ELSE IF (TRIM(KeyWord).eq.'AKP_BAK') THEN
            Npts=load_r(Nval, Rval, Ngrids, Akp_bak)
          ELSE IF (TRIM(KeyWord).eq.'TKENU2') THEN
            Npts=load_r(Nval, Rval, Ngrids, tkenu2)
          ELSE IF (TRIM(KeyWord).eq.'TKENU4') THEN
            Npts=load_r(Nval, Rval, Ngrids, tkenu4)
          ELSE IF (TRIM(KeyWord).eq.'GLS_P') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_p)
          ELSE IF (TRIM(KeyWord).eq.'GLS_M') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_m)
          ELSE IF (TRIM(KeyWord).eq.'GLS_N') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_n)
          ELSE IF (TRIM(KeyWord).eq.'GLS_Kmin') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_Kmin)
          ELSE IF (TRIM(KeyWord).eq.'GLS_Pmin') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_Pmin)
          ELSE IF (TRIM(KeyWord).eq.'GLS_CMU0') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_cmu0)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C1') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c1)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C2') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c2)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C3M') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c3m)
          ELSE IF (TRIM(KeyWord).eq.'GLS_C3P') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_c3p)
          ELSE IF (TRIM(KeyWord).eq.'GLS_SIGK') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_sigk)
          ELSE IF (TRIM(KeyWord).eq.'GLS_SIGP') THEN
            Npts=load_r(Nval, Rval, Ngrids, gls_sigp)
          ELSE IF (TRIM(KeyWord).eq.'RDRG') THEN
            Npts=load_r(Nval, Rval, Ngrids, rdrg)
          ELSE IF (TRIM(KeyWord).eq.'RDRG2') THEN
            Npts=load_r(Nval, Rval, Ngrids, rdrg2)
          ELSE IF (TRIM(KeyWord).eq.'Zob') THEN
            Npts=load_r(Nval, Rval, Ngrids, Zob)
          ELSE IF (TRIM(KeyWord).eq.'Zos') THEN
            Npts=load_r(Nval, Rval, Ngrids, Zos)
          ELSE IF (TRIM(KeyWord).eq.'WTYPE') THEN
            Npts=load_i(Nval, Rval, Ngrids, lmd_Jwt)
          ELSE IF (TRIM(KeyWord).eq.'LEVSFRC') THEN
            Npts=load_i(Nval, Rval, Ngrids, levsfrc)
          ELSE IF (TRIM(KeyWord).eq.'LEVBFRC') THEN
            Npts=load_i(Nval, Rval, Ngrids, levbfrc)
          ELSE IF (TRIM(KeyWord).eq.'THETA_S') THEN
            Npts=load_r(Nval, Rval, Ngrids, theta_s)
          ELSE IF (TRIM(KeyWord).eq.'THETA_B') THEN
            Npts=load_r(Nval, Rval, Ngrids, theta_b)
          ELSE IF (TRIM(KeyWord).eq.'TCLINE') THEN
            Npts=load_r(Nval, Rval, Ngrids, Tcline)
          ELSE IF (TRIM(KeyWord).eq.'RHO0') THEN
            Npts=load_r(Nval, Rval, 1, rho0)
          ELSE IF (TRIM(KeyWord).eq.'DSTART') THEN
            Npts=load_r(Nval, Rval, 1, dstart)
          ELSE IF (TRIM(KeyWord).eq.'TIME_REF') THEN
            Npts=load_r(Nval, Rval, 1, time_ref)
            r_text=ref_att(time_ref,r_date)
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NAT*Ngrids, Rtracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                Tnudg(itrc,ng)=Rtracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ZNUDG') THEN
            Npts=load_r(Nval, Rval, Ngrids, Znudg)
          ELSE IF (TRIM(KeyWord).eq.'M2NUDG') THEN
            Npts=load_r(Nval, Rval, Ngrids, M2nudg)
          ELSE IF (TRIM(KeyWord).eq.'M3NUDG') THEN
            Npts=load_r(Nval, Rval, Ngrids, M3nudg)
          ELSE IF (TRIM(KeyWord).eq.'OBCFAC') THEN
            Npts=load_r(Nval, Rval, Ngrids, obcfac)
          ELSE IF (TRIM(KeyWord).eq.'R0') THEN
            Npts=load_r(Nval, Rval, Ngrids, R0)
          ELSE IF (TRIM(KeyWord).eq.'T0') THEN
            Npts=load_r(Nval, Rval, Ngrids, T0)
          ELSE IF (TRIM(KeyWord).eq.'S0') THEN
            Npts=load_r(Nval, Rval, Ngrids, S0)
          ELSE IF (TRIM(KeyWord).eq.'TCOEF') THEN
            Npts=load_r(Nval, Rval, Ngrids, Tcoef)
          ELSE IF (TRIM(KeyWord).eq.'SCOEF') THEN
            Npts=load_r(Nval, Rval, Ngrids, Scoef)
          ELSE IF (TRIM(KeyWord).eq.'GAMMA2') THEN
            Npts=load_r(Nval, Rval, 1, gamma2)
          ELSE IF (TRIM(KeyWord).eq.'Hout(idFsur)') THEN
            IF (idFsur.eq.0) THEN
              WRITE (out,280) 'idFsur'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idFsur,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbar)') THEN
            IF (idUbar.eq.0) THEN
              WRITE (out,280) 'idUbar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbar,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbar)') THEN
            IF (idVbar.eq.0) THEN
              WRITE (out,280) 'idVbar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbar,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUvel)') THEN
            IF (idUvel.eq.0) THEN
              WRITE (out,280) 'idUvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVvel)') THEN
            IF (idVvel.eq.0) THEN
              WRITE (out,280) 'idVvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idWvel)') THEN
            IF (idWvel.eq.0) THEN
              WRITE (out,280) 'idWvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idWvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idOvel)') THEN
            IF (idOvel.eq.0) THEN
              WRITE (out,280) 'idOvel'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idOvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            IF (maxval(idTvar).eq.0) THEN
              WRITE (out,280) 'idTvar'
              STOP
            END IF
            Npts=load_l(Nval, Cval, NAT*Ngrids, Ltracer)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                i=idTvar(itrc)
                Hout(i,ng)=Ltracer(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUsms)') THEN
            IF (idUsms.eq.0) THEN
              WRITE (out,280) 'idUsms'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUsms,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVsms)') THEN
            IF (idVsms.eq.0) THEN
              WRITE (out,280) 'idVsms'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVsms,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbms)') THEN
            IF (idUbms.eq.0) THEN
              WRITE (out,280) 'idUbms'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbms,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbms)') THEN
            IF (idVbms.eq.0) THEN
              WRITE (out,280) 'idVbms'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbms,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbws)') THEN
            IF (idUbws.eq.0) THEN
              WRITE (out,280) 'idUbws'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbws,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbws)') THEN
            IF (idVbws.eq.0) THEN
              WRITE (out,280) 'idVbws'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbws,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idAbed)') THEN
            IF (idAbed.eq.0) THEN
              WRITE (out,280) 'idAbed'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idAbed,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbed)') THEN
            IF (idUbed.eq.0) THEN
              WRITE (out,280) 'idUbed'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbed,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbed)') THEN
            IF (idVbed.eq.0) THEN
              WRITE (out,280) 'idVbed'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbed,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idUbot)') THEN
            IF (idUbot.eq.0) THEN
              WRITE (out,280) 'idUbot'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idUbot,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVbot)') THEN
            IF (idVbot.eq.0) THEN
              WRITE (out,280) 'idVbot'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVbot,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHrip)') THEN
            IF (idHrip.eq.0) THEN
              WRITE (out,280) 'idHrip'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHrip,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idLrip)') THEN
            IF (idLrip.eq.0) THEN
              WRITE (out,280) 'idLrip'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idLrip,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idZnot)') THEN
            IF (idZnot.eq.0) THEN
              WRITE (out,280) 'idZnot'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idZnot,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idZapp)') THEN
            IF (idZapp.eq.0) THEN
              WRITE (out,280) 'idZapp'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idZapp,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTsur)') THEN
            IF (idTsur(itemp).eq.0) THEN
              WRITE (out,280) 'idTsur'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idTsur(itemp),1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idLhea)') THEN
            IF (idLhea.eq.0) THEN
              WRITE (out,280) 'idLhea'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idLhea,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idShea)') THEN
            IF (idShea.eq.0) THEN
              WRITE (out,280) 'idShea'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idShea,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idLrad)') THEN
            IF (idLrad.eq.0) THEN
              WRITE (out,280) 'idLrad'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idLrad,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSrad)') THEN
            IF (idSrad.eq.0) THEN
              WRITE (out,280) 'idSrad'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSrad,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idDano)') THEN
            IF (idDano.eq.0) THEN
              WRITE (out,280) 'idDano'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idDano,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idVvis)') THEN
            IF (idVvis.eq.0) THEN
              WRITE (out,280) 'idVvis'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idVvis,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTdif)') THEN
            IF (idTdif.eq.0) THEN
              WRITE (out,280) 'idTdif'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idTdif,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idSdif)') THEN
            IF (idSdif.eq.0) THEN
              WRITE (out,280) 'idSdif'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idSdif,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHsbl)') THEN
            IF (idHsbl.eq.0) THEN
              WRITE (out,280) 'idHsbl'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHsbl,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idHbbl)') THEN
            IF (idHbbl.eq.0) THEN
              WRITE (out,280) 'idHbbl'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idHbbl,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idMtke)') THEN
            IF (idMtke.eq.0) THEN
              WRITE (out,280) 'idMtke'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idMtke,1))
          ELSE IF (TRIM(KeyWord).eq.'Hout(idMtls)') THEN
            IF (idMtls.eq.0) THEN
              WRITE (out,280) 'idMtls'
              STOP
            END IF
            Npts=load_l(Nval, Cval, Ngrids, Hout(idMtls,1))
          ELSE IF (TRIM(KeyWord).eq.'NUSER') THEN
            Npts=load_i(Nval, Rval, 1, Nuser)
          ELSE IF (TRIM(KeyWord).eq.'USER') THEN
            Npts=load_r(Nval, Rval, MAX(1,Nuser), user)
          ELSE IF (TRIM(KeyWord).eq.'TITLE') THEN
            IF (Nval.eq.1) THEN
              title=TRIM(ADJUSTL(Cval(Nval)))
            ELSE
              title=TRIM(ADJUSTL(title))//' '//TRIM(ADJUSTL(Cval(Nval)))
            END IF           
          ELSE IF (TRIM(KeyWord).eq.'RSTNAME') THEN
            DO i=1,LEN(rstname(Nval))
              rstname(Nval)(i:i)=blank
            END DO
            rstname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'HISNAME') THEN
            DO i=1,LEN(hisname(Nval))
              hisname(Nval)(i:i)=blank
              hisbase(Nval)(i:i)=blank
            END DO
            hisname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            hisbase(Nval)=TRIM(ADJUSTL(hisname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'AVGNAME') THEN
            DO i=1,LEN(avgname(Nval))
              avgname(Nval)(i:i)=blank
              avgbase(Nval)(i:i)=blank
            END DO
            avgname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
            avgbase(Nval)=TRIM(ADJUSTL(avgname(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'STANAME') THEN
            DO i=1,LEN(staname(Nval))
              staname(Nval)(i:i)=blank
            END DO
            staname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'FLTNAME') THEN
            DO i=1,LEN(fltname(Nval))
              fltname(Nval)(i:i)=blank
            END DO
            fltname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'GRDNAME') THEN
            DO i=1,LEN(grdname(Nval))
              grdname(Nval)(i:i)=blank
            END DO
            grdname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'ININAME') THEN
            DO i=1,LEN(ininame(Nval))
              ininame(Nval)(i:i)=blank
            END DO
            ininame(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'NFFILES') THEN
            Npts=load_i(Nval, Rval, Ngrids, nFfiles)
            DO ng=1,Ngrids
              IF (nFfiles(ng).le.0) THEN
                WRITE (out,260) 'NFFILES',                              &
     &                          'Must be equal or greater than one.'
                STOP
              END IF
            END DO                
            Npts=MAXVAL(nFfiles)           
            allocate ( frcids (Npts,Ngrids) )
            allocate ( frcname(Npts,Ngrids) )
            frcids(1:Npts,1:Ngrids)=-1
            DO ng=1,Ngrids
              DO k=1,Npts
                DO i=1,LEN(frcname(k,ng))
                  frcname(k,ng)(i:i)=blank
                END DO
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'FRCNAME') THEN
            DO ng=1,Ngrids
              DO i=1,nFfiles(ng)
                IF (Nval.eq.(i+(Ngrids-1)*nFfiles(ng))) THEN
                  frcname(i,ng)=TRIM(ADJUSTL(Cval(Nval)))
                END IF
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'CLMNAME') THEN
            DO i=1,LEN(clmname(Nval))
              clmname(Nval)(i:i)=blank
            END DO
            clmname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'BRYNAME') THEN
            DO i=1,LEN(bryname(Nval))
              bryname(Nval)(i:i)=blank
            END DO
            bryname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'FWDNAME') THEN
            DO i=1,LEN(fwdname(Nval))
              fwdname(Nval)(i:i)=blank
            END DO
            fwdname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'APARNAM') THEN
            DO i=1,LEN(aparnam)
              aparnam(i:i)=blank
            END DO
            aparnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'SPOSNAM') THEN
            DO i=1,LEN(sposnam)
              sposnam(i:i)=blank
            END DO
            sposnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'FPOSNAM') THEN
            DO i=1,LEN(fposnam)
              fposnam(i:i)=blank
            END DO
            fposnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'BPARNAM') THEN
            DO i=1,LEN(bparnam)
              bparnam(i:i)=blank
            END DO
            bparnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'SPARNAM') THEN
            DO i=1,LEN(sparnam)
              sparnam(i:i)=blank
            END DO
            sparnam=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'USRNAME') THEN
            DO i=1,LEN(usrname)
              usrname(i:i)=blank
            END DO
            usrname=TRIM(ADJUSTL(Cval(Nval)))
          END IF
        END IF
      END DO
  10  WRITE (out,50) line
      STOP
  20  CLOSE (inp)
!
!  Set switch to create history file.
!
      DO ng=1,Ngrids
#ifdef SOLVE3D
        DO i=1,NV
          IF (Hout(i,ng)) LdefHIS(ng)=.true.
        END DO
#else
        IF (Hout(idFsur,ng).or.Hout(idUbar,ng).or.Hout(idVbar,ng)) THEN
          LdefHIS(ng)=.true.
        END IF
#endif
      END DO
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        WRITE (out,60) TRIM(title)

        DO ng=1,Ngrids
!
!  Report grid size and domain decomposition.  Check for correct tile
!  decomposition.
!
#ifdef DISTRIBUTE
          WRITE (out,70) ng, Lm(ng), Mm(ng), N(ng), numnodes,           &
     &                   NtileI(ng), NtileJ(ng)
          IF ((NtileI(ng)*NtileJ(ng)).ne.numnodes) THEN
            WRITE (out,80) ng
            exit_flag=6
            RETURN
          END IF
#else
          WRITE (out,90) ng, Lm(ng), Mm(ng), N(ng), numthreads,         &
     &                   NtileI(ng), NtileJ(ng)
          IF (NtileI(ng)*NtileJ(ng).le.0) THEN
            WRITE (out,100) ng
            exit_flag=6
            RETURN
          END IF
          IF (MOD(NtileI(ng)*NtileJ(ng),numthreads).ne.0) THEN
            WRITE (out,100) ng
            exit_flag=6
            RETURN
          END IF
#endif
!
!  Report physical parameters.
!
          WRITE (out,110) ng
          WRITE (out,120) ntimes, 'ntimes',                             &
     &          'Number of timesteps for 3-D equations.'
          WRITE (out,140) dt(ng), 'dt',                                 &
     &          'Timestep size (s) for 3-D equations.'
          WRITE (out,130) ndtfast(ng), 'ndtfast',                       &
     &          'Number of timesteps for 2-D equations between',        &
     &          'each 3D timestep.'
          WRITE (out,120) ERstr, 'ERstr',                               &
     &          'Starting ensemble/perturbation run number.'
          WRITE (out,120) ERend, 'ERend',                               &
     &          'Ending ensemble/perturbation run number.'
          WRITE (out,120) nrrec, 'nrrec',                               &
     &          'Number of restart records to read from disk.'
          WRITE (out,170) lcycle(ng), 'lcycle',                         &
     &          'Switch to recycle time-records in restart file.'
          WRITE (out,130) nrst(ng), 'nrst',                             &
     &          'Number of timesteps between the writing of data',      &
     &          'into restart fields.'
          WRITE (out,130) ninfo(ng), 'ninfo',                           &
     &          'Number of timesteps between print of information',     &
     &          'to standard output.'
#ifdef STATIONS
          WRITE (out,130) nsta(ng), 'nsta',                             &
     &          'Number of timesteps between the writing of data',      &
     &          'the stations file.'
#endif
#ifdef FLOATS
          WRITE (out,130) nflt(ng), 'nflt',                             &
     &          'Number of timesteps between the writing of data',      &
     &          'into floats file.'
#endif
          WRITE (out,170) ldefout(ng), 'ldefout',                       &
     &          'Switch to create a new output NetCDF file(s).'
          WRITE (out,130) nhis(ng), 'nhis',                             &
     &          'Number of timesteps between the writing fields',       &
     &          'into history file.'
          IF (ndefhis(ng).gt.0) THEN
            WRITE (out,130) ndefhis(ng), 'ndefhis',                     &
     &            'Number of timesteps between creation of new',        &
     &            'history files.'
          END IF
#ifdef AVERAGES
          WRITE (out,130) ntsavg(ng), 'ntsavg',                         &
     &          'Starting timestep for the accumulation of output',     &
     &          'time-averaged data.'
          WRITE (out,130) navg(ng), 'navg',                             &
     &          'Number of timesteps between the writing of',           &
     &          'time-averaged data into averages file.'
          IF (ndefavg(ng).gt.0) THEN
            WRITE (out,130) ndefavg(ng), 'ndefavg',                     &
     &            'Number of timesteps between creation of new',        &
     &            'time-averaged files.'
          END IF
#endif
#ifdef SOLVE3D
# ifdef TS_DIF2
          DO itrc=1,NAT
            WRITE (out,190) tnu2(itrc,ng), 'tnu2', itrc,                &
     &            'Horizontal, harmonic mixing coefficient (m2/s)',     &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
# ifdef TS_DIF4
          DO itrc=1,NAT
            WRITE (out,190) tnu4(itrc,ng), 'tnu4', itrc,                &
     &            'Horizontal, biharmonic mixing coefficient (m4/s)',   &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# endif
#endif
#ifdef UV_VIS2
          WRITE (out,210) visc2(ng), 'visc2',                           &
     &          'Horizontal, harmonic mixing coefficient (m2/s)',       &
     &          'for momentum.'
#endif
#ifdef UV_VIS4
          WRITE (out,210) visc2(ng), 'visc4',                           &
     &          'Horizontal, biharmonic mixing coefficient (m4/s)',     &
     &          'for momentum.'
#endif
#ifdef SOLVE3D
          DO itrc=1,NAT
            WRITE (out,190) Akt_bak(itrc,ng), 'Akt_bak', itrc,          &
     &            'Background vertical mixing coefficient (m2/s)',      &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
          WRITE (out,210) Akv_bak(ng), 'Akv_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for momentum.'
# if defined MY25_MIXING || defined GLS_MIXING
          WRITE (out,210) Akk_bak(ng), 'Akk_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for turbulent energy.'
          WRITE (out,210) Akp_bak(ng), 'Akp_bak',                       &
     &          'Background vertical mixing coefficient (m2/s)',        &
     &          'for turbulent generic statistical field.'
#  ifdef TKE_DIF2
          WRITE (out,210) tkenu2(ng), 'tkenu2',                         &
     &          'Horizontal, harmonic mixing coefficient (m2/s)',       &
     &          'for turbulent energy.'
#  endif
#  ifdef TKE_DIF4
          WRITE (out,210) tkenu4(ng), 'tkenu4',                         &
     &          'Horizontal, biharmonic mixing coefficient (m4/s)',     &
     &          'for turbulent energy.'
#  endif
# endif
# ifdef GLS_MIXING
          WRITE (out,140) gls_p(ng), 'gls_p',                           &
     &          'GLS stability exponent.'
          WRITE (out,140) gls_m(ng), 'gls_m',                           &
     &          'GLS turbulent kinetic energy exponent.'
          WRITE (out,140) gls_n(ng), 'gls_n',                           &
     &          'GLS turbulent length scale exponent.'
          WRITE (out,200) gls_Kmin(ng), 'gls_Kmin',                     &
     &          'GLS minimum value of turbulent kinetic energy.'
          WRITE (out,200) gls_Pmin(ng), 'gls_Pmin',                     &
     &          'GLS minimum value of dissipation.'
          WRITE (out,200) gls_cmu0(ng), 'gls_cmu0',                     &
     &          'GLS stability coefficient.'
          WRITE (out,200) gls_c1(ng), 'gls_c1',                         &
     &          'GLS shear production coefficient.'
          WRITE (out,200) gls_c2(ng), 'gls_c2',                         &
     &          'GLS dissipation coefficient.'
          WRITE (out,200) gls_c3m(ng), 'gls_c3m',                       &
     &          'GLS stable buoyancy production coefficient.'
          WRITE (out,200) gls_c3p(ng), 'gls_c3p',                       &
     &          'GLS unstable buoyancy production coefficient.'
          WRITE (out,200) gls_sigk(ng), 'gls_sigk',                     &
     &          'GLS constant Schmidt number for TKE.'
          WRITE (out,200) gls_sigp(ng), 'gls_sigp',                     &
     &          'GLS constant Schmidt number for PSI.'
# endif
#endif
          IF (rdrg(ng).ne.0.0_r8) THEN
            WRITE (out,200) rdrg(ng), 'rdrg',                           &
     &            'Linear bottom drag coefficient (m/s).'
          END IF
          IF (rdrg2(ng).ne.0.0_r8) THEN
            WRITE (out,200) rdrg2(ng), 'rdrg2',                         &
     &            'Quadratic bottom drag coefficient.'
          END IF
          IF (Zob(ng).ne.0.0_r8) THEN
            WRITE (out,200) Zob(ng), 'Zob',                             &
     &            'Bottom roughness (m).'
          END IF
#ifdef BBL_MODEL
          IF (Zob(ng).le.0.0_r8) THEN
            WRITE (out,260) 'Zob.',                                     &
     &            'It must be greater than zero when BBL is activated.'
            STOP
          END IF
#endif
#ifdef SOLVE3D
# ifdef GLS_MIXING
          WRITE (out,200) Zos(ng), 'Zos',                               &
     &          'Surface roughness (m).'
# endif
# if defined LMD_SKPP || defined SOLAR_SOURCE
          WRITE (out,120) lmd_Jwt(ng), 'lmd_Jwt',                       &
     &          'Jerlov water type.'
          IF ((lmd_Jwt(ng).lt.1).or.(lmd_Jwt(ng).gt.5)) THEN
            WRITE (out,260) 'lmd_Jwt.',                                 &
     &            'It must between one and five.'
            STOP
          END IF
# endif
# ifdef BODYFORCE
          WRITE (out,130) levsfrc(ng), 'levsfrc',                       &
     &          'Deepest level to apply surface stress as a',           &
     &          'bodyforce.'
          IF ((levsfrc(ng).lt.1).or.(levsfrc(ng).gt.N(ng))) THEN
            WRITE (out,260) 'levsfrc.',                                 &
     &            'Out of range surface bodyforce level.'
            STOP
          END IF
          WRITE (out,130) levbfrc(ng), 'levbfrc',                       &
     &          'Shallowest level to apply bottom stress as a',         &
     &          'bodyforce.'
          IF ((levbfrc(ng).lt.1).or.(levbfrc(ng).gt.N(ng))) THEN
            WRITE (out,260) 'levbfrc.',                                 &
     &            'Out of range bottom bodyforce level.'
            STOP
          END IF
# endif
#endif
#ifdef SOLVE3D
          WRITE (out,200) theta_s(ng), 'theta_s',                       &
     &          'S-coordinate surface control parameter.'
          WRITE (out,200) theta_b(ng), 'theta_b',                       &
     &          'S-coordinate bottom  control parameter.'
          WRITE (out,160) Tcline(ng), 'Tcline',                         &
     &          'S-coordinate surface/bottom layer width (m) used',     &
     &          'in vertical coordinate stretching.'
#endif
          WRITE (out,140) rho0, 'rho0',                                 &
     &          'Mean density (kg/m3) for Boussinesq approximation.'
          WRITE (out,140) dstart, 'dstart',                             &
     &          'Timestamp assigned to model initialization (days).'
          WRITE (out,150) time_ref, 'time_ref',                         &
     &          'Reference time for units attribute (yyyymmdd.dd)'
#ifdef SOLVE3D
          DO itrc=1,NAT
            WRITE (out,190) Tnudg(itrc,ng), 'Tnudg', itrc,              &
     &            'Nudging/relaxation time scale (days)',               &
     &            'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
# ifdef SCORRECTION
          IF (Tnudg(isalt,ng).le.0.0_r8) THEN
            WRITE (out,260) 'Tnudg(isalt).',                            &
     &            'Must be greater than zero for salt flux correction.'
            STOP
          END IF
# endif
#endif
          WRITE (out,210) Znudg(ng), 'Znudg',                           &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for free-surface.'
          WRITE (out,210) M2nudg(ng), 'M2nudg',                         &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for 2D momentum.'
#ifdef SOLVE3D
          WRITE (out,210) M3nudg(ng), 'M3nudg',                         &
     &          'Nudging/relaxation time scale (days)',                 &
     &          'for 3D momentum.'
#endif
          WRITE (out,210) obcfac(ng), 'obcfac',                         &
     &          'Factor between passive and active',                    &
     &          'open boundary conditions.'
#ifdef SOLVE3D
          WRITE (out,140) T0(ng), 'T0',                                 &
     &          'Background potential temperature (C) constant.'
          WRITE (out,140) S0(ng), 'S0',                                 &
     &          'Background salinity (PSU) constant.'
# ifndef NONLIN_EOS
          WRITE (out,160) R0(ng), 'R0',                                 &
     &          'Background density (kg/m3) used in linear Equation',   &
     &          'of State.'
          WRITE (out,200) Tcoef(ng), 'Tcoef',                           &
     &          'Thermal expansion coefficient (1/Celsius).'
          WRITE (out,200) Scoef(ng), 'Scoef',                           &
     &          'Saline contraction coefficient (1/PSU).'
# endif
#endif
          WRITE (out,160) gamma2, 'gamma2',                             &
     &          'Slipperiness variable: free-slip (1.0) or ',           &
     &          '                     no-slip (-1.0).'
          IF (Hout(idFsur,ng)) WRITE (out,170) Hout(idFsur,ng),         &
     &       'Hout(idFsur)',                                            &
     &       'Write out free-surface.'
          IF (Hout(idUbar,ng)) WRITE (out,170) Hout(idUbar,ng),         &
     &       'Hout(idUbar)',                                            &
     &       'Write out 2D U-momentum component.'
          IF (Hout(idVbar,ng)) WRITE (out,170) Hout(idVbar,ng),         &
     &       'Hout(idVbar)',                                            &
     &       'Write out 2D V-momentum component.'
#ifdef SOLVE3D
          IF (Hout(idUvel,ng)) WRITE (out,170) Hout(idUvel,ng),         &
     &       'Hout(idUvel)',                                            &
     &       'Write out 3D U-momentum component.'
          IF (Hout(idVvel,ng)) WRITE (out,170) Hout(idVvel,ng),         &
     &       'Hout(idVvel)',                                            &
     &       'Write out 3D V-momentum component.'
          IF (Hout(idWvel,ng)) WRITE (out,170) Hout(idWvel,ng),         &
     &       'Hout(idWvel)',                                            &
     &       'Write out W-momentum component.'
          IF (Hout(idOvel,ng)) WRITE (out,170) Hout(idOvel,ng),         &
     &       'Hout(idOvel)',                                            &
     &       'Write out omega vertical velocity.'
          DO itrc=1,NAT
            IF (Hout(idTvar(itrc),ng)) WRITE (out,180)                  &
     &          Hout(idTvar(itrc),ng), 'Hout(idTvar)',                  &
     &          'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
          END DO
#endif
          IF (Hout(idUsms,ng)) WRITE (out,170) Hout(idUsms,ng),         &
     &       'Hout(idUsms)',                                            &
     &       'Write out surface U-momentum stress.'
          IF (Hout(idVsms,ng)) WRITE (out,170) Hout(idVsms,ng),         &
     &       'Hout(idVsms)',                                            &
     &       'Write out surface V-momentum stress.'
          IF (Hout(idUbms,ng)) WRITE (out,170) Hout(idUbms,ng),         &
     &       'Hout(idUbms)',                                            &
     &       'Write out bottom U-momentum stress.'
          IF (Hout(idUbms,ng)) WRITE (out,170) Hout(idVbms,ng),         &
     &       'Hout(idVbms)',                                            &
     &       'Write out bottom V-momentum stress.'
#ifdef BBL_MODEL
          IF (Hout(idUbws,ng)) WRITE (out,170) Hout(idUbws,ng),         &
     &       'Hout(idUbws)',                                            &
     &       'Write out wind-induced, bottom U-wave stress.'
          IF (Hout(idVbws,ng)) WRITE (out,170) Hout(idVbws,ng),         &
     &       'Hout(idVbws)',                                            &
     &       'Write out wind-induced, bottom V-wave stress.'
          IF (Hout(idAbed,ng)) WRITE (out,170) Hout(idAbed,ng),         &
     &       'Hout(idAbed)',                                            &
     &       'Write out bed wave excursion amplitude.'
          IF (Hout(idUbed,ng)) WRITE (out,170) Hout(idUbed,ng),         &
     &       'Hout(idUbed)',                                            &
     &       'Write out bed wave orbital U-velocity.'
          IF (Hout(idVbed,ng)) WRITE (out,170) Hout(idVbed,ng),         &
     &       'Hout(idVbed)',                                            &
     &       'Write out bed wave orbital V-velocity.'
          IF (Hout(idUbot,ng)) WRITE (out,170) Hout(idUbot,ng),         &
     &       'Hout(idUbot)',                                            &
     &       'Write out bottom U-momentum above bed.'
          IF (Hout(idVbot,ng)) WRITE (out,170) Hout(idVbot,ng),         &
     &       'Hout(idVbot)',                                            &
     &       'Write out bottom V-momentum above bed.'
          IF (Hout(idHrip,ng)) WRITE (out,170) Hout(idHrip,ng),         &
     &       'Hout(idHrip)',                                            &
     &       'Write out bed ripple height.'
          IF (Hout(idLrip,ng)) WRITE (out,170) Hout(idLrip,ng),         &
     &       'Hout(idLrip)',                                            &
     &       'Write out bed ripple length.'
          IF (Hout(idZnot,ng)) WRITE (out,170) Hout(idZnot,ng),         &
     &       'Hout(idZnot)',                                            &
     &       'Write out bottom roughness.'
          IF (Hout(idZapp,ng)) WRITE (out,170) Hout(idZapp,ng),         &
     &       'Hout(idZapp)',                                            &
     &       'Write out apparent bottom roughness.'
#endif
#ifdef SOLVE3D
          IF (Hout(idTsur(itemp),ng)) WRITE (out,170)                   &
     &        Hout(idTsur(itemp),ng), 'Hout(idTsur)',                   &
     &       'Write out surface net heat flux.'
# ifdef SHORTWAVE
          IF (Hout(idSrad,ng)) WRITE (out,170) Hout(idSrad,ng),         &
     &       'Hout(idSrad)',                                            &
     &       'Write out shortwave radiation flux.'
# endif
# ifdef BULK_FLUXES
          IF (Hout(idLrad,ng)) WRITE (out,170) Hout(idLrad,ng),         &
     &       'Hout(idLrad)',                                            &
     &       'Write out longwave radiation flux.'
          IF (Hout(idLhea,ng)) WRITE (out,170) Hout(idLhea,ng),         &
     &       'Hout(idLhea)',                                            &
     &       'Write out latent heat flux.'
          IF (Hout(idShea,ng)) WRITE (out,170) Hout(idShea,ng),         &
     &       'Hout(idShea)',                                            &
     &       'Write out sensible heat flux.'
# endif
          IF (Hout(idDano,ng)) WRITE (out,170) Hout(idDano,ng),         &
     &       'Hout(idDano)',                                            &
     &       'Write out density anomaly.'
          IF (Hout(idVvis,ng)) WRITE (out,170) Hout(idVvis,ng),         &
     &       'Hout(idVvis)',                                            &
     &       'Write out vertical viscosity coefficient.'
          IF (Hout(idTdif,ng)) WRITE (out,170) Hout(idTdif,ng),         &
     &       'Hout(idTdif)',                                            &
     &       'Write out vertical T-diffusion coefficient.'
          IF (Hout(idSdif,ng)) WRITE (out,170) Hout(idSdif,ng),         &
     &       'Hout(idSdif)',                                            &
     &       'Write out vertical S-diffusion coefficient.'
# ifdef LMD_SKPP
          IF (Hout(idHsbl,ng)) WRITE (out,170) Hout(idHsbl,ng),         &
     &       'Hout(idHsbl)',                                            &
     &       'Write out depth of surface boundary layer.'
# endif
# ifdef LMD_BKPP
          IF (Hout(idHbbl,ng)) WRITE (out,170) Hout(idHbbl,ng),         &
     &       'Hout(idHbbl)',                                            &
     &       'Write out depth of bottom boundary layer.'
# endif
# if defined GLS_MIXING || defined MY25_MIXING
          IF (Hout(idMtke,ng)) WRITE (out,170) Hout(idMtke,ng),         &
     &       'Hout(idMtke)',                                            &
     &       'Write out turbulent kinetic energy.'
          IF (Hout(idMtls,ng)) WRITE (out,170) Hout(idMtls,ng),         &
     &       'Hout(idMtls)',                                            &
     &       'Write out turbulent generic length-scale.'
# endif
#endif
!
!-----------------------------------------------------------------------
!  Report output/input files and check availability of input files.
!-----------------------------------------------------------------------
!
          WRITE (out,220)
          WRITE (out,230) '          Output Restart File:  ',           &
     &                    TRIM(rstname(ng))
          IF (LdefHIS(ng)) THEN
            IF (ndefhis(ng).eq.0) THEN
              WRITE (out,230) '          Output History File:  ',       &
     &                        TRIM(hisname(ng))
            ELSE
              Lstr=LEN_TRIM(hisname(ng))
              WRITE (out,230) '     Prefix for History Files:  ',       &
     &                        hisname(ng)(1:Lstr-3)
            END IF
          END IF
#ifdef AVERAGES
          IF (ndefavg(ng).eq.0) THEN
            WRITE (out,230) '         Output Averages File:  ',         &
     &                      TRIM(avgname(ng))
          ELSE
            Lstr=LEN_TRIM(avgname(ng))
            WRITE (out,230) '    Prefix for Averages Files:  ',         &
     &                      avgname(ng)(1:Lstr-3)
          END IF
#endif
#ifdef STATIONS
          WRITE (out,230) '         Output Stations File:  ',           &
     &                    TRIM(staname(ng))
#endif
#ifdef FLOATS
          WRITE (out,230) '           Output Floats File:  ',           &
     &                    TRIM(fltname(ng))
#endif
#ifndef ANA_GRID
          fname=grdname(ng)
          OPEN (Tunit, FILE=TRIM(fname), FORM='unformatted',            &
     &          STATUS='old', ERR=30)
          WRITE (out,230) '              Input Grid File:  ',           &
     &                    TRIM(fname)
          CLOSE (Tunit)
#endif
#ifndef ANA_INITIAL
          fname=ininame(ng)
          OPEN (Tunit, FILE=TRIM(fname), FORM='unformatted',            &
     &          STATUS='old', ERR=30)
          WRITE (out,230) '           Input Initial File:  ',           &
     &                    TRIM(fname)
          CLOSE (Tunit)
#endif
#ifdef FRC_FILE
          DO i=1,nFfiles(ng)
            fname=frcname(i,ng)
            OPEN (Tunit, FILE=TRIM(fname), FORM='unformatted',          &
     &            STATUS='old', ERR=30)
            WRITE (out,290) '        Input Forcing File ', i, ':  ',    &
     &                      TRIM(fname)
            CLOSE (Tunit)
          END DO
#endif
#if (defined  ZCLIMATOLOGY && !defined ANA_SSH)     || \
    (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA) || \
    (defined  TCLIMATOLOGY && !defined ANA_TCLIMA)  || \
    (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
          fname=clmname(ng)
          OPEN (Tunit, FILE=TRIM(fname), FORM='unformatted',            &
     &          STATUS='old', ERR=30)
          WRITE (out,230) '       Input Climatology File:  ',           &
     &                    TRIM(fname)
          CLOSE (Tunit)
#endif
#ifdef OBC_DATA
          fname=bryname(ng)
          OPEN (Tunit, FILE=TRIM(fname), FORM='unformatted',            &
     &          STATUS='old', ERR=30)
          WRITE (out,230) '          Input Boundary File:  ',           &
     &                    TRIM(fname)
          CLOSE (Tunit)
#endif
#ifdef STATIONS
          fname=sposnam
          OPEN (Tunit, FILE=TRIM(fname), FORM='formatted',              &
     &          STATUS='old', ERR=30)
          WRITE (out,230) '       Station positions File:  ',           &
     &                    TRIM(fname)
          CLOSE (Tunit)
#endif
#if defined ASSIMILATION || defined NUDGING
          fname=aparnam
          OPEN (Tunit, FILE=TRIM(fname), FORM='formatted',              &
     &          STATUS='old', ERR=30)
          WRITE (out,230) ' Assimilation Parameters File:  ',           &
     &                    TRIM(fname)
          CLOSE (Tunit)
#endif
#ifdef FLOATS
          fname=fposnam
          OPEN (Tunit, FILE=TRIM(fname), FORM='formatted',              &
     &          STATUS='old', ERR=30)
          WRITE (out,230) 'Initial Floats Positions File:  ',           &
     &                    TRIM(fname)
          CLOSE (Tunit)
#endif
#ifdef BIOLOGY
          fname=bparnam
          OPEN (Tunit, FILE=TRIM(fname), FORM='formatted',              &
     &          STATUS='old', ERR=30)
          WRITE (out,230) '      Biology Parameters File:  ',           &
     &                    TRIM(fname)
          CLOSE (Tunit)
#endif
          fname=varname
          WRITE (out,230) ' IO Variable Information File:  ',           &
     &                    TRIM(fname)
          GOTO 40
  30      WRITE (out,270) TRIM(fname)
          STOP
  40      CONTINUE
        END DO
!
!-----------------------------------------------------------------------
!  Report generic USER parameters.
!-----------------------------------------------------------------------
!
        IF (Nuser.gt.0) THEN
          WRITE (out,240)
          DO i=1,Nuser
            WRITE (out,250) user(i), i, i
          END DO
        END IF
        WRITE (out,230) '       Input/Output USER File:  ',             &
     &                  TRIM(usrname)

      END IF

  50  FORMAT (/,' READ_PhyPar - Error while processing line: ',/,a)
  60  FORMAT (/,1x,a)
  70  FORMAT (/,' Resolution, Grid ',i2.2,': ',i4.4,'x',i4.4,'x',i3.3,  &
     &        ',',2x,'Parallel Nodes: ',i2,',',2x,'Tiling: ',i2.2,      &
     &        'x',i2.2)
  80  FORMAT (/,' ROMS/TOMS: Wrong choice of domain ',i2.2,1x,          &
     &        'partition or number of parallel threads.',/,7x,          &
     &        'NtileI * NtileJ  must be equal to the number of ',       &
     &        'parallel nodes.',/,7x,'Change -np value to mpirun or',/, &
     &        7x,'change domain partition in input script.')
  90  FORMAT (/,' Resolution, Grid ',i2.2,': ',i4.4,'x',i4.4,'x',i3.3,  &
     &        ',',2x,'Parallel Threads: ',i2,',',2x,'Tiling: ',i2.2,    &
     &        'x',i2.2)
 100  FORMAT (/,' ROMS/TOMS: Wrong choice of domain ',i2.2,1x,          &
     &        'partition or number of parallel threads.',/,7x,          &
     &        'NtileI*NtileJ must be a positive multiple of the ',      &
     &        'number of threads.',/,7x,'Change number of threads ',    &
     &        '(environment variable) or', /,7x,'change domain ',       &
     &        'partition in input script.')
 110  FORMAT (/,/,' Physical Parameters, Grid: ',i2.2,                  &
     &        /,' =============================',/)  
 120  FORMAT (1x,i10,2x,a,t28,a)
 130  FORMAT (1x,i10,2x,a,t28,a,/,t30,a)
 140  FORMAT (f11.3,2x,a,t28,a)
 150  FORMAT (f11.2,2x,a,t28,a)
 160  FORMAT (f11.3,2x,a,t28,a,/,t30,a)
 170  FORMAT (10x,l1,2x,a,t28,a)
 180  FORMAT (10x,l1,2x,a,t28,a,i2.2,':',1x,a)
 190  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t28,a,/,t30,a,i2.2,':',1x,a)
 200  FORMAT (1p,e11.4,2x,a,t28,a)
 210  FORMAT (1p,e11.4,2x,a,t28,a,/,t30,a)
 220  FORMAT (/,' Output/Input Files:',/)
 230  FORMAT (2x,a,a)
 240  FORMAT (/,' Generic User Parameters:',/)
 250  FORMAT (1p,e11.4,2x,'user(',i2.2,')',t28,                         &
     &        'User parameter ',i2.2,'.')
 260  FORMAT (/,' READ_PHYPAR - Invalid input parameter, ',a,/,15x,a)
 270  FORMAT (/,' READ_PHYPAR - could not find input file:  ',a)
 280  FORMAT (/,' READ_PHYPAR - variable info not yet loaded, ', a)
 290  FORMAT (2x,a,i2.2,a,a)

      RETURN
      END SUBROUTINE read_PhyPar

#ifdef BIOLOGY
      SUBROUTINE read_BioPar (inp, out, Lwrite)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine reads in biological model input parameters.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval, i, ng, itrc, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(NBT,Ngrids) :: Lbio

      real(r8), dimension(NBT,Ngrids) :: Rbio
 
      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=80) :: line
      character (len=80), dimension(100) :: Cval

# ifdef BIO_FASHAM
!
!-----------------------------------------------------------------------
!  Read in Fasham biological model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lbiology') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lbiology)
          ELSE IF (TRIM(KeyWord).eq.'BioIter') THEN 
            Npts=load_i(Nval, Rval, Ngrids, BioIter)
          ELSE IF (TRIM(KeyWord).eq.'AttSW') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttSW)
          ELSE IF (TRIM(KeyWord).eq.'AttChl') THEN
            Npts=load_r(Nval, Rval, Ngrids, AttChl)
          ELSE IF (TRIM(KeyWord).eq.'PARfrac') THEN
            Npts=load_r(Nval, Rval, Ngrids, PARfrac)
          ELSE IF (TRIM(KeyWord).eq.'Vp0') THEN
            Npts=load_r(Nval, Rval, Ngrids, Vp0)
          ELSE IF (TRIM(KeyWord).eq.'I_thNH4') THEN
            Npts=load_r(Nval, Rval, Ngrids, I_thNH4)
          ELSE IF (TRIM(KeyWord).eq.'D_p5NH4') THEN
            Npts=load_r(Nval, Rval, Ngrids, D_p5NH4)
          ELSE IF (TRIM(KeyWord).eq.'NitriR') THEN
            Npts=load_r(Nval, Rval, Ngrids, NitriR)
          ELSE IF (TRIM(KeyWord).eq.'K_NO3') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_NO3)
          ELSE IF (TRIM(KeyWord).eq.'K_NH4') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_NH4)
          ELSE IF (TRIM(KeyWord).eq.'K_Phy') THEN
            Npts=load_r(Nval, Rval, Ngrids, K_Phy)
          ELSE IF (TRIM(KeyWord).eq.'Chl2C_m') THEN
            Npts=load_r(Nval, Rval, Ngrids, Chl2C_m)
          ELSE IF (TRIM(KeyWord).eq.'ChlMin') THEN
            Npts=load_r(Nval, Rval, Ngrids, ChlMin)
          ELSE IF (TRIM(KeyWord).eq.'PhyCN') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyCN)
          ELSE IF (TRIM(KeyWord).eq.'PhyIP') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyIP)
          ELSE IF (TRIM(KeyWord).eq.'PhyIS') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyIS)
          ELSE IF (TRIM(KeyWord).eq.'PhyMin') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyMin)
          ELSE IF (TRIM(KeyWord).eq.'PhyMR') THEN
            Npts=load_r(Nval, Rval, Ngrids, PhyMR)
          ELSE IF (TRIM(KeyWord).eq.'ZooAE_C') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooAE_C)
          ELSE IF (TRIM(KeyWord).eq.'ZooAE_N') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooAE_N)
          ELSE IF (TRIM(KeyWord).eq.'ZooBM') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooBM)
          ELSE IF (TRIM(KeyWord).eq.'ZooCN') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooCN)
          ELSE IF (TRIM(KeyWord).eq.'ZooER') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooER)
          ELSE IF (TRIM(KeyWord).eq.'ZooGGE_C') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooGGE_C)
          ELSE IF (TRIM(KeyWord).eq.'ZooGR') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooGR)
          ELSE IF (TRIM(KeyWord).eq.'ZooMin') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooMin)
          ELSE IF (TRIM(KeyWord).eq.'ZooMR') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZooMR)
          ELSE IF (TRIM(KeyWord).eq.'LDeRR') THEN
            Npts=load_r(Nval, Rval, Ngrids, LDeRR)
          ELSE IF (TRIM(KeyWord).eq.'CoagR') THEN
            Npts=load_r(Nval, Rval, Ngrids, CoagR)
          ELSE IF (TRIM(KeyWord).eq.'SDeAR') THEN
            Npts=load_r(Nval, Rval, Ngrids, SDeAR)
          ELSE IF (TRIM(KeyWord).eq.'SDeAR') THEN
            Npts=load_r(Nval, Rval, Ngrids, SDeAR)
          ELSE IF (TRIM(KeyWord).eq.'SDeBR') THEN
            Npts=load_r(Nval, Rval, Ngrids, SDeBR)
          ELSE IF (TRIM(KeyWord).eq.'SDeRR') THEN
            Npts=load_r(Nval, Rval, Ngrids, SDeRR)
          ELSE IF (TRIM(KeyWord).eq.'wPhy') THEN
            Npts=load_r(Nval, Rval, Ngrids, wPhy)
          ELSE IF (TRIM(KeyWord).eq.'wLDeC') THEN
            Npts=load_r(Nval, Rval, Ngrids, wLDeC)
          ELSE IF (TRIM(KeyWord).eq.'wLDeN') THEN
            Npts=load_r(Nval, Rval, Ngrids, wLDeN)
          ELSE IF (TRIM(KeyWord).eq.'wSDeC') THEN
            Npts=load_r(Nval, Rval, Ngrids, wSDeC)
          ELSE IF (TRIM(KeyWord).eq.'wSDeN') THEN
            Npts=load_r(Nval, Rval, Ngrids, wSDeN)
          ELSE IF (TRIM(KeyWord).eq.'TNU2') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu2(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNU4') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                tnu4(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Akt_bak(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'TNUDG') THEN
            Npts=load_r(Nval, Rval, NBT*Ngrids, Rbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idbio(itrc)
                Tnudg(i,ng)=Rbio(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idTvar)') THEN
            Npts=load_l(Nval, Cval, NBT*Ngrids, Lbio)
            DO ng=1,Ngrids
              DO itrc=1,NBT
                i=idTvar(idbio(itrc))
                Hout(i,ng)=Lbio(itrc,ng)
              END DO
            END DO
          END IF
        END IF
      END DO
  10  WRITE (out,30) line
      STOP
  20  CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lbiology(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50) BioIter(ng), 'BioIter',                      &
     &            'Number of iterations for nonlinear convergence.'
            WRITE (out,100) AttSW(ng), 'AttSW',                         &
     &            'Light attenuation of seawater (m-1).'
            WRITE (out,100) AttChl(ng), 'AttChl',                       &
     &            'Light attenuation by chlorophyll (1/(mg_Chl m-2)).'
            WRITE (out,110) PARfrac(ng), 'PARfrac',                     &
     &            'Fraction of shortwave radiation that is',            &
     &            'photosynthetically active (nondimensional).'
            WRITE (out,110) Vp0(ng), 'Vp0',                             &
     &            'Eppley temperature-limited growth parameter',        &
     &            '(nondimensional).'
            WRITE (out,100) I_thNH4(ng), 'I_thNH4',                     &
     &            'Radiation threshold for nitrification (W/m2).'
            WRITE (out,100) D_p5NH4(ng), 'D_p5NH4',                     &
     &            'Half-saturation radiation for nitrification (W/m2).'
            WRITE (out,100) NitriR(ng), 'NitriR',                       &
     &            'Nitrification rate (day-1).'
            WRITE (out,110) K_NO3(ng), 'K_NO3',                         &
     &            'Inverse half-saturation for phytoplankton NO3',      &
     &            'uptake (1/(mmol_N m-3)).'
            WRITE (out,110) K_NH4(ng), 'K_NH4',                         &
     &            'Inverse half-saturation for phytoplankton NH4',      &
     &            'uptake (1/(mmol_N m-3)).'
            WRITE (out,110) K_Phy(ng), 'K_Phy',                         &
     &            'Zooplankton half-saturation constant for ingestion', &
     &            '(day-1).'
            WRITE (out,100) Chl2C_m(ng), 'Chl2C_m',                     &
     &            'Maximum chlorophyll to carbon ratio (mg_Chl/mg_C).'
            WRITE (out,100) ChlMin(ng), 'ChlMin',                       &
     &            'Chlorophyll minimum threshold (mg_Chl/m3).'
            WRITE (out,100) PhyCN(ng), 'PhyCN',                         &
     &            'Phytoplankton Carbon:Nitrogen ratio (mol_C/mol_N).'
            WRITE (out,100) PhyIP(ng), 'PhyIP',                         &
     &            'Phytoplankton NH4 inhibition parameter (1/mmol_N).'
            WRITE (out,110) PhyIS(ng), 'PhyIS',                         &
     &            'Phytoplankton growth, initial slope of P-I curve',   &
     &            '(mg_C/(mg_Chl Watts m-2 day)).'
            WRITE (out,100) PhyMin(ng), 'PhyMin',                       &
     &            'Phytoplankton minimum threshold (mmol_N/m3).'
            WRITE (out,100) PhyMR(ng), 'PhyMR',                         &
     &            'Phytoplankton mortality rate (day-1).'
#  ifdef CARBON
            WRITE (out,110) ZooAE_C(ng), 'ZooAE_C',                     &
     &            'Zooplankton carbon assimilation efficiency',         &
     &            '(nondimensional).'
#  endif
            WRITE (out,110) ZooAE_N(ng), 'ZooAE_N',                     &
     &            'Zooplankton nitrogen assimilation efficiency',       &
     &            '(nondimensional).'
            WRITE (out,100) ZooBM(ng), 'ZooBM',                         &
     &            'Rate for zooplankton basal metabolism (1/day).'
            WRITE (out,100) ZooCN(ng), 'ZooCN',                         &
     &            'Zooplankton Carbon:Nitrogen ratio (mol_C/mol_N).'
            WRITE (out,100) ZooER(ng), 'ZooER',                         &
     &            'Zooplankton specific excretion rate (day-1).'
            WRITE (out,110) ZooGGE_C(ng), 'ZooGGE_C',                   &
     &            'Zooplankton carbon gross growth efficiency',         &
     &            '(nondimensional).'
            WRITE (out,100) ZooGR(ng), 'ZooGR',                         &
     &            'Zooplankton maximum growth rate (day-1).'
            WRITE (out,100) ZooMin(ng), 'ZooMin',                       &
     &            'Zooplankton minimum threshold (mmol_N/m3).'
            WRITE (out,100) ZooER(ng), 'ZooMR',                         &
     &            'Zooplankton mortality rate (day-1).'
            WRITE (out,100) LDeRR(ng), 'LDeRR',                         &
     &            'Large detritus re-mineralization rate (day-1).'
            WRITE (out,100) CoagR(ng), 'CoagR',                         &
     &            'Coagulation rate (day-1).'
            WRITE (out,100) SDeAR(ng), 'SDeAR',                         &
     &            'Small detritus aggregation rate (day-1).'
            WRITE (out,100) SDeBR(ng), 'SDeBR',                         &
     &            'Small detritus breakdown rate (day-1).'
            WRITE (out,100) SDeRR(ng), 'SDeRR',                         &
     &            'Remineralization rate for small detritus (day-1).'
            WRITE (out,100) wPhy(ng), 'wPhy',                           &
     &            'Phytoplankton sinking velocity (m/day).'
#  ifdef CARBON
            WRITE (out,100) wLDeC(ng), 'wLDeC',                         &
     &            'Large carbon detritus sinking velocity (m/day).'
#  endif
            WRITE (out,100) wLDeN(ng), 'wLDeN',                         &
     &            'Large nitrogen detritus sinking velocity (m/day).'
#  ifdef CARBON
            WRITE (out,100) wSDeC(ng), 'wSDeC',                         &
     &            'Small carbon detritus sinking velocity (m/day).'
#  endif
            WRITE (out,100) wSDeN(ng), 'wSDeN',                         &
     &            'Small nitrogen detritus sinking velocity (m/day).'
#  ifdef TS_DIF2
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu2(i,ng), 'tnu2', i,                     &
     &              'Horizontal, harmonic mixing coefficient (m2/s)',   &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
#  ifdef TS_DIF4
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) tnu4(i,ng), 'tnu4', i,                     &
     &              'Horizontal, biharmonic mixing coefficient (m4/s)', &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
#  endif
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE(out,90) Akt_bak(i,ng), 'Akt_bak', i,                &
     &             'Background vertical mixing coefficient (m2/s)',     &
     &             'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              WRITE (out,90) Tnudg(i,ng), 'Tnudg', i,                   &
     &              'Nudging/relaxation time scale (days)',             &
     &              'for tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
            DO itrc=1,NBT
              i=idbio(itrc)
              IF (Hout(idTvar(i),ng)) WRITE (out,60)                    &
     &            Hout(idTvar(i),ng), 'Hout(idTvar)',                   &
     &            'Write out tracer ', i, TRIM(Vname(1,idTvar(i)))
            END DO
          END IF
        END DO
      END IF

  30  FORMAT (/,' READ_BioPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Biology Parameters, Grid: ',i2.2,                   &
     &        /,' ============================',/)  
  50  FORMAT (1x,i10,2x,a,t28,a)
  60  FORMAT (10x,l1,2x,a,t28,a,i2.2,':',1x,a)
  70  FORMAT (f11.3,2x,a,t28,a)
  80  FORMAT (f11.3,2x,a,t28,a,/,t30,a)
  90  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t28,a,/,t30,a,i2.2,':',1x,a)
 100  FORMAT (1p,e11.4,2x,a,t28,a)
 110  FORMAT (1p,e11.4,2x,a,t28,a,/,t30,a)

# endif

      RETURN
      END SUBROUTINE read_BioPar
#endif

#ifdef SEDIMENT
      SUBROUTINE read_SedPar (inp, out, Lwrite)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine reads in cohesive and non-cohesive sediment model      !
!  parameters.                                                         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
      USE mod_sediment
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval, i, ng, itrc, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(Ngrids) :: Lbed
      logical, dimension(NCS,Ngrids) :: Lmud
      logical, dimension(NNS,Ngrids) :: Lsand

      real(r8), dimension(NCS,Ngrids) :: Rmud
      real(r8), dimension(NNS,Ngrids) :: Rsand
 
      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=80) :: line
      character (len=80), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in cohesive and non-cohesive model parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lsediment') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lsediment)
          ELSE IF (TRIM(KeyWord).eq.'MUD_SD50') THEN 
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Sd50(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_CSED') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud )
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Csed(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_SRHO') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Srho(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_WSED') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Wsed(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_ERATE') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                Erate(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAU_CE') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                tau_ce(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TAU_CD') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                tau_cd(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_POROS') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                poros(itrc,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TNU2') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                tnu2(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TNU4') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                tnu4(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                Akt_bak(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'MUD_TNUDG') THEN
            Npts=load_r(Nval, Rval, NCS*Ngrids, Rmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idsed(itrc)
                Tnudg(i,ng)=Rmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idmud)') THEN
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idTvar(idsed(itrc))
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iMfrac)') THEN
            Npts=load_l(Nval, Cval, NCS*Ngrids, Lmud)
            DO ng=1,Ngrids
              DO itrc=1,NCS
                i=idfrac(itrc)
                Hout(i,ng)=Lmud(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_SD50') THEN 
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Sd50(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_CSED') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand )
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Csed(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_SRHO') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Srho(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_WSED') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Wsed(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_ERATE') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                Erate(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TAU_CE') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                tau_ce(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TAU_CD') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                tau_cd(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_POROS') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=NCS+itrc
                poros(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TNU2') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                tnu2(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TNU4') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                tnu4(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_AKT_BAK') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                Akt_bak(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SAND_TNUDG') THEN
            Npts=load_r(Nval, Rval, NNS*Ngrids, Rsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idsed(NCS+itrc)
                Tnudg(i,ng)=Rsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idsand)') THEN
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idTvar(idsed(NCS+itrc))
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iSfrac)') THEN
            Npts=load_l(Nval, Cval, NNS*Ngrids, Lsand)
            DO ng=1,Ngrids
              DO itrc=1,NNS
                i=idfrac(NCS+itrc)
                Hout(i,ng)=Lsand(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(ithck)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(ithck)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iaged)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(iaged)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(iporo)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(iporo)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(idiff)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idSbed(idiff)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(irlen)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idBott(irlen)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(irhgt)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idBott(irhgt)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(ishgt)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idBott(ishgt)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(izNik)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idBott(izNik)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(izbio)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idBott(izbio)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(izbld)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idBott(izbld)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(izbfm)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idBott(izbfm)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Hout(izapp)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Lbed)
            i=idBott(izapp)
            DO ng=1,Ngrids
              Hout(i,ng)=Lbed(ng)
            END DO
          END IF
        END IF
      END DO
  10  WRITE (out,30) line
      STOP
  20  CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lsediment(ng)) THEN
            WRITE (out,40) ng
            WRITE (out,50)
            DO itrc=1,NST
              WRITE (out,60) itrc, Sd50(itrc,ng), Csed(itrc,ng),        &
     &                          Srho(itrc,ng), Wsed(itrc,ng),           &
     &                          Erate(itrc,ng), poros(itrc,ng)
            END DO
            WRITE (out,70)
            DO itrc=1,NST
              i=idsed(itrc)
              WRITE (out,60) itrc, tau_ce(itrc,ng), tau_cd(itrc,ng),    &
     &                          tnu2(i,ng), tnu4(i,ng), Akt_bak(i,ng),  &
     &                          Tnudg(i,ng)
            END DO
            WRITE (out,'(/)')
            DO itrc=1,NST
              i=idTvar(idsed(itrc))
              IF (Hout(i,ng)) WRITE (out,80) Hout(i,ng),                &
     &            'Hout(idTvar)',                                       &
     &            'Write out sediment', itrc, TRIM(Vname(1,i))
            END DO
            DO itrc=1,NST
              i=idfrac(itrc)
              IF (Hout(i,ng)) WRITE (out,80) Hout(i,ng),                &
     &            'Hout(idfrac)',                                       &
     &            'Write out bed fraction, sediment ', itrc,            &
     &            TRIM(Vname(1,i))
            END DO
            DO itrc=1,MBEDP
              i=idSbed(itrc)
              IF (Hout(i,ng)) WRITE (out,80) Hout(i,ng),                &
     &            'Hout(idSbed)',                                       &
     &            'Write out BED property ', itrc, TRIM(Vname(1,i))
            END DO
            DO itrc=1,MBOTP
              i=idBott(itrc)
              IF (Hout(i,ng)) WRITE (out,80) Hout(i,ng),                &
     &            'Hout(idBott)',                                       &
     &            'Write out BOT property ', itrc, TRIM(Vname(1,i))
            END DO
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Scale relevant input parameters
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        DO i=1,NST
          Sd50(i,ng)=Sd50(i,ng)*0.001_r8
          Wsed(i,ng)=Wsed(i,ng)*0.001_r8
          tau_ce(i,ng)=tau_ce(i,ng)/rho0
          tau_cd(i,ng)=tau_cd(i,ng)/rho0
          tnu4(idsed(i),ng)=SQRT(ABS(tnu4(idsed(i),ng)))
          IF (Tnudg(idsed(i),ng).gt.0.0_r8) THEN
            Tnudg(idsed(i),ng)=1.0_r8/(Tnudg(idsed(i),ng)*86400.0_r8)
          ELSE
            Tnudg(idsed(i),ng)=0.0_r8
          END IF
        END DO
      END DO

  30  FORMAT (/,' READ_SedPar - Error while processing line: ',/,a)
  40  FORMAT (/,/,' Sediment Parameters, Grid: ',i2.2,                  &
     &        /,' =============================',/)
  50  FORMAT (/,1x,'Size',5x,'Sd50',8x,'Csed',8x,'Srho',8x,'Wsed',      &
     &        8x,'Erate',7x,'poros',/,1x,'Class',4x,'(mm)',7x,          &
     &        '(mg/l)',6x,'(kg/m3)',5x,'(mm/s)',5x,'(kg/m2/s)',4x,      &
     &        '(nondim)',/)
  60  FORMAT (2x,i2,2x,6(1x,1p,e11.4))
  70  FORMAT (/,9x,'tau_ce',6x,'tau_cd',7x,'tnu2',8x,'tnu4',7x,         &
     &        'Akt_bak',6x,'Tnudg',/,9x,'(N/m2)',6x,'(N/m2)',6x,        &
     &        '(m2/s)',6x,'(m4/s)',7x,'(m2/s)',6x,'(day)',/)
  80  FORMAT (10x,l1,2x,a,t28,a,i2.2,':',1x,a)

      RETURN
      END SUBROUTINE read_SedPar
#endif

#if defined ASSIMILATION || defined NUDGING
      SUBROUTINE read_AssPar (inp, out, Lwrite)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine reads in input model assimilation parameters.       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer, parameter :: Tunit = 2

      integer :: Npts, Nval
      integer :: i, igrid, itrc, k, ng, status

      integer :: decode_line, load_i, load_l, load_r

      logical, dimension(MT,Ngrids) :: Lassi

      real(r8), dimension(MT,Ngrids) :: Rassi

      real(r8), dimension(2,Ngrids) :: pcoef

      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=80) :: line, fname
      character (len=80), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in assimilation parameters. Then, load input data into module.
!  Take into account nested grid configurations.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=10,END=20) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lassimilate') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lassimilate)
          ELSE IF (TRIM(KeyWord).eq.'Emod0') THEN
            Npts=load_r(Nval, Rval, 1, Emod0)
          ELSE IF (TRIM(KeyWord).eq.'Tgrowth') THEN
            Npts=load_r(Nval, Rval, 1, Tgrowth)
          ELSE IF (TRIM(KeyWord).eq.'cor') THEN
            Npts=load_r(Nval, Rval, 1, cor)
          ELSE IF (TRIM(KeyWord).eq.'assi_SSH') THEN
            Npts=load_l(Nval, Cval, Ngrids, assi_SSH)
          ELSE IF (TRIM(KeyWord).eq.'assi_SST') THEN
            Npts=load_l(Nval, Cval, Ngrids, assi_SST)
          ELSE IF (TRIM(KeyWord).eq.'assi_UVsur') THEN
            Npts=load_l(Nval, Cval, Ngrids, assi_UVsur)
          ELSE IF (TRIM(KeyWord).eq.'assi_UV') THEN
            Npts=load_l(Nval, Cval, Ngrids, assi_UV)
          ELSE IF (TRIM(KeyWord).eq.'assi_T') THEN
            Npts=load_l(Nval, Cval, NAT*Ngrids, Lassi)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                assi_T(itrc,ng)=Lassi(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Znudass') THEN
            Npts=load_r(Nval, Rval, Ngrids, Znudass)
          ELSE IF (TRIM(KeyWord).eq.'M2nudass') THEN
            Npts=load_r(Nval, Rval, Ngrids, M2nudass)
          ELSE IF (TRIM(KeyWord).eq.'M3nudass') THEN
            Npts=load_r(Nval, Rval, Ngrids, M3nudass)
          ELSE IF (TRIM(KeyWord).eq.'Tnudass') THEN
            Npts=load_r(Nval, Rval, NAT*Ngrids, Rassi)
            DO ng=1,Ngrids
              DO itrc=1,NAT
                Tnudass(itrc,ng)=Rassi(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'extend_SST') THEN
            Npts=load_l(Nval, Cval, Ngrids, extend_SST)
          ELSE IF (TRIM(KeyWord).eq.'ZmSST') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZmSST)
            DO ng=1,Ngrids
              ZmSST(ng)=-ABS(ZmSST(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ZoSST') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZoSST)
            DO ng=1,Ngrids
              ZoSST(ng)=-ABS(ZoSST(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'npSST') THEN
            Npts=load_i(Nval, Rval, 1, npSST)
          ELSE IF (TRIM(KeyWord).eq.'Pcoef_SST') THEN
            DO ng=1,Ngrids
              DO i=0,npSST(ng)
                read (inp,*) igrid, k, pcoef_SST(k,igrid)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Perr_SST') THEN
            DO ng=1,Ngrids
              DO i=0,npSST(ng)
                read (inp,*) igrid, k, perr_SST(k,igrid)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'extend_UV') THEN
            Npts=load_l(Nval, Cval, Ngrids, extend_UV)
          ELSE IF (TRIM(KeyWord).eq.'ZmUV') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZmUV)
            DO ng=1,Ngrids
              ZmUV(ng)=-ABS(ZmUV(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'ZoUV') THEN
            Npts=load_r(Nval, Rval, Ngrids, ZoUV)
            DO ng=1,Ngrids
              ZoUV(ng)=-ABS(ZoUV(ng))
            END DO
          ELSE IF (TRIM(KeyWord).eq.'npUV') THEN
            Npts=load_i(Nval, Rval, 1, npUV)
          ELSE IF (TRIM(KeyWord).eq.'Pcoef_U') THEN
            DO ng=1,Ngrids
              DO i=0,npUV(ng)
                READ (inp,*) igrid, k, pcoef_U(k,igrid)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Pcoef_V') THEN
            DO ng=1,Ngrids
              DO i=0,npUV(ng)
                READ (inp,*) igrid, k, pcoef_V(k,igrid)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Perr_UV') THEN
            DO ng=1,Ngrids
              DO i=0,npUV(ng)
                READ (inp,*) igrid, k, perr_V(k,igrid)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'SSHname') THEN
            SSHname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'SSTname') THEN
            SSTname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'VSURname') THEN
            VSURname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'VOBSname') THEN
            VOBSname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          ELSE IF (TRIM(KeyWord).eq.'TOBSname') THEN
            TOBSname(Nval)=TRIM(ADJUSTL(Cval(Nval)))
          END IF
        END IF
      END DO
 10   WRITE (out,50) line
      STOP
 20   CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lassimilate(ng)) THEN
            WRITE (out,60) ng
# ifdef ASSIMILATION
            WRITE (out,100) Emod0(ng), 'Emod0',                         &
     &            'Initial model error (percentage).'
            WRITE (out,100) Tgrowth(ng), 'Tgrowth',                     &
     &            'Empirical model error growth scale (days).'
            WRITE (out,100) cor(ng), 'cor',                             &
     &            'Correlation between model and observations.'
            WRITE (out,80) assi_SSH(ng), 'assi_SSH',                    &
     &            'Activate assimilation of SSH.'
#  ifdef SOLVE3D
            WRITE (out,80) assi_SST(ng), 'assi_SST',                    &
     &            'Activate assimilation of SST.'
            DO itrc=1,NAT
              WRITE (out,90) assi_T(itrc,ng), 'assi_T', itrc,           &
     &              'Activate assimilation of tracer', itrc,            &
     &              TRIM(Vname(1,idTvar(itrc)))
            END DO
            WRITE (out,80) assi_UVsur(ng), 'assi_UVsur',                &
     &            'Activate assimilation of surface currents.'
            WRITE (out,80) assi_UV(ng), 'assi_UV',                      &
     &            'Activate assimilation of currents.'
#  endif
# endif
# ifdef SOLVE3D
#  if defined NUDGING_T || defined NUDGING_SST
            DO itrc=1,NAT
              IF (assi_T(itrc,ng)) THEN
                WRITE (out,120) Tnudass(itrc,ng), 'Tnudass', itrc,      &
     &                'Nudging assimilation time scale (days)',         &
     &                'for tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
              END IF
            END DO
#  endif
# endif
# ifdef NUDGING_SSH
            WRITE (out,135) Znudass(ng), 'Znudass',                     &
     &            'Nudging assimilation time scale (days)',             &
     &            'for free-surface.'
# endif
# ifdef SOLVE3D
#  if defined NUDGING_UV || defined NUDGING_UVsur
            WRITE (out,135) M3nudass(ng), 'M3nudass',                   &
     &            'Nudging assimilation time scale (days)',             &
     &            'for 3D momentum.'
#  endif
#  if defined NUDGING_SST || defined ASSIMILATION_SST
            IF (extend_SST(ng)) THEN
              WRITE (out,80) extend_SST(ng), 'extend_SST',              &
     &              'Extend SST vertically.'
              WRITE (out,70) npSST(ng), 'npSST',                        &
     &              'Order of polynomial for SST extension.'
              WRITE (out,100) ZmSST(ng), 'zmSST',                       &
     &              'SST, maximum extension depth (m).'
              WRITE (out,100) ZoSST(ng), 'zoSST',                       &
     &              'E-folding depth (m) to extend SST error variance.'
              IF (KP.lt.npSST(ng)) THEN
                WRITE (out,160) 'KP: ', KP, npSST(ng)
                STOP
              END IF
            END IF
            IF (npSST(ng).gt.0) THEN
              WRITE (out,*)
              DO k=0,npSST(ng)
                WRITE (out,140) k, pcoef_SST(k,ng), perr_SST(k,ng)
              END DO
            END IF
#  endif
#  if defined NUDGING_UVsur || defined ASSIMILATION_UVsur
            IF (extend_UV(ng)) THEN
              WRITE (out,80) extend_UV(ng), 'extend_UV',                &
     &              'Extend surface currents vertically.'
              WRITE (out,70) npUV(ng), 'npUV',                          &
     &             'Order of polynomial for surface currents extension.'
              WRITE (out,100) ZmUV(ng), 'zmUV',                         &
     &              'Surface currents, maximum extension depth (m).'
              WRITE (out,110) ZoSST(ng), 'zoSST',                       &
     &              'E-folding depth (m) to extend surface currents',   &
     &              'error variance.'
              IF (KP.lt.npUV(ng)) THEN
                WRITE (out,160) 'KP: ', KP, npUV(ng)
                STOP
              END IF
            END IF
            IF (npUV(ng).gt.0) THEN
              WRITE (out,*)
              DO k=0,npUV(ng)
                WRITE (out,150) k, pcoef_U(k,ng), pcoef_V(k,ng),        &
     &                             perr_V(k,ng)
              END DO
            END IF
#  endif
# endif
!
!-----------------------------------------------------------------------
!  Report input files and check availability of input files.
!-----------------------------------------------------------------------
!  
            WRITE (out,170)
            WRITE (out,180) ' Assimilation parameters File:  ',         &
     &                      TRIM(aparnam)
# if defined ASSIMILATION_SSH || defined NUDGING_SSH
            fname=SSHname(ng)
            OPEN (Tunit, FILE=TRIM(fname), STATUS='old', ERR=30)
            WRITE (out,180) '      Sea Surface Height File:  ',         &
     &                      TRIM(fname)
            CLOSE (Tunit)
# endif
# ifdef SOLVE3D
#  if defined ASSIMILATION_SST || defined NUDGING_SST
            fname=SSTname(ng)
            OPEN (Tunit, FILE=TRIM(fname), STATUS='old', ERR=30)
            WRITE (out,180) ' Sea Surface Temperature File:  ',         &
     &                      TRIM(fname)
            CLOSE (Tunit)
#  endif
#  if defined ASSIMILATION_T || defined NUDGING_T
            fname=TOBSname(ng)
            OPEN (Tunit, FILE=TRIM(fname), STATUS='old', ERR=30)
            WRITE (out,180) '                 Tracers File:  ',         &
     &                      TRIM(fname)
            CLOSE (Tunit)
#  endif
#  if defined ASSIMILATION_UVsur || defined NUDGING_UVsur
            fname=VSURname(ng)
            OPEN (Tunit, FILE=TRIM(fname), STATUS='old', ERR=30)
            WRITE (out,180) '        Surface Currents File:  ',         &
     &                      TRIM(fname)
            CLOSE (Tunit)
#  endif
#  if defined ASSIMILATION_UV || defined NUDGING_UV
            fname=VOBSname(ng)
            OPEN (Tunit, FILE=TRIM(fname), STATUS='old', ERR=30)
            WRITE (out,180) '     Horizontal Currents File:  ',         &
     &                      TRIM(fname)
            CLOSE (Tunit)
#  endif
# endif
            GOTO 40
  30        WRITE (out,190) TRIM(fname)
            STOP
  40        CONTINUE
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Scale relevant parameters.
!-----------------------------------------------------------------------
!
      DO ng=1,Ngrids
        Emod0(ng)=1.0_r8-0.01_r8*Emod0(ng)
# ifdef NUDGING
        IF (Znudass(ng).gt.0.0_r8) THEN
          Znudass(ng)=1.0_r8/(Znudass(ng)*86400.0_r8)
        ELSE
          Znudass(ng)=0.0_r8
        END IF
        IF (M2nudass(ng).gt.0.0_r8) THEN
          M2nudass(ng)=1.0_r8/(M2nudass(ng)*86400.0_r8)
        ELSE
          M2nudass(ng)=0.0_r8
        END IF
#  ifdef SOLVE3D
        IF (M3nudass(ng).gt.0.0_r8) THEN
          M3nudass(ng)=1.0_r8/(M3nudass(ng)*86400.0_r8)
        ELSE
          M3nudass(ng)=0.0_r8
        END IF
        DO itrc=1,NAT
          IF (Tnudass(itrc,ng).gt.0.0_r8) THEN
            Tnudass(itrc,ng)=1.0_r8/(Tnudass(itrc,ng)*86400.0_r8)
          ELSE
            Tnudass(itrc,ng)=0.0_r8
          END IF
        END DO
#  endif
# endif
      END DO

  50  FORMAT (/,' READ_AssPar - Error while processing line: ',/,a)
  60  FORMAT (/,/,' Assimilation Parameters, Grid: ',i2.2,              &
     &        /,' =================================',/)
  70  FORMAT (1x,i10,2x,a,t28,a)
  80  FORMAT (10x,l1,2x,a,t28,a)
  90  FORMAT (10x,l1,2x,a,'(',i2.2,')',t28,a,1x,i2.2,':',1x,a)
 100  FORMAT (f11.3,2x,a,t28,a)
 110  FORMAT (f11.3,2x,a,t28,a,/,t30,a)
 120  FORMAT (1p,e11.4,2x,a,'(',i2.2,')',t28,a,/,t30,a,i2.2,':',1x,a)
 130  FORMAT (1p,e11.4,2x,a,t28,a)
 135  FORMAT (1p,e11.4,2x,a,t28,a,/,t30,a)
 140  FORMAT (1x,'SST,Err polynomial order ',                           &
     &        i2.2,1x,1p,e15.8,1x,1p,e15.8)
 150  FORMAT (1x,'U,V,Err polynomial order: ',                          &
     &        i2.2,1x,1p,e15.8,1x,1p,e15.8,1x,1p,e15.8)
 160  FORMAT (/,' READ_ASSPAR - too small dimension parameter, ',a,     &
     &        2i4,/,15x,'change file  mod_scalars.F  and recompile.')
 170  FORMAT (/,' Input Assimilation Files:',/)
 180  FORMAT (2x,a,a)
 190  FORMAT (/,' READ_ASSPAR - could not find input file:  ',a)

      RETURN
      END SUBROUTINE read_AssPar
#endif

#ifdef FLOATS
      SUBROUTINE read_FloatsPar (inp, out, Lwrite)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine reads in input station parameters.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_floats
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval
      integer :: i, j, igrid, mc, nc, ng, status

      integer, dimension(Ngrids) :: ncount, nentry

      integer, allocatable :: Fcoor(:,:), Fcount(:,:), Ftype(:,:)

      integer :: decode_line, load_i, load_l, load_r

      real(r8) :: xfloat, yfloat, zfloat

      real(r8), dimension(100) :: Rval

      real(r8), allocatable :: Ft0(:,:), Fx0(:,:), Fy0(:,:), Fz0(:,:)
      real(r8), allocatable :: Fdt(:,:), Fdx(:,:), Fdy(:,:), Fdz(:,:)

      character (len=35) :: frmt
      character (len=40) :: KeyWord
      character (len=80) :: line
      character (len=80), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in initial float locations.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lfloats') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lfloats)
          ELSE IF (TRIM(KeyWord).eq.'FRREC') THEN
            Npts=load_i(Nval, Rval, Ngrids, frrec)
          ELSE IF (TRIM(KeyWord).eq.'NFLOATS') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nfloats)
          ELSE IF (TRIM(KeyWord).eq.'POS') THEN
            Npts=Nfloats(1)
            IF (Ngrids.gt.1) Npts=MAXVAL(Nfloats)
            allocate ( Fcoor (Npts,Ngrids) )
            allocate ( Fcount(Npts,Ngrids) )
            allocate ( Ftype (Npts,Ngrids) )
            allocate ( Ft0(Npts,Ngrids) )
            allocate ( Fx0(Npts,Ngrids) )
            allocate ( Fy0(Npts,Ngrids) )
            allocate ( Fz0(Npts,Ngrids) )
            allocate ( Fdt(Npts,Ngrids) )
            allocate ( Fdx(Npts,Ngrids) )
            allocate ( Fdy(Npts,Ngrids) )
            allocate ( Fdz(Npts,Ngrids) )
            DO ng=1,Ngrids
              CALL allocate_floats (ng)
            END DO
            ncount(1:Ngrids)=0
            nentry(1:Ngrids)=0
            DO WHILE (.true.)
              READ (inp,*,err=30) igrid,                                &
     &                            Fcoor (nentry(igrid)+1,igrid),        &
     &                            Ftype (nentry(igrid)+1,igrid),        &
     &                            Fcount(nentry(igrid)+1,igrid),        &
     &                            Ft0(nentry(igrid)+1,igrid),           &
     &                            Fx0(nentry(igrid)+1,igrid),           &
     &                            Fy0(nentry(igrid)+1,igrid),           &
     &                            Fz0(nentry(igrid)+1,igrid),           &
     &                            Fdt(nentry(igrid)+1,igrid),           &
     &                            Fdx(nentry(igrid)+1,igrid),           &
     &                            Fdy(nentry(igrid)+1,igrid),           &
     &                            Fdz(nentry(igrid)+1,igrid)
              IF (igrid.gt.Ngrids) THEN
                IF (Lwrite) WRITE (out,40) fposnam
                STOP
              END IF
              ncount(igrid)=ncount(igrid)+Fcount(nentry(igrid)+1,igrid)
              nentry(igrid)=nentry(igrid)+1        
            END DO
          END IF
        END IF
      END DO
  20  WRITE (out,50) line
      STOP
  30  CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (ncount(ng).ne.Nfloats(ng)) THEN
            WRITE (stdout,60) ncount(ng), Nfloats(ng)
            STOP
          END IF
          IF (Lfloats(ng)) THEN
            WRITE (out,70) ng
            DO i=1,nentry(ng)
              IF (.not.spherical.and.(Fcoor(i,ng).eq.0)) THEN
                frmt='(i1,i2,i5,f10.4,2f8.2,f8.2,4f9.3)'
              ELSE
                frmt='(i1,i2,i5,f10.4,3f8.2,4f9.3)'
              END IF
              WRITE (out,frmt) Fcoor(i,ng), Ftype(i,ng), Fcount(i,ng),  &
     &                         Ft0(i,ng), Fx0(i,ng), Fy0(i,ng),         &
     &                         Fz0(i,ng), Fdt(i,ng), Fdx(i,ng),         &
     &                         Fdy(i,ng), Fdz(i,ng)
            END DO
            WRITE (out,80) Nfloats(ng),                                 &
     &            'Nfloats',                                            &
     &            'Number of float trajectories to compute.'
          END IF
        END DO
      END IF
!
!-----------------------------------------------------------------------
!  Process initial float locations.
!-----------------------------------------------------------------------
!
!  Set time of float release (seconds after model initialization) and
!  initial float horizontal positions (grid units).  Fill the initial
!  vertical level or depth position.
!
      DO ng=1,Ngrids
        mc=0
        nc=0
        DO i=1,nentry(ng)
          IF (Fcount(i,ng).eq.1) THEN
            nc=nc+1
            FLT(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
            FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
            IF (Fcoor(i,ng).eq.0) THEN
              FLT(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)),        &
     &                                REAL(Lm(ng),r8)+0.5_r8)
              FLT(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)),        &
     &                                REAL(Mm(ng),r8)+0.5_r8)
            ELSE
              mc=mc+1
              FLT(ng)%Flon(mc)=Fx0(i,ng)
              FLT(ng)%Flat(mc)=Fy0(i,ng)
              FLT(ng)%Findex(mc)=nc
            END IF
          ELSE IF (Fcount(i,ng).gt.1) THEN
            DO j=1,Fcount(i,ng)
              nc=nc+1
              IF (Fdt(i,ng).gt.0.0_r8) THEN
                FLT(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng)+              &
     &                                   REAL(j-1,r8)*Fdt(i,ng))*       &
     &                                  day2sec
                FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                IF (Fcoor(i,ng).eq.0) THEN
                  FLT(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,Fx0(i,ng)),    &
     &                                    REAL(Lm(ng),r8)+0.5_r8)
                  FLT(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,Fy0(i,ng)),    &
     &                                    REAL(Mm(ng),r8)+0.5_r8)
                ELSE
                  mc=mc+1
                  FLT(ng)%Flon(mc)=Fx0(i,ng)
                  FLT(ng)%Flat(mc)=Fy0(i,ng)
                  FLT(ng)%Findex(mc)=nc
                END IF
              ELSE
                FLT(ng)%Tinfo(itstr,nc)=(dstart+Ft0(i,ng))*day2sec
                IF (Fdz(i,ng).eq.0.0_r8) THEN
                  FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)
                ELSE
                  IF (Fz0(i,ng).gt.0.0_r8) THEN
                    zfloat=Fz0(i,ng)+REAL(j-1,r8)*Fdz(i,ng)
                    FLT(ng)%Tinfo(izgrd,nc)=MIN(MAX(0.0_r8,zfloat),     &
     &                                      REAL(N(ng),r8))
                  ELSE
                    FLT(ng)%Tinfo(izgrd,nc)=Fz0(i,ng)+                  &
     &                                      REAL(j-1,r8)*Fdz(i,ng)
                  END IF
                END IF
                IF (Fcoor(i,ng).eq.0) THEN
                  xfloat=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                  FLT(ng)%Tinfo(ixgrd,nc)=MIN(MAX(0.5_r8,xfloat),       &
     &                                    REAL(Lm(ng),r8)+0.5_r8)
                  yfloat=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                  FLT(ng)%Tinfo(iygrd,nc)=MIN(MAX(0.5_r8,yfloat),       &
     &                                    REAL(Mm(ng),r8)+0.5_r8)
                ELSE
                  mc=mc+1
                  FLT(ng)%Flon(mc)=Fx0(i,ng)+REAL(j-1,r8)*Fdx(i,ng)
                  FLT(ng)%Flat(mc)=Fy0(i,ng)+REAL(j-1,r8)*Fdy(i,ng)
                END IF
              END IF
            END DO
          END IF
        END DO
        FLT(ng)%Findex(0)=mc
      END DO
!
!  Deallocate local arrays.
!
      deallocate ( Fcoor )
      deallocate ( Fcount )
      deallocate ( Ftype )
      deallocate ( Ft0 )
      deallocate ( Fx0 )
      deallocate ( Fy0 )
      deallocate ( Fz0 )
      deallocate ( Fdt )
      deallocate ( Fdx )
      deallocate ( Fdy )
      deallocate ( Fdz )
!
  40  FORMAT (/,' READ_FloatsPar - Error while reading floats',         &
     &          ' locations in input script: ',a)
  50  FORMAT (/,' READ_FloatsPar - Error while processing line: ',/,a)
  60  FORMAT (/,' READ_FloatsPar - Inconsistent number of floats to',   &
     &          ' process: ', 2i6,/,18x,'change input script.')
  70  FORMAT (/,/,' Floats Initial Locations, Grid: ',i2.2,             &
     &        /,' ==================================',/,/,              &
     &        15x,'Ft0',5x,'Fx0',5x,'Fy0',5x,'Fz0',                     &
     &        6x,'Fdt',6x,'Fdx',6x,'Fdy',6x,'Fdz',/)
 80   FORMAT (/,1x,i10,2x,a,t28,a)

      RETURN
      END SUBROUTINE read_FloatsPar
#endif

#ifdef STATIONS
      SUBROUTINE read_StaPar (inp, out, Lwrite)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine reads in input station parameters.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: inp, out
!
!  Local variable declarations.
!
      integer :: Mstation, Npts, Nval
      integer :: flag, i, igrid, ista, itrc, ng, status

      integer :: decode_line, load_i, load_l, load_r

      real(r8) :: Xpos, Ypos

      logical, dimension(MT,Ngrids) :: Lsta

      integer, dimension(Ngrids) :: is

      real(r8), dimension(100) :: Rval

      character (len=40) :: KeyWord
      character (len=80) :: line
      character (len=80), dimension(100) :: Cval
!
!-----------------------------------------------------------------------
!  Read in stations parameters.
!-----------------------------------------------------------------------
!
      DO WHILE (.true.)
        READ (inp,'(a)',ERR=20,END=30) line
        status=decode_line(line, KeyWord, Nval, Cval, Rval)
        IF (status.gt.0) THEN
          IF (TRIM(KeyWord).eq.'Lstations') THEN 
            Npts=load_l(Nval, Cval, Ngrids, Lstations)
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idWvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idWvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idOvel)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idOvel,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbar)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbar,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbar)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbar,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idFsur)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idFsur,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTvar)') THEN
            Npts=load_l(Nval, Cval, MT*Ngrids, Lsta)
            DO ng=1,Ngrids
              DO itrc=1,NT(ng)
                Sout(idTvar(itrc),ng)=Lsta(itrc,ng)
              END DO
            END DO
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUsms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUsms,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVsms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVsms,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbms,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbms)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbms,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbws)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbws,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbws)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbws,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idAbed)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idAbed,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbed)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbed,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbed)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbed,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idUbot)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idUbot,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVbot)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVbot,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idHrip)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idHrip,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idLrip)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idLrip,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idZnot)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idZnot,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idZapp)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idZapp,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTsur)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idTsur(itemp),1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idLhea)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idLhea,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idShea)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idShea,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idLrad)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idLrad,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idSrad)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idSrad,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idDano)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idDano,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idVvis)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idVvis,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idTdif)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idTdif,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idSdif)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idSdif,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idHsbl)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idHsbl,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idHbbl)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idHbbl,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idMtke)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idMtke,1))
          ELSE IF (TRIM(KeyWord).eq.'Sout(idMtls)') THEN
            Npts=load_l(Nval, Cval, Ngrids, Sout(idMtls,1))
          ELSE IF (TRIM(KeyWord).eq.'NSTATION') THEN
            Npts=load_i(Nval, Rval, Ngrids, Nstation)
          ELSE IF (TRIM(KeyWord).eq.'POS') THEN
            DO ng=1,Ngrids
              allocate ( SCALARS(ng) % Sflag(Nstation(ng)) )            
              allocate ( SCALARS(ng) % SposX(Nstation(ng)) )            
              allocate ( SCALARS(ng) % SposY(Nstation(ng)) )            
            END DO
            is(1:Ngrids)=0
            DO WHILE (.true.)
              READ (inp,*,ERR=10) igrid, flag, Xpos, Ypos
              ng=MAX(1,ABS(igrid))
              is(ng)=is(ng)+1
              SCALARS(ng)%Sflag(is(ng))=flag
              SCALARS(ng)%SposX(is(ng))=Xpos
              SCALARS(ng)%SposY(is(ng))=Ypos
            END DO
 10         DO ng=1,Ngrids
              IF (Nstation(ng).ne.is(ng)) THEN
                WRITE (out,40) Nstation(ng), is(ng)
                STOP
              END IF
            END DO
          END IF
        END IF
      END DO
 20   WRITE (out,50) line
      STOP
 30   CLOSE (inp)
!
!-----------------------------------------------------------------------
!  Report input parameters.
!-----------------------------------------------------------------------
!
      IF (Lwrite) THEN
        DO ng=1,Ngrids
          IF (Lstations(ng)) THEN
            WRITE (out,60) ng
            WRITE (out,70) Nstation(ng), 'Nstation',                    &
     &            'Number of stations to write out into stations file.'
            IF (Sout(idFsur,ng)) WRITE (out,80) Sout(idFsur,ng),        &
     &          'Sout(idFsur)',                                         &
     &          'Write out free-surface.'
            IF (Sout(idUbar,ng)) WRITE (out,80) Sout(idUbar,ng),        &
     &          'Sout(idUbar)',                                         &
     &          'Write out 2D U-momentum component.'
            IF (Sout(idVbar,ng)) WRITE (out,80) Sout(idVbar,ng),        &
     &          'Sout(idVbar)',                                         &
     &          'Write out 2D V-momentum component.'
# ifdef SOLVE3D
            IF (Sout(idUvel,ng)) WRITE (out,80) Sout(idUvel,ng),        &
     &          'Sout(idUvel)',                                         &
     &          'Write out 3D U-momentum component.'
            IF (Sout(idVvel,ng)) WRITE (out,80) Sout(idVvel,ng),        &
     &          'Sout(idVvel)',                                         &
     &          'Write out 3D V-momentum component.'
            IF (Sout(idWvel,ng)) WRITE (out,80) Sout(idWvel,ng),        &
     &          'Sout(idWvel)',                                         &
     &          'Write out W-momentum component.'
            IF (Sout(idOvel,ng)) WRITE (out,80) Sout(idOvel,ng),        &
     &          'Sout(idOvel)',                                         &
     &          'Write out omega vertical velocity.'
            DO itrc=1,NT(ng)
              IF (Sout(idTvar(itrc),ng)) WRITE (out,90)                 &
     &            Sout(idTvar(itrc),ng), 'Hout(idTvar)',                &
     &            'Write out tracer ', itrc, TRIM(Vname(1,idTvar(itrc)))
            END DO
# endif
            IF (Sout(idUsms,ng)) WRITE (out,80) Sout(idUsms,ng),        &
     &          'Sout(idUsms)',                                         &
     &          'Write out surface U-momentum stress.'
            IF (Sout(idVsms,ng)) WRITE (out,80) Sout(idVsms,ng),        &
     &          'Sout(idVsms)',                                         &
     &          'Write out surface V-momentum stress.'
            IF (Sout(idVsms,ng)) WRITE (out,80) Sout(idVsms,ng),        &
     &          'Sout(idVsms)',                                         &
     &          'Write out bottom U-momentum stress.'
            IF (Sout(idVsms,ng)) WRITE (out,80) Sout(idVsms,ng),        &
     &          'Sout(idVsms)',                                         &
     &          'Write out bottom V-momentum stress.'
# ifdef BBL_MODL
            IF (Sout(idUbws,ng)) WRITE (out,80) Sout(idUbws,ng),        &
     &          'Sout(idUbws)',                                         &
     &          'Write out wind-induced, bottom U-wave stress.'
            IF (Sout(idVbws,ng)) WRITE (out,80) Sout(idVbws,ng),        &
     &          'Sout(idVbws)',                                         &
     &          'Write out wind-induced, bottom V-wave stress.'
            IF (Sout(idAbed,ng)) WRITE (out,80) Sout(idAbed,ng),        &
     &          'Sout(idAbed)',                                         &
     &          'Write out bed wave excursion amplitude.'
            IF (Sout(idUbed,ng)) WRITE (out,80) Sout(idUbed,ng),        &
     &          'Sout(idUbed)',                                         &
     &          'Write out bed wave orbital U-velocity.'
            IF (Sout(idVbed,ng)) WRITE (out,80) Sout(idVbed,ng),        &
     &          'Sout(idVbed)',                                         &
     &          'Write out bed wave orbital V-velocity.'
            IF (Sout(idUbot,ng)) WRITE (out,80) Sout(idUbot,ng),        &
     &          'Sout(idUbot)',                                         &
     &          'Write out bottom U-momentum above bed.'
            IF (Sout(idVbot,ng)) WRITE (out,80) Sout(idVbot,ng),        &
     &          'Sout(idVbot)',                                         &
     &          'Write out bottom V-momentum above bed.'
            IF (Sout(idHrip,ng)) WRITE (out,80) Sout(idHrip,ng),        &
     &          'Sout(idHrip)',                                         &
     &          'Write out bed ripple height.'
            IF (Sout(idLrip,ng)) WRITE (out,80) Sout(idLrip,ng),        &
     &          'Sout(idLrip)',                                         &
     &          'Write out bed ripple length.'
            IF (Sout(idZnot,ng)) WRITE (out,80) Sout(idZnot,ng),        &
     &          'Sout(idZnot)',                                         &
     &          'Write out bottom roughness.'
            IF (Sout(idZapp,ng)) WRITE (out,80) Sout(idZapp,ng),        &
     &          'Sout(idZapp)',                                         &
     &          'Write out apparent bottom roughness.'
# endif
# ifdef SOLVE3D
            IF (Sout(idTsur(itemp),ng)) WRITE (out,80)                  &
     &          Sout(idTsur(itemp),ng), 'Sout(idTsur)',                 &
     &          'Write out surface net heat flux.'
#  ifdef SHORTWAVE
            IF (Sout(idSrad,ng)) WRITE (out,80) Sout(idSrad,ng),        &
     &          'Sout(idSrad)',                                         &
     &          'Write out shortwave radiation flux.'
#  endif
#  ifdef BULK_FLUXES
            IF (Sout(idLrad,ng)) WRITE (out,80) Sout(idLrad,ng),        &
     &          'Sout(idLrad)',                                         &
     &          'Write out longwave radiation flux.'
            IF (Sout(idLhea,ng)) WRITE (out,80) Sout(idLhea,ng),        &
     &          'Sout(idLhea)',                                         &
     &          'Write out latent heat flux.'
            IF (Sout(idShea,ng)) WRITE (out,80) Sout(idShea,ng),        &
     &          'Sout(idShea)',                                         &
     &          'Write out sensible heat flux.'
#  endif
            IF (Sout(idDano,ng)) WRITE (out,80) Sout(idDano,ng),        &
     &          'Sout(idDano)',                                         &
     &          'Write out density anomaly.'
            IF (Sout(idVvis,ng)) WRITE (out,80) Sout(idVvis,ng),        &
     &          'Sout(idVvis)',                                         &
     &          'Write out vertical viscosity coefficient.'
            IF (Sout(idTdif,ng)) WRITE (out,80) Sout(idTdif,ng),        &
     &          'Sout(idTdif)',                                         &
     &          'Write out vertical T-diffusion coefficient.'
            IF (Sout(idSdif,ng)) WRITE (out,80) Sout(idSdif,ng),        &
     &          'Sout(idSdif)',                                         &
     &          'Write out vertical S-diffusion coefficient.'
#  ifdef LMD_SKPP
            IF (Sout(idHsbl,ng)) WRITE (out,80) Sout(idHsbl,ng),        &
     &          'Sout(idHsbl)',                                         &
     &          'Write out depth of surface boundary layer.'
#  endif
#  ifdef LMD_BKPP
            IF (Sout(idHbbl,ng)) WRITE (out,80) Sout(idHbbl,ng),        &
     &          'Sout(idHbbl)',                                         &
     &          'Write out depth of bottom boundary layer.'
#  endif
#  if defined GLS_MIXING || defined MY25_MIXING
            IF (Sout(idMtke,ng)) WRITE (out,80) Sout(idMtke,ng),        &
     &          'Sout(idMtke)',                                         &
     &          'Write out turbulent kinetic energy.'
            IF (Sout(idMtls,ng)) WRITE (out,80) Sout(idMtls,ng),        &
     &          'Sout(idMtls)',                                         &
     &          'Write out turbulent generic length-scale.'
#  endif
# endif
            WRITE (out,*)
            DO i=1,Nstation(ng)
              WRITE (out,100) i, SCALARS(ng)%Sflag(i),                  &
     &                           SCALARS(ng)%SposX(i),                  &
     &                           SCALARS(ng)%SposY(i)
            END DO
          END IF
        END DO
      END IF

  40  FORMAT (/,' READ_StaPar - Inconsistent number of stations, ',     &
     &        'Nstation = ',2i8,/,15x,'change input script values.')  
  50  FORMAT (/,' READ_StaPar - Error while processing line: ',/,a)

  60  FORMAT (/,/,' Stations Parameters, Grid: ',i2.2,                  &
     &        /,' =============================',/)  
  70  FORMAT (1x,i10,2x,a,t28,a)
  80  FORMAT (10x,l1,2x,a,t28,a)
  90  FORMAT (10x,l1,2x,a,t28,a,i2.2,':',1x,a)
 100  FORMAT (13x,'Flag and positions for station ',i4.4,':',           &
     &        i3,1x,2f10.4)

      RETURN
      END SUBROUTINE read_StaPar
#endif

      FUNCTION decode_line (line_text, KeyWord, Nval, Cval, Rval)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This function decodes lines of text from input script files.        !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
! Imported variable declarations.
!
      character (len=*), intent(in) :: line_text

      character (len=40), intent(inout) :: KeyWord

      integer, intent(inout) :: Nval

      character (len=80), dimension(100), intent(inout) :: Cval

      real(r8), dimension(100), intent(inout) :: Rval
!
! Local variable declarations
!
      logical :: IsString, Kextract, decode, nested
      integer :: Iblank, Icmm, Kstr, Kend, Linp
      integer :: Lend, LenS, Lstr, Lval, Nmul, Schar
      integer :: copies, i, ic, ie, is, j, status

      integer, dimension(20) :: Imul

      integer :: decode_line

      character (len=1 ), parameter :: blank = ' '
      character (len=80) :: Vstring, line, string
!
!------------------------------------------------------------------------
!  Decode input line.
!------------------------------------------------------------------------
!
!  Initialize.
!
      DO i=1,LEN(line)
        line(i:i)=blank
        Vstring(i:i)=blank
        string(i:i)=blank
      END DO
!
!  Get length of "line", remove leading and trailing blanks.
!
      line=TRIM(ADJUSTL(line_text))
      Linp=LEN_TRIM(line)
!
!  If not a blank or comment line [char(33)=!], decode and extract input
!  values.  Find equal sign [char(61)].
!
      status=-1
      nested=.false.
      IF ((Linp.gt.0).and.(line(1:1).ne.CHAR(33))) THEN
        status=1
        Kstr=1
        Kend=INDEX(line,CHAR(61),BACK=.false.)-1
        Lstr=INDEX(line,CHAR(61),BACK=.true.)+1
!
! Determine if KEYWORD is followed by double equal sign (==) indicating
! nested parameter.
!
        IF ((Lstr-Kend).eq.3) nested=.true.
!
! Extract KEYWORD, trim leading and trailing blanks.
!
        Kextract=.false.
        IF (Kend.gt.0) THEN
          Lend=Linp
          KeyWord=line(Kstr:Kend)
          Nval=0
          Kextract=.true.
        ELSE
          Lstr=1
          Lend=Linp
          Kextract=.true.
        END IF
!
! Extract parameter values string.  Remove comments [char(33)=!] or
! continuation symbol [char(92)=\], if any.  Trim leading trailing
! blanks.
!
        IF (Kextract) THEN
          Icmm=INDEX(line,CHAR(33),BACK=.false.)
          IF (Icmm.gt.0) Lend=Icmm-1
          Icmm=INDEX(line,CHAR(92),BACK=.false.)
          IF (Icmm.gt.0) Lend=Icmm-1
          Vstring=ADJUSTL(line(Lstr:Lend))
          Lval=LEN_TRIM(Vstring)
!
! The TITLE KEYWORD is a special one since it can include strings,
! numbers, spaces, and continuation symbol.
!
          IsString=.false.
          IF (TRIM(KeyWord).eq.'TITLE') THEN
            Nval=Nval+1
            Cval(Nval)=Vstring(1:Lval)
            IsString=.true.
          ELSE
!
! Check if there is a multiplication symbol [char(42)=*] in the variable
! string indicating repection of input values.
!
            Nmul=0
            DO i=1,Lval
              IF (Vstring(i:i).eq.CHAR(42)) THEN
                Nmul=Nmul+1           
                Imul(Nmul)=i
              END IF
            END DO
            ic=1
!
! Check for blank spaces [char(32)=' '] between entries and decode.
!
            is=1
            ie=Lval
            Iblank=0
            decode=.false.
            DO i=1,Lval
              IF (Vstring(i:i).eq.CHAR(32)) THEN
                IF (Vstring(i+1:i+1).ne.CHAR(32)) decode=.true.
                Iblank=i
              ELSE
                ie=i
              ENDIF                                    
              IF (decode.or.(i.eq.Lval)) THEN
                Nval=Nval+1
!
! Processing numeric values.  Check starting character to determine
! if numeric or character values. It is possible to have both when
! processing repetitions via the multiplication symbol.              
!
                Schar=ICHAR(Vstring(is:is))
                IF (((48.le.Schar).and.(Schar.le.57)).or.               &
     &              (Schar.eq.43).or.(Schar.eq.45)) THEN
                  IF ((Nmul.gt.0).and.                                  &
     &                (is.lt.Imul(ic)).and.(Imul(ic).lt.ie)) THEN
                    READ (Vstring(is:Imul(ic)-1),*) copies
                    Schar=ICHAR(Vstring(Imul(ic)+1:Imul(ic)+1))
                    IF ((43.le.Schar).and.(Schar.le.57)) THEN
                      READ (Vstring(Imul(ic)+1:ie),*) Rval(Nval)
                      DO j=1,copies-1
                        Rval(Nval+j)=Rval(Nval)
                      END DO
                    ELSE
                      string=Vstring(Imul(ic)+1:ie)
                      LenS=LEN_TRIM(string)
                      Cval(Nval)=string(1:LenS)
                      DO j=1,copies-1
                        Cval(Nval+j)=Cval(Nval)
                      END DO
                    END IF
                    Nval=Nval+copies-1
                    ic=ic+1
                  ELSE
                    string=Vstring(is:ie)
                    LenS=LEN_TRIM(string)
                    READ (string(1:LenS),*) Rval(Nval)
                  END IF
                ELSE
!
! Processing character values (logicals and strings).
!
                  IF ((Nmul.gt.0).and.                                  &
     &                (is.lt.Imul(ic)).and.(Imul(ic).lt.ie)) THEN
                    READ (Vstring(is:Imul(ic)-1),*) copies
                    Cval(Nval)=Vstring(Imul(ic)+1:ie)
                    DO j=1,copies-1
                      Cval(Nval+j)=Cval(Nval)
                    END DO
                    Nval=Nval+copies-1
                    ic=ic+1
                  ELSE
                    string=Vstring(is:ie)
                    Cval(Nval)=TRIM(ADJUSTL(string))
                  END IF
                  IsString=.true.
                END IF
                is=Iblank+1
                ie=Lval
                decode=.false.
              END IF
            END DO
          END IF
        END IF             
        status=Nval
      END IF
      decode_line=status       
      RETURN
      END FUNCTION decode_line

      FUNCTION load_i (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This function loads input values into a requested model integer     !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output integer variable.                              !
!     load_i     Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      real(r8), intent(in) :: Vinp(Ninp)
      integer, intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: load_i
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vout(i)=INT(Vinp(i))
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          Vout(i)=INT(Vinp(Ninp))
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vout(i)=INT(Vinp(i))
        END DO
      END IF
      load_i=ic

      RETURN
      END FUNCTION load_i

      FUNCTION load_l (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This function loads input values into a requested model logical     !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output integer variable.                              !
!     load_l     Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      character (len=80), intent(in) :: Vinp(Ninp)
      logical, intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: load_l
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.(Vinp(i)(1:1).eq.'t')) THEN
            Vout(i)=.true.
          ELSE
            Vout(i)=.false.
          END IF
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          IF ((Vinp(Ninp)(1:1).eq.'T').or.(Vinp(Ninp)(1:1).eq.'t')) THEN
            Vout(i)=.true.
          ELSE
            Vout(i)=.false.
          END IF
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          IF ((Vinp(i)(1:1).eq.'T').or.(Vinp(i)(1:1).eq.'t')) THEN
            Vout(i)=.true.
          ELSE
            Vout(i)=.false.
          END IF
        END DO
      END IF
      load_l=ic

      RETURN
      END FUNCTION load_l

      FUNCTION load_r (Ninp, Vinp, Nout, Vout)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This function loads input values into a requested model real        !
!  variable.                                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     Ninp       Size of input variable.                               !
!     Vinp       Input values                                          !
!     Nout       Number of output values.                              !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Vout       Output real variable.                                 !
!     load_r     Number of output values processed.                    !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: Ninp, Nout
      real(r8), intent(in) :: Vinp(Ninp)
      real(r8), intent(out) :: Vout(Nout)
!
!  Local variable declarations.
!
      integer :: i, ic
      integer :: load_r
!
!-----------------------------------------------------------------------
!  Load integer variable with input values.
!-----------------------------------------------------------------------
!
!  If not all values are provided for variable, assume the last value
!  for the rest of the array.
!
      ic=0
      IF (Ninp.le.Nout) THEN
        DO i=1,Ninp
          ic=ic+1
          Vout(i)=Vinp(i)
        END DO
        DO i=Ninp+1,Nout
          ic=ic+1
          Vout(i)=Vinp(Ninp)
        END DO
      ELSE
        DO i=1,Nout
          ic=ic+1
          Vout(i)=Vinp(i)
        END DO
      END IF
      load_r=ic

      RETURN
      END FUNCTION load_r
