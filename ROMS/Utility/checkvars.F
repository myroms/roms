#include "cppdefs.h"
      SUBROUTINE checkvars (ng, model, ncname, string, Nrec, Nvar,      &
     &                      tvarnam, get_var, have_var)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine checks if needed state variables are available in      !
!  requested NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     model      Calling model identifier.                             !
!     ncname     NetCDF file name.                                     !
!     string     Identification string.                                !
!     Nrec       Number of time records available.                     !
!     Nvar       Size of logical switches arrays.                      !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     tvarnam    Name of time record variable.                         !
!     get_var    Logical switches (T/F), in terms of variable ID,      !
!                  indicating state variables needed by the model.     !
!     have_var   Logical switches (T/F), in terms of variable ID,      !
!                  indicating state variables available in NetCDF      !
!                  file.                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
#if defined SEDIMENT || defined BBL_MODEL
      USE mod_sediment
#endif
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_bcastl
#endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, Nvar
      integer, intent(inout) :: Nrec

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: string

      character (len=*), intent(inout) :: tvarnam

      logical, dimension(Nvar), intent(out) :: get_var
      logical, dimension(Nvar), intent(out) :: have_var
!
!  Local variable declarations.
!
      integer :: i, itrc, nvd

      integer, dimension(4) :: Vsize

      character (len=80) :: fname
!
!-----------------------------------------------------------------------
!  Determine state variables needed and check if they are available in
!  requested NetCDF file.
!-----------------------------------------------------------------------
!
!  Inquire about the contents of input NetCDF file:  Inquire about
!  the dimensions and variables.  Check for consistency.
!
      IF (InpThread) THEN
        CALL opencdf (ng, 1, ncname, fname, N(ng), 0, Nrec, nvd, Vsize)
        IF (exit_flag.ne.NoError) RETURN
      END IF
!
!  Initialize logical switches.
!
      DO i=1,Nvar
        get_var(i)=.FALSE.
        have_var(i)=.FALSE.
      END DO
!
!  Determine state variables to read from input NetCDF file.
!
#ifdef ANA_INITIAL
      IF (nrrec(ng).ne.0) THEN
        get_var(idFsur)=.TRUE.
        get_var(idUbar)=.TRUE.
        get_var(idVbar)=.TRUE.
# ifdef SOLVE3D
        get_var(idUvel)=.TRUE.
        get_var(idVvel)=.TRUE.
        DO itrc=1,NAT
          get_var(idTvar(itrc))=.TRUE.
        END DO
# endif
      END IF
#else
      get_var(idFsur)=.TRUE.
      get_var(idUbar)=.TRUE.
      get_var(idVbar)=.TRUE.
# ifdef SOLVE3D
      get_var(idUvel)=.TRUE.
      get_var(idVvel)=.TRUE.
      DO itrc=1,NAT
        get_var(idTvar(itrc))=.TRUE.
      END DO
# endif
#endif
#ifdef SOLVE3D
# if defined BIOLOGY
#  if defined ANA_BIOLOGY
      IF (nrrec(ng).ne.0) THEN
        DO itrc=1,NBT
          get_var(idTvar(idbio(itrc)))=.TRUE.
        END DO
      END IF
#  else
      DO itrc=1,NBT
        get_var(idTvar(idbio(itrc)))=.TRUE.
      END DO
#  endif
# endif
# if defined T_PASSIVE
#  if defined ANA_PASSIVE
      IF (nrrec(ng).ne.0) THEN
        DO itrc=1,NPT
          get_var(idTvar(inert(itrc)))=.TRUE.
        END DO
      END IF
#  else
      DO itrc=1,NPT
        get_var(idTvar(inert(itrc)))=.TRUE.
      END DO
#  endif
# endif
# ifdef SEDIMENT
#  ifdef ANA_SEDIMENT
      IF (nrrec(ng).ne.0) THEN
        DO itrc=1,NST
          get_var(idTvar(idsed(itrc)))=.TRUE.
          get_var(idfrac(itrc))=.TRUE.
          get_var(idBmas(itrc))=.TRUE.
        END DO
        DO itrc=1,MBEDP
          get_var(idSbed(itrc))=.TRUE.
        END DO
      END IF
#  else
      DO itrc=1,NST
        get_var(idTvar(idsed(itrc)))=.TRUE.
        get_var(idfrac(itrc))=.TRUE.
        get_var(idBmas(itrc))=.TRUE.
      END DO
      DO itrc=1,MBEDP
        get_var(idSbed(itrc))=.TRUE.
      END DO
#  endif
#  ifdef BEDLOAD
      IF (nrrec(ng).ne.0) THEN
        DO itrc=1,NST
          get_var(idUbld(itrc))=.TRUE.
          get_var(idVbld(itrc))=.TRUE.
        END DO
      END IF
#  endif
#  ifdef SED_MORPH
      IF (nrrec(ng).ne.0) THEN
        get_var(idbath)=.TRUE.
      END IF
#  endif
# endif 
# if defined SEDIMENT || defined BBL_MODEL
      DO itrc=1,MBOTP
        get_var(idBott(itrc))=.TRUE.
      END DO
# endif
#endif
#ifdef PERFECT_RESTART
!
!  Determine perfect restart fields to read from input NetCDF file.
!
       IF ((model.eq.0).and.(nrrec(ng).ne.0)) THEN
# ifdef SOLVE3D
         get_var(idRu3d)=.TRUE.
         get_var(idRv3d)=.TRUE.
# endif
         get_var(idRzet)=.TRUE.
         get_var(idRu2d)=.TRUE.
         get_var(idRv2d)=.TRUE.
# if defined GLS_MIXING || defined MY25_MIXING
         get_var(idMtke)=.TRUE.
         get_var(idMtls)=.TRUE.
         get_var(idVmLS)=.TRUE.
         get_var(idVmKK)=.TRUE.
#  ifdef GLS_MIXING
         get_var(idVmKP)=.TRUE.
#  endif
# endif
       END IF
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
!
!  Determine aditional 4Dvar state varaibles to adjust.
!
       IF ((model.eq.iADM).or.(model.eq.iTLM).or.                       &
     &     (model.eq.5).or.(model.eq.6)) THEN
# ifdef ADJUST_STFLUX
         DO itrc=1,NT(ng)
           get_var(idTsur(itrc))=.TRUE.
         END DO
# endif
# ifdef ADJUST_WSTRESS
         get_var(idUsms)=.TRUE.
         get_var(idVsms)=.TRUE.
# endif
       END IF
#endif
!
!  Scan variable list from input NetCDF and activate switches for
!  model state variables.
!
      IF (InpThread) THEN
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idtime))) THEN
            tvarnam=TRIM(varnam(i))
            have_var(idtime)=.TRUE.
          END IF
#if defined SEDIMENT && defined SED_MORPH
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idbath))) THEN
            have_var(idbath)=.TRUE.
          END IF
#endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idFsur))) THEN
            have_var(idFsur)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRzet))) THEN
            have_var(idRzet)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbar))) THEN
            have_var(idUbar)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRu2d))) THEN
            have_var(idRu2d)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbar))) THEN
            have_var(idVbar)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRv2d))) THEN
            have_var(idRv2d)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idKhor))) THEN
            have_var(idKhor)=.TRUE.
          END IF
#ifdef ADJUST_WSTRESS
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUsms))) THEN
            have_var(idUsms)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVsms))) THEN
            have_var(idVsms)=.TRUE.
          END IF
#endif
#ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idKver))) THEN
            have_var(idKver)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUvel))) THEN
            have_var(idUvel)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRu3d))) THEN
            have_var(idRu3d)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvel))) THEN
            have_var(idVvel)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idRv3d))) THEN
            have_var(idRv3d)=.TRUE.
          END IF
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTvar(itrc)))) THEN
              have_var(idTvar(itrc))=.TRUE.
            END IF
          END DO
# ifdef ADJUST_STFLUX
          DO itrc=1,NT(ng)
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTsur(itrc)))) THEN
              have_var(idTsur(itrc))=.TRUE.
            END IF
          END DO
# endif
# if defined GLS_MIXING || defined MY25_MIXING || defined LMD_MIXING
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVvis))) THEN
            have_var(idVvis)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTdif))) THEN
            have_var(idTdif)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSdif))) THEN
            have_var(idSdif)=.TRUE.
          END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMtke))) THEN
            have_var(idMtke)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idMtls))) THEN
            have_var(idMtls)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVmLS))) THEN
            have_var(idVmLS)=.TRUE.
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVmKK))) THEN
            have_var(idVmKK)=.TRUE.
          END IF
#  ifdef GLS_MIXING
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVmKP))) THEN
            have_var(idVmKP)=.TRUE.
          END IF
#  endif
# endif
# ifdef SEDIMENT
          DO itrc=1,NST
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idfrac(itrc)))) THEN
              have_var(idfrac(itrc))=.TRUE.
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idBmas(itrc)))) THEN
              have_var(idBmas(itrc))=.TRUE.
            END IF
          END DO
          DO itrc=1,MBEDP
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idSbed(itrc)))) THEN
              have_var(idSbed(itrc))=.TRUE.
            END IF
          END DO
#  ifdef BEDLOAD
          DO itrc=1,NST
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idUbld(itrc)))) THEN
              have_var(idUbld(itrc))=.TRUE.
            END IF
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idVbld(itrc)))) THEN
              have_var(idVbld(itrc))=.TRUE.
            END IF
          END DO
#  endif
# endif
# if defined SEDIMENT || defined BBL_MODEL
          DO itrc=1,MBOTP
            IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idBott(itrc)))) THEN
              have_var(idBott(itrc))=.TRUE.
            END IF
          END DO
# endif
#endif
        END DO
!
!  Check if model state variables are available in input NetCDF file.
!
        IF (.not.have_var(idtime)) THEN
          WRITE (stdout,10) string, TRIM(tvarnam), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
#if defined SEDIMENT && defined SED_MORPH
        IF (.not.have_var(idbath).and.get_var(idbath)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idbath)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
#endif
        IF (.not.have_var(idFsur).and.get_var(idFsur)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idFsur)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idRzet).and.get_var(idRzet)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idRzet)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idUbar).and.get_var(idUbar)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idUbar)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idRu2d).and.get_var(idRu2d)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idRu2d)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idVbar).and.get_var(idVbar)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idVbar)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idRv2d).and.get_var(idRv2d)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idRv2d)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
#ifdef ADJUST_WSTRESS
        IF ((model.eq.5).or.(model.eq.6)) THEN
          IF (.not.have_var(idUsms).and.get_var(idUsms)) THEN
            WRITE (stdout,10) string, TRIM(Vname(1,idUsms)),            &
     &                        TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
          IF (.not.have_var(idVsms).and.get_var(idVsms)) THEN
            WRITE (stdout,10) string, TRIM(Vname(1,idVsms)),            &
     &                        TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
        END IF
#endif
#ifdef SOLVE3D
        IF (.not.have_var(idUvel).and.get_var(idUvel)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idUvel)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idRu3d).and.get_var(idRu3d)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idRu3d)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idVvel).and.get_var(idVvel)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idVvel)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idRv3d).and.get_var(idRv3d)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idRv3d)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        DO itrc=1,NT(ng)
          IF (.not.have_var(idTvar(itrc)).and.                          &
     &        get_var(idTvar(itrc))) THEN
            WRITE (stdout,10) string, TRIM(Vname(1,idTvar(itrc))),      &
     &                        TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
        END DO
# ifdef ADJUST_STFLUX
        IF ((model.eq.5).or.(model.eq.6)) THEN
          DO itrc=1,NT(ng)
            IF (.not.have_var(idTsur(itrc)).and.                        &
     &          get_var(idTsur(itrc))) THEN
              WRITE (stdout,10) string, TRIM(Vname(1,idTsur(itrc))),    &
     &                          TRIM(ncname)
              exit_flag=2
              RETURN
            END IF
          END DO
        END IF
# endif
# if defined GLS_MIXING || defined MY25_MIXING
        IF (.not.have_var(idMtke).and.get_var(idMtke)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idMtke)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idMtls).and.get_var(idMtls)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idMtls)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idVmLS).and.get_var(idVmLS)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idVmLS)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
        IF (.not.have_var(idVmKK).and.get_var(idVmKK)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idVmKK)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
#  ifdef GLS_MIXING
        IF (.not.have_var(idVmKP).and.get_var(idVmKP)) THEN
          WRITE (stdout,10) string, TRIM(Vname(1,idVmKP)), TRIM(ncname)
          exit_flag=2
          RETURN
        END IF
#  endif
# endif
# ifdef SEDIMENT
        DO itrc=1,NST
          IF (.not.have_var(idfrac(itrc)).and.                          &
     &        get_var(idfrac(itrc))) THEN
            WRITE (stdout,10) string, TRIM(Vname(1,idfrac(itrc))),      &
     &                        TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
          IF (.not.have_var(idBmas(itrc)).and.                          &
     &        get_var(idBmas(itrc))) THEN
            WRITE (stdout,10) string, TRIM(Vname(1,idBmas(itrc))),      &
     &                        TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
        END DO
        DO itrc=1,MBEDP
          IF (.not.have_var(idSbed(itrc)).and.                          &
     &        get_var(idSbed(itrc))) THEN
            WRITE (stdout,10) string, TRIM(Vname(1,idSbed(itrc))),      &
     &                        TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
        END DO
#  ifdef BEDLOAD
        DO itrc=1,NST
          IF (.not.have_var(idUbld(itrc)).and.                          &
     &        get_var(idUbld(itrc))) THEN
            WRITE (stdout,10) string, TRIM(Vname(1,idUbld(itrc))),      &
     &                        TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
          IF (.not.have_var(idVbld(itrc)).and.                          &
     &        get_var(idVbld(itrc))) THEN
            WRITE (stdout,10) string, TRIM(Vname(1,idVbld(itrc))),      &
     &                        TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
        END DO
#  endif
# endif
# if defined SEDIMENT || defined BBL_MODEL
!
! Check only for the first four required properties, but we looked for
! MBOTP properties above.
!
        DO itrc=1,4
          IF (.not.have_var(idBott(itrc)).and.                          &
     &        get_var(idBott(itrc))) THEN
            WRITE (stdout,10) string, TRIM(Vname(1,idBott(itrc))),      &
     &                        TRIM(ncname)
            exit_flag=2
            RETURN
          END IF
        END DO
# endif
#endif
      END IF

#ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Broadcast variable information to all nodes.
!-----------------------------------------------------------------------
!
      CALL mp_bcastl (ng, model, have_var, Nvar)
      CALL mp_bcasti (ng, model, Nrec, 1)
#endif

  10  FORMAT (/,a,'GET_VARS  - unable to find model variable: ',a,      &
     &        /,17x,'in input NetCDF file: ',a)

      RETURN
      END SUBROUTINE checkvars
