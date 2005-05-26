#include "cppdefs.h"
#ifdef FOUR_DVAR
      SUBROUTINE obs_initial (ng, model, backward)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS adjoint Group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine opens and reads in observations  NetCDF and sets    !
!  various variables needed for processing of the state solution at    !
!  observations locations during variational data assimilation.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_fourdvar
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!    
!  Imported variable declarations.
!
      logical, intent(in) :: backward

      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      logical, dimension(NV) :: got_var(NV)

      integer :: ifield = 0
      integer :: Ifirst, i, recdim, status

      real(r8) :: tend

      character (len=80) :: fname
!
!-----------------------------------------------------------------------
!  Inquire about the contents of observation NetCDF file:  Inquire about
!  the dimensions and variables.
!-----------------------------------------------------------------------
!
      IF (ncOBSid(ng).eq.-1) THEN
        ifield=0
        CALL opencdf (ng, 1, OBSname(ng), fname, N(ng), ifield, recdim)
        IF (exit_flag.ne.NoError) RETURN
!
!  Initialize logical switches.
!
        DO i=1,NV
          got_var(i)=.FALSE.
        END DO
!
!  Scan variable list from observation NetCDF and activate switches for
!  required variables.
!
        DO i=1,nvars
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idOday))) THEN
            got_var(idOday)=.TRUE.
            obsVid(idOday,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idNobs))) THEN
            got_var(idNobs)=.TRUE.
            obsVid(idNobs,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idOtyp))) THEN
            got_var(idOtyp)=.TRUE.
            obsVid(idOtyp,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idObsT))) THEN
            got_var(idObsT)=.TRUE.
            obsVid(idObsT,ng)=i
          END IF
#  ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idObsD))) THEN
            got_var(idObsD)=.TRUE.
            obsVid(idObsD,ng)=i
          END IF
#  endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idObsX))) THEN
            got_var(idObsX)=.TRUE.
            obsVid(idObsX,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idObsY))) THEN
            got_var(idObsY)=.TRUE.
            obsVid(idObsY,ng)=i
          END IF
#  ifdef SOLVE3D
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idObsZ))) THEN
            got_var(idObsZ)=.TRUE.
            obsVid(idObsZ,ng)=i
          END IF
#  endif
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idOvar))) THEN
            got_var(idOvar)=.TRUE.
            obsVid(idOvar,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idOerr))) THEN
            got_var(idOerr)=.TRUE.
            obsVid(idOerr,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idOval))) THEN
            got_var(idOval)=.TRUE.
            obsVid(idOval,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idNLmo))) THEN
            got_var(idNLmo)=.TRUE.
            obsVid(idNLmo,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idTLmo))) THEN
            got_var(idTLmo)=.TRUE.
            obsVid(idTLmo,ng)=i
          END IF
          IF (TRIM(varnam(i)).eq.TRIM(Vname(1,idModH))) THEN
            got_var(idModH)=.TRUE.
            obsVid(idModH,ng)=i
          END IF
        END DO
!
!  Check if needed obsrvation variables are available.
!
        IF (.not.got_var(idOday)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOday)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idNobs)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idNobs)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOtyp)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOtyp)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsT)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idObsT)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idObsD)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idObsD)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
#  endif
        IF (.not.got_var(idObsX)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idObsX)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idObsY)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idObsY)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
#  ifdef SOLVE3D
        IF (.not.got_var(idObsZ)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idObsZ)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
#  endif
        IF (.not.got_var(idOvar)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOvar)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOerr)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOerr)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idOval)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOval)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
        IF (.not.got_var(idNLmo)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idNLmo)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
#  ifdef IS4DVAR
        IF (.not.got_var(idTLmo)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idTLmo)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
#  endif
        IF (.not.got_var(idModH)) THEN
          WRITE (stdout,10) TRIM(Vname(1,idModH)), TRIM(OBSname(ng))
          exit_flag=2
          RETURN
        END IF
!
!  Open observations NetCDF file.
!
        status=nf_open(TRIM(OBSname(ng)), nf_write, ncOBSid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20) TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Read various observation variables from NetCDF file.
!-----------------------------------------------------------------------
!
!  Read in the time of each observation survey.
!
      status=nf_get_vara_TYPE(ncOBSid(ng), obsVid(idOday,ng), 1,        &
     &                        Nsurvey(ng),                              &
     &                        FOURDVAR(ng)%SurveyTime)
      IF (status.ne.nf_noerr) THEN
        WRITE (stdout,30) TRIM(Vname(1,idOday)), TRIM(OBSname(ng))
        exit_flag=2
        ioerror=status
        RETURN
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastf (ng, model, FOURDVAR(ng)%SurveyTime, Nsurvey(ng))
# endif
!
!  Read in observations global variance.
!
      IF (got_var(idOvar)) THEN
        status=nf_get_vara_TYPE(ncOBSid(ng), obsVid(idOvar,ng), 1,      &
     &                          NstateVar(ng),                          &
     &                          FOURDVAR(ng)%ObsVar)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,30) TRIM(Vname(1,idOvar)), TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
      ELSE
        DO i=1,NstateVar(ng)
          FOURDVAR(ng)%ObsVar(i)=1.0_r8
        END DO
      END IF
# ifdef DISTRIBUTE
      CALL mp_bcastf (ng, model, FOURDVAR(ng)%ObsVar, NstateVar(ng))
# endif
!
!  Determine if there is any data available in the model time
!  window.  Set first survey record to process.
!
      IF (backward) THEN
        Ifirst=0
        tend=(time(ng)-(ntstart-1)*dt(ng))*sec2day
        DO i=1,Nsurvey(ng)
          IF ((FOURDVAR(ng)%SurveyTime(i).ge.tend).and.                 &
     &        (FOURDVAR(ng)%SurveyTime(i).le.tdays(ng))) THEN
            Ifirst=MAX(Ifirst,i)
          END IF            
        END DO
        IF (Ifirst.eq.0) THEN
          WRITE (stdout,40) tend, tdays(ng)
          STOP
        END IF
      ELSE
        Ifirst=1
        tend=(time(ng)+ntend*dt(ng))*sec2day
        DO i=1,Nsurvey(ng)
          IF ((FOURDVAR(ng)%SurveyTime(i).ge.tdays(ng)).and.            &
     &        (FOURDVAR(ng)%SurveyTime(i).le.tend)) THEN
            Ifirst=MIN(Ifirst,i)
          END IF            
        END DO
        IF (Ifirst.eq.0) THEN
          WRITE (stdout,40) tdays(ng), tend
          STOP
        END IF
      END IF
      ObsTime(ng)=FOURDVAR(ng)%SurveyTime(Ifirst)*day2sec
!
!  Initialize observation survey counter.  This is the counter of data
!  assimilation cycles or observations survey times on which the model
!  state is extracted (interpolated) at the observation locations.
!
      IF (backward) THEN
        ObsSurvey(ng)=Ifirst+1
      ELSE
        ObsSurvey(ng)=Ifirst-1
      END IF
!
!  Initialize time swith to process model state at observation
!  locations.
!
      ProcessObs(ng)=.FALSE.
!
!  Initialize cost function misfit between model and observations.
!  The IF statement is to avoid rewritting its value before it is
!  written into the initial NetCDF file.
!
      IF (.not.backward.and.(Ipass.eq.1)) THEN
        FOURDVAR(ng)%BackFun(0)=0.0_r8
        FOURDVAR(ng)%CostFun(0)=0.0_r8
        DO i=1,NstateVar(ng)
          FOURDVAR(ng)%BackFun(i)=0.0_r8
          FOURDVAR(ng)%CostFun(i)=0.0_r8
        END DO
      END IF
!
!  Initialize descent algorithm step size factors.
!
      StepSizeB=0.0_r8
      StepSizeT=0.0_r8
!
!  Set staring and ending observation indices.
!
      IF (backward) THEN
        NstrObs(ng)=0
        DO i=1,Ifirst
          NstrObs(ng)=NstrObs(ng)+FOURDVAR(ng)%NobsSurvey(Ifirst)
        END DO
        NstrObs(ng)=NstrObs(ng)+1
        NendObs(ng)=0
      ELSE
        IF (Ifirst.eq.1) THEN
          NstrObs(ng)=0
        ELSE
          NstrObs(ng)=FOURDVAR(ng)%NobsSurvey(Ifirst)-1
        END IF
        NendObs(ng)=0
      END IF
!
  10  FORMAT (/,' OBS_INITIAL - unable to find model variable: ',a,     &
     &        /,12x,'in input NetCDF file: ',a)
  20  FORMAT (/,' OBS_INITIAL - unable to open input NetCDF file: ',a)
  30  FORMAT (/,' OBS_INITIAL - error while reading variable: ',a,2x,   &
     &        /,12x,'in input NetCDF file: ',a)
  40  FORMAT (/,' OBS_INITIAL - No are observations available for',     &
     &          ' time window (days): ',/,12x,f12.4,' - ',f12.4,/)
      RETURN
      END SUBROUTINE obs_initial
#else
      SUBROUTINE obs_initial
      RETURN
      END SUBROUTINE obs_initial
#endif
