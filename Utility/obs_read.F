#include "cppdefs.h"
#ifdef FOUR_DVAR
      SUBROUTINE obs_read (ng, model, backward)
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS adjoint group                          !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine reads in observations data when appropriate from    !
!  observations input NetCDF file.  The observations data is stored    !     
!  for use elsewhere.                                                  ! 
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
      logical :: readNLmod, readTLmod

      integer :: Irec, ObsSum, i, iobs, itrc, status

      integer, dimension(2) :: count, start
!
!---------------------------------------------------------------------
!  Read observation variables needed for interpolating the model
!  state at the observation locations.
!---------------------------------------------------------------------
!
      IF (ProcessObs(ng)) THEN
# if defined S4DVAR
        readNLmod=backward
        Irec=Ipass
# elif defined TLM_OBS
        readNLmod=.TRUE.
        readTLmod=.TRUE.
        Irec=1
# endif
!
!  Initialize observations processing counters.
!
        DO i=1,NstateVar(ng)
          FOURDVAR(ng)%ObsCount(i)=0
        END DO
        IF (backward) THEN
          ObsSurvey(ng)=ObsSurvey(ng)-1
        ELSE
          ObsSurvey(ng)=ObsSurvey(ng)+1
        END IF
!
!  Set number of observations to process.
!
        Nobs(ng)=FOURDVAR(ng)%NobsSurvey(ObsSurvey(ng))
!
!  Set number of datum to process at current time-step.
!
        IF (backward) THEN
          NendObs(ng)=NstrObs(ng)-1
          NstrObs(ng)=NstrObs(ng)-Nobs(ng)
        ELSE
          NstrObs(ng)=NendObs(ng)+1
          NendObs(ng)=NstrObs(ng)+Nobs(ng)-1
        END IF
        start(1)=NstrObs(ng)
        count(1)=Nobs(ng)
!
!  Read in observation type identifier.
!
        status=nf_get_vara_int(ncOBSid(ng), obsVid(idOtyp,ng),          &
     &                         start, count, ObsType)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOtyp)), ObsSurvey(ng),       &
     &                      TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcasti (ng, model, ObsType, Mobs)
# endif
!
!  Read in observation X-location (grid units).
!
        status=nf_get_vara_TYPE(ncOBSid(ng), obsVid(idObsX,ng),         &
     &                          start, count, Xobs)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idObsX)), ObsSurvey(ng),       &
     &                      TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, Xobs, Mobs)
# endif
!
!  Read in observation Y-location (grid units).
!
        status=nf_get_vara_TYPE(ncOBSid(ng), obsVid(idObsY,ng),         &
     &                          start, count, Yobs)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idObsY)), ObsSurvey(ng),       &
     &                      TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, Yobs, Mobs)
# endif

# ifdef SOLVE3D
!
!  Read in observation depth.  If negative, depth is meter. If greater
!  than zero, depth is model level.  The "Zobs" arrays is over-written
!  in "extract_obs3d" with grid units value in order to write grid
!  values into the observation NetCDF file.
!
        status=nf_get_vara_TYPE(ncOBSid(ng), obsVid(idObsD,ng),         &
     &                          start, count, Zobs)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idObsD)), ObsSurvey(ng),       &
     &                      TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
#  ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, Zobs, Mobs)
#  endif
# endif
!
!  Read in observation values.
!
        status=nf_get_vara_TYPE(ncOBSid(ng), obsVid(idOval,ng),         &
     &                          start, count, ObsVal)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOval)), ObsSurvey(ng),       &
     &                      TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, ObsVal, Mobs)
# endif
!
!  Read in observation error.
!
        status=nf_get_vara_TYPE(ncOBSid(ng), obsVid(idOerr,ng),         &
     &                          start, count, ObsErr)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vname(1,idOerr)), ObsSurvey(ng),       &
     &                      TRIM(OBSname(ng))
          exit_flag=2
          ioerror=status
          RETURN
        END IF
# ifdef DISTRIBUTE
        CALL mp_bcastf (ng, model, ObsErr, Mobs)
# endif
!
!  Read in nonlinear model values at observation locations.
!
        IF (readNLmod) THEN
          start(1)=Irec
          count(1)=1
          start(2)=NstrObs(ng)
          count(2)=Nobs(ng)
          status=nf_get_vara_TYPE(ncOBSid(ng), obsVid(idNLmo,ng),       &
     &                            start, count, NLmodVal(1,Irec))
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idNLmo)), ObsSurvey(ng),     &
     &                        TRIM(OBSname(ng))
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_bcastf (ng, model, NLmodVal(:,Irec), Mobs)
# endif
        END IF

# ifdef TLM_OBS
!
!  If adjoint pass and incremental S4DVAR, read in tangent linear model
!  values at observation locations.
!
        IF (readTLmod.and.(Ipass.ne.0)) THEN
          start(1)=Ipass
          count(1)=1
          start(2)=NstrObs(ng)
          count(2)=Nobs(ng)
          status=nf_get_vara_TYPE(ncOBSid(ng), obsVid(idTLmo,ng),       &
     &                            start, count, TLmodVal(1,Ipass))
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idTLmo)), ObsSurvey(ng),     &
     &                        TRIM(OBSname(ng))
            exit_flag=2
            ioerror=status
            RETURN
          END IF
!
!  Reset TLM values to zero at the first pass of the inner loop.
!
          IF ((inner.eq.1).and.(Ipass.eq.1)) THEN
            DO iobs=1,Mobs
              TLmodVal(iobs,1)=0.0_r8
              TLmodVal(iobs,2)=0.0_r8
            END DO
          END IF
#  ifdef DISTRIBUTE
          CALL mp_bcastf (ng, model, TLmodVal, 2*Mobs)
#  endif
        END IF
# endif
!
!  If adjoint pass, read in interpolation weights matrix.
!
        IF (readNLmod) THEN
          start(1)=1
          count(1)=Nweights
          start(2)=NstrObs(ng)
          count(2)=Nobs(ng)
          status=nf_get_vara_TYPE(ncOBSid(ng), obsVid(idModH,ng),       &
     &                            start, count, Hmat)
          IF (status.ne.nf_noerr) THEN
            WRITE (stdout,10) TRIM(Vname(1,idModH)), ObsSurvey(ng),     &
     &                        TRIM(OBSname(ng))
            exit_flag=2
            ioerror=status
            RETURN
          END IF
# ifdef DISTRIBUTE
          CALL mp_bcastf (ng, model, Hmat, Nweights*Mobs)
# endif
        END IF
!
!-----------------------------------------------------------------------
!  Set counters for number of observations to process for each state
!  variable.
!-----------------------------------------------------------------------
!
        DO iobs=1,Nobs(ng)
          IF  (ObsType(iobs).eq.isFsur) THEN
            FOURDVAR(ng)%ObsCount(isFsur)=                              &
     &                           FOURDVAR(ng)%ObsCount(isFsur)+1
          ELSE IF (ObsType(iobs).eq.isUbar) THEN
            FOURDVAR(ng)%ObsCount(isUbar)=                              &
     &                           FOURDVAR(ng)%ObsCount(isUbar)+1
          ELSE IF (ObsType(iobs).eq.isVbar) THEN
            FOURDVAR(ng)%ObsCount(isVbar)=                              &
     &                           FOURDVAR(ng)%ObsCount(isVbar)+1
# ifdef SOLVE3D
          ELSE IF (ObsType(iobs).eq.isUvel) THEN
            FOURDVAR(ng)%ObsCount(isUvel)=                              &
     &                           FOURDVAR(ng)%ObsCount(isUvel)+1
          ELSE IF (ObsType(iobs).eq.isVvel) THEN
            FOURDVAR(ng)%ObsCount(isVvel)=                              &
     &                           FOURDVAR(ng)%ObsCount(isVvel)+1
          ELSE
            DO itrc=1,NT(ng)
              IF (ObsType(iobs).eq.isTvar(itrc)) THEN
                i=isTvar(itrc)
                FOURDVAR(ng)%ObsCount(i)=FOURDVAR(ng)%ObsCount(i)+1
              END IF
            END DO
# endif
          END IF
        END DO
        IF (Master) THEN
          WRITE (stdout,20)
          ObsSum=0.0_r8
          DO i=1,NstateVar(ng)
            IF (FOURDVAR(ng)%ObsCount(i).gt.0) THEN
              WRITE (stdout,30) TRIM(Vname(1,idSvar(i))),               &
     &                          FOURDVAR(ng)%ObsCount(i)
              ObsSum=ObsSum+FOURDVAR(ng)%ObsCount(i)
            END IF
          END DO
          WRITE (stdout,30) 'TOTAL =', ObsSum
          WRITE (stdout,'(/)')
        END IF
!
!-----------------------------------------------------------------------
!  If applicable, set next observation survey time to process.
!-----------------------------------------------------------------------
!
        IF (backward) THEN
          IF ((ObsSurvey(ng)-1).ge.1) THEN
            ObsTime(ng)=FOURDVAR(ng)%SurveyTime(ObsSurvey(ng)-1)*day2sec
          END IF
        ELSE
          IF ((ObsSurvey(ng)+1).le.Nsurvey(ng)) THEN
            ObsTime(ng)=FOURDVAR(ng)%SurveyTime(ObsSurvey(ng)+1)*day2sec
          END IF
        END IF
      END IF
!
  10  FORMAT (/,' OBS_READ - error while reading variable: ',a,         &
     &        ' at survey record = ',i4,/,                              &
     &        11x,'in input NetCDF file: ',a)
  20  FORMAT (/,' Number of State Observations Processed:',/)
  30  FORMAT (42x,a,t52,i10)

      RETURN
      END SUBROUTINE obs_read
#else
      SUBROUTINE obs_read
      RETURN
      END SUBROUTINE obs_read
#endif
