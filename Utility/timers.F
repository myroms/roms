#include "cppdefs.h"
      SUBROUTINE wclock_on (ng, model, region)
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine turns on wall clock to meassure the elapsed time in    !
!  seconds spend by each parallel thread in requested model region.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_strings
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) ::  ng, model, region
!
!  Local variable declarations.
!
      integer :: iregion, MyThread

      integer :: my_getpid
#ifndef DISTRIBUTE
      integer :: my_threadnum
#endif

      real(r8), dimension(2) :: wtime

      real(r8) :: my_wtime
!
!-----------------------------------------------------------------------
! Initialize timing for all threads.
!-----------------------------------------------------------------------
!
      Cstr(region,model,ng)=my_wtime(wtime)
      IF ((region.eq.0).and.(proc(1,model,ng).eq.0)) THEN
        DO iregion=1,Nregion
          Cend(iregion,model,ng)=0.0_r8
          Csum(iregion,model,ng)=0.0_r8
        END DO
        proc(1,model,ng)=1
        proc(0,model,ng)=my_getpid()

!$OMP CRITICAL (START_WCLOCK)
#ifdef DISTRIBUTE
        CALL mp_barrier (ng)
        WRITE (stdout,10) ' Node #', MyRank,                            &
     &                    ' (pid=',proc(0,model,ng),') is active.'
        CALL my_flush (stdout)
#else
        MyThread=my_threadnum()
        WRITE (stdout,10) ' Thread #', MyThread,                        &
     &                    ' (pid=',proc(0,model,ng),') is active.'
#endif
 10     FORMAT (a,i3,a,i8,a)
        thread_count=thread_count+1
        IF (thread_count.eq.numthreads) thread_count=0
!$OMP END CRITICAL (START_WCLOCK)
      END IF
      RETURN
      END SUBROUTINE wclock_on
      SUBROUTINE wclock_off (ng, model, region)
!
!=======================================================================
!                                                                      !
!  This routine turns off wall clock to meassure the elapsed time in   !
!  seconds spend by each parallel thread in requested model region.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_strings
!    
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) ::  ng, model, region
!
!  Local variable declarations.
!
      integer :: iregion, MyThread

      integer :: my_threadnum

      real(r8) :: percent, sumcpu, sumper, total

      real(r8), dimension(2) :: wtime

      real(r8) :: my_wtime

#ifdef DISTRIBUTE
      character (len=3), dimension(0:Nregion) :: op_handle
#endif
!
!-----------------------------------------------------------------------
!  Compute elapsed wall time for all threads.
!-----------------------------------------------------------------------
!
      IF (region.ne.0) THEN
        Cend(region,model,ng)=Cend(region,model,ng)+                    &
     &                        (my_wtime(wtime)-Cstr(region,model,ng))
      END IF
!
!  Report elapsed wall time.
!
      IF ((region.eq.0).and.(proc(1,iNLM,ng).eq.1)) THEN
        Cend(region,iNLM,ng)=Cend(region,iNLM,ng)+                      &
     &                       (my_wtime(wtime)-Cstr(region,iNLM,ng))
        proc(1,iNLM,ng)=0
#ifdef TANGENT
        proc(1,iTLM,ng)=0
#endif
#ifdef ADJOINT
        proc(1,iADM,ng)=0
#endif

!$OMP CRITICAL (FINALIZE_WCLOCK)
!
!  Report total elapsed time (seconds) for each CPU.
!
#ifdef DISTRIBUTE
        CALL mp_barrier (ng)
        WRITE (stdout,10) ' Node   #', MyRank, ' CPU:',                 &
     &                    Cend(region,iNLM,ng)
        CALL my_flush (stdout)
#else
        MyThread=my_threadnum()
        WRITE (stdout,10) ' Thread #', MyThread, ' CPU:',               &
     &                    Cend(region,iNLM,ng)
#endif
 10     FORMAT (a,i3,a,f12.3)
!
! Report elapsed time profile for each region of the code.
!
        thread_count=thread_count+1
        Csum(region,iNLM,ng)=Csum(region,iNLM,ng)+                      &
     &                       Cend(region,iNLM,ng)
#ifdef PROFILE
        DO iregion=1,Nregion
          Csum(iregion,iNLM,ng)=Csum(iregion,iNLM,ng)+                  &
     &                          Cend(iregion,iNLM,ng)
# ifdef TANGENT
          Csum(iregion,iTLM,ng)=Csum(iregion,iTLM,ng)+                  &
     &                          Cend(iregion,iTLM,ng)
# endif
# ifdef ADJOINT
          Csum(iregion,iADM,ng)=Csum(iregion,iADM,ng)+                  &
     &                          Cend(iregion,iADM,ng)
# endif
        END DO
        Csum(35,iNLM,ng)=Csum(35,iNLM,ng)-Csum(15,iNLM,ng)
# ifdef TANGENT
        Csum(35,iTLM,ng)=Csum(35,iTLM,ng)-Csum(15,iTLM,ng)
# endif
# ifdef ADJOINT
        Csum(35,iADM,ng)=Csum(35,iADM,ng)-Csum(15,iADM,ng)
# endif
#endif
        IF (thread_count.eq.numthreads) THEN
          thread_count=0
#ifdef DISTRIBUTE
          op_handle(0:Nregion)='SUM'
          CALL mp_reduce (ng, iNLM, Nregion+1, Csum(0,iNLM,ng),         &
     &                    op_handle)
#endif
          IF (Master) THEN
            total=Csum(region,iNLM,ng)
#ifdef TANGENT
            total=total+Csum(region,iTLM,ng)
#endif
#ifdef ADJOINT
            total=total+Csum(region,iADM,ng)
#endif 
            WRITE (stdout,20) ' Total:', total
 20         FORMAT (a,8x,f14.3)
          END IF
#ifdef PROFILE
!
!  Report nonlinear model profile times.
!
          IF (Master) THEN
            WRITE (stdout,30) ' Nonlinear model elapsed time profile:'
 30         FORMAT (/,a,/)
          END IF
          sumcpu=0.0_r8
          sumper=0.0_r8
          Csum(35,iNLM,ng)=Csum(35,iNLM,ng)-Csum(15,iNLM,ng)
          DO iregion=1,37
            IF (Csum(iregion,iNLM,ng).gt.0.0_r8) THEN
              percent=100.0_r8*Csum(iregion,iNLM,ng)/                   &
     &                         Csum( region,iNLM,ng)
              IF (Master) WRITE (stdout,40) Pregion(iregion),           &
     &                                      Csum(iregion,iNLM,ng),      &
     &                                      percent
              sumcpu=sumcpu+Csum(iregion,iNLM,ng)
              sumper=sumper+percent
            END IF
          END DO
 40       FORMAT (2x,a,t53,f14.3,2x,'(',f7.4,' %)')
          IF (Master) WRITE (stdout,50) sumcpu, sumper
 50       FORMAT (t47,'Total:',f14.3,2x,f8.4)
# ifdef TANGENT
!
!  Report tangent linear model profile times.
!
          IF (Master) THEN
            WRITE (stdout,30)                                           &
     &            ' Tangent linear model elapsed time profile:'
          END IF
          sumcpu=0.0_r8
          sumper=0.0_r8
          Csum(35,iTLM,ng)=Csum(35,iTLM,ng)-Csum(15,iTLM,ng)
          DO iregion=1,37
            IF (Csum(iregion,iTLM,ng).gt.0.0_r8) THEN
              percent=100.0_r8*Csum(iregion,iTLM,ng)/                   &
     &                         Csum( region,iNLM,ng)
              IF (Master) WRITE (stdout,40) Pregion(iregion),           &
     &                                      Csum(iregion,iTLM,ng),      &
     &                                      percent
              sumcpu=sumcpu+Csum(iregion,iTLM,ng)
              sumper=sumper+percent
            END IF
          END DO
          IF (Master) WRITE (stdout,50) sumcpu, sumper
# endif
# ifdef ADJOINT
!
!  Report adjoint model profile times.
!
          IF (Master) THEN
            WRITE (stdout,30) ' Adjoint model elapsed time profile:'
          END IF
          sumcpu=0.0_r8
          sumper=0.0_r8
          Csum(35,iADM,ng)=Csum(35,iADM,ng)-Csum(15,iADM,ng)
          DO iregion=1,37
            IF (Csum(iregion,iADM,ng).gt.0.0_r8) THEN
              percent=100.0_r8*Csum(iregion,iADM,ng)/                   &
     &                         Csum( region,iNLM,ng)
              IF (Master) WRITE (stdout,40) Pregion(iregion),           &
     &                                      Csum(iregion,iADM,ng),      &
     &                                      percent
              sumcpu=sumcpu+Csum(iregion,iADM,ng)
              sumper=sumper+percent
            END IF
          END DO
          IF (Master) WRITE (stdout,50) sumcpu, sumper
# endif
# ifdef DISTRIBUTE
!
!  Report nonlinear model elapsed time for message passage
!  communications.
!
          IF (Master) THEN
            WRITE (stdout,30) ' Nonlinear model message Passage profile:'
          END IF
          sumcpu=0.0_r8
          sumper=0.0_r8
          DO iregion=38,Nregion
            IF (Csum(iregion,iNLM,ng).gt.0.0_r8) THEN
              percent=100.0_r8*Csum(iregion,iNLM,ng)/                   &
     &                         Csum( region,iNLM,ng)
              IF (Master) WRITE (stdout,40) Pregion(iregion),           &
     &                                      Csum(iregion,iNLM,ng),      &
     &                                      percent
              sumcpu=sumcpu+Csum(iregion,iNLM,ng)
              sumper=sumper+percent
            END IF
          END DO
          IF (Master) WRITE (stdout,50) sumcpu, sumper
#  ifdef TANGENT
!
!  Report tangent linear model elapsed time for message passage
!  communications.
!
          IF (Master) THEN
            WRITE (stdout,30)                                           &
     &            ' Tangent linear model message Passage profile:'
          END IF
          sumcpu=0.0_r8
          sumper=0.0_r8
          DO iregion=38,Nregion
            IF (Csum(iregion,iTLM,ng).gt.0.0_r8) THEN
              percent=100.0_r8*Csum(iregion,iTLM,ng)/                   &
     &                         Csum( region,iNLM,ng)
              IF (Master) WRITE (stdout,40) Pregion(iregion),           &
     &                                      Csum(iregion,iTLM,ng),      &
     &                                      percent
              sumcpu=sumcpu+Csum(iregion,iTLM,ng)
              sumper=sumper+percent
            END IF
          END DO
          IF (Master) WRITE (stdout,50) sumcpu, sumper
#  endif
#  ifdef ADJOINT
!
!  Report adjoint model elapsed time for message passage
!  communications.
!
          IF (Master) THEN
            WRITE (stdout,30) 'Adjoint model message Passage profile:'
          END IF
          sumcpu=0.0_r8
          sumper=0.0_r8
          DO iregion=38,Nregion
            IF (Csum(iregion,iADM,ng).gt.0.0_r8) THEN
              percent=100.0_r8*Csum(iregion,iADM,ng)/                   &
     &                         Csum( region,iNLM,ng)
              IF (Master) WRITE (stdout,40) Pregion(iregion),           &
     &                                      Csum(iregion,iADM,ng),      &
     &                                      percent
              sumcpu=sumcpu+Csum(iregion,iADM,ng)
              sumper=sumper+percent
            END IF
          END DO
          IF (Master) WRITE (stdout,50) sumcpu, sumper
#  endif

# endif
#endif
        END IF
!$OMP END CRITICAL (FINALIZE_WCLOCK)
      END IF
      RETURN
      END SUBROUTINE wclock_off
