#include "cppdefs.h"
      SUBROUTINE wclock_on (ng, region)
!
!========================================== Alexander F. Shchepetkin ===
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
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
      integer, intent(in) ::  ng, region
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
      Cstr(region,ng)=my_wtime(wtime)
      IF ((region.eq.0).and.(proc(1,ng).eq.0)) THEN
        DO iregion=1,Nregion
          Cend(iregion,ng)=0.0_r8
          Csum(iregion,ng)=0.0_r8
        END DO
        proc(1,ng)=1
        proc(0,ng)=my_getpid()
!$OMP CRITICAL (START_WCLOCK)
#ifdef DISTRIBUTE
        CALL mp_barrier (ng)
        WRITE (stdout,10) ' Node #', MyRank,                            &
     &                    ' (pid=',proc(0,ng),') is active.'
        CALL my_flush (stdout)
#else
        MyThread=my_threadnum()
        WRITE (stdout,10) ' Thread #', MyThread,                        &
     &                    ' (pid=',proc(0,ng),') is active.'
#endif
 10     FORMAT (a,i3,a,i8,a)
        thread_count=thread_count+1
        IF (thread_count.eq.numthreads) thread_count=0
!$OMP END CRITICAL (START_WCLOCK)
      END IF
      RETURN
      END SUBROUTINE wclock_on
      SUBROUTINE wclock_off (ng, region)
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
      integer, intent(in) ::  ng, region
!
!  Local variable declarations.
!
      integer :: iregion, MyThread

      integer :: my_threadnum

      real(r8) :: percent, sumcpu, sumper

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
      Cend(region,ng)=Cend(region,ng)+(my_wtime(wtime)-Cstr(region,ng))
!
!  Report elapsed wall time.
!
      IF ((region.eq.0).and.(proc(1,ng).eq.1)) THEN
        proc(1,ng)=0
!$OMP CRITICAL (FINALIZE_WCLOCK)
!
!  Report total elapsed time (seconds) for each CPU.
!
#ifdef DISTRIBUTE
        CALL mp_barrier (ng)
        WRITE (stdout,10) ' Node   #',MyRank,' CPU:',Cend(region,ng)
        CALL my_flush (stdout)
#else
        MyThread=my_threadnum()
        WRITE (stdout,10) ' Thread #',MyThread,' CPU:',Cend(region,ng)
#endif
 10     FORMAT (a,i3,a,f12.3)
!
! Report elapsed time profile for each region of the code.
!
        thread_count=thread_count+1
        Csum(region,ng)=Csum(region,ng)+Cend(region,ng)
#ifdef PROFILE
        DO iregion=1,Nregion
          Csum(iregion,ng)=Csum(iregion,ng)+Cend(iregion,ng)
        END DO
        Csum(35,ng)=Csum(35,ng)-Csum(15,ng)
#endif
        IF (thread_count.eq.numthreads) THEN
          thread_count=0
#ifdef DISTRIBUTE
          op_handle(0:Nregion)='SUM'
          CALL mp_reduce (ng, Nregion+1, Csum(0,ng), op_handle)
#endif
          IF (Master) WRITE (stdout,20) ' Total:', Csum(region,ng)
 20       FORMAT (a,8x,f14.3)
#ifdef PROFILE
          IF (Master) WRITE (stdout,30) ' Elapsed time profile:'
 30       FORMAT (/,a,/)
          sumcpu=0.0_r8
          sumper=0.0_r8
          Csum(35,ng)=Csum(35,ng)-Csum(15,ng)
          DO iregion=1,37
            IF (Csum(iregion,ng).gt.0.0_r8) THEN
              percent=100.0_r8*Csum(iregion,ng)/                        &
     &                         Csum( region,ng)
              IF (Master) WRITE (stdout,40) Pregion(iregion),           &
     &                                      Csum(iregion,ng),           &
     &                                      percent
              sumcpu=sumcpu+Csum(iregion,ng)
              sumper=sumper+percent
            END IF
          END DO
 40       FORMAT (2x,a,t53,f14.3,2x,'(',f7.4,' %)')
          IF (Master) WRITE (stdout,50) sumcpu, sumper
 50       FORMAT (t47,'Total:',f14.3,2x,f8.4)
# ifdef DISTRIBUTE
!
! Report elapsed time for message passage communications.
!
          IF (Master) WRITE (stdout,30) ' Message Passage profile:'
          sumcpu=0.0_r8
          sumper=0.0_r8
          DO iregion=38,Nregion
            IF (Csum(iregion,ng).gt.0.0_r8) THEN
              percent=100.0_r8*Csum(iregion,ng)/                        &
     &                         Csum( region,ng)
              IF (Master) WRITE (stdout,40) Pregion(iregion),           &
     &                                      Csum(iregion,ng),           &
     &                                      percent
              sumcpu=sumcpu+Csum(iregion,ng)
              sumper=sumper+percent
            END IF
          END DO
          IF (Master) WRITE (stdout,50) sumcpu, sumper
# endif
#endif
        END IF
!$OMP END CRITICAL (FINALIZE_WCLOCK)
      END IF
      RETURN
      END SUBROUTINE wclock_off
