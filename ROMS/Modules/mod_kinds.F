      MODULE mod_kinds
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!
        implicit none

        integer, parameter :: i1b= selected_int_kind(2)        !  1-byte
        integer, parameter :: i2b= selected_int_kind(4)        !  2-byte
        integer, parameter :: i4b= selected_int_kind(9)        !  4-byte
        integer, parameter :: c8 = selected_real_kind(6,30)    ! 32-byte
        integer, parameter :: r4 = selected_real_kind(6,30)    ! 32-byte
        integer, parameter :: r8 = selected_real_kind(12,300)  ! 64-byte
#if defined SUN || defined AIX || defined NEC || defined SGI || \
    defined CRAYX1 || defined DEC
        integer, parameter :: r16 = selected_real_kind(24,270) !128-byte
#else
        integer, parameter :: r16 = selected_real_kind(15,300) !128-byte
#endif

      END MODULE mod_kinds
