      MODULE mod_kinds
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!
        implicit none

        integer, parameter :: i8 = selected_int_kind(8)
        integer, parameter :: c8 = selected_real_kind(6,30)    ! 32-bit
        integer, parameter :: r4 = selected_real_kind(6,30)    ! 32-bit
        integer, parameter :: r8 = selected_real_kind(12,300)  ! 64-bit
#if defined SUN || defined AIX || defined NEC || defined SGI || \
    defined CRAYX1 || defined DEC
        integer, parameter :: r16 = selected_real_kind(24,270) ! 128-bit
#else
        integer, parameter :: r16 = selected_real_kind(15,300) ! 128-bit
#endif

      END MODULE mod_kinds
