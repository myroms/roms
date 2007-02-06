!-----------------------------------------------------------------------
! CVS $Id$
! CVS $Name: MCT_1_0_12 $ 
!-----------------------------------------------------------------------
! A driver model code for Multi-Process Handshaking utility
! to facilitate a plug & play style programming using single executable.
! each processor only execute one component model once. 
! Written by Yun (Helen) He and Chris Ding, NERSC/LBNL, October 2000.


       program main
       use MPH_all
       implicit none
       integer myProc_global

       external ccm3, cpl, pop2_2

       call MPI_INIT(ierr)
       call MPI_COMM_RANK(MPI_COMM_WORLD,myProc_global,ierr)

! here ccm3.8, pop2.2 etc are subroutine names in component models
! you could list the components in any order or omit any of them
       call MPH_setup_SE (atmosphere=ccm3, coupler=cpl, ocean=pop2_2)

!     write(*,*)'I am proc ', MPH_global_proc_id(), 
!    &  ' of global proc ', MPH_local_proc_id_ME_SE(), ' of ',
!    &  MPH_myName_ME_SE()
!     write(*,*)'=============================================='

       call MPI_FINALIZE(ierr)
       end program

