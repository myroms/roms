#include "cppdefs.h"
      MODULE mod_iounits
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  ESname      Input error subspace file name.                         !
!  Rerror      Running error messages.                                 !
!  SSHname     Input SSH observations file name.                       !
!  SSTname     Input SST observations file name.                       !
!  TOBSname    Input tracer observations file name.                    !
!  VSURname    Input surface currents observations file name.          !
!  VOBSname    Input horizontal currents observations file name.       !
!  aparnam     Input assimilation parameters file name.                !
!  avgbase     Output averages base file name.                         !
!  avgname     Output averages file name.                              !
!  bparnam     Input biology parameters file name.                     !
!  bryname     Input boundary data file name.                          !
!  clmname     Input climatology file name.                            !
!  fltname     Output floats file name.                                !
!  fposnam     Input initial floats positions file name.               !
!  frcids      NetCDF ID descriptors associated with frcname.          ! 
!  frcname     Input forcing fields file name(s).                      !
!  fwdname     Input forward solution file name.                       !
!  grdname     Input grid file name.                                   !
!  hisbase     Output history base file name.                          !
!  hisname     Output history file name.                               !
!  ininame     Input initial conditions file name.                     !
!  nFfiles     Number of forcing files.                                !
!  rstname     Output restart file name.                               !
!  sparnam     Input sediment transport parameters file name.          !
!  sposnam     Input station positions file name.                      !
!  staname     Output station data file name.                          !
!  stdinp      Unit number for standard input (often 5).               !
!  stdout      Unit number for standard output (often 6).              !
!  usrout      Unit number for generic USER output.                    !
!  usrname     USER input/output generic file name.                    !
!  varname     Input IO variables information file name.               !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none

        integer :: stdinp = 5
        integer :: stdout = 6
        integer :: usrout = 10

        integer, dimension(Ngrids) :: nFfiles

        integer, allocatable :: frcids(:,:)

        character (len=50), dimension(6) :: Rerror =                    &
     &       (/ ' ROMS/TOMS - Blows up ................ exit_flag: ',   &
     &          ' ROMS/TOMS - Input error ............. exit_flag: ',   &
     &          ' ROMS/TOMS - Output error ............ exit_flag: ',   &
     &          ' ROMS/TOMS - IO error ................ exit_flag: ',   &
     &          ' ROMS/TOMS - Configuration error ..... exit_flag: ',   &
     &          ' ROMS/TOMS - Partition error ......... exit_flag: ' /)

        character (len=80), allocatable :: frcname(:,:)

        character (len=80), dimension(Ngrids) :: ESname
        character (len=80), dimension(Ngrids) :: SSHname
        character (len=80), dimension(Ngrids) :: SSTname
        character (len=80), dimension(Ngrids) :: TOBSname
        character (len=80), dimension(Ngrids) :: VSURname
        character (len=80), dimension(Ngrids) :: VOBSname
        character (len=80), dimension(Ngrids) :: avgbase
        character (len=80), dimension(Ngrids) :: avgname
        character (len=80), dimension(Ngrids) :: bryname
        character (len=80), dimension(Ngrids) :: clmname
        character (len=80), dimension(Ngrids) :: fltname
        character (len=80), dimension(Ngrids) :: fwdname
        character (len=80), dimension(Ngrids) :: grdname
        character (len=80), dimension(Ngrids) :: hisbase
        character (len=80), dimension(Ngrids) :: hisname
        character (len=80), dimension(Ngrids) :: ininame
        character (len=80), dimension(Ngrids) :: rstname
        character (len=80), dimension(Ngrids) :: staname

        character (len=80) :: aparnam
        character (len=80) :: bparnam
        character (len=80) :: fposnam
        character (len=80) :: sparnam
        character (len=80) :: sposnam
        character (len=80) :: usrname
        character (len=80) :: varname

      END MODULE mod_iounits
