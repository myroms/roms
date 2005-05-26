#include "cppdefs.h"
      MODULE mod_iounits
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  ADJbase     Output adjoint base file name.                          !
!  ADJname     Output adjoint file name.                               !
!  AVGbase     Output averages base file name.                         !
!  AVGname     Output averages file name.                              !
!  BRYname     Input boundary data file name.                          !
!  CLMname     Input climatology file name.                            !
!  DIAbase     Output diagnostics base file name.                      !
!  DIAname     Output diagnostics file name.                           !
!  FLTname     Output floats file name.                                !
!  FRCids      NetCDF ID descriptors associated with FRCname.          ! 
!  FRCname     Input forcing fields file name(s).                      !
!  FWDname     Input forward solution file name.                       !
!  GRDname     Input grid file name.                                   !
!  HISbase     Output history base file name.                          !
!  HISname     Output history file name.                               !
!  INIname     Input nonlinear initial conditions file name.           !
!  ITLname     Input tangent linear initial conditions file name.      !
!  Iname       Physical parameters standard input script file name.    !
!  RSTname     Output restart file name.                               !
!  OBSname     Input/output datum observations file name.              !
!  Rerror      Running error messages.                                 !
!  SSHname     Input SSH observations file name.                       !
!  SSTname     Input SST observations file name.                       !
!  TLMbase     Output tangent linear base file name.                   !
!  TLMname     Output tangent linear file name.                        !
!  TOBSname    Input tracer observations file name.                    !
!  USRname     USER input/output generic file name.                    !
!  VSURname    Input surface currents observations file name.          !
!  VOBSname    Input horizontal currents observations file name.       !
!  aparnam     Input assimilation parameters file name.                !
!  bparnam     Input biology parameters file name.                     !
!  fposnam     Input initial floats positions file name.               !
!  ioerror     IO error flag.                                          !
!  nFfiles     Number of forcing files.                                !
!  sparnam     Input sediment transport parameters file name.          !
!  sposnam     Input station positions file name.                      !
!  STAname     Output station data file name.                          !
!  stdinp      Unit number for standard input (often 5).               !
!  stdout      Unit number for standard output (often 6).              !
!  usrout      Unit number for generic USER output.                    !
!  varname     Input IO variables information file name.               !
!                                                                      !
!=======================================================================
!
        USE mod_param

        implicit none

        integer, parameter :: stdinp = 5
        integer, parameter :: stdout = 6
        integer, parameter :: usrout = 10
        integer :: ioerror = 0

        integer, dimension(Ngrids) :: nFfiles

        integer, allocatable :: FRCids(:,:)

        character (len=50), dimension(6) :: Rerror =                    &
     &       (/ ' ROMS/TOMS - Blows up ................ exit_flag: ',   &
     &          ' ROMS/TOMS - Input error ............. exit_flag: ',   &
     &          ' ROMS/TOMS - Output error ............ exit_flag: ',   &
     &          ' ROMS/TOMS - IO error ................ exit_flag: ',   &
     &          ' ROMS/TOMS - Configuration error ..... exit_flag: ',   &
     &          ' ROMS/TOMS - Partition error ......... exit_flag: ' /)

        character (len=80), allocatable :: FRCname(:,:)

        character (len=80), dimension(Ngrids) :: AVGbase
        character (len=80), dimension(Ngrids) :: AVGname
        character (len=80), dimension(Ngrids) :: ADJbase
        character (len=80), dimension(Ngrids) :: ADJname
        character (len=80), dimension(Ngrids) :: BRYname
        character (len=80), dimension(Ngrids) :: CLMname
        character (len=80), dimension(Ngrids) :: DIAbase
        character (len=80), dimension(Ngrids) :: DIAname
        character (len=80), dimension(Ngrids) :: FLTname
        character (len=80), dimension(Ngrids) :: FWDname
        character (len=80), dimension(Ngrids) :: GRDname
        character (len=80), dimension(Ngrids) :: HISbase
        character (len=80), dimension(Ngrids) :: HISname
        character (len=80), dimension(Ngrids) :: INIname
        character (len=80), dimension(Ngrids) :: ITLname
        character (len=80), dimension(Ngrids) :: OBSname
        character (len=80), dimension(Ngrids) :: SSHname
        character (len=80), dimension(Ngrids) :: SSTname
        character (len=80), dimension(Ngrids) :: TLMbase
        character (len=80), dimension(Ngrids) :: TLMname
        character (len=80), dimension(Ngrids) :: TOBSname
        character (len=80), dimension(Ngrids) :: VSURname
        character (len=80), dimension(Ngrids) :: VOBSname
        character (len=80), dimension(Ngrids) :: RSTname
        character (len=80), dimension(Ngrids) :: STAname

        character (len=80) :: Iname
        character (len=80) :: USRname
        character (len=80) :: aparnam
        character (len=80) :: bparnam
        character (len=80) :: fposnam
        character (len=80) :: sparnam
        character (len=80) :: sposnam
        character (len=80) :: varname

      END MODULE mod_iounits
