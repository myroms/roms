#include "cppdefs.h"
      MODULE mod_strings
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  cdt         F90/F95 compiler used.                                  !
!  fflags      F90/F95 compiler flags.                                 !
!  title       Title of model run.                                     !
!  Coptions    Activated C-preprocessing options.                      !
!  Mstate      Model state processing messages.                        !
!  Pregion     Model regions identifiers used for time profiling.      !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        character (len=80)   :: title
        character (len=1024) :: Coptions

        integer, parameter :: Nregion = 44

        character (len=40), dimension(4) :: Mstate =                    &
     &    (/'Read state initial conditions,          ',                 &
     &      'Read previous state initial conditions, ',                 &
     &      'Read previous adjoint state solution,   ',                 &
     &      'Read latest adjoint state solution,     '/)

        character (len=50), dimension(Nregion) :: Pregion =             &
     &    (/'Initialization ...................................',       &
     &      'OI data assimilation .............................',       &
     &      'Reading of input data ............................',       &
     &      'Processing of input data .........................',       &
     &      'Processing of output time averaged data ..........',       &
     &      'Computation of vertical boundary conditions ......',       &
     &      'Computation of global information integrals ......',       &
     &      'Writing of output data ...........................',       &
     &      'Model 2D kernel ..................................',       &
     &      'Lagrangian floats trajectories ...................',       &
     &      'Tidal forcing ....................................',       &
     &      '2D/3D coupling, vertical metrics .................',       &
     &      'Omega vertical velocity ..........................',       &
     &      'Equation of state for seawater ...................',       &
     &      'Biological module, source/sink terms .............',       &
     &      'Sediment tranport module, source/sink terms ......',       &
     &      'Atmosphere-Ocean bulk flux parameterization ......',       &
     &      'KPP vertical mixing parameterization .............',       &
     &      'GLS vertical mixing parameterization .............',       &
     &      'My2.5 vertical mixing parameterization ...........',       &
     &      '3D equations right-side terms ....................',       &
     &      '3D equations predictor step ......................',       &
     &      'Pressure gradient ................................',       &
     &      'Harmonic mixing of tracers, S-surfaces ...........',       &
     &      'Harmonic mixing of tracers, geopotentials ........',       &
     &      'Harmonic mixing of tracers, isopycnals ...........',       &
     &      'Biharmonic mixing of tracers, S-surfaces .........',       &
     &      'Biharmonic mixing of tracers, geopotentials ......',       &
     &      'Biharmonic mixing of tracers, isopycnals .........',       &
     &      'Harmonic sress tensor, S-surfaces ................',       &
     &      'Harmonic stress tensor, geopotentials ............',       &
     &      'Biharmonic stress tensor, S-surfaces .............',       &
     &      'Biharmonic stress tensor, geopotentials ..........',       &
     &      'Corrector time-step for 3D momentum ..............',       &
     &      'Corrector time-step for tracers ..................',       &
     &      'Two-way Atmosphere-Ocean models coupling .........',       &
     &      'Bottom boundary layer module .....................',       &
     &      'Message Passage: halo exchange ...................',       &
     &      'Message Passage: data broadcast ..................',       &
     &      'Message Passage: data reduction ..................',       &
     &      'Message Passage: data gathering ..................',       &
     &      'Message Passage: data scattering..................',       &
     &      'Message Passage: boundary data gathering .........',       &
     &      'Message Passage: point data gathering ............'/)

        character (len=80) :: my_os = MY_OS
        character (len=80) :: my_cpu = MY_CPU
        character (len=80) :: my_fort = MY_FORT
        character (len=80) :: my_fc = MY_FC
        character (len=160) :: my_fflags = MY_FFLAGS

      END MODULE mod_strings
