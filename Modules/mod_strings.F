#include "cppdefs.h"
      MODULE mod_strings
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  title       Title of model run.                                     !
!  Coptions    Activated C-preprocessing options.                      !
!  Pregion     Model regions identifiers used for time profiling.      !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

        character (len=80)   :: title
        character (len=1024) :: Coptions

#ifdef PROFILE
        character (len=50), dimension(35) :: Pregion =                  &
     &    (/'Initialization ...................................',       &
     &      'OI data assimilation .............................',       &
     &      'Reading of input data ............................',       &
     &      'Processing of input data .........................',       &
     &      'Processing of output time averaged data ..........',       &
     &      'Computation of vertical boundary conditions ......',       &
     &      'Computation of diagnostics .......................',       &
     &      'Writing of output data ...........................',       &
     &      'Model 2D kernel ..................................',       &
     &      'Lagrangian floats trajectories ...................',       &
     &      'Tidal forcing ....................................',       &
     &      'Vertical metrics .................................',       &
     &      'Omega vertical velocity ..........................',       &
     &      'Equation of state for seawater ...................',       &
     &      'Biological MODULE, source/sink terms .............',       &
     &      'Sediment tranport MODULE, source/sink terms ......',       &
     &      'Atmospheric boundary layer coupler ...............',       &
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
     &      'Corrector Time-step for tracers ..................'/)
#endif

      END MODULE mod_strings
