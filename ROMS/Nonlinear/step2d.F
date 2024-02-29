#include "cppdefs.h"
#ifdef NONLINEAR
!
!git $Id$
!svn $Id: step2d.F 1210 2024-01-03 22:03:03Z arango $
!=======================================================================
!  Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license           Hernan G. Arango   !
!    See License_ROMS.md                    Alexander F. Shchepetkin   !
!==================================================== John C. Warner ===
!                                                                      !
!  This module timesteps the Nonlinear Model (NLM) vertically-         !
!  integrated primitive (2D shallow-water) equations for the           !
!  free-surface and 2D momentum. In 3D applications, the ROMS          !
!  numerical kernel is split between baroclinic and barotropic         !
!  dynamics. The barotropic engine uses a smaller timestep in this     !
!  routine to resolve fast gravity wave processes.                     !
# ifdef SOLVE3D
!                                                                      !
!  It also calculates the time filtering variables over all fast-time  !
!  steps to damp high frequency signals in 3D applications.            !
# endif
!                                                                      !
!=======================================================================
!
# if defined STEP2D_FB_AB3_AM4
#  include "step2d_FB.h"
# elif defined STEP2D_FB_LF_AM3
#  include "step2d_FB_LF_AM3.h"
# else
#  include "step2d_LF_AM3.h"
# endif
#else
      MODULE step2d_mod
      END MODULE step2d_mod
#endif
