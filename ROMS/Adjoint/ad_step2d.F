#include "cppdefs.h"
#ifdef ADJOINT
!
!svn $Id: ad_step2d.F 275 2023-02-09 03:21:09Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2023 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module timesteps the Adjoint Model (ADM) vertically-integrated !
!  primitive (2D shallow-water) equations for the free-surface and 2D  !
!  momentum. In 3D applications, the ROMS numerical kernel is split    !
!  between baroclinic and barotropic dynamics. The barotropic engine   !
!  uses a smaller timestep in this routine to resolve fast gravity     !
!  wave processes.                                                     !
# ifdef SOLVE3D
!                                                                      !
!  It also calculates the time filtering variables over all fast-time  !
!  steps to damp high frequency signals in 3D applications.            !
# endif
!                                                                      !
!=======================================================================
!
# if defined STEP2D_FB_AB3_AM4
#  include "ad_step2d_FB.h"
# elif defined STEP2D_FB_LF_AM3
#  include "ad_step2d_FB_LF_AM3.h"
# else
#  include "ad_step2d_LF_AM3.h"
# endif
#else
      MODULE ad_step2d_mod
      END MODULE ad_step2d_mod
#endif
