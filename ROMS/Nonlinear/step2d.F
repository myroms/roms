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
!  This subroutine performs a fast (predictor or corrector) time-step  !
!  for the free-surface  and 2D momentum nonlinear equations.
# ifdef SOLVE3D
!  It also calculates the time filtering variables over all fast-time  !
!  steps  to damp high frequency signals in 3D applications.           !
# endif
!                                                                      !
!=======================================================================
!
# include "step2d_LF_AM3.h"
#else
      MODULE step2d_mod
      END MODULE step2d_mod
#endif
