#include "cppdefs.h"
      MODULE step2d_mod
#ifdef NONLINEAR
/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2007 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This routine performs a fast (predictor or corrector) time-step    **
** for the free-surface and 2D momentum nonlinear equations.          **
**                                                                    **
************************************************************************
*/
      implicit none

      PRIVATE
      PUBLIC  :: step2d

      CONTAINS

# if defined NEARSHORE_MELLOR || defined WET_DRY
#  include "step2d_LF_AM3s.h"
# else
#  include "step2d_LF_AM3.h"
# endif

#endif
      END MODULE step2d_mod
