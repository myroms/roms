#include "cppdefs.h"
      MODULE rp_step2d_mod
#ifdef TL_IOMS
/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2007 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This routine performs a fast (predictor or corrector)  time-step   **
** for the free-surface and 2D momentum representers tangent linear   **
** equations.                                                         **
**                                                                    **
************************************************************************
*/
      implicit none

      PRIVATE
      PUBLIC  :: rp_step2d

      CONTAINS

# include "rp_step2d_LF_AM3.h"

#endif
      END MODULE rp_step2d_mod
