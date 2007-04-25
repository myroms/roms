#include "cppdefs.h"
      MODULE ocean_coupler_mod

#ifdef MODEL_COUPLING
/*
**
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2007 The ROMS/TOMS Group        John C. Warner  **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This module is used to communicate and exchange data between       **
** ROMS/TOMS and other coupled model(s).                              **
**                                                                    **
************************************************************************
*/

# if defined SWAN_COUPLING && defined MCT
#  include "mct_roms_swan.h"
# elif defined REFDIF_COUPLING && defined MCT
!!#  include "mct_roms_refdif.h"
# elif defined WRF_COUPLING && defined MCT
#  include "mct_roms_wrf.h"
# endif

#endif
      END MODULE ocean_coupler_mod
