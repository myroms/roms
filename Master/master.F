#include "cppdefs.h"
/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2008 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Master program to run ROMS/TOMS as single ocean mode or to couple **
**  (sequentially or concurrently) to atmospheric model,  wave model, **
**  or any other modeling system.                                     **
**                                                                    **
************************************************************************
*/

#if defined MODEL_COUPLING
# if defined MCT_LIB
#  if defined AIR_OCEAN
#   include "mct_air_ocean.h"
#  elif defined WAVES_OCEAN
#   include "mct_waves_ocean.h"
#  endif
# elif defined ESMF_LIB
#  include "esmf_coupler.h"
# endif
#else
# if defined ESMF_LIB
#  include "esmf_ocean.h"
# else
#  include "ocean.h"
# endif
#endif
