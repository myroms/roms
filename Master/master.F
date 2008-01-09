#include "cppdefs.h"
/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2008 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Master program to run ROMS/TOMS on single ocean mode or to couple **
**  it (sequentially or concurrently) to an atmospheric model  or any **
**  other Earth modeling system.                                      **
**                                                                    **
************************************************************************
*/

#if defined AIR_OCEAN
# include "air_ocean.h"
#elif defined WAVES_OCEAN
# include "waves_ocean.h"
#else
# include "ocean.h"
#endif
