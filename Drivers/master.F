/*
************************************************************************
*** Copyright (c) 2002 ROMS/TOMS Group, version 2.1                   **
*************************************************** Hernan G. Arango ***
**                                                                    **
**  Master program to run ROMS/TOMS on single ocean mode or to couple **
**  it (sequentially or concurrently) to an atmospheric model  or any **
**  other Earth modeling system.                                      **
**                                                                    **
************************************************************************
*/

#if defined AIR_OCEAN
# include "air_ocean.h"
#else
# include "ocean.h"
#endif
