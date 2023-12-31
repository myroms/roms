#include "cppdefs.h"
/*
** git $Id$
** svn $Id: master.F 1210 2024-01-03 22:03:03Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2024 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Master program to run ROMS/TOMS as single ocean model or coupled  **
**  to other models using the MCT or ESMF libraries.                  **
**                                                                    **
************************************************************************
*/

#if defined MODEL_COUPLING && !defined CMEPS
# if defined MCT_LIB
#  include "mct_driver.h"
# elif defined ESMF_LIB
#  include "esmf_driver.h"
# else
      PROGRAM master         !  ROMS is coupled with other external
      END PROGRAM master     !  software and drivers not written by us
# endif
#else
# include "roms.h"
#endif
