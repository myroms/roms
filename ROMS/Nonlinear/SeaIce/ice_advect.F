#include "cppdefs.h"
/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2023 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
** This module computes the advection of the ice tracer fields.       **
**                                                                    **                                                                  
************************************************************************
*/
#if defined ICE_MODEL && defined ICE_ADVECT
# if defined ICE_SMOLAR
#  include <ice_smolar.h>
# endif
#else
      MODULE ice_advect_mod
      END MODULE ice_advect_mod
#endif
