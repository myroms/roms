#include "cppdefs.h"
/*
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2023 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
** This module computes the ice thermodynamic growth and decay        **
** term. Currently the only option is the Mellor-Kantha (1989)        **
** formulation.                                                       **
**                                                                    **
************************************************************************
*/
#if defined ICE_MODEL && defined ICE_THERMO
# ifdef ICE_MK
#  include <ice_mk.h>
# endif
#else
      MODULE ice_thermo_mod
      END MODULE ice_thermo_mod
#endif

