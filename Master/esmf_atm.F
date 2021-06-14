#include "cppdefs.h"
#if defined MODEL_COUPLING && defined ESMF_LIB
/*
** git $Id$
** svn $Id: esmf_atm.F 1071 2021-06-14 04:22:26Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2021 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Atmospheric model gridded component coupled to ROMS using the     **
**  ESMF library with the NUOPC layer.                                **
**                                                                    **
**  ESMF:   Earth System Modeling Framework (Version 7 or higher)     **
**            https://www.earthsystemcog.org/projects/esmf            **
**                                                                    **
**  NUOPC:  National Unified Operational Prediction Capability        **
**           https://www.earthsystemcog.org/projects/nuopc            **
**                                                                    **
************************************************************************
*/

# if defined ATM_COUPLING
#  if defined COAMPS_COUPLING
#   include "esmf_atm_coamps.h"
#  elif defined REGCM_COUPLING
#   include "esmf_atm_regcm.h"
#  elif defined WRF_COUPLING
#   include "esmf_atm_wrf.h"
#  else
#   include "esmf_atm_void.h"
#  endif
# else
      MODULE esmf_atm_mod
      END MODULE esmf_atm_mod
# endif
#else
      MODULE esmf_atm_mod
      END MODULE esmf_atm_mod
#endif
