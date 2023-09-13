#include "cppdefs.h"
#if defined MODEL_COUPLING && defined ESMF_LIB
/*
** git $Id$
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2023 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  ROMS gridded component ESMF/NUOPC layer for native or third-party **
**  coupling frameworks like CDEPS/CMEPS.                             **
**                                                                    **
**  CDEPS:  Community Data Models for Earth Prediction Systems        **
**            https://github.com/ESCOMP/CDEPS                         **
**                                                                    **
**  CMEPS:  Community Mediator for Earth Prediction Systems           **
**            https://github.com/ESCOMP/CMEPS                         **
**                                                                    **
**  ESMF:   Earth System Modeling Framework (Version 7 or higher)     **
**            https://www.earthsystemcog.org/projects/esmf            **
**                                                                    **
**  NUOPC:  National Unified Operational Prediction Capability        **
**           https://www.earthsystemcog.org/projects/nuopc            **
**                                                                    **
************************************************************************
*/

# if defined CDEPS || defined CMEPS
#  include "cmeps_roms.h"
# else
#  include "esmf_roms.h"
# endif
#else
     MODULE esmf_roms_mod
     END MODULE esmf_roms_mod
#endif
