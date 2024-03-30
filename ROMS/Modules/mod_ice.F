#include "cppdefs.h"
/*
** git $Id$
** svn $Id: mod_ice.F 1210 2024-01-03 22:03:03Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2024 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
** This module declares ice model variables and parameteres.          **
**                                                                    **
** The current design allows the user to have a lot of latitude for   **
** customizing or adding any particular sea-ice model.                **
**                                                                    **
** Each ice model is composed of seven files:                         **
**                                                                    **
**   (1) Internal model parameters declaration:                       **
**                                                                    **
**         ice_mod.h                                                  **
**                                                                    **
**   (2) Model parameters standard input script:                      **
**                                                                    **
**         ice.in                                                     **
**                                                                    **
**   (3) Code to read input model parameters:                         **
**                                                                    **
**         ice_inp.h                                                  **
**                                                                    **
**   (4) Code to assign indices to model variables during the         **
**       reading of metadata information from "varinfo.yaml":         **
**                                                                    **
**         ice_var.h                                                  **
**                                                                    **
** Note that all the *.h files are located in ROMS/Nonlinear/SeaIce   **
** and included within  <...>  to allow the user to customize any of  **
** them in the project directory  while keeping the distributed code  **
** intact (check the build script for details).                       **
**                                                                    **
************************************************************************
*/
#ifdef SEAICE
# if defined ICE_MODEL
#  include <ice_mod.h>
# endif
#else
      MODULE mod_ice
      END MODULE mod_ice
#endif
