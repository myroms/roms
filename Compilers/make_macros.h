/*
** git $Id$
** svn $Id: make_macros.h 1202 2023-10-24 15:36:07Z arango $
********************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2023 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.md                                                     **
*******************************************************************************
**                                                                           **
** This configuration file defines several macros used by the makefile to    **
** select the appropriate code to compile and link. These macros are used    **
** to exclude files from the "modules" and "includes" definitions.           **
**                                                                           **
*******************************************************************************
*/

#include "../ROMS/Include/cppdefs.h"

/*
** Process ROMS model.
*/

  USE_ROMS := on

/*
** Process adjoint model.
*/

#ifdef ADJOINT
  USE_ADJOINT := on
#else
  USE_ADJOINT :=
#endif

/*
** Process tangent linear model.
*/

#ifdef TANGENT
  USE_TANGENT := on
#else
  USE_TANGENT :=
#endif

/*
** Process representers tangent linear model.
*/

#ifdef TL_IOMS
  USE_REPRESENTER := on
#else
  USE_REPRESENTER :=
#endif

/*
** Process ROMS Sea Ice model.
*/

#ifdef ICE_MODEL
  USE_SEAICE := on
#else
  USE_SEAICE :=
#endif

/*
** Process CICE seaice model for coupling.
*/

#ifdef CICE_COUPLING
  USE_CICE := on
#else
  USE_CICE :=
#endif

/*
** Process COAMPS Atmospheric model for coupling.
*/

#ifdef COAMPS_COUPLING
  USE_COAMPS := on
#else
  USE_COAMPS :=
#endif

/*
** Process RegCM Atmospheric model for coupling.
*/

#ifdef REGCM_COUPLING
  USE_REGCM := on
#else
  USE_REGCM :=
#endif

/*
** Process WRF Atmospheric model coupling.
*/

#ifdef WRF_COUPLING
  USE_WRF := on
#else
  USE_WRF :=
#endif

/*
** Process SWAN wave model for coupling.
*/

#ifdef SWAN_COUPLING
  USE_SWAN := on
#else
  USE_SWAN :=
#endif

/*
** Process REFDIF wave model for coupling.
*/

#ifdef REFDIF_COUPLING
  USE_REFDIF := on
#else
  USE_REFDIF :=
#endif

/*
** Process WAM wave model for coupling.
*/

#ifdef WAM_COUPLING
  USE_WAM := on
#else
  USE_WAM :=
#endif

/*
** Determine if the ARPACK library is needed.
*/

#if defined PROPAGATOR
  USE_ARPACK := on
#else
  USE_ARPACK :=
#endif

/*
** Determine if the Model Coupling Tool library is needed.
*/

#ifdef MCT_LIB
  USE_MCT := on
#else
  USE_MCT :=
#endif

/*
** Determine if the Earth System Modeling Framework library is needed.
*/

#ifdef ESMF_LIB
  USE_ESMF := on
#else
  USE_ESMF :=
#endif
