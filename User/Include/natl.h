/*
** git $Id$
*******************************************************************************
** Copyright (c) 2002-2025 The ROMS Group                                    **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.md                                                     **
*******************************************************************************
**
** Options for High Resolution North Atlantic Application.
**
** Application flag:   NATL
** Input script:       ocean_natl.in
*/

#define UV_ADV
#define UV_QDRAG
#define UV_COR
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define STATIONS
#define MASKING
#define WRITE_WATER
#define NO_WRITE_GRID
#define AVERAGES
#define SRELAXATION
#define QCORRECTION
#define SOLAR_SOURCE
#define ANA_BSFLUX
#define ANA_BTFLUX
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif
