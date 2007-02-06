/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for North Pacific Application.
**
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF4
#define MIX_GEO_TS
#define DIFF_GRID
#define NONLIN_EOS
#define SALINITY
#define SOLVE3D
#define CURVGRID
#define MASKING
#define SPLINES
#define QCORRECTION
#define AVERAGES
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#endif
#define TCLIMATOLOGY
#define TCLM_NUDGING
#define ANA_BSFLUX
#define ANA_BTFLUX
