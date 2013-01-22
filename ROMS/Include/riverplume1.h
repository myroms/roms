/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2013 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for River Plume Test (original version).
**
** Application flag:   RIVERPLUME1
** Input script:       ocean_riverplume1.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_PSOURCE
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define TS_PSOURCE
#define NONLIN_EOS
#define SALINITY
#define MASKING
#define SOLVE3D
#define SPLINES
#define AVERAGES

#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
#endif

#define ANA_GRID
#define ANA_INITIAL
#define ANA_PSOURCE
#define ANA_SMFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
