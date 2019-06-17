/*
** git $Id: 560fb376ff8a4576170ebcd4b459de6bcce908f6 $
** svn $Id: riverplume2.h 937 2019-01-28 06:13:04Z arango $
*******************************************************************************
** Copyright (c) 2002-2019 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for River Plume example by Hyatt and Signell described at
** http://smig.usgs.gov/SMIG/features_0300/plumes_inline.html
**
** Application flag:   RIVERPLUME2
** Input script:       roms_riverplume2.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define NONLIN_EOS
#define SALINITY
#define MASKING
#define SOLVE3D
#define AVERAGES

#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
# define RI_SPLINES
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

