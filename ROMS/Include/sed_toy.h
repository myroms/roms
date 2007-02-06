/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for One-Dimensional Sediment Toy.
**
*/

#define UV_ADV
#undef  UV_COR
#undef  UV_LOGDRAG
#define UV_QDRAG
#define BODYFORCE
#define DJ_GRADPS
#undef  TS_U3HADVECTION
#define TS_MPDATA
#undef  NONLIN_EOS
#undef  SALINITY

#undef  FLOATS
#undef  SPLINES
#define SOLVE3D
#define AVERAGES
#undef  AVERAGES_AKV
#undef  AVERAGES_AKT
#undef  AVERAGES_AKS
#define AVERAGES_BEDLOAD
#define OUT_DOUBLE

#define EW_PERIODIC
#define NS_PERIODIC

#define ANA_BPFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SEDIMENT
#define ANA_SMFLUX
#define ANA_SPFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_WWAVE

#undef  MB_BBL
#undef  SG_BBL
#undef  SSW_BBL

#ifdef SG_BBL
# undef  SG_CALC_ZNOT
# undef  SG_LOGINT
#endif

#ifdef MB_BBL
# undef  MB_CALC_ZNOT
# undef  MB_Z0BIO
# undef  MB_Z0BL
# undef  MB_Z0RIP
#endif

#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# undef  SSW_LOGINT
#endif

#define GLS_MIXING
#ifdef GLS_MIXING
# define KANTHA_CLAYSON
# define N2S2_HORAVG
# undef  CRAIG_BANNER
# undef  CHARNOK
#endif

#define SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
# undef  BEDLOAD
# undef  SED_DENS
# undef  SED_MORPH
#endif
