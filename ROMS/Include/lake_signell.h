/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Forced Lake Sediment Test.
**
*/

#define UV_ADV
#undef  UV_COR
#define DJ_GRADPS
#undef  TS_U3HADVECTION
#undef  TS_C4VADVECTION
#define TS_MPDATA
#undef  WJ_GRADP
#define DJ_GRADPS
#define SALINITY
#define SOLVE3D
#define SPLINES
#define FLOATS
#define AVERAGES
#if defined AVERAGES && defined BEDLOAD
# define AVERAGES_BEDLOAD
#endif

#define NORTHERN_WALL
#define SOUTHERN_WALL
#define EASTERN_WALL
#define WESTERN_WALL
#define UV_PSOURCE
#define TS_PSOURCE
#define ANA_PSOURCE

#undef  UV_LOGDRAG
#undef  MB_BBL
#undef  SG_BBL
#define SSW_BBL

#ifdef SG_BBL
# define SG_CALC_ZNOT
# undef  SG_LOGINT
#endif

#ifdef MB_BBL
# define MB_CALC_ZNOT
# undef  MB_Z0BIO
# undef  MB_Z0BL
# undef  MB_Z0RIP
#endif

#ifdef SSW_BBL
# define SSW_CALC_ZNOT
# undef  SSW_LOGINT
#endif

#if defined MB_BBL || defined SG_BBL || defined SSW_BBL
# define ANA_WWAVE
# ifndef ANA_WWAVE
#  define SWAN
# endif
#endif

#define SEDIMENT
#ifdef SEDIMENT
# define SUSPLOAD
# define BEDLOAD
#endif
#if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
# define ANA_SEDIMENT
#endif

#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BPFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_SPFLUX
#define ANA_SRFLUX

#undef ANA_VMIX
#undef MY25_MIXING
#define GLS_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# define  N2S2_HORAVG
# undef CRAIG_BANNER
# undef CHARNOK
#endif
