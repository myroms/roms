/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2015 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Shore Face Planar Beach Test Case.
**
** Application flag:   SHOREFACE
** Input scripts:      ocean_shoreface.h
**                     sediment_shoreface.h
*/

#define UV_VIS2
#define MIX_S_UV
#define DIAGNOSTICS_UV
#define DIAGNOSTICS_TS
#define AVERAGES
#define WET_DRY
#define NEARSHORE_MELLOR08
#define OUT_DOUBLE
#define UV_ADV
#define TS_MPDATA
#define DJ_GRADPS
#define SALINITY
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define SOLVE3D

#define MASKING
#ifdef MASKING
# define ANA_MASK
#endif
#define ANA_GRID
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#define UV_QDRAG

#ifdef SOLVE3D
# undef  SSW_BBL
# ifdef SSW_BBL
#  define SSW_CALC_ZNOT
#  undef  SSW_LOGINT
# endif

# define SEDIMENT
# ifdef SEDIMENT
#  undef  SED_MORPH
#  define SUSPLOAD
#  define BEDLOAD_MPM
#  undef  BEDLOAD_SOULSBY
# endif
# if defined SEDIMENT || defined SG_BBL || defined MB_BBL || defined SSW_BBL
#  define ANA_SEDIMENT
# endif

# define ANA_STFLUX
# define ANA_SSFLUX
# define ANA_BPFLUX
# define ANA_BTFLUX
# define ANA_BSFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
# undef  ANA_VMIX

# define GLS_MIXING
# if defined GLS_MIXING
#  define KANTHA_CLAYSON
#  define N2S2_HORAVG
#  define RI_SPLINES
#  undef CRAIG_BANNER
#  undef CHARNOK
#  undef ZOS_HSIG
#  undef TKE_WAVEDISS
# endif

#endif
