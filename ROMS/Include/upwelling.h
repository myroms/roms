/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Upwelling Test.
**
*/

#define UV_ADV
#define UV_COR
#define UV_LDRAG
#define UV_VIS2
#undef  MIX_GEO_UV
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#undef  TS_MPDATA
#undef  WJ_GRADP
#define DJ_GRADPS
#define TS_DIF2
#undef  TS_DIF4
#undef  MIX_GEO_TS
#define MIX_S_TS

#define SALINITY
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define DIAGNOSTICS_TS
#define DIAGNOSTICS_UV
#define EW_PERIODIC

#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX

#define ANA_VMIX
#undef  GLS_MIXING
#undef  MY25_MIXING
#if defined GLS_MIXING || defined MY25_MIXING
# define KANTHA_CLAYSON
# undef  CANUTO_A
# define N2S2_HORAVG
#endif

#undef  BIO_FASHAM
#undef  NPZD_POWELL

#if defined BIO_FASHAM || defined NPZD_POWELL
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_SRFLUX
#endif

#ifdef BIO_FASHAM
# define CARBON
# define DENITRIFICATION
# define BIO_SEDIMENT
# define DIAGNOSTICS_BIO
#endif

#undef  PERFECT_RESTART
