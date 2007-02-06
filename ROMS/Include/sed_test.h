/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Sediment Tests.
**
*/

#if defined SED_TEST1

/*
**  Suspended Sediment Test in a Channel.
*/

# define UV_ADV
# define UV_PSOURCE
# define UV_LOGDRAG
# define UV_VIS4
# define MIX_S_UV
# define TS_U3HADVECTION
# define TS_C4VADVECTION
# define TS_DIF4
# define MIX_S_TS
# define SALINITY
# define SOLVE3D

# define SEDIMENT
# ifdef SEDIMENT
#  define SUSPLOAD
#  define BEDLOAD
# endif

# undef  SPLINES
# define AVERAGES
# ifdef AVERAGES
#  define AVERAGES_AKV
#  define AVERAGES_AKT
#  define AVERAGES_AKS
# endif

# define NORTHERN_WALL
# define SOUTHERN_WALL
# define WEST_FSRADIATION
# define WEST_M2RADIATION
# define WEST_M3RADIATION
# define WEST_TGRADIENT
# undef  EAST_FSRADIATION
# define EAST_FSCLAMPED
# define EAST_M2RADIATION
# define EAST_M3RADIATION
# undef  EAST_TGRADIENT
# define EAST_TCLAMPED

# undef  ANA_VMIX
# undef  GLS_MIXING
# undef  LMD_MIXING
# define MY25_MIXING

# ifdef MY25_MIXING
#  define KANTHA_CLAYSON
#  undef  N2S2_HORAVG
# endif

# ifdef  LMD_MIXING
#  define LMD_RIMIX
#  define LMD_SKPP
#  define LMD_BKPP
# endif

# define ANA_BPFLUX
# define ANA_BSFLUX
# define ANA_BTFLUX
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SEDIMENT
# define ANA_SMFLUX
# define ANA_SPFLUX
# define ANA_SRFLUX
# define ANA_SSFLUX
# define ANA_STFLUX
# define ANA_PSOURCE
# define ANA_TOBC
# define ANA_FSOBC

#endif
