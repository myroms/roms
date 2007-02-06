/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Coupled Boundary Layer Air-Sea Transfer Application.
**
*/

#define	UV_ADV
#define	UV_COR
#define UV_QDRAG
#undef	UV_VIS2
#undef	MIX_S_UV
#define UV_SADVECTION
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define SOLVE3D
#define	SALINITY
#define	NONLIN_EOS
#define CURVGRID
#define	AVERAGES
#define STATIONS
#undef  FLOATS
#define MASKING
#define SPLINES
#define SOLAR_SOURCE
#undef  LMD_MIXING
#ifdef  LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define ANA_CLOUD
#endif
#define MY25_MIXING
#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
#endif
#undef	SG_BBL
#ifdef SG_BBL
# define SG_CALC_ZNOT
# define ANA_SEDIMENT
# define ANA_WWAVE
#endif
#define  BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE
# define ANA_CLOUD
# define ANA_RAIN
#endif
#define RADIATION_2D
#undef  SPONGE
#undef  EAST_VOLCONS
#undef  WEST_VOLCONS
#undef  SOUTH_VOLCONS
#undef  NORTH_VOLCONS
#define RAMP_TIDES
#define	SSH_TIDES
#ifdef SSH_TIDES
# define ADD_FSOBC
# define EAST_FSCHAPMAN
# define WEST_FSCHAPMAN
# define SOUTH_FSCHAPMAN
# define NORTH_FSCHAPMAN
#else
# define EAST_FSGRADIENT
# define WEST_FSGRADIENT
# define SOUTH_FSGRADIENT
# define NORTH_FSGRADIENT
#endif
#define	UV_TIDES
#ifdef UV_TIDES
# define ADD_M2OBC
# define EAST_M2FLATHER
# define WEST_M2FLATHER
# define SOUTH_M2FLATHER
# define NORTH_M2FLATHER
#else
# define EAST_M2RADIATION
# define WEST_M2RADIATION
# define SOUTH_M2RADIATION
# define NORTH_M2RADIATION
#endif
#define EAST_M3RADIATION
#define EAST_M3NUDGING
#define EAST_TRADIATION
#define EAST_TNUDGING
#define WEST_M3RADIATION
#define WEST_M3NUDGING
#define WEST_TRADIATION
#define WEST_TNUDGING
#define SOUTH_M3RADIATION
#define SOUTH_M3NUDGING
#define SOUTH_TRADIATION
#define SOUTH_TNUDGING
#define NORTH_M3RADIATION
#define NORTH_M3NUDGING
#define NORTH_TRADIATION
#define NORTH_TNUDGING
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_SSFLUX
#undef  ANA_STFLUX
#undef  ANA_SMFLUX
#undef  ANA_SRFLUX
