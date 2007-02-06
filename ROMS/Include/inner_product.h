/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for TLM/ADM Inner Product Test.
**
*/

#define SOLVE3D

#ifdef SOLVE3D                         /* 3D Application */
# define UV_ADV
# define UV_COR
# define UV_LDRAG
# define UV_VIS2
# define MIX_S_UV
# define UV_COR
# define TS_U3HADVECTION
# define TS_DIF2
# define MIX_S_TS
# define WJ_GRADP
# define EAST_M2GRADIENT
# define WEST_M2GRADIENT
# define SOUTH_M2GRADIENT
# define NORTH_M2GRADIENT
# define EAST_FSGRADIENT
# define WEST_FSGRADIENT
# define SOUTH_FSGRADIENT
# define NORTH_FSGRADIENT
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
#else                                  /* 2D Application */
# define UV_ADV
# define UV_LDRAG
# define UV_VIS2
# define UV_COR
# define ANA_GRID
# define ANA_INITIAL
# define EAST_M2GRADIENT
# define WEST_M2GRADIENT
# define SOUTH_M2GRADIENT
# define NORTH_M2GRADIENT
# define EAST_FSGRADIENT
# define WEST_FSGRADIENT
# define SOUTH_FSGRADIENT
# define NORTH_FSGRADIENT
# define ANA_SMFLUX
#endif
