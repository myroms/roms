/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for an Idealized Canyon.
**
*/

#if defined CANYON_A
# define UV_ADV
# define UV_QDRAG
# define UV_VIS2
# define UV_COR
# define EW_PERIODIC
# define BODYFORCE
# define ANA_DIAG
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
#elif defined CANYON_B
# define UV_ADV
# define UV_COR
# define UV_QDRAG
# define UV_VIS2
# define MIX_S_UV
# define DJ_GRADPS
# define TS_A4HADVECTION
# define TS_A4VADVECTION
# define TS_DIF2
# define MIX_GEO_TS
# define SOLVE3D
# define SPLINES
# define EW_PERIODIC
# define ANA_DIAG
# define ANA_GRID
# define ANA_INITIAL
# define ANA_SMFLUX
# define ANA_STFLUX
# define ANA_BTFLUX
# define ANA_VMIX
#endif
