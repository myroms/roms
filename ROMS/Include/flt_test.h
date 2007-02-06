/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Floats Tracking Test.
**
*/

#define UV_ADV
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define FLOATS
#define MASKING
#define EW_PERIODIC
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define SOLVE3D
#ifdef SOLVE3D
# define DJ_GRADPS
# define TS_A4HADVECTION
# define TS_A4VADVECTION
# define BODYFORCE
# define SPLINES
# define ANA_BTFLUX
# define ANA_STFLUX
#endif

