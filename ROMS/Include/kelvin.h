/*
** git $Id: 560fb376ff8a4576170ebcd4b459de6bcce908f6 $
** svn $Id: kelvin.h 937 2019-01-28 06:13:04Z arango $
*******************************************************************************
** Copyright (c) 2002-2019 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Kelvin Wave Test.
**
** Application flag:   KELVIN
** Input script:       roms_kelvin.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_DIF2
#define MIX_S_TS
#define SOLVE3D
#define RADIATION_2D
#define ANA_GRID
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SRFLUX
#define ANA_BTFLUX
