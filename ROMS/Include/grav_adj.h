/*
** git $Id$
*******************************************************************************
** Copyright (c) 2002-2025 The ROMS Group                                    **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.md                                                     **
*******************************************************************************
**
** Options for Gravitational Adjustment Test.
**
** Application flag:   GRAV_ADJ
** Input script:       roms_grav_adj.in
*/

#define UV_ADV
#define UV_VIS2
#define UV_LDRAG
#define MIX_S_UV
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_DIF2
#define MIX_S_TS
#define SOLVE3D
#define AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define OUT_DOUBLE
#define DIAGNOSTICS_TS
#define DIAGNOSTICS_UV
