/*
** git $Id$
*******************************************************************************
** Copyright (c) 2002-2024 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.md                                                     **
*******************************************************************************
**
** Options for North Atlantic DAMEE Application, 3/4 degree resolution
**
** Application flag:   DAMEE_4
** Input script:       roms_damee_4.in
*/

#define UV_ADV
#define UV_COR
#define UV_QDRAG
#define DJ_GRADPS
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define NONLIN_EOS
#define SALINITY
#define SOLVE3D
#define MASKING
#define QCORRECTION
#define SRELAXATION
#define CURVGRID
#define AVERAGES
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
# define RI_SPLINES
#endif
#define ANA_BSFLUX
#define ANA_BTFLUX
