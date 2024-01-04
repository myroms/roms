/*
** git $Id$
** svn $Id: nena.h 1210 2024-01-03 22:03:03Z arango $
*******************************************************************************
** Copyright (c) 2002-2024 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.md                                                     **
*******************************************************************************
**
** Options for North East North America Application.
**
** Application flag:   NENA
** Input script:       ocean_nena.in
*/

#define UV_ADV
#define UV_SADVECTION
#define DJ_GRADPS
#define UV_COR
#define UV_QDRAG
#define UV_VIS2
#define MIX_S_UV
#define SPLINES_VDIFF
#define SPLINES_VVISC
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define MASKING
#define AVERAGES
#define SRELAXATION
#define QCORRECTION
#define SOLAR_SOURCE
#define DIURNAL_SRFLUX
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
# define RI_MIXING
#endif
#define BIO_FENNEL
#ifdef BIO_FENNEL
# define CARBON
# define DENITRIFICATION
# define BIO_SEDIMENT
# define DIAGNOSTICS_BIO
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif

#define ANA_BSFLUX
#define ANA_BTFLUX
