/*
** svn $Id$
*******************************************************************************
** Copyright (c) 2002-2007 The ROMS/TOMS Group                               **
**   Licensed under a MIT/X style license                                    **
**   See License_ROMS.txt                                                    **
*******************************************************************************
**
** Options for Intra-America Sea Application.
**
*/

#undef  AFT_EIGENMODES          /* Adjoint Finite Time Eigenmodes */
#undef  CORRELATION             /* Background-error Correlation Check */
#undef  GRADIENT_CHECK          /* TLM/ADM Gradient Check */
#undef  FORCING_SV              /* Forcing Singular Vectors */
#undef  FT_EIGENMODES           /* Finite Time Eigenmodes */
#undef  IS4DVAR                 /* Incremental, strong constraint 4DVAR */
#define NLM_DRIVER              /* Nonlinear Basic State trajectory */
#undef  OPT_PERTURBATION        /* Optimal perturbations */
#undef  PICARD_TEST             /* Picard Iterations Test */
#undef  R_SYMMETRY              /* Representer Matrix Symmetry Test */
#undef  SANITY_CHECK            /* Sanity Check */
#undef  SO_SEMI                 /* Stochastic Optimals: Semi-norm */
#undef  TLM_CHECK               /* Tangent Linear Model Check */
#undef  W4DPSAS                 /* Weak constraint 4D-PSAS */
#undef  W4DVAR                  /* Weak constraint 4DVAR */

/*
**-----------------------------------------------------------------------------
**  Nonlinear basic state tracjectory.
**-----------------------------------------------------------------------------
*/

#if defined NLM_DRIVER
# define UV_ADV
# define DJ_GRADPS
# define UV_COR
# define UV_QDRAG
# define UV_VIS2
# define MIX_S_UV
# define TS_U3HADVECTION
# define TS_SVADVECTION
# define SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define CURVGRID
# define SPLINES
# define MASKING
# define AVERAGES
# define SRELAXATION
# define QCORRECTION
# define SOLAR_SOURCE
# define DIURNAL_SRFLUX

# define LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  define LMD_NONLOCAL
# endif

# undef  M2CLIMATOLOGY
# undef  M3CLIMATOLOGY
# undef  TCLIMATOLOGY
# undef  ZCLIMATOLOGY
# undef  EAST_VOLCONS
# undef  WEST_VOLCONS
# undef  SOUTH_VOLCONS
# define EAST_FSCHAPMAN
# define EAST_M2FLATHER
# define EAST_M3CLAMPED
# define EAST_TCLAMPED
# define WESTERN_WALL
# define SOUTHERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# define NORTH_M3CLAMPED
# define NORTH_TCLAMPED
# define ANA_BSFLUX
# define ANA_BTFLUX
#endif

/*
**-----------------------------------------------------------------------------
**  Variational Data Assimilation.
**-----------------------------------------------------------------------------
*/

#if defined CORRELATION || defined IS4DVAR      || \
    defined IS4DVAR_OLD || defined SANITY_CHECK || \
    defined W4DPSAS     || defined W4DVAR
# define UV_ADV
# define DJ_GRADPS
# define UV_COR
# define UV_LDRAG
# define UV_VIS2
# define MIX_S_UV
# define TS_U3HADVECTION
# define TS_SVADVECTION
# define SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define CURVGRID
# define SPLINES
# define MASKING
# undef  AVERAGES
# define SRELAXATION
# define QCORRECTION
# define SOLAR_SOURCE
# define DIURNAL_SRFLUX
# undef  LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  define LMD_NONLOCAL
# endif
# undef  GLS_MIXING

# undef  M2CLIMATOLOGY
# undef  M3CLIMATOLOGY
# undef  TCLIMATOLOGY
# undef  ZCLIMATOLOGY
# undef  EAST_VOLCONS
# undef  WEST_VOLCONS
# undef  SOUTH_VOLCONS
# define EAST_FSCHAPMAN
# define EAST_M2FLATHER
# define EAST_M3CLAMPED
# define EAST_TCLAMPED
# define WESTERN_WALL
# define SOUTHERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# define NORTH_M3CLAMPED
# define NORTH_TCLAMPED
# define ANA_BSFLUX
# define ANA_BTFLUX

# ifdef  SANITY_CHECK
#  define ANA_PERTURB
# endif

# define VCONVOLUTION
# define IMPLICIT_VCONV

# define FORWARD_READ
# define FORWARD_MIXING
# define OUT_DOUBLE
#endif

/*
**-----------------------------------------------------------------------------
**  Generalized Stability Theory analysis.
**-----------------------------------------------------------------------------
*/

#if defined AFT_EIGENMODES || defined FT_EIGENMODES    || \
    defined FORCING_SV     || defined OPT_PERTURBATION || \
    defined SO_SEMI
# define UV_ADV
# define DJ_GRADPS
# define UV_COR
# define UV_QDRAG
# define UV_VIS2
# define MIX_S_UV
# undef TS_DIF2
# undef MIX_S_TS
# define TS_U3HADVECTION
# define TS_SVADVECTION
# define SOLVE3D
# define SALINITY
# define NONLIN_EOS
# define CURVGRID
# define SPLINES
# define MASKING
# define SRELAXATION
# define QCORRECTION
# define SOLAR_SOURCE
# define DIURNAL_SRFLUX

# undef LMD_MIXING
# ifdef LMD_MIXING
#  define LMD_RIMIX
#  define LMD_CONVEC
#  define LMD_SKPP
#  define LMD_NONLOCAL
# endif
# undef  GLS_MIXING

# undef  M2CLIMATOLOGY
# undef  M3CLIMATOLOGY
# undef  TCLIMATOLOGY
# undef  ZCLIMATOLOGY
# undef  EAST_VOLCONS
# undef  WEST_VOLCONS
# undef  SOUTH_VOLCONS
# define EAST_FSCHAPMAN
# define EAST_M2FLATHER
# define EAST_M3CLAMPED
# define EAST_TCLAMPED
# define WESTERN_WALL
# define SOUTHERN_WALL
# define NORTH_FSCHAPMAN
# define NORTH_M2FLATHER
# define NORTH_M3CLAMPED
# define NORTH_TCLAMPED
# define ANA_BSFLUX
# define ANA_BTFLUX

# ifdef SO_SEMI
#  undef  SO_SEMI_WHITE
#  define  FULL_GRID
# endif
# define FORWARD_READ
# define FORWARD_MIXING
# undef  OUT_DOUBLE
#endif
