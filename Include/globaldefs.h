/*
** Include file "globaldef.h"
***************************************** Alexander F. Shchepetkin ***
** Copyright (c) 2002 ROMS/TOMS Group                               **
************************************************* Hernan G. Arango ***
**                                                                  **
** WARNING: This  file  contains  a set of  predetermined           **
**          macro definitions which are inserted into the           **
** individual files by the C-preprocessor. It is strongly           **
** recommended to NOT modify any of the definitions below.          **
**                                                                  **
**********************************************************************
*/

/*
** Set assumed-shape array switch.  Imported arrays with dummy 
** arguments that takes the shape of the actual argument passed
** to it.  If off, all the arrays are explicit-shape.  In some
** computer explicit-shape arrays slow down performacnce because
** the arrays are copied when passed by arguments.
*/

#define ASSUMED_SHAPE

/*
** Set internal distributed-memory switch.
*/

#if defined MPI
# define DISTRIBUTE
#endif

/*
** Turn ON/OFF debugging switch. It avoids writing current date and
** CPP options to NetCDF file headers. This is used to compare serial
** and parallel solutions where the UNIX command "diff" is used between
** NetCDF files.  It will only tell us that the binary files are
** different or not.  Finding the parallel bug is complete different
** story.
*/

#define DEBUGGING

/*
** Turn ON/OFF time profiling.
*/

#define PROFILE

/*
** Set default time-averaging filter for barotropic fields.
**
*/

#ifdef SOLVE3D
# undef COSINE2
# define POWER_LAW
#endif

/*
** Turn ON/OFF switch to include/disregard the difference between
** rho0 and surface density in the computation of baroclinic pressure
** term.
*/

#define RHO_SURF

/*
** Activate bacroclinic pressure gradient response due to the
** perturbation of free-surface in the presence of stratification
** and bathymetry.
*/

#ifdef SOLVE3D
# define VAR_RHO_2D
#endif

/*
** Turn ON/OFF double precision for real type variables and
** associated intrinsic functions.
*/

#define DOUBLE_PRECISION

/*
** Define macro for the first 2D time-step.
*/

#ifdef SOLVE3D
# define FIRST_2D_STEP iif(ng).eq.1
#else
# define FIRST_2D_STEP iic(ng).eq.ntstart
#endif

/*
** Define global grid lower and upper bounds in the I- and
** J-directions. These values are a function of periodicity.
** They are used in both shared- and distributed-memory
** configurations.
*/

#ifdef EW_PERIODIC
# ifdef NS_PERIODIC
#  define LOWER_BOUND_I -2
#  define UPPER_BOUND_I Im(ng)+2
#  define LOWER_BOUND_J -2
#  define UPPER_BOUND_J Jm(ng)+2
# else
#  define LOWER_BOUND_I -2
#  define UPPER_BOUND_I Im(ng)+2
#  define LOWER_BOUND_J 0
#  define UPPER_BOUND_J Jm(ng)+1
# endif
#else
# ifdef NS_PERIODIC
#  define LOWER_BOUND_I 0
#  define UPPER_BOUND_I Im(ng)+1
#  define LOWER_BOUND_J -2
#  define UPPER_BOUND_J Jm(ng)+2
# else
#  define LOWER_BOUND_I 0
#  define UPPER_BOUND_I Im(ng)+1
#  define LOWER_BOUND_J 0
#  define UPPER_BOUND_J Jm(ng)+1
# endif
#endif
#define XI_DIM LOWER_BOUND_I:UPPER_BOUND_I
#define ETA_DIM LOWER_BOUND_J:UPPER_BOUND_J
#define GLOBAL_2D_ARRAY XI_DIM,ETA_DIM
#define PRIVATE_1D_SCRATCH_ARRAY Istr-3:Iend+3
#define PRIVATE_2D_SCRATCH_ARRAY Istr-3:Iend+3,Jstr-3:Jend+3

/*
** Set number of ghost-points in the halo region.
*/

#define GHOST_POINTS 2

/*
** Remove OpenMP directives in serial and distributed memory 
** Applications.  This definition will be used in conjunction with
** the pearl script "cpp_clean" to remove the full directive.
*/

#if !defined _OPENMP
# define OMP !
#endif

/*
** Set tile variable for distributed- or shared-memory configurations.
*/

#ifdef DISTRIBUTE
# define TILE MyRank
#else
# define TILE tile
#endif

/*
** The following definitions contain fortran logical expressions
** equivalent to the question: ''Am I the thread working on a tile
** which is adjacent to the WESTERN, EASTERN, SOUTHERN, or NORTHERN
** edges of the model domain?'' These logical expressions are used to
** update domain boundaries and corners.
*/

#define WESTERN_EDGE (Istr.eq.1)
#define EASTERN_EDGE (Iend.eq.Lm(ng))
#define SOUTHERN_EDGE (Jstr.eq.1)
#define NORTHERN_EDGE (Jend.eq.Mm(ng))
#define SOUTH_WEST_CORNER (Istr.eq.1).and.(Jstr.eq.1)
#define NORTH_WEST_CORNER (Istr.eq.1).and.(Jend.eq.Mm(ng))
#define SOUTH_EAST_CORNER (Iend.eq.Lm(ng)).and.(Jstr.eq.1)
#define NORTH_EAST_CORNER (Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))

/*
** The following definitions are fortran logical expressions use
** to update global variables while avoiding mutual overlap between
** threads in shared-memory configurations.
*/

#ifdef DISTRIBUTE
# define SOUTH_WEST_TEST .true.
# define NORTH_WEST_TEST .true.
# define SOUTH_EAST_TEST .true.
# define NORTH_EAST_TEST .true.
#else
# define SOUTH_WEST_TEST (Istr.eq.1).and.(Jstr.eq.1)
# define NORTH_WEST_TEST (Istr.eq.1).and.(Jend.eq.Mm(ng))
# define SOUTH_EAST_TEST (Iend.eq.Lm(ng)).and.(Jstr.eq.1)
# define NORTH_EAST_TEST (Iend.eq.Lm(ng)).and.(Jend.eq.Mm(ng))
#endif

/*
** Choice of double/single precision for real type variables and
** associated intrinsic functions.
*/

#if (defined CRAY || defined CRAYT3E) && !defined CRAYX1
# ifdef  DOUBLE_PRECISION
#  undef  DOUBLE_PRECISION
# endif
#endif

#ifdef DOUBLE_PRECISION
# define nf_get_att_TYPE nf_get_att_double
# define nf_put_att_TYPE nf_put_att_double
# define nf_get_var1_TYPE nf_get_var1_double
# define nf_put_var1_TYPE nf_put_var1_double
# define nf_get_vara_TYPE nf_get_vara_double
# define nf_put_vara_TYPE nf_put_vara_double
#else
# define nf_get_att_TYPE nf_get_att_real
# define nf_put_att_TYPE nf_put_att_real
# define nf_get_var1_TYPE nf_get_var1_real
# define nf_put_var1_TYPE nf_put_var1_real
# define nf_get_vara_TYPE nf_get_vara_real
# define nf_put_vara_TYPE nf_put_vara_real
#endif

/*
** Set output index for multi-time levels variables.
*/

#ifdef SOLVE3D
# define KOUT kstp(ng)
#else
# define KOUT knew(ng)
#endif
#define NOUT nrhs(ng)

/*
** If splines, deactivate horizontal and vertical smoothing of
** Richardson number horizontally and/or vertically.
*/

#ifdef SPLINES
# if defined LMD_MIXING
#  undef RI_HORAVG
#  undef RI_VERAVG
# endif
#endif

/*
** Activate internal switch for the computation of the Brunt-Vaisala
** frequency.
*/

#if defined BVF_MIXING || defined LMD_MIXING  || defined LMD_SKPP    || \
    defined LMD_BKPP   || defined GLS_MIXING  || defined MY25_MIXING
# define BV_FREQUENCY
#endif

/*
** Activate switch for processing climatology data.
*/

#if defined ZCLIMATOLOGY || defined M2CLIMATOLOGY || \
    defined TCLIMATOLOGY || defined M3CLIMATOLOGY || \
    defined ZCLM_NUDGING || defined M2CLM_NUDGING || \
    defined TCLM_NUDGING || defined M3CLM_NUDGING
# define CLIMATOLOGY
#endif

/*
** Activate internal switch for bottom boundary layer closure.
*/

#if defined CS_BBL || defined MB_BBL || defined SG_BBL
# define BBL_MODEL
#endif

/*
** Activate internal switch to set-up nudging coefficients.
*/

#if defined ZCLM_NUDGING    || defined M2CLM_NUDGING   || \
    defined TCLM_NUDGING    || defined M3CLM_NUDGING   || \
    defined WEST_FSNUDGING  || defined EAST_FSNUDGING  || \
    defined SOUTH_FSNUDGING || defined NORTH_FSNUDGING || \
    defined WEST_M2NUDGING  || defined EAST_M2NUDGING  || \
    defined SOUTH_M2NUDGING || defined NORTH_M2NUDGING || \
    defined WEST_TNUDGING   || defined EAST_TNUDGING   || \
    defined SOUTH_TNUDGING  || defined NORTH_TNUDGING  || \
    defined WEST_M3NUDGING  || defined EAST_M3NUDGING  || \
    defined SOUTH_M3NUDGING || defined NORTH_M3NUDGING
# define NUDGING_COFF
#endif

/*
** Activate internal switches requiring open boundary data.
*/

#if (defined WEST_M2RADIATION  && defined WEST_M2NUDGING)  || \
     defined WEST_M2FLATHER    || defined WEST_M2CLAMPED
# define WEST_M2OBC
#endif
#if (defined EAST_M2RADIATION  && defined EAST_M2NUDGING)  || \
     defined EAST_M2FLATHER    || defined EAST_M2CLAMPED
# define EAST_M2OBC
#endif
#if (defined SOUTH_M2RADIATION && defined SOUTH_M2NUDGING) || \
     defined SOUTH_M2FLATHER   || defined SOUTH_M2CLAMPED
# define SOUTH_M2OBC
#endif
#if (defined NORTH_M2RADIATION && defined NORTH_M2NUDGING) || \
     defined NORTH_M2FLATHER   || defined NORTH_M2CLAMPED
# define NORTH_M2OBC
#endif

#if (defined WEST_FSRADIATION  && defined WEST_FSNUDGING)  || \
     defined WEST_M2FLATHER    || defined WEST_FSCLAMPED
# define WEST_FSOBC
#endif
#if (defined EAST_FSRADIATION  && defined EAST_FSNUDGING)  || \
     defined EAST_M2FLATHER    || defined EAST_FSCLAMPED
# define EAST_FSOBC
#endif
#if (defined SOUTH_FSRADIATION && defined SOUTH_FSNUDGING) || \
     defined SOUTH_M2FLATHER   || defined SOUTH_FSCLAMPED
# define SOUTH_FSOBC
#endif
#if (defined NORTH_FSRADIATION && defined NORTH_FSNUDGING) || \
     defined NORTH_M2FLATHER   || defined NORTH_FSCLAMPED
# define NORTH_FSOBC
#endif

#if (defined WEST_M3RADIATION  && defined WEST_M3NUDGING)  || \
     defined WEST_M3CLAMPED
# define WEST_M3OBC
#endif
#if (defined EAST_M3RADIATION  && defined EAST_M3NUDGING)  || \
     defined EAST_M3CLAMPED
# define EAST_M3OBC
#endif
#if (defined SOUTH_M3RADIATION && defined SOUTH_M3NUDGING) || \
     defined SOUTH_M3CLAMPED
# define SOUTH_M3OBC
#endif
#if (defined NORTH_M3RADIATION && defined NORTH_M3NUDGING) || \
     defined NORTH_M3CLAMPED
# define NORTH_M3OBC
#endif

#if (defined WEST_TRADIATION   && defined WEST_TNUDGING)   || \
     defined WEST_TCLAMPED
# define WEST_TOBC
#endif
#if (defined EAST_TRADIATION   && defined EAST_TNUDGING)   || \
     defined EAST_TCLAMPED
# define EAST_TOBC
#endif
#if (defined SOUTH_TRADIATION  && defined SOUTH_TNUDGING)  || \
     defined SOUTH_TCLAMPED
# define SOUTH_TOBC
#endif
#if (defined NORTH_TRADIATION  && defined NORTH_TNUDGING)  || \
     defined NORTH_TCLAMPED
# define NORTH_TOBC
#endif

#ifdef SOLVE3D
# if defined WEST_FSOBC  || defined EAST_FSOBC  || \
     defined SOUTH_FSOBC || defined NORTH_FSOBC || \
     defined WEST_M2OBC  || defined EAST_M2OBC  || \
     defined SOUTH_M2OBC || defined NORTH_M2OBC || \
     defined WEST_M3OBC  || defined EAST_M3OBC  || \
     defined SOUTH_M3OBC || defined NORTH_M3OBC || \
     defined WEST_TOBC   || defined EAST_TOBC   || \
     defined SOUTH_TOBC  || defined NORTH_TOBC
#  define OBC
# endif
#else
# if defined WEST_FSOBC  || defined EAST_FSOBC  || \
     defined SOUTH_FSOBC || defined NORTH_FSOBC || \
     defined WEST_M2OBC  || defined EAST_M2OBC  || \
     defined SOUTH_M2OBC || defined NORTH_M2OBC
#  define OBC
# endif
#endif

/*
** Define internal flag indicating processing of input boundary
** NetCDF file.
*/

#if (!defined ANA_FSOBC && \
     (defined WEST_FSOBC  || defined EAST_FSOBC    || \
      defined SOUTH_FSOBC || defined NORTH_FSOBC)) || \
    (!defined ANA_M2OBC && \
     (defined WEST_M2OBC  || defined EAST_M2OBC    || \
      defined SOUTH_M2OBC || defined NORTH_M2OBC)) || \
    (!defined ANA_M3OBC && \
     (defined WEST_M3OBC  || defined EAST_M3OBC    || \
      defined SOUTH_M3OBC || defined NORTH_M3OBC)) || \
    (!defined ANA_TOBC && \
     (defined WEST_TOBC   || defined EAST_TOBC    || \
      defined SOUTH_TOBC  || defined NORTH_TOBC))
# define OBC_DATA
#endif

/*
** Activate internal switches for volume conservation at open boundary.
*/

#if !defined WEST_M2OBC && defined WEST_VOLCONS
# undef WEST_VOLCONS
#endif
#if !defined EAST_M2OBC && defined EAST_VOLCONS
# undef EAST_VOLCONS
#endif
#if !defined NORTH_M2OBC && defined NORTH_VOLCONS
# undef NORTH_VOLCONS
#endif
#if !defined SOUTH_M2OBC && defined SOUTH_VOLCONS
# undef SOUTH_VOLCONS
#endif

#if defined WEST_VOLCONS  || defined EAST_VOLCONS  || \
    defined NORTH_VOLCONS || defined SOUTH_VOLCONS
# define OBC_VOLCONS
#endif

/*
** Activate assimilation switches.
*/

#if defined RANDOM_ESPERT
# define PERTURBATION
# define ESSE
#endif

#if defined ASSIMILATION_SSH || defined ASSIMILATION_SST   || \
    defined ASSIMILATION_T   || defined ASSIMILATION_UVsur || \
    defined ASSIMILATION_UV
# define ASSIMILATION
#endif
#if defined NUDGING_SSH || defined NUDGING_SST   || \
    defined NUDGING_T   || defined NUDGING_UVsur || \
    defined NUDGING_UV
# define NUDGING
#endif

/*
** Check if it is meaningful to write out time-averaged vertical
** mixing coefficients.
*/

#if !defined LMD_MIXING && !defined MY25_MIXING && !defined GLS_MIXING
# if defined AVERAGES
#  if defined AVERAGES_AKV
#    undef AVERAGES_AKV
#  endif
#  if defined AVERAGES_AKT
#    undef AVERAGES_AKT
#  endif
#  if defined AVERAGES_AKS && !defined SALINITY
#    undef AVERAGES_AKS
#  endif
# endif
#endif

/*
** Define internal flag indicating processing of input forcing
** NetCDF file.
*/

#ifdef SOLVE3D
# ifdef BULK_FLUXES
#  ifdef ANA_SMFLUX
#   undef ANA_SMFLUX
#  endif
#  ifdef ANA_STFLUX
#   undef ANA_STFLUX
#  endif
# endif
# if !defined ANA_BTFLUX   || \
    (!defined AIR_OCEAN    && !defined BULK_FLUXES   && !defined ANA_SMFLUX)   || \
    (!defined BULK_FLUXES  && !defined ANA_STFLUX)   || \
    ( defined SALINITY     && !defined ANA_SSFLUX)   || \
    ( defined BULK_FLUXES  && !defined LONGWAVE)     || \
    ( defined BULK_FLUXES  && !defined ANA_PAIR)     || \
    ( defined BULK_FLUXES  && !defined ANA_TAIR)     || \
    ( defined BULK_FLUXES  && !defined ANA_HUMIDITY) || \
    ( defined BULK_FLUXES  && !defined ANA_CLOUD)    || \
    ( defined BULK_FLUXES  && !defined ANA_RAIN)     || \
    ( defined BULK_FLUXES  && !defined ANA_SRFLUX)   || \
    ( defined LMD_SKPP     && !defined ANA_SRFLUX)   || \
    ( defined SOLAR_SOURCE && !defined ANA_SRFLUX)   || \
    ( defined BBL_MODEL    && !defined ANA_WWAVE)    || \
    ( defined BIOLOGY      && !defined ANA_SPFLUX)   || \
    ( defined BIOLOGY      && !defined ANA_BPFLUX)   || \
    ( defined SEDIMENT     && !defined ANA_SPFLUX)   || \
    ( defined SEDIMENT     && !defined ANA_BPFLUX)
#  define FRC_FILE
# endif
#else
# if !defined ANA_SMFLUX
#  define FRC_FILE
# endif
#endif

/*
** Activate internal biology option when using any type of biological
** module.
*/

#if defined BIO_FASHAM || defined ECOSIM
# define BIOLOGY
#endif

/*
** Define internal shortwave radiation option.  Undefine analytical
** shortwave option if not needed.
*/

#if defined LMD_SKPP     || defined SOLAR_SOURCE   || \
    defined BULK_FLUXES  || defined BIOLOGY
# define SHORTWAVE
#endif
#if !defined SHORTWAVE   && defined ANA_SRFLUX
# undef ANA_SRFLUX
#endif
#if !defined SHORTWAVE   && defined DIURNAL_SRFLUX
# undef DIURNAL_SRFLUX
#endif

/*
** Define internal clouds option.  Undefine analytical
** shortwave option if not needed.
*/

#if (defined BULK_FLUXES && defined LONGWAVE) || defined ECOSIM || \
    (defined ANA_SRFLUX  && defined ALBEDO)
# define CLOUDS
#endif
#if !defined CLOUDS && defined ANA_CLOUD
# undef ANA_CLOUD
#endif

/*
** Check if it is meaningful to write out momentum/tracer diagnostics
** and activate internal diagnostics option.
*/

#if !defined SOLVE3D || defined TS_FIXED
# if defined DIAGNOSTICS_TS
#   undef DIAGNOSTICS_TS
# endif
#endif
#if !defined BIO_FASHAM && defined DIAGNOSTICS_BIO
#  undef DIAGNOSTICS_BIO
#endif
#if defined DIAGNOSTICS_BIO || defined DIAGNOSTICS_TS || \
    defined DIAGNOSTICS_UV
# define DIAGNOSTICS
#endif
