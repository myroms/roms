/*
** Include file "cppdefs.h"
**********************************************************************
** Copyright (c) 2002 ROMS/TOMS Group, version 2.0                  **
************************************************* Hernan G. Arango ***
**                                                                  **
**  Choose the appropriate C-preprocessing options by using the     **
**  command #define to activate option or #undef to deactivate      **
**  option.                                                         **
**                                                                  **
**********************************************************************
*/

/*
**--------------------------------------------------------------------
**  Choose a pre-defined model application.  If building a new
**  application, choose a unique CPP flag for it and select the
**  appropriate configuration.
**--------------------------------------------------------------------
*/

#undef  ADRIATIC1       /* for Adriatic Sea Application, Grid #1 */
#undef  ADRIATIC2       /* for Adriatic Sea Application, Grid #2 */
#undef  BASIN           /* for Big Bad Basin Example */
#undef  BENCHMARK1      /* for Small Benchmark Test */
#undef  BENCHMARK2      /* for Medium Benchmark Test */
#undef  BENCHMARK3      /* for Big Benchmark Test */
#undef  BL_TEST         /* for Boundary Layers Test */
#undef  BISCAY          /* for Bay of Biscay Aplication */
#undef  CALCOFI         /* for Southern California Bight */
#undef  CANYON_A        /* for Canyon_A Example */
#undef  CANYON_B        /* for Canyon_B Example */
#undef  CBLAST          /* for CBLAST application */
#undef  DAMEE_4         /* for North Atlantic DAMEE Application Grid 4 */
#undef  DAMEE_5         /* for North Atlantic DAMEE Application Grid 5 */
#undef  DAMEE_6         /* for North Atlantic DAMEE Application Grid 6 */
#undef  DAMEE_7         /* for North Atlantic DAMEE Application Grid 7 */
#undef  DAMEE_9         /* for North Atlantic DAMEE Application Grid 9 */
#undef  DOUBLE_GYRE     /* for idealized double-gyre Example */
#undef  FLT_TEST        /* for Float Tracking Example */
#undef  GRAV_ADJ        /* for Graviational Adjustment Example */
#undef  KELVIN          /* for Kelvin wave test */
#undef  LAB_CANYON      /* for Lab Canyon, Polar Coordinates Example */
#undef  LMD_TEST        /* Test for LMD and KPP */
#undef  NATL            /* for high resolution North Atlantic Application */
#undef  NENA            /* North East North America Application */
#undef  NJ_BIGHT        /* for New Jersey Bight Application */
#undef  NPACIFIC        /* for North Pacific Application */
#undef  OVERFLOW        /* for Graviational/Overflow Example */
#undef  RIVERPLUME      /* for River plume Example */
#undef  SEAMOUNT        /* for Seamount Example */
#undef  SED_TEST1       /* for Suspended Sediment Test in a Channel */
#undef  SOLITON         /* for Equatorial Rossby Wave Example */
#define UPWELLING       /* for Upwelling Example */
#undef  USWEST          /* for US West Coast Application */
#undef  WEDDELL         /* for Idealized Weddell Sea Shelf Application */

# if defined MY_APPL

/*
**--------------------------------------------------------------------
**  Detailed description of all available CPP options.
**--------------------------------------------------------------------
**
**  Select model dynamics for MOMENTUM equations:
**  (The default advection is third-order upstream bias)
*/

#undef  UV_ADV          /* turn ON or OFF advection terms */
#undef  UV_C2ADVECTION  /* turn ON or OFF 2nd-order centered advection */
#undef  UV_C4ADVECTION  /* turn ON or OFF 4th-order centered advection */
#undef  UV_SADVECTION   /* turn ON or OFF splines vertical advection */
#undef  UV_VIS2         /* turn ON or OFF Laplacian horizontal mixing */
#undef  UV_VIS4         /* turn ON or OFF biharmonic horizontal mixing */
#undef  UV_COR          /* turn ON or OFF Coriolis term */
#undef  UV_PSOURCE      /* turn ON or OFF point Sources/Sinks */

/*
**  Select model dynamics for TRACER equations:
**  (The default horizontal and vertical advection is 4th-order centered)
*/

#undef  TS_A4HADVECTION /* define if 4th-order Akima horiz. advection */
#undef  TS_C2HADVECTION /* define if 2nd-order centered horiz. advection */
#undef  TS_C4HADVECTION /* define if 4th-order centered horiz. advection */
#undef  TS_U3HADVECTION /* define if 3rd-order upstream horiz. advection */
#undef  TS_A4VADVECTION /* define if 4th-order Akima vertical advection */
#undef  TS_C2VADVECTION /* define if 2nd-order centered vertical advection */
#undef  TS_C4VADVECTION /* define if 4th-order centered vertical advection */
#undef  TS_SVADVECTION  /* define if splines vertical advection */
#undef  TS_DIF2         /* turn ON or OFF Laplacian horizontal mixing */
#undef  TS_DIF4         /* turn ON or OFF biharmonic horizontal mixing */
#undef  SALINITY        /* define if using salinity */
#undef  NONLIN_EOS      /* define if using nonlinear equation of state */
#undef  DIAGNOSTIC      /* define if diagnostic simulation */
#undef  QCORRECTION     /* define if net heat flux correction */
#undef  SCORRECTION     /* define if freshwater flux correction */
#undef  SOLAR_SOURCE    /* define solar radiation source term */
#undef  SRELAXATION     /* define if salinity relaxation as a freshwater flux*/
#undef  TS_PSOURCE      /* turn ON or OFF point Sources/Sinks */

/*
**  Select pressure gradient algorithm.  If no option is selected, the
**  pressure gradient term is computed using standard density Jacobian
**  algorithm.  Notice that there are two quartic pressure Jacobian
**  options.  They differ on how the WENO reconsicliation step is done
**  and in the monotonicity constraining algorithms.
*/

#undef  PJ_GRADP        /* Finite volume Pressure Jacobian (Lin, 1997) */
#undef  PJ_GRADPQ2      /* Quartic 2 Pressure Jacobian (Shchepetkin, 2000) */
#undef  PJ_GRADPQ4      /* Quartic 4 Pressure Jacobian (Shchepetkin, 2000) */
#undef  DJ_GRADPS       /* Splines density  Jacobian (Shchepetkin, 2000) */
#undef  WJ_GRADP        /* Weighted density Jacobian (Song, 1998) */

/*
**  Select surface model fluxes formutalion via atmospheric boundary
**  layer (Fairall et al, 1996).
*/

#undef  BULK_FLUXES     /* turn ON or OFF bulk fluxes computation */
#undef  COOL_SKIN       /* turn ON or OFF cool skin correction */
#undef  LONGWAVE        /* Compute net longwave radiation internally */

/*
**  Select options for shortwave radiation.  The shortwave radiation can be
**  computed using the global albedo equation with a cloud correction.
**  Alternatively, input shortwave radiation data computed from average
**  data (with snapshots greater or equal than 24 hours) can be modulated
**  by the local diurnal cycle which is a function longitude, latitude and
**  day-of-year.
*/

#undef  ALBEDO          /* use albedo equation for shortwave radiation */
#undef  DIURNAL_SRFLUX  /* impose shortwave radiation local diurnal cycle */

/*
**  Select model configuration:
*/

#undef  SERIAL          /* define if serial code configuration */
#undef  SOLVE3D         /* define if solving 3D primitive equations */
#undef  CURVGRID        /* define if using  curvilinear coordinate grid*/
#undef  MASKING         /* define if there is land in the domain */
#undef  BODYFORCE       /* define if applying stresses as bodyforces */
#undef  PROFILE         /* define if time profiling */
#undef  AVERAGES        /* define if writing out time-averaged data */
#undef  AVERAGES_AKV    /* define if writing out time-averaged AKv */
#undef  AVERAGES_AKT    /* define if writing out time-averaged AKt */
#undef  AVERAGES_AKS    /* define if writing out time-averaged AKs */
#undef  ICESHELF        /* define if including ice shelf cavities */
#undef  SPHERICAL       /* define if analytical spherical grid */
#undef  STATIONS        /* define if writing out station data */
#undef  STATIONS_CGRID  /* define if extracting data at native C-grid */

/*
** Lagrangian drifters options
*/

#undef  FLOATS          /* define if simulated Lagrangian drifters */
#undef  FLOAT_VWALK     /* define if vertical random walk */

/*
** Activate conservative, parabolic spline reconstruction in the
** vertical.  Notice that there also options (see above) for vertical
** advection of momentum and tracers using splines.
*/

#undef  SPLINES         /* turn ON or OFF parabolic splines reconstruction */

/*
**  Select analytical fields configuration: define if using any of the
**  following options.  Set the appropriate analytical expression in
**  file "analytical.F".
*/

#undef  ANA_BIOLOGY     /* analytical biology initial conditions */
#undef  ANA_BPFLUX      /* analytical bottom passive tracers fluxes */
#undef  ANA_BSEDIM      /* analytical bottom sediment grain size & density */
#undef  ANA_BSFLUX      /* analytical bottom salinity flux */
#undef  ANA_BTFLUX      /* analytical bottom temperature flux */
#undef  ANA_CLOUD       /* analytical cloud fraction */
#undef  ANA_DIAG        /* Customized diagnostics */
#undef  ANA_FSOBC       /* analytical free-surface boundary conditions */
#undef  ANA_GRID        /* analytical model grid set-up */
#undef  ANA_HUMIDITY    /* analytical surface air humidity */
#undef  ANA_INITIAL     /* analytical initial conditions */
#undef  ANA_M2CLIMA     /* analytical 2D momentum climatology */
#undef  ANA_M2OBC       /* analytical 2D momentum boundary conditions */
#undef  ANA_M3CLIMA     /* analytical 3D momentum climatology */
#undef  ANA_M3OBC       /* analytical 3D momentum boundary conditions */
#undef  ANA_PAIR        /* analytical surface air pressure */
#undef  ANA_PSOURCE     /* analytical point Sources/Sinks */
#undef  ANA_RAIN        /* analytical rain fall rate */
#undef  ANA_SEDIMENT    /* analytical sediment initial conditions */
#undef  ANA_SMFLUX      /* analytical surface momentum stress */
#undef  ANA_SPFLUX      /* analytical surface passive tracers fluxes */
#undef  ANA_SPINNING    /* analytical time-varying rotation force */
#undef  ANA_SRFLUX      /* analytical surface shortwave radiation flux */
#undef  ANA_SSFLUX      /* analytical surface salinity flux */
#undef  ANA_SSH         /* analytical sea surface height */
#undef  ANA_SSS         /* analytical sea surface salinity */
#undef  ANA_SST         /* analytical SST and dQdSST */
#undef  ANA_STFLUX      /* analytical surface temperature flux */
#undef  ANA_TAIR        /* analytical surface air temperature */
#undef  ANA_TCLIMA      /* analytical tracers climatology */
#undef  ANA_TOBC        /* analytical tracers boundary conditions */
#undef  ANA_VMIX        /* analytical vertical mixing coefficients */
#undef  ANA_WINDS       /* analytical surface winds */
#undef  ANA_WWAVE       /* analytical wind induced waves */

/*
**  Select options for horizontal mixing of MOMENTUM:
*/

#undef  VISC_GRID       /* viscosity coefficient scaled by grid size */
#undef  MIX_S_UV        /* mixing along constant S-surfaces */
#undef  MIX_GEO_UV      /* mixing on geopotential (constant Z) surfaces */

/*
**  Select options for horizontal mixing of TRACERS:
*/

#undef  DIFF_GRID       /* diffusion coefficient scaled by grid size */
#undef  MIX_S_TS        /* mixing along constant S-surfaces */
#undef  MIX_GEO_TS      /* mixing on geopotential (constant Z) surfaces */
#undef  MIX_ISO_TS      /* mixing on epineutral (constant RHO) surfaces */

/*
**  Select vertical turbulent mixing scheme for MOMENTUM and TRACERS
**  (activate only one closure):
*/

#undef  BVF_MIXING      /* Activate Brunt-Vaisala frequency mixing */
#undef  GLS_MIXING      /* Activate Generic Length-Scale mixing */
#undef  MY25_MIXING     /* Activate Mellor/Yamada Level-2.5 closure */
#undef  LMD_MIXING      /* Activate Large/McWilliams/Doney interior closure */

/*
**  Select options for the Generic Length-Scale closure:
**  (The default advection is third-order upstream bias, G-Scheme)
*/

#undef  N2S2_HORAVG     /* Activate horizontal smoothing of buoyancy/shear */
#undef  CANUTO_A        /* Canuto A-stability function formulation */
#undef  CANUTO_B        /* Canuto B-stability function formulation */
#undef  KANTHA_CLAYSON  /* Kantha and Clayson stability function formulation */
#undef  K_C2ADVECTION   /* turn ON or OFF 2nd-order centered advection */
#undef  K_C4ADVECTION   /* turn ON or OFF 4th-order centered advection */

/*
**  Select options for the Mellor/Yamada level 2.5 closure:
**  (The default advection is third-order upstream bias, G-Scheme)
*/

#undef  N2S2_HORAVG     /* Activate horizontal smoothing of buoyancy/shear */
#undef  KANTHA_CLAYSON  /* Kantha and Clayson stability function formulation */
#undef  K_C2ADVECTION   /* turn ON or OFF 2nd-order centered advection */
#undef  K_C4ADVECTION   /* turn ON or OFF 4th-order centered advection */

/*
**  Select options for the Large/McWilliams/Doney interior mixing:
*/

# ifdef LMD_MIXING
#undef  LMD_RIMIX       /* Add diffusivity due to shear instability */
#undef  LMD_CONVEC      /* Add convective mixing due to shear instability */
#undef  LMD_DDMIX       /* Add double-diffusive mixing */
# endif

/*
**  Select Large/McWilliams/Doney Oceanic Planetary Boundary Layer scheme
**  and its associated options: local K-Profile Parameterization (KPP)
*/

#undef  LMD_SKPP        /* turn ON or OFF surface boundary layer KPP mixing */
#undef  LMD_BKPP        /* turn ON or OFF bottom  boundary layer KPP mixing */
#undef  LMD_NONLOCAL    /* turn ON or OFF nonlocal transport */

/*
**  If not SPLINES, activate smoothing of Richardson number.
*/

#undef  RI_HORAVG       /* Activate horizontal smoothing */
#undef  RI_VERAVG       /* Activate vertical smoothing */

/*
**  Select options for bottom boundary layer closure.
*/

#undef  CS_BBL         /* turn ON or OFF Chris Sherwood BBL closure */
#undef  MB_BBL         /* turn ON or OFF Mainte Blaas BBL closure */
#undef  SG_BBL         /* turn ON or OFF Styles and Glenn (2000) BBL closure */
#undef  SG_ZNOT        /* activate internal computation of bottom roughness */
#undef  SG_LOGINT      /* activate logarithmic interpolation of (Ur,Vr) */

/*
**  Select lateral boundary options.  Usually, select ONE option at each
**  boundary edge for free-surface, 2D momentum, 3D momentum, and tracers.
**  The turbulent kineric energy (TKE) conditions are only activated for
**  the Mellor-Yamada 2.5 closure. If open boundary radiation conditions,
**  an additional option can be activated at each boundary edge to include
**  a passive/active nudging term with weak/strong values for outflow/inflow.
*/

#undef  SPONGE            /* activate areas of enhanced viscosity/diffusion */
#undef  OBC_VOLCONS       /* enforce mass conservation on open boundaries */

#undef  EW_PERIODIC       /* East-West periodic boundaries */
#undef  NS_PERIODIC       /* North-South periodic boundaries */

#undef  EASTERN_WALL      /* Eastern edge, closed wall condition */
#undef  WESTERN_WALL      /* Western edge, closed wall condition */
#undef  NORTHERN_WALL     /* Northern edge, closed wall condition */
#undef  SOUTHERN_WALL     /* Southern edge, closed wall condition */

#undef  RADIATION_2D      /* Tangential phase speed in radiation conditions */

#undef  EAST_FSCHAPMAN    /* Eastern edge, free-surface, Chapman condition */
#undef  EAST_FSGRADIENT   /* Eastern edge, free-surface, gradient condition */
#undef  EAST_FSRADIATION  /* Eastern edge, free-surface, radiation condition */
#undef  EAST_FSNUDGING    /* Eastern edge, free-surface, passive/active term */
#undef  EAST_FSCLAMPED    /* Eastern edge, free-surface, clamped condition */
#undef  EAST_M2FLATHER    /* Eastern edge, 2D momentum, Flather condition */
#undef  EAST_M2GRADIENT   /* Eastern edge, 2D momentum, gradient condition */
#undef  EAST_M2RADIATION  /* Eastern edge, 2D momentum, radiation condition */
#undef  EAST_M2NUDGING    /* Eastern edge, 2D momentum, passive/active term */
#undef  EAST_M2CLAMPED    /* Eastern edge, 2D momentum, clamped condition */
#undef  EAST_M3GRADIENT   /* Eastern edge, 3D momentum, gradient condition */
#undef  EAST_M3RADIATION  /* Eastern edge, 3D momentum, radiation condition */
#undef  EAST_M3NUDGING    /* Eastern edge, 3D momentum, passive/active term */
#undef  EAST_M3CLAMPED    /* Eastern edge, 3D momentum, clamped condition */
#undef  EAST_KGRADIENT    /* Eastern edge, TKE fields, gradient condition */
#undef  EAST_KRADIATION   /* Eastern edge, TKE fields, radiation condition */
#undef  EAST_TGRADIENT    /* Eastern edge, tracers, gradient condition */
#undef  EAST_TRADIATION   /* Eastern edge, tracers, radiation condition */
#undef  EAST_TNUDGING     /* Eastern edge, tracers, passive/active term */
#undef  EAST_TCLAMPED     /* Eastern edge, tracers, clamped condition */

#undef  WEST_FSCHAPMAN    /* Western edge, free-surface, Chapman condition */
#undef  WEST_FSGRADIENT   /* Western edge, free-surface, gradient condition */
#undef  WEST_FSRADIATION  /* Western edge, free-surface, radiation condition */
#undef  WEST_FSNUDGING    /* Western edge, free-surface, passive/active term */
#undef  WEST_FSCLAMPED    /* Western edge, free-surface, clamped condition */
#undef  WEST_M2FLATHER    /* Western edge, 2D momentum, Flather condition */
#undef  WEST_M2GRADIENT   /* Western edge, 2D momentum, gradient condition */
#undef  WEST_M2RADIATION  /* Western edge, 2D momentum, radiation condition */
#undef  WEST_M2NUDGING    /* Western edge, 2D momentum, passive/active term */
#undef  WEST_M2CLAMPED    /* Western edge, 2D momentum, clamped condition */
#undef  WEST_M3GRADIENT   /* Western edge, 3D momentum, gradient condition */
#undef  WEST_M3RADIATION  /* Western edge, 3D momentum, radiation condition */
#undef  WEST_M3NUDGING    /* Western edge, 3D momentum, passive/active term */
#undef  WEST_M3CLAMPED    /* Western edge, 3D momentum, clamped condition */
#undef  WEST_KGRADIENT    /* Western edge, TKE fields, gradient condition */
#undef  WEST_KRADIATION   /* Western edge, TKE fields, radiation condition */
#undef  WEST_TGRADIENT    /* Western edge, tracers, gradient condition */
#undef  WEST_TRADIATION   /* Western edge, tracers, radiation condition */
#undef  WEST_TNUDGING     /* Western edge, tracers, passive/active term */
#undef  WEST_TCLAMPED     /* Western edge, tracers, clamped condition */

#undef  NORTH_FSCHAPMAN   /* Northern edge, 2D momentum, Chapman condition */
#undef  NORTH_FSGRADIENT  /* Northern edge, free-surface, gradient condition */
#undef  NORTH_FSRADIATION /* Northern edge, free-surface, radiation condition*/
#undef  NORTH_FSNUDGING   /* Northern edge, free-surface, passive/active term*/
#undef  NORTH_FSCLAMPED   /* Northern edge, free-surface, clamped condition */
#undef  NORTH_M2FLATHER   /* Northern edge, 2D momentum, Flather condition */
#undef  NORTH_M2GRADIENT  /* Northern edge, 2D momentum, gradient condition */
#undef  NORTH_M2RADIATION /* Northern edge, 2D momentum, radiation condition */
#undef  NORTH_M2NUDGING   /* Northern edge, 2D momentum, passive/active term */
#undef  NORTH_M2CLAMPED   /* Northern edge, 2D momentum, clamped condition */
#undef  NORTH_M3GRADIENT  /* Northern edge, 3D momentum, gradient condition */
#undef  NORTH_M3RADIATION /* Northern edge, 3D momentum, radiation condition */
#undef  NORTH_M3NUDGING   /* Northern edge, 3D momentum, passive/active term */
#undef  NORTH_M3CLAMPED   /* Northern edge, 3D momentum, clamped condition */
#undef  NORTH_KGRADIENT   /* Northern edge, TKE fields, gradient condition */
#undef  NORTH_KRADIATION  /* Northern edge, TKE fields, radiation condition */
#undef  NORTH_TGRADIENT   /* Northern edge, tracers, gradient condition */
#undef  NORTH_TRADIATION  /* Northern edge, tracers, radiation condition */
#undef  NORTH_TNUDGING    /* Northern edge, tracers, passive/active term */
#undef  NORTH_TCLAMPED    /* Northern edge, tracers, clamped condition */

#undef  SOUTH_FSCHAPMAN   /* Southern edge, 2D momentum, Chapman condition */
#undef  SOUTH_FSGRADIENT  /* Southern edge, free-surface, gradient condition */
#undef  SOUTH_FSRADIATION /* Southern edge, free-surface, radiation condition*/
#undef  SOUTH_FSNUDGING   /* Southern edge, free-surface, passive/active term*/
#undef  SOUTH_FSCLAMPED   /* Southern edge, free-surface, clamped condition */
#undef  SOUTH_M2FLATHER   /* Southern edge, 2D momentum, Flather condition */
#undef  SOUTH_M2GRADIENT  /* Southern edge, 2D momentum, gradient condition */
#undef  SOUTH_M2RADIATION /* Southern edge, 2D momentum, radiation condition */
#undef  SOUTH_M2NUDGING   /* Southern edge, 2D momentum, passive/active term */
#undef  SOUTH_M2CLAMPED   /* Southern edge, 2D momentum, clamped condition */
#undef  SOUTH_M3GRADIENT  /* Southern edge, 3D momentum, gradient condition */
#undef  SOUTH_M3RADIATION /* Southern edge, 3D momentum, radiation condition */
#undef  SOUTH_M3NUDGING   /* Southern edge, 3D momentum, passive/active term */
#undef  SOUTH_M3CLAMPED   /* Southern edge, 3D momentum, clamped condition */
#undef  SOUTH_KGRADIENT   /* Southern edge, TKE fields, gradient condition */
#undef  SOUTH_KRADIATION  /* Southern edge, TKE fields, radiation condition */
#undef  SOUTH_TGRADIENT   /* Southern edge, tracers, gradient condition */
#undef  SOUTH_TRADIATION  /* Southern edge, tracers, radiation condition */
#undef  SOUTH_TNUDGING    /* Southern edge, tracers, passive/active term */
#undef  SOUTH_TCLAMPED    /* Southern edge, tracers, clamped condition */

/*
**  Process tidal forcing for desired tidal component, classified by
**  period.  The tidal forcing is computed for the full horizontal grid.
**  If requested, the tidal forcing is added to the processed (read and
**  time-interpolated) open boundary data. Also, if applicable, the tidal
**  forcing is added to the climatological array.
*/

#undef  SSH_TIDES       /* turn on computation of tidal elevation */
#undef  UV_TIDES        /* turn on computation of tidal currents */
#undef  ADD_FSOBC       /* Add tidal elevation to processed OBC data */
#undef  ADD_M2OBC       /* Add tidal currents  to processed OBC data */

/*
**  Turn ON or OFF options for reading and processing of climatological
**  fields.  The nudging of climatology data is primarily used in sponge
**  areas.
*/

#undef  M2CLIMATOLOGY   /* Processing of 2D momentum climatology */
#undef  M3CLIMATOLOGY   /* Processing of 3D momentum climatology */
#undef  TCLIMATOLOGY    /* Processing of tracer climatology */
#undef  ZCLIMATOLOGY    /* Processing of SSH climatology */

#undef  M2CLM_NUDGING   /* Nudging of 2D momentum climatology */
#undef  M3CLM_NUDGING   /* Nudging of 3D momentum climatology */
#undef  TCLM_NUDGING    /* Nudging of tracer climatology */
#undef  ZCLM_NUDGING    /* Nudging of SSH climatology */

/*
**  Select data assimilation options. The assimilation is via OI and
**  intermittent. Nudging is continuous and observations are time-
**  interpolated. If applicable, choose only one option for each
**  field update: assimilation or nudging.
*/

#undef  ASSIMILATION_SSH   /* assimilation of SSH observations */
#undef  ASSIMILATION_SST   /* assimilation of SST observations */
#undef  ASSIMILATION_T     /* assimilation of tracer observations */
#undef  ASSIMILATION_UVsur /* assimilation of surface current observations */
#undef  ASSIMILATION_UV    /* assimilation of horizontal current observations*/
#undef  UV_BAROCLINIC      /* assimilation of only baroclinic currents */
#undef  NUDGING_SSH        /* nudging of SSH observations */
#undef  NUDGING_SST        /* nudging of SST observations */
#undef  NUDGING_T          /* nudging of tracer observations */
#undef  NUDGING_UVsur      /* nudging of surface current observations */
#undef  NUDGING_UV         /* nudging of horizontal currents observations */
#undef  RANDOM_ESPERT      /* random error subspace perturbations */

/*
**  Biological model options.
*/

#undef  BIO_FASHAM         /* Fasham type nitrogen-based model */

/*
**  Sediment transport model options.
*/

#undef  SEDIMENT           /* Activate sediment transport model */

/*
**  NetCDF IO options.
*/

#undef  NO_WRITE_GRID      /* define if not writing grid arrays */
#undef  READ_WATER         /* define if only reading water points data */
#undef  WRITE_WATER        /* define if only writing water points data */
#undef  RST_SINGLE         /* define if single precision restart fields */
#undef  OUT_DOUBLE         /* define if double precision output fields */

/*
**--------------------------------------------------------------------
**  Set pre-defined configuration flags for model test problems.
**--------------------------------------------------------------------
*/

# elif defined ADRIATIC1 || defined ADRIATIC2

/*
**  Adriatic Sea Application, Grid 1.
*/

#define UV_ADV
#define UV_SADVECTION
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_SVADVECTION
#undef  TS_C4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS 
#undef  UV_PSOURCE
#undef  TS_PSOURCE
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define CURVGRID
#define SOLVE3D
#define MASKING
#define SPLINES
#define QCORRECTION
#define SCORRECTION
#define AVERAGES
#ifdef AVERAGES
# define AVERAGES_AKV
# define AVERAGES_AKT
# define AVERAGES_AKS
#endif
#undef LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#endif
#undef  MY25_MIXING
#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
# undef K_C2ADVECTION
# undef K_C4ADVECTION
#endif
#define GLS_MIXING
#ifdef MY25_MIXING
# define N2S2_HORAVG
# undef KANTHA_CLAYSON
#endif
#define EASTERN_WALL
#define WESTERN_WALL
#define NORTHERN_WALL
#define SOUTHERN_WALL
#ifndef SOUTHERN_WALL
# define OBC_VOLCONS
# define RADIATION_2D
# define SOUTH_FSGRADIENT
# define SOUTH_M2RADIATION
# define SOUTH_M3RADIATION
# undef  SOUTH_M3NUDGING
# define SOUTH_TRADIATION
# define SOUTH_TNUDGING
# define SPONGE
#endif
#undef  M3CLIMATOLOGY
#undef  M3CLM_NUDGING
#define TCLIMATOLOGY
#define TCLM_NUDGING
#define ANA_BTFLUX
#define ANA_BSFLUX
#undef  ANA_SMFLUX
#undef  ANA_STFLUX
#undef  ANA_SSFLUX
#define ANA_SRFLUX
#ifdef ANA_SRFLUX
# define ALBEDO
# define ANA_CLOUD
# define ANA_HUMIDITY
# define ANA_TAIR
#endif

# elif defined BASIN

/*
**  Options for Big Bad Basin Example:
*/

#define UV_ADV
#define UV_COR
#define UV_VIS4
#define MIX_S_UV
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define SOLVE3D
#define SPLINES
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL
#define BODYFORCE
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX

# elif defined BENCHMARK1 || defined BENCHMARK2 || defined BENCHMARK3

/*
**  Options Benchmark Test: Idealize Southern Ocean.
*/

#define UV_ADV
#define UV_COR
#define DJ_GRADPS
#define UV_VIS2
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define CURVGRID
#define SOLVE3D
#define SPLINES
#define EW_PERIODIC
#define SOUTHERN_WALL
#define NORTHERN_WALL
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#endif
#define BULK_FLUXES
#ifdef BULK_FLUXES
# define ANA_WINDS
# define ANA_TAIR
# define ANA_PAIR
# define ANA_HUMIDITY
# define ANA_RAIN
# define LONGWAVE
# define ANA_CLOUD
#endif
#define SPHERICAL
#define ANA_GRID
#define ANA_INITIAL
#define ALBEDO
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined BL_TEST

/*
**  Options for Boundary Layers Test.
*/

#define UV_ADV
#define UV_SADVECTION
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define SPLINES
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS
#define STATIONS
#define SOLVE3D
#define WESTERN_WALL
#define NS_PERIODIC
#define RADIATION_2D
#define EAST_FSGRADIENT
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#define EAST_KRADIATION
#define EAST_TRADIATION
#define OBC_VOLCONS
#undef  MY25_MIXING
#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
#endif
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
# define LMD_DDMIX
#endif
#define BULK_FLUXES
#ifdef BULK_FLUXES
# define LONGWAVE
# define ANA_CLOUD
# define ANA_HUMIDITY
# define ANA_PAIR
# define ANA_TAIR
# define ANA_RAIN
# define ANA_WINDS
#else
# define ANA_SMFLUX
# define ANA_STFLUX
#endif
#define SG_BBL
#ifdef SG_BBL
# define SG_ZNOT
# define ANA_BSEDIM
# define ANA_WWAVE
#endif
#define ANA_GRID
#define ANA_INITIAL
#define ALBEDO
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined BISCAY

/*
**  Options for Bay of Biscay Application.
*/

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define VISC_GRID
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define DIFF_GRID
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define SOLVE3D
#define CURVGRID
#define MASKING
#define SPLINES
#define QCORRECTION
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS
#define EASTERN_WALL
#define OBC_VOLCONS
#define SPONGE
#define NORTH_FSGRADIENT
#define NORTH_M2GRADIENT
#define NORTH_M3RADIATION
#define NORTH_TRADIATION
#define SOUTH_FSGRADIENT
#define SOUTH_M2GRADIENT
#define SOUTH_M3RADIATION
#define SOUTH_TRADIATION
#define WEST_FSGRADIENT
#define WEST_M2GRADIENT
#define WEST_M3RADIATION
#define WEST_TRADIATION
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#endif
#define TCLIMATOLOGY
#define TCLM_NUDGING
#define ZCLIMATOLOGY
#define ZCLM_NUDGING
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined CALCOFI

/*
**  Options for Southern California Bight Application.
*/

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define CURVGRID
#define MASKING
#define SOLVE3D
#define SPLINES
#define QCORRECTION
#define EASTERN_WALL
#define OBC_VOLCONS
#define SPONGE
#define OBC_VOLCONS
#define M2CLIMATOLOGY
#define M3CLIMATOLOGY
#define TCLIMATOLOGY
#define M2CLM_NUDGING
#define M3CLM_NUDGING
#define TCLM_NUDGING
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL
#define WEST_FSGRADIENT
#define WEST_M2RADIATION
#define WEST_M2NUDGING
#define WEST_M3RADIATION
#define WEST_M3NUDGING
#define WEST_TRADIATION
#define WEST_TNUDGING
#define SOUTH_FSGRADIENT
#define SOUTH_M2RADIATION
#define SOUTH_M2NUDGING
#define SOUTH_M3RADIATION
#define SOUTH_M3NUDGING
#define SOUTH_TRADIATION
#define SOUTH_TNUDGING
#define NORTH_FSGRADIENT
#define NORTH_M2RADIATION
#define NORTH_M2NUDGING
#define NORTH_M3RADIATION
#define NORTH_M3NUDGING
#define NORTH_TRADIATION
#define NORTH_TNUDGING
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_NONLOCAL
# define LMD_SKPP
#endif
#undef  ANA_SRFLUX
#undef  ALBEDO
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined CANYON_A

/*
**  Options for Canyon A Example.
*/

#define UV_ADV
#define UV_VIS2
#define UV_COR
#define EW_PERIODIC
#define BODYFORCE
#define ANA_DIAG
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX

# elif defined CANYON_B

/*
**  Options for Canyon B Example.
*/

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define SOLVE3D
#define SPLINES
#define EW_PERIODIC
#define ANA_DIAG
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX
#define ANA_VMIX

# elif defined CBLAST

/*
**  Options for CBLAST Application.
*/

#define	UV_ADV
#define	UV_COR
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
# define SG_ZNOT
# define ANA_BSEDIM
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
#undef  OBC_VOLCONS
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

# elif defined DAMEE_4 || defined DAMEE_5 || defined DAMEE_6 || \
       defined DAMEE_7 || defined DAMEE_9

/*
**  Options for North Atlantic DAMEE Application.
*/

#define UV_ADV
#define UV_COR
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4HADVECTION
#define NONLIN_EOS
#define SALINITY
#define SOLVE3D
#define MASKING
#define READ_WATER
#define WRITE_WATER
#define SPLINES
#define QCORRECTION
#define SRELAXATION
#define CURVGRID
#define STATIONS
#define FLOATS
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#endif
#define TCLIMATOLOGY
#define TCLM_NUDGING
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined DOUBLE_GYRE

/*
**  Options for idealized double-gyre example.
*/

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#undef  DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4HADVECTION
#define TS_DIF2
#define MIX_S_TS
#define SOLVE3D
#undef  SPLINES
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL
#define TCLIMATOLOGY
#define TCLM_NUDGING
#define ANA_GRID
#define ANA_INITIAL
#define ANA_TCLIMA
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX

# elif defined FLT_TEST

/*
**  Options for Float Tracking Example.
*/

#define UV_ADV
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
#undef  SOLVE3D
#ifdef SOLVE3D
# define DJ_GRADPS
# define TS_A4HADVECTION
# define TS_A4VADVECTION
# define BODYFORCE
# define SPLINES
# define ANA_BTFLUX
# define ANA_STFLUX
#endif

# elif defined GRAV_ADJ

/*
**  Options for Gravitational Adjustment Example.
*/

#define UV_ADV
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define TS_DIF2
#define MIX_S_TS
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX

# elif defined KELVIN

/*
**  Kelvin wave test.
*/

#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#define UV_VIS2
#define TS_DIF2
#define SOLVE3D
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define RADIATION_2D
#define EAST_FSRADIATION
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#define EAST_TRADIATION
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#undef  WEST_FSCLAMPED
#undef  WEST_M2CLAMPED
#define WEST_M3RADIATION
#define WEST_TRADIATION
#define ANA_GRID
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SRFLUX
#define ANA_BTFLUX

# elif defined LAB_CANYON

/*
**  Options for Lab Canyon Example (Polar Coordinates).
*/

#define UV_COR
#define UV_ADV
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define SOLVE3D
#define CURVGRID
#define AVERAGES
#define SPLINES
#define NS_PERIODIC
#define EASTERN_WALL
#define WESTERN_WALL
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX

# elif defined LMD_TEST

/*
**  Test for LMD and KPP.
*/

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define NONLIN_EOS
#define SALINITY
#define AVERAGES
#define AVERAGES_AKS
#define AVERAGES_AKT
#define AVERAGES_AKV
#define STATIONS
#define SOLVE3D
#define SPLINES
#define EW_PERIODIC
#define NS_PERIODIC
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_DDMIX
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
#endif
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined NATL

/*
**  High Resolution North Atlantic Application
*/

#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPLINES
#define STATIONS
#define MASKING
#define WRITE_WATER
#define NO_WRITE_GRID
#define AVERAGES
#define SRELAXATION
#define QCORRECTION
#define SOLAR_SOURCE
#define ANA_BSFLUX
#define ANA_BTFLUX
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#endif
#define TCLIMATOLOGY
#define TCLM_NUDGING
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL

# elif defined NENA

/*
** North East North America Application.
*/

#define UV_ADV
#define UV_SADVECTION
#define DJ_GRADPS 
#define UV_COR
#define UV_VIS2
#define UV_PSOURCE
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define TS_PSOURCE
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPLINES
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
#endif
#undef  GLS_MIXING
#define BIO_FASHAM
#ifdef BIO_FASHAM
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif
#undef  M2CLIMATOLOGY
#undef  M3CLIMATOLOGY
#undef  TCLIMATOLOGY
#undef  ZCLIMATOLOGY
#undef  OBC_VOLCONS
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3CLAMPED
#define EAST_TCLAMPED
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3CLAMPED
#define WEST_TCLAMPED
#define NORTHERN_WALL
#define SOUTH_FSCHAPMAN
#define SOUTH_M2FLATHER
#define SOUTH_M3CLAMPED
#define SOUTH_TCLAMPED
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined NJ_BIGHT

/*
**  Options for New Jersey Bight Application.
*/

#define UV_ADV
#define UV_SADVECTION
#define UV_COR
#undef  UV_PSOURCE
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#undef  TS_A4HADVECTION
#undef  TS_A4VADVECTION
#undef  TS_DIF2
#undef  MIX_GEO_TS
#undef  TS_PSOURCE
#define SOLAR_SOURCE
#define NONLIN_EOS
#define SALINITY
#define CURVGRID
#define SOLVE3D
#define MASKING
#define SPLINES
#define AVERAGES
#ifdef AVERAGES
# define AVERAGES_AKV
# define AVERAGES_AKT
# define AVERAGES_AKS
#endif
#define STATIONS
#undef  FLOATS
#undef  ASSIMILATION_UVsur
#undef  ASSIMILATION_T
#undef  NUDGING_T
#undef  NUDGING_UVsur
#define WESTERN_WALL
#define NORTHERN_WALL
#define RADIATION_2D
#define EAST_M3RADIATION
#define EAST_TRADIATION
#define SOUTH_M3RADIATION
#define SOUTH_TRADIATION
#define MY25_MIXING
#ifdef MY25_MIXING
# define N2S2_HORAVG
# define KANTHA_CLAYSON
#endif
#undef LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# undef  LMD_BKPP
# define LMD_NONLOCAL
#endif
#define BULK_FLUXES
#ifdef BULK_FLUXES
# define ANA_RAIN
# define LONGWAVE
# ifdef LONGWAVE
#  define ANA_CLOUD
# endif
#endif
#undef SG_BBL
#ifdef SG_BBL
# define SG_ZNOT
# define ANA_BSEDIM
# define ANA_WWAVE
#endif
#define SSH_TIDES
#ifdef SSH_TIDES
# define ANA_FSOBC
# define EAST_FSCHAPMAN
# define SOUTH_FSCHAPMAN
#else
# define EAST_FSGRADIENT
# define SOUTH_M2RADIATION
#endif
#define UV_TIDES
#ifdef UV_TIDES
# define ANA_M2OBC
# define EAST_M2FLATHER
# define SOUTH_M2FLATHER
#else
# define EAST_M2RADIATION
# define SOUTH_FSGRADIENT
#endif
#if defined SSH_TIDES || defined UV_TIDES
# undef  OBC_VOLCONS
#else
# define OBC_VOLCONS
#endif
#undef  BIO_FASHAM
#ifdef BIO_FASHAM
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
#endif
#define ANA_SRFLUX
#define ALBEDO
#define ANA_SMFLUX
#define ANA_SSFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined NPACIFIC

/*
**  Options for North Pacific Application.
*/

#define UV_ADV
#define UV_COR
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF4
#define MIX_GEO_TS
#define DIFF_GRID
#define NONLIN_EOS
#define SALINITY
#define SOLVE3D
#define CURVGRID
#define MASKING
#define SPLINES
#define QCORRECTION
#define AVERAGES
#define EASTERN_WALL
#define WESTERN_WALL
#define SOUTHERN_WALL
#define NORTHERN_WALL
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_NONLOCAL
#endif
#define TCLIMATOLOGY
#define TCLM_NUDGING
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined OVERFLOW

/*
**  Options for Gravitational/Overflow Example.
*/

#define UV_ADV
#define UV_COR
#define UV_VIS2
#define MIX_S_UV
#define DJ_GRADPS
#define TS_U3HADVECTION
#define TS_SVADVECTION
#define TS_DIF2
#define MIX_ISO_TS
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX

# elif defined RIVERPLUME

/*
**  River Plume example.
*/

#define UV_ADV
#define UV_COR
#define UV_PSOURCE
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define TS_PSOURCE
#define NONLIN_EOS
#define SALINITY
#define MASKING
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define AVERAGES_AKV
#define AVERAGES_AKT
#define AVERAGES_AKS
#define NS_PERIODIC
#define WESTERN_WALL
#define EASTERN_WALL
#undef  MY25_MIXING
#define LMD_MIXING
#ifdef LMD_MIXING
# define LMD_RIMIX
# define LMD_CONVEC
# define LMD_SKPP
# define LMD_BKPP
# define LMD_NONLOCAL
#endif
#define ANA_GRID
#define ANA_INITIAL
#define ANA_PSOURCE
#define ANA_SMFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined SEAMOUNT

/*
**  Options for Seamount Example.
*/

#define UV_ADV
#define UV_COR
#define MIX_S_UV
#define DJ_GRADPS
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#define TS_DIF2
#define MIX_GEO_TS
#define SOLVE3D
#define SPLINES
#define EW_PERIODIC
#define ANA_DIAG
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_BTFLUX

# elif defined SED_TEST1

/*
**  Suspended Sediment Test in a Channel.
*/

#define UV_ADV
#define UV_PSOURCE
#define UV_VIS4
#define MIX_S_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#define TS_DIF4
#define MIX_S_TS
#define SALINITY
#define SOLVE3D
#define SEDIMENT
#undef  SPLINES
#define AVERAGES
#ifdef AVERAGES
# define AVERAGES_AKV
# define AVERAGES_AKT
# define AVERAGES_AKS
#endif
#define NORTHERN_WALL
#define SOUTHERN_WALL
#define WEST_FSRADIATION
#define WEST_M2RADIATION
#define WEST_M3RADIATION
#define WEST_TGRADIENT
#undef  EAST_FSRADIATION
#define EAST_FSCLAMPED
#define EAST_M2RADIATION
#define EAST_M3RADIATION
#undef  EAST_TGRADIENT
#define EAST_TCLAMPED
#define MY25_MIXING
#ifdef MY25_MIXING
# define KANTHA_CLAYSON
# undef  N2S2_HORAVG
#endif
#undef  GLS_MIXING
#undef  LMD_MIXING
#ifdef  LMD_MIXING
# define LMD_RIMIX
# define LMD_SKPP
# define LMD_BKPP
#endif
#define ANA_BPFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SEDIMENT
#define ANA_SMFLUX
#define ANA_SPFLUX
#define ANA_SRFLUX
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_PSOURCE
#define ANA_TOBC
#define ANA_FSOBC
#undef  ANA_VMIX

# elif defined SOLITON

/*
**  Options for Equatorial Rossby Example:
*/

#define UV_ADV
#define UV_C4ADVECTION
#define UV_VIS2
#define UV_COR
#define AVERAGES
#define ANA_GRID
#define ANA_INITIAL
#define EW_PERIODIC
#define ANA_SMFLUX

# elif defined UPWELLING

/*
**  Options for Upwelling Example:
*/

#define UV_ADV
#define UV_COR
#undef  UV_VIS2
#undef  UV_VIS4
#undef  MIX_GEO_UV
#define TS_U3HADVECTION
#define TS_C4VADVECTION
#undef  WJ_GRADP
#define DJ_GRADPS
#undef  TS_DIF2
#undef  TS_DIF4
#undef  MIX_GEO_TS
#define SALINITY
#define EW_PERIODIC
#define SOLVE3D
#define SPLINES
#define AVERAGES
#define OUT_DOUBLE
#define ANA_GRID
#define ANA_INITIAL
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_BTFLUX
#define ANA_BSFLUX
#define ANA_VMIX
#undef  BIO_FASHAM
#ifdef BIO_FASHAM
# define ANA_BIOLOGY
# define ANA_SPFLUX
# define ANA_BPFLUX
# define ANA_SRFLUX
#endif

# elif defined USWEST

/*
**  Options for US West Coast Application.
*/

#define UV_ADV
#define UV_VIS2
#define UV_COR
#define DJ_GRADPS
#define NONLIN_EOS
#define SALINITY
#define AVERAGES
#define CURVGRID
#define MASKING
#define SOLVE3D
#define SPLINES
#define SPONGE
#define OBC_VOLCONS
#define M2CLIMATOL0GY
#define M3CLIMATOLOGY
#define TCLIMATOLOGY
#define M2CLM_NUDGING
#define M3CLM_NUDGING
#define TCLM_NUDGING
#define EASTERN_WALL
#define WEST_FSGRADIENT
#define WEST_M2RADIATION
#define WEST_M2NUDGING
#define WEST_M3RADIATION
#define WEST_M3NUDGING
#define WEST_TRADIATION
#define WEST_TNUDGING
#define SOUTH_FSGRADIENT
#define SOUTH_M2RADIATION
#define SOUTH_M2NUDGING
#define SOUTH_M3RADIATION
#define SOUTH_M3NUDGING
#define SOUTH_TRADIATION
#define SOUTH_TNUDGING
#define NORTH_FSGRADIENT
#define NORTH_M2RADIATION
#define NORTH_M2NUDGING
#define NORTH_M3RADIATION
#define NORTH_M3NUDGING
#define NORTH_TRADIATION
#define NORTH_TNUDGING
#define ANA_SSFLUX
#define ANA_STFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

# elif defined WEDDELL

/*
**  Options Idealized Weddell Sea Application: Tides and Ice Shelf Test.
*/

#define UV_ADV
#define DJ_GRADPS
#define UV_COR
#undef  UV_VIS4
#undef  MIX_S_UV
#define TS_A4HADVECTION
#define TS_A4VADVECTION
#undef  TS_DIF4
#undef  MIX_GEO_TS
#define SOLVE3D
#define SALINITY
#define NONLIN_EOS
#define CURVGRID
#define SPLINES
#define ICESHELF
#define AVERAGES
#define NS_PERIODIC
#define OBC_VOLCONS
#define RADIATION_2D
#define EAST_FSCHAPMAN
#define EAST_M2FLATHER
#define EAST_M3RADIATION
#define EAST_TRADIATION
#define WEST_FSCHAPMAN
#define WEST_M2FLATHER
#define WEST_M3RADIATION
#define WEST_TRADIATION
#define ZCLIMATOLOGY
#define M2CLIMATOLOGY
#define ANA_GRID
#define ANA_INITIAL
#define ANA_FSOBC
#define ANA_M2OBC
#define ANA_SMFLUX
#define ANA_STFLUX
#define ANA_SSFLUX
#define ANA_SRFLUX
#define ANA_BSFLUX
#define ANA_BTFLUX

# else

      CPPDEFS - Choose an appropriate ROMS application.

# endif

/*
**--------------------------------------------------------------------
**  Include other internal CPP definitions:
**--------------------------------------------------------------------
*/

#include "globaldefs.h"
