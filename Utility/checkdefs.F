#include "cppdefs.h"
      SUBROUTINE checkdefs
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This subroutine checks activated C-preprocessing options for        !
!  consistency.                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_strings
!
      implicit none
!
!  Local variable declarations.
!
      integer :: FS_east = 0
      integer :: FS_west = 0
      integer :: FS_north = 0
      integer :: FS_south = 0
      integer :: M2_east = 0
      integer :: M2_west = 0
      integer :: M2_north = 0
      integer :: M2_south = 0
      integer :: M3_east = 0
      integer :: M3_west = 0
      integer :: M3_north = 0
      integer :: M3_south = 0
      integer :: K_east = 0
      integer :: K_west = 0
      integer :: K_north = 0
      integer :: K_south = 0
      integer :: T_east = 0
      integer :: T_west = 0
      integer :: T_north = 0
      integer :: T_south = 0
      integer :: ibiology = 0
      integer :: iexample = 0
      integer :: ivmix = 0

      integer :: is, ng
!
!-----------------------------------------------------------------------
!  Report activated C-preprocessing options.
!-----------------------------------------------------------------------
!
      Coptions=' '
      IF (Master) write(stdout,10)
  10  format(/,' Activated C-preprocessing Options:',/)
  20  format(2x,a,t22,a)
#if defined ADD_FSOBC && defined SSH_TIDES
      IF (Master) write(stdout,20) 'ADD_FSOBC',                         &
     &   'Adding tidal elevation to proccesed OBC data'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ADD_FSOBC,'
#endif
#if defined ADD_M2OBC && defined UV_TIDES
      IF (Master) write(stdout,20) 'ADD_M2OBC',                         &
     &   'Adding tidal currents to proccesed OBC data'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ADD_M2OBC,'
#endif
#ifdef ADRIATIC1
      IF (Master) write(stdout,20) 'ADRIATIC1',                         &
     &   'Adriatic Sea Application, Grid 1.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ADRIATIC1,'
      iexample=iexample+1
#endif
#if defined ALBEDO && defined ANA_SRFLUX
      IF (Master) write(stdout,20) 'ALBEDO',                            &
     &   'Shortwave radiation from albedo equation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' ALBEDO,'
#endif
#if defined BIOLOGY && defined ANA_BIOLOGY
      IF (Master) write(stdout,20) 'ANA_BIOLOGY',                       &
     &   'Analytical biology initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_BIOLOGY,'
#endif
#if defined SEDIMENT || defined BIOLOGY
# ifdef ANA_BPFLUX
      IF (Master) write(stdout,20) 'ANA_BPFLUX',                        &
     &   'Analytical bottom passive tracers fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BPFLUX,'
# endif
#endif
#ifdef ANA_BSEDIM
      IF (Master) write(stdout,20) 'ANA_BSEDIM',                        &
     &   'Analytical bottom sediment size and density.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BSEDIM,'
#endif
#ifdef ANA_BSFLUX
      IF (Master) write(stdout,20) 'ANA_BSFLUX',                        &
     &   'Analytical kinematic bottom salinity flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BSFLUX,'
#endif
#ifdef ANA_BTFLUX
      IF (Master) write(stdout,20) 'ANA_BTFLUX',                        &
     &   'Analytical kinematic bottom temperature flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_BTFLUX,'
#endif
#if defined ANA_CLOUD && defined BULK_FLUXES
      IF (Master) write(stdout,20) 'ANA_CLOUD',                         &
     &   'Analytical cloud fraction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_CLOUD,'
#endif
#ifdef ANA_DIAG
      IF (Master) write(stdout,20) 'ANA_DIAG',                          &
     &   'Customized diagnostics.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_DIAG,'
#endif
#ifdef ANA_FSOBC
      IF (Master) write(stdout,20) 'ANA_FSOBC',                         &
     &   'Analytical free-surface boundary conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_FSOBC,'
#endif
#ifdef ANA_GRID
      IF (Master) write(stdout,20) 'ANA_GRID',                          &
     &   'Analytical grid set-up.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_GRID,'
#endif
#if defined ANA_HUMIDITY && defined BULK_FLUXES
      IF (Master) write(stdout,20) 'ANA_HUMIDITY',                      &
     &   'Analytical surface air humidity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ANA_HUMIDITY,'
#endif
#ifdef ANA_INITIAL
      IF (Master) write(stdout,20) 'ANA_INITIAL',                       &
     &   'Analytical initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_INITIAL,'
#endif
#ifdef ANA_M2CLIMA
      IF (Master) write(stdout,20) 'ANA_M2CLIMA',                       &
     &   'Analytical 2D momentum climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_M2CLIMA,'
#endif
#ifdef ANA_M2OBC
      IF (Master) write(stdout,20) 'ANA_M2OBC',                         &
     &   'Analytical 2D momentum boundary conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_M2OBC,'
#endif
#if defined ANA_M3CLIMA && defined SOLVE3D
      IF (Master) write(stdout,20) 'ANA_M3CLIMA',                       &
     &   'Analytical 3D momentum climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_M3CLIMA,'
#endif
#if defined ANA_M3OBC && defined SOLVE3D
      IF (Master) write(stdout,20) 'ANA_M3OBC',                         &
     &   'Analytical 3D momentum boundary conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_M3OBC,'
#endif
#if defined ANA_PAIR && defined BULK_FLUXES
      IF (Master) write(stdout,20) 'ANA_PAIR',                          &
     &   'Analytical surface air pressure.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_PAIR,'
#endif
#ifdef ANA_PSOURCE
      IF (Master) write(stdout,20) 'ANA_PSOURCE',                       &
     &   'Analytical point sources and sinks.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' ANA_PSOURCE,'
#endif
#if defined ANA_RAIN && defined BULK_FLUXES
      IF (Master) write(stdout,20) 'ANA_RAIN',                          &
     &   'Analytical rain fall rate.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_RAIN,'
#endif
#if defined SEDIMENT && defined ANA_SEDIMENT
      IF (Master) write(stdout,20) 'ANA_SEDIMENT',                      &
     &   'Analytical sediment initial conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ANA_SEDIMENT,'
#endif
#ifdef ANA_SMFLUX
      IF (Master) write(stdout,20) 'ANA_SMFLUX',                        &
     &   'Analytical kinematic surface momentum flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SMFLUX,'
#endif
#if defined SEDIMENT || defined BIOLOGY
# ifdef ANA_SPFLUX
      IF (Master) write(stdout,20) 'ANA_SPFLUX',                        &
     &   'Analytical surface passive tracer fluxes.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SPFLUX,'
# endif
#endif
#ifdef ANA_SPINNING
      IF (Master) write(stdout,20) 'ANA_SPINNING',                      &
     &   'Analytical time-varying rotation force.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ANA_SPINNING,'
#endif
#ifdef ANA_SRFLUX
      IF (Master) write(stdout,20) 'ANA_SRFLUX',                        &
     &   'Analytical kinematic shortwave radiation flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SRFLUX,'
#endif
#ifdef ANA_SSH
      IF (Master) write(stdout,20) 'ANA_SSH',                           &
     &   'Analytical sea surface height climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ANA_SSH,'
#endif
#ifdef ANA_SSS
      IF (Master) write(stdout,20) 'ANA_SSS',                           &
     &   'Analytical sea surface salinity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ANA_SSS,'
#endif
#ifdef ANA_SST
      IF (Master) write(stdout,20) 'ANA_SST',                           &
     &   'Analytical SST and dQdSST.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' ANA_SST,'
#endif
#ifdef ANA_SSFLUX
      IF (Master) write(stdout,20) 'ANA_SSFLUX',                        &
     &   'Analytical kinematic surface salinity flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_SSFLUX,'
#endif
#ifdef ANA_STFLUX
      IF (Master) write(stdout,20) 'ANA_STFLUX',                        &
     &   'Analytical kinematic surface temperature flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_STFLUX,'
#endif
#if defined ANA_TCLIMA && defined SOLVE3D
      IF (Master) write(stdout,20) 'ANA_TCLIMA',                        &
     &   'Analytical tracer climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' ANA_TCLIMA,'
#endif
#if defined ANA_TOBC && defined SOLVE3D
      IF (Master) write(stdout,20) 'ANA_TOBC',                          &
     &   'Analytical tracers boundary conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_TOBC,'
#endif
#if defined ANA_VMIX && defined SOLVE3D
      IF (Master) write(stdout,20) 'ANA_VMIX',                          &
     &   'Analytical vertical mixing coefficients.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ANA_VMIX,'
      ivmix=ivmix+1
#endif
#if defined ANA_WINDS && defined BULK_FLUXES
      IF (Master) write(stdout,20) 'ANA_WINDS',                         &
     &   'Analytical surface wind components.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_WINDS,'
#endif
#ifdef ANA_WWAVE
      IF (Master) write(stdout,20) 'ANA_WWAVE',                         &
     &   'Analytical wind induced waves.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' ANA_WWAVE,'
#endif
#ifdef ASSIMILATION_SSH
      IF (Master) write(stdout,20) 'ASSIMILATION_SSH',                  &
     &   'Assimilating SSH observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' ASSIMILATION_SSH,'
#endif
#ifdef ASSIMILATION_SST
      IF (Master) write(stdout,20) 'ASSIMILATION_SST',                  &
     &   'Assimilating SST observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' ASSIMILATION_SST,'
#endif
#ifdef ASSIMILATION_T
      IF (Master) write(stdout,20) 'ASSIMILATION_T',                    &
     &   'Assimilating tracers observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' ASSIMILATION_T,'
#endif
#ifdef ASSIMILATION_UV
      IF (Master) write(stdout,20) 'ASSIMILATION_UV',                   &
     &   'Assimilating horizontal currents observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' ASSIMILATION_UV,'
#endif
#ifdef ASSIMILATION_UVsur
      IF (Master) write(stdout,20) 'ASSIMILATION_UVsur',                &
     &   'Assimilating surface currents observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+20)=' ASSIMILATION_UVsur,'
#endif
#ifdef ASSUMED_SHAPE
      IF (Master) write(stdout,20) 'ASSUMED_SHAPE',                     &
     &   'Using assumed-shape arrays.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' ASSUMED_SHAPE,'
#endif
#ifdef AVERAGES
      IF (Master) write(stdout,20) 'AVERAGES',                          &
     &   'Writing out time-averaged fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' AVERAGES,'
# ifdef AVERAGES_AKS
      IF (Master) write(stdout,20) 'AVERAGES_AKS',                      &
     &   'Writing out time-averaged vertical S-diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' AVERAGES_AKS,'
# endif
# ifdef AVERAGES_AKT
      IF (Master) write(stdout,20) 'AVERAGES_AKT',                      &
     &   'Writing out time-averaged vertical T-diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' AVERAGES_AKT,'
# endif
# ifdef AVERAGES_AKV
      IF (Master) write(stdout,20) 'AVERAGES_AKV',                      &
     &   'Writing out time-averaged vertical viscosity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' AVERAGES_AKV,'
# endif
#endif
#ifdef BASIN
      IF (Master) write(stdout,20) 'BASIN',                             &
     &   'Big Bad Basin Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+7)=' BASIN,'
      iexample=iexample+1
#endif
#ifdef BENCHMARK1
      IF (Master) write(stdout,20) 'BENCHMARK1',                        &
     &   'Small Benchmark, Idealized Southern Ocean.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BENCHMARK1,'
      iexample=iexample+1
#endif
#ifdef BENCHMARK2
      IF (Master) write(stdout,20) 'BENCHMARK2',                        &
     &   'Medium Benchmark, Idealized Southern Ocean.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BENCHMARK2,'
      iexample=iexample+1
#endif
#ifdef BENCHMARK3
      IF (Master) write(stdout,20) 'BENCHMARK3',                        &
     &   'Big Benchmark, Idealized Southern Ocean.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BENCHMARK3,'
      iexample=iexample+1
#endif
#ifdef BIO_FASHAM
      IF (Master) write(stdout,20) 'BIO_FASHAM',                        &
     &   'Fasham type nitrogen-based model.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BIO_FASHAM,'
      ibiology=ibiology+1
#endif
#ifdef BISCAY
      IF (Master) write(stdout,20) 'BISCAY',                            &
     &   'Bay of Biscay Application.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' BISCAY,'
      iexample=iexample+1
#endif
#ifdef BL_TEST
      IF (Master) write(stdout,20) 'BL_TEST',                           &
     &   'Boundary Layers Test.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' BL_TEST,'
      iexample=iexample+1
#endif
#ifdef BODYFORCE
      IF (Master) write(stdout,20) 'BODYFORCE',                         &
     &   'Momentum stresses as body-forces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' BODYFORCE,'
#endif
#ifdef BULK_FLUXES
      IF (Master) write(stdout,20) 'BULK_FLUXES',                       &
     &   'Surface bulk fluxes parametererization.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' BULK_FLUXES,'
#endif
#ifdef BVF_MIXING
      IF (Master) write(stdout,20) 'BVF_MIXING',                        &
     &   'Brunt-Vaisala frequency based vertical mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' BVF_MIXING,'
      ivmix=ivmix+1
#endif
#ifdef CALCOFI
      IF (Master) write(stdout,20) 'CALCOFI',                           &
     &   'Southern California Bight Application.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' CALCOFI,'
      iexample=iexample+1
#endif
#if defined CANUTO_A && defined GLS_MIXING
      IF (Master) write(stdout,20) 'CANUTO_A',                          &
     &   'Canuto A-stability function formulation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' CANUTO_A,'
#endif
#if defined CANUTO_B && defined GLS_MIXING
      IF (Master) write(stdout,20) 'CANUTO_B',                          &
     &   'Canuto B-stability function formulation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' CANUTO_B,'
#endif
#ifdef CANYON_A
      IF (Master) write(stdout,20) 'CANYON_A',                          &
     &   'Canyon A Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' CANYON_A,'
      iexample=iexample+1
#endif
#ifdef CANYON_B
      IF (Master) write(stdout,20) 'CANYON_B',                          &
     &   'Canyon B Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' CANYON_B,'
      iexample=iexample+1
#endif
#if defined COOL_SKIN && defined BULK_FLUXES
      IF (Master) write(stdout,20) 'COOL_SKIN',                         &
     &   'Surface cool skin correction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' COOL_SKIN,'
#endif
#if defined COSINE2 && defined SOLVE3D
      IF (Master) write(stdout,20) 'COSINE2',                           &
     &   'Cosine-squared shape time-averaging barotropic filter.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' COSINE2,'
#endif
#ifdef CS_BBL
      IF (Master) write(stdout,20) 'CS_BBL',                            &
     &   'Sherwood Bottom Boundary Layer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' CS_BBL,'
#endif
#ifdef UV_ADV
# ifdef CURVGRID
      IF (Master) write(stdout,20) 'CURVGRID',                          &
     &   'Orthogonal curvilinear grid.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' CURVGRID,'
# endif
#endif
#ifdef DAMEE_4
      IF (Master) write(stdout,20) 'DAMEE_4',                           &
     &   'North Atlantic DAMEE Application, Grid #4.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' DAMEE_4,'
      iexample=iexample+1
#endif
#ifdef DAMEE_5
      IF (Master) write(stdout,20) 'DAMEE_5',                           &
     &   'North Atlantic DAMEE Application, Grid #5.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' DAMEE_5,'
      iexample=iexample+1
#endif
#ifdef DAMEE_6
      IF (Master) write(stdout,20) 'DAMEE_6',                           &
     &   'North Atlantic DAMEE Application, Grid #6.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' DAMEE_6,'
      iexample=iexample+1
#endif
#ifdef DAMEE_7
      IF (Master) write(stdout,20) 'DAMEE_7',                           &
     &   'North Atlantic DAMEE Application, Grid #7.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' DAMEE_7,'
      iexample=iexample+1
#endif
#ifdef DAMEE_9
      IF (Master) write(stdout,20) 'DAMEE_9',                           &
     &   'North Atlantic DAMEE Application, Grid #9.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' DAMEE_9,'
      iexample=iexample+1
#endif
#ifdef DEBUGGING
      IF (Master) write(stdout,20) 'DEBUGGING',                         &
     &   'Internal debugging switch activated.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DEBUGGING,'
#endif
#ifdef DIAGNOSTIC
      IF (Master) write(stdout,20) 'DIAGNOSTIC',                        &
     &   'Diagnostic configuration, no evolution of tracer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' DIAGNOSTIC,'
#endif
#if defined TS_DIF2 || defined TS_DIF4
# ifdef DIFF_GRID
      IF (Master) write(stdout,20) 'DIFF_GRID',                         &
     &   'Horizontal diffusion coefficient scaled by grid size.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DIFF_GRID,'
# endif
#endif
#ifdef DIURNAL_SRFLUX
      IF (Master) write(stdout,20) 'DIURNAL_SRFLUX',                    &
     &   'Modulate shortwave radiation by the local diurnal cycle.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' DIURNAL_SRFLUX,'
#endif
#ifdef DJ_GRADPS
      IF (Master) write(stdout,20) 'DJ_GRADPS',                         &
     &   'Parabolic Splines density Jacobian (Shchepetkin, 2002).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' DJ_GRADPS,'
#endif
#ifdef DOUBLE_GYRE
      IF (Master) write(stdout,20) 'DOUBLE_GYRE',                       &
     &   'Idealized Double Gyre Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' DOUBLE_GYRE,'
      iexample=iexample+1
#endif
#ifdef DOUBLE_PRECISION
      IF (Master) write(stdout,20) 'DOUBLE_PRECISION',                  &
     &   'Double precision arithmetic.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' DOUBLE_PRECISION,'
#endif
#ifdef EAST_FSCLAMPED
      IF (Master) write(stdout,20) 'EAST_FSCLAMPED',                    &
     &   'Eastern edge, free-surface, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_FSCLAMPED,'
      FS_east=FS_east+1
#endif
#ifdef EAST_FSCHAPMAN
      IF (Master) write(stdout,20) 'EAST_FSCHAPMAN',                    &
     &   'Eastern edge, free-surface, Chapman condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_FSCHAPMAN,'
      FS_east=FS_east+1
#endif
#ifdef EAST_FSGRADIENT
      IF (Master) write(stdout,20) 'EAST_FSGRADIENT',                   &
     &   'Eastern edge, free-surface, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_FSGRADIENT,'
      FS_east=FS_east+1
#endif
#ifdef EAST_FSRADIATION
# ifdef EAST_FSNUDGING
      IF (Master) write(stdout,20) 'EAST_FSNUDGING',                    &
     &   'Eastern edge, free-surface, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_FSNUDGING,'
# endif
      IF (Master) write(stdout,20) 'EAST_FSRADIATION',                  &
     &   'Eastern edge, free-surface, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' EAST_FSRADIATION,'
      FS_east=FS_east+1
#endif
#ifdef EAST_M2CLAMPED
      IF (Master) write(stdout,20) 'EAST_M2CLAMPED',                    &
     &   'Eastern edge, 2D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M2CLAMPED,'
      M2_east=M2_east+1
#endif
#ifdef EAST_M2FLATHER
      IF (Master) write(stdout,20) 'EAST_M2FLATER',                     &
     &   'Eastern edge, 2D momentum, Flather condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M2FLATHER,'
      M2_east=M2_east+1
#endif
#ifdef EAST_M2GRADIENT
      IF (Master) write(stdout,20) 'EAST_M2GRADIENT',                   &
     &   'Eastern edge, 2D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_M2GRADIENT,'
      M2_east=M2_east+1
#endif
#ifdef EAST_M2RADIATION
# ifdef EAST_M2NUDGING
      IF (Master) write(stdout,20) 'EAST_M2NUDGING',                    &
     &   'Eastern edge, 2D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M2NUDGING,'
# endif
      IF (Master) write(stdout,20) 'EAST_M2RADIATION',                  &
     &   'Eastern edge, 2D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' EAST_M2RADIATION,'
      M2_east=M2_east+1
#endif
#ifdef EAST_M3CLAMPED
      IF (Master) write(stdout,20) 'EAST_M3CLAMPED',                    &
     &   'Eastern edge, 3D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M3CLAMPED,'
      M3_east=M3_east+1
#endif
#ifdef EAST_M3GRADIENT
      IF (Master) write(stdout,20) 'EAST_M3GRADIENT',                   &
     &   'Eastern edge, 3D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_M3GRADIENT,'
      M3_east=M3_east+1
#endif
#ifdef EAST_M3RADIATION
# ifdef EAST_M3NUDGING
      IF (Master) write(stdout,20) 'EAST_M3NUDGING',                    &
     &   'Eastern edge, 3D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_M3NUDGING,'
# endif
      IF (Master) write(stdout,20) 'EAST_M3RADIATION',                  &
     &   'Eastern edge, 3D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' EAST_M3RADIATION,'
      M3_east=M3_east+1
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# ifdef EAST_KGRADIENT
      IF (Master) write(stdout,20) 'EAST_KGRADIENT',                    &
     &   'Eastern edge, TKE fields, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_KGRADIENT,'
      K_east=K_east+1
# endif
# ifdef EAST_KRADIATION
      IF (Master) write(stdout,20) 'EAST_KRADIATION',                   &
     &   'Eastern edge, TKE fields, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_KRADIATION,'
      K_east=K_east+1
# endif
#endif
#ifdef EAST_TCLAMPED
      IF (Master) write(stdout,20) 'EAST_TCLAMPED',                     &
     &   'Eastern edge, tracers, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' EAST_TCLAMPED,'
      T_east=T_east+1
#endif
#ifdef EAST_TGRADIENT
      IF (Master) write(stdout,20) 'EAST_TGRADIENT',                    &
     &   'Eastern edge, tracers, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' EAST_TGRADIENT,'
      T_east=T_east+1
#endif
#ifdef EAST_TRADIATION
# ifdef EAST_TNUDGING
      IF (Master) write(stdout,20) 'EAST_TNUDGING',                     &
     &   'Eastern edge, tracers, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' EAST_TNUDGING,'
# endif
      IF (Master) write(stdout,20) 'EAST_TRADIATION',                   &
     &   'Eastern edge, tracers, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' EAST_TRADIATION,'
      T_east=T_east+1
#endif
#ifdef EASTERN_WALL
      IF (Master) write(stdout,20) 'EASTERN_WALL',                      &
     &   'Wall boundary at Eastern edge.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' EASTERN_WALL,'
#endif
#ifdef EW_PERIODIC
      IF (Master) write(stdout,20) 'EW_PERIODIC',                       &
     &   'East-West periodic boundaries.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' EW_PERIODIC,'
      FS_east=FS_east+1
      FS_west=FS_west+1
      M2_east=M2_east+1
      M2_west=M2_west+1
      M3_east=M3_east+1
      M3_west=M3_west+1
      K_east=K_east+1
      K_west=K_west+1
      T_east=T_east+1
      T_west=T_west+1
#endif
#ifdef FLOATS
      IF (Master) write(stdout,20) 'FLOATS',                            &
     &   'Simulated Lagrangian drifters.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' FLOATS,'
#endif
#if defined FLOATS && defined FLOAT_VWALK && defined SOLVE3D
      IF (Master) write(stdout,20) 'FLOAT_VWALK',                       &
     &   'Lagrangian drifter vertical diffusion.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' FLOAT_VWALK,'
#endif
#ifdef FLT_TEST
      IF (Master) write(stdout,20) 'FLT_TEST',                          &
     &   'Float Tracking Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' FLT_TEST,'
      iexample=iexample+1
#endif
#ifdef GLS_MIXING
      IF (Master) write(stdout,20) 'GLS_MIXING',                        &
     &   'Generic Length-Scale turbulence closure.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' GLS_MIXING,'
      ivmix=ivmix+1
#endif
#ifdef GRAV_ADJ
      IF (Master) write(stdout,20) 'GRAV_ADJ',                          &
     &   'Gravitational Adjustment Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' GRAV_ADJ,'
      iexample=iexample+1
#endif
#ifdef ICESHELF
      IF (Master) write(stdout,20) 'ICESHELF',                          &
     &   'Include Ice Shelf Cavities.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' ICESHELF,'
#endif
#if defined KANTHA_CLAYSON && (defined GLS_MIXING || defined MY25_MIXING)
      IF (Master) write(stdout,20) 'KANTHA_CLAYSON',                    &
     &   'Kantha and Clayson stability function formulation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' KANTHA_CLAYSON,'
#endif
#ifdef KELVIN
      IF (Master) write(stdout,20) 'KELVIN',                            &
     &   'Kelvin wave test.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' KELVIN,'
      iexample=iexample+1
#endif
#ifdef LAB_CANYON
      IF (Master) write(stdout,20) 'LAB_CANYON',                        &
     &   'Lab Canyon, an annulus example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' LAB_CANYON,'
      iexample=iexample+1
#endif
#ifdef LMD_BKPP
      IF (Master) write(stdout,20) 'LMD_BKPP',                          &
     &   'KPP bottom boundary layer mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' LMD_BKPP,'
#endif
#ifdef LMD_CONVEC
      IF (Master) write(stdout,20) 'LMD_CONVEC',                        &
     &   'LMD convective mixing due to shear instability.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' LMD_CONVEC,'
#endif
#ifdef LMD_DDMIX
      IF (Master) write(stdout,20) 'LMD_DDMIX',                         &
     &   'LMD double-diffusive mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' LMD_DDMIX,'
#endif
#ifdef LMD_MIXING
      IF (Master) write(stdout,20) 'LMD_MIXING',                        &
     &   'Large/McWilliams/Doney interior mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' LMD_MIXING,'
      ivmix=ivmix+1
#endif
#ifdef LMD_NONLOCAL
      IF (Master) write(stdout,20) 'LMD_NONLOCAL',                      &
     &   'LMD convective nonlocal transport.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' LMD_NONLOCAL,'
#endif
#ifdef LMD_RIMIX
      IF (Master) write(stdout,20) 'LMD_RIMIX',                         &
     &   'LMD diffusivity due to shear instability.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' LMD_RIMIX,'
#endif
#ifdef LMD_SKPP
      IF (Master) write(stdout,20) 'LMD_SKPP',                          &
     &   'KPP surface boundary layer mixing.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' LMD_SKPP,'
#endif
#ifdef LMD_TEST
      IF (Master) write(stdout,20) 'LMD_TEST',                          &
     &                 'Large et al. Inertial Response Test.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' LMD_TEST,'
      iexample=iexample+1
#endif
#if defined LONGWAVE && defined BULK_FLUXES
      IF (Master) write(stdout,20) 'LONGWAVE',                          &
     &   'Compute net longwave radiation internally.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' LONGWAVE,'
#endif
#ifdef M2CLIMATOLOGY
      IF (Master) write(stdout,20) 'M2CLIMATOLOGY',                     &
     &   'Processing 2D momentum climatology data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' M2CLIMATOLOGY,'
#endif
#ifdef M2CLM_NUDGING
      IF (Master) write(stdout,20) 'M2CLM_NUDGING',                     &
     &   'Nudging toward 2D momentum climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' M2CLM_NUDGING,'
#endif
#if defined M3CLIMATOLOGY && defined SOLVE3D
      IF (Master) write(stdout,20) 'M3CLIMATOLOGY',                     &
     &   'Processing 3D momentum climatology data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' M3CLIMATOLOGY,'
#endif
#if defined M3CLM_NUDGING && defined SOLVE3D
      IF (Master) write(stdout,20) 'M3CLM_NUDGING',                     &
     &   'Nudging toward 3D momentum climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' M3CLM_NUDGING,'
#endif
#ifdef MASKING
      IF (Master) write(stdout,20) 'MASKING',                           &
     &   'Land/Sea masking.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' MASKING,'
#endif
#ifdef MB_BBL
      IF (Master) write(stdout,20) 'MB_BBL',                            &
     &   'Blaas Bottom Boundary Layer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' MB_BBL,'
#endif
#if (defined TS_DIF2 || defined TS_DIF4) && defined SOLVE3D
# ifdef MIX_GEO_TS
      IF (Master) write(stdout,20) 'MIX_GEO_TS',                        &
     &   'Mixing of tracers along geopotential surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' MIX_GEO_TS,'
# endif
# ifdef MIX_ISO_TS
      IF (Master) write(stdout,20) 'MIX_ISO_TS',                        &
     &   'Mixing of tracers along isopycnal surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' MIX_ISO_TS,'
# endif
# ifdef MIX_S_TS
      IF (Master) write(stdout,20) 'MIX_S_TS',                          &
     &   'Mixing of tracers along constant S-surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' MIX_S_TS,'
# endif
#endif
#if (defined UV_VIS2 || defined UV_VIS4) && defined SOLVE3D
# ifdef MIX_GEO_UV
      IF (Master) write(stdout,20) 'MIX_GEO_UV',                        &
     &   'Mixing of momentum along geopotential surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' MIX_GEO_UV,'
# endif
# ifdef MIX_S_UV
      IF (Master) write(stdout,20) 'MIX_S_UV',                          &
     &   'Mixing of momentum along constant S-surfaces.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' MIX_S_UV,'
# endif
#endif
#ifdef MPI
      IF (Master) write(stdout,20) 'MPI',                               &
     &   'MPI distributed-memory configuration.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+4)=' MPI,'
#endif
#ifdef MY25_MIXING
      IF (Master) write(stdout,20) 'MY25_MIXING',                       &
     &   'Mellor/Yamada Level-2.5 mixing closure.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' MY25_MIXING,'
      ivmix=ivmix+1
#endif
#ifdef NATL
      IF (Master) write(stdout,20) 'NATL',                              &
     &   'High Resolution North Atlantic DAMEE Application.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+6)=' NATL,'
      iexample=iexample+1
#endif
#ifdef NENA
      IF (Master) write(stdout,20) 'NENA',                              &
     &   'North East North America Application.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+6)=' NENA,'
      iexample=iexample+1
#endif
#ifdef NJ_BIGHT
      IF (Master) write(stdout,20) 'NJ_BIGHT',                          &
     &   'New Jersey Bight Application.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' NJ_BIGHT,'
      iexample=iexample+1
#endif
#if defined NONLIN_EOS && defined SOLVE3D
      IF (Master) write(stdout,20) 'NONLIN_EOS',                        &
     &   'Non-linear Equation of State for seawater.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' NONLIN_EOS,'
#endif
#ifdef NORTHERN_WALL
      IF (Master) write(stdout,20) 'NORTHERN_WALL',                     &
     &   'Wall boundary at Northern edge.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' NORTHERN_WALL,'
#endif
#ifdef NPACIFIC
      IF (Master) write(stdout,20) 'PACIFIC',                           &
     &   'North Pacific Application.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' NPACIFIC,'
      iexample=iexample+1
#endif
#ifdef NORTH_FSCLAMPED
      IF (Master) write(stdout,20) 'NORTH_FSCLAMPED',                   &
     &   'Northern edge, free-surface, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_FSCLAMPED,'
      FS_north=FS_north+1
#endif
#ifdef NORTH_FSCHAPMAN
      IF (Master) write(stdout,20) 'NORTH_FSCHAPMAN',                   &
     &   'Northern edge, free-surface, Chapman condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_FSCHAPMAN,'
      FS_north=FS_north+1
#endif
#ifdef NORTH_FSGRADIENT
      IF (Master) write(stdout,20) 'NORTH_FSGRADIENT',                  &
     &   'Northern edge, free-surface, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NORTH_FSGRADIENT,'
      FS_north=FS_north+1
#endif
#ifdef NORTH_FSRADIATION
# ifdef NORTH_FSNUDGING
      IF (Master) write(stdout,20) 'NORTH_FSNUDGING',                   &
     &   'Northern edge, free-surface, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_FSNUDGING,'
# endif
      IF (Master) write(stdout,20) 'NORTH_FSRADIATION',                 &
     &   'Northern edge, free-surface, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' NORTH_FSRADIATION,'
      FS_north=FS_north+1
#endif
#ifdef NORTH_M2CLAMPED
      IF (Master) write(stdout,20) 'NORTH_M2CLAMPED',                   &
     &   'Northern edge, 2D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M2CLAMPED,'
      M2_north=M2_north+1
#endif
#ifdef NORTH_M2FLATHER
      IF (Master) write(stdout,20) 'NORTH_M2FLATHER',                   &
     &   'Northern edge, 2D momentum, Flather condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M2FLATHER,'
      M2_north=M2_north+1
#endif
#ifdef NORTH_M2GRADIENT
      IF (Master) write(stdout,20) 'NORTH_M2GRADIENT',                  &
     &   'Northern edge, 2D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NORTH_M2GRADIENT,'
      M2_north=M2_north+1
#endif
#ifdef NORTH_M2RADIATION
# ifdef NORTH_M2NUDGING
      IF (Master) write(stdout,20) 'NORTH_M2NUDGING',                   &
     &   'Northern edge, 2D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M2NUDGING,'
# endif
      IF (Master) write(stdout,20) 'NORTH_M2RADIATION',                 &
     &   'Northern edge, 2D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' NORTH_M2RADIATION,'
      M2_north=M2_north+1
#endif
#ifdef NORTH_M3CLAMPED
      IF (Master) write(stdout,20) 'NORTH_M3CLAMPED',                   &
     &   'Northern edge, 3D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M3CLAMPED,'
      M3_north=M3_north+1
#endif
#ifdef NORTH_M3GRADIENT
      IF (Master) write(stdout,20) 'NORTH_M3GRADIENT',                  &
     &   'Northern edge, 3D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NORTH_M3GRADIENT,'
      M3_north=M3_north+1
#endif
#ifdef NORTH_M3RADIATION
# ifdef NORTH_M3NUDGING
      IF (Master) write(stdout,20) 'NORTH_M3NUDGING',                   &
     &   'Northern edge, 3D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_M3NUDGING,'
# endif
      IF (Master) write(stdout,20) 'NORTH_M3RADIATION',                 &
     &   'Northern edge, 3D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' NORTH_M3RADIATION,'
      M3_north=M3_north+1
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# ifdef NORTH_KGRADIENT
      IF (Master) write(stdout,20) 'NORTH_KGRADIENT',                   &
     &   'Northern edge, TKE fields, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_KGRADIENT,'
      K_north=K_north+1
# endif
# ifdef NORTH_KRADIATION
      IF (Master) write(stdout,20) 'NORTH_KRADIATION',                  &
     &                 'Northern edge, TKE fields, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_KRADIATION,'
      K_north=K_north+1
# endif
#endif
#ifdef NORTH_TCLAMPED
      IF (Master) write(stdout,20) 'NORTH_TCLAMPED',                    &
     &   'Northern edge, tracers, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' NORTH_TCLAMPED,'
      T_north=T_north+1
#endif
#ifdef NORTH_TGRADIENT
      IF (Master) write(stdout,20) 'NORTH_TGRADIENT',                   &
     &   'Northern edge, tracers, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' NORTH_TGRADIENT,'
      T_north=T_north+1
#endif
#ifdef NORTH_TRADIATION
# ifdef NORTH_TNUDGING
      IF (Master) write(stdout,20) 'NORTH_TNUDGING',                    &
     &   'Northern edge, tracers, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' NORTH_TNUDGING,'
# endif
      IF (Master) write(stdout,20) 'NORTH_TRADIATION',                  &
     &   'Northern edge, tracers, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' NORTH_TRADIATION,'
      T_north=T_north+1
#endif
#ifdef NO_WRITE_GRID
      IF (Master) write(stdout,20) 'NO_WRITE_GRID',                     &
     &  'Not Writing grid arrays into NetCDF ouput files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' NO_WRITE_GRID,'
#endif
#ifdef NUDGING_SSH
      IF (Master) write(stdout,20) 'NUDGING_SSH',                       &
     &   'Nudging SSH observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' NUDGING_SSH,'
#endif
#ifdef NUDGING_SST
      IF (Master) write(stdout,20) 'NUDGING_SST',                       &
     &   'Nudging SST observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' NUDGING_SST,'
#endif
#ifdef NUDGING_T
      IF (Master) write(stdout,20) 'NUDGING_T',                         &
     &   'Nudging tracers observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' NUDGING_T,'
#endif
#ifdef NUDGING_UV
      IF (Master) write(stdout,20) 'NUDGING_UV',                        &
     &   'Nudging horizontal currents observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' NUDGING_UV,'
#endif
#ifdef NUDGING_UVsur
      IF (Master) write(stdout,20) 'NUDGING_UVsur',                     &
     &   'Nudging surface currents observations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' NUDGING_UVsur,'
#endif
#ifdef NS_PERIODIC
      IF (Master) write(stdout,20) 'NS_PERIODIC',                       &
     &   'North-South periodic boundaries.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' NS_PERIODIC,'
      FS_north=FS_north+1
      FS_south=FS_south+1
      M2_north=M2_north+1
      M2_south=M2_south+1
      M3_north=M3_north+1
      M3_south=M3_south+1
      K_north=K_north+1
      K_south=K_south+1
      T_north=T_north+1
      T_south=T_south+1
#endif
#if defined N2S2_HORAVG && (defined GLS_MIXING || defined MY25_MIXING)
      IF (Master) write(stdout,20) 'N2S2_HORAVG',                       &
     &   'Horizontal smoothing of buoyancy and shear.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' N2S2_HORAVG,'
#endif
#ifdef OUT_DOUBLE
      IF (Master) write(stdout,20) 'OUT_DOUBLE',                        &
     &   'Double precision output fields in NetCDF files.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' OUT_DOUBLE,'
#endif
#ifdef OBC_VOLCONS
      IF (Master) write(stdout,20) 'OBC_VOLCONS',                       &
     &   'Enforce mass conservation on open boundaries.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' OBC_VOLCONS,'
#endif
#ifdef _OPENMP
      IF (Master) write(stdout,20) '_OPENMP',                           &
     &   'OpenMP parallel shared-memory directives.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' _OPENMP,'
#endif
#ifdef OVERFLOW
      IF (Master) write(stdout,20) 'OVERFLOW',                          &
     &   'gravitational/overflow Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' OVERFLOW,'
      iexample=iexample+1
#endif
#ifdef PJ_GRADP
      IF (Master) write(stdout,20) 'PJ_GRADP',                          &
     &   'Finite volume Pressure Jacobian (Lin, 1997).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' PJ_GRADP,'
#endif
#ifdef PJ_GRADPQ2
      IF (Master) write(stdout,20) 'PJ_GRADPQ2',                        &
     &   'Quartic 2 polynomial Pressure Jacobian (Shchepetkin, 2002).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' PJ_GRADPQ2,'
#endif
#ifdef PJ_GRADPQ4
      IF (Master) write(stdout,20) 'PJ_GRADPQ4',                        &
     &   'Quartic 4 polynomial Pressure Jacobian (Shchepetkin, 2002).'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' PJ_GRADPQ4,'
#endif
#if defined POWER_LAW && defined SOLVE3D
      IF (Master) write(stdout,20) 'POWER_LAW',                         &
     &   'Power-law shape time-averaging barotropic filter.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' POWER_LAW,'
#endif
#ifdef PROFILE
      IF (Master) write(stdout,20) 'PROFILE',                           &
     &   'Time profiling activated .'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' PROFILE,'
#endif
#ifdef QCORRECTION
      IF (Master) write(stdout,20) 'QCORRECTION',                       &
     &   'Surface net heat flux correction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' QCORRECTION,'
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# if defined K_C2ADVECTION
      IF (Master) write(stdout,20) 'K_C2ADVECTION',                     &
     &   'Second-order centered differences advection of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' K_C2ADVECTION,'
# elif defined K_C4ADVECTION
      IF (Master) write(stdout,20) 'K_C4ADVECTION',                     &
     &   'Fourth-order centered differences advection of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' K_C4ADVECTION,'
# else
      IF (Master) write(stdout,20) 'K_GSCHEME',                         &
     &   'Third-order upstream bias advection of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' K_GSCHEME,'
# endif
# ifdef TKE_DIF2
      IF (Master) write(stdout,20) 'TKE_DIF2',                          &
     &   'Harmonic mixing of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' TKE_DIF2,'
# endif
# ifdef TKE_DIF4
      IF (Master) write(stdout,20) 'TKE_DIF4',                          &
     &   'Biharmonic mixing of TKE fields.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' TKE_DIF4,'
# endif
#endif
#ifdef RADIATION_2D
      IF (Master) write(stdout,20) 'RADIATION_2D',                      &
     &   'Use tangential phase speed in radiation conditions.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' RADIATION_2D,'
#endif
#ifdef RANDOM_ESPERT
      IF (Master) write(stdout,20) 'RANDOM_ESPERT',                     &
     &   'Random Error Subspace perturbations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' RANDOM_ESPERT,'
#endif
#if defined READ_WATER && defined MASKING
      IF (Master) write(stdout,20) 'READ_WATER',                        &
     &   'Reading data at water points only.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' READ_WATER,'
#endif
#if defined LMD_MIXING
# ifdef RI_HORAVG
      IF (Master) write(stdout,20) 'RI_HORAVG',                         &
     &   'Smooth Richardson number horizontally.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' RI_HORAVG,'
# endif
# ifdef RI_VERAVG
      IF (Master) write(stdout,20) 'RI_VERAVG',                         &
     &   'Smooth Richardson number vertically.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' RI_VERAVG,'
# endif
#endif
#ifdef RIVERPLUME
      IF (Master) write(stdout,20) 'RIVERPLUME',                        &
     &   'River Plume Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' RIVERPLUME,'
      iexample=iexample+1
#endif
#ifdef RST_SINGLE
      IF (Master) write(stdout,20) 'RST_SINGLE',                        &
     &   'Single precision fields in restart NetCDF file.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' RST_SINGLE,'
#else
      IF (Master) write(stdout,20) '!RST_SINGLE',                       &
     &   'Double precision fields in restart NetCDF file.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' !RST_SINGLE,'
#endif
#ifdef SALINITY
      IF (Master) write(stdout,20) 'SALINITY',                          &
     &   'Using salinity.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SALINITY,'
#endif
#ifdef SCORRECTION
      IF (Master) write(stdout,20) 'SCORRECTION',                       &
     &   'Surface salinity flux correction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' SCORRECTION,'
#endif
#ifdef SEAMOUNT
      IF (Master) write(stdout,20) 'SEAMOUNT',                          &
     &   'Seamount Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SEAMOUNT,'
      iexample=iexample+1
#endif
#ifdef SEDIMENT
      IF (Master) write(stdout,20) 'SEDIMENT',                          &
     &   'Cohesive and noncohesive sediments.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' SEDIMENT,'
#endif
#ifdef SED_TEST1
      IF (Master) write(stdout,20) 'SED_TEST1',                         &
     &   'Suspended sediment test in a Channel.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SED_TEST1,'
      iexample=iexample+1
#endif
#ifdef SERIAL
      IF (Master) write(stdout,20) 'SERIAL',                            &
     &   'Serial Code Configuration.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' SERIAL,'
#endif
#ifdef SG_BBL
      IF (Master) write(stdout,20) 'SG_BBL',                            &
     &   'Styles and Glenn Bottom Boundary Layer.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' SG_BBL,'
# ifdef SG_LOGINT
      IF (Master) write(stdout,20) 'SG_LOGINT',                         &
     &   'Bottom currents logarithmic interpolation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SG_LOGINT,'
# endif
# ifdef SG_ZNOT
      IF (Master) write(stdout,20) 'SG_ZNOT',                           &
     &   'Internal computation of bottom roughness.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SG_ZNOT,'
# endif
#endif
#ifdef SRELAXATION
      IF (Master) write(stdout,20) 'SRELAXATION',                       &
     &   'Surface salinity relaxation as surface flux.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' SRELAXATION,'
#endif
#ifdef SOLAR_SOURCE

      IF (Master) write(stdout,20) 'SOLAR_SOURCE',                      &
     &   'Solar Radiation Source Term.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' SOLAR_SOURCE,'
#endif
#ifdef SOLITON
      IF (Master) write(stdout,20) 'SOLITON',                           &
     &   'Equatorial Rossby Soliton Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SOLITON,'
      iexample=iexample+1
#endif
#ifdef SOLVE3D
      IF (Master) write(stdout,20) 'SOLVE3D',                           &
     &   'Solving 3D Primitive Equations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SOLVE3D,'
#endif
#ifdef SOUTH_FSCLAMPED
      IF (Master) write(stdout,20) 'SOUTH_FSCLAMPED',                   &
     &   'Southern edge, free-surface, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_FSCLAMPED,'
      FS_south=FS_south+1
#endif
#ifdef SOUTH_FSCHAPMAN
      IF (Master) write(stdout,20) 'SOUTH_FSCHAPMAN',                   &
     &   'Southern edge, free-surface, Chapman condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_FSCHAPMAN,'
      FS_south=FS_south+1
#endif
#ifdef SOUTH_FSGRADIENT
      IF (Master) write(stdout,20) 'SOUTH_FSGRADIENT',                  &
     &   'Southern edge, free-surface, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SOUTH_FSGRADIENT,'
      FS_south=FS_south+1
#endif
#ifdef SOUTH_FSRADIATION
# ifdef SOUTH_FSNUDGING
      IF (Master) write(stdout,20) 'SOUTH_FSNUDGING',                   &
     &   'Southern edge, free-surface, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_FSNUDGING,'
# endif
      IF (Master) write(stdout,20) 'SOUTH_FSRADIATION',                 &
     &   'Southern edge, free-surface, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' SOUTH_FSRADIATION,'
      FS_south=FS_south+1
#endif
#ifdef SOUTH_M2CLAMPED
      IF (Master) write(stdout,20) 'SOUTH_M2CLAMPED',                   &
     &   'Southern edge, 2D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M2CLAMPED,'
      M2_south=M2_south+1
#endif
#ifdef SOUTH_M2FLATHER
      IF (Master) write(stdout,20) 'SOUTH_M2FLATHER',                   &
     &   'Southern edge, 2D momentum, Flather condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M2FLATHER,'
      M2_south=M2_south+1
#endif
#ifdef SOUTH_M2GRADIENT
      IF (Master) write(stdout,20) 'SOUTH_M2GRADIENT',                  &
     &   'Southern edge, 2D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SOUTH_M2GRADIENT,'
      M2_south=M2_south+1
#endif
#ifdef SOUTH_M2RADIATION
# ifdef SOUTH_M2NUDGING
      IF (Master) write(stdout,20) 'SOUTH_M2NUDGING',                   &
     &   'Southern edge, 2D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M2NUDGING,'
# endif
      IF (Master) write(stdout,20) 'SOUTH_M2RADIATION',                 &
     &   'Southern edge, 2D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' SOUTH_M2RADIATION,'
      M2_south=M2_south+1
#endif
#ifdef SOUTH_M3CLAMPED
      IF (Master) write(stdout,20) 'SOUTH_M3CLAMPED',                   &
     &   'Southern edge, 3D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M3CLAMPED,'
      M3_south=M3_south+1
#endif
#ifdef SOUTH_M3GRADIENT
      IF (Master) write(stdout,20) 'SOUTH_M3GRADIENT',                  &
     &   'Southern edge, 3D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SOUTH_M3GRADIENT,'
      M3_south=M3_south+1
#endif
#ifdef SOUTH_M3RADIATION
# ifdef SOUTH_M3NUDGING
      IF (Master) write(stdout,20) 'SOUTH_M3NUDGING',                   &
     &   'Southern edge, 3D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_M3NUDGING,'
# endif
      IF (Master) write(stdout,20) 'SOUTH_M3RADIATION',                 &
     &   'Southern edge, 3D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+19)=' SOUTH_M3RADIATION,'
      M3_south=M3_south+1
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# ifdef SOUTH_KGRADIENT
      IF (Master) write(stdout,20) 'SOUTH_KGRADIENT',                   &
     &   'Southern edge, TKE fields, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_KGRADIENT,'
      K_south=K_south+1
# endif
# ifdef SOUTH_KRADIATION
      IF (Master) write(stdout,20) 'SOUTH_KRADIATION',                  &
     &   'Southern edge, TKE fields, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_KRADIATION,'
      K_south=K_south+1
# endif
#endif
#ifdef SOUTH_TCLAMPED
      IF (Master) write(stdout,20) 'SOUTH_TCLAMPED',                    &
     &   'Southern edge, tracers, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_TCLAMPED,'
      T_south=T_south+1
#endif
#ifdef SOUTH_TGRADIENT
      IF (Master) write(stdout,20) 'SOUTH_TGRADIENT',                   &
     &   'Southern edge, tracers, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' SOUTH_TGRADIENT,'
      T_south=T_south+1
#endif
#ifdef SOUTH_TRADIATION
# ifdef SOUTH_TNUDGING
      IF (Master) write(stdout,20) 'SOUTH_TNUDGING',                    &
     &   'Southern edge, tracers, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' SOUTH_TNUDGING,'
# endif
      IF (Master) write(stdout,20) 'SOUTH_TRADIATION',                  &
     &   'Southern edge, tracers, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' SOUTH_TRADIATION,'
      T_south=T_south+1
#endif
#ifdef SOUTHERN_WALL
      IF (Master) write(stdout,20) 'SOUTHERN_WALL',                     &
     &   'Wall boundary at Southern edge.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' SOUTHERN_WALL,'
#endif
#ifdef SPLINES
      IF (Master) write(stdout,20) 'SPLINES',                           &
     &   'Conservative parabolic spline reconstruction.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' SPLINES,'
#endif
#ifdef SPHERICAL
      IF (Master) write(stdout,20) 'SPHERICAL',                         &
     &   'Spherical grid configuration.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SPHERICAL,'
#endif
#ifdef SPONGE
      IF (Master) write(stdout,20) 'SPONGE',                            &
     &   'Enhanced horizontal mixing in the sponge areas.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' SPONGE,'
#endif
#if defined SSH_TIDES
      IF (Master) write(stdout,20) 'SSH_TIDES',                         &
     &   'Add tidal elevation to SSH climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' SSH_TIDES,'
#endif
#ifdef STATIONS
      IF (Master) write(stdout,20) 'STATIONS',                          &
     &   'Writing out station data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' STATIONS,'
#endif
#if defined STATIONS_CGRID && defined STATIONS
      IF (Master) write(stdout,20) 'STATIONS_CGRID',                    &
     &   'Extracting station data at native C-grid locations.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' STATIONS_CGRID,'
#endif
#if defined TCLIMATOLOGY && defined SOLVE3D
      IF (Master) write(stdout,20) 'TCLIMATOLOGY',                      &
     &   'Processing tracer climatology data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' TCLIMATOLOGY,'
#endif
#if defined TCLM_NUDGING && defined SOLVE3D
      IF (Master) write(stdout,20) 'TCLM_NUDGING',                      &
     &   'Nudging toward tracer climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' TCLM_NUDGING,'
#endif
#ifdef SOLVE3D
# if defined TS_A4HADVECTION
      IF (Master) write(stdout,20) 'TS_A4HADVECTION',                   &
     &   'Fouth-order Akima horizontal advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_A4HADVECTION,'
# elif defined TS_C2HADVECTION
      IF (Master) write(stdout,20) 'TS_C2HADVECTION',                   &
     &   'Second-order centered horizontal advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C2HADVECTION,'
# elif defined TS_U3HADVECTION
      IF (Master) write(stdout,20) 'TS_U3HADVECTION',                   &
     &   'Third-order upstream bias horizontal advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_U3HADVECTION,'
# else
      IF (Master) write(stdout,20) 'TS_C4HADVECTION',                   &
     &   'Fourth-order centered horizontal advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C4HADVECTION,'
# endif
# if defined TS_A4VADVECTION
      IF (Master) write(stdout,20) 'TS_A4VADVECTION',                   &
     &   'Fouth-order Akima vertical advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_A4VADVECTION,'
# elif defined TS_C2VADVECTION
      IF (Master) write(stdout,20) 'TS_C2VADVECTION',                   &
     &   'Second-order centered vertical advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C2VADVECTION,'
# elif defined TS_SVADVECTION
      IF (Master) write(stdout,20) 'TS_SVADVECTION',                    &
     &   'Parabolic splines vertical advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' TS_SVADVECTION,'
# else
      IF (Master) write(stdout,20) 'TS_C4VADVECTION',                   &
     &   'Fourth-order centered vertical advection of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' TS_C4VADVECTION,'
# endif
#endif
#if defined TS_DIF2 && defined SOLVE3D
      IF (Master) write(stdout,20) 'TS_DIF2',                           &
     &   'Harmonic mixing of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' TS_DIF2,'
#endif
#if defined TS_DIF4 && defined SOLVE3D
      IF (Master) write(stdout,20) 'TS_DIF4',                           &
     &   'Biharmonic mixing of tracers.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' TS_DIF4,'
#endif
#if defined TS_PSOURCE && defined SOLVE3D
      IF (Master) write(stdout,20) 'TS_PSOURCE',                        &
     &   'Tracers point sources and sinks.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' TS_PSOURCE,'
#endif
#ifdef UPWELLING
      IF (Master) write(stdout,20) 'UPWELLING',                         &
     &   'Upwelling/Downwelling Example.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' UPWELLING,'
      iexample=iexample+1
#endif
#ifdef USWEST
      IF (Master) write(stdout,20) 'USWEST',                            &
     &   'US West Coast Application.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' USWEST,'
      iexample=iexample+1
#endif
#ifdef UV_ADV
      IF (Master) write(stdout,20) 'UV_ADV',                            &
     &   'Advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_ADV,'
#endif
#if defined UV_BAROCLINIC && \
   (defined ASSIMILATION_UVsur || defined ASSIMILATION_UV)
      IF (Master) write(stdout,20) 'UV_BAROCLINIC',                     &
     &   'Assimilation of only baroclinic currents.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' UV_BAROCLINIC,'
#endif
#ifdef UV_COR
      IF (Master) write(stdout,20) 'UV_COR',                            &
     &   'Coriolis term.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+8)=' UV_COR,'
#endif
#ifdef UV_ADV
# if defined UV_C2ADVECTION
      IF (Master) write(stdout,20) 'UV_C2ADVECTION',                    &
     &   'Second-order centered differences advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_C2ADVECTION,'
# elif defined UV_C4ADVECTION
      IF (Master) write(stdout,20) 'UV_C4ADVECTION',                    &
     &   'Forth-order centered differences advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_C4ADVECTION,'
# else
#  ifdef SOLVE3D
      IF (Master) write(stdout,20) 'UV_U3HADVECTION',                   &
     &   'Third-order upstream bias advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' UV_U3HADVECTION,'
#  else
      IF (Master) write(stdout,20) 'UV_C4ADVECTION',                    &
     &   'Forth-order centered differences advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' UV_C4ADVECTION,'
#  endif
# endif
# ifdef UV_SADVECTION
      IF (Master) write(stdout,20) 'UV_SADVECTION',                     &
     &   'Parabolic splines vertical advection of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' UV_SADVECTION,'
# endif
#endif
#ifdef UV_PSOURCE
      IF (Master) write(stdout,20) 'UV_PSOURCE',                        &
     &   'Mass point sources and sinks.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' UV_PSOURCE,'
#endif
#if defined UV_TIDES
      IF (Master) write(stdout,20) 'UV_TIDES',                          &
     &   'Add tidal currents to 2D momentum climatologies.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' UV_TIDES,'
#endif
#ifdef UV_VIS2
      IF (Master) write(stdout,20) 'UV_VIS2',                           &
     &   'Harmonic mixing of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' UV_VIS2,'
#endif
#ifdef UV_VIS4
      IF (Master) write(stdout,20) 'UV_VIS4',                           &
     &   'Biharmonic mixing of momentum.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' UV_VIS4,'
#endif
#if defined VAR_RHO_2D && defined SOLVE3D
      IF (Master) write(stdout,20) 'VAR_RHO_2D',                        &
     &   'Variable density barotropic mode.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+12)=' VAR_RHO_2D,'
#endif
#if defined UV_VIS2 || defined UV_VIS4
# ifdef VISC_GRID
      IF (Master) write(stdout,20) 'VISC_GRID',                         &
     &   'Horizontal viscosity coefficient scaled by grid size.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+11)=' VISC_GRID,'
# endif
#endif
#if defined WRITE_WATER && defined MASKING
      IF (Master) write(stdout,20) 'WRITE_WATER',                       &
     &   'Writing data at water points only.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+13)=' WRITE_WATER,'
#endif
#ifdef WEDDELL
      IF (Master) write(stdout,20) 'WEDDELL',                           &
     &   'Idealized Weddell Sea Shelf Application.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+9)=' WEDDELL,'
      iexample=iexample+1
#endif
#ifdef WEST_FSCLAMPED
      IF (Master) write(stdout,20) 'WEST_FSCLAMPED',                    &
     &   'Western edge, free-surface, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_FSCLAMPED,'
      FS_west=FS_west+1
#endif
#ifdef WEST_FSCHAPMAN
      IF (Master) write(stdout,20) 'WEST_FSCHAPMAN',                    &
     &   'Western edge, free-surface, Chapman condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_FSCHAPMAN,'
      FS_west=FS_west+1
#endif
#ifdef WEST_FSGRADIENT
      IF (Master) write(stdout,20) 'WEST_FSGRADIENT',                   &
     &   'Western edge, free-surface, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_FSGRADIENT,'
      FS_west=FS_west+1
#endif
#ifdef WEST_FSRADIATION
# ifdef WEST_FSNUDGING
      IF (Master) write(stdout,20) 'WEST_FSNUDGING',                    &
     &   'Western edge, free-surface, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_FSNUDGING,'
# endif
      IF (Master) write(stdout,20) 'WEST_FSRADIATION',                  &
     &   'Western edge, free-sruface, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' WEST_FSRADIATION,'
      FS_west=FS_west+1
#endif
#ifdef WEST_M2CLAMPED
      IF (Master) write(stdout,20) 'WEST_M2CLAMPED',                    &
     &   'Western edge, 2D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M2CLAMPED,'
      M2_west=M2_west+1
#endif
#ifdef WEST_M2FLATHER
      IF (Master) write(stdout,20) 'WEST_M2FLATHER',                    &
     &   'Western edge, 2D momentum, Flather condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M2FLATHER,'
      M2_west=M2_west+1
#endif
#ifdef WEST_M2GRADIENT
      IF (Master) write(stdout,20) 'WEST_M2GRADIENT',                   &
     &   'Western edge, 2D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_M2GRADIENT,'
      M2_west=M2_west+1
#endif
#ifdef WEST_M2RADIATION
# ifdef WEST_M2NUDGING
      IF (Master) write(stdout,20) 'WEST_M2NUDGING',                    &
     &   'Western edge, 2D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M2NUDGING,'
# endif
      IF (Master) write(stdout,20) 'WEST_M2RADIATION',                  &
     &   'Western edge, 2D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' WEST_M2RADIATION,'
      M2_west=M2_west+1
#endif
#ifdef WEST_M3CLAMPED
      IF (Master) write(stdout,20) 'WEST_M3CLAMPED',                    &
     &   'Western edge, 3D momentum, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M3CLAMPED,'
      M3_west=M3_west+1
#endif
#ifdef WEST_M3GRADIENT
      IF (Master) write(stdout,20) 'WEST_M3GRADIENT',                   &
     &   'Western edge, 3D momentum, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_M3GRADIENT,'
      M3_west=M3_west+1
#endif
#ifdef WEST_M3RADIATION
# ifdef WEST_M3NUDGING
      IF (Master) write(stdout,20) 'WEST_M3NUDGING',                    &
     &   'Western edge, 3D momentum, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_M3NUDGING,'
# endif
      IF (Master) write(stdout,20) 'WEST_M3RADIATION',                  &
     &   'Western edge, 3D momentum, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+18)=' WEST_M3RADIATION,'
      M3_west=M3_west+1
#endif
#if defined GLS_MIXING || defined MY25_MIXING
# ifdef WEST_KGRADIENT
      IF (Master) write(stdout,20) 'WEST_KGRADIENT',                    &
     &   'Western edge, TKE fields, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_KGRADIENT,'
      K_west=K_west+1
# endif
# ifdef WEST_KRADIATION
      IF (Master) write(stdout,20) 'WEST_KRADIATION',                   &
     &   'Western edge, TKE fields, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_KRADIATION,'
      K_west=K_west+1
# endif
#endif
#ifdef WEST_TCLAMPED
      IF (Master) write(stdout,20) 'WEST_TCLAMPED',                     &
     &   'Western edge, tracers, Clamped condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_TCLAMPED,'
      T_west=T_west+1
#endif
#ifdef WEST_TGRADIENT
      IF (Master) write(stdout,20) 'WEST_TGRADIENT',                    &
     &   'Western edge, tracers, gradient condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+16)=' WEST_TGRADIENT,'
      T_west=T_west+1
#endif
#ifdef WEST_TRADIATION
# ifdef WEST_TNUDGING
      IF (Master) write(stdout,20) 'WEST_TNUDGING',                     &
     &   'Western edge, tracers, passive/active outflow/inflow.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+15)=' WEST_TNUDGING,'
# endif
      IF (Master) write(stdout,20) 'WEST_TRADIATION',                   &
     &   'Western edge, tracers, radiation condition.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+17)=' WEST_TRADIATION,'
      T_west=T_west+1
#endif
#ifdef WESTERN_WALL
      IF (Master) write(stdout,20) 'WESTERN_WALL',                      &
     &   'Wall boundary at Western edge.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' WESTERN_WALL,'
#endif
#ifdef WJ_GRADP
      IF (Master) write(stdout,20) 'WJ_GRADP',                          &
     &   'Weighted Jacobians pressure gradient formulation.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+10)=' WJ_GRADP,'
#endif
#ifdef ZCLIMATOLOGY
      IF (Master) write(stdout,20) 'ZCLIMATOLOGY',                      &
     &   'Processing sea surface height climatology data.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ZCLIMATOLOGY,'
#endif
#ifdef ZCLM_NUDGING
      IF (Master) write(stdout,20) 'ZCLM_NUDGING',                      &
     &   'Nudging toward sea surface height climatology.'
      is=LEN_TRIM(Coptions)+1
      Coptions(is:is+14)=' ZCLM_NUDGING,'
#endif
!
!-----------------------------------------------------------------------
!  Check C-preprocessing options.
!-----------------------------------------------------------------------
!
!  Stop if more than one vertical closure scheme is selected.
!
      IF (Master.and.(ivmix.gt.1)) THEN
        write (stdout,30)
  30    format(/,' CHECKDEFS - only one vertical closure scheme',       &
     &                       ' is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more than one biological module is selected.
!
      IF (Master.and.(ibiology.gt.1)) THEN
        write (stdout,40)
  40    format(/,' CHECKDEFS - only one biology MODULE is allowed.')
        exit_flag=5
      END IF

!
!  Stop if more than one model example is selected.
!
      IF (Master.and.(iexample.gt.1)) THEN
        write (stdout,50)
  50    format(/,' CHECKDEFS - only one model example is allowed.')
        exit_flag=5
      END IF
!
!  Stop if more than one boundary conditions option for free-surface
!  is selected at each boundary edge.
!
      IF (Master.and.(FS_east.gt.1)) THEN
        write (stdout,60) 'free-surface','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(FS_west.gt.1)) THEN
        write (stdout,60) 'free-surface','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(FS_north.gt.1)) THEN
        write (stdout,60) 'free-surface','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(FS_south.gt.1)) THEN
        write (stdout,60) 'free-surface','southern edge'
        exit_flag=5
      END IF
  60  format(/,' CHECKDEFS - only one boundary condition option for',   &
     &       1x,a,/,13x,'is allowed at the',1x,a,'.')
!
!  Stop if more than one boundary conditions option for 2D momentum
!  is selected at each boundary edge.
!
      IF (Master.and.(M2_east.gt.1)) THEN
        write (stdout,60) '2D momentum','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M2_west.gt.1)) THEN
        write (stdout,60) '2D momentum','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(M2_north.gt.1)) THEN
        write (stdout,60) '2D momentum','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M2_south.gt.1)) THEN
        write (stdout,60) '2D momentum','southern edge'
        exit_flag=5
      END IF
#ifdef SOLVE3D
!
!  Stop if more than one boundary conditions option for 3D momentum
!  is selected at each boundary edge.
!
      IF (Master.and.(M3_east.gt.1)) THEN
        write (stdout,60) '3D momentum','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M3_west.gt.1)) THEN
        write (stdout,60) '3D momentum','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(M3_north.gt.1)) THEN
        write (stdout,60) '3D momentum','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(M3_south.gt.1)) THEN
        write (stdout,60) '3D momentum','southern edge'
        exit_flag=5
      END IF
!
!  Stop if more than one boundary conditions option for tracers is
!  selected at each boundary edge.
!
      IF (Master.and.(T_east.gt.1)) THEN
        write (stdout,60) 'tracers','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(T_west.gt.1)) THEN
        write (stdout,60) 'tracers','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(T_north.gt.1)) THEN
        write (stdout,60) 'tracers','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(T_south.gt.1)) THEN
        write (stdout,60) 'tracers','southern edge'
        exit_flag=5
      END IF
# if defined GLS_MIXING || defined MY25_MIXING
!
!  Stop if more than one boundary conditions option for tracers is
!  selected at each boundary edge.
!
      IF (Master.and.(K_east.gt.1)) THEN
        write (stdout,60) 'TKE fields','eastern edge'
        exit_flag=5
      END IF
      IF (Master.and.(K_west.gt.1)) THEN
        write (stdout,60) 'TKE fields','western edge'
        exit_flag=5
      END IF
      IF (Master.and.(K_north.gt.1)) THEN
        write (stdout,60) 'TKE fields','northern edge'
        exit_flag=5
      END IF
      IF (Master.and.(K_south.gt.1)) THEN
        write (stdout,60) 'TKE fields','southern edge'
        exit_flag=5
      END IF
# endif
#endif
#ifndef SOLVE3D
!
!  Stop it explicit time-step splitting on shallow water set-up.
!
      DO ng=1,Ngrids
        IF (Master.and.(ndtfast(ng).gt.1)) THEN
          write(stdout,70)
  70      format(/,' CHECKDEFS - explicit time-step splitting is ',     &
     &                         ' inconsistent.',                        &
     &           /,13x,'Change parameter NDTFAST to unity.')
          exit_flag=5
        END IF
      END DO
#endif
      RETURN
      END SUBROUTINE checkdefs
