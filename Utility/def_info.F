#include "cppdefs.h"
      function def_info (ng,ncid,ncname,DimIDs)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine defines information variables in requested NetCDF      !
!  file.                                                               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer).                           !
!     ncid     NetCDF file ID (integer).                               !
!     ncname   NetCDF file name (character).                           !
!     DimIDs   NetCDF dimensions IDs (integer vector):                 !
!                DimIDs( 1) => XI-dimension at RHO-points.             !
!                DimIDs( 2) => XI-dimension at U-points.               !
!                DimIDs( 3) => XI-dimension at V-points.               !
!                DimIDs( 4) => XI-dimension at PSI-points.             !
!                DimIDs( 5) => ETA-dimension at RHO-points.            !
!                DimIDs( 6) => ETA-dimension at U-points.              !
!                DimIDs( 7) => ETA-dimension at V-points.              !
!                DimIDs( 8) => ETA-dimension at PSI-points.            !
!                DimIDs( 9) => S-dimension at RHO-points.              !
!                DimIDs(10) => S-dimension at W-points.                !
!                DimIDs(11) => Number of tracers dimension.            !
!                DimIDs(12) => Unlimited time record dimension.        !
!                DimIDs(13) => Number of stations dimension.           !
!                DimIDs(14) => Boundary dimension.                     !
!                DimIDs(15) => Number of floats dimension.             !
!                DimIDs(16) => Number sediment bed layers dimension.   !
!                DimIDs(17) => Dimension 2D water RHO-points.          !
!                DimIDs(18) => Dimension 2D water U-points.            !
!                DimIDs(19) => Dimension 2D water V-points.            !
!                DimIDs(20) => Dimension 3D water RHO-points.          !
!                DimIDs(21) => Dimension 3D water U-points.            !
!                DimIDs(23) => Dimension 3D water W-points.            !
!                DimIDs(24) => Dimension sediment bed water points.    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     def_info   Error flag (INTEGER_TYPE).                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_strings

      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ncid
      integer, intent(in) :: DimIDs(24)
      character (*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer, parameter :: Natt = 20

      integer :: brydim, i, j, srdim, stadim, status, swdim
      integer :: trcdim, usrdim, varid

      integer :: p2dgrd(2), tbrydim(2)
      integer :: t2dgrd(3), u2dgrd(3), v2dgrd(3)

      integer :: def_info, def_var

      real(r8) :: Aval(4)

      character (len=11) :: frcatt
      character (len=50) :: tiling
      character (len=65) :: Vinfo(Natt)
!
!-----------------------------------------------------------------------
!  Set dimension variables.
!-----------------------------------------------------------------------
!
      p2dgrd(1)=DimIDs(4)
      p2dgrd(2)=DimIDs(8)
      t2dgrd(1)=DimIDs(1)
      t2dgrd(2)=DimIDs(5)
      u2dgrd(1)=DimIDs(2)
      u2dgrd(2)=DimIDs(6)
      v2dgrd(1)=DimIDs(3)
      v2dgrd(2)=DimIDs(7)
      srdim=DimIDs(9)
      swdim=DimIDs(10)
      trcdim=DimIDs(11)
      stadim=DimIDs(13)
      brydim=DimIDs(14)
      tbrydim(1)=DimIDs(11)
      tbrydim(2)=DimIDs(14)
!
!  Set dimension for generic user parameters.
!
      IF (Nuser.gt.0) THEN
        status=nf_def_dim(ncid,'Nuser',25,usrdim)
      END IF
!
!  Initialize local information variable arrays.
!
      DO i=1,Natt
        DO j=1,65
          Vinfo(i)(j:j)=' '
        END DO
      END DO
      DO i=1,4
        Aval(i)=0.0_r8
      END DO
!
!-----------------------------------------------------------------------
!  Define global attributes.
!-----------------------------------------------------------------------
!
!  Define history global attribute.
!
      IF (LEN_TRIM(date_str).gt.0) THEN
        history='ROMS/TOMS, Version '//version//', '//trim(date_str)
      ELSE
        history='ROMS/TOMS, Version '//version
      END IF
!
!  Define tile decomposition global attribute.
!
      WRITE (tiling,10) NtileI(ng), NtileJ(ng)
  10  FORMAT (i3.3,'x',i3.3)
!
!  Define file type global attribute.
!
      IF (ncid.eq.ncavgid(ng)) THEN
        type='ROMS/TOMS averages file'
      ELSE IF (ncid.eq.nchisid(ng)) THEN
        type='ROMS/TOMS hystory file'
      ELSE IF (ncid.eq.ncrstid(ng)) THEN
        type='ROMS/TOMS restart file'
      ELSE IF (ncid.eq.ncstaid(ng)) THEN
        type='ROMS/TOMS station file'
      ELSE IF (ncid.eq.ncfltid(ng)) THEN
        type='ROMS/TOMS floats file'
      END IF
      status=nf_put_att_text(ncid,nf_global,'type',LEN_TRIM(type),type)
!
!  Define other global attributes to NetCDF file.
!
      status=nf_put_att_text(ncid,nf_global,'title',                    &
     &                       LEN_TRIM(title),title)
      status=nf_put_att_text(ncid,nf_global,'var_info',                 &
     &                       LEN_TRIM(varname),varname)
      status=nf_put_att_text(ncid,nf_global,'rst_file',                 &
     &                       LEN_TRIM(rstname(ng)),rstname(ng))
      IF (LdefHIS(ng)) THEN
        IF (ndefhis(ng).gt.0) THEN
          status=nf_put_att_text(ncid,nf_global,'his_base',             &
     &                           LEN_TRIM(hisbase(ng))-3,hisbase(ng))
        ELSE
          status=nf_put_att_text(ncid,nf_global,'his_file',             &
     &                           LEN_TRIM(hisname(ng)),hisname(ng))
        END IF
      END IF
#ifdef AVERAGES
      IF (ndefavg(ng).gt.0) THEN
        status=nf_put_att_text(ncid,nf_global,'avg_base',               &
     &                         LEN_TRIM(avgbase(ng))-3,avgbase(ng))
      ELSE
        status=nf_put_att_text(ncid,nf_global,'avg_file',               &
     &                         LEN_TRIM(avgname(ng)),avgname(ng))
      END IF
#endif
#ifdef STATIONS
      status=nf_put_att_text(ncid,nf_global,'sta_file',                 &
     &                       LEN_TRIM(staname(ng)),staname(ng))
#endif
#ifdef FLOATS
      status=nf_put_att_text(ncid,nf_global,'flt_file',                 &
     &                       LEN_TRIM(fltname(ng)),fltname(ng))
#endif
#ifndef ANA_GRID
      status=nf_put_att_text(ncid,nf_global,'grd_file',                 &
     &                       LEN_TRIM(grdname(ng)),grdname(ng))
#endif
#ifndef ANA_INITIAL
      status=nf_put_att_text(ncid,nf_global,'ini_file',                 &
     &                       LEN_TRIM(ininame(ng)),ininame(ng))
#endif
#ifdef FRC_FILE
      DO i=1,nFfiles(ng)
        WRITE (frcatt,20) i
  20    FORMAT ('frc_file_',i2.2)
        status=nf_put_att_text(ncid,nf_global,frcatt,                   &
     &                       LEN_TRIM(frcname(i,ng)),frcname(i,ng))
      END DO
#endif
#if (defined  ZCLIMATOLOGY && !defined ANA_SSH)     || \
    (defined M2CLIMATOLOGY && !defined ANA_M2CLIMA) || \
    (defined  TCLIMATOLOGY && !defined ANA_TCLIMA)  || \
    (defined M3CLIMATOLOGY && !defined ANA_M3CLIMA)
      status=nf_put_att_text(ncid,nf_global,'clm_file',                 &
     &                       LEN_TRIM(clmname(ng)),clmname(ng))
#endif
#if defined ASSIMILATION || defined NUDGING
      status=nf_put_att_text(ncid,nf_global,'apar_file',                &
     &                       LEN_TRIM(aparnam),aparnam)
#endif
#ifdef BIOLOGY
      status=nf_put_att_text(ncid,nf_global,'bpar_file',                &
     &                       LEN_TRIM(bparnam),bparnam)
#endif
#ifdef FLOATS
      status=nf_put_att_text(ncid,nf_global,'fpos_file',                &
     &                       LEN_TRIM(fposnam),fposnam)
#endif
#ifdef STATIONS
      status=nf_put_att_text(ncid,nf_global,'spos_file',                &
     &                       LEN_TRIM(sposnam),sposnam)
#endif
#ifndef DEBUGGING
!
! Avoid writing this information in order to use "diff" to compare
! NetCDF files during the debugging phase.
!
      status=nf_put_att_text(ncid,nf_global,'history',                  &
     &                       LEN_TRIM(history),history)
      status=nf_put_att_text(ncid,nf_global,'tiling',                   &
     &                       LEN_TRIM(tiling),tiling)
      status=nf_put_att_text(ncid,nf_global,'CPP_options',              &
     &                       LEN_TRIM(Coptions),Coptions)
#endif
!
!-----------------------------------------------------------------------
!  Define running parameters.
!-----------------------------------------------------------------------
!
!  Time stepping parameters.
!
      Vinfo( 1)='ntimes'
      Vinfo( 2)='number of long time-steps'
      status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ndtfast'
      Vinfo( 2)='number of short time-steps'
      status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='dt'
      Vinfo( 2)='size of long time-steps'
      Vinfo( 3)='second'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='dtfast'
      Vinfo( 2)='size of short time-steps'
      Vinfo( 3)='second'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='dstart'
      Vinfo( 2)='time stamp assigned to model initilization'
      IF (INT(time_ref).eq.-2) THEN
        Vinfo( 3)='days since 1968-05-23 00:00:00 GMT'
        Vinfo( 4)='modified Julian day number'
        Vinfo(18)='add_offset'
        Aval(1)=jul_off
      ELSE IF (INT(time_ref).eq.-1) THEN
        Vinfo( 3)='days since 0000-01-01 00:00:00'
        Vinfo( 4)='360 days in every year, 30 days in each month'
      ELSE IF (INT(time_ref).eq.0) THEN
        Vinfo( 3)='days since 0000-01-01 00:00:00'
        Vinfo( 4)='365.25 days in every year'
      ELSE IF (time_ref.gt.0.0_r8) THEN
        Vinfo( 3)='days since '//r_text
      END IF
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='nhis'
      Vinfo( 2)='number of time-steps between history records'
      status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)

      IF (ndefhis(ng).gt.0) THEN
        Vinfo( 1)='ndefhis'
        Vinfo( 2)=                                                      &
     &    'number of time-steps between the creation of history files'
          status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)
      END IF

      Vinfo( 1)='nrst'
      Vinfo( 2)='number of time-steps between restart records'
      IF (lcycle(ng)) THEN
        Vinfo(13)='only latest two records are maintained'
      END IF
      status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)

#ifdef AVERAGES
      Vinfo( 1)='ntsavg'
      Vinfo( 2)=                                                        &
     &   'starting time-step for accumulation of time-averaged fields'
      status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='navg'
      Vinfo( 2)='number of time-steps between time-averaged records'
      status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)

      IF (ndefavg(ng).gt.0) THEN
        Vinfo( 1)='ndefavg'
        Vinfo( 2)=                                                      &
     &    'number of time-steps between the creation of average files'
        status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)
      END IF
#endif

#ifdef STATIONS
      Vinfo( 1)='nsta'
      Vinfo( 2)='number of time-steps between stations records'
      status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)
#endif

#if defined POWER_LAW && defined SOLVE3D
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
      Vinfo( 1)='Falpha'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Fbeta'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Fgamma'
      Vinfo( 2)='Power-law shape barotropic filter parameter'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
#endif
!
!  Horizontal mixing coefficients.
!
#if defined SOLVE3D && defined TS_DIF2
      Vinfo( 1)='tnu2'
      Vinfo( 2)='Laplacian mixing coefficient for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)
#endif

#if defined SOLVE3D && defined TS_DIF4
      Vinfo( 1)='tnu4'
      Vinfo( 2)='biharmonic mixing coefficient for tracers'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)
#endif

#ifdef UV_VIS2
      Vinfo( 1)='visc2'
      Vinfo( 2)='Laplacian mixing coefficient for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
#endif
#ifdef UV_VIS4
      Vinfo( 1)='visc4'
      Vinfo( 2)='biharmonic mixing coefficient for momentum'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
#endif

#if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
# ifdef TKE_DIF2
      Vinfo( 1)='tkenu2'
      Vinfo( 2)='harmonic mixing coefficient for turbulent energy'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
# endif

# ifdef TKE_DIF4
      Vinfo( 1)='tkenu4'
      Vinfo( 2)='biharmonic mixing coefficient for turbulent energy'
      Vinfo( 3)='meter4 second-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
# endif
#endif
#ifdef SOLVE3D
!
!  Background vertical mixing coefficients.
!
      Vinfo( 1)='Akt_bak'
      Vinfo( 2)='background vertical mixing coefficient for tracers'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)

      Vinfo( 1)='Akv_bak'
      Vinfo( 2)='background vertical mixing coefficient for momentum'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

# if defined MY25_MIXING || defined GLS_MIXING
      Vinfo( 1)='Akk_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for turbulent energy'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Akp_bak'
      Vinfo( 2)=                                                        &
     &   'background vertical mixing coefficient for length scale'
      Vinfo( 3)='meter2 second-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
# endif
#endif
!
!  Drag coefficients.
!
      Vinfo( 1)='rdrg'
      Vinfo( 2)='linear drag coefficient'
      Vinfo( 3)='meter second-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='rdrg2'
      Vinfo( 2)='quadratic drag coefficient'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

#ifdef SOLVE3D
      Vinfo( 1)='Zob'
      Vinfo( 2)='bottom roughness'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Zos'
      Vinfo( 2)='surface roughness'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
#endif
#if defined SOLVE3D && defined GLS_MIXING
!
!  Generic length-scale parameters.
!
      Vinfo( 1)='gls_p'
      Vinfo( 2)='stability exponent'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_m'
      Vinfo( 2)='turbulent kinetic energy exponent'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_n'
      Vinfo( 2)='turbulent length scale exponent'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_cmu0'
      Vinfo( 2)='stability coefficient'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_c1'
      Vinfo( 2)='shear production coefficient'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_c2'
      Vinfo( 2)='dissipation coefficient'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_c3m'
      Vinfo( 2)='buoyancy production coefficient (minus)'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_c3p'
      Vinfo( 2)='buoyancy production coefficient (plus)'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_sigk'
      Vinfo( 2)='constant Schmidt number for TKE'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_sigp'
      Vinfo( 2)='constant Schmidt number for PSI'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_Kmin'
      Vinfo( 2)='stability exponent'
      Vinfo( 3)='minimum value of specific turbulent kinetic energy'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='gls_Pmin'
      Vinfo( 2)='minimum Value of dissipation'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
#endif
!
!  Nudging inverse time scales used in various tasks.
!
      Vinfo( 1)='Znudg'
      Vinfo( 2)='free-surface nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='M2nudg'
      Vinfo( 2)='2D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

#ifdef SOLVE3D
      Vinfo( 1)='M3nudg'
      Vinfo( 2)='3D momentum nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Tnudg'
      Vinfo( 2)='Tracers nudging/relaxation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)
#endif
#ifdef NUDGING
!
!  Nudging inverse time scales used in data assimilation.
!
      Vinfo( 1)='Znudass'
      Vinfo( 2)='free-surface nudging assimilation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='M2nudass'
      Vinfo( 2)='2D momentum nudging assimilation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

# ifdef SOLVE3D
      Vinfo( 1)='M3nudass'
      Vinfo( 2)='3D momentum nudging assimilation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Tnudass'
      Vinfo( 2)='Tracers nudging assimilation inverse time scale'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,1,trcdim,Aval,Vinfo,ncname)
# endif
#endif
#ifdef NUDGING_COFF
!
!  Open boundary nudging, inverse time scales.
!
      Vinfo( 1)='FSobc_in'
      Vinfo( 2)='free-surface inflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)

      Vinfo( 1)='FSobc_out'
      Vinfo( 2)='free-surface outflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)

      Vinfo( 1)='M2obc_in'
      Vinfo( 2)='2D momentum inflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)

      Vinfo( 1)='M2obc_out'
      Vinfo( 2)='2D momentum outflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)

# ifdef SOLVE3D
      Vinfo( 1)='Tobc_in'
      Vinfo( 2)='tracers inflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,2,tbrydim,Aval,Vinfo,ncname)

      Vinfo( 1)='Tobc_out'
      Vinfo( 2)='tracers outflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,2,tbrydim,Aval,Vinfo,ncname)

      Vinfo( 1)='M3obc_in'
      Vinfo( 2)='3D momentum inflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)

      Vinfo( 1)='M3obc_out'
      Vinfo( 2)='3D momentum outflow, nudging inverse time scale'
      Vinfo( 3)='second-1'
      status=def_var(ncid,varid,NF_TYPE,1,brydim,Aval,Vinfo,ncname)
# endif
#endif
!
!  Equation of State parameters.
!
      Vinfo( 1)='rho0'
      Vinfo( 2)='mean density used in Boussinesq approximation'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

#if defined SOLVE3D && !defined NONLIN_EOS
      Vinfo( 1)='R0'
      Vinfo( 2)='background density used in linear equation of state'
      Vinfo( 3)='kilogram meter-3'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Tcoef'
      Vinfo( 2)='thermal expansion coefficient'
      Vinfo( 3)='Celsius-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Scoef'
      Vinfo( 2)='Saline contraction coefficient'
      Vinfo( 3)='PSU-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
#endif
#ifdef SOLVE3D
!
!  Various parameters.
!
# ifdef BODYFORCE
      Vinfo( 1)='levsfrc'
      Vinfo( 2)='shallowest level for body-force stress'
      status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='levbfrc'
      Vinfo( 2)='deepest level for body-force stress'
      status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)
# endif
#endif
!
!  Slipperiness parameters.
!
      Vinfo( 1)='gamma2'
      Vinfo( 2)='Slipperiness parameter'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

#if defined BIOLOGY && defined SOLVE3D
# ifdef BIO_FASHAM
!
!  Define Fasham type biological module parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ncid,varid,nf_int,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='AttSW'
      Vinfo( 2)='light attenuation by seawater'
      Vinfo( 3)='meter-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='AttChl'
      Vinfo( 2)='light attenuation by chlorophyll'
      Vinfo( 3)='meter-2 milligram_Chl-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PARfrac'
      Vinfo( 2)='photosynthetically available radiation fraction'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Vp0'
      Vinfo( 2)='Eppley temperature-limited growth parameter'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='I_thNH4'
      Vinfo( 2)='radiation threshold for nitrification'
      Vinfo( 3)='Watts m-2'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='D_p5NH4'
      Vinfo( 2)='half-saturation radiation for nitrification'
      Vinfo( 3)='Watts m-2'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='NitriR'
      Vinfo( 2)='nitrification rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='K_NO3'
      Vinfo( 2)='inverse half-saturation for phytoplankton NO3 uptake'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='K_NH4'
      Vinfo( 2)='inverse half-saturation for phytoplankton NH4 uptake'
      Vinfo( 3)='meter3 millimole_N-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='K_Phy'
      Vinfo( 2)='zooplankton half-saturation constant for ingestion'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Chl2C_m'
      Vinfo( 2)='maximum chlorophyll to carbon ratio'
      Vinfo( 3)='milligram_chl milligram_carbon-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ChlMin'
      Vinfo( 2)='minimum chlorophyll threshold'
      Vinfo( 3)='milligram_chl meter-3'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyCN'
      Vinfo( 2)='phytoplankton Carbon:Nitrogen ratio'
      Vinfo( 3)='mole_C mole_N-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyIP'
      Vinfo( 2)='phytoplankton, NH4 inhibition parameter'
      Vinfo( 3)='millimole_N-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyIS'
      Vinfo( 2)='phytoplankton, initial slope of P-I curve'
      Vinfo( 3)='milligram_C milligram_Chl-1 Watts-1 meter2 day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyMin'
      Vinfo( 2)='minimum phytoplankton threshold'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='PhyMR'
      Vinfo( 2)='phytoplankton mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

#  ifdef CARBON
      Vinfo( 1)='ZooAE_C'
      Vinfo( 2)='zooplankton carbon assimilation efficiency'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
#  endif

      Vinfo( 1)='ZooAE_N'
      Vinfo( 2)='zooplankton mitrogen assimilation efficiency'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooBM'
      Vinfo( 2)='rate for zooplankton basal metabolism'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooCN'
      Vinfo( 2)='zooplankton Carbon:Nitrogen ratio'
      Vinfo( 3)='mole_C mole_N-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooER'
      Vinfo( 2)='zooplankton specific excretion rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooGGE_C'
      Vinfo( 2)='zooplankton carbon gross growth efficiency'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooGR'
      Vinfo( 2)='zooplankton maximum growth rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooMin'
      Vinfo( 2)='minimum zooplankton threshold'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='ZooMR'
      Vinfo( 2)='zooplankton mortality rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='LDeRR'
      Vinfo( 2)='rate of large detritus re-mineralization to NO3'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='CoagR'
      Vinfo( 2)='coagulation rate'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='SDeAR'
      Vinfo( 2)='rate of small detritus aggregation to large detritus'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='SDeBR'
      Vinfo( 2)='rate of small detritus breakdown to NH4'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='SDeRR'
      Vinfo( 2)='remineralization rate for small detritus'
      Vinfo( 3)='day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='wPhy'
      Vinfo( 2)='vertical sinking velocity for phytoplankton'
      Vinfo( 3)='meter day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

#  ifdef CARBON
      Vinfo( 1)='wLDeC'
      Vinfo( 2)='vertical sinking velocity for large carbon detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
#  endif

      Vinfo( 1)='wLDeN'
      Vinfo( 2)='vertical sinking velocity for large nitrogen detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

#  ifdef CARBON
      Vinfo( 1)='wSDeC'
      Vinfo( 2)='vertical sinking velocity for small carbon detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
#  endif

      Vinfo( 1)='wSDeN'
      Vinfo( 2)='vertical sinking velocity for small nitrogen detritus'
      Vinfo( 3)='meter day-1'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

# endif
#endif
!
!-----------------------------------------------------------------------
!  Define grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch: Spherical or Cartesian.
!
      Vinfo( 1)='spherical'
      Vinfo( 2)='grid type logical switch'
      Vinfo( 7)='spherical'
      Vinfo( 8)='Cartesian'
      status=def_var(ncid,varid,nf_char,0,0,Aval,Vinfo,ncname)
!
!  Domain Length.
!
      Vinfo( 1)='xl'
      Vinfo( 2)='domain length in the XI-direction'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='el'
      Vinfo( 2)='domain length in the ETA-direction'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
#ifdef SOLVE3D
!
!  S-coordinate parameters.
!
      Vinfo( 1)='theta_s'
      Vinfo( 2)='S-coordinate surface control parameter'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='theta_b'
      Vinfo( 2)='S-coordinate bottom control parameter'
      Vinfo( 3)='nondimensional'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='Tcline'
      Vinfo( 2)='S-coordinate surface/bottom layer width'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)

      Vinfo( 1)='hc'
      Vinfo( 2)='S-coordinate parameter, critical depth'
      Vinfo( 3)='meter'
      status=def_var(ncid,varid,NF_TYPE,0,0,Aval,Vinfo,ncname)
!
!  S-coordinate non-dimensional independent variable at RHO-points.
!
      Vinfo( 1)='sc_r'
      Vinfo( 2)='S-coordinate at RHO-points'
      Vinfo( 3)='nondimensional'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='sc_r, scalar'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ncid,varid,NF_TYPE,1,srdim,Aval,Vinfo,ncname)
!
!  S-coordinate non-dimensional independent variable at W-points.
!
      Vinfo( 1)='sc_w'
      Vinfo( 2)='S-coordinate at W-points'
      Vinfo( 3)='nondimensional'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='sc_w, scalar'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ncid,varid,NF_TYPE,1,swdim,Aval,Vinfo,ncname)
!
!  S-coordinate non-dimensional stretching curves at RHO-points.
!
      Vinfo( 1)='Cs_r'
      Vinfo( 2)='S-coordinate stretching curves at RHO-points'
      Vinfo( 3)='nondimensional'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='Cs_r, scalar'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ncid,varid,NF_TYPE,1,srdim,Aval,Vinfo,ncname)
!
!  S-coordinate non-dimensional stretching curves at W-points.
!
      Vinfo( 1)='Cs_w'
      Vinfo( 2)='S-coordinate stretching curves at W-points'
      Vinfo( 3)='nondimensional'
      Vinfo( 5)='valid_min'
      Vinfo( 6)='valid_max'
      Vinfo(14)='Cs_w, scalar'
      Aval(2)=-1.0_r8
      Aval(3)=0.0_r8
      status=def_var(ncid,varid,NF_TYPE,1,swdim,Aval,Vinfo,ncname)
#endif
!
!  User generic parameters.
!
      IF (Nuser.gt.0) THEN
        Vinfo( 1)='user'
        Vinfo( 2)='user generic parameters'
        Vinfo(14)='user, scalar'
        status=def_var(ncid,varid,NF_TYPE,1,usrdim,Aval,Vinfo,ncname)
      END IF
#ifdef STATIONS
!
!  Station positions.
!
      IF (ncid.eq.ncstaid(ng)) THEN
        Vinfo( 1)='Ipos'
        Vinfo( 2)='stations I-direction positions'
        status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,          &
     &                 ncname)

        Vinfo( 1)='Jpos'
        Vinfo( 2)='stations J-direction positions'
        status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,          &
     &                 ncname)
      END IF
#endif
#ifndef NO_WRITE_GRID
!
!  Bathymetry.
!
      IF (ncid.ne.ncfltid(ng)) THEN
        Vinfo( 1)='h'
        Vinfo( 2)='bathymetry at RHO-points'
        Vinfo( 3)='meter'
        Vinfo(14)='bath, scalar'
        IF (ncid.eq.ncstaid(ng)) THEN
          status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,        &
     &                   ncname)
        ELSE
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)
        END IF
!
!  Coriolis Parameter.
!
        IF (ncid.ne.ncstaid(ng)) THEN
          Vinfo( 1)='f'
          Vinfo( 2)='Coriolis parameter at RHO-points'
          Vinfo( 3)='second-1'
          Vinfo(14)='coriolis, scalar'
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)
        END IF
!
!  Curvilinear coordinate metrics.
!
        IF (ncid.ne.ncstaid(ng)) THEN
          Vinfo( 1)='pm'
          Vinfo( 2)='curvilinear coordinate metric in XI'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pm, scalar'
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)

          Vinfo( 1)='pn'
          Vinfo( 2)='curvilinear coordinate metric in ETA'
          Vinfo( 3)='meter-1'
          Vinfo(14)='pn, scalar'
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)
        END IF
!
!  Longitude and latitude coordinates of RHO-points, if appropriate.
!
        IF (spherical) THEN
          Vinfo( 1)='lon_rho'
          Vinfo( 2)='longitude of RHO-points'
          Vinfo( 3)='degree_east'
          Vinfo(14)='lon_rho, scalar'
          IF (ncid.eq.ncstaid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,      &
     &                     ncname)
          ELSE
            status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF

          Vinfo( 1)='lat_rho'
          Vinfo( 2)='latitude of RHO-points'
          Vinfo( 3)='degree_north'
          Vinfo(14)='lat_rho, scalar'
          IF (ncid.eq.ncstaid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,      &
     &                     ncname)
          ELSE
            status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
        ELSE
          Vinfo( 1)='x_rho'
          Vinfo( 2)='x-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='x_rho, scalar'
          IF (ncid.eq.ncstaid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,      &
     &                     ncname)
          ELSE
            status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
          Vinfo( 1)='y_rho'
          Vinfo( 2)='y-locations of RHO-points'
          Vinfo( 3)='meter'
          Vinfo(14)='y_rho, scalar'
          IF (ncid.eq.ncstaid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,      &
     &                     ncname)
          ELSE
            status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
        END IF
# ifdef CURVGRID
!
!  Angle between XI-axis and EAST at RHO-points, if appropriate.
!
        IF (spherical) THEN
          Vinfo( 1)='angle'
          Vinfo( 2)='angle between XI-axis and EAST'
          Vinfo( 3)='radians'
          Vinfo(14)='angle, scalar'
          IF (ncid.eq.ncstaid(ng)) THEN
            status=def_var(ncid,varid,NF_TYPE,1,stadim,Aval,Vinfo,      &
     &                     ncname)
          ELSE
            status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,      &
     &                     ncname)
          END IF
        END IF
# endif
# ifdef MASKING
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (ncid.ne.ncstaid(ng)) THEN
          Vinfo( 1)='mask_rho'
          Vinfo( 2)='mask on RHO-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          status=def_var(ncid,varid,NF_TYPE,2,t2dgrd,Aval,Vinfo,        &
     &                   ncname)

          Vinfo( 1)='mask_u'
          Vinfo( 2)='mask on U-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          status=def_var(ncid,varid,NF_TYPE,2,u2dgrd,Aval,Vinfo,        &
     &                   ncname)

          Vinfo( 1)='mask_v'
          Vinfo( 2)='mask on V-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          status=def_var(ncid,varid,NF_TYPE,2,v2dgrd,Aval,Vinfo,        &
     &                   ncname)

          Vinfo( 1)='mask_psi'
          Vinfo( 2)='mask on psi-points'
          Vinfo( 9)='land'
          Vinfo(10)='water'
          status=def_var(ncid,varid,NF_TYPE,2,p2dgrd,Aval,Vinfo,        &
     &                   ncname)
        END IF
# endif
      END IF
#endif
      def_info=status
      RETURN
      end function def_info
