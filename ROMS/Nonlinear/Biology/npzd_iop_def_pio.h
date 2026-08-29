/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2024 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Defines IOP-based, NPZD ecosystem model input parameters in       **
**  output NetCDF files. It is included in routine  "def_info.F".     **
**                                                                    **
************************************************************************
*/

!
!  Define IOP-based, NPZD ecosystem model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ng, model, pioFile, pioVar, PIO_int,               &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='PARfrac'
      Vinfo( 2)='photosynthetically available radiation fraction'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='AttSW'
      Vinfo( 2)='light attenuation due to sea water'
      Vinfo( 3)='meter-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='CDOM_LightAtt'
      Vinfo( 2)='light attenuation due to CDOM'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='CDOM_sigma'
      Vinfo( 2)='light-dependent degradation rate for CDOM'
      Vinfo( 3)='watt-1 m2 day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='K_DIN'
      Vinfo( 2)='Half-saturation for phytoplankton nitrogen uptake'
      Vinfo( 3)='millimole_N meter-3'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='PhyIS'
      Vinfo( 2)='phytoplankton initial slope of the P-I cureve'
      Vinfo( 3)='meter2 watt-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='PhyMRD'
      Vinfo( 2)='phytoplankton mortality rate to the detritus pool'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='PhyMRN'
      Vinfo( 2)='phytoplankton mortality rate to the Nitrogen pool'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='Vm_DIN'
      Vinfo( 2)='nitrogen uptake rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='DetRR'
      Vinfo( 2)='detritus remineralization rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='ThetaM'
      Vinfo( 2)='Maximum ratio of phytoplankton backscatter to '//      &
     &          'absorption'
      Vinfo( 3)='nondimensional'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='BphyMap'
      Vinfo( 2)='Mapping from phytoplankton backscatter to '//          &
     &          'Nitrogen biomass at 440 nm'
      Vinfo( 3)='millimole_N meter-2'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='BdetMap'
      Vinfo( 2)='Mapping from detritus backscatter to Nitrogen '//      &
     &           'biomass at 440 nm'
      Vinfo( 3)='millimole_N meter-2'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='wPhy'
      Vinfo( 2)='phytoplankton sinking rate'
      Vinfo( 3)='m day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN

      Vinfo( 1)='wDet'
      Vinfo( 2)='detrital sinking rate'
      Vinfo( 3)='m day-1'
      status=def_var(ng, model, pioFile, pioVar, PIO_TYPE,              &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
