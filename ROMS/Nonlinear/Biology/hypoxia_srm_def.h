/*
** git $Id: 560fb376ff8a4576170ebcd4b459de6bcce908f6 $
** svn $Id: hypoxia_srm_def.h 937 2019-01-28 06:13:04Z arango $
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2019 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.txt                                             **
************************************************************************
**                                                                    **
**  Defines Hypoxia Simple Respiration Model input parameters in      **
**  output NetCDF files. It is included in routine  "def_info.F".     **
**                                                                    **
************************************************************************
*/

!
!  Define Hypoxia Simple Respiration Model parameters.
!
      Vinfo( 1)='BioIter'
      Vinfo( 2)='number of iterations to achieve convergence'
      status=def_var(ng, model, ncid, varid, nf90_int,                  &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN

      Vinfo( 1)='ResRate'
      Vinfo( 2)='total biological respiration rate'
      Vinfo( 3)='day-1'
      status=def_var(ng, model, ncid, varid, NF_TYPE,                   &
     &               1, (/0/), Aval, Vinfo, ncname,                     &
     &               SetParAccess = .FALSE.)
      IF (FoundError(exit_flag, NoError, __LINE__,                      &
     &               __FILE__)) RETURN
