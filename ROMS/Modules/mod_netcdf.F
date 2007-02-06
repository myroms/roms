#include "cppdefs.h"
      MODULE mod_netcdf
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains all NetCDF variables definitions.              !
!                                                                      !
!  NF_FOUT    Output real array variables type.                        !
!  NF_FRST    Output real array variables type in restart file.        !
!  NF_TYPE    Output real scalar variables type.                       !
!                                                                      !
!=======================================================================
!
        USE mod_kinds

        implicit none

#include "netcdf.inc"

#ifdef OUT_DOUBLE
        integer, parameter :: NF_FOUT = nf_double
#else
        integer, parameter :: NF_FOUT = nf_real
#endif
#ifdef RST_SINGLE
        integer, parameter :: NF_FRST = nf_real
#else
        integer, parameter :: NF_FRST = nf_double
#endif
#ifdef DOUBLE_PRECISION
        integer, parameter :: NF_TYPE = nf_double
#else
        integer, parameter :: NF_TYPE = nf_real
#endif

      END MODULE mod_netcdf
