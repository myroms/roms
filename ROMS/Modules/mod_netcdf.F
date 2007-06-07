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
        USE netcdf

        implicit none

#ifdef OUT_DOUBLE
        integer, parameter :: NF_FOUT = nf90_double
#else
        integer, parameter :: NF_FOUT = nf90_real
#endif
#ifdef RST_SINGLE
        integer, parameter :: NF_FRST = nf90_real
#else
        integer, parameter :: NF_FRST = nf90_double
#endif
#ifdef DOUBLE_PRECISION
        integer, parameter :: NF_TYPE = nf90_double
#else
        integer, parameter :: NF_TYPE = nf90_real
#endif

      END MODULE mod_netcdf
