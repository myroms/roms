#include "cppdefs.h"
      MODULE mod_netcdf
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
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
        integer :: NF_FOUT = nf_double
#else
        integer :: NF_FOUT = nf_real
#endif
#ifdef RST_SINGLE
        integer :: NF_FRST = nf_real
#else
        integer :: NF_FRST = nf_double
#endif
#ifdef DOUBLE_PRECISION
        integer :: NF_TYPE = nf_double
#else
        integer :: NF_TYPE = nf_real
#endif

      END MODULE mod_netcdf
