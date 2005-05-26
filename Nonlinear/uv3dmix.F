#include "cppdefs.h"

      MODULE uv3dmix_mod

#if defined SOLVE3D && defined NONLINEAR && \
   (defined UV_VIS2 || defined UV_VIS4)

!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine computes horizontal viscosity of momentum.             !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
# ifdef UV_VIS2
      PUBLIC uv3dmix2
# endif
# ifdef UV_VIS4
      PUBLIC uv3dmix4
# endif

      CONTAINS

# ifdef UV_VIS2
#  if defined MIX_S_UV
#   include "uv3dmix2_s.h"
#  elif defined MIX_GEO_UV
#   include "uv3dmix2_geo.h"
#  endif
# endif

# ifdef UV_VIS4
#  if defined MIX_S_UV
#   include "uv3dmix4_s.h"
#  elif defined MIX_GEO_UV
#   include "uv3dmix4_geo.h"
#  endif
# endif

#endif

      END MODULE uv3dmix_mod
