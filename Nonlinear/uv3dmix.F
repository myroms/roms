#include "cppdefs.h"

#if defined SOLVE3D && (defined UV_VIS2 || defined UV_VIS4)

/*
**======================================================================
** Copyright (c) 2004 ROMS/TOMS Group                                  !
**================================================= Hernan G. Arango ===
**                                                                     !
** This routine computes horizontal viscosity of momentum.             !
**                                                                     !
**======================================================================
*/

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

