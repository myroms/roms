#include "cppdefs.h"

#if !defined TS_FIXED && defined SOLVE3D &&  \
    (defined TS_DIF2  || defined TS_DIF4)
/*
**======================================================================
** Copyright (c) 2004 ROMS/TOMS Group                                  !
**================================================= Hernan G. Arango ===
**                                                                     !
** This routine computes horizontal mixing of tracers.                 !
**                                                                     !
**======================================================================
*/

# ifdef TS_DIF2
#  if defined MIX_S_TS
#   include "t3dmix2_s.h"
#  elif defined MIX_GEO_TS
#   include "t3dmix2_geo.h"
#  elif defined MIX_ISO_TS
#   include "t3dmix2_iso.h"
#  endif
# endif

# ifdef TS_DIF4
#  if defined MIX_S_TS
#   include "t3dmix4_s.h"
#  elif defined MIX_GEO_TS
#   include "t3dmix4_geo.h"
#  elif defined MIX_ISO_TS
#   include "t3dmix4_iso.h"
#  endif
# endif

#endif
