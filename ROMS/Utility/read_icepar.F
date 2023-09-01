#include "cppdefs.h"
#ifdef SEAICE
!!
!!git $Id$
!!svn $Id: read_icepar.F 1176 2023-07-01 19:23:18Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2023 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!!  This routine reads and reports seaice model input parameters.      !
!!                                                                     !
!!======================================================================
!!
# if defined ICE_MODEL
#  include <ice_inp.h>
# endif
#else
      SUBROUTINE read_IcePar
      END SUBROUTINE read_IcePar
#endif
