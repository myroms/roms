#include "cppdefs.h"
#ifdef SEAICE
!!
!!git $Id$
!!svn $Id: read_icepar.F 1210 2024-01-03 22:03:03Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
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
