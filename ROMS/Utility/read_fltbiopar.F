#include "cppdefs.h"
#if defined FLOATS && defined FLOAT_BIOLOGY
!!
!!git $Id$
!!svn $Id: read_fltbiopar.F 1210 2024-01-03 22:03:03Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2024 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.md                                               !
!!======================================================================
!!                                                                     !
!!  This routine reads and reports biological floats model input       !
!!  parameters.                                                        !
!!                                                                     !
!!======================================================================
!!
# if defined FLOAT_OYSTER
#  include <oyster_floats_inp.h>
# endif
#else
      SUBROUTINE read_FltBioPar
      END SUBROUTINE read_FltBioPar
#endif
