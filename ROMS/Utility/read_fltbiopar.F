#include "cppdefs.h"
#if defined FLOATS && defined FLOAT_BIOLOGY
!!
!!git $Id$
!!svn $Id: read_fltbiopar.F 1099 2022-01-06 21:01:01Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
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
