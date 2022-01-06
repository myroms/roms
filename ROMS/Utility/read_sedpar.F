#include "cppdefs.h"
#ifdef SEDIMENT
!!
!!git $Id$
!!svn $Id: read_sedpar.F 1099 2022-01-06 21:01:01Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2022 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!!                                                                     !
!!  This routine reads and reports sediment model input parameters.    !
!!                                                                     !
!!======================================================================
!!
# include <sediment_inp.h>
#else
      SUBROUTINE read_SedPar
      END SUBROUTINE read_SedPar
#endif
