#include "cppdefs.h"
#ifdef SEDIMENT
!!
!!git $Id$
!!svn $Id: read_sedpar.F 937 2019-01-28 06:13:04Z arango $
!!================================================= Hernan G. Arango ===
!! Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
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
