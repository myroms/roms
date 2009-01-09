#include "cppdefs.h"

      MODULE cgradient_mod
!!
!!svn $Id$
!!================================================== Hernan G. Arango ==
!! Copyright (c) 2002-2009 The ROMS/TOMS Group       Andrew M. Moore   !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!======================================================================
!!                                                                     !
!! This module minimizes  incremental  4Dvar  quadratic cost function  !
!! using a preconditioned version of the conjugate gradient algorithm  !
!! proposed by Mike Fisher (ECMWF) and modified  by  Tshimanga et al.  !
!! (2008) using Limited-Memory Precondtioners (LMP).                   !
!!                                                                     !
!! References:                                                         !
!!                                                                     !
!!   Fisher, M., 1997: Efficient Minimization of Quadratic Penalty     !
!!     funtions, unpublish manuscript, 1-14.                           !
!!                                                                     !
!!   Fisher, M., 1998: Minimization algorithms for variational data    !
!!     data assimilation. 'In Recent Developments in Numerical         !
!!     Methods for Atmospheric Modelling', pp 364-385, ECMWF.          !
!!                                                                     !
!!   Tchimanga, J., S. Gratton, A.T. Weaver, and A. Sartenaer, 2008:   !
!!     Limited-memory preconditioners, with application to incremental !
!!     four-dimensional variational ocean data assimilation, Q.J.R.    !
!!     Meteorol. Soc., 134, 753-771.                                   !
!!                                                                     !
!!======================================================================
!!

#if defined IS4DVAR
#  include "cgradient_lanczos.h"
#endif

      END MODULE cgradient_mod
