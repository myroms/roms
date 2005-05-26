#include "cppdefs.h"
!!
!!================================================== Andrew M. Moore ===
!!  Copyright (c) 2005 ROMS/TOMS Adjoint Group                         !
!!================================================= Hernan G. Arango ===
!!                                                                     !
!!  Regional Ocean Model System (ROMS) Drivers:                        !
!!                                                                     !
!!  nl_ocean.F        ROMS/TOMS nonlinear model.                       !
!!  tl_ocean.F        ROMS/TOMS tangent linear model.                  !
!!  ad_ocean.F        ROMS/TOMS adjoint model.                         !
!!  grad_ocean.F      ROMS/TOMS tangent linear and adjoint models      !
!!                      gradient test.                                 !
!!  tlcheck_ocean.F   ROMS/TOMS tangent linear model linearization     !
!!                      test.                                          !
!!  pert_ocean.F      ROMS/TOMS tangent linear and adjoint models      !
!!                      sanity test.                                   !
!!  picard_ocean.F    ROMS/TOMS Picard test for representers tangent   !
!!                      linear model.                                  ! 
!!  op_ocean.F        ROMS/TOMS optimal perturbations.                 !
!!  s4dvar_ocean.F    ROMS/TOMS strong constraint 4DVAR assimilation.  !
!!  is4dvar_ocean.F   ROMS/TOMS Strong constraint incremental 4DVAR    !
!!                      assimilation.                                  !
!!                                                                     !
!!======================================================================
!!
#if defined OPT_PERTURBATION
# include "op_ocean.h"
#elif defined GRADIENT_CHECK
# include "grad_ocean.h"
#elif defined TLM_CHECK
# include "tlcheck_ocean.h"
#elif defined INNER_PRODUCT || defined SANITY_CHECK
# include "pert_ocean.h"
#elif defined PICARD_TEST
# include "picard_ocean.h"
#elif defined S4DVAR
# include "s4dvar_ocean.h"
#elif defined IS4DVAR
# include "is4dvar_ocean.h"
#else
# include "nl_ocean.h"
#endif
