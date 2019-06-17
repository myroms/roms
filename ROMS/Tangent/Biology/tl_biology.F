#include "cppdefs.h"
      MODULE tl_biology_mod
#if defined TANGENT && defined BIOLOGY
!
!git $Id: 560fb376ff8a4576170ebcd4b459de6bcce908f6 $
!svn $Id: tl_biology.F 937 2019-01-28 06:13:04Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2019 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine computes the tangent linear source and sink terms for  !
!  the selected biology model.                                         !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: tl_biology

      CONTAINS

# if defined BIO_FENNEL
!!#  include <tl_fennel.h>
# elif defined NEMURO
!!#  include <tl_nemuro.h>
# elif defined NPZD_FRANKS
#  include <tl_npzd_Franks.h>
# elif defined NPZD_IRON
#  include <tl_npzd_iron.h>
# elif defined NPZD_POWELL
#  include <tl_npzd_Powell.h>
# endif

#endif
      END MODULE tl_biology_mod
