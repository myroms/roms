#include "cppdefs.h"
      MODULE tracer_index_mod
!
!git $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2025 The ROMS Group                              !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.md                                               !
!=======================================================================
!                                                                      !
!  Given a tracer variable name, this function returns its associated  !
!  array index,                                                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
      USE strings_mod, ONLY : FoundError
!
      implicit none
!
      PUBLIC  :: tracer_index
      PRIVATE
!
      CONTAINS
!
      FUNCTION tracer_index (name) RESULT (TrcIndex)
!
!=======================================================================
!                                                                      !
!  Gets state tracer array index from its variable name.               !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     name         Tracer variable name (string)                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     TrcIndex     Tracer varible index in state array (integer)       !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      character (len=*), intent(in) :: name
!
!  Local variable declarations.
!
      integer :: TrcIndex
!
      character (len=*),  parameter :: MyFile =                         &
     &    __FILE__
!
      SourceFile=MyFile
!
!-----------------------------------------------------------------------
!  Set state tracer array index.
!-----------------------------------------------------------------------
!
      SELECT CASE (TRIM(name))
        CASE ('SST', 'temp',                                            &
     &        'sea_water_temperature',                                  &
     &        'sea_water_potential_temperature',                        &
     &        'waterTemperature')
          TrcIndex=itemp
        CASE ('SSS', 'salt',                                            &
     &        'sea_water_salinity',                                     &
     &        'salinity')
          TrcIndex=isalt
        CASE DEFAULT
          IF (Master) WRITE (stdout,10) TRIM(name)
          exit_flag=5
          IF (FoundError(exit_flag, NoError, __LINE__, MyFile)) RETURN
      END SELECT
!
 10   FORMAT (/,' TRACER_INDEX - Cannot find index for: ',a)
      RETURN
      END FUNCTION tracer_index
!
      END MODULE tracer_index_mod
