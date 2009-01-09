#include "cppdefs.h"
      SUBROUTINE checkadj
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This subroutine checks activated C-preprocessing options for        !
!  consistency with all available algorithms.                          !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE mod_strings
!
      USE strings_mod, ONLY : uppercase
!
      implicit none
!
!  Local variable declarations.
!
      integer :: ic = 0
      integer :: ifound

      character (len=40) :: string
!
!-----------------------------------------------------------------------
!  Report issues with various C-preprocessing options.
!-----------------------------------------------------------------------
!
      string=uppercase('ts_smagorinsky')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string), 'stability problems'
      END IF

      string=uppercase('ts_u3adv_split')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string), 'stability problems'
      END IF

      string=uppercase('uv_smagorinsky')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string), 'stability problems'
      END IF

      string=uppercase('uv_u3adv_split')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        IF (Master) WRITE(stdout,10) TRIM(string), 'stability problems'
      END IF

#if defined TANGENT || defined TL_IOMS || defined ADJOINT
!
!-----------------------------------------------------------------------
!  Stop if unsupported C-preprocessing options are activated for the
!  adjoint-based algorithms.
!-----------------------------------------------------------------------
!
      string=uppercase('assimilation_ssh')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'conflict'
      END IF

      string=uppercase('assimilation_sst')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'conflict'
      END IF

      string=uppercase('assimilation_t')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'conflict'
      END IF

      string=uppercase('assimilation_uv')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'conflict'
      END IF

      string=uppercase('assimilation_uvsur')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'conflict'
      END IF

      string=uppercase('bedload_mpm')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('bedload_soulsby')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('bio_fasham')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('bulk_fluxes')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
!!      ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not working'
      END IF

      string=uppercase('bvf_mixing')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('east_fsradiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('east_m2radiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('east_m3radiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('east_tradiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('ecosim')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('gls_mixing')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
!!      ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string),  'not differentiable'
      END IF

      string=uppercase('lmd_mixing')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
!!      ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not differentiable'
      END IF

      string=uppercase('mb_bbl')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('my25_mixing')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
!!      ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not differentiable'
      END IF

      string=uppercase('nearshore_mellor')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('nemuro')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('north_fsradiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('north_m2radiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('north_m3radiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('north_tradiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('npzd_franks')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
!!      ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not working'
      END IF

      string=uppercase('pj_gradpq2')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('pj_gradpq4')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('refdif_coupling')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not allowed'
      END IF

      string=uppercase('sediment')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('sed_dens')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not tested'
      END IF

      string=uppercase('sed_morph')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('sg_bbl')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('south_fsradiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('south_m2radiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('south_m3radiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('south_tradiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('ssw_bbl')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('suspload')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('swan_coupling')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not allowed'
      END IF

      string=uppercase('ts_mpdata')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('ts_smagorinsky')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('ts_u3adv_split')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('uv_smagorinsky')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('uv_u3adv_split')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF

      string=uppercase('west_fsradiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('west_m2radiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('west_m3radiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('west_tradiation')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not finished'
      END IF

      string=uppercase('wet_dry')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20)  TRIM(string), 'not coded'
      END IF

      string=uppercase('wrf_coupling')
      ifound=INDEX(TRIM(Coptions), TRIM(string))
      IF (ifound.ne.0) THEN
        ic=ic+1
        IF (Master) WRITE(stdout,20) TRIM(string), 'not coded'
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Set execution error flag to stop execution.
!-----------------------------------------------------------------------
!
      IF (ic.gt.0) THEN
        exit_flag=5
      END IF
!
 10   FORMAT (/,' CHECKADJ - use caution when activating: ', a,/,12x,   &
     &        'REASON: ',a,'.')
#if defined TANGENT || defined TL_IOMS || defined ADJOINT
 20   FORMAT (/,' CHECKADJ - unsupported option in adjoint-based',      &
     &        ' algorithms: ',a,/,12x,'REASON: ',a,'.')
#endif

      RETURN
      END SUBROUTINE checkadj
