#include "cppdefs.h"
      FUNCTION wrt_info (ng, ncid, wrtThread, ncname)
!
!=======================================================================
!  Copyright (c) 2002 ROMS/TOMS Group                                  !
!================================================== Hernan G. Arango ===
!                                                                      !
!  This routine writes out information variables into requested        !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng       Nested grid number (integer).                           !
!     ncid     NetCDF file ID (integer).                               !
!     ncname   NetCDF file name (character).                           !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     wrt_info   Error flag (integer).                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
#ifdef BIOLOGY
      USE mod_biology
#endif
      USE mod_grid
      Use mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      USE extract_mod, ONLY : extract2d
!
      implicit none
!
!  Imported variable declarations.
!
      logical, intent(in) :: wrtThread

      integer, intent(in) :: ng, ncid

      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: Cgrid = .true.

      integer :: LBi, UBi, LBj, UBj
      integer :: i, j, k, ibry, ilev, itrc, varid

      integer, dimension(2) :: start, total

      integer :: nf_fwrite, wrt_info

      real(r8) :: scale
#ifdef SOLVE3D
# ifdef TS_DIF4
      real(r8), dimension(NT(ng)) :: diff
# endif
      real(r8), dimension(NT(ng)) :: nudg
      real(r8), dimension(NT(ng),4) :: Tobc
#endif
#ifdef STATIONS
      real(r8), dimension(Nstation(ng)) :: Zpos, wrk
#endif

      character (len=1) :: char1
!
      LBi=lbound(GRID(ng)%h,DIM=1)
      UBi=ubound(GRID(ng)%h,DIM=1)
      LBj=lbound(GRID(ng)%h,DIM=2)
      UBj=ubound(GRID(ng)%h,DIM=2)
!
!-----------------------------------------------------------------------
!  Write out running parameters.
!-----------------------------------------------------------------------
!
      IF (wrtThread) THEN
!
!  Time stepping parameters.
!
        wrt_info=nf_inq_varid(ncid,'ntimes',varid)
        wrt_info=nf_put_var1_int(ncid,varid,1,ntimes)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ntimes', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'ndtfast',varid)
        wrt_info=nf_put_var1_int(ncid,varid,1,ndtfast)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ndtfast', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'dt',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,dt)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'dt', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'dtfast',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,dtfast)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'dtfast', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'dstart',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,dstart)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'dstart', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'nhis',varid)
        wrt_info=nf_put_var1_int(ncid,varid,1,nhis)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'nhis', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        IF (ndefhis(ng).gt.0) THEN
          wrt_info=nf_inq_varid(ncid,'ndefhis',varid)
          wrt_info=nf_put_var1_int(ncid,varid,1,ndefhis(ng))
          IF (wrt_info.ne.nf_noerr) THEN
            WRITE (stdout,10) 'ndefhis', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF

        wrt_info=nf_inq_varid(ncid,'nrst',varid)
        wrt_info=nf_put_var1_int(ncid,varid,1,nrst)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'nrst', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

#ifdef AVERAGES
        wrt_info=nf_inq_varid(ncid,'ntsavg',varid)
        wrt_info=nf_put_var1_int(ncid,varid,1,ntsavg)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ntsavg', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'navg',varid)
        wrt_info=nf_put_var1_int(ncid,varid,1,navg)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'navg', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        IF (ndefavg(ng).gt.0) THEN
          wrt_info=nf_inq_varid(ncid,'ndefavg',varid)
          wrt_info=nf_put_var1_int(ncid,varid,1,ndefavg(ng))
          IF (wrt_info.ne.nf_noerr) THEN
            WRITE (stdout,10) 'ndefavg', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF
#endif

#ifdef STATIONS
        wrt_info=nf_inq_varid(ncid,'nsta',varid)
        wrt_info=nf_put_var1_int(ncid,varid,1,nsta)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'nsta', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif

#if defined POWER_LAW && defined SOLVE3D
!
!  Power-law shape filter parameters for time-averaging of barotropic
!  fields.
!
        wrt_info=nf_inq_varid(ncid,'Falpha',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Falpha)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Falpha', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Fbeta',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Fbeta)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Fbeta', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Fgamma',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Fgamma)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Fgamma', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
!
!  Horizontal mixing coefficients.
!
#if defined SOLVE3D && defined TS_DIF2
        wrt_info=nf_inq_varid(ncid,'tnu2',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,NT(ng),tnu2)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'tnu2', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif

#if defined SOLVE3D && defined TS_DIF4
        wrt_info=nf_inq_varid(ncid,'tnu4',varid)
        DO itrc=1,NT(ng)
          diff(itrc)=tnu4(itrc,ng)*tnu4(itrc,ng)
        END DO
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,NT(ng),diff)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'tnu4', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif

#ifdef UV_VIS2
        wrt_info=nf_inq_varid(ncid,'visc2',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,visc2)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'visc2', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif

#ifdef UV_VIS4
        wrt_info=nf_inq_varid(ncid,'visc4',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,visc4*visc4)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'visc4', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif

#if defined SOLVE3D && (defined MY25_MIXING || defined GLS_MIXING)
# ifdef TKE_DIF2
        wrt_info=nf_inq_varid(ncid,'tkenu2',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,tkenu2(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'tkenu2', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif

# ifdef TKE_DIF4
        wrt_info=nf_inq_varid(ncid,'tkenu4',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,tkenu4(ng)*tkenu4(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'tkenu4', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
#ifdef SOLVE3D
!
!  Background vertical mixing coefficients.
!
        wrt_info=nf_inq_varid(ncid,'Akt_bak',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,NT(ng),Akt_bak(1,ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Akt_bak', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Akv_bak',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Akv_bak(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Akv_bak', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

# if defined MY25_MIXING || defined GLS_MIXING
        wrt_info=nf_inq_varid(ncid,'Akk_bak',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Akk_bak(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Akk_bak', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Akp_bak',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Akp_bak(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Akp_bak', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
!
!  Drag coefficients.
!
        wrt_info=nf_inq_varid(ncid,'rdrg',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,rdrg(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'rdrg', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'rdrg2',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,rdrg2(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'rdrg2', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

#ifdef SOLVE3D
        wrt_info=nf_inq_varid(ncid,'Zob',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Zob(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Zob', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Zos',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Zos(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Zos', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#if defined SOLVE3D && defined GLS_MIXING
!
!  Generic length-scale parameters.
!
        wrt_info=nf_inq_varid(ncid,'gls_p',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_p(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_p', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_m',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_m(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_m', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_n',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_n(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_n', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_cmu0',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_cmu0(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_cmu0', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_c1',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_c1(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_c1', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_c2',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_c2(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_c2', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_c3m',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_c3m(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_c3m', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_c3p',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_c3p(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_c3p', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_sigk',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_sigk(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_sigk', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_sigp',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_sigp(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_sigp', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_Kmin',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_Kmin(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_Kmin', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'gls_Pmin',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gls_Pmin(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gls_Pmin', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
!
!  Nudging inverse time scales used in various tasks.
!
        wrt_info=nf_inq_varid(ncid,'Znudg',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Znudg(ng)/sec2day)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Znudg', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'M2nudg',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,M2nudg(ng)/sec2day)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'M2nudg', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

#ifdef SOLVE3D
        wrt_info=nf_inq_varid(ncid,'M3nudg',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,M3nudg(ng)/sec2day)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'M3nudg', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Tnudg',varid)
        DO itrc=1,NT(ng)
          nudg(itrc)=Tnudg(itrc,ng)/sec2day
        END DO
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,NT(ng),nudg)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Tnudg', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef NUDGING
!
!  Nudging inverse time scales used in data assimilation.
!
        wrt_info=nf_inq_varid(ncid,'Znudass',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Znudass(ng)/sec2day)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Znudass', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'M2nudass',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,M2nudass(ng)/sec2day)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'M2nudass', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

# ifdef SOLVE3D
        wrt_info=nf_inq_varid(ncid,'M3nudass',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,M3nudass(ng)/sec2day)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'M3nudass', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Tnudass',varid)
        DO itrc=1,NT(ng)
          nudg(itrc)=Tnudass(itrc,ng)/sec2day
        END DO
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,NT(ng),nudg)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Tnudass', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
#ifdef NUDGING_COFF
!
!  Open boundary nudging, inverse time scales.
!
        wrt_info=nf_inq_varid(ncid,'FSobc_in',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,4,FSobc_in(ng,1))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'FSobc_in', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'FSobc_out',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,4,FSobc_out(ng,1))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'FSobc_out', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'M2obc_in',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,4,M2obc_in(ng,1))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'M2obc_in', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'M2obc_out',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,4,M2obc_out(ng,1))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'M2obc_out', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

# ifdef SOLVE3D
        start(1)=1
        total(1)=NT(ng)
        start(2)=1
        total(2)=4

        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_in(itrc,ng,ibry)
          END DO
        END DO
        wrt_info=nf_inq_varid(ncid,'Tobc_in',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,start,total,Tobc)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Tobc_in', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        DO ibry=1,4
          DO itrc=1,NT(ng)
            Tobc(itrc,ibry)=Tobc_in(itrc,ng,ibry)
          END DO
        END DO
        wrt_info=nf_inq_varid(ncid,'Tobc_out',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,start,total,Tobc)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Tobc_out', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'M3obc_in',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,4,M3obc_in(ng,1))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'M3obc_in', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'M3obc_out',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,4,M3obc_out(ng,1))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'M3obc_out', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
!
!  Equation of State parameters.
!
        wrt_info=nf_inq_varid(ncid,'rho0',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,rho0)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'rho0', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

#if defined SOLVE3D && !defined NONLIN_EOS
        wrt_info=nf_inq_varid(ncid,'R0',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,R0(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'R0', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Tcoef',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Tcoef(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Tcoef', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Scoef',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Scoef(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Scoef', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
#ifdef SOLVE3D
!
!  Various parameters.
!
# ifdef BODYFORCE
        wrt_info=nf_inq_varid(ncid,'levsfrc',varid)
        wrt_info=nf_put_var1_int(ncid,varid,1,levsfrc(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'levsfrc', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'levbfrc',varid)
        wrt_info=nf_put_var1_int(ncid,varid,1,levbfrc(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'levbfrc', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
# endif
#endif
!
!  Slipperiness parameters.
!
        wrt_info=nf_inq_varid(ncid,'gamma2',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,gamma2)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'gamma2', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

#if defined BIOLOGY && defined SOLVE3D
# ifdef BIO_FASHAM
!
!  Write out Fasham type biological module parameters.
!
        wrt_info=nf_inq_varid(ncid,'BioIter',varid)
        wrt_info=nf_put_var1_int(ncid,varid,1,BioIter(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'BioIter', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'AttSW',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,AttSW(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'AttSW', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'AttChl',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,AttChl(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'AttChl', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'PARfrac',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,PARfrac(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'PARfrac', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Vp0',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Vp0(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Vp0', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'I_thNH4',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,I_thNH4(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'I_thNH4', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'D_p5NH4',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,D_p5NH4(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'D_p5NH4', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'NitriR',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,NitriR(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'NitriR', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'K_NO3',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,K_NO3(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'K_NO3', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'K_NH4',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,K_NH4(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'K_NH4', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'K_Phy',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,K_Phy(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'K_Phy', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Chl2C_m',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Chl2C_m(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Chl2C_m', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'ChlMin',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,ChlMin(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ChlMin', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'PhyCN',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,PhyCN(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'PhyCN', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'PhyIP',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,PhyIP(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'PhyIP', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'PhyIS',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,PhyIS(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'PhyIS', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'PhyMin',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,PhyMin(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'PhyMin', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'PhyMR',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,PhyMR(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'PhyMR', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

#  ifdef CARBON
        wrt_info=nf_inq_varid(ncid,'ZooAE_C',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,ZooAE_C(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ZooAE_C', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif

        wrt_info=nf_inq_varid(ncid,'ZooAE_N',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,ZooAE_N(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ZooAE_N', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'ZooBM',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,ZooBM(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ZooBM', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'ZooCN',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,ZooCN(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ZooCN', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'ZooER',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,ZooER(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ZooER', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'ZooGGE_C',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,ZooGGE_C(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ZooGGE_C', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'ZooGR',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,ZooGR(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ZooGR', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'ZooMin',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,ZooMin(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ZooMin', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'ZooMR',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,ZooMR(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ZooMR', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'LDeRR',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,LDeRR(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'LDeRR', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'CoagR',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,CoagR(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'CoagR', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'SDeAR',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,SDeAR(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'SDeAR', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'SDeBR',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,SDeBR(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'SDeBR', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'SDeRR',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,SDeRR(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'SDeRR', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'wPhy',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,wPhy(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'wPhy', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

#  ifdef CARBON
        wrt_info=nf_inq_varid(ncid,'wLDeC',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,wLDeC(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'wLDeC', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif

        wrt_info=nf_inq_varid(ncid,'wLDeN',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,wLDeN(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'wLDeN', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

#  ifdef CARBON
        wrt_info=nf_inq_varid(ncid,'wSDeC',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,wSDeC(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'wSDeC', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#  endif

        wrt_info=nf_inq_varid(ncid,'wSDeN',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,wSDeN(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'wSDeN', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

# endif
#endif
!
!-----------------------------------------------------------------------
!  Write out grid variables.
!-----------------------------------------------------------------------
!
!  Grid type switch.
!
        wrt_info=nf_inq_varid(ncid,'spherical',varid)
        write(char1,'(l1)') spherical
        wrt_info=nf_put_var1_text(ncid,varid,1,char1)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'spherical', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
!  Domain Length.
!
        wrt_info=nf_inq_varid(ncid,'xl',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,xl)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'xl', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'el',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,el)
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'el', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#ifdef SOLVE3D
!
!  S-coordinate parameters.
!
        wrt_info=nf_inq_varid(ncid,'theta_s',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,theta_s(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'theta_s', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'theta_b',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,theta_b(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'theta_b', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Tcline',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,Tcline(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Tcline', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'hc',varid)
        wrt_info=nf_put_var1_TYPE(ncid,varid,1,hc(ng))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'hc', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
!  S-coordinate non-dimensional independent variables.
!
        wrt_info=nf_inq_varid(ncid,'sc_r',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,N(ng),                   &
     &                            SCALARS(ng)%sc_r(1))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'sc_r', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'sc_w',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,N(ng)+1,                 &
     &                            SCALARS(ng)%sc_w(0))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'sc_w', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
!  S-coordinate non-dimensional stretching curves.
!
        wrt_info=nf_inq_varid(ncid,'Cs_r',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,N(ng),                   &
     &                            SCALARS(ng)%Cs_r(1))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Cs_r', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF

        wrt_info=nf_inq_varid(ncid,'Cs_w',varid)
        wrt_info=nf_put_vara_TYPE(ncid,varid,1,N(ng)+1,                 &
     &                            SCALARS(ng)%Cs_w(0))
        IF (wrt_info.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Cs_w', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
#endif
!
!  User generic parameters.
!
        IF (Nuser.gt.0) THEN
          wrt_info=nf_inq_varid(ncid,'user',varid)
          wrt_info=nf_put_vara_TYPE(ncid,varid,1,Nuser,user)
          IF (wrt_info.ne.nf_noerr) THEN
            WRITE (stdout,10) 'user', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF
#ifdef STATIONS
!
!  Stations positions.
!
        IF (ncid.eq.ncstaid(ng)) THEN
          wrt_info=nf_inq_varid(ncid,'Ipos',varid)
          wrt_info=nf_put_vara_TYPE(ncid,varid,1,Nstation(ng),          &
     &                              SCALARS(ng)%SposX(1))
          IF (wrt_info.ne.nf_noerr) THEN
            WRITE (stdout,10) 'Ipos', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF

          wrt_info=nf_inq_varid(ncid,'Jpos',varid)
          wrt_info=nf_put_vara_TYPE(ncid,varid,1,Nstation(ng),          &
     &                              SCALARS(ng)%SposY(1))
          IF (wrt_info.ne.nf_noerr) THEN
            WRITE (stdout,10) 'Jpos', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF
#endif
      END IF
#ifndef NO_WRITE_GRID
!
!  Bathymetry.
!
      IF (ncid.ne.ncfltid(ng)) THEN
        scale=1.0_r8
        wrt_info=nf_inq_varid(ncid,'h',varid)
        IF (ncid.eq.ncstaid(ng)) THEN
# ifdef STATIONS
          CALL extract2d (ng, Cgrid, varid, r2dvar,                     &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    scale, GRID(ng)%h,                            &
     &                    Nstation(ng), SCALARS(ng)%SposX,              &
     &                    SCALARS(ng)%SposY, wrk)
          IF (wrtThread) THEN
            wrt_info=nf_put_vara_TYPE(ncid,varid,1,Nstation(ng),wrk)
          END IF
# endif
        ELSE
          wrt_info=nf_fwrite(ng, ncid, varid, 0, r2dvar,                &
     &                       LBi, UBi, LBj, UBj, 1, 1, scale,           &
# ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
# endif
     &                       GRID(ng) % h(LBi,LBj))
        END IF
        IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
          WRITE (stdout,10) 'h', TRIM(ncname)
          exit_flag=3
          RETURN
        END IF
!
!  Coriolis parameter.
!
        IF (ncid.ne.ncstaid(ng)) THEN
          scale=1.0_r8
          wrt_info=nf_inq_varid(ncid,'f',varid)
          wrt_info=nf_fwrite(ng, ncid, varid, 0, r2dvar,                &
     &                       LBi, UBi, LBj, UBj, 1, 1, scale,           &
# ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
# endif
     &                       GRID(ng) % f(LBi,LBj))
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'f', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF
!
!  Curvilinear transformation metrics.
!
        IF (ncid.ne.ncstaid(ng)) THEN
          scale=1.0_r8
          wrt_info=nf_inq_varid(ncid,'pm',varid)
          wrt_info=nf_fwrite(ng, ncid, varid, 0, r2dvar,                &
     &                       LBi, UBi, LBj, UBj, 1, 1, scale,           &
# ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
# endif
     &                       GRID(ng) % pm(LBi,LBj))
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'pm', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF

          wrt_info=nf_inq_varid(ncid,'pn',varid)
          wrt_info=nf_fwrite(ng, ncid, varid, 0, r2dvar,                &
     &                       LBi, UBi, LBj, UBj, 1, 1, scale,           &
# ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
# endif
     &                       GRID(ng) % pn(LBi,LBj))
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'pn', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF
!
!  Longitude and latitude of RHO-points, if appropriate.
!
        IF (spherical) THEN
          scale=1.0_r8
          wrt_info=nf_inq_varid(ncid,'lon_rho',varid)
          IF (ncid.eq.ncstaid(ng)) THEN
# ifdef STATIONS
            CALL extract2d (ng, Cgrid, varid, r2dvar,                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%lonr,                       &
     &                      Nstation(ng), SCALARS(ng)%SposX,            &
     &                      SCALARS(ng)%SposY, wrk)
            IF (wrtThread) THEN
              wrt_info=nf_put_vara_TYPE(ncid,varid,1,Nstation(ng),wrk)
            END IF
# endif
          ELSE
            wrt_info=nf_fwrite(ng, ncid, varid, 0, r2dvar,              &
     &                         LBi, UBi, LBj, UBj, 1, 1, scale,         &
# ifdef MASKING
     &                         GRID(ng) % rmask(LBi,LBj),               &
# endif
     &                         GRID(ng) % lonr(LBi,LBj))
          END IF
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'lon_rho', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF

          wrt_info=nf_inq_varid(ncid,'lat_rho',varid)
          IF (ncid.eq.ncstaid(ng)) THEN
# ifdef STATIONS
            CALL extract2d (ng, Cgrid, varid, r2dvar,                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%latr,                       &
     &                      Nstation(ng), SCALARS(ng)%SposX,            &
     &                      SCALARS(ng)%SposY, wrk)
            IF (wrtThread) THEN
              wrt_info=nf_put_vara_TYPE(ncid,varid,1,Nstation(ng),wrk)
            END IF 
# endif
          ELSE
            wrt_info=nf_fwrite(ng, ncid, varid, 0, r2dvar,              &
     &                         LBi, UBi, LBj, UBj, 1, 1, scale,         &
# ifdef MASKING
     &                         GRID(ng) % rmask(LBi,LBj),               &
# endif
     &                         GRID(ng) % latr(LBi,LBj))
          END IF
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'lat_rho', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        ELSE
          scale=1.0_r8
          wrt_info=nf_inq_varid(ncid,'x_rho',varid)
          IF (ncid.eq.ncstaid(ng)) THEN
# ifdef STATIONS
            CALL extract2d (ng, Cgrid, varid, r2dvar,                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%xr,                         &
     &                      Nstation(ng), SCALARS(ng)%SposX,            &
     &                      SCALARS(ng)%SposY, wrk)
            IF (wrtThread) THEN
              wrt_info=nf_put_vara_TYPE(ncid,varid,1,Nstation(ng),wrk)
            END IF
# endif
          ELSE
            wrt_info=nf_fwrite(ng, ncid, varid, 0, r2dvar,              &
     &                         LBi, UBi, LBj, UBj, 1, 1, scale,         &
# ifdef MASKING
     &                         GRID(ng) % rmask(LBi,LBj),               &
# endif
     &                         GRID(ng) % xr(LBi,LBj))
          END IF
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'x_rho', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF

          wrt_info=nf_inq_varid(ncid,'y_rho',varid)
          IF (ncid.eq.ncstaid(ng)) THEN
# ifdef STATIONS
            CALL extract2d (ng, Cgrid, varid, r2dvar,                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%yr,                         &
     &                      Nstation(ng), SCALARS(ng)%SposX,            &
     &                      SCALARS(ng)%SposY, wrk)
            IF (wrtThread) THEN
              wrt_info=nf_put_vara_TYPE(ncid,varid,1,Nstation(ng),wrk)
            END IF
# endif
          ELSE
            wrt_info=nf_fwrite(ng, ncid, varid, 0, r2dvar,              &
     &                         LBi, UBi, LBj, UBj, 1, 1, scale,         &
# ifdef MASKING
     &                         GRID(ng) % rmask(LBi,LBj),               &
# endif
     &                         GRID(ng) % yr(LBi,LBj))
          END IF
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'y_rho', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF
# ifdef CURVGRID
!
!  Angle between XI-axis and EAST at RHO-points, if appropriate.
!
        IF (spherical) THEN
          scale=1.0_r8
          wrt_info=nf_inq_varid(ncid,'angle',varid)
          IF (ncid.eq.ncstaid(ng)) THEN
# ifdef STATIONS
            CALL extract2d (ng, Cgrid, varid, r2dvar,                   &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      scale, GRID(ng)%angler,                     &
     &                      Nstation(ng), SCALARS(ng)%SposX,            &
     &                      SCALARS(ng)%SposY, wrk)
            IF (wrtThread) THEN
              wrt_info=nf_put_vara_TYPE(ncid,varid,1,Nstation(ng),wrk)
            END IF
# endif
          ELSE
            wrt_info=nf_fwrite(ng, ncid, varid, 0, r2dvar,              &
     &                         LBi, UBi, LBj, UBj, 1, 1, scale,         &
# ifdef MASKING
     &                         GRID(ng) % rmask(LBi,LBj),               &
# endif
     &                         GRID(ng) % angler(LBi,LBj))
          END IF
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'angle', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF
# endif
# ifdef MASKING
!
!  Masking fields at RHO-, U-, V-points, and PSI-points.
!
        IF (ncid.ne.ncstaid(ng)) THEN
          scale=1.0_r8
          wrt_info=nf_inq_varid(ncid,'mask_rho',varid)
          wrt_info=nf_fwrite(ng, ncid, varid, 0, r2dvar,                &
     &                       LBi, UBi, LBj, UBj, 1, 1, scale,           &
# ifdef MASKING
     &                       GRID(ng) % rmask(LBi,LBj),                 &
# endif
     &                       GRID(ng) % rmask(LBi,LBj))
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'mask_rho', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF

          wrt_info=nf_inq_varid(ncid,'mask_u',varid)
          wrt_info=nf_fwrite(ng, ncid, varid, 0, u2dvar,                &
     &                       LBi, UBi, LBj, UBj, 1, 1, scale,           &
# ifdef MASKING
     &                       GRID(ng) % umask(LBi,LBj),                 &
# endif
     &                       GRID(ng) % umask(LBi,LBj))
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'mask_u', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF

          wrt_info=nf_inq_varid(ncid,'mask_v',varid)
          wrt_info=nf_fwrite(ng, ncid, varid, 0, v2dvar,                &
     &                       LBi, UBi, LBj, UBj, 1, 1, scale,           &
# ifdef MASKING
     &                       GRID(ng) % vmask(LBi,LBj),                 &
# endif
     &                       GRID(ng) % vmask(LBi,LBj))
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'mask_v', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF

          wrt_info=nf_inq_varid(ncid,'mask_psi',varid)
          wrt_info=nf_fwrite(ng, ncid, varid, 0, p2dvar,                &
     &                       LBi, UBi, LBj, UBj, 1, 1, scale,           &
# ifdef MASKING
     &                       GRID(ng) % pmask(LBi,LBj),                 &
# endif
     &                       GRID(ng) % pmask(LBi,LBj))
          IF (wrtThread.and.(wrt_info.ne.nf_noerr)) THEN
            WRITE (stdout,10) 'mask_psi', TRIM(ncname)
            exit_flag=3
            RETURN
          END IF
        END IF
# endif
      END IF
#endif
!
  10  format(/,' WRT_INFO - error while writing variable: ',a,/,        &
     &       12x,'into NetCDF file: ',a)
      RETURN
      END FUNCTION wrt_info
