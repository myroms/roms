#include "cppdefs.h"
      FUNCTION def_var (ncid,Vid,Vtype,nVdim,Vdim,Aval,Vinfo,ncname)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine defines the requested NetCDF variable.                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ncid     NetCDF file ID (integer).                               !
!     Vtype    NetCDF variable type (integer).                         !
!     mVdim    Number of variable dimensions (integer; 0=scalar).      !
!     Vdim     Dimensions IDs for this variable (integer vector).      !
!     Aval     Attribute values (real vector):                         !
!                Aval(1)   =>  Add offset value.                       !
!                Aval(2)   =>  Valid minimum value.                    !
!                Aval(3)   =>  Valid maximum value.                    !
!                Aval(4)   =>  Missing value.                          !
!                Aval(5)   =>  C-grid variable type.                   !
!                Aval(6)   =>  Fill value.                             !
!     Vinfo    Variable information (character array):                 !
!                Vinfo( 1) =>  Variable name.                          !
!                Vinfo( 2) =>  Variable "longname" attribute.          !
!                Vinfo( 3) =>  Variable "units" attribute.             !
!                Vinfo( 4) =>  Variable "calendar" attribute.          !
!                Vinfo( 5) =>  Variable "valid_min" attribute.         !
!                Vinfo( 6) =>  Variable "valid_max" attribute.         !
!                Vinfo( 7) =>  Variable "option_T" attribute.          !
!                Vinfo( 8) =>  Variable "option_F" attribute.          !
!                Vinfo( 9) =>  Variable "option_0" attribute.          !
!                Vinfo(10) =>  Variable "option_1" attribute.          !
!                Vinfo(11) =>  Variable "negative_value" attribute.    !
!                Vinfo(12) =>  Variable "positive_value" attribute.    !
!                Vinfo(13) =>  Variable "cycle" attribute.             !
!                Vinfo(14) =>  Variable "field" attribute.             !
!                Vinfo(15) =>  Variable "positions" attribute.         !
!                Vinfo(16) =>  Variable "time" attribute.              !
!                Vinfo(17) =>  Variable "missing_value" attribute.     !
!                Vinfo(18) =>  Variable "add_offset" attribute.        !
!                Vinfo(19) =>  Variable "size_class" attribute.        !
!                Vinfo(20) =>  Variable "water_points" attribute.      !
!                Vinfo(21) =>  Variable "standard_name" attribute.     !
!                Vinfo(22) =>  Variable "coordinates" attribute.       !
!                Vinfo(23) =>  Variable "formula_terms" attribute.     !
!                Vinfo(24) =>  Variable "_FillValue" attribute.        !
!     ncname   NetCDF file name.                                       !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     def_var  Error flag (integer).                                   !
!     Vid      NetCDF variable ID (integer).                           !
!                                                                      !
!  Notice that arrays "Aval" and "Vinfo" is destroyed on output to     !
!  facilitate the definition of the next variable.                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent (in) :: ncid, Vtype, nVdim
      integer, dimension(*), intent(in) :: Vdim
      integer, intent (out) :: Vid

      real(r8), dimension(*), intent(inout) :: Aval

      character (len=*), intent(in) :: ncname
      character (len=*), intent(inout) :: Vinfo(24)
!
!  Local variable declarations.
!
      integer :: latt

      integer :: def_var

      character (len=40) text
!
!-----------------------------------------------------------------------
!  Define requested variable and its attributes.
!-----------------------------------------------------------------------
!
!  Define variable.
!
      IF (LEN_TRIM(Vinfo(1)).gt.0) THEN
        def_var=nf_def_var(ncid,TRIM(Vinfo(1)),Vtype,nVdim,Vdim,Vid)
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,10) TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
      END IF
!
!  Define "longname" attribute.
!
      latt=LEN_TRIM(Vinfo(2))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'long_name',                   &
     &                          latt,TRIM(Vinfo(2)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'long_name', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(2)=' '
      END IF
!
!  Define "size_class" attribute.
!
      latt=LEN_TRIM(Vinfo(19))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'size_class',                  &
     &                          latt,TRIM(Vinfo(19)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'size_class', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(19)=' '
      END IF
!
!  If applicable, define "units" attribute.
!
      latt=LEN_TRIM(Vinfo(3))
      IF (latt.gt.0) THEN
        IF (TRIM(Vinfo(3)).ne.'nondimensional') THEN
          def_var=nf_put_att_text(ncid,Vid,'units',                     &
     &                            latt,TRIM(Vinfo(3)))
          IF (def_var.ne.nf_noerr) THEN
            WRITE (stdout,20) 'units', TRIM(Vinfo(1)), TRIM(ncname)
            exit_flag=3
            ioerror=def_var
            RETURN
          END IF
        END IF
        Vinfo(3)=' '
      END IF
!
!  If applicable, define "calendar" attribute.
!
      latt=LEN_TRIM(Vinfo(4))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'calendar',                    &
     &                          latt,TRIM(Vinfo(4)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'calendar', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(4)=' '
      END IF
!
!  If applicable, define "valid_min" attribute.
!
      latt=LEN_TRIM(Vinfo(5))
      IF (latt.gt.0) THEN
        IF (Vtype.eq.nf_int) THEN
          def_var=nf_put_att_int(ncid,Vid,TRIM(Vinfo(5)),               &
     &                           Vtype,1,INT(Aval(2)))
        ELSE
          def_var=nf_put_att_TYPE(ncid,Vid,TRIM(Vinfo(5)),              &
     &                            Vtype,1,Aval(2))
        END IF
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) TRIM(Vinfo(5)), TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Aval(2)=0.0_r8
        Vinfo(5)=' '
      END IF
!
!  If applicable, define "valid_max" attribute.
!
      latt=LEN_TRIM(Vinfo(6))
      IF (latt.gt.0) THEN
        IF (Vtype.eq.nf_int) THEN
          def_var=nf_put_att_int(ncid,Vid,TRIM(Vinfo(6)),               &
     &                           Vtype,1,INT(Aval(3)))
        ELSE
          def_var=nf_put_att_TYPE(ncid,Vid,TRIM(Vinfo(6)),              &
     &                            Vtype,1,Aval(3))
        END IF
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) TRIM(Vinfo(6)), TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Aval(3)=0.0_r8
        Vinfo(6)=' '
      END IF
!
!  If applicable, define "option_T" attribute.
!
      latt=LEN_TRIM(Vinfo(7))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'option_T',                    &
     &                          latt,TRIM(Vinfo(7)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'option_T', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(7)=' '
      END IF
!
!  If applicable, define "option_F" attribute.
!
      latt=LEN_TRIM(Vinfo(8))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'option_F',                    &
     &                          latt,TRIM(Vinfo(8)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'option_F', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(8)=' '
      END IF
!
!  If applicable, define "option_0" attribute.
!
      latt=LEN_TRIM(Vinfo(9))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'option_0',                    &
     &                          latt,TRIM(Vinfo(9)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'option_0', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(9)=' '
      END IF
!
!  If applicable, define "option_1" attribute.
!
      latt=LEN_TRIM(Vinfo(10))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'option_1',                    &
     &                          latt,TRIM(Vinfo(10)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'option_1', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(10)=' '
      END IF
!
!  If applicable, define "negative" attribute.
!
      latt=LEN_TRIM(Vinfo(11))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'negative_value',              &
     &                          latt,TRIM(Vinfo(11)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'negative_value', TRIM(Vinfo(1)),           &
     &                      TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(11)=' '
      END IF
!
!  If applicable, define "option_F" attribute.
!
      latt=LEN_TRIM(Vinfo(12))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'positive_value',              &
     &                          latt,TRIM(Vinfo(12)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'positive_value', TRIM(Vinfo(1)),           &
     &                      TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(12)=' '
      END IF
!
!  If applicable, define "cycle" attribute.
!
      latt=LEN_TRIM(Vinfo(13))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'cycle',                       &
     &                          latt,TRIM(Vinfo(13)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'cycle', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(13)=' '
      END IF
!
!  If applicable, define "positions" attribute.
!
      latt=LEN_TRIM(Vinfo(15))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'positions',                   &
     &                          latt,TRIM(Vinfo(15)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'positions', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(15)=' '
      END IF
!
!  If applicable, define "time" attribute.
!
      latt=LEN_TRIM(Vinfo(16))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'time',                        &
     &                          latt,TRIM(Vinfo(16)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'time', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(16)=' '
      END IF
!
!  If applicable, define "missing_value" attribute.
!
      latt=LEN_TRIM(Vinfo(17))
      IF (latt.gt.0) THEN
        IF (Vtype.eq.nf_int) THEN
          def_var=nf_put_att_int(ncid,Vid,TRIM(Vinfo(17)),              &
     &                           Vtype,1,INT(Aval(4)))
        ELSE
          def_var=nf_put_att_TYPE(ncid,Vid,TRIM(Vinfo(17)),             &
     &                            Vtype,1,Aval(4))
        END IF
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) TRIM(Vinfo(17)), TRIM(Vinfo(1)),            &
     &                      TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Aval(4)=0.0_r8
        Vinfo(17)=' '
      END IF
!
!  If applicable, define "add_offset" attribute.
!
      latt=LEN_TRIM(Vinfo(18))
      IF (latt.gt.0) THEN
        IF (Vtype.eq.nf_int) THEN
          def_var=nf_put_att_int(ncid,Vid,TRIM(Vinfo(18)),              &
     &                           Vtype,1,INT(Aval(1)))
        ELSE
          def_var=nf_put_att_TYPE(ncid,Vid,TRIM(Vinfo(18)),             &
     &                            Vtype,1,Aval(1))
        END IF
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) TRIM(Vinfo(18)), TRIM(Vinfo(1)),            &
     &                      TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Aval(1)=0.0_r8
        Vinfo(18)=' '
      END IF
!
!  If applicable, define "_FillValue" attribute.
!
      latt=LEN_TRIM(Vinfo(24))
      IF (latt.gt.0) THEN
        IF (Vtype.eq.nf_int) THEN
          def_var=nf_put_att_int(ncid,Vid,TRIM(Vinfo(24)),              &
     &                           Vtype,1,INT(Aval(6)))
        ELSE
          def_var=nf_put_att_TYPE(ncid,Vid,TRIM(Vinfo(24)),             &
     &                            Vtype,1,Aval(6))
        END IF
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) TRIM(Vinfo(24)), TRIM(Vinfo(1)),            &
     &                      TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Aval(6)=0.0_r8
        Vinfo(24)=' '
      END IF
!
!  If applicable, define "water_points" attribute.
!
      latt=LEN_TRIM(Vinfo(20))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'water_points',                &
     &                          latt,TRIM(Vinfo(20)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'water_points', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(20)=' '
      END IF
!
!  If applicable, define "standard_name" attribute.
!
      latt=LEN_TRIM(Vinfo(21))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'standard_name',               &
     &                          latt,TRIM(Vinfo(21)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'standard_name', TRIM(Vinfo(1)),            &
     &                      TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(21)=' '
      END IF
!
!  If applicable, define "coordinates" attribute.
!
      latt=LEN_TRIM(Vinfo(22))
      IF (latt.gt.0) THEN
        IF (spherical) THEN
          IF ((INT(Aval(5)).eq.r2dvar).or.                              &
     &        (INT(Aval(5)).eq.b3dvar).or.                              &
     &        (INT(Aval(5)).eq.r3dvar).or.                              &
     &        (INT(Aval(5)).eq.w3dvar)) THEN
            text='lat_rho lon_rho'
          ELSE IF ((INT(Aval(5)).eq.u2dvar).or.                         &
     &             (INT(Aval(5)).eq.u3dvar)) THEN
            text='lat_u lon_u'
          ELSE IF ((INT(Aval(5)).eq.v2dvar).or.                         &
     &             (INT(Aval(5)).eq.v3dvar)) THEN
            text='lat_v lon_v'
          ELSE IF ((INT(Aval(5)).eq.p2dvar).or.                         &
     &             (INT(Aval(5)).eq.p3dvar)) THEN
            text='lat_psi lon_psi'
          END IF
        ELSE
          IF ((INT(Aval(5)).eq.r2dvar).or.                              &
     &        (INT(Aval(5)).eq.b3dvar).or.                              &
     &        (INT(Aval(5)).eq.r3dvar).or.                              &
     &        (INT(Aval(5)).eq.w3dvar)) THEN
            text='x_rho y_rho'
          ELSE IF ((INT(Aval(5)).eq.u2dvar).or.                         &
     &             (INT(Aval(5)).eq.u3dvar)) THEN
            text='x_u y_u'
          ELSE IF ((INT(Aval(5)).eq.v2dvar).or.                         &
     &             (INT(Aval(5)).eq.v3dvar)) THEN
            text='x_v y_v'
          ELSE IF ((INT(Aval(5)).eq.p2dvar).or.                         &
     &             (INT(Aval(5)).eq.p3dvar)) THEN
            text='x_psi y_psi'
          END IF
        END IF
        latt=LEN_TRIM(text)
        def_var=nf_put_att_text(ncid,Vid,TRIM(Vinfo(22)),               &
     &                          latt,TRIM(text))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'coordinates', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(22)=' '
        Aval(5)=0.0_r8
      END IF
!
!  If applicable, define "formula_terms" attribute.
!
      latt=LEN_TRIM(Vinfo(23))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'formula_terms',               &
     &                          latt,TRIM(Vinfo(23)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'formula_terms', TRIM(Vinfo(1)),            &
     &                      TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(23)=' '
      END IF
!
!  If applicable, define "field" attribute (always last).
!
      latt=LEN_TRIM(Vinfo(14))
      IF (latt.gt.0) THEN
        def_var=nf_put_att_text(ncid,Vid,'field',                       &
     &                          latt,TRIM(Vinfo(14)))
        IF (def_var.ne.nf_noerr) THEN
          WRITE (stdout,20) 'field', TRIM(Vinfo(1)), TRIM(ncname)
          exit_flag=3
          ioerror=def_var
          RETURN
        END IF
        Vinfo(14)=' '
      END IF
!
! Clean variable name.
!
      Vinfo(1)=' '
!
 10   FORMAT (/,' DEF_VAR - unable to define variable: ',a,/,           &
     &        11x,'in NetCDF file: ',a)
 20   FORMAT (/,'DEF_VAR - error while defining attribute: ',a,         &
     &        ' for variable: ',a,/,11x,'in NetCDF file: ',a)
      RETURN
      END FUNCTION def_var
