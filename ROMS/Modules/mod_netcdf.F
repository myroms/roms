#include "cppdefs.h"
      MODULE mod_netcdf
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2009 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This MODULE contains all NetCDF variables definitions. It also      !
!  contains several variables and routines to facilitate  generic      !
!  manipulations of NetCDF data.  Usually, several NetCDF library      !
!  calls are required to inquire and read a dimension or variable.     !
!  These routines provide a single interface for such operations.      !
!                                                                      !
!=======================================================================
!
      USE mod_kinds
      USE netcdf
!
      implicit none
!
      INTERFACE netcdf_get_fvar
        MODULE PROCEDURE netcdf_get_fvar_0d
        MODULE PROCEDURE netcdf_get_fvar_1d
        MODULE PROCEDURE netcdf_get_fvar_2d
        MODULE PROCEDURE netcdf_get_fvar_3d
        MODULE PROCEDURE netcdf_get_fvar_4d
      END INTERFACE netcdf_get_fvar

      INTERFACE netcdf_get_ivar
        MODULE PROCEDURE netcdf_get_ivar_0d
        MODULE PROCEDURE netcdf_get_ivar_1d
        MODULE PROCEDURE netcdf_get_ivar_2d
      END INTERFACE netcdf_get_ivar

      INTERFACE netcdf_get_svar
        MODULE PROCEDURE netcdf_get_svar_0d
        MODULE PROCEDURE netcdf_get_svar_1d
      END INTERFACE netcdf_get_svar

      INTERFACE netcdf_put_fvar
        MODULE PROCEDURE netcdf_put_fvar_0d
        MODULE PROCEDURE netcdf_put_fvar_1d
        MODULE PROCEDURE netcdf_put_fvar_2d
        MODULE PROCEDURE netcdf_put_fvar_3d
        MODULE PROCEDURE netcdf_put_fvar_4d
      END INTERFACE netcdf_put_fvar

      INTERFACE netcdf_put_ivar
        MODULE PROCEDURE netcdf_put_ivar_0d
        MODULE PROCEDURE netcdf_put_ivar_1d
        MODULE PROCEDURE netcdf_put_ivar_2d
      END INTERFACE netcdf_put_ivar

      INTERFACE netcdf_put_svar
        MODULE PROCEDURE netcdf_put_svar_0d
        MODULE PROCEDURE netcdf_put_svar_1d
      END INTERFACE netcdf_put_svar
!
!  Local dimension parameters.
!
      integer, parameter :: Mdims = 50  ! maximun number of dimensions
      integer, parameter :: Mvars = 500 ! maximun number of variables
      integer, parameter :: NvarD = 5   ! number of variable dimensions
      integer, parameter :: NvarA = 10  ! number of variable attributes
!
!  Generic information about current NetCDF for all dimensions and
!  all variables.

      integer :: n_dim                  ! number of dimensions
      integer :: n_var                  ! number of variables 
      integer :: n_gatt                 ! number of global attributes
      integer :: rec_id                 ! unlimited dimension ID
      integer :: rec_size               ! unlimited dimension value
      integer :: dim_id(Mdims)          ! dimensions ID
      integer :: dim_size(Mdims)        ! dimensions value
      integer :: var_id(Mvars)          ! variables ID
      integer :: var_natt(Mvars)        ! variables number of attributes
      integer :: var_flag(Mvars)        ! Variables water points flag
      integer :: var_type(Mvars)        ! variables external data type
      integer :: var_ndim(Mvars)        ! variables number of dimensions
      integer :: var_dim(NvarD,Mvars)   ! variables dimensions ID
!
      character (len=40) :: dim_name(Mdims)      ! dimensions name
      character (len=40) :: var_name(Mvars)      ! variables name
!
!  Generic information about requested current variable.
!
      integer :: n_vdim                 ! number of variable dimensions
      integer :: n_vatt                 ! number of variable attributes
      integer :: var_kind               ! external data type
      integer :: var_Dids(NvarD)        ! dimensions ID
      integer :: var_Dsize(NvarD)       ! dimensions values
      integer :: var_Aint(NvarA)        ! attribute integer values
      real(r8) :: var_Afloat(NvarA)     ! attribute float values
!
      character (len=40) :: var_Aname(NvarA)     ! Attribute names
      character (len=40) :: var_Dname(NvarD)     ! dimension names
      character (len=80) :: var_Achar(NvarA)     ! Attribute char values
!
!  External data representation for floating-point variables.
!
#ifdef OUT_DOUBLE
      integer, parameter :: NF_FOUT = nf90_double
#else
      integer, parameter :: NF_FOUT = nf90_real
#endif
#ifdef RST_SINGLE
      integer, parameter :: NF_FRST = nf90_real
#else
      integer, parameter :: NF_FRST = nf90_double
#endif
#ifdef DOUBLE_PRECISION
      integer, parameter :: NF_TYPE = nf90_double
#else
      integer, parameter :: NF_TYPE = nf90_real
#endif
!
!  Netcdf file creation mode flag.
!
#ifdef NETCDF4
      integer :: CMODE = nf90_netcdf4      ! NetCDF-4/HDF5 format file
#else
      integer :: CMODE = nf90_clobber      ! NetCDF classic format file
#endif
!
!  Set NetCDF-4/HDF5 deflate (file compression) parameters.
!
      integer :: shuffle = 1
      integer :: deflate = 1
      integer :: deflate_level = 1

      CONTAINS

      SUBROUTINE netcdf_get_dim (ng, model, ncname, DimName, DimSize,   &
     &                           DimID)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  All the dimension information is stored in the module variables     !
!  declared above.  In addition, if a particular dimension name is     !
!  provided, this routine returns the requested information in the     !
!  optional arguments.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     DimName      Requested dimension name (string, OPTIONAL)         !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     DimSize      Size of requested dimension (integer, OPTIONAL)     !
!     DimID        ID od requested dimension (integer, OPTIONAL)       !
!                                                                      !
!  Other information stored in this module:                            !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     rec_size     Size of unlimited dimension                         !
!     dim_name     Dimensions name (1:n_dim)                           !
!     dim_id       Dimensions ID (1:n_dim)                             !
!     dim_size     Dimensions value (1:n_dim)                          !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in), optional :: DimName

      integer, intent(out), optional :: DimSize
      integer, intent(out), optional :: DimID
!
!  Local variable declarations.
!
      logical :: foundit

      integer :: ncid, i, j, status
      integer :: myID, myValue

#ifdef DISTRIBUTE
      integer, dimension(5) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      n_dim=0
      n_var=0
      n_gatt=0
      rec_id=-1
      rec_size=0
      dim_id=0
      dim_size=0
      DO i=1,Mdims
        DO j=1,LEN(dim_name(1))
          dim_name(i)(j:j)=' '
        END DO
      END DO
!
!  Open NetCDF for reading.
!
      IF (InpThread) THEN
        status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
!
!  Inquire NetCDF file.
!
        IF (status.eq.nf90_noerr) THEN
          status=nf90_inquire(ncid, n_dim, n_var, n_gatt, rec_id)
          IF ((status.eq.nf90_noerr).and.(n_dim.le.Mdims)) THEN   
#ifdef DISTRIBUTE
            ibuffer(1)=n_dim
            ibuffer(2)=n_var
            ibuffer(3)=n_gatt
            ibuffer(4)=rec_id
#endif
!
!  Inquire about dimensions: names, ID, and size.
!
            rec_size=-1
            DO i=1,n_dim
              dim_id(i)=i
              status=nf90_inquire_dimension(ncid, dim_id(i),            &
     &                                      dim_name(i), dim_size(i))
              IF (status.ne.nf90_noerr) THEN
                WRITE (stdout,10) dim_id(i), TRIM(ncname),              &
     &                            TRIM(SourceFile)
                exit_flag=2
                ioerror=status
                EXIT
              END IF
              IF (dim_id(i).eq.rec_id) THEN
                rec_size=dim_size(i)
#ifdef DISTRIBUTE
                ibuffer(5)=rec_size
#endif
              END IF
            END DO
!
!  Close input NetCDF file.
!
            status=nf90_close(ncid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(ncname), TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE      
            IF (n_dim.gt.Mdims) THEN
              WRITE (stdout,30) ' Mdims = ', Mdims, n_dim
              exit_flag=2
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          END IF
        ELSE
          WRITE (stdout,50) TRIM(ncname), TRIM(SourceFile)
          ioerror=status
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast dimension to all processors in the group.
!
      CALL mp_bcasti (ng, model, exit_flag)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, ibuffer)
        n_dim=ibuffer(1)
        n_var=ibuffer(2)
        n_gatt=ibuffer(3)
        rec_id=ibuffer(4)
        rec_size=ibuffer(5)
        CALL mp_bcasti (ng, model, dim_id)
        CALL mp_bcasti (ng, model, dim_size)
        CALL mp_bcasts (ng, model, dim_name)
      END IF
#endif
!
!  Load requested information.
!
      foundit=.FALSE.
      IF (PRESENT(DimName)) THEN
        DO i=1,n_dim
          IF (TRIM(dim_name(i)).eq.TRIM(DimName)) THEN
            foundit=.TRUE.
            myID=dim_id(i)
            myValue=dim_size(i)
          END IF
        END DO
        IF (foundit) THEN
          IF (PRESENT(DimSize)) THEN
            DimSize=myValue
          END IF
          IF (PRESENT(DimID)) THEN
            DimID=myID
          END IF
        ELSE
          WRITE (stdout,60) TRIM(DimName), TRIM(ncname)
          exit_flag=2
          ioerror=status
        END IF
      END IF
!
  10  FORMAT (/,' NETCDF_GET_DIM - error while reading dimension ID:',  &
     &        2x,i3,/,18x,'in input file:',2x,a,/,18x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_DIM - unable to close input NetCDF file:', &
     &        2x,a,/,18x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_DIM - too small dimension parameter,',a,   &
     &        2i5,/,18x,'change file  nc_util.F  and recompile')
  40  FORMAT (/,' NETCDF_GET_DIM - unable to inquire about contents',   &
     &          ' of input NetCDF file:',2x,a,/,18x,'call from:',2x,a)
  50  FORMAT (/,' NETCDF_GET_DIM - unable to open input NetCDF file:',  &
     &        2x,a,/,18x,'call from:',2x,a)
  60  FORMAT (/,' NETCDF_GET_DIM - requested dimension: ',a,/18x,       &
     &        'not found in input file:',2x,a,/,18x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_dim

      SUBROUTINE netcdf_check_dim (ng, model, ncname)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  It checks the file dimensions against model dimension parameters    !
!  for consitency.  All the dimensions information is stored in the    !
!  module variables declared above.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!                                                                      !
!  On output the following information is stored in this module:       !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     rec_size     Size of unlimited dimension                         !
!     dim_name     Dimensions name (1:n_dim)                           !
!     dim_id       Dimensions ID (1:n_dim)                             !
!     dim_size     Dimensions value (1:n_dim)                          !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      logical :: foundit

      integer :: ncid, i, status
      integer :: myID, myValue

#ifdef DISTRIBUTE
      integer, dimension(5) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
      CALL netcdf_get_dim (ng, model, ncname)
      IF (exit_flag.ne.NoError) RETURN
!
!-----------------------------------------------------------------------
!  Check dimensions for consistency.
!-----------------------------------------------------------------------
!
      DO i=1,n_dim
        SELECT CASE (TRIM(ADJUSTL(dim_name(i))))
          CASE ('xi_psi')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_psi) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_psi,        &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_psi')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_psi) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_psi,       &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_rho')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_rho) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_rho,        &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_rho')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_rho) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_rho,       &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_u')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_u) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_u,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_u')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_u) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_u,         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('xi_v')
            IF (dim_size(i).ne.IOBOUNDS(ng)%xi_v) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%xi_v,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('eta_v')
            IF (dim_size(i).ne.IOBOUNDS(ng)%eta_v) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%eta_v,         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
          CASE ('IorJ')
            IF (dim_size(i).ne.IOBOUNDS(ng)%IorJ) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i),                &
     &                                      IOBOUNDS(ng)%IorJ,          &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
#ifdef SOLVE3D
          CASE ('s_rho')
            IF (dim_size(i).ne.N(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i), N(ng),         &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
#endif
#ifdef ADJUST_BOUNDARY
          CASE ('obc_adjust')
            IF (dim_size(i).ne.Nbrec(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i), Nbrec(ng),     &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
#endif
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          CASE ('frc_adjust')
            IF (dim_size(i).ne.Nfrec(ng)) THEN
              IF (Master) WRITE (stdout,10) TRIM(dim_name(i)),          &
     &                                      dim_size(i), Nfrec(ng),     &
     &                                      TRIM(ncname)
              exit_flag=2
              EXIT
            END IF
#endif
        END SELECT
      END DO

 10   FORMAT (/,' NETCDF_CHECK_DIM - inconsistent size of dimension: ', &
     &        a,2x,2i5,/,20x,'in file: ',a)

      RETURN
      END SUBROUTINE netcdf_check_dim

      SUBROUTINE netcdf_inq_var (ng, model, ncname, myVarName,          &
     &                           SearchVar, VarID, nVarDim, nVarAtt)
!
!=======================================================================
!                                                                      !
!  This routine inquires a NetCDF file dimensions names and values.    !
!  All the dimension information is stored in the module variables     !
!  declared above.  In addition, if a particular dimension name is     !
!  provided, this routine returns the requested information in the     !
!  optional arguments.                                                 !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Requested variable name (string, OPTIONAL)          !
!     SearchVar    Switch used when searching a variable over          !
!                    multiple NetCDF files (logical, OPTIONAL)         !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     VarID        ID of requested variable (integer, OPTIONAL)        !
!     nVarDim      Number of variable dimensions (integer, OPTIONAL)   !
!     nVarAtt      Number of variable attributes (integer, OPTIONAL)   !
!                                                                      !
!  Other information stored in this module:                            !
!                                                                      !
!     n_dim        Number of dimensions                                !
!     n_var        Number of variables                                 !
!     n_gatt       Number of global attributes                         !
!     rec_id       Unlimited dimension ID                              !
!     var_name     Variables name (1:n_var)                            !
!     var_id       Variables ID (1:n_var)                              !
!     var_natt     Variables number of attributes (1:n_var)            !
!     var_flag     Variables flag [1=full field, -1=water points only] !
!     var_type     Variables external data type (1:n_var)              !
!     var_ndim     Variables number of dimensions (1:n_var)            !
!     var_dim      Variables dimensions ID (:,1:n_var)                 !
!                                                                      !
!  If the OPTIONAL argument myVarName is provided, the following       !
!  information for requested variable is also stored:                  !
!                                                                      !
!     n_vatt       Number of variable attributes                       !
!     n_vdim       Number of variable dimensions                       !
!     var_kind     Variable external data type                         !
!     var_Aname    Variable attribute names (1:n_vatt)                 !
!     var_Achar    Variable string attribute values (1:n_vatt)         !
!     var_Afloat   Variable float attribute values (1:n_vatt)          !
!     var_Aint     Variable integer attribute values (1:n_vatt)        !
!     var_Dids     Variable dimensions ID (1:n_vdim)                   !
!     var_Dname    Variable dimensions name (1:n_vdim)                 !
!     var_Dsize    Variable dimensions size (1:n_vdim)                 !
!                                                                      !
!     WARNING: This is information is rewritten during each CALL.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_bcastf, mp_bcasts
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in), optional :: myVarName

      logical, intent(out), optional :: SearchVar

      integer, intent(out), optional :: VarID
      integer, intent(out), optional :: nVarDim
      integer, intent(out), optional :: nVarAtt
!
!  Local variable declarations.
!
      logical :: foundit, WriteError

      integer :: ncid, i, j, status
      integer :: my_Alen, my_Atype, my_id

#ifdef DISTRIBUTE
      integer, dimension(4) :: ibuffer
#endif

      character (len=80) :: text
!
!-----------------------------------------------------------------------
!  Inquire about the NetCDF dimensions (names and values).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      n_dim=0
      n_var=0
      n_gatt=0
      rec_id=-1
      var_id=0
      var_natt=0
      var_flag=0
      var_type=0
      var_ndim=0
      var_dim=0
      DO i=1,Mvars
        DO j=1,LEN(var_name(1))
          var_name(i)(j:j)=' '
        END DO
      END DO
!
!  Open NetCDF for reading.
!
      IF (InpThread) THEN
        status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
!
!  Inquire NetCDF file.
!
        IF (status.eq.nf90_noerr) THEN
          status=nf90_inquire(ncid, n_dim, n_var, n_gatt, rec_id)
          IF ((status.eq.nf90_noerr).and.(n_var.le.Mvars)) THEN
#ifdef DISTRIBUTE
            ibuffer(1)=n_dim
            ibuffer(2)=n_var
            ibuffer(3)=n_gatt
            ibuffer(4)=rec_id
#endif
!
!  Inquire about variables: name, ID, dimensions, data type, and number
!  of attributes.
!
            DO i=1,n_var
              var_id(i)=i
              var_flag(i)=1
              status=nf90_inquire_variable(ncid, var_id(i),             &
     &                                     var_name(i), var_type(i),    &
     &                                     var_ndim(i), var_dim(:,i),   &
     &                                     var_natt(i))
              IF (status.eq.nf90_noerr) THEN
                DO j=1,var_natt(i)
                  status=nf90_inq_attname(ncid, var_id(i), j,           &
     &                                    var_Aname(j))
                  IF (status.eq.nf90_noerr) THEN
                    IF (TRIM(var_Aname(j)).eq.'water_points'.and.       &
     &                  (var_ndim(i).gt.0)) THEN
                      var_flag(i)=-1
                    END IF
                  ELSE
                    WRITE (stdout,70) j, TRIM(var_name(i)),             &
     &                                TRIM(ncname), TRIM(SourceFile)
                    exit_flag=2
                    ioerror=status
                    EXIT
                  END IF
                END DO
              ELSE
                WRITE (stdout,10) var_id(i), TRIM(ncname),              &
     &                            TRIM(SourceFile)
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END DO
          ELSE
            IF (n_var.gt.Mvars) THEN
              WRITE (stdout,20) 'Mvars = ', Mvars, n_var
              exit_flag=2
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,30) TRIM(ncname), TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          END IF
        ELSE
          WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
          exit_flag=4
          ioerror=status
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast dimension to all processors in the group.
!
      CALL mp_bcasti (ng, model, exit_flag)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, ibuffer)
        n_dim=ibuffer(1)
        n_var=ibuffer(2)
        n_gatt=ibuffer(3)
        rec_id=ibuffer(4)
        CALL mp_bcasti (ng, model, var_id)
        CALL mp_bcasti (ng, model, var_flag)
        CALL mp_bcasti (ng, model, var_type)
        CALL mp_bcasti (ng, model, var_ndim)
        CALL mp_bcasti (ng, model, var_natt)
        CALL mp_bcasti (ng, model, var_dim)
        CALL mp_bcasts (ng, model, var_name)
      END IF
#endif
!
!  Load requested requested variable information.
!
      IF (exit_flag.eq.NoError) THEN
        foundit=.FALSE.
        IF (PRESENT(myVarName)) THEN
          var_Dids=-1
          var_Dsize=0
          var_Aint=0
          var_Afloat=0.0_r8
          DO i=1,NvarA
            DO j=1,LEN(var_Aname(1))
              var_Aname(i)(j:j)=' '
            END DO
            DO j=1,LEN(var_Achar(1))
              var_Achar(i)(j:j)=' '
            END DO
          END DO
          DO i=1,NvarD
            DO j=1,LEN(var_Dname(1))
              var_Dname(i)(j:j)=' '
            END DO
          END DO
!
          DO i=1,n_var
            IF (TRIM(var_name(i)).eq.TRIM(myVarName)) THEN
              foundit=.TRUE.
              my_id=var_id(i)
              n_vatt=var_natt(i)
              n_vdim=var_ndim(i)
              DO j=1,n_vdim
                var_Dids(j)=var_dim(j,i)
              END DO
              var_kind=var_type(i)
            END IF
          END DO
          IF (foundit) THEN        
            IF (PRESENT(VarID)) THEN
              VarID=my_id
            END IF
            IF (PRESENT(nVarDim)) THEN
              nVarDim=n_vdim
            END IF
            IF (PRESENT(nVarAtt)) THEN
              nVarAtt=n_vatt
            END IF
          END IF
!
!  If founded requested variable, inquire about is dimensions and
!  attributes.
!
          IF (foundit.and.InpThread) THEN
            DO i=1,n_vdim
              status=nf90_inquire_dimension(ncid, var_Dids(i),          &
     &                                      var_Dname(i), var_Dsize(i))
              IF (status.ne.nf90_noerr) THEN
                WRITE (stdout,50) i, TRIM(myVarName), TRIM(ncname),     &
     &                            TRIM(SourceFile)
                exit_flag=2
                ioerror=status
                EXIT
              END IF
            END DO
            IF (status.eq.nf90_noerr) THEN
              DO i=1,n_vatt
                status=nf90_inq_attname(ncid, my_id, i, var_Aname(i))
                IF (status.eq.nf90_noerr) THEN
                  status=nf90_inquire_attribute(ncid, my_id,            &
     &                                          TRIM(var_Aname(i)),     &
     &                                          xtype = my_Atype,       &
     &                                          len = my_Alen)
                  IF (status.eq.nf90_noerr) THEN
                    IF (my_Atype.eq.NF90_INT) THEN
                      status=nf90_get_att(ncid, my_id,                  &
     &                                    TRIM(var_Aname(i)),           &
     &                                    var_Aint(i))
                      IF (status.ne.nf90_noerr) THEN
                        WRITE (stdout,60) 'integer',                    &
     &                                    TRIM(var_Aname(i)),           &
     &                                    TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                    ELSE IF ((my_Atype.eq.NF90_FLOAT ).or.              &
     &                       (my_Atype.eq.NF90_DOUBLE)) THEN
                      status=nf90_get_att(ncid, my_id,                  &
     &                                    TRIM(var_Aname(i)),           &
     &                                    var_Afloat(i))
                      IF (status.ne.nf90_noerr) THEN
                        WRITE (stdout,60) 'float',                      &
     &                                    TRIM(var_Aname(i)),           &
     &                                    TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                    ELSE IF (my_Atype.eq.NF90_CHAR) THEN
                      status=nf90_get_att(ncid, my_id,                  &
     &                                    TRIM(var_Aname(i)),           &
     &                                    text(1:my_Alen))
                      IF (status.ne.nf90_noerr) THEN
                        WRITE (stdout,60) 'string',                     &
     &                                    TRIM(var_Aname(i)),           &
     &                                    TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                      var_Achar(i)=text(1:my_Alen)
                    END IF
                  ELSE
                    WRITE (stdout,70) i, TRIM(myVarName), TRIM(ncname), &
     &                                TRIM(SourceFile)
                    exit_flag=2
                    ioerror=status
                    EXIT
                  END IF
                ELSE
                  WRITE (stdout,80) i, TRIM(myVarName), TRIM(ncname),   &
     &                              TRIM(SourceFile)
                  exit_flag=4
                  ioerror=status
                  EXIT
                END IF
              END DO
            END IF
          END IF

#ifdef DISTRIBUTE
!
!  Broadcast requested variable information to all processors in the
!  group.
!
          IF (foundit) THEN
            CALL mp_bcasti (ng, model, exit_flag)
            IF (exit_flag.eq.NoError) THEN
              CALL mp_bcasts (ng, model, var_Dname)
              IF (n_vdim.gt.0) THEN
                CALL mp_bcasti (ng, model, var_Dsize)
              END IF
              CALL mp_bcasts (ng, model, var_Aname)
              CALL mp_bcasts (ng, model, var_Achar)
              IF (n_vatt.gt.0) THEN
                CALL mp_bcasti (ng, model, var_Aint)
                CALL mp_bcastf (ng, model, var_Afloat)
              END IF
            END IF
          END IF
#endif
!
!  Ignore error message if requested variable not found when searching
!  over multiple input NetCDF files.
!
          IF (PRESENT(SearchVar)) THEN
            SearchVar=foundit
            WriteError=.FALSE.
          ELSE
            WriteError=.TRUE.
          END IF
          IF (.not.foundit.and.WriteError) THEN
            IF (Master) WRITE (stdout,90) TRIM(myVarName),              &
     &                                    TRIM(ncname),                 &
     &                                    TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF
!
!  Close input NetCDF file.
!
      IF (InpThread) THEN
        status=nf90_close(ncid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,100) TRIM(ncname), TRIM(SourceFile)
          exit_flag=2
          ioerror=status
        END IF
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag)
#endif
!
  10  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring variable ID:', &
     &        2x,i3,/,18x,'in input file:',2x,a,/,18x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_INQ_VAR - too small dimension parameter,',a,   &
     &        2i5,/,18x,'change file  nc_util.F  and recompile')
  30  FORMAT (/,' NETCDF_INQ_VAR - unable to inquire about contents',   &
     &          ' of input NetCDF file:',2x,a,/,18x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_INQ_VAR - unable to open input NetCDF file:',  &
     &        2x,a,/,18x,'call from:',2x,a)
  50  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring dimension ',   &
     &        i1,' for variable:',2x,a,/,18x,'in input file:',2x,a,     &
     &        /,18x,'call from:',2x,a)
  60  FORMAT (/,' NETCDF_INQ_VAR - error while reading ',a,             &
     &        'attribute:',1x,a,' for variable ',a,/,18x,               &
     &        'in input file:',2x,a,/,18x,'call from:',2x,a)
  70  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring attribute ',   &
     &        i1,' for variable: ',a,/,18x,'in input file:',2x,a,       &
     &        /,18x,'call from:',2x,a)
  80  FORMAT (/,' NETCDF_INQ_VAR - unable to inquire name of ',         &
     &        'attribute ',i1,' for variable ',a,/,18x,                 &
     &        'in input file:',2x,a,/,18x,'call from:',2x,a)
  90  FORMAT (/,' NETCDF_INQ_VAR - requested variable:',2x,a,/18x,      &
     &        'not found in input file:',2x,a,/,18x,'call from:',2x,a)
 100  FORMAT (/,' NETCDF_INQ_VAR - unable to close input NetCDF file:', &
     &        2x,a,/,18x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_inq_var

      SUBROUTINE netcdf_inq_varid (ng, model, ncname, myVarName,        &
     &                             ncid, VarID)
!
!=======================================================================
!                                                                      !
!  This routine inquires ID of requested NetCDF variable.              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Requested variable name (string)                    !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     VarID        ID of requested variable (integer)                  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      integer, intent(out) :: VarID
!
!  Local variable declarations.
!
      integer :: status
#ifdef DISTRIBUTE
      integer, dimension(3) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Inquire ID of requested variable.
!-----------------------------------------------------------------------
!
      IF (InpThread) THEN
        status=nf90_inq_varid(ncid, TRIM(myVarName), VarID)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(myVarName), TRIM(ncname),              &
     &                      TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
#ifdef DISTRIBUTE
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      ibuffer(3)=VarID
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      VarID=ibuffer(3)
#endif
!
  10  FORMAT (/,' NETCDF_INQ_VARID - error while inquiring ID for ',    &
     &        'variable:',2x,a,/,20x,'in input file:',2x,a,/,           &
     &        20x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_inq_varid

      SUBROUTINE netcdf_get_fvar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point scalar variable from    !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (real)                         !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
      real(r8), dimension(1) :: my_A
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (InpThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=4
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  Read in variable.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
          IF (status.eq.nf90_noerr) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=nf90_get_var(my_ncid, varid, my_A, start, total)
              A=my_A(1)
            ELSE
              status=nf90_get_var(my_ncid, varid, A)
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Compute minimum and maximun values of read variable. Notice that
!  the same read value is assigned since a scalar variable was
!  processed.
!
      IF (PRESENT(min_val)) THEN
        min_val=A
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=A
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_0D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_FVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_FVAR_0D - unable to close NetCDF file:',   &
     &         2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_fvar_0d

      SUBROUTINE netcdf_get_fvar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 1D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:))   !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(:,1))  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

      integer, dimension(1) :: Asize
#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point 1D-array variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)

      IF (InpThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=4
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  Read in variable.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
          IF (status.eq.nf90_noerr) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=nf90_get_var(my_ncid, varid, A, start, total)
            ELSE
              status=nf90_get_var(my_ncid, varid, A)
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Compute minimum and maximun values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_1D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_FVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_FVAR_1D - unable to close NetCDF file:',   &
     &         2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_fvar_1d

      SUBROUTINE netcdf_get_fvar_2d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 2D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:,:)) !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(0:,0:))!
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar(:,:,1))!
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:,:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

      integer, dimension(2) :: Asize
#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)

      IF (InpThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=4
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  Read in variable.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
          IF (status.eq.nf90_noerr) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=nf90_get_var(my_ncid, varid, A, start, total)
            ELSE
              status=nf90_get_var(my_ncid, varid, A)
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname)
            exit_flag=2
            ioerror=status
          END IF
        END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Compute minimum and maximun values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_2D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_2D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_FVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_FVAR_2D - unable to close NetCDF file:',   &
     &         2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_fvar_2d

      SUBROUTINE netcdf_get_fvar_3d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 3D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 3D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:,:,:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

      integer, dimension(3) :: Asize
#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Asize(3)=UBOUND(A, DIM=3)

      IF (InpThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=4
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  Read in variable.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
          IF (status.eq.nf90_noerr) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=nf90_get_var(my_ncid, varid, A, start, total)
            ELSE
              status=nf90_get_var(my_ncid, varid, A)
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Compute minimum and maximun values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_3D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_3D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_FVAR_3D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_FVAR_3D - unable to close NetCDF file:',   &
     &         2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_fvar_3d

      SUBROUTINE netcdf_get_fvar_4d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total,                &
     &                               min_val, max_val)
!
!=======================================================================
!                                                                      !
!  This routine reads requested floating-point 4D-array variable from  !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 4D-array variable (real)                       !
!     min_val      Read data minimum value (real, OPTIONAL)            !
!     max_val      Read data maximum value (real, OPTIONAL)            !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_fvar (ng, iNLM, 'file.nc', 'VarName', fvar)       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      real(r8), intent(out), optional :: min_val
      real(r8), intent(out), optional :: max_val

      real(r8), intent(out) :: A(:,:,:,:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

      integer, dimension(4) :: Asize
#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point 2D-array variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)
      Asize(3)=UBOUND(A, DIM=3)
      Asize(4)=UBOUND(A, DIM=4)

      IF (InpThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=4
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  Read in variable.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
          IF (status.eq.nf90_noerr) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=nf90_get_var(my_ncid, varid, A, start, total)
            ELSE
              status=nf90_get_var(my_ncid, varid, A)
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, A)
      END IF
#endif
!
!  Compute minimum and maximun values of read variable.
!
      IF (PRESENT(min_val)) THEN
        min_val=MINVAL(A)
      END IF
      IF (PRESENT(max_val)) THEN
        max_val=MAXVAL(A)
      END IF
!
  10  FORMAT (/,' NETCDF_GET_FVAR_4D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_FVAR_4D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_FVAR_4D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_FVAR_4D - unable to close NetCDF file:',   &
     &         2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_fvar_4d

      SUBROUTINE netcdf_get_ivar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer scalar variable from specified !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (integer)                      !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      integer, intent(out) :: A
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

      integer, dimension(1) :: my_A
#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in an integer scalar variable.
!-----------------------------------------------------------------------
!
      IF (InpThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=4
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  Read in variable.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
          IF (status.eq.nf90_noerr) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=nf90_get_var(my_ncid, varid, my_A, start, total)
              A=my_A(1)
            ELSE
              status=nf90_get_var(my_ncid, varid, A)
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, A)
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_IVAR_0D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_IVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_IVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_IVAR_0D - unable to close NetCDF file:',   &
     &         2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_ivar_0d

      SUBROUTINE netcdf_get_ivar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer 1D-array variable from         !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (integer)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(0:))   !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(:,1))  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      integer, intent(out) :: A(:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

      integer, dimension(1) :: Asize
#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in an integer 1D-array variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)

      IF (InpThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=4
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  Read in variable.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
          IF (status.eq.nf90_noerr) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=nf90_get_var(my_ncid, varid, A, start, total)
            ELSE
              status=nf90_get_var(my_ncid, varid, A)
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, A)
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_IVAR_1D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_IVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_IVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_IVAR_1D - unable to close NetCDF file:',   &
     &         2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_ivar_1d

      SUBROUTINE netcdf_get_ivar_2d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested integer 2D-array variable from         !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 2D-array variable (integer)                    !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(0:,:)) !
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(0:,0:))!
!    CALL netcdf_get_ivar (ng, iNLM, 'file.nc', 'VarName', ivar(:,:,1))!
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      integer, intent(out) :: A(:,:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

      integer, dimension(2) :: Asize
#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in an integer 2D-array variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)
      Asize(2)=UBOUND(A, DIM=2)

      IF (InpThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=4
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  Read in variable.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
          IF (status.eq.nf90_noerr) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=nf90_get_var(my_ncid, varid, A, start, total)
            ELSE
              status=nf90_get_var(my_ncid, varid, A)
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, A)
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_IVAR_2D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_IVAR_2D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_IVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_IVAR_2D - unable to close NetCDF file:',   &
     &         2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_ivar_2d

      SUBROUTINE netcdf_get_svar_0d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string scalar variable from specified  !
!  NetCDF file.                                                        !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read scalar variable (string)                       !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar(1))    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      character (len=*), intent(out) :: A
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a string scalar variable.
!-----------------------------------------------------------------------
!
      IF (InpThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=4
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  Read in variable.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
          IF (status.eq.nf90_noerr) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=nf90_get_var(my_ncid, varid, A, start, total)
            ELSE
              status=nf90_get_var(my_ncid, varid, A)
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasts (ng, model, A)
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_SVAR_0D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_SVAR_0D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_SVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_SVAR_0D - unable to close NetCDF file:',   &
     &         2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_svar_0d

      SUBROUTINE netcdf_get_svar_1d (ng, model, ncname, myVarName, A,   &
     &                               ncid, start, total)
!
!=======================================================================
!                                                                      !
!  This routine reads requested string 1D-array variable from          !
!  specified NetCDF file.                                              !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     start        Starting index where the first of the data values   !
!                    will be read along each dimension (integer,       !
!                    OPTIONAL)                                         !
!     total        Number of data values to be read along each         !
!                    dimension (integer, OPTIONAL)                     !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     A            Read 1D-array variable (string)                     !
!                                                                      !
!  Examples:                                                           !
!                                                                      !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar)       !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar(0:))   !
!    CALL netcdf_get_svar (ng, iNLM, 'file.nc', 'VarName', ivar(:,1))  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti, mp_bcasts
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in), optional :: ncid
      integer, intent(in), optional :: start(:)
      integer, intent(in), optional :: total(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName

      character (len=*), intent(out) :: A(:)
!
!  Local variable declarations.
!
      integer :: my_ncid, status, varid

      integer, dimension(1) :: Asize
#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a string 1D-array variable.
!-----------------------------------------------------------------------
!
      Asize(1)=UBOUND(A, DIM=1)

      IF (InpThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for reading.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_nowrite, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=4
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  Read in variable.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_inq_varid(my_ncid, TRIM(myVarName), varid)
          IF (status.eq.nf90_noerr) THEN
            IF (PRESENT(start).and.PRESENT(total)) THEN
              status=nf90_get_var(my_ncid, varid, A, start, total)
            ELSE
              status=nf90_get_var(my_ncid, varid, A)
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=2
              ioerror=status
            END IF
          ELSE
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=2
            ioerror=status
          END IF
        END IF
!
!  If NetCDF file ID is not provided, close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname)
            exit_flag=2
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast read variable to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasts (ng, model, A)
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_SVAR_1D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_GET_SVAR_1D - error while reading variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  30  FORMAT (/,' NETCDF_GET_SVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  40  FORMAT (/,' NETCDF_GET_SVAR_1D - unable to close NetCDF file:',   &
     &         2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_get_svar_1d

      SUBROUTINE netcdf_put_fvar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point scalar variable into a NetCDF  !
!  file.  If the NetCDF ID is not provided, it opens the file, writes  !
!  data, and then closes the file.                                     !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE netcdf

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      real(r8), intent(in) :: A

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
      real(r8), dimension(1) :: my_A
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_write, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  If variable ID is not provided, inquire its value.
!
        IF (exit_flag.eq.NoError) THEN
          IF (.not.PRESENT(varid)) THEN
            status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            my_varid=varid
          END IF
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
            status=nf90_put_var(my_ncid, my_varid, A)
          ELSE
            my_A(1)=A
            status=nf90_put_var(my_ncid, my_varid, my_A, start, total)
          END IF
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_0D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_PUT_FVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_PUT_FVAR_0D - unable to close NetCDF file:',   &
     &        2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_put_fvar_0d

      SUBROUTINE netcdf_put_fvar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 1D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE netcdf

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      real(r8), intent(in) :: A(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_write, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  If variable ID is not provided, inquire its value.
!
        IF (exit_flag.eq.NoError) THEN
          IF (.not.PRESENT(varid)) THEN
            status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            my_varid=varid
          END IF
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_1D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_PUT_FVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_PUT_FVAR_1D - unable to close NetCDF file:',   &
     &        2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_put_fvar_1d

      SUBROUTINE netcdf_put_fvar_2d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 2D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE netcdf

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      real(r8), intent(in) :: A(:,:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_write, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  If variable ID is not provided, inquire its value.
!
        IF (exit_flag.eq.NoError) THEN
          IF (.not.PRESENT(varid)) THEN
            status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            my_varid=varid
          END IF
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_2D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_PUT_FVAR_2D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_PUT_FVAR_2D - unable to close NetCDF file:',   &
     &        2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_put_fvar_2d

      SUBROUTINE netcdf_put_fvar_3d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 3D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE netcdf

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      real(r8), intent(in) :: A(:,:,:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_write, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  If variable ID is not provided, inquire its value.
!
        IF (exit_flag.eq.NoError) THEN
          IF (.not.PRESENT(varid)) THEN
            status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            my_varid=varid
          END IF
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_3D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_3D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_PUT_FVAR_3D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_PUT_FVAR_3D - unable to close NetCDF file:',   &
     &        2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_put_fvar_3d

      SUBROUTINE netcdf_put_fvar_4d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a floating-point 4D-array variable into a file. !
!  If the NetCDF ID is not provided,  it opens the file,  writes data, !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (real)                  !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE netcdf

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      real(r8), intent(in) :: A(:,:,:,:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_write, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  If variable ID is not provided, inquire its value.
!
        IF (exit_flag.eq.NoError) THEN
          IF (.not.PRESENT(varid)) THEN
            status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            my_varid=varid
          END IF
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_4D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_4D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_PUT_FVAR_4D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_PUT_FVAR_4D - unable to close NetCDF file:',   &
     &        2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_put_fvar_4d

      SUBROUTINE netcdf_put_ivar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a integer scalar variable into a NetCDF file.   !
!  If the NetCDF ID is not provided, it opens the file, writes data,   !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE netcdf

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      integer, intent(in) :: A

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

      integer, dimension(1) :: my_A
#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_write, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  If variable ID is not provided, inquire its value.
!
        IF (exit_flag.eq.NoError) THEN
          IF (.not.PRESENT(varid)) THEN
            status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            my_varid=varid
          END IF
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          IF ((start(1).eq.0).and.(total(1).eq.0)) THEN
            status=nf90_put_var(my_ncid, my_varid, A)
          ELSE
            my_A(1)=A
            status=nf90_put_var(my_ncid, my_varid, my_A, start, total)
          END IF
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_0D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_PUT_IVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_PUT_IVAR_0D - unable to close NetCDF file:',   &
     &        2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_put_ivar_0d

      SUBROUTINE netcdf_put_ivar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a integer 1D-array variable into a file. If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE netcdf

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      integer, intent(in) :: A(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_write, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  If variable ID is not provided, inquire its value.
!
        IF (exit_flag.eq.NoError) THEN
          IF (.not.PRESENT(varid)) THEN
            status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            my_varid=varid
          END IF
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_1D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_PUT_IVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_PUT_IVAR_1D - unable to close NetCDF file:',   &
     &        2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_put_ivar_1d

      SUBROUTINE netcdf_put_ivar_2d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a integer 2D-array variable into a file. If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (integer)               !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars
      USE netcdf

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      integer, intent(in) :: A(:,:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_write, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  If variable ID is not provided, inquire its value.
!
        IF (exit_flag.eq.NoError) THEN
          IF (.not.PRESENT(varid)) THEN
            status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            my_varid=varid
          END IF
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_put_var(my_ncid, my_varid, A, start, total)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_2D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_2D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_PUT_IVAR_2D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_PUT_IVAR_2D - unable to close NetCDF file:',   &
     &        2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_put_ivar_2d

      SUBROUTINE netcdf_put_svar_0d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a  string scalar variable into a NetCDF file.   !
!  If the NetCDF ID is not provided, it opens the file, writes data,   !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (string)                !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(in) :: start(:), total(:)    ! dummy variables

      integer, intent(in), optional :: ncid, varid

      character (len=*), intent(in) :: A
      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

      integer, dimension(1) :: mystart, mytotal
#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_write, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  If variable ID is not provided, inquire its value.
!
        IF (exit_flag.eq.NoError) THEN
          IF (.not.PRESENT(varid)) THEN
            status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            my_varid=varid
          END IF
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          mystart(1)=1
          mytotal(1)=LEN_TRIM(A)
          status=nf90_put_var(my_ncid, my_varid, A, mystart, mytotal)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),             &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_0D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_0D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_PUT_SVAR_0D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_PUT_SVAR_0D - unable to close NetCDF file:',   &
     &        2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_put_svar_0d

      SUBROUTINE netcdf_put_svar_1d (ng, model, ncname, myVarName, A,   &
     &                               start, total, ncid, varid)
!
!=======================================================================
!                                                                      !
!  This routine writes a string 1D-array variable into a file.  If     !
!  the NetCDF ID is not provided, it opens the file,  writes data,     !
!  and then closes the file.                                           !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       NetCDF file name (string)                           !
!     myVarName    Variable name (string)                              !
!     A            Data value(s) to be written (string)                !
!     start        Starting index where the first of the data values   !
!                    will be written along each dimension (integer)    !
!     total        Number of data values to be written along each      !
!                    dimension (integer)                               !
!     ncid         NetCDF file ID (integer, OPTIONAL)                  !
!     varid        NetCDF variable ID (integer, OPTIONAL)              !
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     exit_flag    Error flag (integer) stored in MOD_SCALARS          !
!     ioerror      NetCDF return code (integer) stored in MOD_IOUNITS  !
!                                                                      !
!  Notice: This routine must be used to write only nontiled variables. !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
 
      integer, intent(in) :: start(:), total(:)

      integer, intent(in), optional :: ncid, varid

      character (len=*), intent(in) :: A(:)

      character (len=*), intent(in) :: ncname
      character (len=*), intent(in) :: myVarName
!
!  Local variable declarations.
!
      integer :: my_ncid, my_varid, status

      integer, dimension(2) :: mystart, mytotal

#ifdef DISTRIBUTE
      integer, dimension(2) :: ibuffer
#endif
!
!-----------------------------------------------------------------------
!  Read in a floating-point scalar variable.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
!
!  If NetCDF file ID is not provided, open NetCDF for writing.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_open(TRIM(ncname), nf90_write, my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        ELSE
          my_ncid=ncid
        END IF
!
!  If variable ID is not provided, inquire its value.
!
        IF (exit_flag.eq.NoError) THEN
          IF (.not.PRESENT(varid)) THEN
            status=nf90_inq_varid(my_ncid, TRIM(myVarName), my_varid)
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname),          &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          ELSE
            my_varid=varid
          END IF
        END IF
!
!  Write out data.
!
        IF (exit_flag.eq.NoError) THEN
          mystart(1)=1
          mytotal(1)=LEN_TRIM(A(1))
          mystart(2)=start(1)
          mytotal(2)=total(1)         
          status=nf90_put_var(my_ncid, my_varid, A, mystart, mytotal)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname),            &
     &                        TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
!
!  Close input NetCDF file.
!
        IF (.not.PRESENT(ncid)) THEN
          status=nf90_close(my_ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,40) TRIM(ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_1D - unable to open NetCDF file:',    &
     &        2x,a,/,22x,'call from:',2x,a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_1D - error while inquiring ID for ',  &
     &        'variable:',2x,a,/,22x,'in input file:',2x,a,/,22x,       &
     &        'call from:',2x,a)
  30  FORMAT (/,' NETCDF_PUT_SVAR_1D - error while writing variable:',  &
     &        2x,a,/,22x,'in input file:',2x,a,/,22x,'call from:',2x,a)
  40  FORMAT (/,' NETCDF_PUT_SVAR_1D - unable to close NetCDF file:',   &
     &        2x,a,/,22x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_put_svar_1d

      SUBROUTINE netcdf_close (ng, model, ncid, ncname)
!
!=======================================================================
!                                                                      !
!  This routine closes requested NetCDF file. If appropriate, it       !
!  also performs additional processing, like updating the global       !
!  attributes, before closing the file.                                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncid         NetCDF file ID (integer)                            !
!     ncname       NetCDF file name (string, OPTIONAL)                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(inout) :: ncid

      character (len=*), intent(in), optional :: ncname
!
!  Local variable declarations.
!
      integer :: i, status
#ifdef BIOLOGY
      integer :: AttLen, is, ie, lstr
#endif
#ifdef DISTRIBUTE
      integer, dimension(3) :: ibuffer
#endif
      character (len=200) :: my_ncname
#ifdef BIOLOGY
      character (len=512) :: bio_file
#endif
!
!-----------------------------------------------------------------------
!  If open, close requested NetCDF file.
!-----------------------------------------------------------------------
!
      IF (OutThread.and.(ncid.ne.-1)) THEN

        DO i=1,LEN(my_ncname)
          my_ncname(i:i)=' '
        END DO

        IF (.not.PRESENT(ncname)) THEN
!
!  Get filename, if any. It will be nice if there is a function in
!  the NetCDF library to do this. Fortunately, the filename is
!  written as a global attribute.
!
          status=nf90_get_att(ncid, nf90_global, 'file', my_ncname)
        ELSE
          my_ncname=TRIM(ncname)
        END IF

#ifdef BIOLOGY
!
!  Update global attribute with the biology header files used.
!
        is=1
        status=nf90_inquire_attribute(ncid, nf90_global, 'bio_file',    &
     &                                len = AttLen)
        IF (status.eq.nf90_noerr) THEN
          DO i=1,512
            bio_file(i:i)=' '
          END DO
          DO i=1,4
            lstr=LEN_TRIM(BIONAME(i))
            IF (lstr.gt.0) THEN
              ie=is+lstr-1
              bio_file(is:ie)=TRIM(BIONAME(i))
              is=ie+1
              bio_file(is:is)=','
              is=is+2
            END IF
          END DO
          lstr=LEN_TRIM(bio_file)-1
          IF (lstr.gt.0) THEN
            status=nf90_put_att(ncid, nf90_global, 'bio_file',          &
     &                          bio_file(1:lstr))
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,10) 'bio_file', TRIM(my_ncname),            &
     &                          TRIM(SourceFile)
              exit_flag=3
              ioerror=status
            END IF
          END IF
        END IF
#endif
!
!  Close requested NetCDF file.
!
        IF (exit_flag.eq.NoError) THEN
          status=nf90_close(ncid)
          IF (status.ne.nf90_noerr) THEN
            WRITE (stdout,20) ncid, TRIM(my_ncname), TRIM(SourceFile)
            exit_flag=3
            ioerror=status
          END IF
          ncid=-1
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast error flags to all processors in the group.
!
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      ibuffer(3)=ncid
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      ncid=ibuffer(3)
#endif
!
  10  FORMAT (/,' NETCDF_CLOSE - error while writing global ',          &
     &        'attribute:',2x,a,/,16x,'file:',2x,a,/,16x,               &
     &        'call from:'2x,a)
  20  FORMAT (/,' NETCDF_CLOSE - error during closing of file, ',       &
     &        'ncid = ',i3,/,16x,'file:',2x,a,/,16x,'call from:'2x,a)

      RETURN
      END SUBROUTINE netcdf_close

      SUBROUTINE netcdf_create (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine creates a new NetCDF file. Currently, it only creates  !
!  file for serial I/O access.  It will be modified in the future for  !
!  parallel I/O access.                                                !
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model

      integer, intent(out) :: ncid

      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: status
#ifdef DISTRIBUTE
      integer :: ibuffer(3)
#endif
!
!-----------------------------------------------------------------------
!  Create requested NetCDF file.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_create(TRIM(ncname), CMODE, ncid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
#ifdef DISTRIBUTE
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      ibuffer(3)=ncid
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      ncid=ibuffer(3)
#endif
!
  10  FORMAT (/,' NETCDF_CREATE - unable to create output NetCDF ',     &
     &        'file:',/,17x,a,/,17x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_create

      SUBROUTINE netcdf_open (ng, model, ncname, omode, ncid)
!
!=======================================================================
!                                                                      !
!  This routine opens an existing NetCDF file for access. Currently,   !
!  it only opens file for serial I/O access.  It will be modified in   !
!  the future for parallel I/O access.                                 !
!                                                                      !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     omode        Open mode flag:                                     !
!                    omode = 0,  read-only access (NF90_NOWRITE)       !
!                    omode = 1,  read and write access (NF90_WRITE)    !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, omode

      integer, intent(out) :: ncid

      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: status
#ifdef DISTRIBUTE
      integer :: ibuffer(3)
#endif
!
!-----------------------------------------------------------------------
!  Create requested NetCDF file.
!-----------------------------------------------------------------------
!
      IF (InpThread) THEN
        SELECT CASE (omode)
          CASE (0)
            status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
          CASE (1)
            status=nf90_open(TRIM(ncname), nf90_write, ncid)
          CASE DEFAULT
            status=nf90_open(TRIM(ncname), nf90_nowrite, ncid)
        END SELECT
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
#ifdef DISTRIBUTE
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      ibuffer(3)=ncid
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      ncid=ibuffer(3)
#endif
!
  10  FORMAT (/,' NETCDF_OPEN - unable to open existing NetCDF file:',  &
     &        /,15x,a,/,15x,'call from:'2x,a)

      RETURN
      END SUBROUTINE netcdf_open

      SUBROUTINE netcdf_openid (nfiles, fname, fid, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine determines association between pathnames and NetCDF    !
!  file descriptors. It is used to get the  ID  of an opened NetCDF    !
!  from file name.  This is done to avoid  opening  too  many files    !
!  simultaneously for the same dataset.                                !
!                                                                      !
!  Input:                                                              !
!                                                                      !
!    nfiles    Number of input NetCDF files.                           !
!    fname     NetCDF file name(s).                                    !
!    fid       NetCDF ID(s) associated with fname.                     !
!    ncname    NetCDF file name to check in association table.         !
!    ncid      NetCDF ID to check in association table.                !
!                                                                      !
!  Output:                                                             !
!                                                                      !
!    fid       If updated, it stores the ID associated with fname.     !
!    ncid      If updated, indicates that ncname has been opened       !
!                previously and this is its associated ID.             !
!                                                                      !
!=======================================================================
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: nfiles

      character (len=*), intent(in) :: fname(nfiles)
      character (len=*), intent(in) :: ncname

      integer, intent(inout) :: fid(nfiles)
      integer, intent(inout) :: ncid
!
!  Local variabel declarations.
!
      integer :: i
!
!-----------------------------------------------------------------------
!  Determine association between NetCDF file name and ID descriptor.
!-----------------------------------------------------------------------
!
!  NetCDF file ID descriptors are set to closed state (-1) during
!  initialization.
!
      DO i=1,nfiles
        IF (TRIM(ncname).eq.TRIM(fname(i))) THEN
!
!  NCFILE was just opened. Store its associated ID for future reference.
!
          IF (ncid.ne.-1) fid(i)=ncid
!
!  NCFILE is already open. Get its ID from reference table.
!
          IF ((fid(i).ne.-1).and.(ncid.eq.-1)) ncid=fid(i)
        END IF
      END DO

      RETURN
      END SUBROUTINE netcdf_openid

      SUBROUTINE netcdf_sync (ng, model, ncname, ncid)
!
!=======================================================================
!                                                                      !
!  This routine synchronize to disk requested NetCDF file with         !
!  in-memory buffer to make data available to other  processes         !
!  immediately after it is written.                                    !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng           Nested grid number (integer)                        !
!     model        Calling model identifier (integer)                  !
!     ncname       Name of the new NetCDF file (string)                !
!     ncid         NetCDF file ID (integer)                            !
!                                                                      !
!  NOTE:                                                               !
!                                                                      !
!  Nowadays, it is recommended to have the writer and readers open     !
!  the NetCDF file with the NF90_SHARE flag to improve performance.    !
!  See NetCDF user manual for more details.                            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_scalars

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcasti
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model, ncid

      character (len=*), intent(in) :: ncname
!
!  Local variable declarations.
!
      integer :: status
#ifdef DISTRIBUTE
      integer :: ibuffer(2)
#endif
!
!-----------------------------------------------------------------------
!  Create requested NetCDF file.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf90_sync(ncid)
        IF (status.ne.nf90_noerr) THEN
          WRITE (stdout,10) TRIM(ncname), TRIM(SourceFile)
          exit_flag=3
          ioerror=status
        END IF
      END IF
#ifdef DISTRIBUTE
      ibuffer(1)=exit_flag
      ibuffer(2)=ioerror
      CALL mp_bcasti (ng, model, ibuffer)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_SYNC - unable to synchronize to disk file:',   &
     &        /,15x,a,/,15x,'call from:',2x,a)

      RETURN
      END SUBROUTINE netcdf_sync

      END MODULE mod_netcdf
