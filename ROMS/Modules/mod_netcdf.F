#include "cppdefs.h"
      MODULE mod_netcdf
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2008 The ROMS/TOMS Group                         !
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
        MODULE PROCEDURE netcdf_get_fvar_s
        MODULE PROCEDURE netcdf_get_fvar_v
        MODULE PROCEDURE netcdf_get_fvar_m
      END INTERFACE netcdf_get_fvar

      INTERFACE netcdf_get_ivar
        MODULE PROCEDURE netcdf_get_ivar_s
        MODULE PROCEDURE netcdf_get_ivar_v
        MODULE PROCEDURE netcdf_get_ivar_m
      END INTERFACE netcdf_get_ivar

      INTERFACE netcdf_get_svar
        MODULE PROCEDURE netcdf_get_svar_s
        MODULE PROCEDURE netcdf_get_svar_v
      END INTERFACE netcdf_get_svar

      INTERFACE netcdf_put_fvar
        MODULE PROCEDURE netcdf_put_fvar_s
        MODULE PROCEDURE netcdf_put_fvar_v
        MODULE PROCEDURE netcdf_put_fvar_m
      END INTERFACE netcdf_put_fvar

      INTERFACE netcdf_put_ivar
        MODULE PROCEDURE netcdf_put_ivar_s
        MODULE PROCEDURE netcdf_put_ivar_v
        MODULE PROCEDURE netcdf_put_ivar_m
      END INTERFACE netcdf_put_ivar

      INTERFACE netcdf_put_svar
        MODULE PROCEDURE netcdf_put_svar_s
        MODULE PROCEDURE netcdf_put_svar_v
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
                WRITE (stdout,10) dim_id(i), TRIM(ncname)
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
              WRITE (stdout,20) TRIM(ncname)
              exit_flag=2
              ioerror=status
            END IF
          ELSE      
            IF (n_dim.gt.Mdims) THEN
              WRITE (stdout,30) ' Mdims = ', Mdims, n_dim
              exit_flag=2
            END IF
            IF (status.ne.nf90_noerr) THEN
              WRITE (stdout,40) TRIM(ncname)
              exit_flag=2
              ioerror=status
            END IF
          END IF
        ELSE
          WRITE (stdout,50) TRIM(ncname)
          ioerror=status
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast dimension to all processors in the group.
!
      CALL mp_bcasti (ng, model, exit_flag, 1)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, ibuffer, 5)
        n_dim=ibuffer(1)
        n_var=ibuffer(2)
        n_gatt=ibuffer(3)
        rec_id=ibuffer(4)
        rec_size=ibuffer(5)
        CALL mp_bcasti (ng, model, dim_id, Mdims)
        CALL mp_bcasti (ng, model, dim_size, Mdims)
        CALL mp_bcasts (ng, model, dim_name, n_dim)
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
     &        1x,i3,/,18x,'in input file: ',a)
  20  FORMAT (/,' NETCDF_GET_DIM - unable to close input NetCDF file: ',&
     &        a)
  30  FORMAT (/,' NETCDF_GET_DIM - too small dimension parameter,',a,   &
     &        2i5,/,18x,'change file  nc_util.F  and recompile')
  40  FORMAT (/,' NETCDF_GET_DIM - unable to inquire about contents',   &
     &          ' of input NetCDF file: ',a)
  50  FORMAT (/,' NETCDF_GET_DIM - unable to open input NetCDF file: ', &
     &        a)
  60  FORMAT (/,' NETCDF_GET_DIM - requested dimension: ',a,/18x,       &
     &        'not found in input file: ',a)

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
#if defined ADJUST_STFLUX || defined ADJUST_WSTRESS
          CASE ('Nadjust')
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
     &                           VarID, nVarDim, nVarAtt)
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
!                                                                      !
!  On Ouput:                                                           !
!                                                                      !
!     VarID        ID of requested dimension (integer, OPTIONAL)       !
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
!     var_Aname    Variable attributes name (1:n_vatt)                 !
!     var_Achar    Variable string attributes value (1:n_vatt)         !
!     var_Afloat   Variable float attributes value (1:n_vatt)          !
!     var_Aint     Variable integer attributes value (1:n_vatt)        !
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

      integer, intent(out), optional :: VarID
      integer, intent(out), optional :: nVarDim
      integer, intent(out), optional :: nVarAtt
!
!  Local variable declarations.
!
      logical :: foundit

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
                    WRITE (stdout,70) j, TRIM(var_name(i)), TRIM(ncname)
                    exit_flag=2
                    ioerror=status
                    EXIT
                  END IF
                END DO
              ELSE
                WRITE (stdout,10) var_id(i), TRIM(ncname)
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
              WRITE (stdout,30) TRIM(ncname)
              exit_flag=2
              ioerror=status
            END IF
          END IF
        ELSE
          WRITE (stdout,40) TRIM(ncname)
          exit_flag=4
          ioerror=status
        END IF
      END IF

#ifdef DISTRIBUTE
!
!  Broadcast dimension to all processors in the group.
!
      CALL mp_bcasti (ng, model, exit_flag, 1)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, ibuffer, 4)
        n_dim=ibuffer(1)
        n_var=ibuffer(2)
        n_gatt=ibuffer(3)
        rec_id=ibuffer(4)
        CALL mp_bcasti (ng, model, var_id, Mvars)
        CALL mp_bcasti (ng, model, var_flag, Mvars)
        CALL mp_bcasti (ng, model, var_type, Mvars)
        CALL mp_bcasti (ng, model, var_ndim, Mvars)
        CALL mp_bcasti (ng, model, var_natt, Mvars)
        CALL mp_bcasti (ng, model, var_dim, NvarD, Mdims)
        CALL mp_bcasts (ng, model, var_name, n_var)
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
                WRITE (stdout,50) i, TRIM(myVarName), TRIM(ncname)
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
     &                                    TRIM(myVarName), TRIM(ncname)
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
     &                                    TRIM(myVarName), TRIM(ncname)
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
     &                                    TRIM(myVarName), TRIM(ncname)
                        exit_flag=2
                        ioerror=status
                        EXIT
                      END IF
                      var_Achar(i)=text(1:my_Alen)
                    END IF
                  ELSE
                    WRITE (stdout,70) i, TRIM(myVarName), TRIM(ncname)
                    exit_flag=2
                    ioerror=status
                    EXIT
                  END IF
                ELSE
                  WRITE (stdout,80) i, TRIM(myVarName), TRIM(ncname)
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
            CALL mp_bcasti (ng, model, exit_flag, 1)
            IF (exit_flag.eq.NoError) THEN
              CALL mp_bcasts (ng, model, var_Dname, n_vdim)
              IF (n_vdim.gt.0) THEN
                CALL mp_bcasti (ng, model, var_Dsize, NvarD)
              END IF
              CALL mp_bcasts (ng, model, var_Aname, n_vatt)
              CALL mp_bcasts (ng, model, var_Achar, n_vatt)
              IF (n_vatt.gt.0) THEN
                CALL mp_bcasti (ng, model, var_Aint, NvarA)
                CALL mp_bcastf (ng, model, var_Afloat, NvarA)
              END IF
            END IF
          END IF
#endif
          IF (.not.foundit) THEN
            IF (Master) WRITE (stdout,90) TRIM(myVarName), TRIM(ncname)
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
          WRITE (stdout,100) TRIM(ncname)
          exit_flag=2
          ioerror=status
        END IF
      END IF
#ifdef DISTRIBUTE
      CALL mp_bcasti (ng, model, exit_flag, 1)
#endif
!
  10  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring variable ID:', &
     &        1x,i3,/,18x,'in input file: ',a)
  20  FORMAT (/,' NETCDF_INQ_VAR - too small dimension parameter,',a,   &
     &        2i5,/,18x,'change file  nc_util.F  and recompile')
  30  FORMAT (/,' NETCDF_INQ_VAR - unable to inquire about contents',   &
     &          ' of input NetCDF file: ',a)
  40  FORMAT (/,' NETCDF_INQ_VAR - unable to open input NetCDF file: ', &
     &        a)
  50  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring dimension ',   &
     &        i1,' for variable: ',a,/,18x,'in input file: ',a)
  60  FORMAT (/,' NETCDF_INQ_VAR - error while reading ',a,             &
     &        'attribute: ',a,' for variable ',a,/,18x,                 &
     &        'in input file: ',a)
  70  FORMAT (/,' NETCDF_INQ_VAR - error while inquiring attribute ',   &
     &        i1,' for variable: ',a,/,18x,'in input file: ',a)
  80  FORMAT (/,' NETCDF_INQ_VAR - unable to inquire name of ',         &
     &        'attribute ',i1,' for variable ',a,/,18x,                 &
     &        'in input file: ',a)
  90  FORMAT (/,' NETCDF_INQ_VAR - requested variable: ',a,/18x,        &
     &        'not found in input file: ',a)
 100  FORMAT (/,' NETCDF_INQ_VAR - unable to close input NetCDF file:', &
     &        1x,a)

      RETURN
      END SUBROUTINE netcdf_inq_var

      SUBROUTINE netcdf_get_fvar_s (ng, model, ncname, myVarName, A,    &
     &                              ncid, start, total)
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, A, 1)
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_FVAR_S - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_GET_FVAR_S - error while reading variable: ',  &
     &        a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_GET_FVAR_S - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_GET_FVAR_S - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_get_fvar_s

      SUBROUTINE netcdf_get_fvar_v (ng, model, ncname, myVarName, A,    &
     &                              ncid, start, total)
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, A, Asize(1))
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_FVAR_V - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_GET_FVAR_V - error while reading variable: ',  &
     &        a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_GET_FVAR_V - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_GET_FVAR_V - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_get_fvar_v

      SUBROUTINE netcdf_get_fvar_m (ng, model, ncname, myVarName, A,    &
     &                              ncid, start, total)
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcastf (ng, model, A, Asize(1), Asize(2))
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_FVAR_M - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_GET_FVAR_M - error while reading variable: ',  &
     &        a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_GET_FVAR_M - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_GET_FVAR_M - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_get_fvar_m

      SUBROUTINE netcdf_get_ivar_s (ng, model, ncname, myVarName, A,    &
     &                              ncid, start, total)
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, A, 1)
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_IVAR_S - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_GET_IVAR_S - error while reading variable: ',  &
     &        a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_GET_IVAR_S - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_GET_IVAR_S - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_get_ivar_s

      SUBROUTINE netcdf_get_ivar_v (ng, model, ncname, myVarName, A,    &
     &                              ncid, start, total)
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, A, Asize(1))
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_IVAR_V - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_GET_IVAR_V - error while reading variable: ',  &
     &        a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_GET_IVAR_V - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_GET_IVAR_V - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_get_ivar_v

      SUBROUTINE netcdf_get_ivar_m (ng, model, ncname, myVarName, A,    &
     &                              ncid, start, total)
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasti (ng, model, A, Asize(1), Asize(2))
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_IVAR_M - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_GET_IVAR_M - error while reading variable: ',  &
     &        a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_GET_IVAR_M - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_GET_IVAR_M - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_get_ivar_m

      SUBROUTINE netcdf_get_svar_s (ng, model, ncname, myVarName, A,    &
     &                              ncid, start, total)
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasts (ng, model, A, 1)
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_SVAR_S - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_GET_SVAR_S - error while reading variable: ',  &
     &        a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_GET_SVAR_S - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_GET_SVAR_S - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_get_svar_s

      SUBROUTINE netcdf_get_svar_v (ng, model, ncname, myVarName, A,    &
     &                              ncid, start, total)
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
      IF (exit_flag.eq.NoError) THEN
        CALL mp_bcasts (ng, model, A, Asize(1))
      END IF
#endif
!
  10  FORMAT (/,' NETCDF_GET_SVAR_V - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_GET_SVAR_V - error while reading variable: ',  &
     &        a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_GET_SVAR_V - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_GET_SVAR_V - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_get_svar_v

      SUBROUTINE netcdf_put_fvar_s (ng, model, ncname, myVarName, A,    &
     &                              start, total, ncid, varid)
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
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,40) TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_S - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_S - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_PUT_FVAR_S - error while writing variable: ',  &
     &        a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_PUT_FVAR_S - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_put_fvar_s

      SUBROUTINE netcdf_put_fvar_v (ng, model, ncname, myVarName, A,    &
     &                              start, total, ncid, varid)
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
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,40) TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_V - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_V - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_PUT_FVAR_V - error while writing variable: ',  &
     &        a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_PUT_FVAR_V - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_put_fvar_v

      SUBROUTINE netcdf_put_fvar_m (ng, model, ncname, myVarName, A,    &
     &                              start, total, ncid, varid)
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
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,40) TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_FVAR_M - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_PUT_FVAR_M - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_PUT_FVAR_M - error while writing variable: ',  &
     &        a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_PUT_FVAR_M - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_put_fvar_m

      SUBROUTINE netcdf_put_ivar_s (ng, model, ncname, myVarName, A,    &
     &                              start, total, ncid, varid)
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
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,40) TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_S - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_S - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_PUT_IVAR_S - error while writing variable: ',  &
     &        a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_PUT_IVAR_S - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_put_ivar_s

      SUBROUTINE netcdf_put_ivar_v (ng, model, ncname, myVarName, A,    &
     &                              start, total, ncid, varid)
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
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,40) TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_V - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_V - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_PUT_IVAR_V - error while writing variable: ',  &
     &        a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_PUT_IVAR_V - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_put_ivar_v

      SUBROUTINE netcdf_put_ivar_m (ng, model, ncname, myVarName, A,    &
     &                              start, total, ncid, varid)
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
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,40) TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_IVAR_M - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_PUT_IVAR_M - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_PUT_IVAR_M - error while writing variable: ',  &
     &        a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_PUT_IVAR_M - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_put_ivar_m

      SUBROUTINE netcdf_put_svar_s (ng, model, ncname, myVarName, A,    &
     &                              start, total, ncid, varid)
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
      USE netcdf

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,40) TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_S - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_S - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_PUT_SVAR_S - error while writing variable: ',  &
     &        a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_PUT_SVAR_S - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_put_svar_s

      SUBROUTINE netcdf_put_svar_v (ng, model, ncname, myVarName, A,    &
     &                              start, total, ncid, varid)
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
      USE netcdf

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_bcastf, mp_bcasti
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
            WRITE (stdout,10) TRIM(ncname)
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
              WRITE (stdout,20) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,30) TRIM(myVarName), TRIM(ncname)
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
            WRITE (stdout,40) TRIM(ncname)
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
      CALL mp_bcasti (ng, model, ibuffer, 2)
      exit_flag=ibuffer(1)
      ioerror=ibuffer(2)
#endif
!
  10  FORMAT (/,' NETCDF_PUT_SVAR_V - unable to open NetCDF file: ',a)
  20  FORMAT (/,' NETCDF_PUT_SVAR_V - error while inquiring ID for ',   &
     &        'variable: ',a,/,21x,'in input file: ',a)
  30  FORMAT (/,' NETCDF_PUT_SVAR_V - error while writing variable: ',  &
     &        a,/,21x,'in input file: ',a)
  40  FORMAT (/,' NETCDF_PUT_SVAR_V - unable to close NetCDF file: ',a)

      RETURN
      END SUBROUTINE netcdf_put_svar_v

      END MODULE mod_netcdf
