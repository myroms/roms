#include "cppdefs.h"
#if defined PROPAGATOR && defined CHECKPOINTING
      SUBROUTINE wrt_gst (ng, model)
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine writes checkpointing fields into GST restart NetCDF    !
!  file.                                                               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_ncparam
      USE mod_netcdf
      USE mod_scalars
      USE mod_storage

#ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_ncwrite
#endif
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
!
!  Local variable declarations.
!
      integer :: i, status, varid

# ifdef DISTRIBUTE
      integer :: vrecord = -1

      real(r8) :: scale = 1.0_r8

      character (len=6) :: var
# else
      integer, dimension(2) :: start, total
# endif
      character (len=1) :: lchar(5)
!
!-----------------------------------------------------------------------
!  Write out checkpointing information variables.
!-----------------------------------------------------------------------
!
      IF (exit_flag.ne.NoError) RETURN
!
!  Write out number of eigenvalues to compute.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'NEV', varid)
        status=nf_put_var1_int(ncGSTid(ng), varid, 1, NEV)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'NEV', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out number of Lanczos vectors to compute.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'NCV', varid)
        status=nf_put_var1_int(ncGSTid(ng), varid, 1, NCV)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'NCV', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out size of the eigenvalue problem.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'Mstate', varid)
        status=nf_put_var1_int(ncGSTid(ng), varid, 1, Mstate(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Mstate', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF

# ifdef DISTRIBUTE
!
!  Write out number of distributed-memory nodes.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'Nnodes', varid)
        status=nf_put_var1_int(ncGSTid(ng), varid, 1, numnodes)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Nnodes', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
# endif
!
!  Write out iteration number.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'iter', varid)
!!      status=nf_put_var1_int(ncGSTid(ng), varid, 1, Nrun+1)
        status=nf_put_var1_int(ncGSTid(ng), varid, 1, Nrun)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'iter', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out reverse communications flag.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'ido', varid)
        status=nf_put_var1_int(ncGSTid(ng), varid, 1, ido)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ido', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out information and error flag.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'info', varid)
        status=nf_put_var1_int(ncGSTid(ng), varid, 1, info)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'info', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out eigenvalue problem type.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'bmat', varid)
        status=nf_put_var1_text(ncGSTid(ng), varid, 1, bmat)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'bmat', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out Ritz eigenvalues to compute.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'which', varid)
        status=nf_put_vara_text(ncGSTid(ng), varid, 1, 2, which)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'which', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out form of basis function.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'howmany', varid)
        status=nf_put_var1_text(ncGSTid(ng), varid, 1, howmany)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'howmay', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out relative accuracy of computed Ritz values.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'Ritz_tol', varid)
        status=nf_put_var1_TYPE(ncGSTid(ng), varid, 1, Ritz_tol)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'Ritz_tol', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out eigenproblem parameters.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'iparam', varid)
        status=nf_put_vara_int(ncGSTid(ng), varid, 1, SIZE(iparam),     &
     &                         iparam)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'iparam', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write out pointers to mark starting location in work arrays.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'ipntr', varid)
        status=nf_put_vara_int(ncGSTid(ng), varid, 1, SIZE(ipntr),      &
     &                         ipntr)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'ipntr', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write ARPACK internal integer parameters to _aupd routines.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'iaupd', varid)
        status=nf_put_vara_int(ncGSTid(ng), varid, 1, SIZE(iaupd),      &
     &                         iaupd)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'iaupd', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write ARPACK internal integer parameters to _aitr routines.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'iaitr', varid)
        status=nf_put_vara_int(ncGSTid(ng), varid, 1, SIZE(iaitr),      &
     &                         iaitr)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'iaitr', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write ARPACK internal integer parameters to _aup2 routines.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'iaup2', varid)
        status=nf_put_vara_int(ncGSTid(ng), varid, 1, SIZE(iaup2),      &
     &                         iaup2)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'iaup2', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write ARPACK internal logical parameters to _aitr routines.
!
      IF (OutThread) THEN
        DO i=1,SIZE(laitr)
          IF (laitr(i)) THEN
            lchar(i)='T'
          ELSE
            lchar(i)='F'
          END IF
        END DO
        status=nf_inq_varid(ncGSTid(ng), 'laitr', varid)
        status=nf_put_vara_text(ncGSTid(ng), varid, 1, SIZE(laitr),     &
     &                          lchar)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'laitr', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Write ARPACK internal logical parameters to _aupd routines.
!
      IF (OutThread) THEN
        DO i=1,SIZE(laup2)
          IF (laup2(i)) THEN
            lchar(i)='T'
          ELSE
            lchar(i)='F'
          END IF
        END DO
        status=nf_inq_varid(ncGSTid(ng), 'laup2', varid)
        status=nf_put_vara_text(ncGSTid(ng), varid, 1, SIZE(laup2),     &
     &                          lchar)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'laup2', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Define ARPACK internal real parameters to _aitr routines.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'raitr', varid)
        status=nf_put_vara_TYPE(ncGSTid(ng), varid, 1, SIZE(raitr),     &
     &                          raitr)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'raitr', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!  Define ARPACK internal real parameters to _aup2 routines.
!
      IF (OutThread) THEN
        status=nf_inq_varid(ncGSTid(ng), 'raup2', varid)
        status=nf_put_vara_TYPE(ncGSTid(ng), varid, 1, SIZE(raup2),     &
     &                          raup2)
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,10) 'raup2', TRIM(GSTname(ng))
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
!
!-----------------------------------------------------------------------
!  Write out checkpointing variables associated with the state vector.
!-----------------------------------------------------------------------
!
!  Write out Lanczos/Arnoldi basis vectors.
!
# ifdef DISTRIBUTE
      var='Bvec'
      status=mp_ncwrite (ng, model, ncGSTid(ng), var, GSTname(ng),      &
     &                   vrecord, Nstr(ng), Nend(ng), 1, NCV, scale,    &
     &                   Bvec(Nstr(ng):,1))
# else
      IF (OutThread) THEN
        start(1)=Nstr(ng)
        total(1)=Nend(ng)-Nstr(ng)+1
        start(2)=1
        total(2)=NCV
        status=nf_inq_varid(ncGSTid(ng), 'Bvec', varid)
        status=nf_put_vara_TYPE(ncGSTid(ng), varid, start, total,       &
     &                          Bvec(Nstr(ng):,1))
      END IF
# endif
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) 'Bvec', TRIM(GSTname(ng))
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out eigenproblem residual vector.
!
# ifdef DISTRIBUTE
      var='resid'
      status=mp_ncwrite(ng, model, ncGSTid(ng), var, GSTname(ng),       &
     &                  vrecord, Nstr(ng), Nend(ng), 1, 1, scale,       &
     &                  resid(Nstr(ng):))
# else
      IF (OutThread) THEN
        start(1)=Nstr(ng)
        total(1)=Nend(ng)-Nstr(ng)+1
        status=nf_inq_varid(ncGSTid(ng), 'resid', varid)
        status=nf_put_vara_TYPE(ncGSTid(ng), varid, start, total,       &
     &                          resid(Nstr(ng):))
      END IF
# endif
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) 'resid', TRIM(GSTname(ng))
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out state reverse communication work array.
# ifdef DISTRIBUTE
!  Notice zero arguments indicating node dimension in NetCDF file.
# endif
!
# ifdef DISTRIBUTE
      var='SworkD'
      status=mp_ncwrite(ng, model, ncGSTid(ng), var, GSTname(ng),       &
     &                  vrecord, 1, 3*Nstate(ng), 0, 0, scale,          &
     &                  SworkD)
# else
      IF (OutThread) THEN
        start(1)=1
        total(1)=3*Nstate(ng)
        status=nf_inq_varid(ncGSTid(ng), 'SworkD', varid)
        status=nf_put_vara_TYPE(ncGSTid(ng), varid, start, total,       &
     &                          SworkD)
      END IF
# endif
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) 'SworkD', TRIM(GSTname(ng))
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!  Write out eigenproblem work array.
# ifdef DISTRIBUTE
!  Notice zero arguments indicating node dimension in NetCDF file.
# endif
!
# ifdef DISTRIBUTE
      var='SworkL'
      status=mp_ncwrite(ng, model, ncGSTid(ng), var, GSTname(ng),       &
     &                  vrecord, 1, LworkL, 0, 0, scale,                &
     &                  SworkL)
# else
      IF (OutThread) THEN
        start(1)=1
        total(1)=LworkL
        status=nf_inq_varid(ncGSTid(ng), 'SworkL', varid)
        status=nf_put_vara_TYPE(ncGSTid(ng), varid, start, total,       &
     &                          SworkL)
      END IF
# endif 
      IF (OutThread.and.(status.ne.nf_noerr)) THEN
        WRITE (stdout,10) 'SworkL', TRIM(GSTname(ng))
        exit_flag=3
        ioerror=status
        RETURN
      END IF
!
!-----------------------------------------------------------------------
!  Synchronize GST checkpointing NetCDF file to disk so the file
!  is available to other processes.
!-----------------------------------------------------------------------
!
      IF (OutThread) THEN
        status=nf_sync(ncGSTid(ng))
        IF (status.ne.nf_noerr) THEN
          WRITE (stdout,20)
          exit_flag=3
          ioerror=status
          RETURN
        END IF
      END IF
      IF (Master) THEN
        WRITE (stdout,30) Nrun+1
      END IF
!
  10  FORMAT (/,' WRT_GST - error while writing variable: ', a,         &
     &        /,11x,'into GST restart NetCDF file: ', a)
  20  FORMAT (/,' WRT_GST - unable to synchronize restart NetCDF to ',  &
     &        'disk.')
  30  FORMAT (6x,'WRT_GST   - wrote GST checkpointing fields at ',      &
     &           'iteration: ', i5.5)

      RETURN
      END SUBROUTINE wrt_gst
#else
      SUBROUTINE wrt_gst
      RETURN
      END SUBROUTINE wrt_gst
#endif
