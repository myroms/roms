#include "cppdefs.h"

#if defined OPT_PERTURBATION
# define ENERGYNORM_SCALE
#endif

#ifdef FULL_GRID
# define IR_RANGE IstrT,IendT
# define IU_RANGE IstrP,IendT
# define JR_RANGE JstrT,JendT
# define JV_RANGE JstrP,JendT
#else
# define IR_RANGE Istr,Iend
# define IU_RANGE IstrU,Iend
# define JR_RANGE Jstr,Jend
# define JV_RANGE JstrV,Jend
#endif

      MODULE packing_mod

#ifdef PROPAGATOR
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2012 The ROMS/TOMS Group       Andrew M. Moore   !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  These routines pack and unpack model state varaibles into/from a    !
!  single vector to interface with  ARPACKs Arnoldi Method  for the    !
!  computation Ritz eigenfunctions.                                    !
!                                                                      !
!=======================================================================
!
      implicit none
!
      PRIVATE
      PUBLIC :: c_norm2
      PUBLIC :: r_norm2
# ifdef ADJOINT
      PUBLIC :: ad_pack
      PUBLIC :: ad_unpack
# endif
# ifdef TANGENT
      PUBLIC :: tl_pack
      PUBLIC :: tl_unpack
# endif
# ifdef SO_SEMI
#  ifdef SO_SEMI_WHITE
      PUBLIC :: so_semi_white
#  else
      PUBLIC :: so_semi_red
#  endif
# endif
!
      CONTAINS
!
      SUBROUTINE c_norm2 (ng, model, Mstr, Mend,                        &
     &                    EvalueR, EvalueI, EvectorR, EvectorI,         &
     &                    state, norm2)
!
!=======================================================================
!                                                                      !
!  This function computes the Euclidean norm between the  propagator   !
!  real/imaginary Ritz eigenvalue (EvalueR, EvalueI) and eigenvector   !
!  (EvectorR, EvectorI) with state vector (state):                     !
!                                                                      !
!     norm2 = Euclidean NORM (state(:) + EvalueR * EvectorR(:) +       !
!                                        EvalueI * EvectorI(:))        !
!                                                                      !
!  WARNING: This function is only intended for serial or distributed   !
!           memory applications.  There is not tiled partitions. All   !
!           quantities are vectors. It replaces the calls to "daxpy"   !
!           and "dnrm2" from the BLAS library. This "legacy" library   !
!           gives  different  results when called inside modules and   !
!           the input arguments are pointers (specially using ifort).  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: Mstr, Mend

      real(r8), intent(in) :: EvalueR
      real(r8), intent(in) :: EvalueI

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: EvectorR(Mstr:)
      real(r8), intent(in) :: EvectorI(Mstr:)
      real(r8), intent(in) :: state(Mstr:)
# else
      real(r8), intent(in) :: EvectorR(Mstr:Mend)
      real(r8), intent(in) :: EvectorI(Mstr:Mend)
      real(r8), intent(in) :: state(Mstr:Mend)
# endif
      real(r8), intent(out) :: norm2
!
!  Local variable declarations.
!
      integer :: NSUB, is

      real(r8) :: cff, my_norm2

# ifdef DISTRIBUTE
      character (len=3), dimension(1) :: op_handle
# endif
!
!-----------------------------------------------------------------------
!  Compute the Euclidean norm of: state(:) + Rvalue * Rvector(:)
!-----------------------------------------------------------------------
!
!  Accumulate squared sum.
!
      my_norm2=0.0_r8
      DO is=Mstr,Mend
        cff=state(is)+EvalueR*EvectorR(is)+                             &
     &                EvalueI*EvectorI(is)
        my_norm2=my_norm2+cff*cff
      END DO
!
!  Take sum squared-root: perform global reduction.
!
      NSUB=1                         ! no horizontal tiles
      IF (tile_count.eq.0) THEN
        norm2=0.0_r8
      END IF
      norm2=norm2+my_norm2
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        op_handle(1)='SUM'
        CALL mp_reduce (ng, model, 1, norm2, op_handle(1))
# endif
      END IF
      norm2=SQRT(norm2)

      RETURN
      END SUBROUTINE c_norm2
!
      SUBROUTINE r_norm2 (ng, model, Mstr, Mend,                        &
     &                    Evalue, Evector, state, norm2)
!
!=======================================================================
!                                                                      !
!  This function computes the Euclidean norm between the propagator    !
!  real Ritz eigenvalue (Evalue) and eigenvector (Evector) with the    !
!  state vector (state):                                               !
!                                                                      !
!     norm2 = Euclidean NORM (state(:) + Evalue * Evector(:))          !
!                                                                      !
!  WARNING: This function is only intended for serial or distributed   !
!           memory applications.  There is not tiled partitions. All   !
!           quantities are vectors. It replaces the calls to "daxpy"   !
!           and "dnrm2" from the BLAS library. This "legacy" library   !
!           gives  different  results when called inside modules and   !
!           the input arguments are pointers (specially using ifort).  !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel

# ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, model
      integer, intent(in) :: Mstr, Mend

      real(r8), intent(in) :: Evalue

# ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Evector(Mstr:)
      real(r8), intent(in) :: state(Mstr:)
# else
      real(r8), intent(in) :: Evector(Mstr:Mend)
      real(r8), intent(in) :: state(Mstr:Mend)
# endif
      real(r8), intent(out) :: norm2
!
!  Local variable declarations.
!
      integer :: NSUB, is

      real(r8) :: cff, my_norm2

# ifdef DISTRIBUTE
      character (len=3), dimension(1) :: op_handle
# endif
!
!-----------------------------------------------------------------------
!  Compute the Euclidean norm of: state(:) + Rvalue * Rvector(:)
!-----------------------------------------------------------------------
!
!  Accumulate squared sum.
!
      my_norm2=0.0_r8
      DO is=Mstr,Mend
        cff=state(is)+Evalue*Evector(is)
        my_norm2=my_norm2+cff*cff
      END DO
!
!  Take sum squared-root: perform global reduction.
!
      NSUB=1                         ! no horizontal tiles
      IF (tile_count.eq.0) THEN
        norm2=0.0_r8
      END IF
      norm2=norm2+my_norm2
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
# ifdef DISTRIBUTE
        op_handle(1)='SUM'
        CALL mp_reduce (ng, model, 1, norm2, op_handle(1))
# endif
      END IF
      norm2=SQRT(norm2)

      RETURN
      END SUBROUTINE r_norm2

# if defined ADJOINT && (defined FORCING_SV || defined SO_SEMI)
!
      SUBROUTINE ad_pack (ng, tile, Mstr, Mend, ad_state)
!
!=======================================================================
!                                                                      !
!  This routine packs the adjoint variables into the state vector.     !
!  The state vector contains only interior water points.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_scatter_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
#  endif

      CALL ad_pack_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   kstp(ng),                                      &
#  ifdef SOLVE3D
     &                   nstp(ng),                                      &
#  endif
#  ifdef MASKING
     &                   GRID(ng) % IJwaterR,                           &
     &                   GRID(ng) % IJwaterU,                           &
     &                   GRID(ng) % IJwaterV,                           &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#  endif
#  ifdef DISTRIBUTE
     &                   1, Mstate(ng), Swork)
#  else
     &                   Mstr, Mend, ad_state)
#  endif

#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Scatter (global to threaded) adjoint state solution to all
!  distributed nodes.
!
      CALL mp_scatter_state (ng, iADM, Mstr, Mend, Mstate(ng),          &
     &                       Swork, ad_state)
#  endif

      RETURN
      END SUBROUTINE ad_pack
!
!***********************************************************************
      SUBROUTINE ad_pack_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         kstp,                                    &
#  ifdef SOLVE3D
     &                         nstp,                                    &
#  endif
#  ifdef MASKING
     &                         IJwaterR, IJwaterU, IJwaterV,            &
     &                         rmask, umask, vmask,                     &
#  endif
     &                         Mstr, Mend, ad_state)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_ncparam
      USE mod_scalars
      USE mod_ocean
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, icount, is, itrc, j, k

      integer, dimension(7+2*NT(ng)) :: offset

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: cff, scale

#  include "set_bounds.h"

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Initialize adjoint state vector with special value (zero) to
!  facilitate gathering/scattering communications between all nodes.
!  This is achieved by summing all the buffers.
!-----------------------------------------------------------------------
!
      DO is=Mstr,Mend
        ad_state(is)=Aspv
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load adjoint STATE and FORCING variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
      icount=0
      offset=0
#  ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=icount
        icount=icount+NwaterR(ng)
      END IF
#   ifdef SOLVE3D
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=icount
        icount=icount+N(ng)*NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=icount
        icount=icount+N(ng)*NwaterV(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=icount
          icount=icount+N(ng)*NwaterR(ng)
        END IF
      END DO
#   else
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=icount
        icount=icount+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=icount
        icount=icount+NwaterV(ng)
      END IF
#   endif
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=icount
        icount=icount+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=icount
        icount=icount+NwaterV(ng)
      END IF
#   ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          offset(isTsur(itrc))=icount
          icount=icount+NwaterR(ng)
        END IF
      END DO
#   endif

#  else

      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=icount
#   ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+2)
#   else
        icount=icount+Lm(ng)*Mm(ng)
#   endif
      END IF
#   ifdef SOLVE3D
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
#    else
        icount=icount+(Lm(ng)-1)*Mm(ng)*N(ng)
#    endif
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
#    else
        icount=icount+Lm(ng)*(Mm(ng)-1)*N(ng)
#    endif
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=icount
#    ifdef FULL_GRID
          icount=icount+(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
#    else
          icount=icount+Lm(ng)*Mm(ng)*N(ng)
#    endif
        END IF
      END DO
#   else
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+1)*(Mm(ng)+2)
#    else
        icount=icount+(Lm(ng)-1)*Mm(ng)
#    endif
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+1)
#    else
        icount=icount+Lm(ng)*(Mm(ng)-1)
#    endif
      END IF
#   endif
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=icount
#   ifdef FULL_GRID
        icount=icount+(Lm(ng)+1)*(Mm(ng)+2)
#   else
        icount=icount+(Lm(ng)-1)*Mm(ng)
#   endif
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=icount
#   ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+1)
#   else
        icount=icount+Lm(ng)*(Mm(ng)-1)
#   endif
      END IF
#   ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          offset(isTsur(itrc))=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+2)
#    else
        icount=icount+Lm(ng)*Mm(ng)
#    endif
        END IF
      END DO
#   endif
#  endif
!
!  Load adjoint of free-surface.
!
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=0
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1
#   endif
#  endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (rmask(i,j).gt.0.0_r8) THEN
              is=IJwaterR(i,j)+offset(isFsur)
              ad_state(is)=scale*OCEAN(ng)%ad_zeta(i,j,kstp)
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
            ad_state(is)=scale*OCEAN(ng)%ad_zeta(i,j,kstp)
#  endif
          END DO
        END DO
      END IF

#  ifndef SOLVE3D
!
!  Load adjoint of 2D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-1
        Ioff=1
        Joff=1
#    endif
#   endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUbar)
              ad_state(is)=scale*OCEAN(ng)%ad_ubar(i,j,kstp)
            END IF
#   else
              is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
              ad_state(is)=scale*OCEAN(ng)%ad_ubar(i,j,kstp)
#   endif
          END DO
        END DO
      END IF
!
!  Load adjoint of 2D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Ioff=0
        Joff=2
#    endif
#   endif
        scale=1.0_r8
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVbar)
              ad_state(is)=scale*OCEAN(ng)%ad_vbar(i,j,kstp)
            END IF
#   else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
            ad_state(is)=scale*OCEAN(ng)%ad_vbar(i,j,kstp)
#   endif
          END DO
        END DO
      END IF

#  else
!
!  Load adjoint of 3D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Jmax=Mm(ng)+2
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-1
        Jmax=Mm(ng)
        Ioff=1
        Joff=1
#    endif
#   endif
        scale=1.0_r8
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
          DO j=JR_RANGE
            DO i=IU_RANGE
#   ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
                is=IJwaterU(i,j)+iadd
                ad_state(is)=scale*OCEAN(ng)%ad_u(i,j,k,nstp)
              END IF
#   else
              is=(i-Ioff)+(j-Joff)*Imax+iadd
              ad_state(is)=scale*OCEAN(ng)%ad_u(i,j,k,nstp)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Load adjoint of 3D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Jmax=Mm(ng)+1
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Jmax=Mm(ng)-1
        Ioff=0
        Joff=2
#    endif
#   endif
        scale=1.0_r8
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
          DO j=JV_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
                is=IJwaterV(i,j)+iadd
                ad_state(is)=scale*OCEAN(ng)%ad_v(i,j,k,nstp)
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              ad_state(is)=scale*OCEAN(ng)%ad_v(i,j,k,nstp)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Load adjoint of tracers variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=1.0_r8
          DO k=1,N(ng)
#   ifdef MASKING
            iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
            iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
            DO j=JR_RANGE
              DO i=IR_RANGE
#   ifdef MASKING
                IF (rmask(i,j).gt.0.0_r8) THEN
                  is=IJwaterR(i,j)+iadd
                  ad_state(is)=scale*OCEAN(ng)%ad_t(i,j,k,nstp,itrc)
                END IF
#   else
                is=(i+Ioff)+(j-Joff)*Imax+iadd
                ad_state(is)=scale*OCEAN(ng)%ad_t(i,j,k,nstp,itrc)
#   endif
              END DO
            END DO
          END DO
        END IF
     END DO
#  endif
!
!  Load adjoint of surface U-stress.
!
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#   else
        Imax=Lm(ng)-1
        Ioff=1
        Joff=1
#   endif
#  endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUstr)
              ad_state(is)=scale*FORCES(ng)%ad_sustr(i,j)
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
            ad_state(is)=scale*FORCES(ng)%ad_sustr(i,j)
#  endif
          END DO
        END DO
      END IF
!
!  Load adjoint of surface V-stress.
!
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=2
#   endif
#  endif
        scale=1.0_r8
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVstr)
              ad_state(is)=scale*FORCES(ng)%ad_svstr(i,j)
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            ad_state(is)=scale*FORCES(ng)%ad_svstr(i,j)
#  endif
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
!
!  Load adjoint of tracers variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=1.0_r8
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+offset(isTsur(itrc))
                ad_state(is)=scale*FORCES(ng)%ad_stflx(i,j,itrc)
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isTsur(itrc))
              ad_state(is)=scale*FORCES(ng)%ad_stflx(i,j,itrc)
#   endif
            END DO
          END DO
        END IF
      END DO
#  endif

      RETURN
      END SUBROUTINE ad_pack_tile

# elif defined ADJOINT
!
      SUBROUTINE ad_pack (ng, tile, Mstr, Mend, ad_state)
!
!=======================================================================
!                                                                      !
!  This routine packs the adjoint variables into the state vector.     !
!  The state vector contains only interior water points.               !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_scatter_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
#  endif

      CALL ad_pack_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   kstp(ng),                                      &
#  ifdef SOLVE3D
     &                   nstp(ng),                                      &
#  endif
#  ifdef DISTRIBUTE
     &                   1, Mstate(ng), Swork,                          &
#  else
     &                   Mstr, Mend, ad_state,                          &
#  endif
#  ifdef MASKING
     &                   GRID(ng) % IJwaterR,                           &
     &                   GRID(ng) % IJwaterU,                           &
     &                   GRID(ng) % IJwaterV,                           &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   GRID(ng) % h,                                  &
#  ifdef SOLVE3D
     &                   GRID(ng) % Hz,                                 &
     &                   OCEAN(ng) % ad_t,                              &
     &                   OCEAN(ng) % ad_u,                              &
     &                   OCEAN(ng) % ad_v,                              &
#  else
     &                   OCEAN(ng) % ad_ubar,                           &
     &                   OCEAN(ng) % ad_vbar,                           &
#  endif
     &                   OCEAN(ng) % ad_zeta)

#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Scatter (global to threaded) adjoint state solution to all
!  distributed nodes.
!
      CALL mp_scatter_state (ng, iADM, Mstr, Mend, Mstate(ng),          &
     &                       Swork, ad_state)
#  endif

      RETURN
      END SUBROUTINE ad_pack
!
!***********************************************************************
      SUBROUTINE ad_pack_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         kstp,                                    &
#  ifdef SOLVE3D
     &                         nstp,                                    &
#  endif
     &                         Mstr, Mend, ad_state,                    &
#  ifdef MASKING
     &                         IJwaterR, IJwaterU, IJwaterV,            &
     &                         rmask, umask, vmask,                     &
#  endif
     &                         h,                                       &
#  ifdef SOLVE3D
     &                         Hz,                                      &
     &                         ad_t, ad_u, ad_v,                        &
#  else
     &                         ad_ubar, ad_vbar,                        &
#  endif
     &                         ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
      real(r8), intent(out) :: ad_state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: cff, scale

#  include "set_bounds.h"

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Initialize adjoint state vector with special value (zero) to
!  facilitate gathering/scattering communications between all nodes.
!  This is achieved by summing all the buffers.
!-----------------------------------------------------------------------
!
      DO is=Mstr,Mend
        ad_state(is)=Aspv
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load adjoint STATE variables into full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
#  ifdef SOLVE3D
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#    else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-1)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#    endif
#   endif
#  else
#    ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#    else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-1)*Mm(ng)
#    endif
#   endif
#  endif
!
!  Load adjoint free-surface.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
      scale=1.0_r8
#  endif
      DO j=JR_RANGE
        DO i=IR_RANGE
#  ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            ad_state(is)=scale*ad_zeta(i,j,kstp)
          END IF
#  else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          ad_state(is)=scale*ad_zeta(i,j,kstp)
#  endif
        END DO
      END DO

#  ifndef SOLVE3D
!
!  Load adjoint 2D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JR_RANGE
        DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            ad_state(is)=scale*ad_ubar(i,j,kstp)
          END IF
#   else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          ad_state(is)=scale*ad_ubar(i,j,kstp)
#   endif
        END DO
      END DO
!
!  Load adjoint 2D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JV_RANGE
        DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            ad_state(is)=scale*ad_vbar(i,j,kstp)
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          ad_state(is)=scale*ad_vbar(i,j,kstp)
#   endif
        END DO
      END DO
#  else
!
!  Load adjoint 3D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-1
      Jmax=Mm(ng)
      Ioff=1
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
              is=IJwaterU(i,j)+iadd
              ad_state(is)=scale*ad_u(i,j,k,nstp)
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            ad_state(is)=scale*ad_u(i,j,k,nstp)
#   endif
          END DO
        END DO
      END DO
!
!  Load adjoint 3D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)-1
      Ioff=0
      Joff=2
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
              is=IJwaterV(i,j)+iadd
              ad_state(is)=scale*ad_v(i,j,k,nstp)
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            ad_state(is)=scale*ad_v(i,j,k,nstp)
#   endif
          END DO
        END DO
      END DO
!
!  Load adjoint tracers variables. For now, use salinity scale for
!  passive tracers.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)
      Ioff=0
      Joff=1
#    endif
#   endif
      DO itrc=1,NT(ng)
#   ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                is=IJwaterR(i,j)+iadd
                ad_state(is)=scale*ad_t(i,j,k,nstp,itrc)
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              ad_state(is)=scale*ad_t(i,j,k,nstp,itrc)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif
      RETURN
      END SUBROUTINE ad_pack_tile

# endif

# if defined ADJOINT && defined SO_SEMI
!
      SUBROUTINE ad_unpack (ng, tile, Mstr, Mend, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the requested adjoint state and/or  surface    !
!  forcing variables used in stochastic optimals.  The state vector    !
!  contains only interior water points.                                !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_forces
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_gather_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef DISTRIBUTE
!
!  Gather (threaded to global) adjoint state solution from all
!  distributed nodes.
!
      CALL mp_gather_state (ng, iNLM, Mstr, Mend, Mstate(ng),           &
     &                      state, Swork)
!
#  endif

#  ifdef PROFILE
      CALL wclock_on (ng, iNLM, 2)
#  endif

      CALL ad_unpack_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     kstp(ng),                                    &
#  ifdef SOLVE3D
     &                     nstp(ng),                                    &
#  endif
#  ifdef MASKING
     &                     GRID(ng) % IJwaterR,                         &
     &                     GRID(ng) % IJwaterU,                         &
     &                     GRID(ng) % IJwaterV,                         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#  endif
#  ifdef DISTRIBUTE
     &                     1, Mstate(ng), Swork)
#  else
     &                     Mstr, Mend, state)
#  endif

#  ifdef PROFILE
      CALL wclock_off (ng, iNLM, 2)
#  endif

      RETURN
      END SUBROUTINE ad_unpack
!
!***********************************************************************
      SUBROUTINE ad_unpack_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           kstp,                                  &
#  ifdef SOLVE3D
     &                           nstp,                                  &
#  endif
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
     &                           Mstr, Mend, state)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: kstp
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
      integer, intent(in) :: Mstr, Mend
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, icount, is, itrc, j, k

#  ifdef SOLVE3D
      integer, dimension(7+2*NT(ng)) :: offset
#  else
      integer, dimension(5) :: offset
#  endif

      real(r8) :: scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract adjoint FORCING variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
      icount=0
      offset=0
#  ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=icount
        icount=icount+NwaterR(ng)
      END IF
#   ifdef SOLVE3D
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=icount
        icount=icount+N(ng)*NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=icount
        icount=icount+N(ng)*NwaterV(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=icount
          icount=icount+N(ng)*NwaterR(ng)
        END IF
      END DO
#   else
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=icount
        icount=icount+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=icount
        icount=icount+NwaterV(ng)
      END IF
#   endif
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=icount
        icount=icount+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=icount
        icount=icount+NwaterV(ng)
      END IF
#   ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          offset(isTsur(itrc))=icount
          icount=icount+NwaterR(ng)
        END IF
      END DO
#   endif

#  else

      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=icount
#   ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+2)
#   else
        icount=icount+Lm(ng)*Mm(ng)
#   endif
      END IF
#   ifdef SOLVE3D
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
#    else
        icount=icount+(Lm(ng)-1)*Mm(ng)*N(ng)
#    endif
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
#    else
        icount=icount+Lm(ng)*(Mm(ng)-1)*N(ng)
#    endif
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=icount
#    ifdef FULL_GRID
          icount=icount+(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
#    else
          icount=icount+Lm(ng)*Mm(ng)*N(ng)
#    endif
        END IF
      END DO
#   else
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+1)*(Mm(ng)+2)
#    else
        icount=icount+(Lm(ng)-1)*Mm(ng)
#    endif
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+1)
#    else
        icount=icount+Lm(ng)*(Mm(ng)-1)
#    endif
      END IF
#   endif
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=icount
#   ifdef FULL_GRID
        icount=icount+(Lm(ng)+1)*(Mm(ng)+2)
#   else
        icount=icount+(Lm(ng)-1)*Mm(ng)
#   endif
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=icount
#   ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+1)
#   else
        icount=icount+Lm(ng)*(Mm(ng)-1)
#   endif
      END IF
#   ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          offset(isTsur(itrc))=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+2)
#    else
        icount=icount+Lm(ng)*Mm(ng)
#    endif
        END IF
      END DO
#   endif
#  endif
!
!  Extract adjoint free-surface.
!
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=0
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1
#   endif
#  endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isFsur)
              OCEAN(ng)%ad_zeta(i,j,kstp)=scale*state(is)
            ELSE
              OCEAN(ng)%ad_zeta(i,j,kstp)=0.0_r8
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isFsur)
            OCEAN(ng)%ad_zeta(i,j,kstp)=scale*state(is)
#  endif
          END DO
        END DO
      END IF

#  ifndef SOLVE3D
!
!  Extract adjoint 2D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-1
        Ioff=1
        Joff=1
#    endif
#   endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUbar)
              OCEAN(ng)%ad_ubar(i,j,kstp)=scale*state(is)
            ELSE
              OCEAN(ng)%ubar(i,j,kstp)=0.0_r8
            END IF
#   else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
            OCEAN(ng)%ad_ubar(i,j,kstp)=scale*state(is)
#   endif
          END DO
        END DO
      END IF
!
!  Extract adjoint 2D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Ioff=0
        Joff=2
#    endif
#   endif
        scale=1.0_r8
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVbar)
              OCEAN(ng)%ad_vbar(i,j,kstp)=scale*state(is)
            ELSE
              OCEAN(ng)%ad_vbar(i,j,kstp)=0.0_r8
            END IF
#   else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            OCEAN(ng)%ad_vbar(i,j,kstp)=scale*state(is)
#   endif
          END DO
        END DO
      END DO

#  else
!
!  Extract adjoint 3D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Jmax=Mm(ng)+2
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-1
        Jmax=Mm(ng)
        Ioff=1
        Joff=1
#    endif
#   endif
        scale=1.0_r8
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
          DO j=JR_RANGE
            DO i=IU_RANGE
#   ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
                is=IJwaterU(i,j)+iadd
                OCEAN(ng)%ad_u(i,j,k,nstp)=scale*state(is)
              ELSE
                OCEAN(ng)%ad_u(i,j,k,nstp)=0.0_r8
              END IF
#   else
              is=(i-Ioff)+(j-Joff)*Imax+iadd
              OCEAN(ng)%ad_u(i,j,k,nstp)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Extract adjoint 3D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Jmax=Mm(ng)+1
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Jmax=Mm(ng)-1
        Ioff=0
        Joff=2
#    endif
#   endif
        scale=1.0_r8
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
          DO j=JV_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
                is=IJwaterV(i,j)+iadd
                OCEAN(ng)%ad_v(i,j,k,nstp)=scale*state(is)
              ELSE
                OCEAN(ng)%ad_v(i,j,k,nstp)=0.0_r8
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              OCEAN(ng)%ad_v(i,j,k,nstp)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Extract adjoint tracers variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=1.0_r8
          DO k=1,N(ng)
#   ifdef MASKING
            iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
            iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
            DO j=JR_RANGE
              DO i=IR_RANGE
#   ifdef MASKING
                IF (rmask(i,j).gt.0.0_r8) THEN
                  is=IJwaterR(i,j)+iadd
                  OCEAN(ng)%ad_t(i,j,k,nstp,itrc)=scale*state(is)
                ELSE
                  OCEAN(ng)%ad_t(i,j,k,nstp,itrc)=0.0_r8
                END IF
#   else
                is=(i+Ioff)+(j-Joff)*Imax+iadd
                OCEAN(ng)%ad_t(i,j,k,nstp,itrc)=scale*state(is)
#   endif
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
!
!  Extract adjoint surface U-momentum stress.
!
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#   else
        Imax=Lm(ng)-1
        Ioff=1
        Joff=1
#   endif
#  endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUstr)
              FORCES(ng)%ad_sustr(i,j)=scale*state(is)
            ELSE
              FORCES(ng)%ad_sustr(i,j)=0.0_r8
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
            FORCES(ng)%ad_sustr(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF
!
!  Extract adjoint surface V-momentum stress.
!
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=2
#   endif
#  endif
        scale=1.0_r8
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVstr)
              FORCES(ng)%ad_svstr(i,j)=scale*state(is)
            ELSE
              FORCES(ng)%ad_svstr(i,j)=0.0_r8
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            FORCES(ng)%ad_svstr(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
!
!  Extract adjoint surface tracer flux variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=1.0_r8
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+offset(isTvar(itrc))
                FORCES(ng)%ad_stflx(i,j,itrc)=scale*state(is)
              ELSE
                FORCES(ng)%ad_stflx(i,j,itrc)=0.0_r8
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isTvar(itrc))
              FORCES(ng)%ad_stflx(i,j,itrc)=scale*state(is)
#   endif
            END DO
          END DO
        END IF
      END DO
#  endif

      RETURN
      END SUBROUTINE ad_unpack_tile

# elif defined ADJOINT
!
      SUBROUTINE ad_unpack (ng, tile, Mstr, Mend, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the adjoint model variables from the state     !
!  vector. If applicable,  the state vector includes only unmasked     !
!  water points.                                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_gather_state
#  endif

!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iADM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Gather (threaded to global) adjoint state solution from all
!  distributed nodes.
!
      CALL mp_gather_state (ng, iTLM, Mstr, Mend, Mstate(ng),           &
     &                      state, Swork)
!
#  endif
      CALL ad_unpack_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     knew(ng),                                    &
#  ifdef SOLVE3D
     &                     nstp(ng),                                    &
#  endif
#  ifdef DISTRIBUTE
     &                     1, Mstate(ng), Swork,                        &
#  else
     &                     Mstr, Mend, state,                           &
#  endif
#  ifdef MASKING
     &                     GRID(ng) % IJwaterR,                         &
     &                     GRID(ng) % IJwaterU,                         &
     &                     GRID(ng) % IJwaterV,                         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     GRID(ng) % h,                                &
#  ifdef SOLVE3D
     &                     GRID(ng) % Hz,                               &
     &                     OCEAN(ng) % ad_t,                            &
     &                     OCEAN(ng) % ad_u,                            &
     &                     OCEAN(ng) % ad_v,                            &
#  else
     &                     OCEAN(ng) % ad_ubar,                         &
     &                     OCEAN(ng) % ad_vbar,                         &
#  endif
     &                     OCEAN(ng) % ad_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iADM, 2)
#  endif
      RETURN
      END SUBROUTINE ad_unpack
!
!***********************************************************************
      SUBROUTINE ad_unpack_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           knew,                                  &
#  ifdef SOLVE3D
     &                           nstp,                                  &
#  endif
     &                           Mstr, Mend, state,                     &
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
     &                           h,                                     &
#  ifdef SOLVE3D
     &                           Hz,                                    &
     &                           ad_t, ad_u, ad_v,                      &
#  else
     &                           ad_ubar, ad_vbar,                      &
#  endif
     &                           ad_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: knew
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: state(Mstr:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: ad_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: ad_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: ad_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: ad_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: state(Mstr:Mend)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: ad_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: ad_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: ad_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(inout) :: ad_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: ad_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: ad_zeta(LBi:UBi,LBj:UBj,3)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8) :: cff, scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract adjoint state variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
#  ifdef SOLVE3D
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#    else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-1)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#    else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-1)*Mm(ng)
#    endif
#   endif
#  endif
!
!  Extract adjoint free-surface.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
      scale=1.0_r8
#  endif
      DO j=JR_RANGE
        DO i=IR_RANGE
#  ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            ad_zeta(i,j,knew)=scale*state(is)
          ELSE
            ad_zeta(i,j,knew)=0.0_r8
          END IF
#  else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          ad_zeta(i,j,knew)=scale*state(is)
#  endif
        END DO
      END DO
#  ifndef SOLVE3D
!
!  Extract adjoint 2D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JR_RANGE
        DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            ad_ubar(i,j,knew)=scale*state(is)
          ELSE
            ad_ubar(i,j,knew)=0.0_r8
          END IF
#   else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          ad_ubar(i,j,knew)=scale*state(is)
#   endif
        END DO
      END DO
!
!  Extract adjoint 2D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JV_RANGE
        DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            ad_vbar(i,j,knew)=scale*state(is)
          ELSE
            ad_vbar(i,j,knew)=0.0_r8
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          ad_vbar(i,j,knew)=scale*state(is)
#   endif
        END DO
      END DO
#  else
!
!  Extract adjoint 3D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-1
      Jmax=Mm(ng)
      Ioff=1
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
              is=IJwaterU(i,j)+iadd
              ad_u(i,j,k,nstp)=scale*state(is)
            ELSE
              ad_u(i,j,k,nstp)=0.0_r8
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            ad_u(i,j,k,nstp)=scale*state(is)
#   endif
          END DO
        END DO
      END DO
!
!  Extract adjoint 3D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)-1
      Ioff=0
      Joff=2
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
              is=IJwaterV(i,j)+iadd
              ad_v(i,j,k,nstp)=scale*state(is)
            ELSE
              ad_v(i,j,k,nstp)=0.0_r8
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            ad_v(i,j,k,nstp)=scale*state(is)
#   endif
          END DO
        END DO
      END DO
!
!  Extract tangent linear tracers variables. For now, use salinity scale
!  for passive tracers.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)
      Ioff=0
      Joff=1
#    endif
#   endif
      DO itrc=1,NT(ng)
#   ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                is=IJwaterR(i,j)+iadd
                ad_t(i,j,k,nstp,itrc)=scale*state(is)
              ELSE
                ad_t(i,j,k,nstp,itrc)=0.0_r8
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              ad_t(i,j,k,nstp,itrc)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif

      RETURN
      END SUBROUTINE ad_unpack_tile
# endif

# ifdef TANGENT
!
      SUBROUTINE tl_pack (ng, tile, Mstr, Mend, tl_state)
!
!=======================================================================
!                                                                      !
!  This routine packs the tangent linear variables into the state      !
!  vetor. The state vector contains only interior water points.        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_scatter_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: tl_state(Mstr:)
#  else
      real(r8), intent(out) :: tl_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
#  endif
      CALL tl_pack_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj,                            &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   krhs(ng), kstp(ng), knew(ng),                  &
#  ifdef SOLVE3D
     &                   nstp(ng),                                      &
#  endif
#  ifdef DISTRIBUTE
     &                   1, Mstate(ng), Swork,                          &
#  else
     &                   Mstr, Mend, tl_state,                          &
#  endif
#  ifdef MASKING
     &                   GRID(ng) % IJwaterR,                           &
     &                   GRID(ng) % IJwaterU,                           &
     &                   GRID(ng) % IJwaterV,                           &
     &                   GRID(ng) % rmask,                              &
     &                   GRID(ng) % umask,                              &
     &                   GRID(ng) % vmask,                              &
#  endif
     &                   GRID(ng) % h,                                  &
#  ifdef SOLVE3D
     &                   GRID(ng) % Hz,                                 &
     &                   OCEAN(ng) % tl_t,                              &
     &                   OCEAN(ng) % tl_u,                              &
     &                   OCEAN(ng) % tl_v,                              &
#  else
     &                   OCEAN(ng) % tl_ubar,                           &
     &                   OCEAN(ng) % tl_vbar,                           &
#  endif
     &                   OCEAN(ng) % tl_zeta)

#  ifdef DISTRIBUTE
!
!  Scatter (global to threaded) tangent linear state solution to all
!  distributed nodes.
!
      CALL mp_scatter_state (ng, iTLM, Mstr, Mend, Mstate(ng),          &
     &                       Swork, tl_state)
#  endif

#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
#  endif
      RETURN
      END SUBROUTINE tl_pack
!
!***********************************************************************
      SUBROUTINE tl_pack_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj,                      &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         krhs, kstp, knew,                        &
#  ifdef SOLVE3D
     &                         nstp,                                    &
#  endif
     &                         Mstr, Mend, tl_state,                    &
#  ifdef MASKING
     &                         IJwaterR, IJwaterU, IJwaterV,            &
     &                         rmask, umask, vmask,                     &
#  endif
     &                         h,                                       &
#  ifdef SOLVE3D
     &                         Hz,                                      &
     &                         tl_t, tl_u, tl_v,                        &
#  else
     &                         tl_ubar, tl_vbar,                        &
#  endif
     &                         tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: krhs, kstp, knew
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)

      real(r8), intent(out) :: tl_state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)

      real(r8), intent(out) :: tl_state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8), parameter :: Aspv = 0.0_r8

      real(r8) :: cff, scale

#  include "set_bounds.h"

#  ifdef DISTRIBUTE
!
!-----------------------------------------------------------------------
!  Initialize tangent linear state vector with special value (zero) to
!  facilitate gathering/scattering communications between all nodes.
!  This is achieved by summing all the buffers.
!-----------------------------------------------------------------------
!
      DO is=Mstr,Mend
        tl_state(is)=Aspv
      END DO
#  endif
!
!-----------------------------------------------------------------------
!  Load tangent linear state variables into full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
#  ifdef SOLVE3D
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#    else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-1)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#    else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-1)*Mm(ng)
#    endif
#   endif
#  endif
!
!  Load tangent linear free-surface.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
      scale=1.0_r8
#  endif
      DO j=JR_RANGE
        DO i=IR_RANGE
#  ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            tl_state(is)=scale*tl_zeta(i,j,knew)
          END IF
#  else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          tl_state(is)=scale*tl_zeta(i,j,knew)
#  endif
        END DO
      END DO
#  ifndef SOLVE3D
!
!  Load tangent linear 2D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JR_RANGE
        DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            tl_state(is)=scale*tl_ubar(i,j,knew)
          END IF
#   else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          tl_state(is)=scale*tl_ubar(i,j,knew)
#   endif
        END DO
      END DO
!
!  Load tangent linear 2D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JV_RANGE
        DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            tl_state(is)=scale*tl_vbar(i,j,knew)
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          tl_state(is)=scale*tl_vbar(i,j,knew)
#   endif
        END DO
      END DO
#  else
!
!  Load tangent linear 3D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-1
      Jmax=Mm(ng)
      Ioff=1
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
              is=IJwaterU(i,j)+iadd
              tl_state(is)=scale*tl_u(i,j,k,nstp)
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            tl_state(is)=scale*tl_u(i,j,k,nstp)
#   endif
          END DO
        END DO
      END DO
!
!  Load tangent linear 3D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)-1
      Ioff=0
      Joff=2
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
              is=IJwaterV(i,j)+iadd
              tl_state(is)=scale*tl_v(i,j,k,nstp)
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            tl_state(is)=scale*tl_v(i,j,k,nstp)
#   endif
          END DO
        END DO
      END DO
!
!  Load tangent linear tracers variables. For now, use salinity scale
!  for passive tracers.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)
      Ioff=0
      Joff=1
#    endif
#   endif
      DO itrc=1,NT(ng)
#   ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                is=IJwaterR(i,j)+iadd
                tl_state(is)=scale*tl_t(i,j,k,nstp,itrc)
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              tl_state(is)=scale*tl_t(i,j,k,nstp,itrc)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif

      RETURN
      END SUBROUTINE tl_pack_tile
# endif

# if defined TANGENT && defined FORCING_SV
!
      SUBROUTINE tl_unpack (ng, tile, Mstr, Mend, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the tangent linear variables from the state    !
!  vector.  If applicable,  the state vector includes only unmasked    !
!  water points.                                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_gather_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Gather (threaded to global) tangent linear state solution from all
!  distributed nodes.
!
      CALL mp_gather_state (ng, iTLM, Mstr, Mend, Mstate(ng),           &
     &                      state, Swork)
!
#  endif
      CALL tl_unpack_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     kstp(ng),                                    &
#  ifdef SOLVE3D
     &                     nstp(ng),                                    &
#  endif
#  ifdef MASKING
     &                     GRID(ng) % IJwaterR,                         &
     &                     GRID(ng) % IJwaterU,                         &
     &                     GRID(ng) % IJwaterV,                         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#  endif
#  ifdef DISTRIBUTE
     &                     1, Mstate(ng), Swork)
#  else
     &                     Mstr, Mend, state)
#  endif

#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
#  endif

      RETURN
      END SUBROUTINE tl_unpack
!
!***********************************************************************
      SUBROUTINE tl_unpack_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           kstp,                                  &
#  ifdef SOLVE3D
     &                           nstp,                                  &
#  endif
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
     &                           Mstr, Mend, state)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_forces
      USE mod_ncparam
      USE mod_ocean
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: state(Mstr:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, icount, is, itrc, j, k

      integer, dimension(7+2*NT(ng)) :: offset

      real(r8) :: cff, scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract tangent linear FORCING variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
      icount=0
      offset=0
#  ifdef MASKING
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=icount
        icount=icount+NwaterR(ng)
      END IF
#   ifdef SOLVE3D
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=icount
        icount=icount+N(ng)*NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=icount
        icount=icount+N(ng)*NwaterV(ng)
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=icount
          icount=icount+N(ng)*NwaterR(ng)
        END IF
      END DO
#   else
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=icount
        icount=icount+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=icount
        icount=icount+NwaterV(ng)
      END IF
#   endif
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=icount
        icount=icount+NwaterU(ng)
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=icount
        icount=icount+NwaterV(ng)
      END IF
#   ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          offset(isTsur(itrc))=icount
          icount=icount+NwaterR(ng)
        END IF
      END DO
#   endif

#  else

      IF (SCALARS(ng)%Fstate(isFsur)) THEN
        offset(isFsur)=icount
#   ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+2)
#   else
        icount=icount+Lm(ng)*Mm(ng)
#   endif
      END IF
#   ifdef SOLVE3D
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
        offset(isUvel)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
#    else
        icount=icount+(Lm(ng)-1)*Mm(ng)*N(ng)
#    endif
      END IF
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
        offset(isVvel)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
#    else
        icount=icount+Lm(ng)*(Mm(ng)-1)*N(ng)
#    endif
      END IF
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
          offset(isTvar(itrc))=icount
#    ifdef FULL_GRID
          icount=icount+(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
#    else
          icount=icount+Lm(ng)*Mm(ng)*N(ng)
#    endif
        END IF
      END DO
#   else
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
        offset(isUbar)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+1)*(Mm(ng)+2)
#    else
        icount=icount+(Lm(ng)-1)*Mm(ng)
#    endif
      END IF
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
        offset(isVbar)=icount
#    ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+1)
#    else
        icount=icount+Lm(ng)*(Mm(ng)-1)
#    endif
      END IF
#   endif
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
        offset(isUstr)=icount
#   ifdef FULL_GRID
        icount=icount+(Lm(ng)+1)*(Mm(ng)+2)
#   else
        icount=icount+(Lm(ng)-1)*Mm(ng)
#   endif
      END IF
      IF (SCALARS(ng)%Fstate(isVstr)) THEN
        offset(isVstr)=icount
#   ifdef FULL_GRID
        icount=icount+(Lm(ng)+2)*(Mm(ng)+1)
#   else
        icount=icount+Lm(ng)*(Mm(ng)-1)
#   endif
      END IF
#   ifdef SOLVE3D
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
          offset(isTsur(itrc))=icount
#    ifdef FULL_GRID
          icount=icount+(Lm(ng)+2)*(Mm(ng)+2)
#    else
          icount=icount+Lm(ng)*Mm(ng)
#    endif
        END IF
      END DO
#   endif
#  endif
!
!  Extract tangent linear free-surface.
!
      IF (SCALARS(ng)%Fstate(isFsur)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=0
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=1
#   endif
#  endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isFsur)
              OCEAN(ng)%tl_zeta(i,j,kstp)=scale*state(is)
            ELSE
              OCEAN(ng)%tl_zeta(i,j,kstp)=0.0_r8
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isFsur)
            OCEAN(ng)%tl_zeta(i,j,kstp)=scale*state(is)
#  endif
          END DO
        END DO
      END IF

#  ifndef SOLVE3D
!
!  Extract tangent linear 2D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-1
        Ioff=1
        Joff=1
#    endif
#   endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUbar)
              OCEAN(ng)%tl_ubar(i,j,kstp)=scale*state(is)
            ELSE
              OCEAN(ng)%tl_ubar(i,j,kstp)=0.0_r8
            END IF
#   else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
            OCEAN(ng)%tl_ubar(i,j,kstp)=scale*state(is)
#   endif
          END DO
        END DO
      END IF
!
!  Extract tangent linear 2D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVbar)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Ioff=0
        Joff=2
#    endif
#   endif
        scale=1.0_r8
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVbar)
              OCEAN(ng)%tl_vbar(i,j,kstp)=scale*state(is)
            ELSE
              OCEAN(ng)%tl_vbar(i,j,kstp)=0.0_r8
            END IF
#   else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            OCEAN(ng)%tl_vbar(i,j,kstp)=scale*state(is)
#   endif
          END DO
        END DO
      END DO

#  else
!
!  Extract tangent linear 3D U-velocity.
!
      IF (SCALARS(ng)%Fstate(isUvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+1
        Jmax=Mm(ng)+2
        Ioff=0
        Joff=0
#    else
        Imax=Lm(ng)-1
        Jmax=Mm(ng)
        Ioff=1
        Joff=1
#    endif
#   endif
        scale=1.0_r8
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
          DO j=JR_RANGE
            DO i=IU_RANGE
#   ifdef MASKING
              IF (umask(i,j).gt.0.0_r8) THEN
                is=IJwaterU(i,j)+iadd
                OCEAN(ng)%tl_u(i,j,k,nstp)=scale*state(is)
              ELSE
                OCEAN(ng)%tl_u(i,j,k,nstp)=0.0_r8
              END IF
#   else
              is=(i-Ioff)+(j-Joff)*Imax+iadd
              OCEAN(ng)%tl_u(i,j,k,nstp)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Extract tangent linear 3D V-velocity.
!
      IF (SCALARS(ng)%Fstate(isVvel)) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
        Imax=Lm(ng)+2
        Jmax=Mm(ng)+1
        Ioff=1
        Joff=1
#    else
        Imax=Lm(ng)
        Jmax=Mm(ng)-1
        Ioff=0
        Joff=2
#    endif
#   endif
        scale=1.0_r8
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
          iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
          DO j=JV_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (vmask(i,j).gt.0.0_r8) THEN
                is=IJwaterV(i,j)+iadd
                OCEAN(ng)%tl_v(i,j,k,nstp)=scale*state(is)
              ELSE
                OCEAN(ng)%tl_v(i,j,k,nstp)=0.0_r8
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              OCEAN(ng)%tl_v(i,j,k,nstp)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END IF
!
!  Extract tangent linear tracers variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTvar(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=1.0_r8
          DO k=1,N(ng)
#   ifdef MASKING
            iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
            iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
            DO j=JR_RANGE
              DO i=IR_RANGE
#   ifdef MASKING
                IF (rmask(i,j).gt.0.0_r8) THEN
                  is=IJwaterR(i,j)+iadd
                  OCEAN(ng)%tl_t(i,j,k,nstp,itrc)=scale*state(is)
                ELSE
                  OCEAN(ng)%tl_t(i,j,k,nstp,itrc)=0.0_r8
                END IF
#   else
                is=(i+Ioff)+(j-Joff)*Imax+iadd
                OCEAN(ng)%tl_t(i,j,k,nstp,itrc)=scale*state(is)
#   endif
              END DO
            END DO
          END DO
        END IF
      END DO
#  endif
!
!  Extract tangent linear surface U-momentum stress.
!
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+1
        Ioff=0
        Joff=0
#   else
        Imax=Lm(ng)-1
        Ioff=1
        Joff=1
#   endif
#  endif
        scale=1.0_r8
        DO j=JR_RANGE
          DO i=IU_RANGE
#  ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
              is=IJwaterU(i,j)+offset(isUstr)
              FORCES(ng)%tl_sustr(i,j)=scale*state(is)
            ELSE
              FORCES(ng)%tl_sustr(i,j)=0.0_r8
            END IF
#  else
            is=(i-Ioff)+(j-Joff)*Imax+offset(isUstr)
            FORCES(ng)%tl_sustr(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF
!
!  Extract tangent linear surface V-momentum stress.
!
      IF (SCALARS(ng)%Fstate(isUstr)) THEN
#  ifndef MASKING
#   ifdef FULL_GRID
        Imax=Lm(ng)+2
        Ioff=1
        Joff=1
#   else
        Imax=Lm(ng)
        Ioff=0
        Joff=2
#   endif
#  endif
        scale=1.0_r8
        DO j=JV_RANGE
          DO i=IR_RANGE
#  ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
              is=IJwaterV(i,j)+offset(isVstr)
              FORCES(ng)%tl_svstr(i,j)=scale*state(is)
            ELSE
              FORCES(ng)%tl_svstr(i,j)=0.0_r8
            END IF
#  else
            is=(i+Ioff)+(j-Joff)*Imax+offset(isVstr)
            FORCES(ng)%tl_svstr(i,j)=scale*state(is)
#  endif
          END DO
        END DO
      END IF

#  ifdef SOLVE3D
!
!  Extract tangent linear surface tracer flux variables.
!
      DO itrc=1,NT(ng)
        IF (SCALARS(ng)%Fstate(isTsur(itrc))) THEN
#   ifndef MASKING
#    ifdef FULL_GRID
          Imax=Lm(ng)+2
          Jmax=Mm(ng)+2
          Ioff=1
          Joff=0
#    else
          Imax=Lm(ng)
          Jmax=Mm(ng)
          Ioff=0
          Joff=1
#    endif
#   endif
          scale=1.0_r8
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
                is=IJwaterR(i,j)+offset(isTvar(itrc))
                FORCES(ng)%tl_stflx(i,j,itrc)=scale*state(is)
              ELSE
                FORCES(ng)%tl_stflx(i,j,itrc)=0.0_r8
              END IF
#   else
              is=(i+Ioff)+(j-Joff)*Imax+offset(isTvar(itrc))
              FORCES(ng)%tl_stflx(i,j,itrc)=scale*state(is)
#   endif
            END DO
          END DO
        END IF
      END DO
#  endif

      RETURN
      END SUBROUTINE tl_unpack_tile

# elif defined TANGENT
!
      SUBROUTINE tl_unpack (ng, tile, Mstr, Mend, state)
!
!=======================================================================
!                                                                      !
!  This routine unpacks the tangent linear variables from the state    !
!  vector.  If applicable,  the state vector includes only unmasked    !
!  water points.                                                       !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ocean
      USE mod_stepping
#  ifdef DISTRIBUTE
      USE mod_storage
#  endif
#  ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_gather_state
#  endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#  ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
#  else
      real(r8), intent(in) :: state(Mstr:Mend)
#  endif
!
!  Local variable declarations.
!
#  include "tile.h"
!
#  ifdef PROFILE
      CALL wclock_on (ng, iTLM, 2)
#  endif

#  ifdef DISTRIBUTE
!
!  Gather (threaded to global) tangent linear state solution from all
!  distributed nodes.
!
      CALL mp_gather_state (ng, iTLM, Mstr, Mend, Mstate(ng),           &
     &                      state, Swork)
!
#  endif
      CALL tl_unpack_tile (ng, tile,                                    &
     &                     LBi, UBi, LBj, UBj,                          &
     &                     IminS, ImaxS, JminS, JmaxS,                  &
     &                     kstp(ng),                                    &
#  ifdef SOLVE3D
     &                     nstp(ng),                                    &
#  endif
#  ifdef DISTRIBUTE
     &                     1, Mstate(ng), Swork,                        &
#  else
     &                     Mstr, Mend, state,                           &
#  endif
#  ifdef MASKING
     &                     GRID(ng) % IJwaterR,                         &
     &                     GRID(ng) % IJwaterU,                         &
     &                     GRID(ng) % IJwaterV,                         &
     &                     GRID(ng) % rmask,                            &
     &                     GRID(ng) % umask,                            &
     &                     GRID(ng) % vmask,                            &
#  endif
     &                     GRID(ng) % h,                                &
#  ifdef SOLVE3D
     &                     GRID(ng) % Hz,                               &
     &                     OCEAN(ng) % tl_t,                            &
     &                     OCEAN(ng) % tl_u,                            &
     &                     OCEAN(ng) % tl_v,                            &
#  else
     &                     OCEAN(ng) % tl_ubar,                         &
     &                     OCEAN(ng) % tl_vbar,                         &
#  endif
     &                     OCEAN(ng) % tl_zeta)
#  ifdef PROFILE
      CALL wclock_off (ng, iTLM, 2)
#  endif
      RETURN
      END SUBROUTINE tl_unpack
!
!***********************************************************************
      SUBROUTINE tl_unpack_tile (ng, tile,                              &
     &                           LBi, UBi, LBj, UBj,                    &
     &                           IminS, ImaxS, JminS, JmaxS,            &
     &                           kstp,                                  &
#  ifdef SOLVE3D
     &                           nstp,                                  &
#  endif
     &                           Mstr, Mend, state,                     &
#  ifdef MASKING
     &                           IJwaterR, IJwaterU, IJwaterV,          &
     &                           rmask, umask, vmask,                   &
#  endif
     &                           h,                                     &
#  ifdef SOLVE3D
     &                           Hz,                                    &
     &                           tl_t, tl_u, tl_v,                      &
#  else
     &                           tl_ubar, tl_vbar,                      &
#  endif
     &                           tl_zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_ncparam
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Mstr, Mend
      integer, intent(in) :: kstp
#  ifdef SOLVE3D
      integer, intent(in) :: nstp
#  endif
!
#  ifdef ASSUMED_SHAPE
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:,LBj:)
      integer, intent(in) :: IJwaterU(LBi:,LBj:)
      integer, intent(in) :: IJwaterV(LBi:,LBj:)

      real(r8), intent(in) :: rmask(LBi:,LBj:)
      real(r8), intent(in) :: umask(LBi:,LBj:)
      real(r8), intent(in) :: vmask(LBi:,LBj:)
#   endif
      real(r8), intent(in) :: state(Mstr:)
      real(r8), intent(in) :: h(LBi:,LBj:)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)

      real(r8), intent(inout) :: tl_t(LBi:,LBj:,:,:,:)
      real(r8), intent(inout) :: tl_u(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: tl_v(LBi:,LBj:,:,:)
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:,LBj:,:)
      real(r8), intent(inout) :: tl_vbar(LBi:,LBj:,:)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:,LBj:,:)
#  else
#   ifdef MASKING
      integer, intent(in) :: IJwaterR(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterU(LBi:UBi,LBj:UBj)
      integer, intent(in) :: IJwaterV(LBi:UBi,LBj:UBj)

      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: umask(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: vmask(LBi:UBi,LBj:UBj)
#   endif
      real(r8), intent(in) :: state(Mstr:Mend)
      real(r8), intent(in) :: h(LBi:UBi,LBj:UBj)
#   ifdef SOLVE3D
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,N(ng))

      real(r8), intent(inout) :: tl_t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
      real(r8), intent(inout) :: tl_u(LBi:UBi,LBj:UBj,N(ng),2)
      real(r8), intent(inout) :: tl_v(LBi:UBi,LBj:UBj,N(ng),2)
#   else
      real(r8), intent(inout) :: tl_ubar(LBi:UBi,LBj:UBj,3)
      real(r8), intent(inout) :: tl_vbar(LBi:UBi,LBj:UBj,3)
#   endif
      real(r8), intent(inout) :: tl_zeta(LBi:UBi,LBj:UBj,3)
#  endif
!
!  Local variable declarations.
!
#  ifndef MASKING
      integer :: Imax, Ioff, Jmax, Joff
#  endif
      integer :: i, iadd, is, itrc, j, k

      integer, dimension(5+NT(ng)) :: offset

      real(r8) :: cff, scale

#  include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Extract tangent linear STATE variables from full 1D state vector.
!-----------------------------------------------------------------------
!
!  Determine the index offset for each variable in the state vector.
#  ifdef MASKING
!  Notice that in Land/Sea masking application the state vector only
!  contains water points to avoid large null space.
#  endif
!
#  ifdef SOLVE3D
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+NwaterR(ng)
      offset(isVvel)=offset(isUvel)+NwaterU(ng)*N(ng)
      iadd=NwaterV(ng)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=NwaterR(ng)*N(ng)
      END DO
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVvel)=offset(isUvel)+(Lm(ng)+1)*(Mm(ng)+2)*N(ng)
      iadd=(Lm(ng)+2)*(Mm(ng)+1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=(Lm(ng)+2)*(Mm(ng)+2)*N(ng)
      END DO
#    else
      offset(isFsur)=0
      offset(isUvel)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVvel)=offset(isUvel)+(Lm(ng)-1)*Mm(ng)*N(ng)
      iadd=Lm(ng)*(Mm(ng)-1)*N(ng)
      DO itrc=1,NT(ng)
        offset(isTvar(itrc))=offset(isTvar(itrc)-1)+iadd
        iadd=Lm(ng)*Mm(ng)*N(ng)
      END DO
#    endif
#   endif
#  else
#   ifdef MASKING
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+NwaterR(ng)
      offset(isVbar)=offset(isUbar)+NwaterU(ng)
#   else
#    ifdef FULL_GRID
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+(Lm(ng)+2)*(Mm(ng)+2)
      offset(isVbar)=offset(isUbar)+(Lm(ng)+1)*(Mm(ng)+2)
#    else
      offset(isFsur)=0
      offset(isUbar)=offset(isFsur)+Lm(ng)*Mm(ng)
      offset(isVbar)=offset(isUbar)+(Lm(ng)-1)*Mm(ng)
#    endif
#   endif
#  endif
!
!  Extract tangent linear free-surface.
!
#  ifndef MASKING
#   ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=0
#   else
      Imax=Lm(ng)
      Ioff=0
      Joff=1
#   endif
#  endif
#  ifdef ENERGYNORM_SCALE
      scale=1.0_r8/SQRT(0.5_r8*g*rho0)
#  else
      scale=1.0_r8
#  endif
      DO j=JR_RANGE
        DO i=IR_RANGE
#  ifdef MASKING
          IF (rmask(i,j).gt.0.0_r8) THEN
            is=IJwaterR(i,j)+offset(isFsur)
            tl_zeta(i,j,kstp)=scale*state(is)
          ELSE
            tl_zeta(i,j,kstp)=0.0_r8
          END IF
#  else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isFsur)
          tl_zeta(i,j,kstp)=scale*state(is)
#  endif
        END DO
      END DO

#  ifndef SOLVE3D
!
!  Extract tangent linear 2D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-1
      Ioff=1
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JR_RANGE
        DO i=IU_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i-1,j)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (umask(i,j).gt.0.0_r8) THEN
            is=IJwaterU(i,j)+offset(isUbar)
            tl_ubar(i,j,kstp)=scale*state(is)
          ELSE
            tl_ubar(i,j,kstp)=0.0_r8
          END IF
#   else
          is=(i-Ioff)+(j-Joff)*Imax+offset(isUbar)
          tl_ubar(i,j,kstp)=scale*state(is)
#   endif
        END DO
      END DO
!
!  Extract tangent linear 2D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Ioff=0
      Joff=2
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO j=JV_RANGE
        DO i=IR_RANGE
#   ifdef ENERGYNORM_SCALE
          scale=1.0_r8/SQRT(cff*(h(i,j-1)+h(i,j)))
#   endif
#   ifdef MASKING
          IF (vmask(i,j).gt.0.0_r8) THEN
            is=IJwaterV(i,j)+offset(isVbar)
            tl_vbar(i,j,kstp)=scale*state(is)
          ELSE
            tl_vbar(i,j,kstp)=0.0_r8
          END IF
#   else
          is=(i+Ioff)+(j-Joff)*Imax+offset(isVbar)
          tl_vbar(i,j,kstp)=scale*state(is)
#   endif
        END DO
      END DO

#  else
!
!  Extract tangent linear 3D U-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+1
      Jmax=Mm(ng)+2
      Ioff=0
      Joff=0
#    else
      Imax=Lm(ng)-1
      Jmax=Mm(ng)
      Ioff=1
      Joff=1
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterU(ng)+offset(isUvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isUvel)
#   endif
        DO j=JR_RANGE
          DO i=IU_RANGE
#   ifdef MASKING
            IF (umask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
              is=IJwaterU(i,j)+iadd
              tl_u(i,j,k,nstp)=scale*state(is)
            ELSE
              tl_u(i,j,k,nstp)=0.0_r8
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i-1,j,k)+Hz(i,j,k)))
#    endif
            is=(i-Ioff)+(j-Joff)*Imax+iadd
            tl_u(i,j,k,nstp)=scale*state(is)
#   endif
          END DO
        END DO
      END DO
!
!  Extract tangent linear 3D V-velocity.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+1
      Ioff=1
      Joff=1
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)-1
      Ioff=0
      Joff=2
#    endif
#   endif
#   ifdef ENERGYNORM_SCALE
      cff=0.25_r8*rho0
#   else
      scale=1.0_r8
#   endif
      DO k=1,N(ng)
#   ifdef MASKING
        iadd=(k-1)*NwaterV(ng)+offset(isVvel)
#   else
        iadd=(k-1)*Imax*Jmax+offset(isVvel)
#   endif
        DO j=JV_RANGE
          DO i=IR_RANGE
#   ifdef MASKING
            IF (vmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
              is=IJwaterV(i,j)+iadd
              tl_v(i,j,k,nstp)=scale*state(is)
            ELSE
              tl_v(i,j,k,nstp)=0.0_r8
            END IF
#   else
#    ifdef ENERGYNORM_SCALE
            scale=1.0_r8/SQRT(cff*(Hz(i,j-1,k)+Hz(i,j,k)))
#    endif
            is=(i+Ioff)+(j-Joff)*Imax+iadd
            tl_v(i,j,k,nstp)=scale*state(is)
#   endif
          END DO
        END DO
      END DO
!
!  Extract tangent linear tracers variables. For now, use salinity scale
!  for passive tracers.
!
#   ifndef MASKING
#    ifdef FULL_GRID
      Imax=Lm(ng)+2
      Jmax=Mm(ng)+2
      Ioff=1
      Joff=0
#    else
      Imax=Lm(ng)
      Jmax=Mm(ng)
      Ioff=0
      Joff=1
#    endif
#   endif
      DO itrc=1,NT(ng)
#   ifdef ENERGYNORM_SCALE
        IF (itrc.eq.itemp) THEN
          cff=0.5_r8*rho0*Tcoef(ng)*Tcoef(ng)*g*g/bvf_bak
        ELSE IF (itrc.eq.isalt) THEN
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        ELSE
          cff=0.5_r8*rho0*Scoef(ng)*Scoef(ng)*g*g/bvf_bak
        END IF
#   else
        scale=1.0_r8
#   endif
        DO k=1,N(ng)
#   ifdef MASKING
          iadd=(k-1)*NwaterR(ng)+offset(isTvar(itrc))
#   else
          iadd=(k-1)*Imax*Jmax+offset(isTvar(itrc))
#   endif
          DO j=JR_RANGE
            DO i=IR_RANGE
#   ifdef MASKING
              IF (rmask(i,j).gt.0.0_r8) THEN
#    ifdef ENERGYNORM_SCALE
                scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
                is=IJwaterR(i,j)+iadd
                tl_t(i,j,k,nstp,itrc)=scale*state(is)
              ELSE
                tl_t(i,j,k,nstp,itrc)=0.0_r8
              END IF
#   else
#    ifdef ENERGYNORM_SCALE
              scale=1.0_r8/SQRT(cff*Hz(i,j,k))
#    endif
              is=(i+Ioff)+(j-Joff)*Imax+iadd
              tl_t(i,j,k,nstp,itrc)=scale*state(is)
#   endif
            END DO
          END DO
        END DO
      END DO
#  endif

      RETURN
      END SUBROUTINE tl_unpack_tile
# endif

# ifdef SO_SEMI
#  ifdef SO_SEMI_WHITE
!
      SUBROUTINE so_semi_white (ng, tile, Mstr, Mend, state, ad_state)
!
!=======================================================================
!                                                                      !
!  This routine computes a new stochastic optimals perturbation vector !
!  (seminorm estimation) assuming white noise forcing using ARPACK.    !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
      USE mod_storage
#   ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#   endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#   ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
      real(r8), intent(out) :: ad_state(Mstr:)
#   else
      real(r8), intent(in) :: state(Mstr:Mend)
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#   endif
!
!  Local variable declarations.
!
      integer :: NSUB, is, rec

      real(r8) :: SOnorm, my_SOnorm, my_TRnorm
      real(r8) :: SOnorm1, my_SOnorm1
      real(r8) :: cff, cff1, cff2

#   ifdef DISTRIBUTE
      real(r8), dimension(3) :: buffer

      character (len=3), dimension(3) :: op_handle
#   endif

#   include "tile.h"
!
!-----------------------------------------------------------------------
!  Compute seminorm, stochastic optimals adjoint perturbation vector.
!-----------------------------------------------------------------------
!
!  Initialize adjoint state vector.
!
      DO is=Mstr,Mend
        ad_state(is)=0.0_r8
      END DO
!
!  Sum over all adjoint surface forcing records available in "so_state'.
!
      IF (Master) THEN
        WRITE (stdout,'(/)')
      END IF
      my_TRnorm=0.0_r8
!
      DO rec=1,Nsemi(ng)
!
!  Compute normalization factor.
!
        cff=REAL((nADJ(ng)-1)*(2*nADJ(ng)-1),r8)/REAL(6*nADJ(ng),r8)
        cff1=1.0_r8+cff
        cff2=0.5_r8*REAL((nADJ(ng)-1))-cff
!
        my_SOnorm=0.0_r8
        my_SOnorm1=0.0_r8
        DO is=Mstr,Mend
          my_SOnorm=my_SOnorm+                                          &
     &              STORAGE(ng)%so_state(is,rec)*state(is)
        END DO
!
        IF (rec.ne.Nsemi(ng)) THEN
          DO is=Mstr,Mend
            my_SOnorm1=my_SOnorm1+                                      &
     &                 STORAGE(ng)%so_state(is,rec+1)*state(is)
            my_TRnorm=my_TRnorm+                                        &
     &                cff1*STORAGE(ng)%so_state(is,rec)*                &
     &                     STORAGE(ng)%so_state(is,rec)+                &
     &                2.0_r8*cff2*STORAGE(ng)%so_state(is,rec  )*       &
     &                            STORAGE(ng)%so_state(is,rec+1)+       &
     &                cff*STORAGE(ng)%so_state(is,rec+1)*               &
     &                    STORAGE(ng)%so_state(is,rec+1)
          END DO
        ELSE
          DO is=Mstr,Mend
            my_TRnorm=my_TRnorm+                                        &
     &                STORAGE(ng)%so_state(is,rec)*                     &
     &                STORAGE(ng)%so_state(is,rec)
          END DO
        END IF
!
!  Global reduction of normalization factor.
!
        IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                      &
     &      DOMAIN(ng)%NorthEast_Corner(tile)) THEN
          NSUB=1                         ! non-tiled application
        ELSE
          NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
        END IF
!$OMP CRITICAL (SO_NORM)
        IF (tile_count.eq.0) THEN
          SOnorm=0.0_r8
          SOnorm1=0.0_r8
          IF (rec.eq.1) THEN
            TRnorm(ng)=0.0_r8
          END IF
        END IF
        SOnorm=SOnorm+my_SOnorm
        SOnorm1=SOnorm1+my_SOnorm1
        tile_count=tile_count+1
        IF (tile_count.eq.NSUB) THEN
          tile_count=0
#   ifdef DISTRIBUTE
          buffer(1)=SOnorm
          buffer(2)=SOnorm1
          op_handle(1)='SUM'
          op_handle(2)='SUM'
          CALL mp_reduce (ng, iADM, 2, buffer, op_handle)
          SOnorm=buffer(1)
          SOnorm1=buffer(2)
#   endif
        END IF
!$OMP END CRITICAL (SO_NORM)
!
!  Report normalization factors.
!
        IF (Master) THEN
          WRITE (stdout,10) rec, SOnorm, SOnorm1
 10       FORMAT (3x,'Rec = ',i2.2,2x,'SOnorm = ',1p,e15.8,0p,          &
     &            2x,'SOnorm1 = ',1p,e15.8)
        END IF
!
!  Compute new perturbation vector.
!
        IF (rec.ne.Nsemi(ng)) THEN
          DO is=Mstr,Mend
            ad_state(is)=ad_state(is)+                                  &
     &                   cff1*SOnorm *STORAGE(ng)%so_state(is,rec  )+   &
     &                   cff2*SOnorm1*STORAGE(ng)%so_state(is,rec  )+   &
     &                   cff2*SOnorm *STORAGE(ng)%so_state(is,rec+1)+   &
     &                   cff *SOnorm1*STORAGE(ng)%so_state(is,rec+1)
          END DO
        ELSE
          DO is=Mstr,Mend
            ad_state(is)=ad_state(is)+                                  &
     &                   SOnorm*STORAGE(ng)%so_state(is,rec)
          END DO
        END IF
      END DO
!
!  Global reduction of normalization factor, TRnorm.
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (TR_NORM)
      IF (tile_count.eq.0) THEN
        TRnorm(ng)=0.0_r8
      END IF
      TRnorm(ng)=TRnorm(ng)+my_TRnorm
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#   ifdef DISTRIBUTE
        op_handle(1)='SUM'
        CALL mp_reduce (ng, iADM, 1, TRnorm(ng), op_handle(1))
#   endif
      END IF
!$OMP END CRITICAL (TR_NORM)

      RETURN
      END SUBROUTINE so_semi_white

#  else
!
      SUBROUTINE so_semi_red (ng, tile, Mstr, Mend, state, ad_state)
!
!=======================================================================
!                                                                      !
!  This routine computes a new stochastic optimals perturbation vector !
!  (seminorm estimation) assuming red noise forcing using ARPACK.      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_iounits
      USE mod_storage
      USE mod_scalars
#   ifdef DISTRIBUTE
!
      USE distribute_mod, ONLY : mp_reduce
#   endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: Mstr, Mend
#   ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: state(Mstr:)
      real(r8), intent(out) :: ad_state(Mstr:)
#   else
      real(r8), intent(in) :: state(Mstr:Mend)
      real(r8), intent(out) :: ad_state(Mstr:Mend)
#   endif
!
!  Local variable declarations.
!
      integer :: NSUB, is, ntAD, ntTL, rec, rec1

      real(r8) :: SOnorm, my_TRnorm

      real(r8), dimension(Nsemi(ng)) :: Bcoef
      real(r8), dimension(Nsemi(ng)) :: SOdotprod
      real(r8), dimension(Nsemi(ng)) :: my_dotprod

#   ifdef DISTRIBUTE
      character (len=3), dimension(Nsemi(ng)) :: op_handle
#   endif

#   include "tile.h"
!
!-----------------------------------------------------------------------
!  Compute seminorm, stochastic optimals adjoint perturbation vector.
!-----------------------------------------------------------------------
!
!  Initialize adjoint state vector.
!
      DO is=Mstr,Mend
        ad_state(is)=0.0_r8
      END DO
!
!  First compute the dot-products.
!
      DO rec=1,Nsemi(ng)
        my_dotprod(rec)=0.0_r8
        DO is=Mstr,Mend
          my_dotprod(rec)=my_dotprod(rec)+                              &
     &                    STORAGE(ng)%so_state(is,rec)*state(is)
        END DO
      END DO
!
!  Global reduction of dot products.
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                         ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)     ! tiled application
      END IF
!$OMP CRITICAL (SO_DOT)
      IF (tile_count.eq.0) THEN
        DO rec=1,Nsemi(ng)
          SOdotprod(rec)=0.0_r8
        END DO
      END IF
      DO rec=1,Nsemi(ng)
        SOdotprod(rec)=SOdotprod(rec)+my_dotprod(rec)
      END DO
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#   ifdef DISTRIBUTE
        DO rec=1,Nsemi(ng)
          op_handle(rec)='SUM'
        END DO
        CALL mp_reduce (ng, iADM, Nsemi(ng), SOdotprod, op_handle)
#   endif
      END IF
!$OMP END CRITICAL (SO_DOT)
!
! Second, loop over time twice allowing for the decorrelation due to the
! red noise AR(1) process with assumed decorrelation time SOdecay.
!
      IF (Master) THEN
        WRITE (stdout,'(/)')
      END IF
      my_TRnorm=0.0_r8
!
      DO rec=1,Nsemi(ng)
        ntAD=(rec-1)*nADJ(ng)+1
        SOnorm=0.0_r8
        DO rec1=1,Nsemi(ng)
          ntTL=(rec1-1)*nADJ(ng)+1
          CALL sp_bcoef (ng, ntAD, ntTL, Bcoef(rec1))
          SOnorm=SOnorm+Bcoef(rec1)*SOdotprod(rec1)
          DO is=Mstr,Mend
            my_TRnorm=my_TRnorm+                                        &
     &                STORAGE(ng)%so_state(is,rec )*Bcoef(rec1)*        &
     &                STORAGE(ng)%so_state(is,rec1)
          END DO
        END DO
!
!  Report normalization factors.
!
        IF (Master) THEN
          WRITE (stdout,10) rec, SOdotprod(rec), Bcoef(rec), SOnorm
 10       FORMAT (1x,'Rec = ',i2.2,1x,'SOdotprod = ',1p,e13.6,0p,       &
     &            1x,'Bcoef = ',1p,e13.6,0p,1x,'SOnorm = ',1p,e13.6)
        END IF
!
!  Compute new perturbation vector.
!
        DO is=Mstr,Mend
          ad_state(is)=ad_state(is)+                                    &
     &                 SOnorm*STORAGE(ng)%so_state(is,rec)
        END DO
      END DO
!
!  Global reduction of normalization factor, TRnorm.
!
      IF (DOMAIN(ng)%SouthWest_Corner(tile).and.                        &
     &    DOMAIN(ng)%NorthEast_Corner(tile)) THEN
        NSUB=1                           ! non-tiled application
      ELSE
        NSUB=NtileX(ng)*NtileE(ng)       ! tiled application
      END IF
!$OMP CRITICAL (TR_NORM)
      IF (tile_count.eq.0) THEN
        TRnorm(ng)=0.0_r8
      END IF
      TRnorm(ng)=TRnorm(ng)+my_TRnorm
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
#   ifdef DISTRIBUTE
        op_handle(1)='SUM'
        CALL mp_reduce (ng, iADM, 1, TRnorm(ng), op_handle(1))
#   endif
      END IF
!$OMP END CRITICAL (TR_NORM)

      RETURN
      END SUBROUTINE so_semi_red
#  endif
# endif

# ifdef SO_SEMI
!
      SUBROUTINE sp_bcoef (ng, ntAD, ntTL, Bcoef)
!
!=======================================================================
!                                                                      !
!  This routine is used to compute  red noise stochastic processes     !
!  time-lagged  coefficient,  Bcoef,  used  to  evaluate  discrete     !
!  double-time integrals.  Currently, a discrete-time Markov chain     !
!  model is assumed with autoregressive order-one processes, AR(1).    !
!  Notice that the routine SP_ACOEF is called to compute the inner     !
!  integral.                                                           !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ntAD, ntTL

      real(r8), intent(out):: Bcoef
!
!  Local variable declarations.
!
      integer :: i, it1, it2

      real(r8) :: Acoef, Acoef1, Acoef2, df1, rov

!
!-----------------------------------------------------------------------
!  Compute red noise stochastic process time-lagged coefficient to
!  evaluate discrete double time-integrals. Currently, only Markov
!  processes, AR(1), are considered.
!-----------------------------------------------------------------------
!
!  Here, ntAD is the current model timestep and ntTL is the timestep
!  associated with forcing.
!
      rov=1.0_r8/REAL(nADJ(ng),r8)
!
      IF ((ntAD.gt.1).and.(ntAD.lt.ntimes(ng)+1)) THEN
        it1=ntAD
        it2=ntAD-nADJ(ng)
        CALL sp_acoef (ng, it1, ntTL, Acoef)
        Bcoef=Acoef
        DO i=1,nADJ(ng)-1
          CALL sp_acoef (ng, it1+i, ntTL, Acoef1)
          CALL sp_acoef (ng, it2+i, ntTL, Acoef2)
          df1=REAL(i,r8)*rov
          Bcoef=Bcoef+(1.0_r8-df1)*Acoef1+df1*Acoef2
        END DO
      ELSE IF (ntAD.eq.1) THEN
        CALL sp_acoef (ng, 1, ntTL, Acoef)
        Bcoef=Acoef
        DO i=1,nADJ(ng)-1
          CALL sp_acoef (ng, 1+i, ntTL, Acoef1)
          df1=REAL(i,r8)*rov
          Bcoef=Bcoef+(1.0_r8-df1)*Acoef1
        END DO
      ELSE IF (ntAD.eq.ntimes(ng)+1) THEN
        CALL sp_acoef (ng, ntimes(ng)+1, ntTL, Acoef)
        Bcoef=Acoef
        DO i=1,nADJ(ng)-1
          CALL sp_acoef (ng, ntimes(ng)+1-nADJ(ng)+i, ntTL, Acoef1)
          df1=REAL(i,r8)*rov
          Bcoef=Bcoef+df1*Acoef1
        END DO
      END IF

      RETURN
      END SUBROUTINE sp_bcoef
!
      SUBROUTINE sp_acoef (ng, ntAD, ntTL, Acoef)
!
!=======================================================================
!                                                                      !
!  This routine is used to compute red noise stochastic processes      !
!  time-lagged coefficient, Acoef,  used  to  evaluate inner time      !
!  integral.  Currently,  a  discrete-time Markov chain model  is      !
!  assumed with autoregressive order-one processes, AR(1). Notice      !
!  that function SP_AUTOC is used to set autocorrelation model.        !
!                                                                      !
!=======================================================================
!
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, ntAD, ntTL

      real(r8), intent(out):: Acoef
!
!  Local variable declarations.
!
      integer :: i, idf1, idf2, idf4

      real(r8) :: df3, rov
!
!-----------------------------------------------------------------------
!  Compute red noise stochastic process time-lagged coefficients to
!  evaluate discrete inner time-integral. Currently, only Markov
!  processes, AR(1), are considered.
!-----------------------------------------------------------------------
!
!  Here, ntAD is the current timestep corresponding to time when
!  solution is saved.
!
      rov=1.0_r8/REAL(nADJ(ng),r8)
      IF ((ntTL.gt.1).and.(ntTL.lt.ntimes(ng)+1)) THEN
        Acoef=0.0_r8
        DO i=1,nADJ(ng)-1
          idf1=IABS(ntAD-ntTL-i)+1
          idf2=IABS(ntAD-(ntTL-nADJ(ng))-i)+1
          df3=REAL(i,r8)*rov
          Acoef=Acoef+sp_autoc(ng,idf1)*(1.0_r8-df3)+                   &
     &                sp_autoc(ng,idf2)*df3
        END DO
        idf4=IABS(ntAD-ntTL)+1
        Acoef=Acoef+sp_autoc(ng,idf4)
      ELSE IF (ntTL.eq.1) THEN
        Acoef=0.0_r8
        DO i=1,nADJ(ng)-1
          idf1=IABS(ntAD-1-i)+1
          df3=REAL(i,r8)*rov
          Acoef=Acoef+sp_autoc(ng,idf1)*(1.0_r8-df3)
        END DO
        idf4=IABS(ntAD-1)+1
        Acoef=Acoef+sp_autoc(ng,idf4)
      ELSE IF (ntTL.eq.ntimes(ng)+1) THEN
        Acoef=0.0_r8
        DO i=1,nADJ(ng)-1
          idf2=IABS(ntAD-ntimes(ng)-1+nADJ(ng)-i)+1
          df3=REAL(i,r8)*rov
          Acoef=Acoef+sp_autoc(ng,idf2)*df3
        END DO
        idf4=IABS(ntAD-ntimes(ng)-1)+1
        Acoef=Acoef+sp_autoc(ng,idf4)
      END IF

      RETURN
      END SUBROUTINE sp_acoef
!
      FUNCTION sp_autoc (ng, idf)
!
!=======================================================================
!                                                                      !
!  This routine is used to compute red noise stochastic processes      !
!  autocorrelation model.  Notice that only  AR(1)  processes are      !
!  considered. However, other models can be easily implemented in      !
!  terms of look tables.                                               !
!
!=======================================================================
!
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, idf
!
!  Function result.
!
      real(r8) :: sp_autoc
!
!-----------------------------------------------------------------------
!  Set autocorrelation model.
!-----------------------------------------------------------------------
#  ifdef SO_NON_AR1
!
!  Use a user-defined temporal decorrelation function such as in the
!  form of a look-up table computed from actual data.
!
      sp_autoc=0.0_r8
#  else
!
!  Assume an AR(1) process with decorrelation time SO_decay.
!
      sp_autoc=EXP(-ABS(REAL(idf-1,r8))*dt(ng)/SO_decay(ng))
#  endif
!
      RETURN
      END FUNCTION sp_autoc
# endif
#endif

#undef IR_RANGE
#undef IU_RANGE
#undef JR_RANGE
#undef JV_RANGE

      END MODULE packing_mod
