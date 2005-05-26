#include "cppdefs.h"
      MODULE utility_mod
!
!=======================================================================
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!=======================================================================
!                                                                      !
!  This module contains several all purpuse generic routines:          !
!                                                                      !
!  Routines:                                                           !
!                                                                      !
!    hindices      Finds model grid cell for any datum.                !
!    try_range     Binary search of model grid cell for any datum.     !
!    inside        Closed polygon datum search.                        !
!    nrng          Gaussian random number generator.                   !
!    urng          Uniform random number generator.                    !
!                                                                      !
!=======================================================================
!
      USE mod_kinds

      implicit none

      CONTAINS

      SUBROUTINE hindices (ng, LBi, UBi, LBj, UBj,                      &
     &                     Istr, Iend, Jstr, Jend,                      &
     &                     angler, Xgrd, Ygrd,                          &
     &                     IJspv, Npos, Xpos, Ypos, Ipos, Jpos)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  Given position vectors Xpos and Ypos of size Npos,  this routine    !
!  finds the corresponding indices Ipos and Jpos of the  model grid    !
!  (Xgrd,Ygrd) cell containing each requested position.                !
!                                                                      !
!  On Input:                                                           !
!                                                                      !
!     ng         Nested grid number.                                   !
!     LBi        I-dimension Lower bound.                              !
!     UBi        I-dimension Upper bound.                              !
!     LBj        J-dimension Lower bound.                              !
!     UBj        J-dimension Upper bound.                              !
!     Istr       Starting model grid I-index to search.                !
!     Iend       Ending   model grid I-index to search.                !
!     Jstr       Starting model grid J-index to search.                !
!     Jend       Ending   model grid J-index to search.                !
!     angler     Angle (radians) between XI-axis and true EAST.        !
!     Xgrd       Geographical X-grid coordinate (usually, longitude).  !
!     Ygrd       Geographycal Y-grid coordinate (usually, latitude).   !
!     IJspv      (I,J) unbounded special value.                        !
!     Npos       Number of points to process.                          !
!     Xpos       Geograpical X-coordinate to process.                  !
!     Ypos       Geograpical Y-coordinate to process.                  !
!                                                                      !
!  On Output:                                                          !
!                                                                      !
!     Ipos       Fractional I-grid coordinate.                         !
!     Jpos       Fractional I-grid coordinate.                         !
!                                                                      !
!  Calls:    Try_Range                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
      integer, intent(in) :: Istr, Iend, Jstr, Jend, Npos

      real(r8), intent(in) :: IJspv

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: Xgrd(LBi:,LBj:)
      real(r8), intent(in) :: Ygrd(LBi:,LBj:)
      real(r8), intent(in) :: Xpos(:), Ypos(:)
      real(r8), intent(out) :: Ipos(:), Jpos(:)
#else
      real(r8), intent(in) :: angler(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Xgrd(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Ygrd(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Xpos(Npos), Ypos(Npos)
      real(r8), intent(out) :: Ipos(Npos), Jpos(Npos)
#endif
!
!  Local variable declarations.
!
      logical :: found

      integer :: Imax, Imin, Jmax, Jmin, i0, j0, k

      real(r8) :: aa2, ang, bb2, diag2, dx, dy, phi
      real(r8) :: xfac, xpp, yfac, ypp
!
!-----------------------------------------------------------------------
!  Determine grid cell indices containing requested position points.
!  Then, interpolate to fractional cell position.
!-----------------------------------------------------------------------
!
!  Initialize all indices.
!
      DO k=1,Npos
        Ipos(k)=IJspv
        Jpos(k)=IJspv
      END DO
!
!  Check each position to find if it falls inside the whole domain.
!  Once it is stablished that it inside, find the exact cell to which
!  it belongs by successively dividing the domain by a half (binary
!  search).
!
      DO k=1,Npos
        found=try_range(ng, LBi, UBi, LBj, UBj,                         &
     &                  Xgrd, Ygrd,                                     &
     &                  Istr, Iend, Jstr, Jend,                         &
     &                  Xpos(k), Ypos(k))
        if (found) THEN
          Imin=Istr
          Imax=Iend
          Jmin=Jstr
          Jmax=Jend
          DO while (((Imax-Imin).gt.1).or.((Jmax-Jmin).gt.1))
            IF ((Imax-Imin).gt.1) THEN
              i0=(Imin+Imax)/2
              found=try_range(ng, LBi, UBi, LBj, UBj,                   &
     &                        Xgrd, Ygrd,                               &
     &                        Imin, i0, Jmin, Jmax,                     &
     &                        Xpos(k), Ypos(k))
              IF (found) THEN
                Imax=i0
              ELSE
                Imin=i0
              END IF
            END IF
            IF ((Jmax-Jmin).gt.1) THEN
              j0=(Jmin+Jmax)/2
              found=try_range(ng, LBi, UBi, LBj, UBj,                   &
     &                        Xgrd, Ygrd,                               &
     &                        Imin, Imax, Jmin, j0,                     &
     &                        Xpos(k), Ypos(k))
              IF (found) THEN
                Jmax=j0
              ELSE
                Jmin=j0
              END IF
            END IF
          END DO
!
!  Knowing the correct cell, calculate the exact indices, accounting
!  for a possibly rotated grid.  If spherical, convert all positions
!  to meters first.
!
          IF (spherical) THEN
            yfac=Eradius*deg2rad
            xfac=yfac*COS(Ypos(k)*deg2rad)
            xpp=(Xpos(k)-Xgrd(Imin,Jmin))*xfac
            ypp=(Ypos(k)-Ygrd(Imin,Jmin))*yfac
          ELSE
            xfac=1.0_r8
            yfac=1.0_r8
            xpp=Xpos(k)-Xgrd(Imin,Jmin)
            ypp=Ypos(k)-Ygrd(Imin,Jmin)
          END IF
!
!  Use Law of Cosines to get cell parallelogram "shear" angle.
!
          diag2=(Xgrd(Imin+1,Jmin)-Xgrd(Imin,Jmin+1))**2+               &
     &          (Ygrd(Imin+1,Jmin)-Ygrd(Imin,Jmin+1))**2
          aa2=(Xgrd(Imin,Jmin)-Xgrd(Imin+1,Jmin))**2+                   &
     &        (Ygrd(Imin,Jmin)-Ygrd(Imin+1,Jmin))**2
          bb2=(Xgrd(Imin,Jmin)-Xgrd(Imin,Jmin+1))**2+                   &
     &        (Ygrd(Imin,Jmin)-Ygrd(Imin,Jmin+1))**2
          phi=ASIN((diag2-aa2-bb2)/(2.0_r8*SQRT(aa2*bb2)))
!
!  Transform float position into curvilinear coordinates. Assume the
!  cell is rectanglar, for now.
!
          ang=angler(Imin,Jmin)
          dx=xpp*COS(ang)+ypp*SIN(ang)
          dy=ypp*COS(ang)-xpp*SIN(ang)
!
!  Correct for parallelogram.
!
          dx=dx+dy*TAN(phi)
          dy=dy/COS(phi)
!
!  Scale with cell side lengths to translate into cell indices.
!
          dx=MIN(MAX(0.0_r8,dx/SQRT(aa2)/xfac),1.0_r8)
          dy=MIN(MAX(0.0_r8,dy/SQRT(bb2)/yfac),1.0_r8)
          Ipos(k)=REAL(Imin,r8)+dx
          Jpos(k)=REAL(Jmin,r8)+dy
        END IF
      END DO
      RETURN
      END SUBROUTINE hindices

      LOGICAL FUNCTION try_range (ng, LBi, UBi, LBj, UBj, Xgrd, Ygrd,   &
     &                            Imin, Imax, Jmin, Jmax, Xo, Yo)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  Given a grided domain with matrix coordinates Xgrd and Ygrd, this   !
!  function finds if the point (Xo,Yo)  is inside the box defined by   !
!  the requested corners (Imin,Jmin) and (Imax,Jmax). It will return   !
!  logical switch  try_range=.TRUE.  if (Xo,Yo) is inside, otherwise   !
!  it will return false.                                               !
!                                                                      !
!  Calls:   inside                                                     !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, LBi, UBi, LBj, UBj
      integer, intent(in) :: Imin, Imax, Jmin, Jmax

#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: Xgrd(LBi:,LBj:)
      real(r8), intent(in) :: Ygrd(LBi:,LBj:)
#else
      real(r8), intent(in) :: Xgrd(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Ygrd(LBi:UBi,LBj:UBj)
#endif

      real(r8), intent(in) :: Xo, Yo
!
!  Local variable declarations.
!
      integer ::  Nb, i, j, shft, ic

      real(r8), dimension(2*(Jmax-Jmin+Imax-Imin)+1) :: Xb, Yb
!
!-----------------------------------------------------------------------
!  Define closed polygon.
!-----------------------------------------------------------------------
!
!  Note that the last point (Xb(Nb),Yb(Nb)) does not repeat first
!  point (Xb(1),Yb(1)).  Instead, in function inside, it is implied
!  that the closing segment is (Xb(Nb),Yb(Nb))-->(Xb(1),Yb(1)). In
!  fact, function inside sets Xb(Nb+1)=Xb(1) and Yb(Nb+1)=Yb(1).
!
      Nb=2*(Jmax-Jmin+Imax-Imin)
      shft=1-Imin
      DO i=Imin,Imax-1
        Xb(i+shft)=Xgrd(i,Jmin)
        Yb(i+shft)=Ygrd(i,Jmin)
      END DO
      shft=1-Jmin+Imax-Imin
      DO j=Jmin,Jmax-1
        Xb(j+shft)=Xgrd(Imax,j)
        Yb(j+shft)=Ygrd(Imax,j)
      END DO
      shft=1+Jmax-Jmin+2*Imax-Imin
      DO i=Imax,Imin+1,-1
        Xb(shft-i)=Xgrd(i,Jmax)
        Yb(shft-i)=Ygrd(i,Jmax)
      END DO
      shft=1+2*Jmax-Jmin+2*(Imax-Imin)
      DO j=Jmax,Jmin+1,-1
        Xb(shft-j)=Xgrd(Imin,j)
        Yb(shft-j)=Ygrd(Imin,j)
      END DO
!
!-----------------------------------------------------------------------
!  Check if point (Xo,Yo) is inside of the defined polygon.
!-----------------------------------------------------------------------
!
      try_range=inside(Nb, Xb, Yb, Xo, Yo)
      RETURN
      END FUNCTION try_range

      LOGICAL FUNCTION inside (Nb, Xb, Yb, Xo, Yo)
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2005 ROMS/TOMS Group                                  !
!========================================== Alexander F. Shchepetkin ===
!                                                                      !
!  Given the vectors Xb and Yb of size Nb, defining the coordinates    !
!  of a closed polygon,  this function find if the point (Xo,Yo) is    !
!  inside the polygon.  If the point  (Xo,Yo)  falls exactly on the    !
!  boundary of the polygon, it still considered inside.                !
!                                                                      !
!  This algorithm does not rely on the setting of  Xb(Nb)=Xb(1) and    !
!  Yb(Nb)=Yb(1).  Instead, it assumes that the last closing segment    !
!  is (Xb(Nb),Yb(Nb)) --> (Xb(1),Yb(1)).                               !
!                                                                      !
!  Reference:                                                          !
!                                                                      !
!    Reid, C., 1969: A long way from Euclid. Oceanography EMR,         !
!      page 174.                                                       !
!                                                                      !
!  Algorithm:                                                          !
!                                                                      !
!  The decision whether the point is  inside or outside the polygon    !
!  is done by counting the number of crossings from the ray (Xo,Yo)    !
!  to (Xo,-infinity), hereafter called meridian, by the boundary of    !
!  the polygon.  In this counting procedure,  a crossing is counted    !
!  as +2 if the crossing happens from "left to right" or -2 if from    !
!  "right to left". If the counting adds up to zero, then the point    !
!  is outside.  Otherwise,  it is either inside or on the boundary.    !
!                                                                      !
!  This routine is a modified version of the Reid (1969) algorithm,    !
!  where all crossings were counted as positive and the decision is    !
!  made  based on  whether the  number of crossings is even or odd.    !
!  This new algorithm may produce different results  in cases where    !
!  Xo accidentally coinsides with one of the (Xb(k),k=1:Nb) points.    !
!  In this case, the crossing is counted here as +1 or -1 depending    !
!  of the sign of (Xb(k+1)-Xb(k)).  Crossings  are  not  counted if    !
!  Xo=Xb(k)=Xb(k+1).  Therefore, if Xo=Xb(k0) and Yo>Yb(k0), and if    !
!  Xb(k0-1) < Xb(k0) < Xb(k0+1),  the crossing is counted twice but    !
!  with weight +1 (for segments with k=k0-1 and k=k0). Similarly if    !
!  Xb(k0-1) > Xb(k0) > Xb(k0+1), the crossing is counted twice with    !
!  weight -1 each time.  If,  on the other hand,  the meridian only    !
!  touches the boundary, that is, for example, Xb(k0-1) < Xb(k0)=Xo    !
!  and Xb(k0+1) < Xb(k0)=Xo, then the crossing is counted as +1 for    !
!  segment k=k0-1 and -1 for segment k=k0, resulting in no crossing.   !
!                                                                      !
!  Note 1: (Explanation of the logical condition)                      !
!                                                                      !
!  Suppose  that there exist two points  (x1,y1)=(Xb(k),Yb(k))  and    !
!  (x2,y2)=(Xb(k+1),Yb(k+1)),  such that,  either (x1 < Xo < x2) or    !
!  (x1 > Xo > x2).  Therefore, meridian x=Xo intersects the segment    !
!  (x1,y1) -> (x2,x2) and the ordinate of the point of intersection    !
!  is:                                                                 !
!                                                                      !
!                 y1*(x2-Xo) + y2*(Xo-x1)                              !
!             y = -----------------------                              !
!                          x2-x1                                       !
!                                                                      !
!  The mathematical statement that point  (Xo,Yo)  either coinsides    !
!  with the point of intersection or lies to the north (Yo>=y) from    !
!  it is, therefore, equivalent to the statement:                      !
!                                                                      !
!         Yo*(x2-x1) >= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 > 0      !
!  or                                                                  !
!         Yo*(x2-x1) <= y1*(x2-Xo) + y2*(Xo-x1),   if   x2-x1 < 0      !
!                                                                      !
!  which, after noting that  Yo*(x2-x1) = Yo*(x2-Xo + Xo-x1) may be    !
!  rewritten as:                                                       !
!                                                                      !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) >= 0,   if   x2-x1 > 0      !
!  or                                                                  !
!        (Yo-y1)*(x2-Xo) + (Yo-y2)*(Xo-x1) <= 0,   if   x2-x1 < 0      !
!                                                                      !
!  and both versions can be merged into  essentially  the condition    !
!  that (Yo-y1)*(x2-Xo)+(Yo-y2)*(Xo-x1) has the same sign as x2-x1.    !
!  That is, the product of these two must be positive or zero.         !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: Nb

      real(r8), intent(in) :: Xo, Yo

#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: Xb(:), Yb(:)
#else
      real(r8), intent(inout) :: Xb(Nb+1), Yb(Nb+1)
#endif
!
!  Local variable declarations.
!
      integer, parameter :: Nstep =128

      integer :: crossings, i, inc, k, kk, nc

      integer, dimension(Nstep) :: Sindex

      real(r8) :: dx1, dx2, dxy
!
!-----------------------------------------------------------------------
!  Find intersections.
!-----------------------------------------------------------------------
!
!  Set crossings counter and close the contour of the polygon.
!
      crossings=0
      Xb(Nb+1)=Xb(1)
      Yb(Nb+1)=Yb(1)
!
!  The search is optimized.  First select the indices of segments
!  where Xb(k) is different from Xb(k+1) and Xo falls between them.
!  Then, further investigate these segments in a separate loop.
!  Doing it in two stages takes less time because the first loop is
!  pipelined.
!
      DO kk=0,Nb-1,Nstep
        nc=0
        DO k=kk+1,MIN(kk+Nstep,Nb)
          IF (((Xb(k+1)-Xo)*(Xo-Xb(k)).ge.0.0_r8).and.                  &
     &        (Xb(k).ne.Xb(k+1))) THEN
            nc=nc+1
            Sindex(nc)=k
          END IF
        END DO
        DO i=1,nc
          k=Sindex(i)
          IF (Xb(k).ne.Xb(k+1)) THEN
            dx1=Xo-Xb(k)
            dx2=Xb(k+1)-Xo
            dxy=dx2*(Yo-Yb(k))-dx1*(Yb(k+1)-Yo)
            inc=0
            IF ((Xb(k).eq.Xo).and.(Yb(k).eq.Yo)) THEN
              crossings=1
              goto 10
            ELSE IF (((dx1.eq.0.0_r8).and.(Yo.ge.Yb(k  ))).or.          &
     &              ((dx2.eq.0.0_r8).and.(Yo.ge.Yb(k+1)))) THEN
              inc=1
            ELSE IF ((dx1*dx2.gt.0.0_r8).and.                           &
     &              ((Xb(k+1)-Xb(k))*dxy.ge.0.0_r8)) THEN  ! see note 1
              inc=2
            END IF
            IF (Xb(k+1).gt.Xb(k)) THEN
              crossings=crossings+inc
            ELSE
              crossings=crossings-inc
            END IF
          END IF
        END DO
      END DO
!
!  Determine if point (Xo,Yo) is inside of closed polygon.
!
  10  IF (crossings.eq.0) THEN
        inside=.FALSE.
      ELSE
        inside=.TRUE.
      END IF
      RETURN
      END FUNCTION inside

      SUBROUTINE nrng (ix, a, n, ierr)
!
!=======================================================================
!                                                                      !
!  Gaussian random-number generator from the NSWC Library. It calls    !
!  the NSWC uniform random-number generator, URNG.                     !
!                                                                      !
!  Modernised and included in ROMS by Mark Hadfield, NIWA.             !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: n

      integer, intent(inout) :: ix

      integer, intent(out) :: ierr

#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: a(:)
#else
      real(r8), intent(out) :: a(n)
#endif
!
!  Local variable declarations.
!
      integer :: i, m

      real(r8), parameter ::pi2 = 6.2831853071796_r8

      real(r8) :: phi, r
      real(r8) :: temp(1)
!
!-----------------------------------------------------------------------
!  Generate Gaussian random numbers.
!-----------------------------------------------------------------------
!
      CALL urng (ix, a, n, ierr)
!
      IF (ierr.ne.0) RETURN
!
      IF (n.gt.1) THEN
        m=n/2
        m=m+m
        DO i=1,m,2
          r=SQRT(-2.0_r8*LOG(a(i)))
          phi=pi2*a(i+1)
          a(i  )=r*COS(phi)
          a(i+1)=r*SIN(phi)
        END DO
         IF (m.eq.n) RETURN
      END IF
!
      CALL urng (ix, temp, 1, ierr)
!
      r=SQRT(-2.0_r8*LOG(a(n)))
!
      a(n)=r*COS(pi2*temp(1))
!
      RETURN
      END SUBROUTINE nrng

      SUBROUTINE urng (ix, x, n, ierr)
!
!=======================================================================
!                                                                      !
!  Uniform random-number generator from the NSWC Library               !
!                                                                      !
!  Uses the recursion ix = ix*a mod p, where 0 < ix < p                !
!                                                                      !
!  Written by Linus Schrage, University of Chicago. Adapted for NSWC   !
!  Library by A. H. Morris. Modernised & included in ROMS by Mark      !
!  Hadfield, NIWA.                                                     !
!                                                                      !
!=======================================================================
!
!  Imported variable declarations.
!
      integer, intent(in) :: n

      integer, intent(inout) :: ix

      integer, intent(out) :: ierr

#ifdef ASSUMED_SHAPE
      real(r8), intent(out) :: x(:)
#else
      real(r8), intent(out) :: x(n)
#endif
!
!  Local variable declarations.
!    
      integer, parameter :: a = 16807          ! 7^5
      integer, parameter :: b15 = 32768        ! 2^15
      integer, parameter :: b16 = 65536        ! 2^16
      integer, parameter :: p = 2147483647     ! 2^31-1

      integer :: fhi, k, l, leftlo, xalo, xhi

      real(r8), parameter :: s = 0.465661E-09_r8
!
!-----------------------------------------------------------------------
!  Generate random numbers.
!-----------------------------------------------------------------------
!
      IF (n.le.0) THEN
        ierr=1
        RETURN
      END IF
      IF ((ix.le.0).or.(ix.ge.p)) THEN
        ierr=2
        RETURN
      END IF
!
      ierr=0
!
      DO l=1,n
!
! Get 15 high order bits of "ix".
!
        xhi=ix/b16
!
! Get 16 lower bits of ix and multiply with "a".
!
        xalo=(ix-xhi*b16)*a
!
! Get 15 high order bits of the product.
!
        leftlo=xalo/b16
!
! Form the 31 highest bits of "a*ix".
!
        fhi=xhi*a+leftlo
!
! Obtain the overflow past the 31st bit of "a*ix".
!
        k=fhi/b15
!
! Assemble all the parts and presubtract "p". The parentheses are
! essential.
!
        ix=(((xalo-leftlo*b16)-p)+(fhi-k*b15)*b16)+k
!
! Add "p" if necessary.
!
        IF (ix.lt.0) ix=ix+p
!
! Rescale "ix", to interpret it as a value between 0 and 1.
! the scale factor "s" is selected to be as near "1/p" as is
! appropriate in order that the floating value for "ix = 1",
! namely "s", be roughly the same distance from 0 as "(p-1)*s"
! is from 1. The current value for "s" assures us that "x(l)"
! is less than 1 for any floating point arithmetic of 6
! or more digits.
!
         x(l)=REAL(ix,r8)*s
      END DO
      RETURN
      END SUBROUTINE urng

      END MODULE utility_mod
