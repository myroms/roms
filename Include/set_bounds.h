/*  Compute derived bounds for the loop indices over a subdomain "tile".
**  The extended bounds (labelled by suffix R) are designed to cover
**  also the outer ghost points, if the subdomain "tile" is adjacent to
**  to the physical boundary. Notice that IstrR, IendR, JstrR, JendR
**  tile bounds computed here DO NOT COVER ghost points associated with
**  periodic boundaries (if any) or then computational margins of MPI
**  subdomains.
**
**  It also computes loop-bounds for U- and V-type variables which
**  belong to the interior of the computational domain. These are
**  labelled by suffixes U,V and they step one grid point inward from
**  the side of the subdomain adjacent to the physical boundary.
**  Conversely, for an internal subdomain which does not include a
**  segments of the physical boundary, all bounds with suffixes R,U,V
**  are set to the same values of corresponding non-suffixed bounds.
*/ 

!
!-----------------------------------------------------------------------
!  Compute lower and upper bounds over a particular domain partition or
!  tile for RHO-, U-, and V-variables.
!-----------------------------------------------------------------------
!
      IF (WESTERN_EDGE) THEN
#ifdef EW_PERIODIC
        IstrR=Istr
        IstrU=Istr
#else
        IstrR=Istr-1
        IstrU=Istr+1
#endif
      ELSE
        IstrR=Istr
        IstrU=Istr
      END IF
      IF (EASTERN_EDGE) THEN
#ifdef EW_PERIODIC
        IendR=Iend
#else
        IendR=Iend+1
#endif
      ELSE
        IendR=Iend
      END IF
      IF (SOUTHERN_EDGE) THEN
#ifdef NS_PERIODIC
        JstrR=Jstr
        JstrV=Jstr
#else
        JstrR=Jstr-1
        JstrV=Jstr+1
#endif
      ELSE
        JstrR=Jstr
        JstrV=Jstr
      END IF
      IF (NORTHERN_EDGE) THEN
#ifdef NS_PERIODIC
        JendR=Jend
#else
        JendR=Jend+1
#endif
      ELSE
        JendR=Jend
      END IF
