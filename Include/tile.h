      integer :: Iend, Istr, Jend, Jstr
      integer :: LBi, UBi, LBj, UBj
      integer :: ChunkSizeI, ChunkSizeJ, MarginI, MarginJ
      integer :: Itile, Jtile
!
!  Set horizontal starting and ending indices for parallel domain
!  partitions in the XI- and ETA-directions.
!
      ChunkSizeI=(Lm(ng)+NtileI(ng)-1)/NtileI(ng)
      ChunkSizeJ=(Mm(ng)+NtileJ(ng)-1)/NtileJ(ng)
      MarginI=(NtileI(ng)*ChunkSizeI-Lm(ng))/2
      MarginJ=(NtileJ(ng)*ChunkSizeJ-Mm(ng))/2
      Jtile=tile/NtileI(ng)
      Itile=tile-Jtile*NtileI(ng)
!
      Istr=1+Itile*ChunkSizeI-MarginI
      Iend=Istr+ChunkSizeI-1
      Istr=MAX(Istr,1)
      Iend=MIN(Iend,Lm(ng))
!
      Jstr=1+Jtile*ChunkSizeJ-MarginJ
      Jend=Jstr+ChunkSizeJ-1
      Jstr=MAX(Jstr,1)
      Jend=MIN(Jend,Mm(ng))
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
#ifdef DISTRIBUTE
      IF (Itile.eq.0) THEN
        LBi=LOWER_BOUND_I
      ELSE
        LBi=Istr-GHOST_POINTS
      END IF
      IF (Itile.eq.(NtileI(ng)-1)) THEN
        UBi=UPPER_BOUND_I
      ELSE
        UBi=Iend+GHOST_POINTS
      END IF
      IF (Jtile.eq.0) THEN
        LBj=LOWER_BOUND_J
      ELSE
        LBj=Jstr-GHOST_POINTS
      END IF
      IF (Jtile.eq.(NtileJ(ng)-1)) THEN
        UBj=UPPER_BOUND_J
      ELSE
        UBj=Jend+GHOST_POINTS
      END IF
#else
      LBi=LOWER_BOUND_I
      UBi=UPPER_BOUND_I
      LBj=LOWER_BOUND_J
      UBj=UPPER_BOUND_J
#endif


