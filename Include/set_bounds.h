!
!---------------------------------------------------------------------
!  Compute bound indices over a particular domain partition or tile
!  for RHO-, U-, and V-variables.
!---------------------------------------------------------------------
!
      if (WESTERN_EDGE) then
#ifdef EW_PERIODIC
        IstrR=Istr-2
        IstrU=Istr
#else
        IstrR=Istr-1
        IstrU=Istr+1
#endif
      else
        IstrR=Istr
        IstrU=Istr
      endif
      if (EASTERN_EDGE) then
#ifdef EW_PERIODIC
        IendR=Iend+2
#else
        IendR=Iend+1
#endif
      else
        IendR=Iend
      endif
      if (SOUTHERN_EDGE) then
#ifdef NS_PERIODIC
        JstrR=Jstr-2
        JstrV=Jstr
#else
        JstrR=Jstr-1
        JstrV=Jstr+1
#endif
      else
        JstrR=Jstr
        JstrV=Jstr
      endif
      if (NORTHERN_EDGE) then
#ifdef NS_PERIODIC
        JendR=Jend+2
#else
        JendR=Jend+1
#endif
      else
        JendR=Jend
      endif
