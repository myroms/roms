/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2025 The ROMS Group             John C. Warner  **
**   Licensed under a MIT/X style license              Neil K. Ganju  **
**   See License_ROMS.txt                              Alexis Beudin  **
************************************************* Tarandeep S. Kalra *** 
**                                                                    **
**  Assigns metadata indices for the vegetation module variables that **
**  are used in input and output NetCDF files.  The metadata          **
**  information is read from "varinfo.dat".                           **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
** Submerged aquatic vegetation model variables.
*/
 
#if defined VEG_DRAG || defined VEG_BIOMASS
          CASE ('idvprp(isDens)')
            idvprp(isDens)=varid
          CASE ('idvprp(isDiam)')
            idvprp(isDiam)=varid
          CASE ('idvprp(isHght)')
            idvprp(isHght)=varid
          CASE ('idvprp(isThck)')
            idvprp(isThck)=varid
!#if defined VEG_BIOMASS
!         CASE ('idvprp(pabbm)')
!           idvprp(pabbm)=varid
!         CASE ('idvprp(pbgbm)')
!           idvprp(pbgbm)=varid
!#endif
#endif
#if defined VEG_STREAMING
          CASE ('idWdvg')
            idWdvg=varid
          CASE ('idCdvg')
            idCdvg=varid
#endif
#if defined MARSH_DYNAMICS
          CASE ('idTims')
            idTims=varid
# if defined MARSH_WAVE_THRUST
          CASE ('idTtot')
            idTtot=varid
# endif
# if defined MARSH_RETREAT
          CASE ('idTmmr')
            idTmmr=varid
# endif
# if defined MARSH_TIDAL_RANGE
          CASE('idTmtr')
            idTmtr=varid
# endif
# if defined MARSH_VERT_GROWTH
          CASE('idTmhw')
            idTmhw=varid
          CASE('idTmlw')
            idTmlw=varid
          CASE('idTmbp')
            idTmbp=varid
          CASE('idTmvg')
            idTmvg=varid
          CASE('idTmvt')
            idTmvt=varid
# endif
#endif
