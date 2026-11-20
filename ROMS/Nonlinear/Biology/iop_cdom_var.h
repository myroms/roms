/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2024 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the IOP-based, CDOM ecosystem model  **
**  variables that are used in input and output NetCDF files.  The    **
**  metadata information is read from "varinfo.dat".                  **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Model state biological tracers.
*/

              CASE ('idTvar(aCDOM(i440n))')
                idTvar(aCDOM(i440n))=varid
              CASE ('idTvar(aCDOM(i510n))')
                idTvar(aCDOM(i510n))=varid

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI

/*
**  Adjoint sensitivity state biological tracers.
*/

              CASE ('idTads(aCDOM(i440n))')
                idTads(aCDOM(i440n))=varid
              CASE ('idTads(aCDOM(i510n))')
                idTads(aCDOM(i510n))=varid
#endif

/*
**  Biological tracers open boundary conditions.
*/

              CASE ('idTbry(iwest,aCDOM(i440n))')
                idTbry(iwest,aCDOM(i440n))=varid
              CASE ('idTbry(ieast,aCDOM(i440n))')
                idTbry(ieast,aCDOM(i440n))=varid
              CASE ('idTbry(isouth,aCDOM(i440n))')
                idTbry(isouth,aCDOM(i440n))=varid
              CASE ('idTbry(inorth,aCDOM(i440n))')
                idTbry(inorth,aCDOM(i440n))=varid

              CASE ('idTbry(iwest,aCDOM(i510n))')
                idTbry(iwest,aCDOM(i510n))=varid
              CASE ('idTbry(ieast,aCDOM(i510n))')
                idTbry(ieast,aCDOM(i510n))=varid
              CASE ('idTbry(isouth,aCDOM(i510n))')
                idTbry(isouth,aCDOM(i510n))=varid
              CASE ('idTbry(inorth,aCDOM(i510n))')
                idTbry(inorth,aCDOM(i510n))=varid

#ifdef TS_PSOURCE

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

              CASE ('idRtrc(aCDOM(i440n))')
                idRtrc(aCDOM(i440n))=varid
              CASE ('idRtrc(aCDOM(i510n))')
                idRtrc(aCDOM(i510n))=varid
#endif
