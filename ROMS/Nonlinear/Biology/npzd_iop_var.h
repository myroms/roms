/*
** git $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2024 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the IOP-based, NPZD ecosystem model  **
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

              CASE ('idTvar(iDIN_)')
                idTvar(iDIN_)=varid
              CASE ('idTvar(iAphy(i440n))')
                idTvar(iAphy(i440n))=varid
              CASE ('idTvar(iAphy(i510n))')
                idTvar(iAphy(i510n))=varid
              CASE ('idTvar(iBphy(i440n))')
                idTvar(iBphy(i440n))=varid
              CASE ('idTvar(iBphy(i510n))')
                idTvar(iBphy(i510n))=varid
              CASE ('idTvar(aCDOM(i440n))')
                idTvar(aCDOM(i440n))=varid
              CASE ('idTvar(aCDOM(i510n))')
                idTvar(aCDOM(i510n))=varid
              CASE ('idTvar(iBdet(i440n))')
                idTvar(iBdet(i440n))=varid
              CASE ('idTvar(iBdet(i510n))')
                idTvar(iBdet(i510n))=varid

#if defined AD_SENSITIVITY   || defined IS4DVAR_SENSITIVITY || \
    defined OPT_OBSERVATIONS || defined SENSITIVITY_4DVAR   || \
    defined SO_SEMI

/*
**  Adjoint sensitivity state biological tracers.
*/

              CASE ('idTads(iDIN_)')
                idTads(iDIN_)=varid
              CASE ('idTads(iAphy(i440n))')
                idTads(iAphy(i440n))=varid
              CASE ('idTads(iAphy(i510n))')
                idTads(iAphy(i510n))=varid
              CASE ('idTads(iBphy(i440n))')
                idTads(iBphy(i440n))=varid
              CASE ('idTads(iBphy(i510n))')
                idTads(iBphy(i510n))=varid
              CASE ('idTads(aCDOM(i440n))')
                idTads(aCDOM(i440n))=varid
              CASE ('idTads(aCDOM(i510n))')
                idTads(aCDOM(i510n))=varid
              CASE ('idTads(iBdet(i440n))')
                idTads(iBdet(i440n))=varid
              CASE ('idTads(iBdet(i510n))')
                idTads(iBdet(i510n))=varid
#endif

/*
**  Biological tracers open boundary conditions.
*/

              CASE ('idTbry(iwest,iDIN_)')
                idTbry(iwest,iDIN_)=varid
              CASE ('idTbry(ieast,iDIN_)')
                idTbry(ieast,iDIN_)=varid
              CASE ('idTbry(isouth,iDIN_)')
                idTbry(isouth,iDIN_)=varid
              CASE ('idTbry(inorth,iDIN_)')
                idTbry(inorth,iDIN_)=varid

              CASE ('idTbry(iwest,iAphy(i440n))')
                idTbry(iwest,iAphy(i440n))=varid
              CASE ('idTbry(ieast,iAphy(i440n))')
                idTbry(ieast,iAphy(i440n))=varid
              CASE ('idTbry(isouth,iAphy(i440n))')
                idTbry(isouth,iAphy(i440n))=varid
              CASE ('idTbry(inorth,iAphy(i440n))')
                idTbry(inorth,iAphy(i440n))=varid

              CASE ('idTbry(iwest,iAphy(i510n))')
                idTbry(iwest,iAphy(i510n))=varid
              CASE ('idTbry(ieast,iAphy(i510n))')
                idTbry(ieast,iAphy(i510n))=varid
              CASE ('idTbry(isouth,iAphy(i510n))')
                idTbry(isouth,iAphy(i510n))=varid
              CASE ('idTbry(inorth,iAphy(i510n))')
                idTbry(inorth,iAphy(i510n))=varid

              CASE ('idTbry(iwest,iBphy(i440n))')
                idTbry(iwest,iBphy(i440n))=varid
              CASE ('idTbry(ieast,iBphy(i440n))')
                idTbry(ieast,iBphy(i440n))=varid
              CASE ('idTbry(isouth,iBphy(i440n))')
                idTbry(isouth,iBphy(i440n))=varid
              CASE ('idTbry(inorth,iBphy(i440n))')
                idTbry(inorth,iBphy(i440n))=varid

              CASE ('idTbry(iwest,iBphy(i510n))')
                idTbry(iwest,iBphy(i510n))=varid
              CASE ('idTbry(ieast,iBphy(i510n))')
                idTbry(ieast,iBphy(i510n))=varid
              CASE ('idTbry(isouth,iBphy(i510n))')
                idTbry(isouth,iBphy(i510n))=varid
              CASE ('idTbry(inorth,iBphy(i510n))')
                idTbry(inorth,iBphy(i510n))=varid

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

              CASE ('idTbry(iwest,iBdet(i440n))')
                idTbry(iwest,iBdet(i440n))=varid
              CASE ('idTbry(ieast,iBdet(i440n))')
                idTbry(ieast,iBdet(i440n))=varid
              CASE ('idTbry(isouth,iBdet(i440n))')
                idTbry(isouth,iBdet(i440n))=varid
              CASE ('idTbry(inorth,iBdet(i440n))')
                idTbry(inorth,iBdet(i440n))=varid

              CASE ('idTbry(iwest,iBdet(i510n))')
                idTbry(iwest,iBdet(i510n))=varid
              CASE ('idTbry(ieast,iBdet(i510n))')
                idTbry(ieast,iBdet(i510n))=varid
              CASE ('idTbry(isouth,iBdet(i510n))')
                idTbry(isouth,iBdet(i510n))=varid
              CASE ('idTbry(inorth,iBdet(i510n))')
                idTbry(inorth,iBdet(i510n))=varid

#ifdef TS_PSOURCE

/*
**  Biological tracers point Source/Sinks (river runoff).
*/

              CASE ('idRtrc(iDIN_)')
                idRtrc(iDIN_)=varid
              CASE ('idRtrc(iAphy(i440n))')
                idRtrc(iAphy(i440n))=varid
              CASE ('idRtrc(iAphy(i510n))')
                idRtrc(iAphy(i510n))=varid
              CASE ('idRtrc(iBphy(i440n))')
                idRtrc(iBphy(i440n))=varid
              CASE ('idRtrc(iBphy(i510n))')
                idRtrc(iBphy(i510n))=varid
              CASE ('idRtrc(aCDOM(i440n))')
                idRtrc(aCDOM(i440n))=varid
              CASE ('idRtrc(aCDOM(i510n))')
                idRtrc(aCDOM(i510n))=varid
              CASE ('idRtrc(iBdet(i440n))')
                idRtrc(iBdet(i440n))=varid
              CASE ('idRtrc(iBdet(i510n))')
                idRtrc(iBdet(i510n))=varid
#endif
