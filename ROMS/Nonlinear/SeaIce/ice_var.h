/*
** git $Id$
** svn $Id$
*************************************************** Hernan G. Arango ***
** Copyright (c) 2002-2024 The ROMS/TOMS Group                        **
**   Licensed under a MIT/X style license                             **
**   See License_ROMS.md                                              **
************************************************************************
**                                                                    **
**  Assigns metadata indices for the sea ice model (Budgell, 2005;    **
**  Durski and Kurapov, 2019, 2020) variables that are used in input  **
**  and output NetCDF files.  The metadata information is read from   **
**  "varinfo.yaml".                                                   **
**                                                                    **
**  This file is included in file "mod_ncparam.F", routine            **
**  "initialize_ncparm".                                              **
**                                                                    **
************************************************************************
*/

/*
**  Ice state variables.
*/

          CASE ('idAice')
            idAice=varid
            iSice(isAice)=idAice
          CASE ('idHage')
            idHage=varid
            iSice(isHage)=idHage
          CASE ('idHice')
            idHice=varid
            iSice(isHice)=idHice
          CASE ('idHmel')
            idHmel=varid
            iSice(isHmel)=idHmel
          CASE ('idHsno')
            idHsno=varid
            iSice(isHsno)=idHsno
          CASE ('idIage')
            idIage=varid
            iSice(isIage)=idIage
          CASE ('idISxx')
            idISxx=varid
            iSice(isISxx)=idISxx
          CASE ('idISxy')
            idISxy=varid
            iSice(isISxy)=idISxy
          CASE ('idISyy')
            idISyy=varid
            iSice(isISyy)=idISyy
          CASE ('idTice')
            idTice=varid
            iSice(isTice)=idTice
          CASE ('idUice')
            idUice=varid
            iSice(isUice)=idUice
          CASE ('idVice')
            idVice=varid
            iSice(isVice)=idVice
          CASE ('idAiCL')
            idAiCL=varid
          CASE ('idHiCL')
            idHiCL=varid
          CASE ('idIOfv')
            idIOfv=varid
            iFice(icIOfv)=idIOfv
          CASE ('idIOmf')
            idIOmf=varid
            iFice(icIOmf)=idIOmf
          CASE ('idIOmt')
            idIOmt=varid
            iFice(icIOmt)=idIOmt
          CASE ('idIsst')
            idIsst=varid
            iFice(icIsst)=idIsst
          CASE ('idS0mk')
            idS0mk=varid
            iFice(icS0mk)=idS0mk
          CASE ('idT0mk')
            idT0mk=varid
            iFice(icT0mk)=idT0mk
          CASE ('idUiCL')
            idUiCL=varid
          CASE ('idUiER')
            idUiER=varid
          CASE ('idViCL')
            idViCL=varid
          CASE ('idViNR')
            idViNR=varid
          CASE ('idWdiv')
            idWdiv=varid
            iFice(icWdiv)=idWdiv
          CASE ('idW_ai')
            idW_ai=varid
            iFice(icW_ai)=idW_ai
          CASE ('idW_ao')
            idW_ao=varid
            iFice(icW_ao)=idW_ao
          CASE ('idW_fr')
            idW_fr=varid
            iFice(icW_fr)=idW_fr
          CASE ('idW_io')
            idW_io=varid
            iFice(icW_io)=idW_io
          CASE ('idW_ro')
            idW_ro=varid
            iFice(icW_io)=idW_io
/*
**  Ice model open boundary conditions.
*/

          CASE ('isAice(iwest)')
            iceOBC(iwest,isAice)=varid
          CASE ('isAice(ieast)')
            iceOBC(ieast,isAice)=varid
          CASE ('isAice(isouth)')
            iceOBC(isouth,isAice)=varid
          CASE ('isAice(inorth)')
            iceOBC(inorth,isAice)=varid

          CASE ('isHice(iwest)')
            iceOBC(iwest,isHice)=varid
          CASE ('isHice(ieast)')
            iceOBC(ieast,isHice)=varid
          CASE ('isHice(isouth)')
            iceOBC(isouth,isHice)=varid
          CASE ('isHice(inorth)')
            iceOBC(inorth,isHice)=varid

          CASE ('isHmel(iwest)')
            iceOBC(iwest,isHmel)=varid
          CASE ('isHmel(ieast)')
            iceOBC(ieast,isHmel)=varid
          CASE ('isHmel(isouth)')
            iceOBC(isouth,isHmel)=varid
          CASE ('isHmel(inorth)')
            iceOBC(inorth,isHmel)=varid

          CASE ('isHsno(iwest)')
            iceOBC(iwest,isHsno)=varid
          CASE ('isHsno(ieast)')
            iceOBC(ieast,isHsno)=varid
          CASE ('isHsno(isouth)')
            iceOBC(isouth,isHsno)=varid
          CASE ('isHsno(inorth)')
            iceOBC(inorth,isHsno)=varid

          CASE ('isIage(iwest)')
            iceOBC(iwest,isIage)=varid
          CASE ('isIage(ieast)')
            iceOBC(ieast,isIage)=varid
          CASE ('isIage(isouth)')
            iceOBC(isouth,isIage)=varid
          CASE ('isIage(inorth)')
            iceOBC(inorth,isIage)=varid

          CASE ('isISxx(iwest)')
            iceOBC(iwest,isISxx)=varid
          CASE ('isISxx(ieast)')
            iceOBC(ieast,isISxx)=varid
          CASE ('isISxx(isouth)')
            iceOBC(isouth,isISxx)=varid
          CASE ('isISxx(inorth)')
            iceOBC(inorth,isISxx)=varid

          CASE ('isISxy(iwest)')
            iceOBC(iwest,isISxy)=varid
          CASE ('isISxy(ieast)')
            iceOBC(ieast,isISxy)=varid
          CASE ('isISxy(isouth)')
            iceOBC(isouth,isISxy)=varid
          CASE ('isISxy(inorth)')
            iceOBC(inorth,isISxy)=varid

          CASE ('isISyy(iwest)')
            iceOBC(iwest,isISyy)=varid
          CASE ('isISyy(ieast)')
            iceOBC(ieast,isISyy)=varid
          CASE ('isISyy(isouth)')
            iceOBC(isouth,isISyy)=varid
          CASE ('isISyy(inorth)')
            iceOBC(inorth,isISyy)=varid

          CASE ('isTice(iwest)')
            iceOBC(iwest,isTice)=varid
          CASE ('isTice(ieast)')
            iceOBC(ieast,isTice)=varid
          CASE ('isTice(isouth)')
            iceOBC(isouth,isTice)=varid
          CASE ('isTice(inorth)')
            iceOBC(inorth,isTice)=varid

          CASE ('isUice(iwest)')
            iceOBC(iwest,isUice)=varid
          CASE ('isUice(ieast)')
            iceOBC(ieast,isUice)=varid
          CASE ('isUice(isouth)')
            iceOBC(isouth,isUice)=varid
          CASE ('isUice(inorth)')
            iceOBC(inorth,isUice)=varid

          CASE ('isVice(iwest)')
            iceOBC(iwest,isVice)=varid
          CASE ('isVice(ieast)')
            iceOBC(ieast,isVice)=varid
          CASE ('isVice(isouth)')
            iceOBC(isouth,isVice)=varid
          CASE ('isVice(inorth)')
            iceOBC(inorth,isVice)=varid
