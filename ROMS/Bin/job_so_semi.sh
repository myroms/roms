#!/bin/csh -f
#
# git $Id: 560fb376ff8a4576170ebcd4b459de6bcce908f6 $
# svn $Id: job_so_semi.sh 937 2019-01-28 06:13:04Z arango $
#######################################################################
# Copyright (c) 2002-2019 The ROMS/TOMS Group                         #
#   Licensed under a MIT/X style license                              #
#   See License_ROMS.txt                                              #
#######################################################################
#                                                                     #
#  Generalized Stability Theory: Stochastic Optimals, seminorm        #
#                                                                     #
#  This script is used to run the ROMS/TOMS stochastic optimals with  #
#  respect to the seminorm of the chosen functional.                  #
#                                                                     #
#######################################################################

# Set ROOT of the directory to run application.  The following
# "dirname" command returns a path by removing any suffix from
# the last slash ('/').  It returns a path above current diretory.

set Dir=`dirname ${PWD}`

# Set basic state trajectory, forward file:

#set HISname=${Dir}/Forward/gyre3d_his_00.nc
 set HISname=${Dir}/Forward/gyre3d_his_01.nc

set FWDname=gyre3d_fwd.nc

if (-e $FWDname) then
  /bin/rm $FWDname
endif
ln -s -v $HISname $FWDname

# Set adjoint model initial conditions file: zero fields.

set IADname=gyre3d_iad.nc

if (-e $IADname) then
  /bin/rm $IADname
endif

ln -s -v ${Dir}/Data/gyre3d_ini_zero.nc $IADname
