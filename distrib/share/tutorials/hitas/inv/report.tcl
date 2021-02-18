#!/usr/bin/env avt_shell
 
#############################################################
# Timing Reporting                                          #
#############################################################
avt_config simToolModel hspice

set fig [ttv_LoadSpecifiedTimingFigure inv_chain]
set clist [ttv_GetPaths $fig * * rf 5 critic path max]

set f [fopen inv_chain.paths "w+"]
ttv_DisplayPathListDetail $f $clist
fclose $f
