## Makefile.am -- Makefile for WHIZARD
##
#
# Copyright (C) 1999-2022 by  
#     Wolfgang Kilian <kilian@physik.uni-siegen.de> 
#     Thorsten Ohl <ohl@physik.uni-wuerzburg.de> 
#     Juergen Reuter <juergen.reuter@desy.de> 
#     with contributions from
#     cf. main AUTHORS file
# 
# WHIZARD is free software; you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by  
# the Free Software Foundation; either version 2, or (at your option) 
# any later version. 
# 
# WHIZARD is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License 
# along with this program; if not, write to the Free Software 
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. 
# 
######################################################################## 
## Process this file with automake to produce Makefile.in

## Distribute the style file and color definition files, and install them here:
feynmptexdir = $(PKGTEXDIR)
dist_feynmptex_DATA = \
  $(FEYNMP_FILES)

## LaTeX style file
FEYNMP_FILES = \
  feynmp.sty feynmp.mp

## Run texhash to include the style file in the TeX search path
install-data-hook:
	test -x `which texhash` && texhash

## Remove backup files
maintainer-clean-local: maintainer-clean-noweb
	-rm -f *~ 
