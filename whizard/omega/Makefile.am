# Makefile.am -- Makefile for O'Mega within and without WHIZARD
##
## Process this file with automake to produce Makefile.in
##
########################################################################
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

## we do our own dependency checking for FORTRAN and OCaml:
AUTOMAKE_OPTIONS = no-dependencies

## Use the m4 directory for local Autoconf macros
ACLOCAL_AMFLAGS = -I m4

## Subdirectories to configure
if OCAML_AVAILABLE 
OPT_OMEGA = bin models share tests extensions
else
OPT_OMEGA = 
endif
SUBDIRS = src lib $(OPT_OMEGA) tools scripts
DIST_SUBDIRS = src lib models bin share tools extensions tests scripts

########################################################################

clean-local:

########################################################################
## The End.
########################################################################

