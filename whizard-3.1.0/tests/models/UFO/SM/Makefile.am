## Makefile.am -- Makefile for models used in WHIZARD test scripts
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
## Definitions

# The UFO model is contained in the OMega test suite.  Copy all files here.
OMEGA_TEST_UFO_DIR = $(top_srcdir)/omega/tests/UFO
MODEL = SM

UFO_FILES = \
  __init__.py \
  coupling_orders.py \
  couplings.py \
  CT_couplings.py \
  decays.py \
  function_library.py \
  lorentz.py \
  object_library.py \
  parameters.py \
  particles.py \
  propagators.py \
  vertices.py \
  write_param_card.py

$(UFO_FILES): $(OMEGA_TEST_UFO_DIR)/$(MODEL)/__init__.py
	cp $(OMEGA_TEST_UFO_DIR)/$(MODEL)/$@ .

########################################################################
## Targets

all: $(UFO_FILES)

distclean-local:
	rm -f $(UFO_FILES)

########################################################################
## The End.
########################################################################
