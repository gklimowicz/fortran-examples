## Makefile.am -- Makefile for WHIZARD
##
## Process this file with automake to produce Makefile.in
#
# Copyright (C) 1999-2020 by
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

########################################################################
## Python Module

if PYTHON_API_AVAILABLE
all-local: build
build:
	$(PYTHON) setup.py build_ext

install-exec-local:
	$(PYTHON) setup.py install \
		--prefix $(DESTDIR)$(prefix) \
		--single-version-externally-managed \
		--record $(DESTDIR)$(PKGLIBDIR)/install_files.txt
uninstall-local:
	xargs rm -rf < $(DESTDIR)$(PKGLIBDIR)/install_files.txt
endif

########################################################################

VPATH = $(srcdir)

EXTRA_DIST = \
  whizard_python.pyx \
  cwhizard.pxd

########################################################################
## Non-standard cleanup tasks

if PYTHON_API_AVAILABLE
clean-local:
	$(PYTHON) setup.py clean \
		--all
	rm -f whizard_python.cpp
	rm -rf PyWHIZARD.egg-info
else
clean-local:
endif

## Remove backup files
maintainer-clean-backup:
	-rm -f *~
.PHONY: maintainer-clean-backup

## Register additional clean targets
maintainer-clean-local: maintainer-clean-backup
