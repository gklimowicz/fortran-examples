# Copyright (c) 2016-2017 Harald Klimach <harald@klimachs.de>
#
# Parts of this file were written by Harald Klimach for University of Siegen.
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
# OR OTHER DEALINGS IN THE SOFTWARE.
# **************************************************************************** #

def set_bench_f_compiler(conf):
	"""Set a specific Fortran compiler to use for compilation.
	   Use this to explicitly set everything for the Fortran compiler.
	"""

	# The Fortran Compiler:
	# ====================

	# Name of the Fortran compiler 'mpif90'
	conf.find_program('mpif90', var='FC')

	# Name of the compiler to use internally for identification.
	conf.env.FC_NAME = 'BENCH'

	# Flags to use for compilation:
	conf.env['FCFLAGS'] = []

	# Flags to use in the linking step:
	conf.env['LINKFLAGS'] = []

	# Which name format does the compiler use for modules?
	# lower = module name in lower case + appended .mod
	# lower.MOD = module name in lower case + appended .MOD
	# UPPER.mod = module name in upper case + appended .mod
	# upper = module name in upper case + appended .MOD
	conf.env.FC_MOD_CAPITALIZATION = 'lower'

	# Version of the compiler to use of the form
	# (major, minor)
	conf.env.FC_VERSION = ('1', '0')

	# Flag to specify source files:
	conf.env['FC_SRC_F']    = []

	# Flag to specify source files in the linking step:
	conf.env['FCLNK_SRC_F'] = []

	# Flag to indicate target objects (NO linking)
	conf.env['FC_TGT_F']    = ['-c', '-o']

	# Flag to indicated target executables (link step)
	conf.env['FCLNK_TGT_F'] = ['-o']

	# Flag to indicate include paths
	conf.env['FCINCPATH_ST']  = '-I%s'

	# Flag to indicate preprocessor defines
	conf.env['FCDEFINES_ST']  = '-D%s'

	# Name pattern to use for static libraries:
	conf.env['fcstlib_PATTERN']   = 'lib%s.a'

	# Flag to provide a library to linkg against:
	conf.env['FCLIB_ST']       = '-l%s'

	# Flag to indicate a path to look for libraries:
	conf.env['FCLIBPATH_ST']   = '-L%s'

	# Corresponding flags for static libraries:
	conf.env['FCSTLIB_ST']     = '-l%s'
	conf.env['FCSTLIBPATH_ST'] = '-L%s'

	conf.add_os_flags('FCFLAGS')

	# The Archiver:
	# ============

	# Find the archiver executable 'ar'.
	conf.find_program('ar', var='AR')

	# Flags for the archiver.
	conf.env.ARFLAGS = ['rcs']

	# ############################################################# #

def set_bench_c_compiler(conf):
	"""Set a specific C compiler to use for compilation.
	   Use this to explicitly set everything for the C compiler
	   and the archiver.
	"""

	# The C Compiler:
	# ==============
	# Name of the C compiler
	conf.find_program('cc', var='CC')

	# Name of the compiler to use internally for identification.
	conf.env.CC_NAME = 'CC'

	# C Compiler flags:
	conf.env['CFLAGS'] = []

	# Flag to specify source files:
	conf.env['CC_SRC_F']            = []

	# Flag to indicate target objects (NO linking)
	conf.env['CC_TGT_F']            = ['-c', '-o']

	# Flag to specify source files in the linking step:
	conf.env['CCLNK_SRC_F']         = []

	# Flag to specify target executables in the linking step:
	conf.env['CCLNK_TGT_F']         = ['-o']

        # Linker
        conf.env['LINK_CC'] = conf.env['CC']

	# Flag to provide include paths
	conf.env['CPPPATH_ST']          = '-I%s'

	# Flag to set defines
	conf.env['DEFINES_ST']          = '-D%s'
                        
	conf.env['LIB_ST']              = '-l%s' # template for adding libs
	conf.env['LIBPATH_ST']          = '-L%s' # template for adding libpaths

	# According templates for static libraries:
	conf.env['STLIB_ST']            = '-l%s'
	conf.env['STLIBPATH_ST']        = '-L%s'

        # Pattern for the program name
	conf.env['cprogram_PATTERN']    = '%s'

        # static lib
	conf.env['LINKFLAGS_cstlib']    = ['-Wl,-Bstatic']
	conf.env['cstlib_PATTERN']      = 'lib%s.a'

	# ############################################################# #


	# The Archiver:
	# ============

	# Find the archiver executable 'ar'.
	conf.find_program('ar', var='AR')

	# Flags for the archiver.
	conf.env.ARFLAGS = ['rcs']

	# ############################################################# #

	# Get flags from the environment, if available:
	conf.add_os_flags('CFLAGS')
	conf.add_os_flags('LINKFLAGS')
