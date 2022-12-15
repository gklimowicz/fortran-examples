# Copyright (c) 2017 Nick Papior <nickpapior@gmail.com>
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

LuaF_DIR = $(TOP_DIR)/LuaFortran

# Append the LuaFortran directory
VPATH := $(VPATH):$(LuaF_DIR)

# Objects in this folder
LuaF_OBJECTS = dump_lua_fif_module.o \
	     flu_binding.o \
	     flu_kinds_module.o \
	     lua_fif.o \
	     lua_parameters.o
LuaF_SOURCES := $(patsubst %.o, %.f90, $(LuaF_OBJECTS))

# Add the C-sources
LuaF_OBJECTS += wrap_lua_dump.o
LuaF_SOURCES += wrap_lua_dump.c


# Add to the global variables
SOURCES += $(LuaF_SOURCES)
OBJECTS += $(LuaF_OBJECTS)


# Local Variables:
#  mode: makefile-gmake
# End:

# Automatic creation of the dependencies may be achieved through:
#   sfmakedepend --depend=obj --modext=o --file=Makefile.inc *.f90

# DO NOT DELETE THIS LINE - used by make depend
flu_binding.o: dump_lua_fif_module.o flu_kinds_module.o lua_fif.o
flu_binding.o: lua_parameters.o
lua_fif.o: lua_parameters.o
