# Copyright (c) 2017-2018 Nick Papior <nickpapior@gmail.com>
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

extdouble_DIR = $(TOP_DIR)/source/extdouble

# Append the source/extdouble directory
VPATH := $(VPATH):$(extdouble_DIR)

# Objects in this folder
ifneq ($(EXTDOUBLE),0)
extdouble_OBJECTS = \
	       aot_extdouble_fun_module.o \
	       aot_extdouble_out_module.o \
	       aot_extdouble_table_module.o \
	       aot_extdouble_top_module.o \
	       aot_extdouble_vector_module.o
aot_top_module.o: aot_extdouble_top_module.o
aot_out_module.o: aot_extdouble_out_module.o
aot_fun_module.o: aot_extdouble_fun_module.o
aot_table_module.o: aot_extdouble_table_module.o
aot_vector_module.o: aot_extdouble_vector_module.o
aot_extdouble_top_module.o: aot_err_module.o
else
extdouble_OBJECTS = \
	       dummy_extdouble_fun_module.o \
	       dummy_extdouble_out_module.o \
	       dummy_extdouble_table_module.o \
	       dummy_extdouble_top_module.o \
	       dummy_extdouble_vector_module.o
aot_top_module.o: dummy_extdouble_top_module.o
aot_out_module.o: dummy_extdouble_out_module.o
aot_fun_module.o: dummy_extdouble_fun_module.o
aot_table_module.o: dummy_extdouble_table_module.o
aot_vector_module.o: dummy_extdouble_vector_module.o
endif
extdouble_SOURCES := $(patsubst %.o, %.f90, $(extdouble_OBJECTS))

# Add dependency for the LuaFortran objects
$(extdouble_OBJECTS): $(LuaF_OBJECTS)

# Add to the global variables
SOURCES += $(extdouble_SOURCES)
OBJECTS += $(extdouble_OBJECTS)


# Local Variables:
#  mode: makefile-gmake
# End:

# Automatic creation of the dependencies may be achieved through:
#   sfmakedepend --depend=obj --modext=o --file=Makefile.inc *.f90

# DO NOT DELETE THIS LINE - used by make depend
aot_extdouble_fun_module.o: aot_extdouble_top_module.o
aot_extdouble_out_module.o: aot_extdouble_top_module.o
aot_extdouble_table_module.o: aot_extdouble_top_module.o
aot_extdouble_vector_module.o: aot_extdouble_top_module.o
