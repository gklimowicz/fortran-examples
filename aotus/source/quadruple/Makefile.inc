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

quadruple_DIR = $(TOP_DIR)/source/quadruple

# Append the source/quadruple directory
VPATH := $(VPATH):$(quadruple_DIR)

# Objects in this folder
ifneq ($(QUADRUPLE),0)
quadruple_OBJECTS = \
	       aot_quadruple_fun_module.o \
	       aot_quadruple_out_module.o \
	       aot_quadruple_table_module.o \
	       aot_quadruple_top_module.o \
	       aot_quadruple_vector_module.o
aot_top_module.o: aot_quadruple_top_module.o
aot_out_module.o: aot_quadruple_out_module.o
aot_fun_module.o: aot_quadruple_fun_module.o
aot_table_module.o: aot_quadruple_table_module.o
aot_vector_module.o: aot_quadruple_vector_module.o
aot_quadruple_top_module.o: aot_err_module.o
else
quadruple_OBJECTS = \
	       dummy_quadruple_fun_module.o \
	       dummy_quadruple_out_module.o \
	       dummy_quadruple_table_module.o \
	       dummy_quadruple_top_module.o \
	       dummy_quadruple_vector_module.o
aot_top_module.o: dummy_quadruple_top_module.o
aot_out_module.o: dummy_quadruple_out_module.o
aot_fun_module.o: dummy_quadruple_fun_module.o
aot_table_module.o: dummy_quadruple_table_module.o
aot_vector_module.o: dummy_quadruple_vector_module.o
endif
quadruple_SOURCES := $(patsubst %.o, %.f90, $(quadruple_OBJECTS))

# Add dependency for the LuaFortran objects
$(quadruple_OBJECTS): $(LuaF_OBJECTS)

# Add to the global variables
SOURCES += $(quadruple_SOURCES)
OBJECTS += $(quadruple_OBJECTS)


# Local Variables:
#  mode: makefile-gmake
# End:


# Automatic creation of the dependencies may be achieved through:
#   sfmakedepend --depend=obj --modext=o --file=Makefile.inc aot*.f90

# DO NOT DELETE THIS LINE - used by make depend
aot_quadruple_fun_module.o: aot_quadruple_top_module.o
aot_quadruple_out_module.o: aot_quadruple_top_module.o
aot_quadruple_table_module.o: aot_quadruple_top_module.o
aot_quadruple_vector_module.o: aot_quadruple_top_module.o
