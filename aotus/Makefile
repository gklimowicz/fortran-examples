
# Define variables for internal Lua source code
EXTERNAL_LUA_V = 5.3.4
# Default platform
PLATFORM ?= linux

# Define VPATH
VPATH ?= $(shell pwd)

# the default target
.PHONY: default
default: lib

# the all target
.PHONY: all static
all: static


# SMEKASETTINGS (DO NOT DELETE)
# DO NOT CHANGE CONTENT IN THIS BLOCK
# IT MAY BE OVERWRITTEN WHEN REINSTALLING SMEKA
#
# This Makefile was created by smeka:
#  github.com/zerothi/smeka

# Top-directory of Makefile/source tree
# If need set, do so ABOVE this block!
TOP_DIR ?= .

# Directory of smeka default Makefiles
SMEKA_DIR = smeka

# Include the smeka settings!
include $(TOP_DIR)/$(SMEKA_DIR)/Makefile.smeka

# SMEKAENDSETTINGS (DO NOT DELETE)


# The linker is a fortran compiler
LINK := $(FC)

# Create targets for the library
.PHONY: lib lua-lib
lib: $(LIBRARIES)

# Determine whether the user has added an external
# Lua header file.
ifdef LUA_DIR
 # User-defined
 INCLUDES += -I$(LUA_DIR)/include
 LUA__LIB = $(LUA_DIR)/lib/liblua.a
lua-lib:
	@echo ""
	@echo "Adding external Lua library to lib$(AOTUS_LIB).*"
	@echo "Extracting $(LUA_DIR)/lib/liblua.a"
	@echo ""
	$(AR) -x $(LUA__LIB)

else

 # lua-sources shipped
 LUA_DIR = external/lua-$(EXTERNAL_LUA_V)
 INCLUDES += -I$(LUA_DIR)/src
 LUA__LIB = $(LUA_DIR)/src/liblua.a
lua-lib:
	mkdir -p external
	cp -rf $(TOP_DIR)/$(LUA_DIR) external/
	$(MAKE) -C $(LUA_DIR) $(PLATFORM)
	$(AR) -x $(LUA__LIB)

# Ensure the Lua-library is built first
$(LIBRARIES): |lua-lib

# Add cleaning target
.PHONY: clean-lua-lib
clean-lua-lib:
	$(MAKE) -C $(LUA_DIR) clean
	-rm -f *.o

clean: clean-lua-lib

endif


# Add the LuaFortran files
include $(TOP_DIR)/LuaFortran/Makefile.inc

# Add source files
include $(TOP_DIR)/source/Makefile.inc
# Add source/extdouble files
include $(TOP_DIR)/source/extdouble/Makefile.inc
# Add source/quadruple files
include $(TOP_DIR)/source/quadruple/Makefile.inc


# Ensure that all objects are required for the libraries
$(LIBRARIES): $(OBJECTS)
