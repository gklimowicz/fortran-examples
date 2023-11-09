
F90 = g95
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS += -DCOMPILER_G95
FFLAGS = -cpp -fno-second-underscore -O2 -fendian=big
F90FLAGS = -cpp -fno-second-underscore -O2 -ffree-form -fendian=big
LFLAGS =
# uncomment next two lines for extensive debugging
# the following switch adds extra debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -fbounds-check -ftrace=full -freal=nan -fpointer=invalid
LFLAGS += #-lefence
endif
