
F90 = pgf90
FMAKEDEP = $(SCRIPTS_DIR)/sfmakedepend
CPPFLAGS += -DCOMPILER_PGI
FFLAGS = -O2  -Mpreprocess -Mbyteswapio
F90FLAGS = -O2  -Mpreprocess -Mbyteswapio
LFLAGS = 
ifeq ($(MP),YES)   # edit for PGI OpenMP compatability???
FFLAGS += -Msmp
LFLAGS += -Msmp
endif
# uncomment next two lines for extensive debugging
ifeq ($(COMPILE_WITH_TRAPS),YES)
FFLAGS += -Mbounds -Mchkfpstk -Mchkptr -Mchkstk
LFLAGS += -Mbounds -Mchkfpstk -Mchkptr -Mchkstk
endif
