
# these libs are supposed to work on "palm"

ifneq (${MPIDIR},)
FFLAGS += -I${MPIDIR}/include
F90FLAGS += -I${MPIDIR}/include
LIBS += -L${MPIDIR}/lib
endif
LIBS += -lmpi -lmpi++ -lstdc++  -lpthread -lrt -lc

#or maybe :
#
#LIBS += -size_lp64 -mp -L${ESMFLIBDIR} -L${MPIDIR}/lib -lesmf -lmpi \
#-lmpi++  -lcprts -limf -lm -lcxa -lunwind -lrt -ldl -threads \
#${NETCDF_STUBS}
