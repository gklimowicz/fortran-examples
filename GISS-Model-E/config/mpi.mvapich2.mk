ifneq (${MPIDIR},)
ifneq ($(wildcard $(MPIDIR)/include/mpi.h),)
  FFLAGS += -I${MPIDIR}/include
  F90FLAGS += -I${MPIDIR}/include
  CPPFLAGS += -I${MPIDIR}/include
else
  $(error MPI distribution not found. Check settings in ~/.modelErc)
endif
LIBS += -L${MPIDIR}/lib
endif

LIBS += -lm -lrt -ldl -lmpich -lfmpich -lstdc++ -lpthread -lrdmacm -libverbs -libumad -lmpl -lnuma

