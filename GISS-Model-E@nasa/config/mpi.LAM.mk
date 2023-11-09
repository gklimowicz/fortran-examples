
# hack for LAM that comes in Linux distrubution
# (compiled with g95 ? )
CPPFLAGS += -DMPILIB_DOUBLE_UNDERSCORE -DMPI_DEFS_HACK

ifneq (${MPIDIR},)
LIBS += -L${MPIDIR}/lib
endif
LIBS += -llammpi++ -llamf77mpi -lmpi \
-llam -lpthread -lrt -lc /usr/lib64/libstdc++.so.6

