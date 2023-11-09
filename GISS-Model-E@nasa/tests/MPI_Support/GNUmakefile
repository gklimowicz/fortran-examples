.PHONY: clean all

ifeq ($(USE_MPI),YES)
SRCS = Test_DistGrid.pf \
	Test_GatherScatter.pf \
	Test_GlobalSum.pf \
	Test_Halo.pf
else
SRCS = Test_MPIPlaceholder.pf
endif

OBJS = $(SRCS:.pf=.o)

LIB = libMPI_Support.a

$(LIB): $(OBJS)
	ar -r $@ $?
	-@ranlib $@

all: $(OBJS) $(LIB)

testSuites.inc: $(SRCS)

%.F90: %.pf
	python $(PFUNIT)/bin/pFUnitParser.py $<  $@

%.o: %.F90
	$(FC) -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a

