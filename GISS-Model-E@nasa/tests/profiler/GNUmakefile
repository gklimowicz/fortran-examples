.PHONY: clean all

ALLFILES:=$(wildcard *.pf)
# We need to exclude MPI test files when running serial tests to 
# avoid compilation errors.
MPIFILES:=Test_TimerParallel.pf

ifeq ($(USE_MPI),YES)
  SRCS := $(ALLFILES)
else
  SRCS := $(filter-out $(MPIFILES),$(ALLFILES))
endif
OBJS := $(SRCS:.pf=.o)

LIB = libprofiler.a

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

