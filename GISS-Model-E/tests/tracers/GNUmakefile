.PHONY: clean tests

SRCS = Test_Tracer.pf 
SRCS += Test_TracerBundle.pf
SRCS += Test_TracerBundleSubset.pf
SRCS += Test_TracerIO.pf
SRCS += Test_TracerReference.pf
SRCS += Test_TracerSurfaceSource.pf
OBJS = $(SRCS:.pf=.o)

LIB = libtracers.a

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

