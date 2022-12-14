.PHONY: clean all

SRCS = $(wildcard *.pf)
OBJS = $(SRCS:.pf=.o FooAssociativeArray.o FooHashMap.o)

LIB = libshared.a

$(LIB): $(OBJS)
	ar -r $@ $?
	-@ranlib $@

# This is needed to find the templates
FFLAGS += -I../../model/include -I../../model/mod

all: FooAssociativeArray.o FooHashMap.o $(OBJS) $(LIB)

testSuites.inc: $(SRCS)

%.F90: %.pf
	python $(PFUNIT)/bin/pFUnitParser.py $<  $@

%.o: %.F90
	$(FC) -c $(FFLAGS) $<

clean:
	$(RM) *.o *.mod *.a

