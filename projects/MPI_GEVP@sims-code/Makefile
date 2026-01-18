#! /usr/bin/make
#
LD        = $(FC)
FC       = mpif90 
SWP       = swplist
RM        = /bin/rm -f
MP        =
ABI       = 
ISA       = 
PROC      = 
ARCH      = $(MP) $(ABI) $(ISA)
OLEVEL    = -O
FOPTS  = 
# For better precision, debugging, and to check unitialized variables, 
# use this also
#FOPTS  =  -g -check bounds -debug -traceback -mieee-fp -check uninit
COPTS     = 
FFLAGS    = $(ARCH) $(OLEVEL) $(FOPTS) 
CFLAGS    = $(ARCH) $(OLEVEL) $(COPTS)
LDFLAGS   = $(ARCH) $(OLEVEL) 
PROF      =

FOBJS  = MPI_GEVP_inviter_qp.o MPI_GEVP_inviter_solve_qp.o

COBJS  = 
OBJS      = $(FOBJS) $(COBJS)

EXEC      = MPI_GEVP_inviter_solve_qp.exe

LIB       = 

all: $(EXEC)

$(EXEC):  $(OBJS) 
	$(LD) -o $@ $(OBJS)

clean:
	$(RM) $(OBJS)

.SUFFIXES:
.SUFFIXES: .o .c .f .for .f0 


.f.o:
	$(FC) -c $(FFLAGS) $<

.for.o:
	$(FC)  -c $(FFLAGS) $<

.c.o:
	$(CC)  -c $(CFLAGS) $(DEFINES) $<
