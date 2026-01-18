# C compiler
#CC = mpicc -O3 -openmp -I/opt/fftw/3.3/include -I./include

# Library path
#LIB = -L/opt/fftw/3.3/lib -lfftw3

# Fortran compiler for Fortran users to compile the Fortran sample
# If (FC = 0) is specified, the Fortran sample will not be compiled
#FC = 0

# Samples of CC and LIB
# Cray XC30 (module load fftw first)
#CC = cc -O3 -openmp -xHOST -xAVX -I./include
#FC = ftn -O3 -openmp -xHOST -xAVX -I./include
#LIB = -L/opt/fftw/3.3.0.4/x86_64/lib -lfftw3

# K computer
#CC = mpifccpx -Kopenmp -Kfast -Dkcomp -I/home/apps/fftw/3.3.3/include -I./include
#LIB = -L/home/apps/fftw/3.3.3/lib64 -lfftw3 
#FC = mpifrtpx -Kopenmp -Kfast -Dkcomp -I/home/apps/fftw/3.3.3/include -I./include

# Fujitsu FX10
#CC = mpifccpx -Kfast,SPARC64IXfx -Kopenmp -Dkcomp-I/usr/local/fftw/3.3/include -I./include
#LIB = -L/usr/local/fftw/3.3/lib64 -lfftw3
#FC = mpifrtpx -Kfast,SPARC64IXfx -Kopenmp -Dkcomp -I/usr/local/fftw/3.3/include -I./include

# SGI Altix UV1000
#CC = icc -O3 -openmp -xHOST -I/opt/intel/mkl/include/fftw -I/opt/sgi/mpt/mpt-2.05/include/ -I./include
#LIB= -L/opt/sgi/mpt/mpt-2.05/lib/ -L/opt/intel/mkl/lib -mkl=parallel -lmpi
#FC = ifort -O3 -openmp -xHOST -I/opt/intel/mkl/include/fftw -I/opt/sgi/mpt/mpt-2.05/include/ -I./include

# SGI Infiniband machine
CC = icc -O3 -openmp -xHOST -xAVX -I/work/duytvt/fftw/3.3.4/include -I./include 
LIB= -L/work/duytvt/fftw-3.3.4/lib -lfftw3 -mkl -lmpi
FC = ifort -O3 -openmp -xHOST -xAVX -I/work/duytvt/fftw/3.3.4/include -I./include

# Source files

SRCDIR = ./source

SRC =  	$(SRCDIR)/openfft_init_c2c_3d.c $(SRCDIR)/openfft_init_r2c_3d.c $(SRCDIR)/openfft_init_c2c_4d.c $(SRCDIR)/openfft_exec_c2c_3d.c $(SRCDIR)/openfft_exec_r2c_3d.c $(SRCDIR)/openfft_exec_c2c_4d.c $(SRCDIR)/openfft_finalize.c $(SRCDIR)/openfft_dtime.c

OBJ = $(SRC:.c=.o)


# OpenFFT library location

OFFTLIB = ./lib/libopenfft.a


all: lib samples

lib: 	libopenfft.a 

libopenfft.a:  $(OFFTLIB)

%.o: %.c
	$(CC) -c $< -o $@

$(OFFTLIB): $(OBJ)
	ar rcs $(OFFTLIB) $(OBJ)


# C sample programs

EXPDIR = ./samples
CEXPDIR = $(EXPDIR)/C
FEXPDIR = $(EXPDIR)/FORTRAN

samples: csamples fsamples

csamples: 

	$(CC) -o $(CEXPDIR)/check_c2c_3d $(CEXPDIR)/check_c2c_3d.c $(OFFTLIB) $(LIB)
	$(CC) -o $(CEXPDIR)/check_r2c_3d $(CEXPDIR)/check_r2c_3d.c $(OFFTLIB) $(LIB)
	$(CC) -o $(CEXPDIR)/check_c2c_4d $(CEXPDIR)/check_c2c_4d.c $(OFFTLIB) $(LIB)
	$(CC) -o $(CEXPDIR)/timing_c2c_3d $(CEXPDIR)/timing_c2c_3d.c $(OFFTLIB) $(LIB)
	$(CC) -o $(CEXPDIR)/timing_r2c_3d $(CEXPDIR)/timing_r2c_3d.c $(OFFTLIB) $(LIB)
	$(CC) -o $(CEXPDIR)/timing_c2c_4d $(CEXPDIR)/timing_c2c_4d.c $(OFFTLIB) $(LIB)
	$(CC) -o $(CEXPDIR)/breaktime_c2c_3d $(CEXPDIR)/breaktime_c2c_3d.c $(OFFTLIB) $(LIB)
	$(CC) -o $(CEXPDIR)/breaktime_r2c_3d $(CEXPDIR)/breaktime_r2c_3d.c $(OFFTLIB) $(LIB)
	$(CC) -o $(CEXPDIR)/breaktime_c2c_4d $(CEXPDIR)/breaktime_c2c_4d.c $(OFFTLIB) $(LIB)

fsamples:

ifneq ($(FC),0)
	$(FC) -o $(FEXPDIR)/check_c2c_3d $(FEXPDIR)/check_c2c_3d.f90 $(OFFTLIB) $(LIB)
	$(FC) -o $(FEXPDIR)/check_r2c_3d $(FEXPDIR)/check_r2c_3d.f90 $(OFFTLIB) $(LIB)
	$(FC) -o $(FEXPDIR)/timing_c2c_3d $(FEXPDIR)/timing_c2c_3d.f90 $(OFFTLIB) $(LIB)
	$(FC) -o $(FEXPDIR)/timing_r2c_3d $(FEXPDIR)/timing_r2c_3d.f90 $(OFFTLIB) $(LIB)
endif


clean:
	rm -f $(OBJ) $(OFFTLIB) $(CEXPDIR)/check_c2c_3d $(CEXPDIR)/check_r2c_3d $(CEXPDIR)/check_c2c_4d $(CEXPDIR)/timing_c2c_3d $(CEXPDIR)/timing_r2c_3d $(CEXPDIR)/timing_c2c_4d $(CEXPDIR)/breaktime_c2c_3d $(CEXPDIR)/breaktime_r2c_3d $(CEXPDIR)/breaktime_c2c_4d $(FEXPDIR)/check_c2c_3d $(FEXPDIR)/check_r2c_3d $(FEXPDIR)/timing_c2c_3d  $(FEXPDIR)/timing_r2c_3d $(CEXPDIR)/*~ $(FEXPDIR)/*~ $(SRCDIR)/*~ ./include/*~ Makefile.bak *~

