#!/bin/bash

# set the number of OpenMP threads per node equal to the number of cores
# (this environment variable does not normally need to be set)
#export OMP_NUM_THREADS=

# no binding of threads to cores
export OMP_PROC_BIND=false

# increase the OpenMP stack size
export OMP_STACKSIZE=256M

# set the soft limit of the stack size equal to the hard limit
ulimit -Ss unlimited

# Elk executable file
~/elk/src/elk
