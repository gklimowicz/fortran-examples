!@sum this file contains macros which are supposed to simplify
!@+   parallelization of HYCOM code.

#ifdef USE_MPI
#define PERIODIC_INDEX(x,y) (x)
#else
#define PERIODIC_INDEX(x,y) (mod(x-1+y,y)+1)
#endif
