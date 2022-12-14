#include "rundeck_opts.h"

      MODULE DOMAIN_DECOMP_ATM
        use dist_grid_mod
        use Halo_mod
        use SpecialIO_mod
        use GatherScatter_mod
        use GlobalSum_mod
        use iso_c_binding
      implicit none
      public
      type(dist_grid), target :: grid
#ifdef GLINT2
      type(c_ptr) :: glint2
#endif
      END MODULE DOMAIN_DECOMP_ATM
