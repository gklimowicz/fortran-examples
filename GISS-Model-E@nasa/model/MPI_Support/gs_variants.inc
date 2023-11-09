#undef _DOMAIN_
#ifdef MULTITILE
#define _DOMAIN_ multitile
#endif
#ifdef TILE
#define _DOMAIN_ 1tile
#endif
#ifdef ROW
#define _DOMAIN_ row
#endif


c
c
c
#undef _ARRDIMSL_
#define _ARRDIMSL_ :,:
      subroutine do_2D
#include "gs_setup.inc"
      return
      end subroutine
c
c
c
#undef _ARRDIMSL_
#define _ARRDIMSL_ :,:,:
      subroutine do_3D
#include "gs_setup.inc"
      return
      end subroutine
c
c
c
#undef _ARRDIMSL_
#define _ARRDIMSL_ :,:,:,:
      subroutine do_4D
#include "gs_setup.inc"
      return
      end subroutine
c
c
c
#undef _ARRDIMSL_
#define _ARRDIMSL_ :,:,:,:,:
      subroutine do_5D
#include "gs_setup.inc"
      return
      end subroutine
