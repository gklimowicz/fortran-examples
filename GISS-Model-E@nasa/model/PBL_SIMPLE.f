#include "rundeck_opts.h"

      MODULE SOCPBL
!@sum  SOCPBL deals with boundary layer physics
!@auth Ye Cheng/G. Hartke (modifications by G. Schmidt)
!@cont pbl,advanc,stars,getl,dflux,simil,griddr,tfix
!@cont ccoeff0,getk,e_eqn,t_eqn,q_eqn,tr_eqn,uv_eqn,level2
!@cont t_eqn_sta,q_eqn_sta,uv_eqn_sta
!@cont inits,tcheck,ucheck,check1,output,rtsafe,fgrid2

      IMPLICIT NONE

      integer, parameter :: n=8  !@param n  no of pbl. layers
      real*8, parameter :: zgs=10. !@var zgs height of surface layer (m)

!input  !@var ZS1    = height of the first model layer (m)
!input !@var TGV    = virtual potential temperature of the ground (K)
!input !@var TKV    = virtual potential temperature of first model layer (K)
!input !@var HEMI   = 1 for northern hemisphere, -1 for southern hemisphere
!input !@var POLE   = .TRUE. if at the north or south pole, .FALSE. otherwise

!@var US     = x component of surface wind, postive eastward (m/s)
!@var VS     = y component of surface wind, positive northward (m/s)
!@var WS     = magnitude of the surface wind (m/s)
!@var WSM    = magnitude of the surface wind - ocean currents (m/s)
!@var WSH    = magnitude of surface wind modified by buoyancy flux(m/s)
!@var TSV    = virtual potential temperature of the surface (K)
!@var QS     = surface value of the specific moisture
!diag !@var PSI    = angular diff. btw geostrophic and surface winds (rads)
!diag !@var DBL    = boundary layer height (m)
!not used @var KMS    = momentum transport coefficient at ZGS (m**2/s)
!@var KHS    = heat transport coefficient at ZGS (m**2/s)
!not used !@var KHQ    = moist transport coefficient at ZGS (m**2/s)
!not used !@var PPBL   = pressure at DBL (mb)
!not used !@var USTAR  = friction speed (square root of momentum flux) (m/s)
!@var CM     = drag coefficient (dimensionless surface momentum flux)
!@var CH     = Stanton number   (dimensionless surface heat flux)
!@var CQ     = Dalton number    (dimensionless surface moisture flux)
!@var z0m   = roughness length for momentum,
!@+           prescribed for itype=3,4 but computed for itype=1,2 (m)
!not used !@var z0h   = roughness length for temperature (m)
!not used !@var z0q   = roughness length for water vapor (m)
!diag !@var UG     = eastward component of the geostrophic wind (m/s)
!diag !@var VG     = northward component of the geostrophic wind (m/s)
!diag !@var WG     = magnitude of the geostrophic wind (m/s)

!input:
      real*8 :: zs1,tgv,tkv,hemi,qg_sat
      logical :: pole
!input compat
      real*8 :: dtsurf
!output
      real*8 :: us,vs,ws,wsm,wsh,tsv,qsrf,khs
     *         ,w2_1
     &     ,wint
!output diag
      real*8 :: psi,dbl,ug,vg,wg
!output compat (not used)
      real*8 :: kms,kqs

      end MODULE SOCPBL


