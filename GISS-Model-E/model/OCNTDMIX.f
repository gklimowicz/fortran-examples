#include "rundeck_opts.h"

      module tdmix_mod
!@auth M. Kelley
!@sum tdmix_mod and associated routines perform mesoscale transport/mixing
!@+   in a dual vertical coordinate which is the union of local isoneutral
!@+   surfaces and model layers.  The duality permits isoneutrally
!@+   stencilled fluxes to be applied directly to model layers.
!@+   Eddy-induced velocities (through the Thickness-Diffusion form
!@+   of GM) are effected via the same flux computations that perform
!@+   slantwise tracer MIXing (hence the module name TDMIX).
!@+
!@+   Implementation and diagnostic notes.
!@+
!@+   Each x- and y-adjacent pair of gridcell columns is decomposed into
!@+   a set of density-space overlap elements (DSOE)s.  Each DSOE
!@+   comprises a left/right pair of water masses having locally
!@+   referenced potential densities which are identical both at
!@+   their upper edges and at their lower edges.  The set of upper-
!@+   and lower-edge DSOE densities is chosen to be the union
!@+   of the two sets of layer-edge potential densities in the
!@+   left and right columns.  This permits the water mass in any
!@+   given model layer to be exactly decomposed as the complete sum of
!@+   a set of (left or right) DSOEs.  This identity also applies to
!@+   tracers, for which conservative vertical remapping is used to
!@+   obtain mean concentrations in each DSOE.
!@+
!@+   Diffusive fluxes are computed pair-by-pair in a loop over DSOEs.
!@+   The bolus component of the fluxes is associated with the difference
!@+   in water masses for each pair, as the diffusion is applied to
!@+   the vertically integrated tracer amount in each DSOE.  Eddy-induced
!@+   horizontal velocities are thus diagnosed from the DSOE-indexed fluxes
!@+   of a "water mass" tracer having unit concentration.
!@+
!@+   The loop over DSOEs effects the entirety of the horizontal fluxes, but
!@+   only part of the vertical fluxes. The DSOE-step convergence of "water mass"
!@+   fluxes induces a vertical displacement of model layer interfaces and thus
!@+   an associated tracer flux in an Eulerian (fixed) coordinate system.
!@+   For economy and accuracy, the calculation of this tracer flux is delegated
!@+   and deferred to the vertical advection step including resolved velocities.
!@+   This combination of velocities requires a post-hoc adjustment of vertical
!@+   flux diagnostics which is inferred periodically on model layers from the
!@+   continuity equation and the saved convergence of water mass fluxes.
!@+
!@+   A literal application of thickness diffusion will generate a nonzero
!@+   barotropic component of eddy-induced flow in the presence of gradients
!@+   of column mass (i.e. bathymetry and surface height); vertical gradients
!@+   of diffusivity also contribute.  Exact baroclinicity is enforced via
!@+   a combination of three mechanisms.  Firstly, the part of the deeper
!@+   column in each pair which is below the bottom of the shallower column
!@+   only participates in the DSOE loop if its density is less than that
!@+   at the bottom of the shallower column.  (In other words, denser fluid
!@+   does not flow upslope.)  Secondly, as outlined below, vertically varying
!@+   diffusivity is presented to the DSOE method via a fractional-in-time
!@+   but instantaneously constant-where-nonzero sampling scheme.  Thirdly, an
!@+   adjustment of mass exchange rates is applied that is small except where
!@+   surface height gradients are large or dense water is bolus-fluxed downslope.
!@+
!@+   The ratio of model timestep to eddy turnover timescale is favorable
!@+   to fractional-in-time application of mesoscale diffusivity.  Vertical
!@+   variations of instantaneous diffusivity are converted into "on/off"
!@+   temporal variations whose "on" frequency at a given depth is equal
!@+   to the ratio of diffusivity at that depth to the peak value in the
!@+   column.  The peak value is used within all depth zones that are "on"
!@+   during a given timestep.  The smoothness and monotonicity of the
!@+   diffusivity profile determines the number and position of such zones.
!@+   At the bottom (or top) of each "on" zone, unpaired DSOEs produce
!@+   a bolus transport which is the counterpart of s*dK/dz in the bolus
!@+   streamfunction psi = d(s*K)/dz.
!@+   For longer timesteps, a variation of this approach could be used
!@+   which generates a set of DSOE mass exchanges by averaging multiple
!@+   applications of the above protocol.
!@+
!@+   Upper boundary conditions and tapering.  The vertical orientation
!@+   of isoneutral surfaces in the mixed layer implies that a literal
!@+   application of thickness diffusion in the ML would produce
!@+   subduction of ML water at a rate dependent on horizontal
!@+   resolution.  Since eddy vertical velocity must vanish at the sea
!@+   surface, DSOE pairings in the ML are modified such that subducting
!@+   and ascending water parcels are re-directed horizontally. In this
!@+   manner, eddy activity in the ML remains at unmodified strength
!@+   but thickness diffusion is transformed into purely horizontal
!@+   mixing. In this context, the ML depth is defined as the vertical
!@+   extent of unpaired DSOEs in the lighter ML.  No post-hoc
!@+   adjustment of the ML diffusivity is performed; it is considered
!@+   the responsibility of the diffusivity calculation scheme to
!@+   prescribe variations of mixing vigor with depth.
!@+
!@+   Definition of potential density for DSOEs.  Since DSOEs are
!@+   constructed independently for each column pair, it is not
!@+   necessary to employ a globally constant reference pressure.
!@+   Instead, each column pair uses a local profile of a "pseudo"
!@+   density obtained via downward integration of the vertical
!@+   gradient of potential density evaluated using in-situ pressure.
!@+   This density coordinate is monotonic as long as local
!@+   stratification remains positive; the pseudo-density is
!@+   monotonized where stratification is negative.
!@+   
!@+   Some comparisons to other numerical forms of mesoscale
!@+   transport in z/p-coordinate models.
!@+   (1)  No column-local vertical transport terms (e.g. those
!@+        that would arise from rotation of the isoneutral mixing
!@+        tensor).  All vertical transport occurs either (a) slantwise
!@+        in DSOE pairs having different layer indices or (b) via the
!@+        induced displacement of model layer interfaces.
!@+   (2)  In constrast to the purely advective form of GM, only
!@+        term (b) in (1) above admits use of a high-order
!@+        and/or large-stencil and/or upwind advection scheme.
!@+        Since the net horizontal mass flux is being modeled via
!@+        diffusion, i.e. bi-directional flow, it is unclear whether
!@+        sharp features that would be preserved by a high-order scheme
!@+        are phenomenologically consistent (also given the presence
!@+        of simultaneous "Redi" mixing which would smear them out).
!@+        However, use of a high-order scheme for effect (b) above is
!@+        justified given that eddy motions are quasi-horizontal.
!@+   (3)  Large isopycnal slopes.  The DSOE method makes no
!@+        explicit reference to slopes since it does not involve
!@+        a rotated mixing tensor.  Large slopes therefore induce
!@+        no numerical difficulties; mass exchanges can "jump" an
!@+        arbitrary number of model layers.  Rather, adaptation of
!@+        the scheme to this limit is associated with the introduction
!@+        of mass-exchange elements that "jump" in density space
!@+        for physical reasons, e.g. as per the discussion of upper
!@+        boundary conditions.
!@+
!@+   To-dos and future areas of exploration:
!@+   (1)  Minor mods for unequal thickness versus tracer diffusivities
!@+        (GM versus Redi) when a compelling model emerges.
!@+   (2)  Accommodate arbitrary layering.
!@+   (3)  Protocol to allow external code to specify profile of
!@+        arbitrarily slanted mass exchange; over multiple angles.
!@+   (4)  Larger stencils for various vertical interpolations/remaps.
!@+   (5)  Use (reconstructed) horizontal subgrid info.

      use ocean, only : lmo
      implicit none
      integer, parameter :: lmaxoverlap=2*lmo+10

!@var dzo_ end-padded version of nominal layer thickness array
      real*8, dimension(0:lmo+1) :: dzo_
!@var by6arr,dzby12,bydzoe geometric factors used for updating
!@+   prognostic subgrid structure
      real*8, dimension(lmo) :: by6arr,dzby12,bydzoe

!@var lm[xy] number of DSOEs between x- and y-adjacent columns
      integer, dimension(:,:), allocatable ::
     &     lmx,
     &     lmy
!@var l[lr][xy] left and right layer indices for x- and y- direction DSOEs
      integer, dimension(:,:,:), allocatable ::
     &     llx,lrx,
     &     lly,lry
!@var dm[lr][xy] left and right water masses of x- and y-direction DSOEs (kg)
!@var c{z,zz}[lr][xy] left and right, x- and y- direction coefficients of
!@+   proportionality relating tracer concentration in each DSOE to the
!@+   1st and 2nd z-derivatives of tracer concentration in the layer
!@+   enclosing the DSOE
      real*8, dimension(:,:,:), allocatable ::
     &     dmlx,dmrx, czlx,czzlx,czrx,czzrx,
     &     dmly,dmry, czly,czzly,czry,czzry

!@var mokgsv pre-tdmix gridbox mass (kg), used in adjustment of subgrid tracer moments
      real*8, dimension(:,:,:), allocatable :: mokgsv
!@var mrel[xyz] relaxation factors for adjustment of subgrid tracer moments
      real*8, dimension(:,:,:), allocatable :: mrelx,mrely,mrelz

#ifdef TRACERS_OCEAN
!@var a[...] accumulation variables corresponding to variable names [...] above.
!@+   Only used when ntrtrans>1.   See description of ntrtrans and long-timestep
!@+   transport in module OCEAN.
!@var almaxoverlap is only an initial DSOE count before on-the-fly resizing
      integer, parameter :: almaxoverlap=3*lmo
      integer, dimension(:,:), allocatable ::
     &     almx,
     &     almy
      integer, dimension(:,:,:), allocatable ::
     &     allx,alrx,
     &     ally,alry
      real*8, dimension(:,:,:), allocatable ::
     &     admlx,admrx, aczlx,aczzlx,aczrx,aczzrx,
     &     admly,admry, aczly,aczzly,aczry,aczzry
#endif

      contains

      subroutine get_dm_intervals(x1,n1,x2,n2,dx1,dx2,l1,l2,
     &     cx1,cxx1, cx2,cxx2,
     &     nxchng,nxupl,nxupr,nxfull,
     &     nmax)
c
c Given two lists of points x1(1:n1+1) and x2(1:n2+1)
c along the "x" axis, merge the two lists and calculate the
c distance increments dx(1:nxchng) separating the points
c in the merged list.  The i_th increment lies in the intervals
c x1(l1(i):l1(i)+1) and x2(l2(i):l2(i)+1).
c Points to the left of max(x1(1),x2(1)) are excluded.
c Unpaired intervals to the right of min(x1(n1+1),x2(n2+1)) are
c associated with the last interval of x1 if x1(n1+1) < x2(n2+1)
c and vice versa.
c
      implicit none
      integer :: n1,n2,nxchng,nxupl,nxupr,nxfull,nmax
      real*8, dimension(n1+1) :: x1
      real*8, dimension(n2+1) :: x2
      real*8, dimension(nmax) :: dx1,dx2
      integer, dimension(nmax) :: l1,l2
      real*8, dimension(nmax) :: cx1,cxx1, cx2,cxx2
      integer :: ls,ld,i,l,nmin
      real*8 :: fac1,fac2
      real*8, allocatable :: x(:)

      nxchng = 0
      nxupr = 0
      nxupl = 0
      nmin = min(n1,n2)
      if(min(x1(n1+1),x2(n2+1)).le.max(x1(1),x2(1))) then
        ! no overlap whatsoever
        if(x1(1).gt.x2(1)) then ! .and. n1.le.n2) then
          nxupr = nmin
          do i=1,nmin
            l1(i) = i
            l2(i) = n2
            dx1(i) = x1(i+1)-x1(i)
            dx2(i) = 0.
            cx1(i) = 0.
            cxx1(i) = 0.
            cx2(i) = 0.
            cxx2(i) = 0.
          enddo
          i = nmin
          nxupl = nmin
          do l=1,nmin
            i = i + 1
            l1(i) = l
            l2(i) = l
            dx1(i) = (x1(l+1)-x1(l))
            dx2(i) = (x2(l+1)-x2(l))*2.
            cx1(i) = 0.
            cxx1(i) = 0.
            cx2(i) = 0.
            cxx2(i) = 0.
          enddo
        elseif(x2(1).gt.x1(1) ) then !.and. n2.le.n1) then
          nxupr = nmin
          do i=1,nmin
            l1(i) = n1
            l2(i) = i
            dx1(i) = 0.
            dx2(i) = x2(i+1)-x2(i)
            cx1(i) = 0.
            cxx1(i) = 0.
            cx2(i) = 0.
            cxx2(i) = 0.
          enddo
          i = nmin
          nxupl = nmin
          do l=1,nmin
            i = i + 1
            l1(i) = l
            l2(i) = l
            dx1(i) = (x1(l+1)-x1(l))*2.
            dx2(i) = (x2(l+1)-x2(l))
            cx1(i) = 0.
            cxx1(i) = 0.
            cx2(i) = 0.
            cxx2(i) = 0.
          enddo
        endif
      else
        allocate(x(nmax+1))
        x(1) = max(x1(1),x2(1))
        ls = 2
        ld = 2
        if(x1(1).gt.x2(1)) then
          do ld=2,n2+1
            if(x2(ld).gt.x(1)) exit
          enddo
        elseif(x1(1).lt.x2(1)) then
          do ls=2,n1+1
            if(x1(ls).gt.x(1)) exit
          enddo
        endif
        do nxchng=1,nmax
          l1(nxchng) = ls-1
          l2(nxchng) = ld-1
          if(x2(ld).le.x1(ls)) then
            if(x2(ld).eq.x1(ls)) ls = ls + 1
            x(nxchng+1) = x2(ld)
            ld = ld + 1
          else
            x(nxchng+1) = x1(ls)
            ls = ls + 1
          endif
          if(ld.gt.n2+1) exit
          if(ls.gt.n1+1) exit
        enddo
        nxfull = nxchng
        dx1(1:nxchng) = x(2:nxchng+1)-x(1:nxchng)
        dx2(1:nxchng) = x(2:nxchng+1)-x(1:nxchng)
        do i=1,nxchng
          if(l1(i).gt.nmin) dx1(i)=0.
          if(l2(i).gt.nmin) dx2(i)=0.
        enddo
        if(    x1(n1+1).gt.x2(n2+1) ) then !.and. n1.le.n2) then
        ! associate unpaired intervals in x1 with the last x2 interval
          nxupr = 0
          do i=nxchng+1,nxchng + n1-l1(nxchng)+1
            l2(i) = n2
            l1(i) = l1(nxchng) + (i-nxchng-1)
            !if(n1.gt.n2) exit
            if(l1(i).gt.n2) exit
            x(i+1) = x1(l1(i)+1)
            dx1(i) = x(i+1)-x(i)
            dx2(i) = 0.
            cx1(i) = 0.
            cxx1(i) = 0.
            cx2(i) = 0.
            cxx2(i) = 0.
            nxupr = nxupr + 1
          enddo
          nxfull = nxchng + nxupr
        elseif(x1(n1+1).lt.x2(n2+1) ) then !.and. n2.le.n1) then
        ! associate unpaired intervals in x2 with the last x1 interval
          nxupr = 0
          do i=nxchng+1,nxchng + n2-l2(nxchng)+1
            l1(i) = n1
            l2(i) = l2(nxchng) + (i-nxchng-1)
            !if(n2.gt.n1) exit
            if(l2(i).gt.n1) exit
            x(i+1) = x2(l2(i)+1)
            dx1(i) = 0.
            dx2(i) = x(i+1)-x(i)
            cx1(i) = 0.
            cxx1(i) = 0.
            cx2(i) = 0.
            cxx2(i) = 0.
            nxupr = nxupr + 1
          enddo
          nxfull = nxchng + nxupr
        endif
c
        i = nxfull
        if(l1(1).ne.l2(1)) then
          if(l1(1).gt.l2(1)) then
            fac1 = 1.
            fac2 = 0.
          else
            fac1 = 0.
            fac2 = 1.
          endif
          nxupl = max(0,min(max(l1(1),l2(1))-1,min(n1,n2)))
          do l=1,nxupl
            i = i + 1
            l1(i) = l
            l2(i) = l
            dx1(i) = (x1(l+1)-x1(l))*fac1
            dx2(i) = (x2(l+1)-x2(l))*fac2
            cx1(i) = 0.
            cxx1(i) = 0.
            cx2(i) = 0.
            cxx2(i) = 0.
          enddo
        else
          nxupl = 0
        endif
        if(x1(1).gt.x2(1)) then
          nxupl = nxupl + 1
          i = i + 1
          l = l2(1)
          dx1(i) = 0.
          dx2(i) = x(1)-x2(l)
          l1(i) = min(l,n1)!1
          l2(i) = l
          cx1(i) = 0.
          cxx1(i) = 0.
          call get_coeffs(x2(l),x2(l+1),x2(l),x(1),cx2(i),cxx2(i))
        elseif(x1(1).lt.x2(1)) then
          nxupl = nxupl + 1
          i = i + 1
          l = l1(1)
          dx1(i) = x(1)-x1(l)
          dx2(i) = 0.
          l1(i) = l
          l2(i) = min(l,n2)!1
          call get_coeffs(x1(l),x1(l+1),x1(l),x(1),cx1(i),cxx1(i))
          cx2(i) = 0.
          cxx2(i) = 0.
        endif
        nxfull = i
c
        do i=1,nxchng
          l = l1(i)
          call get_coeffs(x1(l),x1(l+1),x(i),x(i+1),cx1(i),cxx1(i))
          l = l2(i)
          call get_coeffs(x2(l),x2(l+1),x(i),x(i+1),cx2(i),cxx2(i))
        enddo
        deallocate(x)
      endif
c
      nxfull = nxchng + nxupr + nxupl
      return
      contains
      subroutine get_coeffs(xlb,xrb,xl,xr,cx,cxx)
      real*8 :: xlb,xrb,xl,xr
      real*8 :: cx,cxx
      real*8 :: bydxb,zl,zr
      bydxb = 1d0/(xrb-xlb)
      zl = 2.*(xl-xlb)*bydxb-1d0
      zr = 2.*(xr-xlb)*bydxb-1d0
      cx = .5*(zl+zr)
      !cxx = .5*((zr**3-zl**3)/(zr-zl) -1d0)
      cxx = .5*( (zl+zr)**2 -zl*zr -1d0)
      end subroutine get_coeffs
      end subroutine get_dm_intervals

      subroutine monotonize(q,z,n,dqdzmin)
      implicit none
      integer :: n
      real*8 :: q(n),z(n),dqdzmin
c
      integer :: imin,imax,i,j,k

      real*8 :: bydztot,wt

c
c restore monotonicity at the left edge
c
      imin = sum(minloc(q))
      do j=1,imin-1
        q(j) = q(imin)
      enddo
c
c restore monotonicity at the right edge
c
      imax = sum(maxloc(q))
      do j=imax+1,n
        q(j) = q(imax)
      enddo

c
c restore monotonicity in between
c
      j=imin
      do while(j.lt.imax-1)
        if(q(j+1).lt.q(j)) then
          k=j+2
          do while(q(k).lt.q(j))
            k=k+1
          enddo
          bydztot = 1d0/(z(k)-z(j))
          do i=j+1,k-1
            wt = (z(k)-z(i))*bydztot
            q(i) = wt*q(j) + (1.-wt)*q(k)
          enddo
          j=k
        else
          j=j+1
        endif
      enddo

c
c add background dq/dz
c
      if(dqdzmin.gt.0.) then
        do j=2,n
          q(j) = q(j) + (z(j)-z(1))*dqdzmin
        enddo
      endif

      return
      end subroutine monotonize

      subroutine make_left_right(gl,sl,pl,ml,gr,sr,pr,mr,
     &     rle,rre,dmdrl,dmdrr,lml,lmr)
      use ocean, only : lmo,ze,dzo
      implicit none
      real*8, dimension(lmo) :: gl,sl,pl,ml,gr,sr,pr,mr
      real*8, dimension(0:lmo) :: rle,rre
      real*8, dimension(lmo) :: dmdrl,dmdrr
      integer :: lml,lmr
c
      real*8 :: volgsp ! function
      integer :: l,lmu
      real*8, dimension(lmo) :: g,s,p,r,rl,rr
      real*8 :: drho,wtdn,wtup,pavg
      real*8, parameter :: drdzmin=1d-6
c
      lmu = min(lml,lmr)
c
      do l=1,lmu
        g(l) = .5*(gl(l)+gr(l))
        s(l) = .5*(sl(l)+sr(l))
        p(l) = .5*(pl(l)+pr(l))
      enddo
      r(1) = 1d0/volgsp(g(1),s(1),p(1))
      do l=2,lmu
        pavg = .5d0*(p(l)+p(l-1))
        drho = 
     &       ( 1d0/volgsp(g(l  ),s(l  ),pavg)
     &        -1d0/volgsp(g(l-1),s(l-1),pavg) )
        r(l) = r(l-1) + drho
      enddo
      do l=1,lmu
        drho = 
     &       ( 1d0/volgsp(gr(l),sr(l),p(l))
     &        -1d0/volgsp(gl(l),sl(l),p(l)) )
        rl(l) = r(l) - .5*drho
        rr(l) = r(l) + .5*drho
      enddo
      if(lmr.gt.lml) then
        do l=lmu+1,lmr
          pavg = .5d0*(pr(l)+pr(l-1))
          drho = 
     &       ( 1d0/volgsp(gr(l  ),sr(l  ),pavg)
     &        -1d0/volgsp(gr(l-1),sr(l-1),pavg) )
          rr(l) = rr(l-1) + drho
        enddo
      elseif(lmr.lt.lml) then
        do l=lmu+1,lml
          pavg = .5d0*(pl(l)+pl(l-1))
          drho = 
     &       ( 1d0/volgsp(gl(l  ),sl(l  ),pavg)
     &        -1d0/volgsp(gl(l-1),sl(l-1),pavg) )
          rl(l) = rl(l-1) + drho
        enddo
      endif
      do l=1,lml-1
        wtup = dzo(l+1)/(dzo(l)+dzo(l+1))
        wtdn = 1.-wtup
        rle(l) = wtup*rl(l)+wtdn*rl(l+1)
      enddo
      do l=1,lmr-1
        wtup = dzo(l+1)/(dzo(l)+dzo(l+1))
        wtdn = 1.-wtup
        rre(l) = wtup*rr(l)+wtdn*rr(l+1)
      enddo
      rle(0) = rl(1)-(rle(1)-rl(1))
      rre(0) = rr(1)-(rre(1)-rr(1))
      rle(lml) = rl(lml)+(rl(lml)-rle(lml-1))
      rre(lmr) = rr(lmr)+(rr(lmr)-rre(lmr-1))
      call monotonize(rle,ze,lml+1,drdzmin)
      call monotonize(rre,ze,lmr+1,drdzmin)
      do l=1,lml
        dmdrl(l) = ml(l)/(rle(l)-rle(l-1))
      enddo
      do l=1,lmr
        dmdrr(l) = mr(l)/(rre(l)-rre(l-1))
      enddo
      return
      end subroutine make_left_right

      subroutine get_dsoe_info(
     &     cdir,
     &     lmml,lmmr,idebug,jdebug,
     &     ml,mr,gl,gr,sl,sr,pl,pr,
     &     nlzone,ltb,
     &     lm,
     &     llout,lrout,
     &     dml,dmr,czl,czzl,czr,czzr
     &     )
      use ocean, only : lmo
      implicit none
c input
      character(len=1) :: cdir
      integer :: lmml,lmmr
      integer :: idebug,jdebug
      real*8, dimension(lmo) :: ml,mr,gl,gr,sl,sr,pl,pr
      integer :: nlzone
      integer, dimension(2,lmo/2) :: ltb
c output
      integer :: lm
      integer, dimension(lmaxoverlap) :: llout,lrout
      real*8, dimension(lmaxoverlap) :: dml,dmr,czl,czzl,czr,czzr
c
      integer, parameter :: lmax=2*lmaxoverlap ! just in case
      real*8, dimension(lmax) :: drl,drr
      real*8, dimension(lmo) :: dmdrl,dmdrr
      real*8, dimension(0:lmo) :: rle,rre
      real*8 :: suml,sumr
      integer :: l,n,ll,lr,lmovlap,nlml,nlbbc,lml,lmr,ltop,lbot,lm_
      integer :: lzone
      integer, dimension(lmax) :: ll_,lr_
      real*8, dimension(lmax) :: dml_,dmr_,czl_,czzl_,czr_,czzr_

      real*8 :: delm,redirect_mlsub_eff
      logical :: redirect_mlsub

      redirect_mlsub = .true.
      redirect_mlsub_eff = 1d0

      do lzone=1,nlzone

        ltop = ltb(1,lzone)
        lbot = ltb(2,lzone)

        if(ltop .ge. min(lmml,lmmr)) cycle

        lml = min(lbot,lmml) -ltop + 1
        lmr = min(lbot,lmmr) -ltop + 1
        call make_left_right(
     &       gl(ltop),sl(ltop),pl(ltop),ml(ltop),
     &       gr(ltop),sr(ltop),pr(ltop),mr(ltop),
     &       rle,rre,dmdrl,dmdrr,
     &       lml,lmr )
        call get_dm_intervals(
     &       rle,lml,
     &       rre,lmr,
     &       drl,drr,
     &       ll_,lr_,
     &       czl_,czzl_,czr_,czzr_,
     &       lmovlap,nlml,nlbbc,lm_,lmax)
        if(lm_.gt.lmaxoverlap) then
          write(6,*) 'for '//cdir//'-direction:'
          write(6,*) 'lm>lmaxoverlap at i,j = ',idebug,jdebug
          write(6,*) 'increase lmaxoverlap if reasonable'
          call stop_model('lm>lmaxoverlap',255)
        endif
        do l=1,lm_
          ll = ll_(l)
          lr = lr_(l)
          dml_(l) = dmdrl(ll)*drl(l)
          dmr_(l) = dmdrr(lr)*drr(l)
        enddo

        ! See comments on upper boundary conditions and ML subduction
        ! at beginning of module.
        if(redirect_mlsub .and. (nlml.gt.1 .or. ltop.eq.1)) then
          if(rle(0) .lt. rre(0)) then ! left ml is lighter
            do l=1,lmovlap
              if(lr_(l).gt.nlml) exit ! is this right?
              ll = lr_(l) + lmovlap+nlbbc
              delm = (dmr_(l)-dml_(l))*redirect_mlsub_eff
              if(delm.lt.0.) cycle ! no net subduction
              dmr_(ll) = dmr_(ll) + delm
              dmr_(l ) = dmr_(l ) - delm
            enddo
          elseif(rle(0) .gt. rre(0)) then ! right ml is lighter
            do l=1,lmovlap
              if(ll_(l).gt.nlml) exit ! is this right?
              ll = ll_(l) + lmovlap+nlbbc
              delm = (dml_(l)-dmr_(l))*redirect_mlsub_eff
              if(delm.lt.0.) cycle ! no net subduction
              dml_(ll) = dml_(ll) + delm
              dml_(l ) = dml_(l ) - delm
            enddo
          endif
        endif

        ll = lm
        do l=1,lm_
          ll = ll + 1
          dml(ll) = dml_(l)
          dmr(ll) = dmr_(l)
          llout(ll) = ll_(l) + (ltop - 1)
          lrout(ll) = lr_(l) + (ltop - 1)
          czl(ll) = czl_(l)
          czzl(ll) = czzl_(l)
          czr(ll) = czr_(l)
          czzr(ll) = czzr_(l)
        enddo
        lm = lm + lm_

      enddo ! lzone

      ! Multiplicatively adjust either the left- or right-side
      ! mass exchange elements to enforce zero net vertically
      ! integrated mass flux
      if(lm.gt.0) then
        suml = sum(dml(1:lm))
        sumr = sum(dmr(1:lm))
        if(sumr.gt.suml) then
          dmr(1:lm) = dmr(1:lm)*suml/sumr
        elseif(suml.gt.sumr) then
          dml(1:lm) = dml(1:lm)*sumr/suml
        endif
      endif

      end subroutine get_dsoe_info

      subroutine get_frt_intervals(xrand,lmo,krat,nlzone,ltb)
!@sum get_frt_intervals tablulates the number of intervals in array
!@+   krat which are greater than or equal to xrand, and the
!@+   top and bottom indices of each interval.
      implicit none
      real*8 :: xrand
      integer :: lmo
      real*8, dimension(lmo) :: krat
      integer :: nlzone
      integer, dimension(2,lmo/2) :: ltb
c
      integer :: l,lzone
      if(minval(krat(1:2)).ge.xrand) then
        nlzone = 1
        ltb(1,nlzone) = 1
      else
        nlzone = 0
      endif
      do l=1,lmo-2
        if(krat(l).lt.xrand .and. minval(krat(l+1:l+2)).ge.xrand) then
          nlzone = nlzone + 1
          ltb(1,nlzone) = l+1
        endif
      enddo
      if(minval(krat(lmo-1:lmo)).ge.xrand) ltb(2,nlzone) = lmo
      do lzone=1,nlzone
        do l=2+ltb(1,lzone),lmo
          if(krat(l).lt.xrand) then
            ltb(2,lzone) = l-1
            exit
          endif
        enddo
      enddo
      end subroutine get_frt_intervals

      subroutine get_zmoms(r,rz,rzz,lmom)
      ! low-order version which allows discontinuities at edges
      use ocean, only: lmo
      implicit none
      integer :: lmom
      real*8, dimension(lmom) :: r,rz,rzz
C**** Local variables
      integer :: l
      real*8, dimension(0:lmo) :: rzgrad
      do l=1,lmom-1
        rzgrad(l) = (r(l+1)-r(l))*bydzoe(l)
      enddo
      rzgrad(0) = rzgrad(1)
      rzgrad(lmom) = rzgrad(lmom-1)
      do l=1,lmom
        rz(l) = by6arr(l)*(
     &       dzo_(l+1)*rzgrad(l-1)+dzo_(l-1)*rzgrad(l)
     &       + .5d0*dzo_(l)*(rzgrad(l-1)+rzgrad(l))  )
        rzz(l) = dzby12(l)*(rzgrad(l)-rzgrad(l-1))
      enddo
      return
      end subroutine get_zmoms

      subroutine adj_zmoms(r,rz,rzz,lr,cz,czz,lmo,lmx)
      implicit none
      integer :: lmo,lmx
      real*8, dimension(lmo) :: r,rz,rzz
      integer, dimension(lmx) :: lr
      real*8, dimension(lmx) :: cz,czz
c
      integer :: l,ll
      real*8 :: rr,fac
c
      do l=1,lmx
        ll = lr(l)
        rr = rz(ll)*cz(l) +rzz(ll)*czz(l)
        if(rr + r(ll) .lt. 0.) then
          fac = -max(0d0,r(ll))/(min(0d0,rr)-1d-30)
          rz(ll) = rz(ll)*fac
          rzz(ll) = rzz(ll)*fac
        endif
      enddo
c
      return
      end subroutine adj_zmoms

      end module tdmix_mod

      subroutine alloc_tdmix
      use tdmix_mod
      use domain_decomp_1d, only : getdomainbounds
      use ocean, only : im,lmo,dzo
      use oceanr_dim, only : grid=>ogrid
#ifdef TRACERS_OCEAN
      use ocean, only : ntrtrans
#endif
      implicit none
      integer :: l,j_0h,j_1h

      call getdomainbounds(grid, j_strt_halo=j_0h, j_stop_halo=j_1h)

      call realloc_tdmix_regstep

      allocate(
c
     &        mokgsv(im,j_0h:j_1h,lmo),
c
     &        mrelx(im,j_0h:j_1h,lmo),
     &        mrely(im,j_0h:j_1h,lmo),
     &        mrelz(im,j_0h:j_1h,lmo)
     &        )

      mokgsv = 0.
      mrelx = 0.
      mrely = 0.
      mrelz = 0.

      dzo_(1:lmo) = dzo(1:lmo)
      dzo_(0) = dzo_(1)
      dzo_(lmo+1) = dzo_(lmo)
      do l=1,lmo
        by6arr(l) = .5d0*dzo_(l)/sum(dzo_(l-1:l+1))
        dzby12(l) = .5d0*by6arr(l)*dzo_(l)
      enddo
      do l=1,lmo-1
        bydzoe(l) = 2d0/(dzo(l)+dzo(l+1))
      end do

#ifdef TRACERS_OCEAN
      if(ntrtrans.gt.1) then
        call tdmix_longstep_finish ! to initialize to starting values
      endif
#endif

      return
      end subroutine alloc_tdmix

      subroutine realloc_tdmix_regstep
      use tdmix_mod
      use domain_decomp_1d, only : getdomainbounds
      use ocean, only : im,lmo
      use oceanr_dim, only : grid=>ogrid
      implicit none
      integer :: j_0h,j_1h

      call getdomainbounds(grid, j_strt_halo=j_0h, j_stop_halo=j_1h)

      if(allocated(llx)) then
        deallocate(
     &         lmx,
     &         lmy,
     &         llx,
     &         lrx,
     &         lly,
     &         lry,
c
     &         dmlx,
     &         dmrx,
     &         dmly,
     &         dmry,
c
     &          czlx,
     &         czzlx,
     &          czrx,
     &         czzrx,
     &          czly,
     &         czzly,
     &          czry,
     &         czzry
     &        )
      endif

      allocate(
     &         lmx(im,j_0h:j_1h),
     &         lmy(im,j_0h:j_1h),
     &         llx(lmaxoverlap,im,j_0h:j_1h),
     &         lrx(lmaxoverlap,im,j_0h:j_1h),
     &         lly(lmaxoverlap,im,j_0h:j_1h),
     &         lry(lmaxoverlap,im,j_0h:j_1h),
c
     &         dmlx(lmaxoverlap,im,j_0h:j_1h),
     &         dmrx(lmaxoverlap,im,j_0h:j_1h),
     &         dmly(lmaxoverlap,im,j_0h:j_1h),
     &         dmry(lmaxoverlap,im,j_0h:j_1h),
c
     &          czlx(lmaxoverlap,im,j_0h:j_1h),
     &         czzlx(lmaxoverlap,im,j_0h:j_1h),
     &          czrx(lmaxoverlap,im,j_0h:j_1h),
     &         czzrx(lmaxoverlap,im,j_0h:j_1h),
     &          czly(lmaxoverlap,im,j_0h:j_1h),
     &         czzly(lmaxoverlap,im,j_0h:j_1h),
     &          czry(lmaxoverlap,im,j_0h:j_1h),
     &         czzry(lmaxoverlap,im,j_0h:j_1h)
     &        )

      lmx = 0; lmy = 0; llx = 0; lrx = 0; lly = 0; lry = 0
      dmlx = 0.; dmrx = 0.; dmly = 0.; dmry = 0.

      czlx = 0.; czzlx = 0.; czrx = 0.; czzrx = 0. 
      czly = 0.; czzly = 0.; czry = 0.; czzry = 0.

      end subroutine realloc_tdmix_regstep

      subroutine tdmix_prep(dt,k3dx,k3dy)
!@sum tdmix_prep prepares the DSOE mass-exchange indices/rates as
!@+   described in the comments at the beginning of this file
      use tdmix_mod
      use ocean, only : mo
      use ocean, only : g3d,s3d,p3d
      use ocean, only : im,jm,lmo
      use ocean, only : dxpo,dyvo,dypo,dxvo,cospo,sinpo,dxypo
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use ocean, only : nbyzu,i1yzu,i2yzu,lmu
      use ocean, only : nbyzv,i1yzv,i2yzv,lmv
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds
      use random, only : randu
#ifdef TRACERS_OCEAN
      use ocean, only : ntrtrans
#endif
      implicit none
!@var dt timestep (s)
      real*8 :: dt  ! timestep (s)
!@var k3dx,k3dy (m2/s) diffusvity at x- and y-edges
      real*8, dimension(lmo,im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     k3dx,k3dy
c
      real*8, dimension(lmo) :: krat,k1d
      real*8, dimension(lmo) :: ml,mr,gl,gr,sl,sr,pl,pr
      integer :: i,j,l,n,ll,lr,il,ir,jl,jr
      real*8 :: facj
      real*8 :: xrand,xdum
      integer :: nlzone
      integer, dimension(2,lmo/2) :: ltb

      integer :: j_0,j_1,j_0s,j_1s
      logical :: have_north_pole

      real*8 :: kmax

#ifdef TRACERS_OCEAN
      logical :: found
      integer :: al,lmx_,lmy_
#endif

      call getdomainbounds(grid, j_strt=j_0, j_stop=j_1,
     &               j_strt_skp=j_0s, j_stop_skp=j_1s,
     &               have_north_pole=have_north_pole)


      do l=1,lmo
      do j=j_0,j_1
      do n=1,nbyzm(j,l)
      do i=i1yzm(n,j,l),i2yzm(n,j,l)
        mokgsv(i,j,l) = mo(i,j,l)*dxypo(j)
      enddo
      enddo
      enddo
      enddo

! Obtain random number on [0:1] interval for fractional-in-time treatment
! of vertical variation of diffusivity (see notes at beginning of tdmix_mod).
! Regular cycling over [0:1] is mathematically equivalent but might
! introduce spurious correlations for an unlucky choice of the cycling period.
      !call random_number(xrand)
      xrand = randu(xdum) ! use host PRNG for reproducibility

      do j=j_0s,j_1s
      facj = dt*dypo(j)/dxpo(j)
      if(cospo(j).lt..5d0) ! beta version did this for safety - not clear
     & facj = facj*(cospo(j)/.5d0) ! if still necessary
      do n=1,nbyzu(j,1)
      do i=i1yzu(n,j,1),i2yzu(n,j,1)
        il = i
        if(i.eq.im) then
          ir = 1
        else
          ir = i+1
        endif

        lmx(i,j) = 0

        k1d(:) = k3dx(:,i,j)
        k1d(1) = k1d(2) ! for shallow oceans
        kmax = maxval(k1d)

        krat = k1d / kmax ! normalized K-profile for fractional-in-time
        call get_frt_intervals(xrand,lmo,krat,nlzone,ltb) ! mixing zones
        if(nlzone.eq.0) cycle

        do l = 1,ltb(2,nlzone) ! fill input arrays to maximum depth of any zone
          ml(l) = mo(il,j,l)
          mr(l) = mo(ir,j,l)
          gl(l) = g3d(l,il,j)
          gr(l) = g3d(l,ir,j)
          sl(l) = s3d(l,il,j)
          sr(l) = s3d(l,ir,j)
          pl(l) = p3d(l,il,j)
          pr(l) = p3d(l,ir,j)
        enddo

        call get_dsoe_info(
     &       'x',
     &       lmm(il,j),lmm(ir,j),i,j,
     &       ml,mr,gl,gr,sl,sr,pl,pr,
     &       nlzone,ltb,
     &       lmx(i,j),
     &       llx(1,i,j),lrx(1,i,j),
     &       dmlx(1,i,j),dmrx(1,i,j),
     &       czlx(1,i,j),czzlx(1,i,j),czrx(1,i,j),czzrx(1,i,j)
     &     )

        ! rescale DSOE mass dmlx (kg/m2) -> ((K/dx)*dy*dt)*dmlx = (kg)
        do l=1,lmx(i,j)
          dmlx(l,i,j) = dmlx(l,i,j)*facj*kmax
          dmrx(l,i,j) = dmrx(l,i,j)*facj*kmax
        enddo

#ifdef TRACERS_OCEAN
        ! accumulate DSOE information from this timestep
        if(ntrtrans.gt.1) then
        lmx_ = size(allx,1)
        do l=1,lmx(i,j)
          found = .false.
          do al=1,almx(i,j)
            if(llx(l,i,j).eq.allx(al,i,j) .and.
     &         lrx(l,i,j).eq.alrx(al,i,j)) then
              found = .true.
              exit
            endif
          enddo
          if(.not. found) then
            almx(i,j) = almx(i,j) + 1
            al = almx(i,j)
            if(al.gt.12*lmo) then
              !write(6,*) 'for '//cdir//'-direction:'
              write(6,*) 'lm>12*lmo at i,j = ',i,j
              write(6,*) 'increase limit if reasonable'
              call stop_model('lm>12*lmo (tracers x-dir)',255)
            endif
            if(al.gt.lmx_) then ! resize on the fly
              lmx_ = lmx_ + lmo
              call growi(lmx_,allx)
              call growi(lmx_,alrx)
              call growr(lmx_,admlx)
              call growr(lmx_,admrx)
              call growr(lmx_,aczlx)
              call growr(lmx_,aczzlx)
              call growr(lmx_,aczrx)
              call growr(lmx_,aczzrx)
            endif
            allx(al,i,j) = llx(l,i,j)
            alrx(al,i,j) = lrx(l,i,j)
          endif
          admlx(al,i,j) = admlx(al,i,j) + dmlx(l,i,j)
          admrx(al,i,j) = admrx(al,i,j) + dmrx(l,i,j)
          aczlx(al,i,j) = aczlx(al,i,j) + czlx(l,i,j)*dmlx(l,i,j)
          aczrx(al,i,j) = aczrx(al,i,j) + czrx(l,i,j)*dmrx(l,i,j)
          aczzlx(al,i,j) = aczzlx(al,i,j) + czzlx(l,i,j)*dmlx(l,i,j)
          aczzrx(al,i,j) = aczzrx(al,i,j) + czzrx(l,i,j)*dmrx(l,i,j)
        enddo
        endif
#endif

      enddo
      enddo
      enddo ! j

      do j=max(2,j_0-1),j_1s
      facj = dt*dxvo(j)/dyvo(j)
      do n=1,nbyzv(j,1)
      do i=i1yzv(n,j,1),i2yzv(n,j,1)

        jl = j
        jr = j+1

        lmy(i,j) = 0

        k1d(:) = k3dy(:,i,j)
        k1d(1) = k1d(2) ! for shallow oceans
        kmax = maxval(k1d)

        krat = k1d / kmax ! normalized K-profile for fractional-in-time
        call get_frt_intervals(xrand,lmo,krat,nlzone,ltb) ! mixing zones
        if(nlzone.eq.0) cycle

        do l = 1,ltb(2,nlzone) ! fill input arrays to maximum depth of any zone
          ml(l) = mo(i,jl,l)
          mr(l) = mo(i,jr,l)
          gl(l) = g3d(l,i,jl)
          gr(l) = g3d(l,i,jr)
          sl(l) = s3d(l,i,jl)
          sr(l) = s3d(l,i,jr)
          pl(l) = p3d(l,i,jl)
          pr(l) = p3d(l,i,jr)
        enddo

        call get_dsoe_info(
     &       'y',
     &       lmm(i,jl),lmm(i,jr),i,j,
     &       ml,mr,gl,gr,sl,sr,pl,pr,
     &       nlzone,ltb,
     &       lmy(i,j),
     &       lly(1,i,j),lry(1,i,j),
     &       dmly(1,i,j),dmry(1,i,j),
     &       czly(1,i,j),czzly(1,i,j),czry(1,i,j),czzry(1,i,j)
     &     )

        ! rescale DSOE mass dmly (kg/m2) -> ((K/dy)*dx*dt)*dmly = (kg)
        do l=1,lmy(i,j)
          dmly(l,i,j) = dmly(l,i,j)*facj*kmax
          dmry(l,i,j) = dmry(l,i,j)*facj*kmax
        enddo

#ifdef TRACERS_OCEAN
        ! accumulate DSOE information from this timestep
        if(ntrtrans.gt.1) then
        lmy_ = size(ally,1)
        do l=1,lmy(i,j)
          found = .false.
          do al=1,almy(i,j)
            if(lly(l,i,j).eq.ally(al,i,j) .and.
     &         lry(l,i,j).eq.alry(al,i,j)) then
              found = .true.
              exit
            endif
          enddo
          if(.not. found) then
            almy(i,j) = almy(i,j) + 1
            al = almy(i,j)
            if(al.gt.12*lmo) then
              !write(6,*) 'for '//cdir//'-direction:'
              write(6,*) 'lm>12*lmo at i,j = ',i,j
              write(6,*) 'increase limit if reasonable'
              call stop_model('lm>12*lmo (tracers y-dir)',255)
            endif
            if(al.gt.lmy_) then ! resize on the fly
              lmy_ = lmy_ + lmo
              call growi(lmy_,ally)
              call growi(lmy_,alry)
              call growr(lmy_,admly)
              call growr(lmy_,admry)
              call growr(lmy_,aczly)
              call growr(lmy_,aczzly)
              call growr(lmy_,aczry)
              call growr(lmy_,aczzry)
            endif
            ally(al,i,j) = lly(l,i,j)
            alry(al,i,j) = lry(l,i,j)
          endif
          admly(al,i,j) = admly(al,i,j) + dmly(l,i,j)
          admry(al,i,j) = admry(al,i,j) + dmry(l,i,j)
          aczly(al,i,j) = aczly(al,i,j) + czly(l,i,j)*dmly(l,i,j)
          aczry(al,i,j) = aczry(al,i,j) + czry(l,i,j)*dmry(l,i,j)
          aczzly(al,i,j) = aczzly(al,i,j) + czzly(l,i,j)*dmly(l,i,j)
          aczzry(al,i,j) = aczzry(al,i,j) + czzry(l,i,j)*dmry(l,i,j)
        enddo
        endif
#endif

      enddo
      enddo
      enddo

      call get_tdmix_mrelxyz

      return

#ifdef TRACERS_OCEAN

      contains

      subroutine growi(newlm,arr)
      ! reallocate arr to size newlm while preserving prior contents
      implicit none
      integer :: newlm
      integer, dimension(:,:,:), allocatable :: arr
!
      integer :: j_0h,j_1h,oldlm
      integer, dimension(:,:,:), allocatable :: tmp

      call getdomainbounds(grid, j_strt_halo=j_0h, j_stop_halo=j_1h)

      oldlm = size(arr,1)

      allocate(tmp(newlm,im,j_0h:j_1h))
      tmp = 0
      tmp(1:oldlm,:,:) = arr(:,:,:)
      call move_alloc(tmp,arr)

      end subroutine growi

      subroutine growr(newlm,arr)
      ! reallocate arr to size newlm while preserving prior contents
      implicit none
      integer :: newlm
      real*8, dimension(:,:,:), allocatable :: arr
!
      integer :: j_0h,j_1h,oldlm
      real*8, dimension(:,:,:), allocatable :: tmp

      call getdomainbounds(grid, j_strt_halo=j_0h, j_stop_halo=j_1h)

      oldlm = size(arr,1)

      allocate(tmp(newlm,im,j_0h:j_1h))
      tmp = 0.
      tmp(1:oldlm,:,:) = arr(:,:,:)
      call move_alloc(tmp,arr)

      end subroutine growr

#endif

      end subroutine tdmix_prep

      subroutine get_tdmix_mrelxyz
      use tdmix_mod
      use ocean, only : im,jm,lmo
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds
      implicit none
c
      integer :: i,j,l,n,ll,lr,il,ir

      integer :: j_0,j_1,j_0s,j_1s
      logical :: have_north_pole

      real*8 :: byim,mxfer
      real*8, dimension(lmo) :: mxferx,mxfery,mxferz

      call getdomainbounds(grid, j_strt=j_0, j_stop=j_1,
     &               j_strt_skp=j_0s, j_stop_skp=j_1s,
     &               have_north_pole=have_north_pole)


      ! Weights for relaxing prognostic subgrid tracer profiles
      ! to interpolation-determined "smooth" profiles are equal
      ! to the exchanged water mass divided by the gridbox mass.
      ! The weight in the vertical direction only incorporates
      ! the layer-jumping component of horizontal mass exchange.
      do j=j_0s,j_1s
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        do l=1,lmm(i,j)
          mxferx(l) = 0.
          mxfery(l) = 0.
          mxferz(l) = 0.
        enddo
        if(i.eq.1) then
          il = im
        else
          il = i-1
        endif
        do l=1,lmx(il,j)
          ll = llx(l,il,j)
          lr = lrx(l,il,j)
          mxfer = .5d0*(dmlx(l,il,j)+dmrx(l,il,j))
          mxferx(lr) = mxferx(lr) + mxfer
          if(ll.ne.lr) then ! layer jump
            mxferz(lr) = mxferz(lr) + mxfer
          endif
        enddo
        if(i.eq.im) then
          ir = 1
        else
          ir = i+1
        endif
        do l=1,lmx(i,j)
          ll = llx(l,i,j)
          lr = lrx(l,i,j)
          mxfer = .5d0*(dmlx(l,i,j)+dmrx(l,i,j))
          mxferx(ll) = mxferx(ll) + mxfer
          if(ll.ne.lr) then ! layer jump
            mxferz(ll) = mxferz(ll) + mxfer
          endif
        enddo
        do l=1,lmy(i,j-1)
          ll = lly(l,i,j-1)
          lr = lry(l,i,j-1)
          mxfer = .5d0*(dmly(l,i,j-1)+dmry(l,i,j-1))
          mxfery(lr) = mxfery(lr) + mxfer
          if(ll.ne.lr) then ! layer jump
            mxferz(lr) = mxferz(lr) + mxfer
          endif
        enddo
        do l=1,lmy(i,j)
          ll = lly(l,i,j)
          lr = lry(l,i,j)
          mxfer = .5d0*(dmly(l,i,j)+dmry(l,i,j))
          mxfery(ll) = mxfery(ll) + mxfer
          if(ll.ne.lr) then ! layer jump
            mxferz(ll) = mxferz(ll) + mxfer
          endif
        enddo
        do l=1,lmm(i,j)
          mrelx(i,j,l) = mxferx(l)/mokgsv(i,j,l)
          mrely(i,j,l) = mxfery(l)/mokgsv(i,j,l)
          mrelz(i,j,l) = mxferz(l)/mokgsv(i,j,l)
        enddo
      enddo
      enddo
      enddo

      if(have_north_pole) then
        byim = 1d0/real(im,kind=8)
        j = jm
        do l=1,lmm(1,j)
          mxfery(l) = 0.
          mxferz(l) = 0.
        enddo
        do i=1,im
          do l=1,lmy(i,j-1)
            ll = lly(l,i,j-1)
            lr = lry(l,i,j-1)
            mxfer = .5d0*(dmly(l,i,j-1)+dmry(l,i,j-1))
            mxfery(lr) = mxfery(lr) + mxfer
            if(ll.ne.lr) then ! layer jump
              mxferz(lr) = mxferz(lr) + mxfer
            endif
          enddo
        enddo
        i = 1
        do l=1,lmm(i,j)
          mrely(i,j,l) = byim*mxfery(l)/mokgsv(i,j,l)
          mrelz(i,j,l) = byim*mxferz(l)/mokgsv(i,j,l)
        enddo
      endif

      mrelx = min(mrelx,1d0)
      mrely = min(mrely,1d0)
      mrelz = min(mrelz,1d0)

      return
      end subroutine get_tdmix_mrelxyz

      subroutine tdmix(trm,qlimit,fl3d
#ifdef TDMIX_AUX_DIAGS
     &     ,fl3ds
#endif
     &     )
!@sum tdmix applies mass exchange rates to the tracer field
      use tdmix_mod
      use ocean, only : im,jm,lmo
      use ocean, only : dzo
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use ocean, only : nbyzu,i1yzu,i2yzu,lmu
      use ocean, only : nbyzv,i1yzv,i2yzv,lmv
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds
      use domain_decomp_1d, only : halo_update,halo_update_column,
     &     south
cons      use domain_decomp_1d, only : globalsum
      implicit none
!@var trm input/output tracer mass (kg)
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     trm
!@var qlimit whether to impose non-negativity for trm
      logical :: qlimit
!@var fl3d diagnosed x-y-p fluxes (kg/s)
!@+   Slantwise fluxes are projected into the model x-y-p coordinate
!@+   system by constructing horizontal and vertical fluxes whose
!@+   convergence equals the convergence of the slantwise fluxes for
!@+   each model gridcell.  This is performed separately for each DSOE.
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo,3) ::
     &     fl3d
#ifdef TDMIX_AUX_DIAGS
!@var fl3ds the symmetric part of fl3d
     &    ,fl3ds
#endif
c
      real*8, dimension(lmo,im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     tr0,tz0,tzz0
      real*8 :: facl,facr,byim,tl,tr
      integer :: i,j,l,n,ll,lr,il,ir
      integer :: j_0,j_1,j_0s,j_1s
      logical :: have_north_pole

      real*8, dimension(lmo) :: fl1d,fzxl1d,fzxr1d,fzyl1d,fzyr1d
      integer :: l1,l2
      real*8 :: dmnet
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     fzxl,fzxr,fzyl,fzyr
#ifdef TDMIX_AUX_DIAGS
     &    ,fzxls,fzxrs,fzyls,fzyrs
#endif

cons      real*8, dimension(grid%j_strt_halo:grid%j_stop_halo) ::
cons     &     psumbef,psumaft
cons      real*8 :: sumbef,sumaft
cons      psumbef = 0.; psumaft = 0.

      call getdomainbounds(grid, j_strt=j_0, j_stop=j_1,
     &               j_strt_skp=j_0s, j_stop_skp=j_1s,
     &               have_north_pole=have_north_pole)


      do j=j_0,j_1
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
        if(i.eq.1) then
          il = im
        else
          il = i-1
        endif
        do l=1,lmm(i,j)
!          tr0(l,i,j) = trm(i,j,l)/(mo(i,j,l)*dxypo(j))
          tr0(l,i,j) = trm(i,j,l)/mokgsv(i,j,l)
        enddo
        call get_zmoms(tr0(1,i,j),tz0(1,i,j),tzz0(1,i,j),lmm(i,j))
        if(qlimit) then
          call adj_zmoms(tr0(1,i,j),tz0(1,i,j),tzz0(1,i,j),
     &         lrx(1,il,j),czrx(1,il,j),czzrx(1,il,j),
     &         lmm(i,j),lmx(il,j))
          call adj_zmoms(tr0(1,i,j),tz0(1,i,j),tzz0(1,i,j),
     &         llx(1,i,j),czlx(1,i,j),czzlx(1,i,j),
     &         lmm(i,j),lmx(i ,j))
          call adj_zmoms(tr0(1,i,j),tz0(1,i,j),tzz0(1,i,j),
     &         lry(1,i,j-1),czry(1,i,j-1),czzry(1,i,j-1),
     &         lmm(i,j),lmy(i,j-1))
          call adj_zmoms(tr0(1,i,j),tz0(1,i,j),tzz0(1,i,j),
     &         lly(1,i,j),czly(1,i,j),czzly(1,i,j),
     &         lmm(i,j),lmy(i ,j))
        endif
      enddo
      enddo
      enddo
      if(have_north_pole) then
        do i=2,im
          tr0(:,i,jm) = tr0(:,1,jm)
          tz0(:,i,jm) = tz0(:,1,jm)
          tzz0(:,i,jm) = tzz0(:,1,jm)
        enddo
      endif

      call halo_update_column(grid,tr0)
      call halo_update_column(grid,tz0)
      call halo_update_column(grid,tzz0)

      fl3d = 0.
      fzxl = 0.
      fzxr = 0.
      fzyl = 0.
      fzyr = 0.
#ifdef TDMIX_AUX_DIAGS
      fl3ds = 0.
      fzxls = 0.
      fzxrs = 0.
      fzyls = 0.
      fzyrs = 0.
#endif

      do j=j_0s,j_1s
      do n=1,nbyzm(j,1)
      do i=i1yzm(n,j,1),i2yzm(n,j,1)
cons        psumbef(j) = psumbef(j) + sum(trm(i,j,1:lmm(i,j)))

        ! x-direction
        if(i.eq.1) then
          il = im
        else
          il = i-1
        endif
        do l=1,lmx(il,j)
          ll = llx(l,il,j)
          lr = lrx(l,il,j)
          facl = dmlx(l,il,j)
          facr = dmrx(l,il,j)
          tl = tr0(ll,il,j)
     &       +(tz0(ll,il,j)* czlx(l,il,j)
     &       +tzz0(ll,il,j)*czzlx(l,il,j))
          tr = tr0(lr,i  ,j)
     &       +(tz0(lr,i  ,j)* czrx(l,il,j)
     &       +tzz0(lr,i  ,j)*czzrx(l,il,j))
          trm(i,j,lr) = trm(i,j,lr) + (tl*facl-tr*facr)
        enddo
        if(i.eq.im) then
          ir = 1
        else
          ir = i+1
        endif
        do l=1,lmx(i,j)
          ll = llx(l,i,j)
          lr = lrx(l,i,j)
          facl = dmlx(l,i,j)
          facr = dmrx(l,i,j)
          tl = tr0(ll,i ,j)
     &       +(tz0(ll,i ,j)* czlx(l,i,j)
     &       +tzz0(ll,i ,j)*czzlx(l,i,j))
          tr = tr0(lr,ir,j)
     &       +(tz0(lr,ir,j)* czrx(l,i,j)
     &       +tzz0(lr,ir,j)*czzrx(l,i,j))
          trm(i,j,ll) = trm(i,j,ll) + (tr*facr-tl*facl)

          ! Flux diagnostics.  See above remarks regarding slantwise -> x-y-p
          dmnet = -(tr*facr-tl*facl)
          call iso_to_xyp_diags(
     &         llx(l,i,j),lrx(l,i,j),lmu(i,j),dmnet,
     &         l1,l2,fl1d,fzxl1d,fzxr1d)
          fl3d(i,j,l1:l2,1) = fl3d(i,j,l1:l2,1) + fl1d(l1:l2)
          if(l2.gt.l1) then
            fzxl(i,j,l1:l2-1) = fzxl(i,j,l1:l2-1) + fzxl1d(l1:l2-1)
            fzxr(i,j,l1:l2-1) = fzxr(i,j,l1:l2-1) + fzxr1d(l1:l2-1)
          endif

#ifdef TDMIX_AUX_DIAGS
          dmnet = -(tr-tl)*min(facr,facl)
          call iso_to_xyp_diags(
     &         llx(l,i,j),lrx(l,i,j),lmu(i,j),dmnet,
     &         l1,l2,fl1d,fzxl1d,fzxr1d)
          fl3ds(i,j,l1:l2,1) = fl3ds(i,j,l1:l2,1) + fl1d(l1:l2)
          if(l2.gt.l1) then
            fzxls(i,j,l1:l2-1) = fzxls(i,j,l1:l2-1) + fzxl1d(l1:l2-1)
            fzxrs(i,j,l1:l2-1) = fzxrs(i,j,l1:l2-1) + fzxr1d(l1:l2-1)
          endif
#endif

        enddo

        ! y-direction
        do l=1,lmy(i,j-1)
          ll = lly(l,i,j-1)
          lr = lry(l,i,j-1)
          facl = dmly(l,i,j-1)
          facr = dmry(l,i,j-1)
          tl = tr0(ll,i,j-1)
     &       +(tz0(ll,i,j-1)* czly(l,i,j-1)
     &       +tzz0(ll,i,j-1)*czzly(l,i,j-1))
          tr = tr0(lr,i,j  )
     &       +(tz0(lr,i,j  )* czry(l,i,j-1)
     &       +tzz0(lr,i,j  )*czzry(l,i,j-1))
          trm(i,j,lr) = trm(i,j,lr) + (tl*facl-tr*facr)
        enddo
        do l=1,lmy(i,j)
          ll = lly(l,i,j)
          lr = lry(l,i,j)
          facl = dmly(l,i,j)
          facr = dmry(l,i,j)
          tl = tr0(ll,i,j  )
     &       +(tz0(ll,i,j  )* czly(l,i,j)
     &       +tzz0(ll,i,j  )*czzly(l,i,j))
          tr = tr0(lr,i,j+1)
     &       +(tz0(lr,i,j+1)* czry(l,i,j)
     &       +tzz0(lr,i,j+1)*czzry(l,i,j))
          trm(i,j,ll) = trm(i,j,ll) + (tr*facr-tl*facl)

          ! Flux diagnostics.  See above remarks regarding slantwise -> x-y-p
          dmnet = -(tr*facr-tl*facl)
          call iso_to_xyp_diags(
     &         lly(l,i,j),lry(l,i,j),lmv(i,j),dmnet,
     &         l1,l2,fl1d,fzyl1d,fzyr1d)
          fl3d(i,j,l1:l2,2) = fl3d(i,j,l1:l2,2) + fl1d(l1:l2)
          if(l2.gt.l1) then
            fzyl(i,j,l1:l2-1) = fzyl(i,j,l1:l2-1) + fzyl1d(l1:l2-1)
            fzyr(i,j,l1:l2-1) = fzyr(i,j,l1:l2-1) + fzyr1d(l1:l2-1)
          endif
#ifdef TDMIX_AUX_DIAGS
          dmnet = -(tr-tl)*min(facr,facl)
          call iso_to_xyp_diags(
     &         lly(l,i,j),lry(l,i,j),lmv(i,j),dmnet,
     &         l1,l2,fl1d,fzyl1d,fzyr1d)
          fl3ds(i,j,l1:l2,2) = fl3ds(i,j,l1:l2,2) + fl1d(l1:l2)
          if(l2.gt.l1) then
            fzyls(i,j,l1:l2-1) = fzyls(i,j,l1:l2-1) + fzyl1d(l1:l2-1)
            fzyrs(i,j,l1:l2-1) = fzyrs(i,j,l1:l2-1) + fzyr1d(l1:l2-1)
          endif
#endif

        enddo
cons        psumaft(j) = psumaft(j) + sum(trm(i,j,1:lmm(i,j)))
      enddo
      enddo
      enddo

      if(have_north_pole) then
        byim = 1d0/real(im,kind=8)
        j = jm
cons        psumbef(j) = im*sum(trm(1,j,1:lmm(1,j)))
        do i=1,im
          do l=1,lmy(i,j-1)
            ll = lly(l,i,j-1)
            lr = lry(l,i,j-1)
            facl = dmly(l,i,j-1)
            facr = dmry(l,i,j-1)
            tl = tr0(ll,i,j-1)
     &         +(tz0(ll,i,j-1)* czly(l,i,j-1)
     &         +tzz0(ll,i,j-1)*czzly(l,i,j-1))
            tr = tr0(lr,1,j  )
     &         +(tz0(lr,1,j  )* czry(l,i,j-1)
     &         +tzz0(lr,1,j  )*czzry(l,i,j-1))
            trm(1,j,lr) = trm(1,j,lr) + (tl*facl-tr*facr)*byim
          enddo
        enddo
cons        psumaft(j) = im*sum(trm(1,j,1:lmm(1,j)))
      endif

cons      call globalsum(grid,psumbef,sumbef)
cons      call globalsum(grid,psumaft,sumaft)
cons      if(grid%gid.eq.0) then
cons        write(6,*) 'sumrat ',sumaft/sumbef
cons      endif


      ! Combine x- and y-direction contributions
      ! to z-direction flux diagnostics
      call halo_update(grid,fzyr,from=south)
#ifdef TDMIX_AUX_DIAGS
      call halo_update(grid,fzyrs,from=south)
#endif
      do l=1,lmo-1
        do j=j_0s,j_1s
        i=1
        if(l.lt.lmm(i,j)) then
          fl3d(i,j,l,3) =
     &         fzxr(im,j,l)+fzxl(i,j,l)+fzyr(i,j-1,l)+fzyl(i,j,l)
#ifdef TDMIX_AUX_DIAGS
          fl3ds(i,j,l,3) =
     &         fzxrs(im,j,l)+fzxls(i,j,l)+fzyrs(i,j-1,l)+fzyls(i,j,l)
#endif
        endif
        do n=1,nbyzm(j,l+1)
        do i=max(2,i1yzm(n,j,l+1)),i2yzm(n,j,l+1)
          fl3d(i,j,l,3) =
     &         fzxr(i-1,j,l)+fzxl(i,j,l)+fzyr(i,j-1,l)+fzyl(i,j,l)
#ifdef TDMIX_AUX_DIAGS
          fl3ds(i,j,l,3) =
     &         fzxrs(i-1,j,l)+fzxls(i,j,l)+fzyrs(i,j-1,l)+fzyls(i,j,l)
#endif
        enddo
        enddo
        enddo
      enddo

      if(have_north_pole) then
        byim = 1d0/real(im,kind=8)
        j = jm
        do l=1,lmm(1,j)-1
          fl3d(:,j,l,3) = sum(fzyr(:,j-1,l))*byim
#ifdef TDMIX_AUX_DIAGS
          fl3ds(:,j,l,3) = sum(fzyrs(:,j-1,l))*byim
#endif
        enddo
      endif

      return
      end subroutine tdmix

      subroutine iso_to_xyp_diags(
     &     llx,lrx,lmu,dmnet,l1,l2,fl,fzxl,fzxr
     &     )
      use ocean, only : lmo
      use ocean, only : dzo
      implicit none
      integer :: llx,lrx,lmu
      real*8 :: dmnet
      integer :: l1,l2
      real*8, dimension(lmo) :: fl,fzxl,fzxr
c
      integer :: lll
      real*8 :: fzxl_,fzxr_

      l1 = min(llx,lrx)
      l2 = max(llx,lrx)
      if(l2.eq.l1) then
        fl(l1) = dmnet
      else ! jumping layers
        lll = l2
        if(l2.gt.lmu) then
          fl(lmu+1:l2) = 0.
          lll = lmu
        endif
            !fl(l1:lll) = dmnet/real(lll-l1+1,kind=8)
        fl(l1:lll) = dmnet*(dzo(l1:lll)/sum(dzo(l1:lll)))
        if(lrx.lt.llx) then
          fzxl_ = -sum(fl(l1:l2-1))
          fzxr_ = -fl(l2)
        else
          fzxr_ = sum(fl(l1:l2-1))
          fzxl_ = fl(l2)
        endif
        do lll=l2-1,l1,-1
          fzxl(lll) = fzxl_
          fzxr(lll) = fzxr_
          fzxl_ = fzxl_ + fl(lll)
          fzxr_ = fzxr_ - fl(lll)
        enddo
      endif
      end subroutine iso_to_xyp_diags

      subroutine relax_qusmoms(mnew,trm,
     &     trxm,trym,trzm, trxxm,tryym,trzzm, trxym,tryzm,trzxm
     &     )
!@sum relax_qusmoms relaxes prognostic subgrid profiles toward
!@+   diagnosed subgrid profiles at rates described/computed in
!@+   tdmix_prep
      use ocean, only : im,jm,lmo
      use ocean, only : nbyzm,i1yzm,i2yzm,lmm
      use ocean, only : nbyzu,i1yzu,i2yzu,lmu
      use ocean, only : nbyzv,i1yzv,i2yzv,lmv
      use ocean, only : dzo_orig=>dzo, bydzoe
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds
      use domain_decomp_1d, only : halo_update,halo_update_column,
     &     south,north
      use tdmix_mod, only : mokgsv,mrelx,mrely,mrelz
      implicit none
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     trm,mnew,
     &     trxm,trym,trzm, trxxm,tryym,trzzm, trxym,tryzm,trzxm
c
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,lmo) ::
     &     tr0,try,mrat
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo,0:lmo) ::
     &     trz
c
      integer :: i,j,l,n
      integer :: j_0,j_1,j_0s,j_1s
      logical :: have_north_pole
c
      real*8, dimension(0:im+1) :: tr01d,trx1d
      real*8 :: frel,trxnew,trxxnew,trynew,tryynew,trznew,trzznew

      real*8, dimension(0:lmo+1) :: dzo
      real*8, dimension(lmo) :: by6arr,dzby12

      dzo(1:lmo) = dzo_orig(1:lmo)
      dzo(0) = dzo(1)
      dzo(lmo+1) = dzo(lmo)
      do l=1,lmo
        by6arr(l) = .5d0*dzo(l)/sum(dzo(l-1:l+1))
        dzby12(l) = .5d0*by6arr(l)*dzo(l)
      enddo

      call getdomainbounds(grid, j_strt=j_0, j_stop=j_1,
     &               j_strt_skp=j_0s, j_stop_skp=j_1s,
     &               have_north_pole=have_north_pole)

      mrat = 0.
      do l=1,lmo
      do j=j_0,j_1
      do n=1,nbyzm(j,l)
      do i=i1yzm(n,j,l),i2yzm(n,j,l)
          mrat(i,j,l) = min(1.1d0,mnew(i,j,l)/mokgsv(i,j,l))
          tr0(i,j,l) = trm(i,j,l)/mokgsv(i,j,l)
      enddo
      enddo
      enddo
      enddo

      !if(maxval(mrat).gt.1.01d0) write(6,*) 'mrat ',maxval(mrat)

      if(have_north_pole) then
        j = jm
        do l=1,lmm(1,j)
          tr0(2:im,j,l) = tr0(1,j,l)
        enddo
      endif

      call halo_update(grid,tr0,from=north)

      try = 0.
      do l=1,lmo
        do j=j_0s,j_1s
        do n=1,nbyzv(j,l)
        do i=i1yzv(n,j,l),i2yzv(n,j,l)
          try(i,j,l) = tr0(i,j+1,l)-tr0(i,j,l)
        enddo
        enddo
        enddo
        do j=j_0s,j_1s
          do n=1,nbyzm(j,l)
          do i=i1yzm(n,j,l),i2yzm(n,j,l)
            tr01d(i) = tr0(i,j,l)
          enddo
          enddo
          tr01d(0)    = tr01d(im)
          tr01d(im+1) = tr01d(1)
          trx1d(im) = 0.
          do n=1,nbyzu(j,l)
          do i=i1yzu(n,j,l),i2yzu(n,j,l)
            trx1d(i) = tr01d(i+1)-tr01d(i)
          enddo
          enddo
          trx1d(0) = trx1d(im)
          do n=1,nbyzm(j,l)
            i=i1yzm(n,j,l)
            if(i.gt.1) trx1d(i-1) = 0.
            i=i2yzm(n,j,l)
            if(i.lt.im) trx1d(i) = 0.
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              frel = mrelx(i,j,l)
              trxnew = .25d0*(trx1d(i-1)+trx1d(i))
              trxnew = (1d0-frel)*trxm(i,j,l)*mrat(i,j,l)
     &             + frel*trxnew*mnew(i,j,l)
              trxm(i,j,l) = trxnew
            enddo
            if(i1yzm(n,j,l).ne.i2yzm(n,j,l)) then
              i=i1yzm(n,j,l)
              if(i.gt.1 .or.
     &        (i.eq.1 .and. l.gt.lmm(im,j))
     &           ) then
                trx1d(i-1) = trx1d(i)
              endif
              i=i2yzm(n,j,l)
              if(i.lt.im .or.
     &        (i.eq.im .and. l.gt.lmm(1,j))
     &           ) then
                trx1d(i) = trx1d(i-1)
              endif
            endif
            do i=i1yzm(n,j,l),i2yzm(n,j,l)
              frel = mrelx(i,j,l)
              trxxnew = (trx1d(i)-trx1d(i-1))/12d0
              trxxnew = (1d0-frel)*trxxm(i,j,l)*mrat(i,j,l)
     &             + frel*trxxnew*mnew(i,j,l)
              trxxm(i,j,l) = trxxnew
            enddo
          enddo
        enddo
      enddo

      call halo_update(grid,try,from=south)

      do l=1,lmo
        do j=j_0s,j_1s
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          frel = mrely(i,j,l)
          trynew = .25d0*(try(i,j-1,l)+try(i,j,l))
          trynew = (1d0-frel)*trym(i,j,l)*mrat(i,j,l)
     &         + frel*trynew*mnew(i,j,l)
          trym(i,j,l) = trynew
          if(l.le.minval(lmm(i,j-1:j+1))) then
            tryynew = (try(i,j,l)-try(i,j-1,l))/12d0
          else
            tryynew = 0.
          endif
          tryynew = (1d0-frel)*tryym(i,j,l)*mrat(i,j,l)
     &         + frel*tryynew*mnew(i,j,l)
          tryym(i,j,l) = tryynew
        enddo
        enddo
        enddo
      enddo

      do l=1,lmo-1
      do j=j_0,j_1
        do n=1,nbyzm(j,l+1)
          do i=i1yzm(n,j,l+1),i2yzm(n,j,l+1)
            trz(i,j,l) = (tr0(i,j,l+1)-tr0(i,j,l))*bydzoe(l)
          enddo
        enddo
      enddo
      enddo
      do j=j_0,j_1
        do n=1,nbyzm(j,1)
          do i=i1yzm(n,j,1),i2yzm(n,j,1)
            trz(i,j,0) = trz(i,j,1)
            l = lmm(i,j)
            trz(i,j,l) = trz(i,j,l-1)
          enddo
        enddo
      enddo

      do l=1,lmo
      do j=j_0,j_1
        do n=1,nbyzm(j,l)
          do i=i1yzm(n,j,l),i2yzm(n,j,l)
            frel = mrelz(i,j,l)
            trznew = by6arr(l)*(
     &           dzo(l+1)*trz(i,j,l-1)+dzo(l-1)*trz(i,j,l)
     &           + .5d0*dzo(l)*(trz(i,j,l-1)+trz(i,j,l))  )
            trznew = (1d0-frel)*trzm(i,j,l)*mrat(i,j,l)
     &           + frel*trznew*mnew(i,j,l)
            trzm(i,j,l) = trznew
            trzznew = dzby12(l)*(trz(i,j,l)-trz(i,j,l-1))
            trzznew = (1d0-frel)*trzzm(i,j,l)*mrat(i,j,l)
     &           + frel*trzznew*mnew(i,j,l)
            trzzm(i,j,l) = trzznew
          enddo
        enddo
      enddo
      enddo

      do l=1,lmo
        do j=j_0,j_1
        do n=1,nbyzm(j,l)
        do i=i1yzm(n,j,l),i2yzm(n,j,l)
          trxym(i,j,l) = trxym(i,j,l)*mrat(i,j,l)
          tryzm(i,j,l) = tryzm(i,j,l)*mrat(i,j,l)
          trzxm(i,j,l) = trzxm(i,j,l)*mrat(i,j,l)
        enddo
        enddo
        enddo
      enddo

      end subroutine relax_qusmoms

#ifdef TRACERS_OCEAN
      subroutine tdmix_longstep_prep
      ! Transfer DSOE accumulation arrays a[...] to the [...] used by
      ! transport routines. Instead of copying the a[...] to [...],
      ! move_alloc is used, which deallocates the a[...].
      ! Both [...] and a[...] will be re-allocated the next cycle,
      ! the latter back to the default small size.
      use tdmix_mod
      use ocean, only : nbyzm,i1yzm,i2yzm
      use ocean, only : nbyzu,i1yzu,i2yzu
      use ocean, only : nbyzv,i1yzv,i2yzv
      use ocean, only : ntrtrans,motr,dxypo
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : getdomainbounds,globalmin
      use domain_decomp_1d, only : globalmax ! print stats
      use constant, only : teeny
      implicit none
      integer :: i,j,l,n
      integer :: j_0,j_1,j_0s,j_1s
      integer :: lmloc,lmmin,lmmax

      if(ntrtrans.le.1) return ! nothing to do

      call getdomainbounds(grid, j_strt=j_0, j_stop=j_1,
     &     j_strt_skp=j_0s, j_stop_skp=j_1s)


!      lmloc = size(allx,1)
!      call globalmax(grid,lmloc,lmmax)
!      if(all(nbyzu(j_0:j_1,1).eq.0)) lmloc = 123456
!      call globalmin(grid,lmloc,lmmin)
!      if(grid%gid.eq.0) then
!        write(6,*) 'minlmx,maxlmx ',lmmin,lmmax
!      endif
!      lmloc = size(ally,1)
!      call globalmax(grid,lmloc,lmmax)
!      if(all(nbyzu(j_0:j_1,1).eq.0)) lmloc = 123456
!      call globalmin(grid,lmloc,lmmin)
!      if(grid%gid.eq.0) then
!        write(6,*) 'minlmy,maxlmy ',lmmin,lmmax
!      endif

      call move_alloc(almx,lmx)
      call move_alloc(allx,llx)
      call move_alloc(alrx,lrx)
      call move_alloc(admlx,dmlx)
      call move_alloc(admrx,dmrx)
      call move_alloc(aczlx,czlx)
      call move_alloc(aczzlx,czzlx)
      call move_alloc(aczrx,czrx)
      call move_alloc(aczzrx,czzrx)

      do j=j_0s,j_1s
      do n=1,nbyzu(j,1)
      do i=i1yzu(n,j,1),i2yzu(n,j,1)
        do l=1,lmx(i,j)
          czlx(l,i,j) = czlx(l,i,j)/(dmlx(l,i,j)+teeny)
          czzlx(l,i,j) = czzlx(l,i,j)/(dmlx(l,i,j)+teeny)
          czrx(l,i,j) = czrx(l,i,j)/(dmrx(l,i,j)+teeny)
          czzrx(l,i,j) = czzrx(l,i,j)/(dmrx(l,i,j)+teeny)
        enddo
      enddo
      enddo
      enddo

      call move_alloc(almy,lmy)
      call move_alloc(ally,lly)
      call move_alloc(alry,lry)
      call move_alloc(admly,dmly)
      call move_alloc(admry,dmry)
      call move_alloc(aczly,czly)
      call move_alloc(aczzly,czzly)
      call move_alloc(aczry,czry)
      call move_alloc(aczzry,czzry)

      do j=max(2,j_0-1),j_1s
      do n=1,nbyzv(j,1)
      do i=i1yzv(n,j,1),i2yzv(n,j,1)
        do l=1,lmy(i,j)
          czly(l,i,j) = czly(l,i,j)/(dmly(l,i,j)+teeny)
          czzly(l,i,j) = czzly(l,i,j)/(dmly(l,i,j)+teeny)
          czry(l,i,j) = czry(l,i,j)/(dmry(l,i,j)+teeny)
          czzry(l,i,j) = czzry(l,i,j)/(dmry(l,i,j)+teeny)
        enddo
      enddo
      enddo
      enddo

      do l=1,lmo
      do j=j_0,j_1
      do n=1,nbyzm(j,l)
      do i=i1yzm(n,j,l),i2yzm(n,j,l)
        mokgsv(i,j,l) = motr(i,j,l)*dxypo(j)
      enddo
      enddo
      enddo
      enddo

      call get_tdmix_mrelxyz

      end subroutine tdmix_longstep_prep

      subroutine tdmix_longstep_alloc
      use tdmix_mod
      use domain_decomp_1d, only : getdomainbounds
      use ocean, only : im,lmo
      use oceanr_dim, only : grid=>ogrid
      implicit none
      integer :: j_0h,j_1h

      call getdomainbounds(grid, j_strt_halo=j_0h, j_stop_halo=j_1h)

      allocate(
     &         almx(im,j_0h:j_1h),
     &         almy(im,j_0h:j_1h),
     &         allx(almaxoverlap,im,j_0h:j_1h),
     &         alrx(almaxoverlap,im,j_0h:j_1h),
     &         ally(almaxoverlap,im,j_0h:j_1h),
     &         alry(almaxoverlap,im,j_0h:j_1h),
c
     &         admlx(almaxoverlap,im,j_0h:j_1h),
     &         admrx(almaxoverlap,im,j_0h:j_1h),
     &         admly(almaxoverlap,im,j_0h:j_1h),
     &         admry(almaxoverlap,im,j_0h:j_1h),
c
     &          aczlx(almaxoverlap,im,j_0h:j_1h),
     &         aczzlx(almaxoverlap,im,j_0h:j_1h),
     &          aczrx(almaxoverlap,im,j_0h:j_1h),
     &         aczzrx(almaxoverlap,im,j_0h:j_1h),
     &          aczly(almaxoverlap,im,j_0h:j_1h),
     &         aczzly(almaxoverlap,im,j_0h:j_1h),
     &          aczry(almaxoverlap,im,j_0h:j_1h),
     &         aczzry(almaxoverlap,im,j_0h:j_1h)
     &        )


      almx = 0; almy = 0; allx = 0; alrx = 0; ally = 0; alry = 0
      admlx = 0.; admrx = 0.; admly = 0.; admry = 0.

      aczlx = 0.; aczzlx = 0.; aczrx = 0.; aczzrx = 0. 
      aczly = 0.; aczzly = 0.; aczry = 0.; aczzry = 0.

      end subroutine tdmix_longstep_alloc

      subroutine tdmix_longstep_finish
      use tdmix_mod
      use ocean, only : ntrtrans
      implicit none

      if(ntrtrans.le.1) return ! nothing to do

      call realloc_tdmix_regstep
      call tdmix_longstep_alloc

      end subroutine tdmix_longstep_finish

#endif
