#include "hycom_mpi_hacks.h"
      subroutine tradv0(n,nn)
c
c --- online tracer advection routine designed for intermittent
c --- (i.e., long time step) execution. there are 3 entries:
c
c --- tradv0 - initializes mass flux arrays and saves initial -dp- field
c ---          (should be called immediately  a f t e r  diapfl)
c --- tradv1 - builds up time integral of horizontal mass fluxes
c --- tradv2 - performs the actual transport operation
c ---          (should be called immediately  b e f o r e  diapfl)
c
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : oddev
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      implicit none
      integer i,j,k,l,n,nn
c
c --- initialize arrays
c
      do 1 j=J_0, J_1
      do 1 k=1,kk
c
      do 2 l=1,isu(j)
      do 2 i=ifu(j,l),ilu(j,l)
 2    ufxcum(i,j,k)=0.
c
      do 3 l=1,isv(j)
      do 3 i=ifv(j,l),ilv(j,l)
 3    vfxcum(i,j,k)=0.
c
      do 4 l=1,isp(j)
      do 4 i=ifp(j,l),ilp(j,l)
 4    dpinit(i,j,k)=dp(i,j,k+nn)
 1    continue
      oddev=n
      if(AM_I_ROOT() )
     &  write (*,'(a)') 'tracer transport arrays initialized'
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine tradv1(n,nn)
c
c --- build up time integrals of horiz. mass fluxes
c
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : oddev,delt1
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      implicit none
      integer i,j,k,l,n,nn
c
      if (n.ne.oddev) then
        write (*,'(2(a,i2/))')
     .   'tracer advection time interval begins/ends at n =',oddev,
     .    '=> mass fluxes must be accumulated when n =',oddev
        stop '(n=oddev required in tradv1)'
      end if
c
      do 5 j=J_0, J_1
      do 5 k=1,kk
c
      do 6 l=1,isu(j)
      do 6 i=ifu(j,l),ilu(j,l)
 6    ufxcum(i,j,k)=ufxcum(i,j,k)+uflx(i,j,k)*delt1

      do 7 l=1,isv(j)
      do 7 i=ifv(j,l),ilv(j,l)
 7    vfxcum(i,j,k)=vfxcum(i,j,k)+vflx(i,j,k)*delt1
 5    continue
      if(AM_I_ROOT() )
     &  write (*,'(a)') 'mass fluxes saved for tracer transport'
      return
      end
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      subroutine tradv2(n,nn)
c
c --- advect tracer over 'mixfrq' time steps
c
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : oddev
      USE HYCOM_ARRAYS, only : ufxcum_loc => ufxcum,
     &     vfxcum_loc => vfxcum, p_loc => p, dp_loc => dp,
     &     scp2_loc => scp2,
     &     scp2i_loc => scp2i, util1_loc => util1, util2_loc => util2,
     &     tracer_loc => tracer, dpinit_loc => dpinit
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT, HALO_UPDATE, NORTH
      implicit none
      integer i,j,k,l,n,nn,ib,jb
c
      real vertfx(idm,J_0H:J_1H,kdm),hordiv(idm,J_0H:J_1H,kdm),
     .     coldiv(idm),verdiv,q,fluxdv,thkchg
      real vertfx_glob(idm,jdm,kdm)

!          trcold(jdm,kdm),prold(jdm,kdm+1),
!          trcnew(jdm,kdm),prnew(jdm,kdm+1)
      integer ka,nt
      character string*18
c
      if (n.ne.oddev) then
        write (*,'(2(a,i2))')
     .   'tracer advection interval began at n =',n,
     .    '  and must end at n=',oddev
        stop '(n=oddev required in tradv2)'
      end if
c
c --- compute horizontal flux divergence (units: per transport time step)
c
      do k=1,kk
         call cpy_mJpacJatL(ufxcum_loc(I_0H,J_0H,k) )
! ufxcum(iatls,jatl,k)=-ufxcum(ipacn,jpac,k)
      end do

      CALL HALO_UPDATE(ogrid, vfxcum_loc, FROM=NORTH)
c
      do 9 j=J_0, J_1
      jb = PERIODIC_INDEX(j+1, jj)
! jb=mod(j,jj)+1
      do 9 l=1,isp(j)
c
      do 10 i=ifp(j,l),ilp(j,l)
 10   coldiv(i)=0.
c
      do 11 k=1,kk
      do 11 i=ifp(j,l),ilp(j,l)
      ib=mod(i,ii)+1
      p_loc(i,j,k+1)=p_loc(i,j,k)+dp_loc(i,j,k+nn)
      hordiv(i,j,k)=(ufxcum_loc(ib,j,k)-ufxcum_loc(i,j,k)
     .            +vfxcum_loc(i,jb,k)-vfxcum_loc(i,j,k))*scp2i_loc(i,j)
      coldiv(i)=coldiv(i)+hordiv(i,j,k)
c
cdiag if (i.eq.itest .and. j.eq.jtest) write (*,103) i,j,k,
cdiag. 'mass flux',ufxcum(i,j,k)*scp2i(i,j),vfxcum(i,j,k)*scp2i(i,j),
cdiag.  coldiv(i),vfxcum(i,jb,k)*scp2i(i,j),ufxcum(ib,j,k)*scp2i(i,j)
 11   continue
 103  format (2i5,i3,2x,a/f21.2/f14.2,2f7.2/f21.2)
c
c --- adjust initial layer thickness to cancel effect of column divergence
c
      do 12 i=ifp(j,l),ilp(j,l)
      util1_loc(i,j)=(p_loc(i,j,kk+1)+coldiv(i))/p_loc(i,j,kk+1)
 12   util2_loc(i,j)=1./util1_loc(i,j)
c
      do 13 k=1,kk
      do 13 i=ifp(j,l),ilp(j,l)
      tracer_loc(i,j,k,:)=tracer_loc(i,j,k,:)*util2_loc(i,j)
 13   dpinit_loc(i,j,k)=dpinit_loc(i,j,k)*util1_loc(i,j)
c
c --- compute the various terms in the continuity equation integrated
c --- over time interval since last call to -tradv0-
c --- the continuity eqn is split into horiz. and vert. terms as follows:
c ---        (dpfinl-dpinit) + hordiv + verdiv = 0
c
      do 15 i=ifp(j,l),ilp(j,l)
 15   vertfx(i,j,1)=0.
c
      do 16 k=1,kk
      ka=max(1,k-1)
      do 16 i=ifp(j,l),ilp(j,l)
      verdiv=(dpinit_loc(i,j,k)-dp_loc(i,j,k+nn))-hordiv(i,j,k)
 16   vertfx(i,j,k)=vertfx(i,j,ka)+verdiv	!  flx thru botm of lyr k
 9    continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional: check balance of terms in continuity eqn.
cdiag i=itest
cdiag j=jtest
cdiag ib=mod(i,ii)+1
cdiag jb=mod(j,jj)+1
cdiag write (*,'(2i5,a/a)') i,j,
cdiag. '  trcadv -- time-integrated continuity eqn diagnostics:',
cdiag.  '     thknss_tndcy  horiz.flxdiv   vert.flxdiv      residuum'
cdiag do k=1,kk
cdiag   thkchg=dp(i,j,k+nn)-dpinit(i,j,k)
cdiag   fluxdv=(ufxcum(ib,j,k)-ufxcum(i,j,k)
cdiag.         +vfxcum(i,jb,k)-vfxcum(i,j,k))*scp2i(i,j)
cdiag   if (k.eq.1) then
cdiag     write (*,104) k,thkchg,fluxdv,vertfx(i,j,k),
cdiag.    thkchg+fluxdv+vertfx(i,j,k)
cdiag   else
cdiag     write (*,104) k,thkchg,fluxdv,vertfx(i,j,k)-vertfx(i,j,k-1),
cdiag.    thkchg+fluxdv+vertfx(i,j,k)-vertfx(i,j,k-1)
cdiag   end if
cdiag end do
 104  format (i3,4f14.1)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
cdiag call totals(dpinit,tracer(1,1,1,1),
cdiag.            dpinit,tracer(1,1,1,2),'bfore fct3d')
c
      do nt=1,ntrcr
cdiag   do k=1,kk
cdiag     write (string,'(a,i2,a,i3)') 'trcr',nt,' bfor adv',k
cdiag     call findmx(ip,tracer(1,1,k,nt),idm,ii1,jj,string)
cdiag   end do
c
c       write (string,'(a,i2)') 'trcr_bef',nt
c       call prt9x9(tracer(1,1,1,nt),itest,jtest,0.,100.,string)
c       call prt9x9(u,itest,jtest,0.,100.,'usfc')
c
        call fct3d(2,tracer_loc(I_0H,J_0H,1,nt),ufxcum_loc,vfxcum_loc,
     .    vertfx,scp2_loc,scp2i_loc,dpinit_loc,dp_loc(I_0H,J_0H,1+nn))
c
cdiag   do k=1,kk
cdiag     write (string,'(a,i2,a,i3)') 'trcr',nt,' aftr adv',k
cdiag     call findmx(ip,tracer(1,1,k,nt),idm,ii1,jj,string)
cdiag   end do
c
        do k=1,kk
          call cpy_p_par(tracer_loc(I_0H,J_0H,k,nt))
        end do
c
c       write (string,'(a,i2)') 'trcr_aft',nt
c       call prt9x9(tracer(1,1,1,nt),itest,jtest,0.,100.,string)
      end do
c
cdiag call totals(dp(1,1,1+nn),tracer(1,1,1,1),
cdiag.            dp(1,1,1+nn),tracer(1,1,1,2),'after fct3d')
c
      if (AM_I_ROOT()) write (*,'(a)') 'tracer transport done'
c
      return
      end
c
c
c> Revision history:
c>
c> Dec. 2004 - made code cyclic in both -i- and -j- direction
c> Feb. 2005 - added multiple tracer capability
c> Mar. 2006 - added bering strait exchange logic
c
      subroutine fct3d(iord,fld,u,v,w,scal,scali,fco1,fc1)
c
c --- fully 3-d version of advfct.f
c
c  fld    - transported mixing ratio, e.g., salinity or temperature
c  u,v,w  - mass fluxes (x time step) satisfying continuity equation
c           (w(k) = mass flux per unit area across  b o t t o m  of layer k)
c  scal   - grid cell size
c  scali  - inverse of scal
c  fco,fc - depth of the layer at previous and new time step
c
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : itest, jtest, onemu
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT, HALO_UPDATE, NORTH, SOUTH,
     &                         haveLatitude, GLOBALSUM, broadcast
      implicit none
      integer i,j,k,l,ia,ib,ja,jb
c
      real fld(idm,J_0H:J_1H,kdm),u(idm,J_0H:J_1H,kdm),
     .     v(idm,J_0H:J_1H,kdm),w(idm,J_0H:J_1H,kdm),
     .     scal(idm,J_0H:J_1H),scali(idm,J_0H:J_1H),
     .     fco1(idm,J_0H:J_1H,kdm),fc1(idm,J_0H:J_1H,kdm)

      real fco(idm,J_0H:J_1H,kdm),fc(idm,J_0H:J_1H,kdm),
     .     vertfx(idm,J_0H:J_1H,kdm),vertdv(idm,J_0H:J_1H,kdm)
      real fmx(idm,J_0H:J_1H),fmn(idm,J_0H:J_1H),flp(idm,J_0H:J_1H),
     .     fln(idm,J_0H:J_1H),flx(idm,J_0H:J_1H),fly(idm,J_0H:J_1H),
     .     uan(idm,J_0H:J_1H),van(idm,J_0H:J_1H),
     .     flxdiv(idm,J_0H:J_1H),clipj(J_0H:J_1H),vlumj(J_0H:J_1H)
      real bforej(J_0H:J_1H),afterj(J_0H:J_1H)

! Filling data for 2 HALO positions
! We use 2 HALO operations to obtain data corresponding to 2nd HALO position.
! After 1st HALO operation, the data is copied into a temporary array
! such that data in j location corresponds to j-1 location in the fld array.
! Next, after filling SOUTH HALO position in the temporary array, we copy data
! corrsponding to J_0H position of the temporary array into J_0H-1 position
! of the array with 2 SOUTH HALOs.
      ! tfld = temporary arrays to store fld shifted in ja direction
      ! fld2 array contains fld(i,jaa) data
      ! i.e., fld(i,jaa,k)) = fld2(i,jaa))
      ! k subscript is dropped to save memory

      real fld2(idm,J_0H-1:J_1), tfld(idm,J_0H:J_1H)

      real a(kdm),b(kdm),c(kdm),athird,dx,fcdx,yl,yr
      real q,clip,vlume,amount,bfore,after,slab,dslab,thkchg,
     .     fluxdv,epsil
      integer iord,ip1,im1,jp1,jm1,kp,jaa
      character string*16
      logical wrap,recovr
      data recovr/.true./
      parameter (athird=1./3.)
c
c --- if iord=1, scheme reduces to simple donor cell scheme.
      parameter (epsil=1.e-11)

!========================================================
c
      do 21 j=J_0, J_1
      do 21 k=1,kk
      do 21 l=1,isp(j)
      do 21 i=ifp(j,l),ilp(j,l)
      fco(i,j,k)=max(0.,fco1(i,j,k))
      fc (i,j,k)=max(0.,fc1 (i,j,k))
      if (fco(i,j,k).lt.1.e-30) fco(i,j,k)=0.
      if (fc (i,j,k).lt.1.e-30) fc (i,j,k)=0.
 21   continue
c
c --- optional: check mass conservation
cdiag i=itest
cdiag j=jtest
cdiag jb=mod(j,jj)+1
cdiag write (*,'(2i5,a/a)') i,j,
cdiag. '  fct3d -- time-integrated continuity eqn diagnostics:',
cdiag.  '     thknss_tndcy  horiz.flxdiv   vert.flxdiv      residuum'
cdiag do k=1,kk
cdiag thkchg=fc(i,j,k)-fco(i,j,k)
cdiag fluxdv=(u(i+1,j,k)-u(i,j,k)
cdiag.       +v(i,jb ,k)-v(i,j,k))*scali(i,j)
cdiag if (k.eq.1) then
cdiag   write (*,103) k,thkchg,fluxdv,w(i,j,k),
cdiag.  thkchg+fluxdv+w(i,j,k)
cdiag else
cdiag   write (*,103) k,thkchg,fluxdv,w(i,j,k)-w(i,j,k-1),
cdiag.  thkchg+fluxdv+w(i,j,k)-w(i,j,k-1)
cdiag end if
cdiag end do
 103  format (i3,4f14.1)
cc c
cc       do 6 j=1,jdm
cc       do 6 i=1,idm
cc  6    flxdiv(i,j)=0.
cc       do 5 k=1,kk
cc       do 9 j=1,jj
cc       jb=mod(j,jj)+1
cc       do 9 l=1,isp(j)
cc       if (k.eq.1) then
cc         do 10 i=ifp(j,l),ilp(j,l)
cc  10     flxdiv(i,j)=(u(i+1,j,k)-u(i,j,k)+v(i,jb,k)-v(i,j,k))*scali(i,j)
cc      .    +w(i,j,k)           +fc(i,j,k)-fco(i,j,k)
cc       else
cc         do 12 i=ifp(j,l),ilp(j,l)
cc  12     flxdiv(i,j)=(u(i+1,j,k)-u(i,j,k)+v(i,jb,k)-v(i,j,k))*scali(i,j)
cc      .    +w(i,j,k)-w(i,j,k-1)+fc(i,j,k)-fco(i,j,k)
cc       end if
cc  9    continue
cc       call findmx(ip,flxdiv,idm,ii1,jj,'mass consv')
cc       write (*,*) 'shown below: mass consv. residual in layer',k
cc       call zebra(flxdiv,idm,idm-1,jdm)
cc  5    continue
c
c --- get vertical flux by summing -fld- over upstream slab of thickness -w-
c
      do 26 j=J_0, J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 26 l=1,isp(j)
c
c --- fill massless cells with data from layer above or below
      do 17 k=kk-1,1,-1
      do 17 i=ifp(j,l),ilp(j,l)
 17   fld(i,j,k)=(fld(i,j,k)*fco(i,j,k)+fld(i,j,k+1)*onemu)
     .          /(           fco(i,j,k)+             onemu)
      do 18 k=2,kk
      do 18 i=ifp(j,l),ilp(j,l)
      fld(i,j,k)=(fld(i,j,k)*fco(i,j,k)+fld(i,j,k-1)*onemu)
     .          /(           fco(i,j,k)+             onemu)
      if (fld(i,j,k).le.1.e-30) fld(i,j,k)=0.
      if (k.eq.2.and.fld(i,j,1).le.1.e-30) fld(i,j,1)=0.
 18   continue
c
      do 26 i=ifp(j,l),ilp(j,l)
c
c --- fit 0th, 1st, or 2nd deg. polynomial to tracer in each cell
      a(1 )=fld(i,j,1 )
      b(1 )=0.
      c(1 )=0.
      a(kk)=fld(i,j,kk)
      b(kk)=0.
      c(kk)=0.
      do 16 k=2,kk-1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise constant method:
ccc      a(k)=fld(i,j,k)
ccc      b(k)=0.
ccc      c(k)=0.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise linear method:
c --- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
ccc      a(k)=fld(i,j,k)
ccc      b(k)=0.
ccc      if (fld(i,j,k).le.min(fld(i,j,k-1),fld(i,j,k+1)) .or.
ccc     .    fld(i,j,k).ge.max(fld(i,j,k-1),fld(i,j,k+1))) then
ccc        b(k)=0.
ccc      else if ((fld(i,j,k+1)-fld(i,j,k-1))*(fld(i,j,k-1)+fld(i,j,k+1)
ccc     .  -2.*fld(i,j,k)).gt.0.) then
ccc        b(k)=fld(i,j,k)-fld(i,j,k-1)
ccc      else
ccc        b(k)=fld(i,j,k+1)-fld(i,j,k)
ccc      end if
ccc      c(k)=0.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise parabolic method:
c --- fit parabola a+bx+cx^2 to tracer in each cell (-.5 < x < +.5)
      yl=.5*(fld(i,j,k-1)+fld(i,j,k))
      yr=.5*(fld(i,j,k+1)+fld(i,j,k))
      a(k)=1.5*fld(i,j,k)-.25*(yl+yr)
      b(k)=yr-yl
      c(k)=6.*(.5*(yl+yr)-fld(i,j,k))
      if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-fld(i,j,k))) then
c --- apex of parabola occurs inside interval [-.5,+.5], implying an over-
c --- or undershoot situation. change curve to prevent over/undershoots.
        if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-fld(i,j,k))) then
c --- put apex of parabola on edge of interval [-.5,+.5]
          if ((yr-yl)*(.5*(yl+yr)-fld(i,j,k)) .gt. 0.) then
c --- apex at x=-.5
            a(k)=.25*(3.*fld(i,j,k)+yl)
            c(k)=3.*(fld(i,j,k)-yl)
            b(k)=c(k)
          else
c --- apex at x=+.5
            a(k)=.25*(3.*fld(i,j,k)+yr)
            c(k)=3.*(fld(i,j,k)-yr)
            b(k)=-c(k)
          end if
        else			!  -1/6 < x < +1/6
c --- moving apex won't help. replace parabola by constant.
          a(k)=fld(i,j,k)
          b(k)=0.
          c(k)=0.
        end if
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 16   continue
c
      do 23 k=1,kk-1
      slab=onemu
      if (w(i,j,k).lt.0.) then			! interface moves down
        amount=slab*fld(i,j,k+1)
        kp=k
 24     kp=kp+1
        if (slab.ge.-w(i,j,k)) goto 23
        if (fco(i,j,kp).gt.0.) then
          dslab=min(slab+fco(i,j,kp),-w(i,j,k))
     .         -min(slab            ,-w(i,j,k))
          dx=dslab/fco(i,j,kp)
          fcdx=a(kp)
     .        +b(kp)*.5*(dx-1.)			!  not needed in pcm
     .        +c(kp)*(.25-dx*(.5-dx*athird))	!  not needed in pcm,plm
          amount=amount+fcdx*dslab
          slab=slab+dslab
        end if
        if (kp.lt.kk) go to 24
      else if (w(i,j,k).gt.0.) then		! interface moves up
        amount=slab*fld(i,j,k)
        kp=k+1
 25     kp=kp-1
        if (slab.ge.w(i,j,k)) goto 23
        if (fco(i,j,kp).gt.0.) then
          dslab=min(slab+fco(i,j,kp), w(i,j,k))
     .         -min(slab            , w(i,j,k))
          dx=dslab/fco(i,j,kp)
          fcdx=a(kp)
     .        +b(kp)*.5*(1.-dx)			!  not needed in pcm
     .        +c(kp)*(.25-dx*(.5-dx*athird))	!  not needed in pcm,plm
          amount=amount+fcdx*dslab
          slab=slab+dslab
        end if
        if (kp.gt.2) go to 25
      end if
 23   vertfx(i,j,k)=w(i,j,k)*amount/slab
c
      vertfx(i,j,kk)=0.			!  don't allow flux through bottom
      vertdv(i,j,1)=vertfx(i,j,1)
      do 26 k=2,kk
      vertdv(i,j,k)=vertfx(i,j,k)-vertfx(i,j,k-1)
c
 26   continue
c
      bfore=0.
      after=0.

      CALL HALO_UPDATE(ogrid,ip, FROM=NORTH+SOUTH)
      CALL HALO_UPDATE(ogrid,iv, FROM=SOUTH)
c
      do 4 k=1,kk
c
      do 14 j=J_0, J_1
      bforej(j)=0.
      do 14 l=1,isp(j)
      do 14 i=ifp(j,l),ilp(j,l)
 14   bforej(j)=bforej(j)+fld(i,j,k)*fco(i,j,k)*scal(i,j)
c
c --- compute antidiffusive (high- minus low-order) fluxes
c
      call cpy_p_par(fld(I_0H,J_0H,k))

      CALL HALO_UPDATE(ogrid,fld(:,:,k), FROM=NORTH)
      CALL HALO_UPDATE(ogrid,fld(:,:,k), FROM=SOUTH)
       fld2(:,J_0H:J_1) = fld(:,J_0H:J_1,k)
       tfld(:,J_0:J_1)  = fld(:,J_0H:J_1-1,k)
      CALL HALO_UPDATE(ogrid,tfld, FROM=SOUTH)
       fld2(:,J_0H-1) = tfld(:,J_0H)
c
      do 11 j=J_0, J_1
      ja  = PERIODIC_INDEX(j-1, jj)
      jb  = PERIODIC_INDEX(j+1, jj)
      jaa = PERIODIC_INDEX(j-2, jj)
c
      do 2 l=1,isu(j)
      do 2 i=ifu(j,l),ilu(j,l)
      if (u(i,j,k).ge.0.) then
        q=fld(i-1,j,k)
      else
        q=fld(i  ,j,k)
      end if
      flx(i,j)=u(i,j,k)*q
      q=fld(i,j,k)+fld(i-1,j,k)				!  2nd order
      if (ip(i+1,j)+iu(i-1,j).eq.2)
     .  q=1.125*q-.125*(fld(i+1,j,k)+fld(i-2,j,k))	!  4th order
 2    uan(i,j)=.5*q*u(i,j,k)-flx(i,j)
c
      do 3 l=1,isv(j)
      do 3 i=ifv(j,l),ilv(j,l)
      if (v(i,j,k).ge.0.) then
        q=fld(i,ja ,k)
      else
        q=fld(i,j  ,k)
      end if
      fly(i,j)=v(i,j,k)*q
      q=fld(i,ja ,k)+fld(i,j,k)				!  2nd order
      if (ip(i,jb )+iv(i,ja).eq.2)
     .  q=1.125*q-.125*(fld(i,jb ,k)+fld2(i,jaa))	!  4th order
!    .  q=1.125*q-.125*(fld(i,jb ,k)+fld(i,jaa,k))	!  4th order
 3    van(i,j)=.5*q*v(i,j,k)-fly(i,j)
c
      do 11 l=1,isp(j)
      do 11 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (ip(ib,j).eq.0) ib=i
      !ja=mod(j-2+jj,jj)+1
      ja = PERIODIC_INDEX(j-1, jj)

      if (ip(i,ja).eq.0) ja=j
      !jb=mod(j     ,jj)+1
      jb = PERIODIC_INDEX(j+1, jj)

      if (ip(i,jb).eq.0) jb=j
      fmx(i,j)=max(fld(i,j,k),
     .   fld(ia,j,k),fld(ib,j,k),fld(i,ja,k),fld(i,jb,k))
      fmn(i,j)=min(fld(i,j,k),
     .   fld(ia,j,k),fld(ib,j,k),fld(i,ja,k),fld(i,jb,k))
      if (k.lt.kk) then
        if (w(i,j,k  ).lt.0.) then
          fmx(i,j)=max(fmx(i,j),vertfx(i,j,k  )/w(i,j,k  ))
          fmn(i,j)=min(fmn(i,j),vertfx(i,j,k  )/w(i,j,k  ))
        end if
      end if
      if (k.gt.1) then
        if (w(i,j,k-1).gt.0.) then
          fmx(i,j)=max(fmx(i,j),vertfx(i,j,k-1)/w(i,j,k-1))
          fmn(i,j)=min(fmn(i,j),vertfx(i,j,k-1)/w(i,j,k-1))
        end if
      end if
 11   continue
c
cdiag write (string,'(a6,i2)') 'fmx k=',k
cdiag call findmx(ip,fmx,idm,ii1,jj,string(1:8))
cdiag write (string,'(a6,i2)') 'fmn k=',k
cdiag call findmx(ip,fmx,idm,ii1,jj,string(1:8))

      do 22 j=J_0, J_1
      do 22 l=1,isp(j)
      flx(ifp(j,l)  ,j)=0.
      flx(ilp(j,l)+1,j)=0.
      uan(ifp(j,l)  ,j)=0.
      uan(ilp(j,l)+1,j)=0.
  22  continue
c
      do 33 i=1,ii1
      wrap=jfv(i,1).eq.1	! true if j=1 and j=jj are both water points
      do 33 l=1,jsp(i)
      j=jfp(i,l)
      if (haveLatitude(ogrid, J=j)) then
        if (j.gt.1 .or. .not.wrap) then
        fly(i,j)=0.
        van(i,j)=0.
        end if
      end if
      j=mod(jlp(i,l),jj)+1
      if (haveLatitude(ogrid, J=j)) then
        if (j.gt.1 .or. .not.wrap) then
        fly(i,j)=0.
        van(i,j)=0.
        end if
      end if
   33 continue
c
cdiag i=itest
cdiag j=jtest
cdiag write (*,101) 'advem(1)',i,j,k,fld(i-1,j,k),u(i,j,k),
cdiag. fld(i,j-1,k),v(i,j,k),fld(i,j,k),v(i,j+1,k),fld(i,j+1,k),
cdiag.  u(i+1,j,k),fld(i+1,j,k)
 101  format(a,2i5,i3,f18.3/es39.2/f19.3,es11.2,f9.3,
     .es11.2,f9.3/es39.2/f39.3)

      CALL HALO_UPDATE(ogrid,fly, FROM=NORTH)
c
      do 61 j=J_0, J_1
      jb = PERIODIC_INDEX(j+1, jj)
      if (recovr) vlumj(j)=0.
      if (recovr) clipj(j)=0.
      do 61 l=1,isp(j)
      do 61 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*scali(i,j)
c
cdiag if (i.eq.itest .and. j.eq.jtest)
cdiag. write (*,'(2i5,i3,a,4f10.5,es9.2)') i,j,k,'  fc,fco,divs:',
cdiag.  fc(i,j,k),fco(i,j,k),flxdiv(i,j),vertdv(i,j,k),
cdiag.  fc(i,j,k)-fco(i,j,k)+flxdiv(i,j)+vertdv(i,j,k)
c
      q=fld(i,j,k)*fco(i,j,k)-flxdiv(i,j)-vertdv(i,j,k)
      amount=max(0.,fmn(i,j)*fc(i,j,k),min(q,fmx(i,j)*fc(i,j,k)))
      if (recovr) then
        vlumj(j)=vlumj(j)+scal(i,j)*fc(i,j,k)
        clipj(j)=clipj(j)+(q-amount)*scal(i,j)
      end if
      fld(i,j,k)=(fld(i,j,k)*onemu+amount)/(onemu+fc(i,j,k))
      fld(i,j,k)=max(fmn(i,j),min(fld(i,j,k),fmx(i,j)))	!  just to be sure...
 61   continue
c
cdiag call findmx(ip,fld(1,1,k),idm,ii1,jj,'fld after 61')
c
      if (iord.le.1) go to 100
c
c --- at each grid point, determine the ratio of the largest permissible
c --- pos. (neg.) change in -fld- to the sum of all incoming (outgoing) fluxes

      CALL HALO_UPDATE(ogrid,van, FROM=NORTH)
c
      do 12 j=J_0, J_1
      jb  = PERIODIC_INDEX(j+1, jj)
      do 12 l=1,isp(j)
      do 12 i=ifp(j,l),ilp(j,l)
      flp(i,j)=(fmx(i,j)-fld(i,j,k))*fc(i,j,k)
     ./((max(0.,uan(i,j))-min(0.,uan(i+1,j))
     .  +max(0.,van(i,j))-min(0.,van(i,jb ))+epsil)*scali(i,j))
c
 12   fln(i,j)=(fmn(i,j)-fld(i,j,k))*fc(i,j,k)
     ./((min(0.,uan(i,j))-max(0.,uan(i+1,j))
     .  +min(0.,van(i,j))-max(0.,van(i,jb ))-epsil)*scali(i,j))
c
c---- limit antidiffusive fluxes
c
      call cpy_p_par(flp(I_0H,J_0H))
      call cpy_p_par(fln(I_0H,J_0H))

      CALL HALO_UPDATE(ogrid,flp, FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,fln, FROM=SOUTH)
c
      do 8 j=J_0, J_1
      ja  = PERIODIC_INDEX(j-1, jj)
c
      do 7 l=1,isu(j)
      do 7 i=ifu(j,l),ilu(j,l)
      if (uan(i,j).ge.0.) then
        clip=min(1.,flp(i,j),fln(i-1,j))
      else
        clip=min(1.,fln(i,j),flp(i-1,j))
      end if
 7    flx(i,j)=uan(i,j)*clip
c
      do 8 l=1,isv(j)
      do 8 i=ifv(j,l),ilv(j,l)
      if (van(i,j).ge.0.) then
        clip=min(1.,flp(i,j),fln(i,ja ))
      else
        clip=min(1.,fln(i,j),flp(i,ja ))
      end if
 8    fly(i,j)=van(i,j)*clip
c
cdiag i=itest
cdiag j=jtest
cdiag write (*,101) 'advem(2)',i,j,k,fld(i-1,j,k),u(i,j,k),
cdiag. fld(i,j-1,k),v(i,j,k),fld(i,j,k),v(i,j+1,k),fld(i,j+1,k),
cdiag.  u(i+1,j,k),fld(i+1,j,k)

      CALL HALO_UPDATE(ogrid,fly, FROM=NORTH)
c
      do 62 j=J_0, J_1
      jb  = PERIODIC_INDEX(j+1, jj)
      do 62 l=1,isp(j)
      do 62 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*scali(i,j)
      q=fld(i,j,k)*fc(i,j,k)-flxdiv(i,j)
      amount=max(0.,fmn(i,j)*fc(i,j,k),min(q,fmx(i,j)*fc(i,j,k)))
      if (recovr) clipj(j)=clipj(j)+(q-amount)*scal(i,j)
      fld(i,j,k)=(fld(i,j,k)*onemu+amount)/(onemu+fc(i,j,k))
      fld(i,j,k)=max(fmn(i,j),min(fld(i,j,k),fmx(i,j)))	! just to be sure...
 62   continue
c
cdiag call findmx(ip,fld(1,1,k),idm,ii1,jj,'fld after 62')
c
 100  continue
c
c --- recover 'clipped' amount and return to field layer by layer
c
      if (recovr) then
        vlume=0.
        clip=0.
c
        call GLOBALSUM(ogrid,vlumj,vlume, all=.true.)
        call GLOBALSUM(ogrid,clipj,clip , all=.true.)
!       do 19 j=J_0, J_1
!       vlume=vlume+vlumj(j)
!19     clip=clip+clipj(j)
c
        if (vlume.ne.0.) then
          clip=clip/vlume
cdiag     write (*,'(a,i2,a,es11.3)') 'k=',k,'  tracer drift in fct3d',
cdiag.     -clip
          do 13 j=J_0, J_1
          do 13 l=1,isp(j)
          do 13 i=ifp(j,l),ilp(j,l)
 13       fld(i,j,k)=fld(i,j,k)+clip
        end if
      end if
c
      do 15 j=J_0, J_1
      afterj(j)=0.
      do 15 l=1,isp(j)
      do 15 i=ifp(j,l),ilp(j,l)
 15   afterj(j)=afterj(j)+fld(i,j,k)*fc(i,j,k)*scal(i,j)
c
      call compBforeAfter(bfore,after,bforej,afterj)
c
 4    continue			! k loop

      call broadcast(ogrid, bfore)
      call broadcast(ogrid, after)
c
      if (bfore.ne.0.)
     . write (*,'(a,3es14.6,es11.2)') 'fct3d conservation:',
     .  bfore,after,after-bfore,(after-bfore)/bfore
      q=1.
      if (after.ne.0.) q=bfore/after
      write (*,'(a,f11.6)') 'fct3d: multiply tracer field by',q
ccc   if (q.gt.1.1 .or. q.lt..9) stop '(excessive nonconservation)'
      if (q.gt.2.0 .or. q.lt..5) stop '(excessive nonconservation)'
c
      do 20 j=J_0, J_1
      do 20 k=1,kk
      do 20 l=1,isp(j)
      do 20 i=ifp(j,l),ilp(j,l)
 20   fld(i,j,k)=fld(i,j,k)*q
c
      return
      end
!========================================================
      subroutine compBforeAfter(bfore,after,bforej,afterj)

      USE HYCOM_DIM, only : idm, jdm, kdm, J_0H,  J_1H, jj, ogrid
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
      implicit none

      real bfore,after
      real bforej(J_0H:J_1H),afterj(J_0H:J_1H)
      real bforej_G(jdm),afterj_G(jdm)
      integer j

      call bforejAfterjGath(bforej,afterj,bforej_G,afterj_G)

      !doing this on root for reproducibility
      if (AM_I_ROOT()) then
        do j=1,jj
          bfore=bfore+bforej_G(j)
          after=after+afterj_G(j)
        end do
      endif  ! AM_I_ROOT

      end subroutine compBforeAfter
!========================================================
      subroutine bforejAfterjGath(bforej,afterj,bforej_G,afterj_G)

      USE HYCOM_DIM, only : idm, jdm, kdm, J_0H,  J_1H, ogrid
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA
      implicit none
      real :: bforej(J_0H:J_1H),afterj(J_0H:J_1H)
      real :: bforej_G(jdm),afterj_G(jdm)

      call pack_data( ogrid,  bforej,  bforej_G  )
      call pack_data( ogrid,  afterj,  afterj_G  )

      end subroutine bforejAfterjGath
!========================================================
