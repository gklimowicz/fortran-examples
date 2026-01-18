      subroutine advc1d(kk,fld,z,dz,vrbos,errcon)
c
c --- ppm-based 1-dim transport routine, extracted from fct3d.f
c --- originator: rainer bleck^
c
c --- kk       - number of layers
c --- fld(kk)  - mixing ratio of dependent variable in each layer
c --- z(kk+1)  - interface depth; z(k+1)-z(k) = thickness of layer k
c --- dz(kk+1) - elements of -fld- travel from z(k)-dz(k) to z(k) during
c ---            present time step
c --- vrbos    - if .true., print diagnostic messages
c --- errcon   - set to .true. if error condition encountered
c
      implicit none
      integer, intent(IN) :: kk
      real,    intent(IN) :: z(kk+1),dz(kk+1)
      real, intent(INOUT) :: fld(kk)
      logical, intent(IN) :: vrbos
      logical,intent(OUT) :: errcon
c
      real thk(kk),flux(kk+1),div(kk)
      real a(kk),b(kk),c(kk),dx,fcdx,yl,yr
      real amount,bfore,after,scale,slab,dslab
      integer k,kp
      real, parameter :: athird=1./3.
      real, parameter :: small=1.e-11   !  (for z,dz given in meters)
c
cdiag if (vrbos) write (*,100)
cdiag  do k=1,kk
cdiag  write(*,'(2i5,2e12.4)')k,kk,z(k),dz(k)
cdiag  enddo
cdiag  write (*,100)
cdiag. 'entering advc1d:  old loc''n  new loc''n   variable',
cdiag. (k,z(k),z(k)+dz(k),fld(k),k=1,kk),kk+1,z(kk+1),z(kk+1)+dz(kk+1)
 100  format (a/(i15,2f11.3,es12.3))
c
      bfore=0.
      scale=1.e-33
      do k=1,kk
        thk(k)=z(k+1)-z(k)
        bfore=bfore+fld(k)*thk(k)
        scale=max(scale,abs(fld(k)))
      end do
c
c --- exit if parcels overtake each other
      errcon=.false.
      do k=2,kk
        if (z(k)+dz(k).lt.z(k-1)+dz(k-1) .or.
     .      z(k)+dz(k).gt.z(k+1)+dz(k+1)) then
          write (*,'(a,i3)') 'error advc1d -- crossover at k =',k
          errcon=.true.
          return
          return
        end if
      end do
c
c --- start by filling massless cells with data from layer(s) above or below
c
      do 17 k=kk-1,1,-1
 17   fld(k)=(fld(k)*thk(k)+fld(k+1)*small)
     .      /(       thk(k)+         small)
      do 18 k=2,kk
 18   fld(k)=(fld(k)*thk(k)+fld(k-1)*small)
     .      /(       thk(k)+         small)
c
c --- fit 0th, 1st, or 2nd deg. polynomial to -fld- in each cell
      a(1 )=fld(1 )
      b(1 )=0.
      c(1 )=0.
      a(kk)=fld(kk)
      b(kk)=0.
      c(kk)=0.
c
      do 16 k=2,kk-1
c --- uncomment one of the following 3 options to activate pcm,plm,ppm resp.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise constant method:
ccc      a(k)=fld(k)
ccc      b(k)=0.
ccc      c(k)=0.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise linear method:
c --- fit linear function a+bx to tracer in each cell (-.5 < x < +.5)
ccc      a(k)=fld(k)
ccc      b(k)=0.
ccc      if (fld(k).le.min(fld(k-1),fld(k+1)) .or.
ccc     .    fld(k).ge.max(fld(k-1),fld(k+1))) then
ccc        b(k)=0.
ccc      else if ((fld(k+1)-fld(k-1))*(fld(k-1)+fld(k+1)
ccc     .  -2.*fld(k)).gt.0.) then
ccc        b(k)=fld(k)-fld(k-1)
ccc      else
ccc        b(k)=fld(k+1)-fld(k)
ccc      end if
ccc      c(k)=0.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise parabolic method:
c --- fit parabola a+bx+cx^2 to tracer in each cell (-.5 < x < +.5)
      yl=.5*(fld(k-1)+fld(k))
      yr=.5*(fld(k+1)+fld(k))
      a(k)=1.5*fld(k)-.25*(yl+yr)
      b(k)=yr-yl
      c(k)=6.*(.5*(yl+yr)-fld(k))
      if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-fld(k))) then
c --- apex of parabola lies inside interval [-.5,+.5], implying an over-
c --- or undershoot situation. change curve to prevent over/undershoots.
        if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-fld(k))) then
c --- put apex of parabola on edge of interval [-.5,+.5]
          if ((yr-yl)*(.5*(yl+yr)-fld(k)) .gt. 0.) then
c --- apex at x=-.5
            a(k)=.25*(3.*fld(k)+yl)
            c(k)=3.*(fld(k)-yl)
            b(k)=c(k)
          else
c --- apex at x=+.5
            a(k)=.25*(3.*fld(k)+yr)
            c(k)=3.*(fld(k)-yr)
            b(k)=-c(k)
          end if
        else                    !  -1/6 < x < +1/6
c --- moving apex won't help. replace parabola by constant.
          a(k)=fld(k)
          b(k)=0.
          c(k)=0.
        end if
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 16   continue
c
c --- get flux by summing -fld- over upstream slab of thickness -dz-
c
      do 23 k=1,kk+1
      slab=small
      if (dz(k).lt.0.) then                     ! -fld- moves up
        if (k.eq.kk+1) then
          amount=0.                             ! influx from outside domain
        else
          amount=slab*fld(k)
          kp=k-1
 24       kp=kp+1
          if (kp.le.kk) then
            if (thk(kp).gt.0.) then
              dslab=min(slab+thk(kp),-dz(k),z(kk+1)-z(k))
     .             -min(slab        ,-dz(k),z(kk+1)-z(k))
              dx=dslab/thk(kp)
              fcdx=a(kp)
     .          +b(kp)*.5*(dx-1.)               !  not needed in pcm
     .          +c(kp)*(.25-dx*(.5-dx*athird))  !  not needed in pcm,plm
              amount=amount+fcdx*dslab
              slab=slab+dslab
            end if
          else                                  ! reaching outside domain
            slab=-dz(k)
          end if
          if (slab.lt.-dz(k)) go to 24
        end if
      else if (dz(k).gt.0.) then                ! -fld- moves down
        if (k.eq.1) then
          amount=0.                             ! influx from outside domain
        else
          amount=slab*fld(k-1)
          kp=k
 25       kp=kp-1
          if (kp.ge.1) then
            if (thk(kp).gt.0.) then
              dslab=min(slab+thk(kp), dz(k),z(k)-z(1))
     .             -min(slab        , dz(k),z(k)-z(1))
              dx=dslab/thk(kp)
              fcdx=a(kp)
     .          +b(kp)*.5*(1.-dx)               !  not needed in pcm
     .          +c(kp)*(.25-dx*(.5-dx*athird))  !  not needed in pcm,plm
              amount=amount+fcdx*dslab
              slab=slab+dslab
            end if
          else                                  ! reaching outside domain
            slab=dz(k)
          end if
          if (slab.lt.dz(k)) go to 25
        end if
      else                                      !  dz = 0
        amount=0.
      end if
 23   flux(k)=dz(k)*amount/slab
c
      do 26 k=1,kk
 26   div(k)=flux(k+1)-flux(k)
c
      do 4 k=1,kk
      amount=fld(k)*thk(k)-div(k)
 4    fld(k)=(fld(k)*small+amount)/(small+thk(k))
c
      after=flux(kk+1)-flux(1)                  !  account for outflow loss
      do k=1,kk
        after=after+fld(k)*thk(k)
      end do
      if (abs(bfore-after)*kk.gt.1.e-9*scale*z(kk+1)) then
        write (*,104) '  advc1d - bad column intgl.:',bfore,after
        errcon=.true.
      end if
 104  format (a,2es15.7)
c
cdiag if (vrbos) write (*,100)
cdiag. 'exiting  advc1d:  old loc''n  new loc''n   variable',
cdiag. (k,z(k),z(k)+dz(k),fld(k),k=1,kk),kk+1,z(kk+1),z(kk+1)+dz(kk+1)
c
      return
      end
c
c
c> Revision history:
c>
c> May  2001 - added u/v regridding
c> June 2001 - added interlayer ("diapycnal") mass flux diagnostics
c> Aug. 2001 - added -kn- to list of private variables in loop 13
c> Aug. 2001 - various refinements, incl. multi-layer ingest
c> Dec. 2001 - softened enforcement of minimum layer thickness constraint
c> Feb. 2002 - restricted ingest of multiple layers from below
c> Mar. 2002 - added passive tracer
c> Apr. 2002 - fixed -diaflx- diagnostics by defining -dpold-
c> June 2002 - changed 'slak' logic - now based on gradual restoring
c> June 2002 - curtail downward growth of layers (esp. lowest hybrid layer)
c> May  2003 - introduced variable -dpsum- to make -p_hat- k-dependent
c> Aug. 2003 - added option to maintain 10 cm min.thickness throughout column
c> Sep. 2003 - added logical switch to enable/disable tracer redistribution
c> Sep. 2003 - provided choice between T/S and rho/S conservation
c> Nov. 2003 - replaced dp0 in cushn call by dp0 appropriate for layer k-1
c> Nov. 2003 - accelerated relaxation time for 1st interface (slak x 10)
c> Nov. 2004 - allowed -dplist- values to shrink near equator
c> Nov. 2004 - conversion to piecewise linear advection scheme (PLM)
c> Nov. 2004 - extended PLM advection to velocity field
c> Dec. 2004 - changed from PLM to PPM (piecewise parabolic)
c> Dec. 2004 - added code to update dpu/dpv(n) (loop 9 and dpudpv call)
c> Apr. 2005 - changed sigjmp to 10*sigjmp in formulae computing p_hat
c> May  2005 - ppmadv: bounded wdth by 1 (loop 4)
c> Mar. 2006 - changed pres(k) to pres(k+1) in p_hat calc'n after stmt label 6
c> June 2006 - if layer 1 touches sea floor, don't split it to restore target
c> July 2006 - introduced depth-dependent tapering function for 'slak'
