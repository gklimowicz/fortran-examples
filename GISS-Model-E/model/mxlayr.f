#include "hycom_mpi_hacks.h"
      subroutine mxlayr(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 0.9.6
c
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : dotrcr,theta,onem,onecm,epsil,salmin
     &     ,sigjmp,nstep,delt1,acurcy,time,onemm,huge,thref,g
     &     ,spcifh,thkdff,thkmin,baclin,itest,jtest
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, SOUTH, NORTH, GLOBALSUM
     &     , AM_I_ROOT
      implicit none
c
      integer,intent(IN) :: m,n,mm,nn,k1m,k1n
      integer i,j,k,l,ja,jb,kn,knp1
      real turgen,dpth,ekminv,obuinv,ex,alf1,alf2,cp1,cp3,ape,cc4,spe,
     .     em,en,ea1,ea2,em1,em2,em3,em4,em5,ustar3,thkold,thknew,
     .     pnew,q,temdp,saldp,thermg,small,tem,sal,rho,sup,slo,
     .     ttem(kdm),ssal(kdm),dens(kdm),pres(kdm+1),delp(kdm),
     .     sum1,sum2,siglo,zup,zlo,s1,s2,s3,smax,smin,p1,p2,buoyfl,
     .     totem,tosal,tndcyt,tndcys,athird,trac(kdm,ntrcr),
     .     trcdp(ntrcr),flxhi,flxlo,ctrast,pa,pb
      logical vrbos
      data ctrast/.001/		!  vertical density contrast within mixed layer
      parameter (athird=1./3.)
      real sigocn,dsigdt,dsigds
      external sigocn,dsigdt,dsigds
c
      data ea1, ea2, em1, em2, em3, em4, em5
     .   /0.60,0.30,0.45,2.60,1.90,2.30,0.60/		! Gaspar coefficients
c
c --- advect mixed layer depth with vertically averaged velocity field
c
      do 9 j=J_0,J_1
      do 9 k=1,kk
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
 9    p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
c
      CALL HALO_UPDATE(ogrid,dpmixl,FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,p     ,FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,v     ,FROM=NORTH+SOUTH)
c
      do 2 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
c
      do 3 l=1,isu(j)
      do 3 i=ifu(j,l),ilu(j,l)
      uflux(i,j)=0.
      do 4 k=1,kk
      kn=k+nn
      pb=min(dpmixl(i,j,m)+dpmixl(i-1,j,m),p(i,j,k+1)+p(i-1,j,k+1))
      pa=min(dpmixl(i,j,m)+dpmixl(i-1,j,m),p(i,j,k  )+p(i-1,j,k  ))
      if (pa.eq.pb) exit
 4    uflux(i,j)=uflux(i,j)
     .  +.25*(u(i-1,j,kn)+u(i+1,j,kn)+2.*u(i,j,kn))*(pb-pa)
 3    uflux(i,j)=uflux(i,j)*(dpmixl(i,j,m)-dpmixl(i-1,j,m))*scuy(i,j)/
     .  min(dpmixl(i,j,m)+dpmixl(i-1,j,m),pbot(i,j)+pbot(i-1,j))
c
      do 5 l=1,isv(j)
      do 5 i=ifv(j,l),ilv(j,l)
      vflux(i,j)=0.
      do 6 k=1,kk
      kn=k+nn
      pb=min(dpmixl(i,j,m)+dpmixl(i,ja ,m),p(i,j,k+1)+p(i,ja ,k+1))
      pa=min(dpmixl(i,j,m)+dpmixl(i,ja ,m),p(i,j,k  )+p(i,ja ,k  ))
      if (pa.eq.pb) exit
 6    vflux(i,j)=vflux(i,j)
     .  +.25*(v(i,ja ,kn)+v(i,jb ,kn)+2.*v(i,j,kn))*(pb-pa)
 5    vflux(i,j)=vflux(i,j)*(dpmixl(i,j,m)-dpmixl(i,ja ,m))*scvx(i,j)/
     .  min(dpmixl(i,j,m)+dpmixl(i,ja ,m),pbot(i,j)+pbot(i,ja ))
 2    continue
c
      CALL HALO_UPDATE(ogrid,vflux,FROM=NORTH)
c
      do 13 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 13 l=1,isp(j)
      do 13 i=ifp(j,l),ilp(j,l)
      thkold=dpmixl(i,j,n)
      dpmixl(i,j,n)=min(pbot(i,j),dpmixl(i,j,n)
     .   -.5*(uflux(i+1,j)+uflux(i,j)
     .       +vflux(i,jb )+vflux(i,j))*scp2i(i,j)*delt1)
 13   dpmixl(i,j,m)=.5*dpmixl(i,j,m)+.25*(thkold+dpmixl(i,j,n))
c
      do 1 j=J_0,J_1
      vrbos=.false.
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
cdiag vrbos=i.eq.itest .and. j.eq.jtest
c
c --- extract single column from 3-d fields
      pres(1)=p(i,j,1)
      do 7 k=1,kk
      kn=k+nn
      if (dotrcr) trac(k,:)=tracer(i,j,k,:)
      ttem(k)=temp(i,j,kn)
      ssal(k)=saln(i,j,kn)
      dens(k)=th3d(i,j,kn)
      delp(k)=dp(i,j,kn)
 7    pres(k+1)=pres(k)+delp(k)
c
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f8.2,f8.1))
      if (vrbos) write (*,103) nstep,itest,jtest,
     .'  entering mxlayr:  temp    saln    dens    thkns    dpth',(k,
     .ttem(k),ssal(k),dens(k),delp(k)/onem,pres(k+1)/onem,k=1,kk)
c
c --- store 'old' t/s column integral in totem/tosal (diagnostic use only)
      totem=0.
      tosal=0.
      do k=1,kk
        totem=totem+ttem(k)*delp(k)
        tosal=tosal+ssal(k)*delp(k)
      end do
c
c --- diagnose initial mixed layer depth
c
cc       temdp=ttem(1)*delp(1)
cc       saldp=ssal(1)*delp(1)
cc c
cc       do 11 k=2,kk
cc       if (delp(k).le.0.) go to 11
cc       tem=(temdp+ttem(k)*delp(k))/pres(k+1)
cc       sal=(saldp+ssal(k)*delp(k))/pres(k+1)
cc       rho=sigocn(tem,sal)
cc       if (rho.le.dens(1)+ctrast) then
cc         temdp=temdp+ttem(k)*delp(k)
cc         saldp=saldp+ssal(k)*delp(k)
cc       else
cc c
cc c --- old mixed layer depth (thkold) falls into layer -k-. diagnose -thkold-
cc c --- by requiring that avg. mxlayr density matches surface density + ctrast
cc c
cc         q=delp(k)-pres(k+1)*(rho    -dens(1)-ctrast)/
cc      .                      (dens(k)-dens(1)-ctrast)
cc         q=max(0.,q)
cc         thkold=pres(k)+q
cc         tem=(temdp+ttem(k)*q)/thkold
cc         sal=(saldp+ssal(k)*q)/thkold
cc         dens(1)=sigocn(tem,sal)
cc         go to 12
cc       end if
cc  11   continue
cc       thkold=pres(kk+1)
cc  12   continue
c
      thkold=max(thkmin*onem,dpmixl(i,j,m))
c
c --- ----------------------------------------
c --- slab mixed layer entrainment/detrainment
c --- ----------------------------------------
c
      ustar3=ustar(i,j)**3
c
c --- buoyfl = buoyancy flux, w_prime_buoyancy_prime_bar (m^2/sec^3)
c --- note: surface density increases (column is destabilized) if buoyfl < 0
c
      tem=.5*(temp(i,j,k1m)+temp(i,j,k1n))
      sal=.5*(saln(i,j,k1m)+saln(i,j,k1n))
      buoyfl=-g*thref**2*(dsigds(tem,sal)*salflx(i,j)
     .                   +dsigdt(tem,sal)*surflx(i,j)/spcifh)
c
c --- determine turb.kin.energy generation due to wind stirring (ustar3) and
c --- diabatic forcing (buoyfl). ustar,buoyfl are computed in subr. -thermf-.
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- option 1 :   k r a u s  -  t u r n e r    mixed-layer t.k.e.  closure
c
      em=0.8*exp(-pres(2)/(50.*onem))   !   hadley centre choice (orig.: 1.25)
      en=0.15                           !   hadley centre choice (orig.: 0.4)
      thermg=-.5*((en+1.)*buoyfl+(en-1.)*abs(buoyfl))
      turgen=delt1*(2.*em*g*ustar3/thref+thkold*thermg)/thref**2
c
c --- find monin-obukhov length in case of receding mixed layer (turgen < 0).
c --- the monin-obukhov length is found by stipulating turgen = 0.
c
      if (turgen.lt.0.) then
        thknew=-2.*em*g*ustar3/min(-epsil,thref*thermg)
      else
        thknew=thkold
      end if
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- option 2 :    g a s p a r   mixed-layer t.k.e. closure
c
      dpth=thkold/onem
      ekminv=abs(corio(i,j))/max(epsil,ustar(i,j))
      obuinv=buoyfl/max(epsil,ustar3)
      ex=exp(min(50.,dpth*obuinv))
      alf1=ea1+ea2*max(1.,2.5*dpth*ekminv)*ex
      alf2=ea1+ea2*ex
      cp1=((1.-em5)*(alf1/alf2)+.5*em4)*athird
      cp3=max(0.,(em4*(em2+em3)-(alf1/alf2)*(em2+em3-em3*em5))*athird)
      ape=cp3*ustar3-cp1*dpth*buoyfl
c
      if(ape.lt.0.) then					! detrainment
        turgen=(g*delt1/thref**3)*ape
        thknew=min(thkold,g*cp3/(thref*cp1*max(epsil,obuinv)))
c
      else							! entrainment
        cc4=2.*em4/(em1*em1) * alf1*alf1
        spe=(em2+em3)*ustar3-0.5*dpth*buoyfl
        turgen=(g*delt1/thref**3)*(sqrt((.5*ape-cp1*spe)**2
     .          +2.*cc4*ape*spe)-(.5*ape+cp1*spe))/(cc4-cp1)
        thknew=thkold
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- sum1,sum2 are used to evaluate pot.energy changes during entrainment
      sum1=dens(1)*thkold
      sum2=dens(1)*thkold**2
c
c --- find thknew in case of mixed layer deepening (turgen>0).
c --- entrain as many layers as needed to deplete -turgen-.
c
      if (turgen.ge.0.) then
        do 85 k=2,kk
        if (pres(k+1).gt.thkold) then
          q=max(pres(k),thkold)
          pnew=min(pres(k+1),(2.*turgen+dens(k)*q**2-sum2)/
     .                        max(epsil,dens(k)*q   -sum1))
c --- stop iterating for 'thknew' as soon as pnew < k-th interface pressure
          if (pnew.lt.pres(k)) go to 86
          thknew=pnew
          sum1=sum1+dens(k)*(pres(k+1)   -q   )
          sum2=sum2+dens(k)*(pres(k+1)**2-q**2)
        end if
 85     continue
 86     continue
      end if
c
c --- move mxlyr depth closer to deepest interface embedded in mixed layer.
c --- this is to artificially reduce mixing across mixed layer base.
c
ccc      do 84 k=2,kk
ccc      if (thknew.gt.pres(k) .and. thknew.lt.pres(k+1)) then
ccc        q=(thknew-pres(k))/delp(k)
ccc        q=q*q
ccc        q=q*q
ccc        q=0.			! reduce mxlyr depth to nearest interface
ccc        thknew=pres(k)*(1.-q)+pres(k+1)*q
ccc      end if
ccc 84   continue
c
c --- don't allow mixed layer to get too deep or too shallow.
      thknew=min(pres(kk+1),max(thkmin*onem,delp(1),thknew))
c
c --- integrate t/s over new mixed layer depth
c
      if (dotrcr) trcdp(:)=trac(1,:)*delp(1)
      temdp=ttem(1)*delp(1)
      saldp=ssal(1)*delp(1)
c
      do 15 k=2,kk
      if (pres(k).lt.thknew) then
        q=min(thknew,pres(k+1))-min(thknew,pres(k))
        if (dotrcr) trcdp(:)=trcdp(:)+trac(k,:)*q
        temdp=temdp+ttem(k)*q
        saldp=saldp+ssal(k)*q
      end if
 15   continue
c
      if (vrbos) write (*,'(i9,2i5,a,2f9.3)')
     .  nstep,i,j,'  old/new mixed layer depth:',thkold/onem,thknew/onem
c
c --- distribute thermohaline forcing over new mixed layer depth
c
      ttem(1)=(temdp+surflx(i,j)*delt1*g/spcifh)/thknew
      ssal(1)=(saldp+salflx(i,j)*delt1*g       )/thknew
      dens(1)=sigocn(ttem(1),ssal(1))
      if (dotrcr) trac(1,:)=trcdp(:)/thknew
c
c --- homogenize water mass properties down to new mixed layer depth
c
      do 14 k=2,kk
      if (pres(k+1).le.thknew) then
        if (dotrcr) trac(k,:)=trac(1,:)
        ttem(k)=ttem(1)
        ssal(k)=ssal(1)
        dens(k)=dens(1)
      else if (pres(k).lt.thknew) then
c
        if (vrbos)
     .  write (*,'(i9,2i5,i3,a,3f9.3,25x,2f9.3)') nstep,i,j,k,
     .   '  p_k,thknew,p_k+1,t_1,t_k=',pres(k)/onem,thknew/onem,
     .    pres(k+1)/onem,ttem(1),ttem(k)
c
        q=(thknew-pres(k))/delp(k)
        if (dotrcr) trac(k,:)=trac(1,:)*q+trac(k,:)*(1.-q)
        ttem(k)=ttem(1)*q+ttem(k)*(1.-q)
        ssal(k)=ssal(1)*q+ssal(k)*(1.-q)
        dens(k)=sigocn(ttem(k),ssal(k))
      end if
 14   continue
c
      if (vrbos) write (*,103) nstep,itest,jtest,
     .'  exiting mxlayr:   temp    saln    dens    thkns    dpth',(k,
     .ttem(k),ssal(k),dens(k),delp(k)/onem,pres(k+1)/onem,k=1,kk)
c
c --- compare 'old' with 'new' t/s column integral (diagnostic use only)
c
        tndcyt=-totem
        tndcys=-tosal
        totem=10.*pres(kk+1)
        tosal=35.*pres(kk+1)
        do k=1,kk
          tndcyt=tndcyt+ttem(k)*delp(k)
          tndcys=tndcys+ssal(k)*delp(k)
        end do
        tndcyt=tndcyt-surflx(i,j)*delt1*g/spcifh
        tndcys=tndcys-salflx(i,j)*delt1*g
 101  format (2i5,a,1p,2e16.8,e9.1)
        if (abs(tndcyt).gt.acurcy*totem) write (*,101) i,j,
     .  '  mxlayr - bad temp.intgl.',totem,tndcyt,tndcyt/totem
        if (abs(tndcys).gt.acurcy*tosal) write (*,101) i,j,
     .  '  mxlayr - bad saln.intgl.',tosal,tndcys,tndcys/tosal
ccc        if (max(abs(tndcyt/totem),abs(tndcys/tosal)).gt.
ccc     .  1.e-9) write (*,'(i9,2i5,3x,a,1p,3e10.2/22x,a,3e10.2)')
ccc     .  nstep,i,j,'total saln,srf.flux,tndcy:',tosal/g,
ccc     .  salflx*delt1,tndcys/g,'total temp,srf.flux,tndcy:',
ccc     .  totem/g,surflx*delt1,tndcyt*spcifh/g
c
c --- put single column back into 3-d fields
      do 8 k=1,kk
      kn=k+nn
      if (dotrcr) tracer(i,j,k,:)=trac(k,:)
      temp(i,j,kn)=ttem(k)
      saln(i,j,kn)=ssal(k)
 8    th3d(i,j,kn)=dens(k)

 1    dpmixl(i,j,n)=thknew
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional: smooth mixed layer depth
c
      q=thkdff*baclin
      CALL HALO_UPDATE(ogrid,dpmixl,FROM=SOUTH)
c
      do 51 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
c --- limit fluxes to avoid intersecting the sea floor
      do 50 l=1,isu(j)
      do 50 i=ifu(j,l),ilu(j,l)
      flxhi= .25*(pbot(i  ,j)-dpmixl(i  ,j,n))*scp2(i  ,j)
      flxlo=-.25*(pbot(i-1,j)-dpmixl(i-1,j,n))*scp2(i-1,j)
 50   uflux(i,j)=min(flxhi,max(flxlo,
     .   q*(dpmixl(i-1,j,n)-dpmixl(i,j,n))*scuy(i,j)))
c
      do 51 l=1,isv(j)
      do 51 i=ifv(j,l),ilv(j,l)
      flxhi= .25*(pbot(i,j  )-dpmixl(i,j  ,n))*scp2(i,j  )
      flxlo=-.25*(pbot(i,ja )-dpmixl(i,ja ,n))*scp2(i,ja )
 51   vflux(i,j)=min(flxhi,max(flxlo,
     .   q*(dpmixl(i,ja ,n)-dpmixl(i,j,n))*scvx(i,j)))
c
      CALL HALO_UPDATE(ogrid,vflux,FROM=NORTH)
c
      do 38 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 38 l=1,isp(j)
      do 38 i=ifp(j,l),ilp(j,l)
      dpmixl(i,j,n)=dpmixl(i,j,n)-(uflux(i+1,j)-uflux(i,j)
     .                            +vflux(i,jb )-vflux(i,j))*scp2i(i,j)
 38   dpmxav(i,j)=dpmxav(i,j)+dpmixl(i,j,n)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c --- ---------------
c --- momentum mixing
c --- ---------------
c
c --- homogenize -u- down to new mixed layer depth
c
      small=1.e-4
      CALL HALO_UPDATE(ogrid,dpmixl,FROM=SOUTH)
c
      do 31 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 32 l=1,isu(j)
c
      do 33 i=ifu(j,l),ilu(j,l)
      klist(i,j)=-1
      util1(i,j)=max(dpu(i,j,k1n),.5*(dpmixl(i,j,n)+dpmixl(i-1,j,n)))
      uflux(i,j)=u(i,j,k1n)*dpu(i,j,k1n)
      util2(i,j)=dpu(i,j,k1n)
 33   pu(i,j,2)=dpu(i,j,k1n)
c
      do 34 k=2,kk
      kn=k+nn
      knp1=min(k+1,kk)+nn
      do 34 i=ifu(j,l),ilu(j,l)
      pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,kn)
      if (pu(i,j,k+1).le.util1(i,j)) then
        uflux(i,j)=uflux(i,j)+u(i,j,kn)*dpu(i,j,kn)
        util2(i,j)=util2(i,j)+          dpu(i,j,kn)
      else if (pu(i,j,k).lt.util1(i,j)) then
c --- divide layer k into 2 sublayers. upper one belongs to mixed layer
        zup=util1(i,j)-pu(i,j,k)
        zlo=dpu(i,j,kn)-zup
        s1=u(i,j,kn-1)
        s2=u(i,j,kn  )
        if (k.eq.kk .or. (k.lt.kk .and. dpu(i,j,knp1).lt.onemm)) then
          s3=2.*s2-s1
        else
          s3=u(i,j,kn+1)
        end if
c --- define 'bounding box'
        smax=max(s1,s2,s3)
        smin=min(s1,s2,s3)
        if (s2.lt.smin+small .or. s2.gt.smax-small) then
          sup=s2
          slo=s2
        else
          slo=s3
          sup=(s2*dpu(i,j,kn)-slo*zlo)/zup
          if (sup.gt.smin-small .and. sup.lt.smax+small) go to 36
          sup=s1
          slo=(s2*dpu(i,j,kn)-sup*zup)/zlo
          if (slo.gt.smin-small .and. slo.lt.smax+small) go to 36
          write (*,100) nstep,i,j,'  possible',' error in unmixing u',
     .      dpu(i,j,kn)/onem,zup/onem,zlo/onem,s1,s2,s3,
     .      (s2*dpu(i,j,kn)-slo*zlo)/zup,(s2*dpu(i,j,kn)-sup*zup)/zlo
          sup=s2
          slo=s2
        end if
 36     uflux(i,j)=uflux(i,j)+sup*zup
        util2(i,j)=util2(i,j)+    zup
        util3(i,j)=u(i,j,kn)
        u(i,j,kn)=slo
        klist(i,j)=k
      end if
 34   continue
 100  format (i9,2i5,2a,3f9.3/3f10.4,2(2x,2f10.4))
c
      do 35 i=ifu(j,l),ilu(j,l)
 35   u(i,j,k1n)=uflux(i,j)/util2(i,j)
c
      do 32 k=2,kk
      kn=k+nn
      do 32 i=ifu(j,l),ilu(j,l)
      q=max(0.,min(1.,(util1(i,j)-pu(i,j,k))/(dpu(i,j,kn)+epsil)))
      if (q.eq.0. .and. k.eq.klist(i,j)) then
        u(i,j,kn)=util3(i,j)
      else
        u(i,j,kn)=u(i,j,k1n)*q+u(i,j,kn)*(1.-q)
      end if
 32   continue
c
c --- homogenize -v- down to new mixed layer depth
c
      do 52 l=1,isv(j)
c
      do 53 i=ifv(j,l),ilv(j,l)
      klist(i,j)=-1
      util1(i,j)=max(dpv(i,j,k1n),.5*(dpmixl(i,j,n)+dpmixl(i,ja ,n)))
      vflux(i,j)=v(i,j,k1n)*dpv(i,j,k1n)
      util2(i,j)=dpv(i,j,k1n)
 53   pv(i,j,2)=dpv(i,j,k1n)
c
      do 54 k=2,kk
      kn=k+nn
      knp1=min(k+1,kk)+nn
      do 54 i=ifv(j,l),ilv(j,l)
      pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,kn)
      if (pv(i,j,k+1).le.util1(i,j)) then
        vflux(i,j)=vflux(i,j)+v(i,j,kn)*dpv(i,j,kn)
        util2(i,j)=util2(i,j)+          dpv(i,j,kn)
      else if (pv(i,j,k).lt.util1(i,j)) then
c --- divide layer k into 2 sublayers. upper one belongs to mixed layer
        zup=util1(i,j)-pv(i,j,k)
        zlo=dpv(i,j,kn)-zup
        s1=v(i,j,kn-1)
        s2=v(i,j,kn  )
        if (k.eq.kk .or. (k.lt.kk .and. dpv(i,j,knp1).lt.onemm)) then
          s3=2.*s2-s1
        else
          s3=v(i,j,kn+1)
        end if
c --- define 'bounding box'
        smax=max(s1,s2,s3)
        smin=min(s1,s2,s3)
        if (s2.lt.smin+small .or. s2.gt.smax-small) then
          sup=s2
          slo=s2
        else
          slo=s3
          sup=(s2*dpv(i,j,kn)-slo*zlo)/zup
          if (sup.gt.smin-small .and. sup.lt.smax+small) go to 56
          sup=s1
          slo=(s2*dpv(i,j,kn)-sup*zup)/zlo
          if (slo.gt.smin-small .and. slo.lt.smax+small) go to 56
          write (*,100) nstep,i,j,'  possible',' error in unmixing v',
     .      dpv(i,j,kn)/onem,zup/onem,zlo/onem,s1,s2,s3,
     .      (s2*dpv(i,j,kn)-slo*zlo)/zup,(s2*dpv(i,j,kn)-sup*zup)/zlo
          sup=s2
          slo=s2
        end if
 56     vflux(i,j)=vflux(i,j)+sup*zup
        util2(i,j)=util2(i,j)+    zup
        util3(i,j)=v(i,j,kn)
        v(i,j,kn)=slo
        klist(i,j)=k
      end if
 54   continue
c
      do 55 i=ifv(j,l),ilv(j,l)
 55   v(i,j,k1n)=vflux(i,j)/util2(i,j)
c
      do 52 k=2,kk
      kn=k+nn
      do 52 i=ifv(j,l),ilv(j,l)
      q=max(0.,min(1.,(util1(i,j)-pv(i,j,k))/(dpv(i,j,kn)+epsil)))
      if (q.eq.0. .and. k.eq.klist(i,j)) then
        v(i,j,kn)=util3(i,j)
      else
        v(i,j,kn)=v(i,j,k1n)*q+v(i,j,kn)*(1.-q)
      end if
 52   continue
 31   continue
c
      return
      end
c
c
c> Revision history
c>
c> Mar. 2000 - conversion to SI units
c> Apr. 2000 - conversion of unmixing calculations to double precision
c> Apr. 2000 - changed 'unmixing' calculations to real*8
c> Apr. 2000 - changed dimensions of turgen in light of its use in loop 85
c> May  2000 - added logic to limit mixed layer detrainment rate (-detrmx-)
c> June 2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> June 2000 - total rewrite
c> July 2000 - discarded sublayers in thermodynamic calculations
c> July 2000 - switched sign convention for vertical fluxes (now >0 if down)
c> Nov. 2000 - corrected definition of thkold [now pres(mxbot+1)]
c> Apr. 2001 - eliminated stmt_funcs.h
c> Jun. 2001 - added calculation of -buoyfl- (formerly done in thermf)
c> Jun. 2001 - bypass conv.adjustment (loop 11) in cases of zero layer thknss
c> Dec. 2001 - added logic to reduce mixing across mxlyr bottom (loop 84)
c> Feb. 2002 - fixed bug introduced while rewriting loop 85
c> Apr. 2002 - introduced even-odd time average of t/s in -buoyfl- calculation
c> May  2002 - introduced -ctrast- logic to diagnose old ML depth
c> June 2002 - fixed bug in sum1/sum2 calculation in loop 85
c> Sep. 2003 - added statement basing trac(1) on -trcdp-
c> Sep. 2003 - added logical switch to enable/disable tracer redistribution
c> Feb. 2005 - added multiple tracer capability
c> Aug. 2008 - mixed layer depth (2 time slots) remembered and advected
