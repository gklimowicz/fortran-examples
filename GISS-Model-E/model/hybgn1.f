#include "hycom_mpi_hacks.h"
#include "rundeck_opts.h"
      subroutine hybgen(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 0.9.12
c --- this version allows switching between T/S and rho/S conservation
c --- and between pcm and ppm
c
      USE HYCOM_DIM
      USE HYCOM_SCALARS, only : dotrcr,theta,onem,onecm,epsil,salmin
     &  ,sigjmp,nstep,delt1,acurcy,time,onemm,huge,itest,jtest,dplist
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE, SOUTH, GLOBALSUM,
     &     AM_I_ROOT
      use TimeConstants_mod, only: SECONDS_PER_DAY
      implicit none
c
c --- ---------------------
c --- hybrid grid generator (coordinate restoration exclusively by "dilution")
c --- ---------------------
c
      integer i,j,k,l,m,n,mm,nn,kn,k1m,k1n,ja
c
      real delp,dp0,dp0abv,dpsum,zinteg,tinteg,sinteg,uvintg,
     .     uvscl,phi,plo,pa,pb,dsgdt,dsgds,q1,q2,
     .     rho_lo,tem_lo,sal_lo,rho_up,tem_up,sal_up,tem,sal,p_hat,q,
     .     torho,totem,tosal,totrc,totuv,tndrho,tndtem,tndsal,tndtrc,
     .     tdcyuv,scale,displ(kdm+1),sumrho,sumtem,sumsal
      real targt(kdm+1),dens(kdm),ttem(kdm),ssal(kdm),pres(kdm+1),
     .     uold(kdm),vold(kdm),pold(kdm+1),pnew(kdm+1),
     .     trac(kdm,ntrcr)
      logical abort,tscnsv,vrbos,useppm
      data tscnsv/.true./       ! if true, go with T/S conservation
      data abort/.false./
      data useppm/.true./
      real sigocn,tofsig,dsigdt,dsigds,cushn
      external sigocn,tofsig,dsigdt,dsigds,cushn
      integer lyr,k1,kp,iunit,lpunit,ko,nt
     .        ,ntot2,ntot3,nwrk2,nwrk3
     .        ,ntot2d(J_0H:J_1H),ntot3d(J_0H:J_1H)
     .        ,nwrk2d(J_0H:J_1H),nwrk3d(J_0H:J_1H)
      real :: anwrk, anwrkd(J_0H:J_1H)
      character info*16
      data uvscl/0.02/          !  2 cm/s
      real,parameter :: tfreez=-1.8
      real,parameter :: scalt=-30.,scals=10.   ! oceanic t/s range
c
css   real,parameter :: slak=.5/86400.  ! intfc nudging time scale: 2 days
css   real,parameter :: slak=1./86400.  ! intfc nudging time scale: 1 day
      real,parameter :: slak=2./86400.  ! intfc nudging time scale: 12 hrs
c --- linear taper functions (latitude and depth-dependent) for slak
      real tapr,wgtf,slakf
      tapr(q)=1.+9.*max(0.,1.-.02e-4*q) ! q = pressure (Pa)
      wgtf(q)=(abs(q)-50.)*.1           ! 0->1 for q=50->60
      slakf(q)=min(0.7,max(    tapr(p_hat)*slak*delt1,   ! q = latitude (deg)
     .             0.7*wgtf(q)+tapr(p_hat)*slak*delt1*(1.-wgtf(q))))
c
      do 32 j=J_0,J_1
      do 32 l=1,isp(j)
      do 32 i=ifp(j,l),ilp(j,l)
      vrbos=i.eq.itest .and. j.eq.jtest
      if (vrbos) write (*,103) nstep,i,j,
     .  '  entering hybgen:  temp    saln    dens    thkns    dpth',
     .  (k,temp(i,j,k+nn),saln(i,j,k+nn),
     .  th3d(i,j,k+nn),dp(i,j,k+nn)/onem,
     .  p(i,j,k+1)/onem,k=1,kk)
      if (vrbos) write (*,106) nstep,i,j,
     .  '  entering hybgen:  dpthu      u    dpthv      v',
     .  (k,pu(i,j,k+1)/onem,u(i,j,k+nn),
     .     pv(i,j,k+1)/onem,v(i,j,k+nn),k=1,kk)
 32   continue
 103  format (i9,2i5,a/(33x,i3,2f8.3,f8.3,f8.2,f8.1))
 106  format (i9,2i5,a/(33x,i3,2(f8.1,f8.3)))
!-------------------------------------------------------------------
       delp=0;dp0=0;dp0abv=0;dpsum=0;zinteg=0;tinteg=0;sinteg=0;
       uvintg=0;
       phi=0;plo=0;pa=0;pb=0;dsgdt=0;dsgds=0;
       tem_up=0;sal_up=0;rho_up=0.;tem_lo=0;sal_lo=0;rho_lo=0.;
       tem=0;sal=0;p_hat=0;q=0;q1=0;q2=0;
       torho=0;totem=0;tosal=0;totrc=0;totuv=0;tndrho=0;tndtem=0;
       tndsal=0;tndtrc=0;tdcyuv=0;scale=0;displ=0;
       targt=0;dens=0;ttem=0;ssal=0;pres=0;
       uold=0;vold=0;pold=0;pnew=0;trac=0;
!-------------------------------------------------------------------
c
      do 19 j=J_0,J_1
      do 19 k=1,kk
      do 19 l=1,isp(j)
      do 19 i=ifp(j,l),ilp(j,l)
 19   p(i,j,k+1)=p(i,j,k)+dp(i,j,k+nn)
c
      abort=.false.

      do 12 j=J_0, J_1
      ntot2=0
      ntot3=0
      nwrk2=0
      nwrk3=0
c
      do 2 l=1,isp(j)
      do 2 i=ifp(j,l),ilp(j,l)
c
c --- extract t,s,rho column from 3-d grid
c
      pres(1)=p(i,j,1)
      do 3 k=1,kk
      kn=k+nn
      dens(k)=th3d(i,j,kn)
      ttem(k)=temp(i,j,kn)
      ssal(k)=saln(i,j,kn)
      if (dotrcr) trac(k,:)=tracer(i,j,k,:)
      pres(k+1)=pres(k)+dp(i,j,kn)
 3    targt(k)=theta(k)
c
      vrbos=i.eq.itest .and. j.eq.jtest
      if (vrbos) then
        write (*,99) nstep,i,j,'      o l d   p r o f i l e :'
        do k=1,kk,10
        write (*,100) (pres(k1)/onem,k1=k,min(kk+1,k+10))
        write (*,101) (dens(k1),k1=k,min(kk,k+9))
        write (*,102) (ttem(k1),k1=k,min(kk,k+9))
        write (*,102) (ssal(k1),k1=k,min(kk,k+9))
        end do
      end if
 99   format (i9,2i5,a)
 100  format (11f7.1)
 101  format (4x,10f7.2)
 102  format (4x,10f7.2)
c
      torho=0.
      totem=0.
      tosal=0.
      do k=1,kk
        torho=torho+dens(k)*(pres(k+1)-pres(k))
        totem=totem+ttem(k)*(pres(k+1)-pres(k))
        tosal=tosal+ssal(k)*(pres(k+1)-pres(k))
      end do
c
      kp=1
      do 4 k=2,kk
      if (pres(k).lt.pres(kk+1)-onecm) then
        kp=k
      else
c
c --- absorb near-massless layers on sea floor in layer above
c
        q1=max(epsil,pres(k  )-pres(kp))
        q2=max(   0.,pres(k+1)-pres(k ))
        q=q1/(q1+q2)
        if (q.lt.0. .or. q.gt.1.) then
          write (*,*) 'i,j,q1,q2,q=',i,j,q1,q2,q
ccc          abort=.true.
        end if
        ttem(kp)=ttem(kp)*q+ttem(k)*(1.-q)
        ssal(kp)=ssal(kp)*q+ssal(k)*(1.-q)
        dens(kp)=dens(kp)*q+dens(k)*(1.-q)
        ttem(k)=ttem(kp)
        ssal(k)=max(ssal(kp),salmin(k))
        dens(k)=dens(kp)
        if (dotrcr) then
          trac(kp,:)=trac(kp,:)*q+trac(k,:)*(1.-q)
          trac(k,:)=trac(kp,:)
        end if
c
        if (vrbos)
     .   write (*,'(i9,2i5,a,i3,a,i3,5x,a,f8.3)') nstep,i,j,
     .    '  absorb layer',k,' in',kp,'new rho:',dens(kp)
      end if
 4    continue
c
      do 11 k=kp+1,kk
 11   pres(k)=pres(kk+1)
c
      do 23 k=1,kk
 23   dpold(i,j,k)=pres(k+1)-pres(k)
c
c --- is layer touching sea floor to light?
      if (kp.eq.1) go to 10
      k=kp
      if (dens(k).gt.max(dens(k-1),theta(k-1)) .and.
     .    dens(k).le.theta(k)) then
c
c --- water in layer k is too light. split layer into 2 sublayers
c --- matching target densities of layers k-1 and k respectively.
c --- combine upper sublayer with layer k-1.
c
        tem=ttem(k)
        sal=ssal(k)
ccc        scalt=-abs(ttem(k-1)-tem)
ccc        scals= abs(ssal(k-1)-sal)
ccc        if (scalt+scals.eq.0.) go to 10
        dsgdt=dsigdt(tem,sal)*scalt
        dsgds=dsigds(tem,sal)*scals
        q=1./(dsgdt+dsgds)
c
c --- set properties in lower sublayer:
        sal_lo=min(sal+(theta(k)-dens(k))*q*scals,max(ssal(k-1),sal))
        if (tscnsv) then
          tem_lo=max(tem+(theta(k)-dens(k))*q*scalt,tfreez)
          rho_lo=sigocn(tem_lo,sal_lo)
        else
          rho_lo=theta(k)
          tem_lo=max(tofsig(rho_lo,sal_lo),tfreez)
        end if
c
c --- set properties in upper sublayer:
        rho_up=max(dens(k-1),min(theta(k-1),dens(k)))
        q=(rho_lo-dens(k))/max(epsil,rho_lo-rho_up)
        if (q.ge.0. .and. q.le.1.) then
          p_hat=pres(k)*(1.-q)+pres(k+1)*q
          q=(pres(k+1)-p_hat)/max(epsil,p_hat-pres(k))
          sal_up=sal*(1.+q)-sal_lo*q
          if (tscnsv) then
            tem_up=tem*(1.+q)-tem_lo*q
          else
            tem_up=tofsig(rho_up,sal_up)
          end if
c
          if (vrbos) then
            write (*,'(i9,2i5,i3,a,3f7.3,f8.2)') nstep,i,j,k,
     .      '  t,s,th,dp in upper sblyr:',tem_up,sal_up,
     .      sigocn(tem_up,sal_up),(p_hat-pres(k))/onem
            write (*,'(22x,a,3f7.3,f8.2)')
     .      '  t,s,th,dp in lower sblyr:',tem_lo,sal_lo,
     .      rho_lo,(pres(k+1)-p_hat)/onem
            write (*,'(22x,a,1p,2e11.3)') '  scalt,scals =',scalt,scals
          end if
c
c --- combine upper sublayer with layer k-1
          q=(p_hat-pres(k))/max(p_hat-pres(k-1),epsil)
          if (q.lt.0. .or. q.gt.1.) then
            write (*,*) 'q out of range - i,j,k,p_hat,pres(k),q=',
     .                    i,j,k,p_hat,pres(k),q
ccc            abort=.true.
          end if
          ssal(k-1)=  sal_up*q+ssal(k-1)*(1.-q)
          if (tscnsv) then
            ttem(k-1)=tem_up*q+ttem(k-1)*(1.-q)
          else
            ttem(k-1)=tofsig(rho_up,ssal(k-1))
          end if
          if (dotrcr) trac(k-1,:)=trac(k,:)*q+trac(k-1,:)*(1.-q)
c
          if (vrbos)
     .      write (*,'(22x,a,2f7.3)')
     .     '  old/new th(k-1):',dens(k-1),sigocn(ttem(k-1),ssal(k-1)),
     .     '  old/new th(k  ):',dens(k  ),rho_lo
c
          dens(k-1)=sigocn(ttem(k-1),ssal(k-1))
          ttem(k)=tem_lo
          ssal(k)=sal_lo
          dens(k)=rho_lo
          pres(k)=p_hat
        end if
        nwrk2=nwrk2+1
      end if
 10   continue
      ntot2=ntot2+1
c
      do 29 k=1,kk+1
 29   pold(k)=pres(k)
c
c --- try to restore isopycnic conditions by moving layer interfaces
c
      dpsum=0.
      dp0=huge
      do 8 k=1,kk
      ntot3=ntot3+1
c
c --- set lower limits for layer thknss (dp0) and depth of lower intfc (dpsum)
c
      dp0abv=dp0
      dp0=dplist(k)*onem
c
c --- optional: reduce spacing of z layers near equator, but hide transition
c --- in a subtropical latitude band where z layers are least likely to exist
c     if (k.gt.1) dp0=dp0*max(.6,min(1.,(abs(latij(i,j,3))+5.)*.04))
c
c --- reduce layer thickness in shallow spots, creating sigma coord. effect
      if (4*k.lt.kk) dp0=dp0*min(1.,pbot(i,j)/(200.*onem)+.4)
      dpsum=dpsum+dp0
c
c --- maintain constant thickness in layer 1
      if (k.eq.1) then
        p_hat=dp0
        if (p_hat.gt.pres(2)) then
c --- layer 1 is too thin. entrain water from layers below
          p_hat=min(p_hat,pres(2)+
     .          max(onecm,slakf(latij(i,j,3))*(p_hat-pres(2))))
          info='layer too thin  '
          go to 5
        else if (p_hat.lt.pres(2)) then
c --- layer 1 is too thick. expell layer 1 water into layer 2
          p_hat=max(p_hat,pres(2)+
     .          min(-onecm,slakf(latij(i,j,3))*(p_hat-pres(2))))
          info='layer too thick '
c
          if (vrbos) write (*,105)
     .     i,j,k,info,'lower intfc',pres(k+1)/onem,'=>',p_hat/onem
 105      format (2i5,i3,2x,a,'  try moving',2(1x,a,f9.3))
c
          q=(pres(2)-p_hat)/max(pres(3)-p_hat,epsil)
          if (q.lt.0. .or. q.gt.1.) then
            write (*,*) 'i,j,k,pres(2),p_hat,q=',
     .                    i,j,k,pres(2),p_hat,q
ccc            abort=.true.
          end if
          ssal(2)=ssal(2)*(1.-q)+ssal(1)*q
          if (tscnsv) then
            ttem(2)=ttem(2)*(1.-q)+ttem(1)*q
            dens(2)=sigocn(ttem(2),ssal(2))
          else
            dens(2)=dens(2)*(1.-q)+dens(1)*q
            ttem(2)=tofsig(dens(2),ssal(2))
          end if
          pres(2)=p_hat
          nwrk3=nwrk3+1
        end if
        go to 8
      end if				!  k = 1
c
c --- are we dealing with a near-massless layer on the sea floor?
      if (pres(k).eq.pres(kk+1)) dens(k)=max(targt(k),dens(k))
      if (pres(k).gt.pres(kk+1)-onecm) go to 8
c
c --- is lower intfc too close to the surface?
      p_hat=dpsum
      if (k.lt.kk .and. p_hat.gt.pres(k+1)) then
        p_hat=min(p_hat,pres(k+1)+
     .        max(onecm,slakf(latij(i,j,3))*(p_hat-pres(k+1))))
        info='too close to srf'
        go to 5
      end if
c
c --- is density noticeably different from target value?
      if (abs(dens(k)-targt(k)).lt..1*sigjmp) go to 8
c
      if (dens(k).le.targt(k)) go to 7     !  layer too light
c
c --- water in layer k is too  d e n s e . dilute with water from layer k-1
c                              ^^^^^^^^^
      if (k.eq.2) go to 6             !  don't touch layer 1
      q=(targt(k)-dens(k))/max(targt(k)-dens(k-1),sigjmp*10.)
      p_hat=pres(k)*(1.-q)+pres(k+1)*q
c
c --- maintain minimum layer thickess of layer k-1
      p_hat=pres(k-1)+cushn(p_hat-pres(k-1),dp0abv)
      p_hat=min(p_hat,.5*(pres(k-1)+pres(k+1)))
      info='layer too dense '
c
      if (vrbos) write (*,105)
     . i,j,k,info,'upper intfc',pres(k)/onem,'=>',p_hat/onem
c
      if (p_hat.lt.pres(k)) then
c
c --- upper intfc moves up. entrain layer k-1 water into layer k
c
        p_hat=max(p_hat,pres(k-1),pres(k)+
     .        min(-onecm,slakf(latij(i,j,3))*(p_hat-pres(k))))
        if (useppm .and.
     .    abs(dens(k-1)-targt(k-1)).gt..1*sigjmp) then     !  use ppm
          displ(1)=0.
          displ(2)=0.
          displ(3)=p_hat-pres(k)
          displ(4)=0.
          if (vrbos)
     .    write (*,'(2i5,i3,a)') i,j,k,'  entrain from layer above'
          call ppmad3(pres(k-2),displ,ssal(k-2),ssal(k-2),vrbos)
          if (tscnsv) then
            call ppmad3(pres(k-2),displ,ttem(k-2),ttem(k-2),vrbos)
            dens(k-1)=sigocn(ttem(k-1),ssal(k-1))
            dens(k  )=sigocn(ttem(k  ),ssal(k  ))
          else
            call ppmad3(pres(k-2),displ,dens(k-2),dens(k-2),vrbos)
            ttem(k-1)=tofsig(dens(k-1),ssal(k-1))
            ttem(k  )=tofsig(dens(k  ),ssal(k  ))
          end if
        else        !  use pcm
          q=(pres(k)-p_hat)/max(pres(k+1)-p_hat,epsil)
          if (q.lt.0. .or. q.gt.1.) then
            write (*,*) 'i,j,k,p_hat,pres(k),q=',
     .                    i,j,k,p_hat,pres(k),q
ccc            abort=.true.
          end if
          ssal(k)=ssal(k)*(1.-q)+ssal(k-1)*q
          if (tscnsv) then
            ttem(k)=ttem(k)*(1.-q)+ttem(k-1)*q
            dens(k)=sigocn(ttem(k),ssal(k))
          else
            dens(k)=dens(k)*(1.-q)+dens(k-1)*q
            ttem(k)=tofsig(dens(k),ssal(k))
          end if
        end if
        pres(k)=p_hat
        nwrk3=nwrk3+1
c
      else if (p_hat.gt.pres(k)) then		!  p_hat > pres(k)
c
c --- layer k-1 is too thin for allowing upper intfc to move up.  instead,
c --- move upper interface down and entrain layer k water into layer k-1
c
        p_hat=min(p_hat,pres(k+1),pres(k)+
     .        max(onecm,slakf(latij(i,j,3))*(p_hat-pres(k))))
        if (useppm .and. k.lt.kk) then      !  use ppm
          if (vrbos)
     .    write (*,'(2i5,i3,a)') i,j,k,'  detrain into layer above'
          displ(1)=0.
          displ(2)=p_hat-pres(k)
          displ(3)=0.
          displ(4)=0.
          call ppmad3(pres(k-1),displ,ssal(k-1),ssal(k-1),vrbos)
          if (tscnsv) then
            call ppmad3(pres(k-1),displ,ttem(k-1),ttem(k-1),vrbos)
            dens(k-1)=sigocn(ttem(k-1),ssal(k-1))
            dens(k  )=sigocn(ttem(k  ),ssal(k  ))
          else
            call ppmad3(pres(k-1),displ,dens(k-1),dens(k-1),vrbos)
            ttem(k-1)=tofsig(dens(k-1),ssal(k-1))
            ttem(k  )=tofsig(dens(k  ),ssal(k  ))
          end if
        else         !  use pcm
          q=(p_hat-pres(k))/max(p_hat-pres(k-1),epsil)
          if (q.lt.0. .or. q.gt.1.) then
            write (*,*) 'i,j,k,p_hat,pres(k),q=',
     .                    i,j,k,p_hat,pres(k),q
ccc            abort=.true.
          end if
          ssal(k-1)=ssal(k-1)*(1.-q)+ssal(k)*q
          if (tscnsv) then
            ttem(k-1)=ttem(k-1)*(1.-q)+ttem(k)*q
            dens(k-1)=sigocn(ttem(k-1),ssal(k-1))
          else
            dens(k-1)=dens(k-1)*(1.-q)+dens(k)*q
            ttem(k-1)=tofsig(dens(k-1),ssal(k-1))
          end if
        end if
        pres(k)=p_hat
        nwrk3=nwrk3+1
      end if
c
c --- do we need to inflate layer k by lowering  l o w e r  interface?
c
 6    p_hat=pres(k)+dplist(2)*onem
      if (k.lt.kk .and. pres(k+1).lt.pres(kk+1)-onemm .and.
     .    pres(k+1).lt.p_hat) then
        p_hat=min(p_hat,pres(k+1)+
     .        max(onecm,slakf(latij(i,j,3))*(p_hat-pres(k+1))))
        go to 5
      end if
      go to 8
c
c --- water in layer k is too  l i g h t . dilute with water from layer k+1
c                              ^^^^^^^^^
 7    if (k.ge.kk .or. pres(k+1).gt.pres(kk+1)-onemm) go to 8
c
      q=(dens(k)-targt(k))/max(dens(k+1)-targt(k),sigjmp*10.)
      p_hat=pres(k+1)*(1.-q)+pres(k)*q
c
c --- curtail downward growth of layers (esp. lowest hybrid layer)
      p_hat=max(pres(k+1),min(p_hat,.5*(pres(k)+pres(k+2))))
      p_hat=min(p_hat,pres(k+1)+
     .      max(onecm,slakf(latij(i,j,3))*(p_hat-pres(k+1))))
      info='layer too light '
c
 5    p_hat=min(p_hat,pres(k+2))
      if (p_hat.gt.pres(k+1)+onemm) then
c
        if (vrbos) write (*,105)
     .   i,j,k,info,'lower intfc',pres(k+1)/onem,'=>',p_hat/onem
c
        if (useppm .and. k.lt.kk-1 .and.
     .    abs(dens(k+1)-targt(k+1)).gt..1*sigjmp) then   !  use ppm
          if (vrbos)
     .    write (*,'(2i5,i3,a)') i,j,k,'  entrain from layer below'
          displ(1)=0.
          displ(2)=p_hat-pres(k+1)
          displ(3)=0.
          displ(4)=0.
          call ppmad3(pres(k),displ,ssal(k),ssal(k),vrbos)
          if (tscnsv) then
            call ppmad3(pres(k),displ,ttem(k),ttem(k),vrbos)
            dens(k  )=sigocn(ttem(k  ),ssal(k  ))
            dens(k+1)=sigocn(ttem(k+1),ssal(k+1))
          else
            call ppmad3(pres(k),displ,dens(k),dens(k),vrbos)
            ttem(k  )=tofsig(dens(k  ),ssal(k  ))
            ttem(k+1)=tofsig(dens(k+1),ssal(k+1))
          end if
        else          !  use pcm
          q=(p_hat-pres(k+1))/max(p_hat-pres(k),epsil)
          if (q.lt.0. .or. q.gt.1.) then
            write (*,*) 'i,j,k,p_hat,pres(k+1),q=',
     .                    i,j,k,p_hat,pres(k+1),q
ccc            abort=.true.
          end if
          ssal(k)=ssal(k)*(1.-q)+ssal(k+1)*q
          if (tscnsv) then
            ttem(k)=ttem(k)*(1.-q)+ttem(k+1)*q
            dens(k)=sigocn(ttem(k),ssal(k))
          else
            dens(k)=dens(k)*(1.-q)+dens(k+1)*q
            ttem(k)=tofsig(dens(k),ssal(k))
          end if
        end if
        pres(k+1)=p_hat
        nwrk3=nwrk3+1
      end if
 8    continue
c
      if (dotrcr) then
c
c --- evaluate effect of regridding on tracer field(s)
c
        vrbos=i.eq.itest .and. j.eq.jtest
c
        pnew(1)=pres(1)
c
        do 21 k=1,kk
        pnew(k+1)=pres(k+1)
 21     displ(k+1)=pnew(k+1)-pold(k+1)
        displ(   1)=0.
        displ(kk+1)=0.
c
        do nt=1,ntrcr
          scale=1.e-99
          totrc=0.
          do 22 k=1,kk
          totrc=totrc+trac(k,nt)*(pold(k+1)-pold(k))
 22       scale=scale+abs(trac(k,nt))
c
          do k=1,kk-2
            if (pnew(k+1).lt.pnew(k) .or. pnew(k+1).gt.pold(k+2)) then
ccc              write (*,'(a,3i5)') 'trcr monotonicity problems at',i,j,k
ccc              write (*,'(a/(8f9.0))') 'pold:',pold
ccc              write (*,'(a/(8f9.0))') 'pnew:',pnew
              pnew(k+1)=max(pnew(k),min(pnew(k+1),pold(k+2)))
              displ(k+1)=pnew(k+1)-pold(k+1)
            end if
          end do
c
          call ppmadv(kk,pold,displ,trac(1,nt),trac(1,nt),vrbos)
c
          tndtrc=-totrc
          do 20 k=1,kk
          tracer(i,j,k,nt)=trac(k,nt)
 20       tndtrc=tndtrc+trac(k,nt)*(pnew(k+1)-pnew(k))
c
          if (abs(tndtrc)*kk.gt.acurcy*scale*pnew(kk+1))
     .     write (*,104) i,j,'  hybgen - bad trcr.intgl.:',totrc,
     .      tndtrc,tndtrc/(scale*pnew(kk+1))
        end do          !  ntrcr
      end if            !  dotrcr
c
      tndrho=-torho
      tndtem=-totem
      tndsal=-tosal
      do k=1,kk
        tndrho=tndrho+dens(k)*(pres(k+1)-pres(k))
        tndtem=tndtem+ttem(k)*(pres(k+1)-pres(k))
        tndsal=tndsal+ssal(k)*(pres(k+1)-pres(k))
      end do
      if (tscnsv) then
        if (abs(tndtem).gt.acurcy*10.*pres(kk+1))
     .   write (*,104) i,j,'  hybgen - bad temp.intgl.:',totem,
     .    tndtem,tndtem/(10.*pres(kk+1))
      else
        if (abs(tndrho).gt.acurcy*35.*pres(kk+1))
     .   write (*,104) i,j,'  hybgen - bad dens.intgl.:',torho,
     .    tndrho,tndrho/(35.*pres(kk+1))
      end if
      if (abs(tndsal).gt.acurcy*35.*pres(kk+1))
     . write (*,104) i,j,'  hybgen - bad saln.intgl.:',tosal,
     .  tndsal,tndsal/(35.*pres(kk+1))
 104  format (2i5,a,1p,2e15.7,e9.1)
c
      if (vrbos) then
        write (*,99) nstep,i,j,'      n e w   p r o f i l e :'
        do k=1,kk,10
        write (*,100) (pres(k1)/onem,k1=k,min(kk+1,k+10))
        write (*,101) (dens(k1),k1=k,min(kk,k+9))
        write (*,102) (ttem(k1),k1=k,min(kk,k+9))
        write (*,102) (ssal(k1),k1=k,min(kk,k+9))
        end do
      end if
c
c --- put 1-d column back into 3-d grid
c
      do 2 k=1,kk
      kn=k+nn
      th3d(i,j,kn)=dens(k)
      temp(i,j,kn)=ttem(k)
      saln(i,j,kn)=ssal(k)
      p(i,j,k+1)=pres(k+1)
      dp(i,j,kn)=pres(k+1)-pres(k)
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,kn)-dpold(i,j,k))  !  diapyc.flux
 2    continue
c
      ntot2d(j)=ntot2
      ntot3d(j)=ntot3
      nwrk2d(j)=nwrk2
      nwrk3d(j)=nwrk3
 12   continue

      if (abort) stop '(error in hybgen -- q out of bounds)'
c
      do 1 j=J_0,J_1
      do 1 k=1,kk
      do 1 l=1,isp(j)
      do 1 i=ifp(j,l),ilp(j,l)
 1    p(i,j,k+1)=p(i,j,k)+dpold(i,j,k)

      CALL HALO_UPDATE(ogrid,  p,  FROM=SOUTH)
c
      do 88 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      do 88 k=2,kk+1
c
      do 881 l=1,isu(j)
      do 881 i=ifu(j,l),ilu(j,l)
 881  pu(i,j,k)=min(depthu(i,j),.5*(p(i,j,k)+
     .                                      p(i-1,j,k)))
c
      do 882 l=1,isv(j)
      do 882 i=ifv(j,l),ilv(j,l)
 882  pv(i,j,k)=min(depthv(i,j),.5*(p(i,j,k)+
     .                                      p(i,ja ,k)))
 88   continue
c
      do 9 j=J_0,J_1
      do 9 k=1,kk
      do 9 l=1,isp(j)
      do 9 i=ifp(j,l),ilp(j,l)
 9    p(i,j,k+1)=p(i,j,k)+dp(i,j,k+nn)
c
      call pardpudpv(nn)
c
      do 13 j=J_0,J_1
c
c --- integrate -u- over new depth intervals
c
      do 14 l=1,isu(j)
      do 14 i=ifu(j,l),ilu(j,l)
c
      pold(1)=0.
      pnew(1)=0.
      totuv=0.
c
      do 15 k=1,kk
      kn=k+nn
      uold(k)=u(i,j,kn)
      pold(k+1)=pu(i,j,k+1)
      pnew(k+1)=pnew(k)+dpu(i,j,kn)
 15   totuv=totuv+uold(k)*(pold(k+1)-pold(k))
      tdcyuv=-totuv
c
      do 18 k=1,kk
      phi=pnew(k+1)
      plo=pnew(k  )
      if (phi.gt.plo) then
        uvintg=0.
        pb=plo
        do 16 ko=1,kk
        if (pold(ko+1).le.plo) go to 16
        pa=pb
        pb=min(phi,pold(ko+1))
        uvintg=uvintg+uold(ko)*(pb-pa)
        if (pa.ge.phi) go to 17
 16     continue
 17     tdcyuv=tdcyuv+uvintg
        u(i,j,k+nn)=uvintg/(phi-plo)
      end if
 18   continue
c
      if (abs(tdcyuv).gt.acurcy*uvscl*pold(kk+1))
     . write (*,104) i,j,'  hybgen - bad u intgl.',totuv,
     .  tdcyuv,tdcyuv/(uvscl*pold(kk+1))
 14   continue
c
c --- integrate -v- over new depth intervals
c
      do 24 l=1,isv(j)
      do 24 i=ifv(j,l),ilv(j,l)
c
      pold(1)=0.
      pnew(1)=0.
      totuv=0.
c
      do 25 k=1,kk
      kn=k+nn
      vold(k)=v(i,j,kn)
      pold(k+1)=pv(i,j,k+1)
      pnew(k+1)=pnew(k)+dpv(i,j,kn)
 25   totuv=totuv+vold(k)*(pold(k+1)-pold(k))
      tdcyuv=-totuv
c
      do 28 k=1,kk
      phi=pnew(k+1)
      plo=pnew(k  )
      if (phi.gt.plo) then
        uvintg=0.
        pb=plo
        do 26 ko=1,kk
        if (pold(ko+1).le.plo) go to 26
        pa=pb
        pb=min(phi,pold(ko+1))
        uvintg=uvintg+vold(ko)*(pb-pa)
        if (pa.ge.phi) go to 27
 26     continue
 27     tdcyuv=tdcyuv+uvintg
        v(i,j,k+nn)=uvintg/(phi-plo)
      end if
 28   continue
c
      if (abs(tdcyuv).gt.acurcy*uvscl*pold(kk+1))
     . write (*,104) i,j,'  hybgen - bad v intgl.',totuv,
     .  tdcyuv,tdcyuv/(uvscl*pold(kk+1))
 24   continue
c
 13   continue
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      do 33 j=J_0,J_1
      do 33 l=1,isp(j)
      do 33 i=ifp(j,l),ilp(j,l)
      vrbos=i.eq.itest .and. j.eq.jtest
      if (vrbos) write (*,103) nstep,i,j,
     .  '  exiting  hybgen:  temp    saln    dens    thkns    dpth',
     .  (k,temp(i,j,k+nn),saln(i,j,k+nn),
     .  th3d(i,j,k+nn),dp(i,j,k+nn)/onem,
     .  p(i,j,k+1)/onem,k=1,kk)
      if (vrbos) write (*,106) nstep,i,j,
     .  '  exiting  hybgen:  dpthu      u    dpthv      v',
     .  (k,pu(i,j,k+1)/onem,u(i,j,k+nn),
     .     pv(i,j,k+1)/onem,v(i,j,k+nn),k=1,kk)
 33   continue
c
      if (mod(time+.0001,1.).lt..0002) then
!nwrk2=0
!nwrk3=0
!ntot2=0
!ntot3=0
!do j=1,JDM
!nwrk2=nwrk2+nwrk2d(j)
!nwrk3=nwrk3+nwrk3d(j)
!ntot2=ntot2+ntot2d(j)
!ntot3=ntot3+ntot3d(j)
!end do
        anwrkd(:) = nwrk2d(:)
        call GLOBALSUM(ogrid,anwrkd,anwrk)
        if( AM_I_ROOT() ) nwrk2 = anwrk+0.1
        anwrkd(:) = nwrk3d(:)
        call GLOBALSUM(ogrid,anwrkd,anwrk)
        if( AM_I_ROOT() ) nwrk3 = anwrk+0.1
        anwrkd(:) = ntot2d(:)
        call GLOBALSUM(ogrid,anwrkd,anwrk)
        if( AM_I_ROOT() ) ntot2 = anwrk+0.1
        anwrkd(:) = ntot3d(:)
        call GLOBALSUM(ogrid,anwrkd,anwrk)
        if( AM_I_ROOT() ) ntot3 = anwrk+0.1

        if( AM_I_ROOT() ) then
        write (*,'(a,f6.1,a,i9,a)') 'hybgen - grid restoration at',
     .   100.*float(nwrk3)/float(ntot3),' per cent of',ntot3,' points'
        write (*,'(a,f6.1,a,i9,a)') 'hybgen - new bottom layer at',
     .   100.*float(nwrk2)/float(ntot2),' per cent of',ntot2,' points'
        end if ! AM_I_ROOT
      end if

      return
      end
c
c
      subroutine ppmad3(x,dx,y,ynew,diagno)
c
c --- advection by piecewise parabolic method
c --- this is a special version taylored to nmax = 3
c                                           ^^^^^^^^
c --- input variables:
c --- y(nmax)    - function values at cell midpoints
c --- x(nmax+1)  - cell boundaries
c --- dx(nmax+1) - displacement of cell boundaries during 1 time step
c
c --- output variables:
c --- ynew(nmax) - function values after advection (overwriting of -y- allowed)
c
      implicit none
      integer nmax,n
      parameter (nmax=3)
      real x(nmax+1),dx(nmax+1),y(nmax),ynew(nmax),total,tndcy,scale,
     .     ytmp(nmax),wdth,slab,dxnew,acurcy,a,b,c,athird,onemu,yl,yr
      logical diagno
      data acurcy/1.e-11/,onemu/.098/
      parameter (athird=1./3.)
c
      total=0.
      scale=1.e-99
      do 3 n=1,nmax
      if (x(n).gt.x(n+1)+onemu) then
        write (*,'(a,4f9.0)') 'error: x not monotonic in ppmad3',x
        stop '(ppmad3)'
      end if
c
      ytmp(n)=y(n)*(x(n+1)-x(n))
      ynew(n)=y(n)
      total=total+ytmp(n)
 3    scale=scale+abs(ytmp(n))
      scale=scale/float(nmax)
c
c --- construct parabola whose integral over [-.5,+.5] equals y(2) and
c --- which goes though points [-.5,(y(1)+y(2))/2], [+.5,(y(2)+y(3))/2]
c
      yl=.5*(y(1)+y(2))
      yr=.5*(y(3)+y(2))
      a=1.5*y(2)-.25*(yl+yr)
      b=yr-yl
      c=6.*(.5*(yl+yr)-y(2))
      if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-y(2))) then
c
c --- apex of parabola lies inside interval [-.5,+.5].
c --- => need to change curve to prevent over- and undershoots
c
        if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-y(2))) then
c --- put apex of parabola on edge of interval [-.5,+.5]
          if ((yr-yl)*(.5*(yl+yr)-y(2)) .gt. 0.) then
c --- apex at x=-.5
            a=.25*(3.*y(2)+yl)
            c=3.*(y(2)-yl)
            b=c
          else
c --- apex at x=+.5
            a=.25*(3.*y(2)+yr)
            c=3.*(y(2)-yr)
            b=-c
          end if
        else          !  -1/6 < x < +1/6
c --- can't put apex on edge of interval. only option is to flatten curve
          a=y(2)
          b=0.
          c=0.
        end if
      end if
c
c --- now transport -y- across cell boundaries
c
      if (dx(2).gt.0.) then
c --- velocity at left edge is negative. export slab to cell on left
        wdth= dx(2)/(x(3)-x(2))
        slab=(x(3)-x(2))*
     .   wdth*(a+b*.5*(wdth-1.)+c*(.25-wdth*(.5-wdth*athird)))
        ytmp(2)=ytmp(2)-slab
        ytmp(1)=ytmp(1)+slab
      end if
      if (dx(3).lt.0.) then
c --- velocity at right edge is positive. export slab to cell on right
        wdth=-dx(3)/(x(3)-x(2))
        slab=(x(3)-x(2))*
     .   wdth*(a+b*.5*(1.-wdth)+c*(.25-wdth*(.5-wdth*athird)))
        ytmp(2)=ytmp(2)-slab
        ytmp(3)=ytmp(3)+slab
      end if
c
      if (diagno) write (*,100) 'ppmad3 in: ',x,y
 100  format (a,4f9.0,3f9.4)
c
      tndcy=-total
      do 5 n=1,nmax
      dxnew=x(n+1)-x(n)+dx(n+1)-dx(n)
      if (dxnew.gt.0.) ynew(n)=ytmp(n)/dxnew
 5    tndcy=tndcy+ynew(n)*dxnew
      if (abs(tndcy).gt.acurcy*max(abs(total),scale*x(nmax+1)))
     .  write (*,'(a,1p,2e11.3)') 'ppmad3 - bad intgl.:',total,tndcy
c
      if (diagno) write (*,100) 'ppmad3 out:',dx,ynew
c
      return
      end
c
c
      subroutine ppmadv(nmax,x,dx,y,ynew,diagno)
c
c --- advection by piecewise parabolic method
c --- this version is customized for recursive updating of cell boundaries
c                                    ^^^^^^^^^
c --- input variables:
c --- y(nmax)    - function values at cell midpoints
c --- x(nmax+1)  - cell boundaries
c --- dx(nmax+1) - displacement of cell boundaries during 1 time step
c
c --- output variables:
c --- ynew(nmax) - function values after advection (overwriting of -y- allowed)
c
      implicit none
      integer nmax,n,na,nb,m
      real x(nmax+1),dx(nmax+1),y(nmax),ynew(nmax),
     .     xtmp(nmax+1),ytmp(nmax),yold(nmax),athird,a,b,c,
     .     wdth,slab,dxnew,acurcy,total,tndcy,onemu,scale,yl,yr
      logical diagno
      data acurcy/1.e-12/,onemu/.098/
      parameter (athird=1./3.)
c
      total=0.
      scale=1.e-99
      do 3 n=1,nmax
      if (x(n).gt.x(n+1)+onemu) then
        write (*,'(a/(2(a,i2),a/8f9.0))')
     .   'monotonicity error in ppmadv: x(',n,') > x(',n+1,')',x
        stop '(ppmadv)'
      end if
c
      xtmp(n)=x(n)
      ytmp(n)=y(n)*(x(n+1)-x(n))
      yold(n)=y(n)
      ynew(n)=y(n)
      total=total+ytmp(n)
 3    scale=scale+abs(ytmp(n))
      scale=scale/float(nmax)
      xtmp(nmax+1)=x(nmax+1)
c
      do 4 n=1,nmax
      na=max(   1,n-1)			!  non-cyclic
      nb=min(nmax,n+1)			!  non-cyclic
      if (xtmp(n)+dx(n).lt.xtmp(na )-onemu .or.
     .    xtmp(n)+dx(n).gt.xtmp(n+1)+onemu) then
        write (*,'(a,i3/(i3,3f15.2))')
     .   'ppmadv error: x+dx out of bounds, n =',
     .    n,(m,x(m),dx(m),xtmp(m),m=1,n-1),
     .    n,x(n),dx(n),xtmp(n)+dx(n),n+1,x(n+1)
        stop '(ppmadv)'
      end if
c
      if (xtmp(n+1).gt.xtmp(n)) then
c
c --- construct parabola whose integral over [-.5,+.5] equals y(n) and
c --- which goes though points [-.5,(y(na)+y(n))/2], [+.5,(y(n)+y(nb))/2]
c
      yl=.5*(y(na)+y(n))
      yr=.5*(y(nb)+y(n))
      a=1.5*y(n)-.25*(yl+yr)
      b=yr-yl
      c=6.*(.5*(yl+yr)-y(n))
      if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-y(n))) then
c
c --- apex of parabola lies inside interval [-.5,+.5].
c --- => need to change curve to prevent over- and undershoots
c
        if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-y(n))) then
c --- put apex of parabola on edge of interval [-.5,+.5]
          if ((yr-yl)*(.5*(yl+yr)-y(n)) .gt. 0.) then
c --- apex at x=-.5
            a=.25*(3.*y(n)+yl)
            c=3.*(y(n)-yl)
            b=c
          else
c --- apex at x=+.5
            a=.25*(3.*y(n)+yr)
            c=3.*(y(n)-yr)
            b=-c
          end if
        else			!  -1/6 < x < +1/6
c --- can't put apex on edge of interval. only option is to flatten curve
          a=y(n)
          b=0.
          c=0.
        end if
      end if
c
c --- now transport -y- across cell boundaries
c
      if (dx(n).gt.0.) then
c --- velocity at left edge is negative. export slab to cell on left
        wdth=min(1.,dx(n)/(xtmp(n+1)-xtmp(n)))
        slab=(xtmp(n+1)-xtmp(n))*
     .   wdth*(a+b*.5*(wdth-1.)+c*(.25-wdth*(.5-wdth*athird)))
        ytmp(n )=ytmp(n )-slab
        ytmp(na)=ytmp(na)+slab
        end if
        if (dx(n+1).lt.0.) then
c --- velocity at right edge is positive. export slab to cell on right
        wdth=min(1.,-dx(n+1)/(xtmp(n+1)-xtmp(n)))
        slab=(xtmp(n+1)-xtmp(n))*
     .   wdth*(a+b*.5*(1.-wdth)+c*(.25-wdth*(.5-wdth*athird)))
        ytmp(n )=ytmp(n )-slab
        ytmp(nb)=ytmp(nb)+slab
        end if
      end if
c
c --- before proceding to next cell, update cell boundary
 4    xtmp(n)=xtmp(n)+dx(n)
c
      tndcy=-total
      do 5 n=1,nmax
      dxnew=xtmp(n+1)-xtmp(n)
      if (dxnew.gt.0.) ynew(n)=ytmp(n)/dxnew
 5    tndcy=tndcy+ynew(n)*dxnew
      if (abs(tndcy).gt.acurcy*max(abs(total),scale*x(nmax+1))) then
        write (*,'(a,1p,2e11.3)') 'ppmadv - bad intgl.:',total,tndcy
        write (*,100) 'ppmadv:',
     .   (n,x(n),yold(n),dx(n),xtmp(n),ynew(n),n=1,nmax)
      end if
c
      if (diagno) write (*,100) 'ppmadv:',
     .   (n,x(n),yold(n),dx(n),xtmp(n),ynew(n),n=1,nmax)
 100  format (a/(i3,0pf14.1,1pe17.7,0p,2f14.1,1pe17.7))
c
      return
      end
c
c
      real function cushn(delp,dp0)
c
c --- c u s h i o n   function (from Bleck & Benjamin, 1993):
c
c              /  1                       for x < 2-x1
c              |
c              |      (x + x1 - 2)**2
c --- cushn = <   1 + ---------------     for 2-x1 < x < x1
c              |        4 ( x1 - 1)
c              |
c              \  x                       for x > x1
c
c --- if x = delp/dp0 >>  0, cushn*dp0 returns -delp-
c --- if x = delp/dp0 <<  0, cushn*dp0 returns -dp0-
      implicit none
c
      real delp,dp0,qq,x1,factor
ccc   parameter (x1=4.)                 !  used in Bleck&Benjamin 1993
ccc   parameter (x1=6.)                 !  used in Bleck 2002
ccc   parameter (x1=8.)
      parameter (x1=6.828427125)        !  x1=4+sqrt(8) yields cushn(0)=2
c
      parameter (factor=.25/(x1-1.))
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc   qq=max(-1.,min(2.,delp/(2.*dp0)))
ccc   cushn=dp0*(1.+athird*(qq+1.)**2)
ccc  .            *max(1.,delp/(2.*dp0))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc   qq=max(-.75,min(1.25,delp/(4.*dp0)))
ccc   cushn=dp0*(1.+(qq+.75)**2)
ccc  .            *max(1.,delp/(5.*dp0))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc   qq=max(-1.,min(1.5,delp/(4.*dp0)))
ccc   cushn=dp0*(1+.8*(qq+1.)**2)
ccc  .            *max(1.,delp/(6.*dp0))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      qq=max(2.-x1,min(x1,delp/dp0))
      cushn=dp0*(1.+factor*(qq+x1-2.)**2) * max(1.,delp/(x1*dp0))
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
c> Jan. 2007 - disallow dens(k-1) < theta(k-1) in upper sublayer

