#include "hycom_mpi_hacks.h"
      subroutine cnuity(m,n,mm,nn,k1m,k1n)
c
c --- micom version 2.9
c --- hycom version 0.9
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE,NORTH,SOUTH,GLOBALMAX
      USE HYCOM_DIM, only : ii,jj,isp,ifp,ilp,isu,ifu,ilu,isv,ifv,ilv,kk
     &     ,ip,idm,ii1,JDM
     &     ,jchunk
      USE HYCOM_DIM, only : ogrid,J_0,J_1,J_0H,J_1H
      USE HYCOM_SCALARS, only : acurcy,nstep,delt1,onecm,epsil,thkdff
     &     ,sigjmp,onem,bolus_biharm_constant,bolus_laplc_constant
     &     ,bolus_laplc_exponential
      USE HYCOM_ARRAYS, only : utotn,vtotn,dp,utotm,u,ubavg,scuy,depthu
     &     ,uflux,uflux2,dpu,uflx,vtotm,v,scvx,depthv,vflux
     &     ,vflux2,dpv,vflx,dpold,scp2i,vbavg,util1,util2,scp2,p
     &     ,util3,pbot,diaflx,bolusu,bolusv
      implicit none
c
      integer i,j,k,l,m,n,mm,nn,km,kn,k1m,k1n,ia,ib,ja,jb
c
      integer mask(idm,J_0H:J_1H),jcyc,iz,jz,iter
      real pold(idm,J_0H:J_1H),q,dpmin,dpmn(J_0H:J_1H),clip,flxhi,flxlo,
     &     dtinv,old,dpmin_loc,
     .     hymar,boost,pbot1,pbot2,p1,p2,hyc_pechg1,hyc_pechg2,pa,pb
      external hyc_pechg1,hyc_pechg2
      character text*20
      logical abort

      integer, parameter :: itmax = 15
      !!integer ja_,jb_
cddd      integer my_pet
c
c --- boost = 1 to within 15 m of bottom, then increases linearly to 1.5
c     boost(pbot1,pbot2,p1,p2)=max(1.,1.5-min(pbot1-p1,pbot2-p2)
c    .  /(30.*onem))
c --- boost = 1 to within 30 m of bottom, then increases linearly to 2
c     boost(pbot1,pbot2,p1,p2)=max(1.,2.-min(pbot1-p1,pbot2-p2)
c    .  /(30.*onem))
c --- boost = 1 to within 100 m of bottom, then increases linearly to 3
      boost(pbot1,pbot2,p1,p2)=max(1.,3.-min(pbot1-p1,pbot2-p2)
     .  /(50.*onem))
c
c --- ------------------------------------------------------
c --- continuity equation (flux-corrected transport version)
c --- ------------------------------------------------------
c
      do 41 j=J_0,J_1
c
      !write(0,*) my_pet, j, isp(j)
      do 74 l=1,isp(j)
      do 74 i=ifp(j,l),ilp(j,l)
 74   pold(i,j)=0.
c
      do 40 l=1,isu(j)
      do 40 i=ifu(j,l),ilu(j,l)
 40   utotn(i,j)=0.
c
      do 41 l=1,isv(j)
      do 41 i=ifv(j,l),ilv(j,l)
 41   vtotn(i,j)=0.
c
      do 76 k=1,kk
      km=k+mm
      kn=k+nn
c
c --- uflux/vflux = low-order (diffusive) mass fluxes at old time level.
c --- uflux2/vflux2 = 'antidiffusive' fluxes, defined as high-order minus low-
c --- order fluxes. high-order fluxes are second-order in space, time-centered.
c
      !write(0,*) "ok ",__FILE__,__LINE__
      call cpy_p_par(dp(:,:,kn))
      !write(0,*) "ok ",__FILE__,__LINE__

      call cpy_p_par(pold)
c

      !write(0,*) "ok ",__FILE__,__LINE__
      CALL HALO_UPDATE(ogrid,dp(:,:,kn),SOUTH)
      CALL HALO_UPDATE(ogrid,pold,SOUTH)

!!      do 12 j=1,jj
      do 12 j=J_0,J_1
!!      ja=mod(j-2+jj,jj)+1
        !ja = j-1
        ja = PERIODIC_INDEX(j-1, jj)
c
      do 11 l=1,isu(j)
      do 11 i=ifu(j,l),ilu(j,l)
      utotm(i,j)=(u(i,j,km)+ubavg(i,j,m))*scuy(i,j)
      if (utotm(i,j).ge.0.) then
        q=min(dp(i-1,j,kn),max(0.,depthu(i,j)-pold(i-1,j)))
      else
        q=min(dp(i  ,j,kn),max(0.,depthu(i,j)-pold(i  ,j)))
      end if
      uflux(i,j)=utotm(i,j)*q
      uflux2(i,j)=utotm(i,j)*dpu(i,j,km)-uflux(i,j)
 11   uflx(i,j,k)=uflux(i,j)
c
      do 12 l=1,isv(j)
      do 12 i=ifv(j,l),ilv(j,l)
      vtotm(i,j)=(v(i,j,km)+vbavg(i,j,m))*scvx(i,j)
      if (vtotm(i,j).ge.0.) then
        q=min(dp(i,ja ,kn),max(0.,depthv(i,j)-pold(i,ja )))
      else
        q=min(dp(i,j  ,kn),max(0.,depthv(i,j)-pold(i,j  )))
      end if
      vflux(i,j)=vtotm(i,j)*q
      vflux2(i,j)=vtotm(i,j)*dpv(i,j,km)-vflux(i,j)
 12   vflx(i,j,k)=vflux(i,j)
c
cddd      if (beropn .and.
cddd     .    abs(uflx(ipacs,jpac,k)+uflx(iatln,jatl,k)).gt.
cddd     .max(abs(uflx(ipacs,jpac,k)-uflx(iatln,jatl,k))*acurcy,1.))
cddd     .  write(*,104) nstep,
cddd     . ' cnuity1 WRONG uflx k=',k,uflx(ipacs,jpac,k),uflx(iatln,jatl,k)
 104  format (i9,a,i2,2es15.7)

c
c --- advance -dp- field using low-order (diffusive) flux values
c
      CALL HALO_UPDATE(ogrid,vflux,NORTH)
      !!do 19 j=1,jj
      do 19 j=J_0,J_1
      !!jb=mod(j     ,jj)+1
        !jb = j+1
        jb = PERIODIC_INDEX(j+1, jj)
      dpmn(j)=999.
      do 19 l=1,isp(j)
      do 19 i=ifp(j,l),ilp(j,l)
      dpold(i,j,k)=dp(i,j,kn)
      pold(i,j)=pold(i,j)+dp(i,j,kn)
      dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,jb )-vflux(i,j))*delt1*scp2i(i,j)
 19   dpmn(j)=min(dpmn(j),dp(i,j,kn))
c
ccc      do j=1,jdm
ccc      do i=1,idm
ccc      mask(i,j)=iu(i,j)
ccc      if (i.gt. 1 ) mask(i,j)=mask(i,j)+iu(i-1,j)
ccc      if (i.lt.idm) mask(i,j)=mask(i,j)+iu(i+1,j)
ccc      end do
ccc      end do
ccc      write (text,'(a9,i3,i8)') 'uflux  k=',k,nstep
ccc      call compare(uflux,mask,text)
ccc      do j=1,jdm
ccc      do i=1,idm
ccc      mask(i,j)=iv(i,j)
ccc      if (j.gt. 1 ) mask(i,j)=mask(i,j)+iv(i,ja )
ccc      if (j.lt.jdm) mask(i,j)=mask(i,j)+iv(i,jb )
ccc      end do
ccc      end do
ccc      write (text,'(a9,i3,i8)') 'vflux  k=',k,nstep
ccc      call compare(vflux,mask,text)
c
      dpmin=999.
      do 191 j=J_0,J_1
 191  dpmin=min(dpmin,dpmn(j))
c
      if (dpmin.lt.-onem) then
      do 190 j=J_0,J_1
        do 190 l=1,isp(j)
        do 190 i=ifp(j,l),ilp(j,l)
        if (dp(i,j,kn).eq.dpmin) then
          write (*,100) nstep,i,j,k,19,dpmin/onem
 100      format (i9,' i,j,k=',2i5,i3,' neg. dp (m) in loop ',i3,f9.2)
          iz=i
          jz=j
        end if
 190    continue
        call stencl(iz,jz,k,nn)
      end if
c
cdiag write (*,*) 'time step',nstep,'    layer',k
cdiag do jcyc=jtest-1,jtest+1
cdiag j =mod(jcyc-1+jj,jj)+1
cdiag ja=mod(jcyc-2+jj,jj)+1
cdiag jb=mod(jcyc     ,jj)+1
cdiag do i=itest-1,itest+1
cdiag write (*,101) i,j,k,'old thknss','mid thknss,vel.',
cdiag.'new thknss,fluxes',
cdiag.dpold(i-1,j,k)/onem,u(i,j,km)+ubavg(i,j,m),uflux(i,j),
cdiag.dpold(i,ja,k)/onem,dpold(i,j,k)/onem,dpold(i,jb,k)/onem,
cdiag.v(i,j,km)+vbavg(i,j,m),dp(i,j,km)/onem,v(i,jb,km)+vbavg(i,jb,m),
cdiag.vflux(i,j),dp(i,j,kn)/onem,vflux(i,jb),
cdiag.dpold(i+1,j,k)/onem,u(i,jb,km)+ubavg(i,jb,m),uflux(i+1,j)
cdiag end do
cdiag end do
 101  format (2i5,i3,2x,a,7x,a,10x,a/0p,f14.1,f26.1,1p,e30.2/
     .   0p,3f7.1,f12.1,f7.1,f6.1,1p,e15.2,0p,f8.1,1p,e10.2/
     .   0p,f14.1,f26.1,1p,e30.2)
c
c --- at each grid point, determine the ratio of the largest permissible
c --- pos. (neg.) change in -dp- to the sum of all incoming (outgoing) fluxes
c
      call cpy_p_par(dp(:,:,kn))
c
      CALL HALO_UPDATE(ogrid,dp(:,:,kn))
      CALL HALO_UPDATE(ogrid,vflux2,NORTH)


      !write(922,*) "bounds=",isp,ifp,ilp,ip,j_0,j_1
      !write(910+my_pet,*) "bounds=",j_0,j_1,kn,delt1
      !!do 26 j=1,jj
      do 26 j=J_0,J_1
      do 26 l=1,isp(j)
      do 26 i=ifp(j,l),ilp(j,l)
      ia=max( 1,i-1)
      if (ip(ia,j).eq.0) ia=i
      ib=min(ii,i+1)
      if (ip(ib,j).eq.0) ib=i
      !ja_=mod(j-2+jj,jj)+1
      !ja = j-1
      !ja_ = ja
      ja = PERIODIC_INDEX(j-1, jj)
      if (ip(i,ja).eq.0) ja=j
      !jb_=mod(j     ,jj)+1
      !jb = j+1
      !jb_ = jb
      jb = PERIODIC_INDEX(j+1, jj)
      if (ip(i,jb).eq.0) jb=j
      util1(i,j)=max(dp(i,j,kn),dp(ia,j,kn),dp(ib,j,kn),
     .                          dp(i,ja,kn),dp(i,jb,kn))
      util2(i,j)=max(0.,
     .           min(dp(i,j,kn),dp(ia,j,kn),dp(ib,j,kn),
     .                          dp(i,ja,kn),dp(i,jb,kn)))
c
      !!jb=mod(j     ,jj)+1
      !jb = j+1
      jb = PERIODIC_INDEX(j+1, jj)
      util1(i,j)=(util1(i,j)-dp(i,j,kn))*scp2(i,j)
     ./((max(0.,uflux2(i,j))-min(0.,uflux2(i+1,j))
     .  +max(0.,vflux2(i,j))-min(0.,vflux2(i,jb ))+epsil)*delt1)
c
      util2(i,j)=(util2(i,j)-dp(i,j,kn))*scp2(i,j)
     ./((min(0.,uflux2(i,j))-max(0.,uflux2(i+1,j))
     .  +min(0.,vflux2(i,j))-max(0.,vflux2(i,jb ))-epsil)*delt1)
 26    continue
c
c --- limit antidiffusive fluxes
c --- (keep track in -utotn,vtotn- of discrepancy between high-order
c --- fluxes and the sum of low-order and clipped antidiffusive fluxes.
c --- this will be used later to restore nondivergence of barotropic flow)

      call cpy_p_par(util1)
      call cpy_p_par(util2)
c
      CALL HALO_UPDATE(ogrid,util1,SOUTH)
      CALL HALO_UPDATE(ogrid,util2,SOUTH)

      do 29 j=J_0,j_1
      !!do 29 j=1,jj
      !!ja=mod(j-2+jj,jj)+1
        !ja = j-1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 28 l=1,isu(j)
      do 28 i=ifu(j,l),ilu(j,l)
      if (uflux2(i,j).ge.0.) then
        clip=min(1.,util1(i,j),util2(i-1,j))
      else
        clip=min(1.,util2(i,j),util1(i-1,j))
      end if
      utotn(i,j)=utotn(i,j)+uflux2(i,j)*(1.-clip)
      uflux(i,j)=uflux2(i,j)*clip
 28   uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
c
      do 29 l=1,isv(j)
      do 29 i=ifv(j,l),ilv(j,l)
      if (vflux2(i,j).ge.0.) then
        clip=min(1.,util1(i,j),util2(i,ja ))
      else
        clip=min(1.,util2(i,j),util1(i,ja ))
      end if
      vtotn(i,j)=vtotn(i,j)+vflux2(i,j)*(1.-clip)
      vflux(i,j)=vflux2(i,j)*clip
 29   vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
c
cddd      if (beropn .and.
cddd     .    abs(uflx(ipacs,jpac,k)+uflx(iatln,jatl,k)).gt.
cddd     .max(abs(uflx(ipacs,jpac,k)-uflx(iatln,jatl,k))*acurcy,1.))
cddd     .  write(*,104) nstep,
cddd     . ' cnuity2 WRONG uflx k=',k,uflx(ipacs,jpac,k),uflx(iatln,jatl,k)
c
c --- evaluate effect of antidiffusive fluxes on -dp- field
c
      CALL HALO_UPDATE(ogrid,vflux,NORTH)
      !!do 15 j=1,jj
      do 15 j=J_0,J_1
      !!jb=mod(j     ,jj)+1
        !jb = j+1
        jb = PERIODIC_INDEX(j+1, jj)
      dpmn(j)=999.
      do 15  l=1,isp(j)
      do 15  i=ifp(j,l),ilp(j,l)
      dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,jb )-vflux(i,j))*delt1*scp2i(i,j)
      p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
 15   dpmn(j)=min(dpmn(j),dp(i,j,kn))
c
      dpmin=999.
      do 149 j=J_0,J_1
 149  dpmin=min(dpmin,dpmn(j))
c
      if (dpmin.lt.-onem) then
      do 150 j=J_0,J_1
      do 150 l=1,isp(j)
      do 150 i=ifp(j,l),ilp(j,l)
      if (dp(i,j,kn).eq.dpmin) write (*,100) nstep,i,j,k,15,dpmin/onem
 150  continue
      end if
c
 76   continue
c
c --- restore nondivergence of vertically integrated mass flow by
c --- recovering fluxes lost in the flux limiting process.
c --- treat these fluxes as an 'upstream' barotropic correction to
c --- the sum of diffusive and antidiffusive fluxes obtained so far.
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!

c
      call cpy_p_par(p(:,:,kk+1))
      CALL HALO_UPDATE(ogrid,p(:,:,kk+1),SOUTH)
c
      do 77 k=1,kk
      km=k+mm
      kn=k+nn
c
      call cpy_p_par(dp(:,:,kn))
c
      CALL HALO_UPDATE(ogrid,dp(:,:,kn),SOUTH)
      do 45 j=J_0,J_1
      !!do 45 j=1,jj
      !!ja=mod(j-2+jj,jj)+1
        !ja = j-1
        ja = PERIODIC_INDEX(j-1, jj)
c
      do 44 l=1,isu(j)
      do 44 i=ifu(j,l),ilu(j,l)
      if (utotn(i,j).ge.0.) then
        q=dp(i-1,j,kn)/p(i-1,j,kk+1)
      else
        q=dp(i  ,j,kn)/p(i  ,j,kk+1)
      end if
      uflux(i,j)=utotn(i,j)*q
 44   uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)
c
      do 45 l=1,isv(j)
      do 45 i=ifv(j,l),ilv(j,l)
      if (vtotn(i,j).ge.0.) then
        q=dp(i,ja ,kn)/p(i,ja ,kk+1)
      else
        q=dp(i,j  ,kn)/p(i,j  ,kk+1)
      end if
      vflux(i,j)=vtotn(i,j)*q
 45   vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)
c
      CALL HALO_UPDATE(ogrid,vflux,NORTH)
      !!do 14 j=1,jj
      do 14 j=J_0,J_1
      !!jb=mod(j     ,jj)+1
        !jb = j+1
        jb = PERIODIC_INDEX(j+1, jj)
      dpmn(j)=999.
      do 14 l=1,isp(j)
      do 14 i=ifp(j,l),ilp(j,l)
      dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,jb )-vflux(i,j))*delt1*scp2i(i,j)
      p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
 14   dpmn(j)=min(dpmn(j),dp(i,j,kn))
c
      dpmin=999.
      do 139 j=J_0,J_1
 139  dpmin=min(dpmin,dpmn(j))
c
      if (dpmin.lt.-onem) then
      do 140 j=J_0,J_1
      do 140 l=1,isp(j)
      do 140 i=ifp(j,l),ilp(j,l)
      if (dp(i,j,kn).eq.dpmin) write (*,100) nstep,i,j,k,14,dpmin/onem
 140  continue
      end if
c
cddd      if (beropn .and.
cddd     .    abs(uflx(ipacs,jpac,k)+uflx(iatln,jatl,k)).gt.
cddd     .max(abs(uflx(ipacs,jpac,k)-uflx(iatln,jatl,k))*acurcy,1.))
cddd     .  write(*,104) nstep,
cddd     . ' cnuity3 WRONG uflx k=',k,uflx(ipacs,jpac,k),uflx(iatln,jatl,k)
c
 77   continue
c
c --- add bottom-pressure restoring term arising from split-explicit treatment
c --- of continuity equation (step 4 in appendix B of 1992 BRHS paper).
c --- iterate a few times (mimicking an iterative poisson solver) to spawn
c --- mass fluxes approximately consistent with required column stretching
c
      dtinv=1./delt1
      do 18 iter=1,itmax
c
      do 36 j=J_0,J_1
      dpmn(j)=0.
      do 36 l=1,isp(j)
      do 36 i=ifp(j,l),ilp(j,l)
      util3(i,j)=(pbot(i,j)-p(i,j,kk+1))*scp2(i,j)
      util1(i,j)=1./p(i,j,kk+1)
 36   dpmn(j)=max(dpmn(j),abs(util3(i,j)))
c
      if (iter.eq.itmax) then
        dpmin=0.
      do 37 j=J_0,J_1
 37     dpmin=max(dpmin,dpmn(j))
c
      dpmin_loc = dpmin
      call globalmax(ogrid,dpmin_loc,dpmin)
      do 38 j=J_0,J_1
        do 38 l=1,isp(j)
        do 38 i=ifp(j,l),ilp(j,l)
        if (abs(util3(i,j)).eq.dpmin) write (*,105)
     .   nstep,i,j,'  largest pbot correction after',itmax,
     .    ' iterations:',scp2i(i,j)*util3(i,j)/onecm,' cm'
 105    format (i9,2i5,a,i3,a,f7.1,a)
 38     continue
      end if				!  last iteration
c
      call cpy_p_par(util1)
      call cpy_p_par(util3)
c
      do 20 k=1,kk
      kn=k+nn
c
      call cpy_p_par(dp(:,:,kn))
c
      CALL HALO_UPDATE(ogrid,dp(:,:,kn),SOUTH)
      CALL HALO_UPDATE(ogrid,util1,SOUTH)
      CALL HALO_UPDATE(ogrid,util3,SOUTH)
      !!do 21 j=1,jj
      do 21 j=J_0,J_1
      !!ja=mod(j-2+jj,jj)+1
        !ja = j-1
        ja = PERIODIC_INDEX(j-1, jj)
c
      do 22 l=1,isu(j)
      do 22 i=ifu(j,l),ilu(j,l)
      uflux(i,j)=.25*(util3(i,j)-util3(i-1,j))
      if (uflux(i,j).ge.0.) then
        q=dp(i-1,j,kn)*util1(i-1,j)
      else
        q=dp(i  ,j,kn)*util1(i  ,j)
      end if
      uflux(i,j)=uflux(i,j)*q
 22   uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)*dtinv
c
      do 21 l=1,isv(j)
      do 21 i=ifv(j,l),ilv(j,l)
      vflux(i,j)=.25*(util3(i,j)-util3(i,ja ))
      if (vflux(i,j).ge.0.) then
        q=dp(i,ja ,kn)*util1(i,ja )
      else
        q=dp(i,j  ,kn)*util1(i,j  )
      end if
      vflux(i,j)=vflux(i,j)*q
 21   vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)*dtinv
c
      CALL HALO_UPDATE(ogrid,vflux,NORTH)
      !!do 23 j=1,jj
      do 23 j=J_0,J_1
      !!jb=mod(j     ,jj)+1
        !jb = j+1
        jb = PERIODIC_INDEX(j+1, jj)
      do 23 l=1,isp(j)
      do 23 i=ifp(j,l),ilp(j,l)
      dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,jb )-vflux(i,j))*scp2i(i,j)
 23   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
cddd      if (beropn .and.
cddd     .    abs(uflx(ipacs,jpac,k)+uflx(iatln,jatl,k)).gt.
cddd     .max(abs(uflx(ipacs,jpac,k)-uflx(iatln,jatl,k))*acurcy,1.))
cddd     .  write(*,104) nstep,
cddd     . ' cnuity4 WRONG uflx k=',k,uflx(ipacs,jpac,k),uflx(iatln,jatl,k)
c
 20   continue				!  k loop
 18   continue				!  iter
c
c --- now stretch the grid to take care of the residual pbot discrepancy
c
      do 39 j=J_0,J_1
      do 39 l=1,isp(j)
      do 39 k=1,kk
      kn=k+nn
      do 39 i=ifp(j,l),ilp(j,l)
      old=dp(i,j,kn)
      dp(i,j,kn)=dp(i,j,kn)*pbot(i,j)/p(i,j,kk+1)
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,kn)-old)		! diapyc.flux
 39   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
c --- -----------------------------------
c --- interface smoothing (=> bolus flux)
c --- -----------------------------------
c
      if (thkdff.eq.0.) return
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     if (diagno)
c    .  q=hyc_pechg1(dp(1,1,k1n),th3d(1,1,k1n),32)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      hymar=1./(10.*sigjmp)
c
      do 8 j=J_0,J_1
      do 9 l=1,isu(j)
      do 9 i=ifu(j,l),ilu(j,l)
      utotn (i,j  )=0.
 9    bolusu(i,j,1)=0.
      do 8 l=1,isv(j)
      do 8 i=ifv(j,l),ilv(j,l)
      vtotn (i,j  )=0.
 8    bolusv(i,j,1)=0.
c
      do 13 k=2,kk
      km=k+mm
c
      call cpy_p_par(p(:,:,k))
c
      CALL HALO_UPDATE(ogrid,p(:,:,k))
      CALL HALO_UPDATE(ogrid,pbot)
      !!do 131 j=1,jj
      do 131 j=J_0,J_1
      !!ja=mod(j-2+jj,jj)+1
      !ja_=mod(j-2+jj,jj)+1
        !ja = j-1
        !ja_ = ja
        ja = PERIODIC_INDEX(j-1, jj)
      !!jb=mod(j     ,jj)+1
      !jb_=mod(j     ,jj)+1
      !jb = j+1
      !jb_ = jb
        jb = PERIODIC_INDEX(j+1, jj)
      do 131 l=1,isp(j)
      do 131 i=ifp(j,l),ilp(j,l)
      ia=max(  1,i-1)
      ib=min(ii1,i+1)
      util1(i,j)=p(i,j,k)
      util2(i,j)=p(i,j,k)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- in preparation for biharmonic smoothing, compute 2nd derivatives
c --- of pressure in x,y direction. store in -util1,util2- resp.
c --- don't allow grid points elevated by topography to affect derivatives.
c --- disable next 22 lines to switch from biharm. to laplacian smoothing
      if (bolus_laplc_exponential==0 .and. bolus_laplc_constant==0) then ! use biharm
        pa=p(ia,j,k)
        pb=p(ib,j,k)
        if   (float(ip(ia,j))*pbot(ia,j).lt.p(i,j,k)) then
          pa=p(i,j,k)
          if (float(ip(ib,j))*pbot(ib,j).gt.p(i,j,k)) pa=p(ib,j,k)
        end if
        if   (float(ip(ib,j))*pbot(ib,j).lt.p(i,j,k)) then
          pb=p(i,j,k)
          if (float(ip(ia,j))*pbot(ia,j).gt.p(i,j,k)) pb=p(ia,j,k)
        end if
        util1(i,j)=util1(i,j)-.5*(pa+pb)
        pa=p(i,ja,k)
        pb=p(i,jb,k)
        if   (float(ip(i,ja))*pbot(i,ja).lt.p(i,j,k)) then
          pa=p(i,j,k)
          if (float(ip(i,jb))*pbot(i,jb).gt.p(i,j,k)) pa=p(i,jb,k)
        end if
        if   (float(ip(i,jb))*pbot(i,jb).lt.p(i,j,k)) then
          pb=p(i,j,k)
          if (float(ip(i,ja))*pbot(i,ja).gt.p(i,j,k)) pb=p(i,ja,k)
        end if
        util2(i,j)=util2(i,j)-.5*(pa+pb)
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 131  continue
c
      call cpy_p_par(util1)
      call cpy_p_par(util2)
c
      CALL HALO_UPDATE(ogrid,util2,SOUTH)
      CALL HALO_UPDATE(ogrid,pbot,SOUTH)
      CALL HALO_UPDATE(ogrid,p(:,:,k),SOUTH)
      CALL HALO_UPDATE(ogrid,scp2,SOUTH)
!!! this omp instruction was removed is it correct ? c$OMP. SHARED(k)
      !!do 16 j=1,jj
      do 16 j=J_0,J_1
      !!ja=mod(j-2+jj,jj)+1
        !ja = j-1
        ja = PERIODIC_INDEX(j-1, jj)
      !!jb=mod(j     ,jj)+1
      !!  jb = j+1
c
      do 151 l=1,isu(j)
      do 151 i=ifu(j,l),ilu(j,l)
      if (bolus_laplc_constant==1 .or. bolus_biharm_constant==1) then
        bolusu(i,j,k)=delt1*thkdff*(util1(i,j)-util1(i-1,j))*scuy(i,j) ! thkdff=const
     .   *boost(pbot(i,j),pbot(i-1,j),p(i,j,k),p(i-1,j,k))
      else if (bolus_laplc_exponential==1) then
        bolusu(i,j,k)=delt1*
     .  (thkdff*exp(-(p(i,j,k)+p(i-1,j,k))/(onem*600.))+0.003) !thkdff=3*exp(-z/300m)+0.3 (cm/s)
     .  *(util1(i,j)-util1(i-1,j))*scuy(i,j)
     .   *boost(pbot(i,j),pbot(i-1,j),p(i,j,k),p(i-1,j,k))
      end if
c
c --- suppress uphill fluxes into 'perched' cells
      if (p(i  ,j,k).gt.pbot(i-1,j)) bolusu(i,j,k)=max(0.,bolusu(i,j,k))
      if (p(i-1,j,k).gt.pbot(i  ,j)) bolusu(i,j,k)=min(0.,bolusu(i,j,k))
c
c --- confine interface smoothing to isopycnic coord. subdomain
ccc      bolusu(i,j,k)=bolusu(i,j,k)*min(1.,max(.1,
ccc     .   2.-(max(th3d(i,j,km-1),th3d(i-1,j,km-1))-theta(k-1))*hymar))
ccc     .   2.-(max(th3d(i,j,km  ),th3d(i-1,j,km  ))-theta(k  ))*hymar))
c
c --- keep interfaces from going underground
      flxhi= .25*(pbot(i-1,j)-p(i-1,j,k))*scp2(i-1,j)
      flxlo=-.25*(pbot(i  ,j)-p(i  ,j,k))*scp2(i  ,j)
      bolusu(i,j,k)=min(flxhi,max(flxlo,bolusu(i,j,k)))
c
c --- difference of 2 interface 'pressure fluxes' becomes thickness flux
 151  bolusu(i,j,k-1)=bolusu(i,j,k-1)-bolusu(i,j,k)
c
      do 141 l=1,isv(j)
      do 141 i=ifv(j,l),ilv(j,l)
      if (bolus_laplc_constant==1 .or. bolus_biharm_constant==1) then
         bolusv(i,j,k)=delt1*thkdff*(util2(i,j)-util2(i,ja ))*scvx(i,j)  ! thkdff=const
     .    *boost(pbot(i,j),pbot(i,ja ),p(i,j,k),p(i,ja ,k))
      else if (bolus_laplc_exponential==1) then
        bolusv(i,j,k)=delt1*
     .  (thkdff*exp(-(p(i,j,k)+p(i,ja,k))/(onem*600.))+0.003) !thkdff=3*exp(-z/300m+0.3 (cm/s)
     .  *(util2(i,j)-util2(i,ja ))*scvx(i,j)
     .   *boost(pbot(i,j),pbot(i,ja ),p(i,j,k),p(i,ja ,k))
      end if
c
c --- suppress uphill fluxes into 'perched' cells
      if (p(i,j  ,k).gt.pbot(i,ja )) bolusv(i,j,k)=max(0.,bolusv(i,j,k))
      if (p(i,ja ,k).gt.pbot(i,j  )) bolusv(i,j,k)=min(0.,bolusv(i,j,k))
c
c --- confine interface smoothing to isopycnic coord. subdomain
ccc      bolusv(i,j,k)=bolusv(i,j,k)*min(1.,max(.1,
ccc     .   2.-(max(th3d(i,j,km-1),th3d(i,ja ,km-1))-theta(k-1))*hymar))
ccc     .   2.-(max(th3d(i,j,km  ),th3d(i,ja ,km  ))-theta(k  ))*hymar))
c
c --- keep interfaces from going underground
      flxhi= .25*(pbot(i,ja )-p(i,ja ,k))*scp2(i,ja )
      flxlo=-.25*(pbot(i,j  )-p(i,j  ,k))*scp2(i,j  )
      bolusv(i,j,k)=min(flxhi,max(flxlo,bolusv(i,j,k)))
c
c --- difference of 2 interface 'pressure fluxes' becomes thickness flux
 141  bolusv(i,j,k-1)=bolusv(i,j,k-1)-bolusv(i,j,k)
 16   continue
c
 13   continue
c
      do 10 k=1,kk
      kn=k+nn
c
c --- at each grid point, determine the ratio of the largest permissible
c --- mass loss to the sum of all outgoing bolus fluxes
c
      CALL HALO_UPDATE(ogrid,bolusv(:,:,k),NORTH)
      !!do 261 j=1,jj
      do 261 j=J_0,J_1
      !!jb=mod(j     ,jj)+1
        !jb = j+1
        jb = PERIODIC_INDEX(j+1, jj)
      do 261 l=1,isp(j)
      do 261 i=ifp(j,l),ilp(j,l)
      util2(i,j)=-dp(i,j,kn)*scp2(i,j)
     ./(min(0.,bolusu(i,j,k))-max(0.,bolusu(i+1,j,k))
     . +min(0.,bolusv(i,j,k))-max(0.,bolusv(i,jb ,k))-epsil)
 261  continue
c
c --- limit bolus fluxes (fct-style)
c
      call cpy_p_par(util2)
c
      CALL HALO_UPDATE(ogrid,util2,SOUTH)
      !!do 291 j=1,jj
      do 291 j=J_0,J_1
      !!ja=mod(j-2+jj,jj)+1
        !ja = j-1
        ja = PERIODIC_INDEX(j-1, jj)
c
      do 281 l=1,isu(j)
      do 281 i=ifu(j,l),ilu(j,l)
      if (bolusu(i,j,k).ge.0.) then
        clip=min(1.,util2(i-1,j))
      else
        clip=min(1.,util2(i  ,j))
      end if
c --- clipped part of bolus flux is kept to restore zero column integral later
      utotn(i,j)=utotn(i,j)+bolusu(i,j,k)*(1.-clip)
      bolusu(i,j,k)=bolusu(i,j,k)*clip
c --- add bolus component to total mass flux
 281  uflx(i,j,k)=uflx(i,j,k)+bolusu(i,j,k)*dtinv
c
      do 291 l=1,isv(j)
      do 291 i=ifv(j,l),ilv(j,l)
      if (bolusv(i,j,k).ge.0.) then
        clip=min(1.,util2(i,ja ))
      else
        clip=min(1.,util2(i,j  ))
      end if
c --- clipped part of bolus flux is kept to restore zero column integral later
      vtotn(i,j)=vtotn(i,j)+bolusv(i,j,k)*(1.-clip)
      bolusv(i,j,k)=bolusv(i,j,k)*clip
c --- add bolus component to total mass flux
 291  vflx(i,j,k)=vflx(i,j,k)+bolusv(i,j,k)*dtinv
c
      CALL HALO_UPDATE(ogrid,bolusv(:,:,k),NORTH)
      !!do 181 j=1,jj
      do 181 j=J_0,J_1
      !!jb=mod(j     ,jj)+1
        !jb = j+1
        jb = PERIODIC_INDEX(j+1, jj)
      do 181 l=1,isp(j)
      do 181 i=ifp(j,l),ilp(j,l)
      dp(i,j,kn)=dp(i,j,kn)-(bolusu(i+1,j,k)-bolusu(i,j,k)
     .                      +bolusv(i,jb ,k)-bolusv(i,j,k))*scp2i(i,j)
 181  p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
 10   continue
c
      call cpy_p_par(p(:,:,kk+1))
      CALL HALO_UPDATE(ogrid,p(:,:,kk+1),SOUTH)
c
      do 7 k=1,kk
      kn=k+nn
c
      call cpy_p_par(dp(:,:,kn))
c
c --- restore zero column integral of bolus fluxes by recovering fluxes
c --- lost in the flux limiting process. treat these as an 'upstream'
c --- barotropic correction to the bolus fluxes.
c
      CALL HALO_UPDATE(ogrid,dp(:,:,kn),SOUTH)
      !!do 145 j=1,jj
      do 145 j=J_0,J_1
      !!ja=mod(j-2+jj,jj)+1
        !ja = j-1
        ja = PERIODIC_INDEX(j-1, jj)
c
      do 144 l=1,isu(j)
      do 144 i=ifu(j,l),ilu(j,l)
      if (utotn(i,j).ge.0.) then
        q=dp(i-1,j,kn)/p(i-1,j,kk+1)
      else
        q=dp(i  ,j,kn)/p(i  ,j,kk+1)
      end if
      uflux(i,j)=utotn(i,j)*q
c --- add correction to total mass flux
 144  uflx(i,j,k)=uflx(i,j,k)+uflux(i,j)*dtinv
c
      do 145 l=1,isv(j)
      do 145 i=ifv(j,l),ilv(j,l)
      if (vtotn(i,j).ge.0.) then
        q=dp(i,ja ,kn)/p(i,ja ,kk+1)
      else
        q=dp(i,j  ,kn)/p(i,j  ,kk+1)
      end if
      vflux(i,j)=vtotn(i,j)*q
c --- add correction to total mass flux
 145  vflx(i,j,k)=vflx(i,j,k)+vflux(i,j)*dtinv
c
      CALL HALO_UPDATE(ogrid,vflux,NORTH)
      !!do 182 j=1,jj
      do 182 j=J_0,J_1
      !!jb=mod(j     ,jj)+1
        !jb = j+1
        jb = PERIODIC_INDEX(j+1, jj)
      do 182 l=1,isp(j)
      do 182 i=ifp(j,l),ilp(j,l)
 182  dp(i,j,kn)=dp(i,j,kn)-(uflux(i+1,j)-uflux(i,j)
     .                      +vflux(i,jb )-vflux(i,j))*scp2i(i,j)
c
cddd      if (beropn .and.
cddd     .    abs(uflx(ipacs,jpac,k)+uflx(iatln,jatl,k)).gt.
cddd     .max(abs(uflx(ipacs,jpac,k)-uflx(iatln,jatl,k))*acurcy,1.))
cddd     .  write(*,104) nstep,
cddd     . ' cnuity5 WRONG uflx k=',k,uflx(ipacs,jpac,k),uflx(iatln,jatl,k)
c
 7    continue
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     if (diagno) then
c       q=hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),32)
c       write (*,103) time,'  APE change due to intfc smoothing:',q
c     end if
 103  format (f9.1,a,-12p,f9.3,' TW')
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



      return
      end
c
c
c> Revision history:
c>
c> July 1997 - combined diff. and antidiff.flux calc. (eliminated loops 20,21)
c> Aug. 1997 - set u/vflux=0 before entering thickness smoothing k-loop 13
c> Jul. 1998 - reversed i/j loop nesting in loop 26
c> Apr. 2000 - changed i/j loop nesting to j/i
c> May  2000 - added code to eliminate neg. dp resulting from intfc.smoothing
c> May  2000 - added provisions to enhance -thkdff- in regions of strong flow
c> May  2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> Oct. 2000 - added limiter to intfc smoothing to maintain finite mxlyr thknss
c> Feb. 2001 - further enhanced -thkdff- enhancer by adding -glue-
c> June 2001 - added -dp- change in loop 39 to diapycnal flux diagnostics
c> Mar. 2002 - changed thickness smoothing from biharmonic to laplacian
c> Sep. 2003 - added diagnostics to track magnitude of -pbot- restoration
c> Sep. 2003 - added call to -recast- to improve tracer conservation
c> Dec. 2003 - confined interface smoothing to isopycnic coord. subdomain
c> Jan. 2004 - boosted thickness diffusion near sea floor (loops 141,151)
c> June 2004 - amended def'n of flxlo,flxhi in intfc.smoothing (loops 141,151)
c> Mar. 2006 - added bering strait exchange logic
c> Mar. 2006 - upgraded to biharmonic thickness diffusion
c> Apr. 2006 - iterate to find mass fluxes consistent with botm.pres restor.
c> June 2006 - changed handling of clipped parts of bolus fluxes
c> July 2006 - fixed bug in bolus flux computation (kn undefined in loop 7)
