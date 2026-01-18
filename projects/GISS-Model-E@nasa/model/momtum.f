#include "hycom_mpi_hacks.h"
      subroutine momtum(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 0.9 -- cyclic and noncyclic b.c. combined
      USE HYCOM_DIM
      USE HYCOM_SCALARS
      USE HYCOM_ARRAYS
      USE DOMAIN_DECOMP_1D, only : HALO_UPDATE,SOUTH,NORTH,
     &                          haveLatitude
      implicit none
c
      integer, intent(in) :: m,n,mm,nn
      integer             :: k1m,k1n
c
      integer ::  i,j,k,l,ia,ib,ja,jb,jp,km,kn
c
      real stress(idm,J_0H:J_1H),stresx(idm,J_0H:J_1H),
     .     stresy(idm,J_0H:J_1H),dpmxu(idm,J_0H:J_1H),
     .     dpmxv(idm,J_0H:J_1H),dpmx(idm,J_0H:J_1H),
     .     visc(idm,J_0H:J_1H),vort(idm,J_0H:J_1H),
     .     oneta(idm,J_0H:J_1H),wgtia(idm,J_0H:J_1H),
     .     wgtib(idm,J_0H:J_1H),wgtja(idm,J_0H:J_1H),
     .     wgtjb(idm,J_0H:J_1H),
     .     dpxy,dpia,dpib,dpja,dpjb,visca,viscb,ptopl,pbotl,cutoff,q,
     .     dt1inv,phi,plo,ubot,vbot,thkbop,thk,thka,thkb,avg,slab,
     .     olda,oldb,boost,pbot1,pbot2,p1,p2,botvel,drcoef

      real dl2u(idm,J_0H:J_1H),dl2uja(idm,J_0H:J_1H),
     .     dl2ujb(idm,J_0H:J_1H),dl2v(idm,J_0H:J_1H),
     .     dl2via(idm,J_0H:J_1H),dl2vib(idm,J_0H:J_1H)

      integer kan,jcyc
      character text*20
      real hfharm,sigstar
      external hfharm,sigstar
      data drcoef/.003/
c
      boost(pbot1,pbot2,p1,p2)=max(1.,1.5-min(pbot1-p1,pbot2-p2)
     .  /min(3.*tenm,.125*(pbot1+pbot2)))
c
c --- --------------------
c --- hydrostatic equation
c --- --------------------
c
      do 81 j=J_0,J_1
      do 81 l=1,isp(j)
c
      do 82 k=1,kk
      km=k+mm
      do 82 i=ifp(j,l),ilp(j,l)
c
c --- use upper interface pressure in converting sigma to sigma-star.
c --- this is to avoid density variations in layers intersected by sea floor
c
      if (kappa) then
        thstar(i,j,k)=sigstar(temp(i,j,k),saln(i,j,k),p(i,j,k))
      else
        thstar(i,j,k)=th3d(i,j,k)
      end if
 82   p(i,j,k+1)=p(i,j,k)+dp(i,j,km)
c
      do 80 i=ifp(j,l),ilp(j,l)
c
c --- store (1+eta) (= p_total/p_prime) in -oneta-
      oneta(i,j)=1.+pbavg(i,j,m)/p(i,j,kk+1)
c
c --- m_prime in lowest layer:
 80   montg(i,j,kk)=psikk(i,j)+(p(i,j,kk+1)*(thkk(i,j)-thstar(i,j,kk))
     .              -pbavg(i,j,m)*thstar(i,j,kk))*thref**2
c
c --- m_prime in remaining layers:
      do 81 k=kk-1,1,-1
      km=k+mm
      do 81 i=ifp(j,l),ilp(j,l)
 81   montg(i,j,k)=montg(i,j,k+1)+p(i,j,k+1)*oneta(i,j)
     .          *(thstar(i,j,k+1)-thstar(i,j,k))*thref**2
c
cdiag do j=jtest-1,jtest+1
cdiag do i=itest-1,itest+1
cdiag if (ip(i,j).gt.0) write (*,103) nstep,i,j,
cdiag. '    temp    saln  thstar   thkns    dpth   montg',
cdiag.  (k,temp(i,j,k+mm),saln(i,j,k+mm),thstar(i,j,k),
cdiag.   dp(i,j,k+mm)/onem,p(i,j,k+1)/onem,montg(i,j,k)/g,k=1,kk)
cdiag end do
cdiag end do
 103  format (i9,2i5,a/(i19,3f8.2,2f8.1,f8.3))
c
      call pardpudpv(mm)
c
c +++ ++++++++++++++++++
c +++ momentum equations
c +++ ++++++++++++++++++
c
c --- bottom drag (standard bulk formula)
c
      CALL HALO_UPDATE(ogrid,v, FROM=NORTH)
      CALL HALO_UPDATE(ogrid,vbavg, FROM=NORTH)

      thkbop=thkbot*onem
      do 804 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 804 l=1,isp(j)
c
      do 800 i=ifp(j,l),ilp(j,l)
      util1(i,j)=0.
 800  util2(i,j)=0.
c
      do 801 k=1,kk
      kn=k+nn
      do 801 i=ifp(j,l),ilp(j,l)
      phi=max(p(i,j,k+1),p(i,j,kk+1)-thkbop)
      plo=max(p(i,j,k  ),p(i,j,kk+1)-thkbop)
      util1(i,j)=util1(i,j)+(u(i,j,kn)+u(i+1,j,kn))*(phi-plo)
 801  util2(i,j)=util2(i,j)+(v(i,j,kn)+v(i,jb ,kn))*(phi-plo)
c
      do 804 i=ifp(j,l),ilp(j,l)
      ubot=ubavg(i,j,n)+ubavg(i+1,j,n)+util1(i,j)/thkbop
      vbot=vbavg(i,j,n)+vbavg(i,jb ,n)+util2(i,j)/thkbop
      botvel=.25*sqrt(ubot*ubot+vbot*vbot)+cbar
      ustarb(i,j)=sqrt(drcoef)*botvel
 804  drag(i,j)=min(drcoef*botvel/thkbot,.5/delt1)		! units: 1/s
c
c --- store r.h.s. of barotropic u/v eqn. in -ubrhs,vbrhs-
c --- time-interpolate wind stress

      CALL HALO_UPDATE(ogrid,depthv, FROM=NORTH)
      CALL HALO_UPDATE(ogrid,pvtrop, FROM=NORTH)
      CALL HALO_UPDATE(ogrid,ubavg,  FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,tauy  , FROM=SOUTH)

      CALL HALO_UPDATE(ogrid,montg,  FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,thstar, FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,p,      FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,dp  , FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,dpu,  FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,dpv,  FROM=NORTH)

      CALL HALO_UPDATE(ogrid,depthu, FROM=NORTH+SOUTH)

c
      do 70 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
c
      do 69 l=1,isu(j)
      do 69 i=ifu(j,l),ilu(j,l)
      ubrhs(i,j)=
     . (vbavg(i  ,j,m)*depthv(i  ,j)+vbavg(i  ,jb ,m)*depthv(i  ,jb )
     . +vbavg(i-1,j,m)*depthv(i-1,j)+vbavg(i-1,jb ,m)*depthv(i-1,jb ))
     . *(pvtrop(i,j)+pvtrop(i,jb ))*.125
c
      stresx(i,j)=(taux(i,j)+taux(i-1,j))*.5
     .            *g/min(ekman*onem,depthu(i,j))        !  units: m/s^2
c --- reduce stress under ice
c     stresx(i,j)=stresx(i,j)*(1.-.45*(oice(i,j)+oice(i-1,j)))
 69   continue
c
      do 70 l=1,isv(j)
      do 70 i=ifv(j,l),ilv(j,l)
      vbrhs(i,j)=
     .-(ubavg(i,j  ,m)*depthu(i,j  )+ubavg(i+1,j  ,m)*depthu(i+1,j  )
     . +ubavg(i,ja ,m)*depthu(i,ja )+ubavg(i+1,ja ,m)*depthu(i+1,ja ))
     . *(pvtrop(i,j)+pvtrop(i+1,j))*.125
c
      stresy(i,j)=(tauy(i,j)+tauy(i,ja))*.5
     .            *g/min(ekman*onem,depthv(i,j))        !  units: m/s^2
c --- reduce stress under ice
c     stresy(i,j)=stresy(i,j)*(1.-.45*(oice(i,j)+oice(i,ja )))
 70   continue
c
c --- the old  momeq2.f  starts here
c
      cutoff=5.*onem
c
      do 814 j=J_0,J_1
      do 814 i=1,ii
      dpmxu(i,j)=0.
      dpmxv(i,j)=0.
      dl2u(i,j)=0.
      dl2v(i,j)=0.
      vort(i,j)=huge			!  diagnostic use
c --- spatial weighting function for pressure gradient calculation:
      util1(i,j)=0.
 814  util2(i,j)=0.
c

      do 9 k=1,kk
      km=k+mm
      kn=k+nn
c
c --- store total (barotropic plus baroclinic) flow at old and mid time in
c --- -utotn,vtotn- and -utotm,vtotm- respectively. store minimum thickness
c --- values for use in pot.vort. calculation in -dpmx-.
c
      do 807 j=J_0,J_1
      do 807 l=1,isu(j)
      do 807 i=ifu(j,l),ilu(j,l)
      dpmxu(i,j)=dp(i,j,km)+dp(i-1,j,km)
      utotm(i,j)=u(i,j,km)+ubavg(i,j,m)
      utotn(i,j)=u(i,j,kn)+ubavg(i,j,n)
      uflux(i,j)=utotm(i,j)*max(dpu(i,j,km),cutoff)
 807  pu(i,j,k+1)=pu(i,j,k)+dpu(i,j,km)

c
      do 808 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      do 808 l=1,isv(j)
      do 808 i=ifv(j,l),ilv(j,l)
      dpmxv(i,j)=dp(i,j,km)+dp(i,ja ,km)
      vtotm(i,j)=v(i,j,km)+vbavg(i,j,m)
      vtotn(i,j)=v(i,j,kn)+vbavg(i,j,n)
      vflux(i,j)=vtotm(i,j)*max(dpv(i,j,km),cutoff)
 808  pv(i,j,k+1)=pv(i,j,k)+dpv(i,j,km)

      CALL HALO_UPDATE(ogrid,dpmxu  , FROM=SOUTH)
c
      do 803 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      do 803 i=1,ii
      ia=max(1,i-1)
 803  dpmx(i,j)=max(2.*cutoff,dpmxu(i,j),dpmxu(i  ,ja ),
     .                        dpmxv(i,j),dpmxv(ia ,j  ))
c
c --- define auxiliary velocity fields (via,vib,uja,ujb) to implement
c --- sidewall friction along near-vertical bottom slopes. wgtja,wgtjb,wgtia,
c --- wgtib indicate the extent to which a sidewall is present.

      CALL HALO_UPDATE(ogrid,utotn,  FROM=NORTH+SOUTH)
      CALL HALO_UPDATE(ogrid,vtotn,  FROM=NORTH+SOUTH)
c
      do 806 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
c
      do 805 l=1,isu(j)
      do 805 i=ifu(j,l),ilu(j,l)
      wgtja(i,j)=max(0.,min(1.,(pu(i,j,k+1)-depthu(i,ja))
     .          /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
      wgtjb(i,j)=max(0.,min(1.,(pu(i,j,k+1)-depthu(i,jb))
     .          /max(pu(i,j,k+1)-pu(i,j,k),epsil)))
      uja(i,j)=(1.-wgtja(i,j))*utotn(i,ja)+wgtja(i,j)*slip*utotn(i,j)
      ujb(i,j)=(1.-wgtjb(i,j))*utotn(i,jb)+wgtjb(i,j)*slip*utotn(i,j)
 805  dl2u(i,j)=utotn(i,j)
     .    -.25*(utotn(i+1,j)+utotn(i-1,j)+uja(i,j)+ujb(i,j))
c --- (to switch from biharmonic to laplacian friction, delete previous line)
c
      do 806 l=1,isv(j)
      do 806 i=ifv(j,l),ilv(j,l)
c
c --- if i=1, i-1 must point to zero-filled row (same for i+1 in case i=ii1)
      ia=mod(i-2+ii,ii)+1
      ib=i+1
c
      wgtia(i,j)=max(0.,min(1.,(pv(i,j,k+1)-depthv(ia,j))
     .          /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
      wgtib(i,j)=max(0.,min(1.,(pv(i,j,k+1)-depthv(ib,j))
     .          /max(pv(i,j,k+1)-pv(i,j,k),epsil)))
      via(i,j)=(1.-wgtia(i,j))*vtotn(ia,j)+wgtia(i,j)*slip*vtotn(i,j)
      vib(i,j)=(1.-wgtib(i,j))*vtotn(ib,j)+wgtib(i,j)*slip*vtotn(i,j)
 806  dl2v(i,j)=vtotn(i,j)
     .    -.25*(vtotn(i,jb )+vtotn(i,ja )+via(i,j)+vib(i,j))
c --- (to switch from biharmonic to laplacian friction, delete previous line)

c
c --- vorticity, pot.vort., defor. at lateral boundary points
      do 885 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
c
      do 885 l=1,isv(j)
      i=ifv(j,l)
      vort(i  ,j)= vtotm(i,j)*(1.-slip)*scvy(i,j)*scq2i(i  ,j)
      potvor(i  ,j)=(vort(i  ,j)+corio(i  ,j)) * 8.
     ./max(8.*cutoff,4.*(dp(i,j,km)+dp(i,ja ,km)),dpmx(i,j),dpmx(i+1,j))
      defor2(i  ,j)=(vtotn(i,j)*(1.-slip)*scvy(i,j))**2*scq2i(i  ,j)
      i=ilv(j,l)
      vort(i+1,j)=-vtotm(i,j)*(1.-slip)*scvy(i,j)*scq2i(i+1,j)
      potvor(i+1,j)=(vort(i+1,j)+corio(i+1,j)) * 8.
     ./max(8.*cutoff,4.*(dp(i,j,km)+dp(i,ja ,km)),dpmx(i,j),dpmx(i+1,j))
 885  defor2(i+1,j)=(vtotn(i,j)*(1.-slip)*scvy(i,j))**2*scq2i(i+1,j)

      CALL HALO_UPDATE(ogrid,dpmx,  FROM=NORTH)
      CALL HALO_UPDATE(ogrid,utotm, FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,scux,  FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,dpmx,  FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,utotn, FROM=SOUTH)
c
      do 886 i=1,ii1
c
c --- if i=1, i-1 must point to zero-filled row (same for i+1 in case i=ii1)
      ia=mod(i-2+ii,ii)+1
c
      do 886 l=1,jsu(i)
      j=jfu(i,l)
      if (haveLatitude(ogrid, J=j)) then
        jb = PERIODIC_INDEX(j+1, jj)
        !jb=mod(j     ,jj)+1
        vort(i,j  )=-utotm(i,j)*(1.-slip)*scux(i,j)*scq2i(i,j  )
        potvor(i,j  )=(vort(i,j  )+corio(i,j  )) * 8.
     .  /max(8.*cutoff,4.*(dp(i,j,km)+dp(ia ,j,km)),dpmx(i,j),
     .   dpmx(i,jb ))
        defor2(i,j  )=(utotn(i,j)*(1.-slip)*scux(i,j))**2*scq2i(i,j  )
      endif
      j=jlu(i,l)
      jb=mod(j     ,jj)+1
      if (haveLatitude(ogrid, J=jb)) then
        j = PERIODIC_INDEX(jb-1, jj)
        vort(i,jb )= utotm(i,j)*(1.-slip)*scux(i,j)*scq2i(i,jb )
        potvor(i,jb )=(vort(i,jb )+corio(i,jb )) * 8.
     .  /max(8.*cutoff,4.*(dp(i,j,km)+dp(ia ,j,km)),dpmx(i,j),
     .   dpmx(i,jb ))
        defor2(i,jb )=(utotn(i,j)*(1.-slip)*scux(i,j))**2*scq2i(i,jb )
      endif
 886  continue

      CALL HALO_UPDATE(ogrid,vtotn, FROM=NORTH)
      CALL HALO_UPDATE(ogrid,scvx,  FROM=NORTH)
c
c --- vorticity, pot.vort., defor. at interior points (incl. promontories).
c --- defor1 = du/dx-dv/dy at mass points, defor2 = dv/dx+du/dy at vort. points
c
      do 63 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 63 l=1,isp(j)
      do 63 i=ifp(j,l),ilp(j,l)
      util3(i,j)=.5*(p(i,j,k+1)+p(i,j,k))*oneta(i,j)
 63   defor1(i,j)=((utotn(i+1,j)*scuy(i+1,j)-utotn(i,j)*scuy(i,j))
     .            -(vtotn(i,jb )*scvx(i,jb )-vtotn(i,j)*scvx(i,j)))**2
     .            *scp2i(i,j)

      CALL HALO_UPDATE(ogrid,utotm, FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,scux,  FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,dpmx,  FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,ujb,   FROM=SOUTH)
c
c --- vorticity, pot.vort., defor. at interior points (incl. promontories)
      do 64 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
      do 64 l=1,isq(j)
      do 64 i=ifq(j,l),ilq(j,l)
      vort(i,j)=(vtotm(i,j)*scvy(i,j)-vtotm(i-1,j)*scvy(i-1,j)
     .          -utotm(i,j)*scux(i,j)+utotm(i,ja )*scux(i,ja ))
     .          *scq2i(i,j)
      potvor(i,j)=(vort(i,j)+corio(i,j)) * 8.
     .   /max(8.*cutoff,2.*(dp(i,j  ,km)+dp(i-1,j  ,km)+
     .                      dp(i,ja ,km)+dp(i-1,ja ,km))
     .   ,dpmx(i,j),dpmx(i-1,j),dpmx(i+1,j),dpmx(i,ja ),dpmx(i,jb ))
 64   defor2(i,j)=(vib(i-1,j)*scvy(i,j)-via(i,j)*scvy(i-1,j)
     .            +ujb(i,ja )*scux(i,j)-uja(i,j)*scux(i,ja ))**2
     .            *scq2i(i,j)

      CALL HALO_UPDATE(ogrid,dl2u,  FROM=SOUTH+NORTH)
c
c --- define auxiliary del2 fields (dl2via,dl2vib,dl2uja,dl2ujb) to imple-
c --- ment biharmonic sidewall friction along near-vertical bottom slopes.
c
      do 906 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
c
      do 905 l=1,isu(j)
      do 905 i=ifu(j,l),ilu(j,l)
      dl2uja(i,j)=(1.-wgtja(i,j))*dl2u(i,ja)+wgtja(i,j)*slip*dl2u(i,j)
 905  dl2ujb(i,j)=(1.-wgtjb(i,j))*dl2u(i,jb)+wgtjb(i,j)*slip*dl2u(i,j)
c
      do 906 l=1,isv(j)
      do 906 i=ifv(j,l),ilv(j,l)
c
c --- if i=1, i-1 must point to zero-filled row (same for i+1 in case i=ii1)
      ia=mod(i-2+ii,ii)+1
      ib=i+1
c
      dl2via(i,j)=(1.-wgtia(i,j))*dl2v(ia,j)+wgtia(i,j)*slip*dl2v(i,j)
 906  dl2vib(i,j)=(1.-wgtib(i,j))*dl2v(ib,j)+wgtib(i,j)*slip*dl2v(i,j)
c
ccc      do j=1,jdm
ccc      ja=mod(j-2+jj,jj)+1
ccc      do i=1,idm
ccc      ia=mod(i-2+ii,ii)+1
ccc      mask(i,j)=iq(i,j)+iu(i,j)+iv(i,j)
ccc      mask(i,j)=mask(i,j)+iu(i,ja )
ccc      mask(i,j)=mask(i,j)+iv(ia ,j)
ccc      end do
ccc      end do
ccc      write (text,'(a,i3,i8)') 'potvor k=',k,nstep
ccc      call compare(potvor,mask,text)
c
c --- ----------
c --- u equation
c --- ----------
c
c --- deformation-dependent eddy viscosity coefficient

      CALL HALO_UPDATE(ogrid,defor2, FROM=NORTH)
c
      do 37 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 37 l=1,isu(j)
      do 37 i=ifu(j,l),ilu(j,l)
      visc(i,j)=max(veldff,viscos*
     .sqrt(.5*(defor1(i,j)+defor1(i-1,j)+defor2(i,j)+defor2(i,jb ))))
     .   *boost(pbot(i,j),pbot(i-1,j),p(i,j,k),p(i-1,j,k))
 37   continue

      CALL HALO_UPDATE(ogrid,visc, FROM=SOUTH+NORTH)
ccc      CALL HALO_UPDATE(ogrid,glue, FROM=SOUTH+NORTH)	! not used?
      CALL HALO_UPDATE(ogrid,scqx, FROM=      NORTH)

c
      do 822 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
c
      do 820 l=1,isu(j)
      visc(ifu(j,l)-1,j)=visc(ifu(j,l),j)
 820  visc(ilu(j,l)+1,j)=visc(ilu(j,l),j)
c
c --- longitudinal turb. momentum flux (at mass points)
c
      do 824 l=1,isp(j)
      do 824 i=ifp(j,l),ilp(j,l)
 824  if (iu(i,j)+iu(i+1,j).gt.0)
     . uflux1(i,j)=(visc(i,j)+visc(i+1,j))*(dl2u(i,j)-dl2u(i+1,j))
     .             *hfharm(max(dpu(i  ,j,km),onemm),
     .                     max(dpu(i+1,j,km),onemm))*scpy(i,j)
ccc     .             *glue(i,j)
c
c --- lateral turb. momentum flux (at vorticity points)
c --- (left and right fluxes are evaluated separately because of sidewalls)
c
      do 822 l=1,isu(j)
      do 822 i=ifu(j,l),ilu(j,l)
      dpxy=max(dpu(i,j ,km),onemm)
      dpja=max(dpu(i,ja,km),onemm)
      dpjb=max(dpu(i,jb,km),onemm)
c
c --- check whether variables along coast have been initialized correctly
cdiag if (k.eq.kk) then
cdiag   if (iu(i,ja).eq.0 .and. dpu(i,ja,km).ne.0.) write
cdiag.   (*,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.   '  error - nonzero dpu(ja):',dpu(i,ja,km)
cdiag   if (iu(i,jb).eq.0 .and. dpu(i,jb,km).ne.0.) write
cdiag.   (*,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.   '  error - nonzero dpu(jb):',dpu(i,jb,km)
cdiag end if
c
      if (iu(i,ja).eq.0) then
        visca=visc(i,j )
      else
        visca=visc(i,ja)
      end if
      if (iu(i,jb).eq.0) then
        viscb=visc(i,j )
      else
        viscb=visc(i,jb)
      end if
      uflux2(i,j)=(visc(i,j)+visca)*(dl2uja(i,j)-dl2u(i,j))
     .            *hfharm(dpja+wgtja(i,j)*(dpxy-dpja),dpxy)*scqx(i,j )
ccc     .   *.25*(glue(i,j)+glue(i-1,j)+glue(i,ja )+glue(i-1,ja ))
 822  uflux3(i,j)=(visc(i,j)+viscb)*(dl2u(i,j)-dl2ujb(i,j))
     .            *hfharm(dpjb+wgtjb(i,j)*(dpxy-dpjb),dpxy)*scqx(i,jb)
ccc     .   *.25*(glue(i,j)+glue(i-1,j)+glue(i,jb )+glue(i-1,jb ))
c
c --- pressure force in x direction
c --- ('scheme 2' from appendix -a- in bleck-smith paper)
c
      do 96 j=J_0,J_1
      do 96 l=1,isu(j)
      do 96 i=ifu(j,l),ilu(j,l)
      util1(i,j)=max(0.,min(depthu(i,j)-pu(i,j,k),h1))
 96   pgfx(i,j)=util1(i,j)*
     .    (montg(i,j,k)-montg(i-1,j,k)+(thstar(i,j,k)-thstar(i-1,j,k))
     .      *(p(i,j,k+1)*p(i-1,j,k+1)-p(i,j,k)*p(i-1,j,k))*thref**2
     .      /(dp(i,j,km)+dp(i-1,j,km)+epsil))

      CALL HALO_UPDATE(ogrid,pgfx,  FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,util1, FROM=SOUTH+NORTH)
c
      do 98 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
c
      do 98 l=1,isu(j)
      do 98 i=ifu(j,l),ilu(j,l)
c
c --- check whether variables along coast have been initialized correctly
cdiag if (k.eq.kk) then
cdiag   if (iu(i,ja).eq.0 .and. pgfx(i,ja).ne.0.) write
cdiag.   (*,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.   '  error - nonzero pgfx(ja):',pgfx(i,ja)
cdiag   if (iu(i,jb).eq.0 .and. pgfx(i,jb).ne.0.) write
cdiag.   (*,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.   '  error - nonzero pgfx(jb):',pgfx(i,jb)
cdiag end if
c
 98   gradx(i,j)=(pgfx(i,j)+(h1-util1(i,j))*
     .  (pgfx (i-1,j)+pgfx (i+1,j)+pgfx (i,ja)+pgfx (i,jb))/
     .  (util1(i-1,j)+util1(i+1,j)+util1(i,ja)+util1(i,jb)+epsil))/h1

      CALL HALO_UPDATE(ogrid,vtotm,  FROM=NORTH)
      CALL HALO_UPDATE(ogrid,vflux,  FROM=NORTH)
      CALL HALO_UPDATE(ogrid,potvor, FROM=NORTH)
c
      do 6 j=J_0,J_1
      jb = PERIODIC_INDEX(j+1, jj)
      do 6 l=1,isu(j)
      do 6 i=ifu(j,l),ilu(j,l)
c
      ptopl=min(depthu(i,j),.5*(p(i,j,k  )+p(i-1,j,k  )))
      pbotl=min(depthu(i,j),.5*(p(i,j,k+1)+p(i-1,j,k+1)))
c
c --- top and bottom boundary layer stress. stress profile is assumed linear
      stress(i,j)=(-utotn(i,j)*.5*(drag(i,j)+drag(i-1,j))*
     .     (max(depthu(i,j)-thkbop,          pbotl       )
     .     -max(depthu(i,j)-thkbop,min(ptopl,pbotl-onemm)))
     .     +stresx(i,j)*(min(ekman*onem,pbotl+onemm)
     .                  -min(ekman*onem,ptopl      )))
     .     /max(dpu(i,j,km),onemm)
c
c --- time smoothing of -u- field  (part 1)
      u(i,j,km)=u(i,j,km)*(wuv1*dpu(i,j,km)+onemm)
     .         +u(i,j,kn)* wuv2*dpu(i,j,kn)
c
      util4(i,j)=u(i,j,kn)
 6    u(i,j,kn)=u(i,j,kn)+delt1*(-scuxi(i,j)*(gradx(i,j)
     .+.25*(utotm(i+1,j)**2+vtotm(i  ,j)**2+vtotm(i  ,jb )**2
     .     -utotm(i-1,j)**2-vtotm(i-1,j)**2-vtotm(i-1,jb )**2))
     .+.125*(vflux(i  ,j)+vflux(i  ,jb )+vflux(i-1,j)+vflux(i-1,jb ))
     .     *(potvor(i,j)+potvor(i,jb )) - ubrhs(i,j) + stress(i,j)
     .-(uflux1(i,j)-uflux1(i-1,j)
     . +uflux3(i,j)-uflux2(i,j))/(scu2(i,j)*max(dpu(i,j,km),onemm)))
c
c --- set baroclinic velocity to zero one point away from bering strait seam
      if (haveLatitude(ogrid, J=jpac)) u(ipacs,jpac,kn)=0.
      if (haveLatitude(ogrid, J=jatl)) u(iatln,jatl,kn)=0.
c
      if (itest.gt.0 .and. jtest.gt.0) then
      if (jtest.ge.J_0 .and. jtest.le.J_1) then
      write (*,100) nstep
      do j=max(J_0, jtest-1), min(J_1, jtest+1)
      jb = PERIODIC_INDEX(j+1, jj)
      do i=itest-1,itest+1
      if (iu(i,j).gt.0) then
      write (*,'(2i5,i3,2p,8f8.3)') i,j,k,
     .  util4(i,j),u(i,j,kn),-delt1*gradx(i,j)*scuxi(i,j),
     .  -delt1*scuxi(i,j)*
     . .25*(utotm(i+1,j)**2+vtotm(i  ,j)**2+vtotm(i  ,jb )**2
     .     -utotm(i-1,j)**2-vtotm(i-1,j)**2-vtotm(i-1,jb )**2),
     .   delt1*(vflux(i  ,j)+vflux(i  ,jb )+vflux(i-1,j)+vflux(i-1,jb ))
     .        *(potvor(i,j)+potvor(i,jb ))*.125,
     .  -delt1*ubrhs(i,j),delt1*stress(i,j),
     .  -delt1*(uflux1(i,j)-uflux1(i-1,j)
     .         +uflux3(i,j)-uflux2(i,j))/
     .         (scu2(i,j)*max(dpu(i,j,km),onemm))
      end if
      end do
      end do
 100  format(i9,8x,'uold    unew   gradp  nonlin   corio',
     .3x,'ubrhs  stress    fric')
      end if
      end if
c
c --- ----------
c --- v equation
c --- ----------
c
c --- deformation-dependent eddy viscosity coefficient

      CALL HALO_UPDATE(ogrid,defor1  , FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,pbot  ,   FROM=SOUTH)
c
      do 38 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      do 38 l=1,isv(j)
      do 38 i=ifv(j,l),ilv(j,l)
      visc(i,j)=max(veldff,viscos*
     .sqrt(.5*(defor1(i,j)+defor1(i,ja )+defor2(i,j)+defor2(i+1,j))))
     .   *boost(pbot(i,j),pbot(i,ja ),p(i,j,k),p(i,ja ,k))
 38   continue

      CALL HALO_UPDATE(ogrid,visc, FROM=SOUTH+NORTH)
c
      do 821 i=1,ii1
      do 821 l=1,jsv(i)
      j=jfv(i,l)
      ja=mod(j-2+jj,jj)+1
      if (haveLatitude(ogrid, J=ja)) then
        if (j.ne.1  .or. jlv(i,jsv(i)).ne.jj)
     &     visc(i,ja)=visc( i,PERIODIC_INDEX(ja+1, jj) )
      endif
      j=jlv(i,l)
      jb=mod(j     ,jj)+1
      if (haveLatitude(ogrid, J=jb)) then
        if (j.ne.jj .or. jfv(i,     1).ne.1 )
     &     visc(i,jb)=visc( i,PERIODIC_INDEX(jb-1, jj) )
      endif
 821  continue

      CALL HALO_UPDATE(ogrid,visc, FROM=NORTH)
      CALL HALO_UPDATE(ogrid,dl2v, FROM=NORTH)
ccc   CALL HALO_UPDATE(ogrid,glue, FROM=SOUTH)		! not used?

c
c --- longitudinal turb. momentum flux (at mass points)
c
      do 823 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
c
      do 825 l=1,isp(j)
      do 825 i=ifp(j,l),ilp(j,l)
 825  if (iv(i,j)+iv(i,jb ).gt.0)
     . vflux1(i,j)=(visc(i,j)+visc(i,jb ))*(dl2v(i,j)-dl2v(i,jb ))
     .             *hfharm(max(dpv(i,j  ,km),onemm),
     .                     max(dpv(i,jb ,km),onemm))*scpx(i,j)
ccc     .             *glue(i,j)
c
c --- lateral turb. momentum flux (at vorticity points)
c --- (left and right fluxes are evaluated separately because of sidewalls)
c
      do 823 l=1,isv(j)
      do 823 i=ifv(j,l),ilv(j,l)
c
c --- if i=1, i-1 must point to zero-filled row (same for i+1 in case i=ii1)
      ia=mod(i-2+ii,ii)+1
      ib=i+1
c
      dpxy=max(dpv(i ,j,km),onemm)
      dpia=max(dpv(ia,j,km),onemm)
      dpib=max(dpv(ib,j,km),onemm)
c
c --- check whether variables along coast have been initialized correctly
cdiag if (k.eq.kk) then
cdiag   if (iv(ia,j).eq.0 .and. dpv(ia,j,km).ne.0.) write
cdiag.   (*,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.   '  error - nonzero dpv(ia):',dpv(ia,j,km)
cdiag   if (iv(ib,j).eq.0 .and. dpv(ib,j,km).ne.0.) write
cdiag.   (*,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.   '  error - nonzero dpv(ib):',dpv(ib,j,km)
cdiag end if
c
      if (iv(ia,j).eq.0) then
        visca=visc(i ,j)
      else
        visca=visc(ia,j)
      end if
      if (iv(ib,j).eq.0) then
        viscb=visc(i ,j)
      else
        viscb=visc(ib,j)
      end if
      vflux2(i,j)=(visc(i,j)+visca)*(dl2via(i,j)-dl2v(i,j))
     .            *hfharm(dpia+wgtia(i,j)*(dpxy-dpia),dpxy)*scqy(i ,j)
ccc     .   *.25*(glue(i,j)+glue(i-1,j)+glue(i,ja )+glue(i-1,ja ))
 823  vflux3(i,j)=(visc(i,j)+viscb)*(dl2v(i,j)-dl2vib(i,j))
     .            *hfharm(dpib+wgtib(i,j)*(dpxy-dpib),dpxy)*scqy(ib,j)
ccc     .   *.25*(glue(i,j)+glue(i+1,j)+glue(i,ja )+glue(i+1,ja ))
c
c --- pressure force in y direction
c --- ('scheme 2' from appendix -a- in bleck-smith paper)

c
      do 97 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      do 97 l=1,isv(j)
      do 97 i=ifv(j,l),ilv(j,l)
      util2(i,j)=max(0.,min(depthv(i,j)-pv(i,j,k),h1))
 97   pgfy(i,j)=util2(i,j)*
     .    (montg(i,j,k)-montg(i,ja ,k)+(thstar(i,j,k)-thstar(i,ja ,k))
     .      *(p(i,j,k+1)*p(i,ja ,k+1)-p(i,j,k)*p(i,ja ,k))*thref**2
     .      /(dp(i,j,km)+dp(i,ja ,km)+epsil))

      CALL HALO_UPDATE(ogrid,pgfy,  FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,util2, FROM=SOUTH+NORTH)
c
      do 99 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
      do 99 l=1,isv(j)
      do 99 i=ifv(j,l),ilv(j,l)
c
c --- if i=1, i-1 must point to zero-filled row (same for i+1 in case i=ii1)
      ia=mod(i-2+ii,ii)+1
      ib=i+1
c
c --- check whether variables along coast have been initialized correctly
cdiag if (k.eq.kk) then
cdiag   if (iv(ia,j).eq.0 .and. pgfy(ia,j).ne.0.) write
cdiag.   (*,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.   '  error - nonzero pgfy(ia):',pgfy(ia,j)
cdiag   if (iv(ib,j).eq.0 .and. pgfy(ib,j).ne.0.) write
cdiag.   (*,'(i9,2i5,a,1p,2e9.1)') nstep,i,j,
cdiag.   '  error - nonzero pgfy(ib):',pgfy(ib,j)
cdiag end if
c
 99   grady(i,j)=(pgfy(i,j)+(h1-util2(i,j))*
     .  (pgfy (ia ,j)+pgfy (ib ,j)+pgfy (i,ja )+pgfy (i,jb ))/
     .  (util2(ia ,j)+util2(ib ,j)+util2(i,ja )+util2(i,jb )+epsil))/h1

      CALL HALO_UPDATE(ogrid,drag,   FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,utotm,  FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,uflux,  FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,vflux1, FROM=SOUTH)
      CALL HALO_UPDATE(ogrid,vtotm,  FROM=SOUTH+NORTH)
c
      do 7 j=J_0,J_1
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
      do 7 l=1,isv(j)
      do 7 i=ifv(j,l),ilv(j,l)
c
      ptopl=min(depthv(i,j),.5*(p(i,j,k  )+p(i,ja ,k  )))
      pbotl=min(depthv(i,j),.5*(p(i,j,k+1)+p(i,ja ,k+1)))
c
c --- top and bottom boundary layer stress. stress profile is assumed linear
      stress(i,j)=(-vtotn(i,j)*.5*(drag(i,j)+drag(i,ja ))*
     .     (max(depthv(i,j)-thkbop,          pbotl       )
     .     -max(depthv(i,j)-thkbop,min(ptopl,pbotl-onemm)))
     .     +stresy(i,j)*(min(ekman*onem,pbotl+onemm)
     .                  -min(ekman*onem,ptopl      )))
     .     /max(dpv(i,j,km),onemm)
c
c --- time smoothing of -v- field  (part 1)
      v(i,j,km)=v(i,j,km)*(wuv1*dpv(i,j,km)+onemm)
     .         +v(i,j,kn)* wuv2*dpv(i,j,kn)
c
      util4(i,j)=v(i,j,kn)
 7    v(i,j,kn)=v(i,j,kn)+delt1*(-scvyi(i,j)*(grady(i,j)
     .+.25*(vtotm(i,jb )**2+utotm(i,j  )**2+utotm(i+1,j  )**2
     .     -vtotm(i,ja )**2-utotm(i,ja )**2-utotm(i+1,ja )**2))
     .-.125*(uflux(i,j  )+uflux(i+1,j  )+uflux(i,ja )+uflux(i+1,ja ))
     .     *(potvor(i,j)+potvor(i+1,j)) - vbrhs(i,j) + stress(i,j)
     .-(vflux1(i,j)-vflux1(i,ja )
     . +vflux3(i,j)-vflux2(i,j))/(scv2(i,j)*max(dpv(i,j,km),onemm)))
c
      if (itest.gt.0 .and. jtest.gt.0) then
      if (jtest.ge.J_0 .and. jtest.le.J_1) then
      write (*,101) nstep
      do j=max(J_0, jtest-1), min(J_1, jtest+1)
      ja = PERIODIC_INDEX(j-1, jj)
      jb = PERIODIC_INDEX(j+1, jj)
      do i=itest-1,itest+1
      if (iv(i,j).gt.0) then
      write (*,'(2i5,i3,2p,8f8.3)') i,j,k,
     .  util4(i,j),v(i,j,kn),-delt1*grady(i,j)*scvyi(i,j),
     .  -delt1*scvyi(i,j)*
     . .25*(vtotm(i,jb )**2+utotm(i,j  )**2+utotm(i+1,j  )**2
     .     -vtotm(i,ja )**2-utotm(i,ja )**2-utotm(i+1,ja )**2),
     .  -delt1*(uflux(i,j  )+uflux(i+1,j  )+uflux(i,ja )+uflux(i+1,ja ))
     .        *(potvor(i,j)+potvor(i+1,j))*.125,
     .  -delt1*vbrhs(i,j),delt1*stress(i,j),
     .  -delt1*(vflux1(i,j)-vflux1(i,ja )
     .         +vflux3(i,j)-vflux2(i,j))/
     .         (scv2(i,j)*max(dpv(i,j,km),onemm))
      end if
      end do
      end do
 101  format(i9,8x,'vold    vnew   gradp  nonlin   corio',
     .3x,'vbrhs  stress    fric')
      end if
      end if
c
 9    continue
c
      dt1inv = 1./delt1
c
      do j=J_0,J_1
        do 14 k=1,kk
        kn=k+nn
c
        do 12 l=1,isp(j)
        do 12 i=ifp(j,l),ilp(j,l)
 12     p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
c --- compute new -dpu,dpv- field. save old -dpu,dpv- values in -pu,pv-.
c
        do 13 l=1,isu(j)
        do 13 i=ifu(j,l),ilu(j,l)
 13     pu(i,j,k+1)=dpu(i,j,kn)
c
        do 14 l=1,isv(j)
        do 14 i=ifv(j,l),ilv(j,l)
 14     pv(i,j,k+1)=dpv(i,j,kn)
      end do
c
      call pardpudpv(nn)
c
c --- extract barotropic velocities generated during most recent baroclinic
c --- time step and use them to force barotropic flow field.
c
      slab=onem*vertmx*delt1
      do j=J_0,J_1
c
      do 31 l=1,isu(j)
c
      do 32 i=ifu(j,l),ilu(j,l)
 32   utotn(i,j)=0.
      do 33 k=1,kk
      km=k+mm
      kn=k+nn
      kan=max(1,k-1)+nn
      do 33 i=ifu(j,l),ilu(j,l)
      if (vertmx.eq.0.) then
c --- traditional vertical momentum mixing
        q=min(dpu(i,j,km),dpu(i,j,kn),onem)
        u(i,j,kn)=(u(i,j,kn)*q+u(i,j,kan)*(onem-q))/onem
      else
c --- homogenize velocity over 'slab' straddling each interface
        if (k.lt.kk) then
cdiag     olda=u(i,j,kn  )
cdiag     oldb=u(i,j,kn+1)
          thka=min(dpu(i,j,kn  ),.5*slab)
          thkb=min(dpu(i,j,kn+1),.5*slab)
          thka=          slab-thka
          thkb=max(epsil,slab-thka)
          avg=(u(i,j,kn)*thka+u(i,j,kn+1)*thkb)/(thka+thkb)
          thk=max(0.,dpu(i,j,kn  )-thka)
          u(i,j,kn  )=(u(i,j,kn  )*thk+avg*thka)/(thk+thka)
          thk=max(0.,dpu(i,j,kn+1)-thkb)
          u(i,j,kn+1)=(u(i,j,kn+1)*thk+avg*thkb)/(thk+thkb)
c
cdiag     if (i.eq.itest .and. j.eq.jtest)
cdiag.     write (*,104) nstep,i,j,k,'ua,ub:',dpu(i,j,kn)/onem,
cdiag.      dpu(i,j,kn+1)/onem,olda,oldb,u(i,j,kn),u(i,j,kn+1)
 104       format (i7,2i5,i3,' dpa,dpb,',a,3(f8.2,f7.2))
        end if
      end if
 33   utotn(i,j)=utotn(i,j)+u(i,j,kn)*dpu(i,j,kn)
      do 31 i=ifu(j,l),ilu(j,l)
 31   utotn(i,j)=utotn(i,j)/depthu(i,j)
c
      do 30 l=1,isv(j)
c
      do 34 i=ifv(j,l),ilv(j,l)
 34   vtotn(i,j)=0.
      do 35 k=1,kk
      km=k+mm
      kn=k+nn
      kan=max(1,k-1)+nn
      do 35 i=ifv(j,l),ilv(j,l)
      if (vertmx.eq.0.) then
c --- traditional vertical momentum mixing
        q=min(dpv(i,j,km),dpv(i,j,kn),onem)
        v(i,j,kn)=(v(i,j,kn)*q+v(i,j,kan)*(onem-q))/onem
      else
c --- homogenize velocity over 'slab' straddling each interface
        if (k.lt.kk) then
cdiag     olda=v(i,j,kn  )
cdiag     oldb=v(i,j,kn+1)
          thka=min(dpv(i,j,kn  ),.5*slab)
          thkb=min(dpv(i,j,kn+1),.5*slab)
          thka=          slab-thka
          thkb=max(epsil,slab-thka)
          avg=(v(i,j,kn)*thka+v(i,j,kn+1)*thkb)/(thka+thkb)
          thk=max(0.,dpv(i,j,kn  )-thka)
          v(i,j,kn  )=(v(i,j,kn  )*thk+avg*thka)/(thk+thka)
          thk=max(0.,dpv(i,j,kn+1)-thkb)
          v(i,j,kn+1)=(v(i,j,kn+1)*thk+avg*thkb)/(thk+thkb)
c
cdiag   if (i.eq.itest .and. j.eq.jtest)
cdiag.   write (*,104) nstep,i,j,k,'va,vb:',dpv(i,j,kn)/onem,
cdiag.    dpv(i,j,kn+1)/onem,olda,oldb,v(i,j,kn),v(i,j,kn+1)
        end if
      end if
 35   vtotn(i,j)=vtotn(i,j)+v(i,j,kn)*dpv(i,j,kn)
      do 30 i=ifv(j,l),ilv(j,l)
 30   vtotn(i,j)=vtotn(i,j)/depthv(i,j)
      end do
c
c --- time smoothing of -u,v- fields  (part 2)
c
      do k=1,kk
      km=k+mm
      kn=k+nn
c
      do 22 j=J_0,J_1
c
      do 24 l=1,isu(j)
      do 24 i=ifu(j,l),ilu(j,l)
      u(i,j,kn)=u(i,j,kn)-utotn(i,j)
      u(i,j,km)=(u(i,j,km)+u(i,j,kn)*wuv2*dpu(i,j,kn))/
     .   (wuv1*dpu(i,j,km)+onemm+wuv2*(pu(i,j,k+1)+dpu(i,j,kn)))
c --- build up time integral of velocity field
      uav  (i,j,k)=uav  (i,j,k)+(u(i,j,km)+ubavg(i,j,m))*dpu(i,j,km)
      dpuav(i,j,k)=dpuav(i,j,k)+          dpu(i,j,km)
 24   continue
c
      do 22 l=1,isv(j)
      do 22 i=ifv(j,l),ilv(j,l)
      v(i,j,kn)=v(i,j,kn)-vtotn(i,j)
      v(i,j,km)=(v(i,j,km)+v(i,j,kn)*wuv2*dpv(i,j,kn))/
     .   (wuv1*dpv(i,j,km)+onemm+wuv2*(pv(i,j,k+1)+dpv(i,j,kn)))
c --- build up time integral of velocity field
      vav  (i,j,k)=vav  (i,j,k)+(v(i,j,km)+vbavg(i,j,m))*dpv(i,j,km)
      dpvav(i,j,k)=dpvav(i,j,k)+          dpv(i,j,km)
 22   continue
      end do
c
      do 867 j=J_0,J_1
c
      do 865 l=1,isu(j)
      do 865 i=ifu(j,l),ilu(j,l)
      utotn(i,j)=utotn(i,j)*dt1inv
      ubavg(i,j,n)=ubavg(i,j,m)
 865  continue
c
      do 866 l=1,isv(j)
      do 866 i=ifv(j,l),ilv(j,l)
      vtotn(i,j)=vtotn(i,j)*dt1inv
      vbavg(i,j,n)=vbavg(i,j,m)
 866  continue
c
      do 867 l=1,isp(j)
      do 867 i=ifp(j,l),ilp(j,l)
      pbavg(i,j,n)=pbavg(i,j,m)
      pbavav(i,j)=pbavav(i,j)+pbavg(i,j,m)
      sfhtav(i,j)=sfhtav(i,j)+(montg(i,j,1)+thref*pbavg(i,j,m))/g
 867  continue
c
      return
      end
c
c
      real function hfharm(a,b)
c --- harmonic mean divided by 2
      hfharm=a*b/(a+b)
      return
      end
c
c
c> Revision history
c>
c> Mar. 2000 - conversion to SI units
c> Apr. 2000 - changed 'thkbop' in drag formula (stmt. 804) to 'thkbot'
c> Apr. 2000 - fixed dimension error in stress formula caused by use of 'ekman'
c> June 2000 - modified j-1,j+1 to accomodate both channel & closed basin b.c.
c> Aug. 2000 - separated dpmx(i,j) and dpmx(i,jb) calculation in loop 807
c>             to avoid coincidental execution in multi-threaded runs
c> Sept.2000 - removed erroneous PRIVATE declaration for -thkbop-
c> Oct. 2000 - changed 'thkbot' in drag formula (stmt. 804) to 'thkbop/g'
c> Oct. 2000 - limited bottom drag to  1/(model time step)
c> Apr. 2001 - eliminated stmt_funcs.h
c> Dec. 2001 - added 'glue' to lateral mixing operation
c> Jan. 2004 - boosted viscosity near sea floor (loops 37,38)
c> Mar. 2006 - added bering strait exchange logic
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
