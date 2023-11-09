      subroutine eice(m,n,mm,nn,k1m,k1n)
c
c --- estimate the ice formed (odmsi) in order to keep ocean above freezing
c --- and corresponding saltflux (salflx2)
c
      USE SEAICE, only : fsss,tfrez,Ei
c
      USE HYCOM_DIM,only : isp,ifp,ilp,kk,idm,jchunk,
     *                     J_0,J_1,J_0H,J_1H
      USE HYCOM_SCALARS, only : thkmin,onem,nstep,delt1,g,spcifh
     &     ,equatn,epsil,brntop,brnbot,itest,jtest
      USE HYCOM_ARRAYS
c
      implicit none
c
      integer i,j,k,l,m,n,mm,nn,kn,k1m,k1n
c
      real tmelt,thin,rhoice,kice,fusion,saldif,rate,tmxl,dpth,thkmax,
     .     paybak,borrow,total,totalx,qmax,qmx(J_0H:J_1H),heatfx,dtdflx,
     .     icex(idm,J_0H:J_1H),work(idm,J_0H:J_1H),dm,
     .     salflx2(idm,J_0H:J_1H),top,bot,
     .     thkinv
      integer imx(J_0H:J_1H),jmx(J_0H:J_1H),k1
      logical dosmoo,pump,vrbos
c
c --- tmelt  = melting point (deg)
c --- thin   = min. ice thickness (m)
c --- thkmax = max. ice thickness (m)
c --- rhoice = ice density (kg/m^3)
c --- kice   = heat conductivity in ice (W/m/deg)
c --- fusion = latent heat of fusion (J/kg)
c --- saldif = salinity difference water minus ice
c --- rate   = max. ice freezing and melting rate (m/sec)
c --- dtdflx = d (srf.temperature) / d (heat flux)  (deg m^2/W)
c --- sfxice = net total heat flux between atm and ice (W/m^2)
c --- salflx = salt flux (implied by fresh water flux)
c --- heatfx = heat flux through ice sheet (W/m^2)
c --- covice = ice coverage (rel.units)
c --- thkice = grid-box averaged ice thickness (m)
c --- temice = ice surface temperature
c --- brntop = top of depth interval over which to distribute brine
c --- brnbot = bottom of depth interval over which to distribute brine
c
      data thin/0.2/,rhoice/917./,thkmax/10./
     .    ,kice/2.04/,fusion/334.e3/,rate/5.e-6/,dtdflx/0.05/
c
c --- energy loan: add extra energy to the ocean to keep SST from dropping
c --- below tmelt in winter. return this borrowed energy to the 'energy bank'
c --- in summer as quickly as surflx > 0 allows.
c --- salt loan: analogous to energy loan.
c     
      total=0.
      totalx=0.
      dosmoo=.false.
      do 10 j=J_0,J_1
      do 10 l=1,isp(j)
      do 10 i=ifp(j,l),ilp(j,l)
      vrbos=i.eq.itest .and. j.eq.jtest
      pump=i.gt.equatn
cdiag total=total+thkice(i,j)*scp2(i,j)
c
      saldif=max(0.,saln(i,j,k1n)-10.)
      dpth=max(dp(i,j,k1n),thkmin*onem)
c
c --- calculate hypothetical mixed-layer temp due to diab. forcing 
      if (dp(i,j,k1n).le.0.) then
        write (*,'(i9,2i5,a)') nstep,i,j,'  zero mxlayr thickness'
        stop '(eice error)'
      end if
c
      odmsi(i,j)=0.
      salflx2(i,j)=0.
      tmxl=temp(i,j,k1n)+surflx(i,j)*delt1*g/(spcifh*dpth) 
c 
      if (saln(i,j,k1n).lt.0.) then
        write(*,'(i8,a,2i4,f8.4)')nstep,' neg s=',i,j,saln(i,j,k1n)
        saln(i,j,k1n)=0.
      endif
c
      tmelt=tfrez(saln(i,j,k1n),0.)
      if (tmxl.lt.tmelt) then
        borrow=min((tmelt-tmxl)*spcifh*dpth/(delt1*g),
     .           rate*fusion*rhoice)			! > 0 W/m^2
c
c --- add energy to bring tmxl back to tmelt (only if tmxl < tmelt)
c
        surflx(i,j)=surflx(i,j)+borrow
c --- corrections and generalisation (see Schmidt et al, 2004)
        odmsi(i,j)=-borrow/     ! odmsi: kg/m2/sec
     .       (Ei(tmelt,fsss*saln(i,j,k1n))-tmelt*spcifh) ! > 0
c       odhsi(i,j)=odmsi(i,j)*Ei(tmelt,fsss*saln(i,j,k1n)) ! ohmsi: J/m2/sec
c       odssi(i,j)=odmsi(i,j)*fsss*saln(i,j,k1n)           ! odssi: kg/m2/sec

        salflx2(i,j)=odmsi(i,j)*saln(i,j,k1n)*(1.-fsss)    ! g/m2/sec
c
      endif
c
c --- build up time integrals of surface fluxes
      eminpav(i,j)=eminpav(i,j)+oemnp(i,j)		!m/s
      surflav(i,j)=surflav(i,j)+surflx(i,j)		!W/m2
      salflav(i,j)=salflav(i,j)+salflx(i,j)+salflx2(i,j)!g/m2/s
      brineav(i,j)=brineav(i,j)+osalt(i,j)*1.e3		!g/m2/s 
     .                 -odmsi(i,j)*saln(i,j,k1n)*fsss	! ocn salt gain is +
      tauxav(i,j)=tauxav(i,j)+taux(i,j)	                ! N/m2
      tauyav(i,j)=tauyav(i,j)+tauy(i,j)	                ! N/m2
      oiceav(i,j)=oiceav(i,j)+oice(i,j)
c
c --- deposit brine from brntop to brnbot
c
      if (vrbos .and. pump) write (*,103) nstep,i,j,
     . '  entering eice:  temp    saln    dens    thkns    dpth',
     .  (k,temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn),
     .   dp(i,j,k+nn)/onem,p(i,j,k+1)/onem,k=1,kk)
 103  format (i9,2i5,a/(i34,2f8.3,f8.3,f8.2,f8.1))
c
      if (salflx2(i,j).gt.0. .and. pump) then
        bot=min(brnbot*onem,    pbot(i,j))
        top=min(brntop*onem,.75*pbot(i,j))
        thkinv=1./(bot-top)
        do 11 k=1,kk
        kn=k+nn
        p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
 11     saln(i,j,kn)=saln(i,j,kn)+salflx2(i,j)*delt1*g*thkinv
     .     *max(0.,min(p(i,j,k+1),bot)-max(p(i,j,k),top))
     .       /max(dp(i,j,kn),epsil)
      else
        salflx(i,j)=salflx(i,j)+salflx2(i,j)
      end if
c
      if (vrbos .and. pump) write (*,103) nstep,i,j,
     . '  exiting  eice:  temp    saln    dens    thkns    dpth',
     .  (k,temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn),
     .   dp(i,j,k+nn)/onem,p(i,j,k+1)/onem,k=1,kk)

 10   continue
c
      return
      end
c
c> Revision history
c>
c> June 2000 - conversion to SI units
c> July 2000 - switched sign convention for vertical fluxes (now >0 if down)
c> Aug. 2000 - revised ice surf.temp. (temice) calculation
c> Mar. 2001 - corrected error in -heatfx- calculation
c> Aug. 2001 - introduced mininum thickness -dpth- in loan calculation
c> Sep. 2001 - corrected -surflx- sign error in -tmxl- calculation
c> Feb. 2002 - re-introduced smoother to spread ice thicker than -thkmax-
c> Oct. 2004 - keep tract on ice mass only, add brine subsurface in SH
c> May. 2005 - made choice of brine deposition depth interval more flexible
