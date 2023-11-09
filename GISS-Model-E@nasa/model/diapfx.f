      subroutine diapfl(m,n,mm,nn,k1m,k1n)
c
c --- hycom version 0.9.2
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT
      USE HYCOM_DIM, only : jj,kk,isp,ifp,ilp,idm,kdm,ntrcr
     &     ,jchunk, J_0, J_1
      USE HYCOM_SCALARS, only : diapyc,nstep,dotrcr
     &     ,onemm,g,baclin,onem, onecm ,onemu
     &     ,epsil,mixfrq,sigjmp,thref,acurcy,diapyn
     &     ,itest,jtest
      USE HYCOM_ARRAYS
      implicit none
c
      integer i,j,k,l,m,n,mm,nn,kn,k1m,k1n
      real flxu(kdm),flxl(kdm),pdot(kdm),flngth(kdm),clip(kdm),
     .     ennsq,alfa,beta,q,qmin,qmax,amount,salt,froglp,small,delp,
     .     trflxu(0:kdm+1,ntrcr),trflxl(0:kdm+1,ntrcr),cliptr(ntrcr),
     .     tflxu(0:kdm+1),tflxl(0:kdm+1),clipt,
     .     sflxu(0:kdm+1),sflxl(0:kdm+1),clips,
     .     told(2),sold(2),trold(2,ntrcr),scale,
     .     totem,tosal,totra,tndcyt,tndcys,tndtra
      integer kmin,kmax,ka,kan,nt,iter
      character text*20
      logical event,vrbos
      data small/1.e-22/
      real sigocn,dsigdt,dsigds
      external sigocn,dsigdt,dsigds
c
      if (diapyc.eq.0. .or. mod(nstep,mixfrq).gt.1) return
c
c --- ----------------
c --- diapycnal mixing
c --- ----------------
c
c --- if mixfrq > 1, apply mixing algorithm to both time levels
      froglp=max(2,mixfrq)
c
ccc   salt=0.
c
      do 31 j=J_0, J_1
      do 31 l=1,isp(j)
      do 31 i=ifp(j,l),ilp(j,l)
      vrbos=i.eq.itest .and. j.eq.jtest
c
c --- t/s conservation diagnostics (optional):
      totem=0.
      tosal=0.
      totra=0.
      scale=1.e-99
      do k=1,kk
        kn=k+nn
        totem=totem+temp(i,j,kn)*dp(i,j,kn)
        tosal=tosal+saln(i,j,kn)*dp(i,j,kn)
        if (dotrcr) then
          totra=totra+tracer(i,j,k,1)*dp(i,j,kn)
          scale=scale+abs(tracer(i,j,k,1))
        end if
      end do
c
      do 33 k=1,kk
 33   p(i,j,k+1)=p(i,j,k)+dp(i,j,k+nn)
c
      sold(1)=saln(i,j,kk+nn)
      told(1)=temp(i,j,kk+nn)
      tflxl(   0)=0.
      tflxu(kk+1)=0.
      sflxl(   0)=0.
      sflxu(kk+1)=0.
      if (dotrcr) then
        trold(1,:)=tracer(i,j,kk,:)
        trflxl(   0,:)=0.
        trflxu(kk+1,:)=0.
      end if
c
      if (vrbos) write (*,103) nstep,i,j,
     . '  entering diapfl:  temp    saln    dens    thkns   tracer',
     .  (k,temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn),
     .   dp(i,j,k+nn)/onem,tracer(i,j,k,1),k=1,kk)
 103  format(i9,2i5,a/(33x,i3,3f8.3,f8.2,f9.4))
c
      kmin=kk+1
      kmax=1
c
      do 36 k=2,kk
      kn=k+nn
c
c --- locate lowest mass-containing layer and upper edge of stratified region
      if (p(i,j,k).lt.p(i,j,kk+1)-onecm)  then
        kmax=k
        if (kmin.eq.kk+1 .and. th3d(i,j,kn).gt.th3d(i,j,kn-1)+.1*sigjmp)
     .      kmin=k-1
      end if
 36   continue
c
      if (vrbos) write (*,'(i9,2i5,a,2i5)') nstep,i,j,' kmin,kmax =',
     .  kmin,kmax
c
c --- find buoyancy frequency for each layer
c
      do 43 k=2,kk-1
      kn=k+nn
      if (k.gt.kmin .and. k.lt.kmax) then
c --- ennsq = buoy.freq.^2 / g^2
        ennsq=max(0.,min(th3d(i,j,kn+1)-th3d(i,j,kn  ),
     .                   th3d(i,j,kn  )-th3d(i,j,kn-1)))
     .    /max(p(i,j,k+1)-p(i,j,k),onemu)
c --- store (exch.coeff x buoy.freq.^2 / g x time step) in -flngth-
c --- (dimensions of flngth: length in pressure units)
c -----------------------------------------------------------------------
c --- use the following if exch.coeff. = diapyc / buoyancy frequency
ccc        flngth(k)=diapyc*sqrt(ennsq) * baclin*froglp * onem
c -----------------------------------------------------------------------
c --- use the following if exch.coeff. = diapyc
ccc        flngth(k)=diapyc*ennsq*g * baclin*froglp * onem
c -----------------------------------------------------------------------
        flngth(k)=max(diapyn*sqrt(ennsq),diapyc*ennsq*g)  ! max of two
     .            * baclin*froglp * onem
c
      end if
 43   continue
c
c --- find fluxes at the upper and lower interface of each layer
c --- (compute only the part common to t/s/mass fluxes)
c
      do 37 k=1,kk
      kn=k+nn
      flxu(k)=0.
      flxl(k)=0.
      if (k.gt.kmin .and. k.lt.kmax) then
        alfa=-thref*dsigdt(temp(i,j,kn),saln(i,j,kn))
        beta= thref*dsigds(temp(i,j,kn),saln(i,j,kn))
c
        flxu(k)=flngth(k)/
     .    max(beta*(saln(i,j,kn)-saln(i,j,kn-1))
     .       -alfa*(temp(i,j,kn)-temp(i,j,kn-1)),small)
        flxl(k)=flngth(k)/
     .    max(beta*(saln(i,j,kn+1)-saln(i,j,kn))
     .       -alfa*(temp(i,j,kn+1)-temp(i,j,kn)),small)
      end if
 37   continue
c
c --- determine mass flux -pdot- implied by t/s fluxes.
c
      do 40 k=1,kk
      clip(k)=1.
      pdot(k)=0.
      if (k.gt.kmin .and. k.le.kmax)
     .  pdot(k)=flxu(k)-flxl(k-1)
 40   continue
c
c --- now clip mass fluxes to prevent dp < 0
c
      event=.false.
      do iter=1,5      !  go up and down the column repeatedly
c
      do 42 k=kk*mod(iter,2)+ 2*mod(iter-1,2),
     .         2*mod(iter,2)+kk*mod(iter-1,2),
     .          -mod(iter,2)+   mod(iter-1,2)
      kn=k+nn
      if (k.gt.kmin .and. k.le.kmax) then
        if (pdot(k).gt.0.) then
          q=max(0.,.5*dp(i,j,kn-1))/flxu(k  )
          if (q.lt.clip(k  )) then
            event=.true.
            clip(k  )=q
          end if
c
          if (vrbos .and. clip(k  ).lt.1.)
     .     write (*,'(i3,a,5es10.2)') k,'  pdot,dp,flxu,flxl,clip=',
     .      pdot(k)/onem,dp(i,j,kn-1)/onem,flxu(k)/onem,
     .       flxl(k-1)/onem,clip(k  )
c
        else if (pdot(k).lt.0.) then
          q=max(0.,.5*dp(i,j,kn  ))/flxl(k-1)
          if (q.lt.clip(k-1)) then
            event=.true.
            clip(k-1)=q
          end if
c
          if (vrbos .and. clip(k-1).lt.1.)
     .     write (*,'(i3,a,5es10.2)') k,'  pdot,dp,flxu,flxl,clip=',
     .      pdot(k)/onem,dp(i,j,kn  )/onem,flxu(k)/onem,
     .       flxl(k-1)/onem,clip(k-1)
c
        end if
      end if
 42   continue
c
      do 44 k=1,kk
      kn=k+nn
c
      if (vrbos .and. clip(k  ).lt.1.) then
 101   format (i3,a,(2es10.2,2x))
       write (*,101) k-1,
     .  ' flxu,flxl,clip:',flxu(k-1)/onem,flxl(k-1)/onem,
     .   flxu(k-1)*clip(k-1)/onem,flxl(k-1)*clip(k-1)/onem,clip(k-1)
       write (*,101) k  ,flxu(k  )/onem,flxl(k  )/onem,
     .   flxu(k  )*clip(k  )/onem,flxl(k  )*clip(k  )/onem,clip(k  )
       write (*,101) k+1,flxu(k+1)/onem,flxl(k+1)/onem,
     .   flxu(k+1)*clip(k+1)/onem,flxl(k+1)*clip(k+1)/onem,clip(k+1)
      end if
c
      flxu(k)=flxu(k)*clip(k)
      flxl(k)=flxl(k)*clip(k)
      clip(k)=1.
      if (k.gt.kmin .and. k.le.kmax)
     .  pdot(k)=flxu(k)-flxl(k-1)
 44   continue
c
      if (.not.event) exit
      end do       !  iter
c
c --- convert flxu,flxl into actual t/s (and tracer) fluxes
c
      do 35 k=1,kk
      kn=k+nn
      tflxu(k)=0.
      tflxl(k)=0.
      sflxu(k)=0.
      sflxl(k)=0.
      if (dotrcr) then
        trflxu(k,:)=0.
        trflxl(k,:)=0.
      end if
      if (k.gt.kmin .and. k.lt.kmax) then
        tflxu(k)=flxu(k)*temp(i,j,kn-1)
        tflxl(k)=flxl(k)*temp(i,j,kn+1)
        sflxu(k)=flxu(k)*saln(i,j,kn-1)
        sflxl(k)=flxl(k)*saln(i,j,kn+1)
        if (dotrcr) then
          trflxu(k,:)=flxu(k)*tracer(i,j,k-1,:)
          trflxl(k,:)=flxl(k)*tracer(i,j,k+1,:)
        end if
      end if
 35   continue
c
      if (dotrcr) cliptr(:)=0.
      clipt=0.
      clips=0.
c
c --- update interface pressure and layer temperature/salinity
      do 39 k=kk,1,-1
      kn=k+nn
      ka=max(1,k-1)
      kan=ka+nn
      sold(2)=sold(1)
      sold(1)=saln(i,j,kn)
      told(2)=told(1)
      told(1)=temp(i,j,kn)
      if (dotrcr) then
        trold(2,:)=trold(1,:)
        trold(1,:)=tracer(i,j,k,:)
      end if
c
      dpold(i,j,k)=dp(i,j,kn)
      p(i,j,k)=p(i,j,k)-pdot(k)
      dp(i,j,kn)=p(i,j,k+1)-p(i,j,k)
      if (dp(i,j,kn).lt.-onemm) then
        print '(a,3i5,es10.2)','diapfl: dp<0 at',i,j,k,dp(i,j,kn)/onem
        stop
      end if
c
      if (k.ge.kmin .and. k.le.kmax) then
        delp=dp(i,j,kn)
        if (delp.gt.0.) then
          amount=temp(i,j,kn)*dpold(i,j,k)
     .      -(tflxu(k+1)-tflxu(k)+tflxl(k-1)-tflxl(k))
          q=amount
          qmax=max(temp(i,j,kan),told(1),told(2))
          qmin=min(temp(i,j,kan),told(1),told(2))
          amount=max(qmin*delp,min(amount,qmax*delp))
          clipt=clipt+(q-amount)
          temp(i,j,kn)=amount/delp
c
          amount=saln(i,j,kn)*dpold(i,j,k)
     .      -(sflxu(k+1)-sflxu(k)+sflxl(k-1)-sflxl(k))
          q=amount
          qmax=max(saln(i,j,kan),sold(1),sold(2))
          qmin=min(saln(i,j,kan),sold(1),sold(2))
          amount=max(qmin*delp,min(amount,qmax*delp))
          clips=clips+(q-amount)
          saln(i,j,kn)=amount/delp
c
          if (dotrcr) then
            do nt=1,ntrcr
              amount=tracer(i,j,k,nt)*dpold(i,j,k)
     .            -(trflxu(k+1,nt)-trflxu(k,nt)
     .             +trflxl(k-1,nt)-trflxl(k,nt))
              q=amount
              qmax=max(tracer(i,j,ka,nt),trold(1,nt),trold(2,nt))
              qmin=min(tracer(i,j,ka,nt),trold(1,nt),trold(2,nt))
              amount=max(qmin*delp,min(amount,qmax*delp))
              cliptr(nt)=cliptr(nt)+(q-amount)
              tracer(i,j,k,nt)=amount/delp
            end do
          end if			!  dotrcr
        end if				!  delp > 0
      end if
 39   continue
c
      if (dotrcr) cliptr(:)=cliptr(:)/pbot(i,j)
      clipt=clipt/pbot(i,j)
      clips=clips/pbot(i,j)
c
      do 41 k=1,kk
      kn=k+nn
c
c --- restore 'clipped' t/s amount to column
      temp(i,j,kn)=temp(i,j,kn)+clipt
      saln(i,j,kn)=saln(i,j,kn)+clips
      th3d(i,j,kn)=sigocn(temp(i,j,kn),saln(i,j,kn))
      if (dotrcr) tracer(i,j,k,:)=tracer(i,j,k,:)+cliptr(:)
c
      diaflx(i,j,k)=diaflx(i,j,k)+(dp(i,j,kn)-dpold(i,j,k)) ! diapyc.flux
c --- make sure p is computed from dp, not the other way around (roundoff!)
 41   p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
c
c --- t/s conservation diagnostics (optional):
      tndcys=-tosal
      tndcyt=-totem
      tndtra=-totra
      do k=1,kk
        kn=k+nn
        tndcys=tndcys+saln(i,j,kn)*dp(i,j,kn)
        tndcyt=tndcyt+temp(i,j,kn)*dp(i,j,kn)
        if (dotrcr) tndtra=tndtra+tracer(i,j,k,1)*dp(i,j,kn)
      end do
      if (abs(tndcyt).gt.acurcy*10.*pbot(i,j))
     .  write (*,100) i,j,
     .   '  diapfl - bad temp.intgl.',totem,tndcyt,clipt
        if (abs(tndcys).gt.acurcy*35.*pbot(i,j))
     .  write (*,100) i,j,
     .   '  diapfl - bad saln.intgl.',tosal,tndcys,clips
      if (dotrcr) then
        if (abs(tndtra)*kk.gt.acurcy*scale*pbot(i,j))
     .  write (*,100) i,j,
     .   '  diapfl - bad trcr.intgl.',totra,tndtra,cliptr(1)
      end if
 100  format(2i5,a,1p,e16.8,2e13.5)
c
      if (vrbos) write (*,103) nstep,i,j,
     . '  exiting  diapfl:  temp    saln    dens    thkns   tracer',
     .  (k,temp(i,j,k+nn),saln(i,j,k+nn),th3d(i,j,k+nn),
     .   dp(i,j,k+nn)/onem,tracer(i,j,k,1),k=1,kk)
c
 31   continue
c
ccc   write (*,'(i9,7x,1p,e9.2,a)') nstep,salt*1.e-6/g,
ccc  .  ' kg salt added in diapfl'
c
      call pardpudpv(nn)
c
      if (dotrcr .and. AM_I_ROOT())
     .  write (*,'(a)') 'tracer diapycnal mixing done'
      return
      end
c
c
c> Revision history:
c>
c> Mar. 2000 - conversion to SI units
c> May  2000 - converted T/S advection equations to flux form
c> Apr. 2001 - eliminated stmt_funcs.h
c> Sep. 2003 - added code to return 'clipped' tracer amount to column (loop 41)
c> Sep. 2003 - added logical switch to enable/disable tracer diffusion
c> Feb. 2005 - added multiple tracer capability
c> July 2007 - less restrictive clipping of mass fluxes -pdot-
c> July 2007 - changed order of -i- and -k- loops to reduce indexing complexity
