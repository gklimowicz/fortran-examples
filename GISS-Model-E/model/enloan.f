#include "hycom_mpi_hacks.h"
      subroutine enloan(m,n,mm,nn,k1m,k1n)
c
c --- "energy loan":
c --- heat ocean by forming ice to keep SST from dropping below Tmelt.
c --- this is a stripped version of hycom's original enloan routine.
c
      USE SEAICE,   only : fsss,tfrez,Ei
      USE HYCOM_DIM,only : isp,ifp,ilp,kk,idm,J_0,J_1,J_0H,J_1H
      USE HYCOM_SCALARS, only : thkmin,onem,nstep,delt1,g,spcifh
     & 	,stdsal,equatn,epsil,brntop,brnbot,itest,jtest,baclin
      USE HYCOM_ARRAYS,  only : temp,saln,dp,surflx,p,pbot,salflx
     &  ,brnflx,sqiflx,osalt,odhsi
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
c
      implicit none
      integer i,j,k,l,m,n,mm,nn,kn,k1m,k1n
      real tmelt,tmxl,dpth,transf,capcty,dpgt0,old,amount,top,bot,thkinv
      logical :: vrbos,pump=.false.
c
      do 10 j=J_0,J_1
      do 10 l=1,isp(j)
      do 10 i=ifp(j,l),ilp(j,l)
      vrbos = i==itest .and. j==jtest
      pump=i.gt.equatn		! activate pump in southern hemisphere
c
      if (dp(i,j,k1n).le.0.) then
        write (*,'(i9,2i5,a)') nstep,i,j,'  zero mxlayr thickness'
        stop '(enloan error)'
      end if
c
c --- if saln < 0 due to excessive freshwater input, bring up salt from below
c
      if (saln(i,j,k1n).lt.0.) then
        write(*,'(i8,a,2i4,f8.4)')nstep,' warning: neg S =',i,j
     .      ,saln(i,j,k1n)
        transf=-saln(i,j,k1n)*dp(i,j,k1n)
        saln(i,j,k1n)=0.
        k=1
        do while (transf.gt.0.)
          k=k+1
          if (k.gt.kk) then
            write(*,'(i8,a)') nstep,' enloan error: column S < 0'
          end if
          kn=k+nn
          capcty=min(saln(i,j,kn)*dp(i,j,kn),transf)
          dpgt0=max(epsil,dp(i,j,kn))
          saln(i,j,kn)=saln(i,j,kn)-capcty/dpgt0
          transf=transf-capcty
        end do
      endif
c
      dpth=max(dp(i,j,k1n),thkmin*onem)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- calculate hypothetical mixed-layer temp at end of current time step
      tmxl=temp(i,j,k1n)+surflx(i,j)*delt1*g/(spcifh*dpth)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- use SST to represent mixed-layer temp
c     tmxl=temp(i,j,k1n)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      tmelt=tfrez(saln(i,j,k1n),0.)
c
c --- form new ice if mixed-lyr temp falls below tmelt
c
      odhsi(i,j)=0.
      amount=0.
      brnflx(i,j)=0.
      sqiflx(i,j)=0.
      if (tmxl.lt.tmelt) then
        amount=(tmelt-tmxl)*spcifh*dpth/(delt1*g) ! > 0 W/m^2
        surflx(i,j)=surflx(i,j)+amount            ! > 0 W/m^2
c --- 'odhsi' are thermal energy per atmo time step transmitted to the ocean
        odhsi(i,j)=-amount                        ! < 0 W/m^2
c --- hycom doesn't allow fresh water transfer between ocean and ice
c --- during freezing/melting. the dynamic effect of freezing/melting must
c --- be represented by an "induced" salt flux opposing the freshwater flux.
c --- a newly formed amount of ice 'dmsi' (kg/m^2) sequesters the salt amount
c --- dssi = dmsi*fsss*stdsal/1000. the amount of salt left behind if fsss<1,
c --- (1-fsss)*dmsi*stdsal/1000, becomes the induced salt flux (g/m^2/sec)
c --- (1-fsss)*dmsi*stdsal/dt representing brine ejected during ice fomation.
c --- sign convention: ice mass gain -> dmsi>0
c --- sign convention: ice salt gain -> dssi>0
c --- sign convention: ocn salt gain -> brnflx>0

css     brnflx(i,j)=amount*(1.-fsss)*stdsal/
        brnflx(i,j)=amount*(1.-fsss)*saln(i,j,k1n)/
     .    (tmelt*spcifh-Ei(tmelt,fsss*saln(i,j,k1n)))	! >0  g/m^2/sec
css     sqiflx(i,j)=-amount*   fsss *stdsal/
        sqiflx(i,j)=-amount*   fsss *saln(i,j,k1n)/
     .    (tmelt*spcifh-Ei(tmelt,fsss*saln(i,j,k1n)))	! <0  g/m^2/sec
        if (tmxl.lt.tmelt-1. .or. vrbos)
     .   print 100, nstep,i,j,' (enloan) -borrow- T=',temp(i,j,k1n),
     .     ' W/m2',amount,' brine',brnflx(i,j)*1.e6,' sequst',
     .     sqiflx(i,j)*1.e6
 100    format (i6,2i4,a,f6.2,1x,6(a,"=",es9.2))

        if (pump) then
c --- distribute brine ejected during ice formation in interval (brntop,brnbot)
          bot=min(brnbot*onem,    pbot(i,j))
          top=min(brntop*onem,.75*pbot(i,j))
          thkinv=1./(bot-top)
          do 11 k=1,kk
          kn=k+nn
          p(i,j,k+1)=p(i,j,k)+dp(i,j,kn)
 11       saln(i,j,kn)=saln(i,j,kn)+brnflx(i,j)*delt1*g*thkinv
     .     *max(0.,min(p(i,j,k+1),bot)-max(p(i,j,k),top))
     .       /max(dp(i,j,kn),epsil)
        else
          salflx(i,j)=salflx(i,j)+brnflx(i,j)	! no pump
        end if
      endif
c
      if (vrbos) write (*,103) nstep,i,j,
     .  '  exiting enloan   surflx increment =',amount,
     .  odhsi(i,j),
     .  Ei(tmelt,fsss*saln(i,j,k1n)),
     .  '       temp    saln    thkns    dpth',
     .  (k,temp(i,j,k+nn),saln(i,j,k+nn),
     .  dp(i,j,k+nn)/onem,p(i,j,k+1)/onem,k=1,kk)
 103  format (i8,2i5,a,3es9.2/a/(i3,2f8.3,2f8.1))

 10   continue
c
      return
      end
c
c> Revision history
c>
c> July 2017 - created stripped version of hycom's enloan
