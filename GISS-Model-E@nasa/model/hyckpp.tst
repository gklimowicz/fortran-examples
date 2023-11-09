      subroutine OCEANS
c
c --- ------------------------------
c --- MICOM-based hybrid ocean model
c
c ---     v e r s i o n    0.9
c --- ------------------------------
c
      USE FLUXES, only : e0,prec,eprec,evapor,flowo,eflowo,dmua,dmva
     . ,erunosi,runosi,srunosi,runpsi,srunpsi,dmui,dmvi,dmsi,dhsi,dssi
     . ,gtemp,sss,mlhc,ogeoza,uosurf,vosurf,MELTI,EMELTI,SMELTI
     . ,gmelt,egmelt,solar
      USE SEAICE_COM, only : rsi,msi
      USE SEAICE, only : fsss,tfrez
      USE GEOM, only : dxyp
      USE MODEL_COM, only : focean
      USE CONSTANT, only : lhm,shi,shw
      USE MODEL_COM, only: dtsrc
     *  ,itime,iyear1,nday,jdendofm,jyear,jmon,jday,jdate,jhour,aMON
c
      implicit none
c
#include "dimensions.h"
#include "dimension2.h"
#include "common_blocks.h"
#include "cpl.h"
#include "a2o.h"
#include "kpp.h"
c
      real sum,coord,x,x1,totl,sumice,fusion,saldif,sofsig,tf
     .    ,sigocn,kappaf,check,apehyc,pechg_hyc_bolus
     .    ,hyc_pechg1,hyc_pechg2,q,sum1,sum2
      integer jj1,no,index,nflip,mo0,mo1,mo2,mo3,rename,iatest,jatest
     .       ,OMP_GET_NUM_THREADS,io,jo,maska(iia,jja),nsub
      external rename
      logical master,slave,diag_ape
      character util(idm*jdm+14)*2,charac(20)*1,string*20,
     .          flnm*60
      data charac/'1','2','3','4','5','6','7','8','9','0',
     .            'A','B','C','D','E','F','G','H','I','J'/
      data iatest,jatest/9,7/
c     data iatest,jatest/33,40/           ! Iceland
      data nflip/0/
c
      real*4 user_time,sys_time,total_time,avg_time
      real*4 cnuity_time,cnuity_total_time,cnuity_avg_time
      real*4 tsadvc_time,tsadvc_total_time,tsadvc_avg_time
      real*4 momtum_time,momtum_total_time,momtum_avg_time
      real*4 barotp_time,barotp_total_time,barotp_avg_time
      real*4 thermf_time,thermf_total_time,thermf_avg_time
      real*4 enloan_time,enloan_total_time,enloan_avg_time
      real*4 mxlayr_time,mxlayr_total_time,mxlayr_avg_time
      real*4 hybgen_time,hybgen_total_time,hybgen_avg_time
      real*4 agcm_time
      integer after,before,rate,bfogcm
      data cnuity_total_time,tsadvc_total_time /0.0,0.0/
      data momtum_total_time,barotp_total_time /0.0,0.0/
      data thermf_total_time,enloan_total_time /0.0,0.0/
      data mxlayr_total_time,hybgen_total_time /0.0,0.0/
      data cnuity_time,tsadvc_total_time /0.0,0.0/
      data momtum_time,barotp_total_time /0.0,0.0/
      data thermf_time,enloan_total_time /0.0,0.0/
      data mxlayr_time,hybgen_total_time /0.0,0.0/
c
      real osst(idm,jdm),osss(idm,jdm),osiav(idm,jdm)
     . ,oogeoza(idm,jdm),usf(idm,jdm),vsf(idm,jdm)
     . ,utila(iia,jja)

#include "state_eqn.h"
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- initiate named-pipe comparison utility
c --- (see pipe.f for instructions)
ccc      inquire(file='master',exist=master)
ccc      inquire(file='slave',exist=slave)
ccc      if ((master .and. slave) .or. (.not.master .and. .not.slave))
ccc     .   stop '(master/slave ambiguity)'
ccc      if (master) call pipe_init(.true.)
ccc      if (slave) call pipe_init(.false.)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      call getdte(Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon)
      write(*,'(a,i8,7i5,a)')'chk =',
     .    Itime,Nday,Iyear1,Jyear,Jmon,Jday,Jdate,Jhour,amon
c
      if (mod(jhour,nhr).eq.0.and.mod(itime,nday/24).eq.0) then
c$OMP PARALLEL DO
        do 101 ia=1,iia
        do 101 ja=1,jja
          ataux(ia,ja)=0.
          atauy(ia,ja)=0.
        aflxa2o(ia,ja)=0.
          aemnp(ia,ja)=0.
          asalt(ia,ja)=0.
          aice(ia,ja)=0.
         austar(ia,ja)=0.
         aswflx(ia,ja)=0.
 101    continue
c$OMP END PARALLEL DO
      endif
c
c$OMP PARALLEL DO
c --- accumulate agcm fields and check if ogcm will be called
      do 102 ia=1,iia
      do 102 ja=1,jja
      maska(ia,ja)=0
      if (focean(ia,ja).eq.0.) goto 102
      maska(ia,ja)=1
c --- accumulate 
      aemnp(ia,ja)=aemnp(ia,ja)                               ! kg/m*m => m/s
     .+((prec(ia,ja)-evapor(ia,ja,1))*(1.-rsi(ia,ja))         ! open water
     .+(flowo(ia,ja)+gmelt(ia,ja)+melti(ia,ja))/(dxyp(ja)*focean(ia,ja))!ocn/ice
     .+(runosi(ia,ja)+runpsi(ia,ja))*rsi(ia,ja))*thref        ! ice
     .                                /(3600.*real(nhr))
      aflxa2o(ia,ja)=aflxa2o(ia,ja)                           ! J/m*m => W/m*m
     . +((e0(ia,ja,1)+eprec(ia,ja))*(1.-rsi(ia,ja))           ! ocean water
     . +(eflowo(ia,ja)+egmelt(ia,ja)+emelti(ia,ja))
     .                              /(dxyp(ja)*focean(ia,ja)) ! ocn or ice
css  . +(erunosi(ia,ja)+erunpsi(ia,ja))*rsi(ia,ja))           ! ice
     . + erunosi(ia,ja)*rsi(ia,ja))                           ! ice
     .                                 /(3600.*real(nhr))
      asalt(ia,ja)=asalt(ia,ja)                            ! kg/m*m/sec salt
     .+((srunosi(ia,ja)+srunpsi(ia,ja))*rsi(ia,ja)         !
     .   +smelti(ia,ja)/(dxyp(ja)*focean(ia,ja)))          !
     .                             /(3600.*real(nhr))
       aice(ia,ja)= aice(ia,ja) + rsi(ia,ja)*dtsrc/(real(nhr)*3600.)
c --- dmua on B-grid, dmui on C-grid; Nick aug04
      ataux(ia,ja)=ataux(ia,ja)+(dmua(ia,ja,1)+dmui(ia,ja))    ! scaled by rsi 
     .                                     /(3600.*real(nhr))  ! kg/ms => N/m*m
      atauy(ia,ja)=atauy(ia,ja)+(dmva(ia,ja,1)+dmvi(ia,ja))
     .                                     /(3600.*real(nhr))  ! kg/ms => N/m*m
      austar(ia,ja)=austar(ia,ja)+(
     . sqrt(sqrt((dmua(ia,ja,1)+dmui(ia,ja))**2
     .          +(dmva(ia,ja,1)+dmvi(ia,ja))**2)/dtsrc*thref)) ! sqrt(T/r)=>m/s
     .                               *dtsrc/(real(nhr)*3600.)
      aswflx(ia,ja)=aswflx(ia,ja)+(solar(1,ia,ja)*(1.-rsi(ia,ja))!J/m*m=>W/m*m
     .                            +solar(3,ia,ja)*    rsi(ia,ja))
     .                                 /(3600.*real(nhr))
 102  continue
c$OMP END PARALLEL DO
c
      nsavea=nsavea+1
      if (mod(jhour,nhr).gt.0.or.mod(itime,nday/24).eq.0) return
      print *,' chk agcm saved over hr=',nsavea*24/nday
      nsavea=0
c
      call veca2o(ataux,atauy,taux,tauy)          !wind stress
      call flxa2o(aice,oice)                      !ice coverage
      call flxa2o(aflxa2o,oflxa2o)                !heatflux everywhere
      call flxa2o(asalt,osalt)                    !saltflux from SI
      call flxa2o(aemnp,oemnp)                    !E - P everywhere
      call flxa2o(austar,ustar)                   !friction velocity
      call flxa2o(aswflx,sswflx)                  ! shortwave flux
c
ccc      if (nstep.eq.0) then
cccc$OMP PARALLEL DO
ccc        do 17 j=1,jj
ccc        do 17 l=1,isp(j)
ccc        do 17 i=ifp(j,l),ilp(j,l)
ccc         if (oice(i,j).gt. .95 .and. temp(i,j,1).gt. .0)
ccc     .     write(*,'(f6.2,a,2i4,f6.2,a,f6.2)')
ccc     .   oice(i,j),' bad sst at=',i,j,temp(i,j,1),' new s=',
ccc     . sofsig(th3d(i,j,1)+thbase,-1.8)
ccc        if (oice(i,j).gt. .5) then
ccc          do k=1,kk
ccc          temp(i,j,k   )=min(temp(i,j,k),-1.8)
ccc          saln(i,j,k   )=sofsig(th3d(i,j,k)+thbase,temp(i,j,k))
ccc          temp(i,j,k+kk)=temp(i,j,k)
ccc          saln(i,j,k+kk)=saln(i,j,k)
ccc          enddo
ccc        endif
ccc 17     continue
cccc$OMP END PARALLEL DO
ccc      endif
c
      call system_clock(before)
      call system_clock(count_rate=rate)
      bfogcm=before
      agcm_time = real(bfogcm-afogcm,4)/real(rate,4)
      write (lp,99009) agcm_time,' sec for AGCM step ',nstep
c
      if (nstep.eq.0 .or.nstep.eq.nstep0) then
c
c     Print seconds elapsed since startup (should be almost zero)
c     call system_clock(before)
c     bfogcm=before
      call system_clock(count_rate=rate)
      call system_clock(after)
      user_time = real(after-before,4)/real(rate,4)
css   write (lp,99009) user_time, ' Time 0 HYCOM starting up'
99009 format (f12.3,a,i8)
c
c     agcm_time = real(bfogcm-afogcm,4)/real(rate,4)
c     write (lp,99009) agcm_time,' sec for AGCM step ',nstep
c
c$OMP PARALLEL SHARED(mo0)
      mo0=OMP_GET_NUM_THREADS()
c$OMP END PARALLEL
      call OMP_SET_DYNAMIC(.false.)
ccc   call OMP_SET_NUM_THREADS(max(mo0,16))
c$OMP PARALLEL SHARED(mo1)
ccc   mo1=OMP_GET_NUM_THREADS()
c$OMP END PARALLEL
c     write (lp,'(2(a,i5))') ' hycom thread count',mo0
c     write (lp,'(2(a,i5))') ' hycom thread count',mo0,' changed to',mo1
ccc     . ,' =======  number of threads:',mo1,' ======='
c
      lp=6
c
c --- compute eqn.of state check values
c
      check=36.876506               ! T:[-2:32],S:[16:38]
      if (abs(sigocn(4.,35.)-check).gt..0001) then
       write (lp,'(/2(a,f11.6))') 'error -- sigocn(t=4,s=35) should be',
     . check,', not',sigocn(4.,35.)
css     stop
      end if
c --- reference: t= 0.0, s=34.0, p=0 bar,kap(4.5,34.5,1.e7)=  0.11954594
c     check=.11954594
c     check=0.
c     if (abs(kappaf(4.5,34.5,1.e7)-check).gt..00001) then
c       write (lp,'(/a,2(a,f12.8))') 'error: kappa(4.5,34.5,10^7)',
c    .  '  should be',check,', not',kappaf(4.5,34.5,1.e7)
c       stop
c     end if
      mixfrq=int(43200./baclin+.0001)
c
      write (lp,109) thkdff,temdff,veldff,viscos,diapyc,vertmx
 109  format (' turb. flux parameters:',1p/
     .  ' thkdff,temdff,veldff =',3e9.2/
     .  ' viscos,diapyc,vertmx =',3e9.2)
      write (lp,'(a,1p,e11.3)') 'thbase =',thbase
c
c --- 'lstep' = number of barotropic time steps per baroclinic time step.
c --- lstep   m u s t   be even.
c
      lstep=baclin/batrop
      lstep=2*((lstep+1)/2)
      dlt=baclin/lstep
      nstepi=real(nhr)*3600./baclin + .0001
      write (lp,'(i4,'' barotropic steps per baroclinic time step'')')
     .  lstep
      write (lp,'(''ogcm exchange info w. agcm every step/hr'',2i5)')
     .  nstepi,nhr
c
c --- set up parameters defining the geographic environment
c
css   call geopar                ! moved to agcm for zero start or below
css   call inicon                ! moved to agcm
      call inikpp
c
      watcum=0.
      empcum=0.
c
        write (lp,'(''restart date in ogcm/agcm '',2i5,'' hr '',2i6
     .   /'' iyear1/jyear='',2i5)')
     .  int((nstep0+nstepi-1)*baclin)/(3600*24),itime/nday
     . ,int((nstep0+nstepi-1)*baclin)/3600,itime*24/nday,iyear1,jyear
      if ((nstep0+nstepi-1)*baclin/3600 .ne. itime*24/nday) then
        write (lp,'(/(a,i9,a,i9,a,f7.1,a))')'chk model start step',nstep0
     . ,' changed to: '
     . ,int((itime*24/nday+(iyear1-jyear)*8760)*3600/baclin)-nstepi+1
     . ,' (day '
     .,(int((itime*24/nday+(iyear1-jyear)*8760)*3600/baclin)-nstepi+1)
     .  *baclin/86400,')'
        nstep0=int((itime*24/nday+(iyear1-jyear)*8760)*3600/baclin)
     .                                                  -nstepi+1
        nstep=nstep0
        time0=nstep0*baclin/86400
      else
        write (lp,'(/(a,i9))') 'chk model starts at steps',nstep0
      endif
c
      endif
c
      if (nstepi.le.0) stop 'wrong nstepi'
c
c     print *,' shown below oice'
c     call zebra(oice,iio,iio,jjo)
c     print *,' shown below aflxa2o'
c     call zebra(aflxa2o,iia,iia,jja)
c
c$OMP PARALLEL DO
      do 202 j=1,jj
      do 202 l=1,isp(j)
      do 202 i=ifp(j,l),ilp(j,l)
      osst(i,j)=0.
      osss(i,j)=0.
      osiav(i,j)=0.
 202  continue
c$OMP END PARALLEL DO
c
c --- ---------------------
c --- sub loop starts here
c --- ---------------------
c
c --- letter 'm' refers to mid-time level (example: dp(i,j,km) )
c --- letter 'n' refers to old and new time level
c
      nsaveo=0
      do 15 nsub=1,nstepi
      m=mod(nstep  ,2)+1
      n=mod(nstep+1,2)+1
      mm=(m-1)*kk
      nn=(n-1)*kk
      k1m=1+mm
      k1n=1+nn
      nstep=nstep+1
      time=time0+(nstep-nstep0)*baclin/86400.
c
      diagno=.false.
      if (JDendOfM(jmon).eq.jday.and.Jhour.eq.24.and.nsub.eq.nstepi) 
     .                                    diagno=.true. ! end of month
c
      sum1=0.
      sum2=0.
      do 501 k=1,kk
      km=k+mm
      kn=k+nn
      do 501 j=1,jj
      do 501 l=1,isp(j)
      do 501 i=ifp(j,l),ilp(j,l)
      sum1=sum1+temp(i,j,km)*dp(i,j,km)*scp2(i,j)
      sum2=sum2+temp(i,j,kn)*dp(i,j,kn)*scp2(i,j)
      write(333,'(a,i3,4f7.1,f8.1)') 'i,j,k=',i,j,k
     . ,dp(i,j,k+mm)/onem,temp(i,j,k+mm)
     . ,dp(i,j,k+nn)/onem,temp(i,j,k+nn)
 501  continue
c
      write(*,'(a,i6,2f12.2)')
     .   ' htst1=',nstep-nstep0,sum1*spcifh*1000/area*1.e-8
     .  ,sum2*spcifh*1000/area*1.e-8
c
      diag_ape=diagno
      if (nstep.eq.1) diag_ape=.true.
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- available potential energy diagnostics:
      if (diag_ape)
     . write (501,'(f9.1,a,1p,e15.7)') time,'  APE (J):',
     .  hyc_pechg1(dp(1,1,k1n),th3d(1,1,k1n),31)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 103  format (f9.1,a,-12p,f9.3,' TW')
c
      call system_clock(before)    ! time elapsed since last system_clock
      cnuity_time = 0.0            ! Zero per-call timer
      call cnuity(m,n,mm,nn,k1m,k1n)
      call system_clock(after)     ! time elapsed since last system_clock
      cnuity_time = real(after-before,4)/real(rate,4)  ! time call used
      cnuity_total_time = cnuity_total_time + cnuity_time ! subtotal
      total_time = total_time + cnuity_time               ! total
c     call sstbud(0,'  initialization',temp(1,1,k1n))
c
ccc      write (string,'(a12,i8)') 'cnuity, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     .  q=hyc_pechg1(dp(1,1,k1m),th3d(1,1,k1m),32)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      before = after
      tsadvc_time = 0.0
      call tsadvc(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      tsadvc_time = real(after-before,4)/real(rate,4)
      tsadvc_total_time = tsadvc_total_time + tsadvc_time
      total_time = total_time + tsadvc_time
c     call sstbud(1,' horiz.advection',temp(1,1,k1n))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to time smoothing: ',
     .  hyc_pechg2(dp(1,1,k1m),th3d(1,1,k1m),32)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ccc      write (string,'(a12,i8)') 'tsadvc, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      momtum_time = 0.0
      call momtum(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      momtum_time = real(after-before,4)/real(rate,4)
      momtum_total_time = momtum_total_time + momtum_time
      total_time = total_time + momtum_time
c
ccc      write (string,'(a12,i8)') 'momtum, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      barotp_time = 0.0
      call barotp(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      barotp_time = real(after-before,4)/real(rate,4)
      barotp_total_time = barotp_total_time + barotp_time
      total_time = total_time + barotp_time
c
ccc      write (string,'(a12,i8)') 'barotp, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
c     before = after
c     convec_time = 0.0
c     call convec(m,n,mm,nn,k1m,k1n)
c     call system_clock(after)
c     convec_time = real(after-before,4)/real(rate,4)
c     convec_total_time = convec_total_time + convec_time
c     total_time = total_time + convec_time
c     call sstbud(2,' convec.adjustmt',temp(1,1,k1n))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to mass adjustment:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
ccc      write (string,'(a12,i8)') 'convec, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
c     before = after
c     diapfl_time = 0.0
c     call diapfl(m,n,mm,nn,k1m,k1n)
c     call system_clock(after)
c     diapfl_time = real(after-before,4)/real(rate,4)
c     diapfl_total_time = diapfl_total_time + diapfl_time
c     total_time = total_time + diapfl_time
c     call sstbud(3,'   diapyc.mixing',temp(1,1,k1n))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to diapycn. mixing:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)*2./mixfrq
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
ccc      write (string,'(a12,i8)') 'diapfl, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      thermf_time = 0.0
      call thermf(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      thermf_time = real(after-before,4)/real(rate,4)
      thermf_total_time = thermf_total_time + thermf_time
      total_time = total_time + thermf_time
c
ccc      write (string,'(a12,i8)') 'thermf, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      enloan_time = 0.0
      call eice(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      enloan_time = real(after-before,4)/real(rate,4)
      enloan_total_time = enloan_total_time + enloan_time
      total_time = total_time + enloan_time
c
      before = after
      mxlayr_time = 0.0
c     call mxlayr(m,n,mm,nn,k1m,k1n)
      call mxkpp(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      mxlayr_time = real(after-before,4)/real(rate,4)
      mxlayr_total_time = mxlayr_total_time + mxlayr_time
      total_time = total_time + mxlayr_time
c     call sstbud(4,'  air-sea fluxes',tprime)
c     call sstbud(5,'     entrainment',temp(1,1,k1n))
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to thermal forcing:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
ccc      write (string,'(a12,i8)') 'mxlayr, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c
      before = after
      hybgen_time = 0.0
      call hybgen(m,n,mm,nn,k1m,k1n)
      call system_clock(after)
      hybgen_time = real(after-before,4)/real(rate,4)
      hybgen_total_time = hybgen_total_time + hybgen_time
      total_time = total_time + hybgen_time
c
ccc      write (string,'(a12,i8)') 'hybgrd, step',nstep
ccc      call comparall(m,n,mm,nn,string)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (diag_ape)
     . write (501,103) time,'  APE change due to vert.regridding:',
     .  hyc_pechg2(dp(1,1,k1n),th3d(1,1,k1n),31)
c     call sstbud(6,' vert.regridding',temp(1,1,k1n))
c
c ---------------------------------------------------------------------------
c
c --- output and diagnostic calculations
c
c ---------------------------------------------------------------------------
c
      if (diagno .or. mod(time+.0001,1.).lt..0002) then    !  once a day
c
c --- make line printer plot of mean layer thickness
        index=-1
        nflip=mod(nflip+1,2)
        do 705 k=1+(kk-1)*nflip,kk-(kk-1)*nflip,1-2*nflip
ccc     if (k.eq.kk-(kk-1)*nflip) index=1
        sum=0.
        totl=0.
        do 706 j=1,jj
        do 706 l=1,isp(j)
        do 706 i=ifp(j,l),ilp(j,l)
        totl=totl+oice(i,j)*scp2(i,j)
 706    sum=sum+dp(i,j,k+nn)*scp2(i,j)
        call linout(sum/(area*onem),charac(k),index)
 705    index=0
c --- add ice extend (%) to plot
        call linout(100.*totl/area,'I',1)
c
c --- diagnose mean sea surface height
        sum=0.
        sumice=0.
        do 704 j=1,jj
        do 704 l=1,isp(j)
        do 704 i=ifp(j,l),ilp(j,l)
c --- compute sea surface height
        util1(i,j)=montg(i,j,1)+thref*pbavg(i,j,m)
        sumice=sumice+oice(i,j)*scp2(i,j)
 704    sum=sum+util1(i,j)*scp2(i,j)
        write (lp,'(i9,'' mean sea srf.hgt. (mm):'',f9.2,f12.0)')nstep,
     .  sum/(area*thref*onemm),sumice*1.e-6
c
        call findmx(ip,oflxa2o,idm,ii,jj,'oflxa2o')
        call findmx(maska,asst,iia,iia,jja,'psst')
        call findmx(maska,sss,iia,iia,jja,'osss')
        call findmx(maska,uosurf,iia,iia,jja,'uosurf')
        call findmx(maska,vosurf,iia,iia,jja,'vosurf')
      end if                      ! once a day
c
      before = after
      call system_clock(after)
c
      before = after
      call system_clock(after)
      user_time =  cnuity_time+tsadvc_time+
     . momtum_time+hybgen_time+barotp_time+thermf_time+mxlayr_time+
     . real(after-before,4)/real(rate,4)
      total_time = total_time + real(after-before,4)/real(rate,4)
      write (lp,99009) user_time, '  sec for step ', nstep
      if (mod(nstep,5).eq.0) call flush(lp)
c
      if (.not.diagno) go to 23
c
c     For averaging over nstep steps, subtract times measured in the last
c     step from their corresponding totals and divide by (nstep-1).  (The
c     last step may include time used for diagnostics or other tasks.)
c     avg_time=(total_time-user_time)/(nstep-nstep0-1)
c     cnuity_avg_time=(cnuity_total_time-cnuity_time)/(nstep-nstep0-1)
c     tsadvc_avg_time=(tsadvc_total_time-tsadvc_time)/(nstep-nstep0-1)
c     momtum_avg_time=(momtum_total_time-momtum_time)/(nstep-nstep0-1)
c     barotp_avg_time=(barotp_total_time-barotp_time)/(nstep-nstep0-1)
c     thermf_avg_time=(thermf_total_time-thermf_time)/(nstep-nstep0-1)
c     enloan_avg_time=(enloan_total_time-enloan_time)/(nstep-nstep0-1)
c     mxlayr_avg_time=(mxlayr_total_time-mxlayr_time)/(nstep-nstep0-1)
c     hybgen_avg_time=(hybgen_total_time-hybgen_time)/(nstep-nstep0-1)
c     write (lp,99009) avg_time,' Avg Time at step ', nstep-1
c     write (lp,99009) total_time - user_time, ' Tot Time'
c     write (lp,99009) cnuity_avg_time,' cnuity Avg Time at ', nstep-1
c     write (lp,99009) tsadvc_avg_time,' tsadvc Avg Time at ', nstep-1
c     write (lp,99009) momtum_avg_time,' momtum Avg Time at ', nstep-1
c     write (lp,99009) barotp_avg_time,' barotp Avg Time at ', nstep-1
c     write (lp,99009) thermf_avg_time,' thermf Avg Time at ', nstep-1
c     write (lp,99009) enloan_avg_time,' enloan Avg Time at ', nstep-1
c     write (lp,99009) mxlayr_avg_time,' mxlayr Avg Time at ', nstep-1
c     write (lp,99009) hybgen_avg_time,' hybgen Avg Time at ', nstep-1
c
      write (lp,100) nstep,int((time+.001)/365.),mod(time+.001,365.)
 100  format (' ocn time step',i9,4x,'y e a r',i6,4x,'d a y',f9.2)
c
c --- output to history file
c
      call archiv(n,nn)
c --- diagnose meridional overturning and heat flux
c$OMP PARALLEL DO
        do 3 j=1,jj
        do 3 k=1,kk
        do 3 l=1,isp(j)
        do 3 i=ifp(j,l),ilp(j,l)
 3      p(i,j,k+1)=p(i,j,k)+dp(i,j,k+mm)
c$OMP END PARALLEL DO
      call overtn(mm)
c
      write (lp,105) nstep
 105  format (' step',i9,' -- archiving completed --')
c
c --- test for cyclic-in-j vs. noncyclic domain
c --- (closed-basin coditions are indicated by a land barrier at j=jj)
c
      jj1=jj-1
      do i=1,ii1
        if (ip(i,jj).gt.0) then
          jj1=jj
        end if
      end do
c
c --- output to line printer
c
ccc      call prtmsk(ip,util1,util3,idm,ii1,jj1,0.,1./(thref*onecm),
ccc     .     'sea surface height (cm)')
ccc      call prtmsk(ip,dpmixl(1,1),util3,idm,ii1,jj1,0.,1./onem,
ccc     .     'mixed layer depth (m)')
      call prtmsk(ip,temp(1,1,k1n),util3,idm,ii1,jj1,0.,10.,
     .     'mix.layer temp. (.1 deg)')
ccc      call prtmsk(ip,saln(1,1,k1n),util3,idm,ii1,jj1,35.,100.,
ccc     .     'mx.lay. salin. (.01 mil)')
      call prtmsk(ip,oice,util3,idm,ii1,jj1,0.,100.,
     .     'ice coverage (cm)')
        call prtmsk(maska,asst,util3,iia,iia,jja,0.,1.,
     .     'asst ')
      do 77 j=1,jj1
      do 77 l=1,isu(j)
      do 77 i=ifu(j,l),ilu(j,l)
 77   util1(i,j)=u(i,j,k1n)+ubavg(i,j,n)
ccc      call prtmsk(iu,util1,util3,idm,ii1,jj1,0.,1000.,
ccc     .     'u vel. (mm/s), layer 1')
      do 78 i=1,ii1
      do 78 l=1,jsv(i)
      do 78 j=jfv(i,l),jlv(i,l)
 78   util2(i,j)=v(i,j,k1n)+vbavg(i,j,n)
ccc      call prtmsk(iv,util2,util3,idm,ii1,jj1,0.,1000.,
ccc     .     'v vel. (mm/s), layer 1')
ccc      call prtmsk(iu,ubavg(1,1,n),util3,idm,ii1,jj1,0.,1000.,
ccc     .     'barotrop. u vel. (mm/s)')
ccc      call prtmsk(iv,vbavg(1,1,n),util3,idm,ii1,jj1,0.,1000.,
ccc     .     'barotrop. v vel. (mm/s)')
 23   continue
c --- accumulate fields for agcm
c$OMP PARALLEL DO
      do 201 j=1,jj
      do 201 l=1,isp(j)
      do 201 i=ifp(j,l),ilp(j,l)
      osst(i,j)=osst(i,j)+temp(i,j,k1n)*baclin/(3600.*real(nhr))
      osss(i,j)=osss(i,j)+saln(i,j,k1n)*baclin/(3600.*real(nhr))
      osiav(i,j)=osiav(i,j)+odmsi(i,j)*baclin*dtsrc/(3600.*real(nhr)) !kg/m2=>kg*.5*hr/m2
      omlhc(i,j)=spcifh*max(dp(i,j,k1n)/onem,thkmin)/thref  ! J/(m2*C)
      oogeoza(i,j)=(montg(i,j,1)+thref*pbavg(i,j,m))*g/(thref*onem) ! m^2/s^2
 201  continue
c$OMP END PARALLEL DO
c
      nsaveo=nsaveo+1
 15   continue
      print *,' chk ogcm saved over hr=',nsaveo*baclin/3600
      nsaveo=0
c
c       before = after
c       call system_clock(after)
c       user_time = real(after-before,4)/real(rate,4)
c       write (lp,99009) user_time, ' Time 2 just before STOP'
c       open (unit=no,file='run.status',status='unknown')
c       write (no,*) 'success'
c       close (unit=no)
c
      do 86 i=1,ii
      do 86 j=1,jj
      usf(i,j)=0.
      vsf(i,j)=0.
 86   continue
      do 87 j=1,jj
      do 87 l=1,isu(j)
      do 87 i=ifu(j,l),ilu(j,l)
 87   usf(i,j)=u(i,j,k1n)
      do 88 i=1,ii
      do 88 l=1,jsv(i)
      do 88 j=jfv(i,l),jlv(i,l)
 88   vsf(i,j)=v(i,j,k1n)
c
      call ssto2a(osst,asst)
      call ssto2a(osss,sss)
css   call iceo2a(omlhc,mlhc)
      call ssto2a(oogeoza,ogeoza)
      call ssto2a(osiav,utila)                 !kg/m*m per agcm time step
      call veco2a(usf,vsf,uosurf,vosurf)
      call findmx(ip,osiav,idm,ii,jj,'dmsi')
c$OMP PARALLEL DO PRIVATE(tf)
      do 204 ia=1,iia
      do 204 ja=1,jja
      if (focean(ia,ja).gt.0.) then
        gtemp(1,1,ia,ja)=asst(ia,ja)
        tf=tfrez(sss(ia,ja),0.)
        dmsi(1,ia,ja)=utila(ia,ja)                        !kg/m2 per agcm step
        dhsi(1,ia,ja)=utila(ia,ja)                        !J/m2 per agcm step
     .        *(tf*shi-lhm*(1-0.001*fsss*sss(ia,ja)))
        dssi(1,ia,ja)=1.d-3*dmsi(1,ia,ja)*sss(ia,ja)*fsss !kg/m2 per agcm step
c --- evenly distribute new ice over open water and sea ice
        if (aice(ia,ja).gt.1.e-3) then
          dhsi(2,ia,ja)=dhsi(1,ia,ja)
          dmsi(2,ia,ja)=dmsi(1,ia,ja)
          dssi(2,ia,ja)=dssi(1,ia,ja)
        endif
      endif
 204  continue
c$OMP END PARALLEL DO
      call system_clock(afogcm)
      write(*,'(a,2f10.1)') 'chk agcm/ogcm date at end of ogcm'
     .    ,(itime+1.)/nday,time
c
      return
      end
c
c> Revision history:
c>
c> May  1997 - removed statement "theta(1)=-thbase" after loop 14
c> June 1997 - added loop 60 to fix bug in vertical summation of -diaflx-
c> Mar. 2000 - conversion to SI units
c> June 2000 - added timing logic
c> Aug. 2000 - added code for changing number of OpenMP threads
c> Apr. 2001 - eliminated stmt_funcs.h
c> June 2001 - corrected sign error in diaflx
c> July 2001 - replaced archiving statements by 'call archiv'
c> Oct  2004 - map ice mass to agcm, and then calculate E
