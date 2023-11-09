#include "rundeck_opts.h"

      subroutine obio_ptend(vrbos,kmax,i,j,kdm,nstep,n_co2n,
     &                      DTS,mmo,ddxypo,n_abioDIC,num_tracers,SDIC) 

c  Computes tendencies of biological particles (units/hr) 
c  tracer
c  P(1) = nitrate (uM)
c  P(2) = ammonium (uM)
c  P(3) = silica (uM)
c  P(4) = iron (nM)
c  P(5) = diatoms (mg chl m-3)
c  P(6) = chlorophytes (mg chl m-3)
c  P(7) = cyanobacteria (mg chl m-3)
c  P(8) = coccolithophores (mg chl m-3)
c  P(9) = herbivores (mg chl m-3)
 
      USE obio_dim
      USE obio_incom,only: cnratio,cfratio,remin,obio_wss,bf,cchlratio
     .                    ,wsdeth,rkn,rks,rkf,Rm,phygross,bn,bs,solFe
     .                    ,mgchltouMC,uMtomgm3,NCrrat,O2thr
#ifdef exp_wsdiat
     .                    ,adiat_exp,bdiat_exp
#endif
#ifdef exp_wsdet
     .                    ,adet_exp,bdet_exp
#endif
#ifdef TRACERS_degC
     .                    ,tdegC  !@PL
#endif
      USE obio_forc, only: tirrq
      USE obio_com, only : dp1d,obio_P,obio_ws,P_tend,D_tend,C_tend
     .                    ,gro,rlamz,dratez1,dratez2,rmu3,rmu4 !@PL added rmu3,rmu4 in obio_com
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
     .                    ,O_tend,o21d
#endif
#ifdef TRACERS_abio_O2
     .                    ,Abo_tend
#endif
#endif                       

     .                    ,greff,pnoice,drate,tfac,regen,Fescav
     .                    ,wshc,rikd,rmuplsr,det
     .                    ,gcmax1d,covice_ij,atmFe_ij
     .                    ,temp1d,wsdet,tzoo,p1d
     .                    ,rhs,pp2_1d,flimit,obio_deltat,sday
#ifdef restoreIRON
!AR5 preprocessor option
     .                    ,Iron_BC
#endif
#ifdef TRACERS_degC
     .                    ,ndegC1d,Ndeg_tend
#endif

      implicit none

      integer,intent(in) :: kdm,nstep,n_co2n,n_abioDIC,num_tracers
      real, intent(in) :: dts,mmo,ddxypo,SDIC!(num_tracers) PLdbg

      integer i,j,k,kto
      integer nt,kmax

      real zoo(ntyp)      !herbivores (zooplankton)
      real dphy(ntyp)     !death rate of phytoplankton
      real viscfac(kdm), Sgronfix, SobioP1, ratio



      real ptot
      real Pzoo,Gzoo,Dzoo1,Dzoo2,exc
      real tirrqice,upn,upa,upf,ups

      real rnut2,tnit,tmp,rnut1,rmmn,framm,rmml,rmmlice,rmms
      real rmmf,rlim,rlimice,grate,rlimnfix,rlimrkn,rfix,gratenfix1
      real graterkn,gratenlim,gratenfix,gron,gronfix(kdm)
      real dphyt

      logical vrbos

      real term
      character*7, dimension(12) ::  bio_var = (/
     .  ' NO3   ',' NH4   ',' SiO2  ',' Iron  ','Diatom ','Chlphy ',
     .  'Cyanob ','Coccol ','Herbiv ','N/Cdet ','Silica ','Fe_det '/)

       rhs=0.0
       obio_ws = 0.0
       P_tend = 0.0
       C_tend = 0.0
       D_tend = 0.0
#ifdef TRACERS_Ocean_O2
#ifdef TRACERS_bio_O2
       O_tend = 0.0
#endif
#ifdef TRACERS_abio_O2
       Abo_tend = 0.0
#endif
#endif
#ifdef TRACERS_degC
       Ndeg_tend = 0.0
#endif
       wsdet = 0.0
       rmu4 = 0.0
       rmu3 = 0.0
       zoo  = 0.0
       dphy = 0.0
       gro = 0.0
       viscfac = 0.0
       pp2_1d = 0.0
       gronfix = 0.0


       !define no ice points based on covice (here: covice_ij)
      pnoice(1)=1.-covice_ij
      do k=2,kdm
         pnoice(k)=pnoice(1)
      enddo

      bs = 2.0*bn

!tendency terms are computed in the mid-level (m)

c  Start Model Space Loop
!      m = indext2     !index of "past" (t-1)

! River runoff applied
#ifdef OBIO_RUNOFF
      call obio_rivers(vrbos)
#endif

!Iron + atm iron: disperse in layer and convert to nM
!we do not need to multiply this by pnoice, because iron
!is deposited over ice and presumably when ice melts will enter
!the ocean
       
        k = 1
       term = atmFe_ij*solFe*1.d-3/max(p1d(k+1),1.e-3)
       rhs(k,4,4) = term
       P_tend(k,4) = P_tend(k,4) + term

!change: March 15, 2010
!      do k=2,kmax
!         rhs(k,4,4) = 0.0
!         P_tend(k,4) = 0.0
!      enddo

!Grazing/regeneration of ammonium
       do k = 1,kmax

         ptot = 0.0
         do nt = nnut+1,ntyp-nzoo
          ptot = ptot + obio_P(k,nt)
         enddo
         ptot = max(ptot,1.0E-36)

         Pzoo = obio_P(k,ntyp)
         gzoo = tzoo*Rm*(1.0-exp(-rlamz*ptot))*Pzoo
         dzoo1 = dratez1*Pzoo
         dzoo2 = dratez2*Pzoo*Pzoo

!!!!     term = ((1.0-greff)*gzoo-dzoo1-dzoo2) * pnoice(k)
         term = gzoo * pnoice(k)           !herbivore growth
         rhs(k,9,9) = term
         P_tend(k,ntyp) = P_tend(k,ntyp) + term

         term = -greff*gzoo * pnoice(k)    !grazing efficiency
         rhs(k,9,10) = term
         P_tend(k,ntyp) = P_tend(k,ntyp) + term

         term = -dzoo1 * pnoice(k)        !death of herbiv
         rhs(k,9,11) = term
         P_tend(k,ntyp) = P_tend(k,ntyp) + term

         term = -dzoo2 * pnoice(k)        !death of herbiv
         rhs(k,9,12) = term
         P_tend(k,ntyp) = P_tend(k,ntyp) + term

         do nt = nnut+1,ntyp-nzoo
          !fraction of grazing for this group
          zoo(nt) = gzoo*obio_P(k,nt)/ptot
          dphy(nt) = drate*obio_P(k,nt)

          term = -zoo(nt) * pnoice(k)
          rhs(k,nt,9) = term
          P_tend(k,nt) = P_tend(k,nt) + term

          term = -dphy(nt) * pnoice(k)           !death of phytoplankton
          rhs(k,nt,nt) = term
          P_tend(k,nt) = P_tend(k,nt) + term

         enddo

         !remineralization 
         !@PL ratio of o2/nitrate remineralization, not applied to Silica
         !@PL define delta function for O2
          NCrrat = 1.d0
#ifdef TRACERS_bio_O2
         if (o21d(k).le.O2thr) then
          NCrrat = 0.4d0
         else
          NCrrat = 1.d0
         endif
#endif
         term = NCrrat*tfac(k)*remin(1)*det(k,1)/cnratio * pnoice(k)
         rhs(k,1,10) = term
         P_tend(k,1) = P_tend(k,1) + term

         !regeneration from zooplankton
         exc = greff*gzoo
!!!!!    term = bn*(exc + regen*dzoo2) * pnoice(k)
         term = bn * exc * pnoice(k)
         rhs(k,2,9) = term
         P_tend(k,2) = P_tend(k,2) + term

         term = bn*regen*dzoo2 * pnoice(k)
         rhs(k,2,10) = term
         P_tend(k,2) = P_tend(k,2) + term

         !remineralization 
         term = tfac(k)*remin(2)*det(k,2) * pnoice(k)
         rhs(k,3,11) = term    !put this in diff column
         P_tend(k,3) = P_tend(k,3) + term

         !regeneration from zooplankton
         term = bf*exc * pnoice(k)
         rhs(k,4,9) = term
         P_tend(k,4) = P_tend(k,4) + term

         term = bf*regen*dzoo2 * pnoice(k)
         rhs(k,4,11) = term
         P_tend(k,4) = P_tend(k,4) + term

         !remineralization 
         term = NCrrat*tfac(k)*remin(3)*det(k,3) * pnoice(k)
         rhs(k,4,12) = term                        !put this in diff column
         P_tend(k,4) = P_tend(k,4) + term

         term = -Fescav(k)
         rhs(k,4,13) = term             !this should actually be rhs(4,4) but we
                                        !have already defined it so let it be rhs(4,13)
         P_tend(k,4) = P_tend(k,4) + term

!!!!     term = exc*mgchltouMC*pnoice(k)
!!!! *                  + regen*dzoo2*mgchltouMC*pnoice(k)
         term = exc*mgchltouMC*pnoice(k)
         rhs(k,13,9) = term
         C_tend(k,1) = C_tend(k,1) + term

         term = regen*dzoo2*mgchltouMC*pnoice(k)
         rhs(k,13,10) = term
         C_tend(k,1) = C_tend(k,1) + term

         dphyt = 0.0
         do nt = nnut+1,ntyp-nzoo
          dphyt = dphyt + dphy(nt)
         enddo
 
!1st detrital fraction is carbon
         term = dphyt*cchlratio * pnoice(k)
         rhs(k,10,5) = term
         D_tend(k,1) = D_tend(k,1) + term

!!!!     term = dzoo1*cchlratio * pnoice(k)
!!!! .        + (1.0-regen)*dzoo2*cchlratio * pnoice(k)
         term = dzoo1*cchlratio * pnoice(k)       !death of herbiv
         rhs(k,10,7) = term
         D_tend(k,1) = D_tend(k,1) + term

         term = dzoo2*cchlratio * pnoice(k)
         rhs(k,10,8) = term
         D_tend(k,1) = D_tend(k,1) + term

         term = -regen*dzoo2*cchlratio * pnoice(k)
         rhs(k,10,9) = term
         D_tend(k,1) = D_tend(k,1) + term

         term = -NCrrat*tfac(k)*remin(1)*det(k,1) * pnoice(k)
         rhs(k,10,10) = term
         D_tend(k,1) = D_tend(k,1) + term

!@auth PL: remove a portion of degradable carbon to become non-degradable
!@ added Jan 19, 2020.
#ifdef TRACERS_degC
         term = -tdegC*det(k,1) * pnoice(k)
         rhs(k,10,11) = term
         D_tend(k,1) = D_tend(k,1) + term

         term = tdegC*det(k,1)*pnoice(k)
         rhs(k,ndimndegC,11) = term
         Ndeg_tend(k) = Ndeg_tend(k) + term
#endif

!2nd detrital fraction is silica
         term = bs*dphy(nnut+1) * pnoice(k)
         rhs(k,11,5) = term
         D_tend(k,2) = D_tend(k,2) + term

         term = bs*zoo(nnut+1) * pnoice(k)
         rhs(k,11,9) = term
         D_tend(k,2) = D_tend(k,2) + term

         term = -tfac(k)*remin(2)*det(k,2) * pnoice(k)
         rhs(k,11,11) = term
         D_tend(k,2) = D_tend(k,2) + term

!3rd detrital fraction is iron
         term = bf*dphyt * pnoice(k)
         rhs(k,12,5) = term
         D_tend(k,3) = D_tend(k,3) + term

         term = bf*dzoo1 * pnoice(k)
         rhs(k,12,9) = term
         D_tend(k,3) = D_tend(k,3) + term

         term = bf*dzoo2 * pnoice(k)
         rhs(k,12,6) = term
         D_tend(k,3) = D_tend(k,3) + term

         term = -bf*regen*dzoo2 * pnoice(k)
         rhs(k,12,7) = term
         D_tend(k,3) = D_tend(k,3) + term

         term = -NCrrat*tfac(k)*remin(3)*det(k,3) * pnoice(k)
         rhs(k,12,12) = term
         D_tend(k,3) = D_tend(k,3) + term

!change June 1, 2010
         term = 0.1 * Fescav(k)
#ifdef DETSCAV     
!this is an AR5 preprocessor option
         term = Fescav(k)
#endif
!endofchange
         rhs(k,12,4) = term
         D_tend(k,3) = D_tend(k,3) + term

cdiag    if (vrbos) then
cdiag     if(k.eq.1) write(103,'(a,a)')
cdiag.    'nstep    k  dp1d(k)   dphyt   dzoo1     dzoo2   dphy(5) '
cdiag.   ,'zoo(5)   tfac(k)  D_t(1)     D_t(2) D_t(3)'
cdiag    write(103,'(2i5,10(1x,es8.2))')nstep,k,dp1d(k),
cdiag.   dphyt,dzoo1,dzoo2,dphy(nnut+1),zoo(nnut+1),tfac(k),
cdiag.   D_tend(k,1),D_tend(k,2),D_tend(k,3)
cdiag    endif
cdiag    if (vrbos) write(*,*)'ptend5 :',
cdiag.   nstep,i,j,k,
cdiag.   dp1d(k),dphyt,dzoo1,dzoo2,dphy(nnut+1),zoo(nnut+1),tfac(k),
cdiag.   D_tend(k,1),D_tend(k,2),D_tend(k,3)

       enddo !kmax

c Day: Grow
      flimit(:,:,:)=0.0
      do k = 1,kmax

      if (tirrq(k) .gt. 0.0)then
        tirrqice = tirrq(k)*0.01  !reduce light in ice by half

c Light-regulated growth
!!#if NCHL_DEFINED > 0
      if (nchl > 0) then

        ! Diatoms
        nt = 1

        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rmml =0.d0; rmmlice=0.d0; rmmn=0.d0; rmms=0.d0; rmmf=0.d0;
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
         tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))     !nitrate
          tmp = 1.0 - rnut2
        ! Enforce preferential utilization of ammonium
          rnut1 = min(tmp,tnit)
           rmmn = rnut1 + rnut2         !nitrate limitation
          if ( rmmn .ne. 0.d0 ) then
            framm = rnut2/rmmn
          else
            framm = 0.d0
          endif
           rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))         !light limitation
        rmmlice = tirrqice/(tirrqice+0.5*rikd(k,nt))

        rmms = obio_P(k,3)/(rks(nt)+obio_P(k,3))      !silicate limitation
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron limitation
        rlim = min(rmml,rmmn,rmms,rmmf)
        rlimice = min(rmmlice,rmmn,rmms,rmmf)

        !compute limiting factor array
        flimit(k,nt,1) = rmml*pnoice(k)
        flimit(k,nt,2) = rmmlice*(1.0-pnoice(k))
        flimit(k,nt,3) = rmmn*pnoice(k)  
        flimit(k,nt,4) = rmms*pnoice(k)  
        flimit(k,nt,5) = rmmf*pnoice(k)  

        grate = rmuplsr(k,nt)*rlim*pnoice(k)
     .        + rmuplsr(k,nt)*rlimice*(1.0-pnoice(k))
        rmu4(k,nt) = grate*framm
        rmu3(k,nt) = grate*(1.0-framm)
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt)
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

        !Net primary production  in mgC/m2/day because:
        ! [gro]= mg,chl/m3/s, [dp]= m,                     !July 2016
        ! [cchlratio]= mgl/mgl,[phygross]=no units
        pp2_1d(k,nt) = gro(k,nt) * phygross 
     .               * dp1d(k) * cchlratio *sday    !July 2016

      endif

!!#endif


!!#if NCHL_DEFINED > 1
      if (nchl > 1) then

! Chlorophytes
        nt = 2

        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rmml =0.d0; rmmlice=0.d0; rmmn=0.d0; rmms=0.d0; rmmf=0.d0;
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit  = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp   = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2         !nitrate limitation                
        if ( rmmn .ne. 0.d0 ) then
          framm = rnut2/rmmn                           
        else
          framm = 0.d0
        endif
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))     !light limitation
        rmmlice = tirrqice/(tirrqice+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron limitation
        rlim = min(rmml,rmmn,rmmf)
        rlimice = min(rmmlice,rmmn,rmmf)

        !compute limiting factor array
        flimit(k,nt,1) = rmml
        flimit(k,nt,2) = rmmlice*(1.0-pnoice(k))
        flimit(k,nt,3) = rmmn   
        flimit(k,nt,4) = rmms   
        flimit(k,nt,5) = rmmf

        grate = rmuplsr(k,nt)*rlim * pnoice(k)
     .        + rmuplsr(k,nt)*rlimice * (1.0-pnoice(k))
        rmu4(k,nt) = grate*framm
        rmu3(k,nt) = grate*(1.0-framm)
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt)
        rhs(k,nt+nnut,13) = term  
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

        !Net primary production  in mgC/m2/day because:
        ! [gro]= mg,chl/m3/s, [dp]= m,
        ! [cchlratio]= mgl/mgl,[phygross]=no units
        pp2_1d(k,nt) = gro(k,nt) * phygross
     .               * dp1d(k) * cchlratio * sday    !July 2016
      endif
!!#endif


!!#if NCHL_DEFINED > 2
      if (nchl > 2) then
! Cyanobacteria
        nt = 3
        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rmml =0.d0; rmmlice=0.d0; rmmn=0.d0; rmms=0.d0; rmmf=0.d0;
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2
        if ( rmmn .ne. 0.d0 ) then
          framm = rnut2/rmmn
        else
          framm = 0.d0
        endif
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
        rlim = min(rmml,rmmn,rmmf)
        rlimnfix = min(rmml,rmmf)         !limitation for N2 fixation
        rlimrkn = min(rmml,rkn(nt),rmmf)   !limitation at kn

        !compute limiting factor array
        flimit(k,nt,1) = rmml
        flimit(k,nt,2) = rmmlice*(1.0-pnoice(k))
        flimit(k,nt,3) = rmmn   
        flimit(k,nt,4) = rmms   
        flimit(k,nt,5) = rmmf

        grate = rmuplsr(k,nt)*rlim*pnoice(k)
        rmu4(k,nt) = grate*framm
        rmu3(k,nt) = grate*(1.0-framm)
        rfix = 0.25*exp(-(75.0*obio_P(k,nt+nnut)))
        rfix = max(rfix,0.0)
c        rfix = min(rfix,0.2)

        gratenfix1 = rmuplsr(k,nt)*rlimnfix*rfix  !N fix
        graterkn = rmuplsr(k,nt)*rlimrkn
        gratenlim = graterkn - grate
        gratenfix = min(gratenlim,gratenfix1)  !N fix cannot exceed kn
        gratenfix = max(gratenfix,0.0) * pnoice(k)

        gron = grate*obio_P(k,nt+nnut)
        gronfix(k) = gratenfix*obio_P(k,nt+nnut)
!       gro(k,nt) = gron + gronfix(k)
        gro(k,nt) = gron             !add the gronfix later

        term = gro(k,nt)
        rhs(k,nt+nnut,13) = gron
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

        !Net primary production  in mgC/m2/day because:
        ! [gro]= mg,chl/m3/s, [dp]= m,
        ! [cchlratio]= mgl/mgl,[phygross]=no units
        pp2_1d(k,nt) = gro(k,nt) * phygross
     .               * dp1d(k) * cchlratio * sday     !July 2016

      endif
!!#endif


!!#if NCHL_DEFINED > 3
      if (nchl > 3) then
! Coccolithophores
        nt = 4

        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rmml =0.d0; rmmlice=0.d0; rmmn=0.d0; rmms=0.d0; rmmf=0.d0;
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2
        if ( rmmn .ne. 0.d0 ) then
          framm = rnut2/rmmn
        else
          framm = 0.d0
        endif
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
        rlim = min(rmml,rmmn,rmmf)

        !compute limiting factor array
        flimit(k,nt,1) = rmml
        flimit(k,nt,2) = rmmlice*(1.0-pnoice(k))
        flimit(k,nt,3) = rmmn   
        flimit(k,nt,4) = rmms   
        flimit(k,nt,5) = rmmf

        grate = rmuplsr(k,nt)*rlim*pnoice(k)
        rmu4(k,nt) = grate*framm
        rmu3(k,nt) = grate*(1.0-framm)
        gro(k,nt) = grate*obio_P(k,nt+nnut)

!!!     term = gro(k,nt) * pnoice(k)
        term = gro(k,nt)
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term
        gcmax1d(k) = max(gcmax1d(k),grate)

        !Net primary production  in mgC/m2/day because:
        ! [gro]= mg,chl/m3/s, [dp]= m,
        ! [cchlratio]= mgl/mgl,[phygross]=no units
        pp2_1d(k,nt) = gro(k,nt) * phygross
     .               * dp1d(k) * cchlratio * sday    !July 2016

      endif
!!#endif


!!#if NCHL_DEFINED > 4
      if (nchl > 4) then
! Dinoflagellates
        nt = 5

        ! Nutrient-regulated growth; Michaelis-Menton uptake kinetics
        rmml =0.d0; rmmlice=0.d0; rmmn=0.d0; rmms=0.d0; rmmf=0.d0;
        rnut2 = obio_P(k,2)/(rkn(nt)+obio_P(k,2))     !ammonium
        tnit = obio_P(k,1)/(rkn(nt)+obio_P(k,1))      !nitrate
        tmp = 1.0 - rnut2

        ! Enforce preferential utilization of ammonium
        rnut1 = min(tmp,tnit)
        rmmn = rnut1 + rnut2
        if ( rmmn .ne. 0.d0 ) then
          framm = rnut2/rmmn
        else
          framm = 0.d0
        endif
        rmml = tirrq(k)/(tirrq(k)+0.5*rikd(k,nt))
        rmmf = obio_P(k,4)/(rkf(nt)+obio_P(k,4))      !iron
        rlim = min(rmml,rmmn,rmmf)

        !compute limiting factor array
        flimit(k,nt,1) = rmml
        flimit(k,nt,2) = rmmlice*(1.0-pnoice(k))
        flimit(k,nt,3) = rmmn   
        flimit(k,nt,4) = rmms   
        flimit(k,nt,5) = rmmf

        grate = rmuplsr(k,nt)*rlim
        rmu4(k,nt) = grate*framm * pnoice(k)
        rmu3(k,nt) = grate*(1.0-framm) * pnoice(k)
        gro(k,nt) = grate*obio_P(k,nt+nnut)

        term = gro(k,nt)
        rhs(k,nt+nnut,13) = term
        P_tend(k,nt+nnut) = P_tend(k,nt+nnut) + term

        !Net primary production  in mgC/m2/day because:
        ! [gro]= mg,chl/m3/s, [dp]= m,
        ! [cchlratio]= mgl/mgl,[phygross]=no units
        pp2_1d(k,nt) = gro(k,nt) * phygross
     .               * dp1d(k) * cchlratio * sday    !July 2016

      endif
!!#endif
       endif !tirrq(k) .gt. 0.0
      enddo  !kmax

#ifdef new_NFIXATION
!this is an AR5 preprocessor option
      Sgronfix = 0.d0
      SobioP1  = 0.d0
      do k=1,kmax
        Sgronfix = Sgronfix + gronfix(k)  * dp1d(k)  !integrate gronfix
         SobioP1 =  SobioP1 + obio_P(k,1) * dp1d(k)  !integrate nitrogen
      enddo

      !New fixation takes place only if the Sgronfix does not exceed 50% of SobioP1
      if (bn*Sgronfix*obio_deltat .le. 0.5*SobioP1) then

      !cyanobacteria growth due to N-fixation
      do k=1,kmax
      if (tirrq(k) .gt. 0.0)then
        term = gronfix(k)
        rhs(k,7,12) = term
        P_tend(k,7) = P_tend(k,7) + term    !cyanobacteria only     !Ctest1
        gro(k,3) = gro(k,3) + gronfix(k)
        pp2_1d(k,3) = pp2_1d(k,3)
     .              +gronfix(k)*phygross*dp1d(k)*cchlratio * sday   !July 2016
      endif
      enddo

      !nitrate reduction due to uptake via N-fixation
      do k=1,kmax
       !we distribute total_gronfix into each layer according to its layer thickness
       kto = (kmax-k)+1
       if (dp1d(kto).gt.1.) then    !avoid massless layers
       ratio = obio_P(kto,1)*dp1d(kto) / SobioP1
       term = -bn* ratio * Sgronfix/dp1d(kto)
       rhs(kto,1,11) = term
       P_tend(kto,1) = P_tend(kto,1) + term     !Ctest1

c      if (vrbos)write(*,'(a,5i6,8e18.8)')'Nfixation_new diag:',
c    .     nstep,i,j,k,kto,dp1d(k),dp1d(kto),
c    .     obio_P(kto,1),SobioP1,ratio,Sgronfix,
c    .     rhs(kto,1,11),rhs(k,7,12)*bn
       endif
      enddo
      endif  !Sgronfix*dt <= 0.5*SobioP1
#endif

! Nutrient uptake
      do k=1,kmax
      if (tirrq(k) .gt. 0.0)then
        upn = 0.0
        upa = 0.0
        upf = 0.0
        do nt = 1,nchl
         term = -bn*(rmu3(k,nt)*obio_P(k,nnut+nt))
         upn = upn + term
         rhs(k,1,nnut+nt) = rhs(k,1,nnut+nt) + term

         term = -bn*(rmu4(k,nt)*obio_P(k,nnut+nt))
         upa = upa + term
         rhs(k,2,nnut+nt) = rhs(k,2,nnut+nt) + term

         term = -bf*gro(k,nt)
         upf = upf + term
         rhs(k,4,nnut+nt) = rhs(k,4,nnut+nt) + term
        enddo

        term = -bs*gro(k,1)
        ups = term
        rhs(k,3,nnut+1) = term

        P_tend(k,1) = P_tend(k,1) + upn
        P_tend(k,2) = P_tend(k,2) + upa
        P_tend(k,3) = P_tend(k,3) + ups
        P_tend(k,4) = P_tend(k,4) + upf

        endif !tirrq(k) .gt. 0.0
      enddo  !kmax

!#ifndef new_NFIXATION
!!  Fix up nitrogen uptake from N-fixation -- give back to water
!!  column in reverse order of layers
!      do k = 1,kmax
!       kto = (kmax-k)+1
!       term = - bn*gronfix(k)
!!    .      * dp1d(k)/max(1.,dp1d(kto))  !protect against zero, in case dp1d(kto) is vanishing
!       rhs(kto,1,11) = term
!       P_tend(kto,1) = P_tend(kto,1) + term
!       if (vrbos)write(*,'(a,5i6,5e18.8)')'Nfixation diag:',
!     .     nstep,i,j,k,kto,dp1d(k),dp1d(kto),gronfix(k),
!     .     rhs(kto,1,11),rhs(k,7,12)
!      enddo
!#endif


#ifdef restoreIRON
!this is an AR5 preprocessor option
!iron bottom sink/source
!whether sink or source is determined from Iron_BC
!Iron_BC > 0 sink of iron through sedimentation
!Iron_BC < 0 source of iron through sediment resuspension (Moore et al, 2004) 
      k = kmax
        if (p1d(kmax) >= 3700.) then        !for deep regions, bottom cell (lower 200m)
            term =  - obio_P(k,4) / (200/(Iron_BC*3700))
            term = term/(365.d0*sday)     !convert to per s   July 2016
            rhs(k,4,14) = term
            P_tend(k,4) = P_tend(k,4) + term
        endif
!alternative
!       if (p1d(kmax) <= 1100.) then    !for shelf regions, bottom cell
!         !!do k = 1,kmax
!           k=kmax
!           term =  50.d0/(24.d0*dp1d(k))     ! in nano-moleFe/m3/hr
!           rhs(k,4,14) = term
!           P_tend(k,4) = P_tend(k,4) + term
!         !!enddo
!       write(*,'(a,5i5,3e12.4)')'obio_ptend, iron source:',
!    .  nstep,i,j,k,kmax,dp1d(kmax),P_tend(k,4)-term,term
!       endif
#endif


#ifdef TRACERS_Alkalinity
#ifdef TOPAZ_params
      call obio_alkalinity_topaz(vrbos,kmax,i,j)
#else
      call obio_alkalinity(kmax,i,j,nstep)
#endif
#endif

!NOT FOR HYCOM:dts,mo,dxypo, n_abioDIC, num_tracers,trmo and oij passed
!to the subroutine.  
      call obio_carbon(gro,vrbos,kmax,i,j,nstep,kdm,n_co2n,
     &                 DTS,mmo,ddxypo,n_abioDIC,num_tracers,SDIC)
#ifdef TRACERS_Ocean_O2
!@PL sum tendencies of O2 in obio_O2
      call obio_o2(gro,vrbos,kmax,i,j,nstep,kdm,
     &                       DTS,mmo,ddxypo)
#endif

 107  format (/'lyr',i3,4x,'amount   tndcy   ',
     .    9(2x,a7)/(a7,2es9.1,2x,6es9.1))

!compute here rates but all detritus update done inside the update routine

c Sinking rate temperature (viscosity) dependence (also convert to /hr) -> convert to /s July 2016
      do k = 1,kmax
        viscfac(k) = 0.451 + 0.0178*temp1d(k)
      enddo
      do nt = 1,nchl
       do k = 1,kmax
!        obio_ws(k,nt) = obio_wsh(nt)*viscfac(k)*pnoice(k)
         obio_ws(k,nt) = obio_wss(nt)*viscfac(k)*pnoice(k)    !July 2016
       enddo
      enddo
#ifdef exp_wsdiat
! exponential profile coefficients for diatoms
      nt = 1  !diatoms ONLY
      do k = 1,kmax
        obio_ws(k,nt) = viscfac(k)*pnoice(k)
     .                * adiat_exp *exp(obio_P(k,nnut+nt)*bdiat_exp)   !m/s
     .                / sday
      obio_ws(k,nt) = min(obio_ws(k,nt), 20./3600./24.)    !max ws for diat 20m/day
      enddo
#endif
      nt = 4
      do k = 1,kmax
        obio_ws(k,nt) = wshc(k)*viscfac(k)*pnoice(k)
      enddo
      do nt = 1,nchl
        obio_ws(kmax+1,nt)=obio_ws(kmax,nt)
      enddo
      do nt = 1,ndet
       do k = 1,kmax
         wsdet(k,nt) = wsdeth(nt)*viscfac(k)*pnoice(k)
#ifdef exp_wsdet
! exponential profile coefficients for detritus
         wsdet(k,nt) = viscfac(k)*pnoice(k)
     .               * adet_exp(nt)* exp(det(k,nt)*bdet_exp(nt))   !m/s 
     .               / sday
      wsdet(k,nt) = min(wsdet(k,nt), 50./3600./24.)   !max wsdet 50 m/day
#endif
       enddo
      enddo

!     if(vrbos)then
!       do k=1,kmax
!          write(*,'(a,4i5,8e12.4)')
!    .       'obio_ptend, ws:',
!    .        nstep,i,j,k,temp1d(k),viscfac(k),pnoice(k)
!    .       ,obio_P(k,nnut+nt),(obio_ws(k,nt),nt=1,nchl)
!       enddo
!       do k=1,kmax
!          write(*,'(a,4i5,9e12.4)')
!    .       'obio_ptend, wsdet:',
!    .        nstep,i,j,k,temp1d(k),viscfac(k),pnoice(k)
!    .       ,(det(k,nt),nt=1,ndet),(wsdet(k,nt),nt=1,ndet)
!       enddo
!     endif  !vrbos

c Save method for hard boundary condition (no flux)
c      srate = 0.0 - obio_wsh(n)*tracer(i,k-1,m,n)

      return
      end subroutine obio_ptend
