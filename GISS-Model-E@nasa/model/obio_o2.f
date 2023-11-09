#include "rundeck_opts.h"
!@sum tendencey terms for oxygen. 
!@auth: Paul Lerner, A. Romanou, D.R. Nicholson
!NOT FOR HYCOM:dts,mo,dxypo, n_abioDIC, num_tracers,trmo and oij passed to the subroutine.
      subroutine obio_o2(gro,vrbos,kmax,i,j,nstep,kdm,
     &                       DTS,mmo,ddxypo)

      USE MODEL_COM, only: dtsrc
      USE CONSTANT, only: tf
      USE obio_dim
      USE obio_incom, only : cnratio,rlamdoc,rkdoc1,rkdoc2
     .                      ,rlampoc,uMtomgm3,Pzo,stdslp
     .                      ,excz,resz,remin,excp,resp,bn,cchlratio
     .                      ,mgchltouMC,bf,ko2,HvO2,O2thr
     .                      ,ro2c_DET,ro2c_NH4,ro2c_NO3,NCrrat
      USE obio_forc, only: wind,tirrq
      USE obio_com, only : obio_P,P_tend
     .                    ,tfac,det,D_tend,tzoo,pnoice,pHsfc
     .                    ,temp1d,saln1d,dp1d,rhs,alk1d,trmo_unit_factor
     .                    ,rho_water,docbac,dicresp
     .                    ,rmu3,rmu4,rho1d,p1d,gronfix
#ifdef TRACERS_bio_O2
     .                    ,O_tend,o21d,pO2_ij,o2flux
#endif
#ifdef TRACERS_abio_O2
     .                    ,Abo_tend,abo21d,pabO2_ij,abo2flux
#endif

!    .                    ,dic_river_sink,p1d


      USE OFLUXES, only: oAPRESS
#ifdef pH2O_Ocean

     .                  ,ocnatm
#endif

      use TimeConstants_mod, only: SECONDS_PER_HOUR, DAYS_PER_YEAR,
     &                             HOURS_PER_DAY
      

      use dictionary_mod, only: get_param
      use domain_decomp_1d, only: am_i_root

      implicit none

      integer,intent(in) :: kdm,nstep
      real, intent(in) :: dts,mmo,ddxypo

!     real, parameter :: awan=0.337d0/(3.6d5) !piston vel coeff., from
!                                             !Wanninkof 1992, but adjusted
!                                             !by OCMIP, and converted from
!                                             !cm/hr to m/s
      real, parameter :: awan=0.251d0/(3.6d5) !piston vel coeff., OCMIP2016
                                              !converted from
                                              !cm/hr to m/s
      integer :: i,j,k

      integer :: nt,kmax
      real  :: rmmzoo,docexcz,rndep,docdep
      real  :: docdet,o2resz,sumdoc,sumutk,sumres,totgro
      real  :: docexcp(nchl),sco2,sco2arg,wssq,rkwo2
      real  :: Ts,tk,tk100,tk1002,O2sat0,xo2,flxmolm3,flxmolm3h,
     &         flxmolkg,flxmolkgh

#ifdef TRACERS_bio_O2
      real  :: delto2
#endif
#ifdef TRACERS_abio_O2
      real  :: deltabo2
#endif
      real  :: ta,ta2,ta3,ta4,ta5,ps2,ps,Ksol,SLP,pH2O,O2sat
      real  :: gro(kdm,nchl)
      real  :: termb4(kdm,nchl),termb5(kdm,nchl)
      real  :: term,termab,dummy
      real  :: termb1(kdm),termb2(kdm),termb3(kdm),termb6(kdm)
      real bs


      logical vrbos


#ifdef TRACERS_bio_O2
      do k = 1,kmax
!---------------------------------------------------------------

!@PL define delta function for O2

        if (o21d(k).le.O2thr) then
         HvO2 = 0.d0
        else
         HvO2 = 1.d0
        endif

         
!O2 RRR  O2 tracer index == 15??
        o2resz =HvO2*ro2c_DET*tzoo*resz*obio_P(k,ntyp) !zoopl O2 consump (resp)
        term = o2resz*mgchltouMC * pnoice(k)/rho1d(k) !@PL mg/(m3 s) -> mmol/(kg s)
        rhs(k,ndimo2,15) = term
        termb1(k)=term
        O_tend(k) = O_tend(k) + term

! respiration of DOC
!@PL the o2-dependant factor is used to prevent negative O2 values,
! and represents transfer from aerobic remin to denitr below 20mmol/kg O2
        term = HvO2* !
     .          ro2c_DET*docbac(k) * pnoice(k)/rho1d(k) !@PL mg/(m3 s) -> mmol/(kg s)
        rhs(k,ndimo2,14) = term
        termb2(k)=term
        O_tend(k) = O_tend(k) + term
! remineralization of detritus
        term = HvO2*
     .         ro2c_DET*tfac(k)*remin(1)*det(k,1)/uMtomgm3 * pnoice(k)
     .         /rho1d(k) !@PL uM/s -> mmol/(kg s)
        rhs(k,ndimo2,10) = term
        termb3(k)=term
        O_tend(k) = O_tend(k) + term
      enddo  !k=1,kmax
      ! Phytoplankton components related to growth
      do k = 1,kmax
      if (tirrq(k) .gt. 0.d0)then
      do nt = 1,nchl
        ! totgro = gro(k,nt)
        ! RRR: adding O2 production here - need to check units
        ! b/c of different O2:C ratios, need to accumulate for each PFT
        ! still need to add O2 prod during n-fix growth

        ! ammonium supported O2 production
        term = -ro2c_NH4*rmu4(k,nt)*obio_P(k,nt+nnut)
     .    * mgchltouMC/rho1d(k) !@PL mg/(m3 s) -> mmol/(kg s)
        rhs(k,ndimo2,6) = rhs(k,ndimo2,6) + term    !accumulate
        termb4(k,nt)=term
        O_tend(k) = O_tend(k) + term
        ! NO3 supported O2 production
        term =  -ro2c_NO3*rmu3(k,nt) *obio_P(k,nt+nnut)
     .   * mgchltouMC/rho1d(k) !@PL mg/(m3 s) -> mmol/(kg s)
        rhs(k,ndimo2,7) = rhs(k,ndimo2,7) + term    !accumulate
        termb5(k,nt)=term
        O_tend(k) = O_tend(k) + term
#ifdef new_NFIXATION

        ! for nitrogen fixation, add O2 only for cyanobacteria
      if (nt.eq.3) then
        term = -ro2c_NO3*gronfix(k)
     .           *mgchltouMC/rho1d(k) !@PL mg/m3 s) -> mmol/(kg s)
        rhs(k,ndimo2,8) = rhs(k,ndimo2,8) + term !accumulate
        O_tend(k) = O_tend(k) + term
      endif
#endif

        ! autotrophic resp O2 consumption
        term =  ro2c_DET * gro(k,nt)*resp
     .    * mgchltouMC/rho1d(k) !@PL mg/(m3 s) -> mmol/(kg s)
        rhs(k,ndimo2,5) = rhs(k,ndimo2,5) + term    !accumulate
       termb6(k)=term
       O_tend(k) = O_tend(k) + term


      enddo !nt
      endif !tirrq>0
      enddo !k

!@PLdbg check values at o2 at test lat and lon
      do k=1,kmax

!       if (vrbos) then
!          write(6,'(a,3i7,13e12.4)')'obio_o2(bio):',
!     .      nstep,i,j,p1d(k),termb1,termb2,termb3,termb4,termb5
!     .      ,termb6,ro2c_NH4,ro2c_NO3,rmu4(k,nt),rmu3(k,nt),o21d(k)
!     .      ,O_tend(k)
!        endif

!@PL check term components

!       if (vrbos) then
!          write(6,'(a,3i7,12e12.4)')'obio_o2(bio):',
!     .    nstep,i,j,termb1(k),termb2(k),termb3(k)
!     .           ,termb4(k,1),termb4(k,2),termb4(k,3),termb4(k,4)
!     .           ,termb5(k,1),termb5(k,2),termb5(k,3),termb5(k,4)
!     .           ,termb6(k)
!        endif

      if (ISNAN(o21d(k))) then
          write(6,'(a,3i7,9e12.4)')'obio_o2(bionan):',
     .      nstep,i,j,p1d(k),termb1,termb2,termb3,termb4,termb5
     .      ,termb6,o21d(k),O_tend(k)
        endif
      enddo
!@PLdbg
#endif
        
!@PL: for gas exchange on ocean grid, calculate O2 tendencies from O2 ocean/atm gas exchange
!@PL: xO2 is a constant mole fraction
      
      


        k = 1
        Ts = temp1d(k)
        ps = saln1d(k) !@PL salinity in PSS-78
!       scco2 = 2073.1 - 125.62*Ts + 3.6276*Ts**2 - 0.043219*Ts**3
        !new OCMIP2016 values
        sco2 = 1920.4d0 - 135.6d0*Ts + 5.2122d0*Ts**2 - 0.10939d0*Ts**3 !@PL Schmidt number
     &          + 0.00094743d0*Ts**4

        wssq = wind*wind
        if (sco2.lt.0.) then
          sco2arg=1.d-10
          rkwo2=1.d-10
        else
          sco2arg = (sco2/660.D0)**(-0.5)      !@PL function of schmidt number
          rkwo2 = awan*wssq*sco2arg           !transfer coeff (units of m/s)
        endif


      ta = log((298.15d0-Ts)/(Ts + tf)) !@PL scaled T from Garcia and Gordon, 1992
      ta2 = ta*ta
      ta3 = ta2*ta
      ta4 = ta3*ta
      ta5 = ta4*ta
      ps2 = ps*ps
      tk100 = 100.d0/(Ts + tf)
!@Pldbg
      SLP = ((oAPRESS(i,j)/100.d0)+stdslp) !@PL oAPRESS is pressure anomoly in Pa, stdslp in hPa
!      SLP = stdslp   

      !@ PL
      !@sum following Orr 2017, compute O2sat0, reference O2sat from Garcia and Gordan 1992
      !@+ and convert to solubility function using eq. (8) in Garcia and
      !@+ Gordon 1992, and eq. (16) in Orr 2017
      !@auth Paul Lerner. var list in obio_bioinit
      !@ PL

     

    
   
       O2sat0 = (1.d-06)*exp(5.80818d0 + 3.20684d0*ta
     .         + 4.11890*ta2 + 4.93845d0*ta2 + 4.93845d0*ta3
     .         + 1.01567d0*ta4 + 1.41575d0*ta5 - ps * (0.00701211d0
     .         - 0.00725958d0*ta + 0.00793334d0*ta2 - 0.000554491d0*ta3)
     .         -0.000000132412d0*ps2) !@PL mol/kg


        xo2 = 2.0946d-1 !@PL units mol fraction O2 in dry air

#ifdef pH2O_Ocean
        pH2O = ocnatm%qsavg(i,j) * SLP/1013.25d0
#else
        pH2O = exp(24.4543d0 - 67.4509d0*(tk100) 
     .          - 4.8489d0*log(1.d0/tk100) 
     .        - 0.000544d0 * ps) ! @PL humidity used in correction term, units atm
#endif
        Ksol = O2sat0 / (xO2*((stdslp/1013.25d0)-pH2O)) !@PL units mol/kg/atm
        O2sat = Ksol*((SLP/1013.25d0) - pH2O)*xo2   !@PL mol/kg

#ifdef TRACERS_bio_O2
!@PL this block is for gas exchange with O2 affected by physics and bio/chemistry
        pO2_ij = o21d(1)/Ksol  !@PL partial pressure of O2 at surface
        delto2 = (O2sat-(o21d(1)/1000.d0)) !@PL units  mol/kg
        flxmolkg = (rkwo2*delto2/dp1d(k))   !unit:qs of mol/m3/s
        flxmolkgh = flxmolkg*SECONDS_PER_HOUR !@PL units of mol/kg/hr
!       flxmolm3h = flxmolm3*SECONDS_PER_HOUR !units of mol/m3/hr       July 2016
        term = flxmolkg*1000.D0*pnoice(k)    !@PL units of mmol/kg/s (=mili-mol/kg/s)
        rhs(k,ndimo2,16) = term
        O_tend(k) = O_tend(k) + term

      !flux sign is (atmos-ocean)>0, i.e. positive flux is INTO the ocean
!     @PLdbg
        o2flux= rkwo2*(O2sat-(o21d(1)/1000.d0))*rho1d(k)*pnoice(k)! air-sea o2 flux
     .            *SECONDS_PER_HOUR                             ! mol/m2/hr
     .            *HOURS_PER_DAY*DAYS_PER_YEAR            ! mol/m2/yr


#endif

#ifdef TRACERS_abio_O2
!@PL this block is for gas exchange with abiotic O2
        pabO2_ij = abo21d(1)/Ksol  !@PL partial pressure of O2 at surface
        deltabo2 = (O2sat-(abo21d(1)/1000.d0)) !@PL units  mol/kg
        flxmolkg = (rkwo2*deltabo2/dp1d(k))   !unit:qs of mol/m3/s
        flxmolkgh = flxmolkg*SECONDS_PER_HOUR !@PL units of mol/kg/hr
!       flxmolm3h = flxmolm3*SECONDS_PER_HOUR !units of mol/m3/hr       July 2016
        termab = flxmolkg*1000.D0*pnoice(k)    !@PL units of mmol/kg/s (=mili-mol/kg/s)
        rhs(k,ndimabo2,16) = termab
        Abo_tend(k) = Abo_tend(k) + termab


!      @PLdbg
       abo2flux= rkwo2*(O2sat-(abo21d(1)/1000.d0))*rho1d(k)*pnoice(k)! air-sea o2 flux
     .            *SECONDS_PER_HOUR                             ! mol/m2/hr
     .            *HOURS_PER_DAY*DAYS_PER_YEAR            ! mol/m2/yr
#endif

        if (vrbos) then
            write(*,'(a,3i7,13e12.4)')'obio_o2, fluxdiag:',
     .      nstep,i,j,dp1d(k),Ts,saln1d(k),sco2arg,wssq,rkwo2,
     .      O2sat
#ifdef TRACERS_bio_O2
     .      ,o21d(1),o2flux,O_tend(k)
#endif
#ifdef TRACERS_abio_O2
     .      ,abo21d(1),abo2flux,Abo_tend(k)
#endif
        endif
!@PLdbg
        if (vrbos) then
          write(6,'(a,3i7,13e12.4)')'obio_o2(watson):',
     .      nstep,i,j,Ts,sco2arg,wssq,rkwo2,O2sat
#ifdef TRACERS_bio_O2
     .      ,o21d(1),
     .      rkwo2*(o2sat-(o21d(1)/1000.d0))*rho1d(k),O_tend(k),term
#endif
#ifdef TRACERS_abio_O2
     . ,abo21d(1),
     . rkwo2*(o2sat-(abo21d(1)/1000.d0))*rho1d(k),Abo_tend(k),termab   !this flux should have units mol,o2/m2/s
#endif
        endif





      end subroutine obio_o2



