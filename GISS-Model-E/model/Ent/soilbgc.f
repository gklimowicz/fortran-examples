      module soilbgc
!@sum Routines to simulate soil biogeochemistry:
!@sum microbial dynamics, C & N pools, respiration and N fluxes.

      use ent_const
      use ent_types
      use ent_pfts
      implicit none

      public soil_bgc

      contains
      
!***********************************************************************      
      subroutine soil_bgc(dtsec, pp, soilmoist_in, soiltemp_in)
!@sum soil_bgc  Main routine to interface with driver to calculate 
!@+   soil respiration and update soil carbon pools. 
!@+   Sets up drivers then calls physics routine.
!@auth P.Kharecha, N.Y.Kiang

      use patches, only : print_Tpool
      implicit none

      real*8, intent(in) :: dtsec !main ent time step (s)
      type(patch),pointer :: pp
      real*8, optional :: soilmoist_in, soiltemp_in
      !----Local----------
      real*8 :: Soilmoist(N_CASA_LAYERS) !(vol fraction) by CASA layers.
      integer :: ivt                     !ivt = pft
      real*8 :: Soiltemp(N_CASA_LAYERS)  !(C) by CASA layers.
      real*8 :: clayfrac                 !fractional clay content in soil
      real*8 :: sandfrac                 !fractional clay content in soil
      real*8 :: siltfrac                 !fractional clay content in soil
      real*8 :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS)  !total plant and soil C,N pools (gC/m2)
      real*8 :: Cflux                    !total soil C flux to atm (gC/m2/s) **C flux, NOT CO2!**

      ! do nothing if no vegetation
      if ( .not. ASSOCIATED(pp%tallest) ) return

      ivt = pp%tallest%pft     
!      Soilmoist = pp%Soilmoist    !**soil moist eventually should vary by patch**
!      Soilmoist(:) = pp%cellptr%Soilmoist(:) 
!      Soiltemp(:) = pp%cellptr%Soiltemp(:)  !soil temp, texture vary by cell
      if ( .not. present(soiltemp_in) ) then
        call Soillayer_convert_Ent(pp%cellptr%Soilmoist(:), SOILDEPTH_m, 
     &       Soilmoist) 
        Soilmoist(:) = Soilmoist(:)*pp%cellptr%soil_Phi !Convert from rel.sat. to vol fraction. 
        call Soillayer_convert_Ent(pp%cellptr%Soiltemp(:), SOILDEPTH_m, 
     &       Soiltemp) 

        pp%acc_sbgc_Soilmoist=pp%acc_sbgc_Soilmoist + Soilmoist(1)*dtsec
        pp%acc_sbgc_Soiltemp=pp%acc_sbgc_Soiltemp + Soiltemp(1)*dtsec
      else
        !print *,"called soil_bgc", soilmoist_in, soiltemp_in
        Soilmoist(1) = soilmoist_in
        Soiltemp(1) = soiltemp_in
      endif

      clayfrac = pp%cellptr%soil_texture(3) !in GHY.f, 1-sand,2-loam,3-clay,4-peat(+bedrock) 
      sandfrac = pp%cellptr%soil_texture(1)
!     use siltfrac = 1 - (clayfrac + sandfrac) or 0.4*"loam" for now -PK 6/14/06
      siltfrac = 0.4d0*pp%cellptr%soil_texture(2) !Approximate "silt" from loam.
      Tpool(:,:,:) = pp%Tpool(:,:,:) !Added - NYK 7/27/06

      call casa_bgfluxes(dtsec, Soilmoist,  ivt    
     &     , Soiltemp, clayfrac, sandfrac, siltfrac
     &     , Tpool, Cflux)
      
      !* Convert Cflux from gC/m2/s to kgC/m2/s 
      !* and assign patch soil_resp, Tpool 
      pp%Soil_resp = Cflux*1.d-3 
      pp%Tpool(:,:,:) = Tpool(:,:,:)

      end subroutine soil_bgc

!***********************************************************************
      subroutine casa_bgfluxes(dtsec, Soilmoist, ivt
     i     ,Soiltemp, clayfrac, sandfrac, siltfrac 
     o     ,Tpool,Cflux)
!@sum casa_bgfluxes  Main physics routine to calculate soil respiration
!@+   and update soil carbon pools. 
!@+   From CASA with modified moisture and temperature responses.
!@auth P.Kharecha, N.Y.Kiang
      implicit none

      real*8,intent(in) :: dtsec !main ent time step (s)
      real*8,intent(in) :: Soilmoist(N_CASA_LAYERS) !(vol frac)
      integer,intent(in) :: ivt !ivt = pp%tallest%pft (assigned in soilbgc)
      real*8,intent(in) :: Soiltemp(N_CASA_LAYERS) !(Celsius)
      real*8,intent(in) :: clayfrac !fractional clay content in soil
      real*8,intent(in) :: sandfrac !fractional sand content in soil
      real*8,intent(in) :: siltfrac !fractional silt content in soil
      real*8,intent(inout) :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS) !total plant and soil C,N pools
      real*8,intent(out) :: Cflux !total respiration flux to atm (gC/m2/s)

! ------------------------ local variables ------------------------
      integer ::  n,m
      integer ::  ipool, i
      real*8 :: Cdead(N_CASA_LAYERS) 
!decomp coefs      
      real*8 :: fact_soilmic(N_PFT)
      real*8 :: fact_slow(N_PFT)
      real*8 :: fact_passive(N_PFT)
      real*8 :: eff(NRESP_PATHS)          
      real*8 :: frac_donor(NRESP_PATHS)   
      !real*8 :: kdt(N_PFT,NPOOLS)
      real*8 :: kdt_ivt(NLIVE+1:NPOOLS)
!met variables and functions
      real*8,dimension(N_CASA_LAYERS) :: atmp, bgmoist, bgtemp
      real*8,dimension(N_CASA_LAYERS) :: Wlim !water limitation function for CASA (func of next 2)
#ifdef SBGC_THIS_CODE_IS_UNUSED
      real*8 :: watopt          !"optimal" water content for ET
      real*8 :: watdry          !water content when ET stops (~wilting point)
      !watopt, watdry are funcs of next 3 (which are funcs of clayfrac,sandfrac -- see lsmtci.F) 
      real*8 :: watsat          !saturated volumetric soil water content (porosity)
      real*8 :: smpsat          !soil matric potential at saturation (mm)
      real*8 :: bch             !clapp and hornberger "b"
#endif
      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !amt C lost per pool (gC/m2)
      real*8 :: Resp(PTRACE,NPOOLS,N_CASA_LAYERS) !amt C lost to atm per pool (gC/m2), used to calculate Cflux
      real*8 :: poolsum         !accumulates Resp

!!! HACK  set Closs to something so that 1:NLIVE are initialized
!!! please fix! ??
      Closs(:,:,:) = 0.d0
      
! define rate coefs fact_soilmic, fact_slow, fact_passive (from casatci.F) -PK 5/25/06
      if(ivt.eq.CROPS)then      
         fact_soilmic(ivt) = 1.25d0
         fact_slow(ivt)    = 1.5d0
         fact_passive(ivt) = 1.5d0
      else
         fact_soilmic(ivt) = 1.d0
         fact_slow(ivt)    = 1.d0
         fact_passive(ivt) = 1.d0
      end if
      if(ivt.eq.0)then
         fact_soilmic(ivt) = 0.d0
         fact_slow(ivt)    = 0.d0
         fact_passive(ivt) = 0.d0
      end if

* Maximum RATE CONSTANTS FOR EACH POOL SCALED TO LENGTH OF TIME STEP
* For small delta_t, kdt   is the same as annK*delta_t
*iyf:  Consider dM/dt =  - M/tau
*iyf:  Analytic solution:  M(t) = M0 exp (-t/tau)
*iyf:  Integrate Flux*dt from t=(n-1)*dt to t=n*dt:
*iyf:  integral = M(n*dt) - M[(n-1)*dt]
*iyf:           = - M[(n-1)*dt] {1 - exp [-dt/tau]}
*iyf:    approx = M[(n-1)*dt] * [dt/tau]
*iyf:  variable kdt was previously Krate in CASA

** NOTE: kdt will be used for WOOD and dead pools only
**       so no need to worry about adding stressCD to annK here
      !! computing only those kdt which are needed
      !do n = 1, NPOOLS
      !   do m = 1, N_PFT    
      !      kdt(m,n)= 1.d0 - (exp(-annK(m,n)*dtsec))
      !   enddo
      !enddo

      do n = NLIVE+1, NPOOLS
        kdt_ivt(n)= 1.d0 - (exp(-annK(ivt,n)*dtsec))
      enddo


*---Step 1: heterotrophic (microbial) respiration -------------------------
C.. Initialize respiration fluxes each timestep
C.. Note: these should be over dead pools only (see resp_path_index)
*iyf:  Resp in unit of gC/m2/timestep

      Resp(:,:,:) = 0.d0

*---Step 1a: TEMPERATURE AND MOISTURE CONSTRAINTS ON DECOMP 
! TEMPERATURE DEPENDENCE
      !* Original CASA.
!      bgtemp(:) = (Q10 ** ((Soiltemp(:) - 30.d0) / 10.d0)) !original CASA function -PK
      !* Exponential Arrhenius temperature response (close to Q10 but smaller at 0 Celsius).
!      do i = 1,N_CASA_LAYERS
!        if (Soiltemp(i).le.-33.d0) then
!          bgtemp(i) = 0.d0
!        else
!          bgtemp(i) = exp(308.56d0*(1/63.15d0 - 1/(Soiltemp(i)+33.15))) !exp(308.56d0*(1/(KELVIN+30.d0-240.d0) - 1/(Soiltemp(:)+KELVIN-240.d0)
!        endif
!      enddo
      !* Function f(Tsoil) from DelGrosso et al. (Biogeoch. 73, 2005)**  -PK 2/07
        !allows for variable Q10 rather than fixed 
!            bgtemp(:) = 0.56d0+
!     &      (1.46d0*atan(PI*0.0309d0*(Soiltemp(:)-15.7d0))) / PI

      !* Linear fit to Del Grosso ftemp. Min=0.125 @Soiltemp=0, Max=1.0 @Soiltemp=30 degC- NK
!      bgtemp = max(0.125d0, 
!     &     0.125d0 + (1.d0-0.125d0)/(30.d0-0.d0)*Soiltemp(:)) !linear, no upper cap -NK
!      bgtemp = max(0.125d0, min(1.d0,
!     &     0.125d0 + (1.d0-0.125d0)/(30.d0-0.d0)*Soiltemp(:))) !linear, Del Grosso intercept -NK
!      bgtemp = max(0.d0, min(1.d0,
!     &     0.d0 + (1.d0-0.d0)/(30.d0-0.d0)*Soiltemp(:))) !linear, zero intercept at freezing - NK
      bgtemp = max(0.046d0, min(1.d0,
     &     0.046d0 + (1.d0-0.046d0)/(30.d0-0.d0)*Soiltemp(:))) !linear, log-fitted intercept to original Del Grosso data -NK
      
      !* S-function fit to Del Grosso. - NK
!      bgtemp=1.15d0*(1.d0/(1.d0+EXP(-0.14d0*(Soiltemp(:)-17.d0))))

! MOISTURE DEPENDENCE
* mimic calculation of bevap in surphy.F to get Wlim
* but use Soilmoist,Soiltemp instead of h2osoi,tsoi 
*   watdry = water content when evapotranspiration stops = wp
*   watdry x 0.5d0 = rough estimate of hygroscopic point for microbial respiration.
! Equations take soil textures in percents rather than fractions.
#ifdef SBGC_THIS_CODE_IS_UNUSED
      watsat = 0.489d0 - 0.00126d0*(sandfrac*100.d0)
      smpsat = -10.d0 * ( 10.d0**(1.88d0-
     &     (0.0131d0*sandfrac*100.d0)) )
      bch = 2.91d0 + 0.159d0*(clayfrac*100.d0)
!      watdry = 0.5d0*watsat * (-316230.d0/smpsat) ** (-1.d0/bch)
      watdry = watsat * (-316230.d0/smpsat) ** (-1.d0/bch) !Original WP.
      watopt = watsat * (-158490.d0/smpsat) ** (-1.d0/bch)
#endif
!! 03/11/21 note: there are no limits on Wlim, except if Soiltemp < 0
!!Wlim ultimately used to get total C loss per pool,  
!!to both atm and other pools (see step 1b below) -PK 6/8/06
      do n=1,N_CASA_LAYERS
         if (Soiltemp(n) .gt. 0.d0) then
!     Wlim(n) = min( max(Soilmoist(n)-watdry,0.d0) /  !original CASA function -PK
!     &                   (watopt-watdry), 1.d0)
        !**function RWC from DelGrosso et al., 2005** -PK 2/07
!               Wlim(n) = (Soilmoist(n)-watdry)/(watopt - watdry)
!SOILMOIST_OLD
!            Wlim(n) = max(0.0d0,
!     &           (Soilmoist(n)-watdry)/(watsat - watdry)) !Made this REW instead of Wlim - NK
            Wlim(n) = Soilmoist(n) !Soilmoist is the saturated fraction. - 08/2009 - KIM
         else
            Wlim(n) = 0.01d0
         end if
      end do

!     bgmoist(:) = 0.25d0 + 0.75d0*Wlim(:)  !original CASA function -PK
!**functions f(RWC), Rh from DelGrosso et al. (Biogeoch. 73, 2005)** -PK 2/07
!     bgmoist(:) = 5.d0 *
!     &                 (0.287d0+(atan(PI*0.009d0*(Wlim(:)-17.47d0)))/PI)
      bgmoist(:) = min(1.d0,0.01d0 + 
     &     (1.d0-0.01d0)/(0.7-0.d0)*Wlim(:)) !linear - NK
      atmp(:) = bgtemp(:) * bgmoist(:)
!#DEBUG
!      write(990,*) bgmoist, bgtemp

*---Step 1b: DETERMINE loss of C FROM EACH DEAD POOL (donor) PER TIMESTEP
*iyf:  Closs is the amount of carbon each pool loses in gC/m2/timestep.
*iyf:  A fraction of Closs is transferred to another pool (receiver pool), 
*iyf:  the remainder to the atm (Resp).
*iyf:  The distribution of Closs is done in subroutine casa_respire, Step 1c. 
      do i=1,N_CASA_LAYERS
         do n = 1, NDEAD
            ipool = NLIVE + n
            !Cdead(i) = Tpool(CARBON,ipool,i) * kdt(ivt,ipool) * atmp(i)
            Cdead(i) = Tpool(CARBON,ipool,i) * kdt_ivt(ipool) * atmp(i)
            Closs(CARBON,ipool,i) = Cdead(i)
         end do  

** adjust pools
         Closs(CARBON,SURFSTR,i) = Closs(CARBON,SURFSTR,i) 
     &        * lignineffect(ivt)
         Closs(CARBON,SOILSTR,i) = Closs(CARBON,SOILSTR,i)
     &        * lignineffect(ivt)
         Closs(CARBON,SOILMIC,i) = Closs(CARBON,SOILMIC,i)
     &        * (1.d0-(0.75d0*(siltfrac+clayfrac))) !see Potter et al 1993 for ref -PK 6/8/06
     &        * fact_soilmic(ivt)
         Closs(CARBON,SLOW,i) = Closs(CARBON,SLOW,i)
     &        * fact_slow(ivt)
         Closs(CARBON,PASSIVE,i) = Closs(CARBON,PASSIVE,i)
     &                  * fact_passive(ivt)

*iyf:  limits on loss from dead pools.
ciyf - No need to track limits on loss rate
ciyf (only need to track limits on inventories)
         !do i=1,N_CASA_LAYERS
         do n = NLIVE+1,NPOOLS
            Closs(CARBON,n,i)=MIN(Closs(CARBON,n,i),Tpool(CARBON,n,i))
         end do
         
      end do                    !N_CASA_LAYERS
 
*---Step 1c:  SOM C DECOMPOSITION
*iyf:  update Tpool's:  C is transferred between donor and receiver pools 
*iyf:  Resp is the amount of C to atm, in gC/m2/timestep

* MICROBIAL EFFICIENCIES FOR PARTICULAR FLOWS
      eff( 1) =  0.45           ! SLOW,PASSIVE
      eff( 2) =  0.45           ! SLOW,SOILMIC
      eff( 3) =  0.40           ! SURFMET,SURFMIC
      eff( 4) =  0.40           ! SURFSTR,SURFMIC
      eff( 5) =  0.70           ! SURFSTR,SLOW
      eff( 6) =  0.45           ! SOILMET,SOILMIC
      eff( 7) =  0.45           ! SOILSTR,SOILMIC
      eff( 8) =  0.70           ! SOILSTR,SLOW
      eff( 9) =  0.40           ! CWD,SURFMIC
      eff(10) =  0.70           ! CWD,SLOW
      eff(11) =  0.40           ! SURFMIC,SLOW
      eff(12) =  0.85 - (0.68 * (siltfrac+clayfrac)) ! SOILMIC,PASSIVE 
      eff(13) =  0.85 - (0.68 * (siltfrac+clayfrac)) ! SOILMIC,SLOW 
      eff(14) =  0.45           ! PASSIVE,SOILMIC
* EXTRA RESPIRATION TRANSFER EFFICIENCIES
      frac_donor( 1) =  0.003 + (0.009*clayfrac)
      frac_donor( 2) =  1.0 - frac_donor(1)
      frac_donor( 3) =  1.0
      frac_donor( 4) =  1.0 - structurallignin(ivt)
      frac_donor( 5) =  structurallignin(ivt)
      frac_donor( 6) =  1.0
      frac_donor( 7) =  1.0 - structurallignin(ivt)
      frac_donor( 8) =  structurallignin(ivt)
      frac_donor( 9) =  1.0 - woodligninfract
      frac_donor(10) =  woodligninfract
      frac_donor(11) =  1.0
      frac_donor(12) =  0.003 + (0.032*clayfrac)
      frac_donor(13) =  1.0 - frac_donor(12)
      frac_donor(14) =  1.0

! Determine amount of C transferred from each "donor" pool and where 
! it goes ("receiver" pool), in addition to total C respired (Resp-->Cflux)
      call casa_respire(ivt ,eff ,frac_donor
     &     ,Closs, Resp, Tpool)
     
*---Step 2-----------------------------------------------------------
* CALCULATE NITROGEN POOLS !***keep this for future use*** -PK 8/23/06                                 
!            do n = 1,NPOOLS
!              Tpool(Nitrogen,n)=Tpool(Carbon,n)/CNratio(n)
!            enddo   

!--- Step 2.5 ------------------------------------------------------
! **calculate vertical C transport from each pool (assume closed at top&bottom)
! and adjust each pool accordingly** -PK 5/07  -NO GOOD, DO NOT USE -NK
      !if (N_CASA_LAYERS == 2)  call vertCtransport(dtsec, Tpool) !only applies when multiple layers present

*--- Step 3 --------------------------------------------------------
*  Get Total C Fluxes to atm = Sum over respiring (dead) pools/dtsec
*  Note: Resp for live pools should be zero
*iyf:  Cflux in gC/m2/s
      poolsum = 0.0
      do i=1,N_CASA_LAYERS
         do n = NLIVE+1,NPOOLS
            poolsum=poolsum+Resp(CARBON,n,i)
         end do  
      end do                    !over dead pools -PK
      Cflux=poolsum/dtsec       !total C flux to atm (gC/m2/s)


      Closs(CARBON,:,:)=Closs(CARBON,:,:)/dtsec !g/m2/sec
!     Closs(NITROGEN,:,:)=Closs(NITROGEN,:,:)/dtsec  !not used (yet)

      end subroutine casa_bgfluxes
      
!***********************************************************************
      subroutine casa_respire(ivt ,eff ,frac_donor
     &     ,Closs, Resp, Tpool)
!@sum casa_respire  Calculate soil respiration and transfers between soil pools
      implicit none

! ------------------------ input/output variables -----------------
! input
      integer,intent(in) :: ivt !ivt = pp%tallest%pft (assigned in soilbgc)
      real*8,intent(in) :: eff(NRESP_PATHS) !decomp coefs
      real*8,intent(in) :: frac_donor(NRESP_PATHS) !decomp coefs
      real*8,intent(in) :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !amt C lost per pool (gC/m2)
      real*8,intent(inout) :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS) !total plant and soil C,N pools
      real*8,intent(out) :: Resp(PTRACE,NPOOLS,N_CASA_LAYERS) !C lost to atm per pool (gC/m2)-->Cflux in casa_bgfluxes (step 3)
      
! ------------------------ local variables ------------------------
      integer ::  resp_path_index(2,NRESP_PATHS)   
      integer ::  irtype,n
      integer ::  donor_pool
      integer ::  recvr_pool
      real*8 :: Out(N_CASA_LAYERS)
      
! ----------------------------------------------------------
! set values of parameters/constants used in CASA Respiration
! these are in order in which respiration is called
      resp_path_index = reshape ( 
     1     (/SLOW     ,PASSIVE,
     2     SLOW     ,SOILMIC,
     3     SURFMET  ,SURFMIC,
     4     SURFSTR  ,SURFMIC,
     5     SURFSTR  ,SLOW   ,
     6     SOILMET  ,SOILMIC,
     7     SOILSTR  ,SOILMIC,
     8     SOILSTR  ,SLOW   ,
     9     CWD      ,SURFMIC,
     a     CWD      ,SLOW   ,
     b     SURFMIC  ,SLOW   ,
     c     SOILMIC  ,PASSIVE,
     d     SOILMIC  ,SLOW   ,
     e     PASSIVE  ,SOILMIC/),
     &     (/2,NRESP_PATHS/)   )

! Loop over all respiring pools
      do n=1,N_CASA_LAYERS
         do irtype = 1, NRESP_PATHS
            donor_pool = resp_path_index(1,irtype)
            recvr_pool = resp_path_index(2,irtype)
            Out(n)  = Closs(CARBON,donor_pool,n) * frac_donor(irtype)
            Tpool(CARBON,donor_pool,n) = Tpool(CARBON,donor_pool,n)
     &           - Out(n)
            Tpool(CARBON,recvr_pool,n) = Tpool(CARBON,recvr_pool,n)  
     &           + (Out(n) * eff(irtype))
            Resp(CARBON,donor_pool,n) =  Resp(CARBON,donor_pool,n) 
     &           + Out(n) * (1.d0 - eff(irtype))

            !** 01/10/01 make sure donor pool does not fall below zero
            if(Tpool(CARBON,donor_pool,n).le.0.0d0)then
               Tpool(CARBON,donor_pool,n)=0.d0
            end if
         end do                 !irtype - over all 14 respiring pools
      end do                    !N_CASA_LAYERS

      end subroutine casa_respire

!***********************************************************************
      subroutine Soillayer_convert_Ent(s, depthm, Savglayer) 
!@sum Calculates layer-weighted average of soil variables to produce 
!@+   average values for 0-30 cm and, optional, 30-100 cm.
      implicit none
      real*8 :: s(N_DEPTH)            !Any intensive (non-extensive) variable by LSM soil layer.
      real*8,intent(in) :: depthm(N_DEPTH) !(m) Depths of bottoms of LSM soil layers, increasing.
      real*8, intent(out) :: Savglayer(N_CASA_LAYERS)
      !---Local-----
#ifdef NCASA2
      real*8, parameter :: ENTD(N_CASA_LAYERS) = ( / 0.30d0, 1.0d0 /)
#else
      real*8, parameter :: ENTD(N_CASA_LAYERS) = ( / 0.30d0 /)
#endif
      integer :: k,ek
      real*8 :: sk,ddk,depthup
      real*8 :: savg, dsum
      logical :: entdpassed !Flag for if Ent soil layer depth was passed.

      k=1
      do ek = 1,N_CASA_LAYERS
         entdpassed = .false.
         savg = 0.d0
         dsum = 0.d0
         ddk = depthm(k)
         do while ((.not.entdpassed).and.(k.le.N_DEPTH))
            if (k.eq.1) then !Top layer
               depthup = 0.d0
            else !Lower layers
               depthup = depthm(k-1)
            endif
            if (depthm(k).le.ENTD(ek)) then
               ddk = depthm(k)- depthup
            else    
                ddk = ENTD(ek) - depthup
                entdpassed = .true.
            endif
            savg = savg + ddk*s(k)
            dsum = dsum + ddk
            k = k+1
         end do
         Savglayer(ek) = savg/dsum
      end do

      end subroutine Soillayer_convert_Ent

!***********************************************************************
c Commented out due to lack of parameterization.
c      subroutine vertCtransport(dtsec, Tpool)
c      !calculates inter-layer C transport (based on Baisden et al., GBC 16, 2002) -PK 5/07
c      implicit none
c      
! input 
c      real*8,intent(in) :: dtsec           !main ent time step (s)
! i/o
c      real*8,intent(inout) :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS)  !total plant and soil C,N pools
c! ------------------------ local variables ------------------------
c      integer ::  i, n, numlayers
c      !coefs for vertical C transport (based on Baisden et al., GBC 16, 2002) -PK 5/07        
c      !layer depths(mm) for use w/vertical transport coefs. **need to know these a priori**  -PK 5/07
c      real*8 :: dz(N_CASA_LAYERS) !move to ent_const at some point? -PK 7/07   
c      !rate coefficients 
c      real*8, dimension(N_CASA_LAYERS) :: vfast, vslow, vpass
c      
c      !**temporary workaround to avoid compile trap and allow for N_CASA_LAYERS=2** -PK 7/9/07
c        numlayers = 2
c      !might increase dz(2) to 1700 mm, i.e. put bottom of column at 2 m -PK
c        dz(1) = 300.d0  !layer 1 = 300 mm
c#ifdef NCASA2
c        dz(2) = 700.d0  !layer 2 = 700 mm   
c#endif
c        
c      !calculate 1st-order transport coefs in appropriate units (unitless) 
c      do i=1,N_CASA_LAYERS
c        vfast(i) = 4.d0/dz(i)/SECPY*dtsec   !v1=4.0 mm/yr for fast pool of annual grassland from Baisden et al., 2002 
c        vslow(i) = 0.5d0/dz(i)/SECPY*dtsec  !v2=0.5 mm/yr for slow pool ...
c        vpass(i) = 0.4d0/dz(i)/SECPY*dtsec  !v3=0.4 mm/yr for passive pool ...
c      end do
c      
c      !adjust TPOOL to reflect vertical (inter-layer) C transport
c      do n = NLIVE+1,NPOOLS
c        !top(1st) layer -- loss term only
c        if (n<=SOILMIC .AND. n.ne.CWD) then  !for SURFMET,SURFSTR,SOILMET,SOILSTR,SURFMIC,SOILMIC  
c          TPOOL(CARBON,n,1) =  (1.d0-vfast(1)) * TPOOL(CARBON,n,1) 
c        else if (n==CWD .OR. n==SLOW) then  !for CWD,SLOW
c          TPOOL(CARBON,n,1) =  (1.d0-vslow(1)) * TPOOL(CARBON,n,1)
c        else if (n==PASSIVE) then  !for PASSIVE
c          TPOOL(CARBON,n,1) =  (1.d0-vpass(1)) * TPOOL(CARBON,n,1)
c        end if
c        !second layer -- gain from layer 1 only
c        if (n==6 .OR. n==7 .OR. n==10) then  !for SOILMET,SOILSTR,SOILMIC  
c          TPOOL(CARBON,n,numlayers) = TPOOL(CARBON,n,numlayers) 
c     &                       + vfast(1) * TPOOL(CARBON,n,1)
c        else if (n==8 .OR. n==11) then  !for CWD,SLOW
c          TPOOL(CARBON,n,numlayers) = TPOOL(CARBON,n,numlayers)
c     &                       + vslow(1) * TPOOL(CARBON,n,1)
c        else if (n==12) then  !for PASSIVE
c          TPOOL(CARBON,n,numlayers) = TPOOL(CARBON,n,numlayers)
c     &                       + vpass(1) * TPOOL(CARBON,n,1)
c        end if
c      end do  !loop through all 9 soil C ('dead') pools
c       
c      end subroutine vertCtransport
      
!***********************************************************************

      end module soilbgc
