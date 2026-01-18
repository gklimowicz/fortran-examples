!#define DEBUG 1

      module respauto_physio
!@sum Routines to calculate autotrophic respiration and physiological status
!@+   like water stress, and diagnostics like ci.
!@+   This module is entirely physics-level and has no knowledge of
!@+   ent data structures.
!@+   Respiration routines include those by plant component and plant 
!@+   carbon pool size, and fluxes from:
!@+     leaf dark respiration,
!@+     photosynthetic activity growth respiration,
!@+     maintenance respiration,
!@+     tissue growth respiration (from phenology/growth module, distributed
!@+   over the day at physical time step), 
!@+     and estimation of a plant's daily respiration requirement.
!@+   Physiological routines include water stress routines and other
!@+   miscellaneous.
!@+   Called by module photcondmod, canopy biophysics, and phenology
!@auth N.Y.Kiang

      
      use ent_const
      use ent_pfts
      use FarquharBBpspar !eventually need to consolidate with ent_pfts

      implicit none

      public Rdark
      public Resp_cpool_maint, Resp_plant_day, Resp_can_growth
      public water_stress3

      contains
!################# AUTOTROPHIC RESPIRATION ######################################

      real*8 function Rdark(Vcmax)
!@sum Rdark  Leaf dark respiration, Rd (umol m-2_leaf s-1)
!@+   From S. von Caemmerer (2000) Biochemical Models of Leaf Photosynthesis,
!@+   CSIRO book.
!@auth N.Y.Kiang

      real*8, intent(in) :: Vcmax

      Rdark = 0.015d0 * Vcmax !von Caemmerer book.
      
      end function Rdark
!-----------------------------------------------------------------------------

!      OLD Friend & Kiang canopy foliage respiration.
!      Rd = Respveg(pspar%Nleaf,Tl)  !Old F&K Respveg is not only leaf respir.

!---------------------------------------------------------------------!
      real*8 function Resp_cpool_maint(pft,C,CN,T_k,T_k_10d,facclim) 
     &     Result(R_maint)
!@sum Maintenance respiration for a plant carbon pool of size C (umol/plant/s)
!@+   Based on biomass amount (total N in pool). From CLM3.0.
!@auth N.Y.Kiang
      !C3 vs. C4:  Byrd et al. (1992) showed no difference in maintenance
      ! respiration costs between C3 and C4 leaves in a lab growth study.
      ! Also, maintenance (dark) respiration showed no relation to
      ! leaf nitrogen content (assimilation and growth respiration did 
      ! respond to leaf N content). In lab conditions, leaf dark respiration
      ! was about 1 umol-CO2 m-2 s-1 for an N range of ~70 to 155 mmol-N m-2.
! Re CLM parameters: 
!       CLM calculates this per individual*population/area_fraction
!       to give flux per area of pft cover rather than per ground area.

!Other versions:
!        !*Original CASA *!
!     &       exp(308.56d0*(1/56.02d0 - (1/(T_k-227.13d0)))) * 
!     &       ugBiomass_per_gC/ugBiomass_per_umolCO2
!        !*Acclimation vertical shift*! 56.02 = 10+273.15-227.13.  76.02 = 30+273.15-227.13.
!        !*Acclimation horizontal shift.
!     &       exp(308.56d0*                                     
!     &       (1/56.02d0 
!     &       - (1/(T_k-min(30.d0,max(10.d0,T_k_10d))+10.d0-227.13d0))))
!     &       * ugBiomass_per_gC/ugBiomass_per_umolCO2

      implicit none
      integer :: pft            !Plant functional type.
      real*8 :: C               !g-C/individual 
                                !Can be leaf, stem, or root pools.
      real*8 :: CN              !C:N ratio of the respective pool
      real*8 :: T_k             !Temperature of canopy (Kelvin)
      real*8 :: T_k_10d         !Temperature of air 10-day average (Kelvin)
                                !  Should be canopy temp, but 10-day avg. okay.
      real*8 :: facclim         !frost-hardiness stress factor

      !---Local-------
      real*8,parameter :: k_CLM = 6.34d-07 !(s-1) rate from CLM.
      real*8,parameter :: ugBiomass_per_gC = 2.d6
      real*8,parameter :: ugBiomass_per_umolCO2 = 28.5
!      real*8 :: k_pft !Factor for different PFT respiration rates.

!      if (pfpar(pft)%leaftype.eq.NEEDLELEAF) then
!        k_pft = 2.d0
!      else
!        k_pft = 1.d0
!      endif

      if (T_k>228.15d0) then    ! set to cut-off at 45 deg C 
        R_maint = facclim * pfpar(pft)%r * k_CLM * (C/CN) * 
     &       exp(308.56d0*                                     
     &       (1/min(max(56.02d0,T_k_10d-227.13d0),76.02d0)
     &       - (1/(T_k-227.13d0))))
     &       * ugBiomass_per_gC/ugBiomass_per_umolCO2
      else 
         R_maint = 0.d0
      endif

      end function Resp_cpool_maint
!---------------------------------------------------------------------!
      real*8 function Resp_plant_maint(pft,cpools,TcanopyK,TsoilK,
     &     TairK_10d, TsoilK_10d, facclim, rpools) Result(Rmaintp)
!@sum (kgC/s/plant) Maintenance respiration of a plant.
!@auth N.Y.Kiang
      !NOTES:  Because of module hierarchy, this function cannot access
      !  pspar in FBBphotosynthesis (photconmod), so some biochemical 
      !  parameters are not scaled with current temperature here but are the
      !  maximum, like Vcmax.  Good enough for calculating thresholds.
      !Sapwood:   330 C:N mass ratio from CLM, 
      !   factor 0.5/7=0.0714 relative to foliage from Ruimy et al (1996);
      !    58 from Tatarinov & Cienciala (2006) BIOME-BGC pine live wood,
      !        range 42-73.5 kg-C/kg-N

      implicit none
      integer,intent(in) :: pft
      real*8 :: cpools(N_BPOOLS) !plant carbon pools (gC)
      real*8,intent(in) :: TcanopyK
      real*8,intent(in) :: TsoilK
      real*8,intent(in) :: TairK_10d
      real*8,intent(in) :: TsoilK_10d
      real*8,intent(in) :: facclim
      real*8 :: rpools(N_BPOOLS) !plant respiration by pools (kgC/s/plant)
      !---Local---
      real*8 :: Rd
      real*8 :: Resp_fol, Resp_sw, Resp_lab, Resp_froot
      real*8 :: C2N
      real*8,parameter :: umols_to_kgCs = 0.012d-6

      C2N = 1/(pftpar(pft)%Nleaf*.001d0*pfpar(pft)%SLA)
      Rd = Rdark(pftpar(pft)%Vcmax)*cpools(FOL)*(pfpar(pft)%SLA*.001d0) !Rd*LAI

      !* Maintenance respiration - leaf + sapwood + storage
      Resp_fol = umols_to_kgCs *
     &     (Rd + Resp_cpool_maint(pft, cpools(FOL),C2N,
     &     TcanopyK, TairK_10d,facclim)) !Foliage
      Resp_sw = umols_to_kgCs *
     &     Resp_cpool_maint(pft,cpools(SW), 
     &     330.d0,TcanopyK,TairK_10d,facclim) 
      Resp_lab = 0.d0           !kg-C/m2/s - Storage - NON-RESPIRING
      !* Assume fine root C:N same as foliage C:N
      Resp_froot = umols_to_kgCs * Resp_cpool_maint(pft
     &     ,cpools(FR),C2N,TsoilK,TsoilK_10d,facclim) 

      rpools(FOL) = Resp_fol
      rpools(SW) = Resp_sw
      rpools(HW) = 0.d0
      rpools(FR) = Resp_froot
      rpools(CR) = 0.d0

      Rmaintp =  Resp_froot + Resp_fol + Resp_sw + Resp_lab

#ifdef DEBUG
      write(206,*) cpools(FOL),cpools(SW),cpools(HW),cpools(FR)
     &     ,cpools(CR)
      write(204,*) Rd*umols_to_kgCs, rpools(FOL),rpools(SW),rpools(HW)
     &     ,rpools(FR),rpools(CR), Rmaintp,umols_to_kgCs
#endif

      end function Resp_plant_maint
!---------------------------------------------------------------------!

      real*8 function Resp_plant_day(pft,cpools,TcanopyK,TsoilK,
     &     TairK_10d, TsoilK_10d,facclim) Result(Rauto_day)
!@sum (kgC/day/plant) Estimate of total respiration required by
!@+    a plant for one day, excluding tissue growth respiration.
!@auth N.Y.Kiang
      !  Commented out, because estimate too big:
      !  Light growth respiration is estimated based on GPP at half Vcmax
      !  over the whole day (accounts for tundra plants getting 24-hr light).
!      use photcondmod, only :  pspar
      implicit none
      integer,intent(in) :: pft
      real*8 :: cpools(N_BPOOLS) !plant carbon pools (gC)
      real*8,intent(in) :: TcanopyK
      real*8,intent(in) :: TsoilK
      real*8,intent(in) :: TairK_10d
      real*8,intent(in) :: TsoilK_10d
      real*8,intent(in) :: facclim
      !---Local----
      real*8 :: Rpools(N_BPOOLS) !plant respiration by pools (kgC/s/plant)
      real*8 :: Rmaint
      real*8 :: Rgrowth
      real*8 :: GPPplant !plant GPP (kgC/s/plant)
      real*8 :: LAplant
      real*8, parameter :: s2day = 24.d0*60.d0*60.d0

      Rmaint = Resp_plant_maint(pft,cpools,TcanopyK,TsoilK,
     &     TairK_10d, TsoilK_10d,facclim, Rpools)

!      LAplant = cpools(FOL)*1.d-3*pfpar(pft)%sla !1d-3*gC*m2/kgC = m2
!      GPPplant = facclim*0.5d0*pspar%Vcmax *0.012D-6*LAplant !kgC/s/plant
!      Rgrowth = 0.012D-6*Rdark()*LAplant + 
!     &     Resp_can_growth(pft,GPPplant,Rmaint,0.d0)
      Rgrowth = 0.d0  !Just do maintenance respiration requirement.

      Rauto_day = (Rmaint + Rgrowth)*s2day

#ifdef DEBUG
      write(205,*) Rmaint, Rgrowth, Rauto_day
#endif
      
      end function Resp_plant_day
!---------------------------------------------------------------------!

      real*8 function Resp_root(Tcelsius,froot_kgCm2) Result(Rootresp)
!@sum Frootresp = fine root respiration (kgC/s/m2)
!@+   NOT USED.
      !From ED model.  Not used.
      real*8 :: Tcelsius, froot_kgCm2
      
      Rootresp = OptCurve(Tcelsius,1.0d0,3000.d0) * froot_kgCm2/SECPY
     &     /((1.d0 + exp(0.4d0*(5.0d0-Tcelsius)))
     &     *(1.d0 + exp(0.4d0*(Tcelsius-45.d0))))

      end function Resp_root
!---------------------------------------------------------------------!

      real*8 function OptCurve(Tcelsius,x,y) Result(OptCurveResult)
!@sum Optimum curve, where OptCurveResult=x at 15 Celsius.
      !From ED model.
      real*8 :: Tcelsius, x, y
      
      OptCurveResult = x * exp(y*(1/288.15d0 - 1/(Tcelsius+KELVIN)))
      end function OptCurve

!---------------------------------------------------------------------!
      real*8 function Resp_can_growth(pft,Acan,Rmaint,Rtgrowth) 
     &     Result(R_growth)
!@sum Growth (light) respiration (units as input for Acan and Rmaint).
!@auth N.Y.Kiang
      !Based on photosynthetic activity. See Amthor (2000) review of
      ! Mcree - de Wit - Penning de Vries - Thornley respiration paradigms.
      ! See also Ruimy et al. (1996) analysis of growth_r.
      !Fixed to min 0.d0 like ED2. - NYK
      integer :: pft
      real*8 :: Acan !Canopy photosynthesis rate (mass/m2/s)(or any units)
      real*8 :: Rmaint !Canopy maintenance respiration rate (mass/m2/s)
      real*8 :: Rtgrowth !Growth respiration from tissue growth (mass/m2/s)
      real*8 :: growth_r !pft-dependent. E.g.CLM3.0-0.25, ED2 conifer-0.53, ED2 hw-0.33

!      if (pfpar(pft)%leaftype.eq.NEEDLELEAF) then
      if (pfpar(pft)%woody) then
        growth_r = 0.40d0 !Default 0.4d0, Amthor (2000) range 0.39-0.77
      else
        growth_r = 0.28d0       !0.28 Value from Ruimy et al. (1996)
      endif

      R_growth = max(0.d0, growth_r*(Acan - Rmaint) - Rtgrowth)
      end function Resp_can_growth

!---------------------------------------------------------------------!
      function water_stress(nlayers, soilmp, fracroot, fice,
     &     hwilt, betadl) Result(betad)
!@sum Rosensweig & Abramopoulos (1997) plant water stress function.

      implicit none
      integer,intent(in) :: nlayers !Number of soil layers
      real*8,intent(in) ::  soilmp(:) !Soil matric potential (m)
      real*8,intent(in) :: fracroot(:) !Fraction of roots in layer
      real*8,intent(in) :: fice(:)  !Fraction of ice in layer
      real*8,intent(in) :: hwilt  !Wilting point of pft, matric pot. (m)
      real*8,intent(out) :: betadl(:) !Water stress in layers
      real*8 :: betad !Stress value, 0-1, 1=no stress
      !---Local-----------
      integer :: k
      
      betad = 0.d0
      do k = 1,nlayers
        betadl(k) = (1.d0-fice(k))*fracroot(k)
!     &       *max((hwilt-soilmp(k))/hwilt,0.d0) !R&A original
     &       *min(1.d0,max((hwilt-soilmp(k))/(hwilt + 25.d0),0.d0))  !With unstressed range to h=-25 m.
        betad = betad + betadl(k) 
      end do
      if (betad < EPS2) betad=0.d0

      end function water_stress

!----------------------------------------------------------------------!
      function water_stress2(pft, nlayers, thetas, thetasat, thetamin, 
     &     fracroot, fice, betadl) Result(betad)
!@sum Rodriguez-Iturbe et al. (2001) water stress function.
!@+   Version if input is volumetric soil water content.
!@auth N.Y.Kiang
      !  thetasat = watsat = 0.489d0 - 0.00126d0*sandfrac  !From soilbgc.f

      implicit none
      integer,intent(in) :: pft  !Plant functional type number.
      integer,intent(in) :: nlayers !Number of soil layers
      real*8,intent(in) ::  thetas(:) !Soil vol. water (vol.water/vol.soil)
      real*8,intent(in) :: thetasat  !Saturated soil water (vol.water/vol.soil)
                                !Equals porosity
      real*8,intent(in) :: thetamin !Hygroscopic H2O cont(vol.water/vol.soil)
      real*8,intent(in) :: fracroot(:) !Fraction of roots in layer
      real*8,intent(in) :: fice(:)  !Fraction of ice in layer
      real*8,intent(out) :: betadl(:) !Water stress in layers
      real*8 :: betad !Stress value, 0-1, 1=no stress, weighted by layers

      !Local vars
      integer :: k
      real*8 :: s  !Normalized soil moisture, s=thetas/thetasat
      real*8 :: betak  !Stress value for layer k

      !2. Rodriguez-Iturbe, Laio, & Porporato (2001 set) water stress fn 
      betad = 0.d0
      do k = 1,nlayers
        s = thetas(k)/thetasat
        if (s.ge.pfpar(pft)%sstar) then
          betak = 1.d0  !No stress
        else if ((s.lt.pfpar(pft)%sstar).and.
     &         (s.gt.(pfpar(pft)%swilt))) then
          betak = (s-pfpar(pft)%swilt)/
     &         (pfpar(pft)%sstar-pfpar(pft)%swilt) !Just linear
        else
          betak = 0.d0
        end if
        betadl(k) = (1.d0-fice(k))*fracroot(k)*betak
        betad = betad +  (1.d0-fice(k))*fracroot(k)*betak
      end do
      if (betad < EPS2) betad=0.d0

      end function water_stress2
!----------------------------------------------------------------------!
      function water_stress3(pft, nlayers, thetarel, 
     &     fracroot, fice, betadl) Result(betad)
!@sum Rodriguez-Iturbe et al. (2001) water stress function.
!@+   Version if input is relative soil water content (saturated fraction).
!@auth N.Y.Kiang
      implicit none
      integer,intent(in) :: pft  !Plant functional type number.
      integer,intent(in) :: nlayers !Number of soil layers
!      real*8,intent(in) ::  thetas(:) !Soil vol. water (vol.water/vol.soil)
      real*8,intent(in) ::  thetarel(:) !Relative soil vol. water (vol.water/vol. saturated)
      real*8,intent(in) :: fracroot(:) !Fraction of roots in layer
      real*8,intent(in) :: fice(:)  !Fraction of ice in layer
      real*8,intent(out) :: betadl(:) !Water stress in layers
      real*8 :: betad !Stress value, 0-1, 1=no stress, weighted by layers

      !Local vars
      integer :: k
      real*8 :: s  !Normalized soil moisture, s=thetas/thetasat
      real*8 :: betak  !Stress value for layer k

      !2. Rodriguez-Iturbe, Laio, & Porporato (2001 set) water stress fn 
      betad = 0.d0
      do k = 1,nlayers
        s = thetarel(k)
        if (s.ge.pfpar(pft)%sstar) then
          betak = 1.d0  !No stress
        else if ((s.lt.pfpar(pft)%sstar).and.
     &         (s.gt.(pfpar(pft)%swilt))) then
          betak = (s-pfpar(pft)%swilt)/
     &         (pfpar(pft)%sstar-pfpar(pft)%swilt) !Just linear
        else
          betak = 0.d0
        end if
        betadl(k) = (1.d0-fice(k))*fracroot(k)*betak
        betad = betad +  (1.d0-fice(k))*fracroot(k)*betak
      end do
      if (betad < EPS2) betad=0.d0

      end function water_stress3

!----------------------------------------------------------------------!

      function water_stress4(pft, nlayers, thetarel,
     &     fracroot, fice, betadl) Result(betad)
!@sum Rodriguez-Iturbe et al. (2001) water stress function.
!@+   Version if input is relative soil water content (saturated fraction).
!@+   Version if water access is equal throughout the root zone (not waited by root mass distr.).
!@+      Distribution betadl(:) scaled so that sum betad equals betak stress of least stress layer.
!@auth N.Y.Kiang
      implicit none
      integer,intent(in) :: pft  !Plant functional type number.
      integer,intent(in) :: nlayers !Number of soil layers
!      real*8,intent(in) ::  thetas(:) !Soil vol. water (vol.water/vol.soil)
      real*8,intent(in) ::  thetarel(:) !Relative soil vol. water (vol.water/vol. saturated)
      real*8,intent(in) :: fracroot(:) !Fraction of roots in layer
      real*8,intent(in) :: fice(:)  !Fraction of ice in layer
      real*8,intent(out) :: betadl(:) !Water stress in layers
      real*8 :: betad !Stress value, 0-1, 1=no stress, weighted by layers

      !Local vars
      integer :: k
      real*8 :: s  !Normalized soil moisture, s=thetas/thetasat
      real*8 :: betak  !Stress value for layer k
      real*8 :: rootzone(N_DEPTH)
      real*8 :: betadlsum

      !4. Rodriguez-Iturbe, Laio, & Porporato (2001 set) water stress fn 
      !   with uniform depth access, and stress equal to the least stressed soil layer.
      betad = 0.d0
      do k=1,nlayers
        if (fracroot(k).gt.0.d0) then
           rootzone(k) = 1.d0
        else
           rootzone(k) = 0.d0
        endif
      end do

      do k = 1,nlayers
        s = thetarel(k)
        if (s.ge.pfpar(pft)%sstar) then
          betak = 1.d0  !No stress
        else if ((s.lt.pfpar(pft)%sstar).and.
     &         (s.gt.(pfpar(pft)%swilt))) then
          betak = (s-pfpar(pft)%swilt)/
     &         (pfpar(pft)%sstar-pfpar(pft)%swilt) !Just linear
        else
          betak = 0.d0
        end if
        betadl(k) = (1.d0-fice(k))*rootzone(k)*betak
        betad = max( betad, betadl(k))  !Stress is of least stressed layer
      end do
      betadlsum = sum(betadl(:))
      if ( betadlsum < 1.d-12 ) then  ! to avoid 0/0 divisions
        betadl(:) = 0.d0
        betad = 0.d0
      else
        betadl(:) = betadl(:) * betad/betadlsum  !Scale so that betadl(:) sums to betad.
      endif
      if (betad < EPS2) betad=0.d0

      end function water_stress4

!---------------------------------------------------------------------------

      function calc_Ci_canopy(Ca,Gb, Gs,Anet,LAI,IPAR) Result(ci)
!@sum Foliage internal CO2 conc (mol mol-3) assuming diffusive flux of CO2
!@sum is at steady-state with biochemical uptake by photosynthesis and
!@sum that there is zero leaf boundary layer resistance (infinite gb),
!@sum and that there is no leaf cuticular conductance of CO2.
!@sum Full equation:  ci = ca - Anet*(1.37/gb + 1.65/gs)
!@sum 1.37 = ratio of diffusivities of CO2 and water vapor in laminar flow
!@sum       in the leaf boundary layer
!@sum 1.65 = ratio of diffusivities of CO2 and water vapor in still air at
!@sum       the leaf surface
!@sum (Monteith, 1995;  Kiang, 2003;  Collatz 1991)
!@sum ##### NOTE: Parameter Ball_b should actually be passed in with pspar.
!@sum #####       Have to set up for generic pspar for different photosynthesis
!@sum #####       routines.  (NK)
!@+   FOR DIAGNOSTIC/TESTING, NOT USED FOR REGULAR RUNS.

      implicit none

      real*8,intent(in) :: ca !CO2 mole fraction at surface reference height (umol mol-1)
      real*8,intent(in) :: Gb !Canopy boundary layer conductance of water vapor (mol m-2 s-1)
      real*8,intent(in) :: Gs !Canopy Stomatal conductance of water vapor(mol m-2 s-1)
      real*8,intent(in) :: Anet !Leaf net assimilation of CO2 (umol m-2 s-1)
      real*8,intent(in) :: LAI !Leaf area index 
      real*8,intent(in) :: IPAR !Incident PAR (umol m-2 s-1)
      real*8 :: ci              !Leaf internal CO2 concentration (umol mol-1)
      !----Local------
      real*8,parameter :: MINPARMOL=50  !umol m-2 s-1
      real*8,parameter :: Ball_b = 0.01 !mol m-2 s-1

      if (IPAR.lt.MINPARMOL) then  !Stomates closed
        ci = ca - Anet*1.37/Ball_b
      else
        ci = ca - Anet*(1.37/Gb + 1.65/Gs) !LAI cancels in numerator and denominator.
      endif

      end function calc_Ci_canopy

!***************************************************************************
      end module respauto_physio
