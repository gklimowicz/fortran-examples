!#define  DEBUG  1

      module biophysics !canopygort
!@sum GORT and two-stream canopy radiative transfer (sunlit/shaded 
!@sum leaves), and non-equidistant multiple cohort canopy layering.
!@+   Analytical Clumped Two-Stream (ACTS) model of Ni-Meister et al. (2010).
!@+   Accounts for foliage clumping.  Alternative to canopyspitters.f.
!@+   Scales up leaf-level fluxes.
!@+   UNDER DEVELOPMENT.
!@auth W.Yang, N.Y.Kiang

      use ent_types
      use ent_const
      use ent_pfts
      use photcondmod, only : pscondleaf, ciMIN
      use canopyrad, only : get_canopy_rad
      use FarquharBBpspar
      use ent_debug_mod, only : ent_d 

      implicit none
      
      public photosynth_cond !This is interface for Ent.

!      real*8,parameter :: pi = 3.1415926535897932d0 !@param pi    pi
!      real*8,parameter :: EPS=1.d-8   !Small, to prevent div zero.
      real*8,parameter :: IPARMINMOL=LOW_PAR_LIMIT  !umol m-2 s-1
      real*8,parameter :: O2frac=.20900 !fraction of Pa, partial pressure.


      contains
!################## MAIN SUBROUTINE #########################################
      subroutine photosynth_cond(dtsec, pp)
!@sum photosynth_cond  Main routine to set up drivers and calculate 
!@sum canopy scale fluxes.
!@+   Version that calls Farquhar-Ball-Berry leaf biophysics.
!@+   Calculates photosynthesis, autotrophic respiration, H2O conductance,
!@+   looping through cohorts.
!@+   Inputs:  met drivers, radiation, Ca, Tcanopy
!@+   Outputs:  GPP, NPP, autotrophic respiration components
 
      use ent_const
      use ent_types
      use FarquharBBpspar !pspartype, psdrvtype
      use photcondmod, only : biophysdrv_setup, calc_Pspar,pspar
      use respauto_physio, only : Rdark, water_stress3
      use patches, only : patch_print
      use physutil, only : QSAT
      implicit none

      real*8, intent(in) :: dtsec
      type(patch),pointer :: pp
      !----Local----------------!
      type(cohort),pointer :: cop
      type(psdrvtype) :: psdrvpar !Met biophysics drivers, except for radiation.
      real*8 :: ci_umol !umol/mol, Leaf internal CO2 
      real*8 :: ca_umol !umol/mol, Ambient air CO2
      real*8 :: TsurfK, TcanK, TsoilK, Pa !,rh
      real*8 :: CosZen !,betad
      real*8 :: IPAR            !Incident PAR 400-700 nm (W m-2)
      real*8 :: fdir            !Fraction of IPAR that is direct
      real*8 :: Gb !Leaf boundary layer conductance of water vapor(mol m-2 s-1)
      real*8 :: fdry_pft_eff ! pft-specific effective dry canoopy fraction   
      real*8 :: Anet,Atot,Rd    !umol m-2 s-1
      real*8 :: Iemis ! Nadine's isoprene emission umol m-2 s-1
      real*8 :: GCANOPY,TRANS_SW ! Ci,NPP !,R_auto
      real*8 :: GCANOPYsum, Ciavg, GPPsum, NPPsum, R_autosum,C_labsum,
     &          R_rootsum  !PK 5/15/07
      real*8 :: IPPsum
      real*8 :: molconc_to_umol

      if ( .NOT.ASSOCIATED(pp%tallest)) then ! bare soil
        pp%TRANS_SW = 1.d0
        return
      endif

      if (( pp%tallest%pft.eq.0).or.(pp%tallest%pft > N_PFT)) then
        print *,"photosynth_cond: wrong pft = ", pp%tallest%pft
        call patch_print(6,pp,"ERROR ")
        call stop_model("photosynth_cond: wrong pft",255)
      endif

      !* ZERO SOME OUTPUT VARIABLES AT PATCH LEVEL
      pp%TRANS_SW = 1.d0 !Case of zero LAI.
      !* Time-stepped outputs:  CNC, Ci, Qf.

      !* INITIALIZE SUMMARY OUTPUT VARIABLES *!
      GCANOPYsum = 0.d0
      Ciavg = 0.d0
      GPPsum = 0.d0
      NPPsum = 0.d0
      IPPsum = 0.d0
      R_autosum = 0.d0
      R_rootsum = 0.d0
      C_labsum = 0.d0


      !* SET UP DRIVERS *!
      !* Patch-level water stress only needed for Friend&Kiang conductance.
      !* Cohort-level water stress is used for cohort-level photosynthesis.
!      pp%betad = water_stress(N_DEPTH, pp%cellptr%Soilmp(:)
!     i     ,pp%fracroot(:)
!     i     ,pp%cellptr%fice(:), pfpar(pp%tallest%pft)%hwilt
!     o     , pp%betadl(:))

      !* Radiation drivers *!
      IPAR = pp%cellptr%IPARdir + pp%cellptr%IPARdif
      if (pp%cellptr%IPARdir.eq.0.d0) then
        fdir = 0.d0
      else
        fdir = pp%cellptr%IPARdir / IPAR
      endif
      CosZen = pp%cellptr%CosZen
      Pa = pp%cellptr%P_mbar * 100.d0

      !* Other photosynthesis drivers *!
      !Set up psdrvpar - pack in met drivers.
      Gb = pp%cellptr%Ch*pp%cellptr%U* Pa/
     &     (gasc*(pp%cellptr%TairC+KELVIN)) !m/s * N/m2 * mol K/J * 1/K = mol/m2/s
      molconc_to_umol = gasc * (pp%cellptr%TcanopyC + KELVIN)/Pa * 1d6
      ca_umol = pp%cellptr%Ca * molconc_to_umol  !Convert to umol/mol or ppm.
      ci_umol = 0.7d0*ca_umol !pp%cellptr%Ci * molconc_to_umol  !This is solved for is pscubic in FBBphotosynthesis.f.  Replace with dummy initialization.
      TcanK = pp%cellptr%TcanopyC + KELVIN
      TsurfK = pp%cellptr%TairC + KELVIN
      TsoilK = pp%cellptr%Soiltemp(1) + KELVIN
      call biophysdrv_setup(ca_umol,ci_umol,
     &     pp%cellptr%TcanopyC,Pa,
     &     min(1.d0,  !RH
     &     max(pp%cellptr%Qf,0.d0)/Qsat(TcanK,
     &     2500800.d0 - 2360.d0*(TsurfK-KELVIN),Pa/100.d0)),
     &     psdrvpar)  !Equation for latent heat of vaporizaiton comes from ..?
!#DEBUG
!      write(994,*) ca_umol, ci_umol, pp%cellptr%TcanopyC, Pa, TsurfK,
!     &     pp%cellptr%Qf, QSAT(TcanK,
!     &     2500800.d0 - 2360.d0*(TsurfK-KELVIN),Pa/100.d0)
!----------------------
! STOMATAL SUICIDE TEST
!      call biophysdrv_setup(ca_umol,ci_umol,
!     &     pp%cellptr%TcanopyC,Pa,
!     &     1.d0,!RH for stomotal suicide test
!     &     psdrvpar)
!----------------------
!
!     call get_canopy_rad to calculate
!     1. albedo for the patch
!     2. transmittance for the patch
!     3. some profiles, as foliage, fraction and absorption 
!        using specified layering scheme

      ! print *, 'before get_canopy_rad...'
      if (IPAR*4.05d0.gt.LOW_PAR_LIMIT) then
         call get_canopy_rad(pp, IPAR*4.05d0, fdir)
      else 
         pp%TRANS_SW = 0.d0
      endif
      !* LOOP THROUGH COHORTS *!
      cop => pp%tallest
      do while (ASSOCIATED(cop))
      !* Assign vegpar
         !print *, 'cop%LAI = ', cop%LAI
         if ((cop%LAI.gt.0.d0).and.(IPAR*4.05d0.gt.LOW_PAR_LIMIT))
     &        then
!SOILMOIST_OLD
!          cop%stressH2O = water_stress(N_DEPTH, pp%cellptr%Soilmp(:)
!     i         ,cop%fracroot(:)
!     i         ,pp%cellptr%fice(:), pfpar(cop%pft)%hwilt
!     o         , cop%stressH2Ol(:))

!          cop%stressH2O = 1.d0 !### TEMP RUN HACK FOR MMSF ####!
          !* Can't use water_stress2 until have Soilmoist all layers.
!          cop%stressH2O = water_stress2(cop%pft, N_DEPTH, 
!     i         pp%cellptr%Soilmoist(:), pp%cellptr%soil_Phi, 
!     i         pp%cellptr%soil_dry, 
!     &         cop%fracroot, pp%cellptr%fice(:), cop%stressH2Ol(:))
!          betad = cop%stressH2O

          !KIM - water_stress3 uses Soilmoist as a saturated fraction
          cop%stressH2O = water_stress3(cop%pft, N_DEPTH, 
     i         pp%cellptr%Soilmoist(:), 
     &         cop%fracroot, pp%cellptr%fice(:), cop%stressH2Ol(:))

          !#HACK TEMPORARY
          cop%stressH2O = 1.d0

          call calc_Pspar(dtsec,cop%pft,psdrvpar%Pa,psdrvpar%Tc
     i         ,O2frac*psdrvpar%Pa
     i         ,cop%stressH2O,cop%Sacclim,cop%llspan)

          call canopyfluxes(dtsec, pp, cop
     &         ,Gb
     &         ,psdrvpar
     &         ,GCANOPY,Anet,Atot,Rd !NOTE: Ci should be cohort level
     &         ,Iemis)
!     &         ,TRANS_SW)       !NOTE:  Should include stressH2O.
!     &       ,if_ci)  

          if (pfpar(cop%pft)%leaftype.eq.BROADLEAF) then
            ! stomata on underside of leaves so max stomatal blocking = 0
            fdry_pft_eff = 1.d0
          elseif(pfpar(cop%pft)%leaftype.eq.NEEDLELEAF) then
            ! Bosveld & Bouten (2003) max stomatal blocking = 1/3
            fdry_pft_eff = 1.d0 - min(pp%cellptr%fwet_canopy, 0.333d0)
          else
             fdry_pft_eff = 1.d0
          endif

         !* Assign outputs to cohort *!
         !* Account for wet leaves with pp%cellptr%fwet_canopy.
          cop%GCANOPY = GCANOPY*fdry_pft_eff*(gasc*TsurfK)/Pa  !Convert mol-H2O m-2 s-1 to m/s
          cop%Ci = psdrvpar%ci * !ci is in mole fraction
     &         psdrvpar%Pa/(gasc * (psdrvpar%Tc+KELVIN)) !mol m-3
          cop%GPP = Atot * fdry_pft_eff * 0.012d-6 !umol m-2 s-1 to kg-C/m2-ground/s
          cop%IPP = Iemis * 0.012d-6 !umol m-2 s-1 to kg-C/m2-ground/s

          Anet = cop%GPP - Rd 
       else                     !Zero LAI or night
          cop%GCANOPY=0.d0      !May want minimum conductance for stems.
          cop%Ci = EPS
          cop%GPP = 0.d0
          cop%IPP = 0.d0
          Rd = Rdark(pspar%Vcmax)*cop%LAI
       endif
        !* Update cohort respiration components, NPP, C_lab
       call Respauto_NPP_Clabile(dtsec, TcanK,TsoilK,
     &       pp%cellptr%airtemp_10d+KELVIN,
     &       pp%cellptr%soiltemp_10d+KELVIN, Rd, cop)

!        call Allocate_NPP_to_labile(dtsec, cop) !(Previous)

        ! update total carbon
        cop%C_total = cop%C_total + cop%NPP*dtsec

        !* pp cohort flux summaries
        GCANOPYsum = GCANOPYsum + cop%GCANOPY
        Ciavg = Ciavg + cop%Ci*cop%LAI
        GPPsum = GPPsum + cop%GPP
        IPPsum = IPPsum + cop%IPP
        NPPsum = NPPsum + cop%NPP
        R_autosum = R_autosum + cop%R_auto
        R_rootsum = R_rootsum + cop%R_root  !PK 5/15/07
        C_labsum = C_labsum + cop%C_lab * cop%n !Sum for cohort.

        cop => cop%shorter
      end do

      !* Patch-level OUTPUTS *!
      pp%GCANOPY = GCANOPYsum
      if ( pp%LAI > 0.d0 ) then
        pp%Ci = Ciavg/pp%LAI
      else
        pp%Ci = 0.d0
      endif
      pp%GPP = GPPsum
      pp%IPP = IPPsum
      pp%NPP = NPPsum
      pp%R_auto = R_autosum
      pp%R_root = R_rootsum

      !* Accumulate uptake. 
      !* Respiration should be from leaves and not draw down C_lab. ## Need to allocate respiration to leaves.##
!      pp%C_lab = pp%C_lab + max(C_labsum, 0.d0)  !(kg/m2) ###Eventually need to convert to kg/individual.
      pp%C_lab = C_labsum!(kg/m2) ###Eventually need to convert to kg/individual.


      end subroutine photosynth_cond


!---------------------------------------------------------------------------
      subroutine canopyfluxes(dt, pptr, cop
     i     ,Gb,psp
     o     ,Gs,Anet,Atot,Rd,Iemis)
!     i     ,if_ci)
!@sum canopyfluxes Calculates photosynthesis and conductance with
!@sum Farqhuar et al. (1980) photosynthesis, Ball-Berry stomatal conductance,
!@sum and ACTS (Ni-Meister et al. 2010) canopy radiation (clumped, sunlit, 
!@sum shaded leaves).
!@sum Integrates vertically over the canopy.
!@sum Ci is updated at the canopy level using the canopy boundary layer
!@sum conductance as in Friend and Kiang (2005). 

!@sum If PAR is not directly available, the following conversions may be used:
      !  From total shortwave (W m-2) to PAR (umol m-2 s-1) (Monteith & Unsworth):
      !          PAR(umol m-2 s-1) = 2.3(umol/J)*SW(W m-2)
      !  From PAR (W m-2) to PAR (umol m-2 s-1) (U.Maryland, Dept. of Met., PAR Project),
      !  suggest nominal 485 nm for conversion, which gives:
      !          PAR(umol m-2 s-1) = 4.05(umol/J) * PAR(W m-2)
      !  Dye, D.G. (2004) suggests a slightly different conversion:
      !          PAR(umol m-2 s-1) = 4.56(umol/J) * PAR(W m-2)
      
      implicit none
      real*8,intent(in) :: dt   !time step (seconds)
      type(patch),pointer :: pptr    ! current patch 
      type(cohort),pointer :: cop    ! current cohort
      real*8,intent(in) :: Gb   !Canopy boundary layer conductance of water vapor (mol m-2 s-1)
      type(psdrvtype) :: psp !Photosynthesis drivers, except for radiation.
      real*8,intent(inout) :: Gs !Canopy stomatal conductance of water vapor (mol m-2 s-1)
      real*8,intent(out) :: Anet !Leaf net photosynthesis (micromol m-2 s-1)
      real*8,intent(out) :: Atot !Leaf gross photosynthesis (micromol m-2 s-1)
      real*8,intent(out) :: Rd  !Leaf respiration (umol/m2/s)
      ! real*8,intent(out) :: TRANS_SW !Transmittance of shortwave to ground surface.
      real*8,intent(out) :: Iemis ! Leaf isoprene emission Nadine (micromol m-2 s-1)
      
!Passed parameters
!      type(psdrvtype) :: psp

!----------------------------------------------------------------------!
!Local variables
!nu   real*8, parameter :: EPS=1.D-3
      integer :: L, layers
      real*8 :: Aleaf1 !Mean net photosynthesis at layer bottom (umol[CO2]/m2/s).
      real*8 :: Aleaf2 !Mean net photosynthesis at layer top (umol[CO2]/m2/s).
      real*8 :: gleaf1, gleaf2 !Mean conductance at L1, L2 (umol[H2O]/m2/s).
      real*8 :: Rdleaf1, Rdleaf2 !Mean leaf respiration at L1, L2 (umol[CO2]/m2/s)
      real*8 :: Ileaf1,Ileaf2 ! mean isop emis at L1,L2 (umol[C]/m2/s)
      real*8 :: SUM,SUMg,SUMr,SUMi

      ! print *, 'pptr%crad%LAI(1)=', pptr%crad%LAI(1) 
      if (.NOT.ASSOCIATED(pptr%crad%LAI)) then ! no LAI
         call stop_model('not associated LAI',255)
      endif 

      layers=size(pptr%crad%LAI)
         SUM=0.D0
         SUMg=0.d0
         SUMr=0.d0
         SUMi=0.d0
         do 11 L=2,layers
           call photosynth_sunshd(pptr,cop,L,psp,Gb,Aleaf1,gleaf1,
     &          Rdleaf1,Ileaf1)
!           SUM = SUM + 0.5D0*(cop%fp(L)-cop%fp(L-1))*Aleaf1
!           SUMg = SUMg + 0.5D0*(cop%fp(L)-cop%fp(L-1))*gleaf1
!           SUMr = SUMr + 0.5D0*(cop%fp(L)-cop%fp(L-1))*Rdleaf1
!           SUMi = SUMi + 0.5D0*(cop%fp(L)-cop%fp(L-1))*Ileaf1
           SUM = SUM + (cop%fp(L)-cop%fp(L-1))*Aleaf1
           SUMg = SUMg + (cop%fp(L)-cop%fp(L-1))*gleaf1
           SUMr = SUMr + (cop%fp(L)-cop%fp(L-1))*Rdleaf1
           SUMi = SUMi + (cop%fp(L)-cop%fp(L-1))*Ileaf1
   11    continue

         Atot = SUM
         Rd = SUMr
         Gs = SUMg
         Anet = Atot - Rd
         Iemis = SUMi

!#ifdef DEBUG        
!        write(992,*) CosZen,IPAR,cradpar,psdrvpar
!     &       ,Gb,Gsint,Gs,Atot,Anet,Rd,TRANS_SW
!#endif
      end subroutine canopyfluxes

!################## PHOTOSYNTHESIS #########################################

      subroutine photosynth_sunshd(
!@sum photosynth_sunshd  Calculates sunlit and shaded leaf fluxes.
!@auth N.Y.Kiang
      !new scheme parameters
     i     pptr                 !patch pointer
     i     ,cop                 !cohort pointer
     i     ,L                   !layer (top to bottum, from 2) 
     i     ,psd                 !Photosynthesis met drivers
     i     ,Gb                  !Leaf boundary layer conductance (mol/m2/s)
     o     ,Aleaf              !Leaf Net assimilation of CO2 in layer (umol m-2 s-1)
     o     ,gsleaf            !Leaf Conductance of water vapor in layer (mol m-2 s-1)
     o     ,Rdleaf             !Leaf respiration (umol m-2 s-1)
     o     ,Ileaf)           ! Leaf isoprene emission (umol m-2 s-1)
      implicit none
      type(patch),pointer :: pptr
      type(cohort),pointer :: cop
      integer :: L
      type(psdrvtype) :: psd
      real*8,intent(in) :: Gb
      !real*8,intent(in):: cs,Tl,Pa,rh
      !real*8,intent(inout) :: ci
      real*8,intent(out) :: Aleaf !Flux for single leaf
      real*8,intent(out) :: gsleaf !Conductance for single leaf
      real*8,intent(out) :: Rdleaf !Respiration for single leaf
      real*8,intent(out) :: Ileaf ! Isop emis for single leaf
      !------Local---------
      real*8 fsl  !Fraction of sunlit foliage in layer at Lc (unitless).
      real*8 Isl !PAR absorbed by sunlit foliage at Lc (umol/m[foliage]2/s).
      real*8 Ish !PAR absorbed by shaded foliage at Lc (umol/m[foliage]2/s).
      real*8 Asl !Anet from sunlit leaves (umol/m2/s)
      real*8 Ash !Anet from shaded leaves (umol/m2/s)
      real*8 gssl,gssh !Leaf conductance, sunlit,shaded (mol-H2O/m2/s)
      real*8 Rdsl, Rdsh
      real*8 Iemisl  ! Nadine's isop emis from sunlit leaves (umol/m2/s)
      real*8 Iemiss  ! Nadine's isop emis from shaded leaves (umol/m2/s)
      integer :: sunlitshaded !1-sunlit, 2-shaded

!      call canopy_rad(Lcum,crp,Isl,Ish,fsl)
      Isl = pptr%crad%I_sun(L) 
      Ish = pptr%crad%I_sha(L) 
      fsl = pptr%crad%f_sun(L)
!      write(997,*) 'Lcum,sigma, sqrtexpr,kdf,rhor,kbl,pft,canalbedo,
!     & LAI,solarzen,I0df,I0dr,Isl,Ish,fsl',Lcum,crp,Isl,Ish,fsl

      sunlitshaded = 1
      call pscondleaf(cop%pft,Isl,psd,Gb,gssl,Asl,Rdsl,sunlitshaded,
     & Iemisl)
!      write(992,*) 'shaded'
      sunlitshaded = 2
      call pscondleaf(cop%pft,Ish,psd,Gb,gssh,Ash,Rdsh,sunlitshaded,
     & Iemiss)
      !call Collatz(crp%pft, Isl,cs,Tl,rh, Pa,ci,gssl,Asl)
      !call Collatz(crp%pft, Ish,cs,Tl,rh, Pa,ci,gssh,Ash)
                   
      Aleaf = fsl*Asl + (1.0d0 - fsl)*Ash
      gsleaf = fsl*gssl + (1.0d0 - fsl)*gssh
      Rdleaf = fsl*Rdsl + (1.0d0 - fsl)*Rdsh
      Ileaf = fsl*Iemisl + (1.0d0 - fsl)*Iemiss
!#ifdef DEBUG
!      write(998,*) Lcum,crp,Isl,Ish,fsl,psd,gssl,gssh,Asl,Ash,Rdsl,Rdsh
!#endif
      end subroutine photosynth_sunshd

!################# CANOPY CONDUCTANCE ######################################

      subroutine Gs_bound(dt, LAI, Gsnew, Gsinout)  
!@sum Gs_bound Limit rate of change of Gs (umol m-2 s-1)
!@+   Call to this was commented out - Igor?
!@auth A.Friend
      ! Required change in canopy conductance to reach equilibrium (m/s).
      implicit none
      real*8,intent(in) :: dt, LAI
      real*8,intent(in) :: Gsnew !Canopy conductance of curr time step (mol m-2 s-1)
      real*8,intent(inout) :: Gsinout !Bounded canopy conductance (mol m-2 s-1)
      !---Local----!
      real*8 :: Gsold           !Canopy conductance at prev time step (mol m-2 s-1)
      real*8, parameter :: rhoH2O = 998.2 !Density of water (1000 kg/m^3)
      real*8, parameter :: MW_H2O = 18.015 !Molecular weight of water (g/mol)
      !real*8, parameter :: ghi = 0.006d0*rhoh2o*1000./MW_H2O !Upper limit of gs leaf (mol m-2 s-1)
      !real*8, parameter :: glo = 0.000001d0*rhoH2O*1000./MW_H2O !Lower limit of gs leaf (mol m-2 s-1), See Ball and Berry paper.
      real*8, parameter :: ghi = 333.0 !Conversion from 6 mm s-1 upper limit.(mol m-2 s-1)
      real*8, parameter :: glo = .015 !Temperature grassland. Korner (1994) (mol m-2 s-1)

      real*8 :: dGs, dGs_max 

      dGs=Gsnew-Gsinout
      Gsold = Gsinout
      Gsinout = Gsnew
      !nu Limit Gs change over timestep because of guard cell mechanics (m/s)
      dGs_max=dt*LAI*(ghi-glo)/1800.0D0
      if( dGs.gt.dGs_max) Gsinout = Gsold + dGs_max
      if(-dGs.gt.dGs_max) Gsinout = Gsold - dGs_max
      ! Biological limits of absolute Gs (m/s).
      if(Gsinout.gt.ghi*LAI) Gsinout=ghi*LAI
      if(Gsinout.lt.glo*LAI) Gsinout=glo*LAI

      end subroutine Gs_bound

!---------------------------------------------------------------------------
      function Gs_from_Ci(Anet,Ca,Gb,Gs,Ci,IPAR) Result(gsout)
!@sum Inversion of calc_Ci_canopy
      implicit none
      real*8 :: Anet !Leaf net assimilation of CO2 (umol m-2 s-1)
      real*8 :: ca !Ambient air CO2 mole fraction at surface reference height (umol mol-1)
      real*8 :: gb !Canopy boundary layer conductance of water vapor (mol m-2 s-1)
      real*8 :: gs !Canopy Stomatal conductance of water vapor(mol m-2 s-1)
      real*8 :: ci !Leaf internal CO2 concentration (umol mol-1)
      real*8 :: IPAR !Incident PAR (umol m-2 s-1)
      real*8 :: gsout !Adjusted Gs (mol m-2 s-1)
      !----- Local ------

      if (IPAR.lt.IPARMINMOL) then  !Stomates closed
        gsout = gs
      else
        gsout =  1.65/((Ca - Ci)/(Anet) - 1.37/gb )
      endif

      end function Gs_from_Ci

!################# NPP STORAGE ALLOCATION ######################################
!      subroutine Allocate_NPP_to_labile(dtsec,cop)
!@sum Allocate_NPP_storage.  Allocates C_lab.
!@sum This allows C_lab to go negative at the half-hourly time scale, 
!@sum This does not check for C_lab going negative, because either:
!@sum 1) Respiration_autotrophic adjusts respiration to prevent this from 
!@sum happening (check current setup); or 
!@sum 2) C_lab may go negative for less than 24 hours, and then
!@sum the phenology growth module will do senescence and retranslocation
!@sum to compensate for the loss.
!@sum For prescribed LAI, C_lab provides a measure of the imbalance between 
!@sum the biophysics and prescribed LAI.
!
!      implicit none
!      real*8 :: dtsec
!      type(cohort) :: cop
!
!      cop%C_lab = cop%C_lab + 1000.d0*cop%NPP*dtsec/cop%n !(g-C/individual)
!
!      end subroutine Allocate_NPP_to_labile
!################# AUTOTROPHIC RESPIRATION ###################################

      subroutine Respauto_NPP_Clabile(dtsec,TcanopyK,TsoilK,
     &     TairK_10d, TsoilK_10d, Rd, cop)
!@sum Updates cohort-level autotrophic respiration, NPP, and C_lab
!@sum (All in kg-C/m^2/s)
!@auth N.Y.Kiang
!@+   Note:  This does not check for C_lab going negative, because
!@+   the phenology/growth module compensates on a daily basis
!@+   for negative C_lab by senescence and retranslocation.
!@+   For prescribed LAI, C_lab provides a measure of the imbalance between 
!@+   the biophysics and prescribed LAI.
!@+   For sapwood:
!@+    Resp_cpool_maint(cop%pft,0.0714d0*cop%C_sw, !Sapwood - 330 C:N from CLM, factor 0.5/7=0.0714 relative to foliage from Ruimy et al (1996); 58 from Tatarinov & Cienciala (2006) BIOME-BGC pine live wood range 42-73.5 kg-C/kg-N; 

      use photcondmod, only:  frost_hardiness
      use respauto_physio
      implicit none
      real*8,intent(in) :: dtsec
      real*8,intent(in) :: TcanopyK
      real*8,intent(in) :: TsoilK
      real*8,intent(in) :: TairK_10d
      real*8,intent(in) :: TsoilK_10d
      real*8,intent(in) :: Rd !umol m-2 s-1 Calculated at leaf level in photosynthesis module.
      type(cohort),pointer :: cop
      !----Local-----
      real*8 :: Resp_fol, Resp_sw, Resp_lab, Resp_froot, Resp_maint
      real*8 :: Resp_growth, C2N, Resp_growth_1
      real*8 :: facclim
      real*8 :: umols_to_kgCm2s
!      real*8 :: Rauto, adj

      umols_to_kgCm2s = 0.012D-6 * cop%n !Convert individ flux to canopy.
      facclim = frost_hardiness(cop%Sacclim)
      C2N = 1/(pftpar(cop%pft)%Nleaf*1d-3*pfpar(cop%pft)%SLA)

      !* Maintenance respiration - leaf + sapwood + storage
      Resp_fol = umols_to_kgCm2s *
!!     &       Canopy_resp(vegpar%Ntot, TcanopyC+KELVIN) !Foliage
     &     (Rd + Resp_cpool_maint(cop%pft, cop%C_fol,C2N,
     &     TcanopyK, TairK_10d, facclim)) !Foliage
      Resp_sw = umols_to_kgCm2s *
     &     Resp_cpool_maint(cop%pft,cop%C_sw, 
     &     100.d0,TcanopyK,TairK_10d, facclim) 
      Resp_lab = 0.d0           !kg-C/m2/s - Storage - NON-RESPIRING
      !* Assume fine root C:N same as foliage C:N
      Resp_froot = umols_to_kgCm2s * Resp_cpool_maint(
     &     cop%pft,cop%C_froot,C2N,TsoilK,TsoilK_10d,facclim) 

      Resp_maint = Resp_froot + Resp_fol + Resp_sw + Resp_lab
!     &       Canopy_resp(vegpar%Ntot, TcanopyC+KELVIN))

      !* Growth respiration from biomass tissue growth.
      !Resp_growth_1 = cop%C_growth/(24.d0*3600.d0) !Convert from d-1 to s-1.
      Resp_growth_1 = cop%C_growth_flux
      cop%C_growth = cop%C_growth - Resp_growth_1*dtsec 

      !* Growth respiration tied to GPP; compensates for tissue growth respir.
      Resp_growth = Resp_can_growth(cop%pft, 
     &     cop%GPP,Resp_maint, Resp_growth_1)

      !* NK - commented out -- should be done in phenology with senescence.
      !*** Check for Rauto greater than C_labile + GPP:
      !* If Rauto > C_lab+GPP, then adjust Resp_maint and Resp_growth, but
      !* cannot adjust Resp_growth_1, because biomass growth was already
      !* calculated and checked for adequate C_lab for C_growth.
!      Rauto = Resp_maint + Resp_growth + Resp_growth_1
!      if ((Rauto).gt.(cop%C_lab+cop%GPP)) then !Reduce respiration
!         !if (Rauto.gt.0.d0) then !No div by zero - should not have to check.
!         if ((cop%C_lab +cop%GPP -Resp_growth_1).gt.0.d0) then
!            adj = (cop%C_lab +cop%GPP -Resp_growth_1)/
!     &           (Resp_maint + Resp_growth + EPS)
!         else
!            !call stop_model("Rauto le zero.", 255)
!            !C_lab might go below zero with Resp_growth_1.
!            adj = 0
!         endif
!         Resp_maint = adj * Resp_maint
!         Resp_growth = adj * Resp_growth
!      endif

      !* Update cop respiration, NPP, C_lab.
      cop%R_auto =  Resp_maint + Resp_growth + Resp_growth_1
      cop%R_root = Resp_froot
      cop%NPP = cop%GPP - cop%R_auto !kg-C/m2-ground/s
      cop%C_lab = cop%C_lab + 1000.d0*cop%NPP*dtsec/cop%n !(g-C/individual)

      !set values for debugging
        ent_d%Resp_fol(cop%pft) = Resp_fol
        ent_d%Resp_sw(cop%pft) = Resp_sw
        ent_d%Resp_lab(cop%pft) = Resp_lab
        ent_d%Resp_root(cop%pft) = Resp_froot
        ent_d%Resp_maint(cop%pft) = Resp_maint
        ent_d%Resp_growth_1(cop%pft) = Resp_growth_1
        ent_d%Resp_growth(cop%pft) = Resp_growth

C#define OFFLINE 1
#ifdef OFFLINE
      write(998,*) cop%C_lab,cop%GPP,cop%NPP,Resp_fol,Resp_sw,Resp_lab,
     &     Resp_froot,Resp_maint,Resp_growth, Resp_growth_1
!      write(997,*) cop%C_fol,cop%C_froot,cop%C_sw,cop%C_hw,cop%C_croot
#endif

      end subroutine Respauto_NPP_Clabile
   

!============================================================================

      end module biophysics !canopyspitters
