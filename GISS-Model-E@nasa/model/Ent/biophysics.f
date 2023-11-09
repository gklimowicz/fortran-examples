!@sum Photosynthesis and canopy conductance.
      module biophysics
!@sum Photosynthesis and canopy conductance.
!@auth A. Friend, N. Kiang

      !Ent MODULES TO USE
      use ent_const
      use ent_types
      use ent_pfts

      implicit none
      private
      save

!******************************************************************************
!*  Calculates canopy conductance of water vapor and net assimilation of CO2
!*  by photosynthesis (and related diagnostic, GPP).  
!*
!*  1) Driver program should first call subroutine veg_set_struct to set up canopy
!*  structural parameters.
!*
!*  2) At the surface flux time step, driver program will then call 
!* subroutine veg_conductance:
!*     Canopy conductance, CNC, is given by:
!*           CNC=betad*(1.0D0-0.0075D0*vegpar%vh)*650.0D-6*Anet_max*
!*             ((Ci+0.004D0)/(Ci+EPS))*2.8D0**(-80.0D0*dQs)
!*     where:
!*     betad = water stress factor (0 to 1, 1 unstressed)
!*     vegpar%vh = canopy height (m)
!*     Anet_max = photosynthesis at saturating Ci (leaf internal CO2 conc.)
!*     Ci = leaf internal CO2 concentration, simulated
!*     EPS = dummy small number to prevent divide-by-zero
!*     dQs = vapor pressure deficit (kg[H2O]/kg[air])
!*
!*  2a) The subroutine veg_conductance as of 5/30/04 performs calculations of
!*  radiative transfer, photosynthesis, and canopy conductance.  It calls a generic
!*  integration routine qsimp to calculate photosynthesis (subroutine phot) in layers.
!*  THIS RADIATIVE TRANSFER SCHEME CAN BE REPLACED WITH ANOTHER (see more notes below).
!*
!*  3) Vegetation conductance variables that must be saved/updated between time steps:
!*  (These would be saved into the restart file for a GCM):
!*     CNC is checked from the previous time step to place bounds on how quickly
!*  CNC can change.
!*     Ci is updated given the new CNC  and is saved for the next time step.
!*     Qf (foliage surface vapor mixing ratio, needed for dQs) must be updated 
!*  by the GCM land surface model, since it requires calculation of the total 
!*  surface evapotranspiration.
!*  
!*
!*     CNC is provided back to the surface hydrology and planetary boundary 
!*  layer modules.  The latter are responsible for performing the surface 
!*  energy balance and water balance, calculating total surface fluxes, 
!*  providing temperature (soil, canopy, air), and updating soil moisture 
!*  (to determine betad) and surface air humidity (determine dQs).
!*
!*    5/30/04:  Meteorological drivers that are supplied by other modules are:
!*  PAR, shortwave and longwave incident radiation (top of canopy), wind speed,
!*  atmospheric CO2 concentration.
!*
!*    5/30/04:  CANOPY RADIATIVE TRANSFER.  Currently, photosynthesis is 
!*  calculated by assuming a Beer's Law vertical extinction of radiation 
!*  through the canopy, as well as a similar vertical decline in leaf nitrogen 
!*  dependent on incident light at that canopy level (after Kull and Kruijt, 
!*  1999). 
!*  Integration over canopy layers is performed by increasing the number of 
!*  canopy layers until the total photosynthesis calculated converges.
!******************************************************************************

      public photosynth_cond

      type :: veg_par_type
         !from veg_set_cell: VEGETATION SPECIFICATIONS
         sequence
         !Canopy radiative transfer
         !@var sigma Leaf scattering coefficient (?unitless).
         real*8 :: sigma        !=0.2D0
         real*8  sqrtexpr
         !@var kdf Canopy extinction coeff. for diffuse radiation (unitless).
         real*8 :: kdf          !=0.71D0
         !@var rhor Canopy reflectivity (?unitless).
         real*8 rhor
         !@var kbl Canopy extinction coeff. for direct radiation (unitless).
         real*8 kbl

         !Vegetation geometry, biology
         !@var alai  LAI over whole grid cell
         real*8 :: alai
         !@var nm   Mean canopy nitrogen (g/m2[leaf]) over whole grid cell
         real*8 :: nm
         !@var vh   Mean canopy height (m) over whole grid cell
         real*8 :: vh
         !@var Ntot Total canopy nitrogen (g/m[ground]2).
         real*8 Ntot
         !@var alait  Array of LAI by pft, alait*vfraction=grid LAI contrib.
         !real*8, dimension(8) :: alait 
         !@var vfraction  Array of cover fraction by pft in grid cell
         !real*8, dimension(8) :: vfraction 

         !@var vegalbedo  Vegetation canopy albedo averaged over the grid cell.
         !UPDATE WITH CANOPY RADIATIVE TRANSFER SCHEME
         real*8 :: vegalbedo  
      end type veg_par_type



      type photosynth_par_type
         sequence
         !@var pcp Photorespiratory compensation point (Pa).
         real*8 :: pcp
         !@var Kc Rubisco Michaelis-Menten constant for CO2 (Pa).
         real*8 :: Kc
         !@var Ko Rubisco Michaelis-Menten constant for O2 (Pa).
         real*8 :: Ko
         !@var n1 Proportionality between Jmax and N (umol[CO2]/mmol[N]/s).
         real*8 :: n1
         !@var n2 Proportionality between Vcmax and N (umol[CO2]/mmol[N]/s).
         real*8 :: n2
         !@var m1 CO2 dependence of light-limited photosynthesis (?units).
         real*8 :: m1
         !@var m2 CO2 dependence of light-saturated photosynthesis (?units).
         real*8 :: m2
         !@var msat Nitrogen dependence of photosynthesis (?units).
         real*8 :: msat
         !@var N0 Foliage nitrogen at top of canopy (mmol/m2).
         real*8 :: N0
         !@var Oi Internal foliage O2 partial pressure (kPa).
         real*8 :: Oi !=20.9D0
         !@var k Canopy nitrogen extinction coefficient (?unitless).
         real*8 :: k  !=0.11D0
      end type photosynth_par_type

!     METEOROLOGICAL DRIVERS
      type met_drv_type 
!@var tcan   Canopy temperature (Celsius)
      real*8 :: tcan
!@var qfol  Foliage surface specific humidity (kg vapor/ kg air)
      real*8 :: qfol
!@var dt  Model time step (s)
      real*8 :: dt
!@var pres  Surface air pressure (mbar)
      real*8 :: pres  !Atmospheric pressure (mb)
!@var ch   Ground to surface heat transfer coefficient 
      real*8 :: ch    
!@var U  Surface layer wind speed (m s-1)
      real*8 :: U
!@var parinc Incident photosynthetically active (visible solar, dir+dif)
!@+   radiation on the surface (W/m2) (nyk)
      real*8 :: parinc
!@var fdir Fraction of surface visible radiation that is direct (adf)
      real*8 :: fdir
!@var CosZen cos of Solar zenith angle
      real*8 :: CosZen
!@var Ca Atmospheric CO2 concentration at surface height (mol/m3).
      real*8 :: Ca 
!@var betad  Vegetation water stress (0-1, 1=unstressed)
      real*8 :: betad
      end type met_drv_type

!**** PUBLIC FUNCTIONS*********************************************************
!      public veg_conductance, update_veg_locals

!GLOBAL TO MODULE ONLY:
!      type(veg_par_type) :: vegpar
!----------------------------------------------------------------------------


!Physics variables specific to vegetetation
!@var Ci Original internal foliage CO2 (mol/m3) (adf)
!      real*8 :: Ci
!@var Qf Original foliage surface mixing ratio (kg/kg) (adf)
!      real*8 :: Qf

!----------------------------------------------------------------------------
!OUTPUT:
!@var CNC  Canopy conductance of water vapor (m/s).  
!@+ Must be saved for next time step.  Global so that a GCM driver can
!@+ re-assign a saved value for a grid cell.
!      real*8, public :: CNC

!OUTPUT (via var parameter):
!      real*8 GPP 

!Common block is needed in case of parallelization.--------------------------
!      common /veg_private/
!     &     vegpar
!           ,pres,ch,vsm,parinc,fdir,sbeta
!     &     ,Ci,Qf,Ca,betad,tcan,qv
!     &     ,CNC
!- - - - - - - - - - - - - - - - - - - - -
!     &      alai,nm,nf,vh
!     &     ,pres,ch,vsm
!    &     ,parinc,fdir,vegalbedo,sbeta,Ci,Qf,Ca
!     &     ,betad,tcan,qv
!     &     ,CNC
!$OMP  THREADPRIVATE (/veg_private/)

      contains
      !************************************************************************

      subroutine photosynth_cond(dtsec, pp)
      !Computes photosynthesis and conductance for a patch
      ! by summing through cohorts.
      !1. Canopy radiative transfer is calculate for the entire patch:
      !   canopy albedo, TRANS_SW, light profile
      !   (clumping is calculated outside this routine when canopy
      !    structure is updated).
      !2. Calls veg() to update cohort (cp%) conductance & CO2 uptake:
      !   gcanopy, GPP, Ci
      !3. Calculates autotrophic respiration for each cohort:
      !   NPP, R_auto = R_growth + R_maint(leaf, stem root)
      !4. Update cohort C_lab labile carbon pool.
      !5. Updates the same patch (pp%) variables as above.
                                !
      use ent_const
      use ent_types
      use patches, only : patch_print
      implicit none

      real*8, intent(in) :: dtsec
      type(patch),pointer :: pp  
      type(cohort),pointer :: cop
      !-----Local----------!

      integer :: pft
      real*8 :: TcanopyC,P_mbar,Ch,U,CosZen,Ca,betad,Qf
      real*8 :: GCANOPY,Ci,TRANS_SW,GPP ! ,NPP !,R_auto
      real*8 :: GCANOPYsum, Ciavg, GPPsum, NPPsum, R_autosum,C_labsum,
     &          R_rootsum !, C_reprosum
      real*8 :: IPAR            !Incident PAR 400-700 nm (W m-2)
      real*8 :: fdir            !Fraction of IPAR that is direct
      type(veg_par_type) :: vegpar !Vegetation parameters

#ifdef DEBUG
      print *,"Started photosynth_cond FK" ! with patch:"
      !call patch_print(6,pp," ")
#endif


      if ( .NOT.ASSOCIATED(pp%tallest)) then ! bare soil
        pp%TRANS_SW = 1.d0
        return
      endif

      if (( pp%tallest%pft.eq.0).or.(pp%tallest%pft > N_PFT)) then
        print *,"photosynth_cond: wrong pft = ", pp%tallest%pft
        call patch_print(6,pp,"ERROR ")
        call stop_model("photosynth_cond: wrong pft",255)
      endif

      !* SET UP DRIVERS *!
      !* Patch-level water stress only needed for Friend&Kiang conductance.
      !* Cohort-level water stress is used for cohort-level photosynthesis.
      pp%betad = water_stress(N_DEPTH, pp%cellptr%Soilmp(:)
     i     ,pp%fracroot(:)
     i     ,pp%cellptr%fice(:), pfpar(pp%tallest%pft)%hwilt
     o     , pp%betadl(:))

      IPAR = pp%cellptr%IPARdir + pp%cellptr%IPARdif
      if (pp%cellptr%IPARdir.eq.0.d0) then
        fdir = 0.d0
      else
        fdir = pp%cellptr%IPARdir / IPAR
      endif

      TcanopyC = pp%cellptr%TcanopyC
      P_mbar = pp%cellptr%P_mbar
      Ch = pp%cellptr%Ch
      U = pp%cellptr%U
      CosZen = pp%cellptr%CosZen
      Ca = pp%cellptr%Ca
      betad = pp%betad
      Qf = pp%cellptr%Qf

      !* ZERO SOME OUTPUT VARIABLES AT PATCH LEVEL
      pp%TRANS_SW = 1.d0 !Case of zero LAI.
      !* Time-stepped outputs:  CNC, Ci, Qf.

      !* INITIALIZE SUMMARY OUTPUT VARIABLES *!
      GCANOPYsum = 0.d0
      Ciavg = 0.d0
      GPPsum = 0.d0
      NPPsum = 0.d0
      R_autosum = 0.d0
      R_rootsum = 0.d0
      C_labsum = 0.d0
!      C_reprosum = 0.d0

      !* LOOP THROUGH COHORTS *!
      cop => pp%tallest
      do while (ASSOCIATED(cop))
        !* Assign vegpar

        !* PHOTOSYNTHESIS *!
        if (cop%LAI.gt.0.d0) then
          cop%stressH2O = water_stress(N_DEPTH, pp%cellptr%Soilmp(:)
     i         ,cop%fracroot(:)
     i         ,pp%cellptr%fice(:), pfpar(cop%pft)%hwilt
     o         , cop%stressH2Ol(:))
          betad = cop%stressH2O
          
          vegpar%alai = cop%LAI
!        write(777,*) __FILE__,__LINE__,cop%LAI
          vegpar%nm = cop%nm
          vegpar%vh = cop%h
          vegpar%vegalbedo = pp%albedo(1) !Visible band. NOTE: Patch level.
          pft =  cop%pft

          !GCM land surface saved variables.
          GCANOPY = pp%GCANOPY  !##HACK FOR FRIEND SCHEME.  This only works if GCANOPY is not partitioned over different cohorts.
         !GCANOPY = pp%GCANOPY *cop%LAI/pp%LAI !Fraction contrib. by cohort. ##HACK as it is not necessarily proportional by LAI.
          Ci = pp%Ci

         !print *,"Calling veg..."
          call veg(
     i         dtsec, pft,
     i         TcanopyC,
     i         P_mbar,Ch,U,
     i         IPAR,fdir,CosZen,
     i         Ca,
     i         betad,           !NOTE:  betad is stressH2O of cohort
     i         Qf, 
     &         vegpar,
     &         GCANOPY, Ci, 
     o         TRANS_SW, GPP)   !, NPP )
          
          !* Assign outputs to cohort *!
          cop%GCANOPY = GCANOPY
          cop%Ci = Ci
          cop%GPP = GPP         !kg-C/m2-ground/s
          if (GPP.lt.0.d0) then
            print *,"BAD GPP:",cop%GPP,vegpar%alai,cop%lai
          endif
        else !Zero LAI, no photosynthesis
          vegpar%alai = cop%LAI
          vegpar%Ntot = 0.d0
          cop%GCANOPY = 0.d0 !May want minimum conductance for stems & cuticle.
          cop%Ci = EPS
          GPP = 0.d0
          cop%GPP = GPP
        end if

        !* RESPIRATION FLUXES *!
        call Respiration_autotrophic(dtsec,TcanopyC,cop,vegpar)
        call Allocate_NPP_to_labile(dtsec, cop)

        !* pp cohort flux summaries
        GCANOPYsum = GCANOPYsum + cop%GCANOPY
        Ciavg = Ciavg + Ci*cop%LAI
        GPPsum = GPPsum + cop%GPP
        NPPsum = NPPsum + cop%NPP
        R_autosum = R_autosum + cop%R_auto
        R_rootsum = R_rootsum + cop%R_root  !PK 5/15/07
        C_labsum = C_labsum + cop%C_lab * cop%n !Sum for cohort.
!        C_reprosum = C_reprosum + cop%C_repro * cop%n
        cop => cop%shorter

      end do

      !* Patch-level OUTPUTS *!
      pp%GCANOPY = GCANOPYsum
      pp%Ci = Ciavg/pp%LAI
      pp%GPP = GPPsum
      pp%NPP = NPPsum
      pp%R_auto = R_autosum
      pp%R_root = R_rootsum
      pp%TRANS_SW = TRANS_SW 

      !* Accumulate uptake. 
      !* Respiration should be from leaves and not draw down C_lab. ## Need to allocate respiration to leaves.##
!      pp%C_lab = pp%C_lab + max(C_labsum, 0.d0)  !(kg/m2) ###Eventually need to convert to kg/individual.
      pp%C_lab = C_labsum!(kg/m2) ###Eventually need to convert to kg/individual.
!      pp%C_repro = C_reprosum
      end subroutine photosynth_cond

      !************************************************************************
!NOT USED
!      subroutine veg_set_struct( 
!     i     rsfile)              !Name of restart file
!     o     alaip,               !LAI
!     o     nmp,                 !nitrogen per leaf area (g-N/m^2-leaf)
!     o     nfp,                 !nitrogen parameter
!     o     vhp,                 !canopy height, m
!     o     vegalbedop          !vegetation canopy albedo
!     o )

      ! Called by driver module to set up vegetation structural variables.
      ! Need this for restarts or to set up a grid cell.

!      implicit none
      
      !Set up vegetation structural parameters
!      character*60, intent(in) :: rsfile
      !real*8, intent(out) ::  alaip, nmp, nfp, vhp, vegalbedop

!      open(10,File=rsfile)
!      read(10,*) !Skip header line
!      read(10,*)  vegpar%alai,vegpar%nm,vegpar%vh,vegpar%vegalbedo
!      close(10)

!      end subroutine veg_set_struct

      !************************************************************************

      subroutine veg_set_init(
     i     rsfile,              !Name of restart file
     o     Tc,                  !canopy temperature, Kelvin
     o     Ci,                  !Leaf internal CO2 concentration, mol/m3
     o     Qf,                  !Foliage-to-surface vapor mixing ratio, kg/kg
     o     CNC)                !Conductance water vapor, m/s

      ! Called by driver module to initialize vegetation physics simulated
      ! variables.  Need this for restarts or to set up a grid cell.

      implicit none
      character*60, intent(in) :: rsfile
      real*8, intent(out) :: Tc, Ci, Qf, CNC

      open(10,File=rsfile)
      read(10,*) !Skip header line
      read(10,*) Tc, Ci, Qf, CNC
      close(10)

      !ALTERNATE DEFAULT VALUES
      !Ci = 0.7 * 380e-06 * Pa/(R*TKelvin)  !Gives mol m-3
      !CNC = 0.0006*LAI

      !* DEBUG *!
      !write(10,*) Tc, Ci, Qf, CNC
      end subroutine veg_set_init

      !************************************************************************

      subroutine veg(
     i     dt, pft,tcan,pres,ch,U,parinc,fdir,CosZen,Ca,
     i     betad, 
     i     Qf_IN, 
     &     vegpar,
     &     CNC_INOUT, Ci_INOUT, 
     o     TRANS_SW_OUT, GPP_OUT)    !, NPP_OUT )

!----------------------------------------------------------------------!
! Vegetation dynamics routine (adf, nyk). (May 2004)
! Calculates canopy stomatal conductance (CNC, m/s) using a
! biologically-based formulation, along with canopy CO2 fluxes (GPP;
! kg[C]/m2/s).  Vegetation structure and meteorological drivers for the 
! grid cell should be previously set.  (Eventually should pass these in
! as data structures).
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!@var dt  Time step (seconds)
      real*8, intent(in) :: dt
!@var pft Plant functional type number
      integer, intent(in) :: pft
!@var tcan   Canopy temperature (Celsius)
      real*8, intent(in) :: tcan
!@var pres  Surface air pressure (mbar)
      real*8, intent(in) :: pres  !Atmospheric pressure (mb)
!@var ch   Ground to surface heat transfer coefficient 
      real*8, intent(in) :: ch    
!@var U  Surface layer wind speed (m s-1)
      real*8, intent(in) :: U
!@var parinc Incident photosynthetically active (visible solar, dir+dif)
!@+   radiation on the surface (W/m2) (nyk)
      real*8, intent(in) :: parinc
!@var fdir Fraction of surface visible radiation that is direct (adf)
      real*8, intent(in) :: fdir
!@var CosZen cos of Solar zenith angle
      real*8, intent(in) :: CosZen
!@var Ca Atmospheric CO2 concentration at surface height (mol/m3).
      real*8, intent(in) :: Ca 
!@var betad  Vegetation water stress (0-1, 1=unstressed)
      real*8, intent(in) :: betad

!var Qf_IN Foliage surface H2O vapor mixing ratio (kg[H2O]/kg[air])
      real*8, intent(in) :: Qf_IN
      type(veg_par_type) :: vegpar
!---OUTPUT VARIABLES--------------
!@var TRANS_SW_OUT  Transmittance of shortwave through canopy to soil surface.
      real*8, intent(out) :: TRANS_SW_OUT
!var CNC_inout  Canopy conductance of water vapor (m s-1). 
      real*8, intent(inout) :: CNC_INOUT
!var Ci_INOUT Internal foliage CO2 concentration (mol/m3)
      real*8, intent(inout) :: Ci_INOUT
!@var GPP_OUT  Gross primary productivity (kg[C]/m2/s).
      real*8, intent(out) :: GPP_OUT
!@var NPP_OUT Net primary productivity (kg[C]/m2/s)
!      real*8, intent(out) :: NPP_OUT
!----Local------------------------
!var ps  Photosynthetic parameters (nyk)
      type(photosynth_par_type) :: ps

!################## RADIATIVE TRANSFER ########################################
!@var PAR Incident photosynthetically active radiation (umol/m2/s).
!      real*8 PAR
!@var I0dr Direct beam PAR incident on canopy (umol/m2/s).
      real*8 I0dr
!@var I0df Diffuse PAR incident on canopy (umol/m2/s).
      real*8 I0df
!@var rho Canopy reflectivity (?unitless).
!      real*8 rhor
!@var kbl Canopy extinction coeff. for direct radiation (unitless).
!      real*8 kbl
!#############################################################################
!@var tk Canopy temperature (K).
      real*8 tk
!@var prpa Atmospheric pressure (Pa).
      real*8 prpa
!@var CiPa Internal foliage CO2 partial pressure (Pa).
      real*8 CiPa
!@var Acan Canopy A (umol[CO2]/m[ground]2/s).
      real*8 Acan
!@var Anet Canopy net A (umol[CO2]/m[ground]2/s).
      real*8 Anet
!@var Amax Canopy A at saturating CiPa (umol[CO2]/m[ground]2/s).
      real*8 Amax
!@var Acan Canopy net A at saturating CiPa (umol[CO2]/m[ground]2/s).
      real*8 Anet_max
!@var Rcan Canopy mitochondrial respiration (umol/m2/s).
      real*8 Rcan
!@var dqs Humidity deficit across canopy surface (kg/kg).
      real*8 dQs
!@var CNCN New equilibrium canopy conductance to water (m/s).
      real*8 CNCN
!nu@var dCNC Change in canopy conductance to reach equilibrium (m/s).
      real*8 dCNC
!nu@var dCNC_max Maximum possible change in CNC over timestep (m/s).
      real*8 dCNC_max
!@var gt Conductance from inside foliage to surface height at 30m (m/s).
! NOTE:  Should be changed to surface height at 10 m with new PBL.
      real*8 gt
!----------------------------------------------------------------------------
      real*8 :: sbeta  !Gets cos zenith angle
!@var qv  Canopy saturated specific humidity (kg vapor/ kg air)
      real*8 :: qvsat
      real*8 :: Ci_old, Ci
      real*8 :: dts,dtt
      integer :: N

!## DEBUG  ##!
#ifdef DEBUG
      write(951,*)  dt, pft,tcan,pres,ch,U,parinc,fdir,CosZen,Ca,
     i     betad,
     i     Qf_IN, 
     &     vegpar,
     &     CNC_INOUT, Ci_INOUT, 
     o     TRANS_SW_OUT, GPP_OUT !, NPP_OUT 
#endif

!      write(666,'(e16.5,i5,100e16.5)')      
!     i     dt, pft,tcan,pres,ch,U,parinc,fdir,CosZen,Ca,
!     i     betad,
!     i     Qf_IN,
!     i     vegpar%alai
!      write(667,*)
!     &     vegpar
!     &     CNC_INOUT, Ci_INOUT,
!     o     TRANS_SW_OUT, GPP_OUT)    !, NPP_OUT )


      N = 0
      dts = dt
      dtt = 0.d0
      Ci = Ci_INOUT
      Ci_old = Ci_INOUT
!----------------------------------------------------------------------------
! Make sure there is some leaf area (m2/m2).
      if(vegpar%alai.le.EPS) vegpar%alai=EPS
! Convert canopy temperature from oC to K.
      tk=tcan+tfrz
! Convert atmospheric pressure from mbar to Pa.
      prpa=100.0D0*pres
! Internal foliage CO2 partial pressure from concentration (Pa).
!      CiPa=Ci_INOUT*(gasc*tk)
 10   CiPa=Ci*(gasc*tk)
      ps%Oi = 20.9D0
      ps%k = 0.11D0

!################## RADIATIVE TRANSFER########################################
      !Get incident diffuse and direction PAR radiation 
      !and set up canopy radiative transfer parameters.
      sbeta = CosZen
      call canopy_rad_setup(sbeta, fdir, parinc,
     o     vegpar,I0df, I0dr)
!#############################################################################

! Photorespiratory compensation point (Pa).
      ps%pcp=exp(19.02D0-37.83D3/(gasc*tk))*prpa/1.0D6
! Rubisco Michaelis-Menten constant for CO2 (Pa).
      ps%Kc=exp(38.05D0-79.43D3/(gasc*tk))*prpa/1.0D6
! Rubisco Michaelis-Menten constant for O2 (Pa).
      ps%Ko=exp(20.30D0-36.38D3/(gasc*tk))*prpa/1.0D6
! Proportionality coefficient between Jmax and N (umol[CO2]/mmol[N]/s).
      ps%n1=pfpar(pft)%nf*0.12D0*0.95D+14*exp(-79.5D3/(gasc*tk))/
     1    (1.0D0+exp((650.0D0*tk-199.0D3)/(gasc*tk)))
! Proportionality coefficient between Vcmax and N (umol[CO2]/mmol[N]/s).
      ps%n2=pfpar(pft)%nf*0.23D0*exp(26.35D0-65.33D3/(gasc*tk))
! CO2 dependence of light-limited photosynthesis (?units).
      ps%m1=CiPa/(CiPa+2.0D0*ps%pcp)
! CO2 dependence of light-saturated photosynthesis (?units).
      ps%m2=CiPa/(CiPa+ps%kc*(1.0D0+ps%Oi/ps%ko))
! Nitrogen dependence of photosynthesis (?units).
      ps%msat=min(ps%m1*ps%n1,ps%m2*ps%n2)
! Total canopy nitrogen (g/m[ground]2).
      vegpar%Ntot=vegpar%nm*vegpar%alai
! Top of canopy nitrogen (mmol/m[foliage]2).
      ps%N0=(1000.0D0/14.0D0)*ps%k*vegpar%Ntot/
     &     (1.0D0-exp(-ps%k*vegpar%alai))
! Canopy photosynthesis (Acan: umol/m2/s).
      call qsimp(vegpar,ps,Acan,sbeta,I0dr,I0df,CiPa,tk,prpa)
! Saturating Ci to calculate canopy conductance (mol/m3).
      CiPa=1.0D6
! CO2 dependence of light-limited photosynthesis (?units).
      ps%m1=CiPa/(CiPa+2.0D0*ps%pcp)
! CO2 dependence of light-saturated photosynthesis (?units).
      ps%m2=CiPa/(CiPa+ps%kc*(1.0D0+ps%Oi/ps%ko))
! Nitrogen dependence of photosynthesis (?units).
      ps%msat=min(ps%m1*ps%n1,ps%m2*ps%n2)
! Canopy photosynthesis at saturating CiPa (Amax: umol/m[ground]2/s).
      call qsimp(vegpar,ps,Amax,sbeta,I0dr,I0df,CiPa,tk,prpa)
!----------------------------------------------------------------------!
! Canopy mitochondrial respiration (umol/m[ground]2/s). Above-ground.
!      Rcan=0.20D0*vegpar%Ntot*exp(18.72D0-46.39D3/(gasc*tk))
      Rcan = Canopy_Resp(vegpar%Ntot, tk)
! Net canopy photosynthesis (umol/m[ground]2/s).
      Anet=Acan-Rcan
! Net canopy photosynthesis at saturating CiPa (umol/m[ground]2/s).
      Anet_max=max(Amax-Rcan,0.d0)
!----------------------------------------------------------------------!
! Humidity deficit across canopy surface (inside leaf to leaf surface) (kg/kg).
!      qvsat = QSAT(tk,2500800.-2360.*tk*(101325./prpa)**(gasc/cp),pres) !Wrong previous -NK 9/14/07
!      qvsat = QSAT(tk,LHE0+(SHV100-SHW25)*(tk-KELVIN),pres)
      qvsat = QSAT(tk,LHEW,pres)
      dQs= qvsat-Qf_IN
      if(dQs.lt.zero)dQs=zero
!----------------------------------------------------------------------!
! New equilibrium canopy conductance to moisture (m/s). 
      CNCN=betad*(1.0D0-0.0075D0*vegpar%vh)*650.0D-6*Anet_max*
     &   ((Ci+0.004D0)/(Ci+EPS))*2.8D0**(-80.0D0*dQs)
! Required change in canopy conductance to reach equilibrium (m/s).
      dCNC=CNCN-CNC_INOUT
#ifdef DEBUG      !## DEBUG ##
      write(94,*) betad, vegpar%vh, Ci,dQs,CNCN,dCNC,CNC_INOUT
#endif
!nu Limit CNC change over timestep because of guard cell mechanics (m/s)
!      dCNC_max=dt*vegpar%alai*(0.006D0-0.00006D0)/1800.0D0
 20   N = N + 1
      dCNC_max=dts*vegpar%alai*(0.006D0-0.00006D0)/1800.0D0
      if( dCNC.gt.dCNC_max)CNCN=CNC_INOUT+dCNC_max
      IF(-dCNC.gt.dCNC_max)CNCN=CNC_INOUT-dCNC_max
! Biological limits of absolute CNC (m/s).
      if(CNCN.gt.0.006*vegpar%alai)CNCN=0.006*vegpar%alai
!NOTE:  Water balance issue due to not modeling canopy water content
!       explicitly.  This needs to be considered at some point. -nyk
      if(CNCN.lt.0.00006*vegpar%alai)CNCN=0.00006*vegpar%alai

!----------------------------------------------------------------------!
! Update Ci for next time step:
! Total conductance from inside foliage to surface height at 30m (m/s),
! where CO2 is assumed at Ca.
      gt=1.0D0/(1.42D0/CNCN+1.65D0/(ch*U+EPS))
! Ci update.
      Ci=Ca-1.0D-6*Anet/gt !(mol/m3)

! Limit Cin to physical realism (mol/m3). It is possible that
! oscillations could occur due to the Ci<>CNC feedback. Also, something
! to watch out for is that setting Ci like this does not conserve CO2.
!#ifdef DEBUG
      !* Iterate to match Ci and CNCN
!      if ((Ci.lt.(0.5*Ca)).AND.(dts.gt.1.d0)) then 
!      if ((abs(Ci-Ci_old)/Ci_old.gt.(0.5*Ci_old)).AND.
!      if ((Ci.lt.EPS).AND.(dts.gt.1.d0)) then 
      if ((Ci.lt.EPS).AND.(N.lt.50)) then 
        dts = 0.5*dts
        Ci = Ci_old !Reduce time step until new Ci remains positive.
!        dtt = dtt + dts
        go to 20
      else
        if(Ci.lt.EPS) Ci=EPS 
        Ci_old = Ci
      end if

!#ifdef DEBUG
!      write(94,*) dts,Ci_old,Ci, CNCN, dCNC, dCNC_max,gt,
!     &     Rcan,Acan,Amax, betad,vegpar%vh, dQs,Ca
!#endif

      if ((dtt<dt).AND.(N.lt.50)) then
        dtt = dtt + dts
        dts = dts + dts
!        dts = max(dts, (dt-dtt)/20.d0)
        go to 10  !Recalculate Acan with new Ci
      end if
!#endif

      if(Ci.lt.EPS) Ci=EPS  
!        if ((I0df+I0dr).gt.50.d0) then
!          write(101,*) "1",Anet, Ci, CNCN, gt
!          Anet = Anet - (EPS-Ci)*CNCN/1.42d0
!          Ci = EPS
!          CNCN=1.42d0/((Ca-Ci)/(1.0d-6*Anet) - 1.65d0/(ch*U+EPS))
!          write(101,*) "2",Anet, Ci, CNCN, gt
!        else
!          Ci = EPS
!        end if
!      end if
!      Ci = 0.7*Ca  !## DEBUG HACK
!......................................................................
! NOTE:  Land surface ground hydrology module must update Qf immediately 
!        following this subroutine call, because Qf requires evap_tot 
!        calculated by ground hydrology.  Therefore updated by driver through 
!        update_veg_locals (GISS GCM).
!......................................................................
! Transmission of shortwave through canopy.
      call canopy_transmittance(TRANS_SW_OUT,sbeta,fdir,vegpar)
      if ( TRANS_SW_OUT < 0.d0 .or. TRANS_SW_OUT > 1.d0 ) then
        write(99,*)"veg: unphysical TRANS_SW_OUT", TRANS_SW_OUT
      endif
!......................................................................
! Gross primary productivity (kg[C]/m2/s).
      GPP_OUT=0.012D-6*Acan !Positive for GPP
      !write(*,*) parinc, GPP, Acan
! Net primary productivity (kg[C]/m2/s).  Not available yet.
!     NPP = GPP - 0.012D-6*(Rcan + BoleRespir + RootRespir)
!      NPP_OUT = 0.012D06*Anet
!     &     - Resp_root(Tcan, pp%Tpool(CARBON,FROOT)*1e-3)
! Net ecosystem exchange (kg[C]/m2/s).  Not available yet.
!     NEE = NPP - 0.012D-6* SoilRespir
! Update canopy conductance for next timestep (m/s).
      CNC_INOUT=CNCN
! Update leaf CO2 concentration for next timestep (mol/m3).
      Ci_INOUT=Ci
!----------------------------------------------------------------------!
      return
      end subroutine veg

!***********************************************************************
      subroutine veg_C4(
     i     dt, pft,tcan,pres,ch,U,parinc,fdir,CosZen,Ca,
     i     betad,
     i     Qf_IN,
     &     vegpar,
     &     CNC_INOUT, Ci_INOUT,
     o     TRANS_SW_OUT, GPP_OUT)  !, NPP_OUT )

!----------------------------------------------------------------------!
! Vegetation dynamics routine (adf, nyk). (May 2004, July 2006).
! For C4 vegetation.
! Calculates canopy stomatal conductance (CNC, m/s) using a
! biologically-based formulation, along with canopy CO2 fluxes (GPP;
! kg[C]/m2/s).  Vegetation structure and meteorological drivers for the
! grid cell should be previously set.  (Eventually should pass these in
! as data structures).
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!@var dt  Time step (seconds)
      real*8, intent(in) :: dt
!@var pft Plant functional type number
      integer, intent(in) :: pft
!@var tcan   Canopy temperature (Celsius)
      real*8, intent(in) :: tcan
!@var pres  Surface air pressure (mbar)
      real*8, intent(in) :: pres  !Atmospheric pressure (mb)
!@var ch   Ground to surface heat transfer coefficient
      real*8, intent(in) :: ch
!@var U  Surface layer wind speed (m s-1)
      real*8, intent(in) :: U
!@var parinc Incident photosynthetically active (visible solar, dir+dif)
!@+   radiation on the surface (W/m2) (nyk)
      real*8, intent(in) :: parinc
!@var fdir Fraction of surface visible radiation that is direct (adf)
      real*8, intent(in) :: fdir
!@var CosZen cos of Solar zenith angle
      real*8, intent(in) :: CosZen
!@var Ca Atmospheric CO2 concentration at surface height (mol/m3).
      real*8, intent(in) :: Ca
!@var betad  Vegetation water stress (0-1, 1=unstressed)
      real*8, intent(in) :: betad

!var Qf_IN Foliage surface H2O vapor mixing ratio (kg[H2O]/kg[air])
      real*8, intent(in) :: Qf_IN
      type(veg_par_type) :: vegpar
!---OUTPUT VARIABLES--------------
!@var TRANS_SW_OUT  Transmittance of shortwave through canopy to soil surface.
      real*8, intent(out) :: TRANS_SW_OUT
!var CNC_inout  Canopy conductance of water vapor (m s-1).
      real*8, intent(inout) :: CNC_INOUT
!var Ci_INOUT Internal foliage CO2 concentration (mol/m3)
      real*8, intent(inout) :: Ci_INOUT
!@var GPP_OUT  Gross primary productivity (kg[C]/m2/s).
      real*8, intent(out) :: GPP_OUT
!@var NPP_OUT Net primary productivity (kg[C]/m2/s)
!      real*8, intent(out) :: NPP_OUT
!----Local------------------------
!var ps  Photosynthetic parameters (nyk)
      type(photosynth_par_type) :: ps

!################## RADIATIVE TRANSFER ########################################
!@var PAR Incident photosynthetically active radiation (umol/m2/s).
!      real*8 PAR
!@var I0dr Direct beam PAR incident on canopy (umol/m2/s).
      real*8 I0dr
!@var I0df Diffuse PAR incident on canopy (umol/m2/s).
      real*8 I0df
!@var rho Canopy reflectivity (?unitless).
!      real*8 rhor
!@var kbl Canopy extinction coeff. for direct radiation (unitless).
!      real*8 kbl
!#############################################################################
!@var tk Canopy temperature (K).
      real*8 tk
!@var prpa Atmospheric pressure (Pa).
      real*8 prpa
!@var CiPa Internal foliage CO2 partial pressure (Pa).
      real*8 CiPa
!@var Acan Canopy A (umol[CO2]/m[ground]2/s).
      real*8 Acan
!@var Anet Canopy net A (umol[CO2]/m[ground]2/s).
!
!@var Anet Canopy net A (umol[CO2]/m[ground]2/s).
      real*8 Anet
!@var Amax Canopy A at saturating CiPa (umol[CO2]/m[ground]2/s).
      real*8 Amax
!@var Acan Canopy net A at saturating CiPa (umol[CO2]/m[ground]2/s).
      real*8 Anet_max
!@var Rcan Canopy mitochondrial respiration (umol/m2/s).
      real*8 Rcan
!@var dqs Humidity deficit across canopy surface (kg/kg).
      real*8 dQs
!@var CNCN New equilibrium canopy conductance to water (m/s).
      real*8 CNCN
!nu@var dCNC Change in canopy conductance to reach equilibrium (m/s).
      real*8 dCNC
!nu@var dCNC_max Maximum possible change in CNC over timestep (m/s).
      real*8 dCNC_max
!@var gt Conductance from inside foliage to surface height at 30m (m/s).
! NOTE:  Should be changed to surface height at 10 m with new PBL.
      real*8 gt
!----------------------------------------------------------------------------
      real*8 :: sbeta  !Gets cos zenith angle
!@var qv  Canopy saturated specific humidity (kg vapor/ kg air)
      real*8 :: qvsat
      real*8 :: Ci_old, Ci
      real*8 :: dts,dtt
      integer :: N

!hack!!!
      !Ci_INOUT = 0.02d0

!## DEBUG  ##!
#ifdef DEBUG
      write(95,*)  dt, pft,tcan,pres,ch,U,parinc,fdir,CosZen,Ca,
     i     betad,
     i     Qf_IN,
     &     vegpar,
     &     CNC_INOUT, Ci_INOUT,
     o     TRANS_SW_OUT, GPP_OUT)  !, NPP_OUT
#endif

      N = 0
      dts = dt
      dtt = 0.d0
      Ci = Ci_INOUT
      Ci_old = Ci_INOUT
!----------------------------------------------------------------------------
! Make sure there is some leaf area (m2/m2).
      if(vegpar%alai.le.EPS) vegpar%alai=EPS
! Convert canopy temperature from oC to K.
      tk=tcan+tfrz
! Convert atmospheric pressure from mbar to Pa.
      prpa=100.0D0*pres
! Internal foliage CO2 partial pressure from concentration (Pa).
!      CiPa=Ci_INOUT*(gasc*tk)
 10   CiPa=Ci*(gasc*tk)
      ps%Oi = 20.9D0
      ps%k = 0.11D0

!################## RADIATIVE TRANSFER########################################
      !Get incident diffuse and direction PAR radiation 
      !and set up canopy radiative transfer parameters.
      sbeta = CosZen
      call canopy_rad_setup(sbeta, fdir, parinc,
     o     vegpar,I0df, I0dr)
!#############################################################################

! Photorespiratory compensation point (Pa).
      ps%pcp=exp(19.02D0-37.83D3/(gasc*tk))*prpa/1.0D6
! Rubisco Michaelis-Menten constant for CO2 (Pa).
      ps%Kc=exp(38.05D0-79.43D3/(gasc*tk))*prpa/1.0D6
! Rubisco Michaelis-Menten constant for O2 (Pa).
      ps%Ko=exp(20.30D0-36.38D3/(gasc*tk))*prpa/1.0D6
! Proportionality coefficient between Jmax and N (umol[CO2]/mmol[N]/s).
      ps%n1=pfpar(pft)%nf*0.12D0*0.95D+14*exp(-79.5D3/(gasc*tk))/
     1    (1.0D0+exp((650.0D0*tk-199.0D3)/(gasc*tk)))
! Proportionality coefficient between Vcmax and N (umol[CO2]/mmol[N]/s).
      ps%n2=pfpar(pft)%nf*0.23D0*exp(26.35D0-65.33D3/(gasc*tk))
! CO2 dependence of light-limited photosynthesis (?units).
      ps%m1=CiPa/(CiPa+2.0D0*ps%pcp)
! CO2 dependence of light-saturated photosynthesis (?units).
      ps%m2=CiPa/(CiPa+ps%kc*(1.0D0+ps%Oi/ps%ko))
! Nitrogen dependence of photosynthesis (?units).
      ps%msat=min(ps%m1*ps%n1,ps%m2*ps%n2)
! Total canopy nitrogen (g/m[ground]2).
      vegpar%Ntot=vegpar%nm*vegpar%alai
! Top of canopy nitrogen (mmol/m[foliage]2).
      ps%N0=(1000.0D0/14.0D0)*ps%k*vegpar%Ntot/
     &     (1.0D0-exp(-ps%k*vegpar%alai))
! Canopy photosynthesis (Acan: umol/m2/s).
      call qsimp(vegpar,ps,Acan,sbeta,I0dr,I0df,CiPa,tk,prpa)
! Saturating Ci to calculate canopy conductance (mol/m3).
      !CiPa=1.0D6 ! C3
      CiPa = 1.5D3 !### CHECK UNITS, MAGNITUDE
! CO2 dependence of light-limited photosynthesis (?units).
      ps%m1=CiPa/(CiPa+2.0D0*ps%pcp)
! CO2 dependence of light-saturated photosynthesis (?units).
      ps%m2=CiPa/(CiPa+ps%kc*(1.0D0+ps%Oi/ps%ko))
! Nitrogen dependence of photosynthesis (?units).
      ps%msat=min(ps%m1*ps%n1,ps%m2*ps%n2)
! Canopy photosynthesis at saturating CiPa (Amax: umol/m[ground]2/s).
      call qsimp(vegpar,ps,Amax,sbeta,I0dr,I0df,CiPa,tk,prpa)
!----------------------------------------------------------------------!
! Canopy mitochondrial respiration (umol/m[ground]2/s).- C4
      Rcan=0.33D0*vegpar%Ntot*exp(18.72D0-46.39D3/(gasc*tk))
! Net canopy photosynthesis (umol/m[ground]2/s).
      Anet=Acan-Rcan
! Net canopy photosynthesis at saturating CiPa (umol/m[ground]2/s).
      Anet_max=max(Amax-Rcan,0.d0)
!----------------------------------------------------------------------!
! Humidity deficit across canopy surface (kg/kg).
!      qvsat = QSAT(tk,2500800.-2360.*tk*(101325./prpa)**(gasc/cp),pres) !Wrong previous - NK 09/14/07
      qvsat = QSAT(tk,2500800.-2360.*(tk-KELVIN),pres)
      dQs=qvsat-Qf_IN
      if(dQs.lt.zero)dQs=zero
!----------------------------------------------------------------------!
! New equilibrium canopy conductance to moisture (m/s).
      CNCN=betad*(1.0D0-0.0075D0*vegpar%vh)*650.0D-6*Anet_max*
     &   ((Ci+0.004D0)/(Ci+EPS))*2.8D0**(-80.0D0*dQs)
! Required change in canopy conductance to reach equilibrium (m/s).
      dCNC=CNCN-CNC_INOUT
#ifdef DEBUG      !## DEBUG ##
      write(94,*) betad, vegpar%vh, Ci,dQs,CNCN,dCNC,CNC_INOUT
#endif
!nu Limit CNC change over timestep because of guard cell mechanics (m/s)
!      dCNC_max=dt*vegpar%alai*(0.006D0-0.00006D0)/1800.0D0
 20   N = N + 1
      dCNC_max=dts*vegpar%alai*(0.006D0-0.00006D0)/1800.0D0
      if( dCNC.gt.dCNC_max)CNCN=CNC_INOUT+dCNC_max
      IF(-dCNC.gt.dCNC_max)CNCN=CNC_INOUT-dCNC_max
! Biological limits of absolute CNC (m/s).
      if(CNCN.gt.0.006*vegpar%alai)CNCN=0.006*vegpar%alai
!NOTE:  Water balance issue due to not modeling canopy water content
!       explicitly.  This needs to be considered at some point. -nyk
      if(CNCN.lt.0.00006*vegpar%alai)CNCN=0.00006*vegpar%alai

!----------------------------------------------------------------------!
! Update Ci for next time step:
! Total conductance from inside foliage to surface height at 30m (m/s),
! where CO2 is assumed at Ca.
      gt=1.0D0/(1.42D0/CNCN+1.65D0/(ch*U+EPS))
! Ci update.
      Ci=Ca-1.0D-6*Anet/gt !(mol/m3)

! Limit Cin to physical realism (mol/m3). It is possible that
! oscillations could occur due to the Ci<>CNC feedback. Also, something
! to watch out for is that setting Ci like this does not conserve CO2.
!#ifdef DEBUG
      !* Iterate to match Ci and CNCN
!      if ((Ci.lt.(0.5*Ca)).AND.(dts.gt.1.d0)) then
!      if ((abs(Ci-Ci_old)/Ci_old.gt.(0.5*Ci_old)).AND.
!      if ((Ci.lt.EPS).AND.(dts.gt.1.d0)) then
      if ((Ci.lt.EPS).AND.(N.lt.50)) then
        dts = 0.5*dts
        Ci = Ci_old !Reduce time step until new Ci remains positive.
!        dtt = dtt + dts
        go to 20
      else
        if(Ci.lt.EPS) Ci=EPS
        Ci_old = Ci
      end if

!#ifdef DEBUG
!      write(94,*) dts,Ci_old,Ci, CNCN, dCNC, dCNC_max,gt,
!     &     Rcan,Acan,Amax, betad,vegpar%vh, dQs,Ca
!#endif

      if ((dtt<dt).AND.(N.lt.50)) then
        dtt = dtt + dts
        dts = dts + dts
!        dts = max(dts, (dt-dtt)/20.d0)
        go to 10  !Recalculate Acan with new Ci
      end if
!#endif

      if(Ci.lt.EPS) Ci=EPS  
!        if ((I0df+I0dr).gt.50.d0) then
!          write(101,*) "1",Anet, Ci, CNCN, gt
!          Anet = Anet - (EPS-Ci)*CNCN/1.42d0
!          Ci = EPS
!          CNCN=1.42d0/((Ca-Ci)/(1.0d-6*Anet) - 1.65d0/(ch*U+EPS))
!          write(101,*) "2",Anet, Ci, CNCN, gt
!        else
!          Ci = EPS 
!        end if
!      end if 
!      Ci = 0.7*Ca  !## DEBUG HACK 
!......................................................................
! NOTE:  Land surface ground hydrology module must update Qf immediately
!        following this subroutine call, because Qf requires evap_tot 
!        calculated by ground hydrology.  Therefore updated by driver through
!        update_veg_locals (GISS GCM).
!......................................................................
! Transmission of shortwave through canopy.
      call canopy_transmittance(TRANS_SW_OUT,sbeta,fdir,vegpar)
      if ( TRANS_SW_OUT < 0.d0 .or. TRANS_SW_OUT > 1.d0 ) then
        write(99,*)"veg: unphysical TRANS_SW_OUT", TRANS_SW_OUT
      endif 
!......................................................................
! Gross primary productivity (kg[C]/m2/s). 
      GPP_OUT=0.012D-6*Acan !Acan[umol/m2/s] to kg-C/m2/s, (umol/m2/s) x MW-C(12g/mol) * (1e-6 mol/umol) * (1e-3 kg/g)
      !write(*,*) parinc, GPP, Acan
! Net primary productivity (kg[C]/m2/s).  Not available yet.
!     NPP = GPP - 0.012D-6*(Rcan + BoleRespir + RootRespir)
!      NPP_OUT = 0.012D-6*Anet
!     &     - Resp_root(tcan, pp%Tpool(CARBON,FROOT)*1e-3) !kg-C/m2/s
! Updata canopy conductance for next timestep (m/s).
      CNC_INOUT=CNCN 
! Updata leaf CO2 concentration for next timestep (mol/m3).
      Ci_INOUT=Ci
!----------------------------------------------------------------------!
      return
      end subroutine veg_C4

!---------------------------------------------------------------------!
      real*8 function Canopy_Resp(N_gm2,T_kelvin) Result(MaintResp)
      !Canopy maintenance (mitochondrial) respiration (umol-C/m2/s)
      !Based on biomass (or canopy nitrogen). Friend & Kiang (2005).
      real*8 :: N_gm2 !Canopy foliage nitrogen (g-N/m2-ground)
      real*8 :: T_kelvin !Canopy temperature (Kelvin)

      MaintResp=0.20D0*N_gm2*exp(18.72D0-46.39D3/(gasc*T_kelvin))
      end function Canopy_Resp

!---------------------------------------------------------------------!
!################# NPP STORAGE ALLOCATION ######################################
      subroutine Allocate_NPP_to_labile(dtsec,cop)
!@sum Allocate_NPP_storage.  Allocates C_lab.
!@sum If NPP is negative then C_lab is reduced by all of NPP, and if C_lab is
!@sum is not enough, then the other biomass pools are reduced to compensate.
!@sum Does not adjust LAI - TEMPORARY HACK.
      implicit none
      real*8 :: dtsec
      type(cohort) :: cop
      !---Local-----
      real*8 :: Cdiff

      !* Accumulate uptake.*!
      if ( (cop%NPP*dtsec/cop%n.lt.0.d0).and.
     &     (abs(cop%NPP*dtsec/cop%n).ge.cop%C_lab) ) then 
        !Don't let C_lab go below zero.
        cop%C_lab = EPS
        !* NYK - TEMPORARY HACK - SHOULD CALL STRESS SUBROUTINE FOR SENESCENCE.
        !*       Needs checks for negative NPP loss greater than C pools.
        Cdiff = cop%NPP*dtsec/cop%n + (cop%C_lab - EPS)
        cop%C_fol = cop%C_fol + 0.33d0*Cdiff
        cop%C_froot = cop%C_froot + 0.33d0*Cdiff
        cop%C_sw = cop%C_sw + 0.33d0*Cdiff
        !cop%LAI = !should update
      else
        if (cop%NPP.gt.0.d0) then
          cop%C_lab = cop%C_lab + 0.8d0*cop%NPP*dtsec/cop%n !(kg/individual)
          cop%pptr%Reproduction(cop%pft) = 
     &         cop%pptr%Reproduction(cop%pft) + 0.2d0*cop%NPP*dtsec !(kg/m2-patch) !Reprod. fraction in ED is 0.3, in CLM-DGVM 0.1, so take avg=0.2.
        else !Negative NPP is taken only from C_lab
          cop%C_lab = cop%C_lab + cop%NPP*dtsec/cop%n !(kg/individual)          
        endif
      endif
      end subroutine Allocate_NPP_to_labile
!################# AUTOTROPHIC RESPIRATION ######################################

      subroutine Respiration_autotrophic(dtsec,TcanopyC,cop,vegpar)
      !@sum Autotrophic respiration - updates cohort respiration,NPP,C_lab
      !@sum Returns kg-C/m^2/s

      implicit none
      real*8,intent(in) :: dtsec
      real*8,intent(in) :: TcanopyC
      type(cohort),pointer :: cop
      type(veg_par_type) :: vegpar
      !----Local-----
      real*8 :: Resp_fol, Resp_sw, Resp_lab, Resp_root, Resp_maint
      real*8 ::Resp_growth, C2N, TcanopyK

      !NOTE: NEED TO FIX Canopy maintenance respiration for different
      !C:N ratios for the different pools.

      TcanopyK = TcanopyC + Kelvin
      !C2N = 1/(pftpar(cop%pft)%Nleaf*1d-3*pfpar(cop%pft)%SLA)
      C2N = 1/(vegpar%nm*1e-3 * pfpar(cop%pft)%SLA)

      !* Maintenance respiration - leaf + sapwood + storage
      Resp_fol = 0.012D-6 * !kg-C/m2/s
!     &       Canopy_resp(vegpar%Ntot, TcanopyC+KELVIN) !Foliage
     &     Resp_can_maint(cop%pft, cop%C_fol,C2N,
     &     TcanopyK,cop%n) !Foliage
      Resp_sw = 0.012D-6 *  !kg-C/m2/s
     &     Resp_can_maint(cop%pft,0.0714d0*cop%C_sw, !Sapwood - 330 C:N from CLM, factor 0.5/7=0.0714 relative to foliage from Ruimy et al (1996)
     &     330.d0,TcanopyK,cop%n) 
      Resp_lab = 0.012D-6 * !kg-C/m2/s
     &     Resp_can_maint(cop%pft,cop%C_lab, !Storage
     &     C2N,TcanopyK,cop%n)
      !* Assume root C:N same as foliage C:N
      Resp_root = 0.012D-6 * Resp_can_maint(cop%pft,cop%C_froot,
     &     C2N,TcanopyK,cop%n)
      Resp_maint = Resp_root + Resp_fol + Resp_sw + Resp_lab
!     &       Canopy_resp(vegpar%Ntot, TcanopyC+KELVIN))

      !* Growth respiration
      Resp_growth = 0.012D-6 * Resp_can_growth(cop%pft, 
     &     cop%GPP/0.012D-6,(Resp_maint)/0.012D-6 )

      !* Total respiration : maintenance + growth
      cop%R_auto =  Resp_maint + Resp_growth
      cop%R_root = Resp_root
      cop%NPP = cop%GPP - cop%R_auto !kg-C/m2-ground/s

!      write(998,*) Resp_fol,Resp_sw,Resp_lab,Resp_root,Resp_maint
!     &     ,Resp_growth
      end subroutine Respiration_autotrophic

!---------------------------------------------------------------------!
      real*8 function Resp_can_maint(pft,C,CN,T_k,n) 
     &     Result(R_maint)
      !Canopy maintenance respiration (umol/m2-ground/s)
      !Based on biomass amount (total N in pool). From CLM3.0.
      implicit none
      integer :: pft            !Plant functional type.
      real*8 :: C               !kg-C/individual  !not per m2 -PK 5/15/07
                                !Can be leaf, stem, or root pools.
      real*8 :: CN              !C:N ratio of the respective pool
!      real*8 :: R_maint !Canopy maintenance respiration rate (umol/m2/s)
      real*8 :: T_k             !Temperature of canopy (Kelvin)
      real*8 :: n               !Density of individuals (no./m2)
      !---Local-------
      real*8,parameter :: k_CLM = 6.34d-07 !(s-1) rate from CLM.
      
      if (T_k>228.15d0) then ! set to a cut-off at 45 deg C
         R_maint = pfpar(pft)%r * k_CLM * (C/CN) * 
     &        exp(308.56d0*(1/56.02d0 - (1/(T_k-227.13d0)))) *
     &        2.d6*n/28.5d0
      else
         R_maint = 0.d0
      endif
      !Note:  CLM calculates this per individual*population/area_fraction
      !      to give flux per area of pft cover rather than per ground area.
      end function Resp_can_maint
!---------------------------------------------------------------------!
      real*8 function Resp_root(Tcelsius,froot_kgCm2) Result(Rootresp)
!@sum Frootresp = fine root respiration (kgC/s/m2)
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
      real*8 function Resp_can_growth(pft,Acan,Rmaint) Result(R_growth)
      !Canopy growth respiration (umol/m2/s)
      !Based on photosynthetic activity. From CLM3.0.
      !Fixed to min 0.d0 like ED2. - NYK
      integer :: pft
      real*8 :: Acan !Canopy photosynthesis rate (umol/m2/s)
      real*8 :: Rmaint !Canopy maintenance respiration rate (umol/m2/s)
      real*8 :: growth_r !pft-dependent. E.g.CLM3.0-0.25, ED2 conifer-0.53, ED2 hw-0.33

      growth_r = pfpar(pft)%r * 0.5d0    !Factor 0.5 gives growth_r=0.3 for C3 grass.
      R_growth = max(0.d0, growth_r * (Acan - Rmaint))
      end function Resp_can_growth


!---------------------------------------------------------------------!
      function water_stress(nlayers, soilmp, fracroot, fice,
     &     hwilt, betadl) Result(betad)
      !1. Rosensweig & Abramopoulos water stress fn.

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
     &       *max((hwilt-soilmp(k))/hwilt,0.d0)
        betad = betad + betadl(k) 
      end do
      if (betad < EPS2) betad=0.d0

      end function water_stress
!----------------------------------------------------------------------!
      function water_stress2(pft, nlayers, thetas, thetasat, thetamin, 
     &     fracroot, fice, hwilt) Result(betad)

      implicit none
      integer,intent(in) :: pft  !Plant functional type number.
      integer,intent(in) :: nlayers !Number of soil layers
      real*8,intent(in) ::  thetas(:) !Soil vol. water (vol.water/vol.soil)
      real*8,intent(in) :: thetasat  !Saturated soil water (vol.water/vol.soil)
                                !Equals porosity
      real*8,intent(in) :: thetamin !Hygroscopic H2O cont(vol.water/vol.soil)
      real*8,intent(in) :: fracroot(:) !Fraction of roots in layer
      real*8,intent(in) :: fice(:)  !Fraction of ice in layer
      real*8,intent(in) :: hwilt  !Wilting point of pft, matric pot. (m)
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
        betad = betad +  (1.d0-fice(k))*fracroot(k)*betak
      end do
      if (betad < EPS2) betad=0.d0

      end function water_stress2

!----------------------------------------------------------------------!
      function water_stress3(nlayers, thetas, thetafc, thetamin, 
     &     fracroot, fice, hwilt) Result(betad)

      implicit none
      integer,intent(in) :: nlayers !Number of soil layers
      real*8,intent(in) ::  thetas(:) !Soil vol. water (vol.water/vol.soil)
      real*8,intent(in) :: thetafc  !Field capacity (vol.water/vol.soil)
      real*8,intent(in) :: thetamin !Hygroscopic H2O cont(vol.water/vol.soil)
      real*8,intent(in) :: fracroot(:) !Fraction of roots in layer
      real*8,intent(in) :: fice(:)  !Fraction of ice in layer
      real*8,intent(in) :: hwilt  !Wilting point of pft, matric pot. (m)
      real*8 :: betad !Stress value, 0-1, 1=no stress
      integer :: k
      
      !3. Relative extractable water,REW
      betad = 0.d0
      do k = 1,nlayers
        betad = betad +  (1.d0-fice(k))*fracroot(k)*
     &       (thetas(k)-thetamin)/(thetafc-thetamin)
      end do
      if (betad < EPS2) betad=0.d0

      end function water_stress3


      !************************************************************************
      subroutine qsimp(vpar,ppar,S,sbeta,I0dr,I0df,Ci,T,P)
!----------------------------------------------------------------------!
! qsimp calculates canopy photosynthesis by increasing the number of
! layers in the integration until the result (S) changes by less than
! 0.1 umol[CO2]/m2/s.
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!Passed parameters
      type(veg_par_type) :: vpar
      type(photosynth_par_type) :: ppar
      real*8 :: S,sbeta,I0dr, I0df, Ci,T,P
!----------------------------------------------------------------------!
!Local variables
!nu   real*8, parameter :: EPS=1.D-3
      integer, parameter :: MAXIT=6
      real*8,parameter :: ERRLIM=0.1d0
      integer IT, canopylayers
      real*8 A, B, OST,OS,ST,ERROR
!----------------------------------------------------------------------!
!################## RADIATIVE TRANSFER########################################
! Calculate canopy radiative transfer and photosynthesis in increasing 
!   number of layers.
      A=0.0D0
      B=vpar%alai
      OST=-1.D30
      OS= -1.D30
      canopylayers=1
      do 11 IT=1,MAXIT
         CALL TRAPZD(A,B,ST,IT,
     $       sbeta,I0dr,I0df,Ci,T,P,canopylayers,vpar,ppar)
         S=(4.D0*ST-OST)/3.D0
         ERROR=ABS(S-OS)
         IF (ERROR.lt.ERRLIM) RETURN
         OS=S
         OST=ST
         if(IT.gt.1) canopylayers=canopylayers*2
   11 enddo
!      write(99,*) 'Too many steps.'
!      write(99,*) S,ERROR,100.0D0*ERROR/S
!#############################################################################
      return
      end subroutine qsimp

      !************************************************************************
      subroutine trapzd(A,B,S,N,
     $     sbeta,I0dr,I0df,Ci,T,P,canopylayers,vpar,ppar)
!----------------------------------------------------------------------!
! Integrates canopy photosynthesis over canopy layers using Simpson's
! Rule (Press et al., 19??).
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
      integer :: N,canopylayers
      real*8 :: A,B,S,sbeta,I0dr,I0df,Ci,T,P 
      !nf,pcp,Kc,Oi,Ko,N0,k,n1,n2,m1,msat
      !real*8 sigma,temp,rhor,kdf,kbl
      type(veg_par_type) :: vpar
      type(photosynth_par_type) :: ppar
      !----Local----------------------
      integer :: L
      real*8 :: DEL,X,SUM,RCL
!@var func1 Mean net photosynthesis at Lc (umol[CO2]/m2/s).
      real*8 :: func1
!@var func2 Mean net photosynthesis at Lc (umol[CO2]/m2/s).
      real*8 :: func2

      if(N.eq.1)then
         call phot(A,sbeta,I0dr,I0df,Ci,T,P,B,func1,vpar,ppar)
         call phot(B,sbeta,I0dr,I0df,Ci,T,P,B,func2,vpar,ppar)
         S=0.5D0*(B-A)*(func1+func2)
      else
         RCL=canopylayers  ! convert to real*8
         DEL=(B-A)/RCL
         X=A+0.5D0*DEL
         SUM=0.D0
         do 11 L=1,canopylayers
            call phot(X,sbeta,I0dr,I0df,Ci,T,P,B,func1,vpar,ppar)
           SUM=SUM+func1
           X=X+DEL
   11    continue
         S=0.5D0*(S+(B-A)*SUM/RCL)
      endif
      return
      end subroutine trapzd

      !************************************************************************
      subroutine phot(Lc,sbeta,I0dr,I0df,Ci,T,P,alai,func,vpar,ppar)
!----------------------------------------------------------------------!
! Calculate mean leaf photosynthesis at cumulative leaf area index Lc
! (from top of canopy), returning result as func (umol[CO2]/m2/s).
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
      real*8 :: func,sbeta,I0dr,I0df
      real*8 :: Ci  !(Pa)
      real*8 :: T,P,alai
      type(veg_par_type) :: vpar
      type(photosynth_par_type) :: ppar
!@var Lc Cumulative LAI from top of canopy (m2/m2).
      real*8 Lc
!----------------------------------------------------------------------!
!@var alpha Intrinsic quantum efficiency (?units).
      real*8, parameter :: alpha=0.08D0
!@var ka Chlorophyll extinction coefficient (?units).
      real*8, parameter :: ka=0.005D0
!@var n3 Ratio of foliage chlorophyll to N (?units).
      real*8 n3
!@var Np Foliage nitrogen content (mmol/m2).
      real*8 Np
!@var Isha PAR penetrating shaded foliage at Lc (umol/m[foliage]2/s).
      real*8 Isha
!@var Isla PAR penetrating sunlit foliage at Lc (umol/m[foliage]2/s).
      real*8 Isla
!@var fsl Fraction of sunlit foliage in layer at Lc (unitless).
      real*8 fsl
!@var Nlim Theoretical cumulative N for light-saturated A (mmol/m2).
      real*8 Nlim
!@var Nsat Actual cumulative N for light-saturated A (mmol/m2).
      real*8 Nsat
!@var fabs Absorbed radiation in light-limited chloroplasts (fraction).
      real*8 fabs
!@var FUNCsl Photosynthesis in sunlit foliage (umol/m2/s).
      real*8 FUNCsl
!@var FUNCsh Photosynthesis in shaded foliage (umol/m2/s).
      real*8 FUNCsh
!@var To stop dbzs.
      real*8, parameter :: EPS = 1.d-8
!----------------------------------------------------------------------!
! Check following OK. n3 is ratio of chlorophyll to foliage N (units?).
      n3=6.0D0-3.6D0*exp(-0.7D0*Lc)
! Foliage N at canopy depth Lc (mmol/m2).
      Np=ppar%N0*exp(-ppar%k*Lc)

!############################ RADIATIVE TRANSFER #############################
! Get incident radiation on layer, diffuse/direct, sunlit/shaded.
      call canopy_rad(sbeta, Lc, I0df, I0dr, vpar,Isla, Isha, fsl)
!#############################################################################
!----------------------------------------------------------------------!
! Cumulative nitrogen concentration at which photosynthesis becomes
! light-limited in sunlit foliage (mmol/m2).
      Nlim=-dlog(ppar%msat/(alpha*Isla*ka*n3*ppar%m1+EPS))/(ka*n3)
! Impose limits on Nsat based on actual foliage nitrogen (mmol/m2).
      if(Nlim.lt.zero)then
        Nsat=zero
      elseif(Nlim.gt.Np)then
        Nsat=Np
      else
        Nsat=Nlim
      endif
! Absorbed radiation in light-limited chloroplasts (fraction).
      fabs=exp(-ka*n3*Nsat)-exp(-ka*n3*Np)
! Photosynthesis in sunlit foliage (umol/m2/s).
      FUNCsl=(1.0D0-ppar%PCP/(Ci+EPS))*
     $     (ppar%msat*Nsat+alpha*ppar%m1*Isla*fabs)
!----------------------------------------------------------------------!
! Cumulative nitrogen concentration at which photosynthesis becomes
! light-limited in shaded foliage (mmol/m2).
      Nlim=-log(ppar%msat/(alpha*Isha*ka*n3*ppar%m1+EPS))/(ka*n3)
! Impose limits on Nsat based on actual foliage nitrogen (mmol/m2).
      if(Nlim.lt.0.0D0)then
        Nsat=0.0D0
      elseif(Nlim.gt.Np)then
        Nsat=Np
      else
        Nsat=Nlim
      endif
! Absorbed radiation in light-limited chloroplasts (fraction).
      fabs=exp(-ka*n3*Nsat)-exp(-ka*n3*Np)
! Photosynthesis in shaded foliage (umol/m2/s).
      FUNCsh=(1.0D0-ppar%PCP/(Ci+EPS))*
     &  (ppar%msat*Nsat+alpha*ppar%m1*Isha*fabs)
!----------------------------------------------------------------------!
! Mean photosynthesis in layer (umol/m2/s).
      func=fsl*FUNCsl+(1.0D0-fsl)*FUNCsh
!----------------------------------------------------------------------!
#ifdef DEBUG      !## DEBUG ##
      if (Ci.ne.1.0d6) then
      !write(97,*)"sbeta,Lc,I0df,I0dr,vpar.sigma,sqrtexpr,kdf,rhor,kbl,alai,nm,vh,Ntot,vegalbedo,Isla,Isha,fsl,func",
        write(97,*) sbeta,Lc,I0df,I0dr,vpar,ppar,Isla,Isha,fsl,Ci,
     &       FUNCsl,FUNCsh,func
      end if
#endif
      return
      end subroutine phot


!##################RADIATIVE TRANSFER#########################################
!# THIS WILL BE REPLACED WITH ROUTINES IN ENT canopyradiation.f
!#############################################################################

      subroutine canopy_rad_setup(
     i     sbeta, fdir, parinc_W_m2,vegpar,
     o     I0df, I0dr)

      implicit none

! NOTES FOR GISS GCM ONLY:
! If only incident shortwave is available, convert 
!  total shortwave (W/m2) to PAR (umol/m2/s).
!  2.3 umol/J = shortwave to fraction that is PAR (Monteith & Unsworth).
!      PAR=2.3d0*max(srht,zero)/(1.0D0-0.08D0)
! Current: Replaced back-calculation with actual incident PAR at surface.

! GENERAL PAR CONVERSIONS.
! nyk  For 400-700 nm, energy is 3.3-6 umol photons/J.
!      Top-of-atmosphere flux-weighted average of energy is 4.54 umol/J.
!      Univ. Maryland, Dept. of Met., PAR Project, suggests nominal
!      485 nm for conversion, which gives 4.05 umol/J.

      ! CONSTANTS
      real*8,parameter :: sigma=0.2d0  !Leaf scattering coefficient
      real*8,parameter :: kdf=0.71d0  !Canopy extinction coefficient

      !Input parameters
      real*8 :: sbeta, fdir, parinc_W_m2
      !Output parameters
      real*8 :: I0df, I0dr !(umol/m2/s)
      type(veg_par_type) :: vegpar
      !----Local var---------------------
      real*8 :: PAR

      PAR=4.05d0*parinc_W_m2          !W/m2 to umol/m2/s
! Direct beam PAR incident on canopy (umol/m2/s).
      I0dr=fdir*PAR
! Diffuse PAR incident on canopy (umol/m2/s).
      I0df=(1.0D0-fdir)*PAR
! Canopy reflectivity depends on sun angle. Or take prescribed albedos.
      vegpar%sigma = sigma
      vegpar%kdf = kdf
      vegpar%sqrtexpr=sqrt(1.0D0-vegpar%sigma)
      !Hack canopy reflectivity
      if (vegpar%vegalbedo.eq.0) then !because radiation gets initialized late
        vegpar%rhor=((1.0D0-vegpar%sqrtexpr)/(1.0D0+vegpar%sqrtexpr))*
     &       (2.0D0/(1.0D0+1.6D0*sbeta))
        !write (99,*) 'rhor', rhor
      else
       !Temporary: keep prescribed veg albedos until have canopy scheme.
        vegpar%rhor = vegpar%vegalbedo
        !write (99,*) 'vegalbedo', vegalbedo
      end if
! Canopy extinction coefficient for direct beam radiation depends on
! sun angle (unitless).
      vegpar%kbl=0.5D0*vegpar%kdf/(0.8D0*vegpar%sqrtexpr*sbeta+EPS)
      end subroutine canopy_rad_setup

!-----------------------------------------------------------------------------
      subroutine canopy_rad(
     i     sbeta,               !cos of solar zenith angle
     i     Lc,                  !Cumulative LAI at layer from top of canopy
     i     I0df,                !Incident diffuse PAR, top of layer Lc
     i     I0dr,                !Incident direct PAR, top of layer Lc
     i     vpar,                !veg parameters
     o     Isla,                !PAR on shaded foliage in layer Lc
     o     Isha,                 !PAR on sunlit foliage in layer Lic
     o     fsl)                 !Fraction of sunlit foliage at Lc (unitless).


      implicit none
! Get incident radiation on layer, diffuse/direct, sunlit/shaded.
      !Input parameters
      real*8 :: sbeta, Lc, I0df, I0dr
      type(veg_par_type) :: vpar
      !Output var parameters
!@var Isla PAR penetrating sunlit foliage at Lc (umol/m[foliage]2/s).
      real*8 :: Isla
!@var Isha PAR penetrating shaded foliage at Lc (umol/m[foliage]2/s).
      real*8 :: Isha
      !Local vars------------------------------------------------------
!@var Idfa Diffuse PAR in canopy at Lc (umol/m[ground]2/s).
      real*8 Idfa
!@var Idra Direct PAR in canopy at Lc (umol/m[ground]2/s).
      real*8 Idra
!@var Idrdra Direct PAR at Lc that remains direct (umol/m[ground]2/s).
      real*8 Idrdra
!@var fsl Fraction of sunlit foliage in layer at Lc (unitless).
      real*8 fsl

! When sun is above the horizon.
      if(sbeta.gt.zero)then
! Diffuse PAR in canopy at Lc (umol/m[ground]2/s).
        Idfa=(1.0D0-vpar%rhor)*I0df* vpar%kdf *exp(-vpar%kdf*Lc)
! Direct PAR in canopy at Lc (umol/m[ground]2/s).
        Idra=(1.0D0-vpar%rhor)*I0dr* vpar%sqrtexpr * vpar%kbl *
     $        exp(-vpar%sqrtexpr*vpar%kbl*Lc)
! Direct PAR at Lc that remains direct (umol/m[ground]2/s).
        Idrdra=(1.0D0-vpar%sigma)*I0dr * vpar%kbl *
     $        exp(-vpar%sqrtexpr*vpar%kbl*Lc)
! PAR penetrating shaded foliage at Lc (umol/m[foliage]2/s).
        Isha=Idfa+(Idra-Idrdra)
! PAR penetrating sunlit foliage at Lc (umol/m[foliage]2/s)._
        Isla=Isha+(1.0D0-vpar%sigma)* vpar%kbL * I0dr
! Fraction of sunlit foliage in layer at Lc (unitless).
        fsl=exp(-vpar%kbL*Lc)
        if(Isha.lt.0.1D0)Isha=0.1D0
        if(Isla.lt.0.1D0)Isla=0.1D0
        if(fsl.lt.zero)fsl=zero 
      else
        Isha=0.1D0
        Isla=0.1D0
        fsl=zero
      end if
!      write(99,*) I0dr, I0df,Isha, Isla, kbl

      end subroutine canopy_rad

!-----------------------------------------------------------------------------

      subroutine canopy_transmittance(TRANS_SW, sbeta,fdir,vp)
! Calculates the transmittance (fraction) of radiation through the canopy
! to the soil surface.      
! Transmission of shortwave radiation through canopy. This has errors.
!----------------------------------------------------------------------!
      implicit none
!----------------------------------------------------------------------!
!Parameters -------------------------------------------------------------
!@var trans_sw Total canopy transmittance of shortwave (fraction)
      real*8, intent(out) :: TRANS_SW  
!@var sbeta Cos of solar zenith angle
      real*8, intent(in) :: sbeta
!@var fdir Fraction of surface visible radiation that is direct
      real*8, intent(in) :: fdir
! Vegetation type cell parameters
      type(veg_par_type) :: vp

!Local variables --------------------------------------------------------
! Get assigned from vp passed in.
      real*8 sigma,sqrtexpr,rhor,kbl,kdf,alai
!@var absdf Absorbance of diffuse SW radiation (fraction)
      real*8 absdf
!@var absdr Absorbance of direct SW radiation (fraction)
      real*8 absdr
!@var absdrdr Absorbance of direct SW that remains direct (fraction)
      real*8 absdrdr
!@var abssh Absorbance of SW through shaded foliage (fraction)
      real*8 abssh
!@var abssl Absorbance of SW through sunlit foliage (fraction)
      real*8 abssl
!@var fracsl Fraction of leaves in canopy that are sunlit (fraction)
      real*8 fracsl
!------------------------------------------------------------------------
      sigma=vp%sigma
      sqrtexpr=vp%sqrtexpr
      rhor=vp%rhor
      kbl=vp%kbl
      kdf=vp%kdf
      alai=vp%alai

! Diffuse shortwave in canopy (umol/m[ground]2/s).
      if(sbeta.gt.zero)then
        absdf=(1-fdir)*(1.0D0-rhor)*kdf*exp(-kdf*alai)
! Direct shortwave in canopy (umol/m[ground]2/s).
        absdr=fdir*(1.0D0-rhor)*sqrtexpr*kbl*exp(-sqrtexpr*kbl*alai)
! Direct shortwave that remains direct (umol/m[ground]2/s).
        absdrdr=absdr*(1.0D0-sigma)*kbl*exp(-sqrtexpr*kbl*alai)
! Shortwave penetrating shaded foliage (umol/m[foliage]2/s).
        abssh=absdf + (absdr-absdrdr) 
! Shortwave penetrating sunlit foliage (umol/m[foliage]2/s).
        abssl=abssh+(1.0D0-sigma)*kbL*fdir
        if(abssh.lt.0.0D0) abssh=0.0D0
        if(abssl.lt.0.0D0) abssl=0.0D0
        if(abssh.gt.1.0D0) abssh=1.0D0
        if(abssl.gt.0.0D0) abssl=1.0D0
        fracsl=exp(-kbl*alai)
      else
        abssh=0.001D0
        abssl=0.0D0
        fracsl=0
      endif
!      TRANS_SW = 1-((1-fracsl)*abssh + fracsl*abssl)
      TRANS_SW = ((1.d0-fracsl)*abssh + fracsl*abssl)
!      write(110,*) TRANS_SW,I0dr,I0df,abssh, abssl, sbeta, sigma, kbl


      end subroutine canopy_transmittance

!##############################################################################

      FUNCTION QSAT (TM,LH,PR)
!@sum  QSAT calculates saturation vapour mixing ratio
!@auth Gary Russell
!@ver  1.0
!      USE CONSTANT, only : mrat,rvap,tf
      IMPLICIT NONE
!@var Physical constants from GISS GCM CONST.f
      real*8, parameter :: MWAT = 18.015d0 !molecular weight of water vapour (g/mol)
      real*8, parameter :: MAIR = 28.9655d0 !molecular weight of dry air (28.9655 g/mol)
      real*8, parameter :: MRAT = MWAT/MAIR 
      real*8, parameter :: RVAP = 1d3 * GASC/MWAT !gas constant for water vapour (461.5 J/K kg)
!@var A,B,C   expansion coefficients for QSAT
      REAL*8, PARAMETER :: A=6.108d0*MRAT    !3.797915d0
      REAL*8, PARAMETER :: B= 1./(RVAP*TF)   !7.93252d-6
      REAL*8, PARAMETER :: C= 1./RVAP        !2.166847d-3
C**** Note that if LH is considered to be a function of temperature, the
C**** correct argument in QSAT is the average LH from t=0 (C) to TM, ie.
C**** LH = 0.5*(LH(0)+LH(t))
      REAL*8, INTENT(IN) :: TM  !@var TM   temperature (K)
      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
      REAL*8, INTENT(IN) :: LH  !@var LH   lat. heat of vap./sub. (J/kg)
      REAL*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
      QSAT = A*EXP(LH*(B-C/max(130.d0,TM)))/PR
      RETURN
      END FUNCTION QSAT

!======================================================================!

!     real*8 FUNCTION QSATold (TM,QL,PR)
!      implicit none
!!@sum  QSAT calculates saturation vapour mixing ratio (kg/kg)
!!@auth Gary Russell
!!@ver  1.0
!!      USE CONSTANT, only : mrat,rvap,tf
!!      IMPLICIT NONE
!!@var A,B,C   expansion coefficients for QSAT
!      REAL*8, PARAMETER :: A=3.797915e0    !3.797915d0
!      REAL*8, PARAMETER :: B=7.93252e-6    !7.93252d-6
!      REAL*8, PARAMETER :: C=2.166847e-3         !2.166847d-3
!      real*8 :: TM, QL, PR
!!**** Note that if QL is considered to be a function of temperature, the
!!**** correct argument in QSAT is the average QL from t=0 (C) to TM, ie.
!!**** QL = 0.5*(QL(0)+QL(t))
!!      REAL*8, INTENT(IN) :: TM  !@var TM   potential temperature (K)
!!      REAL*8, INTENT(IN) :: PR  !@var PR   air pressure (mb)
!!     REAL*8, INTENT(IN) :: QL  !@var QL   lat. heat of vap./sub. (J/kg)
!!      REAL*8 :: QSAT            !@var QSAT sat. vapour mixing ratio
!      QSAT = A*EXP(QL*(B-C/TM))/PR
!      RETURN
!      END FUNCTION QSATold
!======================================================================!


!**   OPTIONAL SUBROUTINES **!

      subroutine update_veg_locals(evap_tot, rho, ch, U,Qsurf, Qf)
      !Update vegetation input variables that require values external
      !to vegetation module to update.  For values that change within
      !the smaller time step of the GHY modules.
      !---------------------------------------------------------------!
      implicit none
      !---------------------------------------------------------------!
      real*8,parameter:: rhow=1d3 !CONSTANT, density of water (1000 kg/m^3)

      real*8,intent(in):: evap_tot !Total evaporation over veg+soil (m/s)
      real*8,intent(in):: rho   !Density of air at the surface (kg/m^3)
      real*8,intent(in):: ch  !Ground to surface heat transfer coefficient
      real*8,intent(in):: U !Surface wind speed (m/s)
      real*8,intent(in):: Qsurf !Surface humidity (kg/kg)
      real*8,intent(out) :: Qf !Foliage surface H2O vapor mixing ratio (kg[H2O]/kg[air])
      
      !Get new Qf
      Qf=evap_tot*rhow/(rho*ch*U)+Qsurf ! New mixing ratio for next timestep
      end subroutine update_veg_locals


      !----------------------------------------------------------------


      subroutine veg_accum
      !Accumulate GPP, R_growth, R_maint, R_root, Soil_resp
      !---------------------------------------------------------------!
      implicit none
      !---------------------------------------------------------------!
      ! agpp = agpp + gpp*(1.d0-fr_snow(2)*fm)*fv*dts

      entry veg_accm0
      
      !agpp=0.d0   ! new accumulator, nyk 4/25/03
      

      end subroutine veg_accum



      !************************************************************************
      
      end module biophysics
