#include "rundeck_opts.h"
      module phenology
!@sum Routines to calculate phenological change in an entcell:
!@sum budburst/leafout, albedo change, senescence
!@auth Y. Kim
#ifdef ENT_STANDALONE_DIAG
#define PHENOLOGY_DIAG
#define DEBUG
#endif

      use ent_types
      use ent_const
      use ent_pfts
      use allometryfn

      implicit none
!      public veg_init
      public clim_stats
      public pheno_update
      public veg_update !may change the name into veg_update 
      public litter_cohort, litter_patch   !Now called from veg_update

      private pheno_update_coldwoody
      private pheno_update_coldherbaceous
      private pheno_update_drought
      private growth_cpools_active
      private growth_cpools_structural
!      private senesce_cpools
      private recruit_annual
      private photosyn_acclim
      private phenology_diag

      !**********************************************************************

      !** GROWTH MODEL CONSTANTS - phenology & carbon allocation 
      !*l_fract: fraction of leaf C retained after leaf fall (unitless) (value from ED) 
      real*8, parameter :: l_fract = 0.50d0 
      !*q: ratio of root to leaf biomass (unitless) (value from ED)
      !real*8, parameter :: q=1.0d0  !Moved to allometryfn.f
      !*iqsw: sapwood biomass per (leaf area x wood height) (kgC/m2/m) 
      !3900.0: leaf area per sapwood area (m2/m2) (value from ED)
      !1000.0: sapwood density (kg/m3)
      !2.0:  biomass per carbon (kg/kgC)
      !(qsw)=(iqsw*sla) (1/m) & (qsw*h): ratio of sapwood to leaf biomass (unitless)
      !(iqsw)=1000.0d0/3900.0d0/2.0d0=0.1282 !NOTE: This value corrects an error in the coefficient in Moorcroft et al. (2001) Appendix D, which had the value too small by a factor of 100, at 0.00128.
      !real*8, parameter :: iqsw=1000.0d0/3900.0d0/2.0d0 !Moved to allometryfn.f

      !*hw_fract: ratio of above ground stem to total stem (stem plus structural roots) (value from ED)
      !real*8, parameter :: hw_fract = 0.70d0 !Moved to allometryfn.f
      !*C2B: ratio of biomass to carbon (kg-Biomass/kg-Carbon) 
      real*8, parameter :: C2B = 2.0d0 
      !*temperature constraint for cold-deciduous PFTs (Botta et al. 1997)
      !*airtemp_par !base temperature to calculate the growing degree days (gdd)
      !*gdd_par1/2/3: paramters to estimate the threshold for gdd
      !*gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)    
      real*8, parameter :: airtemp_par = 5.d0 
      real*8, parameter :: gdd_par1 = -68.d0 
      real*8, parameter :: gdd_par2 = 638.d0
      real*8, parameter :: gdd_par3 = -0.01d0 
      !*gdd_length - tuning parameters (tuned for HF & MMSF)
      real*8, parameter :: gdd_length = 200.d0 
      !*airt_threshold - tuning parameters (tuned for HF & MMSF)
      real*8, parameter :: airt_max_w = 15.d0
      real*8, parameter :: airt_min_w = 5.d0
      !*soilt_threshold - tuning parameters 
      real*8, parameter :: soilt_max = 5.d0  !10.d0 - old repository value (tuned for Barrow)
      real*8, parameter :: soilt_min = 0.d0
      real*8, parameter :: soilt_base = -5.d0 !Added from YK
      !*sgdd_threshold & length - tuning parameters (tuned for Barrow)
      real*8, parameter :: soiltemp_par = 0.d0
      real*8, parameter :: sgdd_threshold = 100.d0
      real*8, parameter :: sgdd_length=50.d0   
      !*ld_threshold (minute): light length constraint for cold-deciduous woody PFTs (White et al. 1997)
      real*8, parameter :: ld_threshold = 655.d0
      real*8, parameter :: ld_min =540.d0
      real*8, parameter :: ld_max =550.d0
      !*tsoil_threshold1, tsoil_threshold2 : soil temperature constraint for cold-deciduous woody PFTs (White et al. 1997)
!      real*8, parameter :: tsoil_threshold1 = 11.15d0
!      real*8, parameter :: tsoil_threshold2 = 2.d0
      !*ddfacu: the rate of leaf fall (1/day) (value from IBIS)
      real*8, parameter :: ddfacu = 1.d0/15.d0
      !*betad : water_stress3  - tunning parameters (sstar/swilt determines betad, then check those first and tune these)
      !_w for woody & _h  (tunned for MMSF) - max/min/resistance parameters for woody
      real*8, parameter :: betad_max_w = 0.1d0
      real*8, parameter :: betad_min_w = 0.d0 
      real*8, parameter :: betad_res_w = 0.25d0
      ! _h for herbaceous  (tunned for Vaira/Tonzi) - max/min/resistance prameters for woody
      real*8, parameter :: betad_max_h = 0.9d0 
      real*8, parameter :: betad_min_h = 0.4d0
      real*8, parameter :: betad_res_h = 1.0d0      
      !*light-controll phenology model (Kim et al; originally for ED2)
      !*different from the original implementation, as it cannot be directly implemented due to model difference.
      !*e.g.) PAR instead of Rshort & other differences in parameterization requires the model to be tunned for Ent.
      !*par_turnover_int & par_turnover_slope  - tunning parameters (tunned for TNF; not finalized)
      real*8, parameter :: par_turnover_int = -12.d0 
      real*8, parameter :: par_turnover_slope = 0.18d0  
      !*r_fract: fraction of excess c going to seed reproduction (value from ED)
      real*8, parameter :: r_fract = 0.3d0
      !*c_fract: fraction of excess c going to clonal reproduction - only for herbaceous (value from ED)
      real*8, parameter :: c_fract = 0.7d0
      !*mort_seedling: mortality rate for seedling (value from ED)
      real*8, parameter :: mort_seedling = 0.90d0 

      contains


      !*********************************************************************
      subroutine clim_stats(dtsec, ecp, config,dailyupdate)
!@sum Calculate climate statistics such as 10 day running mean 
!@+   Called by ent.f every physical time step.
      use soilbgc, only : Soillayer_convert_Ent 
      real*8,intent(in) :: dtsec           !dt in seconds
      type(entcelltype) :: ecp      
      type(ent_config) :: config
      logical, intent(in) :: dailyupdate  
      !-----local--------
      type(patch), pointer :: pp
      type(cohort), pointer :: cop 
      !*local variables for entcell-level envrionment variable
      real*8 :: airtemp        !air temperature degC 
      real*8 :: soiltemp       !soil temperature degC
      real*8 :: par            !photosynthetic active radiation (PAR)
      !*local variables for entcell-level variables, 
      !*updated in this subroutine
      real*8 :: airtemp_10d    !10 day running mean of air temperature
      real*8 :: soiltemp_10d   !10 day running mean of soil temperature
      real*8 :: par_10d        !10 day running mean of PAR
      real*8 :: gdd            !growing degree days, based on air temperature
      real*8 :: ncd            !number of chilling days, based on air temperature
      real*8 :: sgdd           !growing degree days, based on soil temperature
      !*PAR-limited phenology parameters
      real*8 :: par_crit       !PAR threshold
      logical :: par_limit     !logical whether PAR-limited phenology parameterization is applied or not for certain PFTs
      real*8 :: turnover0      !turnover amplitude, calculated with the phenology parameterization
      real*8 ::  llspan0       !leaf life span, calculated with the phenology parameterization
      !*soil temperature for CASA layers
      real*8 :: Soiltemp2layer(N_CASA_LAYERS)  
     
  
      airtemp = ecp%TairC
      call Soillayer_convert_Ent(ecp%Soiltemp(:), SOILDEPTH_m, 
     &     Soiltemp2layer) 
      soiltemp = Soiltemp2layer(1)

      soiltemp_10d = ecp%soiltemp_10d
      airtemp_10d = ecp%airtemp_10d
      par_10d = ecp%par_10d  
      gdd = ecp%gdd
      ncd = ecp%ncd
      sgdd = ecp%sgdd

      !*10-day running mean of Air Temperature
      airtemp_10d = running_mean(dtsec, 10.d0, airtemp, airtemp_10d)

      !*10-day running mean of Soil Temperature
      soiltemp_10d = running_mean(dtsec, 10.d0, soiltemp, soiltemp_10d)

      !*10-day running mean of PAR
      par = ecp%IPARdif + ecp%IPARdir !total PAR is the sum of diffused and direct PARs
      par_10d = running_mean(dtsec, 10.d0, par, par_10d)

      !*daylength
      if (ecp%CosZen > 0.d0) then
         ecp%daylength(2) = ecp%daylength(2) + dtsec/60.d0
      end if

      !*GDD & NCD - Update Once a day 
      if (dailyupdate) then
         !*Calculate Growing degree days
         if (airtemp_10d .ge. airtemp_par) then
            gdd = gdd + ( airtemp_10d - airtemp_par )
         end if
         !*Calculate Growing degree days for soil temperature 
         if (soiltemp_10d .ge. soiltemp_par) then
            sgdd = sgdd + ( soiltemp_10d -soiltemp_par )
         end if
         !*Number of chilling days 
         !number of days below airtemp_par - chilling requirements
         if (airtemp_10d .lt. airtemp_par) then
            ncd = ncd +  1.d0
         end if
         !*If the season is fall or not (if fall, it's 1; else, it's 0)
         !1) it is to control the phenological status      
         !2) it is determined according to whether the daylength is decreasing (i.e., fall) or not. 
         if (NInt(ecp%daylength(2)) .lt. NInt(ecp%daylength(1)) ) then
            ecp%fall = 1
         else if (NInt(ecp%daylength(2)).gt.NInt(ecp%daylength(1))) then
            ecp%fall = 0
         end if 
       end if

      pp => ecp%oldest 
      do while (ASSOCIATED(pp)) 

        cop => pp%tallest
        do while(ASSOCIATED(cop))
          
          !*10-day running mean of stressH2O (betad)
          cop%betad_10d = running_mean(dtsec, 10.d0, 
     &                    cop%stressH2O, cop%betad_10d)
     &                     
          !*Daily carbon balance
          !it is used for the carbon allocation
          cop%CB_d =  cop%CB_d + cop%NPP*dtsec/cop%n*1000.d0

          !*********************************************
          !* evergreen broadleaf - PAR limited - not finalized yet!
          !********************************************* 

          !if it is evergreen & broadleaf, radiation-limited phenology is working    
          par_limit = ((pfpar(cop%pft)%phenotype.eq.EVERGREEN).and.  
     &                (pfpar(cop%pft)%leaftype.eq.BROADLEAF))
          !par_limit = .false. !temp. suppress
          !raidation-limited phenology model 
          if (par_limit) then
             par_crit = - par_turnover_int/par_turnover_slope 

             !calculate the turnover amplitude 
             !(relative ratio of turnover compared to its intrinsic turnover rate)
             !based on PAR
             turnover0 = min(100.d0, max(0.01d0, 
     &          par_turnover_slope*par_10d + par_turnover_int))

             if (par_10d .lt. par_crit) turnover0 = 0.01d0

             !calculate 10 day running mean of turnover amplitude
             cop%turnover_amp = running_mean(dtsec,10.d0, 
     &                          turnover0, cop%turnover_amp)

             !calculate the leaf life span based on turnover amplitude
             !lrage is in year, llspan is in month, and then 12 is used to convert the units
             llspan0 = pfpar(cop%pft)%lrage*12.d0/cop%turnover_amp 

             !calculate 90 day running mean of llspan
             cop%llspan = running_mean(dtsec, 90.d0,llspan0,cop%llspan)

          else

             cop%turnover_amp = 1.d0 
             cop%llspan = undef !-999.d0

          endif
          

          !**************************************************************
          !* Update photosynthetic acclimation factor for evergreen veg
          !**************************************************************
          if (config%do_frost_hardiness) then 
             if (((pfpar(cop%pft)%phenotype.eq.EVERGREEN).and.  
     &           (pfpar(cop%pft)%leaftype.eq.NEEDLELEAF)).or.
     &           (pfpar(cop%pft)%phenotype.eq.COLDDECID).or.
     &           (pfpar(cop%pft)%phenotype.eq.COLDDROUGHTDECID)) then
                call photosyn_acclim(dtsec,airtemp_10d,cop%Sacclim) 
             else
                cop%Sacclim = 25.d0 !Force no cold hardening, mild temperature.
             endif
	  else
              cop%Sacclim = 25.d0 !Force no cold hardening, mild temperature.
          endif

          cop => cop%shorter  
        end do        
        pp => pp%younger 
      end do 

      ecp%soiltemp_10d = soiltemp_10d
      ecp%airtemp_10d = airtemp_10d
      ecp%par_10d = par_10d
      ecp%gdd = gdd
      ecp%ncd = ncd
      ecp%sgdd = sgdd

      end subroutine clim_stats

      !*********************************************************************   
      subroutine photosyn_acclim(dtsec,Ta,Sacc)
!@sum Calculate Sacclim accumulator state used to calculate 
!@+   acclimation/frost hardiness for boreal coniferous forests.  
!@+   Called by clim_stats during update of phenology climate stats.
!@+   Based on Repo et al (1990), 
!@+   Hanninen & Kramer (2007),and  Makela et al (2006)
      implicit none
      real*8,intent(in) :: dtsec ! time step size [sec]
      real*8,intent(in) :: Ta ! air temperature [deg C]
      real*8,intent(inout) :: Sacc ! state of acclimation [deg C]

      !----Local-----
      real*8,parameter :: tau_inv = 2.22222e-6 
                          ! inverse of time constant of delayed
                          ! response to ambient temperature [sec] = 125 hr 
                          ! Makela et al (2004) for Scots pine

!      Use a first-order Euler scheme
       Sacc = Sacc + dtsec*(tau_inv)*(Ta - Sacc) 

!     Predictor-corrector method requires temperature from next timestep
!       Sacc_old = Sacc
!       Sacc = Sacc_old + dtsec*(1/tau_acclim)*(Ta - Sacc ) 
!       Sacc = Sacc_old + ((1/tau_acclim)*(Ta - Sacc_old)+
!     &                     (1/tau_acclim)*(Ta_next - Sacc))*0.5d*dtsec

      end subroutine photosyn_acclim


      !*********************************************************************   
      subroutine pheno_update(pp)
!@sum Update statstics for phneology_update    
!@sum DAILY time step.
      use ent_const

      type(patch) :: pp
      !--Local-----
      type(cohort), pointer :: cop
      integer :: pft
      integer :: phenotype
      real*8 :: airtemp_10d    !10 day running mean of air temperature
      real*8 :: soiltemp_10d   !10 day running mean of soil temperature    
      real*8::  betad_10d      !10 day running mean of betad (calculated with stressH2O)  
      real*8 :: ld             !day length in minutes
      real*8 :: gdd            !growing degree days, based on air temperature
      real*8 :: ncd            !number of chilling days, based on air temperature
      real*8 :: sgdd           !growing degree days, based on soil temperature 
      real*8 :: airt_adj       !adjustment for air temperature threshold (airt_max & airt_min)
      real*8 :: soilt_adj      !adjustment for soil temperature threshold (soilt_max & soilt_min)
      !*phenofactor : phenological elongation factor, ranging 0 (no leaf) to 1 (full leaf)
      !_c for the cold-deciduous & _d for the drought-deciduous
      real*8 :: phenofactor
      real*8 :: phenofactor_c
      real*8 :: phenofactor_d
      !*phenostatus :  define phenological status
      !1 - no leaf (phenofactor, equal to 0; after leaf-off in the  fall, until leaf green-up in the next spring)
      !2 - growing leaf (phenofactor, increasing from 0 to 1 in the spring) 
      !3 - leaf in full growth (phenofactor, equal to 1)
      !4 - leaf senescence (phenofactor, decreasing from 1 to 0 in the fall)
      real*8 :: phenostatus
      real*8 :: phenostatus_c
      real*8 :: phenostatus_d
      !*phenogy type 
      logical :: cold_limit    !.true.=cold deciduous
      logical :: drought_limit  !.true.=drought deciduous
      !*GDD treshold (White et al.)      
      real*8 :: gdd_threshold
      !*whther the season is fall or not (determined in clim_stats)
      logical :: fall         !.true. for fall
      !*whether  PFT is wood or not
      logical :: woody
      !*X day to mature: mature=1+X/1000 
      !1.1 means 100 days to mature.
      !it is devised to prevent abrupt grass green-up 
      !in the middle of winter due to a couple of warm days. 
      real*8 :: mature = 1.1d0 
 
      soiltemp_10d = pp%cellptr%soiltemp_10d
      airtemp_10d = pp%cellptr%airtemp_10d
      gdd = pp%cellptr%gdd
      ncd = pp%cellptr%ncd
      if ( pp%cellptr%fall == 1 ) then
        fall = .true.
      else
        fall = .false.
      endif
      sgdd = pp%cellptr%sgdd
      ld = pp%cellptr%daylength(2)

      cop => pp%tallest
      do while(ASSOCIATED(cop))                
         phenofactor_c=cop%phenofactor_c
         phenofactor_d=cop%phenofactor_d
         phenofactor=cop%phenofactor
         phenostatus_c=cop%phenostatus
         phenostatus_d=cop%phenostatus
         betad_10d=cop%betad_10d
         pft=cop%pft
         phenotype=pfpar(pft)%phenotype
         woody = pfpar(pft)%woody

         !***********************************************
         !*Determine whther PFT cold or drought deciduous
         !***********************************************
         if (phenotype .eq. COLDDECID) then 
            cold_limit = .true.
            drought_limit = .false.
         else if (phenotype .eq. DROUGHTDECID) then 
            cold_limit = .true.
            drought_limit = .true.
         else if (phenotype .eq. EVERGREEN) then
            cold_limit = .false.
            drought_limit = .false.
         else !any of cold and drought deciduous
            cold_limit = .true.
            drought_limit = .true.            
         end if

         !*Set the air/temperature adjustment for specific PFTs
         airt_adj=0.d0
         if (pft.eq.DROUGHTDECIDBROAD)airt_adj=10.d0
         soilt_adj=0.d0
         if (pft .eq. GRASSC3ARCTIC)soilt_adj=-5.d0      


         !*******************************************
         !*Update the phenology for Cold-deciduous
         !*******************************************
         if (cold_limit)then

           !*Cold-deciduous Woody
           if (woody) then  
             !*Update phenofactor and phenostatus
             call pheno_update_coldwoody(cop%phenostatus, 
     i            fall, airtemp_10d, airt_adj, ld,
     o            gdd, ncd, 
     o            phenofactor_c, phenostatus_c)

          !*Cold-decidous Herbaceous
          else 
             !*Update phenofactor and phenostatus
             call pheno_update_coldherbaceous(cop%phenostatus,
     i            fall, soiltemp_10d, soilt_adj, 
     o            sgdd, 
     o            phenofactor_c, phenostatus_c)
           end if

         else !cold_limit=.false.
           phenofactor_c = 1.d0
           phenostatus_c  = 3.d0
         end if          
                 
         !*******************************************
         !*Update the phenology for Drought-deciduous
         !*******************************************
         if (drought_limit)then

           !*Drought-deciduous Woody
           if (woody) then 
             !*Update phenofactor and phenostatus
             call pheno_update_drought(cop%phenostatus, 
     i            mature, 
     i            betad_10d, betad_min_w, betad_max_w, betad_res_w, 
     o            phenofactor_d, phenostatus_d)

         !*Drought-decidous Herbaceous
          else 
             !*Update phenofactor and phenostatus
             call pheno_update_drought(cop%phenostatus, 
     i            mature, 
     i            betad_10d, betad_min_h, betad_max_h, betad_res_h, 
     o            phenofactor_d, phenostatus_d)
           end if   
    
         else !drought_limit=.false.
           phenofactor_d = 1.d0
           phenostatus_d = 3.d0
         end if  

         !*******************************************
         !*Update the phenology according to PFTs 
         !*******************************************
         if (phenotype .eq. COLDDECID) then 
            phenofactor = phenofactor_c
            phenostatus = phenostatus_c
         else if (phenotype .eq. EVERGREEN) then !leaf in full growth
            phenofactor = 1.d0 
            phenostatus = 3.d0
         else if (phenotype .eq. DROUGHTDECIDBROAD .and. 
     &            .not.phenostatus_c.lt.3.d0) then
             phenofactor = phenofactor_c
             phenostatus = phenostatus_c
         else !any of cold and drought deciduous
            phenofactor = phenofactor_c * phenofactor_d   
           if((phenostatus_c.ge.4.d0.and.phenostatus_d.ge.2.d0).or.
     &       (phenostatus_d.ge.4.d0.and.phenostatus_c.ge.2.d0))then
             phenostatus = max(phenostatus_c, phenostatus_d)
           else
             phenostatus = min(phenostatus_c, phenostatus_d)
           end if
          end if
         
         !*increment phenostatus by 0.001 
         !to track how many days after phenostatus has been changed
         if (aint(cop%phenostatus).eq.aint(phenostatus))then
            phenostatus = cop%phenostatus + 1.d0/1000.d0  
         end if
    
#ifdef DEBUG
!             write(202,'(3(i5),100(1pe16.8))') pp%cellptr%fall
             write(202,*) pp%cellptr%fall
     &      ,phenofactor
     &      ,phenofactor_c,phenofactor_d
     &      ,phenostatus, cop%phenostatus
     &      ,phenostatus_c, phenostatus_d
     &      ,betad_10d,soiltemp_10d, mature
     &      ,gdd, ncd
#endif
         
         cop%phenofactor_c=phenofactor_c
         cop%phenofactor_d=phenofactor_d
         cop%phenofactor=phenofactor
         cop%phenostatus=phenostatus
   
         cop => cop%shorter 
      
      end do   

      pp%cellptr%gdd = gdd
      pp%cellptr%ncd = ncd
      pp%cellptr%sgdd = sgdd
      

      end subroutine pheno_update
      !*********************************************************************  
      subroutine pheno_update_coldwoody(phenostatus, 
     i            fall, airtemp_10d, airt_adj, ld,
     o            gdd, ncd, 
     o            phenofactor_c, phenostatus_c)
!@sum Update phenology for cold-decidous woody PFTs
!@sum Called from pheno_update

      use ent_const

      !input variables
      real*8, intent(in) :: phenostatus !phenological status (refer pheno_update for details) 
      logical,intent(in) :: fall        !.true. if the season is fall
      real*8, intent(in) :: airtemp_10d !10 day running mean of air temperature
      real*8, intent(in) :: airt_adj    !adjustment for air temperature threshold (airt_max & airt_min)
      real*8, intent(in) :: ld          !day length in minutes
      !in/output variables
      real*8, intent(inout) :: gdd      !growing degree days, based on air temperature
      real*8, intent(inout) :: ncd      !number of chilling days, based on air temperature  
      !output variables
      real*8, intent(out) :: phenofactor_c !phenological factor for cold deciduous (refer pheno_update for details) 
      real*8, intent(out) :: phenostatus_c !phenological status for cold deciduous (refer pheno_update for details)
      !local variables
      real*8 :: gdd_threshold   !GDD treshold (White et al.)      

      !*GDD threshold for leaf green-up 
      gdd_threshold = gdd_par1 + gdd_par2*exp(gdd_par3*ncd)  

      !*Leaf-on in the spring, triggered by thermal sum
      !if gdd is larger than its threshold 
      !in the spring (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if ((.not. fall) .and.
     &   (phenostatus.lt.3.d0).and.(gdd.gt.gdd_threshold))then  
         !determine phenofactor by scaling gdd with gdd_threshold and gdd_length
         phenofactor_c = min (1.d0,(gdd-gdd_threshold)/gdd_length)
         if (phenofactor_c .lt. 1.d0) then
            phenostatus_c = 2.d0 !growing leaf
         else 
            phenostatus_c = 3.d0 !leaf in full grwoth
         end if
      end if

      !*Leaf-off in the fall 
      !*Leaf-off triggered by air temperature
      !if air temperature is falling below its maximum 
      !in the fall (when it's full-leaf (phenostatus=3.X) or leaf is senescening (phenostatus=4.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if (fall .and. 
     &   (phenostatus.ge.3.d0).and.
     &   (airtemp_10d.lt.airt_max_w+airt_adj)) then
         !determine phenofactor by scaling air temperature
         !with its minimum(min+adj) and maximum(max+adj).
         phenofactor_c = min(phenofactor_c,max(0.d0,
     &      (airtemp_10d-airt_min_w-airt_adj)/
     &      (airt_max_w-airt_min_w)))
         if (phenofactor_c .eq. 0.d0) then
            phenostatus_c = 1.d0   !no leaf
            ncd = 0.d0             !zero-out ncd once complete leaf-off occurs.
            gdd = 0.d0             !zero-out gdd once complete leaf-off occurs.
         else  
            phenostatus_c = 4.d0   !leaf senescence
         end if 
      end if
      !*Leaf-off triggered by day-length
      !if day length is falling shorter than its maximum 
      !in the fall (when it's full-leaf (phenostatus=3.X) or leaf is senescening (phenostatus=4.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if (fall .and.
     &   (phenostatus.ge.3.d0).and.(ld.lt.ld_max)) then
         !dtermine phenofactor by scaling the day length with its min and max.
         phenofactor_c = min(phenofactor_c, max(0.d0,
     &       (ld - ld_min)/(ld_max-ld_min)))
         if (phenofactor_c .eq. 0.0d0) then
            phenostatus_c =1.d0    !no leaf
            ncd = 0.d0             !zero-out ncd once complete leaf-off occurs.
            gdd = 0.d0             !zero-out gdd once complete leaf-off occurs.
         else
            phenostatus_c = 4.d0   !leaf senescence
         end if
      end if    
       
      end subroutine pheno_update_coldwoody
      !*********************************************************************  
      subroutine pheno_update_coldherbaceous(phenostatus, 
     i            fall, soiltemp_10d, soilt_adj,
     o            sgdd,  
     o            phenofactor_c, phenostatus_c)
!@sum Update phenology for cold-decidous herbaceous PFTs
!@sum Called from pheno_update

      use ent_const

      !input variables
      real*8, intent(in) :: phenostatus !phenological status (refer pheno_update for details) 
      logical,intent(in) :: fall         !.true. if the season is fall
      real*8, intent(in) :: soiltemp_10d !10 day running mean of soil temperature
      real*8, intent(in) :: soilt_adj    !adjustment for soil temperature threshold (soilt_max & soilt_min)
      !in/output variables
      real*8, intent(inout) :: sgdd      !growing degree days, based on soil temperature
      !output variables
      real*8, intent(out) :: phenofactor_c !phenological factor for cold deciduous (refer pheno_update for details) 
      real*8, intent(out) :: phenostatus_c !phenological status for cold deciduous (refer pheno_update for details)
        
      !*Leaf-on in the spring, triggered by thermal sum, based on soil temperature
      !if sgdd is larger than its threshold 
      !in the spring (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if ((.not. fall) .and.
     &   (phenostatus.lt.3.d0) .and.
     &   (sgdd.gt.sgdd_threshold)) then
        !determine phenofactor by scaling sgdd with gsdd_threshold and sgdd_length
         phenofactor_c  
     &      = min (1.d0,(sgdd-sgdd_threshold)/sgdd_length)
         if (phenofactor_c .lt. 1.d0) then
            phenostatus_c = 2.d0    !growing leaf
         else 
            phenostatus_c = 3.d0    !leaf in full growth
         end if
      end if

      !*Leaf-off in the fall 
      !*Leaf-off triggered by soil temperature
      !if soil temperature is falling below its maximum 
      !in the fall (when it's full-leaf (phenostatus=3.X) or leaf is senescening (phenostatus=4.X)), 
      !determine the phenofactor and corresponding phenostatus.
      if (fall .and.
     &   (phenostatus.ge.3.d0).and.
     &   (soiltemp_10d.lt.soilt_max+soilt_adj)) then
         !determine phenofactor by scaling soil temperature
         !with its minimum(min+adj) and maximum(max+adj).
         phenofactor_c = min(phenofactor_c,max(0.d0,
     &      (soiltemp_10d-soilt_min-soilt_adj)/(soilt_max-soilt_min)))
         if (phenofactor_c .eq. 0.d0) then
            phenostatus_c = 1.d0    !no leaf
            sgdd = 0.d0             !zero-out sgdd once complete leaf-off occurs.
         else  
            phenostatus_c = 4.d0    !leaf senescence
         end if 
      end if
    

      end subroutine pheno_update_coldherbaceous
      !*********************************************************************  
      subroutine pheno_update_drought(phenostatus, 
     i            mature, betad_10d, betad_min, betad_max, betad_res,
     o            phenofactor_d, phenostatus_d)
!@sum Update phenology for drought-deciduous PFTs
!@sum Called from pheno_update

      use ent_const

      !input variables
      real*8, intent(in) :: phenostatus !phenological status (refer pheno_update for details) 
      real*8, intent(in) :: mature      !X day to mature: mature=1+X/1000  (refer pheno_update for details) 
      real*8, intent(in) :: betad_10d   !10 day running mean of betad (calculated with stressH2O)  
      real*8, intent(in) :: betad_min   !betad minimum
      real*8, intent(in) :: betad_max   !betad maximum
      real*8, intent(in) :: betad_res   !betad resistance
      !output variables
      real*8, intent(out) :: phenofactor_d !phenological factor for drought deciduous (refer pheno_update for details) 
      real*8, intent(out) :: phenostatus_d !phenological status for drought deciduous (refer pheno_update for details)

      !*Leaf-on in the spring, triggered by water stress 
      !if betad is larger than its minimum 
      !in the spring (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X), 
      !and after long enough after the leaf-off), 
      !determine the phenofactor and corresponding phenostatus.
      if ((phenostatus .gt. mature) .and.
     &   (phenostatus.lt.3.d0).and. (betad_10d.gt.betad_min))then
         !determine phenofactor by scaling betad with its mim, max and resistance factor
         phenofactor_d = min(1.d0,
     &      ((betad_10d-betad_min)/(betad_max-betad_min))**betad_res)
         if (phenofactor_d .ge. 0.95d0) then
            phenostatus_d = 3.d0    !leaf in full growth
         else
            phenostatus_d = 2.d0    !growing leaf
         end if

      !*Leaf-off in the fall 
      !*Leaf-off triggered by water stress
      !if betad is falling smaller than its maximum 
      !in the fall (when there's no leaf (phenostatus=1.X) or leaf is growing (phenostatus=2.X)), 
      !determine the phenofactor and corresponding phenostatus.
      else if ((phenostatus.ge.3.d0).and. (betad_10d.lt.betad_max))then
         !determine phenofactor by scaling betad with its mim, max and resistance factor
         phenofactor_d = max(0.d0,
     &      ((betad_10d-betad_min)/(betad_max-betad_min))**betad_res)
         if (phenofactor_d .le. EPS) then
            phenostatus_d = 1.d0    !no leaf
         else
            phenostatus_d = 4.d0    !leaf senescence
         end if 
      end if    
      
      end subroutine pheno_update_drought

      !*********************************************************************   

      real*8 function turnover_acclim(cop)
!@sum  turnover_acclim  Seasonal acclimation factor (0-1) for turnover rate of plant tissues.
!@+    This function uses pft,phenotype, dbh, h, lai, and Sacclim if evergreen.
!@+    If evergreen, then assign the frost-hardiness factor (only affects cold climates).
!@+    If deciduous, then assign Cfol/Cfolmax (no turnover when dormant)
      use photcondmod, only : frost_hardiness
      type(cohort), pointer :: cop 
      !---Local-----
      integer :: pft
      
      pft = cop%pft
      if (pfpar(cop%pft)%phenotype.eq.EVERGREEN) then
         turnover_acclim = frost_hardiness(cop%Sacclim) 
!!      else
!!         turnover_acclim = cop%lai/(pfpar(pft)%sla*cop%n*
!!     &        Cfol_fn(pft,cop%dbh,cop%h))*1.d3
      else if (pfpar(cop%pft)%woody) then
         turnover_acclim = cop%lai/(pfpar(pft)%sla*cop%n)
     &        /(Cfol_fn(pft,cop%dbh,cop%h)*1.d-3)
      else
         turnover_acclim = cop%lai/(pfpar(pft)%sla*cop%n)
     &        /(Cfol_fn(pft,0.d0,5.d0)*1.d-3)  !Yeonjoo 5.d0 hi value
!!     &        /(Cfol_fn(pft,0.d0,1.5d0)*1.d-3)    !Ent GVSD perennial grass max height
                                                 !## Need to replace with ent_pft pft parameter.
                                                 !## Safest if hmax > hinit.
      endif
      end function turnover_acclim

      !*********************************************************************   

      subroutine veg_update(pp,config)
!@sum Update the vegetation state and carbon pools:
!@+   DAILY call.
!@+   Updates LAI, senescefrac, DBH, height 
!@+   carbon pools of foliage, sapwood, fineroot, hardwood, coarseroot
!@+   AND growth respiration from growth and tissue turnover.
!@+   Algorithms for this module are still subject to refinement - NK.
!@+   Current subroutine steps:
!@+   1. Calculate allocation ratios from the phenofactor.
!@+   2. Calculate litter from turnover, update Clabile and available CB_d.
!@+   3. Calculate Cactive_max given dbh, h, and alloc and phenofactor from 1.
!@+   4. Allocate CB_d to Cactive and update Clab.  However, if Clab<0 (I will
!@+   change this to daily maintenance respiration requirement), then senesce
!@+   to return Clab above the threshold value.  This can result in
!@+   Cactive>Cactive_max.  Why?
!@+   5. Update the active pools, Cfol, Csw, Cfroot
!@+   6. Given new Clab,  Csw and Cfol, re-calculate change in Cdead, Cactive,
!@+   and Crepro.  There are two options here:  SGrowthModel=1, update the
!@+   active pool given increase in Cdead;  or SGrowthModel=2, growth Cdead
!@+   only but not Cactive.
!@+   7. Update the patch-level Reproduction pool.
!@+   8. Update the cop active pools (again), capped by Cactive_max (new
!@+   Cactive_max or old?).
!@+   9. Update the structural pools, or allocate to litter if
!@+   do_structuralgrowth is FALSE.
!@+   10. HACK fix:  prevent cop cpools from going negative.
!@+   11. Calculate litter from growth.
!@+   12. Update the senescence fraction.
!@+   13. Update h and dbh (need conditional for if do_structuralgrowth=TRUE;
!@+   probably should do this same time as 9).
!@+   14. Update LAI.
!@+   15. Zero out CB_d.
!@+   16. After all cohorts are done, update Tpools from Clossacc for patch.
!@+   17. Update daily tissue growth respiration accumulated to distribute
!@+   over half-hourly fluxes. 

      use ent_const
      use cohorts, only : cohort_carbon
      use patches, only : patch_carbon
      use photcondmod, only :  frost_hardiness
      use respauto_physio, only : Resp_plant_day
      implicit none
      type(ent_config) :: config 
      type(patch),pointer :: pp 
      type(cohort), pointer :: cop
      integer :: pft
      logical :: woody
      logical :: is_annual 
      real*8 :: C_fol_old,C_froot_old,C_sw_old,C_hw_old,C_croot_old
      real*8 :: Cactive_old
      real*8 :: C_fol, C_froot, C_croot, C_sw, C_hw
      real*8 :: C_lab
      real*8 :: phenofactor !phenological elongation factor [0,1] (unitless)
      real*8 :: Cactive     !active carbon pool: foliage, sapwood, fine root (gC/pool/individual)
      real*8 :: Cactive_max !maximum active carbon pool allowed by the allometric constraint 
      real*8 :: Cdead       !dead carbon pool, including hardwood and coarse root (gC/pool/individual)
      real*8 :: qsw         !allometric factor for sapwood as a function of 
      real*8 :: dbh         !diameter at the breast height (cm)
      real*8 :: h           !plant height (m
      real*8 :: nplant      !plant population (#-individual/m2-ground)
      real*8 :: alloc,ialloc
      real*8 :: senescefrac !This is now net fraction of foliage that is litter.
      real*8 :: CB_d        !daily carbon balance (gC/individual)
      real*8 :: laipatch
      real*8 :: qf
      real*8 :: dCrepro
      real*8 :: dC_lab
      real*8 :: C_lab_old
      logical :: dormant
      real*8 :: dC_litter_hw
      real*8 :: dC_litter_croot
      real*8 :: Cfol_half
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      real*8 :: resp_auto_patch, resp_root_patch !kg-C/m/s
      integer :: cohortnum
      real*8 :: cpool(N_BPOOLS) 
      real*8 :: turn_leaf
      real*8 :: resp_growth
      real*8 :: resp_growth1
      real*8 :: resp_growth2
      integer, parameter :: irecruit=1
      real*8 :: dummy
      real*8 :: Rauto_day_gC ! Plant respiration estimate for one day.-NK
      real*8 :: facclim
      real*8 :: adj

      !Check carbon conservation
      real*8 :: tot_c_old, tot_c
      real*8 :: patch_tot_c_old, patch_tot_c
      real*8 :: d_tot_c(20), tot_closs(20)
      real*8 :: tot_closs_acc, tot_closs_acc_old
      real*8 :: cop_n_old, cop_n

      !Initialize
      laipatch = 0.d0
      Clossacc(:,:,:) = 0.d0 
      resp_auto_patch = 0.d0
      resp_root_patch = 0.d0
      cohortnum = 0
      resp_growth=0.d0
      resp_growth1=0.d0
      resp_growth2=0.d0
      cpool(:) = 0.d0
      dCrepro = 0.d0

      !C conservation check
      tot_closs_acc_old = 0.d0
      tot_closs_acc = 0.d0
      patch_tot_c_old = patch_carbon(pp)

      cop => pp%tallest

      do while(ASSOCIATED(cop))

         pft = cop%pft
         phenofactor = cop%phenofactor   

         tot_c_old = cop%n*cohort_carbon(cop)
         cop_n_old = cop%n
               
         cohortnum = cohortnum + 1

         is_annual = .false.

         is_annual = (pfpar(pft)%phenotype .eq. ANNUAL) !?Was commented out?
cddd
cddd        if (is_annual) then
cddd            if (phenofactor .gt. 0.d0 .AND. cop%C_fol .eq. 0.d0) then
cddd               cop%h = 0.05d0 !min. height = 0.05m
cddd	       select case (irecruit)
cddd	       case(1)
cddd               call recruit_annual(cop%pptr%Reproduction(pft),cop%h,
cddd     o              cop%C_fol, cop%C_froot,cop%n, cop%LAI)    
cddd               cop%pptr%Reproduction(pft) = 0.d0
cddd               case(2)
cddd               cop%C_fol = height2Cfol(pft,cop%h)
cddd!               cop%LAI=cop%n*pfpar(pft)%sla*
cddd               cop%LAI=cop%n*sla(pft,cop%llspan)*
cddd     &              (height2Cfol(pft,cop%h)/1000.0d0) 
cddd               cop%C_froot = q*cop%C_fol
cddd               end select 
cddd           end if
cddd         end if

 
         dbh = cop%dbh
         h = cop%h
         nplant = cop%n
         woody = pfpar(pft)%woody
         C_lab = cop%C_lab 
         C_fol = cop%C_fol
         C_froot = cop%C_froot
         C_sw = cop%C_sw
         C_hw = cop%C_hw
         C_croot = cop%C_croot 
         !Set array for passing values in - NK
         cpool(LABILE) = cop%C_lab
         cpool(FOL) = cop%C_fol
         cpool(FR) = cop%C_froot
         cpool(SW) = cop%C_sw
         cpool(HW) = cop%C_hw
         cpool(CR) = cop%C_croot

#ifdef DEBUG
!         write(207,*) cpool(LABILE),cpool(FOL),cpool(FR)
!     &        ,cpool(SW),cpool(HW),cpool(CR)
         write(207,*) cop%C_lab,cop%C_fol,cop%C_froot,cop%C_sw
     &        ,cop%C_hw,cop%C_croot
#endif

         !*************************************************
         !*Allometric relation - qsw, qf, ialloc
         !*************************************************
         !*calculate qsw 
         !qsw*h: ratio of sapwood to leaf biomass
         qsw  = sla(pft,cop%llspan)*iqsw_fn(pft)
         !for herbaceous, no allocation to the sap wood   
         if (.not.woody) qsw = 0.0d0      
         
         !*calculate qf
         !q: ratio of root to leaf biomass  
         qf = q 
         !for annual grasses, fine roots are assumed to be proportional to the foliage
         if (is_annual .and. pfpar(pft)%leaftype.eq.MONOCOT)
     &      qf=q*phenofactor 
         
         !*calculate total allocation relative to C_fol
         !*multiplying by C_fol gives total C for fol, froot, and sw.
         !phenofactor - for foliage, qf - for root, qsw*h - for sapwood
         alloc = phenofactor+qf+h*qsw
         !ialloc: inverse of alloc
         if (alloc .ne. 0.0d0) then
           ialloc = 1.d0/alloc
         else
           ialloc = 0.d0
         end if

         !**********************
         !*Litter from turnover
         !**********************
         !*save existing plant carbon pools 
         Cactive = C_froot + C_fol + C_sw !active = fine root + fol + sapw
         Cdead = C_hw + C_croot           !dead = heartwood + coarse root
         C_fol_old = C_fol
         C_froot_old = C_froot
         C_croot_old = C_croot
         C_sw_old = C_sw
         C_hw_old = C_hw
         Cactive_old =Cactive
         
         !*calculate the litter from turnover when mature growth
         if (cop%phenostatus.ge.3.d0 .and. cop%phenostatus.lt.4.d0) then
            call litter_turnover_cohort(SDAY,
     i           C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &           cop,Clossacc,
     &           turn_leaf,resp_growth1)
         endif

#ifdef DEBUG
         write(209,*) Clossacc(CARBON,:,:)
#endif

         !*update labile carbon from litter_turnover_cohort,
         !and resulting daily carbon balance
         C_lab_old =C_lab
         C_lab = cop%C_lab 
         CB_d = cop%CB_d - (C_lab_old-C_lab)

         !****************************************************
         !*Active growth: increment Cactive and decrease C_lab
         !****************************************************
         !*determine the potential max. of active pool
         Cactive_max=Cfol_fn(pft,dbh,h)*(alloc+(1.d0-phenofactor)) 

         !Get threshold Cfol_half 
         !call init_Clab(pft,dbh,h,Cfol_half)
         Cfol_half = 0.5d0 * Cfol_fn(pft,dbh,h)

         !*calculate growth of active carbon pools  
         facclim = frost_hardiness(cop%Sacclim)
         Rauto_day_gC = Resp_plant_day(pft,cpool
     &        ,cop%pptr%cellptr%TcanopyC + KELVIN
     &        ,cop%pptr%cellptr%Soiltemp(1) + KELVIN
     &        ,cop%pptr%cellptr%airtemp_10d + KELVIN
     &        ,cop%pptr%cellptr%soiltemp_10d + KELVIN
     &        ,facclim)*1.d3   !kgC/day to gC/day

!         if (.not.dormant) then
            call growth_cpools_active(pft,phenofactor,cop%phenostatus
     i           ,ialloc
     i           ,Rauto_day_gC,Cactive_max,Cfol_half,CB_d,C_fol
     o           ,Cactive,C_lab)
!         endif
#ifdef DEBUG
         write(203,*) Cactive_max, Cactive, C_lab, C_fol, dbh !##NK DEBUG
     &        ,phenofactor, ialloc, qf, h, qsw, h*qsw, Rauto_day_gC
#endif
      
         !*update the active carbon pools 
         cop%C_fol = phenofactor * Cactive *ialloc
!!!         if (cop%C_fol>Cfol_fn(pft,dbh,h)) then  !##NK dEBUG
!!!            print *,"C_fol>Cfol_max", Cactive_max, Cactive, pft, dbh
!!!     &           ,phenofactor, ialloc, qf, h, qsw, h*qsw
!            print *,(phenofactor + qf + h*qsw)
!!!         endif
         cop%C_froot = qf * Cactive * ialloc !qf = q * phenofactor
!         if (cop%phenostatus.ge.4.d0 .or. cop%phenostatus.eq.1.d0) then 
!            !senescing
!#ifdef DEBUG
!            write(208,*) 1, cop%phenostatus
!     &           ,C_sw_old, Cactive * h *qsw * ialloc, C_lab
!#endif            
!            cop%C_sw = min(C_sw_old, Cactive * h *qsw * ialloc)
!            adj = max(0.d0,Cactive * h *qsw * ialloc-cop%C_sw)
!            C_lab = C_lab + adj
!            Cactive = Cactive - adj 
!         else
            cop%C_sw = Cactive * h *qsw * ialloc
!         endif

         !*********************************************************
         !*Structural growth (& corresponding active, reproductive) 
         !*********************************************************
         !*determine whether the structural growth is activated or not.
         !for most of cases, it is activated,
         !but the structural growth of cold deciduous trees
         !only in the spring (when leaf is growing).
         dormant = .false.
         dormant =
     &      (pfpar(pft)%phenotype .eq. COLDDECID .and. 
     &      pfpar(pft)%woody .and. 
     &      cop%phenostatus .lt. 2.d0 .and. cop%phenostatus .ge. 3.d0 )
!     &      cop%phenostatus .lt. 2.d0 .or. cop%phenostatus .ge. 3.d0 )

         !*calculate the structural growth, 
          !and corresponding changes in active & reproductive pools
            if ((.not.dormant)) then
!            if ((.not.dormant).and.(cop%phenostatus.lt.4.d0)) then
               call growth_cpools_structural(pft,dbh,h,qsw,qf
     i              ,phenofactor, cop%phenostatus,dormant
     i              ,Rauto_day_gC,C_sw,Cactive_max,C_fol,CB_d
     i              ,Cactive_old
     o              ,Cactive,C_lab,Cdead,dCrepro) 
            cop%C_lab = C_lab !Update in litter_growth_cohort after checking do_structuralgrowth
            endif

         !*update Cactive_max in case there was structural growth.
!         Cactive_max=Cfol_fn(pft,dbh,h)*(alloc+(1.d0-phenofactor)) !NK
         
         !****************************************
         !*Update the active and structural pools
         !****************************************
         !*allocate the active carbon into foliage, fine root and sapwood, 
         !according to the allocation ratio
         !*foliage = phenofactor * active / (phenofactor + qf + h*qsw)

         !*First update alloc given possibly new h.
!!!         alloc = phenofactor+qf+h*qsw
!!!         !ialloc: inverse of alloc
!!!         if (alloc .ne. 0.0d0) then
!!!           ialloc = 1.d0/alloc
!!!         else
!!!           ialloc = 0.d0
!!!         end if
!         cop%C_fol = phenofactor *ialloc  * Cactive
!         !*froot = qf * active / (phenofactor + qf + h*qsw)
!         cop%C_froot =  qf * ialloc * Cactive
!         !*sw = (h * qsw) * active / (phenofactor + qf + h*qsw)  
!         cop%C_sw = min(Cactive,Cactive_max) * h *qsw * ialloc
!         cop%C_sw = h *qsw * ialloc * Cactive

         !********************
         !*Litter from growth 
         !********************
         !phenology + turnover + C_lab change + growth respiration
         !senescefrac returned is fraction of foliage that is litter.

!!! HACK !!!
! just to keep the things going reset negavive pools to zero
!!!         cop%C_fol = max( 0.d0, cop%C_fol)
!!!         cop%C_froot = max( 0.d0, cop%C_froot)
!!!         cop%C_croot = max( 0.d0, cop%C_croot)
!!!         cop%C_sw = max( 0.d0, cop%C_sw)
!!!         cop%C_hw = max( 0.d0, cop%C_hw)

         !*allocate structural growth and calculate the litter from growth
         !* and update cop%C_lab

         call litter_growth_cohort(config%do_structuralgrowth
     i        ,dCrepro, 
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     i        dC_litter_hw,dC_litter_croot,
     o        Cdead,cop,Clossacc,resp_growth2)
#ifdef DEBUG
         write(210,*) Clossacc(CARBON,:,:)
#endif

         !*update the fraction of senescence 
         if (C_fol_old.eq.0.d0) then
           cop%senescefrac = 0.d0
         else
           cop%senescefrac = l_fract *
     &          (max(0.d0,C_fol_old - cop%C_fol) + turn_leaf)/C_fol_old
         endif

         !****************************************
         !*Update the plant size, LAI & nitrogen
         !****************************************
         !*update dbh & height  
         if (woody) then !.and.(config%do_structuralgrowth)) then
            cop%dbh = Cdead2dbh(pft,Cdead)
            cop%h = dbh2height(pft,cop%dbh)
         else
            cop%dbh = 0.0d0
            cop%h = Cfol2height_herb(pft,cop%C_fol)
         end if
         
         !*update LAI
         !cop%LAI=cop%C_fol/1000.0d0*pfpar(pft)%sla*cop%n
         cop%LAI=cop%C_fol/1000.0d0*sla(pft,cop%llspan)*cop%n 
         !print *,'sla:',pfpar(pft)%sla,sla(pft,cop%llspan) !##NK DEBUG
         if (cop%LAI .lt. EPS) cop%LAI=EPS
         laipatch = laipatch + cop%lai  

         !*update Ntot
         cop%Ntot = cop%nm * cop%LAI 

         !* Summarize for patch level *!
         !Total respiration flux including growth increment.
         resp_auto_patch = resp_auto_patch  + cop%R_auto 
         resp_root_patch = resp_root_patch + cop%R_root

         !* Tissue growth respiration is distributed over day 
         !* subtracted at physical time step in canopy biophysics 
         !* module with R_auto.
         resp_growth = resp_growth1 + resp_growth2
         cop%C_growth = cop%C_growth + resp_growth*cop%n*1.d-3 !kg-C m-2 day-1
         cop%C_growth_flux = cop%C_growth/(24.d0*3600.d0) ! resp flux, kg-C m-2 s-1
          !*update the patch-level reproductive pool
          cop%pptr%Reproduction(cop%pft) = 
     &        cop%pptr%Reproduction(cop%pft)+ dCrepro*cop%n


         !*Zero-out the daily accumulated carbon 
         cop%CB_d = 0.d0   

         d_tot_c(cohortnum) = cop%n*cohort_carbon(cop) - tot_c_old
         tot_closs_acc = tot_closs_acc + Clossacc(CARBON,LEAF,1)
     &       +Clossacc(CARBON,FROOT,1)
     &       + Clossacc(CARBON,WOOD,1)
         tot_closs(cohortnum) = tot_closs_acc - tot_closs_acc_old
         tot_closs_acc_old = tot_closs_acc
         cop_n = cop%n
         cop%N_up = tot_closs_acc !### HACK ## TEMP USE OF UNUSED VARIABLE #-NK

cddd         write(901,*)  
cddd         write(901,*) "pft ", cop%pft
cddd         write(901,*) "deltaC ", tot_c - tot_c_old, tot_c_old
cddd         write(901,*) "deltaC*n ", (tot_c - tot_c_old)*cop%n
         !write(901,*) "Clossacc ", Clossacc


#ifdef PHENOLOGY_DIAG
         call phenology_diag(cohortnum,cop)  
#endif
            
         cop => cop%shorter 
      end do !looping through cohorts
  
      !*Update Tpool from all litter.
      call litter_patch(pp, Clossacc) 

      !*Update patch-level LAI
      pp%LAI = laipatch  

      !* Update patch fluxes with growth respiration. *!
      !* The daily respiration fluxes are accumulated in C_growth to
      !* be distributed over the course of a day to avoid pulses at night.
      pp%R_root = pp%R_root + resp_root_patch !Total flux including growth increment.


      patch_tot_c = patch_carbon(pp)

#ifdef DEBUG
      if( abs(patch_tot_c - patch_tot_c_old) > 1d-10 ) then
        write(903,*) "P ",patch_tot_c - patch_tot_c_old, patch_tot_c_old
     &      ,"C ",d_tot_c(1:cohortnum)*1000, "S ",tot_closs(1:cohortnum)
     &       ,"N ", cop_n_old, cop_n
      endif
#endif
      end subroutine veg_update

      !*********************************************************************
      subroutine growth_cpools_active(pft,phenofactor,phenostatus
     &     ,ialloc
     &     ,Rauto_day_gC,Cactive_max,Cfol_half,CB_d,C_fol,Cactive,C_lab)
!@sum Update Cactive and C_lab, reserving fraction for Reproduction.
!@+   This routine does NOT partition the individual active pools.
      use respauto_physio, only : Resp_plant_day
      !input variables
      integer, intent(in) :: pft
      real*8, intent(in) :: phenofactor
      real*8, intent(in) :: phenostatus
      real*8, intent(in) :: ialloc
      real*8,intent(in) :: Rauto_day_gC   !gC/plant/day estimated
      real*8, intent(in) :: Cactive_max
      real*8, intent(in) :: Cfol_half
      real*8, intent(in) :: CB_d
      real*8, intent(in) :: C_fol
      !output variables updated in the subroutine 
      real*8, intent(inout) :: Cactive
      real*8, intent(inout) :: C_lab
      !---Local variables----------------
      real*8 :: Cactive_pot !potential maximum of active carbon pools
      real*8 :: dC_lab      !change in labile [g-C/individual],  negative for reduction of C_lab for growth.
      real*8 :: dCactive    !change in active 
      real*8 :: dCavail     !available carbon for active pool
      !AGrowthModel options for active growth:
      !1 grass growth with no storage; 
      !2 grass growth with storage; 
      !3 grass growth with storage after the certain size; 
      !4 grass/tree growth with storage
      integer, parameter :: CPotModel = 2
      integer, parameter :: AGrowthModel= 3 

      !*******************************
      !*Calculate the change in C_lab 
      !*******************************
      !*1) Allocate the carbon to active pool, 
      !if both labile pool and daily carbon balance are positive             
!      if (C_lab .gt.0.d0 .and. CB_d .gt. 0.d0) then
#ifdef DEBUG      
      write(300,*) phenofactor
#endif
      if (C_lab .gt.(Rauto_day_gC) .and. CB_d .gt. 0.d0) then
        
         !*1-1) if there's no leaf, 
         !no change in labile and active pool;
         !i.e., accumulated labile carbon stays in the labile. 
         if (phenofactor .eq. 0.d0) then
!         if (phenostatus .lt. 2.d0) then
            dC_lab = 0.d0 !store the carbon in the labile
            dCactive = 0.d0

         !*1-2) if leaves are existing, 
         !certain amount of the labile carbon is relocated into the active carbon.
         else                   !growing
            !Cactive_max (max. allowed pool size according to the DBH)
            !Cactive_pot (current size + daily accumulated carbon)  
            !Cactive (current size)
            !Cactive_pot = Cactive + C_lab
            Cactive_pot = min(Cactive_max,Cactive + min(C_lab, CB_d)) !only new carbon is used for growth.
!            dCavail = min(Cactive_max, Cactive_pot) - Cactive
            dCavail = Cactive_pot - Cactive

            select case (AGrowthModel)
            case(1) !no storage - default
               dCactive = dCavail
            case(2) !storage and grass-only reproduction
               if (.not.pfpar(pft)%woody) then !herbaceous
                  dCactive = (1.d0-r_fract) * dCavail 
               else !woody
                  dCactive = dCavail
               end if
            case(3) !storage, with minimum C_fol required for herb repro.
               if (.not.pfpar(pft)%woody) then !herbaceous
                  if (C_fol .gt. Cfol_half) then
                     dCactive = (1.d0-r_fract) * dCavail 
                  else
                     dCactive = dCavail
                  end if
               else !woody
                  dCactive = dCavail
               end if
            case(4) !storage and repro, both grass and tree
               if (.not.pfpar(pft)%woody) then !herbaceous
                  dCactive = (1.d0-r_fract) * dCavail 
               else !woody
                  dCactive = (1.d0-r_fract) * dCavail
               end if
            end select

            !Update dC_lab according to dCactive
            if (dCactive .lt. 0.d0) then
               
               dC_lab = - dCactive * l_fract !translocated from active to labile
               !### NEED TO PUT IN LITTER HERE FOR SENESCED ACTIVE.

            else 
                dC_lab = - dCactive  !labile carbon allocated to active growth
            end if
         end if

      !*2) Translocate the active carbon to labile, 
      !if the labile carbon is negative.
!y      else if (C_lab .lt. 0.d0 ) then
!y         dCactive = C_lab / l_fract
!y         dC_lab = - C_lab 
      else if (C_lab .lt. Rauto_day_gC ) then !C_lab could be + or -.
         !* Senesce if C_lab drops low.
         dCactive = min(0.d0,(C_lab-Rauto_day_gC)) / (1.d0-l_fract) !negative change
         !dC_lab = - C_lab 
         dC_lab = -dCactive * l_fract  !translocation to C_lab
         !## NEED TO PUT IN LITTER FOR SENESCED ACTIVE.

      !*3) Otherwise, nothing happens, only respiration in biophysics module.
      else  
         dCactive = 0.d0
         dC_lab = 0.d0
      end if

      !******************************************
      !*Update the labile and active carbon pools
      !******************************************
      C_lab = C_lab + dC_lab
      Cactive = Cactive +dCactive

#ifdef DEBUG
      write(200,'(100(1pe16.8))') CB_d,C_lab,dC_lab,dCactive,Cactive, 
     &     Cactive_max,Cactive_pot, Rauto_day_gC,dCavail
#endif


#ifdef DEBUG
      if (Cactive>Cactive_max) then !## NK DEBGU
         print *,'Cactive>Cactive_max',Cactive_max,Cactive
     &        ,CB_d,Rauto_day_gC
      endif
#endif

      end subroutine growth_cpools_active

      !*********************************************************************
      subroutine growth_cpools_structural(
     i     pft,dbh,h,qsw,qf,phenofactor,
     i     phenostatus,dormant,
     i     Rauto_day, C_sw,Cactive_max,C_fol,CB_d,Cactive_old,
     o     Cactive, C_lab, Cdead, dCrepro)
!@sum Grow structural carbon: update Cactive, C_lab, Cdead, and Reproduction.
      
      implicit none
      integer, intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8, intent(in) :: h
      real*8, intent(in) :: qsw
      real*8, intent(in) :: qf
      real*8, intent(in) :: phenofactor
      real*8, intent(in) :: phenostatus
      logical, intent(in) :: dormant
      real*8, intent(in) :: Rauto_day
      real*8, intent(in) :: C_fol
      real*8, intent(in) :: CB_d
      real*8, intent(in) :: Cactive_max
      real*8, intent(in) :: C_sw
      real*8, intent(inout) :: Cactive,Cactive_old
      real*8, intent(inout) :: C_lab
      real*8, intent(inout) :: Cdead
      real*8, intent(out) :: dCrepro
      !---Local variables--------------
      real*8 :: Cavail        !available carbon for growth
      real*8 :: dCdead        
      real*8 :: dCactive
      real*8 :: gr_fract      !fraction of excess C going to structural growth
      real*8 :: rp_fract      !fraction of C to reproductive pool
      real*8 :: qs            
      real*8 :: qsprime
      real*8 :: dCfoldCdead
      real*8 :: dCfrootdCdead
      real*8 :: dHdCdead
      real*8 :: dCswdCdead
      !option for structural growth
      !1 - based on ED1 (Moorcroft et al.); 2 - based on ED2 (Medvigy et al.)
      integer :: SGrowthModel=1

      !**************************************************
      !*Calculate the growth fraction for different pools
      !**************************************************
      !*Herbaceous   
      if (.not.pfpar(pft)%woody) then 
         !all remained labile carbon can be allocated to different pools. 
         Cavail = C_lab
         !when labile carbon exists
         if (C_lab .gt. Rauto_day )then
            !if there's no leaf, 
            !carbon is allocated btw reproductive and active pools.
            if (phenofactor .eq. 0.d0) then
               qs = 0.d0 !no structural pools for grass 
               rp_fract = r_fract + c_fract  !THIS SUMS TO 1.0?? - NK
               gr_fract = 1.d0 - rp_fract
            !otherwise, labile carbon stays in the pool.
            else
               qs = 0.d0 !no structural pools for grass 
               rp_fract = 0.d0
               gr_fract = 0.d0
            end if
         !if labile pool is negative, 
         !carbon is translocated from the active pool.
         else 
            qs = 0.d0    !no structural pools for grass 
            rp_fract = 0.d0
            gr_fract = 1.d0 / l_fract
         end if

      !*Woody
      else 
         !C used for growth is limited by both size of avaiable labile storage 
         !& size of sapwood pool.
         !Cavail = min(C_lab,C_sw) !Yeonjoo's
         !2017-10-06: C for structural growth is limited by 3-year turnover time
         !            of Clabile pool. -NK
         
!         if ((dormant).or.(phenostatus.ge.4.d0)) then
!         if ((phenostatus.lt.2.d0).or.(phenostatus.ge.3.d0)) then
!            Cavail = 0.d0 !no growth when dormant or senescing
!         else
            !Cavail = min(max(0.d0,C_lab-Rauto_day),C_sw) !YKim
            !2 years = 2*365 days * secperday/secperday !NK daily alloc
            !3 years = 3*365 days * secperday/secperday !NK daily alloc
             Cavail = min(max(0.d0,C_lab-Rauto_day),C_sw)/(2.d0*365.d0)
!         endif
!         Cavail = min(C_lab, C_sw)

         !*1)Option 1: Growth, based on ED1
         !in this option, the active growth is accounted
         !corresponding to the structural growth
         if (SGrowthModel.eq.1) then 
            if (C_fol .gt. 0.d0 .and. 
     &           Cactive .ge. Cactive_max .and. C_lab .gt. 0.d0 ) then
               if (dbh .le. maxdbh(pft))then
                  dCfoldCdead = dDBHdCdead(pft,Cdead)
     &                          /dDBHdCfol(pft,C_fol)
                  dCfrootdCdead = qf *dCfoldCdead
                  dHdCdead = dHdDBH(pft, dbh)  * dDBHdCdead(pft,Cdead) 
                  dCswdCdead = qsw*
     &                 (h*dCfoldCdead + C_fol*dHdCdead)
                  qsprime
     &               = 1.d0 / (dCfoldCdead + dCfrootdCdead +
     &                 dCswdCdead)
                  qs=qsprime/(1.d0+qsprime)
                  rp_fract = r_fract
                  gr_fract = 1.d0 - rp_fract
               else
                  qs = 1.d0
                  rp_fract = r_fract
                  gr_fract = 1.d0 - rp_fract 
               end if
            else if (C_lab .le. Rauto_day)then
               qs = 0.d0
               rp_fract = 0.d0
               gr_fract = 1.d0 / l_fract
            else
               qs = 0.d0
               rp_fract = 0.d0
               gr_fract = 1.d0
            end if
        
         !*2) Option 2: Growth, based on ED2
         !in this option, available carbon goes to the structural pool 
         !without corresponding growth in the active pool
         else if (SGrowthModel.eq.2) then 
!            if (C_lab .gt. Rauto_day .and. CB_d .gt.0.d0 )then
            if (Cavail .gt. Rauto_day .and. CB_d .gt.0.d0 ) then
               qs = 1.d0
               rp_fract = r_fract
               gr_fract = 1.d0 - rp_fract 
            else
               qs = 0.d0
               rp_fract = r_fract
               gr_fract = 1.d0 - rp_fract
            end if
         end if
      end if
       
      dCdead = gr_fract * qs  * Cavail 
      dCactive = gr_fract *(1.d0 - qs) * Cavail !This can go >Cactive_max
      dCrepro =  rp_fract  * Cavail


      !************************
      !*update the carbon pools
      !************************
      C_lab = C_lab - (dCdead + dCactive + dCrepro)
      Cdead = Cdead + dCdead
      Cactive = Cactive + dCactive

#ifdef DEBUG
      write(201,'(100(1pe16.8))') Cavail,C_lab, Rauto_day
     &     ,dCactive, dCdead,dCrepro
#endif

      end subroutine growth_cpools_structural
!*************************************************************************
      subroutine recruit_annual(reproduction,height,
     &           C_fol,C_froot,nplant,lai)
!@sum Recruitment for annual grasses. 
!@+   Initializes new seedling carbon pools, density, LAI, given initial height
      real*8, intent(in) :: reproduction
      real*8, intent(in) :: height
      real*8, intent(out) :: C_fol
      real*8, intent(out) :: C_froot
      real*8, intent(out) :: nplant
      real*8, intent(out) :: lai
      real*8 :: recruit
      integer :: pft
      
      !this subroutine is sepcifically written for C3 annual grass temporarily
      pft=GRASSC3
      
      !amount of carbon, used for recruit
      recruit = reproduction * (1.d0 - mort_seedling) ! per patch-area

      !properties for new seedling
!      C_fol=height2Cfol(pft,height)  !Yeonjoo's
      C_fol = Cfol_fn(pft, 0.d0, height)
      C_froot = q *C_fol
      nplant = recruit / (C_fol + C_froot)
      LAI = nplant*pfpar(pft)%sla*(C_fol/1000.d0)
      
      end subroutine recruit_annual
!*************************************************************************

      subroutine assign_Closs(n, copfracroot, dC_fol, dC_froot, dC_sw
     i     ,dC_hw, dC_croot
     o     ,Closs, dC_total)
!@sum assign_Closs Generic routine to put litter fluxes into array.
      use cohorts, only : calc_CASArootfrac 
      real*8, intent(in) :: n   !cohort density (#/m^2)
      real*8, pointer :: copfracroot(:)
      real*8,  intent(in) ::    !g-C/plant
     &     dC_fol, dC_froot,dC_sw ,dC_hw, dC_croot
      real*8,intent(inout) :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !gC/m^2 Litter per cohort
      real*8,intent(out) :: dC_total !kg-C/m^2 !For carbonb balance check
      !---Local----
      integer :: i
      real*8 :: fracrootCASA(N_CASA_LAYERS)

      !assign root fractions for CASA layers -PK
      call calc_CASArootfrac(copfracroot,fracrootCASA)
      
      do i=1,N_CASA_LAYERS   
         if (i.eq.1) then       !only top CASA layer has leaf and wood litter
            Closs(CARBON,LEAF,i) = n * (1.d0-l_fract) * 
     &           max(0.d0,-dC_fol)
            Closs(CARBON,WOOD,i) = n * (
     &           max(0.d0,-dC_hw)  + 
     &           fracrootCASA(i) * 
     &           (max(0.d0,-dC_croot)))
         else    
            Closs(CARBON,LEAF,i) = 0.d0 
            Closs(CARBON,WOOD,i) = n * 
     &           fracrootCASA(i) * max(0.d0,-dC_croot)
         end if

         ! both layers have fine root litter. sapwood treated like fine roots?
         Closs(CARBON,FROOT,i) = n * (1.d0-l_fract)
     &        * fracrootCASA(i) * ( max(0.d0,-dC_froot)
     &        + max(0.d0,-dC_sw) )
      enddo
      
      dC_total = 0.d0
      do i=1,N_CASA_LAYERS 
         dC_total = dC_total - Closs(CARBON,LEAF,i) -
     &        Closs(CARBON,FROOT,i) - Closs(CARBON,WOOD,i)
      enddo
      dC_total = dC_total*1.d-3 ! convert it to kg
      
      end subroutine assign_Closs
!*************************************************************************

      subroutine accumulate_Clossacc(pft,Closs, Clossacc)
!@sum Accumulate litter into soil carbon pools.
      integer, intent(in) :: pft
      real*8,intent(in) :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort by depth.
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      !---Local-----
      integer :: i

      !loop through CASA layers-->cumul litter per pool per layer -PK
      do i=1,N_CASA_LAYERS

        !* Accumulate *!
        Clossacc(CARBON,LEAF,i) = Clossacc(CARBON,LEAF,i)
     &       + Closs(CARBON,LEAF,i)
        Clossacc(CARBON,FROOT,i) = Clossacc(CARBON,FROOT,i) 
     &       + Closs(CARBON,FROOT,i)
        Clossacc(CARBON,WOOD,i) = Clossacc(CARBON,WOOD,i) 
     &       + Closs(CARBON,WOOD,i)
      
        !* NDEAD POOLS *!
        Clossacc(CARBON,SURFMET,i) = Clossacc(CARBON,SURFMET,i) 
     &       + Closs(CARBON,LEAF,i) * solubfract(pft)
        Clossacc(CARBON,SOILMET,i) = Clossacc(CARBON,SOILMET,i) 
     &       + Closs(CARBON,FROOT,i) * solubfract(pft)
        Clossacc(CARBON,SURFSTR,i) = Clossacc(CARBON,SURFSTR,i)
     &       + Closs(CARBON,LEAF,i) * (1-solubfract(pft))
        Clossacc(CARBON,SOILSTR,i) = Clossacc(CARBON,SOILSTR,i) 
     &       + Closs(CARBON,FROOT,i) * (1-solubfract(pft))
        Clossacc(CARBON,CWD,i) = Clossacc(CARBON,CWD,i) 
     &       + Closs(CARBON,WOOD,i)
      end do   

      !* Return Clossacc *!
      end subroutine accumulate_Clossacc

      !*********************************************************************
      subroutine litter_turnover_cohort(dt,
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &        cop,Clossacc,turn_leaf,resp_growth)
!@sum litter_turnover_cohort. 
!@sum CALLED BY phenology veg_update.
!@sum DAILY TIME STEP.
!@sum Calculates at daily time step litterfall from cohort 
!@sum     to soil, tissue growth,growth respiration, and updates the following
!@sum     variables:
!@sum     cohort: C_lab
!@sum             C_growth (daily total tissue growth respiration),
!@sum             senescefrac
!@sum     patch:  Clossacc
      !* NOTES:
      !* Determine litter from cohort carbon pools and accumulate litter into
      !* Clossacc array.  
      !* Active pool loss to litter from turnover is replenished by same amount
      !* from C_lab, so no change to standing pools except for C_lab.
      !* Turnover tissue provides retranslocated carbon back to C_lab.
      !* No litter from sapwood.
      !* Dead pool loss to litter from turnover is replenished by same amount
      !* from C_lab, but without retranslocation. ## MAY WANT TO EXPERIMENT.
      !* Tissue growth respiration in C_growth is allocated by canopy
      !* biophysics module to fluxes over the course of the whole (next) day.
      !* After CASA, but called at daily time step. - NYK 7/27/06
      !* Update cohort pools - SUMMARY *!
      !C_fol replenished from C_lab: no change
      !C_froot replenished from C_lab: no change
      !C_sw =  No litter from sapwood
      !C_hw replenished from C_lab: no change
      !C_croot replenished from C_lab: no change

      use photcondmod, only : frost_hardiness
      real*8,intent(in) :: dt !seconds, time since last call
      real*8,intent(in) ::C_fol_old,C_froot_old,C_hw_old,C_croot_old,
     &     C_sw_old
      type(cohort),pointer :: cop
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      real*8,intent(inout) :: resp_growth
      !--Local-----------------
      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort.  !explicitly depth-structured -PK 7/07
      integer :: pft,i
      real*8 :: fracrootCASA(N_CASA_LAYERS)
      real*8 :: turnoverdtleaf !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtfroot !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtwood !Closs amount from intrinsic turnover of biomass pool.
!      real*8 :: turnoverdttotal!Total
      real*8 :: turn_leaf, turn_froot, turn_sw, turn_hw,turn_croot!g-C/plant
      real*8 :: dC_fol, dC_froot, dC_hw, dC_sw, dC_croot,dC_lab !g-C/individual
      real*8 :: adj !Adjustment to keep loss less than C_lab
      real*8 :: resp_growth_root          !g-C/individ/ms/s
      real*8 :: resp_turnover, resp_newgrowth !g-C/individ
      real*8 :: i2a !1d-3*cop%n -- Convert g-C/individual to kg-C/m^2
!      real*8 :: Csum
      real*8 :: dC_total, dClab_dbiomass
      real*8 :: tacclim !turnover rate factor
      real*8 :: dC_lab_corr_litter ! kg/m^2

      Closs(:,:,:) = 0.d0
      !Clossacc(:,:,:) = 0.d0 !Initialized outside of this routine

      !* Calculate fresh litter from a cohort *!
      pft = cop%pft
        
      !calculate turnover rates
      tacclim = turnover_acclim(cop)
      turnoverdtleaf = tacclim*cop%turnover_amp*annK(pft,LEAF)*SDAY !s^-1 * s/day = day^-1
!      turnoverdtleaf = tacclim*annK(pft,LEAF)*SDAY !s^-1 * s/day = day^-1
      turnoverdtfroot = tacclim*annK(pft,FROOT)*SDAY
      turnoverdtwood = 0.3d0/pfpar(pft)%lrage*
     &     (1.d0-exp(-annK(pft,WOOD)*SDAY)) !SW not HW, tuning factor.

      !turnover fluxes
      if (C_fol_old .gt. 0.d0) then
         turn_leaf = C_fol_old * turnoverdtleaf 
         turn_froot =  C_froot_old * turnoverdtfroot
         turn_sw = 0.d0
         turn_hw = C_hw_old * turnoverdtwood 
         turn_croot = C_croot_old * turnoverdtwood
      else
      !no turnover during the winter
      	turn_leaf = 0.d0
      	turn_froot= 0.d0
        turn_sw = 0.d0
      	turn_hw = 0.d0
      	turn_croot = 0.d0
      end if

      !respiration from turnover and regrowth.
      !* Coefficients from Amthor (2000) Table 3
      !* turn_hw is regrown as though sapwood to maintain a carbon balance. 
      resp_turnover = 0.16d0*turn_froot + 0.014d0*turn_leaf
     &     + 0.d0*turn_sw
      resp_newgrowth = 0.16d0*(max(0.d0,turn_hw)+max(0.d0,turn_croot)) 

      !turnover draws down C_lab from regrowth and retranslocation.
      !* (but not respiration, which is distributed over the day)
      dClab_dbiomass =
     &     max(0.d0,turn_leaf) + max(0.d0,turn_froot)
     &     + max(0.d0,turn_sw)
     &     + max(0.d0,turn_hw) + max(0.d0,turn_croot)
      dC_lab = - dClab_dbiomass  !Turnover growth
     &     + (1-l_fract)*(turn_leaf + turn_froot) !Retranslocation

      !* Calculate adjustment factor if loss amount is too large for C_lab.          !* Limit turnover litter if regrowth and respiration exceed C_lab.*!
      if (cop%C_lab+dC_lab-resp_turnover-resp_newgrowth.lt.0.d0) then
         if ((0.5d0*cop%C_lab -resp_newgrowth).lt.0.d0)
     &        then
            adj = 0.d0   !No turnover litter to preserve C_lab for growth.
                !C_lab will probably go negative here, but only a short while.
!        else                    !Reduce rate of turnover litter.
        else if ((turn_leaf + turn_froot).ne.0.d0)then
          adj = (0.5d0*cop%C_lab - dClab_dbiomass - resp_newgrowth)/
     &         ((1-l_fract)*(turn_leaf + turn_froot)
     &         + resp_turnover)
        else
          adj = 1.d0
        endif
      else
        adj = 1.d0
      endif

      !* Adjust turnover losses to accommodate low C_lab. *!
      turn_leaf = adj*turn_leaf
      turn_froot = adj*turn_froot
      turn_hw = adj*turn_hw
      turn_croot = adj*turn_croot

      if (adj.lt.1.d0) then
      !* Update C_lab again due to growth and retranslocation
      !*  (but not day respiration, which is distributed over the day).
      dC_lab = 
     &      (1-l_fract)*(turn_leaf + turn_froot) !retranslocation
     &     - max(0.d0,turn_leaf) - max(0.d0,turn_froot) !active regrowth
     &        - max(0.d0,turn_sw)
     &     - (max(0.d0,turn_hw)+max(0.d0,turn_croot)) !structural regrowth
      end if

      resp_growth_root = 0.16d0 * turn_froot + 0.16d0*turn_croot  
      resp_growth = 0.d0 !Turned off resp_growth_root + 0.14d0*turn_leaf+0.16d0*turn_hw  

      !* Calculate litter from turnover
      !* Change from senescence is calculated as max(0.d0, C_pool_old-C_pool).

      call assign_Closs(cop%n, cop%fracroot
     &     , -turn_leaf, -turn_froot, -turn_sw
     &     ,-turn_hw, -turn_croot, Closs, dC_total)


!#define RESTRICT_LITTER_FLUX
#ifdef RESTRICT_LITTER_FLUX
      call do_restrict_litter_flux(Closs, cop%C_lab*cop%n*1.d-3,
     &     dC_total, dC_lab_corr_litter)
      dC_lab = dC_lab + dC_lab_corr_litter/cop%n*1.d3
#endif
      cop%C_total = cop%C_total + dC_total


      !* Update C_lab *!
      cop%C_lab = cop%C_lab + dC_lab

      !* Return Clossacc *!
      call accumulate_Clossacc(pft, Closs, Clossacc)
!      Csum = 0.d0
!      do i=1,NPOOLS
!        Csum = Csum + Clossacc(CARBON,i,1)
!      enddo

      !* Return resp_growth for the day for C_growth *!

      end subroutine litter_turnover_cohort
      !*********************************************************************
      subroutine litter_growth_cohort(do_structuralgrowth
     i     ,dC_repro, !dt,dCrepro,
     i     C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &     dC_litter_hw,dC_litter_croot,
     &     Cdead,cop,Clossacc,resp_growth)
!@sum litter_growth cohort for prognostic growth.
!@+   CALLED BY phenology veg_update.
!@+   DAILY TIME STEP.
!@+   Calculates cohort litterfall and tissue growth respiration
!@+   Updates the following cohort variables:
!@+     cohort: C_lab
!@+             senescefrac
!@+   Outputs the following variables 
!@+     cohort:  resp_growth for C_growth (daily total tissue growth respir),
!@+     patch:  Clossacc (array for litter component fluxes)

      !* NOTES:
      !* Determine litter from cohort carbon pools due to pool size changes
      !*  and accumulate litter into Clossacc array.  
      !* No litter from sapwood.
      !* Since establishment from reproduction for woody plants is not
      !*  implemented, yet, dCrepro is put into foliage litter to conserve
      !*  carbon. TEMPORARY
      !* Tissue growth respiration in C_growth is allocated by canopy
      !* biophysics module to fluxes over the course of the whole (next) day.

      use photcondmod, only : frost_hardiness

      logical,intent(in) :: do_structuralgrowth
      real*8,intent(inout) ::dC_repro  !TEMPORARY inout.
      real*8,intent(in) ::C_fol_old,C_froot_old,C_hw_old,C_croot_old,
     &     C_sw_old
      real*8, intent(in) :: dC_litter_hw, dC_litter_croot
      real*8, intent(inout) :: Cdead !g-C/plant, structural carbon
      type(cohort),pointer :: cop
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !gC/m^2 Litter accumulator.
      real*8,intent(out)  :: resp_growth !g-C/individ (mass total per day)
      !--Local-----------------
      integer :: pft,i
      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !gC/m^2 Litter per cohort. 
      real*8 :: dC_fol, dC_froot, dC_hw, dC_sw, dC_croot,dC_lab !g-C/individual
!      real*8 :: adj !Adjustment to keep loss less than C_lab
      real*8 :: resp_growth_root !g-C/individ (mass total per day)
      real*8 :: dC_total, dClab_dbiomass
      real*8 :: Csum
      real*8 :: dC_lab_corr_litter ! kg/m^2

      Closs(:,:,:) = 0.d0
      !Clossacc(:,:,:) = 0.d0 !Initialized outside of this routine

      pft = cop%pft
        
      !* Structural plant pools change and structural pool litter.
      if (.not.do_structuralgrowth) then
      !* If structural growth is off (set in the rundeck),
      !* then any estimated changes in hw & croot becomes litter.        
         dC_hw = cop%C_hw - C_hw_old !Save difference
         dC_croot = cop%C_croot - C_croot_old
         cop%C_hw=C_hw_old
         cop%C_croot=C_croot_old
         Cdead = cop%C_hw+cop%C_croot !Keep Cdead static.
      else
      !* If structural growth is on (set in the rundeck),
      !* then change in dead pool is allocated into hw & croot.         
!         cop%C_hw = Cdead * hw_fract  !OLD DELETE!
!         cop%C_croot = Cdead * (1-hw_fract) !OLD DELETE!
         cop%C_hw = Cdead2Chw(pft, Cdead)
         cop%C_croot = Chw2Ccroot(pft,cop%C_hw)
         dC_hw = cop%C_hw - C_hw_old
         dC_croot = cop%C_croot - C_croot_old
      endif

      !*  Active plant tissue pools change 
      dC_fol = cop%C_fol-C_fol_old
      dC_froot = cop%C_froot - C_froot_old
      dC_sw = cop%C_sw - C_sw_old

      !* C_lab change for biomass growth or senescence (not turnover)
      dClab_dbiomass = +max(0.d0, dC_fol) + max(0.d0,dC_froot)!New tissue
     &     + max(0.d0,dC_sw)    
     &     + max(0.d0,dC_hw) + max(0.d0,dC_croot)
     &     + max(0.d0,dC_repro)
     &     - l_fract*( max(0.d0,-dC_fol) + max(0.d0,-dC_froot) !Retranslocation
     &     + max(0.d0,-dC_sw))

      !* Growth and retranslocation.
      !* NOTE: Respiration is distributed over the day by canopy module,
      !*       so does not decrease C_lab here.
      dC_lab = 
     &      -dClab_dbiomass 

      !* Respiration from biomass growth
      resp_growth_root = 0.16d0*(max(0.d0,dC_froot)+max(0.d0,dC_croot))
      resp_growth = resp_growth_root + 0.14d0 * 
     &      (max(0.d0,dC_fol)+max(0.d0,dC_sw)) + 0.16d0* max(0.d0,dC_hw) 
     
      !* Calculate litter from senescence and when do_structuralgrowth=FALSE 
      !* Change from senescence is calculated as max(0.d0, C_pool_old-C_pool).
      !* For non-annual plants, dump dCrepro into foliage litter until
      !*   ecological dynamics are implemented.  TEMPORARY
      if (pfpar(pft)%phenotype.eq.ANNUAL) then
         call assign_Closs(cop%n, cop%fracroot
     i        ,dC_fol, dC_froot, dC_sw,dC_hw, dC_croot
     o        ,Closs, dC_total)
      else
            call assign_Closs(cop%n, cop%fracroot
     i           ,dC_fol-dC_repro, dC_froot, dC_sw,dC_hw, dC_croot
     o           ,Closs, dC_total)
            dC_repro = 0.d0     !All put in litter.
      endif

!#define RESTRICT_LITTER_FLUX
#ifdef RESTRICT_LITTER_FLUX
      call do_restrict_litter_flux(Closs, cop%C_lab*cop%n*1.d-3,
     &     dC_total, dC_lab_corr_litter)
      dC_lab = dC_lab + dC_lab_corr_litter/cop%n*1.d3
#endif
      cop%C_total = cop%C_total + dC_total

      call accumulate_Clossacc(pft, Closs, Clossacc)

       !* Update C_lab *!
      cop%C_lab = cop%C_lab + dC_lab
      !at this point, C_lab<0 comes from the rounding errors...

      !* Return Clossacc and resp_growth *!
!      Csum = 0.d0
!      do i=1,NPOOLS
!        Csum = Csum + Clossacc(CARBON,i,1)
!      enddo

      !* Return resp_growth (plant) for cop%C_growth *!
      !* Tissue growth respiration is subtracted at physical time step
      !* distributed over day in canopy biophysics module with R_auto.


      !----------------Error check--------------------------------
#ifdef DEBUG        
!      if ( abs( (dC_fol+dC_froot+dC_hw+dC_sw+dC_croot
!     &     +dC_lab)*cop%n + Closs(CARBON,LEAF,1)+Closs(CARBON,FROOT,1)
!     &     + Closs(CARBON,WOOD,1) ) > 1d-10 ) then

         write(901,*) "Closs ", Closs(CARBON,LEAF,1)
     &       +Closs(CARBON,FROOT,1)
     &       + Closs(CARBON,WOOD,1)
     &       ,"dC ", dC_fol, dC_froot, dC_hw, dC_sw, dC_croot
     &       ,dC_lab
     &       ,"dC*n ", (dC_fol+dC_froot+dC_hw+dC_sw+dC_croot
     &       +dC_lab)*cop%n
     &       ,"dC_repro ", dC_repro, dC_repro*cop%n
!       endif

      if (cop%C_lab < -1.d-8) then
        write(902,*) "WARNING: Clab ", cop%C_lab, dC_lab, cop%pft
     &       ,"dC ", dC_fol, dC_froot, dC_hw, dC_sw, dC_croot
     &       ,dC_lab
     &       ,"dC*n ", (dC_fol+dC_froot+dC_hw+dC_sw+dC_croot
     &       +dC_lab)*cop%n
     &       ,"dCrepro ", dC_repro, dC_repro*cop%n
      endif
!      if (cop%C_lab < 0.d0) then
!         print*,dC_fol,cop%C_fol,dC_sw,cop%C_sw
!         print*,dC_lab,cop%C_lab, dC_froot, cop%C_froot
!         print*,dC_hw,cop%C_hw,dC_croot,cop%C_croot
!         stop
!      endif
#endif
      end subroutine litter_growth_cohort
      !*********************************************************************
      subroutine litter_cohort_fff(dt,
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &        cop,Clossacc)
!@sum litter_cohort for static woody structure 
!@sum CALLED BY ent_prescribed_updates.
!@sum DAILY TIME STEP.
!@sum Calculates litterfall from cohort to soil, tissue growth,
!@sum growth respiration, and updates the following
!@sum     variables:
!@sum     cohort: C_lab
!@sum             C_growth (daily total tissue growth respiration),
!@sum             senescefrac
!@sum     patch:  Clossacc
      !* NOTES:
      !* Determine litter from cohort carbon pools and accumulate litter into
      !* Clossacc array.  
      !* Active pool loss to litter from turnover is replenished by same amount
      !* from C_lab, so no change to standing pools except for C_lab.
      !* Turnover tissue provides retranslocated carbon back to C_lab.
      !* No litter from sapwood.
      !* Dead pool loss to litter from turnover is replenished by same amount
      !* from C_lab, but without retranslocation. ## MAY WANT TO EXPERIMENT.
      !* Tissue growth respiration in C_growth is allocated by canopy
      !* biophysics module to fluxes over the course of the whole (next) day.
      !* After CASA, but called at daily time step. - NYK 7/27/06
      !* Update cohort pools - SUMMARY *!
      !C_fol replenished from C_lab: no change
      !C_froot replenished from C_lab: no change
      !C_sw =  No litter from sapwood
      !C_hw replenished from C_lab: no change
      !C_croot replenished from C_lab: no change

      use cohorts, only : calc_CASArootfrac 
      use photcondmod, only : frost_hardiness
      real*8,intent(in) :: dt !seconds, time since last call
      real*8,intent(in) ::C_fol_old,C_froot_old,C_hw_old,C_croot_old,
     &     C_sw_old
      type(cohort),pointer :: cop
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      !--Local-----------------
      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort.  !explicitly depth-structured -PK 7/07
      integer :: pft,i
      real*8 :: fracrootCASA(N_CASA_LAYERS)
      real*8 :: turnoverdtleaf !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtfroot !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtwood !Closs amount from intrinsic turnover of biomass pool.
!      real*8 :: turnoverdttotal!Total
      real*8 :: turn_leaf, turn_froot, turn_hw,turn_croot, turn_live !g-C/individual
      real*8 :: dC_fol, dC_froot, dC_hw, dC_sw, dC_croot,dC_lab !g-C/individual
      real*8 :: adj !Adjustment to keep turnover less than C_lab
      real*8 :: resp_growth,resp_growth_root !g-C/individ/ms/s
      real*8 :: resp_turnover, resp_newgrowth !g-C/individ
      real*8 :: i2a !1d-3*cop%n -- Convert g-C/individual to kg-C/m^2
      real*8 :: Csum
      real*8 :: dC_total, dClab_dbiomass
      real*8 :: tacclim !turnover parameter, linear with Cfol/Cfolmax
      real*8 :: dC_lab_corr_litter ! kg/m^2

      Closs(:,:,:) = 0.d0
      !Clossacc(:,:,:) = 0.d0 !Initialized outside of this routine

      !* Calculate fresh litter from a cohort *!
      pft = cop%pft
        
      !assign root fractions for CASA layers -PK
      call calc_CASArootfrac(cop%fracroot,fracrootCASA)

      !* NLIVE POOLS *! 
      tacclim = turnover_acclim(cop)
      turnoverdtleaf = tacclim*cop%turnover_amp*annK(pft,LEAF)*SDAY !s^-1 * s/day = day^-1
!      turnoverdtleaf = tacclim*annK(pft,LEAF)*SDAY !s^-1 * s/day = day^-1
      turnoverdtfroot = tacclim*annK(pft,FROOT)*SDAY
      turnoverdtwood = 0.32d0/pfpar(pft)%lrage*
     &     (1.d0-exp(-annK(pft,WOOD)*SDAY)) !Sapwood not hardwood,coefficient is a tuning factor.

      !* Turnover draws down C_lab. *!
      !* Calculate adjustment factor if loss amount is too large for C_lab.
      turn_leaf = C_fol_old * turnoverdtleaf 
      turn_froot =  C_froot_old * turnoverdtfroot
      turn_live = turn_leaf + turn_froot

      !Wood losses:  
      turn_hw = C_hw_old * turnoverdtwood 
      turn_croot = C_croot_old * turnoverdtwood
      !## Need to add hw turn corresponding to woody allocation that goes to litter for static woody structure.

      !* Change in plant tissue pools. *!
      dC_fol = cop%C_fol-C_fol_old
      dC_froot = cop%C_froot - C_froot_old
      dC_hw = cop%C_hw - C_hw_old
      dC_sw = cop%C_sw - C_sw_old
      dC_croot = cop%C_croot - C_croot_old

      !* Distinguish respiration from turnover vs. from new growth.
      !* ### With constant prescribed structural tissue, dC_sw=0.d0,but
      !* ### there must still be regrowth of sapwood to replace that converted
      !* ### to dead heartwood.  For a hack, turn_hw is regrown as sapwood to
      !* ### maintain a carbon balance. 
      resp_turnover = 0.16d0*turn_froot + 0.014d0*turn_leaf !Coefficients from Amthor (2000) Table 3
      resp_newgrowth = 0.d0  !Turned off
!      resp_newgrowth = 0.16d0*max(0.d0,dC_froot) + 
!     &     0.14d0*(max(0.d0,dC_fol)+max(0.d0,dC_sw))
!     &     +0.16d0*(max(0.d0,turn_hw)+max(0.d0,turn_croot)) !##THIS IS RESPIRATION FOR REGROWTH OF SAPWOOD TO ACCOUNT FOR CONVERSION TO HEARTWOOD WITH CONSTANT PLANT STRUCTURE.

      !* C_lab required for biomass growth or senescence (not turnover)
      dClab_dbiomass = -max(0.d0, dC_fol) - max(0.d0,dC_froot) !Growth of new tissue
     &     - max(0.d0,dC_sw)    !For constant structural tissue, dC_sw=0, but still need to account for sapwood growth.
     &     - max(0.d0,turn_hw)-max(0.d0,turn_croot) !### Sapwood growth to replace that converted to heartwood.
     &     + l_fract*( max(0.d0,-dC_fol) + max(0.d0,-dC_froot) !Retranslocated carbon from senescence.
     &     - max(0.d0,-dC_sw))

      !* Growth and retranslocation.
      !* NOTE: Respiration is distributed over the day by canopy module,
      !*       so does not decrease C_lab here.
      dC_lab = 
     &     (l_fract)*(turn_leaf + turn_froot) !Retranslocated carbon from turnover
     &     + dClab_dbiomass       !Growth (new growth or senescence)
          !- resp_growth          !Distrib resp_growth in cop%C_growth over day.

      !* Limit turnover litter if losses and respiration exceed C_lab.*!
c      if (cop%C_lab+dC_lab-resp_turnover-resp_newgrowth.lt.0.d0) then
        if ((0.5d0*cop%C_lab + dClab_dbiomass-resp_newgrowth).lt.0.d0)
     &       then
         adj = 0.5         !Reduce turnover litter to preserve C_lab for growth.
c         print *,'DEBUG#: adj=0.5*adj',adj !C_lab will probably go negative here, but only a short while.

        else                    !Reduce rate of turnover litter.
c          adj = (0.5d0*cop%C_lab - dClab_dbiomass - resp_newgrowth)/
c     &         ((1-l_fract)*(turn_leaf + turn_froot)
c     &         + resp_turnover)
           adj = 1.d0
        endif
c      else
c        adj = 1.d0
c      endif


      !* Adjust turnover losses to accommodate low C_lab. *!
      if (adj < 1.d0) then
        turn_leaf = adj*turn_leaf
        turn_froot = adj*turn_froot
        turn_hw = adj*turn_hw
        turn_croot = adj*turn_croot
        
        !* Recalculate dClab_dbiomass *!
        dClab_dbiomass = max(0.d0, dC_fol) + max(0.d0,dC_froot) !Growth of new tissue
     &       + max(0.d0,dC_sw)  !For constant structural tissue, dC_sw=0, but still need to account for sapwood growth.
     &       + max(0.d0,turn_hw)+max(0.d0,turn_croot) !### This is sapwood growth to replace that converted to heartwood.
     &       - l_fract*( max(0.d0,-dC_fol) + max(0.d0,-dC_froot) !Retranslocated carbon from senescence.
     &       + max(0.d0,-dC_sw))
      endif

      !* Recalculate respiration. 
      !  Distinguish below- vs. above-ground autotrophic respiration.
        resp_growth_root = 0.16d0 * ( !Coefficient from Amthor (2000) Table 3
     &       turn_froot            !Turnover growth
     &       + max(0.d0,dC_froot)) !New biomass growth
     &       + 0.16d0*turn_croot   !# Hack for regrowth of sapwood converted to replace senesced coarse root.
        resp_growth = resp_growth_root + 
     &       0.14d0 *              !Coefficient from Amthor (2000) Table 3
     &       ( turn_leaf           !Turnover growth    
     &       +max(0.d0,dC_fol)+max(0.d0,dC_sw)) !New biomass growth
     &       + 0.16d0*turn_hw      !# Hack for regrowth of sapwood converted to replace senesced hw.
      !* Recalculate dC_lab in case adj < 1.0.
      dC_lab = 
     &     (l_fract)*(turn_leaf + turn_froot) !Retranslocated carbon from turnover
     &     + dClab_dbiomass       !Growth (new growth or senescence) 
          !- resp_growth          !Distrib resp_growth in cop%C_growth over day.

#ifdef DEBUG
      write(991,*)  tacclim,cop%turnover_amp,adj
     &       ,turn_froot,turn_croot,max(0.d0,dC_froot)
     &       ,max(0.d0,dC_croot), turn_leaf,turn_hw
     &       ,max(0.d0,dC_fol), max(0.d0,dC_sw) !New biomass growth
     &       ,resp_growth_root, resp_growth
#endif

      !* Calculate litter to soil from turnover and from senescence*!
      !* Change from senescence is calculated as max(0.d0, C_pool_old-C_pool).
      ! Senescefrac factor diagnostic also calculated.
      do i=1,N_CASA_LAYERS   
        if (i.eq.1) then        !only top CASA layer has leaf and wood litter -PK   
          Closs(CARBON,LEAF,i) = cop%n * (1.d0-l_fract) * (turn_leaf +
     &         max(0.d0,-dC_fol))
          Closs(CARBON,WOOD,i) = cop%n * (turn_hw + 
     &         max(0.d0,-dC_hw) +
     &         fracrootCASA(i)*
     &         (turn_croot+max(0.d0,-dC_croot)))
        else    
          Closs(CARBON,LEAF,i) = 0.d0 
          Closs(CARBON,WOOD,i) = cop%n * 
     &       (fracrootCASA(i)
     &         *(turn_croot+max(0.d0,-dC_croot)))
        end if
        ! both layers have fine root litter 
        Closs(CARBON,FROOT,i) = cop%n * (1.d0-l_fract)
     &       * fracrootCASA(i) 
     &       * (turn_froot + max(0.d0,-dC_froot))
      enddo

      !* Diagnostic
      dC_total = 0.d0
      do i=1,N_CASA_LAYERS 
        dC_total = dC_total - Closs(CARBON,LEAF,i) -
     &       Closs(CARBON,FROOT,i) - Closs(CARBON,WOOD,i)
      enddo
      dC_total = dC_total*1.d-3  ! convert it to kg

!#define RESTRICT_LITTER_FLUX
#ifdef RESTRICT_LITTER_FLUX
      call do_restrict_litter_flux(Closs, cop%C_lab*cop%n*1.d-3,
     &     dC_total, dC_lab_corr_litter)
      dC_lab = dC_lab + dC_lab_corr_litter/cop%n*1.d3
#endif
      cop%C_total = cop%C_total + dC_total

      call accumulate_Clossacc(pft,Closs, Clossacc)

!      write(992,*) C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
!     &     cop%C_lab,cop%C_fol,cop%C_froot,cop%C_hw,cop%C_sw,
!     &     cop%C_croot, cop%dbh,turn_leaf,turn_froot,turn_hw,turn_croot,
!     &     dC_fol,dC_froot,dC_hw,dC_sw,dC_croot,
!     &     Closs(CARBON,:,:), Clossacc(CARBON,:,:),adj,cop%turnover_amp,
!     &     tacclim,turnoverdtleaf,turnoverdtfroot, turnoverdtwood

      !################ ###################################################
      !#### DUE TO TIMING OF LAI UPDATE IN GISS GCM AT THE DAILY TIME STEP,
      !#### GROWTH RESPIRATION FROM CHANGE IN LAI NEEDS TO BE SAVED AS 
      !#### A RESTART VARIABLE IN ORDER TO SEND THAT FLUX TO THE ATMOSPHERE.
      !#### Igor has put in code to distribute C_growth over the day.
      !####################################################################

      cop%C_lab = cop%C_lab + dC_lab

      !* Tissue growth respiration is subtracted at physical time step
      !* distributed over day in canopy biophysics module with R_auto.
      cop%C_growth = cop%C_growth + resp_growth*cop%n*1.d-3
      cop%C_growth_flux = cop%C_growth/(24.d0*3600.d0) ! resp flux, C s-1

      !Cactive = Cactive - turn_leaf - turn_froot !No change in active

      !* Diagnostic
      if ( C_fol_old > 0.d0 ) then
        cop%senescefrac = l_fract *
     &     (max(0.d0,C_fol_old - cop%C_fol) + turn_leaf)/C_fol_old
        !cop%senescefrac = max(0.d0,-dC_fol/C_fol_old)
      else
        cop%senescefrac = 0.d0
      endif

      !* Return Clossacc *!
      Csum = 0.d0
      do i=1,NPOOLS
        Csum = Csum + Clossacc(CARBON,i,1)
      enddo

      end subroutine litter_cohort_fff


      subroutine litter_cohort(dt,
     i        C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &        cop,Clossacc)
!@sum litter_cohort for static woody structure 
!@sum CALLED BY ent_prescribed_updates.
!@sum DAILY TIME STEP.
!@sum Calculates litterfall from cohort to soil, tissue growth,
!@sum growth respiration, and updates the following
!@sum     variables:
!@sum     cohort: C_lab
!@sum             C_growth (daily total tissue growth respiration),
!@sum             senescefrac
!@sum     patch:  Clossacc
      !* NOTES:
      !* Determine litter from cohort carbon pools and accumulate litter into
      !* Clossacc array.  
      !* Active pool loss to litter from turnover is replenished by same amount
      !* from C_lab, so no change to standing pools except for C_lab.
      !* Turnover tissue provides retranslocated carbon back to C_lab.
      !* No litter from sapwood.
      !* Dead pool loss to litter from turnover is replenished by same amount
      !* from C_lab, but without retranslocation. ## MAY WANT TO EXPERIMENT.
      !* Tissue growth respiration in C_growth is allocated by canopy
      !* biophysics module to fluxes over the course of the whole (next) day.
      !* After CASA, but called at daily time step. - NYK 7/27/06
      !* Update cohort pools - SUMMARY *!
      !C_fol replenished from C_lab: no change
      !C_froot replenished from C_lab: no change
      !C_sw =  No litter from sapwood
      !C_hw replenished from C_lab: no change
      !C_croot replenished from C_lab: no change

      use cohorts, only : calc_CASArootfrac 
      !use biophysics, only: Resp_can_growth
      real*8,intent(in) :: dt !seconds, time since last call
      real*8,intent(in) ::C_fol_old,C_froot_old,C_hw_old,C_croot_old,
     &     C_sw_old
      type(cohort),pointer :: cop
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      !--Local-----------------
      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort.  !explicitly depth-structured -PK 7/07
      integer :: pft,i
      real*8 :: fracrootCASA(N_CASA_LAYERS)
      real*8 :: turnoverdtleaf !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtfroot !Closs amount from intrinsic turnover of biomass pool.
      real*8 :: turnoverdtwood !Closs amount from intrinsic turnover of biomass pool.
!      real*8 :: turnoverdttotal!Total
      real*8 :: turn_leaf, turn_froot, turn_hw,turn_croot, turn_live !g-C/individual
      real*8 :: dC_fol, dC_froot, dC_hw, dC_sw, dC_croot,dC_lab !g-C/individual
      real*8 :: adj !Adjustment to keep turnover less than C_lab
      real*8 :: resp_growth,resp_growth_root !g-C/individ/ms/s
      real*8 :: resp_turnover, resp_newgrowth !g-C/individ
      real*8 :: i2a !1d-3*cop%n -- Convert g-C/individual to kg-C/m^2
      real*8 :: Csum
      real*8 :: dC_total, dClab_dbiomass
      real*8 :: facclim !Frost hardiness parameter - affects turnover rates in winter.
      real*8 :: dC_lab_corr
      real*8 :: clab_init ! value of clab from init_Clab (g per plant)
      real*8 :: clab_max ! max clab computed from clab_init (g/m^2)
      real*8 :: dC_lab_corr_litter ! kg/m^2

      Closs(:,:,:) = 0.d0
      !Clossacc(:,:,:) = 0.d0 !Initialized outside of this routine

      !* Calculate fresh litter from a cohort *!
      pft = cop%pft
        
      !assign root fractions for CASA layers -PK
      call calc_CASArootfrac(cop%fracroot,fracrootCASA)

      !* NLIVE POOLS *! 
      facclim = frost_hardiness(cop%Sacclim)
      turnoverdtleaf = facclim*cop%turnover_amp*annK(pft,LEAF)*SDAY !s^-1 * s/day = day^-1
!      turnoverdtleaf = facclim*annK(pft,LEAF)*SDAY !s^-1 * s/day = day^-1
      turnoverdtfroot = facclim*annK(pft,FROOT)*SDAY
      turnoverdtwood = 0.32d0/pfpar(pft)%lrage*
     &     (1.d0-exp(-annK(pft,WOOD)*SDAY)) !Sapwood not hardwood.  0.08d0 is a tuning factor.

      !* Turnover draws down C_lab. *!
      !* Calculate adjustment factor if loss amount is too large for C_lab.
      turn_leaf = C_fol_old * turnoverdtleaf 
      turn_froot =  C_froot_old * turnoverdtfroot
      turn_live = turn_leaf + turn_froot

      !Wood losses:  
      turn_hw = C_hw_old * turnoverdtwood 
      turn_croot = C_croot_old * turnoverdtwood
      !## Need to add hw turn corresponding to woody allocation that goes to litter for static woody structure.

      !* Change in plant tissue pools. *!
      dC_fol = cop%C_fol-C_fol_old
      dC_froot = cop%C_froot - C_froot_old
      dC_hw = cop%C_hw - C_hw_old
      dC_sw = cop%C_sw - C_sw_old
      dC_croot = cop%C_croot - C_croot_old

      !* Distinguish respiration from turnover vs. from new growth.
      !* ### With constant prescribed structural tissue, dC_sw=0.d0,but
      !* ### there must still be regrowth of sapwood to replace that converted
      !* ### to dead heartwood.  For a hack, turn_hw is regrown as sapwood to
      !* ### maintain a carbon balance. 
      resp_turnover = 0.16d0*turn_froot + 0.014d0*turn_leaf !Coefficients from Amthor (2000) Table 3
      resp_newgrowth = 0.16d0*max(0.d0,dC_froot) + 
     &     0.14d0*(max(0.d0,dC_fol))+  0 * ( max(0.d0,dC_sw)) ! hack : sw was .14
     &     +0.16d0*(max(0.d0,turn_hw)+max(0.d0,turn_croot)) !##THIS IS RESPIRATION FOR REGROWTH OF SAPWOOD TO ACCOUNT FOR CONVERSION TO HEARTWOOD WITH CONSTANT PLANT STRUCTURE.

      !* C_lab required for biomass growth or senescence (not turnover)
      dClab_dbiomass = max(0.d0, dC_fol) + max(0.d0,dC_froot) !Growth of new tissue
     &     + max(0.d0,dC_sw)    !For constant structural tissue, dC_sw=0, but still need to account for sapwood growth.
     &     + max(0.d0,turn_hw)+max(0.d0,turn_croot) !### Sapwood growth to replace that converted to heartwood.
     &     - l_fract*( max(0.d0,-dC_fol) + max(0.d0,-dC_froot) !Retranslocated carbon from senescence.
     &     + max(0.d0,-dC_sw))

      !* Growth and retranslocation.
      !* NOTE: Respiration is distributed over the day by canopy module,
      !*       so does not decrease C_lab here.
      dC_lab = 
     &     - (1-l_fract)*(turn_leaf + turn_froot) !Retranslocated carbon from turnover
     &     - dClab_dbiomass       !Growth (new growth or senescence)
          !- resp_growth          !Distrib resp_growth in cop%C_growth over day.

      !* Limit turnover litter if losses and respiration exceed C_lab.*!
      if (cop%C_lab+dC_lab-resp_turnover-resp_newgrowth.lt.0.d0) then
        if ((0.5d0*cop%C_lab - dClab_dbiomass-resp_newgrowth).lt.0.d0)
     &       then
          adj = 0.d0            !No turnover litter to preserve C_lab for growth.
                                !C_lab will probably go negative here, but only a short while.
        else                    !Reduce rate of turnover litter.
          adj = (0.5d0*cop%C_lab - dClab_dbiomass - resp_newgrowth)/
     &         ((1-l_fract)*(turn_leaf + turn_froot)
     &         + resp_turnover)
        endif
      else
        adj = 1.d0
      endif

      !* Adjust turnover losses to accommodate low C_lab. *!
      if (adj < 1.d0) then
        turn_leaf = adj*turn_leaf
        turn_froot = adj*turn_froot
        turn_hw = adj*turn_hw
        turn_croot = adj*turn_croot
        
        !* Recalculate dClab_dbiomass *!
        dClab_dbiomass = max(0.d0, dC_fol) + max(0.d0,dC_froot) !Growth of new tissue
     &       + max(0.d0,dC_sw)  !For constant structural tissue, dC_sw=0, but still need to account for sapwood growth.
     &       + max(0.d0,turn_hw)+max(0.d0,turn_croot) !### This is sapwood growth to replace that converted to heartwood.
     &       - l_fract*( max(0.d0,-dC_fol) + max(0.d0,-dC_froot) !Retranslocated carbon from senescence.
     &       + max(0.d0,-dC_sw))
      endif

      !* Recalculate respiration. 
      !  Distinguish below- vs. above-ground autotrophic respiration.
        resp_growth_root = 0.16d0 * ( !Coefficient from Amthor (2000) Table 3
     &       turn_froot            !Turnover growth
     &       + max(0.d0,dC_froot)) !New biomass growth
     &       + 0.16d0*turn_croot   !# Hack for regrowth of sapwood converted to replace senesced coarse root.
        resp_growth = resp_growth_root + 
     &       0.14d0 *              !Coefficient from Amthor (2000) Table 3
     &       ( turn_leaf           !Turnover growth    
     &       +max(0.d0,dC_fol) )
     &       + 0.d0*(max(0.d0,dC_sw)) !New biomass growth (sw was 0.14d0)
     &       + 0.16d0*turn_hw      !# Hack for regrowth of sapwood converted to replace senesced hw.

!      write(991,*)  facclim,turn_froot,turn_croot,max(0.d0,dC_froot),
!     &     max(0.d0,dC_croot), turn_leaf,turn_hw,
!     &     max(0.d0,dC_fol), max(0.d0,dC_sw) !New biomass growth

      !* Recalculate dC_lab in case adj < 1.0.
      dC_lab = 
     &     - (1-l_fract)*(turn_leaf + turn_froot) !Retranslocated carbon from turnover
     &     - dClab_dbiomass       !Growth (new growth or senescence)
          !- resp_growth          !Distrib resp_growth in cop%C_growth over day.

      !* Calculate litter to soil from turnover and from senescence*!
      !* Change from senescence is calculated as max(0.d0, C_pool_old-C_pool).
      ! Senescefrac factor diagnostic also calculated.
      do i=1,N_CASA_LAYERS   
        if (i.eq.1) then        !only top CASA layer has leaf and wood litter -PK   
          Closs(CARBON,LEAF,i) = cop%n * (1.d0-l_fract) * (turn_leaf +
     &         max(0.d0,-dC_fol))
          Closs(CARBON,WOOD,i) = cop%n * (turn_hw + 
     &         max(0.d0,-dC_hw) +
     &         fracrootCASA(i)*
     &         (turn_croot+max(0.d0,-dC_croot)))
        else    
          Closs(CARBON,LEAF,i) = 0.d0 
          Closs(CARBON,WOOD,i) = cop%n * 
     &       (fracrootCASA(i)
     &         *(turn_croot+max(0.d0,-dC_croot)))
        end if
        ! both layers have fine root litter 
        Closs(CARBON,FROOT,i) = cop%n * (1.d0-l_fract)
     &       * fracrootCASA(i) 
     &       * (turn_froot + max(0.d0,-dC_froot))
      enddo

#ifdef ENT_DISABLE_CLAB_DUMPING_TO_SOIL
      continue ! do nothing
#else

#ifdef ENT_PROPER_CLAB_MAX
      call init_Clab(cop%pft, cop%dbh, cop%h, clab_init)
      clab_max = clab_init*2.d0*cop%n
#ifdef ENT_SMALLER_CLAB_MAX_FOR_GRASSC3PER
      if ( cop%pft == GRASSC3PER .or. cop%pft == GRASSC4 ) then
        clab_max = clab_init*1.d0*cop%n
      endif
#endif
      if ( cop%C_lab*cop%n > clab_max ) then
        dC_lab_corr = cop%C_lab*cop%n - clab_max
        dC_lab = dC_lab - dC_lab_corr/cop%n
        Closs(CARBON,LEAF,1) = Closs(CARBON,LEAF,1) + dC_lab_corr
      endif
#else
!!! hack to prevent C_lab from growing infinitely
      if ( cop%C_lab*cop%n > 3000.d0 ) then
        dC_lab_corr = cop%C_lab*cop%n - 3000.d0
        dC_lab = dC_lab - dC_lab_corr/cop%n
        Closs(CARBON,LEAF,1) = Closs(CARBON,LEAF,1) + dC_lab_corr
      endif
#endif

#endif

      !* Diagnostic
      dC_total = 0.d0
      do i=1,N_CASA_LAYERS 
        dC_total = dC_total - Closs(CARBON,LEAF,i) -
     &       Closs(CARBON,FROOT,i) - Closs(CARBON,WOOD,i)
      enddo
      dC_total = dC_total*1.d-3  ! convert it to kg

!#define RESTRICT_LITTER_FLUX
#ifdef RESTRICT_LITTER_FLUX
      call do_restrict_litter_flux(Closs, cop%C_lab*cop%n*1.d-3,
     &     dC_total, dC_lab_corr_litter)
      dC_lab = dC_lab + dC_lab_corr_litter/cop%n*1.d3
#endif
      cop%C_total = cop%C_total + dC_total

      call accumulate_Clossacc(pft,Closs, Clossacc)

!      write(992,*) C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
!     &     cop%C_lab,cop%C_fol,cop%C_froot,cop%C_hw,cop%C_sw,
!     &     cop%C_croot, cop%dbh,turn_leaf,turn_froot,turn_hw,turn_croot,
!     &     dC_fol,dC_froot,dC_hw,dC_sw,dC_croot,
!     &     Closs(CARBON,:,:), Clossacc(CARBON,:,:),adj,cop%turnover_amp,
!     &     facclim,turnoverdtleaf,turnoverdtfroot, turnoverdtwood

      !################ ###################################################
      !#### DUE TO TIMING OF LAI UPDATE IN GISS GCM AT THE DAILY TIME STEP,
      !#### GROWTH RESPIRATION FROM CHANGE IN LAI NEEDS TO BE SAVED AS 
      !#### A RESTART VARIABLE IN ORDER TO SEND THAT FLUX TO THE ATMOSPHERE.
      !#### Igor has put in code to distribute C_growth over the day.
      !####################################################################

      cop%C_lab = cop%C_lab + dC_lab

      !* Tissue growth respiration is subtracted at physical time step
      !* distributed over day in canopy biophysics module with R_auto.
      cop%C_growth = cop%C_growth + resp_growth*cop%n*1.d-3
      cop%C_growth_flux = cop%C_growth/(24.d0*3600.d0) ! resp flux, C s-1

      !Cactive = Cactive - turn_leaf - turn_froot !No change in active

      !* Diagnostic
      if ( C_fol_old > 0.d0 ) then
        cop%senescefrac = l_fract *
     &     (max(0.d0,C_fol_old - cop%C_fol) + turn_leaf)/C_fol_old
        !cop%senescefrac = max(0.d0,-dC_fol/C_fol_old)
      else
        cop%senescefrac = 0.d0
      endif

      !* Return Clossacc *!
      Csum = 0.d0
      do i=1,NPOOLS
        Csum = Csum + Clossacc(CARBON,i,1)
      enddo

      end subroutine litter_cohort


!*************************************************************************
      real*8 function frost_hardiness(Sacclim) Result(facclim)
!@sum frost_hardiness.  Calculate factor for adjusting photosynthetic capacity
!@sum  due to frost hardiness phenology.
      real*8,intent(in) :: Sacclim 
      !----Local-----
      real*8,parameter :: Tacclim=-5.93d0 ! threshold temperature for photosynthesis [deg C]
                        ! Site specific thres. temp.: state of photosyn.acclim
                        ! Hyytiala Scots Pine, -5.93 deg C Makela et al (2006)
      real*8,parameter :: a_const=0.0595 ! factor to convert from Sacclim [degC] to facclim [-]
                        ! Site specific; conversion (1/Sacclim_max)=1/16.8115
                        ! estimated by using the max S from Hyytiala 1998
!      real*8 :: facclim ! acclimation/frost hardiness factor [-]

      if (Sacclim > Tacclim) then ! photosynthesis occurs 
         facclim = a_const * (Sacclim-Tacclim) 
         if (facclim > 1.d0) facclim = 1.d0
!      elseif (Sacclim < -1E10)then !UNDEFINED
      elseif (Sacclim.eq.UNDEF)then !UNDEFINED
         facclim = 1.d0   ! no acclimation for this pft and/or simualtion
      else
         facclim = 0.01d0 ! arbitrary min value so that photosyn /= zero
      endif

      end function frost_hardiness



!**********************************************************************
      subroutine litter_patch(pp, Clossacc)
!@sum litter_dead.  Update soil Tpools following litterfall.
      type(patch),pointer :: pp 
      real*8,intent(in) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !gC/m^2 Litter accumulator.
      !----Local------
      integer :: i

      !* NDEAD POOLS *!
       do i=1,N_CASA_LAYERS
        pp%Tpool(CARBON,SURFMET,i) = pp%Tpool(CARBON,SURFMET,i) 
     &     + Clossacc(CARBON,SURFMET,i)
        pp%Tpool(CARBON,SOILMET,i) = pp%Tpool(CARBON,SOILMET,i) 
     &     + Clossacc(CARBON,SOILMET,i)
        pp%Tpool(CARBON,SURFSTR,i) = pp%Tpool(CARBON,SURFSTR,i)
     &     + Clossacc(CARBON,SURFSTR,i)
        pp%Tpool(CARBON,SOILSTR,i) = pp%Tpool(CARBON,SOILSTR,i) 
     &     + Clossacc(CARBON,SOILSTR,i)
        pp%Tpool(CARBON,CWD,i) = pp%Tpool(CARBON,CWD,i) 
     &     + Clossacc(CARBON,CWD,i)
       end do   !loop through CASA layers-->total C per pool per layer -PK

       ! update accumulators for external soil spinup
       i = 1
        pp%acc_litter_SURFMET = pp%acc_litter_SURFMET
     &     + Clossacc(CARBON,SURFMET,i)
        pp%acc_litter_SOILMET = pp%acc_litter_SOILMET
     &     + Clossacc(CARBON,SOILMET,i)
        pp%acc_litter_SURFSTR = pp%acc_litter_SURFSTR
     &     + Clossacc(CARBON,SURFSTR,i)
        pp%acc_litter_SOILSTR = pp%acc_litter_SOILSTR
     &     + Clossacc(CARBON,SOILSTR,i)
        pp%acc_litter_CWD = pp%acc_litter_CWD
     &     + Clossacc(CARBON,CWD,i)


       end subroutine litter_patch

!*********************************************************************
c      subroutine litter_old( pp)
c      !* Determine litter from live carbon pools and update Tpool.
c      !* After CASA, but called at daily time step. - NYK 7/27/06
c      
c      use cohorts, only : calc_CASArootfrac 
c
c      !real*8 :: dtsec           !dt in seconds
c      !type(timestruct) :: tt    !Greenwich Mean Time
c      type(patch),pointer :: pp
c      !--Local-----------------
c      type(cohort),pointer :: cop
c      real*8 :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter per cohort.  !explicitly depth-structured -PK 7/07
c      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
c      integer :: pft,i
c      real*8 :: fracrootCASA(N_CASA_LAYERS)
c      real*8 :: turnoverdtleaf !Closs amount from intrinsic turnover of biomass pool.
c      real*8 :: turnoverdtfroot !Closs amount from intrinsic turnover of biomass pool.
c      real*8 :: turnoverdtwood !Closs amount from intrinsic turnover of biomass pool.
c      real*8 :: turnoverdttotal, adj !Total, adjustment factor if larger than C_lab.
c
c      Closs(:,:,:) = 0.d0
c      Clossacc(:,:,:) = 0.d0
c
c      !* Calculate fresh litter from each cohort *!
c      cop => pp%tallest  !changed to => (!) -PK 7/11/06
c      do while(ASSOCIATED(cop)) 
c        pft = cop%pft
c        
c      !assign root fractions for CASA layers -PK
c      call calc_CASArootfrac(cop,fracrootCASA)
c!      print *, 'from litter(pheno*.f): fracrootCASA(:) =', fracrootCASA !***test*** -PK 11/27/06  
c
c       do i=1,N_CASA_LAYERS  !do this over all CASA layers -PK
c        !* NLIVE POOLS *! 
c        turnoverdtleaf = annK(pft,LEAF)*SDAY
c        turnoverdtfroot = annK(pft,FROOT)*SDAY
c        turnoverdtwood = 1.d0-exp(-annK(pft,WOOD)*SDAY) !Sapwood not hardwood
c
c        !* UPDATE C_LAB: Turnover should draw down C_lab. *!
c        ! Check that amount not too large 
c        turnoverdttotal = turnoverdtleaf+turnoverdtfroot+turnoverdtwood
c        !if (turnoverdttotal.lt.cop%C_lab) then
c          adj = 1.d0  !No adjustment, enough C_lab
c        !else
c        !  adj = (cop%C_lab - EPS)/turnoverdttotal !Turnover can reduce C_lab only to EPS.
c        !endif
c        !cop%C_lab = cop%C_lab - adj*turnoverdttotal
c        !* NEED TO PUT RETRANSLOCATION IN HERE, TOO *!
c
c        !* Calculate litter *!
c        ! Senescefrac factor can be calculated by either prescribed or prognostic phenology: ****** NYK!
c        if (i.eq.1) then  !only top CASA layer has leaf and wood litter -PK   
c         Closs(CARBON,LEAF,i) = 
c     &       cop%C_fol * cop%n *
c     &       (adj*turnoverdtleaf + cop%senescefrac) !* x tune factor
c         Closs(CARBON,WOOD,i) = 
c     &       (cop%C_hw + cop%C_croot) * cop%n
c     &       *(adj*turnoverdtwood) !* expr is kdt; x tune factor
c        else    
c         Closs(CARBON,LEAF,i) = 0.d0 
c         Closs(CARBON,WOOD,i) = 0.d0
c        end if
c         Closs(CARBON,FROOT,i) =  !both layers have root litter -PK 
c     &       fracrootCASA(i)*cop%C_froot * cop%n * 
c     &       (adj*turnoverdtfroot + cop%senescefrac) !* x tune factor
c
!        write(98,*) 'In litter: ',dtsec
!        write(98,*) cop%pft, cop%C_fol, cop%n, annK(pft,LEAF),pp%betad
!        write(98,*) cop%pft, cop%C_froot, cop%n, annK(pft,FROOT)
!        write(98,*) cop%pft, cop%C_hw, cop%n, annK(pft,WOOD)
!        write(98,*) 'solubfract(pft)', solubfract(pft)
!        write(98,*) 'Closs(CARBON,LEAF)',Closs(CARBON,LEAF)
!        write(98,*) 'Closs(CARBON,FROOT)',Closs(CARBON,FROOT)
c         Clossacc(CARBON,LEAF,i) = Clossacc(CARBON,LEAF,i)
c     &        + Closs(CARBON,LEAF,i)
c         Clossacc(CARBON,FROOT,i) = Clossacc(CARBON,FROOT,i) 
c     &        + Closs(CARBON,FROOT,i)
c         Clossacc(CARBON,WOOD,i) = Clossacc(CARBON,WOOD,i) 
c     &        + Closs(CARBON,WOOD,i)
c
c        !* NDEAD POOLS *!
c         Clossacc(CARBON,SURFMET,i) = Clossacc(CARBON,SURFMET,i) 
c     &        + Closs(CARBON,LEAF,i) * solubfract(pft)
c         Clossacc(CARBON,SOILMET,i) = Clossacc(CARBON,SOILMET,i) 
c     &        + Closs(CARBON,FROOT,i) * solubfract(pft)
c         Clossacc(CARBON,SURFSTR,i) = Clossacc(CARBON,SURFSTR,i)
c     &        + Closs(CARBON,LEAF,i) * (1-solubfract(pft))
c         Clossacc(CARBON,SOILSTR,i) = Clossacc(CARBON,SOILSTR,i) 
c     &        + Closs(CARBON,FROOT,i) * (1-solubfract(pft))
c         Clossacc(CARBON,CWD,i) = Clossacc(CARBON,CWD,i) 
c     &        + Closs(CARBON,WOOD,i)
c        end do  !loop through CASA layers-->cumul litter per pool per layer -PK
c     
c        cop => cop%shorter  !added -PK 7/12/06
c      end do  !loop through cohorts

      !* NDEAD POOLS *!
c       do i=1,N_CASA_LAYERS
c        pp%Tpool(CARBON,SURFMET,i) = pp%Tpool(CARBON,SURFMET,i) 
c     &     + Clossacc(CARBON,SURFMET,i)
c        pp%Tpool(CARBON,SOILMET,i) = pp%Tpool(CARBON,SOILMET,i) 
c     &     + Clossacc(CARBON,SOILMET,i)
c        pp%Tpool(CARBON,SURFSTR,i) = pp%Tpool(CARBON,SURFSTR,i)
c     &     + Clossacc(CARBON,SURFSTR,i)
c        pp%Tpool(CARBON,SOILSTR,i) = pp%Tpool(CARBON,SOILSTR,i) 
c     &     + Clossacc(CARBON,SOILSTR,i)
c        pp%Tpool(CARBON,CWD,i) = pp%Tpool(CARBON,CWD,i) 
c     &     + Clossacc(CARBON,CWD,i)
c       end do   !loop through CASA layers-->total C per pool per layer -PK
c!       print *, __FILE__,__LINE__,'pp%Tpool=',pp%Tpool(CARBON,:,:) !***test*** -PK 7/24/07  
c
c      end subroutine litter_old

!*********************************************************************
      real*8 function running_mean(dtsec,numd,var,var_mean) 
!@sum Function for expiring running mean for phenological climate statistics.
!@+   Daily average, taking into account time step of function call.
      real*8, intent(in) :: dtsec
      real*8, intent(in) :: numd !number of days for running mean
      real*8, intent(in) :: var
      real*8, intent(in) :: var_mean
      real*8 :: zweight

      zweight=exp(-1.d0/(numd*86400.d0/dtsec))
      running_mean=zweight*var_mean+(1.d0-zweight)*var  
      
      end function running_mean
!*************************************************************************

      subroutine phenology_diag(cohortnum, cop)
!@sum Debugging routine to print out phenology diagnostics to fort.990.
      implicit none
      type(cohort), pointer ::cop
      integer :: cohortnum
      real*8 :: fall_real
      if (cop%pptr%cellptr%fall==1) then
        fall_real=1.d0
      else
         fall_real=0.d0
      end if
         write(990,'(2(i5),27(1pe16.8))')
     &        cohortnum, !1
     &        cop%pft,  
     &        cop%phenofactor_c,
     &        cop%phenofactor_d,
     &        cop%phenostatus,
     &        cop%LAI,
     &        cop%C_fol,
     &        cop%C_lab,
     &        cop%C_sw,
     &        cop%C_hw,
     &        cop%C_froot,!11
     &        cop%C_croot,
     &        cop%NPP,
     &        cop%dbh,
     &        cop%h,
     &        cop%CB_d,
     &        cop%senescefrac,
     &        cop%llspan,
     &        cop%turnover_amp,
     &        cop%pptr%cellptr%airtemp_10d,
     &        cop%pptr%cellptr%soiltemp_10d, !21
     &        cop%betad_10d,
     &        cop%pptr%cellptr%par_10d,
     &        cop%pptr%cellptr%gdd,
     &        cop%pptr%cellptr%ncd,
     &        cop%pptr%cellptr%CosZen, 
     &        cop%pptr%cellptr%daylength(1),
     &        cop%pptr%cellptr%daylength(2),
     &        fall_real
!     &        cop%pptr%cellptr%fall

      end subroutine phenology_diag
!*************************************************************************

      subroutine do_restrict_litter_flux(Closs, C_lab_kg_m2,
     &     dC_total, dC_lab_corr_litter)
      real*8, intent(inout) :: Closs(PTRACE,NPOOLS,N_CASA_LAYERS) ! g/m^2
      real*8, intent(in) :: C_lab_kg_m2 ! kg/m^2
      real*8, intent(inout) :: dC_total ! kg/m^2
      real*8, intent(out) :: dC_lab_corr_litter ! kg/m^2
      !---
      real*8 :: scale

      if ( dC_total < 0.d0 .and. C_lab_kg_m2 < 0.d0 ) then
        scale = 1.d0 - min(1.d0, C_lab_kg_m2/dC_total)
        dC_lab_corr_litter = -dC_total * (1.d0 - scale)
        Closs(CARBON,:,:) = Closs(CARBON,:,:) * scale
        dC_total = dC_total * scale
      endif

      end subroutine do_restrict_litter_flux
      end module phenology
