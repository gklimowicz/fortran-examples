      module entcells
!@sum entcells contains routines to perform summaries at the ent grid cell
!@+   level.


      !Ent MODULES TO USE
      use ent_const
      use ent_types
      use ent_pfts
      use patches
      use cohorts

      implicit none
      private
      save

      public zero_entcell, summarize_entcell, entcell_print
      public entcell_print_diag
      public assign_entcell, assign_entcell_soilcarbon
      public init_simple_entcell, entcell_construct, entcell_destruct
      public entcell_extract_pfts, entcell_extract_heights
      public entcell_carbon

      contains
!**************************************************************************

      subroutine zero_entcell_patchsum(ecp)
!@sum Zeros entcell patch summary variables.
      implicit none
      type(entcelltype) :: ecp

      !--Cohort-level properties--!
      ecp%nm = 0.d0
      ecp%Ntot = 0.d0
      ecp%LMA = 0.d0
      ecp%LAI = 0.d0
      ecp%LAIpft(:) = 0.d0
      ecp%h = 0.d0
!      ecp%crown_dx = 0.d0
!      ecp%crown_dy = 0.d0
!      ecp%clump = 0.d0
      ecp%fracroot(:) = 0.d0
      ecp%C_fol = 0.d0
      ecp%N_fol = 0.d0
      ecp%C_w = 0.d0
      ecp%N_w = 0.d0
      ecp%C_lab = 0.d0
      ecp%N_lab = 0.d0
      ecp%C_froot = 0.d0
      ecp%N_froot = 0.d0
      ecp%C_root = 0.d0
      ecp%N_root = 0.d0
      ecp%Ci = 0.0127D0         !Internal foliage CO2 (mol/m3) non-zero
      ecp%GCANOPY = 0.d0
      ecp%GPP = 0.d0
      ecp%IPP = 0.d0
      ecp%NPP = 0.d0
      ecp%R_auto = 0.d0
      ecp%R_root = 0.d0  !PK 5/15/07
      ecp%N_up = 0.d0

      !--Patch-level properties, to be calculated from patch averages --!
      ecp%z0 = 0.d0
      ecp%albedo(:) = 0.d0
      ecp%betad = 0.d0
      ecp%betadl(:) = 0.d0
      ecp%TRANS_SW = 0.d0 
      ecp%CO2flux = 0.d0
      ecp%Soil_resp = 0.d0
      ecp%Tpool(:,:,:) = 0.d0

!      ecp%Soilmoist = 0.d0  !commented out -- it's a forcing, so shouldn't be altered -PK 7/24/07
!      ecp%N_deposit = 0.d0

!      nullify(ecp%crad%heights)
!      nullify(ecp%crad%lai)
!      ecp%crad%gortclump = 0.d0

      ecp%fuel = 0.d0
      ecp%ignition_rate = 0.d0
      ecp%lambda1(:) = 0.d0
      ecp%disturbance_rate(:) = 0.d0

      ecp%C_total = 0.d0
      ecp%C_growth = 0.d0

!      ecp%soil_type - meaning at this level?

      end subroutine zero_entcell_patchsum
!**************************************************************************

      subroutine zero_entcell(ecp)
!@sum Zeros import/export variables of entdata type
      implicit none
      type(entcelltype) :: ecp

      ecp%area = 0.d0

      call zero_entcell_patchsum(ecp)
      
      !VEGETATION - EXPORT STATE
      ecp%fv = 0.d0
      ecp%heat_capacity = 0.d0
      ecp%fwet_canopy = 0.d0
      !SOIL - IMPORT
      ecp%soil_Phi = 0.0         !Soil porosity (m3/m3)
      ecp%soildepth = 0.0        !Soil depth (m)
      ecp%theta_max = 0.0        !Saturated soil water volume (m/m)
      ecp%k_sat = 0.0            !Saturated hydraulic conductivity
      ecp%root_Phi = 0.0         !Infiltration factor promoted by roots (units?)

      !-----
      !METEOROLOGICAL - IMPORT STATE VARIABLES

      !VEGETATION - PRIVATE - Initial values not zero.
!      ecp%Ci = 0.0127D0         !Internal foliage CO2 (mol/m3) !!Cohort level or patch??
      ecp%Qf = 3.D-6            !Foliage surface vapor mixing ratio (kg/kg)
      
      !Cell-level summary values - CALCULATED BY GCM/EWB OR OFF-LINE FILE
      ecp%TairC = 0.0           !Air temperature (Celsius) 
      ecp%TcanopyC = 0.d0       !Canopy temperatue (Celsius)
!      ecp%Qv = 0.0               !Canopy air specif humidity (kg vapor/ kg air)
      ecp%P_mbar = 0.d0         !Atmospheric pressure (mb)
      ecp%Ca = 0.d0             !@Atmos CO2 conc at surface height (mol/m3).
      ecp%Soilmoist(:) = 0.d0   !Soil moisture (volumetric fraction), depth-structured  -PK 6/28/06
      ecp%Soiltemp(:) = 0.d0    !Soil temp (Celsius), depth-structured 
      ecp%fice = 0.d0           !Fraction of soil layer that is ice
      ecp%Ch = 0.d0             !Ground to surface heat transfer coefficient 
      ecp%U = 0.d0              !Surface layer wind speed (m s-1)

      !Radiation - IMPORT STATE VARIABLE
      !may later be broken down into hyperspectral increments in an array
      ecp%IPARdif = 0.0         !Incident diffuse PAR 400-700 nm (W m-2)
      ecp%IPARdir = 0.0         !Incident direct PAR 400-700 nm (W m-2)
      ecp%CosZen = 0.0          !cos of solar zenith angle

      !Phenology/Growth
      ecp%soiltemp_10d = 0.0d0 !0.7d0! for Hyytiala10-day running avergeage of soil temp (degC)
      ecp%airtemp_10d = 0.0d0 !-3.17!for Hyytiala  !10-day running average of air temp (degC) 
      ecp%paw_10d = 0.5d0!10-day running average of soil moisture (-)
      ecp%par_10d = 100.d0
      ecp%gdd = 0.0d0 !growing degree day
      ecp%ncd = 0.0d0 !number of chilling day
      ecp%sgdd = 0.d0
      ecp%daylength(:) = 0.0d0  !day length (min)
      ecp%fall = 1 !KIM- now it's integer...true. !KIM - starting in the winter
!      ecp%fall = .false. !KIM - starting in the middle of growing season

      ecp%C_total = 0.d0
      ecp%C_growth = 0.d0
      ecp%Soilmp(:) = 0.d0

      end subroutine zero_entcell
!**************************************************************************
      subroutine summarize_entcell(ecp)
!@sum Summarize patch properties within a grid cell.
      type(entcelltype) :: ecp
      !-----Local variables------------
      type(patch),pointer :: pp
      !real*8 :: albedo(N_BANDS)
      !real*8 :: z0              !Roughness length, average over patch
      !real*8 :: GCANOPY         !Canopy conductance of water vapor (mm/s)
      !real*8 :: CO2flux         !Net CO2 flux up (umol-CO2/m2-gnd/s)
      !real*8 :: Soilmoist(N_DEPTH) !Available soil moisture by depth (mm)
      !real*8 :: N_deposit       !N deposition (kgN/m2)
      !real*8 :: LAI
      !real*8 :: fracroot(N_DEPTH)
      integer :: ip             !#patches
      integer :: ia             !counter variable
      real*8 :: fa, laifa, laifasum
!      real*8 ::  vfraction(N_COVERTYPES) ! needed for a hack to compute canopy
                                 ! heat capacity exactly as in GISS GCM

      ecp%fv = 0.d0
      call zero_entcell_patchsum(ecp)

      ip = 0
      fa = 0.d0
      laifasum = 0.d0
      pp => ecp%oldest
      do while (ASSOCIATED(pp)) 
        ip = ip + 1

        call summarize_patch(pp) ! make sure patch is summarized, or assume

        fa = fa + pp%area       !For doing wtd avgs
        laifa = pp%area * pp%LAI !For doing wtd avgs
        laifasum = laifasum + laifa

        !- - - - Cohort - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ecp%nm = ecp%nm + pp%nm * laifa !wtd avg by LAI
        ecp%Ntot = ecp%Ntot + pp%Ntot * pp%area !wtd avg by area
        ecp%LMA = ecp%LMA + ecp%LMA * laifa !wtd avg by LAI

        ecp%LAI = ecp%LAI + pp%LAI*pp%area
        do ia=1,N_COVERTYPES
!          ecp%LAI = ecp%LAI + pp%LAIpft(ia) * pp%area !wtd avg by area
          ecp%LAIpft(ia) = ecp%LAIpft(ia) + pp%LAIpft(ia) * pp%area !wtd avg by area
        end do

        ecp%h = ecp%h + pp%h * laifa
!        ecp%crown_dx =  ecp%crown_dx + pp%crown_dx * laifa
!        ecp%crown_dy = ecp%crown_dy +  pp%crown_dy * laifa
!        ecp%clump = ecp%clump + pp%clump * laifa
        do ia=1,N_DEPTH
          ecp%fracroot(ia) = ecp%fracroot(ia) + pp%fracroot(ia)*pp%area
        end do

        ecp%C_fol = ecp%C_fol + pp%C_fol * pp%area
        ecp%N_fol = ecp%N_fol + pp%N_fol * pp%area
        ecp%C_w =ecp%C_w + pp%C_w * pp%area
        ecp%N_w = ecp%N_w + pp%N_w * pp%area
        ecp%C_lab = ecp%C_lab + pp%C_lab * pp%area
        ecp%N_lab = ecp%N_lab + pp%N_lab * pp%area
        ecp%C_froot = ecp%C_froot + pp%C_froot * pp%area
        ecp%N_froot = ecp%N_froot + pp%N_froot * pp%area
        ecp%C_root = ecp%C_root + pp%C_root * pp%area
        ecp%N_root = ecp%N_root + pp%N_root * pp%area

        !* IMPORT/EXPORT
        ecp%Ci = ecp%Ci + pp%Ci*pp%area
        ecp%GCANOPY = ecp%GCANOPY + pp%GCANOPY*pp%area

        !* DIAGNOSTICS
        ecp%GPP = ecp%GPP + pp%GPP*pp%area
        ecp%IPP = ecp%IPP + pp%IPP*pp%area
        ecp%NPP = ecp%NPP + pp%NPP*pp%area
        ecp%R_auto = ecp%R_auto + pp%R_auto*pp%area
        ecp%R_root = ecp%R_root + pp%R_root*pp%area  !PK 5/15/07
        ecp%N_up = ecp%N_up + pp%N_up*pp%area

        !- - - - Patch - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        !* IMPORT-PRESCRIBED, EXPORT-SIMULATED
        ecp%z0 = ecp%z0 + pp%z0*pp%area !Area-weighted average
        !* EXPORT
        do ia=1,N_BANDS  !wtd avg by area
          ecp%albedo(ia) = ecp%albedo(ia) + pp%albedo(ia)*pp%area
        end do

        ecp%betad = ecp%betad + pp%betad*pp%area
        do ia=1,N_DEPTH
          ecp%betadl(ia) = ecp%betadl(ia) + pp%betadl(ia)*pp%area
        end do
        ecp%TRANS_SW = ecp%TRANS_SW + pp%TRANS_SW*pp%area !Area-weighted average
        ecp%CO2flux = ecp%CO2flux + pp%CO2flux*pp%area
        
        !* DIAGNOSTICS
        ecp%Soil_resp = ecp%Soil_resp + pp%Soil_resp*pp%area   
        ecp%Tpool(:,:,:) = ecp%Tpool(:,:,:) + pp%Tpool(:,:,:)*pp%area                                                

        !* IMPORT Variables calculated by GCM/EWB - downscaled from grid cell
        ecp%Soilmoist(:) = ecp%Soilmoist(:) + pp%Soilmoist(:)*pp%area !##
        !ecp%N_deposit     !N deposition (kgN/m2)

        !* Variables for biophysics and biogeochemistry
        !type(canradtype) :: crad !Data structure for light profile

        !* Disturbance values
        ecp%fuel = ecp%fuel + pp%fuel*pp%area
        ecp%ignition_rate = ecp%ignition_rate+pp%ignition_rate*pp%area
        !ecp%lambda1(T_SUB)      !Site-averaged fire dist. rate during year
        !ecp%disturbance_rate(N_DIST_TYPES)

        !* Rates of change - patch total
        !ecp%dadt                !Rate of change of patch age = 1
        !ecp%dpdt                !Rate of change of patch area
        !ecp%dwdt                !Rate of change of available soil water

        !- - - - - Entcell - - - - - - - - - - - - - - - - - - - - - - - - - -

        if ( associated( pp%tallest ) ) ecp%fv = ecp%fv + pp%area
        !ecp%heat_capacity = !Currently imported

        !* daigs and hacks
        ecp%C_total = ecp%C_total + pp%C_total*pp%area
        ecp%C_growth = ecp%C_growth + pp%C_growth*pp%area
       
        pp => pp%younger
      end do
      
!!!!
    !!! hack to avoid divizion by zero
      !laifasum = max( laifasum, 1.d-10 )
      !------DO AVERAGES-------------------------------------------------
      !!!CHECK IF ECP%AREA IS ZERO!
      if (ASSOCIATED(ecp%oldest)) then
        !- - - - Cohort - - - - - - - - - - - - - - - - - - - - - - - - - 
        !laifasum weighting
         if (laifasum.gt.0.) then
            ecp%nm = ecp%nm/laifasum
            ecp%LMA = ecp%LMA/laifasum
            ecp%h = ecp%h/laifasum
         else
            ecp%nm = 0.d0
            ecp%LMA = 0.d0
            ecp%h = 0.d0
         endif

         ecp%Ntot = ecp%Ntot/fa
         ecp%LAI = ecp%LAI/fa
         
         do ia=1,N_COVERTYPES
            ecp%LAIpft(ia) = ecp%LAIpft(ia)/fa
         end do

         ecp%fracroot = ecp%fracroot/fa
         ecp%C_fol = ecp%C_fol/fa
         ecp%N_fol = ecp%N_fol/fa
         ecp%C_w = ecp%C_w/fa
         ecp%N_w = ecp%N_w/fa
         ecp%C_lab = ecp%C_lab/fa
         ecp%N_lab = ecp%N_lab/fa
         ecp%C_froot = ecp%C_froot/fa
         ecp%N_froot = ecp%N_froot/fa
         ecp%C_root = ecp%C_root/fa
         ecp%N_root = ecp%N_root/fa

        !* Flux variables for GCM/EWB - patch total wtd averages
        ecp%Ci = ecp%Ci/fa
        ecp%GCANOPY = ecp%GCANOPY/fa
        ecp%GPP = ecp%GPP/fa
C NADINE - IS THIS CORRECT?
        ecp%IPP = ecp%IPP/fa
        ecp%NPP = ecp%NPP/fa
        ecp%R_auto = ecp%R_auto/fa
        ecp%R_root = ecp%R_root/fa  !PK 5/15/07
        ecp%N_up = ecp%N_up/fa

        !- - - - - -  Patch-level summary values - PHYSICAL ------------------

        ecp%z0 = ecp%z0/fa !Area-weighted average
        do ia=1,N_BANDS         !Area-weighted average
          ecp%albedo(ia) = ecp%albedo(ia)/fa
        end do 
        ecp%betad = ecp%betad/fa
        ecp%TRANS_SW = ecp%TRANS_SW/fa !Area-weighted average
        ecp%CO2flux = ecp%CO2flux/fa
        ecp%Soil_resp = ecp%Soil_resp/fa  
        ecp%Tpool(:,:,:) = ecp%Tpool(:,:,:)/fa       
        
        !* Variables calculated by GCM/EWB - up/downscaled to/from grid cell
        ecp%Soilmoist(:) = ecp%Soilmoist(:)/fa
        do ia=1,N_DEPTH
          ecp%betadl(ia) = ecp%betadl(ia)/fa
        end do
        !ecp%N_deposit     !N deposition (kgN/m2)

        !* Variables for biophysics and biogeochemistry
        !type(canradtype) :: crad !Data structure for light profile

        !* Disturbance values
        ecp%fuel = ecp%fuel/fa
        ecp%ignition_rate = ecp%ignition_rate/fa

       !- - - - - Entcell - - - - - - - - - - - - - - - - - - - - - - - - - -
        !ecp%heat_capacity = !Currently imported
        !ecp%fv = no averaging needed

        ecp%area = fa !## NK used in ent_diag.
      end if

      ! use original formula for canopy heat cpacity
      !lai = ecp%LAI
      !ecp%heat_capacity=(.010d0+.002d0*lai+.001d0*lai**2)*shw*rhow
      !extract vfraction exactly like in GISS GCM
      !!vfraction(:) = 0.d0
      !!call entcell_extract_pfts(ecp, vfraction(2:) )
      !!ecp%heat_capacity=GISS_calc_shc(vfraction)

#ifndef MIXED_CANOPY
      call entcell_update_shc_mosaicveg(ecp)
#else
      call entcell_update_shc(ecp)
#endif

      end subroutine summarize_entcell
!**************************************************************************

      subroutine entcell_update_shc( ecp )
!@sum entcell_update_shc. 
!     Mixed canopies calculation of shc, or generic any veg structure.
!     Old GISS GCM version: shc(avg(lai*vfraction)).
!     Correct version: avg(shc(lai)*vfraction).
      use patches, only : shc_patch
!      use ent_prescr_veg, only : GISS_shc
      implicit none
      type(entcelltype) :: ecp
      !-----Local---------
      type(patch),pointer :: pp
      real*8 :: shc, pfrac

      shc = 0.d0
      pfrac = 0.d0
      pp => ecp%oldest
      do while (associated(pp))
         shc = shc + shc_patch(pp)*pp%area
         pfrac = pfrac + pp%area
         pp=>pp%younger
      enddo

      if (pfrac>EPS) then 
         shc = shc/pfrac        !Correct way.
      else
         !shc = GISS_shc(0.d0)  !Non-zero shc for zero lai from zero patch area.
         shc = 0.d0
      endif

      ecp%heat_capacity = shc
      end subroutine entcell_update_shc
!******************************************************************

      subroutine entcell_update_shc_mosaicveg( ecp )
!@sum Returns GISS GCM specific heat capacity for entcell.
!     This version preserves old, slightly incorrect lai averaging.
!     shc_entcell = shc(mean lai*vfraction)
!Old entcell_update_shc      
      use ent_const
      use ent_pfts, only: COVEROFFSET, alamax, alamin
      use allometryfn, only: do_geo,GISS_shc
      type(entcelltype) :: ecp
      !-----Local---------
      real*8 vfraction(N_COVERTYPES) ! needed for a hack to compute canopy
      type(patch),pointer :: pp      
      real*8 lai, fsum  !lai is mean annual entcell lai.
      integer pft, anum

      lai = 0.d0
      fsum = 0.d0
      vfraction(:) = 0.d0

      !Cover-weighted average
      if (.not.do_geo) then
         !Matthews mean annual LAI
         call entcell_extract_pfts( ecp, vfraction )
         do pft=1,N_PFT
            anum = pft+COVEROFFSET
            lai = lai + .5d0*(alamax(anum) 
     &           + alamin(anum))*vfraction(anum)
            fsum = fsum + vfraction(anum)
         enddo
      else
         !Geographic LAI - NYK
         !This calculation is just to preserve stability with R&A shc scheme.
         !Eventually, shc should be based on canopy biomass and water, not LAI.
         pp = ecp%oldest
         do while (ASSOCIATED(pp))
            lai = lai + laimean_annual_patch(pp)*pp%area
            fsum = fsum + pp%area
            pp => pp%younger
         enddo
      endif
      if ( fsum > EPS ) then 
         lai = lai/fsum
      else
         lai = 0.d0
      endif
         
      !shc = (.010d0+.002d0*lai+.001d0*lai**2)*shw*rhow
      
      ecp%heat_capacity=GISS_shc(lai)

      end subroutine entcell_update_shc_mosaicveg

!**************************************************************************

#ifdef SUMROOTSCELL
      subroutine sum_roots_patches2cell(ecp)
!@sum Calculate grid-averaged depth-, mass-, and cover-weighted average
!@sum of fine roots.
      type(entcelltype),pointer :: ecp
      !-----Local variables-------
      type(patch),pointer :: pp
      integer :: ia
      real*8 :: fracroot(N_DEPTH)
      real*8 :: frootC_total
      real*8 :: cfrac, tcfrac !cover fraction, total cover fraction

      !* Re-zero summary variable.
      do ia=1,N_DEPTH
        fracroot(ia) = 0.0
      end do
      frootC_total = 0.0
      cfrac = 0.0
      tcfrac = 0.0

      pp = ecp%oldest
      do while (ASSOCIATED(pp))
        cfrac = pp%area !This is area fraction.
        tcfrac = tcfrac + cfrac
        frootC_total = frootC_total + pp%C_froot
        do ia=1,N_DEPTH
          fracroot(ia) = fracroot(ia) + cfrac*pp%fracroot(ia)*pp%C_froot  
        end do
        pp = pp%younger
      end do

      ecp%fracroot = fracroot/(tcfrac*frootC_total)

      end subroutine sum_roots_patches2cell
#endif

!*************************************************************************
      subroutine init_simple_entcell( ecp,
     i     vegdata,popdens,laidata,hdata,dbhdata,craddata,
     i     cpooldata,nmdata,
     i     fracrootdata,soildata,albedodata,soil_texture,
     i     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpool_ini,
     i     reinitialize)
!@sum Initializes an entcell assuming one cohort per patch.
      use patches, only : summarize_patch
      use allometryfn, only : crown_radius_vert !Need to pass this in as array values
      type(entcelltype), target :: ecp
      real*8,intent(in) :: vegdata(N_COVERTYPES) !Veg cover fractions.
      real*8,intent(in) :: popdens(N_COVERTYPES) !Veg population density
      real*8,intent(in) :: laidata(N_COVERTYPES) !LAI
      real*8,intent(in) :: hdata(N_COVERTYPES) !Height
      real*8,intent(in) :: dbhdata(N_COVERTYPES) !Woody plant diameter (cm)
      real*8,intent(in) :: craddata(N_COVERTYPES)
      real*8,intent(in) :: cpooldata(N_COVERTYPES,N_BPOOLS)
      real*8,intent(in) :: nmdata(N_COVERTYPES) !Nitrogen parameter
      real*8,intent(in) :: fracrootdata(N_COVERTYPES,N_DEPTH) !Root profile.
      integer,intent(in) :: soildata(N_COVERTYPES)
      real*8,intent(in) :: albedodata(N_BANDS,N_COVERTYPES) !patch, NOTE:snow
      real*8,intent(in) :: soil_texture(N_SOIL_TEXTURES) !soil texture fractions.
      real*8 :: Ci_ini, CNC_ini, Tcan_ini, Qf_ini
      real*8,intent(in) ::
     &      Tpool_ini(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS)  !prescribed soil pools, g/m2 -PK
      logical,intent(in) :: reinitialize
      !-----Local---------
      integer :: ncov, pft
      type(patch),pointer :: pp, pp_tmp, pp_ncov
      real*8 :: sandfrac,clayfrac,smpsat,bch,watsat,watdry

      if ( reinitialize ) then
        ! destroy all existing patches since we are going to 
        ! re-initialize the cell
        pp => ecp%oldest      
        do while ( associated(pp) )
          pp_tmp => pp
          pp => pp%younger
          call patch_destruct(pp_tmp)
        enddo
        nullify(ecp%oldest)
        nullify(ecp%youngest)
      else
        ! just set all areas to 0 since we will reset them
        pp => ecp%oldest      
        do while ( associated(pp) )
          pp%area = 0.d0
          pp => pp%younger
        enddo

      endif

      do ncov=1,N_COVERTYPES           !One patch with one cohort per pft or bare
        pft = ncov - COVEROFFSET
        !### Get from GISS GCM ## vfraction of grid cell and area.

        if (vegdata(ncov)>0.0) then

          call get_patch_by_cover(ecp,ncov,pp_ncov)
          if ( associated(pp_ncov) ) then
            ! if patch is present - just resize it
            pp_ncov%area = vegdata(ncov)
            cycle
          endif
          
         !call insert_patch(ecp,GCMgridareas(j)*vegdata(pnum))
          call insert_patch(ecp,vegdata(ncov),soildata(ncov))
          pp => ecp%youngest

          !## Supply also geometry, clumping index

          ! insert cohort only if population density > 0 (i.e. skip bare soil)
          if ( popdens(ncov) > EPS ) then 
            if ( pft < 1 .or. pft > N_PFT ) then
              print *,"init_simple_entcell: wrong pft:", pft
              call patch_print(6,pp,"ERROR ")
              call stop_model("init_simple_entcell: wrong pft",255)
            endif
            call assign_patch(pp,Ci_ini, CNC_ini, pft, Tpool_ini)
            call insert_cohort(pp,pft,popdens(ncov),hdata(ncov),
     &           nmdata(ncov),laidata(ncov),
     &           craddata(ncov), 
!     &           min(0.45*hdata(ncov),craddata(ncov)*2.7d0), !old crad_dy vert NK
     &           crown_radius_vert(pft
     &             , hdata(ncov),craddata(ncov)), !dependency NK
     &           dbhdata(ncov),0.d0,0.d0, 
!     &           craddata(ncov),0.d0,dbhdata(ncov),0.d0,0.d0,
     &           0.d0, fracrootdata(ncov,:),
     &           cpooldata(ncov,FOL),0.d0,cpooldata(ncov,SW),0.d0,
     &           cpooldata(ncov,HW),0.d0,
     &           cpooldata(ncov,LABILE),0.d0,
     &           cpooldata(ncov,FR),0.d0,cpooldata(ncov,CR),0.d0,
     &           Ci_ini, CNC_ini,0.d0,0.d0,0.d0,0.d0, !added last 0 for R_root -PK 5/15/07
     &           0.d0, 0.d0, 
     &           0.d0, 1.d0,0.d0,1.d0, !KIM - starting in the winter for cold-dec.
!     &           1.d0, 1.d0,1.d0,3.d0, !KIM - starting in the middle of growing season
     &           1.d0, 0.d0,
     &           1.d0, -999.d0) !KIM-added for phenology/growth
          endif
          call summarize_patch(pp)

          !CALL CALC_ALBEDO HERE
          pp%albedo = albedodata(:,ncov) !##GISS HACK
        end if
      end do

      ! get rid of patches with 0 area
      pp => ecp%oldest      
      do while ( associated(pp) )
        pp_tmp => pp
        pp => pp%younger
        if( pp_tmp%area == 0.d0 ) call delete_patch(ecp, pp_tmp)
      enddo


      if ( reinitialize ) then
        !Initialize canopy met variables.
        ecp%TcanopyC = Tcan_ini
        ecp%Qf = Qf_ini

        ! soil textures for CASA
        ecp%soil_texture(:) = soil_texture(:)

      !Soil porosity and wilting? hygroscopic? point for water stress2 calculation. From soilbgc.f.
!??      sandfrac = pp%cellptr%soil_texture(1)
!??      clayfrac = pp%cellptr%soil_texture(3)
! is it supposed to be
        sandfrac = soil_texture(1)
        clayfrac = soil_texture(3)

        watsat =  0.489d0 - 0.00126d0*sandfrac !porosity, saturated soil fraction
        smpsat = -10.d0*(10.d0**(1.88d0-0.0131d0*sandfrac))
        bch = 2.91d0 + 0.159d0*clayfrac
        watdry = watsat * (-316230.d0/smpsat) ** (-1.d0/bch)
!      watopt = watsat * (-158490.d0/smpsat) ** (-1.d0/bch)
        ecp%soil_Phi = watsat
        ecp%soil_dry = watdry
      endif

#ifdef OFFLINE
      write(*,*) "soil_Phi, soil_dry",ecp%soil_Phi, ecp%soil_dry
#endif

      call summarize_entcell(ecp)

      
      !print *,"leaving init_simple_entcell:"
      !call entcell_print(6,ecp)

      end subroutine init_simple_entcell

!*************************************************************************
      subroutine assign_entcell( ecp,
     i     soil_texture,        !soil_type,
!     i     albedodata,
     i     Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpool_ini, 
     i     reinitialize)
!@sum assign_entcell. Assigns entcell level values as passed in parameters.
!+    NOTE:  soil_type is a patch-level variable by cover type in Matthews for
!+     the purpose of calculating albedo, whereas soil_texture is an 
!+     entcell-level variable. 

      use patches, only : summarize_patch
      implicit none
      type(entcelltype) :: ecp
      real*8,intent(in) :: soil_texture(N_SOIL_TEXTURES)!soil texture fractions.
!      real*8,intent(in) :: albedodata(N_BANDS,N_COVERTYPES) !#patch, NOTE:snow
!      integer,intent(in) :: soil_type(N_COVERTYPES) !#patch, from file or by vegtype.

      real*8 :: Ci_ini, CNC_ini, Tcan_ini, Qf_ini
      real*8,intent(in) ::
     &      Tpool_ini(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS)  
      logical,intent(in) :: reinitialize
      !-----Local---------
      integer :: ncov, pft
      type(patch),pointer :: pp, pp_tmp, pp_ncov
      real*8 :: sandfrac,clayfrac,smpsat,bch,watsat,watdry

      !Assign patch values.
!      pp => ecp%oldest      
!      do while ( associated(pp) )
!         pp%area = 0.d0
!         call assign_patch(pp,Ci_ini, CNC_ini, pft, Tpool_ini)
!         call summarize_patch(pp) 

         !CALL CALC_ALBEDO HERE 
!        !*GISS HACK VERSION for initialization.
         !* These values should get overwritten by prognostic 
         !* canopy radiative transfer. This provides non-zero default values.
!         if (associated(pp%tallest)) then
!            ncov = pp%tallest%pft + COVEROFFSET
!         else  !Bare soil type needed to be assigned from struct csv file.
!            ncov = SAND  !##HACK - NK - NEED TO ASSIGN FROM FILE.
!         endif
!        pp%albedo = ALBVND(ncov,1,:) !Winter albedo default. Should check hemi.
         !* Prognostic canopy radiation requires solar zenith angle
         !call get_patchalbedo(pp) !Currently dummy routine.  Use Two-stream.
!    
!      pp => pp%younger
!      enddo

      !Assign entcell variables.
      ecp%TcanopyC = Tcan_ini
      ecp%Qf = Qf_ini

      ! soil textures for CASA
      ecp%soil_texture(:) = soil_texture(:)

      !Soil porosity and wilting? hygroscopic? point for water stress2 calculation. From soilbgc.f.
!??      sandfrac = pp%cellptr%soil_texture(1)
!??      clayfrac = pp%cellptr%soil_texture(3)
! is it supposed to be
      sandfrac = soil_texture(1)
      clayfrac = soil_texture(3)

      watsat =  0.489d0 - 0.00126d0*sandfrac !porosity, saturated soil fraction
      smpsat = -10.d0*(10.d0**(1.88d0-0.0131d0*sandfrac))
      bch = 2.91d0 + 0.159d0*clayfrac
      watdry = watsat * (-316230.d0/smpsat) ** (-1.d0/bch)
!      watopt = watsat * (-158490.d0/smpsat) ** (-1.d0/bch)

      ecp%soil_Phi = watsat
      ecp%soil_dry = watdry


#ifdef OFFLINE
      write(*,*) "soil_Phi, soil_dry",ecp%soil_Phi, ecp%soil_dry
#endif

      call summarize_entcell(ecp)

      !print *,"leaving assign_entcell"
      !call entcell_print(6,ecp)

      end subroutine assign_entcell

  !*********************************************************************
      subroutine assign_entcell_soilcarbon(ecp,Tpools)
!@sum Distribute entcell grid soil carbon (e.g. from file) to subgrid pathes.
      use patches, only : assign_patch, summarize_patch
      implicit none
      type(entcelltype) :: ecp
      real*8,intent(in) ::
     &      Tpools(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS)  
      !--- Local ---
      type(patch),pointer :: pp

      pp => ecp%oldest      
      do while ( associated(pp) )
         pp%Tpool(NITROGEN,(NLIVE+1):NPOOLS,:) = 0.d0 !## HACK WHILE NO PROGNOSTIC N #
            if (ASSOCIATED(pp%tallest)) then
                call assign_patch_soilcarbon(pp,pp%tallest%pft,Tpools)
                call summarize_patch(pp) 
            endif
          pp => pp%younger
      enddo

      end subroutine assign_entcell_soilcarbon
  !*********************************************************************

      subroutine entcell_construct(ecp)
!@sum entcell_construct  Allocate memory for an entcell and nullify pointers.
      implicit none
      type(entcelltype), pointer :: ecp

      ! allocate memory
      allocate( ecp )
      allocate( ecp%LAIpft(N_COVERTYPES) )
      allocate( ecp%Soilmp(N_DEPTH) )
      allocate( ecp%fice(N_DEPTH) )
      allocate( ecp%fracroot(N_DEPTH) )
      allocate( ecp%betadl(N_DEPTH) )

      ! set pointers
      nullify( ecp%youngest )
      nullify( ecp%oldest   )

      ! for now set all values to zero or defaults
      call zero_entcell(ecp)

      ! maybe the following is not needed, but I guess we wanted
      ! to initialize the cells by default to bare soil (sand)
      call insert_patch(ecp,1.d0, soil_color_prescribed(SAND))

      end subroutine entcell_construct

 !*********************************************************************

      subroutine entcell_destruct(ecp)
!@sum entcell_destruct  Deallocate memory for an entcell.
      implicit none
      type(entcelltype), pointer :: ecp
      !---
      type(patch), pointer :: pp, pp_tmp

      ! destroy other patches
      pp => ecp%oldest      
      do while ( associated(pp) )
        pp_tmp => pp
        pp => pp%younger
        call patch_destruct(pp_tmp)
      enddo

      ! deallocate memory
      deallocate( ecp )
      nullify( ecp )

      end subroutine entcell_destruct

 !*********************************************************************
      
      subroutine entcell_print(iu, ecp)
!@sum entcell_print  Print contents of an entcell.
      use patches, only : patch_print
      integer, intent(in) :: iu
      type(entcelltype), intent(in) :: ecp
      !---
      type(patch), pointer :: pp
      character*1 :: prefix=" "
      character*8 prefix_p
      integer np

      write(iu, '(a,"entcell:")') prefix
      !print '(a," = ",f10.7)',"GCANOPY ",ecp%GCANOPY
      write(iu, 1) prefix,"soiltemp_10d",
     &     ecp%soiltemp_10d
      write(iu, 1) prefix,"airtemp_10d",
     &     ecp%airtemp_10d
      write(iu, 1) prefix,"paw_10d", ecp%paw_10d
      write(iu, 1) prefix,"par_10d", ecp%par_10d
      write(iu, 1) prefix,"gdd", ecp%gdd
      write(iu, 1) prefix,"ncd", ecp%ncd
      write(iu, 1) prefix,"ld", ecp%daylength(:)


      write(iu, '(a,a," = ",i7)') prefix,"fall", ecp%fall
      write(iu, 1) prefix,"nm	", ecp%nm	    
      write(iu, 1) prefix,"Ntot	", ecp%Ntot	    
      write(iu, 1) prefix,"LMA	", ecp%LMA	    
      write(iu, 1) prefix,"LAI      ", ecp%LAI      
      write(iu, 1) prefix,"h        ", ecp%h        
      write(iu, 1) prefix,"C_fol    ", ecp%C_fol    
      write(iu, 1) prefix,"N_fol    ", ecp%N_fol    
      write(iu, 1) prefix,"C_w      ", ecp%C_w      
      write(iu, 1) prefix,"N_w      ", ecp%N_w      
      write(iu, 1) prefix,"C_lab    ", ecp%C_lab    
      write(iu, 1) prefix,"N_lab    ", ecp%N_lab    
      write(iu, 1) prefix,"C_froot  ", ecp%C_froot  
      write(iu, 1) prefix,"N_froot  ", ecp%N_froot  
      write(iu, 1) prefix,"C_root   ", ecp%C_root   
      write(iu, 1) prefix,"N_root   ", ecp%N_root   
      write(iu, 1) prefix,"Ci       ", ecp%Ci       
      write(iu, 1) prefix,"GCANOPY  ", ecp%GCANOPY  
      write(iu, 1) prefix,"GPP      ", ecp%GPP      
      write(iu, 1) prefix,"IPP      ", ecp%IPP      
      write(iu, 1) prefix,"NPP      ", ecp%NPP      
      write(iu, 1) prefix,"R_auto   ", ecp%R_auto   
      write(iu, 1) prefix,"R_root   ", ecp%R_root   
      write(iu, 1) prefix,"N_up     ", ecp%N_up     
      write(iu, 1) prefix,"z0       ", ecp%z0       
      write(iu, 1) prefix,"albedo   ", ecp%albedo(:)
      write(iu, 1) prefix,"betad    ", ecp%betad    
      write(iu, 1) prefix,"TRANS_SW ", ecp%TRANS_SW 
      write(iu, 1) prefix,"CO2flux  ", ecp%CO2flux  
      write(iu, 1) prefix,"Soil_resp", ecp%Soil_resp
      write(iu, 1) prefix,"Tpool    ", ecp%Tpool
      write(iu, 1) prefix,"fuel	", ecp%fuel	    
      write(iu, 1) prefix,"ignit", ecp%ignition_rate
      write(iu, 1) prefix,"lambda1  ",ecp%lambda1(:)
      write(iu, 1) prefix,"dis",ecp%disturbance_rate
      write(iu, 1) prefix,"fv       ", ecp%fv       
      write(iu, 1) prefix,"heat_", ecp%heat_capacity
      write(iu, 1) prefix,"fwet_ca", ecp%fwet_canopy
      write(iu, 1) prefix,"soil_Phi ", ecp%soil_Phi 
      write(iu, 1) prefix,"soil_dry ", ecp%soil_dry 
      write(iu, 1) prefix,"soildepth", ecp%soildepth
      write(iu, 1) prefix,"theta_max", ecp%theta_max
      write(iu, 1) prefix,"k_sat    ", ecp%k_sat    
      write(iu, 1) prefix,"root_Phi ", ecp%root_Phi 
      write(iu, 1) prefix,"soil_t", ecp%soil_texture
      write(iu, 1) prefix,"TairC ! A", ecp%TairC
      write(iu, 1) prefix,"TcanopyC ", ecp%TcanopyC 
      write(iu, 1) prefix,"Qf       ", ecp%Qf       
      write(iu, 1) prefix,"P_mbar   ", ecp%P_mbar   
      write(iu, 1) prefix,"Ca       ", ecp%Ca       
      write(iu, 1) prefix,"Soilmoist", ecp%Soilmoist
      write(iu, 1) prefix,"Soiltem", ecp%Soiltemp(:)
      write(iu, 1) prefix,"Ch       ", ecp%Ch       
      write(iu, 1) prefix,"U        ", ecp%U        
      write(iu, 1) prefix,"IPARdir  ", ecp%IPARdir  
      write(iu, 1) prefix,"IPARdif  ", ecp%IPARdif  
      write(iu, 1) prefix,"CosZen   ", ecp%CosZen   
      write(iu, 1) prefix,"soilte", ecp%soiltemp_10d
      write(iu, 1) prefix,"airtemp", ecp%airtemp_10d
      write(iu, 1) prefix,"paw_10d  ", ecp%paw_10d  
      write(iu, 1) prefix,"par_10d  ", ecp%par_10d  
      write(iu, 1) prefix,"gdd      ", ecp%gdd	    
      write(iu, 1) prefix,"ncd	", ecp%ncd	    
      write(iu, 1) prefix,"daylength", ecp%daylength
      write(iu, 1) prefix,"C_total  ", ecp%C_total  
      write(iu, 1) prefix,"C_growth ", ecp%C_growth 

      if ( associated(ecp%LAIpft) )
     & write(iu,1) prefix,"LAIpft",ecp%LAIpft
      if ( associated(ecp%fracroot) )
     & write(iu,1) prefix,"fracroot",ecp%fracroot
      if ( associated(ecp%betadl) )
     & write(iu,1) prefix,"betadl",ecp%betadl
      if ( associated(ecp%Soilmp) )
     & write(iu,1) prefix,"Soilmp",ecp%Soilmp
      if ( associated(ecp%fice) )
     & write(iu,1) prefix,"fice",ecp%fice


      write(iu, '(a,"patches:")') prefix
      pp => ecp%oldest
      np = 0
      do while( associated(pp) )
        np = np + 1
        write( prefix_p, '(i2,"      ")' ) np
        call patch_print(iu, pp, prefix//prefix_p)
        pp => pp%younger
      enddo

 1    format(a,a," = ",99e23.16)  ! e12.5

      end subroutine entcell_print


 !*********************************************************************
      subroutine entcell_print_diag(iu, ecp)
!@sum Prints out other calculated diagnostics not printed by entcell_print
!@auth - NK
      use patches, only : patch_print_diag
      integer, intent(in) :: iu
      type(entcelltype), intent(in) :: ecp
      !---
      type(patch), pointer :: pp
      character*1 :: prefix=" "
      character*8 prefix_p
      integer np

      write(iu, '(a,"entcell:")') prefix

      write(iu, 1) prefix,"LMA	", ecp%LMA	    
      write(iu, 1) prefix,"LAI      ", ecp%LAI      
      write(iu, 1) prefix,"h        ", ecp%h        
      write(iu, 1) prefix,"C_fol    ", ecp%C_fol    
      write(iu, 1) prefix,"C_w      ", ecp%C_w      
      write(iu, 1) prefix,"C_lab    ", ecp%C_lab    
      write(iu, 1) prefix,"C_froot  ", ecp%C_froot  
      write(iu, 1) prefix,"C_root   ", ecp%C_root   
      write(iu, 1) prefix,"albedo   ", ecp%albedo(:)
      write(iu, 1) prefix,"Tpool    ", ecp%Tpool
      write(iu, 1) prefix,"fv       ", ecp%fv       
      write(iu, 1) prefix,"C_total  ", ecp%C_total  
      write(iu, 1) prefix,"C_growth ", ecp%C_growth 

      if ( associated(ecp%LAIpft) )
     & write(iu,1) prefix,"LAIpft",ecp%LAIpft

      write(iu, '(a,"patches:")') prefix
      pp => ecp%oldest
      np = 0
      do while( associated(pp) )
        np = np + 1
        write( prefix_p, '(i2,"      ")' ) np
        call patch_print_diag(iu, pp, prefix//prefix_p)
        pp => pp%younger
      enddo

 1    format(a,a," = ",99e23.16)  ! e12.5

      end subroutine entcell_print_diag

 !*********************************************************************

      subroutine entcell_extract_pfts(ecp, vfraction)
!@sum entcell_extract_pfts Extract cover fraction of subgrid patches.
      type(entcelltype) :: ecp
      real*8 :: vfraction(:)
      !---
      type(patch), pointer :: pp
      real*8 :: vfraction_patch(size(vfraction))

      vfraction(:) = 0.d0
      pp => ecp%oldest
      do while( associated(pp) )
        call patch_extract_pfts(pp, vfraction_patch)
        vfraction(:) = vfraction(:) + vfraction_patch(:)*pp%area
        pp => pp%younger
      enddo

      end subroutine entcell_extract_pfts

 !*********************************************************************

      subroutine entcell_extract_heights(ecp, h)
!@sum return maximum height per pft
      type(entcelltype) :: ecp
      real*8 :: h(:)
      !---
      type(patch), pointer :: pp
      type(cohort),pointer :: cop

      h(:) = 0.d0
      pp => ecp%oldest
      do while( associated(pp) )
        cop => pp%tallest
        if( associated(cop) ) then
          h(cop%pft + COVEROFFSET) = max(h(cop%pft + COVEROFFSET),cop%h)
        endif
        pp => pp%younger
      enddo

      end subroutine entcell_extract_heights

 !*********************************************************************

      subroutine get_patch_by_cover(ecp, ncov, pp_ncov)
!@sum Point to a patch of a given cover type in the given entcell.
      type(entcelltype) :: ecp
      integer :: ncov
      type(patch), pointer :: pp_ncov
      !---
      type(patch), pointer :: pp

      !write(0,*) "entered get_patch_by_cover", ncov
      nullify( pp_ncov )
      pp => ecp%oldest
      do while( associated(pp) )
        !write(0,*) "inside loop"
        if ( associated(pp%tallest) ) then
          !write(0,*) pp%tallest%pft+COVEROFFSET
          if ( pp%tallest%pft == ncov - COVEROFFSET ) then
            pp_ncov => pp
            return
          endif
        else
          if ( ncov == COVER_SAND .and. pp%soil_type == 1 ) then
            !write(0,*) 1
            pp_ncov => pp
            return
          else if ( ncov == COVER_DIRT .and. pp%soil_type == 2 ) then
            !write(0,*) 10
            pp_ncov => pp
            return
          endif
        endif
        pp => pp%younger
      enddo
      
      !write(0,*) "null"

      end subroutine get_patch_by_cover

!----------------------------------------------------------------------
      real*8 function entcell_carbon(ecp, 
     &     ecp_Cfol,ecp_Cstem,ecp_Croot, ecp_Csoil)
!@sum entcell_carbon (kg-C m-2). Return total carbon in entcell per land
!@+       area (exclude any water area in entcell grid).  
!@+       Optionally return component carbon pools.
      type(entcelltype) :: ecp
      real*8, optional :: ecp_Cfol
      real*8, optional :: ecp_Cstem
      real*8, optional :: ecp_Croot
      real*8, optional :: ecp_Csoil
      !--Local----
      type(patch), pointer :: pp
      real*8 :: ecp_kgCm2
      real*8 :: pp_Cfol, pp_Cstem, pp_Croot, pp_Csoil
      real*8 :: sumarea

      if (present(ecp_Cfol)) ecp_Cfol = 0.d0
      if (present(ecp_Cstem)) ecp_Cstem = 0.d0
      if (present(ecp_Croot)) ecp_Croot = 0.d0
      if (present(ecp_Csoil)) ecp_Csoil = 0.d0

      ecp_kgCm2 = 0.d0
      pp_Cfol = 0.d0
      pp_Cstem = 0.d0
      pp_Croot = 0.d0
      pp_Csoil = 0.d0
      sumarea = 0.d0

      pp => ecp%oldest

      do while (associated(pp))
         sumarea = sumarea + pp%area 
         ecp_kgCm2 = ecp_kgCm2 + pp%area*patch_carbon(pp,
     &        pp_Cfol, pp_Cstem, pp_Croot, pp_Csoil)

         if (present(ecp_Cfol)) ecp_Cfol = ecp_Cfol + pp%area*pp_Cfol
         if (present(ecp_Cstem)) ecp_Cstem = ecp_Cstem +pp%area*pp_Cstem
         if (present(ecp_Croot)) ecp_Croot = ecp_Croot +pp%area*pp_Croot
         if (present(ecp_Csoil)) ecp_Csoil = ecp_Csoil +pp%area*pp_Csoil

         pp => pp%younger
      end do

      if (sumarea.eq.0.d0) then
         entcell_carbon = 0.d0
      else
         entcell_carbon = ecp_kgCm2/sumarea
         if (present(ecp_Cfol)) ecp_Cfol = ecp_Cfol/sumarea
         if (present(ecp_Cstem)) ecp_Cstem = ecp_Cstem/sumarea
         if (present(ecp_Croot)) ecp_Croot = ecp_Croot/sumarea
         if (present(ecp_Csoil)) ecp_Csoil = ecp_Csoil/sumarea
      endif

      end function entcell_carbon
!----------------------------------------------------------------------

      end module entcells
