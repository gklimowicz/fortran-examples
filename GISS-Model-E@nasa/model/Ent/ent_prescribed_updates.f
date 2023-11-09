#define CHECK_CARBON_CONSERVATION

      module ent_prescribed_updates
!@sum Routines for updating prescribed (Matthews) vegetation. 
!@+   These routines work on the entcell level or lower.

!#define DEBUG TRUE  !NYK

      use ent_types

      implicit none
      private

      public entcell_vegupdate
      !public entcell_update_lai, entcell_update_albedo
      !public entcell_update_shc

      contains

!************************************************************************

      subroutine entcell_vegupdate(ecp, hemi, jday
     &     ,do_giss_phenology, do_giss_lai, do_giss_albedo!,mixed_veg
     &     ,laidata, hdata, albedodata, cropsdata, init )
!@sum Main routine for prescribed vegetation structure updates
!@+   at the entcell level (and down). 
!@+   DAILY TIME STEP ASSUMED. Coordinate-dependent
!@+   All var parameters except entcell are optional. 
!@+   Var arrays have pointer attribute to provide a way to tell the 
!@+   program that an argument is actually optional and not missing
!@+   (see how it is used in ent_prescribe_vegupdate)
      use ent_pfts, only: COVEROFFSET
      use patches, only : summarize_patch
      use entcells,only : summarize_entcell
      use ent_prescr_veg, only : prescr_veg_albedo
      implicit none
      type(entcelltype),pointer :: ecp
      integer,intent(in) :: jday
      integer,intent(in) :: hemi
      logical, intent(in) :: do_giss_phenology
      logical, intent(in) :: do_giss_lai
      logical, intent(in) :: do_giss_albedo
!      logical, intent(in) :: mixed_veg
      real*8,  pointer :: laidata(:)  !Array of length N_PFT
      real*8,  pointer :: hdata(:)  !Array of length N_PFT
      real*8,  pointer :: albedodata(:,:)
      real*8,  pointer :: cropsdata
      logical, intent(in) :: init
      !----Local------
      type(patch),pointer :: pp
      
      !* 1. Update crops to get right patch/cover distribution.
      !*     NOTE:  CARBON CONSERVATION NEEDS TO BE CALCULATED FOR CHANGING VEG/CROP COVER ##
      !* 2. Update height to get any height growth (with GISS veg, height is static)
      !* 3. Update LAI, and accumulate litter from new LAI and growth/senescence.
      !*       Cohort litter is accumulated to the patch level.
      !*    3a. If external LAI, then litter is calculated based on that external LAI change.
      !*    3b. If GISS prescribed LAI, then new LAI is calculated, and then litter.
      !* 4. Update albedo based on new vegetation structure.

      !* VEGETATION STRUCTURE AND LITTER *!
      ! veg structure update with external data if provided
      if (.not.do_giss_lai) then

        if ( associated(cropsdata) )
     &       call entcell_update_crops(ecp, cropsdata)
        !print *,"update cropsdata: ", associated(cropsdata) !##debug
        
        if ( associated(hdata) )
     &       call entcell_update_height(ecp, hdata, init)!, mixed_veg)
        !print *, "update hdata: ", associated(hdata) !##debug

        if ( associated(laidata) )
     &       call entcell_update_lai_poolslitter(ecp,laidata,
     &       init)!,mixed_veg)
        !print *, "update laidata: ",associated(laidata) !##debug

      endif
      ! or veg structure from prescribed GISS LAI phenology 
      if ( do_giss_phenology ) then !do_giss_phenology is redundant with do_giss_lai.
        if ( hemi<-2 .or. jday <-2 )
     &       call stop_model("entcell_vegupdate: needs hemi,jday",255)
        pp => ecp%oldest
        do while (ASSOCIATED(pp))
        !* LAI, SENESCEFRAC *!
          if (do_giss_lai) 
     &         call prescr_phenology(jday,hemi, pp, do_giss_lai)
         !call summarize_patch(pp) !* Redundant because summarize_entcell is called.
          pp => pp%younger
        end do
      endif


      !* ALBEDO *!
      if ( associated(albedodata) ) then
        call entcell_update_albedo(ecp, albedodata)
        !print *, "update albedodata from array"
      else
        pp => ecp%oldest
        do while (ASSOCIATED(pp))
          ! update if have vegetation or not prognostic albedo
          if ( ASSOCIATED(pp%tallest).and.do_giss_albedo )
     &         call prescr_veg_albedo(hemi, pp%tallest%pft+COVEROFFSET, 
     &         jday, pp%albedo)
          pp => pp%younger
        end do
        !print *, "update albedodata hemi"
      endif

      call summarize_entcell(ecp)
!      print *,'ecp%LAI', ecp%LAI
!      print *,'ecp%oldest%LAI', ecp%oldest%LAI

!Moved this into summarize_entcell
!#ifndef MIXED_CANOPY
!      call entcell_update_shc_mosaicveg(ecp)
!#else
!      call entcell_update_shc(ecp)
!#endif

      !if (YEAR_FLAG.eq.0) call ent_GISS_init(entcellarray,im,jm,jday,year)
      !!!### REORGANIZE WTIH ent_prog.f ####!!!

      end subroutine entcell_vegupdate

!************************************************************************
      subroutine entcell_update_lai_poolslitter( ecp,
     i    laidata,init)!,mixed_VEG)
!@sum sets prescribed LAI over the cell and update  C_fol, C_froot, senescefrac
      use phenology, only : litter_patch
      implicit none
      type(entcelltype) :: ecp  !?pointer? ##
      real*8,intent(in) :: laidata(N_PFT) !@var LAI for all PFT's 
      logical,intent(in) :: init!,mixed_VEG
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch
      type(cohort), pointer :: cop !@var current cohort
      real*8 :: laipatch, lai_old,lai_new
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      integer :: i

      laipatch = 0.d0
      lai_new = 0.d0
      lai_old = 0.d0
      i = 0
!      Clossacc(:,:,:) = 0.d0  !?Moved to inside patch loop - NK
      pp => ecp%oldest      
      do while ( associated(pp) )
        Clossacc(:,:,:) = 0.d0
        laipatch = 0.d0
        cop => pp%tallest
        do while ( associated(cop) )
          i = i + 1
          lai_old = cop%LAI
!          if (.not.mixed_VEG) then !* LAI array is by PFT.
             lai_new = laidata(cop%pft)
!          else !* mixed_VEG, LAI array is in cohort order.
#ifdef MIXED_CANOPY
             lai_new = laidata(i)
             print *,i,'lai_new',lai_new
             print *,'Bug? Need above print stmt to update laipatch(??)'
             !### For some reason laipatch is zero if this print statement
             !### is not done.  Print statement can even go after 
             !### laipatch is assigned.?? -NK
!          endif
#endif
!MIXED_CANOPY
          !* Update biomass pools, senescefrac, accumulate litter
          call prescr_veglitterupdate_cohort(cop,lai_new,Clossacc,init)

          laipatch = laipatch + cop%lai
          cop => cop%shorter
        enddo
        call litter_patch(pp,Clossacc) 
        pp%LAI = laipatch
        !print *, 'laipatch',laipatch,pp%LAI
        pp => pp%younger
      enddo

      end subroutine entcell_update_lai_poolslitter


      subroutine entcell_update_height( ecp,
     i    hdata,init)!,mixed_VEG)
!@sum sets prescribed canopy height for entcell subgrid fractions.
      use ent_prescr_veg, only : popdensity
      use allometryfn, only : crown_radius_horiz_allom,crown_radius_vert
     &     ,height2dbh, allom_plant_cpools
      use phenology, only : litter_growth_cohort, litter_patch
      use ent_pfts
      type(entcelltype) :: ecp
      real*8,intent(in) :: hdata(N_PFT) !@var LAI for all PFT's 
      logical,intent(in) :: init!, mixed_VEG
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch
      type(cohort), pointer :: cop !@var current cohort
      real*8 :: cpool(N_BPOOLS)
      real*8 :: Cdead,dC_repro
      real*8 :: C_sw_old, C_hw_old,C_croot_old,C_fol_old,C_froot_old
      real*8 :: h_old
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      real*8 :: resp_growth2
      integer :: i

      cpool(:) = 0.d0
      pp => ecp%oldest      
      i = 0

      do while ( associated(pp) )
        Clossacc(:,:,:) = 0.d0
        cop => pp%tallest
        do while ( associated(cop) )
          i = i + 1
          h_old = cop%h
          C_fol_old = cop%C_fol
          C_sw_old = cop%C_sw
          C_hw_old = cop%C_hw
          C_froot_old = cop%C_froot
          C_croot_old = cop%C_croot
!          if (.not.mixed_VEG) then !hdata array is by PFT
#ifndef MIXED_CANOPY
             cop%h = hdata(cop%pft)
!          else !mixed_VEG, hdata array is by cohort order
#else  /* MIXED_CANOPY */
             cop%h = hdata(i)
!          endif
#endif /* MIXED_CANOPY */
          if (pfpar(cop%pft)%woody) then !update dbh
            cop%dbh = height2dbh(cop%pft,cop%h)
!            if (init) then !Set population density and crown geometry
!            !This routine should get called after initialization, so
!            !density should not be updated, only plant geometry.
!              cop%n = popdensity(cop%pft,cop%dbh,
!     &              alamax(cop%PFT+COVEROFFSET))
!              cop%crown_dx = crown_radius_horiz_HF(cop%pft,cop%dbh,cop%n)
              cop%crown_dx = crown_radius_horiz_allom(
     &           cop%pft,cop%h,cop%n)
              cop%crown_dy = crown_radius_vert(cop%pft
     &             ,cop%h,cop%crown_dx)
!            endif
#ifdef ENT_STANDALONE_DIAG
            print *,'pft,n,h,dbh,crown_dx,crown_dy',
     &           cop%pft,cop%n,cop%h,cop%dbh,cop%crown_dx,cop%crown_dy
#endif            
          endif

          !* Update biomass pools except for C_lab.
          call allom_plant_cpools(cop%pft, cop%lai, cop%h, 
     &         cop%dbh, cop%n, cpool )

          cop%C_fol = cpool(FOL)
          cop%C_sw = cpool(SW)
          cop%C_hw = cpool(HW)
          cop%C_froot = cpool(FR)
          cop%C_croot = cpool(CR)

          !* If h loss, dump to litter.  If h growth, withdraw from C_lab.-NK
          if (.not.init) then
             if (cop%h<h_old) then !Decrease in height -> put in litter
                Cdead = cop%C_hw + cop%C_croot 
                dC_repro = 0.d0
                call litter_growth_cohort(.true. !prescribed structural growth
     i             ,dC_repro
     i             ,C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old
     &             ,cop%C_hw-C_hw_old,cop%C_croot-C_croot_old
     &             ,Cdead
     &             ,cop,Clossacc, resp_growth2)
             elseif (cop%h>h_old) then  !Increase in height -> withdraw from C_lab
                cop%C_lab = cop%C_lab
     &               - (cop%C_fol - C_fol_old)
     &               - (cop%C_sw - C_sw_old)
     &               - (cop%C_hw - C_hw_old)
     &               - (cop%C_froot - C_froot_old)
     &               - (cop%C_croot - C_croot_old)
                !###NK - Need to add tissue growth respiration here.
             endif
          endif
          cop => cop%shorter
          print *,'After cop=>cop%shorter'
        enddo
        call litter_patch(pp, Clossacc)
        pp => pp%younger
      enddo
      !summarize_patch - called outside this routine
      end subroutine entcell_update_height


      subroutine entcell_update_albedo( ecp,
     i    albedodata)
!@sum sets prescribed albedo in vegetated patches of the cell (skips bare soil)
!@+   This subroutine assumes one cohort per patch !!!
      type(entcelltype) :: ecp
      real*8,intent(in) :: albedodata(N_BANDS,N_PFT) !@var albedo for all PFTs 
      !-----Local---------
      type(patch), pointer :: pp  !@var p current patch

      pp => ecp%oldest      
      do while ( associated(pp) )
        ! update albedo of vegetated patches only
        if ( associated(pp%tallest) ) then ! have vegetation
          if ( pp%tallest%pft > N_PFT .or.  pp%tallest%pft < 1 )
     &         call stop_model("entcell_update_albedo: bad pft",255)
          pp%albedo(1:N_BANDS) = albedodata(1:N_BANDS, pp%tallest%pft)
        endif
        pp => pp%younger
      enddo

      end subroutine entcell_update_albedo

      subroutine entcell_update_crops( ecp,
     i    cropsdata)
!@sum sets prescribed crop cover given input crop cover array.
!@+   This subroutine assumes one cohort per patch !!!
!@+  
!@+   - Calculate vegetation, v, and soil, s, carbon/area, c, for non-crop patches, k.
!@+       c_v,k + c_s,k  where k ranges over non-crop patches.
!@+   - Calculate vegetation and soil carbon/area for crop patch, cr.
!@+       c_v,cr + c_s,cr
!@+
!@+   Crops area changes by df_c: 
!@+     (if df_c>0, the crops area increase; if df_c<0, then crops area decreases)
!@+   - Decrease non-crop patch fractional areas,f_k, proportionally by crop fraction
!@+       area increase, df_cr.
!@+       f_k,new = f_k - (f_k/sum_f_k)*df_cr
!@+       s.th.  sum_k((f_k/sum_f_k)*df_c) = df_cr
!@+   - Increase crop area (or if new, insert crop patch).
!@+       f_c,new = f_c + df_cr
!@+   - If df_c>0:
!@+     Update crop soil carbon to be sum of acquired non-crop patch soils:
!@+       I.e. soil carbon from non-crop patches is merged into the new crop area.
!@+       (sum each soil carbon pool, p, then divide by new crop area).
!@+       c_s,c,p,new = { sum_k(c_s,k,p * df_cr*f_k/sum_f_k) + c_s,cr,p*f_c }/f_c,new
!@+   - If df_c<0:
!@+     Update non-crop patch soil carbon to merge in crop patch increment soil carbon.
!@+       (sum each soil carbon pool, p, then divide by new patch area).
!@+       c_s,k,p,new = { c_s,k,p*f_k + c_s,c,p*f_k/sum_f_k*df_cr }/f_k,new
!@+   - Net carbon to atmosphere is:
!@+       C_to_atm = loss from non-crop patch veg - gain from crop patch increment
!@+                = sum_k( c_v,k * df_cr*f_k/sum_f_k) - c_v,cr*df_cr
!@+
!@+   - Where to put a net loss in biomass carbon (all soil carbon stays in the soil):
!@+       Option 1:  Put it all as litter into the crop patch.
!@+                  This will require partitioning the lost biomass pools into litter
!@+                  pools with relevant ligning content and C:N ratios.
!@+       Option 2:  Put it into a slow-release pool to send it all to the atmosphere
!@+                  over a few months or the next year.
!@+       Option 1 is more strictly correct, but involves a little more computation of
!@+       litter pools, and allows also an option to sequester, e.g. wood harvest.
!@+
!@+   - Where to put a net gain in biomass carbon:
!@+      Option 1:   There's only one option:  to suck it all from the atmosphere at
!@+                  once.  Given the net sum of gains and losses over the entire
!@+                  land surface, hopefully this shouldn't cause a dramatic dent in
!@+                  atmospheric CO2 over one time step.

      use entcells, only : entcell_extract_pfts, summarize_entcell
      use patches, only : patch_has_pft, patch_split, patch_merge,
     &     delete_patch, patch_set_pft, patch_delete_cohort
      type(entcelltype) :: ecp
      real*8,intent(in) :: cropsdata !@var crop cover fraction.
      !-----Local---------
      integer, parameter :: PFT_CROPS = 8
      type(patch), pointer :: pp, pp_crops, pp_tmp, pp_end
      real*8 :: crops_old, dcrops, dfr
      real*8 :: vfraction(12)

!#define UNFINISHED_CROPS_CODE
#ifdef UNFINISHED_CROPS_CODE


      ! for debug 

      call entcell_extract_pfts(ecp, vfraction)
      print *,"updating crops, old vfraction:", vfraction

      ! find fraction of old crops
      crops_old = 0.d0
      pp => ecp%oldest      
      do while ( associated(pp) )
        if ( patch_has_pft(pp, PFT_CROPS ) )
     &       crops_old = crops_old + pp%area
        pp => pp%younger
      enddo
      dcrops = cropsdata - crops_old
      if ( abs(dcrops) < 1.d-2 ) return ! ignore very small changes 1d-6


      call summarize_entcell(ecp )
      write(997,*) "before"
      write(997,*) "ecp%C_fol   ", ecp%C_fol 
      write(997,*) "ecp%N_fol   ", ecp%N_fol 
      write(997,*) "ecp%C_w     ", ecp%C_w 
      write(997,*) "ecp%N_w     ", ecp%N_w 
      write(997,*) "ecp%C_lab   ", ecp%C_lab 
      write(997,*) "ecp%N_lab   ", ecp%N_lab 
      write(997,*) "ecp%C_froot ", ecp%C_froot
      write(997,*) "ecp%N_froot ", ecp%N_froot
      write(997,*) "ecp%C_root  ", ecp%C_root 
      write(997,*) "ecp%N_root  ", ecp%N_root 
      write(997,*) "ecp%nm ", ecp%nm
      write(997,*) "ecp%Ntot ", ecp%Ntot
      write(997,*) "ecp%LMA ", ecp%LMA
      write(997,*) "ecp%h ", ecp%h
      write(997,*) "ecp%LAI ", ecp%LAI
      write(997,*) "ecp%LAIpft ", ecp%LAIpft
      write(997,*) "ecp%fracroot ", ecp%fracroot

      call summarize_entcell(ecp )
      write(6,*) __LINE__,"ecp%C_froot ", ecp%C_froot
      call entcell_extract_pfts(ecp, vfraction)
      write(6,*) "sum= ", sum(vfraction(:))


      if ( crops_old > 1.d0 - 1.d-6 .or.
     &     cropsdata > 1.d0 - 1.d-6) call stop_model(
     &     "Fractions with 100% of crops are not supported",255)
      dfr = dcrops/(1.d0 - crops_old)
      if ( dcrops > 0 ) then
        ! split each patch to create crops patch of area = pp%area*dfr 
        print *,"here ",__FILE__,__LINE__
        pp_end => ecp%youngest
        pp => ecp%oldest
        do while ( associated(pp) )
          if ( .not. patch_has_pft(pp, PFT_CROPS ) ) then
            call patch_split(pp, pp%area*dfr, pp_tmp)
      !debug...........
      call summarize_entcell(ecp )
      write(6,*) __LINE__,"ecp%C_froot ", ecp%C_froot
      call entcell_extract_pfts(ecp, vfraction)
      write(6,*) "sum= ", sum(vfraction(:))
      !............
            call patch_set_pft(pp_tmp, PFT_CROPS)
          endif
          if ( associated(pp, pp_end) ) exit
          pp => pp%younger
        enddo
      !debug...........
      call summarize_entcell(ecp )
      write(6,*) __LINE__,"ecp%C_froot ", ecp%C_froot
      call entcell_extract_pfts(ecp, vfraction)
      write(6,*) "sum= ", sum(vfraction(:))
      !............
        ! find first crops patch
        print *,"here ",__FILE__,__LINE__
        pp => ecp%oldest      
        do while ( associated(pp) )
          if ( patch_has_pft(pp, PFT_CROPS ) ) exit
          pp => pp%younger
        enddo
        pp_crops => pp
        ! merge all crops patches into one patch
        print *,"here ",__FILE__,__LINE__
        pp => ecp%oldest      
        do while ( associated(pp) )
          pp_tmp => pp%younger
          if ( patch_has_pft(pp, PFT_CROPS )
     &         .and. .not. associated(pp, pp_crops) )
     &         call patch_merge(pp_crops, pp)
          pp => pp_tmp
        enddo
      !debug...........
      call summarize_entcell(ecp )
      write(6,*) __LINE__,"ecp%C_froot ", ecp%C_froot
      call entcell_extract_pfts(ecp, vfraction)
      write(6,*) "sum= ", sum(vfraction(:))
      !............
        print *,"here" ,__FILE__,__LINE__
      else
        ! find first crops patch
        print *,"here ",__FILE__,__LINE__
        pp => ecp%oldest      
        do while ( associated(pp) )
          if ( patch_has_pft(pp, PFT_CROPS ) ) exit
          pp => pp%younger
        enddo
        print *,"here ",__FILE__,__LINE__
        pp_crops => pp
        if ( abs(crops_old - pp_crops%area) > 1.d-6 ) call stop_model(
     &       "More than 1 crops patches per cell not supported",255)
        ! split crops patch into pieces and merge them with other patches
        print *,"here ",__FILE__,__LINE__
      !debug...........
      call summarize_entcell(ecp )
      write(6,*) __LINE__,"ecp%C_froot ", ecp%C_froot
      call entcell_extract_pfts(ecp, vfraction)
      write(6,*) "sum= ", sum(vfraction(:))
      !............
        pp => ecp%oldest      
        do while ( associated(pp) )
          if ( .not. patch_has_pft(pp, PFT_CROPS ) ) then
            call patch_split(pp_crops, -pp%area*dfr, pp_tmp)
            ! hack to deal with bare land (do not conserve C here)
            if ( .not. associated( pp%tallest ) ) then
              call patch_delete_cohort(pp_tmp,pp_tmp%tallest)
            else
              call patch_set_pft(pp_tmp, pp%tallest%pft)
            endif
            call patch_merge(pp, pp_tmp)
          endif
          pp => pp%younger
        enddo
      !debug...........
      call summarize_entcell(ecp )
      write(6,*) __LINE__,"ecp%C_froot ", ecp%C_froot
      call entcell_extract_pfts(ecp, vfraction)
      write(6,*) "sum= ", sum(vfraction(:))
      !............
        print *,"here ",__FILE__,__LINE__
        if ( pp_crops%area < 1.d-6 ) call delete_patch(ecp, pp_crops)
      endif

        print *,"here ",__FILE__,__LINE__
      call entcell_extract_pfts(ecp, vfraction)
      print *,"updating crops, new vfraction:", vfraction

      call summarize_entcell(ecp )
      write(997,*) "after"
      write(997,*) "before"
      write(997,*) "ecp%C_fol   ", ecp%C_fol 
      write(997,*) "ecp%N_fol   ", ecp%N_fol 
      write(997,*) "ecp%C_w     ", ecp%C_w 
      write(997,*) "ecp%N_w     ", ecp%N_w
      write(997,*) "ecp%C_lab   ", ecp%C_lab 
      write(997,*) "ecp%N_lab   ", ecp%N_lab 
      write(997,*) "ecp%C_froot ", ecp%C_froot
      write(997,*) "ecp%N_froot ", ecp%N_froot
      write(997,*) "ecp%C_root  ", ecp%C_root 
      write(997,*) "ecp%N_root  ", ecp%N_root 
      write(997,*) "ecp%nm ", ecp%nm
      write(997,*) "ecp%Ntot ", ecp%Ntot
      write(997,*) "ecp%LMA ", ecp%LMA
      write(997,*) "ecp%h ", ecp%h
      write(997,*) "ecp%LAI ", ecp%LAI
      write(997,*) "ecp%LAIpft ", ecp%LAIpft
      write(997,*) "ecp%fracroot ", ecp%fracroot


#endif

      end subroutine entcell_update_crops


!******************************************************************

      subroutine prescr_phenology(jday,hemi,pp,do_giss_lai)
!@sum Prescribed (Matthews 1983) phenology (LAI seasonality)
!@+   and associated litterfall.
!@+   DAILY TIME STEP.
!@+   Calculate new LAI, biomass poos, and senescefrac 
!@+      for given jday, for prescribed vegetation structure.
      use ent_pfts
      use ent_prescr_veg, only : prescr_calc_lai
     &     ,prescr_veg_albedo
      use phenology, only : litter_patch
      implicit none
      integer,intent(in) :: jday !Day of year.
      integer,intent(in) :: hemi !@var hemi -1: S.hemisphere, 1: N.hemisphere
      type(patch),pointer :: pp
      logical, intent(in) :: do_giss_lai
      !-------local-----
      type(cohort),pointer :: cop
      real*8 :: laipatch
      real*8 :: cpool(N_BPOOLS)
      real*8 :: lai_new
      real*8 :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.

      if (ASSOCIATED(pp)) then
        !* LAI *AND* BIOMASS - carbon pools *!
        laipatch = 0.d0         !Initialize for summing
        cpool(:) = 0.d0
        Clossacc(:,:,:) = 0.d0
        
        cop => pp%tallest
        do while (ASSOCIATED(cop))
          if ( do_giss_lai ) then  !This if statement is now redundant
            lai_new = prescr_calc_lai(cop%pft+COVEROFFSET, jday, hemi)
            call prescr_veglitterupdate_cohort(cop,
     &           lai_new,Clossacc,.false.)
            laipatch = laipatch + cop%LAI
         endif
          cop => cop%shorter
        end do
        call litter_patch(pp,Clossacc) !Update Tpools following litter accumulation. Daily time step
        pp%LAI = laipatch

        !* ALBEDO *! - Moved these lines to entcell_vegupdate.
!        if ( ASSOCIATED(pp%tallest) ) then ! update if have vegetation
!          call prescr_veg_albedo(hemi, pp%tallest%pft+COVEROFFSET, 
!     &         jday, pp%albedo)
!        endif
      endif
      end subroutine prescr_phenology

!******************************************************************************
      subroutine prescr_veglitterupdate_cohort(
     &     cop,lai_new,Clossacc,init)
!@sum Given new LAI, update cohort biomass pools, litter, senescefrac.
      use phenology, only : litter_cohort
      use allometryfn, only :allom_plant_cpools
      implicit none
      type(cohort),pointer :: cop
      real*8,intent(in) :: lai_new
      real*8,intent(inout) :: Clossacc(PTRACE,NPOOLS,N_CASA_LAYERS) !Litter accumulator.
      logical,intent(in) :: init
       !-------local-----
      real*8 :: cpool(N_BPOOLS)
      real*8 :: lai_old
      real*8 :: C_fol_old,C_froot_old,C_hw_old,C_croot_old,C_sw_old
      real*8 :: C_lab_old, Clossacc_old, dC
!#ifdef DEBUG
      integer :: i
      real*8 :: Csum
!#endif

      lai_old = cop%LAI
      C_fol_old = cop%C_fol
      C_froot_old = cop%C_froot
      C_hw_old = cop%C_hw
      C_sw_old = cop%C_sw
      C_croot_old = cop%C_croot
      cop%LAI = lai_new
#ifdef CHECK_CARBON_CONSERVATION
      C_lab_old = cop%C_lab
      Clossacc_old = sum(Clossacc(CARBON,SURFMET:CWD,1))
#endif
 
      !* Get new biomass pools except for C_lab.
      call allom_plant_cpools(cop%pft, cop%lai, cop%h, 
     &     cop%dbh, cop%n, cpool )
      cop%C_fol = cpool(FOL)
      cop%C_sw = cpool(SW)
      cop%C_hw = cpool(HW)
      cop%C_froot = cpool(FR)
      cop%C_croot = cpool(CR)

      !* Update litter, C_lab, senescefrac
      if (.not.init) then
        call litter_cohort(SDAY,
     i       C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
     &       cop,Clossacc)
!        write(992,*)C_fol_old,C_froot_old,C_hw_old,C_sw_old,C_croot_old,
!     &       cop%C_fol,cop%C_froot,cop%C_hw,cop%C_sw,cop%C_croot
      endif

#ifdef CHECK_CARBON_CONSERVATION
      dC = cop%C_fol-C_fol_old + cop%C_sw-C_sw_old +
     &     cop%C_hw-C_hw_old + cop%C_froot-C_froot_old +
     &     cop%C_croot-C_croot_old + cop%C_lab-C_lab_old +
     &     (sum(Clossacc(CARBON,SURFMET:CWD,1))-Clossacc_old)/cop%n

      if( abs(dC) > 1.d-8 ) then
        write(99,*) "litter", dC,
     &       (sum(Clossacc(CARBON,SURFMET:CWD,1))-Clossacc_old)/cop%n,
     &       cop%C_fol-C_fol_old, cop%C_sw-C_sw_old, cop%C_hw-C_hw_old,
     &       cop%C_froot-C_froot_old, cop%C_croot-C_croot_old,
     &       cop%C_lab-C_lab_old

      endif

!!!! HACK to conserve carbon (add error to C_lab)
!!!!      cop%C_lab = cop%C_lab - dC  !! seems like not needed any more ...
#endif

        !*## DEBUG ##*!
!#ifdef DEBUG
      Csum = 0.d0
      do i=1,NPOOLS
        Csum = Csum + Clossacc(CARBON,i,1)
      enddo
      !write(999,*) Csum
      cop%N_up = Csum !####### HACK ## TEMPORARY USE OF UNUSED VARIABLE ###-NK
!#endif
        
      
      cop%Ntot = cop%nm * cop%LAI !This should eventually go into N allocation routine if dynamic nm.
      cop%LAI = lai_new  
 
      end subroutine prescr_veglitterupdate_cohort
!******************************************************************************

      end module ent_prescribed_updates
