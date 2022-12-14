      module patches
!@sum Routines to organize patches in an entcell.

      use ent_types

      implicit none


      contains
      !*********************************************************************
      subroutine insert_patch(gp, area, soil_type)
!@sum insert_patch Insert patch at youngest end of patch list.
!@+   Blank patch with no cohorts.
      implicit none
      type(entcelltype), target :: gp
      real*8, intent(in) :: area
      integer, intent(in) :: soil_type
      !----
      type(patch),  pointer :: pp !Would be good to return pointer to patch.

      !* If !ASSOCIATED(gp) ERROR *!

      ! create patch
      call patch_construct(pp, gp, area, soil_type)

      !* If there exist patches, then insert.  If no patches, then allocate.
      if (ASSOCIATED(gp%youngest)) then
        !allocate(gp%youngest%younger)
        gp%youngest%younger => pp
        pp%older => gp%youngest
        gp%youngest => pp
      else !first patch on entcell
        gp%oldest => pp
        gp%youngest => pp
      end if

      end subroutine insert_patch

      !*********************************************************************
      subroutine assign_patch(pp,
     &     Ci_ini, CNC_ini, pft, Tpool_ini)
!@sum assign_patch   Initialize some canopy and soil variables for a patch.
      use ent_const  
      
      !Eventually may want to include all patch variables as optional vars.
      type(patch),pointer :: pp
      integer, intent(in) :: pft
      real*8 :: Ci_ini, CNC_ini
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS) ::
     &         Tpool_ini  !in g/m2 -PK
      !----Local-------
      integer :: i, n

      pp%Ci = Ci_ini
      pp%GCANOPY = CNC_ini

      !Assign soil (dead) pools.
      call assign_patch_soilcarbon(pp, pft, Tpool_ini)
!      do n=1,N_CASA_LAYERS 
!       do i=NLIVE+1,NPOOLS
!        pp%Tpool(CARBON,i,n) = Tpool_ini(pft,CARBON,i-NLIVE,n)  
!       end do
!      end do
c!#else
c!      pp%Tpool(:,:,:) = 0.d0  
c!#endif

#ifdef OFFLINE  
!      print*, 'soil C pools: ', pp%Tpool(:,CARBON,:,:)  !optional test -PK
#endif
      end subroutine assign_patch

      !*********************************************************************
      subroutine assign_patch_soilcarbon(pp, pft, Tpool_ini)
!@sum assign_patch_soilcarbon  Initialize.
      !Eventually may want to include all patch variables as optional vars.
      !in assign_patch and get rid of this subroutine
     
      use ent_const  
      
      type(patch),pointer :: pp
      integer, intent(in) :: pft
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS) ::
     &         Tpool_ini  !gC/m2, soil pools only
      !----Local-------
      integer :: i, n

      !Assign soil (dead) pools.
      do n=1,N_CASA_LAYERS 
       do i=NLIVE+1,NPOOLS
        pp%Tpool(CARBON,i,n) = Tpool_ini(pft,CARBON,i-NLIVE,n)  
       end do
      end do
!      pp%Tpool(CARBON,(NLIVE+1):NPOOLS,1:N_CASA_LAYERS) =
!     &     Tpool_init(pft,CARBON,1:(NPOOLS-NLIVE),1:N_CASA_LAYERS)

      end subroutine assign_patch_soilcarbon
      !*********************************************************************
      subroutine delete_patch(gp, pp)
!@sum Delete patch pointed to by pp
      !* NOTE:  THIS DOES NOT AUTOMATICALLY UPDATE ENTCELL SUMMMARY VALUES
      implicit none
      type(entcelltype) :: gp
      type(patch), pointer :: pp

      if (.NOT.ASSOCIATED(pp%older)) then !is oldest
        gp%oldest => gp%oldest%younger
        nullify(gp%oldest%older)
      else if(.NOT.ASSOCIATED(pp%younger)) then !is youngest
        gp%youngest => gp%youngest%older
        nullify(gp%youngest%younger)
      else !pp points somewhere in the middle of patch list
        pp%older%younger => pp%younger
        pp%younger%older => pp%older
      end if
      call patch_destruct(pp)

      end subroutine delete_patch

 !**************************************************************************
      real*8 function shc_patch(pp) Result(shc)
!@sum Return GISS GCM specific heat capacity for patch, but with subgrid
!@+   calculation of shc.  See comments in entcell_update_shc_mosaicveg.
!@+   Awkward legacy code.
      use ent_const
      use ent_pfts, only: COVEROFFSET, alamax, alamin
      !use ent_prescr_veg, only : GISS_shc
      use allometryfn, only : GISS_shc
      implicit none
      type(patch),pointer :: pp
      !-----
      type(cohort),pointer :: cop
      real*8 :: lai
      integer :: anum

      shc = 0.d0
      lai = 0.d0
      cop => pp%shortest
      do while (associated(cop))
         anum = cop%pft + COVEROFFSET
         lai = lai + .5d0*(alamax(anum) + alamin(anum)) 
         !shc = shc + GISS_shc(cop%pft) !Alternative way.
         cop => cop%shorter
      enddo

      !shc = (.010d0+.002d0*lai+.001d0*lai**2)*shw*rhow
      shc = GISS_shc(lai)  !Note: original GISS shc is at cell level, not patch.

      end function shc_patch
!**************************************************************************

      real*8 function laimean_annual_patch(pp) Result(laires)
!@sum Calculate lai annual mean ONLY for old R&A shc scheme.
      !NYK:  This function called by entcell_update_shc/_mosaicveg,
      !      made just to preserve R&A shc scheme, which uses 
      !      mean annual entcell LAI, constant for Matthews veg,
      !      but which needs to be replaced by canopy biomass and water, and
      !      should seasonally vary!
      use allometryfn, only : Cfol_fn
      use ent_pfts, only : pfpar
      implicit none
      type(patch),pointer :: pp
      !-----
      type(cohort),pointer :: cop
      integer :: pft
      real*8 :: laimax

      laires = 0.d0
      cop => pp%tallest
      do while (ASSOCIATED(cop))
         pft = cop%pft
         laimax = .001d0 * Cfol_fn(pft,cop%dbh,cop%h)
     &        * pfpar(pft)%sla
         if (pfpar(pft)%phenotype.eq.EVERGREEN) then
            laires = laires + 0.8d0*laimax  !Arbitrary guess, annual variation
         else !cold/drought/deciduous & annual
            laires = laires + 0.5d0*laimax
            !May need to set minimum of 1.5 for stability, but try actual.
         endif
         cop => cop%shorter
      enddo

      end function laimean_annual_patch


      

      subroutine summarize_patch(pp)
!@sum summarize_patch Calculate patch-level summary values of cohort pools.
!@+   Intensive properties (e.g. geometry, LMA) are averages weighted by
!@+   total number of individuals (may want to do by biomass of cohort)
!@+   Extensive properties (e.g. biomass, Ntot) are totals per m2 ground
      use cohorts, only: calc_CASArootfrac  !PK 7/07
      !use canopyrad
      implicit none
      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: cop
      real*8 :: nc, nsum  !density, sum
      integer :: ia  !array index
      integer :: pft
      real*8 :: fracrootCASA(N_CASA_LAYERS)  !to map fracroot to fracrootCASA
      real*8 :: gCindiv_to_kgCm2

      !* Zero out cohort summary variables *!
      call zero_patch_cohortsum(pp)
      pp%Tpool(:,1:NLIVE,:) = 0.d0 !###

      if ( .not. associated(pp%tallest) ) return ! no cohorts in this patch

      do ia=1,N_COVERTYPES
        pp%LAIpft(ia) = 0.d0
      end do

      nsum = 0
      cop => pp%tallest
      do while(ASSOCIATED(cop))
        pft = cop%pft
        
      !assign root fractions for CASA layers -PK 11/06 
      call calc_CASArootfrac(cop%fracroot,fracrootCASA)

        !*- - - - - - COHORT SUMMARY VARIABLES - - - - - - - - - -
        nc = cop%n  
        nsum = nsum + nc !Number of individuals for summing.

        !* PFT PARAMETERS
        ! Only need to index array of pftypes.

        !* NITROGEN status and LEAF */
        pp%nm = pp%nm + cop%nm*cop%LAI !wtd average
        pp%Ntot = pp%Ntot + cop%Ntot  !Total
        pp%LAI = pp%LAI + cop%LAI  !Total
        pp%LAIpft(pft) = pp%LAIpft(pft) + cop%LAI  !Total
        pp%LMA = pp%LMA + cop%LMA * cop%LAI  !wtd avg

         !* GEOMETRY - WEIGHTED AVERAGES
        pp%h = pp%h + cop%h * nc  !wtd avg
        pp%crown_dx = pp%crown_dx + cop%crown_dx * nc !wtd avg
        pp%crown_dy = pp%crown_dy + cop%crown_dy * nc !wtd avg
        pp%clump = pp%clump + cop%clump * nc !wtd avg - ##SHOULD COME FROM RADIATIVE TRANSFER
        !* Do fracroot(:) outside this loop, below.

         !* BIOMASS POOLS - TOTALS
        gCindiv_to_kgCm2 = 1d-3 * nc !Convert gC/individ to kg-C/m^2
        pp%C_fol = pp%C_fol + cop%C_fol * gCindiv_to_kgCm2
        pp%N_fol = pp%N_fol + cop%N_fol * gCindiv_to_kgCm2
        pp%C_w = pp%C_w + (cop%C_sw  + cop%C_hw) * gCindiv_to_kgCm2
        pp%N_w = pp%N_w + (cop%N_sw  + cop%N_hw) * gCindiv_to_kgCm2
        pp%C_lab = pp%C_lab + cop%C_lab * gCindiv_to_kgCm2
        pp%N_lab = pp%N_lab + cop%N_lab * gCindiv_to_kgCm2
        pp%C_froot = pp%C_froot + (cop%C_froot) * gCindiv_to_kgCm2
        pp%N_froot = pp%N_froot + (cop%N_froot) * gCindiv_to_kgCm2
        pp%C_root = pp%C_root + 
     &       (cop%C_froot + cop%C_croot) * gCindiv_to_kgCm2
        pp%N_root = pp%N_root + 
     &       (cop%N_froot + cop%N_croot) * gCindiv_to_kgCm2

         !* FLUXES - TOTALS
        pp%Ci = pp%Ci + cop%Ci * cop%LAI !wtd average
        pp%GCANOPY = pp%GCANOPY + cop%GCANOPY !total
        pp%GPP = pp%GPP + cop%GPP     !total
        pp%IPP = pp%IPP + cop%IPP     !total
        pp%NPP = pp%NPP + cop%NPP     !total
        pp%R_auto = pp%R_auto + cop%R_auto        !total
        pp%R_root = pp%R_root + cop%R_root        !total 
        pp%N_up = pp%N_up + cop%N_up  !total
!        pp%C_litter = pp%C_litter + cop%C_litter
!        pp%N_litter = pp%N_litter + cop%N_litter  
!        pp%C_to_Nfix = pp%C_to_Nfix + cop%C_to_Nfix

        pp%betad = pp%betad + cop%stressH2O*cop%LAI
        do ia=1,N_DEPTH
          pp%betadl(ia) = pp%betadl(ia) + cop%stressH2Ol(ia)*cop%LAI
        end do
        !*- - - - - end cohort summary variables - - - - - - - - - - - - -

        !* CASA pools * These are redundant but not 1-1 with the cohort pools.
       do ia=1,N_CASA_LAYERS  !loop over all CASA layers -PK 7/07
        if (ia.eq.1) then  !only top CASA layer has leaves & wood
         pp%Tpool(CARBON,LEAF,ia) = pp%Tpool(CARBON,LEAF,ia)
     &        + cop%C_fol * cop%n    
         pp%Tpool(CARBON,WOOD,ia) = pp%Tpool(CARBON,WOOD,ia)  !note: phenology.f uses C_croot instead of C_sw -PK 7/07
     &       + (cop%C_sw + cop%C_hw + cop%C_croot)* cop%n 
        else    
         pp%Tpool(CARBON,LEAF,ia) = 0.d0  
         pp%Tpool(CARBON,WOOD,ia) = 0.d0
        end if 
         pp%Tpool(CARBON,FROOT,ia) = pp%Tpool(CARBON,FROOT,ia)
     &        + fracrootCASA(ia)*cop%C_froot * cop%n 
       end do  !N_CASA_LAYERS
        !* Tpool(NITROGEN,:,:) gets updated in casa_bgfluxes.f

       pp%C_total = pp%C_total + cop%C_total
       pp%C_growth = pp%C_growth + cop%C_growth
       cop => cop%shorter
      end do    !loop through cohorts

      !* ------- DO AVERAGES ----------------------------------------------*!

      !* CHECK IF NC = 0. If ASSOCIATE(pp%tallest) then should not be a problem.*!
!      pp%n = nsum                          !* Total individuals *!

      !* NITROGEN status and LEAF */
      if (pp%LAI.gt.0.d0) then
        pp%nm = pp%nm / pp%LAI  !wtd average
        pp%LMA = pp%LMA / pp%LAI !wtd avg
      else
        pp%nm = 0.d0
        pp%LMA = 0.d0
      endif

      !* GEOMETRY - WEIGHTED AVERAGES
      if ( nsum > 0.d0 ) then
        pp%h = pp%h / nsum      !wtd avg
        pp%crown_dx = pp%crown_dx / nsum !wtd avg
        pp%crown_dy = pp%crown_dy / nsum !wtd avg
        pp%clump = pp%clump / nsum !wtd avg - ##SHOULD COME FROM RADIATIVE TRANSFER
        CALL  sum_roots_cohorts2patch(pp) !froot and C_froot
      endif

      !* FLUXES 
      if (pp%LAI.gt.0.d0) then
        pp%Ci = pp%Ci / pp%LAI  !wtd average
        
        pp%betad = pp%betad/pp%LAI
        pp%betadl = pp%betadl/pp%LAI
      else
        pp%Ci = 0.d0
        pp%betad = 0.d0
        pp%betadl(:) = 0.0d0
      endif
      
      !* Flux variables for GCM/EWB - patch total
      !pp%z0 =0.d0               !## Dummy ##!
      !##### call get_patchalbedo(jday,pp)######
      !##### Calculate TRANS_SW ################

      !* Variables calculated by GCM/EWB - downscaled from grid cell
      !Soilmoist(:)        !Calculated by GCM/EWB
      !pp%N_deposit = 0.0  !Dummy until we have N cycle
      
      !* Variables for biophysics and biogeochemistry
      !nullify(pp%crad)     !Dummy until we have GORT clumping.

      !* Disturbance values - UPDATE WHEN WE INCLUDE DISTURBANCE
      !pp%fuel = 0.d0          !## Dummy ##!
      !pp%ignition_rate =0.d0  !## Dummy ##!
      !pp%lambda1(T_SUB) = 0.d0!## Dummy ##!
      !pp%disturbance_rate(N_DIST_TYPES)=0.d0 !## Dummy ##!

      end subroutine summarize_patch

      !*********************************************************************

      subroutine zero_patch_cohortsum(pp)
!@sum Zero only the patch's cohort summary variables
      implicit none
      type(patch),pointer :: pp

      pp%nm = 0.d0
      pp%Ntot = 0.d0
      pp%LAI = 0.d0
      pp%LMA = 0.d0
      pp%h = 0.d0
      pp%crown_dx = 0.d0
      pp%crown_dy = 0.d0
      pp%clump = 0.d0
      pp%fracroot(:) = 0.d0
      pp%C_fol = 0.d0
      pp%N_fol = 0.d0
      pp%C_w = 0.d0
      pp%N_w = 0.d0
      pp%C_lab = 0.d0
      pp%N_lab = 0.d0
      pp%C_froot = 0.d0
      pp%N_froot = 0.d0
      pp%C_root = 0.d0
      pp%N_root = 0.d0
      pp%Ci = 0.d0
      pp%GCANOPY = 0.d0
      pp%GPP = 0.d0
      pp%IPP = 0.d0
      pp%NPP = 0.d0
      pp%R_auto = 0.d0
      pp%R_root = 0.d0  !PK 5/15/07
      pp%N_up = 0.d0
      pp%betad = 0.d0
      pp%betadl(:) = 0.d0
      ! diags and hacks
      pp%C_total = 0.d0
      pp%C_growth = 0.d0
      end subroutine zero_patch_cohortsum
      !*********************************************************************

      subroutine zero_patch(pp)
!@sum Initialize patch, zeroing variables.
      implicit none
      type(patch),pointer :: pp
      !----Local----

            pp%age=0.d0
            pp%area = 0.d0
            pp%Reproduction(:) = 0.d0

            !* Structural variables *!
            call zero_patch_cohortsum(pp)

            !* Flux variables for GCM/EWB - patch total
            pp%z0 = 0.0d0    !## Get GISS z0 ######!
            pp%albedo = 0.0d0!## Get GISS albveg ##!
            pp%TRANS_SW = 1.0d0 !## Calculate for zero veg.##!
            pp%GCANOPY = 0.d0 !
            pp%CO2flux = 0.d0 
            pp%Ci = 0.0127D0  !Initial value not zero.
            pp%betad = 0.d0
            pp%betadl = 0.d0
            pp%Soil_resp = 0.d0

            !* Variables calculated by GCM/EWB - downscaled from grid cell
            pp%Soilmoist(:) = 0.0d0 !## Get GISS soil moisture layers ##!
!            pp%N_deposit = 0.d0

            !* Variables for biophysics and biogeochemistry
            !pp%crad%##### = !## Get GORT canopy radiation params ##!
            ! initializing it at least to something 
            nullify(pp%crad%heights)
            nullify(pp%crad%lai)
            pp%crad%gortclump = 0.d0
            nullify(pp%crad%crad_heights)

            !* Disturbance values
            pp%fuel = 0.d0      !## Dummy ##!
            pp%ignition_rate =0.d0 !## Dummy ##!
            pp%lambda1 = 0.d0 !## Dummy ##!
            pp%disturbance_rate = 0.d0 !## Dummy ##!

            !* DIAGNOSTIC SUMMARIES
            !* Biomass pools - patch total
            pp%LAIpft(:) = 0.d0
            pp%Tpool(:,:,:) = 0.d0

            !* Soil type
            pp%soil_type = -1    ! set to undefined soil type (maybe use -1?)

#ifdef NEWDIAG
            pp%plant_ag_Cp = 0.d0 !## Dummy ##!
            pp%plant_bg_Cp = 0.d0 !## Dummy ##!
            pp%plant_ag_Np = 0.d0 !## Dummy ##!
            pp%plant_bg_Np = 0.d0 !## Dummy ##!

            !* Biomass pools - by pft
            pp%plant_ag_C(pnum) = 0.d0 !## Dummy ##!
            pp%plant_bg_C(pnum) = 0.d0 !## Dummy ##!
            pp%plant_ag_N(pnum) = 0.d0 !## Dummy ##!
            pp%plant_bg_N(pnum) = 0.d0 !## Dummy ##!

            !* Soil pools - patch total
            pp%REW = 0.0!## CALCULATE FROM GISS SOIL MOISTURE
            pp%soil_labile_C = 0.d0 !## Dummy ##!
            pp%soil_slow_C = 0.d0 !## Dummy ##!
            pp%soil_labile_N = 0.d0 !## Dummy ##!
            pp%soil_slow_N = 0.d0 !## Dummy ##!
            pp%mineral_N = 0.d0 !## Dummy ##!

            !* Rates of change - patch total
            pp%dadt = 0.d0 !## Dummy ##!
            pp%dpdt = 0.d0 !## Dummy ##!
            pp%dwdt = 0.d0 !## Dummy ##!
#endif

            !* Activity diagnostics - can be summed by month, year, etc.

      end subroutine zero_patch

      !*********************************************************************
      !*********************************************************************
      subroutine sum_roots_cohorts2patch(pp)
!@sum Calculate patch-level depth- and mass-weighted average
!@sum of fine roots depth fractions, fracroot.
      implicit none
      type(patch),pointer :: pp
      !-----Local variables-------
      type(cohort),pointer :: cop
      integer :: n
      real*8 :: fracroot(N_DEPTH)
      real*8 :: frootC_total

      !Initialize
      do n=1,N_DEPTH  
        fracroot(n) = 0.0
      end do 
      frootC_total = 0.0

      cop => pp%tallest
      do while(ASSOCIATED(cop)) 
        frootC_total = frootC_total + cop%n*cop%C_froot !Mass sum per area
        !* NYK - haven't set up fracroot yet.
        do n=1,N_DEPTH
          !fracroot(n) = fracroot(n) + cop%n*cop%C_froot*cop%fracroot(n) !Mass wtd avg.
          fracroot(n) = fracroot(n) + cop%n*cop%fracroot(n) !Population wtd avg.##HACK
        end do
        cop => cop%shorter
      end do
      if (frootC_total > 0.0) then
        pp%fracroot = fracroot/frootC_total
      end if
!      pp%C_froot = frootC_total !Already calculated outside of this routine.

      end subroutine sum_roots_cohorts2patch
      !*********************************************************************

      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      !*********************************************************************
      
      subroutine reorganize_patches(entcell)
!@sum Place holder.
      implicit none
      type(entcelltype) :: entcell


      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !Not needed in GISS replication test.
      
      end subroutine reorganize_patches


      !*********************************************************************

      subroutine patch_construct(pp, parent_entcell, area, soil_type)
!@sum patch_construct  Allocate memory for a new patch and point to entcell
!@+   but nullify older and younger and cohort pointers.
      use cohorts, only : cohort_construct
      implicit none
      type(patch), pointer :: pp
      type(entcelltype), intent(in), target :: parent_entcell
      real*8, intent(in) :: area
      integer, intent(in) :: soil_type
      !---

      if ( area > 1.d0 ) then 
        print *, "patch area > 1, area=",area
        call stop_model("patch_construct: area >1",255)
      endif
      ! allocate memory
      allocate( pp )
      allocate( pp%LAIpft(N_COVERTYPES) )  !remove -- have made into scalar -PK 6/28/06
      allocate( pp%betadl(N_DEPTH) )
      allocate( pp%fracroot(N_DEPTH) )
!      allocate( pp%LAIpft(N_COVERTYPES))
      allocate( pp%Reproduction(N_PFT) ) 

      ! set pointers
      pp%cellptr => parent_entcell
      nullify( pp%older )
      nullify( pp%younger )
      nullify( pp%tallest )
      nullify( pp%shortest )

      ! set variables
      call zero_patch(pp)

      pp%area = area
      pp%soil_type = soil_type
      
      end subroutine patch_construct

      !*********************************************************************

      subroutine patch_destruct(pp)
!@sum patch_destruct  Deallocate memory for a patch.
      use cohorts, only : cohort_destruct
      implicit none
      type(patch), pointer :: pp
      !---
      type(cohort),pointer :: cop, cop_tmp

      ! destroy other cohorts
      cop => pp%tallest
      do while( associated(cop) )
        cop_tmp => cop
        cop => cop%shorter
        call cohort_destruct(cop_tmp)
      enddo
      !nullify(pp%tallest)
      !nullify(pp%shortest)

      ! deallocate memory
!      deallocate( pp%Soilmoist )
      deallocate( pp )
      nullify( pp )

      end subroutine patch_destruct


      !*********************************************************************
      subroutine patch_copy(ppin,ppout)
!@sum patch_copy  Copy values in ppin into ppout, only patch-level variables.
!@+   Does not copy cohorts or entcell pointer.
!@+   Patch ppout should already be allocated.
      type(patch),intent(in) :: ppin
      type(patch),intent(out) :: ppout

!      if (associated(ppin)) then
!         if (associated(ppout)) then
         ppout%age = ppin %age
         ppout%area = ppin%area
         ppout%soil_type = ppin%soil_type
            ppout%Tpool(:,:,:) = ppin%Tpool(:,:,:)
         !* Other: nm, albedo, fluxes if desired.
!         else
!            call stop_model("patch_copy: ppout not allocated",255)
!         endif
         !If blank ppin, then do nothing.  May want to stop_model.
 !     endif
      end subroutine patch_copy

      !*********************************************************************
      subroutine patch_print(iu,pp,prefix)
!@sum patch_print  Print contents of a patch.
      use cohorts, only : cohort_print
      integer, intent(in) :: iu
      type(patch), intent(in) :: pp
      character*(*), optional, intent(in) :: prefix !Optional text
      !---
      integer n, nc, m, i
      type(cohort),pointer :: cop
      character*8 prefix_c

      prefix_c = "        "

      write(iu,1) prefix,"area",pp%area
      write(iu,1) prefix,"age ",pp%age
      write(iu,1) prefix,"fracroot",pp%fracroot
      write(iu,1) prefix,"GCANOPY ",
     &     pp%GCANOPY
      write(iu,1) prefix,"Ci ",pp%Ci
      write(iu,1) prefix,"nm ",pp%nm
      write(iu,1) prefix,"soil H2O ",pp%Soilmoist
      write(iu,'(a,"albedo:")') prefix
      do n=1,N_BANDS
        write(iu,'(a,"      ",99e23.16)') prefix,pp%albedo(n)
      enddo
      write(iu,'(a,"Tpool:")') prefix
      do m=1,PTRACE
        do n=1,NPOOLS,3
         do i=1,N_CASA_LAYERS
          write(iu,'(a,"      ",i1,"  ",e23.16,e23.16,e23.16)') prefix,m
     &         ,pp%Tpool(m,n,i),pp%Tpool(m,n+1,i),pp%Tpool(m,n+2,i)
         end do
        enddo
      enddo
      write(iu,'(a,a," = ",i7)') prefix,"soil type",pp%soil_type


      write(iu,1) prefix,"nm	 ",pp%nm		   
      write(iu,1) prefix,"Ntot	",pp%Ntot		   
      write(iu,1) prefix,"LMA   ",pp%LMA		   
      write(iu,1) prefix,"LAI   ",pp%  LAI             
      write(iu,1) prefix,"h       ",pp%  h              !
      write(iu,1) prefix,"crown_dx  ",pp%  crown_dx       !
      write(iu,1) prefix,"crown_dy  ",pp%  crown_dy       !
      write(iu,1) prefix,"clump      ",pp%  clump          !
      write(iu,1) prefix,"C_fol       ",pp%  C_fol          !
      write(iu,1) prefix,"N_fol        ",pp%  N_fol          !
      write(iu,1) prefix,"C_w         ",pp%  C_w           ! 
      write(iu,1) prefix,"N_w       ",pp%  N_w           ! 
      write(iu,1) prefix,"C_lab      ",pp%  C_lab          !
      write(iu,1) prefix,"N_lab       ",pp%  N_lab          !
      write(iu,1) prefix,"C_froot      ",pp% C_froot        !
      write(iu,1) prefix,"N_froot   ",pp%  N_froot        !
      write(iu,1) prefix,"C_root    ",pp%  C_root        ! 
      write(iu,1) prefix,"N_root    ",pp%  N_root        ! 
      write(iu,1) prefix,"Ci        ",pp%  Ci             !
      write(iu,1) prefix,"GCANOPY   ",pp%  GCANOPY        !
      write(iu,1) prefix,"GPP        ",pp%  GPP            !
      write(iu,1) prefix,"IPP        ",pp%  IPP            !
      write(iu,1) prefix,"NPP        ",pp%  NPP            !
      write(iu,1) prefix,"R_auto     ",pp%  R_auto         !
      write(iu,1) prefix,"R_root     ",pp%  R_root         !
      write(iu,1) prefix,"N_up        ",pp%  N_up           !
      write(iu,1) prefix,"z0          ",pp%  z0              
      write(iu,1) prefix,"albedo(N" ,pp%  albedo(:) 
      write(iu,1) prefix,"betad     ",pp%  betad           
      write(iu,1) prefix,"TRANS_SW  ",pp%  TRANS_SW        
      write(iu,1) prefix,"CO2flux    ",pp%  CO2flux         
      write(iu,1) prefix,"Soil_resp",pp%Soil_resp       
      write(iu,1) prefix,"Soilmoist",pp%Soilmoist(:)
      write(iu,1) prefix,"fuel ",pp%fuel		   
      write(iu,1) prefix,"ign_rate",pp%ignition_rate   
      write(iu,1) prefix,"lambda1(",pp%  lambda1(:) !
      write(iu,1) prefix,"d_r",pp%  disturbance_rate
      write(iu,1) prefix,"C_total	",pp%C_total	   
      write(iu,1) prefix,"C_growth  ",pp%  C_growth        

      if ( associated(pp%LAIpft) )
     & write(iu,1) prefix,"LAIpft",pp%LAIpft
      if ( associated(pp%betadl) )
     & write(iu,1) prefix,"betadl",pp%betadl
      if ( associated(pp%Reproduction) )
     & write(iu,1) prefix,"Reproduction",
     &     pp%Reproduction

      write(iu,'(a,"cohorts:")') prefix
      cop => pp%tallest
      nc = 0
      do while( associated(cop) )
        nc = nc + 1
        write( prefix_c, '(i2,"      ")' ) nc
        call cohort_print(iu, cop, prefix//prefix_c)
        cop => cop%shorter
      enddo

 1    format(a,a," = ",99e23.16)  ! e12.5
      end subroutine patch_print

!*************************************************************************

      subroutine patch_print_diag(iu,pp,prefix)
!@sum Print other calculated diagnostics not output by patch_print
!@auth N.Y.Kiang
      use cohorts, only : cohort_print, cohort_print_diag
      integer, intent(in) :: iu
      type(patch), intent(in) :: pp
      character*(*), optional, intent(in) :: prefix !Optional text
      !---
      integer n, nc, m, i
      type(cohort),pointer :: cop
      character*8 prefix_c

      prefix_c = "        "

      write(iu,'(a,"albedo:")') prefix
      do n=1,N_BANDS
        write(iu,'(a,"      ",99e23.16)') prefix,pp%albedo(n)
      enddo
      write(iu,'(a,"Tpool:")') prefix
      do m=1,PTRACE
        do n=1,NPOOLS,3
         do i=1,N_CASA_LAYERS
          write(iu,'(a,"      ",i1,"  ",e23.16,e23.16,e23.16)') prefix,m
     &         ,pp%Tpool(m,n,i),pp%Tpool(m,n+1,i),pp%Tpool(m,n+2,i)
         end do
        enddo
      enddo

      write(iu,1) prefix,"LMA   ",pp%LMA		   
      write(iu,1) prefix,"LAI   ",pp%  LAI             
      write(iu,1) prefix,"h       ",pp%  h              !
      write(iu,1) prefix,"crown_dx  ",pp%  crown_dx       !
      write(iu,1) prefix,"crown_dy  ",pp%  crown_dy       !
      write(iu,1) prefix,"clump      ",pp%  clump          !
      write(iu,1) prefix,"C_fol       ",pp%  C_fol          !
      write(iu,1) prefix,"C_w         ",pp%  C_w           ! 
      write(iu,1) prefix,"C_lab      ",pp%  C_lab          !
      write(iu,1) prefix,"C_froot      ",pp% C_froot        !
      write(iu,1) prefix,"C_root    ",pp%  C_root        ! 
      write(iu,1) prefix,"albedo(N" ,pp%  albedo(:) 
      write(iu,1) prefix,"C_total	",pp%C_total	   

      if ( associated(pp%LAIpft) )
     & write(iu,1) prefix,"LAIpft",pp%LAIpft
      if ( associated(pp%Reproduction) )
     & write(iu,1) prefix,"Reproduction",
     &     pp%Reproduction

      write(iu,'(a,"cohorts:")') prefix
      cop => pp%tallest
      nc = 0
      do while( associated(cop) )
        nc = nc + 1
        write( prefix_c, '(i2,"      ")' ) nc
        call cohort_print_diag(iu, cop, prefix//prefix_c)
        cop => cop%shorter
      enddo

 1    format(a,a," = ",99e23.16)  ! e12.5
      end subroutine patch_print_diag
!*************************************************************************

      subroutine patch_extract_pfts(pp, vfraction)
!@sum patch_extract_pfts  Extract cohort density fraction from a patch.
!@+   This is a place holder routine, as the calculation is meaningless.
      use ent_pfts
      use ent_const
      use cohorts, only : cohort_print
      type(patch) :: pp
      real*8 :: vfraction(:)
      !---
      type(cohort),pointer :: cop
      real*8 :: density

      vfraction(:) = 0.d0
      density = 0.d0
      cop => pp%tallest
      do while( associated(cop) )
        vfraction( cop%pft + COVEROFFSET ) = 
     &       vfraction( cop%pft + COVEROFFSET ) + cop%n
        density = density + cop%n
        cop => cop%shorter
      enddo

      if ( density>0.d0 ) then
        vfraction(:) = vfraction(:)/density
      else
        if ( pp%soil_type == 1 ) then
          vfraction(COVER_SAND) = 1.d0
        else
          vfraction(COVER_DIRT) = 1.d0
        endif
      endif

      end subroutine patch_extract_pfts
     

!**************************************************************************

      subroutine print_Tpool(Tpool)
!@sum Print carbon pools in Tpool.
      real*8 :: Tpool(PTRACE,NPOOLS,N_CASA_LAYERS)
      print*,'Tpool(CARBON,LEAF,:)',Tpool(CARBON,LEAF,:)
      print*,'Tpool(CARBON,FROOT,:)',Tpool(CARBON,FROOT,:)
      print*,'Tpool(CARBON,WOOD,:)',Tpool(CARBON,WOOD,:)
      print*,'Tpool(CARBON,SURFMET,:)',Tpool(CARBON,SURFMET,:)
      print*,'Tpool(CARBON,SURFSTR,:)',Tpool(CARBON,SURFSTR,:)
      print*,'Tpool(CARBON,SOILMET,:)',Tpool(CARBON,SOILMET,:)
      print*,'Tpool(CARBON,SOILSTR,:)',Tpool(CARBON,SOILSTR,:)
      print*,'Tpool(CARBON,CWD,:)',Tpool(CARBON,CWD,:)
      print*,'Tpool(CARBON,SOILMIC,:)',Tpool(CARBON,SURFMIC,:)
      print*,'Tpool(CARBON,SURFMIC,:)',Tpool(CARBON,SOILMIC,:)
      print*,'Tpool(CARBON,SLOW,:)',Tpool(CARBON,SLOW,:)
      print*,'Tpool(CARBON,PASSIVE,:)',Tpool(CARBON,PASSIVE,:)
      end subroutine print_Tpool

!**************************************************************************


      subroutine patch_split(pp, area, pp_new)
!@sum patch_split.  Split a 0- or 1-cohort patch into two cohorts.

      use cohorts, only : insert_cohort
      type(patch), pointer :: pp, pp_new
      real*8 area
      !---
      type(entcelltype), pointer :: gp

      gp => pp%cellptr
      call insert_patch(gp, area, pp%soil_type)
      pp_new => gp%youngest

      if ( associated(pp%tallest) ) then ! have cohort: will insert similar
        call insert_cohort(pp_new, pp%tallest%pft, pp%tallest%n)
      endif

      pp%area = pp%area - area
      if ( pp%area < -1.d-6 )
     &     call stop_model("patch_split: area > pp%area",255)
      pp%area = max( 0.d0, pp%area )

      call patch_merge_data( pp_new, 0.d0, pp, 1.d0 )

      end subroutine patch_split


      subroutine patch_set_pft(pp, pft)
!@sum patch_set_pft  For mosaicked veg, insert a cohort of pft into a patch.
!!!   this subroutine works only with 0 or 1 coherts per patch !
      use cohorts, only : insert_cohort
      type(patch), pointer :: pp
      integer :: pft
      !---

      if ( .not. associated(pp%tallest) ) then
        ! no cohorts, will insert one
        call insert_cohort(pp, pft, 0.d0)
      else if ( associated(pp%tallest, pp%shortest) ) then
        ! one cohort - just change pft
        pp%tallest%pft = pft
      else
        ! more then one cohort
        call stop_model("patch_set_pft: called for >1 cohorts",255)
      endif
       
      
      end subroutine patch_set_pft


      subroutine patch_merge(pp1, pp2)
!@sum patch_merge  Merge two patches into one.
!@+   Merges areas and calls another routine to merge patch data.
      implicit none
      type(patch), pointer :: pp1, pp2
      !---
      type(entcelltype), pointer :: gp
      real*8 :: w1, w2, tot_area

      gp => pp1%cellptr

      tot_area = pp1%area + pp2%area
      w1 = pp1%area/tot_area
      w2 = pp2%area/tot_area
      call patch_merge_data( pp1, w1, pp2, w2 )

      pp1%area = tot_area
      pp1%area = min( 1.d0, pp1%area )
      call delete_patch( gp, pp2 )

      end subroutine patch_merge

      subroutine patch_merge_data( pp1, w1, pp2, w2 )
!@sum patch_merge_data.  Merge the data in two patches into one patch.
      use cohorts, only : cohort_merge_data
      implicit none
      type(patch) :: pp1, pp2
      real*8 :: w1, w2
      !----Local----

      ! do we need to reset age? - Ideally, the area that is converted would become
      !      a new patch entirely, rather than being merged with another patch.
      !      In this case, it's okay to keep the original patch ages.
      !      If the crop patch is completely new, then age=0.
      !pp%age
      ! area is set outside
      !pp%area

      !nullify(pp%crad%heights)
      !nullify(pp%crad%lai)
      !pp%crad%gortclump = 0.d0

      !pp%soil_type = -1         ! set to undefined soil type (maybe use -1?)

      !* Weighted average for Reproduction okay.
      pp1%Reproduction(:) =w1*pp1%Reproduction    +w2*pp2%Reproduction

      !* Cohort variables below should be updated with summarize_patch AFTER
      !* calling cohort_merge_data.
      pp1%nm              =w1*pp1%nm              +w2*pp2%nm              
      pp1%Ntot            =w1*pp1%Ntot            +w2*pp2%Ntot            
      pp1%LAI             =w1*pp1%LAI             +w2*pp2%LAI             
      pp1%LMA             =w1*pp1%LMA             +w2*pp2%LMA             
      pp1%h               =w1*pp1%h               +w2*pp2%h               
      pp1%crown_dx        =w1*pp1%crown_dx        +w2*pp2%crown_dx        
      pp1%crown_dy        =w1*pp1%crown_dy        +w2*pp2%crown_dy        
      pp1%clump           =w1*pp1%clump           +w2*pp2%clump           
      pp1%fracroot        =w1*pp1%fracroot        +w2*pp2%fracroot        
      pp1%C_fol           =w1*pp1%C_fol           +w2*pp2%C_fol           
      pp1%N_fol           =w1*pp1%N_fol           +w2*pp2%N_fol           
      pp1%C_w             =w1*pp1%C_w             +w2*pp2%C_w             
      pp1%N_w             =w1*pp1%N_w             +w2*pp2%N_w             
      pp1%C_lab           =w1*pp1%C_lab           +w2*pp2%C_lab           
      pp1%N_lab           =w1*pp1%N_lab           +w2*pp2%N_lab           
      pp1%C_froot         =w1*pp1%C_froot         +w2*pp2%C_froot         
      pp1%N_froot         =w1*pp1%N_froot         +w2*pp2%N_froot         
      pp1%C_root          =w1*pp1%C_root          +w2*pp2%C_root          
      pp1%N_root          =w1*pp1%N_root          +w2*pp2%N_root          
      !* Value of Ci doesn't matter -- saved only for diagnostic testing.
      pp1%Ci              =w1*pp1%Ci              +w2*pp2%Ci              
 
      !* Okay to do weighted average for fluxes.
      pp1%GCANOPY         =w1*pp1%GCANOPY         +w2*pp2%GCANOPY         
      pp1%GPP             =w1*pp1%GPP             +w2*pp2%GPP             
      pp1%IPP             =w1*pp1%IPP             +w2*pp2%IPP             
      pp1%NPP             =w1*pp1%NPP             +w2*pp2%NPP             
      pp1%R_auto          =w1*pp1%R_auto          +w2*pp2%R_auto          
      pp1%R_root          =w1*pp1%R_root          +w2*pp2%R_root          
      pp1%N_up            =w1*pp1%N_up            +w2*pp2%N_up            
      
      !* Best to summarize_patch after calling cohort_merge_data, but this has
      !* to be done first by biophysics, anyway, so no harm.
      pp1%betad           =w1*pp1%betad           +w2*pp2%betad           
      pp1%betadl(:)       =w1*pp1%betadl(:)       +w2*pp2%betadl(:)       
                                                      
      !* Weighted average okay.
      pp1%C_total         =w1*pp1%C_total         +w2*pp2%C_total         
      pp1%C_growth        =w1*pp1%C_growth        +w2*pp2%C_growth        
                                                        
      !* Need to recalculate radiative transfer vars after calling cohort_merge_data.
      pp1%z0              =w1*pp1%z0              +w2*pp2%z0
      pp1%albedo          =w1*pp1%albedo          +w2*pp2%albedo
      pp1%TRANS_SW        =w1*pp1%TRANS_SW        +w2*pp2%TRANS_SW

      !* Weighted average okay.
      pp1%CO2flux         =w1*pp1%CO2flux         +w2*pp2%CO2flux         
      pp1%Soil_resp       =w1*pp1%Soil_resp       +w2*pp2%Soil_resp       
                                                     
      pp1%Soilmoist(:)    =w1*pp1%Soilmoist(:)    +w2*pp2%Soilmoist(:)    
                                                         
      pp1%fuel            =w1*pp1%fuel            +w2*pp2%fuel            
      pp1%ignition_rate   =w1*pp1%ignition_rate   +w2*pp2%ignition_rate   
      pp1%lambda1         =w1*pp1%lambda1         +w2*pp2%lambda1         
      pp1%disturbance_rate=w1*pp1%disturbance_rate
     &     +w2*pp2%disturbance_rate
                                                            
      !* Need to summarize_patch after calling cohort_merge_data.
      pp1%LAIpft(:)       =w1*pp1%LAIpft(:)       +w2*pp2%LAIpft(:)       
      pp1%Tpool(:,:,:)    =w1*pp1%Tpool(:,:,:)    +w2*pp2%Tpool(:,:,:)

      ! process cohorts
      if ( (.not. associated(pp1%tallest)) .and.
     &     (.not. associated(pp2%tallest)) ) return ! no cohorts

      if ( (.not. associated(pp1%tallest)) .or.
     &     (.not. associated(pp2%tallest)) ) call stop_model(
     &     "Can't merge patches with and without cohorts",255)

      if ( .not. associated(pp1%tallest,pp1%shortest) ) call stop_model(
     &     "Can't merge patches with more than 1 cohort",255)

      if ( .not. associated(pp2%tallest,pp2%shortest) ) call stop_model(
     &     "Can't merge patches with more than 1 cohort",255)

      !* CALL THIS AT BEGINNING OF patch_merge_data ROUTINE then summarize_patch *!
      !* This should loop through cohorts to merge all cohorts in both patches
      call cohort_merge_data( pp1%tallest, w1, pp2%tallest, w2 )

      end subroutine patch_merge_data


      function patch_has_pft( pp, pft )
!@sum Logical if a pft exists in a patch.
      logical patch_has_pft
      type(patch), intent(in) :: pp
      integer, intent(in) :: pft
      !---
      type(cohort),pointer :: cop
      
      patch_has_pft = .false.
      cop => pp%tallest
      do while( associated(cop) )
        if ( cop%pft == pft ) then
          patch_has_pft = .true.
          return
        endif
        cop => cop%shorter
      enddo

      end function patch_has_pft


      subroutine patch_delete_cohort( pp, cop_del )
!@sum patch_delete_cohort  Delete a cohort from a patch
      use cohorts, only : cohort_destruct
      type(patch), pointer :: pp
      type(cohort), pointer :: cop_del
      !---
      type(cohort),pointer :: cop, cop_next
      
      cop => pp%tallest
      do while( associated(cop) )

        cop_next => cop%shorter
        if ( associated( cop, cop_del ) ) then
          if ( associated( cop%taller ) ) then
            cop%taller%shorter => cop%shorter
          else
            pp%tallest => cop%shorter
          endif
          if ( associated( cop%shorter ) ) then
            cop%shorter%taller => cop%taller
          else
            pp%shortest => cop
          endif
          call cohort_destruct( cop )
        endif

        cop => cop_next
      enddo


      end subroutine patch_delete_cohort

!----------------------------------------------------------------------
      real*8 function patch_carbon(pp,
     &     pp_Cfol, pp_Cstem, pp_Croot, pp_Csoil)
!@sum patch_carbon (kg-C m-2) - Return biomass + soil carbon per land area
!@+   Also optionally return component carbon pools (kg-C m-2).      
      use cohorts, only : cohort_carbon
      type(patch), pointer :: pp
      real*8,optional :: pp_Cfol
      real*8,optional :: pp_Cstem
      real*8,optional :: pp_Croot
      real*8,optional :: pp_Csoil
      !--Local-----------------
      type(cohort), pointer :: cop
      real*8 :: cohort_kgC, cohort_kgCm2

      integer :: i,k
      real*8 :: soil_gCm2

      cohort_kgC = 0.d0
      cohort_kgCm2 = 0.d0
      soil_gCm2 = 0.d0

      if (present(pp_Cfol)) pp_Cfol = 0.d0
      if (present(pp_Cstem)) pp_Cstem = 0.d0
      if (present(pp_Croot)) pp_Croot = 0.d0

      !* Biomass.
      cop => pp%tallest
      do while (associated(cop))
         cohort_kgCm2 = cohort_kgCm2 + cohort_carbon(cop)*cop%n

         if (present(pp_Cfol)) pp_Cfol = pp_Cfol + cop%C_fol
         if (present(pp_Cstem)) pp_Cstem = pp_Cstem + 
     &        cop%C_sw + cop%C_hw
         if (present(pp_Croot)) pp_Croot = pp_Croot + 
     &        cop%C_froot + cop%C_croot

         cop => cop%shorter
      end do

      !* Soil.
      do i=1,N_CASA_LAYERS
         do k=(NLIVE+1),NPOOLS
            soil_gCm2 = soil_gCm2 + pp%Tpool(CARBON,k,i)
         enddo
      enddo
      if (present(pp_Csoil)) pp_Csoil = soil_gCm2*1.d-3

      !* Total patch carbon per m2-ground.
      patch_carbon = cohort_kgCm2 + soil_gCm2*1.d-3

      end function patch_carbon
!----------------------------------------------------------------------

      end module patches
