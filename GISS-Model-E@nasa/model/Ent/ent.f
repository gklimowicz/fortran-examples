      module ent
      
!@sum Contains main routines to perform single time step Ent model simulation
!@sum on a single grid cell, entcell, whose parameters are set previously by
!@sum ent_driver routines and the subroutine ent_model called by a 
!@sum main program.

!@auth N.Y. Kiang

!#define ENT_STANDALONE_DIAG TRUE
!#define DEBUG

      use ent_const
      use ent_types

      implicit none
      private
      save

      public ent_integrate 
      public update_veg_structure
      public ent_biophysics
      public ent_ecosystem_dynamics
      public spinup_soil_bgc_step

      contains
      !*********************************************************************

      subroutine spinup_soil_bgc_step(dtsec, ecp, buf)
      use soilbgc, only : soil_bgc
      use entcells, only : summarize_entcell
      real*8 :: dtsec  !dt in seconds
      type(entcelltype) :: ecp
      real*8, intent(in) :: buf(:,:)
      !---
      type(patch),pointer :: pp
      integer :: n, i

      pp => ecp%oldest
      n = 0
      do while (ASSOCIATED(pp)) 
        n = n + 1

cddd        buf(1,n) = pp%acc_litter_SURFMET
cddd        buf(2,n) = pp%acc_litter_SOILMET
cddd        buf(3,n) = pp%acc_litter_SURFSTR
cddd        buf(4,n) = pp%acc_litter_SOILSTR
cddd        buf(5,n) = pp%acc_litter_CWD

        i = 1
        pp%Tpool(CARBON,SURFMET,i) = pp%Tpool(CARBON,SURFMET,i) 
     &     + buf(1,n)
        pp%Tpool(CARBON,SOILMET,i) = pp%Tpool(CARBON,SOILMET,i) 
     &     + buf(2,n)
        pp%Tpool(CARBON,SURFSTR,i) = pp%Tpool(CARBON,SURFSTR,i)
     &     + buf(3,n)
        pp%Tpool(CARBON,SOILSTR,i) = pp%Tpool(CARBON,SOILSTR,i) 
     &     + buf(4,n)
        pp%Tpool(CARBON,CWD,i) = pp%Tpool(CARBON,CWD,i) 
     &     + buf(5,n)

cddd        if ( associated( pp%tallest ) ) then
cddd          write(638,*) pp%tallest%pft, buf(1:5,n)
cddd        endif
        call soil_bgc( dtsec, pp, buf(6,n)/dtsec, buf(7,n)/dtsec )

        pp => pp%younger 
      end do

      !call summarize_entcell(ecp)

      end subroutine spinup_soil_bgc_step


      subroutine ent_integrate(dtsec, ecp, update_day, config)
!@sum Main routine to control Ent biophysics/biogeochemistry. 
!@+   (Patch ecological dynamics TBA)
      use cohorts
      use patches
      use biophysics, only : photosynth_cond
      use soilbgc, only : soil_bgc
      use phenology, only : clim_stats, pheno_update, veg_update
!      use canopyrad, only : recalc_radpar
      use entcells, only : summarize_entcell, entcell_print

      implicit none
      real*8 :: dtsec  !dt in seconds
      !type(timestruct),pointer :: tt !Time in year.fraction, Greenwich Mean Time
      type(entcelltype) :: ecp
      logical :: update_day
      type(ent_config) :: config 

      !-----local--------
      integer :: patchnum
      type(patch),pointer :: pp


      call clim_stats(dtsec,ecp,config,update_day)

      !* Patch dynamics
      !if (config%do_patchdynamics) then
      !  call patch_dynamics(pp,monthlyupdate)
      ! call summarize_entcell(ecp)
      !endif

      !* Dynamic phenology
      if (update_day) then 
        call update_veg_structure(ecp, config)
      endif

      !* Biophysics
      patchnum = 0
      pp => ecp%oldest 
      do while (ASSOCIATED(pp)) 
        patchnum = patchnum + 1
        !call patch_print(771,pp," ff ")
        call photosynth_cond(dtsec, pp)

        if (config%do_soilresp) then
          call soil_bgc(dtsec, pp)
        endif
        pp%CO2flux = -pp%NPP + pp%Soil_resp
        
        pp%age = pp%age + dtsec

          !*********** DIAGNOSTICS FOR PLOTTING ********************!
#ifdef  ENT_STANDALONE_DIAG         
        call summarize_patch(pp)
        call ent_diagnostics(patchnum, pp)
#endif
          !*********************************************************!

        pp => pp%younger 
      end do

      call summarize_entcell(ecp)

      call debug_diags(ecp)

#ifdef DEBUG
      print *,"End of ent_integrate"
      call entcell_print(6, ecp)
#endif

#ifdef ENT_STANDALONE_DIAG
      call ent_diagnostics_entcell(ecp)
!      call entcell_print(6, ecp)
#endif

      end subroutine ent_integrate



      subroutine update_veg_structure(ecp, config)
!@sum Update prognostic vegetation structure (seasonal) at the end of day
      use phenology, only : pheno_update, veg_update
      use entcells, only : summarize_entcell, entcell_print
     &     ,entcell_carbon
      implicit none
      type(entcelltype) :: ecp
      type(ent_config) :: config 
      !-----local--------
      type(patch),pointer :: pp
      real*8 c_before, c_after


      c_before = entcell_carbon(ecp)
      !* Loop through patches
      pp => ecp%oldest 
      do while (ASSOCIATED(pp)) 
        
        if (config%do_phenology_activegrowth) then
          call pheno_update(pp)
          !call patch_print(771,pp," bb ")
          call veg_update(pp,config)
        endif

        pp => pp%younger 
      end do

      call summarize_entcell(ecp)

      ecp%daylength(1) = ecp%daylength(2)
      ecp%daylength(2) = 0.d0

      c_after = entcell_carbon(ecp)

#ifdef DEBUG
      if ( abs(c_after-c_before) > 1.d-10 ) then
        write(904,*) "dC_cell ", c_after-c_before, c_before
      endif
#endif

      end subroutine update_veg_structure

      !*********************************************************************











      !*********************************************************************
      subroutine ent_biophysics(dtsec, ecp, config)
!@sum  Photosynthesis CO2 uptake and conductance of water vapor.
!@+    This optional routine may be used by a GCM or land surface driver
!@+    to calculate only surface conductance without biogeochemistry.
!@+    Useful for implicit schemes, but not used in the default setup.
!@+    If do_soilresp, then also  soil respiration for net CO2 fluxes.
      use biophysics, only : photosynth_cond
      use soilbgc, only : soil_bgc
      use patches, only : summarize_patch
      use entcells, only : summarize_entcell, entcell_print
      implicit none
      real*8 :: dtsec  !dt in seconds
      type(entcelltype) :: ecp
      type(ent_config) :: config
      !---Local--------
      type(patch),pointer :: pp
      integer :: patchnum

      call stop_model("ent_biophysics: obsolete code",255)
       
cddd      patchnum = 0
cddd      pp => ecp%oldest
cddd      do while(ASSOCIATED(pp))
cddd        patchnum = patchnum + 1
cddd        !print*,'NEXT PATCH'
cddd        !print*,'Calling photosynth_cond'
cddd        call photosynth_cond(dtsec, pp)
cddd        if (config%do_soilresp) then 
cddd          !print*,'Calling soil_bgc'
cddd          call soil_bgc(dtsec,pp)
cddd          pp%CO2flux = -pp%NPP + pp%Soil_resp
cddd        ! Litter is updated daily in ent_prescribe_vegupdates.
cddd        ! is do_soilresp flag ok or different flag is needed ?
cddd        ! if ( dailyupdate ) call litter(pp) 
cddd
cddd          !*********** DIAGNOSTICS FOR PLOTTING ********************!
cddd#ifdef ENT_STANDALONE_DIAG         
cddd          call summarize_patch(pp)
cddd          call ent_diagnostics(patchnum,pp)
cddd#endif
cddd          !*********************************************************!
cddd        else 
cddd          pp%CO2flux = UNDEF
cddd        endif
cddd        pp%age = pp%age + dtsec
cddd        !call summarize_patch(pp)
cddd
cddd        pp => pp%younger
cddd
cddd      end do
cddd      call summarize_entcell(ecp)
cddd
cddd#ifdef DEBUG
cddd      print *,"End of ent_biophysics"
cddd      call entcell_print(6, ecp)
cddd      print *,"*"
cddd#endif
      end subroutine ent_biophysics
      !*********************************************************************

      subroutine ent_ecosystem_dynamics(dtsec,tt,ecp,ALBEDO_FLAG)
!@sum Ent ecosystem dynamics. UNDER DEVELOPMENT
!@auth N.Y.Kiang
      use phenology
!      use canopyrad
      use disturbance
      use cohorts, only : reorganize_cohorts
      use patches, only : reorganize_patches, summarize_patch
      use entcells, only : summarize_entcell

      real*8,intent(in) :: dtsec
      type(timestruct),pointer :: tt
      type(entcelltype) :: ecp
      logical, optional, intent(in) :: ALBEDO_FLAG
      !---
      type(patch), pointer :: pp

      call stop_model("ent_ecosystem_dynamics: not supportd in gcm",255)

      pp => ecp%oldest
      do while (ASSOCIATED(pp)) 
      !#### THIS LOOP: NEED TO REPLACE ALL CALLS WITH ecp TO CALLS WITH pp ##

        if (ALBEDO_FLAG) then
!          call get_patchalbedo(pp)  !Place holder
        end if

!        call ent_integrate(dtsec,ecp,0.0) !Biophysics, respiration

        if (STRUCT_FLAG(tt,ecp)) then
          call reorganize_cohorts(pp)
!          call phenology_update (dtsec,tt, pp) !UPDATE LAI
!          call recalc_radpar (pp) !UPDATE canopy radiative transfer
        end if
        call summarize_patch(pp)
      
        !Flag when it's time to update disturbance.
        !May be at set time intervals, or function of biomass accumulation, etc.
        !For now, monthly update as a place holder.
        if (STRUCT_FLAG_MONTH(tt,ecp)) then
        !* Update phenology and disturbance
        !call phenology_update (dtsec,tt, pp) !UPDATE LAI - put in ent_integrate
          call fire_frequency_cell (dtsec,tt, ecp) !DUMMY
!          call recalc_radpar_cell (ecp) !
          call reorganize_patches(ecp)
          call calc_cell_disturbance_rates(dtsec,tt,ecp)
        else
          call calc_cell_disturbance_rates(dtsec,tt,ecp)
        end if
        pp => pp%younger
      end do

      call summarize_entcell(ecp)

        end subroutine ent_ecosystem_dynamics

      !*********************************************************************


      function STRUCT_FLAG(tt, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
        type(timestruct) :: tt
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        !------local------

        call stop_model("STRUCT_FLAG: not allowed in gcm",255)

        update_struct = STRUCT_FLAG_DAY(tt,ecp)

      end function STRUCT_FLAG
      !*********************************************************************

      function STRUCT_FLAG_DAY(tt, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
!@sum Below is a simple end-of-day flag, but can make more
!@sum sophisticated as a function of biomass increment, etc.

        type(timestruct) :: tt
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        !-----local----------
        real*8 :: hourfrac

        call stop_model("STRUCT_FLAG_DAY: not allowed in gcm",255)
        
        hourfrac = tt%hour + tt%minute/60.0 + tt%seconds/3600.0
!        if (hourfrac.le.dtsec) then !Midnight
        if (hourfrac.eq.0.0) then  !Midnight
          update_struct = .true.
        else
           update_struct = .false.
        end if
      end function STRUCT_FLAG_DAY

      !*********************************************************************

      function STRUCT_FLAG_MONTH(tt, ecp) Result(update_struct)
!@sum Flag to determine if it's time to update vegetation structure.
!@sum Below is a simple beginning-of-the-month flag, but can make more
!@sum sophisticated as a function of biomass increment, etc.

        type(timestruct) :: tt
        type(entcelltype) :: ecp !Not needed this version, but will be.
        logical :: update_struct
        
        real*8 :: hourfrac

        call stop_model("STRUCT_FLAG_MONTH: not allowed in gcm",255)
        
        hourfrac = tt%hour + tt%minute/60.0 + tt%seconds/3600.0
        if ((tt%day.eq.1).and.
     &       (hourfrac.eq.0.0)) then
          update_struct = .true.
        else
           update_struct = .false.
        end if
      end function STRUCT_FLAG_MONTH


!*****************************************************************************
      subroutine ent_diagnostics(patchnum, pp)
      !*********** DIAGNOSTICS FOR PLOTTING ********************!
      !use patches, only : patch_print
      implicit none
      integer :: patchnum
      type(patch),pointer :: pp
      !---Local------
      integer :: tmp_pft
      real*8 :: tmp_n,tmp_senescefrac,tmp_Sacclim
      type(cohort),pointer :: cop
!#ifdef DEBUG

      tmp_pft = -1
      tmp_n = -1
      tmp_senescefrac = 0.d0
      tmp_Sacclim = -1

      !call patch_print(6,pp)
      if ( ASSOCIATED(pp%tallest) ) then
        tmp_pft = pp%tallest%pft
        tmp_n = pp%tallest%n
        tmp_senescefrac = pp%tallest%senescefrac
        tmp_Sacclim = pp%tallest%Sacclim
        cop => pp%tallest
        do while (ASSOCIATED(cop))
          cop => cop%shorter
        end do
      endif

      write(995,'(i5,3(1pe16.8),i5,100(1pe16.8))') !Fluxes are positive up.
     &     patchnum,pp%cellptr%IPARdir,pp%cellptr%IPARdif, 
     &     pp%cellptr%coszen,
     &     tmp_pft,tmp_n,pp%lai, pp%h, pp%Tpool(CARBON,:,:), 
     &     pp%C_fol, pp%C_w, pp%C_froot, pp%C_root, pp%C_lab,
     &     pp%Reproduction(tmp_pft),
     &     pp%TRANS_SW,
     &     pp%Ci, pp%GPP,pp%R_auto,pp%Soil_resp,
     &     pp%NPP,pp%CO2flux,pp%GCANOPY,pp%IPP,
     &     tmp_senescefrac,tmp_Sacclim,pp%c_total,
!### HACK: c_growth is Igor's hack to store daily growth respiration somewhere
!### HACK: N_up is temporarily litterfall, using unused variable -NK
     &     pp%c_growth,pp%N_up,pp%betad

      if (pp%GPP.lt.0.d0) then
        print *,"ent.f: BAD GPP:",pp%lai, pp%GPP
      endif
!      write(999,*) pp%cellptr%Soilmp, pp%Soilmoist
!     &     ,pp%cellptr%betad, pp%cellptr%betadl
!      write(994,*) pp%cellptr%GCANOPY
!#endif
      end subroutine ent_diagnostics
          !*********************************************************!

      subroutine ent_diagnostics_entcell(ecp)
      implicit none
      type(entcelltype) :: ecp

      write(996,*) ecp%area,ecp%IPARdir,ecp%IPARdif,ecp%LAI,ecp%fv
     &     ,ecp%Tpool(CARBON,:,:)
     &     ,ecp%C_fol, ecp%C_w, ecp%C_froot, ecp%C_root, ecp%C_lab
     &     ,ecp%TRANS_SW
     &     ,ecp%Ci, ecp%GPP,ecp%R_auto,ecp%Soil_resp
     &     ,ecp%NPP,ecp%CO2flux,ecp%GCANOPY

      end subroutine ent_diagnostics_entcell


      subroutine debug_diags(ecp)
      use ent_debug_mod
      type(entcelltype), intent(in) :: ecp
      !---
      type(patch),pointer :: pp
      type(cohort), pointer :: cop

      real*8 :: defacc_scale= 1.d0/(3600*24*1000.d0)  !!! hack
      real*8 :: area, n, scale
      integer i, k, pft


      ent_d%total(:) = 0.d0
      ent_d%C_lab(:) = 0.d0
      ent_d%C_fol(:) = 0.d0
      ent_d%C_sw(:) = 0.d0
      ent_d%C_hw(:) = 0.d0
      ent_d%C_froot(:) = 0.d0
      ent_d%C_croot(:) = 0.d0
      ent_d%C_soil(:) = 0.d0
      ent_d%phenofactor(:) = 0.d0
      ent_d%betad(:) = 0.d0

      pp => ecp%oldest

      do while (associated(pp))
        area = pp%area

        pft = 0
        cop => pp%tallest
        do while (associated(cop))
          pft = cop%pft
          n = cop%n
          scale = n*1.d-3 * defacc_scale      ! was : n*area*1.d-3

          ent_d%C_lab(pft) = ent_d%C_lab(pft) + cop%C_lab*scale
          ent_d%C_fol(pft) = ent_d%C_fol(pft) + cop%C_fol*scale
          ent_d%C_sw(pft) = ent_d%C_sw(pft) + cop%C_sw*scale
          ent_d%C_hw(pft) = ent_d%C_hw(pft) + cop%C_hw*scale
          ent_d%C_froot(pft) = ent_d%C_froot(pft) + cop%C_froot*scale
          ent_d%C_croot(pft) = ent_d%C_croot(pft) + cop%C_croot*scale
          ent_d%phenofactor(pft) = ent_d%phenofactor(pft)
     &         + cop%phenofactor*defacc_scale
          ent_d%betad(pft) = ent_d%betad(pft)
     &         + cop%stressH2O*defacc_scale
          !ent_d%C_(pft) = ent_d%C_(pft) + cop%C_*scale

          cop => cop%shorter
        end do

                                !!! assume 1 cohort per patch
        if ( pft > 0) then      ! skip cells with no vegetation
          do i=1,N_CASA_LAYERS
            do k=(NLIVE+1),NPOOLS
              ent_d%C_soil(pft)=ent_d%C_soil(pft)
     &             + pp%Tpool(CARBON,k,i)*1.d-3 * defacc_scale ! *area
            enddo
          enddo
          ent_d%Resp_soil(pft)=ent_d%Resp_soil(pft) + pp%Soil_resp
        endif

        pp => pp%younger
      end do

      ent_d%total(:) = ent_d%C_lab(:) + ent_d%C_fol(:)
     &	   + ent_d%C_sw(:) + ent_d%C_hw(:)
     &     + ent_d%C_froot(:) + ent_d%C_croot(:)
     &     + ent_d%C_soil(:)

      end subroutine debug_diags

      end module ent
