      module ent_make_struct
!@sum ent_make_struct  Off-line module for generating an Ent vegetation
!@+   data structure, given ascii input file of entcell-patch-cohort structure.
      use ent_types
      use ent_const
      use cohorts
      use patches
      use entcells

      public ent_struct_readcsv
      
      contains

!************************************************************************
      subroutine skipstar (iu_entstruct)
      !* Skip a comment line in an Ent structure csv file.
      integer :: iu_entstruct
      !---
      character :: check

      read(iu_entstruct,*) check
      !write(*,*) 'ss',check
      do while (check.eq.'*') 
        read(iu_entstruct,*) check
      end do
      backspace(iu_entstruct)
      end subroutine skipstar

!************************************************************************


      subroutine read_entcell_struct ( ecp, iu_entstruct )
      type(entcelltype) :: ecp
      integer :: iu_entstruct
      !---
      !real*8 :: stext1,stext2,stext3,stext4,stext5 !soil_texture
      
      !read(iu_entstruct,*) stext1,stext2,stext3,stext4,stext5 !soil_texture
      !ecp%soil_texture(1) = stext1
      !ecp%soil_texture(2) = stext2
      !ecp%soil_texture(3) = stext3
      !ecp%soil_texture(4) = stext4
      !ecp%soil_texture(5) = stext5
      call skipstar(iu_entstruct)
      read(iu_entstruct,*) ecp%soil_texture 
      end subroutine read_entcell_struct

!************************************************************************
     
      subroutine read_patch_struct (iu_entstruct
     &     , age, area, soil_type, Tpool)
      use ent_prescribed_drv, only : read_soilcarbon_patch
      implicit none
      integer :: iu_entstruct
      real*8, intent(out) :: age, area
      integer, intent(out) :: soil_type      
      real*8, intent(out) :: Tpool(PTRACE,NPOOLS-NLIVE
     &     ,N_CASA_LAYERS)!prescribed soil pools, g/m2

      call skipstar(iu_entstruct)
      read(iu_entstruct,*) age, area, soil_type

      call read_soilcarbon_patch(iu_entstruct,Tpool)
      end subroutine read_patch_struct

!************************************************************************

      subroutine read_patch_struct_old (pp,iu_entstruct )
      type(patch) :: pp
      integer :: iu_entstruct
      !---
      integer :: layer
!      real*8 :: age, area, soil_type
!     real*8 :: May also add soil texture
!      real*8 :: tc1,tc2,tc3,tc4,tc5,tc6,tc7,tc8,tc9,tc10,tc11,tc12 !TpoolC
!      real*8 :: tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10,tn11,tn12 !TpoolN

      call skipstar(iu_entstruct)
      read(iu_entstruct,*) pp%age, pp%area, pp%soil_type
      do layer = 1,N_CASA_LAYERS
        call skipstar(iu_entstruct)
        read(iu_entstruct,*) pp%Tpool(CARBON,:,layer)
        read(iu_entstruct,*) pp%Tpool(NITROGEN,:,layer)
      end do
      end subroutine read_patch_struct_old

!************************************************************************

      subroutine get_patch_struct (iu_entstruct, option, ec)
!@sum get_patch_struct  For mixed canopies, read csv file for initialization
!@+      values for patch, insert patch, write intialization values.      
      implicit none
      integer :: iu_entstruct
      integer :: option !1-overwrite first dummy patch ec%youngest=ec%oldest
                        !2-insert new patch
      type(entcelltype) :: ec
      !---------------
      real*8 :: age, area
      integer :: soil_type 
      real*8 :: Tpool_soil(PTRACE,NPOOLS-NLIVE
     &     ,N_CASA_LAYERS)!prescribed soil pools, g/m2
     
      call read_patch_struct(iu_entstruct
     &     , age, area, soil_type, Tpool_soil)

      if (option.eq.1) then
         ec%youngest%age = age
         ec%youngest%area = area
         ec%youngest%soil_type = soil_type
         ec%youngest%Tpool(:,(NLIVE+1):NPOOLS,:) = Tpool_soil(:,:,:)
      else
         call insert_patch(ec, area, soil_type)
         ec%youngest%Tpool(:,(NLIVE+1):NPOOLS,:) = Tpool_soil(:,:,:)
      endif
      write(*,*) 'patch age area soil_type Tpool'
     &     ,age,area,soil_type,Tpool_soil
      
      end subroutine get_patch_struct

!************************************************************************

      subroutine read_cohort_struct ( iu_entstruct
     o     , pft, LAImax, hm, LAIinit )
      integer :: iu_entstruct
      integer, intent(out) :: pft
      real*8, intent(out) :: LAImax, hm, LAIinit
      !---
      
      call skipstar(iu_entstruct)
      read(iu_entstruct,*) pft
      call skipstar(iu_entstruct)
      read(iu_entstruct,*) LAImax,hm,LAIinit

      end subroutine read_cohort_struct

!************************************************************************

      subroutine read_cohort_struct_old ( cop, iu_entstruct )
      type(cohort) :: cop
      integer :: iu_entstruct
      !---

      call skipstar(iu_entstruct)
      read(iu_entstruct,*) cop%pft
      call skipstar(iu_entstruct)
      read(iu_entstruct,*) cop%n,cop%h,cop%crown_dx,cop%crown_dy,
     &     cop%dbh,cop%root_d,cop%clump
!      write(*,*) cop%n,cop%h,cop%crown_dx,cop%crown_dy,
!     &     cop%dbh,cop%root_d,cop%clump
!      call skipstar(iu_entstruct)
!      read(iu_entstruct,*) cop%C_fol,cop%N_fol,cop%C_sw,cop%N_sw,
!     &     cop%C_hw,cop%N_hw,cop%C_lab,cop%N_lab,cop%
!     &     C_froot,cop%N_froot,cop%C_croot,cop%N_croot

      end subroutine read_cohort_struct_old

!************************************************************************

      subroutine get_cohort_struct ( iu_entstruct, pp )
!@sum get_cohort_struct  For mixed canpopes, read csv file for initialization
!@+     of cohort structure, insert cohort, write initialization values.
      use ent_prescribed_drv, only:  init_canopy_physical_single
      implicit none
      integer :: iu_entstruct
      type(patch),pointer :: pp 
      !-----------
      integer :: pft
      real*8 :: LAImax, hm, LAIinit
      real*8 :: dbh,npop,cradx, crady
      real*8 :: nm, fracroot(N_DEPTH)
      real*8 :: cpool(N_BPOOLS) !g-C/pool/plant
      real*8 :: Ci, GCOHORT
      real*8 :: Tcan, Qf
      real*8 :: albedo(N_BANDS)

      if (.not.ASSOCIATED(pp)) then
         call stop_model("ent_make_struct: null patch",255)
      else
         call read_cohort_struct( iu_entstruct, pft,LAImax,hm,LAIinit )

         call calc_cohort_allometry_single (pft,LAImax,hm,LAIinit
     o        ,dbh,npop,cradx,crady,cpool,nm,fracroot)

         call init_canopy_physical_single(Ci, GCOHORT, Tcan, Qf)

         call insert_cohort(pp,pft,npop, hm,
     &        nm,LAIinit,
     &        cradx, crady,dbh, 0.d0, 0.d0, 0.d0,fracroot,
     &        cpool(FOL), 0.d0, cpool(SW), 0.d0, cpool(HW), 0.d0, 
     &        cpool(LABILE), 0.d0,
     &        cpool(FR), 0.d0, cpool(CR), 0.d0,
     &        Ci, GCOHORT,0.d0, 0.d0, 0.d0, 0.d0, 
     &        0.d0, 0.d0,
              !phenofactor_c, phenofactor_d, phenofactor, phenostatus,
     &        0.d0, 1.d0,0.d0,1.d0, !KIM - starting in the winter for cold-dec.
     &        0.d0, 0.d0,       !betad_10d, CB_d
     &        1.d0, -999.d0)    !turnover_amp, llspan

         write(*,*) 'cohort pft, height',pft, hm 
      endif

      end subroutine get_cohort_struct

!************************************************************************

      subroutine calc_cohort_allometry_single (pft,LAImax,hinit,LAIinit
     o     ,dbh,npop,cradx,crady,cpool,nm,fracroot)
!@sum calc_cohort_allometry.  
!+    Calculates cohort density, allometry and biomass pools.      
      !use ent_prescr_veg, only :  popdensity
      use ent_prescr_veg, only :  prescr_calc_rootprof
      use allometryfn, only : height2dbh, Cfol_fn, nplant
     &    , crown_radius_horiz_allom, crown_radius_vert
     &    , allom_plant_cpools, init_Clab
     &     ,Crown_rad_max_from_density, Crown_rad_allom
      use ent_pfts, only : pfpar,COVEROFFSET,nmv,form
      implicit none
      integer, intent(in) :: pft
      real*8, intent(in) :: LAImax, hinit, LAIinit
      real*8, intent(out) :: dbh,npop,cradx, crady
      real*8, intent(out) :: cpool(N_BPOOLS) !g-C/pool/plant
      real*8, intent(out) :: nm, fracroot(N_DEPTH)
      !---Local------
      !real*8 :: albedo(N_BANDS)
      real*8 :: LAmax, htop

!!      !*Follows ent_prescribed_drv.f:prescr_get_ent_plant, non-Matthews, non-geo
!!      dbh = height2dbh(pft,hm)
!!      !cop%n = popdensity(pft,dbh,LAImax) !Used in prescr_calc_canopy_geometry Matthews
!!      npop = nplant(pft,dbh,hm,LAImax)
!!      !cop%dx = Crown_rad_max_from_density(n)  !This is used in prescr_calc_canopy_geometry Matthews
!!      cradx = min(Crown_rad_max_from_density(npop)
!!      ,Crown_rad_allom(cop%pft,hm))
!!      !crady = crown_radius_vert(hm,cradx) !Not assigned in ent_prescribed_drv.f!
!!      call allom_plant_cpools(pft,LAIinit,hm,dbh,npop,cpool(:))
!!      call init_Clab(pft,dbh,hm,cpool(LABILE))

      if (form(pft).eq.HERB) then
         htop = 1.5d0 !Max height for grasses to calculate npop.
                     !## When mortality/growth are run, then npop will be by
                     !actual height and obey self-thinning law.
                     !## Should eventually be pft-specific. Current Ent GVSD:
                     !annual grass 0.5 m, C4 and perennial grass 1.5 m, crops
                     !herb 0.5 m
                     !Literature: wheat 2 m, corn 3+ m (record is 10 m).
      else !tree and shrubs
         htop = hinit
      endif

      !*Follows ent_prescribed_drv.f: prescr_vegdata, init_entvegdata_geo
      if ((htop.gt.0.d0).and.(LAImax.gt.0.d0)) then
        dbh = height2dbh(pft,htop)
        LAmax =  0.001d0 * Cfol_fn(pft,dbh,htop) *pfpar(pft)%sla !gC to kgC for sla in m2/kgC
        if (LAmax.gt.0.d0) then
           npop = LAImax/LAmax
        else
           npop = 0.d0
        endif
        !cradx = min(Crown_rad_max_from_density(npop)
        !   ,Crown_rad_allom(pft,hm))
        write(*,*) 'Crown_rad_max_from dens,LAImax, LAmax'
     &       ,Crown_rad_max_from_density(npop),LAImax, LAmax
        write(*,*) 'Crown_rad_allom',Crown_rad_allom(pft,hinit)
        cradx = crown_radius_horiz_allom(pft,hinit,npop)
        crady = crown_radius_vert(pft,hinit,cradx)
        call allom_plant_cpools(pft,LAIinit,hinit,dbh,npop,cpool)
        call init_Clab(pft,dbh,hinit,cpool(LABILE))
        nm = nmv(pft + COVEROFFSET)
        call prescr_calc_rootprof(fracroot,pft + COVEROFFSET)
      else 
        dbh = 0.d0
        npop = nplant(pft,dbh,htop,LAImax) !If no height, default.
        cradx = 0.d0
        crady = 0.d0
        cpool(:) = 0.d0
        nm = 0.d0
        fracroot(:) = 0.d0
      endif
       !cop%LAI = No need to assign LAI here, because max is used for allometry.

      end subroutine calc_cohort_allometry_single
!************************************************************************
      subroutine calc_cohort_allometry_old (cop)
!@sum calc_cohort_allometry.  cop comes initialized with pft, n, h.
!+    This subroutine calculates other allometry and biomass pools.      
      !use ent_prescr_veg, only : crown_radius_hw
      use ent_prescr_veg, only :  prescr_calc_rootprof
      use allometryfn, only : Crown_rad_allom, height2dbh
     &    , allom_plant_cpools, init_Clab
      use ent_pfts, only : COVEROFFSET,nmv
      implicit none
      type(cohort),pointer :: cop
      !---Local------
      real*8 :: cpool(N_BPOOLS) !g-C/pool/plant

      cop%dbh = height2dbh(cop%pft,cop%h)
      !cop%crown_dx = crown_radius_hw(cop%dbh) !## Eventually need to make pft-specific
      cop%crown_dx = Crown_rad_allom(cop%pft,cop%h)
      !cop%crown_dy = 0.66d0*cop%h !* Temporary estimate
      cop%nm = nmv(cop%pft+COVEROFFSET)
      !if (.not.FORCE_VEG)
      !cop%LAI = No need to assign LAI here, because max is used for allometry.
      call allom_plant_cpools(cop%pft,0.d0,cop%h,cop%dbh,cop%n,cpool)
      !if .not.FORCE_INIT_CLAB 
      call init_Clab(cop%pft,cop%dbh,cop%h,cpool(LABILE))
      !else *ent_struct_readcsv should provide Clab
      call prescr_calc_rootprof(cop%fracroot,cop%pft + COVEROFFSET)
      cop%C_fol = cpool(FOL)
      cop%C_sw = cpool(SW)
      cop%C_hw = cpool(HW)
      cop%C_lab = cpool(LABILE) 
      cop%C_froot = cpool(FR)
      cop%C_croot = cpool(CR)
      
      end subroutine calc_cohort_allometry_old
!************************************************************************

      subroutine ent_struct_readcsv (ec, iu_entstruct)
!@sum Read ent vegetation structure from ASCII CSV file
!@sum Called by do_ent_struct or ent_prog:
!@sum    - do_ent_struct loops through entcells and checks ij bounds
!@sum    - ent_struct_read_csv loops through reading patches and cohorts
!@sum Order in CSV file must be:
!@sum First line:  N_CASA_LAYERS
!@sum e, entcells in any order spanning dimension IM,JM
!@sum p, patches in order from oldest to youngest
!@sum c, cohorts any order (insert_cohort sorts by height)
!@sum *, comment line
!@sum $, end of entcell
!@sum #, end of file

      implicit none
      type(entcelltype) :: ec
      integer :: iu_entstruct      
      !----
      type(patch),pointer :: pp, tpp
      type(cohort),pointer :: cop,tcp
      character :: next
      integer :: counter, pk
      logical :: end_of_entcell

      !* Set up pointers.
!      call entcell_construct(ec)
!      call zero_entcell( ec )
      !!call patch_construct(pp,null(),1.d0,2)!Blank patch to hold values.
      !call patch_construct(pp,ec,1.d0,2)!Blank patch to hold values.
      !call zero_patch(pp)
      !call cohort_construct(cop)
      !call zero_cohort(cop)
      
      nullify(pp)  !define init
      nullify(cop) !define init
    
      counter = 0
      pk = 0
      end_of_entcell = .false.
      do while (.not.end_of_entcell) 
        counter = counter + 1
        if ( counter >= 1000) 
     &       call stop_model("ent_readcsv: infinite loop",255)
        read(iu_entstruct,*) next
        if (next.eq.'$') then
           write(*,*) 'End of entcell.'
           end_of_entcell = .true.
           !exit
        else if (next.eq.'*') then !skip comment
        else if (next.eq.'p') then !new patch
          write(*,*) 'p'
          pk = pk + 1
          write(*,*) 'pk = ', pk
          if (pk.eq.1) then !First patch overwrites initial dummy.
             write(*,*) 'ent_make_struct: first patch'
             call get_patch_struct( iu_entstruct, 1, ec)
          else !New patch
             write(*,*) 'ent_make_struct: new patch'
             call get_patch_struct( iu_entstruct, 2, ec)
             write(*,*) 'inserted patch'
          endif
             !write(*,*) 'patch age area soil',pp%age,pp%area,pp%soil_type
          pp=>ec%youngest
        else if (next.eq.'c') then !new cohort
          write(*,*) 'c'
          !## Assumes last patch is youngest 
          !## (insert_patch does not currently sort by age).
          call get_cohort_struct( iu_entstruct, pp )
          write(*,*) 'inserted cohort'
        end if
      end do

      call summarize_entcell(ec)

      !## Check entcell structure ----------
      write(*,*) 'ENTSTRUCT entcell summary -'
      write(*,*) 'entcell: area sum C_w =',ec%area,ec%C_w
      tpp => ec%oldest
      do while (associated(tpp))
         write(*,*) '  patch: area soil_type =',tpp%area,tpp%soil_type
         tcp=>ec%oldest%tallest
         do while (associated(tcp))
            write(*,*) '  cohort: pft=', tcp%pft,' h=',tcp%h
            tcp=>tcp%shorter
         enddo
         tpp=>tpp%younger
      enddo
      !##-----------------------------------


      end subroutine ent_struct_readcsv


      end module ent_make_struct
