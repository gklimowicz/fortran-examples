      module cohorts
!@sum Routines to organize cohorts within an entcell.

      use ent_const
      use ent_types

      implicit none


      contains
      !*********************************************************************
      subroutine insert_cohort(pp,pft,n, h,
     &     nm,LAI,
     &     crown_dx, crown_dy,dbh, clump,LMA, root_d,fracroot,
     &     C_fol, N_fol, C_sw, N_sw, C_hw, N_hw, C_lab, N_lab,
     &     C_froot, N_froot, C_croot, N_croot,
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     N_up, C_to_Nfix, 
     &     phenofactor_c, phenofactor_d, phenofactor, phenostatus, 
     &     betad_10d, CB_d,
     &     turnover_amp, llspan) !KIM -7 vars for phenology
!@sum insert_cohort Insert new cohort into a canopy patch.

      type(patch),pointer :: pp
      integer :: pft
      real*8, intent(in) :: n
      real*8, optional, intent(in) :: h, nm, LAI,
     &     crown_dx, crown_dy,dbh, clump
      real*8, optional, intent(in) :: root_d,fracroot(N_DEPTH)
      real*8, optional, intent(in) :: LMA, C_fol, N_fol, C_sw,
     &     N_sw, C_hw, N_hw,
     &     C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     N_up, C_to_Nfix,
     &     phenofactor_c, phenofactor_d, phenofactor, phenostatus, 
     &     betad_10d, CB_d,
     &     turnover_amp, llspan
!     &     stressH2O, stressH2O(N_DEPTH) !No need to assign biophysical values initialized in cohort_construct.
      !------------------
      type(cohort),pointer :: cop, csp, newc

      !If !ASSOCIATED(pp) then ERROR

      !* If pp has cohorts, then insert, else allocate.
      !* Insert at correct height in list, starting from shortest,
      !* since most new cohorts will be shortest.
        !!ALLOCATE(newc)
        !!newc%pptr = pp
      call cohort_construct(newc, pp, pft)
      newc%n = n

      if ( present(h) ) then
        call assign_cohort(newc,pft,n, h, nm, LAI,
     &       crown_dx, crown_dy,dbh, clump, LMA, root_d,fracroot,
     &       C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &       C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &       Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &       N_up, C_to_Nfix,
     &       phenofactor_c, phenofactor_d, phenofactor,phenostatus,
     &       betad_10d, CB_d,
     &       turnover_amp, llspan)

        newc%Ntot = nm*LAI
      endif

        if (ASSOCIATED(pp%shortest)) then !A. There are other cohorts.
          if (pp%shortest%h.ge.newc%h) then !newc is shortest
            pp%shortest%shorter => newc  !changed = to => -PK 9/28/07 
            newc%taller => pp%shortest
            pp%shortest => newc
          else if (pp%tallest%h.lt.newc%h) then !newc is tallest
            pp%tallest%taller => newc  
            newc%shorter => pp%tallest
            pp%tallest => newc
          else !newc is neither tallest nor shortest
            cop => pp%shortest
            do while (cop%h.lt.newc%h) !find next taller cohort
              cop => cop%taller
            end do
            newc%taller => cop
            newc%shorter => cop%shorter
            newc%shorter%taller => newc
            cop%shorter => newc
          end if
          !Now parse through csp's
          csp => newc%taller
          if (ASSOCIATED(csp)) then
            do while (csp%pft.ne.newc%pft)
              if (ASSOCIATED(csp%taller).and.(csp%pft.ne.newc%pft)) then
                csp => csp%taller
              else
                exit !exit loop
              end if
            end do
            if (csp%pft.eq.newc%pft) then !1.newc is not tallest csp
              newc%csptaller => csp
              newc%cspshorter => csp%cspshorter
              csp%cspshorter => newc
            else !2. no taller con-specifics
              nullify(newc%csptaller)
            end if
          else  !3. no taller con-specifics
            nullify(newc%csptaller)
          end if
          if (.NOT.ASSOCIATED(newc%cspshorter)) then !Case 1 did not hold
            csp => newc%shorter
            if (ASSOCIATED(csp)) then
              do while (csp%pft.ne.newc%pft)
                if (ASSOCIATED(csp%shorter).and.
     &               (csp%pft.ne.newc%pft)) then
                  csp => csp%shorter
                else
                  exit
                end if
              end do
              if (csp%pft.eq.newc%pft) then !4. newc is not shortest csp
                newc%cspshorter => csp
                newc%csptaller => csp%csptaller
                csp%csptaller => newc
              else !5. no shorter con-specifics
                nullify(newc%cspshorter)
              end if
            else !6. no shorter con-specifics
              nullify(newc%cspshorter)
            end if
          end if
        else !B. newc is the only cohort
          pp%tallest => newc
          pp%shortest => newc
          ! actually following are already nullified in cohort_construct
          !nullify(newc%taller) 
          !nullify(newc%shorter)
          !nullify(newc%csptaller)
          !nullify(newc%cspshorter)
        end if

      end subroutine insert_cohort
      !*********************************************************************
      
      subroutine assign_cohort(cop,pft,n, h, nm, LAI,
     &     crown_dx, crown_dy,dbh, clump,LMA,root_d,fracroot,
     &     C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &     C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     N_up, C_to_Nfix,
     &     phenofactor_c, phenofactor_d, phenofactor, phenostatus, 
     &     betad_10d, CB_d,
     &     turnover_amp, llspan)
!     &     stressH2O, stressH2Ol)
!@sum assign_cohort  Assign values to a cohort.

      use ent_pfts
      type(cohort) :: cop
      integer :: pft
      real*8,optional :: n, h, nm, LAI,
     &     crown_dx, crown_dy,dbh, root_d, fracroot(:),clump,
     &     LMA, C_fol, N_fol, C_sw, N_sw, C_hw, N_hw,
     &     C_lab, N_lab, C_froot, N_froot, C_croot, N_croot,
     &     Ci, GCANOPY, GPP, NPP, R_auto, R_root,
     &     N_up, C_to_Nfix,
     &     phenofactor_c, phenofactor_d, phenofactor, phenostatus,
     &     betad_10d, CB_d,
     &     turnover_amp, llspan
!     &     stressH2O, stressH2Ol(:)
      cop%pft = pft
      cop%n = n
      cop%nm = nm
      cop%LAI = LAI
!      write(777,*) __FILE__,__LINE__,cop%LAI
      cop%h = h
      cop%crown_dx =  crown_dx 
      cop%crown_dy =  crown_dy
      cop%dbh =  dbh 
      cop%root_d = root_d
      cop%fracroot(:) = fracroot(:)
      cop%clump = clump
      cop%LMA =  LMA 
      cop%C_fol =  C_fol 
      cop%N_fol = N_fol
      cop%C_sw =  C_sw 
      cop%N_sw =  N_sw 
      cop%C_hw =  C_hw 
      cop%N_hw = N_hw
      cop%C_lab =  C_lab 
      cop%N_lab =  N_lab 
      cop%C_froot = C_froot
      cop%N_froot =  N_froot 
      cop%C_croot =  C_croot 
      cop%N_croot = N_croot
      cop%Ci = Ci
      cop%GCANOPY =  GCANOPY 
      cop%GPP =  GPP 
      cop%NPP =  NPP 
      cop%R_auto = 0.0
      cop%R_root = 0.0 !PK 5/15/07
      cop%N_up =  N_up 
!      cop%C_litter = C_litter
!      cop%N_litter = N_litter
      cop%C_to_Nfix = C_to_Nfix
      cop%phenofactor_c = phenofactor_c
      cop%phenofactor_d = phenofactor_d
      cop%phenofactor = phenofactor
      cop%phenostatus = phenostatus
      cop%betad_10d = betad_10d
      cop%CB_d = CB_d
      cop%turnover_amp = turnover_amp
      cop%llspan = llspan
      if (cop%llspan.eq.-999.d0 .and.
     &   (pfpar(cop%pft)%phenotype.eq.EVERGREEN).and.  
     &   (pfpar(cop%pft)%leaftype.eq.BROADLEAF))
     &   cop%llspan=pfpar(cop%pft)%lrage*12.d0
!      cop%stressH2O = stressH2O      !Calculated in biophysics
!      cop%stressH2Ol = stressH2Ol    !Calculated in biophysics
      !* diags and hacks
      cop%C_total = 0.d0
      cop%C_growth = 0.d0
      cop%C_growth_flux = 0.d0

      end subroutine assign_cohort
      !*********************************************************************
!!! basically replaced with cohort_construct()
cddd      subroutine init_cohort_defaults(cop,pnum)
cddd!     @sum Initialize a cohort with default values
cddd      type(cohort),pointer :: cop
cddd      integer :: pnum
cddd
cddd      cop%pft = pnum
cddd      cop%n = 1                 !## Dummy ##!
cddd!     cop%pptr = pp
cddd!     call nullify(cop%taller) !Only one cohort
cddd!     call nullify(cop%shorter)
cddd!     call nullify(cop%csptaller)
cddd!     call nullify(cop%cspshorter)            
cddd
cddd      call zero_cohort(cop)
cddd      
cddd      end subroutine init_cohort_defaults

      subroutine zero_cohort(cop)
!@sum Zero all real variables in cohort record.      
      use ent_pfts
      type(cohort),pointer :: cop

      cop%nm = 0.0
      cop%Ntot = 0.0

      !* Individual plant properties *!
      !* GEOMETRY *!
      cop%h = 0.0
      cop%crown_dx = 0.0
      cop%crown_dy = 0.0
      cop%dbh = 0.0
      cop%root_d = 0.0
      cop%LAI = 0.0
      cop%clump = 0.0
      cop%fracroot(:) = 0.0

      !* BIOMASS POOLS *!
      cop%LMA = 0.0
      cop%C_fol = 0.0
      cop%N_fol = 0.0
      cop%C_sw = 0.0
      cop%N_sw = 0.0
      cop%C_hw = 0.0
      cop%N_hw = 0.0
      cop%C_lab = 0.0
      cop%N_lab = 0.0
      cop%C_froot = 1.0         !Dummy
      cop%N_froot = 0.0
      cop%C_croot = 0.0
      cop%N_croot = 0.0
      
      !* FLUXES *!
      cop%Ci =  0.0127D0        !Initial value not zero.
      cop%GCANOPY = 0.0
      cop%GPP = 0.0
      cop%IPP = 0.0
      cop%NPP = 0.0
      cop%R_auto = 0.0
      cop%R_root = 0.0 !PK 5/15/07
      cop%N_up = 0.0
!      cop%C_litter = 0.0
!      cop%N_litter = 0.0
      cop%C_to_Nfix = 0.0

      !* PHENOLOGY/GROWTH *!
      !KIM - starting in the middle of winter for cold-dec.
      cop%phenofactor_c = 0.d0
      cop%phenofactor_d = 1.d0
      cop%phenofactor = 0.d0
      cop%phenostatus = 1.d0
      !KIM - starting in the middle of growing season
!      cop%phenofactor_c = 1.d0
!      cop%phenofactor_d = 1.d0
!      cop%phenofactor = 1.d0
!      cop%phenostatus = 3.d0
      cop%betad_10d = 1.d0
      cop%CB_d = 0.d0
      cop%turnover_amp = 1.d0
      cop%llspan = -999.d0
!!!   cop%pft is not known here !!!
!!      if ((pfpar(cop%pft)%phenotype.eq.EVERGREEN).and.  
!!     &   (pfpar(cop%pft)%leaftype.eq.BROADLEAF))
!!     &   cop%llspan=pfpar(cop%pft)%lrage*12.d0
      cop%Sacclim = 25.d0 !NK - force mild average temperatures default.

      !* PHYSIOLOGICAL STATUS *!
      cop%stressH2O = 1.d0 !Default no stress.
      cop%stressH2Ol(:) = 1.d0 !Default no stress.
      cop%senescefrac = 0.d0

      !* diags and hacks
      cop%C_total = 0.d0
      cop%C_growth = 0.d0
      cop%C_growth_flux = 0.d0
      end subroutine zero_cohort


      !*********************************************************************
      !*********************************************************************
       
      subroutine reorganize_cohorts(pp)
!@sum Place holder.
      type(patch),pointer :: pp

      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !Not need in GISS replication test.
      
      end subroutine reorganize_cohorts
      !*********************************************************************


      subroutine cohort_construct(cop, parent_patch, pnum)
!@sum cohort_construct  Create a cohort with default values. if optional values
!@+ are provided - set them
!@auth I.Aleinov.
      ! This function may eventually be combined with assign_cohort
      ! for better performance
      type(cohort),pointer :: cop
      integer, optional :: pnum
      type(patch), optional, target :: parent_patch

      ! allocate memory
      allocate( cop )
      allocate( cop%fracroot(N_DEPTH) )
      allocate( cop%stressH2Ol(N_DEPTH) )

      ! just in case...
      nullify( cop%height_dz )
      nullify( cop%fp_dz )
      nullify( cop%height )
      nullify( cop%fp )

      ! set pointers if any
      nullify(cop%pptr )
      nullify(cop%taller )
      nullify(cop%shorter )
      nullify(cop%csptaller )
      nullify(cop%cspshorter )
      if ( present(parent_patch) ) then
        cop%pptr => parent_patch
      endif

      ! set variables
      cop%pft = -1              ! = -1 if pft not set
      if ( present(pnum) ) cop%pft = pnum
      cop%n = 0.0

      call zero_cohort(cop)

      end subroutine cohort_construct
      !*********************************************************************


      subroutine cohort_destruct(cop)
!@sum cohort_destruct Deallocate memory used by cohort
!@auth I.Aleinov.
      type(cohort),pointer :: cop

      ! we may want ot collapse hole between "taller" and "shorter"
      ! here if this functionality is needed

      ! deallocate all memory
      deallocate( cop%stressH2Ol )
      deallocate( cop%fracroot )
      deallocate( cop )
      nullify( cop )

      end subroutine cohort_destruct


      subroutine cohort_print_old(iu, cop, prefix)
!@sum Debug routine to print contents of a cohort.
      integer, intent(in) :: iu
      type(cohort), intent(in) :: cop
      character*(*), optional, intent(in) :: prefix
      !---
      integer n

      write(iu, '(a,a," = ",i7)') prefix,"pft ",cop%pft
      write(iu, '(a,a," = ",f10.7)') prefix,"n   ",cop%n
      write(iu, '(a,a," = ",f10.7)') prefix,"nm  ",cop%nm
! Ntot doesn't seem terribly important ( always LAI*nm? )
      !write(iu, '(a,a," = ",f10.7)') prefix,"Ntot",cop%Ntot
      write(iu, '(a,a," = ",f10.7)') prefix,"LAI ",cop%LAI
      write(iu, '(a,a," = ",f10.7)') prefix,"h   ",cop%h
      write(iu, '(a,a," = ",f10.7)') prefix,"dbh ",cop%dbh
      write(iu, '(a,"fracroot = " )') prefix
      do n=1,N_DEPTH
        write(iu, '(a,"      ",f10.7)') prefix,cop%fracroot(n)
      enddo
      write(iu, '(a,a," = ",f10.7)') prefix,"Gcan",cop%gcanopy
      write(iu, '(a,a," = ",f10.7)') prefix,"GPP ",cop%GPP

      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_fol",  cop%C_fol
      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_froot",cop%C_froot
      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_hw",   cop%C_hw

      write(iu, '(a,a," = ",f10.7)') prefix,"c%Ci",  cop%Ci

      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_fol",   cop%C_fol   
      write(iu, '(a,a," = ",f10.7)') prefix,"c%N_fol",   cop%N_fol   
      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_sw",    cop%C_sw    
      write(iu, '(a,a," = ",f10.7)') prefix,"c%N_sw",    cop%N_sw    
      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_hw",    cop%C_hw    
      write(iu, '(a,a," = ",f10.7)') prefix,"c%N_hw",    cop%N_hw    
      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_lab",   cop%C_lab   
      write(iu, '(a,a," = ",f10.7)') prefix,"c%N_lab",   cop%N_lab   
      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_froot", cop%C_froot 
      write(iu, '(a,a," = ",f10.7)') prefix,"c%N_froot", cop%N_froot 
      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_croot", cop%C_croot 
      write(iu, '(a,a," = ",f10.7)') prefix,"c%N_croot", cop%N_croot 
                                                             
      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_growth",cop%C_growth
      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_growth_flux"
     &     ,cop%C_growth_flux
      write(iu, '(a,a," = ",f10.7)') prefix,"c%C_total", cop%C_total 

      write(iu, '(a,a," = ",f10.7)') prefix,"c%llspan",  cop%llspan
      write(iu, '(a,a," = ",f10.7)') prefix,"turnover_amp",
     &     cop%turnover_amp


      end subroutine cohort_print_old
      
      subroutine cohort_print(iu, cop, prefix)
!@sum Debug routine to print contents of a cohort.
      integer, intent(in) :: iu
      type(cohort), intent(in) :: cop
      character*(*), optional, intent(in) :: prefix
      !---
      integer n

      write(iu,'(a,a," = ",i7)') prefix,"pft ",cop%pft
      write(iu,1)prefix,"n   ",cop%n
      write(iu,1)prefix,"nm          ",cop%nm            
      write(iu,1)prefix,"Ntot        ",cop%Ntot          
      write(iu,1)prefix,"LAI   	     ",cop%LAI 
      write(iu,1)prefix,"LMA         ",cop%LMA           
      write(iu,1)prefix,"h      	 ",cop%h      	 
      write(iu,1)prefix,"crown_dx  	 ",cop%crown_dx  	 
      write(iu,1)prefix,"crown_dy  	 ",cop%crown_dy  	 
      write(iu,1)prefix,"dbh         ",cop%dbh           
      write(iu,1)prefix,"root_d      ",cop%root_d        
      write(iu,1)prefix,"clump       ",cop%clump         
      write(iu,1)prefix,"fracroot(:) ",
     &     cop%fracroot(:)	    
      if (associated(cop%height_dz) )
     &     write(iu,1)prefix,"height_dz(:)",
     &     cop%height_dz(:)  
      if ( associated(cop%fp_dz) )
     & write(iu,1)prefix,"fp_dz(:)    ",cop%fp_dz(:)      
      if ( associated(cop%height) )
     &write(iu,1)prefix,"height(:)   ",cop%height(:)     
      if ( associated(cop%fp) )
     &   write(iu,1)prefix,"fp(:)       ",cop%fp(:)         
      write(iu,1)prefix,"C_fol       ",cop%C_fol         
      write(iu,1)prefix,"N_fol       ",cop%N_fol         
      write(iu,1)prefix,"C_sw        ",cop%C_sw          
      write(iu,1)prefix,"N_sw        ",cop%N_sw          
      write(iu,1)prefix,"C_hw        ",cop%C_hw          
      write(iu,1)prefix,"N_hw        ",cop%N_hw          
      write(iu,1)prefix,"C_lab       ",cop%C_lab
      write(iu,1)prefix,"N_lab       ",cop%N_lab         
      write(iu,1)prefix,"C_froot     ",cop%C_froot       
      write(iu,1)prefix,"N_froot     ",cop%N_froot       
      write(iu,1)prefix,"C_croot     ",cop%C_croot       
      write(iu,1)prefix,"N_croot     ",cop%N_croot       
      write(iu,1)prefix,"Ci          ",cop%Ci            
      write(iu,1)prefix,"gcanopy     ",cop%gcanopy       
      write(iu,1)prefix,"GPP         ",cop%GPP           
      write(iu,1)prefix,"IPP         ",cop%IPP           
      write(iu,1)prefix,"NPP         ",cop%NPP           
      write(iu,1)prefix,"R_auto      ",cop%R_auto        
      write(iu,1)prefix,"R_root      ",cop%R_root        
      write(iu,1)prefix,"N_up        ",cop%N_up          
      write(iu,1)prefix,"C_to_Nfix   ",cop%C_to_Nfix     
      write(iu,1)prefix,"phenofactor ",
     &     cop%phenofactor   
      write(iu,1)prefix,"phenofactorc",
     &     cop%phenofactor_c  
      write(iu,1)prefix,"phenofactord",
     &     cop%phenofactor_d  
      write(iu,1)prefix,"phenostatus ",
     &     cop%phenostatus   
      write(iu,1)prefix,"betad_10d	 ",cop%betad_10d	 
      write(iu,1)prefix,"CB_d	 ",cop%CB_d	 
      write(iu,1)prefix,"turnover_amp",
     &     cop%turnover_amp  
      write(iu,1)prefix,"llspan      ",cop%llspan        
      write(iu,1)prefix,"Sacclim     ",cop%Sacclim    	         
      write(iu,1)prefix,"stressH2O   ",cop%stressH2O
      if ( associated(cop%stressH2Ol) )
     &     write(iu,1)prefix,"stressH2Ol(:",
     &     cop%stressH2Ol(:)
      write(iu,1)prefix,"senescefrac ",
     &     cop%senescefrac	         
      write(iu,1)prefix,"C_growth    ",cop%C_growth      
      write(iu,1)prefix,"C_growth_flux",cop%C_growth_flux    
      write(iu,1)prefix,"C_total     ",cop%C_total       
   
 1    format(a,a," = ",99e23.16)  ! e12.5

      end subroutine cohort_print
       
      !*********************************************************************
      subroutine cohort_print_diag(iu, cop, prefix)
!@sum Debug routine to print other calculated diagnostics not output 
!@+   by cohort_print
      !@auth NK
      use ent_pfts
      use allometryfn, only : Cfol_fn
      implicit none
      integer, intent(in) :: iu
      type(cohort), intent(in) :: cop
      character*(*), optional, intent(in) :: prefix
      !---
      integer :: n

      write(iu,'(a,a," = ",i7)') prefix,"pft ",cop%pft
      write(iu,1)prefix,"n   ",cop%n
      write(iu,1)prefix,"LAI   	     ",cop%LAI 
      write(iu,1)prefix,"LMA         ",cop%LMA           
      write(iu,1)prefix,"h      	 ",cop%h      	 
      write(iu,1)prefix,"crown_dx  	 ",cop%crown_dx  	 
      write(iu,1)prefix,"crown_dy  	 ",cop%crown_dy  	 
      write(iu,1)prefix,"dbh         ",cop%dbh           
      write(iu,1)prefix,"Cfolmax     ",Cfol_fn(cop%pft,cop%dbh,cop%h)
      write(iu,1)prefix,"LA max      ",
     &     Cfol_fn(cop%pft,cop%dbh,cop%h)/(pfpar(cop%pft)%sla)
      write(iu,1)prefix,"LAI max     ",
     &     Cfol_fn(cop%pft,cop%dbh,cop%h)/(pfpar(cop%pft)%sla) * cop%n
      
!      write(iu,1)prefix,"fracroot(:) ",
!     &     cop%fracroot(:)	    
      if (associated(cop%height_dz) )
     &     write(iu,1)prefix,"height_dz(:)",
     &     cop%height_dz(:)  
      if ( associated(cop%fp_dz) )
     & write(iu,1)prefix,"fp_dz(:)    ",cop%fp_dz(:)      
      if ( associated(cop%height) )
     &write(iu,1)prefix,"height(:)   ",cop%height(:)     
      if ( associated(cop%fp) )
     &   write(iu,1)prefix,"fp(:)       ",cop%fp(:)         
      write(iu,1)prefix,"C_fol       ",cop%C_fol         
      write(iu,1)prefix,"C_sw        ",cop%C_sw          
      write(iu,1)prefix,"C_hw        ",cop%C_hw          
      write(iu,1)prefix,"C_lab       ",cop%C_lab
      write(iu,1)prefix,"C_froot     ",cop%C_froot       
      write(iu,1)prefix,"C_croot     ",cop%C_croot       

      write(iu,1)prefix,"C_total     ",cop%C_total       
   
 1    format(a,a," = ",99e23.16)  ! e12.5

      end subroutine cohort_print_diag
       
      !*********************************************************************
      subroutine calc_CASArootfrac(copfracroot,fracrootCASA)  !PK 11/06
!@sum calc_CASArootfrac  Maps fracroot(N_DEPTH) to fracrootCASA(N_CASA_LAYERS)
!@+   ifdef customization required dependent on thicknesses of CASA layers
!@+   and GCM layers 
      !type(cohort),intent(in) :: cop
      real*8,pointer :: copfracroot(:)
      real*8,intent(out) :: fracrootCASA(N_CASA_LAYERS)

      if (N_CASA_LAYERS == 1) then
         fracrootCASA = 1.d0  !if there is no explicit depth structure
      else
#ifdef NCASA2
      !***scheme for N_CASA_LAYERS=2 (layers: 0-30, 30-100 cm)*** 
         fracrootCASA(1) = copfracroot(1) + copfracroot(2)  !CASA layer 1 --> GISS GCM layers 1,2
         fracrootCASA(2) = copfracroot(3) + copfracroot(4)  !CASA layer 2 --> GISS layers 3,4
     &                + copfracroot(5)                    !need to add 5th GISS layer (mainly for trees) -PK 6/26/07
#endif
      end if
                               
      end subroutine calc_CASArootfrac
      !*********************************************************************


      subroutine cohort_merge_data( cop1, wp1, cop2, wp2 )
!@sum Cohort management, to merge similar cohorts into one.
      type(cohort)  :: cop1, cop2
      real*8 :: wp1, wp2
      !---
      real*8 :: w1, w2, w_tot
      
      if ( cop1%pft .ne. cop2%pft ) then
        print *,"cohort_merge_data: prf1, pft2=", cop1%pft, cop2%pft
        call stop_model("Can't merge cohorts with different pfts",255)
      endif

      w_tot = wp1*cop1%n + wp2*cop2%n
      w1 = wp1*cop1%n/w_tot
      w2 = wp2*cop2%n/w_tot

      ! values weighted with patch area
      cop1%n            =wp1*cop1%n           +wp2*cop2%n
      
      ! values weighted with (patch area)*n
      cop1%nm           =w1*cop1%nm           +w2*cop2%nm           
      cop1%Ntot         =w1*cop1%Ntot         +w2*cop2%Ntot         
                                                                           
      cop1%h            =w1*cop1%h            +w2*cop2%h            
      cop1%crown_dx     =w1*cop1%crown_dx     +w2*cop2%crown_dx     
      cop1%crown_dy     =w1*cop1%crown_dy     +w2*cop2%crown_dy     
      cop1%dbh          =w1*cop1%dbh          +w2*cop2%dbh          
      cop1%root_d       =w1*cop1%root_d       +w2*cop2%root_d       
      cop1%LAI          =w1*cop1%LAI          +w2*cop2%LAI          
      cop1%clump        =w1*cop1%clump        +w2*cop2%clump        
      cop1%fracroot(:)  =w1*cop1%fracroot(:)  +w2*cop2%fracroot(:)  
                                                                           
      cop1%LMA          =w1*cop1%LMA          +w2*cop2%LMA          
      cop1%C_fol        =w1*cop1%C_fol        +w2*cop2%C_fol        
      cop1%N_fol        =w1*cop1%N_fol        +w2*cop2%N_fol        
      cop1%C_sw         =w1*cop1%C_sw         +w2*cop2%C_sw         
      cop1%N_sw         =w1*cop1%N_sw         +w2*cop2%N_sw         
      cop1%C_hw         =w1*cop1%C_hw         +w2*cop2%C_hw         
      cop1%N_hw         =w1*cop1%N_hw         +w2*cop2%N_hw         
      cop1%C_lab        =w1*cop1%C_lab        +w2*cop2%C_lab        
      cop1%N_lab        =w1*cop1%N_lab        +w2*cop2%N_lab        
      cop1%C_froot      =w1*cop1%C_froot      +w2*cop2%C_froot      
      cop1%N_froot      =w1*cop1%N_froot      +w2*cop2%N_froot      
      cop1%C_croot      =w1*cop1%C_croot      +w2*cop2%C_croot      
      cop1%N_croot      =w1*cop1%N_croot      +w2*cop2%N_croot      
                                                         
      cop1%Ci           =w1*cop1%Ci           +w2*cop2%Ci           
      cop1%GCANOPY      =w1*cop1%GCANOPY      +w2*cop2%GCANOPY      
      cop1%GPP          =w1*cop1%GPP          +w2*cop2%GPP          
      cop1%IPP          =w1*cop1%IPP          +w2*cop2%IPP          
      cop1%NPP          =w1*cop1%NPP          +w2*cop2%NPP          
      cop1%R_auto       =w1*cop1%R_auto       +w2*cop2%R_auto       
      cop1%R_root       =w1*cop1%R_root       +w2*cop2%R_root       
      cop1%N_up         =w1*cop1%N_up         +w2*cop2%N_up         
      cop1%C_to_Nfix    =w1*cop1%C_to_Nfix    +w2*cop2%C_to_Nfix    
                                                                           
      cop1%phenofactor_c=w1*cop1%phenofactor_c+w2*cop2%phenofactor_c
      cop1%phenofactor_d=w1*cop1%phenofactor_d+w2*cop2%phenofactor_d
      cop1%phenofactor  =w1*cop1%phenofactor  +w2*cop2%phenofactor  

      ! phenostatus is int, not sere how to deal with it
      !cop1%phenostatus  =w1*cop1%phenostatus  +w2*cop2%phenostatus  
      cop1%betad_10d    =w1*cop1%betad_10d    +w2*cop2%betad_10d  
      cop1%CB_d         =w1*cop1%CB_d         +w2*cop2%CB_d         
      cop1%turnover_amp =w1*cop1%turnover_amp +w2*cop2%turnover_amp 
      cop1%llspan       =w1*cop1%llspan       +w2*cop2%llspan       
      cop1%Sacclim      =w1*cop1%Sacclim      +w2*cop2%Sacclim      
                                                                           
      cop1%stressH2O    =w1*cop1%stressH2O    +w2*cop2%stressH2O    
      cop1%stressH2Ol(:)=w1*cop1%stressH2Ol(:)+w2*cop2%stressH2Ol(:)
      cop1%senescefrac  =w1*cop1%senescefrac  +w2*cop2%senescefrac  
                                                                           
      cop1%C_total      =w1*cop1%C_total      +w2*cop2%C_total      
      cop1%C_growth     =w1*cop1%C_growth     +w2*cop2%C_growth    
      cop1%C_growth_flux=w1*cop1%C_growth_flux+w2*cop2%C_growth_flux    
 
      end subroutine cohort_merge_data

!----------------------------------------------------------------------
      real*8 function cohort_carbon(cop) 
!@sum Return kg-C/individual in a cohort.
      type(cohort),pointer :: cop
      !---Local------
      real*8 :: kgC_indiv             !kg-C per individual

      kgC_indiv = 0.d0
      if (ASSOCIATED(cop)) then
         kgC_indiv = 0.001d0 * (cop%C_fol + cop%C_sw + cop%C_hw + 
     &       cop%C_lab + cop%C_froot + cop%C_croot )
      endif

      cohort_carbon = kgC_indiv

      end function cohort_carbon
!----------------------------------------------------------------------

      end module cohorts
