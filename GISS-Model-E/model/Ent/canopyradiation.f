!This file contains two modules, qsort_c_module and canopyrad.

      module qsort_c_module
!@sum Module with routines to sort a series in descending order
!@auth W.Yang

#define DEBUG

      implicit none
      public :: QsortC
      private :: Partition

      contains

      recursive subroutine QsortC(A)
      real*8, intent(in out), dimension(:) :: A
      integer :: iq

      if(size(A) > 1) then
         call Partition(A, iq)
         call QsortC(A(:iq-1))
         call QsortC(A(iq:))
      endif
      end subroutine QsortC

      subroutine Partition(A, marker)
      real*8, intent(in out), dimension(:) :: A
      integer, intent(out) :: marker
      integer :: i, j
      real*8 :: temp
      real*8 :: x      ! pivot point
      x = A(1)
      i= 0
      j= size(A) + 1

      do
        j = j-1
        do
          if (A(j) >= x) exit
          j = j-1
        end do
        i = i+1
        do
          if (A(i) <= x) exit
          i = i+1
        end do
        if (i < j) then
          ! exchange A(i) and A(j)
          temp = A(i)
          A(i) = A(j)
          A(j) = temp
        elseif (i == j) then
          marker = i+1
          return
        else
          marker = i
          return
        endif
      end do

      end subroutine Partition

      end module qsort_c_module

      module canopyrad
!@sum Routines for calculating canopy radiation and albedo with the
!@+   Analytical Clumped Two-Stream (ACTS) model (Ni-Meister et al. 2010).
!@+   Is called by canopygort module to provide light for photosynthesis.
!@auth W.Ni-Meister
!@+   UNDER DEVELOPMENT

      !Ent MODULES TO USE
      use ent_const
      use ent_types
      use patches, only : patch_print
      implicit none
      private
      save

      public recalc_radpar_cell, recalc_radpar
      public get_canopy_rad, GORT_clumping, TwoStream       
      public :: gort_input
      public :: pft_pstate_type

      ! type definition
      type gort_input
        real*8 :: dens_tree         ! crown count density (m-2)
        real*8 :: dens_foliage      ! foliage area volume density of a single crown(m-1)
        real*8 :: h1, h2            ! lower and upper bound of crown centers (m)
        real*8 :: delta_z           ! height increase for each height level
        real*8 :: horz_radius       ! horizontal crown radius (m)
        real*8 :: vert_radius       ! vertical crown radius (m)
        real*8 :: zenith            ! solar zenith angle (rad)
        real*8 :: LAI               ! leaf area index (m2m-2)
        real*8 :: dbh               ! stem diameter at breast height (m)
        integer :: pft              ! plant function type
      end type gort_input

      type profile_params
        real*8, dimension(:), pointer :: height_levels ! height level for each profile
        real*8, dimension(:), pointer :: fp            ! foliage profile, no clump involved
        real*8, dimension(:), pointer :: rdfp          ! clumping factor * foliage profile
        real*8, dimension(:), pointer :: rifp          ! clumping factor * foliage profile
        real*8, dimension(:), pointer :: efp           ! effective foliage profile for each profile
        real*8, dimension(:), pointer :: vz            ! trunk volumn
      end type profile_params

      !----------------------------------------------------
      ! pft physical state variables structure
      !----------------------------------------------------
      type pft_pstate_type
        real*8, pointer :: albd(:,:)    !surface albedo (direct)                       (N_BANDS)
        real*8, pointer :: albi(:,:)    !surface albedo (indirect)                     (N_BANDS)
        real*8, pointer :: fabd(:,:)    !flux absorbed by veg per unit direct flux     (N_BANDS)
        real*8, pointer :: fabi(:,:)    !flux absorbed by veg per unit diffuse flux    (N_BANDS)
        real*8, pointer :: ftdd(:,:)    !down direct flux below veg per unit dir flx   (N_BANDS)
        real*8, pointer :: ftid(:,:)    !down diffuse flux below veg per unit dir flx  (N_BANDS)
        real*8, pointer :: ftii(:,:)    !down diffuse flux below veg per unit dif flx  (N_BANDS)
        real*8, pointer :: frid(:,:)    !up diffuse flux below veg per unit dir flx    (N_BANDS)
        real*8, pointer :: frii(:,:)    !up diffuse flux below veg per unit dif flx    (N_BANDS) 
        real*8, pointer :: gdir(:)      !leaf projection in solar direction (0 to 1)
        real*8, pointer :: omega(:,:)   !fraction of intercepted radiation that is scattered (0 to 1)
        real*8, pointer :: isha(:,:)    !shaded flux per unit incident flux            (N_BANDS)         
        real*8, pointer :: isun(:,:)    !sunlit flux per unit incident flux            (N_BANDS)         
      end type pft_pstate_type

      real*8 :: K
      !real*8, parameter :: PI=3.14159265
      real*8, parameter :: ang_dif = 1.0071   ! 57.7 deg in rad
      real*8, parameter :: dz_gin = 0.1       ! init delta z for profile
      real*8, parameter :: lai_thres = 0.001    ! don't calculate if lower
      ! integer, parameter :: N_BANDS      =   2   ! number of solar radiation bands: vis, nir
      integer, parameter :: lbp = 1
      integer :: ubp 
      integer :: num_vegsol 
      ! integer, parameter :: N_PFT = 16
 
      contains
      !*********************************************************************
      
      subroutine recalc_radpar_cell(pptr)
!@sum Calculate canopy radiation geometrical parameters following structural
!@sum changes.  At patch level.  PLACE HOLDER.

      type(entcelltype) :: pptr

!!! passing actual structure instead of pointer, change if necessary
!!!      if (ASSOCIATED(pptr)) then 
        !call get_patchalbedo(tt,pptr) !*This is called by summarize_patch
        !call GORT_clumping(pptr)
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
!!!      end if

      end subroutine recalc_radpar_cell


      subroutine recalc_radpar(pptr)
!@sum Calculate canopy radiation geometrical parameters following structural
!@sum changes.  At patch level.

      type(patch) :: pptr

!!! passing actual structure instead of pointer, change if necessary
!!!      if (ASSOCIATED(pptr)) then 
        !call get_patchalbedo(tt,pptr) !*This is called by summarize_patch
        !call GORT_clumping(pptr)
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
!!!      end if

      end subroutine recalc_radpar


      !*********************************************************************
      subroutine recalc_radpar_OLD(entcell)
!@sum Calculate canopy radiation geometrical parameters following structural
!@sum changes.  At entcell level.

      type(entcelltype) :: entcell
      type(patch),pointer :: pptr

      pptr = entcell%youngest

      do while (ASSOCIATED(pptr)) 
        !call GORT_clumping(pptr)
        pptr = pptr%older
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      end do

      end subroutine recalc_radpar_OLD

      !*********************************************************************

      subroutine calc_canopy_rad(pptr, h)
!@sum Get incident light profiles in canopy and return in pptr%crad
      type(patch),pointer :: pptr
      real*8 :: h               !Height in canopy
      !--------------------------


      !Get incident light profiles in canopy and return in crad
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

      end subroutine calc_canopy_rad
      
      !*********************************************************************
      subroutine get_canopy_rad(pptr, IPAR, Id)
      !ACTS canopy radiation.
      !UNDER DEVELOPMENT
      use ent_pfts, only : rhol, taul, rhos, taus
      type(patch),pointer :: pptr
      ! real*8 :: h               !Height in canopy
      
      real*8 :: IPAR, tempCos 
      real*8 :: Id, Ii
      real*8 :: dz
      real*8, parameter :: mpe = 1.e-06     ! prevents overflow for division by zero

      integer :: i, n_height_level, N2
      integer :: fp,g,c,p                ! indices
      integer :: ib                      ! band index
      integer :: len

!      real*8 :: rhol(N_PFT,N_BANDS), taul(N_PFT,N_BANDS) 
!      real*8 :: rhos(N_PFT,N_BANDS), taus(N_PFT,N_BANDS)
      real*8 :: tran(2)

      real*8, dimension(:), pointer :: ffp, rdfp, rifp, sunlit, shaded
      real*8, dimension(:), pointer :: height_levels, vz, crad_heights
      real*8, dimension(:), pointer :: T_sun, T_sha, I_sun, I_sha
      real*8, dimension(:), pointer :: h_coh
      
      real*8, dimension(:), pointer :: vai, coszen, esai, elai, wl, ws
      real*8, dimension(:,:), pointer :: rho, tau
      integer, dimension(:), pointer :: filter_vegsol, ivt, pcolumn

!      integer :: filter_vegsol(num_vegsol), ivt(num_vegsol) 
!      integer :: pcolumn(num_vegsol)
!      real*8 :: vai(num_vegsol) 
!      real*8 :: rho(num_vegsol,N_BANDS), tau(num_vegsol,N_BANDS)
!      real*8 :: coszen(num_vegsol), esai(num_vegsol), elai(num_vegsol)    
!      real*8 :: wl(lbp:ubp)                 ! fraction of LAI+SAI that is LAI
!      real*8 :: ws(lbp:ubp)                 ! fraction of LAI+SAI that is SAI

      type(pft_pstate_type) :: sout
      type(gort_input), dimension(:), pointer :: gin
    
!      data (ivt(i),i=1,num_vegsol) /1/
!      data (pcolumn(i),i=1,num_vegsol) /num_vegsol*1/    ! 1-soil 2-snow
!      data (elai(i),i=1,num_vegsol) /num_vegsol*1.2/
!      data (esai(i),i=1,num_vegsol) /num_vegsol*1/

!      if (.NOT.ASSOCIATED(pptr)) then
!         print *,'In GORT_clumping: pptr=NULL'
!         stop
!      endif

      Ii = 1 - Id
      if (Ii .eq. 1.d0) then
         tempCos = pptr%cellptr%CosZen !Save input CosZen
         pptr%cellptr%CosZen = 0.5343  !Substitute diffusive angle
      end if 
      if (IPAR.ge.LOW_PAR_LIMIT) then  !Only adjust CosZen here
         !Only call canopy rad if there's light.
         !Check for horizon or night, and set bound, because clumping 
         !calculation uses tan(zenith). Only near-horizon times can pose
         !a problem if there is light when the sun is at or below the horizon.
         !Otherwise, should check for zenith>=PI/2 and zenith<=1.5PI.
      ! More operations are as below
         if (pptr%cellptr%CosZen.lt.0.001d0) then
            pptr%cellptr%CosZen = 0.001d0  
         endif
      endif

      call geo_to_gin(pptr, gin, h_coh)
      if (.NOT.ASSOCIATED(gin)) then ! bare soil or LAI=0
          !print *, 'lai zero'
         call stop_model('GORT_clumping: Bare soil',255)
      endif 

      ubp = size(gin)
      num_vegsol = ubp - lbp + 1

      allocate(filter_vegsol(num_vegsol))
      allocate(ivt(num_vegsol)) 
      allocate(pcolumn(num_vegsol))
      allocate(vai(num_vegsol)) 
      allocate(rho(num_vegsol,N_BANDS))
      allocate(tau(num_vegsol,N_BANDS))
      allocate(coszen(num_vegsol))
      allocate(esai(num_vegsol))
      allocate(elai(num_vegsol))
      allocate(wl(num_vegsol))                 ! fraction of LAI+SAI that is LAI
      allocate(ws(num_vegsol))                 ! fraction of LAI+SAI that is SAI
 
      ! need further update of elai, esai, pcolumn
      do i=1, num_vegsol
         ivt(i) = gin(i)%pft
         pcolumn(i) = 1
         elai(i) = 1.2
         esai(i) = 1
      end do

      ! print *, 'before GORT_clumping...'
      if (IPAR.ge.LOW_PAR_LIMIT) then  !Only call canopy rad if there's light.
         !Check for horizon or night, and set bound, because clumping 
         !calculation uses tan(zenith). Only near-horizon times can pose
         !a problem if there is light when the sun is at or below the horizon.
         !Otherwise, should check for zenith>=PI/2 and zenith<=1.5PI.
         call GORT_clumping(pptr, height_levels, ffp, rdfp, rifp
     &        , sunlit, shaded, vz, gin, h_coh)
         if (.NOT.ASSOCIATED(height_levels)) then ! bare soil or LAI=0
            call stop_model( 'LAI=0',255)  !Trunks
         endif 
         ! print *, 'after GORT_clumping...'

         ! Weight reflectance/transmittance by lai and sai
         ! Only perform on vegetated pfts where coszen > 0

         do i=1, num_vegsol
            filter_vegsol(i) = i
            coszen(i) = pptr%cellptr%CosZen
            p = filter_vegsol(i)
            vai(p) = elai(p) + esai(p)
            wl(p) = elai(p) / max( vai(p), mpe )
            ws(p) = esai(p) / max( vai(p), mpe )
         end do

         do ib = 1, N_BANDS
            do fp = 1,num_vegsol
               p = filter_vegsol(fp)
               rho(p,ib) = max( rhol(ivt(p),ib)*wl(p) + 
     &              rhos(ivt(p),ib)*ws(p), mpe )
               tau(p,ib) = max( taul(ivt(p),ib)*wl(p) + 
     &              taus(ivt(p),ib)*ws(p), mpe )
            end do
         end do

#ifdef ENT_STANDALONE_DIAG
         write(1079,*) lbp,ubp,filter_vegsol(1),num_vegsol,pcolumn(1)
     &        ,coszen(1),Id,size(height_levels),rho(1,1),tau(1,1)
#endif

         call TwoStream(lbp, ubp, filter_vegsol, num_vegsol, ivt, 
     &        pcolumn, coszen, Id, height_levels, rdfp, rifp, rho, tau,
     &        sout)
         
         N2 = size(height_levels)
         do ib = 1, N_BANDS
            pptr%albedo(ib) = Id * sout%albd(1,ib) 
     &           + Ii * sout%albi(1,ib)
         end do
         
         tran(1) = Id * (sout%ftdd(N2,1) + sout%ftid(N2,1))
     &        + Ii * sout%ftii(N2,1)
         tran(2) = Id * (sout%ftdd(N2,2) + sout%ftid(N2,2))
     &        + Ii * sout%ftii(1,2)
         ! pptr%TRANS_SW = tran(1)
         pptr%TRANS_SW = tran(1) * exp(-1 * vz(1))
         pptr%crad%LAI => ffp
         allocate(T_sun(N2))
         allocate(T_sha(N2))
         allocate(I_sun(N2))
         allocate(I_sha(N2))
         
         do i = 1, N2
            T_sun(i) = IPAR*(sout%ftdd(i,1) + sout%ftid(i,1))
            T_sha(i) = IPAR*sout%ftii(i,1)
            I_sun(i) = IPAR*sout%isun(i,1)
            I_sha(i) = IPAR*sout%isha(i,1)
         end do

      else  !IPAR is less than LOW_LIGHT_LIMIT
         sunlit = 0.d0
         shaded = 0.d0
         N2 = 1
         do ib = 1, N_BANDS
            pptr%albedo(ib) = 0.d0
         end do
         pptr%TRANS_SW = tran(1) * exp(-1 * vz(1))
         pptr%crad%LAI => ffp
         allocate(T_sun(N2))
         allocate(T_sha(N2))
         allocate(I_sun(N2))
         allocate(I_sha(N2))
         do i = 1, N2
            T_sun(i) = 0.d0
            T_sha(i) = 0.d0
            I_sun(i) = 0.d0
            I_sha(i) = 0.d0
         end do
      endif !if IPAR.ge.LOW_LIGHT_LIMIT

#ifdef ENT_STANDALONE_DIAG
      write(1080,*) Id,Ii,pptr%albedo(1),pptr%albedo(2),
     & tran(1),tran(2),I_sun(1),I_sha(1)
      write(1081,*) 'T_sun', T_sun
      write(1081,*) 'T_sha', T_sha
#endif
      
      pptr%crad%f_sun => sunlit
      pptr%crad%f_sha => shaded
      pptr%crad%T_sun => T_sun
      pptr%crad%T_sha => T_sha
      pptr%crad%I_sun => I_sun
      pptr%crad%I_sha => I_sha

      !deallocate(T_sun)
      !deallocate(T_sha)
      !deallocate(I_sun)
      !deallocate(I_sha)
      !deallocate(sunlit)
      !deallocate(shaded)
      !deallocate(ffp)
      !deallocate(rdfp)
      !deallocate(rifp)
      !deallocate(height_levels)
      !deallocate(vz)
      !deallocate(crad_heights)
#ifdef ENT_STANDALONE_DIAG
      write(1081,*) 'crad%T_sun', pptr%crad%T_sun
      write(1081,*) 'crad%T_sha', pptr%crad%T_sha
#endif


      if (Ii .eq. 1.d0) then
         pptr%cellptr%CosZen = tempCos
      end if

      end subroutine get_canopy_rad

      !*********************************************************************
      subroutine GORT_clumping(pptr, height_levels, ffp, rdfp, rifp, 
     &   sunlit, shaded, vz, gin, h_coh)
!@sum Calculate the GORT clumping index in canopy layers and save into
!@sum variable ppt%crad
!@sum Requirements:  pptr cannot be NULL or bare soil -- check before call.
      type(patch), pointer :: pptr
      real*8, dimension(:), pointer :: height_levels
      real*8, dimension(:), pointer :: ffp, rdfp, rifp
      real*8, dimension(:), pointer :: sunlit, shaded, vz 
      !--- Local -------
      type(gort_input) ::  gin1  ! input structure for gort model
      type(gort_input), dimension(:), pointer :: gin
      real*8, dimension(:), pointer :: h_coh
      integer :: num_profiles
      
!      if (.NOT.ASSOCIATED(pptr)) then
!         print *,'In GORT_clumping: pptr=NULL'
!         stop
!      endif

      ! moved to get_canopy_rad()
!      call geo_to_gin(pptr, gin, h_coh)
!      if (.NOT.ASSOCIATED(gin)) then ! bare soil or LAI=0
!          !print *, 'lai zero'
!         call stop_model('GORT_clumping: Bare soil',255)
!      endif 

      num_profiles = size(gin)
      if (num_profiles .eq. 1.d0) then
          gin1%dens_tree = gin(1)%dens_tree
          gin1%h1 = gin(1)%h1
          gin1%h2 = gin(1)%h2
          gin1%horz_radius = gin(1)%horz_radius
          gin1%vert_radius = gin(1)%vert_radius
          gin1%zenith = gin(1)%zenith
          gin1%delta_z = gin(1)%delta_z
          gin1%dens_foliage = gin(1)%dens_foliage
          gin1%pft = gin(1)%pft
          gin1%LAI = gin(1)%LAI
          gin1%dbh = gin(1)%dbh
          ! print *, 'before run_single_gort ...'
          ! print *, "gin(1)%dens_foliage is ", gin(1)%dens_foliage 
          call run_single_gort(gin1, pptr, height_levels, ffp, rdfp, 
     &        rifp, sunlit, shaded, vz, h_coh)
      else
          ! print *, 'before run_conveolute_gort ...'
          call run_convolute_gort(gin, pptr, height_levels, ffp, rdfp, 
     &        rifp, sunlit, shaded, vz, h_coh)
      endif
      
      ! height_levels = pptr%crad%heights
      ! deallocate(h_coh)
      
      end subroutine GORT_clumping

      !*********************************************************************

      subroutine run_single_gort(gin, pptr, crad_heights, fp, rdfp, 
     &    rifp, sunlit, shaded, vz, h_coh)
!      subroutine run_single_gort(gin, direct_light_ratio, ivt, height_levels,
!     &        fp, rdfp, rifp, transmit, sunlit, shaded)

      use ent_pfts, only : pfpar
      ! Input
      type(gort_input) ::  gin  ! input structure for gort model
      type(patch), pointer :: pptr
      ! Output
      real*8, dimension(:), pointer :: crad_heights 
      real*8, dimension(:), pointer :: fp, rdfp, rifp
      real*8, dimension(:), pointer :: sunlit     ! sunlit leaf area fraction
      real*8, dimension(:), pointer :: shaded     ! shaded leaf area fraction
      real*8, dimension(:), pointer :: vz
      real*8, dimension(:), pointer :: h_coh
      !---- Local ---
      real*8, dimension(:), pointer :: efp

      integer :: N_height_level, N2
      integer :: ilevel, jlevel, tmpl
      real*8 :: clumpd, clumpi
      real*8, dimension(:), pointer :: height_levels
      real*8, dimension(:), pointer :: fpt, vzt
      real*8 :: tmp_zenith
      type(cohort),pointer :: cop

      ! ----------------------------------------------------------
      ! DEBUG
      !print *, 'I am in run_single_gort ..........'
   
      cop => pptr%tallest
      K = get_K(cos(gin%zenith),gin%pft)

      ! Check the validity of the inputs
      if (pfpar(gin%pft)%woody) then  !## NEED TO CHECK FOR NON-WOODY, TOO
         call check_inputs(gin)
      endif

      ! define the vertical layers, 
      ! first at interval of 1m, then rescale to input height level 
      call get_height_level(gin, height_levels)
      N_height_level = size(height_levels)
      allocate(fpt(N_height_level))
      allocate(vzt(N_height_level)) 
      !print *, 'before foliage profile' 
      call get_foliage_profile(gin,height_levels,fpt,vzt)
      call layering(pptr%crad,height_levels,fpt,gin%delta_z,h_coh,
     &               N2)
#ifdef ENT_STANDALONE_DIAG
      print *, 'N2=', N2
#endif
      ! N2 = size(pptr%crad%heights)
      !print *, 'before allocat crad'
      allocate(crad_heights(N2))
      !print *, 'after allocat crad'
      crad_heights = pptr%crad%heights
      ! N2 = size(crad_heights)
      !print *, 'after set crad'

      ! get intermediate variables, clump and fp
      allocate(fp(N2))
      allocate(rdfp(N2))
      allocate(rifp(N2))
      allocate(efp(N2))
      allocate(vz(N2))
      !print *, 'after allocat vz'
 
      ! rescale to input height level
      do ilevel = 1, N2
         do jlevel = 1, N_height_level
            if (crad_heights(ilevel) .le. height_levels(jlevel)) then
               tmpl = jlevel
               !print *, 'jlevel=', jlevel
               exit
            end if
         end do
         fp(ilevel) = sum(fpt(tmpl:N_height_level))*gin%delta_z
         vz(ilevel) = vzt(tmpl)
      end do
      !print *, 'N_height_level=', N_height_level
      clumpd = get_analytical_clump(gin)
      tmp_zenith=gin%zenith
      gin%zenith = ang_dif
      clumpi = get_analytical_clump(gin)
      gin%zenith = tmp_zenith

      ! add clump factor into foliage profile data
      efp = K * clumpd * fp
      rdfp = clumpd * fp
      rifp = clumpi * fp

      allocate(sunlit(N2))
      allocate(shaded(N2))
      sunlit = T(efp, gin%zenith)
      shaded = 1 - sunlit
#ifdef ENT_STANDALONE_DIAG
      write(1078,*) "heights",size(crad_heights),crad_heights
      write(1078,*) "ffp", size(fp), fp
      write(1078,*) "K,clump=",K,clumpd,clumpi
      write(1078,*) "rdfp", size(rdfp),rdfp
      write(1078,*) "efp", size(efp),efp
      write(1078,*) "sunlit", size(sunlit),sunlit
#endif
            
      pptr%crad%GORTclump = clumpd
      cop%height => crad_heights
      cop%fp => fp

      end subroutine run_single_gort

      !*********************************************************************

      subroutine run_convolute_gort (gin, pptr, crad_heights, fp, rdfp, 
     &    rifp, sunlit, shaded, vz, h_coh)         ! convoluted profile values
!      subroutine run_convolute_gort (gin, direct_light_ratio, ivt,     &   ! All inputs for gort
!          height_levels, fp, rdfp, rifp, transmit, sunlit, shaded)            ! convoluted profile values

      use ent_pfts, only : pfpar
      ! Input
      type(gort_input), dimension(:), pointer :: gin
      type(patch), pointer :: pptr

      ! Output
      real*8, dimension(:), pointer   :: crad_heights
      real*8, dimension(:), pointer   :: fp
      real*8, dimension(:), pointer   :: rdfp       ! clumping factor * foliage profile
      real*8, dimension(:), pointer   :: rifp       ! clumping factor * foliage profile
      real*8, dimension(:), pointer   :: sunlit     ! sunlit leaf area fraction
      real*8, dimension(:), pointer   :: shaded     ! shaded leaf area fraction
      real*8, dimension(:), pointer   :: vz         ! trunk volumn
      real*8, dimension(:), pointer   :: h_coh

      ! local variables
      real*8, dimension(:), pointer   :: height_levels  ! Convolute profile height levels
      integer :: num_profiles, iprofile, N_height_level, N2, i
      integer :: ilevel, jlevel, tmpl, iszn
      real*8    :: szn, delta_z, dszn
      real*8    :: tmp_zenith
      real*8    :: clumpd, clumpi
      type(profile_params), dimension(:), pointer :: all_convolute_input
      real*8, dimension(:), pointer   :: efp  
      real*8, dimension(:), pointer   :: fpt, rdfpt, rifpt, efpt, vzt
      real*8, dimension(:), pointer   :: tmp_fp, tmp_rdfp, tmp_rifp  
      real*8, dimension(:), pointer   :: tmp_height_levels, tmp_efp
      real*8, dimension(:), pointer   :: tmp_vz

      real*8, dimension(:), pointer   :: d_fp, d_rdfp, d_rifp, d_efp
      real*8, dimension(:), pointer   :: d_height_levels, d_vz 
      real*8, dimension(:,:), pointer :: convoluted_szn_transmit
      type(cohort),pointer :: cop
      !-----------------------------------------------------------------

      ! DEBUG
      ! print *, 'I am in run_convolute_gort.........'
      num_profiles = size(gin)
      K = get_K2(gin)
      ! print *, 'K = ', K


      ! check to verify that all profiles in all_gort_in have the same zenith angle
      szn = gin(1)%zenith
      do iprofile=2, num_profiles
        if (gin(iprofile)%zenith .ne. szn) then
          call stop_model( 'The sun zenith angle of all
     &	     profiles should be the same. Exit.....',255)
        end if
      end do

      !print *, 'The sun zenith angle for the convoluted profile is ', szn


      ! Get effective efp for each profile

      nullify(all_convolute_input)
      allocate(all_convolute_input(num_profiles))

      do iprofile=1, num_profiles

        if (pfpar(gin(iprofile)%pft)%woody) then   !## NEED TO CHECK FOR NON-WOODY, TOO
           call check_inputs(gin(iprofile))
        endif

        ! print *, iprofile, 'iprof after check inputs' 
        ! define the vertical layers
        call get_height_level(gin(iprofile), tmp_height_levels) !tmp_height_levels is allocated in the subroutine
        N_height_level = size(tmp_height_levels)
        !print *, 'N_height_level = ', N_height_level
        !print *, tmp_height_levels

        ! get intermediate variables, clump and fp
        allocate(tmp_fp(N_height_level))
        allocate(tmp_rdfp(N_height_level))
        allocate(tmp_rifp(N_height_level))
        allocate(tmp_efp(N_height_level))
        allocate(tmp_vz(N_height_level))
   
        ! consider putting the following two subroutines as function instead
        call get_foliage_profile(gin(iprofile), tmp_height_levels, 
     &      tmp_fp, tmp_vz) 
        clumpd = get_analytical_clump(gin(iprofile))
        tmp_zenith=gin(iprofile)%zenith
        gin(iprofile)%zenith = ang_dif
        clumpi = get_analytical_clump(gin(iprofile))
        gin(iprofile)%zenith = tmp_zenith

        ! add clump factor into foliage profile data 
        tmp_efp = K * clumpd * tmp_fp 
        tmp_rdfp = clumpd * tmp_fp
        tmp_rifp = clumpi * tmp_fp

        !DEBUG
        !print *, tmp_fp(10), clumpd, tmp_efp(10)
        !print *, minval(tmp_height_levels), maxval(tmp_height_levels)
        !print *, 'Assinging single_gort output to convolute_input'

        ! fill in convolute_input
        all_convolute_input(iprofile)%fp => tmp_fp
        all_convolute_input(iprofile)%rdfp => tmp_rdfp
        all_convolute_input(iprofile)%rifp => tmp_rifp
        all_convolute_input(iprofile)%efp => tmp_efp
        all_convolute_input(iprofile)%vz => tmp_vz
        all_convolute_input(iprofile)%height_levels => tmp_height_levels
     
        !DEBUG
        !print *, all_convolute_input(iprofile)%efp(10)

        ! deallocate and null the tmp pointers
        nullify(tmp_fp)
        nullify(tmp_rdfp)
        nullify(tmp_rifp)
        nullify(tmp_efp)
        nullify(tmp_vz)
        nullify(tmp_height_levels)

      end do
      ! print *, '2 fp: ', all_convolute_input(2)%fp
      ! print *, '1 height_levles: ', all_convolute_input(1)%height_levels
      ! print *, '2 height_levels: ', all_convolute_input(2)%height_levels


      ! Convolute fp, efp by calling convolute, but here only fp is final output for run_gort
      ! efp will be used for getting transmittance for direct light
      ! intermediate results will be recorded in cohort series

      ! print *, 'Calling convolute ......' 
      call convolute(all_convolute_input, pptr,      ! input
     &   fpt, rdfpt, rifpt, efpt, vzt, height_levels)    ! convoluted and run_gort output
      ! print *, 'c height_levels: ', height_levels

      ! print *, 'Calling layering ......' 
      call layering(pptr%crad,height_levels,fpt, gin(1)%delta_z,h_coh,
     &              N2)
      N_height_level = size(height_levels)

      allocate(crad_heights(N2))
      crad_heights = pptr%crad%heights
      ! N2 = size(crad_heights)
      allocate(fp(N2))
      allocate(rdfp(N2))
      allocate(rifp(N2))
      allocate(efp(N2))
      allocate(vz(N2))
      ! print *, 'After allocating fp, rdfp .... N2 =', N2 

      ! rescale to input height level for each cohort
      cop => pptr%tallest
      do while (ASSOCIATED(cop))

      ! print *, 'crad_heights: ', crad_heights
      ! print *, 'height_levels: ', height_levels
      ! print *, 'cop%fp_dz: ', cop%fp_dz
      do ilevel = 1, N2
         do jlevel = 1, N_height_level
            if (crad_heights(ilevel) .le. height_levels(jlevel)) then
              tmpl = jlevel
              exit
            end if
         end do
         fp(ilevel) = sum(cop%fp_dz(tmpl:N_height_level))*gin(1)%delta_z
      end do
      ! print *, 'After associate cohort ......' 
 
        cop%height => crad_heights
        cop%fp => fp
        cop => cop%shorter
      end do
      ! print *, 'After associate cop ......' 

      ! rescale to input height level for whole patch
      do ilevel = 1, N2
         do jlevel = 1, N_height_level
            if (crad_heights(ilevel) .le. height_levels(jlevel)) then
               tmpl = jlevel
               exit
            end if
         end do
         fp(ilevel) = sum(fpt(tmpl:N_height_level))*gin(1)%delta_z
         rdfp(ilevel) = sum(rdfpt(tmpl:N_height_level))*gin(1)%delta_z
         rifp(ilevel) = sum(rifpt(tmpl:N_height_level))*gin(1)%delta_z
         efp(ilevel) = sum(efpt(tmpl:N_height_level))*gin(1)%delta_z
         vz(ilevel) = vzt(tmpl)
      end do
      ! print *, 'After rescale ......' 


      ! Compute the transmit, sunlit, shaded from efp
      ! Sunlit and Shaded will be the run_gort output, but transmit will be depending on 
      ! if transmit from diffuse light needs to be considered

      allocate(sunlit(N2))
      allocate(shaded(N2))
      sunlit = T(efp, szn)
      shaded = 1 - sunlit

      do iprofile=1, num_profiles
         deallocate(all_convolute_input(iprofile)%fp)
         deallocate(all_convolute_input(iprofile)%rdfp)
         deallocate(all_convolute_input(iprofile)%rifp)
         deallocate(all_convolute_input(iprofile)%efp)
         deallocate(all_convolute_input(iprofile)%vz)
         deallocate(all_convolute_input(iprofile)%height_levels)
      end do
      deallocate(all_convolute_input)

      ! pptr%crad%heights = height_levels
      ! pptr%crad%LAI = rdfp
      pptr%crad%GORTclump = clumpd

      end subroutine run_convolute_gort

      !*********************************************************************

!      subroutine TwoStream(lbp, ubp, filter_vegsol, num_vegsol, pcolumn, 
!     &    coszen, Id, height_levels, rdfp, rifp, rho, tau, sout)
      subroutine TwoStream (lbp, ubp, filter_vegsol, num_vegsol, ivt, 
     &   pcolumn, coszen, Id, height_levels, rdfp, rifp, rho, tau, sout)
!
! !DESCRIPTION:
! Two-stream fluxes for canopy radiative transfer
! Use two-stream approximation of Dickinson (1983) Adv Geophysics
! 25:305-353 and Sellers (1985) Int J Remote Sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse
! flux given an underlying surface with known albedo.
!
! !USES:
!    use clmtype
!    use clm_varpar, only : N_BANDS
!    use clm_varcon, only : omegas, tfrz, betads, betais
!
! !ARGUMENTS:
!    implicit none
      use ent_pfts, only : xl
      type(patch),pointer :: pptr
      integer , intent(in)  :: lbp, ubp                 ! pft bounds
      integer , intent(in)  :: filter_vegsol(ubp-lbp+1) ! filter for vegetated pfts with coszen>0
      integer , intent(in)  :: num_vegsol               ! number of vegetated pfts where coszen>0

      integer :: ivt(num_vegsol)
      ! integer , intent(in)  :: ivt(lbp:ubp)             ! pft vegetation type
      integer :: pcolumn(num_vegsol)
      ! integer , intent(in)  :: pcolumn(lbp:ubp)         ! column of corresponding pft

      ! real*8 :: coszen(num_vegsol)
      real*8, intent(in), dimension(:)  :: coszen        ! cosine solar zenith angle for next time step
      real*8, intent(in)  :: Id                       ! direct light ratio, used for isha and isun
      real*8, intent(in), dimension(:)  :: height_levels
      real*8, intent(in), dimension(:)  :: rdfp, rifp
      real*8, intent(in)  :: rho(lbp:ubp,N_BANDS)    ! leaf/stem refl weighted by fraction LAI and SAI
      real*8, intent(in)  :: tau(lbp:ubp,N_BANDS)    ! leaf/stem tran weighted by fraction LAI and SAI

      type(pft_pstate_type) :: sout    ! output structure for 2stream model

! !Constant from CLM code to avoid using 'use' model    
      real*8, parameter :: betads  = 0.5            ! two-stream parameter betad for snow
      real*8, parameter :: betais  = 0.5            ! two-stream parameter betai for snow
      real*8 :: omegas(N_BANDS)           ! two-stream parameter omega for snow by band
      integer :: i
      data (omegas(i),i=1,N_BANDS) /0.8, 0.4, 0.4, 0.4, 0.4, 0.4/
      real*8, parameter :: SHR_CONST_TKFRZ   = 273.15       ! freezing T of fresh water          ~ K 
      real*8, parameter :: tfrz   = SHR_CONST_TKFRZ !freezing temperature [K]

      ! integer, parameter :: N_PFT      =   16   ! number of plant function types
C      real :: xl(N_PFT)         ! ecophys const - leaf/stem orientation index
C      data (xl(i),i=1,N_PFT) /0.01, 0.01, 0.01, 0.10, 0.10, 0.01, 0.25, 
C     &    0.25, 0.01, 0.25, 0.25, -0.30, -0.30, -0.30, -0.30, -0.30/
      real*8 :: albgrd(2,N_BANDS)   ! ground albedo (direct) (column-level)
      real*8 :: albgri(2,N_BANDS)   ! ground albedo (diffuse)(column-level)
      data (albgrd(1,i),i=1,N_BANDS) /0.075, 0.314, 0.314, 0.314, 0.314,
     &    0.314/
      data (albgri(1,i),i=1,N_BANDS) /0.075, 0.314, 0.314, 0.314, 0.314,
     &    0.314/
      data (albgrd(2,i),i=1,N_BANDS) /0.939, 0.787, 0.787, 0.787, 0.787,
     &    0.787/
      data (albgri(2,i),i=1,N_BANDS) /0.939, 0.787, 0.787, 0.787, 0.787,
     &    0.787/
      real*8 :: Ii                 !diffuse light ratio
!
! !CALLED FROM:
! subroutine SurfaceAlbedo in this module
!
! !REVISION HISTORY:
! Author: Gordon Bonan
! Modified for speedup: Mariana Vertenstein, 8/26/02
! Vectorized routine: Mariana Vertenstein:  8/20/03
!
! !LOCAL VARIABLES:
!
! local pointers to implicit in scalars
!
!    integer , pointer :: pcolumn(:)    ! column of corresponding pft
!    real, pointer :: albgrd(:,:)   ! ground albedo (direct) (column-level)
!    real, pointer :: albgri(:,:)   ! ground albedo (diffuse)(column-level)
!    real, pointer :: t_veg(:)      ! vegetation temperature (Kelvin)
!    real, pointer :: fwet(:)       ! fraction of canopy that is wet (0 to 1)
!    integer , pointer :: ivt(:)        ! pft vegetation type
!    real, pointer :: xl(:)         ! ecophys const - leaf/stem orientation index
!
! local pointers to implicit out scalars
!
      real*8, pointer :: albd(:,:)     ! surface albedo (direct)
      real*8, pointer :: albi(:,:)     ! surface albedo (diffuse)
      real*8, pointer :: fabd(:,:)     ! flux absorbed by veg per unit direct flux
      real*8, pointer :: fabi(:,:)     ! flux absorbed by veg per unit diffuse flux
      real*8, pointer :: ftdd(:,:)     ! down direct flux below veg per unit dir flx
      real*8, pointer :: ftid(:,:)     ! down diffuse flux below veg per unit dir flx
      real*8, pointer :: ftii(:,:)     ! down diffuse flux below veg per unit dif flx
      real*8, pointer :: frid(:,:)     ! up diffuse flux below veg per unit dir flx
      real*8, pointer :: frii(:,:)     ! up diffuse flux below veg per unit dif flx
      real*8, pointer :: gdir(:)       ! leaf projection in solar direction (0 to 1)
      real*8, pointer :: omega(:,:)    ! fraction of intercepted radiation that is scattered (0 to 1)
      real*8, pointer :: isha(:,:)     ! shaded flux per unit incident flux                     
      real*8, pointer :: isun(:,:)     ! sunlit flux per unit incident flux                     
!
!EOP
!
! !OTHER LOCAL VARIABLES:
!
      integer  :: fp, p, c         ! array indices
      !integer  :: ic               ! 0=unit incoming direct; 1=unit incoming diffuse
      integer  :: ib !, i          ! waveband number !Have to move declaration of i to above.
      real*8 :: cosz             ! 0.001 <= coszen <= 1.000
      real*8 :: asu              ! single scattering albedo
      real*8 :: chil(lbp:ubp)    ! -0.4 <= xl <= 0.6
      real*8 :: twostext(lbp:ubp)! optical depth of direct beam per unit leaf area
      real*8 :: avmu(lbp:ubp)    ! average diffuse optical depth
      real*8 :: omegal           ! omega for leaves
      real*8 :: betai            ! upscatter parameter for diffuse radiation
      real*8 :: betail           ! betai for leaves
      real*8 :: betad            ! upscatter parameter for direct beam radiation
      real*8 :: betadl           ! betad for leaves
      real*8 :: tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9 ! temporary
      real*8 :: p1,p2,p3,p4,s1,s2,u1,u2,u3                        ! temporary
      real*8 :: b,c1,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10   ! temporary
      real*8 :: phi1,phi2,sigma                                   ! temporary
      real*8 :: temp0(lbp:ubp),temp1,temp2(lbp:ubp)               ! temporary
      real*8 :: t1
      real*8, pointer :: evai(:)
      real*8 :: delta_z, vai(lbp:ubp)

      integer :: N_height_level, ilevel
!      integer , intent(in)  :: ivt(lbp:ubp)             ! pft vegetation type

!-----------------------------------------------------------------------
      N_height_level = size(height_levels)
      Ii = 1 - Id
#ifdef ENT_STANDALONE_DIAG
      print *, "N_height_level, Ii=", N_height_level, Ii
#endif

!    if (ubp-lbp > 0) then 
!       allocate(vai(lbp:ubp))
!    end if
      allocate(evai(N_height_level))
      allocate(sout%albd(lbp:ubp,N_BANDS))
      allocate(sout%albi(lbp:ubp,N_BANDS))
      allocate(sout%fabd(N_height_level,N_BANDS))
      allocate(sout%fabi(N_height_level,N_BANDS))
      allocate(sout%ftdd(N_height_level,N_BANDS))
      allocate(sout%ftid(N_height_level,N_BANDS))
      allocate(sout%ftii(N_height_level,N_BANDS))
      allocate(sout%frid(N_height_level,N_BANDS))
      allocate(sout%frii(N_height_level,N_BANDS))
      allocate(sout%gdir(lbp:ubp))
      allocate(sout%omega(lbp:ubp,N_BANDS))
      allocate(sout%isha(N_height_level,N_BANDS))
      allocate(sout%isun(N_height_level,N_BANDS))

    ! Assign local pointers to derived subtypes components (column-level)

!    albgrd  => clm3%g%l%c%cps%albgrd
!    albgri  => clm3%g%l%c%cps%albgri

    ! Assign local pointers to derived subtypes components (pft-level)

!      pcolumn => clm3%g%l%c%p%column
!      fwet    => sout%fwet
!      t_veg   => clm3%g%l%c%p%pes%t_veg
!      ivt     => clm3%g%l%c%p%itype
      albd    => sout%albd
      albd    => sout%albd
      albi    => sout%albi
      fabd    => sout%fabd
      fabi    => sout%fabi
      ftdd    => sout%ftdd
      ftid    => sout%ftid
      ftii    => sout%ftii
      frid    => sout%frid
      frii    => sout%frii
      gdir    => sout%gdir
      omega   => sout%omega
      isha    => sout%isha
      isun    => sout%isun
!      xl      => pftcon%xl

      ! Calculate two-stream parameters omega, betad, betai, avmu, gdir, twostext.
      ! Omega, betad, betai are adjusted for snow. Values for omega*betad
      ! and omega*betai are calculated and then divided by the new omega
      ! because the product omega*betai, omega*betad is used in solution.
      ! Also, the transmittances and reflectances (tau, rho) are linear
      ! weights of leaf and stem values.

      do fp = 1,num_vegsol
         p = filter_vegsol(fp)
       
         ! note that the following limit only acts on cosz values > 0 and less than 
         ! 0.001, not on values cosz = 0, since these zero have already been filtered
         ! out in filter_vegsol
         cosz = max(0.001, coszen(p))
       
         chil(p) = min( max(xl(ivt(p)), -0.4), 0.6 )
         if (abs(chil(p)) <= 0.01) chil(p) = 0.01
         phi1 = 0.5 - 0.633*chil(p) - 0.330*chil(p)*chil(p)
         phi2 = 0.877 * (1.-2.*phi1)
         ! gdir(p) = phi1 + phi2*cosz
         ! twostext(p) = gdir(p)/cosz
         gdir(p) = K
         twostext(p) = gdir(p)/cosz
	 
         !write(1076,*) 'size(cosz),p,cosz,K=',size(coszen),p,cosz,K
         avmu(p) = (1. - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
         temp0(p) = gdir(p) + phi2*cosz
         temp1 = phi1*cosz
         temp2(p) = (1. - temp1/temp0(p) * log((temp1+temp0(p))/temp1))
      end do

      do ib = 1, N_BANDS
         do fp = 1,num_vegsol
            p = filter_vegsol(fp)
            c = pcolumn(p)

            omegal = rho(p,ib) + tau(p,ib)
            asu = 0.5*omegal*K/temp0(p) *temp2(p)
            betadl = (1.+avmu(p)*twostext(p))/(omegal*avmu(p)*
     &         twostext(p))*asu
            betail = 0.5 *((rho(p,ib)+tau(p,ib)) + (rho(p,ib)-tau(p,ib)) 
     &         * ((1.+chil(p))/2.)**2) / omegal

            ! Adjust omega, betad, and betai for intercepted snow

            !if (t_veg(p) > tfrz) then                             !no snow
               tmp0 = omegal
               tmp1 = betadl
               tmp2 = betail
            !else
            !   tmp0 =   (1.-fwet(p))*omegal        + fwet(p)*omegas(ib)
            !   tmp1 = ( (1.-fwet(p))*omegal*betadl + fwet(p)*omegas(ib)*betads ) / tmp0
            !   tmp2 = ( (1.-fwet(p))*omegal*betail + fwet(p)*omegas(ib)*betais ) / tmp0
            !end if
            omega(p,ib) = tmp0           
            betad = tmp1 
            betai = tmp2  

            ! Absorbed, reflected, transmitted fluxes per unit incoming radiation

            b = 1. - omega(p,ib) + omega(p,ib)*betai
            c1 = omega(p,ib)*betai
            tmp0 = avmu(p)*twostext(p)
            d = tmp0 * omega(p,ib)*betad
            f = tmp0 * omega(p,ib)*(1.-betad)
            tmp1 = b*b - c1*c1
            h = sqrt(tmp1) / avmu(p)
            sigma = tmp0*tmp0 - tmp1
            p1 = b + avmu(p)*h
            p2 = b - avmu(p)*h
            p3 = b + tmp0
            p4 = b - tmp0
          
            ! Determine fluxes for vegetated pft for unit incoming direct 
            ! Loop over incoming direct and incoming diffuse
            ! 0=unit incoming direct; 1=unit incoming diffuse

            ! ic = 0 unit incoming direct flux
            ! ========================================

            u1 = b - c1/albgrd(c,ib)
            u2 = b - c1*albgrd(c,ib)
            u3 = f + c1*albgrd(c,ib)

            tmp2 = u1 - avmu(p)*h
            tmp3 = u1 + avmu(p)*h
            tmp4 = u2 + avmu(p)*h
            tmp5 = u2 - avmu(p)*h
            h1 = -d*p4 - c1*f
            tmp6 = d - h1*p3/sigma
            h4 = -f*p3 - c1*d
            tmp8 = h4/sigma

            ! delta_z = height_levels(2)-height_levels(1)
            vai(p) = rdfp(N_height_level)
            ! vai(p) = rdfp(1)
            ! PET, 3/1/04: added this test to avoid floating point errors in exp()
            t1 = min(h*vai(p), 40.d0)
            s1 = exp(-t1)
            t1 = min(twostext(p)*vai(p), 40.d0)
            s2 = exp(-t1)
          
            d1 = p1*tmp2/s1 - p2*tmp3*s1
            d2 = tmp4/s1 - tmp5*s1          
            tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
            h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
            h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
            tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
            h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
            h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
            h7 = (c1*tmp2) / (d1*s1)
            h8 = (-c1*tmp3*s1) / d1
            h9 = tmp4 / (d2*s1)
            h10 = (-tmp5*s1) / d2


            ! Flux reflected by vegetation (ic = 0)
            albd(p,ib) = h1/sigma + h2 + h3

            do ilevel=1, N_height_level
               evai(ilevel) = rdfp(ilevel)   ! cumulative clumping factor * foliage profile

               t1 = min(h*evai(ilevel), 40.d0)
               s1 = exp(-t1)
               t1 = min(twostext(p)*evai(ilevel), 40.d0)
               s2 = exp(-t1)

               ! Downward direct and diffuse fluxes below vegetation (ic = 0)

               ftdd(ilevel,ib) = s2
               ftid(ilevel,ib) = h4*s2/sigma + h5*s1 + h6/s1
               frid(ilevel,ib) = h1*s2/sigma + h2*s1 + h3/s1
               ! Flux absorbed by vegetation (ic = 0)

               fabd(ilevel,ib) = 1. - albd(p,ib) 
     &             - (1.-albgrd(c,ib))*ftdd(ilevel,ib) 
     &             - (1.-albgri(c,ib))*ftid(ilevel,ib)
          
            end do
          

            ! ic = 1 unit incoming diffuse
            ! ========================================

            u1 = b - c1/albgri(c,ib)
            u2 = b - c1*albgri(c,ib)
            u3 = f + c1*albgri(c,ib)

            tmp2 = u1 - avmu(p)*h
            tmp3 = u1 + avmu(p)*h
            tmp4 = u2 + avmu(p)*h
            tmp5 = u2 - avmu(p)*h
            h1 = -d*p4 - c1*f
            tmp6 = d - h1*p3/sigma
            h4 = -f*p3 - c1*d
            tmp8 = h4/sigma

            vai(p) = rifp(N_height_level)
            ! vai(p) = rifp(1)
            t1 = min(h*vai(p), 40.d0)
            s1 = exp(-t1)
            t1 = min(twostext(p)*vai(p), 40.d0)
            s2 = exp(-t1)

            d1 = p1*tmp2/s1 - p2*tmp3*s1
            d2 = tmp4/s1 - tmp5*s1
            tmp7 = ( d - c1 - h1/sigma*(u1+tmp0) ) * s2
            h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
            h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
            tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
            h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
            h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
            h7 = (c1*tmp2) / (d1*s1)
            h8 = (-c1*tmp3*s1) / d1
            h9 = tmp4 / (d2*s1)
            h10 = (-tmp5*s1) / d2

            ! Flux reflected by vegetation
  
            albi(p,ib) = h7 + h8

            do ilevel=1, N_height_level
               evai(ilevel) = rifp(ilevel)   ! cumulative clumping factor * foliage profile

               t1 = min(h*evai(ilevel), 40.d0)
               s1 = exp(-t1)
               t1 = min(twostext(p)*evai(ilevel), 40.d0)
               s2 = exp(-t1)

               ! Downward direct and diffuse fluxes below vegetation

               ftii(ilevel,ib) = h9*s1 + h10/s1
               frii(ilevel,ib) = h7*s1 + h8/s1

               ! Flux absorbed by vegetation

               fabi(ilevel,ib) = 1. - albi(p,ib) 
     &            - (1.-albgri(c,ib))*ftii(ilevel,ib)
	     
	       ! Sunlit/shaded intensity
	     
               isha(ilevel,ib) = (1. - omegal) * (Id * (ftid(ilevel,ib) 
     &            + frid(ilevel,ib)) + Ii * (ftii(ilevel,ib) 
     &            + frii(ilevel,ib)))
               isun(ilevel,ib) = (1. - omegal) * twostext(p) * Id 
     &            + isha(ilevel,ib)
             
            end do   ! end of foliage level loop
         end do   ! end of pft loop
      end do   ! end of radiation band loop
    
      end subroutine TwoStream

! ----------------------------------------------------------

      function T(fp, szn)

      ! return a transmittance array based on the size of input fp

      real*8, dimension(:), pointer :: fp
      real*8 :: szn
      real*8, dimension(size(fp)) :: T
   
      ! local vars
      integer :: N_height_level, ilevel

      N_height_level = size(fp)
      do ilevel=1, N_height_level
        T(ilevel) = exp(-fp(ilevel) / cos(szn))
      end do

      end function 
! ----------------------------------------------------------
      subroutine check_inputs(gin)
      implicit none
      ! check the validity of inputs
      type(gort_input) :: gin                  ! gort input

      real*8 :: dens_tree            ! Crown count density (m-2)
      real*8 :: dens_foliage         ! foliage area volume density of a single crown(m-1)
      real*8 :: h1, h2               ! lower and upper bound of crown centers (m)
      real*8 :: delta_z              ! height increase for each height level
      real*8 :: horz_radius          ! horizontal crown radius (m)
      real*8 :: vert_radius          ! vertical crown radius (m)
      real*8 :: zenith               ! solar zenith angle (degree)

      dens_tree = gin%dens_tree
      dens_foliage = gin%dens_foliage
      h1 = gin%h1
      h2 = gin%h2
      delta_z = gin%delta_z
      horz_radius = gin%horz_radius
      vert_radius = gin%vert_radius
      zenith = gin%zenith
#ifdef ENT_STANDALONE_DIAG
      print *,"dens_tree,dens_fol,h1,h2,dz,hrad,vertrad,zenith=",
     &     dens_tree,dens_foliage,h1,h2,delta_z,horz_radius,vert_radius,
     &     zenith
#endif      
      ! check density
      !if ( (dens_tree <= 0.) .or. (dens_foliage <= 0.)) then
      if ( (dens_tree <= 0.) ) then
        !write(*,*) 'Density must be greater than 0. STOPPING ...' !WZ
         call stop_model('Bare soil, check before canrad',255) !NK
         !This check should have been done before calling canopy radiation.-NK
      end if

      ! check the height
      if ((h1 <= 0.) .or. (h2 <= 0.)) then
        call stop_model('Center of crown height must be >0.',255)
      end if

      if (h1 > h2) then
        call stop_model('h2 must be no less than h1.',255)
      end if

      if ( h1 < vert_radius ) then
        write(*,*) 'canopyradiation.f: dens_tree, dens_foliage,
     & h1, h2, delta_z, horz_radius, vert_radius, zenith',
     &    dens_tree, dens_foliage, h1, h2, delta_z, 
     &    horz_radius, vert_radius, zenith
       call stop_model('h1 must be greater than vert_radius.',255)
      end if

      ! check radius
      if ((horz_radius <= 0.) .or. (vert_radius <= 0.)) then
        call stop_model('Radius must be greater than 0.',255)
      end if

      ! check angle
      if ((zenith < 0.) .or. (zenith > PI/2.)) then
      !if ((zenith.ge.0.5*PI).and.(zenith.le.1.5*PI)) !Night
        !call stop_model('Solar zenith angle must be within [0.,1.]. CONTINUING ...'
      end if

      end subroutine check_inputs

! ----------------------------------------------------------
      subroutine get_height_level(gin, height_levels)

        ! compute the height levels between h1 and h2
        ! NOTE : here last layer could be less than delta_z.
        !        It may be necessary to return the dz at the last layer to
        !        ensure the accuracy of calculation, otherwise, we can make
        !        sure (h2-h1) can be integer divided by delta_z.

      type(gort_input), intent(in) :: gin
      real*8, dimension(:), pointer :: height_levels

      ! local
      real*8 :: h1, h2, vert_radius, delta_z
      integer :: n_levels, i, temp1
      real*8 :: jh1, jh2, temp2, dif

      ! -------------------------------------------


      h1 = gin%h1
      h2 = gin%h2
      vert_radius = gin%vert_radius
      delta_z = gin%delta_z

      jh1 = h1 - vert_radius
      jh2 = h2 + vert_radius

      ! Added the rounding checking for different machines.

      temp1 = ceiling((jh2-jh1)/delta_z)
      temp2 = (jh2-jh1)/delta_z

      if ((temp1-temp2) > 0.999) then
        n_levels = ceiling((jh2-jh1)/delta_z)
      else
        n_levels = ceiling((jh2-jh1)/delta_z) + 1
      end if


      allocate(height_levels(n_levels))

      do i =1, n_levels-1
        height_levels(i) = jh1 + delta_z * (i-1)
      end do

      height_levels(n_levels) = jh2

      end subroutine get_height_level

! ----------------------------------------------------------

      subroutine get_foliage_profile(gin, height_levels, fp, vz)

      use ent_pfts, only : pfpar
      type(gort_input) :: gin
      real*8, dimension(:), intent(in) :: height_levels
      real*8, dimension(:), intent(inout) :: fp, vz

      ! local
      real*8 :: h1, h2, vert_radius, tmp, tmpv
      real*8 :: z, iz, jh11, jh12, jh21, jh22, hdif
      real*8 :: maxv
      integer :: n_levels, i

      real*8, dimension(:), allocatable :: xx

      ! ---------------------------------
      ! print *, 'entering get_foliage_profile...'
      n_levels = size(height_levels)
      allocate(xx(n_levels))

      if (pfpar(gin%pft)%woody) then
      h1 = gin%h1
      h2 = gin%h2
      vert_radius = gin%vert_radius

      tmp=gin%dens_tree * gin%dens_foliage * PI * (gin%horz_radius ** 2)
#ifdef ENT_STANDALONE_DIAG
      write(1077,*) "LAI, dbh, lumda, fa, R, b:", gin%LAI, gin%dbh,
     &  gin%dens_tree, gin%dens_foliage,
     &  gin%horz_radius, gin%vert_radius
#endif
      jh11 = h1 - vert_radius
      jh12 = h1 + vert_radius
      jh21 = h2 - vert_radius
      jh22 = h2 + vert_radius
      hdif = h2 - h1

      ! allocate(vz(n_levels))
      ! allocate(fp(n_levels))
      ! print *, 'n_levels=', n_levels
      ! define xx based on the location of each height level
      do i = 1, n_levels

        z = height_levels(i)
        if (hdif < gin%delta_z) then
           if ((z >= jh21) .and. (z <= jh22))
     &         xx(i) = 1.0 - ((z-h2) / vert_radius) ** 2
        else ! (h1 != h2) 
        if (jh12 < jh21) then
           if ((z >= jh11) .and. (z <= jh12))   
     &         xx(i) = height_function1(z, h1, vert_radius, hdif)
           if ((z > jh12) .and. (z <= jh21))   
     &         xx(i) = height_function3(hdif, vert_radius)
           if ((z > jh21) .and. (z <= jh22))   
     &         xx(i) = height_function2(z, h2, vert_radius, hdif)
        else ! (jh12 >= jh21)
           if ((z >= jh11) .and. (z <= jh21))   
     &         xx(i) = height_function1(z, h1, vert_radius, hdif)
           if ((z > jh21) .and. (z <= jh12))   
     &         xx(i) = height_function4(z, h1, h2, vert_radius)
           if ((z > jh12) .and. (z <= jh22))   
     &         xx(i) = height_function2(z, h2, vert_radius, hdif)
        end if
        end if

        fp(i) = tmp * xx(i)
        
        xx(i) = 0.0

        if (hdif < gin%delta_z) then
          maxv = h2 - z
          if (maxv < 0) maxv = 0
          xx(i) = gin%dbh * gin%dens_tree * tan(gin%zenith) * maxv
        else
          tmpv = gin%dbh * gin%dens_tree * tan(gin%zenith) / hdif
          iz =gin%h1 + gin%delta_z/2 
          do while (iz <= gin%h2) 
          ! do iz = gin%h1 + gin%delta_z/2, gin%h2, gin%delta_z
             maxv = iz - z
             if (maxv < 0) maxv = 0
             xx(i) = xx(i) + gin%delta_z * maxv
             ! vz(i) = vz(i) + gin%delta_z * (z2/2-min(z,iz))
             iz = iz + gin%delta_z 
          end do
          xx(i) = xx(i) * tmpv
        end if
        ! print *, 'xx(i)=', xx(i) 
        ! vz(i) = xx(i)
      end do

      vz = xx

      else  ! herbacious
        do i = 1, n_levels
           fp(i) = gin%LAI / (n_levels * gin%delta_z)
           vz(i) = 0
        end do 
      endif

      ! deallocate(xx)

      end subroutine get_foliage_profile

! ----------------------------------------------------
 
      subroutine convolute(all_convolute_input, pptr,     ! input
     &   fp, rdfp, rifp, efp, vz, height_levelC)         ! convoluted output


      ! NOTE: This subroutine convolute multiple profiles into one.
      !       The algorithm used here is mapping each profile to the desired height levels
      !       and then convolute (either sum or multiply) profiles together.
      !       The height_levelC will be based upon the input for desired height, however,
      !       it will be adjusted to cover all the participating profiles.
      !
      ! note the following params are going to be decided by the input profiles
      !  real :: hc_low       ! Low end for the convoluted profile
      !  real :: hc_high      ! High end for the convoluted profile
      !  real :: delta_z      ! Delta z for the convoluted profile
      ! --------------------------------------------------------------------------


      ! Input
      type(profile_params), dimension(:), pointer :: all_convolute_input
      type(patch),pointer :: pptr


      ! Output
      real*8, dimension(:), pointer   :: fp         ! foliage profile
      real*8, dimension(:), pointer   :: rdfp       ! clumping factor * foliage profile
      real*8, dimension(:), pointer   :: rifp       ! clumping factor * foliage profile
      real*8, dimension(:), pointer   :: efp        ! effective foliage profile
      real*8, dimension(:), pointer   :: vz         ! stem volumn
      real*8, dimension(:), pointer   :: height_levelC  ! Convolute profile height levels


      ! Local variables
      real*8, dimension(:), allocatable :: arrthl, arrthh, arrdz
      real*8 :: hc_low, hc_high, delta_z
      integer :: n_profiles, iprofile, n_convolute_levels
      type(profile_params), pointer :: interpolate_out
      type(cohort),pointer :: cop

      !---------------------------------------------------------------------

      ! print *, 'I am in convolute now....'

      n_profiles = size(all_convolute_input)

      ! print *, all_convolute_input(1)%fp(10)
      ! print *, all_convolute_input(2)%fp(10)


      ! Get the height levels for the convoluted profile
      !----------------------------------------------------------------------------------------------
      allocate(arrthl(n_profiles))
      allocate(arrthh(n_profiles))
      allocate(arrdz(n_profiles))

      ! compare with all profiles and adjust the low end and high end

      do iprofile = 1, n_profiles
        arrthl(iprofile) = minval(
     &      all_convolute_input(iprofile)%height_levels)
        arrthh(iprofile) = maxval(
     &      all_convolute_input(iprofile)%height_levels)
        arrdz(iprofile) = all_convolute_input(iprofile)%height_levels(2) 
     &      - all_convolute_input(iprofile)%height_levels(1)
      end do

      ! print *, arrthl
      ! print *, arrthh
      ! print *, arrdz

      hc_low=minval(arrthl)
      hc_high=maxval(arrthh)
      delta_z=maxval(arrdz)
      ! print *, 'hc_low, hc_high, delta_z = ', hc_low, hc_high, delta_z

      deallocate(arrthl, arrthh, arrdz)

      ! print *, 'Getting convolute height levels ...'
      call get_convolute_height_level(hc_low, hc_high, delta_z,
     &   height_levelC)
      n_convolute_levels=size(height_levelC)


      ! allocate output variables
      ! print *, 'Allocating for the convoluted outputs ....'

      allocate(fp(n_convolute_levels))
      allocate(rdfp(n_convolute_levels))
      allocate(rifp(n_convolute_levels))
      allocate(efp(n_convolute_levels))
      allocate(vz(n_convolute_levels))

      ! allocate temporary interpolate_out
      ! print *, 'Allocating for the interpolation output ....'
      allocate(interpolate_out)
      allocate(interpolate_out%fp(n_convolute_levels))
      allocate(interpolate_out%rdfp(n_convolute_levels))
      allocate(interpolate_out%rifp(n_convolute_levels))
      allocate(interpolate_out%efp(n_convolute_levels))
      allocate(interpolate_out%vz(n_convolute_levels))
      allocate(interpolate_out%height_levels(n_convolute_levels))

      interpolate_out%height_levels = height_levelC


      ! interpolate profiles to the convoluted height level
      !---------------------------------------------------------

      cop => pptr%tallest
      ! interpolate profile1 to the convoluted height levels
      do iprofile = 1, n_profiles


         ! print *, 'Interpolating for profile No. ', iprofile
         call interpolate_profile (all_convolute_input(iprofile),   
     &      interpolate_out)
         ! print *, interpolate_out%fp(10)

         ! print *, 'Reassign or convolute the output'

         if (iprofile == 1) then
            fp = interpolate_out%fp
            rdfp = interpolate_out%rdfp
            rifp = interpolate_out%rifp
            efp = interpolate_out%efp
            vz = interpolate_out%vz
         else
            fp = fp + interpolate_out%fp
            rdfp = rdfp + interpolate_out%rdfp
            rifp = rifp + interpolate_out%rifp
            efp = efp + interpolate_out%efp
            vz = vz + interpolate_out%vz
         end if

         cop%height_dz => height_levelC
         cop%fp_dz => interpolate_out%fp
         cop => cop%shorter 
C      print *, iprofile, 'convolute input fp ', 
C     &        all_convolute_input(iprofile)%fp
C      print *, iprofile, 'convolute output fp ', 
C     &        interpolate_out%fp 

      end do

      !deallocate(interpolate_out%fp)
      !deallocate(interpolate_out%rdfp)
      !deallocate(interpolate_out%rifp)
      !deallocate(interpolate_out%efp)
      !deallocate(interpolate_out%vz)
      !deallocate(interpolate_out%height_levels)
      !deallocate(interpolate_out)


      return


      end subroutine convolute


 ! ----------------------------------------------------------------------------------

      subroutine interpolate_profile (interpolate_in, interpolate_out)

      ! Input
      type(profile_params), intent(in) :: interpolate_in          ! data for input profile

      ! Output
      type(profile_params), pointer :: interpolate_out            ! interpolated profile

      ! local variables
      integer :: i, j
      integer :: nlevel_in, nlevel_out
      real*8 :: h, hin_low, hin_high
      integer :: id1, id2
      integer, dimension(1) :: id
      real*8, dimension(:), pointer :: height_in, height_out      ! height levels for input, output profile
      real*8, dimension(:), pointer :: diff

      ! ---------------------------------------
      ! print *, interpolate_in%fp(10)

      height_in => interpolate_in%height_levels
      height_out => interpolate_out%height_levels

      nlevel_in = size(height_in)
      nlevel_out = size (height_out)
      hin_low = height_in(1)
      hin_high = height_in(nlevel_in)
      allocate(diff(nlevel_in))


      ! loop through each height level for the convoluted height levels.
      ! Based upon the height comparison to get the interpolated profile

      do i=1, nlevel_out

        h = height_out(i)

        if ( (h >= hin_low) .and. (h <= hin_high)) then

          ! find the closest two height level
          do j=1, nlevel_in
            diff(j) = (height_in(j) - h)
          end do
          id = minloc ( abs(diff) ) ! find the location that has the height closes to h
          if (diff(id(1)) > 0) then
            id1 = id(1) - 1
            id2 = id(1)
          else
            id1 = id(1)
            id2 = id(1) + 1
          end if

          ! use the data from the above two points to interpolate for the desired height
       

          interpolate_out%fp(i) = interpolate_in%fp(id1) + (h - 
     &       height_in(id1)) * (interpolate_in%fp(id2) - 
     &       interpolate_in%fp(id1)) / (height_in(id2) - height_in(id1))
          interpolate_out%rdfp(i) = interpolate_in%rdfp(id1) + (h - 
     &       height_in(id1)) * (interpolate_in%rdfp(id2) -  
     &       interpolate_in%rdfp(id1)) / (height_in(id2) - 
     &       height_in(id1))
          interpolate_out%rifp(i) = interpolate_in%rifp(id1) + (h - 
     &       height_in(id1)) * (interpolate_in%rifp(id2) -
     &       interpolate_in%rifp(id1)) / (height_in(id2) - 
     &       height_in(id1))
          interpolate_out%efp(i) = interpolate_in%efp(id1) + (h - 
     &       height_in(id1)) * (interpolate_in%efp(id2) -   
     &       interpolate_in%efp(id1)) /(height_in(id2) - height_in(id1))
          interpolate_out%vz(i) = interpolate_in%vz(id1) + (h - 
     &       height_in(id1)) * (interpolate_in%vz(id2) -   
     &       interpolate_in%vz(id1)) /(height_in(id2) - height_in(id1))
       

        else

          ! just fill in the output profile with 0

          if (h < hin_low) then
            interpolate_out%fp(i) = 0
            interpolate_out%rdfp(i) = 0
            interpolate_out%rifp(i) = 0
            interpolate_out%efp(i) = 0
            interpolate_out%vz(i) = 0
          end if

          if (h > hin_high) then
            interpolate_out%fp(i) = 0
            interpolate_out%rdfp(i) = 0
            interpolate_out%rifp(i) = 0
            interpolate_out%efp(i) = 0
            interpolate_out%vz(i) = 0
          end if

          ! just fill in the output profile with the ceiling or bottom value

          !if (h <= hin_low) then
          !  interpolate_out%fp(i) = interpolate_in%fp(1)
          !  interpolate_out%rdfp(i) = interpolate_in%rdfp(1)
          !  interpolate_out%rifp(i) = interpolate_in%rifp(1)
          !  interpolate_out%efp(i) = interpolate_in%efp(1)
          !  interpolate_out%vz(i) = interpolate_in%vz(1)
          !end if

          !if (h >= hin_high) then
          !  interpolate_out%fp(i) = interpolate_in%fp(nlevel_in)
          !  interpolate_out%rdfp(i) = interpolate_in%rdfp(nlevel_in)
          !  interpolate_out%rifp(i) = interpolate_in%rifp(nlevel_in)
          !  interpolate_out%efp(i) = interpolate_in%efp(nlevel_in)
          !  interpolate_out%vz(i) = interpolate_in%vz(nlevel_in)
          !end if

        end if

      end do

      nullify(height_in)
      nullify(height_out)

      deallocate(diff)
      return

      end subroutine interpolate_profile


 !-------------------------------------------------------------------------------

      subroutine get_convolute_height_level(h1, h2, delta_z, 
     &   height_levels)

      ! compute the height levels between h1 and h2
      ! NOTE : here last layer could be less than delta_z.
      !        It may be necessary to return the dz at the last layer to
      !        ensure the accuracy of calculation, otherwise, we can make
      !        sure (h2-h1) can be integer divided by delta_z.

      real*8, intent(in) :: h1, h2, delta_z
      real*8, dimension(:), pointer :: height_levels

      ! local
      integer :: n_levels, i, temp1
      real*8 :: temp2, dif

      ! -------------------------------------------

      ! Added the rounding checking for different machines.

      temp1 = ceiling((h2-h1)/delta_z)
      temp2 = (h2-h1)/delta_z

      if ((temp1-temp2) > 0.999) then
        n_levels = ceiling((h2-h1)/delta_z)
      else
        n_levels = ceiling((h2-h1)/delta_z) + 1
      end if


      allocate(height_levels(n_levels))

      do i =1, n_levels-1
        height_levels(i) = h1 + delta_z * (i-1)
      end do

      height_levels(n_levels) = h2
    
      end subroutine get_convolute_height_level

! ----------------------------------------------------

      subroutine geo_to_gin(pp, gin, h_coh)
!@sum Calculate gin from Ent input, get mean patch geometry from cohort value,
!@sum then get the gort input at patch level, note the whole Entcell share the 
!@sum same input radiation.
!@sum Note in calculating clumping indices, dbh is not needed
!@sum Required inputs:  Check for bare soil before calling this routine.
      
      use ent_pfts, only : pfpar
      type(patch),pointer :: pp
      type(gort_input),dimension(:),pointer :: gin
      real*8,pointer :: h_coh(:)
      !----Local----------------!
      type(cohort),pointer :: cop
      integer :: count, i
      real*8 :: vc
 
#ifdef ENT_STANDALONE_DIAG
      print *,'Started geo_to_gin' ! with patch
#endif

!      if ( .NOT.ASSOCIATED(pp%tallest)) then ! bare soil
!        pp%crad%heights = 0.d0
!        pp%crad%LAI = 0.d0
!        pp%crad%GORTclump = 0.d0
!     return
!      endif

      if (( pp%tallest%pft.eq.0).or.(pp%tallest%pft > N_PFT)) then
        print *, 'GORT_clumping: wrong pft = ', pp%tallest%pft
        call patch_print(6,pp,"ERROR ")
        call stop_model("GORT_clumping: wrong pft",255)
      endif

      !* LOOP THROUGH COHORTS *!
      count = 0
      cop => pp%tallest
      do while (ASSOCIATED(cop))
        !* Assign vegpar
         if (cop%n.gt.0.d0) then
         !if (cop%LAI.gt.lai_thres) then
         ! if (cop%LAI.gt.0.d0) then
            count = count + 1
         endif

         cop => cop%shorter
      end do

    !  if (count == 0) then
    !     print *, 'no vegetation,L,b=', cop%LAI,
    ! & cop%crown_dy
    !     pp%crad%heights = 0.d0
    !     pp%crad%LAI = 0.d0
    !     pp%crad%GORTclump = 0.d0
    !     return
    !  endif
 
      if (count.gt.0.d0) then
        if (count.eq.1.d0) then ! single cohort
          allocate(gin(1))
          allocate(h_coh(1))   
          cop => pp%tallest
          gin(1)%dens_tree = cop%n
          gin(1)%vert_radius = cop%crown_dy  ! debugging
          gin(1)%h1 = cop%h - gin(1)%vert_radius
          gin(1)%h2 = gin(1)%h1
          gin(1)%horz_radius = cop%crown_dx
          !print *, 'r,b=',cop%crown_dx, cop%crown_dy
          gin(1)%zenith = acos(pp%cellptr%CosZen)
          gin(1)%delta_z = dz_gin
!          vc = 1.33333 * PI * gin(1)%horz_radius ** 2 * 
!     &      gin(1)%vert_radius
          if (pfpar(cop%pft)%woody) then
          vc = (4.d0/3.d0) * PI * (gin(1)%horz_radius ** 2) * 
     &      gin(1)%vert_radius
          gin(1)%dens_foliage = cop%LAI / (vc * 
     &      gin(1)%dens_tree)
          else
            gin(1)%dens_foliage = cop%LAI
          endif
          gin(1)%pft = cop%pft
          gin(1)%LAI = cop%LAI
          gin(1)%dbh = cop%dbh
          h_coh(1) = cop%h 
        else ! multi cohort
          allocate(gin(count)) 
          allocate(h_coh(count))
  
          !* LOOP THROUGH COHORTS *!
          cop => pp%tallest
          i = 0;
          do while (ASSOCIATED(cop))
            if (cop%n.gt.0.d0) then
            ! if (cop%LAI.gt.lai_thres) then
            ! if (cop%LAI.gt.0.d0) then
              i = i + 1
              gin(i)%dens_tree = cop%n
              gin(i)%h1 = cop%h - cop%crown_dy
              gin(i)%h2 = gin(i)%h1
              gin(i)%horz_radius = cop%crown_dx
              gin(i)%vert_radius = cop%crown_dy
              gin(i)%zenith = acos(pp%cellptr%CosZen)
              gin(i)%delta_z = dz_gin

              if (pfpar(cop%pft)%woody) then
              vc = 1.33333 * PI * gin(i)%horz_radius ** 2 * 
     &           gin(i)%vert_radius
              gin(i)%dens_foliage = cop%LAI / (vc * gin(i)%dens_tree)
              else
                 gin(i)%dens_foliage = cop%LAI
              endif
              gin(i)%pft = cop%pft
              gin(i)%LAI = cop%LAI
              gin(i)%dbh = cop%dbh
              h_coh(i) = cop%h 
            endif 
	
            cop => cop%shorter
          end do
        endif
      else
         !* This should never be called, as check for bare soil should
         !* have been done outside of this module.
        !print *,'no vegetation,L,b='!, cop%LAI, cop%crown_dy
         print *,'GORT: No vegetation, bare soil'
        !call patch_print(6,pp,"ERROR ")
        !call stop_model("GORT_clumping: no vegetation",255)
      endif
 
      end subroutine geo_to_gin

! ----------------------------------------------------
      function get_analytical_clump(gin)
!@sum calculate clumping index from canopy geometry 

      use ent_pfts, only : pfpar
      type(gort_input),intent(in) :: gin
      real*8 :: get_analytical_clump

      ! local variables:
      real*8 :: tr, a, b, c, d, G

      if (.not. pfpar(gin%pft)%woody .or. 
     &      gin%dens_foliage.lt.lai_thres) then
         get_analytical_clump = 1.d0 !No foliage clumping, no foliage -NK
      else
         a = (tan(gin%zenith)) ** 2
         b = (gin%vert_radius/gin%horz_radius) ** 2
         c = sqrt((1+a)/(1+a*b))
         G = get_K(cos(gin%zenith),gin%pft)
         tr = G * gin%dens_foliage * gin%horz_radius * c

         d = (1 - (2*tr+1) * exp(-2*tr)) / (2*tr*tr)
         
         get_analytical_clump = 3*(1-d)/(4*tr)
      endif
      end function get_analytical_clump

! ------------------------------------------------------

      function height_function1(z, h, b, hdif)
!@sum the four height functions describes 3D GORT geometry, 
!@sum from Ni-Meister et al., 2001

      real*8,intent(in) :: z, h, b, hdif
      real*8 :: height_function1

      height_function1 = ((b+z-h)**2) * (2*b + h - z) / (3*hdif*b**2)

      end function height_function1

! ------------------------------------------------------

      function height_function2(z, h, b, hdif)

      real*8,intent(in) :: z, h, b, hdif
      real*8 :: height_function2

      height_function2 = ((b+h-z)**2) * (2*b + z - h) / (3*hdif*b**2)

      end function height_function2

! -------------------------------------------------------

      function height_function3(hdif, b)

      real*8,intent(in) :: b, hdif
      real*8 :: height_function3

      height_function3 = (4./3.) * b / hdif

      end function height_function3

! -------------------------------------------------------

      function height_function4(z, h1, h2, b)

      real*8,intent(in) :: z, h1, h2, b
      real*8 :: height_function4

      height_function4 = 1 - ((h2-h1)**2 + 3*(z-h1)*(z-h2)) / (3*b**2)

      end function height_function4

! --------------------------------------------------------

      function get_K(CosZen, pft)
!@sum G is the relative projected area of leaf elements in the direction zenith
!@sum get_K here is actually G, following the definition of Li et al., 1995
!@sum G function comes from Sellers, 1985 

      use ent_pfts, only : xl
      real*8,intent(in) :: CosZen
      integer,intent(in) :: pft
      real*8 :: get_K

      ! integer, parameter :: N_PFT = 16   ! number of plant function types
      integer :: i 
C      real*8 :: xl(N_PFT)       ! ecophys const - leaf/stem orientation index
C      data (xl(i),i=1,N_PFT) /0.01, 0.01, 0.01, 0.10, 0.10, 0.01, 0.25, 
C     &    0.25, 0.01, 0.25, 0.25, -0.30, -0.30, -0.30, -0.30, -0.30/
      real*8 :: phi1,phi2,chil   
                                    
      chil = min( max(xl(pft), -0.4d0), 0.6d0 )
      if (abs(chil) <= 0.01) chil = 0.01
      phi1 = 0.5 - 0.633*chil - 0.330*chil*chil
      phi2 = 0.877 * (1.-2.*phi1)
      get_K = phi1 + phi2*CosZen

      end function get_K

! --------------------------------------------------------

      function get_K2(gin)
!@sum G is the relative projected area of leaf elements in the direction zenith
!@sum get_K2 here is actually G for single or multiple layers

      type(gort_input), dimension(:), pointer :: gin  ! input structure for gort model
      real*8 :: get_K2
      real*8 :: szn

      ! local variables
      integer :: num_profiles, iprofile
      real*8, dimension(:), pointer :: clump, G
      real*8 :: eLAI, rLAI
      
      num_profiles = size(gin)
      allocate(clump(num_profiles))
      allocate(G(num_profiles))
      eLAI = 0
      rLAI = 0
      
      do iprofile = 1, num_profiles
         szn = gin(iprofile)%zenith
C         print *, 'dens_foliage, pft = ', gin(iprofile)%dens_foliage, 
C     &         gin(iprofile)%pft
         G(iprofile) = get_K(cos(szn),gin(iprofile)%pft)
         clump(iprofile) = get_analytical_clump(gin(iprofile))
         !print *, 'G, clump =', G(iprofile), clump(iprofile)
         eLAI = eLAI + clump(iprofile) * G(iprofile) * gin(iprofile)%LAI
         rLAI = rLAI + clump(iprofile) * gin(iprofile)%LAI
      end do

      !print *, 'eLAI, rLAI = ', eLAI, rLAI
      if (rLAI .eq. 0.d0) then
         get_K2 = 0.5
      else
         get_K2 = eLAI / rLAI
      endif

      deallocate(clump)
      deallocate(G)

      end function get_K2

      subroutine layering(crad, height_levels, fpt, dz, h_coh,
     &                    n2)
!@sum Calculate the needed height levels according to LAI profile and 
!@sum cohort heights, and save to variable ppt%crad

      use qsort_c_module

      !type(patch), pointer :: pptr
      type(canradtype) :: crad
      real*8, dimension(:), pointer :: height_levels, fpt
      real*8 :: dz
      real*8, dimension(:), pointer :: h_tmp, ha 
      real*8, dimension(:), pointer :: dlai, clai, tdlai 
      real*8, dimension(:), pointer :: h_coh, lai_bound
      integer :: len_h, len_l, i, j, len, n2
      integer, dimension(:), pointer :: index

      ! need to put height and foliage profile upside down

      len = size(height_levels)
      allocate(h_tmp(len))
      allocate(dlai(len))
      allocate(clai(len))
      allocate(tdlai(len))

      do i=1,len ! need to reverse heights, so h(1)=top
         h_tmp(len+1-i)=height_levels(i)
         dlai(len+1-i)=fpt(i)
      end do

      ! cumulative lai from top of canopy
      clai(1)=dlai(1)
      do i=2,len
         clai(i)=clai(i-1)+dlai(i)
      end do
      clai = clai * dz

      !write(1077,*) "h_coh in layering:", h_coh
      ! write(1077,*) "pptr%crad%h_coh", pptr%crad%h_coh
      len_h = ceiling(clai(len))+3             ! add 0, 0.5, 1.5
      len_l = size(h_coh)-1
      !len_l = size(pptr%crad%h_coh)-1          ! don't count the tallest cohort

      allocate(lai_bound(len_h))
      allocate(index(len_h)) 
      allocate(ha(len_h+len_l))
      n2 = len_h+len_l

      do i=1,len_h
         if (i .le. 4) then
            lai_bound(i) = (i-1)*0.5
         else
            lai_bound(i) = i - 3
         end if 
         tdlai = clai - lai_bound(i)
         do j=1,len
            tdlai(j) = abs(tdlai(j))
         end do
         index(i) = minloc(tdlai,dim=1)
         ha(i) = h_tmp(index(i))
      end do
 
      crad%h_lai => ha
   
      if (len_l.gt.0) then  ! multiple cohorts
      do i=len_h+1, len_h+len_l
         ha(i)=h_coh(i-len_h+1)
         ! ha(i)=pptr%crad%h_coh(i-len_h+1)
      end do
      endif

      call QsortC(ha)

      !write(1077,*) "ha before", ha
      !ha(1) = height_levels(len)
      !write(1077,*) "ha after", ha

      crad%heights => ha

      deallocate(h_tmp)
      deallocate(dlai)
      deallocate(clai)
      deallocate(tdlai)
      deallocate(lai_bound)
      deallocate(index)
      !deallocate(ha)
      !write(1077,*) "crad%heights", crad%heights

      end subroutine layering

      end module canopyrad
