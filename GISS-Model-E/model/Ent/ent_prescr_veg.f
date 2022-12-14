      module ent_prescr_veg
!@sum ent_prescr_veg - This module contains vegetation structure routines 
!@+   that are primarily specific to the Matthews (1983) prescribed vegetation
!@+   arrays.
!@+   Uses same allometry functions as ent_prescribed_drv_geo, but
!@+   with Matthews arrays for alamax, height.
!@+   Contains soil color initialization by veg type.
      
      use ent_const
!      use ent_pfts

      implicit none
      private
      save

      public 
     &     init_params,!prescr_calcconst renamed init_params 
     &     prescr_veg_albedo,prescr_calc_rootprof,
     &     prescr_calc_lai
!      public prescr_plant_cpools !Moved to allometryfn.f
      public prescr_calc_hdata
     &     ,prescr_calc_woodydiameter,prescr_get_pop,prescr_get_crownrad
     &     ,prescr_calc_initnm, prescr_calc_rootprof_all
     &     ,prescr_calc_soilcolor
      public popdensity

#ifdef ENT_STANDALONE_DIAG
      public print_ent_pfts
#endif

      real*8,parameter :: EDPERY=365. !GISS CONST.f

      contains

!***************************************************************************

      subroutine init_params()
!@sum Initialize array constants for Ent, mostly for soil biogeochemistry.
      use ent_const
      use ent_pfts
      !--Local------
      integer :: n
      real*8 :: lnscl

      !* annK - Turnover time of litter and soil carbon *!
      do n = 1, N_PFT
        if(pfpar(n)%lrage.gt.0.0d0)then
          annK(n,LEAF)  = (1.0d0/(pfpar(n)%lrage*secpy))
          annK(n,FROOT) = (1.0d0/(pfpar(n)%lrage*secpy))
        else
          annK(n,LEAF)  = 1.0d-40  !CASA originally 1.0e-40
          annK(n,FROOT) = 1.0d-40
        end if

        if(pfpar(n)%woodage.gt.0.0d0)then
          annK(n,WOOD)  = 1.0d0/(pfpar(n)%woodage*secpy)
        else
          annK(n,WOOD)  = 1.0d-40
        end if
       !* iyf: 1/(turnover times) for dead pools.  Want annK in sec-1.
        annK(n,SURFMET)    = 14.8d0  /secpy
        annK(n,SURFMIC)    = 6.0d0   /secpy
        annK(n,SURFSTR)    = 3.9d0   /secpy
        annK(n,SOILMET)    = 18.5d0  /secpy
        annK(n,SOILMIC)    = 7.3d0   /secpy
        annK(n,SOILSTR)    = 4.9d0   /secpy        ! 4.8 in casa v3.0
        annK(n,CWD)        = 0.2424d0/secpy
        annK(n,SLOW)       = 0.2d0   /secpy
        annK(n,PASSIVE)    = 0.1d0 * 0.02d0  /secpy
      enddo

#ifdef DEBUG
      do n = 1,N_PFT
        write(98,*) 'pft',n
        write(98,*) 'annK(n,LEAF)',annK(n,LEAF)
        write(98,*) 'annK(n,FROOT)',annK(n,FROOT)
        write(98,*) 'annK(n,WOOD)' ,annK(n,WOOD) 
        write(98,*) 'annK(n,SURFMET)',annK(n,SURFMET)
        write(98,*) 'annK(n,SURFMIC)',annK(n,SURFMIC)
        write(98,*) 'annK(n,SURFSTR)',annK(n,SURFSTR)
        write(98,*) 'annK(n,SOILMET)',annK(n,SOILMET)
        write(98,*) 'annK(n,SOILMIC)',annK(n,SOILMIC)
        write(98,*) 'annK(n,SOILSTR)',annK(n,SOILSTR)
        write(98,*) 'annK(n,CWD)',annK(n,CWD)    
        write(98,*) 'annK(n,SLOW)',annK(n,SLOW)   
        write(98,*) 'annK(n,PASSIVE)',annK(n,PASSIVE)
      enddo
#endif

      !* solubfrac - Soluble fraction of litter*!
      !structurallignin, lignineffect - frac of structural C from lignin, effect of lignin on decomp -PK 6/29/06 
      do n = 1,N_PFT
        lnscl = pfpar(n)%lit_C2N * pfpar(n)%lignin * 2.22 !lignin:nitrogen scalar        
        solubfract(n) = 0.85 - (0.018 * lnscl)
        structurallignin(n) = (pfpar(n)%lignin * 0.65 * 2.22) 
     &                      / (1. - solubfract(n))
        lignineffect(n) = exp(-3.0 * structuralLignin(n))
      end do

      end subroutine init_params

!**************************************************************************
      
      real*8 function prescr_calc_lai(pnum,jday,hemi ) RESULT(lai)
!@sum Returns GISS GCM leaf area index for given vegetation type, julian day
!@+   and hemisphere
      use ent_const
      use ent_pfts, only: alamax, alamin, laday
      !real*8, intent(out) :: lai !@var lai leaf area index - returned
      integer, intent(in) :: pnum !@var pnum cover type
      integer, intent(in) :: jday !@var jday julian day
      integer, intent(in) :: hemi !@var hemi =1 in N. hemisphere, =-1 S.hemi
      !-----Local variables------
      real*8 dphi

      dphi = 0
      if ( hemi < 0 ) dphi = 2d0*pi*.5d0

      !* Return lai *!
      lai =  .5d0 * (alamax(pnum) + alamin(pnum))
     $     + .5d0 * (alamax(pnum) - alamin(pnum))
     $     * cos( 2d0*pi*(laday(pnum)-jday)/dble(EDPERY) + dphi )

      end function prescr_calc_lai


!*************************************************************************

      subroutine prescr_veg_albedo(hemi, ncov, jday, albedo)
!@sum Returns prescribed (Matthews, 1983) albedo for vegetation of type pft 
      !use ent_pfts, only:  COVEROFFSET, albvnd
      use ent_pfts, only:  ALBVND
      integer, intent(in) :: hemi !@hemi hemisphere (-1 south, +1 north)
      !integer, intent(in) :: pft !@var pftlike iv, plant functional type
      integer, intent(in) :: ncov !@var cover type, soil or pft, if pft then ncov=pft+COVEROFFSET
      integer, intent(in) :: jday !@jday julian day
      real*8, intent(out) :: albedo(N_BANDS) !@albedo returned albedo
      !----------Local----------
      integer, parameter :: NV=N_COVERTYPES
      !@var SEASON julian day for start of season (used for veg albedo calc)
C                      1       2       3       4
C                    WINTER  SPRING  SUMMER  AUTUMN
      real*8, parameter, dimension(4)::
     *     SEASON=(/ 15.00,  105.0,  196.0,  288.0/)
C**** parameters used for vegetation albedo
!@var albvnd veg alb by veg type, season and band
!@+   albvnd has been moved to ent_pfts_ENT.f - NK

ccc or pass k-vegetation type, L-band and 1 or 2 for Hemisphere
      integer k,kh1,kh2,l
      real*8 seasn1,seasn2,wt2,wt1
c
c                      define seasonal albedo dependence
c                      ---------------------------------
c
      seasn1=-77.0d0
      do k=1,4
        seasn2=SEASON(k)
        if(jday.le.seasn2) go to 120
        seasn1=seasn2
      end do
      k=1
      seasn2=380.0d0
  120 continue
      wt2=(jday-seasn1)/(seasn2-seasn1)
      wt1=1.d0-wt2
      if ( hemi == -1 ) then    ! southern hemisphere
        kh1=1+mod(k,4)
        kh2=1+mod(k+1,4)
      else                      ! northern hemisphere
        kh1=1+mod(k+2,4)
        kh2=k
      endif

      do l=1,6

!        albedo(l)=wt1*ALBVND(pft+COVEROFFSET,kh1,l)
!     &        +wt2*ALBVND(pft+COVEROFFSET,kh2,l)
        albedo(l)=wt1*ALBVND(ncov,kh1,l)
     &        +wt2*ALBVND(ncov,kh2,l)
      enddo

      end subroutine prescr_veg_albedo

!**************************************************************************

      subroutine prescr_calc_rootprof(rootprof, ncov)
!@sum Return prescribed array rootprof of fractions of roots in soil layer
!@+   for single cover type. (Rosenzweig & Abrampoulos, 1997)
!@+   !Cohort/patch level.
      !This could be moved to allometryfn, but is kept here since it
      !is the original Rosenzweig & Abrampoulos (1997) formulation.
      use ent_pfts, only: COVEROFFSET, aroot, broot
      real*8 :: rootprof(:)
      integer :: ncov !plant functional type + COVEROFFSET
      !-----Local variables------------------
      real*8,parameter :: dz_soil(1:6)=  !N_DEPTH
     &     (/  0.99999964d-01,  0.17254400d+00,
     &     0.29771447d+00,  0.51368874d+00,  0.88633960d+00,
     &     0.15293264d+01 /)
      integer :: n,l
      real*8 :: z, frup,frdn

c**** calculate root fraction afr averaged over vegetation types
      !Initialize zero
      do l=1,N_DEPTH
        rootprof(l) = 0.0
      end do
      do n=1,N_DEPTH
        if (dz_soil(n) <= 0.0) exit !Get last layer w/roots in it.
      end do
      n=n-1
      z=0.
      frup=0.
      do l=1,n
        z=z+dz_soil(l)
        frdn=aroot(ncov)*z**broot(ncov) !cumulative root distrib.
        !frdn=min(frdn,one)
        frdn=min(frdn,1d0)
        if(l.eq.n)frdn=1.
        rootprof(l) = frdn-frup
        frup=frdn
      end do
      !Return rootprof(:)
      end subroutine prescr_calc_rootprof

!**************************************************************************

      subroutine prescr_calc_rootprof_all(rootprofdata)
!@sum Return array of prescribed root fraction profiles by soil depth 
!@+   for all cover types.
      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH) 
      !---Local--------
      integer :: ncov !plant functional type + COVEROFFSET     

      do ncov=1,N_COVERTYPES
        call prescr_calc_rootprof(rootprofdata(ncov,:), ncov)
        !Return array rootprof of fractions of roots in soil layer
        !by vegetation type.
      end do
      end subroutine prescr_calc_rootprof_all

!**************************************************************************

      subroutine prescr_calc_hdata(hdata)
!@sum Return array of prescribed (Matthews, 1983) vegetation heights (m)
!@+   by vegetation type.
      use ent_pfts, only : vhght
      real*8 :: hdata(N_COVERTYPES) 
      !------

      ! For prescr Model E replication, don't need to fill in an
      ! i,j array of vegetation height, but just can use
      ! constant arry for each vegetation pft type.
      ! For full-fledged Ent, will need to read in a file
      ! containing (i,j,N_PFT) matrix of vegetation heights.

      !* Return hdata heights for all vegetation types
      hdata = vhght
      end subroutine prescr_calc_hdata

!**************************************************************************

      subroutine prescr_calc_initnm(nmdata)
!@sum Return PFT array of mean canopy nitrogen (nmv; g/m2[leaf])
      use ent_pfts, only : nmv
      real*8 :: nmdata(N_COVERTYPES)
      !-------

      !* Return intial nm for all vegetation and cover types
      nmdata = nmv
      end subroutine prescr_calc_initnm

!*************************************************************************

      subroutine prescr_get_pop(dbhdata,laimaxdata,popdata)
!@sum Return PFT array of vegetation population density (#/m2) calculated
!@+   from Matthews (1983) veg type heights and alamax.
      !* Derived from Moorcroft, et al. (2001)
      use ent_pfts, only : COVEROFFSET, alamax
      real*8,intent(in) :: dbhdata(N_COVERTYPES)
      real*8,intent(in) :: laimaxdata(N_COVERTYPES)
      real*8,intent(out) :: popdata(N_COVERTYPES)
      !---Local-----------
      integer :: ncov,pft

      popdata(:) = 0.0 !Zero initialize, and zero bare soil.
      do pft=1,N_PFT
        ncov = pft + COVEROFFSET
        !popdata(ncov) = popdensity(pft,dbhdata(ncov),alamax(ncov))
        popdata(ncov) = popdensity(pft,dbhdata(ncov),laimaxdata(ncov))
      enddo
      end subroutine prescr_get_pop

!*************************************************************************
      real*8 function popdensity(pft,dbh,LAImax) Result(popdens)
!@sum (#plants/m^2) Return plant population density for a PFT based on
!@+   allometric functions from Moorcroft et al. (2001).
!@+   
      use ent_pfts, only: pfpar, COVEROFFSET, alamax
      use allometryfn, only : wooddensity_gcm3
       integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8, intent(in) :: LAImax
      !---Local-----------
      real*8 :: Blmax !Max foliage mass per tree (kg-C/individ)
      real*8 :: wooddens !Wood density (g/cm^3)

      if (.not.pfpar(pft)%woody) then
        popdens = 10.d0       !Grass ##HACK See Stampfli et al 2008 (~25 seedlings/m2 for cover %1-10, but big range)
      else if (dbh > 0.0001) then
        wooddens = wooddensity_gcm3(pft)
        Blmax = 0.0419d0 * (dbh**1.56d0) * (wooddens**0.55d0)
        !popdens = (alamax(pft+COVEROFFSET)/pfpar(pft)%sla)/Blmax 
        popdens = (LAImax/pfpar(pft)%sla)/Blmax 
        !leaf area/ground area/(leaf area/kg-C) / (kg-C/individ)
      else
        popdens = 0.d0
      endif
      end function popdensity

!*************************************************************************

      subroutine prescr_calc_woodydiameter(hdata, wddata)
      !* Return array of woody plant diameters at breast height (dbh, cm)
      use ent_pfts, only : COVEROFFSET
      use allometryfn, only : height2dbh
      real*8,intent(in) :: hdata(N_COVERTYPES)
      real*8,intent(out) :: wddata(N_COVERTYPES)
      !----Local---------
      integer :: ncov,pft

      wddata(:) = 0.0           !Zero initialize.
      do pft = 1,N_PFT
         ncov = pft + COVEROFFSET
         wddata(ncov) = height2dbh(pft,hdata(ncov))
      enddo
      end subroutine prescr_calc_woodydiameter

!*************************************************************************

      subroutine prescr_get_crownrad(popdata,craddata)
!@sum prescr_get_crownrad - assumes closed-canopy packing of crowns
!@+      in rows and columns (not staggered).
      use ent_pfts, only : COVEROFFSET
      use allometryfn, only : Crown_rad_max_from_density
      real*8,intent(in) :: popdata(N_COVERTYPES)
      real*8,intent(out) :: craddata(N_COVERTYPES)
      !---Local----
      integer :: n, pft

      craddata(:) = 0.0 !Zero initialize.
      do pft=1,N_PFT
        n = pft + COVEROFFSET
        !craddata(n) = crown_radius_closed(popdata(n))
        craddata(n) = Crown_rad_max_from_density(popdata(n))
      end do

      end subroutine prescr_get_crownrad

!*************************************************************************

!      subroutine prescr_init_Clab(pft,n,laimax,Clabile)
!@sum prescr_init_Clab - Initializes labile carbon pool
!@sum Deciduous woody: Clab = 4 x max Cfol of plant.
!@sum Evergreen woody: Clab = 0.5 x max Cfol of plant.
!@sum 4x requirement is from Bill Parton (personal communication).
!@sum 7/1/2012 - Replaced with different subroutine - NK

!! This subroutine has been revised and moved to phenology.f.      
!      subroutine prescr_init_Clab_old(pft,n,cpool)
!!@sum prescr_init_Clab - Initializes labile carbon pool to 4x mass of 
!!@sum (alamax - alamin) for woody and perennial plants. 
!!@sum 4x requirement is from Bill Parton (personal communication).
!!@sum For herbaceous annuals, assume seed provides 0.5 of alamax mass (guess).
!!@sum 7/1/2012 - Replaced with different subroutine - NK
!      use ent_pfts, only : COVEROFFSET, pfpar !, alamax, alamin
!      implicit none
!      integer, intent(in) :: pft
!      real*8, intent(in) :: n !Density (#/m^2)
!      real*8, intent(in) :: laimax, laimin
!      real*8, intent(inout) :: cpool(N_BPOOLS) !g-C/pool/plant
!      
!!      cpool(LABILE) = 0.5d0*alamax(pft+COVEROFFSET)/pfpar(pft)%sla/n*1d3 !g-C/individ.
!
!      if (pfpar(pft)%phenotype.ne.ANNUAL) then
!        !Enough to grow peak foliage and fine roots.
!        cpool(LABILE) = (alamax(pft+COVEROFFSET)-alamin(pft+COVEROFFSET)
!     &       )*4.d0/pfpar(pft)%sla/n*1d3 !g-C/individ.
!      else
!       if( n > 0.d0 ) then
!        cpool(LABILE) = 0.5d0*alamax(pft+COVEROFFSET)/
!     &       pfpar(pft)%sla/n*1d3 !g-C/individ.
!       else
!        cpool(LABILE) = 0.d0
!       endif
!      endif
!
!      end subroutine prescr_init_Clab_old
!*************************************************************************


      subroutine prescr_calc_soilcolor(soil_color)
!@sum Return arrays of Matthews (1983) soil color and texture.
      !## Can get rid of this subroutine and replace with array assignment.
      use ent_pfts, only : soil_color_prescribed
      integer, intent(out) :: soil_color(N_COVERTYPES)
      !------

      soil_color(:) = soil_color_prescribed(:)

      end subroutine prescr_calc_soilcolor
!*************************************************************************
#ifdef ENT_STANDALONE_DIAG
      !Assume PS_MODEL=FBB if running Ent_standalone, so can print out FBBpfts.f.
      subroutine print_ent_pfts()
!@sum Print all parameter sets used on initialization of Ent.
      use ent_const
      use ent_types
      use ent_pfts
      use FarquharBBpspar
      integer pft

      write(*,*) "ent_pfts pfpar:"
      write(*,*) "pst,woody,leaftype,hwilt,sstar,swilt,nf,sla,
     &r,lrage,woodage,lit_C2N,lignin,croot_ratio,phenotype 
     &b1Cf, b2Cf, b1Cd, b2Cd, b1Ht, b2Ht"
       do pft = 1,N_PFT
         write(*,*) pft, pfpar(pft)
      enddo

      write(*,*) "FarquharBBpspar pftpar: "
      write(*,*) "pst, PARabsorb,Vcmax,m,b,Nleaf"
      do pft = 1,N_PFT
         write(*,*) pft, pftpar(pft)
      enddo
      end subroutine print_ent_pfts
#endif
!*************************************************************************
      end module ent_prescr_veg

