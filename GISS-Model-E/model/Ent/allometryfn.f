      module allometryfn
!@sum Physics-level functions to calculate plant physiological state, 
!@+   allometry (geometry and biomass), and tissue properties, shc, 
!@+   and initial C_labile.
!@+   Individual level, NO knowledge of entcells, patches, or cohort struct.
!@+   Used for initialization of vegetation structure, updating phenology and
!@+   growth, and calculating biophysical state.
!@auth NYK  July 2012

      use ent_const
      use ent_pfts

      implicit none

      public do_geo
      public GRASS, HERB, SHRUB, TREE, BARE
      public GISS_shc
      public wooddensity_gcm3
      public dbh2Cfol,dbh2height, height2dbh
      public Cfol2height_herb !, height2Cfol_herb
      public dbh2Cdead, Cdead2dbh, Cdead2Chw, Chw2Ccroot
      public dDBHdCdead, dDBHdCfol, dHdDBH
      public maxdbh, Cfol_fn, Csw_fn
      public Crown_rad_allom, Crown_rad_max_from_density
      public crown_radius_horiz_allom, crown_radius_vert
      public allom_plant_cpools, init_Clab, nplant

      logical, save :: do_geo = .false.!Can change on initialization

      !*q: ratio of root to leaf biomass (unitless) (value from ED)
      real*8, parameter :: q=1.0d0 

      !*hw_fract: ratio of above ground stem to total stem (stem plus structural roots) (value from ED)
      real*8, parameter :: hw_fract = 0.70d0 
      real*8, parameter :: B2C = 0.5d0  !Inverse, to avoid having to divide.

      contains
!*************************************************************************

      subroutine init_do_geo_flag
!@sum Flag to do calculations for geographic version instead of Matthews.
!@auth N.Y.Kiang
!@+   Created this for modules that call GISS_shc and that need to
!@+     calculate mean annual LAI, which is currently really clunky.

      do_geo = .true.
      end subroutine init_do_geo_flag


      integer function gform(PFT)
!@sum Plant growth form
!@auth N.Y.Kiang
      integer :: PFT

      if ((PFT.ge.1).and.(PFT.le.8)) then
         gform = TREE
      elseif ((PFT.eq.9).or.(PFT.eq.10)) then
         gform = SHRUB
      elseif ((PFT.ge.11).and.(PFT.le.14)) then
         gform = GRASS
      elseif ((PFT.eq.15)) then
         gform = HERB
      elseif (PFT.eq.16) then
         gform = TREE
      else
         gform = BARE
      endif
      end function gform

!*************************************************************************

      real*8 function GISS_shc(meanLAI) Result(shc)
!@sum Returns GISS GCM specific heat capacity for cohort as function of LAI.
!     meanLAI = (sum over pfts) 
!           {.5d0*(alamax(anum) + alamin(anum)) * vfraction }/sum(vfraction)
!     I.e. GISS GCM computes shc_entcell = shc(mean entcell LAI*vfaction)
!     instead of shc_entcell = mean(shc(patch LAI)*vfraction)
      use ent_const
      real*8,intent(in) :: meanlai    

      !Seems like this ought to use actual LAI, too? - NK
      !-----Local----

      shc = (.010d0+.002d0*meanLAI+.001d0*meanLAI**2)*shw*rhow

      end function GISS_shc

!*************************************************************************

      real*8 function wooddensity_gcm3(pft) Result(wooddens)
!@sum Wood density (g-Biomass cm-3). Moorcroft et al. (2001), corrected
!@+   for long-lived leaves (evergreen needleleaf and tundra) by averages 
!@+   from the Wood Density Database (WDD)
!@+   http://www.worldagroforestry.org/sea/products/afdbases/wd/
!@+   Returns wood density in total mass per volume 
!@+   (total mass = dry mass = C + N + everything else)
!@+   The Moorcroft et al. (2001) formulation of wood density as a function
!@+   of leaf lifespan (lrage used as proxy) predicts too high density for
!@+   long-lived leaves like evergreen needeleaf and tundra shrubs.
!@+   This formulation probably is only valid for hardwood and not ring-porous.
!@+   Cold shrub density of 1.05 is temporary guess, for lack of data.
!@+   This function might be replaced by array wdens_g_cm3.
      use ent_pfts, only : pfpar
      integer,intent(in) :: pft

      if (pfpar(pft)%leaftype.eq.NEEDLELEAF) then
        wooddens = 0.5d0
      else
        wooddens = min(
     &       max(0.5d0, 0.5d0 + 0.2d0*(pfpar(pft)%lrage-1.d0)), 1.05d0)       
      endif

      end function wooddensity_gcm3

!*************************************************************************
      real*8 function sla(pft,llspan)
!@sum Specific leaf area (m^2/kgC) as function of leaf lifespan
      !Used only for tropical evergreen broadleaf trees.
      integer, intent(in) :: pft
      real*8, intent(in) :: llspan !Leaf lifespan (years)
      
      if (pfpar(pft)%phenotype .eq. EVERGREEN .and. 
     &    pfpar(pft)%leaftype .eq. BROADLEAF .and.
     &    llspan .gt. 0.d0) then 
         sla = 10.0**(1.6923-0.3305*log10(llspan))
      else 
         sla = pfpar(pft)%sla
      endif
      end function sla
!*************************************************************************

      real*8 function maxdbh(pft)
!@sum Maximum dbh (cm) for woody plant of given PFT.
      integer,intent(in) :: pft       

      maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+a0h(pft))-a0h(pft))  
     &     /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht

c$$$      if (.not.pfpar(pft)%woody) then !grasses/crops 
c$$$         maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+1.3)-1.3)  
c$$$     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
c$$$      else !woody   
c$$$         maxdbh=log(1.0-(0.999*pfpar(pft)%b1Ht-1.3)  
c$$$     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
c$$$      end if
 
      end function maxdbh
!*************************************************************************

      real*8 function dbh2Cdead(pft,dbh)
!@sum (gC/plant) Carbon in dead of structural tissue as a function of dbh
!@+   Structural tissue = dead stem (HW) + coarse roots (CR)
!@+   Follows formulation for cold deciduous trees from Albani et al. (2006).
!@+   Other PFTs for Ent full documentation.
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh !(cm)
  
!      !Old original ED formulation (Moorcroft et al., 2001).
!       Bs = .069d0*(h**0.572d0)*(dbh**1.94d0) * 
!     &       (wooddensity_gcm3(pft)**0.931d0) *1d3

      !Albani et al. (2006) formulation
      dbh2Cdead = B2C * pfpar(pft)%b1Cd * 
     &     (dbh**pfpar(pft)%b2Cd) * 1000.0d0 !kgC to gC

      end function  dbh2Cdead

!*************************************************************************
      real*8 function Cdead2dbh(pft,Cdead)
!@sum (cm) Diameter at breast height (DBH) given carbon in structural
!@+   tissue.  Inverse of dbh2Cdead.
      integer,intent(in) :: pft
      real*8, intent(in) :: Cdead !gC/plant

!      Cdead2dbh = (Cdead*B2C/1000.0d0/pfpar(pft)%b1Cd)
!     &            **(1.0d0/pfpar(pft)%b2Cd)

      Cdead2dbh = (Cdead/(B2C*pfpar(pft)%b1Cd*1000.d0))
     &     **(1.0d0/pfpar(pft)%b2Cd)

      end function Cdead2dbh

!*************************************************************************
      real*8 function Cdead2Chw(pft,Cdead)
      integer :: pft
      real*8 :: Cdead
!@sum (gC/plant) Return portion of Cdead (structural carbon) that is 
!@+     Chw (heartwood stem carbon)

      Cdead2Chw = Cdead /(pfpar(pft)%croot_ratio+1.d0)
      end function Cdead2Chw

!*************************************************************************
      real*8 function Chw2Ccroot(pft,Chw)
      integer :: pft
      real*8 :: Chw
!@sum (gC/plant) Return Croot (coarse root carbon) given Chw (heartwood C)
      Chw2Ccroot = Chw * pfpar(pft)%croot_ratio
      end function Chw2Ccroot

!*************************************************************************
      real*8 function dDBHdCdead(pft,Cdead)
!@sum (cm/kg)  Derivative w.r.t. Cdead of Cdead2dbh.
      integer,intent(in) :: pft
      real*8, intent(in) :: Cdead

!      dDBHdCdead=(C2B/1000.0d0/pfpar(pft)%b1Cd)**(1.0d0/pfpar(pft)%b2Cd)
!     &          *Cdead**((1.0d0/pfpar(pft)%b2Cd)-1.0d0)
!     &          /pfpar(pft)%b2Cd

      dDBHdCdead=(1/(B2C*1000.0d0*pfpar(pft)%b1Cd))**
     &     (1.0d0/pfpar(pft)%b2Cd)
     &     *Cdead**((1.0d0/pfpar(pft)%b2Cd)-1.0d0)
     &     /pfpar(pft)%b2Cd

      end function dDBHdCdead

!*************************************************************************
      real*8 function dbh2Cfol(pft,dbh)
!@sum (gC/plant) Maximum carbon in foliage as a function of dbh
!@+   After formulation for cold-deciduous trees in Albani et al. (2006).
!@+   Other Ent PFTs, see full Ent documentation.
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8 :: dbhmax

      dbhmax=maxdbh(pft)

c$$$      if (.not.pfpar(pft)%woody) then !herbaceous
c$$$         maxdbh=log(1.0-(0.999*(pfpar(pft)%b1Ht+1.3)-1.3)  
c$$$     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
c$$$      else !woody
c$$$         maxdbh=log(1.0-(0.999*pfpar(pft)%b1Ht-1.3)  
c$$$     &        /pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
c$$$      end if

      dbh2Cfol=B2C *pfpar(pft)%b1Cf *1000.0d0 !kgC to gC
     &        * min(dbh, dbhmax)**pfpar(pft)%b2Cf

      end function dbh2Cfol
!*************************************************************************
 
      real*8 function dDBHdCfol(pft,Cfol)
!@sum Change in DBH (cm) given change in Cfol. Derivative of dbh2Cfol.
!@+   For woody plants only, of course.
      integer,intent(in) :: pft
      real*8, intent(in) :: Cfol !(gC/plant) foliage carbon

      dDBHdCfol=
     &   (1.d0/(B2C*1000.0d0*pfpar(pft)%b1Cf))**(1.0d0/pfpar(pft)%b2Cf)
     &     *Cfol**((1.0d0/pfpar(pft)%b2Cf)-1.0d0)
     &     /pfpar(pft)%b2Cf

      end function dDBHdCfol
!*************************************************************************

      real*8 function iqsw_fn(pft)
!@sum iqsw: sapwood biomass per (leaf area x wood height) (kgC/m2/m) 
!@+   Values from ED):
!@+      !3900.0: leaf area per sapwood area (m2/m2) 
!@+      !1000.0: sapwood density (kgB/m3)
!@+      !2.0:  biomass per carbon (kgB/kgC)
      !(qsw)=(iqsw*sla) (1/m) & (qsw*h): ratio of sapwood to leaf biomass (unitless)
      !(iqsw)=1000.0d0/3900.0d0/2.0d0=0.1282 !NOTE: This value corrects an error in the coefficient in Moorcroft et al. (2001) Appendix D, which had the value too small by a factor of 100, at 0.00128.

      !use ent_prescr_veg, only : wooddensity_gcm3
      integer, intent(in) :: pft
      !--- Local ---
      real*8, parameter :: iqsw0=1000.0d0/3900.0d0/2.0d0
      real*8 :: rho_wood        !kgC/m3
      real*8 :: Csw

      !Yeonjoo's
      !iqsw = iqsw0*wooddensity_gcm3(pft)

      !wdens_g_cm3 is static array
!      rho_wood = wdens_g_cm3(pft) * 1000.d0 !gB/cm3 to kgB/m3
      !wooddensity_gcm3 is f(lrage or llspan), but prescr. needle/tundra
      rho_wood = wooddensity_gcm3(pft) * 1000.d0 !gB/cm3 to kgB/m3
      iqsw_fn = B2C*rho_wood/3900.d0 !kgCwood/m3wood * m2sw/m2fol = kgCsw/m2fol/msw

      end function iqsw_fn
!*************************************************************************
      real*8 function Csw_fn(pft, DBH_cm, h_m, la)
!@sum (gC/plant) Carbon in sapwood.
      !Csw = iqsw (kgCsw/m2fol/msw) * SLA (m2fol/kgCfol) * h (msw) * Cfol (kgC)
      !    = kgCsw or 1000.*gCsw
      !0.5 converts density from dry biomass to C, both rho_wood and Cfol
      !From data of Cleary et al (2008) for sagebrush, it appears that
      ! the piperatio 3900 = h*factor.  The factor varies with h
      ! to keep h*factor constant.  So I will look at cutting h out of eqn.-nk
      integer :: pft
      real*8 :: DBH_cm, h_m, la !la = current lai, not max
      !-------
      real*8 :: rho_wood  !kgC/m3
      real*8 :: Csw
      real*8 :: iqsw

      iqsw = iqsw_fn(pft) !kgCsw/m2fol/msw
      !rho_wood = wdens_g_cm3(pft) * 1000. !dry mass, convert gB/cm3 to kgB/m3
      !iqsw = B2C*rho_wood/3900. !kgCwood/m3wood * m2sw/m2fol = kgCsw/m2fol/msw
	
      if (pfpar(pft)%woody) then
!         if (DBH_cm.lt.10.) then

!         else
!           Csw = h * qsw* Cfolmax !bad old way to calculate, qsw=sla*iqsw
!           Csw = iqsw * pfpar(pft)%sla * h_m
!     &        * B2C * pfpar(pft)%b1Cf*(DBH_cm**pfpar(pft)%b2Cf)
!     &        *1.d3             !kgC to gC
         Csw = iqsw * h_m * la * 1.d3 !kgC to gC
!         endif
      else 
         Csw = (0.0)
      endif

      Csw_fn = Csw
      end function Csw_fn


!*************************************************************************
      real*8 function Cfol_fn(pft,DBH_cm, h_m) 
!@sum (gC/plant) Maximum carbon in foliage for woody or herbaceous plants,
!@+   given plant geometry.
      integer :: pft
      real*8 :: DBH_cm, h_m
      !-------
      integer :: vform
      real*8 :: rho_wood
      real*8  :: Cfol !gC/plant

      vform = gform(pft)

      if (pfpar(pft)%woody) then
         rho_wood = wdens_g_cm3(pft) * 1000. !convert from g/cm3 to kg/m3
         !Csw = Csw_fn(pft,DBH_cm,h_m)
         if ((h_m.gt.0.0).and.(DBH_cm.gt.0.0)) then
            !* Based on foliage dependence on sapwood.- NYK
            !Cfol = (Csw*3900./(0.5*rho_wood*sla_m2_kgC(pft)*h_m))
            ! Yeonjoo's
            Cfol = dbh2Cfol(pft,DBH_cm)
         else
            Cfol = 0.0
         endif
      elseif ((vform.eq.HERB).or.(vform.eq.GRASS)) then
         if (h_m.gt.0.0) then
            !This departs from Yeonjoo's height2Cfol, instead uses
            !Nancy's formulation.
            !Cfol = 0.5*(ag(pft)*h_m**bg(pft))
            Cfol = B2C*(pfpar(pft)%b1Cf*h_m**pfpar(pft)%b2Cf)*1000.d0
         else
            Cfol = 0
         endif
      else
         Cfol = (0.0)
      endif
      if (isNaN(Cfol)) then
         write(*,*) "Cfol=NaN",pft,vform, rho_wood,DBH_cm,h_m
     &        ,pfpar(pft)%sla,pfpar(pft)%woody
     &        ,pfpar(pft)%b1Cf,pfpar(pft)%b2Cf!,ag(pft),bg(pft)
         !STOP
      endif
      Cfol_fn = Cfol
      end function Cfol_fn

!*************************************************************************
      real*8 function Cfol2height_herb(pft,Cfol)
!@sum Update herb height (m), given new foliage carbon.
      integer,intent(in) :: pft
      real*8, intent(in) :: Cfol !gC/pool/plant
      real*8,parameter :: h1Cf = 1.66d0 
      real*8,parameter :: h2Cf = 1.5d0
      real*8,parameter :: nplant = 2500.d0

      if (Cfol.gt.0.0d0) then
         Cfol2height_herb=
     &        exp(log(Cfol*1.d0/(B2C*h1Cf)*nplant)/h2Cf)/100.0d0
      else
         Cfol2height_herb=0.d0
      end if
      return
      end function Cfol2height_herb

!*************************************************************************

      real*8 function dbh2height(pft,dbh)
!@sum (m) Height of woody plant given dbh (for shrubs, dbh~basal diameter).
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      !----
      real*8 :: hcrit !max height

      hcrit=a0h(pft)+pfpar(pft)%b1Ht-EPS
      dbh2height = min(a0h(pft) + pfpar(pft)%b1Ht * 
     &             (1.0-exp(pfpar(pft)%b2Ht*dbh)), hcrit)

      end function dbh2height
!*************************************************************************
      real*8 function dHdDBH(pft,dbh)
!@sum Change in height (m) given change in DBH (cm). Derivative of dbh2height.
!@+   For woody plants only, of course.
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh !Diameter at breast height (cm)

      dHdDBH = - pfpar(pft)%b1Ht*pfpar(pft)%b2Ht
     &         * exp(pfpar(pft)%b2Ht * dbh)
      end function dHdDBH

!*************************************************************************

      real*8 function height2dbh(pft,h)
!@sum (cm) Woody plant diameter (cm) given height (m).
!@+        trees:  diameter at breast height (dbh), shrubs: basal diameter.
!@+   From Albani et al. (2006)
!@+   Fixed to avoid negative dbh for h=0.
      !Note: a0h>0 is problematic, since get non-zero height for dbh=0.
      !   !!!Need to fix allometry for pfts with a0h>0!!!
      !NK -07/02/2012
      !Original relation:  
      !  h_m = a0h + pfpar(pft)%b1Ht*(1-exp(pfpar(pft)%b2Ht*DBH_cm))
      !Inverse relation:   
      !  DBH_cm = log(1-(h_m - a0h)/pfpar(pft)%b1Ht)/pfpar(pft)%b2Ht
      !This places a cap on DBH if h_m goes above pfpar(pft)%b1Ht, 
      ! which is really wrong, since it is h that is limited while
      ! DBH gets fatter -
      !  but we don't know what the upper limit of h is.
      !  Need to change function unlimited form h = a*dbh^b
  
      integer,intent(in) :: pft
      real*8, intent(in) :: h
      real*8 :: hcrit, hin

      if (pfpar(pft)%woody) then
         !ED version, Moorcroft, et al. (2001)
         !dbh = ((1/2.34d0)*h)**(1/0.64d0)

         hcrit=a0h(pft)+pfpar(pft)%b1Ht-EPS
         hin=min(h,hcrit)
         
         height2dbh = max(0.d0, 
     &        (1.d0/pfpar(pft)%b2Ht)
     &        *log(1.d0-(hin-a0h(pft))/pfpar(pft)%b1Ht))
      else
         height2dbh = 0.d0 !herb
      endif

      end function height2dbh

!*************************************************************************
!Keep for development purposes - NK
!      real*8 function height2Cfol_herb(pft,height)
!!@sum (gC/plant) Maximum carbon in foliage for herbaceous plant
!      !Yeonjoo's.  No documentation or units.
!      !?? Is this only for annual herbaceous recruitment?
!      integer,intent(in) :: pft
!      real*8, intent(in) :: height
!      real*8,parameter :: h1Cf = 1.66d0 
!      real*8,parameter :: h2Cf = 1.50d0 
!      real*8,parameter :: nplant = 2500.d0 
!      height2Cfol= B2C *h1Cf 
!     &        *((height*100.0d0)**h2Cf)/nplant
!
!      end function height2Cfol_herb
!*************************************************************************

      real*8 function Crown_rad_allom(pft, h_m)
!@sum Return crown radius (m) calculate from allometry.
!@auth NYK
      integer :: pft
      real*8 :: h_m
      Crown_rad_allom = acr(pft)*h_m

      end function Crown_rad_allom

!*************************************************************************
      real*8 function Crown_rad_max_from_density(ndens) 
!@sum Return potential crown horizontal radius (m) based on maximum packing:
!@+   n = pi/sqrt(12)/(2*(a+da))^2  #Lagrange hexagonal maximum packing density
!@+   where n = density, a=crown radius, da = distance between crowns
!@+      a = 0.5*sqrt(pi/sqrt(12)/n) - da
      real*8 :: ndens !Plant density in #/m^2

      if (ndens.gt.0.d0) then
         Crown_rad_max_from_density = 0.5d0*sqrt(pi/sqrt(12.d0)/ndens)
      else
         Crown_rad_max_from_density = 0.d0
      endif

      end function Crown_rad_max_from_density

!*************************************************************************
      real*8 function crown_radius_horiz_allom(pft,h,popdensity)
!@sum Horizontal crown radius (m), minimum from allometry or max packing.
      integer,intent(in) :: pft
      real*8, intent(in) :: h  !m
      real*8, intent(in) :: popdensity !#/m^2

      if (pfpar(pft)%woody) then
         crown_radius_horiz_allom = min(
     &        Crown_rad_max_from_density(popdensity)
     &        ,Crown_rad_allom(pft,h))
      else !Herbs based on density
         crown_radius_horiz_allom =
     &        Crown_rad_max_from_density(popdensity)
      endif

      end function crown_radius_horiz_allom
!*************************************************************************
      real*8 function crown_radius_horiz_HF(pft,dbh,popdensity)
!@sum Horizontal crown radius (m), minimum from allometry or max packing.
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh !cm
      real*8, intent(in) :: popdensity !#/m^2

      crown_radius_horiz_HF = min(Crown_rad_max_from_density(popdensity)
     &     ,crown_radius_pft_HF(pft,dbh)) !!Yeonjoo Harvard Forest hack

      end function crown_radius_horiz_HF
!*************************************************************************
      real*8 function crown_radius_pft_HF(pft,dbh) Result(cradm)
!@sum Horizontal crown radius from Harvard Forest late successional 
!@+   hardwood allometry with mean conifer dbh_max limit.
!@+   Coefficient 0.107 for late-succ hw is approx. mean for all types.
      use ent_pfts, only : is_conifer,is_hw
      integer, intent(in) :: pft
      real*8, intent(in) :: dbh
      !------
      !dbh_max parameter values from Harvard Forest allometry
      real*8, parameter :: dbh_max_conifer = 42.d0 !cm.  Mean of for pine and late successional conifer.  
      real*8, parameter :: dbh_max_hw = 150.d0 !cm.  
      real*8 :: dbh_max

      if (is_conifer(pft)) then
         dbh_max = dbh_max_conifer
      elseif (is_hw(pft)) then
         dbh_max = dbh_max_hw
      else
         dbh_max = dbh
      endif

      cradm = .107d0 * min(dbh, dbh_max) 

      end function crown_radius_pft_HF
!*************************************************************************
!      real*8 function crown_radius_closed(popdensity) Result(cradm)
!      !* Return plant crown radius (m).
!      !* Assumes closed canopy packing given popdensity, square regular distrib.
!      real*8, intent(in) :: popdensity !#/m^2
!      
!      cradm = 0.5d0*sqrt(1/popdensity)
!      end function crown_radius_closed
!*************************************************************************

      real*8 function crown_radius_vert(pft, h,crx)
!@sum Vertical crown radius (m) from allometry.
!@+   Subject to change.  Currently allows tall ellipsoid to spherical but not
!@+   oblate crowns. May want to allow oblate for understory crowns.
      integer :: pft
      real*8 :: h, crx !Tree height, crown horizontal radius

      if (form(pft).eq.HERB) then
         crown_radius_vert = 0.5d0*h  !Entire height of plant
      else !TREE, SHRUB
         !crown_radius_vert = min(0.45*h,crx*2.7d0)
         crown_radius_vert = max(0.45*h,crx) !##This does not allow oblate shape.
      endif
      end function crown_radius_vert

!*************************************************************************

      subroutine allom_plant_cpools(pft,lai,h,dbh,popdens,cpool)
!@sum Calculate plant carbon pools for single plant (g-C/plant)
!@+   based on geometry allometry and LAI.
!@+   Does NOT update LABILE, which is prognostic, not static allometry.
      !* Coarse root fraction estimated from Zerihun (2007) for Pinus radiata
      !*  at h=20 m for evergr, 50% wood is carbon, and their relation:
      !*  CR(kg-C/tree) = 0.5*CR(kg/tree) = 0.5*exp(-4.4835)*(dbh**2.5064).
      !*  The ratio of CR(kg-C/tree)/HW(kg-C/tree) at h=20 m is ~0.153.
      !*  This is about the same as the ratio of their CR biomass/AG biomass.  
      use ent_pfts, only: pfpar !COVEROFFSET, alamax
!      use allometryfn, only : Cfol_fn, Csw_fn ,dbh2Cdead
      integer,intent(in) :: pft !plant functional type
      !real*8, intent(in) :: laimax !max lai should not be needed here
      real*8, intent(in) ::lai  !lai
      real*8, intent(in) :: h,dbh,popdens !h(m), dbh(cm),popd(#/m2)
      real*8, intent(out) :: cpool(N_BPOOLS) !g-C/pool/plant
      !----Local-----
      real*8 :: LA,LAmax

      !* Initialize
      cpool(:) = 0.d0 

      if (popdens.eq.0.d0) then
         LA = 0.d0
         cpool(FOL) = 0.d0
      else
         LA = lai/popdens
         cpool(FOL) = LA/pfpar(pft)%sla*1.d3 !gC/plant
      endif
      LAmax = Cfol_fn(pft,dbh,h)*(pfpar(pft)%sla*.001d0)
      !max_cpoolFOL = alamax(pft+COVEROFFSET)/pfpar(pft)%sla/popdens*1d3
      !max_cpoolFOL = laimax/pfpar(pft)%sla/popdens*1.d3
      cpool(FR) = q * cpool(FOL)   !Br
      if (pfpar(pft)%woody) then !Woody
         !cpool(SW) = Csw_fn(pft,dbh,h,max(LA,0.5d0*LAmax))
         cpool(SW) = Csw_fn(pft,dbh,h,LAmax)
         cpool(HW) = Cdead2Chw(pft, dbh2Cdead(pft,dbh))
        !cpool(HW) = dbh2Cdead(pft,dbh) * hw_fract
        !cpool(CR) = cpool(HW) * (1-hw_fract)/hw_fract !=dbh2Cdead*(1-hw_fract)
         cpool(CR) =  Chw2Ccroot(pft,cpool(HW))
      else
         cpool(SW) = 0.d0
         cpool(HW) = 0.d0
         cpool(CR) = 0.d0
      endif

      end subroutine allom_plant_cpools

!**************************************************************************
      subroutine init_Clab(pft,dbh,h,Clabile)
!@sum init_Clab - renamed from prescr_init_Clab (only depends on allometry)
!@sum - Initializes labile carbon pool
!@sum Deciduous woody: Clab = 4 x max Cfol of plant.
!@sum Evergreen woody: Clab = 0.5 x max Cfol of plant.
!@sum 4x requirement is from Bill Parton (personal communication).
!@sum 7/1/2012 - Replaced with different subroutine - NK
      use ent_pfts, only : COVEROFFSET, pfpar !, alamax, alamin
      implicit none
      integer, intent(in) :: pft
      real*8, intent(in) :: dbh !dbh (cm) !0 for herbs
      real*8, intent(in) :: h !height (m)
      !real*8, intent(in) :: n   !Density (#/m^2)
      !real*8, intent(in) :: laimax
      real*8, intent(inout) :: Clabile !g-C/pool/plant
      !----- Local ------
      real*8 :: Cfolmax

      !Cfolmax = laimax/n/pfpar(pft)%sla !kgC/plant - old way based on laimax
      Cfolmax = Cfol_fn(pft,dbh,h) !g-C/individ -- CHECK Cfol_fn!!!

      if (pfpar(pft)%phenotype.eq.EVERGREEN) then
        !Enough to grow peak foliage and fine roots, and some reproduction.
        Clabile = Cfolmax*2.d0
      else
        Clabile = Cfolmax*4.d0
      endif

      end subroutine init_Clab
!*************************************************************************

      real*8 function nplant(pft,dbh,h,laimax)
!@sum (#plants/m^2) Plant density given individual geometry and patch LAImax.
      integer,intent(in) :: pft
      real*8, intent(in) :: dbh
      real*8, intent(in) :: h
      real*8, intent(in) :: laimax

      if (form(pft).eq.HERB) then !grasses/crops/non-woody 
!         nplant = laimax/pfpar(pft)%sla/(height2Cfol(pft,h)/1000.0d0) !YK's
         nplant = laimax/pfpar(pft)%sla/(Cfol_fn(pft,0.d0,h)/1000.0d0) 
      else
         nplant = laimax/pfpar(pft)%sla/(dbh2Cfol(pft,dbh)/1000.0d0)
      end if

      end function nplant


!**************************************************************************
!## This function does not appear to be used - NK
!      real*8 function prescr_calc_shoot(pft,hdata,dbhdata)
!     &     Result(Bshoot)
!!@sum Returns GISS GCM veg shoot kg-C per plant for given vegetation type.
!!@+   From Moorcroft, et al. (2001), who takes allometry data from
!!@+   Saldarriaga et al. (1998).
!     use ent_pfts, only : COVEROFFSET
!      integer,intent(in) :: pft !@var pft vegetation type
!      real*8,intent(in) :: hdata(N_COVERTYPES), dbhdata(N_COVERTYPES)
!      !-----Local-------
!      real*8 :: wooddens
!      integer :: ncov !covertypes index
!
!      ncov = pft + COVEROFFSET
!      wooddens = wooddensity_gcm3(pft)
!      Bshoot = 0.069 * (hdata(ncov))**0.572
!     &     * (dbhdata(ncov))**1.94 * (wooddens**0.931)
!
!      end function prescr_calc_shoot
!**************************************************************************

      end module allometryfn

