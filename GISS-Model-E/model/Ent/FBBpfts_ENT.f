      module FarquharBBpspar
!@sum Module  Ent plant functional type (PFT) parameters for Ent 16 PFTs default.
!@+   for Farqhuar-von Caemmerer (1982) photosynthesis and 
!@+   Ball-Berry (1985) stomatal conductance.
      use ent_const, only : N_PFT

      implicit none

      !======DECLARED TYPES====== !
      type pspartype
      integer :: pst            !Photosynth type.  1-C3, 2=C4
      real*8 :: PARabsorb       !Leaf PAR absorptance (fraction)

      !Photosynthesis/Conductance - Farquhar/Ball-Berry parameters
      real*8 :: Vcmax           !Maximum photosynthetic capacity (umol m-2 s-1)
!      real*8 :: Kc              !Michaelis-Menten constant for CO2 (Pa)
!      real*8 :: Ko              !Michaelis-Menten constant for O2 (Pa)
!      real*8 :: KcQ10           !Kc Q10 exponent
!      real*8 :: KoQ10           !Ko Q10 exponent
      !real*8 :: GammastarQ10    !CO2 compensation point Q10 (Pa)
      real*8 :: m               !Slope of Ball-Berry equation
      real*8 :: b               !Intercept of Ball-Berry equation (mol m-2 s-1)
!      real*8 :: Rdc             !Dark respiration c scaling factor (Harley&Tenhunen, 1991)
!      real*8 :: RdH             !Dark respiration deltaH (J mol-1) (Harley&Tenhunen, 1991)
      real*8 :: Nleaf           !g-N/m2[leaf]
      end type pspartype


      type psdrvtype
      real*8 :: ca              !Surface CO2 mole fraction (umol mol-1)
      real*8 :: ci              !Leaf internal CO2 mole fraction (umol mol-1)
      real*8 :: Tc              !Canopy (foliage) temperature (Celsius)
      real*8 :: Pa              !Atmospheric pressure (Pa)
      real*8 :: rh              !Relative humidity (fraction)
      end type psdrvtype


      !=======CONSTANTS========!
!      integer,parameter :: N_PFT = 16 !In ent_const.f

!*********************************************************************
!* Ent PFTs
!* 1.  evergreen broadleaf early successional
!* 2.  evergreen broadleaf late successional
!* 3.  evergreen needleleaf early successional
!* 4.  evergreen needleleaf late successional
!* 5.  cold deciduous broadleaf early successional
!* 6.  cold deciduous broadleaf late successional
!* 7.  drought deciduous broadleaf
!* 8.  decidous needleleaf
!* 9.  cold adapted shrub
!* 10.  arid adapted shrub
!* 11.  C3 grass - perennial
!* 12.  C4 grass - perennial
!* 13.  C3 grass - annual
!* 14.  arctic C3 grass
!* 15.  crops - C4 herbaceous
!* 16.  crops - broadleaf woody
!*********************************************************************

      type(pspartype),parameter :: pftpar(N_PFT) = !PFT parameters for Ent veg types
     &!     pft PARabsorb Vcmax Kc Ko KcQ10 KoQ10 Gammastar  m b !Rdc RdH
     &     (/
     &     pspartype(1          !1. EVERGREEN BROADLEAF EARLY SUCCESSIONAL
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BET temperate & tropical, Table 3.1 (Oleson, et al 2004)
     &     ,75.d0               !Vmax25, CLM BET tropical, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0,.002d0         !m,b, CLM
!     &     ,6.0d0               !Nleaf (gN/m2-leaf). Low Est. from Reich 1997 (big range).
     &     ,2.7d0),               !Nleaf (gN/m2-leaf). Friend&Kiang (2005), Table 1.

     &     pspartype(1          !2. EVERGREEN BROADLEAF LATE SUCCESSIONAL
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BET & BDT temperate & tropical, Table 3.1 (Oleson, et al 2004)
!     &     ,69.d0               !Vmax25, CLM BET temperate, Table 8.2 (Oleson, et al 2004)
!     &     ,54.d0               !Tapajo KM67
     &     ,50.d0               !
     &     ,9.d0,.002d0          !m,b, CLM
!     &     ,8.0d0               !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ,2.5d0),               !Nleaf (gN/m2-leaf). fort.995=read.table("../Code/GISSClim/GISSClim/WorkSpace/r_ent_fbb_tnfg/fort.995", header=TRUE)

     &     pspartype(1          !3. EVERGREEN NEEDLELEAF EARLY SUCCESSIONAL
     &     ,.93d0               !from leaf VIS 1-albedo,CLM NET & NDT temperate & boreal, Table 3.1 (Oleson, et al 2004)
     &     ,51.d0               !Vmax25, CLM NET temperate, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0,.002d0         !m,b, CLM
!     &     ,2.8d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
!     &     ,2.9d0               !Nleaf (gN/m2-leaf). Friend&Kiang (2005), Table 1.
!     &     ,0.66d0),                 !Nleaf (gN/m2-leaf). Reich (1997) for SLA of 18 m2/kg. This give C:N of 84!  Biome-BGC is 42, CLM 29.
     &     ,1.8d0),             !Nleaf (gN/m2-leaf). Pinus sylvestris, Pensa and Sellin (2002).

     &     pspartype(1          !4. EVERGREEN NEEDLELEAF LATE SUCCESSIONAL
     &     ,.93d0               !from leaf VIS 1-albedo,CLM NET & NDT temperate & boreal, Table 3.1 (Oleson, et al 2004)
     &     ,43.d0               !Vmax25, CLM NET boreal, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0,.002d0              !m,b, CLM
!     &     ,3.0d0               !Nleaf (gN/m2-leaf). High Est. from Reich 1997 (big range).
!     &     ,0.66d0),            !Nleaf (gN/m2-leaf). Reich (1997) for SLA of 18 m2/kg. This give C:N of 84!  Biome-BGC is 42, CLM 29.
     &     ,1.8d0),             !Nleaf (gN/m2-leaf). Pinus sylvestris, Pensa and Sellin (2002).

     &     pspartype(1          !5. COLD DECIDUOUS BROADLEAF EARLY SUCCESSIONAL
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDT temperate, Table 3.1 (Oleson, et al 2004)
     &     ,60.d0               !Vmax25, CLM BDT temperate, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0  ,.002d0              !m,b, CLM
!     &     ,6.0d0               !Nleaf (gN/m2-leaf). Low Est. from Reich 1997 (big range).
!     &     ,3.1d0               !Nleaf (gN/m2-leaf), Q. ilex, Mediavilla & Escudero(2003)
     &     ,1.07d0              !Nleaf (gN/m2-leaf), Oak, derived from C:N 27.2 and SLA 34.5 m2/kg-C in Tatarinov & Cienciala (2006) for BIOME-BGC.
     &     ),                   

     &     pspartype(1          !6. COLD DECIDUOUS BROADLEAF LATE SUCCESSIONAL
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDT temperate, Table 3.1 (Oleson, et al 2004)
     &     ,51.d0               !Vmax25, CLM BDT boreal, Table 8.2 (Oleson, et al 2004). Wilson et al (2001) Walker Branch species average to 51.
!     &     ,40.d0               !Vmax25, tuned to MMSF, was 40.0; changed to 39.0 8/9/11-NK
     &     ,9.d0 ,.002d0              !m,b, CLM
!     &     ,6.7d0                !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ,1.25d0),               !Nleaf (gN/m2-leaf). Friend&Kiang (2005), Table 1.

     &     pspartype(1          !7. DROUGHT DECIDUOUS BROADLEAF 
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDT temperate & tropical, Table 3.1 (Oleson, et al 2004)
     &     ,56.4d0              !Vcmax25, Wang et al (2007) GCB Fluxnet inversion avg of grass+tree 50.1 & 56.4 umol m-2 s-1
!     &     ,100.d0              !Vcmax, Liukang Xu's Tonzi Ranch blue oak leaf measurements go as high as 100 umol m-2 s-1, highly seasonal. See Kiang dissertation (2002) Fig. 2.34.
!     &     ,17.d0               !Vmax25, CLM BES & BDS temperate, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0 ,.002d0              !m,b, CLM
!     &     ,6.7d0              !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
!     &     ,1.5d0),               !Nleaf (gN/m2-leaf). Increment Friend&Kiang (2005), Table 1.
!     &     ,3.1d0                  !Nleaf (gN/m2-leaf), Q. ilex, Mediavilla & Escudero(2003)
     &     ,1.07d0              !Nleaf (gN/m2-leaf), Oak, derived from C:N 27.2 and SLA 34.5 m2/kg-C in Tatarinov & Cienciala (2006) for BIOME-BGC.
     &     ), 

     &     pspartype(1          !8. DECIDUOUS NEEDLELEAF
     &     ,.93d0               !from leaf VIS 1-albedo,CLM NDT boreal, Table 3.1 (Oleson, et al 2004)
     &     ,43.d0               !Vmax25, CLM NDT boreal, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0  ,.002d0              !m,b, CLM
     &     ,1.2d0),               !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range). Avg 2.8 (evergr) and 6.7 (decid).

     &     pspartype(1          !9. COLD ADAPTED SHRUB (TUNDRA)
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDS boreal, Table 3.1 (Oleson, et al 2004)
     &     ,33.d0               !Vmax25, CLM BDS boreal, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0 ,.002d0              !m,b, CLM
!     &     ,7d0                 !Nleaf (gN/m2-leaf). Est. from Reich 1997 (big range).
     &     ,1.6d0),               !Nleaf (gN/m2-leaf). Friend&Kiang (2005), Table 1.

     &     pspartype(1          !10. ARID ADAPTED SHRUB
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDS temperate, Table 3.1 (Oleson, et al 2004)
     &     ,17.d0               !Vmax25, CLM BES & BDS temperate, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0  ,.002d0              !m,b, CLM
!     &     ,7d0                 !Nleaf (gN/m2-leaf). =tundra. Est. from Reich 1997 (big range).
     &     ,2.38d0),               !Nleaf (gN/m2-leaf). Friend&Kiang (2005), Table 1.

     &     pspartype(1          !11. GRASSC3 - perennial
     &    ,.86d0                !PARabsorb, Collatz et al. (1991)
!     &     ,.89d0               !from leaf VIS 1-albedo,CLM C3 grass, Table 3.1 (Oleson, et al 2004)
!     &     ,56.4d0              !Vcmax25, Wang et al (2007) GCB Fluxnet inversion avg of grass+tree 50.1 & 56.4 umol m-2 s-1
!     &     ,60d0                !Vcmax, von Caemmerer, CSIRO 2000, VARIOUS VALUES
!     &     ,93.d0               !Vcmax, S. Verma and J. Berry,http://nigec.ucdavis.edu/publications/annual97/greatplains/project86.html
     &     ,43.d0   ,11.d0               !m, b, X.Mo, et al. (2001)
     &     ,.008d0              !b, X.Mo, et al. (2001)
     &     ,2.460),              !Nleaf (gN/m2-leaf). Ponca winter wheat, from C:N 18.85, SLA 21.6 m2/kg-C.

     &     pspartype(2          !12. GRASSC4 - perennial
     &     ,.9d0               !leaf VIS albedo,CLM C4 grass, Table 3.1 (Oleson, et al 2004)
     &     ,24.d0               !Vmax25, CLM C4 grass, Table 8.2 (Oleson, et al 2004)
     &     ,5d0                !m, CLM C4 grass, Table 8.2 (Oleson, et al 2004)
     &     ,.002d0              !b, CLM (Oleson, et al 2004, Section 8, p. 129)
     &     ,2.1d0),             !Nleaf (gN/m2-leaf), from Shidler, C4 grass.  Nleaf=1.0 from Reich (1997) C:N of ~24 for SLA of 41.1 m2/kg-C (leaf longevity 0.5 yr).  Ponca has 3.27d0 g-N/m2-green leaf (seasonal curve), avg C:N of 20 kg/kg, avg SLA 9.4 m2/kg-C (fat leaves).

     &     pspartype(1          !13. GRASSC3 - annual
     &    ,.86d0                !PARabsorb, Collatz et al. (1991)
!     &     ,.89d0               !leaf VIS albedo,CLM C3 grass, Table 3.1 (Oleson, et al 2004)
!     &     ,43.d0               !Vmax25, CLM C3 grass, Table 8.2 (Oleson, et al 2004)
     &     ,60.d0               !Vcmax, best fit guess, see plots of 08/10/2007.
!     &     ,93.d0               !Vcmax for Vaira, S. Verma and J. Berry,http://nigec.ucdavis.edu/publications/annual97/greatplains/project86.html
!     &     ,9d0                !m, CLM C3 grass, Table 8.2 (Oleson, et al 2004)
!     &     ,.002d0              !b, CLM (Oleson, et al 2004, Section 8, p. 129)
     &     ,11.d0               !m, X.Mo, et al. (2001)
     &     ,.008d0              !b, X.Mo, et al. (2001)
     &     ,2.46d0),              !Nleaf (gN/m2-leaf) from Ponca. Ponca has 3.27d0 g-N/m2-green leaf (seasonal curve), avg C:N of 20 kg/kg, avg SLA 9.4 m2/kg-C dry wt, or SLA 21.6 m2/kg-C +- 9.6.

     &     pspartype(1          !14. GRASSC3 - arctic
     &     ,.89d0               !from leaf VIS 1-albedo,CLM C3 grass, Table 3.1 (Oleson, et al 2004)
     &     ,43d0               !Vmax25, CLM C3 arctic grass, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0               !m, CLM C3 arctic grass, Table 8.2 (Oleson, et al 2004)
     &     ,.002d0             !b, CLM (Oleson, et al 2004, Section 8, p. 129)
!     &     ,1.15d0),              !Nleaf (gN/m2-leaf) Increment from Reich (1997) C:N of ~24 for SLA of 41.1 m2/kg-C (leaf longevity 0.5 yr).  Ponca has 3.27d0 g-N/m2-green leaf (seasonal curve), avg C:N of 20 kg/kg, avg SLA 9.4 m2/kg-C (fat leaves).
     &     ,2.46d0),              !Ponca Ntot/LA average (actually seasonal curve); Rd not large deviance from direct daily Ntot 

     &     pspartype(2          !15. CROPS - C4
     &     ,.89d0               !from leaf VIS 1-albedo,CLM Crop1 & Crop2, Table 3.1 (Oleson, et al 2004)
     &     ,50d0               !Vmax25, CLM Crop1, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0               !m, CLM Crop1, Table 8.2 (Oleson, et al 2004)
     &     ,.002d0              !b, CLM (Oleson, et al 2004, Section 8, p. 129)
     &     ,2.5d0),                 !Nleaf (gN/m2-leaf). Friend&Kiang (2005), Table 1.
!     &     ,4.0d0),              !Round up high from Ponca Ntot/LA average (actually seasonal curve); Rd not large deviance from direct daily Ntot 

     &     pspartype(1          !16. CROPS - BROADLEAF WOODY
     &     ,.90d0               !from leaf VIS 1-albedo,CLM BDT, Table 3.1 (Oleson, et al 2004)
     &     ,51.d0               !Vmax25, CLM BDT temperate, Table 8.2 (Oleson, et al 2004)
     &     ,9.d0               !m, CLM BDT temperate, Table 8.2 (Oleson, et al 2004)
     &     ,.002d0              !b, CLM (Oleson, et al 2004, Section 8, p. 129)
     &     ,3.1d0)               !Nleaf (gN/m2-leaf), Q. ilex, Mediavilla & Escudero(2003)
     &/)


      !NOTES:
        !--------Collatz, et al. (1991) Farquhar parameters-----
        !Vmax values:
        !Collatz, C3 grass. Vcmax = 200.
        !Harley, et al. (1992), cotton Vcmax=51-127 umol m-2 s-1
        !Ponca, Oklahoma (Fluxnet), winter wheat, Vcmax fit ~30.
!        pspar%Vmax = 200./(1 + exp((-220.e03+703.*(Tl+Kelvin))



!****************************************************************************
      end module FarquharBBpspar
