      module ent_pfts
!@sum Parameter set for Ent default supported 16 plant functional types
!@+   tailored to Fluxnet site Ent_standalone runs.

      !use ent_pftconst
      use ent_const
      use ent_types
      implicit none

      !***************************************************
      !*      ENT PLANT FUNCTIONAL TYPES                 *
      !***************************************************

      character*50, parameter :: Ent_title(N_COVERTYPES) =
     &     (/
     &     '1 - evergreen broadleaf early succ               ',
     &     '2 - evergreen broadleaf late succ                ',
     &     '3 - evergreen needleleaf early succ              ',
     &     '4 - evergreen needleleaf late succ               ',
     &     '5 - cold deciduous broadleaf early succ          ',
     &     '6 - cold deciduous broadleaf late succ           ',
     &     '7 - drought deciduous broadleaf                  ',
     &     '8 - deciduous needleleaf                         ',
     &     '9 - cold adapted shrub                           ',
     &     '10 - arid adapted shrub                          ',
     &     '11 - C3 grass perennial                          ',
     &     '12 - C4 grass                                    ',
     &     '13 - C3 grass - annual                           ',
     &     '14 - arctic C3 grass                             ',
     &     '15 - crops herb                                  ',
     &     '16 - crops woody                                 ',
     &     '17 - Permanent snow/ice                          ',
     &     '18 - Bare or sparsely vegetated, urban           '
     &     /)


      !* 1 - evergreen broadleaf early successional
      !* 2 - evergreen broadleaf late successional
      !* 3 - evergreen needleleaf early successional
      !* 4 - evergreen needleleaf late successional
      !* 5 - cold deciduous broadleaf early successional
      !* 6 - cold deciduous broadleaf late successional
      !* 7 - drought deciduous broadleaf
      !* 8 - deciduous needleleaf
      !* 9 - cold adapted shrub
      !* 10 - arid adapted shrub
      !* 11- C3 grass perennial
      !* 12 - C4 grass
      !* 13 - C3 grass - annual
      !* 14- arctic C3 grass
      !* 15- C4 crops
      !* 16 - crops broadleaf woody

      integer, parameter :: EVGRBROADEARLY = 1
      integer, parameter :: EVGRBROADLATE = 2
      integer, parameter :: EVGRNEEDLEEARLY = 3
      integer, parameter :: EVGRNEEDLELATE = 4
      integer, parameter :: COLDDECIDBROADEARLY = 5
      integer, parameter :: COLDDECIDBROADLATE = 6
      integer, parameter :: DROUGHTDECIDBROAD =7
      integer, parameter :: DECIDNEEDLE = 8
      integer, parameter :: COLDSHRUB = 9
      integer, parameter :: TUNDRA = 9 !Hack for same name as Matthews - NK
      integer, parameter :: ARIDSHRUB = 10
      integer, parameter :: GRASSC3PER = 11
      integer, parameter :: GRASSC4 = 12
      integer, parameter :: GRASSC3 = 13
      integer, parameter :: GRASSC3ARCTIC = 14
      integer, parameter :: CROPSC4 = 15
      integer, parameter :: CROPSWOODY = 16
      integer, parameter :: SAND = 17
      integer, parameter :: BDIRT = 18
!##### TEMPORARY HACK - YK #####
!to avoid the conflict in ent_prescribed_drv.f90, using CROPS!
      integer, parameter :: CROPS = 15
!##### END OF TEMPORARY HACK #####


      !* netcdf names *!
      character(len=13), parameter :: ent_cover_names(N_COVERTYPES) = (/
     &     "ever_br_early",
     &     "ever_br_late ",
     &     "ever_nd_early",
     &     "ever_nd_late ",
     &     "cold_br_early",
     &     "cold_br_late ",
     &     "drought_br   ",
     &     "decid_nd     ",
     &     "cold_shrub   ",
     &     "arid_shrub   ",
     &     "c3_grass_per ",
     &     "c4_grass     ",
     &     "c3_grass_ann ",
     &     "c3_grass_arct",
     &     "crops_herb   ",
     &     "crops_woody  ",
     &     "bare_bright  ",
     &     "bare_dark    "
     &     /)

      
      ! other parameters needed for Ent to compile
      integer, parameter :: COVEROFFSET = 0

      !pst - photosynthetic pathway, 1-C3, 2-C4
      !woody - Woody, FALSE=NO, TRUE=YES.
      !leaftype - 1=broadleaf, 2=needleleaf, 3=monocot (not crops)
      !hwilt - wilting point (m)
      !sstar - soil moisture stress onset point (fraction of soil volumetric saturation)
      !swilt - wilting point (fraction of soil volumetric saturation)
      !nf - canopy nitrogen factor (dimensionless)
      !* Parameters from CASA:
      !sla - specific leaf area (m2 leaf area/kg leaf C)(CASA)
      !lrage - leaf and root litter age (years) (CASA)
      !woodage - stem litter age (years) (CASA)
      !lit_C2N - litter C:N (CASA)
      !lignin - lignin content
      !* Parameter for phenology
      !phenotype - phenological types
      !          - evergreen (1), 
      !            cold deciduous (2), 
      !            drought deciduous (3), 
      !            cold/drought deciduous (4),
      !            annual  (5)
      !annual - Annual, FALSE = NO, TRUE=YES.
      !* Parameters for plant allomteries 
      !* (Albani et al. Global Change Biology 2006 & 
      !* estimated from KM67 (Santarem, Amazon) tree survey.)
      !b1Cf - para 1 for allometric relation between DBH & foliage C 
      !b2Cf - para 2 for allometric relation between DBH & foliage C
      !b1Cd - para 1 for allometric relation between DBH & structural(dead) C
      !b2Cd - para 2 for allometric relation between DBH & structural(dead) C
      !b1Ht - para 1 for allometric relation between DBH & height
      !b2Ht - para 2 for allometric relation between DBH & height

      !SOURCES:
      !soil  moisture points:  7 & 10 from Kiang (2002) dissertation.
      !                       Other values are guesses.
      !nf:  Kull&Kruijt ps cap/leaf N param. Guesses from Friend&Kiang (2005) 
      !sla:  from LSM parameterizations.
      !lrage, woodage,lit_C2N,lignin: from CASA parameterizations.
      
      !***************************************************
      !Temp values for Ent pfts (See ent_const.f for types)
      type(pftype),parameter :: pfpar(N_PFT) =         !PFT parameters
      !pst, woody,leaftype, hwilt, sstar, swilt,nf,sla,r,
      !lrage,woodage,lit_C2N,lignin,croot_ratio,phenotype, 
      !b1Cf, b2Cf, b1Cd, b2Cd, b1Ht, b2Ht
     &     (/      
     ! !* 1 - evergreen broadleaf early successional
     &     pftype(1,.true.,1,-153.d0, .60d0, .29d0, 1.2d0, 8.8d0, 0.5d0, 
     &     1.8d0, 41.0d0, 40.d0, 0.2d0, 0.075d0, 1,
     &     0.0347d0, 1.560d0, 0.0816d0, 2.306d0, 34.62d0, -0.0232d0),
     ! !* 2 - evergreen broadleaf late successional
!YK - KM67, Tapajo National Forest, Santarem, Brazil
!YK - lrage=3.0 & sla=9.7
!YK     &     pftype(1,.true.,1,-153.d0, .60d0, .29d0, 1.1d0, 8.8d0, 0.5d0, 
!YK     &     1.8d0,41.0d0, 40.d0, 0.2d0, 0.075d0, 1,
     &     pftype(1,.true.,1,-153.d0, .60d0, .29d0, 1.1d0, 9.7d0, 0.5d0, 
     &     3.0d0,41.0d0, 40.d0, 0.2d0, 0.075d0, 1,
     &     0.0395d0, 1.560d0, 0.1017d0, 2.306d0, 34.62d0, -0.0232d0),
     ! !* 3 - evergreen needleleaf early successional
     &     pftype(1,.true.,2,-153.d0, .50d0, .25d0, 0.9d0, 9.5d0, 1.2d0, 
     &     4.0d0, 42.0d0, 80.d0, 0.25d0, 0.184d0, 1,
     &     0.0240d0, 1.899d0, 0.1470d0, 2.238d0, 27.14d0, -0.0388d0),
     ! !* 4 - evergreen needleleaf late successional
     &     pftype(1,.true.,2,-153.d0, .50d0, .25d0, 0.85d0,9.5d0, 1.2d0, 
     &     4.0d0, 42.0d0, 80.d0, 0.25d0, 0.184d0, 1,
     &     0.0450d0, 1.683d0, 0.1617d0, 2.1536d0, 22.79d0, -0.0445d0),
     ! !* 5 - cold deciduous broadleaf early successional
     &     pftype(1,.true.,1,-500.d0, .50d0, .29d0, 1.5d0,34.5d0,0.6d0, !34.5, 0.6d0
!     &     pftype(1,.true.,1,-500.d0, .50d0, .29d0, 1.5d0, 12.4d0,0.6d0, 
     &     1.2d0, 58.0d0, 57.d0, 0.2d0, 0.093d0, 2,
     &     0.0240d0, 1.860d0, 0.1480d0, 2.411d0, 25.18d0, -0.0496d0),
     ! !* 6 - cold deciduous broadleaf late successional
     &     pftype(1,.true.,1,-500.d0, .50d0, .29d0, 1.4d0,34.0d0,0.6d0, !old SLA 11.5, 0.6
     &     0.75d0, 58.0d0, 57.d0, 0.2d0, 0.093d0, 2,
     &     0.0170d0, 1.731d0, 0.2350d0, 2.252d0, 27.9d0, -0.0540d0), !##YK
!     &     0.0170d0, 1.731d0, 0.2350d0, 2.252d0, 36.8d0, -0.0540d0),   !##NK
     ! !* 7 - drought deciduous broadleaf
!     &     pftype(1,.true.,1,-500.d0, .45d0, .22d0, 1.4d0, 8.3d0, 0.5d0,
!YKIM - Oaks at Tonzi
!YKIM - sstar = 0.34 & swilt 0.28
!     &     pftype(1,.true.,1,-500.d0, .34d0, .28d0, 1.4d0, 8.3d0, 0.5d0, 
!NK - SLA=29 gives right tree density for Tonzi given Cfol_fn
     &     pftype(1,.true.,1,-500.d0, .34d0, .28d0, 1.4d0, 29.d0, 0.5d0, 
     &     1.2d0,25.0d0, 60.d0, 0.2d0, 0.153d0, 3,
     &     0.0296d0, 1.560d0, 0.0621d0, 2.306d0, 34.62d0, -0.0232d0),
     ! !* 8 - deciduous needleleaf !## SLA from Reich (1997) leaf longev. 1 yr
     &     pftype(1,.true.,2,-100.d0,.55d0, .25d0, 0.9d0, 10.0d0, 0.9d0, 
     &     1.8d0, 27.0d0, 50.d0, 0.2d0, 0.2d0, 2,
     &     0.0240d0, 1.899d0, 0.1470d0, 2.238d0, 27.14d0, -0.0388d0),
     ! !* 9 - cold adapted shrub
     &     pftype(1,.true.,1,-153.d0,.50d0, .30d0, 1.4d0, 2.25d0, 0.6d0, 
     &     2.8d0, 5.5d0, 50.d0, 0.15d0, 1.40d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 10 - arid adapted shrub
     &     pftype(1,.true.,1,-2030.d0,.40d0,.22d0, 1.3d0, 3.25d0, 0.6d0, 
     &     1.0d0, 5.5d0, 65.d0, 0.2d0, 0.32d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 11 - C3 grass perennial
     &     pftype(1,.false.,3,-2030.d0,.30d0,.10d0,1.5d0, 21.6d0, 1.2d0, 
     &     1.5d0, UNDEF, 50.d0, 0.1d0, 0.0d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 12 - C4 grass
     &     pftype(2,.false.,3,-2030.d0,.65d0,.27d0, 1.3d0,9.4d0, 0.6d0, 
     &     1.5d0, UNDEF, 50.d0, 0.1d0, 0.0d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 13 - C3 grass - annual 
!YKIM - Grasses at Tonzi
!YKIM - sstar = 0.65 & swilt 0.27
     &     pftype(1,.false.,3,-2030.d0,.65d0,.27d0,1.5d0, 21.6d0, 1.2d0, !10->15  
     &     1.5d0, UNDEF, 50.d0, 0.1d0, 0.0d0, 5,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 14 - arctic C3 grass
     &     pftype(1,.false.,3,-153.d0,.60d0, .27d0, 1.4d0, 9.0d0, 0.6d0, 
     &     1.5d0, UNDEF, 50.d0, 0.1d0, 0.0d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 15 - C4 crops herbaceous
     &     pftype(2,.false., 1,-153.d0,.45d0,.27d0,1.3d0, 15.d0, 0.6d0, 
     &     1.1d0, UNDEF, 52.5d0, 0.16d0, 0.0d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     ! !* 16 - crops - broadleaf woody !## COPIED FROM BROAD COLDDECID LATE ##
     &     pftype(1,.true.,1,-153.d0,.50d0, .29d0, 1.4d0, 8.3d0, 0.9d0, 
     &     1.2d0, 58.0d0, 50.d0, 0.2d0, 0.093d0, 4,
     &     0.0170d0, 1.731d0, 0.2350d0, 2.252d0, 23.39d0, -0.0540d0)
     &     /)

      !*********************************************************
      !* Prescribed albedoes for Matthews prescribed seasonality
      !*********************************************************
      real*8, parameter :: ALBVND(N_COVERTYPES,4,6) = RESHAPE( (/
C
!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-updated: 
! decid_needle=decid fall-winter & needle spring-summer, 
! c3grass_ann =c3grass opposite seasons
!, cropstree = cold_broad
C
C     (1)  >SRBALB(6) = VIS  (300 - 770 nm)
     1 .061,.061,.067,.067,.100,.100,.078, .100,.067,.089,.089,.089,
     &     .091, .089, .089, .100, .500,.000,
     2 .061,.061,.067,.067,.055,.055,.073, .067,.100,.100,.100,.100,
     &     .089, .100, .100, .055, .500,.000,
     3 .061,.061,.083,.083,.058,.058,.085, .083,.139,.139,.091,.091,
     &     .100, .091, .091, .058, .500,.000,
     4 .061,.061,.061,.061,.055,.055,.064, .061, .111,.111,.090,.090,
     &     .100,.090, .090, .055, .500,.000,
C
C     (2)  >SRBALB(5) = NIR  (770 - 860 nm)    (ANIR=Ref)
     1 .183,.183,.200,.200,.300,.300,.233,.300,.200,.267,.267,.267,
     &     .350,.267,.267,.300,.500,.000,
     2 .183,.183,.200,.200,.218,.218,.241,.200,.300,.300,.350,.350,
     &     .364,.350,.350,.218,.500,.000,
     3 .183,.183,.250,.250,.288,.288,.297,.250,.417,.417,.364,.364,
     &     .267,.364,.364,.288,.500,.000,
     4 .183,.183,.183,.183,.218,.218,.204,.183,.333,.333,.315,.315,
     &     .315,.315,.315,.218,.500,.000,
C
C     (3)  >SRBALB(4) = NIR  (860 -1250 nm)    (ANIR*1.0)
     1 .183,.183,.200,.200,.300,.300,.233,.300,.200,.267,.267,.267,
     &     .350,.267,.267,.300,.500,.000,
     2 .183,.183,.200,.200,.218,.218,.241,.200,.300,.300,.350,.350,
     &     .364,.350,.350,.218,.500,.000,
     3 .183,.183,.250,.250,.288,.288,.297,.250,.417,.417,.364,.364,
     &     .267,.364,.364,.288,.500,.000,
     4 .183,.183,.183,.183,.218,.218,.204,.218,.333,.333,.315,.315,
     &     .315,.315,.315,.218,.500,.000,
C
C     (4)  >SRBALB(3) = NIR  (1250-1500 nm)    (ANIR*0.4)
     1 .073,.073,.080,.080,.120,.120,.093,.120,.080,.107,.107,.107,
     &     .140,.107,.107,.120,.500,.000,
     2 .073,.073,.080,.080,.083,.083,.096,.080,.120,.120,.140,.140,
     &     .145,.140,.140,.083,.500,.000,
     3 .073,.073,.100,.100,.115,.115,.119,.100,.167,.167,.145,.145,
     &     .107,.145,.145,.115,.500,.000,
     4 .073,.073,.073,.073,.087,.087,.081,.087,.132,.132,.126,.126,
     &     .126,.126,.126,.087,.500,.000,
C
C     (5)  >SRBALB(2) = NIR  (1500-2200 nm)    (ANIR*0.5)
     1 .091,.091,.100,.100,.150,.150,.116,.150,.100,.133,.133,.133,
     &     .175,.133,.133,.150,.500,.000,
     2 .091,.091,.100,.100,.109,.109,.120,.100,.150,.150,.175,.175,
     &     .182,.175,.175,.109,.500,.000,
     3 .091,.091,.125,.125,.144,.144,.148,.125,.208,.208,.182,.182,
     &     .133,.182,.182,.144,.500,.000,
     4 .091,.091,.091,.091,.109,.109,.102,.109,.166,.166,.157,.157,
     &     .157,.157,.157,.109,.500,.000,
C
C     (6)  >SRBALB(1) = NIR  (2200-4000 nm)    (ANIR*0.1)

     1 .018,.018,.020,.020,.030,.030,.023,.030,.020,.027,.027,.027,
     &     .035,.027,.027,.030,.500,.000,
     2 .018,.018,.020,.020,.022,.022,.024,.020,.030,.030,.035,.035,
     &     .036,.035,.035,.022,.500,.000,
     3 .018,.018,.025,.025,.029,.029,.030,.025,.042,.042,.036,.036,
     &     .027,.036,.036,.029,.500,.000,
     4 .018,.018,.018,.018,.022,.022,.020,.022,.033,.033,.032,.032,
     &     .032,.032,.032,.022,.500,.000
     *     /),(/N_COVERTYPES,4,6/) )

      real*8, parameter :: rhol(N_PFT,N_BANDS) = RESHAPE( (/
     1   0.07, 0.07, 0.07, 0.10, 0.10, 0.10,
     &   0.10, 0.10, 0.07, 0.10, 0.10, 0.11, 0.11, 0.11, 0.11, 0.11,
     2   0.35, 0.35, 0.35, 0.45, 0.45, 0.45,
     &   0.45, 0.45, 0.35, 0.45, 0.45, 0.58, 0.58, 0.58, 0.58, 0.58,
     3   0.35, 0.35, 0.35, 0.45, 0.45, 0.45,
     &   0.45, 0.45, 0.35, 0.45, 0.45, 0.58, 0.58, 0.58, 0.58, 0.58,
     4   0.35, 0.35, 0.35, 0.45, 0.45, 0.45,
     &   0.45, 0.45, 0.35, 0.45, 0.45, 0.58, 0.58, 0.58, 0.58, 0.58,
     5   0.35, 0.35, 0.35, 0.45, 0.45, 0.45,
     &   0.45, 0.45, 0.35, 0.45, 0.45, 0.58, 0.58, 0.58, 0.58, 0.58,
     6   0.35, 0.35, 0.35, 0.45, 0.45, 0.55,
     &   0.45, 0.45, 0.35, 0.45, 0.45, 0.58, 0.58, 0.58, 0.58, 0.58
     *     /), (/N_PFT, N_BANDS/) )

      real*8, parameter :: rhos(N_PFT,N_BANDS) = RESHAPE( (/
     1    0.16, 0.16, 0.16, 0.16, 0.16, 0.16,
     &    0.16, 0.16, 0.16, 0.16, 0.16, 0.36, 0.36, 0.36, 0.36, 0.36,
     2    0.39, 0.39, 0.39, 0.39, 0.39, 0.39,
     &    0.39, 0.39, 0.39, 0.39, 0.39, 0.58, 0.58, 0.58, 0.58, 0.58,
     3    0.39, 0.39, 0.39, 0.39, 0.39, 0.39,
     &    0.39, 0.39, 0.39, 0.39, 0.39, 0.58, 0.58, 0.58, 0.58, 0.58,
     4    0.39, 0.39, 0.39, 0.39, 0.39, 0.39,
     &    0.39, 0.39, 0.39, 0.39, 0.39, 0.58, 0.58, 0.58, 0.58, 0.58,
     5    0.39, 0.39, 0.39, 0.39, 0.39, 0.39,
     &    0.39, 0.39, 0.39, 0.39, 0.39, 0.58, 0.58, 0.58, 0.58, 0.58,
     6    0.39, 0.39, 0.39, 0.39, 0.39, 0.39,
     &    0.39, 0.39, 0.39, 0.39, 0.39, 0.58, 0.58, 0.58, 0.58, 0.58
     *     /), (/N_PFT, N_BANDS/) )

      real*8, parameter :: taul(N_PFT,N_BANDS) = RESHAPE( (/
     1   0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
     &   0.05, 0.05, 0.05, 0.05, 0.05, 0.07, 0.07, 0.07, 0.07, 0.07,
     2   0.10, 0.10, 0.10, 0.25, 0.25, 0.25,
     &   0.25, 0.25, 0.10, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
     3   0.10, 0.10, 0.10, 0.25, 0.25, 0.25,
     &   0.25, 0.25, 0.10, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
     4   0.10, 0.10, 0.10, 0.25, 0.25, 0.25,
     &   0.25, 0.25, 0.10, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
     5   0.10, 0.10, 0.10, 0.25, 0.25, 0.25,
     &   0.25, 0.25, 0.10, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
     6   0.10, 0.10, 0.10, 0.25, 0.25, 0.25,
     &   0.25, 0.25, 0.10, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25
     *     /), (/N_PFT, N_BANDS/) )

      real*8, parameter :: taus(N_PFT,N_BANDS) = RESHAPE( (/
     1   0.001, 0.001, 0.001, 0.001, 0.001,
     &   0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.220, 0.220,
     &   0.220, 0.220, 0.220,
     2   0.001, 0.001, 0.001, 0.001, 0.001,
     &   0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.380, 0.380,
     &   0.380, 0.380, 0.380,
     3   0.001, 0.001, 0.001, 0.001, 0.001,
     &   0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.380, 0.380,
     &   0.380, 0.380, 0.380,
     4   0.001, 0.001, 0.001, 0.001, 0.001,
     &   0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.380, 0.380,
     &   0.380, 0.380, 0.380,
     5   0.001, 0.001, 0.001, 0.001, 0.001,
     &   0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.380, 0.380,
     &   0.380, 0.380, 0.380,
     6   0.001, 0.001, 0.001, 0.001, 0.001,
     &   0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.380, 0.380,
     &   0.380, 0.380, 0.380
     *     /), (/N_PFT, N_BANDS/) )

      !***************************************************
      !* ecophys const - leaf/stem orientation index
      !***************************************************
      real, parameter :: xl(N_PFT) =
     &    (/0.01, 0.01, 0.01, 0.10, 0.10, 0.01, 0.25, 0.25,
     &    0.01, 0.25, 0.25, -0.30, -0.30, -0.30, -0.30, -0.30/)

      !***************************************************
      !* PFT categories
      !***************************************************

      logical, parameter :: is_crop(N_PFT) =
     &     (/ .false., .false., .false., .false.,
     &     .false., .false., .false., .false.,
     &     .false., .false., .false., .false.,
     &     .false., .false., .true., .true. /)

      logical, parameter :: is_hw(N_PFT) = 
     &     (/ .true., .true., .false., .false.,
     &     .true., .true., .true., .false.,
     &     .false., .false., .false., .false.,
     &     .false., .false., .false., .true. /)

      logical, parameter :: is_conifer(N_PFT) =
     &     (/ .false., .false., .true., .true.,
     &     .false., .false., .false., .true.,
     &     .false., .false., .false., .false.,
     &     .false., .false., .false., .false. /)

      logical, parameter :: is_grass(N_PFT) =
     &     (/ .false., .false., .false., .false., 
     &     .false., .false., .false., .false.,
     &     .false., .false., .true., .true.,
     &     .true., .true., .false., .false. /)


      !***************************************************
      !* Prescribed max and min LAI
      !***************************************************
!#ifdef FLUXNETINIT  - No need for this flag here any more - NK
!      !* FLUXNET LAI *!
      real*8, parameter :: alamax(N_COVERTYPES) =
     $     (/ 6.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0 ,6.0d0, 4.0d0
     &     ,6.d0, 1.5d0, 2.5d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0
     &     ,4.5d0, 6.0d0, 0.d0, 0.d0/)
      real*8, parameter :: alamin(N_COVERTYPES) =
     $     (/ 5.0d0, 4.0d0, 6.0d0, 2.25d0, 1.0d0, 1.0d0,1.0d0
     &     ,1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.1d0, 1.0d0
     &     ,1.0d0, 1.0d0, 0.d0, 0.d0 /)

            integer, parameter :: laday(N_COVERTYPES) =
     $     (/ 196, 196, 196, 196, 196, 196, 196
     &     ,196, 196, 196, 196, 196, 196, 196
     &     , 196, 196, 0, 0 /)

      !***************************************************
      !* Rosenzweig & Abramopoulos root profile parameters
      !***************************************************
 !--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-updated
      real*8, parameter :: aroot(N_COVERTYPES) = 
     $     (/ 1.1d0, 1.1d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.25d0
     &     ,0.25d0, 0.8d0, 0.8d0, 0.9d0, 0.9d0, 0.9d0, 0.9d0
     &     ,0.9d0, 0.25d0, 0.d0, 0.d0 /)
      real*8, parameter :: broot(N_COVERTYPES) = 
     $     (/ 0.4d0, 0.4d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0, 2.0d0
     &     ,2.0d0, 0.4d0, 0.4d0, 0.9d0, 0.9d0, 0.9d0, 0.9d0
     &     , 0.9d0, 2.0d0, 0.0d0, 0.0d0 /)


      !***************************************************
      !* Prescribed vegetation height (m)
      !***************************************************
      real*8, parameter :: vhght(N_COVERTYPES) =
!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-updated
     $     (/25d0, 25d0, 30d0, 30d0, 7.1d0, 7.1d0, 5d0
     &     , 30d0, 0.1d0, 5d0, 1.5d0, 1.5d0, 1.5d0, 1.5d0
     &     ,1.75d0, 7.1d0, 0.d0, 0.d0 /)

       !***************************************************
      !*  Mean canopy nitrogen (nmv; g/m2[leaf])
      !***************************************************
      real*8, parameter :: nmv(N_COVERTYPES) =
!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-updated
     $     (/2.7d0, 2.7d0, 2.9d0, 2.9d0, 1.25d0, 1.25d0, 1.03d0
     &     , 2.9d0, 1.6d0, 2.38d0, 3.27d0, 0.82d0, 3.27d0
     &     , 3.27d0, 2.50d0, 1.25d0, 0.d0, 0.d0 /)

      !***************************************************
      !* Soil color by cover type (GISS GCM)
      !***************************************************
!--- ever_ES_broad ever_LS_broad ever_ES_needle ever_LS_needle 
!----cold_ES_broad cold_LS_broad drought_broad decid_needle shrub_cold 
!----shrub_arid c3grass c4grass c3grass_ann c3grass_arctic 
!----cropsc4 cropstree
!----sand bdirt
!KIM - temp. values, NK-updated
      integer, parameter :: soil_color_prescribed(N_COVERTYPES) =
     $     (/ 2, 2,  2, 2, 2, 2, 2
     &     , 2, 2, 2, 2, 2 ,2, 2, 2, 2, 1, 2 /)


      !***************************************************
      !* Soil carbon by cover type (GISS GCM) 
      !* (fractions of total soil carbon in layer) - 1 layer only, modeled
      !***************************************************
!YK - temp. values, modified from 8 GISS pfts below
!      real*8,parameter :: Cpool_fracs(N_PFT,NPOOLS-NLIVE,N_CASA_LAYERS)=
      real*8,parameter :: Cpool_fracs1(N_PFT,NPOOLS-NLIVE,1)=
     &     RESHAPE( (/
        !1. ever_ES_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349,
        !2. ever_LS_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349,
        !3. ever_ES_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349,
        !4. ever_LS_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349,
        !5.  cold_ES_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349,
        !6.  cold_LS_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349,
        !7.  drought_broad
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989,
        !8.  decid_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349,
        !9.  shrub_cold
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822,
        !10.  shrub_arid
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989,
        !11.  c3grass
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822,
        !12.  c4grass
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822,
        !13.  c3grass_ann
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822,
        !14.  c3grass_arctic
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822,
        !15.  cropsc4
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822,
        !16.  cropstree
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /),
     &     ( /N_PFT,NPOOLS-NLIVE,1/ ))
!     &     ( /N_PFT,NPOOLS-NLIVE,N_CASA_LAYERS/ )
      
      !*********************************************************************

      !wdens_g_cm3 is not calculated from wooddensity_gcm3 but is from data.
      real*8, DIMENSION(N_PFT), parameter :: wdens_g_cm3 =
     &     (/ 0.66d0,0.7d0,0.5d0,0.5d0,0.54d0,0.54d0,0.6d0
     &     ,0.54d0,0.6d0,0.6d0
     &     ,undef,undef,undef,undef,undef,0.54d0 /)
      real*8, DIMENSION(N_PFT), parameter :: a0h = !Later move to ent_pfts_ENT.f
     &     (/ 1.3d0,0.d0,1.3d0,0.d0,1.3d0,1.3d0,0.d0,1.3d0,0.d0,0.d0
     &     ,undef,undef,undef,undef,undef,0.d0 /)
      real*8, DIMENSION(N_PFT), parameter :: acr =
     &     (/ 0.1407d0,0.1407d0,0.2855d0,0.2570d0,0.3070d0,0.2773d0
     &     ,0.3868d0,0.3046d0,0.500d0,0.500d0
     &     ,undef,undef,undef,undef,undef,0.3237d0 /)
      real*8, DIMENSION(N_PFT), parameter :: bcr =
     &     (/ 1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0,1.d0
     &     ,undef,undef,undef,undef,undef ,1.d0 /)
      real*8, DIMENSION(N_PFT), parameter :: bR =
     &     (/ 1.50d0,1.50d0,2.00d0,2.00d0,1.00d0,1.00d0,0.75d0
     &     ,2.0d0,1.0d0,1.0d0,40.0d0,40.0d0,40.0d0,40.0d0,10.0d0,1.d0 /)
      integer, DIMENSION(N_PFT), parameter  :: form =
     &     (/ TREE,TREE,TREE,TREE,TREE,TREE,TREE,TREE,SHRUB,SHRUB
     &     ,HERB,HERB,HERB,HERB,HERB,TREE /)
      real*8, DIMENSION(N_PFT), parameter  :: DBHBAmax_cm =
     &     (/ 150.d0,150.d0,150.d0,150.d0,150.d0,150.d0,150.d0,150.d0
     &     ,10.d0,10.d0,undef,undef,undef,undef,undef,50.d0 /)
      logical, DIMENSION(N_PFT), parameter  :: crop =
     &     (/ .false.,.false.,.false.,.false.,.false.,.false.
     &     ,.false.,.false.,.false.,.false.,.false.,.false.
     &     ,.false.,.false.,.true.,.true. /)

!**********************************************************************
      end module ent_pfts
