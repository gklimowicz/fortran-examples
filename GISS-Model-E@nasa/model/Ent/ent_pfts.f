      module ent_pfts
!@sum Parameter sets for GISS 8 vegetation biome types for Ent model.

      use ent_const
      use ent_types
      implicit none


      !***************************************************
      !*            GISS VEGETATION TYPES                *
      !***************************************************
      !Types: These are GISS GCM types until we get full data sets for Ent.
      !1-tundra, 2-grassland, 3-shrubland, 4-savanna, 5-deciduous forest,
      !6-evergreen needleleaf forest, 7-tropical rainforest, 8-crops
      !**** Correspondence to CASA/LSM vegetation types for sla from 
      !**** Dickinson et al. (J. Clim., Nov. 1998).  See casa_bgfluxes.f and
      !**** 
      !#-GISS   = Dickinson = CASA/LSM
      !1-tundra = tundra and semidesert
      !2-grassland = avg short grass & tall grass = LSM cool C3, warm C4 grasses
      !3-shrubland = evergreen shrubs, deciduous shrubs
      !4-savanna = shrubs & broadleaf evergreen & grasses (avg 25,35,40)
      !5-deciduous forest = deciduous broadleaf & decid needleleaf
      !6-evergreen needleleaf = needleleaf evergreen
      !7-tropical rainforest = deciduous broadleaf = LSM tropical seasonal tree 
      !8-crops = crops

      integer, parameter :: COVEROFFSET = 1 !SAND in first position in GISS array
      !PFT number
      integer, parameter :: TUNDRA = 1
      integer, parameter :: GRASSC3 = 2
      integer, parameter :: SHRUBGRASS = 3
      integer, parameter :: SAVANNA = 4
      integer, parameter :: DECIDFOREST = 5
      integer, parameter :: EVERGRNEEDLE = 6
      integer, parameter :: TROPRAINF = 7
      integer, parameter :: CROPS = 8
      !Array position
      integer, parameter :: SAND = 1
      integer, parameter :: BDIRT = 10

!##### TEMPORARY HACK - YK #####
!to avoid the conflict in phenology.f
      integer, parameter :: DROUGHTDECIDBROAD = 20
      integer, parameter :: GRASSC3ARCTIC = 21

      character(len=15), parameter ::
     &     ent_cover_names(N_COVERTYPES-N_OTHER) = (/
     &     'brightsoil     ','tundra         ','grass          ',
     &     'shrub_and_grass','tree_and_grass ','deciduous      ',
     &     'evergreen      ','rainforest     ','cultivation    ',
     &     'darksoil       '
     &     /)

!      !* netcdf names *!
!      character(len=13), parameter :: ent_cover_names(N_COVERTYPES) = (/
!     &     "BSAND        ",
!     &     "TNDRA        ",
!     &     "GRASS        ",
!     &     "SHRUBGRASS   ",
!     &     "SAVANNA      ",
!     &     "DECIDFOREST  ",
!     &     "EVERGRNEEDLE ",
!     &     "TROPRAINF    ",
!     &     "CROPS        ",
!     &     "BDIRT        ",
!     &     "ALGAE        ",
!     &     "GRAC4        "
!     &     /)


!##### TEMPORARY HACK - YK #####
      !*-----------------------------------------
      !* Veg types correspondence between models:
*       LSM:   1  2    3  4  5  6  7  8  9 10 11 12 13 14^M
*      CASA:   4  5    1  2  6  7  9 11 10 10 12 12  6  8^M
*      GISS:   6  5,6  7  5  4  2  3 x  1  1  8  8   4  x
*     NOTE:  Need to distinguish: - needleleaf decid vs. evergreen
*                                 - bare soil
*                                 - savanna mix
*            Check for:  sla, lrage, woodage, albedo
      !*-----------------------------------------
      !For dependencies to work, these constants must go in ent_const.f
!      integer,parameter :: N_PFT = 8
!      integer,parameter :: N_SOILCOV = 2 !light sand, dark dirt (GISS)
!      integer,parameter :: N_OTHER = 1
!      integer,parameter :: N_COVERTYPES = N_PFT + N_SOILCOV + N_OTHER
       !pst - photosynthetic pathway, 1-C3, 2-C4
      !hwilt - wilting point (m)
      !sstar - soil moisture stress onset point (fraction of soil volumetric saturation)
      !swilt - wilting point (fraction of soil volumetric saturation)
      !nf - canopy nitrogen factor (dimensionless)
      !* Parameters from CASA:
      !sla - specific leaf area (m2 leaf area/kg leaf C)(CASA)
      !r - CASA respiration parameter (gC/gN)(CASA)
      !lrage - leaf and root litter turnover time (years) (CASA)
      !woodage - stem litter turnover time (years) (CASA)
      !lit_C2N - litter C:N (CASA)
      !lignin - lignin content
      !*for GISS pft, 
      !*phenotype, b1Cf, b2Cf, b1Cd, b2Cd, b1Ht, b2Ht are dummies!
      !phenotype - phenological type 
      !* Parameters for plant allomteries 
      !* (Albani et al. Global Change Biology 2006 & 
      !* estimated from KM67 (Santarem, Amazon) tree survey.)
      !b1Cf - para 1 for allometric relation between DBH & foliage C 
      !b2Cf - para 2 for allometric relation between DBH & foliage C
      !b1Cd - para 1 for allometric relation between DBH & structural(dead) C
      !b2Cd - para 2 for allometric relation between DBH & structural(dead) C
      !b1Ht - para 1 for allometric relation between DBH & height
      !b2Ht - para 2 for allometric relation between DBH & height

      type(pftype),parameter :: pfpar(N_PFT) =          !PFT parameters
     &!        pst,woody,leaftype,hwilt,sstar,swilt,nf,sla,
     &!        r,lrage,woodage,lit_C2N,lignin,croot_ratio,phenotype ! 
     &!        b1Cf, b2Cf, b1Cd, b2Cd, b1Ht, b2Ht
     &     (/
     &     pftype(1,.true., 1,-153.d0,  .50d0, .30d0,  1.4d0, !tundra
     &     2.25d0, 0.6d0, 2.8d0, 5.5d0, 50.0d0,0.15d0,1.4d0,2,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     &     pftype(1,.false., 3,-2030.d0,  .30d0, .10d0,  1.5d0, !grassC3 !hwilt=Vaira grassland final senescence soilmp layer2
     &     11.7d0, 1.2d0, 1.5d0, UNDEF, 50.0d0, 0.1d0,0.d0,4, ! Vaira SLA from adjustment for C_labile.
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     &     pftype(1,.true., 1,-153.d0,  .40d0, .22d0,  1.3d0, !shrub
     &     3.25d0, 0.6d0, 1.0d0, 5.5d0, 57.5d0, 0.15d0,0.32d0,3,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0),
     &     pftype(1,.true., 1,-2030.d0,  .65d0, .22d0,  1.3d0, !savanna
     &     5.1d0, 1.d0, 1.8d0, 25.d0, 50.0d0, 0.15d0,0.153d0,3,
     &     0.0296d0, 1.560d0, 0.0621d0, 2.306d0, 27.d0, -0.0232d0),
     &     pftype(1,.true., 1,-500.d0,  .50d0, .29d0,  1.5d0, !decidforest
!     &     8.3d0, 0.6d0, 1.2d0, 58.0d0, 50.0d0, 0.2d0,0.093,2, !SLA for Quercus ilex, Mediavilla & Escudero(2003)
     &     34.5d0, 0.6d0, 1.2d0, 58.0d0, 57.0d0, 0.2d0,0.093,2, !SLA for oak, Tatarinov & Cienciala (2006)
     &     0.0170d0, 1.731d0, 0.2350d0, 2.252d0, 23.39d0, -0.0540d0), !late successional
     &     pftype(1,.true., 2,-153.d0,  .50d0, .25d0,  0.9d0, !evergrneedle
     &     5.9d0, 1.2d0, 4.d0,42.0d0, 80.0d0,0.25d0,0.184d0,1, !SLA for Pinus sylvestris, Pensa and Sellin (2002).SLA-CLM 10. lrage WAS 5.0!!!! Pinus sylvestris 2-4 years!!!-NYK
!     &     0.0240d0, 1.899d0, 0.1470d0, 2.238d0, 27.14d0, -0.0388d0), !early succ
     &     0.0450d0, 1.683d0, 0.1617d0, 2.1536d0, 22.79d0, -0.0445d0), !late succ
     &     pftype(1,.true., 1, -153.d0,  .60d0, .29d0,  1.1d0, !troprainf
     &     9.9d0, 0.5d0, 1.8d0, 41.0d0, 40.0d0, 0.2d0,0.075d0,1,
     &     0.0347d0, 1.560d0, 0.0816d0, 2.306d0, 34.62d0, -0.0232d0),
     &     pftype(1,.false., 1,-153.d0,  .45d0, .27d0,  1.3d0, !crops
     &     6.36d0, 0.6d0, 1.1d0, UNDEF, 52.5d0, 0.16d0,0.0d0, 4,
     &     0.0800d0, 1.000d0, 0.00001d0, 1.000d0, 0.4778d0, -0.75d0)
!     &     pftype(1,   -100.d0,  .50d0, .30d0,  0.76d0)&
!     &     0.0d0,0.0d0,0.0d0,0.0d0)
     &     /)

      !**NOTE:  Above, sstar and swilt are guesses for all except grassland
      ! and savanna.

      !*********************************************************
      !* Prescribed albedoes for Matthews prescribed seasonality
      !*********************************************************
      real*8, parameter :: ALBVND(N_COVERTYPES,4,6) = RESHAPE( (/
C     (1)  >SRBALB(6) = VIS  (300 - 770 nm)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.067,.089,.089,.078,.100,.067,.061,.089,.000,.200,.089,
     2 .500,.062,.100,.100,.073,.055,.067,.061,.100,.000,.200,.100,
     3 .500,.085,.091,.139,.085,.058,.083,.061,.091,.000,.200,.091,
     4 .500,.080,.090,.111,.064,.055,.061,.061,.090,.000,.200,.090,
C
C     (2)  >SRBALB(5) = NIR  (770 - 860 nm)    (ANIR=Ref)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.200,.267,.267,.233,.300,.200,.183,.267,.000,.200,.267,
     2 .500,.206,.350,.300,.241,.218,.200,.183,.350,.000,.200,.350,
     3 .500,.297,.364,.417,.297,.288,.250,.183,.364,.000,.200,.364,
     4 .500,.255,.315,.333,.204,.218,.183,.183,.315,.000,.200,.315,
C
C     (3)  >SRBALB(4) = NIR  (860 -1250 nm)    (ANIR*1.0)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.200,.267,.267,.233,.300,.200,.183,.267,.000,.200,.267,
     2 .500,.206,.350,.300,.241,.218,.200,.183,.350,.000,.200,.350,
     3 .500,.297,.364,.417,.297,.288,.250,.183,.364,.000,.200,.364,
     4 .500,.255,.315,.333,.204,.218,.183,.183,.315,.000,.200,.315,
C
C     (4)  >SRBALB(3) = NIR  (1250-1500 nm)    (ANIR*0.4)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.080,.107,.107,.093,.120,.080,.073,.107,.000,.200,.107,
     2 .500,.082,.140,.120,.096,.083,.080,.073,.140,.000,.200,.140,
     3 .500,.119,.145,.167,.119,.115,.100,.073,.145,.000,.200,.145,
     4 .500,.102,.126,.132,.081,.087,.073,.073,.126,.000,.200,.126,
C
C     (5)  >SRBALB(2) = NIR  (1500-2200 nm)    (ANIR*0.5)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.100,.133,.133,.116,.150,.100,.091,.133,.000,.200,.133,
     2 .500,.103,.175,.150,.120,.109,.100,.091,.175,.000,.200,.175,
     3 .500,.148,.182,.208,.148,.144,.125,.091,.182,.000,.200,.182,
     4 .500,.127,.157,.166,.102,.109,.091,.091,.157,.000,.200,.157,
C
C     (6)  >SRBALB(1) = NIR  (2200-4000 nm)    (ANIR*0.1)
C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
     1 .500,.020,.027,.027,.023,.030,.020,.018,.027,.000,.200,.027,
     2 .500,.021,.035,.030,.024,.022,.020,.018,.035,.000,.200,.035,
     3 .500,.030,.036,.042,.030,.029,.025,.018,.036,.000,.200,.036,
     4 .500,.026,.032,.033,.020,.022,.018,.018,.032,.000,.200,.032
     *     /),(/N_COVERTYPES,4,6/) )

      !***************************************************
      !* PFT categories
      !***************************************************

      logical, parameter :: is_crop(N_PFT) =
     &     (/ .false., .false., .false., .false.,
     &     .false., .false., .false., .true. /)

      logical, parameter :: is_hw(N_PFT) = 
     &     (/ .false., .false., .false., .true.,
     &     .true., .false., .true., .false. /)

      logical, parameter :: is_conifer(N_PFT) =
     &     (/ .false., .false., .false., .false.,
     &     .false., .true., .false., .false. /)

      logical, parameter :: is_grass(N_PFT) =
     &     (/ .false., .true., .false., .false.,
     &     .false., .false., .false., .false. /)


      !***************************************************
      !* Prescribed max and min LAI
      !***************************************************
!--- sand tundr grass shrub trees decid evrgr rainf crops bdirt algae c4grass
      real*8, parameter :: alamax(N_COVERTYPES) =
      !* Matthews LAI *!
!     $     (/ 0.d0, 1.5d0, 2.0d0, 2.5d0, 4.0d0, 6.0d0,10.0d0,8.0d0,4.5d0
!     &     ,0.d0, 0.d0, 2.d0 /)
      !* Revised Matthews LAI *!
     $     (/ 0.d0, 1.5d0, 2.0d0, 2.5d0, 4.0d0, 6.0d0,8.0d0,7.0d0,3.0d0
     &     ,0.d0, 0.d0, 2.d0 /)

      real*8, parameter :: alamin(N_COVERTYPES) =
      !* Matthews LAI *!
!     $     (/ 0.d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 8.0d0,6.0d0,1.0d0
!     &     ,0.d0, 0.d0, 1.d0 /)
      !* Revised Matthews LAI *!
     $     (/ 0.d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 6.0d0,6.0d0,1.0d0
     &     ,0.d0, 0.d0, 1.d0 /)

      integer, parameter :: laday(N_COVERTYPES) =
     $     (/ 0, 196,  196,  196,  196,  196,  196,  196,  196
     &     ,0, 0, 196 /)

      !***************************************************
      !* Rosenzweig & Abramopoulos root profile parameters
      !***************************************************
 !--- sand tundr grass shrub trees decid evrgr rainf crops bdirt algae c4grass
      real*8, parameter :: aroot(N_COVERTYPES) = 
     $     (/ 0.d0,12.5d0, 0.9d0, 0.8d0,0.25d0,0.25d0,0.25d0,1.1d0,0.9d0
     &     ,0.d0, 0.d0, 0.9d0 /)
      real*8, parameter :: broot(N_COVERTYPES) = 
     $     (/ 0.d0, 1.0d0, 0.9d0, 0.4d0,2.00d0,2.00d0,2.00d0,0.4d0,0.9d0
     &     ,0.d0, 0.d0, 0.9d0 /)

      !***************************************************
      !* Prescribed vegetation height
      !***************************************************
      real*8, parameter :: vhght(N_COVERTYPES) =
!* bsand tundrv grass shrub trees  decid evrgr rainf crops bdirt algae c4grass
      !GISS ORIGINAL
!     $     (/0.d0, 0.1d0, 1.5d0,   5d0,  15d0,  20d0,  30d0, 25d0,1.75d0
!     &     ,0.d0, 0.d0, 1.5d0 /)
      !NYK ADJUSTED VALUES
     $     (/0.d0, 0.1d0, 1.5d0,  5d0,  7.1d0,  20d0,  30d0, 25d0,1.75d0
     &     ,0.d0, 0.d0, 1.5d0 /)

      !***************************************************
      !*  Mean canopy nitrogen (nmv; g/m2[leaf])
      !***************************************************
      real*8, parameter :: nmv(N_COVERTYPES) =
     $     (/0.d0,1.6d0,3.27d0,2.38d0,1.03d0,1.25d0,2.9d0,2.7d0,2.50d0
     &     ,0.d0, 0.d0, 0.82d0 /)

      !***************************************************
      !* Soil color by cover type (GISS GCM)
      !***************************************************
!* bsand tundrv grass shrub trees  decid evrgr rainf crops bdirt algae c4grass
      integer, parameter :: soil_color_prescribed(N_COVERTYPES) =
     $     (/1, 2, 2,  2, 2, 2, 2, 2
     &     ,2, 2, 2, 2 /)


      !***************************************************
      !* Soil carbon by cover type (GISS GCM) 
      !* (fractions of total soil carbon in layer) - 1 layer only, modeled
      !***************************************************
!      real*8,parameter :: Cpool_fracs(N_PFT,NPOOLS-NLIVE,N_CASA_LAYERS)=
      real*8,parameter :: Cpool_fracs1(N_PFT,NPOOLS-NLIVE,1)=
     &     RESHAPE( (/
       !1.  tundra (for now=C3 grass)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822,
       !2.  C3 grass (Vaira)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822,
       !3.  shrub (for now=savanna)
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989,
       !4.  savanna (Tonzi)
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989,
       !5.  decid broadl (MMSF)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349,
       !6.  evergr needl (for now=decid broadl)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349,
       !7.  trop rainf (for now=decid broadl)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349,
       !8.  crops (for now=C3 grass)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /),
     &     ( /N_PFT,NPOOLS-NLIVE,1/ ) )
!     &     ( /N_PFT,NPOOLS-NLIVE,N_CASA_LAYERS/ )



       !***************************************************

C        1    2    3    4    5    6    7    8    9   10   11    12
C      BSAND     GRASS     TREES     EVERG     CROPS     ALGAE
C           TNDRA     SHRUB     DECID     RAINF     BDIRT     GRAC4
      real*8, DIMENSION(N_PFT), parameter :: wdens_g_cm3 =
     &     (/ 0.6d0,undef,0.6d0,0.6d0,0.54d0,0.5d0,0.7d0,undef /)
      real*8, DIMENSION(N_PFT), parameter :: a0h = !Later move to ent_pfts_ENT.f
     &     (/ 0.d0,undef,0.d0,0.d0,1.3d0,0.d0,0.d0,undef /)
      real*8, DIMENSION(N_PFT), parameter :: acr =
     &     (/ 0.500d0,undef,0.500d0,0.3868d0,0.2773d0
     &        ,0.2570d0,0.1407d0,0.500d0 /)
      real*8, DIMENSION(N_PFT), parameter :: bcr =
     &     (/ 1.d0,undef,1.d0,1.d0,1.d0,1.d0,1.d0,undef /)
      real*8, DIMENSION(N_PFT), parameter :: bR =
     &     (/ 1.d0,40.d0,1.d0,0.75d0,1.0d0,2.0d0,1.5d0,10.0d0 /)
      integer, DIMENSION(N_PFT), parameter  :: form =
     &     (/ SHRUB,HERB,SHRUB,TREE,TREE,TREE,TREE,HERB /)
      real*8, DIMENSION(N_PFT), parameter  :: DBHBAmax_cm =
     &     (/ 10.d0,undef,10.d0,150.d0,150.d0,150.d0,150.d0,undef /)
      logical, DIMENSION(N_PFT), parameter  :: crop =
     &     (/ .false.,.false.,.false.,.false.
     &       ,.false.,.false.,.false.,.true. /)

      !***************************************************
      !*         END - GISS VEGETATION TYPES             *
      !***************************************************
      end module ent_pfts
