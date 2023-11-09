      module ent_const
!@sum  CONSTANT definitions for physical constants
!@sum  These are constants that would be common to the GCM/EWB, so should
!@sum  have the GCM/EWB constants substituted in for coupled runs.
!@auth N.Kiang
!@ver  1.0

      !* COUPLED RUNS - Replace with values from GCM constants*!
      !use Name_of_GCM_constants_module  

      use TimeConstants_mod, only: sday=>SECONDS_PER_DAY
      implicit none
      save

#define OFFLINE 1
      !************************************************************************
      !* INTEGRATION - TIME STEPS
      integer,parameter :: T_SUB=12  !Number of sub-time steps in a year
      !************************************************************************
      !* NUMERICAL CONSTANTS
      real*8,parameter :: PI = 3.1415926535897932d0 !@param pi    pi
      real*8,parameter :: zero = 0.d0

      !************************************************************************
      !* PHYSICAL CONSTANTS - may set these equal to constants of GCM.
      real*8,parameter :: stbo =5.67051d-8 !Stefan-Boltzman (W/m2/K4)
      real*8,parameter :: lhe0 = 2.5008d6 !Latent heat evap at 0 oC (2.5008d6 J/kg)
      real*8,parameter :: lhew = 2.260d6 !Latent heat evap at 100 oC (J/kg)
      real*8,parameter :: shv100 = 2080.d0 !Specific heat of water vapor @ 100 oC (J kg-1 K-1)
      real*8,parameter :: shw25 = 4181.3d0 !Specific heat of water liquid @ 25 oC (J kg-1 K-1) (Wiki)
      real*8,parameter :: shw20 = 4185.d0 !Specific heat of water liquid @ 20 oC (J kg-1 K-1) (GISS GCM)
      real*8,parameter :: rhow = 1d3  !Density of pure water (1000 kg/m^3)
      real*8,parameter :: tfrz = 273.15d0 !freezing pt of H2O at 1 atm (Kelvin)
      real*8,parameter :: KELVIN = tfrz
      real*8,parameter :: gasc = 8.314510d0 !gas constant (8.314510 J/mol K)
      real*8,parameter :: Avogadro=6.023d23 !Avogadro's constant (atmos/mole)
      real*8,parameter :: cp=1012. !Heat capacity of dry air (J kg-1 K-1)
!@param shw heat capacity of water (at 20 C) (4185 J/kg C)
      real*8,parameter :: shw  = 4185.

      !************************************************************************
      !* CONSTANTS FROM GISS CONST.f
      !@param mwat molecular weight of water vapour
      real*8,parameter :: mwat = 18.015d0
      !@param mair molecular weight of dry air (28.9655 g/mol)
      real*8,parameter :: mair = 28.9655d0
      !@param mrat  mass ratio of air to water vapour (0.62197)
      real*8,parameter :: mrat = mwat/mair    ! = 0.62197....
      !@param rvap  gas constant for water vapour (461.5 J/K kg)
      !**** defined as R/M_W = 1000* 8.314510 J/mol K /18.015 g/mol
      real*8,parameter :: rvap = 1d3 * gasc / mwat ! = 461.5...
      !@param tf freezing point of water at 1 atm (273.15 K)
      real*8,parameter :: tf = tfrz

      !************************************************************************
       !* ASTRONOMICAL CONSTANTS

      !************************************************************************
      !* RUN CONTROL
      real*8,parameter :: undef=-1.d30 ! Missing value
      real*8,parameter :: teeny=1.d-30 ! Small positive value to avoid 0/0
      real*8, parameter :: EPS = 1.d-8  !Small error
      real*8, parameter :: EPS2 = 1.d-12 !Smaller error



      
      !************************************************************************
      !* Ent CONSTANTS
      !******************
      !* PATCH DYNAMICS *
      !******************
      real, parameter ::  F_AREA=.01 !* min area of patch as fraction of total
      real, parameter ::  BTOL=.00001!* min cohort biomass for termination (kgC/m2)
      real, parameter ::  NTOL=.001  !* min plant density for whatever
      integer,parameter :: PATCH_DYNAMICS = 0 ! 0-No, 1=Yes

      !************************************************************************
       !* ASTRONOMICAL CONSTANTS
      ! real*8,parameter :: sday = 86400.d0! sec per day (s)
      real*8,parameter :: SECPY = 31536000.d0  ! sec per year (s)

      !************************************************************************
      !********************
      !* SOIL / HYDROLOGY *
      !********************
!      integer :: N_DEPTH        !Number of soil layers.  SET IN ENT_INIT
      integer, parameter :: N_DEPTH = 6 !Number of soil layers. SET AS PAR BUT LATER CAN BE VARIABLE.
      real*8, parameter :: SOILDEPTH_m(N_DEPTH) = !Bottom depths of soil layers (m)
     &     (/ 0.1,0.27,0.57,1.08,1.97,3.5 /) !GISS GCM

      !**********************
      !* RADIATIVE TRANSFER *
      !**********************
      integer,parameter :: N_BANDS = 6 !Number of spectral bands (GISS 6)
                                !Expect to adjust to hyperspectral

      !***********************
      !* ECOLOGICAL DYNAMICS *
      !***********************
      integer,parameter :: N_DIST_TYPES = 2 !Number of disturbance types
      real*8,parameter :: LOW_PAR_LIMIT = 2.5d0 !umol m-2 s-1.  Nobel 1999, lower light limit for green plants is 0.7 W m-2 ~ 3 umol m-2 s-1.

!#define PFT_MODEL_ENT
#ifdef PFT_MODEL_ENT
      !************************************************
      !*  ENT PLANT FUNCTIONAL TYPE CONSTANTS         *
      !************************************************
      integer,parameter :: N_PFT = 16
      integer,parameter :: N_OTHER = 0
      integer,parameter :: COVER_SAND = 17
      integer,parameter :: COVER_DIRT = 18
#else
      !************************************************
      !*  GISS VEGETATION CONSTANTS                   *
      !************************************************
      integer,parameter :: N_PFT = 8 !8
      integer,parameter :: N_OTHER = 2 !2 algae, c4 grass
      integer,parameter :: COVER_SAND = 1
      integer,parameter :: COVER_DIRT = 10
#endif
      !************************************************
      !* GENERIC VEGETATION CONSTANTS                 *
      !************************************************
      !* photosynthetic pathway type pst in pftype
      integer,parameter :: C3 = 1
      integer,parameter :: C4 = 2

      !* leaftype in pftype *!
      integer,parameter :: BROADLEAF = 1
      integer,parameter :: NEEDLELEAF = 2
      integer,parameter :: MONOCOT = 3 !NOT CROPS

      !* phenotype in pftype *!
      integer,parameter :: EVERGREEN = 1
      integer,parameter :: COLDDECID = 2
      integer,parameter :: DROUGHTDECID = 3
      integer,parameter :: COLDDROUGHTDECID = 4
      integer,parameter :: ANNUAL = 5 

      !* growth form *!
      integer, parameter :: GRASS = 1
      integer, parameter :: HERB = 2
      integer, parameter :: SHRUB = 3
      integer, parameter :: TREE = 4
      integer, parameter :: BARE = 5


      !************************************************
      !*  COVER SUMMARY CONSTANTS                     *
      !************************************************

      integer,parameter :: N_SOILCOV = 2 !2-light sand, dark dirt (GISS) 
      integer,parameter :: N_COVERTYPES = N_PFT + N_SOILCOV + N_OTHER


      !************************************************************************
      !* COHORT Biomass pools
      integer,parameter :: N_BPOOLS = 7 !foliage,sapwood,hardwood,labile, fine root,coarse root
      integer,parameter :: FOL = 1   !FOLIAGE Array indices for array(N_BPOOLS)
      integer,parameter :: SW = 2    !SAPWOOD
      integer,parameter :: HW = 3    !HARDWOOD
      integer,parameter :: LABILE = 4 !LABILE
      integer,parameter :: FR = 5    !FINE ROOT
      integer,parameter :: CR = 6 !COARSE ROOT
      integer,parameter :: RP = 7 !REPRODUCTION
      !************************************************************************
      !************************************************************************
      !* CASA SOIL CONSTANTS
      !* THESE SHOULD GO IN THE CASA MODULE BUT NEED TO BE USED IN ENT_TYPES.F
      !* Array sizes
      integer,parameter :: N_SOIL_TEXTURES = 5  !in GHY.f, textures are sand,loam,clay,peat(+bedrock) -PK 7/13/06
      integer,parameter :: NLIVE = 3
      integer,parameter :: NDEAD = 9
      integer,parameter :: NPOOLS = NLIVE + NDEAD
      !layers for soil bgc -- need to make this into variable if
      !affected arrays later allocated dynamically -PK 7/23/07   
#ifdef NCASA2
      integer,parameter :: N_CASA_LAYERS = 2  !Option for 2 soil biogeochemistry layers, 0-30 cm, 30-100 cm.
#else 
      integer,parameter :: N_CASA_LAYERS = 1   !Default option for 1 soil biogeochemistry layer, 0-30 cm.
#endif
      !* Total pool array indices
      integer,parameter :: PTRACE = 2 !Trace elements in Tpools, C and N
      integer,parameter :: Carbon = 1
      integer,parameter :: Nitrogen = 2
!      integer,parameter :: ptrace = 2  !num. nutrient pools used in CASA resp. routine -PK
      integer,parameter :: NRESP_PATHS = 14  !num. pathways between pools in CASA soil respiration.
      real*8,parameter :: Q10 = 2.d0        !Q10 used in belowground calculations --> value from lit -PK 5/25/06

      !* Live pool array indices
      integer,parameter :: LEAF = 1  !Array index
      integer,parameter :: FROOT = 2  !Array index
      integer,parameter :: WOOD = 3  !Array index

      !* Dead pool array indices
      integer,parameter :: SURFMET = 4 !sfc metabolic
      integer,parameter :: SURFSTR = 5 !sfc structural
      integer,parameter :: SOILMET = 6 !soil metabolic
      integer,parameter :: SOILSTR = 7 !soil structural
      integer,parameter :: CWD = 8     !coarse woody debris
      integer,parameter :: SURFMIC = 9 !sfc microbial
      integer,parameter :: SOILMIC = 10 !soil microbial
      integer,parameter :: SLOW = 11 !slowly decomposing soil o.m. pool (up to a decade)
      integer,parameter :: PASSIVE = 12 !very slowly decomposing soil o.m. pool (decades-centuries)

      !************************************************************************
      real*8 :: CNratio(NPOOLS)
            data CNratio/
     1            30.0,       ! C:N ratio of leaf pool
     2           130.0,       ! C:N ratio of wood pool
     3            55.0,       ! C:N ratio of froot pool
     4            30.0,       ! C:N ratio of surfmet pool
     5            50.0,       ! C:N ratio of surfstr pool
     6            25.0,       ! C:N ratio of soilmet pool
     7            50.0,       ! C:N ratio of soilstr pool
     8           135.0,       ! C:N ratio of cwd pool
     9            12.5,       ! C:N ratio of surfmic pool
     a            12.5,       ! C:N ratio of soilmic pool
     b            12.5,       ! C:N ratio of slow pool
     c             8.5/       ! C:N ratio of passive pool

      real*8,dimension(N_PFT,NPOOLS) :: annK !CASA turnover times (sec-1)
      !real*8,dimension(N_PFT,NPOOLS) :: kdt !CASA turnover times for wood & dead (yr-1)
      real*8,dimension(N_PFT) :: solubfract !Soluble ("metabolic") fraction of litter 
      real*8,dimension(N_PFT) :: structurallignin !fraction of structural C from lignin -PK 7/5/06 
      real*8,dimension(N_PFT) :: lignineffect !effect of lignin on decomp -PK 7/5/06
      real*8,parameter :: woodligninfract = 0.40 !amt lignin in wood C -PK 7/5/06
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            character*13, parameter :: Ent_cpool_title(N_BPOOLS) =
     &     (/
     &     'FOLIAGE      '
     &     ,'SAPWOOD      '
     &     ,'HARDWOOD     '
     &     ,'LABILE       '
     &     ,'FINE ROOT    '
     &     ,'COARSE ROOT  '
     &     ,'REPRODUCTION '
     &     /)

      !************************************************************************
      end module ent_const
