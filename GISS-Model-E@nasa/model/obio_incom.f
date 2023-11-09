#include "rundeck_opts.h"

      MODULE obio_incom

! parameters and arrays neccessary for obio_init and obio_bioinit

      USE obio_dim
      USE obio_com, only: sday
      USE ocalbedo_mod, only: nlt

      implicit none

      real :: rmumax(nchl)          !max phyto growth rate at 20oC, d/
      real :: rik(3,nchl)           !light saturation parameter umol quanta/m2/s
      real :: obio_wsd(nchl)        !phyto sinking rate m/d
      real :: obio_wsh(nchl)        !phyto sinking rate m/h
      real :: obio_wss(nchl)        !phyto sinking rate m/s

      real :: rkn(nchl)   !half-saturation constants for nitrogen (uM)
      real :: rks(nchl)   !half-saturation constants for silica (uM)
      real :: rkf(nchl)   !half-saturation constant for iron (nM)

      real rad, pi2

      real Pdeep(ntyp)            !deep BC
      real detdeep(ndet)          !detrital deep BC
      real cardeep(ncar)          !carbon deep BC
 
      real :: cchl(3)             !C:chl ratio g:g for 3 light states

      real :: cnratio             !C:N ratio
      real :: csratio             !C:Si ratio
      real :: cfratio             !C:Fe ratio
      real :: ro2c_DET            !O2:C ratio from detritus
      real :: ro2c_NH4            !O2:C ratio from ammonium production
      real :: ro2c_NO3            !O2:C rati from NO3 productio    
      
      real :: mgchltouMC

      real :: wsdeth(ndet)        !sinking rate of detritus (m/d): changed to m/s    !July 2016
      real :: remin(ndet)         !detrital remineralization rate /d: changed to m/s !July 2016
#ifdef TRACERS_degC
      real :: tdegC               ! transfer rate constant from degradable to nondegradable carbon !Jan 2020
#endif 
      
      real :: Fescavrate(2)       !scavenging rate for dissolved iron
      real:: HvO2                !@PL dirac delta function for O2 consumption processes
      real:: NCrrat              !@PL ratio of anaerobioc/aerobic remin
C if CARBON == 1
      real, parameter :: excp=0.05              !excretion of DOC by phyto growth
      real, parameter :: resp=0.05              !respiration of DIC by phyto growth
!     real, parameter :: excz=0.05/24.0         !excretion of DOC by zoopl/hr
      real, parameter :: excz=0.05/sday         !excretion of DOC by zoopl/s        !July 2016
!     real, parameter :: resz=0.05/24.0         !respiration of DIC by zoopl/hr
      real, parameter :: resz=0.05/sday         !respiration of DIC by zoopl/s      !July 2016
      real, parameter :: phygross=1.0-(excp+resp) !factor to derive gross PP!

!change: March 15, 2010
!     real, parameter :: rlamdoc=0.017/24.0    !N-dependent DOC remin/hr
!     real, parameter :: rlamdoc=0.005/24.0    !N-dependent DOC remin/hr
      real, parameter :: rlamdoc=0.005/sday    !N-dependent DOC remin/s    !July 2016

      real, parameter :: rkdoc1=0.3*10.0       !N-dep half sat uM(PO4) modified
                                               !for nitrate by mult*10.0,
                                               !based on Conkright et al. 1994
      real, parameter :: rkdoc2=15.0           !DOC-dep half-sat uM(C)
!     real, parameter :: rlampoc=0.05/24.0     !detrital breakdown/hr
      real, parameter :: rlampoc=0.05/sday     !detrital breakdown/s      !July 2016
      real, parameter :: uMtomgm3=12.0         !conversion uM to mg/m3 C
      real, parameter :: Pzo=1.0*uMtomgm3/50.0 !zoopl half-sat for
                                               !DOC excretion mg/m3(chl,assuming
                                               !C:chl ratio of 50))
      real, parameter :: stdslp=1013.25        !standard sea level pressure in mb
      real, parameter :: ko2=20e-3             !@PL half-saturation constant for transfer from aerobic remin to denit
      real, parameter :: O2thr=2e-3            !@PL O2 concentration threshold for cessation of aerobic remin/respiration

!     real, parameter :: Rm=1.20/24.0          !max zoopl. growth rate/hr
      real, parameter :: Rm=1.20/sday          !max zoopl. growth rate/s    !July 1016
                                               !increase to account for excretion
                                               !and respiration
C if CARBON /=1    parameter(Rm=1.0/24.0)      !max zoopl. growth rate/hr


      !array of factors to compute mean irradiance w/in water column
      real facirr(nh,nch,5,ncd)  

      !absorption and scattering coefficients of chlorophyll
      real ac(nchl,nlt),bc(nchl,nlt)

!     real, parameter :: solFe=0.02    !solubility of iron: this is the default
!     real, parameter :: solFe=0.05    !solubility of iron
      real solFe                       !now defined in rundeck

c     parameter(bn=0.5,bs=0.5)        !N/chl and Si/chl ratios
      
      real bn,bs,bf,cchlratio


      integer nl450

      real excdom(nlt),bbw,Dmax,rd,ru,rmus,rmuu

#ifdef exp_wsdiat
      real :: adiat_exp, bdiat_exp
#endif
#ifdef exp_wsdet
      real :: adet_exp(3), bdet_exp(3)
#endif

!define compensation depth
      real, parameter ::  zc = 75. ! in meters (from OCMIP)

#ifdef TRACERS_Alkalinity
      real, parameter ::  rain_ratio=0.07
      real, parameter ::  npratio   =16     ! N:P Redfield ratio
      real, parameter ::  cpratio   =117    ! C:P Redfield ratio,
                                            ! these values are from OCMIP 
                                            ! protocol (Najjar and Orr, 1998)
                   ! Yamanaka and Tajika(1996) suggest R=0.08, rCP=106

! sigma  fraction of net downward flux of organic matter across the
!        compensation depth that is in dissolved form, ie form of
!        export production that is in dissolved form
      real, parameter ::  sigma_Ca  = 0.67 ! (Yamanaka and Tajika, 1997)

      real, parameter ::  d_Ca = 3500.    ! in meters (Yamanaka and Tajika, 1996)
      real, parameter ::  kappa_Ca = 2.    ! (1/0.5years)^-1   !OCMIP
#endif

#ifdef OBIO_RUNOFF
      real, parameter ::  estFe = 0.01   ! estuarine retention rate, can vary between 0.2 and 0.01 (daCunha 2007)
#endif


      END MODULE obio_incom
