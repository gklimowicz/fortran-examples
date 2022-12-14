#include "rundeck_opts.h"

      module mo_bulk2m_driver_gcm
!@sum 2-moment bulk cloud microphysics scheme                                
!@contains routines to calculate cloud and ice crystal number, autoconversion, phase-interactions 
!@auth Surabi Menon and Igor Sednev ; based on Morrison 2008 scheme in CCSM
!@+    as in Morrison and Gettelman, 2008, J Clim 21, 15, 3642-3659 and Gettelman et al. 2008, J Clim,21,3660-3679
      USE resolution,ONLY: im,jm,lm
      USE CONSTANT, only : tf,rgas,rvap,pi,grav
      IMPLICIT NONE
      PRIVATE
      integer,parameter :: mx=1 !,mx=lm
        INTEGER      :: il0,jl0,kl0,nmodes,i,j,k,ktau
        INTEGER      :: ier,ios,ou=6
        REAL*8       :: dt,rhw,rhi
        REAL*8       :: tp0,tk0,qk0,pk0
        LOGICAL      :: lcheck,ldummy,ltermo
        CHARACTER*80 :: msg
C==============================================================================
C Include file: blk_add_declare.f
C========================================================================
        REAL*8,ALLOCATABLE,DIMENSION(:,:) :: na4d    !#/m3 
        REAL*8,ALLOCATABLE,DIMENSION(:) :: 
     1       tk3d,           ! temperature (K)
     1       qv3d,           ! water vapor mixing ratio (kg/kg)
     1       pp3d,           ! atmospheric pressure (pa)
     1       pdel,           ! difference in pressure across vertical level 
     1       sw3d,           ! supersaturation with respect to water
     1       si3d,           ! supersaturation with respect to ice
c    1       rttd,           ! radiative heating rate (K/s)
     1       ww3d,           ! grid-scale vertical velocity (m/s)
     1       wvar,           ! sub-grid vertical velocity (m/s)
     1       tk3dten,        ! temperature tendency (K/s)
     1       qv3dten         ! water vapor mixing ratio tendency (kg/kg/s)
C Allocatable 1d arrays: 
C input/output parameters           ! description (units)
	REAL*8,ALLOCATABLE,DIMENSION(:)  :: ! (mkx)
     1       qc3dten,        ! cloud water mixing ratio tendency (kg/kg/s)
     1       qi3dten,        ! cloud ice mixing ratio tendency (kg/kg/s)
     1       qs3dten,        ! snow mixing ratio tendency (kg/kg/s)
     1       qr3dten,        ! rain mixing ratio tendency (kg/kg/s)
     1       nc3dten,        ! cloud droplet number concentration (1/kg/s)
     1       ni3dten,        ! cloud ice number concentration (1/kg/s)
     1       ns3dten,        ! snow number concentration (1/kg/s)
     1       nr3dten         ! rain number concentration (1/kg/s)
	REAL*8,ALLOCATABLE,DIMENSION(:)  :: ! (mkx)
     1       qc3d,           ! cloud water mixing ratio (kg/kg)
     1       qi3d,           ! cloud ice mixing ratio (kg/kg)
     1       qs3d,           ! snow mixing ratio (kg/kg)
     1       qr3d,           ! rain mixing ratio (kg/kg)
     1       nc3d,           ! cloud droplet number concentration (1/kg)
     1       ni3d,           ! cloud ice number concentration (1/kg)
     1       ns3d,           ! snow number concentration (1/kg)
     1       nr3d            ! rain number concentration (1/kg)
C output variables
	REAL*8,ALLOCATABLE,DIMENSION(:)  :: ! (mkx)
     1       effc,           ! droplet effective radius (micron)
     1       effi,           ! cloud ice effective radius (micron)
     1       effs,           ! snow effective radius (micron)
     1       effr            ! rain effective radius (micron)
	REAL*8
     1       precrt,         ! precip rate (mm/hr) (rain plus snow)
     1       snowrt          ! snowfall rate (mm/hr)
C size parameter variables
	REAL*8,ALLOCATABLE,DIMENSION(:)  :: ! (mkx)
     1       lamc,           ! slope parameter for droplets (m^-1)
     1       lami,           ! slope parameter for cloud ice (m^-1)
     1       lams,           ! slope parameter for snow (m^-1)
     1       lamr,           ! slope parameter for rain (m^-1)
     1       cdist1,         ! PSD parameter for droplets
     1       n0i,            ! intercept parameter for cloud ice (kg-1 m-1)
     1       n0s,            ! intercept parameter for snow (kg-1 m-1)
     1       n0rr,           ! intercept parameter for rain (kg-1 m-1)
     1       pgam            ! spectral shape parameter for droplets
c microphysical processes
	REAL*8,ALLOCATABLE,DIMENSION(:)  ::     ! (mkx)
     1       nsubc,              ! loss of nc during evap
     1       nsubi,              ! loss of ni during sub.
     1       nsubs,              ! loss of ns during sub.
     1       nsubr,              ! loss of nr during evap
c
c     1       prd,                ! dep/sub cloud ice
c     1       pre,                ! cond/evap of rain
c     1       prds,               ! dep/sub snow
c
     1       ncondc,             ! change n cond/evap droplets
     1       ncondi,             ! change n dep/sub cloud ice
     1       ncondr,             ! change n cond/evap of rain
     1       nconds,             ! change n dep/sub snow
     1       mcondc,             ! change q cond/evap droplets
     1       mcondi,             ! change q dep/sub cloud ice
     1       mcondr,             ! change q cond/evap of rain
     1       mconds,             ! change q dep/sub snow
c     1       npcc,               ! cond/evap droplets number
c     1       mpcc,               ! cond/evap droplets mass
     1       nnuccc,             ! change n due to con. droplets freez 
     1       mnuccc,             ! change q due to con. droplets freez 
     1       nnucci,             ! change n due to imm. droplets freez 
     1       mnucci,             ! change q due to imm. droplets freez
     1       npccn,              ! change n droplets activation
     1       mpccn,              ! change q droplet activation
     1       nnuccd,             ! change n freezing aerosol (prim ice nuc)
     1       mnuccd,             ! change q freezing aerosol (prim ice nuc) 
     1       nnucmd,             ! change n cond freezing Meyers (prim ice nuc)
     1       mnucmd,             ! change q cond freezing Meyers (prim ice nuc) 
     1       nnucmt,             ! change n cont freezing Meyers (prim ice nuc)
     1       mnucmt,             ! change q cont freezing Meyers (prim ice nuc) 
     1       nnuccr,             ! change n due to con rain freez 
     1       mnuccr,             ! change q due to con rain freez 
     1       nnucir,             ! change n due to imm rain freez 
     1       mnucir,             ! change q due to imm rain freez  
     1       ncactv,             ! activated fraction of aerosols for all modes 
c............................................................................
c water-water interactions
c............................................................................
     1       nprc,               ! change n autoconversion droplets
     1       mprc,               ! change q autoconversion droplets
     1       ncagg,              ! change n self-collection droplets
     1       mcagg,              ! change q self-collection droplets is ZERO
     1       npra,               ! change n droplets accretion by rain
     1       mpra,               ! change q droplets accretion by rain
     1       nragg,              ! change n self-collection rain
     1       mragg,              ! change q self-collection rain is ZERO
c
c ice-ice interactions
c crystal-crystal
     1       nprci,              ! change n autoconversion cloud ice
     1       mprci,              ! change q autoconversion cloud ice
     1       niagg,              ! change n self-collection cloud ice
     1       miagg,              ! change q self-collection cloud ice is ZERO
c crystal-snow
     1       nprai,              ! change n accretion cloud ice by snow
     1       mprai,              ! change q accretion cloud ice by snow
c snow-snow
     1       nsagg,              ! change n self-collection of snow
     1       msagg,              ! change q self-collection of snow
c
c water-ice interactions
c drop-crystal
c drop-snow
     1       npsacws,            ! change n droplet acc droplets by snow
     1       mpsacws,            ! change q droplet acc droplets by snow
     1       nmults,             ! change n ice mult acc droplets by snow
     1       mmults,             ! change q ice mult acc droplets by snow
c rain-crystal
c rain-snow
     1       npracs,             ! change n collection rain by snow
     1       mpracs,             ! change q collection rain by snow
     1       nmultr,             ! change n due to ice mult acc rain by snow
     1       mmultr,             ! change q due to ice mult acc rain by snow
     1       nsmlts,             ! change n melting snow evaporating
     1       msmlts,             ! chnage q melting snow evaporating
     1       nsmltr,             ! change n melting snow to rain
     1       msmltr,             ! change q melting snow to rain
     1       nmltii,             ! change n melting cloud ice evaporating
     1       mmltii,             ! change q melting cloud ice evaporating
     1       nmltic,             ! change n melting cloud ice to droplets
     1       mmltic,             ! change q melting cloud ice to droplets
     1       nhfrc,              ! change n homog freezing droplets to ice
     1       mhfrc,              ! change q homog freezing droplets to ice
     1       nhfrr,              ! change n homog freezing rain to graupel
     1       mhfrr,              ! change q homog freezing rain to graupel
     1       nsedimc,            ! change n sedimentation drop
     1       msedimc,            ! change q sedimentation drop
     1       nsedimi,            ! change n sedimentation ice
     1       msedimi,            ! change q sedimentation ice
     1       nsedimr,            ! change n sedimentation rain
     1       msedimr,            ! change q sedimentation rain
     1       nsedims,            ! change n sedimentation snow
     1       msedims,            ! change q sedimentation snow
     1       nvtermc,            ! n-weighted sedim vel drop
     1       mvtermc,            ! q-weighted sedim vel drop
     1       nvtermi,            ! n-weighted sedim vel ice
     1       mvtermi,            ! q-weighted sedim vel ice
     1       nvtermr,            ! n-weighted sedim vel rain
     1       mvtermr,            ! q-weighted sedim vel rain
     1       nvterms,            ! n-weighted sedim vel snow
     1       mvterms             ! q-weighted sedim vel snow
c
	REAL*8,ALLOCATABLE,DIMENSION(:)  ::     ! (mkx)
     1       wnuc                ! vertical velocity for primary nucleation
     1      ,anuc                ! "a" coefficient for primary crys nucleation
     1      ,bnuc                ! "b" coefficient for primary crys nucleation
c time-varying atmospheric parameters
	REAL*8,ALLOCATABLE,DIMENSION(:)  ::     ! (mkx)
     1       kap,                ! thermal conductivity of air
     1       evs,                ! saturation vapor pressure
     1       eis,                ! ice saturation vapor pressure
     1       qvs,                ! saturation mixing ratio
     1       qvi,                ! ice saturation mixing ratio
     1       qvqvs,              ! watre satration ratio
     1       qvqvsi,             ! ice   saturaion ratio
     1       Dv,                 ! diffusivity of water vapor in air
     1       Dap,                ! diffusivity of aerosol
     1       xxls,               ! latent heat of sublimation
     1       xxlv,               ! latent heat of vaporization
     1       mu,                 ! viscocity of air
     1       sc,                 ! Schmidt number
     1       xlf,                ! latent heat of freezing
     1       rho,                ! air density
     1       ab,                 ! correction to cond rate due to latent heating
     1       abi,                ! correction to depo rate due to latent heating
     1       fdum,               ! mean free path 
     1       eii                 ! collection efficiency 
	REAL*8,ALLOCATABLE,DIMENSION(:)  ::     ! (mkx)
     1       epsc,               ! 1/phase rel. time (see M2005), drop
     1       epsi,               ! 1/phase rel. time (see M2005), crys
     1       epss,               ! 1/phase rel. time (see M2005), snow
     1       epsr                ! 1/phase rel. time (see M2005), rain
C
C fall speed working variables (defined in code)
c	REAL*8,ALLOCATABLE,DIMENSION(:)  ::     ! (mkx)
c     1       dumfmi,dumfmr,dumfni,
c     1       fmr,
c     1       fmi,fni,faloutmr,faloutmi,
c     1       faloutni,
c     1       dumqs,dumfns,fms,fns,
c     1       faloutms,faloutns,
c     1       dumfmc,dumfnc,fmc,faloutmc,
c     1       faloutnc,fnc,
c     1       dumfnr,faloutnr,fnr
C
C fall speed working variables (defined in code)
	REAL*8,ALLOCATABLE,DIMENSION(:)  ::     ! (mkx)
     1       fnc,dumfnc,faloutnc,
     1       fmc,dumfmc,faloutmc,
     1       fni,dumfni,faloutni,
     1       fmi,dumfmi,faloutmi,
     1       fnr,dumfnr,faloutnr,
     1       fmr,dumfmr,faloutmr,
     1       fns,dumfns,faloutns,
     1       fms,dumfms,faloutms
C
C fall-speed parameter 'a' with air density correction
	REAL*8,ALLOCATABLE,DIMENSION(:)  ::     ! (mkx)
     1       ain,arn,asn,acn
C
C microphysics parameters and physical constants, defined during first time step
	REAL*8 
     1     ai,ac,as,ar  ! 'a' parameter in fallspeed-diam relationship
	REAL*8 
     1     bi,bc,bs,br  ! 'b' parameter in fallspeed-diam relationship
	REAL*8 
     1     r,           ! gas constant for dry air
     1     rv,          ! gas constant for water vapor
     1     cp,          ! specific heat at constant pressure for dry air
     1     rhosu,       ! standard air density at 850 mb
     1     rhow,        ! density of liquid water
     1     rhoi,        ! bulk density of cloud ice
     1     rhosn,       ! bulk density of snow
     1     aimm,        ! parameter in Bigg immersion freezing
     1     bimm,        ! parameter in Bigg immersion freezing
     1     ecr,         ! collection efficiency droplets/rain and snow/rain  
     1     dcs,         ! threshold size for cloud ice autoconversion  
     1     rw0,         ! initial size of nucleated droplet, [m]   
     1     mw0,         ! initial mass of nucleated droplet, [kg]   
     1     ri0,         ! initial size of nucleated crystal, [m]   
     1     mi0,         ! initial mass of nucleated crystal, [kg]   
     1     f1s,         ! ventilation parameter for snow
     1     f2s,         ! ventilation parameter for snow
     1     f1r,         ! ventilation parameter for rain
     1     f2r,         ! ventilation parameter for rain
     1     g,           ! gravitational acceleration
     1     qsmall,      ! smallest allowed hydrometeor mixing ratio
     1     nsmall       ! smallest allowed hydrometeor mixing concentration

	REAL*8 
     1     ci,di,cs,ds, ! size distribution parameters for cloud ice and snow
     1     dr,          ! size distribution parameters for rain
     1     rin          ! radius of contact nuclei (m)

c aerosol parameters
	REAL*8 
     1     mw,          ! molecular weight water (kg/mol)
     1     osm,         ! osmotic coefficient
     1     vi,          ! number of ion dissociated in solution
     1     epsm,        ! aerosol soluble fraction
     1     rhoa,        ! aerosol bulk density (kg/m3)
     1     mapp,        ! molecular weight aerosol (kg/mol)
     1     ma,          ! molecular weight of 'air' (kg/mol)
     1     rr,          ! universal gas constant
     1     bact,        ! activation parameter
     1     rm1,         ! geometric mean radius, mode 1 (m)
     1     rm2,         ! geometric mean radius, mode 2 (m)
     1     nanew1,      ! total aerosol concentration, mode 1 (m^-3)
     1     nanew2,      ! total aerosol concentration, mode 2 (m^-3)
     1     sig1,        ! standard deviation of aerosol s.d., mode 1
     1     sig2,        ! standard deviation of aerosol s.d., mode 2
     1     f11,         ! correction factor for activation, mode 1
     1     f12,         ! correction factor for activation, mode 1
     1     f21,         ! correction factor for activation, mode 2
     1     f22          ! correction factor for activation, mode 2
c
	LOGICAL,ALLOCATABLE,DIMENSION(:) ::     ! (mkx)
     1     lcold,lwarm,
     1     lqc3d,lqi3d,lqs3d,lqr3d,
     1     lnc3d,lni3d,lns3d,lnr3d

C==============================================================================
        INTERFACE init_bulk2m_driver
!         module procedure init_bulk2m_driver_gcm
          module procedure init_bulk2m_driver_gcm_mat
        END INTERFACE init_bulk2m_driver
C
        INTERFACE make_distribution
          module procedure make_drop_dist
          module procedure make_rain_dist
          module procedure make_crys_dist
          module procedure make_snow_dist
          module procedure make_xxxx_dist
          module procedure make_distribution0
        END INTERFACE make_distribution
C
        INTERFACE execute_bulk2m_driver
          module procedure get_rate_bulk2m
          module procedure get_value_bulk2m
          module procedure set_termo_bulk2m
          module procedure set_micro_bulk2m
          module procedure set_micro_bulk2m_mat
          module procedure set_micro_bulk2m_tomas
          module procedure update_micro_tendencies_bulk2m
          module procedure execute_bulk2m_sedimentation_gcm
c          module procedure execute_bulk2m_driver0
          module procedure execute_bulk2m_driver_gcm
c          module procedure execute_bulk2m_driver_gcm0
c          module procedure execute_bulk2m_driver_gcm1
c          module procedure hugh_ice_nucleation
c          module procedure hugh_ice_nucleation0
c          module procedure hugh_drop_nucleation
c          module procedure hugh_autoconversion
        END INTERFACE execute_bulk2m_driver
C
      PUBLIC    :: init_bulk2m_driver
      PUBLIC    :: execute_bulk2m_driver
      PUBLIC    :: cleanup_bulk2m_driver

C
C==============================================================================
C
        CONTAINS
!
C==============================================================================
!     LOGICAL FUNCTION init_bulk2m_driver_gcm(dt0,il0,jl0,kl0,name) 
!    *  RESULT (la)
!     IMPLICIT NONE
!       REAL*8,             INTENT (in)       :: dt0
!       INTEGER,OPTIONAL,   INTENT (in)       :: il0,jl0,kl0
!       CHARACTER (len = *),INTENT (in),optional :: name
!       CHARACTER (len = 120)                 :: fname
!       INTEGER                               :: ier,i,ios,ou,k
!       integer :: il,jl,kl
!       CHARACTER*24 :: sname='init_bulk2m_driver_gcm: '
!       dt = dt0
!       if(PRESENT(il0)) then 
!         il = il0
!       else 
!         il=im
!       endif
!       if(PRESENT(jl0)) then 
!         jl = jl0
!       else 
!         jl=jm
!       endif
!       if(PRESENT(kl0)) then 
!         kl = kl0
!       else 
!         kl=lm
!       endif
!       if(PRESENT(name)) then 
!         fname = trim(name)
!       else 
!         fname='blk_input.dat'
!       endif
!       if(il .eq. 1) ltermo=.true.
!       lcheck=.false.
!       ou = 6; ktau = 0; i=1; j = 1; rhw=0.95; rhi=0.95
! Allocate TERMO arrays
!        ldummy = allocate_arrays_bulk2m_driver(mx)
!       ldummy = allocate_arrays_bulk2m_driver(kl)
! Hugh 
! constants and parameters
!       ldummy=hugh_constants_init()
!       ldummy=hugh_drop_nucleation_init()
!       if(ldummy) then
!         print *,sname,'Initialization has been completed...'
!       else
!         lcheck = .true.
!         print *,sname,'BULK2M initialization has not been done ...'
!         ldummy=cleanup_bulk2m_driver()
!       endif
!       la = .NOT.lcheck
!     END FUNCTION init_bulk2m_driver_gcm
C==============================================================================
      LOGICAL FUNCTION init_bulk2m_driver_gcm_mat(dt0,il0,jl0,kl0,nm0,
     *name) 
     *  RESULT (la)
      IMPLICIT NONE
        REAL*8,             INTENT (in)       :: dt0
        INTEGER,OPTIONAL,   INTENT (in)       :: il0,jl0,kl0,nm0
        CHARACTER (len = *),INTENT (in),optional :: name
        CHARACTER (len = 120)                 :: fname
        INTEGER                               :: ier,i,ios,ou,k
        integer :: il,jl,kl
        CHARACTER*28 :: sname='init_bulk2m_driver_gcm_mat: '
        dt = dt0
        if(PRESENT(il0)) then 
          il = il0
        else 
          il=im
        endif
        if(PRESENT(jl0)) then 
          jl = jl0
        else 
          jl=jm
        endif
        if(PRESENT(kl0)) then 
          kl = kl0
        else 
          kl=lm
        endif
        if(PRESENT(nm0)) then 
          nmodes = nm0
        else 
	nmodes=2
        endif
        if(PRESENT(name)) then 
          fname = trim(name)
        else 
          fname='blk_input.dat'
        endif
        if(il .eq. 1) ltermo=.true.
        lcheck=.false.
        ou = 6; ktau = 0; i=1; j = 1; rhw=0.95; rhi=0.95
c Allocate TERMO arrays
c        ldummy = allocate_arrays_bulk2m_driver(mx)
        ldummy = allocate_arrays_bulk2m_driver(kl)
c Hugh 
c constants and parameters
        ldummy=hugh_constants_init()
        ldummy=hugh_drop_nucleation_init()
        if(ldummy) then
          print *,sname,'Initialization has been completed...'
        else
          lcheck = .true.
          print *,sname,'BULK2M initialization has not been done ...'
          ldummy=cleanup_bulk2m_driver()
        endif
        la = .NOT.lcheck
      END FUNCTION init_bulk2m_driver_gcm_mat

C==============================================================================
C Include file: blk_add01.f
C========================================================================
      LOGICAL FUNCTION hugh_constants_init() RESULT (la)
      IMPLICIT NONE
! Need to make sure this works with the CONSTant values set in CONSTANT
!     USE CONSTANT, only: pi,rhow,r->rgas,rv->rvap,cp->sha,g->grav 
        CHARACTER*21 :: sname='hugh_constants_init: '
        la     = .true.
C Set constants
c initialize parameters for initial call to scheme
c fallspeed parameters (V=aD^b)
	   ai = 700.
	   ac = 3.e7
	   as = 11.72
	   ar = 841.99667
	   bi = 1.
	   bc = 2.
	   bs = 0.41
	   br = 0.8
c constants and parameters
! 	   pi = 4.0d0*datan(1.0d0)
 	   r = rgas
 	   rv = rvap
 	   cp = 1005.0d0
           rhosu = 85000./(rgas*tf)
 	   rhow = 997.
	   rhoi = 500.
	   rhosn = 100.
	   aimm = 0.66
	   bimm = 100.
	   ecr = 1.
	   dcs = 250.e-6
	   ri0 = 1.e-6
	   mi0 = 4./3.*pi*rhoi*(1.e-6)**3
	   mi0 = 4./3.*pi*rhoi*ri0**3
Cigs
	   rw0 = 1.e-6
	   mw0 = 4./3.*pi*rhow*rw0**3
Cigs
	   f1s = 0.86
	   f2s = 0.28 
	   f1r = 0.78
	   f2r = 0.32
!	   g = grav
	   qsmall = 1.e-20
	   nsmall = 1.e-20

c size distribution parameters
	   ci = rhoi*pi/6.0d0
	   di = 3.0d0
	   cs = rhosn*pi/6.0d0
	   ds = 3.0d0
	   dr = 1.0d0/3.0d0

c radius of contact nuclei
         rin = 0.1e-6
        la=.NOT.lcheck
      END FUNCTION hugh_constants_init
C==============================================================================
      logical function make_coefs(mx0) result(la)
      implicit none
        integer,            INTENT (in)       :: mx0
      real*8  :: dqsdt,dqsidt
c counting variables
      integer :: k
      character*12 :: sname='make_coefs: '
        la     = .true.
        do k=1,mx0
c atmospheric parameters that vary in time and height
          rho(k)=pp3d(k)/(r*tk3d(k))
c thermal conductivity for air
	      kap(k) = 1.414d3*1.496d-6*tk3d(k)**
     1          1.5d0/(tk3d(k)+120.0d0) 
c Saturation vapor pressure and mixing ratio
            evs(k) = polysvp(tk3d(k),0)   ! pa
            eis(k) = polysvp(tk3d(k),1)   ! pa
c make sure ice saturation does not exceed water sat. near freezing
	    if (eis(k).gt.evs(k)) eis(k) = evs(k)
c Saturation mixing ratio
            qvs(k) = .622d0*evs(k)/(pp3d(k)-evs(k))
            qvi(k) = .622d0*eis(k)/(pp3d(k)-eis(k))
            qvqvs(k)  = qv3d(k)/qvs(k)
            qvqvsi(k) = qv3d(k)/qvi(k)
c diffusivity of water vapor
	      Dv(k) = 8.794d-5*tk3d(k)**1.81/pp3d(k)
c latent heat of vaporation
            xxlv(k) = 3.1484d6-2370.0d0*tk3d(k)
c latent heat of sublimation
            xxls(k) = 3.15d6-2370.0d0*tk3d(k)+0.3337d6
c heat of fusion
	    xlf(k)  = xxls(k)-xxlv(k)
c viscosity of air
c schmit number
	      mu(k) = 1.496E-6*tk3d(k)**1.5/(tk3d(k)+120.)/rho(k)	
	      sc(k) = mu(k)/Dv(k)
c psychometic corrections
c rate of change sat. mix. ratio with temperature
	      dqsdt = xxlv(k)*qvs(k)/(rv*tk3d(k)**2)
             dqsidt = xxls(k)*qvi(k)/(rv*tk3d(k)**2)		
            abi(k) = 1.+dqsidt*xxls(k)/cp
            ab(k)  = 1.+dqsdt *xxlv(k)/cp
c..................................................................
c microphysics parameters varying in time/height
c fall speed with density correction (heymsfield and benssemer 2006)
	      ain(k) = (rhosu/rho(k))**0.54*ai
	      arn(k) = (rhosu/rho(k))**0.54*ar
	      asn(k) = (rhosu/rho(k))**0.54*as
	      acn(k) = (rhosu/rho(k))**0.54*ac
c mean free path
            fdum(k) = 7.37*tk3d(k)/(288.*10.*pp3d(k))/100.
c effective diffusivity of contact nuclei
c based on brownian diffusion
            dap(k) = 4.*pi*1.38e-23*tk3d(k)*(1.+fdum(k)/rin)/
     1                   (6.*pi*rin*mu(k))
c calculate collection efficiency as function of temperature
c loosely following Pruppacher and Klett 1997
        eii(k) = 0.01
        if (tk3d(k).ge.253.15 .and.tk3d(k).lt.263.15) then
           eii(k) = 0.4
        else if (tk3d(k).ge.263.15) then
           eii(k) = 0.7
        end if
        enddo
        la=.NOT.lcheck
      end function make_coefs
C==============================================================================
      logical function make_crys_dist(dt,ci,di,pi,qsmall,mx0) result(la)
      implicit none 
        real*8, intent(in)  :: dt,ci,di,pi,qsmall
        integer, intent(in) :: mx0
c local variables
      real*8 lammax,lammin,dum
c external function call
c	real*8 gamma
c counting variables
      integer i,j,k
      character*16 :: sname='make_crys_dist: '
        la     = .true.
c cloud ice
!      do i = 1,mix
!         do j = 1,mjx
       do k = 1,mx0
c If mixing ratio < 1.e-20 set mixing ratio and number conc to zero
	 if (qi3d(k).lt.qsmall) then
	   qi3d(k) = 0.
           ni3d(k) = 0.
	 end if
c make sure number concentrations are not negative
	ni3d(k) = max(0.0d0,ni3d(k))
c calculate size distribution parameters
	if (qi3d(k).ge.qsmall) then
	  lami(k) = (gamma(1.+di)*ci*
     1         ni3d(k)/qi3d(k))**(1./di)
	  n0i(k) = ni3d(k)*lami(k)
c check for slope
	lammax = 1.0d0/1.d-6
	lammin = 1.0d0/(2.0d0*dcs+100.0d-6)
c adjust vars
	if (lami(k).lt.lammin) then
	  lami(k) = lammin
	  n0i(k) = lami(k)**4*qi3d(k)/(ci*gamma(1.0d0+di))
c      ni3dten(k) = ni3dten(k)+
c     1       (n0i(k)/lami(k)-ni3d(k))/dtmic
	  ni3d(k) = n0i(k)/lami(k)
	else if (lami(k).gt.lammax) then
	  lami(k) = lammax
  	  n0i(k) = lami(k)**4*qi3d(k)/(ci*gamma(1.0d0+di))
c        ni3dten(k) = ni3dten(k)+
c     1       (n0i(k)/lami(k)-ni3d(k))/dtmic
	  ni3d(k) = n0i(k)/lami(k)
	end if
	else
	  lami(k) = 0.0d0
	  n0i(k) = 0.0d0
	end if
c calculate effective radius: default value is 25 micron
c assume spherical ice habit

	if (qi3d(k).ge.qsmall) then
	   effi(k) = 3.0d0/lami(k)/2.0*1.d6
	else
	   effi(k) = 25.0d0
	end if
c put limits of values for pass to radiation code
        effi(k) = max(effi(k),2.0d0)
        effi(k) = min(effi(k),150.0d0)
       end do  ! k loop
!          end do ! j loop
!        end do  ! i loop
        la=.NOT.lcheck
      end function make_crys_dist
C
C==============================================================================
      logical function make_snow_dist(dt,cs,ds,qsmall,mx0) result(la)
      implicit none 
        real*8,  intent(in) :: dt,cs,ds,qsmall
        integer, intent(in) :: mx0
c local variables
      real*8 lammax,lammin,dum
c external function call
c	real*8 gamma
c counting variables
      integer i,j,k
      character*16 :: sname='make_snow_dist: '
        la     = .true.
c snow
!      do i = 1,mix
!         do j = 1,mjx
       do k = 1,mx0
c If mixing ratio < 1.e-20 set mixing ratio and number conc to zero
	 if (qs3d(k).lt.qsmall) then
         qs3d(k) = 0.
         ns3d(k) = 0.
         end if
c make sure number concentrations are not negative
	ns3d(k) = max(0.0d0,ns3d(k))
c calculate size distribution parameters
	if (qs3d(k).ge.qsmall) then
	  lams(k) = (gamma(1.0d0+ds)*cs*ns3d(k)/
     1       qs3d(k))**(1.0d0/ds)
	  n0s(k) = ns3d(k)*lams(k)
c check for slope
	lammax = 1.0d0/5.0d-6
	lammin = 1.0d0/5000.0d-6
c adjust vars
	if (lams(k).lt.lammin) then
	  lams(k) = lammin
	  n0s(k) = lams(k)**4*qs3d(k)/(cs*gamma(1.0d0+ds))
c        ns3dten(k) = ns3dten(k)+
c     1      (n0s(k)/lams(k)-ns3d(k))/dtmic
	  ns3d(k) = n0s(k)/lams(k)
	else if (lams(k).gt.lammax) then
	  lams(k) = lammax
	  n0s(k) = lams(k)**4*qs3d(k)/(cs*gamma(1.0d0+ds))
c        ns3dten(k) = ns3dten(k)+
c     1      (n0s(k)/lams(k)-ns3d(k))/dtmic
	  ns3d(k) = n0s(k)/lams(k)
	end if
	else
	  lams(k) = 0.0d0
	  n0s(k) = 0.0d0
	end if
c calculate effective radius: default value is 25 micron
c assume spherical snow habit
	if (qs3d(k).ge.qsmall) then
	   effs(k) = 3.0d0/lams(k)/2.0*1.d6
	else
	   effs(k) = 25.0d0
	end if	
       end do  ! k loop
!          end do ! j loop
!        end do  ! i loop
        la=.NOT.lcheck
      end function make_snow_dist
C==============================================================================
      logical function make_rain_dist(dt,dr,pi,rhow,qsmall
     *  ,mx0,ind) result(la)
      implicit none 
        real*8, intent(in)  :: dt,dr,pi,rhow,qsmall
        integer, intent(in) :: mx0,ind
c local variables
      real*8 lammax,lammin,dum
c external function call
c	real*8 gamma
c counting variables
      integer i,j,k
      character*16 :: sname='make_rain_dist: '
        la     = .true.
c rain
!      do i = 1,mix
!         do j = 1,mjx
       do k = 1,mx0
c If mixing ratio < 1.e-20 set mixing ratio and number conc to zero
	 if (qr3d(k).lt.qsmall) then
	   qr3d(k) = 0.
           nr3d(k) = 0.
	 end if
c make sure number concentrations are not negative
        nr3d(k) = max(0.0d0,nr3d(k))
c calculate lamr
	if (qr3d(k).ge.qsmall) then
	  lamr(k) = (pi*rhow*nr3d(k)/qr3d(k))**(1.0d0/3.0d0)
	  n0rr(k) = nr3d(k)*lamr(k)
c check for slope
	lammax = 1.0d0/5.d-6
	lammin = 1.0d0/5000.d-6
c adjust vars
	if (lamr(k).lt.lammin) then
	  lamr(k) = lammin
	  n0rr(k) = lamr(k)**4*qr3d(k)/(pi*rhow)
c        nr3dten(k) = nr3dten(k)+
c     1       (n0rr(k)/lamr(k)-nr3d(k))/dtmic
	  nr3d(k) = n0rr(k)/lamr(k)
	else if (lamr(k).gt.lammax) then
	  lamr(k) = lammax
	  n0rr(k) = lamr(k)**4*qr3d(k)/(pi*rhow)
c        nr3dten(k) = nr3dten(k)+
c     1       (n0rr(k)/lamr(k)-nr3d(k))/dtmic
   	  nr3d(k) = n0rr(k)/lamr(k)
	end if
	else
	  lamr(k) = 0.0d0
	  n0rr(k) = 0.0d0
	end if
c calculate effective radius: default value is 25 micron
	if (qr3d(k).ge.qsmall) then
	   effr(k) = 3.0d0/lamr(k)/2.0*1.d6
	else
	   effr(k) = 25.0d0
	end if	
       end do  ! k loop
!          end do ! j loop
!        end do  ! i loop
        la=.NOT.lcheck
      end function make_rain_dist
C==============================================================================
      logical function make_drop_dist(dt,pi,rhow,tk,pk,qsmall,mx0) 
     *  result(la)
      implicit none 
        integer, intent(in) :: mx0
        real*8,  intent(in) :: tk(mx0),pk(mx0)
        real*8,  intent(in) :: dt,pi,rhow,qsmall
!       integer, intent(in) :: mx0
c local variables
      real*8 lammax,lammin,dum
c external function call
c	real*8 gamma
c counting variables
      integer i,j,k
      character*16 :: sname='make_drop_dist: '
        la     = .true.
c cloud droplets
!      do i = 1,mix
!         do j = 1,mjx
       do k = 1,mx0
c If mixing ratio < 1.e-20 set mixing ratio and number conc to zero
	 if (qc3d(k).lt.qsmall) then
	   qc3d(k) = 0.
           nc3d(k) = 0.
         end if
c make sure number concentrations are not negative
	 nc3d(k) = max(0.0d0,nc3d(k))
c Martin et al. (1994) formula for pgam
         pgam(k) = 0.0d0
	if (qc3d(k).ge.qsmall) then
         dum = pk(k)/(287.15*tk(k))
         pgam(k)=0.0005714*(nc3d(k)/1.0d6/dum)+0.2714d0
         pgam(k)=1.0d0/(pgam(k)**2)-1.0d0
         pgam(k)=max(pgam(k), 2.0d0)
         pgam(k)=min(pgam(k),15.0d0)
c calculate lamc
	lamc(k) = (pi/6.*rhow*nc3d(k)*gamma(pgam(k)+4.0d0)/
     1            (qc3d(k)*gamma(pgam(k)+1.0d0)))**(1./3.)
c       print*,"Here",lamc(k),k 
c
c lammin, 60 micron diameter
c lammax, 1 micron
	lammin = (pgam(k)+1.)/60.0d-6
	lammax = (pgam(k)+1.)/ 1.0d-6
	if (lamc(k).lt.lammin) then
	  lamc(k) = lammin
c        nc3dten(k) = nc3dten(k)+
c     1     (6.0d0*lamc(k)**3*qc3d(k)*
c     1           gamma(pgam(k)+1.0d0)/
c     1          (pi*rhow*gamma(pgam(k)+4.0d0))-nc3d(k))/dtmic
	  nc3d(k) = 6.0d0*lamc(k)**3*qc3d(k)*
     1           gamma(pgam(k)+1.0d0)/
     1          (pi*rhow*gamma(pgam(k)+4.0d0))
	else if (lamc(k).gt.lammax) then
	  lamc(k) = lammax
c        nc3dten(k) = nc3dten(k)+
c     1     (6.0d0*lamc(k)**3*qc3d(k)*
c     1           gamma(pgam(k)+1.0d0)/
c     1          (pi*rhow*gamma(pgam(k)+4.0d0))-nc3d(k))/dtmic
  	         nc3d(k) = 6.0d0*lamc(k)**3*qc3d(k)*
     1             gamma(pgam(k)+1.0d0)/
     1               (pi*rhow*gamma(pgam(k)+4.0d0))
	       end if
c to calculate droplet freezing
                 cdist1(k) = nc3d(k)/gamma(pgam(k)+1.0d0) 
	       else
	         lamc(k) = 0.0d0
                 cdist1(k) = 0.0d0
               end if
c calculate drop effective radius default value is 25 micron
	       if (qc3d(k).ge.qsmall) then
	         effc(k) = gamma(pgam(k)+4.0d0)/
     1             gamma(pgam(k)+3.0d0)/lamc(k)/2.0*1.d6
	       else
                  effc(k) = 25.0d0
	       end if
c put limits of values for pass to radiation code
               effc(k) = max(effc(k),2.0d0)
               effc(k) = min(effc(k),40.0d0)
       end do  ! k loop
!          end do ! j loop
!        end do  ! i loop
        la=.NOT.lcheck
      end function make_drop_dist
C
C==============================================================================
      logical function make_xxxx_dist(dt0,ci,di,pi,rhow,cs,ds,dr
     *  ,tk,pk,qsmall,tag,mx0) result(la)
      implicit none 
        integer, intent(in) :: mx0
        real*8,  intent(in) :: tk(mx0),pk(mx0)
        real*8,  intent(in) :: dt0,ci,di,pi,rhow,cs,ds,dr,qsmall
        character (len = *),INTENT (in)       :: tag
c       integer, intent(in) :: mx0
c local variables
      real*8 lammax,lammin,dum,dtmic
c external function call
c	real*8 gamma
c counting variables
      integer      :: i,j,k
c dummy integer
      integer      :: ind=0
      character*16 :: sname='make_xxxx_dist: '
        la     = .true.
        dtmic=dt0
        case_tag: select case(trim(tag))
c Make droplet distribution
          case('drop')
            la= make_drop_dist(dtmic,pi,rhow,tk,pk,qsmall,mx0) 
c Make rain distribution
          case('rain')
            la= make_rain_dist(dtmic,dr,pi,rhow,qsmall,mx0,ind) 
c Make crystal distribution
          case('crys')
            la= make_crys_dist(dtmic,ci,di,pi,qsmall,mx0) 
c Make rain distribution
          case('snow')
            la= make_snow_dist(dtmic,cs,ds,qsmall,mx0)
c Make liquid phase distributions
          case('liquid')
            la= make_drop_dist(dtmic,pi,rhow,tk,pk,qsmall,mx0) 
            la= make_rain_dist(dtmic,dr,pi,rhow,qsmall,mx0,ind) 
c Make solid phase distributions
          case('solid')
            la= make_crys_dist(dtmic,ci,di,pi,qsmall,mx0) 
            la= make_snow_dist(dtmic,cs,ds,qsmall,mx0)
c Make all distributions
          case('all')
            la= make_drop_dist(dtmic,pi,rhow,tk,pk,qsmall,mx0) 
            la= make_rain_dist(dtmic,dr,pi,rhow,qsmall,mx0,ind) 
            la= make_crys_dist(dtmic,ci,di,pi,qsmall,mx0) 
            la= make_snow_dist(dtmic,cs,ds,qsmall,mx0)
          case default
            msg="Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select case_tag
        la=.NOT.lcheck
      end function make_xxxx_dist
C
C==============================================================================
      logical function make_distribution0(tk,pk,dt,ci,di,pi,rhow
     *  ,cs,ds,dr,qsmall,mx0) result (la)
      implicit none 
        integer, intent(in) :: mx0
        real*8,  intent(in) :: tk(mx0),pk(mx0)
        real*8,  intent(in) :: dt,ci,di,pi,rhow,cs,ds,dr,qsmall
c       integer, intent(in) :: mx0
c local variables
      real*8 lammax,lammin,dum
c external function call
c	real*8 gamma
c counting variables
      integer i,j,k
      character*20 :: sname='make_distribution0: '
c calculate size distribution parameters
c..................................................................
!      do i = 1,mix
!         do j = 1,mjx
            do k = 1,mx0
               pgam(k) = 0.0
c calculate size distribution parameters
c..................................................................
c If mixing ratio < 1.e-20 set mixing ratio and number conc to zero
	 if (qi3d(k).lt.qsmall) then
	   qi3d(k) = 0.
           ni3d(k) = 0.
	 end if
	 if (qs3d(k).lt.qsmall) then
         qs3d(k) = 0.
         ns3d(k) = 0.
       end if
	 if (qc3d(k).lt.qsmall) then
	   qc3d(k) = 0.
           nc3d(k) = 0.
       end if
	 if (qr3d(k).lt.qsmall) then
	   qr3d(k) = 0.
           nr3d(k) = 0.
	 end if
c make sure number concentrations are not negative
	ni3d(k) = max(0.0d0,ni3d(k))
	ns3d(k) = max(0.0d0,ns3d(k))
	nc3d(k) = max(0.0d0,nc3d(k))
        nr3d(k) = max(0.0d0,nr3d(k))
c......................................................................
c cloud ice
	if (qi3d(k).ge.qsmall) then
	  lami(k) = (gamma(1.+di)*ci*
     1         ni3d(k)/qi3d(k))**(1./di)
	  n0i(k) = ni3d(k)*lami(k)
c check for slope
	lammax = 1.0d0/1.d-6
	lammin = 1.0d0/(2.0d0*dcs+100.0d-6)
c adjust vars
	if (lami(k).lt.lammin) then
	  lami(k) = lammin
	  n0i(k) = lami(k)**4*qi3d(k)/(ci*gamma(1.0d0+di))
c      ni3dten(k) = ni3dten(k)+
c     1       (n0i(k)/lami(k)-ni3d(k))/dtmic
	  ni3d(k) = n0i(k)/lami(k)
	else if (lami(k).gt.lammax) then
	  lami(k) = lammax
  	  n0i(k) = lami(k)**4*qi3d(k)/(ci*gamma(1.0d0+di))
c        ni3dten(k) = ni3dten(k)+
c     1       (n0i(k)/lami(k)-ni3d(k))/dtmic
	  ni3d(k) = n0i(k)/lami(k)
	end if
	else
	  lami(k) = 0.0d0
	  n0i(k) = 0.0d0
	end if
c......................................................................
c rain
	if (qr3d(k).ge.qsmall) then
	  lamr(k) = (pi*rhow*nr3d(k)/qr3d(k))**(1.0d0/3.0d0)
	  n0rr(k) = nr3d(k)*lamr(k)
c check for slope
	lammax = 1.0d0/5.d-6
	lammin = 1.0d0/5000.d-6
c adjust vars
	if (lamr(k).lt.lammin) then
	  lamr(k) = lammin
	  n0rr(k) = lamr(k)**4*qr3d(k)/(pi*rhow)
c        nr3dten(k) = nr3dten(k)+
c     1       (n0rr(k)/lamr(k)-nr3d(k))/dtmic

	  nr3d(k) = n0rr(k)/lamr(k)
	else if (lamr(k).gt.lammax) then
	  lamr(k) = lammax
	  n0rr(k) = lamr(k)**4*qr3d(k)/(pi*rhow)
c        nr3dten(k) = nr3dten(k)+
c     1       (n0rr(k)/lamr(k)-nr3d(k))/dtmic
   	  nr3d(k) = n0rr(k)/lamr(k)
	end if
	else
	  lamr(k) = 0.0d0
	  n0rr(k) = 0.0d0
	end if
c......................................................................
c cloud droplets
c Martin et al. (1994) formula for pgam
	if (qc3d(k).ge.qsmall) then
         dum = pp3d(k)/(287.15*tk3d(k))
         pgam(k)=0.0005714*(nc3d(k)/1.0d6/dum)+0.2714d0
         pgam(k)=1.0d0/(pgam(k)**2)-1.0d0
         pgam(k)=max(pgam(k), 2.0d0)
         pgam(k)=min(pgam(k),15.0d0)
c calculate lamc
	lamc(k) = (pi/6.*rhow*nc3d(k)*gamma(pgam(k)+4.0d0)/
     1            (qc3d(k)*gamma(pgam(k)+1.0d0)))**(1./3.)
c lammin, 60 micron diameter
c lammax, 1 micron
	lammin = (pgam(k)+1.)/60.0d-6
	lammax = (pgam(k)+1.)/ 1.0d-6
	if (lamc(k).lt.lammin) then
	  lamc(k) = lammin
c        nc3dten(k) = nc3dten(k)+
c     1     (6.0d0*lamc(k)**3*qc3d(k)*
c     1           gamma(pgam(k)+1.0d0)/
c     1          (pi*rhow*gamma(pgam(k)+4.0d0))-nc3d(k))/dtmic
	  nc3d(k) = 6.0d0*lamc(k)**3*qc3d(k)*
     1           gamma(pgam(k)+1.0d0)/
     1          (pi*rhow*gamma(pgam(k)+4.0d0))
	else if (lamc(k).gt.lammax) then
	  lamc(k) = lammax
c        nc3dten(k) = nc3dten(k)+
c     1     (6.0d0*lamc(k)**3*qc3d(k)*
c     1           gamma(pgam(k)+1.0d0)/
c     1          (pi*rhow*gamma(pgam(k)+4.0d0))-nc3d(k))/dtmic
  	  nc3d(k) = 6.0d0*lamc(k)**3*qc3d(k)*
     1           gamma(pgam(k)+1.0d0)/
     1          (pi*rhow*gamma(pgam(k)+4.0d0))
	end if
c to calculate droplet freezing
        cdist1(k) = nc3d(k)/gamma(pgam(k)+1.0d0) 
	else
	  lamc(k) = 0.0d0
          cdist1(k) = 0.0d0
	end if

c......................................................................
c snow
	if (qs3d(k).ge.qsmall) then
	  lams(k) = (gamma(1.0d0+ds)*cs*ns3d(k)/
     1       qs3d(k))**(1.0d0/ds)
	  n0s(k) = ns3d(k)*lams(k)
c check for slope
	lammax = 1.0d0/5.0d-6
	lammin = 1.0d0/5000.0d-6
c adjust vars
	if (lams(k).lt.lammin) then
	  lams(k) = lammin
	  n0s(k) = lams(k)**4*qs3d(k)/(cs*gamma(1.0d0+ds))
c        ns3dten(k) = ns3dten(k)+
c     1      (n0s(k)/lams(k)-ns3d(k))/dtmic
	  ns3d(k) = n0s(k)/lams(k)
	else if (lams(k).gt.lammax) then
	  lams(k) = lammax
	  n0s(k) = lams(k)**4*qs3d(k)/(cs*gamma(1.0d0+ds))
c        ns3dten(k) = ns3dten(k)+
c     1      (n0s(k)/lams(k)-ns3d(k))/dtmic
	  ns3d(k) = n0s(k)/lams(k)
	end if
	else
	  lams(k) = 0.0d0
	  n0s(k) = 0.0d0
	end if
c......................................................................
c calculate effective radius
c assume spherical ice habit
c default value (if no cloud/precip water is present) is 25 micron
	if (qi3d(k).ge.qsmall) then
	   effi(k) = 3.0d0/lami(k)/2.0*1.d6
	else
	   effi(k) = 25.0d0
	end if

	if (qs3d(k).ge.qsmall) then
	   effs(k) = 3.0d0/lams(k)/2.0*1.d6
	else
	   effs(k) = 25.0d0
	end if	

	if (qr3d(k).ge.qsmall) then
	   effr(k) = 3.0d0/lamr(k)/2.0*1.d6
	else
	   effr(k) = 25.0d0
	end if	
	
	if (qc3d(k).ge.qsmall) then
	effc(k) = gamma(pgam(k)+4.0d0)/
     1        gamma(pgam(k)+3.0d0)/lamc(k)/2.0*1.d6
	else
          effc(k) = 25.0d0
	end if

c put limits of values for pass to radiation code
        effc(k) = max(effc(k),2.0d0)
        effi(k) = max(effi(k),2.0d0)
        effc(k) = min(effc(k),40.0d0)
        effi(k) = min(effi(k),150.0d0)

            end do  ! k loop
!          end do ! j loop
!        end do  ! i loop

        la=.NOT.lcheck
       end function make_distribution0
C
C==============================================================================
      LOGICAL FUNCTION execute_bulk2m_driver_gcm(tag,action,dt0,k00,
     *r1,r2)     
     *  RESULT (la)
      IMPLICIT NONE
        CHARACTER (len = *),INTENT (in)       :: tag,action
        REAL*8,             INTENT (in)       :: dt0
        REAL*8,optional,    INTENT (in)       :: r1,r2  
        INTEGER,optional,   INTENT (in)       :: k00
c Local
        integer      :: k0
        CHARACTER*26 :: sname='execute_bulk2m_driver_gcm: '
        la     = .true.
        if(PRESENT(k00)) then
          k0=k00
        else
          k0=1
        endif
C
C**********This is to use the old Gultepe/Isaac scheme for CDNC
        case_tag: select case(trim(tag))
          case ('gult')
            case_action3: select case (trim(action))
              case ("drop_nucl")
                ldummy=gult_drop_activation(dt0,k0,r1,r2)
              case default
                msg=trim(tag)//
     *             " Case "//trim(action)//" is not implemented"
                call stop_model(msg,255)
              end select case_action3
C************************************************************************************
C**********This is to use the MATRIX activated fraction 
c       case_tag: select case(trim(tag))
          case ('matr')
            case_action4: select case (trim(action))
              case ("drop_nucl")
                ldummy=matr_drop_activation(dt0,k0,.false.)
              case default
                msg=trim(tag)//
     *             " Case "//trim(action)//" is not implemented"
                call stop_model(msg,255)
              end select case_action4
C************************************************************************************
C**********This is to use the TOMAS activated fraction 
c       case_tag: select case(trim(tag))
          case ('toma')
            case_action5: select case (trim(action))
              case ("drop_nucl")
                ldummy=toma_drop_activation(dt0,k0,.false.)
              case default
                msg=trim(tag)//
     *             " Case "//trim(action)//" is not implemented"
                call stop_model(msg,255)
              end select case_action5
C************************************************************************************
          case ('surabi')
            case_action1: select case (trim(action))
              case ("GET_CDNC_UPD")
                ldummy=surabi_drop_activation(dt0,k0)
              case ("drop_nucl")
                ldummy=hugh_drop_nucleation(dt0,k0)
              case ("drop_auto")
                ldummy=hugh_drop_autoconversion(dt0,k0)
              case ("drop_rain")
                ldummy= hugh_drop_by_rain_accretion(dt0,k0)
              case ("drop_frzn")
                ldummy=hugh_drop_freezing(dt0,k0)
              case ("rain_frzn")
                ldummy=hugh_rain_freezing(dt0,k0)
              case ("crys_nucl")
                ldummy=hugh_crys_nucleation(dt0,k0)
              case ("crys_auto")
                ldummy=hugh_crys_autoconversion(dt0,k0)
              case ("crys_snow")
                ldummy= hugh_crys_by_snow_accretion(dt0,k0)
              case default
                msg=trim(tag)//
     *             " Case "//trim(action)//" is not implemented"
                call stop_model(msg,255)
              end select case_action1
          case ('hugh')
            case_action2: select case (trim(action))
              case ("drop_nucl")
                ldummy=hugh_drop_nucleation(dt0,k0)
              case ("drop_auto")
                ldummy=hugh_drop_autoconversion(dt0,k0)
              case ("drop_rain")
                ldummy= hugh_drop_by_rain_accretion(dt0,k0)
              case ("drop_snow")
                ldummy= hugh_drop_by_snow_accretion(dt0,k0)
              case ("drop_cond")
                ldummy=hugh_drop_condensation(dt0,k0)
              case ("rain_snow")
                ldummy= hugh_rain_by_snow_accretion(dt0,k0)
              case ("drop_frzn")
                ldummy=hugh_drop_freezing(dt0,k0)
              case ("drop_sedim")
                ldummy=hugh_drop_sedimentation(dt0,k0)
              case ("rain_frzn")
                ldummy=hugh_rain_freezing(dt0,k0)
              case ("crys_nucl")
                ldummy=hugh_crys_nucleation(dt0,k0)
              case ("crys_cond")
                ldummy=hugh_crys_condensation(dt0,k0)
              case ("crys_auto")
                ldummy=hugh_crys_autoconversion(dt0,k0)
              case ("crys_snow")
                ldummy= hugh_crys_by_snow_accretion(dt0,k0)
              case ("crys_melt")
                ldummy=hugh_crys_melting(dt0,k0)
              case ("snow_cond")
                ldummy=hugh_crys_condensation(dt0,k0)
              case ("snow_melt")
                ldummy=hugh_snow_melting(dt0,k0)
              case default
                msg=trim(tag)//
     *             " Case "//trim(action)//" is not implemented"
                call stop_model(msg,255)
              end select case_action2
          case default
            msg="Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select case_tag
C
        la=.NOT.lcheck
      END FUNCTION execute_bulk2m_driver_gcm
C
C==============================================================================
C Include file: blk_add_execute_interface.f
C========================================================================
C
      LOGICAL FUNCTION update_micro_tendencies_bulk2m(tag,tag0,k0) 
     *  RESULT (la)
c     *  ,nc0,qc0,ni0,qi0,tag0,nr0,qr0,ns0,qs0) RESULT (la)
      IMPLICIT NONE
        CHARACTER (len = *),INTENT (in)       :: tag,tag0
        integer,            intent(in)        :: k0
c Local
        integer      :: mx0!,k0
        character*32 :: sname='update_micro_tendencies_bulk2m: '
        real*8                :: dtmic
        la     = .true.
        mx0=k0
        dtmic = dt
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
c          ldummy=hugh_balance_tendencies(dtmic,mx0)
        which_tendency: select case (trim(tag))
          case ('tkqv')
            action_tkqv: select case (trim(tag0))
              case ('tk_qv')
c                ldummy=hugh_balance_tendencies(dtmic,mx0)
                ldummy=hugh_update_tkqv_tendencies0(dtmic,mx0)
                ldummy=hugh_update_tkqv_tendencies(dtmic,mx0)
              case ('2zero')
                ldummy=hugh_update_tkqv_2zero(dtmic,mx0)
              case default
                msg="Case "//trim(tag)//":"//trim(tag0)
     *            //" is not implemented"
                call stop_model(msg,255)             
            end select action_tkqv
c
          case ('drop')
            action_drop: select case (trim(tag0))
              case ('01') 
                ldummy=hugh_update_drop_tendencies0(dtmic,mx0)
              case ('2zero')
                ldummy=hugh_update_drop_2zero(dtmic,mx0)
              case default
                msg="Case "//trim(tag)//":"//trim(tag0)
     *            //" is not implemented"
                call stop_model(msg,255)             
            end select action_drop
c
          case ('crys')
            action_crys: select case (trim(tag0))
              case ('01') 
                ldummy=hugh_update_crys_tendencies0(dtmic,mx0)
              case ('2zero')
                ldummy=hugh_update_crys_2zero(dtmic,mx0)
              case default
                msg="Case "//trim(tag)//":"//trim(tag0)
     *            //" is not implemented"
                call stop_model(msg,255)             
            end select action_crys
c
          case ('rain')
            action_rain: select case (trim(tag0))
              case ('01') 
                ldummy=hugh_update_rain_tendencies0(dtmic,mx0)
              case ('2zero')
                ldummy=hugh_update_rain_2zero(dtmic,mx0)
              case default
                msg="Case "//trim(tag)//":"//trim(tag0)
     *            //" is not implemented"
                call stop_model(msg,255)             
            end select action_rain
c
          case ('snow')
            action_snow: select case (trim(tag0))
              case ('01') 
                ldummy=hugh_update_snow_tendencies0(dtmic,mx0)
              case ('2zero')
                ldummy=hugh_update_snow_2zero(dtmic,mx0)
              case default
                msg="Case "//trim(tag)//":"//trim(tag0)
     *            //" is not implemented"
                call stop_model(msg,255)             
            end select action_snow
c
          case ('all')
            action_all: select case (trim(tag0))
              case ('make_balance') 
                ldummy=hugh_balance_tendencies(dtmic,mx0)
                ldummy=srbm_balance_tendencies(dtmic,mx0)
                ldummy=hugh_update_drop_tendencies0(dtmic,mx0)
                ldummy=hugh_update_crys_tendencies0(dtmic,mx0)
                ldummy=hugh_update_rain_tendencies0(dtmic,mx0)
                ldummy=hugh_update_snow_tendencies0(dtmic,mx0)
                ldummy=hugh_update_tkqv_tendencies0(dtmic,mx0)
                ldummy=hugh_update_tkqv_tendencies(dtmic,mx0)
              case ('02') 
                ldummy=hugh_update_tkqv_tendencies(dtmic,mx0)
              case ('2zero')
                ldummy=hugh_update_all_2zero(dtmic,mx0)
              case default
                msg="Case "//trim(tag)//":"//trim(tag0)
     *            //" is not implemented"
            end select action_all
c
          case default
            msg=sname//"Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select which_tendency
c 
        la=.NOT.lcheck
      END FUNCTION update_micro_tendencies_bulk2m
C
C==============================================================================
      LOGICAL FUNCTION set_micro_bulk2m(tag,nc0,qc0,ni0,qi0,tag0
     *  ,nr0,qr0,ns0,qs0) RESULT (la)
      IMPLICIT NONE
        CHARACTER (len = *),INTENT (in)       :: tag,tag0
        real*8,dimension(:),intent(in)          :: nc0,qc0,ni0,qi0
        real*8,dimension(:),intent(in),optional :: nr0,qr0,ns0,qs0
c Local
        integer      :: mx0
        character*18 :: sname='set_micro_bulk2m: '
        real*8, dimension(size(nc0,1)) :: nr00,qr00,ns00,qs00
        real*8                :: dtmic
        la     = .true.
        mx0=SIZE(nc0,1)
c        stop 890
        dtmic = dt
        if(PRESENT(nr0)) then
          nr00=nr0/rho
        else
          nr00=0.0d0
        endif
        if(PRESENT(qr0)) then
          qr00=qr0
        else
          qr00=0.0d0
        endif
        if(PRESENT(ns0)) then
          ns00=ns0/rho
        else
          ns00=0.0d0
        endif
        if(PRESENT(qs0)) then
          qs00=qs0
        else
          qs00=0.0d0
        endif
        lcold =  tk3d.le.tf
        lwarm =  .NOT. lcold
        select case (trim(tag))
           case ('drop')
             nc3d=nc0/rho;qc3d=qc0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'drop',mx0)
             lqc3d = qc3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'drop',mx0)
           case ('crys')
             ni3d=ni0/rho;qi3d=qi0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'crys',mx0)
             lqi3d = qi3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'crys',mx0)
           case ('all')
             nc3d=nc0/rho;qc3d=qc0;ni3d=ni0/rho;qi3d=qi0
             nr3d=nr00/rho;qr3d=qr00;ns3d=ns00/rho;qs3d=qs00
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'all',mx0)
             lqc3d = qc3d.ge.qsmall
             lqi3d = qi3d.ge.qsmall
             lqs3d = qs3d.ge.qsmall
             lqr3d = qr3d.ge.qsmall
           case ('nc')
             nc3d=nc0/rho
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'drop',mx0)
             lqc3d = qc3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'drop',mx0)
           case ('qc')
             qc3d=qc0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'drop',mx0)
             lqc3d = qc3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'drop',mx0)
           case ('ni')
             ni3d=ni0/rho
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'crys',mx0)
             lqi3d = qi3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'crys',mx0)
           case ('qi')
             qi3d=qi0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'crys',mx0)
             lqi3d = qi3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'crys',mx0)
           case ('nr'); nr3d=nr00/rho
           case ('qr'); qr3d=qr00
           case ('ns'); ns3d=ns00/rho
           case ('qs'); qs3d=qs00
          case default
            msg=sname//"Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select
c
c Make supsat timescale 
        ldummy=make_xxxx_timescale(dtmic,'all',mx0)
        la=.NOT.lcheck
      END FUNCTION set_micro_bulk2m
C
C==============================================================================
      LOGICAL FUNCTION set_micro_bulk2m_mat(tag,nc0,qc0,ni0,qi0 
     * ,na0,nm0,tag0
     * ,nr0,qr0,ns0,qs0) RESULT (la)
      IMPLICIT NONE
        CHARACTER (len = *),INTENT (in)       :: tag,tag0
        real*8,dimension(:),intent(in)          :: nc0,qc0,ni0,qi0
        real*8,dimension(:,:),intent(in)        :: na0
        real*8,dimension(:),intent(in),optional :: nr0,qr0,ns0,qs0
        integer      , INTENT (in)   ::nm0 
c Local
        integer      :: mx0
        character*22 :: sname='set_micro_bulk2m_mat: '
        real*8, dimension(size(nc0,1)) :: nr00,qr00,ns00,qs00
        real*8                :: dtmic
        la     = .true.
        mx0=SIZE(nc0,1)
c        stop 890
        dtmic = dt

        ncactv(1:mx0) = sum(na0(1:mx0,:)) 
c       write(6,*)"IN MODES",nm0
c       do i =1,nm0
c        write(6,*)"MODES",na0(:,i),nm0
c       enddo
c       write(6,*)"NUMBER",ncactv(1:mx0),mx0
        if(PRESENT(nr0)) then
          nr00=nr0/rho
        else
          nr00=0.0d0
        endif
        if(PRESENT(qr0)) then
          qr00=qr0
        else
          qr00=0.0d0
        endif
        if(PRESENT(ns0)) then
          ns00=ns0/rho
        else
          ns00=0.0d0
        endif
        if(PRESENT(qs0)) then
          qs00=qs0
        else
          qs00=0.0d0
        endif
        lcold =  tk3d.le.tf
        lwarm =  .NOT. lcold
        select case (trim(tag))
           case ('drop')
             nc3d=nc0/rho;qc3d=qc0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'drop',mx0)
             lqc3d = qc3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'drop',mx0)
           case ('crys')
             ni3d=ni0/rho;qi3d=qi0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'crys',mx0)
             lqi3d = qi3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'crys',mx0)
           case ('all')
             nc3d=nc0/rho;qc3d=qc0;ni3d=ni0/rho;qi3d=qi0
             nr3d=nr00/rho;qr3d=qr00;ns3d=ns00/rho;qs3d=qs00
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'all',mx0)
             lqc3d = qc3d.ge.qsmall
             lqi3d = qi3d.ge.qsmall
             lqs3d = qs3d.ge.qsmall
             lqr3d = qr3d.ge.qsmall
           case ('nc')
             nc3d=nc0/rho
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'drop',mx0)
             lqc3d = qc3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'drop',mx0)
           case ('qc')
             qc3d=qc0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'drop',mx0)
             lqc3d = qc3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'drop',mx0)
           case ('ni')
             ni3d=ni0/rho
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'crys',mx0)
             lqi3d = qi3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'crys',mx0)
           case ('qi')
             qi3d=qi0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'crys',mx0)
             lqi3d = qi3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'crys',mx0)
           case ('nr'); nr3d=nr00/rho
           case ('qr'); qr3d=qr00
           case ('ns'); ns3d=ns00/rho
           case ('qs'); qs3d=qs00
          case default
            msg=sname//"Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select
c
c Make supsat timescale 
        ldummy=make_xxxx_timescale(dtmic,'all',mx0)
        la=.NOT.lcheck
      END FUNCTION set_micro_bulk2m_mat
      LOGICAL FUNCTION set_micro_bulk2m_tomas(tag,nc0,qc0,ni0,qi0 
     * ,na0,tag0
     * ,nr0,qr0,ns0,qs0) RESULT (la)
      IMPLICIT NONE
        CHARACTER (len = *),INTENT (in)       :: tag,tag0
        real*8,dimension(:),intent(in)          :: nc0,qc0,ni0,qi0
        real*8,dimension(:),intent(in)        :: na0
        real*8,dimension(:),intent(in),optional :: nr0,qr0,ns0,qs0
!        integer      , INTENT (in)   ::nm0 
c Local
        integer      :: mx0
        character*22 :: sname='set_micro_bulk2m_tom: '
        real*8, dimension(size(nc0,1)) :: nr00,qr00,ns00,qs00
        real*8                :: dtmic
        la     = .true.
        mx0=SIZE(nc0,1)
c        stop 890
        dtmic = dt

        ncactv(1:mx0) = na0(1:mx0)
c       write(6,*)"IN MODES",nm0
!         print*,'bulk2m_tomas',mx0,na0(mx0),ncactv(1:mx0)
c       write(6,*)"NUMBER",ncactv(1:mx0),mx0
        if(PRESENT(nr0)) then
          nr00=nr0/rho
        else
          nr00=0.0d0
        endif
        if(PRESENT(qr0)) then
          qr00=qr0
        else
          qr00=0.0d0
        endif
        if(PRESENT(ns0)) then
          ns00=ns0/rho
        else
          ns00=0.0d0
        endif
        if(PRESENT(qs0)) then
          qs00=qs0
        else
          qs00=0.0d0
        endif
        lcold =  tk3d.le.tf
        lwarm =  .NOT. lcold
        select case (trim(tag))
           case ('drop')
             nc3d=nc0/rho;qc3d=qc0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'drop',mx0)
             lqc3d = qc3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'drop',mx0)
           case ('crys')
             ni3d=ni0/rho;qi3d=qi0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'crys',mx0)
             lqi3d = qi3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'crys',mx0)
           case ('all')
             nc3d=nc0/rho;qc3d=qc0;ni3d=ni0/rho;qi3d=qi0
             nr3d=nr00/rho;qr3d=qr00;ns3d=ns00/rho;qs3d=qs00
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'all',mx0)
             lqc3d = qc3d.ge.qsmall
             lqi3d = qi3d.ge.qsmall
             lqs3d = qs3d.ge.qsmall
             lqr3d = qr3d.ge.qsmall
           case ('nc')
             nc3d=nc0/rho
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'drop',mx0)
             lqc3d = qc3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'drop',mx0)
           case ('qc')
             qc3d=qc0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'drop',mx0)
             lqc3d = qc3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'drop',mx0)
           case ('ni')
             ni3d=ni0/rho
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'crys',mx0)
             lqi3d = qi3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'crys',mx0)
           case ('qi')
             qi3d=qi0
             ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *         ,tk3d,pp3d,qsmall,'crys',mx0)
             lqi3d = qi3d.ge.qsmall
             ldummy=make_xxxx_timescale(dtmic,'crys',mx0)
           case ('nr'); nr3d=nr00/rho
           case ('qr'); qr3d=qr00
           case ('ns'); ns3d=ns00/rho
           case ('qs'); qs3d=qs00
          case default
            msg=sname//"Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select
c
c Make supsat timescale 
        ldummy=make_xxxx_timescale(dtmic,'all',mx0)
        la=.NOT.lcheck
      END FUNCTION set_micro_bulk2m_tomas
C
C==============================================================================
      LOGICAL FUNCTION set_termo_bulk2m(tag,tk0,qk0,pk0,w0
     *  ,v0,r0) RESULT (la)
      IMPLICIT NONE
        CHARACTER (len = *),INTENT (in)       :: tag
        real*8,dimension(:)                  :: tk0,qk0,pk0
        real*8,dimension(:),optional         :: w0,v0,r0
c Local
        integer      :: mx0
        character*18 :: sname='set_termo_bulk2m: '
        la     = .true.
        mx0=SIZE(tk0,1)
        select case (trim(tag))
           case ('all')
             tk3d=tk0;qv3d=qk0;pp3d=1.0d2*pk0;rho=pp3d/(r*tk3d)
             ww3d=w0;wvar=v0     !;rttd=r0
           case ('tk')
             tk3d=tk0;rho=pp3d/(r*tk3d)
             ldummy=make_coefs(mx0)
           case ('qv')
             qv3d=qk0
             ldummy=make_coefs(mx0)
           case ('pp')
             pp3d=1.0d2*pk0;rho=pp3d/(r*tk3d)
             ldummy=make_coefs(mx0)
           case ('ww'); ww3d=w0
           case ('wv'); wvar=v0
c          case ('rt'); rttd=r0
          case default
            msg=sname//"Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select
c 
        la=.NOT.lcheck
      END FUNCTION set_termo_bulk2m
C
C==============================================================================
      FUNCTION get_value_bulk2m(tag,dummy,field) RESULT (ra)
      IMPLICIT NONE
        CHARACTER (len = *),INTENT (in)       :: tag,dummy,field
        real*8                                :: ra(size(nc3d,1))
c Local
        character*18 :: sname='get_value_bulk2m: '
        integer      :: k
        ra=0.0d0
        k=size(ra,1)
        case_tag: select case(trim(field))
c Micro values
c       nc3d, [No/kq]  ra, [No/m^3]  ! drop concentration 
          case ('nc'); ra(1:k)=nc3d(1:k)*rho(1:k)
c       qc3d, [kq/kq]  ra, [kq/kq]   ! drop content 
          case ('qc'); ra(1:k)=qc3d(1:k)
c       ni3d, [No/kq]  ra, [No/m^3]  ! crys concentration 
          case ('ni'); ra(1:k)=ni3d(1:k)*rho(1:k)
c       qc3d, [kq/kq]  ra, [kq/kq]   ! crys content 
          case ('qi'); ra(1:k)=qi3d(1:k)
c       nc3d, [No/kq]  ra, [No/m^3]  ! rain concentration 
          case ('nr'); ra(1:k)=nr3d(1:k)*rho(1:k)
c       qr3d, [kq/kq]  ra, [kq/kq]   ! rain content 
          case ('qr'); ra(1:k)=qr3d(1:k)
c       ns3d, [No/kq]  ra, [No/m^3]  ! snow concentration 
          case ('ns'); ra(1:k)=ns3d(1:k)*rho(1:k)
c       qs3d, [kq/kq]  ra, [kq/kq]   ! snow content 
          case ('qs'); ra(1:k)=qs3d(1:k)
c       effc, [m]  ra, [micron]      ! drop effective radius 
          case ('ec'); ra(1:k)=effc(1:k)*1.0d-0
c       effi, [m]  ra, [micron]      ! crys effective radius 
          case ('ei'); ra(1:k)=effi(1:k)*1.0d-0
c       effr, [m]  ra, [micron]      ! rain effective radius 
          case ('er'); ra(1:k)=effr(1:k)*1.0d-0
c       effs, [m]  ra, [micron]      ! snow effective radius 
          case ('es'); ra(1:k)=effs(1:k)*1.0d-0
c       lamc, [m]  ra, [micron]      ! drop reverse mean diameter 
          case ('lc'); ra(1:k)=lamc(1:k)*1.0d-0
c       lami, [m]  ra, [micron]      ! crys reverse mean diameter 
          case ('li'); ra(1:k)=lami(1:k)*1.0d-0
c       lamr, [m]  ra, [micron]      ! rain reverse mean diameter 
          case ('lr'); ra(1:k)=lamr(1:k)*1.0d-0
c       lams, [m]  ra, [micron]      ! snow reverse mean diameter
          case ('ls'); ra(1:k)=lams(1:k)*1.0d-0
c
c Thermo values
c       tk3d, [K]  ra, [K]           ! temperature 
          case ('tk'); ra(1:k)=tk3d(1:k)
c       qv3d, [kq/kq]  ra, [kq/kq]   ! water vapor content
          case ('qv'); ra(1:k)=qv3d(1:k)
c       pp3d, [Pa]  ra, [mB]   ! water vapor content
          case ('pk'); ra(1:k)=pp3d(1:k)*1.0d-2
c       ww3d, [m/s]  ra, [m/c]       ! resolved vertical velocity
          case ('ww'); ra(1:k)=ww3d(1:k)
c       wvar, [m/s]  ra, [m/c]       ! sub-grid vertical velocity
          case ('wv'); ra(1:k)=wvar(1:k)
c       rttd, [K/s]  ra, [K/c]       ! radiative heating
c         case ('rt'); ra(1:k)=rttd(1:k)
c       sw, [%]  ra, [%]             ! supersaturation wrt water
          case ('sw')
            ra(1:k)=1.0d2*hugh_make_supw(pp3d(1:k)
     *,tk3d(1:k)+dt*tk3dten(1:k),qv3d(1:k)+dt*qv3dten(1:k))
c       si, [%]  ra, [%]             ! supersaturation wrt ice
          case ('si')
            ra(1:k)=1.0d2*hugh_make_supi(pp3d(1:k)
     *,tk3d(1:k)+dt*tk3dten(1:k),qv3d(1:k)+dt*qv3dten(1:k))
c
c Tendency values
c       tk3dten, [K/s]  ra, [K/c]    ! temperature tendencies
          case ('tk_tnd'); ra(1:k)=tk3dten(1:k)
c       qv3dten, [kg/kg/s]  ra, [kg/kg/c]  ! vapor content tendencies
          case ('qv_tnd'); ra(1:k)=qv3dten(1:k)
c       qc3dten, [kg/kg/s]  ra, [kg/kg/c]  ! drop cont tendencies
          case ('qc_tnd'); ra(1:k)=qc3dten(1:k)
c       qi3dten, [kg/kg/s]  ra, [kg/kg/c]  ! crys cont tendencies
          case ('qi_tnd'); ra(1:k)=qi3dten(1:k)
c       qr3dten, [kg/kg/s]  ra, [kg/kg/c]  ! rain cont tendencies
          case ('qr_tnd'); ra(1:k)=qr3dten(1:k)
c       qs3dten, [kg/kg/s]  ra, [kg/kg/c]  ! snow cont tendencies
          case ('qs_tnd'); ra(1:k)=qs3dten(1:k)
c       nc3dten, [No/kg/s]  ra, [No/m^3/c] ! drop conc tendencies
          case ('nc_tnd'); ra(1:k)=nc3dten(1:k)*rho(1:k)
c       ni3dten, [No/kg/s]  ra, [No/m^3/c] ! crys conc tendencies
          case ('ni_tnd'); ra(1:k)=ni3dten(1:k)*rho(1:k)
c       nr3dten, [No/kg/s]  ra, [No/m^3/c]  ! rain conc tendencies
          case ('nr_tnd'); ra(1:k)=nr3dten(1:k)*rho(1:k)
c       ns3dten, [No/kg/s]  ra, [No/m^3/c]  ! snow conc tendencies
          case ('ns_tnd'); ra(1:k)=ns3dten(1:k)*rho(1:k)
c
          case default
            msg=sname//"Case "//trim(field)//" is not implemented"
            call stop_model(msg,255)             
        end select case_tag
c
      END FUNCTION get_value_bulk2m
C
C==============================================================================
      FUNCTION get_rate_bulk2m(tag,field) RESULT (ra)
      IMPLICIT NONE
        CHARACTER (len = *),INTENT (in)       :: tag,field
        real*8                                :: ra(size(mpccn,1))
c Local
        character*17 :: sname='get_rate_bulk2m: '
        integer      :: k
        ra=0.0d0
        k=mx
        k=size(ra,1)
c Output: concentration tendencies in [No/m^3/s]
c Output: content tendencies in [kq/kg]
c Exception: tendencies due to Meyers are given in [No/m^3] and [kq/kq]
        case_tag: select case(trim(field))
c       nnuccc             ! change n due to con droplets freez 
          case ('nnuccc'); ra(1:k)=nnuccc(1:k)*rho(1:k)
c       mnuccc             ! change q due to con droplets freez 
          case ('mnuccc'); ra(1:k)=mnuccc(1:k)
c       nnucci             ! change n due to imm droplets freez 
          case ('nnucci'); ra(1:k)=nnucci(1:k)*rho(1:k)
c       mnucci             ! change q due to imm droplets freez 
          case ('mnucci'); ra(1:k)=mnucci(1:k)
c       nnuccd             ! change n freezing aerosol (prim ice nuc)
          case ('nnuccd'); ra(1:k)=nnuccd(1:k)*rho(1:k)
c       mnuccd             ! change q freezing aerosol (prim ice nuc) 
          case ('mnuccd'); ra(1:k)=mnuccd(1:k)
c       nnucmd             ! change n cond freezing Meyers (prim ice nuc)
          case ('nnucmd'); ra(1:k)=nnucmd(1:k)*rho(1:k)
c       mnucmd             ! change q cond freezing Meyers (prim ice nuc) 
          case ('mnucmd'); ra(1:k)=mnucmd(1:k)
c       nnucmt             ! change n cont freezing Meyers (prim ice nuc)
          case ('nnucmt'); ra(1:k)=nnucmt(1:k)*rho(1:k)
c       mnucmt             ! change q cont freezing Meyers (prim ice nuc) 
          case ('mnucmt'); ra(1:k)=mnucmt(1:k)
c       nnuccr             ! change n due to con rain freez 
          case ('nnuccr'); ra(1:k)=nnuccr(1:k)*rho(1:k)
c       mnuccr             ! change q due to con rain freez 
          case ('mnuccr'); ra(1:k)=mnuccr(1:k)
c       nnucir             ! change n due to imm rain freez 
          case ('nnucir'); ra(1:k)=nnucir(1:k)*rho(1:k)
c       mnucir             ! change q due to imm rain freez 
          case ('mnucir'); ra(1:k)=mnucir(1:k)
c       nmults             ! ice mult due to acc droplets by snow
          case ('nmults'); ra(1:k)=nmults(1:k)*rho(1:k)
c       mmults             ! change q due to ice mult droplets/snow
          case ('mmults'); ra(1:k)=mmults(1:k)
c       nmultr             ! ice mult due to acc rain by snow
          case ('nmultr'); ra(1:k)=nmultr(1:k)*rho(1:k)
c       mmultr             ! change q due to ice mult droplets/rain
          case ('mmultr'); ra(1:k)=mmultr(1:k)
c       npccn              ! change n droplets activation
          case ('npccn');  ra(1:k)=npccn(1:k)*rho(1:k)
c       pccn               ! change q droplet activation
          case ('mpccn');  ra(1:k)=mpccn(1:k)
c       nhfrc              ! change n homog freezing droplets to ice
          case ('nhfrc');  ra(1:k)=nhfrc(1:k)*rho(1:k)
c       mhfrc              ! change q homog freezing droplets to ice
          case ('mhfrc');  ra(1:k)=mhfrc(1:k)
c       nhfrr               ! change n homog freezing rain to graupel
          case ('nhfrr');  ra(1:k)=nhfrr(1:k)*rho(1:k)
c       mhfrr              ! change q homog freezing rain to graupel
          case ('mhfrr');  ra(1:k)=mhfrr(1:k)
c
c       nmltii             ! change n melting cloud ice evaporating
          case ('nmltii'); ra(1:k)=nmltii(1:k)*rho(1:k)
c       mmltii             ! change q melting cloud ice evaporating
          case ('mmltii'); ra(1:k)=mmltii(1:k)
c       nmltic             ! change n melting cloud ice to droplets
          case ('nmltic'); ra(1:k)=nmltic(1:k)*rho(1:k)
c       mmltic             ! change q melting cloud ice to droplets
          case ('mmltic'); ra(1:k)=mmltic(1:k)
c
c       nsmlts             ! change n melting snow evaporating
          case ('nsmlts'); ra(1:k)=nsmlts(1:k)*rho(1:k)
c       msmlts             ! chnage q melting snow evaporating
          case ('msmlts'); ra(1:k)=msmlts(1:k)
c       nsmltr             ! change n melting snow to rain
          case ('nsmltr'); ra(1:k)=nsmltr(1:k)*rho(1:k)
c       msmltr             ! change q melting snow to rain
          case ('msmltr'); ra(1:k)=msmltr(1:k)
c
c       nprc               ! change n autoconversion droplets
          case ('nprc');   ra(1:k)=nprc(1:k)*rho(1:k)
c       mprc               ! autoconversion droplets
          case ('mprc');   ra(1:k)=mprc(1:k)
c       ncagg              ! change n self-collection droplets
          case ('ncagg');  ra(1:k)=ncagg(1:k)*rho(1:k)
c       mcagg              ! self-collection droplets is ZERO
          case ('mcagg');  ra(1:k)=mcagg(1:k)
c       npra               ! change in n due to droplet acc by rain
          case ('npra');   ra(1:k)=npra(1:k)*rho(1:k)
c       mpra               ! accretion droplets by rain
          case ('mpra');   ra(1:k)=mpra(1:k)
c       nragg              ! change n self-collection rain
          case ('nragg');  ra(1:k)=nragg(1:k)*rho(1:k)
c       mragg              ! self-collection rain is ZERO
          case ('mragg');  ra(1:k)=mragg(1:k)
c       ncondc             ! cond/evap droplets number
          case ('ncondc'); ra(1:k)=ncondc(1:k)*rho(1:k)
c       ncondi             ! dep/sub cloud ice number
          case ('ncondi'); ra(1:k)=ncondi(1:k)*rho(1:k)
c       ncondr             ! cond/evap of rain number
          case ('ncondr'); ra(1:k)=ncondr(1:k)*rho(1:k)
c       nconds             ! dep/sub snow mass number
          case ('nconds'); ra(1:k)=nconds(1:k)*rho(1:k)
c       mcondc             ! cond/evap droplets mass
          case ('mcondc'); ra(1:k)=mcondc(1:k)
c       mcondi             ! dep/sub cloud ice mass
          case ('mcondi'); ra(1:k)=mcondi(1:k)
c       mcondr             ! cond/evap of rain mass
          case ('mcondr'); ra(1:k)=mcondr(1:k)
c       mconds             ! dep/sub snow mass
          case ('mconds'); ra(1:k)=mconds(1:k)
c
          case default
            msg=sname//"Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select case_tag
c
      END FUNCTION get_rate_bulk2m
C
C==============================================================================
      FUNCTION get_micro_bulk2m(tag,field) RESULT (ra)
      IMPLICIT NONE
        CHARACTER (len = *),INTENT (in)       :: tag,field
        real*8                                :: ra(size(mpccn,1))
c Local
        character*18 :: sname='get_micro_bulk2m: '
        integer      :: k
        ra=0.0d0
        k=mx
        k=size(ra,1)
c Output: concentration tendencies in [No/m^3/s]
c Output: content tendencies in [kq/kg]
c Exception: tendencies due to Meyers are given in [No/m^3] and [kq/kq]
c Rate values
        case_tag: select case(trim(field))
c       nnuccc             ! change n due to con droplets freez 
          case ('nnuccc'); ra(1:k)=nnuccc(1:k)*rho(1:k)
c       mnuccc             ! change q due to con droplets freez 
          case ('mnuccc'); ra(1:k)=mnuccc(1:k)
c       nnucci             ! change n due to imm droplets freez 
          case ('nnucci'); ra(1:k)=nnucci(1:k)*rho(1:k)
c       mnucci             ! change q due to imm droplets freez 
          case ('mnucci'); ra(1:k)=mnucci(1:k)
c       nnuccd             ! change n freezing aerosol (prim ice nuc)
          case ('nnuccd'); ra(1:k)=nnuccd(1:k)*rho(1:k)
c       mnuccd             ! change q freezing aerosol (prim ice nuc) 
          case ('mnuccd'); ra(1:k)=mnuccd(1:k)
c       nnucmd             ! change n cond freezing Meyers (prim ice nuc)
          case ('nnucmd'); ra(1:k)=nnucmd(1:k)*rho(1:k)
c       mnucmd             ! change q cond freezing Meyers (prim ice nuc) 
          case ('mnucmd'); ra(1:k)=mnucmd(1:k)
c       nnucmt             ! change n cont freezing Meyers (prim ice nuc)
          case ('nnucmt'); ra(1:k)=nnucmt(1:k)*rho(1:k)
c       mnucmt             ! change q cont freezing Meyers (prim ice nuc) 
          case ('mnucmt'); ra(1:k)=mnucmt(1:k)
c       nnuccr             ! change n due to con rain freez 
          case ('nnuccr'); ra(1:k)=nnuccr(1:k)*rho(1:k)
c       mnuccr             ! change q due to con rain freez 
          case ('mnuccr'); ra(1:k)=mnuccr(1:k)
c       nnucir             ! change n due to imm rain freez 
          case ('nnucir'); ra(1:k)=nnucir(1:k)*rho(1:k)
c       mnucir             ! change q due to imm rain freez 
          case ('mnucir'); ra(1:k)=mnucir(1:k)
c       nmults             ! ice mult due to acc droplets by snow
          case ('nmults'); ra(1:k)=nmults(1:k)*rho(1:k)
c       mmults             ! change q due to ice mult droplets/snow
          case ('mmults'); ra(1:k)=mmults(1:k)
c       nmultr             ! ice mult due to acc rain by snow
          case ('nmultr'); ra(1:k)=nmultr(1:k)*rho(1:k)
c       mmultr             ! change q due to ice mult droplets/rain
          case ('mmultr'); ra(1:k)=mmultr(1:k)
c       npccn              ! change n droplets activation
          case ('npccn');  ra(1:k)=npccn(1:k)*rho(1:k)
c       pccn               ! change q droplet activation
          case ('mpccn');  ra(1:k)=mpccn(1:k)
c       nhfrc              ! change n homog freezing droplets to ice
          case ('nhfrc');  ra(1:k)=nhfrc(1:k)*rho(1:k)
c       mhfrc              ! change q homog freezing droplets to ice
          case ('mhfrc');  ra(1:k)=mhfrc(1:k)
c       nhfrr               ! change n homog freezing rain to graupel
          case ('nhfrr');  ra(1:k)=nhfrr(1:k)*rho(1:k)
c       mhfrr              ! change q homog freezing rain to graupel
          case ('mhfrr');  ra(1:k)=mhfrr(1:k)
c       ncondc             ! cond/evap droplets number
          case ('ncondc'); ra(1:k)=ncondc(1:k)*rho(1:k)
c       ncondi             ! dep/sub cloud ice number
          case ('ncondi'); ra(1:k)=ncondi(1:k)*rho(1:k)
c       ncondr             ! cond/evap of rain number
          case ('ncondr'); ra(1:k)=ncondr(1:k)*rho(1:k)
c       nconds             ! dep/sub snow mass number
          case ('nconds'); ra(1:k)=nconds(1:k)*rho(1:k)
c       mcondc             ! cond/evap droplets mass
          case ('mcondc'); ra(1:k)=mcondc(1:k)
c       mcondi             ! dep/sub cloud ice mass
          case ('mcondi'); ra(1:k)=mcondi(1:k)
c       mcondr             ! cond/evap of rain mass
          case ('mcondr'); ra(1:k)=mcondr(1:k)
c       mconds             ! dep/sub snow mass
          case ('mconds'); ra(1:k)=mconds(1:k)
c
c............................................................................
c water-water interactions
c............................................................................
c       nprc               ! change n autoconversion droplets
          case ('nprc');   ra(1:k)=nprc(1:k)*rho(1:k)
c       mprc               ! autoconversion droplets
          case ('mprc');   ra(1:k)=mprc(1:k)
c       ncagg              ! change n self-collection droplets
          case ('ncagg');  ra(1:k)=ncagg(1:k)*rho(1:k)
c       mcagg              ! self-collection droplets is ZERO
          case ('mcagg');  ra(1:k)=mcagg(1:k)
c       npra               ! change in n due to droplet acc by rain
          case ('npra');   ra(1:k)=npra(1:k)*rho(1:k)
c       mpra               ! accretion droplets by rain
          case ('mpra');   ra(1:k)=mpra(1:k)
c       nragg              ! change n self-collection rain
          case ('nragg');  ra(1:k)=nragg(1:k)*rho(1:k)
c       mragg              ! self-collection rain is ZERO
          case ('mragg');  ra(1:k)=mragg(1:k)
c
c............................................................................
c crystal-crystal interactions
c............................................................................
c       nprci              ! change n autoconversion cloud ice
          case ('nprci');  ra(1:k)=nprci(1:k)*rho(1:k)
c       mprci              ! change q autoconversion cloud ice
          case ('mprci');  ra(1:k)=mprci(1:k)
c       niagg              ! change n self-collection cloud ice
          case ('niagg');  ra(1:k)=niagg(1:k)*rho(1:k)
c       miagg              ! self-collection cloud ice is ZERO
          case ('miagg');  ra(1:k)=miagg(1:k)
c
c............................................................................
c crystal-snow interactions
c............................................................................
c       nprai              ! change n accretion cloud ice by snow
          case ('nprai');  ra(1:k)=nprai(1:k)*rho(1:k)
c       mprai              ! change q accretion cloud ice by snow
          case ('mprai');  ra(1:k)=mprai(1:k)
c
c............................................................................
c snow-snow interactions
c............................................................................
c       nsagg              ! change n self-collection of snow
          case ('nsagg');  ra(1:k)=nsagg(1:k)*rho(1:k)
c       msagg,              ! self-collection of snow
          case ('msagg');  ra(1:k)=msagg(1:k)
c
c Micro values
c       nc3d, [No/kq]  ra, [No/m^3]  ! drop concentration 
          case ('nc'); ra(1:k)=nc3d(1:k)*rho(1:k)
c       qc3d, [kq/kq]  ra, [kq/kq]   ! drop content 
          case ('qc'); ra(1:k)=qc3d(1:k)
c       ni3d, [No/kq]  ra, [No/m^3]  ! crys concentration 
          case ('ni'); ra(1:k)=ni3d(1:k)*rho(1:k)
c       qc3d, [kq/kq]  ra, [kq/kq]   ! crys content 
          case ('qi'); ra(1:k)=qi3d(1:k)
c       nc3d, [No/kq]  ra, [No/m^3]  ! rain concentration 
          case ('nr'); ra(1:k)=nr3d(1:k)*rho(1:k)
c       qr3d, [kq/kq]  ra, [kq/kq]   ! rain content 
          case ('qr'); ra(1:k)=qr3d(1:k)
c       ns3d, [No/kq]  ra, [No/m^3]  ! snow concentration 
          case ('ns'); ra(1:k)=ns3d(1:k)*rho(1:k)
c       qs3d, [kq/kq]  ra, [kq/kq]   ! snow content 
          case ('qs'); ra(1:k)=qs3d(1:k)
c       effc, [m]  ra, [micron]      ! drop effective radius 
          case ('ec'); ra(1:k)=effc(1:k)*1.0d6
c       effi, [m]  ra, [micron]      ! crys effective radius 
          case ('ei'); ra(1:k)=effi(1:k)*1.0d6
c       effr, [m]  ra, [micron]      ! rain effective radius 
          case ('er'); ra(1:k)=effr(1:k)*1.0d6
c       effs, [m]  ra, [micron]      ! snow effective radius 
          case ('es'); ra(1:k)=effs(1:k)*1.0d6
c       lamc, [m]  ra, [micron]      ! drop reverse mean diameter 
          case ('lc'); ra(1:k)=lamc(1:k)*1.0d6
c       lami, [m]  ra, [micron]      ! crys reverse mean diameter 
          case ('li'); ra(1:k)=lami(1:k)*1.0d6
c       lamr, [m]  ra, [micron]      ! rain reverse mean diameter 
          case ('lr'); ra(1:k)=lamr(1:k)*1.0d6
c       lams, [m]  ra, [micron]      ! snow reverse mean diameter
          case ('ls'); ra(1:k)=lams(1:k)*1.0d6
c
c Thermo values
c       tk3d, [K]  ra, [K]           ! temperature 
          case ('tk'); ra(1:k)=tk3d(1:k)
c       qv3d, [kq/kq]  ra, [kq/kq]   ! water vapor content
          case ('qv'); ra(1:k)=qv3d(1:k)
c       pp3d, [Pa]  ra, [mB]   ! water vapor content
          case ('pk'); ra(1:k)=pp3d(1:k)*1.0d-2
c       ww3d, [m/s]  ra, [m/c]       ! resolved vertical velocity
          case ('ww'); ra(1:k)=ww3d(1:k)
c       wvar, [m/s]  ra, [m/c]       ! sub-grid vertical velocity
          case ('wv'); ra(1:k)=wvar(1:k)
c       rttd, [K/s]  ra, [m/c]       ! radiative heating
c         case ('rt'); ra(1:k)=rttd(1:k)
c
c Tendency value
c       tk3dten, [K/s]  ra, [K/c]    ! temperature tendencies
          case ('tk_tnd'); ra(1:k)=tk3dten(1:k)
c       qv3dten, [kg/kg/s]  ra, [kg/kg/c]  ! vapor content tendencies
          case ('qv_tnd'); ra(1:k)=qv3dten(1:k)
c       qc3dten, [kg/kg/s]  ra, [kg/kg/c]  ! drop cont tendencies
          case ('qc_tnd'); ra(1:k)=qc3dten(1:k)
c       qi3dten, [kg/kg/s]  ra, [kg/kg/c]  ! crys cont tendencies
          case ('qi_tnd'); ra(1:k)=qi3dten(1:k)
c       qr3dten, [kg/kg/s]  ra, [kg/kg/c]  ! rain cont tendencies
          case ('qr_tnd'); ra(1:k)=qr3dten(1:k)
c       qs3dten, [kg/kg/s]  ra, [kg/kg/c]  ! snow cont tendencies
          case ('qs_tnd'); ra(1:k)=qs3dten(1:k)
c       nc3dten, [No/kg/s]  ra, [No/m^3/c] ! drop conc tendencies
          case ('nc_tnd'); ra(1:k)=nc3dten(1:k)*rho(1:k)
c       ni3dten, [No/kg/s]  ra, [No/m^3/c] ! crys conc tendencies
          case ('ni_tnd'); ra(1:k)=ni3dten(1:k)*rho(1:k)
c       nr3dten, [No/kg/s]  ra, [No/m^3/c]  ! rain conc tendencies
          case ('nr_tnd'); ra(1:k)=nr3dten(1:k)*rho(1:k)
c       ns3dten, [No/kg/s]  ra, [No/m^3/c]  ! snow conc tendencies
          case ('ns_tnd'); ra(1:k)=ns3dten(1:k)*rho(1:k)
c
          case default
            msg=sname//"Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select case_tag
c
      END FUNCTION get_micro_bulk2m
C
C==============================================================================
      logical function execute_bulk2m_sedimentation_gcm(dt0,tag,mx0) 
     *  result(la)
      implicit none 
        real*8,  intent(in) :: dt0
        character (len = *),INTENT (in)       :: tag
        integer, intent(in) :: mx0
c Local
        real*8     :: dtmic
c counting variables
      integer      :: i,j,k
      character*25 :: sname='execute_bulk2m_sedimentation_gcm: '
        la     = .true.
        dtmic=dt0
        case_tag: select case(trim(tag))
c Make droplet distribution
          case('drop')
            la= hugh_drop_sedimentation(dtmic,mx0) 
c Make rain distribution
          case('rain')
c            la= hugh_rain_sedimentation(dtmic,mx0) 
c Make crystal distribution
          case('crys')
            la= hugh_crys_sedimentation(dtmic,mx0) 
c Make rain distribution
          case('snow')
c            la= hugh_rain_sedimentation(dtmic,mx0) 
c Make liquid phase distributions
          case('liquid')
            la= hugh_drop_sedimentation(dtmic,mx0) 
c            la= hugh_rain_sedimentation(dtmic,mx0) 
c Make solid phase distributions
          case('solid')
            la= hugh_crys_sedimentation(dtmic,mx0) 
c            la= hugh_snow_sedimentation(dtmic,mx0) 
c Make cloud water and cloud ice 
          case('cloud')
            la= hugh_drop_sedimentation(dtmic,mx0) 
            la= hugh_crys_sedimentation(dtmic,mx0) 
c Make precipitating water and  precipitating ice 
          case('precip')
c            la= hugh_rain_sedimentation(dtmic,mx0) 
c            la= hugh_snow_sedimentation(dtmic,mx0) 
c Make all distributions
          case('all')
            la= hugh_drop_sedimentation(dtmic,mx0) 
            la= hugh_crys_sedimentation(dtmic,mx0) 
c            la= hugh_rain_sedimentation(dtmic,mx0) 
c            la= hugh_snow_sedimentation(dtmic,mx0) 
          case default
            msg="Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select case_tag
c
        la=.NOT.lcheck
c
      end function execute_bulk2m_sedimentation_gcm
C
C========================================================================
C Include file: blk_add_hugh_interfaces.f
C========================================================================
      LOGICAL FUNCTION hugh_drop_condensation(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8       :: dtmic
        integer      :: mx0,k=1
        character*25 :: sname='hugh_drop_condensation: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
cigs
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'drop',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
c Make supsat timescale 
        ldummy=make_drop_timescale(dtmic,mx0)
cigs
        ldummy=hugh_make_water_condensation(dtmic,mx0)
c Check distribution
          if((maxval(nc3d) .eq. 0) .and. (maxval(nc3d) .eq. 0) .and.
     *      (maxval(mcondc).ne.0) ) then
            qc3d=qc3d+mcondc*dt0
            ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'drop',mx0)
          endif
c
        la=.NOT.lcheck
      END FUNCTION hugh_drop_condensation
C
C==============================================================================
      LOGICAL FUNCTION hugh_crys_condensation(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8       :: dtmic
        integer      :: mx0,k=1
        character*25 :: sname='hugh_crys_condensation: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
cigs
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'crys',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
c Make supsat timescale 
        ldummy=make_crys_timescale(dtmic,mx0)
cigs
        ldummy=hugh_make_ice_condensation(dtmic,mx0)
c
        la=.NOT.lcheck
c
      END FUNCTION hugh_crys_condensation
C
C==============================================================================
      LOGICAL FUNCTION hugh_snow_condensation(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8       :: dtmic
        integer      :: mx0,k=1
        character*25 :: sname='hugh_snow_condensation: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
cigs
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'snow',mx0)
c
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
c Make supsat timescale 
        ldummy=make_snow_timescale(dtmic,mx0)
        ldummy=hugh_make_ice_condensation(dtmic,mx0)
c
        la=.NOT.lcheck
      END FUNCTION hugh_snow_condensation
C
C==============================================================================
      LOGICAL FUNCTION hugh_crys_melting(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8       :: dtmic
        integer      :: mx0,k=1
        character*20 :: sname='hugh_crys_melting: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
cigs
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'crys',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
c Make supsat timescale 
        ldummy=make_crys_timescale(dtmic,mx0)
cigs
        ldummy=hugh_make_crys_melting(dtmic,mx0)
        la=.NOT.lcheck
      END FUNCTION hugh_crys_melting
C
C==============================================================================
      LOGICAL FUNCTION hugh_snow_melting(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8       :: dtmic
        integer      :: mx0,k=1
        character*20 :: sname='hugh_snow_melting: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
cigs
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'snow',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
c Make supsat timescale 
        ldummy=make_snow_timescale(dtmic,mx0)
        ldummy=hugh_make_crys_melting(dtmic,mx0)
        la=.NOT.lcheck
      END FUNCTION hugh_snow_melting
C
C==============================================================================
      LOGICAL FUNCTION hugh_rain_by_snow_accretion(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8                     :: dtmic
        integer      :: mx0,k=1
        CHARACTER*30 :: sname='hugh_rain_by_snow_accretion: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'rain',mx0)
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'snow',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall;lqs3d = qs3d.ge.1.0d-10
          lqr3d = qr3d.ge.qsmall;lqr3d = qr3d.ge.1.0d-10
        ldummy=hugh_make_rain_by_snow_accretion(dtmic,mx0)
c
        la=.NOT.lcheck
      END FUNCTION hugh_rain_by_snow_accretion
C
C==============================================================================
      LOGICAL FUNCTION hugh_drop_by_snow_accretion(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8                     :: dtmic
        integer      :: mx0,k=1
        CHARACTER*30 :: sname='hugh_drop_by_snow_accretion: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'drop',mx0)
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'snow',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=hugh_make_drop_by_snow_accretion(dtmic,mx0)
        la=.NOT.lcheck

      END FUNCTION hugh_drop_by_snow_accretion
C
C==============================================================================
      LOGICAL FUNCTION hugh_crys_by_snow_accretion(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8                     :: dtmic
        integer      :: mx0,k=1
        CHARACTER*30 :: sname='hugh_crys_by_snow_accretion: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'crys',mx0)
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'snow',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall;lqc3d = qc3d.ge.1.0d-8
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=hugh_make_crys_by_snow_accretion(dtmic,mx0)
        la=.NOT.lcheck
c
      END FUNCTION hugh_crys_by_snow_accretion
C
C==============================================================================
      LOGICAL FUNCTION hugh_crys_autoconversion(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8                     :: dtmic
        integer      :: mx0,k=1
        CHARACTER*26 :: sname='hugh_crys_autoconversion: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'crys',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=hugh_make_crys_autoconversion(dtmic,mx0)
        la=.NOT.lcheck
c
      END FUNCTION hugh_crys_autoconversion
C
C==============================================================================
      LOGICAL FUNCTION hugh_drop_by_rain_accretion(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8                     :: dtmic
        integer      :: mx0,k=1
        CHARACTER*30 :: sname='hugh_drop_by_rain_accretion: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'drop',mx0)
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'rain',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall;lqc3d = qc3d.ge.1.0d-8
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=hugh_make_drop_by_rain_accretion(dtmic,mx0)
        la=.NOT.lcheck
c
      END FUNCTION hugh_drop_by_rain_accretion
C
C==============================================================================
      LOGICAL FUNCTION hugh_drop_autoconversion(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8                     :: dtmic
        integer      :: mx0,k=1
        CHARACTER*26 :: sname='hugh_drop_autoconversion: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'drop',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall;lqc3d = qc3d.ge.1.0d-8
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=hugh_make_drop_autoconversion(dtmic,mx0)
        la=.NOT.lcheck
c
      END FUNCTION hugh_drop_autoconversion
C
C==============================================================================
      LOGICAL FUNCTION hugh_drop_freezing(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8  :: dtmic
        integer      :: mx0,k=1
        CHARACTER*20 :: sname='hugh_drop_freezing: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'drop',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=hugh_make_drop_freezing(dtmic,mx0)
        la=.NOT.lcheck
c
      END FUNCTION hugh_drop_freezing
C
C==============================================================================
      LOGICAL FUNCTION hugh_rain_freezing(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8  :: dtmic
        integer      :: mx0,k=1
        CHARACTER*20 :: sname='hugh_rain_freezing: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'rain',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=hugh_make_rain_freezing(dtmic,mx0)
        la=.NOT.lcheck
c
      END FUNCTION hugh_rain_freezing
C
C==============================================================================
      LOGICAL FUNCTION hugh_crys_nucleation(dt0,k0)  RESULT (la)     
      IMPLICIT NONE
        REAL*8,             INTENT (in)       :: dt0
        INTEGER,            INTENT (in)       :: k0
c my Local
        real*8       :: dtmic
        integer      :: mx0,k=1
        logical      :: lnucl
        character*22 :: sname='hugh_crys_nucleation: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'crys',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=hugh_make_crys_nucleation(dtmic,mx0)
        la=.NOT.lcheck
c
      END FUNCTION hugh_crys_nucleation
C
C==============================================================================
      LOGICAL FUNCTION hugh_drop_nucleation(dt0,k0)  RESULT (la)     
      IMPLICIT NONE
        REAL*8,             INTENT (in)       :: dt0
        INTEGER,            INTENT (in)       :: k0
c Local
        real*8       :: dtmic
        integer      :: mx0,k=1
        LOGICAL      :: lnucl
        character*22 :: sname='hugh_drop_nucleation: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'drop',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=hugh_make_drop_nucleation(dtmic,mx0)
        la=.NOT.lcheck
c
      END FUNCTION hugh_drop_nucleation
C
C========================================================================
C Include file: blk_add_ice_nucleation.f
C==============================================================================
      logical function hugh_make_drop_freezing(dtmic,mx0) result(la)
      implicit none
        real*8,             intent (in)       :: dtmic
        integer,            intent(in)        :: mx0
c Local
      real*8  :: Nacnt
c counting variables
      integer :: k
      logical :: lfreezC,lfreezI,lfreezH
      character*25 :: sname='hugh_make_drop_freezing: '
        la     = .true.
        lfreezC = .true.!; lfreezC = .false.
        lfreezI = .true.!; lfreezI = .false.
        lfreezH = .true.!; lfreezH = .false.
c
c freezing of cloud droplets only allowed below -4 C
        do k=1,mx0
          nnuccc(k)=0.0d0; mnuccc(k)=0.0d0
          nnucci(k)=0.0d0; mnucci(k)=0.0d0
          nhfrc(k) =0.0d0; mhfrc(k) =0.0d0
          Nacnt = 0.0d0
	  if_drop_exists: if (lqc3d(k) .and. tk3d(k).lt.269.15) then

c number of contact nuclei (m^-3) from Meyers et al., 1992
c factor of 1000 is to convert from L^-1 to m^-3
c meyers curve
c     Nacnt = exp(-2.80d0+0.262d0*(tf-tk3d(k)))*1.0d3
c cooper curve
         Nacnt =  5.*exp(0.304*(tf-tk3d(k)))
         if(Nacnt.gt.209000.) Nacnt=209000.
c flecther
c	nacnt = 0.01*exp(0.6*(tf-tk3d(k)))
c drop contact freezing
           drop_freez_contact: if(lfreezC) then
	     mnuccc(k) = pi*pi/3.*Dap(k)*Nacnt*rhow*cdist1(k)*
     1              gamma(pgam(k)+5.0d0)/lamc(k)**4
	     nnuccc(k) = 2.0d0*pi*Dap(k)*Nacnt*cdist1(k)*
     1               gamma(pgam(k)+2.0d0)/lamc(k)
c        write(6,*)"2M_QAUT1",nnuccc(k),k,Nacnt,lamc(k)  
           endif drop_freez_contact
c drop immersion freezing (bigg 1953)
           drop_freez_immersion: if(lfreezI) then
             mnucci(k) = pi*pi/36.*rhow*
     1             cdist1(k)*gamma(7.0d0+pgam(k))*
     1              bimm*exp(aimm*(tf-tk3d(k)))/
     1              lamc(k)**3/lamc(k)**3
             nnucci(k) = pi/6.*
     1             cdist1(k)*gamma(pgam(k)+4.0d0)*bimm*
     1            exp(aimm*(tf-tk3d(k)))/lamc(k)**3
c        write(6,*)"2M_QAUT2",nnucci(k),k,lamc(k),tk3d(k)  
           endif drop_freez_immersion
c drop homogeneous freezing
           drop_freez_homogeneous: if(lfreezH) then
	      mhfrc(k)  = 0.0d0
              nhfrc(k)  = 0.0d0
           endif drop_freez_homogeneous
c put in a catch here to prevent divergence between number conc. and
c mixing ratio
	     if (mnuccc(k).eq.0.) nnuccc(k) = 0.0d0
             if (nnuccc(k).eq.0.) mnuccc(k) = 0.0d0
	     if (mnucci(k).eq.0.) nnucci(k) = 0.0d0
             if (nnucci(k).eq.0.) mnucci(k) = 0.0d0
	     if (mhfrc(k) .eq.0.) nhfrc(k)  = 0.0d0
	     if (nhfrc(k) .eq.0.) mhfrc(k)  = 0.0d0
          else
	  endif if_drop_exists
        enddo
        la=.NOT.lcheck
c
      end function hugh_make_drop_freezing
C
C==============================================================================
      logical function hugh_make_rain_freezing(dtmic,mx0) result(la)
      implicit none
        real*8,             intent (in)       :: dtmic
        integer,            intent(in)        :: mx0
c Local
      real*8  :: Nacnt
c counting variables
      integer :: k
      logical :: lfreezC,lfreezI,lfreezH
      character*25 :: sname='hugh_make_rain_freezing: '
        la     = .true.
        lfreezC = .true.!; lfreezC = .false.
        lfreezI = .true.!; lfreezI = .false.
        lfreezH = .true.!; lfreezH = .false.
c
c freezing of rain drops
c number of contact aerosol given by meyers et al 1992
c Freezing time of drop assumed to be instantaneous once freezing occurs
c freezing allowed below -4 C
        do k=1,mx0
          nnuccr(k)=0.0d0; mnuccr(k)=0.0d0
          nnucir(k)=0.0d0; mnucir(k)=0.0d0
          nhfrr(k) =0.0d0; mhfrr(k) =0.0d0
	   if_rain_exists: if (tk3d(k) .lt. 269.15 .and. lqr3d(k)) then
c meyers curve
	      Nacnt = exp(-2.80d0+0.262d0*(tf-tk3d(k)))*1.0d3
c cooper curve
c        Nacnt =  5.*exp(0.304*(tf-tk3d(k)))
c flecther curve
c	nacnt = 0.01*exp(0.6*(tf-tk3d(k)))
c rain contact freezing
           rain_freez_contact: if(lfreezC) then
	      mnuccr(k) = 8.0d0*pi*pi*Dap(k)*Nacnt*rhow*
     1              n0rr(k)/lamr(k)**5	
	      nnuccr(k)=2.0d0*pi*Dap(k)*Nacnt*n0rr(k)/lamr(k)**2
           endif rain_freez_contact
c rain immersion freezing
           rain_freez_immersion: if(lfreezI) then
	      mnucir(k) = 20.0d0*pi*pi*rhow*nr3d(k)*bimm*
     1             exp(aimm*(tf-tk3d(k)))/lamr(k)**3/lamr(k)**3

	      nnucir(k) = pi*nr3d(k)*bimm*
     1              exp(aimm*(tf-tk3d(k)))/lamr(k)**3
           endif rain_freez_immersion
c rain homogeneous freezing
           rain_freez_homogeneous: if(lfreezH) then
	      mhfrr(k)  = 0.0d0
              nhfrr(k)  = 0.0d0
           endif rain_freez_homogeneous
c prevent divergence between mixing ratio and number conc
	    if (mnuccr(k).eq.0.) nnuccr(k) = 0.0d0
            if (nnuccr(k).eq.0.) mnuccr(k) = 0.0d0
	    if (mnucir(k).eq.0.) nnucir(k) = 0.0d0
            if (nnucir(k).eq.0.) mnucir(k) = 0.0d0
	    if (mhfrr(k) .eq.0.) nhfrr(k)  = 0.0d0
	    if (nhfrr(k) .eq.0.) mhfrr(k)  = 0.0d0
           else
	   end if if_rain_exists
        enddo
        la=.NOT.lcheck
c
      end function hugh_make_rain_freezing
C
C==============================================================================
      logical function hugh_make_crys_nucleation(dtmic,mx0) result(la)
      implicit none
      real*8, parameter  :: temp1= -5.0d0
      real*8, parameter  :: temp2= -2.0d0
      real*8, parameter  :: temp3=-20.0d0
      real*8, parameter  :: as_meyers=-0.639d0   ! 
      real*8, parameter  :: bs_meyers=0.1296d0   ! [1/K]
      real*8, parameter  :: cs_meyers=1.0d3      ! [No/m^3]
      real*8, parameter  :: at_meyers=-2.8d0     ! 
      real*8, parameter  :: bt_meyers=0.2620d0   ! [1/K]
      real*8, parameter  :: ct_meyers=1.0d3      ! [No/m^3]
        real*8,             intent (in)       :: dtmic
        integer,            intent(in)        :: mx0
      real*8  :: dum,dumt,dumqv,dumqsi,dumqss,kc2
      real*8  :: wef,kc0,kc1,supw,supi,tpc,Nacnt,Nacns
      real*8, dimension(size(nc3d,1))  :: sw,si
      logical,dimension(size(nc3d,1))  :: lcnucl
c
c      real*8,dimension(mx) :: wnuc,anuc,bnuc
c     
c counting variables
      integer :: k
c
      logical :: lmeyers,lmeyersS,lmeyersT,lhugh
      character*27 :: sname='hugh_make_crys_nucleation: '
        la     = .true.
        lhugh    = .true.!; lhugh    = .false.
        lmeyers  = .true.!; lmeyers  = .false.
        lmeyersS = .true.!; lmeyersS = .false.
        lmeyersT = .true.!; lmeyersT = .false.
c
c nucleation of cloud ice from homogeneous and heterogeneous freezing on 
c aerosol calculated from curve fit to parcel simulations using rates from
c khvorostyanov and curry, 2005
        do k=1,mx0
           nnuccd(k) = 0.0d0; mnuccd(k) = 0.0d0
           nnucmd(k) = 0.0d0; mnucmd(k) = 0.0d0
           nnucmt(k) = 0.0d0; mnucmt(k) = 0.0d0
           wnuc(k) = 0.0d0; anuc(k) = 0.0d0; bnuc(k) = 0.0d0
           wef = 0.0d0; kc0=0.0d0; kc1=0.0d0; kc2=0.0d0
           Nacnt=0.0d0;Nacns=0.0d0;dum = 0.0d0
           lcnucl(k)=.false.
	   dumt = tk3d(k)
	   dumqv = qv3d(k)
c updated ice, water saturation mixing ratio
	   dumqsi = 0.622*polysvp(dumt,1)/
     1          (pp3d(k)-polysvp(dumt,1))
	   dumqss = 0.622*polysvp(dumt,0)/
     1          (pp3d(k)-polysvp(dumt,0))
c updated saturation ratio with respect to water and ice
	   supw = dumqv/dumqss
	   supi = dumqv/dumqsi
           sw(k) = 1.0d2*(supw-1.0d0)
           si(k) = 1.0d2*(supi -1.0d0)
           tpc= dumt - tf
c Meyers formulation 
         if_meyers: if(lmeyers) then
           ice_saturation: if(si(k) .gt. 0 ) then
             meyers_deposition: if(tpc .lt. temp1 .AND. lmeyersS) then
               si(k) = min(25.0d0, si(k))
c meyers supersaturation curve, [No/kq]
c	       Nacns = exp(-0.639d0+0.262d0*supi)*1.0d3
	       Nacns = cs_meyers*exp(as_meyers+bs_meyers*si(k))/rho(k)
               nnucmd(k) = max(Nacns-ni3d(k),0.0d0)
               mnucmd(k) = nnucmd(k)*mi0
             endif meyers_deposition
             meyers_contact: if(tpc .lt. temp2 .AND. lmeyersT) then
               tpc = max(temp3, tpc)
c meyers temperature curve, [no/kq]
c	     Nacnt = exp(-2.80d0+0.262d0*(tf-tk3d(k)))*1.0d3
	     Nacnt = ct_meyers*exp(at_meyers-bt_meyers*tpc)/rho(k)
               nnucmt(k) = max(Nacnt-ni3d(k),0.0d0)
               mnucmt(k) = nnucmt(k)*mi0
             endif meyers_contact
           endif ice_saturation
           kc0 = nnucmd(k) + nnucmt(k)
           kc0 = Nacns + Nacnt
           kc1 = max(kc0-ni3d(k),0.0d0)
           if(kc1 .ne. 0) then
             nnucmd(k) = NacnS*kc1/kc0
             mnucmd(k) = nnucmd(k)*mi0
             nnucmt(k) = NacnT*kc1/kc0
             mnucmt(k) = nnucmt(k)*mi0
           endif
         endif if_meyers
c freezing of aerosol only allowed below -5 C
c and above deliquescence threshold of 80%
c and above ice saturation 
         if_hugh: if(lhugh) then
	   if (dumt.ge.251. .or.
c	   if (dumt.ge.271. .or.
     1        dumqv.le.dumqsi .or. dumqv/dumqss.le.0.75) go to 10
	     lcnucl(k)=.true.
c           print *,dumt,dumqv,dumqsi,dumqv/dumqss
c determine if greater than threshold RH
c qv = 0.25
c	if (dumt.lt.237.) then
c	dum = 0.34987+0.0022953*dumt
c	else if (dumt.ge.237.) then
c	dum = -0.49692+0.0058683*dumt
c	end if
c qv = 0.75

           if (dumt.lt.210.) then
             dum = 0.63335+0.00088638*dumt
           else if (dumt.ge.210.) then  
             dum = -0.012998+0.0039642*dumt
           end if
c	   if_supsat: if (sat1.ge.dum) then
           if(supw .lt. dum) cycle
c nucleation occurs
c determine effective vertical velocity
c             wef = -cp/g*rttd(k)+ww3d(k)
              wef = ww3d(k)
c set minimum value of 10 cm/s for sub-grid eff. vertical vel.
              wef = max(wef,0.1d0)
              wnuc(k) = wef
c qv = 0.75
            dumt=251.
            if (dumt.lt.225.) then
            anuc(k) = -1.689502e10+2.551999e8*dumt
     1        -1.26874e6*dumt**2+2082.19*dumt**3
            bnuc(k) = -132.729+1.977388*dumt
     1        -9.66254e-3*dumt**2+1.56567e-5*dumt**3
            else if (dumt.ge.225.and.dumt.lt.230.) then
              anuc(k) = 2.90801e8-1.244498e6*dumt
              bnuc(k) = -0.986563+0.01043599*dumt
            else if (dumt.ge.230.) then
              anuc(k) = 3.43141e9-4.13828e7*dumt
     1          +1.66606e5*dumt**2-223.907*dumt**3
              bnuc(k) = -4257.98+53.0047*dumt
     1          -0.21965*dumt**2+3.03093e-4*dumt**3
            end if
c below for qv = 0.25
c        if (dumt.lt.223.) then
c        anuc(k) = -1.6842018e10+2.677928e8*dumt
c     1  -1.389664e6*dumt**2+2366.984*dumt**3
c        bnuc(k) = -174.99+2.6159123*dumt
c     1  -1.2881e-2*dumt**2+2.106594e-5*dumt**3
c        else if (dumt.ge.223.and.dumt.lt.227.) then
c        anuc(k) = 4.97688e8-2.164938e6*dumt
c        bnuc(k) = 1.44143-1.75774e-4*dumt
c        else if (dumt.ge.227.) then
c        anuc(k) = 2.812037e9-3.41905e7*dumt
c     1   +1.38705e5*dumt**2-187.728*dumt**3
c        bnuc(k) = 175.1624-2.1298*dumt
c     1   +8.74645e-3*dumt**2-1.205017e-5*dumt**3
c        end if
	    kc2 = anuc(k)*wef**bnuc(k)
c            kc2 = abs(kc2)
            kc2 = max(kc2,0.d0)
            if (kc2.gt.ni3d(k)) then
              nnuccd(k) = (kc2-ni3d(k)) !/ dtmic
              mnuccd(k) = nnuccd(k)*mi0
            end if
c          endif if_supsat
 10       continue	
         endif if_hugh
        enddo
c Storage: concentration in [No/kq]
           nnucmd = nnucmd !* rho
           nnucmt = nnucmt ! *rho
c Storage: content in [kq/kq]
           mnucmd = mnucmd ! *rho
           mnucmt = mnucmt ! *rho
        la=.NOT.lcheck
      end function hugh_make_crys_nucleation
C
C==============================================================================
      logical function hugh_make_snow_freezing(dtmic,mx0) result(la)
      implicit none
        real*8,             intent (in)       :: dtmic
        integer,            intent(in)        :: mx0
c Local
      real*8  :: dqsdt,dqsidt
c counting variables
      integer :: k
      character*25 :: sname='hugh_make_rain_freezing: '
        la     = .true.
        do k=1,mx0
c
        enddo
        la=.NOT.lcheck
      end function hugh_make_snow_freezing
C
C==============================================================================
C Include file: blk_add_drop_nucleation.f
C==============================================================================
      LOGICAL FUNCTION hugh_make_drop_nucleation(dt0,mx0) RESULT (la)     
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: mx0
c droplet activation/freezing aerosol
       real*8  
     1       ct,      ! droplet activation parameter
     1       temp1,   ! dummy temperature
     1       sat1,    ! dummy saturation
     1       sigvl,   ! surface tension liq/vapor
     1       kel,     ! Kelvin parameter
     1       kc2      ! total ice nucleation rate
       real*8 
     1       cry,kry   ! aerosol activation parameters
c
c working parameters for aerosol activation
        real*8 
     1   aact,gamm,gg,psi,eta1,eta2,sm1,sm2,smax,uu1,uu2,
     1   alpha
c dummy variables
	real*8 
     1       dum,dum1,dum2,dum3,dum4,dumt,qdumt,dumqv,
     1       dumqss,dumqsi,dums,dume
c my Local
        real*8       :: p0,w0,v0,r0,ra,dtmic
        integer :: k=1
        LOGICAL :: lnucl                            
        CHARACTER*27 :: sname='hugh_make_drop_nucleation: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
c        ldummy=make_coefs(mx0)

c update variables
       altitude_loop: do k=1,mx0
c      do k=1,mx0
        npccn(k)=0.0d0; mpccn(k)=0.0d0
	dumt = tk3d(k)
	dumqv = qv3d(k)
	dumqss = 0.622*polysvp(dumt,0)/
     1          (pp3d(k)-polysvp(dumt,0))
        dume  = pp3d(k)*dumqv / (0.622d0 + dumqv)
	lnucl=.false.
	if (dumqv/dumqss.ge.0.95) then
	lnucl=.true.
c droplet activation from Abdul-Razzak and Ghan (2000)
c convert heating rate to effective vertical velocity
c effective vertical velocity (m/s)
c          dum = ww3d(k)-cp/g*rttd(k)
           dum = ww3d(k)  !-cp/g*rttd(k)
c add sub-grid vertical velocity
	   dum = dum+wvar(k)
c assume minimum eff. sub-grid velocity 0.20 m/s
c          dum = max(dum,0.40d0)                 
           dum = max(dum,0.20d0)                 
           sigvl = 0.0761-1.55e-4*(tk3d(k)-tf)
           aact = 2.*mw/(rhow*rr)*sigvl/tk3d(k)
           alpha = g*mw*xxlv(k)/(cp*rr*tk3d(k)**2)-
     1       g*ma/(rr*tk3d(k))
           gamm = rr*tk3d(k)/(evs(k)*mw)+mw*xxlv(k)**2/
     1      (cp*pp3d(k)*ma*tk3d(k))
           gg = 1./(rhow*rr*tk3d(k)/(evs(k)*Dv(k)*mw)+
     1       xxlv(k)*rhow/(kap(k)*tk3d(k))*(xxlv(k)*mw/
     1         (tk3d(k)*rr)-1.))
           psi = 2./3.*(alpha*dum/gg)**0.5*aact
           eta1 = (alpha*dum/gg)**1.5/
     1          (2.*pi*rhow*gamm*nanew1)
           eta2 = (alpha*dum/gg)**1.5/
     1          (2.*pi*rhow*gamm*nanew2)
           sm1 = 2./bact**0.5*(aact/(3.*rm1))**1.5
           sm2 = 2./bact**0.5*(aact/(3.*rm2))**1.5
           dum1 = 1./sm1**2*(f11*(psi/eta1)**1.5+
     1            f21*(sm1**2/(eta1+3.*psi))**0.75)
           dum2 = 1./sm2**2*(f12*(psi/eta2)**1.5+
     1            f22*(sm2**2/(eta2+3.*psi))**0.75)
           smax = 1./(dum1+dum2)**0.5

           uu1 = 2.*log(sm1/smax)/(4.242*log(sig1))
           uu2 = 2.*log(sm2/smax)/(4.242*log(sig2))
           dum1 = nanew1/2.*(1.-derf1(uu1))
           dum2 = nanew2/2.*(1.-derf1(uu2))  
           dum3 = (dum1+dum2)/rho(k)  !convert to kg-1
c make sure this value is not greater than total number of aerosol
	      dum3 = min((nanew1+nanew2)/rho(k),dum3)
	      dum4 = (dum3-nc3d(k)) ! / dtmic
	      dum4 = max(0.d0,dum4)
              dum4=MERGE(dum4,0.0d0, (dum4 .gt. 0)) /dtmic
c              dum4=MERGE(dum3-nc3d(k),0.0d0,(dum3-nc3d(k)) .gt. 0)
	      npccn(k) = dum4/2.
c      write(6,*)"2M_BLK",k,sigvl,aact,alpha,gamm,gg,psi,eta1,eta2,sm1
c    &,sm2,smax,uu1,uu2,tk3d,qv3d,pp3d,dume,dum,dumt,dumqv,dum1,dum2,
c    &dum3,dum4,nc3d(k),npccn(k),qc3d,mpccn,xxlv(k),xxls(k),
c    &xlf(k),rho(k)
c      if(nc3d(k).ne.0..and.npccn(k).ne.0.) write(6,*)"2M_BLK",
c    &nc3d(k),npccn(k),dume,dumt,pp3d(k)
Cigs
Cigs all new droplets have the same mass: 
              mpccn(k) = mw0*npccn(k)
        end if  ! qv/qvs

      enddo altitude_loop
c
      la=.NOT.lcheck
C
      END FUNCTION hugh_make_drop_nucleation
C
C==============================================================================
C
      LOGICAL FUNCTION hugh_drop_nucleation_init() RESULT (la)
      IMPLICIT NONE
        CHARACTER*27 :: sname='hugh_drop_nucleation_init: '
        la     = .true.
c aerosol activation parameters
c parameters currently set for ammonium sulfate

           mw = 0.018
           osm = 1.
           vi = 3.
           epsm = 1.0
           rhoa = 1777.
           mapp = 0.132
           ma = 0.0284
           rr = 8.3187
           bact = vi*osm*epsm*mw*rhoa/(mapp*rhow)
c
	   rhow = 997.

c aerosol size distribution parameters currently set for mpace
c mode 1

           rm1 = 0.052e-6
           sig1 = 2.04
c          nanew1 = 72.2e6
           nanew1 = 172.2e6
           f11 = 0.5*exp(2.5*(log(sig1))**2)
           f21 = 1.+0.25*log(sig1)

c mode 2
           rm2 = 1.3e-6
           sig2 = 2.5
c          nanew2 = 1.8e6
           nanew2 = 11.8e6
           f12 = 0.5*exp(2.5*(log(sig2))**2)
           f22 = 1.+0.25*log(sig2)
C
        la=.NOT.lcheck
      END FUNCTION hugh_drop_nucleation_init
C========================================================================
C Include file: blk_add_water_water_interactions.f
C========================================================================
      LOGICAL FUNCTION hugh_make_drop_autoconversion(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: ldrop_auto,ldrop_self
        CHARACTER*31 :: sname='hugh_make_drop_autoconversion: '
        la     = .true.
        ldrop_auto = .true.!; ldrop_auto = .false.
        ldrop_self = .true.; ldrop_self = .false.

c Process description
        altitude_loop: do k=1,mx0
          nprc(k)=0.0d0; mprc(k)=0.0d0
          ncagg(k)=0.0d0; mcagg(k)=0.0d0
c
c use minimum value of 1.e-8 to prevent floating point error
	  if_drop_exists: if (lqc3d(k)) then
c.................................................................
c autoconversion of cloud liquid water to rain
c formula from Beheng (1994)
c using numerical simulation of stochastic collection equation
c and initial cloud droplet size distribution specified
c as a gamma distribution THAT depends on MASSES but not on RADII !!!
	   if_drop_auto: if (ldrop_auto) then
	      nprc(k) = 7.7d9
     1                *(6.0d28/rho(k)*pgam(k)**(-1.7)*
     1           (1.0d-6*rho(k)*nc3d(k))**(-3.3)*
     1           (1.0d-3*rho(k)*qc3d(k))**( 4.7))
	      mprc(k) = 6.0d28/rho(k)*pgam(k)**(-1.7)*
     1           (1.0d-6*rho(k)*nc3d(k))**(-3.3)*     ! conc in [No/cc]
     1           (1.0d-3*rho(k)*qc3d(k))**( 4.7)      ! cont in [g/cc]
c        write(6,*)"2M_QAUT",nprc(k),k,nc3d(k),qc3d(k)
           else
	   end if if_drop_auto
c
c.......................................................................
c self-collection of cloud droplets
c from Beheng(1994)
c from numerical simulation of the stochastic collection equation
c as described above for autoconversion
	   if_drop_self: if (ldrop_self) then
	     ncagg(k) = -5.5d16/rho(k)*pgam(k)**(-0.63)*
     1             (1.0d-3*rho(k)*qc3d(k))**2
             mcagg(k)=0.0d0
          else
	   end if if_drop_self
	  endif if_drop_exists
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_make_drop_autoconversion
C
C==============================================================================
      LOGICAL FUNCTION hugh_make_drop_by_rain_accretion(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: ldrop_rain,lrain_self
        CHARACTER*34 :: sname='hugh_make_drop_by_rain_accretion: '
        la     = .true.
        ldrop_rain = .true.!; ldrop_rain = .false.
        lrain_self = .true.; lrain_self = .false.
c Process description
        altitude_loop: do k=1,mx0
          npra(k)=0.0d0; mpra(k)=0.0d0
          nragg(k)=0.0d0; mragg(k)=0.0d0
c use minimum value of 1.d-8 to prevent floating point error
c
	  if_drop_rain_exists: if (lqc3d(k) .AND. lqr3d(k)) then
c.......................................................................
c accretion of cloud liquid water by rain
c continuous collection equation with
c gravitational collection kernel, droplet fall speed neglected
c use formula from Beheng
c
	    if_drop_rain: if (ldrop_rain) then
              mpra(k) = 6.0d0*rho(k)*(qc3d(k)*qr3d(k))
              npra(k) = mpra(k)*rho(k)/1.0d3*(nc3d(k)*rho(k)/1.0d6)/
     1        (qc3d(k)*rho(k)/1.0d3)*1.0d6/rho(k)

c continuous collection
c	mpra(k) = pi/4.*rho(k)*qc3d(k)*Ecr*arn(k)*n0rr(k)*
c     1           gamma(br+3.0d0)/lamr(k)**(br+3.0d0)
c	npra(k) = pi/4.*rho(k)*nc3d(k)*Ecr*arn(k)*n0rr(k)*
c     1           gamma(br+3.0d0)/lamr(k)**(br+3.0d0)
            else
	    endif if_drop_rain
	  endif if_drop_rain_exists
	  if_rain_exists: if (lqr3d(k)) then
c
c.......................................................................
c Self-collection of rain drops
c from Beheng(1994)
c from numerical simulation of the stochastic collection equation
c as descrined above for autoconversion
	    if_rain_self: if (lrain_self) then
	      nragg(k) = -8.0d0*nr3d(k)*qr3d(k)*rho(k)
	      mragg(k) =  0.0d0
            else
	    endif if_rain_self
	  endif if_rain_exists
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_make_drop_by_rain_accretion
C
C
C========================================================================
C Include file: blk_add_ice_ice_interactions.f
C==============================================================================
      LOGICAL FUNCTION hugh_make_crys_by_snow_accretion(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: lcrys_snow,lsnow_self
        CHARACTER*29 :: sname='hugh_make_crys_by_snow_accretion: '
        la     = .true.
        lcrys_snow = .true.; lcrys_snow = .false.
        lsnow_self = .true.; lsnow_self = .false.

c Process description
        altitude_loop: do k=1,mx0
          nprai(k)=0.0d0; mprai(k)=0.0d0
          nsagg(k)=0.0d0; msagg(k)=0.0d0
c
	  if_crys_snow_exists: 
     *      if (lqi3d(k) .AND. lqs3d(k) .and. lcold(k)) then
c.......................................................................
c Accretion of cloud ice by snow
c For this calculation, it is assumed that the Vs >> Vi
c and Ds >> Di for continuous collection
	    if_crys_snow: if (lcrys_snow) then
	      mprai(k) = pi/4.0d0*asn(k)*qi3d(k)*rho(k)*
     1                n0s(k)*Eii(k)*gamma(bs+3.0d0)/
     1                lams(k)**(bs+3.0d0)	
	      Nprai(k) = pi/4.0d0*asn(k)*ni3d(k)*rho(k)*
     1                n0s(k)*Eii(k)*gamma(bs+3.0d0)/
     1                lams(k)**(bs+3.0d0)
            else
	    endif if_crys_snow
	  endif if_crys_snow_exists
	  if_snow_exists: if (lqs3d(k) .and. lcold(k)) then
c
c.......................................................................
c snow aggregation from passarelli, 1978, used by reisner, 1998
c this is hard-wired for bs = 0.4 for now
	    if_snow_self: if (lsnow_self) then
	       Nsagg(k) = -1108.*asn(k)*Eii(k)*
     1           pi**((1.0d0-bs)/3.)*rhosn**((-2.0d0-bs)/3.)*rho(k)**
     1           ((2.0d0+bs)/3.)*qs3d(k)**((2.0d0+bs)/3.)*
     1           (ns3d(k)*rho(k))**((4.0d0-bs)/3.)/
     1           (4.0d0*720.0d0*rho(k))
               msagg(k) = 0.0d0
            else
	    endif if_snow_self
	  endif if_snow_exists
        enddo altitude_loop
        la=.NOT.lcheck
      END FUNCTION hugh_make_crys_by_snow_accretion
C
C
C==============================================================================
      LOGICAL FUNCTION hugh_make_crys_autoconversion(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      real*8  :: coffi
      logical :: lcrys_auto,lcrys_self
        CHARACTER*31 :: sname='hugh_make_crys_autoconversion: '
        la     = .true.
        lcrys_auto = .true.; lcrys_auto = .false.
        lcrys_self = .true.; lcrys_self = .false.
c Process description
        altitude_loop: do k=1,mx0
          nprci(k)=0.0d0; mprci(k)=0.0d0
          niagg(k)=0.0d0; miagg(k)=0.0d0
c use minimum value of 1.e-8 to prevent floating point error
	  if_crys_exists: if (lqi3d(k)) then
c.......................................................................
c Autoconversion of cloud ice to snow
c following Harrington et al. (1995) with modification
c here it is assumed that autoconversion can only occur when the
c ice is growing, i.e. in conditions of ice supersaturation
c only allow autoconversion when mean cloud ice diameter > 0.5*dcs
	   if_crys_auto: if (lcrys_auto) then
	     if_cold_supi: if(lcold(k) .and. qvqvsi(k).ge.1.) then
               coffi = 2./lami(k)
               if_size: if (coffi.ge.dcs) then
	         nprci(k) = 2.*pi*(qv3d(k)-qvi(k))*rho(k)/
     1             (dcs*rhoi)*n0i(k)*exp(-lami(k)*dcs)*Dv(k)/abi(k)
	         mprci(k) = pi*rhoi*dcs**3/6.*nprci(k)
               else
               endif if_size
             endif if_cold_supi
           else
	   endif if_crys_auto
c.......................................................................
c self-collection of ice crystals
	   if_crys_self: if (lcrys_self) then
	     niagg(k) = 0.0d0
             miagg(k) = 0.0d0
          else
	   end if if_crys_self
	  endif if_crys_exists
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_make_crys_autoconversion
C
C========================================================================
C Include file: blk_add_water_ice_interactions.f
C========================================================================
      LOGICAL FUNCTION hugh_make_rain_by_snow_accretion(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
c      real*8  :: dc0,ds0,eci,dum,uns,dum1,fmult
      real*8  :: uns,ums,unr,umr,dum1,fmult
      logical :: lrain_snow,lsnow_mult
        CHARACTER*29 :: sname='hugh_make_rain_by_snow_accretion: '
        la     = .true.
        lrain_snow = .true.!; lrain_snow = .false.
        lsnow_mult = .true.; lsnow_mult = .false.

c Process description
        altitude_loop: do k=1,mx0
          npracs(k)=0.0d0; mpracs(k)=0.0d0
          nmultr(k)=0.0d0; nmultr(k)=0.0d0
c
	  if_rain_snow_exists: 
     *      if (lqr3d(k) .AND. lqs3d(k) .and. lcold(k)) then
c.......................................................................
c accretion of cloud droplets onto snow/graupel
c here use continuous collection equation with
c simple gravitational collection kernel ignoring
c collisions between droplets/cloud ice
c since minimum size ice particle for accretion is 50 - 250 micron
c
	    if_rain_snow: if (lrain_snow) then
	      ums = asn(k)*gamma(4.0d0+bs)/(6.0d0*lams(k)**bs)
	      umr = arn(k)*gamma(4.0d0+br)/(6.0d0*lamr(k)**br)
	      uns = asn(k)*gamma(1.0d0+bs)/lams(k)**bs
	      unr = arn(k)*gamma(1.0d0+br)/lamr(k)**br 

c set reaslistic limits on fallspeeds
	      ums=min(ums,3.d0)
	      uns=min(uns,3.d0)
	      umr=min(umr,9.1d0)
	      unr=min(unr,9.1d0)

	      mpracs(k) = pi*pi*ecr*(((1.2*umr-0.95*ums)**2+
     1             0.08*ums*umr)**0.5*rhosn*rho(k)*
     1            n0rr(k)*n0s(k)*
     1             (5.0d0/(lams(k)**6*lamr(k))+
     1             2./(lams(k)**5*lamr(k)**2)+
     1             0.5/(lams(k)**4*lamr(k)**3)))+
     1             pi*pi*ecr*(((1.2*umr-0.95*ums)**2+
     1             0.08*ums*umr)**0.5*rhow*rho(k)*
     1             n0rr(k)*n0s(k)*
     1             (5./(lamr(k)**6*lams(k))+
     1             2./(lamr(k)**5*lams(k)**2)+
     1             0.5/(lamr(k)**4*lams(k)**3)))
	      npracs(k) = pi/2.*rho(k)*ecr*(1.7*(unr-uns)**2+
     1           0.3*unr*uns)**0.5*n0rr(k)*n0s(k)*
     1           (1./(lamr(k)**3*lams(k))+
     1            1./(lamr(k)**2*lams(k)**2)+
     1            1./(lamr(k)*lams(k)**3))
c make sure pracs does not exceed total rain mixing ratio
c as this may otherwise result in too much transfer of water during
c rime-splintering
            mpracs(k) = min(mpracs(k),qr3d(k)/dt)
            else
	    endif if_rain_snow
	  endif if_rain_snow_exists
	  if_snow_exists: if (lqs3d(k) .and. lcold(k)) then
c.......................................................................
c.......................................................................
c rime-splintering
c hallet-mossop (1974)
c number of splinters formed is based on mass of rimed water
c dum1 = mass of individual splinters
	    if_snow_mult: if (lsnow_mult .and. (mpsacws(k).gt.0) ) then
              if_cotton:
     1          if (tk3d(k).lt.270.16 .and. tk3d(k).gt.265.16) then
	         dum1 = 4.0d0/3.0d0*pi*(5.0d-6)**3*rhoi
	         if (tk3d(k).gt.270.16d0) then
	            fmult = 0.0d0
	         else if (tk3d(k).le.270.16.and.tk3d(k).gt.268.16) 
     1		 then 
	            fmult = (270.16-tk3d(k))/2.0d0
	         else if (tk3d(k).ge.265.16.and.tk3d(k).le.268.16) 
     1		 then
	            fmult = (tk3d(k)-265.16)/3.0d0
	         else if (tk3d(k).lt.265.16) then
	            fmult = 0.d0
	      end if

c 1000 is to convert from kg to g
c riming and splintering from accreted raindrops
	      if (mpracs(k).gt.0.) then
                nmultr(k) = 35.0d4*mpracs(k)*fmult*1.0d3
                mmultr(k) = nmultr(k)*dum1 
c constrain so that transfer of mass from snow to ice cannot be more mass
c than was rimed onto snow
                mmultr(k) = min(mmultr(k),mpracs(k))
	        mpracs(k) = mpracs(k)-mmultr(k)
              endif
              endif if_cotton
            else
	    endif if_snow_mult
	  endif if_snow_exists
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_make_rain_by_snow_accretion
C
C==============================================================================
      LOGICAL FUNCTION hugh_make_drop_by_snow_accretion(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      real*8  :: dc0,ds0,eci,dum,uns,dum1,fmult
      logical :: ldrop_snow,lsnow_mult
        CHARACTER*34 :: sname='hugh_make_drop_by_snow_accretion: '
        la     = .true.
        ldrop_snow = .true.!; ldrop_snow = .false.
        lsnow_mult = .true.; lsnow_mult = .false.

c Process description
        altitude_loop: do k=1,mx0
          npsacws(k)=0.0d0; mpsacws(k)=0.0d0
          nmults(k)=0.0d0; nmults(k)=0.0d0
	  if_drop_snow_exists: 
     *      if (lqc3d(k) .AND. lqs3d(k) .and. lcold(k)) then
c.......................................................................
c accretion of cloud droplets onto snow/graupel
c here use continuous collection equation with
c simple gravitational collection kernel ignoring
c collisions between droplets/cloud ice
c since minimum size ice particle for accretion is 50 - 250 micron
c
	    if_drop_snow: if (ldrop_snow) then

              dc0 = (pgam(k)+1.0d0)/lamc(k)
              ds0 = 2.0d0/lams(k)
	      uns = asn(k)*gamma(1.0d0+bs)/lams(k)**bs
              dum = dc0*dc0*uns*rhow/(9.0d0*mu(k)*ds0)
              eci = dum*dum/((dum+0.4d0)*(dum+0.4d0))
              eci = max(eci,0.d0)
              eci = min(eci,1.d0)
	      mpsacws(k) = pi/4.0d0*asn(k)*qc3d(k)*rho(k)*
     1             n0s(k)*Eci*gamma(bs+3.0d0)/
     1             lams(k)**(bs+3.0d0)	
	      npsacws(k) = pi/4.0d0*asn(k)*nc3d(k)*rho(k)*
     1             n0s(k)*Eci*gamma(bs+3.0d0)/
     1             lams(k)**(bs+3.0d0)
            else
	    endif if_drop_snow
	  endif if_drop_snow_exists
c
	  if_snow_exists: if (lqs3d(k) .and. lcold(k)) then
c.......................................................................
c.......................................................................
c rime-splintering
c hallet-mossop (1974)
c number of splinters formed is based on mass of rimed water
c
c dum1 = mass of individual splinters
c
	    if_snow_mult: if (lsnow_mult .and. (mpsacws(k).gt.0) ) then
              if_cotton:
     1          if (tk3d(k).lt.270.16 .and. tk3d(k).gt.265.16) then

	         dum1 = 4.0d0/3.0d0*pi*(5.0d-6)**3*rhoi
 
	         if (tk3d(k).gt.270.16d0) then
	            fmult = 0.0d0
	         else if (tk3d(k).le.270.16.and.tk3d(k).gt.268.16) 
     1		 then 
	            fmult = (270.16-tk3d(k))/2.0d0
	         else if (tk3d(k).ge.265.16.and.tk3d(k).le.268.16) 
     1		 then
	            fmult = (tk3d(k)-265.16)/3.0d0
	         else if (tk3d(k).lt.265.16) then
	            fmult = 0.d0
	      end if

c 1000 is to convert from kg to g

c splintering from droplets accreted onto snow
	      if (mpsacws(k).gt.0.) then
                nmults(k) = 35.0d4*mpsacws(k)*fmult*1.0d3
                mmults(k) = nmults(k)*dum1 

c constrain so that transfer of mass from snow to ice cannot be more mass
c than was rimed onto snow
                mmults(k) = min(mmults(k),mpsacws(k))
	        mpsacws(k) = mpsacws(k)-mmults(k)
              endif
              endif if_cotton
            else
	    endif if_snow_mult
	  endif if_snow_exists
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_make_drop_by_snow_accretion
C
C========================================================================
C Include file: blk_add_ice_melting.f
C========================================================================
C==============================================================================
      LOGICAL FUNCTION hugh_make_snow_melting(dtmic,mx0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: lmeltSR,lmeltSE
        CHARACTER*30 :: sname='hugh_make_snow_melting: '
        la     = .true.
        lmeltSR = .true.; lmeltSR = .false.
        lmeltSE = .true.; lmeltSE = .false.
c.......................................................................
c melting of snow
c snow may persits above freezing, formula from rutledge and hobbs, 1984
        altitude_loop: do k=1,mx0
          nsmltr(k)=0.0d0; msmltr(k)=0.0d0
          nsmlts(k)=0.0d0; msmlts(k)=0.0d0
	  if_snow_exists: if (lqs3d(k) .and. lwarm(k) ) then
c if water supersaturation, snow melts to form rain
            if_meltSR: if(lmeltSR) then
	      msmltr(k)=2.*pi*n0s(k)*kap(k)*(tf-tk3d(k))/
     1               xlf(k)*rho(k)*(f1s/(lams(k)*lams(k))+
     1               f2s*(asn(k)*rho(k)/mu(k))**0.5*
     1               sc(k)**(1./3.)*gamma(5./2.+bs/2.)/
     1              (lams(k)**(5./2.+bs/2.)))*4./(2.*pi)
              nsmltr(k) = 0.0d0
            endif if_meltSR
            if_meltSE: if(lmeltSE .and. qvqvs(k) .lt. 1.0d0) then
	      msmlts(k) = 2.*pi*n0s(k)*(qv3d(k)-qvs(k))*epss(k)/ab(k)
     1             *4./(2.*pi)
              msmlts(k) = max(msmlts(k),msmltr(k))
              nsmlts(k) = 0.0d0
            endif if_meltSE
            msmltr(k) = msmltr(k)-msmlts(k)
          else
	  endif if_snow_exists
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_make_snow_melting
C
C==============================================================================
C
      LOGICAL FUNCTION hugh_make_crys_melting(dtmic,mx0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: lmeltC,lmeltS
        CHARACTER*30 :: sname='hugh_make_crys_melting: '
        la     = .true.
        lmeltC = .true.; lmeltC = .false.
        lmeltS = .true.; lmeltS = .false.
c Process description
        altitude_loop: do k=1,mx0
          nmltic(k)=0.0d0; mmltic(k)=0.0d0
          nmltii(k)=0.0d0; mmltii(k)=0.0d0
	  if_crys_exists: if (lqc3d(k) .and. lwarm(k)) then
            if_meltC: if(lmeltC) then
            endif if_meltC
            if_meltS: if(lmeltS) then
            endif if_meltS
	  endif if_crys_exists
        enddo altitude_loop
c
        la=.NOT.lcheck
      END FUNCTION hugh_make_crys_melting
C
C========================================================================
C Include file: blk_add_supsat_timescale.f
C========================================================================
      logical function make_drop_timescale(dt0,mx0) 
     *  result(la)
      implicit none 
        real*8,  intent(in) :: dt0
        integer, intent(in) :: mx0
c Local
c counting variables
      integer k !,i,j
      character*21 :: sname='make_drop_timescale: '
        la     = .true.
c drop
!      do i = 1,mix
!         do j = 1,mjx
       altitude_loop: do k = 1,mx0
	   epsc(k) = 0.0d0
	 if (.not. lqc3d(k)) cycle  altitude_loop
	   epsc(k) = 0.0d0
       enddo  altitude_loop
!          end do ! j loop
!        end do  ! i loop
c
        la=.NOT.lcheck
      end function make_drop_timescale
C
C==============================================================================
      logical function make_rain_timescale(dt0,mx0) 
     *  result(la)
      implicit none 
        real*8,  intent(in) :: dt0
        integer, intent(in) :: mx0
c Local
c counting variables
      integer k !,i,j
      character*21 :: sname='make_rain_timescale: '
        la     = .true.
c rain
c ventilation for rain
!      do i = 1,mix
!         do j = 1,mjx
       altitude_loop: do k = 1,mx0
	 epsr(k) = 0.0d0
	 if (.not. lqr3d(k)) cycle  altitude_loop

           epsr(k) = 2.*pi*n0rr(k)*rho(k)*Dv(k)*
     1              (f1r/(lamr(k)*lamr(k))+
     1               f2r*(arn(k)*rho(k)/mu(k))**0.5*
     1               sc(k)**(1./3.)*gamma(5./2.+br/2.)/
     1              (lamr(k)**(5./2.+br/2.)))

       enddo  altitude_loop
!          end do ! j loop
!        end do  ! i loop
c
        la=.NOT.lcheck
c
      end function make_rain_timescale
C
C==============================================================================
      logical function make_crys_timescale(dt0,mx0) 
     *  result(la)
      implicit none 
        real*8,  intent(in) :: dt0
        integer, intent(in) :: mx0
c Local
c counting variables
      integer k !,i,j
      character*21 :: sname='make_crys_timescale: '
        la     = .true.
c crys
c no ventilation for cloud ice
!      do i = 1,mix
!         do j = 1,mjx
       altitude_loop: do k = 1,mx0
	 epsi(k) = 0.0d0
	 if (.not. lqi3d(k)) cycle  altitude_loop

	   epsi(k) = 2.*pi*n0i(k)*rho(k)*Dv(k)
     1              /(lami(k)*lami(k))

       enddo  altitude_loop
!          end do ! j loop
!        end do  ! i loop
c
        la=.NOT.lcheck
c
      end function make_crys_timescale
C
C==============================================================================
      logical function make_snow_timescale(dt0,mx0) 
     *  result(la)
      implicit none 
        real*8,  intent(in) :: dt0
        integer, intent(in) :: mx0
c Local
c counting variables
      integer k !,i,j
      character*21 :: sname='make_snow_timescale: '
        la     = .true.
c snow
c ventilation for snow
c
!      do i = 1,mix
!         do j = 1,mjx
       altitude_loop: do k = 1,mx0
	 epss(k) = 0.0d0
	 if (.not. (lqs3d(k).and.lcold(k))) cycle  altitude_loop

           epss(k) = 2.*pi*n0s(k)*rho(k)*Dv(k)*
     1              (f1s/(lams(k)*lams(k))+
     1               f2s*(asn(k)*rho(k)/mu(k))**0.5*
     1               sc(k)**(1./3.)*gamma(5./2.+bs/2.)/
     1              (lams(k)**(5./2.+bs/2.)))

       enddo  altitude_loop
!          end do ! j loop
!        end do  ! i loop
c
        la=.NOT.lcheck
      end function make_snow_timescale
C
C==============================================================================
      logical function make_xxxx_timescale(dt0,tag,mx0) result(la)
      implicit none 
        real*8,  intent(in) :: dt0
        character (len = *),INTENT (in)       :: tag
        integer, intent(in) :: mx0
c Local
        real*8     :: dtmic
c counting variables
      integer      :: i,j,k
      character*21 :: sname='make_xxxx_timescale: '
        la     = .true.
        dtmic=dt0
        case_tag: select case(trim(tag))
c Make droplet distribution
          case('drop')
            la= make_drop_timescale(dtmic,mx0) 
c Make rain distribution
          case('rain')
            la= make_rain_timescale(dtmic,mx0) 
c Make crystal distribution
          case('crys')
            la= make_crys_timescale(dtmic,mx0) 
c Make rain distribution
          case('snow')
            la= make_rain_timescale(dtmic,mx0) 
c Make liquid phase distributions
          case('liquid')
            la= make_drop_timescale(dtmic,mx0) 
            la= make_rain_timescale(dtmic,mx0) 
c Make solid phase distributions
          case('solid')
            la= make_crys_timescale(dtmic,mx0) 
            la= make_snow_timescale(dtmic,mx0) 
c Make cloud water and cloud ice 
          case('cloud')
            la= make_drop_timescale(dtmic,mx0) 
            la= make_crys_timescale(dtmic,mx0) 
c Make precipitating water and  precipitating ice 
          case('precip')
            la= make_rain_timescale(dtmic,mx0) 
            la= make_snow_timescale(dtmic,mx0) 
c Make all distributions
          case('all')
            la= make_drop_timescale(dtmic,mx0) 
            la= make_rain_timescale(dtmic,mx0) 
            la= make_crys_timescale(dtmic,mx0) 
            la= make_snow_timescale(dtmic,mx0) 
          case default
            msg="Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select case_tag
c
        la=.NOT.lcheck
      end function make_xxxx_timescale
C
C========================================================================
C Include file: blk_add_ice_condensation.f
C==============================================================================
      LOGICAL FUNCTION hugh_make_ice_condensation(dtmic,mx0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      real*8  :: dum,fudgef,sum_dep
      logical :: lcondI,lcondS
        CHARACTER*34 :: sname='hugh_make_ice_condensation: '
        la     = .true.
        lcondI = .true.; lcondI = .false.
        lcondS = .true.; lcondS = .false.
c.......................................................................
c sublimation/deposition of crys/snow below freezing
        altitude_loop: do k=1,mx0
          ncondi(k)=0.0d0; mcondi(k)=0.0d0
          nconds(k)=0.0d0; mconds(k)=0.0d0
c sublimation/deposition of crys
	  if_condI: if (lcondI .and. lcold(k) ) then
	    mcondi(k) = epsi(k)*(qv3d(k)-qvi(k))/abi(k)
          else

          endif if_condI

c sublimation/deposition of snow
	  if_condS: if (lcondS .and. lcold(k) ) then
	    mconds(k) = epss(k)*(qv3d(k)-qvi(k))/abi(k)
          endif if_condS
c make sure not pushed into ice supersat/subsat
c formula from reisner 2 scheme
           dum = (qv3d(k)-qvi(k))/dt
           fudgef = 0.9999
           sum_dep = mcondi(k)+mconds(k)+mnuccd(k)
           IF( (dum.GT.0. .AND. SUM_DEP.GT.dum*FUDGEF) . 
     +      OR. (dum.LT.0. .AND. SUM_DEP.LT.dum*FUDGEF) ) THEN
               mnuccd(K) = FUDGEF*mnuccd(K)*dum/SUM_DEP
               mcondi(K) = FUDGEF*mcondi(K)*dum/SUM_DEP
               mconds(K) = FUDGEF*mconds(K)*dum/SUM_DEP
           ENDIF                                         
c make sure mixing ratio not negative
	    mcondi(k) = max(-qi3d(k)/dt,mcondi(k))
	    mconds(k) = max(-qs3d(k)/dt,mconds(k))
        enddo altitude_loop
        la=.NOT.lcheck
      END FUNCTION hugh_make_ice_condensation
C
C========================================================================
C Include file: blk_add_water_condensation.f
C========================================================================
      LOGICAL FUNCTION hugh_make_water_condensation(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      real*8  :: dumt,dumqv,dumqc,dumqss,dums,dumres
      logical :: lcondC,lcondR
        CHARACTER*36 :: sname='hugh_make_water_condensation: '
        la     = .true.
        lcondC = .true.!; lcondC = .false.
        lcondR = .true.!; lcondR = .false.
c.......................................................................
c saturation adjustment 
        altitude_loop: do k=1,mx0
          ncondc(k)=0.0d0; mcondc(k)=0.0d0
          ncondr(k)=0.0d0; mcondr(k)=0.0d0
c no condensation onto rain, only evap
	  if_condR: if (lcondR .and. qv3d(k).lt.qvs(k) ) then
	    mcondr(k) = epsr(k)*(qv3d(k)-qvs(k))/ab(k)
	    mcondr(k) = max(-qr3d(k)/dt,mcondr(k))
          endif if_condR

c now calculate saturation adjustment to condense extra vapor above
c water saturation: mcondr(k) <= 0 !!!
	  if_condC: if (lcondC) then
            dumt   = tk3d(k) + dtmic*(tk3dten(k)+mcondr(k)*xxlv(k)/cp)
            dumqv  = qv3d(k) + dtmic*(qv3dten(k)-mcondr(k))
            dumqss = 0.622*polysvp(dumt,0)/
     1          (pp3d(k)-polysvp(dumt,0))
            dumqc  = qc3d(k)+dtmic*qc3dten(k)
            dumqc  = max(dumqc,0.0d0)
c saturation adjustment for liquid
            dums   = dumqv-dumqss
c            mcondc(k)=dums/(1.+xxlv(k)**2*dumqss/(cp*rv*dumt**2))/dtmic
            dumres=dums/(1.+xxlv(k)**2*dumqss/(cp*rv*dumt**2))	   
c            if(mcondc(k)*dtmic+dumqc.lt.0.) then
            if(dumres+dumqc.lt.0.) then
              mcondc(k) = -dumqc/dtmic
            else
              mcondc(k) = dumres/dtmic
            endif
          else
          endif if_condC
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_make_water_condensation
C
C========================================================================
C Include file: blk_add_sedimentation.f
C========================================================================
      LOGICAL FUNCTION hugh_make_drop_sedimentation(dtmic,mx0) 
     *  RESULT (la)
      use TimeConstants_mod, only: SECONDS_PER_HOUR
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k,n,nstep=1
      real*8  :: dtsed,lammax,lammin,rgvm
      real*8  :: unc,umc,faltndmc,faltndnc,dlamc
      logical :: lsedimC
        CHARACTER*35 :: sname='hugh_make_drop_sedimentation: '
        la     = .true.
        lsedimC = .true.; lsedimC = .false.
        nstep=1
        nvtermc=0.0d0; mvtermc=0.0d0
        nsedimc=0.0d0; msedimc=0.0d0
c
c Drop sedimentation
        if(.not. lsedimC) return
          altitude_loop0: do k=1,mx0
            nvtermc(k)=0.0d0; mvtermc(k)=0.0d0
            nsedimc(k)=0.0d0; msedimc(k)=0.0d0
            unc=0.0d0;umc=0.0d0
	    fmc(k)=0.0d0;fnc(k)=0.0d0
	    dumfmc(k)=0.0d0;dumfnc(k)=0.0d0
	    if_drop_exists: if (lqc3d(k)) then
	      dumfmc(k) = qc3d(k)+qc3dten(k)*dtmic
	      dumfnc(k) = nc3d(k)+nc3dten(k)*dtmic

c make sure number concentrations are positive
	      dumfnc(k) = max(0.0d0,dumfnc(k))
	      if (dumfmc(k).ge.qsmall) then
c get dummy lamda
	        dlamc = (pi/6.*rhow*dumfnc(k)*gamma(pgam(k)+4.)/
     1            (dumfmc(k)*gamma(pgam(k)+1.)))**(1./3.)
	        lammin = (pgam(k)+1.)/60.e-6
	        lammax = (pgam(k)+1.)/1.e-6
                dlamc=max(dlamc,lammin)
                dlamc=min(dlamc,lammax)
c calculate number-weighted and mass-weighted terminal fall speeds
	        unc =  acn(k)*gamma(1.+bc+pgam(k))/
     1            (dlamc**bc*gamma(pgam(k)+1.))
	        umc = acn(k)*gamma(4.+bc+pgam(k))/
     1            (dlamc**bc*gamma(pgam(k)+4.))
              endif

              nvtermc(k)  = unc        ! [m/s]
              mvtermc(k)  = umc        ! [m/s]
	      fmc(k) = g*rho(k)*umc   ! [Pa/s]
	      fnc(k) = g*rho(k)*unc   ! [Pa/s]
c
	    endif if_drop_exists
c calculate number of split time steps

	    rgvm = max(fmc(k),fnc(k))
! VVT changed IFIX -> INT (generic function)
	    nstep = max(int(rgvm*dtmic/pdel(k)+1.),nstep) 
c
          enddo altitude_loop0
c
c Time splitting 
c
        dtsed=dtmic/nstep
	time_splitting: do n = 1,nstep
	  faloutmc = fmc*dumfmc
	  faloutnc = fnc*dumfnc

c top of model
	  k = mx0
	  faltndmc = faloutmc(k)/pdel(k)
          faltndnc = faloutnc(k)/pdel(k)
c add fallout terms to eulerian tendencies	
c	  qc3dten(k) = qc3dten(k)-faltndmc/nstep
c          nc3dten(k) = nc3dten(k)-faltndnc/nstep
          nsedimc(k) = nsedimc(k)-faltndnc/nstep
          msedimc(k) = msedimc(k)-faltndnc/nstep

	  dumfmc(k) = dumfmc(k)-faltndmc*dtsed
          dumfnc(k) = dumfnc(k)-faltndnc*dtsed
c	
          altitude_loop: do k=mx0-1,1,-1
	    faltndmc = (faloutmc(k)-faloutmc(k+1))/pdel(k)
	    faltndnc = (faloutnc(k)-faloutnc(k+1))/pdel(k)
c add fallout terms to eulerian tendencies
c	    qc3dten(k) = qc3dten(k)-faltndmc/nstep
c           nc3dten(k) = nc3dten(k)-faltndnc/nstep
            nsedimc(k) = nsedimc(k)-faltndnc/nstep
            msedimc(k) = msedimc(k)-faltndnc/nstep

	    dumfmc(k) = dumfmc(k)-faltndmc*dtsed
            dumfnc(k) = dumfnc(k)-faltndnc*dtsed

            fnc(k)=max(fnc(k)/pdel(k),fnc(k+1)/
     1        pdel(k+1))*pdel(k)
            fmc(K)=MAX(fmc(K)/pdel(k),fmc(K+1)/
     1        pdel(K+1))*pdel(k)
c
          enddo altitude_loop
c
c get precipitation and snowfall rate at the surface (mm/hour)
c
c bottom of model
	  k = 1
	  precrt = precrt+faloutmc(k)
     1                /g/10./nstep*10.*SECONDS_PER_HOUR
c
        enddo time_splitting
c
        la=.NOT.lcheck
      END FUNCTION hugh_make_drop_sedimentation
C
C==============================================================================
      LOGICAL FUNCTION hugh_make_crys_sedimentation(dtmic,mx0) 
     *  RESULT (la)
      use TimeConstants_mod, only: SECONDS_PER_HOUR
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k,n,nstep=1
      real*8  :: dtsed,lammax,lammin,rgvm,pgami
      real*8  :: uni,umi,faltndmi,faltndni,dlami
      logical :: lsedimI
        CHARACTER*35 :: sname='hugh_make_crys_sedimentation: '
        la     = .true.
        lsedimI = .true.; lsedimI = .false.
        nstep=1;pgami=5.
        nvtermi=0.0d0; mvtermi=0.0d0
        nsedimi=0.0d0; msedimi=0.0d0
c Crys sedimentation
        if(.not. lsedimI) return
          altitude_loop0: do k=1,mx0
            nvtermi(k)=0.0d0; mvtermi(k)=0.0d0
            nsedimi(k)=0.0d0; msedimi(k)=0.0d0
            uni=0.0d0;umi=0.0d0
	    fmi(k)=0.0d0;fni(k)=0.0d0
	    dumfmi(k)=0.0d0;dumfni(k)=0.0d0
	    if_crys_exists: if (lqi3d(k)) then
	      dumfmi(k) = qi3d(k)+qi3dten(k)*dtmic
	      dumfni(k) = ni3d(k)+ni3dten(k)*dtmic

c make sure number concentrations are positive
	      dumfni(k) = max(0.0d0,dumfni(k))
	      if (dumfmi(k).ge.qsmall) then
c get dummy lamda
	        dlami = (gamma(1.+di)*ci*
     1            dumfni(k)/dumfmi(k))**(1./di)
	        lammax = 1./1.e-6
	        lammin = 1./(2.*dcs)
                dlami=max(dlami,lammin)
                dlami=min(dlami,lammax)

c calculate number-weighted and mass-weighted terminal fall speeds
	        uni =  ain(k)*gamma(1.+bi+pgami)/
     1            (dlami**bi*gamma(pgami+1.))
	        umi = ain(k)*gamma(4.+bi+pgami)/
     1            (dlami**bi*gamma(pgami+4.))
              endif

c set realistic limits on fallspeed
              umi=min(umi,1.2d0)
              uni=min(uni,1.2d0)
              nvtermi(k)  = uni        ! [m/s]
              mvtermi(k)  = umi        ! [m/s]
	      fmi(k) = g*rho(k)*umi    ! [Pa/s]
	      fni(k) = g*rho(k)*uni    ! [Pa/s]
c
	    endif if_crys_exists
c calculate number of split time steps

	    rgvm = max(fmi(k),fni(k))
! VVT changed IFIX -> INT (generic function)
	    nstep = max(int(rgvm*dtmic/pdel(k)+1.),nstep) 
c
          enddo altitude_loop0
c
c Time splitting 
        dtsed=dtmic/nstep
	time_splitting: do n = 1,nstep
	  faloutmi = fmi*dumfmi
	  faloutni = fni*dumfni

c top of model
	  k = mx0
	  faltndmi = faloutmi(k)/pdel(k)
          faltndni = faloutni(k)/pdel(k)
c add fallout terms to eulerian tendencies	
c	  qi3dten(k) = qi3dten(k)-faltndmi/nstep
c         ni3dten(k) = ni3dten(k)-faltndni/nstep
          nsedimi(k) = nsedimi(k)-faltndni/nstep
          msedimi(k) = msedimi(k)-faltndni/nstep

	  dumfmi(k) = dumfmi(k)-faltndmi*dtsed
          dumfni(k) = dumfni(k)-faltndni*dtsed
c	
          altitude_loop: do k=mx0-1,1,-1
	    faltndmi = (faloutmi(k)-faloutmi(k+1))/pdel(k)
	    faltndni = (faloutni(k)-faloutni(k+1))/pdel(k)
c add fallout terms to eulerian tendencies
c	    qi3dten(k) = qi3dten(k)-faltndmi/nstep
c           ni3dten(k) = ni3dten(k)-faltndni/nstep
            nsedimi(k) = nsedimi(k)-faltndni/nstep
            msedimi(k) = msedimi(k)-faltndni/nstep

	    dumfmi(k) = dumfmi(k)-faltndmi*dtsed
            dumfni(k) = dumfni(k)-faltndni*dtsed

            fni(k)=max(fni(k)/pdel(k),fni(k+1)/
     1        pdel(k+1))*pdel(k)
            fmi(K)=MAX(fmi(K)/pdel(k),fmi(K+1)/
     1        pdel(K+1))*pdel(k)
c
          enddo altitude_loop
c
c get precipitation and snowfall rate at the surface (mm/hour)
c bottom of model
	  k = 1
	  precrt = precrt+faloutmi(k)
     1                /g/10./nstep*10.*SECONDS_PER_HOUR
c
	  snowrt = snowrt+faloutmi(k)
     1                /g/10./nstep*10.*SECONDS_PER_HOUR
c
        enddo time_splitting
        la=.NOT.lcheck
      END FUNCTION hugh_make_crys_sedimentation
C
C==============================================================================
      LOGICAL FUNCTION hugh_make_rain_sedimentation(dtmic,mx0) 
     *  RESULT (la)
      use TimeConstants_mod, only: SECONDS_PER_HOUR
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k,n,nstep=1
      real*8  :: dtsed,lammax,lammin,rgvm
      real*8  :: unr,umr,faltndmr,faltndnr,dlamr
      logical :: lsedimR
        CHARACTER*35 :: sname='hugh_make_rain_sedimentation: '
        la     = .true.
        lsedimR = .true.; lsedimR = .false.
        nstep=1
        nvtermr=0.0d0; mvtermr=0.0d0
        nsedimr=0.0d0; msedimr=0.0d0
c Crys sedimentation
        if(.not. lsedimR) return
          altitude_loop0: do k=1,mx0
            nvtermr(k)=0.0d0; mvtermr(k)=0.0d0
            nsedimr(k)=0.0d0; msedimr(k)=0.0d0
            unr=0.0d0;umr=0.0d0
	    fmr(k)=0.0d0;fnr(k)=0.0d0
	    dumfmr(k)=0.0d0;dumfnr(k)=0.0d0
	    if_rain_exists: if (lqr3d(k)) then
	      dumfmr(k) = qr3d(k)+qr3dten(k)*dtmic
	      dumfnr(k) = nr3d(k)+nr3dten(k)*dtmic

c make sure number concentrations are positive
	      dumfnr(k) = max(0.0d0,dumfnr(k))
	      if (dumfmr(k).ge.qsmall) then
c get dummy lamda
	        dlamr = (pi*rhow*dumfnr(k)/dumfmr(k))**(1./3.)
	        lammax = 1./5.d-6
	        lammin = 1./5000.d-6
                dlamr=max(dlamr,lammin)
                dlamr=min(dlamr,lammax)

c calculate number-weighted and mass-weighted terminal fall speeds
	        unr = arn(k)*gamma(1.+br)/dlamr**br
	        umr = arn(k)*gamma(4.+br)/(6.*dlamr**br)
              endif

c set realistic limits on fallspeed
              umr=min(umr,9.1d0)
              unr=min(unr,9.1d0)
              nvtermr(k)  = unr        ! [m/s]
              mvtermr(k)  = umr        ! [m/s]
	      fmr(k) = g*rho(k)*umr    ! [Pa/s]
	      fnr(k) = g*rho(k)*unr    ! [Pa/s]
c
	    endif if_rain_exists
c calculate number of split time steps
	    rgvm = max(fmr(k),fnr(k))
! VVT changed IFIX -> INT (generic function)
	    nstep = max(int(rgvm*dtmic/pdel(k)+1.),nstep) 
c
          enddo altitude_loop0
c
c Time splitting 
        dtsed=dtmic/nstep
	time_splitting: do n = 1,nstep
	  faloutmr = fmr*dumfmr
	  faloutnr = fnr*dumfnr

c top of model
	  k = mx0
	  faltndmr = faloutmr(k)/pdel(k)
          faltndnr = faloutnr(k)/pdel(k)
c add fallout terms to eulerian tendencies	
c	  qr3dten(k) = qr3dten(k)-faltndmr/nstep
c         nr3dten(k) = nr3dten(k)-faltndnr/nstep
          nsedimr(k) = nsedimr(k)-faltndnr/nstep
          msedimr(k) = msedimr(k)-faltndnr/nstep

	  dumfmr(k) = dumfmr(k)-faltndmr*dtsed
          dumfnr(k) = dumfnr(k)-faltndnr*dtsed
c	
          altitude_loop: do k=mx0-1,1,-1
	    faltndmr = (faloutmr(k)-faloutmr(k+1))/pdel(k)
	    faltndnr = (faloutnr(k)-faloutnr(k+1))/pdel(k)
c add fallout terms to eulerian tendencies
c	    qr3dten(k) = qr3dten(k)-faltndmr/nstep
c           nr3dten(k) = nr3dten(k)-faltndnr/nstep
            nsedimr(k) = nsedimr(k)-faltndnr/nstep
            msedimr(k) = msedimr(k)-faltndnr/nstep

	    dumfmr(k) = dumfmr(k)-faltndmr*dtsed
            dumfnr(k) = dumfnr(k)-faltndnr*dtsed

            fnr(k)=max(fnr(k)/pdel(k),fnr(k+1)/
     1        pdel(k+1))*pdel(k)
            fmr(K)=MAX(fmr(K)/pdel(k),fmr(K+1)/
     1        pdel(K+1))*pdel(k)
c
          enddo altitude_loop
c
c get precipitation and snowfall rate at the surface (mm/hour)
c bottom of model
	  k = 1
	  precrt = precrt+faloutmr(k)
     1                /g/10./nstep*10.*SECONDS_PER_HOUR
c
        enddo time_splitting
        la=.NOT.lcheck
      END FUNCTION hugh_make_rain_sedimentation
C
C==============================================================================
      LOGICAL FUNCTION hugh_make_snow_sedimentation(dtmic,mx0) 
     *  RESULT (la)
      use TimeConstants_mod, only: SECONDS_PER_HOUR
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k,n,nstep=1
      real*8  :: dtsed,lammax,lammin,rgvm
      real*8  :: uns,ums,faltndms,faltndns,dlams
      logical :: lsedimS
        CHARACTER*35 :: sname='hugh_make_snow_sedimentation: '
        la     = .true.
        lsedimS = .true.; lsedimS = .false.
        nstep=1
        nvterms=0.0d0; mvterms=0.0d0
        nsedims=0.0d0; msedims=0.0d0
c Crys sedimentation
        if(.not. lsedimS) return
          altitude_loop0: do k=1,mx0
            nvterms(k)=0.0d0; mvterms(k)=0.0d0
            nsedims(k)=0.0d0; msedims(k)=0.0d0
            uns=0.0d0;ums=0.0d0
	    fms(k)=0.0d0;fns(k)=0.0d0
	    dumfms(k)=0.0d0;dumfns(k)=0.0d0
	    if_snow_exists: if (lqs3d(k)) then
	      dumfms(k) = qs3d(k)+qs3dten(k)*dtmic
	      dumfns(k) = ns3d(k)+ns3dten(k)*dtmic
c make sure number concentrations are positive
	      dumfns(k) = max(0.0d0,dumfns(k))
	      if (dumfms(k).ge.qsmall) then
c get dummy lamda
	        dlams = (gamma(1.+ds)*cs*dumfns(k)/
     1            dumfms(k))**(1./ds)
	        lammax = 1./5.e-6
	        lammin = 1./5000.e-6
                dlams=max(dlams,lammin)
                dlams=min(dlams,lammax)
c calculate number-weighted and mass-weighted terminal fall speeds
	        ums = asn(k)*gamma(4.+bs)/(6.*dlams**bs)
	        uns = asn(k)*gamma(1.+bs)/dlams**bs
              endif
c set realistic limits on fallspeed
              ums=min(ums,1.2d0)
              uns=min(uns,1.2d0)
              nvterms(k)  = uns        ! [m/s]
              mvterms(k)  = ums        ! [m/s]
	      fms(k) = g*rho(k)*ums    ! [Pa/s]
	      fns(k) = g*rho(k)*uns    ! [Pa/s]
	    endif if_snow_exists
c calculate number of split time steps
	    rgvm = max(fms(k),fns(k))
! VVT changed IFIX -> INT (generic function)
	    nstep = max(int(rgvm*dtmic/pdel(k)+1.),nstep) 
          enddo altitude_loop0
c Time splitting 
        dtsed=dtmic/nstep
	time_splitting: do n = 1,nstep
	  faloutms = fms*dumfms
	  faloutns = fns*dumfns
c top of model
	  k = mx0
	  faltndms = faloutms(k)/pdel(k)
          faltndns = faloutns(k)/pdel(k)
c add fallout terms to eulerian tendencies	
c	  qs3dten(k) = qs3dten(k)-faltndms/nstep
c         ns3dten(k) = ns3dten(k)-faltndns/nstep
          nsedims(k) = nsedims(k)-faltndns/nstep
          msedims(k) = msedims(k)-faltndns/nstep
	  dumfms(k) = dumfms(k)-faltndms*dtsed
          dumfns(k) = dumfns(k)-faltndns*dtsed
          altitude_loop: do k=mx0-1,1,-1
	    faltndms = (faloutms(k)-faloutms(k+1))/pdel(k)
	    faltndns = (faloutns(k)-faloutns(k+1))/pdel(k)
c add fallout terms to eulerian tendencies
c	    qs3dten(k) = qs3dten(k)-faltndms/nstep
c           ns3dten(k) = ns3dten(k)-faltndns/nstep
            nsedims(k) = nsedims(k)-faltndns/nstep
            msedims(k) = msedims(k)-faltndns/nstep

	    dumfms(k) = dumfms(k)-faltndms*dtsed
            dumfns(k) = dumfns(k)-faltndns*dtsed

            fns(k)=max(fns(k)/pdel(k),fns(k+1)/
     1        pdel(k+1))*pdel(k)
            fms(K)=MAX(fms(K)/pdel(k),fms(K+1)/
     1        pdel(K+1))*pdel(k)
          enddo altitude_loop
c
c get precipitation and snowfall rate at the surface (mm/hour)
c bottom of model
	  k = 1
	  precrt = precrt+faloutms(k)
     1                /g/10./nstep*10.*SECONDS_PER_HOUR
	  snowrt = snowrt+faloutms(k)
     1                /g/10./nstep*10.*SECONDS_PER_HOUR
        enddo time_splitting
        la=.NOT.lcheck
      END FUNCTION hugh_make_snow_sedimentation
C
C==============================================================================
      logical function execute_bulk2m_sedimentation_gcm0(dt0,tag,k0) 
     *  result(la)
      implicit none 
        real*8,  intent(in) :: dt0
        character (len = *),INTENT (in)       :: tag
        integer, intent(in) :: k0
c Local
        real*8     :: dtmic
c counting variables
      integer      :: k,mx0!,i,j
      character*25 :: sname='execute_bulk2m_sedimentation_gcm: '
        la     = .true.
        dtmic=dt0
        mx0=k0
        case_tag: select case(trim(tag))
c Make droplet distribution
          case('drop')
            la= hugh_drop_sedimentation(dtmic,mx0) 
c Make rain distribution
          case('rain')
c            la= hugh_rain_sedimentation(dtmic,mx0) 
c Make crystal distribution
          case('crys')
            la= hugh_crys_sedimentation(dtmic,mx0) 
c Make rain distribution
          case('snow')
c            la= hugh_rain_sedimentation(dtmic,mx0) 
c Make liquid phase distributions
          case('liquid')
            la= hugh_drop_sedimentation(dtmic,mx0) 
c            la= hugh_rain_sedimentation(dtmic,mx0) 
c Make solid phase distributions
          case('solid')
c            la= hugh_crys_sedimentation(dtmic,mx0) 
c            la= hugh_snow_sedimentation(dtmic,mx0) 
c Make cloud water and cloud ice 
          case('cloud')
            la= hugh_drop_sedimentation(dtmic,mx0) 
c            la= hugh_crys_sedimentation(dtmic,mx0) 
c Make precipitating water and  precipitating ice 
          case('precip')
c            la= hugh_rain_sedimentation(dtmic,mx0) 
c            la= hugh_snow_sedimentation(dtmic,mx0) 
c Make all distributions
          case('all')
            la= hugh_drop_sedimentation(dtmic,mx0) 
c            la= hugh_rain_sedimentation(dtmic,mx0) 
c            la= hugh_crys_sedimentation(dtmic,mx0) 
c            la= hugh_snow_sedimentation(dtmic,mx0) 
          case default
            msg="Case "//trim(tag)//" is not implemented"
            call stop_model(msg,255)             
        end select case_tag
c
        la=.NOT.lcheck
      end function execute_bulk2m_sedimentation_gcm0
C
C==============================================================================
      LOGICAL FUNCTION hugh_drop_sedimentation(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8       :: dtmic
        integer      :: mx0,k=1
        character*25 :: sname='hugh_drop_sedimentation: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
c        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
c     *  ,tk3d,pp3d,qsmall,'drop',mx0)
c Make time & altitude dependent coefficients
c        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
c Make supsat timescale 
c        ldummy=make_drop_timescale(dtmic,mx0)
        ldummy=hugh_make_drop_sedimentation(dtmic,mx0)
        la=.NOT.lcheck
      END FUNCTION hugh_drop_sedimentation
C
C==============================================================================
      LOGICAL FUNCTION hugh_crys_sedimentation(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8       :: dtmic
        integer      :: mx0,k=1
        character*25 :: sname='hugh_crys_sedimentation: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
c        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
c     *  ,tk3d,pp3d,qsmall,'crys',mx0)
c Make time & altitude dependent coefficients
c        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
c
c Make supsat timescale 
c        ldummy=make_drop_timescale(dtmic,mx0)
        ldummy=hugh_make_crys_sedimentation(dtmic,mx0)
        la=.NOT.lcheck
c
      END FUNCTION hugh_crys_sedimentation
C
C========================================================================
C Include file: blk_add_update_tendencies.f
C========================================================================
C
      LOGICAL FUNCTION hugh_balance_tendencies(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      real*8  :: dum,ratio
      logical :: ldrop0
        CHARACTER*25 :: sname='hugh_balance_tendencies: '
        la = .true.
c Balance tendencies
        altitude_loop: do k=1,mx0
         if(lqc3d(k).and.lqi3d(k).and.lqr3d(k)
     1     .and.lqs3d(k)) cycle altitude_loop
         warm_or_cold_balance: if(lcold(k)) then
c.......................................................................
c conservation of water
c This is adopted loosely from MM5 resiner code. However, here we
c only adjust processes that are negative, rather than all processes.
c this section is separated into two parts, if T < 0 C, T > 0 C
c due to different processes that act depending on freezing/above freezing
c if mixing ratios less than qsmall, then no depletion of water
c through microphysical processes, skip conservation
c conservation of qc
	dum = (mprc(k)+mpra(k)+mnuccc(k)+
     1         mpsacws(k)+mmults(k))*dtmic
	if (dum.gt.qc3d(k).and.qc3d(k).ge.qsmall) then
          ratio = qc3d(k)/dum
	  mprc(k) = mprc(k)*ratio
	  mpra(k) = mpra(k)*ratio
	  mnuccc(k)  = mnuccc(k)*ratio
	  mpsacws(k) = mpsacws(k)*ratio
	  mmults(k)  = mmults(k)*ratio
        end if
c conservation of qi
	dum = (-mcondi(k)-mnuccc(k)+mprci(k)+
     1         mprai(k)-mmults(k)-mmultr(k)-mnuccd(k))*dtmic
	if (dum.gt.qi3d(k).and.qi3d(k).ge.qsmall) then
        ratio = qi3d(k)/dum
        mcondi(k) = mcondi(k)*ratio	
        mnuccc(k) = mnuccc(k)*ratio
        mprci(k)  = mprci(k)*ratio
        mprai(k)  = mprai(k)*ratio
        mmults(k) = mmults(k)*ratio
        mmultr(k) = mmultr(k)*ratio
        mnuccd(k) = mnuccd(k)*ratio
        end if
	! this point always works
c conservation of qr
	dum=((mpracs(k)-mcondr(k))+(mmultr(k)-mprc(k))
     1    +(mnuccr(k)+mnucir(k)-mpra(k)))*dtmic
	if (dum.gt.qr3d(k).and.qr3d(k).ge.qsmall) then
        ratio = qr3d(k)/dum
        mcondr(k) = mcondr(k)*ratio
        mprc(k)   = mprc(k)*ratio
        mpra(k)   = mpra(k)*ratio
        mpracs(k) = mpracs(k)*ratio
        mmultr(k) = mmultr(k)*ratio
        mnuccr(k) = mnuccr(k)*ratio
        mnucir(k) = mnucir(k)*ratio

        end if
c conservation of qs
	dum = (-mconds(k)-mpsacws(k)-mprai(k)-mprci(k)-
     1          mpracs(k)-mnuccr(k)-mnucir(k))*dtmic

	if (dum.gt.qs3d(k).and.qs3d(k).ge.qsmall) then
        ratio = qs3d(k)/dum
        mconds(k)  = mconds(k)*ratio
        mpsacws(k) = mpsacws(k)*ratio
        mprai(k)   = mprai(k)*ratio
        mprci(k)   = mprci(k)*ratio
        mpracs(k)  = mpracs(k)*ratio
        mnuccr(k)  = mnuccr(k)*ratio
        mnucir(k)  = mnucir(k)*ratio

       end if
         else ! warm cloud: temperature above freezing
c for cloud ice and snow, only processes operating at T > tf is
c melting/evaporating, which is already conserved during process
c calculation
c conservation of qc
	dum = (mprc(k)+mpra(k))*dtmic
	if (dum.gt.qc3d(k).and.qc3d(k).ge.qsmall) then
        ratio = qc3d(k)/dum
        mprc(k) = mprc(k)*ratio
        mpra(k) = mpra(k)*ratio
        end if
c conservation of qr
c all terms are positive, do not need conservation of rain
c conservation of snow
        dum = (-msmltr(k)-msmlts(k))*dtmic
        if (dum.gt.qs3d(k).and.qs3d(k).ge.qsmall) then
        ratio = qs3d(k)/dum
        msmltr(k) = msmltr(k)*ratio
        msmlts(k) = msmlts(k)*ratio
        endif
         endif warm_or_cold_balance
        enddo altitude_loop
c
        la=.NOT.lcheck
      END FUNCTION hugh_balance_tendencies
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_tkqv_tendencies(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      real*8  :: dum
      logical :: ldrop0,lsnow0
        CHARACTER*30 :: sname='hugh_update_tkqv_tendencies0: '
        la     = .true.
        ldrop0 = .true.!; ldrop0 = .false.
        lsnow0 = .true.; lsnow0 = .false.
c Update tendensies
        altitude_loop: do k=1,mx0
        if_drop0: if(ldrop0) then
          qv3dten(k) = qv3dten(k)-mcondc(k)
          tk3dten(k) = tk3dten(k)+mcondc(k)*xxlv(k)/cp
          qc3dten(k) = qc3dten(k)+mcondc(k)
        endif if_drop0
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c sublimate, melt, or evaporate number concentration
c this formulation assumes 1:1 ratio between mass loss and
c loss of number concentration
        if_snow0: if(lsnow0) then
c	if (mcondc(k).lt.0.) then
c	   dum = mcondc(k)*dtmic/qc3d(k)
c           dum = max(-1.0d0,dum)
c	   nsubc(k) = dum*nc3d(k)/dtmic
c	end if
c	if (mcondi(k).lt.0.) then
c	   dum = mcondi(k)*dtmic/qi3d(k)
c           dum = max(-1.0d0,dum)
c	   nsubi(k) = dum*ni3d(k)/dtmic
c	end if
c
	if (mconds(k).lt.0.) then
	   dum = mconds(k)*dtmic/qs3d(k)
           dum = max(-1.0d0,dum)
	   nsubs(k) = dum*ns3d(k)/dtmic
	end if
	if (mcondr(k).lt.0.) then
	   dum = mcondr(k)*dtmic/qr3d(k)
           dum = max(-1.0d0,dum)
	   nsubr(k) = dum*nr3d(k)/dtmic
	end if

        if (tk3d(k).ge.tf) then
        if (msmlts(k)+msmltr(k).lt.0.) then
	   dum = (msmlts(k)+msmltr(k))*dtmic/qs3d(k)
           dum = max(-1.0d0,dum)
	   nsmlts(k) = dum*ns3d(k)/dtmic
        end if
        if (msmltr(k).lt.0.) then
          dum = msmltr(k)*dtmic/qs3d(k)
          dum = max(-1.0d0,dum)
          nsmltr(k) = dum*ns3d(k)/dtmic
        end if
        end if

c update tendencies
c        nc3dten(k) = nc3dten(k)+nsubc(k)
c        ni3dten(k) = ni3dten(k)+nsubi(k)
        ns3dten(k) = ns3dten(k)+(nsubs(k)+nsmlts(k))
        nr3dten(k) = nr3dten(k)+(nsubr(k)-nsmltr(k))
        endif if_snow0
        enddo altitude_loop
        la=.NOT.lcheck
      END FUNCTION hugh_update_tkqv_tendencies
C
C==============================================================================
C
      LOGICAL FUNCTION hugh_update_tkqv_tendencies0(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: ldrop0
        CHARACTER*30 :: sname='hugh_update_tkqv_tendencies0: '
        la     = .true.
        ldrop0 = .true.; ldrop0 = .false.
c Update qc3dten
        altitude_loop: do k=1,mx0
         warm_or_cold_tkqv: if(lcold(k)) then
          qv3dten(k)=qv3dten(k)+(
c ! change q evap rain:                        mcondr==>pre
     1              -mcondr(k)
c ! change q dep/sub cloud ice:                mcondi==>prd
     1              -mcondi(k)
c ! change q dep/sub snow                      mcondr==>prds
     1              -mconds(k)
c ! change q freezing aerosol (prim ice nuc)   mnuccd
     1              -mnuccd(k)
c igs
c ! change q droplets activation
     1              -mpccn(k)     
     1                          )
          tk3dten(k)=tk3dten(k)+(
c ! change q evap rain:                        mcondr==>pre
     1              -mcondr(k)*xxlv(k)  ! should be *xxls(k)
c ! change q dep/sub cloud ice:                mcondi==>prd
     1              -mcondi(k)*xxls(k)
c ! change q dep/sub snow                      mcondr==>prds
     1              -mconds(k)*xxls(k)
c ! change q freezing aerosol (prim ice nuc)   mnuccd
     1              +mnuccd(k)*xxls(k)
c ! change q droplets accretion by snow:       mpsacws==>psacws
     1              +mpsacws(k)*xlf(k)
c ! change q collection rain by snow:          mpracs==>pracs          
     1              +mpracs(k)*xlf(k)
c ! change q ice mult acc droplets by snow     mmults==>qmults
     1              +mmults(k)*xlf(k)
c change q due to ice mult acc rain by snow    mmultr==>qmultr
     1              +mmultr(k)*xlf(k)
c ! change q due to con. freez droplets:       mnuccc+mnucci==>mnuccc
     1              +mnuccc(k)*xlf(k)
c ! change q due to imm. freez droplets:
     1              +mnucci(k)*xlf(k)
c ! change q due to con rain freez:            mnuccr
     1             +mnuccr(k)*xlf(k)
c igs
c
c ! change q droplets activation
     1              +mpccn(k)*xxlv(k)     
c ! change q due to imm rain freez:            mnucir
     1             +mnucir(k)*xlf(k)
c ! change n cond freezing Meyers (prim ice nuc):new set to zero for M2005
     1              +mnucmd(k)*xxlv(k)  ! because of condensational nucl
c ! change n cont freezing Meyers (prim ice nuc):new set to zero for M2005
     1              +mnucmt(k)*xlf(k)
     1                          )/cp
c
         else ! warm cloud: temperature above freezing
c
          qv3dten(k)=qv3dten(k)+(
c ! change q evap rain:                        mcondr==>pre
     1              -mcondr(k)
c ! chnage q melting snow evaporating:         msmlts==>evpms          
     1              -msmlts(k)                 
c igs
c
c ! change q droplets activation
     1              -mpccn(k)     
     1                          )

          tk3dten(k)=tk3dten(k)+(
c ! change q evap rain:                        mcondr==>pre
     1              +mcondr(k)*xxlv(k)
c ! chnage q melting snow evaporating:         msmlts==>evpms          
     1              +msmlts(k)*xxls(k)                 
c ! change q melting snow to rain              msmltr==>psmlt
     1              +msmltr(k)*xlf(k)
c igs
c
c ! change q droplets activation
     1              +mpccn(k)*xxlv(k)     
     1                          )/cp
c
         endif warm_or_cold_tkqv
c
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_update_tkqv_tendencies0
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_snow_tendencies0(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: lsnow0
        CHARACTER*30 :: sname='hugh_update_snow_tendencies0: '
        la     = .true.
        lsnow0 = .true.; lsnow0 = .false.
c Update qc3dten
        altitude_loop: do k=1,mx0
         warm_or_cold_snow: if(lcold(k)) then
	 qs3dten(k) = qs3dten(k)+(
c ! change q autoconversion cloud ice:         mprci==>prci
     1              +mprci(k)
c ! change q self-collection of snow is ZERO:  msagg
     1              +msagg(k)
c ! change q accretion cloud ice by snow:      mprai==>prai
     1              +mprai(k)
c ! change q droplets accretion by snow:       mpsacws==>psacws
     1              +mpsacws(k)
c ! change q collection rain by snow:          mpracs==>pracs          
     1              +mpracs(k)
c ! change q dep/sub snow                      mcondr==>prds
     1              +mconds(k)
c ! change q due to con rain freez:            mnuccr
     1              +mnuccr(k)
c ! change q homog freezing rain to graupel:   mhfrr==>qhfrr is ZERO
     1              +mhfrr(k)
     1                           )
          qs3dten(k)=qs3dten(k)+(
     1     +0.0d0
     1                          )
c
          ns3dten(k)=ns3dten(k)+(
c ! change n autoconversion cloud ice:         nprci==>nprci
     1              +nprci(k)
c ! change n self-collection of snow:          nsagg
     1              +nsagg(k)
c ! change n due to con rain freez:            nnuccr
     1              +nnuccr(k)
     1                          )
c
         else ! warm cloud: temperature above freezing
c
          qs3dten(k)=qs3dten(k)+(
c ! change q melting snow to rain              msmltr==>psmlt
     1              +msmltr(k)
c ! chnage q melting snow evaporating:         msmlts==>evpms          
     1              +msmlts(k)                 
     1                          )
c
          ns3dten(k)=ns3dten(k)+(
     1     +0.0d0
     1                          )
c
         endif warm_or_cold_snow
c
        enddo altitude_loop

c
        la=.NOT.lcheck
c
      END FUNCTION hugh_update_snow_tendencies0
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_rain_tendencies0(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: lrain0
        CHARACTER*30 :: sname='hugh_update_rain_tendencies0: '
        la     = .true.
        lrain0 = .true.; lrain0 = .false.
c Update qc3dten
        altitude_loop: do k=1,mx0
         warm_or_cold_rain: if(lcold(k)) then
	 qr3dten(k) = qr3dten(k)+(
c ! change q autoconversion droplets:          mprc ==> prc
     1              +mprc(k)
c ! change q self-collection rain is ZERO:     mragg
     1              +mragg(k)
c ! change q droplets accretion by rain:       mpra==>pra
     1              +mpra(k)
c ! change q collection rain by snow:          mpracs==>pracs          
     1             -mpracs(k)
c ! change q evap rain:                        mcondr==>pre
     1             +mcondr(k)
c change q due to ice mult acc rain by snow    mmultr==>qmultr
     1             -mmultr(k)
c ! change q due to con rain freez:            mnuccr
     1             -mnuccr(k)
c ! change q due to imm rain freez:            mnuccr
     1             -mnucir(k)
     1                          )

c
	 nr3dten(k) = nr3dten(k)+(
c ! change n autoconversion droplets:          nprc
     1              +nprc(k)*0.5
c ! change n self-collection droplets:         nragg
     1              +nragg(k)
c ! change n collection rain by snow:          npracs
     1              -npracs(k)
c ! change n due to con rain freez:            nnuccr
     1              -nnuccr(k)
c ! change n due to imm rain freez:            nnuccr
     1              -nnucir(k)
     1                          )
c
         else ! warm cloud: temperature above freezing
c
          qr3dten(k)=qr3dten(k)+(
c ! change q autoconversion droplets:          mprc ==> prc
     1              +mprc(k)
c ! change q self-collection droplets is ZERO: mragg
     1              +mragg(k)
c ! change q droplets accretion by rain:       mpra==>pra
     1              +mpra(k)
c ! change q evap rain:                        mcondr==>pre
     1             +mcondr(k)
c ! change q melting snow to rain              msmltr==>psmlt
     1             -msmltr(k)
     1                          )
c
          nr3dten(k)=nr3dten(k)+(
c ! change n autoconversion droplets:          nprc
     1              +nprc(k)*0.5
c ! change n self-collection droplets:         nragg
     1              +nragg(k)
     1                          )
c
         endif warm_or_cold_rain
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_update_rain_tendencies0
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_crys_tendencies0(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: lcrys0
        CHARACTER*30 :: sname='hugh_update_crys_tendencies0: '
        la     = .true.
        lcrys0 = .true.; lcrys0 = .false.
c Update qc3dten
        altitude_loop: do k=1,mx0
         warm_or_cold_crys: if(lcold(k)) then
          qi3dten(k)=qi3dten(k)+(
c ! change q autoconversion ice    :           mprci ==> prca
     1              -mprci(k)
c ! change q self-collection ice is ZERO:      miagg
     1              +miagg(k)
c ! change q ice accretion by snow:            mprai==>prci
     1              -mprai(k)
c ! change q due to con. freez droplets:       mnuccc+mnucci==>mnuccc
     1              +mnuccc(k)
c ! change q due to imm. freez droplets:
     1              +mnucci(k)
c ! change n cond freezing Meyers (prim ice nuc):new set to zero for M2005
     1              +mnucmd(k)
c ! change n cont freezing Meyers (prim ice nuc):new set to zero for M2005
     1              +mnucmt(k)
c ! change q dep/sub cloud ice:                mcondi==>prd
     1              +mcondi(k)
c ! change q ice mult acc droplets by snow     mmults==>qmults
     1              +mmults(k)
c change q due to ice mult acc rain by snow    mmultr==>qmultr
     1              +mmultr(k)
c ! change q freezing aerosol (prim ice nuc)
     1              +mnuccd(k)
     1                          )
c
          ni3dten(k)=ni3dten(k)+(
c
c ! change n autoconversion ice    :           nprci ==> nprca
     1              -nprci(k)
c ! change n self-collection ice is ZERO:      niagg
     1              +niagg(k)
c ! change n ice accretion by snow:            nprai==>nprci
     1              -nprai(k)
c ! change n due to con. freez droplets:       nnuccc+nnucci==>nnuccc
     1              +nnuccc(k)
c ! change n due to imm. freez droplets:
     1              +nnucci(k)
c ! change n cond freezing Meyers (prim ice nuc):new set to zero for M2005
     1              +nnucmd(k)
c ! change n cont freezing Meyers (prim ice nuc):new set to zero for M2005
     1              +nnucmt(k)
c ! change n dep/sub cloud ice:                ncondi==>nsubi ????
     1              +ncondi(k)
c ! change n ice mult acc droplets by snow     nmults
     1              +nmults(k)
c ! change n due to ice mult acc rain by snow  nmultr
     1              +nmultr(k)
c ! change n freezing aerosol (prim ice nuc)
     1              +nnuccd(k)
     1                          )
c
         else ! warm cloud: temperature above freezing
c
          qi3dten(k)=qi3dten(k)
     1     +0.0d0
c
          ni3dten(k)=ni3dten(k)
     1     +0.0d0
c
         endif warm_or_cold_crys
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_update_crys_tendencies0
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_drop_tendencies0(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: ldrop0
        CHARACTER*30 :: sname='hugh_update_drop_tendencies0: '
        la     = .true.
        ldrop0 = .true.; ldrop0 = .false.
c Update qc3dten
        altitude_loop: do k=1,mx0
         warm_or_cold_drop: if(lcold(k)) then
          qc3dten(k)=qc3dten(k)+(
c ! change q autoconversion droplets:          mprc ==> prc
     1              -mprc(k)
c ! change q self-collection droplets is ZERO: mcagg
     1              +mcagg(k)
c ! change q droplets accretion by rain:       mpra==>pra
     1              -mpra(k)
c ! change q due to con. freez droplets:       mnuccc+mnucci==>mnuccc
     1              -mnuccc(k)
c ! change q due to imm. freez droplets:
     1              -mnucci(k)
c ! change q cond/evap droplets:               mcondc==>pcc
c     1              +mcondc(k)   ! updated in tk_tendencies
c ! change q droplets accretion by snow:       mpsacws==>psacws
     1              -mpsacws(k)
c ! change q due to ice mult droplets/snow:    mmults==>qmults
     1              -mmults(k)
c ! change q droplets activation
     1              +mpccn(k)     
     1                          )
c
          nc3dten(k)=nc3dten(k)+(
c ! change n autoconversion droplets:
     1              -nprc(k)
c ! change n self-collection droplets:
     1              +ncagg(k)
c ! change n droplets accretion by rain:
     1              -npra(k)
c ! change n due to con. freez droplets:       nnuccc+nnucci==>nnuccc
     1              -nnuccc(k)
c ! change n due to imm. freez droplets:
     1              -nnucci(k)
c ! change n cond/evap droplets:               ncondc==>nsubc ????
     1              +ncondc(k)
c ! change n droplets accretion by snow:
     1              -npsacws(k)
c ! change n due to ice mult droplets/snow is ZERO:
     1              -nmults(k)
c ! change n droplets activation
     1              +npccn(k)     
     1                          )
c
         else ! warm cloud: temperature above freezing
c
          qc3dten(k)=qc3dten(k)+(
c
c ! change q autoconversion droplets:          mprc ==> prc
     1              -mprc(k)
c ! change q self-collection droplets is ZERO: mcagg
     1              +mcagg(k)
c ! change q droplets accretion by rain:       mpra==>pra
     1              -mpra(k)
c ! change q droplets activation
     1              +mpccn(k)     
     1                          )
c
          nc3dten(k)=nc3dten(k)+(
c ! change n autoconversion droplets:
     1              -nprc(k)
c ! change n self-collection droplets:
     1              +ncagg(k)
c ! change n droplets accretion by rain:
     1              -npra(k)
c ! change n droplets activation
     1              +npccn(k)     
     1                          )
c
         endif warm_or_cold_drop
c
        enddo altitude_loop
c
c
        la=.NOT.lcheck
c
      END FUNCTION hugh_update_drop_tendencies0
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_drop_tendencies00(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: ldrop0
        CHARACTER*31 :: sname='hugh_update_drop_tendencies00: '
        la     = .true.
        ldrop0 = .true.; ldrop0 = .false.
c Update qc3dten
        altitude_loop: do k=1,mx0
          qc3dten(k)=qc3dten(k)
     1     +0.0d0
          nc3dten(k)=nc3dten(k)
     1     +0.0d0
        enddo altitude_loop
        la=.NOT.lcheck
c
      END FUNCTION hugh_update_drop_tendencies00
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_tkqv_2zero(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: ltkqv2zero
        CHARACTER*24 :: sname='hugh_update_tkqv_2zero: '
        la     = .true.
        ltkqv2zero = .true.; ltkqv2zero = .false.
c Set tk3dten & qv3dten
        altitude_loop: do k=1,mx0
          tk3dten(k)=0.0d0; qv3dten(k)=0.0d0
        enddo altitude_loop
        la=.NOT.lcheck
      END FUNCTION hugh_update_tkqv_2zero
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_drop_2zero(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: ldrop2zero
        CHARACTER*24 :: sname='hugh_update_drop_2zero: '
        la     = .true.
        ldrop2zero = .true.; ldrop2zero = .false.
c Update qc3dten
        altitude_loop: do k=1,mx0
          qc3dten(k)=0.0d0; nc3dten(k)=0.0d0
        enddo altitude_loop
        la=.NOT.lcheck
      END FUNCTION hugh_update_drop_2zero
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_crys_2zero(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: lcrys2zero
        CHARACTER*24 :: sname='hugh_update_crys_2zero: '
        la     = .true.
        lcrys2zero = .true.; lcrys2zero = .false.
c Update qc3dten
        altitude_loop: do k=1,mx0
          qi3dten(k)=0.0d0; ni3dten(k)=0.0d0
        enddo altitude_loop
        la=.NOT.lcheck
      END FUNCTION hugh_update_crys_2zero
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_rain_2zero(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: lrain2zero
        CHARACTER*24 :: sname='hugh_update_rain_2zero: '
        la     = .true.
        lrain2zero = .true.; lrain2zero = .false.
c Update qc3dten
        altitude_loop: do k=1,mx0
          qr3dten(k)=0.0d0; nr3dten(k)=0.0d0
        enddo altitude_loop
        la=.NOT.lcheck
      END FUNCTION hugh_update_rain_2zero
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_snow_2zero(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: lsnow2zero
        CHARACTER*24 :: sname='hugh_update_snow_2zero: '
        la     = .true.
        lsnow2zero = .true.; lsnow2zero = .false.
c Update qc3dten
        altitude_loop: do k=1,mx0
          qs3dten(k)=0.0d0; ns3dten(k)=0.0d0
        enddo altitude_loop
        la=.NOT.lcheck
      END FUNCTION hugh_update_snow_2zero
C
C==============================================================================
      LOGICAL FUNCTION hugh_update_all_2zero(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: mkx
      logical :: lall2zero
        CHARACTER*24 :: sname='hugh_update_all_2zero: '
        la     = .true.
        lall2zero = .true.!; lall2zero = .false.
c Set to zero:  qc3dten,qi3dten,qr3dten,qs3dten
c Set to zero:  nc3dten,ni3dten,nr3dten,ns3dten
c Set to zero:  tk3dten,qv3dten
c Set to zero:  all microphysical tendencies
c
        if(.not. lall2zero) return
        la = hugh_update_tkqv_2zero(dtmic,mx0)
        la = hugh_update_drop_2zero(dtmic,mx0)
        la = hugh_update_crys_2zero(dtmic,mx0)
        la = hugh_update_rain_2zero(dtmic,mx0)
        la = hugh_update_snow_2zero(dtmic,mx0)
        altitude_loop: do mkx=1,mx0
             nsubc(mkx)=0.0d0         ! loss of nc during evap
             nsubi(mkx)=0.0d0         ! loss of ni during sub.
             nsubs(mkx)=0.0d0         ! loss of ns during sub.
             nsubr(mkx)=0.0d0         ! loss of nr during evap
             ncondc(mkx)=0.0d0        ! cond/evap droplets number
             ncondi(mkx)=0.0d0        ! dep/sub cloud ice number
             ncondr(mkx)=0.0d0        ! cond/evap of rain number
             nconds(mkx)=0.0d0        ! dep/sub snow mass number

             mcondc(mkx)=0.0d0        ! cond/evap droplets mass
             mcondi(mkx)=0.0d0        ! dep/sub cloud ice mass
             mcondr(mkx)=0.0d0        ! cond/evap of rain mass
             mconds(mkx)=0.0d0        ! dep/sub snow mass
             nnuccc(mkx)=0.0d0        ! change n due to con. droplets freez 
             mnuccc(mkx)=0.0d0        ! change q due to con. droplets freez 
             nnucci(mkx)=0.0d0        ! change n due to imm. droplets freez 
             mnucci(mkx)=0.0d0        ! change q due to imm. droplets freez 
             npccn(mkx)=0.0d0         ! change n droplets activation
             mpccn(mkx)=0.0d0         ! change q droplets activation
             nnuccd(mkx)=0.0d0        ! change n freezing aerosol
             mnuccd(mkx)=0.0d0        ! change q freezing aerosol
             nnucmd(mkx)=0.0d0        ! change n cond freezing Meyers
             mnucmd(mkx)=0.0d0        ! change q cond freezing Meyers
             nnucmt(mkx)=0.0d0        ! change n cont freezing Meyers
             mnucmt(mkx)=0.0d0        ! change q cont freezing Meyers
             nnuccr(mkx)=0.0d0        ! change n due to con rain freez 
             mnuccr(mkx)=0.0d0        ! change q due to con rain freez 
             nnucir(mkx)=0.0d0        ! change n due to imm rain freez 
             mnucir(mkx)=0.0d0        ! change q due to imm rain freez 
c............................................................................
c water-water interactions
c............................................................................
             nprc(mkx)=0.0d0          ! change n autoconversion droplets
             mprc(mkx)=0.0d0          ! autoconversion droplets
             ncagg(mkx)=0.0d0         ! change n self-collection droplets
             mcagg(mkx)=0.0d0         ! self-collection droplets is ZERO
             npra(mkx)=0.0d0          ! change in n due to droplet acc by rain
             mpra(mkx)=0.0d0          ! accretion droplets by rain
             nragg(mkx)=0.0d0         ! change n self-collection rain
             mragg(mkx)=0.0d0         ! self-collection rain is ZERO
c ice-ice interactions
c crystal-crystal
             nprci(mkx)=0.0d0         ! change n autoconversion cloud ice
             mprci(mkx)=0.0d0         ! change q autoconversion cloud ice
             niagg(mkx)=0.0d0         ! change n self-collection cloud ice
             miagg(mkx)=0.0d0         ! self-collection cloud ice is ZERO
c crystal-snow
             nprai(mkx)=0.0d0         ! change n accretion cloud ice by snow
             mprai(mkx)=0.0d0         ! change q accretion cloud ice by snow
c snow-snow
             nsagg(mkx)=0.0d0         ! change n self-collection of snow
             msagg(mkx)=0.0d0         ! self-collection of snow is ZERO
c water-ice interactions
c drop-crystal
c drop-snow
             npsacws(mkx)=0.0d0       ! change n droplet accretion by snow
             mpsacws(mkx)=0.0d0       ! change q droplet accretion by snow
             nmults(mkx)=0.0d0        ! ice mult due to acc droplets by snow
             mmults(mkx)=0.0d0        ! change q due to ice mult droplets/snow
c rain-crystal
c rain-snow
             npracs(mkx)=0.0d0        ! change n collection rain by snow
             mpracs(mkx)=0.0d0        ! change q collection rain by snow
             nmultr(mkx)=0.0d0        ! ice mult due to acc rain by snow
             mmultr(mkx)=0.0d0        ! change q due to ice mult rain/snow
c melting  snow to rain
             nsmltr(mkx)=0.0d0        ! change n melting snow to rain
             msmltr(mkx)=0.0d0        ! change q melting snow to rain
             nsmlts(mkx)=0.0d0        ! change n melting snow
             msmlts(mkx)=0.0d0        ! chnage q melting snow evaporating
c melting  crys to drop
             nmltii(mkx)=0.0d0        ! change n melting cloud ice evaporating
             mmltii(mkx)=0.0d0        ! change q melting cloud ice evaporating
             nmltic(mkx)=0.0d0        ! change n melting cloud ice to droplets
             mmltic(mkx)=0.0d0        ! change q melting cloud ice to droplets
c homogenous freezing
             nhfrc(mkx)=0.0d0         ! change n homog freezing droplets to ice
             mhfrc(mkx)=0.0d0         ! change q homog freezing droplets to ice
             nhfrr(mkx)=0.0d0         ! change n homog freezing rain to graupel
             mhfrr(mkx)=0.0d0         ! change q homog freezing rain to graupel
             nsedimc(mkx)=0.0d0       ! change n sedimentation drop
             msedimc(mkx)=0.0d0       ! change q sedimentation drop
             nsedimi(mkx)=0.0d0       ! change n sedimentation ice
             msedimi(mkx)=0.0d0       ! change q sedimentation ice
             nsedimr(mkx)=0.0d0       ! change n sedimentation rain
             msedimr(mkx)=0.0d0       ! change q sedimentation rain
             nsedims(mkx)=0.0d0       ! change n sedimentation snow
             msedims(mkx)=0.0d0       ! change q sedimentation snow
             nvtermc(mkx)=0.0d0       ! n-weighted sedim vel drop
             mvtermc(mkx)=0.0d0       ! q-weighted sedim vel drop
             nvtermi(mkx)=0.0d0       ! n-weighted sedim vel ice
             mvtermi(mkx)=0.0d0       ! q-weighted sedim vel ice
             nvtermr(mkx)=0.0d0       ! n-weighted sedim vel rain
             mvtermr(mkx)=0.0d0       ! q-weighted sedim vel rain
             nvterms(mkx)=0.0d0       ! n-weighted sedim vel snow
             mvterms(mkx)=0.0d0       ! q-weighted sedim vel snow
        enddo altitude_loop
        la=.NOT.lcheck
      END FUNCTION hugh_update_all_2zero
C
C========================================================================
C Include file: blk_add_surabi_interfaces.f
C========================================================================
C==============================================================================
      LOGICAL FUNCTION surabi_drop_activation(dt0,k0) RESULT (la)
      IMPLICIT NONE
        REAL*8,             INTENT (in)       :: dt0
        INTEGER,            INTENT (in)       :: k0
c Local
        real*8       :: dtmic
        integer      :: mx0!,k=1
        LOGICAL      :: lnucl!,make_distribution
        character*24 :: sname='surabi_drop_activation: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'drop',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=surabi_make_drop_activation(dtmic,mx0)
        la=.NOT.lcheck
      END FUNCTION surabi_drop_activation
C
C========================================================================
C Include file: blk_add_surabi_drop_activation.f
C========================================================================
C==============================================================================
      LOGICAL FUNCTION surabi_make_drop_activation(dt0,mx0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: mx0
c droplet activation/freezing aerosol

       real*8  
     1       ct,      ! droplet activation parameter
     1       temp1,   ! dummy temperature
     1       sat1,    ! dummy saturation
     1       sigvl,   ! surface tension liq/vapor
     1       kel,     ! Kelvin parameter
     1       kc2      ! total ice nucleation rate

       real*8 
     1       cry,kry   ! aerosol activation parameters
c working parameters for aerosol activation
        real*8 
     1   aact,gamm,gg,psi,eta1,eta2,sm1,sm2,smax,uu1,uu2,
     1   alpha
c dummy variables
	real*8 
     1       dum,dum1,dum2,dum3,dum4,dumt,qdumt,dumqv,
     1       dumqss,dumqsi,dums,dume
c my Local
        real*8       :: p0,w0,v0,r0,ra,dtmic
        integer :: k=1
        LOGICAL :: lnucl,make_distribution
        character*29 :: sname='surabi_make_drop_activation: '
        la     = .true.
        dtmic=dt0!; dtmic=1.0d0
        ldummy=make_coefs(mx0)

c update variables
       altitude_loop: do k=1,mx0
       npccn(k)=0.0d0; mpccn(k)=0.0d0
	dumt = tk3d(k)
	dumqv = qv3d(k)
	dumqss = 0.622*polysvp(dumt,0)/
     1          (pp3d(k)-polysvp(dumt,0))
        dume  = pp3d(k)*dumqv / (0.622d0 + dumqv)
	lnucl=.false.
	if (dumqv/dumqss.ge.0.95) then
           lnucl=.true.
c droplet activation from Abdul-Razzak and Ghan (2000)
c convert heating rate to effective vertical velocity
c effective vertical velocity (m/s)
           dum = ww3d(k) !-cp/g*rttd(k)
c add sub-grid vertical velocity
	   dum = dum+wvar(k)
c assume minimum eff. sub-grid velocity 0.20 m/s
           dum = max(dum,0.40d0)                 
           sigvl = 0.0761-1.55e-4*(tk3d(k)-tf)
           aact = 2.*mw/(rhow*rr)*sigvl/tk3d(k)
           alpha = g*mw*xxlv(k)/(cp*rr*tk3d(k)**2)-
     1       g*ma/(rr*tk3d(k))
           gamm = rr*tk3d(k)/(evs(k)*mw)+mw*xxlv(k)**2/
     1      (cp*pp3d(k)*ma*tk3d(k))
           gg = 1./(rhow*rr*tk3d(k)/(evs(k)*Dv(k)*mw)+
     1       xxlv(k)*rhow/(kap(k)*tk3d(k))*(xxlv(k)*mw/
     1         (tk3d(k)*rr)-1.))
           psi = 2./3.*(alpha*dum/gg)**0.5*aact
           eta1 = (alpha*dum/gg)**1.5/
     1          (2.*pi*rhow*gamm*nanew1)
           eta2 = (alpha*dum/gg)**1.5/
     1          (2.*pi*rhow*gamm*nanew2)
           sm1 = 2./bact**0.5*(aact/(3.*rm1))**1.5
           sm2 = 2./bact**0.5*(aact/(3.*rm2))**1.5
           dum1 = 1./sm1**2*(f11*(psi/eta1)**1.5+
     1            f21*(sm1**2/(eta1+3.*psi))**0.75)
           dum2 = 1./sm2**2*(f12*(psi/eta2)**1.5+
     1            f22*(sm2**2/(eta2+3.*psi))**0.75)
           smax = 1./(dum1+dum2)**0.5
           uu1 = 2.*log(sm1/smax)/(4.242*log(sig1))
           uu2 = 2.*log(sm2/smax)/(4.242*log(sig2))
           dum1 = nanew1/2.*(1.-derf1(uu1))
           dum2 = nanew2/2.*(1.-derf1(uu2))  
           dum3 = (dum1+dum2)/rho(k)  !convert to kg-1
c make sure this value is not greater than total number of aerosol
	      dum3 = min((nanew1+nanew2)/rho(k),dum3)
	      dum4 = (dum3-nc3d(k)) ! / dtmic
	      dum4 = max(0.d0,dum4)
              dum4=MERGE(dum4,0.0d0, (dum4 .gt. 0)) /dtmic
c              dum4=MERGE(dum3-nc3d(k),0.0d0,(dum3-nc3d(k)) .gt. 0)
	      npccn(k) = dum4/2.
Cigs all new droplets have the same mass: 
              mpccn(k) = mw0*npccn(k)
        endif
      enddo altitude_loop
c  
      la=.NOT.lcheck
c
      END FUNCTION surabi_make_drop_activation
C
C==============================================================================
      LOGICAL FUNCTION hugh_template(dt0,k0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: k0
c Local
        real*8       :: dtmic
        integer      :: mx0,k=1
        character*16 :: sname='hugh_template: '
        la     = .true.
        dtmic=dt0; dtmic=1.0d0
        mx0=k0
c Make hydrometeors distributions
        ldummy=make_distribution(dtmic,ci,di,pi,rhow,cs,ds,dr
     *  ,tk3d,pp3d,qsmall,'all',mx0)
c Make time & altitude dependent coefficients
        ldummy=make_coefs(mx0)
          lcold =  tk3d.le.tf
          lwarm =  .NOT. lcold
          lqc3d = qc3d.ge.qsmall
          lqi3d = qi3d.ge.qsmall
          lqs3d = qs3d.ge.qsmall
          lqr3d = qr3d.ge.qsmall
        ldummy=hugh_make_template(dtmic,mx0)
        la=.NOT.lcheck
      END FUNCTION hugh_template
C
C==============================================================================
      LOGICAL FUNCTION hugh_make_template(dtmic,mx0) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      logical :: lfreezC,lfreezI
        CHARACTER*26 :: sname='hugh_make_template: '
        la     = .true.
        lfreezC = .true.; lfreezC = .false.
        lfreezI = .true.; lfreezI = .false.

c Process description
        altitude_loop: do k=1,mx0
          nnuccr(k)=0.0d0; mnuccr(k)=0.0d0
          nnucir(k)=0.0d0; mnucir(k)=0.0d0
	  if_xxxx_exists: if (lqc3d(k) .and. tk3d(k).lt.269.15) then
            if_freezC: if(lfreezC) then
            endif if_freezC
            if_freezI: if(lfreezI) then
            endif if_freezI
	  endif if_xxxx_exists
c
        enddo altitude_loop
        la=.NOT.lcheck
      END FUNCTION hugh_make_template
C
C==============================================================================
C Include file: blk_add_allocate.f
C========================================================================
      LOGICAL FUNCTION allocate_arrays_bulk2m_driver(mkx) RESULT (la)
      IMPLICIT NONE
        INTEGER,            INTENT (in)       :: mkx
        CHARACTER*31 :: sname='allocate_arrays_bulk2m_driver: '
        la     = .true.
        print *,sname,'mkx     = ',mkx
	ALLOCATE(
     1       tk3dten(mkx), 
     1       qv3dten(mkx))
	ALLOCATE(
     1       tk3d(mkx),    
     1       qv3d(mkx),   
     1       pp3d(mkx),   
     1       pdel(mkx),   
     1       sw3d(mkx),  
     1       si3d(mkx),  
c    1       rttd(mkx),  
     1       ww3d(mkx),    
     1       wvar(mkx))   
C
	ALLOCATE(
     1       qc3d(mkx),
     1       qi3d(mkx),
     1       qs3d(mkx),
     1       qr3d(mkx), 
     1       nc3d(mkx), 
     1       ni3d(mkx), 
     1       ns3d(mkx), 
     1       nr3d(mkx)) 
C
	ALLOCATE(
     1       na4d(mkx,nmodes))
	ALLOCATE(
     1       qc3dten(mkx),
     1       qi3dten(mkx),
     1       qs3dten(mkx),
     1       qr3dten(mkx),
     1       nc3dten(mkx),
     1       ni3dten(mkx),
     1       ns3dten(mkx),
     1       nr3dten(mkx))
C
	ALLOCATE(
     1       effc(mkx),  
     1       effi(mkx),  
     1       effs(mkx),  
     1       effr(mkx))  
C
	ALLOCATE(
     1       lamc(mkx),
     1       lami(mkx),
     1       lams(mkx),
     1       lamr(mkx),
     1       cdist1(mkx),
     1       n0i(mkx),
     1       n0s(mkx),
     1       n0rr(mkx),
     1       pgam(mkx))
c
c Allocate 1D microphysics arrays
        la=allocate_arrays_bulk2m_1d(mkx)
        la=.NOT.lcheck
      END FUNCTION allocate_arrays_bulk2m_driver
C
C==============================================================================
      LOGICAL FUNCTION allocate_arrays_bulk2m_1d(mkx) RESULT (la)
      IMPLICIT NONE
        INTEGER,            INTENT (in)       :: mkx
        CHARACTER*27 :: sname='allocate_arrays_bulk2m_1d: '
        la     = .true.
        print *,sname,'mkx     = ',mkx
c microphysical processes
	ALLOCATE(
     1       nsubc(mkx),         ! loss of nc during evap
     1       nsubi(mkx),         ! loss of ni during sub.
     1       nsubs(mkx),         ! loss of ns during sub.
     1       nsubr(mkx),         ! loss of nr during evap
c
c     1       npcc(mkx),          ! cond/evap droplets
c     1       mpcc(mkx),          ! cond/evap droplets
c
c     1       prd(mkx),           ! dep/sub cloud ice
c     1       pre(mkx),           ! cond/evap of rain
c     1       prds(mkx),          ! dep/sub snow
c
     1       ncondc(mkx),        ! cond/evap droplets number
     1       ncondi(mkx),        ! dep/sub cloud ice number
     1       ncondr(mkx),        ! cond/evap of rain number
     1       nconds(mkx),        ! dep/sub snow mass number

     1       mcondc(mkx),        ! cond/evap droplets mass
     1       mcondi(mkx),        ! dep/sub cloud ice mass
     1       mcondr(mkx),        ! cond/evap of rain mass
     1       mconds(mkx),        ! dep/sub snow mass
     1       nnuccc(mkx),        ! change n due to con. droplets freez 
     1       mnuccc(mkx),        ! change q due to con. droplets freez 
     1       nnucci(mkx),        ! change n due to imm. droplets freez 
     1       mnucci(mkx),        ! change q due to imm. droplets freez 
     1       npccn(mkx),         ! change n droplets activation
     1       mpccn(mkx),         ! change q droplets activation
     1       nnuccd(mkx),        ! change n freezing aerosol (prim ice nuc)
     1       mnuccd(mkx),        ! change q freezing aerosol (prim ice nuc) 
     1       nnucmd(mkx),        ! change n cond freezing Meyers (prim ice nuc)
     1       mnucmd(mkx),        ! change q cond freezing Meyers (prim ice nuc)
     1       nnucmt(mkx),        ! change n cont freezing Meyers (prim ice nuc)
     1       mnucmt(mkx),        ! change q cont freezing Meyers (prim ice nuc)
     1       nnuccr(mkx),        ! change n due to con rain freez 
     1       mnuccr(mkx),        ! change q due to con rain freez 
     1       nnucir(mkx),        ! change n due to imm rain freez 
     1       mnucir(mkx),        ! change q due to imm rain freez 
     1       ncactv(mkx),        !activated fraction for all modes
c............................................................................
c water-water interactions
c............................................................................
     1       nprc(mkx),          ! change n autoconversion droplets
     1       mprc(mkx),          ! autoconversion droplets
     1       ncagg(mkx),         ! change n self-collection droplets
     1       mcagg(mkx),         ! self-collection droplets is ZERO
     1       npra(mkx),          ! change in n due to droplet acc by rain
     1       mpra(mkx),          ! accretion droplets by rain
     1       nragg(mkx),         ! change n self-collection rain
     1       mragg(mkx),         ! self-collection rain is ZERO
c ice-ice interactions
c crystal-crystal
     1       nprci(mkx),         ! change n autoconversion cloud ice
     1       mprci(mkx),         ! change q autoconversion cloud ice
     1       niagg(mkx),         ! change n self-collection cloud ice
     1       miagg(mkx),         ! self-collection cloud ice is ZERO
c crystal-snow
     1       nprai(mkx),         ! change n accretion cloud ice by snow
     1       mprai(mkx),         ! change q accretion cloud ice by snow
c snow-snow
     1       nsagg(mkx),         ! change n self-collection of snow
     1       msagg(mkx),         ! self-collection of snow is ZERO
c water-ice interactions
c drop-crystal
c drop-snow
     1       npsacws(mkx),       ! change n droplet accretion by snow
     1       mpsacws(mkx),       ! change q droplet accretion by snow
     1       nmults(mkx),        ! ice mult due to acc droplets by snow
     1       mmults(mkx),        ! change q due to ice mult droplets/snow
c rain-crystal
c rain-snow
     1       npracs(mkx),        ! change n collection rain by snow
     1       mpracs(mkx),        ! change q collection rain by snow
     1       nmultr(mkx),        ! ice mult due to acc rain by snow
     1       mmultr(mkx),        ! change q due to ice mult rain/snow
c melting  snow to rain
     1       nsmltr(mkx),        ! change n melting snow to rain
     1       msmltr(mkx),        ! change q melting snow to rain
     1       nsmlts(mkx),        ! change n melting snow
     1       msmlts(mkx),        ! chnage q melting snow evaporating
c melting  crys to drop
     1       nmltii(mkx),        ! change n melting cloud ice evaporating
     1       mmltii(mkx),        ! change q melting cloud ice evaporating
     1       nmltic(mkx),        ! change n melting cloud ice to droplets
     1       mmltic(mkx),        ! change q melting cloud ice to droplets
c homogenous freezing
     1       nhfrc(mkx),         ! change n homog freezing droplets to ice
     1       mhfrc(mkx),         ! change q homog freezing droplets to ice
     1       nhfrr(mkx),         ! change n homog freezing rain to graupel
     1       mhfrr(mkx))         ! change q homog freezing rain to graupel
c
	ALLOCATE(
     1       nsedimc(mkx),       ! change n sedimentation drop
     1       msedimc(mkx),       ! change q sedimentation drop
     1       nsedimi(mkx),       ! change n sedimentation ice
     1       msedimi(mkx),       ! change q sedimentation ice
     1       nsedimr(mkx),       ! change n sedimentation rain
     1       msedimr(mkx),       ! change q sedimentation rain
     1       nsedims(mkx),       ! change n sedimentation snow
     1       msedims(mkx),       ! change q sedimentation snow
     1       nvtermc(mkx),       ! n-weighted sedim vel drop
     1       mvtermc(mkx),       ! q-weighted sedim vel drop
     1       nvtermi(mkx),       ! n-weighted sedim vel ice
     1       mvtermi(mkx),       ! q-weighted sedim vel ice
     1       nvtermr(mkx),       ! n-weighted sedim vel rain
     1       mvtermr(mkx),       ! q-weighted sedim vel rain
     1       nvterms(mkx),       ! n-weighted sedim vel snow
     1       mvterms(mkx))       ! q-weighted sedim vel snow
	ALLOCATE(
     1       wnuc(mkx),          ! vertical velocity for primary crys nucleation
     1       anuc(mkx),          ! "a" coefficient for primary crys nucleation 
     1       bnuc(mkx))          ! "b" coefficient for primary crys nucleation 
c time-varying atmospheric parameters
	ALLOCATE(
     1       kap(mkx),           ! thermal conductivity of air
     1       evs(mkx),           ! saturation vapor pressure
     1       eis(mkx),           ! ice saturation vapor pressure
     1       qvs(mkx),           ! saturation mixing ratio
     1       qvi(mkx),           ! ice saturation mixing ratio
     1       qvqvs(mkx),         ! watre satration ratio
     1       qvqvsi(mkx),        ! ice   saturaion ratio
     1       Dv(mkx),            ! diffusivity of water vapor in air
     1       Dap(mkx),           ! diffusivity of aerosol
     1       xxls(mkx),          ! latent heat of sublimation
     1       xxlv(mkx),          ! latent heat of vaporization
     1       mu(mkx),            ! viscocity of air
     1       sc(mkx),            ! Schmidt number
     1       xlf(mkx),           ! latent heat of freezing
     1       rho(mkx),           ! air density
     1       ab(mkx),            ! correction to cond rate due to latent heating
     1       abi(mkx),           ! correction to depo rate due to latent heating
     1       fdum(mkx),          ! mean free path
     1       eii(mkx))           ! collection efficiency
c
	ALLOCATE(
     1       epsc(mkx),          ! 1/phase rel. time (see M2005), drop
     1       epsi(mkx),          ! 1/phase rel. time (see M2005), crys
     1       epss(mkx),          ! 1/phase rel. time (see M2005), snow
     1       epsr(mkx))          ! 1/phase rel. time (see M2005), rain
C
C fall speed working variables (defined in code)
C
 	ALLOCATE(
     1       fnc(mkx),dumfnc(mkx),faloutnc(mkx),
     1       fmc(mkx),dumfmc(mkx),faloutmc(mkx),
     1       fni(mkx),dumfni(mkx),faloutni(mkx),
     1       fmi(mkx),dumfmi(mkx),faloutmi(mkx),
     1       fnr(mkx),dumfnr(mkx),faloutnr(mkx),
     1       fmr(mkx),dumfmr(mkx),faloutmr(mkx),
     1       fns(mkx),dumfns(mkx),faloutns(mkx),
     1       fms(mkx),dumfms(mkx),faloutms(mkx))
c
c     1       dumi(mkx),dumr(mkx),dumfni(mkx),
c     1       fr(mkx),
c     1       fi(mkx),fni(mkx),faloutr(mkx),falouti(mkx),
c     1       faloutni(mkx),
c     1       dumqs(mkx),dumfns(mkx),fs(mkx),fns(mkx),
c     1       falouts(mkx),faloutns(mkx),
c     1       dumc(mkx),dumfnc(mkx),fc(mkx),faloutc(mkx),
c     1       faloutnc(mkx),fnc(mkx),
c     1       dumfnr(mkx),faloutnr(mkx),fnr(mkx))
C
C fall-speed parameter 'a' with air density correction
C
 	ALLOCATE(ain(mkx),arn(mkx),asn(mkx),acn(mkx))
 	ALLOCATE(
     1       lwarm(mkx),lcold(mkx),
     1       lqc3d(mkx),lqi3d(mkx),lqr3d(mkx),lqs3d(mkx),
     1       lnc3d(mkx),lni3d(mkx),lnr3d(mkx),lns3d(mkx))
        la=.NOT.lcheck
      END FUNCTION allocate_arrays_bulk2m_1d
C
C==============================================================================
      LOGICAL FUNCTION deallocate_arrays_bulk2m_driver() RESULT (la)
      IMPLICIT NONE
        CHARACTER*33 :: sname='deallocate_arrays_bulk2m_driver: '
        la     = .true.
	DEALLOCATE(
     1       qc3dten,
     1       qi3dten,
     1       qs3dten,
     1       qr3dten,
     1       nc3dten,
     1       ni3dten,
     1       ns3dten,
     1       nr3dten)
	DEALLOCATE(
     1       qc3d,
     1       qi3d,
     1       qs3d,
     1       qr3d,
     1       nc3d,
     1       ni3d,
     1       ns3d,
     1       na4d,
     1       nr3d) 
	DEALLOCATE(
     1       tk3dten,
     1       qv3dten,
     1       tk3d,
     1       qv3d,
     1       pp3d,
     1       pdel,
     1       sw3d, 
     1       si3d, 
c    1       rttd,  
     1       ww3d, 
     1       wvar)   
	DEALLOCATE(
     1       effc,  
     1       effi, 
     1       effs,  
     1       effr)  
C
C input/output parameters           ! description (units)
	DEALLOCATE(
     1       lamc,
     1       lami,
     1       lams,
     1       lamr,
     1       cdist1,
     1       n0i,
     1       n0s,
     1       n0rr,
     1       pgam)
        la=deallocate_arrays_bulk2m_1d()
        la=.NOT.lcheck
      END FUNCTION deallocate_arrays_bulk2m_driver
C
C==============================================================================
      LOGICAL FUNCTION deallocate_arrays_bulk2m_1d() RESULT (la)
      IMPLICIT NONE
        CHARACTER*29 :: sname='deallocate_arrays_bulk2m_1d: '
        la     = .true.
        print *,sname
c microphysical processes
	DEALLOCATE(
     1       nsubc,
     1       nsubi,
     1       nsubs,
     1       nsubr,
     1       niagg,
c
c     1       npcc,
c     1       mpcc,
c
c     1       prd,
c     1       pre,
c     1       prds,
     1       ncondc,               ! cond/evap droplets number
     1       ncondi,             ! dep/sub cloud ice number
     1       ncondr,             ! cond/evap of rain number
     1       nconds,             ! dep/sub snow mass number
     1       mcondc,               ! cond/evap droplets mass
     1       mcondi,             ! dep/sub cloud ice mass
     1       mcondr,             ! cond/evap of rain mass
     1       mconds,             ! dep/sub snow mass
     1       nnuccc,
     1       mnuccc,
     1       nnucci,
     1       mnucci,
     1       nnuccd,
     1       mnuccd,
     1       nnucmd,
     1       mnucmd,
     1       nnucmt,
     1       mnucmt,
     1       nnuccr,
     1       mnuccr,
     1       nnucir,
     1       mnucir,
     1       nprc,
     1       mprc,
     1       ncagg,
     1       mcagg,
     1       npra,
     1       mpra,
     1       nragg,
     1       mragg,
     1       npccn,
     1       mpccn,
     1       nsagg,
     1       mprai,
     1       mprci,
     1       mpsacws,
     1       npsacws,
     1       nprci,
     1       nprai,
     1       nmults,
     1       nmultr,
     1       mmults,
     1       mmultr,
     1       mpracs,
     1       npracs,
     1       nsmltr,
     1       msmltr,
     1       nsmlts,
     1       msmlts,
     1       nmltii,
     1       mmltii,
     1       nmltic,
     1       mmltic,
     1       nhfrc,
     1       mhfrc,
     1       nhfrr,
     1       mhfrr,
     1       nsedimc,            ! change n sedimentation drop
     1       msedimc,            ! change q sedimentation drop
     1       nsedimi,            ! change n sedimentation ice
     1       msedimi,            ! change q sedimentation ice
     1       nsedimr,            ! change n sedimentation rain
     1       msedimr,            ! change q sedimentation rain
     1       nsedims,            ! change n sedimentation snow
     1       msedims,            ! change q sedimentation snow
     1       nvtermc,            ! n-weighted sedim vel drop
     1       mvtermc,            ! q-weighted sedim vel drop
     1       nvtermi,            ! n-weighted sedim vel ice
     1       mvtermi,            ! q-weighted sedim vel ice
     1       nvtermr,            ! n-weighted sedim vel rain
     1       mvtermr,            ! q-weighted sedim vel rain
     1       nvterms,            ! n-weighted sedim vel snow
     1       mvterms,            ! q-weighted sedim vel snow
     1       wnuc,
     1       anuc,
     1       ncactv,
     1       bnuc)

c time-varying atmospheric parameters
	DEALLOCATE(
     1       kap,
     1       evs,
     1       eis,
     1       qvs,
     1       qvi,
     1       qvqvs,
     1       qvqvsi,
     1       Dv,
     1       Dap,
     1       xxls,
     1       xxlv,
     1       mu,
     1       sc,
     1       xlf,
     1       rho,
     1       ab,
     1       abi,
     1       fdum,
     1       eii)
	DEALLOCATE(
     1       epsc,
     1       epsi,
     1       epss,
     1       epsr)
C
C fall speed working variables (defined in code)
C
c 	DEALLOCATE(
c     1       dumi,dumr,dumfni,
c     1       fr,
c     1       fi,fni,faloutr,falouti,
c     1       faloutni,
c     1       dumqs,dumfns,fs,fns,
c     1       falouts,faloutns,
c     1       dumc,dumfnc,fc,faloutc,
c     1       faloutnc,fnc,
c     1       dumfnr,faloutnr,fnr)

 	DEALLOCATE(
     1       fnc,dumfnc,faloutnc,
     1       fmc,dumfmc,faloutmc,
     1       fni,dumfni,faloutni,
     1       fmi,dumfmi,faloutmi,
     1       fnr,dumfnr,faloutnr,
     1       fmr,dumfmr,faloutmr,
     1       fns,dumfns,faloutns,
     1       fms,dumfms,faloutms)

C
C fall-speed parameter 'a' with air density correction
 	DEALLOCATE(ain,arn,asn,acn)
 	DEALLOCATE(
     1       lwarm,lcold,
     1       lqc3d,lqi3d,lqr3d,lqs3d,
     1       lnc3d,lni3d,lnr3d,lns3d)
        la=.NOT.lcheck
      END FUNCTION deallocate_arrays_bulk2m_1d
C
C==============================================================================
      LOGICAL FUNCTION cleanup_bulk2m_driver() RESULT (la)
      IMPLICIT NONE
        INTEGER                               :: i
        CHARACTER*23 :: sname='cleanup_bulk2m_driver: '
c        ldummy=cleanup_bulk2m()
        ldummy = deallocate_arrays_bulk2m_driver()
        la=.NOT.lcheck
      END FUNCTION cleanup_bulk2m_driver
C===================================================================
C
C========================================================================
C Include file: blk_add_hugh_functions.f
C========================================================================
c===================================================================
      function hugh_make_supi(p,t,q) RESULT (ra)
      implicit none
      integer,parameter            :: iu6=6
      real*8,intent(   in),dimension(:) :: p     ! pressure,      [dyn/cm**2]
      real*8,intent(   in),dimension(:) :: t     ! temperature,   [K]
      real*8,intent(   in),dimension(:) :: q     ! vapor content, [g/g]
c Output ra
      real*8,dimension(size(p,1))  :: ra
c Local
      integer      :: k,kl
      real*8,dimension(size(p,1))  :: sw,si
      real*8       :: dum,dumw,dumi,tk,qk,pk
      character*16 :: sname='hugh_make_supi: '
      kl=size(p,1)
      do k=1,kl
            qk = q(k)
            tk = t(k)
            pk = p(k)
            dum=qk*pk/(0.622d0+qk)
            dumw=polysvp(tk,0)
            dumi=polysvp(tk,1)
c            print *,sname,k,dumw,dumi,dum
            sw(k)=dum/dumw-1.d0        
            si(k)=dum/dumi-1.d0
            if(tk.gt. tf) si(k)=sw(k)
      enddo
      ra = si
      end function hugh_make_supi
c===================================================================
      function hugh_make_supw(p,t,q) RESULT (ra)
      implicit none
      integer,parameter            :: iu6=6
      real*8,intent(   in),dimension(:) :: p     ! pressure,      [dyn/cm**2]
      real*8,intent(   in),dimension(:) :: t     ! temperature,   [K]
      real*8,intent(   in),dimension(:) :: q     ! vapor content, [g/g]
c Output ra
      real*8,dimension(size(p,1))  :: ra
c Local
      integer      :: k,kl
      real*8,dimension(size(p,1))  :: sw,si
      real*8       :: dum,dumw,dumi,tk,qk,pk
      character*16 :: sname='hugh_make_supw: '
      kl=size(p,1)
      do k=1,kl
            qk = q(k)
            tk = t(k)
            pk = p(k)
            dum=qk*pk/(0.622d0+qk)
            dumw=polysvp(tk,0)
            dumi=polysvp(tk,1)
            sw(k)=dum/dumw-1.d0        
            si(k)=dum/dumi-1.d0
            if(tk.gt. tf) si(k)=sw(k)
      enddo
      ra = sw
      end function hugh_make_supw
c===================================================================
      subroutine brm_make_supsat1d(p,t,q,sw,si,kl)
      implicit none
      integer, intent(in)                :: kl
      integer,parameter            :: iu6=6
      real*8,intent(   in),dimension(kl) :: p    ! pressure,      [dyn/cm**2]
      real*8,intent(   in),dimension(kl) :: t    ! temperature,   [K]
      real*8,intent(   in),dimension(kl) :: q    ! vapor content, [g/g]
      real*8,intent(inout),dimension(kl) :: sw   ! water supersaturation
      real*8,intent(inout),dimension(kl) :: si   ! ice   supersaturation
c Local
      integer      :: k
      real*8       :: dum,dumw,dumi,tk,qk,pk
      character*19 :: sname='brm_make_supsat1d: '
      do k=1,kl
            qk = q(k)
            tk = t(k)
            pk = p(k)
            dum=qk*pk/(0.622d0+0.378d0*qk)
            dumw=polysvp(tk,0)
            dumi=polysvp(tk,1)
            sw(k)=dum/dumw-1.d0        
            si(k)=dum/dumi-1.d0
            if(tk.gt. tf) si(k)=sw(k)
      enddo
      return
      end subroutine brm_make_supsat1d
C==============================================================================
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! error function in single precision
c    Copyright(C) 1996 Takuya OOURA (email: ooura@mmm.t.u-tokyo.ac.jp).
c    You may use, copy, modify this code for any purpose and 
c    without fee. You may distribute this ORIGINAL package.
      real*8 function derf1(xx)
	implicit none
	real*8,intent(in) ::  xx
	real x
      double precision a(0 : 64), b(0 : 64)
	double precision w,t,y
	integer k,i
      data (a(i), i = 0, 12) / 
     &    0.00000000005958930743d0, -0.00000000113739022964d0, 
     &    0.00000001466005199839d0, -0.00000016350354461960d0, 
     &    0.00000164610044809620d0, -0.00001492559551950604d0, 
     &    0.00012055331122299265d0, -0.00085483269811296660d0, 
     &    0.00522397762482322257d0, -0.02686617064507733420d0, 
     &    0.11283791670954881569d0, -0.37612638903183748117d0, 
     &    1.12837916709551257377d0 / 
      data (a(i), i = 13, 25) / 
     &    0.00000000002372510631d0, -0.00000000045493253732d0, 
     &    0.00000000590362766598d0, -0.00000006642090827576d0, 
     &    0.00000067595634268133d0, -0.00000621188515924000d0, 
     &    0.00005103883009709690d0, -0.00037015410692956173d0, 
     &    0.00233307631218880978d0, -0.01254988477182192210d0, 
     &    0.05657061146827041994d0, -0.21379664776456006580d0, 
     &    0.84270079294971486929d0 / 
      data (a(i), i = 26, 38) / 
     &    0.00000000000949905026d0, -0.00000000018310229805d0, 
     &    0.00000000239463074000d0, -0.00000002721444369609d0, 
     &    0.00000028045522331686d0, -0.00000261830022482897d0, 
     &    0.00002195455056768781d0, -0.00016358986921372656d0, 
     &    0.00107052153564110318d0, -0.00608284718113590151d0, 
     &    0.02986978465246258244d0, -0.13055593046562267625d0, 
     &    0.67493323603965504676d0 / 
      data (a(i), i = 39, 51) / 
     &    0.00000000000382722073d0, -0.00000000007421598602d0, 
     &    0.00000000097930574080d0, -0.00000001126008898854d0, 
     &    0.00000011775134830784d0, -0.00000111992758382650d0, 
     &    0.00000962023443095201d0, -0.00007404402135070773d0, 
     &    0.00050689993654144881d0, -0.00307553051439272889d0, 
     &    0.01668977892553165586d0, -0.08548534594781312114d0, 
     &    0.56909076642393639985d0 / 
      data (a(i), i = 52, 64) / 
     &    0.00000000000155296588d0, -0.00000000003032205868d0, 
     &    0.00000000040424830707d0, -0.00000000471135111493d0, 
     &    0.00000005011915876293d0, -0.00000048722516178974d0, 
     &    0.00000430683284629395d0, -0.00003445026145385764d0, 
     &    0.00024879276133931664d0, -0.00162940941748079288d0, 
     &    0.00988786373932350462d0, -0.05962426839442303805d0, 
     &    0.49766113250947636708d0 / 
      data (b(i), i = 0, 12) / 
     &    -0.00000000029734388465d0, 0.00000000269776334046d0, 
     &    -0.00000000640788827665d0, -0.00000001667820132100d0, 
     &    -0.00000021854388148686d0, 0.00000266246030457984d0, 
     &    0.00001612722157047886d0, -0.00025616361025506629d0, 
     &    0.00015380842432375365d0, 0.00815533022524927908d0, 
     &    -0.01402283663896319337d0, -0.19746892495383021487d0, 
     &    0.71511720328842845913d0 / 
      data (b(i), i = 13, 25) / 
     &    -0.00000000001951073787d0, -0.00000000032302692214d0, 
     &    0.00000000522461866919d0, 0.00000000342940918551d0, 
     &    -0.00000035772874310272d0, 0.00000019999935792654d0, 
     &    0.00002687044575042908d0, -0.00011843240273775776d0, 
     &    -0.00080991728956032271d0, 0.00661062970502241174d0, 
     &    0.00909530922354827295d0, -0.20160072778491013140d0, 
     &    0.51169696718727644908d0 / 
      data (b(i), i = 26, 38) / 
     &    0.00000000003147682272d0, -0.00000000048465972408d0, 
     &    0.00000000063675740242d0, 0.00000003377623323271d0, 
     &    -0.00000015451139637086d0, -0.00000203340624738438d0, 
     &    0.00001947204525295057d0, 0.00002854147231653228d0, 
     &    -0.00101565063152200272d0, 0.00271187003520095655d0, 
     &    0.02328095035422810727d0, -0.16725021123116877197d0, 
     &    0.32490054966649436974d0 / 
      data (b(i), i = 39, 51) / 
     &    0.00000000002319363370d0, -0.00000000006303206648d0, 
     &    -0.00000000264888267434d0, 0.00000002050708040581d0, 
     &    0.00000011371857327578d0, -0.00000211211337219663d0, 
     &    0.00000368797328322935d0, 0.00009823686253424796d0, 
     &    -0.00065860243990455368d0, -0.00075285814895230877d0, 
     &    0.02585434424202960464d0, -0.11637092784486193258d0, 
     &    0.18267336775296612024d0 / 
      data (b(i), i = 52, 64) / 
     &    -0.00000000000367789363d0, 0.00000000020876046746d0, 
     &    -0.00000000193319027226d0, -0.00000000435953392472d0, 
     &    0.00000018006992266137d0, -0.00000078441223763969d0, 
     &    -0.00000675407647949153d0, 0.00008428418334440096d0, 
     &    -0.00017604388937031815d0, -0.00239729611435071610d0, 
     &    0.02064129023876022970d0, -0.06905562880005864105d0, 
     &    0.09084526782065478489d0 / 
c
      x=xx
c
      w = abs(x)
      if (w .lt. 2.2d0) then
          t = w * w
          k = int(t)
          t = t - k
          k = k * 13
          y = ((((((((((((a(k) * t + a(k + 1)) * t + 
     &        a(k + 2)) * t + a(k + 3)) * t + a(k + 4)) * t + 
     &        a(k + 5)) * t + a(k + 6)) * t + a(k + 7)) * t + 
     &        a(k + 8)) * t + a(k + 9)) * t + a(k + 10)) * t + 
     &        a(k + 11)) * t + a(k + 12)) * w
      else if (w .lt. 6.9d0) then
          k = int(w)
          t = w - k
          k = 13 * (k - 2)
          y = (((((((((((b(k) * t + b(k + 1)) * t + 
     &        b(k + 2)) * t + b(k + 3)) * t + b(k + 4)) * t + 
     &        b(k + 5)) * t + b(k + 6)) * t + b(k + 7)) * t + 
     &        b(k + 8)) * t + b(k + 9)) * t + b(k + 10)) * t + 
     &        b(k + 11)) * t + b(k + 12)
          y = y * y
          y = y * y
          y = y * y
          y = 1 - y * y
      else
          y = 1
      end if
      if (x .lt. 0) y = -y
      derf1 = y
      end function derf1
C
C==============================================================================
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      real*8 function polysvp (T,type)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Compute saturation vapor pressure by using 
c function from Goff and Gatch (1946)
c  Polysvp returned in units of pa.
c  T is input in units of K.
c  type refers to saturation with respect to liquid (0) or ice (1)
      implicit none
      real*8, intent(in) :: T
      integer,intent(in) :: type
      real*8 dum
c ice
      if (type.eq.1) then
         polysvp = 10.**(-9.09718*(tf/t-1.)-3.56654*
     1     log10(tf/t)+0.876793*(1.-t/tf)+
     1     log10(6.1071))*100.
      end if
c liquid
      if (type.eq.0) then
         polysvp = 10.**(-7.90298*(373.16/t-1.)+
     1        5.02808*log10(373.16/t)-
     1        1.3816e-7*(10**(11.344*(1.-t/373.16))-1.)+
     1        8.1328e-3*(10**(-3.49149*(373.16/t-1.))-1.)+
     1        log10(1013.246))*100.                      
         end if

      end function polysvp
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      REAL*8 FUNCTION GAMMA(X)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

CD    DOUBLE PRECISION FUNCTION DGAMMA(X)
C----------------------------------------------------------------------
C
C THIS ROUTINE CALCULATES THE GAMMA FUNCTION FOR A REAL ARGUMENT X.
C   COMPUTATION IS BASED ON AN ALGORITHM OUTLINED IN REFERENCE 1.
C   THE PROGRAM USES RATIONAL FUNCTIONS THAT APPROXIMATE THE GAMMA
C   FUNCTION TO AT LEAST 20 SIGNIFICANT DECIMAL DIGITS.  COEFFICIENTS
C   FOR THE APPROXIMATION OVER THE INTERVAL (1,2) ARE UNPUBLISHED.
C   THOSE FOR THE APPROXIMATION FOR X .GE. 12 ARE FROM REFERENCE 2.
C   THE ACCURACY ACHIEVED DEPENDS ON THE ARITHMETIC SYSTEM, THE
C   COMPILER, THE INTRINSIC FUNCTIONS, AND PROPER SELECTION OF THE
C   MACHINE-DEPENDENT CONSTANTS.
C*******************************************************************
C*******************************************************************
C EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
C BETA   - RADIX FOR THE FLOATING-POINT REPRESENTATION
C MAXEXP - THE SMALLEST POSITIVE POWER OF BETA THAT OVERFLOWS
C XBIG   - THE LARGEST ARGUMENT FOR WHICH GAMMA(X) IS REPRESENTABLE
C          IN THE MACHINE, I.E., THE SOLUTION TO THE EQUATION
C                  GAMMA(XBIG) = BETA**MAXEXP
C XINF   - THE LARGEST MACHINE REPRESENTABLE FLOATING-POINT NUMBER;
C          APPROXIMATELY BETA**MAXEXP
C EPS    - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
C          1.0+EPS .GT. 1.0
C XMININ - THE SMALLEST POSITIVE FLOATING-POINT NUMBER SUCH THAT
C          1/XMININ IS MACHINE REPRESENTABLE
C
C     APPROXIMATE VALUES FOR SOME IMPORTANT MACHINES ARE:
C
C                            BETA       MAXEXP        XBIG
C
C CRAY-1         (S.P.)        2         8191        966.961
C CYBER 180/855
C   UNDER NOS    (S.P.)        2         1070        177.803
C IEEE (IBM/XT,
C   SUN, ETC.)   (S.P.)        2          128        35.040
C IEEE (IBM/XT,
C   SUN, ETC.)   (D.P.)        2         1024        171.624
C IBM 3033       (D.P.)       16           63        57.574
C VAX D-FORMAT   (D.P.)        2          127        34.844
C VAX G-FORMAT   (D.P.)        2         1023        171.489
C
C                            XINF         EPS        XMININ
C
C CRAY-1         (S.P.)   5.45E+2465   7.11E-15    1.84E-2466
C CYBER 180/855
C   UNDER NOS    (S.P.)   1.26E+322    3.55E-15    3.14E-294
C IEEE (IBM/XT,
C   SUN, ETC.)   (S.P.)   3.40E+38     1.19E-7     1.18E-38
C IEEE (IBM/XT,
C   SUN, ETC.)   (D.P.)   1.79D+308    2.22D-16    2.23D-308
C IBM 3033       (D.P.)   7.23D+75     2.22D-16    1.39D-76
C VAX D-FORMAT   (D.P.)   1.70D+38     1.39D-17    5.88D-39
C VAX G-FORMAT   (D.P.)   8.98D+307    1.11D-16    1.12D-308
C
C*******************************************************************
C*******************************************************************
C ERROR RETURNS
C  THE PROGRAM RETURNS THE VALUE XINF FOR SINGULARITIES OR
C     WHEN OVERFLOW WOULD OCCUR.  THE COMPUTATION IS BELIEVED
C     TO BE FREE OF UNDERFLOW AND OVERFLOW.
C
C  INTRINSIC FUNCTIONS REQUIRED ARE:
C     INT, DBLE, EXP, LOG, REAL, SIN
C REFERENCES:  AN OVERVIEW OF SOFTWARE DEVELOPMENT FOR SPECIAL
C              FUNCTIONS   W. J. CODY, LECTURE NOTES IN MATHEMATICS,
C              506, NUMERICAL ANALYSIS DUNDEE, 1975, G. A. WATSON
C              (ED.), SPRINGER VERLAG, BERLIN, 1976.
C
C              COMPUTER APPROXIMATIONS, HART, ET. AL., WILEY AND
C              SONS, NEW YORK, 1968.
C  LATEST MODIFICATION: OCTOBER 12, 1989
C  AUTHORS: W. J. CODY AND L. STOLTZ
C           APPLIED MATHEMATICS DIVISION
C           ARGONNE NATIONAL LABORATORY
C           ARGONNE, IL 60439
C----------------------------------------------------------------------
      implicit none
      REAL*8,intent(in) :: X
      INTEGER I,N
      LOGICAL PARITY
      REAL
     1    C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE,
     2    TWO,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO
      DIMENSION C(7),P(8),Q(8)
C----------------------------------------------------------------------
C  MATHEMATICAL CONSTANTS
C----------------------------------------------------------------------
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0E0,0.5E0,12.0E0,2.0E0,0.0E0/,
     1     SQRTPI/0.9189385332046727417803297E0/,
     2     PI/3.1415926535897932384626434E0/
CD    DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,
CD   1     SQRTPI/0.9189385332046727417803297D0/,
CD   2     PI/3.1415926535897932384626434D0/
C----------------------------------------------------------------------
C  MACHINE DEPENDENT PARAMETERS
C----------------------------------------------------------------------
      DATA XBIG,XMININ,EPS/35.040E0,1.18E-38,1.19E-7/,
     1     XINF/3.4E38/
CD    DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,
CD   1     XINF/1.79D308/
C----------------------------------------------------------------------
C  NUMERATOR AND DENOMINATOR COEFFICIENTS FOR RATIONAL MINIMAX
C     APPROXIMATION OVER (1,2).
C----------------------------------------------------------------------
      DATA P/-1.71618513886549492533811E+0,2.47656508055759199108314E+1,
     1       -3.79804256470945635097577E+2,6.29331155312818442661052E+2,
     2       8.66966202790413211295064E+2,-3.14512729688483675254357E+4,
     3       -3.61444134186911729807069E+4,6.64561438202405440627855E+4/
      DATA Q/-3.08402300119738975254353E+1,3.15350626979604161529144E+2,
     1      -1.01515636749021914166146E+3,-3.10777167157231109440444E+3,
     2        2.25381184209801510330112E+4,4.75584627752788110767815E+3,
     3      -1.34659959864969306392456E+5,-1.15132259675553483497211E+5/
CD    DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,
CD   1       -3.79804256470945635097577D+2,6.29331155312818442661052D+2,
CD   2       8.66966202790413211295064D+2,-3.14512729688483675254357D+4,
CD   3       -3.61444134186911729807069D+4,6.64561438202405440627855D+4/
CD    DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,
CD   1      -1.01515636749021914166146D+3,-3.10777167157231109440444D+3,
CD   2        2.25381184209801510330112D+4,4.75584627752788110767815D+3,
CD   3      -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/
C----------------------------------------------------------------------
C  COEFFICIENTS FOR MINIMAX APPROXIMATION OVER (12, INF).
C----------------------------------------------------------------------
      DATA C/-1.910444077728E-03,8.4171387781295E-04,
     1     -5.952379913043012E-04,7.93650793500350248E-04,
     2     -2.777777777777681622553E-03,8.333333333333333331554247E-02,
     3      5.7083835261E-03/
CD    DATA C/-1.910444077728D-03,8.4171387781295D-04,
CD   1     -5.952379913043012D-04,7.93650793500350248D-04,
CD   2     -2.777777777777681622553D-03,8.333333333333333331554247D-02,
CD   3      5.7083835261D-03/
C----------------------------------------------------------------------
C  STATEMENT FUNCTIONS FOR CONVERSION BETWEEN INTEGER AND FLOAT
C----------------------------------------------------------------------
      CONV(I) = REAL(I)
CD    CONV(I) = DBLE(I)
      PARITY=.FALSE.
      FACT=ONE
      N=0
      Y=X
      IF(Y.LE.ZERO)THEN
C----------------------------------------------------------------------
C  ARGUMENT IS NEGATIVE
C----------------------------------------------------------------------
        Y=-X
        Y1=AINT(Y)
        RES=Y-Y1
        IF(RES.NE.ZERO)THEN
          IF(Y1.NE.AINT(Y1*HALF)*TWO)PARITY=.TRUE.
          FACT=-PI/SIN(PI*RES)
          Y=Y+ONE
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
C----------------------------------------------------------------------
C  ARGUMENT IS POSITIVE
C----------------------------------------------------------------------
      IF(Y.LT.EPS)THEN
C----------------------------------------------------------------------
C  ARGUMENT .LT. EPS
C----------------------------------------------------------------------
        IF(Y.GE.XMININ)THEN
          RES=ONE/Y
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ELSEIF(Y.LT.TWELVE)THEN
        Y1=Y
        IF(Y.LT.ONE)THEN
C----------------------------------------------------------------------
C  0.0 .LT. ARGUMENT .LT. 1.0
C----------------------------------------------------------------------
          Z=Y
          Y=Y+ONE
        ELSE
C----------------------------------------------------------------------
C  1.0 .LT. ARGUMENT .LT. 12.0, REDUCE ARGUMENT IF NECESSARY
C----------------------------------------------------------------------
          N=INT(Y)-1
          Y=Y-CONV(N)
          Z=Y-ONE
        ENDIF
C----------------------------------------------------------------------
C  EVALUATE APPROXIMATION FOR 1.0 .LT. ARGUMENT .LT. 2.0
C----------------------------------------------------------------------
        XNUM=ZERO
        XDEN=ONE
        DO 260 I=1,8
          XNUM=(XNUM+P(I))*Z
          XDEN=XDEN*Z+Q(I)
  260   CONTINUE
        RES=XNUM/XDEN+ONE
        IF(Y1.LT.Y)THEN
C----------------------------------------------------------------------
C  ADJUST RESULT FOR CASE  0.0 .LT. ARGUMENT .LT. 1.0
C----------------------------------------------------------------------
          RES=RES/Y1
        ELSEIF(Y1.GT.Y)THEN
C----------------------------------------------------------------------
C  ADJUST RESULT FOR CASE  2.0 .LT. ARGUMENT .LT. 12.0
C----------------------------------------------------------------------
          DO 290 I=1,N
            RES=RES*Y
            Y=Y+ONE
  290     CONTINUE
        ENDIF
      ELSE
C----------------------------------------------------------------------
C  EVALUATE FOR ARGUMENT .GE. 12.0,
C----------------------------------------------------------------------
        IF(Y.LE.XBIG)THEN
          YSQ=Y*Y
          SUM=C(7)
          DO 350 I=1,6
            SUM=SUM/YSQ+C(I)
  350     CONTINUE
          SUM=SUM/Y-Y+SQRTPI
          SUM=SUM+(Y-HALF)*LOG(Y)
          RES=EXP(SUM)
        ELSE
          RES=XINF
          GOTO 900
        ENDIF
      ENDIF
C----------------------------------------------------------------------
C  FINAL ADJUSTMENTS AND RETURN
C----------------------------------------------------------------------
      IF(PARITY)RES=-RES
      IF(FACT.NE.ONE)RES=FACT/RES
  900 GAMMA=RES
CD900 DGAMMA = RES
C ---------- LAST LINE OF GAMMA ----------
      END FUNCTION GAMMA
C
C===================================================================
C
C
C========================================================================
C Include file: blk_add_gult_drop_activation.f
C========================================================================
C==============================================================================
      LOGICAL FUNCTION gult_drop_activation(dt0,mx0,r1,r2) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0,r1,r2
        integer, intent(in)                   :: mx0
       altitude_loop: do k=1,mx0
c      do k=1,mx0
        npccn(k)=max(0.d0,((r2-r1)/dt0/rho(k)))   
c       if(r2.gt.1.d7)print*, "gult_drop",r2,r1,dt0,npccn(k) 
C** gives new - old CDNC/timestep = rate 
        mpccn(k)=mw0*npccn(k)

       enddo altitude_loop

        la=.NOT.lcheck
      END FUNCTION gult_drop_activation

C========================================================================
C Include file: blk_add_matr_drop_activation.f
C========================================================================
C==============================================================================
      LOGICAL FUNCTION matr_drop_activation(dt0,mx0,l1) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: mx0
        logical, intent(in) :: l1
       altitude_loop: do k=1,mx0

c       npccn(k)=max(0.d0, ( (ncactv(k)/dt0/rho(k)) )! change to 2*dt0 if needed
        npccn(k)=max(0.d0, ( (ncactv(k)/rho(k) - nc3d(k))/(dt0) ) )! change to 2*dt0 if needed
c       write(6,*)"MATRIX_HM",nc3d(k)*1d-6,ncactv(k)*1d-6,npccn(k)*1d-6
        mpccn(k)=mw0*npccn(k)

       enddo altitude_loop

        la=.NOT.lcheck
      END FUNCTION matr_drop_activation
C========================================================================
C Include file: blk_add_matr_drop_activation.f
C========================================================================
C==============================================================================
      LOGICAL FUNCTION toma_drop_activation(dt0,mx0,l1) RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dt0
        integer, intent(in)                   :: mx0
        logical, intent(in) :: l1
       altitude_loop: do k=1,mx0

        npccn(k)=max(0.d0, ( (ncactv(k)/rho(k) - nc3d(k))/(dt0) ) )! change to 2*dt0 if needed
        mpccn(k)=mw0*npccn(k)

       enddo altitude_loop

        la=.NOT.lcheck
      END FUNCTION toma_drop_activation
C
C========================================================================
C
      LOGICAL FUNCTION srbm_balance_tendencies(dtmic,mx0) 
     *  RESULT (la)
      IMPLICIT NONE
        real*8,             INTENT (in)       :: dtmic
        integer, intent(in)                   :: mx0
c Local
c counting variables
      integer :: k
      real*8  :: dum,ratio
      logical :: ldrop0
        CHARACTER*25 :: sname='srbm_balance_tendencies: '
        la = .true.
c Balance tendencies
        altitude_loop: do k=1,mx0
         if(lqc3d(k).and.lqi3d(k).and.lqr3d(k)
     1     .and.lqs3d(k)) cycle altitude_loop
         warm_or_cold_balance: if(lcold(k)) then
c.......................................................................
c positiveness of concentration
c this section is separated into two parts, if T < 0 C, T > 0 C
c due to different processes that act depending on freezing/above freezing
c if concentrations less than nsmall, then no depletion of water
c through microphysical processes, skip conservation
c positiveness of nc
	dum = (nprc(k)+npra(k)+nnuccc(k)+
     1         npsacws(k)+nmults(k))*dtmic
	if (dum.gt.nc3d(k).and.nc3d(k).ge.nsmall) then
          ratio      = nc3d(k)/dum
	  nprc(k)    = nprc(k)*ratio
	  npra(k)    = npra(k)*ratio
	  nnuccc(k)  = nnuccc(k)*ratio
	  npsacws(k) = npsacws(k)*ratio
	  nmults(k)  = nmults(k)*ratio
        end if
c positiveness of ni
	dum = (-ncondi(k)-nnuccc(k)+nprci(k)+
     1         nprai(k)-nmults(k)-nmultr(k)-nnuccd(k))*dtmic
	if (dum.gt.ni3d(k).and.ni3d(k).ge.nsmall) then
          ratio = ni3d(k)/dum
          ncondi(k)  = ncondi(k)*ratio	
          nnuccc(k)  = nnuccc(k)*ratio
          nprci(k)   = nprci(k)*ratio
          nprai(k)   = nprai(k)*ratio
          nmults(k)  = nmults(k)*ratio
          nmultr(k)  = nmultr(k)*ratio
          nnuccd(k)  = nnuccd(k)*ratio
        end if
	! this point always works
c positiveness of nr
	dum=((npracs(k)-ncondr(k))+(nmultr(k)-nprc(k))
     1    +(nnuccr(k)+nnucir(k)-npra(k)))*dtmic
	if (dum.gt.nr3d(k).and.nr3d(k).ge.nsmall) then
          ratio = nr3d(k)/dum
          ncondr(k)  = ncondr(k)*ratio
          nprc(k)    = nprc(k)*ratio
          npra(k)    = npra(k)*ratio
          npracs(k)  = npracs(k)*ratio
          nmultr(k)  = nmultr(k)*ratio
          nnuccr(k)  = nnuccr(k)*ratio
          nnucir(k)  = nnucir(k)*ratio

        end if
c positiveness of ns
	dum = (-nconds(k)-npsacws(k)-nprai(k)-nprci(k)-
     1          npracs(k)-nnuccr(k)-nnucir(k))*dtmic

	if (dum.gt.ns3d(k).and.ns3d(k).ge.nsmall) then
          ratio = ns3d(k)/dum
          nconds(k)  = nconds(k)*ratio
          npsacws(k) = npsacws(k)*ratio
          nprai(k)   = nprai(k)*ratio
          nprci(k)   = nprci(k)*ratio
          npracs(k)  = npracs(k)*ratio
          nnuccr(k)  = nnuccr(k)*ratio
          nnucir(k)  = nnucir(k)*ratio

       end if
         else ! warm cloud: temperature above freezing
c for cloud ice and snow, only processes operating at T > tf is
c melting/evaporating, which is already conserved during process
c calculation
c positiveness of nc
	dum = (nprc(k)+npra(k))*dtmic
	if (dum.gt.nc3d(k).and.nc3d(k).ge.nsmall) then
          ratio   = nc3d(k)/dum
          nprc(k) = nprc(k)*ratio
          npra(k) = npra(k)*ratio
        end if
c positiveness of nr
c all terms are positive, do not need conservation of rain
c positiveness of snow
        dum = (-nsmltr(k)-nsmlts(k))*dtmic
        if (dum.gt.ns3d(k).and.ns3d(k).ge.nsmall) then
          ratio     = ns3d(k)/dum
          nsmltr(k) = nsmltr(k)*ratio
          nsmlts(k) = nsmlts(k)*ratio
        endif
         endif warm_or_cold_balance
        enddo altitude_loop
c
        la=.NOT.lcheck
      END FUNCTION srbm_balance_tendencies

      END module mo_bulk2m_driver_gcm
