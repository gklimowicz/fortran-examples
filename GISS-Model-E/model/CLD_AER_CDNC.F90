#include "rundeck_opts.h"

!@sum module containing (CLD_AER_CDNC || BLK_2MOM) code blocks
!@+   formerly embedded in LSCOND which have been wrapped in
!@+   subroutines.  No further structure has been introduced.
!@+   The never-used capability to run BLK_2MOM without CLD_AER_CDNC
!@+   was removed during the code transfer.

module cld_aer_cdnc_mod
  use mo_bulk2m_driver_gcm, only: execute_bulk2m_driver
#if defined(TRACERS_AMP)
  use AERO_CONFIG, only: NMODES
#endif
  implicit none
  private
  public :: sntm,cld_aer_cdnc_block1,cld_aer_cdnc_block2
  integer, parameter :: SNTM=31

#if defined(TRACERS_AEROSOLS_Koch) || defined(TRACERS_AEROSOLS_SEASALT) || \
    defined(TRACERS_DUST) || defined(TRACERS_NITRATE) || \
    defined(TRACERS_HETCHEM) || defined(TRACERS_SOA) || \
    defined(TRACERS_AEROSOLS_OCEAN) || defined(TRACERS_AEROSOLS_VBS)
  public :: cld_aer_cdnc_block0
#endif

    integer,parameter         :: mkx=1   ! lm
    real*8,parameter          :: mw0 = 2.094395148947515E-15
    real*8,parameter          :: mi0 = 2.094395148947515E-15
    logical, parameter        :: lSCM=.false.
    logical, parameter        :: wSCM=.false.
    real(8), parameter        :: tiny = 1.0D-30
    real*8,dimension(mkx)     :: tk0,qk0,pk0,w0,v0,r0,rablk
    real*8,dimension(mkx)     :: tk0new,qk0new
    real*8,dimension(mkx)     :: ndrop,mdrop,ncrys,mcrys
    real*8,dimension(mkx)     :: nrain,mrain,mtau,rtau,ctau,dtau
    real*8                    :: qrArray(5),row,mr0,piby6
    real*8,dimension(mkx)     :: vrain
    real*8,dimension(mkx)     :: ndrop_old,mdrop_old
    real*8,dimension(mkx)     :: ndrop_new,mdrop_new
    real*8,dimension(mkx)     :: ndrop_blk,mdrop_blk
    real*8,dimension(mkx)     :: ndrop_res,mdrop_res
    real*8,dimension(mkx)     :: ncrys_old,mcrys_old
    real*8,dimension(mkx)     :: ncrys_new,mcrys_new
    real*8,dimension(mkx)     :: ncrys_blk,mcrys_blk
    real*8,dimension(mkx)     :: ncrys_res,mcrys_res
    real*8,dimension(mkx)     :: npccn,nprc,nnuccc,nnucci
    real*8,dimension(mkx)     :: mpccn,mprc,mnuccc,mnucci
    real*8,dimension(mkx)     :: nnuccd,nnucmd,nnucmt
    real*8,dimension(mkx)     :: mnuccd,mnucmd,mnucmt
    real*8,dimension(mkx)     :: nc_tnd,qc_tnd,ni_tnd,qi_tnd
    real*8,dimension(mkx)     :: nc_tot,qc_tot,ni_tot,qi_tot
    real*8                    :: DTB2M,QAUT_B2M
    real*8                    :: NEWCDNC,OLDCDNC
    real*8                    :: DCLD
    logical                   :: ldummy=.false.
    integer                   :: nm,iuo=801

contains

#if defined(TRACERS_AEROSOLS_Koch) || defined(TRACERS_AEROSOLS_SEASALT) || \
    defined(TRACERS_DUST) || defined(TRACERS_NITRATE) || \
    defined(TRACERS_HETCHEM) || defined(TRACERS_SOA) || \
    defined(TRACERS_AEROSOLS_OCEAN) || defined(TRACERS_AEROSOLS_VBS)

  subroutine cld_aer_cdnc_block0( &
       ntx,ntix, &
       lhx,fcld,vvel,dxypij, &
       qclx,qcix,cleara,cldsavl,pl,tl,ncll,sme,airm, &
       tm, &
       dsu, &
       oldcdn,newcdn &
       )
    use constant, only : lhe
    use TRACER_COM, only: NTM
    use OldTracer_mod, only: trname
    implicit none

    real*8 :: lhx,fcld,vvel,oldcdn,newcdn,dxypij
    real*8 DSS(SNTM),DSU(SNTM)
    real*8 :: qclx,qcix,cleara,cldsavl,pl,tl,ncll,sme,airm
    integer :: ntx
    integer, dimension(ntm) :: ntix
    real*8, dimension(ntm) :: tm

    integer :: n
    real*8 :: qcx,cdnl0,cdnl1,snd

!@auth Menon  - storing var for cloud droplet number

!@auth Menon  saving aerosols mass for CDNC prediction
      do N=1,SNTM
        DSS(N)=1.d-10
      end do
      do N=1,NTX
        select case (trname(ntix(n)))
#ifdef TRACERS_AEROSOLS_Koch
        case('SO4')
          DSS(1)=tm(n)     !n=4
        case('OCIA')
          DSS(4)=tm(n)     !n=12
        case('OCB')
          DSS(5)=tm(n)     !n=13
        case('BCIA')
          DSS(6)=tm(n)     !n=9
        case('BCB')
          DSS(7)=tm(n)     !n=10
        case('OCII')
          DSS(8)=tm(n)     !n=11
        case('BCII')
          DSS(9)=tm(n)     !n=8
#endif  /* TRACERS_AEROSOLS_Koch */
#ifdef TRACERS_AEROSOLS_SEASALT
        case('seasalt1')
          DSS(2)=tm(n)     !n=6
        case('seasalt2')
          DSS(3)=tm(n)     !n=7
#endif  /* TRACERS_AEROSOLS_SEASALT */
#ifdef TRACERS_DUST
        case('Clay')
          DSS(10)=tm(n)    !n=23
        case('Silt1')
          DSS(11)=tm(n)    !n=23
        case('Silt2')
          DSS(12)=tm(n)    !n=23
        case('Silt3')
          DSS(13)=tm(n)    !n=23
#endif  /* TRACERS_DUST */
#ifdef TRACERS_NITRATE
        case('NO3p')
          DSS(14)=tm(n)    !n=23
#endif
#ifdef TRACERS_HETCHEM
          !**** Here are dust particles coated with sulfate
        case('SO4_d1')
          DSS(15)=tm(n)    !n=20
        case('SO4_d2')
          DSS(16)=tm(n)    !n=21
        case('SO4_d3')
          DSS(17)=tm(n)    !n=22
#endif  /* TRACERS_HETCHEM */
#ifdef TRACERS_AEROSOLS_SOA
        case('isopp1a')
          DSS(18)=tm(n)
        case('isopp2a')
          DSS(19)=tm(n)
#ifdef TRACERS_TERP
        case('apinp1a')
          DSS(20)=tm(n)
        case('apinp2a')
          DSS(21)=tm(n)
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_OCEAN
        case('OCocean')
          DSS(22)=tm(n)
#endif  /* TRACERS_AEROSOLS_OCEAN */
#ifdef TRACERS_AEROSOLS_VBS
        case('vbsAm2')
          DSS(23)=tm(n)
        case('vbsAm1')
          DSS(24)=tm(n)
        case('vbsAz')
          DSS(25)=tm(n)
        case('vbsAp1')
          DSS(26)=tm(n)
        case('vbsAp2')
          DSS(27)=tm(n)
        case('vbsAp3')
          DSS(28)=tm(n)
        case('vbsAp4')
          DSS(29)=tm(n)
        case('vbsAp5')
          DSS(30)=tm(n)
        case('vbsAp6')
          DSS(31)=tm(n)
#endif  /* TRACERS_AEROSOLS_VBS */
        end select
      end do      !end of n loop for tracers

      if( LHX.eq.LHE )then
        QCX = QCLX
      else
        QCX = QCIX
      endif
      call GET_CDNC(LHX,AIRM,QCX,DXYPIJ, &
           FCLD,CLEARA,CLDSAVL,DSS,PL,TL, &
           NCLL,VVEL,SME,DSU,CDNL0,CDNL1)
      !     write(6,*)"Where is",DSU(L),l
      SNd=CDNL1
      !C** Pass old and new cloud droplet number
      NEWCDN=SNd
      OLDCDN=CDNL0
      !     if(SNd.gt.20.) write(6,*)"SM 11 CDNC",NEWCDN,OLDCDN,L

  end subroutine cld_aer_cdnc_block0

#endif

  subroutine cld_aer_cdnc_block1( &
       i_debug,j_debug, &
       oldcdn,newcdn, &
       dtsrc,vvel,lhx,fcld,dxypij,pearth, &
       prebar,tl,ql,pl,wturb,cldsavl,qclx,qcix,ncll,ncil,airm &
#if defined(TRACERS_AMP)
       ,nactc &
#endif
#ifdef TRACERS_TOMAS
       ,tm &
       ,cdnc_tomas &
#endif
       ,SNd &
       ,scdncw,scdnci &
       )
  use constant, only : lhe,rgas,teeny,mb2kg
  use tracer_com, only: ntm
  implicit none

  integer :: i_debug,j_debug
  real*8 :: oldcdn,newcdn,SNd
  real*8 :: dtsrc,vvel,lhx,fcld,scdncw,scdnci,dxypij,pearth
  real*8 :: prebar,tl,ql,pl,wturb,cldsavl,qclx,qcix,ncll,ncil,airm
#ifdef TRACERS_AMP
  real*8 :: nactc(nmodes)
#endif
#ifdef TRACERS_TOMAS
  real*8, dimension(ntm) :: tm
  real*8 :: cdnc_tomas
#endif

! local vars
  real*8 :: svlhxl,wmxice,rho,sndi
#ifdef TRACERS_AMP
  real*8 :: naero (mkx,nmodes)
#endif

#ifdef TRACERS_TOMAS
      REAL*8,dimension(mkx)     :: nactl
!Can
!Can Droplet parameterization quantities
!Can
      INTEGER :: NCCNMx,NCC,NSECi
      PARAMETER (NCCNMx=100, NCC=10)
      REAL*8 SULFI, BOXVL, TOTi,TOT_MI, &
           TPi(NCCNMx), MLi(NCCNMx), SLVL(NCC), CCON(NCC), &
           NACTEarth, NACTOcean, NACT, NACTBL, &
           SMAXEarth, SMAXOcean, SMAX, &
           REFFEarth, REFFOcean, REFF, REFFBL, REFFGISS, &
           CLDTAUEarth, CLDTAUOcean, CLDTAU, CLDTAUBL, CLDTAULIQ, &
           CLDTAUICE, TPARC,PPARC, &
           WPARCOcean, WPARCEarth, WPARC, &
           RHOSI,QLWC,EPSILON,AUTO(6),DIFFLWMR,DIFFEPS

      INTEGER ITYP
      LOGICAL EX

#endif


      ! Microphysical time step
      dtB2M=DTsrc
      ! Set all tendencies to zero
      ldummy=execute_bulk2m_driver('all','2zero',mkx)
      ! Set thermodynamics
      tk0=TL                 ! Temperature, [K]
      qk0=QL                 ! Water vapor mixing ratio, [kq/kq]
      pk0=PL                 ! Pressure, [hPa]
      w0=VVEL*1.d-02            ! Large-scale velocity, [m/s]
      v0=WTURB !; v0=w3d(k,i,j) ! Sub-grid velocity, [m/s]
      r0=0.0  !RDTNDL(L)        ! tendency due to radiation, [K/s], not needed
      !      print*,"rad tendency",w0,v0,r0,L
      piby6=4.*atan(1.0)/6.0; row=1.e+03
      qrArray(1)=842.0e+00                      ! [m^(1-b)s]
      qrArray(2)=  0.8e+00                      ! [unitless]
      qrArray(3)=500.0e-06                      ! assumed mean rain diameter [m]
      qrArray(4)=piby6*row*qrArray(3)**3        ! mean rain mass [kg]
      qrArray(5)=1000.0e-06                     ! assumed rain conc [No/m^3]
      !      write(ou,*) 'k,   prebar(k)   vrain(k)    mrain(k)    nrain(k)'
      vrain=min(qrArray(1)*qrArray(3)**qrArray(2),9.2d0)  ! [m/s]
      !      write(6,*)"VRAIN",vrain,l
      mrain=100.d0*prebar/vrain             ! [kq/m^3]
      nrain=mrain/qrArray(4)                ! [No/m^3]
      RHO=1d5*PL/(RGAS*TL)
      mrain=1.d3*mrain/RHO                  !kg water/kg air
      nrain=1.d3*nrain/RHO                  !number/kg air
      ldummy=execute_bulk2m_driver('all' &
           ,tk0,qk0,pk0,w0,v0,r0)
      ! Set microphysics
      if(LHX.eq.LHE)  then
        mdrop=QCLX            ! drop content, [kg water/kg air]
        ndrop =NCLL*1.d6   !convert from cm-3 to m-3
        if (QCLX.eq.0.) ndrop=0.d0
        ncrys=0.d0;mcrys=0.0d0
      else
        WMXICE = QCIX
        mcrys=WMXICE         ! crys content, [kg water/kg air]
        ncrys=NCIL*1.d6     ! convert cm-3 to m-3; set at 0.1 l-1 = 1.d-4 cm-3
        if (QCIX.eq.0.) ncrys=0.d0
        ndrop=0.0d0;mdrop=0.0d0
      endif
      !
      ndrop_old=ndrop;mdrop_old=mdrop;ncrys_old=ncrys;mcrys_old=mcrys
      ndrop_new=0.0d0;mdrop_new=0.0d0;ncrys_new=0.0d0;mcrys_new=0.0d0
      nc_tnd=0.0d0;qc_tnd=0.0d0;ni_tnd=0.0d0;qi_tnd=0.0d0
      nc_tot=0.0d0;qc_tot=0.0d0;ni_tot=0.0d0;qi_tot=0.0d0
      !
      !** Convert from l-1 to cm-3====>  1 l^-1 = 10^-3 cm^-3 = 10^3 m^-3
      !      ldummy=execute_bulk2m_driver('all'
      !    *           ,ndrop,mdrop,ncrys,mcrys,'end')
#ifdef TRACERS_AMP
      do nm=1,nmodes
        naero(mkx,nm)=nactc(nm)
        !       if(l.eq.1) then
        !       if(nactc(l,nm).gt.1.)write(6,*)"Callmatrix",naero(mkx,nm)*1.e-6
        !        endif
      enddo
      ldummy=execute_bulk2m_driver('all' &
           ,ndrop,mdrop,ncrys,mcrys,naero,nmodes,'end',qr0=mrain, &
           nr0=nrain)
#endif
#ifdef TRACERS_TOMAS
!CCC
!Can *************************************************************************
!Can      CLOUD DROPLET CALCULATION
!Can *************************************************************************
!CCC
!CCC *** Input properties for parameterization
!CCC
      TOT_MI    = 0d0
      WPARC      = 0d0
      SMAX       = 0d0
      NACT       = 0d0
      REFF       = 0d0
!      CLDTAU     = 0d0
!      CLDTAUBL   = 0d0 ! I don't account BL case- yhl
!c$$$      QautP6     = 0d0
!c$$$      QautKK     = 0d0
!c$$$      QautMC     = 0d0
!c$$$      QautBH     = 0d0
!c$$$      QautGI     = 0d0
!c$$$      QautNS     = 0d0
!c$$$C
!C Get CCN properties
!C
!      avol(l) = axyp(i_debug,j_debug)*MA(i_debug,j_debug,l)/mair*
!!!   byam(l) = [m2/kg of air]
      boxvl = DXYPIJ*airm*mb2kg*rgas*TL  &
          /100./PL

      CALL getCCN (TM,BOXVL,TOT_MI,TOTi,TPi,MLi, &
         NCCNMx,NSECi)        ! Get CCN properties
!C
!C Call cloud microphysics
!C
      IF (LHX.EQ.LHE) THEN         ! Liquid clouds present
!CCC
!CCC *** Nenes & Seinfeld parameterization - calcuilate droplet number
!CCC

!Two options for Updrate velocity

!1. fixed as a constant

!         WPARCOcean = 0.15d0       ! Fix Ocean and terrestrial updrafts for now
!         WPARCEarth = 0.3d0
!         WPARC      = (1.d0-PEARTH)*WPARCOcean + PEARTH*WPARCEarth

!2. computed using EGCM
! TOMAS (Nov 2011) WPARC results in too high CDNC. So it reduced by 7 times (arbitrary)
!        WPARC=v0(mkx)/7. !wturb=sqrt(0.6667*EGCM(l,i,j))
! TOMAS (NOV 2013) WPARC now use large-scale vertical velocity (v0 is sub-grid scale velocity)
         WPARC=(VVEL*1.d-02+v0(mkx))/7.           !wturb=sqrt(0.6667*EGCM(l,i,j))
         WPARC=MAX(WPARC,0.0) ! make sure it is positive
         WPARC=MIN(WPARC,0.4) ! arbituary max
!End of updrate velocity option.


!        QLWC=WMX(L)/(FCLD + teeny)     !in-cloud dimensionless LWC
         QLWC=(QCLX+QCIX)/(FCLD + teeny)     !in-cloud dimensionless LWC or IWC
         QLWC=MIN(QLWC, 3.d-03)    !(upper limit for the QLWC)

         RHO=1.d5*PL/(RGAS*TL)
         RHOSI = RHO*1.d-3
         if(rhosi.eq.0.) print*,'zero rho',rho,pl,tl

         TPARC=tk0(mkx)
         PPARC=pk0(mkx)*100.d0  ! mbar to Pa
         IF (TOTi.GT.6.d7.and.WPARC.gt.0.) THEN  ! more than 60 particles per cc, call droplet
                                                 ! activation
            CALL CALCNd (TPARC,PPARC,TPi,MLi,NSECi,WPARC,NACT & ! Activate droplets
                 ,SMAX ,RHOSI,QLWC,EPSILON,AUTO,DIFFLWMR,DIFFEPS,pearth)
         ELSE
!YUNHA- The minimum NACT is set to 1 instead of 40.d6, which is used for old GISS-TOMAS model.
            NACT = 1.0 !ndrop(mkx) ! 40.d6      ! Minimum droplet number [#/m3]
            SMAX = 0.0001    ! Minimum supersaturation
         ENDIF

       ENDIF

       NACTL(mkx)=NACT
       CDNC_TOMAS=nactl(mkx)*1.e-6 !m-3 to cm-3

       ldummy=execute_bulk2m_driver('all' &
            ,ndrop,mdrop,ncrys,mcrys,nactl,'end',qr0=mrain, &
           nr0=nrain)

#endif
      ! Make calls to calculate hydrometeors' growth rates due to microphysical
      ! processes
      ! content      :       [kg water/kg air/s]
      ! concentration:       [No/kg air/s]
      ! Activation of cloud droplets: prescribed AP spectrum
      !*** For the originial HM scheme with fixed distributions for amm. sulfate
      !       ldummy=execute_bulk2m_driver('hugh','drop_nucl',dtB2M,mkx)
#if defined(TRACERS_AMP)
      !*** Using the AMP_actv interface from MATRIX
      ldummy=execute_bulk2m_driver('matr','drop_nucl',dtB2M,mkx)
#elif defined(TRACERS_TOMAS)
      !*** Using the TOMAS_actv interface from TOMAS
      ldummy=execute_bulk2m_driver('toma','drop_nucl',dtB2M,mkx)
#else
      !*** Use this if using the Lohmann or Gultepe scheme for mass to number
      ldummy=execute_bulk2m_driver('all' &
           ,ndrop,mdrop,ncrys,mcrys,'end',qr0=mrain,nr0=nrain)
      OLDCDNC=OLDCDN*1.d6  !convert from cm-3 to m-3
      NEWCDNC=NEWCDN*1.d6  !convert from cm-3 to m-3
      ldummy=execute_bulk2m_driver('gult','drop_nucl',dtB2M,mkx, &
           OLDCDNC,NEWCDNC)
#endif
      ! Droplets' autoconversion: Beheng (concentration and content)
      !       ldummy=execute_bulk2m_driver('hugh','drop_auto',dtB2M,mkx)
      ! Droplets' autoconversion: Seifert and Beheng (concentration and content)
      !       ldummy=execute_bulk2m_driver('beheng','drop_auto',dtB2M,mkx)
      ! Freezing of cloud droplets (contact and immersion)
      ldummy=execute_bulk2m_driver('hugh','drop_frzn',dtB2M,mkx)
      ! Crystal nucleation: Ncrys=anuc(k)*wef**bnuc(k)
      ldummy=execute_bulk2m_driver('hugh','crys_nucl',dtB2M,mkx)
      ! Numerous processes of water-water,water-ice, ice-ice interaction,
      ! condensation/evaporation/deposition/sublimation, ice multiplication
      ! and sedimention are ready to be called. There are only a few examples:
      !        ldummy=execute_bulk2m_driver('hugh','drop_rain',dtB2M,mkx)
      !        ldummy=execute_bulk2m_driver('hugh','drop_snow',dtB2M,mkx)
      !        ldummy=execute_bulk2m_driver('hugh','crys_auto',dtB2M,mkx)
      !        ldummy=execute_bulk2m_driver('hugh','crys_snow',dtB2M,mkx)
      !        ldummy=execute_bulk2m_driver('hugh','crys_cond',dtB2M,mkx)
      !        ldummy=execute_bulk2m_driver('hugh','snow_melt',dtB2M,mkx)
      !
      ! In this chain of events the very last call, which applies saturation
      ! adjustment to keep environment about water saturation, is supposed to be
      !        ldummy=execute_bulk2m_driver('hugh','drop_cond',dtB2M,mkx)
      !
      ! To ensure calculated growth rates don't lead to negative contents:
      ! Previous call ('drop_cond') HAS to be used.
      !        ldummy=execute_bulk2m_driver('all','make_balance',mkx)
      ! Otherwise, the caller is responsible to handle the problem.
      !
      ! To get particular growth rate:
      !        rablk = execute_bulk2m_driver('get','rate','mprc')
      ! return value "rablk" is real*8 array whose dimension is equal to "mkx"
      !
      ! To get tendencies of temperature, water vapor mixing ratio or
      ! tendencies of concentration/content of particular hydrometeor:
      !        rablk = execute_bulk2m_driver('get','tnd','qc_tnd')
      ! return value "rablk" is real*8 array whose dimension is equal to "mkx"
      !
      ! To get parameters of hydrometeors size distributions:
      !        rablk = execute_bulk2m_driver('get','val','ec')
      ! return value "rablk" is real*8 array whose dimension is equal to "mkx"
      !
      ! To update concentration and contents due to uncommented processes
      ! If "make_balance" was used:
      !        ndrop=ndrop+dtB2M*execute_bulk2m_driver('get','tnd','nc_tnd')
      !        mdrop=mdrop+dtB2M*execute_bulk2m_driver('get','tnd','qc_tnd')
      !        ncrys=ncrys+dtB2M*execute_bulk2m_driver('get','tnd','ni_tnd')
      !        mcrys=mcrys+dtB2M*execute_bulk2m_driver('get','tnd','qi_tnd')

      !
      ! This is our case
      ! Without "make_balance":
      ! Droplet concentration
      !       if(l.eq.1) then
      !        npccn=execute_bulk2m_driver('get','npccn')
      !        nprc=execute_bulk2m_driver('get','nprc')
      !        nnuccc=execute_bulk2m_driver('get','nnuccc')
      !        nnucci=execute_bulk2m_driver('get','nnucci')
      !         if(npccn(l).gt.1.)write(6,*)"check BLK",ndrop*1.e-6,
      !    *npccn*1.e-6,nprc*1.e-6,nnuccc*1.e-6,nnucci*1.e-6
      !       endif

      ndrop=ndrop+( &
                                !       npccn              ! change n droplets activation
           +execute_bulk2m_driver('get','npccn') &
                                !       nprc               ! change n autoconversion of droplets:
           -execute_bulk2m_driver('get','nprc') &
                                !       nnuccc             ! change n due to con droplets freez
           -execute_bulk2m_driver('get','nnuccc') &
                                !       nnucci             ! change n due to imm droplets freez
           -execute_bulk2m_driver('get','nnucci') &
                                !
           )*dtB2M
      !
      ! Droplet content
      mdrop=mdrop+( &
                                !       mpccn              ! change q droplets activation
           +execute_bulk2m_driver('get','mpccn') &
                                !       mprc               ! change q autoconversion of droplets:
           -execute_bulk2m_driver('get','mprc') &
                                !       mnuccc             ! change q due to con droplets freez
           -execute_bulk2m_driver('get','mnuccc') &
                                !       mnucci             ! change q due to imm droplets freez
           -execute_bulk2m_driver('get','mnucci') &
                                !
           )*dtB2M
      !
      ! Crystal concentration
      ncrys=ncrys+( &
                                !       nnuccc             ! change n due to contact droplets freez
           +execute_bulk2m_driver('get','nnuccc') &
                                !       nnucci            ! change n due to immersion droplets freez
           +execute_bulk2m_driver('get','nnucci') &
                                !       nnuccd            ! change n freezing aerosol (prim ice nuc)
           +execute_bulk2m_driver('get','nnuccd') &
           )*dtB2M &
                                !      nnucmd        ! change n cond freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','nnucmd') &
                                !      nnucmt        ! change n cont freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','nnucmt')
      !
      ! Crystal content
      mcrys=mcrys+( &
                                !       mnuccc             ! change q due to con droplets freez
           +execute_bulk2m_driver('get','mnuccc') &
                                !       mnucci             ! change q due to imm droplets freez
           +execute_bulk2m_driver('get','mnucci') &
                                !       mnuccd            ! change q freezing aerosol (prim ice nuc)
           +execute_bulk2m_driver('get','mnuccd') &
           )*dtB2M &
                                !      mnucmd        ! change q cond freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','mnucmd') &
                                !      mnucmt        ! change q cont freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','mnucmt')
      !
      if_balance: if( (ndrop(mkx) .lt. 0) .or. (mdrop(mkx) .lt. 0)  .or. &
           (ncrys(mkx) .lt. 0) .or. (mcrys(mkx) .lt. 0)) then
        !
        if(lSCM) then
          write(6,*)"stop BLK: ndrop_old,mdrop_old,ncrys_old,mcrys_old" &
               ,ndrop_old*1.e-6,mdrop_old*1.e+3,ncrys_old*1.e-3,mcrys_old*1.e+3
          !
          write(6,*)"stop BLK: ndrop,mdrop,ncrys,mcrys" &
               ,ndrop*1.e-6,mdrop*1.e+3,ncrys*1.e-3,mcrys*1.e+3
        endif
        !
        ! No/m^3
        !
        npccn  =             & ! change n droplets activation
             +execute_bulk2m_driver('get','npccn')*dtB2M
        nprc   =             & ! change n autoconversion of droplets:
             -execute_bulk2m_driver('get','nprc')*dtB2M
        nnuccc =             & ! change n due to con droplets freez
             -execute_bulk2m_driver('get','nnuccc')*dtB2M
        nnucci =             & ! change n due to imm droplets freez
             -execute_bulk2m_driver('get','nnucci')*dtB2M

        nc_tot = npccn + nprc + nnuccc + nnucci
        !
        ! No/cc
        !
        if(lSCM) then
          write(6,*)"stop BLK: ndrop_old,nc_tot,ndrop" &
               ,ndrop_old*1.e-6,nc_tot*1.e-6,ndrop*1.e-6
          write(6,*)"stop BLK: npccn,nprc,nnuccc,nnucci" &
               ,npccn*1.e-6,nprc*1.e-6,nnuccc*1.e-6,nnucci*1.e-6
        endif
        !
        ! kg/kg
        !
        mpccn  =              & ! change q droplets activation
             execute_bulk2m_driver('get','mpccn')*dtB2M
        mprc   =              & ! change q autoconversion of droplets:
             -execute_bulk2m_driver('get','mprc')*dtB2M
        mnuccc =              & ! change q due to con droplets freez
             -execute_bulk2m_driver('get','mnuccc')*dtB2M
        mnucci =              & ! change q due to imm droplets freez
             -execute_bulk2m_driver('get','mnucci')*dtB2M

        qc_tot = mpccn + mprc + mnuccc + mnucci
        !
        ! g/kg
        !
        if(lSCM) then
          write(6,*)"stop BLK: mdrop_old,qc_tot,mdrop" &
               ,mdrop_old*1.e+3,qc_tot*1.e+3,mdrop*1.e+3
          write(6,*)"stop BLK: mpccn,mprc,mnuccc,mnucci" &
               ,mpccn*1.e+3,mprc*1.e+3,mnuccc*1.e+3,mnucci*1.e+3
        endif
        !
        ! No/m^3
        !
        nnuccc  =             & ! change n due to contact droplets freez
             +execute_bulk2m_driver('get','nnuccc')*dtB2M
        nnucci  =             & ! change n due to immersion droplets freez
             +execute_bulk2m_driver('get','nnucci')*dtB2M
        nnuccd  =             & ! change n freezing aerosol (prim ice nuc)
             +execute_bulk2m_driver('get','nnuccd') ! *dtB2M
        nnucmd =              & ! change n cond freezing Meyer's (prim ice nuc)
             +execute_bulk2m_driver('get','nnucmd') ! *dtB2M
        nnucmt =              & ! change n cont freezing Meyer's (prim ice nuc)
             +execute_bulk2m_driver('get','nnucmt') ! *dtB2M

        ni_tot = nnuccc + nnucci + nnuccd + nnucmd + nnucmt
        !
        ! No/l
        !
        if(lSCM) then
          write(6,*)"stop BLK: ncrys_old,ni_tot,ncrys" &
               ,ncrys_old*1.e-3,ni_tot*1.e-3,ncrys*1.e-3
          write(6,*)"stop BLK: nnuccc,nnucci,nnuccd,nnucmd,nnucmt" &
               ,nnuccc*1.e-3,nnucci*1.e-3,nnuccd*1.e-3,nnucmd*1.e-3 &
               ,nnucmt*1.e-3
        endif
        !
        ! kg/kg
        !
        mnuccc  =             & ! change q due to contact droplets freez
             +execute_bulk2m_driver('get','mnuccc')*dtB2M
        mnucci  =             & ! change q due to immersion droplets freez
             +execute_bulk2m_driver('get','mnucci')*dtB2M
        mnuccd  =             & ! change q freezing aerosol (prim ice nuc)
             +execute_bulk2m_driver('get','mnuccd') ! *dtB2M
        mnucmd =              & ! change q cond freezing Meyer's (prim ice nuc)
             +execute_bulk2m_driver('get','mnucmd') ! *dtB2M
        mnucmt =              & ! change q cont freezing Meyer's (prim ice nuc)
             +execute_bulk2m_driver('get','mnucmt') ! *dtB2M

        qi_tot = mnuccc + mnucci + mnuccd + mnucmd + mnucmt
        !
        ! g/m^3
        !
        if(lSCM) then
          write(6,*)"stop BLK: mcrys_old,qi_tot,mcrys" &
               ,mcrys_old*1.e+3,qi_tot*1.e+3,mcrys*1.e+3
          write(6,*)"stop BLK: mnuccc,mnucci,mnuccd,mnucmd,mnucmt" &
               ,mnuccc*1.e+3,mnucci*1.e+3,mnuccd*1.e+3,mnucmd*1.e+3 &
               ,mnucmt*1.e+3
        endif
        !
        ! balanced tendecies:
        !
        ldummy=execute_bulk2m_driver('all','make_balance',mkx)
        !
        nc_tnd=dtB2M*execute_bulk2m_driver('get','tnd','nc_tnd')
        if(lSCM) then
          write(6,*)"stop BLK:00: nc_tnd",nc_tnd*1.e-6
        endif
        !
        npccn  =             & ! change n droplets activation
             +execute_bulk2m_driver('get','npccn')*dtB2M
        nprc   =             & ! change n autoconversion of droplets:
             -execute_bulk2m_driver('get','nprc')*dtB2M
        nnuccc =             & ! change n due to con droplets freez
             -execute_bulk2m_driver('get','nnuccc')*dtB2M
        nnucci =             & ! change n due to imm droplets freez
             -execute_bulk2m_driver('get','nnucci')*dtB2M

        nc_tnd = npccn + nprc + nnuccc + nnucci
        if(lSCM) then
          write(6,*)"stop BLK:01: nc_tnd",nc_tnd*1.e-6
        endif
        !
        qc_tnd=dtB2M*execute_bulk2m_driver('get','tnd','qc_tnd')
        if(lSCM) then
          write(6,*)"stop BLK:00: qc_tnd",qc_tnd*1.e+3
        endif
        !
        mpccn  =              & ! change q droplets activation
             execute_bulk2m_driver('get','mpccn')*dtB2M
        mprc   =              & ! change q autoconversion of droplets:
             -execute_bulk2m_driver('get','mprc')*dtB2M
        mnuccc =              & ! change q due to con droplets freez
             -execute_bulk2m_driver('get','mnuccc')*dtB2M
        mnucci =              & ! change q due to imm droplets freez
             -execute_bulk2m_driver('get','mnucci')*dtB2M

        qc_tnd = mpccn + mprc + mnuccc + mnucci
        if(lSCM) then
          write(6,*)"stop BLK:01: qc_tnd",qc_tnd*1.e+3
        endif
        !
        ni_tnd=dtB2M*execute_bulk2m_driver('get','tnd','ni_tnd')
        if(lSCM) then
          write(6,*)"stop BLK:01: ni_tnd",ni_tnd*1.e-3
        endif
        !
        qi_tnd=dtB2M*execute_bulk2m_driver('get','tnd','qi_tnd')
        if(lSCM) then
          write(6,*)"stop BLK:01: qi_tnd",qi_tnd*1.e+3
        endif
        !
        ndrop_new = ndrop_old + nc_tnd
        mdrop_new = mdrop_old + qc_tnd
        ncrys_new = ncrys_old + ni_tnd
        mcrys_new = mcrys_old + qi_tnd
        !
        ndrop_blk = execute_bulk2m_driver('get','val','nc')
        mdrop_blk = execute_bulk2m_driver('get','val','qc')
        ncrys_blk = execute_bulk2m_driver('get','val','ni')
        mcrys_blk = execute_bulk2m_driver('get','val','qi')
        !
        ndrop_res = ndrop_blk + nc_tnd
        mdrop_res = mdrop_blk + qc_tnd
        ncrys_res = ncrys_blk + ni_tnd
        mcrys_res = mcrys_blk + qi_tnd
        !
        if(wSCM) then
          write(6,*) &
               "stop BLK: ndrop_old,nc_tnd,ndrop_new" &
               ,ndrop_old*1.e-6,nc_tnd*1.e-6,ndrop_new*1.e-6
          !
          write(6,*) &
               "stop BLK: mdrop_old,qc_tnd,mdrop_new" &
               ,mdrop_old*1.e+3,qc_tnd*1.e3,mdrop_new*1.e+3
          !
          write(6,*) &
               "stop BLK: ncrys_old,ni_tnd,ncrys_new" &
               ,ncrys_old*1.e-3,ni_tnd*1.e-3,ncrys_new*1.e-3
          !
          write(6,*) &
               "stop BLK: mcrys_old,qi_tnd,mcrys_new" &
               ,mcrys_old*1.e+3,qi_tnd*1.e3,mcrys_new*1.e+3
          !
          write(6,*) &
               "stop BLK: ndrop_blk,nc_tnd,ndrop_res" &
               ,ndrop_blk*1.e-6,nc_tnd*1.e-6,ndrop_res*1.e-6
          !
          write(6,*) &
               "stop BLK: mdrop_blk,qc_tnd,mdrop_res" &
               ,mdrop_blk*1.e+3,qc_tnd*1.e3,mdrop_res*1.e+3
          !
          write(6,*) &
               "stop BLK: ndrop_old,ndrop,ndrop_new,nc_tot,nc_tnd" &
               ,ndrop_old*1.e-6,ndrop_new*1.e-6,ndrop*1.e-6 &
               ,nc_tot*1.e-6,nc_tnd*1.e-6
          !
          write(6,*) &
               "stop BLK: mdrop_old,mdrop,mdrop_new,qc_tot,qc_tnd" &
               ,mdrop_old*1.e+3,mdrop_new*1.e+3,mdrop*1.e+3 &
               ,qc_tot*1.e+3,qc_tnd*1.e+3
        endif
        !
        ! output for standalone internal variables
        !
        if(lSCM) then
          write(iuo,*) 'dtB2M'
          write(iuo,*) dtB2M

          write(iuo,*) 'wmx wmxice tl ql pl svlhxl lhx '
          write(iuo,*) qclx
          write(iuo,*) wmxice
          write(iuo,*) tl
          write(iuo,*) ql
          write(iuo,*) pl
          write(iuo,*) svlhxl
          write(iuo,*) lhx

          write(iuo,*) 'ndrop_old mdrop_old ncrys_old mcrys_old'
          write(iuo,*) ndrop_old
          write(iuo,*) mdrop_old
          write(iuo,*) ncrys_old
          write(iuo,*) mcrys_old

          write(iuo,*) 'nc_tnd qc_tnd ni_tnd qi_tnd'
          write(iuo,*) nc_tnd
          write(iuo,*) qc_tnd
          write(iuo,*) ni_tnd
          write(iuo,*) qi_tnd

          write(iuo,*) 'ndrop_new mdrop_new ncrys_new mcrys_new'
          write(iuo,*) ndrop_new
          write(iuo,*) mdrop_new
          write(iuo,*) ncrys_new
          write(iuo,*) mcrys_new

          write(iuo,*) 'ndrop mdrop ncrys mcrys'
          write(iuo,*) ndrop
          write(iuo,*) mdrop
          write(iuo,*) ncrys
          write(iuo,*) mcrys

          write(iuo,*) 'ndrop_res mdrop_res ncrys_res mcrys_res'
          write(iuo,*) ndrop_res
          write(iuo,*) mdrop_res
          write(iuo,*) ncrys_res
          write(iuo,*) mcrys_res

          write(iuo,*) 'npccn nprc nnuccc nnucci nc_tot'
          write(iuo,*) npccn
          write(iuo,*) nprc
          write(iuo,*) nnuccc
          write(iuo,*) nnucci
          write(iuo,*) nc_tot

          write(iuo,*) 'mpccn mprc mnuccc mnucci qc_tot'
          write(iuo,*) mpccn
          write(iuo,*) mprc
          write(iuo,*) mnuccc
          write(iuo,*) mnucci
          write(iuo,*) qc_tot

          write(iuo,*) 'nnuccd nnucmd nnucmt ni_tot'
          write(iuo,*) nnuccd
          write(iuo,*) nnucmd
          write(iuo,*) nnucmt
          write(iuo,*) ni_tot

          write(iuo,*) 'mnuccd mnucmd mnucmt qi_tot'
          write(iuo,*) mnuccd
          write(iuo,*) mnucmd
          write(iuo,*) mnucmt
          write(iuo,*) qi_tot
          !
        endif
        !
        if( (ndrop_res(mkx) .lt. 0) .or. (mdrop_res(mkx) .lt. 0)  .or. &
             (ncrys_res(mkx) .lt. 0) .or. (mcrys_res(mkx) .lt. 0)) then
          !        call stop_model("BLK2MOM: Negative conc/cont...", 255)
          !         write(6,*)"We reached -ve con.",ndrop_res(mkx),mdrop_res(mkx),
          !     * ncrys_res(mkx), mcrys_res(mkx),l
          ndrop_res(mkx)=20.*1.d06
          mdrop_res(mkx)=1*1.d-06
          ncrys_res(mkx)=1*1.d-06
          mcrys_res(mkx)=1*1.d02
        else
          ndrop=ndrop_res;mdrop=mdrop_res;ncrys=ncrys_res;mcrys=mcrys_res
        endif
        !
      endif if_balance
      !
      ! To calculate "new" temperature and vapor mixing ratio:
      ldummy=execute_bulk2m_driver('tkqv','tk_qv',mkx)
      tk0new=tk0+dtB2M*execute_bulk2m_driver('get','tnd','tk_tnd')
      qk0new=qk0+dtB2M*execute_bulk2m_driver('get','tnd','qv_tnd')
      !
      ! At this point you have 2 phases separately.
      ! Almost all processes are switched off, but you can calculate also
      ! accreation of droplets/crystal  by rain/snow, for example, and use
      ! rain/snow as diagnostic variables. But you need one additional
      ! long-storage array to keep ice crystal hydrometeor content as a minimum
      !
      !     IF(LHX.EQ.LHE)  THEN
      !        WMX(L)=mdrop(mkx)
      !      ELSE
      !        WMX(L)=mcrys(mkx)
      !      ENDIF
      !
      ! GCM logics...........  SNd0, SNdL [No/cc; ]SNdI Units also in /cc
      !
      !      SNdI=ncrys(mkx)*1.0d-6          ! ncrys, [No/m^3]
      SNdI = 0.06417127d0
      !      if(SNdI.gt.0.) write(6,*)"ICE CRY",SNdI, SNdI/dtB2M
      if(SNdI.gt.1.d0) SNdI=1.d0      !try to limit to 1000 /l
      SNd=ndrop(mkx)*1.d-6                 ! ndrop, [No/m^3]
      !      if(SNd.gt.20.) write(6,*)"SM 12 CDNC",SNd   ,l
      !**** Old treatment to get CDNC for cloud changes within time steps
      DCLD = FCLD-CLDSAVL ! cloud fraction change
      !** If previous time step is clear sky
      if(CLDSAVL.eq.0.) then
        SNd=SNd
      elseif (DCLD.le.0.d0) then
        SNd=NCLL
      elseif(DCLD.gt.0.d0) then
        SNd=( (NCLL*CLDSAVL) + (SNd*DCLD) )/FCLD
      endif
      !* If using an alternate definition for QAUT
      rablk=execute_bulk2m_driver('get','mprc')
      QAUT_B2M=rablk(mkx)

      SCDNCW=SNd      ! we have already passed the grid box value
      SCDNCI=SNdI
      if (SCDNCI.le.0.0d0) SCDNCI=teeny         !set min ice crystal, do we need this, please check
      if (SCDNCW.le.20.d0) SCDNCW=20.d0         !set min CDNC, sensitivity test
      !     if(SCDNCW.gt.2000.) write(6,*)"PROBLEM",SCDNCW,L
      if (SCDNCW.ge.1400.d0) SCDNCW=1400.d0     !set max CDNC, sensitivity test
      !     if (SNd.gt.20.) write(6,*)"CDNC LSS",SCDNCW,SNd,L

  end subroutine cld_aer_cdnc_block1

  subroutine cld_aer_cdnc_block2( &
         lhx,fcld,dtb2m, &
         dsu, &
         cleara,vvel,CLDSAV0,qclx,qcix,sme,ncll,ncil,tl,ql,pl,wturb,wmxice, &
         snd_l, & ! only for tracers_amp
         rbeta, & ! output
         scdncw,scdnci &
      )
      use constant, only : lhe,teeny
      implicit none

      real*8 :: vvel,lhx,fcld,scdncw,scdnci,dtb2m,snd_l,rbeta
      real*8 :: cleara,CLDSAV0,qclx,qcix,sme,ncll,ncil,tl,ql,pl,wturb,wmxice
      real*8 :: dsu(sntm)

      ! local variables
      real*8 CDNL1,CDNL0,NEWCDN,OLDCDN,SNd,qcx,SNdi
      real*8 NEWCLD,SAVCLD
      real*8 :: Repsis,Repsi
#ifdef TRACERS_AMP
    real*8 :: naero (mkx,nmodes)
#endif


      NEWCLD = 1.-CLEARA
      SAVCLD = CLDSAV0 ! from prev. timestep

!@auth Menon for CDNC prediction
#if defined(TRACERS_AMP)
      NCLL=SNd_L
      !NCIL=SNdi
#elif defined(TRACERS_TOMAS)
       NCLL=SNd_L
       !NCIL=SNdi
#else
      if( LHX.eq.LHE )then
        QCX = QCLX
      else
        QCX = QCIX
      endif
      call GET_CDNC_UPD(LHX,QCX,FCLD,NEWCLD, &
           SAVCLD,VVEL,SME,DSU,NCLL, &
           CDNL0,CDNL1)
      NCLL = CDNL1
      SNd=CDNL1
      !** Pass old and new cloud droplet number
      NEWCDN=SNd
      OLDCDN=CDNL0
      !     if (L.eq.1)write(6,*)"BLK_2M NUPD",NEWCDN,OLDCDN
#endif

      ! Update thermo if environment was changed
      tk0=TL                 ! Temperature, [K]
      qk0=QL                 ! Water vapor mixing ratio, [kq/kq]
      pk0=PL                 ! Pressure, [hPa]
      w0=VVEL*1.d-02           ! Large-scale velocity, [m/s]
      v0=WTURB !; v0=w3d(k,i,j) ! Sub-grid velocity, [m/s]
      r0= 0.0  !RDTNDL(L)               ! T tendency due to radiation, [K/s]
      ldummy=execute_bulk2m_driver('all' &
           ,tk0,qk0,pk0,w0,v0,r0)
      ! Update micro if contents were changed
      !        mdrop=WMX(L)            ! drop content, [kg water/kg air]
      !        ndrop=mdrop/mw0         ! drop concent, [No/m3]
      !        mcrys=WMXICE(L)         ! crys content, [kg water/kg air]
      !        ncrys=mcrys/mi0         ! crys concent, [No/m3]
      if(LHX.eq.LHE)  then
        mdrop =QCLX
        ndrop= NCLL*1.d6  !mdrop/mw0         ! drop concent, [No/m3]
        if(QCLX.eq.0.) ndrop=0.0
        ncrys = 0.
        mcrys = 0.
      else
        mcrys =QCIX
        WMXICE = QCIX
        ncrys= NCIL*1.d6  !mcrys/mi0         ! crystal concent, [No/m3]
        if(QCIX.eq.0.) ncrys=0.0
        ndrop = 0.
        mdrop = 0.
      endif
      !      if(L.eq.1)write(6,*)"5th check BLK_2M",
      !    *WMX(L),NCLL(L),NCIL(L)
      !
      ldummy=execute_bulk2m_driver('all' &
           ,ndrop,mdrop,ncrys,mcrys,'end')
      ! Get new drop & crys concentration
      !     ldummy=execute_bulk2m_driver('surabi','GET_CDNC_UPD',dtB2M,mkx)
#if defined(TRACERS_AMP)
      !*** Using the AMP_actv interface from MATRIX
      ldummy=execute_bulk2m_driver('matr','drop_nucl',dtB2M,mkx)
#elif defined(TRACERS_TOMAS)
      !*** Using the AMP_actv interface from MATRIX
        ldummy=execute_bulk2m_driver('toma','drop_nucl',dtB2M,mkx)
#else
      !*** Call Lohmann's or Gultepe's scheme for CDNC
      OLDCDNC=OLDCDN*1.d6  !convert from cm-3 to m-3
      NEWCDNC=NEWCDN*1.d6  !convert from cm-3 to m-3
      ldummy=execute_bulk2m_driver('gult','drop_nucl',dtB2M,mkx, &
           OLDCDNC,NEWCDNC)
#endif
      rablk=execute_bulk2m_driver('get','value','nc') + ( &
           +execute_bulk2m_driver('get','npccn') &
           )*dtB2M
      SNd=rablk(mkx)*1.0d-6            ! ndrop, [No/m^3], SNdL, [No/cc]
      !     if(SNd.gt.20.) write(6,*)"Finally out",SNd, l
      !**** Old treatment to get CDNC for cloud changes within time steps
      DCLD = NEWCLD-SAVCLD ! cloud fraction change
      !** If previous time step is clear sky
      if(SAVCLD.eq.0.) then
        SNd=SNd
      elseif (DCLD.le.0.d0) then
        SNd=NCLL
      elseif(DCLD.gt.0.d0) then
        SNd=( (NCLL*SAVCLD) + (SNd*DCLD) )/NEWCLD
      endif
      rablk=execute_bulk2m_driver('get','value','ni') + ( &
                                !       nnuccc             ! change n due to contact droplets freez
           +execute_bulk2m_driver('get','nnuccc') &
                                !       nnucci            ! change n due to immersion droplets freez
           +execute_bulk2m_driver('get','nnucci') &
                                !       nnuccd            ! change n freezing aerosol (prim ice nuc)
           +execute_bulk2m_driver('get','nnuccd') &
           )*dtB2M &
                                !      nnucmd        ! change n cond freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','nnucmd') &
                                !      nnucmt        ! change n cont freezing Meyer's (prim ice nuc)
           +execute_bulk2m_driver('get','nnucmt')

      !      SNdI=rablk(mkx)*1.0d-6             ! from ncrys [No/m^3] to SNdI in [No/cc]
      SNdI = 0.06417127d0
      if(SNdI.gt.1.d0) SNdI=1.d0      !try to limit to 1000 /l
      NCLL = SNd
      NCIL = SNdI
#ifdef TRACERS_AMP
      !nactc(l,1:nmodes) =  naero(mkx,1:nmodes)
      !      do nm=1,nmodes
      !        if(nactc(l,nm).gt.0.)
      !    *   write(6,*)"NMOD1",nactc(l,nm),l,nm
      !       enddo
#endif
      !      if(L.eq.1) write(6,*)"6_LO check BLK_2M",SNd,SNdI
      ! To get effective radii in micron
      rablk=execute_bulk2m_driver('get','value','ec')  ! [micron]

      SCDNCW=SNd
      SNdI = 0.06417127d0
      SCDNCI=SNdI
      if (SCDNCW.le.20.d0) SCDNCW=20.d0   !set min CDNC sensitivity test
      !     If (SCDNCI.le.0.06d0) SCDNCI=0.06417127d0   !set min ice crystal
      if (SCDNCI.le.0.0d0) SCDNCI=teeny           !set min ice crystal
      if(SCDNCW.gt.1400.d0) SCDNCw=1400.d0
      !     if (SCDNCW.gt.20.) write(6,*) "SCND CDNC",SCDNCW,NCLL(l),l

      if(lhx.eq.lhe) then
        !** Using the Liu and Daum paramet
        !** for spectral dispersion effects on droplet size distribution
        Repsi=1.d0 - 0.7d0*exp(-0.003d0*SCDNCW)
        Repsis=Repsi*Repsi
        Rbeta=(((1.d0+2.d0*Repsis)**0.667d0))/((1.d0+Repsis)**0.333d0)
      endif

  end subroutine cld_aer_cdnc_block2

end module cld_aer_cdnc_mod
