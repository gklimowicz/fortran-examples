#include "rundeck_opts.h"

      MODULE CH4_SOURCES
      USE TRACER_COM
      implicit none
!@var CH4_src CH4 surface sources and sinks (kg/s)
      integer, parameter :: nch4src=14
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: CH4_src
!@var frqlos chemical loss rate for methane in troposphere
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: frqlos
      END MODULE CH4_SOURCES

      MODULE PRATHER_CHEM_COM
!@sum Variables for chemical tracer routines that were provided by
!@+    Michael Prather.  The chemistry is parameterized as frequencies
C**** Stratospheric chemistry: consolidated commons
C**** These variables are used by both ozone and strat chem routines
      USE RESOLUTION, only: jm,lm
      INTEGER, ALLOCATABLE, DIMENSION(:) :: jlatmd
      real*8 p0l(lm+1)
!@dbparam NSTRTC # of strat chem layers (counting top down) (=12 for M23)
      integer :: nstrtc=12,lmtc

      contains
      subroutine set_prather_constants
      Use RESOLUTION, Only: JM,LM
      USE ATM_COM, only: pednl00
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE CH4_SOURCES
      USE Dictionary_mod
      implicit none
      real*8 yedge(GRID%J_STRT_HALO:GRID%J_STOP_HALO+1)
      real*8 yedge1,yedgen,xlatmd
      integer j,jxxx,l,lr

      INTEGER :: J_1,  J_0,  I_1,  I_0
      INTEGER :: J_1H, J_0H, I_1H, I_0H
      INTEGER :: IER, lmtc

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,  J_STOP     =J_1,
     *               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      call getDomainBounds(grid, I_STRT     =I_0,  I_STOP     =I_1,
     *               I_STRT_HALO=I_0H, I_STOP_HALO=I_1H)

      call sync_param("NSTRTC",NSTRTC)
C**** ESMF: This array is read in only
      lmtc = lm-nstrtc
      ALLOCATE(   frqlos(I_0:I_1,J_0:J_1,lmtc),
     *          STAT=IER)

C---calculate nearest latitude to std lats
      yedge1=-90.
      yedgen= 90.
c GISS-ESMF EXCEPTIONAL CASE
      do j=J_0,J_1+1
        yedge(j)=yedge1*float(jm+1-j)/float(jm) +
     *           yedgen*float(j-1)/float(jm)
      end do
      do j=J_0,J_1
        xlatmd = 0.5*(yedge(j)+yedge(j+1))
        jxxx = xlatmd/10. + 10.
        jlatmd(j) = min(18,max(1,jxxx))
      end do
C---just calculate average P(mbar) at edge of each level
C---Calculate average P(mbar) at edge of each level PLEVL(1)=Psurf
      do lr=1,lm+1
        p0l(lr) = pednl00(lm+2-lr)
      end do
      return

      end subroutine set_prather_constants
      end MODULE PRATHER_CHEM_COM


      MODULE TRACERS_MPchem_COM
!@sum Variables for Prather's Stratospheric chemistry loss model
      USE RESOLUTION, only: jm,lm
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
      implicit none
!@var nMPtable:  Number of tracers that will use the frequency
!@+    tables and share strat chem code
      integer :: nMPtable
!@param lz_schem Number of heights in stratchem tables
      integer, parameter :: lz_schem=20,lz_sx=lz_schem+7
!@var tscparm: Contains mean loss prequency in grid box
      real*8, allocatable, dimension(:,:,:,:) :: tscparm
!@var TLtrm,TLtzm,TLtzzm: loss freq and moments of loss freq from tables
      real*8, ALLOCATABLE, DIMENSION (:,:,:) :: tltrm,tltzm,tltzzm
!@var PS Used in STRT2M
      real*8 ps(lz_sx+1)

      contains

      SUBROUTINE STRATCHEM_SETUP
C**** Prather stratospheric chemistry
      USE FILEMANAGER, only: openunit,closeunit
      USE PRATHER_CHEM_COM, only: set_prather_constants
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
      use TRACER_COM, only: ntm
      use OldTracer_mod, only: trname,iMPtable
      implicit none

      integer n,j,k,m,iu
      character*80 titlch
      character*16 filein
      integer l,nl
      real*8    XPSD,XPSLM1,XPSL

      do n=1,ntm
        if (iMPtable(n) == 0) cycle
        filein = trim(trname(n))//'_TABLE'
        call openunit(filein,iu,.false.,.true.)
        read (iu,'(a)')   titlch
        if (AM_I_ROOT()) write(6,'(1x,a)') titlch
        do m=1,INT_MONTHS_PER_YEAR
          do j=1,18
            read(iu,'(20x,6e10.3/(8e10.3))')
     *          (tscparm(k,j,m,iMPtable(n)),k=lz_schem,1,-1)
          end do
        end do
        call closeunit(iu)
        if (AM_I_ROOT()) write(6,'(2A)') 'STRATCHEM TABLES READ for ',
     *                                   trim(trname(n))
      enddo

      call set_prather_constants

C**** This code moved from STRT2M to go faster
c-----------------------------------------------------------------------
c       set up std z* atmosphere: p = 1000 * 10**(-z*/16 km)
c       assume that stratospheric chemical parameters always start at
cc       52 km (N=27) scan downward from 52 km to 14 km (NX=20) by 2 km
c       58 km (N=30) scan downward from 58 km to 10 km (NX=25) by 2 km
c       intervals, constant >58km
c-------- N.B. F(@30km) assumed to be constant from 29-31 km (by mass)
      nl = lz_sx
      XPSD       = 10.D0 **(-0.125D0)      !Z=2 km
      XPSLM1     = 1000.D0                 !Z=0 km, P=1000 mb
      PS(1)      = 1000.D0
      DO L = 2,NL
        XPSL     = XPSLM1 *XPSD
        PS(L)    = 0.5D0 *(XPSLM1 +XPSL) !Z=2*I km, P=1000*10**(-Z/16)
        XPSLM1   = XPSL
      ENDDO
      PS(NL+1)   = 0.D0
      return

      end SUBROUTINE STRATCHEM_SETUP
      end MODULE TRACERS_MPchem_COM


      SUBROUTINE Strat_chem_Prather(ns,n)
!@sum Strat_chem_Prather calculates stratospheric chemistry for
!@+        N2O, CFC and CH4
!@auth Michael Prather (J.Lerner adopted code)
!@var nsc index for chemical tracer number (to access tables)
!@var T0L is the amount 'lost' to chemistry this time step
!@var facbb: APPLY AN AD-HOC FACTOR TO BRING CH4 INTO BALANCE
      USE CONSTANT, only: by3
      USE RESOLUTION, only: im,jm,lm
      USE MODEL_COM, only: dtsrc
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds
      USE GEOM, only: imaxj
      USE QUSDEF, only : mz,mzz
      use OldTracer_mod, only: itime_tr0,trname,tcscale,iMPtable
      USE TRACER_COM, only: trm, trmom
      USE TRACERS_MPchem_COM, only: tltrm,tltzm,tltzzm
      USE PRATHER_CHEM_COM, only: nstrtc
      USE FLUXES, only: tr3Dsource
      implicit none
      integer i,j,l,lr,n,ns,najl,nsc,lmtc
      real*8, parameter :: by7=1./7.d0
      real*8 told(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &            GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm)
      real*8 f0l,f1l,f2l,g0l,g1l,g2l,t0l,t1l,t2l,facbb

      INTEGER :: J_1, J_0, I_0, I_1
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      facbb = 1.
      if (trname(n).eq.'CH4') facbb = (40.d0/25.73d0)*(40.d0/35.177d0)
C-----STRATOSPHERIC LOSS occurs in top NSTRTC layers
C-----ALLOW FOR ONLY ONE SET OF TSCPARMs FOR STRATOSPHERIC LOSS,
C-----uses TCSCALE for different tracers to scale loss
C-----uses S.O.M. formulation for vertical losses
C-----NOTE that TLTRM(J,LR,N) stored from top (=LM) down

      told(:,:,:) = trm(:,:,:,n)
      lmtc = lm-nstrtc
      do 150 l=lm,lmtc+1,-1
      lr = lm+1-l
      do 140 j=J_0,J_1
C-----TSCPARM->TLtrm contains mean loss freq in grid box:
        f0l = max(TLtrm(j,lr,iMPtable(n)),0d0)
        if (f0l.le.0.) go to 140
        f1l = tltzm(j,lr,iMPtable(n))  ! FOM of loss freq from tables
        f2l = tltzzm(j,lr,iMPtable(n))  ! SOM of loss freq from tables
        do 130 i=I_0,imaxj(j)
          if (trm(i,j,l,n).le.0.) go to 130
C------Couple the moments of the loss freq with moments of the tracer:
C       dSo  = dt*( So*Lo + Sz*Lz/3 + Szz*Lzz/5)
C       dSz  = dt*( Sz*Lo + So*Lz + 2(Sz*Lzz + Szz*Lz)/5)
C       dSzz = dt*( Szz*Lo + SoLzz +2Sz*Lz/3 + 2Szz*Lzz/7)
          t0l = f0l*trm(i,j,l,n) +
     *           f1l*trmom(mz,i,j,l,n)*by3+f2l*trmom(mzz,i,j,l,n)*.2d0
          if (t0l.lt.0.) go to 130
          t1l = f0l*trmom( mz,i,j,l,n) + f1l*trm(i,j,l,n) +
     *      0.4*(f2l*trmom(mz,i,j,l,n)    +f1l*trmom(mzz,i,j,l,n))
          t2l = f0l*trmom(mzz,i,j,l,n) + f2l*trm(i,j,l,n) +
     *      2.0*(f1l*trmom(mz,i,j,l,n)*by3+f2l*trmom(mzz,i,j,l,n)*by7)
C---calculate e-folding of tracer mass (trm) & scale moment losses
C---  to this change (T0L):
          g0l = t0l/trm(i,j,l,n)
          g1l = t1l/t0l
          g2l = t2l/t0l
          t0l = (1.0-exp(-g0l*dtsrc*tcscale(n)))*trm(i,j,l,n)
          t0l = t0l*facbb  ! APPLY AN AD-HOC FACTOR
          tr3Dsource(i,j,l,ns,n)=-t0l/dtsrc
cc          trm(i,j,l,n) = trm(i,j,l,n) - t0l
C**** moments are modified here since they are calculated specially
C**** moments ARE NOT modified in apply_tracer_3Dsource
          trmom( mz,i,j,l,n) = trmom( mz,i,j,l,n) - t0l*g1l
          trmom(mzz,i,j,l,n) = trmom(mzz,i,j,l,n) - t0l*g2l
  130     CONTINUE
  140     CONTINUE
  150   CONTINUE
cc      najl = jls_3Dsource(ns,n)
cc      do l=1,lm
cc      do j=J_0,J_1
cc      do i=I_0,imaxj(j)
cc        call inc_tajls(i,j,l,najl,trm(i,j,l,n)-told(i,j,l))
cc      end do
cc      end do
cc      end do
      return
      END SUBROUTINE Strat_chem_Prather


      SUBROUTINE STRTL
C**** This is called at the beginning of each month
C**** Prather strat chem
      USE RESOLUTION, only: jm,lm
      USE MODEL_COM, only: modelEclock
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds
      USE PRATHER_CHEM_COM, only: nstrtc,jlatmd,p0l
      USE TRACERS_MPchem_COM, only: tscparm,nMPtable,
     *    tltrm,tltzm,tltzzm,lz_schem,lz_sx,ps
      use OldTracer_mod, only: iMPtable
      implicit none
C-----------------------------------------------------------------------
C---monthly set up of chemical loss parameters
      real*8 strt0l(lm),strt1l(lm),strt2l(lm),strtx(lz_schem)
      real*8 f(lz_sx)
C-----------------------------------------------------------------------
C--tscparm(lz_schem,18,12,N) defined for 18 lats (85S, 75S, ...85N)
C --                   & 12 months
C----  do NOT interpolate, just pick nearest latitude
C---assume given MONTH = month #, JM=#lats, etc.
      integer n,j,jj,k,lr,jmon

      INTEGER :: J_1, J_0
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      jmon = modelEclock%getMonth()

      DO 800 N=1,nMPtable
        DO 700 J=J_0,J_1
          JJ = JLATMD(J)
          DO K=1,lz_schem
            STRTX(K) = tscparm(K,JJ,jmon,n)
          END DO
          CALL STRT2M(STRTX,lz_schem,STRT0L,STRT1L,STRT2L,P0L,NSTRTC
     *      ,ps,f,lz_sx)
C----store loss freq & moments in TLtrm/tzm/tzzm for exact
C---- CTM layers LM down
          DO 600 LR=1,NSTRTC
            TLtrm(J,LR,n) = STRT0L(LR)
            TLtzm(J,LR,n) = STRT1L(LR)
            TLtzzm(J,LR,n) = STRT2L(LR)
600       CONTINUE
700     CONTINUE
800   CONTINUE
      RETURN
      END SUBROUTINE STRTL


      SUBROUTINE Trop_chem_CH4(ns,n)
!@sum Trop_chem_CH4 calculates tropospheric chemistry for CH4
!@+     by applying a pre-determined chemical loss rate
!@auth Jean Lerner
      USE RESOLUTION, only: im,jm,lm
      USE MODEL_COM, only: nday,itime,dtsrc,modelEclock
      use TimeConstants_mod, only: INT_DAYS_PER_YEAR
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds, AM_I_ROOT,
     *  readt8_parallel,haveLatitude,broadcast
      USE GEOM, only: imaxj,byim
      USE PRATHER_CHEM_COM, only: nstrtc
      USE TRACER_COM
      USE CH4_SOURCES, only : frqlos
      USE FLUXES, only: tr3Dsource
      USE FILEMANAGER, only: openunit,closeunit,nameunit
      implicit none
      integer n,ns,i,j,l,FRQfile,lmtc,it_from_Dec29hr12,nrec,nr
      real*8 tune
      real*8, save :: taux=-1. ! hr of latest record read in (0.->8640.)
      parameter (tune = 445./501.)
      integer, save :: ifirst=1, nrec_cur=-1
      character*80 title
      real*8, dimension(:,:,:), allocatable :: arr_dummy_3d
      save FRQfile
      INTEGER :: J_1, J_0, I_0, I_1
      INTEGER :: J_1H, J_0H, I_1H, I_0H
      INTEGER :: IER
      integer :: jyear, jhour , jdate , month
      character*16, save :: FRQname='OHCH4_FRQ_interp'

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO

      lmtc = lm-nstrtc

      call modelEclock%get(hour=jhour, date=jdate, year=jyear ,
     &                     month=month)

      it_from_Dec29hr12 = mod(itime+5*NDAY/2,Nday*INT_DAYS_PER_YEAR)
      nrec = 1 + it_from_Dec29hr12/(5*Nday) ! record for current 5-day period

C**** Check whether chem.loss rate is up-to-date (updated every 5 days)
      if(nrec == nrec_cur) go to 550

C**** Create interpolated table for this resolution, position input file
      if (ifirst==1) then
        IF (AM_I_ROOT()) THEN
          call get_Trop_chem_CH4_freq(FRQname)
          call openunit(FRQname,FRQfile,.true.,.true.)
          do nr=1,nrec-1
            read(FRQfile)
          end do
          ifirst = 0
        END IF
        ! the following call is actually serving as MPI_Barrier
        ! do not remove it unless you know what you are doing
        call broadcast(grid, ifirst)
      end if

C**** Read chemical loss rate dataset (1 year of data, 5-day frequency)
      ALLOCATE(arr_dummy_3d(I_0:I_1,J_0:J_1,lmtc), STAT=IER)
      IF (AM_I_ROOT()) THEN
        read(FRQfile) title
        read (title,'(f10.0)') taux
        backspace(FRQfile)
      END IF
      CALL READT8_PARALLEL(grid,FRQfile,FRQname,arr_dummy_3d,0)
      nrec_cur = nrec

      frqlos = arr_dummy_3d
      DEALLOCATE (arr_dummy_3d)

      IF (AM_I_ROOT()) THEN
C**** near END OF YEAR, cycle to FIRST RECORD
        IF (taux==8640.) rewind FRQfile
        WRITE(6,'(2A,F10.0,I10,I7,a1,I2.2,a1,I2.2,a3,I3)')
     *    ' *** Chemical Loss Rates in Trop_chem_CH4 read for',
     *    ' taux,itime,yr/month/day/hour=',
     *     taux,itime,jyear,'/',month,'/',jdate,' hr',jhour
      END IF
C**** AVERAGE POLES
      if(haveLatitude(grid,J=1)) then
        do l=1,lmtc
          frqlos(1, 1,l) = sum(frqlos(:, 1,l))*byim
        end do
      end if
      if(haveLatitude(grid,J=JM)) then
        do l=1,lmtc
          frqlos(1,jm,l) = sum(frqlos(:,jm,l))*byim
        end do
      end if
C**** APPLY AN AD-HOC FACTOR TO BRING INTO BALANCE
      frqlos(:,:,:) = frqlos(:,:,:)*tune

C**** Apply the chemistry
  550 continue
      do l=1,lmtc
      do j=J_0,J_1
        do i=I_0,imaxj(j)
          tr3Dsource(i,j,l,ns,n) = -frqlos(i,j,l)*trm(i,j,l,n)
        end do
      end do
      end do
      return
      END SUBROUTINE Trop_chem_CH4


      MODULE LINOZ_CHEM_COM
C**** linoz with sol variability
!@sum Variables for linoz chemistry.  Original code was provided by
!@+    Michael Prather and Chris McLinden.
C**** Ozone tracer; linoz chemistry
C     NCTABLE=# of linoz tables (includes solar UV)
C     n_O3=tracer number for linoz O3
!@dbparam dsol describes portion of solar cycle being modeled for linoz
!@+      +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
      USE CONSTANT, only: mair
      USE RESOLUTION, only: im,jm,lm
      USE MODEL_COM, only: dtsrc
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
      USE ATM_COM, only: pednl00
      USE TRACER_COM, only: tr_mm
      USE PRATHER_CHEM_COM, only: set_prather_constants,nstrtc
      implicit none
      integer lmtc    !=11 for lm=23
!@param lz_linoz Number of heights in linoz tables
      integer, PARAMETER :: lz_linoz=25,nctable=7,lz_lx=lz_linoz+5
C****    lz_linoz heights, 18 lats, 12 months, nctable parameters
      real*8 TLPARM(lz_linoz,18,INT_MONTHS_PER_YEAR,nctable)
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: TLT0M, TLTZM, TLTZZM
      real*8 dsol
!@var PS,F Used in STRT2M
      real*8 PS(lz_lx+1)
C**** Harvard troposphere production and loss rates, deposition vel
!@var O3trop_Loss, O3trop_Prod, O3_DepVel: production and loss
!@+       rates, deposition vel from L. Mickley
      real*8, dimension(:,:,:,:), allocatable ::
     &     O3trop_Loss,O3trop_Prod ! (im,jm,lm,12)
      real*8, dimension(:,:,:), allocatable ::
     &     O3_DepVel !(im,jm,12)

      contains
      SUBROUTINE LINOZ_SETUP(n_O3)
C**** Needed for linoz chemistry
      USE FILEMANAGER, only: openunit,closeunit,nameunit
      USE DOMAIN_DECOMP_ATM, only: AM_I_ROOT,grid,readt_parallel
      use TimeConstants_mod, only: INT_MONTHS_PER_YEAR
      implicit none
      integer iu,i,j,k,l,m,n,n_O3,nl
      character*80 titlch
      real*8    XPSD,XPSLM1,XPSL
      real*8, dimension(:,:,:), allocatable :: arr_dummy_3d
      real*8, dimension(:,:), allocatable :: arr_dummy_2d

      call set_prather_constants
      lmtc = lm-nstrtc

      allocate(arr_dummy_3d(grid%i_strt_halo:grid%i_stop_halo,
     &                      grid%j_strt_halo:grid%j_stop_halo,lmtc))
      allocate(arr_dummy_2d(grid%i_strt_halo:grid%i_stop_halo,
     &                      grid%j_strt_halo:grid%j_stop_halo))

      call openunit('LINOZ_TABLE',iu,.false.,.true.)
      read (iu,'(a)')   titlch
      if (AM_I_ROOT()) write(6,'(1x,a)') titlch
      do n=1,nctable
        read (iu,'(a)')   titlch
        if (AM_I_ROOT()) write(6,'(1x,a)') titlch
        do m=1,INT_MONTHS_PER_YEAR
          do j=1,18
            read(iu,'(20x,6e10.3/(8e10.3))')
     *           (tlparm(k,j,m,n),k=lz_linoz,1,-1)
          end do
        end do
      end do
      if (AM_I_ROOT()) write(6,'(a)') ' linoz tables read'
      call closeunit(iu)

C****
C**** Harvard troposphere rate data (L.Mickley)
C****
C     Loss rates
      call openunit('LO3_Trop_loss',iu,.true.,.true.)
      do m = 1,12
        CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),arr_dummy_3d,0)
        O3trop_Loss(:,:,1:lmtc,m) = arr_dummy_3d
      enddo
      call closeunit(iu)

C     Production rates
      call openunit('LO3_Trop_prod',iu,.true.,.true.)
      do m = 1,12
        CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),arr_dummy_3d,0)
        O3trop_Prod(:,:,1:lmtc,m) = arr_dummy_3d
      enddo
      call closeunit(iu)

C     Deposition velocities
      call openunit('LINOZ_Dep_vel',iu,.true.,.true.)
      do m = 1,12
        CALL READT_PARALLEL(grid,iu,NAMEUNIT(iu),arr_dummy_2d,0)
        O3_DepVel(:,:,m) = arr_dummy_2d
      enddo
      call closeunit(iu)

      deallocate(arr_dummy_3d,arr_dummy_2d)

C**** This code moved from STRT2M to go faster
c-----------------------------------------------------------------------
c       set up std z* atmosphere: p = 1000 * 10**(-z*/16 km)
c       assume that stratospheric chemical parameters always start at
cc       52 km (N=27) scan downward from 52 km to 14 km (NX=20) by 2 km
c       58 km (N=30) scan downward from 58 km to 10 km (NX=25) by 2 km
c       intervals, constant >58km
c-------- N.B. F(@30km) assumed to be constant from 29-31 km (by mass)
      nl = lz_lx
      XPSD       = 10.D0 **(-0.125D0)      !Z=2 km
      XPSLM1     = 1000.D0                 !Z=0 km, P=1000 mb
      PS(1)      = 1000.D0
      DO L = 2,NL
        XPSL     = XPSLM1 *XPSD
        PS(L)    = 0.5D0 *(XPSLM1 +XPSL) !Z=2*I km, P=1000*10**(-Z/16)
        XPSLM1   = XPSL
      ENDDO
      PS(NL+1)   = 0.D0
      return

      end SUBROUTINE LINOZ_SETUP
      end MODULE LINOZ_CHEM_COM


      SUBROUTINE Trop_chem_O3(nsp,nsl,n)
c
c-----------------------------------------------------------------------
c   Troposphere is forced by Harvard tables
c-----------------------------------------------------------------------
c
      USE RESOLUTION, only: jm
      USE MODEL_COM, only: modelEclock,itime,dtsrc
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds
      USE CONSTANT, only : grav,rgas
      USE GEOM, only: imaxj,axyp
      USE ATM_COM, only: t,pmid,pk,pdsig
      USE TRACER_COM
      USE LINOZ_CHEM_COM, only: O3trop_Prod,O3trop_Loss,lmtc
      USE FLUXES, only: tr3Dsource
      implicit none
      integer i,j,l,n,nsp,nsl,jmon
      real*8 rprod,rloss,factor,tk
      real*8 dz != -dP/rhoG; rho=PRT; Deposition velocity
      INTEGER :: J_1, J_0, I_0, I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      jmon = modelEclock%getMonth()

C**** Convert from kg/cm3/s to kg
        do l=1,lmtc
        do j=J_0,J_1
        do i=I_0,imaxj(j)
          tk = t(i,j,l)*pk(l,i,j)            ! Temp in kelvin
          dz = pdsig(l,i,j)*rgas*tk/(pmid(l,i,j)*grav)   ! meters
          factor = dtsrc*axyp(i,j)*dz*1.d6    !/cm3->/m3
          rprod = O3trop_Prod(i,j,l,jmon)*factor     ! unit=kg
          rloss = O3trop_Loss(i,j,l,jmon)*factor*trm(i,j,l,n)
          if(trm(i,j,l,n) +(rprod-rloss).lt.0.) then
            write(6,'(a,3i3,4e14.3)') ' Negative O3 due to trop chem',
     *             i,j,l,trm(i,j,l,n),rprod,rloss,itime
            rloss = trm(i,j,l,n)+rprod
cc          scalmom = max(1.d0-rloss/(trm(i,j,l,n)+1.d-40),0.d0)
cc          trm(i,j,l,n) = 0.d0 !trm=0 here avoids possible roundoff
          else
cc          scalmom = max(1.d0-rloss/(trm(i,j,l,n)+1.d-40),0.d0)
cc          trm(i,j,l,n) = trm(i,j,l,n) + (rprod-rloss)
          end if
          tr3Dsource(i,j,l,nsp,n) = rprod/dtsrc
          tr3Dsource(i,j,l,nsl,n) = -rloss/dtsrc
        enddo
        enddo
        enddo
      RETURN
      END SUBROUTINE Trop_chem_O3


      SUBROUTINE linoz_depo(ns,n)
C****
C**** Deposition from layer 1
C**** Deposition Velocity is in cm/sec.  Convert to kg
C****
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds
      USE LINOZ_CHEM_COM, only: O3_DepVel
      USE RESOLUTION, only: im,jm
      USE MODEL_COM, only: modelEclock,itime,dtsrc
      USE ATM_COM, only: t,pmid,pk,pdsig
      USE TRACER_COM
      USE CONSTANT, only : grav,rgas
      USE GEOM, only: imaxj
      USE QUSDEF, only : mz,mzz
      USE FLUXES, only: trsource
      implicit none
      integer i,j,l,n,ns,jmon
      INTEGER :: J_1, J_0, I_0, I_1
      real*8 tmsurf,dmass,tk
      real*8 dz ! = -dP/rhoG; rho=PRT; Deposition velocity

C**** Extract useful local domain parameters from "grid"
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP
      jmon = modelEclock%getMonth()

      l=1
      do j=j_0,j_1
        do i=i_0,imaxj(j)
c         write(*,*)' pdsig',l,i,j,pdsig(l,i,j)
          tk = t(i,j,l)*pk(l,i,j)            ! Temp in kelvin
          dz = pdsig(l,i,j)*rgas*tk/(pmid(l,i,j)*grav)   ! meters
          tmsurf = max(0.d0,trm(i,j,l,n)
     *                -trmom(mz,I,J,L,n)+trmom(mzz,I,J,L,n))
          dmass = -O3_DepVel(i,j,jmon)*tmsurf*dtsrc/(dz*100.)
          if(trm(i,j,l,n) + dmass .lt.0.) then
            write(6,'(a,3i3,2f12.0,3E12.3,i9,e12.3,2f7.1)')
     *      ' Negative O3 from deposition', i,j,l,trm(i,j,l,n),
     *      dmass,O3_DepVel(i,j,jmon),dz,tmsurf,itime
     *      ,pdsig(l,i,j),tk,pmid(l,i,j)
            dmass = -trm(i,j,l,n)
c           trm(i,j,l,n) = 0.d0
          else
c           trm(i,j,l,n) = trm(i,j,l,n) + dmass
          end if
          trsource(i,j,ns,n) = dmass/dtsrc
        enddo
      enddo
      RETURN
      END


      SUBROUTINE Strat_chem_O3(ns,n)
!@vers 2013/03/26
c-----------------------------------------------------------------------
c   Strat_chem_O3 applies linearized chemistry based on tables from
c    PRATMO model using climatological T, O3, time of year
c-----------------------------------------------------------------------
c  stratospheric chem occurs in top NSTRTC layers of CTM
c  TLT0M(J,LR,N) is stored LR from top (=LM) down (=LM+1-NCSTRT)
c
c Stratospheric Chemistry Tables for O3:
c ======================================
c   7 tables, each a function of month (12), latitude
c   (18, -85 to 85 in 10 deg. increments) and altitude
c   (25, z*=10-58 km in 2 km increments).
c  1- ozone (Logan climatology), v/v
c  2- Temperature climatology, K
c  3- Column ozone climatology, Logan ozone integrated above box, DU
c  4- ozone (P-L) for climatological ozone, v/v/s
c  5- d(P-L) / dO3, 1/s
c  6- d(P-L) / dT, v/v/s/K
c  7- d(P-L) / d(column O3), v/v/s/DU
cXX8- d(P-L) / d(sol.flx.) (-)   <<<XXX not in this version>>>
cXXXXX DSOL NOT USED XXXXX
c!@var dsol variable describes portion of solar cycle being modeled
c!@+ +1.0 = solar max, 0.0 = neutral, -1.0 = solar min
c!@+ can go beyond +1,-1 for more extreme situations
c!@+ (but remember this is a linear approximation)
c!@+ use dsol=0.0 for 'standard linoz' runs
cXXXXX DSOL NOT USED XXXXX

      USE CONSTANT, only : avog
      USE RESOLUTION, only: im,jm,lm
      USE MODEL_COM, only: itime,dtsrc
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds
      USE ATM_COM, only: t,pk,MA! Air mass of each box (kg/m^2)
      USE GEOM, only: imaxj,axyp
      USE TRACER_COM
      USE PRATHER_CHEM_COM, only: nstrtc
      USE LINOZ_CHEM_COM, only: tlT0M,TLTZM,TLTZZM,dsol
      USE FLUXES, only: tr3Dsource
      implicit none
      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm) ::
     &     dcolo3,colo3
      real*8
     &  dero3,scalmom,pmltot,dertmp,dtmp,derco3,dco3,sso3,
     &  climo3,climpml,dersol
      real*8 dmass,T0Mold
      integer i,j,l,lr,n,ns,najl   ,kx
      INTEGER :: J_1, J_0, I_0, I_1
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

cc      najl = jls_3Dsource(ns,n)
c start at top layer and continue to lowest layer for strat. chem
      DO 330 l = lm,lm+1-nstrtc,-1
        LR = LM+1-L
        DO 320 J=J_0,J_1
            if (tlT0M(j,lr,5) == 0.) go to 320
          do 310 i=I_0,imaxj(j)
            if (trm(i,j,l,n).le.0.d0) goto 310

c calculate ozone column above box (and save)
c   dcolo3 = ozone column (in DU) in given layer
c   colo3 =  ozone column above layer + half of column in layer
            if (l.eq.lm) then           !top model layer
              dcolo3(i,j,l) = trm(i,j,l,n) / axyp(i,j) *
     &          avog/(tr_mm(n)*1d-3)/ 2.687d16 * 1d-4
              colo3(i,j,l) = dcolo3(i,j,l)*0.5
            else
              dcolo3(i,j,l) = trm(i,j,l,n)/ axyp(I,J) *
     &          avog/(tr_mm(n)*1d-3)/ 2.687d16 * 1d-4
              colo3(i,j,l) = colo3(i,j,l+1) +
     &          (dcolo3(i,j,l)+dcolo3(i,j,l+1))*0.5
            endif

c ****** O3 Chemistry  ******
c store tracer mass before chemistry
            T0Mold=trm(i,j,l,n)
c climatological P-L:
            climpml = tlT0M(j,lr,4)/mass2vol(n)*MA(l,i,j)*axyp(i,j)
c local ozone feedback:
            dero3=tlT0M(j,lr,5)
            climo3 = tlT0M(j,lr,1)/mass2vol(n)*MA(l,i,j)*axyp(i,j)
c column ozone feedback:
            derco3 = tlT0M(j,lr,7)/mass2vol(n)*MA(l,i,j)*axyp(i,j)
            dco3=(colo3(i,j,l)-tlT0M(j,lr,3))
c temperature feedback: T is potential temp, need to convert
            dertmp = tlT0M(j,lr,6)/mass2vol(n)*MA(l,i,j)*axyp(i,j)
            dtmp=(t(i,j,l)*PK(L,I,J)-tlT0M(j,lr,2))
c define sol.flux. derivative and convert from mixing ratio to mass
CXXX        dersol = tlT0M(j,lr,8)/mass2vol(n)*MA(l,i,j)*axyp(i,j)
c calulate steady-state ozone:
            sso3=climo3 - (climpml+dco3*derco3+dtmp*dertmp)/dero3
CXXX        sso3=climo3 -
CXXX *        (climpml+dco3*derco3+dtmp*dertmp+dsol*dersol)/dero3
c change in ozone mass due to chemistry:
            dmass=(sso3-T0Mold)*(1.0-exp(dero3*dtsrc))
c update ozone mass
            if (T0Mold+dmass.lt.0.) then
               write(6,'(a,3i4,2f20.2,i9)')
     *           ' Negative tracer in Strat_chem_O3',
     *                i,j,l,T0Mold,dmass,itime
               dmass = -T0Mold
cc               trm(i,j,l,n) = 0.d0
cc            else
cc               trm(i,j,l,n) = T0Mold + dmass
            end if
            tr3Dsource(i,j,l,ns,n) = dmass/dtsrc
cc            call inc_tajls(i,j,l,najl,dmass)
c scale moments by fractional change in total tracer mass
cc        if (dmass.lt.0.d0) then
cc          scalmom = trm(i,j,l,n)/T0Mold
cc          trmom(1:nmom,I,J,L,n) = trmom(1:nmom,I,J,L,n) * scalmom
cc        end if
  310 continue
  320 continue
  330 continue
      return
      end SUBROUTINE Strat_chem_O3


c-----------------------------------------------------------------------
      SUBROUTINE linoz_STRATL
c-----------------------------------------------------------------------
c-------- monthly fixup of chemistry PARAM'S
c
      USE RESOLUTION, only: jm,lm
      USE MODEL_COM, only: modelEclock
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds
      USE PRATHER_CHEM_COM, only: jlatmd,p0l,NSTRTC
      USE LINOZ_CHEM_COM, only: nctable,TLPARM,
     *    tlt0m,tltzm,tltzzm,lz_linoz,lz_lx,ps
      implicit none
      real*8  STRT0L(LM),STRT1L(LM),STRT2L(LM),STRTX(lz_linoz)
      real*8 f(lz_lx)
      integer j,jj,k,lr,n,jmon

      INTEGER :: J_1, J_0
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      jmon = modelEclock%getMonth()

c-------- TLPARM(25,18,12,N) defined for -----------------------------
c lz_linoz  25 layers from 58 km to 10 km by 2 km intervals
c            18 LATS (85S, 75S, ...85N)
c            12 months
c            N tables = NCTABLE
c-------- skip interpolating, pick nearest latitude --------------------

      DO N = 1,NCTABLE
      DO J = J_0,J_1
        JJ = JLATMD(J)
        DO K = 1,lz_linoz
          STRTX(K) = TLPARM(K,JJ,jmon,N)
        ENDDO

c-------- stratospheric chem occurs in top NSTRTC layers ---------------
c
        CALL STRT2M(STRTX,lz_linoz,STRT0L,STRT1L,STRT2L,P0L,NSTRTC
     *      ,ps,f,lz_lx)
c
c-------store loss freq/yields & moments in TLT0M/TLTZM/TLTZZM
c-------             for exact CTM layers LM down
          DO LR = 1,NSTRTC
            TLT0M(J,LR,N) = STRT0L(LR)
            TLTZM(J,LR,N) = STRT1L(LR)
            TLTZZM(J,LR,N) = STRT2L(LR)
          ENDDO
        ENDDO   ! J
      ENDDO     ! N
      return
      END SUBROUTINE linoz_STRATL


      SUBROUTINE STRT2M (STRTX,NX,STRT0L,STRT1L,STRT2L,P0L,NSTRT,
     *  ps,f,nl)
c-----------------------------------------------------------------------
      implicit none
      integer  NX,NSTRT,nl,l,k
      real*8   P0L(*),STRT0L(*),STRT1L(*),STRT2L(*),STRTX(*)
      real*8   P1,P2,F0,F1,F2,PS(*),F(*)
c-----------------------------------------------------------------------
c       set up std z* atmosphere: p = 1000 * 10**(-z*/16 km)
c       assume that stratospheric chemical parameters always start at
cc       52 km (N=27) scan downward from 52 km to 14 km (NX=20) by 2 km
c       58 km (N=30) scan downward from 58 km to 10 km (NX=25) by 2 km
c       intervals, constant >58km
c-------- N.B. F(@30km) assumed to be constant from 29-31 km (by mass)
c
      DO L = 1,NL-NX
        F(L)     = 0.D0
      ENDDO
      DO K = 1,NX    ! K=1 is at the top of atmosphere
        F(NL+1-K)= STRTX(K)
      ENDDO
      DO K = 1,NSTRT
        P1       = P0L(K+1)
        P2       = P0L(K)
        CALL SOMLFQ(P1,P2,F0,F1,F2,PS,F,NL)
        STRT0L(K)= F0
        STRT1L(K)= F1
        STRT2L(K)= F2
      ENDDO
      RETURN
      END SUBROUTINE STRT2M


      SUBROUTINE SOMLFQ(P1,P2,F0,F1,F2,PS,F,NL)
c-----------------------------------------------------------------------
c--- calculate loss freq moments from a set of loss freq's at std z*
c--  given a CTM model interval pressure range: P1 > P2 (decreasing up)
c-----  the pressure levels BETWEEN z* values are:
c                      PS(i) > PS(i+1) bounds z*(i)
c----- NL:  z* levels, ==> PS(NL+1) = 0
c-----            (extrapolate chemical loss to top)
c      Z1 = 16.D0*LOG10(1000.D0/P1)
c      Z2 = 16.D0*LOG10(1000.D0/P2)
c
c----- The MOMENTS for a square-wave or 'bar':
c-----             F(x)=f0  b<=x<=c, =0.0 else
c---    S0 =   f0 (x)                      [from x=b to x=c]
c---    S1 = 3 f0 (x^2 - x)                [from x=b to x=c]
c---    S2 = 5 f0 (2x^3 - 3x^2 + x)        [from x=b to x=c]
c-----------------------------------------------------------------------
      USE CONSTANT, only: by3
      implicit none
      integer  NL,I
      real*8   P1,P2,F0,F1,F2,PS(NL+1),F(NL),sgnf0
      real*8   XB,XC,PC,PB
c-----------------------------------------------------------------------
      F0     = 0.D0
      F1     = 0.D0
      F2     = 0.D0
      DO I = 1,NL
        PC   = MIN(P1,PS(I))
        PB   = MAX(P2,PS(I+1))
        IF (PC .GT. PB)  THEN
C------ have condition:  P1>=PC > PB>=P2, 0<=XB < XC<=1 --------------
          XC = (PC-P2)/(P1-P2)
          XB = (PB-P2)/(P1-P2)
c
c------ assume that the loss freq, F, is constant over interval [XB,XC],
c------ F0: (c-b),  F1: 6((c2-c)-(b2-b)), F2: 5((2c3-3c2+c)-(2b3-3b2+b))
c------ calculate its contribution to the moments in the interval [0,1]
c
          F0 = F0 +F(I) *(XC-XB)
          F1 = F1 +F(I) *3.D0 *((XC*XC-XC) - (XB*XB-XB))
          F2 = F2 +F(I) *5.D0 *
     &         ((XC+XC-1.D0)*(XC*XC-XC) - (XB+XB-1.D0)*(XB*XB-XB))
        ENDIF
      ENDDO
c
c-------- RESTRAIN moments: force monotonicity & positive at min end pt
c
c -=-=- cam: tables can be + or -
      if (f0.ne.0.0) then
        sgnf0=f0 / abs(f0)
      else
        sgnf0=1.0
      endif
      f0=abs(f0)

      IF (F2 .GT. 0.D0)  THEN
c
c-------- do not allow reversal of curvature: F2 > 0 -------------------
        F2   = MIN(F2, ABS(F1)*by3, 5.D-1*F0)
        IF (F1 .LT .0.D0)  THEN
          F1 = MAX(-(F0+F2), F1)
         ELSE
          F1 = MIN(+(F0+F2), F1)
        ENDIF
      ELSE
c
c-------- F2 < 0 = curved down at ends, allow if F1 < F0 ---------------
        F1  = MIN(F0,MAX(-F0,F1))
        F2  = MAX(F2,(ABS(F1)-F0),(-ABS(F1)*by3))
      ENDIF
c
c -=-=- cae: apply sign
      f0=sgnf0 * f0
      f1=sgnf0 * f1
      f2=sgnf0 * f2
      RETURN
      END SUBROUTINE SOMLFQ


      subroutine read_CH4_sources(nt)
!@sum reads in CH4 sources and sinks
!@auth Jean Lerner
C****
C**** There are 3 monthly sources and 11 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE RESOLUTION, only: im,jm
      USE MODEL_COM, only: itime,modelEclock
      use TimeConstants_mod, only: SECONDS_PER_DAY, INT_DAYS_PER_YEAR
      use TimeConstants_mod, only: JDPERY
      USE FLUXES, only: focean,fearth0,flake0
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds,
     *  readt_parallel, AM_I_ROOT
      use OldTracer_mod, only: itime_tr0,trname
      USE FILEMANAGER, only: openunit,closeunit
      USE FILEMANAGER, only: nameunit
      USE CH4_SOURCES, only: src=>ch4_src,nsrc=>nch4src
      implicit none
      character*80 title
!@var adj Factors that tune the total amount of individual sources
      real*8 adj(nsrc)
      data adj/1.3847,1.0285,3.904,1.659,1.233,1.194,0.999,
!    *  7.2154, 3.7247d-5,3.1399d-4,5.4838d-5,   !Model II prime
     *  7.2154, 3.5997d-5,17.330d-4,5.3558d-5,
     *  0.4369,0.7533,0.9818/
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=11,nmons=3
      integer ann_units(nanns-3),mon_units(nmons)
      character*12 :: ann_files(nanns-3) =
     *  (/'CH4_ANIMALS ','CH4_COALMINE','CH4_GASLEAK ','CH4_GASVENT ',
     *    'CH4_CITYDUMP','CH4_SOIL_ABS','CH4_TERMITES','CH4_COALBURN'/)
      logical :: ann_bins=.true.
      character*8 :: mon_files(nmons) =
     *   (/'CH4_BURN','CH4_RICE','CH4_WETL'/)
      real*8 adj_wet(GRID%J_STRT_HALO:GRID%J_STOP_HALO)
      integer :: kwet=14         !!! position of wetlands array in src
      logical :: mon_bins=.true.

c GISS-ESMF EXCEPTIONAL CASE - SAVE variable, I/O
      real*8, allocatable :: tlca(:,:,:), tlcb(:,:,:)   ! for monthly sources
      real*8 frac
      integer i,j,nt,iu,k,imon(nmons)
      logical :: ifirst=.true.
      integer :: jdlast=0
      integer :: jday
      save ifirst,jdlast,tlca,tlcb,mon_units,imon
      INTEGER :: J_1, J_0, J_0H, J_1H, I_0H, I_1H

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      jday = modelEclock%getDayOfYear()

      adj_wet = 0
      do J=J_0,J_1
         if ( (J>JM/3) .AND. (J<=JM-JM/3) ) then   !!! low latitudes
            adj_wet(J) = 1.761
         else
            adj_wet(J) = 0.6585                  !!! high latitudes
         endif
      end do

      if (itime.lt.itime_tr0(nt)) return
C****
C**** Annual Sources and sinks
C**** Apply adjustment factors to bring sources into balance
C**** Annual sources are in KG C/M2/Y
C**** Sources need to be kg/m^2 s; convert /year to /s
C****
      if (ifirst) then
         call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
         I_0H=GRID%I_STRT_HALO
         I_1H=GRID%I_STOP_HALO
         Allocate(tlca(i_0H:i_1H,j_0H:j_1H,nmons))
         Allocate(tlcb(i_0H:i_1H,j_0H:j_1H,nmons))
        k = 0
        call openunit(ann_files,ann_units,ann_bins)
        do iu = 1,nanns-3
          k = k+1
          call readt_parallel (grid,
     &         ann_units(iu),nameunit(ann_units(iu)),src(:,:,k),1)
!ref      call readt          (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,:,k) = src(:,:,k)*adj(k)/
     &                 (SECONDS_PER_DAY*INT_DAYS_PER_YEAR)
        end do
        call closeunit(ann_units)
        ! 3 miscellaneous sources
          k = k+1
          src(:,:,k) = focean(:,:)*adj(k)/
     &                 (SECONDS_PER_DAY*INT_DAYS_PER_YEAR)
          k = k+1
          src(:,:,k) = flake0(:,:)*adj(k)/
     &                 (SECONDS_PER_DAY*INT_DAYS_PER_YEAR)
          k = k+1
          src(:,:,k) = fearth0(:,:)*adj(k)/
     &                 (SECONDS_PER_DAY*INT_DAYS_PER_YEAR)

      call openunit(mon_files,mon_units,mon_bins)
      endif
C****
C**** Monthly sources are interpolated to the current day
C****
C**** Also, Apply adjustment factors to bring sources into balance
C**** Monthly sources are in KG C/M2/S => src in kg/m^2 s
C****
      ifirst = .false.
      j = 0
      do k = nanns+1,nsrc
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(:,:,j),tlcb(:,:,j),src(:,:,k),frac,imon(j))
        src(:,J_0:J_1,k) = src(:,J_0:J_1,k)*adj(k)
      end do
      jdlast = jday
      if (AM_I_ROOT())
     *  write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
C**** Zonal adjustment for combined wetlands and tundra
C****
      do j=J_0,J_1
        src(:,j,kwet) = src(:,j,kwet)*adj_wet(j)
      end do
      return
      end subroutine read_CH4_sources


      MODULE CO2_SOURCES
      USE TRACER_COM
!@var co2_src C02 surface sources and sinks (kg/s)
      integer, parameter :: nco2src=6
      real*8, ALLOCATABLE, DIMENSION(:,:,:) :: co2_src
      END MODULE CO2_SOURCES

      subroutine read_CO2_sources(nt)
!@sum reads in CO2 sources and sinks
!@auth Jean Lerner
C****
C**** There are two monthly sources and 4 annual sources
C**** Annual sources are read in at start and re-start of run only
C**** Monthly sources are interpolated each day
      USE RESOLUTION, only: im,jm
      USE MODEL_COM, only: itime,modelEclock
      use TimeConstants_mod, only: JDPERY
      USE DOMAIN_DECOMP_ATM, only : grid, getDomainBounds, AM_I_ROOT
      use TimeConstants_mod, only: SECONDS_PER_DAY, INT_DAYS_PER_YEAR
      USE DOMAIN_DECOMP_ATM, only : READT_PARALLEL
      use OldTracer_mod, only: itime_tr0,trname
      USE CO2_SOURCES, only: src=>co2_src,nsrc=>nco2src
      USE FILEMANAGER, only: openunit,closeunit
      USE FILEMANAGER, only: nameunit
      implicit none
      character*80 title
!@var adj Factors that tune the total amount of individual sources
      real*8 adj(nsrc)
      data adj/3.81d0,3.67d0,3.67d0,19.54d0,3.67d0,6.42d0/
!@var nanns,nmons: number of annual and monthly input files
      integer, parameter :: nanns=4,nmons=2
      integer ann_units(nanns),mon_units(nmons)
      character*12 :: ann_files(nanns) =
     *  (/'CO2_FOS_FUEL','CO2_FERT    ','CO2_REGROWTH','CO2_LAND_USE'/)
      logical :: ann_bins=.true.
      character*9 :: mon_files(nmons) = (/'CO2_VEG  ','CO2_OCEAN'/)
      logical :: mon_bins=.true.

c GISS-ESMF EXCEPTIONAL CASE - SAVE and I/O issues
      real*8, Allocatable, DIMENSION(:,:,:) :: tlca, tlcb ! for monthly sources
      real*8 frac
      integer i,j,nt,iu,k,imon(nmons)
      logical :: ifirst=.true.
      integer :: jdlast=0
      save ifirst,jdlast,tlca,tlcb,mon_units,imon
      integer :: J_0, J_1, J_0H, J_1H, I_0H, I_1H
      integer :: jday

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      jday = modelEclock%getDayOfYear()

      if (itime.lt.itime_tr0(nt)) return
C****
C**** Annual Sources and sink
C**** Apply adjustment factors to bring sources into balance
C**** Annual sources are in KG C/M2/Y
C**** Sources need to be kg/m^2 s; convert /year to /s
C****
      if (ifirst) then
        I_0H=GRID%I_STRT_HALO
        I_1H=GRID%I_STOP_HALO
        J_0H=GRID%J_STRT_HALO
        J_1H=GRID%J_STOP_HALO
        Allocate(tlca(i_0H:i_1H,j_0H:j_1H,nmons))
        Allocate(tlcb(i_0H:i_1H,j_0H:j_1H,nmons))
        call openunit(ann_files,ann_units,ann_bins)
        k = 0
        do iu = 1,nanns
          k = k+1
          call readt_parallel (grid,
     &         ann_units(iu),nameunit(ann_units(iu)),src(:,:,k),1)
!ref      call readt (iu,0,src(1,1,k),im*jm,src(1,1,k),1)
          src(:,J_0:J_1,k) = src(:,J_0:J_1,k)*adj(k)/
     &                       (SECONDS_PER_DAY*INT_DAYS_PER_YEAR)
        end do
        call closeunit(ann_units)

      call openunit(mon_files,mon_units,mon_bins)
      endif
C****
C**** Monthly sources are interpolated to the current day
C****
C**** Also, Apply adjustment factors to bring sources into balance
C**** Monthly sources are in KG C/M2/S => src in kg/m^2 s
      ifirst = .false.
      j = 0
      do k=nanns+1,nsrc
        j = j+1
        call read_monthly_sources(mon_units(j),jdlast,
     *    tlca(:,:,j),tlcb(:,:,j),src(:,:,k),frac,imon(j))
        src(:,J_0:J_1,k) = src(:,J_0:J_1,k)*adj(k)
      end do
      jdlast = jday
      if (AM_I_ROOT())
     *  write(6,*) trname(nt),'Sources interpolated to current day',frac
      call sys_flush(6)
C****
      return
      end subroutine read_CO2_sources


      SUBROUTINE get_14CO2_IC(CO2IJL)
!@sum GET_14CO2_IC Calculates initial distribution for 14CO2 tracer
!@auth J.Lerner (modified from program by G. Russell)
C**** NOTE: tracer is supposed to start on 10/16
C**** October 1963 14CO2 Concentrations for GCM  2/26/99
C**** 2/2/2: generalized code for modelE
C****
      USE RESOLUTION, ONLY: im,jm,lm
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds
      USE ATM_COM, only: pedn
      USE FILEMANAGER, only: openunit,closeunit
      USE GEOM, only : DLAT_DG
      IMPLICIT NONE
      integer,PARAMETER :: kmwco2=60
      REAL*4 CO2W(37,0:30)
      real*8 p(0:60)
      real*8 CO2JK(GRID%J_STRT_HALO:GRID%J_STOP_HALO,0:kmwco2)
      real*8 CO2IJL(GRID%I_STRT:GRID%I_STOP,
     &              GRID%J_STRT:GRID%J_STOP,LM)
      CHARACTER*80 TITLE
      integer i,j,jw,k,l,n,iu_in,iu_out
      Real*8 :: pup,cup,pdn,cdn,psum,csum,w,zk !,stratm

      INTEGER :: I_0, I_1, J_1, J_0
      INTEGER :: J_0S, J_1S
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT     =J_0,    J_STOP     =J_1,
     &               J_STRT_SKP =J_0S,   J_STOP_SKP =J_1S,
     &               HAVE_SOUTH_POLE = HAVE_SOUTH_POLE,
     &               HAVE_NORTH_POLE = HAVE_NORTH_POLE)
      I_0 = GRID%I_STRT
      I_1 = GRID%I_STOP

C****
C**** Read in CO2 concentrations from workshop
C****
c     OPEN (1,FILE='workshop.14co2',STATUS='OLD')
      call openunit('14CO2_IC_DATA',iu_in,.false.)
      DO 110 N=1,10
  110 READ (iu_in,911)
      READ (iu_in,911) ((CO2W(J,K),J=1,13),K=0,30)
      READ (iu_in,911)
      READ (iu_in,911)
      READ (iu_in,912) ((CO2W(J,K),J=14,25),K=0,30)
      READ (iu_in,912)
      READ (iu_in,912)
      READ (iu_in,912) ((CO2W(J,K),J=26,37),K=0,30)
      call closeunit(iu_in)
C**** Calculate workshop pressure levels (Pa)
      DO K=0,kmwco2
        ZK = 2.*K
        P(K) = 100000.*10**(-ZK/16.)
      end do
C****
C**** Interpolate workshop CO2W to CO2JK on GCM latitudes
C****
      DO K=31,kmwco2    ! above 17 Pa set all equal
        CO2JK(:,K) = 250.
      end do
      DO K=0,30
        IF (HAVE_SOUTH_POLE) CO2JK( 1,K) = CO2W( 1,K)
        IF (HAVE_NORTH_POLE) CO2JK(JM,K) = CO2W(37,K)
        DO J=J_0S,J_1S
          W = 1. + (J-1)*0.2*DLAT_DG
          JW=W
          CO2JK(J,K) = CO2W(JW,K)*(JW+1-W) + CO2W(JW+1,K)*(W-JW)
      end do; end do
C****
C**** Interpoate CO2J to CO2IJL on GCM grid conserving vertical means
C****
C**** P(K), PUP, PDN, PSUM in (Pa) = (mb*100)
      DO 440 J=J_0,J_1
      DO 440 I=I_0,I_1
      PDN = pedn(1,i,j)*100.
      CDN = CO2JK(J,0)
      K=1
      DO 430 L=1,LM
      PSUM = 0.
      CSUM = 0.
      PUP = PEDN(L+1,I,J)*100
  410 IF(P(K).LE.PUP)  GO TO 420
      PSUM = PSUM +  PDN-P(K)
      CSUM = CSUM + (PDN-P(K))*(CDN+CO2JK(J,K))/2.
      PDN  = P(K)
      CDN  = CO2JK(J,K)
      K=K+1
      if (k.gt.kmwco2) call stop_model(
     &     ' Please increase kmwco2 in get_14CO2_IC',255)
      GO TO 410
C****
  420 CUP  = CO2JK(J,K) + (CO2JK(J,K-1)-CO2JK(J,K))*(PUP-P(K))/
     /       (P(K-1)-P(K))
      PSUM = PSUM +  PDN-PUP
      CSUM = CSUM + (PDN-PUP)*(CDN+CUP)/2.
      CO2IJL(I,J,L) = CSUM/PSUM
      PDN = PUP
  430 CDN = CUP
  440 CONTINUE
C**** Scale data to proper units (10**-18 kg 14CO2/kg air)
      CO2IJL(:,:,:) = CO2IJL(:,:,:)*4.82d0*(14.+16.*2.)/29.029d0
      RETURN
C****
  911 FORMAT (5X,13F5.0)
  912 FORMAT (5X,12F5.0)
      END SUBROUTINE get_14CO2_IC


      module lhntr_com
!@sum LHNTR routines from LHNTR.GRD for Linear Horizontal Interpolation
!@auth G. Russell, modified by J. Lerner for modelE
      IMPLICIT NONE
      REAL*8 OFFIA,DIVJA,OFFIB,DIVJB,SKIB,SKIP
      REAL*8 WWEST(720),WSOUTH(361)
      INTEGER*4 IWEST(720),IEAST(720),JSOUTH(361),JNORTH(361)
      LOGICAL*4 QMPOLE
      contains

C**** LHNTR.GRD    Linear Horizontal Interpolation REAL*8    1/25/93
C****
      SUBROUTINE LHNTR0 (IMA,JMA,OFFIA,DIVJA,
     *                   IMB,JMB,OFFIB,DIVJB,SKIB)
C****
C**** LHNTR performs a linear horizontal interpolation of per unit
C**** mass or per unit area quantities defined on grid A, calculating
C**** the quantity on grid B.  B grid values, that cannot be
C**** calculated because the four surrounding A grid points have
C**** weight zero, are set to the value of SKIP.  This scheme does
C**** not conserve the area weighted integral of the quantity.
C****
C**** The center of grid point values are calculated as follows:
C****   Longitude = -180 + 360*(I-.5+OFFI)/IM
C****   Latitude  =  -90 + 180*(J-.5+OFFJ)/DIVJ
C****
C****
      INTEGER :: IMA,JMA,IMB,JMB,IB,JB
      REAL*8 :: OFFIA,DIVJA,OFFIB,DIVJB,SKIB,RIMA,RIA,RJA,FJAEQ,FJBEQ

      SKIP = SKIB
      IF(IMB.LT.1 .OR. IMB.GT.720 .OR.
     *   JMB.LT.1 .OR. JMB.GT.361)  GO TO 300
C****
C**** Determine location and weighting of B grid points with
C**** respect to A grid points in the West-East direction.
C**** IWEST(IB) = A grid point just West of B grid point
C**** IEAST(IB) = A grid point just East of B grid point
C**** WWEST(IB) = weight of A grid point IWEST(IB)
C****
      RIMA = IMA
      DO 110 IB=1,IMB
      RIA = DMOD(.5D0-OFFIA+(IB-.5D0+OFFIB)*IMA/IMB+2.*RIMA,RIMA)
      IWEST(IB) = RIA
      IEAST(IB) = IWEST(IB)+1
      IF(IWEST(IB).LT.1)  IWEST(IB) = IMA
  110 WWEST(IB) = IEAST(IB)-RIA
C****
C**** Determine location and weighting of B grid points with
C**** respect to A grid points in the South-North direction.
C**** JSOUTH(JB) = A grid point just South of B grid point
C**** JNORTH(JB) = A grid point just North of B grid point
C**** WSOUTH(JB) = weight of A grid point JSOUTH(JB)
C****
      FJAEQ = (1+JMA)/2.
      FJBEQ = (1+JMB)/2.
      DO 230 JB=1,JMB
      RJA = FJAEQ + (JB-FJBEQ)*DIVJA/DIVJB
      IF(RJA.LT.1.)   GO TO 210
      IF(RJA.GE.JMA)  GO TO 220
      JSOUTH(JB) = RJA
      JNORTH(JB) = JSOUTH(JB)+1
      WSOUTH(JB) = JNORTH(JB)-RJA
      GO TO 230
C**** B grid point is South of most Southward A grid point
  210 JSOUTH(JB) = 1
      JNORTH(JB) = 2
      WSOUTH(JB) = 1.
      GO TO 230
C**** B grid point is North of most Northward A grid point
  220 JSOUTH(JB) = JMA-1
      JNORTH(JB) = JMA
      WSOUTH(JB) = 0.
  230 CONTINUE
      RETURN
C****
C**** Invalid arguments or B dimensions are out of range
C****
  300 WRITE (6,930) IMA,JMA,OFFIA,DIVJA,
     *              IMB,JMB,OFFIB,DIVJB,SKIP
      call stop_model('LHNTR: 300',255)
  930 FORMAT ('0Arguments received by LHNTR0 in order:'/
     *   2I12,' = IMA,JMA = array dimensions for A grid'/
     *  E24.8,' = OFFIA   = fractional number of grid boxes from',
     *                    ' IDL to left edge of grid box I=1'/
     *  E24.8,' = DIVJA   = number of whole grid boxes from SP to NP'/
     *   2I12,' = IMB,JMB = array dimensions for B grid'/
     *  E24.8,' = OFFIB   = fractional number of grid boxes from',
     *                    ' IDL to left edge of grid box I=1'/
     *  E24.8,' = DIVJB   = number of whole grid boxes from SP to NP'/
     *  E24.8,' = SKIP    = value to be put in B array when B',
     *  ' grid box is subset of A grid boxes with WTA = 0'/
     *  '0These arguments are invalid or out of range.')
      end SUBROUTINE LHNTR0

      subroutine LHNTR (WTA,A,B,IMA,JMA,IMB,JMB)
C****
C**** LHNTR performs the linear interpolation
C**** Input: WTA = weighting array for values on the A grid
C****          A = per unit mass or per unit area quantity
C**** Output:  B = linearly interpolated quantity on B grid
C****
      INTEGER :: IMA,JMA,IMB,JMB,JB,JS,JN,IB,IW,IE
      REAL*8 :: WEIGHT,VALUE,BMEAN
      REAL*8 WTA(IMA,JMA),A(IMA,JMA),B(IMB,JMB)
      QMPOLE = .FALSE.
      DO 510 JB=1,JMB
      JS = JSOUTH(JB)
      JN = JNORTH(JB)
      DO 510 IB=1,IMB
      IW = IWEST(IB)
      IE = IEAST(IB)
      B(IB,JB) = SKIP
      WEIGHT = (WTA(IW,JS)*    WWEST(IB)
     *        + WTA(IE,JS)*(1.-WWEST(IB)))*    WSOUTH(JB)
     *       + (WTA(IW,JN)*    WWEST(IB)
     *        + WTA(IE,JN)*(1.-WWEST(IB)))*(1.-WSOUTH(JB))
      IF(WEIGHT.EQ.0.)  GO TO 510
      VALUE = (WTA(IW,JS)*A(IW,JS)*    WWEST(IB)
     *       + WTA(IE,JS)*A(IE,JS)*(1.-WWEST(IB)))*    WSOUTH(JB)
     *      + (WTA(IW,JN)*A(IW,JN)*    WWEST(IB)
     *       + WTA(IE,JN)*A(IE,JN)*(1.-WWEST(IB)))*(1.-WSOUTH(JB))
      B(IB,JB) = VALUE/WEIGHT
  510 CONTINUE
C****
C**** Replace individual values near the poles by longitudinal mean
C****
      IF(.NOT.QMPOLE)  RETURN
      DO 630 JB=1,JMB,JMB-1
      WEIGHT = 0.
      VALUE  = 0.
      DO 610 IB=1,IMB
      IF(B(IB,JB).EQ.SKIP)  GO TO 610
      WEIGHT = WEIGHT + 1.
      VALUE  = VALUE  + B(IB,JB)
  610 CONTINUE
      BMEAN  = SKIP
      IF(WEIGHT.NE.0.)  BMEAN = VALUE/WEIGHT
      DO 620 IB=1,IMB
  620 B(IB,JB) = BMEAN
  630 CONTINUE
      RETURN
      end subroutine LHNTR
      end module lhntr_com


      subroutine get_Trop_chem_CH4_freq(FRQname)
!@sum get_Trop_chem_CH4_freq interpolates troposphereic chemical
!@+     rates for CH4 from n-grid, 9 layers to 4X5, lm layers
!@+     The resulting rate file is written to disk for use throughout
!@+     the run.  The rates are changed every 5 days.
!@auth Jean Lerner
C****  Input: CLIM.RUN.OHCH4.FRQ
C**** Output: temporary file is for this vertical resolution
C**** WARNING: RESULTS ARE INTENDED FOR USE TO ABOUT 26.5 mb ONLY
C**** However, we'll stop at LMTC, which is lower, and let strat
C****   chem pick up from there.
      USE RESOLUTION, only: im,jm,lm
      USE RESOLUTION, only: psf
      USE ATM_COM, only: pmidl00
      USE GEOM, only : DLAT_DG
      USE FILEMANAGER, only: openunit,closeunit
      USE lhntr_com, only: LHNTR,LHNTR0
      implicit none
      integer l,i,j,km,imo,jmo,lmo,kmo,InFile,interp_file,ltopx,it
      parameter (km=im*jm, imo=36,jmo=24,lmo=9,kmo=imo*jmo)
      character*80 title
      character*16 :: FRQname
      logical :: debug=.true.
c     real*4 fold(kmo,lmo)   ,rlat(jm)
c     real*8 wta(kmo),foldlm(kmo,lm),
c    *  fnew(km,lm),pold(lmo),pnew(lm),ain(lmo),aout(lm)
      real*4 tau
      real*4, ALLOCATABLE, DIMENSION(:,:) :: fold
      real*8, ALLOCATABLE, DIMENSION(:,:) :: foldlm,fnew
      real*8, ALLOCATABLE, DIMENSION(:) :: wta
      real*8 pold(lmo),pnew(lm),ain(lmo),aout(lm),divj
      real*8 :: sigo(lmo) = (/.974264d0,.907372d0,.796957d0,.640124d0,
     *    .470418d0,.318899d0,.195759d0,.094938d0,.016897d0/)

      ALLOCATE (fold(kmo,lmo),foldlm(kmo,lm),fnew(km,lm),wta(kmo))

!     initialize
      pold(:) = sigo(:)*(psf-10.)+10.
      pnew(:) = pmidl00(1:lm)
      wta = 1.
!     find a top for the output data (a drop sloppy!)
!     Note: ltopx is always higher than lmtc so it doesn't matter..
      do 5 l=1,lm
        ltopx = l
        if (pnew(l).lt.pold(lmo)) go to 10
    5 continue
   10 continue
      divj = 180./dlat_dg !divj is number of whole grid boxes from SP to NP
      call LHNTR0(imo,jmo,-.25d0,22.5d0, im,jm,0.d0,divj,0.d0)

!     Open input and output files
      call openunit('CH4_TROP_FRQ',InFile,.true.,.true.) !uninterpolated
      call openunit(FRQname,interp_file,.true.,.false.)

C**** outer loop over tau
      do 500 it=0,8640,120
      read (InFile) tau,fold
C**** interpolate vertically  fold-->foldlm
      do i=1,kmo
        ain(:) = fold(i,:)
        debug = .false.
        aout = 0.
        call v_int(debug,ltopx,lmo,pold,ain,lm,pnew,aout)
        foldlm(i,:) = aout(:)
      end do
C**** interpolate horizontally  foldlm-->fnew
      do l=1,lm
        call LHNTR(wta,foldlm(1,l),fnew(1,l),imo,jmo,im,jm)
      end do
C**** Write interpolated file in more standard format
      write(title,'(f10.0,a)') tau,' CLIM.RUN.OHCH4.FRQ interpolated'
      write(interp_file) title,fnew
  500 continue
      call closeunit(InFile)
      call closeunit(interp_file)
      DEALLOCATE (fold,foldlm,fnew,wta)
      write(6,*) ' SUBROUTINE get_Trop_chem_CH4_freq executed'
      return
      end subroutine get_Trop_chem_CH4_freq


      subroutine v_int(debug,ltopx,lm_old,pold,ain,lm_new,pnew,aout)
!@sum v_int vertical interpolation for CH4 tables
C**** We want to interpolate from an irregular to an irregular grid.
C**** This is a little sloppy because value at ltopx is not accurate
C****    But what WOULD be correct???
      implicit none
      real*8 ain(*),aout(*),pnew(*),pold(*)
      real*8 dist,pint,plbot,pltop,step
      integer l,lx,lm_new,lm_old,lbot,ltop,ltopx
      logical debug

      step = -1.
      aout(1) = ain(1)
      if (debug) write(6,*)'data in:',(ain(l),l=1,lm_old)

      do 190 l=2,ltopx
      pint = pnew(l)
      plbot = pold(1)
      do lx = 1,lm_old-1     ! avoid reference to pold(lm_old+1)
        lbot = lx
        ltop = lbot+1
        pltop = pold(ltop)
        if (step.gt.0 .and. pint.le.pltop) go to 180
        if (step.lt.0 .and. pint.ge.pltop) go to 180
        if (lx.lt.lm_old-1) plbot = pltop   ! don't end with plbot = pltop
      end do
  180 continue
      dist = (pint-pltop)/(plbot-pltop) !distance from upper boundary
      aout(l) = ain(lbot)*dist+ain(ltop)*(1.-dist)
      if (.not.(aout(l).gt.0..or.aout(l).le.0.)) aout(l) = 0.
      if (debug) write(6,*)l,lbot,ltop,pint,plbot,pltop,
     *  ain(lbot),ain(ltop),aout(l)
  190 continue
      aout(ltopx) = ain(lm_old)   !!! fudgy, but it doesn't matter
      if (debug)write(6,*) 'data out',(aout(l),l=1,lm_new)
      return
      end subroutine v_int


      SUBROUTINE get_wofsy_gas_IC(cgas,gasjl)
C**** Input data are from Wofsy
C**** 1995 CH4 Concentrations in ppb; 1995 CO2 Concentrations in ppm
      USE RESOLUTION, ONLY: im,jm,lm
      USE MODEL_COM, ONLY: amonth,jmon0
      USE DOMAIN_DECOMP_ATM, only: GRID, getDomainBounds, AM_I_ROOT
      USE ATM_COM, only: pednl00
      USE FILEMANAGER, only: openunit,closeunit
      USE GEOM, only : DLAT_DG,lat_dg
      implicit none
      integer j,jw,k,kstart,l,n,iu
      integer, parameter :: kmw=200
      real*4 GASX(3,62),xj,t
      REAL*8 GASW(3,62),
     *  GASJK(GRID%J_STRT_HALO:GRID%J_STOP_HALO,0:kmw),
     *  GASJL(GRID%J_STRT_HALO:GRID%J_STOP_HALO,lm),
     *  P(0:kmw),pup,cup,pdn,cdn,psum,csum,w,zk,scale
      CHARACTER*80 card,titlew,dfile*24
      character*(*) cgas

      INTEGER :: J_1, J_0
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)

C****
C**** Read in GAS concentrations from Wofsy
C**** Data are in jm latitudinal bands, lm layers, km time periods
C**** lat: 1=tropics +/-15, 2=N mid/high, 3=S Mid./high
C**** Data start at 14 km, end at 45 km and are at edges.  Center
C**** of lowest box is at 14.5 km.  There are 28 undefined boxes below,
C**** and 62 defined ones from 14 km up.
C****
      dfile = trim(cgas)//'_IC'
      call openunit(dfile,iu,.false.,.true.)
C**** skip over lines at top
      do n=1,11
        read(iu,'(a)') card
      end do
C**** There are 73 time periods at 5-day intervals
      kstart = nint((jmon0-1.)*73./12.+1.)
C**** process
      do k=1,kstart
        read(iu,*) xj,t,(GASx(2,l),l=1,62) !tropics
        read(iu,*) xj,t,(GASx(3,l),l=1,62) !north
        read(iu,*) xj,t,(GASx(1,l),l=1,62) !south
      end do
      if (AM_I_ROOT())
     *   write(6,'(a,i3,2f9.2)') ' Read Wofsy data at month', jmon0,xj,t
      call closeunit(iu)
C**** Scale to a particular value
      scale = 1.
      if (trim(cgas) == 'CH4') then
        scale = 1.
      else if (trim(cgas) == 'CO2') then
        scale = 334./GASx(1,1)
      end if
      GASW = GASx*scale
C****
C**** Interpolate data GASW to GASJK on GCM latitudes
C****
C**** Below 14 km set all equal

      GASJK(:,0:28) = GASW(1,1)
C**** Keep step function except at two transition points (+/- 15 deg)
      do j=j_0,j_1
        if (lat_dg(j,1).gt.15+0.5*DLAT_DG) then ! nh
          GASJK(J,29:90) = GASW(3,1:62)
        elseif (lat_dg(j,1).lt.-15-0.5*DLAT_DG) then ! sh
          GASJK(J,29:90) = GASW(1,1:62)
        elseif (lat_dg(j,1).gt.-15+0.5*DLAT_DG .and. lat_dg(j,1).lt.15
     *         -0.5*DLAT_DG) then ! tropics
          GASJK(J,29:90) = GASW(2,1:62)
        else                    ! edge points
          W = 1. + (J-1)*DLAT_DG/90.
          JW=W
          GASJK(J,29:90)=GASW(JW,1:62)*(JW+1-W)+GASW(JW+1,1:62)*(W-JW)
        end if
C**** Above, extend
        GASJK(j,91:kmw) = GASJK(j,90)
      end do
C**** Calculate data pressure levels (mb)
C**** z* in km (=7 ln (1000/p)
      zk = 0.
      DO 120 K=0,kmw
      P(K) = 1000.*exp(-ZK/7.)
      ZK = zk+.5
  120 continue
C****
C**** Interpoate GASJ to GASJL on GCM grid boxes conserving vertical
C**** means
C****
      GASJL = 0.
      DO 440 J=J_0,J_1
      PDN = pednl00(1)
      CDN = GASJK(J,0)
      K=1
      DO 430 L=1,lm
      PSUM = 0.
      CSUM = 0.
      PUP  = pednl00(l+1)
  410 IF(P(K).LE.PUP)  GO TO 420
      PSUM = PSUM +  PDN-P(K)
      CSUM = CSUM + (PDN-P(K))*(CDN+GASJK(J,K))/2.
      PDN  = P(K)
      CDN  = GASJK(J,K)
      K=K+1
      if (k.gt.kmw) call stop_model(
     &     ' Please incrase kmw in get_wofsy_gas_IC',255)
      GO TO 410
C****
  420 CUP  = GASJK(J,K) + (GASJK(J,K-1)-GASJK(J,K))*(PUP-P(K))/
     /       (P(K-1)-P(K))
      PSUM = PSUM +  PDN-PUP
      CSUM = CSUM + (PDN-PUP)*(CDN+CUP)/2.
      GASJL(J,L) = CSUM/PSUM
      PDN = PUP
      CDN = CUP
  430 continue
  440 CONTINUE
C****
C**** Write GAS concentration on GCM grid boxes to disk (to check)
C****
c     call openunit('CH4check',iu,.true.)
c     TITLEW =
c    * 'Wolfsy CH4 1995 CONCENTRATION for '//amonth(jmon0)//' 1'
c     WRITE (iu) TITLEW,sngl(GASJL)
c     call closeunit(iu)

      RETURN
      END SUBROUTINE get_wofsy_gas_IC


      SUBROUTINE ALLOC_TRACER_SPECIAL_Lerner_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE PRATHER_CHEM_COM
      USE TRACERS_MPchem_COM
      USE CO2_SOURCES
      USE CH4_SOURCES
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, getDomainBounds
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid
      INTEGER :: J_1H, J_0H, I_1H, I_0H
      INTEGER :: IER
C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H=GRID%I_STRT_HALO
      I_1H=GRID%I_STOP_HALO

      ALLOCATE( jlatmd(J_0H:J_1H),
     *          STAT=IER )
      ALLOCATE(  tltrm(J_0H:J_1H,lm,nMPtable),
     *           tltzm(J_0H:J_1H,lm,nMPtable),
     *          tltzzm(J_0H:J_1H,lm,nMPtable),
     *          STAT=IER )
      ALLOCATE(  CH4_src(I_0H:I_1H,J_0H:J_1H,nch4src),
     *           CO2_src(I_0H:I_1H,J_0H:J_1H,nco2src),
     *           STAT=IER )
      allocate(tscparm(lz_schem,18,INT_MONTHS_PER_YEAR,nMPtable))

      END SUBROUTINE ALLOC_TRACER_SPECIAL_Lerner_COM


      SUBROUTINE ALLOC_LINOZ_CHEM_COM(grid)
!@sum  To allocate arrays whose sizes now need to be determined at
!@+    run time
!@auth NCCS (Goddard) Development Team
      USE LINOZ_CHEM_COM
      USE DOMAIN_DECOMP_ATM, ONLY : DIST_GRID, getDomainBounds
      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H, I_1H, I_0H
      INTEGER :: IER

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = GRID%I_STRT_HALO
      I_1H = GRID%I_STOP_HALO

      ALLOCATE(  TLT0M(J_0H:J_1H,lm,nctable),
     *           TLTZM(J_0H:J_1H,lm,nctable),
     *          TLTZZM(J_0H:J_1H,lm,nctable),
     *          STAT=IER )

      allocate(O3trop_Loss(I_0H:I_1H,J_0H:J_1H,lm,12),
     &         O3trop_Prod(I_0H:I_1H,J_0H:J_1H,lm,12),
     &           O3_DepVel(I_0H:I_1H,J_0H:J_1H,12))

      END SUBROUTINE ALLOC_LINOZ_CHEM_COM


