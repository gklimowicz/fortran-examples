C****  
C**** OCNGISSVM.f  GISS ocean vertical mixing scheme, 12/28/2012
C****
#include "rundeck_opts.h"

      MODULE GISSMIX_COM
!@sum  GISSMIX_COM holds variables related to the GISS mixing scheme
!@ref Canuto et al. 2010, Ocean Modelling, 34, 70-91 (C2010)
!@auth AHoward/YCheng
#ifdef TRACERS_OCEAN
c     USE OCN_TRACER_COM, only : tracerlist
#endif
      USE OCEAN, only : im,jm,lmo
      USE CONSTANT, only : omega,by3,grav

      IMPLICIT NONE
      SAVE

      integer, parameter :: mt0=107,mt=2*mt0  !@var mt dim of ri in table
      integer, parameter :: nt0=54, nt=2*nt0  !@var nt dim of rr in table
      real*8, parameter :: kmin=1d-3    !@var kmin min of diffusivities, (m^2/s)
      real*8, parameter :: kmax=100.    !@var kmax max of diffusivities, (m^2/s)
      real*8, parameter :: osocb1=21.6,kappa=0.4
      real*8 ::
     &    ria(mt)       !@var ria ri 1d array of richardson #, for 2d tables
     &   ,rra(nt)       !@var rra rr 1d array of density ratio,for 2d tables
     &   ,gma(mt,nt)    !@var gma 2d table for gm=(tau*shear)^2
     &   ,sma(mt,nt)    !@var sma 2d table for sm=structure fuction for momentum
     &   ,sha(mt,nt)    !@var sha 2d table for sh=structure fuction for heat
     &   ,ssa(mt,nt)    !@var ssa 2d table for ss=structure fuction for salinity
     &   ,sca(mt,nt)    !@var sca 2d table for sc=structure fuction for passive scalars
     &   ,phim2a(mt,nt) !@var phim2a 2d table for phim2, used for bottom shear
      real*8 :: 
     &    rimin         !@var rimin min of richardson # ri
     &   ,rimax         !@var rimax max of richardson # ri
     &   ,rrmin         !@var rrmin min of density ratio rr
     &   ,rrmax         !@var rrmax max of density ratio rr

!@var otke turbulent kinetic energy in ocean (m/s)^2
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: otke
!@var otke_init_max maximum initial value of otke (m/s)^2
      real*8, parameter :: otke_init_max=0.5d0/800.
!@var emin minimum value of otke (m/s)^2
!@var emax maximum value of otke (m/s)^2
      real*8, parameter :: emin=1d-6,emax=1000.
#ifdef TRACERS_OCEAN
c     REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: TRMO1,TXMO1,TYMO1
#endif
      ! bottom drag over lon,lat used by bottom drag routine
      ! C2010 section 7.2 equation (72)
!@var rhobot(:,:) in-situ density at ocean bottom (kg/m^3)
!@+   which is calculated in subroutine OCONV and passed from here
!@+   and will be used to calculate the ocean bottom drag in the GCM host
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: rhobot
!@var taubx x component of tau_b, kinematic bottom drag in (m/s)^2
!@var tauby y component of tau_b, kinematic bottom drag in (m/s)^2
!@+   taubx and tauby are calculated in subroutine gissmix
!@+   and will be used to calculate the ocean bottom drag in the GCM host
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: taubx
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: tauby
!@var exya internal tidal energy (w/m^2)
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: exya
!@var ut2a unresolved bottom velocity squared (m/s)^2
      REAL*8, ALLOCATABLE, DIMENSION(:,:) :: ut2a
!@var idrag =1: tides produce bottom drag;
!@+         =0: tides do not produce bottom drag
!@+   idrag is used in subroutine OBDRAG in OCNDYN.f and OCNDYN2.f
      integer, parameter :: idrag=1

      END MODULE GISSMIX_COM

      MODULE GISS_OTURB
!@sum GISS_OTURB contains variables and routines for GISS mixing scheme
!@ref Canuto et al. 2010, Ocean Modelling, 34, 70-91 (C2010)
!@+   Canuto et al. 2011, Ocean Modelling, 36, 198-207 (C2011)
!@+   Cheng et al. 2002, JAS, 59,1550-1565 (C2002)
!@auth AHoward/YCheng

      USE GISSMIX_COM
      implicit none

      ! Variables used for background diffusivities
      real*8 :: bv0, f30, bv0byf30, byden, epsbyn2, q, byzet 

      CONTAINS

      subroutine gissmix_init(iniOCEAN)
!@sum creates tables for the giss diffusivities(ri,rr) 
!@+   gissmix_init is called from routine init_OCEAN in OCNDYN.f
!@auth AHoward/YCheng
!@ref Canuto et al. 2010, Ocean Modelling, 34, 70-91 (C2010)
!@+   Canuto et al. 2011, Ocean Modelling, 36, 198-207 (C2011)
!@+   Cheng et al. 2002, JAS, 59,1550-1565 (C2002)
      USE FILEMANAGER
      USE DOMAIN_DECOMP_1D, only : READT_PARALLEL
      USE OCEANR_DIM, only : ogrid

      implicit none

      ! in:
      logical, intent(in) :: iniOCEAN

      ! local:

      integer j,k,l,iu_TIDES
      real*8 tmpm(mt),tmpn(nt),rini,rend,ff

      ! establish ri grids for look-up table
      ! the grids extend from -rend to rend, with more grids near zero

      rini=1d-4
      rend=1d4
      ff=(mt0-1)/(dlog(rend/rini)/dlog(2.d0))
      do j=1,mt0-1
         tmpm(j)=rini*2**(dble(j-1)/ff)
      end do
      tmpm(mt0)=rend
      do j=1,mt
         if(j.le.mt0) then
            ria(j)=-tmpm(mt0-j+1)
         else
            ria(j)=tmpm(j-mt0)
         endif
      end do

      ! establish rr grids for look-up table
      ! the grids extend from -1 to rend, with more grids near 0
      ! abs(rr)>=rini
      rini=1d-2
      rend=1.-1d-2
      ff=(nt0-1)/(dlog(rend/rini)/dlog(2.d0))
      do k=1,nt0-1
         tmpn(k)=rini*2**(dble(k-1)/ff)
      end do
      tmpn(nt0)=rend
      do k=1,nt
          if(k.le.nt0) then
             rra(k)=-tmpn(nt0-k+1)
          else
             rra(k)=tmpn(k-nt0)
          endif
       end do
       rra(1)=-1.

      ! 2d tables for gm,sm,sh,ss,sc,rf,phim2
      ! phim2 is defined in C2002, (36), with l=kz; for bottom shear

      rimin=ria(1) ! -1.d4
      rimax=ria(mt) !  1.d4
      rrmin=rra(1) ! -1.
      rrmax=rra(nt) !  0.99
      do k=1,nt
         do j=1,mt
            call gissdiffus(
            ! In:
     &         ria(j),rra(k)
            ! Out:
     &        ,gma(j,k),sma(j,k),sha(j,k),ssa(j,k),sca(j,k))
            phim2a(j,k)=2./osocb1**2*gma(j,k)**.5d0/sma(j,k)
         end do ! loop j
      end do ! loop k

C**** initialize exya, ut2a
c tidally induced diffusivities: C2010, (69), (75)-(77)
c the data file opened below, TIDAL_e_v2_1QX1, is of format
c character*80::title, real*4::data(288,180)
c title1: Wave Bottom Tidal Dissipation in Watts/(meter^2)
c title2: (Tidal velocity)^2 from Wave Data in (meters/second)^2
  
      call openunit("TIDES",iu_TIDES,.true.,.true.)
      CALL READT_PARALLEL(ogrid,iu_TIDES,NAMEUNIT(iu_TIDES),exya ,1)
      CALL READT_PARALLEL(ogrid,iu_TIDES,NAMEUNIT(iu_TIDES),ut2a ,1)

C**** initialize rhobot, taubx, tauby
      rhobot(:,:)=0.
      taubx(:,:)=0.
      tauby(:,:)=0.

C**** initialize otke
      if(iniOCEAN) then
         do l=1,lmo
            otke(l,:,:)=min(max(otke_init_max/(float(l)**2),emin),emax)
         end do
      endif

      bv0 = 5.24d-3                     ! (1/s)
      f30 = omega                       ! (1/s)
      bv0byf30 = bv0/f30                ! (1)
      byden = 1./(f30*acosh(bv0byf30))  ! (1)
      epsbyn2 = .288d-4                 ! (m^2/s)
      q = .7d0        ! fraction of baroclinic energy into creating mixing
      byzet = 1./500.d0                 ! upward decaying factor (1/m)

      return
      end subroutine gissmix_init

      subroutine gissdiffus(
      ! In:
     &    ri,rr
      ! Out:
     &   ,gm,sm,sh,ss,sc)
!@sum finds structure functions sm,sh,ss,sc using giss model
!@ref Canuto et al. 2010, Ocean Modelling, 34, 70-91 (C2010)
!@auth AHoward/YCheng
      
      implicit none

      ! in:
      real*8 :: ri,rr
      ! out:
      real*8 :: gm,sm,sh,ss,sc
      ! local:
      real*8 :: pi1,pi2,pi3,pi4,pi5,a3,a2,a1,a0,p1,p2
      real*8 :: tmp,a,b,c
      real*8 :: x1r,x2r,x3r,x1i,x2i,x3i,omrr
      real*8 :: gh,gr,q,p,bygam,xx,am,ah,as,ac,ar,w2byk
      complex*16 x1,x2,x3

      ! prepare for solving follwoing cubic eqn for gm
      ! a3*gm**3+a2*gm**2+a1*gm+a0=0

      call cubic_coeff(
      ! In:
     &    ri,rr
      ! Out:
     &   ,pi1,pi2,pi3,pi4,pi5
     &   ,a3,a2,a1,a0)
      if(abs(a3).lt.1d-10.and.abs(a2).lt.1d-10) then
         gm=-a0/a1
      elseif(abs(a3).lt.1d-10) then
         ! solve quadratic a2*x**2+a1*x+a0=0
         ! x=gm
         tmp=sqrt(a1**2-4*a2*a0)
         gm=(-a1-tmp)/(2*a2)
         if(gm.lt.0.) then
            gm=(-a1+tmp)/(2*a2)
         endif
      else
         ! solve cubic x**3+a*x**2+b*x+c=0
         ! x=gm
         a=a2/a3
         b=a1/a3
         c=a0/a3
         call cubic_r(a,b,c,x1,x2,x3)
         x1r=dble(x1); x1i=aimag(x1);
         x2r=dble(x2); x2i=aimag(x2);
         x3r=dble(x3); x3i=aimag(x3);
         if(x1r.lt.0.) x1r=1d40
         if(x2r.lt.0.) x2r=2d40
         if(x3r.lt.0.) x3r=3d40
         if(abs(x1i).gt.1d-30) x1r=1d30
         if(abs(x2i).gt.1d-30) x2r=2d30
         if(abs(x3i).gt.1d-30) x3r=3d30
         gm=min(x1r,min(x2r,x3r))
         if(gm.ge.1d30) then
            write(*,*) "gm=",gm
            write(*,*) "why gm >= 1d30; stop and check"
            stop
         endif
      endif
      if(gm.lt.0.) then
         write(*,*) "gm < 0", gm
         stop
      endif
      omrr=1-rr
      gh=ri*gm/omrr
      q=pi1*(pi2*(1+rr)-pi3*rr)
      p=pi4*(pi5-pi2*(1+rr))
      bygam=rr/( (pi4/pi1)*(1+q*gh)/(1+p*gh) )
      ah=pi4/(1+p*gh+pi2*pi4*gh*(1-bygam))
      xx=(1-bygam)*gh*ah
      ar=xx/(ri*gm)
      p1=(1+rr**2)/(1/pi1+rr**2/pi4)
      p2=pi2
      gr=ri*gm
      am=2/gm*(15./7.+xx)
      w2byk=2. / (30./7.+xx)
      as=ah/( (pi4/pi1)*(1+q*gh)/(1+p*gh) )
      ac=p1/(1+p1*p2*gr)*(1-p2*gr*ar)
      sm=am*w2byk
      sh=ah*w2byk
      ss=as*w2byk
      sc=ac*w2byk
      if((gm.le.0.).or.(sm.le.0.).or.(sh.le.0.).or.(ss.le.0.).
     &  or.(sc.le.0.)) then
         write(100,*) "gm,sm,sh,ss or sc <= 0., stop and check"
         write(100,*)  gm,sm,sh,ss,sc
         stop
      endif
      !sr=(sh-rr*ss)/(1-rr)
      !rf=ri*sr/sm

      return
      end subroutine gissdiffus

      subroutine cubic_coeff(
         ! In:
     &      ri,rr
         ! Out:
     &     ,pi1,pi2,pi3,pi4,pi5
     &     ,a3,a2,a1,a0
     &     )
!@sum calculates the coefficients of the cubic eqn for gm
!@+   a3*gm**3+a2*gm**2+a1*gm+a0=0
!@ref Canuto et al. 2010, Ocean Modelling, 34, 70-91 (C2010)
!@ref Canuto et al. 2013, in preparation
!@auth AHoward/YCheng

      implicit none

      real*8 ri,rr
      real*8 pi1,pi2,pi3,pi4,pi5
      real*8 a3,a2,a1,a0

      real*8 sgmt0,b,rrsy,tmp,aa1,aa2,aa3,aa4,aa5,aa6
      real*8 pi10,pi20,pi30,pi40,pi50

      !@ pi's are the time scale ratios

      sgmt0=.72d0
      pi20=1./3.d0
      pi40=1./5.*(1./(1.+1/sgmt0))
      pi10=pi40
      pi50=sgmt0
      pi30=pi50
      pi3=pi30
      pi5=pi50
      if(ri.lt.0.) then
         pi1=pi10
         pi2=pi20
         pi4=pi40
      else ! ri >= 0
         ! the pi's have been updated by Canuto et al. 2013 
         if(rr.gt.0.) then
            b=2.*ri/(.1+ri)
            rrsy = 2./(rr + 1/rr)
            tmp=atan2(.5*rrsy**2*(1-rrsy**2),2*rrsy**2-1)/3.1415927
            pi1 = pi10/(1.+ri**tmp*tmp)
            pi4 = pi1
            pi2=pi20*( 1-b*rrsy*(1-rrsy) )
         else ! rr < 0
            pi1=pi10/(1+ri)
            pi2=pi20
            pi4=pi40/(1+ri)
         endif
      endif

!@    a3*gm**3+a2*gm**2+a1*gm+a0=0

      aa1 = pi1 * pi4 *
     &  (    pi2 * (15 * pi3 + 7) *(rr**2 +1) + (14 * pi2
     &  - 14 * pi3 - 15 * pi3 ** 2) * rr    )
     & * (pi1 *  rr - pi4) / (-1 + rr) ** 3 / 150

      aa2 = pi1 * pi4 * 
     &  (    pi2 * (-150 * pi3 + 210 * pi1 + 7) *( rr **2+1)
     &    +( 14*(pi2-pi3)*(1+15*(pi1+pi4))+150*pi3**2 )*rr
     &    + 210* pi2 * ( pi4-pi1)    )
     &  / (-1 + rr) ** 2  / 9000

      aa3 = (pi1 * ( 5*pi2*pi4*(30 * pi3  + 17 )
     &      + pi1*(15 * pi3 + 7) ) * (rr ** 2+1)
     &   - ( 10*pi1*pi3*pi4*(15*pi3+17) + 15*pi2*(pi4**2+pi1**2)
     &    + 14*pi1*pi4*(1-10*pi2) ) * rr
     &       - (15*pi3+7)*(pi1**2-pi4**2))
     &      / (-1 + rr) ** 2 / 150

      aa4 = (( 150 *(pi1*pi3+pi4*pi2)-7 *pi1*(1+30*pi1) )*rr
     &    - 150*(pi1*pi2+pi4*pi3)+7 *pi4*(1 + 30 * pi4)  )
     &      / (1 - rr) / 9000

      aa5 =  ((-17 * pi1 - 30 *( pi1 * pi3 +  pi4 * pi2)) * rr  
     &          +17 * pi4 + 30 *( pi1 * pi2 +  pi4 * pi3))
     &    / (1 - rr) / 30

      aa6 = -1. / 60.

      a3 = aa1*ri**3+aa2*ri**2
      a2 = aa3*ri**2+aa4*ri
      a1 = aa5*ri+aa6
      a0 = 1.

      return
      end subroutine cubic_coeff
      
      subroutine cubic_r(a,b,c,x1,x2,x3)
!@sum solves the cubic eqn x**3+a*x**2+b*x+c=0
!@+   for a, b, c are real
      implicit none
      real*8 a,b,c
      complex*16 x1,x2,x3
      real*8 pi,q,r,the,temp,aa,bb
      pi=acos(-1.d0)
      q=(a**2-3.d0*b)/9.d0
      r=(2.d0*a**3-9.d0*a*b+27.d0*c)/54.d0
      if (r**2.lt.q**3) then
          the=acos(r/sqrt(q**3))
          temp=2.d0*sqrt(q)
          x1=-temp*cos(the/3.d0)-a/3.d0
          x2=-temp*cos((the+2.d0*pi)/3.d0)-a/3.d0
          x3=-temp*cos((the-2.d0*pi)/3.d0)-a/3.d0
      else
          aa=-sign(1.d0,r)*(abs(r)+sqrt(r**2-q**3))**(1.d0/3.d0)
          if(aa.ne.0.d0) then
              bb=q/aa
          else
              bb=0.d0
          endif
          x1=(aa+bb)-a/3.d0
          x2=cmplx(-0.5d0*(aa+bb)-a/3.d0, sqrt(3.d0)/2.d0*(aa-bb))
          x3=conjg(x2)
      endif
      return
      end subroutine cubic_r

      subroutine locatex(n,xa,x,jl,ju,a,b)
!@sum locatex finds the grids that embrace x
!@+   after call locatex, the interpreted value can be calculated as
!@+   y=a*ya(jl)+b*ya(ju)
      implicit none
      ! in:
      integer n
      real*8 xa(n),x
      ! out:
      integer jl,ju
      real*8 a,b ! a=(xa(ju)-x)/h, b=(x-xa(jl))/h
      ! local:
      integer j
      real*8 h
      jl=1
      ju=n
      do while (ju-jl.gt.1)
        j=(ju+jl)/2
        if(xa(j).gt.x)then
          ju=j
        else
          jl=j
        endif
      end do
      h=xa(ju)-xa(jl)
      if (h.eq.0.) then
         write(*,*) 'h=0 in locatex'
         stop
      endif
      a=(xa(ju)-x)/h
      b=(x-xa(jl))/h
      return
      end subroutine locatex


      END MODULE GISS_OTURB


      subroutine gissmix( 
      ! in:
     &    n,ze,zg,db,dv2,adt,bds,rho,ustarb2,exy,fc,hbl,strait
      ! out:
     &   ,ri,rr,bv2,km,kh,ks,kc,e) 

!@sum giss turbulence model for ocean
!@ref Canuto et al. 2004, GRL, 31, L16305 (C2004)
!@+   Canuto et al. 2010, Ocean Modelling, 34, 70-91 (C2010)
!@+   Canuto et al. 2011, Ocean Modelling, 36, 198-207 (C2011)
!@auth AHoward and YCheng
!@var grav acceleration of gravity m/s^2
!@var lmo max. number of vertical layers
!@var taubx(IM,J_0H:J_1H) x component of tau_b = kinematic bottom drag (m/s)^2
!@var tauby(IM,J_0H:J_1H) y component of tau_b = kinematic bottom drag (m/s)^2
!@var exya internal tidal energy (w/m^2)
!@var ut2a unresolved bottom shear squared (m/s)^2

      USE GISS_OTURB
!@ MODULE OTURB's variables will be used and/or updated
!@ for example, taubx,tauby will be updated
      implicit none

      ! in:
      integer n          !@var n number of vert. layers on this column
      real*8 ze(0:lmo)   !@var ze vertical grid-edge depth (m), > 0
      real*8 zg(0:lmo+1) !@var zg vertical grid depth (m), < 0
      real*8 db(lmo)     !@var db -grav/rho*d(rho) (m/s^2)
      real*8 dv2(lmo)    !@var dv2 vel. diff. squared btw layers (m/s)^2
      real*8 adt(lmo)    !@var adt rho*alpha*DT (kg/m^3)
      real*8 bds(lmo)    !@var bds rho*beta*DS  (kg/m^3)
      real*8 rho(lmo)    !@var rho density
      real*8 ustarb2     !@var ustarb2 velocity squared at zg(n)
      real*8 exy         !@var exy tidal power input
      real*8 fc          !@var fc Coriolis parameter=2*omega*sin(lat) (1/s)
      real*8 hbl         !@var hbl pbl depth (m)
      integer strait     !@var strait 0: not in strait; 1: in strait
      intent (in) n,ze,zg,db,dv2,adt,bds,rho,ustarb2,exy,fc,hbl,strait

      ! out:
      real*8 ri(0:lmo+1) !@var ri local richardson number
      real*8 rr(0:lmo+1) !@var rr local alpha*dTdz/(beta*dSdz)
      real*8 bv2(0:lmo+1)!@var bv2 Brunt Vaisala frequency squared (1/s**2)
      real*8 km(0:lmo+1) !@var km vertical momentun diffusivity (m**2/s)
      real*8 kh(0:lmo+1) !@var kh vertical heat diffusivity (m**2/s)
      real*8 ks(0:lmo+1) !@var ks vertical salinity diffusivity (m**2/s)
      real*8 kc(0:lmo+1) !@var kc vertical passive scalar diffusivity (m**2/s)
      real*8 e(lmo)      !@var e ocean turbulent kinetic energy (m/s)**2
      intent (out) ri,rr,bv2,km,kh,ks,kc,e

      ! local:
      integer :: flag !@var flag =0 if abs(rr)<=1; =1 if abs(rr)>1
      ! num_smooth=1: to smooth bv2,vs2 and rr; num_smooth=0: no smooth
      integer, parameter :: num_smooth=0
      !@var l0 constant length scale within obl 
      !@var l1 length scale throughout the vertical ocean
      !@var l1min minimum of l1 below obl 
      !@var lr length scale reduction factor by stable buoyancy
      !@var l2 length scale after reduced by stable buoyancy 
      !@var l2min minimum of l2
      !@var len final length scale
      real*8, parameter :: l1min=3.d0,l2min=.05d0 ! (m)
      integer iter, mr
      real*8 vs2(0:lmo+1)!@var vs2 velocity shear squared (1/s**2)
      real*8 len         !@var len turbulence length scale (m)

      integer l,jlo,jhi,klo,khi
      real*8 a1,a2,b1,b2,c1,c2,c3,c4
      real*8 ril,rrl,gm,sm,sh,ss,sc,kml,khl,ksl,kcl,lr,etau
      real*8 l0,l1,l2,kz,zbyh,bydz,zl,tmp,tmp1
      real*8 fbyden,afc,ltn,bvbyf,fac,kmbg,khbg,ksbg
      real*8 den,fz,epstd_byn2,kmtd,khtd,kstd
      real*8 phim2,zb,unr20

      ! Vertical grid diagram
      !
      !             grid levels                       interface levels
      !
      !                   -----------------------------  surf
      !                1    - - - - - - - - - - - - -
      !                   -----------------------------  1
      !                2    - - - - - - - - - - - - -
      !                   -----------------------------  2
      !               l-1   - - - - - - - - - - - - -
      !                   -----------------------------  l-1
      !zg(l)<0,tracer  l    - - - - - - - - - - - - -
      !                   -----------------------------  l  ze(l)>=0,edge
      !               l+1   - - - - - - - - - - - - -
      !                   -----------------------------  l+1
      !               lm    - - - - - - - - - - - - -
      !                   -----------------------------  lm  
      !             lm+1    - - - - - - - - - - - - -
      !                   -----------------------------  lm+1

      !-----------------------------------------------------------------
      ! more details of giss model are included in subroutine gissdiffus
      ! which is called from within the subroutine gissmix_init, the
      ! latter is called only once from routine init_OCEAN in OCNDYN.f
      !-----------------------------------------------------------------

      do l=1,n-1
         bydz=1./(zg(l)-zg(l+1))
         bv2(l)=db(l)*bydz ! N^2
         vs2(l)=dv2(l)*bydz*bydz
      end do
      do mr = 1,num_smooth ! no smoothing if num_smooth=0
         call z121(bv2,n-1,lmo)
         call z121(vs2,n-1,lmo)
      end do
      do l=1,n-1
         ri(l)=bv2(l)/(vs2(l)+1d-30)
         rr(l)=bds(l)/(adt(l)+1d-30)
      end do
      do mr = 1,num_smooth
        call z121(rr,n-1,lmo)
      end do

      l0=.15*hbl

      ! modify ri at the interface nearest to ocean bottom
      ! due to unresolved bottom shear, C2010, eqs,(73)-(77)
      ! using C2002, (35)-(36), find phim(ri), make a table for phim2(ri)
      ! (generalized to include rr dependence)
      ! phim=f/(1-a*f*rf), f=sqrt(2.)/b1*(gm/sm**2)**.25, a=0 or 2.7
      ! iterate for ri=n2/(sig2+unr2)

      if(strait.eq.0) then
         l=n-1
         zb=ze(n)-ze(l)
         unr20=ustarb2/(kappa*zb)**2
         rrl=rr(l)
         if(abs(rrl).gt.1.) rrl=1/(rr(l)+1d-30)
         rrl=min(max(rrl,rrmin),rrmax)
        
         ! locatex uses bisection method to lookup tables
         call locatex(nt,rra,rrl,klo,khi,a2,b2)
         ril=ri(l)
         do iter=1,10
            ril=min(max(ril,rimin),rimax)
            call locatex(mt,ria,ril,jlo,jhi,a1,b1)
            phim2=a2*(a1*phim2a(jlo,klo)+b1*phim2a(jhi,klo))
     &           +b2*(b1*phim2a(jhi,khi)+a1*phim2a(jlo,khi))
            tmp=ril
            ril=bv2(l)/(vs2(l)+unr20*phim2+1d-30)
            ril=.9d0*ril+.1d0*tmp ! best
            if(abs((ril-tmp)/(ril+tmp)).le.1d-2) exit
         end do
         ri(l)=ril
      endif
      
      den=1.-exp(-ze(n)*byzet)
      afc=abs(fc)
      fbyden=afc*byden

      do l=1,n-1

         ril=ri(l)
         rrl=rr(l)
         flag=0
         if(abs(rrl).gt.1.) then
            rrl=1/(rr(l)+1d-30)
            flag=1
         endif

         ! find foreground diffusivities from 2d lookup tables

         ril=min(max(ril,rimin),rimax)
         rrl=min(max(rrl,rrmin),rrmax)
         call locatex(mt,ria,ril,jlo,jhi,a1,b1)
         call locatex(nt,rra,rrl,klo,khi,a2,b2)
         c1=a1*a2
         c2=b1*a2
         c3=b1*b2
         c4=a1*b2
         gm=c1*gma(jlo,klo)+c2*gma(jhi,klo)
     &     +c3*gma(jhi,khi)+c4*gma(jlo,khi)
         sm=c1*sma(jlo,klo)+c2*sma(jhi,klo)
     &     +c3*sma(jhi,khi)+c4*sma(jlo,khi)
         sh=c1*sha(jlo,klo)+c2*sha(jhi,klo)
     &     +c3*sha(jhi,khi)+c4*sha(jlo,khi)
         ss=c1*ssa(jlo,klo)+c2*ssa(jhi,klo)
     &     +c3*ssa(jhi,khi)+c4*ssa(jlo,khi)
         sc=c1*sca(jlo,klo)+c2*sca(jhi,klo)
     &     +c3*sca(jhi,khi)+c4*sca(jlo,khi)
         ! symmetry of giss model: sh <-> ss if rr<->1/rr
         if(flag.eq.1) then
            tmp=sh
            sh=ss
            ss=tmp
         endif

         ! length scale:
         zl=ze(l)
         zbyh=zl/hbl
         kz=kappa*zl
         !@var lr length scale reduction factor by stable buoyancy
         !@var l0 constant length scale within obl 
         !@var l1 length scale before reduced by stable buoyancy 
         !@var l1min minimum of l1 below obl 
         !@var l2 length scale after reduced by stable buoyancy 
         !@var l2min minimum of l2
         !@var len final length scale
         lr=1./(1.+max(ril,0.d0))
         if(zl.le.hbl) then   ! within obl
            l1=l0
         else
            l1=l1min+max(l0-l1min,0.d0)*exp(1.-zbyh)
         endif
         l2=max(l1*lr,l2min)
         len=l2*kz/(l2+kz)
         tmp=(osocb1*len)**2*vs2(l)
         e(l)=.5*tmp/(gm+1.d-20)
         e(l)=min(max(e(l),emin),emax)
         etau=.5*osocb1*sqrt(2.*e(l))*len
         kml=etau*sm
         khl=etau*sh
         ksl=etau*ss
         kcl=etau*sc

         ! background and tidally induced diffusivities

         if(ril.gt.0.) then
            ! C2010, eqs.(65a)-(66); C2011, near end of Sec 4, p.203
            ! afc = abs(Coriol)
            bvbyf=sqrt(max(bv2(l),1d-8))/(afc+1.d-30) 
            if(bvbyf.gt.1.) then
               ltn=acosh(bvbyf)*fbyden  ! dimensionless
               ltn=max(ltn,7.d-2)       ! limit described in C2004
            else
               ltn=7.d-2                ! dimensionless
            endif
            fac=epsbyn2*ltn      ! in m^2/s, Km=GAMMAm*fac
            kmbg=fac
            khbg=by3*fac
            ksbg=khbg

            ! tidally induced diffusivities, C2010, eqs.(69)-(71)
            ! use the following fz instead of (70) of C10:
            fz=(exp((-zg(l+1)-ze(n))*byzet)
     &         -exp((-zg(l)  -ze(n))*byzet))/(den*(zg(l)-zg(l+1)))
            epstd_byn2=q*exy*fz*2./(rho(l)+rho(l+1))/max(bv2(l),1d-8)
            kmtd=epstd_byn2
            khtd=by3*epstd_byn2
            kstd=khtd
         else
            kmbg=0.
            khbg=0.
            ksbg=0.
            kmtd=0.
            khtd=0.
            kstd=0.
         endif      

         ! foreground, background and tidal diffusivities added up
         ! C2010, eqs.(69)-(71)

         km(l)=min(kml+kmbg+kmtd,kmax)
         kh(l)=min(khl+khbg+khtd,kmax)
         ks(l)=min(ksl+ksbg+kstd,kmax)
         kc(l)=min(kcl+ksbg+kstd,kmax)
         
      end do  
      km(0)=0.; kh(0)=0.; ks(0)=0.; kc(0)=0.
      km(1)=max(km(1),kmin);kh(1)=max(kh(1),kmin);ks(1)=max(ks(1),kmin)
      kc(1)=max(kc(1),kmin)
      km(n:lmo+1)=0.; kh(n:lmo+1)=0.; ks(n:lmo+1)=0.; kc(n:lmo+1)=0.
      ri(0)=0.; rr(0)=0.; bv2(0)=0.
      ri(n:lmo+1)=0.; rr(n:lmo+1)=0.; bv2(n:lmo+1)=0.; e(n:lmo)=emin
      
      return
      end subroutine gissmix


#ifdef OCN_GISSMIX /*alloc_gissmix_com,def_rsf_gissmix,new_io_gissmix*/
      ! the three routines below are for allocation and input/output

      SUBROUTINE alloc_gissmix_com(grid)
!@sum  To allocate arrays who sizes now need to be determined at
!@+    run-time; called from UBROUTINE alloc_ocean in OCEAN_COM.f
!@auth Reto Ruedy/Ye Cheng

      USE DOMAIN_DECOMP_1D, only : dist_grid,getDomainBounds
!      USE OCEANR_DIM

      USE GISSMIX_COM

      IMPLICIT NONE
      TYPE (DIST_GRID), INTENT(IN) :: grid

      INTEGER :: J_1H, J_0H, L
      INTEGER :: IER

      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      ALLOCATE( otke(lmo,IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( rhobot(IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( taubx(IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( tauby(IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( exya(IM,J_0H:J_1H) , STAT = IER)
      ALLOCATE( ut2a(IM,J_0H:J_1H) , STAT = IER)

#ifdef TRACERS_OCEAN
c     ALLOCATE( TRMO1(tracerlist%getsize(),IM,J_0H:J_1H),
c    *          TXMO1(tracerlist%getsize(),IM,J_0H:J_1H),
c    *          TYMO1(tracerlist%getsize(),IM,J_0H:J_1H),
c    *   STAT = IER)
#endif

      END SUBROUTINE alloc_gissmix_com

      subroutine def_rsf_gissmix(fid)
!@sum  this subroutine is called at the end of
!@+    subroutine def_rsf_ocean, the latter is in OCNDYN.f
!@+    remember to put #ifdef OCN_GISSMIX around the call 
!@date Nov 29,2011
      use pario, only : defvar
      use oceanr_dim, only : grid=>ogrid
      use gissmix_com, only : otke
      use straits, only : otkest
      implicit none
      integer fid   !@var fid file id
      call defvar(grid,fid,otke,'otke(lmo,dist_imo,dist_jmo)')
      call defvar(grid,fid,otkest,'otkest(lmo,nmst)')
      return
      end subroutine def_rsf_gissmix

      subroutine new_io_gissmix(fid,iaction)
!@sum  this subroutine is called at the end of
!@+    subroutine new_io_ocean, the latter is in OCNDYN.f
!@+    remember to put #ifdef OCN_GISSMIX around the call 
!@date Nov 29,2011
      use model_com, only : ioread,iowrite
      use pario, only : read_dist_data,write_dist_data
     &                 ,read_data,write_data
      use oceanr_dim, only : grid=>ogrid
      use gissmix_com, only : otke
      use straits, only : otkest
      implicit none
      integer fid     !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      select case (iaction)
      case (iowrite)            ! output to restart file
        call write_dist_data(grid,fid,'otke',otke,jdim=3)
        call write_data(grid,fid,'otkest',otkest)
      case (ioread)            ! input from restart file
        call read_dist_data(grid,fid,'otke',otke,jdim=3)
        call read_data(grid,fid,'otkest',otkest,bcast_all=.true.)
      end select
      return
      end subroutine new_io_gissmix

#endif /*alloc_gissmix_com,def_rsf_gissmix,new_io_gissmix*/
