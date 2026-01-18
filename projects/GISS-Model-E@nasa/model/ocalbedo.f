#include "rundeck_opts.h"      

      module ocalbedo_mod

      implicit none

      integer, parameter :: nlt=33                ! number of spectral channels
      real*8, dimension(nlt), private :: wfac
      real*8, dimension(nlt), protected :: aw, bw ! absorption,scattering coefficients of water

      integer, dimension(nlt), protected :: lam   ! wavelength in nm
      logical, private :: initialized=.false.

      contains

      subroutine ocalbedo(wind,solz,bocvn,xocvn,chl,
     .                         rod,ros,hycgr,i,j)

***********************************************************************
***** this routine is used by both atmosphere and ocean at each (i,j)
***** where implicitly it is assumed that
***** wind,solz,bocvn,xocvn,chl are in the atmos gird if hycgr=.false. 
***** wind,solz,                are in the ocean gird if hycgr=.true.
***** when this routine is called from within ocean, it does not compute
***** albedo coefficients. Those are computed on the atmos grid.
***** when this routine is called from within ocean, it does not pass
***** reflectances
***********************************************************************

c  Computes ocean surface albedo from solar zenith angle (solz)
c  and wind speed (wind, m/s).
c  Albedo is provided as direct (albd) and diffuse (albs).
c  Derive surface reflectance as a function of solz and wind
c  Includes spectral dependence of foam reflectance derived from Frouin
c  et al., 1996 (JGR)
      USE CONSTANT, only : radian

      implicit none

      real*8, intent(in)  :: wind, solz, chl
      logical, intent(in) :: hycgr
      integer, intent(in) :: i,j
      real*8, dimension(6), intent(out) :: bocvn, xocvn
      real*8, dimension(:), intent(out) :: rod, ros

      integer nl
      real*8 cn,rof,rosps,rospd,rtheta
      real*8 sintr,rthetar,rmin,rpls,sinrmin,sinrpls,tanrmin
      real*8 tanrpls,sinp,tanp,a,b

      real*8 :: sunz
 
      real*8 :: sum1, sum2, part_sum
      logical :: res
      real*8, dimension(nlt) :: refl

!!!!!!!!!!  Boris' part !!!!!!!!!!!!!!!!!
!@sum Those are the weights which were obtained by normalizing solar flux
!@sum for the Lamda's range 0 - 4000 nm. We used Landau fitting function for it.
!@sum They are used to be used for getting 6 band approximation based on 33
!@sum Watson Gregg band calculations
!@sum band_6(j) = Sum(band_33(i)*weight(i))/Sum(part_sum(i))
C**** WHY IS WEIGHT ONLY DECLARED TO BE 31 AND NOT 33?
      real*8 :: weight(31) = (/0.0158378,0.0201205,0.0241885,0.0277778,
     .                       0.0307124,0.0329082,0.0343586,0.0351143,
     .                       0.0352609,0.0349008,0.0341389,0.0330742,
     .                       0.0317941,0.0303725,0.0288696,0.0273329,
     .                       0.0500921,0.0643897,0.0686573,0.0532013,
     .                       0.0416379,0.0330341,0.0265929,0.0217156,
     .                       0.0179725,0.0150596,0.0127618,0.0158128,
     .                       0.0232875,0.0313132,0.0184843/)

C**** gband is the distribution of the 31 bands for the 6 band GISS code
      integer :: gband(7) = (/ 1, 18, 19, 23, 26, 30, 32 /)
c      real*8 :: part_sum(6) = (/0.526854,0.0643897,0.196531,0.066281,
c     .                        0.066922,0.0497974/)
        real*8 :: lam8(nlt)
!!!!!!!!!! end Boris' part !!!!!!!!!!!!!!!!!

      real*8 :: roair, rn
      integer :: ngiss

      call init
      rn = 1.341d0  ! index of refraction of pure seawater
      roair = 1.2D3 ! density of air g/m3  SHOULD BE INTERACTIVE?

      sunz=acos(solz)/radian  !in degs

c  Foam and diffuse reflectance
      if (wind .gt. 4.0) then
        if (wind .le. 7.0) then
          cn = 6.2D-4 + 1.56D-3/wind
          rof = roair*cn*2.2D-5*wind*wind - 4.0D-4
        else
          cn = 0.49D-3 + 0.065D-3*wind
          rof = (roair*cn*4.5D-5 - 4.0D-5)*wind*wind
        endif
        rosps = 0.057d0
      else
        rof = 0.0
        rosps = 0.066d0
      endif
      
c  Direct
c   Fresnel reflectance for sunz < 40, wind < 2 m/s
      if (sunz .lt. 40.0 .or. wind .lt. 2.0) then
        if (sunz .eq. 0.0) then
          rospd = 0.0211d0
        else
          rtheta = sunz*radian
          sintr = sin(rtheta)/rn
          rthetar = asin(sintr)
          rmin = rtheta - rthetar
          rpls = rtheta + rthetar
          sinrmin = sin(rmin)
          sinrpls = sin(rpls)
          tanrmin = tan(rmin)
          tanrpls = tan(rpls)
          sinp = (sinrmin*sinrmin)/(sinrpls*sinrpls)
          tanp = (tanrmin*tanrmin)/(tanrpls*tanrpls)
          rospd = 0.5*(sinp + tanp)
        endif
      else
       !Empirical fit otherwise
        a = 0.0253d0
        b = -7.14D-4*wind + 0.0618d0
        rospd = a*exp(b*(sunz-40.0))
      endif

c  Reflectance totals
      do nl = 1,nlt
        ros(nl) = rosps + rof*wfac(nl)
        rod(nl) = rospd + rof*wfac(nl)
      enddo

      !lam is integer, lam8 is real8
      lam8=float(lam)

      if (hycgr) return  !we do not compute albedo coefs
                         !from within ocean, but from atmos
      
!!!!!!!!!! Boris' part !!!!!!!!!!!!!!!!!
C**** get chlorophyll term
      ! function reflectance calculates reflectance 
      ! as a function of chl and wavelength (lam)

      res = reflectance(refl,chl,lam8,nlt,i,j)
 
!  transition between band33 and band6 approximation

! loop over giss radiation bands
      do ngiss=1,6
        sum1 = 0.0
        sum2 = 0.0
        do nl=gband(ngiss), gband(ngiss+1)-1 
          if (refl(nl).lt.0) refl(nl) = 0.0
          ros(nl) = ros(nl) + refl(nl)
          sum1 = sum1+weight(nl)*rod(nl)
          sum2 = sum2+weight(nl)*ros(nl)
        enddo
        part_sum=sum(weight(gband(ngiss):gband(ngiss+1)-1))
        xocvn(ngiss) = sum1/part_sum
        bocvn(ngiss) = sum2/part_sum

      if (xocvn(ngiss).ge.1. .or. bocvn(ngiss).ge.1.) then
         print*, 'XOCVN/BOCVN greater than 1 at ngiss,i,j=',
     .                ngiss,i,j
         do nl=gband(ngiss), gband(ngiss+1)-1
           write(*,'(i5,4e12.4)')
     .     nl,rod(nl),ros(nl),refl(nl),weight(gband(ngiss))
         enddo
         stop
      endif
      enddo
    
!!!!!!!!!! end Boris' part !!!!!!!!!!!!!!!!!
      return
      end subroutine ocalbedo

!=======================================================================
      subroutine init

      use filemanager, only: openunit, closeunit
      implicit none

      real*8 :: saw, sbw
      real*8 :: b0, b1, b2, b3, a0, a1, a2, a3, expterm, tlog, fac, rlam
      integer :: nl,ic , iu_bio, lambda
      character title*50
      data a0,a1,a2,a3 /0.9976d0, 0.2194d0,  5.554d-2,  6.7d-3 /
      data b0,b1,b2,b3 /5.026d0, -0.01138d0, 9.552d-6, -2.698d-9/

      if (initialized) then
        return
      else
        initialized=.true.
      endif
      call openunit('cfle1',iu_bio,.false.,.true.)
      do ic = 1,6
        read(iu_bio,'(a50)')title
      enddo
      do nl = 1,nlt
        read(iu_bio,20) lambda,saw,sbw
        lam(nl) = lambda
        aw(nl) = saw
        bw(nl) = sbw
        if (lam(nl) .lt. 900) then
          expterm = exp(-(aw(nl)+0.5*bw(nl)))
          tlog = dlog(1.0D-36+expterm)
          fac = a0 + a1*tlog + a2*tlog*tlog + a3*tlog*tlog*tlog
          wfac(nl) = max(0d0,min(fac,1d0))
        else
          rlam = float(lam(nl))
          fac = b0 + b1*rlam + b2*rlam*rlam + b3*rlam*rlam*rlam
          wfac(nl) = max(fac,0d0)
        endif
      enddo
      call closeunit(iu_bio)
 20   format(i5,f15.4,f10.4)
      return

      end subroutine init
!=======================================================================

      logical function reflectance(r2,chl_in,lam,nlt,ilon,jlat)

!  r2:  the reflectance calculated as a function of chl and wavelength
!  chl: chloryphyl
!
! compute Irradiance reflectance just below the surface (R = Eu(0-)/Ed(0-))
! described in Morel & Maritorena (2001):  Bio-optical properties of
! oceanic waters : a reappraisal. J. Geophys. Res., 106(C4) :
! 7163-7180. 
! Input : [Chl] in mg/m3 (self generated or as input), wavelengths ; ;
! Output : R = Eu(0-)/Ed(0-) or Rrs = Lu(0+)/Ed(0+) at all wavelengths between
! 350 and 700 nm with 5 nm step  
! The program computes Kd (diffuse attenuation coefficient for
! downwelling irradiance) as modeled in Morel & Maritorena (2001) and 
! uses Kw with Pope & Fry (1997) for aw, Morel (1974) 
! for bbw and Loisel & Morel (1998) for bp(550). 
!
!       wavelength is lam

! bbw_out backscattering coefficient in pure water
! bbt     backscattering coefficient due to suspended particles
! C       suspended particles (here chlorophyll) in mg/m^3
! lam     in nm
! Aw      diffuse albedo of seawater
! aw    absorption coefficient of optically pure sea water
! bw    scattering coefficient of optically pure sea water
! Kw_out  vertical diffuse attenuation coefficient for pure water
! Kd      diffuse attenuation coefficient for downward irradiance
!         (spectral attenuation for downward irradiance)



      USE MODEL_COM, only: nstep=>itime

      implicit  none
    
      integer :: nlt 
      real*8, dimension(nlt) :: lam, bbw_out, aw_out, kw_out, Kd, a, 
     .   r2, mud
      integer  :: i, j
      logical key
      integer :: ilon,jlat
  
c:       Created below.  -rjh 9/18/2008 
c:       Values interpolated/extrapolated from Morel 2001 Table 2
c:        to lam using Matlab linear interp 
c:       More changes, aromanou 3/4/2009
c: corrections regarding out-of-bounds calculations of bbw_out, kw_out, 
c: Kd and r2 (iterative procedure replaced by analytic solution)
c: r2<0 are missing values

       real*8, dimension(33) :: X = (/0.234333,0.173333,0.15300, 
     &             0.13100,0.11748,0.12086,0.10165,0.08287,
     &             0.06579,0.05072,0.04111,0.03400,0.03400,
     &             0.04000,0.04500,0.05200,0.03000,
     &             0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)

       real*8, dimension(33) :: e = (/0.99800,0.83300,0.77800,
     &             0.70000,0.64358,0.66259,0.67692,0.68657,
     &             0.68880,0.67649,0.64927,0.61500,0.62600,
     &             0.64700,0.67200,0.69700,0.60000,0.50000,
     &             0.30000,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0./)

      real*8, dimension(nlt) :: u2, bbt, bb
      real*8 :: chl_in, bp550, var_exp, c, discr
      integer :: count_loops, iterator
   

      key =.false.
   
      ! Backscattering
      bp550 = 0.416 * chl_in**0.766 ! Loisel & Morel (1998)
	 
      ! Varying backscattering exponent
      if(chl_in < 2.0 .AND. chl_in > 0 ) then
         var_exp = 0.5 * (log10(chl_in) - 0.3)
      else 
         var_exp = 0.0
      endif
	 
! Change wave references  to lam, num=nlt=33
      key =  bbw(lam,nlt,bbw_out)
      do j = 1, nlt
        !Morel and Maritorena (2001)
!       bbt(j) = 0.002 + 0.01 * (0.5 - 0.25 * dlog10(chl_in))*(lam(j)/
!    .                        550.0)**var_exp
! correction related to the Huot et al2008 paper (eqn 5)
        if (chl_in.gt.0.d0) then
             bbt(j) = (0.002 + 0.01 * (0.5 - 0.25 * log10(chl_in)))
     .              * (lam(j)/550.0d0)**var_exp
        else
             bbt(j) = 0.d0
        endif
        if (bbw_out(j) > 0.) then
           bb(j) = bbw_out(j) + bbt(j) * bp550
        else
           bb(j) = bbt(j) * bp550
        endif
      enddo
	  
      if (chl_in <= 0.001) then 
        key =  get_aw(lam,nlt,aw_out)
        key =  kw(lam,nlt,kw_out)
        key =  bbw(lam,nlt,bbw_out)
        do j = 1, nlt
           Kd(j) = kw_out(j)
           bb(j) = bbw_out(j)
           a(j) = aw_out(j)
        enddo
      else 
        key =  kw(lam,nlt,kw_out)
        do j = 1, nlt
          if (kw_out(j) > 0.) then
            Kd(j) = kw_out(j) + X(j)*chl_in**e(j) ! Morel (1988)
          else
            Kd(j) = X(j)*chl_in**e(j) ! Morel (1988)
          endif
        enddo
      endif

      ! Bilinear interpolation (chl and lambda) in the mu_d look-up table
      key = bilin_mud(chl_in,lam,nlt,mud)

      ! Loop on wavelengths
      do j = 1, nlt
        if (chl_in <= 0.001) then
          r2(j) = 0.3*bb(j)/a(j)
        else if(Kd(j).le.0.) then
          r2(j) = -1.
        else
          c = 0.33*bb(j)/Kd(j)/mud(j)
          discr = (2.25*c - 1)**2 - 4*c
          r2(j)=-1.
          if(discr >= 0d0) r2(j) = .5*(1-2.25*c - sqrt(discr))
          if(r2(j)>1) then     ! probably temporary monitoring ???
            write(0,*) 'chlorophyll albedo too large',nstep,ilon,jlat,j
            r2(j)=-1 
          end if
        end if
      enddo  

      reflectance = .true.
 
      end function reflectance
  
      logical function get_virtual_index(int_point,data_array,
     &                                   dim,ridx,vidx)

! Function returns the index of high bound of original coordinate data (ridx) 
! and virtual index (vidx, between 0 and 1.0) of the interpolated point.
! Input :
! 	int_point     : coordinate of the point for which interpolation is
!                       required (ex. wavelength)
! 	data_array    : original coordinate data array from which to interpolate

      implicit none
      integer, intent(in) :: dim
      real*8, intent(in)   :: int_point
      real*8, dimension(dim), intent(in) :: data_array
      real*8, intent(out) :: vidx
      integer, intent(out) :: ridx

      real*8 :: min, max
   
      min = find_min(data_array,dim)
      max = find_max(data_array,dim)
    
      if (int_point  < min) then
         vidx = 0.0
         ridx = 1
      else if (int_point >= max) then
         vidx = 1.0 
         ridx = size(data_array)
      else
         ridx = 1
      	 do while (data_array(ridx) <= int_point) 
               ridx = ridx + 1
         enddo
       vidx = (int_point - data_array(ridx - 1))/
     &        (data_array(ridx) - data_array(ridx - 1))
      endif
   
      get_virtual_index = .true.
   
      end function get_virtual_index

      real*8 function find_min(data_arr, size_arr)

      implicit none
  
      integer, intent(in)  :: size_arr
      real*8, dimension(size_arr), intent(in) :: data_arr
  
      integer :: i
  
      find_min = data_arr(1)
	  
      do i=2, size_arr
	   if (find_min >= data_arr(i)) then
	       find_min = data_arr(i)
	   endif
      enddo   
      end function find_min

      real*8 function find_max(data_arr, size_arr)

      implicit none
  
      integer, intent(in)  :: size_arr
      real*8, dimension(size_arr), intent(in) :: data_arr
  
      integer :: i
  
         find_max = data_arr(1)
	  
      do i=2, size_arr
   	   if (find_max < data_arr(i)) then
	       find_max = data_arr(i)
	   endif
      enddo
      end function find_max

      logical function  bbw(data_arr, size_arr, bbw_out)

      integer, intent(in) :: size_arr
      real*8, dimension(size_arr), intent(in) :: data_arr
      real*8, dimension(size_arr), intent(out) :: bbw_out

      real*8, dimension(61) :: bw_sb81 = (/0.151, 0.119, 0.0995, 0.082,  
     &    0.0685, 0.0575, 0.0485, 0.0415, 0.0353, 0.0305, 0.0262, 
     &    0.0229, 0.02,   0.0175, 0.0153, 0.0134, 0.012,  0.0106, 
     &    0.0094, 0.0084, 0.0076, 0.0068, 0.0061, 0.0055, 0.0049, 
     &    0.0045, 0.0041, 0.0037, 0.0034, 0.0031, 0.0029, 0.0026, 
     &    0.0024, 0.0022, 0.0021, 0.0019, 0.0018, 0.0017, 0.0016, 
     &    0.0015, 0.0014, 0.0013, 0.0012, 0.0011, 0.001,  0.001,  
     &    0.0008, 0.0008, 0.0007, 0.0007, 0.0007, 0.0007, 0.0006, 
     &    0.0006, 0.0006, 0.0005, 0.0005, 0.0005, 0.0004, 0.0004, 
     &    0.0004/)
	    
      integer, parameter :: wl_sb_size = 61
      real*8, dimension(wl_sb_size) :: wl_sb
   
      integer :: i, ri_sb81
      real*8 :: vi_sb81
      logical key 
      real*8 :: min, max, int_bbw
   

      ! Smith & Baker, 1981 lambdas: 200-800 nm with 10nm steps
      do i = 1, wl_sb_size
         wl_sb(i) = (i-1)*10 + 200.0  
      enddo
   
      do i = 1, size_arr
       bbw_out(i) = 0.0
      enddo
   

      min = find_min(wl_sb, wl_sb_size)
      max = find_max(wl_sb, wl_sb_size)
  
      do i = 1, size_arr
       key = get_virtual_index(data_arr(i),wl_sb, 
     &          wl_sb_size,ri_sb81,vi_sb81)
   
       if ((data_arr(i) < min) .or. (data_arr(i) >= max)) then
            int_bbw = -1.0 
       else
         int_bbw = 0.5 * interpolate(bw_sb81(ri_sb81-1), 
     &                   bw_sb81(ri_sb81), vi_sb81)
       endif
	
	bbw_out(i) = int_bbw    
	
      enddo
      bbw = .true.
      end function bbw

      real*8 function interpolate(low,high,vi)

      real*8, intent(in) ::low, high, vi
      real*8 :: temp
  
      temp = low + vi*(high - low)
      interpolate = temp
  
      end function interpolate
  
      logical function  get_aw(data_arr, size_arr, aw_out)

       integer, intent(in) :: size_arr
       real*8, dimension(size_arr), intent(in)  :: data_arr
       real*8, dimension(size_arr), intent(out) :: aw_out
   
       integer, parameter :: wl_pf_size = 153
   
       integer :: i, ri_aw
       real*8 :: vi_aw
       logical key
       real*8 :: min, max, int_aw_pf97
   
       real*8, dimension(158) :: aw_pf97 = (/0.0325,  0.0204,  0.0156,  
     &               0.0114,  
     &   0.01137, 0.01044, 0.00941, 0.00917, 0.00851, 0.00829, 0.00813,
     &   0.00775, 0.00663, 0.00579, 0.0053,  0.00503, 0.00473, 0.00452,
     &   0.00444, 0.00442, 0.00454, 0.00474, 0.00478, 0.00482, 0.00495,
     &   0.00504, 0.0053,  0.0058,  0.00635, 0.00696, 0.00751, 0.0083, 
     &   0.00922, 0.00969, 0.00962, 0.00957, 0.00979, 0.01005, 0.01011,
     &   0.0102,  0.0106,  0.0109,  0.0114,  0.0121,  0.0127,  0.0131, 
     &   0.0136,  0.0144,  0.015,   0.0162,  0.0173,  0.0191,  0.0204, 
     &   0.0228,  0.0256,  0.028,   0.0325,  0.0372,  0.0396,  0.0399, 
     &   0.0409,  0.0416,  0.0417,  0.0428,  0.0434,  0.0447,  0.0452, 
     &   0.0466,  0.0474,  0.0489,  0.0511,  0.0537,  0.0565,  0.0593, 
     &   0.0596,  0.0606,  0.0619,  0.064,   0.0642,  0.0672,  0.0695, 
     &   0.0733,  0.0772,  0.0836,  0.0896,  0.0989,  0.11,    0.122,  
     &   0.1351,  0.1516,  0.1672,  0.1925,  0.2224,  0.247,   0.2577, 
     &   0.2629,  0.2644,  0.2665,  0.2678,  0.2707,  0.2755,  0.281,  
     &   0.2834,  0.2904,  0.2916,  0.2995,  0.3012,  0.3077,  0.3108, 
     &   0.322,   0.325,   0.335,   0.34,    0.358,   0.371,   0.393,  
     &   0.41,    0.424,   0.429,   0.436,   0.439,   0.448,   0.448,  
     &   0.461,   0.465,   0.478,   0.486,   0.502,   0.516,   0.538,  
     &   0.559,   0.592,   0.624,   0.663,   0.704,   0.756,   0.827,  
     &   0.914,   1.007,   1.119,   1.231,   1.356,   1.489,   1.678,  
     &   -1.0,    -1.0,    -1.0,    -1.0,    -1.0,    -1.0,    -1.0,   
     &   -1.0,    -1.0,    -1.0,    -1.0,    -1.0,    -1.0,    -1.0/)

       real*8, dimension(wl_pf_size) :: wl_pf

       wl_pf(1)=340.0
       wl_pf(2)=350.0
       wl_pf(3)=360.0
       wl_pf(4)=370.0

  
       do i = 1, 149
          wl_pf(i+4) = (i-1)*2.5 + 380.0  ! 380.-750. with 2.5nm steps
       enddo
 
   
       min = find_min(wl_pf, wl_pf_size)
       max = find_max(wl_pf, wl_pf_size)
  
       do i = 1, size_arr
        key = get_virtual_index(data_arr(i),wl_pf,
     &                     wl_pf_size,ri_aw,vi_aw)
    
        if ((data_arr(i) < min) .or. (data_arr(i) >= max)) then
         int_aw_pf97 = -1.0 
        else
         int_aw_pf97 = interpolate(aw_pf97(ri_aw-1), 
     &                    aw_pf97(ri_aw), vi_aw)
        endif
	
	aw_out(i) = int_aw_pf97    
	
       enddo
      get_aw = .true.

      end function get_aw

      logical function  kw(data_arr, size_arr, kw_out)

       integer, intent(in) :: size_arr
       real*8, dimension(size_arr), intent(in) :: data_arr
       real*8, dimension(size_arr), intent(out) :: kw_out

	    
       real*8, dimension(size_arr) :: water_abs, water_bb
       real*8, dimension(size_arr) :: aw_out, bbw_out
    
       integer :: i
       logical key
   
   
        key =  get_aw(data_arr,size_arr,aw_out)
        key =  bbw(data_arr,size_arr,bbw_out)
   
        do i = 1, size_arr
          water_abs(i) = aw_out(i)
    	 water_bb(i) = bbw_out(I)
	 kw_out(i) = water_abs(i)+water_bb(i)
	 if (kw_out(i) < 0.0) then 
	     kw_out(i) = -1;
	 endif
        enddo
    
        kw = .true.
	
      end function kw

      logical function  bilin_mud(chlin, data_arr, size_arr, intmud)

! mu_d look-up table (Morel & Maritorena, 2001)
! Bilinear interpolation (chl and wavelength) in the mu_d look-up table

      implicit none

       integer, intent(in) :: size_arr
       real*8, dimension(size_arr), intent(in) :: data_arr
       real*8, dimension(size_arr), intent(out) :: intmud
       real*8, intent(in) :: chlin
   
       real*8, dimension(size_arr) :: all_vil
   
       integer, parameter :: wl_size = 10	
       real*8, dimension(wl_size) :: wl = (/350.0, 400.0, 412.0, 
     &                        443.0, 490.0,
     &                        510.0, 555.0, 620.0, 670.0, 700.0/)
					    
       integer, parameter :: chlut_size = 6					   
       real*8, dimension(6) :: chlut = 
     &                         (/0.03, 0.1, 0.3, 1.0, 3.0, 10.0/)
 
       real*8, dimension(6) :: ch = (/0.03, 0.1, 0.3, 1.0, 3.0, 10.0/)
   
       real*8, dimension(60) :: lutmud = 
     &                       (/0.770, 0.769, 0.766, 0.767, 0.767, 0.767,
     &                       0.770, 0.769, 0.766, 0.767, 0.767, 0.767,
     &                       0.765, 0.770, 0.774, 0.779, 0.782, 0.782,
     &                       0.800, 0.797, 0.796, 0.797, 0.799, 0.799,
     &                       0.841, 0.824, 0.808, 0.797, 0.791, 0.791,
     &                       0.872, 0.855, 0.834, 0.811, 0.796, 0.796,
     &                       0.892, 0.879, 0.858, 0.827, 0.795, 0.795,
     &                       0.911, 0.908, 0.902, 0.890, 0.871, 0.871,
     &                       0.914, 0.912, 0.909, 0.901, 0.890, 0.890,
     &                       0.914, 0.912, 0.909, 0.901, 0.890, 0.890/)
						 					 
	real*8, dimension(10,6) :: lut
	real*8    :: min, max, min_wl, max_wl, vic, vil	
	integer :: ic=0, il=0, k				 

	 integer :: row,col,i


   
  	do row=1, 10
	  do col=1, 6
	   lut(row,col)=lutmud(col+6*(row-1))
	   enddo
	enddo					 
  

      min = minval(ch)
      max = maxval(ch)
	
   
      if (chlin < min) then
	      vic = 1 
      else if (chlin >= max) then
	      vic = chlut_size
	 else
	    ic = 1
	      do while (ch(ic) <= chlin)
		   ic = ic + 1
		 enddo
          vic = (ic - 1) + (chlin - ch(ic - 1))/(ch(ic) - ch(ic - 1))
      endif
	 
   
        min_wl = minval(wl)
        max_wl = maxval(wl)
	 
   
       do k = 1, size_arr

	 if (data_arr(k) < min_wl) then
	   vil = 1.d0 
	 else if (data_arr(k) >= max_wl) then
	   vil = wl_size 
         else
	   il = 1
	    do while (wl(il) <= data_arr(k)) 
	      il = il + 1
	    enddo 
	   vil =  float(il - 1) + (data_arr(k) - wl(il - 1))
     .         /  (wl(il) - wl(il - 1))
	endif
	
	all_vil(k) = vil
	
	
      enddo
   
   
       do i=1, size_arr
       intmud(i) = bilin(lut, 6, 10, NINT(vic), NINT(all_vil(i)))
       enddo

!   do i=1, 3
!    result = bilin(arr, 4, 4, x(i), y(i))
!    write(*,*) result
!   enddo
   
   
      bilin_mud = .true.
  
      end function bilin_mud  

      real*8 function bilin(data_arr, x_size, y_size, ind_x, ind_y)

      integer, intent(in) :: x_size, y_size
      real*8, dimension(y_size,x_size), intent(in) :: data_arr
      integer,  intent(in) :: ind_x,ind_y
      integer :: low, high
      integer :: low_i, high_i, low_j, high_j
      real*8 :: R1, R2, f
   
       low_i = ind_x
       high_i = ind_x + 1
    
       low_j = ind_y
       high_j = ind_y + 1 
	
!	write(*,*) low_i,high_i,low_j,high_j
	 
       if(high_j > y_size) then
           high_j = y_size
	   low_j = y_size -1
       endif
    
       if(high_i > x_size) then
         high_i = x_size
	 low_i = x_size -1
       endif

      if(low_i < 1) then
        low_i  = 1
        high_i = 2
       endif
    
       if(low_j < 1) then
           low_j  = 1
	   high_j = 2
       endif


       R1 = (high_i-ind_x)/(high_i-low_i)*data_arr(low_j,low_i) + 
     &       (ind_x-low_i)/(high_i-low_i)*data_arr(low_j,high_i)
       R2 = (high_i-ind_x)/(high_i-low_i)*data_arr(high_j,low_i) + 
     &       (ind_x-low_i)/(high_i-low_i)*data_arr(high_j,high_i)
       f = (high_j - ind_y)/(high_j-low_j)*R1 + (ind_y-low_j)/   
     &                     (high_j-low_j)*R2
    
!    write(*,*) R1, R2
!       write(*,*) data_arr(low_j,low_i), data_arr(low_j,high_i),data_arr(high_j,high_i),data_arr(high_j,low_i)
       bilin = f

       end function bilin

      end module ocalbedo_mod
