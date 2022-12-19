! parameters.UED.f90
!
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module parameters_ued
  use kinds
  use constants
  use sm_physics !NODEP!
  implicit none
  private

  real(default), dimension(39), public :: mass, width
  real(default), public :: e, g, gp, sinthw, costhw, sin2thw, tanthw, e_em
  complex(default), public :: gs, igs
  real(default), public :: qelep, qeup, qedwn
  real(default), public :: ttop, tbot, tch, ttau, tw
  real(default), public :: ltop, lbot, lc, ltau, lw
  complex(default), public :: qlep, qup, qdwn, gcc, qw, &
       gzww, gwww, ghww, ghhww, ghzz, ghhzz, &
       ghbb, ghtt, ghcc, ghtautau, gh3, gh4, &
       ghgaga, ghgaz, ghgg, ghmm, & 		
       iqw, igzww, igwww, gw4, gzzww, gazww, gaaww
  complex(default), public :: ggrav, iqwk
  real(default), public :: vev
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn
  real(default), public :: xi0 = 0, xipm = 0
  
  integer, parameter, public :: &
       n0 = 5, nloop = 2 
  real(default), parameter :: &
       acc = 1.e-12_default
  real(default), parameter :: &
       asmz = 0.118_default


  public :: import_from_whizard, model_update_alpha_s

contains

  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(134), intent(in) :: par_array
    integer, intent(in) :: scheme
    !!! This corresponds to 1/alpha = 137.03598949333
    real(default), parameter :: &
         alpha = 1.0_default/137.03598949333_default
    real(default), parameter :: &
         asmz = 0.118_default
    type :: parameter_set
       real(default) :: gf
       real(default) :: mZ
       real(default) :: mW
       real(default) :: mH
       real(default) :: alphas
       real(default) :: me
       real(default) :: mmu
       real(default) :: mtau
       real(default) :: ms
       real(default) :: mc
       real(default) :: mb
       real(default) :: mtop
       real(default) :: wtop
       real(default) :: wZ
       real(default) :: wW
       real(default) :: wH
       real(default) :: mGG
       real(default) :: wGG
       real(default) :: khgaz
       real(default) :: khgaga
       real(default) :: khgg
       real(default) :: ggrav
       real(default) :: mgk1
       real(default) :: wgk1
       real(default) :: mgk2
       real(default) :: wgk2
       real(default) :: mb1
       real(default) :: wb1
       real(default) :: mb2
       real(default) :: wb2
       real(default) :: mz1
       real(default) :: wz1
       real(default) :: mz2
       real(default) :: wz2
       real(default) :: mw1
       real(default) :: ww1
       real(default) :: mw2
       real(default) :: ww2
       real(default) :: mek1l
       real(default) :: mmuk1l
       real(default) :: mtauk1l
       real(default) :: mek1r
       real(default) :: mmuk1r
       real(default) :: mtauk1r
       real(default) :: mnuek1
       real(default) :: mnumuk1
       real(default) :: mnutauk1
       real(default) :: mek2l
       real(default) :: mmuk2l
       real(default) :: mtauk2l
       real(default) :: mek2r
       real(default) :: mmuk2r
       real(default) :: mtauk2r
       real(default) :: mnuek2
       real(default) :: mnumuk2
       real(default) :: mnutauk2
       real(default) :: wek1l
       real(default) :: wmuk1l
       real(default) :: wtauk1l
       real(default) :: wek1r
       real(default) :: wmuk1r
       real(default) :: wtauk1r
       real(default) :: wnuek1
       real(default) :: wnumuk1
       real(default) :: wnutauk1
       real(default) :: wek2l
       real(default) :: wmuk2l
       real(default) :: wtauk2l
       real(default) :: wek2r
       real(default) :: wmuk2r
       real(default) :: wtauk2r
       real(default) :: wnuek2
       real(default) :: wnumuk2
       real(default) :: wnutauk2
       real(default) :: muk1l
       real(default) :: mdk1l
       real(default) :: mck1l
       real(default) :: msk1l
       real(default) :: mtk1l
       real(default) :: mbk1l
       real(default) :: muk1r
       real(default) :: mdk1r
       real(default) :: mck1r
       real(default) :: msk1r
       real(default) :: mtk1r
       real(default) :: mbk1r
       real(default) :: muk2l
       real(default) :: mdk2l
       real(default) :: mck2l
       real(default) :: msk2l
       real(default) :: mtk2l
       real(default) :: mbk2l
       real(default) :: muk2r
       real(default) :: mdk2r
       real(default) :: mck2r
       real(default) :: msk2r
       real(default) :: mtk2r
       real(default) :: mbk2r
       real(default) :: wuk1l
       real(default) :: wdk1l
       real(default) :: wck1l
       real(default) :: wsk1l
       real(default) :: wtk1l
       real(default) :: wbk1l
       real(default) :: wuk1r
       real(default) :: wdk1r
       real(default) :: wck1r
       real(default) :: wsk1r
       real(default) :: wtk1r
       real(default) :: wbk1r
       real(default) :: wuk2l
       real(default) :: wdk2l
       real(default) :: wck2l
       real(default) :: wsk2l
       real(default) :: wtk2l
       real(default) :: wbk2l
       real(default) :: wuk2r
       real(default) :: wdk2r
       real(default) :: wck2r
       real(default) :: wsk2r
       real(default) :: wtk2r
       real(default) :: wbk2r
       real(default) :: mh1u
       real(default) :: wh1u
       real(default) :: mh1d
       real(default) :: wh1d
       real(default) :: mh2u
       real(default) :: wh2u
       real(default) :: mh2d
       real(default) :: wh2d
       real(default) :: v
       real(default) :: cw
       real(default) :: sw
       real(default) :: ee
    end type parameter_set
    type(parameter_set) :: par
    e_em = sqrt(4.0_default * PI * alpha)
    par%gf      = par_array(1)
    par%mZ      = par_array(2)
    par%mW      = par_array(3)
    par%mH      = par_array(4)
    par%alphas  = par_array(5)
    par%me      = par_array(6)
    par%mmu     = par_array(7)
    par%mtau    = par_array(8)
    par%ms      = par_array(9)
    par%mc      = par_array(10)
    par%mb      = par_array(11)
    par%mtop    = par_array(12)
    par%wtop    = par_array(13)
    par%wZ      = par_array(14)
    par%wW      = par_array(15)
    par%wH      = par_array(16)
    par%mGG     = par_array(17)
    par%wGG     = par_array(18)
    par%khgaz   = par_array(19)
    par%khgaga  = par_array(20)
    par%khgg    = par_array(21)
    par%ggrav   = par_array(22)
    par%mgk1    = par_array(23)
    par%wgk1    = par_array(24)
    par%mgk2    = par_array(25)
    par%wgk2    = par_array(26)
    par%mb1     = par_array(27)
    par%wb1     = par_array(28)
    par%mb2     = par_array(29)
    par%wb2     = par_array(30)
    par%mz1     = par_array(31)
    par%wz1     = par_array(32)
    par%mz2     = par_array(33)
    par%wz2     = par_array(34)
    par%mw1     = par_array(35)
    par%ww1     = par_array(36)
    par%mw2     = par_array(37)
    par%ww2     = par_array(38)
    par%mek1l   = par_array(39)
    par%mmuk1l  = par_array(40)
    par%mtauk1l = par_array(41)
    par%mek1r   = par_array(42)
    par%mmuk1r  = par_array(43)
    par%mtauk1r = par_array(44)
    par%mnuek1  = par_array(45)
    par%mnumuk1 = par_array(46)
    par%mnutauk1= par_array(47)
    par%mek2l   = par_array(48)
    par%mmuk2l  = par_array(49)
    par%mtauk2l = par_array(50)
    par%mek2r   = par_array(51)
    par%mmuk2r  = par_array(52)
    par%mtauk2r = par_array(53)
    par%mnuek2  = par_array(54)
    par%mnumuk2 = par_array(55)
    par%mnutauk2= par_array(56)
    par%wek1l   = par_array(57)
    par%wmuk1l  = par_array(58)
    par%wtauk1l = par_array(59)
    par%wek1r   = par_array(60)
    par%wmuk1r  = par_array(61)
    par%wtauk1r = par_array(62)
    par%wnuek1  = par_array(63)
    par%wnumuk1 = par_array(64)
    par%wnutauk1= par_array(65)
    par%wek2l   = par_array(66)
    par%wmuk2l  = par_array(67)
    par%wtauk2l = par_array(68)
    par%wek2r   = par_array(69)
    par%wmuk2r  = par_array(70)
    par%wtauk2r = par_array(71)
    par%wnuek2  = par_array(72)
    par%wnumuk2 = par_array(73)
    par%wnutauk2= par_array(74)
    par%muk1l   = par_array(75)
    par%mdk1l   = par_array(76)
    par%mck1l   = par_array(77)
    par%msk1l   = par_array(78)
    par%mtk1l   = par_array(79)
    par%mbk1l   = par_array(80)
    par%muk1r   = par_array(81)
    par%mdk1r   = par_array(82)
    par%mck1r   = par_array(83)
    par%msk1r   = par_array(84)
    par%mtk1r   = par_array(85)
    par%mbk1r   = par_array(86)
    par%muk2l   = par_array(87)
    par%mdk2l   = par_array(88)
    par%mck2l   = par_array(89)
    par%msk2l   = par_array(90)
    par%mtk2l   = par_array(91)
    par%mbk2l   = par_array(92)
    par%muk2r   = par_array(93)
    par%mdk2r   = par_array(94)
    par%mck2r   = par_array(95)
    par%msk2r   = par_array(96)
    par%mtk2r   = par_array(97)
    par%mbk2r   = par_array(98)
    par%wuk1l   = par_array(90)
    par%wdk1l   = par_array(100)
    par%wck1l   = par_array(101)
    par%wsk1l   = par_array(102)
    par%wtk1l   = par_array(103)
    par%wbk1l   = par_array(104)
    par%wuk1r   = par_array(105)
    par%wdk1r   = par_array(106)
    par%wck1r   = par_array(107)
    par%wsk1r   = par_array(108)
    par%wtk1r   = par_array(109)
    par%wbk1r   = par_array(110)
    par%wuk2l   = par_array(111)
    par%wdk2l   = par_array(112)
    par%wck2l   = par_array(113)
    par%wsk2l   = par_array(114)
    par%wtk2l   = par_array(115)
    par%wbk2l   = par_array(116)
    par%wuk2r   = par_array(117)
    par%wdk2r   = par_array(118)
    par%wck2r   = par_array(119)
    par%wsk2r   = par_array(120)
    par%wtk2r   = par_array(121)
    par%wbk2r   = par_array(122)
    par%mh1u    = par_array(123)
    par%wh1u    = par_array(124)
    par%mh1d    = par_array(125)
    par%wh1d    = par_array(126)
    par%mh2u    = par_array(127)
    par%wh2u    = par_array(128)
    par%mh2d    = par_array(129)
    par%wh2d    = par_array(130)
    par%v       = par_array(131)
    par%cw      = par_array(132)
    par%sw      = par_array(133)
    par%ee      = par_array(134)
    mass(1:39) = 0
    width(1:39) = 0
    mass(3) = par%ms
    mass(4) = par%mc
    mass(5) = par%mb
    mass(6) = par%mtop
    width(6) = par%wtop
    mass(11) = par%me
    mass(13) = par%mmu
    mass(15) = par%mtau
    mass(23) = par%mZ
    width(23) = par%wZ
    mass(24) = par%mW
    width(24) = par%wW
    mass(25) = par%mH
    width(25) = par%wH
    mass(26) =  xi0 * mass(23)
    width(26) =  0
    mass(27) =  xipm * mass(24)
    width(27) =  0
    mass(39) = par%mgg
    width(39) = par%wgg
    ttop = 4.0_default * mass(6)**2 / mass(25)**2
    tbot = 4.0_default * mass(5)**2 / mass(25)**2
    tch  = 4.0_default * mass(4)**2 / mass(25)**2
    ttau = 4.0_default * mass(15)**2 / mass(25)**2
    tw   = 4.0_default * mass(24)**2 / mass(25)**2  
    ltop = 4.0_default * mass(6)**2 / mass(23)**2
    lbot = 4.0_default * mass(5)**2 / mass(23)**2  
    lc   = 4.0_default * mass(4)**2 / mass(23)**2
    ltau = 4.0_default * mass(15)**2 / mass(23)**2
    lw   = 4.0_default * mass(24)**2 / mass(23)**2
    vev = 2 * mass(24) * par%sw / par%ee 
    e = par%ee
    sinthw = par%sw
    sin2thw = sinthw**2
    costhw = par%cw
    tanthw = sinthw/costhw
    qelep = - 1
    qeup = 2.0_default / 3.0_default
    qedwn = - 1.0_default / 3.0_default
    g = e / sinthw
    gp = e / costhw
    gcc = - g / 2 / sqrt (2.0_default)
    gncneu(1) = - g / 2 / costhw * ( + 0.5_default)
    gnclep(1) = - g / 2 / costhw * ( - 0.5_default - 2 * qelep * sin2thw)
    gncup(1)  = - g / 2 / costhw * ( + 0.5_default - 2 * qeup  * sin2thw)
    gncdwn(1) = - g / 2 / costhw * ( - 0.5_default - 2 * qedwn * sin2thw)
    gncneu(2) = - g / 2 / costhw * ( + 0.5_default)
    gnclep(2) = - g / 2 / costhw * ( - 0.5_default)
    gncup(2)  = - g / 2 / costhw * ( + 0.5_default)
    gncdwn(2) = - g / 2 / costhw * ( - 0.5_default)
    qlep = - e * qelep
    qup = - e * qeup
    qdwn = - e * qedwn
    qw = e
    iqw = (0,1)*qw
    gzww = g * costhw
    igzww = (0,1)*gzww
    gwww = g
    igwww = (0,1)*gwww
    gw4 = gwww**2
    gzzww = gzww**2
    gazww = gzww * qw
    gaaww = qw**2
    ghww = mass(24) * g
    ghhww = g**2 / 2.0_default
    ghzz = mass(23) * g / costhw
    ghhzz = g**2 / 2.0_default / costhw**2
    ghtt = - mass(6) / vev
    ghbb = - mass(5) / vev
    ghcc = - mass(4) / vev
    ghtautau = - mass(15) / vev
    gh3 = - 3 * mass(25)**2 / vev
    gh4 = - 3 * mass(25)**2 / vev**2
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt(2.0_default * PI * par%alphas)
    igs = cmplx(0.0_default, 1.0_default, kind=default) * gs    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Higgs anomaly couplings
    !!! SM LO loop factor (top,bottom,W)
    ghgaga = (-1._default) * alpha / vev / 2.0_default / PI * &
         (( 4.0_default * (fonehalf(ttop) + fonehalf(tch)) &
         + fonehalf(tbot)) / 3.0_default + fonehalf(ttau) + fone(tw)) &
         * sqrt(par%khgaga)
    !!! asymptotic limit:
    !!! ghgaga = (par%ee)**2 / vev / &
    !!!      9.0_default / pi**2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SM LO loop factor (only top and W)
    ghgaz = e * e_em / 8.0_default / PI**2 / vev * ( &
          ( - 2.0_default + &
          16.0_default/3.0_default * sin2thw) * &
          (tri_i1(ttop,ltop) - tri_i2(ttop,ltop)) / costhw & 
          + ( - 1.0_default + &
          4.0_default/3.0_default * sin2thw) & 
          * (tri_i1(tbot,lbot) - tri_i2(tbot,lbot)) / costhw &
          + (-1.0_default + 4.0_default * sin2thw) &
          * (tri_i1(ttau,ltau) - tri_i2(ttau,ltau)) / costhw &
           - costhw * ( 4.0_default * (3.0_default - tanthw**2) * &
           tri_i2(tw,lw) + ((1 + 2.0_default/tw) * tanthw**2 - ( &
           5.0_default + 2.0_default/tw)) * tri_i1(tw,lw)) &
          )/sinthw * sqrt(par%khgaz)
    !!! SM LO order loop factor with 
    !!! N(N)LO K factor = 2.1 (only top)
    !!! Limit of infinite top quark mass:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! We use par%gg because of sqrt(2) above
    ghgg = (-1._double) * par%alphas / vev / 4.0_default / PI * &
         (fonehalf(ttop) + fonehalf(tbot) + fonehalf(tch)) * &
         sqrt(par%khgg)
    !!! ghgg   = par%alphas / 3.0_default &
    !!!      / vev / pi * 2.1_default
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!   GRAVITATIONAL COUPLING      !!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ggrav = par%ggrav
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!   SPECIAL UED COUPLINGS       !!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    iqwk = g*gp / sqrt(g**2 + gp**2)
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx(0.0_default, 1.0_default, kind=default) * gs     
    !!! The Hgg coupling should not get a running alpha_s
  end subroutine model_update_alpha_s
end module parameters_ued
