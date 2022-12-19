! parameters.Simplest_univ.f90
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
module parameters_simplest_univ
  use kinds
  use constants
  use sm_physics !NODEP!
  implicit none
  private

  real(default), dimension(45), public :: mass, width
  real(default), public :: as
  complex(default), public :: gs, igs
  real(default), public :: e, g, sinthw, costhw, sin2thw, tanthw, e_em
  real(default), public :: qelep, qeup, qedwn
  real(default), public :: ttop, tbot, tch, ttau, tw
  real(default), public :: tqh1, tqh2, tqh3
  real(default), public :: ltop, lbot, lc, ltau, lw
  complex(default), public :: qlep, qup, qdwn, gcc, qw, &
       gzww, gwww, ghww, ghhww, ghzz, ghhzz, &
       ghbb, ghtt, ghcc, ghtautau, gh3, gh4, &
       ghgaga, ghgaz, ghgg, ghmm, & 		
       iqw, igzww, igwww, gw4, gzzww, gazww, gaaww
  real(default), public :: vev
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn
  real(default), parameter, public :: xi0 = 0.0, xipm = 0.0

  complex(default), public :: &
       gh0ww, gh0zz, &
       gh0tt, gh0bb, gh0cc, gh0tautau, gh0mumu, &
       iga0tt, iga0bb, iga0cc, iga0tautau, iga0mumu, &
       gahh, gzhh, igzha, igzh0a
  complex(default), dimension(2), public :: &
       ghptb, ghpcs, ghptaunu, ghpmunu
  !!! Specific simple group parameters
  complex(default), public :: &
       gncx, gncxt, gncy, gncyt
  complex(default), dimension(2), public :: &
       gnchlep, gnchneu, gnchup, gnchdwn, &
       gnchtop, gnchbot, gnchn, &
       gnchu, gnchd
  complex(default), public :: &
       ghhzzh, iqzh, igz1, igz2, igz3, & 
       igz4, igz5, igz6, i_gcc, gnch
  complex(default), public :: &
        ghtht, ghhthth, gncht, ghqhq, getht
  complex(default), public :: &
        gzeh, gzheh, gebb, gett
  complex(default), public :: & 
       ghyhvv, ghyhww
  integer, parameter, public :: &
       n0 = 5, nloop = 2 
  real(default), parameter :: &
       acc = 1.e-12_default
  real(default), parameter :: &
       asmz = 0.118_default

  public :: import_from_whizard, model_update_alpha_s

contains

  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(48), intent(in) :: par_array
    integer, intent(in) :: scheme
    !!! This corresponds to 1/alpha = 137.03598949333
    real(default), parameter :: &
         alpha = 1.0_default/137.03598949333_default
    real(default), parameter :: &
         asmz = 0.118_default
    real(default), parameter :: & 
         one = 1.0_default, two = 2.0_default, three = 3.0_default, &
         four = 4.0_default, five = 5.0_default
    complex(default), parameter :: &
         imago = (0.0_default, 1.0_default)
    real(default) :: f, d_nu, d_top, tanb, sinb, cosb, &
         xlam, xlamp, gnzh
    real(default) :: n1, n2, n3, gethth, gegg
    type :: parameter_set
       real(default) :: gf
       real(default) :: mZ
       real(default) :: mW
       real(default) :: mH
       real(default) :: meta
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
       real(default) :: weta
       real(default) :: khgaz
       real(default) :: khgaga
       real(default) :: khgg
       real(default) :: mtoph
       real(default) :: muh
       real(default) :: mch
       real(default) :: mnh1
       real(default) :: mnh2
       real(default) :: mnh3
       real(default) :: mzh
       real(default) :: mwh
       real(default) :: mxh
       real(default) :: myh
       real(default) :: wtoph
       real(default) :: wuh
       real(default) :: wch
       real(default) :: wnh1
       real(default) :: wnh2
       real(default) :: wnh3
       real(default) :: wzh
       real(default) :: wwh
       real(default) :: wxh
       real(default) :: wyh
       real(default) :: f
       real(default) :: tanb
       real(default) :: xlam
       real(default) :: v
       real(default) :: cw
       real(default) :: sw
       real(default) :: ee
    end type parameter_set
    type(parameter_set) :: par
    par%gf     = par_array(1)
    par%mZ     = par_array(2)
    par%mW     = par_array(3)
    par%mH     = par_array(4)
    par%meta   = par_array(5)
    par%alphas = par_array(6)
    par%me     = par_array(7)
    par%mmu    = par_array(8)
    par%mtau   = par_array(9)
    par%ms     = par_array(10)
    par%mc     = par_array(11)
    par%mb     = par_array(12)
    par%mtop   = par_array(13)
    par%wtop   = par_array(14)
    par%wZ     = par_array(15)
    par%wW     = par_array(16)
    par%wH     = par_array(17)
    par%weta   = par_array(18)
    par%khgaz  = par_array(19)
    par%khgaga = par_array(20)
    par%khgg   = par_array(21)
    par%mtoph  = par_array(22)
    par%muh    = par_array(23)
    par%mch    = par_array(24)
    par%mnh1   = par_array(25)
    par%mnh2   = par_array(26)
    par%mnh3   = par_array(27)
    par%mzh    = par_array(28)
    par%mwh    = par_array(29)
    par%mxh    = par_array(30)
    par%myh    = par_array(31)
    par%wtoph  = par_array(32)
    par%wuh    = par_array(33)
    par%wch    = par_array(34)
    par%wnh1   = par_array(35)
    par%wnh2   = par_array(36)
    par%wnh3   = par_array(37)
    par%wzh    = par_array(38)
    par%wwh    = par_array(39)
    par%wxh    = par_array(40)
    par%wyh    = par_array(41)
    par%f      = par_array(42)
    par%tanb   = par_array(43)
    par%xlam   = par_array(44)
    par%v      = par_array(45)
    par%cw     = par_array(46)
    par%sw     = par_array(47)
    par%ee     = par_array(48)
    e_em = sqrt(4.0_default * PI * par%alphas)
    mass(1:45) = 0
    width(1:45) = 0
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
    mass(36) = par%meta
    width(36) = par%weta
    mass(26) =  xi0 * mass(23)
    width(26) =  0
    mass(27) =  xipm * mass(24)
    width(27) =  0
    mass(32) = par%mzh
    width(32) = par%wzh
    mass(33) = par%mxh
    width(33) = par%wxh
    mass(34) = par%mwh
    width(34) = par%wwh
    mass(38) = par%myh
    width(38) = par%wyh
    mass(40) = par%muh
    width(40) = par%wuh
    mass(41) = par%mnh1
    width(41) = par%wnh1
    mass(42) = par%mch
    width(42) = par%wch
    mass(43) = par%mnh2
    width(43) = par%wnh2
    mass(44) = par%mtoph
    width(44) = par%wtoph
    mass(45) = par%mnh3
    width(45) = par%wnh3
    f = par%f
    tanb = par%tanb
    sinb = (tanb / sqrt ((1.0_default + (tanb**2))))
    cosb = (1.0_default / sqrt ((1.0_default + (tanb**2))))
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
    tqh1 = 4.0_default * mass(40)**2 / mass(25)**2
    tqh2 = 4.0_default * mass(42)**2 / mass(25)**2
    tqh3 = 4.0_default * mass(44)**2 / mass(25)**2
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
    gcc = - g / 2 / sqrt (2.0_double)
    i_gcc = imago * gcc
    n1 = sinb**2 - cosb**2
    n2 = tanb - 1.0_default/tanb
    n3 = 4.0_default * mass(6)**2/vev**2 * (sinb**4 + cosb**4)
    gncneu(1) = - g / 2 / costhw * ( + 0.5_double)
    gnclep(1) = - g / 2 / costhw * ( - 0.5_double - 2 * qelep * sin2thw)
    gncup(1)  = - g / 2 / costhw * ( + 0.5_double - 2 * qeup  * sin2thw)
    gncdwn(1) = - g / 2 / costhw * ( - 0.5_double - 2 * qedwn * sin2thw)
    gncneu(2) = - g / 2 / costhw * ( + 0.5_double)
    gnclep(2) = - g / 2 / costhw * ( - 0.5_double)
    gncup(2)  = - g / 2 / costhw * ( + 0.5_double)
    gncdwn(2) = - g / 2 / costhw * ( - 0.5_double)
    qlep = - e * qelep
    qup = - e * qeup
    qdwn = - e * qedwn
    qw = e
    iqw = imago*qw
    gzww = g * costhw
    igzww = imago*gzww
    gwww = g
    igwww = imago*gwww
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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Specific Simplest Little Higgs parameters
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !!! We are always considering the case that
    !!! lambda1/lambda2 = f2/f1 minimizing the mtoph !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Yukawa couplings
    !!!
    ghtht = (cosb**2 - sinb**2) * ghtt / two
    ghhthth = ghtt**2 * (one + (cosb - sinb)**2) 
    ghqhq = - mass(40)/two/sqrt(two)/f/tanb
    gett = n2*mass(6)/sqrt(two)/F - mass(6)**2*n1/vev/mass(44)
    getht = ghtt / two
    gebb =  - imago * mass(5) * (tanb - one/tanb) / f / sqrt(two)
    gethth = n1*mass(6)/vev + n1*n3*vev*mass(6)/two/mass(44)**2
    gh4 = - 3 * mass(25)**2 / vev**2
    !!! Color flow basis, divided by sqrt(2)
    gs = sqrt(2.0_default * PI * par%alphas)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs    
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
    gegg = (par%ee)**2 / 16.0_default / PI**2 * &
         Abs( fonehalf_pseudo(four*mass(6)**2/mass(36)**2) * gett/mass(6) &
         + fonehalf_pseudo(four*mass(44)**2/mass(36)**2) * gethth/mass(44)  &
         + fonehalf_pseudo(four*mass(40)**2/mass(36)**2) / tanb/F/sqrt(two) &
         + fonehalf_pseudo(four*mass(42)**2/mass(36)**2) / tanb/F/sqrt(two))
    xlam = mass(6) / vev * (cosb**2 - sinb**2)
    xlamp = mass(6)**2 / vev**2 * (one + (cosb**2 - sinb**2)**2)
    d_nu = - vev / sqrt(two) / f / tanb 
    d_top = xlam * vev / mass(44)
    gnzh = g / costhw / sqrt(three - four * sin2thw) / two
    !!!
    gncx = gcc * d_nu
    gncxt = - gcc * d_top        
    gncy = imago * gncx
    gncyt = imago * gncxt
    gnchneu(1) = gnzh * ( one/two + qelep * sin2thw)
    gnchlep(1) = gnchneu(1)
    gnchup(1)  = gnzh * ( - one/two + qeup * sin2thw)
    gnchdwn(1) = gnchup(1)
    gnchtop(1) = gnzh * ( one/two + qedwn * sin2thw)
    gnchbot(1) = gnchtop(1)
    gnchneu(2) = 0
    gnchlep(2) = gnzh * qelep * sin2thw
    gnchup(2)  = gnzh * qeup * sin2thw
    gnchdwn(2) = gnzh * qedwn * sin2thw
    gnchtop(2) = gnchup(2)
    gnchbot(2) =  gnchdwn(2)
    gnch = g / four / costhw * d_nu
    gncht = g / four / costhw * d_top
    !!!
    gnchn(1) = - gnzh * ( - one + sin2thw)
    gnchu(1) = - gnzh * ( - one + five / three * sin2thw)
    gnchd(1) = - gnzh * ( - one + five / three * sin2thw)
    gnchn(2) = 0
    gnchu(2) = - gnzh * qeup * sin2thw
    gnchd(2) = - gnzh * qedwn * sin2thw
    !!! 
    ghhzzh = g**2 * (one - sin2thw/costhw**2) / four / costhw / &
         sqrt(three - sin2thw/costhw**2)
    !!!
    iqzh = g * (one - two * sin2thw) / two / costhw
    igz1 = g / four / costhw
    igz2 = g / four * sqrt( one - sin2thw / three / costhw**2) * sqrt(three)
    igz3 = imago * g / two
    igz4 = g / two
    igz5 = imago * igz2
    igz6 = - imago * g * (one - sin2thw/costhw**2) * sqrt(three - sin2thw / &
         costhw**2) * vev**2 / 8.0_default / f**2
    !!! 
    gzeh = mass(23) / sqrt(two) / f * (tanb - one/tanb)
    gzheh = gzeh * costhw * sqrt(one - three * tanthw**2)
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default * PI * alpha_s)
    igs = cmplx (0.0_default, 1.0_default,kind=default) * gs   
  end subroutine model_update_alpha_s
end module parameters_simplest_univ
