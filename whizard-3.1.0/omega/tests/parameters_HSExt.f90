!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

module parameters_hsext
  use kinds
  use constants
  implicit none
  private

  real(default), dimension(35), public :: mass, width
  real(default), public :: as
  complex(default), public :: gs, igs

  real(default), public :: e, g, e_em
  real(default), public :: sinthw, costhw, sin2thw, tanthw
  real(default), public :: qelep, qeup, qedwn
  real(default), public :: ttop, tbot, tch, ttau, tw
  real(default), public :: ltop, lbot, lc, ltau, lw
  complex(default), public :: qlep, qup, qdwn, gcc, qw, &
       gzww, gwww, ghww, gsww, ghhww, ghsww, gssww, &
       ghzz, ghhzz, gszz, ghszz, gsszz, &
       ghbb, ghtt, ghcc, ghtautau, gh3, gh4, ghhs, ghss, gsss, &
       gxgaga, gxgaz, gxgg, ghgaga, ghgaz, ghgg, gsgaga, gsgaz, gsgg, &
       ghmm, gsbb, gstt, gstautau, gscc, &
       gsmm, iqw, igzww, igwww, gw4, gzzww, gazww, gaaww, &
       gh4_1, gh4_2, gh4_3, gh4_4, gh4_5
  real(default), public :: vev
  real(default) :: sina, cosa
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn
  real(default), public :: xi0 = 0, xipm = 0

  public :: import_from_whizard, model_update_alpha_s, init_parameters

contains

  subroutine init_parameters
    real(default), dimension(27) :: vars
    vars(1)  = 1.16639E-5   ! Fermi constant
    vars(2)  = 91.1882      ! Z-boson mass
    vars(3)  = 80.419       ! W-boson mass
    vars(4)  = 125          ! Higgs mass
    vars(5)  = 200    ! Singlet mass
    vars(6)  = 0.1178       ! Strong coupling constant (Z point)
    vars(7)  = 0.000510997  ! electron mass
    vars(8)  = 0.105658389  ! muon mass
    vars(9)  = 1.77705      ! tau-lepton mass
    vars(10)  = 0.095        ! s-quark mass
    vars(11) = 1.2          ! c-quark mass
    vars(12) = 4.2          ! b-quark mass
    vars(13) = 171.9        ! t-quark mass
    vars(14) = 1.523        ! t-quark width
    vars(15) = 2.443        ! Z-boson width
    vars(16) = 2.049        ! W-boson width
    vars(17) = 0.004143     ! Higgs width
    vars(18) = 0.001  ! Singlet width
    vars(19) = 0.1    ! Mixing parameter for the singlet
    vars(20) = 0.1    ! Ratio of vacuum expectation values
    vars(21) = 0.100        ! anomaly Higgs couplings K factors
    vars(22) = 0.200        ! anomaly Higgs couplings K factors
    vars(23) = 0.300        ! anomaly Higgs couplings K factors
    vars(24) = 1 / sqrt (sqrt (2.) * vars(1))    ! v (Higgs vev)
    vars(25) = vars(3) / vars(2)                 ! cos(theta-W)
    vars(26) = sqrt (1-vars(25)**2)              ! sin(theta-W)
    vars(27) = 2 * vars(26) * vars(3) / vars(24) ! em-coupling (GF scheme)
    call import_from_whizard (vars)
  end subroutine init_parameters

  subroutine import_from_whizard (par_array)
    real(default), dimension(27), intent(in) :: par_array
    type :: parameter_set
       real(default) :: gf
       real(default) :: mZ
       real(default) :: mW
       real(default) :: mH
       real(default) :: msing
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
       real(default) :: wsing
       real(default) :: sinal
       real(default) :: tanbe
       real(default) :: khgaz
       real(default) :: khgaga
       real(default) :: khgg
       real(default) :: v
       real(default) :: cw
       real(default) :: sw
       real(default) :: ee
    end type parameter_set
    type(parameter_set) :: par
    !!! This corresponds to 1/alpha = 137.03598949333
    real(default), parameter :: &
         alpha = 1.0_default/137.03598949333_default
    e_em = sqrt(4.0_default * PI * alpha)
    par%gf     = par_array(1)
    par%mZ     = par_array(2)
    par%mW     = par_array(3)
    par%mH     = par_array(4)
    par%msing  = par_array(5)
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
    par%wsing  = par_array(18)
    par%sinal  = par_array(19)
    par%tanbe  = par_array(20)
    par%khgaz  = par_array(21)
    par%khgaga = par_array(22)
    par%khgg   = par_array(23)
    par%v      = par_array(24)
    par%cw     = par_array(25)
    par%sw     = par_array(26)
    par%ee     = par_array(27)
    mass(1:35) = 0
    width(1:35) = 0
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
    mass(35) = par%msing
    width(35) = par%wsing
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
    vev = par%v
    e = par%ee
    sinthw = par%sw
    sin2thw = par%sw**2
    costhw = par%cw
    tanthw = sinthw/costhw
    sina = par%sinal
    cosa = sqrt(1 - sina**2)
    qelep = - 1
    qeup = 2.0_default / 3.0_default
    qedwn = - 1.0_default / 3.0_default
    g = e / sinthw
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
    ghww = mass(24) * g * cosa
    gsww = mass(24) * g * sina
    ghhww = g**2 / 2.0_default * cosa**2
    ghsww = g**2 / 2.0_default * cosa * sina
    gssww = g**2 / 2.0_default * sina**2
    ghzz = mass(23) * g / costhw * cosa
    gszz = mass(23) * g / costhw * sina
    ghhzz = g**2 / 2.0_default / costhw**2 * cosa**2
    ghszz = g**2 / 2.0_default / costhw**2 * cosa * sina
    gsszz = g**2 / 2.0_default / costhw**2 * sina**2
    ghtt = - mass(6) / vev * cosa
    ghbb = - mass(5) / vev * cosa
    ghcc = - mass(4) / vev * cosa
    ghtautau = - mass(15) / vev * cosa
    ghmm = - mass(13) / vev * cosa
    gstt = - mass(6) / vev * sina
    gsbb = - mass(5) / vev * sina
    gscc = - mass(4) / vev * sina
    gstautau = - mass(15) / vev * sina
    gsmm = - mass(13) / vev * sina
    gh3 = - 3 * mass(25)**2 / vev * sina**3
    ghhs = - (mass(25)**2 + mass(35)**2) / (2 * vev) * &
         (sina * par%tanbe + cosa) * sina * cosa
    !!! These couplings are not checked, and probably wrong
    ghss = ghhs
    gsss = - 3 * mass(25)**2 / vev * cosa**3
    gh4 = - 3 * mass(25)**2 / vev**2
    gh4_1 = gh4
    gh4_2 = gh4
    gh4_3 = gh4
    gh4_4 = gh4
    gh4_5 = gh4
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt(2.0_default*PI*par%alphas)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Higgs anomaly couplings
    !!! SM LO loop factor (top,bottom,W)
    gxgaga = (-1._default) * alpha / vev / 2.0_default / PI * &
         (( 4.0_default * (fonehalf(ttop) + fonehalf(tch)) &
         + fonehalf(tbot)) / 3.0_default + fonehalf(ttau) + fone(tw)) &
         * sqrt(par%khgaga)
    ghgaga = gxgaga * cosa
    gsgaga = gxgaga * sina
    !!! asymptotic limit:
    !!! ghgaga = (par%ee)**2 / vev / &
    !!!      9.0_default / pi**2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SM LO loop factor (only top and W)
    gxgaz = e * e_em / 8.0_default / PI**2 / vev * ( &
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
    ghgaz = gxgaz * cosa
    gsgaz = gxgaz * sina
    !!! SM LO order loop factor with
    !!! N(N)LO K factor = 2.1 (only top)
    !!! Limit of infinite top quark mass:
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! We use par%gg because of sqrt(2) above
    gxgg = (-1._double) * par%alphas / vev / 4.0_default / PI * &
         (fonehalf(ttop) + fonehalf(tbot) + fonehalf(tch)) * &
         sqrt(par%khgg)
    ghgg = gxgg * cosa
    gsgg = gxgg * sina
    !!! ghgg   = par%alphas / 3.0_default &
    !!!      / vev / pi * 2.1_default
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs
    !!! The Hgg coupling should not get a running alpha_s
  end subroutine model_update_alpha_s

  elemental function faux (x) result (y)
    real(default), intent(in) :: x
    complex(default) :: y
    if (1 <= x) then
       y = asin(sqrt(1/x))**2
    else
       y = - 1/4.0_default * (log((1 + sqrt(1 - x))/ &
            (1 - sqrt(1 - x))) - cmplx (0.0_default, pi, kind=default))**2
    end if
  end function faux

  elemental function fone (x) result  (y)
    real(default), intent(in) :: x
    complex(default) :: y
    if (x==0) then
       y = 2.0_default
    else
       y = 2.0_default + 3.0_default * x + &
            3.0_default * x * (2.0_default - x) * &
            faux(x)
    end if
  end function fone

  elemental function fonehalf (x) result (y)
     real(default), intent(in) :: x
     complex(default) :: y
     if (x==0) then
        y = 0
     else
        y = - 2.0_default * x * (1 + (1 - x) * faux(x))
     end if
  end function fonehalf

  elemental function tri_i1 (a,b) result (y)
    real(default), intent(in) :: a,b
    complex(default) :: y
    if (a < epsilon(a) .or. b < epsilon (b)) then
       y = 0
    else
       y = a*b/2.0_default/(a-b) + a**2 * b**2/2.0_default/(a-b)**2 * &
            (faux(a) - faux(b)) + &
            a**2 * b/(a-b)**2 * (gaux(a) - gaux(b))
    end if
  end function tri_i1

  elemental function tri_i2 (a,b) result (y)
    real(default), intent(in) :: a,b
    complex(default) :: y
    if (a < epsilon (a) .or. b < epsilon(b)) then
       y = 0
    else
       y = - a * b / 2.0_default / (a-b) * (faux(a) - faux(b))
    end if
  end function tri_i2

  elemental function gaux (x) result (y)
    real(default), intent(in) :: x
    complex(default) :: y
    if (1 <= x) then
       y = sqrt(x - 1) * asin(sqrt(1/x))
    else
       y = sqrt(1 - x) * (log((1 + sqrt(1 - x)) / &
            (1 - sqrt(1 - x))) - &
            cmplx (0.0_default, pi, kind=default)) / 2.0_default
    end if
  end function gaux

end module parameters_hsext
