! parameters.SM_ac_CKM.f90
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
module parameters_sm_ac_ckm
  use kinds
  use constants 
  use sm_physics !NODEP!
  implicit none
  private

  real(default), dimension(27), public :: mass, width
  real(default), public :: as
  complex(default), public :: gs, igs

  real(default), public :: e, g, e_em
  real(default), public :: sinthw, costhw, sin2thw, tanthw
  real(default), public :: qelep, qeup, qedwn
  real(default), public :: ttop, tbot, tch, ttau, tw
  real(default), public :: ltop, lbot, lc, ltau, lw
  complex(default), public :: qlep, qup, qdwn, gcc, qw, &
       gzww, gwww, ghww, ghhww, ghzz, ghhzz, &
       ghbb, ghtt, ghcc, ghtautau, gh3, gh4, ghmm, & 		
       iqw, igzww, igwww, gw4, gzzww, gazww, gaaww
  complex(default), public :: &
       gccq11 = 0, gccq12 = 0, gccq13 = 0, gccq21 = 0, &
       gccq22 = 0, gccq23 = 0, gccq31 = 0, gccq32 = 0, gccq33 = 0     
  real(default), public :: vev
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn
  real(default), public :: a4 = 0, a5 = 0, a6 = 0, a7 = 0, a10 = 0
  real(default), public :: tau4, tau5
  complex(default), public :: ig1a, ig1z, rg5a, rg5z, &
       ig1pkpg4a, ig1pkpg4z, ig1pkmg4a, ig1pkmg4z, &
       ig1mkpg4a, ig1mkpg4z, ig1mkmg4a, ig1mkmg4z, &
       ila, ilz, il5a, il5z, ik5a, ik5z, &
       alww0, alww2, alzw0, alzw1, alzz
  complex(default), private :: ghgaga_sm, ghgaz_sm
  complex(default), public :: ghgaga_ac, ghgaz_ac, ghzz_ac, ghww_ac
  complex(default), public :: ghgaz_u, ghzz_u, ghww_u
  complex(default), public :: lambda_h, fw, fww, fb, fbb
  complex(default), private :: h_anom

  public :: import_from_whizard, model_update_alpha_s

contains
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(63), intent(in) :: par_array
    integer, intent(in) :: scheme
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
       real(default) :: xi0
       real(default) :: xipm
       real(default) :: vckm11
       real(default) :: vckm12
       real(default) :: vckm13
       real(default) :: vckm21
       real(default) :: vckm22
       real(default) :: vckm23
       real(default) :: vckm31
       real(default) :: vckm32
       real(default) :: vckm33
       real(default) :: a4
       real(default) :: a5
       real(default) :: a6
       real(default) :: a7
       real(default) :: a10
       real(default) :: g1a
       real(default) :: g1z
       real(default) :: g4a
       real(default) :: g4z
       real(default) :: g5a
       real(default) :: g5z
       real(default) :: ka
       real(default) :: kz
       real(default) :: la
       real(default) :: lz
       real(default) :: k5a
       real(default) :: k5z
       real(default) :: l5a
       real(default) :: l5z
       real(default) :: v
       real(default) :: cw
       real(default) :: sw
       real(default) :: ee
       real(default) :: csw
       real(default) :: az
       real(default) :: awz1
       real(default) :: awz2
       real(default) :: fac_gh3
       real(default) :: fac_gh4
       real(default) :: fghgaga
       real(default) :: fghgaz
       real(default) :: lambdah
       real(default) :: fw
       real(default) :: fww
       real(default) :: fb
       real(default) :: fbb
    end type parameter_set
    type(parameter_set) :: par
    !!! This corresponds to 1/alpha = 137.03598949333
    real(default), parameter :: &
         alpha = 1.0_default/137.03598949333_default
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
    par%xi0     = par_array(17)
    par%xipm    = par_array(18)
    par%vckm11  = par_array(19)
    par%vckm12  = par_array(20)
    par%vckm13  = par_array(21)
    par%vckm21  = par_array(22)
    par%vckm22  = par_array(23)
    par%vckm23  = par_array(24)
    par%vckm31  = par_array(25)
    par%vckm32  = par_array(26)
    par%vckm33  = par_array(27)
    par%a4      = par_array(28)         
    par%a5      = par_array(29)
    par%a6      = par_array(30)
    par%a7      = par_array(31)
    par%a10     = par_array(32)
    par%g1a     = par_array(33)
    par%g1z     = par_array(34)
    par%g4a     = par_array(35)
    par%g4z     = par_array(36)
    par%g5a     = par_array(37)
    par%g5z     = par_array(38)
    par%ka      = par_array(39)
    par%kz      = par_array(40)
    par%la      = par_array(41)
    par%lz      = par_array(42)
    par%k5a     = par_array(43)
    par%k5z     = par_array(44)
    par%l5a     = par_array(45)
    par%l5z     = par_array(46)
    par%fac_gh3 = par_array(47)
    par%fac_gh4 = par_array(48)
    par%fghgaga = par_array(49)
    par%fghgaz  = par_array(50)
    par%lambdah = par_array(51)
    par%fw      = par_array(52)
    par%fww     = par_array(53)
    par%fb      = par_array(54)
    par%fbb     = par_array(55)
    par%v       = par_array(56)
    par%cw      = par_array(57)
    par%sw      = par_array(58)
    par%ee      = par_array(59)
    par%csw     = par_array(60)
    par%aZ      = par_array(61)
    par%aWZ1    = par_array(62)
    par%aWZ2    = par_array(63)
    mass(1:27) = 0
    width(1:27) = 0
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
    mass(26) =  par%xi0 * mass(23)
    width(26) =  0
    mass(27) =  par%xipm * mass(24)
    width(27) =  0
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
    qelep = - 1
    qeup = 2.0_default / 3.0_default
    qedwn = - 1.0_default / 3.0_default
    g = e / sinthw
    gcc = - g / 2 / sqrt (2.0_default)
    gccq11 = gcc * par%vckm11
    gccq12 = gcc * par%vckm12
    gccq13 = gcc * par%vckm13
    gccq21 = gcc * par%vckm21
    gccq22 = gcc * par%vckm22
    gccq23 = gcc * par%vckm23
    gccq31 = gcc * par%vckm31
    gccq32 = gcc * par%vckm32
    gccq33 = gcc * par%vckm33
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
    ghmm = - mass(13) / vev
    gh3 = - par%fac_gh3 * 3 * mass(25)**2 / vev
    gh4 = - par%fac_gh4 * 3 * mass(25)**2 / vev**2
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt(2.0_default*PI*par%alphas)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs    
    a4 = par%a4
    a5 = par%a5
    a6 = par%a6
    a7 = par%a7
    a10 = par%a10
    ig1a = iqw * par%g1a
    ig1z = igzww * par%g1z
    ig1pkpg4a = iqw   * (par%g1a + par%ka + par%g4a) / 2
    ig1pkpg4z = igzww * (par%g1z + par%kz + par%g4z) / 2
    ig1pkmg4a = iqw   * (par%g1a + par%ka - par%g4a) / 2
    ig1pkmg4z = igzww * (par%g1z + par%kz - par%g4z) / 2
    ig1mkpg4a = iqw   * (par%g1a - par%ka + par%g4a) / 2
    ig1mkpg4z = igzww * (par%g1z - par%kz + par%g4z) / 2
    ig1mkmg4a = iqw   * (par%g1a - par%ka - par%g4a) / 2
    ig1mkmg4z = igzww * (par%g1z - par%kz - par%g4z) / 2
    ila = iqw   * par%la / (mass(24)*mass(24))
    ilz = igzww * par%lz / (mass(24)*mass(24))
    rg5a = qw   * par%g5a
    rg5z = gzww * par%g5z
    ik5a = iqw   * par%k5a
    ik5z = igzww * par%k5z
    il5a = iqw   * par%l5a / (mass(24)*mass(24))
    il5z = igzww * par%l5z / (mass(24)*mass(24))
    alww0 = g**4 * (a4 + 2 * a5)
    alww2 = g**4 * 2 * a4
    alzw1 = g**4 / costhw**2 * (a4 + a6)
    alzw0 = g**4 / costhw**2 * 2 * (a5 + a7)
    alzz = g**4 / costhw**4 * 2 * (a4 + a5 + (a6+a7+a10)*2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Higgs anomaly couplings
    !!! SM LO loop factor (top,bottom,W)
    ghgaga_sm = (-1._default) * alpha / vev / 2.0_default / PI * &
         (( 4.0_default * (fonehalf(ttop) + fonehalf(tch)) &
         + fonehalf(tbot)) / 3.0_default + fonehalf(ttau) + fone(tw)) 
    !!! asymptotic limit:
    !!! ghgaga_sm = (par%ee)**2 / vev / &
    !!!      9.0_default / pi**2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! SM LO loop factor (only top and W)
    ghgaz_sm = e * e_em / 8.0_default / PI**2 / vev * ( &
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
          )/sinthw
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    h_anom = g * mass(25) / par%lambdah**2
    ghgaga_ac = par%fghgaga * ghgaga_sm - h_anom * &
         sin2thw * (par%fbb + par%fww) / 2.0_default
    ghgaz_ac  = par%fghgaz * ghgaz_sm + h_anom * &
         sinthw * (sinthw**2 * par%fbb - costhw**2 * par%fww) / costhw
    ghzz_ac = - h_anom * (sinthw**4 * par%fbb + costhw**4 * par%fww) / &
         2.0_default / costhw**2
    ghww_ac = - h_anom * par%fww
    ghgaz_u = h_anom * sinthw * (par%fw - par%fb) / 2.0_default / costhw
    ghzz_u = h_anom * (costhw**2 * par%fw + sinthw**2 * par%fb) / &
         2.0_default / costhw**2
    ghww_u = h_anom * par%fw / 2.0_default
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx (0.0_default, 1.0_default, kind=default) * gs     
  end subroutine model_update_alpha_s
end module parameters_sm_ac_ckm
