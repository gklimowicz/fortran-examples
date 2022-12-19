! parameters.SM_dim6.f90
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
module parameters_sm_dim6
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
  real(default), public :: vev
  complex(default), dimension(2), public :: &
       gncneu, gnclep, gncup, gncdwn
  real(default), public :: tau4, tau5
  complex(default), private :: ghgaga_sm, ghgaz_sm
  complex(default), public :: ghgaga_ac, ghgaz_ac, ghzz_ac, ghww_ac
  complex(default), public :: ghgaz_u, ghzz_u, ghww_u
  complex(default), public :: c66, c6p, c6t, c6dw, c6db, c6dpw, c6dpb, &
       c6pb, c6pg, c6g, c6w
  complex(default), public :: lambdac6
  complex(default), public :: ghzz6v3, ghzz6d, ghzz6dp, ghzz6pb
  complex(default), public :: ghaz6d, ghaz6dp, ghaz6pb, dim6awwdp, &
       dim6awwdw, dim6awwgauge
  complex(default), public :: dim6awwggg, dim6cphi, dim6vev3, ghgaga6, &
       ghww6d, ghww6dp
  complex(default), public :: dim6gggg ,dim6gggcg, dim6wwzw, &
       dim6wwzdpwdw, dim6wwzdw, dim6wwzd
  complex(default), public :: adim6h4v2, adim6h4p2, adim6ahwwdpb, &
       adim6ahwwdpw, adim6ahwwdw
  complex(default), public :: adim6aawwdw, adim6aawww, &
       adim6hhwwdw, adim6hhwwdpw 
  complex(default), public :: adim6hwwzdw, adim6hwwzddpw, &
       adim6hwwzdpw, adim6hwwzdpb
  complex(default), public :: adim6ahhzd, adim6ahhzdp, adim6ahhzpb, &
       adim6azwww, adim6azwwdwdpw
  complex(default), public :: adim6wwwww, adim6wwwwdwdpw, &
       adim6wwzzw, adim6wwzzdwdpw, adim6hhaa
  complex(default), public :: adim6hhzzd, adim6hhzzdp, adim6hhzzpb, adim6hhzzt
  public :: import_from_whizard, model_update_alpha_s

contains
  subroutine import_from_whizard (par_array, scheme)
    real(default), dimension(37), intent(in) :: par_array
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
       real(default) :: fac_gh3
       real(default) :: fghgaga
       real(default) :: fghgaz
       real(default) :: c66
       real(default) :: c6p
       real(default) :: c6t
       real(default) :: c6dw
       real(default) :: c6db
       real(default) :: c6dpw
       real(default) :: c6dpb
       real(default) :: c6pb
       real(default) :: c6pg
       real(default) :: c6g
       real(default) :: c6w
       real(default) :: lambdac6
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
    par%fac_gh3 = par_array(19)
    par%fghgaga = par_array(20)
    par%fghgaz  = par_array(21)
    par%c66     = par_array(22)
    par%c6p     = par_array(23)
    par%c6t     = par_array(24)
    par%c6dw    = par_array(25)
    par%c6db    = par_array(26)
    par%c6dpw   = par_array(27)
    par%c6dpb   = par_array(28)
    par%c6pb    = par_array(29)
    par%c6pg    = par_array(30)
    par%c6g     = par_array(31)
    par%c6w     = par_array(32)
    par%lambdac6= par_array(33)
    par%v       = par_array(34)
    par%cw      = par_array(35)
    par%sw      = par_array(36)
    par%ee      = par_array(37)
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
    ghww =  mass(24) * g
    ghhww = g**2 / 2.0_default
    ghzz = mass(23) * g / costhw
    ghhzz = g**2 / 2.0_default / costhw**2
    ghtt = - mass(6) / vev
    ghbb = - mass(5) / vev
    ghcc = - mass(4) / vev
    ghtautau = - mass(15) / vev
    ghmm = - mass(13) / vev
    gh3 = - par%fac_gh3 * 3 * mass(25)**2 / vev
    gh4 = - 3 * mass(25)**2 / vev**2
    !!! Color flow basis, divide by sqrt(2)
    gs = sqrt(2.0_default*PI*par%alphas)
    igs = cmplx(0.0_default, 1.0_default, kind=default) * gs    
    c66 = par%c66  
    c6p = par%c6p  
    c6t = par%c6t  
    c6dw = par%c6dw
    c6db = par%c6db 
    c6dpw = par%c6dpw
    c6dpb = par%c6dpb
    c6pb = par%c6pb  
    c6pg = par%c6pg  
    c6g = par%c6g    
    c6w = par%c6w
    lambdac6 = par%lambdac6    
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
    ghgaga_ac = par%fghgaga * ghgaga_sm   !!! Not yet implemented
    ghgaz_ac  = par%fghgaz * ghgaz_sm     !!! Not yet implemented
    ghzz_ac = 0.0_default                 !!! Not yet implemented
    ghww_ac = 0.0_default                 !!! Not yet implemented
    ghgaz_u = 0.0_default                 !!! Not yet implemented
    ghzz_u  = 0.0_default                 !!! Not yet implemented 
    ghww_u  = 0.0_default                 !!! Not yet implemented 
    ghaz6d = - e / costhw * vev * c6dw/(lambdac6**2) + &
         g * vev * c6db/(lambdac6**2)
    ghaz6dp = 0.5_default * e / costhw * vev * c6dpw/(lambdac6**2) - &
         0.5_default * g * vev * c6dpb/(lambdac6**2)
    ghaz6pb = 4.0_default * costhw * sinthw * vev * c6pb/(lambdac6**2)
    ghgaga6 = 4.0_default * costhw * costhw * vev * c6pb/(lambdac6**2)
    dim6awwggg = 6.0_default * sinthw * c6w/(lambdac6**2)
    dim6awwdp = 0.25_default * g * costhw * vev**2 *  &
         (g * c6dpb/(lambdac6**2) + e / costhw * c6dpw/(lambdac6**2))
    dim6awwdw = 0.5_default * e**2 / sinthw * vev**2 * c6dw/(lambdac6**2)
    dim6vev3 = 15 * vev**3 * c66/(lambdac6**2)
    dim6cphi = 4.0_default * vev * c6p/(lambdac6**2)
    ghww6d = g * vev * c6dw/(lambdac6**2)
    ghww6dp = 0.5_default * g * vev * c6dpw/(lambdac6**2)
    dim6wwzw = 6 * costhw * c6w/(lambdac6**2)
    dim6wwzdpwdw = 0.25_default * g**2 / costhw * vev**2 * c6dpw/(lambdac6**2) &
         + 0.5_default * g * vev**2 * (sinthw*e/costhw + 3*costhw*g) &
         * c6dw/(lambdac6**2)
    dim6wwzdw = 0.5_default * e**2 / costhw * vev**2 * c6dw/(lambdac6**2)
    dim6wwzd = 0.25_default * g**2 * costhw * vev**2 * c6dpw/(lambdac6**2) &
         + 1.5 * costhw * g**2 * vev**2 * c6dw/(lambdac6**2) &
         - 0.25_default * g**2 * sinthw * vev**2 * c6dpb/(lambdac6**2) 
    ghzz6v3 = - 2 * g**2 / costhw**2 * vev**3 * c6t/(lambdac6**2) 
    ghzz6d = g * vev * c6dw/(lambdac6**2) + e / costhw * vev &
         * c6db/(lambdac6**2)
    ghzz6dp = 0.5_default * vev * (g * c6dpw/(lambdac6**2) + &
         (e / costhw) * c6dpb/(lambdac6**2))
    ghzz6pb = - 4.0_default * sinthw**2 * vev * c6pb/(lambdac6**2)
    ! dim6gggg = 6.0_default * c6g    
    ! dim6gggcg = 2.0_default * gs * vev**2 * c6pg
    adim6h4v2 = 45.0_default * vev**2 * c66/(lambdac6**2)
    adim6h4p2 = 4.0_default * c6p/(lambdac6**2)
    adim6ahwwdpb = 0.5_default * g**2 * costhw * vev * c6dpb/(lambdac6**2)
    adim6ahwwdpw = 0.5_default * (e**2 / sinthw) * vev * c6dpw/(lambdac6**2)
    adim6ahwwdw = (e**2 / sinthw) * vev * c6dw/(lambdac6**2)
    adim6aawwdw = (e**3 / sinthw) * vev**2 * c6dw/(lambdac6**2)
    adim6aawww =  6.0_default * e * sinthw * c6w/(lambdac6**2)
    adim6hhwwdpw = 0.5_default * g * c6dpw/(lambdac6**2) 
    adim6hhwwdw = g * c6dw/(lambdac6**2) 
    adim6hwwzdw = (e**2 / costhw) * vev *c6dw/(lambdac6**2) 
    adim6hwwzddpw = g**2 * costhw * vev * ( 3*c6dw/(lambdac6**2) + &
         0.5_default*c6dpw/(lambdac6**2) )
    adim6hwwzdpw = 0.5_default * (e**2 / costhw) * vev * c6dpw/(lambdac6**2) 
    adim6hwwzdpb = 0.5_default * (e**2 / sinthw) * vev * c6dpb/(lambdac6**2) 
    adim6ahhzd = g * c6db/(lambdac6**2)  - e/costhw * c6dw/(lambdac6**2) 
    adim6ahhzdp = 0.5_default * (e/costhw * c6dpw/(lambdac6**2) - &
         g * c6dpb/(lambdac6**2))
    adim6ahhzpb = 4.0_default * costhw * sinthw * c6pb/(lambdac6**2) 
    adim6azwww = 6.0_default * e * costhw * c6w/(lambdac6**2) 
    adim6azwwdwdpw = (e**3/costhw + costhw*(e**3)/(sinthw**2)) &
         * vev**2 * 0.25_default * c6dpw/(lambdac6**2) &
         + (e**3/costhw + 3*costhw*(e**3)/(sinthw**2))* vev**2 * 0.5_default &
         * c6dw/(lambdac6**2) 
    adim6wwwww = 6.0_default * e/sinthw * c6w/(lambdac6**2) 
    adim6wwwwdwdpw = (e**3)/(sinthw**3) * vev**2 * &
         (0.5_default * c6dpw/(lambdac6**2) + 2.0_default * c6dw/(lambdac6**2))
    adim6wwzzw = 6.0_default * costhw**2 * e/sinthw * c6w/(lambdac6**2) 
!!! Recalculated WK 2018-09-13: sign flip
    adim6wwzzdwdpw = -0.5_default * (costhw**2 * e**3 / sinthw**3 +  &
         e**3 / sinthw) * vev**2 * c6dpw/(lambdac6**2) &
         + (2.0_default * costhw**2 * e**3 / sinthw**3 + e**3 / sinthw) * &
         vev**2 *c6dw/(lambdac6**2) 
!!! Original version
!     adim6wwzzdwdpw = 0.5_default * (costhw**2 * e**3 / sinthw**3 +  &
!          e**3 / sinthw) * vev**2 * c6dpw/(lambdac6**2) &
!          + (2.0_default * costhw**2 * e**3 / sinthw**3 + e**3 / sinthw) * &
!          vev**2 *c6dw/(lambdac6**2) 
    adim6hhaa = 4.0_default * costhw**2 * c6pb/(lambdac6**2)
    adim6hhzzd = (costhw**2 * e / sinthw + e * sinthw ) * c6dw/(lambdac6**2) &
         + (costhw * e + sinthw**2 * e/costhw) *c6db/(lambdac6**2)
    adim6hhzzdp = 0.5_default*(costhw**2 * e / sinthw + e * sinthw ) * &
         c6dpw/(lambdac6**2) &
         + 0.5_default*(costhw * e + sinthw**2 * e/costhw) *c6dpb/(lambdac6**2)
    adim6hhzzpb = 4.0_default * sinthw**2 * c6pb/(lambdac6**2)
    adim6hhzzt = - 6.0_default * vev**2 * (g/costhw)**2 * c6t/(lambdac6**2)
  end subroutine import_from_whizard

  subroutine model_update_alpha_s (alpha_s)
    real(default), intent(in) :: alpha_s
    gs = sqrt(2.0_default*PI*alpha_s)
    igs = cmplx(0.0_default, 1.0_default, kind=default) * gs     
  end subroutine model_update_alpha_s
end module parameters_sm_dim6
