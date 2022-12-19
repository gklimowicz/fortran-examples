! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
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
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'

module sm_physics

  use kinds, only: default, double
  use constants
  use physics_defs
  use lorentz

  implicit none
  private

  public :: zeta2, zeta3, zeta4, zeta5
  public :: eulerc
  public :: beta0, beta1, beta2
  public :: coeff_b0, coeff_b1, coeff_b2, coeffqed_b0, coeffqed_b1
  public :: running_as, running_as_lam, running_alpha, running_alpha_num
  public :: lambda_qcd
  public :: gamma_g, k_g
  public :: Li2
  public :: psic
  public :: psir
  public :: psim
  public :: psimr
  public :: cnielsen
  public :: nielsen
  public :: polylog
  public :: dilog
  public :: trilog
  public :: faux
  public :: fonehalf
  public :: fonehalf_pseudo
  public :: fone
  public :: gaux
  public :: tri_i1
  public :: tri_i2
  public :: run_b0
  public :: run_b1
  public :: run_aa
  public :: ff_dipole
  public :: fi_dipole
  public :: if_dipole
  public :: ii_dipole
  public :: delta
  public :: plus_distr
  public :: pqq
  public :: pgq
  public :: pqg
  public :: pgg
  public :: pqq_reg
  public :: pgg_reg
  public :: kbarqg
  public :: kbargq
  public :: kbarqq
  public :: kbargg
  public :: ktildeqq
  public :: ktildeqg
  public :: ktildegq
  public :: ktildegg
  public :: insert_q
  public :: insert_g
  public :: k_q_al, k_g_al
  public :: plus_distr_al
  public :: kbarqg_al
  public :: kbargq_al
  public :: kbarqq_al
  public :: kbargg_al
  public :: ktildeqq_al
  public :: log_plus_distr
  public :: log2_plus_distr
  public :: log2_plus_distr_al
  public :: p_qqg
  public :: p_gqq
  public :: p_ggg
  public :: integral_over_p_qqg
  public :: integral_over_p_gqq
  public :: integral_over_p_ggg
  public :: p_qqg_pol
  public :: pqqm
  public :: top_width_sm_lo
  public :: g_mu_from_alpha
  public :: alpha_from_g_mu
  public :: top_width_sm_qcd_nlo_massless_b
  public :: f0
  public :: f1
  public :: top_width_sm_qcd_nlo_jk
  public :: top_width_sm_qcd_nlo_ce
  public :: ff0
  public :: ff_f0
  public :: ff_lambda
  public :: ff1

  real(default), parameter :: &
       zeta2 = 1.64493406684822643647241516665_default, &
       zeta3 = 1.20205690315959428539973816151_default, &
       zeta4 = 1.08232323371113819151600369654_default, &
       zeta5 = 1.03692775514336992633136548646_default

  real(default), parameter :: &
       eulerc =0.5772156649015328606065120900824024310422_default

  real(default), parameter, public ::  gamma_q = three/two * CF, &
     k_q = (7.0_default/two - pi**2/6.0_default) * CF


  interface
    pure module function beta0 (nf)
      real(default), intent(in) :: nf
      real(default) :: beta0
    end function beta0
    pure module function beta1 (nf)
      real(default), intent(in) :: nf
      real(default) :: beta1
    end function beta1
    pure module function beta2 (nf)
      real(default), intent(in) :: nf
      real(default) :: beta2
    end function beta2
    pure module function coeff_b0 (nf)
      real(default), intent(in) :: nf
      real(default) :: coeff_b0
    end function coeff_b0
    pure module function coeff_b1 (nf)
      real(default), intent(in) :: nf
      real(default) :: coeff_b1
    end function coeff_b1
    pure module function coeff_b2 (nf)
      real(default), intent(in) :: nf
      real(default) :: coeff_b2
    end function coeff_b2
    pure module function coeffqed_b0 (nf, nlep)
      integer, intent(in) :: nf, nlep
      real(default) :: n_lep, coeffqed_b0
    end function coeffqed_b0
    pure module function coeffqed_b1 (nf, nlep)
      integer, intent(in) :: nf, nlep
      real(default) :: n_lep, coeffqed_b1
    end function coeffqed_b1
    pure module function running_as (scale, al_mz, mz, order, nf) result (ascale)
      real(default), intent(in) :: scale
      real(default), intent(in), optional :: al_mz, nf, mz
      integer, intent(in), optional :: order
      real(default) :: ascale
    end function running_as
    pure module function running_as_lam (nf, scale, lambda, order) result (ascale)
      real(default), intent(in) :: nf, scale
      real(default), intent(in), optional :: lambda
      integer, intent(in), optional :: order
      real(default) :: ascale
     end function running_as_lam
    pure module function running_alpha &
         (scale, al_me, me, order, nf, nlep) result (ascale)
      real(default), intent(in) :: scale
      real(default), intent(in), optional :: al_me, me
      integer, intent(in), optional :: order, nf, nlep
      real(default) :: ascale
    end function running_alpha
    pure module function running_alpha_num &
         (scale, al_me, me, order, nf, nlep) result (ascale)
      real(default), intent(in) :: scale
      real(default), intent(in), optional :: al_me, me
      integer, intent(in), optional :: order, nf, nlep
      real(default) :: ascale
    end function running_alpha_num
    module function lambda_qcd (as_q, q, nf, order) result (lambda)
      real(default), intent(in) :: as_q, q
      integer, intent(in) :: order, nf
      real(default) :: lambda
    end function lambda_qcd
    elemental module function gamma_g (nf) result (gg)
      real(default), intent(in) :: nf
      real(default) :: gg
    end function gamma_g
    elemental module function k_g (nf) result (kg)
      real(default), intent(in) :: nf
      real(default) :: kg
    end function k_g
    elemental module function Li2 (x)
      real(default), intent(in) :: x
      real(default) :: Li2
    end function Li2
    elemental module function psic (z) result (psi)
      complex(default), intent(in) :: z
      complex(default) :: psi
    end function psic
    elemental module function psir (x) result (psi)
      real(default), intent(in) :: x
      real(default) :: psi
    end function psir
    elemental module function psim (z, m) result (psi)
      complex(default), intent(in) :: z
      integer, intent(in) :: m
      complex(default) :: psi
    end function psim
    elemental module function psimr (x, m) result (psi)
      real(default), intent(in) :: x
      integer, intent(in) :: m
      real(default) :: psi
    end function psimr
    module function cnielsen (n, m, x) result (nplog)
      integer, intent(in) :: n, m
      real(default), intent(in) :: x
      complex(default) :: nplog
    end function cnielsen
    module function nielsen (n, m, x) result (nplog)
      integer, intent(in) :: n, m
      real(default), intent(in) :: x
      real(default) :: nplog
    end function nielsen
    module function polylog (n, x) result (plog)
      integer, intent(in) :: n
      real(default), intent(in) :: x
      real(default) :: plog
    end function polylog
    module function dilog (x) result (dlog)
      real(default), intent(in) :: x
      real(default) :: dlog
    end function dilog
    module function trilog (x) result (tlog)
      real(default), intent(in) :: x
      real(default) :: tlog
    end function trilog
    elemental module function faux (x) result (y)
      real(default), intent(in) :: x
      complex(default) :: y
    end function faux
    elemental module function fonehalf (x) result (y)
      real(default), intent(in) :: x
      complex(default) :: y
    end function fonehalf
    module function fonehalf_pseudo (x) result (y)
      real(default), intent(in) :: x
      complex(default) :: y
    end function fonehalf_pseudo
    elemental module function fone (x) result  (y)
      real(default), intent(in) :: x
      complex(default) :: y
    end function fone
    elemental module function gaux (x) result (y)
      real(default), intent(in) :: x
      complex(default) :: y
    end function gaux
    elemental module function tri_i1 (a,b) result (y)
      real(default), intent(in) :: a,b
      complex(default) :: y
    end function tri_i1
    elemental module function tri_i2 (a,b) result (y)
      real(default), intent(in) :: a,b
      complex(default) :: y
    end function tri_i2
    elemental module function run_b0 (nf) result (bnull)
      integer, intent(in) :: nf
      real(default) :: bnull
    end function run_b0
    elemental module function run_b1 (nf) result (bone)
      integer, intent(in) :: nf
      real(default) :: bone
    end function run_b1
    elemental module function run_aa (nf) result (aaa)
        integer, intent(in) :: nf
        real(default) :: aaa
    end function run_aa
    pure module subroutine ff_dipole (v_ijk, y_ijk, p_ij, pp_k, p_i, p_j, p_k)
      type(vector4_t), intent(in) :: p_i, p_j, p_k
      type(vector4_t), intent(out) :: p_ij, pp_k
      real(default), intent(out) :: y_ijk
      real(default), intent(out) :: v_ijk
    end subroutine ff_dipole
    pure module subroutine fi_dipole (v_ija, x_ija, p_ij, pp_a, p_i, p_j, p_a)
      type(vector4_t), intent(in) :: p_i,  p_j,  p_a
      type(vector4_t), intent(out) :: p_ij, pp_a
      real(default), intent(out) :: x_ija
      real(default), intent(out) :: v_ija
    end subroutine fi_dipole
    pure module subroutine if_dipole (v_kja, u_j, p_aj, pp_k, p_k, p_j, p_a)
      type(vector4_t), intent(in) :: p_k, p_j, p_a
      type(vector4_t), intent(out) :: p_aj, pp_k
      real(default), intent(out) :: u_j
      real(default), intent(out) :: v_kja
    end subroutine if_dipole
    pure module subroutine ii_dipole (v_jab, v_j, p_in, p_out, flag_1or2)
      type(vector4_t), dimension(:), intent(in) :: p_in
      type(vector4_t), dimension(size(p_in)-1), intent(out) :: p_out
      logical, intent(in) :: flag_1or2
      real(default), intent(out) :: v_j
      real(default), intent(out) :: v_jab
    end subroutine ii_dipole
    elemental module function delta (x,eps) result (z)
       real(default), intent(in) :: x, eps
       real(default) :: z
    end function delta
    elemental module function plus_distr (x,eps) result (plusd)
      real(default), intent(in) :: x, eps
      real(default) :: plusd
    end function plus_distr
    elemental module function pqq (x,eps) result (pqqx)
      real(default), intent(in) :: x, eps
      real(default) :: pqqx
    end function pqq
    elemental module function pgq (x) result (pgqx)
      real(default), intent(in) :: x
      real(default) :: pgqx
    end function pgq
    elemental module function pqg (x) result (pqgx)
      real(default), intent(in) :: x
      real(default) :: pqgx
    end function pqg
    elemental module function pgg (x, nf, eps) result (pggx)
      real(default), intent(in) :: x, nf, eps
      real(default) :: pggx
    end function pgg
    elemental module function pqq_reg (x) result (pqqregx)
       real(default), intent(in) :: x
       real(default) :: pqqregx
    end function pqq_reg
    elemental module function pgg_reg (x) result (pggregx)
       real(default), intent(in) :: x
       real(default) :: pggregx
    end function pgg_reg
    module function kbarqg (x) result (kbarqgx)
      real(default), intent(in) :: x
      real(default) :: kbarqgx
    end function kbarqg
    module function kbargq (x) result (kbargqx)
      real(default), intent(in) :: x
      real(default) :: kbargqx
    end function kbargq
    module function kbarqq (x,eps) result (kbarqqx)
      real(default), intent(in) :: x, eps
      real(default) :: kbarqqx
    end function kbarqq
    module function kbargg (x,eps,nf) result (kbarggx)
      real(default), intent(in) :: x, eps, nf
      real(default) :: kbarggx
    end function kbargg
    module function ktildeqq (x,eps) result (ktildeqqx)
      real(default), intent(in) :: x, eps
      real(default) :: ktildeqqx
    end function ktildeqq
    module function ktildeqg (x,eps) result (ktildeqgx)
      real(default), intent(in) :: x, eps
      real(default) :: ktildeqgx
    end function ktildeqg
    module function ktildegq (x,eps) result (ktildegqx)
      real(default), intent(in) :: x, eps
      real(default) :: ktildegqx
    end function ktildegq
    module function ktildegg (x,eps) result (ktildeggx)
      real(default), intent(in) :: x, eps
      real(default) :: ktildeggx
    end function ktildegg
    pure module function insert_q () result (i_q)
      real(default), dimension(0:2) :: i_q
    end function insert_q
    pure module function insert_g (nf) result (i_g)
      real(default), intent(in) :: nf
      real(default), dimension(0:2) :: i_g
    end function insert_g
    pure module function k_q_al (alpha)
      real(default), intent(in) :: alpha
      real(default) :: k_q_al
    end function k_q_al
    pure module function k_g_al (alpha, nf)
      real(default), intent(in) :: alpha, nf
      real(default) :: k_g_al
    end function k_g_al
    module function plus_distr_al (x,alpha,eps) result (plusd_al)
      real(default), intent(in) :: x,  eps, alpha
      real(default) :: plusd_al
    end function plus_distr_al
    module function kbarqg_al (x,alpha,eps) result (kbarqgx)
      real(default), intent(in) :: x, alpha, eps
      real(default) :: kbarqgx
    end function kbarqg_al
    module function kbargq_al (x,alpha,eps) result (kbargqx)
      real(default), intent(in) :: x, alpha, eps
      real(default) :: kbargqx
    end function kbargq_al
    module function kbarqq_al (x,alpha,eps) result (kbarqqx)
      real(default), intent(in) :: x, alpha, eps
      real(default) :: kbarqqx
    end function kbarqq_al
    module function kbargg_al (x,alpha,eps,nf) result (kbarggx)
      real(default), intent(in) :: x, alpha, eps, nf
      real(default) :: kbarggx
    end function kbargg_al
    module function ktildeqq_al (x,alpha,eps) result (ktildeqqx)
      real(default), intent(in) :: x, eps, alpha
      real(default) :: ktildeqqx
    end function ktildeqq_al
    module function log_plus_distr (x,eps) result (lpd)
       real(default), intent(in) :: x, eps
       real(default) :: lpd, eps2
    end function log_plus_distr
    module function log2_plus_distr (x,eps) result (lpd)
      real(default), intent(in) :: x, eps
      real(default) :: lpd
    end function log2_plus_distr
    module function log2_plus_distr_al (x,alpha,eps) result (lpd_al)
      real(default), intent(in) :: x, eps, alpha
      real(default) :: lpd_al
    end function log2_plus_distr_al
    elemental module function p_qqg (z) result (P)
      real(default), intent(in) :: z
      real(default) :: P
    end function p_qqg
    elemental module function p_gqq (z) result (P)
      real(default), intent(in) :: z
      real(default) :: P
    end function p_gqq
    elemental module function p_ggg (z) result (P)
      real(default), intent(in) :: z
      real(default) :: P
    end function p_ggg
    pure module function integral_over_p_qqg (zmin, zmax) result (integral)
      real(default), intent(in) :: zmin, zmax
      real(default) :: integral
    end function integral_over_p_qqg
    pure module function integral_over_p_gqq (zmin, zmax) result (integral)
      real(default), intent(in) :: zmin, zmax
      real(default) :: integral
    end function integral_over_p_gqq
    pure module function integral_over_p_ggg (zmin, zmax) result (integral)
      real(default), intent(in) :: zmin, zmax
      real(default) :: integral
    end function integral_over_p_ggg
    elemental module function p_qqg_pol (z, l_a, l_b, l_c) result (P)
      real(default), intent(in) :: z
      integer, intent(in) :: l_a, l_b, l_c
      real(default) :: P
    end function p_qqg_pol
    module function pqqm (n, c_f) result (pqq_m)
      integer, intent(in) :: n
      real(default), intent(in) :: c_f
      complex(default) :: pqq_m
    end function pqqm
    elemental module function top_width_sm_lo (alpha, sinthw, vtb, mtop, mw, mb) &
           result (gamma)
      real(default) :: gamma
      real(default), intent(in) :: alpha, sinthw, vtb, mtop, mw, mb
    end function top_width_sm_lo
    elemental module function g_mu_from_alpha (alpha, mw, sinthw) result (g_mu)
      real(default) :: g_mu
      real(default), intent(in) :: alpha, mw, sinthw
    end function g_mu_from_alpha
    elemental module function alpha_from_g_mu (g_mu, mw, sinthw) result (alpha)
      real(default) :: alpha
      real(default), intent(in) :: g_mu, mw, sinthw
    end function alpha_from_g_mu
    elemental module function top_width_sm_qcd_nlo_massless_b &
           (alpha, sinthw, vtb, mtop, mw, alphas) result (gamma)
      real(default) :: gamma
      real(default), intent(in) :: alpha, sinthw, vtb, mtop, mw, alphas
    end function top_width_sm_qcd_nlo_massless_b
    elemental module function f0 (w2) result (f)
      real(default) :: f
      real(default), intent(in) :: w2
    end function f0
    elemental module function f1 (w2) result (f)
      real(default) :: f
      real(default), intent(in) :: w2
    end function f1
    elemental module function top_width_sm_qcd_nlo_jk &
           (alpha, sinthw, vtb, mtop, mw, mb, alphas) result (gamma)
      real(default) :: gamma
      real(default), intent(in) :: alpha, sinthw, vtb, mtop, mw, mb, alphas
    end function top_width_sm_qcd_nlo_jk
    elemental module function top_width_sm_qcd_nlo_ce &
         (alpha, sinthw, vtb, mtop, mw, mb, alpha_s) result (gamma)
      real(default) :: gamma
      real(default), intent(in) :: alpha, sinthw, vtb, mtop, mw, mb, alpha_s
    end function top_width_sm_qcd_nlo_ce
    elemental module function ff0 (eps2, w2) result (f)
      real(default) :: f
      real(default), intent(in) :: eps2, w2
    end function ff0
    elemental module function ff_f0 (eps2, w2) result (f)
      real(default) :: f
      real(default), intent(in) :: eps2, w2
    end function ff_f0
    elemental module function ff_lambda (eps2, w2) result (l)
      real(default) :: l
      real(default), intent(in) :: eps2, w2
    end function ff_lambda
    elemental module function ff1 (eps2, w2) result (f)
      real(default) :: f
      real(default), intent(in) :: eps2, w2
    end function ff1
  end interface

end module sm_physics
