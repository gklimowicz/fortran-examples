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

module sm_physics_uti

  use kinds, only: default
  use numeric_utils
  use format_defs, only: FMT_15
  use constants

  use sm_physics

  implicit none
  private

  public :: sm_physics_1
  public :: sm_physics_2
  public :: sm_physics_3

contains

  subroutine sm_physics_1 (u)
    integer, intent(in) :: u
    real(default) :: z = 0.75_default

    write (u, "(A)")  "* Test output: sm_physics_1"
    write (u, "(A)")  "*   Purpose: check analytic properties"
    write (u, "(A)")

    write (u, "(A)")  "* Splitting functions:"
    write (u, "(A)")

    call assert (u, vanishes (p_qqg_pol (z, +1, -1, +1)), "+-+")
    call assert (u, vanishes (p_qqg_pol (z, +1, -1, -1)), "+--")
    call assert (u, vanishes (p_qqg_pol (z, -1, +1, +1)), "-++")
    call assert (u, vanishes (p_qqg_pol (z, -1, +1, -1)), "-+-")

    !call assert (u, nearly_equal ( &
         !p_qqg_pol (z, +1, +1, -1) + p_qqg_pol (z, +1, +1, +1), &
         !p_qqg (z)), "pol sum")

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sm_physics_1"

  end subroutine sm_physics_1

  subroutine sm_physics_2 (u)
    integer, intent(in) :: u
    real(default) :: mtop, mw, mz, mb, g_mu, sinthw, alpha, vtb, gamma0
    real(default) :: w2, alphas, alphas_mz, gamma1
    write (u, "(A)")  "* Test output: sm_physics_2"
    write (u, "(A)")  "*   Purpose: Check different top width computations"
    write (u, "(A)")

    write (u, "(A)")  "*   Values from [[1207.5018]] (massless b)"
    mtop = 172.0
    mw = 80.399
    mz = 91.1876
    mb = zero
    mb = 0.00001
    g_mu = 1.16637E-5
    sinthw = sqrt(one - mw**2 / mz**2)
    alpha = alpha_from_g_mu (g_mu, mw, sinthw)
    vtb = one
    w2 = mw**2 / mtop**2

    write (u, "(A)")  "*   Check Li2 implementation"
    call assert_equal (u, Li2(w2), 0.2317566263959552_default, &
         "Li2(w2)", rel_smallness=1.0E-6_default)
    call assert_equal (u, Li2(one - w2), 1.038200378935867_default, &
         "Li2(one - w2)", rel_smallness=1.0E-6_default)

    write (u, "(A)")  "*   Check LO Width"
    gamma0 = top_width_sm_lo (alpha, sinthw, vtb, mtop, mw, mb)
    call assert_equal (u, gamma0, 1.4655_default, &
         "top_width_sm_lo", rel_smallness=1.0E-5_default)
    alphas = zero
    gamma0 = top_width_sm_qcd_nlo_massless_b &
         (alpha, sinthw, vtb, mtop, mw, alphas)
    call assert_equal (u, gamma0, 1.4655_default, &
         "top_width_sm_qcd_nlo_massless_b", rel_smallness=1.0E-5_default)
    gamma0 = top_width_sm_qcd_nlo_jk &
         (alpha, sinthw, vtb, mtop, mw, mb, alphas)
    call assert_equal (u, gamma0, 1.4655_default, &
         "top_width_sm_qcd_nlo", rel_smallness=1.0E-5_default)

    write (u, "(A)")  "*   Check NLO Width"
    alphas_mz = 0.1202      ! MSTW2008 NLO fit
    alphas = running_as (mtop, alphas_mz, mz, 1, 5.0_default)
    gamma1 = top_width_sm_qcd_nlo_massless_b &
         (alpha, sinthw, vtb, mtop, mw, alphas)
    call assert_equal (u, gamma1, 1.3376_default, rel_smallness=1.0E-4_default)
    gamma1 = top_width_sm_qcd_nlo_jk &
         (alpha, sinthw, vtb, mtop, mw, mb, alphas)
    ! It would be nice to get one more significant digit but the
    ! expression is numerically rather unstable for mb -> 0
    call assert_equal (u, gamma1, 1.3376_default, rel_smallness=1.0E-3_default)

    write (u, "(A)")  "*   Values from threshold validation (massive b)"
    alpha = one / 125.924
    ! ee = 0.315901
    ! cw = 0.881903
    ! v = 240.024
    mtop = 172.0 ! This is the value for M1S !!!
    mb = 4.2
    sinthw = 0.47143
    mz = 91.188
    mw = 80.419
    call assert_equal (u, sqrt(one - mw**2 / mz**2), sinthw, &
         "sinthw", rel_smallness=1.0E-6_default)

    write (u, "(A)")  "*   Check LO Width"
    gamma0 = top_width_sm_lo (alpha, sinthw, vtb, mtop, mw, mb)
    call assert_equal (u, gamma0, 1.5386446_default, &
         "gamma0", rel_smallness=1.0E-7_default)
    alphas = zero
    gamma0 = top_width_sm_qcd_nlo_jk &
         (alpha, sinthw, vtb, mtop, mw, mb, alphas)
    call assert_equal (u, gamma0, 1.5386446_default, &
         "gamma0", rel_smallness=1.0E-7_default)

    write (u, "(A)")  "*   Check NLO Width"
    alphas_mz = 0.118 !(Z pole, NLL running to mu_h)
    alphas = running_as (mtop, alphas_mz, mz, 1, 5.0_default)
    write (u, "(A," // FMT_15 // ")")  "*   alphas = ", alphas
    gamma1 = top_width_sm_qcd_nlo_jk &
         (alpha, sinthw, vtb, mtop, mw, mb, alphas)
    write (u, "(A," // FMT_15 // ")")  "*   Gamma1 = ", gamma1

    mb = zero
    gamma1 = top_width_sm_qcd_nlo_massless_b &
         (alpha, sinthw, vtb, mtop, mw, alphas)
    alphas = running_as (mtop, alphas_mz, mz, 1, 5.0_default)
    write (u, "(A," // FMT_15 // ")")  "*   Gamma1(mb=0) = ", gamma1

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sm_physics_2"
  end subroutine sm_physics_2

  subroutine sm_physics_3 (u)
    integer, intent(in) :: u
    complex(default) :: z1 = (0.75_default, 1.25_default)
    complex(default) :: z2 = (1.33_default, 11.25_default)
    complex(default) :: psiz
    real(default) :: x1 = 0.045847700_default
    real(default) :: psir

    write (u, "(A)")  "* Test output: sm_physics_3"
    write (u, "(A)")  "*   Purpose: check special functions"
    write (u, "(A)")

    write (u, "(A)")  "* Complex digamma function:"
    write (u, "(A)")

    psiz = psic (z1)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z1      = ", &
         real(z1), aimag(z1)
    write (u, "(1x,A,'(',F7.5,',',F7.5,')')")  " psi(z1) = ", &
         real(psiz), aimag(psiz)
    psiz = psic (z2)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z2      = ", &
         real(z2), aimag(z2)
    write (u, "(1x,A,'(',F7.5,',',F7.5,')')")  " psi(z2) = ", &
         real(psiz), aimag(psiz)

    write (u, "(A)")
    write (u, "(A)")  "* Complex polygamma function:"
    write (u, "(A)")

    psiz = psim (z1,1)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z1        = ", &
         real(z1), aimag(z1)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " psi(z1,1) = ", &
         real(psiz), aimag(psiz)
    psiz = psim (z2,1)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z2        = ", &
         real(z2), aimag(z2)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " psi(z2,1) = ", &
         real(psiz), aimag(psiz)

    write (u, "(A)")

    psiz = psim (z1,2)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z1        = ", &
         real(z1), aimag(z1)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " psi(z1,2) = ", &
         real(psiz), aimag(psiz)
    psiz = psim (z2,2)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z2        = ", &
         real(z2), aimag(z2)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " psi(z2,2) = ", &
         real(psiz), aimag(psiz)

    write (u, "(A)")

    psiz = psim (z1,3)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z1        = ", &
         real(z1), aimag(z1)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " psi(z1,3) = ", &
         real(psiz), aimag(psiz)
    psiz = psim (z2,3)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z2        = ", &
         real(z2), aimag(z2)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " psi(z2,3) = ", &
         real(psiz), aimag(psiz)

    write (u, "(A)")

    psiz = psim (z1,4)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z1        = ", &
         real(z1), aimag(z1)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " psi(z1,4) = ", &
         real(psiz), aimag(psiz)
    psiz = psim (z2,4)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z2        = ", &
         real(z2), aimag(z2)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " psi(z2,4) = ", &
         real(psiz), aimag(psiz)

    write (u, "(A)")

    psiz = psim (z1,5)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z1        = ", &
         real(z1), aimag(z1)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " psi(z1,5) = ", &
         real(psiz), aimag(psiz)
    psiz = psim (z2,5)
    write (u, "(1x,A,'(',F5.2,',',F5.2,')')")  " z2        = ", &
         real(z2), aimag(z2)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " psi(z2,5) = ", &
         real(psiz), aimag(psiz)

    write (u, "(A)")
    write (u, "(A)")  "* Real polygamma function:"
    write (u, "(A)")

    psir = psimr (x1,1)
    write (u, "(1x,A,'(',F8.5,')')")  " x1        = ", x1
    write (u, "(1x,A,'(',F8.4,')')")  " psir      = ", psir

    write (u, "(A)")
    write (u, "(A)")  "* Generalized Nielsen polylogarithm:"
    write (u, "(A)")

    write (u, "(1x,A,F8.5)")  " S(1,1,0) = ", &
         nielsen(1,1,0._default)
    write (u, "(1x,A,F8.5)")  " S(1,1,-1) = ", &
         nielsen(1,1,-1._default)
    write (u, "(1x,A,F8.5)")  " S(1,2,-1) = ", &
         nielsen(1,2,-1._default)
    write (u, "(1x,A,F8.5)")  " S(2,1,-1) = ", &
         nielsen(2,1,-1._default)
    write (u, "(1x,A,F8.5)")  " S(1,3,-1) = ", &
         nielsen(1,3,-1._default)
    write (u, "(1x,A,F8.5)")  " S(2,2,-1) = ", &
         nielsen(2,2,-1._default)
    write (u, "(1x,A,F8.5)")  " S(3,1,-1) = ", &
         nielsen(3,1,-1._default)
    write (u, "(1x,A,F8.5)")  " S(1,4,-1) = ", &
         nielsen(1,4,-1._default)
    write (u, "(1x,A,F8.5)")  " S(2,3,-1) = ", &
         nielsen(2,3,-1._default)
    write (u, "(1x,A,F8.5)")  " S(3,2,-1) = ", &
         nielsen(3,2,-1._default)
    write (u, "(1x,A,F8.5)")  " S(4,1,-1) = ", &
         nielsen(4,1,-1._default)
    write (u, "(1x,A,F8.5)")  " S(1,1,0.2) = ", &
         nielsen(1,1,0.2_default)
    write (u, "(1x,A,F8.5)")  " S(1,2,0.2) = ", &
         nielsen(1,2,0.2_default)
    write (u, "(1x,A,F8.5)")  " S(2,1,0.2) = ", &
         nielsen(2,1,0.2_default)
    write (u, "(1x,A,F8.5)")  " S(1,3,0.2) = ", &
         nielsen(1,3,0.2_default)
    write (u, "(1x,A,F8.5)")  " S(2,2,0.2) = ", &
         nielsen(2,2,0.2_default)
    write (u, "(1x,A,F8.5)")  " S(3,1,0.2) = ", &
         nielsen(3,1,0.2_default)
    write (u, "(1x,A,F8.5)")  " S(1,4,0.2) = ", &
         nielsen(1,4,0.2_default)
    write (u, "(1x,A,F8.5)")  " S(2,3,0.2) = ", &
         nielsen(2,3,0.2_default)
    write (u, "(1x,A,F8.5)")  " S(3,2,0.2) = ", &
         nielsen(3,2,0.2_default)
    write (u, "(1x,A,F8.5)")  " S(4,1,0.2) = ", &
         nielsen(4,1,0.2_default)
    write (u, "(1x,A,F8.5)")  " S(1,1,1) = ", &
         nielsen(1,1,1._default)
    write (u, "(1x,A,F8.5)")  " S(1,2,1) = ", &
         nielsen(1,2,1._default)
    write (u, "(1x,A,F8.5)")  " S(2,1,1) = ", &
         nielsen(2,1,1._default)
    write (u, "(1x,A,F8.5)")  " S(1,3,1) = ", &
         nielsen(1,3,1._default)
    write (u, "(1x,A,F8.5)")  " S(2,2,1) = ", &
         nielsen(2,2,1._default)
    write (u, "(1x,A,F8.5)")  " S(3,1,1) = ", &
         nielsen(3,1,1._default)
    write (u, "(1x,A,F8.5)")  " S(1,4,1) = ", &
         nielsen(1,4,1._default)
    write (u, "(1x,A,F8.5)")  " S(2,3,1) = ", &
         nielsen(2,3,1._default)
    write (u, "(1x,A,F8.5)")  " S(3,2,1) = ", &
         nielsen(3,2,1._default)
    write (u, "(1x,A,F8.5)")  " S(4,1,1) = ", &
         nielsen(4,1,1._default)
    write (u, "(1x,A,F8.5)")  " S(1,1,0.75) = ", &
         nielsen(1,1,0.75_default)
    write (u, "(1x,A,F8.5)")  " S(1,3,0.75) = ", &
         nielsen(1,3,0.75_default)
    write (u, "(1x,A,F8.5)")  " S(1,4,0.75) = ", &
         nielsen(1,4,0.75_default)
    write (u, "(1x,A,F8.5)")  " S(2,2,0.75) = ", &
         nielsen(2,2,0.75_default)
    write (u, "(1x,A,'(',F8.5,',',F8.5,')')")  " S(1,1,2) = ", &
         real(cnielsen(1,1,3._default)), &
         aimag(cnielsen(1,1,3._default))

    write (u, "(A)")
    write (u, "(A)")  "* Dilog, trilog, polylog:"
    write (u, "(A)")

    write (u, "(1x,A,F8.5)")  " Li2(0.66)    = ", &
         dilog(0.66_default)
    write (u, "(1x,A,F8.5)")  " Li3(0.66)    = ", &
         trilog(0.66_default)
    write (u, "(1x,A,F8.5)")  " Poly(4,0.66) = ", &
         polylog(4,0.66_default)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sm_physics_3"

  end subroutine sm_physics_3


end module sm_physics_uti
