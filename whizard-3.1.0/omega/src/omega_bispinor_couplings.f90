!  omegalib.nw --
!
!  Copyright (C) 1999-2022 by
!      Wolfgang Kilian <kilian@physik.uni-siegen.de>
!      Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!      Juergen Reuter <juergen.reuter@desy.de>
!      with contributions from                                                                                                                                    
!      Fabian Bach <fabian.bach@t-online.de>                                                                                                                 
!      Bijan Chokoufe Nejad <bijan.chokoufe@desy.de>                                                                                                              
!      Christian Speckner <cnspeckn@googlemail.com>     
!
!  WHIZARD is free software; you can redistribute it and/or modify it
!  under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2, or (at your option)
!  any later version.
!
!  WHIZARD is distributed in the hope that it will be useful, but
!  WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module omega_bispinor_couplings
  use kinds
  use constants
  use omega_bispinors
  use omega_vectorspinors
  use omega_vectors
  use omega_couplings
  implicit none
  private
  public :: u, v, ghost
  public :: brs_u, brs_v
  public :: va_ff, v_ff, a_ff, vl_ff, vr_ff, vlr_ff, va2_ff, tva_ff, tvam_ff, &
            tlr_ff, tlrm_ff
  public :: f_vaf, f_vf, f_af, f_vlf, f_vrf, f_vlrf, f_va2f, &
            f_tvaf, f_tlrf, f_tvamf, f_tlrmf
  public :: sp_ff, s_ff, p_ff, sl_ff, sr_ff, slr_ff
  public :: f_spf, f_sf, f_pf, f_slf, f_srf, f_slrf
  private :: vv_ff, f_vvf
  public :: vmom_ff, mom_ff, mom5_ff, moml_ff, momr_ff, lmom_ff, rmom_ff
  public :: f_vmomf, f_momf, f_mom5f, f_momlf, f_momrf, f_lmomf, f_rmomf
  public :: v2_ff, sv1_ff, sv2_ff, pv1_ff, pv2_ff, svl1_ff, svl2_ff, &
       svr1_ff, svr2_ff, svlr1_ff, svlr2_ff
  public :: f_v2f, f_svf, f_pvf, f_svlf, f_svrf, f_svlrf
  public :: pot_grf, pot_fgr, s_grf, s_fgr, p_grf, p_fgr, &
       sl_grf, sl_fgr, sr_grf, sr_fgr, slr_grf, slr_fgr
  private :: fgvgr, fgvg5gr, fggvvgr, grkgf, grkggf, grkkggf, &
       fgkgr, fg5gkgr, grvgf, grg5vgf, grkgggf, fggkggr
  public :: f_potgr, f_sgr, f_pgr, f_vgr, f_vlrgr, f_slgr, f_srgr, f_slrgr
  public :: gr_potf, gr_sf, gr_pf, gr_vf, gr_vlrf, gr_slf, gr_srf, gr_slrf
  public :: v_grf, v_fgr
  public :: vlr_grf, vlr_fgr
  public :: f_s2gr, f_svgr, f_slvgr, f_srvgr, f_slrvgr, f_pvgr, f_v2gr, f_v2lrgr
  public :: gr_s2f, gr_svf, gr_pvf, gr_slvf, gr_srvf, gr_slrvf, gr_v2f, gr_v2lrf
  public :: s2_grf, s2_fgr, sv1_grf, sv2_grf, sv1_fgr, sv2_fgr, &
            slv1_grf, slv2_grf, slv1_fgr, slv2_fgr, &
            srv1_grf, srv2_grf, srv1_fgr, srv2_fgr, &
            slrv1_grf, slrv2_grf, slrv1_fgr, slrv2_fgr, &
            pv1_grf, pv2_grf, pv1_fgr, pv2_fgr, v2_grf, v2_fgr, &
            v2lr_grf, v2lr_fgr
  public :: pr_psi, pr_grav
  public :: pj_psi, pg_psi
  integer, parameter, public :: omega_bispinor_cpls_2010_01_A = 0
contains
  pure function u (mass, p, s) result (psi)
    type(bispinor) :: psi
    real(kind=default), intent(in) :: mass
    type(momentum), intent(in) :: p
    integer, intent(in) :: s
    complex(kind=default), dimension(2) :: chip, chim
    real(kind=default) :: pabs, norm, delta, m
    m = abs(mass)
    pabs = sqrt (dot_product (p%x, p%x))
    if (m < epsilon (m) * pabs) then
        delta = 0
    else
        delta = sqrt (max (p%t - pabs, 0._default))
    end if
    if (pabs + p%x(3) <= 1000 * epsilon (pabs) * pabs) then
       chip = (/ cmplx ( 0.0, 0.0, kind=default), &
                 cmplx ( 1.0, 0.0, kind=default) /)
       chim = (/ cmplx (-1.0, 0.0, kind=default), &
                 cmplx ( 0.0, 0.0, kind=default) /)
    else
       norm = 1 / sqrt (2*pabs*(pabs + p%x(3)))
       chip = norm * (/ cmplx (pabs + p%x(3), kind=default), &
                        cmplx (p%x(1), p%x(2), kind=default) /)
       chim = norm * (/ cmplx (-p%x(1), p%x(2), kind=default), &
                        cmplx (pabs + p%x(3), kind=default) /)
    end if
    if (s > 0) then
       psi%a(1:2) = delta * chip
       psi%a(3:4) = sqrt (p%t + pabs) * chip
    else
       psi%a(1:2) = sqrt (p%t + pabs) * chim
       psi%a(3:4) = delta * chim
    end if
    pabs = m ! make the compiler happy and use m
    if (mass < 0) then
       psi%a(1:2) = - imago * psi%a(1:2)
       psi%a(3:4) = + imago * psi%a(3:4)
    end if
  end function u
  pure function v (mass, p, s) result (psi)
    type(bispinor) :: psi
    real(kind=default), intent(in) :: mass
    type(momentum), intent(in) :: p
    integer, intent(in) :: s
    complex(kind=default), dimension(2) :: chip, chim
    real(kind=default) :: pabs, norm, delta, m
    pabs = sqrt (dot_product (p%x, p%x))
    m = abs(mass)
    if (m < epsilon (m) * pabs) then
        delta = 0
    else
        delta = sqrt (max (p%t - pabs, 0._default))
    end if
    if (pabs + p%x(3) <= 1000 * epsilon (pabs) * pabs) then
       chip = (/ cmplx ( 0.0, 0.0, kind=default), &
                 cmplx ( 1.0, 0.0, kind=default) /)
       chim = (/ cmplx (-1.0, 0.0, kind=default), &
                 cmplx ( 0.0, 0.0, kind=default) /)
    else
       norm = 1 / sqrt (2*pabs*(pabs + p%x(3)))
       chip = norm * (/ cmplx (pabs + p%x(3), kind=default), &
                        cmplx (p%x(1), p%x(2), kind=default) /)
       chim = norm * (/ cmplx (-p%x(1), p%x(2), kind=default), &
                        cmplx (pabs + p%x(3), kind=default) /)
    end if
    if (s > 0) then
       psi%a(1:2) = - sqrt (p%t + pabs) * chim
       psi%a(3:4) = delta * chim
    else
       psi%a(1:2) = delta * chip
       psi%a(3:4) = - sqrt (p%t + pabs) * chip
    end if
    pabs = m ! make the compiler happy and use m
    if (mass < 0) then
       psi%a(1:2) = - imago * psi%a(1:2)
       psi%a(3:4) = + imago * psi%a(3:4)
    end if
  end function v
  pure function ghost (m, p, s) result (psi)
      type(bispinor) :: psi
      real(kind=default), intent(in) :: m
      type(momentum), intent(in) :: p
      integer, intent(in) :: s
      psi%a(:) = 0
      select case (s)
      case (1)
         psi%a(1)   = 1
         psi%a(2:4) = 0
      case (2)
         psi%a(1)   = 0
         psi%a(2)   = 1
         psi%a(3:4) = 0
      case (3)
         psi%a(1:2) = 0
         psi%a(3)   = 1
         psi%a(4)   = 0
      case (4)
         psi%a(1:3) = 0
         psi%a(4)   = 1
      case (5)
         psi%a(1) =    1.4
         psi%a(2) = -  2.3
         psi%a(3) = - 71.5
         psi%a(4) =    0.1
      end select
  end function ghost
  pure function brs_u (m, p, s) result (dpsi)
      type(bispinor) :: dpsi, psi
      real(kind=default), intent(in) :: m
      type(momentum), intent(in) :: p
      integer, intent(in) :: s
      type (vector)::vp
      complex(kind=default), parameter :: one = (1, 0)
      vp=p
      psi=u(m,p,s)
      dpsi=cmplx(0.0,-1.0)*(f_vf(one,vp,psi)-m*psi)
  end function brs_u
  pure function brs_v (m, p, s) result (dpsi)
      type(bispinor) :: dpsi, psi
      real(kind=default), intent(in) :: m
      type(momentum), intent(in) :: p
      integer, intent(in) ::   s
      type (vector)::vp
      complex(kind=default), parameter :: one = (1, 0)
      vp=p
      psi=v(m,p,s)
      dpsi=cmplx(0.0,1.0)*(f_vf(one,vp,psi)+m*psi)
  end function brs_v
  pure function va_ff (gv, ga, psil, psir) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gv, ga
    type(bispinor), intent(in) :: psil, psir
    complex(kind=default) :: gl, gr
    complex(kind=default) :: g13, g14, g23, g24, g31, g32, g41, g42
    gl = gv + ga
    gr = gv - ga
    g13 = psil%a(1)*psir%a(3)
    g14 = psil%a(1)*psir%a(4)
    g23 = psil%a(2)*psir%a(3)
    g24 = psil%a(2)*psir%a(4)
    g31 = psil%a(3)*psir%a(1)
    g32 = psil%a(3)*psir%a(2)
    g41 = psil%a(4)*psir%a(1)
    g42 = psil%a(4)*psir%a(2)
    j%t    =  gr * (   g14 - g23) + gl * ( - g32 + g41)
    j%x(1) =  gr * (   g13 - g24) + gl * (   g31 - g42)
    j%x(2) = (gr * (   g13 + g24) + gl * (   g31 + g42)) * (0, 1)
    j%x(3) =  gr * ( - g14 - g23) + gl * ( - g32 - g41)
  end function va_ff
  pure function va2_ff (gva, psil, psir) result (j)
    type(vector) :: j
    complex(kind=default), intent(in), dimension(2) :: gva
    type(bispinor), intent(in) :: psil, psir
    complex(kind=default) :: gl, gr
    complex(kind=default) :: g13, g14, g23, g24, g31, g32, g41, g42
    gl = gva(1) + gva(2)
    gr = gva(1) - gva(2)
    g13 = psil%a(1)*psir%a(3)
    g14 = psil%a(1)*psir%a(4)
    g23 = psil%a(2)*psir%a(3)
    g24 = psil%a(2)*psir%a(4)
    g31 = psil%a(3)*psir%a(1)
    g32 = psil%a(3)*psir%a(2)
    g41 = psil%a(4)*psir%a(1)
    g42 = psil%a(4)*psir%a(2)
    j%t    =  gr * (   g14 - g23) + gl * ( - g32 + g41)
    j%x(1) =  gr * (   g13 - g24) + gl * (   g31 - g42)
    j%x(2) = (gr * (   g13 + g24) + gl * (   g31 + g42)) * (0, 1)
    j%x(3) =  gr * ( - g14 - g23) + gl * ( - g32 - g41)
  end function va2_ff
  pure function v_ff (gv, psil, psir) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gv
    type(bispinor), intent(in) :: psil, psir
    complex(kind=default) :: g13, g14, g23, g24, g31, g32, g41, g42
    g13 = psil%a(1)*psir%a(3)
    g14 = psil%a(1)*psir%a(4)
    g23 = psil%a(2)*psir%a(3)
    g24 = psil%a(2)*psir%a(4)
    g31 = psil%a(3)*psir%a(1)
    g32 = psil%a(3)*psir%a(2)
    g41 = psil%a(4)*psir%a(1)
    g42 = psil%a(4)*psir%a(2)
    j%t    =   gv * (   g14 - g23 - g32 + g41)
    j%x(1) =   gv * (   g13 - g24 + g31 - g42)
    j%x(2) =   gv * (   g13 + g24 + g31 + g42) * (0, 1)
    j%x(3) =   gv * ( - g14 - g23 - g32 - g41)
  end function v_ff
  pure function a_ff (ga, psil, psir) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: ga
    type(bispinor), intent(in) :: psil, psir
    complex(kind=default) :: g13, g14, g23, g24, g31, g32, g41, g42
    g13 = psil%a(1)*psir%a(3)
    g14 = psil%a(1)*psir%a(4)
    g23 = psil%a(2)*psir%a(3)
    g24 = psil%a(2)*psir%a(4)
    g31 = psil%a(3)*psir%a(1)
    g32 = psil%a(3)*psir%a(2)
    g41 = psil%a(4)*psir%a(1)
    g42 = psil%a(4)*psir%a(2)
    j%t    =  -ga * (   g14 - g23 + g32 - g41)
    j%x(1) =  -ga * (   g13 - g24 - g31 + g42)
    j%x(2) =  -ga * (   g13 + g24 - g31 - g42) * (0, 1)
    j%x(3) =  -ga * ( - g14 - g23 + g32 + g41)
  end function a_ff
  pure function vl_ff (gl, psil, psir) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl
    type(bispinor), intent(in) :: psil, psir
    complex(kind=default) :: gl2
    complex(kind=default) :: g31, g32, g41, g42
    gl2 = 2 * gl
    g31 = psil%a(3)*psir%a(1)
    g32 = psil%a(3)*psir%a(2)
    g41 = psil%a(4)*psir%a(1)
    g42 = psil%a(4)*psir%a(2)
    j%t    =   gl2 * ( - g32 + g41)
    j%x(1) =   gl2 * (   g31 - g42)
    j%x(2) =   gl2 * (   g31 + g42) * (0, 1)
    j%x(3) =   gl2 * ( - g32 - g41)
  end function vl_ff
  pure function vr_ff (gr, psil, psir) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gr
    type(bispinor), intent(in) :: psil, psir
    complex(kind=default) :: gr2
    complex(kind=default) :: g13, g14, g23, g24
    gr2 = 2 * gr
    g13 = psil%a(1)*psir%a(3)
    g14 = psil%a(1)*psir%a(4)
    g23 = psil%a(2)*psir%a(3)
    g24 = psil%a(2)*psir%a(4)
    j%t    = gr2 * (   g14 - g23)
    j%x(1) = gr2 * (   g13 - g24)
    j%x(2) = gr2 * (   g13 + g24) * (0, 1)
    j%x(3) = gr2 * ( - g14 - g23)
  end function vr_ff
  pure function vlr_ff (gl, gr, psibar, psi) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(bispinor), intent(in) :: psibar
    type(bispinor), intent(in) :: psi
    j = va_ff (gl+gr, gl-gr, psibar, psi)
  end function vlr_ff
  pure function tva_ff (gv, ga, psibar, psi) result (t)
    type(tensor2odd) :: t
    complex(kind=default), intent(in) :: gv, ga
    type(bispinor), intent(in) :: psibar
    type(bispinor), intent(in) :: psi
    complex(kind=default) :: gl, gr
    complex(kind=default) :: g11, g22, g33, g44, g1p2, g3p4
    gr     = gv + ga
    gl     = gv - ga
    g11    = psibar%a(1)*psi%a(1)
    g22    = psibar%a(2)*psi%a(2)
    g1p2   = psibar%a(1)*psi%a(2) + psibar%a(2)*psi%a(1)
    g3p4   = psibar%a(3)*psi%a(4) + psibar%a(4)*psi%a(3)
    g33    = psibar%a(3)*psi%a(3)
    g44    = psibar%a(4)*psi%a(4)
    t%e(1) = (gl * ( - g11 + g22) + gr * ( - g33 + g44)) * (0, 1)
    t%e(2) =  gl * (   g11 + g22) + gr * (   g33 + g44)  
    t%e(3) = (gl * (   g1p2     ) + gr * (   g3p4     )) * (0, 1)  
    t%b(1) =  gl * (   g11 - g22) + gr * ( - g33 + g44)  
    t%b(2) = (gl * (   g11 + g22) + gr * ( - g33 - g44)) * (0, 1)
    t%b(3) =  gl * ( - g1p2     ) + gr * (   g3p4     )
  end function tva_ff
  pure function tlr_ff (gl, gr, psibar, psi) result (t)
    type(tensor2odd) :: t
    complex(kind=default), intent(in) :: gl, gr
    type(bispinor), intent(in) :: psibar
    type(bispinor), intent(in) :: psi
    t = tva_ff (gr+gl, gr-gl, psibar, psi)
  end function tlr_ff
  pure function tvam_ff (gv, ga, psibar, psi, p) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gv, ga
    type(bispinor), intent(in) :: psibar
    type(bispinor), intent(in) :: psi
    type(momentum), intent(in) :: p
    j = (tva_ff(gv, ga, psibar, psi) * p) * (0,1)
  end function tvam_ff
  pure function tlrm_ff (gl, gr, psibar, psi, p) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(bispinor), intent(in) :: psibar
    type(bispinor), intent(in) :: psi
    type(momentum), intent(in) :: p
    j = tvam_ff (gr+gl, gr-gl, psibar, psi, p)
  end function tlrm_ff
  pure function f_vaf (gv, ga, v, psi) result (vpsi)
    type(bispinor) :: vpsi
    complex(kind=default), intent(in) :: gv, ga
    type(vector), intent(in) :: v
    type(bispinor), intent(in) :: psi
    complex(kind=default) :: gl, gr
    complex(kind=default) :: vp, vm, v12, v12s
    gl = gv + ga
    gr = gv - ga
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = gr * (   vm  * psi%a(3) - v12s * psi%a(4))
    vpsi%a(2) = gr * ( - v12 * psi%a(3) + vp   * psi%a(4))
    vpsi%a(3) = gl * (   vp  * psi%a(1) + v12s * psi%a(2))
    vpsi%a(4) = gl * (   v12 * psi%a(1) + vm   * psi%a(2))
  end function f_vaf
  pure function f_va2f (gva, v, psi) result (vpsi)
    type(bispinor) :: vpsi
    complex(kind=default), intent(in), dimension(2) :: gva
    type(vector), intent(in) :: v
    type(bispinor), intent(in) :: psi
    complex(kind=default) :: gl, gr
    complex(kind=default) :: vp, vm, v12, v12s
    gl = gva(1) + gva(2)
    gr = gva(1) - gva(2)
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = gr * (   vm  * psi%a(3) - v12s * psi%a(4))
    vpsi%a(2) = gr * ( - v12 * psi%a(3) + vp   * psi%a(4))
    vpsi%a(3) = gl * (   vp  * psi%a(1) + v12s * psi%a(2))
    vpsi%a(4) = gl * (   v12 * psi%a(1) + vm   * psi%a(2))
  end function f_va2f
  pure function f_vf (gv, v, psi) result (vpsi)
    type(bispinor) :: vpsi
    complex(kind=default), intent(in) :: gv
    type(vector), intent(in) :: v
    type(bispinor), intent(in) :: psi
    complex(kind=default) :: vp, vm, v12, v12s
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = gv * (   vm  * psi%a(3) - v12s * psi%a(4))
    vpsi%a(2) = gv * ( - v12 * psi%a(3) + vp   * psi%a(4))
    vpsi%a(3) = gv * (   vp  * psi%a(1) + v12s * psi%a(2))
    vpsi%a(4) = gv * (   v12 * psi%a(1) + vm   * psi%a(2))
  end function f_vf
  pure function f_af (ga, v, psi) result (vpsi)
    type(bispinor) :: vpsi
    complex(kind=default), intent(in) :: ga
    type(vector), intent(in) :: v
    type(bispinor), intent(in) :: psi
    complex(kind=default) :: vp, vm, v12, v12s
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = ga * ( - vm  * psi%a(3) + v12s * psi%a(4))
    vpsi%a(2) = ga * (   v12 * psi%a(3) - vp   * psi%a(4))
    vpsi%a(3) = ga * (   vp  * psi%a(1) + v12s * psi%a(2))
    vpsi%a(4) = ga * (   v12 * psi%a(1) + vm   * psi%a(2))
  end function f_af
  pure function f_vlf (gl, v, psi) result (vpsi)
    type(bispinor) :: vpsi
    complex(kind=default), intent(in) :: gl
    type(vector), intent(in) :: v
    type(bispinor), intent(in) :: psi
    complex(kind=default) :: gl2
    complex(kind=default) :: vp, vm, v12, v12s
    gl2 = 2 * gl
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = 0
    vpsi%a(2) = 0
    vpsi%a(3) = gl2 * (   vp  * psi%a(1) + v12s * psi%a(2))
    vpsi%a(4) = gl2 * (   v12 * psi%a(1) + vm   * psi%a(2))
  end function f_vlf
  pure function f_vrf (gr, v, psi) result (vpsi)
    type(bispinor) :: vpsi
    complex(kind=default), intent(in) :: gr
    type(vector), intent(in) :: v
    type(bispinor), intent(in) :: psi
    complex(kind=default) :: gr2
    complex(kind=default) :: vp, vm, v12, v12s
    gr2 = 2 * gr
    vp = v%t + v%x(3)
    vm = v%t - v%x(3)
    v12  =  v%x(1) + (0,1)*v%x(2)
    v12s =  v%x(1) - (0,1)*v%x(2)
    vpsi%a(1) = gr2 * (   vm  * psi%a(3) - v12s * psi%a(4))
    vpsi%a(2) = gr2 * ( - v12 * psi%a(3) + vp   * psi%a(4))
    vpsi%a(3) = 0
    vpsi%a(4) = 0
  end function f_vrf
  pure function f_vlrf (gl, gr, v, psi) result (vpsi)
    type(bispinor) :: vpsi
    complex(kind=default), intent(in) :: gl, gr
    type(vector), intent(in) :: v
    type(bispinor), intent(in) :: psi
    vpsi = f_vaf (gl+gr, gl-gr, v, psi)
  end function f_vlrf
  pure function f_tvaf (gv, ga, t, psi) result (tpsi)
    type(bispinor) :: tpsi
    complex(kind=default), intent(in) :: gv, ga
    type(tensor2odd), intent(in) :: t
    type(bispinor), intent(in) :: psi
    complex(kind=default) :: gl, gr
    complex(kind=default) :: e21, e21s, b12, b12s, be3, be3s
    gr   = gv + ga
    gl   = gv - ga
    e21  = t%e(2) + t%e(1)*(0,1)
    e21s = t%e(2) - t%e(1)*(0,1)
    b12  = t%b(1) + t%b(2)*(0,1)
    b12s = t%b(1) - t%b(2)*(0,1)
    be3  = t%b(3) + t%e(3)*(0,1)
    be3s = t%b(3) - t%e(3)*(0,1)
    tpsi%a(1) =   2*gl * (   psi%a(1) * be3  + psi%a(2) * ( e21 +b12s))
    tpsi%a(2) =   2*gl * ( - psi%a(2) * be3  + psi%a(1) * (-e21s+b12 ))
    tpsi%a(3) =   2*gr * (   psi%a(3) * be3s + psi%a(4) * (-e21 +b12s))
    tpsi%a(4) =   2*gr * ( - psi%a(4) * be3s + psi%a(3) * ( e21s+b12 ))
  end function f_tvaf
  pure function f_tlrf (gl, gr, t, psi) result (tpsi)
    type(bispinor) :: tpsi
    complex(kind=default), intent(in) :: gl, gr
    type(tensor2odd), intent(in) :: t
    type(bispinor), intent(in) :: psi
    tpsi = f_tvaf (gr+gl, gr-gl, t, psi)
  end function f_tlrf
  pure function f_tvamf (gv, ga, v, psi, k) result (vpsi)
    type(bispinor) :: vpsi
    complex(kind=default), intent(in) :: gv, ga
    type(vector), intent(in) :: v
    type(bispinor), intent(in) :: psi
    type(momentum), intent(in) :: k
    type(tensor2odd) :: t
    t = (v.wedge.k) * (0, 0.5)
    vpsi = f_tvaf(gv, ga, t, psi)
  end function f_tvamf
  pure function f_tlrmf (gl, gr, v, psi, k) result (vpsi)
    type(bispinor) :: vpsi
    complex(kind=default), intent(in) :: gl, gr
    type(vector), intent(in) :: v
    type(bispinor), intent(in) :: psi
    type(momentum), intent(in) :: k
    vpsi = f_tvamf (gr+gl, gr-gl, v, psi, k)
  end function f_tlrmf
  pure function sp_ff (gs, gp, psil, psir) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gs, gp
    type(bispinor), intent(in) :: psil, psir
    j =    (gs - gp) * (psil%a(1)*psir%a(2) - psil%a(2)*psir%a(1)) &
         + (gs + gp) * (- psil%a(3)*psir%a(4) + psil%a(4)*psir%a(3))
  end function sp_ff
  pure function s_ff (gs, psil, psir) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gs
    type(bispinor), intent(in) :: psil, psir
    j = gs * (psil * psir)
  end function s_ff
  pure function p_ff (gp, psil, psir) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gp
    type(bispinor), intent(in) :: psil, psir
    j = gp * (- psil%a(1)*psir%a(2) + psil%a(2)*psir%a(1) &
              - psil%a(3)*psir%a(4) + psil%a(4)*psir%a(3))
  end function p_ff
  pure function sl_ff (gl, psil, psir) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gl
    type(bispinor), intent(in) :: psil, psir
    j =  2 * gl * (psil%a(1)*psir%a(2) - psil%a(2)*psir%a(1))
  end function sl_ff
  pure function sr_ff (gr, psil, psir) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gr
    type(bispinor), intent(in) :: psil, psir
    j = 2 * gr * (- psil%a(3)*psir%a(4) + psil%a(4)*psir%a(3))
  end function sr_ff
  pure function slr_ff (gl, gr, psibar, psi) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(bispinor), intent(in) :: psibar
    type(bispinor), intent(in) :: psi
    j = sp_ff (gr+gl, gr-gl, psibar, psi)
  end function slr_ff
  pure function f_spf (gs, gp, phi, psi) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: gs, gp
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    phipsi%a(1:2) = ((gs - gp) * phi) * psi%a(1:2)
    phipsi%a(3:4) = ((gs + gp) * phi) * psi%a(3:4)
  end function f_spf
  pure function f_sf (gs, phi, psi) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: gs
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    phipsi%a = (gs * phi) * psi%a
  end function f_sf
  pure function f_pf (gp, phi, psi) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: gp
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    phipsi%a(1:2) = (- gp * phi) * psi%a(1:2)
    phipsi%a(3:4) = (  gp * phi) * psi%a(3:4)
  end function f_pf
  pure function f_slf (gl, phi, psi) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: gl
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    phipsi%a(1:2) = (2 * gl * phi) * psi%a(1:2)
    phipsi%a(3:4) = 0
  end function f_slf
  pure function f_srf (gr, phi, psi) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: gr
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    phipsi%a(1:2) = 0
    phipsi%a(3:4) = (2 * gr * phi) * psi%a(3:4)
  end function f_srf
  pure function f_slrf (gl, gr, phi, psi) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: gl, gr
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    phipsi =  f_spf (gr+gl, gr-gl, phi, psi)
  end function f_slrf
  pure function vv_ff (psibar, psi, k) result (psibarpsi)
    type(vector) :: psibarpsi
    type(bispinor), intent(in) :: psibar, psi
    type(vector), intent(in) :: k
    complex(kind=default) :: kp, km, k12, k12s
    type(bispinor) :: kgpsi1, kgpsi2, kgpsi3, kgpsi4
    kp = k%t + k%x(3)
    km = k%t - k%x(3)
    k12  =  k%x(1) + (0,1)*k%x(2)
    k12s =  k%x(1) - (0,1)*k%x(2)
    kgpsi1%a(1) = -k%x(3) * psi%a(1) - k12s * psi%a(2)
    kgpsi1%a(2) = -k12 * psi%a(1) + k%x(3) * psi%a(2)
    kgpsi1%a(3) = k%x(3) * psi%a(3) + k12s * psi%a(4)
    kgpsi1%a(4) = k12 * psi%a(3) - k%x(3) * psi%a(4)
    kgpsi2%a(1) = ((0,-1) * k%x(2)) * psi%a(1) - km * psi%a(2)
    kgpsi2%a(2) = - kp * psi%a(1) + ((0,1) * k%x(2)) * psi%a(2)
    kgpsi2%a(3) = ((0,-1) * k%x(2)) * psi%a(3) + kp * psi%a(4)
    kgpsi2%a(4) = km * psi%a(3) + ((0,1) * k%x(2)) * psi%a(4)
    kgpsi3%a(1) = (0,1) * (k%x(1) * psi%a(1) + km * psi%a(2))
    kgpsi3%a(2) = (0,-1) * (kp * psi%a(1) + k%x(1) * psi%a(2))
    kgpsi3%a(3) = (0,1) * (k%x(1) * psi%a(3) - kp * psi%a(4))
    kgpsi3%a(4) = (0,1) * (km * psi%a(3) - k%x(1) * psi%a(4))
    kgpsi4%a(1) = -k%t * psi%a(1) - k12s * psi%a(2)
    kgpsi4%a(2) = k12 * psi%a(1) + k%t * psi%a(2)
    kgpsi4%a(3) = k%t * psi%a(3) - k12s * psi%a(4)
    kgpsi4%a(4) = k12 * psi%a(3) - k%t * psi%a(4)
    psibarpsi%t    = 2 * (psibar * kgpsi1)
    psibarpsi%x(1) = 2 * (psibar * kgpsi2)
    psibarpsi%x(2) = 2 * (psibar * kgpsi3)
    psibarpsi%x(3) = 2 * (psibar * kgpsi4)
  end function vv_ff
  pure function f_vvf (v, psi, k) result (kvpsi)
    type(bispinor) :: kvpsi
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: k, v
      complex(kind=default) :: kv30, kv21, kv01, kv31, kv02, kv32
    complex(kind=default) :: ap, am, bp, bm, bps, bms
    kv30 = k%x(3) * v%t - k%t * v%x(3)
    kv21 = (0,1) * (k%x(2) * v%x(1) - k%x(1) * v%x(2))
    kv01 = k%t * v%x(1) - k%x(1) * v%t
    kv31 = k%x(3) * v%x(1) - k%x(1) * v%x(3)
    kv02 = (0,1) * (k%t * v%x(2) - k%x(2) * v%t)
    kv32 = (0,1) * (k%x(3) * v%x(2) - k%x(2) * v%x(3))
    ap  = 2 * (kv30 + kv21)
    am  = 2 * (-kv30 + kv21)
    bp  = 2 * (kv01 + kv31 + kv02 + kv32)
    bm  = 2 * (kv01 - kv31 + kv02 - kv32)
    bps = 2 * (kv01 + kv31 - kv02 - kv32)
    bms = 2 * (kv01 - kv31 - kv02 + kv32)
    kvpsi%a(1) = am * psi%a(1) + bms * psi%a(2)
    kvpsi%a(2) = bp * psi%a(1) - am * psi%a(2)
    kvpsi%a(3) = ap * psi%a(3) - bps * psi%a(4)
    kvpsi%a(4) = -bm * psi%a(3) - ap * psi%a(4)
  end function f_vvf
  pure function vmom_ff (g, psibar, psi, k) result (psibarpsi)
    type(vector) :: psibarpsi
    complex(kind=default), intent(in) :: g
    type(bispinor), intent(in) :: psibar, psi
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    psibarpsi = g * vv_ff (psibar, psi, vk)
  end function vmom_ff
  pure function mom_ff (g, m, psibar, psi, k) result (psibarpsi)
    complex(kind=default) :: psibarpsi
    type(bispinor), intent(in) :: psibar, psi
    type(momentum), intent(in) :: k
    complex(kind=default), intent(in) :: g, m
    type(bispinor) :: kmpsi
    complex(kind=default) :: kp, km, k12, k12s
    kp = k%t + k%x(3)
    km = k%t - k%x(3)
    k12  =  k%x(1) + (0,1)*k%x(2)
    k12s =  k%x(1) - (0,1)*k%x(2)
    kmpsi%a(1) = km * psi%a(3) - k12s * psi%a(4)
    kmpsi%a(2) = kp * psi%a(4) - k12 * psi%a(3)
    kmpsi%a(3) = kp * psi%a(1) + k12s * psi%a(2)
    kmpsi%a(4) = k12 * psi%a(1) + km * psi%a(2)
    psibarpsi = g * (psibar * kmpsi) + s_ff (m, psibar, psi)
  end function mom_ff
  pure function mom5_ff (g, m, psibar, psi, k) result (psibarpsi)
    complex(kind=default) :: psibarpsi
    type(bispinor), intent(in) :: psibar, psi
    type(momentum), intent(in) :: k
    complex(kind=default), intent(in) :: g, m
    type(bispinor) :: g5psi
    g5psi%a(1:2) = - psi%a(1:2)
    g5psi%a(3:4) = psi%a(3:4)
    psibarpsi = mom_ff (g, m, psibar, g5psi, k)
  end function mom5_ff
  pure function moml_ff (g, m, psibar, psi, k) result (psibarpsi)
    complex(kind=default) :: psibarpsi
    type(bispinor), intent(in) :: psibar, psi
    type(momentum), intent(in) :: k
    complex(kind=default), intent(in) :: g, m
    type(bispinor) :: leftpsi
    leftpsi%a(1:2) = 2 * psi%a(1:2)
    leftpsi%a(3:4) = 0
    psibarpsi = mom_ff (g, m, psibar, leftpsi, k)
  end function moml_ff
  pure function momr_ff (g, m, psibar, psi, k) result (psibarpsi)
    complex(kind=default) :: psibarpsi
    type(bispinor), intent(in) :: psibar, psi
    type(momentum), intent(in) :: k
    complex(kind=default), intent(in) :: g, m
    type(bispinor) :: rightpsi
    rightpsi%a(1:2) = 0
    rightpsi%a(3:4) = 2 * psi%a(3:4)
    psibarpsi = mom_ff (g, m, psibar, rightpsi, k)
  end function momr_ff
  pure function lmom_ff (g, m, psibar, psi, k) result (psibarpsi)
    complex(kind=default) :: psibarpsi
    type(bispinor), intent(in) :: psibar, psi
    type(momentum), intent(in) :: k
    complex(kind=default), intent(in) :: g, m
    psibarpsi = mom_ff  (g, m, psibar, psi, k) + &
                mom5_ff (g,-m, psibar, psi, k)
  end function lmom_ff
  pure function rmom_ff (g, m, psibar, psi, k) result (psibarpsi)
    complex(kind=default) :: psibarpsi
    type(bispinor), intent(in) :: psibar, psi
    type(momentum), intent(in) :: k
    complex(kind=default), intent(in) :: g, m
    psibarpsi = mom_ff  (g, m, psibar, psi, k) - &
                mom5_ff (g,-m, psibar, psi, k)
  end function rmom_ff
  pure function f_vmomf (g, v, psi, k) result (kvpsi)
    type(bispinor) :: kvpsi
    type(bispinor), intent(in) :: psi
    complex(kind=default), intent(in) :: g
    type(momentum), intent(in) :: k
    type(vector), intent(in) :: v
    type(vector) :: vk
    vk = k
    kvpsi = g * f_vvf (v, psi, vk)
  end function f_vmomf
  pure function f_momf (g, m, phi, psi, k) result (kmpsi)
    type(bispinor) :: kmpsi
    type(bispinor), intent(in) :: psi
    complex(kind=default), intent(in) :: phi, g, m
    type(momentum), intent(in) :: k
    complex(kind=default) :: kp, km, k12, k12s
    kp = k%t + k%x(3)
    km = k%t - k%x(3)
    k12  =  k%x(1) + (0,1)*k%x(2)
    k12s =  k%x(1) - (0,1)*k%x(2)
    kmpsi%a(1) = km * psi%a(3) - k12s * psi%a(4)
    kmpsi%a(2) = -k12 * psi%a(3) + kp * psi%a(4)
    kmpsi%a(3) = kp * psi%a(1) + k12s * psi%a(2)
    kmpsi%a(4) = k12 * psi%a(1) + km * psi%a(2)
    kmpsi = g * (phi * kmpsi) + f_sf (m, phi, psi)
  end function f_momf
  pure function f_mom5f (g, m, phi, psi, k) result (kmpsi)
    type(bispinor) :: kmpsi
    type(bispinor), intent(in) :: psi
    complex(kind=default), intent(in) :: phi, g, m
    type(momentum), intent(in) :: k
    type(bispinor) :: g5psi
    g5psi%a(1:2) = - psi%a(1:2)
    g5psi%a(3:4) =   psi%a(3:4)
    kmpsi = f_momf (g, m, phi, g5psi, k)
  end function f_mom5f
  pure function f_momlf (g, m, phi, psi, k) result (kmpsi)
    type(bispinor) :: kmpsi
    type(bispinor), intent(in) :: psi
    complex(kind=default), intent(in) :: phi, g, m
    type(momentum), intent(in) :: k
    type(bispinor) :: leftpsi
    leftpsi%a(1:2) = 2 * psi%a(1:2)
    leftpsi%a(3:4) = 0
    kmpsi = f_momf (g, m, phi, leftpsi, k)
  end function f_momlf
  pure function f_momrf (g, m, phi, psi, k) result (kmpsi)
    type(bispinor) :: kmpsi
    type(bispinor), intent(in) :: psi
    complex(kind=default), intent(in) :: phi, g, m
    type(momentum), intent(in) :: k
    type(bispinor) :: rightpsi
    rightpsi%a(1:2) = 0
    rightpsi%a(3:4) = 2 * psi%a(3:4)
    kmpsi = f_momf (g, m, phi, rightpsi, k)
  end function f_momrf
  pure function f_lmomf (g, m, phi, psi, k) result (kmpsi)
    type(bispinor) :: kmpsi
    type(bispinor), intent(in) :: psi
    complex(kind=default), intent(in) :: phi, g, m
    type(momentum), intent(in) :: k
    kmpsi = f_momf  (g, m, phi, psi, k) + &
            f_mom5f (g,-m, phi, psi, k)
  end function f_lmomf
  pure function f_rmomf (g, m, phi, psi, k) result (kmpsi)
    type(bispinor) :: kmpsi
    type(bispinor), intent(in) :: psi
    complex(kind=default), intent(in) :: phi, g, m
    type(momentum), intent(in) :: k
    kmpsi = f_momf  (g, m, phi, psi, k) - &
            f_mom5f (g,-m, phi, psi, k)
  end function f_rmomf
  pure function v2_ff (g, psibar, v, psi) result (v2)
    type(vector) :: v2
    complex (kind=default), intent(in) :: g
    type(bispinor), intent(in) :: psibar, psi
    type(vector), intent(in) :: v
    v2 = (-g) * vv_ff (psibar, psi, v)
  end function v2_ff
  pure function sv1_ff (g, psibar, v, psi) result (phi)
    complex(kind=default) :: phi
    type(bispinor), intent(in) :: psibar, psi
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: g
    phi = psibar * f_vf (g, v, psi)
  end function sv1_ff
  pure function sv2_ff (g, psibar, phi, psi) result (v)
    type(vector) :: v
    complex(kind=default), intent(in) :: phi, g
    type(bispinor), intent(in) :: psibar, psi
    v = phi * v_ff (g, psibar, psi)
  end function sv2_ff
  pure function pv1_ff (g, psibar, v, psi) result (phi)
    complex(kind=default) :: phi
    type(bispinor), intent(in) :: psibar, psi
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: g
    phi = - (psibar * f_af (g, v, psi))
  end function pv1_ff
  pure function pv2_ff (g, psibar, phi, psi) result (v)
    type(vector) :: v
    complex(kind=default), intent(in) :: phi, g
    type(bispinor), intent(in) :: psibar, psi
    v = -(phi * a_ff (g, psibar, psi))
  end function pv2_ff
  pure function svl1_ff (g, psibar, v, psi) result (phi)
    complex(kind=default) :: phi
    type(bispinor), intent(in) :: psibar, psi
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: g
    phi = psibar * f_vlf (g, v, psi)
  end function svl1_ff
  pure function svl2_ff (g, psibar, phi, psi) result (v)
    type(vector) :: v
    complex(kind=default), intent(in) :: phi, g
    type(bispinor), intent(in) :: psibar, psi
    v = phi * vl_ff (g, psibar, psi)
  end function svl2_ff
  pure function svr1_ff (g, psibar, v, psi) result (phi)
    complex(kind=default) :: phi
    type(bispinor), intent(in) :: psibar, psi
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: g
    phi = psibar * f_vrf (g, v, psi)
  end function svr1_ff
  pure function svr2_ff (g, psibar, phi, psi) result (v)
    type(vector) :: v
    complex(kind=default), intent(in) :: phi, g
    type(bispinor), intent(in) :: psibar, psi
    v = phi * vr_ff (g, psibar, psi)
  end function svr2_ff
  pure function svlr1_ff (gl, gr, psibar, v, psi) result (phi)
    complex(kind=default) :: phi
    type(bispinor), intent(in) :: psibar, psi
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: gl, gr
    phi = psibar * f_vlrf (gl, gr, v, psi)
  end function svlr1_ff
  pure function svlr2_ff (gl, gr, psibar, phi, psi) result (v)
    type(vector) :: v
    complex(kind=default), intent(in) :: phi, gl, gr
    type(bispinor), intent(in) :: psibar, psi
    v = phi * vlr_ff (gl, gr, psibar, psi)
  end function svlr2_ff
  pure function f_v2f (g, v1, v2, psi) result (vpsi)
    type(bispinor) :: vpsi
    complex(kind=default), intent(in) :: g
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v1, v2
    vpsi = g * f_vvf (v2, psi, v1)
  end function f_v2f
  pure function f_svf (g, phi, v, psi) result (pvpsi)
    type(bispinor) :: pvpsi
    complex(kind=default), intent(in) :: g, phi
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    pvpsi = phi * f_vf (g, v, psi)
  end function f_svf
  pure function f_pvf (g, phi, v, psi) result (pvpsi)
    type(bispinor) :: pvpsi
    complex(kind=default), intent(in) :: g, phi
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    pvpsi = -(phi * f_af (g, v, psi))
  end function f_pvf
  pure function f_svlf (g, phi, v, psi) result (pvpsi)
    type(bispinor) :: pvpsi
    complex(kind=default), intent(in) :: g, phi
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    pvpsi = phi * f_vlf (g, v, psi)
  end function f_svlf
  pure function f_svrf (g, phi, v, psi) result (pvpsi)
    type(bispinor) :: pvpsi
    complex(kind=default), intent(in) :: g, phi
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    pvpsi = phi * f_vrf (g, v, psi)
  end function f_svrf
  pure function f_svlrf (gl, gr, phi, v, psi) result (pvpsi)
    type(bispinor) :: pvpsi
    complex(kind=default), intent(in) :: gl, gr, phi
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    pvpsi = phi * f_vlrf (gl, gr, v, psi)
  end function f_svlrf
  pure function pot_grf (g, gravbar, psi) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(vectorspinor) :: gamma_psi
    gamma_psi%psi(1)%a(1) = psi%a(3)
    gamma_psi%psi(1)%a(2) = psi%a(4)
    gamma_psi%psi(1)%a(3) = psi%a(1)
    gamma_psi%psi(1)%a(4) = psi%a(2)
    gamma_psi%psi(2)%a(1) = psi%a(4)
    gamma_psi%psi(2)%a(2) = psi%a(3)
    gamma_psi%psi(2)%a(3) = - psi%a(2)
    gamma_psi%psi(2)%a(4) = - psi%a(1)
    gamma_psi%psi(3)%a(1) = (0,-1) * psi%a(4)
    gamma_psi%psi(3)%a(2) = (0,1) * psi%a(3)
    gamma_psi%psi(3)%a(3) = (0,1) * psi%a(2)
    gamma_psi%psi(3)%a(4) = (0,-1) * psi%a(1)
    gamma_psi%psi(4)%a(1) = psi%a(3)
    gamma_psi%psi(4)%a(2) = - psi%a(4)
    gamma_psi%psi(4)%a(3) = - psi%a(1)
    gamma_psi%psi(4)%a(4) = psi%a(2)
    j = g * (gravbar * gamma_psi)
  end function pot_grf
  pure function pot_fgr (g, psibar, grav) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g
    type(bispinor), intent(in) :: psibar
    type(vectorspinor), intent(in) :: grav
    type(bispinor) :: gamma_grav
    gamma_grav%a(1) = grav%psi(1)%a(3) - grav%psi(2)%a(4) + &
           ((0,1)*grav%psi(3)%a(4)) - grav%psi(4)%a(3)
    gamma_grav%a(2) = grav%psi(1)%a(4) - grav%psi(2)%a(3) - &
           ((0,1)*grav%psi(3)%a(3)) + grav%psi(4)%a(4)
    gamma_grav%a(3) = grav%psi(1)%a(1) + grav%psi(2)%a(2) - &
           ((0,1)*grav%psi(3)%a(2)) + grav%psi(4)%a(1)
    gamma_grav%a(4) = grav%psi(1)%a(2) + grav%psi(2)%a(1) + &
           ((0,1)*grav%psi(3)%a(1)) - grav%psi(4)%a(2)
    j = g * (psibar * gamma_grav)
  end function pot_fgr
  pure function grvgf (gravbar, psi, k) result (j)
    complex(kind=default) :: j
    complex(kind=default) :: kp, km, k12, k12s
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: k
    type(vectorspinor) :: kg_psi
    kp = k%t + k%x(3)
    km = k%t - k%x(3)
    k12  =  k%x(1) + (0,1)*k%x(2)
    k12s =  k%x(1) - (0,1)*k%x(2)
    !!! Since we are taking the spinor product here, NO explicit
    !!! charge conjugation matrix is needed!
    kg_psi%psi(1)%a(1) = km * psi%a(1) - k12s * psi%a(2)
    kg_psi%psi(1)%a(2) = (-k12) * psi%a(1) + kp * psi%a(2)
    kg_psi%psi(1)%a(3) = kp * psi%a(3) + k12s * psi%a(4)
    kg_psi%psi(1)%a(4) = k12 * psi%a(3) + km * psi%a(4)
    kg_psi%psi(2)%a(1) = k12s * psi%a(1) - km * psi%a(2)
    kg_psi%psi(2)%a(2) = (-kp) * psi%a(1) + k12 * psi%a(2)
    kg_psi%psi(2)%a(3) = k12s * psi%a(3) + kp * psi%a(4)
    kg_psi%psi(2)%a(4) = km * psi%a(3) + k12 * psi%a(4)
    kg_psi%psi(3)%a(1) = (0,1) * (k12s * psi%a(1) + km * psi%a(2))
    kg_psi%psi(3)%a(2) = (0,1) * (- kp * psi%a(1) - k12 * psi%a(2))
    kg_psi%psi(3)%a(3) = (0,1) * (k12s * psi%a(3) - kp * psi%a(4))
    kg_psi%psi(3)%a(4) = (0,1) * (km * psi%a(3) - k12 * psi%a(4))
    kg_psi%psi(4)%a(1) = (-km) * psi%a(1) - k12s * psi%a(2)
    kg_psi%psi(4)%a(2) = k12 * psi%a(1) + kp * psi%a(2)
    kg_psi%psi(4)%a(3) = kp * psi%a(3) - k12s * psi%a(4)
    kg_psi%psi(4)%a(4) = k12 * psi%a(3) - km * psi%a(4)
    j = gravbar * kg_psi
  end function grvgf
  pure function grg5vgf (gravbar, psi, k) result (j)
    complex(kind=default) :: j
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: k
    type(bispinor) :: g5_psi
    g5_psi%a(1:2) = - psi%a(1:2)
    g5_psi%a(3:4) =   psi%a(3:4)
    j = grvgf (gravbar, g5_psi, k)
  end function grg5vgf
  pure function s_grf (g, gravbar, psi, k) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    j = g * grvgf (gravbar, psi, vk)
  end function s_grf
  pure function sl_grf (gl, gravbar, psi, k) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gl
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_l
    type(momentum), intent(in) :: k
    psi_l%a(1:2) = psi%a(1:2)
    psi_l%a(3:4) = 0
    j = s_grf (gl, gravbar, psi_l, k)
  end function sl_grf
  pure function sr_grf (gr, gravbar, psi, k) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gr
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_r
    type(momentum), intent(in) :: k
    psi_r%a(1:2) = 0
    psi_r%a(3:4) = psi%a(3:4)
    j = s_grf (gr, gravbar, psi_r, k)
  end function sr_grf
  pure function slr_grf (gl, gr, gravbar, psi, k) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(momentum), intent(in) :: k
    j = sl_grf (gl, gravbar, psi, k) + sr_grf (gr, gravbar, psi, k)
  end function slr_grf
  pure function fgkgr (psibar, grav, k) result (j)
    complex(kind=default) :: j
    complex(kind=default) :: kp, km, k12, k12s
    type(bispinor), intent(in) :: psibar
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: k
    type(bispinor) :: gk_grav
    kp = k%t + k%x(3)
    km = k%t - k%x(3)
    k12  =  k%x(1) + (0,1)*k%x(2)
    k12s =  k%x(1) - (0,1)*k%x(2)
    !!! Since we are taking the spinor product here, NO explicit
    !!! charge conjugation matrix is needed!
    gk_grav%a(1) =  kp * grav%psi(1)%a(1) + k12s * grav%psi(1)%a(2) &
                 - k12 * grav%psi(2)%a(1) - km * grav%psi(2)%a(2) &
                 + (0,1) * k12 * grav%psi(3)%a(1)   &
                 + (0,1) * km * grav%psi(3)%a(2) &
                 - kp * grav%psi(4)%a(1) - k12s * grav%psi(4)%a(2)
    gk_grav%a(2) = k12 * grav%psi(1)%a(1) + km * grav%psi(1)%a(2) &
                 - kp * grav%psi(2)%a(1) - k12s * grav%psi(2)%a(2) &
                 - (0,1) * kp * grav%psi(3)%a(1) &
                 - (0,1) * k12s * grav%psi(3)%a(2)  &
                 + k12 * grav%psi(4)%a(1) + km * grav%psi(4)%a(2)
    gk_grav%a(3) = km * grav%psi(1)%a(3) - k12s * grav%psi(1)%a(4) &
                 - k12 * grav%psi(2)%a(3) + kp * grav%psi(2)%a(4) &
                 + (0,1) * k12 * grav%psi(3)%a(3)   &
                 - (0,1) * kp * grav%psi(3)%a(4) &
                 + km * grav%psi(4)%a(3) - k12s * grav%psi(4)%a(4)
    gk_grav%a(4) = - k12 * grav%psi(1)%a(3) + kp * grav%psi(1)%a(4) &
                 + km * grav%psi(2)%a(3) - k12s * grav%psi(2)%a(4) &
                 + (0,1) * km * grav%psi(3)%a(3) &
                 - (0,1) * k12s * grav%psi(3)%a(4)  &
                 + k12 * grav%psi(4)%a(3) - kp * grav%psi(4)%a(4)
    j = psibar * gk_grav
  end function fgkgr
  pure function fg5gkgr (psibar, grav, k) result (j)
    complex(kind=default) :: j
    type(bispinor), intent(in) :: psibar
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: k
    type(bispinor) :: psibar_g5
    psibar_g5%a(1:2) = - psibar%a(1:2)
    psibar_g5%a(3:4) =   psibar%a(3:4)
    j = fgkgr (psibar_g5, grav, k)
  end function fg5gkgr
  pure function s_fgr (g, psibar, grav, k) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g
    type(bispinor), intent(in) :: psibar
    type(vectorspinor), intent(in) :: grav
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    j = g * fgkgr (psibar, grav, vk)
  end function s_fgr
  pure function sl_fgr (gl, psibar, grav, k) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gl
    type(bispinor), intent(in) :: psibar
    type(bispinor) :: psibar_l
    type(vectorspinor), intent(in) :: grav
    type(momentum), intent(in) :: k
    psibar_l%a(1:2) = psibar%a(1:2)
    psibar_l%a(3:4) = 0
    j = s_fgr (gl, psibar_l, grav, k)
  end function sl_fgr
  pure function sr_fgr (gr, psibar, grav, k) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gr
    type(bispinor), intent(in) :: psibar
    type(bispinor) :: psibar_r
    type(vectorspinor), intent(in) :: grav
    type(momentum), intent(in) :: k
    psibar_r%a(1:2) = 0
    psibar_r%a(3:4) = psibar%a(3:4)
    j = s_fgr (gr, psibar_r, grav, k)
  end function sr_fgr
  pure function slr_fgr (gl, gr, psibar, grav, k) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(bispinor), intent(in) :: psibar
    type(vectorspinor), intent(in) :: grav
    type(momentum), intent(in) :: k
    j = sl_fgr (gl, psibar, grav, k) + sr_fgr (gr, psibar, grav, k)
  end function slr_fgr
  pure function p_grf (g, gravbar, psi, k) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    j = g * grg5vgf (gravbar, psi, vk)
  end function p_grf
  pure function p_fgr (g, psibar, grav, k) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g
    type(bispinor), intent(in) :: psibar
    type(vectorspinor), intent(in) :: grav
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    j = g * fg5gkgr (psibar, grav, vk)
  end function p_fgr
  pure function f_potgr (g, phi, psi) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: phi
    type(vectorspinor), intent(in) :: psi
    phipsi%a(1) = (g * phi) * (psi%psi(1)%a(3) - psi%psi(2)%a(4) + &
                  ((0,1)*psi%psi(3)%a(4)) - psi%psi(4)%a(3))
    phipsi%a(2) = (g * phi) * (psi%psi(1)%a(4) - psi%psi(2)%a(3) - &
                  ((0,1)*psi%psi(3)%a(3)) + psi%psi(4)%a(4))
    phipsi%a(3) = (g * phi) * (psi%psi(1)%a(1) + psi%psi(2)%a(2) - &
                  ((0,1)*psi%psi(3)%a(2)) + psi%psi(4)%a(1))
    phipsi%a(4) = (g * phi) * (psi%psi(1)%a(2) + psi%psi(2)%a(1) + &
                  ((0,1)*psi%psi(3)%a(1)) - psi%psi(4)%a(2))
  end function f_potgr
  pure function fgvgr (psi, k) result (kpsi)
    type(bispinor) :: kpsi
    complex(kind=default) :: kp, km, k12, k12s
    type(vector), intent(in) :: k
    type(vectorspinor), intent(in) :: psi
    kp = k%t + k%x(3)
    km = k%t - k%x(3)
    k12  =  k%x(1) + (0,1)*k%x(2)
    k12s =  k%x(1) - (0,1)*k%x(2)
    kpsi%a(1) = kp * psi%psi(1)%a(1) + k12s * psi%psi(1)%a(2) &
              - k12 * psi%psi(2)%a(1) - km * psi%psi(2)%a(2) &
              + (0,1) * k12 * psi%psi(3)%a(1) + (0,1) * km * psi%psi(3)%a(2) &
              - kp * psi%psi(4)%a(1) - k12s * psi%psi(4)%a(2)
    kpsi%a(2) = k12 * psi%psi(1)%a(1) + km * psi%psi(1)%a(2) &
              - kp * psi%psi(2)%a(1) - k12s * psi%psi(2)%a(2) &
              - (0,1) * kp * psi%psi(3)%a(1) - (0,1) * k12s * psi%psi(3)%a(2) &
              + k12 * psi%psi(4)%a(1) + km * psi%psi(4)%a(2)
    kpsi%a(3) = km * psi%psi(1)%a(3) - k12s * psi%psi(1)%a(4) &
              - k12 * psi%psi(2)%a(3) + kp * psi%psi(2)%a(4) &
              + (0,1) * k12 * psi%psi(3)%a(3) - (0,1) * kp * psi%psi(3)%a(4) &
              + km * psi%psi(4)%a(3) - k12s * psi%psi(4)%a(4)
    kpsi%a(4) = - k12 * psi%psi(1)%a(3) + kp * psi%psi(1)%a(4) &
              + km * psi%psi(2)%a(3) - k12s * psi%psi(2)%a(4) &
              + (0,1) * km * psi%psi(3)%a(3) - (0,1) * k12s * psi%psi(3)%a(4) &
              + k12 * psi%psi(4)%a(3) - kp * psi%psi(4)%a(4)
  end function fgvgr
  pure function f_sgr (g, phi, psi, k) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: phi
    type(momentum), intent(in) :: k
    type(vectorspinor), intent(in) :: psi
    type(vector) :: vk
    vk = k
    phipsi = (g * phi) * fgvgr (psi, vk)
  end function f_sgr
  pure function f_slgr (gl, phi, psi, k) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: gl
    complex(kind=default), intent(in) :: phi
    type(momentum), intent(in) :: k
    type(vectorspinor), intent(in) :: psi
    phipsi = f_sgr (gl, phi, psi, k)
    phipsi%a(3:4) = 0
  end function f_slgr
  pure function f_srgr (gr, phi, psi, k) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: gr
    complex(kind=default), intent(in) :: phi
    type(momentum), intent(in) :: k
    type(vectorspinor), intent(in) :: psi
    phipsi = f_sgr (gr, phi, psi, k)
    phipsi%a(1:2) = 0
  end function f_srgr
  pure function f_slrgr (gl, gr, phi, psi, k) result (phipsi)
    type(bispinor) :: phipsi, phipsi_l, phipsi_r
    complex(kind=default), intent(in) :: gl, gr
    complex(kind=default), intent(in) :: phi
    type(momentum), intent(in) :: k
    type(vectorspinor), intent(in) :: psi
    phipsi_l = f_slgr (gl, phi, psi, k)
    phipsi_r = f_srgr (gr, phi, psi, k)
    phipsi%a(1:2) = phipsi_l%a(1:2)
    phipsi%a(3:4) = phipsi_r%a(3:4)
  end function f_slrgr
  pure function fgvg5gr (psi, k) result (kpsi)
    type(bispinor) :: kpsi
    type(vector), intent(in) :: k
    type(vectorspinor), intent(in) :: psi
    type(bispinor) :: kpsi_dum
    kpsi_dum = fgvgr (psi, k)
    kpsi%a(1:2) = - kpsi_dum%a(1:2)
    kpsi%a(3:4) =   kpsi_dum%a(3:4)
  end function fgvg5gr
  pure function f_pgr (g, phi, psi, k) result (phipsi)
    type(bispinor) :: phipsi
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: phi
    type(momentum), intent(in) :: k
    type(vectorspinor), intent(in) :: psi
    type(vector) :: vk
    vk = k
    phipsi = (g * phi) * fgvg5gr (psi, vk)
  end function f_pgr
  pure function fggvvgr (v, psi, k) result (psikv)
    type(bispinor) :: psikv
    type(vectorspinor), intent(in) :: psi
    type(vector), intent(in) :: v, k
    complex(kind=default) :: kv30, kv21, kv01, kv31, kv02, kv32
    complex(kind=default) :: ap, am, bp, bm, bps, bms
    kv30 = k%x(3) * v%t - k%t * v%x(3)
    kv21 = (0,1) * (k%x(2) * v%x(1) - k%x(1) * v%x(2))
    kv01 = k%t * v%x(1) - k%x(1) * v%t
    kv31 = k%x(3) * v%x(1) - k%x(1) * v%x(3)
    kv02 = (0,1) * (k%t * v%x(2) - k%x(2) * v%t)
    kv32 = (0,1) * (k%x(3) * v%x(2) - k%x(2) * v%x(3))
    ap  = 2 * (kv30 + kv21)
    am  = 2 * (-kv30 + kv21)
    bp  = 2 * (kv01 + kv31 + kv02 + kv32)
    bm  = 2 * (kv01 - kv31 + kv02 - kv32)
    bps = 2 * (kv01 + kv31 - kv02 - kv32)
    bms = 2 * (kv01 - kv31 - kv02 + kv32)
    psikv%a(1) = (-ap) * psi%psi(1)%a(3) + bps * psi%psi(1)%a(4) &
               + (-bm) * psi%psi(2)%a(3) + (-ap) * psi%psi(2)%a(4) &
               + (0,1) * (bm * psi%psi(3)%a(3) + ap * psi%psi(3)%a(4)) &
               + ap * psi%psi(4)%a(3) + (-bps) * psi%psi(4)%a(4)
    psikv%a(2) =  bm * psi%psi(1)%a(3) + ap * psi%psi(1)%a(4) &
               + ap * psi%psi(2)%a(3) + (-bps) * psi%psi(2)%a(4) &
               + (0,1) * (ap * psi%psi(3)%a(3) - bps * psi%psi(3)%a(4)) &
               + bm * psi%psi(4)%a(3) + ap * psi%psi(4)%a(4)
    psikv%a(3) =  am * psi%psi(1)%a(1) + bms * psi%psi(1)%a(2) &
               + bp * psi%psi(2)%a(1) + (-am) * psi%psi(2)%a(2) &
               + (0,-1) * (bp * psi%psi(3)%a(1) + (-am) * psi%psi(3)%a(2)) &
               + am * psi%psi(4)%a(1) + bms * psi%psi(4)%a(2)
    psikv%a(4) =  bp * psi%psi(1)%a(1) + (-am) * psi%psi(1)%a(2) &
               + am * psi%psi(2)%a(1) + bms * psi%psi(2)%a(2) &
               + (0,1) * (am * psi%psi(3)%a(1) + bms * psi%psi(3)%a(2)) &
               + (-bp) * psi%psi(4)%a(1) + am * psi%psi(4)%a(2)
  end function fggvvgr
  pure function f_vgr (g, v, psi, k) result (psikkkv)
    type(bispinor) :: psikkkv
    type(vectorspinor), intent(in) :: psi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: k
    complex(kind=default), intent(in) :: g
    type(vector) :: vk
    vk = k
    psikkkv = g * (fggvvgr (v, psi, vk))
  end function f_vgr
  pure function f_vlrgr (gl, gr, v, psi, k) result (psikv)
    type(bispinor) :: psikv
    type(vectorspinor), intent(in) :: psi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: k
    complex(kind=default), intent(in) :: gl, gr
    type(vector) :: vk
    vk = k
    psikv = fggvvgr (v, psi, vk)
    psikv%a(1:2) = gl * psikv%a(1:2)
    psikv%a(3:4) = gr * psikv%a(3:4)
  end function f_vlrgr
  pure function gr_potf (g, phi, psi) result (phipsi)
    type(vectorspinor) :: phipsi
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    phipsi%psi(1)%a(1) = (g * phi) * psi%a(3)
    phipsi%psi(1)%a(2) = (g * phi) * psi%a(4)
    phipsi%psi(1)%a(3) = (g * phi) * psi%a(1)
    phipsi%psi(1)%a(4) = (g * phi) * psi%a(2)
    phipsi%psi(2)%a(1) = (g * phi) * psi%a(4)
    phipsi%psi(2)%a(2) = (g * phi) * psi%a(3)
    phipsi%psi(2)%a(3) = ((-g) * phi) * psi%a(2)
    phipsi%psi(2)%a(4) = ((-g) * phi) * psi%a(1)
    phipsi%psi(3)%a(1) = ((0,-1) * g * phi) * psi%a(4)
    phipsi%psi(3)%a(2) = ((0,1) * g * phi) * psi%a(3)
    phipsi%psi(3)%a(3) = ((0,1) * g * phi) * psi%a(2)
    phipsi%psi(3)%a(4) = ((0,-1) * g * phi) * psi%a(1)
    phipsi%psi(4)%a(1) = (g * phi) * psi%a(3)
    phipsi%psi(4)%a(2) = ((-g) * phi) * psi%a(4)
    phipsi%psi(4)%a(3) = ((-g) * phi) * psi%a(1)
    phipsi%psi(4)%a(4) = (g * phi) * psi%a(2)
  end function gr_potf
  pure function grkgf (psi, k) result (kpsi)
    type(vectorspinor) :: kpsi
    complex(kind=default) :: kp, km, k12, k12s
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: k
    kp = k%t + k%x(3)
    km = k%t - k%x(3)
    k12  =  k%x(1) + (0,1)*k%x(2)
    k12s =  k%x(1) - (0,1)*k%x(2)
    kpsi%psi(1)%a(1) = km * psi%a(1) - k12s * psi%a(2)
    kpsi%psi(1)%a(2) = (-k12) * psi%a(1) + kp * psi%a(2)
    kpsi%psi(1)%a(3) = kp * psi%a(3) + k12s * psi%a(4)
    kpsi%psi(1)%a(4) = k12 * psi%a(3) + km * psi%a(4)
    kpsi%psi(2)%a(1) = k12s * psi%a(1) - km * psi%a(2)
    kpsi%psi(2)%a(2) = (-kp) * psi%a(1) + k12 * psi%a(2)
    kpsi%psi(2)%a(3) = k12s * psi%a(3) + kp * psi%a(4)
    kpsi%psi(2)%a(4) = km * psi%a(3) + k12 * psi%a(4)
    kpsi%psi(3)%a(1) = (0,1) * (k12s * psi%a(1) + km * psi%a(2))
    kpsi%psi(3)%a(2) = (0,-1) * (kp * psi%a(1) + k12 * psi%a(2))
    kpsi%psi(3)%a(3) = (0,1) * (k12s * psi%a(3) - kp * psi%a(4))
    kpsi%psi(3)%a(4) = (0,1) * (km * psi%a(3) - k12 * psi%a(4))
    kpsi%psi(4)%a(1) = -(km * psi%a(1) + k12s * psi%a(2))
    kpsi%psi(4)%a(2) = k12 * psi%a(1) + kp * psi%a(2)
    kpsi%psi(4)%a(3) = kp * psi%a(3) - k12s * psi%a(4)
    kpsi%psi(4)%a(4) = k12 * psi%a(3) - km * psi%a(4)
  end function grkgf
  pure function gr_sf (g, phi, psi, k) result (phipsi)
    type(vectorspinor) :: phipsi
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    phipsi = (g * phi) * grkgf (psi, vk)
  end function gr_sf
  pure function gr_slf (gl, phi, psi, k) result (phipsi)
    type(vectorspinor) :: phipsi
    complex(kind=default), intent(in) :: gl
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_l
    type(momentum), intent(in) :: k
    psi_l%a(1:2) = psi%a(1:2)
    psi_l%a(3:4) = 0
    phipsi = gr_sf (gl, phi, psi_l, k)
  end function gr_slf
  pure function gr_srf (gr, phi, psi, k) result (phipsi)
    type(vectorspinor) :: phipsi
    complex(kind=default), intent(in) :: gr
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_r
    type(momentum), intent(in) :: k
    psi_r%a(1:2) = 0
    psi_r%a(3:4) = psi%a(3:4)
    phipsi = gr_sf (gr, phi, psi_r, k)
  end function gr_srf
  pure function gr_slrf (gl, gr, phi, psi, k) result (phipsi)
    type(vectorspinor) :: phipsi
    complex(kind=default), intent(in) :: gl, gr
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    type(momentum), intent(in) :: k
    phipsi = gr_slf (gl, phi, psi, k) + gr_srf (gr, phi, psi, k)
  end function gr_slrf
  pure function grkggf (psi, k) result (kpsi)
    type(vectorspinor) :: kpsi
    complex(kind=default) :: kp, km, k12, k12s
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: k
    kp = k%t + k%x(3)
    km = k%t - k%x(3)
    k12  =  k%x(1) + (0,1)*k%x(2)
    k12s =  k%x(1) - (0,1)*k%x(2)
    kpsi%psi(1)%a(1) = (-km) * psi%a(1) + k12s * psi%a(2)
    kpsi%psi(1)%a(2) = k12 * psi%a(1) - kp * psi%a(2)
    kpsi%psi(1)%a(3) = kp * psi%a(3) + k12s * psi%a(4)
    kpsi%psi(1)%a(4) = k12 * psi%a(3) + km * psi%a(4)
    kpsi%psi(2)%a(1) = (-k12s) * psi%a(1) + km * psi%a(2)
    kpsi%psi(2)%a(2) = kp * psi%a(1) - k12 * psi%a(2)
    kpsi%psi(2)%a(3) = k12s * psi%a(3) + kp * psi%a(4)
    kpsi%psi(2)%a(4) = km * psi%a(3) + k12 * psi%a(4)
    kpsi%psi(3)%a(1) = (0,-1) * (k12s * psi%a(1) + km * psi%a(2))
    kpsi%psi(3)%a(2) = (0,1) * (kp * psi%a(1) + k12 * psi%a(2))
    kpsi%psi(3)%a(3) = (0,1) * (k12s * psi%a(3) - kp * psi%a(4))
    kpsi%psi(3)%a(4) = (0,1) * (km * psi%a(3) - k12 * psi%a(4))
    kpsi%psi(4)%a(1) = km * psi%a(1) + k12s * psi%a(2)
    kpsi%psi(4)%a(2) = -(k12 * psi%a(1) + kp * psi%a(2))
    kpsi%psi(4)%a(3) = kp * psi%a(3) - k12s * psi%a(4)
    kpsi%psi(4)%a(4) = k12 * psi%a(3) - km * psi%a(4)
  end function grkggf
  pure function gr_pf (g, phi, psi, k) result (phipsi)
    type(vectorspinor) :: phipsi
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: phi
    type(bispinor), intent(in) :: psi
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    phipsi = (g * phi) * grkggf (psi, vk)
  end function gr_pf
  pure function grkkggf (v, psi, k) result (psikv)
    type(vectorspinor) :: psikv
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v, k
    complex(kind=default) :: kv30, kv21, kv01, kv31, kv02, kv32
    complex(kind=default) :: ap, am, bp, bm, bps, bms, imago
    imago = (0.0_default,1.0_default)
    kv30 = k%x(3) * v%t - k%t * v%x(3)
    kv21 = imago * (k%x(2) * v%x(1) - k%x(1) * v%x(2))
    kv01 = k%t * v%x(1) - k%x(1) * v%t
    kv31 = k%x(3) * v%x(1) - k%x(1) * v%x(3)
    kv02 = imago * (k%t * v%x(2) - k%x(2) * v%t)
    kv32 = imago * (k%x(3) * v%x(2) - k%x(2) * v%x(3))
    ap  = 2 * (kv30 + kv21)
    am  = 2 * ((-kv30) + kv21)
    bp  = 2 * (kv01 + kv31 + kv02 + kv32)
    bm  = 2 * (kv01 - kv31 + kv02 - kv32)
    bps = 2 * (kv01 + kv31 - kv02 - kv32)
    bms = 2 * (kv01 - kv31 - kv02 + kv32)
    psikv%psi(1)%a(1) = am * psi%a(3) + bms * psi%a(4)
    psikv%psi(1)%a(2) = bp * psi%a(3) + (-am) * psi%a(4)
    psikv%psi(1)%a(3) = (-ap) * psi%a(1) + bps * psi%a(2)
    psikv%psi(1)%a(4) = bm * psi%a(1) + ap * psi%a(2)
    psikv%psi(2)%a(1) = bms * psi%a(3) + am * psi%a(4)
    psikv%psi(2)%a(2) = (-am) * psi%a(3) + bp * psi%a(4)
    psikv%psi(2)%a(3) = (-bps) * psi%a(1) + ap * psi%a(2)
    psikv%psi(2)%a(4) = (-ap) * psi%a(1) + (-bm) * psi%a(2)
    psikv%psi(3)%a(1) = imago * (bms * psi%a(3) - am * psi%a(4))
    psikv%psi(3)%a(2) = (-imago) * (am * psi%a(3) + bp * psi%a(4))
    psikv%psi(3)%a(3) = (-imago) * (bps * psi%a(1) + ap * psi%a(2))
    psikv%psi(3)%a(4) = imago * ((-ap) * psi%a(1) + bm * psi%a(2))
    psikv%psi(4)%a(1) = am * psi%a(3) + (-bms) * psi%a(4)
    psikv%psi(4)%a(2) = bp * psi%a(3) + am * psi%a(4)
    psikv%psi(4)%a(3) = ap * psi%a(1) + bps * psi%a(2)
    psikv%psi(4)%a(4) = (-bm) * psi%a(1) + ap * psi%a(2)
  end function grkkggf
  pure function gr_vf (g, v, psi, k) result (psikv)
    type(vectorspinor) :: psikv
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: k
    complex(kind=default), intent(in) :: g
    type(vector) :: vk
    vk = k
    psikv = g * (grkkggf (v, psi, vk))
  end function gr_vf
  pure function gr_vlrf (gl, gr, v, psi, k) result (psikv)
    type(vectorspinor) :: psikv
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_l, psi_r
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: k
    complex(kind=default), intent(in) :: gl, gr
    type(vector) :: vk
    vk = k
    psi_l%a(1:2) = psi%a(1:2)
    psi_l%a(3:4) = 0
    psi_r%a(1:2) = 0
    psi_r%a(3:4) = psi%a(3:4)
    psikv = gl * grkkggf (v, psi_l, vk) + gr * grkkggf (v, psi_r, vk)
  end function gr_vlrf
  pure function grkgggf (psil, psir, k) result (j)
    type(vector) :: j
    type(vectorspinor), intent(in) :: psil
    type(bispinor), intent(in) :: psir
    type(vector), intent(in) :: k
    type(vectorspinor) :: c_psir0, c_psir1, c_psir2, c_psir3
    complex(kind=default) :: kp, km, k12, k12s, ik2
    kp = k%t + k%x(3)
    km = k%t - k%x(3)
    k12  =  (k%x(1) + (0,1)*k%x(2))
    k12s =  (k%x(1) - (0,1)*k%x(2))
    ik2 = (0,1) * k%x(2)
    !!! New version:
    c_psir0%psi(1)%a(1) = (-k%x(3)) * psir%a(3) + (-k12s) * psir%a(4)
    c_psir0%psi(1)%a(2) = (-k12) * psir%a(3) + k%x(3) * psir%a(4)
    c_psir0%psi(1)%a(3) = (-k%x(3)) * psir%a(1) + (-k12s) * psir%a(2)
    c_psir0%psi(1)%a(4) = (-k12) * psir%a(1) + k%x(3) * psir%a(2)
    c_psir0%psi(2)%a(1) = (-k12s) * psir%a(3) + (-k%x(3)) * psir%a(4)
    c_psir0%psi(2)%a(2) = k%x(3) * psir%a(3) + (-k12) * psir%a(4)
    c_psir0%psi(2)%a(3) = k12s * psir%a(1) + k%x(3) * psir%a(2)
    c_psir0%psi(2)%a(4) = (-k%x(3)) * psir%a(1) + k12 * psir%a(2)
    c_psir0%psi(3)%a(1) = (0,1) * ((-k12s) * psir%a(3) + k%x(3) * psir%a(4))
    c_psir0%psi(3)%a(2) = (0,1) * (k%x(3) * psir%a(3) + k12 * psir%a(4))
    c_psir0%psi(3)%a(3) = (0,1) * (k12s * psir%a(1) + (-k%x(3)) * psir%a(2))
    c_psir0%psi(3)%a(4) = (0,1) * ((-k%x(3)) * psir%a(1) + (-k12) * psir%a(2))
    c_psir0%psi(4)%a(1) = (-k%x(3)) * psir%a(3) + k12s * psir%a(4)
    c_psir0%psi(4)%a(2) = (-k12) * psir%a(3) + (-k%x(3)) * psir%a(4)
    c_psir0%psi(4)%a(3) = k%x(3) * psir%a(1) + (-k12s) * psir%a(2)
    c_psir0%psi(4)%a(4) = k12 * psir%a(1) + k%x(3) * psir%a(2)
    !!!
    c_psir1%psi(1)%a(1) = (-ik2) * psir%a(3) + (-km) * psir%a(4)
    c_psir1%psi(1)%a(2) = (-kp) * psir%a(3) + ik2 * psir%a(4)
    c_psir1%psi(1)%a(3) = ik2 * psir%a(1) + (-kp) * psir%a(2)
    c_psir1%psi(1)%a(4) = (-km) * psir%a(1) + (-ik2) * psir%a(2)
    c_psir1%psi(2)%a(1) = (-km) * psir%a(3) + (-ik2) * psir%a(4)
    c_psir1%psi(2)%a(2) = ik2 * psir%a(3) + (-kp) * psir%a(4)
    c_psir1%psi(2)%a(3) = kp * psir%a(1) + (-ik2) * psir%a(2)
    c_psir1%psi(2)%a(4) = ik2 * psir%a(1) + km * psir%a(2)
    c_psir1%psi(3)%a(1) = ((0,-1) * km) * psir%a(3) + (-k%x(2)) * psir%a(4)
    c_psir1%psi(3)%a(2) = (-k%x(2)) * psir%a(3) + ((0,1) * kp) * psir%a(4)
    c_psir1%psi(3)%a(3) = ((0,1) * kp) * psir%a(1) + (-k%x(2)) * psir%a(2)
    c_psir1%psi(3)%a(4) = (-k%x(2)) * psir%a(1) + ((0,-1) * km) * psir%a(2)
    c_psir1%psi(4)%a(1) = (-ik2) * psir%a(3) + km * psir%a(4)
    c_psir1%psi(4)%a(2) = (-kp) * psir%a(3) + (-ik2) * psir%a(4)
    c_psir1%psi(4)%a(3) = (-ik2) *  psir%a(1) + (-kp) * psir%a(2)
    c_psir1%psi(4)%a(4) = km * psir%a(1) + (-ik2) * psir%a(2)
    !!!
    c_psir2%psi(1)%a(1) = (0,1) * (k%x(1) * psir%a(3) + km * psir%a(4))
    c_psir2%psi(1)%a(2) = (0,-1) * (kp * psir%a(3) + k%x(1) * psir%a(4))
    c_psir2%psi(1)%a(3) = (0,1) * ((-k%x(1)) * psir%a(1) + kp * psir%a(2))
    c_psir2%psi(1)%a(4) = (0,1) * ((-km) * psir%a(1) + k%x(1) * psir%a(2))
    c_psir2%psi(2)%a(1) = (0,1) * (km * psir%a(3) + k%x(1) * psir%a(4))
    c_psir2%psi(2)%a(2) = (0,-1) * (k%x(1) * psir%a(3) + kp * psir%a(4))
    c_psir2%psi(2)%a(3) = (0,-1) * (kp * psir%a(1) + (-k%x(1)) * psir%a(2))
    c_psir2%psi(2)%a(4) = (0,-1) * (k%x(1) * psir%a(1) + (-km) * psir%a(2))
    c_psir2%psi(3)%a(1) = (-km) * psir%a(3) + k%x(1) * psir%a(4)
    c_psir2%psi(3)%a(2) = k%x(1) * psir%a(3) + (-kp) * psir%a(4)
    c_psir2%psi(3)%a(3) = kp * psir%a(1) + k%x(1) * psir%a(2)
    c_psir2%psi(3)%a(4) = k%x(1) * psir%a(1) + km * psir%a(2)
    c_psir2%psi(4)%a(1) = (0,1) * (k%x(1) * psir%a(3) + (-km) * psir%a(4))
    c_psir2%psi(4)%a(2) = (0,1) * ((-kp) * psir%a(3) + k%x(1) * psir%a(4))
    c_psir2%psi(4)%a(3) = (0,1) * (k%x(1) * psir%a(1) + kp * psir%a(2))
    c_psir2%psi(4)%a(4) = (0,1) * (km * psir%a(1) + k%x(1) * psir%a(2))
    !!!
    c_psir3%psi(1)%a(1) = (-k%t) * psir%a(3) - k12s * psir%a(4)
    c_psir3%psi(1)%a(2) = k12 * psir%a(3) + k%t * psir%a(4)
    c_psir3%psi(1)%a(3) = (-k%t) * psir%a(1) + k12s * psir%a(2)
    c_psir3%psi(1)%a(4) = (-k12) * psir%a(1) + k%t * psir%a(2)
    c_psir3%psi(2)%a(1) = (-k12s) * psir%a(3) + (-k%t) * psir%a(4)
    c_psir3%psi(2)%a(2) = k%t * psir%a(3) + k12 * psir%a(4)
    c_psir3%psi(2)%a(3) = (-k12s) * psir%a(1) + k%t * psir%a(2)
    c_psir3%psi(2)%a(4) = (-k%t) * psir%a(1) + k12 * psir%a(2)
    c_psir3%psi(3)%a(1) = (0,-1) * (k12s * psir%a(3) + (-k%t) * psir%a(4))
    c_psir3%psi(3)%a(2) = (0,1) * (k%t * psir%a(3) + (-k12) * psir%a(4))
    c_psir3%psi(3)%a(3) = (0,-1) * (k12s * psir%a(1) + k%t * psir%a(2))
    c_psir3%psi(3)%a(4) = (0,-1) * (k%t * psir%a(1) + k12 * psir%a(2))
    c_psir3%psi(4)%a(1) = (-k%t) * psir%a(3) + k12s * psir%a(4)
    c_psir3%psi(4)%a(2) = k12 * psir%a(3) + (-k%t) * psir%a(4)
    c_psir3%psi(4)%a(3) = k%t * psir%a(1) + k12s * psir%a(2)
    c_psir3%psi(4)%a(4) = k12 * psir%a(1) + k%t * psir%a(2)
    j%t    =   2 * (psil * c_psir0)
    j%x(1) =   2 * (psil * c_psir1)
    j%x(2) =   2 * (psil * c_psir2)
    j%x(3) =   2 * (psil * c_psir3)
  end function grkgggf
  pure function v_grf (g, psil, psir, k) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: g
    type(vectorspinor), intent(in) :: psil
    type(bispinor), intent(in) :: psir
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    j = g * grkgggf (psil, psir, vk)
  end function v_grf
  pure function vlr_grf (gl, gr, psil, psir, k) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(vectorspinor), intent(in) :: psil
    type(bispinor), intent(in) :: psir
    type(bispinor) :: psir_l, psir_r
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    psir_l%a(1:2) = psir%a(1:2)
    psir_l%a(3:4) = 0
    psir_r%a(1:2) = 0
    psir_r%a(3:4) = psir%a(3:4)
    j = gl * grkgggf (psil, psir_l, vk) + gr * grkgggf (psil, psir_r, vk)
  end function vlr_grf
  pure function fggkggr (psil, psir, k) result (j)
    type(vector) :: j
    type(vectorspinor), intent(in) :: psir
    type(bispinor), intent(in) :: psil
    type(vector), intent(in) :: k
    type(bispinor) :: c_psir0, c_psir1, c_psir2, c_psir3
    complex(kind=default) :: kp, km, k12, k12s, ik1, ik2
    kp = k%t + k%x(3)
    km = k%t - k%x(3)
    k12  =  k%x(1) + (0,1)*k%x(2)
    k12s =  k%x(1) - (0,1)*k%x(2)
    ik1 = (0,1) * k%x(1)
    ik2 = (0,1) * k%x(2)
    c_psir0%a(1) = k%x(3) * (psir%psi(1)%a(4) + psir%psi(4)%a(4) &
                 + psir%psi(2)%a(3) + (0,1) * psir%psi(3)%a(3)) &
                 - k12 * (psir%psi(1)%a(3) + psir%psi(4)%a(3)) &
                 + k12s * (psir%psi(2)%a(4) + (0,1) * psir%psi(3)%a(4))
    c_psir0%a(2) = k%x(3) * (psir%psi(1)%a(3) - psir%psi(4)%a(3) + &
                   psir%psi(2)%a(4) - (0,1) * psir%psi(3)%a(4)) + &
                   k12s * (psir%psi(1)%a(4) - psir%psi(4)%a(4)) - &
                   k12 * (psir%psi(2)%a(3) - (0,1) * psir%psi(3)%a(3))
    c_psir0%a(3) = k%x(3) * (-psir%psi(1)%a(2) + psir%psi(4)%a(2) + &
                   psir%psi(2)%a(1) + (0,1) * psir%psi(3)%a(1)) + &
                   k12 * (psir%psi(1)%a(1) - psir%psi(4)%a(1)) + &
                   k12s * (psir%psi(2)%a(2) + (0,1) * psir%psi(3)%a(2))
    c_psir0%a(4) = k%x(3) * (-psir%psi(1)%a(1) - psir%psi(4)%a(1) + &
                   psir%psi(2)%a(2) - (0,1) * psir%psi(3)%a(2)) -  &
                   k12s * (psir%psi(1)%a(2) + psir%psi(4)%a(2)) - &
                   k12 * (psir%psi(2)%a(1) - (0,1) * psir%psi(3)%a(1))
    !!!
    c_psir1%a(1) = ik2 * (-psir%psi(1)%a(4) - psir%psi(4)%a(4) - &
                   psir%psi(2)%a(3) - (0,1) * psir%psi(3)%a(3)) - &
                   km * (psir%psi(1)%a(3) + psir%psi(4)%a(3)) + &
                   kp * (psir%psi(2)%a(4) + (0,1) * psir%psi(3)%a(4))
    c_psir1%a(2) = ik2 * (-psir%psi(1)%a(3) - psir%psi(2)%a(4) + &
                   psir%psi(4)%a(3) + (0,1) * psir%psi(3)%a(4)) + &
                   kp * (psir%psi(1)%a(4) - psir%psi(4)%a(4)) - &
                   km * (psir%psi(2)%a(3) - (0,1) * psir%psi(3)%a(3))
    c_psir1%a(3) = ik2 * (-psir%psi(1)%a(2) + psir%psi(2)%a(1) + &
                   psir%psi(4)%a(2) + (0,1) * psir%psi(3)%a(1)) + &
                   kp * (psir%psi(1)%a(1) - psir%psi(4)%a(1)) + &
                   km * (psir%psi(2)%a(2) + (0,1) * psir%psi(3)%a(2))
    c_psir1%a(4) = ik2 * (-psir%psi(1)%a(1) + psir%psi(2)%a(2) - &
                   psir%psi(4)%a(1) - (0,1) * psir%psi(3)%a(2)) - &
                   km * (psir%psi(1)%a(2) + psir%psi(4)%a(2)) - &
                   kp * (psir%psi(2)%a(1) - (0,1) * psir%psi(3)%a(1))
    !!!
    c_psir2%a(1) = ik1 * (psir%psi(2)%a(3) + psir%psi(1)%a(4) &
                   + psir%psi(4)%a(4) + (0,1) * psir%psi(3)%a(3)) - &
                   ((0,1)*km) * (psir%psi(1)%a(3) + psir%psi(4)%a(3)) &
                   + kp * (psir%psi(3)%a(4) - (0,1) * psir%psi(2)%a(4))
    c_psir2%a(2) = ik1 * (psir%psi(1)%a(3) + psir%psi(2)%a(4) - &
                   psir%psi(4)%a(3) - (0,1) * psir%psi(3)%a(4)) - &
                   ((0,1)*kp) * (psir%psi(1)%a(4) - psir%psi(4)%a(4)) &
                   - km * (psir%psi(3)%a(3) + (0,1) * psir%psi(2)%a(3))
    c_psir2%a(3) = ik1 * (psir%psi(1)%a(2) - psir%psi(2)%a(1) - &
                   psir%psi(4)%a(2) - (0,1) * psir%psi(3)%a(1)) + &
                   ((0,1)*kp) * (psir%psi(1)%a(1) - psir%psi(4)%a(1)) &
                   + km * (psir%psi(3)%a(2) - (0,1) * psir%psi(2)%a(2))
    c_psir2%a(4) = ik1 * (psir%psi(1)%a(1) - psir%psi(2)%a(2) + &
                   psir%psi(4)%a(1) + (0,1) * psir%psi(3)%a(2)) + &
                   ((0,1)*km) * (psir%psi(1)%a(2) + psir%psi(4)%a(2)) - &
                   kp * (psir%psi(3)%a(1) + (0,1) * psir%psi(2)%a(1))
    !!!
    c_psir3%a(1) = k%t * (psir%psi(1)%a(4) + psir%psi(4)%a(4) + &
                   psir%psi(2)%a(3) + (0,1) * psir%psi(3)%a(3)) - &
                   k12 * (psir%psi(1)%a(3) + psir%psi(4)%a(3)) - &
                   k12s * (psir%psi(2)%a(4) + (0,1) * psir%psi(3)%a(4))
    c_psir3%a(2) = k%t * (psir%psi(1)%a(3) - psir%psi(4)%a(3) + &
                   psir%psi(2)%a(4) - (0,1) * psir%psi(3)%a(4)) - &
                   k12s * (psir%psi(1)%a(4) - psir%psi(4)%a(4)) - &
                   k12 * (psir%psi(2)%a(3) - (0,1) * psir%psi(3)%a(3))
    c_psir3%a(3) = k%t * (-psir%psi(1)%a(2) + psir%psi(2)%a(1) + &
                   psir%psi(4)%a(2) + (0,1) * psir%psi(3)%a(1)) - &
                   k12 * (psir%psi(1)%a(1) - psir%psi(4)%a(1)) + &
                   k12s * (psir%psi(2)%a(2) + (0,1) * psir%psi(3)%a(2))
    c_psir3%a(4) = k%t * (-psir%psi(1)%a(1) + psir%psi(2)%a(2) - &
                   psir%psi(4)%a(1) - (0,1) * psir%psi(3)%a(2)) - &
                   k12s * (psir%psi(1)%a(2) + psir%psi(4)%a(2)) + &
                   k12 * (psir%psi(2)%a(1) - (0,1) * psir%psi(3)%a(1))
    !!! Because we explicitly multiplied the charge conjugation matrix
    !!! we have to omit it from the spinor product and take the
    !!! ordinary product!
    j%t    =   2 * dot_product (conjg (psil%a), c_psir0%a)
    j%x(1) =   2 * dot_product (conjg (psil%a), c_psir1%a)
    j%x(2) =   2 * dot_product (conjg (psil%a), c_psir2%a)
    j%x(3) =   2 * dot_product (conjg (psil%a), c_psir3%a)
  end function fggkggr
  pure function v_fgr (g, psil, psir, k) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: g
    type(vectorspinor), intent(in) :: psir
    type(bispinor), intent(in) :: psil
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    j = g * fggkggr (psil, psir, vk)
  end function v_fgr
  pure function vlr_fgr (gl, gr, psil, psir, k) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(vectorspinor), intent(in) :: psir
    type(bispinor), intent(in) :: psil
    type(bispinor) :: psil_l
    type(bispinor) :: psil_r
    type(momentum), intent(in) :: k
    type(vector) :: vk
    vk = k
    psil_l%a(1:2) = psil%a(1:2)
    psil_l%a(3:4) = 0
    psil_r%a(1:2) = 0
    psil_r%a(3:4) = psil%a(3:4)
    j = gl * fggkggr (psil_l, psir, vk) + gr * fggkggr (psil_r, psir, vk)
  end function vlr_fgr
  pure function f_s2gr (g, phi1, phi2, psi) result (phipsi)
    type(bispinor) :: phipsi
    type(vectorspinor), intent(in) :: psi
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: phi1, phi2
    phipsi = phi2 * f_potgr (g, phi1, psi)
  end function f_s2gr
  pure function f_svgr (g, phi, v, grav) result (phigrav)
    type(bispinor) :: phigrav
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: g, phi
    phigrav = (g * phi) * fgvg5gr (grav, v)
  end function f_svgr
  pure function f_slvgr (gl, phi, v, grav) result (phigrav)
    type(bispinor) :: phigrav, phidum
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: gl, phi
    phidum = (gl * phi) * fgvg5gr (grav, v)
    phigrav%a(1:2) = phidum%a(1:2)
    phigrav%a(3:4) = 0
  end function f_slvgr
  pure function f_srvgr (gr, phi, v, grav) result (phigrav)
    type(bispinor) :: phigrav, phidum
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: gr, phi
    phidum = (gr * phi) * fgvg5gr (grav, v)
    phigrav%a(1:2) = 0
    phigrav%a(3:4) = phidum%a(3:4)
  end function f_srvgr
  pure function f_slrvgr (gl, gr, phi, v, grav) result (phigrav)
    type(bispinor) :: phigrav
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: gl, gr, phi
    phigrav = f_slvgr (gl, phi, v, grav) + f_srvgr (gr, phi, v, grav)
  end function f_slrvgr
  pure function f_pvgr (g, phi, v, grav) result (phigrav)
    type(bispinor) :: phigrav
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: g, phi
    phigrav = (g * phi) * fgvgr (grav, v)
  end function f_pvgr
  pure function f_v2gr (g, v1, v2, grav) result (psi)
    type(bispinor) :: psi
    complex(kind=default), intent(in) :: g
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v1, v2
    psi = g * fggvvgr (v2, grav, v1)
  end function f_v2gr
  pure function f_v2lrgr (gl, gr, v1, v2, grav) result (psi)
    type(bispinor) :: psi
    complex(kind=default), intent(in) :: gl, gr
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v1, v2
    psi = fggvvgr (v2, grav, v1)
    psi%a(1:2) = gl * psi%a(1:2)
    psi%a(3:4) = gr * psi%a(3:4)
  end function f_v2lrgr
  pure function gr_s2f (g, phi1, phi2, psi) result (phipsi)
    type(vectorspinor) :: phipsi
    type(bispinor), intent(in) :: psi
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: phi1, phi2
    phipsi = phi2 * gr_potf (g, phi1, psi)
  end function gr_s2f
  pure function gr_svf (g, phi, v, psi) result (phipsi)
    type(vectorspinor) :: phipsi
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: g, phi
    phipsi = (g * phi) * grkggf (psi, v)
  end function gr_svf
  pure function gr_slvf (gl, phi, v, psi) result (phipsi)
    type(vectorspinor) :: phipsi
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_l
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: gl, phi
    psi_l%a(1:2) = psi%a(1:2)
    psi_l%a(3:4) = 0
    phipsi = (gl * phi) * grkggf (psi_l, v)
  end function gr_slvf
  pure function gr_srvf (gr, phi, v, psi) result (phipsi)
    type(vectorspinor) :: phipsi
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_r
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: gr, phi
    psi_r%a(1:2) = 0
    psi_r%a(3:4) = psi%a(3:4)
    phipsi = (gr * phi) * grkggf (psi_r, v)
  end function gr_srvf
  pure function gr_slrvf (gl, gr, phi, v, psi) result (phipsi)
    type(vectorspinor) :: phipsi
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: gl, gr, phi
    phipsi = gr_slvf (gl, phi, v, psi) + gr_srvf (gr, phi, v, psi)
  end function gr_slrvf
  pure function gr_pvf (g, phi, v, psi) result (phipsi)
    type(vectorspinor) :: phipsi
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    complex(kind=default), intent(in) :: g, phi
    phipsi = (g * phi) * grkgf (psi, v)
  end function gr_pvf
  pure function gr_v2f (g, v1, v2, psi) result (vvpsi)
    type(vectorspinor) :: vvpsi
    complex(kind=default), intent(in) :: g
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v1, v2
    vvpsi = g * grkkggf (v2, psi, v1)
  end function gr_v2f
  pure function gr_v2lrf (gl, gr, v1, v2, psi) result (vvpsi)
    type(vectorspinor) :: vvpsi
    complex(kind=default), intent(in) :: gl, gr
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_l, psi_r
    type(vector), intent(in) :: v1, v2
    psi_l%a(1:2) = psi%a(1:2)
    psi_l%a(3:4) = 0
    psi_r%a(1:2) = 0
    psi_r%a(3:4) = psi%a(3:4)
    vvpsi = gl * grkkggf (v2, psi_l, v1) + gr * grkkggf (v2, psi_r, v1)
  end function gr_v2lrf
  pure function s2_grf (g, gravbar, phi, psi) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g, phi
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    j = phi * pot_grf (g, gravbar, psi)
  end function s2_grf
  pure function s2_fgr (g, psibar, phi, grav) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g, phi
    type(bispinor), intent(in) :: psibar
    type(vectorspinor), intent(in) :: grav
    j = phi * pot_fgr (g, psibar, grav)
  end function s2_fgr
  pure function sv1_grf (g, gravbar, v, psi) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    j = g * grg5vgf (gravbar, psi, v)
  end function sv1_grf
  pure function slv1_grf (gl, gravbar, v, psi) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gl
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_l
    type(vector), intent(in) :: v
    psi_l%a(1:2) = psi%a(1:2)
    psi_l%a(3:4) = 0
    j = gl * grg5vgf (gravbar, psi_l, v)
  end function slv1_grf
  pure function srv1_grf (gr, gravbar, v, psi) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gr
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_r
    type(vector), intent(in) :: v
    psi_r%a(1:2) = 0
    psi_r%a(3:4) = psi%a(3:4)
    j = gr * grg5vgf (gravbar, psi_r, v)
  end function srv1_grf
  pure function slrv1_grf (gl, gr, gravbar, v, psi) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_l, psi_r
    type(vector), intent(in) :: v
    psi_l%a(1:2) = psi%a(1:2)
    psi_l%a(3:4) = 0
    psi_r%a(1:2) = 0
    psi_r%a(3:4) = psi%a(3:4)
    j = gl * grg5vgf (gravbar, psi_l, v) + gr * grg5vgf (gravbar, psi_r, v)
  end function slrv1_grf
  pure function sv2_grf (g, gravbar, phi, psi) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: g, phi
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(vectorspinor) :: g0_psi, g1_psi, g2_psi, g3_psi
    g0_psi%psi(1)%a(1:2) = - psi%a(1:2)
    g0_psi%psi(1)%a(3:4) = psi%a(3:4)
    g0_psi%psi(2)%a(1) = psi%a(2)
    g0_psi%psi(2)%a(2) = psi%a(1)
    g0_psi%psi(2)%a(3) = psi%a(4)
    g0_psi%psi(2)%a(4) = psi%a(3)
    g0_psi%psi(3)%a(1) = (0,-1) * psi%a(2)
    g0_psi%psi(3)%a(2) = (0,1) * psi%a(1)
    g0_psi%psi(3)%a(3) = (0,-1) * psi%a(4)
    g0_psi%psi(3)%a(4) = (0,1) * psi%a(3)
    g0_psi%psi(4)%a(1) = psi%a(1)
    g0_psi%psi(4)%a(2) = - psi%a(2)
    g0_psi%psi(4)%a(3) = psi%a(3)
    g0_psi%psi(4)%a(4) = - psi%a(4)
    g1_psi%psi(1)%a(1:4) = - g0_psi%psi(2)%a(1:4)
    g1_psi%psi(2)%a(1:4) = - g0_psi%psi(1)%a(1:4)
    g1_psi%psi(3)%a(1) = (0,1) * psi%a(1)
    g1_psi%psi(3)%a(2) = (0,-1) * psi%a(2)
    g1_psi%psi(3)%a(3) = (0,-1) * psi%a(3)
    g1_psi%psi(3)%a(4) = (0,1) * psi%a(4)
    g1_psi%psi(4)%a(1) = - psi%a(2)
    g1_psi%psi(4)%a(2) = psi%a(1)
    g1_psi%psi(4)%a(3) = psi%a(4)
    g1_psi%psi(4)%a(4) = - psi%a(3)
    g2_psi%psi(1)%a(1:4) = - g0_psi%psi(3)%a(1:4)
    g2_psi%psi(2)%a(1:4) = - g1_psi%psi(3)%a(1:4)
    g2_psi%psi(3)%a(1:4) = - g0_psi%psi(1)%a(1:4)
    g2_psi%psi(4)%a(1) = (0,1) * psi%a(2)
    g2_psi%psi(4)%a(2) = (0,1) * psi%a(1)
    g2_psi%psi(4)%a(3) = (0,-1) * psi%a(4)
    g2_psi%psi(4)%a(4) = (0,-1) * psi%a(3)
    g3_psi%psi(1)%a(1:4) = - g0_psi%psi(4)%a(1:4)
    g3_psi%psi(2)%a(1:4) = - g1_psi%psi(4)%a(1:4)
    g3_psi%psi(3)%a(1:4) = - g2_psi%psi(4)%a(1:4)
    g3_psi%psi(4)%a(1:4) = - g0_psi%psi(1)%a(1:4)
    j%t    =   (g * phi) * (gravbar * g0_psi)
    j%x(1) =   (g * phi) * (gravbar * g1_psi)
    j%x(2) =   (g * phi) * (gravbar * g2_psi)
    j%x(3) =   (g * phi) * (gravbar * g3_psi)
  end function sv2_grf
  pure function slv2_grf (gl, gravbar, phi, psi) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl, phi
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_l
    psi_l%a(1:2) = psi%a(1:2)
    psi_l%a(3:4) = 0
    j = sv2_grf (gl, gravbar, phi, psi_l)
  end function slv2_grf
  pure function srv2_grf (gr, gravbar, phi, psi) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gr, phi
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_r
    psi_r%a(1:2) = 0
    psi_r%a(3:4) = psi%a(3:4)
    j = sv2_grf (gr, gravbar, phi, psi_r)
  end function srv2_grf
  pure function slrv2_grf (gl, gr, gravbar, phi, psi) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl, gr, phi
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_l, psi_r
    psi_l%a(1:2) = psi%a(1:2)
    psi_l%a(3:4) = 0
    psi_r%a(1:2) = 0
    psi_r%a(3:4) = psi%a(3:4)
    j = sv2_grf (gl, gravbar, phi, psi_l) + sv2_grf (gr, gravbar, phi, psi_r)
  end function slrv2_grf
  pure function sv1_fgr (g, psibar, v, grav) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g
    type(bispinor), intent(in) :: psibar
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v
    j = g * fg5gkgr (psibar, grav, v)
  end function sv1_fgr
  pure function slv1_fgr (gl, psibar, v, grav) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gl
    type(bispinor), intent(in) :: psibar
    type(bispinor) :: psibar_l
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v
    psibar_l%a(1:2) = psibar%a(1:2)
    psibar_l%a(3:4) = 0
    j = gl * fg5gkgr (psibar_l, grav, v)
  end function slv1_fgr
  pure function srv1_fgr (gr, psibar, v, grav) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gr
    type(bispinor), intent(in) :: psibar
    type(bispinor) :: psibar_r
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v
    psibar_r%a(1:2) = 0
    psibar_r%a(3:4) = psibar%a(3:4)
    j = gr * fg5gkgr (psibar_r, grav, v)
  end function srv1_fgr
  pure function slrv1_fgr (gl, gr, psibar, v, grav) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(bispinor), intent(in) :: psibar
    type(bispinor) :: psibar_l, psibar_r
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v
    psibar_l%a(1:2) = psibar%a(1:2)
    psibar_l%a(3:4) = 0
    psibar_r%a(1:2) = 0
    psibar_r%a(3:4) = psibar%a(3:4)
    j = gl * fg5gkgr (psibar_l, grav, v)  + gr * fg5gkgr (psibar_r, grav, v)
  end function slrv1_fgr
  pure function sv2_fgr (g, psibar, phi, grav) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: g, phi
    type(bispinor), intent(in) :: psibar
    type(vectorspinor), intent(in) :: grav
    type(bispinor) :: g0_grav, g1_grav, g2_grav, g3_grav
    g0_grav%a(1) = -grav%psi(1)%a(1) +  grav%psi(2)%a(2) - &
                  (0,1) * grav%psi(3)%a(2) + grav%psi(4)%a(1)
    g0_grav%a(2) = -grav%psi(1)%a(2) + grav%psi(2)%a(1) + &
                  (0,1) * grav%psi(3)%a(1) - grav%psi(4)%a(2)
    g0_grav%a(3) = grav%psi(1)%a(3) + grav%psi(2)%a(4) - &
                  (0,1) * grav%psi(3)%a(4) + grav%psi(4)%a(3)
    g0_grav%a(4) = grav%psi(1)%a(4) + grav%psi(2)%a(3) + &
                  (0,1) * grav%psi(3)%a(3) - grav%psi(4)%a(4)
    !!!
    g1_grav%a(1) = grav%psi(1)%a(2) - grav%psi(2)%a(1) + &
                  (0,1) * grav%psi(3)%a(1) - grav%psi(4)%a(2)
    g1_grav%a(2) = grav%psi(1)%a(1) - grav%psi(2)%a(2) - &
                  (0,1) * grav%psi(3)%a(2) + grav%psi(4)%a(1)
    g1_grav%a(3) = grav%psi(1)%a(4) + grav%psi(2)%a(3) - &
                  (0,1) * grav%psi(3)%a(3) + grav%psi(4)%a(4)
    g1_grav%a(4) = grav%psi(1)%a(3) + grav%psi(2)%a(4) + &
                  (0,1) * grav%psi(3)%a(4) - grav%psi(4)%a(3)
    !!!
    g2_grav%a(1) = (0,1) * (-grav%psi(1)%a(2) - grav%psi(2)%a(1) + &
                  grav%psi(4)%a(2)) - grav%psi(3)%a(1)
    g2_grav%a(2) = (0,1) * (grav%psi(1)%a(1) + grav%psi(2)%a(2) + &
                  grav%psi(4)%a(1)) - grav%psi(3)%a(2)
    g2_grav%a(3) = (0,1) * (-grav%psi(1)%a(4) + grav%psi(2)%a(3) - &
                  grav%psi(4)%a(4)) + grav%psi(3)%a(3)
    g2_grav%a(4) = (0,1) * (grav%psi(1)%a(3) - grav%psi(2)%a(4) - &
                  grav%psi(4)%a(3)) + grav%psi(3)%a(4)
    !!!
    g3_grav%a(1) = -grav%psi(1)%a(2) + grav%psi(2)%a(2) - &
                  (0,1) * grav%psi(3)%a(2) - grav%psi(4)%a(1)
    g3_grav%a(2) = grav%psi(1)%a(1) - grav%psi(2)%a(1) - &
                  (0,1) * grav%psi(3)%a(1) - grav%psi(4)%a(2)
    g3_grav%a(3) = -grav%psi(1)%a(2) - grav%psi(2)%a(4) + &
                  (0,1) * grav%psi(3)%a(4) + grav%psi(4)%a(3)
    g3_grav%a(4) = -grav%psi(1)%a(4) + grav%psi(2)%a(3) + &
                  (0,1) * grav%psi(3)%a(3) + grav%psi(4)%a(4)
    j%t    =   (g * phi) * (psibar * g0_grav)
    j%x(1) =   (g * phi) * (psibar * g1_grav)
    j%x(2) =   (g * phi) * (psibar * g2_grav)
    j%x(3) =   (g * phi) * (psibar * g3_grav)
  end function sv2_fgr
  pure function slv2_fgr (gl, psibar, phi, grav) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl, phi
    type(bispinor), intent(in) :: psibar
    type(bispinor) :: psibar_l
    type(vectorspinor), intent(in) :: grav
    psibar_l%a(1:2) = psibar%a(1:2)
    psibar_l%a(3:4) = 0
    j = sv2_fgr (gl, psibar_l, phi, grav)
  end function slv2_fgr
  pure function srv2_fgr (gr, psibar, phi, grav) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gr, phi
    type(bispinor), intent(in) :: psibar
    type(bispinor) :: psibar_r
    type(vectorspinor), intent(in) :: grav
    psibar_r%a(1:2) = 0
    psibar_r%a(3:4) = psibar%a(3:4)
    j = sv2_fgr (gr, psibar_r, phi, grav)
  end function srv2_fgr
  pure function slrv2_fgr (gl, gr, psibar, phi, grav) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl, gr, phi
    type(bispinor), intent(in) :: psibar
    type(bispinor) :: psibar_l, psibar_r
    type(vectorspinor), intent(in) :: grav
    psibar_l%a(1:2) = psibar%a(1:2)
    psibar_l%a(3:4) = 0
    psibar_r%a(1:2) = 0
    psibar_r%a(3:4) = psibar%a(3:4)
    j = sv2_fgr (gl, psibar_l, phi, grav) + sv2_fgr (gr, psibar_r, phi, grav)
  end function slrv2_fgr
  pure function pv1_grf (g, gravbar, v, psi) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    j = g * grvgf (gravbar, psi, v)
  end function pv1_grf
  pure function pv2_grf (g, gravbar, phi, psi) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: g, phi
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(bispinor) :: g5_psi
    g5_psi%a(1:2) = - psi%a(1:2)
    g5_psi%a(3:4) = psi%a(3:4)
    j = sv2_grf (g, gravbar, phi, g5_psi)
  end function pv2_grf
  pure function pv1_fgr (g, psibar, v, grav) result (j)
    complex(kind=default) :: j
    complex(kind=default), intent(in) :: g
    type(bispinor), intent(in) :: psibar
    type(vectorspinor), intent(in) :: grav
    type(vector), intent(in) :: v
    j = g * fgkgr (psibar, grav, v)
  end function pv1_fgr
  pure function pv2_fgr (g, psibar, phi, grav) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: g, phi
    type(vectorspinor), intent(in) :: grav
    type(bispinor), intent(in) :: psibar
    type(bispinor) :: psibar_g5
    psibar_g5%a(1:2) = - psibar%a(1:2)
    psibar_g5%a(3:4) = psibar%a(3:4)
    j = sv2_fgr (g, psibar_g5, phi, grav)
  end function pv2_fgr
  pure function v2_grf (g, gravbar, v, psi) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: g
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(vector), intent(in) :: v
    j = -g * grkgggf (gravbar, psi, v)
  end function v2_grf
  pure function v2lr_grf (gl, gr, gravbar, v, psi) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(vectorspinor), intent(in) :: gravbar
    type(bispinor), intent(in) :: psi
    type(bispinor) :: psi_l, psi_r
    type(vector), intent(in) :: v
    psi_l%a(1:2) = psi%a(1:2)
    psi_l%a(3:4) = 0
    psi_r%a(1:2) = 0
    psi_r%a(3:4) = psi%a(3:4)
    j = -(gl * grkgggf (gravbar, psi_l, v) + gr * grkgggf (gravbar, psi_r, v))
  end function v2lr_grf
  pure function v2_fgr (g, psibar, v, grav) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: g
    type(vectorspinor), intent(in) :: grav
    type(bispinor), intent(in) :: psibar
    type(vector), intent(in) :: v
    j = -g * fggkggr (psibar, grav, v)
  end function v2_fgr
  pure function v2lr_fgr (gl, gr, psibar, v, grav) result (j)
    type(vector) :: j
    complex(kind=default), intent(in) :: gl, gr
    type(vectorspinor), intent(in) :: grav
    type(bispinor), intent(in) :: psibar
    type(bispinor) :: psibar_l, psibar_r
    type(vector), intent(in) :: v
    psibar_l%a(1:2) = psibar%a(1:2)
    psibar_l%a(3:4) = 0
    psibar_r%a(1:2) = 0
    psibar_r%a(3:4) = psibar%a(3:4)
    j = -(gl * fggkggr (psibar_l, grav, v) + gr * fggkggr (psibar_r, grav, v))
  end function v2lr_fgr
  pure function pr_psi (p, m, w, cms, psi) result (ppsi)
    type(bispinor) :: ppsi
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    type(bispinor), intent(in) :: psi
    logical, intent(in) :: cms
    type(vector) :: vp
    complex(kind=default), parameter :: one = (1, 0)
    complex(kind=default) :: num_mass
    vp = p
    if (cms) then
       num_mass = sqrt(cmplx(m**2, -m*w, kind=default))
    else
       num_mass = cmplx (m, 0, kind=default)
    end if  
    ppsi = (1 / cmplx (p*p - m**2, m*w, kind=default)) &
         * (- f_vf (one, vp, psi) + num_mass * psi)
  end function pr_psi
  pure function pj_psi (p, m, w, psi) result (ppsi)
    type(bispinor) :: ppsi
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    type(bispinor), intent(in) :: psi
    type(vector) :: vp
    complex(kind=default), parameter :: one = (1, 0)
    vp = p
    ppsi = (0, -1) * sqrt (PI / m / w) * (- f_vf (one, vp, psi) + m * psi)
  end function pj_psi
  pure function pg_psi (p, m, w, psi) result (ppsi)
    type(bispinor) :: ppsi
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    type(bispinor), intent(in) :: psi
    type(vector) :: vp
    complex(kind=default), parameter :: one = (1, 0)
    vp = p
    ppsi = gauss (p*p, m, w) * (- f_vf (one, vp, psi) + m * psi)
  end function pg_psi
  pure function pr_grav (p, m, w, grav) result (propgrav)
    type(vectorspinor) :: propgrav
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    type(vectorspinor), intent(in) :: grav
    type(vector) :: vp
    type(bispinor) :: pgrav, ggrav, ggrav1, ggrav2, ppgrav
    type(vectorspinor) :: etagrav_dum, etagrav, pppgrav, &
                          gg_grav_dum, gg_grav
    complex(kind=default), parameter :: one = (1, 0)
    real(kind=default) :: minv
    integer :: i
    vp = p
    minv = 1/m
    pgrav = p%t    * grav%psi(1) - p%x(1) * grav%psi(2) - &
            p%x(2) * grav%psi(3) - p%x(3) * grav%psi(4)
    ggrav%a(1) = grav%psi(1)%a(3) - grav%psi(2)%a(4) + (0,1) * &
                 grav%psi(3)%a(4) - grav%psi(4)%a(3)
    ggrav%a(2) = grav%psi(1)%a(4) - grav%psi(2)%a(3) - (0,1) * &
                 grav%psi(3)%a(3) + grav%psi(4)%a(4)
    ggrav%a(3) = grav%psi(1)%a(1) + grav%psi(2)%a(2) - (0,1) * &
                 grav%psi(3)%a(2) + grav%psi(4)%a(1)
    ggrav%a(4) = grav%psi(1)%a(2) + grav%psi(2)%a(1) + (0,1) * &
                 grav%psi(3)%a(1) - grav%psi(4)%a(2)
    ggrav1 = ggrav - minv * pgrav
    ggrav2 = f_vf (one, vp, ggrav1) + m * ggrav - pgrav
    ppgrav = (-minv**2) * f_vf (one, vp, pgrav) + minv * pgrav
    do i = 1, 4
    etagrav_dum%psi(i) = f_vf (one, vp, grav%psi(i))
    end do
    etagrav = etagrav_dum - m * grav
    pppgrav%psi(1) = p%t    * ppgrav
    pppgrav%psi(2) = p%x(1) * ppgrav
    pppgrav%psi(3) = p%x(2) * ppgrav
    pppgrav%psi(4) = p%x(3) * ppgrav
    gg_grav_dum%psi(1) = p%t    * ggrav2
    gg_grav_dum%psi(2) = p%x(1) * ggrav2
    gg_grav_dum%psi(3) = p%x(2) * ggrav2
    gg_grav_dum%psi(4) = p%x(3) * ggrav2
    gg_grav = gr_potf (one, one, ggrav2) - minv * gg_grav_dum
    propgrav = (1 / cmplx (p*p - m**2, m*w, kind=default)) * &
         (etagrav + pppgrav + (1/3.0_default) * gg_grav)
  end function pr_grav
end module omega_bispinor_couplings
