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
module omega_couplings
  use kinds
  use constants
  use omega_vectors
  use omega_tensors
  implicit none
  private
  public :: g_gg
  public :: x_gg, g_gx
  public :: v_ss, s_vs
  public :: s_vv_t, v_sv_t
  public :: tphi_vv, tphi_vv_cf, v_tphiv, v_tphiv_cf
  public :: tkv_vv, lkv_vv, tv_kvv, lv_kvv, kg_kgkg
  public :: t5kv_vv, l5kv_vv, t5v_kvv, l5v_kvv, kg5_kgkg, kg_kg5kg
  public :: dv_vv, v_dvv, dv_vv_cf, v_dvv_cf
  public :: dv_phi2,phi_dvphi, dv_phi2_cf, phi_dvphi_cf
  public :: phi_vv, v_phiv, phi_u_vv, v_u_phiv
  public :: s_vv_6D, v_sv_6D, s_vv_6DP, v_sv_6DP, a_hz_D, h_az_D, z_ah_D, &
       a_hz_DP, h_az_DP, z_ah_DP, h_hh_6
  public :: g_gg_13, g_gg_23, g_gg_6, kg_kgkg_i
  public ::a_ww_DP, w_aw_DP, a_ww_DW
  public :: w_wz_DPW, z_ww_DPW, w_wz_DW, z_ww_DW, w_wz_D, z_ww_D
  public :: hhhh_p2, a_hww_DPB, h_aww_DPB, w_ahw_DPB, a_hww_DPW, h_aww_DPW, &
       w_ahw_DPW, a_hww_DW, h_aww_DW, w3_ahw_DW, w4_ahw_DW
  public ::a_aww_DW, w_aaw_DW, a_aww_W, w_aaw_W
  public :: h_hww_D, w_hhw_D, h_hww_DP, w_hhw_DP, h_hvv_PB, v_hhv_PB
  public :: a_hhz_D, h_ahz_D, z_ahh_D, a_hhz_DP, h_ahz_DP, z_ahh_DP, &
       a_hhz_PB, h_ahz_PB, z_ahh_PB
  public :: h_wwz_DW, w_hwz_DW, z_hww_DW, h_wwz_DPB, w_hwz_DPB, z_hww_DPB 
  public :: h_wwz_DDPW, w_hwz_DDPW, z_hww_DDPW, h_wwz_DPW, w_hwz_DPW, z_hww_DPW
  public :: phi_dim5s2
  public :: tphi_ss, tphi_ss_cf, s_tphis, s_tphis_cf
  public :: phi_phi2v_1, v_phi2v_1, phi_phi2v_2, v_phi2v_2
  public :: s_dim8s3
  public :: phi_phi2v_m_0, v_phi2v_m_0, phi_phi2v_m_1, v_phi2v_m_1, phi_phi2v_m_7, v_phi2v_m_7
  public :: g_dim8g3_t_0, g_dim8g3_t_1, g_dim8g3_t_2
  public :: g_dim8g3_m_0, g_dim8g3_m_1, g_dim8g3_m_7
  public :: s_gravs, v_gravv, grav_ss, grav_vv
  public :: t2_vv, v_t2v, t2_vv_cf, v_t2v_cf, & 
         t2_vv_1, v_t2v_1, t2_vv_t, v_t2v_t, &
         t2_phi2, phi_t2phi, t2_phi2_cf, phi_t2phi_cf
  public :: t2_vv_d5_1, v_t2v_d5_1
  public :: t2_vv_d5_2, v_t2v_d5_2
  public :: t2_vv_d7, v_t2v_d7
  public :: wd_tl
  public :: wd_run
  public :: gauss
  public :: pr_phi, pr_unitarity, pr_feynman, pr_gauge, pr_rxi
  public :: pr_vector_pure
  public :: pj_phi, pj_unitarity
  public :: pg_phi, pg_unitarity
  public :: pr_tensor, pr_tensor_pure
  integer, parameter, public :: omega_couplings_2010_01_A = 0
contains
  pure function g_gg (g, a1, k1, a2, k2) result (a)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: a
    a = (0, -1) * g * ((k1 - k2) * (a1 * a2) &
                        + ((2*k2 + k1) * a1) * a2 - a1 * ((2*k1 + k2) * a2))
  end function g_gg
  pure function x_gg (g, a1, a2) result (x)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(tensor2odd) :: x
    x = g * (a1 .wedge. a2)
  end function x_gg
  pure function g_gx (g, a1, x) result (a)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: a1
    type(tensor2odd), intent(in) :: x
    type(vector) :: a
    a = g * (a1 * x)
  end function g_gx
  pure function v_ss (g, phi1, k1, phi2, k2) result (v)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    v = (k1 - k2) * (g * phi1 * phi2)
  end function v_ss
  pure function s_vs (g, v1, k1, phi2, k2) result (phi)
    complex(kind=default), intent(in) :: g, phi2
    type(vector), intent(in) :: v1
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: phi
    phi = g * ((k1 + 2*k2) * v1) * phi2
  end function s_vs
  pure function s_vv_t (g, v1, k1, v2, k2) result (phi)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: phi
    phi = g * ((v1*v2) * (k1*k2) - (v1*k2) * (v2*k1))
  end function s_vv_t
  pure function v_sv_t (g, phi, kphi,v, kv) result (vout)
    complex(kind=default), intent(in) :: g, phi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: kv, kphi
    type(momentum) :: kout
    type(vector)  :: vout
    kout = - (kv + kphi)
    vout = g * phi * ((kout*kv) * v - (v * kout) * kv)
  end function v_sv_t
  pure function tphi_vv (g, v1, k1, v2, k2) result (phi)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: phi
    type(momentum) :: k
    k = - (k1 + k2)
    phi = 2 * g * (v1*k) * (v2*k)
  end function tphi_vv
  pure function tphi_vv_cf (g, v1, k1, v2, k2) result (phi)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: phi
    type(momentum) :: k
    k = - (k1 + k2)
    phi = - g/2 * (v1*v2) * (k*k)
  end function tphi_vv_cf
  pure function v_tphiv (g, phi, kphi,v, kv) result (vout)
    complex(kind=default), intent(in) :: g, phi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: kv, kphi
    type(momentum) :: kout
    type(vector)  :: vout
    kout = - (kv + kphi)
    vout = 2 * g * phi * ((v * kout) * kout)
  end function v_tphiv
  pure function v_tphiv_cf (g, phi, kphi,v, kv) result (vout)
    complex(kind=default), intent(in) :: g, phi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: kv, kphi
    type(momentum) :: kout
    type(vector)  :: vout
    kout = - (kv + kphi)
    vout = -g/2 * phi * (kout*kout) * v
  end function v_tphiv_cf
  pure function tkv_vv (g, v1, k1, v2, k2) result (v)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    v = (k1 - k2) * ((0, 1) * g * (v1*v2))
  end function tkv_vv
  pure function t5kv_vv (g, v1, k1, v2, k2) result (v)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    type(vector) :: k
    k = k1 - k2
    v = (0, 1) * g * pseudo_vector (k, v1, v2)
  end function t5kv_vv
  pure function lkv_vv (g, v1, k1, v2, k2) result (v)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    v = (k1 + k2) * ((0, 1) * g * (v1*v2))
  end function lkv_vv
  pure function l5kv_vv (g, v1, k1, v2, k2) result (v)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    type(vector) :: k
    k = k1 + k2
    v = (0, 1) * g * pseudo_vector (k, v1, v2)
  end function l5kv_vv
  pure function tv_kvv (g, v1, k1, v2, k2) result (v)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    v = v2 * ((0, 1) * g * ((2*k2 + k1)*v1))
  end function tv_kvv
  pure function t5v_kvv (g, v1, k1, v2, k2) result (v)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    type(vector) :: k
    k = k1 + 2*k2
    v = (0, 1) * g * pseudo_vector (k, v1, v2)
  end function t5v_kvv
  pure function lv_kvv (g, v1, k1, v2) result (v)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1
    type(vector) :: v
    v = v2 * ((0, -1) * g * (k1*v1))
  end function lv_kvv
  pure function l5v_kvv (g, v1, k1, v2) result (v)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1
    type(vector) :: v
    type(vector) :: k
    k = k1
    v = (0, -1) * g * pseudo_vector (k, v1, v2)
  end function l5v_kvv
  pure function kg_kgkg (g, a1, k1, a2, k2) result (a)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: a
    real(kind=default) :: k1k1, k2k2, k1k2, kk1, kk2
    complex(kind=default) :: a1a2, k2a1, ka1, k1a2, ka2
    k1k1 = k1 * k1
    k1k2 = k1 * k2
    k2k2 = k2 * k2
    kk1 = k1k1 + k1k2
    kk2 = k1k2 + k2k2
    k2a1 = k2 * a1
    ka1 = k2a1 + k1 * a1
    k1a2 = k1 * a2
    ka2 = k1a2 + k2 * a2
    a1a2 = a1 * a2
    a = (0, -1) * g * (   (kk2  * k1a2 - k1k2 * ka2 ) * a1 &
                        + (k1k2 * ka1  - kk1  * k2a1) * a2 &
                        + (ka2  * k2a1 - kk2  * a1a2) * k1 &
                        + (kk1  * a1a2 - ka1  * k1a2) * k2 )
  end function kg_kgkg
  pure function kg5_kgkg (g, a1, k1, a2, k2) result (a)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: a
    type(vector) :: kv, k1v, k2v
    kv = - k1 - k2
    k1v = k1
    k2v = k2
    a = (0, -2) * g * (   (k2*A1) * pseudo_vector (kv, k1v, a2 ) &
                        + (k1*A2) * pseudo_vector (kv, A1 , k2v) &
                        - (A1*A2) * pseudo_vector (kv, k1v, k2v) &
                        - (k1*k2) * pseudo_vector (kv, a1 , a2 ) )
  end function kg5_kgkg
  pure function kg_kg5kg (g, a1, k1, a2, k2) result (a)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: a
    type(vector) :: kv, k1v, k2v
    kv = - k1 - k2
    k1v = k1
    k2v = k2
    a = (0, -1) * g * (   (kv*k2v) * pseudo_vector (a2 , k1v, a1) &
                        - (kv*a2 ) * pseudo_vector (k2v, k1v, a1) &
                        -  k2v * pseudo_scalar (kv, a2,  k1v, a1) &
                        +  a2  * pseudo_scalar (kv, k2v, k1v, a1) )
  end function kg_kg5kg
  pure function dv_vv (g, v1, k1, v2, k2) result (v)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    type(vector) :: k
    k = -(k1 + k2)
    v = g * ((k * v1) * v2 + (k * v2) * v1)
  end function dv_vv
  pure function dv_vv_cf (g, v1, k1, v2, k2) result (v)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    type(vector) :: k
    k = -(k1 + k2)
    v = - g/2 * (v1 * v2) * k
  end function dv_vv_cf
  pure function v_dvv (g, v, k, v2) result (v1)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v, v2
    type(momentum), intent(in) :: k
    type(vector) :: v1
    v1 = g * ((v * v2) * k + (k * v2) * v)
  end function v_dvv
  pure function v_dvv_cf (g, v, k, v2) result (v1)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) ::  v, v2
    type(momentum), intent(in) :: k
    type(vector) :: v1
    v1 = - g/2 * (v * k) * v2
  end function v_dvv_cf
  pure function dv_phi2 (g, phi1, k1, phi2, k2) result (v)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    v = g * phi1 * phi2 * ( &
        (k1 * k2 + k2 * k2 ) * k1 + &
        (k1 * k2 + k1 * k1 ) * k2 )
  end function dv_phi2
  pure function dv_phi2_cf (g, phi1, k1, phi2, k2) result (v)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: v
    v = - g/2 * phi1 * phi2 * (k1 * k2) * (k1 + k2)
  end function dv_phi2_cf
  pure function phi_dvphi (g, v, k, phi2, k2) result (phi1)
    complex(kind=default), intent(in) :: g, phi2
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: k, k2
    complex(kind=default) :: phi1
    type(momentum) :: k1
    k1 = - (k + k2)
    phi1 = g * phi2 * ( &
        (k1 * k2 + k2 * k2 ) * ( k1 * V ) + &
        (k1 * k2 + k1 * k1 ) * ( k2 * V ) )
  end function phi_dvphi
  pure function phi_dvphi_cf (g, v, k, phi2, k2) result (phi1)
    complex(kind=default), intent(in) :: g, phi2
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: k, k2
    complex(kind=default) :: phi1
    type(momentum) :: k1
    k1 = -(k + k2)
    phi1 = - g/2 * phi2 * (k1 * k2)  * ((k1 + k2) * v)
  end function phi_dvphi_cf
  pure function phi_vv (g, k1, k2, v1, v2) result (phi)
    complex(kind=default), intent(in) :: g
    type(momentum), intent(in) :: k1, k2
    type(vector), intent(in) :: v1, v2
    complex(kind=default) :: phi
    phi = g * pseudo_scalar (k1, v1, k2, v2)
  end function phi_vv
  pure function v_phiv (g, phi, k1, k2, v) result (w)
    complex(kind=default), intent(in) :: g, phi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: k1, k2
    type(vector) :: w
    w = g * phi * pseudo_vector (k1, k2, v)
  end function v_phiv
  pure function phi_u_vv (g, k1, k2, v1, v2) result (phi)
    complex(kind=default), intent(in) :: g
    type(momentum), intent(in) :: k1, k2
    type(vector), intent(in) :: v1, v2
    complex(kind=default) :: phi
    phi = g * ((k1*v2)*((-(k1+k2))*v1) + &
               (k2*v1)*((-(k1+k2))*v2) + &
               (((k1+k2)*(k1+k2)) * (v1*v2)))
  end function phi_u_vv
  pure function v_u_phiv (g, phi, k1, k2, v) result (w)
    complex(kind=default), intent(in) :: g, phi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: k1, k2
    type(vector) :: w
    w = g * phi * ((k1*v)*k2 + &
         ((-(k1+k2))*v)*k1 + &
         ((k1*k1)*v))
  end function v_u_phiv
  pure function s_vv_6D (g, v1, k1, v2, k2) result (phi)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: phi
    phi =  g * (-(k1 * v1) * (k1 * v2) - (k2 * v1) * (k2 * v2) &
         + ((k1 * k1) + (k2 * k2)) * (v1 * v2))
  end function s_vv_6D
  pure function v_sv_6D (g, phi, kphi, v, kv) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: phi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: kphi, kv
    type(vector) :: vout
    vout = g * ( - phi * (kv * v) * kv - phi * ((kphi + kv) * v) * (kphi + kv) &
         + phi * (kv * kv) * v + phi * ((kphi + kv)*(kphi + kv)) * v)
  end function v_sv_6D
  pure function s_vv_6DP (g, v1, k1, v2, k2) result (phi)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: phi
    phi = g * ( (-(k1+k2)*v1) * (k1*v2) - ((k1+k2)*v2) * (k2*v1) + &
         ((k1+k2)*(k1+k2))*(v1*v2) )
  end function s_vv_6DP
  pure function v_sv_6DP (g, phi, kphi, v, kv) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: phi
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: kphi, kv
    type(vector) :: vout
    vout = g * phi * ((-(kphi + kv)*v) * kphi + (kphi * v) * kv + &
         (kphi*kphi) * v )
  end function v_sv_6DP 
  pure function a_hz_D (g, h1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * h1 * (((k1 + k2) * v2) * (k1 + k2) + &
         ((k1 + k2) * (k1 + k2)) * v2)
  end function a_hz_D
  pure function h_az_D (g, v1, k1, v2, k2) result (hout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: hout
    hout = g * ((k1 * v1) * (k1 * v2) + (k1 * k1) * (v1 * v2))
  end function h_az_D
  pure function z_ah_D (g, v1, k1, h2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h2
    type(vector), intent(in) :: v1
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * h2 * ((k1 * v1) * k1 + ((k1 * k1)) *v1)
  end function z_ah_D
  pure function a_hz_DP (g, h1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * ((- h1 * (k1 + k2) * v2) * (k1) &
         + h1 * ((k1 + k2) * (k1)) *v2)
  end function a_hz_DP
  pure function h_az_DP (g, v1, k1, v2, k2) result (hout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: hout
    hout = g * (- (k1 * v2) * ((k1 + k2) * v1) + (k1 * (k1 + k2)) * (v1 * v2))
  end function h_az_DP
  pure function z_ah_DP (g, v1, k1, h2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h2
    type(vector), intent(in) :: v1
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * h2* ((k2 * v1) * k1 - (k1 * k2) * v1)
  end function z_ah_DP
  pure function h_hh_6 (g, h1, k1, h2, k2) result (hout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1, h2
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: hout
    hout =  g * ((k1* k1) + (k2 * k2) + (k1* k2)) * h1 * h2
  end function h_hh_6
  pure function g_gg_23 (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * (v1 * (-2*(k1*v2)) + v2 * (2*k2 * v1) + (k1 - k2) * (v1*v2))
  end function g_gg_23
  pure function g_gg_13 (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * (v1 * (2*(k1 + k2)*v2) - v2 * ((k1 + 2*k2) * v1) + 2*k2 * (v1 * v2))
  end function g_gg_13
  pure function g_gg_6 (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * &
         ( k1 * ((-(k1 + k2) * v2) * (k2 * v1) + ((k1 + k2) * k2) * (v1 * v2)) &
         + k2 * (((k1 + k2) * v1) * (k1 * v2) - ((k1 + k2) * k1) * (v1 * v2)) &
         + v1 * (-((k1 + k2) * k2) * (k1 * v2) + (k1 * k2) * ((k1 + k2) * v2)) &
         + v2 * (((k1 + k2) * k1) * (k2 * v1) - (k1 * k2) * ((k1 + k2) * v1)))
  end function g_gg_6
  pure function kg_kgkg_i (g, a1, k1, a2, k2) result (a)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: a1, a2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: a
    real(kind=default) :: k1k1, k2k2, k1k2, kk1, kk2
    complex(kind=default) :: a1a2, k2a1, ka1, k1a2, ka2
    k1k1 = k1 * k1
    k1k2 = k1 * k2
    k2k2 = k2 * k2
    kk1 = k1k1 + k1k2
    kk2 = k1k2 + k2k2
    k2a1 = k2 * a1
    ka1 = k2a1 + k1 * a1
    k1a2 = k1 * a2
    ka2 = k1a2 + k2 * a2
    a1a2 = a1 * a2
    a = (-1) * g * (   (kk2  * k1a2 - k1k2 * ka2 ) * a1 &
         + (k1k2 * ka1  - kk1  * k2a1) * a2 &
         + (ka2  * k2a1 - kk2  * a1a2) * k1 &
         + (kk1  * a1a2 - ka1  * k1a2) * k2 )
  end function kg_kgkg_i
  pure function a_ww_DP (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * ( - ((k1 + k2) * v2) * v1 + ((k1 + k2) * v1) * v2)
  end function a_ww_DP
  pure function w_aw_DP (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * ((k1 * v2) * v1 - (v1 * v2) * k1)
  end function w_aw_DP
  pure function a_ww_DW (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * (v1 * (- (4*k1 + 2*k2) * v2) &
         + v2 * ( (2*k1 + 4*k2) * v1) &
         + (k1 - k2) * (2*v1*v2))
  end function a_ww_DW
  pure function w_wz_DPW (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * (v1 * (-(k1+k2)*v2 - k1*v2) + v2 * ((k1+k2)*v1) + k1 * (v1*v2))
  end function w_wz_DPW
  pure function z_ww_DPW (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * (k1*(v1*v2) - k2*(v1*v2) - v1*(k1*v2) + v2*(k2*v1))
  end function z_ww_DPW
  pure function w_wz_DW (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * (v2 * (v1 * k2) - k2 * (v1 * v2))
  end function w_wz_DW
  pure function z_ww_DW (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * (v1 * ((-1)*(k1+k2) * v2) + v2 * ((k1+k2) * v1))
  end function z_ww_DW
  pure function w_wz_D (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * (v2 * (k2*v1) - k2 * (v1*v2))
  end function w_wz_D
  pure function z_ww_D (g, v1, k1, v2, k2) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: vout
    vout = g * (v1 * (- (k1 + k2) * v2) + v2 * ((k1 + k2) * v1))
  end function z_ww_D

  pure function hhhh_p2 (g, h1, k1, h2, k2, h3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1, h2, h3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * h1*h2*h3* (k1*k1 + k2*k2 +k3*k3 + k1*k3 + k1*k2 + k2*k3)
  end function hhhh_p2
  pure function a_hww_DPB (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * (v3*((k1+k2+k3)*v2) - v2*((k1+k2+k3)*v3))
  end function a_hww_DPB
  pure function h_aww_DPB (g, v1, k1, v2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * ((k1 * v3) * (v1 * v2) - (k1 * v2) * (v1 * v3))
  end function h_aww_DPB
  pure function w_ahw_DPB (g, v1, k1, h2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h2
    type(vector), intent(in) :: v1, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h2 * (v1 * (k1 * v3) - k1 * (v1 * v3))
  end function w_ahw_DPB
  pure function a_hww_DPW (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * (v3 * ((2*k1+k2+k3)*v2) - v2 * ((2*k1+k2+k3)*v3))
  end function a_hww_DPW
  pure function h_aww_DPW (g, v1, k1, v2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * ((-(2*k1+k2+k3)*v2)*(v1*v3)+((2*k1+k2+k3)*v3)*(v1*v2))
  end function h_aww_DPW
  pure function w_ahw_DPW (g, v1, k1, h2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h2
    type(vector), intent(in) :: v1, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h2 * ((k2 - k1) * (v1 * v3) + v1 * ((k1 - k2) * v3))
  end function w_ahw_DPW
  pure function a_hww_DW (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * ( v2 * (-(3*k1 + 4*k2 + 4*k3) * v3) &
         + v3 * ((3*k1 + 2*k2 + 4*k3) * v2)  &
         + (k2 - k3) *2*(v2 * v3))
  end function a_hww_DW
  pure function h_aww_DW (g, v1, k1, v2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * ((v1*v2) * ((3*k1 - k2 - k3)*v3) &
         + (v1*v3) * ((-3*k1 - k2 + k3)*v2) &
         + (v2*v3) * (2*(k2-k3)*v1))
  end function h_aww_DW
  pure function w3_ahw_DW (g, v1, k1, h2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h2
    type(vector), intent(in) :: v1, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h2 * (v1 * ((4*k1 + k2) * v3) &
         +v3 * (-2*(k1 + k2 + 2*k3) * v1) &
         +(-2*k1 + k2 + 2*k3) * (v1*v3))
  end function w3_ahw_DW
  pure function w4_ahw_DW (g, v1, k1, h2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h2
    type(vector), intent(in) :: v1, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h2 * (v1 * (-(4*k1 + k2 + 2*k3) * v3) &
         + v3 * (2*(k1 + k2 + 2*k3) * v1) &
         +(4*k1 + k2) * (v1*v3))
  end function w4_ahw_DW
  pure function a_aww_DW (g, v1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * (2*v1*(v2*v3) - v2*(v1*v3) - v3*(v1*v2))
  end function a_aww_DW
  pure function w_aaw_DW (g, v1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * (2*v3*(v1*v2) - v2*(v1*v3) - v1*(v2*v3))
  end function w_aaw_DW
  pure function a_aww_W (g, v1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
  !!! Recalculated WK 2018-08-24
    type(momentum) :: k4
    k4 = -(k1+k2+k3)
  !!! negative sign (-g) causes expected gauge cancellation
    vout = (-g) * ( &
         + (k1*v3)*(k3*v2)*v1 - (k3*v2)*(v1*v3)*k1 &
         - (k1*k3)*(v2*v3)*v1 + (k3*v1)*(v2*v3)*k1 &
         - (k1*v3)*(v1*v2)*k3 + (k1*v2)*(v1*v3)*k3 &
         + (k1*k3)*(v1*v2)*v3 - (k3*v1)*(k1*v2)*v3 &
         + (k3*v2)*(k4*v3)*v1 - (k3*v2)*(k4*v1)*v3 &
         - (k3*k4)*(v2*v3)*v1 + (k4*v1)*(v2*v3)*k3 &
         - (k3*v1)*(k4*v3)*v2 + (k3*v1)*(k4*v2)*v3 &
         + (k3*k4)*(v1*v3)*v2 - (k4*v2)*(v1*v3)*k3 &
         + (k1*v2)*(k2*v3)*v1 - (k2*v3)*(v1*v2)*k1 &
         - (k1*k2)*(v2*v3)*v1 + (k2*v1)*(v2*v3)*k1 &
         - (k1*v2)*(v1*v3)*k2 + (k1*v3)*(v1*v2)*k2 &
         + (k1*k2)*(v1*v3)*v2 - (k2*v1)*(k1*v3)*v2 &
         + (k2*v3)*(k4*v2)*v1 - (k2*v3)*(k4*v1)*v2 &
         - (k2*k4)*(v2*v3)*v1 + (k4*v1)*(v2*v3)*k2 &
         - (k2*v1)*(k4*v2)*v3 + (k2*v1)*(k4*v3)*v2 &
         + (k2*k4)*(v1*v2)*v3 - (k4*v3)*(v1*v2)*k2 &
         )
  !!! Original Version
  !   vout = g * (v1*((-(k2+k3)*v2)*(k2*v3) + (-(k2+k3)*v3)*(k3*v2)) &
  !        +v2*((-((k2-k3)*v1)*(k1+k2+k3)*v3) - (k1*v3)*(k2*v1) &
  !        + ((k1+k2+k3)*v1)*(k2*v3)) &
  !        +v3*(((k2-k3)*v1)*((k1+k2+k3)*v2) - (k1*v2)*(k3*v1) &
  !        + ((k1+k2+k3)*v1)*(k3*v2)) &
  !        +(v1*v2)*(((2*k1+k2+k3)*v3)*k2 - (k2*v3)*k1 -(k1*v3)*k3) &
  !        +(v1*v3)*(((2*k1+k2+k3)*v2)*k3 - (k3*v2)*k1 - (k1*v2)*k3) &
  !        +(v2*v3)*((-(k1+k2+k3)*v1)*(k2+k3) + ((k2+k3)*v1)*k1) &
  !        +(-(k1+k2+k3)*k3 +k1*k2)*((v1*v3)*v2 - (v2*v3)*v1) &
  !        +(-(k1+k2+k3)*k2 + k1*k3)*((v1*v2)*v3 - (v2*v3)*v1))
  end function a_aww_W
  pure function w_aaw_W (g, v1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
  !!! Recalculated WK 2018-08-25
    type(momentum) :: k4
    k4 = -(k1+k2+k3)
  !!! negative sign (-g) causes expected gauge cancellation
    vout = (-g) * ( &
         + (k3*v1)*(k1*v2)*v3 - (k1*v2)*(v3*v1)*k3 &
         - (k3*k1)*(v2*v1)*v3 + (k1*v3)*(v2*v1)*k3 &
         - (k3*v1)*(v3*v2)*k1 + (k3*v2)*(v3*v1)*k1 &
         + (k3*k1)*(v3*v2)*v1 - (k1*v3)*(k3*v2)*v1 &
         + (k1*v2)*(k4*v1)*v3 - (k1*v2)*(k4*v3)*v1 &
         - (k1*k4)*(v2*v1)*v3 + (k4*v3)*(v2*v1)*k1 &
         - (k1*v3)*(k4*v1)*v2 + (k1*v3)*(k4*v2)*v1 &
         + (k1*k4)*(v3*v1)*v2 - (k4*v2)*(v3*v1)*k1 &
         + (k3*v2)*(k2*v1)*v3 - (k2*v1)*(v3*v2)*k3 &
         - (k3*k2)*(v2*v1)*v3 + (k2*v3)*(v2*v1)*k3 &
         - (k3*v2)*(v3*v1)*k2 + (k3*v1)*(v3*v2)*k2 &
         + (k3*k2)*(v3*v1)*v2 - (k2*v3)*(k3*v1)*v2 &
         + (k2*v1)*(k4*v2)*v3 - (k2*v1)*(k4*v3)*v2 &
         - (k2*k4)*(v2*v1)*v3 + (k4*v3)*(v2*v1)*k2 &
         - (k2*v3)*(k4*v2)*v1 + (k2*v3)*(k4*v1)*v2 &
         + (k2*k4)*(v3*v2)*v1 - (k4*v1)*(v3*v2)*k2 &
         )
  !!! Original Version
  !   vout = g * (v1*((k1*v3)*(-(k1+k2+2*k3)*v2) + (k2*v3)*((k1+k2+k3)*v2) &
  !        + (k1*v2)*((k1+k2+k3)*v3)) &
  !        + v2*(((k1-k2)*v3)*((k1+k2+k3)*v1) - (k2*v3)*(k3*v1) &
  !        + (k2*v1)*((k1+k2+k3)*v3)) &
  !        + v3*((k1*v2)*(-(k1+k2)*v1) + (k2*v1)*(-(k1+k2)*v2)) &
  !        + (v1*v2)*((k1+k2)*(-(k1+k2+k3)*v3) + k3*((k1+k2)*v3))&
  !        + (v1*v3)*(-k2*(k3*v2) - k3*(k1*v2) + k1*((k1+k2+2*k3)*v2)) &
  !        + (v2*v3)*(-k1*(k3*v1) - k3*(k2*v1) + k2*((k1+k2+2*k3)*v1)) &
  !        + (-k2*(k1+k2+k3) + k1*k3)*(v1*(v2*v3) - v3*(v1*v2)) &
  !        + (-k1*(k1+k2+k3) + k2*k3)*(v2*(v1*v3) - v3*(v1*v2)) )
  end function w_aaw_W
  pure function h_hww_D (g, h1, k1, v2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1 
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * h1 * ((v2*v3)*((k2*k2)+(k3*k3)) - (k2*v2)*(k2*v3) &
         - (k3*v2)*(k3*v3))
  end function h_hww_D
  pure function w_hhw_D (g, h1, k1, h2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1, h2 
    type(vector), intent(in) :: v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * h2 * (v3 * ((k1+k2+k3)*(k1+k2+k3)+(k3*k3)) &
         - (k1+k2+k3) * ((k1+k2+k3)*v3) - k3 * (k3*v3))
  end function w_hhw_D
  pure function h_hww_DP (g, h1, k1, v2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1 
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * h1 * (-((k2+k3)*v2)*(k2*v3) - &
         ((k2+k3)*v3)*(k3*v2)+ (v2*v3)*((k2+k3)*(k2+k3)))
  end function h_hww_DP
  pure function w_hhw_DP (g, h1, k1, h2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1, h2 
    type(vector), intent(in) :: v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * h2 * (k3*((k1+k2)*v3) + (k1+k2)*(-(k1+k2+k3)*v3) &
         + v3*((k1+k2)*(k1+k2)))
  end function w_hhw_DP
  pure function h_hvv_PB (g, h1, k1, v2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1 
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * h1 * ((k2*v3)*(k3*v2) - (k2*k3)*(v2*v3))
  end function h_hvv_PB
  pure function v_hhv_PB (g, h1, k1, h2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1, h2 
    type(vector), intent(in) :: v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * h2 * ((-(k1+k2+k3)*v3)*k3 + ((k1+k2+k3)*k3)*v3)
  end function v_hhv_PB
  pure function a_hhz_D (g, h1, k1, h2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1, h2
    type(vector), intent(in) :: v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * h2 * ((k1+k2+k3) * ((k1+k2+k3)*v3) &
         - v3 * ((k1+k2+k3)*(k1+k2+k3)))
  end function a_hhz_D
  pure function h_ahz_D (g, v1, k1, h2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h2 
    type(vector), intent(in) :: v1, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * h2 * ((k1*v1)*(k1*v3) - (k1*k1)*(v1*v3))
  end function h_ahz_D
  pure function z_ahh_D (g, v1, k1, h2, k2, h3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1
    complex(kind=default), intent(in) :: h2, h3 
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h2 * h3 * ((k1*v1)*k1 - (k1*k1)*v1)
  end function z_ahh_D
  pure function a_hhz_DP (g, h1, k1, h2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1, h2
    type(vector), intent(in) :: v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * h2 * ((-(k1+k2+k3)*v3)*(k1+k2) + ((k1+k2+k3)*(k1+k2))*v3)
  end function a_hhz_DP
  pure function h_ahz_DP (g, v1, k1, h2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h2 
    type(vector), intent(in) :: v1, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * h2 * ( (k1*v3)*(-(k1+k3)*v1) + (k1*(k1+k3))*(v1*v3) )
  end function h_ahz_DP
  pure function z_ahh_DP (g, v1, k1, h2, k2, h3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1
    complex(kind=default), intent(in) :: h2, h3 
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h2 * h3 * (k1*((k2+k3)*v1) - v1*(k1*(k2+k3)))
  end function z_ahh_DP
  pure function a_hhz_PB (g, h1, k1, h2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1, h2
    type(vector), intent(in) :: v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * h2 * (k3*((k1+k2+k3)*v3) - v3*((k1+k2+k3)*k3))
  end function a_hhz_PB
  pure function h_ahz_PB (g, v1, k1, h2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h2 
    type(vector), intent(in) :: v1, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * h2 * ((-k1*v3)*(k3*v1) + (k1*k3)*(v1*v3))
  end function h_ahz_PB
  pure function z_ahh_PB (g, v1, k1, h2, k2, h3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1
    complex(kind=default), intent(in) :: h2, h3 
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h2 * h3 * (k1*((k1+k2+k3)*v1) - v1*(k1*(k1+k2+k3)))
  end function z_ahh_PB
  pure function h_wwz_DW (g, v1, k1, v2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * (((k1-k2)*v3)*(v1*v2)-((2*k1+k2)*v2)*(v1*v3) + &
         ((k1+2*k2)*v1)*(v2*v3))
  end function h_wwz_DW
  pure function w_hwz_DW (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * ( v2*(-(k1+2*k2+k3)*v3) + v3*((2*k1+k2+2*k3)*v2) - &
         (k1 - k2 + k3)*(v2*v3))
  end function w_hwz_DW
  pure function z_hww_DW (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * ((k2-k3)*(v2*v3) - v2*((2*k2+k3)*v3) + v3*((k2+2*k3)*v2))
  end function z_hww_DW
  pure function h_wwz_DPB (g, v1, k1, v2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * ((k3*v1)*(v2*v3) - (k3*v2)*(v1*v3))
  end function h_wwz_DPB
  pure function w_hwz_DPB (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * (k3*(v2*v3) - v3*(k3*v2))
  end function w_hwz_DPB
  pure function z_hww_DPB (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * (((k1+k2+k3)*v3)*v2 - ((k1+k2+k3)*v2)*v3)
  end function z_hww_DPB
  pure function h_wwz_DDPW (g, v1, k1, v2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * (((k1-k2)*v3)*(v1*v2)-((k1-k3)*v2)*(v1*v3)+((k2-k3)*v1)*(v2*v3))
  end function h_wwz_DDPW
  pure function w_hwz_DDPW (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * ((-(k1+2*k2+k3)*v3)*v2 + ((k1+k2+2*k3)*v2)*v3 + &
         (v2*v3)*(k2-k3))
  end function w_hwz_DDPW
  pure function z_hww_DDPW (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * ((v2*v3)*(k2-k3) - ((k1+2*k2+k3)*v3) *v2 + &
         ((k1+k2+2*k3)*v2)*v3 )
  end function z_hww_DDPW
  pure function h_wwz_DPW (g, v1, k1, v2, k2, v3, k3) result (hout)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    complex(kind=default) :: hout
    hout = g * (((k1-k2)*v3)*(v1*v2) + (-(2*k1+k2+k3)*v2)*(v1*v3) + &
         ((k1+2*k2+k3)*v1)*(v2*v3))
  end function h_wwz_DPW
  pure function w_hwz_DPW (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * ((-(k1+2*k2+k3)*v3)*v2 + ((2*k1+k2+k3)*v2)*v3 + &
         (v2*v3)*(k2-k1))
  end function w_hwz_DPW
  pure function z_hww_DPW (g, h1, k1, v2, k2, v3, k3) result (vout)
    complex(kind=default), intent(in) :: g
    complex(kind=default), intent(in) :: h1
    type(vector), intent(in) :: v2, v3
    type(momentum), intent(in) :: k1, k2, k3
    type(vector) :: vout
    vout = g * h1 * ((v2*v3)*(k2-k3) + ((k1-k2)*v3)*v2 + ((k3-k1)*v2)*v3)
  end function z_hww_DPW

  pure function phi_dim5s2 (g, phi2, k2, phi3, k3) result (phi1)
    complex(kind=default), intent(in) :: g, phi2, phi3
    type(momentum), intent(in) :: k2, k3
    complex(kind=default) :: phi1
    phi1 = g * phi2 * phi3 * (k2 * k3)
  end function phi_dim5s2
  pure function tphi_ss (g, phi1, k1, phi2, k2) result (phi)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: phi
    phi = 2 * g * phi1 * phi2 * &
             ((k1 * k2)+ (k1 * k1)) * &
             ((k1 * k2)+ (k2 * k2))
  end function tphi_ss
  pure function tphi_ss_cf (g, phi1, k1, phi2, k2) result (phi)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2
    complex(kind=default) :: phi
    phi = - g/2 * phi1 * phi2 * &
             (k1 * k2) * &
             ((k1 + k2) * (k1 + k2))
  end function tphi_ss_cf
  pure function s_tphis (g, phi, k, phi2, k2) result (phi1)
    complex(kind=default), intent(in) :: g, phi, phi2
    type(momentum), intent(in) :: k, k2
    complex(kind=default) :: phi1
    type(momentum) :: k1
    k1 = - ( k + k2)
    phi1 = 2 * g * phi * phi2 * &
             ((k1 * k2)+ (k1 * k1)) * &
             ((k1 * k2)+ (k2 * k2))
  end function s_tphis
  pure function s_tphis_cf (g, phi, k, phi2, k2) result (phi1)
    complex(kind=default), intent(in) :: g, phi, phi2
    type(momentum), intent(in) :: k, k2
    complex(kind=default) :: phi1
    type(momentum) :: k1
    k1 = - ( k + k2)
    phi1 = - g/2 * phi * phi2 * &
             (k1 * k2) * &
             ((k1 + k2) * (k1 + k2))
  end function s_tphis_cf
  pure function phi_phi2v_1 (g, phi1, k1, v1, k_v1, v2, k_v2) result (phi2)
    complex(kind=default), intent(in) :: g, phi1
    type(momentum), intent(in) :: k1, k_v1, k_v2
    type(momentum) :: k2
    type(vector), intent(in) :: v1, v2
    complex(kind=default) :: phi2
    k2 = - k1 - k_v1 - k_v2
    phi2 = g * phi1 * &
          ( (k1 * v1) * (k2 * v2) + (k1 * v2) * (k2 * v1) )
  end function phi_phi2v_1
  pure function v_phi2v_1 (g, phi1, k1, phi2, k2, v1) result (v2)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2
    type(vector), intent(in) :: v1
    type(vector) :: v2
    v2 = g * phi1 * phi2 * &
          ( k1  * (k2 * v1) + k2 * (k1 * v1) )
  end function v_phi2v_1
  pure function phi_phi2v_2 (g, phi1, k1, v1,k_v1, v2, k_v2) result (phi2)
    complex(kind=default), intent(in) :: g, phi1
    type(momentum), intent(in) :: k1, k_v1, k_v2
    type(vector), intent(in) :: v1, v2
    type(momentum) :: k2
    complex(kind=default) :: phi2
    k2 = - k1 - k_v1 - k_v2
    phi2 = g * phi1 * (k1 * k2) * (v1 * v2)
  end function phi_phi2v_2
  pure function v_phi2v_2 (g, phi1, k1, phi2, k2, v1) result (v2)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2
    type(vector), intent(in) :: v1
    type(vector) :: v2
    v2 = g * phi1 * phi2 * &
          ( k1  * k2 ) * v1
  end function v_phi2v_2
  pure function s_dim8s3 (g, phi2, k2, phi3, k3, phi4, k4) result (phi1)
    complex(kind=default), intent(in) :: g, phi2, phi3, phi4
    type(momentum), intent(in) :: k2, k3, k4
    type(momentum) :: k1
    complex(kind=default) :: phi1
    k1 = - k2 - k3 - k4
    phi1 = g * ( (k1 * k2) * (k3 * k4) + (k1 * k3) * (k2 * k4) &
            + (k1 * k4) * (k2 * k3) ) * phi2 * phi3 * phi4
  end function s_dim8s3
  pure function phi_phi2v_m_0 (g, phi1, k1, v1, k_v1, v2, k_v2) result (phi2)
    complex(kind=default), intent(in) :: g, phi1
    type(momentum), intent(in) :: k1, k_v1, k_v2
    type(momentum) :: k2
    type(vector), intent(in) :: v1, v2
    complex(kind=default) :: phi2
    k2 = - k1 - k_v1 - k_v2
    phi2 = g * phi1 * &
              ( (v1 * k_v2) * (v2 * k_v1) * (k1 * k2) &
              - (v1 * v2) * (k_v1 * k_v2) * (k1 * k2) )
  end function phi_phi2v_m_0 
  pure function v_phi2v_m_0 (g, phi1, k1, phi2, k2, v1, k_v1) result (v2)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2, k_v1
    type(vector), intent(in) :: v1
    type(momentum) :: k_v2
    type(vector) :: v2
    k_v2 = - k_v1 - k1 - k2
    v2 = g * phi1 * phi2 * &
            ( k_v1 * (v1 *  k_v2) * (k1 * k2) &
             - v1 * (k_v2 * k_v1) * (k1 * k2) )
  end function v_phi2v_m_0
  pure function phi_phi2v_m_1 (g, phi1, k1, v1, k_v1, v2, k_v2) result (phi2)
    complex(kind=default), intent(in) :: g, phi1
    type(momentum), intent(in) :: k1, k_v1, k_v2
    type(momentum) :: k2
    type(vector), intent(in) :: v1, v2
    complex(kind=default) :: phi2
    k2 = - k1 - k_v1 - k_v2
    phi2 = g * phi1 * &
              ( (v1 * v2) * (k1 * k_v2) * (k2 * k_v1) &
              + (v1 * v2) * (k1 * k_v1) * (k2 * k_v2) &
              + (v1 * k2) * (v2 * k1) * (k_v1 * k_v2) &
              + (v1 * k1) * (v2 * k2) * (k_v1 * k_v2) &
              - (v1 * k_v2) * (v2 * k2) * (k1 * k_v1) &
              - (v1 * k2) * (v2 * k_v1) * (k1 * k_v2) &
              - (v1 * k_v2) * (v2 * k1) * (k2 * k_v1) &
              - (v1 * k1) * (v2 * k_v1) * (k2 * k_v2) )
  end function phi_phi2v_m_1
  pure function v_phi2v_m_1 (g, phi1, k1, phi2, k2, v1, k_v1) result (v2)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2, k_v1
    type(vector), intent(in) :: v1
    type(momentum) :: k_v2
    type(vector) :: v2
    k_v2 = - k_v1 - k1 - k2
    v2 = g * phi1 * phi2 * &
            ( k1 * (v1 * k2) * (k_v1 * k_v2) &
            + k2 * (v1 * k1) * (k_v1 * k_v2) &
            + v1 * (k_v1 * k1) * (k_v2 * k2) &
            + v1 * (k_v1 * k2) * (k_v2 * k1) &
            - k1 * (v1 * k_v2) * (k_v1 * k2) &
            - k2 * (v1 * k_v2) * (k_v1 * k1) &
            - k_v1 * (v1 * k1) * (k_v2 * k2) &
            - k_v1 * (v1 * k2) * (k_v2 * k1) )
  end function v_phi2v_m_1
  pure function phi_phi2v_m_7 (g, phi1, k1, v1, k_v1, v2, k_v2) result (phi2)
    complex(kind=default), intent(in) :: g, phi1
    type(momentum), intent(in) :: k1, k_v1, k_v2
    type(momentum) :: k2
    type(vector), intent(in) :: v1, v2
    complex(kind=default) :: phi2
    k2 = - k1 - k_v1 - k_v2
    phi2 = g * phi1 * &
              ( (v1 * k_v2) * (k1 * v2) * (k2 * k_v1) &
              + (v1 * k_v2) * (k1 * k_v1) * (k2 * v2) &
              + (v1 * k1) * (v2 * k_v1) * (k2 * k_v2) &
              + (v1 * k2) * (v2 * k_v1) * (k1 * k_v2) &
              - (v1 * v2) * (k1 * k_v2) * (k2 * k_v1) &
              - (v1 * v2) * (k1 * k_v1) * (k2 * k_v2) &
              - (v1 * k2) * (v2 * k1) * (k_v1 * k_v2) &
              - (v1 * k1) * (v2 * k2) * (k_v1 * k_v2) )
  end function phi_phi2v_m_7
  pure function v_phi2v_m_7 (g, phi1, k1, phi2, k2, v1, k_v1) result (v2)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2, k_v1
    type(vector), intent(in) :: v1
    type(momentum) :: k_v2
    type(vector) :: v2
    k_v2 = - k_v1 - k1 - k2
    v2 = g * phi1 * phi2 * &
            ( k1 * (v1 * k_v2) * (k2 * k_v1) &
            + k2 * (v1 * k_v2) * (k1 * k_v1) &
            + k_v1 * (v1 * k1) * (k2 * k_v2) &
            + k_v1 * (v1 * k2) * (k1 * k_v2) &
            - k1 * (v1 * k2) * (k_v1 * k_v2) &
            - k2 * (v1 * k1) * (k_v1 * k_v2) &
            - v1 * (k1 * k_v2) * (k2 * k_v1) &
            - v1 * (k1 * k_v1) * (k2 * k_v2) )
  end function v_phi2v_m_7
  pure function g_dim8g3_t_0 (g, v2, k2, v3, k3, v4, k4) result (v1)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v2, v3, v4
    type(momentum), intent(in) :: k2, k3, k4
    type(vector) :: v1
    type(momentum) :: k1
    k1 = - k2 - k3 - k4
    v1 = g * (k2 * (k1 * v2) - v2 * (k1 * k2)) & 
           * ((k3 * v4) * (k4 * v3) - (v3 * v4) * (k3 * k4))
  end function g_dim8g3_t_0
  pure function g_dim8g3_t_1 (g, v2, k2, v3, k3, v4, k4) result (v1)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v2, v3, v4
    type(momentum), intent(in) :: k2, k3, k4
    type(vector) :: v1
    type(momentum) :: k1
    k1 = - k2 - k3 - k4
    v1 = g * (v3 * (v2 * k4) * (k1 * k3) * (k2 * v4) &
            + v4 * (v2 * k3) * (k1 * k4) * (k2 * v3) &
            + k3 * (v2 * v4) * (k1 * v3) * (k2 * k4) &
            + k4 * (v2 * v3) * (k1 * v4) * (k2 * k3) &
            - v3 * (v2 * v4) * (k1 * k3) * (k2 * k4) &
            - v4 * (v2 * v3) * (k1 * k4) * (k2 * k3) &
            - k3 * (v2 * k4) * (k1 * v3) * (k2 * v4) &
            - k4 * (v2 * k3) * (k1 * v4) * (k2 * v3))
  end function g_dim8g3_t_1
  pure function g_dim8g3_t_2 (g, v2, k2, v3, k3, v4, k4) result (v1)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v2, v3, v4
    type(momentum), intent(in) :: k2, k3, k4
    type(vector) :: v1
    type(momentum) :: k1
    k1 = - k2 - k3 - k4
    v1 = g * (k2 * (v2 * k3) * (v3 * k4) * (v4 * k1) &
            + k3 * (v2 * k1) * (v3 * k4) * (v4 * k2) &
            + k2 * (v2 * k4) * (v3 * k1) * (v4 * k3) &
            + k4 * (v2 * k1) * (v3 * k2) * (v4 * k3) &
            + k4 * (v2 * k3) * (v3 * v4) * (k1 * k2) &
            + k3 * (v2 * k4) * (v3 * v4) * (k1 * k2) &
            - k3 * (v2 * v4) * (v3 * k4) * (k1 * k2) & 
            - v4 * (v2 * k3) * (v3 * k4) * (k1 * k2) &
            - k4 * (v2 * v3) * (v4 * k3) * (k1 * k2) &
            - v3 * (v2 * k4) * (v4 * k3) * (k1 * k2) &
            - k2 * (v2 * k4) * (v3 * v4) * (k1 * k3) &
            + k2 * (v2 * v4) * (v3 * k4) * (k1 * k3) &
            - v2 * (v3 * k4) * (v4 * k2) * (k1 * k3) &
            - k2 * (v2 * k3) * (v3 * v4) * (k1 * k4) &
            + k2 * (v2 * v3) * (v4 * k3) * (k1 * k4) &
            - v2 * (v3 * k2) * (v4 * k3) * (k1 * k4) &
            - k4 * (v2 * k1) * (v3 * v4) * (k2 * k3) &
            + v4 * (v2 * k1) * (v3 * k4) * (k2 * k3) &
            - v2 * (v3 * k4) * (v4 * k1) * (k2 * k3) &
            + v2 * (v3 * v4) * (k1 * k4) * (k2 * k3) &
            - k3 * (v2 * k1) * (v3 * v4) * (k2 * k4) &
            + v3 * (v2 * k1) * (v4 * k3) * (k2 * k4) &
            - v2 * (v3 * k1) * (v4 * k3) * (k2 * k4) &
            + v2 * (v3 * v4) * (k1 * k3) * (k2 * k4) &
            - k2 * (v2 * v4) * (v3 * k1) * (k3 * k4) &
            - v4 * (v2 * k1) * (v3 * k2) * (k3 * k4) &
            - k2 * (v2 * v3) * (v4 * k1) * (k3 * k4) &
            + v2 * (v3 * k2) * (v4 * k1) * (k3 * k4) &
            - v3 * (v2 * k1) * (v4 * k2) * (k3 * k4) &
            + v2 * (v3 * k1) * (v4 * k2) * (k3 * k4) &
            + v4 * (v2 * v3) * (k1 * k2) * (k3 * k4) &
            + v3 * (v2 * v4) * (k1 * k2) * (k3 * k4))
  end function g_dim8g3_t_2
  pure function g_dim8g3_m_0 (g1, g2, v2, k2, v3, k3, v4, k4) result (v1)
    complex(kind=default), intent(in) :: g1, g2
    type(vector), intent(in) :: v2, v3, v4
    type(momentum), intent(in) :: k2, k3, k4
    type(vector) :: v1
    type(momentum) :: k1
    k1 = - k2 - k3 - k4
    v1 = g1 * (v2 * (v3 * v4) * (k1 * k2)  &
             - k2 * (v2 * k1) * (v3 * v4)) &
       + g2 * (v2 * (v3 * v4) * (k3 * k4)  &
             - v2 * (v3 * k4) * (v4 * k3))
  end function g_dim8g3_m_0
  pure function g_dim8g3_m_1 (g1, g2, v2, k2, v3, k3, v4, k4) result (v1)
    complex(kind=default), intent(in) :: g1, g2
    type(vector), intent(in) :: v2, v3, v4
    type(momentum), intent(in) :: k2, k3, k4
    type(vector) :: v1
    type(momentum) :: k1
    k1 = - k2 - k3 - k4
    v1 = g1 * (k2 * (v2 * v4) * (v3 * k1)  &
             + v4 * (v2 * k1) * (v3 * k2)  &
             + k2 * (v2 * v3) * (v4 * k1)  &
             + v3 * (v2 * k1) * (v4 * k2)  &
             - v2 * (v3 * k2) * (v4 * k1)  &          
             - v2 * (v3 * k1) * (v4 * k2)  &
             - v4 * (v2 * v3) * (k1 * k2)  &
             - v3 * (v2 * v4) * (k1 * k2)) &
       + g2 * (k3 * (v2 * v4) * (v3 * k4)  &   
             - k4 * (v2 * k3) * (v3 * v4)  &
             - k3 * (v2 * k4) * (v3 * v4)  &
             + v4 * (v2 * k3) * (v3 * k4)  &
             + k4 * (v2 * v3) * (v4 * k3)  &
             + v3 * (v2 * k4) * (v4 * k3)  &
             - v4 * (v2 * v3) * (k3 * k4)  &
             - v3 * (v2 * v4) * (k3 * k4))
  end function g_dim8g3_m_1
  pure function g_dim8g3_m_7 (g1, g2, g3, v2, k2, v3, k3, v4, k4) result (v1)
    complex(kind=default), intent(in) :: g1, g2, g3
    type(vector), intent(in) :: v2, v3, v4
    type(momentum), intent(in) :: k2, k3, k4
    type(vector) :: v1
    type(momentum) :: k1
    k1 = - k2 - k3 - k4
    v1 = g1 * (v2 * (v3 * k2) * (v4 * k1)  &
             + v2 * (v3 * k1) * (v4 * k2)  & 
             + v4 * (v2 * v3) * (k1 * k2)  &
             + v3 * (v2 * v4) * (k1 * k2)  &
             - k2 * (v2 * v4) * (v3 * k1)  &
             - v4 * (v2 * k1) * (v3 * k2)  &
             - k2 * (v2 * v3) * (v4 * k1)  &
             - v3 * (v2 * k1) * (v4 * k2)) &
       + g2 * (k3 * (v2 * k1) * (v3 * v4)  & 
             + k4 * (v2 * k1) * (v3 * v4)  &
             + k2 * (v2 * k3) * (v3 * v4)  &
             + k2 * (v2 * k4) * (v3 * v4)  &
             + v4 * (v2 * k4) * (v3 * k1)  &
             + k4 * (v2 * v4) * (v3 * k2)  &
             + v3 * (v2 * k3) * (v4 * k1)  &
             + v2 * (v3 * k4) * (v4 * k1)  &
             + k3 * (v2 * v3) * (v4 * k2)  &
             + v2 * (v3 * k4) * (v4 * k2)  &
             + v2 * (v3 * k1) * (v4 * k3)  &
             + v2 * (v3 * k2) * (v4 * k3)  &
             + v4 * (v2 * v3) * (k1 * k3)  &
             + v3 * (v2 * v4) * (k1 * k4)  &
             + v3 * (v2 * v4) * (k2 * k3)  &
             + v4 * (v2 * v3) * (k2 * k4)  &
             - k4 * (v2 * v4) * (v3 * k1)  &
             - v4 * (v2 * k3) * (v3 * k1)  &
             - k3 * (v2 * v4) * (v3 * k2)  &
             - v4 * (v2 * k4) * (v3 * k2)  &
             - k2 * (v2 * v4) * (v3 * k4)  &
             - v4 * (v2 * k1) * (v3 * k4)  &
             - k3 * (v2 * v3) * (v4 * k1)  &
             - v3 * (v2 * k4) * (v4 * k1)  &
             - k4 * (v2 * v3) * (v4 * k2)  &
             - v3 * (v2 * k3) * (v4 * k2)  &
             - k2 * (v2 * v3) * (v4 * k3)  &
             - v3 * (v2 * k1) * (v4 * k3)  &
             - v2 * (v3 * v4) * (k1 * k3)  &
             - v2 * (v3 * v4) * (k1 * k4)  &
             - v2 * (v3 * v4) * (k2 * k3)  &
             - v2 * (v3 * v4) * (k2 * k4)) &
       + g3 * (k4 * (v2 * k3) * (v3 * v4)  &
             + k3 * (v2 * k4) * (v3 * v4)  &
             + v4 * (v2 * v3) * (k3 * k4)  &
             + v3 * (v2 * v4) * (k3 * k4)  &
             - k3 * (v2 * v4) * (v3 * k4)  &
             - v4 * (v2 * k3) * (v3 * k4)  &
             - k4 * (v2 * v3) * (v4 * k3)  &
             - v3 * (v2 * k4) * (v4 * k3)) 
  end function g_dim8g3_m_7
  pure function s_gravs (g, m, k1, k2, t, s) result (phi)
    complex(kind=default), intent(in) :: g, s
    real(kind=default), intent(in) :: m
    type(momentum), intent(in) :: k1, k2
    type(tensor), intent(in) :: t
    complex(kind=default) :: phi, t_tr
    t_tr = t%t(0,0) - t%t(1,1) - t%t(2,2) - t%t(3,3)
    phi = g * s * (((t*k1)*k2) + ((t*k2)*k1) &
        - g * (m**2 + (k1*k2))*t_tr)/2.0_default
  end function s_gravs
  pure function grav_ss (g, m, k1, k2, s1, s2) result (t)
    complex(kind=default), intent(in) :: g, s1, s2
    real(kind=default), intent(in) :: m
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: t_metric, t
    t_metric%t = 0
    t_metric%t(0,0) = 1.0_default
    t_metric%t(1,1) = - 1.0_default
    t_metric%t(2,2) = - 1.0_default
    t_metric%t(3,3) = - 1.0_default
    t = g*s1*s2/2.0_default * (-(m**2 + (k1*k2)) * t_metric &
      + (k1.tprod.k2) + (k2.tprod.k1))
  end function grav_ss
  pure function v_gravv (g, m, k1, k2, t, v) result (vec)
    complex(kind=default), intent(in) :: g
    real(kind=default), intent(in) :: m
    type(momentum), intent(in) :: k1, k2
    type(vector), intent(in) :: v
    type(tensor), intent(in) :: t
    complex(kind=default) :: t_tr
    real(kind=default) :: xi
    type(vector) :: vec
    xi = 1.0_default
    t_tr = t%t(0,0) - t%t(1,1) - t%t(2,2) - t%t(3,3)
    vec = (-g)/ 2.0_default * (((k1*k2) + m**2) * &
         (t*v + v*t - t_tr * v) + t_tr * (k1*v) * k2 &
         - (k1*v) * ((k2*t) + (t*k2)) &
         - ((k1*(t*v)) + (v*(t*k1))) * k2 &
         + ((k1*(t*k2)) + (k2*(t*k1))) * v)
  !!!       Unitarity gauge: xi -> Infinity
  !!!       + (1.0_default/xi) * (t_tr * ((k1*v)*k2) + &
  !!!       (k2*v)*k2 + (k2*v)*k1 - (k1*(t*v))*k1 + &
  !!!       (k2*v)*(k2*t) - (v*(t*k1))*k1 - (k2*v)*(t*k2)))
  end function v_gravv
  pure function grav_vv (g, m, k1, k2, v1, v2) result (t)
    complex(kind=default), intent(in) :: g
    type(momentum), intent(in) :: k1, k2
    real(kind=default), intent(in) :: m
    real(kind=default) :: xi
    type(vector), intent (in) :: v1, v2
    type(tensor) :: t_metric, t
    xi = 0.00001_default
    t_metric%t = 0
    t_metric%t(0,0) = 1.0_default
    t_metric%t(1,1) = - 1.0_default
    t_metric%t(2,2) = - 1.0_default
    t_metric%t(3,3) = - 1.0_default
    t = (-g)/2.0_default * ( &
         ((k1*k2) + m**2) * ( &
         (v1.tprod.v2) +  (v2.tprod.v1) - (v1*v2) * t_metric) &
         + (v1*k2)*(v2*k1)*t_metric &
         - (k2*v1)*((v2.tprod.k1) + (k1.tprod.v2)) &
         - (k1*v2)*((v1.tprod.k2) + (k2.tprod.v1)) &
         + (v1*v2)*((k1.tprod.k2) + (k2.tprod.k1)))
  !!!       Unitarity gauge: xi -> Infinity
  !!!       + (1.0_default/xi) * ( &
  !!!       ((k1*v1)*(k1*v2) + (k2*v1)*(k2*v2) + (k1*v1)*(k2*v2))* &
  !!!       t_metric) - (k1*v1) * ((k1.tprod.v2) + (v2.tprod.k1)) &
  !!!       - (k2*v2) * ((k2.tprod.v1) + (v1.tprod.k2)))
  end function grav_vv
  pure function t2_vv (g, v1, v2) result (t)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(tensor) :: t
    type(tensor) :: tmp
    tmp = v1.tprod.v2
    t%t = g * (tmp%t + transpose (tmp%t))
  end function t2_vv
  pure function v_t2v (g, t, v) result (tv)
    complex(kind=default), intent(in) :: g
    type(tensor), intent(in) :: t
    type(vector), intent(in) :: v
    type(vector) :: tv
    type(tensor) :: tmp
    tmp%t = t%t + transpose (t%t)
    tv = g * (tmp * v)
  end function v_t2v
  pure function t2_vv_cf (g, v1, v2) result (t)
    complex(kind=default), intent(in) :: g
    complex(kind=default) :: tmp_s
    type(vector), intent(in) :: v1, v2
    type(tensor) :: t_metric, t
    t_metric%t = 0
    t_metric%t(0,0) =   1.0_default
    t_metric%t(1,1) = - 1.0_default
    t_metric%t(2,2) = - 1.0_default
    t_metric%t(3,3) = - 1.0_default
    tmp_s = v1 * v2
    t%t = - (g /2.0_default) * tmp_s * t_metric%t 
  end function t2_vv_cf
  pure function v_t2v_cf (g, t, v) result (tv)
    complex(kind=default), intent(in) :: g
    type(tensor), intent(in) :: t
    type(vector), intent(in) :: v
    type(vector) :: tv, tmp_tv
    tmp_tv =  ( t%t(0,0)-t%t(1,1)-t%t(2,2)-t%t(3,3) ) * v
    tv = - ( g /2.0_default) * tmp_tv
  end function v_t2v_cf
  pure function t2_phi2 (g, phi1, k1, phi2, k2) result (t)
    complex(kind=default), intent(in) :: g, phi1, phi2
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: t
    type(tensor) :: tmp
    tmp = k1.tprod.k2
    t%t = g * (tmp%t + transpose (tmp%t)) * phi1 * phi2
  end function t2_phi2
  pure function phi_t2phi (g, t, kt, phi2, k2) result (phi1)
    complex(kind=default), intent(in) :: g, phi2
    type(tensor), intent(in) :: t
    type(momentum), intent(in) :: kt, k2
    type(momentum) :: k1
    complex(kind=default) :: phi1
    type(tensor) :: tmp
    k1 = -kt - k2 
    tmp%t = t%t + transpose (t%t)
    phi1 = g * ( (tmp * k2) * k1) * phi2
  end function phi_t2phi
  pure function t2_phi2_cf (g, phi1, k1, phi2, k2) result (t)
    complex(kind=default), intent(in) :: g, phi1, phi2
    complex(kind=default) :: tmp_s
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: t_metric, t
    t_metric%t = 0
    t_metric%t(0,0) =   1.0_default
    t_metric%t(1,1) = - 1.0_default
    t_metric%t(2,2) = - 1.0_default
    t_metric%t(3,3) = - 1.0_default
    tmp_s = (k1 * k2) * phi1 * phi2
    t%t = - (g /2.0_default) * tmp_s * t_metric%t 
  end function t2_phi2_cf
  pure function phi_t2phi_cf (g, t, kt, phi2, k2) result (phi1)
    complex(kind=default), intent(in) :: g, phi2
    type(tensor), intent(in) :: t
    type(momentum), intent(in) :: kt, k2
    type(momentum) :: k1
    complex(kind=default) ::  tmp_ts, phi1
    k1 = - kt - k2
    tmp_ts =  ( t%t(0,0)-t%t(1,1)-t%t(2,2)-t%t(3,3) ) 
    phi1 = - ( g /2.0_default) * tmp_ts * (k1 * k2) * phi2
  end function phi_t2phi_cf
  pure function t2_vv_1 (g, v1, v2) result (t)
    complex(kind=default), intent(in) :: g
    complex(kind=default) :: tmp_s
    type(vector), intent(in) :: v1, v2
    type(tensor) :: tmp
    type(tensor) :: t_metric, t
    t_metric%t = 0
    t_metric%t(0,0) =   1.0_default
    t_metric%t(1,1) = - 1.0_default
    t_metric%t(2,2) = - 1.0_default
    t_metric%t(3,3) = - 1.0_default
    tmp = v1.tprod.v2
    tmp_s = v1 * v2
    t%t = g * (tmp%t + transpose (tmp%t) - tmp_s * t_metric%t )
  end function t2_vv_1
  pure function v_t2v_1 (g, t, v) result (tv)
    complex(kind=default), intent(in) :: g
    type(tensor), intent(in) :: t
    type(vector), intent(in) :: v
    type(vector) :: tv, tmp_tv
    type(tensor) :: tmp
    tmp_tv =  ( t%t(0,0)-t%t(1,1)-t%t(2,2)-t%t(3,3) ) * v
    tmp%t = t%t + transpose (t%t)
    tv = g * (tmp * v - tmp_tv)
  end function v_t2v_1
  pure function t2_vv_t (g, v1, k1, v2, k2) result (t)
    complex(kind=default), intent(in) :: g
    complex(kind=default) :: tmp_s
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: tmp, tmp_v1k2, tmp_v2k1, tmp_k1k2, tmp2
    type(tensor) :: t_metric, t
    t_metric%t = 0
    t_metric%t(0,0) =   1.0_default
    t_metric%t(1,1) = - 1.0_default
    t_metric%t(2,2) = - 1.0_default
    t_metric%t(3,3) = - 1.0_default
    tmp = v1.tprod.v2
    tmp_s = v1 * v2
    tmp_v1k2 = (v2 * k1) * (v1.tprod.k2)
    tmp_v2k1 = (v1 * k2) * (v2.tprod.k1)
    tmp_k1k2 = tmp_s * (k1.tprod.k2)
    tmp2%t = tmp_v1k2%t + tmp_v2k1%t - tmp_k1k2%t
    t%t = g * ( (k1*k2) * (tmp%t + transpose (tmp%t) - tmp_s * t_metric%t ) &
         + ((v1 * k2) * (v2 * k1)) * t_metric%t &
         - tmp2%t - transpose(tmp2%t))
  end function t2_vv_t
  pure function v_t2v_t (g, t, kt, v, kv) result (tv)
    complex(kind=default), intent(in) :: g
    type(tensor), intent(in) :: t
    type(vector), intent(in) :: v
    type(momentum), intent(in) :: kt, kv
    type(momentum) :: kout
    type(vector) :: tv, tmp_tv
    type(tensor) :: tmp
    kout = - (kt + kv)
    tmp_tv =  ( t%t(0,0)-t%t(1,1)-t%t(2,2)-t%t(3,3) ) * v
    tmp%t = t%t + transpose (t%t)
    tv = g * ( (tmp * v - tmp_tv) * (kv * kout )&
         + ( t%t(0,0)-t%t(1,1)-t%t(2,2)-t%t(3,3) ) * (kout * v ) * kv &
         - (kout * v) * ( tmp * kv) &
         - (v* (t * kout) + kout * (t * v)) * kv &
         + (kout* (t * kv) + kv * (t * kout)) * v)
  end function v_t2v_t
  pure function t2_vv_d5_1 (g, v1, k1, v2, k2) result (t)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: t
    t = (g * (v1 * v2)) * (k1-k2).tprod.(k1-k2)
  end function t2_vv_d5_1
  pure function v_t2v_d5_1 (g, t1, k1, v2, k2) result (tv)
    complex(kind=default), intent(in) :: g
    type(tensor), intent(in) :: t1
    type(vector), intent(in) :: v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: tv
    tv = (g * ((k1+2*k2).tprod.(k1+2*k2) * t1)) * v2
  end function v_t2v_d5_1
  pure function t2_vv_d5_2 (g, v1, k1, v2, k2) result (t)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: t
    t = (g * (k2 * v1)) * (k2-k1).tprod.v2
    t%t = t%t + transpose (t%t)
  end function t2_vv_d5_2
  pure function v_t2v_d5_2 (g, t1, k1, v2, k2) result (tv)
    complex(kind=default), intent(in) :: g
    type(tensor), intent(in) :: t1
    type(vector), intent(in) :: v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: tv
    type(tensor) :: tmp
    type(momentum) :: k1_k2, k1_2k2
    k1_k2 = k1 + k2
    k1_2k2 = k1_k2 + k2
    tmp%t = t1%t + transpose (t1%t)
    tv = (g * (k1_k2 * v2)) * (k1_2k2 * tmp)
  end function v_t2v_d5_2
  pure function t2_vv_d7 (g, v1, k1, v2, k2) result (t)
    complex(kind=default), intent(in) :: g
    type(vector), intent(in) :: v1, v2
    type(momentum), intent(in) :: k1, k2
    type(tensor) :: t
    t = (g * (k2 * v1) * (k1 * v2)) * (k1-k2).tprod.(k1-k2)
  end function t2_vv_d7
  pure function v_t2v_d7 (g, t1, k1, v2, k2) result (tv)
    complex(kind=default), intent(in) :: g
    type(tensor), intent(in) :: t1
    type(vector), intent(in) :: v2
    type(momentum), intent(in) :: k1, k2
    type(vector) :: tv
    type(vector) :: k1_k2, k1_2k2
    k1_k2 = k1 + k2
    k1_2k2 = k1_k2 + k2
    tv = (- g * (k1_k2 * v2) * (k1_2k2.tprod.k1_2k2 * t1)) * k2
  end function v_t2v_d7
  pure function wd_tl (p, w) result (width)
    real(kind=default) :: width
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: w
    if (p*p > 0) then
       width = w
    else
       width = 0
    end if
  end function wd_tl
  pure function wd_run (p, m, w) result (width)
    real(kind=default) :: width
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m
    real(kind=default), intent(in) :: w
    if (p*p > 0) then
       width = w * (p*p) / m**2
    else
       width = 0
    end if
  end function wd_run
  pure function gauss (x, mu, w) result (gg)
    real(kind=default) :: gg
    real(kind=default), intent(in) :: x, mu, w
    if (w > 0) then
      gg = exp(-(x - mu**2)**2/4.0_default/mu**2/w**2) * &
           sqrt(sqrt(PI/2)) / w / mu
      else
      gg = 1.0_default
    end if
  end function gauss
  pure function pr_phi (p, m, w, phi) result (pphi)
    complex(kind=default) :: pphi
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    complex(kind=default), intent(in) :: phi
    pphi = (1 / cmplx (p*p - m**2, m*w, kind=default)) * phi
  end function pr_phi
  pure function pj_phi (m, w, phi) result (pphi)
    complex(kind=default) :: pphi
    real(kind=default), intent(in) :: m, w
    complex(kind=default), intent(in) :: phi
    pphi = (0, -1) * sqrt (PI / m / w) * phi
  end function pj_phi
  pure function pg_phi (p, m, w, phi) result (pphi)
    complex(kind=default) :: pphi
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    complex(kind=default), intent(in) :: phi
    pphi = ((0, 1) * gauss (p*p, m, w)) * phi
  end function pg_phi
  pure function pr_unitarity (p, m, w, cms, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    type(vector), intent(in) :: e
    logical, intent(in) :: cms
    type(vector) :: pv
    complex(kind=default) :: c_mass2
    pv = p
    if (cms) then
       c_mass2 = cmplx (m**2, -m*w, kind=default)
    else
       c_mass2 = m**2
    end if
    pe = - (1 / cmplx (p*p - m**2, m*w, kind=default)) &
         * (e - (p*e / c_mass2) * pv)
  end function pr_unitarity
  pure function pj_unitarity (p, m, w, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    type(vector), intent(in) :: e
    type(vector) :: pv
    pv = p
    pe = (0, 1) * sqrt (PI / m / w) * (e - (p*e / m**2) * pv)
  end function pj_unitarity
  pure function pg_unitarity (p, m, w, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    type(vector), intent(in) :: e
    type(vector) :: pv
    pv = p
    pe = - gauss (p*p, m, w) &
         * (e - (p*e / m**2) * pv)
  end function pg_unitarity
  pure function pr_feynman (p, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    type(vector), intent(in) :: e
    pe = - (1 / (p*p)) * e
  end function pr_feynman
  pure function pr_gauge (p, xi, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: xi
    type(vector), intent(in) :: e
    real(kind=default) :: p2
    type(vector) :: pv
    p2 = p*p
    pv = p
    pe = - (1 / p2) * (e - ((1 - xi) * (p*e) / p2) * pv)
  end function pr_gauge
  pure function pr_rxi (p, m, w, xi, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w, xi
    type(vector), intent(in) :: e
    real(kind=default) :: p2
    type(vector) :: pv
    p2 = p*p
    pv = p
    pe = - (1 / cmplx (p2 - m**2, m*w, kind=default)) &
         * (e - ((1 - xi) * (p*e) / (p2 - xi * m**2)) * pv)
  end function pr_rxi
  pure function pr_vector_pure (p, m, w, e) result (pe)
    type(vector) :: pe
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    type(vector), intent(in) :: e
    real(kind=default) :: p2
    type(vector) :: pv
    p2 = p*p
    pv = p
    pe = - (1 / cmplx (p2 - m**2, m*w, kind=default)) * e 
  end function pr_vector_pure
  pure function pr_tensor (p, m, w, t) result (pt)
    type(tensor) :: pt
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    type(tensor), intent(in) :: t
    complex(kind=default) :: p_dd_t
    real(kind=default), dimension(0:3,0:3) :: p_uu, p_ud, p_du, p_dd
    integer :: i, j
    p_uu(0,0) = 1 - p%t * p%t / m**2
    p_uu(0,1:3) = - p%t * p%x / m**2
    p_uu(1:3,0) = p_uu(0,1:3)
    do i = 1, 3
       do j = 1, 3
          p_uu(i,j) = - p%x(i) * p%x(j) / m**2
       end do
    end do
    do i = 1, 3
       p_uu(i,i) = - 1 + p_uu(i,i)
    end do
    p_ud(:,0) = p_uu(:,0)
    p_ud(:,1:3) = - p_uu(:,1:3)
    p_du = transpose (p_ud)
    p_dd(:,0) = p_du(:,0)
    p_dd(:,1:3) = - p_du(:,1:3)
    p_dd_t = 0
    do i = 0, 3
       do j = 0, 3
          p_dd_t = p_dd_t + p_dd(i,j) * t%t(i,j)
       end do
    end do
    pt%t = matmul (p_ud, matmul (0.5_default * (t%t + transpose (t%t)), p_du)) &
         - (p_dd_t / 3.0_default) * p_uu
    pt%t = pt%t / cmplx (p*p - m**2, m*w, kind=default)
  end function pr_tensor
  pure function pr_tensor_pure (p, m, w, t) result (pt)
    type(tensor) :: pt
    type(momentum), intent(in) :: p
    real(kind=default), intent(in) :: m, w
    type(tensor), intent(in) :: t
    complex(kind=default) :: p_dd_t
    real(kind=default), dimension(0:3,0:3) :: g_uu
    integer :: i, j
    g_uu(0,0) = 1 
    g_uu(0,1:3) = 0
    g_uu(1:3,0) = g_uu(0,1:3)
    do i = 1, 3
       do j = 1, 3
          g_uu(i,j) = 0
       end do
    end do
    do i = 1, 3
       g_uu(i,i) = - 1
    end do
    p_dd_t = t%t(0,0) - t%t(1,1) - t%t(2,2) - t%t(3,3) 
    pt%t =  0.5_default * ((t%t + transpose (t%t)) &
         - p_dd_t * g_uu )
    pt%t = pt%t / cmplx (p*p - m**2, m*w, kind=default)
  end function pr_tensor_pure
end module omega_couplings
