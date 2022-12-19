! $Id$
!
! Copyright (C) 1999-2022 by
!    Wolfgang Kilian <kilian@physik.uni-siegen.de>
!    Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!    Juergen Reuter <juergen.reuter@desy.de>
!    Christian Speckner <christian.speckner@physik.uni-wuerzburg.de>
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file has been generated from the tangled threeshl sources,
! passed through the C preprocessor and has been stripped of all
! in this process. The documented sources are available from
! http://theorie.physik.uni-wuerzburg.de/~cnspeckn
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module tdefs
use iso_fortran_env
use kinds, only:default
use constants, only:pi
implicit none
public

integer, parameter :: double=default
integer, parameter :: slength=256

logical, parameter :: threeshl_lt_enabled = .true.




contains

function pad_string (string) result(res)
character(len=*), intent(in) :: string
character(len=slength) :: res
        res = string
end function pad_string

subroutine exit (v)
integer, intent(in) :: v
end subroutine exit

end module

module threeshl
use tdefs

use nlowidth

implicit none
save
private

integer :: i
integer, parameter :: errstack_size=5
integer :: errstack_pos=0
character(len=slength), dimension(errstack_size) :: errstack=""
integer, parameter :: err_invalid_parameters=1, err_invalid_root=2, &
 err_sign_error=3, err_range_error=4, err_num_ids=4
character(len=slength), dimension(err_num_ids) :: err_messages= &
 (/ "invalid parameters in function call - invalid model parameters?        ", &
  "square root of negative value - invalid model parameters?              ", &
  "wrong sign found - invalid model parameters?                           ", &
  "result of function exceeds parameter range - invalid model parameters? "/)
logical, public :: threeshl_quit_on_panic=.true., threeshl_error=.false., &
 threeshl_print_errors=.true.
integer, public :: threeshl_errunit = output_unit
integer, parameter, public :: e_bcd=8, nue_bcd=0, mu_bcd=24, numu_bcd=16, tau_bcd=40, &
 nutau_bcd=32, he_bcd=9, hnue_bcd=1, hmu_bcd=25, hnumu_bcd=17, htau_bcd=41, &
 hnutau_bcd=33, u_bcd=4, d_bcd=12, c_bcd=20, s_bcd=28, t_bcd=36, b_bcd=44, hu_bcd=5, &
 hd_bcd=13, hc_bcd=21, hs_bcd=29, ht_bcd=37, hb_bcd=45, w_bcd=2, hw_bcd=3, z_bcd=6, hz_bcd=7, &
 a_bcd=63
integer, parameter, public :: l_chir=100, r_chir=101, l_mode=110, h_mode=111, &
 lh_mode=112, iso_up=120, iso_down=121, lat_0=130, lat_1=131, lat_2=132, &
 gen_0=140, gen_1=141, gen_2=142, ftype_l=150, ftype_q=151, &
 ptype_b=160, ptype_f=161, btype_w=170, btype_z=171, btype_a=172
real(kind=double), public, target :: mass_array(0:63)=0., width_array(0:63)=0.
real(kind=double), public, target :: &
 g_w_lep(l_mode:h_mode, l_mode:h_mode, gen_0:gen_2, &
  l_mode:h_mode, gen_0:gen_2, l_chir:r_chir)= 0., &
 g_w_quark(l_mode:h_mode, l_mode:h_mode, gen_0:gen_2, &
 l_mode:h_mode, gen_0:gen_2, l_chir:r_chir)= 0.
real(kind=double), public, target :: &
 g_z_lep(l_mode:h_mode, l_mode:lh_mode, gen_0:gen_2, iso_up:iso_down, l_chir:r_chir)= 0.,&
 g_z_quark(l_mode:h_mode, l_mode:lh_mode, gen_0:gen_2, iso_up:iso_down, l_chir:r_chir)= 0.
real(kind=double), public, target :: &
 g_wwz(l_mode:lh_mode, l_mode:h_mode)= 0.
real(kind=double), public, target :: &
 g_wwww(0:4)=(/(0.,i=1,5)/), &
 g_wwzz(l_mode:lh_mode, l_mode:lh_mode)= 0., &
 g_wwza(l_mode:lh_mode, l_mode:h_mode)= 0.
real(kind=double), target :: sigma_vev=0., g0=0., g1=0., g2=0., x=0., lambda=0., &
 t=0., eps_l, eps_r(ftype_l:ftype_q, gen_0:gen_2, iso_up:iso_down)= 0., e
real(kind=double) :: &
 wfunct_w(l_mode:h_mode, lat_0:lat_1)= 0., &
 wfunct_z(l_mode:h_mode, lat_0:lat_2)= 0.
real(kind=double) :: &
 wfunct_lep_l(l_mode:h_mode, gen_0:gen_2, iso_up:iso_down, lat_0:lat_1)= 0.,&
 wfunct_lep_r(l_mode:h_mode, gen_0:gen_2, iso_up:iso_down, lat_1:lat_2)= 0.,&
 wfunct_quark_l(l_mode:h_mode, gen_0:gen_2, iso_up:iso_down, lat_0:lat_1)= 0.,&
 wfunct_quark_r(l_mode:h_mode, gen_0:gen_2, iso_up:iso_down, lat_1:lat_2)= 0.
real(kind=double), public :: &
 me_pdg=0._double, mmu_pdg=0.106_double, mtau_pdg=1.78_double, &
 muq_pdg=0._double, mdq_pdg=0._double, mcq_pdg=1.25_double, &
 msq_pdg=0.95_double, mtq_pdg=174._double, mbq_pdg=4.2_double, &
 mw_pdg=80.403_double, mz_pdg=91.188_double, e_pdg=0.313329_double, &
 ww_pdg=2.048_double, wz_pdg=2.443_double, wt_pdg=1.523_double
character(len=slength), public :: particle_names(0:63)= &
 (/("[INVALID PARTICLE BCD OR ARRAY NOT INITIALIZED]",i=0,63)/)
real(kind=double), pointer, public :: threeshl_sigma_vev, threeshl_g0, &
 threeshl_g1, threeshl_g2, threeshl_x, threeshl_lambda, threeshl_t, &
 threeshl_eps_l, threeshl_e

logical, public :: threeshl_use_nlow=.true.




interface threeshl_retrieve_bcd_f
 module procedure retrieve_bcd_f
end interface threeshl_retrieve_bcd_f
public :: threeshl_retrieve_bcd_f
interface threeshl_retrieve_bcd_b
 module procedure retrieve_bcd_b
end interface threeshl_retrieve_bcd_b
public :: threeshl_retrieve_bcd_b
interface threeshl_pdg_init_wgap_bmass
 module procedure pdg_init_wgap_bmass
end interface threeshl_pdg_init_wgap_bmass
public :: threeshl_pdg_init_wgap_bmass
interface threeshl_calculate_widths
 module procedure calculate_widths
end interface threeshl_calculate_widths
public :: threeshl_calculate_widths
interface threeshl_init_ward
 module procedure init_ward
end interface threeshl_init_ward
public :: threeshl_init_ward
interface threeshl_init
 module procedure init
end interface threeshl_init
public :: threeshl_init
interface threeshl_finalize
 module procedure finalize
end interface threeshl_finalize
public :: threeshl_finalize
public :: threeshl_eps_r
interface threeshl_print_particles
 module procedure print_particles
end interface threeshl_print_particles
public :: threeshl_print_particles
interface threeshl_print_gwff
 module procedure print_gwff
end interface threeshl_print_gwff
public :: threeshl_print_gwff
interface threeshl_print_gzff
 module procedure print_gzff
end interface threeshl_print_gzff
public :: threeshl_print_gzff
interface threeshl_print_gauge_coup
 module procedure print_gauge_coup
end interface threeshl_print_gauge_coup
public :: threeshl_print_gauge_coup
interface threeshl_print_parameters
 module procedure print_parameters
end interface threeshl_print_parameters
public :: threeshl_print_parameters

contains

subroutine set_names
 particle_names(e_bcd) = "electron"
 particle_names(nue_bcd) = "electron neutrino"
 particle_names(mu_bcd) = "muon"
 particle_names(numu_bcd) = "muon neutrino"
 particle_names(tau_bcd) = "tauon"
 particle_names(nutau_bcd) = "tauon neutrino"
 particle_names(he_bcd) = "heavy electron"
 particle_names(hnue_bcd) = "heavy electron neutrino"
 particle_names(hmu_bcd) = "heavy muon"
 particle_names(hnumu_bcd) = "heavy muon neutrino"
 particle_names(htau_bcd) = "heavy tauon"
 particle_names(hnutau_bcd) = "heavy tauon neutrino"
 particle_names(u_bcd) = "up quark"
 particle_names(d_bcd) = "down quark"
 particle_names(c_bcd) = "charm quark"
 particle_names(s_bcd) = "strange quark"
 particle_names(t_bcd) = "top quark"
 particle_names(b_bcd) = "bottom quark"
 particle_names(hu_bcd) = "heavy up quark"
 particle_names(hd_bcd) = "heavy down quark"
 particle_names(hc_bcd) = "heavy charm quark"
 particle_names(hs_bcd) = "heavy strange quark"
 particle_names(ht_bcd) = "heavy top quark"
 particle_names(hb_bcd) = "heavy bottom quark"
 particle_names(w_bcd) = "W boson"
 particle_names(hw_bcd) = "heavy W boson"
 particle_names(z_bcd) = "Z boson"
 particle_names(hz_bcd) = "heavy Z boson"
end subroutine set_names
function retrieve_bcd_f (kkmode, ftype, generation, isospin) result(bcd)
integer, intent(in) :: ftype, kkmode, isospin, generation
integer :: bcd
character(len=slength), parameter :: fname="retrieve_bcd_f"
 call errstack_push(fname)
 if ((ftype < ftype_l) .or. (ftype > ftype_q) .or. &
  (kkmode < l_mode) .or. (kkmode > h_mode) .or. (isospin < iso_up) .or. &
  (isospin > iso_down) .or. (generation < gen_0) .or. (generation > gen_2)) &
 call panic(err_invalid_parameters, 0)
 bcd = kkmode - l_mode + &
  ishft(ftype - ftype_l, 2) + &
  ishft(isospin - iso_up, 3) + &
  ishft(generation - gen_0, 4)
 call errstack_pop
end function retrieve_bcd_f
function retrieve_bcd_b (kkmode, btype) result(bcd)
integer, intent(in) :: kkmode, btype
integer :: bcd
character(len=slength), parameter :: fname="retrieve_bcd_b"
 call errstack_push(fname)
 if ((kkmode < l_mode) .or. (kkmode > h_mode) .or. &
  (btype < btype_w) .or. (btype > btype_z)) call panic(err_invalid_parameters, 0)
 bcd = kkmode - l_mode + 2 + ishft(btype - btype_w, 2)
 call errstack_pop
end function retrieve_bcd_b
subroutine decode_BCD(part, ptype_stat, kkmode, ptype, gen, isospin)
integer, intent(in) :: part
integer, intent(out) :: ptype_stat, kkmode, ptype, gen, isospin
character(len=slength), parameter :: fname="decode_BCD"
 call errstack_push(fname)
 if (ishft(part, -6) .ne. 0) call panic(err_invalid_parameters, 0)
 if (iand(part, 1) == 1) then
  kkmode = h_mode
 else
  kkmode = l_mode
 end if
 if (iand(part, 2) == 0) then
  ptype_stat = ptype_f
  if (iand(part, 4) == 0) then
   ptype = ftype_l
  else
   ptype = ftype_q
  end if
  if (iand(part, 8) == 0) then
   isospin = iso_up
  else
   isospin = iso_down
  end if
  if (ishft(part, -4) .gt. 2) call panic(err_invalid_parameters, 1)
  gen = ishft(part, -4) + gen_0
 else
  ptype_stat = ptype_b
  if (ishft(part, -3) .ne. 0) call panic(err_invalid_parameters, 2)
  if (iand(part, 4) == 0) then
   ptype = btype_w
  else
   ptype = btype_z
  end if
 end if
 call errstack_pop
end subroutine decode_BCD
function x2_of_mratio2 (mratio2) result (x2)
real(kind=double), intent(in) :: mratio2
real(kind=double) :: x2
character(len=slength), parameter :: fname="x2_of_mratio"
 call errstack_push(fname)
 if ((mratio2 > 1.) .or. (mratio2 < 0.)) &
  call panic(err_invalid_parameters, 0)
 x2 = (1._double - 2._double*mratio2 + mratio2**2 - &
  (1._double + mratio2)*msqrt(1._double - 6._double*mratio2 + mratio2**2)) / &
  2._double / mratio2
 if (x2 < 0.) call panic(err_range_error, 0)
 call errstack_pop
end function x2_of_mratio2
function c2_of_ctw2 (ctw2, x) result(c2)
real(kind=double), intent(in) :: ctw2, x
real(kind=double) :: c2, a
character(len=slength), parameter :: fname="c2_of_ctw2"
 call errstack_push(fname)
 if ((ctw2 > 1.) .or. (ctw2 < -1.) .or. (x < 0.)) &
  call panic(err_invalid_parameters, 0)
 a = sqrt(4._double + x**4)
 c2 =(ctw2*(x**2*(2._double + a + x**2) + 2._double*ctw2**3*(x**4 + x**6) -&
 ctw2*(2._double*(2._double + a) + (10._double + 3._double*a) *&
 x**2 + (3._double + a)*x**4 + x**6) + &
   ctw2**2*(4._double + 6._double*x**2 - x**6 + a*(2._double + 2._double*x**2 + x**4))))/&
   (2._double*(x**2 + ctw2**4*x**6 - 2._double*ctw2*(2._double + x**2) +&
 2._double*ctw2**3*x**2*(2._double + x**2) - ctw2**2 *&
 (-4._double + 4._double*x**2 + 2._double*x**4 + x**6)))
 if ((c2 > 1.) .or. (c2 < 0.)) call panic(err_range_error, 0)
 call errstack_pop
end function c2_of_ctw2
function eps_l_of_x (x) result(el)
real(kind=double), intent(in) :: x
real(kind=double) :: el, a, c
character(len=slength), parameter :: fname="eps_l_of_x"
 call errstack_push(fname)
 if (x < 0.) call panic(err_invalid_parameters, 0)
 a = sqrt(4._double + x**4)
 c = sqrt(8._double * (2._double + a) * (1._double + x**2) +&
  2._double * (5._double + 2._double*a) * x**4 +&
  (4._double + a) * x**6 + x**8)
 el = (16._double + 8._double*a - sqrt(2._double)*a*c -&
  2._double*sqrt(2._double)*c + (-4._double - 2._double*a +&
  3._double*sqrt(2._double)*c) * x**2 + 2._double*(3._double + a) *&
  x**4 + 2._double * x**6) / 8._double / (2._double + a) / x**2
 el = msqrt(el)
 call errstack_pop
end function eps_l_of_x
subroutine gauge_cpl_from_sm_wgap (mw, mz, e_, mwprime)
real(kind=double), intent(in) :: mw, mz, e_, mwprime
character(len=slength), parameter :: fname="gauge_cpl_from_sm_wgap"
 call errstack_push(fname)
 if ((mw < 0.) .or. (mz < 0.) .or. (mw > mz) .or. (e_ < 0.) &
  .or. (mwprime < 0.) .or. (mwprime < mw)) call panic(err_invalid_parameters, 0)
 x = msqrt(x2_of_mratio2(mw**2 / mwprime**2))
 t = msqrt(1._double/c2_of_ctw2(mw**2 / mz**2, x) - 1._double)
 e = e_
 g1 = e * msqrt(1._double + 1._double / x**2 * (1._double + 1._double / t**2))
 g0 = x*g1
 g2 = t*g0
 sigma_vev = 2._double * mw / g1 / &
  msqrt(2._double + x**2 - sqrt(4._double + x**4))
 call errstack_pop
end subroutine gauge_cpl_from_sm_wgap
function eps_r_of_m (m) result(er)
real(kind=double), intent(in) :: m
real(kind=double) :: er, a
character(len=slength), parameter :: fname="eps_r_of_m"
 call errstack_push(fname)
 if (m < 0.) call panic(err_invalid_parameters, 0)
 a = m**2 / 2._double / lambda**2 / sigma_vev**2
 er = msqrt(a * (1._double + eps_l**2 - a) / (eps_l**2 - a))
 call errstack_pop
end function eps_r_of_m
subroutine translate_fermion_masses(masses_l, masses_qu, masses_qd)
real(kind=double), intent(in), dimension(3) :: &
 masses_l, masses_qu, masses_qd
character(len=slength), parameter :: fname="translate_fermion_masses"
 call errstack_push(fname)
 if (minval((/masses_l, masses_qu, masses_qd/)) < 0.) &
  call panic(err_invalid_parameters, 0)
 do i = 1, 3
  eps_r(ftype_l, gen_0+i-1, iso_down) =&
   eps_r_of_m(masses_l(i))
  eps_r(ftype_q, gen_0+i-1, iso_up) =&
   eps_r_of_m(masses_qu(i))
  eps_r(ftype_q, gen_0+i-1, iso_down) =&
   eps_r_of_m(masses_qd(i))
 end do
 call errstack_pop
end subroutine translate_fermion_masses
subroutine calculate_fermion(er, wf_LL, wf_LR, wf_HL, wf_HR, m_L, m_H)
real(kind=double), intent(in) :: er
real(kind=double), intent(out), dimension(2) ::&
 wf_LL, wf_LR, wf_HL, wf_HR
real(kind=double), intent(out) :: m_L, m_H
real(kind=double) :: a
real(kind=double), dimension(2,2) :: ut, v, m, mdiag
character(len=slength), parameter :: fname="calculate_fermion"
 call errstack_push(fname)
 if (er < 0.) call panic(err_invalid_parameters, 0)
 if (er == 0.) then
  a = sqrt(eps_l**2 + 1._double)
  wf_LL = (/-1._double / a, eps_l / a/)
  wf_LR = (/0._double, 1._double/)
  wf_HL = (/eps_l / a, 1._double / a/)
  wf_HR = (/1._double, 0._double/)
 else
  a = msqrt (1._double + 2._double*eps_l**2 + eps_l**4 + &
   2._double*er**2 - 2._double*eps_l**2*er**2 + er**4)
  wf_LL = (/-1._double * (1._double - eps_l**2 + er**2 + a) /&
   2._double / eps_l, 1._double/)
  wf_LR = (/(1._double + eps_l**2 - er**2 - a) / 2._double / er, 1._double/)
  wf_HL = (/-1._double * (1._double - eps_l**2 + er**2 - a) /&
   2._double / eps_l, 1._double/)
  wf_HR = (/(1._double + eps_l**2 - er**2 + a) / 2._double / er, 1._double/)
  wf_LL = wf_LL / sqrt(wf_LL(1)**2 + wf_LL(2)**2)
  wf_LR = wf_LR / sqrt(wf_LR(1)**2 + wf_LR(2)**2)
  wf_HL = wf_HL / sqrt(wf_HL(1)**2 + wf_HL(2)**2)
  wf_HR = wf_HR / sqrt(wf_HR(1)**2 + wf_HR(2)**2)
 end if
 ut(:,1) = wf_LL
 ut(:,2) = wf_HL
 v(1,:) = wf_LR
 v(2,:) = wf_HR
 m = reshape((/eps_L, 1._double, 0._double, er/), (/2,2/))
 mdiag = matmul(ut, matmul(m, v))
 m_L = sqrt(2._double) * lambda * sigma_vev * mdiag(1,1)
 m_H = sqrt(2._double) * lambda * sigma_vev * mdiag(2,2)
 if ((m_L < 0.) .or. (m_H < 0.)) call panic(err_range_error, 0)
 call errstack_pop
end subroutine calculate_fermion
function get_gwff_l (wf_w, wf_fu, wf_fd) result(coup)
real(kind=double), intent(in), dimension(2) :: wf_w, wf_fu, wf_fd
real(kind=double) :: coup
 coup = g0 * (wf_fu(1) * wf_fd(1) * wf_w(1) + &
  wf_fu(2) * wf_fd(2) * wf_w(2) / x)
end function get_gwff_l
function get_gwff_r (wf_w, wf_fu, wf_fd) result(coup)
real(kind=double), intent(in), dimension(2) :: wf_w, wf_fu, wf_fd
real(kind=double) :: coup
 coup = g0 / x * wf_w(2) * wf_fu(1) * wf_fd(1)
end function get_gwff_r
function get_gzff_l (wf_z, kk, ftype, gen, isospin) result(coup)
real(kind=double), intent(in), dimension(3) :: wf_z
integer, intent(in) :: kk, ftype, gen, isospin
real(kind=double) :: coup, y, t3
real(kind=double), dimension(2) :: wf_f1, wf_f2
character(len=slength), parameter :: fname="get_gzff_l"
 call errstack_push(fname)
 if ((kk < l_mode) .or. (kk > lh_mode) .or. (ftype < ftype_l) &
  .or. (ftype > ftype_q) .or. (gen < gen_0) .or. (gen > gen_2) &
  .or. (isospin < iso_up) .or. (isospin > iso_down)) &
  call panic(err_invalid_parameters, 0)
 if (isospin == iso_up) then
  t3 = 0.5_double
 else
  t3 = -0.5_double
 end if
 if (ftype == ftype_l) then
  y = -0.5_double
  select case(kk)
   case(l_mode, h_mode)
    wf_f1 = wfunct_lep_l(kk, gen, isospin, :)
    wf_f2 = wf_f1
   case(lh_mode)
    wf_f1 = wfunct_lep_l(l_mode, gen, isospin, :)
    wf_f2 = wfunct_lep_l(h_mode, gen, isospin, :)
   end select
 else
  y = 1._double / 6._double
  select case(kk)
   case(l_mode, h_mode)
    wf_f1 = wfunct_quark_l(kk, gen, isospin, :)
    wf_f2 = wf_f1
   case(lh_mode)
    wf_f1 = wfunct_quark_l(l_mode, gen, isospin, :)
    wf_f2 = wfunct_quark_l(h_mode, gen, isospin, :)
   end select
 end if
 coup = &
  wf_f1(1) * wf_f2(1) * (g0 * wf_z(1) * t3 + g2 * wf_z(3) * y) +&
  wf_f1(2) * wf_f2(2) * (g1 * wf_z(2) * t3 + g2 * wf_z(3) * y)
 call errstack_pop
end function get_gzff_l
function get_gzff_r (wf_z, kk, ftype, gen, isospin) result(coup)
real(kind=double), intent(in), dimension(3) :: wf_z
integer, intent(in) :: kk, ftype, gen, isospin
real(kind=double) :: coup, y, t3
real(kind=double), dimension(2) :: wf_f1, wf_f2
character(len=slength), parameter :: fname="get_gzff_r"
 call errstack_push(fname)
 if ((kk < l_mode) .or. (kk > lh_mode) .or. (ftype < ftype_l) &
  .or. (ftype > ftype_q) .or. (gen < gen_0) .or. (gen > gen_2) &
  .or. (isospin < iso_up) .or. (isospin > iso_down)) &
  call panic(err_invalid_parameters, 0)
 if (isospin == iso_up) then
  t3 = 0.5_double
 else
  t3 = -0.5_double
 end if
 if (ftype == ftype_l) then
  y = -0.5_double
  select case(kk)
   case(l_mode, h_mode)
    wf_f1 = wfunct_lep_r(kk, gen, isospin, :)
    wf_f2 = wf_f1
   case(lh_mode)
    wf_f1 = wfunct_lep_r(l_mode, gen, isospin, :)
    wf_f2 = wfunct_lep_r(h_mode, gen, isospin, :)
   end select
 else
  y = 1._double / 6._double
  select case(kk)
   case(l_mode, h_mode)
    wf_f1 = wfunct_quark_r(kk, gen, isospin, :)
    wf_f2 = wf_f1
   case(lh_mode)
    wf_f1 = wfunct_quark_r(l_mode, gen, isospin, :)
    wf_f2 = wfunct_quark_r(h_mode, gen, isospin, :)
   end select
 end if
 coup = &
  wf_f1(1) * wf_f2(1) * (g1 * wf_z(2) - g2 * wf_z(3)) * t3 +&
  g2 * wf_z(3) * (y + t3) * (wf_f1(1)*wf_f2(1) + wf_f1(2)*wf_f2(2))
 call errstack_pop
end function get_gzff_r
function get_gwwzz(kk1, kk2) result(coup)
integer, intent(in) :: kk1, kk2
real(kind=double) :: coup
real(kind=double), dimension(2) :: wf_w1, wf_w2
real(kind=double), dimension(3) :: wf_z1, wf_z2
character(len=slength), parameter :: fname="get_gwwzz"
 call errstack_push(fname)
 select case(kk1)
  case(l_mode, h_mode)
   wf_w1 = wfunct_w(kk1, :)
   wf_w2 = wf_w1
  case(lh_mode)
   wf_w1 = wfunct_w(l_mode, :)
   wf_w2 = wfunct_w(h_mode, :)
  case default
   call panic(err_invalid_parameters, 0)
 end select
 select case(kk2)
  case(l_mode, h_mode)
   wf_z1 = wfunct_z(kk2, :)
   wf_z2 = wf_z1
  case(lh_mode)
   wf_z1 = wfunct_z(l_mode, :)
   wf_z2 = wfunct_z(h_mode, :)
  case default
   call panic(err_invalid_parameters, 0)
 end select
 coup = g0**2 * wf_w1(1) * wf_w2(1) * wf_z1(1) * wf_z2(1) +&
  g1**2 * wf_w1(2) * wf_w2(2) * wf_z1(2) * wf_z2(2)
 call errstack_pop
end function get_gwwzz
subroutine calculate_couplings
integer :: isospin, gen, kk1, kk2, kk3
 call errstack_push( pad_string( "calculate_couplings" ))
 do gen = gen_0, gen_2
 do kk1 = l_mode, h_mode
 do kk2 = l_mode, h_mode
 do kk3 = l_mode, h_mode
  g_w_lep(kk1, kk2, gen, kk3, gen, l_chir) = get_gwff_l( &
   wfunct_w(kk1, :), wfunct_lep_l(kk2, gen, iso_up, :), &
   wfunct_lep_l(kk3, gen, iso_down, :))
  g_w_lep(kk1, kk2, gen, kk3, gen, r_chir) = get_gwff_r( &
   wfunct_w(kk1, :), wfunct_lep_r(kk2, gen, iso_up, :), &
   wfunct_lep_r(kk3, gen, iso_down, :))
  g_w_quark(kk1, kk2, gen, kk3, gen, l_chir) = get_gwff_l( &
   wfunct_w(kk1, :), wfunct_quark_l(kk2, gen, iso_up, :), &
   wfunct_quark_l(kk3, gen, iso_down, :))
  g_w_quark(kk1, kk2, gen, kk3, gen, r_chir) = get_gwff_r( &
   wfunct_w(kk1, :), wfunct_quark_r(kk2, gen, iso_up, :), &
   wfunct_quark_r(kk3, gen, iso_down, :))
 end do
 end do
 end do
 end do
 do kk1 = l_mode, h_mode
 do kk2 = l_mode, lh_mode
 do isospin = iso_up, iso_down
 do gen = gen_0, gen_2
  g_z_lep(kk1, kk2, gen, isospin, l_chir) = get_gzff_l( &
   wfunct_z(kk1, :), kk2, ftype_l, gen, isospin)
  g_z_lep(kk1, kk2, gen, isospin, r_chir) = get_gzff_r( &
   wfunct_z(kk1, :), kk2, ftype_l, gen, isospin)
  g_z_quark(kk1, kk2, gen, isospin, l_chir) = get_gzff_l( &
   wfunct_z(kk1, :), kk2, ftype_q, gen, isospin)
  g_z_quark(kk1, kk2, gen, isospin, r_chir) = get_gzff_r( &
   wfunct_z(kk1, :), kk2, ftype_q, gen, isospin)
 end do
 end do
 end do
 end do
 do kk1 = l_mode, h_mode
  do kk2 = l_mode, h_mode
   g_wwz(kk2, kk1) = &
    g0 * wfunct_w(kk2, lat_0)**2 * wfunct_z(kk1, lat_0) +&
    g1 * wfunct_w(kk2, lat_1)**2 * wfunct_z(kk1, lat_1)
   g_wwza(kk2, kk1) = &
    g0 * e * wfunct_w(kk2, lat_0)**2 * wfunct_z(kk1, lat_0) +&
    g1 * e * wfunct_w(kk2, lat_1)**2 * wfunct_z(kk1, lat_1)
  end do
  g_wwz(lh_mode, kk1) = &
   g0 * wfunct_w(l_mode, lat_0) * wfunct_w(h_mode, lat_0) &
    * wfunct_z(kk1, lat_0) +&
   g1 * wfunct_w(l_mode, lat_1) * wfunct_w(h_mode, lat_1) &
    * wfunct_z(kk1, lat_1)
  g_wwza(lh_mode, kk1) = &
   g0 * e * wfunct_w(l_mode, lat_0) * wfunct_w(h_mode, lat_0) *&
    wfunct_z(kk1, lat_0) +&
   g1 * e * wfunct_w(l_mode, lat_1) * wfunct_w(h_mode, lat_1) *&
    wfunct_z(kk1, lat_1)
 end do
 do i = 0, 4
  g_wwww(i) = &
   g0**2 * wfunct_w(l_mode, lat_0)**(4-i) * wfunct_w(h_mode, lat_0)**i +&
   g1**2 * wfunct_w(l_mode, lat_1)**(4-i) * wfunct_w(h_mode, lat_1)**i
 end do
 do kk1 = l_mode, lh_mode
 do kk2 = l_mode, lh_mode
  g_wwzz(kk1, kk2) = get_gwwzz(kk1, kk2)
 end do
 end do
 call errstack_pop
end subroutine calculate_couplings
subroutine diagonalize
real(kind=double) :: a
integer :: isospin, gen
character(len=slength), parameter :: fname="diagonalize"
 call errstack_push(fname)
 mass_array(w_bcd) = sigma_vev * 0.5_double * g1 *&
  msqrt(2._double + x**2 - sqrt(4._double + x**4))
 mass_array(hw_bcd) = sigma_vev * 0.5_double * g1 *&
  msqrt(2._double + x**2 + sqrt(4._double + x**4))
 mass_array(z_bcd) = sigma_vev * 0.5_double * g1 *&
 msqrt(2._double + (1._double + t**2)*x**2 -&
  msqrt(4._double + (t**2 - 1._double)**2 * x**4))
 mass_array(hz_bcd) = sigma_vev * 0.5_double * g1 *&
  msqrt(2._double + (1._double + t**2)*x**2 +&
  msqrt(4._double + (t**2 - 1._double)**2 * x**4))
 wfunct_w(L_mode,:) = (/ (2._double - x**2 +&
  sqrt(4._double + x**4)) / 2._double / x, 1._double /)
 wfunct_w(L_mode,:) = wfunct_w(L_mode,:) / &
  sqrt(wfunct_w(L_mode, lat_0)**2 + wfunct_w(L_mode, lat_1)**2)
 wfunct_w(H_mode,:) = (/ (2._double - x**2 -&
  sqrt(4._double + x**4)) / 2._double / x, 1._double /)
 wfunct_w(H_mode,:) = wfunct_w(H_mode,:) / &
  sqrt(wfunct_w(H_mode, lat_0)**2 + wfunct_w(H_mode, lat_1)**2)
 a = msqrt(4._double + x**4 - 2._double*t**2*x**4 + t**4*x**4)
 wfunct_z(L_mode,:) = (/ (x**2 - t**2*x**2 - a) / 2._double / t, &
  (-2._double - x**2 + t**2*x**2 + a) / 2._double / t / x, 1._double/)
 wfunct_z(L_mode,:) = wfunct_z(L_mode,:) / sqrt(wfunct_z(L_mode, lat_0)**2 + &
  wfunct_z(L_mode, lat_1)**2 + wfunct_z(L_mode, lat_2)**2)
 wfunct_z(H_mode,:) = (/ (x**2 - t**2*x**2 + a) / 2._double / t, &
  (-2._double - x**2 + t**2*x**2 - a) / 2._double / t / x, 1._double/)
 wfunct_z(H_mode,:) = wfunct_z(H_mode,:) / sqrt(wfunct_z(H_mode, lat_0)**2 + &
  wfunct_z(H_mode, lat_1)**2 + wfunct_z(H_mode, lat_2)**2)
 do isospin = iso_up, iso_down
  do gen = gen_0, gen_2
   call calculate_fermion(eps_r(ftype_l, gen, isospin), &
    wfunct_lep_l(l_mode, gen, isospin, :), &
    wfunct_lep_r(l_mode, gen, isospin, :), &
    wfunct_lep_l(h_mode, gen, isospin, :), &
    wfunct_lep_r(h_mode, gen, isospin, :), &
    mass_array(retrieve_bcd_f(l_mode, ftype_l, gen, isospin)), &
    mass_array(retrieve_bcd_f(h_mode, ftype_l, gen, isospin)))
   call calculate_fermion(eps_r(ftype_q, gen, isospin), &
    wfunct_quark_l(l_mode, gen, isospin, :), &
    wfunct_quark_r(l_mode, gen, isospin, :), &
    wfunct_quark_l(h_mode, gen, isospin, :), &
    wfunct_quark_r(h_mode, gen, isospin, :), &
    mass_array(retrieve_bcd_f(l_mode, ftype_q, gen, isospin)), &
    mass_array(retrieve_bcd_f(h_mode, ftype_q, gen, isospin)))
  end do
 end do
 call calculate_couplings
 call errstack_pop
end subroutine diagonalize
subroutine pdg_init_wgap_bmass (mhw, bmass, el)
real(kind=double), intent(in) :: mhw, bmass
real(kind=double), optional, intent(in) :: el
character(len=slength), parameter :: fname="pdb_init_wgap_bmass"
 call errstack_push(fname)
 if ((mhw < mw_pdg) .or. (bmass < 0.)) call panic(err_invalid_parameters, 0)
 call gauge_cpl_from_sm_wgap(mw_pdg, mz_pdg, e_pdg, mhw)
 if (present (el)) then
  eps_l = el
 else
  eps_l = eps_l_of_x(x)
 end if
 lambda = bmass / sqrt(2._double) / sigma_vev
 call translate_fermion_masses( (/me_pdg, mmu_pdg, mtau_pdg/), &
  (/muq_pdg, mcq_pdg, mtq_pdg/), (/mdq_pdg, msq_pdg, mbq_pdg/))
 call diagonalize
 call errstack_pop
end subroutine pdg_init_wgap_bmass
subroutine calculate_widths
character(len=slength), parameter :: fname="calculate_widths"
integer :: kk1, kk2, kk3, gen1, gen2, iso1, iso2, bcd, ftype
call errstack_push (fname)
width_array = 0
width_array(w_bcd) = ww_pdg
width_array(z_bcd) = wz_pdg
width_array(t_bcd) = wt_pdg
width_array (hw_bcd) = &
 wd_ggg (mass_array(hw_bcd), mass_array(w_bcd), mass_array(z_bcd), &
  g_wwz(lh_mode, l_mode)) / 3._double + &
 wd_ggg (mass_array(hw_bcd), mass_array(w_bcd), mass_array(hz_bcd), &
  g_wwz(lh_mode, h_mode)) / 3._double
 do kk1 = l_mode, h_mode; do kk2 = l_mode, h_mode
 do gen1 = gen_0, gen_2; do gen2 = gen_0, gen_2
  width_array (hw_bcd) = width_array(hw_bcd) + &
   wd_gff (mass_array(hw_bcd), &
   mass_array(retrieve_bcd_f (kk1, ftype_l, gen1, iso_up)), &
   mass_array(retrieve_bcd_f (kk2, ftype_l, gen2, iso_down)), &
   g_w_lep(h_mode, kk1, gen1, kk2, gen2, :) / sqrt(2._double)) / 3._double &
  + &
   wd_gff (mass_array(hw_bcd), &
   mass_array(retrieve_bcd_f (kk1, ftype_q, gen1, iso_up)), &
   mass_array(retrieve_bcd_f (kk2, ftype_q, gen2, iso_down)), &
   g_w_quark(h_mode, kk1, gen1, kk2, gen2, :) / sqrt(2._double))
 end do; end do; end do; end do
call errstack_pop
do kk2 = l_mode, h_mode; do kk3 = l_mode, h_mode
 if (kk2 == kk3) then
  kk1 = kk2
 else
  kk1 = lh_mode
 end if
 width_array(hz_bcd) = width_array(hz_bcd) + &
  wd_ggg (mass_array(hz_bcd), &
   mass_array(retrieve_bcd_b (kk2, btype_w)), mass_array(retrieve_bcd_b (kk3, btype_w)), &
   g_wwz (kk1, h_mode)) /3._double
 do gen1 = gen_0, gen_2; do iso1 = iso_up, iso_down
  width_array(hz_bcd) = width_array(hz_bcd) + wd_gff (mass_array(hz_bcd), &
   mass_array(retrieve_bcd_f (kk2, ftype_l, gen1, iso1)), &
   mass_array(retrieve_bcd_f (kk3, ftype_l, gen1, iso1)), &
   g_z_lep (h_mode, kk1, gen1, iso1, :)) / 3._double &
  + wd_gff (mass_array(hz_bcd), &
   mass_array(retrieve_bcd_f (kk2, ftype_q, gen1, iso1)), &
   mass_array(retrieve_bcd_f (kk3, ftype_q, gen1, iso1)), &
   g_z_quark (h_mode, kk1, gen1, iso1, :))
 end do; end do
end do; end do

nlow_mc = mass_array(c_bcd)
nlow_mb = mass_array(b_bcd)
nlow_mt = mass_array(t_bcd)

do gen1 = gen_0, gen_2; do iso1 = iso_up, iso_down
 bcd = retrieve_bcd_f (h_mode, ftype_l, gen1, iso1)
 do kk1 = l_mode, h_mode
  width_array(bcd) = width_array(bcd) + wd_fgf (mass_array(bcd), &
    mass_array(retrieve_bcd_b (kk1, btype_z)), &
    mass_array(retrieve_bcd_f (l_mode, ftype_l, gen1, iso1)), &
    cpl_hfzf (ftype_l, gen1, iso1, kk1, l_mode)) / 2._double
  do gen2 = gen_0, gen_2; do kk2 = l_mode, h_mode
   width_array(bcd) = width_array(bcd) + wd_fgf (mass_array(bcd), &
    mass_array(retrieve_bcd_b (kk1, btype_w)), &
    mass_array(retrieve_bcd_f (kk2, ftype_l, gen2, flip_isospin(iso1))), &
    cpl_hfwf (ftype_l, gen1, iso1, kk1, kk2, gen2) / sqrt(2._double)) / 2._double
  end do; end do
 end do
 bcd = retrieve_bcd_f (h_mode, ftype_q, gen1, iso1)
 do kk1 = l_mode, h_mode
  width_array(bcd) = width_array(bcd) + wd_qgq (mass_array(bcd), &
    mass_array(retrieve_bcd_b (kk1, btype_z)), &
    mass_array(retrieve_bcd_f (l_mode, ftype_q, gen1, iso1)), &
    cpl_hfzf (ftype_q, gen1, iso1, kk1, l_mode)) / 2._double
  do gen2 = gen_0, gen_2; do kk2 = l_mode, h_mode
   width_array(bcd) = width_array(bcd) + wd_qgq (mass_array(bcd), &
    mass_array(retrieve_bcd_b (kk1, btype_w)), &
    mass_array(retrieve_bcd_f (kk2, ftype_q, gen2, flip_isospin(iso1))), &
    cpl_hfwf (ftype_q, gen1, iso1, kk1, kk2, gen2) / sqrt(2._double)) / 2._double
  end do; end do
 end do
end do; end do
contains
function width_normal (m1, m2, m3) result(res)
real(kind=double), intent(in) :: m1, m2, m3
real(kind=double) :: res
 if (m1 < (m2 + m3)) then
  res = 0
 else
  res = msqrt((m1**2 - (m2 + m3)**2) * (m1**2 - (m2 - m3)**2)) / &
   16._double / pi / m1**3
 end if
end function width_normal
function wd_ggg (m1, m2, m3, g) result(amp)
real(kind=double), intent(in) :: m1, m2, m3, g
real(kind=double) :: amp
 if (min (m1, m2, m3) < 0) call panic (err_invalid_parameters, 2)
 if (m1 < (m2 + m3)) then
  amp = 0
 else
  amp = g**2 * ( -8._double*(m1**2 + m2**2 + m3**2) + &
   2._double*((m3**4 + m2**4)/m1**2 + (m1**4 + m3**4)/m2**2 + &
   (m1**4 + m2**4)/m3**2) - 4.5_double*( &
   m2**2*m3**2/m1**2 + m1**2*m3**2/m2**2 + m1**2*m2**2/m3**2) +&
   0.25_double*(m1**6/m2**2/m3**2 + m2**6/m1**2/m3**2 + m3**6/m1**2/m2**2)) *&
   width_normal (m1, m2, m3)
 end if
end function wd_ggg
function wd_gff (m1, m2, m3, glr) result(amp)
real(kind=double), intent(in) :: m1, m2, m3, glr(2)
real(kind=double) :: amp, gva(2)
 if (min (m1, m2, m3) < 0) call panic (err_invalid_parameters, 3)
 if (m1 < (m2 + m3)) then
  amp = 0
 else
  gva = (/glr(1) + glr(2), glr(1) - glr(2)/)/2._double
  amp = ((gva(1)**2 + gva(2)**2) * ( &
   2._double*(2._double*m1**2 - m2**2 - m3**2) -&
   2._double*(m2**4 + m3**4)/m1**2 + 4._double*m2**2*m3**2/m1**2) +&
   12._double * (gva(1)**2 - gva(2)**2) *m2*m3) *&
   width_normal (m1, m2, m3)
 end if
end function wd_gff
function wd_fgf (m1, m2, m3, glr) result(amp)
real(kind=double), intent(in) :: m1, m2, m3, glr(2)
real(kind=double) :: amp, gva(2)
 if (min (m1, m2, m3) < 0) call panic(err_invalid_parameters, 4)
 if (m1 < (m2 + m3)) then
  amp = 0
 else
  gva = (/glr(1) + glr(2), glr(1) - glr(2)/)/2._double
  amp = ((gva(1)**2 + gva(2)**2) * ( &
   2._double*(m1**2 + m3**2 - 2._double*m2**2) +&
   2._double*(m1**4 + m3**4)/m2**2 - 4._double*(m1**2*m3**2)/m2**2) -&
   12._double*(gva(1)**2 - gva(2)**2)*m1*m3) *&
   width_normal (m1, m2, m3)
 end if
end function wd_fgf
function wd_qgq (m1, m2, m3, glr) result(amp)
real(kind=double), intent(in) :: m1, m2, m3, glr(2)
real(kind=double) :: amp, gva(2)
logical, save :: warned=.false.
 if (min (m1, m2, m3) < 0) call panic(err_invalid_parameters, 4)
 if (m1 < (m2 + m3)) then
  amp = 0
 else
  gva = (/glr(1) + glr(2), glr(1) - glr(2)/)/2._double

  if (threeshl_use_nlow .and. m3 > 0.0001_double) then
   amp = nlow_width (m1, m3, m2, gva(1), gva(2), &
    nlow_alfas (sqrt (2._double) * sigma_vev * lambda), 1._double) * 2._double
  elseif (threeshl_use_nlow) then
   amp = nlow_width (m1, m2, gva(1), gva(2), &
    nlow_alfas (sqrt (2._double) * sigma_vev * lambda), 1._double) * 2._double
  else

   amp = wd_fgf (m1, m2, m3, glr)

  end if
 end if
end function wd_qgq
function cpl_hfzf (ftype_hf, gen_hf, iso_hf, kk_b, kk_f) result(cpl)
integer, intent(in) :: gen_hf, iso_hf, kk_b, kk_f, ftype_hf
real(kind=double), dimension(2) :: cpl
integer :: kkcomb
 select case(kk_f)
  case(h_mode); kkcomb = h_mode
  case(l_mode); kkcomb = lh_mode
 end select
 select case(ftype_hf)
  case(ftype_l)
   cpl = g_z_lep (kk_b, kkcomb, gen_hf, iso_hf, :)
  case(ftype_q)
   cpl = g_z_quark (kk_b, kkcomb, gen_hf, iso_hf, :)
  case default
   call panic (err_invalid_parameters, 5)
 end select
end function cpl_hfzf
function cpl_hfwf (ftype_hf, gen_hf, iso_hf, kk_b, kk_f, gen_f) result(cpl)
integer, intent(in) :: ftype_hf, gen_hf, iso_hf, kk_b, kk_f, gen_f
real(kind=double), dimension(2) :: cpl
integer :: kk_up, kk_dwn, gen_up, gen_dwn
 select case(iso_hf)
  case(iso_up)
   kk_up = h_mode; kk_dwn = kk_f
   gen_up = gen_hf; gen_dwn = gen_f
  case(iso_down)
   kk_up = kk_f; kk_dwn = h_mode
   gen_up = gen_f; gen_dwn = gen_hf
  case default
   call panic (err_invalid_parameters, 6)
 end select
 select case(ftype_hf)
  case(ftype_l)
   cpl = g_w_lep (kk_b, kk_up, gen_up, kk_dwn, gen_dwn, :)
  case(ftype_q)
   cpl = g_w_quark (kk_b, kk_up, gen_up, kk_dwn, gen_dwn, :)
  case default
   call panic (err_invalid_parameters, 6)
 end select
end function cpl_hfwf
function flip_isospin (ispin) result(res)
integer, intent(in) :: ispin
integer :: res
 select case(ispin)
  case(iso_up); res = iso_down
  case(iso_down); res = iso_up
  case default; call panic (err_invalid_parameters, 7)
 end select
end function flip_isospin
end subroutine calculate_widths
subroutine init_ward (mx, ct, ph)
real(kind=double), intent(in) :: mx, ct, ph
integer :: gen, iso
character(len=slength), parameter :: fname="init_ward"
 call errstack_push(fname)
 if ( (ct .le. -1._double) .or. (ct .ge. 1._double) .or. (mx == 0.)) &
  call panic (err_invalid_parameters, 0)
 mass_array = 0.
 t = msqrt(1._double/ct**2 - 1._double)
 x = mx
 e = e_pdg
 g0 = e * msqrt( 1._double + x**2 + 1._double/t**2 )
 g1 = g0 / x
 g2 = g0 * t
 sigma_vev = 0.
 lambda = 0.
 eps_l = 0.
 eps_r = 0.
 wfunct_w(l_mode, :) = (/cos(ph), sin(ph)/)
 wfunct_w(h_mode, :) = (/-sin(ph), cos(ph)/)
 wfunct_z(l_mode, :) = (/-g2/2._double/g1 - g1/g2 , &
  g2/2._double/g0 - g0/2._double/g2 , g1/g0 + g0/2._double/g1 /)
 wfunct_z(l_mode, :) = wfunct_z(l_mode, :) / &
  msqrt(wfunct_z(l_mode, lat_0)**2 + wfunct_z(l_mode, lat_1)**2 + &
   wfunct_z(l_mode, lat_2)**2)
 wfunct_z(h_mode, :) = &
  (/-g0/2._double, g1, -g2/2._double/) / &
  msqrt(g0**2/4._double + g1**2 + g2**2/4._double)
 do gen = gen_0, gen_2 ; do iso = iso_up, iso_down
  wfunct_lep_l(l_mode, gen, iso, :) = (/cos(ph), sin(ph)/)
  wfunct_lep_l(h_mode, gen, iso, :) = (/-sin(ph), cos(ph)/)
  wfunct_lep_r(l_mode, gen, iso, :) = (/sin(ph), cos(ph)/)
  wfunct_lep_r(h_mode, gen ,iso, :) = (/cos(ph), -sin(ph)/)
  wfunct_quark_l(l_mode, gen, iso, :) = (/cos(ph), sin(ph)/)
  wfunct_quark_l(h_mode, gen, iso, :) = (/-sin(ph), cos(ph)/)
  wfunct_quark_r(l_mode, gen, iso, :) = (/sin(ph), cos(ph)/)
  wfunct_quark_r(h_mode, gen ,iso, :) = (/cos(ph),-sin(ph)/)
 end do; end do
 call calculate_couplings
 call errstack_pop
end subroutine init_ward
subroutine init_pointers
 threeshl_sigma_vev => sigma_vev
 threeshl_g0 => g0
 threeshl_g1 => g1
 threeshl_g2 => g2
 threeshl_x => x
 threeshl_t => t
 threeshl_lambda => lambda
 threeshl_eps_l => eps_l
 threeshl_e => e
end subroutine init_pointers
subroutine init
 call nlow_init
 call set_names
 call init_pointers
end subroutine init
subroutine finalize
 call nlow_finalize
end subroutine finalize
function threeshl_eps_r (ftype, gen, isospin) result(eps)
integer, intent(in) :: ftype, gen, isospin
real(kind=double), pointer :: eps
 eps => eps_r(ftype, gen, isospin)
end function threeshl_eps_r
function msqrt (x) result(res)
real(kind=double), intent(in) :: x
real(kind=double) :: res
 if (x < 0) call panic(err_invalid_root, 0)
 res = sqrt(x)
end function msqrt
subroutine panic (err_id, err_supplement)
integer, intent(in) :: err_id, err_supplement
 if (threeshl_print_errors) then
  write (threeshl_errunit, *) "-----"
  write (threeshl_errunit, *) "error occured during diagonalization; message:"
  if ((err_id > err_num_ids) .or. (err_id < 1)) then
   write (threeshl_errunit, *) "[invalid error message]"
  else
   write (threeshl_errunit, *) trim(err_messages(err_id))
  end if
  write (threeshl_errunit, *)
  write (threeshl_errunit, *) "supplemental ID: ", err_supplement
  write (threeshl_errunit, *)
  write (threeshl_errunit, *) "function history:"
  if (errstack_pos == 0) then
   write (threeshl_errunit, *) "[empty]"
  else
   do i = errstack_pos, 1, -1
    write (threeshl_errunit, *) trim(errstack(i))
   end do
  end if
  write (threeshl_errunit, *) "-----"
 end if
 threeshl_error = .true.
 if (threeshl_quit_on_panic) call exit(1)
end subroutine panic
subroutine errstack_push(str)
character(len=slength), intent(in) :: str
 if (errstack_pos == errstack_size) then
  write (threeshl_errunit, *) "FATAL: history stack overflow!"
  call exit(1)
 else
  errstack_pos = errstack_pos + 1
  errstack(errstack_pos) = str
 end if
end subroutine errstack_push
subroutine errstack_pop
 if (errstack_pos == 0) then
  write (threeshl_errunit, *) "FATAL: history stack underflow!"
  call exit(1)
 else
  errstack_pos = errstack_pos - 1
 end if
end subroutine errstack_pop
function kkmode_to_text (mode) result(str)
integer, intent(in) :: mode
character(len=slength) :: str
character(len=slength), parameter :: fname="kkmode_to_text"
 call errstack_push(fname)
 if ((mode > h_mode) .or. (mode < l_mode)) call panic(err_invalid_parameters, 0)
 if (mode == l_mode) then
  str = "light"
 else
  str = "heavy"
 end if
 call errstack_pop
end function kkmode_to_text
subroutine print_particle_properties (part, unit)
integer, intent(in) :: part, unit
integer :: ptype_stat, ptype, kkmode, gen, isospin
character(len=slength) :: ptype_char, kkmode_char, isospin_char
character(len=slength), parameter :: fname="print_particle_properties"
 call errstack_push(fname)
 call decode_BCD(part, ptype_stat, kkmode, ptype, gen, isospin)
 if (kkmode == h_mode) then
  kkmode_char = "heavy"
 else
  kkmode_char = "light"
 end if
 if (ptype_stat == ptype_f) then
  if (ptype == ftype_l) then
   ptype_char = "lepton"
  else
   ptype_char = "quark"
  end if
  if (isospin == iso_up) then
   isospin_char = "+1/2"
  else
   isospin_char = "-1/2"
  end if
  write (unit, '(A," : ",A," ",A," with isospin ",A,&
   &" belonging to generation ",I1.1)') trim(particle_names(part)), &
   trim(kkmode_char), trim(ptype_char), trim(isospin_char), gen-gen_0
  write (unit, '(" eps_r: ",F15.7)') eps_r(ptype, gen, isospin)
  if (ptype == ftype_l) then
   write (unit, '(" wavefunction L :",F15.7,5X,F15.7)') &
    wfunct_lep_L(kkmode, gen, isospin, lat_0), &
    wfunct_lep_L(kkmode, gen, isospin, lat_1)
   write (unit, '(" wavefunction R :",F15.7,5X,F15.7)') &
    wfunct_lep_R(kkmode, gen, isospin, lat_1), &
    wfunct_lep_R(kkmode, gen, isospin, lat_2)
  else
   write (unit, '(" wavefunction L :",F15.7,5X,F15.7)') &
    wfunct_quark_L(kkmode, gen, isospin, lat_0), &
    wfunct_quark_L(kkmode, gen, isospin, lat_1)
   write (unit, '(" wavefunction R :",F15.7,5X,F15.7)') &
    wfunct_quark_R(kkmode, gen, isospin, lat_1), &
    wfunct_quark_R(kkmode, gen, isospin, lat_2)
  end if
 else
  if (ptype == btype_w) then
   ptype_char = "W boson"
  else
   ptype_char = "Z boson"
  end if
  write (unit, '(A," : ",A," ",A)') trim(particle_names(part)), &
   trim(kkmode_char), trim(ptype_char)
  if (ptype == btype_w) then
   write (unit, '(" wavefuntion: ",F15.7,5X,F15.7)') &
    wfunct_w(kkmode, lat_0), wfunct_w(kkmode, lat_1)
  else
   write (unit, '(" wavefuntion: ",F15.7,5X,F15.7,5X,F15.7)') &
    wfunct_z(kkmode, lat_0), wfunct_z(kkmode, lat_1), &
    wfunct_z(kkmode, lat_2)
  end if
 end if
 write (unit, '(" mass: ",F15.7)') mass_array(part)
 write (unit, '(" width:",F15.7)') width_array(part)
 call errstack_pop
end subroutine print_particle_properties
subroutine print_particles (unit)
integer, intent(in), optional :: unit
integer :: uunit
character(len=slength), parameter :: fname="output_particles"
 call errstack_push(fname)
 if (present (unit)) then
  uunit = unit
 else
  uunit = output_unit
 end if
 write (uunit, *)
 do i = 0, 11
  call print_particle_properties(ishft(i, 2), uunit)
  write (uunit, *)
  call print_particle_properties(ishft(i, 2) + 1, uunit)
  write (uunit, *)
 end do
 call print_particle_properties(w_BCD, uunit)
 write (uunit, *)
 call print_particle_properties(hw_BCD, uunit)
 write (uunit, *)
 call print_particle_properties(z_BCD, uunit)
 write (uunit, *)
 call print_particle_properties(hz_BCD, uunit)
 call errstack_pop
end subroutine print_particles
subroutine print_gwff (unit)
integer, intent(in), optional :: unit
integer :: kk1, kk2, kk3, gen, uunit
character(len=slength), parameter :: fname="output_gwff"
 call errstack_push(fname)
 if (present (unit)) then
  uunit = unit
 else
  uunit = output_unit
 end if
 do kk1 = l_mode, h_mode
 do kk2 = l_mode, h_mode
 do kk3 = l_mode, h_mode
  write (uunit, *)
  write (uunit, '("-------------------------------------------")')
  write (uunit, '("W boson:",T25,A)') trim(kkmode_to_text(kk1))
  write (uunit, '("isospin +1/2 lepton:",T25,A)') trim(kkmode_to_text(kk2))
  write (uunit, '("isospin -1/2 lepton:",T25,A)') trim(kkmode_to_text(kk3))
  write (uunit, *)
  write (uunit, '("left-handed coupling matrix (rows <-> iso-up, &
   & columns <-> iso-down):")')
  write (uunit, *)
  do gen = gen_0, gen_2
   write (uunit, '(5X,F15.7,5X,F15.7,5X,F15.7)') &
    g_w_lep(kk1, kk2, gen, kk3, gen_0, l_chir), &
    g_w_lep(kk1, kk2, gen, kk3, gen_1, l_chir), &
    g_w_lep(kk1, kk2, gen, kk3, gen_2, l_chir)
  end do
  write (uunit, *)
  write (uunit, '("right-handed coupling matrix (rows <-> iso-up, &
   & columns <-> iso-down):")')
  write (uunit, *)
  do gen = gen_0, gen_2
   write (uunit, '(5X,F15.7,5X,F15.7,5X,F15.7)') &
    g_w_lep(kk1, kk2, gen, kk3, gen_0, r_chir), &
    g_w_lep(kk1, kk2, gen, kk3, gen_1, r_chir), &
    g_w_lep(kk1, kk2, gen, kk3, gen_2, r_chir)
  end do
  write (uunit, *)
  write (uunit, '("-------------------------------------------")')
  write (uunit, *)
  write (uunit, *)
  write (uunit, '("-------------------------------------------")')
  write (uunit, '("W boson:",T25,A)') trim(kkmode_to_text(kk1))
  write (uunit, '("isospin +1/2 quark:",T25,A)') trim(kkmode_to_text(kk2))
  write (uunit, '("isospin -1/2 quark:",T25,A)') trim(kkmode_to_text(kk3))
  write (uunit, *)
  write (uunit, '("left-handed coupling matrix (rows <-> iso-up, &
   & columns <-> iso-down):")')
  write (uunit, *)
  do gen = gen_0, gen_2
   write (uunit, '(5X,F15.7,5X,F15.7,5X,F15.7)') &
    g_w_quark(kk1, kk2, gen, kk3, gen_0, l_chir), &
    g_w_quark(kk1, kk2, gen, kk3, gen_1, l_chir), &
    g_w_quark(kk1, kk2, gen, kk3, gen_2, l_chir)
  end do
  write (uunit, *)
  write (uunit, '("right-handed coupling matrix (rows <-> iso-up, &
   & columns <-> iso-down):")')
  write (uunit, *)
  do gen = gen_0, gen_2
   write (uunit, '(5X,F15.7,5X,F15.7,5X,F15.7)') &
    g_w_quark(kk1, kk2, gen, kk3, gen_0, r_chir), &
    g_w_quark(kk1, kk2, gen, kk3, gen_1, r_chir), &
    g_w_quark(kk1, kk2, gen, kk3, gen_2, r_chir)
  end do
  write (uunit, *)
  write (uunit, '("-------------------------------------------")')
  write (uunit, *)
 end do
 end do
 end do
 call errstack_pop
end subroutine print_gwff
subroutine print_gzff (unit)
integer, intent(in), optional :: unit
integer :: isospin, gen, kkb, kkf, uunit
character(len=slength), parameter :: fname="print_gzff"
 call errstack_push(fname)
 if (present (unit)) then
  uunit = unit
 else
  uunit = output_unit
 end if
 write (uunit, *)
 do kkb = l_mode, h_mode
  if (kkb == l_mode) then
   write (uunit, '("couplings to the light Z")')
   write (uunit, '("------------------------")')
  else
   write (uunit, '("couplings to the heavy Z")')
   write (uunit, '("------------------------")')
  end if
  do gen = gen_0, gen_2
  do isospin = iso_up, iso_down
   do kkf = l_mode, h_mode
    write (uunit, '("2x ",A,":")') trim(particle_names(retrieve_bcd_f( &
     kkf, ftype_l, gen, isospin)))
    write (uunit, '(" L: ",F15.7,5X,"R: ",F15.7)') &
     g_z_lep(kkb, kkf, gen, isospin, l_chir), &
     g_z_lep(kkb, kkf, gen, isospin, r_chir)
    write (uunit, *)
    write (uunit, '("2x ",A,":")') trim(particle_names(retrieve_bcd_f( &
     kkf, ftype_q, gen, isospin)))
    write (uunit, '(" L: ",F15.7,5X,"R: ",F15.7)') &
     g_z_quark(kkb, kkf, gen, isospin, l_chir), &
     g_z_quark(kkb, kkf, gen, isospin, r_chir)
    write (uunit, *)
   end do
   write (uunit, '(A," , ",A,":")') trim(particle_names(retrieve_bcd_f( &
    l_mode, ftype_l, gen, isospin))), trim(particle_names(retrieve_bcd_f( &
    h_mode, ftype_l, gen, isospin)))
   write (uunit, '(" L: ",F15.7,5X,"R: ",F15.7)') &
    g_z_lep(kkb, lh_mode, gen, isospin, l_chir), &
    g_z_lep(kkb, lh_mode, gen, isospin, r_chir)
   write (uunit, *)
   write (uunit, '(A," , ",A,":")') trim(particle_names(retrieve_bcd_f( &
    l_mode, ftype_q, gen, isospin))), trim(particle_names(retrieve_bcd_f( &
    h_mode, ftype_q, gen, isospin)))
   write (uunit, '(" L: ",F15.7,5X,"R: ",F15.7)') &
    g_z_quark(kkb, lh_mode, gen, isospin, l_chir), &
    g_z_quark(kkb, lh_mode, gen, isospin, r_chir)
   write (uunit, *)
  end do
  end do
 end do
 call errstack_pop
end subroutine print_gzff
subroutine print_gauge_coup (unit)
integer, intent(in), optional :: unit
integer :: kk1, kk2, uunit
character(len=slength), parameter :: fname="print_gauge_coup"
 call errstack_push(fname)
 if (present (unit)) then
  uunit = unit
 else
  uunit = output_unit
 end if
 write (uunit, *)
 write (uunit, '("couplings between gauge bosons")')
 write (uunit, '("------------------------------")')
 write (uunit, *)
 do kk1 = l_mode, lh_mode
 do kk2 = l_mode, h_mode
  if (kk1 == lh_mode) then
   write (uunit, '(A," , ",A," , ", A)', advance="no") trim(particle_names(w_bcd)), &
    trim(particle_names(hw_bcd)), trim(particle_names(retrieve_bcd_b(kk2, btype_z)))
  else
   write (uunit, '("2x",A," , ",A)', advance="no") trim(particle_names(retrieve_bcd_b(kk1, btype_w))), &
    trim(particle_names(retrieve_bcd_b(kk2, btype_z)))
  end if
  write (uunit, '(" : ",F15.7)') g_wwz(kk1, kk2)
  write (uunit, *)
 end do
 end do
 do i = 0, 4
  write (uunit, '(I1.1,"x ",A," , ",I1.1,"x ",A)', advance="no") (4-i), &
   trim(particle_names(w_bcd)), i, trim(particle_names(hw_bcd))
  write (uunit, '(" : ",F15.7)') g_wwww(i)
  write (uunit, *)
 end do
 do kk1 = l_mode, lh_mode
 do kk2 = l_mode, h_mode
  if (kk1 == lh_mode) then
   write (uunit, '(A," , ",A," , ",A," , photon")', advance="no") trim(particle_names(w_bcd)), &
    trim(particle_names(hw_bcd)), trim(particle_names(retrieve_bcd_b(kk2, btype_z)))
  else
   write (uunit, '("2x ",A," , ",A," , photon")', advance="no") &
    trim(particle_names(retrieve_bcd_b(kk1, btype_w))), &
    trim(particle_names(retrieve_bcd_b(kk2, btype_z)))
  end if
  write (uunit, '(" : ",F15.7)') g_wwza(kk1, kk2)
  write (uunit, *)
 end do
 end do
 do kk1 = l_mode, lh_mode
 do kk2 = l_mode, lh_mode
  if (kk1 == lh_mode) then
   write (uunit, '(A," , ",A," , ")', advance="no") trim(particle_names(w_bcd)), &
    trim(particle_names(hw_bcd))
  else
   write (uunit, '("2x ",A," , ")', advance="no") trim(particle_names(retrieve_bcd_b( &
    kk1, btype_w)))
  end if
  if (kk2 == lh_mode) then
   write (uunit, '(A," , ",A)', advance="no") trim(particle_names(z_bcd)), &
    trim(particle_names(hz_bcd))
  else
   write (uunit, '("2x ",A)', advance="no") trim(particle_names(retrieve_bcd_b( &
    kk2, btype_z)))
  end if
  write (uunit, '(" : ",F15.7)') g_wwzz(kk1, kk2)
  write (uunit, *)
 end do
 end do
 call errstack_pop
end subroutine print_gauge_coup
subroutine print_parameters (unit)
integer, intent(in), optional :: unit
integer :: uunit
 if (present (unit)) then
  uunit = unit
 else
  uunit = output_unit
 end if
 write (uunit, *)
 write (uunit, '("g0 = ",F10.7,5X,"g1 = ",F10.7,5X,"g2 = ",F10.7)') g0, g1, g2
 write (uunit, '(5X,"-> x = ",F10.7,5X,"t = ",F10.7,5X,"e = ",F10.7)') &
  x, t, e
 write (uunit, '(5X,"-> eps_L = ",F10.7)') eps_l
 write (uunit, '("v = ",F15.7,5X,"lambda =",F15.7)') sigma_vev, lambda
 write (uunit, *)
end subroutine print_parameters
end module threeshl
module tscript
use threeshl
use tdefs
implicit none
save
private
type tscript_tokenize_object
 character(len=slength) :: string
 integer :: next_char
 logical :: fully_tokenized=.false.
end type tscript_tokenize_object
public :: tscript_tokenize_object
real(kind=double), target :: zero=0._double, one=1._double
real(kind=double), target :: alpha_t, rwidths(0:63)
logical, public :: tscript_show_syntax=.false., tscript_error=.false.
interface tscript_calculate
 module procedure calculate
end interface tscript_calculate
public :: tscript_calculate
interface tscript_print_syntax
 module procedure print_syntax
end interface tscript_print_syntax
public :: tscript_print_syntax
interface tscript_tokenize
 module procedure tokenize
end interface tscript_tokenize
public :: tscript_tokenize
interface tscript_create_tobject
 module procedure create_tobject
end interface tscript_create_tobject
public :: tscript_create_tobject
interface tscript_decode_fspec
 module procedure decode_fspec
end interface tscript_decode_fspec
public :: tscript_decode_fspec
contains
subroutine calculate
integer :: i
 alpha_t = threeshl_eps_r(ftype_q, gen_2, iso_up)**4 * &
  threeshl_lambda**2 / 8._double / pi**2
 do i = 0, 63
  if (width_array(i) > 0.) then
   rwidths(i) = width_array(i) / mass_array(i)
  else
   rwidths(i) = 0.
  end if
 end do
end subroutine calculate
subroutine panic (msg)
character(len=*), intent(in), optional :: msg
 if (present (msg) .and. threeshl_print_errors) write (threeshl_errunit, *) msg
 if (tscript_show_syntax) then
  call print_syntax (quit=.false.)
 end if
 if (threeshl_quit_on_panic) call exit(1)
 tscript_error = .true.
end subroutine panic
subroutine print_syntax (quit)
logical, intent(in), optional :: quit
 write (threeshl_errunit, *)
 write (threeshl_errunit, *) "a 3shl function may be specified via one of the three possibilities:"
 write (threeshl_errunit, *)
 write (threeshl_errunit, *) "parameter( para )"
 write (threeshl_errunit, *) "lhcoupling( coup ), rhcoupling( coup )"
 write (threeshl_errunit, *) "mass( partspec )"
 write (threeshl_errunit, *) "width( partspec )"
 write (threeshl_errunit, *) "function( funspec )"
 write (threeshl_errunit, *) "rwidth( partspec )"
 write (threeshl_errunit, *)
 write (threeshl_errunit, *) "Model parameters are specified via ""parameter"", where ""para"" has to be one of:"
 write (threeshl_errunit, *) "g0, g1, g2, x, t, lambda, eps_l, eps_r (which has to be followed by a"
 write (threeshl_errunit, *) "fermion identifier)"
 write (threeshl_errunit, *)
 write (threeshl_errunit, *) "Vertex factors are specified via ""lhcoupling"" or ""rhcoupling""; the chirality"
 write (threeshl_errunit, *) "is ignored for vertices coupling only gauge bosons. The parameter ""coup"" is a"
 write (threeshl_errunit, *) "list of the particles meeting at the vertex seperated by spaces. Only couplings"
 write (threeshl_errunit, *) "between 3 or 4 particles with either no or 2 fermions are accepted as valid."
 write (threeshl_errunit, *)
 write (threeshl_errunit, *) "Masses are specified via ""mass"" with ""partspec"" being a particle identifier"
 write (threeshl_errunit, *)
 write (threeshl_errunit, *) "Widths are specified similarly to masses via ""mass"" and ""partspec""."
 write (threeshl_errunit, *)
 write (threeshl_errunit, *) """function"" is used to calculate functions of the threeshl parameters. At the"
 write (threeshl_errunit, *) "moment, these are:"
 write (threeshl_errunit, *) """alpha_t""   : electroweak precision observable as estimated by Chivukula et al."
 write (threeshl_errunit, *)
 write (threeshl_errunit, *) "Possible light fermions are e, nue, mu, numu, tau, nutau, u, d, c, s, t, b ; the"
 write (threeshl_errunit, *) "light gauge bosons are idenitified by w, z, a. The heavy KK modes are select by"
 write (threeshl_errunit, *) "the same specifiers with the first letter capitalized."
 write (threeshl_errunit, *)
 if (present (quit)) then
  if (quit) call exit (1)
 else
  call exit (1)
 end if
end subroutine print_syntax
function tokenize (input) result(token)
type(tscript_tokenize_object), intent(inout) :: input
character(len=slength) :: token
integer :: i
character(len=1) :: current_char
logical :: wait_for_token, token_complete
 if (input%next_char > len_trim(input%string)) then
  write (threeshl_errunit, *) "INTERNAL ERROR: tokenize called on already fully parsed or empty string"
  call exit(1)
 end if
 token = ""
 i = 1
 input%fully_tokenized = .false.
 wait_for_token = .true.
 token_complete = .false.
 do while (.not. token_complete)
  current_char = input%string(input%next_char:input%next_char)
  select case(current_char)
   case(" ")
    if (.not. wait_for_token) token_complete = .true.
    input%next_char = input%next_char + 1
   case("(",")","#")
    if (wait_for_token) then
     token(i:i) = current_char
     input%next_char = input%next_char + 1
    end if
    token_complete = .true.
   case default
    token(i:i) = current_char
    input%next_char = input%next_char +1
    wait_for_token = .false.
    i = i + 1
  end select
  if (input%next_char > len_trim(input%string)) then
   input%fully_tokenized = .true.
   token_complete = .true.
  end if
 end do
end function tokenize
function decode_parameter (input, full_parse) result(ptr)
type(tscript_tokenize_object), intent(inout) :: input
logical, optional, intent(in) :: full_parse
real(kind=double), pointer :: ptr
integer, parameter :: state_start=1, state_get_para=2, state_end=3, state_stop=4, &
 state_get_ferm=5
integer :: state
character(len=slength) :: token, emesg
logical :: parse
 nullify (ptr)
 emesg = "ERROR: malformed parameter specification in """ // trim(input%string) // """"
 state = state_start
 parse = .true.
 do while((.not. input%fully_tokenized) .and. parse)
  token = tokenize(input)
  select case(trim (token))
   case("(")
    if (state == state_start) then
     state = state_get_para
    else
     call panic (msg=trim (emesg))
     return
    end if
   case(")")
    if (state == state_end) then
     state = state_stop
     if (present(full_parse)) parse = full_parse
    else
     call panic (msg=trim (emesg))
     return
    end if
   case("eps_r")
    if (state == state_get_para) then
     state = state_get_ferm
    else
     call panic (msg=trim (emesg))
     return
    end if
   case default
    select case(state)
     case(state_get_para)
      state = state_end
      ptr => dec_token(token)
     case(state_get_ferm)
      state = state_end
      ptr => dec_ferm(token)
     case default
      call panic (msg=trim (emesg))
      return
    end select
  end select
 end do
 if (state .ne. state_stop) call panic (msg=trim (emesg))
contains
function dec_token(tok) result(res)
character(len=slength), intent(in) :: tok
real(kind=double), pointer :: res
 nullify (res)
 select case(trim(tok))
  case("sigma_vev")
   res => threeshl_sigma_vev
  case("g0")
   res => threeshl_g0
  case("g1")
   res => threeshl_g1
  case("g2")
   res => threeshl_g2
  case("lambda")
   res => threeshl_lambda
  case("eps_l")
   res => threeshl_eps_l
  case("x")
   res => threeshl_x
  case("t")
   res => threeshl_t
  case default
   call panic (msg=trim (emesg))
 end select
end function dec_token
function dec_ferm(tok) result(res)
character(len=slength), intent(in) :: tok
real(kind=double), pointer :: res
integer :: kk, ftype, gen, isospin
 nullify (res)
 call decode_fermion_name(tok, kk, ftype, gen, isospin)
 if (ftype >= 0) res => threeshl_eps_r(ftype, gen, isospin)
end function dec_ferm
end function decode_parameter
function decode_genericpart (input, nm, full_parse) result(bcd)
type(tscript_tokenize_object), intent(inout) :: input
character(len=*), intent(in) :: nm
logical, intent(in), optional :: full_parse
integer :: bcd
integer, parameter :: state_start=1, state_get_part=2, state_end=3, state_stop=4
integer :: state
character(len=slength) :: token, emesg
logical :: parse
 emesg = "ERROR: malformed " // nm // " specification in """ // trim(input%string) // """"
 parse = .true.
 state = state_start
 do while ((.not. input%fully_tokenized) .and. parse)
  token = tokenize(input)
  select case(trim(token))
   case("(")
    if (state == state_start) then
     state = state_get_part
    else
     call panic (msg=trim (emesg))
     return
    end if
   case(")")
    if (state == state_end) then
     state = state_stop
     if (present(full_parse)) parse = full_parse
    else
     call panic (msg=trim (emesg))
     return
    end if
   case default
    if (state == state_get_part) then
     bcd = bcd_from_name(token)
     state = state_end
    else
     call panic (msg=trim (emesg))
    end if
  end select
 end do
 if (state .ne. state_stop) call panic (msg=trim (emesg))
end function decode_genericpart
function decode_function (input, full_parse) result(ptr)
type(tscript_tokenize_object), intent(inout) :: input
logical, intent(in), optional :: full_parse
real(kind=double), pointer :: ptr
integer, parameter :: state_start=1, state_get_fun=2, state_end=3, state_stop=4
integer :: state
character(len=slength) :: token, emesg
logical :: parse
 emesg = "ERROR: malformed function specification in """ // trim(input%string) // """"
 nullify (ptr)
 parse = .true.
 state = state_start
 do while ((.not. input%fully_tokenized) .and. parse)
  token = tokenize(input)
  select case(trim(token))
   case("(")
    if (state == state_start) then
     state = state_get_fun
    else
     call panic (msg=trim (emesg))
     return
    end if
   case(")")
    if (state == state_end) then
     state = state_stop
     if (present(full_parse)) parse = full_parse
    else
     call panic (msg=trim (emesg))
     return
    end if
   case("alpha_t")
    if (state == state_get_fun) then
     ptr => alpha_t
     state = state_end
    else
     call panic (msg=trim (emesg))
     return
    end if
   case default
    call panic (msg=trim (emesg))
    return
  end select
 end do
 if (state .ne. state_stop) call panic (msg=trim(emesg))
end function decode_function
function decode_coupling (input, chir, full_parse) result(ptr)
type(tscript_tokenize_object), intent(inout) :: input
integer, intent(in) :: chir
logical, intent(in), optional :: full_parse
real(kind=double), pointer :: ptr
type fermion
 integer :: kk, ftype, gen, isospin
end type fermion
type(fermion), dimension(2) :: flist
type(fermion) :: tmp
integer :: nf, nw, nz, na, nhw, nhz, tmpi, btype, bkk
integer, parameter :: state_start=1, state_get_pname=2, state_stop=4
integer :: state
character(len=slength) :: token, emesg
logical :: parse
 nullify (ptr)
 emesg = "ERROR: malformed coupling specification in """ // trim(input%string) // """"
 if ((chir < l_chir) .or. (chir > r_chir)) call internal_error
 state = state_start
 parse = .true.
 nf = 0; nw = 0; nz = 0; na = 0; nhw = 0; nhz = 0
 do while ((.not. input%fully_tokenized) .and. parse)
  token = tokenize(input)
  select case(token)
   case ("(")
    if (state == state_start) then
     state = state_get_pname
    else
     call panic (msg=trim (emesg))
     return
    end if
   case(")")
    if (state == state_get_pname) then
     state = state_stop
     if (present(full_parse)) parse = full_parse
    else
     call panic (msg=trim (emesg))
     return
    end if
   case default
    if (state == state_stop) then
     call panic (msg=trim (emesg))
     return
    end if
    select case(token)
     case("w", "W", "z", "Z", "a")
      if (nb() .ge. 4) then
       call panic (msg="ERROR: too many gauge bosons in coupling specification in """ &
        // trim(input%string) // """")
       return
      else
       call decode_boson_name(token, bkk, btype)
       if (bkk < 0) return
       select case(btype)
        case(btype_w)
         nw = nw + 1
         if (bkk == h_mode) nhw = nhw + 1
        case(btype_z)
         nz = nz + 1
         if (bkk == h_mode) nhz = nhz + 1
        case(btype_a)
         na = na + 1
        case default
         call internal_error
       end select
      end if
     case default
      if (nf .ge. 2) then
       call panic (msg="ERROR: too many fermions in coupling specification&
        & (our illegal particle name)")
       return
      else
       nf = nf + 1
       call decode_fermion_name(token, flist(nf)%kk, flist(nf)%ftype,&
        flist(nf)%gen, flist(nf)%isospin)
        if (flist(nf)%kk < 0) return
      end if
    end select
  end select
 end do
 if ((state .ne. state_stop) .or. &
  (.not. ((nf == 0) .or. (nf == 2))) .or. &
  (nf + nb() > 4) .or. (nf + nb() < 3)) then
   call panic (msg=trim (emesg))
   return
 end if
 if ((nf == 2) .and. (nb() == 1) .and. (flist(1)%ftype == flist(2)%ftype)) then
  if ((flist(1)%isospin == iso_down) .and. (flist(2)%isospin == iso_up)) then
   tmp = flist(1)
   flist(1) = flist(2)
   flist(2) = tmp
  end if
  select case(btype)
   case(btype_w)
    if ((flist(1)%isospin == iso_up) .and. (flist(2)%isospin == iso_down)) then
     select case(flist(1)%ftype)
      case(ftype_l)
       ptr => g_w_lep(bkk, flist(1)%kk, flist(1)%gen, &
        flist(2)%kk, flist(2)%gen, chir)
      case(ftype_q)
       ptr => g_w_quark(bkk, flist(1)%kk, flist(1)%gen, &
        flist(2)%kk, flist(2)%gen, chir)
      case default
       call internal_error
     end select
    else
     ptr => zero
    end if
   case(btype_z)
    if ((flist(1)%gen == flist(2)%gen) .and. &
     (flist(1)%isospin == flist(2)%isospin)) then
     if (flist(1)%kk ==flist(2)%kk) then
      tmpi = flist(1)%kk
     else
      tmpi = lh_mode
     end if
     select case(flist(1)%ftype)
      case(ftype_l)
       ptr => g_z_lep(bkk, tmpi, flist(1)%gen, flist(1)%isospin, chir)
      case(ftype_q)
       ptr => g_z_quark(bkk, tmpi, flist(1)%gen, flist(1)%isospin, chir)
      case default
       call internal_error
     end select
    else
     ptr => zero
    end if
   case(btype_a)
    if ((flist(1)%kk == flist(2)%kk) .and. (flist(1)%gen == flist(2)%gen)&
     .and. (flist(1)%isospin == flist(2)%isospin)) then
     ptr => threeshl_e
    else
     ptr => zero
    end if
   case default
    call internal_error
  end select
 else if ((nb() == 3) .and. (nw == 2.)) then
  if ((nz == 1) .and. (na == 0)) then
   ptr => g_wwz(nhx_to_modex(nhw), l_mode + nhz)
  else if ((nz == 0) .and. (na == 1)) then
   if ((nhw == 0) .or. (nhw == 2)) then
    ptr => threeshl_e
   else
    ptr => zero
   end if
  else
   call internal_error
  end if
 else if (nw == 4) then
  ptr => g_wwww(nhw)
 else if ((nw == 2) .and. (nz == 2)) then
  ptr => g_wwzz(nhx_to_modex(nhw), nhx_to_modex(nhz))
 else if ((nw == 2) .and. (nz == 1) .and. (na == 1)) then
  ptr => g_wwza(nhx_to_modex(nhw), l_mode + nhz)
 else
  ptr => zero
 end if
contains
function nb () result (res)
integer :: res
 res = na + nw + nz
end function nb
subroutine internal_error
 write (threeshl_errunit, *) "internal error in decode_coupling; very bad"
 write (threeshl_errunit, *) "string parsed: """, trim(input%string), """"
 call exit(1)
end subroutine internal_error
function nhx_to_modex (n) result(res)
integer, intent(in) :: n
integer :: res
 select case(n)
  case(0); res = l_mode
  case(1); res = lh_mode
  case(2); res = h_mode
  case default; call internal_error
 end select
end function nhx_to_modex
end function decode_coupling
subroutine decode_fermion_name(str, kk, ftype, gen, isospin)
character(len=slength), intent(in) :: str
integer, intent(out) :: kk, ftype, gen, isospin
 if (lgt(str(1:1), "Z")) then
  kk = l_mode
 else
  kk = h_mode
 end if
 select case(trim(str))
  case("nue", "Nue")
   ftype = ftype_l; gen = gen_0; isospin = iso_up
  case("e", "E")
   ftype = ftype_l; gen = gen_0; isospin = iso_down
  case("numu", "Numu")
   ftype = ftype_l; gen = gen_1; isospin = iso_up
  case("mu", "Mu")
   ftype = ftype_l; gen = gen_1; isospin = iso_down
  case("nutau", "Nutau")
   ftype = ftype_l; gen = gen_2; isospin = iso_up
  case("tau", "Tau")
   ftype = ftype_l; gen = gen_2; isospin = iso_down
  case("u", "U")
   ftype = ftype_q; gen = gen_0; isospin = iso_up
  case("d", "D")
   ftype = ftype_q; gen = gen_0; isospin = iso_down
  case("c", "C")
   ftype = ftype_q; gen = gen_1; isospin = iso_up
  case("s", "S")
   ftype = ftype_q; gen = gen_1; isospin = iso_down
  case("t", "T")
   ftype = ftype_q; gen = gen_2; isospin = iso_up
  case("b", "B")
   ftype = ftype_q; gen = gen_2; isospin = iso_down
  case default
   call panic (msg="ERROR: invalid fermion identifier """ // trim(str) // """")
   kk = -1; ftype = -1; gen = -1; isospin = -1
 end select
end subroutine decode_fermion_name
subroutine decode_boson_name (str, kk, btype)
character(len=slength), intent(in) :: str
integer, intent(out) :: kk, btype
 select case(trim(str))
  case("w")
   kk = l_mode; btype = btype_w
  case("W")
   kk = h_mode; btype = btype_w
  case("z")
   kk = l_mode; btype = btype_z
  case("Z")
   kk = h_mode; btype = btype_z
  case("a")
   kk = l_mode; btype = btype_a
  case default
   call panic (msg="ERROR: invalid gauge boson identifier """ // trim(str) // """")
   kk = -1; btype = -1
 end select
end subroutine decode_boson_name
function bcd_from_name (str) result(bcd)
character(len=slength), intent(in) :: str
integer :: bcd, kk, ftype, btype, gen, isospin
 select case(trim(str))
  case("w", "W", "z", "Z", "a")
   call decode_boson_name(str, kk, btype)
   if (kk < 0) then
    bcd = -1
   else
    bcd = threeshl_retrieve_bcd_b(kk, btype)
   end if
  case default
   call decode_fermion_name(str, kk, ftype, gen, isospin)
   if (kk < 0) then
    bcd = -1
   else
    bcd = threeshl_retrieve_bcd_f(kk, ftype, gen, isospin)
   end if
 end select
end function bcd_from_name
function create_tobject (string) result(tobject)
character(len=slength), intent(in) :: string
type(tscript_tokenize_object) :: tobject
 tobject%string = string
 tobject%next_char = 1
 if (len_trim(string) == 0) then
  tobject%fully_tokenized = .true.
 else
  tobject%fully_tokenized = .false.
 end if
end function create_tobject
function decode_fspec (input, full_parse, allow_comment) result(ptr)
type(tscript_tokenize_object), intent(inout) :: input
logical, intent(in), optional :: full_parse, allow_comment
real(kind=double), pointer :: ptr
character(len=slength) :: token
logical :: comment_allowed
 nullify (ptr)
 if (present(allow_comment)) then
  comment_allowed = allow_comment
 else
  comment_allowed = .false.
 end if
 if (input%fully_tokenized) then
  print *, "ERROR: empty 3shl function definition!"
  call print_syntax
 end if
 token = tokenize(input)
 select case(trim(token))
  case("parameter")
   ptr => decode_parameter(input, full_parse=full_parse)
  case("mass")
   ptr => mass_array(decode_genericpart (input, "mass", full_parse=full_parse))
  case("width")
   ptr => width_array(decode_genericpart (input, "width", full_parse=full_parse))
  case("rwidth")
   ptr => rwidths(decode_genericpart (input, "rwidth", full_parse=full_parse))
  case("lhcoupling")
   ptr => decode_coupling(input, l_chir, full_parse=full_parse)
  case("rhcoupling")
   ptr => decode_coupling(input, r_chir, full_parse=full_parse)
  case("function")
   ptr => decode_function(input, full_parse=full_parse)
  case("#")
   if (comment_allowed) then
    nullify(ptr)
   else
    call panic (msg="ERROR: comment not allowed at this point in """ // &
     trim(input%string) // """")
    return
   end if
  case default
   call panic (msg="ERROR: illegal 3shl function specification in """ // &
    trim(input%string) //"""")
   return
 end select
end function decode_fspec
end module tscript
module tglue
use tdefs
use threeshl
use tscript
implicit none
save
private
complex(kind=double), public :: g_a_lep, g_a_quark (iso_up:iso_down), g_aaww, &
 ig_aww, ig_wwz (l_mode:lh_mode, l_mode:h_mode)
complex(kind=double), public :: &
 g_w_lep_va (1:2, l_mode:h_mode, l_mode:h_mode, gen_0:gen_2, l_mode:h_mode, gen_0:gen_2), &
 g_w_quark_va (1:2, l_mode:h_mode, l_mode:h_mode, gen_0:gen_2, l_mode:h_mode, gen_0:gen_2), &
 g_z_lep_va (1:2, l_mode:h_mode, l_mode:lh_mode, gen_0:gen_2, iso_up:iso_down), &
 g_z_quark_va (1:2, l_mode:h_mode, l_mode:lh_mode, gen_0:gen_2, iso_up:iso_down)
real(kind=double), public :: g_s=1.218_double
complex(kind=double), public :: ig_s_norm, g_s_norm, g_s_norm2
public :: tglue_set_alphas
public :: tglue_init
public :: tglue_init_ward
public :: tglue_finalize
public :: tglue_mcid_to_bcd
contains
subroutine init_couplings
integer :: kk1, kk2, kk3, gen1, gen2, iso
 g_a_lep = (-1._double) * threeshl_e
 g_a_quark (iso_up) = threeshl_e * 2._double / 3._double
 g_a_quark (iso_down) = threeshl_e / (-3._double)
 ig_aww = (0, 1._double) * threeshl_e
 g_aaww = threeshl_e**2
 do kk1 = l_mode, lh_mode
 do kk2 = l_mode, h_mode
  ig_wwz(kk1, kk2) = (0., 1._double) * g_wwz(kk1, kk2)
 end do; end do
 do kk1 = l_mode, h_mode
 do kk2 = l_mode, h_mode
 do kk3 = l_mode, h_mode
 do gen1 = gen_0, gen_2
 do gen2 = gen_0, gen_2
  g_w_lep_va(:, kk1, kk2, gen1, kk3, gen2) = (/ &
   g_w_lep(kk1, kk2, gen1, kk3, gen2, r_chir) + g_w_lep(kk1, kk2, gen1, kk3, gen2, l_chir), &
   g_w_lep(kk1, kk2, gen1, kk3, gen2, l_chir) - g_w_lep(kk1, kk2, gen1, kk3, gen2, r_chir) /) &
   / sqrt(8._double)
  g_w_quark_va(:, kk1, kk2, gen1, kk3, gen2) = (/ &
   g_w_quark(kk1, kk2, gen1, kk3, gen2, r_chir) + g_w_quark(kk1, kk2, gen1, kk3, gen2, l_chir), &
   g_w_quark(kk1, kk2, gen1, kk3, gen2, l_chir) - g_w_quark(kk1, kk2, gen1, kk3, gen2, r_chir) /) &
   / sqrt(8._double)
 end do; end do; end do; end do; end do
 do kk1 = l_mode, h_mode
 do kk2 = l_mode, lh_mode
 do gen1 = gen_0, gen_2
 do iso = iso_up, iso_down
  g_z_lep_va(:, kk1, kk2, gen1, iso) = (/ &
   g_z_lep(kk1, kk2, gen1, iso, r_chir) + g_z_lep(kk1, kk2, gen1, iso, l_chir), &
   g_z_lep(kk1, kk2, gen1, iso, l_chir) - g_z_lep(kk1, kk2, gen1, iso, r_chir) /) / 2._double
  g_z_quark_va(:, kk1, kk2, gen1, iso) = (/ &
   g_z_quark(kk1, kk2, gen1, iso, r_chir) + g_z_quark(kk1, kk2, gen1, iso, l_chir), &
   g_z_quark(kk1, kk2, gen1, iso, l_chir) - g_z_quark(kk1, kk2, gen1, iso, r_chir) /) / 2._double
 end do; end do; end do; end do
 g_s_norm = g_s / sqrt (2._double)
 g_s_norm2 = g_s_norm**2
 ig_s_norm = (0., 1._double) * g_s_norm
end subroutine init_couplings
subroutine tglue_set_alphas (as)
real(kind=double), intent(in) :: as
 g_s = sqrt (4._double * pi * as)
 g_s_norm = g_s / sqrt (2._double)
 g_s_norm2 = g_s_norm**2
 ig_s_norm = (0., 1._double) * g_s_norm
end subroutine tglue_set_alphas
subroutine tglue_init (mhw, bmass, el)
real(kind=double), intent(in) :: mhw, bmass
real(kind=double), intent(in), optional :: el
 call threeshl_init
 call threeshl_pdg_init_wgap_bmass (mhw, bmass, el)
 if (.not. threeshl_error) then
  call threeshl_calculate_widths
  call init_couplings
  call tscript_calculate
 end if
end subroutine tglue_init
subroutine tglue_init_ward (x, ct, ph)
real(kind=double), intent(in) :: x, ct, ph
 call threeshl_init
 call threeshl_init_ward (x, ct, ph)
 width_array = 0._double
 call init_couplings
end subroutine tglue_init_ward
subroutine tglue_finalize
 call threeshl_finalize
end subroutine tglue_finalize
function tglue_mcid_to_bcd (mcid) result(bcd)
integer, intent(in) :: mcid
integer :: bcd, t
 if (abs (mcid) .ge. 9900) then
  t = mcid - sign (9900, mcid)
 else
  t = mcid
 end if
 select case (abs(t))
  case(1); bcd = d_bcd
  case(2); bcd = u_bcd
  case(3); bcd = s_bcd
  case(4); bcd = c_bcd
  case(5); bcd = b_bcd
  case(6); bcd = t_bcd
  case(11); bcd = e_bcd
  case(12); bcd = nue_bcd
  case(13); bcd = mu_bcd
  case(14); bcd = numu_bcd
  case(15); bcd = tau_bcd
  case(16); bcd = nutau_bcd
  case(22); bcd = a_bcd
  case(23); bcd = z_bcd
  case(24); bcd = w_bcd
  case default
   print *, "Internal error in tglue_mcid_to_bcd: invalid MC id ", t
   call exit(1)
 end select
 if (t /= mcid) then; if (bcd == a_bcd) then
  print *, "Internal error in tglue_mcid_to_bcd:&
   &invalid MC id (there is no heavy photon)."
  call exit(1)
 else
  bcd = bcd + 1
 end if; end if
end function tglue_mcid_to_bcd
end module
