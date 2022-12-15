
! Copyright (C) 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin and Lars Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Stub routines for Wannier90 library.

subroutine wannier_setup(seed__name,mp_grid_loc,num_kpts_loc, &
 real_lattice_loc,recip_lattice_loc,kpt_latt_loc,num_bands_tot, &
 num_atoms_loc,atom_symbols_loc,atoms_cart_loc, gamma_only_loc,spinors_loc, &
 nntot_loc,nnlist_loc,nncell_loc,num_bands_loc,num_wann_loc, &
 proj_site_loc,proj_l_loc,proj_m_loc,proj_radial_loc,proj_z_loc, &
 proj_x_loc,proj_zona_loc,exclude_bands_loc,proj_s_loc,proj_s_qaxis_loc)
implicit none
! arguments
character(*) seed__name
integer mp_grid_loc(3)
integer num_kpts_loc
real(8) real_lattice_loc(3,3)
real(8) recip_lattice_loc(3,3)
real(8) kpt_latt_loc(3,num_kpts_loc)
integer num_bands_tot
integer num_atoms_loc
character(*) atom_symbols_loc(num_atoms_loc)
real(8) atoms_cart_loc(3,num_atoms_loc)
logical gamma_only_loc
logical spinors_loc
integer nntot_loc
integer nnlist_loc(num_kpts_loc,*)
integer nncell_loc(3,num_kpts_loc,*)
integer num_bands_loc
integer num_wann_loc
real(8) proj_site_loc(3,num_bands_tot)
integer proj_l_loc(num_bands_tot)
integer proj_m_loc(num_bands_tot)
integer proj_radial_loc(num_bands_tot)
real(8) proj_z_loc(3,num_bands_tot)
real(8) proj_x_loc(3,num_bands_tot)
real(8) proj_zona_loc(num_bands_tot)
integer exclude_bands_loc(num_bands_tot)
integer, optional :: proj_s_loc(num_bands_tot)
real(8), optional :: proj_s_qaxis_loc(3,num_bands_tot)
write(*,*)
write(*,'("Error(wannier_setup): libwannier not or improperly installed")')
write(*,*)
error stop
end subroutine
!EOC

