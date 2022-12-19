!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     Christian Speckner <christian.speckner@physik.uni-wuerzburg.de>
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

subroutine init_external_parameters (par) bind (C)
use iso_c_binding
use kinds
use io_units
use threeshl
use tglue
use diagnostics
implicit none

real(c_default_float), dimension(*), intent(inout) :: par
real(default) :: alphas, mhw, mbulk, eps_l, nideloc, dump, nlow
integer :: u

	me_pdg =    par(1)
	mmu_pdg =   par(2)
	mtau_pdg =  par(3)
	muq_pdg =   par(4)
	mdq_pdg =   par(5)
	mcq_pdg =   par(6)
	msq_pdg =   par(7)
	mtq_pdg =   par(8)
	mbq_pdg =   par(9)
	mw_pdg =    par(10)
	mz_pdg =    par(11)
	wt_pdg =    par(12)
	ww_pdg =    par(13)
	wz_pdg =    par(14)
	e_pdg =     par(15)
	alphas =    par(16)
	mhw =       par(17)
	mbulk =     par(18)
	eps_l =     par(19)
	nideloc=    par(20)
	dump =      par(21)
	nlow =      par(22)

	threeshl_error = .false.
	threeshl_quit_on_panic = .false.
	threeshl_print_errors = .false.
	call tglue_set_alphas (alphas)
	if (nideloc > 0.) then
		call tglue_init (mhw, mbulk, eps_l)
	else
		call tglue_init (mhw, mbulk)
	end if
	if (threeshl_error) then
		call msg_warning ( &
			"3SHL initialized with invalid parameters - now in undefined state!")
	else
		call msg_message ("3SHL initialized sucessfully.")
	end if
	if (dump > 0.) then
		u = free_unit ()
		open (unit=u, file="threeshl.whizard.dump", action="write", status="replace")
		call threeshl_print_parameters (unit=u)
		call threeshl_print_particles (unit=u)
		call threeshl_print_gwff (unit=u)
		call threeshl_print_gzff (unit=u)
		call threeshl_print_gauge_coup (unit=u)
		close (u)
	end if

	par(23) = mass_array(hz_bcd)
	par(24) = mass_array(hu_bcd)
	par(25) = mass_array(hd_bcd)
	par(26) = mass_array(hc_bcd)
	par(27) = mass_array(hs_bcd)
	par(28) = mass_array(ht_bcd)
	par(29) = mass_array(hb_bcd)
	par(30) = mass_array(he_bcd)
	par(31) = mass_array(hnue_bcd)
	par(32) = mass_array(hmu_bcd)
	par(33) = mass_array(hnumu_bcd)
	par(34) = mass_array(htau_bcd)
	par(35) = mass_array(hnutau_bcd)
	par(36) = width_array(hw_bcd)
	par(37) = width_array(hz_bcd)
	par(38) = width_array(hu_bcd)
	par(39) = width_array(hd_bcd)
	par(40) = width_array(hc_bcd)
	par(41) = width_array(hs_bcd)
	par(42) = width_array(ht_bcd)
	par(43) = width_array(hb_bcd)
	par(44) = width_array(he_bcd)
	par(45) = width_array(hnue_bcd)
	par(46) = width_array(hmu_bcd)
	par(47) = width_array(hnumu_bcd)
	par(48) = width_array(htau_bcd)
	par(49) = width_array(hnutau_bcd)

end subroutine init_external_parameters
