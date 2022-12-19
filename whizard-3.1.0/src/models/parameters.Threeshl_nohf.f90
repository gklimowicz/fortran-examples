! parameters.Threeshl.f90 
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
module parameters_threeshl_nohf
use threeshl
use tglue
use kinds, only:default
use io_units, only: free_unit
implicit none
public

private :: default, free_unit

contains

subroutine import_from_whizard (par, scheme)
  real(default), dimension(49), intent(in) :: par
  integer, intent(in) :: scheme
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
	nideloc =   par(20)
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
	if (dump > 0.) then
		u = free_unit ()
		open (unit=u, file="threeshl.omega.dump", action="write", status="replace")
		call threeshl_print_parameters (unit=u)
		call threeshl_print_particles (unit=u)
		call threeshl_print_gwff (unit=u)
		call threeshl_print_gzff (unit=u)
		call threeshl_print_gauge_coup (unit=u)
		close (u)
	end if

end subroutine import_from_whizard

subroutine model_update_alpha_s (alpha_s)
real(default), intent(in) :: alpha_s
	call tglue_set_alphas (alpha_s)
end subroutine model_update_alpha_s

end module parameters_threeshl_nohf
