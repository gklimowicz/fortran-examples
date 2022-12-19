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

module dispatch_mci

  use iso_varying_string, string_t => varying_string
  use variables
  use mci_base

  implicit none
  private

  public :: dispatch_mci_setup
  public :: setup_grid_path

  interface
    module subroutine dispatch_mci_setup (mci, var_list, process_id, is_nlo)
      class(mci_t), allocatable, intent(out) :: mci
      type(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: process_id
      logical, intent(in), optional :: is_nlo
    end subroutine dispatch_mci_setup
    module subroutine setup_grid_path (grid_path)
      type(string_t), intent(in) :: grid_path
    end subroutine setup_grid_path
  end interface

end module dispatch_mci
