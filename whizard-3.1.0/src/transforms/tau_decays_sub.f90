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

submodule (tau_decays) tau_decays_s

  use io_units
  use format_utils, only: write_separator

  implicit none

contains

  module subroutine evt_tau_decays_write_name (evt, unit)
    class(evt_tau_decays_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Event transform: tau decays"
  end subroutine evt_tau_decays_write_name

  module subroutine evt_tau_decays_write &
       (evt, unit, verbose, more_verbose, testflag)
    class(evt_tau_decays_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose, more_verbose, testflag
    integer :: u
    u = given_output_unit (unit)
    call write_separator (u, 2)
    call evt%write_name (u)
    call write_separator (u)
    call evt%base_write (u, testflag = testflag, show_set = .false.)
    if (evt%particle_set_exists)  &
         call evt%particle_set%write &
         (u, summary = .true., compressed = .true., testflag = testflag)
    call write_separator (u)
  end subroutine evt_tau_decays_write

  module subroutine evt_tau_decays_generate_weighted (evt, probability)
    class(evt_tau_decays_t), intent(inout) :: evt
    real(default), intent(inout) :: probability
    logical :: valid
    evt%particle_set = evt%previous%particle_set
    !!! To be checked or expanded
    probability = 1
    valid = .true.
    evt%particle_set_exists = valid
  end subroutine evt_tau_decays_generate_weighted

  module subroutine evt_tau_decays_make_particle_set &
       (evt, factorization_mode, keep_correlations, r)
    class(evt_tau_decays_t), intent(inout) :: evt
    integer, intent(in) :: factorization_mode
    logical, intent(in) :: keep_correlations
    real(default), dimension(:), intent(in), optional :: r
    logical :: valid
    !!! to be checked and expanded
    valid = .true.
    evt%particle_set_exists = evt%particle_set_exists .and. valid
  end subroutine evt_tau_decays_make_particle_set

  module subroutine evt_tau_decays_prepare_new_event (evt, i_mci, i_term)
    class(evt_tau_decays_t), intent(inout) :: evt
    integer, intent(in) :: i_mci, i_term
    call evt%reset ()
  end subroutine evt_tau_decays_prepare_new_event


end submodule tau_decays_s

