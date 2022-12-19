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

module shower_base_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use format_utils, only: write_separator
  use variables
  use shower_base

  implicit none
  private

  public :: shower_base_1

contains

  subroutine shower_base_1 (u)
    integer, intent(in) :: u
    type(var_list_t) :: var_list
    type(shower_settings_t) :: shower_settings

    write (u, "(A)")  "* Test output: shower_base_1"
    write (u, "(A)")  "*   Purpose: setting ISR/FSR shower"
    write (u, "(A)")

    write (u, "(A)")  "* Default settings"
    write (u, "(A)")

    call var_list%init_defaults (0)
    call var_list%set_log (var_str ("?alphas_is_fixed"), &
         .true., is_known = .true.)
    call shower_settings%init (var_list)
    call write_separator (u)
    call shower_settings%write (u)
    call write_separator (u)

    write (u, "(A)")
    write (u, "(A)")  "* Switch on ISR/FSR showers, hadronization"
    write (u, "(A)")  "      and MLM matching"
    write (u, "(A)")

    call var_list%set_string (var_str ("$shower_method"), &
         var_str ("PYTHIA6"), is_known = .true.)
    call var_list%set_log (var_str ("?ps_fsr_active"), &
         .true., is_known = .true.)
    call var_list%set_log (var_str ("?ps_isr_active"), &
         .true., is_known = .true.)
    call var_list%set_log (var_str ("?hadronization_active"), &
         .true., is_known = .true.)
    call var_list%set_log (var_str ("?mlm_matching"), &
         .true., is_known = .true.)
    call var_list%set_int &
         (var_str ("ps_max_n_flavors"), 4, is_known = .true.)
    call var_list%set_real &
         (var_str ("ps_isr_z_cutoff"), 0.1234_default, &
          is_known=.true.)
    call var_list%set_real (&
         var_str ("mlm_etamax"), 3.456_default, is_known=.true.)
    call var_list%set_string (&
         var_str ("$ps_PYTHIA_PYGIVE"), var_str ("abcdefgh"), is_known=.true.)
    call shower_settings%init (var_list)
    call write_separator (u)
    call shower_settings%write (u)
    call write_separator (u)

    call var_list%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: shower_base_1"

  end subroutine shower_base_1


end module shower_base_uti
