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

module dispatch_phase_space

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use variables, only: var_list_t
  use os_interface, only: os_data_t

  use sf_mappings, only: sf_channel_t
  use beam_structures, only: beam_structure_t
  use dispatch_beams, only: sf_prop_t, strfun_mode

  use mappings
  use phs_forests, only: phs_parameters_t
  use phs_base

  implicit none
  private

  public :: dispatch_phs
  public :: dispatch_sf_channels

  interface
    module subroutine dispatch_phs (phs, var_list, os_data, process_id, &
           mapping_defaults, phs_par, phs_method_in)
      class(phs_config_t), allocatable, intent(inout) :: phs
      type(var_list_t), intent(in) :: var_list
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in) :: process_id
      type(mapping_defaults_t), intent(in), optional :: mapping_defaults
      type(phs_parameters_t), intent(in), optional :: phs_par
      type(string_t), intent(in), optional :: phs_method_in
    end subroutine dispatch_phs
    module subroutine dispatch_sf_channels (sf_channel, sf_string, sf_prop, &
         coll, var_list, sqrts, beam_structure)
      type(sf_channel_t), dimension(:), allocatable, intent(out) :: sf_channel
      type(string_t), intent(out) :: sf_string
      type(sf_prop_t), intent(in) :: sf_prop
      type(phs_channel_collection_t), intent(in) :: coll
      type(var_list_t), intent(in) :: var_list
      real(default), intent(in) :: sqrts
      type(beam_structure_t), intent(in) :: beam_structure
    end subroutine dispatch_sf_channels
  end interface

end module dispatch_phase_space
