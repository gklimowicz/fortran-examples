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

submodule (dispatch_beams) dispatch_beams_s

  implicit none

contains

  module function strfun_mode (name) result (n)
    type(string_t), intent(in) :: name
    integer :: n
    select case (char (name))
    case ("none")
       n = 0
    case ("sf_test_0", "sf_test_1")
       n = 1
    case ("pdf_builtin","pdf_builtin_photon", &
          "lhapdf","lhapdf_photon")
       n = 1
    case ("isr","epa","ewa")
       n = 1
    case ("circe1", "circe2")
       n = 2
    case ("gaussian")
       n = 2
    case ("beam_events")
       n = 2
    case ("energy_scan")
       n = 2
    case default
       n = -1
       call msg_bug ("Structure function '" // char (name) &
            // "' not supported yet")
    end select
  end function strfun_mode

  module subroutine dispatch_sf_config (sf_config, sf_prop, beam_structure, &
         var_list, var_list_global, model, os_data, sqrts, pdg_prc)
    type(sf_config_t), dimension(:), allocatable, intent(out) :: sf_config
    type(sf_prop_t), intent(out) :: sf_prop
    type(beam_structure_t), intent(inout) :: beam_structure
    type(var_list_t), intent(in) :: var_list
    type(var_list_t), intent(inout) :: var_list_global
    class(model_data_t), target, intent(in) :: model
    type(os_data_t), intent(in) :: os_data
    real(default), intent(in) :: sqrts
    class(sf_data_t), allocatable :: sf_data
    type(beam_structure_t) :: beam_structure_tmp
    type(pdg_array_t), dimension(:,:), intent(in) :: pdg_prc
    type(string_t), dimension(:), allocatable :: prt_in
    type(pdg_array_t), dimension(:), allocatable :: pdg_in
    type(flavor_t) :: flv_in
    integer :: n_beam, n_record, i
    beam_structure_tmp = beam_structure
    call beam_structure_tmp%expand (strfun_mode)
    n_record = beam_structure_tmp%get_n_record ()
    allocate (sf_config (n_record))
    n_beam = beam_structure_tmp%get_n_beam ()
    if (n_beam > 0) then
       allocate (prt_in (n_beam), pdg_in (n_beam))
       prt_in = beam_structure_tmp%get_prt ()
       do i = 1, n_beam
          call flv_in%init (prt_in(i), model)
          pdg_in(i) = flv_in%get_pdg ()
       end do
    else
       n_beam = size (pdg_prc, 1)
       allocate (pdg_in (n_beam))
       pdg_in = pdg_prc(:,1)
    end if
    do i = 1, n_record
       call dispatch_sf_data (sf_data, &
            beam_structure_tmp%get_name (i), &
            beam_structure_tmp%get_i_entry (i), &
            sf_prop, var_list, var_list_global, model, os_data, sqrts, &
            pdg_in, pdg_prc, &
            beam_structure_tmp%polarized ())
       call sf_config(i)%init (beam_structure_tmp%get_i_entry (i), sf_data)
       deallocate (sf_data)
    end do
  end subroutine dispatch_sf_config


end submodule dispatch_beams_s

