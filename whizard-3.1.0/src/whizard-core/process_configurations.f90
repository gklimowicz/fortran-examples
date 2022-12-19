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

module process_configurations

  use iso_varying_string, string_t => varying_string
  use models
  use particle_specifiers
  use process_libraries
  use rt_data
  use variables, only: var_list_t

  implicit none
  private

  public :: process_configuration_t

  type :: process_configuration_t
     type(process_def_entry_t), pointer :: entry => null ()
     type(string_t) :: id
     integer :: num_id = 0
   contains
     procedure :: write => process_configuration_write
     procedure :: init => process_configuration_init
     procedure :: setup_component => process_configuration_setup_component
     procedure :: set_fixed_emitter => process_configuration_set_fixed_emitter
     procedure :: set_coupling_powers => process_configuration_set_coupling_powers
     procedure :: set_component_associations => &
          process_configuration_set_component_associations
     procedure :: record => process_configuration_record
  end type process_configuration_t


  interface
    module subroutine process_configuration_write (config, unit)
      class(process_configuration_t), intent(in) :: config
      integer, intent(in), optional :: unit
    end subroutine process_configuration_write
    module subroutine process_configuration_init &
         (config, prc_name, n_in, n_components, model, var_list, &
         nlo_process, negative_sf)
      class(process_configuration_t), intent(out) :: config
      type(string_t), intent(in) :: prc_name
      integer, intent(in) :: n_in
      integer, intent(in) :: n_components
      type(model_t), intent(in), pointer :: model
      type(var_list_t), intent(in) :: var_list
      logical, intent(in), optional :: nlo_process, negative_sf
    end subroutine process_configuration_init
    module subroutine process_configuration_setup_component &
         (config, i_component, prt_in, prt_out, model, var_list, &
          nlo_type, can_be_integrated)
      class(process_configuration_t), intent(inout) :: config
      integer, intent(in) :: i_component
      type(prt_spec_t), dimension(:), intent(in) :: prt_in
      type(prt_spec_t), dimension(:), intent(in) :: prt_out
      type(model_t), pointer, intent(in) :: model
      type(var_list_t), intent(in) :: var_list
      integer, intent(in), optional :: nlo_type
      logical, intent(in), optional :: can_be_integrated
    end subroutine process_configuration_setup_component
    module subroutine process_configuration_set_fixed_emitter &
         (config, i, emitter)
       class(process_configuration_t), intent(inout) :: config
       integer, intent(in) :: i, emitter
    end subroutine process_configuration_set_fixed_emitter
    module subroutine process_configuration_set_coupling_powers &
         (config, alpha_power, alphas_power)
      class(process_configuration_t), intent(inout) :: config
      integer, intent(in) :: alpha_power, alphas_power
    end subroutine process_configuration_set_coupling_powers
    module subroutine process_configuration_set_component_associations &
           (config, i_list, remnant, use_real_finite, mismatch)
      class(process_configuration_t), intent(inout) :: config
      integer, dimension(:), intent(in) :: i_list
      logical, intent(in) :: remnant, use_real_finite, mismatch
    end subroutine process_configuration_set_component_associations
    module subroutine process_configuration_record (config, global)
      class(process_configuration_t), intent(inout) :: config
      type(rt_data_t), intent(inout) :: global
    end subroutine process_configuration_record
  end interface

end module process_configurations
