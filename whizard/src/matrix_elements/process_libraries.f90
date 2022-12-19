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

module process_libraries

  use, intrinsic :: iso_c_binding !NODEP!

  use iso_varying_string, string_t => varying_string
  use physics_defs
  use os_interface
  use model_data
  use particle_specifiers
  use process_constants
  use prclib_interfaces
  use prc_core_def

  implicit none
  private

  public :: strip_equation_lhs
  public :: process_component_def_t
  public :: process_def_t
  public :: process_def_entry_t
  public :: process_def_list_t
  public :: process_library_t

  integer, parameter, public :: STAT_UNKNOWN = 0
  integer, parameter, public :: STAT_OPEN = 1
  integer, parameter, public :: STAT_CONFIGURED = 2
  integer, parameter, public :: STAT_SOURCE = 3
  integer, parameter, public :: STAT_COMPILED = 4
  integer, parameter, public :: STAT_LINKED = 5
  integer, parameter, public :: STAT_ACTIVE = 6

  integer, parameter, public :: ASSOCIATED_BORN = 1
  integer, parameter, public :: ASSOCIATED_REAL = 2
  integer, parameter, public :: ASSOCIATED_VIRT = 3
  integer, parameter, public :: ASSOCIATED_SUB = 4
  integer, parameter, public :: ASSOCIATED_PDF = 5
  integer, parameter, public :: ASSOCIATED_REAL_SING = 6
  integer, parameter, public :: ASSOCIATED_REAL_FIN = 7
  integer, parameter, public :: N_ASSOCIATED_COMPONENTS = 7

  character, dimension(0:6), parameter :: STATUS_LETTER = &
       ["?", "o", "f", "s", "c", "l", "a"]


  type :: process_component_def_t
     private
     type(string_t) :: basename
     logical :: initial = .false.
     integer :: n_in = 0
     integer :: n_out = 0
     integer :: n_tot = 0
     type(prt_spec_t), dimension(:), allocatable :: prt_in
     type(prt_spec_t), dimension(:), allocatable :: prt_out
     type(string_t) :: method
     type(string_t) :: description
     class(prc_core_def_t), allocatable :: core_def
     character(32) :: md5sum = ""
     integer :: nlo_type = BORN
     integer, dimension(N_ASSOCIATED_COMPONENTS) :: associated_components = 0
     logical :: active
     integer :: fixed_emitter = -1
     integer :: alpha_power = 0
     integer :: alphas_power = 0
   contains
     procedure :: write => process_component_def_write
     procedure :: read => process_component_def_read
     procedure :: show => process_component_def_show
     procedure :: compute_md5sum => process_component_def_compute_md5sum
     procedure :: get_def_type_string => process_component_def_get_def_type_string
     procedure :: allocate_driver => process_component_def_allocate_driver
     procedure :: needs_code => process_component_def_needs_code
     procedure :: get_writer_ptr => process_component_def_get_writer_ptr
     procedure :: get_features => process_component_def_get_features
     procedure :: connect => process_component_def_connect
     procedure :: get_core_def_ptr => process_component_get_core_def_ptr
     procedure :: get_n_in  => process_component_def_get_n_in
     procedure :: get_n_out => process_component_def_get_n_out
     procedure :: get_n_tot => process_component_def_get_n_tot
     procedure :: get_prt_in => process_component_def_get_prt_in
     procedure :: get_prt_out => process_component_def_get_prt_out
     procedure :: get_prt_spec_in => process_component_def_get_prt_spec_in
     procedure :: get_prt_spec_out => process_component_def_get_prt_spec_out
     procedure :: get_pdg_in => process_component_def_get_pdg_in
     procedure :: get_md5sum => process_component_def_get_md5sum
     procedure :: get_nlo_type => process_component_def_get_nlo_type
     procedure :: get_associated_born &
          => process_component_def_get_associated_born
     procedure :: get_associated_real_fin &
          => process_component_def_get_associated_real_fin
     procedure :: get_associated_real_sing &
          => process_component_def_get_associated_real_sing
     procedure :: get_associated_subtraction &
          => process_component_def_get_associated_subtraction
     procedure :: get_association_list &
          => process_component_def_get_association_list
     procedure :: can_be_integrated &
          => process_component_def_can_be_integrated
     procedure :: get_associated_real => process_component_def_get_associated_real
     procedure :: get_me_method => process_component_def_get_me_method
     procedure :: get_fixed_emitter => process_component_def_get_fixed_emitter
     procedure :: get_coupling_powers => process_component_def_get_coupling_powers
  end type process_component_def_t

  type :: process_def_t
     private
     type(string_t) :: id
     integer :: num_id = 0
     class(model_data_t), pointer :: model => null ()
     type(string_t) :: model_name
     integer :: n_in  = 0
     integer :: n_initial = 0
     integer :: n_extra = 0
     type(process_component_def_t), dimension(:), allocatable :: initial
     type(process_component_def_t), dimension(:), allocatable :: extra
     character(32) :: md5sum = ""
     logical :: nlo_process = .false.
     logical :: negative_sf = .false.
     logical :: requires_resonances = .false.
   contains
     procedure :: write => process_def_write
     procedure :: read => process_def_read
     procedure :: show => process_def_show
     procedure :: init => process_def_init
     procedure :: set_model_name => process_def_set_model_name
     procedure :: import_component => process_def_import_component
     procedure :: get_n_components => process_def_get_n_components
     procedure :: set_fixed_emitter => process_def_set_fixed_emitter
     procedure :: set_coupling_powers => process_def_set_coupling_powers
     procedure :: set_associated_components => &
          process_def_set_associated_components
     procedure :: compute_md5sum => process_def_compute_md5sum
     procedure :: get_md5sum => process_def_get_md5sum
     procedure :: get_core_def_ptr => process_def_get_core_def_ptr
     procedure :: needs_code => process_def_needs_code
     procedure :: get_pdg_in_1 => process_def_get_pdg_in_1
     procedure :: is_nlo => process_def_is_nlo
     procedure :: get_nlo_type => process_def_get_nlo_type
     procedure :: get_negative_sf => process_def_get_negative_sf
     procedure :: get_n_in => process_def_get_n_in
     procedure :: get_component_def_ptr => process_def_get_component_def_ptr
  end type process_def_t

  type, extends (process_def_t) :: process_def_entry_t
     private
     type(process_def_entry_t), pointer :: next => null ()
  end type process_def_entry_t

  type :: process_def_list_t
     private
     type(process_def_entry_t), pointer :: first => null ()
     type(process_def_entry_t), pointer :: last => null ()
   contains
     procedure :: final => process_def_list_final
     procedure :: write => process_def_list_write
     procedure :: show => process_def_list_show
     procedure :: read => process_def_list_read
     procedure :: append => process_def_list_append
     procedure :: get_n_processes => process_def_list_get_n_processes
     procedure :: get_process_id_list => process_def_list_get_process_id_list
     procedure :: get_process_id_req_resonant => &
          process_def_list_get_process_id_req_resonant
     procedure :: get_process_def_ptr => process_def_list_get_process_def_ptr
     procedure :: contains => process_def_list_contains
     procedure :: get_entry_index => process_def_list_get_entry_index
     procedure :: get_num_id => process_def_list_get_num_id
     procedure :: get_model_name => process_def_list_get_model_name
     procedure :: get_n_in => process_def_list_get_n_in
     procedure :: get_pdg_in_1 => process_def_list_get_pdg_in_1
     procedure :: get_component_list => process_def_list_get_component_list
     procedure :: get_component_description_list => &
          process_def_list_get_component_description_list
     procedure :: req_resonant => process_def_list_req_resonant
  end type process_def_list_t

  type :: process_library_entry_t
     private
     integer :: status = STAT_UNKNOWN
     type(process_def_t), pointer :: def => null ()
     integer :: i_component = 0
     integer :: i_external = 0
     class(prc_core_driver_t), allocatable :: driver
   contains
     procedure :: to_string => process_library_entry_to_string
     procedure :: init => process_library_entry_init
     procedure :: connect => process_library_entry_connect
     procedure :: fill_constants => process_library_entry_fill_constants
  end type process_library_entry_t

  type, extends (process_def_list_t) :: process_library_t
     private
     type(string_t) :: basename
     integer :: n_entries = 0
     logical :: external = .false.
     integer :: status = STAT_UNKNOWN
     logical :: static = .false.
     logical :: driver_exists = .false.
     logical :: makefile_exists = .false.
     integer :: update_counter = 0
     type(process_library_entry_t), dimension(:), allocatable :: entry
     class(prclib_driver_t), allocatable :: driver
     character(32) :: md5sum = ""
   contains
     procedure :: write => process_library_write
     procedure :: show => process_library_show
     procedure :: init => process_library_init
     procedure :: init_static => process_library_init_static
     procedure :: configure => process_library_configure
     procedure :: allocate_entries => process_library_allocate_entries
     procedure :: init_entry => process_library_init_entry
     procedure :: compute_md5sum => process_library_compute_md5sum
     procedure :: write_makefile => process_library_write_makefile
     procedure :: write_driver => process_library_write_driver
     procedure :: update_status => process_library_update_status
     procedure :: make_source => process_library_make_source
     procedure :: make_compile => process_library_make_compile
     procedure :: make_link => process_library_make_link
     procedure :: load => process_library_load
     procedure :: load_entries => process_library_load_entries
     procedure :: unload => process_library_unload
     procedure :: clean => process_library_clean
     procedure :: open => process_library_open
     procedure :: get_name => process_library_get_name
     procedure :: is_active => process_library_is_active
     procedure :: get_status => process_library_get_status
     procedure :: get_update_counter => process_library_get_update_counter
     procedure :: set_status => process_library_set_status
     procedure :: is_loaded => process_library_is_loaded
     procedure :: fill_constants => process_library_fill_constants
     procedure :: connect_process => process_library_connect_process
     procedure :: test_transfer_md5sum => process_library_test_transfer_md5sum
     procedure :: get_nlo_type => process_library_get_nlo_type
     procedure :: get_modellibs_ldflags => process_library_get_modellibs_ldflags
     procedure :: get_static_modelname => process_library_get_static_modelname
  end type process_library_t


  interface
    module subroutine strip_equation_lhs (buffer)
      character(*), intent(inout) :: buffer
    end subroutine strip_equation_lhs
    module subroutine process_component_def_write (object, unit)
      class(process_component_def_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine process_component_def_write
    module subroutine process_component_def_read (component, unit, core_def_templates)
      class(process_component_def_t), intent(out) :: component
      integer, intent(in) :: unit
      type(prc_template_t), dimension(:), intent(in) :: core_def_templates
    end subroutine process_component_def_read
    module subroutine process_component_def_show (object, unit)
      class(process_component_def_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine process_component_def_show
    module subroutine process_component_def_compute_md5sum (component, model)
      class(process_component_def_t), intent(inout) :: component
      class(model_data_t), intent(in), optional, target :: model
    end subroutine process_component_def_compute_md5sum
    module function process_component_def_get_def_type_string (component) result (type_string)
      type(string_t) :: type_string
      class(process_component_def_t), intent(in) :: component
    end function process_component_def_get_def_type_string
    module subroutine process_component_def_allocate_driver (component, driver)
      class(process_component_def_t), intent(in) :: component
      class(prc_core_driver_t), intent(out), allocatable :: driver
    end subroutine process_component_def_allocate_driver
    module function process_component_def_needs_code (component) result (flag)
      class(process_component_def_t), intent(in) :: component
      logical :: flag
    end function process_component_def_needs_code
    module function process_component_def_get_writer_ptr (component) result (writer)
      class(process_component_def_t), intent(in), target :: component
      class(prc_writer_t), pointer :: writer
    end function process_component_def_get_writer_ptr
    module function process_component_def_get_features (component) result (features)
      class(process_component_def_t), intent(in) :: component
      type(string_t), dimension(:), allocatable :: features
    end function process_component_def_get_features
    module subroutine process_component_def_connect &
         (component, lib_driver, i, proc_driver)
      class(process_component_def_t), intent(in) :: component
      class(prclib_driver_t), intent(in) :: lib_driver
      integer, intent(in) :: i
      class(prc_core_driver_t), intent(inout) :: proc_driver
    end subroutine process_component_def_connect
    module function process_component_get_core_def_ptr (component) result (ptr)
      class(process_component_def_t), intent(in), target :: component
      class(prc_core_def_t), pointer :: ptr
    end function process_component_get_core_def_ptr
    module function process_component_def_get_n_in (component) result (n_in)
      class(process_component_def_t), intent(in) :: component
      integer :: n_in
    end function process_component_def_get_n_in
    module function process_component_def_get_n_out (component) result (n_out)
      class(process_component_def_t), intent(in) :: component
      integer :: n_out
    end function process_component_def_get_n_out
    module function process_component_def_get_n_tot (component) result (n_tot)
      class(process_component_def_t), intent(in) :: component
      integer :: n_tot
    end function process_component_def_get_n_tot
    module subroutine process_component_def_get_prt_in (component, prt)
      class(process_component_def_t), intent(in) :: component
      type(string_t), dimension(:), intent(out), allocatable :: prt
    end subroutine process_component_def_get_prt_in
    module subroutine process_component_def_get_prt_out (component, prt)
      class(process_component_def_t), intent(in) :: component
      type(string_t), dimension(:), intent(out), allocatable :: prt
    end subroutine process_component_def_get_prt_out
    module function process_component_def_get_prt_spec_in (component) result (prt)
      class(process_component_def_t), intent(in) :: component
      type(prt_spec_t), dimension(:), allocatable :: prt
    end function process_component_def_get_prt_spec_in
    module function process_component_def_get_prt_spec_out (component) result (prt)
      class(process_component_def_t), intent(in) :: component
      type(prt_spec_t), dimension(:), allocatable :: prt
    end function process_component_def_get_prt_spec_out
    module subroutine process_component_def_get_pdg_in (component, model, pdg)
      class(process_component_def_t), intent(in) :: component
      class(model_data_t), intent(in), target :: model
      integer, intent(out), dimension(:) :: pdg
    end subroutine process_component_def_get_pdg_in
    pure module function process_component_def_get_md5sum (component) result (md5sum)
      class(process_component_def_t), intent(in) :: component
      character(32) :: md5sum
    end function process_component_def_get_md5sum
    elemental module function process_component_def_get_nlo_type &
         (component) result (nlo_type)
      integer :: nlo_type
      class(process_component_def_t), intent(in) :: component
    end function process_component_def_get_nlo_type
    elemental module function process_component_def_get_associated_born &
         (component) result (i_born)
      integer :: i_born
      class(process_component_def_t), intent(in) :: component
    end function process_component_def_get_associated_born
    elemental module function process_component_def_get_associated_real_fin &
         (component) result (i_rfin)
      integer :: i_rfin
      class(process_component_def_t), intent(in) :: component
    end function process_component_def_get_associated_real_fin
    elemental module function process_component_def_get_associated_real_sing &
         (component) result (i_rsing)
      integer :: i_rsing
      class(process_component_def_t), intent(in) :: component
    end function process_component_def_get_associated_real_sing
    elemental module function process_component_def_get_associated_subtraction &
         (component) result (i_sub)
      integer :: i_sub
      class(process_component_def_t), intent(in) :: component
    end function process_component_def_get_associated_subtraction
    elemental module function process_component_def_can_be_integrated &
         (component) result (active)
      logical :: active
      class(process_component_def_t), intent(in) :: component
    end function process_component_def_can_be_integrated
    module function process_component_def_get_association_list &
         (component, i_skip_in) result (list)
      integer, dimension(:), allocatable :: list
      class(process_component_def_t), intent(in) :: component
      integer, intent(in), optional :: i_skip_in
    end function process_component_def_get_association_list
    module function process_component_def_get_associated_real &
         (component) result (i_real)
      integer :: i_real
      class(process_component_def_t), intent(in) :: component
    end function process_component_def_get_associated_real
    elemental module function process_component_def_get_me_method (component) result (method)
      type(string_t) :: method
      class(process_component_def_t), intent(in) :: component
    end function process_component_def_get_me_method
    module function process_component_def_get_fixed_emitter (component) result (emitter)
      integer :: emitter
      class(process_component_def_t), intent(in) :: component
    end function process_component_def_get_fixed_emitter
    pure module subroutine process_component_def_get_coupling_powers &
         (component, alpha_power, alphas_power)
      class(process_component_def_t), intent(in) :: component
      integer, intent(out) :: alpha_power, alphas_power
    end subroutine process_component_def_get_coupling_powers
    module subroutine process_def_write (object, unit)
      class(process_def_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine process_def_write
    module subroutine process_def_read (object, unit, core_def_templates)
      class(process_def_t), intent(out) :: object
      integer, intent(in) :: unit
      type(prc_template_t), dimension(:), intent(in) :: core_def_templates
    end subroutine process_def_read
    module subroutine process_def_show (object, unit)
      class(process_def_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine process_def_show
    module subroutine process_def_init (def, id, &
         model, model_name, n_in, n_components, num_id, &
         nlo_process, negative_sf, requires_resonances)
      class(process_def_t), intent(out) :: def
      type(string_t), intent(in), optional :: id
      class(model_data_t), intent(in), optional, target :: model
      type(string_t), intent(in), optional :: model_name
      integer, intent(in), optional :: n_in
      integer, intent(in), optional :: n_components
      integer, intent(in), optional :: num_id
      logical, intent(in), optional :: nlo_process
      logical, intent(in), optional :: negative_sf
      logical, intent(in), optional :: requires_resonances
    end subroutine process_def_init
    module subroutine process_def_set_model_name (def, model_name)
      class(process_def_t), intent(inout) :: def
      type(string_t), intent(in) :: model_name
    end subroutine process_def_set_model_name
    module subroutine process_def_import_component (def, &
         i, n_out, prt_in, prt_out, method, variant, &
         nlo_type, can_be_integrated)
      class(process_def_t), intent(inout) :: def
      integer, intent(in) :: i
      integer, intent(in), optional :: n_out
      type(prt_spec_t), dimension(:), intent(in), optional :: prt_in
      type(prt_spec_t), dimension(:), intent(in), optional :: prt_out
      type(string_t), intent(in), optional :: method
      integer, intent(in), optional :: nlo_type
      logical, intent(in), optional :: can_be_integrated
      class(prc_core_def_t), &
           intent(inout), allocatable, optional :: variant
    end subroutine process_def_import_component
    module function process_def_get_n_components (def) result (n)
      class(process_def_t), intent(in) :: def
      integer :: n
    end function process_def_get_n_components
    module subroutine process_def_set_fixed_emitter (def, i, emitter)
      class(process_def_t), intent(inout) :: def
      integer, intent(in) :: i, emitter
    end subroutine process_def_set_fixed_emitter
    module subroutine process_def_set_coupling_powers (def, alpha_power, alphas_power)
      class(process_def_t), intent(inout) :: def
      integer, intent(in) :: alpha_power, alphas_power
    end subroutine process_def_set_coupling_powers
    module subroutine process_def_set_associated_components (def, i, &
         i_list, remnant, real_finite, mismatch)
      class(process_def_t), intent(inout) :: def
      logical, intent(in) :: remnant, real_finite, mismatch
      integer, intent(in) :: i
      integer, dimension(:), intent(in) :: i_list
    end subroutine process_def_set_associated_components
    module subroutine process_def_compute_md5sum (def, model)
      class(process_def_t), intent(inout) :: def
      class(model_data_t), intent(in), optional, target :: model
    end subroutine process_def_compute_md5sum
    module function process_def_get_md5sum (def, i_component) result (md5sum)
      class(process_def_t), intent(in) :: def
      integer, intent(in), optional :: i_component
      character(32) :: md5sum
    end function process_def_get_md5sum
    module function process_def_get_core_def_ptr (def, i_component) result (ptr)
      class(process_def_t), intent(in), target :: def
      integer, intent(in) :: i_component
      class(prc_core_def_t), pointer :: ptr
    end function process_def_get_core_def_ptr
    module function process_def_needs_code (def, i_component) result (flag)
      class(process_def_t), intent(in) :: def
      integer, intent(in) :: i_component
      logical :: flag
    end function process_def_needs_code
    module subroutine process_def_get_pdg_in_1 (def, pdg)
      class(process_def_t), intent(in), target :: def
      integer, dimension(:), intent(out) :: pdg
    end subroutine process_def_get_pdg_in_1
    elemental module function process_def_is_nlo (def) result (flag)
      logical :: flag
      class(process_def_t), intent(in) :: def
    end function process_def_is_nlo
    elemental module function process_def_get_nlo_type (def, i_component) result (nlo_type)
      integer :: nlo_type
      class(process_def_t), intent(in) :: def
      integer, intent(in) :: i_component
    end function process_def_get_nlo_type
    elemental module function process_def_get_negative_sf (def) result (neg_sf)
      logical :: neg_sf
      class(process_def_t), intent(in) :: def
    end function process_def_get_negative_sf
    module function process_def_get_n_in (def) result (n_in)
      class(process_def_t), intent(in) :: def
      integer :: n_in
    end function process_def_get_n_in
    module function process_def_get_component_def_ptr (def, i) result (component)
      type(process_component_def_t), pointer :: component
      class(process_def_t), intent(in), target :: def
      integer, intent(in) :: i
    end function process_def_get_component_def_ptr
    module subroutine process_def_list_final (list)
      class(process_def_list_t), intent(inout) :: list
    end subroutine process_def_list_final
    module subroutine process_def_list_write (object, unit, libpath)
      class(process_def_list_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: libpath
    end subroutine process_def_list_write
    module subroutine process_def_list_show (object, unit)
      class(process_def_list_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine process_def_list_show
    module subroutine process_def_list_read (object, unit, core_def_templates)
      class(process_def_list_t), intent(out) :: object
      integer, intent(in) :: unit
      type(prc_template_t), dimension(:), intent(in) :: core_def_templates
    end subroutine process_def_list_read
    module subroutine process_def_list_append (list, entry)
      class(process_def_list_t), intent(inout) :: list
      type(process_def_entry_t), intent(inout), pointer :: entry
    end subroutine process_def_list_append
    module function process_def_list_get_n_processes (list) result (n)
      integer :: n
      class(process_def_list_t), intent(in) :: list
    end function process_def_list_get_n_processes
    module subroutine process_def_list_get_process_id_list (list, id)
      class(process_def_list_t), intent(in) :: list
      type(string_t), dimension(:), allocatable, intent(out) :: id
    end subroutine process_def_list_get_process_id_list
    module subroutine process_def_list_get_process_id_req_resonant (list, id)
      class(process_def_list_t), intent(in) :: list
      type(string_t), dimension(:), allocatable, intent(out) :: id
    end subroutine process_def_list_get_process_id_req_resonant
    module function process_def_list_get_process_def_ptr (list, id) result (entry)
      type(process_def_entry_t), pointer :: entry
      class(process_def_list_t), intent(in) :: list
      type(string_t), intent(in) :: id
    end function process_def_list_get_process_def_ptr
    module function process_def_list_contains (list, id) result (flag)
      logical :: flag
      class(process_def_list_t), intent(in) :: list
      type(string_t), intent(in) :: id
    end function process_def_list_contains
    module function process_def_list_get_entry_index (list, id) result (n)
      integer :: n
      class(process_def_list_t), intent(in) :: list
      type(string_t), intent(in) :: id
    end function process_def_list_get_entry_index
    module function process_def_list_get_num_id (list, id) result (num_id)
      integer :: num_id
      class(process_def_list_t), intent(in) :: list
      type(string_t), intent(in) :: id
    end function process_def_list_get_num_id
    module function process_def_list_get_model_name (list, id) result (model_name)
      type(string_t) :: model_name
      class(process_def_list_t), intent(in) :: list
      type(string_t), intent(in) :: id
    end function process_def_list_get_model_name
    module function process_def_list_get_n_in (list, id) result (n)
      integer :: n
      class(process_def_list_t), intent(in) :: list
      type(string_t), intent(in) :: id
    end function process_def_list_get_n_in
    module subroutine process_def_list_get_pdg_in_1 (list, id, pdg)
      class(process_def_list_t), intent(in) :: list
      type(string_t), intent(in) :: id
      integer, dimension(:), intent(out) :: pdg
    end subroutine process_def_list_get_pdg_in_1
    module subroutine process_def_list_get_component_list (list, id, cid)
      class(process_def_list_t), intent(in) :: list
      type(string_t), intent(in) :: id
      type(string_t), dimension(:), allocatable, intent(out) :: cid
    end subroutine process_def_list_get_component_list
    module subroutine process_def_list_get_component_description_list &
         (list, id, description)
      class(process_def_list_t), intent(in) :: list
      type(string_t), intent(in) :: id
      type(string_t), dimension(:), allocatable, intent(out) :: description
    end subroutine process_def_list_get_component_description_list
    module function process_def_list_req_resonant (list, id) result (flag)
      class(process_def_list_t), intent(in) :: list
      type(string_t), intent(in) :: id
      logical :: flag
    end function process_def_list_req_resonant
    module function process_library_entry_to_string (object) result (string)
      type(string_t) :: string
      class(process_library_entry_t), intent(in) :: object
    end function process_library_entry_to_string
    module subroutine process_library_entry_init (object, &
         status, def, i_component, i_external, driver_template)
      class(process_library_entry_t), intent(out) :: object
      integer, intent(in) :: status
      type(process_def_t), target, intent(in) :: def
      integer, intent(in) :: i_component
      integer, intent(in) :: i_external
      class(prc_core_driver_t), intent(inout), allocatable, optional &
           :: driver_template
    end subroutine process_library_entry_init
    module subroutine process_library_entry_connect (entry, lib_driver, i)
      class(process_library_entry_t), intent(inout) :: entry
      class(prclib_driver_t), intent(in) :: lib_driver
      integer, intent(in) :: i
    end subroutine process_library_entry_connect
    module subroutine process_library_write (object, unit, libpath)
      class(process_library_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: libpath
    end subroutine process_library_write
    module subroutine process_library_show (object, unit)
      class(process_library_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine process_library_show
    module subroutine process_library_init (lib, basename)
      class(process_library_t), intent(out) :: lib
      type(string_t), intent(in) :: basename
    end subroutine process_library_init
    module subroutine process_library_init_static (lib, basename)
      class(process_library_t), intent(out) :: lib
      type(string_t), intent(in) :: basename
    end subroutine process_library_init_static
    module subroutine process_library_configure (lib, os_data)
      class(process_library_t), intent(inout) :: lib
      type(os_data_t), intent(in) :: os_data
    end subroutine process_library_configure
    module subroutine process_library_allocate_entries (lib, n_entries)
      class(process_library_t), intent(inout) :: lib
      integer, intent(in) :: n_entries
    end subroutine process_library_allocate_entries
    module subroutine process_library_init_entry (lib, i, &
         status, def, i_component, i_external, driver_template)
      class(process_library_t), intent(inout) :: lib
      integer, intent(in) :: i
      integer, intent(in) :: status
      type(process_def_t), target, intent(in) :: def
      integer, intent(in) :: i_component
      integer, intent(in) :: i_external
      class(prc_core_driver_t), intent(inout), allocatable, optional &
           :: driver_template
    end subroutine process_library_init_entry
    module subroutine process_library_compute_md5sum (lib, model)
      class(process_library_t), intent(inout) :: lib
      class(model_data_t), intent(in), optional, target :: model
    end subroutine process_library_compute_md5sum
    module subroutine process_library_write_makefile &
         (lib, os_data, force, verbose, testflag, workspace)
      class(process_library_t), intent(inout) :: lib
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: force, verbose
      logical, intent(in), optional :: testflag
      type(string_t), intent(in), optional :: workspace
    end subroutine process_library_write_makefile
    module subroutine process_library_write_driver (lib, force, workspace)
      class(process_library_t), intent(inout) :: lib
      logical, intent(in) :: force
      type(string_t), intent(in), optional :: workspace
    end subroutine process_library_write_driver
    module subroutine process_library_update_status (lib, os_data, workspace)
      class(process_library_t), intent(inout) :: lib
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in), optional :: workspace
    end subroutine process_library_update_status
    module subroutine process_library_make_source &
         (lib, os_data, keep_old_source, workspace)
      class(process_library_t), intent(inout) :: lib
      type(os_data_t), intent(in) :: os_data
      logical, intent(in), optional :: keep_old_source
      type(string_t), intent(in), optional :: workspace
    end subroutine process_library_make_source
    module subroutine process_library_make_compile &
         (lib, os_data, keep_old_source, workspace)
      class(process_library_t), intent(inout) :: lib
      type(os_data_t), intent(in) :: os_data
      logical, intent(in), optional :: keep_old_source
      type(string_t), intent(in), optional :: workspace
    end subroutine process_library_make_compile
    module subroutine process_library_make_link &
         (lib, os_data, keep_old_source, workspace)
      class(process_library_t), intent(inout) :: lib
      type(os_data_t), intent(in) :: os_data
      logical, intent(in), optional :: keep_old_source
      type(string_t), intent(in), optional :: workspace
    end subroutine process_library_make_link
    module subroutine process_library_load (lib, os_data, keep_old_source, workspace)
      class(process_library_t), intent(inout) :: lib
      type(os_data_t), intent(in) :: os_data
      logical, intent(in), optional :: keep_old_source
      type(string_t), intent(in), optional :: workspace
    end subroutine process_library_load
    module subroutine process_library_load_entries (lib)
      class(process_library_t), intent(inout) :: lib
    end subroutine process_library_load_entries
    module subroutine process_library_unload (lib)
      class(process_library_t), intent(inout) :: lib
    end subroutine process_library_unload
    module subroutine process_library_clean (lib, os_data, distclean, workspace)
      class(process_library_t), intent(inout) :: lib
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: distclean
      type(string_t), intent(in), optional :: workspace
    end subroutine process_library_clean
    module subroutine process_library_open (lib)
      class(process_library_t), intent(inout) :: lib
    end subroutine process_library_open
    module function process_library_get_name (lib) result (name)
      class(process_library_t), intent(in) :: lib
      type(string_t) :: name
    end function process_library_get_name
    module function process_library_is_active (lib) result (flag)
      logical :: flag
      class(process_library_t), intent(in) :: lib
    end function process_library_is_active
    module function process_library_get_status (lib, i) result (status)
      class(process_library_t), intent(in) :: lib
      integer, intent(in), optional :: i
      integer :: status
    end function process_library_get_status
    module function process_library_get_update_counter (lib) result (counter)
      class(process_library_t), intent(in) :: lib
      integer :: counter
    end function process_library_get_update_counter
    module subroutine process_library_set_status (lib, status, entries)
      class(process_library_t), intent(inout) :: lib
      integer, intent(in) :: status
      logical, intent(in), optional :: entries
    end subroutine process_library_set_status
    module function process_library_is_loaded (lib) result (flag)
      class(process_library_t), intent(in) :: lib
      logical :: flag
    end function process_library_is_loaded
    module subroutine process_library_entry_fill_constants (entry, driver, data)
      class(process_library_entry_t), intent(in) :: entry
      class(prclib_driver_t), intent(in) :: driver
      type(process_constants_t), intent(out) :: data
    end subroutine process_library_entry_fill_constants
    module subroutine process_library_fill_constants (lib, id, i_component, data)
      class(process_library_t), intent(in) :: lib
      type(string_t), intent(in) :: id
      integer, intent(in) :: i_component
      type(process_constants_t), intent(out) :: data
    end subroutine process_library_fill_constants
    module subroutine process_library_connect_process &
         (lib, id, i_component, data, proc_driver)
      class(process_library_t), intent(in) :: lib
      type(string_t), intent(in) :: id
      integer, intent(in) :: i_component
      type(process_constants_t), intent(out) :: data
      class(prc_core_driver_t), allocatable, intent(out) :: proc_driver
    end subroutine process_library_connect_process
    module subroutine process_library_test_transfer_md5sum (lib, r, e, c)
      class(process_library_t), intent(inout) :: lib
      integer, intent(in) :: r, e, c
    end subroutine process_library_test_transfer_md5sum
    module function process_library_get_nlo_type (lib, id, i_component) result (nlo_type)
      integer :: nlo_type
      class(process_library_t), intent(in) :: lib
      type(string_t), intent(in) :: id
      integer, intent(in) :: i_component
    end function process_library_get_nlo_type
    module function process_library_get_modellibs_ldflags (prc_lib, os_data) result (flags)
      class(process_library_t), intent(in) :: prc_lib
      type(os_data_t), intent(in) :: os_data
      type(string_t) :: flags
    end function process_library_get_modellibs_ldflags
    module function process_library_get_static_modelname (prc_lib, os_data) result (name)
      class(process_library_t), intent(in) :: prc_lib
      type(os_data_t), intent(in) :: os_data
      type(string_t) :: name
    end function process_library_get_static_modelname
  end interface

end module process_libraries

