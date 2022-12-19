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

module sf_base

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use numeric_utils, only: pacify
  use lorentz
  use quantum_numbers
  use pdg_arrays
  use interactions
  use evaluators
  use beams
  use sf_aux
  use sf_mappings

  implicit none
  private

  public :: sf_rescale_t
  public :: sf_rescale_collinear_t
  public :: sf_rescale_real_t
  public :: sf_rescale_dglap_t
  public :: sf_data_t
  public :: sf_config_t
  public :: sf_int_t
  public :: sf_chain_t
  public :: sf_chain_instance_t

  integer, parameter, public :: SF_UNDEFINED = 0
  integer, parameter, public :: SF_INITIAL = 1
  integer, parameter, public :: SF_DONE_LINKS = 2
  integer, parameter, public :: SF_FAILED_MASK = 3
  integer, parameter, public :: SF_DONE_MASK = 4
  integer, parameter, public :: SF_FAILED_CONNECTIONS = 5
  integer, parameter, public :: SF_DONE_CONNECTIONS = 6
  integer, parameter, public :: SF_SEED_KINEMATICS = 10
  integer, parameter, public :: SF_FAILED_KINEMATICS = 11
  integer, parameter, public :: SF_DONE_KINEMATICS = 12
  integer, parameter, public :: SF_FAILED_EVALUATION = 13
  integer, parameter, public :: SF_EVALUATED = 20


  type, abstract :: sf_rescale_t
    integer :: i_beam = 0
  contains
    procedure (sf_rescale_apply), deferred :: apply
    procedure :: set_i_beam => sf_rescale_set_i_beam
  end type sf_rescale_t

  type, extends (sf_rescale_t) :: sf_rescale_collinear_t
     real(default) :: xi_tilde
  contains
    procedure :: apply => sf_rescale_collinear_apply
    procedure :: set => sf_rescale_collinear_set
  end type sf_rescale_collinear_t

  type, extends (sf_rescale_t) :: sf_rescale_real_t
     real(default) :: xi, y
  contains
    procedure :: apply => sf_rescale_real_apply
    procedure :: set => sf_rescale_real_set
  end type sf_rescale_real_t

  type, extends(sf_rescale_t) :: sf_rescale_dglap_t
     real(default), dimension(:), allocatable :: z
   contains
     procedure :: apply => sf_rescale_dglap_apply
     procedure :: set => sf_rescale_dglap_set
  end type sf_rescale_dglap_t

  type, abstract :: sf_data_t
   contains
     procedure (sf_data_write), deferred :: write
     procedure :: is_generator => sf_data_is_generator
     procedure (sf_data_get_int), deferred :: get_n_par
     procedure (sf_data_get_pdg_out), deferred :: get_pdg_out
     procedure (sf_data_allocate_sf_int), deferred :: allocate_sf_int
     procedure :: get_pdf_set => sf_data_get_pdf_set
     procedure :: get_beam_file => sf_data_get_beam_file
  end type sf_data_t

  type :: sf_config_t
     integer, dimension(:), allocatable :: i
     class(sf_data_t), allocatable :: data
   contains
     procedure :: write => sf_config_write
     procedure :: init => sf_config_init
     procedure :: get_pdf_set => sf_config_get_pdf_set
     procedure :: get_beam_file => sf_config_get_beam_file
  end type sf_config_t

  type, abstract, extends (interaction_t) :: sf_int_t
     integer :: status = SF_UNDEFINED
     real(default), dimension(:), allocatable :: mi2
     real(default), dimension(:), allocatable :: mr2
     real(default), dimension(:), allocatable :: mo2
     integer :: on_shell_mode = KEEP_ENERGY
     logical :: qmin_defined = .false.
     logical :: qmax_defined = .false.
     real(default), dimension(:), allocatable :: qmin
     real(default), dimension(:), allocatable :: qmax
     integer, dimension(:), allocatable :: beam_index
     integer, dimension(:), allocatable :: incoming
     integer, dimension(:), allocatable :: radiated
     integer, dimension(:), allocatable :: outgoing
     integer, dimension(:), allocatable :: par_index
     integer, dimension(:), allocatable :: par_primary
   contains
     procedure :: base_write => sf_int_base_write
     procedure (sf_int_type_string), deferred :: type_string
     procedure (sf_int_write), deferred :: write
     procedure :: base_init => sf_int_base_init
     procedure :: set_incoming => sf_int_set_incoming
     procedure :: set_radiated => sf_int_set_radiated
     procedure :: set_outgoing => sf_int_set_outgoing
     procedure (sf_int_init), deferred :: init
     procedure :: setup_constants => sf_int_setup_constants
     procedure :: set_beam_index => sf_int_set_beam_index
     procedure :: set_par_index => sf_int_set_par_index
     generic :: seed_kinematics => sf_int_receive_momenta
     generic :: seed_kinematics => sf_int_seed_momenta
     generic :: seed_kinematics => sf_int_seed_energies
     procedure :: sf_int_receive_momenta
     procedure :: sf_int_seed_momenta
     procedure :: sf_int_seed_energies
     procedure :: is_generator => sf_int_is_generator
     procedure :: generate_free => sf_int_generate_free
     procedure (sf_int_complete_kinematics), deferred :: complete_kinematics
     procedure (sf_int_inverse_kinematics), deferred :: inverse_kinematics
     procedure :: split_momentum => sf_int_split_momentum
     procedure :: split_momenta => sf_int_split_momenta
     procedure :: reduce_momenta => sf_int_reduce_momenta
     procedure :: recover_x => sf_int_recover_x
     procedure :: base_recover_x => sf_int_recover_x
     procedure (sf_int_apply), deferred :: apply
     procedure :: get_n_in => sf_int_get_n_in
     procedure :: get_n_rad => sf_int_get_n_rad
     procedure :: get_n_out => sf_int_get_n_out
     procedure :: get_n_states => sf_int_get_n_states
     procedure :: get_state => sf_int_get_state
     procedure :: get_values => sf_int_get_values
     procedure :: compute_values => sf_int_compute_values
     procedure :: compute_value => sf_int_compute_value
  end type sf_int_t

  type :: sf_instance_t
     class(sf_int_t), allocatable :: int
     type(evaluator_t) :: eval
     real(default), dimension(:,:), allocatable :: r
     real(default), dimension(:,:), allocatable :: rb
     real(default), dimension(:), allocatable :: f
     logical, dimension(:), allocatable :: m
     real(default), dimension(:), allocatable :: x
     real(default), dimension(:), allocatable :: xb
  end type sf_instance_t

  type, extends (beam_t) :: sf_chain_t
     type(beam_data_t), pointer :: beam_data => null ()
     integer :: n_in = 0
     integer :: n_strfun = 0
     integer :: n_par = 0
     integer :: n_bound = 0
     type(sf_instance_t), dimension(:), allocatable :: sf
     logical :: trace_enable = .false.
     integer :: trace_unit = 0
   contains
     procedure :: final => sf_chain_final
     procedure :: write => sf_chain_write
     procedure :: init => sf_chain_init
     procedure :: receive_beam_momenta => sf_chain_receive_beam_momenta
     procedure :: set_beam_momenta => sf_chain_set_beam_momenta
     procedure :: set_strfun => sf_chain_set_strfun
     procedure :: get_n_par => sf_chain_get_n_par
     procedure :: get_n_bound => sf_chain_get_n_bound
     procedure :: get_beam_int_ptr => sf_chain_get_beam_int_ptr
     procedure :: setup_tracing => sf_chain_setup_tracing
     procedure :: final_tracing => sf_chain_final_tracing
     procedure :: write_trace_header => sf_chain_write_trace_header
     procedure :: trace => sf_chain_trace
  end type sf_chain_t

  type, extends (beam_t) :: sf_chain_instance_t
     type(sf_chain_t), pointer :: config => null ()
     integer :: status = SF_UNDEFINED
     type(sf_instance_t), dimension(:), allocatable :: sf
     integer, dimension(:), allocatable :: out_sf
     integer, dimension(:), allocatable :: out_sf_i
     integer :: out_eval = 0
     integer, dimension(:), allocatable :: out_eval_i
     integer :: selected_channel = 0
     real(default), dimension(:,:), allocatable :: p, pb
     real(default), dimension(:,:), allocatable :: r, rb
     real(default), dimension(:), allocatable :: f
     real(default), dimension(:), allocatable :: x, xb
     logical, dimension(:), allocatable :: bound
     real(default) :: x_free = 1
     type(sf_channel_t), dimension(:), allocatable :: channel
   contains
     procedure :: final => sf_chain_instance_final
     procedure :: write => sf_chain_instance_write
     procedure :: init => sf_chain_instance_init
     procedure :: select_channel => sf_chain_instance_select_channel
     procedure :: set_channel => sf_chain_instance_set_channel
     procedure :: link_interactions => sf_chain_instance_link_interactions
     procedure :: exchange_mask => sf_chain_exchange_mask
     procedure :: init_evaluators => sf_chain_instance_init_evaluators
     procedure :: write_interaction => sf_chain_instance_write_interaction
     procedure :: compute_kinematics => sf_chain_instance_compute_kinematics
     procedure :: inverse_kinematics => sf_chain_instance_inverse_kinematics
     procedure :: recover_kinematics => sf_chain_instance_recover_kinematics
     procedure :: return_beam_momenta => sf_chain_instance_return_beam_momenta
     procedure :: evaluate => sf_chain_instance_evaluate
     procedure :: get_out_momenta => sf_chain_instance_get_out_momenta
     procedure :: get_out_int_ptr => sf_chain_instance_get_out_int_ptr
     procedure :: get_out_i => sf_chain_instance_get_out_i
     procedure :: get_out_mask => sf_chain_instance_get_out_mask
     procedure :: get_mcpar => sf_chain_instance_get_mcpar
     procedure :: get_f => sf_chain_instance_get_f
     procedure :: get_status => sf_chain_instance_get_status
     procedure :: get_matrix_elements => sf_chain_instance_get_matrix_elements
     procedure :: get_beam_int_ptr => sf_chain_instance_get_beam_int_ptr
     procedure :: get_n_sub => sf_chain_instance_get_n_sub
  end type sf_chain_instance_t


  abstract interface
     subroutine sf_rescale_apply (func, x)
       import
       class(sf_rescale_t), intent(in) :: func
       real(default), intent(inout) :: x
     end subroutine sf_rescale_apply
  end interface

  abstract interface
     subroutine sf_data_write (data, unit, verbose)
       import
       class(sf_data_t), intent(in) :: data
       integer, intent(in), optional :: unit
       logical, intent(in), optional :: verbose
     end subroutine sf_data_write
  end interface

  abstract interface
     function sf_data_get_int (data) result (n)
       import
       class(sf_data_t), intent(in) :: data
       integer :: n
     end function sf_data_get_int
  end interface

  abstract interface
     subroutine sf_data_get_pdg_out (data, pdg_out)
       import
       class(sf_data_t), intent(in) :: data
       type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
     end subroutine sf_data_get_pdg_out
  end interface

  abstract interface
     subroutine sf_data_allocate_sf_int (data, sf_int)
       import
       class(sf_data_t), intent(in) :: data
       class(sf_int_t), intent(inout), allocatable :: sf_int
     end subroutine sf_data_allocate_sf_int
  end interface

  abstract interface
     function sf_int_type_string (object) result (string)
       import
       class(sf_int_t), intent(in) :: object
       type(string_t) :: string
     end function sf_int_type_string
  end interface

  abstract interface
     subroutine sf_int_write (object, unit, testflag)
       import
       class(sf_int_t), intent(in) :: object
       integer, intent(in), optional :: unit
       logical, intent(in), optional :: testflag
     end subroutine sf_int_write
  end interface

  abstract interface
     subroutine sf_int_init (sf_int, data)
       import
       class(sf_int_t), intent(out) :: sf_int
       class(sf_data_t), intent(in), target :: data
     end subroutine sf_int_init
  end interface

  abstract interface
     subroutine sf_int_complete_kinematics (sf_int, x, xb, f, r, rb, map)
       import
       class(sf_int_t), intent(inout) :: sf_int
       real(default), dimension(:), intent(out) :: x
       real(default), dimension(:), intent(out) :: xb
       real(default), intent(out) :: f
       real(default), dimension(:), intent(in) :: r
       real(default), dimension(:), intent(in) :: rb
       logical, intent(in) :: map
     end subroutine sf_int_complete_kinematics
  end interface

  abstract interface
     subroutine sf_int_inverse_kinematics (sf_int, x, xb, f, r, rb, map, set_momenta)
       import
       class(sf_int_t), intent(inout) :: sf_int
       real(default), dimension(:), intent(in) :: x
       real(default), dimension(:), intent(in) :: xb
       real(default), intent(out) :: f
       real(default), dimension(:), intent(out) :: r
       real(default), dimension(:), intent(out) :: rb
       logical, intent(in) :: map
       logical, intent(in), optional :: set_momenta
     end subroutine sf_int_inverse_kinematics
  end interface

  abstract interface
     subroutine sf_int_apply (sf_int, scale, negative_sf, rescale, i_sub)
       import
       class(sf_int_t), intent(inout) :: sf_int
       real(default), intent(in) :: scale
       logical, intent(in), optional :: negative_sf
       class(sf_rescale_t), intent(in), optional :: rescale
       integer, intent(in), optional :: i_sub
     end subroutine sf_int_apply
  end interface


  interface
    module subroutine sf_rescale_set_i_beam (func, i_beam)
      class(sf_rescale_t), intent(inout) :: func
      integer, intent(in) :: i_beam
    end subroutine sf_rescale_set_i_beam
    module subroutine sf_rescale_collinear_apply (func, x)
      class(sf_rescale_collinear_t), intent(in) :: func
      real(default), intent(inout) :: x
    end subroutine sf_rescale_collinear_apply
    module subroutine sf_rescale_collinear_set (func, xi_tilde)
      class(sf_rescale_collinear_t), intent(inout) :: func
      real(default), intent(in) :: xi_tilde
    end subroutine sf_rescale_collinear_set
    module subroutine sf_rescale_real_apply (func, x)
      class(sf_rescale_real_t), intent(in) :: func
      real(default), intent(inout) :: x
    end subroutine sf_rescale_real_apply
    module subroutine sf_rescale_real_set (func, xi, y)
      class(sf_rescale_real_t), intent(inout) :: func
      real(default), intent(in) :: xi, y
    end subroutine sf_rescale_real_set
    module subroutine sf_rescale_dglap_apply (func, x)
      class(sf_rescale_dglap_t), intent(in) :: func
      real(default), intent(inout) :: x
    end subroutine sf_rescale_dglap_apply
    module subroutine sf_rescale_dglap_set (func, z)
      class(sf_rescale_dglap_t), intent(inout) :: func
      real(default), dimension(:), intent(in) :: z
    end subroutine sf_rescale_dglap_set
    module function sf_data_is_generator (data) result (flag)
      class(sf_data_t), intent(in) :: data
      logical :: flag
    end function sf_data_is_generator
    elemental module function sf_data_get_pdf_set (data) result (pdf_set)
      class(sf_data_t), intent(in) :: data
      integer :: pdf_set
    end function sf_data_get_pdf_set
    module function sf_data_get_beam_file (data) result (file)
      class(sf_data_t), intent(in) :: data
      type(string_t) :: file
    end function sf_data_get_beam_file
    module subroutine sf_config_write (object, unit, verbose)
      class(sf_config_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine sf_config_write
    module subroutine sf_config_init (sf_config, i_beam, sf_data)
      class(sf_config_t), intent(out) :: sf_config
      integer, dimension(:), intent(in) :: i_beam
      class(sf_data_t), intent(in) :: sf_data
    end subroutine sf_config_init
    elemental module function sf_config_get_pdf_set (sf_config) result (pdf_set)
      class(sf_config_t), intent(in) :: sf_config
      integer :: pdf_set
    end function sf_config_get_pdf_set
    module function sf_config_get_beam_file (sf_config) result (file)
      class(sf_config_t), intent(in) :: sf_config
      type(string_t) :: file
    end function sf_config_get_beam_file
    module subroutine sf_int_base_write (object, unit, testflag)
      class(sf_int_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine sf_int_base_write
    module subroutine sf_int_base_init &
         (sf_int, mask, mi2, mr2, mo2, qmin, qmax, hel_lock)
      class(sf_int_t), intent(out) :: sf_int
      type (quantum_numbers_mask_t), dimension(:), intent(in) :: mask
      real(default), dimension(:), intent(in) :: mi2, mr2, mo2
      real(default), dimension(:), intent(in), optional :: qmin, qmax
      integer, dimension(:), intent(in), optional :: hel_lock
    end subroutine sf_int_base_init
    module subroutine sf_int_set_incoming (sf_int, incoming)
      class(sf_int_t), intent(inout) :: sf_int
      integer, dimension(:), intent(in) :: incoming
    end subroutine sf_int_set_incoming
    module subroutine sf_int_set_radiated (sf_int, radiated)
      class(sf_int_t), intent(inout) :: sf_int
      integer, dimension(:), intent(in) :: radiated
    end subroutine sf_int_set_radiated
    module subroutine sf_int_set_outgoing (sf_int, outgoing)
      class(sf_int_t), intent(inout) :: sf_int
      integer, dimension(:), intent(in) :: outgoing
    end subroutine sf_int_set_outgoing
    module subroutine sf_int_setup_constants (sf_int)
      class(sf_int_t), intent(inout), target :: sf_int
    end subroutine sf_int_setup_constants
    module subroutine sf_int_set_beam_index (sf_int, beam_index)
      class(sf_int_t), intent(inout) :: sf_int
      integer, dimension(:), intent(in) :: beam_index
    end subroutine sf_int_set_beam_index
    module subroutine sf_int_set_par_index (sf_int, par_index)
      class(sf_int_t), intent(inout) :: sf_int
      integer, dimension(:), intent(in) :: par_index
    end subroutine sf_int_set_par_index
    module subroutine sf_int_receive_momenta (sf_int)
      class(sf_int_t), intent(inout) :: sf_int
    end subroutine sf_int_receive_momenta
    module subroutine sf_int_seed_momenta (sf_int, k)
      class(sf_int_t), intent(inout) :: sf_int
      type(vector4_t), dimension(:), intent(in) :: k
    end subroutine sf_int_seed_momenta
    module subroutine sf_int_seed_energies (sf_int, E)
      class(sf_int_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: E
      type(vector4_t), dimension(:), allocatable :: k
    end subroutine sf_int_seed_energies
    module function sf_int_is_generator (sf_int) result (flag)
      class(sf_int_t), intent(in) :: sf_int
      logical :: flag
    end function sf_int_is_generator
    module subroutine sf_int_generate_free (sf_int, r, rb,  x_free)
      class(sf_int_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(inout) :: x_free
    end subroutine sf_int_generate_free
    module subroutine sf_int_split_momentum (sf_int, x, xb)
      class(sf_int_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
    end subroutine sf_int_split_momentum
    module subroutine sf_int_split_momenta (sf_int, x, xb)
      class(sf_int_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
    end subroutine sf_int_split_momenta
    module subroutine sf_int_reduce_momenta (sf_int, x)
      class(sf_int_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
    end subroutine sf_int_reduce_momenta
    module subroutine sf_int_recover_x (sf_int, x, xb, x_free)
      class(sf_int_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_int_recover_x
    pure module function sf_int_get_n_in (object) result (n_in)
      class(sf_int_t), intent(in) :: object
      integer :: n_in
    end function sf_int_get_n_in
    pure module function sf_int_get_n_rad (object) result (n_rad)
      class(sf_int_t), intent(in) :: object
      integer :: n_rad
    end function sf_int_get_n_rad
    pure module function sf_int_get_n_out (object) result (n_out)
      class(sf_int_t), intent(in) :: object
      integer :: n_out
    end function sf_int_get_n_out
    module function sf_int_get_n_states (sf_int) result (n_states)
      class(sf_int_t), intent(in) :: sf_int
      integer :: n_states
    end function sf_int_get_n_states
    module function sf_int_get_state (sf_int, i) result (qn)
      class(sf_int_t), intent(in) :: sf_int
      type(quantum_numbers_t), dimension(:), allocatable :: qn
      integer, intent(in) :: i
    end function sf_int_get_state
    module subroutine sf_int_get_values (sf_int, value)
      class(sf_int_t), intent(in) :: sf_int
      real(default), dimension(:), intent(out) :: value
    end subroutine sf_int_get_values
    module subroutine sf_int_compute_values (sf_int, value, x, xb, scale, E)
      class(sf_int_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: value
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(in) :: scale
      real(default), dimension(:), intent(in), optional :: E
    end subroutine sf_int_compute_values
    module subroutine sf_int_compute_value &
         (sf_int, i_state, value, x, xb, scale, E)
      class(sf_int_t), intent(inout) :: sf_int
      integer, intent(in) :: i_state
      real(default), intent(out) :: value
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(in) :: scale
      real(default), dimension(:), intent(in), optional :: E
    end subroutine sf_int_compute_value
    module subroutine sf_chain_final (object)
      class(sf_chain_t), intent(inout) :: object
    end subroutine sf_chain_final
    module subroutine sf_chain_write (object, unit)
      class(sf_chain_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_chain_write
    module subroutine sf_chain_init (sf_chain, beam_data, sf_config)
      class(sf_chain_t), intent(out) :: sf_chain
      type(beam_data_t), intent(in), target :: beam_data
      type(sf_config_t), dimension(:), intent(in), optional, target :: sf_config
    end subroutine sf_chain_init
    module subroutine sf_chain_receive_beam_momenta (sf_chain)
      class(sf_chain_t), intent(inout), target :: sf_chain
      type(interaction_t), pointer :: beam_int
    end subroutine sf_chain_receive_beam_momenta
    module subroutine sf_chain_set_beam_momenta (sf_chain, p)
      class(sf_chain_t), intent(inout) :: sf_chain
      type(vector4_t), dimension(:), intent(in) :: p
    end subroutine sf_chain_set_beam_momenta
    module subroutine sf_chain_set_strfun (sf_chain, i, beam_index, data)
      class(sf_chain_t), intent(inout) :: sf_chain
      integer, intent(in) :: i
      integer, dimension(:), intent(in) :: beam_index
      class(sf_data_t), intent(in), target :: data
    end subroutine sf_chain_set_strfun
    module function sf_chain_get_n_par (sf_chain) result (n)
      class(sf_chain_t), intent(in) :: sf_chain
      integer :: n
    end function sf_chain_get_n_par
    module function sf_chain_get_n_bound (sf_chain) result (n)
      class(sf_chain_t), intent(in) :: sf_chain
      integer :: n
    end function sf_chain_get_n_bound
    module function sf_chain_get_beam_int_ptr (sf_chain) result (int)
      type(interaction_t), pointer :: int
      class(sf_chain_t), intent(in), target :: sf_chain
    end function sf_chain_get_beam_int_ptr
    module subroutine sf_chain_setup_tracing (sf_chain, file)
      class(sf_chain_t), intent(inout) :: sf_chain
      type(string_t), intent(in) :: file
    end subroutine sf_chain_setup_tracing
    module subroutine sf_chain_final_tracing (sf_chain)
      class(sf_chain_t), intent(inout) :: sf_chain
    end subroutine sf_chain_final_tracing
    module subroutine sf_chain_write_trace_header (sf_chain)
      class(sf_chain_t), intent(in) :: sf_chain
    end subroutine sf_chain_write_trace_header
    module subroutine sf_chain_trace (sf_chain, c_sel, p, x, f, sf_sum)
      class(sf_chain_t), intent(in) :: sf_chain
      integer, intent(in) :: c_sel
      real(default), dimension(:,:), intent(in) :: p
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: f
      real(default), intent(in) :: sf_sum
    end subroutine sf_chain_trace
    module subroutine sf_chain_instance_final (object)
      class(sf_chain_instance_t), intent(inout) :: object
    end subroutine sf_chain_instance_final
    module subroutine sf_chain_instance_write (object, unit, col_verbose)
      class(sf_chain_instance_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: col_verbose
    end subroutine sf_chain_instance_write
    module subroutine sf_chain_instance_init (chain, config, n_channel)
      class(sf_chain_instance_t), intent(out), target :: chain
      type(sf_chain_t), intent(in), target :: config
      integer, intent(in) :: n_channel
    end subroutine sf_chain_instance_init
    module subroutine sf_chain_instance_select_channel (chain, channel)
      class(sf_chain_instance_t), intent(inout) :: chain
      integer, intent(in), optional :: channel
    end subroutine sf_chain_instance_select_channel
    module subroutine sf_chain_instance_set_channel (chain, c, channel)
      class(sf_chain_instance_t), intent(inout) :: chain
      integer, intent(in) :: c
      type(sf_channel_t), intent(in) :: channel
    end subroutine sf_chain_instance_set_channel
    module subroutine sf_chain_instance_link_interactions (chain)
      class(sf_chain_instance_t), intent(inout), target :: chain
    end subroutine sf_chain_instance_link_interactions
    module subroutine sf_chain_exchange_mask (chain)
      class(sf_chain_instance_t), intent(inout), target :: chain
    end subroutine sf_chain_exchange_mask
    module subroutine sf_chain_instance_init_evaluators (chain, extended_sf)
      class(sf_chain_instance_t), intent(inout), target :: chain
      logical, intent(in), optional :: extended_sf
    end subroutine sf_chain_instance_init_evaluators
    module subroutine sf_chain_instance_write_interaction &
         (chain, i_sf, i_int, unit)
      class(sf_chain_instance_t), intent(in) :: chain
      integer, intent(in) :: i_sf, i_int
      integer, intent(in) :: unit
    end subroutine sf_chain_instance_write_interaction
    module subroutine sf_chain_instance_compute_kinematics (chain, c_sel, p_in)
      class(sf_chain_instance_t), intent(inout), target :: chain
      integer, intent(in) :: c_sel
      real(default), dimension(:), intent(in) :: p_in
    end subroutine sf_chain_instance_compute_kinematics
    module subroutine sf_chain_instance_inverse_kinematics (chain, x, xb)
      class(sf_chain_instance_t), intent(inout), target :: chain
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
    end subroutine sf_chain_instance_inverse_kinematics
    module subroutine sf_chain_instance_recover_kinematics (chain, c_sel)
      class(sf_chain_instance_t), intent(inout), target :: chain
      integer, intent(in) :: c_sel
    end subroutine sf_chain_instance_recover_kinematics
    module subroutine sf_chain_instance_return_beam_momenta (chain)
      class(sf_chain_instance_t), intent(in), target :: chain
      type(interaction_t), pointer :: int
    end subroutine sf_chain_instance_return_beam_momenta
    module subroutine sf_chain_instance_evaluate &
         (chain, scale, negative_sf, sf_rescale)
      class(sf_chain_instance_t), intent(inout), target :: chain
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(inout), optional :: sf_rescale
    end subroutine sf_chain_instance_evaluate
    module subroutine sf_chain_instance_get_out_momenta (chain, p)
      class(sf_chain_instance_t), intent(in), target :: chain
      type(vector4_t), dimension(:), intent(out) :: p
    end subroutine sf_chain_instance_get_out_momenta
    module function sf_chain_instance_get_out_int_ptr (chain) result (int)
      class(sf_chain_instance_t), intent(in), target :: chain
      type(interaction_t), pointer :: int
    end function sf_chain_instance_get_out_int_ptr
    module function sf_chain_instance_get_out_i (chain, j) result (i)
      class(sf_chain_instance_t), intent(in) :: chain
      integer, intent(in) :: j
      integer :: i
    end function sf_chain_instance_get_out_i
    module function sf_chain_instance_get_out_mask (chain) result (mask)
      class(sf_chain_instance_t), intent(in), target :: chain
      type(quantum_numbers_mask_t), dimension(:), allocatable :: mask
    end function sf_chain_instance_get_out_mask
    module subroutine sf_chain_instance_get_mcpar (chain, c, r)
      class(sf_chain_instance_t), intent(in) :: chain
      integer, intent(in) :: c
      real(default), dimension(:), intent(out) :: r
    end subroutine sf_chain_instance_get_mcpar
    module function sf_chain_instance_get_f (chain, c) result (f)
      class(sf_chain_instance_t), intent(in) :: chain
      integer, intent(in) :: c
      real(default) :: f
    end function sf_chain_instance_get_f
    module function sf_chain_instance_get_status (chain) result (status)
      class(sf_chain_instance_t), intent(in) :: chain
      integer :: status
    end function sf_chain_instance_get_status
    module subroutine sf_chain_instance_get_matrix_elements (chain, i, ff)
       class(sf_chain_instance_t), intent(in) :: chain
       integer, intent(in) :: i
       real(default), intent(out), dimension(:), allocatable :: ff
    end subroutine sf_chain_instance_get_matrix_elements
    module function sf_chain_instance_get_beam_int_ptr (chain) result (int)
      type(interaction_t), pointer :: int
      class(sf_chain_instance_t), intent(in), target :: chain
    end function sf_chain_instance_get_beam_int_ptr
    module function sf_chain_instance_get_n_sub (chain) result (n_sub)
      type(interaction_t), pointer :: int
      class(sf_chain_instance_t), intent(in), target :: chain
      integer :: n_sub
    end function sf_chain_instance_get_n_sub
  end interface

end module sf_base
