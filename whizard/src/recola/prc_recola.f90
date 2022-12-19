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

module prc_recola

  use kinds
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use diagnostics
  use lorentz
  use physics_defs
  use variables, only: var_list_t
  use os_interface, only: os_data_t
  use sm_qcd, only: qcd_t
  use model_data, only: model_data_t

  use prc_core, only: prc_core_state_t
  use prc_core_def, only: prc_core_driver_t, prc_core_def_t
  use prc_external
  use process_libraries, only: process_library_t

  implicit none
  private

  public :: abort_if_recola_not_active
  public :: recola_def_t
  public :: prc_recola_t

  integer, parameter :: RECOLA_UNDEFINED = 0, RECOLA_QCD = 1, &
       RECOLA_EW = 2, RECOLA_FULL = 3


  type, extends (prc_external_def_t) :: recola_def_t
     type(string_t) :: suffix
     type(string_t) :: order
     integer :: alpha_power = 0
     integer :: alphas_power = 0
     integer :: corr = RECOLA_UNDEFINED
  contains
    procedure, nopass :: type_string => recola_def_type_string
    procedure :: write => recola_def_write
    procedure :: read => recola_def_read
    procedure :: init => recola_def_init
    procedure :: allocate_driver => recola_def_allocate_driver
  end type recola_def_t

  type, extends (prc_external_writer_t) :: recola_writer_t
     private
     type(string_t) :: id
     type(string_t) :: order
     integer :: alpha_power = 0
     integer :: alphas_power = 0
  contains
    procedure, nopass :: type_name => recola_writer_type_name
    procedure :: set_id => recola_writer_set_id
    procedure :: set_order => recola_writer_set_order
    procedure :: set_coupling_powers => recola_writer_set_coupling_powers
    procedure :: write_makefile_code => recola_writer_write_makefile_code
    procedure :: register_processes => prc_recola_register_processes
  end type recola_writer_t

  type, extends (prc_external_driver_t) :: recola_driver_t
  contains
    procedure, nopass :: type_name => recola_driver_type_name
  end type recola_driver_t

  type, extends (prc_external_t) :: prc_recola_t
     integer, dimension(:), allocatable :: recola_ids
     integer, dimension(:,:), allocatable :: color_state
     integer :: n_f = 0
     logical :: helicity_and_color_arrays_are_replaced = .false.
  contains
    procedure :: write_name => prc_recola_write_name
    procedure :: has_matrix_element => prc_recola_has_matrix_element
    procedure :: write => prc_recola_write
    procedure :: allocate_workspace => prc_recola_allocate_workspace
    procedure :: get_alpha_power => prc_recola_get_alpha_power
    procedure :: get_alphas_power => prc_recola_get_alphas_power
    procedure :: compute_alpha_s => prc_recola_compute_alpha_s
    procedure :: includes_polarization => prc_recola_includes_polarization
    procedure :: prepare_external_code => &
         prc_recola_prepare_external_code
    procedure :: set_parameters => prc_recola_set_parameters
    procedure :: init => prc_recola_init
    procedure :: replace_helicity_and_color_arrays => &
         prc_recola_replace_helicity_and_color_arrays
    procedure :: compute_amplitude => prc_recola_compute_amplitude
    procedure :: compute_sqme => prc_recola_compute_sqme
    procedure :: compute_sqme_virt => prc_recola_compute_sqme_virt
    procedure :: get_alpha_qed => prc_recola_get_alpha_qed
    procedure :: compute_sqme_color_c_raw => prc_recola_compute_sqme_color_c_raw
  end type prc_recola_t

  type, extends (prc_external_state_t) :: recola_state_t
  contains
    procedure :: write => recola_state_write
  end type recola_state_t




  interface
    module subroutine abort_if_recola_not_active ()
    end subroutine abort_if_recola_not_active
    module function recola_def_type_string () result (string)
      type(string_t) :: string
    end function recola_def_type_string
    module subroutine recola_def_write (object, unit)
      class(recola_def_t), intent(in) :: object
      integer, intent(in) :: unit
    end subroutine recola_def_write
    module subroutine recola_def_read (object, unit)
      class(recola_def_t), intent(out) :: object
      integer, intent(in) :: unit
    end subroutine recola_def_read
    module function recola_writer_type_name () result (string)
      type(string_t) :: string
    end function recola_writer_type_name
    module subroutine recola_writer_set_id (writer, id)
      class(recola_writer_t), intent(inout) :: writer
      type(string_t), intent(in) :: id
    end subroutine recola_writer_set_id
    module subroutine recola_writer_set_order (writer, order)
      class(recola_writer_t), intent(inout) :: writer
      type(string_t), intent(in) :: order
    end subroutine recola_writer_set_order
    module subroutine recola_writer_set_coupling_powers &
         (writer, alpha_power, alphas_power)
      class(recola_writer_t), intent(inout) :: writer
      integer, intent(in) :: alpha_power
      integer, intent(in) :: alphas_power
    end subroutine recola_writer_set_coupling_powers
    module subroutine recola_writer_write_makefile_code &
         (writer, unit, id, os_data, verbose, testflag)
      class(recola_writer_t), intent(in) :: writer
      integer, intent(in) :: unit
      type(string_t), intent(in) :: id
      type(os_data_t), intent(in) :: os_data
      logical, intent(in) :: verbose
      logical, intent(in), optional :: testflag
    end subroutine recola_writer_write_makefile_code
    module subroutine prc_recola_register_processes (writer, recola_ids)
      class(recola_writer_t), intent(in) :: writer
      integer, dimension (:), intent(inout) :: recola_ids
    end subroutine prc_recola_register_processes
    module function recola_driver_type_name () result (type)
      type(string_t) :: type
    end function recola_driver_type_name
    module subroutine prc_recola_write_name (object, unit)
      class(prc_recola_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_recola_write_name
    module function prc_recola_has_matrix_element (object) result (flag)
      logical :: flag
      class(prc_recola_t), intent(in) :: object
    end function prc_recola_has_matrix_element
    module subroutine prc_recola_write (object, unit)
      class(prc_recola_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine prc_recola_write
    module subroutine recola_state_write (object, unit)
      class(recola_state_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine recola_state_write
    module function prc_recola_get_alpha_power (object) result (p)
      class(prc_recola_t), intent(in) :: object
      integer :: p
    end function prc_recola_get_alpha_power
    module function prc_recola_get_alphas_power (object) result (p)
      class(prc_recola_t), intent(in) :: object
      integer :: p
    end function prc_recola_get_alphas_power
    module subroutine prc_recola_compute_alpha_s (object, core_state, ren_scale)
      class(prc_recola_t), intent(in) :: object
      class(prc_external_state_t), intent(inout) :: core_state
      real(default), intent(in) :: ren_scale
    end subroutine prc_recola_compute_alpha_s
    module function prc_recola_includes_polarization (object) result (polarized)
      logical :: polarized
      class(prc_recola_t), intent(in) :: object
    end function prc_recola_includes_polarization
    module subroutine prc_recola_prepare_external_code &
         (core, flv_states, var_list, os_data, libname, model, i_core, is_nlo)
      class(prc_recola_t), intent(inout) :: core
      integer, intent(in), dimension(:,:), allocatable :: flv_states
      type(var_list_t), intent(in) :: var_list
      type(os_data_t), intent(in) :: os_data
      type(string_t), intent(in) :: libname
      type(model_data_t), intent(in), target :: model
      integer, intent(in) :: i_core
      logical, intent(in) :: is_nlo
    end subroutine prc_recola_prepare_external_code
    module subroutine prc_recola_set_parameters (object, qcd, model)
      class(prc_recola_t), intent(inout) :: object
      type(qcd_t), intent(in) :: qcd
      class(model_data_t), intent(in), target, optional :: model
    end subroutine prc_recola_set_parameters
    module subroutine prc_recola_init (object, def, lib, id, i_component)
      class(prc_recola_t), intent(inout) :: object
      class(prc_core_def_t), intent(in), target :: def
      type(process_library_t), intent(in), target :: lib
      type(string_t), intent(in) :: id
      integer, intent(in) :: i_component
    end subroutine prc_recola_init
    module subroutine prc_recola_replace_helicity_and_color_arrays (object)
      class(prc_recola_t), intent(inout) :: object
    end subroutine prc_recola_replace_helicity_and_color_arrays
    module function prc_recola_compute_amplitude &
       (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
       core_state) result (amp)
      complex(default) :: amp
      class(prc_recola_t), intent(in) :: object
      integer, intent(in) :: j
      type(vector4_t), intent(in), dimension(:) :: p
      integer, intent(in) :: f, h, c
      real(default), intent(in) :: fac_scale, ren_scale
      real(default), intent(in), allocatable :: alpha_qcd_forced
      class(prc_core_state_t), intent(inout), allocatable, optional :: &
           core_state
    end function prc_recola_compute_amplitude
    module subroutine prc_recola_compute_sqme (object, i_flv, i_hel, p, &
         ren_scale, sqme, bad_point)
      class(prc_recola_t), intent(in) :: object
      integer, intent(in) :: i_flv, i_hel
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: ren_scale
      real(default), intent(out) :: sqme
      logical, intent(out) :: bad_point
    end subroutine prc_recola_compute_sqme
    module subroutine prc_recola_compute_sqme_virt (object, i_flv, i_hel,  &
         p, ren_scale, es_scale, loop_method, sqme, bad_point)
      class(prc_recola_t), intent(in) :: object
      integer, intent(in) :: i_flv, i_hel
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: ren_scale, es_scale
      integer, intent(in) :: loop_method
      real(default), dimension(4), intent(out) :: sqme
      real(default) :: amp
      logical, intent(out) :: bad_point
    end subroutine prc_recola_compute_sqme_virt
    module function prc_recola_get_alpha_qed &
         (object, core_state) result (alpha_qed)
      class(prc_recola_t), intent(in) :: object
      class(prc_core_state_t), intent(in), allocatable :: core_state
      real(default) :: alpha_qed
    end function prc_recola_get_alpha_qed
    module subroutine prc_recola_compute_sqme_color_c_raw (object, &
         i_flv, i_hel, p, ren_scale, sqme_color_c, bad_point)
      class(prc_recola_t), intent(in) :: object
      integer, intent(in) :: i_hel, i_flv
      type(vector4_t), dimension(:), intent(in) :: p
      real(default), intent(in) :: ren_scale
      real(default), dimension(:), intent(out) :: sqme_color_c
      logical, intent(out) :: bad_point
    end subroutine prc_recola_compute_sqme_color_c_raw
  end interface

contains

  subroutine recola_def_init (object, basename, model_name, &
       prt_in, prt_out, nlo_type, alpha_power, alphas_power, &
       correction_type, restrictions)
    class(recola_def_t), intent(inout) :: object
    type(string_t), intent(in) :: basename, model_name
    type(string_t), dimension(:), intent(in) :: prt_in, prt_out
    integer, intent(in) :: nlo_type
    integer, intent(in) :: alpha_power
    integer, intent(in) :: alphas_power
    type(string_t), intent(in) :: correction_type
    type(string_t), intent(in), optional :: restrictions
    if (debug_on) call msg_debug (D_ME_METHODS, "recola_def_init: " &
         // char (basename) // ", nlo_type", nlo_type)
    object%basename = basename
    object%alpha_power = alpha_power
    object%alphas_power = alphas_power
    select case (char (correction_type))
    case ("QCD")
       object%corr = RECOLA_QCD
    case ("EW")
       object%corr = RECOLA_EW
    case ("Full")
       object%corr = RECOLA_FULL
    end select
    allocate (recola_writer_t :: object%writer)
    select case (nlo_type)
    case (BORN)
       object%suffix = '_BORN'
       object%order = "LO"
    case (NLO_REAL)
       object%suffix = '_REAL'
       object%order = "LO"
       if (object%corr == RECOLA_QCD)  object%alphas_power = alphas_power + 1
       if (object%corr == RECOLA_EW)  object%alpha_power = alpha_power + 1
    case (NLO_VIRTUAL)
       object%suffix = '_LOOP'
       object%order = "NLO"
    case (NLO_SUBTRACTION)
       object%suffix = '_SUB'
       object%order = "LO"
    case (NLO_MISMATCH)
       object%suffix = '_MISMATCH'
       object%order = "LO"
    case (NLO_DGLAP)
       object%suffix = '_DGLAP'
       object%order = "LO"
    end select
    select type (writer => object%writer)
    class is (recola_writer_t)
       call writer%init (model_name, prt_in, prt_out, restrictions)
       call writer%set_id (basename // object%suffix)
       call writer%set_order (object%order)
       call writer%set_coupling_powers (object%alpha_power, object%alphas_power)
    end select
  end subroutine recola_def_init

  subroutine recola_def_allocate_driver (object, driver, basename)
    class(recola_def_t), intent(in) :: object
    class(prc_core_driver_t), intent(out), allocatable :: driver
    type(string_t), intent(in) :: basename
    if (debug_on) call msg_debug2 (D_ME_METHODS, "recola_def_allocate_driver")
    allocate (recola_driver_t :: driver)
  end subroutine recola_def_allocate_driver

  subroutine prc_recola_allocate_workspace (object, core_state)
    class(prc_recola_t), intent(in) :: object
    class(prc_core_state_t), intent(inout), allocatable :: core_state
    allocate (recola_state_t :: core_state)
  end subroutine prc_recola_allocate_workspace


end module prc_recola
