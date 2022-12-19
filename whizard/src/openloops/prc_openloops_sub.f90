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

submodule (prc_openloops) prc_openloops_s

  use debug_master, only: debug_on
  use io_units
  use constants
  use numeric_utils
  use diagnostics
  use system_dependencies
  use sm_physics, only: top_width_sm_lo, top_width_sm_qcd_nlo_jk
  use sm_qcd
  use prclib_interfaces

  implicit none

  real(default), parameter :: openloops_default_bmass = 0._default
  real(default), parameter :: openloops_default_topmass = 172._default
  real(default), parameter :: openloops_default_topwidth = 0._default
  real(default), parameter :: openloops_default_wmass = 80.399_default
  real(default), parameter :: openloops_default_wwidth = 0._default
  real(default), parameter :: openloops_default_zmass = 91.1876_default
  real(default), parameter :: openloops_default_zwidth = 0._default
  real(default), parameter :: openloops_default_higgsmass = 125._default
  real(default), parameter :: openloops_default_higgswidth = 0._default


contains

  module function openloops_threshold_data_compute_top_width &
       (data, mtop, alpha_s) result (wtop)
    real(default) :: wtop
    class(openloops_threshold_data_t), intent(in) :: data
    real(default), intent(in) :: mtop, alpha_s
    if (data%nlo) then
       wtop = top_width_sm_qcd_nlo_jk (data%alpha_ew, data%sinthw, &
              data%vtb, mtop, data%m_W, data%m_b, alpha_s)
    else
       wtop = top_width_sm_lo (data%alpha_ew, data%sinthw, data%vtb, &
              mtop, data%m_W, data%m_b)
    end if
  end function openloops_threshold_data_compute_top_width

  module subroutine openloops_state_init_threshold (object, model)
    class(openloops_state_t), intent(inout) :: object
    type(model_data_t), intent(in) :: model
    if (model%get_name () == "SM_tt_threshold") then
       allocate (object%threshold_data)
       associate (data => object%threshold_data)
          data%nlo = btest (int (model%get_real (var_str ('offshell_strategy'))), 0)
          data%alpha_ew = one / model%get_real (var_str ('alpha_em_i'))
          data%sinthw = model%get_real (var_str ('sw'))
          data%m_b = model%get_real (var_str ('mb'))
          data%m_W = model%get_real (var_str ('mW'))
          data%vtb = model%get_real (var_str ('Vtb'))
       end associate
    end if
  end subroutine openloops_state_init_threshold

  module function openloops_writer_type_name () result (string)
    type(string_t) :: string
    string = "openloops"
  end function openloops_writer_type_name

  module function openloops_def_type_string () result (string)
    type(string_t) :: string
    string = "openloops"
  end function openloops_def_type_string

  module subroutine openloops_def_write (object, unit)
    class(openloops_def_t), intent(in) :: object
    integer, intent(in) :: unit
    select type (writer => object%writer)
    type is (openloops_writer_t)
       call writer%write (unit)
    end select
  end subroutine openloops_def_write

  module subroutine openloops_driver_init_dlaccess_to_library &
     (object, os_data, dlaccess, success)
    class(openloops_driver_t), intent(in) :: object
    type(os_data_t), intent(in) :: os_data
    type(dlaccess_t), intent(out) :: dlaccess
    logical, intent(out) :: success
    type(string_t) :: ol_library, msg_buffer
    ol_library = OPENLOOPS_DIR // '/lib/libopenloops.' // &
         os_data%shrlib_ext
    msg_buffer = "One-Loop-Provider: Using OpenLoops"
    call msg_message (char(msg_buffer))
    msg_buffer = "Loading library: " // ol_library
    call msg_message (char(msg_buffer))
    if (os_file_exist (ol_library)) then
       call dlaccess_init (dlaccess, var_str (""), ol_library, os_data)
    else
       call msg_fatal ("Link OpenLoops: library not found")
    end if
    success = .not. dlaccess_has_error (dlaccess)
  end subroutine openloops_driver_init_dlaccess_to_library

  module subroutine openloops_driver_set_alpha_s (driver, alpha_s)
    class(openloops_driver_t), intent(in) :: driver
    real(default), intent(in) :: alpha_s
    integer :: ierr
    if (associated (driver%blha_olp_set_parameter)) then
       call driver%blha_olp_set_parameter &
            (c_char_'alphas'//c_null_char, &
             dble (alpha_s), 0._double, ierr)
    else
       call msg_fatal ("blha_olp_set_parameter not associated!")
    end if
    if (ierr == 0) call parameter_error_message (var_str ('alphas'), &
         var_str ('openloops_driver_set_alpha_s'))
  end subroutine openloops_driver_set_alpha_s

  module subroutine openloops_driver_set_alpha_qed (driver, alpha)
    class(openloops_driver_t), intent(inout) :: driver
    real(default), intent(in) :: alpha
    integer :: ierr
    call driver%blha_olp_set_parameter &
       (c_char_'alpha_qed'//c_null_char, &
        dble (alpha), 0._double, ierr)
    if (ierr == 0) call ew_parameter_error_message (var_str ('alpha_qed'))
  end subroutine openloops_driver_set_alpha_qed

  module subroutine openloops_driver_set_GF (driver, GF)
    class(openloops_driver_t), intent(inout) :: driver
    real(default), intent(in) :: GF
    integer :: ierr
    call driver%blha_olp_set_parameter &
       (c_char_'Gmu'//c_null_char, &
        dble(GF), 0._double, ierr)
    if (ierr == 0) call ew_parameter_error_message (var_str ('Gmu'))
  end subroutine openloops_driver_set_GF

  module subroutine openloops_driver_set_weinberg_angle (driver, sw2)
    class(openloops_driver_t), intent(inout) :: driver
    real(default), intent(in) :: sw2
    integer :: ierr
    call driver%blha_olp_set_parameter &
       (c_char_'sw2'//c_null_char, &
        dble(sw2), 0._double, ierr)
    if (ierr == 0) call ew_parameter_error_message (var_str ('sw2'))
  end subroutine openloops_driver_set_weinberg_angle

  module subroutine openloops_driver_print_alpha_s (object)
    class(openloops_driver_t), intent(in) :: object
    call object%blha_olp_print_parameter (c_char_'alphas'//c_null_char)
  end subroutine openloops_driver_print_alpha_s

  module function openloops_driver_type_name () result (type)
    type(string_t) :: type
    type = "OpenLoops"
  end function openloops_driver_type_name

  module subroutine openloops_driver_load_procedures &
       (object, os_data, success)
    class(openloops_driver_t), intent(inout) :: object
    type(os_data_t), intent(in) :: os_data
    logical, intent(out) :: success
    type(dlaccess_t) :: dlaccess
    type(c_funptr) :: c_fptr
    logical :: init_success

    call object%init_dlaccess_to_library (os_data, dlaccess, init_success)

    c_fptr = dlaccess_get_c_funptr (dlaccess, var_str ("ol_evaluate_scpowheg"))
    call c_f_procpointer (c_fptr, object%evaluate_spin_correlations_powheg)
    call check_for_error (var_str ("ol_evaluate_scpowheg"))
    c_fptr = dlaccess_get_c_funptr (dlaccess, var_str ("ol_getparameter_double"))
    call c_f_procpointer (c_fptr, object%get_parameter_double)
    call check_for_error (var_str ("ol_getparameter_double"))
    success = .true.

    contains
      subroutine check_for_error (function_name)
        type(string_t), intent(in) :: function_name
        if (dlaccess_has_error (dlaccess)) &
           call msg_fatal (char ("Loading of " // function_name // " failed!"))
     end subroutine check_for_error
  end subroutine openloops_driver_load_procedures

  module subroutine openloops_def_read (object, unit)
    class(openloops_def_t), intent(out) :: object
    integer, intent(in) :: unit
  end subroutine openloops_def_read

  module subroutine openloops_state_write (object, unit)
    class(openloops_state_t), intent(in) :: object
    integer, intent(in), optional :: unit
  end subroutine openloops_state_write

  module subroutine prc_openloops_init_driver (object, os_data)
    class(prc_openloops_t), intent(inout) :: object
    type(os_data_t), intent(in) :: os_data
    type(string_t) :: olp_file, olc_file
    type(string_t) :: suffix

    select type (def => object%def)
    type is (openloops_def_t)
       suffix = def%suffix
       olp_file = def%basename // suffix // '.olp'
       olc_file = def%basename // suffix // '.olc'
    class default
       call msg_bug ("prc_openloops_init_driver: core_def should be openloops-type")
    end select

    select type (driver => object%driver)
    type is (openloops_driver_t)
       driver%olp_file = olp_file
       driver%contract_file = olc_file
       driver%nlo_suffix = suffix
    end select
  end subroutine prc_openloops_init_driver

  module subroutine prc_openloops_write (object, unit)
    class(prc_openloops_t), intent(in) :: object
    integer, intent(in), optional :: unit
    call msg_message (unit = unit, string = "OpenLoops")
  end subroutine prc_openloops_write

  module subroutine prc_openloops_write_name (object, unit)
    class(prc_openloops_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u,"(1x,A)") "Core: OpenLoops"
  end subroutine prc_openloops_write_name

  module subroutine prc_openloops_prepare_library (object, os_data, model)
    class(prc_openloops_t), intent(inout) :: object
    type(os_data_t), intent(in) :: os_data
    type(model_data_t), intent(in), target :: model
    call object%load_driver (os_data)
    call object%reset_parameters ()
    call object%set_particle_properties (model)
    select type(def => object%def)
    type is (openloops_def_t)
       call object%set_verbosity (def%verbosity)
    end select
  end subroutine prc_openloops_prepare_library

  module subroutine prc_openloops_load_driver (object, os_data)
    class(prc_openloops_t), intent(inout) :: object
    type(os_data_t), intent(in) :: os_data
    logical :: success
    select type (driver => object%driver)
    type is (openloops_driver_t)
       call driver%load (os_data, success)
       call driver%load_procedures (os_data, success)
    end select
  end subroutine prc_openloops_load_driver

  module subroutine prc_openloops_start (object)
    class(prc_openloops_t), intent(inout) :: object
    integer :: ierr
    select type (driver => object%driver)
    type is (openloops_driver_t)
       call driver%blha_olp_start (char (driver%olp_file)//c_null_char, ierr)
    end select
  end subroutine prc_openloops_start

  module subroutine prc_openloops_set_n_external (object, n)
    class(prc_openloops_t), intent(inout) :: object
    integer, intent(in) :: n
    N_EXTERNAL = n
  end subroutine prc_openloops_set_n_external

  module subroutine prc_openloops_reset_parameters (object)
    class(prc_openloops_t), intent(inout) :: object
    integer :: ierr
    select type (driver => object%driver)
    type is (openloops_driver_t)
       call driver%blha_olp_set_parameter ('mass(5)'//c_null_char, &
            dble(openloops_default_bmass), 0._double, ierr)
       call driver%blha_olp_set_parameter ('mass(6)'//c_null_char, &
            dble(openloops_default_topmass), 0._double, ierr)
       call driver%blha_olp_set_parameter ('width(6)'//c_null_char, &
            dble(openloops_default_topwidth), 0._double, ierr)
       call driver%blha_olp_set_parameter ('mass(23)'//c_null_char, &
            dble(openloops_default_zmass), 0._double, ierr)
       call driver%blha_olp_set_parameter ('width(23)'//c_null_char, &
            dble(openloops_default_zwidth), 0._double, ierr)
       call driver%blha_olp_set_parameter ('mass(24)'//c_null_char, &
            dble(openloops_default_wmass), 0._double, ierr)
       call driver%blha_olp_set_parameter ('width(24)'//c_null_char, &
            dble(openloops_default_wwidth), 0._double, ierr)
       call driver%blha_olp_set_parameter ('mass(25)'//c_null_char, &
            dble(openloops_default_higgsmass), 0._double, ierr)
       call driver%blha_olp_set_parameter ('width(25)'//c_null_char, &
            dble(openloops_default_higgswidth), 0._double, ierr)
    end select
  end subroutine prc_openloops_reset_parameters

  module subroutine prc_openloops_set_verbosity (object, verbose)
    class(prc_openloops_t), intent(inout) :: object
    integer, intent(in) :: verbose
    integer :: ierr
    select type (driver => object%driver)
    type is (openloops_driver_t)
       call driver%blha_olp_set_parameter ('verbose'//c_null_char, &
            dble(verbose), 0._double, ierr)
    end select
  end subroutine prc_openloops_set_verbosity

  module subroutine prc_openloops_prepare_external_code &
       (core, flv_states, var_list, os_data, libname, model, i_core, is_nlo)
    class(prc_openloops_t), intent(inout) :: core
    integer, intent(in), dimension(:,:), allocatable :: flv_states
    type(var_list_t), intent(in) :: var_list
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in) :: libname
    type(model_data_t), intent(in), target :: model
    integer, intent(in) :: i_core
    logical, intent(in) :: is_nlo
    integer :: ierr
    core%sqme_tree_pos = 1
    call core%set_n_external (core%data%get_n_tot ())
    call core%prepare_library (os_data, model)
    call core%start ()
    call core%set_electroweak_parameters (model)
    select type (driver => core%driver)
    type is (openloops_driver_t)
       !!! We have to set the external vector boson wavefunction to the MadGraph convention
       !!! in order to be consistent with our calculation of the spin correlated contributions.
       call driver%blha_olp_set_parameter ('wf_v_select'//c_null_char, &
            3._double, 0._double, ierr)
       if (ierr == 0) call parameter_error_message (var_str ('wf_v_select'), &
            var_str ('prc_openloops_prepare_external_code'))
    end select
    call core%read_contract_file (flv_states)
    call core%print_parameter_file (i_core)
  end subroutine prc_openloops_prepare_external_code

  module subroutine prc_openloops_compute_sqme_spin_c (object, &
       i_flv, i_hel, em, p, ren_scale, sqme_spin_c, bad_point)
    class(prc_openloops_t), intent(inout) :: object
    integer, intent(in) :: i_flv, i_hel
    integer, intent(in) :: em
    type(vector4_t), intent(in), dimension(:) :: p
    real(default), intent(in) :: ren_scale
    real(default), intent(out), dimension(6) :: sqme_spin_c
    logical, intent(out) :: bad_point
    real(default), dimension(16) :: sqme_spin_c_tmp
    real(double), dimension(5*N_EXTERNAL) :: mom
    real(double) :: res
    real(double), dimension(16) :: res_munu
    real(default) :: alpha_s
    if (object%i_spin_c(i_flv, i_hel) >= 0) then
       mom = object%create_momentum_array (p)
       if (vanishes (ren_scale)) call msg_fatal &
            ("prc_openloops_compute_sqme_spin_c: ren_scale vanishes")
       alpha_s = object%qcd%alpha%get (ren_scale)
       select type (driver => object%driver)
       type is (openloops_driver_t)
          call driver%set_alpha_s (alpha_s)
          call driver%evaluate_spin_correlations_powheg &
               (object%i_spin_c(i_flv, i_hel), mom, em, res, res_munu)
       end select
       sqme_spin_c_tmp = res_munu
       bad_point = .false.
       if (object%includes_polarization ()) &
            sqme_spin_c_tmp = object%n_hel * sqme_spin_c_tmp
       if (debug_on) then
          if (sum(sqme_spin_c_tmp) == 0) then
             call msg_debug(D_SUBTRACTION,'Spin-correlated matrix elements provided by OpenLoops are zero!')
          end if
       end if
    else
       sqme_spin_c_tmp = zero
    end if
    !!! Using symmetry of the 4x4 matrix of spin correlated squared Born MEs and
    !!! the fact that we multiply only with vectors with E=0. We thus store the
    !!! upper triangle of the lower 3x3 matrix as columns in a 1-dim array
    sqme_spin_c(1:2) = sqme_spin_c_tmp(6:7)
    sqme_spin_c(3) = sqme_spin_c_tmp(11)
    sqme_spin_c(4) = sqme_spin_c_tmp(8)
    sqme_spin_c(5) = sqme_spin_c_tmp(12)
    sqme_spin_c(6) = sqme_spin_c_tmp(16)
  end subroutine prc_openloops_compute_sqme_spin_c

  module function prc_openloops_get_alpha_qed &
       (object, core_state) result (alpha_qed)
    class(prc_openloops_t), intent(in) :: object
    class(prc_core_state_t), intent(in), allocatable :: core_state
    real(double) :: value
    real(default) :: alpha_qed
    select type (driver => object%driver)
    type is (openloops_driver_t)
       call driver%get_parameter_double ('alpha_qed'//c_null_char, value)
       alpha_qed = value
       return
    end select
    call msg_fatal ("prc_openloops_get_alpha_qed: " // &
         "called by wrong driver, only supported for OpenLoops!")
  end function prc_openloops_get_alpha_qed


end submodule prc_openloops_s

