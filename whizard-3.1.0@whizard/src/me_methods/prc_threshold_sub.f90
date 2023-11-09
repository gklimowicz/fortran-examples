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

submodule (prc_threshold) prc_threshold_s

  use constants
  use numeric_utils
  use string_utils, only: lower_case
  use io_units
  use system_defs, only: TAB
  use sm_qcd

  implicit none

contains

  pure module subroutine threshold_writer_init &
       (writer, model_name, prt_in, prt_out, restrictions)
    class(threshold_writer_t), intent(inout) :: writer
    type(string_t), intent(in) :: model_name
    type(string_t), dimension(:), intent(in) :: prt_in, prt_out
    type(string_t), intent(in), optional :: restrictions
    call writer%base_init (model_name, prt_in, prt_out, restrictions)
    writer%amp_triv = .false.
  end subroutine threshold_writer_init

  module subroutine threshold_writer_write_makefile_extra &
       (writer, unit, id, os_data, verbose, nlo_type)
    class(threshold_writer_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    integer, intent(in) :: nlo_type
    type(string_t) :: f90in, f90, lo, extra
    if (debug_on) call msg_debug &
         (D_ME_METHODS, "threshold_writer_write_makefile_extra")
    if (nlo_type /= BORN) then
       extra = "_" // component_status (nlo_type)
    else
       extra = var_str ("")
    end if
    f90 = id // "_threshold" // extra //".f90"
    f90in = f90 // ".in"
    lo = id // "_threshold" // extra // ".lo"
    write (unit, "(A)") "OBJECTS += " // char (lo)
    write (unit, "(A)") char (f90in) // ":"
    write (unit, "(A)") char (TAB // "if ! test -f " // f90in // &
         "; then cp " // os_data%whizard_sharepath // &
         "/SM_tt_threshold_data/threshold" // extra // ".f90 " // &
         f90in // "; fi")
    write (unit, "(A)") char(f90) // ": " // char (f90in)
    write (unit, "(A)") TAB // "sed 's/@ID@/" // char (id) // "/' " // &
         char (f90in) // " > " // char (f90)
    write (unit, "(5A)")  "CLEAN_SOURCES += ", char (f90)
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (f90in)
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), "_threshold.mod"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (lo)
    write (unit, "(A)") char(lo) // ": " // char (f90) // " " // &
         char(id) // ".f90"
    write (unit, "(5A)")  TAB, "$(LTFCOMPILE) $<"
    if (.not. verbose) then
       write (unit, "(5A)")  TAB // '@echo  "  FC       " $@'
    end if
  end subroutine threshold_writer_write_makefile_extra

  module subroutine threshold_writer_write_makefile_code &
       (writer, unit, id, os_data, verbose, testflag)
    class(threshold_writer_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    if (debug_on) &
         call msg_debug (D_ME_METHODS, "threshold_writer_write_makefile_code")
    call writer%base_write_makefile_code &
         (unit, id, os_data, verbose, testflag = testflag)
    call writer%write_makefile_extra (unit, id, os_data, verbose, BORN)
    if (writer%nlo_type == NLO_VIRTUAL .and. writer%active) &
         call writer%write_makefile_extra (unit, id, os_data, verbose, writer%nlo_type)
  end subroutine threshold_writer_write_makefile_code

  module function threshold_writer_type_name () result (string)
    type(string_t) :: string
    string = "Threshold"
  end function threshold_writer_type_name

  module function threshold_driver_type_name () result (type)
    type(string_t) :: type
    type = "Threshold"
  end function threshold_driver_type_name

  module subroutine threshold_driver_load (threshold_driver, dlaccess)
    class(threshold_driver_t), intent(inout) :: threshold_driver
    type(dlaccess_t), intent(inout) :: dlaccess
    type(c_funptr) :: c_fptr
    type(string_t) :: lower_case_id
    if (debug_on) call msg_debug (D_ME_METHODS, "threshold_driver_load")
    lower_case_id = lower_case (threshold_driver%id)
    c_fptr = dlaccess_get_c_funptr (dlaccess, lower_case_id // "_set_process_mode")
    call c_f_procpointer (c_fptr, threshold_driver%set_process_mode)
    call check_for_error (lower_case_id // "_set_process_mode")
    c_fptr = dlaccess_get_c_funptr (dlaccess, lower_case_id // "_get_amp_squared")
    call c_f_procpointer (c_fptr, threshold_driver%get_amp_squared)
    call check_for_error (lower_case_id // "_get_amp_squared")
    c_fptr = dlaccess_get_c_funptr (dlaccess, lower_case_id // "_threshold_init")
    call c_f_procpointer (c_fptr, threshold_driver%init)
    call check_for_error (lower_case_id // "_threshold_init")
    select type (threshold_driver)
    type is (threshold_driver_t)
       if (threshold_driver%nlo_type == NLO_VIRTUAL) then
          c_fptr = dlaccess_get_c_funptr &
               (dlaccess, lower_case_id // "_start_openloops")
          call c_f_procpointer (c_fptr, threshold_driver%start_openloops)
          call check_for_error (lower_case_id // "_start_openloops")
          c_fptr = dlaccess_get_c_funptr (dlaccess, lower_case_id // "_olp_eval2")
          call c_f_procpointer (c_fptr, threshold_driver%olp_eval2)
          call check_for_error (lower_case_id // "_olp_eval2")
       end if
    end select
    call msg_message ("Loaded extra threshold functions")
    contains
      subroutine check_for_error (function_name)
        type(string_t), intent(in) :: function_name
        if (dlaccess_has_error (dlaccess))  call msg_fatal &
             (char ("Loading of " // function_name // " failed!"))
     end subroutine check_for_error
  end subroutine threshold_driver_load

  module function threshold_def_type_string () result (string)
    type(string_t) :: string
    string = "threshold computation"
  end function threshold_def_type_string

  module subroutine threshold_def_write (object, unit)
    class(threshold_def_t), intent(in) :: object
    integer, intent(in) :: unit
  end subroutine threshold_def_write

  module subroutine threshold_def_read (object, unit)
    class(threshold_def_t), intent(out) :: object
    integer, intent(in) :: unit
  end subroutine threshold_def_read

  module subroutine threshold_def_connect (def, lib_driver, i, proc_driver)
    class(threshold_def_t), intent(in) :: def
    class(prclib_driver_t), intent(in) :: lib_driver
    integer, intent(in) :: i
    class(prc_core_driver_t), intent(inout) :: proc_driver
    type(dlaccess_t) :: dlaccess
    logical :: skip
    if (debug_on) call msg_debug (D_ME_METHODS, "threshold_def_connect")
    call def%omega_connect (lib_driver, i, proc_driver)
    select type (lib_driver)
    class is (prclib_driver_dynamic_t)
       dlaccess = lib_driver%dlaccess
    end select
    select type (proc_driver)
    class is (threshold_driver_t)
       select type (writer => def%writer)
       type is (threshold_writer_t)
          skip = writer%nlo_type == NLO_VIRTUAL .and. .not. writer%active
          if (.not. skip) call proc_driver%load (dlaccess)
       end select
    end select
  end subroutine threshold_def_connect

  module subroutine threshold_state_write (object, unit)
    class(threshold_state_t), intent(in) :: object
    integer, intent(in), optional :: unit
  end subroutine threshold_state_write

  module subroutine prc_threshold_write (object, unit)
    class(prc_threshold_t), intent(in) :: object
    integer, intent(in), optional :: unit
    call msg_message ("Supply amplitudes squared for threshold computation")
  end subroutine prc_threshold_write

  module subroutine prc_threshold_write_name (object, unit)
    class(prc_threshold_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u,"(1x,A)") "Core: Threshold"
  end subroutine prc_threshold_write_name

  module subroutine prc_threshold_set_beam_pol (object, has_beam_pol)
    class(prc_threshold_t), intent(inout) :: object
    logical, intent(in), optional :: has_beam_pol
    if (present (has_beam_pol)) then
       object%has_beam_pol = has_beam_pol
    end if
  end subroutine prc_threshold_set_beam_pol

  module function prc_threshold_compute_amplitude &
       (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
       core_state)  result (amp)
    class(prc_threshold_t), intent(in) :: object
    integer, intent(in) :: j
    type(vector4_t), dimension(:), intent(in) :: p
    integer, intent(in) :: f, h, c
    real(default), intent(in) :: fac_scale, ren_scale
    real(default), intent(in), allocatable :: alpha_qcd_forced
    class(prc_core_state_t), intent(inout), allocatable, optional :: core_state
    complex(default) :: amp
    select type (core_state)
    class is (prc_external_test_state_t)
       core_state%alpha_qcd = object%qcd%alpha%get (ren_scale)
    end select
    amp = 0
  end function prc_threshold_compute_amplitude

  module subroutine prc_threshold_set_offshell_momenta (object, p)
    class(prc_threshold_t), intent(inout) :: object
    type(vector4_t), intent(in), dimension(:) :: p
    integer :: i
    do i = 1, size(p)
       object%parray_ofs(:,i) = p(i)%p
    end do
  end subroutine prc_threshold_set_offshell_momenta

  module subroutine prc_threshold_set_onshell_momenta (object, p)
    class(prc_threshold_t), intent(inout) :: object
    type(vector4_t), intent(in), dimension(:) :: p
    integer :: i
    do  i = 1, size(p)
       object%parray_ons(:,i) = p(i)%p
    end do
  end subroutine prc_threshold_set_onshell_momenta

  module subroutine prc_threshold_set_leg (object, leg)
    class(prc_threshold_t), intent(inout) :: object
    integer, intent(in) :: leg
    object%leg = leg
  end subroutine prc_threshold_set_leg

  module subroutine prc_threshold_set_process_mode (object, mode)
    class(prc_threshold_t), intent(in) :: object
    integer(kind = c_int), intent(in) :: mode
    select type (driver => object%driver)
    class is (threshold_driver_t)
       if (associated (driver%set_process_mode)) &
            call driver%set_process_mode (mode)
    end select
  end subroutine prc_threshold_set_process_mode

  module subroutine prc_threshold_compute_sqme (object, i_flv, i_hel, p, &
         ren_scale, sqme, bad_point)
    class(prc_threshold_t), intent(in) :: object
    integer, intent(in) :: i_flv, i_hel
    type(vector4_t), intent(in), dimension(:) :: p
    real(default), intent(in) :: ren_scale
    real(default), intent(out) :: sqme
    logical, intent(out) :: bad_point
    integer :: n_tot
    if (debug_on) call msg_debug2 (D_ME_METHODS, "prc_threshold_compute_sqme")
    n_tot = size (p)
    select type (driver => object%driver)
    class is (threshold_driver_t)
       if (object%has_beam_pol) then
          call driver%get_amp_squared (sqme, object%parray_ofs, &
               object%parray_ons, object%leg, n_tot, i_flv - 1)
       else
          call driver%get_amp_squared (sqme, object%parray_ofs, &
               object%parray_ons, object%leg, n_tot, -1)
       end if
    end select
    bad_point = .false.
  end subroutine prc_threshold_compute_sqme

  module subroutine prc_threshold_compute_sqme_virt (object, i_flv, i_hel, &
         p, ren_scale, es_scale, loop_method, sqme, bad_point)
    class(prc_threshold_t), intent(in) :: object
    integer, intent(in) :: i_flv, i_hel
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(in) :: ren_scale, es_scale
    integer, intent(in) :: loop_method
    real(default), dimension(4), intent(out) :: sqme
    real(c_default_float), dimension(:,:), allocatable, save :: parray
    logical, intent(out) :: bad_point
    integer :: n_tot, i
    real(default) :: mu
    real(c_default_float), dimension(4) :: sqme_c
    real(c_default_float) :: mu_c, acc_c, alpha_s_c
    integer(c_int) :: i_flv_c
    if (debug_on) call msg_debug2 (D_ME_METHODS, "prc_threshold_compute_sqme_virt")
    n_tot = size (p)
    if (allocated (parray)) then
       if (size(parray) /= n_tot) deallocate (parray)
    end if
    if (.not. allocated (parray))  allocate (parray (0:3, n_tot))
    forall (i = 1:n_tot)  parray(:,i) = p(i)%p

    if (vanishes (ren_scale)) then
      mu = sqrt (2* (p(1)*p(2)))
    else
      mu = ren_scale
    end if
    mu_c = mu
    alpha_s_c = object%qcd%alpha%get (mu)
    i_flv_c = i_flv
    select type (driver => object%driver)
    class is (threshold_driver_t)
       if (associated (driver%olp_eval2)) then
          if (object%has_beam_pol) then
             call driver%olp_eval2 (1, alpha_s_c, &
                  parray, mu_c, i_flv_c - 1, sqme_c, acc_c)
          else
             call driver%olp_eval2 (i_flv_c, alpha_s_c, &
                  parray, mu_c, -1, sqme_c, acc_c)
          end if
          bad_point = real(acc_c, kind=default) > object%maximum_accuracy
          sqme = sqme_c
       else
          sqme = 0._default
          bad_point = .true.
       end if
    end select
  end subroutine prc_threshold_compute_sqme_virt

  module subroutine prc_threshold_init (object, def, lib, id, i_component)
    class(prc_threshold_t), intent(inout) :: object
    class(prc_core_def_t), intent(in), target :: def
    type(process_library_t), intent(in), target :: lib
    type(string_t), intent(in) :: id
    integer, intent(in) :: i_component
    integer :: n_tot
    call object%base_init (def, lib, id, i_component)
    n_tot = object%data%n_in + object%data%n_out
    allocate (object%parray_ofs (0:3,n_tot), object%parray_ons (0:3,n_tot))
    if (n_tot == 4) then
       call object%set_process_mode (PROC_MODE_TT)
    else
       call object%set_process_mode (PROC_MODE_WBWB)
    end if
    call object%activate_parameters ()
  end subroutine prc_threshold_init

  module subroutine prc_threshold_activate_parameters (object)
    class (prc_threshold_t), intent(inout) :: object
    if (debug_on) &
         call msg_debug (D_ME_METHODS, "prc_threshold_activate_parameters")
    if (allocated (object%driver)) then
       if (allocated (object%par)) then
          select type (driver => object%driver)
          type is (threshold_driver_t)
             if (associated (driver%init)) then
                call driver%init (object%par, object%scheme)
             end if
          end select
       else
          call msg_bug &
               ("prc_threshold_activate: parameter set is not allocated")
       end if
    else
       call msg_bug ("prc_threshold_activate: driver is not allocated")
    end if
  end subroutine prc_threshold_activate_parameters

  module subroutine prc_threshold_prepare_external_code &
       (core, flv_states, var_list, os_data, libname, model, i_core, is_nlo)
    class(prc_threshold_t), intent(inout) :: core
    integer, intent(in), dimension(:,:), allocatable :: flv_states
    type(var_list_t), intent(in) :: var_list
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in) :: libname
    type(model_data_t), intent(in), target :: model
    integer, intent(in) :: i_core
    logical, intent(in) :: is_nlo
    if (debug_on) call msg_debug (D_ME_METHODS, &
         "prc_threshold_prepare_external_code")
    if (allocated (core%driver)) then
       select type (driver => core%driver)
       type is (threshold_driver_t)
          if (driver%nlo_type == NLO_VIRTUAL) call driver%start_openloops ()
       end select
    else
       call msg_bug ("prc_threshold_prepare_external_code: " &
            // "driver is not allocated")
    end if
  end subroutine prc_threshold_prepare_external_code

  module function prc_threshold_includes_polarization &
       (object) result (polarized)
    logical :: polarized
    class(prc_threshold_t), intent(in) :: object
    polarized = object%has_beam_pol
  end function prc_threshold_includes_polarization


end submodule prc_threshold_s

