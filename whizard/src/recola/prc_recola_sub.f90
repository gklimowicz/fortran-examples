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

submodule (prc_recola) prc_recola_s

  use constants, only: pi, zero
  use string_utils, only: str
  use system_defs, only: TAB
  use io_units
  use recola_wrapper !NODEP!

  implicit none

contains

  module subroutine abort_if_recola_not_active ()
    if (.not. rclwrap_is_active) call msg_fatal ("You want to use Recola, ", &
       [var_str("but either the compiler with which Whizard has been build "), &
        var_str("is not supported by it, or you have not linked Recola "), &
        var_str("correctly to Whizard. Either reconfigure Whizard with a path to "), &
        var_str("a valid Recola installation (for details consult the manual), "), &
        var_str("or choose a different matrix-element method.")])
  end subroutine abort_if_recola_not_active

  module function recola_def_type_string () result (string)
    type(string_t) :: string
    string = "recola"
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
    string = "recola"
  end function recola_writer_type_name

  module subroutine recola_writer_set_id (writer, id)
    class(recola_writer_t), intent(inout) :: writer
    type(string_t), intent(in) :: id
    if (debug_on)  call msg_debug2 &
         (D_ME_METHODS, "Recola writer: id = " // char (id))
    writer%id = id
  end subroutine recola_writer_set_id

  module subroutine recola_writer_set_order (writer, order)
    class(recola_writer_t), intent(inout) :: writer
    type(string_t), intent(in) :: order
    if (debug_on)  call msg_debug2 &
         (D_ME_METHODS, "Recola writer: order = " // char (order))
    writer%order = order
  end subroutine recola_writer_set_order

  module subroutine recola_writer_set_coupling_powers &
       (writer, alpha_power, alphas_power)
    class(recola_writer_t), intent(inout) :: writer
    integer, intent(in) :: alpha_power
    integer, intent(in) :: alphas_power
    if (debug_on)  call msg_debug2 &
         (D_ME_METHODS, "Recola writer: alphas_power", alphas_power)
    if (debug_on)  call msg_debug2 &
         (D_ME_METHODS, "Recola writer: alpha_power", alpha_power)
    writer%alpha_power = alpha_power
    writer%alphas_power = alphas_power
  end subroutine recola_writer_set_coupling_powers

  function flv_file_name (id)
    type(string_t), intent(in) :: id
    type(string_t) :: flv_file_name
    flv_file_name = id // ".flv.dat"
  end function flv_file_name

  module subroutine recola_writer_write_makefile_code &
       (writer, unit, id, os_data, verbose, testflag)
    class(recola_writer_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    type(string_t) :: src_file
    type(string_t) :: flv_file
    call writer%base_write_makefile_code (unit, id, os_data, verbose, testflag)
    src_file = trim (char(id)) // ".f90"
    flv_file = flv_file_name (writer%id)
    write (unit, *)
    write (unit, "(5A)")  "# Flavor state listing for RECOLA process generation"
    write (unit, "(5A)")  char (flv_file), ": ", char (src_file)
    if (verbose) then
       write (unit, "(5A)", advance="no")  TAB
    else
       write (unit, "(5A)")  TAB, '@echo  "  MAKE      ', char (flv_file), '"'
       write (unit, "(5A)", advance="no")  TAB, "@"
    end if
    write (unit, "(5A)")  &
         "grep 'data table_flavor_states' $< ", &
         "| sed -e 's/.*\/\(.*\)\/.*/\1/' -e 's/,//g' > $@"
    write (unit, "(5A)")  "SOURCES += ", char (flv_file)
    write (unit, "(5A)")  "CLEAN_SOURCES += ", char (flv_file)
  end subroutine recola_writer_write_makefile_code

  module subroutine prc_recola_register_processes (writer, recola_ids)
    class(recola_writer_t), intent(in) :: writer
    integer, dimension (:), intent(inout) :: recola_ids
    integer :: recola_id
    integer :: i_flv
    integer :: n_tot
    integer :: unit, iostat
    integer, dimension(:), allocatable :: pdg
    type(string_t), dimension(:), allocatable :: particle_names
    type(string_t) :: process_string
    integer :: i_part
    !!! TODO (cw-2016-08-08): Include helicities
    call msg_message ("Recola: registering processes for '" // char (writer%id) // "'")
    i_flv = 0
    n_tot = writer%n_in + writer%n_out
    allocate (pdg (n_tot))
    allocate (particle_names (n_tot))
    call open_flv_list (writer%id, unit)
    call rclwrap_request_generate_processes ()
    SCAN_FLV_LIST: do
       read (unit, *, iostat = iostat)  pdg
       if (iostat < 0) then
          exit SCAN_FLV_LIST
       else if (iostat > 0) then
          call err_flv_list (writer%id)
       end if
       i_flv = i_flv + 1
       call rclwrap_get_new_recola_id (recola_id)
       recola_ids(i_flv) = recola_id
       particle_names(:) = get_recola_particle_string (pdg)
       process_string = var_str ("")
       do i_part = 1, n_tot
          process_string = process_string // &
               particle_names (i_part) // var_str (" ")
          if (i_part == writer%n_in) then
             process_string = process_string // var_str ("-> ")
          end if
       end do
       call msg_message ("Recola: " &
            // "process #" // char (str (recola_id)) &
            // ": " // char (process_string) &
            // "(" // char (writer%order) // ")")
       call rclwrap_add_process (recola_id, process_string, writer%order)
       call rclwrap_define_processes ()
    end do SCAN_FLV_LIST
    call close_flv_list (unit)
    if (debug_on) call msg_debug (D_ME_METHODS, "RECOLA: processes for '" &
         // char (writer%id) // "' registered")
  end subroutine prc_recola_register_processes

  subroutine open_flv_list (id, unit)
    type(string_t), intent(in) :: id
    integer, intent(out) :: unit
    type(string_t) :: flv_file
    integer :: iostat
    flv_file = flv_file_name (id)
    open (file = char (flv_file), newunit = unit, &
         status = "old", action = "read", &
         iostat = iostat)
    if (iostat /= 0) then
       call msg_fatal ("Recola: attempt to open flavor-list file '" &
            // char (flv_file) // "' failed")
    end if
  end subroutine open_flv_list

  subroutine err_flv_list (id)
    type(string_t), intent(in) :: id
    type(string_t) :: flv_file
    flv_file = flv_file_name (id)
    call msg_fatal ("Recola: error while reading from flavor-list file '" &
            // char (flv_file) // "'")
  end subroutine err_flv_list

  subroutine close_flv_list (unit)
    integer, intent(in) :: unit
    close (unit)
  end subroutine close_flv_list

  module function recola_driver_type_name () result (type)
    type(string_t) :: type
    type = "Recola"
  end function recola_driver_type_name

  module subroutine prc_recola_write_name (object, unit)
    class(prc_recola_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u,"(1x,A)") "Core: Recola"
  end subroutine prc_recola_write_name

  module function prc_recola_has_matrix_element (object) result (flag)
    logical :: flag
    class(prc_recola_t), intent(in) :: object
    flag = .true.
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
    p = 0
    if (associated (object%def)) then
       select type (def => object%def)
       type is (recola_def_t)
          p = def%alpha_power
       end select
    end if
  end function prc_recola_get_alpha_power

  module function prc_recola_get_alphas_power (object) result (p)
    class(prc_recola_t), intent(in) :: object
    integer :: p
    p = 0
    if (associated (object%def)) then
       select type (def => object%def)
       type is (recola_def_t)
          p = def%alphas_power
       end select
    end if
  end function prc_recola_get_alphas_power

  module subroutine prc_recola_compute_alpha_s (object, core_state, ren_scale)
    class(prc_recola_t), intent(in) :: object
    class(prc_external_state_t), intent(inout) :: core_state
    real(default), intent(in) :: ren_scale
    core_state%alpha_qcd = object%qcd%alpha%get (ren_scale)
  end subroutine prc_recola_compute_alpha_s

  module function prc_recola_includes_polarization (object) result (polarized)
    logical :: polarized
    class(prc_recola_t), intent(in) :: object
    polarized = .false.
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
    if (debug_on)  call msg_debug &
         (D_ME_METHODS, "prc_recola_prepare_external_code (no-op)")
  end subroutine prc_recola_prepare_external_code

  module subroutine prc_recola_set_parameters (object, qcd, model)
    class(prc_recola_t), intent(inout) :: object
    type(qcd_t), intent(in) :: qcd
    class(model_data_t), intent(in), target, optional :: model

    if (debug_on) call msg_debug (D_ME_METHODS, "RECOLA: set_parameters")
    object%qcd = qcd
    call rclwrap_set_dynamic_settings ()
    call rclwrap_set_pole_mass &
         (11, dble(model%get_real (var_str ('me'))), 0._double)
    call rclwrap_set_pole_mass &
         (13, dble(model%get_real (var_str ('mmu'))), 0._double)
    call rclwrap_set_pole_mass &
         (15, dble(model%get_real (var_str ('mtau'))), 0._double)

    call rclwrap_set_pole_mass (1, 0._double, 0._double)
    call rclwrap_set_pole_mass (2, 0._double, 0._double)

    call rclwrap_set_pole_mass (3, dble(model%get_real (var_str ('ms'))), 0._double)
    call rclwrap_set_pole_mass (4, dble(model%get_real (var_str ('mc'))), 0._double)
    call rclwrap_set_pole_mass (5, dble(model%get_real (var_str ('mb'))), 0._double)
    call rclwrap_set_pole_mass (6, dble(model%get_real (var_str ('mtop'))), &
         dble(model%get_real (var_str ('wtop'))))

    call rclwrap_set_pole_mass (23, dble(model%get_real (var_str ('mZ'))), &
         dble(model%get_real (var_str ('wZ'))))
    call rclwrap_set_pole_mass (24, dble(model%get_real (var_str ('mW'))), &
         dble(model%get_real (var_str ('wW'))))
    call rclwrap_set_pole_mass (25, dble(model%get_real (var_str ('mH'))), &
         dble(model%get_real (var_str ('wH'))))

    !!! TODO PB 03-03-2022: Automatize EW input schemes
    call rclwrap_use_gfermi_scheme (dble(model%get_real (var_str ('GF'))))
    !!! TODO PB 03-03-2022: Automatize mass threshold for light fermions
    call rclwrap_set_light_fermions (0._double)
    call rclwrap_set_delta_ir (0._double, dble(pi**2 / 6))
  end subroutine prc_recola_set_parameters

  module subroutine prc_recola_init (object, def, lib, id, i_component)
    class(prc_recola_t), intent(inout) :: object
    class(prc_core_def_t), intent(in), target :: def
    type(process_library_t), intent(in), target :: lib
    type(string_t), intent(in) :: id
    integer, intent(in) :: i_component
    integer :: n_flv
    if (debug_on) call msg_debug (D_ME_METHODS, "RECOLA: init process object")
    call object%base_init (def, lib, id, i_component)
    n_flv = size (object%data%flv_state, 2)
    allocate (object%recola_ids(n_flv))
    select type (writer => object%def%writer)
    type is (recola_writer_t)
       call writer%register_processes (object%recola_ids)
    end select
    call rclwrap_generate_processes ()
    call object%replace_helicity_and_color_arrays ()
  end subroutine prc_recola_init

  module subroutine prc_recola_replace_helicity_and_color_arrays (object)
    class(prc_recola_t), intent(inout) :: object
    integer, dimension(:,:), allocatable :: col_recola
    integer :: i
    if (debug_on)  call msg_debug &
         (D_ME_METHODS, "RECOLA: replace_helicity_and_color_arrays")
    deallocate (object%data%hel_state)
    call rclwrap_get_helicity_configurations &
         (object%recola_ids(1), object%data%hel_state)
    call rclwrap_get_color_configurations (object%recola_ids(1), col_recola)
    allocate (object%color_state (object%data%n_in + object%data%n_out, &
           size (col_recola, dim = 2)))
    do i = 1, size (col_recola, dim = 2)
       object%color_state (:, i) = col_recola (:, i)
    end do
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
    real(double), dimension(0:3, object%data%n_in + object%data%n_out) :: &
         p_recola
    integer :: i
    logical :: new_event
    complex(double) :: amp_dble

    if (debug_on) call msg_debug2 (D_ME_METHODS, "prc_recola_compute_amplitude")
    if (present (core_state)) then
       if (allocated (core_state)) then
          select type (core_state)
          type is (recola_state_t)
             new_event = core_state%new_kinematics
             core_state%new_kinematics = .false.
          end select
       end if
    end if
    if (new_event) then
       do i = 1, object%data%n_in + object%data%n_out
          p_recola(:, i) = dble(p(i)%p)
       end do
       call rclwrap_compute_process (object%recola_ids(f), p_recola, 'LO')
    end if
    call rclwrap_get_amplitude (object%recola_ids(f), 0, 'LO', &
         object%color_state (:, c), object%data%hel_state (h, :), amp_dble)
    amp = amp_dble
  end function prc_recola_compute_amplitude

  module subroutine prc_recola_compute_sqme (object, i_flv, i_hel, p, &
       ren_scale, sqme, bad_point)
    class(prc_recola_t), intent(in) :: object
    integer, intent(in) :: i_flv, i_hel
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(in) :: ren_scale
    real(default), intent(out) :: sqme
    logical, intent(out) :: bad_point
    real(double) :: sqme_dble
    real(double), dimension(0:3, object%data%n_in + object%data%n_out) :: &
         p_recola
    real(default) :: alpha_s
    integer :: i
    integer :: alphas_power
    ! TODO sbrass: Helicity for RECOLA
    if (debug_on) call msg_debug2 (D_ME_METHODS, "prc_recola_compute_sqme")
    do i = 1, object%data%n_in + object%data%n_out
       p_recola(:, i) = dble(p(i)%p)
    end do
    alpha_s = object%qcd%alpha%get (ren_scale)
    if (debug_on) call msg_debug2 (D_ME_METHODS, "alpha_s", alpha_s)
    if (debug_on) call msg_debug2 (D_ME_METHODS, "ren_scale", ren_scale)
    call rclwrap_set_alpha_s (dble (alpha_s), dble (ren_scale), object%qcd%n_f)
    call rclwrap_set_mu_ir (dble (ren_scale))
    call rclwrap_compute_process (object%recola_ids(i_flv), p_recola, 'LO')
    call rclwrap_get_squared_amplitude &
         (object%recola_ids(i_flv), object%get_alphas_power (), 'LO', sqme_dble)
    sqme = real(sqme_dble, kind=default)
    bad_point = .false.
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
    real(double), dimension(0:3, object%data%n_in + object%data%n_out) :: &
         p_recola
    real(double) :: sqme_dble
    real(default) :: alpha_s
    integer :: i, as_coupling_power
    logical :: ew_corr
    ! TODO sbrass Helicity for RECOLA
    if (debug_on) call msg_debug2 (D_ME_METHODS, "prc_recola_compute_sqme_virt")
    sqme = zero
    do i = 1, object%data%n_in + object%data%n_out
       p_recola(:, i) = dble(p(i)%p)
    end do
    call rclwrap_set_mu_ir (dble (ren_scale))
    alpha_s = object%qcd%alpha%get (ren_scale)
    call rclwrap_set_alpha_s (dble (alpha_s), dble (ren_scale), object%qcd%n_f)
    call rclwrap_compute_process (object%recola_ids(i_flv), p_recola, 'NLO')
    if (associated (object%def)) then
       select type (def => object%def)
       type is (recola_def_t)
          ew_corr = def%corr == RECOLA_EW
       end select
    end if
    if (ew_corr) then
       as_coupling_power = object%get_alphas_power ()
    else
       as_coupling_power = object%get_alphas_power () + 1
    end if
    call rclwrap_get_squared_amplitude (object%recola_ids(i_flv), &
         as_coupling_power, 'NLO', sqme_dble)
    sqme(3) = sqme_dble 
    call rclwrap_get_squared_amplitude &
         (object%recola_ids(i_flv), object%get_alphas_power (), 'LO', sqme_dble)
    sqme(4) = sqme_dble
    bad_point = .false.
  end subroutine prc_recola_compute_sqme_virt

  module function prc_recola_get_alpha_qed &
       (object, core_state) result (alpha_qed)
    class(prc_recola_t), intent(in) :: object
    class(prc_core_state_t), intent(in), allocatable :: core_state
    real(double) :: value
    real(default) :: alpha_qed
       alpha_qed = rclwrap_get_alpha ()
  end function prc_recola_get_alpha_qed

  module subroutine prc_recola_compute_sqme_color_c_raw (object, &
       i_flv, i_hel, p, ren_scale, sqme_color_c, bad_point)
    class(prc_recola_t), intent(in) :: object
    integer, intent(in) :: i_hel, i_flv
    type(vector4_t), dimension(:), intent(in) :: p
    real(double), dimension(0:3, object%data%n_in + object%data%n_out) :: &
         p_recola
    real(default), intent(in) :: ren_scale
    real(default), dimension(:), intent(out) :: sqme_color_c
    logical, intent(out) :: bad_point
    integer :: i1, i2, i, n_tot
    real(double) :: sqme_dble
    do i = 1, object%data%n_in + object%data%n_out
       p_recola(:, i) = dble(p(i)%p)
    end do
    n_tot = object%data%n_in + object%data%n_out
    i = 0
    do i1 = 1, n_tot
       do i2 = 1, i1-1
          i = i + 1
          call rclwrap_compute_color_correlation &
               (object%recola_ids(i_flv), p_recola, i1, i2, sqme_dble)
          sqme_color_c(i) = real (sqme_dble, kind=default)
          select case (abs (object%data%flv_state (i1, i_flv)))
          case (1:6)
             sqme_color_c(i) = CF * sqme_color_c(i)
          case (9,21)
             sqme_color_c(i) = CA * sqme_color_c(i)
          end select
       end do
    end do
  end subroutine prc_recola_compute_sqme_color_c_raw


end submodule prc_recola_s

