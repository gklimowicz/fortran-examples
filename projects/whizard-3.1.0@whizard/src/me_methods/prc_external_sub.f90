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

submodule (prc_external) prc_external_s

  use debug_master, only: debug_on
  use io_units
  use pdg_arrays, only: is_gluon, is_quark
  use system_defs, only: TAB
  use physics_defs, only: CF
  use diagnostics
  use prc_omega, only: omega_state_t

  implicit none

  integer, parameter :: LEPTONS = 1
  integer, parameter :: HADRONS = 2

contains

  module subroutine sf_handler_init (sf_handler, sf_chain)
    class(sf_handler_t), intent(out) :: sf_handler
    type(sf_chain_instance_t), intent(in) :: sf_chain
    integer :: i
    sf_handler%n_sf = size (sf_chain%sf)
    if (sf_handler%n_sf == 0) then
       sf_handler%initial_state_type = LEPTONS
    else
       do i = 1, sf_handler%n_sf
          select type (int => sf_chain%sf(i)%int)
          type is (pdf_builtin_t)
             sf_handler%initial_state_type = HADRONS
          type is (lhapdf_t)
             sf_handler%initial_state_type = HADRONS
          class default
             sf_handler%initial_state_type = LEPTONS
          end select
       end do
     end if
  end subroutine sf_handler_init

  module subroutine sf_handler_init_dummy (sf_handler)
    class(sf_handler_t), intent(out) :: sf_handler
    sf_handler%n_sf = 0
    sf_handler%initial_state_type = LEPTONS
  end subroutine sf_handler_init_dummy

  module subroutine sf_handler_apply_structure_functions &
       (sf_handler, sf_chain, flavors)
     class(sf_handler_t), intent(inout) :: sf_handler
     type(sf_chain_instance_t), intent(in) :: sf_chain
     integer, intent(in), dimension(2) :: flavors
     integer :: i
     real(default), dimension(:), allocatable :: f
     if (sf_handler%n_sf < 0) call msg_fatal ("sf_handler not initialized")
     sf_handler%val = one
     do i = 1, sf_handler%n_sf
        select case (sf_handler%initial_state_type)
        case (HADRONS)
           sf_handler%val = sf_handler%val * &
                sf_handler%get_pdf (sf_chain, i, flavors(i))
        case (LEPTONS)
           call sf_chain%get_matrix_elements (i, f)
           sf_handler%val = sf_handler%val * f(1)
        case default
           call msg_fatal ("sf_handler not initialized")
        end select
     end do
  end subroutine sf_handler_apply_structure_functions

  module function sf_handler_get_pdf &
       (sf_handler, sf_chain, i, flavor) result (f)
     real(default) :: f
     class(sf_handler_t), intent(in) :: sf_handler
     type(sf_chain_instance_t), intent(in) :: sf_chain
     integer, intent(in) :: i, flavor
     integer :: k
     real(default), dimension(:), allocatable :: ff
     integer, parameter :: n_flv_light = 6

     call sf_chain%get_matrix_elements (i, ff)

     if (is_gluon (flavor)) then
        k = n_flv_light + 1
     else if (is_quark (abs(flavor))) then
        k = n_flv_light + 1 + flavor
     else
        call msg_fatal ("Not a colored particle")
     end if

     f = ff(k)
  end function sf_handler_get_pdf

  module subroutine prc_external_state_reset_new_kinematics (object)
    class(prc_external_state_t), intent(inout) :: object
    object%new_kinematics = .true.
  end subroutine prc_external_state_reset_new_kinematics

  module function prc_external_needs_external_code () result (flag)
    logical :: flag
    flag = .true.
  end function prc_external_needs_external_code

  pure module function prc_external_get_n_flvs (object, i_flv) result (n)
    integer :: n
    class(prc_external_t), intent(in) :: object
    integer, intent(in) :: i_flv
    n = size (object%data%flv_state (:,i_flv))
  end function prc_external_get_n_flvs

  module function prc_external_get_flv_state (object, i_flv) result (flv)
    integer, dimension(:), allocatable :: flv
    class(prc_external_t), intent(in) :: object
    integer, intent(in) :: i_flv
    allocate (flv (size (object%data%flv_state (:,i_flv))))
    flv = object%data%flv_state (:,i_flv)
  end function prc_external_get_flv_state

  module subroutine prc_external_compute_sqme (object, i_flv, i_hel, p, &
         ren_scale, sqme, bad_point)
     class(prc_external_t), intent(in) :: object
     integer, intent(in) :: i_flv, i_hel
     type(vector4_t), dimension(:), intent(in) :: p
     real(default), intent(in) :: ren_scale
     real(default), intent(out) :: sqme
     logical, intent(out) :: bad_point
     sqme = one
     bad_point = .false.
  end subroutine prc_external_compute_sqme

  module subroutine prc_external_compute_sqme_virt (object, i_flv, i_hel, &
     p, ren_scale, es_scale, loop_method, sqme, bad_point)
    class(prc_external_t), intent(in) :: object
    integer, intent(in) :: i_flv, i_hel
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(in) :: ren_scale, es_scale
    integer, intent(in) :: loop_method
    logical, intent(out) :: bad_point
    real(default), dimension(4), intent(out) :: sqme
    if (debug_on) call msg_debug2 &
         (D_ME_METHODS, "prc_external_compute_sqme_virt")
    sqme(1) = 0.001_default
    sqme(2) = 0.001_default
    sqme(3) = 0.001_default
    sqme(4) = 0.0015_default
    bad_point = .false.
  end subroutine prc_external_compute_sqme_virt

  module subroutine prc_external_compute_sqme_color_c (object, i_flv, &
       i_hel, p, ren_scale, born_color_c, bad_point, born_out)
    class(prc_external_t), intent(inout) :: object
    integer, intent(in) :: i_flv, i_hel
    type(vector4_t), intent(in), dimension(:) :: p
    real(default), intent(in) :: ren_scale
    real(default), intent(inout), dimension(:,:) :: born_color_c
    logical, intent(out) :: bad_point
    real(default), intent(out), optional :: born_out
    if (debug_on) call msg_debug2 (D_ME_METHODS, "prc_external_compute_sqme_color_c")
    if (size (p) == 4) then
       if (present (born_out)) then
          born_out = 0.0015_default
          born_color_c = zero
          born_color_c(3,3) = - CF * born_out
          born_color_c(4,4) = - CF * born_out
          born_color_c(3,4) = CF * born_out
          born_color_c(4,3) = born_color_c(3,4)
          bad_point = .false.
       end if
    else
       if (present (born_out)) born_out = zero
       born_color_c = zero
    end if
  end subroutine prc_external_compute_sqme_color_c

  module subroutine prc_external_compute_alpha_s &
       (object, core_state, ren_scale)
    class(prc_external_t), intent(in) :: object
    class(prc_external_state_t), intent(inout) :: core_state
    real(default), intent(in) :: ren_scale
    core_state%alpha_qcd = object%qcd%alpha%get (ren_scale)
  end subroutine prc_external_compute_alpha_s

  module function prc_external_get_alpha_s &
       (object, core_state) result (alpha_qcd)
    class(prc_external_t), intent(in) :: object
    class(prc_core_state_t), intent(in), allocatable :: core_state
    real(default) :: alpha_qcd
    if (allocated (core_state)) then
      select type (core_state)
      class is (prc_external_state_t)
         alpha_qcd = core_state%alpha_qcd
      type is (omega_state_t)
         alpha_qcd = core_state%alpha_qcd
      class default
         alpha_qcd = zero
      end select
    else
       alpha_qcd = zero
    end if
  end function prc_external_get_alpha_s

  module function prc_external_get_alpha_qed &
       (object, core_state) result (alpha_qed)
    class(prc_external_t), intent(in) :: object
    class(prc_core_state_t), intent(in), allocatable :: core_state
    real(default) :: alpha_qed
    if (allocated (core_state)) then
       select type (core_state)
       class is (prc_external_state_t)
          alpha_qed = core_state%alpha_qed
       type is (omega_state_t)
          alpha_qed = core_state%alpha_qed
       class default
          alpha_qed = zero
       end select
    else
       alpha_qed = zero
    end if
  end function prc_external_get_alpha_qed

  module function prc_external_is_allowed &
       (object, i_term, f, h, c) result (flag)
    class(prc_external_t), intent(in) :: object
    integer, intent(in) :: i_term, f, h, c
    logical :: flag
    logical(c_bool) :: cflag
    select type (driver => object%driver)
    class is (prc_external_driver_t)
       call driver%is_allowed (f, h, c, cflag)
       flag = cflag
    class default
       call msg_fatal &
          ("Driver does not fit to prc_external_t")
    end select
  end function prc_external_is_allowed

  module function prc_external_get_nflv (object) result (n_flv)
    class(prc_external_t), intent(in) :: object
    integer :: n_flv
    n_flv = object%n_flv
  end function prc_external_get_nflv

  module subroutine prc_external_compute_hard_kinematics &
       (object, p_seed, i_term, int_hard, core_state)
    class(prc_external_t), intent(in) :: object
    type(vector4_t), dimension(:), intent(in) :: p_seed
    integer, intent(in) :: i_term
    type(interaction_t), intent(inout) :: int_hard
    class(prc_core_state_t), intent(inout), allocatable :: core_state
    call int_hard%set_momenta (p_seed)
    if (allocated (core_state)) then
      select type (core_state)
      class is (prc_external_state_t); core_state%new_kinematics = .true.
      end select
    end if
  end subroutine prc_external_compute_hard_kinematics

  module subroutine prc_external_compute_eff_kinematics &
       (object, i_term, int_hard, int_eff, core_state)
    class(prc_external_t), intent(in) :: object
    integer, intent(in) :: i_term
    type(interaction_t), intent(in) :: int_hard
    type(interaction_t), intent(inout) :: int_eff
    class(prc_core_state_t), intent(inout), allocatable :: core_state
  end subroutine prc_external_compute_eff_kinematics

  module subroutine prc_external_update_alpha_s (object, core_state, scale)
    class(prc_external_t), intent(in) :: object
    class(prc_core_state_t), intent(inout), allocatable :: core_state
    real(default), intent(in) :: scale
    real(default) :: alpha_qcd
    if (allocated (object%qcd%alpha)) then
       alpha_qcd = object%qcd%alpha%get (scale)
       select type (driver => object%driver)
       class is (prc_external_driver_t)
          call driver%update_alpha_s (alpha_qcd)
       end select
       select type (core_state)
       class is (prc_external_state_t)
          core_state%alpha_qcd = alpha_qcd
       type is (omega_state_t)
          core_state%alpha_qcd = alpha_qcd
       end select
    end if
  end subroutine prc_external_update_alpha_s

  module subroutine prc_external_init_sf_handler (core, sf_chain)
     class(prc_external_t), intent(inout) :: core
     type(sf_chain_instance_t), intent(in) :: sf_chain
     if (allocated (sf_chain%sf)) then
        call core%sf_handler%init (sf_chain)
     else
        call core%sf_handler%init_dummy ()
     end if
  end subroutine prc_external_init_sf_handler

  module subroutine prc_external_init_sf_handler_dummy (core)
     class(prc_external_t), intent(inout) :: core
     call core%sf_handler%init_dummy ()
  end subroutine prc_external_init_sf_handler_dummy

  module subroutine prc_external_apply_structure_functions &
       (core, sf_chain, flavors)
    class(prc_external_t), intent(inout) :: core
    type(sf_chain_instance_t), intent(in) :: sf_chain
    integer, dimension(2), intent(in) :: flavors
    call core%sf_handler%apply_structure_functions (sf_chain, flavors)
  end subroutine prc_external_apply_structure_functions

  module function prc_external_get_sf_value (core) result (val)
    real(default) :: val
    class(prc_external_t), intent(in) :: core
    val = core%sf_handler%val
  end function prc_external_get_sf_value

  module subroutine prc_external_def_set_active_writer (def, active)
    class(prc_external_def_t), intent(inout) :: def
    logical, intent(in) :: active
    select type (writer => def%writer)
    class is (prc_external_writer_t)
       writer%active = active
    end select
  end subroutine prc_external_def_set_active_writer

  module subroutine prc_external_def_get_features (features)
    type(string_t), dimension(:), allocatable, intent(out) :: features
    allocate (features (6))
    features = [ &
         var_str ("init"), &
         var_str ("update_alpha_s"), &
         var_str ("reset_helicity_selection"), &
         var_str ("is_allowed"), &
         var_str ("new_event"), &
         var_str ("get_amplitude")]
  end subroutine prc_external_def_get_features

  module subroutine prc_external_def_connect (def, lib_driver, i, proc_driver)
    class(prc_external_def_t), intent(in) :: def
    class(prclib_driver_t), intent(in) :: lib_driver
    integer, intent(in) :: i
    class(prc_core_driver_t), intent(inout) :: proc_driver
    integer :: pid, fid
    type(c_funptr) :: fptr
    select type (proc_driver)
    class is (prc_external_driver_t)
       pid = i
       fid = 2
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%update_alpha_s)
       fid = 4
       call lib_driver%get_fptr (pid, fid, fptr)
       call c_f_procpointer (fptr, proc_driver%is_allowed)
    end select
  end subroutine prc_external_def_connect

  module function prc_external_def_needs_code () result (flag)
    logical :: flag
    flag = .true.
  end function prc_external_def_needs_code

  pure module subroutine prc_external_writer_init &
       (writer, model_name, prt_in, prt_out, restrictions)
    class(prc_external_writer_t), intent(inout) :: writer
    type(string_t), intent(in) :: model_name
    type(string_t), dimension(:), intent(in) :: prt_in, prt_out
    type(string_t), intent(in), optional :: restrictions
    integer :: i
    writer%model_name = model_name
    if (present (restrictions)) then
       writer%restrictions = restrictions
    else
       writer%restrictions = ""
    end if
    writer%n_in = size (prt_in)
    writer%n_out = size (prt_out)
    select case (size (prt_in))
       case(1); writer%process_mode = " -decay"
       case(2); writer%process_mode = " -scatter"
    end select
    associate (s => writer%process_string)
      s = " '"
      do i = 1, size (prt_in)
         if (i > 1) s = s // " "
         s = s // prt_in(i)
      end do
      s = s // " ->"
      do i = 1, size (prt_out)
         s = s // " " // prt_out(i)
      end do
      s = s // "'"
    end associate
  end subroutine prc_external_writer_init

  module function prc_external_writer_get_module_name (id) result (name)
    type(string_t) :: name
    type(string_t), intent(in) :: id
    name = "opr_" // id
  end function prc_external_writer_get_module_name

  module subroutine prc_external_writer_write_wrapper &
       (writer, unit, id, feature)
    class(prc_external_writer_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id, feature
    type(string_t) :: name
    name = writer%get_c_procname (id, feature)
    write (unit, *)
    select case (char (feature))
    case ("init")
       write (unit, "(9A)")  "subroutine ", char (name), &
            " (par, scheme) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use opr_", char (id)
       write (unit, "(2x,9A)")  "real(c_default_float), dimension(*), &
            &intent(in) :: par"
       write (unit, "(2x,9A)")  "integer(c_int), intent(in) :: scheme"
       if (c_default_float == default .and. c_int == kind (1)) then
          write (unit, "(2x,9A)")  "call ", char (feature), " (par, scheme)"
       end if
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("update_alpha_s")
       write (unit, "(9A)")  "subroutine ", char (name), " (alpha_s) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use opr_", char (id)
       if (c_default_float == default) then
          write (unit, "(2x,9A)")  "real(c_default_float), intent(in) &
               &:: alpha_s"
          write (unit, "(2x,9A)")  "call ", char (feature), " (alpha_s)"
       end if
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("reset_helicity_selection")
       write (unit, "(9A)")  "subroutine ", char (name), &
            " (threshold, cutoff) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use opr_", char (id)
       if (c_default_float == default) then
          write (unit, "(2x,9A)")  "real(c_default_float), intent(in) &
               &:: threshold"
          write (unit, "(2x,9A)")  "integer(c_int), intent(in) :: cutoff"
          write (unit, "(2x,9A)")  "call ", char (feature), &
               " (threshold, int (cutoff))"
       end if
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("is_allowed")
       write (unit, "(9A)")  "subroutine ", char (name), &
            " (flv, hel, col, flag) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use opr_", char (id)
       write (unit, "(2x,9A)")  "integer(c_int), intent(in) :: flv, hel, col"
       write (unit, "(2x,9A)")  "logical(c_bool), intent(out) :: flag"
       write (unit, "(2x,9A)")  "flag = ", char (feature), &
            " (int (flv), int (hel), int (col))"
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("new_event")
       write (unit, "(9A)")  "subroutine ", char (name), " (p) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use opr_", char (id)
       if (c_default_float == default) then
          write (unit, "(2x,9A)")  "real(c_default_float), dimension(0:3,*), &
               &intent(in) :: p"
          write (unit, "(2x,9A)")  "call ", char (feature), " (p)"
       end if
       write (unit, "(9A)")  "end subroutine ", char (name)
    case ("get_amplitude")
       write (unit, "(9A)")  "subroutine ", char (name), &
            " (flv, hel, col, amp) bind(C)"
       write (unit, "(2x,9A)")  "use iso_c_binding"
       write (unit, "(2x,9A)")  "use kinds"
       write (unit, "(2x,9A)")  "use opr_", char (id)
       write (unit, "(2x,9A)")  "integer(c_int), intent(in) :: flv, hel, col"
       write (unit, "(2x,9A)")  "complex(c_default_complex), intent(out) &
            &:: amp"
       write (unit, "(2x,9A)")  "amp = ", char (feature), &
            " (int (flv), int (hel), int (col))"
       write (unit, "(9A)")  "end subroutine ", char (name)
    end select

  end subroutine prc_external_writer_write_wrapper

  module subroutine prc_external_writer_write_interface &
       (writer, unit, id, feature)
    class(prc_external_writer_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(string_t), intent(in) :: feature
    type(string_t) :: name
    name = writer%get_c_procname (id, feature)
    write (unit, "(2x,9A)")  "interface"
    select case (char (feature))
    case ("init")
       write (unit, "(5x,9A)")  "subroutine ", char (name), &
            " (par, scheme) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "real(c_default_float), dimension(*), &
            &intent(in) :: par"
       write (unit, "(7x,9A)")  "integer(c_int), intent(in) :: scheme"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    case ("update_alpha_s")
       write (unit, "(5x,9A)")  "subroutine ", char (name), " (alpha_s) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "real(c_default_float), intent(in) :: alpha_s"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    case ("reset_helicity_selection")
       write (unit, "(5x,9A)")  "subroutine ", char (name), " &
            &(threshold, cutoff) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "real(c_default_float), intent(in) :: threshold"
       write (unit, "(7x,9A)")  "integer(c_int), intent(in) :: cutoff"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    case ("is_allowed")
       write (unit, "(5x,9A)")  "subroutine ", char (name), " &
            &(flv, hel, col, flag) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "integer(c_int), intent(in) :: flv, hel, col"
       write (unit, "(7x,9A)")  "logical(c_bool), intent(out) :: flag"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    case ("new_event")
       write (unit, "(5x,9A)")  "subroutine ", char (name), " (p) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "real(c_default_float), dimension(0:3,*), &
            &intent(in) :: p"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    case ("get_amplitude")
       write (unit, "(5x,9A)")  "subroutine ", char (name), " &
            &(flv, hel, col, amp) bind(C)"
       write (unit, "(7x,9A)")  "import"
       write (unit, "(7x,9A)")  "integer(c_int), intent(in) :: flv, hel, col"
       write (unit, "(7x,9A)")  "complex(c_default_complex), intent(out) &
            &:: amp"
       write (unit, "(5x,9A)")  "end subroutine ", char (name)
    end select
    write (unit, "(2x,9A)")  "end interface"
  end subroutine prc_external_writer_write_interface

  module subroutine prc_external_writer_write_source_code (writer, id)
    class(prc_external_writer_t), intent(in) :: writer
    type(string_t), intent(in) :: id
    if (debug_on) call msg_debug (D_ME_METHODS, &
         "prc_external_writer_write_source_code (no-op)")
    !!! This is a dummy
  end subroutine prc_external_writer_write_source_code

  module subroutine prc_external_writer_before_compile (writer, id)
    class(prc_external_writer_t), intent(in) :: writer
    type(string_t), intent(in) :: id
    if (debug_on) call msg_debug (D_ME_METHODS, &
         "prc_external_writer_before_compile (no-op)")
    !!! This is a dummy
  end subroutine prc_external_writer_before_compile

  module subroutine prc_external_writer_after_compile (writer, id)
    class(prc_external_writer_t), intent(in) :: writer
    type(string_t), intent(in) :: id
    if (debug_on) call msg_debug (D_ME_METHODS, &
         "prc_external_writer_after_compile (no-op)")
    !!! This is a dummy
  end subroutine prc_external_writer_after_compile

  module subroutine prc_external_writer_write_makefile_code &
       (writer, unit, id, os_data, verbose, testflag)
    class(prc_external_writer_t), intent(in) :: writer
    integer, intent(in) :: unit
    type(string_t), intent(in) :: id
    type(os_data_t), intent(in) :: os_data
    logical, intent(in) :: verbose
    logical, intent(in), optional :: testflag
    type(string_t) :: omega_binary, omega_path
    type(string_t) :: restrictions_string, amp_triv_string
    omega_binary = "omega_" // writer%model_name // ".opt"
    omega_path = os_data%whizard_omega_binpath // "/" // omega_binary
    if (.not. verbose)  omega_path = "@" // omega_path
    if (writer%restrictions /= "") then
       restrictions_string = " -cascade '" // writer%restrictions // "'"
    else
       restrictions_string = ""
    end if
    amp_triv_string = ""
    if (writer%amp_triv)  amp_triv_string = " -target:amp_triv"
    write (unit, "(5A)")  "OBJECTS += ", char (id), ".lo"
    write (unit, "(5A)")  char (id), ".f90:"
    if (.not. verbose) then
       write (unit, "(5A)")  TAB // '@echo  "  OMEGA     ', trim (char (id)), '.f90"'
    end if
    write (unit, "(99A)")  TAB, char (omega_path), &
         " -o ", char (id), ".f90", &
         " -target:whizard", char (amp_triv_string), &
         " -target:parameter_module parameters_", char (writer%model_name), &
         " -target:module opr_", char (id), &
         " -target:md5sum '", writer%md5sum, "'", &
         char (writer%process_mode), char (writer%process_string), &
         char (restrictions_string)
    write (unit, "(5A)")  "clean-", char (id), ":"
    write (unit, "(5A)")  TAB, "rm -f ", char (id), ".f90"
    write (unit, "(5A)")  TAB, "rm -f opr_", char (id), ".mod"
    write (unit, "(5A)")  TAB, "rm -f ", char (id), ".lo"
    write (unit, "(5A)")  "CLEAN_SOURCES += ", char (id), ".f90"
    write (unit, "(5A)")  "CLEAN_OBJECTS += opr_", char (id), ".mod"
    write (unit, "(5A)")  "CLEAN_OBJECTS += ", char (id), ".lo"
    write (unit, "(5A)")  char (id), ".lo: ", char (id), ".f90"
    if (.not. verbose) then
       write (unit, "(5A)")  TAB // '@echo  "  FC       " $@'
    end if
    write (unit, "(5A)")  TAB, "$(LTFCOMPILE) $<"

  end subroutine prc_external_writer_write_makefile_code

  module function prc_external_writer_writer_get_procname &
       (feature) result (name)
    type(string_t) :: name
    type(string_t), intent(in) :: feature
    select case (char (feature))
    case ("n_in");   name = "number_particles_in"
    case ("n_out");  name = "number_particles_out"
    case ("n_flv");  name = "number_flavor_states"
    case ("n_hel");  name = "number_spin_states"
    case ("n_col");  name = "number_color_flows"
    case ("n_cin");  name = "number_color_indices"
    case ("n_cf");   name = "number_color_factors"
    case ("flv_state");  name = "flavor_states"
    case ("hel_state");  name = "spin_states"
    case ("col_state");  name = "color_flows"
    case default
       name = feature
    end select
  end function prc_external_writer_writer_get_procname

  module function prc_external_test_writer_type_name () result (string)
    type(string_t) :: string
    string = "External matrix element dummy"
  end function prc_external_test_writer_type_name

  module subroutine prc_external_test_state_write (object, unit)
    class(prc_external_test_state_t), intent(in) :: object
    integer, intent(in), optional :: unit
  end subroutine prc_external_test_state_write

  module function prc_external_test_driver_type_name () result (type)
    type(string_t) :: type
    type = "External matrix element dummy"
  end function prc_external_test_driver_type_name

  module function prc_external_test_def_type_string () result (string)
    type(string_t) :: string
    string = "external test dummy"
  end function prc_external_test_def_type_string

  module subroutine prc_external_test_def_write (object, unit)
    class(prc_external_test_def_t), intent(in) :: object
    integer, intent(in) :: unit
  end subroutine prc_external_test_def_write

  module subroutine prc_external_test_def_read (object, unit)
    class(prc_external_test_def_t), intent(out) :: object
    integer, intent(in) :: unit
  end subroutine prc_external_test_def_read

  module subroutine prc_external_test_write (object, unit)
    class(prc_external_test_t), intent(in) :: object
    integer, intent(in), optional :: unit
    call msg_message ("Test external matrix elements")
  end subroutine prc_external_test_write

  module subroutine prc_external_test_write_name (object, unit)
    class(prc_external_test_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u,"(1x,A)") "Core: external test"
  end subroutine prc_external_test_write_name

  module function prc_external_test_compute_amplitude &
       (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
       core_state)  result (amp)
    class(prc_external_test_t), intent(in) :: object
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
    amp = 0.0
  end function prc_external_test_compute_amplitude

  module function prc_external_test_includes_polarization &
       (object) result (polarized)
    logical :: polarized
    class(prc_external_test_t), intent(in) :: object
    polarized = .false.
  end function prc_external_test_includes_polarization

  module subroutine prc_external_test_prepare_external_code &
       (core, flv_states, var_list, os_data, libname, model, i_core, is_nlo)
    class(prc_external_test_t), intent(inout) :: core
    integer, intent(in), dimension(:,:), allocatable :: flv_states
    type(var_list_t), intent(in) :: var_list
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in) :: libname
    type(model_data_t), intent(in), target :: model
    integer, intent(in) :: i_core
    logical, intent(in) :: is_nlo
  end subroutine prc_external_test_prepare_external_code


end submodule prc_external_s

