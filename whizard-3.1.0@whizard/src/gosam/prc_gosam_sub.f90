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

submodule (prc_gosam) prc_gosam_s

  use io_units
  use constants
  use numeric_utils
  use system_defs, only: TAB
  use system_dependencies
  use file_utils
  use string_utils
  use diagnostics
  use sm_qcd
  use flavors
  use pdg_arrays
  use process_constants
  use prclib_interfaces
  use process_libraries

  implicit none

contains

  module subroutine gosam_writer_write_config (gosam_writer)
    class(gosam_writer_t), intent(in) :: gosam_writer
    integer :: unit
    unit = free_unit ()
    open (unit, file = "golem.in", status = "replace", action = "write")
    call gosam_writer%generate_configuration_file (unit)
    close(unit)
  end subroutine gosam_writer_write_config

  module function gosam_def_type_string () result (string)
    type(string_t) :: string
    string = "gosam"
  end function gosam_def_type_string

  module subroutine gosam_def_write (object, unit)
    class(gosam_def_t), intent(in) :: object
    integer, intent(in) :: unit
    select type (writer => object%writer)
    type is (gosam_writer_t)
      call writer%write (unit)
    end select
  end subroutine gosam_def_write

  module subroutine gosam_def_read (object, unit)
    class(gosam_def_t), intent(out) :: object
    integer, intent(in) :: unit
  end subroutine gosam_def_read

  module function gosam_writer_type_name () result (string)
    type(string_t) :: string
    string = "gosam"
  end function gosam_writer_type_name

  pure module subroutine gosam_writer_init &
       (writer, model_name, prt_in, prt_out, restrictions)
    class(gosam_writer_t), intent(inout) :: writer
    type(string_t), intent(in) :: model_name
    type(string_t), dimension(:), intent(in) :: prt_in, prt_out
    type(string_t), intent(in), optional :: restrictions
    writer%gosam_dir = GOSAM_DIR
    writer%golem_dir = GOLEM_DIR
    writer%samurai_dir = SAMURAI_DIR
    writer%ninja_dir = NINJA_DIR
    writer%form_dir = FORM_DIR
    writer%qgraf_dir = QGRAF_DIR
    call writer%base_init (model_name, prt_in, prt_out)
  end subroutine gosam_writer_init

  module function gosam_driver_type_name () result (string)
    type(string_t) :: string
    string = "gosam"
  end function gosam_driver_type_name

  module subroutine gosam_driver_init_gosam (object, os_data, olp_file, &
       olc_file, olp_dir, olp_lib)
    class(gosam_driver_t), intent(inout) :: object
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in) :: olp_file, olc_file, olp_dir, olp_lib
    object%gosam_dir = GOSAM_DIR
    object%olp_file = olp_file
    object%contract_file = olc_file
    object%olp_dir = olp_dir
    object%olp_lib = olp_lib
  end subroutine gosam_driver_init_gosam

  module subroutine gosam_driver_init_dlaccess_to_library &
     (object, os_data, dlaccess, success)
    class(gosam_driver_t), intent(in) :: object
    type(os_data_t), intent(in) :: os_data
    type(dlaccess_t), intent(out) :: dlaccess
    logical, intent(out) :: success
    type(string_t) :: libname, msg_buffer
    libname = object%olp_dir // '/.libs/libgolem_olp.' // &
      os_data%shrlib_ext
    msg_buffer = "One-Loop-Provider: Using Gosam"
    call msg_message (char(msg_buffer))
    msg_buffer = "Loading library: " // libname
    call msg_message (char(msg_buffer))
    call dlaccess_init (dlaccess, var_str ("."), libname, os_data)
    success = .not. dlaccess_has_error (dlaccess)
  end subroutine gosam_driver_init_dlaccess_to_library

  module subroutine gosam_writer_generate_configuration_file &
       (object, unit)
    class(gosam_writer_t), intent(in) :: object
    integer, intent(in) :: unit
    type(string_t) :: fc_bin
    type(string_t) :: form_bin, qgraf_bin, haggies_bin
    type(string_t) :: fcflags_golem, ldflags_golem
    type(string_t) :: fcflags_samurai, ldflags_samurai
    type(string_t) :: fcflags_ninja, ldflags_ninja
    type(string_t) :: ldflags_avh_olo, ldflags_qcdloop
    fc_bin = DEFAULT_FC
    form_bin = object%form_dir // '/bin/tform'
    qgraf_bin = object%qgraf_dir // '/bin/qgraf'
    if (object%gosam_dir /= "") then
       haggies_bin = '/usr/bin/java -jar ' // object%gosam_dir // &
            '/share/golem/haggies/haggies.jar'
    else
       call msg_fatal ("generate_configuration_file: At least " // &
            "the GoSam Directory has to be specified!")
    end if
    if (object%golem_dir /= "") then
       fcflags_golem = "-I" // object%golem_dir // "/include/golem95"
       ldflags_golem = "-L" // object%golem_dir // "/lib -lgolem"
    end if
    if (object%samurai_dir /= "") then
       fcflags_samurai = "-I" // object%samurai_dir // "/include/samurai"
       ldflags_samurai = "-L" // object%samurai_dir // "/lib -lsamurai"
       ldflags_avh_olo = "-L" // object%samurai_dir // "/lib -lavh_olo"
       ldflags_qcdloop = "-L" // object%samurai_dir // "/lib -lqcdloop"
    end if
    if (object%ninja_dir /= "") then
       fcflags_ninja = "-I" // object%ninja_dir // "/include/ninja " &
            // "-I" // object%ninja_dir // "/include"
       ldflags_ninja = "-L" // object%ninja_dir // "/lib -lninja"
    end if
    write (unit, "(A)") "#+avh_olo.ldflags=" &
         // char (ldflags_avh_olo)
    write (unit, "(A)") "reduction_programs=golem95, samurai, ninja"
    write (unit, "(A)") "extensions=autotools"
    write (unit, "(A)") "#+qcdloop.ldflags=" &
         // char (ldflags_qcdloop)
    write (unit, "(A)") "#+zzz.extensions=qcdloop, avh_olo"
    write (unit, "(A)") "#fc.bin=" // char (fc_bin)
    write (unit, "(A)") "form.bin=" // char (form_bin)
    write (unit, "(A)") "qgraf.bin=" // char (qgraf_bin)
    write (unit, "(A)") "#golem95.fcflags=" // char (fcflags_golem)
    write (unit, "(A)") "#golem95.ldflags=" // char (ldflags_golem)
    write (unit, "(A)") "haggies.bin=" // char (haggies_bin)
    write (unit, "(A)") "#samurai.fcflags=" // char (fcflags_samurai)
    write (unit, "(A)") "#samurai.ldflags=" // char (ldflags_samurai)
    write (unit, "(A)") "#ninja.fcflags=" // char (fcflags_ninja)
    write (unit, "(A)") "#ninja.ldflags=" // char (ldflags_ninja)
    !!! This might collide with the mass-setup in the order-file
    !!! write (unit, "(A)") "zero=mU,mD,mC,mS,mB"
    !!! This is covered by the BLHA2 interface
    write (unit, "(A)") "PSP_check=False"
      if (char (object%filter_lo) /= "") &
         write (unit, "(A)") "filter.lo=" // char (object%filter_lo)
      if (char (object%filter_nlo) /= "") &
         write (unit, "(A)") "filter.nlo=" // char (object%filter_nlo)
      if (char (object%symmetries) /= "") &
         write (unit, "(A)") "symmetries=" // char(object%symmetries)
      write (unit, "(A,I0)") "form.threads=", object%form_threads
      write (unit, "(A,I0)") "form.workspace=", object%form_workspace
      if (char (object%fc) /= "") &
         write (unit, "(A)") "fc.bin=" // char(object%fc)
  end subroutine gosam_writer_generate_configuration_file

  module subroutine gosam_driver_write_makefile (object, unit, libname)
    class(gosam_driver_t), intent(in) :: object
    integer, intent(in) :: unit
    type(string_t), intent(in) :: libname
    write (unit, "(2A)")  "OLP_FILE = ", char (object%olp_file)
    write (unit, "(2A)")  "OLP_DIR = ", char (object%olp_dir)
    write (unit, "(A)")
    write (unit, "(A)")   "all: $(OLP_DIR)/config.log"
    write (unit, "(2A)")  TAB, "make -C $(OLP_DIR) install"
    write (unit, "(A)")
    write (unit, "(3A)")  "$(OLP_DIR)/config.log: "
    write (unit, "(4A)")  TAB, char (object%gosam_dir // "/bin/gosam.py "), &
                             "--olp $(OLP_FILE) --destination=$(OLP_DIR)", &
                             " -f -z"
    write (unit, "(3A)")  TAB, "cd $(OLP_DIR); ./autogen.sh --prefix=", &
         "$(dir $(abspath $(lastword $(MAKEFILE_LIST))))"
  end subroutine gosam_driver_write_makefile
  module subroutine gosam_driver_set_alpha_s (driver, alpha_s)
     class(gosam_driver_t), intent(in) :: driver
     real(default), intent(in) :: alpha_s
     integer :: ierr
     call driver%blha_olp_set_parameter &
              (c_char_'alphaS'//c_null_char, &
               dble (alpha_s), 0._double, ierr)
  end subroutine gosam_driver_set_alpha_s

  module subroutine gosam_driver_set_alpha_qed (driver, alpha)
    class(gosam_driver_t), intent(inout) :: driver
    real(default), intent(in) :: alpha
    integer :: ierr
    call driver%blha_olp_set_parameter &
       (c_char_'alpha'//c_null_char, &
        dble (alpha), 0._double, ierr)
    if (ierr == 0) call ew_parameter_error_message (var_str ('alpha'))
  end subroutine gosam_driver_set_alpha_qed

  module subroutine gosam_driver_set_GF (driver, GF)
    class(gosam_driver_t), intent(inout) :: driver
    real(default), intent(in) :: GF
    integer :: ierr
    call driver%blha_olp_set_parameter &
       (c_char_'GF'//c_null_char, &
        dble(GF), 0._double, ierr)
    if (ierr == 0) call ew_parameter_error_message (var_str ('GF'))
  end subroutine gosam_driver_set_GF

  module subroutine gosam_driver_set_weinberg_angle (driver, sw2)
    class(gosam_driver_t), intent(inout) :: driver
    real(default), intent(in) :: sw2
    integer :: ierr
    call driver%blha_olp_set_parameter &
       (c_char_'sw2'//c_null_char, &
        dble(sw2), 0._double, ierr)
    if (ierr == 0) call ew_parameter_error_message (var_str ('sw2'))
  end subroutine gosam_driver_set_weinberg_angle

  module subroutine gosam_driver_print_alpha_s (object)
    class(gosam_driver_t), intent(in) :: object
    call object%blha_olp_print_parameter (c_char_'alphaS'//c_null_char)
  end subroutine gosam_driver_print_alpha_s

  module subroutine prc_gosam_prepare_library (object, os_data, libname)
    class(prc_gosam_t), intent(inout) :: object
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in) :: libname
    select type (writer => object%def%writer)
    type is (gosam_writer_t)
       call writer%write_config ()
    end select
    call object%create_olp_library (libname)
    call object%load_driver (os_data)
  end subroutine prc_gosam_prepare_library

  module subroutine prc_gosam_prepare_external_code &
       (core, flv_states, var_list, os_data, libname, model, i_core, is_nlo)
    class(prc_gosam_t), intent(inout) :: core
    integer, intent(in), dimension(:,:), allocatable :: flv_states
    type(var_list_t), intent(in) :: var_list
    type(os_data_t), intent(in) :: os_data
    type(string_t), intent(in) :: libname
    type(model_data_t), intent(in), target :: model
    integer, intent(in) :: i_core
    logical, intent(in) :: is_nlo
    core%sqme_tree_pos = 4
    call core%prepare_library (os_data, libname)
    call core%start ()
    call core%read_contract_file (flv_states)
    call core%set_particle_properties (model)
    call core%set_electroweak_parameters (model)
    call core%print_parameter_file (i_core)
  end subroutine prc_gosam_prepare_external_code

  module subroutine prc_gosam_write_makefile (object, unit, libname)
    class(prc_gosam_t), intent(in) :: object
    integer, intent(in) :: unit
    type(string_t), intent(in) :: libname
    select type (driver => object%driver)
    type is (gosam_driver_t)
       call driver%write_makefile (unit, libname)
    end select
  end subroutine prc_gosam_write_makefile

  module subroutine prc_gosam_execute_makefile (object, libname)
    class(prc_gosam_t), intent(in) :: object
    type(string_t), intent(in) :: libname
    select type (driver => object%driver)
    type is (gosam_driver_t)
       call os_system_call ("make -f " // &
          libname // "_gosam.makefile")
    end select
  end subroutine prc_gosam_execute_makefile

  module subroutine prc_gosam_create_olp_library (object, libname)
    class(prc_gosam_t), intent(inout) :: object
    type(string_t), intent(in) :: libname
    integer :: unit
    select type (driver => object%driver)
    type is (gosam_driver_t)
       unit = free_unit ()
       open (unit, file = char (libname // "_gosam.makefile"), &
            status = "replace", action= "write")
       call object%write_makefile (unit, libname)
       close (unit)
       call object%execute_makefile (libname)
    end select
  end subroutine prc_gosam_create_olp_library

  module subroutine prc_gosam_load_driver (object, os_data)
    class(prc_gosam_t), intent(inout) :: object
    type(os_data_t), intent(in) :: os_data
    logical :: dl_success

    select type (driver => object%driver)
    type is (gosam_driver_t)
       call driver%load (os_data, dl_success)
       if (.not. dl_success) &
          call msg_fatal ("GoSam Libraries could not be loaded")
    end select
  end subroutine prc_gosam_load_driver

  module subroutine prc_gosam_start (object)
    class(prc_gosam_t), intent(inout) :: object
    integer :: ierr
    if (object%includes_polarization())  &
         call msg_fatal ('GoSam does not support polarized beams!')
    select type (driver => object%driver)
    type is (gosam_driver_t)
       call driver%blha_olp_start (string_f2c (driver%contract_file), ierr)
    end select
  end subroutine prc_gosam_start

  module subroutine prc_gosam_write (object, unit)
    class(prc_gosam_t), intent(in) :: object
    integer, intent(in), optional :: unit
    call msg_message (unit = unit, string = "GOSAM")
  end subroutine prc_gosam_write

  module subroutine prc_gosam_write_name (object, unit)
    class(prc_gosam_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u,"(1x,A)") "Core: GoSam"
  end subroutine prc_gosam_write_name

  module subroutine prc_gosam_init_driver (object, os_data)
    class(prc_gosam_t), intent(inout) :: object
    type(os_data_t), intent(in) :: os_data
    type(string_t) :: olp_file, olc_file, olp_dir
    type(string_t) :: suffix

    select type (def => object%def)
    type is (gosam_def_t)
       suffix = def%suffix
       olp_file = def%basename // suffix // '.olp'
       olc_file = def%basename // suffix // '.olc'
       olp_dir = def%basename // suffix // '_olp_modules'
    class default
       call msg_bug ("prc_gosam_init_driver: core_def should be of gosam-type")
    end select

    select type(driver => object%driver)
    type is (gosam_driver_t)
       driver%nlo_suffix = suffix
       call driver%init_gosam (os_data, olp_file, olc_file, olp_dir, &
            var_str ("libgolem_olp"))
    end select
  end subroutine prc_gosam_init_driver

  module subroutine prc_gosam_set_initialized (prc_gosam)
    class(prc_gosam_t), intent(inout) :: prc_gosam
    prc_gosam%initialized = .true.
  end subroutine prc_gosam_set_initialized

  module subroutine prc_gosam_compute_sqme_spin_c (object, &
       i_flv, i_hel, em, p, ren_scale, me_sc, bad_point)
    class(prc_gosam_t), intent(inout) :: object
    integer, intent(in) :: i_flv, i_hel
    integer, intent(in) :: em
    type(vector4_t), intent(in), dimension(:) :: p
    real(default), intent(in) :: ren_scale
    complex(default), intent(out) :: me_sc
    logical, intent(out) :: bad_point
    real(double), dimension(5*object%n_particles) :: mom
    real(double), dimension(OLP_RESULTS_LIMIT) :: r
    real(double) :: ren_scale_dble
    integer :: i, igm1, n
    integer :: pos_real, pos_imag
    real(double) :: acc_dble
    real(default) :: acc, alpha_s
    if (object%i_spin_c(i_flv, i_hel) >= 0) then
       me_sc = cmplx (zero ,zero, kind=default)
       mom = object%create_momentum_array (p)
       if (vanishes (ren_scale)) &
          call msg_fatal ("prc_gosam_compute_sqme_spin_c: ren_scale vanishes")
       alpha_s = object%qcd%alpha%get (ren_scale)
       ren_scale_dble = dble (ren_scale)
       select type (driver => object%driver)
       type is (gosam_driver_t)
          call driver%set_alpha_s (alpha_s)
          call driver%blha_olp_eval2 (object%i_spin_c(i_flv, i_hel), &
               mom, ren_scale_dble, r, acc_dble)
       end select
       igm1 = em - 1
       n = size(p)
       do i = 0, n - 1
          pos_real = 2 * igm1 + 2 * n * i + 1
          pos_imag = pos_real + 1
          me_sc = me_sc + cmplx (r(pos_real), r(pos_imag), default)
       end do

       me_sc = - conjg(me_sc) / CA

       acc = acc_dble
       if (acc > object%maximum_accuracy) bad_point = .true.
    else
       r = 0._double
    end if
  end subroutine prc_gosam_compute_sqme_spin_c

  module subroutine gosam_state_write (object, unit)
    class(gosam_state_t), intent(in) :: object
    integer, intent(in), optional :: unit
    call msg_warning (unit = unit, string = "gosam_state_write: What to write?")
  end subroutine gosam_state_write

  module subroutine prc_gosam_set_particle_properties (object, model)
    class(prc_gosam_t), intent(inout) :: object
    class(model_data_t), intent(in), target :: model
    integer :: i, i_pdg
    type(flavor_t) :: flv
    real(default) :: mass, width
    integer :: ierr
    real(default) :: top_yukawa
    do i = 1, OLP_N_MASSIVE_PARTICLES
       i_pdg = OLP_MASSIVE_PARTICLES(i)
       if (i_pdg < 0) cycle
       call flv%init (i_pdg, model)
       mass = flv%get_mass (); width = flv%get_width ()
       select type (driver => object%driver)
       class is (blha_driver_t)
          if (i_pdg == 13) then
             call driver%set_mass_and_width (i_pdg, mass = mass)
          else
             call driver%set_mass_and_width (i_pdg, mass = mass, width = width)
          end if
          if (i_pdg == 5) call driver%blha_olp_set_parameter &
             ('yuk(5)'//c_null_char, dble(mass), 0._double, ierr)
          if (i_pdg == 6) then
             if (driver%external_top_yukawa > 0._default) then
                top_yukawa = driver%external_top_yukawa
             else
                top_yukawa = mass
             end if
             call driver%blha_olp_set_parameter &
                ('yuk(6)'//c_null_char, dble(top_yukawa), 0._double, ierr)
          end if
          if (driver%switch_off_muon_yukawas) then
             if (i_pdg == 13) call driver%blha_olp_set_parameter &
                ('yuk(13)' //c_null_char, 0._double, 0._double, ierr)
          end if
       end select
    end do
  end subroutine prc_gosam_set_particle_properties


end submodule prc_gosam_s

