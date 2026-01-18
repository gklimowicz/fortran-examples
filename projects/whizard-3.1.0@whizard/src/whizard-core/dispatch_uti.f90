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

module dispatch_uti
  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface, only: os_data_t
  use physics_defs, only: ELECTRON, PROTON
  use sm_qcd, only: qcd_t
  use flavors, only: flavor_t
  use interactions, only: reset_interaction_counter
  use pdg_arrays, only: pdg_array_t, assignment(=)
  use prc_core_def, only: prc_core_def_t
  use prc_test_core, only: test_t
  use prc_core, only: prc_core_t
  use prc_test, only: prc_test_def_t
  use prc_omega, only: omega_def_t, prc_omega_t
  use sf_mappings, only: sf_channel_t
  use sf_base, only: sf_data_t, sf_config_t
  use phs_base, only: phs_channel_collection_t
  use variables, only: var_list_t
  use model_data, only: model_data_t
  use models, only: syntax_model_file_init, syntax_model_file_final
  use rt_data, only: rt_data_t

  use dispatch_phase_space, only: dispatch_sf_channels
  use dispatch_beams, only: sf_prop_t, dispatch_qcd
  use dispatch_beams, only: dispatch_sf_config, dispatch_sf_data
  use dispatch_me_methods, only: dispatch_core_def, dispatch_core
  use dispatch_me_methods, only: dispatch_core_update, dispatch_core_restore

  use sf_base_ut, only: sf_test_data_t

  implicit none
  private

  public :: dispatch_sf_data_test

  public :: dispatch_1
  public :: dispatch_2
  public :: dispatch_7
  public :: dispatch_8
  public :: dispatch_10
  public :: dispatch_11

contains

  subroutine dispatch_1 (u)
    integer, intent(in) :: u
    type(string_t), dimension(2) :: prt_in, prt_out
    type(rt_data_t), target :: global
    class(prc_core_def_t), allocatable :: core_def

    write (u, "(A)")  "* Test output: dispatch_1"
    write (u, "(A)")  "*   Purpose: select process configuration method"
    write (u, "(A)")

    call global%global_init ()

    call global%set_log (var_str ("?omega_openmp"), &
         .false., is_known = .true.)

    prt_in = [var_str ("a"), var_str ("b")]
    prt_out = [var_str ("c"), var_str ("d")]

    write (u, "(A)")  "* Allocate core_def as prc_test_def"

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call dispatch_core_def (core_def, prt_in, prt_out, global%model, global%var_list)
    select type (core_def)
    type is (prc_test_def_t)
       call core_def%write (u)
    end select

    deallocate (core_def)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate core_def as omega_def"
    write (u, "(A)")

    call global%set_string (var_str ("$method"), &
         var_str ("omega"), is_known = .true.)
    call dispatch_core_def (core_def, prt_in, prt_out, global%model, global%var_list)
    select type (core_def)
    type is (omega_def_t)
       call core_def%write (u)
    end select

    call global%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_1"

  end subroutine dispatch_1

  subroutine dispatch_2 (u)
    integer, intent(in) :: u
    type(string_t), dimension(2) :: prt_in, prt_out
    type(rt_data_t), target :: global
    class(prc_core_def_t), allocatable :: core_def
    class(prc_core_t), allocatable :: core

    write (u, "(A)")  "* Test output: dispatch_2"
    write (u, "(A)")  "*   Purpose: select process configuration method"
    write (u, "(A)")  "             and allocate process core"
    write (u, "(A)")

    call syntax_model_file_init ()
    call global%global_init ()

    prt_in = [var_str ("a"), var_str ("b")]
    prt_out = [var_str ("c"), var_str ("d")]

    write (u, "(A)")  "* Allocate core as test_t"
    write (u, "(A)")

    call global%set_string (var_str ("$method"), &
         var_str ("unit_test"), is_known = .true.)
    call dispatch_core_def (core_def, prt_in, prt_out, global%model, global%var_list)
    call dispatch_core (core, core_def)
    select type (core)
    type is (test_t)
       call core%write (u)
    end select

    deallocate (core)
    deallocate (core_def)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate core as prc_omega_t"
    write (u, "(A)")

    call global%set_string (var_str ("$method"), &
         var_str ("omega"), is_known = .true.)
    call dispatch_core_def (core_def, prt_in, prt_out, global%model, global%var_list)

    call global%select_model (var_str ("Test"))

    call global%set_log (&
         var_str ("?helicity_selection_active"), &
         .true., is_known = .true.)
    call global%set_real (&
         var_str ("helicity_selection_threshold"), &
         1e9_default, is_known = .true.)
    call global%set_int (&
         var_str ("helicity_selection_cutoff"), &
         10, is_known = .true.)

    call dispatch_core (core, core_def, &
         global%model, &
         global%get_helicity_selection ())
    call core_def%allocate_driver (core%driver, var_str (""))

    select type (core)
    type is (prc_omega_t)
       call core%write (u)
    end select

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_2"

  end subroutine dispatch_2

  subroutine dispatch_7 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global
    type(os_data_t) :: os_data
    type(string_t) :: prt, sf_method
    type(sf_prop_t) :: sf_prop
    class(sf_data_t), allocatable :: data
    type(pdg_array_t), dimension(1) :: pdg_in
    type(pdg_array_t), dimension(1,1) :: pdg_prc
    type(pdg_array_t), dimension(1) :: pdg_out
    integer, dimension(:), allocatable :: pdg1

    write (u, "(A)")  "* Test output: dispatch_7"
    write (u, "(A)")  "*   Purpose: select and configure &
         &structure function data"
    write (u, "(A)")

    call global%global_init ()

    call os_data%init ()
    call syntax_model_file_init ()
    call global%select_model (var_str ("QCD"))

    call reset_interaction_counter ()
    call global%set_real (var_str ("sqrts"), &
         14000._default, is_known = .true.)
    prt = "p"
    call global%beam_structure%init_sf ([prt, prt], [1])
    pdg_in = 2212

    write (u, "(A)")  "* Allocate data as sf_pdf_builtin_t"
    write (u, "(A)")

    sf_method = "pdf_builtin"
    call dispatch_sf_data (data, sf_method, [1], sf_prop, &
         global%get_var_list_ptr (), global%var_list, &
         global%model, global%os_data, global%get_sqrts (), &
         pdg_in, pdg_prc, .false.)
    call data%write (u)

    call data%get_pdg_out (pdg_out)
    pdg1 = pdg_out(1)
    write (u, "(A)")
    write (u, "(1x,A,99(1x,I0))")  "PDG(out) = ", pdg1

    deallocate (data)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate data for different PDF set"
    write (u, "(A)")

    pdg_in = 2212

    call global%set_string (var_str ("$pdf_builtin_set"), &
         var_str ("CTEQ6M"), is_known = .true.)
    sf_method = "pdf_builtin"
    call dispatch_sf_data (data, sf_method, [1], sf_prop, &
         global%get_var_list_ptr (), global%var_list, &
         global%model, global%os_data, global%get_sqrts (), &
         pdg_in, pdg_prc, .false.)
    call data%write (u)

    call data%get_pdg_out (pdg_out)
    pdg1 = pdg_out(1)
    write (u, "(A)")
    write (u, "(1x,A,99(1x,I0))")  "PDG(out) = ", pdg1

    deallocate (data)

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_7"

  end subroutine dispatch_7

  subroutine dispatch_8 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global
    type(os_data_t) :: os_data
    type(flavor_t), dimension(2) :: flv
    type(sf_config_t), dimension(:), allocatable :: sf_config
    type(sf_prop_t) :: sf_prop
    type(sf_channel_t), dimension(:), allocatable :: sf_channel
    type(phs_channel_collection_t) :: coll
    type(string_t) :: sf_string
    integer :: i
    type(pdg_array_t), dimension (2,1) :: pdg_prc

    write (u, "(A)")  "* Test output: dispatch_8"
    write (u, "(A)")  "*   Purpose: configure a structure-function chain"
    write (u, "(A)")

    call global%global_init ()

    call os_data%init ()
    call syntax_model_file_init ()
    call global%select_model (var_str ("QCD"))

    write (u, "(A)")  "* Allocate LHC beams with PDF builtin"
    write (u, "(A)")

    call flv(1)%init (PROTON, global%model)
    call flv(2)%init (PROTON, global%model)

    call reset_interaction_counter ()
    call global%set_real (var_str ("sqrts"), &
         14000._default, is_known = .true.)

    call global%beam_structure%init_sf (flv%get_name (), [1])
    call global%beam_structure%set_sf (1, 1, var_str ("pdf_builtin"))

    call dispatch_sf_config (sf_config, sf_prop, global%beam_structure, &
         global%get_var_list_ptr (), global%var_list, &
         global%model, global%os_data, global%get_sqrts (), pdg_prc)
    do i = 1, size (sf_config)
       call sf_config(i)%write (u)
    end do

    call dispatch_sf_channels (sf_channel, sf_string, sf_prop, coll, &
         global%var_list, global%get_sqrts(), global%beam_structure)
    write (u, "(1x,A)")  "Mapping configuration:"
    do i = 1, size (sf_channel)
       write (u, "(2x)", advance = "no")
       call sf_channel(i)%write (u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Allocate ILC beams with CIRCE1"
    write (u, "(A)")

    call global%select_model (var_str ("QED"))
    call flv(1)%init ( ELECTRON, global%model)
    call flv(2)%init (-ELECTRON, global%model)

    call reset_interaction_counter ()
    call global%set_real (var_str ("sqrts"), &
         500._default, is_known = .true.)
    call global%set_log (var_str ("?circe1_generate"), &
         .false., is_known = .true.)

    call global%beam_structure%init_sf (flv%get_name (), [1])
    call global%beam_structure%set_sf (1, 1, var_str ("circe1"))

    call dispatch_sf_config (sf_config, sf_prop, global%beam_structure, &
         global%get_var_list_ptr (), global%var_list, &
         global%model, global%os_data, global%get_sqrts (), pdg_prc)
    do i = 1, size (sf_config)
       call sf_config(i)%write (u)
    end do

    call dispatch_sf_channels (sf_channel, sf_string, sf_prop, coll, &
         global%var_list, global%get_sqrts(), global%beam_structure)
    write (u, "(1x,A)")  "Mapping configuration:"
    do i = 1, size (sf_channel)
       write (u, "(2x)", advance = "no")
       call sf_channel(i)%write (u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_8"

  end subroutine dispatch_8

  subroutine dispatch_10 (u)
    integer, intent(in) :: u
    type(string_t), dimension(2) :: prt_in, prt_out
    type(rt_data_t), target :: global
    class(prc_core_def_t), allocatable :: core_def
    class(prc_core_t), allocatable :: core, saved_core
    type(var_list_t), pointer :: model_vars

    write (u, "(A)")  "* Test output: dispatch_10"
    write (u, "(A)")  "*   Purpose: select process configuration method,"
    write (u, "(A)")  "             allocate process core,"
    write (u, "(A)")  "             temporarily reset parameters"
    write (u, "(A)")

    call syntax_model_file_init ()
    call global%global_init ()

    prt_in = [var_str ("a"), var_str ("b")]
    prt_out = [var_str ("c"), var_str ("d")]

    write (u, "(A)")  "* Allocate core as prc_omega_t"
    write (u, "(A)")

    call global%set_string (var_str ("$method"), &
         var_str ("omega"), is_known = .true.)
    call dispatch_core_def (core_def, prt_in, prt_out, global%model, global%var_list)

    call global%select_model (var_str ("Test"))

    call dispatch_core (core, core_def, global%model)
    call core_def%allocate_driver (core%driver, var_str (""))

    select type (core)
    type is (prc_omega_t)
       call core%write (u)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Update core with modified model and helicity selection"
    write (u, "(A)")

    model_vars => global%model%get_var_list_ptr ()

    call model_vars%set_real (var_str ("gy"), 2._default, &
         is_known = .true.)
    call global%model%update_parameters ()

    call global%set_log (&
         var_str ("?helicity_selection_active"), &
         .true., is_known = .true.)
    call global%set_real (&
         var_str ("helicity_selection_threshold"), &
         2e10_default, is_known = .true.)
    call global%set_int (&
         var_str ("helicity_selection_cutoff"), &
         5, is_known = .true.)

    call dispatch_core_update (core, &
         global%model, &
         global%get_helicity_selection (), &
         saved_core = saved_core)
    select type (core)
    type is (prc_omega_t)
       call core%write (u)
    end select

    write (u, "(A)")
    write (u, "(A)")  "* Restore core from save"
    write (u, "(A)")

    call dispatch_core_restore (core, saved_core)
    select type (core)
    type is (prc_omega_t)
       call core%write (u)
    end select

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_10"

  end subroutine dispatch_10

  subroutine dispatch_11 (u)
    integer, intent(in) :: u
    type(rt_data_t), target :: global
    type(var_list_t), pointer :: model_vars
    type(qcd_t) :: qcd

    write (u, "(A)")  "* Test output: dispatch_11"
    write (u, "(A)")  "*   Purpose: select QCD coupling formula"
    write (u, "(A)")

    call syntax_model_file_init ()
    call global%global_init ()
    call global%select_model (var_str ("SM"))
    model_vars => global%get_var_list_ptr ()

    write (u, "(A)")  "* Allocate alpha_s as fixed"
    write (u, "(A)")

    call global%set_log (var_str ("?alphas_is_fixed"), &
         .true., is_known = .true.)
    call dispatch_qcd (qcd, global%get_var_list_ptr (), global%os_data)
    call qcd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate alpha_s as running (built-in)"
    write (u, "(A)")

    call global%set_log (var_str ("?alphas_is_fixed"), &
         .false., is_known = .true.)
    call global%set_log (var_str ("?alphas_from_mz"), &
         .true., is_known = .true.)
    call global%set_int &
         (var_str ("alphas_order"), 1, is_known = .true.)
    call model_vars%set_real (var_str ("alphas"), 0.1234_default, &
          is_known=.true.)
    call model_vars%set_real (var_str ("mZ"), 91.234_default, &
          is_known=.true.)
    call dispatch_qcd (qcd, global%get_var_list_ptr (), global%os_data)
    call qcd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate alpha_s as running (built-in, Lambda defined)"
    write (u, "(A)")

    call global%set_log (var_str ("?alphas_from_mz"), &
         .false., is_known = .true.)
    call global%set_log (&
         var_str ("?alphas_from_lambda_qcd"), &
         .true., is_known = .true.)
    call global%set_real &
         (var_str ("lambda_qcd"), 250.e-3_default, &
          is_known=.true.)
    call global%set_int &
         (var_str ("alphas_order"), 2, is_known = .true.)
    call global%set_int &
         (var_str ("alphas_nf"), 4, is_known = .true.)
    call dispatch_qcd (qcd, global%get_var_list_ptr (), global%os_data)
    call qcd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Allocate alpha_s as running (using builtin PDF set)"
    write (u, "(A)")

    call global%set_log (&
         var_str ("?alphas_from_lambda_qcd"), &
         .false., is_known = .true.)
    call global%set_log &
         (var_str ("?alphas_from_pdf_builtin"), &
         .true., is_known = .true.)
    call dispatch_qcd (qcd, global%get_var_list_ptr (), global%os_data)
    call qcd%write (u)

    call global%final ()
    call syntax_model_file_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: dispatch_11"

  end subroutine dispatch_11


  subroutine dispatch_sf_data_test (data, sf_method, i_beam, sf_prop, &
       var_list, var_list_global, model, os_data, sqrts, pdg_in, pdg_prc, polarized)
    class(sf_data_t), allocatable, intent(inout) :: data
    type(string_t), intent(in) :: sf_method
    integer, dimension(:), intent(in) :: i_beam
    type(var_list_t), intent(in) :: var_list
    type(var_list_t), intent(inout) :: var_list_global
    class(model_data_t), target, intent(in) :: model
    type(os_data_t), intent(in) :: os_data
    real(default), intent(in) :: sqrts
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_in
    type(pdg_array_t), dimension(:,:), intent(in) :: pdg_prc
    type(sf_prop_t), intent(inout) :: sf_prop
    logical, intent(in) :: polarized
    select case (char (sf_method))
    case ("sf_test_0", "sf_test_1")
       allocate (sf_test_data_t :: data)
       select type (data)
       type is (sf_test_data_t)
          select case (char (sf_method))
          case ("sf_test_0");  call data%init (model, pdg_in(i_beam(1)))
          case ("sf_test_1");  call data%init (model, pdg_in(i_beam(1)),&
               mode = 1)
          end select
       end select
    end select
  end subroutine dispatch_sf_data_test


end module dispatch_uti

