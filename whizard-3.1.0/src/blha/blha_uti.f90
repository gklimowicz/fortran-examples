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

module blha_uti

  use iso_varying_string, string_t => varying_string
  use format_utils, only: write_separator
  use variables, only: var_list_t
  use os_interface
  use models
  use blha_config


  implicit none
  private

  public :: blha_1
  public :: blha_2
  public :: blha_3

contains

  subroutine setup_and_write_blha_configuration (u, single, polarized)
    integer, intent(in) :: u
    logical, intent(in), optional :: single
    logical, intent(in), optional :: polarized
    logical :: polrzd, singl
    type(blha_master_t) :: blha_master
    integer :: i
    integer :: n_in, n_out
    integer :: alpha_power, alphas_power
    integer, dimension(:,:), allocatable :: flv_born, flv_real
    type(string_t) :: proc_id, method, correction_type
    type(os_data_t) :: os_data
    type(model_list_t) :: model_list
    type(var_list_t) :: var_list
    type(model_t), pointer :: model => null ()
    integer :: openloops_phs_tolerance

    polrzd = .false.; if (present (polarized)) polrzd = polarized
    singl = .true.; if (present (single)) singl = single

    if (singl) then
       write (u, "(A)") "* Process: e+ e- -> W+ W- b b~"
       n_in = 2; n_out = 4
       alpha_power = 4; alphas_power = 0
       allocate (flv_born (n_in + n_out, 1))
       allocate (flv_real (n_in + n_out + 1, 1))
       flv_born(1,1) = 11; flv_born(2,1) = -11
       flv_born(3,1) = 24; flv_born(4,1) = -24
       flv_born(5,1) = 5; flv_born(6,1) = -5
       flv_real(1:6,1) = flv_born(:,1)
       flv_real(7,1) = 21
    else
       write (u, "(A)") "* Process: e+ e- -> u:d:s U:D:S"
       n_in = 2; n_out = 2
       alpha_power = 2; alphas_power = 0
       allocate (flv_born (n_in + n_out, 3))
       allocate (flv_real (n_in + n_out + 1, 3))
       flv_born(1,:) = 11; flv_born(2,:) = -11
       flv_born(3,1) = 1; flv_born(4,1) = -1
       flv_born(3,2) = 2; flv_born(4,2) = -2
       flv_born(3,3) = 3; flv_born(4,3) = -3
       flv_real(1:4,:) = flv_born
       flv_real(5,:) = 21
    end if
    proc_id = var_str ("BLHA_Test")

    call syntax_model_file_init ()
    call os_data%init ()
    call model_list%read_model &
         (var_str ("SM"), var_str ("SM.mdl"), os_data, model)

    write (u, "(A)") "* BLHA matrix elements assumed for all process components"
    write (u, "(A)") "* Mode: GoSam"

    method = var_str ("gosam")
    correction_type = var_str ("QCD")
    call var_list%append_string (var_str ("$born_me_method"), method)
    call var_list%append_string (var_str ("$real_tree_me_method"), method)
    call var_list%append_string (var_str ("$loop_me_method"), method)
    call var_list%append_string (var_str ("$correlation_me_method"), method)
    call blha_master%set_ew_scheme (var_str ("GF"))
    call blha_master%set_methods (.true., var_list)
    call blha_master%allocate_config_files ()
    call blha_master%set_correction_type (correction_type)
    call blha_master%generate (proc_id, model, n_in, &
         alpha_power, alphas_power, flv_born, flv_real)

    call test_output (u)

    call blha_master%final ()
    call var_list%final ()
    write (u, "(A)") "* Switch to OpenLoops"
    openloops_phs_tolerance = 7

    method = var_str ("openloops")
    correction_type = var_str ("QCD")
    call var_list%append_string (var_str ("$born_me_method"), method)
    call var_list%append_string (var_str ("$real_tree_me_method"), method)
    call var_list%append_string (var_str ("$loop_me_method"), method)
    call var_list%append_string (var_str ("$correlation_me_method"), method)
    call blha_master%set_methods (.true., var_list)
    call blha_master%allocate_config_files ()
    call blha_master%set_correction_type (correction_type)
    call blha_master%generate (proc_id, model, n_in, &
         alpha_power, alphas_power, flv_born, flv_real)

    if (polrzd) then
       do i = 1, 4
          call blha_master%set_polarization (i)
       end do
    end if
    call blha_master%setup_additional_features &
         (openloops_phs_tolerance, .false., 0)

    call test_output (u)

  contains

    subroutine test_output (u)
      integer, intent(in) :: u
      do i = 1, 4
         call write_separator (u)
         call write_component_type (i, u)
         call write_separator (u)
         call blha_configuration_write &
              (blha_master%blha_cfg(i), blha_master%suffix(i), u, no_version = .true.)
      end do
    end subroutine test_output

    subroutine write_component_type (i, u)
      integer, intent(in) :: i, u
      type(string_t) :: message, component_type
      message = var_str ("OLP-File content for ")
      select case (i)
      case (1)
         component_type = var_str ("loop")
      case (2)
         component_type = var_str ("subtraction")
      case (3)
         component_type = var_str ("real")
      case (4)
         component_type = var_str ("born")
      end select
      message = message // component_type // " matrix elements"
      write (u, "(A)") char (message)
    end subroutine write_component_type

  end subroutine setup_and_write_blha_configuration


  subroutine blha_1 (u)
    integer, intent(in) :: u
    write (u, "(A)") "* Test output: blha_1"
    write (u, "(A)") "* Purpose: Test the creation of olp-files for single "&
         &"and unpolarized flavor structures"
    write (u, "(A)")
    call setup_and_write_blha_configuration (u, single = .true., polarized = .false.)
  end subroutine blha_1

  subroutine blha_2 (u)
    integer, intent(in) :: u
    write (u, "(A)") "* Test output: blha_2"
    write (u, "(A)") "* Purpose: Test the creation of olp-files for multiple "&
         &"and unpolarized flavor structures"
    write (u, "(A)")
    call setup_and_write_blha_configuration (u, single = .false., polarized = .false.)
  end subroutine blha_2

  subroutine blha_3 (u)
    integer, intent(in) :: u
    write (u, "(A)") "* Test output: blha_3"
    write (u, "(A)") "* Purpose: Test the creation of olp-files for single "&
         &"and polarized flavor structures"
    write (u, "(A)")
    call setup_and_write_blha_configuration (u, single = .true., polarized = .true.)
  end subroutine blha_3


end module blha_uti

