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

module cascades2_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use numeric_utils

  use cascades2
  use flavors
  use phs_forests, only: phs_parameters_t
  use model_data

  implicit none
  private

  public :: cascades2_1
  public :: cascades2_2

contains

  subroutine cascades2_1 (u)
    integer, intent(in) :: u
    type(feyngraph_set_t) :: feyngraph_set
    type(model_data_t) :: model
    integer :: n_in = 1
    integer :: n_out = 6
    type(flavor_t), dimension(7,1) :: flv
    type(phs_parameters_t) :: phs_par
    logical :: fatal_beam_decay = .true.
    integer :: u_in = 8

    write (u, "(A)")  "* Test output: cascades2_1"
    write (u, "(A)")  "*   Purpose: create a test phs file (decay) with the forest"
    write (u, "(A)")  "*            output of O'Mega"
    write (u, "(A)")

    write (u, "(A)")  "* Initializing"
    write (u, "(A)")

    call init_sm_full_test (model)

    call flv(1,1)%init (6, model)
    call flv(2,1)%init (5, model)
    call flv(3,1)%init (-11, model)
    call flv(4,1)%init (12, model)
    call flv(5,1)%init (21, model)
    call flv(6,1)%init (22, model)
    call flv(7,1)%init (21, model)

    phs_par%sqrts = 173.1_default
    phs_par%m_threshold_s = 50._default
    phs_par%m_threshold_t = 100._default
    phs_par%keep_nonresonant = .true.
    phs_par%off_shell = 2

    open (unit=u_in, file="cascades2_1.fds", status='old', action='read')

    write (u, "(A)")
    write (u, "(A)")  "* Generating phase-space parametrizations"
    write (u, "(A)")

    call feyngraph_set_generate (feyngraph_set, model, n_in, n_out, &
         flv, phs_par, fatal_beam_decay, u_in, use_dag = .false., &
         vis_channels = .false.)
    call feyngraph_set_write_process_bincode_format (feyngraph_set, u)
    call feyngraph_set_write_file_format (feyngraph_set, u)

    write (u, "(A)")  "* Cleanup"
    write (u, "(A)")

    close (u_in)
    call feyngraph_set%final ()
    call model%final ()

    write (u, *)
    write (u, "(A)")  "* Test output end: cascades2_1"
  end subroutine cascades2_1

  subroutine cascades2_2 (u)
    integer, intent(in) :: u
    type(feyngraph_set_t) :: feyngraph_set
    type(model_data_t) :: model
    integer :: n_in = 2
    integer :: n_out = 5
    type(flavor_t), dimension(7,1) :: flv
    type(phs_parameters_t) :: phs_par
    logical :: fatal_beam_decay = .true.
    integer :: u_in = 8

    write (u, "(A)")  "* Test output: cascades2_2"
    write (u, "(A)")  "*   Purpose: create a test phs file (scattering) with the"
    write (u, "(A)")  "*            parsable DAG output of O'Mega"
    write (u, "(A)")

    write (u, "(A)")  "* Initializing"
    write (u, "(A)")

    call init_sm_full_test (model)

    call flv(1,1)%init (-11, model)
    call flv(2,1)%init (11, model)
    call flv(3,1)%init (-11, model)
    call flv(4,1)%init (12, model)
    call flv(5,1)%init (1, model)
    call flv(6,1)%init (-2, model)
    call flv(7,1)%init (22, model)

    phs_par%sqrts = 500._default
    phs_par%m_threshold_s = 50._default
    phs_par%m_threshold_t = 100._default
    phs_par%keep_nonresonant = .true.
    phs_par%off_shell = 2
    phs_par%t_channel = 6

    open (unit=u_in, file="cascades2_2.fds", &
         status='old', action='read')

    write (u, "(A)")
    write (u, "(A)")  "* Generating phase-space parametrizations"
    write (u, "(A)")

    call feyngraph_set_generate (feyngraph_set, model, n_in, n_out, &
         flv, phs_par, fatal_beam_decay, u_in, use_dag = .true., &
         vis_channels = .false.)
    call feyngraph_set_write_process_bincode_format (feyngraph_set, u)
    call feyngraph_set_write_file_format (feyngraph_set, u)

    write (u, "(A)")  "* Cleanup"
    write (u, "(A)")

    close (u_in)
    call feyngraph_set%final ()
    call model%final ()

    write (u, *)
    write (u, "(A)")  "* Test output end: cascades2_2"
  end subroutine cascades2_2


end module cascades2_uti
