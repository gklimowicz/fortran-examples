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

module parton_states_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use constants, only: zero
  use numeric_utils
  use flavors
  use colors
  use helicities
  use quantum_numbers
  use sf_base, only: sf_chain_instance_t
  use state_matrices, only: state_matrix_t
  use prc_template_me, only: prc_template_me_t
  use interactions, only: interaction_t
  use models, only: model_t, create_test_model
  use parton_states

  implicit none
  private

  public :: parton_states_1

contains

  subroutine parton_states_1 (u)
    integer, intent(in) :: u
    type(state_matrix_t), allocatable :: state
    type(flavor_t), dimension(2) :: flv_in
    type(flavor_t), dimension(2) :: flv_out1, flv_out2
    type(flavor_t), dimension(4) :: flv_tot
    type(helicity_t), dimension(4) :: hel
    type(color_t), dimension(4) :: col
    integer :: h1, h2, h3, h4
    integer :: f
    integer :: i
    type(quantum_numbers_t), dimension(4) :: qn
    type(prc_template_me_t) :: core
    type(sf_chain_instance_t), target :: sf_chain
    type(interaction_t), target :: int
    type(isolated_state_t) :: isolated_state
    integer :: n_states = 0
    integer, dimension(:), allocatable :: col_flow_index
    type(quantum_numbers_mask_t), dimension(2) :: qn_mask
    integer, dimension(8) :: i_allowed_states
    complex(default), dimension(8) :: me
    complex(default) :: me_check_tot, me_check_1, me_check_2, me2
    logical :: tmp1, tmp2
    type(model_t), pointer :: test_model => null ()

    write (u, "(A)") "* Test output: parton_states_1"
    write (u, "(A)") "* Purpose: Test the standard parton states"
    write (u, "(A)")

    call flv_in%init ([11, -11])
    call flv_out1%init ([1, -1])
    call flv_out2%init ([2, -2])

    write (u, "(A)") "* Using incoming flavors: "
    call flavor_write_array (flv_in, u)
    write (u, "(A)") "* Two outgoing flavor structures: "
    call flavor_write_array (flv_out1, u)
    call flavor_write_array (flv_out2, u)


    write (u, "(A)") "* Initialize state matrix"
    allocate (state)
    call state%init ()

    write (u, "(A)") "* Fill state matrix"
    call col(3)%init ([1])
    call col(4)%init ([-1])
    do f = 1, 2
       do h1 = -1, 1, 2
          do h2 = -1, 1, 2
             do h3 = -1, 1, 2
                do h4 = -1, 1, 2
                   n_states = n_states + 1
                   call hel%init ([h1, h2, h3, h4], [h1, h2, h3, h4])
                   if (f == 1) then
                      flv_tot = [flv_in, flv_out1]
                   else
                      flv_tot = [flv_in, flv_out2]
                   end if
                   call qn%init (flv_tot, col, hel)
                   call state%add_state (qn)
                end do
             end do
          end do
       end do
    end do

    !!! Two flavors, one color flow, 2 x 2 x 2 x 2 helicity configurations
    !!! -> 32 states.
    write (u, "(A)")
    write (u, "(A,I2)") "* Generated number of states: ", n_states

    call state%freeze ()

    !!! Indices of the helicity configurations which are non-zero
    i_allowed_states = [6, 7, 10, 11, 22, 23, 26, 27]
    me = [cmplx (-1.89448E-5_default,  9.94456E-7_default, default), &
          cmplx (-8.37887E-2_default,  4.30842E-3_default, default), &
          cmplx (-1.99997E-1_default, -1.01985E-2_default, default), &
          cmplx ( 1.79717E-5_default,  9.27038E-7_default, default), &
          cmplx (-1.74859E-5_default,  8.78819E-7_default, default), &
          cmplx ( 1.67577E-1_default, -8.61683E-3_default, default), &
          cmplx ( 2.41331E-1_default,  1.23306E-2_default, default), &
          cmplx (-3.59435E-5_default, -1.85407E-6_default, default)]
    me_check_tot = cmplx (zero, zero, default)
    me_check_1 = cmplx (zero, zero, default)
    me_check_2 = cmplx (zero, zero, default)
    do i = 1, 8
       me2 = me(i) * conjg (me(i))
       me_check_tot = me_check_tot + me2
       if (i < 5) then
          me_check_1 = me_check_1 + me2
       else
          me_check_2 = me_check_2 + me2
       end if
       call state%set_matrix_element (i_allowed_states(i), me(i))
    end do

    !!! Do not forget the color factor
    me_check_tot = 3._default * me_check_tot
    me_check_1 = 3._default * me_check_1
    me_check_2 = 3._default * me_check_2
    write (u, "(A)")

    write (u, "(A)") "* Setup interaction"
    call int%basic_init (2, 0, 2, set_relations = .true.)
    call int%set_state_matrix (state)

    core%data%n_in = 2; core%data%n_out = 2
    core%data%n_flv = 2
    allocate (core%data%flv_state (4, 2))
    core%data%flv_state (1, :) = [11, 11]
    core%data%flv_state (2, :) = [-11, -11]
    core%data%flv_state (3, :) = [1, 2]
    core%data%flv_state (4, :) = [-1, -2]
    core%use_color_factors = .false.
    core%nc = 3

    write (u, "(A)") "* Init isolated state"
    call isolated_state%init (sf_chain, int)
    !!! There is only one color flow.
    allocate (col_flow_index (n_states)); col_flow_index = 1
    call qn_mask%init (.false., .false., .true., mask_cg = .false.)
    write (u, "(A)") "* Give a trace to the isolated state"
    call isolated_state%setup_square_trace (core, qn_mask, col_flow_index, .false.)
    call isolated_state%evaluate_trace ()
    write (u, "(A)")
    write (u, "(A)", advance = "no") "* Squared matrix element correct: "
    write (u, "(L1)") nearly_equal (me_check_tot, &
       isolated_state%trace%get_matrix_element (1), rel_smallness = 0.00001_default)


    write (u, "(A)") "* Give a matrix to the isolated state"
    call create_test_model (var_str ("SM"), test_model)
    call isolated_state%setup_square_matrix (core, test_model, qn_mask, col_flow_index)
    call isolated_state%evaluate_matrix ()

    write (u, "(A)") "* Sub-matrixelements correct: "
    tmp1 = nearly_equal (me_check_1, &
       isolated_state%matrix%get_matrix_element (1), rel_smallness = 0.00001_default)
    tmp2 = nearly_equal (me_check_2, &
       isolated_state%matrix%get_matrix_element (2), rel_smallness = 0.00001_default)
    write (u, "(A,L1,A,L1)") "* 1: ", tmp1, ", 2: ", tmp2

    write (u, "(A)") "* Test output end: parton_states_1"
  end subroutine parton_states_1


end module parton_states_uti
