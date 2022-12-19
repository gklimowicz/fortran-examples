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

module evaluators_uti

  use kinds, only: default
  use lorentz
  use flavors
  use colors
  use helicities
  use quantum_numbers
  use interactions
  use model_data

  use evaluators

  implicit none
  private

  public :: evaluator_1
  public :: evaluator_2
  public :: evaluator_3
  public :: evaluator_4

contains

  subroutine evaluator_1 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(interaction_t), target :: int_qqtt, int_tbw, int1, int2
    type(flavor_t), dimension(:), allocatable :: flv
    type(color_t), dimension(:), allocatable :: col
    type(helicity_t), dimension(:), allocatable :: hel
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    integer :: f, c, h1, h2, h3
    type(vector4_t), dimension(4) :: p
    type(vector4_t), dimension(2) :: q
    type(quantum_numbers_mask_t) :: qn_mask_conn
    type(quantum_numbers_mask_t), dimension(:), allocatable :: qn_mask2
    type(evaluator_t), target :: eval, eval2, eval3

    call model%init_sm_test ()

    write (u, "(A)")   "*** Evaluator for matrix product"
    write (u, "(A)")   "***   Construct interaction for qq -> tt"
    write (u, "(A)")
    call int_qqtt%basic_init (2, 0, 2, set_relations=.true.)
    allocate (flv (4), col (4), hel (4), qn (4))
    allocate (qn_mask2 (4))
    do c = 1, 2
       select case (c)
       case (1)
          call col%init_col_acl ([1, 0, 1, 0], [0, 2, 0, 2])
       case (2)
          call col%init_col_acl ([1, 0, 2, 0], [0, 1, 0, 2])
       end select
       do f = 1, 2
          call flv%init ([f, -f, 6, -6], model)
          do h1 = -1, 1, 2
             call hel(3)%init (h1)
             do h2 = -1, 1, 2
                call hel(4)%init (h2)
                call qn%init (flv, col, hel)
                call int_qqtt%add_state (qn)
             end do
          end do
       end do
    end do
    call int_qqtt%freeze ()
    deallocate (flv, col, hel, qn)
    write (u, "(A)")  "***   Construct interaction for t -> bW"
    call int_tbw%basic_init (1, 0, 2, set_relations=.true.)
    allocate (flv (3), col (3), hel (3), qn (3))
    call flv%init ([6, 5, 24], model)
    call col%init_col_acl ([1, 1, 0], [0, 0, 0])
    do h1 = -1, 1, 2
       call hel(1)%init (h1)
       do h2 = -1, 1, 2
          call hel(2)%init (h2)
          do h3 = -1, 1
             call hel(3)%init (h3)
             call qn%init (flv, col, hel)
             call int_tbw%add_state (qn)
          end do
       end do
    end do
    call int_tbw%freeze ()
    deallocate (flv, col, hel, qn)
    write (u, "(A)")  "***   Link interactions"
    call int_tbw%set_source_link (1, int_qqtt, 3)
    qn_mask_conn = quantum_numbers_mask (.false.,.false.,.true.)
    write (u, "(A)")
    write (u, "(A)")  "***   Show input"
    call int_qqtt%basic_write (unit = u)
    write (u, "(A)")
    call int_tbw%basic_write (unit = u)
    write (u, "(A)")
    write (u, "(A)")  "***   Evaluate product"
    call eval%init_product (int_qqtt, int_tbw, qn_mask_conn)
    call eval%write (unit = u)

     call int1%basic_init (2, 0, 2, set_relations=.true.)
     call int2%basic_init (1, 0, 2, set_relations=.true.)
     p(1) = vector4_moving (1000._default, 1000._default, 3)
     p(2) = vector4_moving (200._default, 200._default, 2)
     p(3) = vector4_moving (100._default, 200._default, 1)
     p(4) = p(1) - p(2) - p(3)
     call int1%set_momenta (p)
     q(1) = vector4_moving (50._default,-50._default, 3)
     q(2) = p(2) + p(4) - q(1)
     call int2%set_momenta (q, outgoing=.true.)
     call int1%set_matrix_element ([(2._default,0._default), &
          (4._default,1._default), (-3._default,0._default)])
     call int2%set_matrix_element ([(-3._default,0._default), &
          (0._default,1._default), (1._default,2._default)])
     call eval%receive_momenta ()
     call eval%evaluate ()
     call int1%basic_write (unit = u)
     write (u, "(A)")
     call int2%basic_write (unit = u)
     write (u, "(A)")
     call eval%write (unit = u)
     write (u, "(A)")
     call int1%final ()
     call int2%final ()
     call eval%final ()

     write (u, "(A)")
     write (u, "(A)")   "*** Evaluator for matrix square"
     allocate (flv(4), col(4), qn(4))
     call int1%basic_init (2, 0, 2, set_relations=.true.)
     call flv%init ([1, -1, 21, 21], model)
     call col(1)%init ([1])
     call col(2)%init ([-2])
     call col(3)%init ([2, -3])
     call col(4)%init ([3, -1])
     call qn%init (flv, col)
     call int1%add_state (qn)
     call col(3)%init ([3, -1])
     call col(4)%init ([2, -3])
     call qn%init (flv, col)
     call int1%add_state (qn)
     call col(3)%init ([2, -1])
     call col(4)%init (.true.)
     call qn%init (flv, col)
     call int1%add_state (qn)
     call int1%freeze ()
     ! [qn_mask2 not set since default is false]
     call eval%init_square (int1, qn_mask2, nc=3)
     call eval2%init_square_nondiag (int1, qn_mask2)
     qn_mask2 = quantum_numbers_mask (.false., .true., .true.)
     call eval3%init_square_diag (eval, qn_mask2)
     call int1%set_matrix_element &
          ([(2._default,0._default), &
          (4._default,1._default), (-3._default,0._default)])
     call int1%set_momenta (p)
     call int1%basic_write (unit = u)
     write (u, "(A)")
     call eval%receive_momenta ()
     call eval%evaluate ()
     call eval%write (unit = u)
     write (u, "(A)")
     call eval2%receive_momenta ()
     call eval2%evaluate ()
     call eval2%write (unit = u)
     write (u, "(A)")
     call eval3%receive_momenta ()
     call eval3%evaluate ()
     call eval3%write (unit = u)
     call int1%final ()
     call eval%final ()
     call eval2%final ()
     call eval3%final ()

     call model%final ()
  end subroutine evaluator_1

  subroutine evaluator_2 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(interaction_t), target :: int
    integer :: h1, h2, h3, h4
    type(helicity_t), dimension(4) :: hel
    type(color_t), dimension(4) :: col
    type(flavor_t), dimension(4) :: flv
    type(quantum_numbers_t), dimension(4) :: qn
    type(vector4_t), dimension(4) :: p
    type(evaluator_t) :: eval
    integer :: i

    call model%init_sm_test ()

    write (u, "(A)") "*** Creating interaction for e+ e- -> W+ W-"
    write (u, "(A)")

    call flv%init ([11, -11, 24, -24], model)
    do i = 1, 4
       call col(i)%init ()
    end do
    call int%basic_init (2, 0, 2, set_relations=.true.)
    do h1 = -1, 1, 2
       call hel(1)%init (h1)
       do h2 = -1, 1, 2
          call hel(2)%init (h2)
          do h3 = -1, 1
             call hel(3)%init (h3)
             do h4 = -1, 1
                call hel(4)%init (h4)
                call qn%init (flv, col, hel)
                call int%add_state (qn)
             end do
          end do
       end do
    end do
    call int%freeze ()
    call int%set_matrix_element &
       ([(cmplx (i, kind=default), i = 1, 36)])
    p(1) = vector4_moving (1000._default, 1000._default, 3)
    p(2) = vector4_moving (1000._default, -1000._default, 3)
    p(3) = vector4_moving (1000._default, &
       sqrt (1E6_default - 80._default**2), 3)
    p(4) = p(1) + p(2) - p(3)
    call int%set_momenta (p)
    write (u, "(A)") "*** Setting up evaluator"
    write (u, "(A)")

    call eval%init_identity (int)
    write (u, "(A)") "*** Transferring momenta and evaluating"
    write (u, "(A)")

    call eval%receive_momenta ()
    call eval%evaluate ()
    write (u, "(A)")  "*******************************************************"
    write (u, "(A)")  "   Interaction dump"
    write (u, "(A)")  "*******************************************************"
    call int%basic_write (unit = u)
    write (u, "(A)")
    write (u, "(A)")  "*******************************************************"
    write (u, "(A)")  "   Evaluator dump"
    write (u, "(A)")  "*******************************************************"
    call eval%write (unit = u)
    write (u, "(A)")
    write (u, "(A)")   "*** cleaning up"
    call int%final ()
    call eval%final ()

    call model%final ()
  end subroutine evaluator_2

  subroutine evaluator_3 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(interaction_t), target :: int
    integer :: h1, h2, h3, h4
    type(helicity_t), dimension(4) :: hel
    type(color_t), dimension(4) :: col
    type(flavor_t), dimension(4) :: flv1, flv2
    type(quantum_numbers_t), dimension(4) :: qn
    type(vector4_t), dimension(4) :: p
    type(evaluator_t) :: eval1, eval2, eval3
    type(quantum_numbers_mask_t), dimension(4) :: qn_mask
    integer :: i

    call model%init_sm_test ()

    write (u, "(A)")  "*** Creating interaction for e+/mu+ e-/mu- -> W+ W-"
    call flv1%init ([11, -11, 24, -24], model)
    call flv2%init ([13, -13, 24, -24], model)
    do i = 1, 4
       call col (i)%init ()
    end do
    call int%basic_init (2, 0, 2, set_relations=.true.)
    do h1 = -1, 1, 2
       call hel(1)%init (h1)
       do h2 = -1, 1, 2
          call hel(2)%init (h2)
          do h3 = -1, 1
             call hel(3)%init (h3)
             do h4 = -1, 1
                call hel(4)%init (h4)
                call qn%init (flv1, col, hel)
                call int%add_state (qn)
                call qn%init (flv2, col, hel)
                call int%add_state (qn)
             end do
          end do
       end do
    end do
    call int%freeze ()
    call int%set_matrix_element &
       ([(cmplx (1, kind=default), i = 1, 72)])
    p(1) = vector4_moving (1000._default, 1000._default, 3)
    p(2) = vector4_moving (1000._default, -1000._default, 3)
    p(3) = vector4_moving (1000._default, &
       sqrt (1E6_default - 80._default**2), 3)
    p(4) = p(1) + p(2) - p(3)
    call int%set_momenta (p)
    write (u, "(A)")  "*** Setting up evaluators"
    call qn_mask%init (.false., .true., .true.)
    call eval1%init_qn_sum (int, qn_mask)
    call qn_mask%init (.true., .true., .true.)
    call eval2%init_qn_sum (int, qn_mask)
    call qn_mask%init (.false., .true., .false.)
    call eval3%init_qn_sum (int, qn_mask, &
      [.false., .false., .false., .true.])
    write (u, "(A)")  "*** Transferring momenta and evaluating"
    call eval1%receive_momenta ()
    call eval1%evaluate ()
    call eval2%receive_momenta ()
    call eval2%evaluate ()
    call eval3%receive_momenta ()
    call eval3%evaluate ()
    write (u, "(A)")  "*******************************************************"
    write (u, "(A)")  "   Interaction dump"
    write (u, "(A)")  "*******************************************************"
    call int%basic_write (unit = u)
    write (u, "(A)")
    write (u, "(A)")  "*******************************************************"
    write (u, "(A)")  "   Evaluator dump --- spin sum"
    write (u, "(A)")  "*******************************************************"
    call eval1%write (unit = u)
    call eval1%basic_write (unit = u)
    write (u, "(A)")  "*******************************************************"
    write (u, "(A)")  "   Evaluator dump --- spin / flavor sum"
    write (u, "(A)")  "*******************************************************"
    call eval2%write (unit = u)
    call eval2%basic_write (unit = u)
    write (u, "(A)")  "*******************************************************"
    write (u, "(A)")  "   Evaluator dump --- flavor sum, drop last W"
    write (u, "(A)")  "*******************************************************"
    call eval3%write (unit = u)
    call eval3%basic_write (unit = u)
    write (u, "(A)")
    write (u, "(A)")  "*** cleaning up"
    call int%final ()
    call eval1%final ()
    call eval2%final ()
    call eval3%final ()

    call model%final ()
  end subroutine evaluator_3

  subroutine evaluator_4 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(interaction_t), target :: int1, int2
    integer :: h1, h2, h3
    type(helicity_t), dimension(3) :: hel
    type(color_t), dimension(3) :: col
    type(flavor_t), dimension(2) :: flv1, flv2
    type(flavor_t), dimension(3) :: flv3, flv4
    type(quantum_numbers_t), dimension(3) :: qn
    type(evaluator_t) :: eval1, eval2, eval3, eval4
    type(quantum_numbers_mask_t) :: qn_mask
    type(flavor_t) :: flv_filter
    type(helicity_t) :: hel_filter
    type(color_t) :: col_filter
    type(quantum_numbers_t) :: qn_filter
    integer :: i

    write (u, "(A)")  "* Test output: evaluator_4"
    write (u, "(A)")  "*   Purpose: test evaluator products &
         &with mask and filter"
    write (u, "(A)")

    call model%init_sm_test ()

    write (u, "(A)")  "* Creating interaction for e- -> W+/Z"
    write (u, "(A)")

    call flv1%init ([11, 24], model)
    call flv2%init ([11, 23], model)
    do i = 1, 3
       call col(i)%init ()
    end do
    call int1%basic_init (1, 0, 1, set_relations=.true.)
    do h1 = -1, 1, 2
       call hel(1)%init (h1)
       do h2 = -1, 1
          call hel(2)%init (h2)
          call qn(:2)%init (flv1, col(:2), hel(:2))
          call int1%add_state (qn(:2))
          call qn(:2)%init (flv2, col(:2), hel(:2))
          call int1%add_state (qn(:2))
       end do
    end do
    call int1%freeze ()
    call int1%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Creating interaction for W+/Z -> u ubar/dbar"
    write (u, "(A)")

    call flv3%init ([24, 2, -1], model)
    call flv4%init ([23, 2, -2], model)

    call int2%basic_init (1, 0, 2, set_relations=.true.)
    do h1 = -1, 1
       call hel(1)%init (h1)
       do h2 = -1, 1, 2
          call hel(2)%init (h2)
          do h3 = -1, 1, 2
             call hel(3)%init (h3)
             call qn(:3)%init (flv3, col(:3), hel(:3))
             call int2%add_state (qn(:3))
             call qn(:3)%init (flv4, col(:3), hel(:3))
             call int2%add_state (qn(:3))
          end do
       end do
    end do
    call int2%freeze ()

    call int2%set_source_link (1, int1, 2)
    call int2%basic_write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Product evaluator"
    write (u, "(A)")

    call qn_mask%init (.false., .false., .false.)
    call eval1%init_product (int1, int2, qn_mask_conn = qn_mask)
    call eval1%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Product evaluator with helicity mask"
    write (u, "(A)")

    call qn_mask%init (.false., .false., .true.)
    call eval2%init_product (int1, int2, qn_mask_conn = qn_mask)
    call eval2%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Product with flavor filter and helicity mask"
    write (u, "(A)")

    call qn_mask%init (.false., .false., .true.)
    call flv_filter%init (24, model)
    call hel_filter%init ()
    call col_filter%init ()
    call qn_filter%init (flv_filter, col_filter, hel_filter)
    call eval3%init_product (int1, int2, &
         qn_mask_conn = qn_mask, qn_filter_conn = qn_filter)
    call eval3%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Product with helicity filter and mask"
    write (u, "(A)")

    call qn_mask%init (.false., .false., .true.)
    call flv_filter%init ()
    call hel_filter%init (0)
    call col_filter%init ()
    call qn_filter%init (flv_filter, col_filter, hel_filter)
    call eval4%init_product (int1, int2, &
         qn_mask_conn = qn_mask, qn_filter_conn = qn_filter)
    call eval4%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call eval1%final ()
    call eval2%final ()
    call eval3%final ()
    call eval4%final ()

    call int1%final ()
    call int2%final ()

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: evaluator_4"

  end subroutine evaluator_4


end module evaluators_uti
