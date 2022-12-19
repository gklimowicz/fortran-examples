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

module eval_trees_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string

  use ifiles
  use lexers
  use lorentz
  use syntax_rules, only: syntax_write
  use pdg_arrays
  use subevents
  use variables
  use observables

  use eval_trees

  implicit none
  private

  public :: expressions_1
  public :: expressions_2
  public :: expressions_3
  public :: expressions_4

contains

  subroutine expressions_1 (u)
    integer, intent(in) :: u
    type(var_list_t), pointer :: var_list => null ()
    type(eval_node_t), pointer :: node => null ()
    type(prt_t), pointer :: prt => null ()
    type(string_t) :: var_name

    write (u, "(A)")  "* Test output: Expressions"
    write (u, "(A)")  "*   Purpose: test simple observable and node evaluation"
    write (u, "(A)")

    write (u, "(A)")  "* Setting a unary observable:"
    write (u, "(A)")

    allocate (var_list)
    allocate (prt)
    call var_list_set_observables_unary (var_list, prt)
    call var_list%write (u)

    write (u, "(A)")  "* Evaluating the observable node:"
    write (u, "(A)")

    var_name = "PDG"

    allocate (node)
    call node%test_obs (var_list, var_name)
    call node%write (u)

    write (u, "(A)")  "* Cleanup"
    write (u, "(A)")

    call node%final_rec ()
    deallocate (node)
    !!! Workaround for NAGFOR 6.2 
    ! call var_list%final ()
    deallocate (var_list)
    deallocate (prt)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: expressions_1"

  end subroutine expressions_1

  subroutine expressions_2 (u)
    integer, intent(in) :: u
    type(ifile_t) :: ifile
    type(stream_t) :: stream
    type(eval_tree_t) :: eval_tree
    type(string_t) :: expr_text
    type(var_list_t), pointer :: var_list => null ()

    write (u, "(A)")  "* Test output: Expressions"
    write (u, "(A)")  "*   Purpose: test parse routines"
    write (u, "(A)")

    call syntax_expr_init ()
    call syntax_write (syntax_expr, u)
    allocate (var_list)
    call var_list%append_real (var_str ("tolerance"), 0._default)
    call var_list%append_real (var_str ("x"), -5._default)
    call var_list%append_int  (var_str ("foo"), -27)
    call var_list%append_real (var_str ("mb"), 4._default)
    expr_text = &
         "let real twopi = 2 * pi in" // &
         "  twopi * sqrt (25.d0 - mb^2)" // &
         "  / (let int mb_or_0 = max (mb, 0) in" // &
         "       1 + (if -1 TeV <= x < mb_or_0 then abs(x) else x endif))"
    call ifile_append (ifile, expr_text)
    call stream_init (stream, ifile)
    call var_list%write (u)
    call eval_tree%init_stream (stream, var_list=var_list)
    call eval_tree%evaluate ()
    call eval_tree%write (u)

    write (u, "(A)")  "* Input string:"
    write (u, "(A,A)")  "     ", char (expr_text)
    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call stream_final (stream)
    call ifile_final (ifile)
    call eval_tree%final ()
    call var_list%final ()
    deallocate (var_list)
    call syntax_expr_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: expressions_2"

  end subroutine expressions_2

  subroutine expressions_3 (u)
    integer, intent(in) :: u
    type(subevt_t) :: subevt

    write (u, "(A)")  "* Test output: Expressions"
    write (u, "(A)")  "*   Purpose: test subevent expressions"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize subevent:"
    write (u, "(A)")

    call subevt_init (subevt)
    call subevt%reset (1)
    call subevt%set_incoming (1, 22, &
         vector4_moving (1.e3_default, 1.e3_default, 1), 0._default, [2])
    call subevt%write (u)
    call subevt%reset (4)
    call subevt%reset (3)
    call subevt%set_incoming (1, 21, &
         vector4_moving (1.e3_default, 1.e3_default, 3), 0._default, [1])
    call subevt_polarize (subevt, 1, -1)
    call subevt%set_outgoing (2, 1, &
         vector4_moving (0._default, 1.e3_default, 3), &
         -1.e6_default, [7])
    call subevt%set_composite (3, &
         vector4_moving (-1.e3_default, 0._default, 3), [2, 7])
    call subevt%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: expressions_3"

  end subroutine expressions_3

  subroutine expressions_4 (u)
    integer, intent(in) :: u
    type(subevt_t), target :: subevt
    type(string_t) :: expr_text
    type(ifile_t) :: ifile
    type(stream_t) :: stream
    type(eval_tree_t) :: eval_tree
    type(var_list_t), pointer :: var_list => null ()
    type(pdg_array_t) :: aval

    write (u, "(A)")  "* Test output: Expressions"
    write (u, "(A)")  "*   Purpose: test pdg array expressions"
    write (u, "(A)")

    write (u, "(A)")  "* Initialization:"
    write (u, "(A)")

    call syntax_pexpr_init ()
    call syntax_write (syntax_pexpr, u)
    allocate (var_list)
    call var_list%append_real (var_str ("tolerance"), 0._default)
    aval = 0
    call var_list%append_pdg_array (var_str ("particle"), aval)
    aval = [11,-11]
    call var_list%append_pdg_array (var_str ("lepton"), aval)
    aval = 22
    call var_list%append_pdg_array (var_str ("photon"), aval)
    aval = 1
    call var_list%append_pdg_array (var_str ("u"), aval)
    call subevt_init (subevt)
    call subevt%reset (6)
    call subevt%set_incoming (1, 1, &
         vector4_moving (1._default, 1._default, 1), 0._default)
    call subevt%set_incoming (2, -1, &
         vector4_moving (2._default, 2._default, 1), 0._default)
    call subevt%set_outgoing (3, 22, &
         vector4_moving (3._default, 3._default, 1), 0._default)
    call subevt%set_outgoing (4, 22, &
         vector4_moving (4._default, 4._default, 1), 0._default)
    call subevt%set_outgoing (5, 11, &
         vector4_moving (5._default, 5._default, 1), 0._default)
    call subevt%set_outgoing (6, -11, &
         vector4_moving (6._default, 6._default, 1), 0._default)
    write (u, "(A)")
    write (u, "(A)")  "* Expression:"
    expr_text = &
         "let alias quark = pdg(1):pdg(2):pdg(3) in" // &
         "  any E > 3 GeV " // &
         "    [sort by - Pt " // &
         "       [select if Index < 6 " // &
         "          [photon:pdg(-11):pdg(3):quark " // &
         "           & incoming particle]]]" // &
         "  and" // &
         "  eval Theta [extract index -1 [photon]] > 45 degree" // &
         "  and" // &
         "  count [incoming photon] * 3 > 0"
    write (u, "(A,A)")  "     ", char (expr_text)
    write (u, "(A)")

    write (u, "(A)")
    write (u, "(A)")  "* Extract the evaluation tree:"
    write (u, "(A)")

    call ifile_append (ifile, expr_text)
    call stream_init (stream, ifile)
    call eval_tree%init_stream (stream, var_list, subevt, V_LOG)
    call eval_tree%write (u)
    call eval_tree%evaluate ()

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate the tree:"
    write (u, "(A)")

    call eval_tree%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"
    write (u, "(A)")

    call stream_final (stream)
    call ifile_final (ifile)
    call eval_tree%final ()
    call var_list%final ()
    deallocate (var_list)
    call syntax_pexpr_final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: expressions_4"

  end subroutine expressions_4


end module eval_trees_uti
