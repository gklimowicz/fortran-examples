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

module parser

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lexers
  use syntax_rules

  implicit none
  private

  public :: parse_node_t
  public :: parse_node_p
  public :: parse_node_write_rec
  public :: parse_node_write
  public :: parse_node_final
  public :: parse_node_create_key
  public :: parse_node_create_value
  public :: parse_node_set_value
  public :: parse_node_create_branch
  public :: parse_node_append_sub
  public :: parse_node_freeze_branch
  public :: parse_node_replace_rule
  public :: parse_node_replace_last_sub
  public :: parse_node_get_rule_ptr
  public :: parse_node_get_n_sub
  public :: parse_node_get_sub_ptr
  public :: parse_node_get_next_ptr
  public :: parse_node_get_last_sub_ptr
  public :: parse_node_check
  public :: parse_node_mismatch
  public :: parse_node_get_logical
  public :: parse_node_get_integer
  public :: parse_node_get_real
  public :: parse_node_get_cmplx
  public :: parse_node_get_string
  public :: parse_node_get_key
  public :: parse_node_get_rule_key
  public :: parse_node_get_md5sum
  public :: parse_tree_t
  public :: parse_tree_init
  public :: parse_tree_final
  public :: parse_tree_write
  public :: parse_tree_bug
  public :: parse_tree_reduce
  public :: parse_tree_get_process_ptr

  type :: token_t
     private
     integer :: type = S_UNKNOWN
     logical, pointer :: lval => null ()
     integer, pointer :: ival => null ()
     real(default), pointer :: rval => null ()
     complex(default), pointer :: cval => null ()
     type(string_t), pointer :: sval => null ()
     type(string_t), pointer :: kval => null ()
     type(string_t), dimension(:), pointer :: quote => null ()
  end type token_t

  type :: parse_node_t
     private
     type(syntax_rule_t), pointer :: rule => null ()
     type(token_t) :: token
     integer :: n_sub = 0
     type(parse_node_t), pointer :: sub_first => null ()
     type(parse_node_t), pointer :: sub_last => null ()
     type(parse_node_t), pointer :: next => null ()
   contains
     procedure :: write => parse_node_write_rec
     procedure :: copy => parse_node_copy
     procedure :: append_sub => parse_node_append_sub
     procedure :: get_rule_ptr => parse_node_get_rule_ptr
     procedure :: get_n_sub => parse_node_get_n_sub
     procedure :: get_sub_ptr => parse_node_get_sub_ptr
     procedure :: get_next_ptr => parse_node_get_next_ptr
     procedure :: get_logical => parse_node_get_logical
     procedure :: get_integer => parse_node_get_integer
     procedure :: get_real => parse_node_get_real
     procedure :: get_cmplx => parse_node_get_cmplx
     procedure :: get_string => parse_node_get_string
     procedure :: get_key => parse_node_get_key
     procedure :: get_rule_key => parse_node_get_rule_key
  end type parse_node_t

  type :: parse_node_p
    type(parse_node_t), pointer :: ptr => null ()
  end type parse_node_p

  type :: parse_tree_t
     private
     type(parse_node_t), pointer :: root_node => null ()
   contains
     procedure :: parse => parse_tree_init
     procedure :: final => parse_tree_final
     procedure :: write => parse_tree_write
     procedure :: get_root_ptr => parse_tree_get_root_ptr
  end type parse_tree_t


  interface assignment(=)
     module procedure token_assign
     module procedure token_assign_integer
     module procedure token_assign_real
     module procedure token_assign_complex
     module procedure token_assign_logical
     module procedure token_assign_string
  end interface


  interface
    module subroutine token_assign (token, token_in)
      type(token_t), intent(out) :: token
      type(token_t), intent(in) :: token_in
    end subroutine token_assign
    module subroutine token_assign_integer (token, ival)
      type(token_t), intent(out) :: token
      integer, intent(in) :: ival
    end subroutine token_assign_integer
    module subroutine token_assign_real (token, rval)
      type(token_t), intent(out) :: token
      real(default), intent(in) :: rval
    end subroutine token_assign_real
    module subroutine token_assign_complex (token, cval)
      type(token_t), intent(out) :: token
      complex(default), intent(in) :: cval
    end subroutine token_assign_complex
    module subroutine token_assign_logical (token, lval)
      type(token_t), intent(out) :: token
      logical, intent(in) :: lval
    end subroutine token_assign_logical
    module subroutine token_assign_string (token, sval)
      type(token_t), intent(out) :: token
      type(string_t), intent(in) :: sval
    end subroutine token_assign_string
    recursive module subroutine parse_node_write_rec (node, unit, short, depth)
      class(parse_node_t), intent(in), target :: node
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: short
      integer, intent(in), optional :: depth
    end subroutine parse_node_write_rec
    module subroutine parse_node_write (node, unit, short)
      class(parse_node_t), intent(in) :: node
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: short
    end subroutine parse_node_write
    recursive module subroutine parse_node_final (node, recursive)
      type(parse_node_t), intent(inout) :: node
      logical, intent(in), optional :: recursive
    end subroutine parse_node_final
    module subroutine parse_node_create_key (node, rule)
      type(parse_node_t), intent(out) :: node
      type(syntax_rule_t), intent(in), target :: rule
    end subroutine parse_node_create_key
    module subroutine parse_node_create_value (node, rule, ival, rval, cval, sval, lval)
      type(parse_node_t), intent(out) :: node
      type(syntax_rule_t), intent(in), target :: rule
      integer, intent(in), optional :: ival
      real(default), intent(in), optional :: rval
      complex(default), intent(in), optional :: cval
      type(string_t), intent(in), optional :: sval
      logical, intent(in), optional :: lval
    end subroutine parse_node_create_value
    module subroutine parse_node_set_value (node, ival, rval, cval, sval, lval)
      type(parse_node_t), intent(inout) :: node
      integer, intent(in), optional :: ival
      real(default), intent(in), optional :: rval
      complex(default), intent(in), optional :: cval
      type(string_t), intent(in), optional :: sval
      logical, intent(in), optional :: lval
    end subroutine parse_node_set_value
    module subroutine parse_node_create_branch (node, rule)
      type(parse_node_t), pointer :: node
      type(syntax_rule_t), intent(in), target :: rule
    end subroutine parse_node_create_branch
    module subroutine parse_node_copy (node, copy)
      class(parse_node_t), intent(in) :: node
      type(parse_node_t), pointer, intent(out) :: copy
    end subroutine parse_node_copy
    module subroutine parse_node_append_sub (node, sub)
      class(parse_node_t), intent(inout) :: node
      type(parse_node_t), pointer :: sub
    end subroutine parse_node_append_sub
    module subroutine parse_node_freeze_branch (node)
      type(parse_node_t), pointer :: node
    end subroutine parse_node_freeze_branch
    module subroutine parse_node_replace_rule (node, rule)
      type(parse_node_t), pointer :: node
      type(syntax_rule_t), intent(in), target :: rule
    end subroutine parse_node_replace_rule
    module subroutine parse_node_replace_last_sub (node, pn_target)
      type(parse_node_t), intent(inout), target :: node
      type(parse_node_t), intent(in), target :: pn_target
    end subroutine parse_node_replace_last_sub
    module function parse_node_get_rule_ptr (node) result (rule)
      class(parse_node_t), intent(in) :: node
      type(syntax_rule_t), pointer :: rule
    end function parse_node_get_rule_ptr
    module function parse_node_get_last_sub_ptr (node, tag, required) result (sub)
      type(parse_node_t), pointer :: sub
      type(parse_node_t), intent(in), target :: node
      character(*), intent(in), optional :: tag
      logical, intent(in), optional :: required
    end function parse_node_get_last_sub_ptr
    module function parse_node_get_n_sub (node) result (n)
      class(parse_node_t), intent(in) :: node
      integer :: n
    end function parse_node_get_n_sub
    module function parse_node_get_sub_ptr (node, n, tag, required) result (sub)
      class(parse_node_t), intent(in), target :: node
      type(parse_node_t), pointer :: sub
      integer, intent(in), optional :: n
      character(*), intent(in), optional :: tag
      logical, intent(in), optional :: required
    end function parse_node_get_sub_ptr
    module function parse_node_get_next_ptr (sub, n, tag, required) result (next)
      class(parse_node_t), intent(in), target :: sub
      type(parse_node_t), pointer :: next
      integer, intent(in), optional :: n
      character(*), intent(in), optional :: tag
      logical, intent(in), optional :: required
    end function parse_node_get_next_ptr
    module subroutine parse_node_check (node, tag, required)
      type(parse_node_t), pointer :: node
      character(*), intent(in), optional :: tag
      logical, intent(in), optional :: required
    end subroutine parse_node_check
    module subroutine parse_node_mismatch (string, parse_node)
      character(*), intent(in) :: string
      type(parse_node_t), intent(in) :: parse_node
    end subroutine parse_node_mismatch
    module function parse_node_get_logical (node) result (lval)
      class(parse_node_t), intent(in), target :: node
      logical :: lval
    end function parse_node_get_logical
    module function parse_node_get_integer (node) result (ival)
      class(parse_node_t), intent(in), target :: node
      integer :: ival
    end function parse_node_get_integer
    module function parse_node_get_real (node) result (rval)
      class(parse_node_t), intent(in), target :: node
      real(default) :: rval
    end function parse_node_get_real
    module function parse_node_get_cmplx (node) result (cval)
      class(parse_node_t), intent(in), target :: node
      complex(default) :: cval
    end function parse_node_get_cmplx
    module function parse_node_get_string (node) result (sval)
      class(parse_node_t), intent(in), target :: node
      type(string_t) :: sval
    end function parse_node_get_string
    module function parse_node_get_key (node) result (kval)
      class(parse_node_t), intent(in), target :: node
      type(string_t) :: kval
    end function parse_node_get_key
    module function parse_node_get_rule_key (node) result (kval)
      class(parse_node_t), intent(in), target :: node
      type(string_t) :: kval
    end function parse_node_get_rule_key
    module function parse_node_get_token_ptr (node) result (token)
      type(token_t), pointer :: token
      type(parse_node_t), intent(in), target :: node
    end function parse_node_get_token_ptr
    module function parse_node_get_md5sum (pn) result (md5sum_pn)
      character(32) :: md5sum_pn
      type(parse_node_t), intent(in) :: pn
    end function parse_node_get_md5sum
    module subroutine parse_tree_init &
         (parse_tree, syntax, lexer, key, check_eof)
      class(parse_tree_t), intent(inout) :: parse_tree
      type(lexer_t), intent(inout) :: lexer
      type(syntax_t), intent(in), target :: syntax
      type(string_t), intent(in), optional :: key
      logical, intent(in), optional :: check_eof
    end subroutine parse_tree_init
    module subroutine parse_tree_final (parse_tree)
      class(parse_tree_t), intent(inout) :: parse_tree
    end subroutine parse_tree_final
    module subroutine parse_tree_write (parse_tree, unit, verbose)
      class(parse_tree_t), intent(in) :: parse_tree
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine parse_tree_write
    module subroutine parse_tree_bug (node, keys)
      type(parse_node_t), intent(in) :: node
      character(*), intent(in) :: keys
    end subroutine parse_tree_bug
    module function parse_tree_get_root_ptr (parse_tree) result (node)
      class(parse_tree_t), intent(in) :: parse_tree
      type(parse_node_t), pointer :: node
    end function parse_tree_get_root_ptr
    module subroutine parse_tree_reduce (parse_tree, rule_key)
      type(parse_tree_t), intent(inout) :: parse_tree
      type(string_t), dimension(:), intent(in) :: rule_key
    end subroutine parse_tree_reduce
    module function parse_tree_get_process_ptr (parse_tree, process) result (node)
      type(parse_node_t), pointer :: node
      type(parse_tree_t), intent(in), target :: parse_tree
      type(string_t), intent(in) :: process
    end function parse_tree_get_process_ptr
  end interface

end module parser
