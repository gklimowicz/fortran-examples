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

module eval_trees

  use, intrinsic :: iso_c_binding !NODEP!
  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use ifiles
  use lexers
  use syntax_rules
  use parser
  use pdg_arrays
  use subevents
  use var_base
  use expr_base
  use variables

  implicit none
  private

  public :: eval_node_t
  public :: syntax_expr
  public :: syntax_pexpr
  public :: syntax_expr_init
  public :: syntax_pexpr_init
  public :: syntax_expr_final
  public :: syntax_pexpr_final
  public :: syntax_pexpr_write
  public :: define_expr_syntax
  public :: parse_tree_init_expr
  public :: parse_tree_init_lexpr
  public :: parse_tree_init_pexpr
  public :: parse_tree_init_cexpr
  public :: parse_tree_init_sexpr
  public :: eval_tree_t
  public :: eval_log
  public :: eval_int
  public :: eval_real
  public :: eval_cmplx
  public :: eval_subevt
  public :: eval_pdg_array
  public :: eval_string
  public :: eval_numeric
  public :: eval_tree_factory_t

  integer, parameter :: EN_UNKNOWN = 0, EN_UNARY = 1, EN_BINARY = 2
  integer, parameter :: EN_CONSTANT = 3, EN_VARIABLE = 4
  integer, parameter :: EN_CONDITIONAL = 5, EN_BLOCK = 6
  integer, parameter :: EN_RECORD_CMD = 7
  integer, parameter :: EN_OBS1_INT = 11, EN_OBS2_INT = 12
  integer, parameter :: EN_OBSEV_INT = 13
  integer, parameter :: EN_OBS1_REAL = 21, EN_OBS2_REAL = 22
  integer, parameter :: EN_OBSEV_REAL = 23
  integer, parameter :: EN_PRT_FUN_UNARY = 101, EN_PRT_FUN_BINARY = 102
  integer, parameter :: EN_EVAL_FUN_UNARY = 111, EN_EVAL_FUN_BINARY = 112
  integer, parameter :: EN_LOG_FUN_UNARY = 121, EN_LOG_FUN_BINARY = 122
  integer, parameter :: EN_INT_FUN_UNARY = 131, EN_INT_FUN_BINARY = 132
  integer, parameter :: EN_REAL_FUN_UNARY = 141, EN_REAL_FUN_BINARY = 142
  integer, parameter :: EN_REAL_FUN_CUM = 151
  integer, parameter :: EN_FORMAT_STR = 161

  type :: eval_node_t
     private
     type(string_t) :: tag
     integer :: type = EN_UNKNOWN
     integer :: result_type = V_NONE
     type(var_list_t), pointer :: var_list => null ()
     type(string_t) :: var_name
     logical, pointer :: value_is_known => null ()
     logical,           pointer :: lval => null ()
     integer,           pointer :: ival => null ()
     real(default),     pointer :: rval => null ()
     complex(default),  pointer :: cval => null ()
     type(subevt_t),  pointer :: pval => null ()
     type(pdg_array_t), pointer :: aval => null ()
     type(string_t),    pointer :: sval => null ()
     type(eval_node_t), pointer :: arg0 => null ()
     type(eval_node_t), pointer :: arg1 => null ()
     type(eval_node_t), pointer :: arg2 => null ()
     type(eval_node_t), pointer :: arg3 => null ()
     type(eval_node_t), pointer :: arg4 => null ()
     procedure(obs_unary_int),   nopass, pointer :: obs1_int  => null ()
     procedure(obs_unary_real),  nopass, pointer :: obs1_real => null ()
     procedure(obs_binary_int),  nopass, pointer :: obs2_int  => null ()
     procedure(obs_binary_real), nopass, pointer :: obs2_real => null ()
     procedure(obs_sev_int), nopass, pointer :: obsev_int => null ()
     procedure(obs_sev_real), nopass, pointer :: obsev_real => null ()
     integer, pointer :: prt_type => null ()
     integer, pointer :: index => null ()
     real(default), pointer :: tolerance => null ()
     integer, pointer :: jet_algorithm => null ()
     real(default), pointer :: jet_r => null ()
     real(default), pointer :: jet_p => null ()
     real(default), pointer :: jet_ycut => null ()
     real(default), pointer :: jet_dcut => null ()
     real(default), pointer :: photon_iso_eps => null ()
     real(default), pointer :: photon_iso_n => null ()
     real(default), pointer :: photon_iso_r0 => null ()
     real(default), pointer :: photon_rec_r0 => null ()
     type(prt_t), pointer :: prt1 => null ()
     type(prt_t), pointer :: prt2 => null ()
     procedure(unary_log),  nopass, pointer :: op1_log  => null ()
     procedure(unary_int),  nopass, pointer :: op1_int  => null ()
     procedure(unary_real), nopass, pointer :: op1_real => null ()
     procedure(unary_cmplx), nopass, pointer :: op1_cmplx => null ()
     procedure(unary_pdg),  nopass, pointer :: op1_pdg  => null ()
     procedure(unary_sev),  nopass, pointer :: op1_sev  => null ()
     procedure(unary_str),  nopass, pointer :: op1_str  => null ()
     procedure(unary_cut),  nopass, pointer :: op1_cut  => null ()
     procedure(unary_evi),  nopass, pointer :: op1_evi  => null ()
     procedure(unary_evr),  nopass, pointer :: op1_evr  => null ()
     procedure(binary_log),  nopass, pointer :: op2_log  => null ()
     procedure(binary_int),  nopass, pointer :: op2_int  => null ()
     procedure(binary_real), nopass, pointer :: op2_real => null ()
     procedure(binary_cmplx), nopass, pointer :: op2_cmplx => null ()
     procedure(binary_pdg),  nopass, pointer :: op2_pdg  => null ()
     procedure(binary_sev),  nopass, pointer :: op2_sev  => null ()
     procedure(binary_str),  nopass, pointer :: op2_str  => null ()
     procedure(binary_cut),  nopass, pointer :: op2_cut  => null ()
     procedure(binary_evi),  nopass, pointer :: op2_evi  => null ()
     procedure(binary_evr),  nopass, pointer :: op2_evr  => null ()
     procedure(cum_evi), nopass, pointer :: opcum_evi => null ()
     procedure(cum_evr), nopass, pointer :: opcum_evr => null ()
   contains
     procedure :: final_rec => eval_node_final_rec
     procedure :: write => eval_node_write
     procedure :: test_obs => eval_node_test_obs
  end type eval_node_t

  type, extends (expr_t) :: eval_tree_t
     private
     type(parse_node_t), pointer :: pn => null ()
     type(var_list_t) :: var_list
     type(eval_node_t), pointer :: root => null ()
   contains
     procedure :: init_stream => eval_tree_init_stream
     procedure :: init_expr  => eval_tree_init_expr
     procedure :: init_lexpr => eval_tree_init_lexpr
     procedure :: init_pexpr => eval_tree_init_pexpr
     procedure :: init_cexpr => eval_tree_init_cexpr
     procedure :: init_sexpr => eval_tree_init_sexpr
     procedure :: setup_expr  => eval_tree_setup_expr
     procedure :: setup_lexpr => eval_tree_setup_lexpr
     procedure :: setup_pexpr => eval_tree_setup_pexpr
     procedure :: setup_cexpr => eval_tree_setup_cexpr
     procedure :: setup_sexpr => eval_tree_setup_sexpr
     procedure :: init_numeric_value => eval_tree_init_numeric_value
     procedure :: final => eval_tree_final
     procedure :: evaluate => eval_tree_evaluate
     procedure :: is_known => eval_tree_result_is_known
     procedure :: get_log => eval_tree_get_log
     procedure :: get_int => eval_tree_get_int
     procedure :: get_real => eval_tree_get_real
     procedure :: get_cmplx => eval_tree_get_cmplx
     procedure :: get_pdg_array => eval_tree_get_pdg_array
     procedure :: get_subevt => eval_tree_get_subevt
     procedure :: get_string => eval_tree_get_string
     procedure :: write => eval_tree_write
  end type eval_tree_t

  type, extends (expr_factory_t) :: eval_tree_factory_t
     private
     type(parse_node_t), pointer :: pn => null ()
   contains
    procedure :: write => eval_tree_factory_write
    procedure :: init => eval_tree_factory_init
    procedure :: build => eval_tree_factory_build
  end type eval_tree_factory_t


  abstract interface
     logical function unary_log (arg)
       import eval_node_t
       type(eval_node_t), intent(in) :: arg
     end function unary_log
  end interface
  abstract interface
     integer function unary_int (arg)
       import eval_node_t
       type(eval_node_t), intent(in) :: arg
     end function unary_int
  end interface
  abstract interface
     real(default) function unary_real (arg)
       import default
       import eval_node_t
       type(eval_node_t), intent(in) :: arg
     end function unary_real
  end interface
  abstract interface
     complex(default) function unary_cmplx (arg)
       import default
       import eval_node_t
       type(eval_node_t), intent(in) :: arg
     end function unary_cmplx
  end interface
  abstract interface
     subroutine unary_pdg (pdg_array, arg)
       import pdg_array_t
       import eval_node_t
       type(pdg_array_t), intent(out) :: pdg_array
       type(eval_node_t), intent(in) :: arg
     end subroutine unary_pdg
  end interface
  abstract interface
     subroutine unary_sev (subevt, arg, arg0)
       import subevt_t
       import eval_node_t
       type(subevt_t), intent(inout) :: subevt
       type(eval_node_t), intent(in) :: arg
       type(eval_node_t), intent(inout), optional :: arg0
     end subroutine unary_sev
  end interface
  abstract interface
     subroutine unary_str (string, arg)
       import string_t
       import eval_node_t
       type(string_t), intent(out) :: string
       type(eval_node_t), intent(in) :: arg
     end subroutine unary_str
  end interface
  abstract interface
     logical function unary_cut (arg1, arg0)
       import eval_node_t
       type(eval_node_t), intent(in) :: arg1
       type(eval_node_t), intent(inout) :: arg0
     end function unary_cut
  end interface
  abstract interface
     subroutine unary_evi (ival, arg1, arg0)
       import eval_node_t
       integer, intent(out) :: ival
       type(eval_node_t), intent(in) :: arg1
       type(eval_node_t), intent(inout), optional :: arg0
     end subroutine unary_evi
  end interface
  abstract interface
     subroutine unary_evr (rval, arg1, arg0)
       import eval_node_t, default
       real(default), intent(out) :: rval
       type(eval_node_t), intent(in) :: arg1
       type(eval_node_t), intent(inout), optional :: arg0
     end subroutine unary_evr
  end interface
  abstract interface
     logical function binary_log (arg1, arg2)
       import eval_node_t
       type(eval_node_t), intent(in) :: arg1, arg2
     end function binary_log
  end interface
  abstract interface
     integer function binary_int (arg1, arg2)
       import eval_node_t
       type(eval_node_t), intent(in) :: arg1, arg2
     end function binary_int
  end interface
  abstract interface
     real(default) function binary_real (arg1, arg2)
       import default
       import eval_node_t
       type(eval_node_t), intent(in) :: arg1, arg2
     end function binary_real
  end interface
  abstract interface
     complex(default) function binary_cmplx (arg1, arg2)
       import default
       import eval_node_t
       type(eval_node_t), intent(in) :: arg1, arg2
     end function binary_cmplx
  end interface
  abstract interface
     subroutine binary_pdg (pdg_array, arg1, arg2)
       import pdg_array_t
       import eval_node_t
       type(pdg_array_t), intent(out) :: pdg_array
       type(eval_node_t), intent(in) :: arg1, arg2
     end subroutine binary_pdg
  end interface
  abstract interface
     subroutine binary_sev (subevt, arg1, arg2, arg0)
       import subevt_t
       import eval_node_t
       type(subevt_t), intent(inout) :: subevt
       type(eval_node_t), intent(in) :: arg1, arg2
       type(eval_node_t), intent(inout), optional :: arg0
     end subroutine binary_sev
  end interface
  abstract interface
     subroutine binary_str (string, arg1, arg2)
       import string_t
       import eval_node_t
       type(string_t), intent(out) :: string
       type(eval_node_t), intent(in) :: arg1, arg2
     end subroutine binary_str
  end interface
  abstract interface
     logical function binary_cut (arg1, arg2, arg0)
       import eval_node_t
       type(eval_node_t), intent(in) :: arg1, arg2
       type(eval_node_t), intent(inout) :: arg0
     end function binary_cut
  end interface
  abstract interface
     subroutine binary_evi (ival, arg1, arg2, arg0)
       import eval_node_t
       integer, intent(out) :: ival
       type(eval_node_t), intent(in) :: arg1, arg2
       type(eval_node_t), intent(inout), optional :: arg0
     end subroutine binary_evi
  end interface
  abstract interface
     subroutine binary_evr (rval, arg1, arg2, arg0)
       import eval_node_t, default
       real(default), intent(out) :: rval
       type(eval_node_t), intent(in) :: arg1, arg2
       type(eval_node_t), intent(inout), optional :: arg0
     end subroutine binary_evr
  end interface
  abstract interface
     integer function cum_evi (arg1, arg0)
       import eval_node_t
       type(eval_node_t), intent(in) :: arg1
       type(eval_node_t), intent(inout) :: arg0
     end function cum_evi
  end interface
  abstract interface
     real(default) function cum_evr (arg1, arg0)
       import eval_node_t, default
       type(eval_node_t), intent(in) :: arg1
       type(eval_node_t), intent(inout) :: arg0
     end function cum_evr
  end interface


  type(syntax_t), target, save :: syntax_expr
  type(syntax_t), target, save :: syntax_pexpr


  interface
    recursive module subroutine eval_node_final_rec (node)
      class(eval_node_t), intent(inout) :: node
    end subroutine eval_node_final_rec
    module subroutine eval_node_write (node, unit, indent)
      class(eval_node_t), intent(in) :: node
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: indent
    end subroutine eval_node_write
    module subroutine eval_node_test_obs (node, var_list, var_name)
      class(eval_node_t), intent(inout) :: node
      type(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: var_name
    end subroutine eval_node_test_obs
    module subroutine syntax_expr_init ()
    end subroutine syntax_expr_init
    module subroutine syntax_pexpr_init ()
    end subroutine syntax_pexpr_init
    module subroutine syntax_expr_final ()
    end subroutine syntax_expr_final
    module subroutine syntax_pexpr_final ()
    end subroutine syntax_pexpr_final
    module subroutine syntax_pexpr_write (unit)
      integer, intent(in), optional :: unit
    end subroutine syntax_pexpr_write
    module subroutine define_expr_syntax (ifile, particles, analysis)
      type(ifile_t), intent(inout) :: ifile
      logical, intent(in) :: particles, analysis
    end subroutine define_expr_syntax
    module subroutine parse_tree_init_expr (parse_tree, stream, particles)
      type(parse_tree_t), intent(out) :: parse_tree
      type(stream_t), intent(inout), target :: stream
      logical, intent(in) :: particles
    end subroutine parse_tree_init_expr
    module subroutine parse_tree_init_lexpr (parse_tree, stream, particles)
      type(parse_tree_t), intent(out) :: parse_tree
      type(stream_t), intent(inout), target :: stream
      logical, intent(in) :: particles
    end subroutine parse_tree_init_lexpr
    module subroutine parse_tree_init_pexpr (parse_tree, stream)
      type(parse_tree_t), intent(out) :: parse_tree
      type(stream_t), intent(inout), target :: stream
    end subroutine parse_tree_init_pexpr
    module subroutine parse_tree_init_cexpr (parse_tree, stream)
      type(parse_tree_t), intent(out) :: parse_tree
      type(stream_t), intent(inout), target :: stream
    end subroutine parse_tree_init_cexpr
    module subroutine parse_tree_init_sexpr (parse_tree, stream, particles)
      type(parse_tree_t), intent(out) :: parse_tree
      type(stream_t), intent(inout), target :: stream
      logical, intent(in) :: particles
    end subroutine parse_tree_init_sexpr
    module subroutine eval_tree_init_stream &
         (eval_tree, stream, var_list, subevt, result_type)
      class(eval_tree_t), intent(out), target :: eval_tree
      type(stream_t), intent(inout), target :: stream
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), target, optional :: subevt
      integer, intent(in), optional :: result_type
    end subroutine eval_tree_init_stream
    module subroutine eval_tree_init_expr &
        (expr, parse_node, var_list, subevt)
      class(eval_tree_t), intent(out), target :: expr
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
    end subroutine eval_tree_init_expr
    module subroutine eval_tree_init_lexpr &
        (expr, parse_node, var_list, subevt)
      class(eval_tree_t), intent(out), target :: expr
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
    end subroutine eval_tree_init_lexpr
    module subroutine eval_tree_init_pexpr &
        (expr, parse_node, var_list, subevt)
      class(eval_tree_t), intent(out), target :: expr
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
    end subroutine eval_tree_init_pexpr
    module subroutine eval_tree_init_cexpr &
        (expr, parse_node, var_list, subevt)
      class(eval_tree_t), intent(out), target :: expr
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
    end subroutine eval_tree_init_cexpr
    module subroutine eval_tree_init_sexpr &
        (expr, parse_node, var_list, subevt)
      class(eval_tree_t), intent(out), target :: expr
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
    end subroutine eval_tree_init_sexpr
    module subroutine eval_tree_setup_expr (expr, vars)
      class(eval_tree_t), intent(inout), target :: expr
      class(vars_t), intent(in), target :: vars
    end subroutine eval_tree_setup_expr
    module subroutine eval_tree_setup_lexpr (expr, vars)
      class(eval_tree_t), intent(inout), target :: expr
      class(vars_t), intent(in), target :: vars
    end subroutine eval_tree_setup_lexpr
    module subroutine eval_tree_setup_pexpr (expr, vars)
      class(eval_tree_t), intent(inout), target :: expr
      class(vars_t), intent(in), target :: vars
    end subroutine eval_tree_setup_pexpr
    module subroutine eval_tree_setup_cexpr (expr, vars)
      class(eval_tree_t), intent(inout), target :: expr
      class(vars_t), intent(in), target :: vars
    end subroutine eval_tree_setup_cexpr
    module subroutine eval_tree_setup_sexpr (expr, vars)
      class(eval_tree_t), intent(inout), target :: expr
      class(vars_t), intent(in), target :: vars
    end subroutine eval_tree_setup_sexpr
    module subroutine eval_tree_init_numeric_value (eval_tree, parse_node)
      class(eval_tree_t), intent(out), target :: eval_tree
      type(parse_node_t), intent(in), target :: parse_node
    end subroutine eval_tree_init_numeric_value
    module subroutine eval_tree_final (expr)
      class(eval_tree_t), intent(inout) :: expr
    end subroutine eval_tree_final
    module subroutine eval_tree_evaluate (expr)
      class(eval_tree_t), intent(inout) :: expr
    end subroutine eval_tree_evaluate
    module function eval_tree_get_result_type (expr) result (type)
      integer :: type
      class(eval_tree_t), intent(in) :: expr
    end function eval_tree_get_result_type
    module function eval_tree_result_is_known (expr) result (flag)
      logical :: flag
      class(eval_tree_t), intent(in) :: expr
    end function eval_tree_result_is_known
    module function eval_tree_result_is_known_ptr (expr) result (ptr)
      logical, pointer :: ptr
      class(eval_tree_t), intent(in) :: expr
    end function eval_tree_result_is_known_ptr
    module function eval_tree_get_log (expr) result (lval)
      logical :: lval
      class(eval_tree_t), intent(in) :: expr
    end function eval_tree_get_log
    module function eval_tree_get_int (expr) result (ival)
      integer :: ival
      class(eval_tree_t), intent(in) :: expr
    end function eval_tree_get_int
    module function eval_tree_get_real (expr) result (rval)
      real(default) :: rval
      class(eval_tree_t), intent(in) :: expr
    end function eval_tree_get_real
    module function eval_tree_get_cmplx (expr) result (cval)
      complex(default) :: cval
      class(eval_tree_t), intent(in) :: expr
    end function eval_tree_get_cmplx
    module function eval_tree_get_pdg_array (expr) result (aval)
      type(pdg_array_t) :: aval
      class(eval_tree_t), intent(in) :: expr
    end function eval_tree_get_pdg_array
    module function eval_tree_get_subevt (expr) result (pval)
      type(subevt_t) :: pval
      class(eval_tree_t), intent(in) :: expr
    end function eval_tree_get_subevt
    module function eval_tree_get_string (expr) result (sval)
      type(string_t) :: sval
      class(eval_tree_t), intent(in) :: expr
    end function eval_tree_get_string
    module subroutine eval_tree_write (expr, unit, write_vars)
      class(eval_tree_t), intent(in) :: expr
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: write_vars
    end subroutine eval_tree_write
    module function eval_log &
         (parse_node, var_list, subevt, is_known) result (lval)
      logical :: lval
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
      logical, intent(out), optional :: is_known
    end function eval_log
    module function eval_int &
         (parse_node, var_list, subevt, is_known) result (ival)
      integer :: ival
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
      logical, intent(out), optional :: is_known
    end function eval_int
    module function eval_real &
         (parse_node, var_list, subevt, is_known) result (rval)
      real(default) :: rval
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
      logical, intent(out), optional :: is_known
    end function eval_real
    module function eval_cmplx &
         (parse_node, var_list, subevt, is_known) result (cval)
      complex(default) :: cval
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
      logical, intent(out), optional :: is_known
    end function eval_cmplx
    module function eval_subevt &
         (parse_node, var_list, subevt, is_known) result (pval)
      type(subevt_t) :: pval
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
      logical, intent(out), optional :: is_known
    end function eval_subevt
    module function eval_pdg_array &
         (parse_node, var_list, subevt, is_known) result (aval)
      type(pdg_array_t) :: aval
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
      logical, intent(out), optional :: is_known
    end function eval_pdg_array
    module function eval_string &
         (parse_node, var_list, subevt, is_known) result (sval)
      type(string_t) :: sval
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
      logical, intent(out), optional :: is_known
    end function eval_string
    module subroutine eval_numeric &
         (parse_node, var_list, subevt, ival, rval, cval, &
          is_known, result_type)
      type(parse_node_t), intent(in), target :: parse_node
      type(var_list_t), intent(in), target :: var_list
      type(subevt_t), intent(in), optional, target :: subevt
      integer, intent(out), optional :: ival
      real(default), intent(out), optional :: rval
      complex(default), intent(out), optional :: cval
      logical, intent(out), optional :: is_known
      integer, intent(out), optional :: result_type
    end subroutine eval_numeric
    module subroutine eval_tree_factory_write (expr_factory, unit)
      class(eval_tree_factory_t), intent(in) :: expr_factory
      integer, intent(in), optional :: unit
    end subroutine eval_tree_factory_write
    module subroutine eval_tree_factory_init (expr_factory, pn)
      class(eval_tree_factory_t), intent(out) :: expr_factory
      type(parse_node_t), intent(in), pointer :: pn
    end subroutine eval_tree_factory_init
  end interface

contains

  subroutine eval_tree_factory_build (expr_factory, expr)
    class(eval_tree_factory_t), intent(in) :: expr_factory
    class(expr_t), intent(out), allocatable :: expr
    if (associated (expr_factory%pn)) then
       allocate (eval_tree_t :: expr)
       select type (expr)
       type is (eval_tree_t)
          expr%pn => expr_factory%pn
       end select
    end if
  end subroutine eval_tree_factory_build


end module eval_trees
