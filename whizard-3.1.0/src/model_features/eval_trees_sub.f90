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

submodule (eval_trees) eval_trees_s

  use io_units
  use constants, only: DEGREE, IMAGO, PI
  use format_defs, only: FMT_19
  use numeric_utils, only: nearly_equal
  use diagnostics
  use lorentz
  use md5
  use formats
  use sorting
  use analysis
  use jets
  use observables

  implicit none

contains

  recursive module subroutine eval_node_final_rec (node)
    class(eval_node_t), intent(inout) :: node
    select case (node%type)
    case (EN_UNARY)
       call eval_node_final_rec (node%arg1)
    case (EN_BINARY)
       call eval_node_final_rec (node%arg1)
       call eval_node_final_rec (node%arg2)
    case (EN_CONDITIONAL)
       call eval_node_final_rec (node%arg0)
       call eval_node_final_rec (node%arg1)
       call eval_node_final_rec (node%arg2)
    case (EN_BLOCK)
       call eval_node_final_rec (node%arg0)
       call eval_node_final_rec (node%arg1)
    case (EN_PRT_FUN_UNARY, EN_EVAL_FUN_UNARY, &
          EN_LOG_FUN_UNARY, EN_INT_FUN_UNARY, EN_REAL_FUN_UNARY)
       if (associated (node%arg0))  call eval_node_final_rec (node%arg0)
       call eval_node_final_rec (node%arg1)
       deallocate (node%index)
       deallocate (node%prt1)
    case (EN_REAL_FUN_CUM)
       if (associated (node%arg0))  call eval_node_final_rec (node%arg0)
       call eval_node_final_rec (node%arg1)
       deallocate (node%index)
       deallocate (node%prt1)
    case (EN_PRT_FUN_BINARY, EN_EVAL_FUN_BINARY, &
          EN_LOG_FUN_BINARY, EN_INT_FUN_BINARY, EN_REAL_FUN_BINARY)
       if (associated (node%arg0))  call eval_node_final_rec (node%arg0)
       call eval_node_final_rec (node%arg1)
       call eval_node_final_rec (node%arg2)
       deallocate (node%index)
       deallocate (node%prt1)
       deallocate (node%prt2)
    case (EN_FORMAT_STR)
       if (associated (node%arg0))  call eval_node_final_rec (node%arg0)
       if (associated (node%arg1))  call eval_node_final_rec (node%arg1)
       deallocate (node%ival)
    case (EN_RECORD_CMD)
       if (associated (node%arg0))  call eval_node_final_rec (node%arg0)
       if (associated (node%arg1))  call eval_node_final_rec (node%arg1)
       if (associated (node%arg2))  call eval_node_final_rec (node%arg2)
       if (associated (node%arg3))  call eval_node_final_rec (node%arg3)
       if (associated (node%arg4))  call eval_node_final_rec (node%arg4)
    end select
    select case (node%type)
    case (EN_UNARY, EN_BINARY, EN_CONDITIONAL, EN_CONSTANT, EN_BLOCK, &
          EN_PRT_FUN_UNARY, EN_PRT_FUN_BINARY, &
          EN_EVAL_FUN_UNARY, EN_EVAL_FUN_BINARY, &
          EN_LOG_FUN_UNARY, EN_LOG_FUN_BINARY, &
          EN_INT_FUN_UNARY, EN_INT_FUN_BINARY, &
          EN_REAL_FUN_UNARY, EN_REAL_FUN_BINARY, &
          EN_REAL_FUN_CUM, &
          EN_FORMAT_STR, EN_RECORD_CMD)
       select case (node%result_type)
       case (V_LOG);  deallocate (node%lval)
       case (V_INT);  deallocate (node%ival)
       case (V_REAL); deallocate (node%rval)
       case (V_CMPLX); deallocate (node%cval)
       case (V_SEV);  deallocate (node%pval)
       case (V_PDG);  deallocate (node%aval)
       case (V_STR);  deallocate (node%sval)
       end select
       deallocate (node%value_is_known)
    end select
  end subroutine eval_node_final_rec

  subroutine eval_node_init_log (node, lval)
    type(eval_node_t), intent(out) :: node
    logical, intent(in) :: lval
    node%type = EN_CONSTANT
    node%result_type = V_LOG
    allocate (node%lval, node%value_is_known)
    node%lval = lval
    node%value_is_known = .true.
  end subroutine eval_node_init_log

  subroutine eval_node_init_int (node, ival)
    type(eval_node_t), intent(out) :: node
    integer, intent(in) :: ival
    node%type = EN_CONSTANT
    node%result_type = V_INT
    allocate (node%ival, node%value_is_known)
    node%ival = ival
    node%value_is_known = .true.
  end subroutine eval_node_init_int

  subroutine eval_node_init_real (node, rval)
    type(eval_node_t), intent(out) :: node
    real(default), intent(in) :: rval
    node%type = EN_CONSTANT
    node%result_type = V_REAL
    allocate (node%rval, node%value_is_known)
    node%rval = rval
    node%value_is_known = .true.
  end subroutine eval_node_init_real

  subroutine eval_node_init_cmplx (node, cval)
    type(eval_node_t), intent(out) :: node
    complex(default), intent(in) :: cval
    node%type = EN_CONSTANT
    node%result_type = V_CMPLX
    allocate (node%cval, node%value_is_known)
    node%cval = cval
    node%value_is_known = .true.
  end subroutine eval_node_init_cmplx

  subroutine eval_node_init_subevt (node, pval)
    type(eval_node_t), intent(out) :: node
    type(subevt_t), intent(in) :: pval
    node%type = EN_CONSTANT
    node%result_type = V_SEV
    allocate (node%pval, node%value_is_known)
    node%pval = pval
    node%value_is_known = .true.
  end subroutine eval_node_init_subevt

  subroutine eval_node_init_pdg_array (node, aval)
    type(eval_node_t), intent(out) :: node
    type(pdg_array_t), intent(in) :: aval
    node%type = EN_CONSTANT
    node%result_type = V_PDG
    allocate (node%aval, node%value_is_known)
    node%aval = aval
    node%value_is_known = .true.
  end subroutine eval_node_init_pdg_array

  subroutine eval_node_init_string (node, sval)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: sval
    node%type = EN_CONSTANT
    node%result_type = V_STR
    allocate (node%sval, node%value_is_known)
    node%sval = sval
    node%value_is_known = .true.
  end subroutine eval_node_init_string

  subroutine eval_node_init_log_ptr (node, name, lval, is_known)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    logical, intent(in), target :: lval
    logical, intent(in), target :: is_known
    node%type = EN_VARIABLE
    node%tag = name
    node%result_type = V_LOG
    node%lval => lval
    node%value_is_known => is_known
  end subroutine eval_node_init_log_ptr

  subroutine eval_node_init_int_ptr (node, name, ival, is_known)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    integer, intent(in), target :: ival
    logical, intent(in), target :: is_known
    node%type = EN_VARIABLE
    node%tag = name
    node%result_type = V_INT
    node%ival => ival
    node%value_is_known => is_known
  end subroutine eval_node_init_int_ptr

  subroutine eval_node_init_real_ptr (node, name, rval, is_known)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    real(default), intent(in), target :: rval
    logical, intent(in), target :: is_known
    node%type = EN_VARIABLE
    node%tag = name
    node%result_type = V_REAL
    node%rval => rval
    node%value_is_known => is_known
  end subroutine eval_node_init_real_ptr

  subroutine eval_node_init_cmplx_ptr (node, name, cval, is_known)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    complex(default), intent(in), target :: cval
    logical, intent(in), target :: is_known
    node%type = EN_VARIABLE
    node%tag = name
    node%result_type = V_CMPLX
    node%cval => cval
    node%value_is_known => is_known
  end subroutine eval_node_init_cmplx_ptr

  subroutine eval_node_init_subevt_ptr (node, name, pval, is_known)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    type(subevt_t), intent(in), target :: pval
    logical, intent(in), target :: is_known
    node%type = EN_VARIABLE
    node%tag = name
    node%result_type = V_SEV
    node%pval => pval
    node%value_is_known => is_known
  end subroutine eval_node_init_subevt_ptr

  subroutine eval_node_init_pdg_array_ptr (node, name, aval, is_known)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    type(pdg_array_t), intent(in), target :: aval
    logical, intent(in), target :: is_known
    node%type = EN_VARIABLE
    node%tag = name
    node%result_type = V_PDG
    node%aval => aval
    node%value_is_known => is_known
  end subroutine eval_node_init_pdg_array_ptr

  subroutine eval_node_init_string_ptr (node, name, sval, is_known)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    type(string_t), intent(in), target :: sval
    logical, intent(in), target :: is_known
    node%type = EN_VARIABLE
    node%tag = name
    node%result_type = V_STR
    node%sval => sval
    node%value_is_known => is_known
  end subroutine eval_node_init_string_ptr

  subroutine eval_node_init_obs1_int_ptr (node, name, obs1_iptr, p1)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    procedure(obs_unary_int), intent(in), pointer :: obs1_iptr
    type(prt_t), intent(in), target :: p1
    node%type = EN_OBS1_INT
    node%tag = name
    node%result_type = V_INT
    node%obs1_int => obs1_iptr
    node%prt1 => p1
    allocate (node%ival, node%value_is_known)
    node%value_is_known = .false.
  end subroutine eval_node_init_obs1_int_ptr

  subroutine eval_node_init_obs2_int_ptr (node, name, obs2_iptr, p1, p2)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    procedure(obs_binary_int), intent(in), pointer :: obs2_iptr
    type(prt_t), intent(in), target :: p1, p2
    node%type = EN_OBS2_INT
    node%tag = name
    node%result_type = V_INT
    node%obs2_int => obs2_iptr
    node%prt1 => p1
    node%prt2 => p2
    allocate (node%ival, node%value_is_known)
    node%value_is_known = .false.
  end subroutine eval_node_init_obs2_int_ptr

  subroutine eval_node_init_obsev_int_ptr (node, name, obsev_iptr, pval)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    procedure(obs_sev_int), intent(in), pointer :: obsev_iptr
    type(subevt_t), intent(in), target :: pval
    node%type = EN_OBSEV_INT
    node%tag = name
    node%result_type = V_INT
    node%obsev_int => obsev_iptr
    node%pval => pval
    allocate (node%rval, node%value_is_known)
    node%value_is_known = .false.
  end subroutine eval_node_init_obsev_int_ptr

  subroutine eval_node_init_obs1_real_ptr (node, name, obs1_rptr, p1)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    procedure(obs_unary_real), intent(in), pointer :: obs1_rptr
    type(prt_t), intent(in), target :: p1
    node%type = EN_OBS1_REAL
    node%tag = name
    node%result_type = V_REAL
    node%obs1_real => obs1_rptr
    node%prt1 => p1
    allocate (node%rval, node%value_is_known)
    node%value_is_known = .false.
  end subroutine eval_node_init_obs1_real_ptr

  subroutine eval_node_init_obs2_real_ptr (node, name, obs2_rptr, p1, p2)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    procedure(obs_binary_real), intent(in), pointer :: obs2_rptr
    type(prt_t), intent(in), target :: p1, p2
    node%type = EN_OBS2_REAL
    node%tag = name
    node%result_type = V_REAL
    node%obs2_real => obs2_rptr
    node%prt1 => p1
    node%prt2 => p2
    allocate (node%rval, node%value_is_known)
    node%value_is_known = .false.
  end subroutine eval_node_init_obs2_real_ptr

  subroutine eval_node_init_obsev_real_ptr (node, name, obsev_rptr, pval)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: name
    procedure(obs_sev_real), intent(in), pointer :: obsev_rptr
    type(subevt_t), intent(in), target :: pval
    node%type = EN_OBSEV_REAL
    node%tag = name
    node%result_type = V_REAL
    node%obsev_real => obsev_rptr
    node%pval => pval
    allocate (node%rval, node%value_is_known)
    node%value_is_known = .false.
  end subroutine eval_node_init_obsev_real_ptr

  subroutine eval_node_init_branch (node, tag, result_type, arg1, arg2)
    type(eval_node_t), intent(out) :: node
    type(string_t), intent(in) :: tag
    integer, intent(in) :: result_type
    type(eval_node_t), intent(in), target :: arg1
    type(eval_node_t), intent(in), target, optional :: arg2
    if (present (arg2)) then
       node%type = EN_BINARY
    else
       node%type = EN_UNARY
    end if
    node%tag = tag
    node%result_type = result_type
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    if (present (arg2))  node%arg2 => arg2
  end subroutine eval_node_init_branch

  subroutine eval_node_allocate_value (node)
    type(eval_node_t), intent(inout) :: node
    select case (node%result_type)
    case (V_LOG);  allocate (node%lval)
    case (V_INT);  allocate (node%ival)
    case (V_REAL); allocate (node%rval)
    case (V_CMPLX); allocate (node%cval)
    case (V_PDG);  allocate (node%aval)
    case (V_SEV);  allocate (node%pval)
       call subevt_init (node%pval)
    case (V_STR);  allocate (node%sval)
    end select
    allocate (node%value_is_known)
    node%value_is_known = .false.
  end subroutine eval_node_allocate_value

  subroutine eval_node_init_block (node, name, type, var_def, var_list)
    type(eval_node_t), intent(out), target :: node
    type(string_t), intent(in) :: name
    integer, intent(in) :: type
    type(eval_node_t), intent(in), target :: var_def
    type(var_list_t), intent(in), target :: var_list
    node%type = EN_BLOCK
    node%tag = "var_def"
    node%var_name = name
    node%arg1 => var_def
    allocate (node%var_list)
    call node%var_list%link (var_list)
    if (var_def%type == EN_CONSTANT) then
       select case (type)
       case (V_LOG)
          call node%var_list%append_log (name, var_def%lval)
       case (V_INT)
          call node%var_list%append_int (name, var_def%ival)
       case (V_REAL)
          call node%var_list%append_real (name, var_def%rval)
       case (V_CMPLX)
          call node%var_list%append_cmplx (name, var_def%cval)
       case (V_PDG)
          call node%var_list%append_pdg_array (name, var_def%aval)
       case (V_SEV)
          call node%var_list%append_subevt (name, var_def%pval)
       case (V_STR)
          call node%var_list%append_string (name, var_def%sval)
       end select
    else
       select case (type)
       case (V_LOG);  call node%var_list%append_log_ptr &
            (name, var_def%lval, var_def%value_is_known)
       case (V_INT);  call node%var_list%append_int_ptr &
            (name, var_def%ival, var_def%value_is_known)
       case (V_REAL); call node%var_list%append_real_ptr &
            (name, var_def%rval, var_def%value_is_known)
       case (V_CMPLX); call node%var_list%append_cmplx_ptr &
            (name, var_def%cval, var_def%value_is_known)
       case (V_PDG);  call node%var_list%append_pdg_array_ptr &
            (name, var_def%aval, var_def%value_is_known)
       case (V_SEV); call node%var_list%append_subevt_ptr &
            (name, var_def%pval, var_def%value_is_known)
       case (V_STR); call node%var_list%append_string_ptr &
            (name, var_def%sval, var_def%value_is_known)
       end select
    end if
  end subroutine eval_node_init_block

  subroutine eval_node_set_expr (node, arg, result_type)
    type(eval_node_t), intent(inout) :: node
    type(eval_node_t), intent(in), target :: arg
    integer, intent(in), optional :: result_type
    if (present (result_type)) then
       node%result_type = result_type
    else
       node%result_type = arg%result_type
    end if
    call eval_node_allocate_value (node)
    node%arg0 => arg
  end subroutine eval_node_set_expr

  subroutine eval_node_init_conditional (node, result_type, cond, arg1, arg2)
    type(eval_node_t), intent(out) :: node
    integer, intent(in) :: result_type
    type(eval_node_t), intent(in), target :: cond, arg1, arg2
    node%type = EN_CONDITIONAL
    node%tag = "cond"
    node%result_type = result_type
    call eval_node_allocate_value (node)
    node%arg0 => cond
    node%arg1 => arg1
    node%arg2 => arg2
  end subroutine eval_node_init_conditional

  subroutine eval_node_init_record_cmd &
      (node, event_weight, id, arg1, arg2, arg3, arg4)
    type(eval_node_t), intent(out) :: node
    real(default), pointer :: event_weight
    type(eval_node_t), intent(in), target :: id
    type(eval_node_t), intent(in), optional, target :: arg1, arg2, arg3, arg4
    call eval_node_init_log (node, .true.)
    node%type = EN_RECORD_CMD
    node%rval => event_weight
    node%tag = "record_cmd"
    node%arg0 => id
    if (present (arg1)) then
       node%arg1 => arg1
       if (present (arg2)) then
          node%arg2 => arg2
          if (present (arg3)) then
             node%arg3 => arg3
             if (present (arg4)) then
                node%arg4 => arg4
             end if
          end if
       end if
    end if
  end subroutine eval_node_init_record_cmd

  subroutine eval_node_init_prt_fun_unary (node, arg1, name, proc)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1
    type(string_t), intent(in) :: name
    procedure(unary_sev) :: proc
    node%type = EN_PRT_FUN_UNARY
    node%tag = name
    node%result_type = V_SEV
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    allocate (node%index, source = 0)
    allocate (node%prt1)
    node%op1_sev => proc
  end subroutine eval_node_init_prt_fun_unary

  subroutine eval_node_init_prt_fun_binary (node, arg1, arg2, name, proc)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1, arg2
    type(string_t), intent(in) :: name
    procedure(binary_sev) :: proc
    node%type = EN_PRT_FUN_BINARY
    node%tag = name
    node%result_type = V_SEV
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    node%arg2 => arg2
    allocate (node%index, source = 0)
    allocate (node%prt1)
    allocate (node%prt2)
    node%op2_sev => proc
  end subroutine eval_node_init_prt_fun_binary

  subroutine eval_node_init_eval_fun_unary (node, arg1, name)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1
    type(string_t), intent(in) :: name
    node%type = EN_EVAL_FUN_UNARY
    node%tag = name
    node%result_type = V_REAL
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    allocate (node%index, source = 0)
    allocate (node%prt1)
  end subroutine eval_node_init_eval_fun_unary

  subroutine eval_node_init_eval_fun_binary (node, arg1, arg2, name)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1, arg2
    type(string_t), intent(in) :: name
    node%type = EN_EVAL_FUN_BINARY
    node%tag = name
    node%result_type = V_REAL
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    node%arg2 => arg2
    allocate (node%index, source = 0)
    allocate (node%prt1)
    allocate (node%prt2)
  end subroutine eval_node_init_eval_fun_binary

  subroutine eval_node_init_log_fun_unary (node, arg1, name, proc)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1
    type(string_t), intent(in) :: name
    procedure(unary_cut) :: proc
    node%type = EN_LOG_FUN_UNARY
    node%tag = name
    node%result_type = V_LOG
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    allocate (node%index, source = 0)
    allocate (node%prt1)
    node%op1_cut => proc
  end subroutine eval_node_init_log_fun_unary

  subroutine eval_node_init_log_fun_binary (node, arg1, arg2, name, proc)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1, arg2
    type(string_t), intent(in) :: name
    procedure(binary_cut) :: proc
    node%type = EN_LOG_FUN_BINARY
    node%tag = name
    node%result_type = V_LOG
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    node%arg2 => arg2
    allocate (node%index, source = 0)
    allocate (node%prt1)
    allocate (node%prt2)
    node%op2_cut => proc
  end subroutine eval_node_init_log_fun_binary

  subroutine eval_node_init_int_fun_unary (node, arg1, name, proc)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1
    type(string_t), intent(in) :: name
    procedure(unary_evi) :: proc
    node%type = EN_INT_FUN_UNARY
    node%tag = name
    node%result_type = V_INT
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    allocate (node%index, source = 0)
    allocate (node%prt1)
    node%op1_evi => proc
  end subroutine eval_node_init_int_fun_unary

  subroutine eval_node_init_int_fun_binary (node, arg1, arg2, name, proc)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1, arg2
    type(string_t), intent(in) :: name
    procedure(binary_evi) :: proc
    node%type = EN_INT_FUN_BINARY
    node%tag = name
    node%result_type = V_INT
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    node%arg2 => arg2
    allocate (node%index, source = 0)
    allocate (node%prt1)
    allocate (node%prt2)
    node%op2_evi => proc
  end subroutine eval_node_init_int_fun_binary

  subroutine eval_node_init_real_fun_unary (node, arg1, name, proc)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1
    type(string_t), intent(in) :: name
    procedure(unary_evr) :: proc
    node%type = EN_REAL_FUN_UNARY
    node%tag = name
    node%result_type = V_REAL
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    allocate (node%index, source = 0)
    allocate (node%prt1)
    node%op1_evr => proc
  end subroutine eval_node_init_real_fun_unary

  subroutine eval_node_init_real_fun_binary (node, arg1, arg2, name, proc)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1, arg2
    type(string_t), intent(in) :: name
    procedure(binary_evr) :: proc
    node%type = EN_REAL_FUN_BINARY
    node%tag = name
    node%result_type = V_REAL
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    node%arg2 => arg2
    allocate (node%index, source = 0)
    allocate (node%prt1)
    allocate (node%prt2)
    node%op2_evr => proc
  end subroutine eval_node_init_real_fun_binary

  subroutine eval_node_init_real_fun_cum (node, arg1, name, proc)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), intent(in), target :: arg1
    type(string_t), intent(in) :: name
    procedure(cum_evr) :: proc
    node%type = EN_REAL_FUN_CUM
    node%tag = name
    node%result_type = V_REAL
    call eval_node_allocate_value (node)
    node%arg1 => arg1
    allocate (node%index, source = 0)
    allocate (node%prt1)
    node%opcum_evr => proc
  end subroutine eval_node_init_real_fun_cum

  subroutine eval_node_init_format_string (node, fmt, arg, name, n_args)
    type(eval_node_t), intent(out) :: node
    type(eval_node_t), pointer :: fmt, arg
    type(string_t), intent(in) :: name
    integer, intent(in) :: n_args
    node%type = EN_FORMAT_STR
    node%tag = name
    node%result_type = V_STR
    call eval_node_allocate_value (node)
    node%arg0 => fmt
    node%arg1 => arg
    allocate (node%ival)
    node%ival = n_args
  end subroutine eval_node_init_format_string

  subroutine eval_node_set_observables (node, var_list)
    type(eval_node_t), intent(inout) :: node
    type(var_list_t), intent(in), target :: var_list
    logical, save, target :: known = .true.
    allocate (node%var_list)
    call node%var_list%link (var_list)
    allocate (node%index, source = 0)
    call node%var_list%append_int_ptr &
         (var_str ("Index"), node%index, known, intrinsic=.true.)
    if (.not. associated (node%prt2)) then
       call var_list_set_observables_unary &
            (node%var_list, node%prt1)
       if (associated (node%pval)) then
          call var_list_set_observables_sev &
               (node%var_list, node%pval)
       end if
    else
       call var_list_set_observables_binary &
            (node%var_list, node%prt1, node%prt2)
    end if
  end subroutine eval_node_set_observables

  module subroutine eval_node_write (node, unit, indent)
    class(eval_node_t), intent(in) :: node
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: indent
    integer :: u, ind
    u = given_output_unit (unit);  if (u < 0)  return
    ind = 0;  if (present (indent)) ind = indent
    write (u, "(A)", advance="no")  repeat ("|  ", ind) // "o "
    select case (node%type)
    case (EN_UNARY, EN_BINARY, EN_CONDITIONAL, &
          EN_PRT_FUN_UNARY, EN_PRT_FUN_BINARY, &
          EN_EVAL_FUN_UNARY, EN_EVAL_FUN_BINARY, &
          EN_LOG_FUN_UNARY, EN_LOG_FUN_BINARY, &
          EN_INT_FUN_UNARY, EN_INT_FUN_BINARY, &
          EN_REAL_FUN_UNARY, EN_REAL_FUN_BINARY, &
          EN_REAL_FUN_CUM)
       write (u, "(A)", advance="no")  "[" // char (node%tag) // "] ="
    case (EN_CONSTANT)
       write (u, "(A)", advance="no")  "[const] ="
    case (EN_VARIABLE)
       write (u, "(A)", advance="no")  char (node%tag) // " =>"
    case (EN_OBS1_INT, EN_OBS2_INT, EN_OBS1_REAL, EN_OBS2_REAL)
       write (u, "(A)", advance="no")  char (node%tag) // " ="
    case (EN_BLOCK)
       write (u, "(A)", advance="no")  "[" // char (node%tag) // "]" // &
            char (node%var_name) // " [expr] = "
    case default
       write (u, "(A)", advance="no")  "[???] ="
    end select
    select case (node%result_type)
    case (V_LOG)
       if (node%value_is_known) then
          if (node%lval) then
             write (u, "(1x,A)") "true"
          else
             write (u, "(1x,A)") "false"
          end if
       else
          write (u, "(1x,A)") "[unknown logical]"
       end if
    case (V_INT)
       if (node%value_is_known) then
          write (u, "(1x,I0)")  node%ival
       else
          write (u, "(1x,A)") "[unknown integer]"
       end if
    case (V_REAL)
       if (node%value_is_known) then
          write (u, "(1x," // FMT_19 // ")") node%rval
       else
          write (u, "(1x,A)") "[unknown real]"
       end if
   case (V_CMPLX)
       if (node%value_is_known) then
          write (u, "(1x,'('," // FMT_19 // ",','," // &
               FMT_19 // ",')')") node%cval
       else
          write (u, "(1x,A)") "[unknown complex]"
       end if
    case (V_SEV)
       if (char (node%tag) == "@evt") then
          write (u, "(1x,A)") "[event subevent]"
       else if (node%value_is_known) then
          call node%pval%write (unit, prefix = repeat ("|  ", ind + 1))
       else
          write (u, "(1x,A)") "[unknown subevent]"
       end if
    case (V_PDG)
       write (u, "(1x)", advance="no")
       call node%aval%write (u);  write (u, *)
    case (V_STR)
       if (node%value_is_known) then
          write (u, "(A)")  '"' // char (node%sval) // '"'
       else
          write (u, "(1x,A)") "[unknown string]"
       end if
    case default
       write (u, "(1x,A)") "[empty]"
    end select
    select case (node%type)
    case (EN_OBS1_INT, EN_OBS1_REAL)
       write (u, "(A,6x,A)", advance="no")  repeat ("|  ", ind), "prt1 ="
       call prt_write (node%prt1, unit)
    case (EN_OBS2_INT, EN_OBS2_REAL)
       write (u, "(A,6x,A)", advance="no")  repeat ("|  ", ind), "prt1 ="
       call prt_write (node%prt1, unit)
       write (u, "(A,6x,A)", advance="no")  repeat ("|  ", ind), "prt2 ="
       call prt_write (node%prt2, unit)
    end select
  end subroutine eval_node_write

  recursive subroutine eval_node_write_rec (node, unit, indent)
    type(eval_node_t), intent(in) :: node
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: indent
    integer :: u, ind
    u = given_output_unit (unit);  if (u < 0)  return
    ind = 0;  if (present (indent))  ind = indent
    call eval_node_write (node, unit, indent)
    select case (node%type)
    case (EN_UNARY)
       if (associated (node%arg0)) &
            call eval_node_write_rec (node%arg0, unit, ind+1)
       call eval_node_write_rec (node%arg1, unit, ind+1)
    case (EN_BINARY)
       if (associated (node%arg0)) &
            call eval_node_write_rec (node%arg0, unit, ind+1)
       call eval_node_write_rec (node%arg1, unit, ind+1)
       call eval_node_write_rec (node%arg2, unit, ind+1)
    case (EN_BLOCK)
       call eval_node_write_rec (node%arg1, unit, ind+1)
       call eval_node_write_rec (node%arg0, unit, ind+1)
    case (EN_CONDITIONAL)
       call eval_node_write_rec (node%arg0, unit, ind+1)
       call eval_node_write_rec (node%arg1, unit, ind+1)
       call eval_node_write_rec (node%arg2, unit, ind+1)
    case (EN_PRT_FUN_UNARY, EN_EVAL_FUN_UNARY, &
          EN_LOG_FUN_UNARY, EN_INT_FUN_UNARY, EN_REAL_FUN_UNARY, &
          EN_REAL_FUN_CUM)
       if (associated (node%arg0)) &
            call eval_node_write_rec (node%arg0, unit, ind+1)
       call eval_node_write_rec (node%arg1, unit, ind+1)
    case (EN_PRT_FUN_BINARY, EN_EVAL_FUN_BINARY, &
          EN_LOG_FUN_BINARY, EN_INT_FUN_BINARY, EN_REAL_FUN_BINARY)
       if (associated (node%arg0)) &
            call eval_node_write_rec (node%arg0, unit, ind+1)
       call eval_node_write_rec (node%arg1, unit, ind+1)
       call eval_node_write_rec (node%arg2, unit, ind+1)
    case (EN_RECORD_CMD)
       if (associated (node%arg1)) then
          call eval_node_write_rec (node%arg1, unit, ind+1)
          if (associated (node%arg2)) then
             call eval_node_write_rec (node%arg2, unit, ind+1)
             if (associated (node%arg3)) then
                call eval_node_write_rec (node%arg3, unit, ind+1)
                if (associated (node%arg4)) then
                   call eval_node_write_rec (node%arg4, unit, ind+1)
                end if
             end if
          end if
       end if
    end select
  end subroutine eval_node_write_rec

  subroutine eval_node_set_op1_log (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(unary_log) :: op
    en%op1_log => op
    end subroutine eval_node_set_op1_log

  subroutine eval_node_set_op1_int (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(unary_int) :: op
    en%op1_int => op
  end subroutine eval_node_set_op1_int

  subroutine eval_node_set_op1_real (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(unary_real) :: op
    en%op1_real => op
  end subroutine eval_node_set_op1_real

  subroutine eval_node_set_op1_cmplx (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(unary_cmplx) :: op
    en%op1_cmplx => op
  end subroutine eval_node_set_op1_cmplx

  subroutine eval_node_set_op1_pdg (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(unary_pdg) :: op
    en%op1_pdg => op
  end subroutine eval_node_set_op1_pdg

  subroutine eval_node_set_op1_sev (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(unary_sev) :: op
    en%op1_sev => op
  end subroutine eval_node_set_op1_sev

  subroutine eval_node_set_op1_str (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(unary_str) :: op
    en%op1_str => op
  end subroutine eval_node_set_op1_str

  subroutine eval_node_set_op2_log (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(binary_log) :: op
    en%op2_log => op
  end subroutine eval_node_set_op2_log

  subroutine eval_node_set_op2_int (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(binary_int) :: op
    en%op2_int => op
  end subroutine eval_node_set_op2_int

  subroutine eval_node_set_op2_real (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(binary_real) :: op
    en%op2_real => op
  end subroutine eval_node_set_op2_real

  subroutine eval_node_set_op2_cmplx (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(binary_cmplx) :: op
    en%op2_cmplx => op
  end subroutine eval_node_set_op2_cmplx

  subroutine eval_node_set_op2_pdg (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(binary_pdg) :: op
    en%op2_pdg => op
  end subroutine eval_node_set_op2_pdg

  subroutine eval_node_set_op2_sev (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(binary_sev) :: op
    en%op2_sev => op
  end subroutine eval_node_set_op2_sev

  subroutine eval_node_set_op2_str (en, op)
    type(eval_node_t), intent(inout) :: en
    procedure(binary_str) :: op
    en%op2_str => op
  end subroutine eval_node_set_op2_str

  integer function add_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival + en2%ival
  end function add_ii
  real(default) function add_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival + en2%rval
  end function add_ir
  complex(default) function add_ic (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival + en2%cval
  end function add_ic
  real(default) function add_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval + en2%ival
  end function add_ri
  complex(default) function add_ci (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval + en2%ival
  end function add_ci
  complex(default) function add_cr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval + en2%rval
  end function add_cr
  complex(default) function add_rc (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval + en2%cval
  end function add_rc
  real(default) function add_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval + en2%rval
  end function add_rr
  complex(default) function add_cc (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval + en2%cval
  end function add_cc

  integer function sub_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival - en2%ival
  end function sub_ii
  real(default) function sub_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival - en2%rval
  end function sub_ir
  real(default) function sub_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval - en2%ival
  end function sub_ri
  complex(default) function sub_ic (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival - en2%cval
  end function sub_ic
  complex(default) function sub_ci (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval - en2%ival
  end function sub_ci
  complex(default) function sub_cr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval - en2%rval
  end function sub_cr
  complex(default) function sub_rc (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval - en2%cval
  end function sub_rc
  real(default) function sub_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval - en2%rval
  end function sub_rr
  complex(default) function sub_cc (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval - en2%cval
  end function sub_cc

  integer function mul_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival * en2%ival
  end function mul_ii
  real(default) function mul_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival * en2%rval
  end function mul_ir
  real(default) function mul_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval * en2%ival
  end function mul_ri
  complex(default) function mul_ic (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival * en2%cval
  end function mul_ic
  complex(default) function mul_ci (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval * en2%ival
  end function mul_ci
  complex(default) function mul_rc (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval * en2%cval
  end function mul_rc
  complex(default) function mul_cr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval * en2%rval
  end function mul_cr
  real(default) function mul_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval * en2%rval
  end function mul_rr
  complex(default) function mul_cc (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval * en2%cval
  end function mul_cc

  integer function div_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (en2%ival == 0) then
       if (en1%ival >= 0) then
          call msg_warning ("division by zero: " // int2char (en1%ival) // &
             " / 0 ; result set to 0")
       else
          call msg_warning ("division by zero: (" // int2char (en1%ival) // &
             ") / 0 ; result set to 0")
       end if
       y = 0
       return
    end if
    y = en1%ival / en2%ival
  end function div_ii
  real(default) function div_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival / en2%rval
  end function div_ir
  real(default) function div_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval / en2%ival
  end function div_ri
  complex(default) function div_ic (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival / en2%cval
  end function div_ic
  complex(default) function div_ci (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval / en2%ival
  end function div_ci
  complex(default) function div_rc (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval / en2%cval
  end function div_rc
  complex(default) function div_cr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval / en2%rval
  end function div_cr
  real(default) function div_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval / en2%rval
  end function div_rr
  complex(default) function div_cc (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval / en2%cval
  end function div_cc

  integer function pow_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    integer :: a, b
    real(default) :: rres
    a = en1%ival
    b = en2%ival
    if ((a == 0) .and. (b < 0)) then
       call msg_warning ("division by zero: " // int2char (a) // &
          " ^ (" // int2char (b) // ") ; result set to 0")
       y = 0
       return
    end if
    rres = real(a, default) ** b
    y = rres
    if (real(y, default) /= rres) then
       if (b < 0) then
          call msg_warning ("result of all-integer operation " // &
                int2char (a) // " ^ (" // int2char (b) // &
                ") has been trucated to "// int2char (y), &
             [  var_str ("Chances are that you want to use " // &
                "reals instead of integers at this point.") ])
       else
          call msg_warning ("integer overflow in " // int2char (a) // &
                " ^ " // int2char (b) // " ; result is " // int2char (y), &
             [  var_str ("Using reals instead of integers might help.")])
       end if
    end if
  end function pow_ii
  real(default) function pow_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval ** en2%ival
  end function pow_ri
  complex(default) function pow_ci (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval ** en2%ival
  end function pow_ci
  real(default) function pow_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival ** en2%rval
  end function pow_ir
  real(default) function pow_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval ** en2%rval
  end function pow_rr
  complex(default) function pow_cr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval ** en2%rval
  end function pow_cr
  complex(default) function pow_ic (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival ** en2%cval
  end function pow_ic
  complex(default) function pow_rc (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval ** en2%cval
  end function pow_rc
  complex(default) function pow_cc (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%cval ** en2%cval
  end function pow_cc

  integer function max_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = max (en1%ival, en2%ival)
  end function max_ii
  real(default) function max_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = max (real (en1%ival, default), en2%rval)
  end function max_ir
  real(default) function max_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = max (en1%rval, real (en2%ival, default))
  end function max_ri
  real(default) function max_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = max (en1%rval, en2%rval)
  end function max_rr
  integer function min_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = min (en1%ival, en2%ival)
  end function min_ii
  real(default) function min_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = min (real (en1%ival, default), en2%rval)
  end function min_ir
  real(default) function min_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = min (en1%rval, real (en2%ival, default))
  end function min_ri
  real(default) function min_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = min (en1%rval, en2%rval)
  end function min_rr

  integer function mod_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = mod (en1%ival, en2%ival)
  end function mod_ii
  real(default) function mod_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = mod (real (en1%ival, default), en2%rval)
  end function mod_ir
  real(default) function mod_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = mod (en1%rval, real (en2%ival, default))
  end function mod_ri
  real(default) function mod_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = mod (en1%rval, en2%rval)
  end function mod_rr
  integer function modulo_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = modulo (en1%ival, en2%ival)
  end function modulo_ii
  real(default) function modulo_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = modulo (real (en1%ival, default), en2%rval)
  end function modulo_ir
  real(default) function modulo_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = modulo (en1%rval, real (en2%ival, default))
  end function modulo_ri
  real(default) function modulo_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = modulo (en1%rval, en2%rval)
  end function modulo_rr

  real(default) function real_i (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = en%ival
  end function real_i
  real(default) function real_c (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = en%cval
  end function real_c
  integer function int_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = en%rval
  end function int_r
  complex(default) function cmplx_i (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = en%ival
  end function cmplx_i
  integer function int_c (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = en%cval
  end function int_c
  complex(default) function cmplx_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = en%rval
  end function cmplx_r
  integer function nint_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = nint (en%rval)
  end function nint_r
  integer function floor_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = floor (en%rval)
  end function floor_r
  integer function ceiling_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = ceiling (en%rval)
  end function ceiling_r

  integer function neg_i (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = - en%ival
  end function neg_i
  real(default) function neg_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = - en%rval
  end function neg_r
  complex(default) function neg_c (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = - en%cval
  end function neg_c
  integer function abs_i (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = abs (en%ival)
  end function abs_i
  real(default) function abs_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = abs (en%rval)
  end function abs_r
  real(default) function abs_c (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = abs (en%cval)
  end function abs_c
  integer function conjg_i (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = en%ival
  end function conjg_i
  real(default) function conjg_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = en%rval
  end function conjg_r
  complex(default) function conjg_c (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = conjg (en%cval)
  end function conjg_c
  integer function sgn_i (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = sign (1, en%ival)
  end function sgn_i
  real(default) function sgn_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = sign (1._default, en%rval)
  end function sgn_r

  real(default) function sqrt_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = sqrt (en%rval)
  end function sqrt_r
  real(default) function exp_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = exp (en%rval)
  end function exp_r
  real(default) function log_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = log (en%rval)
  end function log_r
  real(default) function log10_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = log10 (en%rval)
  end function log10_r

  complex(default) function sqrt_c (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = sqrt (en%cval)
  end function sqrt_c
  complex(default) function exp_c (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = exp (en%cval)
  end function exp_c
  complex(default) function log_c (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = log (en%cval)
  end function log_c

  real(default) function sin_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = sin (en%rval)
  end function sin_r
  real(default) function cos_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = cos (en%rval)
  end function cos_r
  real(default) function tan_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = tan (en%rval)
  end function tan_r
  real(default) function asin_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = asin (en%rval)
  end function asin_r
  real(default) function acos_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = acos (en%rval)
  end function acos_r
  real(default) function atan_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = atan (en%rval)
  end function atan_r

  complex(default) function sin_c (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = sin (en%cval)
  end function sin_c
  complex(default) function cos_c (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = cos (en%cval)
  end function cos_c

  real(default) function sinh_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = sinh (en%rval)
  end function sinh_r
  real(default) function cosh_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = cosh (en%rval)
  end function cosh_r
  real(default) function tanh_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = tanh (en%rval)
  end function tanh_r
  real(default) function asinh_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = asinh (en%rval)
  end function asinh_r
  real(default) function acosh_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = acosh (en%rval)
  end function acosh_r
  real(default) function atanh_r (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = atanh (en%rval)
  end function atanh_r

  logical function ignore_first_ll (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en2%lval
  end function ignore_first_ll
  logical function or_ll (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%lval .or. en2%lval
  end function or_ll
  logical function and_ll (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%lval .and. en2%lval
  end function and_ll

  logical function comp_lt_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival < en2%ival
  end function comp_lt_ii
  logical function comp_lt_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival < en2%rval
  end function comp_lt_ir
  logical function comp_lt_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval < en2%ival
  end function comp_lt_ri
  logical function comp_lt_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval < en2%rval
  end function comp_lt_rr

  logical function comp_gt_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival > en2%ival
  end function comp_gt_ii
  logical function comp_gt_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival > en2%rval
  end function comp_gt_ir
  logical function comp_gt_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval > en2%ival
  end function comp_gt_ri
  logical function comp_gt_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval > en2%rval
  end function comp_gt_rr

  logical function comp_le_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival <= en2%ival
  end function comp_le_ii
  logical function comp_le_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival <= en2%rval
  end function comp_le_ir
  logical function comp_le_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval <= en2%ival
  end function comp_le_ri
  logical function comp_le_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval <= en2%rval
  end function comp_le_rr

  logical function comp_ge_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival >= en2%ival
  end function comp_ge_ii
  logical function comp_ge_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival >= en2%rval
  end function comp_ge_ir
  logical function comp_ge_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval >= en2%ival
  end function comp_ge_ri
  logical function comp_ge_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval >= en2%rval
  end function comp_ge_rr

  logical function comp_eq_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival == en2%ival
  end function comp_eq_ii
  logical function comp_eq_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival == en2%rval
  end function comp_eq_ir
  logical function comp_eq_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval == en2%ival
  end function comp_eq_ri
  logical function comp_eq_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval == en2%rval
  end function comp_eq_rr
  logical function comp_eq_ss (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%sval == en2%sval
  end function comp_eq_ss

  logical function comp_ne_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival /= en2%ival
  end function comp_ne_ii
  logical function comp_ne_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%ival /= en2%rval
  end function comp_ne_ir
  logical function comp_ne_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval /= en2%ival
  end function comp_ne_ri
  logical function comp_ne_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%rval /= en2%rval
  end function comp_ne_rr
  logical function comp_ne_ss (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    y = en1%sval /= en2%sval
  end function comp_ne_ss

  logical function comp_se_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = abs (en1%ival - en2%ival) <= en1%tolerance
    else
       y = en1%ival == en2%ival
    end if
  end function comp_se_ii
  logical function comp_se_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = abs (en1%rval - en2%ival) <= en1%tolerance
    else
       y = en1%rval == en2%ival
    end if
  end function comp_se_ri
  logical function comp_se_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = abs (en1%ival - en2%rval) <= en1%tolerance
    else
       y = en1%ival == en2%rval
    end if
  end function comp_se_ir
  logical function comp_se_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = abs (en1%rval - en2%rval) <= en1%tolerance
    else
       y = en1%rval == en2%rval
    end if
  end function comp_se_rr
  logical function comp_ns_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = abs (en1%ival - en2%ival) > en1%tolerance
    else
       y = en1%ival /= en2%ival
    end if
  end function comp_ns_ii
  logical function comp_ns_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = abs (en1%rval - en2%ival) > en1%tolerance
    else
       y = en1%rval /= en2%ival
    end if
  end function comp_ns_ri
  logical function comp_ns_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = abs (en1%ival - en2%rval) > en1%tolerance
    else
       y = en1%ival /= en2%rval
    end if
  end function comp_ns_ir
  logical function comp_ns_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = abs (en1%rval - en2%rval) > en1%tolerance
    else
       y = en1%rval /= en2%rval
    end if
  end function comp_ns_rr

  logical function comp_ls_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%ival <= en2%ival + en1%tolerance
    else
       y = en1%ival <= en2%ival
    end if
  end function comp_ls_ii
  logical function comp_ls_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%rval <= en2%ival + en1%tolerance
    else
       y = en1%rval <= en2%ival
    end if
  end function comp_ls_ri
  logical function comp_ls_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%ival <= en2%rval + en1%tolerance
    else
       y = en1%ival <= en2%rval
    end if
  end function comp_ls_ir
  logical function comp_ls_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%rval <= en2%rval + en1%tolerance
    else
       y = en1%rval <= en2%rval
    end if
  end function comp_ls_rr

  logical function comp_ll_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%ival < en2%ival - en1%tolerance
    else
       y = en1%ival < en2%ival
    end if
  end function comp_ll_ii
  logical function comp_ll_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%rval < en2%ival - en1%tolerance
    else
       y = en1%rval < en2%ival
    end if
  end function comp_ll_ri
  logical function comp_ll_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%ival < en2%rval - en1%tolerance
    else
       y = en1%ival < en2%rval
    end if
  end function comp_ll_ir
  logical function comp_ll_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%rval < en2%rval - en1%tolerance
    else
       y = en1%rval < en2%rval
    end if
  end function comp_ll_rr

  logical function comp_gs_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%ival >= en2%ival - en1%tolerance
    else
       y = en1%ival >= en2%ival
    end if
  end function comp_gs_ii
  logical function comp_gs_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%rval >= en2%ival - en1%tolerance
    else
       y = en1%rval >= en2%ival
    end if
  end function comp_gs_ri
  logical function comp_gs_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%ival >= en2%rval - en1%tolerance
    else
       y = en1%ival >= en2%rval
    end if
  end function comp_gs_ir
  logical function comp_gs_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%rval >= en2%rval - en1%tolerance
    else
       y = en1%rval >= en2%rval
    end if
  end function comp_gs_rr

  logical function comp_gg_ii (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%ival > en2%ival + en1%tolerance
    else
       y = en1%ival > en2%ival
    end if
  end function comp_gg_ii
  logical function comp_gg_ri (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%rval > en2%ival + en1%tolerance
    else
       y = en1%rval > en2%ival
    end if
  end function comp_gg_ri
  logical function comp_gg_ir (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%ival > en2%rval + en1%tolerance
    else
       y = en1%ival > en2%rval
    end if
  end function comp_gg_ir
  logical function comp_gg_rr (en1, en2) result (y)
    type(eval_node_t), intent(in) :: en1, en2
    if (associated (en1%tolerance)) then
       y = en1%rval > en2%rval + en1%tolerance
    else
       y = en1%rval > en2%rval
    end if
  end function comp_gg_rr

  logical function not_l (en) result (y)
    type(eval_node_t), intent(in) :: en
    y = .not. en%lval
  end function not_l

  subroutine pdg_i (pdg_array, en)
    type(pdg_array_t), intent(out) :: pdg_array
    type(eval_node_t), intent(in) :: en
    pdg_array = en%ival
  end subroutine pdg_i

  subroutine concat_cc (pdg_array, en1, en2)
    type(pdg_array_t), intent(out) :: pdg_array
    type(eval_node_t), intent(in) :: en1, en2
    pdg_array = en1%aval // en2%aval
  end subroutine concat_cc

  subroutine collect_p (subevt, en1, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask1
    integer :: n, i
    n = en1%pval%get_length ()
    allocate (mask1 (n))
    if (present (en0)) then
       do i = 1, n
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          call eval_node_evaluate (en0)
          mask1(i) = en0%lval
       end do
    else
       mask1 = .true.
    end if
    call subevt_collect (subevt, en1%pval, mask1)
  end subroutine collect_p

  subroutine cluster_p (subevt, en1, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask1
    integer :: n, i
    !!! Should not be initialized for every event
    type(jet_definition_t) :: jet_def
    logical :: keep_jets, exclusive
    call jet_def%init (en1%jet_algorithm, en1%jet_r, en1%jet_p, en1%jet_ycut)
    n = en1%pval%get_length ()
    allocate (mask1 (n))
    if (present (en0)) then
       do i = 1, n
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          call eval_node_evaluate (en0)
          mask1(i) = en0%lval
       end do
    else
       mask1 = .true.
    end if
    if (associated (en1%var_list)) then
       keep_jets = en1%var_list%get_lval (var_str("?keep_flavors_when_clustering"))
    else
       keep_jets = .false.
    end if
    exclusive = .false.
    select case (en1%jet_algorithm)
    case (ee_kt_algorithm)
       exclusive = .true.
    case (ee_genkt_algorithm)
       if (en1%jet_r > Pi)  exclusive = .true.
    end select
    call subevt_cluster (subevt, en1%pval, en1%jet_dcut, mask1, &
         jet_def, keep_jets, exclusive)
    call jet_def%final ()
  end subroutine cluster_p

  subroutine photon_recombination_p (subevt, en1, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask1
    type(prt_t), dimension(:), allocatable :: prt
    integer :: n, i
    real(default) :: reco_r0
    logical :: keep_flv
    reco_r0 = en1%photon_rec_r0
    n = en1%pval%get_length ()
    allocate (prt (n))
    do i = 1, n
       prt(i) = en1%pval%get_prt (i)
       if (.not. prt_is_recombinable (prt (i))) then
          call msg_fatal ("Only charged leptons, quarks, and " //&
               "photons can be included in photon recombination.")
       end if
    end do
    if (count (prt_is_photon (prt)) > 1) &
         call msg_fatal ("Photon recombination is supported " // &
            "only for single photons.")
    allocate (mask1 (n))
    if (present (en0)) then
       do i = 1, n
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          call eval_node_evaluate (en0)
          mask1(i) = en0%lval
       end do
    else
       mask1 = .true.
    end if
    if (associated (en1%var_list)) then
       keep_flv = en1%var_list%get_lval &
            (var_str("?keep_flavors_when_recombining"))
    else
       keep_flv = .false.
    end if
    call subevt_recombine &
         (subevt, en1%pval, mask1, reco_r0, keep_flv)
  end subroutine photon_recombination_p

  subroutine select_p (subevt, en1, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask1
    integer :: n, i
    n = en1%pval%get_length ()
    allocate (mask1 (n))
    if (present (en0)) then
       do i = 1, n
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          call eval_node_evaluate (en0)
          mask1(i) = en0%lval
       end do
    else
       mask1 = .true.
    end if
    call subevt_select (subevt, en1%pval, mask1)
  end subroutine select_p

  subroutine select_b_jet_p (subevt, en1, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask1
    integer :: n, i
    n = en1%pval%get_length ()
    allocate (mask1 (n))    
    do i = 1, n
       mask1(i) = prt_is_b_jet (en1%pval%get_prt (i))
       if (present (en0)) then
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          call eval_node_evaluate (en0)
          mask1(i) = en0%lval .and. mask1(i)
        end if
    end do
    call subevt_select (subevt, en1%pval, mask1)
  end subroutine select_b_jet_p

  subroutine select_non_b_jet_p (subevt, en1, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask1
    integer :: n, i
    n = en1%pval%get_length ()
    allocate (mask1 (n))    
    do i = 1, n
       mask1(i) = .not. prt_is_b_jet (en1%pval%get_prt (i))
       if (present (en0)) then
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          call eval_node_evaluate (en0)
          mask1(i) = en0%lval .and. mask1(i)
        end if
    end do
    call subevt_select (subevt, en1%pval, mask1)
  end subroutine select_non_b_jet_p

  subroutine select_c_jet_p (subevt, en1, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask1
    integer :: n, i
    n = en1%pval%get_length ()
    allocate (mask1 (n))    
    do i = 1, n
       mask1(i) = .not. prt_is_b_jet (en1%pval%get_prt (i)) &
            .and. prt_is_c_jet (en1%pval%get_prt (i))
       if (present (en0)) then
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          call eval_node_evaluate (en0)
          mask1(i) = en0%lval .and. mask1(i)
        end if
    end do
    call subevt_select (subevt, en1%pval, mask1)
  end subroutine select_c_jet_p

  subroutine select_light_jet_p (subevt, en1, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask1
    integer :: n, i
    n = en1%pval%get_length ()
    allocate (mask1 (n))    
    do i = 1, n
       mask1(i) = .not. prt_is_b_jet (en1%pval%get_prt (i)) &
            .and. .not. prt_is_c_jet (en1%pval%get_prt (i))
       if (present (en0)) then
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          call eval_node_evaluate (en0)
          mask1(i) = en0%lval .and. mask1(i)
        end if
    end do
    call subevt_select (subevt, en1%pval, mask1)
  end subroutine select_light_jet_p

  subroutine extract_p (subevt, en1, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    integer :: index
    if (present (en0)) then
       call eval_node_evaluate (en0)
       select case (en0%result_type)
       case (V_INT);  index = en0%ival
       case default
          call eval_node_write (en0)
          call msg_fatal (" Index parameter of 'extract' must be integer.")
       end select
    else
       index = 1
    end if
    call subevt_extract (subevt, en1%pval, index)
  end subroutine extract_p

  subroutine sort_p (subevt, en1, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    integer, dimension(:), allocatable :: ival
    real(default), dimension(:), allocatable :: rval
    integer :: i, n
    n = en1%pval%get_length ()
    if (present (en0)) then
       select case (en0%result_type)
       case (V_INT);  allocate (ival (n))
       case (V_REAL); allocate (rval (n))
       end select
       do i = 1, n
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          call eval_node_evaluate (en0)
          select case (en0%result_type)
          case (V_INT);  ival(i) = en0%ival
          case (V_REAL); rval(i) = en0%rval
          end select
       end do
       select case (en0%result_type)
       case (V_INT);  call subevt_sort (subevt, en1%pval, ival)
       case (V_REAL); call subevt_sort (subevt, en1%pval, rval)
       end select
    else
       call subevt_sort (subevt, en1%pval)
    end if
  end subroutine sort_p

  function all_p (en1, en0) result (lval)
    logical :: lval
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout) :: en0
    integer :: i, n
    n = en1%pval%get_length ()
    lval = .true.
    do i = 1, n
       en0%index = i
       en0%prt1 = en1%pval%get_prt (i)
       call eval_node_evaluate (en0)
       lval = en0%lval
       if (.not. lval)  exit
    end do
  end function all_p

  function any_p (en1, en0) result (lval)
    logical :: lval
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout) :: en0
    integer :: i, n
    n = en1%pval%get_length ()
    lval = .false.
    do i = 1, n
       en0%index = i
       en0%prt1 = en1%pval%get_prt (i)
       call eval_node_evaluate (en0)
       lval = en0%lval
       if (lval)  exit
    end do
  end function any_p

  function no_p (en1, en0) result (lval)
    logical :: lval
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout) :: en0
    integer :: i, n
    n = en1%pval%get_length ()
    lval = .true.
    do i = 1, n
       en0%index = i
       en0%prt1 = en1%pval%get_prt (i)
       call eval_node_evaluate (en0)
       lval = .not. en0%lval
       if (lval)  exit
    end do
  end function no_p

  subroutine count_a (ival, en1, en0)
    integer, intent(out) :: ival
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout), optional :: en0
    integer :: i, n, count
    n = en1%pval%get_length ()
    if (present (en0)) then
       count = 0
       do i = 1, n
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          call eval_node_evaluate (en0)
          if (en0%lval)  count = count + 1
       end do
       ival = count
    else
       ival = n
    end if
  end subroutine count_a

  function sum_a (en1, en0) result (rval)
    real(default) :: rval
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout) :: en0
    integer :: i, n
    n = en1%pval%get_length ()
    rval = 0._default
    do i = 1, n
       en0%index = i
       en0%prt1 = en1%pval%get_prt (i)
       call eval_node_evaluate (en0)
       rval = rval + en0%rval
    end do
  end function sum_a

  function prod_a (en1, en0) result (rval)
    real(default) :: rval
    type(eval_node_t), intent(in) :: en1
    type(eval_node_t), intent(inout) :: en0
    integer :: i, n
    n = en1%pval%get_length ()
    rval = 1._default
    do i = 1, n
       en0%index = i
       en0%prt1 = en1%pval%get_prt (i)
       call eval_node_evaluate (en0)
       rval = rval * en0%rval
    end do
  end function prod_a

  subroutine join_pp (subevt, en1, en2, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask2
    integer :: i, j, n1, n2
    n1 = en1%pval%get_length ()
    n2 = en2%pval%get_length ()
    allocate (mask2 (n2))
    mask2 = .true.
    if (present (en0)) then
       do i = 1, n1
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          do j = 1, n2
             en0%prt2 = en2%pval%get_prt (j)
             call eval_node_evaluate (en0)
             mask2(j) = mask2(j) .and. en0%lval
          end do
       end do
    end if
    call subevt_join (subevt, en1%pval, en2%pval, mask2)
  end subroutine join_pp

  subroutine combine_pp (subevt, en1, en2, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:,:), allocatable :: mask12
    integer :: i, j, n1, n2
    n1 = en1%pval%get_length ()
    n2 = en2%pval%get_length ()
    if (present (en0)) then
       allocate (mask12 (n1, n2))
       do i = 1, n1
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          do j = 1, n2
             en0%prt2 = en2%pval%get_prt (j)
             call eval_node_evaluate (en0)
             mask12(i,j) = en0%lval
          end do
       end do
       call subevt_combine (subevt, en1%pval, en2%pval, mask12)
    else
       call subevt_combine (subevt, en1%pval, en2%pval)
    end if
  end subroutine combine_pp

  subroutine collect_pp (subevt, en1, en2, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask1
    integer :: i, j, n1, n2
    n1 = en1%pval%get_length ()
    n2 = en2%pval%get_length ()
    allocate (mask1 (n1))
    mask1 = .true.
    if (present (en0)) then
       do i = 1, n1
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          do j = 1, n2
             en0%prt2 = en2%pval%get_prt (j)
             call eval_node_evaluate (en0)
             mask1(i) = mask1(i) .and. en0%lval
          end do
       end do
    end if
    call subevt_collect (subevt, en1%pval, mask1)
  end subroutine collect_pp

  subroutine select_pp (subevt, en1, en2, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout), optional :: en0
    logical, dimension(:), allocatable :: mask1
    integer :: i, j, n1, n2
    n1 = en1%pval%get_length ()
    n2 = en2%pval%get_length ()
    allocate (mask1 (n1))
    mask1 = .true.
    if (present (en0)) then
       do i = 1, n1
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          do j = 1, n2
             en0%prt2 = en2%pval%get_prt (j)
             call eval_node_evaluate (en0)
             mask1(i) = mask1(i) .and. en0%lval
          end do
       end do
    end if
    call subevt_select (subevt, en1%pval, mask1)
  end subroutine select_pp

  subroutine sort_pp (subevt, en1, en2, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout), optional :: en0
    integer, dimension(:), allocatable :: ival
    real(default), dimension(:), allocatable :: rval
    integer :: i, n1
    n1 = en1%pval%get_length ()
    if (present (en0)) then
       select case (en0%result_type)
       case (V_INT);  allocate (ival (n1))
       case (V_REAL); allocate (rval (n1))
       end select
       do i = 1, n1
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          en0%prt2 = en2%pval%get_prt (1)
          call eval_node_evaluate (en0)
          select case (en0%result_type)
          case (V_INT);  ival(i) = en0%ival
          case (V_REAL); rval(i) = en0%rval
          end select
       end do
       select case (en0%result_type)
       case (V_INT);  call subevt_sort (subevt, en1%pval, ival)
       case (V_REAL); call subevt_sort (subevt, en1%pval, rval)
       end select
    else
       call subevt_sort (subevt, en1%pval)
    end if
  end subroutine sort_pp

  function all_pp (en1, en2, en0) result (lval)
    logical :: lval
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout) :: en0
    integer :: i, j, n1, n2
    n1 = en1%pval%get_length ()
    n2 = en2%pval%get_length ()
    lval = .true.
    LOOP1: do i = 1, n1
       en0%index = i
       en0%prt1 = en1%pval%get_prt (i)
       do j = 1, n2
          en0%prt2 = en2%pval%get_prt (j)
          if (are_disjoint (en0%prt1, en0%prt2)) then
             call eval_node_evaluate (en0)
             lval = en0%lval
             if (.not. lval)  exit LOOP1
          end if
       end do
    end do LOOP1
  end function all_pp

  function any_pp (en1, en2, en0) result (lval)
    logical :: lval
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout) :: en0
    integer :: i, j, n1, n2
    n1 = en1%pval%get_length ()
    n2 = en2%pval%get_length ()
    lval = .false.
    LOOP1: do i = 1, n1
       en0%index = i
       en0%prt1 = en1%pval%get_prt (i)
       do j = 1, n2
          en0%prt2 = en2%pval%get_prt (j)
          if (are_disjoint (en0%prt1, en0%prt2)) then
             call eval_node_evaluate (en0)
             lval = en0%lval
             if (lval)  exit LOOP1
          end if
       end do
    end do LOOP1
  end function any_pp

  function no_pp (en1, en2, en0) result (lval)
    logical :: lval
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout) :: en0
    integer :: i, j, n1, n2
    n1 = en1%pval%get_length ()
    n2 = en2%pval%get_length ()
    lval = .true.
    LOOP1: do i = 1, n1
       en0%index = i
       en0%prt1 = en1%pval%get_prt (i)
       do j = 1, n2
          en0%prt2 = en2%pval%get_prt (j)
          if (are_disjoint (en0%prt1, en0%prt2)) then
             call eval_node_evaluate (en0)
             lval = .not. en0%lval
             if (lval)  exit LOOP1
          end if
       end do
    end do LOOP1
  end function no_pp

  function photon_isolation_pp (en1, en2, en0) result (lval)
    logical :: lval
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout) :: en0
    type(prt_t) :: prt
    type(prt_t), dimension(:), allocatable :: prt_gam0, prt_lep
    type(vector4_t), dimension(:), allocatable :: &
         p_gam0, p_lep0, p_lep, p_par
    integer :: i, j, n1, n2, n_par, n_lep, n_gam, n_delta
    real(default), dimension(:), allocatable :: delta_r, et_sum
    integer, dimension(:), allocatable :: index
    real(default) :: eps, iso_n, r0, pt_gam
    logical, dimension(:,:), allocatable :: photon_mask
    n1 = en1%pval%get_length ()
    n2 = en2%pval%get_length ()
    allocate (p_gam0 (n1), prt_gam0 (n1))
    eps = en1%photon_iso_eps
    iso_n = en1%photon_iso_n
    r0 = en1%photon_iso_r0
    lval = .true.
    do i = 1, n1
       en0%index = i
       prt = en1%pval%get_prt (i)
       prt_gam0(i) = prt
       if (.not. prt_is_photon (prt_gam0(i))) &
            call msg_fatal ("Photon isolation can only " // &
            "be applied to photons.")
       p_gam0(i) = prt_get_momentum (prt_gam0(i))
       en0%prt1 = prt
       call eval_node_evaluate (en0)
       lval = en0%lval
       if (.not. lval) return
    end do
    if (n1 == 0) then
       call msg_fatal ("Photon isolation applied on empty photon sample.")
    end if
    n_par = 0
    n_lep = 0
    n_gam = 0
    do i = 1, n2
       prt = en2%pval%get_prt (i)
       if (prt_is_parton (prt) .or. prt_is_clustered (prt)) then
          n_par = n_par + 1
       end if
       if (prt_is_lepton (prt)) then
          n_lep = n_lep + 1
       end if
       if (prt_is_photon (prt)) then
          n_gam = n_gam + 1
       end if
    end do
    if (n_lep > 0 .and. n_gam == 0) then
       call msg_fatal ("Photon isolation from EM energy: photons " // &
            "have to be included.")
    end if
    if (n_lep > 0 .and. n_gam /= n1) then
       call msg_fatal ("Photon isolation: photon samples do not match.")
    end if
    allocate (p_par (n_par))
    allocate (p_lep0 (n_gam+n_lep), prt_lep(n_gam+n_lep))
    n_par = 0
    n_lep = 0
    do i = 1, n2
       prt = en2%pval%get_prt (i)
       if (prt_is_parton (prt) .or. prt_is_clustered (prt)) then
          n_par = n_par + 1
          p_par(n_par) = prt_get_momentum (prt)
       end if
       if (prt_is_lepton (prt) .or. prt_is_photon(prt)) then
          n_lep = n_lep + 1
          prt_lep(n_lep) = prt
          p_lep0(n_lep) = prt_get_momentum (prt_lep(n_lep))
       end if       
    end do
    if (n_par > 0) then
       allocate (delta_r (n_par), index (n_par))
       HADRON_ISOLATION: do i = 1, n1
          pt_gam = transverse_part (p_gam0(i))
          delta_r(1:n_par) = sort (eta_phi_distance (p_gam0(i), p_par(1:n_par)))
          index(1:n_par) = order (eta_phi_distance (p_gam0(i), p_par(1:n_par)))
          n_delta = count (delta_r < r0)
          allocate (et_sum(n_delta))
          do j = 1, n_delta
             et_sum(j) = sum (transverse_part (p_par (index (1:j))))
             if (.not. et_sum(j) <= &
                  iso_chi_gamma (delta_r(j), r0, iso_n, eps, pt_gam)) then
                lval = .false.
                return
             end if
          end do
          deallocate (et_sum)
       end do HADRON_ISOLATION
       deallocate (delta_r)
       deallocate (index)
    end if
    if (n_lep > 0) then
       allocate (photon_mask(n1,n_lep))
       do i = 1, n1
          photon_mask(i,:) = .not. (prt_gam0(i) .match. prt_lep(:))
       end do
       allocate (delta_r (n_lep-1), index (n_lep-1), p_lep(n_lep-1))
       EM_ISOLATION: do i = 1, n1
          pt_gam = transverse_part (p_gam0(i))
          p_lep = pack (p_lep0, photon_mask(i,:))          
          delta_r(1:n_lep-1) = sort (eta_phi_distance (p_gam0(i), p_lep(1:n_lep-1)))
          index(1:n_lep-1) = order (eta_phi_distance (p_gam0(i), p_lep(1:n_lep-1)))
          n_delta = count (delta_r < r0)
          allocate (et_sum(n_delta))
          do j = 1, n_delta
             et_sum(j) = sum (transverse_part (p_lep (index(1:j))))
             if (.not. et_sum(j) <= &
                  iso_chi_gamma (delta_r(j), r0, iso_n, eps, pt_gam)) then
                lval = .false.
                return
             end if
          end do
          deallocate (et_sum)
       end do EM_ISOLATION
       deallocate (delta_r)
       deallocate (index)
    end if
  contains
    function iso_chi_gamma (dr, r0_gam, n_gam, eps_gam, pt_gam) result (iso)
      real(default) :: iso
      real(default), intent(in) :: dr, r0_gam, n_gam, eps_gam, pt_gam
      iso = eps_gam * pt_gam
      if (.not. nearly_equal (abs(n_gam), 0._default)) then
         iso = iso * ((1._default - cos(dr)) / &
              (1._default - cos(r0_gam)))**abs(n_gam)
      end if
    end function iso_chi_gamma
  end function photon_isolation_pp

  subroutine eval_pp (en1, en2, en0, rval, is_known)
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout) :: en0
    real(default), intent(out) :: rval
    logical, intent(out) :: is_known
    integer :: i, j, n1, n2
    n1 = en1%pval%get_length ()
    n2 = en2%pval%get_length ()
    rval = 0
    is_known = .false.
    LOOP1: do i = 1, n1
       en0%index = i
       en0%prt1 = en1%pval%get_prt (i)
       do j = 1, n2
          en0%prt2 = en2%pval%get_prt (j)
          if (are_disjoint (en0%prt1, en0%prt2)) then
             call eval_node_evaluate (en0)
             rval = en0%rval
             is_known = .true.
             exit LOOP1
          end if
       end do
    end do LOOP1
  end subroutine eval_pp

  subroutine count_pp (ival, en1, en2, en0)
    integer, intent(out) :: ival
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout), optional :: en0
    integer :: i, j, n1, n2, count
    n1 = en1%pval%get_length ()
    n2 = en2%pval%get_length ()
    if (present (en0)) then
       count = 0
       do i = 1, n1
          en0%index = i
          en0%prt1 = en1%pval%get_prt (i)
          do j = 1, n2
             en0%prt2 = en2%pval%get_prt (j)
             if (are_disjoint (en0%prt1, en0%prt2)) then
                call eval_node_evaluate (en0)
                if (en0%lval)  count = count + 1
             end if
          end do
       end do
    else
       count = 0
       do i = 1, n1
          do j = 1, n2
             if (are_disjoint (en1%pval%get_prt (i), &
                               en2%pval%get_prt (j))) then
                count = count + 1
             end if
          end do
       end do
    end if
    ival = count
  end subroutine count_pp

  subroutine select_pdg_ca (subevt, en1, en2, en0)
    type(subevt_t), intent(inout) :: subevt
    type(eval_node_t), intent(in) :: en1, en2
    type(eval_node_t), intent(inout), optional :: en0
    if (present (en0)) then
       call subevt_select_pdg_code (subevt, en1%aval, en2%pval, en0%ival)
    else
       call subevt_select_pdg_code (subevt, en1%aval, en2%pval)
    end if
  end subroutine select_pdg_ca

  subroutine concat_ss (string, en1, en2)
    type(string_t), intent(out) :: string
    type(eval_node_t), intent(in) :: en1, en2
    string = en1%sval // en2%sval
  end subroutine concat_ss

  recursive subroutine eval_node_compile_genexpr &
       (en, pn, var_list, result_type)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    integer, intent(in), optional :: result_type
    if (debug_active (D_MODEL_F)) then
       print *, "read genexpr";  call parse_node_write (pn)
    end if
    if (present (result_type)) then
       select case (result_type)
       case (V_INT, V_REAL, V_CMPLX)
          call eval_node_compile_expr  (en, pn, var_list)
       case (V_LOG)
          call eval_node_compile_lexpr (en, pn, var_list)
       case (V_SEV)
          call eval_node_compile_pexpr (en, pn, var_list)
       case (V_PDG)
          call eval_node_compile_cexpr (en, pn, var_list)
       case (V_STR)
          call eval_node_compile_sexpr (en, pn, var_list)
       end select
    else
       call eval_node_compile_expr  (en, pn, var_list)
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done genexpr"
    end if
  end subroutine eval_node_compile_genexpr

  recursive subroutine eval_node_compile_expr (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_term, pn_addition, pn_op, pn_arg
    type(eval_node_t), pointer :: en1, en2
    type(string_t) :: key
    integer :: t1, t2, t
    if (debug_active (D_MODEL_F)) then
       print *, "read expr";  call parse_node_write (pn)
    end if
    pn_term => parse_node_get_sub_ptr (pn)
    select case (char (parse_node_get_rule_key (pn_term)))
    case ("term")
       call eval_node_compile_term (en, pn_term, var_list)
       pn_addition => parse_node_get_next_ptr (pn_term, tag="addition")
    case ("addition")
       en => null ()
       pn_addition => pn_term
    case default
       call parse_node_mismatch ("term|addition", pn)
    end select
    do while (associated (pn_addition))
       pn_op => parse_node_get_sub_ptr (pn_addition)
       pn_arg => parse_node_get_next_ptr (pn_op, tag="term")
       call eval_node_compile_term (en2, pn_arg, var_list)
       t2 = en2%result_type
       if (associated (en)) then
          en1 => en
          t1 = en1%result_type
       else
          allocate (en1)
          select case (t2)
          case (V_INT);  call eval_node_init_int  (en1, 0)
          case (V_REAL); call eval_node_init_real (en1, 0._default)
          case (V_CMPLX); call eval_node_init_cmplx (en1, cmplx &
                         (0._default, 0._default, kind=default))
          end select
          t1 = t2
       end if
       t = numeric_result_type (t1, t2)
       allocate (en)
       key = parse_node_get_key (pn_op)
       if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
          select case (char (key))
          case ("+")
             select case (t1)
             case (V_INT)
                select case (t2)
                case (V_INT);  call eval_node_init_int  (en, add_ii (en1, en2))
                case (V_REAL); call eval_node_init_real (en, add_ir (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, add_ic (en1, en2))
                end select
             case (V_REAL)
                select case (t2)
                case (V_INT);  call eval_node_init_real (en, add_ri (en1, en2))
                case (V_REAL); call eval_node_init_real (en, add_rr (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, add_rc (en1, en2))
                end select
              case (V_CMPLX)
                select case (t2)
                case (V_INT);  call eval_node_init_cmplx (en, add_ci (en1, en2))
                case (V_REAL); call eval_node_init_cmplx (en, add_cr (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, add_cc (en1, en2))
                end select
             end select
          case ("-")
             select case (t1)
             case (V_INT)
                select case (t2)
                case (V_INT);  call eval_node_init_int  (en, sub_ii (en1, en2))
                case (V_REAL); call eval_node_init_real (en, sub_ir (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, sub_ic (en1, en2))
                end select
             case (V_REAL)
                select case (t2)
                case (V_INT);  call eval_node_init_real (en, sub_ri (en1, en2))
                case (V_REAL); call eval_node_init_real (en, sub_rr (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, sub_rc (en1, en2))
                end select
             case (V_CMPLX)
                select case (t2)
                case (V_INT);  call eval_node_init_cmplx (en, sub_ci (en1, en2))
                case (V_REAL); call eval_node_init_cmplx (en, sub_cr (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, sub_cc (en1, en2))
                end select
             end select
          end select
          call eval_node_final_rec (en1)
          call eval_node_final_rec (en2)
          deallocate (en1, en2)
       else
          call eval_node_init_branch (en, key, t, en1, en2)
          select case (char (key))
          case ("+")
             select case (t1)
             case (V_INT)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_int  (en, add_ii)
                case (V_REAL);  call eval_node_set_op2_real (en, add_ir)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, add_ic)
                end select
             case (V_REAL)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_real (en, add_ri)
                case (V_REAL);  call eval_node_set_op2_real (en, add_rr)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, add_rc)
                end select
             case (V_CMPLX)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_cmplx (en, add_ci)
                case (V_REAL);  call eval_node_set_op2_cmplx (en, add_cr)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, add_cc)
                end select
             end select
          case ("-")
             select case (t1)
             case (V_INT)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_int  (en, sub_ii)
                case (V_REAL);  call eval_node_set_op2_real (en, sub_ir)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, sub_ic)
                end select
             case (V_REAL)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_real (en, sub_ri)
                case (V_REAL);  call eval_node_set_op2_real (en, sub_rr)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, sub_rc)
                end select
             case (V_CMPLX)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_cmplx (en, sub_ci)
                case (V_REAL);  call eval_node_set_op2_cmplx (en, sub_cr)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, sub_cc)
                end select
             end select
          end select
       end if
       pn_addition => parse_node_get_next_ptr (pn_addition)
    end do
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done expr"
    end if
  end subroutine eval_node_compile_expr

  recursive subroutine eval_node_compile_term (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_factor, pn_multiplication, pn_op, pn_arg
    type(eval_node_t), pointer :: en1, en2
    type(string_t) :: key
    integer :: t1, t2, t
    if (debug_active (D_MODEL_F)) then
       print *, "read term";  call parse_node_write (pn)
    end if
    pn_factor => parse_node_get_sub_ptr (pn, tag="factor")
    call eval_node_compile_factor (en, pn_factor, var_list)
    pn_multiplication => &
         parse_node_get_next_ptr (pn_factor, tag="multiplication")
    do while (associated (pn_multiplication))
       pn_op => parse_node_get_sub_ptr (pn_multiplication)
       pn_arg => parse_node_get_next_ptr (pn_op, tag="factor")
       en1 => en
       call eval_node_compile_factor (en2, pn_arg, var_list)
       t1 = en1%result_type
       t2 = en2%result_type
       t = numeric_result_type (t1, t2)
       allocate (en)
       key = parse_node_get_key (pn_op)
       if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
          select case (char (key))
          case ("*")
             select case (t1)
             case (V_INT)
                select case (t2)
                case (V_INT);  call eval_node_init_int  (en, mul_ii (en1, en2))
                case (V_REAL); call eval_node_init_real (en, mul_ir (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, mul_ic (en1, en2))
                end select
             case (V_REAL)
                select case (t2)
                case (V_INT);  call eval_node_init_real (en, mul_ri (en1, en2))
                case (V_REAL); call eval_node_init_real (en, mul_rr (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, mul_rc (en1, en2))
                end select
              case (V_CMPLX)
                select case (t2)
                case (V_INT);  call eval_node_init_cmplx (en, mul_ci (en1, en2))
                case (V_REAL); call eval_node_init_cmplx (en, mul_cr (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, mul_cc (en1, en2))
                end select
             end select
          case ("/")
             select case (t1)
             case (V_INT)
                select case (t2)
                case (V_INT);  call eval_node_init_int  (en, div_ii (en1, en2))
                case (V_REAL); call eval_node_init_real (en, div_ir (en1, en2))
                case (V_CMPLX); call eval_node_init_real (en, div_ir (en1, en2))
                end select
             case (V_REAL)
                select case (t2)
                case (V_INT);  call eval_node_init_real (en, div_ri (en1, en2))
                case (V_REAL); call eval_node_init_real (en, div_rr (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, div_rc (en1, en2))
                end select
             case (V_CMPLX)
                select case (t2)
                case (V_INT);  call eval_node_init_cmplx (en, div_ci (en1, en2))
                case (V_REAL); call eval_node_init_cmplx (en, div_cr (en1, en2))
                case (V_CMPLX); call eval_node_init_cmplx (en, div_cc (en1, en2))
                end select
             end select
          end select
          call eval_node_final_rec (en1)
          call eval_node_final_rec (en2)
          deallocate (en1, en2)
       else
          call eval_node_init_branch (en, key, t, en1, en2)
          select case (char (key))
          case ("*")
             select case (t1)
             case (V_INT)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_int  (en, mul_ii)
                case (V_REAL);  call eval_node_set_op2_real (en, mul_ir)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, mul_ic)
                end select
             case (V_REAL)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_real (en, mul_ri)
                case (V_REAL);  call eval_node_set_op2_real (en, mul_rr)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, mul_rc)
                end select
             case (V_CMPLX)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_cmplx (en, mul_ci)
                case (V_REAL);  call eval_node_set_op2_cmplx (en, mul_cr)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, mul_cc)
                end select
             end select
          case ("/")
             select case (t1)
             case (V_INT)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_int  (en, div_ii)
                case (V_REAL);  call eval_node_set_op2_real (en, div_ir)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, div_ic)
                end select
             case (V_REAL)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_real (en, div_ri)
                case (V_REAL);  call eval_node_set_op2_real (en, div_rr)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, div_rc)
                end select
             case (V_CMPLX)
                select case (t2)
                case (V_INT);   call eval_node_set_op2_cmplx (en, div_ci)
                case (V_REAL);  call eval_node_set_op2_cmplx (en, div_cr)
                case (V_CMPLX);  call eval_node_set_op2_cmplx (en, div_cc)
                end select
             end select
          end select
       end if
       pn_multiplication => parse_node_get_next_ptr (pn_multiplication)
    end do
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done term"
    end if
  end subroutine eval_node_compile_term

  recursive subroutine eval_node_compile_factor (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_value, pn_exponentiation, pn_op, pn_arg
    type(eval_node_t), pointer :: en1, en2
    type(string_t) :: key
    integer :: t1, t2, t
    if (debug_active (D_MODEL_F)) then
       print *, "read factor";  call parse_node_write (pn)
    end if
    pn_value => parse_node_get_sub_ptr (pn)
    call eval_node_compile_signed_value (en, pn_value, var_list)
    pn_exponentiation => &
         parse_node_get_next_ptr (pn_value, tag="exponentiation")
    if (associated (pn_exponentiation)) then
       pn_op => parse_node_get_sub_ptr (pn_exponentiation)
       pn_arg => parse_node_get_next_ptr (pn_op)
       en1 => en
       call eval_node_compile_signed_value (en2, pn_arg, var_list)
       t1 = en1%result_type
       t2 = en2%result_type
       t = numeric_result_type (t1, t2)
       allocate (en)
       key = parse_node_get_key (pn_op)
       if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);   call eval_node_init_int   (en, pow_ii (en1, en2))
             case (V_REAL);  call eval_node_init_real  (en, pow_ir (en1, en2))
             case (V_CMPLX); call eval_node_init_cmplx (en, pow_ic (en1, en2))
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);   call eval_node_init_real  (en, pow_ri (en1, en2))
             case (V_REAL);  call eval_node_init_real  (en, pow_rr (en1, en2))
             case (V_CMPLX); call eval_node_init_cmplx (en, pow_rc (en1, en2))
             end select
          case (V_CMPLX)
             select case (t2)
             case (V_INT);   call eval_node_init_cmplx (en, pow_ci (en1, en2))
             case (V_REAL);  call eval_node_init_cmplx (en, pow_cr (en1, en2))
             case (V_CMPLX); call eval_node_init_cmplx (en, pow_cc (en1, en2))
             end select
          end select
          call eval_node_final_rec (en1)
          call eval_node_final_rec (en2)
          deallocate (en1, en2)
       else
          call eval_node_init_branch (en, key, t, en1, en2)
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);   call eval_node_set_op2_int  (en, pow_ii)
             case (V_REAL,V_CMPLX);  call eval_type_error (pn, "exponentiation", t1)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);   call eval_node_set_op2_real (en, pow_ri)
             case (V_REAL);  call eval_node_set_op2_real (en, pow_rr)
             case (V_CMPLX);  call eval_type_error (pn, "exponentiation", t1)
             end select
          case (V_CMPLX)
             select case (t2)
             case (V_INT);   call eval_node_set_op2_cmplx (en, pow_ci)
             case (V_REAL);  call eval_node_set_op2_cmplx (en, pow_cr)
             case (V_CMPLX);  call eval_node_set_op2_cmplx (en, pow_cc)
             end select
          end select
       end if
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done factor"
    end if
  end subroutine eval_node_compile_factor

  recursive subroutine eval_node_compile_signed_value (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_arg
    type(eval_node_t), pointer :: en1
    integer :: t
    if (debug_active (D_MODEL_F)) then
       print *, "read signed value";  call parse_node_write (pn)
    end if
    select case (char (parse_node_get_rule_key (pn)))
    case ("signed_value")
       pn_arg => parse_node_get_sub_ptr (pn, 2)
       call eval_node_compile_value (en1, pn_arg, var_list)
       t = en1%result_type
       allocate (en)
       if (en1%type == EN_CONSTANT) then
          select case (t)
          case (V_INT);  call eval_node_init_int  (en, neg_i (en1))
          case (V_REAL); call eval_node_init_real (en, neg_r (en1))
          case (V_CMPLX); call eval_node_init_cmplx (en, neg_c (en1))
          end select
          call eval_node_final_rec (en1)
          deallocate (en1)
       else
          call eval_node_init_branch (en, var_str ("-"), t, en1)
          select case (t)
          case (V_INT);  call eval_node_set_op1_int  (en, neg_i)
          case (V_REAL); call eval_node_set_op1_real (en, neg_r)
          case (V_CMPLX); call eval_node_set_op1_cmplx (en, neg_c)
          end select
       end if
    case default
       call eval_node_compile_value (en, pn, var_list)
    end select
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done signed value"
    end if
  end subroutine eval_node_compile_signed_value

  recursive subroutine eval_node_compile_value (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    if (debug_active (D_MODEL_F)) then
       print *, "read value";  call parse_node_write (pn)
    end if
    select case (char (parse_node_get_rule_key (pn)))
    case ("integer_value", "real_value", "complex_value")
       call eval_node_compile_numeric_value (en, pn)
    case ("pi")
       call eval_node_compile_constant (en, pn)
    case ("I")
       call eval_node_compile_constant (en, pn)
    case ("variable")
       call eval_node_compile_variable (en, pn, var_list)
    case ("result")
       call eval_node_compile_result (en, pn, var_list)
    case ("expr")
       call eval_node_compile_expr (en, pn, var_list)
    case ("block_expr")
       call eval_node_compile_block_expr (en, pn, var_list)
    case ("conditional_expr")
       call eval_node_compile_conditional (en, pn, var_list)
    case ("unary_function")
       call eval_node_compile_unary_function (en, pn, var_list)
    case ("binary_function")
       call eval_node_compile_binary_function (en, pn, var_list)
    case ("eval_fun")
       call eval_node_compile_eval_function (en, pn, var_list)
    case ("count_fun")
       call eval_node_compile_count_function (en, pn, var_list)
    case ("sum_fun", "prod_fun")
       call eval_node_compile_numeric_function (en, pn, var_list)
    case default
       call parse_node_mismatch &
            ("integer|real|complex|constant|variable|" // &
             "expr|block_expr|conditional_expr|" // &
             "unary_function|binary_function|numeric_pexpr", pn)
    end select
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done value"
    end if
  end subroutine eval_node_compile_value

  subroutine eval_node_compile_numeric_value (en, pn)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in), target :: pn
    type(parse_node_t), pointer :: pn_val, pn_unit
    allocate (en)
    pn_val => parse_node_get_sub_ptr (pn)
    pn_unit => parse_node_get_next_ptr (pn_val)
    select case (char (parse_node_get_rule_key (pn)))
    case ("integer_value")
       if (associated (pn_unit)) then
          call eval_node_init_real (en, &
               parse_node_get_integer (pn_val) * parse_node_get_unit (pn_unit))
       else
          call eval_node_init_int (en, parse_node_get_integer (pn_val))
       end if
    case ("real_value")
       if (associated (pn_unit)) then
          call eval_node_init_real (en, &
               parse_node_get_real (pn_val) * parse_node_get_unit (pn_unit))
       else
          call eval_node_init_real (en, parse_node_get_real (pn_val))
       end if
    case ("complex_value")
       if (associated (pn_unit)) then
          call eval_node_init_cmplx (en, &
               parse_node_get_cmplx (pn_val) * parse_node_get_unit (pn_unit))
       else
          call eval_node_init_cmplx (en, parse_node_get_cmplx (pn_val))
       end if
    case ("neg_real_value")
       pn_val => parse_node_get_sub_ptr (parse_node_get_sub_ptr (pn, 2))
       pn_unit => parse_node_get_next_ptr (pn_val)
       if (associated (pn_unit)) then
          call eval_node_init_real (en, &
               - parse_node_get_real (pn_val) * parse_node_get_unit (pn_unit))
       else
          call eval_node_init_real (en, - parse_node_get_real (pn_val))
       end if
    case ("pos_real_value")
       pn_val => parse_node_get_sub_ptr (parse_node_get_sub_ptr (pn, 2))
       pn_unit => parse_node_get_next_ptr (pn_val)
       if (associated (pn_unit)) then
          call eval_node_init_real (en, &
               parse_node_get_real (pn_val) * parse_node_get_unit (pn_unit))
       else
          call eval_node_init_real (en, parse_node_get_real (pn_val))
       end if
    case default
       call parse_node_mismatch &
       ("integer_value|real_value|complex_value|neg_real_value|pos_real_value", pn)
    end select
  end subroutine eval_node_compile_numeric_value

  function parse_node_get_unit (pn) result (factor)
    real(default) :: factor
    real(default) :: unit
    type(parse_node_t), intent(in) :: pn
    type(parse_node_t), pointer :: pn_unit, pn_unit_power
    type(parse_node_t), pointer :: pn_frac, pn_num, pn_int, pn_div, pn_den
    integer :: num, den
    pn_unit => parse_node_get_sub_ptr (pn)
    select case (char (parse_node_get_key (pn_unit)))
    case ("TeV");  unit = 1.e3_default
    case ("GeV");  unit = 1
    case ("MeV");  unit = 1.e-3_default
    case ("keV");  unit = 1.e-6_default
    case ("eV");   unit = 1.e-9_default
    case ("meV");   unit = 1.e-12_default
    case ("nbarn");  unit = 1.e6_default
    case ("pbarn");  unit = 1.e3_default
    case ("fbarn");  unit = 1
    case ("abarn");  unit = 1.e-3_default
    case ("rad");     unit = 1
    case ("mrad");    unit = 1.e-3_default
    case ("degree");  unit = degree
    case ("%");  unit = 1.e-2_default
    case default
       call msg_bug (" Unit '" // &
            char (parse_node_get_key (pn)) // "' is undefined.")
    end select
    pn_unit_power => parse_node_get_next_ptr (pn_unit)
    if (associated (pn_unit_power)) then
       pn_frac => parse_node_get_sub_ptr (pn_unit_power, 2)
       pn_num => parse_node_get_sub_ptr (pn_frac)
       select case (char (parse_node_get_rule_key (pn_num)))
       case ("neg_int")
          pn_int => parse_node_get_sub_ptr (pn_num, 2)
          num = - parse_node_get_integer (pn_int)
       case ("pos_int")
          pn_int => parse_node_get_sub_ptr (pn_num, 2)
          num = parse_node_get_integer (pn_int)
       case ("integer_literal")
          num = parse_node_get_integer (pn_num)
       case default
          call parse_node_mismatch ("neg_int|pos_int|integer_literal", pn_num)
       end select
       pn_div => parse_node_get_next_ptr (pn_num)
       if (associated (pn_div)) then
          pn_den => parse_node_get_sub_ptr (pn_div, 2)
          den = parse_node_get_integer (pn_den)
       else
          den = 1
       end if
    else
       num = 1
       den = 1
    end if
    factor = unit ** (real (num, default) / den)
  end function parse_node_get_unit

  subroutine eval_node_compile_constant (en, pn)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    if (debug_active (D_MODEL_F)) then
       print *, "read constant";  call parse_node_write (pn)
    end if
    allocate (en)
    select case (char (parse_node_get_key (pn)))
    case ("pi");     call eval_node_init_real (en, pi)
    case ("I");      call eval_node_init_cmplx (en, imago)
    case default
       call parse_node_mismatch ("pi or I", pn)
    end select
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done constant"
    end if
  end subroutine eval_node_compile_constant

  recursive subroutine eval_node_compile_variable (en, pn, var_list, var_type)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in), target :: pn
    type(var_list_t), intent(in), target :: var_list
    integer, intent(in), optional :: var_type
    type(parse_node_t), pointer :: pn_name
    type(string_t) :: var_name
    logical, target, save :: no_lval
    real(default), target, save :: no_rval
    type(subevt_t), target, save :: no_pval
    type(string_t), target, save :: no_sval
    logical, target, save :: unknown = .false.
    integer :: type
    logical :: defined
    logical, pointer :: known
    logical, pointer :: lptr
    integer, pointer :: iptr
    real(default), pointer :: rptr
    complex(default), pointer :: cptr
    type(subevt_t), pointer :: pptr
    type(string_t), pointer :: sptr
    procedure(obs_unary_int), pointer :: obs1_iptr
    procedure(obs_unary_real), pointer :: obs1_rptr
    procedure(obs_binary_int), pointer :: obs2_iptr
    procedure(obs_binary_real), pointer :: obs2_rptr
    procedure(obs_sev_int), pointer :: obsev_iptr
    procedure(obs_sev_real), pointer :: obsev_rptr
    type(prt_t), pointer :: p1, p2
    if (debug_active (D_MODEL_F)) then
       print *, "read variable";  call parse_node_write (pn)
    end if
    if (present (var_type)) then
       select case (var_type)
       case (V_REAL, V_OBS1_REAL, V_OBS2_REAL, V_INT, V_OBS1_INT, &
                V_OBS2_INT, V_CMPLX)
          pn_name => pn
       case default
          pn_name => parse_node_get_sub_ptr (pn, 2)
       end select
    else
       pn_name => pn
    end if
    select case (char (parse_node_get_rule_key (pn_name)))
    case ("expr")
       call eval_node_compile_expr (en, pn_name, var_list)
    case ("lexpr")
       call eval_node_compile_lexpr (en, pn_name, var_list)
    case ("sexpr")
       call eval_node_compile_sexpr (en, pn_name, var_list)
    case ("pexpr")
       call eval_node_compile_pexpr (en, pn_name, var_list)
    case ("variable")
       var_name = parse_node_get_string (pn_name)
       if (present (var_type)) then
          select case (var_type)
          case (V_LOG);  var_name = "?" // var_name
          case (V_SEV);  var_name = "@" // var_name
          case (V_STR);  var_name = "$" // var_name   ! $ sign
          end select
       end if
       call var_list%get_var_properties &
            (var_name, req_type=var_type, type=type, is_defined=defined)
       allocate (en)
       if (defined) then
          select case (type)
          case (V_LOG)
             call var_list%get_lptr (var_name, lptr, known)
             call eval_node_init_log_ptr (en, var_name, lptr, known)
          case (V_INT)
             call var_list%get_iptr (var_name, iptr, known)
             call eval_node_init_int_ptr (en, var_name, iptr, known)
          case (V_REAL)
             call var_list%get_rptr (var_name, rptr, known)
             call eval_node_init_real_ptr (en, var_name, rptr, known)
          case (V_CMPLX)
             call var_list%get_cptr (var_name, cptr, known)
             call eval_node_init_cmplx_ptr (en, var_name, cptr, known)
          case (V_SEV)
             call var_list%get_pptr (var_name, pptr, known)
             call eval_node_init_subevt_ptr (en, var_name, pptr, known)
          case (V_STR)
             call var_list%get_sptr (var_name, sptr, known)
             call eval_node_init_string_ptr (en, var_name, sptr, known)
          case (V_OBS1_INT)
             call var_list%get_obs1_iptr (var_name, obs1_iptr, p1)
             call eval_node_init_obs1_int_ptr (en, var_name, obs1_iptr, p1)
          case (V_OBS2_INT)
             call var_list%get_obs2_iptr (var_name, obs2_iptr, p1, p2)
             call eval_node_init_obs2_int_ptr (en, var_name, obs2_iptr, p1, p2)
          case (V_OBSEV_INT)
             call var_list%get_obsev_iptr (var_name, obsev_iptr, pptr)
             call eval_node_init_obsev_int_ptr (en, var_name, obsev_iptr, pptr)
          case (V_OBS1_REAL)
             call var_list%get_obs1_rptr (var_name, obs1_rptr, p1)
             call eval_node_init_obs1_real_ptr (en, var_name, obs1_rptr, p1)
          case (V_OBS2_REAL)
             call var_list%get_obs2_rptr (var_name, obs2_rptr, p1, p2)
             call eval_node_init_obs2_real_ptr (en, var_name, obs2_rptr, p1, p2)
          case (V_OBSEV_REAL)
             call var_list%get_obsev_rptr (var_name, obsev_rptr, pptr)
             call eval_node_init_obsev_real_ptr (en, var_name, obsev_rptr, pptr)
          case default
             call parse_node_write (pn)
             call msg_fatal ("Variable of this type " // &
                  "is not allowed in the present context")
             if (present (var_type)) then
                select case (var_type)
                case (V_LOG)
                   call eval_node_init_log_ptr (en, var_name, no_lval, unknown)
                case (V_SEV)
                   call eval_node_init_subevt_ptr &
                        (en, var_name, no_pval, unknown)
                case (V_STR)
                   call eval_node_init_string_ptr &
                        (en, var_name, no_sval, unknown)
                end select
             else
                call eval_node_init_real_ptr (en, var_name, no_rval, unknown)
             end if
          end select
       else
          call parse_node_write (pn)
          call msg_error ("This variable is undefined at this point")
          if (present (var_type)) then
             select case (var_type)
             case (V_LOG)
                call eval_node_init_log_ptr (en, var_name, no_lval, unknown)
             case (V_SEV)
                call eval_node_init_subevt_ptr &
                     (en, var_name, no_pval, unknown)
             case (V_STR)
                call eval_node_init_string_ptr (en, var_name, no_sval, unknown)
             end select
          else
             call eval_node_init_real_ptr (en, var_name, no_rval, unknown)
          end if
       end if
    end select
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done variable"
    end if
  end subroutine eval_node_compile_variable

  subroutine check_var_type (pn, ok, type_actual, type_requested)
    type(parse_node_t), intent(in) :: pn
    logical, intent(out) :: ok
    integer, intent(in) :: type_actual
    integer, intent(in), optional :: type_requested
    if (present (type_requested)) then
       select case (type_requested)
       case (V_LOG)
          select case (type_actual)
          case (V_LOG)
          case default
             call parse_node_write (pn)
             call msg_fatal ("Variable type is invalid (should be logical)")
             ok = .false.
          end select
       case (V_SEV)
          select case (type_actual)
          case (V_SEV)
          case default
             call parse_node_write (pn)
             call msg_fatal &
                  ("Variable type is invalid (should be particle set)")
             ok = .false.
          end select
       case (V_PDG)
          select case (type_actual)
          case (V_PDG)
          case default
             call parse_node_write (pn)
             call msg_fatal &
                  ("Variable type is invalid (should be PDG array)")
             ok = .false.
          end select
       case (V_STR)
          select case (type_actual)
          case (V_STR)
          case default
             call parse_node_write (pn)
             call msg_fatal &
                  ("Variable type is invalid (should be string)")
             ok = .false.
          end select
       case default
          call parse_node_write (pn)
          call msg_bug ("Variable type is unknown")
       end select
    else
       select case (type_actual)
       case (V_REAL, V_OBS1_REAL, V_OBS2_REAL, V_INT, V_OBS1_INT, &
                V_OBS2_INT, V_CMPLX)
       case default
          call parse_node_write (pn)
          call msg_fatal ("Variable type is invalid (should be numeric)")
          ok = .false.
       end select
    end if
    ok = .true.
  end subroutine check_var_type

  subroutine eval_node_compile_result (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in), target :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_key, pn_prc_id
    type(string_t) :: key, prc_id, var_name
    integer, pointer :: iptr
    real(default), pointer :: rptr
    logical, pointer :: known
    if (debug_active (D_MODEL_F)) then
       print *, "read result";  call parse_node_write (pn)
    end if
    pn_key => parse_node_get_sub_ptr (pn)
    pn_prc_id => parse_node_get_next_ptr (pn_key)
    key = parse_node_get_key (pn_key)
    prc_id = parse_node_get_string (pn_prc_id)
    var_name = key // "(" // prc_id // ")"
    if (var_list%contains (var_name)) then
       allocate (en)
       select case (char(key))
       case ("num_id", "n_calls")
          call var_list%get_iptr (var_name, iptr, known)
          call eval_node_init_int_ptr (en, var_name, iptr, known)
       case ("integral", "error")
          call var_list%get_rptr (var_name, rptr, known)
          call eval_node_init_real_ptr (en, var_name, rptr, known)
       end select
    else
       call msg_fatal ("Result variable '" // char (var_name) &
            // "' is undefined (call 'integrate' before use)")
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done result"
    end if
  end subroutine eval_node_compile_result

  recursive subroutine eval_node_compile_unary_function (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_fname, pn_arg
    type(eval_node_t), pointer :: en1
    type(string_t) :: key
    integer :: t
    if (debug_active (D_MODEL_F)) then
       print *, "read unary function";  call parse_node_write (pn)
    end if
    pn_fname => parse_node_get_sub_ptr (pn)
    pn_arg => parse_node_get_next_ptr (pn_fname, tag="function_arg1")
    call eval_node_compile_expr &
         (en1, parse_node_get_sub_ptr (pn_arg, tag="expr"), var_list)
    t = en1%result_type
    allocate (en)
    key = parse_node_get_key (pn_fname)
    if (en1%type == EN_CONSTANT) then
       select case (char (key))
       case ("complex")
          select case (t)
          case (V_INT);  call eval_node_init_cmplx (en, cmplx_i (en1))
          case (V_REAL); call eval_node_init_cmplx (en, cmplx_r (en1))
          case (V_CMPLX); deallocate (en);  en => en1;  en1 => null ()
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("real")
          select case (t)
          case (V_INT);  call eval_node_init_real (en, real_i (en1))
          case (V_REAL); deallocate (en);  en => en1;  en1 => null ()
          case (V_CMPLX); call eval_node_init_real (en, real_c (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("int")
          select case (t)
          case (V_INT);  deallocate (en);  en => en1;  en1 => null ()
          case (V_REAL); call eval_node_init_int  (en, int_r (en1))
          case (V_CMPLX); call eval_node_init_int  (en, int_c (en1))
          end select
       case ("nint")
          select case (t)
          case (V_INT);  deallocate (en);  en => en1;  en1 => null ()
          case (V_REAL); call eval_node_init_int  (en, nint_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("floor")
          select case (t)
          case (V_INT);  deallocate (en);  en => en1;  en1 => null ()
          case (V_REAL); call eval_node_init_int  (en, floor_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("ceiling")
          select case (t)
          case (V_INT);  deallocate (en);  en => en1;  en1 => null ()
          case (V_REAL); call eval_node_init_int  (en, ceiling_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("abs")
          select case (t)
          case (V_INT);  call eval_node_init_int  (en, abs_i (en1))
          case (V_REAL); call eval_node_init_real (en, abs_r (en1))
          case (V_CMPLX); call eval_node_init_real (en, abs_c (en1))
          end select
       case ("conjg")
          select case (t)
          case (V_INT);  call eval_node_init_int  (en, conjg_i (en1))
          case (V_REAL); call eval_node_init_real (en, conjg_r (en1))
          case (V_CMPLX); call eval_node_init_cmplx (en, conjg_c (en1))
          end select
       case ("sgn")
          select case (t)
          case (V_INT);  call eval_node_init_int  (en, sgn_i (en1))
          case (V_REAL); call eval_node_init_real (en, sgn_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("sqrt")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, sqrt_r (en1))
          case (V_CMPLX); call eval_node_init_cmplx (en, sqrt_c (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("exp")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, exp_r (en1))
          case (V_CMPLX); call eval_node_init_cmplx (en, exp_c (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("log")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, log_r (en1))
          case (V_CMPLX); call eval_node_init_cmplx (en, log_c (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("log10")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, log10_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("sin")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, sin_r (en1))
          case (V_CMPLX); call eval_node_init_cmplx (en, sin_c (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("cos")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, cos_r (en1))
          case (V_CMPLX); call eval_node_init_cmplx (en, cos_c (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("tan")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, tan_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("asin")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, asin_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("acos")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, acos_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("atan")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, atan_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("sinh")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, sinh_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("cosh")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, cosh_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("tanh")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, tanh_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("asinh")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, asinh_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("acosh")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, acosh_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("atanh")
          select case (t)
          case (V_REAL); call eval_node_init_real (en, atanh_r (en1))
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case default
          call parse_node_mismatch ("function name", pn_fname)
       end select
       if (associated (en1)) then
          call eval_node_final_rec (en1)
          deallocate (en1)
       end if
    else
       select case (char (key))
       case ("complex")
          call eval_node_init_branch (en, key, V_CMPLX, en1)
       case ("real")
          call eval_node_init_branch (en, key, V_REAL, en1)
       case ("int", "nint", "floor", "ceiling")
          call eval_node_init_branch (en, key, V_INT, en1)
       case default
          call eval_node_init_branch (en, key, t, en1)
       end select
       select case (char (key))
       case ("complex")
          select case (t)
          case (V_INT);  call eval_node_set_op1_cmplx (en, cmplx_i)
          case (V_REAL); call eval_node_set_op1_cmplx (en, cmplx_r)
          case (V_CMPLX); deallocate (en);  en => en1
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("real")
          select case (t)
          case (V_INT);  call eval_node_set_op1_real (en, real_i)
          case (V_REAL); deallocate (en);  en => en1
          case (V_CMPLX); call eval_node_set_op1_real (en, real_c)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("int")
          select case (t)
          case (V_INT);  deallocate (en);  en => en1
          case (V_REAL); call eval_node_set_op1_int (en, int_r)
          case (V_CMPLX); call eval_node_set_op1_int (en, int_c)
          end select
       case ("nint")
          select case (t)
          case (V_INT);  deallocate (en);  en => en1
          case (V_REAL); call eval_node_set_op1_int (en, nint_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("floor")
          select case (t)
          case (V_INT);  deallocate (en);  en => en1
          case (V_REAL); call eval_node_set_op1_int (en, floor_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("ceiling")
          select case (t)
          case (V_INT);  deallocate (en);  en => en1
          case (V_REAL); call eval_node_set_op1_int (en, ceiling_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("abs")
          select case (t)
          case (V_INT);  call eval_node_set_op1_int  (en, abs_i)
          case (V_REAL); call eval_node_set_op1_real (en, abs_r)
          case (V_CMPLX);
             call eval_node_init_branch (en, key, V_REAL, en1)
             call eval_node_set_op1_real (en, abs_c)
          end select
       case ("conjg")
          select case (t)
          case (V_INT);  call eval_node_set_op1_int  (en, conjg_i)
          case (V_REAL); call eval_node_set_op1_real (en, conjg_r)
          case (V_CMPLX);  call eval_node_set_op1_cmplx (en, conjg_c)
          end select
       case ("sgn")
          select case (t)
          case (V_INT);  call eval_node_set_op1_int  (en, sgn_i)
          case (V_REAL); call eval_node_set_op1_real (en, sgn_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("sqrt")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, sqrt_r)
          case (V_CMPLX); call eval_node_set_op1_cmplx (en, sqrt_c)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("exp")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, exp_r)
          case (V_CMPLX); call eval_node_set_op1_cmplx (en, exp_c)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("log")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, log_r)
          case (V_CMPLX); call eval_node_set_op1_cmplx (en, log_c)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("log10")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, log10_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("sin")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, sin_r)
          case (V_CMPLX); call eval_node_set_op1_cmplx (en, sin_c)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("cos")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, cos_r)
          case (V_CMPLX); call eval_node_set_op1_cmplx (en, cos_c)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("tan")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, tan_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("asin")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, asin_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("acos")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, acos_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("atan")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, atan_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("sinh")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, sinh_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("cosh")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, cosh_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("tanh")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, tanh_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("asinh")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, asinh_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("acosh")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, acosh_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case ("atanh")
          select case (t)
          case (V_REAL); call eval_node_set_op1_real (en, atanh_r)
          case default;  call eval_type_error (pn, char (key), t)
          end select
       case default
          call parse_node_mismatch ("function name", pn_fname)
       end select
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done function"
    end if
  end subroutine eval_node_compile_unary_function

  recursive subroutine eval_node_compile_binary_function (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_fname, pn_arg, pn_arg1, pn_arg2
    type(eval_node_t), pointer :: en1, en2
    type(string_t) :: key
    integer :: t1, t2
    if (debug_active (D_MODEL_F)) then
       print *, "read binary function";  call parse_node_write (pn)
    end if
    pn_fname => parse_node_get_sub_ptr (pn)
    pn_arg => parse_node_get_next_ptr (pn_fname, tag="function_arg2")
    pn_arg1 => parse_node_get_sub_ptr (pn_arg, tag="expr")
    pn_arg2 => parse_node_get_next_ptr (pn_arg1, tag="expr")
    call eval_node_compile_expr (en1, pn_arg1, var_list)
    call eval_node_compile_expr (en2, pn_arg2, var_list)
    t1 = en1%result_type
    t2 = en2%result_type
    allocate (en)
    key = parse_node_get_key (pn_fname)
    if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
       select case (char (key))
       case ("max")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_init_int  (en, max_ii (en1, en2))
             case (V_REAL); call eval_node_init_real (en, max_ir (en1, en2))
             case default;  call eval_type_error (pn, char (key), t2)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_init_real (en, max_ri (en1, en2))
             case (V_REAL); call eval_node_init_real (en, max_rr (en1, en2))
             case default;  call eval_type_error (pn, char (key), t2)
             end select
           case default;  call eval_type_error (pn, char (key), t1)
         end select
       case ("min")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_init_int  (en, min_ii (en1, en2))
             case (V_REAL); call eval_node_init_real (en, min_ir (en1, en2))
             case default;  call eval_type_error (pn, char (key), t2)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_init_real (en, min_ri (en1, en2))
             case (V_REAL); call eval_node_init_real (en, min_rr (en1, en2))
             case default;  call eval_type_error (pn, char (key), t2)
             end select
           case default;  call eval_type_error (pn, char (key), t1)
         end select
       case ("mod")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_init_int  (en, mod_ii (en1, en2))
             case (V_REAL); call eval_node_init_real (en, mod_ir (en1, en2))
             case default;  call eval_type_error (pn, char (key), t2)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_init_real (en, mod_ri (en1, en2))
             case (V_REAL); call eval_node_init_real (en, mod_rr (en1, en2))
             case default;  call eval_type_error (pn, char (key), t2)
             end select
           case default;  call eval_type_error (pn, char (key), t1)
          end select
       case ("modulo")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_init_int  (en, modulo_ii (en1, en2))
             case (V_REAL); call eval_node_init_real (en, modulo_ir (en1, en2))
             case default;  call eval_type_error (pn, char (key), t2)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_init_real (en, modulo_ri (en1, en2))
             case (V_REAL); call eval_node_init_real (en, modulo_rr (en1, en2))
             case default;  call eval_type_error (pn, char (key), t2)
             end select
           case default;  call eval_type_error (pn, char (key), t2)
         end select
       case default
          call parse_node_mismatch ("function name", pn_fname)
       end select
       call eval_node_final_rec (en1)
       deallocate (en1)
    else
       call eval_node_init_branch (en, key, t1, en1, en2)
       select case (char (key))
       case ("max")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_int  (en, max_ii)
             case (V_REAL); call eval_node_set_op2_real (en, max_ir)
             case default;  call eval_type_error (pn, char (key), t2)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_real (en, max_ri)
             case (V_REAL); call eval_node_set_op2_real (en, max_rr)
             case default;  call eval_type_error (pn, char (key), t2)
             end select
           case default;  call eval_type_error (pn, char (key), t2)
         end select
       case ("min")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_int  (en, min_ii)
             case (V_REAL); call eval_node_set_op2_real (en, min_ir)
             case default;  call eval_type_error (pn, char (key), t2)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_real (en, min_ri)
             case (V_REAL); call eval_node_set_op2_real (en, min_rr)
             case default;  call eval_type_error (pn, char (key), t2)
             end select
           case default;  call eval_type_error (pn, char (key), t2)
         end select
       case ("mod")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_int  (en, mod_ii)
             case (V_REAL); call eval_node_set_op2_real (en, mod_ir)
             case default;  call eval_type_error (pn, char (key), t2)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_real (en, mod_ri)
             case (V_REAL); call eval_node_set_op2_real (en, mod_rr)
             case default;  call eval_type_error (pn, char (key), t2)
             end select
           case default;  call eval_type_error (pn, char (key), t2)
        end select
       case ("modulo")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_int  (en, modulo_ii)
             case (V_REAL); call eval_node_set_op2_real (en, modulo_ir)
             case default;  call eval_type_error (pn, char (key), t2)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_real (en, modulo_ri)
             case (V_REAL); call eval_node_set_op2_real (en, modulo_rr)
             case default;  call eval_type_error (pn, char (key), t2)
             end select
           case default;  call eval_type_error (pn, char (key), t2)
         end select
       case default
          call parse_node_mismatch ("function name", pn_fname)
       end select
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done function"
    end if
  end subroutine eval_node_compile_binary_function

  recursive subroutine eval_node_compile_block_expr &
       (en, pn, var_list, result_type)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    integer, intent(in), optional :: result_type
    type(parse_node_t), pointer :: pn_var_spec, pn_var_subspec
    type(parse_node_t), pointer :: pn_var_type, pn_var_name, pn_var_expr
    type(parse_node_t), pointer :: pn_expr
    type(string_t) :: var_name
    type(eval_node_t), pointer :: en1, en2
    integer :: var_type
    logical :: new
    if (debug_active (D_MODEL_F)) then
       print *, "read block expr";  call parse_node_write (pn)
    end if
    new = .false.
    pn_var_spec => parse_node_get_sub_ptr (pn, 2)
    select case (char (parse_node_get_rule_key (pn_var_spec)))
    case ("var_num");      var_type = V_NONE
       pn_var_name => parse_node_get_sub_ptr (pn_var_spec)
    case ("var_int");      var_type = V_INT
       new = .true.
       pn_var_name => parse_node_get_sub_ptr (pn_var_spec, 2)
    case ("var_real");     var_type = V_REAL
       new = .true.
       pn_var_name => parse_node_get_sub_ptr (pn_var_spec, 2)
    case ("var_cmplx");     var_type = V_CMPLX
       new = .true.
       pn_var_name => parse_node_get_sub_ptr (pn_var_spec, 2)
    case ("var_logical_new");  var_type = V_LOG
       new = .true.
       pn_var_subspec => parse_node_get_sub_ptr (pn_var_spec, 2)
       pn_var_name => parse_node_get_sub_ptr (pn_var_subspec, 2)
    case ("var_logical_spec");  var_type = V_LOG
       pn_var_name => parse_node_get_sub_ptr (pn_var_spec, 2)
    case ("var_plist_new");    var_type = V_SEV
       new = .true.
       pn_var_subspec => parse_node_get_sub_ptr (pn_var_spec, 2)
       pn_var_name => parse_node_get_sub_ptr (pn_var_subspec, 2)
    case ("var_plist_spec");    var_type = V_SEV
       new = .true.
       pn_var_name => parse_node_get_sub_ptr (pn_var_spec, 2)
    case ("var_alias");    var_type = V_PDG
       new = .true.
       pn_var_name => parse_node_get_sub_ptr (pn_var_spec, 2)
    case ("var_string_new");   var_type = V_STR
       new = .true.
       pn_var_subspec => parse_node_get_sub_ptr (pn_var_spec, 2)
       pn_var_name => parse_node_get_sub_ptr (pn_var_subspec, 2)
    case ("var_string_spec");   var_type = V_STR
       pn_var_name => parse_node_get_sub_ptr (pn_var_spec, 2)
    case default
       call parse_node_mismatch &
            ("logical|int|real|plist|alias", pn_var_type)
    end select
    pn_var_expr => parse_node_get_next_ptr (pn_var_name, 2)
    pn_expr => parse_node_get_next_ptr (pn_var_spec, 2)
    var_name = parse_node_get_string (pn_var_name)
    select case (var_type)
    case (V_LOG);  var_name = "?" // var_name
    case (V_SEV);  var_name = "@" // var_name
    case (V_STR);  var_name = "$" // var_name    ! $ sign
    end select
    call var_list%check_user_var (var_name, var_type, new)
    call eval_node_compile_genexpr (en1, pn_var_expr, var_list, var_type)
    call insert_conversion_node (en1, var_type)
    allocate (en)
    call eval_node_init_block (en, var_name, var_type, en1, var_list)
    call eval_node_compile_genexpr (en2, pn_expr, en%var_list, result_type)
    call eval_node_set_expr (en, en2)
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done block expr"
    end if
  end subroutine eval_node_compile_block_expr

  subroutine insert_conversion_node (en, result_type)
    type(eval_node_t), pointer :: en
    integer, intent(in) :: result_type
    type(eval_node_t), pointer :: en_conv
    select case (en%result_type)
    case (V_INT)
       select case (result_type)
       case (V_REAL)
          allocate (en_conv)
          call eval_node_init_branch (en_conv, var_str ("real"), V_REAL, en)
          call eval_node_set_op1_real (en_conv, real_i)
          en => en_conv
       case (V_CMPLX)
          allocate (en_conv)
          call eval_node_init_branch (en_conv, var_str ("complex"), V_CMPLX, en)
          call eval_node_set_op1_cmplx (en_conv, cmplx_i)
          en => en_conv
       end select
    case (V_REAL)
       select case (result_type)
       case (V_INT)
          allocate (en_conv)
          call eval_node_init_branch (en_conv, var_str ("int"), V_INT, en)
          call eval_node_set_op1_int (en_conv, int_r)
          en => en_conv
       case (V_CMPLX)
          allocate (en_conv)
          call eval_node_init_branch (en_conv, var_str ("complex"), V_CMPLX, en)
          call eval_node_set_op1_cmplx (en_conv, cmplx_r)
          en => en_conv
       end select
    case (V_CMPLX)
       select case (result_type)
       case (V_INT)
          allocate (en_conv)
          call eval_node_init_branch (en_conv, var_str ("int"), V_INT, en)
          call eval_node_set_op1_int (en_conv, int_c)
          en => en_conv
       case (V_REAL)
          allocate (en_conv)
          call eval_node_init_branch (en_conv, var_str ("real"), V_REAL, en)
          call eval_node_set_op1_real (en_conv, real_c)
          en => en_conv
       end select
     case default
     end select
  end subroutine insert_conversion_node

  recursive subroutine eval_node_compile_conditional &
       (en, pn, var_list, result_type)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    integer, intent(in), optional :: result_type
    type(parse_node_t), pointer :: pn_condition, pn_expr
    type(parse_node_t), pointer :: pn_maybe_elsif, pn_elsif_branch
    type(parse_node_t), pointer :: pn_maybe_else, pn_else_branch, pn_else_expr
    type(eval_node_t), pointer :: en0, en1, en2
    integer :: restype
    if (debug_active (D_MODEL_F)) then
       print *, "read conditional";  call parse_node_write (pn)
    end if
    pn_condition => parse_node_get_sub_ptr (pn, 2, tag="lexpr")
    pn_expr => parse_node_get_next_ptr (pn_condition, 2)
    call eval_node_compile_lexpr (en0, pn_condition, var_list)
    call eval_node_compile_genexpr (en1, pn_expr, var_list, result_type)
    if (present (result_type)) then
       restype = major_result_type (result_type, en1%result_type)
    else
       restype = en1%result_type
    end if
    pn_maybe_elsif => parse_node_get_next_ptr (pn_expr)
    select case (char (parse_node_get_rule_key (pn_maybe_elsif)))
    case ("maybe_elsif_expr", &
          "maybe_elsif_lexpr", &
          "maybe_elsif_pexpr", &
          "maybe_elsif_cexpr", &
          "maybe_elsif_sexpr")
       pn_elsif_branch => parse_node_get_sub_ptr (pn_maybe_elsif)
       pn_maybe_else => parse_node_get_next_ptr (pn_maybe_elsif)
       select case (char (parse_node_get_rule_key (pn_maybe_else)))
       case ("maybe_else_expr", &
          "maybe_else_lexpr", &
          "maybe_else_pexpr", &
          "maybe_else_cexpr", &
          "maybe_else_sexpr")
          pn_else_branch => parse_node_get_sub_ptr (pn_maybe_else)
          pn_else_expr => parse_node_get_sub_ptr (pn_else_branch, 2)
       case default
          pn_else_expr => null ()
       end select
       call eval_node_compile_elsif &
            (en2, pn_elsif_branch, pn_else_expr, var_list, restype)
    case ("maybe_else_expr", &
          "maybe_else_lexpr", &
          "maybe_else_pexpr", &
          "maybe_else_cexpr", &
          "maybe_else_sexpr")
       pn_maybe_else => pn_maybe_elsif
       pn_maybe_elsif => null ()
       pn_else_branch => parse_node_get_sub_ptr (pn_maybe_else)
       pn_else_expr => parse_node_get_sub_ptr (pn_else_branch, 2)
       call eval_node_compile_genexpr &
            (en2, pn_else_expr, var_list, restype)
    case ("endif")
       call eval_node_compile_default_else (en2, restype)
    case default
       call msg_bug ("Broken conditional: unexpected " &
            // char (parse_node_get_rule_key (pn_maybe_elsif)))
    end select
    call eval_node_create_conditional (en, en0, en1, en2, restype)
    call conditional_insert_conversion_nodes (en, restype)
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done conditional"
    end if
  end subroutine eval_node_compile_conditional

  recursive subroutine eval_node_compile_elsif &
       (en, pn, pn_else_expr, var_list, result_type)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in), target :: pn
    type(parse_node_t), pointer :: pn_else_expr
    type(var_list_t), intent(in), target :: var_list
    integer, intent(inout) :: result_type
    type(parse_node_t), pointer :: pn_next, pn_condition, pn_expr
    type(eval_node_t), pointer :: en0, en1, en2
    pn_condition => parse_node_get_sub_ptr (pn, 2, tag="lexpr")
    pn_expr => parse_node_get_next_ptr (pn_condition, 2)
    call eval_node_compile_lexpr (en0, pn_condition, var_list)
    call eval_node_compile_genexpr (en1, pn_expr, var_list, result_type)
    result_type = major_result_type (result_type, en1%result_type)
    pn_next => parse_node_get_next_ptr (pn)
    if (associated (pn_next)) then
       call eval_node_compile_elsif &
            (en2, pn_next, pn_else_expr, var_list, result_type)
       result_type = major_result_type (result_type, en2%result_type)
    else if (associated (pn_else_expr)) then
       call eval_node_compile_genexpr &
            (en2, pn_else_expr, var_list, result_type)
       result_type = major_result_type (result_type, en2%result_type)
    else
       call eval_node_compile_default_else (en2, result_type)
    end if
    call eval_node_create_conditional (en, en0, en1, en2, result_type)
  end subroutine eval_node_compile_elsif

  subroutine eval_node_compile_default_else (en, result_type)
    type(eval_node_t), pointer :: en
    integer, intent(in) :: result_type
    type(subevt_t) :: pval_empty
    type(pdg_array_t) :: aval_undefined
    allocate (en)
    select case (result_type)
    case (V_LOG);  call eval_node_init_log (en, .false.)
    case (V_INT);  call eval_node_init_int (en, 0)
    case (V_REAL);  call eval_node_init_real (en, 0._default)
    case (V_CMPLX)
         call eval_node_init_cmplx (en, (0._default, 0._default))
    case (V_SEV)
       call subevt_init (pval_empty)
       call eval_node_init_subevt (en, pval_empty)
    case (V_PDG)
       call eval_node_init_pdg_array  (en, aval_undefined)
    case (V_STR)
       call eval_node_init_string (en, var_str (""))
    case default
       call msg_bug ("Undefined type for 'else' branch in conditional")
    end select
  end subroutine eval_node_compile_default_else

  subroutine eval_node_create_conditional (en, en0, en1, en2, result_type)
    type(eval_node_t), pointer :: en, en0, en1, en2
    integer, intent(in) :: result_type
    if (en0%type == EN_CONSTANT) then
       if (en0%lval) then
          en => en1
          call eval_node_final_rec (en2)
          deallocate (en2)
       else
          en => en2
          call eval_node_final_rec (en1)
          deallocate (en1)
       end if
    else
       allocate (en)
       call eval_node_init_conditional (en, result_type, en0, en1, en2)
    end if
  end subroutine eval_node_create_conditional

  function major_result_type (t1, t2) result (t)
    integer :: t
    integer, intent(in) :: t1, t2
    select case (t1)
    case (V_INT)
       select case (t2)
       case (V_INT, V_REAL, V_CMPLX)
          t = t2
       case default
          call type_mismatch ()
       end select
    case (V_REAL)
       select case (t2)
       case (V_INT)
          t = t1
       case (V_REAL, V_CMPLX)
          t = t2
       case default
          call type_mismatch ()
       end select
    case (V_CMPLX)
       select case (t2)
       case (V_INT, V_REAL, V_CMPLX)
          t = t1
       case default
          call type_mismatch ()
       end select
    case default
       if (t1 == t2) then
          t = t1
       else
          call type_mismatch ()
       end if
    end select
  contains
    subroutine type_mismatch ()
      call msg_bug ("Type mismatch in branches of a conditional expression")
    end subroutine type_mismatch
  end function major_result_type

  recursive subroutine conditional_insert_conversion_nodes (en, result_type)
    type(eval_node_t), intent(inout), target :: en
    integer, intent(in) :: result_type
    select case (result_type)
    case (V_INT, V_REAL, V_CMPLX)
       call insert_conversion_node (en%arg1, result_type)
       if (en%arg2%type == EN_CONDITIONAL) then
          call conditional_insert_conversion_nodes (en%arg2, result_type)
       else
          call insert_conversion_node (en%arg2, result_type)
       end if
    end select
  end subroutine conditional_insert_conversion_nodes

  recursive subroutine eval_node_compile_lexpr (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_term, pn_sequel, pn_arg
    type(eval_node_t), pointer :: en1, en2
    if (debug_active (D_MODEL_F)) then
       print *, "read lexpr";  call parse_node_write (pn)
    end if
    pn_term => parse_node_get_sub_ptr (pn, tag="lsinglet")
    call eval_node_compile_lsinglet (en, pn_term, var_list)
    pn_sequel => parse_node_get_next_ptr (pn_term, tag="lsequel")
    do while (associated (pn_sequel))
       pn_arg => parse_node_get_sub_ptr (pn_sequel, 2, tag="lsinglet")
       en1 => en
       call eval_node_compile_lsinglet (en2, pn_arg, var_list)
       allocate (en)
       if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
          call eval_node_init_log (en, ignore_first_ll (en1, en2))
          call eval_node_final_rec (en1)
          call eval_node_final_rec (en2)
          deallocate (en1, en2)
       else
          call eval_node_init_branch &
               (en, var_str ("lsequel"), V_LOG, en1, en2)
          call eval_node_set_op2_log (en, ignore_first_ll)
       end if
       pn_sequel => parse_node_get_next_ptr (pn_sequel)
    end do
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done lexpr"
    end if
  end subroutine eval_node_compile_lexpr

  recursive subroutine eval_node_compile_lsinglet (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_term, pn_alternative, pn_arg
    type(eval_node_t), pointer :: en1, en2
    if (debug_active (D_MODEL_F)) then
       print *, "read lsinglet";  call parse_node_write (pn)
    end if
    pn_term => parse_node_get_sub_ptr (pn, tag="lterm")
    call eval_node_compile_lterm (en, pn_term, var_list)
    pn_alternative => parse_node_get_next_ptr (pn_term, tag="alternative")
    do while (associated (pn_alternative))
       pn_arg => parse_node_get_sub_ptr (pn_alternative, 2, tag="lterm")
       en1 => en
       call eval_node_compile_lterm (en2, pn_arg, var_list)
       allocate (en)
       if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
          call eval_node_init_log (en, or_ll (en1, en2))
          call eval_node_final_rec (en1)
          call eval_node_final_rec (en2)
          deallocate (en1, en2)
       else
          call eval_node_init_branch &
               (en, var_str ("alternative"), V_LOG, en1, en2)
          call eval_node_set_op2_log (en, or_ll)
       end if
       pn_alternative => parse_node_get_next_ptr (pn_alternative)
    end do
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done lsinglet"
    end if
  end subroutine eval_node_compile_lsinglet

  recursive subroutine eval_node_compile_lterm (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_term, pn_coincidence, pn_arg
    type(eval_node_t), pointer :: en1, en2
    if (debug_active (D_MODEL_F)) then
       print *, "read lterm";  call parse_node_write (pn)
    end if
    pn_term => parse_node_get_sub_ptr (pn)
    call eval_node_compile_lvalue (en, pn_term, var_list)
    pn_coincidence => parse_node_get_next_ptr (pn_term, tag="coincidence")
    do while (associated (pn_coincidence))
       pn_arg => parse_node_get_sub_ptr (pn_coincidence, 2)
       en1 => en
       call eval_node_compile_lvalue (en2, pn_arg, var_list)
       allocate (en)
       if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
          call eval_node_init_log (en, and_ll (en1, en2))
          call eval_node_final_rec (en1)
          call eval_node_final_rec (en2)
          deallocate (en1, en2)
       else
          call eval_node_init_branch &
               (en, var_str ("coincidence"), V_LOG, en1, en2)
          call eval_node_set_op2_log (en, and_ll)
       end if
       pn_coincidence => parse_node_get_next_ptr (pn_coincidence)
    end do
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done lterm"
    end if
  end subroutine eval_node_compile_lterm

  recursive subroutine eval_node_compile_lvalue (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    if (debug_active (D_MODEL_F)) then
       print *, "read lvalue";  call parse_node_write (pn)
    end if
    select case (char (parse_node_get_rule_key (pn)))
    case ("true")
       allocate (en)
       call eval_node_init_log (en, .true.)
    case ("false")
       allocate (en)
       call eval_node_init_log (en, .false.)
    case ("negation")
       call eval_node_compile_negation (en, pn, var_list)
    case ("lvariable")
       call eval_node_compile_variable (en, pn, var_list, V_LOG)
    case ("lexpr")
       call eval_node_compile_lexpr (en, pn, var_list)
    case ("block_lexpr")
       call eval_node_compile_block_expr (en, pn, var_list, V_LOG)
    case ("conditional_lexpr")
       call eval_node_compile_conditional (en, pn, var_list, V_LOG)
    case ("compared_expr")
       call eval_node_compile_compared_expr (en, pn, var_list, V_REAL)
    case ("compared_sexpr")
       call eval_node_compile_compared_expr (en, pn, var_list, V_STR)
    case ("all_fun", "any_fun", "no_fun", "photon_isolation_fun")
       call eval_node_compile_log_function (en, pn, var_list)
    case ("record_cmd")
       call eval_node_compile_record_cmd (en, pn, var_list)
    case default
       call parse_node_mismatch &
            ("true|false|negation|lvariable|" // &
             "lexpr|block_lexpr|conditional_lexpr|" // &
             "compared_expr|compared_sexpr|logical_pexpr", pn)
    end select
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done lvalue"
    end if
  end subroutine eval_node_compile_lvalue

  recursive subroutine eval_node_compile_negation (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_arg
    type(eval_node_t), pointer :: en1
    if (debug_active (D_MODEL_F)) then
       print *, "read negation";  call parse_node_write (pn)
    end if
    pn_arg => parse_node_get_sub_ptr (pn, 2)
    call eval_node_compile_lvalue (en1, pn_arg, var_list)
    allocate (en)
    if (en1%type == EN_CONSTANT) then
       call eval_node_init_log (en, not_l (en1))
       call eval_node_final_rec (en1)
       deallocate (en1)
    else
       call eval_node_init_branch (en, var_str ("not"), V_LOG, en1)
       call eval_node_set_op1_log (en, not_l)
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done negation"
    end if
  end subroutine eval_node_compile_negation

  recursive subroutine eval_node_compile_compared_expr (en, pn, var_list, type)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    integer, intent(in) :: type
    type(parse_node_t), pointer :: pn_comparison, pn_expr1
    type(eval_node_t), pointer :: en0, en1, en2
    if (debug_active (D_MODEL_F)) then
       print *, "read comparison";  call parse_node_write (pn)
    end if
    select case (type)
    case (V_INT, V_REAL)
       pn_expr1 => parse_node_get_sub_ptr (pn, tag="expr")
       call eval_node_compile_expr (en1, pn_expr1, var_list)
       pn_comparison => parse_node_get_next_ptr (pn_expr1, tag="comparison")
    case (V_STR)
       pn_expr1 => parse_node_get_sub_ptr (pn, tag="sexpr")
       call eval_node_compile_sexpr (en1, pn_expr1, var_list)
       pn_comparison => parse_node_get_next_ptr (pn_expr1, tag="str_comparison")
    end select
    call eval_node_compile_comparison &
         (en, en1, en2, pn_comparison, var_list, type)
    pn_comparison => parse_node_get_next_ptr (pn_comparison)
    SCAN_FURTHER: do while (associated (pn_comparison))
       if (en%type == EN_CONSTANT) then
          if (en%lval) then
             en1 => en2
             call eval_node_final_rec (en);  deallocate (en)
             call eval_node_compile_comparison &
                  (en, en1, en2, pn_comparison, var_list, type)
          else
             exit SCAN_FURTHER
          end if
       else
          allocate (en1)
          if (en2%type == EN_CONSTANT) then
             select case (en2%result_type)
             case (V_INT);  call eval_node_init_int    (en1, en2%ival)
             case (V_REAL); call eval_node_init_real   (en1, en2%rval)
             case (V_STR);  call eval_node_init_string (en1, en2%sval)
             end select
          else
             select case (en2%result_type)
             case (V_INT);  call eval_node_init_int_ptr &
                  (en1, var_str ("(previous)"), en2%ival, en2%value_is_known)
             case (V_REAL); call eval_node_init_real_ptr &
                  (en1, var_str ("(previous)"), en2%rval, en2%value_is_known)
             case (V_STR);  call eval_node_init_string_ptr &
                  (en1, var_str ("(previous)"), en2%sval, en2%value_is_known)
             end select
          end if
          en0 => en
          call eval_node_compile_comparison &
               (en, en1, en2, pn_comparison, var_list, type)
          if (en%type == EN_CONSTANT) then
             if (en%lval) then
                call eval_node_final_rec (en);  deallocate (en)
                en => en0
             else
                call eval_node_final_rec (en0);  deallocate (en0)
                exit SCAN_FURTHER
             end if
          else
             en1 => en
             allocate (en)
             call eval_node_init_branch (en, var_str ("and"), V_LOG, en0, en1)
             call eval_node_set_op2_log (en, and_ll)
          end if
       end if
       pn_comparison => parse_node_get_next_ptr (pn_comparison)
    end do SCAN_FURTHER
    if (en%type == EN_CONSTANT .and. associated (en2)) then
       call eval_node_final_rec (en2);  deallocate (en2)
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done compared_expr"
    end if
  end subroutine eval_node_compile_compared_expr

  recursive subroutine eval_node_compile_comparison &
       (en, en1, en2, pn, var_list, type)
    type(eval_node_t), pointer :: en, en1, en2
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    integer, intent(in) :: type
    type(parse_node_t), pointer :: pn_op, pn_arg
    type(string_t) :: key
    integer :: t1, t2
    real(default), pointer :: tolerance_ptr
    pn_op => parse_node_get_sub_ptr (pn)
    key = parse_node_get_key (pn_op)
    select case (type)
    case (V_INT, V_REAL)
       pn_arg => parse_node_get_next_ptr (pn_op, tag="expr")
       call eval_node_compile_expr (en2, pn_arg, var_list)
    case (V_STR)
       pn_arg => parse_node_get_next_ptr (pn_op, tag="sexpr")
       call eval_node_compile_sexpr (en2, pn_arg, var_list)
    end select
    t1 = en1%result_type
    t2 = en2%result_type
    allocate (en)
    if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
       call var_list%get_rptr (var_str ("tolerance"), tolerance_ptr)
       en1%tolerance => tolerance_ptr
       select case (char (key))
       case ("<")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_lt_ii (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_ll_ir (en1, en2))
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_ll_ri (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_ll_rr (en1, en2))
             end select
          end select
       case (">")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_gt_ii (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_gg_ir (en1, en2))
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_gg_ri (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_gg_rr (en1, en2))
             end select
          end select
       case ("<=")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_le_ii (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_ls_ir (en1, en2))
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_ls_ri (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_ls_rr (en1, en2))
             end select
          end select
       case (">=")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_ge_ii (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_gs_ir (en1, en2))
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_gs_ri (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_gs_rr (en1, en2))
             end select
          end select
       case ("==")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_eq_ii (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_se_ir (en1, en2))
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_se_ri (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_se_rr (en1, en2))
             end select
          case (V_STR)
             select case (t2)
             case (V_STR);  call eval_node_init_log (en, comp_eq_ss (en1, en2))
             end select
          end select
       case ("<>")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_ne_ii (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_ns_ir (en1, en2))
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_init_log (en, comp_ns_ri (en1, en2))
             case (V_REAL); call eval_node_init_log (en, comp_ns_rr (en1, en2))
             end select
          case (V_STR)
             select case (t2)
             case (V_STR);  call eval_node_init_log (en, comp_ne_ss (en1, en2))
             end select
          end select
       end select
       call eval_node_final_rec (en1)
       deallocate (en1)
    else
       call eval_node_init_branch (en, key, V_LOG, en1, en2)
       select case (char (key))
       case ("<")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_lt_ii)
             case (V_REAL); call eval_node_set_op2_log (en, comp_ll_ir)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_ll_ri)
             case (V_REAL); call eval_node_set_op2_log (en, comp_ll_rr)
             end select
          end select
       case (">")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_gt_ii)
             case (V_REAL); call eval_node_set_op2_log (en, comp_gg_ir)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_gg_ri)
             case (V_REAL); call eval_node_set_op2_log (en, comp_gg_rr)
             end select
          end select
       case ("<=")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_le_ii)
             case (V_REAL); call eval_node_set_op2_log (en, comp_ls_ir)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_ls_ri)
             case (V_REAL); call eval_node_set_op2_log (en, comp_ls_rr)
             end select
          end select
       case (">=")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_ge_ii)
             case (V_REAL); call eval_node_set_op2_log (en, comp_gs_ir)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_gs_ri)
             case (V_REAL); call eval_node_set_op2_log (en, comp_gs_rr)
             end select
          end select
       case ("==")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_eq_ii)
             case (V_REAL); call eval_node_set_op2_log (en, comp_se_ir)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_se_ri)
             case (V_REAL); call eval_node_set_op2_log (en, comp_se_rr)
             end select
          case (V_STR)
             select case (t2)
             case (V_STR);  call eval_node_set_op2_log (en, comp_eq_ss)
             end select
          end select
       case ("<>")
          select case (t1)
          case (V_INT)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_ne_ii)
             case (V_REAL); call eval_node_set_op2_log (en, comp_ns_ir)
             end select
          case (V_REAL)
             select case (t2)
             case (V_INT);  call eval_node_set_op2_log (en, comp_ns_ri)
             case (V_REAL); call eval_node_set_op2_log (en, comp_ns_rr)
             end select
          case (V_STR)
             select case (t2)
             case (V_STR);  call eval_node_set_op2_log (en, comp_ne_ss)
             end select
          end select
       end select
       call var_list%get_rptr (var_str ("tolerance"), tolerance_ptr)
       en1%tolerance => tolerance_ptr
    end if
  end subroutine eval_node_compile_comparison

  recursive subroutine eval_node_compile_record_cmd (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_key, pn_tag, pn_arg
    type(parse_node_t), pointer :: pn_arg1, pn_arg2, pn_arg3, pn_arg4
    type(eval_node_t), pointer :: en0, en1, en2, en3, en4
    real(default), pointer :: event_weight
    if (debug_active (D_MODEL_F)) then
       print *, "read record_cmd";  call parse_node_write (pn)
    end if
    pn_key => parse_node_get_sub_ptr (pn)
    pn_tag => parse_node_get_next_ptr (pn_key)
    pn_arg => parse_node_get_next_ptr (pn_tag)
    select case (char (parse_node_get_key (pn_key)))
    case ("record")
       call var_list%get_rptr (var_str ("event_weight"), event_weight)
    case ("record_unweighted")
       event_weight => null ()
    case ("record_excess")
       call var_list%get_rptr (var_str ("event_excess"), event_weight)
    end select
    select case (char (parse_node_get_rule_key (pn_tag)))
    case ("analysis_id")
       allocate (en0)
       call eval_node_init_string (en0, parse_node_get_string (pn_tag))
    case default
       call eval_node_compile_sexpr (en0, pn_tag, var_list)
    end select
    allocate (en)
    if (associated (pn_arg)) then
       pn_arg1 => parse_node_get_sub_ptr (pn_arg)
       call eval_node_compile_expr (en1, pn_arg1, var_list)
       if (en1%result_type == V_INT) &
            call insert_conversion_node (en1, V_REAL)
       pn_arg2 => parse_node_get_next_ptr (pn_arg1)
       if (associated (pn_arg2)) then
          call eval_node_compile_expr (en2, pn_arg2, var_list)
          if (en2%result_type == V_INT) &
               call insert_conversion_node (en2, V_REAL)
          pn_arg3 => parse_node_get_next_ptr (pn_arg2)
          if (associated (pn_arg3)) then
             call eval_node_compile_expr (en3, pn_arg3, var_list)
             if (en3%result_type == V_INT) &
                  call insert_conversion_node (en3, V_REAL)
             pn_arg4 => parse_node_get_next_ptr (pn_arg3)
             if (associated (pn_arg4)) then
                call eval_node_compile_expr (en4, pn_arg4, var_list)
                if (en4%result_type == V_INT) &
                     call insert_conversion_node (en4, V_REAL)
                call eval_node_init_record_cmd &
                     (en, event_weight, en0, en1, en2, en3, en4)
             else
                call eval_node_init_record_cmd &
                     (en, event_weight, en0, en1, en2, en3)
             end if
          else
             call eval_node_init_record_cmd (en, event_weight, en0, en1, en2)
          end if
       else
          call eval_node_init_record_cmd (en, event_weight, en0, en1)
       end if
    else
       call eval_node_init_record_cmd (en, event_weight, en0)
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done record_cmd"
    end if
  end subroutine eval_node_compile_record_cmd

  recursive subroutine eval_node_compile_pexpr (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_pterm, pn_concatenation, pn_op, pn_arg
    type(eval_node_t), pointer :: en1, en2
    type(subevt_t) :: subevt
    if (debug_active (D_MODEL_F)) then
       print *, "read pexpr";  call parse_node_write (pn)
    end if
    pn_pterm => parse_node_get_sub_ptr (pn)
    call eval_node_compile_pterm (en, pn_pterm, var_list)
    pn_concatenation => &
         parse_node_get_next_ptr (pn_pterm, tag="pconcatenation")
    do while (associated (pn_concatenation))
       pn_op => parse_node_get_sub_ptr (pn_concatenation)
       pn_arg => parse_node_get_next_ptr (pn_op)
       en1 => en
       call eval_node_compile_pterm (en2, pn_arg, var_list)
       allocate (en)
       if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
          call subevt_join (subevt, en1%pval, en2%pval)
          call eval_node_init_subevt (en, subevt)
          call eval_node_final_rec (en1)
          call eval_node_final_rec (en2)
          deallocate (en1, en2)
       else
          call eval_node_init_branch &
               (en, var_str ("join"), V_SEV, en1, en2)
          call eval_node_set_op2_sev (en, join_pp)
       end if
       pn_concatenation => parse_node_get_next_ptr (pn_concatenation)
    end do
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done pexpr"
    end if
  end subroutine eval_node_compile_pexpr

  recursive subroutine eval_node_compile_pterm (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_pvalue, pn_combination, pn_op, pn_arg
    type(eval_node_t), pointer :: en1, en2
    type(subevt_t) :: subevt
    if (debug_active (D_MODEL_F)) then
       print *, "read pterm";  call parse_node_write (pn)
    end if
    pn_pvalue => parse_node_get_sub_ptr (pn)
    call eval_node_compile_pvalue (en, pn_pvalue, var_list)
    pn_combination => &
         parse_node_get_next_ptr (pn_pvalue, tag="pcombination")
    do while (associated (pn_combination))
       pn_op => parse_node_get_sub_ptr (pn_combination)
       pn_arg => parse_node_get_next_ptr (pn_op)
       en1 => en
       call eval_node_compile_pvalue (en2, pn_arg, var_list)
       allocate (en)
       if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
          call subevt_combine (subevt, en1%pval, en2%pval)
          call eval_node_init_subevt (en, subevt)
          call eval_node_final_rec (en1)
          call eval_node_final_rec (en2)
          deallocate (en1, en2)
       else
          call eval_node_init_branch &
               (en, var_str ("combine"), V_SEV, en1, en2)
          call eval_node_set_op2_sev (en, combine_pp)
       end if
       pn_combination => parse_node_get_next_ptr (pn_combination)
    end do
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done pterm"
    end if
  end subroutine eval_node_compile_pterm

  recursive subroutine eval_node_compile_pvalue (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_prefix_cexpr
    type(eval_node_t), pointer :: en1, en2, en0
    type(string_t) :: key
    type(subevt_t), pointer :: evt_ptr
    logical, pointer :: known
    if (debug_active (D_MODEL_F)) then
       print *, "read pvalue";  call parse_node_write (pn)
    end if
    select case (char (parse_node_get_rule_key (pn)))
    case ("pexpr_src")
       call eval_node_compile_prefix_cexpr (en1, pn, var_list)
       allocate (en2)
       if (var_list%contains (var_str ("@evt"))) then
          call var_list%get_pptr (var_str ("@evt"), evt_ptr, known)
          call eval_node_init_subevt_ptr (en2, var_str ("@evt"), evt_ptr, known)
          allocate (en)
          call eval_node_init_branch &
               (en, var_str ("prt_selection"), V_SEV, en1, en2)
          call eval_node_set_op2_sev (en, select_pdg_ca)
          allocate (en0)
          pn_prefix_cexpr => parse_node_get_sub_ptr (pn)
          key = parse_node_get_rule_key (pn_prefix_cexpr)
          select case (char (key))
          case ("beam_prt")
             call eval_node_init_int (en0, PRT_BEAM)
             en%arg0 => en0
          case ("incoming_prt")
             call eval_node_init_int (en0, PRT_INCOMING)
             en%arg0 => en0
          case ("outgoing_prt")
             call eval_node_init_int (en0, PRT_OUTGOING)
             en%arg0 => en0
          case ("unspecified_prt")
             call eval_node_init_int (en0, PRT_OUTGOING)
             en%arg0 => en0
          end select
       else
          call parse_node_write (pn)
          call msg_bug (" Missing event data while compiling pvalue")
       end if
    case ("pvariable")
       call eval_node_compile_variable (en, pn, var_list, V_SEV)
    case ("pexpr")
       call eval_node_compile_pexpr (en, pn, var_list)
    case ("block_pexpr")
       call eval_node_compile_block_expr (en, pn, var_list, V_SEV)
    case ("conditional_pexpr")
       call eval_node_compile_conditional (en, pn, var_list, V_SEV)
    case ("join_fun", "combine_fun", "collect_fun", "cluster_fun", &
          "select_fun", "extract_fun", "sort_fun", "select_b_jet_fun", &
          "select_non_bjet_fun", "select_c_jet_fun", &
          "select_light_jet_fun", "photon_reco_fun")
       call eval_node_compile_prt_function (en, pn, var_list)
    case default
       call parse_node_mismatch &
            ("prefix_cexpr|pvariable|" // &
             "grouped_pexpr|block_pexpr|conditional_pexpr|" // &
             "prt_function", pn)
    end select
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done pvalue"
    end if
  end subroutine eval_node_compile_pvalue

  recursive subroutine eval_node_compile_prt_function (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_clause, pn_key, pn_cond, pn_args
    type(parse_node_t), pointer :: pn_arg0, pn_arg1, pn_arg2
    type(eval_node_t), pointer :: en0, en1, en2
    type(string_t) :: key
    if (debug_active (D_MODEL_F)) then
       print *, "read prt_function";  call parse_node_write (pn)
    end if
    pn_clause => parse_node_get_sub_ptr (pn)
    pn_key  => parse_node_get_sub_ptr (pn_clause)
    pn_cond => parse_node_get_next_ptr (pn_key)
    if (associated (pn_cond)) &
         pn_arg0 => parse_node_get_sub_ptr (pn_cond, 2)
    pn_args => parse_node_get_next_ptr (pn_clause)
    pn_arg1 => parse_node_get_sub_ptr (pn_args)
    pn_arg2 => parse_node_get_next_ptr (pn_arg1)
    key = parse_node_get_key (pn_key)
    call eval_node_compile_pexpr (en1, pn_arg1, var_list)
    allocate (en)
    if (.not. associated (pn_arg2)) then
       select case (char (key))
       case ("collect")
          call eval_node_init_prt_fun_unary (en, en1, key, collect_p)
       case ("cluster")
          if (fastjet_available ()) then
             call fastjet_init ()
          else
             call msg_fatal &
               ("'cluster' function requires FastJet, which is not enabled")
          end if
          en1%var_list => var_list
          call eval_node_init_prt_fun_unary (en, en1, key, cluster_p)
          call var_list%get_iptr (var_str ("jet_algorithm"), en1%jet_algorithm)
          call var_list%get_rptr (var_str ("jet_r"), en1%jet_r)
          call var_list%get_rptr (var_str ("jet_p"), en1%jet_p)
          call var_list%get_rptr (var_str ("jet_ycut"), en1%jet_ycut)
          call var_list%get_rptr (var_str ("jet_dcut"), en1%jet_dcut)
       case ("photon_recombination")
          en1%var_list => var_list
          call eval_node_init_prt_fun_unary &
               (en, en1, key, photon_recombination_p)
          call var_list%get_rptr (var_str ("photon_rec_r0"), en1%photon_rec_r0)
       case ("select")
          call eval_node_init_prt_fun_unary (en, en1, key, select_p)
       case ("extract")
          call eval_node_init_prt_fun_unary (en, en1, key, extract_p)
       case ("sort")
          call eval_node_init_prt_fun_unary (en, en1, key, sort_p)
       case ("select_b_jet")
          call eval_node_init_prt_fun_unary (en, en1, key, select_b_jet_p)
       case ("select_non_b_jet")
          call eval_node_init_prt_fun_unary (en, en1, key, select_non_b_jet_p)
       case ("select_c_jet")
          call eval_node_init_prt_fun_unary (en, en1, key, select_c_jet_p)
       case ("select_light_jet")
          call eval_node_init_prt_fun_unary (en, en1, key, select_light_jet_p)
       case default
          call msg_bug (" Unary particle function '" // char (key) // &
               "' undefined")
       end select
    else
       call eval_node_compile_pexpr (en2, pn_arg2, var_list)
       select case (char (key))
       case ("join")
          call eval_node_init_prt_fun_binary (en, en1, en2, key, join_pp)
       case ("combine")
          call eval_node_init_prt_fun_binary (en, en1, en2, key, combine_pp)
       case ("collect")
          call eval_node_init_prt_fun_binary (en, en1, en2, key, collect_pp)
       case ("select")
          call eval_node_init_prt_fun_binary (en, en1, en2, key, select_pp)
       case ("sort")
          call eval_node_init_prt_fun_binary (en, en1, en2, key, sort_pp)
       case default
          call msg_bug (" Binary particle function '" // char (key) // &
               "' undefined")
       end select
    end if
    if (associated (pn_cond)) then
       call eval_node_set_observables (en, var_list)
       select case (char (key))
       case ("extract", "sort")
          call eval_node_compile_expr (en0, pn_arg0, en%var_list)
       case default
          call eval_node_compile_lexpr (en0, pn_arg0, en%var_list)
       end select
       en%arg0 => en0
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done prt_function"
    end if
  end subroutine eval_node_compile_prt_function

  recursive subroutine eval_node_compile_eval_function (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_key, pn_arg0, pn_args, pn_arg1, pn_arg2
    type(eval_node_t), pointer :: en0, en1, en2
    type(string_t) :: key
    if (debug_active (D_MODEL_F)) then
       print *, "read eval_function";  call parse_node_write (pn)
    end if
    pn_key => parse_node_get_sub_ptr (pn)
    pn_arg0 => parse_node_get_next_ptr (pn_key)
    pn_args => parse_node_get_next_ptr (pn_arg0)
    pn_arg1 => parse_node_get_sub_ptr (pn_args)
    pn_arg2 => parse_node_get_next_ptr (pn_arg1)
    key = parse_node_get_key (pn_key)
    call eval_node_compile_pexpr (en1, pn_arg1, var_list)
    allocate (en)
    if (.not. associated (pn_arg2)) then
       call eval_node_init_eval_fun_unary (en, en1, key)
    else
       call eval_node_compile_pexpr (en2, pn_arg2, var_list)
       call eval_node_init_eval_fun_binary (en, en1, en2, key)
    end if
    en%pval => en1%pval
    call eval_node_set_observables (en, var_list)
    call eval_node_compile_expr (en0, pn_arg0, en%var_list)
    if (en0%result_type == V_INT)  &
         call insert_conversion_node (en0, V_REAL)
    if (en0%result_type /= V_REAL) &
         call msg_fatal (" 'eval' function does not result in real value")
    call eval_node_set_expr (en, en0)
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done eval_function"
    end if
  end subroutine eval_node_compile_eval_function

  recursive subroutine eval_node_compile_log_function (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_clause, pn_key, pn_str, pn_cond
    type(parse_node_t), pointer :: pn_arg0, pn_args, pn_arg1, pn_arg2
    type(eval_node_t), pointer :: en0, en1, en2
    type(string_t) :: key
    if (debug_active (D_MODEL_F)) then
       print *, "read log_function";  call parse_node_write (pn)
    end if
    select case (char (parse_node_get_rule_key (pn)))
    case ("all_fun", "any_fun", "no_fun")
       pn_key => parse_node_get_sub_ptr (pn)
       pn_arg0 => parse_node_get_next_ptr (pn_key)
       pn_args => parse_node_get_next_ptr (pn_arg0)
    case ("photon_isolation_fun")
       pn_clause => parse_node_get_sub_ptr (pn)
       pn_key => parse_node_get_sub_ptr (pn_clause)
       pn_cond => parse_node_get_next_ptr (pn_key)
       if (associated (pn_cond)) then
          pn_arg0 => parse_node_get_sub_ptr (pn_cond, 2)
       else
          pn_arg0 => null ()
       end if
       pn_args => parse_node_get_next_ptr (pn_clause)
    case default
       call parse_node_mismatch ("all_fun|any_fun|" // &
            "no_fun|photon_isolation_fun", pn)
    end select
    pn_arg1 => parse_node_get_sub_ptr (pn_args)
    pn_arg2 => parse_node_get_next_ptr (pn_arg1)
    key = parse_node_get_key (pn_key)
    call eval_node_compile_pexpr (en1, pn_arg1, var_list)
    allocate (en)
    if (.not. associated (pn_arg2)) then
       select case (char (key))
       case ("all")
          call eval_node_init_log_fun_unary (en, en1, key, all_p)
       case ("any")
          call eval_node_init_log_fun_unary (en, en1, key, any_p)
       case ("no")
          call eval_node_init_log_fun_unary (en, en1, key, no_p)
       case default
          call msg_bug ("Unary logical particle function '" // char (key) // &
               "' undefined")
       end select
    else
       call eval_node_compile_pexpr (en2, pn_arg2, var_list)
       select case (char (key))
       case ("all")
          call eval_node_init_log_fun_binary (en, en1, en2, key, all_pp)
       case ("any")
          call eval_node_init_log_fun_binary (en, en1, en2, key, any_pp)
       case ("no")
          call eval_node_init_log_fun_binary (en, en1, en2, key, no_pp)
       case ("photon_isolation")
          en1%var_list => var_list
          call var_list%get_rptr (var_str ("photon_iso_eps"), en1%photon_iso_eps)
          call var_list%get_rptr (var_str ("photon_iso_n"), en1%photon_iso_n)
          call var_list%get_rptr (var_str ("photon_iso_r0"), en1%photon_iso_r0)
          call eval_node_init_log_fun_binary (en, en1, en2, key, photon_isolation_pp)
       case default
          call msg_bug ("Binary logical particle function '" // char (key) // &
               "' undefined")
       end select
    end if
    if (associated (pn_arg0)) then
       call eval_node_set_observables (en, var_list)
       select case (char (key))
       case ("all", "any", "no", "photon_isolation")
          call eval_node_compile_lexpr (en0, pn_arg0, en%var_list)
       case default
          call msg_bug ("Compiling logical particle function: missing mode")
       end select
       call eval_node_set_expr (en, en0, V_LOG)
    else
       select case (char (key))
       case ("photon_isolation")
          allocate (en0)
          call eval_node_init_log (en0, .true.)
          call eval_node_set_expr (en, en0, V_LOG)
       case default
          call msg_bug ("Only photon isolation can be called unconditionally")
       end select
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done log_function"
    end if
  end subroutine eval_node_compile_log_function

  recursive subroutine eval_node_compile_count_function (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_clause, pn_key, pn_cond, pn_args
    type(parse_node_t), pointer :: pn_arg0, pn_arg1, pn_arg2
    type(eval_node_t), pointer :: en0, en1, en2
    type(string_t) :: key
    if (debug_active (D_MODEL_F)) then
       print *, "read count_function";  call parse_node_write (pn)
    end if
    select case (char (parse_node_get_rule_key (pn)))
    case ("count_fun")
       pn_clause => parse_node_get_sub_ptr (pn)
       pn_key => parse_node_get_sub_ptr (pn_clause)
       pn_cond => parse_node_get_next_ptr (pn_key)
       if (associated (pn_cond)) then
          pn_arg0 => parse_node_get_sub_ptr (pn_cond, 2)
       else
          pn_arg0 => null ()
       end if
       pn_args => parse_node_get_next_ptr (pn_clause)
    end select
    pn_arg1 => parse_node_get_sub_ptr (pn_args)
    pn_arg2 => parse_node_get_next_ptr (pn_arg1)
    key = parse_node_get_key (pn_key)
    call eval_node_compile_pexpr (en1, pn_arg1, var_list)
    allocate (en)
    if (.not. associated (pn_arg2)) then
       select case (char (key))
       case ("count")
          call eval_node_init_int_fun_unary (en, en1, key, count_a)
       case default
          call msg_bug ("Unary subevent function '" // char (key) // &
               "' undefined")
       end select
    else
       call eval_node_compile_pexpr (en2, pn_arg2, var_list)
       select case (char (key))
       case ("count")
          call eval_node_init_int_fun_binary (en, en1, en2, key, count_pp)
       case default
          call msg_bug ("Binary subevent function '" // char (key) // &
               "' undefined")
       end select
    end if
    if (associated (pn_arg0)) then
       call eval_node_set_observables (en, var_list)
       select case (char (key))
       case ("count")
          call eval_node_compile_lexpr (en0, pn_arg0, en%var_list)
          call eval_node_set_expr (en, en0, V_INT)
       end select
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done count_function"
    end if
  end subroutine eval_node_compile_count_function

  recursive subroutine eval_node_compile_numeric_function (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_key, pn_args
    type(parse_node_t), pointer :: pn_arg0, pn_arg1, pn_arg2
    type(eval_node_t), pointer :: en0, en1
    type(string_t) :: key
    if (debug_active (D_MODEL_F)) then
       print *, "read numeric_function";  call parse_node_write (pn)
    end if
    select case (char (parse_node_get_rule_key (pn)))
    case ("sum_fun", "prod_fun")
       if (debug_active (D_MODEL_F)) then
          print *, "read sum_fun"; call parse_node_write (pn)
       end if
       pn_key => parse_node_get_sub_ptr (pn)
       pn_arg0 => parse_node_get_next_ptr (pn_key)
       pn_args => parse_node_get_next_ptr (pn_arg0)
    end select
    pn_arg1 => parse_node_get_sub_ptr (pn_args)
    pn_arg2 => parse_node_get_next_ptr (pn_arg1)
    key = parse_node_get_key (pn_key)
    call eval_node_compile_pexpr (en1, pn_arg1, var_list)
    if (associated (pn_arg2)) then
       call msg_fatal ("The " // char (key) // &
            " function can only be used for unary observables.")
    end if
    allocate (en)
    select case (char (key))
    case ("sum")
       call eval_node_init_real_fun_cum (en, en1, key, sum_a)
    case ("prod")
       call eval_node_init_real_fun_cum (en, en1, key, prod_a)
    case default
       call msg_bug ("Unary subevent function '" // char (key) // &
            "' undefined")
    end select
    call eval_node_set_observables (en, var_list)
    call eval_node_compile_expr (en0, pn_arg0, en%var_list)
    if (en0%result_type == V_INT) &
         call insert_conversion_node (en0, V_REAL)
    call eval_node_set_expr (en, en0, V_REAL)
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done numeric_function"
    end if
  end subroutine eval_node_compile_numeric_function

  recursive subroutine eval_node_compile_prefix_cexpr (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_avalue, pn_prt
    type(string_t) :: key
    if (debug_active (D_MODEL_F)) then
       print *, "read prefix_cexpr";  call parse_node_write (pn)
    end if
    pn_avalue => parse_node_get_sub_ptr (pn)
    key = parse_node_get_rule_key (pn_avalue)
    select case (char (key))
    case ("beam_prt")
       pn_prt => parse_node_get_sub_ptr (pn_avalue, 2)
       call eval_node_compile_cexpr (en, pn_prt, var_list)
    case ("incoming_prt")
       pn_prt => parse_node_get_sub_ptr (pn_avalue, 2)
       call eval_node_compile_cexpr (en, pn_prt, var_list)
    case ("outgoing_prt")
       pn_prt => parse_node_get_sub_ptr (pn_avalue, 2)
       call eval_node_compile_cexpr (en, pn_prt, var_list)
    case ("unspecified_prt")
       pn_prt => parse_node_get_sub_ptr (pn_avalue, 1)
       call eval_node_compile_cexpr (en, pn_prt, var_list)
    case default
       call parse_node_mismatch &
            ("beam_prt|incoming_prt|outgoing_prt|unspecified_prt", &
             pn_avalue)
    end select
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done prefix_cexpr"
    end if
  end subroutine eval_node_compile_prefix_cexpr

   recursive subroutine eval_node_compile_cexpr (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_prt, pn_concatenation
    type(eval_node_t), pointer :: en1, en2
    type(pdg_array_t) :: aval
    if (debug_active (D_MODEL_F)) then
       print *, "read cexpr";  call parse_node_write (pn)
    end if
    pn_prt => parse_node_get_sub_ptr (pn)
    call eval_node_compile_avalue (en, pn_prt, var_list)
    pn_concatenation => parse_node_get_next_ptr (pn_prt)
    do while (associated (pn_concatenation))
       pn_prt => parse_node_get_sub_ptr (pn_concatenation, 2)
       en1 => en
       call eval_node_compile_avalue (en2, pn_prt, var_list)
       allocate (en)
       if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
          call concat_cc (aval, en1, en2)
          call eval_node_init_pdg_array (en, aval)
          call eval_node_final_rec (en1)
          call eval_node_final_rec (en2)
          deallocate (en1, en2)
       else
          call eval_node_init_branch (en, var_str (":"), V_PDG, en1, en2)
          call eval_node_set_op2_pdg (en, concat_cc)
       end if
       pn_concatenation => parse_node_get_next_ptr (pn_concatenation)
    end do
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done cexpr"
    end if
  end subroutine eval_node_compile_cexpr

  recursive subroutine eval_node_compile_avalue (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    if (debug_active (D_MODEL_F)) then
       print *, "read avalue";  call parse_node_write (pn)
    end if
    select case (char (parse_node_get_rule_key (pn)))
    case ("pdg_code")
       call eval_node_compile_pdg_code (en, pn, var_list)
    case ("cvariable", "variable", "prt_name")
       call eval_node_compile_cvariable (en, pn, var_list)
    case ("cexpr")
       call eval_node_compile_cexpr (en, pn, var_list)
    case ("block_cexpr")
       call eval_node_compile_block_expr (en, pn, var_list, V_PDG)
    case ("conditional_cexpr")
       call eval_node_compile_conditional (en, pn, var_list, V_PDG)
    case default
       call parse_node_mismatch &
            ("grouped_cexpr|block_cexpr|conditional_cexpr|" // &
             "pdg_code|cvariable|prt_name", pn)
    end select
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done avalue"
    end if
  end subroutine eval_node_compile_avalue

  subroutine eval_node_compile_pdg_code (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in), target :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_arg
    type(eval_node_t), pointer :: en1
    type(string_t) :: key
    type(pdg_array_t) :: aval
    integer :: t
    if (debug_active (D_MODEL_F)) then
       print *, "read PDG code";  call parse_node_write (pn)
    end if
    pn_arg => parse_node_get_sub_ptr (pn, 2)
    call eval_node_compile_expr &
         (en1, parse_node_get_sub_ptr (pn_arg, tag="expr"), var_list)
    t = en1%result_type
    allocate (en)
    key = "PDG"
    if (en1%type == EN_CONSTANT) then
       select case (t)
       case (V_INT)
          call pdg_i (aval, en1)
          call eval_node_init_pdg_array (en, aval)
       case default;  call eval_type_error (pn, char (key), t)
       end select
       call eval_node_final_rec (en1)
       deallocate (en1)
    else
       select case (t)
       case (V_INT);  call eval_node_set_op1_pdg (en, pdg_i)
       case default;  call eval_type_error (pn, char (key), t)
       end select
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done function"
    end if
  end subroutine eval_node_compile_pdg_code

  subroutine eval_node_compile_cvariable (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in), target :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_name
    type(string_t) :: var_name
    type(pdg_array_t), pointer :: aptr
    type(pdg_array_t), target, save :: no_aval
    logical, pointer :: known
    logical, target, save :: unknown = .false.
    if (debug_active (D_MODEL_F)) then
       print *, "read cvariable";  call parse_node_write (pn)
    end if
    pn_name => pn
    var_name = parse_node_get_string (pn_name)
    allocate (en)
    if (var_list%contains (var_name)) then
       call var_list%get_aptr (var_name, aptr, known)
       call eval_node_init_pdg_array_ptr (en, var_name, aptr, known)
    else
       call parse_node_write (pn)
       call msg_error ("This PDG-array variable is undefined at this point")
       call eval_node_init_pdg_array_ptr (en, var_name, no_aval, unknown)
    end if
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done cvariable"
    end if
  end subroutine eval_node_compile_cvariable

  recursive subroutine eval_node_compile_sexpr (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_svalue, pn_concatenation, pn_op, pn_arg
    type(eval_node_t), pointer :: en1, en2
    type(string_t) :: string
    if (debug_active (D_MODEL_F)) then
       print *, "read sexpr";  call parse_node_write (pn)
    end if
    pn_svalue => parse_node_get_sub_ptr (pn)
    call eval_node_compile_svalue (en, pn_svalue, var_list)
    pn_concatenation => &
         parse_node_get_next_ptr (pn_svalue, tag="str_concatenation")
    do while (associated (pn_concatenation))
       pn_op => parse_node_get_sub_ptr (pn_concatenation)
       pn_arg => parse_node_get_next_ptr (pn_op)
       en1 => en
       call eval_node_compile_svalue (en2, pn_arg, var_list)
       allocate (en)
       if (en1%type == EN_CONSTANT .and. en2%type == EN_CONSTANT) then
          call concat_ss (string, en1, en2)
          call eval_node_init_string (en, string)
          call eval_node_final_rec (en1)
          call eval_node_final_rec (en2)
          deallocate (en1, en2)
       else
          call eval_node_init_branch &
               (en, var_str ("concat"), V_STR, en1, en2)
          call eval_node_set_op2_str (en, concat_ss)
       end if
       pn_concatenation => parse_node_get_next_ptr (pn_concatenation)
    end do
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done sexpr"
    end if
  end subroutine eval_node_compile_sexpr

  recursive subroutine eval_node_compile_svalue (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    if (debug_active (D_MODEL_F)) then
       print *, "read svalue";  call parse_node_write (pn)
    end if
    select case (char (parse_node_get_rule_key (pn)))
    case ("svariable")
       call eval_node_compile_variable (en, pn, var_list, V_STR)
    case ("sexpr")
       call eval_node_compile_sexpr (en, pn, var_list)
    case ("block_sexpr")
       call eval_node_compile_block_expr (en, pn, var_list, V_STR)
    case ("conditional_sexpr")
       call eval_node_compile_conditional (en, pn, var_list, V_STR)
    case ("sprintf_fun")
       call eval_node_compile_sprintf (en, pn, var_list)
    case ("string_literal")
       allocate (en)
       call eval_node_init_string (en, parse_node_get_string (pn))
    case default
       call parse_node_mismatch &
            ("svariable|" // &
             "grouped_sexpr|block_sexpr|conditional_sexpr|" // &
             "string_function|string_literal", pn)
    end select
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done svalue"
    end if
  end subroutine eval_node_compile_svalue

  recursive subroutine eval_node_compile_sprintf (en, pn, var_list)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_clause, pn_key, pn_args
    type(parse_node_t), pointer :: pn_arg0
    type(eval_node_t), pointer :: en0, en1
    integer :: n_args
    type(string_t) :: key
    if (debug_active (D_MODEL_F)) then
       print *, "read sprintf_fun";  call parse_node_write (pn)
    end if
    pn_clause => parse_node_get_sub_ptr (pn)
    pn_key  => parse_node_get_sub_ptr (pn_clause)
    pn_arg0 => parse_node_get_next_ptr (pn_key)
    pn_args => parse_node_get_next_ptr (pn_clause)
    call eval_node_compile_sexpr (en0, pn_arg0, var_list)
    if (associated (pn_args)) then
       call eval_node_compile_sprintf_args (en1, pn_args, var_list, n_args)
    else
       n_args = 0
       en1 => null ()
    end if
    allocate (en)
    key = parse_node_get_key (pn_key)
    call eval_node_init_format_string (en, en0, en1, key, n_args)
    if (debug_active (D_MODEL_F)) then
       call eval_node_write (en)
       print *, "done sprintf_fun"
    end if
  end subroutine eval_node_compile_sprintf

  subroutine eval_node_compile_sprintf_args (en, pn, var_list, n_args)
    type(eval_node_t), pointer :: en
    type(parse_node_t), intent(in) :: pn
    type(var_list_t), intent(in), target :: var_list
    integer, intent(out) :: n_args
    type(parse_node_t), pointer :: pn_arg
    integer :: i
    type(eval_node_t), pointer :: en1, en2
    n_args = parse_node_get_n_sub (pn)
    en => null ()
    do i = n_args, 1, -1
       pn_arg => parse_node_get_sub_ptr (pn, i)
       select case (char (parse_node_get_rule_key (pn_arg)))
       case ("lvariable")
          call eval_node_compile_variable (en1, pn_arg, var_list, V_LOG)
       case ("svariable")
          call eval_node_compile_variable (en1, pn_arg, var_list, V_STR)
       case ("expr")
          call eval_node_compile_expr (en1, pn_arg, var_list)
       case default
          call parse_node_mismatch ("variable|svariable|lvariable|expr", pn_arg)
       end select
       if (associated (en)) then
          en2 => en
          allocate (en)
          call eval_node_init_branch &
               (en, var_str ("sprintf_arg"), V_NONE, en1, en2)
       else
          allocate (en)
          call eval_node_init_branch &
               (en, var_str ("sprintf_arg"), V_NONE, en1)
       end if
    end do
  end subroutine eval_node_compile_sprintf_args

  subroutine evaluate_sprintf (string, n_args, en_fmt, en_arg)
    type(string_t), intent(out) :: string
    integer, intent(in) :: n_args
    type(eval_node_t), pointer :: en_fmt
    type(eval_node_t), intent(in), optional, target :: en_arg
    type(eval_node_t), pointer :: en_branch, en_var
    type(sprintf_arg_t), dimension(:), allocatable :: arg
    type(string_t) :: fmt
    logical :: autoformat
    integer :: i, j, sprintf_argc
    autoformat = .not. associated (en_fmt)
    if (autoformat) fmt = ""
    if (present (en_arg)) then
       sprintf_argc = 0
       en_branch => en_arg
       do i = 1, n_args
          select case (en_branch%arg1%result_type)
             case (V_CMPLX); sprintf_argc = sprintf_argc + 2
             case default  ; sprintf_argc = sprintf_argc + 1
          end select
          en_branch => en_branch%arg2
       end do
       allocate (arg (sprintf_argc))
       j = 1
       en_branch => en_arg
       do i = 1, n_args
          en_var => en_branch%arg1
          select case (en_var%result_type)
          case (V_LOG)
             call sprintf_arg_init (arg(j), en_var%lval)
             if (autoformat) fmt = fmt // "%s "
          case (V_INT);
             call sprintf_arg_init (arg(j), en_var%ival)
             if (autoformat) fmt = fmt // "%i "
          case (V_REAL);
             call sprintf_arg_init (arg(j), en_var%rval)
             if (autoformat) fmt = fmt // "%g "
          case (V_STR)
             call sprintf_arg_init (arg(j), en_var%sval)
             if (autoformat) fmt = fmt // "%s "
          case (V_CMPLX)
             call sprintf_arg_init (arg(j), real (en_var%cval, default))
             j = j + 1
             call sprintf_arg_init (arg(j), aimag (en_var%cval))
             if (autoformat) fmt = fmt // "(%g + %g * I) "
          case default
             call eval_node_write (en_var)
             call msg_error ("sprintf is implemented " &
                  // "for logical, integer, real, and string values only")
          end select
          j = j + 1
          en_branch => en_branch%arg2
       end do
    else
       allocate (arg(0))
    end if
    if (autoformat) then
       string = sprintf (trim (fmt), arg)
    else
       string = sprintf (en_fmt%sval, arg)
    end if
  end subroutine evaluate_sprintf

  subroutine eval_type_error (pn, string, t)
    type(parse_node_t), intent(in) :: pn
    character(*), intent(in) :: string
    integer, intent(in) :: t
    type(string_t) :: type
    select case (t)
    case (V_NONE); type = "(none)"
    case (V_LOG);  type = "'logical'"
    case (V_INT);  type = "'integer'"
    case (V_REAL); type = "'real'"
    case (V_CMPLX); type = "'complex'"
    case default;  type = "(unknown)"
    end select
    call parse_node_write (pn)
    call msg_fatal (" The " // string // &
         " operation is not defined for the given argument type " // &
         char (type))
  end subroutine eval_type_error

  function numeric_result_type (t1, t2) result (t)
    integer, intent(in) :: t1, t2
    integer :: t
    if (t1 == V_INT .and. t2 == V_INT) then
       t = V_INT
    else if (t1 == V_INT .and. t2 == V_REAL) then
       t = V_REAL
    else if (t1 == V_REAL .and. t2 == V_INT) then
       t = V_REAL
    else if (t1 == V_REAL .and. t2 == V_REAL) then
       t = V_REAL
    else
       t = V_CMPLX
    end if
  end function numeric_result_type

  recursive subroutine eval_node_evaluate (en)
    type(eval_node_t), intent(inout) :: en
    logical :: exist
    select case (en%type)
    case (EN_UNARY)
       if (associated (en%arg1)) then
          call eval_node_evaluate (en%arg1)
          en%value_is_known = en%arg1%value_is_known
       else
          en%value_is_known = .false.
       end if
       if (en%value_is_known) then
          select case (en%result_type)
          case (V_LOG);  en%lval = en%op1_log  (en%arg1)
          case (V_INT);  en%ival = en%op1_int  (en%arg1)
          case (V_REAL); en%rval = en%op1_real (en%arg1)
          case (V_CMPLX); en%cval = en%op1_cmplx (en%arg1)
          case (V_PDG);
             call en%op1_pdg  (en%aval, en%arg1)
          case (V_SEV)
             if (associated (en%arg0)) then
                call en%op1_sev (en%pval, en%arg1, en%arg0)
             else
                call en%op1_sev (en%pval, en%arg1)
             end if
          case (V_STR)
             call en%op1_str (en%sval, en%arg1)
          end select
       end if
    case (EN_BINARY)
       if (associated (en%arg1) .and. associated (en%arg2)) then
          call eval_node_evaluate (en%arg1)
          call eval_node_evaluate (en%arg2)
          en%value_is_known = &
               en%arg1%value_is_known .and. en%arg2%value_is_known
       else
          en%value_is_known = .false.
       end if
       if (en%value_is_known) then
          select case (en%result_type)
          case (V_LOG);  en%lval = en%op2_log  (en%arg1, en%arg2)
          case (V_INT);  en%ival = en%op2_int  (en%arg1, en%arg2)
          case (V_REAL); en%rval = en%op2_real (en%arg1, en%arg2)
          case (V_CMPLX); en%cval = en%op2_cmplx (en%arg1, en%arg2)
          case (V_PDG)
             call en%op2_pdg  (en%aval, en%arg1, en%arg2)
          case (V_SEV)
             if (associated (en%arg0)) then
                call en%op2_sev (en%pval, en%arg1, en%arg2, en%arg0)
             else
                call en%op2_sev (en%pval, en%arg1, en%arg2)
             end if
          case (V_STR)
             call en%op2_str (en%sval, en%arg1, en%arg2)
          end select
       end if
    case (EN_BLOCK)
       if (associated (en%arg1) .and. associated (en%arg0)) then
          call eval_node_evaluate (en%arg1)
          call eval_node_evaluate (en%arg0)
          en%value_is_known = en%arg0%value_is_known
       else
          en%value_is_known = .false.
       end if
       if (en%value_is_known) then
          select case (en%result_type)
          case (V_LOG);  en%lval = en%arg0%lval
          case (V_INT);  en%ival = en%arg0%ival
          case (V_REAL); en%rval = en%arg0%rval
          case (V_CMPLX); en%cval = en%arg0%cval
          case (V_PDG);  en%aval = en%arg0%aval
          case (V_SEV);  en%pval = en%arg0%pval
          case (V_STR);  en%sval = en%arg0%sval
          end select
       end if
    case (EN_CONDITIONAL)
       if (associated (en%arg0)) then
          call eval_node_evaluate (en%arg0)
          en%value_is_known = en%arg0%value_is_known
       else
          en%value_is_known = .false.
       end if
       if (en%arg0%value_is_known) then
          if (en%arg0%lval) then
             call eval_node_evaluate (en%arg1)
             en%value_is_known = en%arg1%value_is_known
             if (en%value_is_known) then
                select case (en%result_type)
                case (V_LOG);  en%lval = en%arg1%lval
                case (V_INT);  en%ival = en%arg1%ival
                case (V_REAL); en%rval = en%arg1%rval
                case (V_CMPLX); en%cval = en%arg1%cval
                case (V_PDG);  en%aval = en%arg1%aval
                case (V_SEV);  en%pval = en%arg1%pval
                case (V_STR);  en%sval = en%arg1%sval
                end select
             end if
          else
             call eval_node_evaluate (en%arg2)
             en%value_is_known = en%arg2%value_is_known
             if (en%value_is_known) then
                select case (en%result_type)
                case (V_LOG);  en%lval = en%arg2%lval
                case (V_INT);  en%ival = en%arg2%ival
                case (V_REAL); en%rval = en%arg2%rval
                case (V_CMPLX); en%cval = en%arg2%cval
                case (V_PDG);  en%aval = en%arg2%aval
                case (V_SEV);  en%pval = en%arg2%pval
                case (V_STR);  en%sval = en%arg2%sval
                end select
             end if
          end if
       end if
    case (EN_RECORD_CMD)
       exist = .true.
       en%lval = .false.
       call eval_node_evaluate (en%arg0)
       if (en%arg0%value_is_known) then
          if (associated (en%arg1)) then
             call eval_node_evaluate (en%arg1)
             if (en%arg1%value_is_known) then
                if (associated (en%arg2)) then
                   call eval_node_evaluate (en%arg2)
                   if (en%arg2%value_is_known) then
                      if (associated (en%arg3)) then
                         call eval_node_evaluate (en%arg3)
                         if (en%arg3%value_is_known) then
                            if (associated (en%arg4)) then
                               call eval_node_evaluate (en%arg4)
                               if (en%arg4%value_is_known) then
                                  if (associated (en%rval)) then
                                     call analysis_record_data (en%arg0%sval, &
                                          en%arg1%rval, en%arg2%rval, &
                                          en%arg3%rval, en%arg4%rval, &
                                          weight=en%rval, exist=exist, &
                                          success=en%lval)
                                  else
                                     call analysis_record_data (en%arg0%sval, &
                                          en%arg1%rval, en%arg2%rval, &
                                          en%arg3%rval, en%arg4%rval, &
                                          exist=exist, success=en%lval)
                                  end if
                               end if
                            else
                               if (associated (en%rval)) then
                                  call analysis_record_data (en%arg0%sval, &
                                       en%arg1%rval, en%arg2%rval, &
                                       en%arg3%rval, &
                                       weight=en%rval, exist=exist, &
                                       success=en%lval)
                               else
                                  call analysis_record_data (en%arg0%sval, &
                                       en%arg1%rval, en%arg2%rval, &
                                       en%arg3%rval, &
                                       exist=exist, success=en%lval)
                               end if
                            end if
                         end if
                      else
                         if (associated (en%rval)) then
                            call analysis_record_data (en%arg0%sval, &
                                 en%arg1%rval, en%arg2%rval, &
                                 weight=en%rval, exist=exist, &
                                 success=en%lval)
                         else
                            call analysis_record_data (en%arg0%sval, &
                                 en%arg1%rval, en%arg2%rval, &
                                 exist=exist, success=en%lval)
                         end if
                      end if
                   end if
                else
                   if (associated (en%rval)) then
                      call analysis_record_data (en%arg0%sval, &
                           en%arg1%rval, &
                           weight=en%rval, exist=exist, success=en%lval)
                   else
                      call analysis_record_data (en%arg0%sval, &
                           en%arg1%rval, &
                           exist=exist, success=en%lval)
                   end if
                end if
             end if
          else
             if (associated (en%rval)) then
                call analysis_record_data (en%arg0%sval, 1._default, &
                     weight=en%rval, exist=exist, success=en%lval)
             else
                call analysis_record_data (en%arg0%sval, 1._default, &
                     exist=exist, success=en%lval)
             end if
          end if
          if (.not. exist) then
             call msg_error ("Analysis object '" // char (en%arg0%sval) &
                  // "' is undefined")
             en%arg0%value_is_known = .false.
          end if
       end if
    case (EN_OBS1_INT)
       en%ival = en%obs1_int (en%prt1)
       en%value_is_known = .true.
    case (EN_OBS2_INT)
       en%ival = en%obs2_int (en%prt1, en%prt2)
       en%value_is_known = .true.
    case (EN_OBSEV_INT)
       en%ival = en%obsev_int (en%pval)
       en%value_is_known = .true.
    case (EN_OBS1_REAL)
       en%rval = en%obs1_real (en%prt1)
       en%value_is_known = .true.
    case (EN_OBS2_REAL)
       en%rval = en%obs2_real (en%prt1, en%prt2)
       en%value_is_known = .true.
    case (EN_OBSEV_REAL)
       en%rval = en%obsev_real (en%pval)
       en%value_is_known = .true.
    case (EN_PRT_FUN_UNARY)
       call eval_node_evaluate (en%arg1)
       en%value_is_known = en%arg1%value_is_known
       if (en%value_is_known) then
          if (associated (en%arg0)) then
             en%arg0%index => en%index
             en%arg0%prt1 => en%prt1
             call en%op1_sev (en%pval, en%arg1, en%arg0)
          else
             call en%op1_sev (en%pval, en%arg1)
          end if
       end if
    case (EN_PRT_FUN_BINARY)
       call eval_node_evaluate (en%arg1)
       call eval_node_evaluate (en%arg2)
       en%value_is_known = &
            en%arg1%value_is_known .and. en%arg2%value_is_known
       if (en%value_is_known) then
          if (associated (en%arg0)) then
             en%arg0%index => en%index
             en%arg0%prt1 => en%prt1
             en%arg0%prt2 => en%prt2
             call en%op2_sev (en%pval, en%arg1, en%arg2, en%arg0)
          else
             call en%op2_sev (en%pval, en%arg1, en%arg2)
          end if
       end if
    case (EN_EVAL_FUN_UNARY)
       call eval_node_evaluate (en%arg1)
       en%value_is_known = en%arg1%pval%is_nonempty ()
       if (en%value_is_known) then
          en%arg0%index => en%index
          en%index = 1
          en%arg0%prt1 => en%prt1
          en%prt1 = en%arg1%pval%get_prt (1)
          call eval_node_evaluate (en%arg0)
          en%rval = en%arg0%rval
       end if
    case (EN_EVAL_FUN_BINARY)
       call eval_node_evaluate (en%arg1)
       call eval_node_evaluate (en%arg2)
       en%value_is_known = &
            en%arg1%pval%is_nonempty () .and. en%arg2%pval%is_nonempty ()
       if (en%value_is_known) then
          en%arg0%index => en%index
          en%arg0%prt1 => en%prt1
          en%arg0%prt2 => en%prt2
          en%index = 1
          call eval_pp (en%arg1, en%arg2, en%arg0, en%rval, en%value_is_known)
       end if
    case (EN_LOG_FUN_UNARY)
       call eval_node_evaluate (en%arg1)
       en%value_is_known = .true.
       if (en%value_is_known) then
          en%arg0%index => en%index
          en%arg0%prt1 => en%prt1
          en%lval = en%op1_cut (en%arg1, en%arg0)
       end if
    case (EN_LOG_FUN_BINARY)
       call eval_node_evaluate (en%arg1)
       call eval_node_evaluate (en%arg2)
       en%value_is_known = .true.
       if (en%value_is_known) then
          en%arg0%index => en%index
          en%arg0%prt1 => en%prt1
          en%arg0%prt2 => en%prt2
          en%lval = en%op2_cut (en%arg1, en%arg2, en%arg0)
       end if
    case (EN_INT_FUN_UNARY)
       call eval_node_evaluate (en%arg1)
       en%value_is_known = en%arg1%value_is_known
       if (en%value_is_known) then
          if (associated (en%arg0)) then
             en%arg0%index => en%index
             en%arg0%prt1 => en%prt1
             call en%op1_evi (en%ival, en%arg1, en%arg0)
          else
             call en%op1_evi (en%ival, en%arg1)
          end if
       end if
    case (EN_INT_FUN_BINARY)
       call eval_node_evaluate (en%arg1)
       call eval_node_evaluate (en%arg2)
       en%value_is_known = &
            en%arg1%value_is_known .and. &
            en%arg2%value_is_known
       if (en%value_is_known) then
          if (associated (en%arg0)) then
             en%arg0%index => en%index
             en%arg0%prt1 => en%prt1
             en%arg0%prt2 => en%prt2
             call en%op2_evi (en%ival, en%arg1, en%arg2, en%arg0)
          else
             call en%op2_evi (en%ival, en%arg1, en%arg2)
          end if
       end if
    case (EN_REAL_FUN_UNARY)
       call eval_node_evaluate (en%arg1)
       en%value_is_known = en%arg1%value_is_known
       if (en%value_is_known) then
          if (associated (en%arg0)) then
             en%arg0%index => en%index
             en%arg0%prt1 => en%prt1
             call en%op1_evr (en%rval, en%arg1, en%arg0)
          else
             call en%op1_evr (en%rval, en%arg1)
          end if
       end if
    case (EN_REAL_FUN_BINARY)
       call eval_node_evaluate (en%arg1)
       call eval_node_evaluate (en%arg2)
       en%value_is_known = &
            en%arg1%value_is_known .and. &
            en%arg2%value_is_known
       if (en%value_is_known) then
          if (associated (en%arg0)) then
             en%arg0%index => en%index
             en%arg0%prt1 => en%prt1
             en%arg0%prt2 => en%prt2
             call en%op2_evr (en%rval, en%arg1, en%arg2, en%arg0)
          else
             call en%op2_evr (en%rval, en%arg1, en%arg2)
          end if
       end if
     case (EN_REAL_FUN_CUM)
        call eval_node_evaluate (en%arg1)
        en%value_is_known = .true.
        if (en%value_is_known) then
           en%arg0%index => en%index
           en%arg0%prt1 => en%prt1
           en%rval = en%opcum_evr (en%arg1, en%arg0)
       end if
    case (EN_FORMAT_STR)
       if (associated (en%arg0)) then
          call eval_node_evaluate (en%arg0)
          en%value_is_known = en%arg0%value_is_known
       else
          en%value_is_known = .true.
       end if
       if (associated (en%arg1)) then
          call eval_node_evaluate (en%arg1)
          en%value_is_known = &
               en%value_is_known .and. en%arg1%value_is_known
          if (en%value_is_known) then
             call evaluate_sprintf (en%sval, en%ival, en%arg0, en%arg1)
          end if
       else
          if (en%value_is_known) then
             call evaluate_sprintf (en%sval, en%ival, en%arg0)
          end if
       end if
    end select
    if (debug2_active (D_MODEL_F)) then
       print *, "eval_node_evaluate"
       call eval_node_write (en)
    end if
  end subroutine eval_node_evaluate

  module subroutine eval_node_test_obs (node, var_list, var_name)
    class(eval_node_t), intent(inout) :: node
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: var_name
    procedure(obs_unary_int), pointer :: obs1_iptr
    type(prt_t), pointer :: p1
    call var_list%get_obs1_iptr (var_name, obs1_iptr, p1)
    call eval_node_init_obs1_int_ptr (node, var_name, obs1_iptr, p1)
  end subroutine eval_node_test_obs

  module subroutine syntax_expr_init ()
    type(ifile_t) :: ifile
    call define_expr_syntax (ifile, particles=.false., analysis=.false.)
    call syntax_init (syntax_expr, ifile)
    call ifile_final (ifile)
  end subroutine syntax_expr_init

  module subroutine syntax_pexpr_init ()
    type(ifile_t) :: ifile
    call define_expr_syntax (ifile, particles=.true., analysis=.false.)
    call syntax_init (syntax_pexpr, ifile)
    call ifile_final (ifile)
  end subroutine syntax_pexpr_init

  module subroutine syntax_expr_final ()
    call syntax_final (syntax_expr)
  end subroutine syntax_expr_final

  module subroutine syntax_pexpr_final ()
    call syntax_final (syntax_pexpr)
  end subroutine syntax_pexpr_final

  module subroutine syntax_pexpr_write (unit)
    integer, intent(in), optional :: unit
    call syntax_write (syntax_pexpr, unit)
  end subroutine syntax_pexpr_write

  module subroutine define_expr_syntax (ifile, particles, analysis)
    type(ifile_t), intent(inout) :: ifile
    logical, intent(in) :: particles, analysis
    type(string_t) :: numeric_pexpr
    type(string_t) :: var_plist, var_alias
    if (particles) then
       numeric_pexpr = " | numeric_pexpr"
       var_plist = " | var_plist"
       var_alias = " | var_alias"
    else
       numeric_pexpr = ""
       var_plist = ""
       var_alias = ""
    end if
    call ifile_append (ifile, "SEQ expr = subexpr addition*")
    call ifile_append (ifile, "ALT subexpr = addition | term")
    call ifile_append (ifile, "SEQ addition = plus_or_minus term")
    call ifile_append (ifile, "SEQ term = factor multiplication*")
    call ifile_append (ifile, "SEQ multiplication = times_or_over factor")
    call ifile_append (ifile, "SEQ factor = value exponentiation?")
    call ifile_append (ifile, "SEQ exponentiation = to_the value")
    call ifile_append (ifile, "ALT plus_or_minus = '+' | '-'")
    call ifile_append (ifile, "ALT times_or_over = '*' | '/'")
    call ifile_append (ifile, "ALT to_the = '^' | '**'")
    call ifile_append (ifile, "KEY '+'")
    call ifile_append (ifile, "KEY '-'")
    call ifile_append (ifile, "KEY '*'")
    call ifile_append (ifile, "KEY '/'")
    call ifile_append (ifile, "KEY '^'")
    call ifile_append (ifile, "KEY '**'")
    call ifile_append (ifile, "ALT value = signed_value | unsigned_value")
    call ifile_append (ifile, "SEQ signed_value = '-' unsigned_value")
    call ifile_append (ifile, "ALT unsigned_value = " // &
         "numeric_value | constant | variable | " // &
         "result | " // &
         "grouped_expr | block_expr | conditional_expr | " // &
         "unary_function | binary_function" // &
         numeric_pexpr)
    call ifile_append (ifile, "ALT numeric_value = integer_value | " &
         // "real_value | complex_value")
    call ifile_append (ifile, "SEQ integer_value = integer_literal unit_expr?")
    call ifile_append (ifile, "SEQ real_value = real_literal unit_expr?")
    call ifile_append (ifile, "SEQ complex_value = complex_literal unit_expr?")
    call ifile_append (ifile, "INT integer_literal")
    call ifile_append (ifile, "REA real_literal")
    call ifile_append (ifile, "COM complex_literal")
    call ifile_append (ifile, "SEQ unit_expr = unit unit_power?")
    call ifile_append (ifile, "ALT unit = " // &
         "TeV | GeV | MeV | keV | eV | meV | " // &
         "nbarn | pbarn | fbarn | abarn | " // &
         "rad | mrad | degree | '%'")
    call ifile_append (ifile, "KEY TeV")
    call ifile_append (ifile, "KEY GeV")
    call ifile_append (ifile, "KEY MeV")
    call ifile_append (ifile, "KEY keV")
    call ifile_append (ifile, "KEY eV")
    call ifile_append (ifile, "KEY meV")
    call ifile_append (ifile, "KEY nbarn")
    call ifile_append (ifile, "KEY pbarn")
    call ifile_append (ifile, "KEY fbarn")
    call ifile_append (ifile, "KEY abarn")
    call ifile_append (ifile, "KEY rad")
    call ifile_append (ifile, "KEY mrad")
    call ifile_append (ifile, "KEY degree")
    call ifile_append (ifile, "KEY '%'")
    call ifile_append (ifile, "SEQ unit_power = '^' frac_expr")
    call ifile_append (ifile, "ALT frac_expr = frac | grouped_frac")
    call ifile_append (ifile, "GRO grouped_frac = ( frac_expr )")
    call ifile_append (ifile, "SEQ frac = signed_int div?")
    call ifile_append (ifile, "ALT signed_int = " &
         // "neg_int | pos_int | integer_literal")
    call ifile_append (ifile, "SEQ neg_int = '-' integer_literal")
    call ifile_append (ifile, "SEQ pos_int = '+' integer_literal")
    call ifile_append (ifile, "SEQ div = '/' integer_literal")
    call ifile_append (ifile, "ALT constant = pi | I")
    call ifile_append (ifile, "KEY pi")
    call ifile_append (ifile, "KEY I")
    call ifile_append (ifile, "IDE variable")
    call ifile_append (ifile, "SEQ result = result_key result_arg")
    call ifile_append (ifile, "ALT result_key = " // &
         "num_id | integral | error")
    call ifile_append (ifile, "KEY num_id")
    call ifile_append (ifile, "KEY integral")
    call ifile_append (ifile, "KEY error")
    call ifile_append (ifile, "GRO result_arg = ( process_id )")
    call ifile_append (ifile, "IDE process_id")
    call ifile_append (ifile, "SEQ unary_function = fun_unary function_arg1")
    call ifile_append (ifile, "SEQ binary_function = fun_binary function_arg2")
    call ifile_append (ifile, "ALT fun_unary = " // &
         "complex | real | int | nint | floor | ceiling | abs | conjg | sgn | " // &
         "sqrt | exp | log | log10 | " // &
         "sin | cos | tan | asin | acos | atan | " // &
         "sinh | cosh | tanh | asinh | acosh | atanh")
    call ifile_append (ifile, "KEY complex")
    call ifile_append (ifile, "KEY real")
    call ifile_append (ifile, "KEY int")
    call ifile_append (ifile, "KEY nint")
    call ifile_append (ifile, "KEY floor")
    call ifile_append (ifile, "KEY ceiling")
    call ifile_append (ifile, "KEY abs")
    call ifile_append (ifile, "KEY conjg")
    call ifile_append (ifile, "KEY sgn")
    call ifile_append (ifile, "KEY sqrt")
    call ifile_append (ifile, "KEY exp")
    call ifile_append (ifile, "KEY log")
    call ifile_append (ifile, "KEY log10")
    call ifile_append (ifile, "KEY sin")
    call ifile_append (ifile, "KEY cos")
    call ifile_append (ifile, "KEY tan")
    call ifile_append (ifile, "KEY asin")
    call ifile_append (ifile, "KEY acos")
    call ifile_append (ifile, "KEY atan")
    call ifile_append (ifile, "KEY sinh")
    call ifile_append (ifile, "KEY cosh")
    call ifile_append (ifile, "KEY tanh")
    call ifile_append (ifile, "KEY asinh")
    call ifile_append (ifile, "KEY acosh")
    call ifile_append (ifile, "KEY atanh")
    call ifile_append (ifile, "ALT fun_binary = max | min | mod | modulo")
    call ifile_append (ifile, "KEY max")
    call ifile_append (ifile, "KEY min")
    call ifile_append (ifile, "KEY mod")
    call ifile_append (ifile, "KEY modulo")
    call ifile_append (ifile, "ARG function_arg1 = ( expr )")
    call ifile_append (ifile, "ARG function_arg2 = ( expr, expr )")
    call ifile_append (ifile, "GRO grouped_expr = ( expr )")
    call ifile_append (ifile, "SEQ block_expr = let var_spec in expr")
    call ifile_append (ifile, "KEY let")
    call ifile_append (ifile, "ALT var_spec = " // &
         "var_num | var_int | var_real | var_complex | " // &
         "var_logical" // var_plist // var_alias // " | var_string")
    call ifile_append (ifile, "SEQ var_num = var_name '=' expr")
    call ifile_append (ifile, "SEQ var_int = int var_name '=' expr")
    call ifile_append (ifile, "SEQ var_real = real var_name '=' expr")
    call ifile_append (ifile, "SEQ var_complex = complex var_name '=' complex_expr")
    call ifile_append (ifile, "ALT complex_expr = " // &
         "cexpr_real | cexpr_complex")
    call ifile_append (ifile, "ARG cexpr_complex = ( expr, expr )")
    call ifile_append (ifile, "SEQ cexpr_real = expr")
    call ifile_append (ifile, "IDE var_name")
    call ifile_append (ifile, "KEY '='")
    call ifile_append (ifile, "KEY in")
    call ifile_append (ifile, "SEQ conditional_expr = " // &
         "if lexpr then expr maybe_elsif_expr maybe_else_expr endif")
    call ifile_append (ifile, "SEQ maybe_elsif_expr = elsif_expr*")
    call ifile_append (ifile, "SEQ maybe_else_expr = else_expr?")
    call ifile_append (ifile, "SEQ elsif_expr = elsif lexpr then expr")
    call ifile_append (ifile, "SEQ else_expr = else expr")
    call ifile_append (ifile, "KEY if")
    call ifile_append (ifile, "KEY then")
    call ifile_append (ifile, "KEY elsif")
    call ifile_append (ifile, "KEY else")
    call ifile_append (ifile, "KEY endif")
    call define_lexpr_syntax (ifile, particles, analysis)
    call define_sexpr_syntax (ifile)
    if (particles) then
       call define_pexpr_syntax (ifile)
       call define_cexpr_syntax (ifile)
       call define_var_plist_syntax (ifile)
       call define_var_alias_syntax (ifile)
       call define_numeric_pexpr_syntax (ifile)
       call define_logical_pexpr_syntax (ifile)
    end if

  end subroutine define_expr_syntax

  subroutine define_lexpr_syntax (ifile, particles, analysis)
    type(ifile_t), intent(inout) :: ifile
    logical, intent(in) :: particles, analysis
    type(string_t) :: logical_pexpr, record_cmd
    if (particles) then
       logical_pexpr = " | logical_pexpr"
    else
       logical_pexpr = ""
    end if
    if (analysis) then
       record_cmd = " | record_cmd"
    else
       record_cmd = ""
    end if
    call ifile_append (ifile, "SEQ lexpr = lsinglet lsequel*")
    call ifile_append (ifile, "SEQ lsequel = ';' lsinglet")
    call ifile_append (ifile, "SEQ lsinglet = lterm alternative*")
    call ifile_append (ifile, "SEQ alternative = or lterm")
    call ifile_append (ifile, "SEQ lterm = lvalue coincidence*")
    call ifile_append (ifile, "SEQ coincidence = and lvalue")
    call ifile_append (ifile, "KEY ';'")
    call ifile_append (ifile, "KEY or")
    call ifile_append (ifile, "KEY and")
    call ifile_append (ifile, "ALT lvalue = " // &
         "true | false | lvariable | negation | " // &
         "grouped_lexpr | block_lexpr | conditional_lexpr | " // &
         "compared_expr | compared_sexpr" // &
         logical_pexpr //  record_cmd)
    call ifile_append (ifile, "KEY true")
    call ifile_append (ifile, "KEY false")
    call ifile_append (ifile, "SEQ lvariable = '?' alt_lvariable")
    call ifile_append (ifile, "KEY '?'")
    call ifile_append (ifile, "ALT alt_lvariable = variable | grouped_lexpr")
    call ifile_append (ifile, "SEQ negation = not lvalue")
    call ifile_append (ifile, "KEY not")
    call ifile_append (ifile, "GRO grouped_lexpr = ( lexpr )")
    call ifile_append (ifile, "SEQ block_lexpr = let var_spec in lexpr")
    call ifile_append (ifile, "ALT var_logical = " // &
         "var_logical_new | var_logical_spec")
    call ifile_append (ifile, "SEQ var_logical_new = logical var_logical_spec")
    call ifile_append (ifile, "KEY logical")
    call ifile_append (ifile, "SEQ var_logical_spec = '?' var_name = lexpr")
    call ifile_append (ifile, "SEQ conditional_lexpr = " // &
         "if lexpr then lexpr maybe_elsif_lexpr maybe_else_lexpr endif")
    call ifile_append (ifile, "SEQ maybe_elsif_lexpr = elsif_lexpr*")
    call ifile_append (ifile, "SEQ maybe_else_lexpr = else_lexpr?")
    call ifile_append (ifile, "SEQ elsif_lexpr = elsif lexpr then lexpr")
    call ifile_append (ifile, "SEQ else_lexpr = else lexpr")
    call ifile_append (ifile, "SEQ compared_expr = expr comparison+")
    call ifile_append (ifile, "SEQ comparison = compare expr")
    call ifile_append (ifile, "ALT compare = " // &
         "'<' | '>' | '<=' | '>=' | '==' | '<>'")
    call ifile_append (ifile, "KEY '<'")
    call ifile_append (ifile, "KEY '>'")
    call ifile_append (ifile, "KEY '<='")
    call ifile_append (ifile, "KEY '>='")
    call ifile_append (ifile, "KEY '=='")
    call ifile_append (ifile, "KEY '<>'")
    call ifile_append (ifile, "SEQ compared_sexpr = sexpr str_comparison+")
    call ifile_append (ifile, "SEQ str_comparison = str_compare sexpr")
    call ifile_append (ifile, "ALT str_compare = '==' | '<>'")
    if (analysis) then
       call ifile_append (ifile, "SEQ record_cmd = " // &
            "record_key analysis_tag record_arg?")
       call ifile_append (ifile, "ALT record_key = " // &
            "record | record_unweighted | record_excess")
       call ifile_append (ifile, "KEY record")
       call ifile_append (ifile, "KEY record_unweighted")
       call ifile_append (ifile, "KEY record_excess")
       call ifile_append (ifile, "ALT analysis_tag = analysis_id | sexpr")
       call ifile_append (ifile, "IDE analysis_id")
       call ifile_append (ifile, "ARG record_arg = ( expr+ )")
    end if
  end subroutine define_lexpr_syntax

  subroutine define_sexpr_syntax (ifile)
    type(ifile_t), intent(inout) :: ifile
    call ifile_append (ifile, "SEQ sexpr = svalue str_concatenation*")
    call ifile_append (ifile, "SEQ str_concatenation = '&' svalue")
    call ifile_append (ifile, "KEY '&'")
    call ifile_append (ifile, "ALT svalue = " // &
         "grouped_sexpr | block_sexpr | conditional_sexpr | " // &
         "svariable | string_function | string_literal")
    call ifile_append (ifile, "GRO grouped_sexpr = ( sexpr )")
    call ifile_append (ifile, "SEQ block_sexpr = let var_spec in sexpr")
    call ifile_append (ifile, "SEQ conditional_sexpr = " // &
         "if lexpr then sexpr maybe_elsif_sexpr maybe_else_sexpr endif")
    call ifile_append (ifile, "SEQ maybe_elsif_sexpr = elsif_sexpr*")
    call ifile_append (ifile, "SEQ maybe_else_sexpr = else_sexpr?")
    call ifile_append (ifile, "SEQ elsif_sexpr = elsif lexpr then sexpr")
    call ifile_append (ifile, "SEQ else_sexpr = else sexpr")
    call ifile_append (ifile, "SEQ svariable = '$' alt_svariable")
    call ifile_append (ifile, "KEY '$'")
    call ifile_append (ifile, "ALT alt_svariable = variable | grouped_sexpr")
    call ifile_append (ifile, "ALT var_string = " // &
         "var_string_new | var_string_spec")
    call ifile_append (ifile, "SEQ var_string_new = string var_string_spec")
    call ifile_append (ifile, "KEY string")
    call ifile_append (ifile, "SEQ var_string_spec = '$' var_name = sexpr") ! $
    call ifile_append (ifile, "ALT string_function = sprintf_fun")
    call ifile_append (ifile, "SEQ sprintf_fun = sprintf_clause sprintf_args?")
    call ifile_append (ifile, "SEQ sprintf_clause = sprintf sexpr")
    call ifile_append (ifile, "KEY sprintf")
    call ifile_append (ifile, "ARG sprintf_args = ( sprintf_arg* )")
    call ifile_append (ifile, "ALT sprintf_arg = " &
         // "lvariable | svariable | expr")
    call ifile_append (ifile, "QUO string_literal = '""'...'""'")
  end subroutine define_sexpr_syntax

  subroutine define_pexpr_syntax (ifile)
    type(ifile_t), intent(inout) :: ifile
    call ifile_append (ifile, "SEQ pexpr = pterm pconcatenation*")
    call ifile_append (ifile, "SEQ pconcatenation = '&' pterm")
    ! call ifile_append (ifile, "KEY '&'")   !!! (Key exists already)
    call ifile_append (ifile, "SEQ pterm = pvalue pcombination*")
    call ifile_append (ifile, "SEQ pcombination = '+' pvalue")
    ! call ifile_append (ifile, "KEY '+'")   !!! (Key exists already)
    call ifile_append (ifile, "ALT pvalue = " // &
         "pexpr_src | pvariable | " // &
         "grouped_pexpr | block_pexpr | conditional_pexpr | " // &
         "prt_function")
    call ifile_append (ifile, "SEQ pexpr_src = prefix_cexpr")
    call ifile_append (ifile, "ALT prefix_cexpr = " // &
         "beam_prt | incoming_prt | outgoing_prt | unspecified_prt")
    call ifile_append (ifile, "SEQ beam_prt = beam cexpr")
    call ifile_append (ifile, "KEY beam")
    call ifile_append (ifile, "SEQ incoming_prt = incoming cexpr")
    call ifile_append (ifile, "KEY incoming")
    call ifile_append (ifile, "SEQ outgoing_prt = outgoing cexpr")
    call ifile_append (ifile, "KEY outgoing")
    call ifile_append (ifile, "SEQ unspecified_prt = cexpr")
    call ifile_append (ifile, "SEQ pvariable = '@' alt_pvariable")
    call ifile_append (ifile, "KEY '@'")
    call ifile_append (ifile, "ALT alt_pvariable = variable | grouped_pexpr")
    call ifile_append (ifile, "GRO grouped_pexpr = '[' pexpr ']'")
    call ifile_append (ifile, "SEQ block_pexpr = let var_spec in pexpr")
    call ifile_append (ifile, "SEQ conditional_pexpr = " // &
         "if lexpr then pexpr maybe_elsif_pexpr maybe_else_pexpr endif")
    call ifile_append (ifile, "SEQ maybe_elsif_pexpr = elsif_pexpr*")
    call ifile_append (ifile, "SEQ maybe_else_pexpr = else_pexpr?")
    call ifile_append (ifile, "SEQ elsif_pexpr = elsif lexpr then pexpr")
    call ifile_append (ifile, "SEQ else_pexpr = else pexpr")
    call ifile_append (ifile, "ALT prt_function = " // &
         "join_fun | combine_fun | collect_fun | cluster_fun | " // &
         "photon_reco_fun | " // &
         "select_fun | extract_fun | sort_fun | " // &
         "select_b_jet_fun | select_non_b_jet_fun | " // &
         "select_c_jet_fun | select_light_jet_fun")
    call ifile_append (ifile, "SEQ join_fun = join_clause pargs2")
    call ifile_append (ifile, "SEQ combine_fun = combine_clause pargs2")
    call ifile_append (ifile, "SEQ collect_fun = collect_clause pargs1")
    call ifile_append (ifile, "SEQ cluster_fun = cluster_clause pargs1")
    call ifile_append (ifile, "SEQ photon_reco_fun = photon_reco_clause pargs1")
    call ifile_append (ifile, "SEQ select_fun = select_clause pargs1")
    call ifile_append (ifile, "SEQ extract_fun = extract_clause pargs1")
    call ifile_append (ifile, "SEQ sort_fun = sort_clause pargs1")
    call ifile_append (ifile, "SEQ select_b_jet_fun = " // &
          "select_b_jet_clause pargs1")
    call ifile_append (ifile, "SEQ select_non_b_jet_fun = " // &
          "select_non_b_jet_clause pargs1")
    call ifile_append (ifile, "SEQ select_c_jet_fun = " // &
          "select_c_jet_clause pargs1")
    call ifile_append (ifile, "SEQ select_light_jet_fun = " // &
          "select_light_jet_clause pargs1")
    call ifile_append (ifile, "SEQ join_clause = join condition?")
    call ifile_append (ifile, "SEQ combine_clause = combine condition?")
    call ifile_append (ifile, "SEQ collect_clause = collect condition?")
    call ifile_append (ifile, "SEQ cluster_clause = cluster condition?")
    call ifile_append (ifile, "SEQ photon_reco_clause = photon_recombination condition?")
    call ifile_append (ifile, "SEQ select_clause = select condition?")
    call ifile_append (ifile, "SEQ extract_clause = extract position?")
    call ifile_append (ifile, "SEQ sort_clause = sort criterion?")
    call ifile_append (ifile, "SEQ select_b_jet_clause = " // &
         "select_b_jet condition?")
    call ifile_append (ifile, "SEQ select_non_b_jet_clause = " // &
         "select_non_b_jet condition?")
    call ifile_append (ifile, "SEQ select_c_jet_clause = " // &
         "select_c_jet condition?")
    call ifile_append (ifile, "SEQ select_light_jet_clause = " // &
         "select_light_jet condition?")
    call ifile_append (ifile, "KEY join")
    call ifile_append (ifile, "KEY combine")
    call ifile_append (ifile, "KEY collect")
    call ifile_append (ifile, "KEY cluster")
    call ifile_append (ifile, "KEY photon_recombination")
    call ifile_append (ifile, "KEY select")
    call ifile_append (ifile, "SEQ condition = if lexpr")
    call ifile_append (ifile, "KEY extract")
    call ifile_append (ifile, "SEQ position = index expr")
    call ifile_append (ifile, "KEY sort")
    call ifile_append (ifile, "KEY select_b_jet")
    call ifile_append (ifile, "KEY select_non_b_jet")
    call ifile_append (ifile, "KEY select_c_jet")
    call ifile_append (ifile, "KEY select_light_jet")
    call ifile_append (ifile, "SEQ criterion = by expr")
    call ifile_append (ifile, "KEY index")
    call ifile_append (ifile, "KEY by")
    call ifile_append (ifile, "ARG pargs2 = '[' pexpr, pexpr ']'")
    call ifile_append (ifile, "ARG pargs1 = '[' pexpr, pexpr? ']'")
  end subroutine define_pexpr_syntax

  subroutine define_cexpr_syntax (ifile)
    type(ifile_t), intent(inout) :: ifile
    call ifile_append (ifile, "SEQ cexpr = avalue concatenation*")
    call ifile_append (ifile, "SEQ concatenation = ':' avalue")
    call ifile_append (ifile, "KEY ':'")
    call ifile_append (ifile, "ALT avalue = " // &
         "grouped_cexpr | block_cexpr | conditional_cexpr | " // &
         "variable | pdg_code | prt_name")
    call ifile_append (ifile, "GRO grouped_cexpr = ( cexpr )")
    call ifile_append (ifile, "SEQ block_cexpr = let var_spec in cexpr")
    call ifile_append (ifile, "SEQ conditional_cexpr = " // &
         "if lexpr then cexpr maybe_elsif_cexpr maybe_else_cexpr endif")
    call ifile_append (ifile, "SEQ maybe_elsif_cexpr = elsif_cexpr*")
    call ifile_append (ifile, "SEQ maybe_else_cexpr = else_cexpr?")
    call ifile_append (ifile, "SEQ elsif_cexpr = elsif lexpr then cexpr")
    call ifile_append (ifile, "SEQ else_cexpr = else cexpr")
    call ifile_append (ifile, "SEQ pdg_code = pdg pdg_arg")
    call ifile_append (ifile, "KEY pdg")
    call ifile_append (ifile, "ARG pdg_arg = ( expr )")
    call ifile_append (ifile, "QUO prt_name = '""'...'""'")
  end subroutine define_cexpr_syntax

  subroutine define_var_plist_syntax (ifile)
    type(ifile_t), intent(inout) :: ifile
    call ifile_append (ifile, "ALT var_plist = var_plist_new | var_plist_spec")
    call ifile_append (ifile, "SEQ var_plist_new = subevt var_plist_spec")
    call ifile_append (ifile, "KEY subevt")
    call ifile_append (ifile, "SEQ var_plist_spec = '@' var_name '=' pexpr")
  end subroutine define_var_plist_syntax

  subroutine define_var_alias_syntax (ifile)
    type(ifile_t), intent(inout) :: ifile
    call ifile_append (ifile, "SEQ var_alias = alias var_name '=' cexpr")
    call ifile_append (ifile, "KEY alias")
  end subroutine define_var_alias_syntax

  subroutine define_numeric_pexpr_syntax (ifile)
    type(ifile_t), intent(inout) :: ifile
    call ifile_append (ifile, "ALT numeric_pexpr = " &
         // "eval_fun | count_fun | sum_fun | " &
         // "prod_fun")
    call ifile_append (ifile, "SEQ eval_fun = eval expr pargs1")
    call ifile_append (ifile, "SEQ count_fun = count_clause pargs1")
    call ifile_append (ifile, "SEQ count_clause = count condition?")
    call ifile_append (ifile, "SEQ sum_fun = sum expr pargs1")
    call ifile_append (ifile, "SEQ prod_fun = prod expr pargs1")
    call ifile_append (ifile, "KEY eval")
    call ifile_append (ifile, "KEY count")
    call ifile_append (ifile, "KEY sum")
    call ifile_append (ifile, "KEY prod")
  end subroutine define_numeric_pexpr_syntax

  subroutine define_logical_pexpr_syntax (ifile)
    type(ifile_t), intent(inout) :: ifile
    call ifile_append (ifile, "ALT logical_pexpr = " // &
         "all_fun | any_fun | no_fun | " // &
         "photon_isolation_fun")
    call ifile_append (ifile, "SEQ all_fun = all lexpr pargs1")
    call ifile_append (ifile, "SEQ any_fun = any lexpr pargs1")
    call ifile_append (ifile, "SEQ no_fun = no lexpr pargs1")
    call ifile_append (ifile, "SEQ photon_isolation_fun = " // &
         "photon_isolation_clause pargs2")
    call ifile_append (ifile, "SEQ photon_isolation_clause = " // &
         "photon_isolation condition?")
    call ifile_append (ifile, "KEY all")
    call ifile_append (ifile, "KEY any")
    call ifile_append (ifile, "KEY no")
    call ifile_append (ifile, "KEY photon_isolation")
  end subroutine define_logical_pexpr_syntax

  subroutine lexer_init_eval_tree (lexer, particles)
    type(lexer_t), intent(out) :: lexer
    logical, intent(in) :: particles
    type(keyword_list_t), pointer :: keyword_list
    if (particles) then
       keyword_list => syntax_get_keyword_list_ptr (syntax_pexpr)
    else
       keyword_list => syntax_get_keyword_list_ptr (syntax_expr)
    end if
    call lexer_init (lexer, &
         comment_chars = "#!", &
         quote_chars = '"', &
         quote_match = '"', &
         single_chars = "()[],;:&%?$@", &
         special_class = [ "+-*/^", "<>=~ " ] , &
         keyword_list = keyword_list)
  end subroutine lexer_init_eval_tree

  module subroutine parse_tree_init_expr (parse_tree, stream, particles)
    type(parse_tree_t), intent(out) :: parse_tree
    type(stream_t), intent(inout), target :: stream
    logical, intent(in) :: particles
    type(lexer_t) :: lexer
    call lexer_init_eval_tree (lexer, particles)
    call lexer_assign_stream (lexer, stream)
    if (particles) then
       call parse_tree_init &
            (parse_tree, syntax_pexpr, lexer, var_str ("expr"))
    else
       call parse_tree_init &
            (parse_tree, syntax_expr, lexer, var_str ("expr"))
    end if
    call lexer_final (lexer)
  end subroutine parse_tree_init_expr

  module subroutine parse_tree_init_lexpr (parse_tree, stream, particles)
    type(parse_tree_t), intent(out) :: parse_tree
    type(stream_t), intent(inout), target :: stream
    logical, intent(in) :: particles
    type(lexer_t) :: lexer
    call lexer_init_eval_tree (lexer, particles)
    call lexer_assign_stream (lexer, stream)
    if (particles) then
       call parse_tree_init &
            (parse_tree, syntax_pexpr, lexer, var_str ("lexpr"))
    else
       call parse_tree_init &
            (parse_tree, syntax_expr, lexer, var_str ("lexpr"))
    end if
    call lexer_final (lexer)
  end subroutine parse_tree_init_lexpr

  module subroutine parse_tree_init_pexpr (parse_tree, stream)
    type(parse_tree_t), intent(out) :: parse_tree
    type(stream_t), intent(inout), target :: stream
    type(lexer_t) :: lexer
    call lexer_init_eval_tree (lexer, .true.)
    call lexer_assign_stream (lexer, stream)
    call parse_tree_init &
         (parse_tree, syntax_pexpr, lexer, var_str ("pexpr"))
    call lexer_final (lexer)
  end subroutine parse_tree_init_pexpr

  module subroutine parse_tree_init_cexpr (parse_tree, stream)
    type(parse_tree_t), intent(out) :: parse_tree
    type(stream_t), intent(inout), target :: stream
    type(lexer_t) :: lexer
    call lexer_init_eval_tree (lexer, .true.)
    call lexer_assign_stream (lexer, stream)
    call parse_tree_init &
         (parse_tree, syntax_pexpr, lexer, var_str ("cexpr"))
    call lexer_final (lexer)
  end subroutine parse_tree_init_cexpr

  module subroutine parse_tree_init_sexpr (parse_tree, stream, particles)
    type(parse_tree_t), intent(out) :: parse_tree
    type(stream_t), intent(inout), target :: stream
    logical, intent(in) :: particles
    type(lexer_t) :: lexer
    call lexer_init_eval_tree (lexer, particles)
    call lexer_assign_stream (lexer, stream)
    if (particles) then
       call parse_tree_init &
            (parse_tree, syntax_pexpr, lexer, var_str ("sexpr"))
    else
       call parse_tree_init &
            (parse_tree, syntax_expr, lexer, var_str ("sexpr"))
    end if
    call lexer_final (lexer)
  end subroutine parse_tree_init_sexpr

  module subroutine eval_tree_init_stream &
       (eval_tree, stream, var_list, subevt, result_type)
    class(eval_tree_t), intent(out), target :: eval_tree
    type(stream_t), intent(inout), target :: stream
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), target, optional :: subevt
    integer, intent(in), optional :: result_type
    type(parse_tree_t) :: parse_tree
    type(parse_node_t), pointer :: nd_root
    integer :: type
    type = V_REAL;  if (present (result_type))  type = result_type
    select case (type)
    case (V_INT, V_REAL, V_CMPLX)
       call parse_tree_init_expr (parse_tree, stream, present (subevt))
    case (V_LOG)
       call parse_tree_init_lexpr (parse_tree, stream, present (subevt))
    case (V_SEV)
       call parse_tree_init_pexpr (parse_tree, stream)
    case (V_PDG)
       call parse_tree_init_cexpr (parse_tree, stream)
    case (V_STR)
       call parse_tree_init_sexpr (parse_tree, stream, present (subevt))
    end select
    nd_root => parse_tree%get_root_ptr ()
    if (associated (nd_root)) then
       select case (type)
       case (V_INT, V_REAL, V_CMPLX)
          call eval_tree_init_expr (eval_tree, nd_root, var_list, subevt)
       case (V_LOG)
          call eval_tree_init_lexpr (eval_tree, nd_root, var_list, subevt)
       case (V_SEV)
          call eval_tree_init_pexpr (eval_tree, nd_root, var_list, subevt)
       case (V_PDG)
          call eval_tree_init_cexpr (eval_tree, nd_root, var_list, subevt)
       case (V_STR)
          call eval_tree_init_sexpr (eval_tree, nd_root, var_list, subevt)
       end select
    end if
    call parse_tree_final (parse_tree)
  end subroutine eval_tree_init_stream

  module subroutine eval_tree_init_expr &
      (expr, parse_node, var_list, subevt)
    class(eval_tree_t), intent(out), target :: expr
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    call eval_tree_link_var_list (expr, var_list)
    if (present (subevt))  call eval_tree_set_subevt (expr, subevt)
    call eval_node_compile_expr &
         (expr%root, parse_node, expr%var_list)
  end subroutine eval_tree_init_expr

  module subroutine eval_tree_init_lexpr &
      (expr, parse_node, var_list, subevt)
    class(eval_tree_t), intent(out), target :: expr
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    call eval_tree_link_var_list (expr, var_list)
    if (present (subevt))  call eval_tree_set_subevt (expr, subevt)
    call eval_node_compile_lexpr &
         (expr%root, parse_node, expr%var_list)
  end subroutine eval_tree_init_lexpr

  module subroutine eval_tree_init_pexpr &
      (expr, parse_node, var_list, subevt)
    class(eval_tree_t), intent(out), target :: expr
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    call eval_tree_link_var_list (expr, var_list)
    if (present (subevt))  call eval_tree_set_subevt (expr, subevt)
    call eval_node_compile_pexpr &
         (expr%root, parse_node, expr%var_list)
  end subroutine eval_tree_init_pexpr

  module subroutine eval_tree_init_cexpr &
      (expr, parse_node, var_list, subevt)
    class(eval_tree_t), intent(out), target :: expr
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    call eval_tree_link_var_list (expr, var_list)
    if (present (subevt))  call eval_tree_set_subevt (expr, subevt)
    call eval_node_compile_cexpr &
         (expr%root, parse_node, expr%var_list)
  end subroutine eval_tree_init_cexpr

  module subroutine eval_tree_init_sexpr &
      (expr, parse_node, var_list, subevt)
    class(eval_tree_t), intent(out), target :: expr
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    call eval_tree_link_var_list (expr, var_list)
    if (present (subevt))  call eval_tree_set_subevt (expr, subevt)
    call eval_node_compile_sexpr &
         (expr%root, parse_node, expr%var_list)
  end subroutine eval_tree_init_sexpr

  module subroutine eval_tree_setup_expr (expr, vars)
    class(eval_tree_t), intent(inout), target :: expr
    class(vars_t), intent(in), target :: vars
    call eval_tree_link_var_list (expr, vars)
    call eval_node_compile_expr (expr%root, expr%pn, expr%var_list)
  end subroutine eval_tree_setup_expr

  module subroutine eval_tree_setup_lexpr (expr, vars)
    class(eval_tree_t), intent(inout), target :: expr
    class(vars_t), intent(in), target :: vars
    call eval_tree_link_var_list (expr, vars)
    call eval_node_compile_lexpr (expr%root, expr%pn, expr%var_list)
  end subroutine eval_tree_setup_lexpr

  module subroutine eval_tree_setup_pexpr (expr, vars)
    class(eval_tree_t), intent(inout), target :: expr
    class(vars_t), intent(in), target :: vars
    call eval_tree_link_var_list (expr, vars)
    call eval_node_compile_pexpr (expr%root, expr%pn, expr%var_list)
  end subroutine eval_tree_setup_pexpr

  module subroutine eval_tree_setup_cexpr (expr, vars)
    class(eval_tree_t), intent(inout), target :: expr
    class(vars_t), intent(in), target :: vars
    call eval_tree_link_var_list (expr, vars)
    call eval_node_compile_cexpr (expr%root, expr%pn, expr%var_list)
  end subroutine eval_tree_setup_cexpr

  module subroutine eval_tree_setup_sexpr (expr, vars)
    class(eval_tree_t), intent(inout), target :: expr
    class(vars_t), intent(in), target :: vars
    call eval_tree_link_var_list (expr, vars)
    call eval_node_compile_sexpr (expr%root, expr%pn, expr%var_list)
  end subroutine eval_tree_setup_sexpr

  module subroutine eval_tree_init_numeric_value (eval_tree, parse_node)
    class(eval_tree_t), intent(out), target :: eval_tree
    type(parse_node_t), intent(in), target :: parse_node
    call eval_node_compile_numeric_value (eval_tree%root, parse_node)
  end subroutine eval_tree_init_numeric_value

  subroutine eval_tree_link_var_list (eval_tree, vars)
    type(eval_tree_t), intent(inout), target :: eval_tree
    class(vars_t), intent(in), target :: vars
    call eval_tree%var_list%link (vars)
  end subroutine eval_tree_link_var_list

  subroutine eval_tree_set_subevt (eval_tree, subevt)
    type(eval_tree_t), intent(inout), target :: eval_tree
    type(subevt_t), intent(in), target :: subevt
    logical, save, target :: known = .true.
    call eval_tree%var_list%append_subevt_ptr &
         (var_str ("@evt"), subevt, known, intrinsic=.true.)
  end subroutine eval_tree_set_subevt

  module subroutine eval_tree_final (expr)
    class(eval_tree_t), intent(inout) :: expr
    call expr%var_list%final ()
    if (associated (expr%root)) then
       call eval_node_final_rec (expr%root)
       deallocate (expr%root)
    end if
  end subroutine eval_tree_final

  module subroutine eval_tree_evaluate (expr)
    class(eval_tree_t), intent(inout) :: expr
    if (associated (expr%root)) then
       call eval_node_evaluate (expr%root)
    end if
  end subroutine eval_tree_evaluate

  function eval_tree_is_defined (eval_tree) result (flag)
    logical :: flag
    type(eval_tree_t), intent(in) :: eval_tree
    flag = associated (eval_tree%root)
  end function eval_tree_is_defined

  function eval_tree_is_constant (eval_tree) result (flag)
    logical :: flag
    type(eval_tree_t), intent(in) :: eval_tree
    if (associated (eval_tree%root)) then
       flag = eval_tree%root%type == EN_CONSTANT
    else
       flag = .false.
    end if
  end function eval_tree_is_constant

  subroutine eval_tree_convert_result (eval_tree, result_type)
    type(eval_tree_t), intent(inout) :: eval_tree
    integer, intent(in) :: result_type
    if (associated (eval_tree%root)) then
       call insert_conversion_node (eval_tree%root, result_type)
    end if
  end subroutine eval_tree_convert_result

  module function eval_tree_get_result_type (expr) result (type)
    integer :: type
    class(eval_tree_t), intent(in) :: expr
    if (associated (expr%root)) then
       type = expr%root%result_type
    else
       type = V_NONE
    end if
  end function eval_tree_get_result_type

  module function eval_tree_result_is_known (expr) result (flag)
    logical :: flag
    class(eval_tree_t), intent(in) :: expr
    if (associated (expr%root)) then
       select case (expr%root%result_type)
       case (V_LOG, V_INT, V_REAL)
          flag = expr%root%value_is_known
       case default
          flag = .true.
       end select
    else
       flag = .false.
    end if
  end function eval_tree_result_is_known

  module function eval_tree_result_is_known_ptr (expr) result (ptr)
    logical, pointer :: ptr
    class(eval_tree_t), intent(in) :: expr
    logical, target, save :: known = .true.
    if (associated (expr%root)) then
       select case (expr%root%result_type)
       case (V_LOG, V_INT, V_REAL)
          ptr => expr%root%value_is_known
       case default
          ptr => known
       end select
    else
       ptr => null ()
    end if
  end function eval_tree_result_is_known_ptr

  module function eval_tree_get_log (expr) result (lval)
    logical :: lval
    class(eval_tree_t), intent(in) :: expr
    if (associated (expr%root))  lval = expr%root%lval
  end function eval_tree_get_log

  module function eval_tree_get_int (expr) result (ival)
    integer :: ival
    class(eval_tree_t), intent(in) :: expr
    if (associated (expr%root)) then
       select case (expr%root%result_type)
       case (V_INT);  ival = expr%root%ival
       case (V_REAL); ival = expr%root%rval
       case (V_CMPLX); ival = expr%root%cval
       end select
    end if
  end function eval_tree_get_int

  module function eval_tree_get_real (expr) result (rval)
    real(default) :: rval
    class(eval_tree_t), intent(in) :: expr
    if (associated (expr%root)) then
       select case (expr%root%result_type)
       case (V_REAL); rval = expr%root%rval
       case (V_INT);  rval = expr%root%ival
       case (V_CMPLX);  rval = expr%root%cval
       end select
    end if
  end function eval_tree_get_real

  module function eval_tree_get_cmplx (expr) result (cval)
    complex(default) :: cval
    class(eval_tree_t), intent(in) :: expr
    if (associated (expr%root)) then
       select case (expr%root%result_type)
       case (V_CMPLX); cval = expr%root%cval
       case (V_REAL); cval = expr%root%rval
       case (V_INT);  cval = expr%root%ival
       end select
    end if
  end function eval_tree_get_cmplx

  module function eval_tree_get_pdg_array (expr) result (aval)
    type(pdg_array_t) :: aval
    class(eval_tree_t), intent(in) :: expr
    if (associated (expr%root)) then
       aval = expr%root%aval
    end if
  end function eval_tree_get_pdg_array

  module function eval_tree_get_subevt (expr) result (pval)
    type(subevt_t) :: pval
    class(eval_tree_t), intent(in) :: expr
    if (associated (expr%root)) then
       pval = expr%root%pval
    end if
  end function eval_tree_get_subevt

  module function eval_tree_get_string (expr) result (sval)
    type(string_t) :: sval
    class(eval_tree_t), intent(in) :: expr
    if (associated (expr%root)) then
       sval = expr%root%sval
    end if
  end function eval_tree_get_string

  function eval_tree_get_log_ptr (eval_tree) result (lval)
    logical, pointer :: lval
    type(eval_tree_t), intent(in) :: eval_tree
    if (associated (eval_tree%root)) then
       lval => eval_tree%root%lval
    else
       lval => null ()
    end if
  end function eval_tree_get_log_ptr

  function eval_tree_get_int_ptr (eval_tree) result (ival)
    integer, pointer :: ival
    type(eval_tree_t), intent(in) :: eval_tree
    if (associated (eval_tree%root)) then
       ival => eval_tree%root%ival
    else
       ival => null ()
    end if
  end function eval_tree_get_int_ptr

  function eval_tree_get_real_ptr (eval_tree) result (rval)
    real(default), pointer :: rval
    type(eval_tree_t), intent(in) :: eval_tree
    if (associated (eval_tree%root)) then
       rval => eval_tree%root%rval
    else
       rval => null ()
    end if
  end function eval_tree_get_real_ptr

  function eval_tree_get_cmplx_ptr (eval_tree) result (cval)
    complex(default), pointer :: cval
    type(eval_tree_t), intent(in) :: eval_tree
    if (associated (eval_tree%root)) then
       cval => eval_tree%root%cval
    else
       cval => null ()
    end if
  end function eval_tree_get_cmplx_ptr

  function eval_tree_get_subevt_ptr (eval_tree) result (pval)
    type(subevt_t), pointer :: pval
    type(eval_tree_t), intent(in) :: eval_tree
    if (associated (eval_tree%root)) then
       pval => eval_tree%root%pval
    else
       pval => null ()
    end if
  end function eval_tree_get_subevt_ptr

  function eval_tree_get_pdg_array_ptr (eval_tree) result (aval)
    type(pdg_array_t), pointer :: aval
    type(eval_tree_t), intent(in) :: eval_tree
    if (associated (eval_tree%root)) then
       aval => eval_tree%root%aval
    else
       aval => null ()
    end if
  end function eval_tree_get_pdg_array_ptr

  function eval_tree_get_string_ptr (eval_tree) result (sval)
    type(string_t), pointer :: sval
    type(eval_tree_t), intent(in) :: eval_tree
    if (associated (eval_tree%root)) then
       sval => eval_tree%root%sval
    else
       sval => null ()
    end if
  end function eval_tree_get_string_ptr

  module subroutine eval_tree_write (expr, unit, write_vars)
    class(eval_tree_t), intent(in) :: expr
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: write_vars
    integer :: u
    logical :: vl
    u = given_output_unit (unit);  if (u < 0)  return
    vl = .false.;  if (present (write_vars))  vl = write_vars
    write (u, "(1x,A)") "Evaluation tree:"
    if (associated (expr%root)) then
       call eval_node_write_rec (expr%root, unit)
    else
       write (u, "(3x,A)") "[empty]"
    end if
    if (vl)  call expr%var_list%write (unit)
  end subroutine eval_tree_write

  function eval_tree_get_md5sum (eval_tree) result (md5sum_et)
    character(32) :: md5sum_et
    type(eval_tree_t), intent(in) :: eval_tree
    integer :: u
    u = free_unit ()
    open (unit = u, status = "scratch", action = "readwrite")
    call eval_tree_write (eval_tree, unit=u)
    rewind (u)
    md5sum_et = md5sum (u)
    close (u)
  end function eval_tree_get_md5sum

  module function eval_log &
       (parse_node, var_list, subevt, is_known) result (lval)
    logical :: lval
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    logical, intent(out), optional :: is_known
    type(eval_tree_t), target :: eval_tree
    call eval_tree_init_lexpr &
         (eval_tree, parse_node, var_list, subevt)
    call eval_tree_evaluate (eval_tree)
    if (eval_tree_result_is_known (eval_tree)) then
       if (present (is_known))  is_known = .true.
       lval = eval_tree_get_log (eval_tree)
    else if (present (is_known)) then
       is_known = .false.
    else
       call eval_tree_unknown (eval_tree, parse_node)
       lval = .false.
    end if
    call eval_tree_final (eval_tree)
  end function eval_log

  module function eval_int &
       (parse_node, var_list, subevt, is_known) result (ival)
    integer :: ival
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    logical, intent(out), optional :: is_known
    type(eval_tree_t), target :: eval_tree
    call eval_tree_init_expr &
         (eval_tree, parse_node, var_list, subevt)
    call eval_tree_evaluate (eval_tree)
    if (eval_tree_result_is_known (eval_tree)) then
       if (present (is_known))  is_known = .true.
       ival = eval_tree_get_int (eval_tree)
    else if (present (is_known)) then
       is_known = .false.
    else
       call eval_tree_unknown (eval_tree, parse_node)
       ival = 0
    end if
    call eval_tree_final (eval_tree)
  end function eval_int

  module function eval_real &
       (parse_node, var_list, subevt, is_known) result (rval)
    real(default) :: rval
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    logical, intent(out), optional :: is_known
    type(eval_tree_t), target :: eval_tree
    call eval_tree_init_expr &
         (eval_tree, parse_node, var_list, subevt)
    call eval_tree_evaluate (eval_tree)
    if (eval_tree_result_is_known (eval_tree)) then
       if (present (is_known))  is_known = .true.
       rval = eval_tree_get_real (eval_tree)
    else if (present (is_known)) then
       is_known = .false.
    else
       call eval_tree_unknown (eval_tree, parse_node)
       rval = 0
    end if
    call eval_tree_final (eval_tree)
  end function eval_real

  module function eval_cmplx &
       (parse_node, var_list, subevt, is_known) result (cval)
    complex(default) :: cval
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    logical, intent(out), optional :: is_known
    type(eval_tree_t), target :: eval_tree
    call eval_tree_init_expr &
         (eval_tree, parse_node, var_list, subevt)
    call eval_tree_evaluate (eval_tree)
    if (eval_tree_result_is_known (eval_tree)) then
       if (present (is_known))  is_known = .true.
       cval = eval_tree_get_cmplx (eval_tree)
    else if (present (is_known)) then
       is_known = .false.
    else
       call eval_tree_unknown (eval_tree, parse_node)
       cval = 0
    end if
    call eval_tree_final (eval_tree)
  end function eval_cmplx

  module function eval_subevt &
       (parse_node, var_list, subevt, is_known) result (pval)
    type(subevt_t) :: pval
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    logical, intent(out), optional :: is_known
    type(eval_tree_t), target :: eval_tree
    call eval_tree_init_pexpr &
         (eval_tree, parse_node, var_list, subevt)
    call eval_tree_evaluate (eval_tree)
    if (eval_tree_result_is_known (eval_tree)) then
       if (present (is_known))  is_known = .true.
       pval = eval_tree_get_subevt (eval_tree)
    else if (present (is_known)) then
       is_known = .false.
    else
       call eval_tree_unknown (eval_tree, parse_node)
    end if
    call eval_tree_final (eval_tree)
  end function eval_subevt

  module function eval_pdg_array &
       (parse_node, var_list, subevt, is_known) result (aval)
    type(pdg_array_t) :: aval
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    logical, intent(out), optional :: is_known
    type(eval_tree_t), target :: eval_tree
    call eval_tree_init_cexpr &
         (eval_tree, parse_node, var_list, subevt)
    call eval_tree_evaluate (eval_tree)
    if (eval_tree_result_is_known (eval_tree)) then
       if (present (is_known))  is_known = .true.
       aval = eval_tree_get_pdg_array (eval_tree)
    else if (present (is_known)) then
       is_known = .false.
    else
       call eval_tree_unknown (eval_tree, parse_node)
    end if
    call eval_tree_final (eval_tree)
  end function eval_pdg_array

  module function eval_string &
       (parse_node, var_list, subevt, is_known) result (sval)
    type(string_t) :: sval
    type(parse_node_t), intent(in), target :: parse_node
    type(var_list_t), intent(in), target :: var_list
    type(subevt_t), intent(in), optional, target :: subevt
    logical, intent(out), optional :: is_known
    type(eval_tree_t), target :: eval_tree
    call eval_tree_init_sexpr &
         (eval_tree, parse_node, var_list, subevt)
    call eval_tree_evaluate (eval_tree)
    if (eval_tree_result_is_known (eval_tree)) then
       if (present (is_known))  is_known = .true.
       sval = eval_tree_get_string (eval_tree)
    else if (present (is_known)) then
       is_known = .false.
    else
       call eval_tree_unknown (eval_tree, parse_node)
       sval = ""
    end if
    call eval_tree_final (eval_tree)
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
    type(eval_tree_t), target :: eval_tree
    call eval_tree_init_expr &
         (eval_tree, parse_node, var_list, subevt)
    call eval_tree_evaluate (eval_tree)
    if (eval_tree_result_is_known (eval_tree)) then
       if (present (ival))  ival = eval_tree_get_int (eval_tree)
       if (present (rval))  rval = eval_tree_get_real (eval_tree)
       if (present (cval))  cval = eval_tree_get_cmplx (eval_tree)
       if (present (is_known))  is_known = .true.
    else
       call eval_tree_unknown (eval_tree, parse_node)
       if (present (ival))  ival = 0
       if (present (rval))  rval = 0
       if (present (cval))  cval = 0
       if (present (is_known))  is_known = .false.
    end if
    if (present (result_type))  &
         result_type = eval_tree_get_result_type (eval_tree)
    call eval_tree_final (eval_tree)
  end subroutine eval_numeric

  subroutine eval_tree_unknown (eval_tree, parse_node)
    type(eval_tree_t), intent(in) :: eval_tree
    type(parse_node_t), intent(in) :: parse_node
    call parse_node_write_rec (parse_node)
    call eval_tree_write (eval_tree)
    call msg_error ("Evaluation yields an undefined result, inserting default")
  end subroutine eval_tree_unknown

  module subroutine eval_tree_factory_write (expr_factory, unit)
    class(eval_tree_factory_t), intent(in) :: expr_factory
    integer, intent(in), optional :: unit
    if (associated (expr_factory%pn)) then
       call parse_node_write_rec (expr_factory%pn, unit)
    end if
  end subroutine eval_tree_factory_write

  module subroutine eval_tree_factory_init (expr_factory, pn)
    class(eval_tree_factory_t), intent(out) :: expr_factory
    type(parse_node_t), intent(in), pointer :: pn
    expr_factory%pn => pn
  end subroutine eval_tree_factory_init


end submodule eval_trees_s

