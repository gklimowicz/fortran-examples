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

submodule (slha_interface) slha_interface_s

  use io_units
  use constants
  use string_utils, only: upper_case
  use system_defs, only: VERSION_STRING
  use system_defs, only: EOF
  use diagnostics
  use ifiles
  use syntax_rules

  implicit none

  integer, parameter :: MODE_SKIP = 0, MODE_DATA = 1, MODE_INFO = 2

  integer, parameter :: MDL_MSSM = 0
  integer, parameter :: MDL_NMSSM = 1
  integer, parameter :: MSSM_GENERIC = 0
  integer, parameter :: MSSM_SUGRA = 1
  integer, parameter :: MSSM_GMSB = 2
  integer, parameter :: MSSM_AMSB = 3


  type(syntax_t), target :: syntax_slha


contains

  subroutine slha_preprocess (stream, custom_block_name, ifile)
    type(stream_t), intent(inout), target :: stream
    type(string_t), dimension(:), intent(in) :: custom_block_name
    type(ifile_t), intent(out) :: ifile
    type(string_t) :: buffer, line, item
    integer :: iostat
    integer :: mode
    mode = MODE
    SCAN_FILE: do
       call stream_get_record (stream, buffer, iostat)
       select case (iostat)
       case (0)
          call split (buffer, line, "#")
          if (len_trim (line) == 0)  cycle SCAN_FILE
          select case (char (extract (line, 1, 1)))
          case ("B", "b")
             call check_block_handling (line, custom_block_name, mode)
             call ifile_append (ifile, line // "$")
          case ("D", "d")
             mode = MODE_DATA
             call ifile_append (ifile, line // "$")
          case (" ")
             select case (mode)
             case (MODE_DATA)
                call ifile_append (ifile, "DATA" // line // "$")
             case (MODE_INFO)
                line = adjustl (line)
                call split (line, item, " ")
                call ifile_append (ifile, "INFO" // " " // item // " " &
                     // '"' // trim (adjustl (line)) // '" $')
             end select
          case default
             call msg_message (char (line))
             call msg_fatal ("SLHA: Incomprehensible line")
          end select
       case (EOF)
          exit SCAN_FILE
       case default
          call msg_fatal ("SLHA: I/O error occured while reading SLHA input")
       end select
    end do SCAN_FILE
  end subroutine slha_preprocess

  subroutine check_block_handling (line, custom_block_name, mode)
    type(string_t), intent(inout) :: line
    type(string_t), dimension(:), intent(in) :: custom_block_name
    integer, intent(out) :: mode
    type(string_t) :: buffer, key, block_name
    integer :: i
    buffer = trim (line)
    call split (buffer, key, " ")
    buffer = adjustl (buffer)
    call split (buffer, block_name, " ")
    buffer = adjustl (buffer)
    block_name = trim (adjustl (upper_case (block_name)))
    select case (char (block_name))
    case ("MODSEL", "MINPAR", "SMINPUTS")
       mode = MODE_DATA
    case ("MASS")
       mode = MODE_DATA
    case ("NMIX", "UMIX", "VMIX", "STOPMIX", "SBOTMIX", "STAUMIX")
       mode = MODE_DATA
    case ("NMHMIX", "NMAMIX", "NMNMIX", "NMSSMRUN")
       mode = MODE_DATA
    case ("ALPHA", "HMIX")
       mode = MODE_DATA
    case ("AU", "AD", "AE")
       mode = MODE_DATA
    case ("SPINFO", "DCINFO")
       mode = MODE_INFO
    case default
       mode = MODE_SKIP
       CHECK_CUSTOM_NAMES: do i = 1, size (custom_block_name)
          if (block_name == custom_block_name(i)) then
             mode = MODE_DATA
             call mangle_keywords (block_name)
             line = key // " " // block_name // " " // trim (buffer)
             exit CHECK_CUSTOM_NAMES
          end if
       end do CHECK_CUSTOM_NAMES
    end select
  end subroutine check_block_handling

  subroutine mangle_keywords (name)
    type(string_t), intent(inout) :: name
    select case (char (name))
    case ("BLOCK", "DATA", "INFO", "DECAY")
       name = name // "_"
    end select
  end subroutine mangle_keywords
  
  subroutine demangle_keywords (name)
    type(string_t), intent(inout) :: name
    select case (char (name))
    case ("BLOCK_", "DATA_", "INFO_", "DECAY_")
       name = extract (name, 1, len(name)-1)
    end select
  end subroutine demangle_keywords
  
  module subroutine syntax_slha_init ()
    type(ifile_t) :: ifile
    call define_slha_syntax (ifile)
    call syntax_init (syntax_slha, ifile)
    call ifile_final (ifile)
  end subroutine syntax_slha_init

  module subroutine syntax_slha_final ()
    call syntax_final (syntax_slha)
  end subroutine syntax_slha_final

  module subroutine syntax_slha_write (unit)
    integer, intent(in), optional :: unit
    call syntax_write (syntax_slha, unit)
  end subroutine syntax_slha_write

  subroutine define_slha_syntax (ifile)
    type(ifile_t), intent(inout) :: ifile
    call ifile_append (ifile, "SEQ slha = chunk*")
    call ifile_append (ifile, "ALT chunk = block_def | decay_def")
    call ifile_append (ifile, "SEQ block_def = " &
         // "BLOCK blockgen '$' block_line*")
    call ifile_append (ifile, "ALT blockgen = block_spec | q_spec")
    call ifile_append (ifile, "KEY BLOCK")
    call ifile_append (ifile, "SEQ q_spec = QNUMBERS pdg_code")
    call ifile_append (ifile, "KEY QNUMBERS")
    call ifile_append (ifile, "SEQ block_spec = block_name qvalue?")
    call ifile_append (ifile, "IDE block_name")
    call ifile_append (ifile, "SEQ qvalue = qname '=' qitem")
    call ifile_append (ifile, "IDE qname")
    call ifile_append (ifile, "ALT qitem = signed_real | real")
    call ifile_append (ifile, "SEQ signed_real = sign real")
    call ifile_append (ifile, "KEY '='")
    call ifile_append (ifile, "REA real")
    call ifile_append (ifile, "KEY '$'")
    call ifile_append (ifile, "ALT block_line = block_data | block_info")
    call ifile_append (ifile, "SEQ block_data = DATA data_line '$'")
    call ifile_append (ifile, "KEY DATA")
    call ifile_append (ifile, "SEQ data_line = data_item+")
    call ifile_append (ifile, "ALT data_item = signed_number | number")
    call ifile_append (ifile, "SEQ signed_number = sign number")
    call ifile_append (ifile, "ALT sign = '+' | '-'")
    call ifile_append (ifile, "ALT number = integer | real")
    call ifile_append (ifile, "INT integer")
    call ifile_append (ifile, "KEY '-'")
    call ifile_append (ifile, "KEY '+'")
    call ifile_append (ifile, "SEQ block_info = INFO info_line '$'")
    call ifile_append (ifile, "KEY INFO")
    call ifile_append (ifile, "SEQ info_line = integer string_literal")
    call ifile_append (ifile, "QUO string_literal = '""'...'""'")
    call ifile_append (ifile, "SEQ decay_def = " &
         // "DECAY decay_spec '$' decay_data*")
    call ifile_append (ifile, "KEY DECAY")
    call ifile_append (ifile, "SEQ decay_spec = pdg_code data_item")
    call ifile_append (ifile, "ALT pdg_code = signed_integer | integer")
    call ifile_append (ifile, "SEQ signed_integer = sign integer")
    call ifile_append (ifile, "SEQ decay_data = DATA decay_line '$'")
    call ifile_append (ifile, "SEQ decay_line = data_item integer pdg_code+")
  end subroutine define_slha_syntax

  module subroutine lexer_init_slha (lexer)
    type(lexer_t), intent(out) :: lexer
    call lexer_init (lexer, &
         comment_chars = "#", &
         quote_chars = '"', &
         quote_match = '"', &
         single_chars = "+-=$", &
         special_class = [ "" ], &
         keyword_list = syntax_get_keyword_list_ptr (syntax_slha), &
         upper_case_keywords = .true.)  ! $
  end subroutine lexer_init_slha

  function slha_get_block_ptr &
       (parse_tree, block_name, required) result (pn_block)
    type(parse_node_t), pointer :: pn_block
    type(parse_tree_t), intent(in) :: parse_tree
    type(string_t), intent(in) :: block_name
    type(string_t) :: block_def
    logical, intent(in) :: required
    type(parse_node_t), pointer :: pn_root, pn_block_spec, pn_block_name
    pn_root => parse_tree%get_root_ptr ()
    pn_block => parse_node_get_sub_ptr (pn_root)
    do while (associated (pn_block))
       select case (char (parse_node_get_rule_key (pn_block)))
       case ("block_def")
          pn_block_spec => parse_node_get_sub_ptr (pn_block, 2)
          pn_block_name => parse_node_get_sub_ptr (pn_block_spec)
          select case (char (pn_block_name%get_rule_key ()))
          case ("block_name")
             block_def = trim (adjustl (upper_case &
                  (pn_block_name%get_string ())))
          case ("QNUMBERS")
             block_def = "QNUMBERS"
          end select
          if (block_def == block_name) then
             return
          end if
       end select
       pn_block => parse_node_get_next_ptr (pn_block)
    end do
    if (required) then
       call msg_fatal ("SLHA: block '" // char (block_name) // "' not found")
    end if
  end function slha_get_block_ptr

  function slha_get_first_decay_ptr (parse_tree) result (pn_decay)
    type(parse_node_t), pointer :: pn_decay
    type(parse_tree_t), intent(in) :: parse_tree
    type(parse_node_t), pointer :: pn_root
    pn_root => parse_tree%get_root_ptr ()
    pn_decay => parse_node_get_sub_ptr (pn_root)
    do while (associated (pn_decay))
       select case (char (parse_node_get_rule_key (pn_decay)))
       case ("decay_def")
          return
       end select
       pn_decay => parse_node_get_next_ptr (pn_decay)
    end do
  end function slha_get_first_decay_ptr

  function slha_get_next_decay_ptr (pn_block) result (pn_decay)
    type(parse_node_t), pointer :: pn_decay
    type(parse_node_t), intent(in), target :: pn_block
    pn_decay => parse_node_get_next_ptr (pn_block)
    do while (associated (pn_decay))
       select case (char (parse_node_get_rule_key (pn_decay)))
       case ("decay_def")
          return
       end select
       pn_decay => parse_node_get_next_ptr (pn_decay)
    end do
  end function slha_get_next_decay_ptr

  subroutine slha_find_index_ptr (pn_block, pn_data, pn_item, code)
    type(parse_node_t), intent(in), target :: pn_block
    type(parse_node_t), intent(out), pointer :: pn_data
    type(parse_node_t), intent(out), pointer :: pn_item
    integer, intent(in) :: code
    pn_data => parse_node_get_sub_ptr (pn_block, 4)
    call slha_next_index_ptr (pn_data, pn_item, code)
  end subroutine slha_find_index_ptr

  subroutine slha_find_index_pair_ptr (pn_block, pn_data, pn_item, code1, code2)
    type(parse_node_t), intent(in), target :: pn_block
    type(parse_node_t), intent(out), pointer :: pn_data
    type(parse_node_t), intent(out), pointer :: pn_item
    integer, intent(in) :: code1, code2
    pn_data => parse_node_get_sub_ptr (pn_block, 4)
    call slha_next_index_pair_ptr (pn_data, pn_item, code1, code2)
  end subroutine slha_find_index_pair_ptr

  subroutine slha_next_index_ptr (pn_data, pn_item, code)
    type(parse_node_t), intent(inout), pointer :: pn_data
    integer, intent(in) :: code
    type(parse_node_t), intent(out), pointer :: pn_item
    type(parse_node_t), pointer :: pn_line, pn_code
    do while (associated (pn_data))
       pn_line => parse_node_get_sub_ptr (pn_data, 2)
       pn_code => parse_node_get_sub_ptr (pn_line)
       select case (char (parse_node_get_rule_key (pn_code)))
       case ("integer")
          if (parse_node_get_integer (pn_code) == code) then
             pn_item => parse_node_get_next_ptr (pn_code)
             return
          end if
       end select
       pn_data => parse_node_get_next_ptr (pn_data)
    end do
    pn_item => null ()
  end subroutine slha_next_index_ptr

  subroutine slha_next_index_pair_ptr (pn_data, pn_item, code1, code2)
    type(parse_node_t), intent(inout), pointer :: pn_data
    integer, intent(in) :: code1, code2
    type(parse_node_t), intent(out), pointer :: pn_item
    type(parse_node_t), pointer :: pn_line, pn_code1, pn_code2
    do while (associated (pn_data))
       pn_line => parse_node_get_sub_ptr (pn_data, 2)
       pn_code1 => parse_node_get_sub_ptr (pn_line)
       select case (char (parse_node_get_rule_key (pn_code1)))
       case ("integer")
          if (parse_node_get_integer (pn_code1) == code1) then
             pn_code2 => parse_node_get_next_ptr (pn_code1)
             if (associated (pn_code2)) then
                select case (char (parse_node_get_rule_key (pn_code2)))
                case ("integer")
                   if (parse_node_get_integer (pn_code2) == code2) then
                      pn_item => parse_node_get_next_ptr (pn_code2)
                      return
                   end if
                end select
             end if
          end if
       end select
       pn_data => parse_node_get_next_ptr (pn_data)
    end do
    pn_item => null ()
  end subroutine slha_next_index_pair_ptr

  subroutine retrieve_strings_in_block (pn_block, code, str_array)
    type(parse_node_t), intent(in), target :: pn_block
    integer, intent(in) :: code
    type(string_t), dimension(:), allocatable, intent(out) :: str_array
    type(parse_node_t), pointer :: pn_data, pn_item
    type(str_entry_t), pointer :: first => null ()
    type(str_entry_t), pointer :: current => null ()
    integer :: n
    n = 0
    call slha_find_index_ptr (pn_block, pn_data, pn_item, code)
    if (associated (pn_item)) then
       n = n + 1
       allocate (first)
       first%str = parse_node_get_string (pn_item)
       current => first
       do while (associated (pn_data))
          pn_data => parse_node_get_next_ptr (pn_data)
          call slha_next_index_ptr (pn_data, pn_item, code)
          if (associated (pn_item)) then
             n = n + 1
             allocate (current%next)
             current => current%next
             current%str = parse_node_get_string (pn_item)
          end if
       end do
       allocate (str_array (n))
       n = 0
       do while (associated (first))
          n = n + 1
          current => first
          str_array(n) = current%str
          first => first%next
          deallocate (current)
       end do
    else
       allocate (str_array (0))
    end if
  end subroutine retrieve_strings_in_block

  function get_parameter_in_block (pn_block, code, name, var_list) result (var)
    real(default) :: var
    type(parse_node_t), intent(in), target :: pn_block
    integer, intent(in) :: code
    type(string_t), intent(in) :: name
    type(var_list_t), intent(in), target :: var_list
    type(parse_node_t), pointer :: pn_data, pn_item
    call slha_find_index_ptr (pn_block, pn_data, pn_item, code)
    if (associated (pn_item)) then
       var = get_real_parameter (pn_item)
    else
       var = var_list%get_rval (name)
    end if
  end function get_parameter_in_block

  subroutine set_data_item (pn_block, code, name, var_list)
    type(parse_node_t), intent(in), target :: pn_block
    integer, intent(in) :: code
    type(string_t), intent(in) :: name
    type(var_list_t), intent(inout), target :: var_list
    type(parse_node_t), pointer :: pn_data, pn_item
    call slha_find_index_ptr (pn_block, pn_data, pn_item, code)
    if (associated (pn_item)) then
       call var_list%set_real (name,  get_real_parameter (pn_item), &
             is_known=.true., ignore=.true.)
    end if
  end subroutine set_data_item

  subroutine set_matrix_element (pn_block, code1, code2, name, var_list)
    type(parse_node_t), intent(in), target :: pn_block
    integer, intent(in) :: code1, code2
    type(string_t), intent(in) :: name
    type(var_list_t), intent(inout), target :: var_list
    type(parse_node_t), pointer :: pn_data, pn_item
    call slha_find_index_pair_ptr (pn_block, pn_data, pn_item, code1, code2)
    if (associated (pn_item)) then
       call var_list%set_real (name, get_real_parameter (pn_item), &
             is_known=.true., ignore=.true.)
    end if
  end subroutine set_matrix_element

  subroutine write_integer_data_item (u, code, name, var_list, comment)
    integer, intent(in) :: u
    integer, intent(in) :: code
    type(string_t), intent(in) :: name
    type(var_list_t), intent(in) :: var_list
    character(*), intent(in) :: comment
    integer :: item
    if (var_list%contains (name)) then
       item = nint (var_list%get_rval (name))
       call write_integer_parameter (u, code, item, comment)
    end if
  end subroutine write_integer_data_item

  subroutine write_real_data_item (u, code, name, var_list, comment)
    integer, intent(in) :: u
    integer, intent(in) :: code
    type(string_t), intent(in) :: name
    type(var_list_t), intent(in) :: var_list
    character(*), intent(in) :: comment
    real(default) :: item
    if (var_list%contains (name)) then
       item = var_list%get_rval (name)
       call write_real_parameter (u, code, item, comment)
    end if
  end subroutine write_real_data_item

  subroutine write_matrix_element (u, code1, code2, name, var_list, comment)
    integer, intent(in) :: u
    integer, intent(in) :: code1, code2
    type(string_t), intent(in) :: name
    type(var_list_t), intent(in) :: var_list
    character(*), intent(in) :: comment
    real(default) :: item
    if (var_list%contains (name)) then
       item = var_list%get_rval (name)
       call write_real_matrix_element (u, code1, code2, item, comment)
    end if
  end subroutine write_matrix_element

  subroutine write_block_header (u, name, comment)
    integer, intent(in) :: u
    character(*), intent(in) :: name, comment
    write (u, "(A,1x,A,3x,'#',1x,A)")  "BLOCK", name, comment
  end subroutine write_block_header

  function get_real_parameter (pn_item) result (var)
    real(default) :: var
    type(parse_node_t), intent(in), target :: pn_item
    type(parse_node_t), pointer :: pn_sign, pn_var
    integer :: sign
    select case (char (parse_node_get_rule_key (pn_item)))
    case ("signed_number")
       pn_sign => parse_node_get_sub_ptr (pn_item)
       pn_var  => parse_node_get_next_ptr (pn_sign)
       select case (char (parse_node_get_key (pn_sign)))
       case ("+");  sign = +1
       case ("-");  sign = -1
       end select
    case default
       sign = +1
       pn_var => pn_item
    end select
    select case (char (parse_node_get_rule_key (pn_var)))
    case ("integer");  var = sign * parse_node_get_integer (pn_var)
    case ("real");     var = sign * parse_node_get_real (pn_var)
    end select
  end function get_real_parameter

  function get_integer_parameter (pn_item) result (var)
    integer :: var
    type(parse_node_t), intent(in), target :: pn_item
    type(parse_node_t), pointer :: pn_sign, pn_var
    integer :: sign
    select case (char (parse_node_get_rule_key (pn_item)))
    case ("signed_integer")
       pn_sign => parse_node_get_sub_ptr (pn_item)
       pn_var  => parse_node_get_next_ptr (pn_sign)
       select case (char (parse_node_get_key (pn_sign)))
       case ("+");  sign = +1
       case ("-");  sign = -1
       end select
    case ("integer")
       sign = +1
       pn_var => pn_item
    case default
       call parse_node_write (pn_var)
       call msg_error ("SLHA: Integer parameter expected")
       var = 0
       return
    end select
    var = sign * parse_node_get_integer (pn_var)
  end function get_integer_parameter

  subroutine write_integer_parameter (u, code, item, comment)
    integer, intent(in) :: u
    integer, intent(in) :: code
    integer, intent(in) :: item
    character(*), intent(in) :: comment
1   format (1x, I9, 3x, 3x, I9, 4x, 3x, '#', 1x, A)
    write (u, 1)  code, item, comment
  end subroutine write_integer_parameter

  subroutine write_real_parameter (u, code, item, comment)
    integer, intent(in) :: u
    integer, intent(in) :: code
    real(default), intent(in) :: item
    character(*), intent(in) :: comment
1   format (1x, I9, 3x, 1P, E16.8, 0P, 3x, '#', 1x, A)
    write (u, 1)  code, item, comment
  end subroutine write_real_parameter

  subroutine write_real_matrix_element (u, code1, code2, item, comment)
    integer, intent(in) :: u
    integer, intent(in) :: code1, code2
    real(default), intent(in) :: item
    character(*), intent(in) :: comment
1   format (1x, I2, 1x, I2, 3x, 1P, E16.8, 0P, 3x, '#', 1x, A)
    write (u, 1)  code1, code2, item, comment
  end subroutine write_real_matrix_element

  module subroutine slha_interpret_parse_tree &
       (parse_tree, model, input, spectrum, decays)
    type(parse_tree_t), intent(in) :: parse_tree
    type(model_t), intent(inout), target :: model
    logical, intent(in) :: input, spectrum, decays
    logical :: errors
    integer :: mssm_type
    if (model%supports_custom_slha ()) then
       call slha_handle_custom_file (parse_tree, model)
    else
       call slha_handle_MODSEL (parse_tree, model, mssm_type)
       if (input) then
          call slha_handle_SMINPUTS (parse_tree, model)
          call slha_handle_MINPAR (parse_tree, model, mssm_type)
       end if
       if (spectrum) then
          call slha_handle_info_block (parse_tree, "SPINFO", errors)
          if (errors)  return
          call slha_handle_MASS (parse_tree, model)
          call slha_handle_matrix_block (parse_tree, "NMIX", "mn_", 4, 4, model)
          call slha_handle_matrix_block (parse_tree, "NMNMIX", "mixn_", 5, 5, model)
          call slha_handle_matrix_block (parse_tree, "UMIX", "mu_", 2, 2, model)
          call slha_handle_matrix_block (parse_tree, "VMIX", "mv_", 2, 2, model)
          call slha_handle_matrix_block (parse_tree, "STOPMIX", "mt_", 2, 2, model)
          call slha_handle_matrix_block (parse_tree, "SBOTMIX", "mb_", 2, 2, model)
          call slha_handle_matrix_block (parse_tree, "STAUMIX", "ml_", 2, 2, model)
          call slha_handle_matrix_block (parse_tree, "NMHMIX", "mixh0_", 3, 3, model)
          call slha_handle_matrix_block (parse_tree, "NMAMIX", "mixa0_", 2, 3, model)
          call slha_handle_ALPHA (parse_tree, model)
          call slha_handle_HMIX (parse_tree, model)
          call slha_handle_NMSSMRUN (parse_tree, model)
          call slha_handle_matrix_block (parse_tree, "AU", "Au_", 3, 3, model)
          call slha_handle_matrix_block (parse_tree, "AD", "Ad_", 3, 3, model)
          call slha_handle_matrix_block (parse_tree, "AE", "Ae_", 3, 3, model)
       end if
    end if
    if (decays) then
       call slha_handle_info_block (parse_tree, "DCINFO", errors)
       if (errors)  return
       call slha_handle_decays (parse_tree, model)
    end if
  end subroutine slha_interpret_parse_tree

  subroutine slha_handle_info_block (parse_tree, block_name, errors)
    type(parse_tree_t), intent(in) :: parse_tree
    character(*), intent(in) :: block_name
    logical, intent(out) :: errors
    type(parse_node_t), pointer :: pn_block
    type(string_t), dimension(:), allocatable :: msg
    integer :: i
    pn_block => slha_get_block_ptr &
         (parse_tree, var_str (block_name), required=.true.)
    if (.not. associated (pn_block)) then
       call msg_error ("SLHA: Missing info block '" &
            // trim (block_name) // "'; ignored.")
       errors = .true.
       return
    end if
    select case (block_name)
    case ("SPINFO")
       call msg_message ("SLHA: SUSY spectrum program info:")
    case ("DCINFO")
       call msg_message ("SLHA: SUSY decay program info:")
    end select
    call retrieve_strings_in_block (pn_block, 1, msg)
    do i = 1, size (msg)
       call msg_message ("SLHA: " // char (msg(i)))
    end do
    call retrieve_strings_in_block (pn_block, 2, msg)
    do i = 1, size (msg)
       call msg_message ("SLHA: " // char (msg(i)))
    end do
    call retrieve_strings_in_block (pn_block, 3, msg)
    do i = 1, size (msg)
       call msg_warning ("SLHA: " // char (msg(i)))
    end do
    call retrieve_strings_in_block (pn_block, 4, msg)
    do i = 1, size (msg)
       call msg_error ("SLHA: " // char (msg(i)))
    end do
    errors = size (msg) > 0
  end subroutine slha_handle_info_block

  subroutine slha_handle_MODSEL (parse_tree, model, mssm_type)
    type(parse_tree_t), intent(in) :: parse_tree
    type(model_t), intent(in), target :: model
    integer, intent(out) :: mssm_type
    type(parse_node_t), pointer :: pn_block, pn_data, pn_item
    type(string_t) :: model_name
    pn_block => slha_get_block_ptr &
         (parse_tree, var_str ("MODSEL"), required=.false.)
    if (.not. associated (pn_block)) then
       mssm_type = -1
       return
    end if
    call slha_find_index_ptr (pn_block, pn_data, pn_item, 1)
    if (associated (pn_item)) then
       mssm_type = get_integer_parameter (pn_item)
    else
       mssm_type = MSSM_GENERIC
    end if
    call slha_find_index_ptr (pn_block, pn_data, pn_item, 3)
    if (associated (pn_item)) then
       select case (parse_node_get_integer (pn_item))
       case (MDL_MSSM);  model_name = "MSSM"
       case (MDL_NMSSM); model_name = "NMSSM"
       case default
          call msg_fatal ("SLHA: unknown model code in MODSEL")
          return
       end select
    else
       model_name = "MSSM"
    end if
    call slha_find_index_ptr (pn_block, pn_data, pn_item, 4)
    if (associated (pn_item)) then
      call msg_fatal (" R-parity violation is currently not supported by WHIZARD.")
    end if
    call slha_find_index_ptr (pn_block, pn_data, pn_item, 5)
    if (associated (pn_item)) then
      call msg_fatal (" CP violation is currently not supported by WHIZARD.")
    end if
    select case (char (model_name))
    case ("MSSM")
       select case (char (model%get_name ()))
       case ("MSSM","MSSM_CKM","MSSM_Grav","MSSM_Hgg")
          model_name = model%get_name ()
       case default
          call msg_fatal ("Selected model '" &
               // char (model%get_name ()) // "' does not match model '" &
               // char (model_name) // "' in SLHA input file.")
          return
       end select
    case ("NMSSM")
       select case (char (model%get_name ()))
       case ("NMSSM","NMSSM_CKM","NMSSM_Hgg")
          model_name = model%get_name ()
       case default
          call msg_fatal ("Selected model '" &
               // char (model%get_name ()) // "' does not match model '" &
               // char (model_name) // "' in SLHA input file.")
          return
       end select
    case default
       call msg_bug ("SLHA model name '" &
            // char (model_name) // "' not recognized.")
       return
    end select
    call msg_message ("SLHA: Initializing model '" // char (model_name) // "'")
  end subroutine slha_handle_MODSEL

  subroutine slha_write_MODSEL (u, model, mssm_type)
    integer, intent(in) :: u
    type(model_t), intent(in), target :: model
    integer, intent(out) :: mssm_type
    type(var_list_t), pointer :: var_list
    integer :: model_id
    type(string_t) :: mtype_string
    var_list => model%get_var_list_ptr ()
    if (var_list%contains (var_str ("mtype"))) then
       mssm_type = nint (var_list%get_rval (var_str ("mtype")))
    else
       call msg_error ("SLHA: parameter 'mtype' (SUSY breaking scheme) " &
            // "is unknown in current model, no SLHA output possible")
       mssm_type = -1
       return
    end if
    call write_block_header (u, "MODSEL", "SUSY model selection")
    select case (mssm_type)
    case (0);  mtype_string = "Generic MSSM"
    case (1);  mtype_string = "SUGRA"
    case (2);  mtype_string = "GMSB"
    case (3);  mtype_string = "AMSB"
    case default
       mtype_string = "unknown"
    end select
    call write_integer_parameter (u, 1, mssm_type, &
         "SUSY-breaking scheme: " // char (mtype_string))
    select case (char (model%get_name ()))
    case ("MSSM");  model_id = MDL_MSSM
    case ("NMSSM"); model_id = MDL_NMSSM
    case default
       model_id = 0
    end select
    call write_integer_parameter (u, 3, model_id, &
         "SUSY model type: " // char (model%get_name ()))
  end subroutine slha_write_MODSEL

  subroutine slha_handle_SMINPUTS (parse_tree, model)
    type(parse_tree_t), intent(in) :: parse_tree
    type(model_t), intent(inout), target :: model
    type(parse_node_t), pointer :: pn_block
    real(default) :: alpha_em_i, GF, alphas, mZ
    real(default) :: ee, vv, cw_sw, cw2, mW
    real(default) :: mb, mtop, mtau
    type(var_list_t), pointer :: var_list
    var_list => model%get_var_list_ptr ()
    pn_block => slha_get_block_ptr &
         (parse_tree, var_str ("SMINPUTS"), required=.true.)
    if (.not. (associated (pn_block)))  return
    alpha_em_i = &
         get_parameter_in_block (pn_block, 1, var_str ("alpha_em_i"), var_list)
    GF = get_parameter_in_block (pn_block, 2, var_str ("GF"), var_list)
    alphas = &
         get_parameter_in_block (pn_block, 3, var_str ("alphas"), var_list)
    mZ   = get_parameter_in_block (pn_block, 4, var_str ("mZ"), var_list)
    mb   = get_parameter_in_block (pn_block, 5, var_str ("mb"), var_list)
    mtop = get_parameter_in_block (pn_block, 6, var_str ("mtop"), var_list)
    mtau = get_parameter_in_block (pn_block, 7, var_str ("mtau"), var_list)
    ee = sqrt (4 * pi / alpha_em_i)
    vv = 1 / sqrt (sqrt (2._default) * GF)
    cw_sw = ee * vv / (2 * mZ)
    if (2*cw_sw <= 1) then
       cw2 = (1 + sqrt (1 - 4 * cw_sw**2)) / 2
       mW = mZ * sqrt (cw2)
       call var_list%set_real (var_str ("GF"), GF, .true.)
       call var_list%set_real (var_str ("mZ"), mZ, .true.)
       call var_list%set_real (var_str ("mW"), mW, .true.)
       call var_list%set_real (var_str ("mtau"), mtau, .true.)
       call var_list%set_real (var_str ("mb"), mb, .true.)
       call var_list%set_real (var_str ("mtop"), mtop, .true.)
       call var_list%set_real (var_str ("alphas"), alphas, .true.)
    else
       call msg_fatal ("SLHA: Unphysical SM parameter values")
       return
    end if
  end subroutine slha_handle_SMINPUTS

  subroutine slha_write_SMINPUTS (u, model)
    integer, intent(in) :: u
    type(model_t), intent(in), target :: model
    type(var_list_t), pointer :: var_list
    var_list => model%get_var_list_ptr ()
    call write_block_header (u, "SMINPUTS", "SM input parameters")
    call write_real_data_item (u, 1, var_str ("alpha_em_i"), var_list, &
         "Inverse electromagnetic coupling alpha (Z pole)")
    call write_real_data_item (u, 2, var_str ("GF"), var_list, &
         "Fermi constant")
    call write_real_data_item (u, 3, var_str ("alphas"), var_list, &
         "Strong coupling alpha_s (Z pole)")
    call write_real_data_item (u, 4, var_str ("mZ"), var_list, &
         "Z mass")
    call write_real_data_item (u, 5, var_str ("mb"), var_list, &
         "b running mass (at mb)")
    call write_real_data_item (u, 6, var_str ("mtop"), var_list, &
         "top mass")
    call write_real_data_item (u, 7, var_str ("mtau"), var_list, &
         "tau mass")
  end subroutine slha_write_SMINPUTS

  subroutine slha_handle_MINPAR (parse_tree, model, mssm_type)
    type(parse_tree_t), intent(in) :: parse_tree
    type(model_t), intent(inout), target :: model
    integer, intent(in) :: mssm_type
    type(var_list_t), pointer :: var_list
    type(parse_node_t), pointer :: pn_block
    var_list => model%get_var_list_ptr ()
    call var_list%set_real &
         (var_str ("mtype"),  real(mssm_type, default), is_known=.true.)
    pn_block => slha_get_block_ptr &
         (parse_tree, var_str ("MINPAR"), required=.true.)
    select case (mssm_type)
    case (MSSM_SUGRA)
       call set_data_item (pn_block, 1, var_str ("m_zero"), var_list)
       call set_data_item (pn_block, 2, var_str ("m_half"), var_list)
       call set_data_item (pn_block, 3, var_str ("tanb"), var_list)
       call set_data_item (pn_block, 4, var_str ("sgn_mu"), var_list)
       call set_data_item (pn_block, 5, var_str ("A0"), var_list)
    case (MSSM_GMSB)
       call set_data_item (pn_block, 1, var_str ("Lambda"), var_list)
       call set_data_item (pn_block, 2, var_str ("M_mes"), var_list)
       call set_data_item (pn_block, 3, var_str ("tanb"), var_list)
       call set_data_item (pn_block, 4, var_str ("sgn_mu"), var_list)
       call set_data_item (pn_block, 5, var_str ("N_5"), var_list)
       call set_data_item (pn_block, 6, var_str ("c_grav"), var_list)
    case (MSSM_AMSB)
       call set_data_item (pn_block, 1, var_str ("m_zero"), var_list)
       call set_data_item (pn_block, 2, var_str ("m_grav"), var_list)
       call set_data_item (pn_block, 3, var_str ("tanb"), var_list)
       call set_data_item (pn_block, 4, var_str ("sgn_mu"), var_list)
    case default
       call set_data_item (pn_block, 3, var_str ("tanb"), var_list)
    end select
  end subroutine slha_handle_MINPAR

  subroutine slha_write_MINPAR (u, model, mssm_type)
    integer, intent(in) :: u
    type(model_t), intent(in), target :: model
    integer, intent(in) :: mssm_type
    type(var_list_t), pointer :: var_list
    var_list => model%get_var_list_ptr ()
    call write_block_header (u, "MINPAR", "Basic SUSY input parameters")
    select case (mssm_type)
    case (MSSM_SUGRA)
       call write_real_data_item (u, 1, var_str ("m_zero"), var_list, &
            "Common scalar mass")
       call write_real_data_item (u, 2, var_str ("m_half"), var_list, &
            "Common gaugino mass")
       call write_real_data_item (u, 3, var_str ("tanb"), var_list, &
            "tan(beta)")
       call write_integer_data_item (u, 4, &
            var_str ("sgn_mu"), var_list, &
            "Sign of mu")
       call write_real_data_item (u, 5, var_str ("A0"), var_list, &
            "Common trilinear coupling")
    case (MSSM_GMSB)
       call write_real_data_item (u, 1, var_str ("Lambda"), var_list, &
            "Soft-breaking scale")
       call write_real_data_item (u, 2, var_str ("M_mes"), var_list, &
            "Messenger scale")
       call write_real_data_item (u, 3, var_str ("tanb"), var_list, &
            "tan(beta)")
       call write_integer_data_item (u, 4, &
            var_str ("sgn_mu"), var_list, &
            "Sign of mu")
       call write_integer_data_item (u, 5, var_str ("N_5"), var_list, &
            "Messenger index")
       call write_real_data_item (u, 6, var_str ("c_grav"), var_list, &
            "Gravitino mass factor")
    case (MSSM_AMSB)
       call write_real_data_item (u, 1, var_str ("m_zero"), var_list, &
            "Common scalar mass")
       call write_real_data_item (u, 2, var_str ("m_grav"), var_list, &
            "Gravitino mass")
       call write_real_data_item (u, 3, var_str ("tanb"), var_list, &
            "tan(beta)")
       call write_integer_data_item (u, 4, &
            var_str ("sgn_mu"), var_list, &
            "Sign of mu")
    case default
       call write_real_data_item (u, 3, var_str ("tanb"), var_list, &
            "tan(beta)")
    end select
  end subroutine slha_write_MINPAR

  subroutine slha_handle_MASS (parse_tree, model)
    type(parse_tree_t), intent(in) :: parse_tree
    type(model_t), intent(inout), target :: model
    type(parse_node_t), pointer :: pn_block, pn_data, pn_line, pn_code
    type(parse_node_t), pointer :: pn_mass
    integer :: pdg
    real(default) :: mass
    pn_block => slha_get_block_ptr &
         (parse_tree, var_str ("MASS"), required=.true.)
    if (.not. (associated (pn_block)))  return
    pn_data => parse_node_get_sub_ptr (pn_block, 4)
    do while (associated (pn_data))
       pn_line => parse_node_get_sub_ptr (pn_data, 2)
       pn_code => parse_node_get_sub_ptr (pn_line)
       if (associated (pn_code)) then
          pdg = get_integer_parameter (pn_code)
          pn_mass => parse_node_get_next_ptr (pn_code)
          if (associated (pn_mass)) then
             mass = get_real_parameter (pn_mass)
             call model%set_field_mass (pdg, mass)
          else
             call msg_error ("SLHA: Block MASS: Missing mass value")
          end if
       else
          call msg_error ("SLHA: Block MASS: Missing PDG code")
       end if
       pn_data => parse_node_get_next_ptr (pn_data)
    end do
  end subroutine slha_handle_MASS

  subroutine slha_handle_decays (parse_tree, model)
    type(parse_tree_t), intent(in) :: parse_tree
    type(model_t), intent(inout), target :: model
    type(parse_node_t), pointer :: pn_decay, pn_decay_spec, pn_code, pn_width
    integer :: pdg
    real(default) :: width
    pn_decay => slha_get_first_decay_ptr (parse_tree)
    do while (associated (pn_decay))
       pn_decay_spec => parse_node_get_sub_ptr (pn_decay, 2)
       pn_code => parse_node_get_sub_ptr (pn_decay_spec)
       pdg = get_integer_parameter (pn_code)
       pn_width => parse_node_get_next_ptr (pn_code)
       width = get_real_parameter (pn_width)
       call model%set_field_width (pdg, width)
       pn_decay => slha_get_next_decay_ptr (pn_decay)
    end do
  end subroutine slha_handle_decays

  subroutine slha_handle_matrix_block &
       (parse_tree, block_name, var_prefix, dim1, dim2, model)
    type(parse_tree_t), intent(in) :: parse_tree
    character(*), intent(in) :: block_name, var_prefix
    integer, intent(in) :: dim1, dim2
    type(model_t), intent(inout), target :: model
    type(parse_node_t), pointer :: pn_block
    type(var_list_t), pointer :: var_list
    integer :: i, j
    character(len=len(var_prefix)+2) :: var_name
    var_list => model%get_var_list_ptr ()
    pn_block => slha_get_block_ptr &
         (parse_tree, var_str (block_name), required=.false.)
    if (.not. (associated (pn_block)))  return
    do i = 1, dim1
       do j = 1, dim2
          write (var_name, "(A,I1,I1)")  var_prefix, i, j
          call set_matrix_element (pn_block, i, j, var_str (var_name), var_list)
       end do
    end do
  end subroutine slha_handle_matrix_block

  subroutine slha_handle_ALPHA (parse_tree, model)
    type(parse_tree_t), intent(in) :: parse_tree
    type(model_t), intent(inout), target :: model
    type(parse_node_t), pointer :: pn_block, pn_line, pn_data, pn_item
    type(var_list_t), pointer :: var_list
    real(default) :: al_h
    var_list => model%get_var_list_ptr ()
    pn_block => slha_get_block_ptr &
         (parse_tree, var_str ("ALPHA"), required=.false.)
    if (.not. (associated (pn_block)))  return
    pn_data => parse_node_get_sub_ptr (pn_block, 4)
    pn_line => parse_node_get_sub_ptr (pn_data, 2)
    pn_item => parse_node_get_sub_ptr (pn_line)
    if (associated (pn_item)) then
       al_h = get_real_parameter (pn_item)
       call var_list%set_real (var_str ("al_h"), al_h, &
            is_known=.true., ignore=.true.)
    end if
  end subroutine slha_handle_ALPHA

  subroutine slha_handle_HMIX (parse_tree, model)
    type(parse_tree_t), intent(in) :: parse_tree
    type(model_t), intent(inout), target :: model
    type(parse_node_t), pointer :: pn_block
    type(var_list_t), pointer :: var_list
    var_list => model%get_var_list_ptr ()
    pn_block => slha_get_block_ptr &
         (parse_tree, var_str ("HMIX"), required=.false.)
    if (.not. (associated (pn_block)))  return
    call set_data_item (pn_block, 1, var_str ("mu_h"), var_list)
    call set_data_item (pn_block, 2, var_str ("tanb_h"), var_list)
  end subroutine slha_handle_HMIX

  subroutine slha_handle_NMSSMRUN (parse_tree, model)
    type(parse_tree_t), intent(in) :: parse_tree
    type(model_t), intent(inout), target :: model
    type(parse_node_t), pointer :: pn_block
    type(var_list_t), pointer :: var_list
    var_list => model%get_var_list_ptr ()
    pn_block => slha_get_block_ptr &
         (parse_tree, var_str ("NMSSMRUN"), required=.false.)
    if (.not. (associated (pn_block)))  return
    call set_data_item (pn_block, 1, var_str ("ls"), var_list)
    call set_data_item (pn_block, 2, var_str ("ks"), var_list)
    call set_data_item (pn_block, 3, var_str ("a_ls"), var_list)
    call set_data_item (pn_block, 4, var_str ("a_ks"), var_list)
    call set_data_item (pn_block, 5, var_str ("nmu"), var_list)
    end subroutine slha_handle_NMSSMRUN

  subroutine slha_handle_custom_file (parse_tree, model)
    type(parse_tree_t), intent(in) :: parse_tree
    type(model_t), intent(inout), target :: model

    type(parse_node_t), pointer :: pn_root, pn_block
    type(parse_node_t), pointer :: pn_block_spec, pn_block_name
    type(parse_node_t), pointer :: pn_data, pn_line, pn_code, pn_item
    type(string_t) :: block_name
    integer, dimension(:), allocatable :: block_index
    integer :: n_index, i
    real(default) :: value
    
    pn_root => parse_tree%get_root_ptr ()
    pn_block => pn_root%get_sub_ptr ()
    HANDLE_BLOCKS: do while (associated (pn_block))
       select case (char (pn_block%get_rule_key ()))
       case ("block_def")
          call slha_handle_custom_block (pn_block, model)
       end select
       pn_block => pn_block%get_next_ptr ()
    end do HANDLE_BLOCKS

  end subroutine slha_handle_custom_file
  
  subroutine slha_handle_custom_block (pn_block, model)
    type(parse_node_t), intent(in), target :: pn_block
    type(model_t), intent(inout), target :: model
    
    type(parse_node_t), pointer :: pn_block_spec, pn_block_name
    type(parse_node_t), pointer :: pn_data, pn_line, pn_code, pn_item
    type(string_t) :: block_name
    integer, dimension(:), allocatable :: block_index
    integer :: n_index, i
    real(default) :: value

    pn_block_spec => parse_node_get_sub_ptr (pn_block, 2)
    pn_block_name => parse_node_get_sub_ptr (pn_block_spec)
    select case (char (parse_node_get_rule_key (pn_block_name)))
    case ("block_name")
       block_name = trim (adjustl (upper_case (pn_block_name%get_string ())))
    case ("QNUMBERS")
       block_name = "QNUMBERS"
    end select
    call demangle_keywords (block_name)
    pn_data => pn_block%get_sub_ptr (4)
    HANDLE_LINES: do while (associated (pn_data))
       select case (char (pn_data%get_rule_key ()))
       case ("block_data")
          pn_line => pn_data%get_sub_ptr (2)
          n_index = pn_line%get_n_sub () - 1
          allocate (block_index (n_index))
          pn_code => pn_line%get_sub_ptr ()
          READ_LINE: do i = 1, n_index
             select case (char (pn_code%get_rule_key ()))
             case ("integer");  block_index(i) = pn_code%get_integer ()
             case default
                pn_code => null ()
                exit READ_LINE
             end select
             pn_code => pn_code%get_next_ptr ()
          end do READ_LINE
          if (associated (pn_code)) then
             value = get_real_parameter (pn_code)
             call model%slha_set_par (block_name, block_index, value)
          end if
          deallocate (block_index)
       end select
       pn_data => pn_data%get_next_ptr ()
    end do HANDLE_LINES
    
  end subroutine slha_handle_custom_block
  
  subroutine slha_parse_stream (stream, custom_block_name, parse_tree)
    type(stream_t), intent(inout), target :: stream
    type(string_t), dimension(:), intent(in) :: custom_block_name
    type(parse_tree_t), intent(out) :: parse_tree
    type(ifile_t) :: ifile
    type(lexer_t) :: lexer
    type(stream_t), target :: stream_tmp
    call slha_preprocess (stream, custom_block_name, ifile)
    call stream_init (stream_tmp, ifile)
    call lexer_init_slha (lexer)
    call lexer_assign_stream (lexer, stream_tmp)
    call parse_tree_init (parse_tree, syntax_slha, lexer)
    call lexer_final (lexer)
    call stream_final (stream_tmp)
    call ifile_final (ifile)
  end subroutine slha_parse_stream

  module subroutine slha_parse_file &
       (file, custom_block_name, os_data, parse_tree)
    type(string_t), intent(in) :: file
    type(string_t), dimension(:), intent(in) :: custom_block_name
    type(os_data_t), intent(in) :: os_data
    type(parse_tree_t), intent(out) :: parse_tree
    logical :: exist
    type(string_t) :: filename
    type(stream_t), target :: stream
    call msg_message ("Reading SLHA input file '" // char (file) // "'")
    filename = file
    inquire (file=char(filename), exist=exist)
    if (.not. exist) then
       filename = os_data%whizard_susypath // "/" // file
       inquire (file=char(filename), exist=exist)
       if (.not. exist) then
          call msg_fatal ("SLHA input file '" // char (file) // "' not found")
          return
       end if
    end if
    call stream_init (stream, char (filename))
    call slha_parse_stream (stream, custom_block_name, parse_tree)
    call stream_final (stream)
  end subroutine slha_parse_file

  module subroutine slha_read_file &
       (file, os_data, model, input, spectrum, decays)
    type(string_t), intent(in) :: file
    type(os_data_t), intent(in) :: os_data
    type(model_t), intent(inout), target :: model
    logical, intent(in) :: input, spectrum, decays
    type(string_t), dimension(:), allocatable :: custom_block_name
    type(parse_tree_t) :: parse_tree
    call model%get_custom_slha_blocks (custom_block_name)
    call slha_parse_file (file, custom_block_name, os_data, parse_tree)
    if (associated (parse_tree%get_root_ptr ())) then
       call slha_interpret_parse_tree &
            (parse_tree, model, input, spectrum, decays)
       call parse_tree_final (parse_tree)
       call model%update_parameters ()
    end if
  end subroutine slha_read_file

  module subroutine slha_write_file (file, model, input, spectrum, decays)
    type(string_t), intent(in) :: file
    type(model_t), target, intent(in) :: model
    logical, intent(in) :: input, spectrum, decays
    integer :: mssm_type
    integer :: u
    u = free_unit ()
    call msg_message ("Writing SLHA output file '" // char (file) // "'")
    open (unit=u, file=char(file), action="write", status="replace")
    write (u, "(A)")  "# SUSY Les Houches Accord"
    write (u, "(A)")  "# Output generated by " // trim (VERSION_STRING)
    call slha_write_MODSEL (u, model, mssm_type)
    if (input) then
       call slha_write_SMINPUTS (u, model)
       call slha_write_MINPAR (u, model, mssm_type)
    end if
    if (spectrum) then
       call msg_bug ("SLHA: spectrum output not supported yet")
    end if
    if (decays) then
       call msg_bug ("SLHA: decays output not supported yet")
    end if
    close (u)
  end subroutine slha_write_file

  module subroutine dispatch_slha (var_list, input, spectrum, decays)
    type(var_list_t), intent(inout), target :: var_list
    logical, intent(out) :: input, spectrum, decays
    input = var_list%get_lval (var_str ("?slha_read_input"))
    spectrum = var_list%get_lval (var_str ("?slha_read_spectrum"))
    decays = var_list%get_lval (var_str ("?slha_read_decays"))
  end subroutine dispatch_slha


end submodule slha_interface_s

