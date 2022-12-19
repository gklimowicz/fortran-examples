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

module slha_interface

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use os_interface
  use lexers
  use parser
  use variables
  use models

  implicit none
  private

  public :: syntax_slha_init
  public :: syntax_slha_final
  public :: syntax_slha_write
  public :: lexer_init_slha
  public :: slha_interpret_parse_tree
  public :: slha_parse_file
  public :: slha_read_file
  public :: slha_write_file
  public :: dispatch_slha

  type :: str_entry_t
     type(string_t) :: str
     type(str_entry_t), pointer :: next => null ()
  end type str_entry_t


  interface
    module subroutine syntax_slha_init ()
    end subroutine syntax_slha_init
    module subroutine syntax_slha_final ()
    end subroutine syntax_slha_final
    module subroutine syntax_slha_write (unit)
      integer, intent(in), optional :: unit
    end subroutine syntax_slha_write
    module subroutine lexer_init_slha (lexer)
      type(lexer_t), intent(out) :: lexer
    end subroutine lexer_init_slha
    module subroutine slha_interpret_parse_tree &
         (parse_tree, model, input, spectrum, decays)
      type(parse_tree_t), intent(in) :: parse_tree
      type(model_t), intent(inout), target :: model
      logical, intent(in) :: input, spectrum, decays
    end subroutine slha_interpret_parse_tree
    module subroutine slha_parse_file &
         (file, custom_block_name, os_data, parse_tree)
      type(string_t), intent(in) :: file
      type(string_t), dimension(:), intent(in) :: custom_block_name
      type(os_data_t), intent(in) :: os_data
      type(parse_tree_t), intent(out) :: parse_tree
    end subroutine slha_parse_file
    module subroutine slha_read_file &
         (file, os_data, model, input, spectrum, decays)
      type(string_t), intent(in) :: file
      type(os_data_t), intent(in) :: os_data
      type(model_t), intent(inout), target :: model
      logical, intent(in) :: input, spectrum, decays
    end subroutine slha_read_file
    module subroutine slha_write_file (file, model, input, spectrum, decays)
      type(string_t), intent(in) :: file
      type(model_t), target, intent(in) :: model
      logical, intent(in) :: input, spectrum, decays
    end subroutine slha_write_file
    module subroutine dispatch_slha (var_list, input, spectrum, decays)
      type(var_list_t), intent(inout), target :: var_list
      logical, intent(out) :: input, spectrum, decays
    end subroutine dispatch_slha
  end interface

  save

end module slha_interface
