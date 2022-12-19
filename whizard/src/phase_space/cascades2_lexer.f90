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

module cascades2_lexer

  use kinds, only: default
  use kinds, only: TC, i8

  implicit none
  private

  public :: dag_token_t
  public :: dag_string_t
  public :: dag_chain_t
  public :: assignment (=)
  public :: operator (//)
  public :: operator (==)
  public :: operator (/=)
  public :: char

  integer, parameter :: PRT_NAME_LEN = 20
  character(len=1), parameter, public :: BACKSLASH_CHAR = "\\"
  character(len=1), parameter :: BLANC_CHAR = " "
  integer, parameter, public :: NEW_LINE_TK = -2
  integer, parameter :: BLANC_SPACE_TK = -1
  integer, parameter :: EMPTY_TK = 0
  integer, parameter, public :: NODE_TK = 1
  integer, parameter, public :: DAG_NODE_TK = 2
  integer, parameter, public :: DAG_OPTIONS_TK = 3
  integer, parameter, public :: DAG_COMBINATION_TK = 4
  integer, parameter, public :: COLON_TK = 11
  integer, parameter, public :: COMMA_TK = 12
  integer, parameter, public :: VERTICAL_BAR_TK = 13
  integer, parameter, public :: OPEN_PAR_TK = 21
  integer, parameter, public :: CLOSED_PAR_TK = 22
  integer, parameter, public :: OPEN_CURLY_TK = 31
  integer, parameter, public :: CLOSED_CURLY_TK = 32


  type :: dag_token_t
     integer :: type = EMPTY_TK
     integer :: char_len = 0
     integer(TC) :: bincode = 0
     character(len=PRT_NAME_LEN) :: particle_name=""
     integer :: index = 0
   contains
       procedure :: init_dag_object_token => dag_token_init_dag_object_token
  end type dag_token_t

  type :: dag_string_t
     integer :: char_len = 0
     type(dag_token_t), dimension(:), allocatable :: t
     type(dag_string_t), pointer :: next => null ()
   contains
       procedure :: clean => dag_string_clean
       procedure :: update_char_len => dag_string_update_char_len
       procedure :: final => dag_string_final
  end type dag_string_t

  type :: dag_chain_t
     integer :: char_len = 0
     integer :: t_size = 0
     type(dag_string_t), pointer :: first => null ()
     type(dag_string_t), pointer :: last => null ()
   contains
       procedure :: append => dag_chain_append_string
       procedure :: compress => dag_chain_compress
       procedure :: final => dag_chain_final
  end type dag_chain_t


  interface assignment (=)
     module procedure dag_token_assign_from_char_string
     module procedure dag_token_assign_from_dag_token
     module procedure dag_string_assign_from_dag_token
     module procedure dag_string_assign_from_char_string
     module procedure dag_string_assign_from_dag_string
     module procedure dag_string_assign_from_dag_token_array
  end interface assignment (=)

  interface operator (//)
     module procedure concat_dag_token_dag_token
     module procedure concat_dag_string_dag_token
     module procedure concat_dag_token_dag_string
     module procedure concat_dag_string_dag_string
  end interface operator (//)

  interface operator (==)
     module procedure dag_token_eq_dag_token
     module procedure dag_string_eq_dag_string
     module procedure dag_token_eq_dag_string
     module procedure dag_string_eq_dag_token
     module procedure dag_token_eq_char_string
     module procedure char_string_eq_dag_token
     module procedure dag_string_eq_char_string
     module procedure char_string_eq_dag_string
  end interface operator (==)

  interface operator (/=)
     module procedure dag_token_ne_dag_token
     module procedure dag_string_ne_dag_string
     module procedure dag_token_ne_dag_string
     module procedure dag_string_ne_dag_token
     module procedure dag_token_ne_char_string
     module procedure char_string_ne_dag_token
     module procedure dag_string_ne_char_string
     module procedure char_string_ne_dag_string
  end interface operator (/=)

  interface char
     module procedure char_dag_token
     module procedure char_dag_string
  end interface char


  interface
    module subroutine dag_token_init_dag_object_token (dag_token, type, index)
      class(dag_token_t), intent(out) :: dag_token
      integer, intent(in) :: index
      integer :: type
    end subroutine dag_token_init_dag_object_token
    elemental module subroutine dag_token_assign_from_char_string &
         (dag_token, char_string)
      type(dag_token_t), intent(out) :: dag_token
      character(len=*), intent(in) :: char_string
    end subroutine dag_token_assign_from_char_string
    elemental module subroutine dag_token_assign_from_dag_token &
         (token_out, token_in)
      type(dag_token_t), intent(out) :: token_out
      type(dag_token_t), intent(in) :: token_in
    end subroutine dag_token_assign_from_dag_token
    elemental module subroutine dag_string_assign_from_dag_token &
         (dag_string, dag_token)
      type(dag_string_t), intent(out) :: dag_string
      type(dag_token_t), intent(in) :: dag_token
    end subroutine dag_string_assign_from_dag_token
    module subroutine dag_string_assign_from_dag_token_array &
         (dag_string, dag_token)
      type(dag_string_t), intent(out) :: dag_string
      type(dag_token_t), dimension(:), intent(in) :: dag_token
    end subroutine dag_string_assign_from_dag_token_array
    elemental module subroutine dag_string_assign_from_char_string &
         (dag_string, char_string)
      type(dag_string_t), intent(out) :: dag_string
      character(len=*), intent(in) :: char_string
    end subroutine dag_string_assign_from_char_string
    elemental module subroutine dag_string_assign_from_dag_string &
         (string_out, string_in)
      type(dag_string_t), intent(out) :: string_out
      type(dag_string_t), intent(in) :: string_in
    end subroutine dag_string_assign_from_dag_string
    module function concat_dag_token_dag_token &
         (token1, token2) result (res_string)
      type(dag_token_t), intent(in) :: token1, token2
      type(dag_string_t) :: res_string
    end function concat_dag_token_dag_token
    module function concat_dag_string_dag_token &
         (dag_string, dag_token) result (res_string)
      type(dag_string_t), intent(in) :: dag_string
      type(dag_token_t), intent(in) :: dag_token
      type(dag_string_t) :: res_string
    end function concat_dag_string_dag_token
    module function concat_dag_token_dag_string &
         (dag_token, dag_string) result (res_string)
      type(dag_token_t), intent(in) :: dag_token
      type(dag_string_t), intent(in) :: dag_string
      type(dag_string_t) :: res_string
      integer :: t_size
    end function concat_dag_token_dag_string
    module function concat_dag_string_dag_string &
         (string1, string2) result (res_string)
      type(dag_string_t), intent(in) :: string1, string2
      type(dag_string_t) :: res_string
    end function concat_dag_string_dag_string
    elemental module function dag_token_eq_dag_token &
         (token1, token2) result (flag)
      type(dag_token_t), intent(in) :: token1, token2
      logical :: flag
    end function dag_token_eq_dag_token
    elemental module function dag_string_eq_dag_string &
         (string1, string2) result (flag)
      type(dag_string_t), intent(in) :: string1, string2
      logical :: flag
    end function dag_string_eq_dag_string
    elemental module function dag_token_eq_dag_string &
         (dag_token, dag_string) result (flag)
      type(dag_token_t), intent(in) :: dag_token
      type(dag_string_t), intent(in) :: dag_string
      logical :: flag
    end function dag_token_eq_dag_string
    elemental module function dag_string_eq_dag_token &
         (dag_string, dag_token) result (flag)
      type(dag_token_t), intent(in) :: dag_token
      type(dag_string_t), intent(in) :: dag_string
      logical :: flag
    end function dag_string_eq_dag_token
    elemental module function dag_token_eq_char_string &
         (dag_token, char_string) result (flag)
      type(dag_token_t), intent(in) :: dag_token
      character(len=*), intent(in) :: char_string
      logical :: flag
    end function dag_token_eq_char_string
    elemental module function char_string_eq_dag_token &
         (char_string, dag_token) result (flag)
      type(dag_token_t), intent(in) :: dag_token
      character(len=*), intent(in) :: char_string
      logical :: flag
    end function char_string_eq_dag_token
    elemental module function dag_string_eq_char_string &
         (dag_string, char_string) result (flag)
      type(dag_string_t), intent(in) :: dag_string
      character(len=*), intent(in) :: char_string
      logical :: flag
    end function dag_string_eq_char_string
    elemental module function char_string_eq_dag_string &
         (char_string, dag_string) result (flag)
      type(dag_string_t), intent(in) :: dag_string
      character(len=*), intent(in) :: char_string
      logical :: flag
    end function char_string_eq_dag_string
    elemental module function dag_token_ne_dag_token &
         (token1, token2) result (flag)
      type(dag_token_t), intent(in) :: token1, token2
      logical :: flag
    end function dag_token_ne_dag_token
    elemental module function dag_string_ne_dag_string &
         (string1, string2) result (flag)
      type(dag_string_t), intent(in) :: string1, string2
      logical :: flag
    end function dag_string_ne_dag_string
    elemental module function dag_token_ne_dag_string &
         (dag_token, dag_string) result (flag)
      type(dag_token_t), intent(in) :: dag_token
      type(dag_string_t), intent(in) :: dag_string
      logical :: flag
    end function dag_token_ne_dag_string
    elemental module function dag_string_ne_dag_token &
         (dag_string, dag_token) result (flag)
      type(dag_token_t), intent(in) :: dag_token
      type(dag_string_t), intent(in) :: dag_string
      logical :: flag
    end function dag_string_ne_dag_token
    elemental module function dag_token_ne_char_string &
         (dag_token, char_string) result (flag)
      type(dag_token_t), intent(in) :: dag_token
      character(len=*), intent(in) :: char_string
      logical :: flag
    end function dag_token_ne_char_string
    elemental module function char_string_ne_dag_token &
         (char_string, dag_token) result (flag)
      type(dag_token_t), intent(in) :: dag_token
      character(len=*), intent(in) :: char_string
      logical :: flag
    end function char_string_ne_dag_token
    elemental module function dag_string_ne_char_string &
         (dag_string, char_string) result (flag)
      type(dag_string_t), intent(in) :: dag_string
      character(len=*), intent(in) :: char_string
      logical :: flag
    end function dag_string_ne_char_string
    elemental module function char_string_ne_dag_string &
         (char_string, dag_string) result (flag)
      type(dag_string_t), intent(in) :: dag_string
      character(len=*), intent(in) :: char_string
      logical :: flag
    end function char_string_ne_dag_string
    pure module function char_dag_token (dag_token) result (char_string)
      type(dag_token_t), intent(in) :: dag_token
      character (dag_token%char_len) :: char_string
    end function char_dag_token
    pure module function char_dag_string (dag_string) result (char_string)
      type(dag_string_t), intent(in) :: dag_string
      character (dag_string%char_len) :: char_string
    end function char_dag_string
    module subroutine dag_string_clean (dag_string)
      class(dag_string_t), intent(inout) :: dag_string
    end subroutine dag_string_clean
    module subroutine dag_string_update_char_len (dag_string)
      class(dag_string_t), intent(inout) :: dag_string
    end subroutine dag_string_update_char_len
    module subroutine dag_chain_append_string (dag_chain, char_string)
      class(dag_chain_t), intent(inout) :: dag_chain
      character(len=*), intent(in) :: char_string
    end subroutine dag_chain_append_string
    module subroutine dag_chain_compress (dag_chain)
      class(dag_chain_t), intent(inout) :: dag_chain
    end subroutine dag_chain_compress
    module subroutine dag_string_final (dag_string)
      class(dag_string_t), intent(inout) :: dag_string
    end subroutine dag_string_final
    module subroutine dag_chain_final (dag_chain)
      class(dag_chain_t), intent(inout) :: dag_chain
    end subroutine dag_chain_final
  end interface

end module cascades2_lexer

