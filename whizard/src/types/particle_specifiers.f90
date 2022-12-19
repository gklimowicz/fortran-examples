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

module particle_specifiers

  use iso_varying_string, string_t => varying_string

  implicit none
  private

  public :: prt_expr_t
  public :: prt_spec_t
  public :: prt_spec_write
  public :: prt_spec_read
  public :: new_prt_spec
  public :: prt_spec_list_t
  public :: prt_spec_sum_t

  type, abstract :: prt_spec_expr_t
   contains
     procedure (prt_spec_expr_to_string), deferred :: to_string
     procedure (prt_spec_expr_expand_sub), deferred :: expand_sub
  end type prt_spec_expr_t

  type :: prt_expr_t
     class(prt_spec_expr_t), allocatable :: x
   contains
     procedure :: to_string => prt_expr_to_string
     procedure :: init_spec => prt_expr_init_spec
     procedure :: init_list => prt_expr_init_list
     procedure :: init_sum => prt_expr_init_sum
     procedure :: get_n_terms => prt_expr_get_n_terms
     procedure :: term_to_array => prt_expr_term_to_array
     procedure :: expand => prt_expr_expand
  end type prt_expr_t

  type, extends (prt_spec_expr_t) :: prt_spec_t
     private
     type(string_t) :: name
     logical :: polarized = .false.
     type(string_t), dimension(:), allocatable :: decay
   contains
     procedure :: get_name => prt_spec_get_name
     procedure :: to_string => prt_spec_to_string
     procedure :: is_polarized => prt_spec_is_polarized
     procedure :: is_unstable => prt_spec_is_unstable
     procedure :: get_n_decays => prt_spec_get_n_decays
     procedure :: get_decays => prt_spec_get_decays
     procedure :: expand_sub => prt_spec_expand_sub
  end type prt_spec_t

  type, extends (prt_spec_expr_t) :: prt_spec_list_t
     type(prt_expr_t), dimension(:), allocatable :: expr
   contains
     procedure :: to_string => prt_spec_list_to_string
     procedure :: flatten => prt_spec_list_flatten
     procedure :: expand_sub => prt_spec_list_expand_sub
  end type prt_spec_list_t

  type, extends (prt_spec_expr_t) :: prt_spec_sum_t
     type(prt_expr_t), dimension(:), allocatable :: expr
   contains
     procedure :: to_string => prt_spec_sum_to_string
     procedure :: flatten => prt_spec_sum_flatten
     procedure :: expand_sub => prt_spec_sum_expand_sub
  end type prt_spec_sum_t


  abstract interface
     function prt_spec_expr_to_string (object) result (string)
       import
       class(prt_spec_expr_t), intent(in) :: object
       type(string_t) :: string
     end function prt_spec_expr_to_string
  end interface

  abstract interface
     subroutine prt_spec_expr_expand_sub (object)
       import
       class(prt_spec_expr_t), intent(inout) :: object
     end subroutine prt_spec_expr_expand_sub
  end interface

  interface prt_spec_write
     module procedure prt_spec_write1
     module procedure prt_spec_write2
  end interface prt_spec_write
  interface prt_spec_read
     module procedure prt_spec_read1
     module procedure prt_spec_read2
  end interface prt_spec_read
  interface new_prt_spec
     module procedure new_prt_spec_
     module procedure new_prt_spec_polarized
     module procedure new_prt_spec_unstable
  end interface new_prt_spec

  interface
    recursive module function prt_expr_to_string (object) result (string)
      class(prt_expr_t), intent(in) :: object
      type(string_t) :: string
    end function prt_expr_to_string
    module function prt_expr_get_n_terms (object) result (n)
      class(prt_expr_t), intent(in) :: object
      integer :: n
    end function prt_expr_get_n_terms
    recursive module subroutine prt_expr_term_to_array (object, array, i)
      class(prt_expr_t), intent(in) :: object
      type(prt_spec_t), dimension(:), intent(inout), allocatable :: array
      integer, intent(in) :: i
    end subroutine prt_expr_term_to_array
    module subroutine prt_spec_write1 (object, unit, advance)
      type(prt_spec_t), intent(in) :: object
      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: advance
    end subroutine prt_spec_write1
    module subroutine prt_spec_write2 (prt_spec, unit, advance)
      type(prt_spec_t), dimension(:), intent(in) :: prt_spec
      integer, intent(in), optional :: unit
      character(len=*), intent(in), optional :: advance
    end subroutine prt_spec_write2
    pure module subroutine prt_spec_read1 (prt_spec, string)
      type(prt_spec_t), intent(out) :: prt_spec
      type(string_t), intent(in) :: string
    end subroutine prt_spec_read1
    pure module subroutine prt_spec_read2 (prt_spec, string)
      type(prt_spec_t), dimension(:), intent(out), allocatable :: prt_spec
      type(string_t), intent(in) :: string
    end subroutine prt_spec_read2
    elemental module function new_prt_spec_ (name) result (prt_spec)
      type(string_t), intent(in) :: name
      type(prt_spec_t) :: prt_spec
    end function new_prt_spec_
    elemental module function new_prt_spec_polarized (name, polarized) result (prt_spec)
      type(string_t), intent(in) :: name
      logical, intent(in) :: polarized
      type(prt_spec_t) :: prt_spec
    end function new_prt_spec_polarized
    pure module function new_prt_spec_unstable (name, decay) result (prt_spec)
      type(string_t), intent(in) :: name
      type(string_t), dimension(:), intent(in) :: decay
      type(prt_spec_t) :: prt_spec
    end function new_prt_spec_unstable
    elemental module function prt_spec_get_name (prt_spec) result (name)
      class(prt_spec_t), intent(in) :: prt_spec
      type(string_t) :: name
    end function prt_spec_get_name
    module function prt_spec_to_string (object) result (string)
      class(prt_spec_t), intent(in) :: object
      type(string_t) :: string
    end function prt_spec_to_string
    elemental module function prt_spec_is_polarized (prt_spec) result (flag)
      class(prt_spec_t), intent(in) :: prt_spec
      logical :: flag
    end function prt_spec_is_polarized
    elemental module function prt_spec_is_unstable (prt_spec) result (flag)
      class(prt_spec_t), intent(in) :: prt_spec
      logical :: flag
    end function prt_spec_is_unstable
    elemental module function prt_spec_get_n_decays (prt_spec) result (n)
      class(prt_spec_t), intent(in) :: prt_spec
      integer :: n
    end function prt_spec_get_n_decays
    module subroutine prt_spec_get_decays (prt_spec, decay)
      class(prt_spec_t), intent(in) :: prt_spec
      type(string_t), dimension(:), allocatable, intent(out) :: decay
    end subroutine prt_spec_get_decays
    module subroutine prt_spec_expand_sub (object)
      class(prt_spec_t), intent(inout) :: object
    end subroutine prt_spec_expand_sub
    recursive module function prt_spec_list_to_string (object) result (string)
      class(prt_spec_list_t), intent(in) :: object
      type(string_t) :: string
    end function prt_spec_list_to_string
    module subroutine prt_spec_list_flatten (object)
      class(prt_spec_list_t), intent(inout) :: object
    end subroutine prt_spec_list_flatten
    recursive module subroutine prt_spec_list_expand_sub (object)
      class(prt_spec_list_t), intent(inout) :: object
    end subroutine prt_spec_list_expand_sub
    recursive module function prt_spec_sum_to_string (object) result (string)
      class(prt_spec_sum_t), intent(in) :: object
      type(string_t) :: string
    end function prt_spec_sum_to_string
    module subroutine prt_spec_sum_flatten (object)
      class(prt_spec_sum_t), intent(inout) :: object
    end subroutine prt_spec_sum_flatten
    recursive module subroutine prt_spec_sum_expand_sub (object)
      class(prt_spec_sum_t), intent(inout) :: object
    end subroutine prt_spec_sum_expand_sub
  end interface

contains

  subroutine prt_expr_init_spec (object, spec)
    class(prt_expr_t), intent(out) :: object
    type(prt_spec_t), intent(in) :: spec
    allocate (prt_spec_t :: object%x)
    select type (x => object%x)
    type is (prt_spec_t)
       x = spec
    end select
  end subroutine prt_expr_init_spec

  subroutine prt_expr_init_list (object, n)
    class(prt_expr_t), intent(out) :: object
    integer, intent(in) :: n
    allocate (prt_spec_list_t :: object%x)
    select type (x => object%x)
    type is (prt_spec_list_t)
       allocate (x%expr (n))
    end select
  end subroutine prt_expr_init_list

  subroutine prt_expr_init_sum (object, n)
    class(prt_expr_t), intent(out) :: object
    integer, intent(in) :: n
    allocate (prt_spec_sum_t :: object%x)
    select type (x => object%x)
    type is (prt_spec_sum_t)
       allocate (x%expr (n))
    end select
  end subroutine prt_expr_init_sum

  subroutine distribute_prt_spec_list (object)
    class(prt_spec_expr_t), intent(inout), allocatable :: object
    class(prt_spec_expr_t), allocatable :: new_object
    integer, dimension(:), allocatable :: n, ii
    integer :: k, n_expr, n_terms, i_term
    select type (object)
    type is (prt_spec_list_t)
       n_expr = size (object%expr)
       allocate (n (n_expr), source = 1)
       allocate (ii (n_expr), source = 1)
       do k = 1, size (object%expr)
          select type (y => object%expr(k)%x)
          type is (prt_spec_sum_t)
             n(k) = size (y%expr)
          end select
       end do
       n_terms = product (n)
       if (n_terms > 1) then
          allocate (prt_spec_sum_t :: new_object)
          select type (new_object)
          type is (prt_spec_sum_t)
             allocate (new_object%expr (n_terms))
             do i_term = 1, n_terms
                allocate (prt_spec_list_t :: new_object%expr(i_term)%x)
                select type (x => new_object%expr(i_term)%x)
                type is (prt_spec_list_t)
                   allocate (x%expr (n_expr))
                   do k = 1, n_expr
                      select type (y => object%expr(k)%x)
                      type is (prt_spec_sum_t)
                         x%expr(k) = y%expr(ii(k))
                      class default
                         x%expr(k) = object%expr(k)
                      end select
                   end do
                end select
                INCR_INDEX: do k = n_expr, 1, -1
                   if (ii(k) < n(k)) then
                      ii(k) = ii(k) + 1
                      exit INCR_INDEX
                   else
                      ii(k) = 1
                   end if
                end do INCR_INDEX
             end do
          end select
       end if
    end select
    if (allocated (new_object)) call move_alloc (from = new_object, to = object)
  end subroutine distribute_prt_spec_list

  recursive subroutine prt_expr_expand (expr)
    class(prt_expr_t), intent(inout) :: expr
    if (allocated (expr%x)) then
       call distribute_prt_spec_list (expr%x)
       call expr%x%expand_sub ()
       select type (x => expr%x)
       type is (prt_spec_list_t)
          call x%flatten ()
       type is (prt_spec_sum_t)
          call x%flatten ()
       end select
    end if
  end subroutine prt_expr_expand


end module particle_specifiers
