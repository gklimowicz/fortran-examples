! Copyright (c) 2012-2013, 2016 Harald Klimach <harald@klimachs.de>
!
! Parts of this file were written by Harald Klimach for
! German Research School of Simulation Sciences and University of Siegen.
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in
! all copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
! IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
! DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
! OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE
! OR OTHER DEALINGS IN THE SOFTWARE.
! **************************************************************************** !

!> A module providing quadruple number input to Lua functions
!!
!! Note that Lua actually only handles double precision, and the numbers are
!! converted accordingly. Thus this is merely a convenience interface, to allow
!! the usage of the functions from this module with quadruple precision numbers.
module aot_quadruple_fun_module
  use flu_binding
  use flu_kinds_module, only: double_k
  use aot_quadruple_top_module, only: quad_k
  use aot_fun_declaration_module, only: aot_fun_type
  use aot_table_module, only: aot_table_from_1Darray

  implicit none

  private

  public :: aot_fun_put

  !> Put an argument into the lua function.
  !!
  !! Arguments have to be in order, first put the first argument then the second
  !! and so on.
  !! Here we add support for quadruple precision numbers
  interface aot_fun_put
    module procedure aot_fun_put_quadruple
    module procedure aot_fun_put_quadruple_v
  end interface aot_fun_put

contains

  !> Put an argument of type double into the list of arguments for the function.
  subroutine aot_fun_put_quadruple(L, fun, arg)
    type(flu_state) :: L !< Handle for the Lua script.

    !> Handle of the function, this argument should be put into.
    type(aot_fun_type) :: fun

    !> Actual argument to hand over to the Lua function.
    real(kind=quad_k), intent(in) :: arg

    real(kind=double_k) :: locarg

    ! Only do something, if the function is actually properly defined.
    if (fun%handle /= 0) then

      locarg = real(arg, kind=double_k)

      ! If the function was executed before this call, it has to be
      ! reset.
      if (fun%arg_count == -1) then
        ! Set the top of the stack to the reference of the function.
        ! Discarding anything above it.
        call flu_settop(L, fun%handle)
        ! Push a copy of the function itself on the stack again, before
        ! adding arguments, to savely survive popping of the function
        ! upon execution.
        call flu_pushvalue(L, fun%handle)
        ! Increase the argument count to 0 again (really start counting
        ! arguments afterwards.
        fun%arg_count = fun%arg_count+1
      end if

      call flu_pushNumber(L, locarg)
      fun%arg_count = fun%arg_count+1
    end if

  end subroutine aot_fun_put_quadruple


  !> Put an array of quadruples into the list of arguments for the
  !! function.
  subroutine aot_fun_put_quadruple_v(L, fun, arg)
    type(flu_state) :: L !< Handle for the Lua script.

    !> Handle of the function, this argument should be put into.
    type(aot_fun_type) :: fun

    !> Actual argument to hand over to the Lua function.
    real(kind=quad_k), intent(in) :: arg(:)

    real(kind=double_k) :: locarg(size(arg))
    integer :: thandle

    ! Only do something, if the function is actually properly defined.
    if (fun%handle /= 0) then

      locarg = real(arg, kind=double_k)

      ! If the function was executed before this call, it has to be
      ! reset.
      if (fun%arg_count == -1) then
        ! Set the top of the stack to the reference of the function.
        ! Discarding anything above it.
        call flu_settop(L, fun%handle)
        ! Push a copy of the function itself on the stack again, before
        ! adding arguments, to savely survive popping of the function
        ! upon execution.
        call flu_pushvalue(L, fun%handle)
        ! Increase the argument count to 0 again (really start counting
        ! arguments afterwards.
        fun%arg_count = fun%arg_count+1
      end if

      call aot_table_from_1Darray(L, thandle, locarg)
      fun%arg_count = fun%arg_count+1
    end if

  end subroutine aot_fun_put_quadruple_v

end module aot_quadruple_fun_module
