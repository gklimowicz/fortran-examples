! Copyright (c) 2011-2016 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2012 James Spencer <j.spencer@imperial.ac.uk>
! Copyright (c) 2018 Raphael Haupt <Raphael.Haupt@student.uni-siegen.de>
!
! Parts of this file were written by Harald Klimach for
! German Research School of Simulation Sciences and University of Siegen.
! Parts of this file were written by Raphael Haupt for University of Siegen.
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

!> This module provides access to Lua functions
!!
!! Intented usage:
!!
!! - First open a function with [[aot_fun_open]].
!! - Then put required parameters into it with [[aot_fun_put]].
!! - Execute the function with [[aot_fun_do]].
!! - Retrieve the possibly multiple results with [[aot_top_get_val]].
!!   If there are multiple results to be retrieved from the function
!!   repeat calling [[aot_top_get_val]] for each of them. Keep in mind that they
!!   will be in reversed order on the stack!
!! - Repeat putting and retrieving as needed (for multiple function
!!   evaluations).
!! - Close the function finally with [[aot_fun_close]].
module aot_fun_module
  use flu_binding
  use flu_kinds_module, only: double_k, single_k
  use aot_fun_declaration_module, only: aot_fun_type
  use aot_table_module, only: aot_table_push, aot_table_from_1Darray
  use aot_top_module, only: aot_err_handler
  use aot_references_module, only: aot_reference_to_top

  ! Include quadruple precision interfaces if available
  use aot_quadruple_fun_module

  ! Support for extended double precision
  use aot_extdouble_fun_module

  implicit none

  private

  public :: aot_fun_type, aot_fun_open, aot_fun_close, aot_fun_put, aot_fun_do
  public :: aot_fun_top
  public :: aot_fun_id

  !> Open a Lua function for evaluation.
  !!
  !! After it is opened, arguments might be put into the function, and it might
  !! be executed.
  !! Execution might be repeated for an arbitrary number of iterations, to
  !! retrieve more than one evaluation of a single function, before closing it
  !! again with [[aot_fun_close]].
  interface aot_fun_open
    module procedure aot_fun_table
    module procedure aot_fun_ref
  end interface aot_fun_open

  !> Put an argument into the lua function.
  !!
  !! Arguments have to be in order, first put the first argument then the second
  !! and so on.
  !! Currently only real number arguments are supported.
  interface aot_fun_put
    module procedure aot_fun_put_top
    module procedure aot_fun_put_double
    module procedure aot_fun_put_single
    module procedure aot_fun_put_double_v
    module procedure aot_fun_put_single_v
  end interface aot_fun_put


contains


  !> Return the stack of the top as a function.
  !!
  !! If it actually is not a Lua function, the returned handle will be 0.
  function aot_fun_top(L) result(fun)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle to the function on the top of the stack.
    type(aot_fun_type) :: fun

    fun%handle = 0
    fun%arg_count = 0
    if (flu_isFunction(L, -1)) then
      ! Keep a handle to this function.
      fun%handle = flu_gettop(L)
      fun%id = flu_topointer(L, -1)
      ! Push a copy of the function right after it, the function will
      ! be popped from the stack upon execution. Thus, this copy is
      ! used to ensure the reference to the function is kept across
      ! several executions of the function.
      call flu_pushvalue(L, -1)
    end if
  end function aot_fun_top

  !> Get a function defined as component of a table.
  !!
  !! Functions in tables might be retrieved by position or key.
  !! If both optional parameters are provided, the key is attempted to be read
  !! first. Only when that fails, the position will be tested.
  subroutine aot_fun_table(L, parent, fun, key, pos)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle to the table to look in for the function.
    integer, intent(in), optional :: parent

    !> Returned handle, providing access to the function.
    type(aot_fun_type), intent(out) :: fun

    !> Name of the function to look up in the table.
    character(len=*), intent(in), optional :: key

    !> Position of the function to look up in the table.
    integer, intent(in), optional :: pos
    call aot_table_push(L, parent, key, pos)
    fun = aot_fun_top(L)
  end subroutine aot_fun_table


  !> Get a function from a previously defned Lua reference.
  !!
  !! Use a previously (with [[aot_reference_for]]) defined reference to get a
  !! function.
  subroutine aot_fun_ref(L, fun, ref)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Returned handle, providing access to the function.
    type(aot_fun_type), intent(out) :: fun

    !> Lua reference to the function.
    integer, intent(in) :: ref

    call aot_reference_to_top(L, ref)
    fun = aot_fun_top(L)
  end subroutine aot_fun_ref


  !> Close the function again (pop everything above from the stack).
  subroutine aot_fun_close(L, fun)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle to the function to close.
    type(aot_fun_type) :: fun

    if (fun%handle > 0) call flu_settop(L, fun%handle-1)
    fun%handle = 0
    fun%id = 0
    fun%arg_count = 0
  end subroutine aot_fun_close


  !> Put the top of the stack as argument into the list of arguments for the
  !! function.
  subroutine aot_fun_put_top(L, fun)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle of the function, this argument should be put into.
    type(aot_fun_type) :: fun

    integer :: curtop

    ! Only do something, if the function is actually properly defined.
    if (fun%handle /= 0) then

      ! Get position of current top of the stack.
      curtop = flu_gettop(L)

      ! If the function was executed before this call, it has to be
      ! reset.
      if (fun%arg_count == -1) then
        ! Only procede, if curtop is exactly one above the function reference,
        ! that is after executing the function previously, only one item was
        ! put into the stack, which should now be used as an argument.
        if (curtop == fun%handle+1) then
          ! Push a copy of the function itself on the stack again, before
          ! adding arguments, to savely survive popping of the function
          ! upon execution. (insert this copy before the already added argument)
          call flu_insert(L, fun%handle+1)
          ! Increase the argument count to 0 again (really start counting
          ! arguments afterwards.
          fun%arg_count = fun%arg_count+1
          curtop = curtop + 1
        end if
      end if

      ! Only proceed, if the current top is actually a new argument (that is, it
      ! is especially not the function copy at fun%handle + 1 itself).
      if ((curtop - fun%arg_count) == (fun%handle + 2)) then
        fun%arg_count = fun%arg_count+1
      end if
    end if

  end subroutine aot_fun_put_top


  !> Put an argument of type double into the list of arguments for the function.
  subroutine aot_fun_put_double(L, fun, arg)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle of the function, this argument should be put into.
    type(aot_fun_type) :: fun

    !> Actual argument to hand over to the Lua function.
    real(kind=double_k), intent(in) :: arg

    ! Only do something, if the function is actually properly defined.
    if (fun%handle /= 0) then

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

      call flu_pushNumber(L, arg)
      fun%arg_count = fun%arg_count+1
    end if

  end subroutine aot_fun_put_double


  !> Put an argument of type single into the list of arguments for the function.
  subroutine aot_fun_put_single(L, fun, arg)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle of the function, this argument should be put into.
    type(aot_fun_type) :: fun

    !> Actual argument to hand over to the Lua function.
    real(kind=single_k), intent(in) :: arg

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

  end subroutine aot_fun_put_single


  !> Put an array of doubles into the list of arguments for the function.
  subroutine aot_fun_put_double_v(L, fun, arg)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle of the function, this argument should be put into.
    type(aot_fun_type) :: fun

    !> Actual argument to hand over to the Lua function.
    real(kind=double_k), intent(in) :: arg(:)

    integer :: thandle

    ! Only do something, if the function is actually properly defined.
    if (fun%handle /= 0) then

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

      call aot_table_from_1Darray(L, thandle, arg)
      fun%arg_count = fun%arg_count+1
    end if

  end subroutine aot_fun_put_double_v


  !> Put an array of singles into the list of arguments for the function.
  subroutine aot_fun_put_single_v(L, fun, arg)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle of the function, this argument should be put into.
    type(aot_fun_type) :: fun

    !> Actual argument to hand over to the Lua function.
    real(kind=single_k), intent(in) :: arg(:)

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

  end subroutine aot_fun_put_single_v




  !> Execute a given function and put its results on the stack, where it is
  !! retrievable with [[aot_top_get_val]].
  !!
  !! The optional arguments ErrCode and ErrString provide some feedback on the
  !! success of the function execution.
  !! If none of them are in the argument list, the execution of the application
  !! will be stopped, and the error will be printed to the standard output.
  !! You have to provide the number of results to obtain in nresults. Keep in
  !! mind, that multiple results have to obtained in reverse order from the
  !! stack.
  !!
  !! @note You might want to return multiple values as a single argument in a
  !!       table instead of several single values.
  subroutine aot_fun_do(L, fun, nresults, ErrCode, ErrString)
    type(flu_state) :: L !! Handle for the Lua script.

    !> Handle to the function to execute.
    type(aot_fun_type) :: fun

    !> Number of resulting values the caller wants to obtain from the Lua
    !! function.
    integer, intent(in) :: nresults

    !> Error code returned by Lua during execution of the function.
    integer, intent(out), optional :: ErrCode

    !> Obtained error string from the Lua stack if an error occured.
    character(len=*), intent(out), optional :: ErrString

    integer :: err

    if (fun%handle /= 0) then
      err = flu_pcall(L, fun%arg_count, nresults, 0)
      call aot_err_handler(L=L, err=err, msg="Failed aot_fun_do! ", &
        &                  ErrCode = ErrCode, ErrString = ErrString)
      fun%arg_count = -1
    end if
  end subroutine aot_fun_do


  !> A string identifying the function uniquely in the Lua script.
  function aot_fun_id(fun) result(id)
    !> Function to identify.
    type(aot_fun_type), intent(in) :: fun

    !> Identification of the function as a string.
    character(len=32) :: id

    character(len=32) :: tmp

    write(tmp,'(i0)') fun%id
    id = adjustl(tmp)
  end function aot_fun_id

end module aot_fun_module
