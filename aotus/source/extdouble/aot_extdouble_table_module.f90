! Copyright (c) 2012-2013, 2016 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2018 Nick Papior <nickpapior@gmail.com>
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

module aot_extdouble_table_module
  use flu_binding
  use flu_kinds_module, only: double_k
  use aot_err_module, only: aoterr_Fatal, aoterr_NonExistent
  use aot_extdouble_top_module, only: xdble_k
  use aot_top_module, only: aot_top_get_val
  use aot_table_ops_module, only: aot_table_open, aot_table_close, &
    &                             aot_table_length, aot_table_first, &
    &                             aot_table_top, aot_table_push
  implicit none

  private

  public :: aot_table_get_val, aot_table_set_val, aot_table_from_1Darray, &
    &       aot_get_val

  !> Get a value from a table.
  !!
  !! First the given key is looked up, if this fails, the value
  !! at the given position is looked up, and if this also fails,
  !! the default value is returned.
  !! Positional addressing is only valid, as long,
  !! as no value was provided by an explicit key
  !! in the list before the entry in question.
  interface aot_table_get_val
    module procedure get_table_extdouble
  end interface

  !> Set a value in a table.
  !!
  !! The given value will be put at the entry named by key into the table
  !! provided in thandle.
  !! Alternatively you can also put the value by position into the table by
  !! providing the pos argument.
  !! If both, pos and key are provided, the key will be used.
  !! Though, both of them are optional, at least one of them has to be provided.
  interface aot_table_set_val
    module procedure set_table_extdouble
  end interface

  !> Get a value from a table.
  !!
  !! First the given key is looked up, if this fails, the value
  !! at the given position is looked up, and if this also fails,
  !! the default value is returned.
  !! Positional addressing is only valid, as long,
  !! as no value was provided by an explicit key
  !! in the list before the entry in question.
  !!
  !! The interface to access table values looks like:
  !! `call aot_get_val(val, errCode, L, thandle, key, pos, default)`.
  !! Position pos and key are both optional, but one of them has to be provided.
  !! If both are provided the key takes precedence over the pos, and the pos
  !! will only be tried if the access to the key fails.
  !! See for example get_table_real() for a more detailed
  !! description of the parameters.
  !!
  !! Note that positional addressing only works intuitively as long as there
  !! have been no entries specified by keys in the table.
  !! This kind of resembles the behavior of Fortran interfaces with named or
  !! unnamed arguments, as soon as you provide a name, all following arguments
  !! have to be given by key also.
  !! Just stick to this rule for the Lua tables as well to avoid too much
  !! headache.
  !!
  !! The reason for this is, that positional addressing in Lua refers only to
  !! the unnamed entries of the tables.
  interface aot_get_val
    module procedure get_table_extdouble
  end interface

  !> This interface enables the simple creation of uniform one dimensional
  !! arrays as tables in the Lua context.
  !!
  !! It takes an one dimensional array of values and returns a thandle to
  !! identify the newly generated table.
  interface aot_table_from_1Darray
    module procedure create_1Darray_extdouble
  end interface

contains

  !> Retrieve a extdouble precision real value from a table.
  !!
  !! NOTE that Lua actually only provides double precision numbers, and this
  !! interface is merely a convenience for Fortran implementations with
  !! extdouble precision real numbers.
  subroutine get_table_extdouble(val, ErrCode, L, thandle, key, pos, &
    &                            default)
    type(flu_State) :: L !< Handle to the Lua script.

    !> Handle to the table to look the value up in.
    integer, intent(in), optional :: thandle

    !> Value of the table entry if it exists.
    real(kind=xdble_k), intent(out) :: val

    !> Error code to indicate what kind of problem might have occured.
    integer, intent(out) :: ErrCode

    !> Name of the entry to look for.
    !!
    !! Key and pos are both optional, however at least one of them has to be
    !! supplied.
    !! The key takes precedence over the pos if both are given.
    character(len=*), intent(in), optional :: key

    !> Position of the entry to look for in the table.
    !!
    !! It allows the access to unnamed arrays in the Lua tables.
    integer, intent(in), optional :: pos

    !> Some default value, that should be used, if the variable is not set in
    !! the Lua script.
    real(kind=xdble_k), intent(in), optional :: default

    logical :: valid_args
    integer :: toptype

    valid_args = .true.
    if (present(thandle)) then
      call aot_table_push(L=L, thandle=thandle, &
        &                 key=key, pos=pos      )
    else
      if (present(key)) then
        toptype = flu_getglobal(L, key)
      else
        valid_args = .false.
      end if
    end if
    if (valid_args) then
      call aot_top_get_val(val, ErrCode, L, default)
    else
      ErrCode = ibSet(0, aoterr_NonExistent)
      ErrCode = ibSet(ErrCode, aoterr_Fatal)
    end if

  end subroutine get_table_extdouble


  !> Put a extdouble precision real value into a table.
  subroutine set_table_extdouble(val, L, thandle, key, pos)
    type(flu_State) :: L !< Handle to the Lua script.

    !> Handle to the table to look the value up in.
    integer, intent(in) :: thandle

    !> Value of the table entry if it exists.
    real(kind=xdble_k), intent(in) :: val

    !> Name of the entry to look for.
    !!
    !! Key and pos are both optional, however at least one of them has to be
    !! supplied.
    !! The key takes precedence over the pos if both are given.
    character(len=*), intent(in), optional :: key

    !> Position of the entry to look for in the table.
    !!
    !! It allows the access to unnamed arrays in the Lua tables.
    integer, intent(in), optional :: pos

    real(kind=double_k) :: locval

    locval = real(val, kind=double_k)

    if (thandle > 0) then
      if (present(key)) then
        ! If there is a key, use that.
        ! First put the value on the top of the stack
        call flu_pushNumber(L, locval)
        ! Now put it into the table
        call flu_setField(L, thandle, trim(key))
      else
        ! No key given, try to put the value by position
        if (present(pos)) then
          ! First put the index, where to write the value into the table, on the
          ! stack.
          call flu_pushInteger(L, pos)
          ! Now put the actual value on the top of the stack.
          call flu_pushNumber(L, locval)
          ! Get the two entries from the stack into the table.
          call flu_setTable(L, thandle)
        end if
      end if
    end if

  end subroutine set_table_extdouble


  !> This subroutine takes a one dimensional array, and puts it as a table
  !! into the Lua context.
  !!
  !! The returned thandle provides the index to access this newly created
  !! table.
  subroutine create_1Darray_extdouble(L, thandle, val)
    type(flu_State) :: L !< Handle to the Lua script.

    !> Handle to access the newly created table.
    integer, intent(out) :: thandle

    !> Values to put into the new table.
    real(kind=xdble_k), intent(in) :: val(:)

    integer :: tab
    integer :: nvals
    integer :: i
    real(kind=double_k), allocatable :: locval(:)

    nVals = size(val)
    allocate(locVal(nVals))
    locVal(:) = real(val, kind=double_k)
    call flu_createtable(L, nVals, 0)
    thandle = flu_gettop(L)
    tab = thandle

    do i=1,nVals
      call flu_pushInteger(L, i)
      call flu_pushNumber(L, locval(i))
      call flu_settable(L, tab)
    end do

    deallocate(locval)

  end subroutine create_1Darray_extdouble


end module aot_extdouble_table_module
