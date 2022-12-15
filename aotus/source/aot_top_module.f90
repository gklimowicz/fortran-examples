! Copyright (c) 2012-2016, 2018 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2013 James Spencer <j.spencer@imperial.ac.uk>
!
! Parts of this file were written by Harald Klimach for
! German Research School of Simulation Sciences and University of
! Siegen.
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

!> Module for interaction with topmost element of the Lua stack.
!!
!! This is a basic module which provides the fundamental functionality to
!! access the topmost element in the stack of the Lua API.
!! All intrinsic variables except complex numbers can be accessed this way.
module aot_top_module
  use flu_binding
  use flu_kinds_module, only: double_k, single_k, int_k, long_k
  use aot_err_module, only: aoterr_Fatal, aoterr_NonExistent, &
    &                       aoterr_WrongType, aot_err_handler

  ! The following module enables an interface for quadruple precision numbers,
  ! if the compiler supports them. However, you should be aware, that this is
  ! merely a convenience interface, as the values provided by Lua are only
  ! double precision.
  use aot_quadruple_top_module

  ! Support for extended double precision.
  use aot_extdouble_top_module

  implicit none

  private

  public :: aoterr_Fatal, aoterr_NonExistent, aoterr_WrongType
  public :: aot_top_get_val
  public :: aot_err_handler

  !> Get the value on top of the Lua API stack
  !!
  !! This is the most basic operation to retrieve a value.
  !! It is also most flexible in the sense, that it does not matter how the
  !! value actually gets on top of the stack by previous Lua operations.
  !!
  !! The interface looks like this:
  !! `call aot_top_get_val(val, errCode, L, default)`.
  !! See for example [[aot_top_get_real]] for a more detailed description of the
  !! parameters.
  !!
  !! aot_top_get_val can not be in the same generic interface as the other
  !! [[aot_get_val]] routines, as it results in ambiguities of the interfaces.
  !!
  !! @note The retrieved value will be popped from the Lua API stack.
  interface aot_top_get_val
    module procedure aot_top_get_real
    module procedure aot_top_get_double
    module procedure aot_top_get_integer
    module procedure aot_top_get_long
    module procedure aot_top_get_string
    module procedure aot_top_get_logical
    module procedure aot_top_get_userdata
  end interface

contains

  !> Interpret topmost entry on Lua stack as a single precision real.
  subroutine aot_top_get_real(val, ErrCode, L, default)
    type(flu_State) :: L !! Handle to the Lua script

    !> Value of the Variable in the script
    real(kind=single_k), intent(out) :: val

    !> Error code to indicate what kind of problem might have occured.
    integer, intent(out) :: ErrCode

    !> Some default value, that should be used, if the variable is not set in
    !! the Lua script.
    real(kind=single_k), optional, intent(in) :: default

    logical :: not_retrievable

    ErrCode = 0
    not_retrievable = .false.

    if (flu_isNoneOrNil(L, -1)) then
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      not_retrievable = .true.
    else
      if (flu_isNumber(L, -1)) then
        val = flu_toNumber(L, -1)
      else
        ErrCode = ibSet(ErrCode, aoterr_WrongType)
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
        not_retrievable = .true.
      end if
    end if

    if (not_retrievable) then
      if (present(default)) then
        val = default
      else
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call flu_pop(L)

  end subroutine aot_top_get_real


  !> Interpret topmost entry on Lua stack as a double precision real.
  subroutine aot_top_get_double(val, ErrCode, L, default)
    type(flu_State) :: L !! Handle to the Lua script

    !> Value of the Variable in the script
    real(kind=double_k), intent(out) :: val

    !> Error code to indicate what kind of problem might have occured.
    integer, intent(out) :: ErrCode

    !> Some default value, that should be used, if the variable is not set in
    !! the Lua script.
    real(kind=double_k), optional, intent(in) :: default

    logical :: not_retrievable

    ErrCode = 0
    not_retrievable = .false.

    if (flu_isNoneOrNil(L, -1)) then
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      not_retrievable = .true.
    else
      if (flu_isNumber(L, -1)) then
        val = flu_toDouble(L, -1)
      else
        ErrCode = ibSet(ErrCode, aoterr_WrongType)
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
        not_retrievable = .true.
      end if
    end if

    if (not_retrievable) then
      if (present(default)) then
        val = default
      else
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call flu_pop(L)

  end subroutine aot_top_get_double


  !> Interpret topmost entry on Lua stack as a default integer.
  subroutine aot_top_get_integer(val, ErrCode, L, default)
    type(flu_State) :: L !! Handle to the Lua script

    !> Value of the Variable in the script
    integer(kind=int_k), intent(out) :: val

    !> Error code to indicate what kind of problem might have occured.
    integer, intent(out) :: ErrCode

    !> Some default value, that should be used, if the variable is not set in
    !! the Lua script.
    integer(kind=int_k), optional, intent(in) :: default

    logical :: not_retrievable

    ErrCode = 0
    not_retrievable = .false.

    if (flu_isNoneOrNil(L, -1)) then
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      not_retrievable = .true.
    else
      if (flu_isNumber(L, -1)) then
        val = int(flu_toDouble(L, -1))
      else
        ErrCode = ibSet(ErrCode, aoterr_WrongType)
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
        not_retrievable = .true.
      end if
    end if

    if (not_retrievable) then
      if (present(default)) then
        val = default
      else
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call flu_pop(L)

  end subroutine aot_top_get_integer


  !> Interpret topmost entry on Lua stack as a single precision real.
  subroutine aot_top_get_long(val, ErrCode, L, default)
    type(flu_State) :: L !! Handle to the Lua script

    !> Value of the Variable in the script
    integer(kind=long_k), intent(out) :: val

    !> Error code to indicate what kind of problem might have occured.
    integer, intent(out) :: ErrCode

    !> Some default value, that should be used, if the variable is not set in
    !! the Lua script.
    integer(kind=long_k), optional, intent(in) :: default

    logical :: not_retrievable

    ErrCode = 0
    not_retrievable = .false.

    if (flu_isNoneOrNil(L, -1)) then
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      not_retrievable = .true.
    else
      if (flu_isNumber(L, -1)) then
        val = int(flu_toDouble(L, -1),kind=long_k)
      else
        ErrCode = ibSet(ErrCode, aoterr_WrongType)
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
        not_retrievable = .true.
      end if
    end if

    if (not_retrievable) then
      if (present(default)) then
        val = default
      else
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call flu_pop(L)

  end subroutine aot_top_get_long


  !> Interpret topmost entry on Lua stack as a single precision real.
  subroutine aot_top_get_logical(val, ErrCode, L, default)
    type(flu_State) :: L !! Handle to the Lua script

    !> Value of the Variable in the script
    logical, intent(out) :: val

    !> Error code to indicate what kind of problem might have occured.
    integer, intent(out) :: ErrCode

    !> Some default value, that should be used, if the variable is not set in
    !! the Lua script.
    logical, optional, intent(in) :: default

    logical :: not_retrievable

    ErrCode = 0
    not_retrievable = .false.

    if (flu_isNoneOrNil(L, -1)) then
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      not_retrievable = .true.
    else
      if (flu_isBoolean(L, -1)) then
        val = flu_toBoolean(L, -1)
      else
        ErrCode = ibSet(ErrCode, aoterr_WrongType)
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
        not_retrievable = .true.
      end if
    end if

    if (not_retrievable) then
      if (present(default)) then
        val = default
      else
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call flu_pop(L)

  end subroutine aot_top_get_logical


  !> Interpret topmost entry on Lua stack as a single precision real.
  subroutine aot_top_get_string(val, ErrCode, L, default)
    type(flu_State) :: L !! Handle to the Lua script

    !> Value of the Variable in the script
    character(len=*) :: val

    !> Error code to indicate what kind of problem might have occured.
    integer, intent(out) :: ErrCode

    !> Some default value, that should be used, if the variable is not set in
    !! the Lua script.
    character(len=*), optional, intent(in) :: default

    logical :: not_retrievable
    character, pointer :: cstring(:)
    integer :: i, StrLen, StrLimit

    ErrCode = 0
    not_retrievable = .false.

    if (flu_isNoneOrNil(L, -1)) then
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      not_retrievable = .true.
    else
      if (flu_isString(L, -1)) then
        cstring => flu_toLString(L, -1, StrLen)
        StrLimit = min(StrLen, len(val))
        val = ''
        do i=1,StrLimit
          val(i:i) = cstring(i)
        end do
      else
        ErrCode = ibSet(ErrCode, aoterr_WrongType)
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
        not_retrievable = .true.
      end if
    end if

    if (not_retrievable) then
      if (present(default)) then
        val = default
      else
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call flu_pop(L)

  end subroutine aot_top_get_string


  !> Interpret topmost entry on Lua stack as userdata.
  subroutine aot_top_get_userdata(val, ErrCode, L, default)
    use, intrinsic :: iso_c_binding
    type(flu_State) :: L !! Handle to the Lua script

    !> Value of the Variable in the script
    type(c_ptr), intent(out) :: val

    !> Error code to indicate what kind of problem might have occured.
    integer, intent(out) :: ErrCode

    !> Some default value, that should be used, if the variable is not set in
    !! the Lua script.
    type(c_ptr), optional, intent(in) :: default

    logical :: not_retrievable

    ErrCode = 0
    not_retrievable = .false.

    if (flu_isNoneOrNil(L, -1)) then
      ErrCode = ibSet(ErrCode, aoterr_NonExistent)
      not_retrievable = .true.
    else
      if (flu_islightuserdata(L, -1)) then
        val = flu_touserdata(L, -1)
      else
        ErrCode = ibSet(ErrCode, aoterr_WrongType)
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
        not_retrievable = .true.
      end if
    end if

    if (not_retrievable) then
      if (present(default)) then
        val = default
      else
        ErrCode = ibSet(ErrCode, aoterr_Fatal)
      end if
    end if
    call flu_pop(L)

  end subroutine aot_top_get_userdata


end module aot_top_module
