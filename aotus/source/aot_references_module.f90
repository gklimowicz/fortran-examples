! Copyright (c) 2016 Harald Klimach <harald@klimachs.de>
!
! Parts of this file were written by Harald Klimach for University of Siegen.
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

!> This module provides procedures to deal with Lua references.
!!
!! See for example Programming in Lua
!! [27.3.2](http://www.lua.org/pil/27.3.2.html).
!!
module aot_references_module
  use flu_binding
  use lua_parameters, only: LUA_REGISTRYINDEX
  use aot_table_ops_module, only: aot_table_push

  implicit none

  private

  public :: aot_reference_for, aot_reference_to_top


contains


  !> Get a reference for the entry defined by thandle, key and pos, or
  !! the current top entry in the stack.
  !!
  !! The reference can be used to refer to a given object in the Lua
  !! table by storing a reference to it in the LUA_REGISTRYINDEX table.
  !!
  !! The object can then be put onto the stack again by using this reference.
  function aot_reference_for(L, thandle, key, pos) result(ref)
    type(flu_State) :: L !! Handle to the Lua script

    !> Handle to the table containing the object to get a reference for.
    integer, intent(in), optional :: thandle

    !> Name of the object to look up, if thandle is not present, this is
    !! a global definition.
    !!
    !! If neither thandle nor key is provided, a reference to the top of
    !! the stack is returned.
    character(len=*), intent(in), optional :: key

    !> Positional index of the object inside thandle to get the reference
    !! for. If thandle is not provided, this argument is ignored.
    !!
    !! If a key is provided, that takes precedent over pos.
    integer, intent(in), optional :: pos

    integer :: toptype
    integer :: ref

    if (present(thandle)) then
      call aot_table_push(L=L, thandle=thandle, &
        &                 key=key, pos=pos      )
    else if (present(key)) then
      toptype = flu_getglobal(L, key)
    end if

    ref = fluL_ref(L, LUA_REGISTRYINDEX)
  end function aot_reference_for


  !> Put a given reference (ref) in the Lua script (L) to the top of the stack.
  subroutine aot_reference_to_top(L, ref)
    type(flu_State) :: L !! Handle to the Lua script

    !> Reference retrieved by [[aot_reference_for]] to put on the top
    !! of the stack
    integer :: ref

    integer :: luatype

    luatype = flu_rawgeti(L, LUA_REGISTRYINDEX, ref)
  end subroutine aot_reference_to_top

end module aot_references_module
