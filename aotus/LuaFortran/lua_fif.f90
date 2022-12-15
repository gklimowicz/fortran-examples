! Copyright (c) 2011-2012, 2015-2016 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2012-2013 James Spencer <j.spencer@imperial.ac.uk>
! Copyright (c) 2015 Peter Vitt <peter.vitt2@uni-siegen.de>
! Copyright (c) 2016 Ruth Franklin
!
! Parts of this file were written by Harald Klimach for
! German Research School of Simulation Sciences and University of
! Siegen.
! Parts of this file were written by Peter Vitt for University of
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

!> This module provides a direct translation of some
!! Lua 5.3.2
!! C-Interfaces to Fortran 2003 interfaces using the
!! ISO_C_BINDING facilities.
module lua_fif
  use, intrinsic :: iso_c_binding
  use lua_parameters

  implicit none

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
  ! Lua API interfaces
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
  interface

    subroutine lua_close(L) bind(c, name="lua_close")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
    end subroutine lua_close

    subroutine lua_createtable(L, narr, nrec) bind(c, name="lua_createtable")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: narr
      integer(kind=c_int), value :: nrec
    end subroutine lua_createtable

    function lua_getglobal(L, k) bind(c, name="lua_getglobal")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      character(kind=c_char), dimension(*) :: k
      integer(kind=c_int) :: lua_getglobal
    end function lua_getglobal

    function lua_getfield(L, index, k) bind(c, name="lua_getfield")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      character(kind=c_char), dimension(*) :: k
      integer(kind=c_int) :: lua_getfield
    end function lua_getfield

    function lua_gettable(L, index) bind(c, name="lua_gettable")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      integer(kind=c_int) :: lua_gettable
    end function lua_gettable

    function lua_gettop(L) bind(c, name="lua_gettop")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int) :: lua_gettop
    end function lua_gettop

    function lua_isNumber(L, index) bind(c, name="lua_isnumber")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      integer(kind=c_int) :: lua_isnumber
    end function lua_isnumber

    function lua_isString(L, index) bind(c, name="lua_isstring")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      integer(kind=c_int) :: lua_isString
    end function lua_isString

    function lua_next(L, index) bind(c, name="lua_next")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      integer(kind=c_int) :: lua_next
    end function lua_next

    function lua_pcallk(L, nargs, nresults, errfunc, ctx, k) bind(c, name="lua_pcallk")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: nargs
      integer(kind=c_int), value :: nresults
      integer(kind=c_int), value :: errfunc
      integer(kind=c_int), value :: ctx
      type(c_ptr), value :: k
      integer(kind=c_int) :: lua_pcallk
    end function lua_pcallk

    subroutine lua_pushinteger(L, n) bind(c, name="lua_pushinteger")
      use, intrinsic :: iso_c_binding
      use lua_parameters, only: lua_int
      type(c_ptr), value :: L
      integer(kind=lua_int), value :: n
    end subroutine lua_pushinteger

    subroutine lua_pushboolean(L, n) bind(c, name="lua_pushboolean")
      use, intrinsic :: iso_c_binding
      use lua_parameters, only: lua_int
      type(c_ptr), value :: L
      integer(kind=lua_int), value :: n
    end subroutine lua_pushboolean

    subroutine lua_pushnil(L) bind(c, name="lua_pushnil")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
    end subroutine lua_pushnil

    subroutine lua_pushnumber(L, n) bind(c, name="lua_pushnumber")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      real(kind=c_double), value :: n
    end subroutine lua_pushnumber

    function lua_pushlstring(L, s, len) bind(c, name="lua_pushlstring")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      character(kind=c_char), dimension(*) :: s
      integer(kind=c_size_t), value :: len
      type(c_ptr) :: lua_pushlstring
    end function lua_pushlstring

    subroutine lua_pushvalue(L, index) bind(c, name="lua_pushvalue")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
    end subroutine lua_pushvalue

    function lua_rawgeti(L, index, n) bind(c, name="lua_rawgeti")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      integer(kind=c_int), value :: n
      integer(kind=c_int) :: lua_rawgeti
    end function lua_rawgeti

    subroutine lua_rotate(L, idx, n) bind(c, name="lua_rotate")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: idx
      integer(kind=c_int), value :: n
    end subroutine lua_rotate

    subroutine lua_setfield(L, index, k) bind(c, name="lua_setfield")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      character(kind=c_char), dimension(*) :: k
    end subroutine lua_setfield

    subroutine lua_setglobal(L, k) bind(c, name="lua_setglobal")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      character(kind=c_char), dimension(*) :: k
    end subroutine lua_setglobal

    subroutine lua_settable(L, index) bind(c, name="lua_settable")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
    end subroutine lua_settable

    subroutine lua_settop(L, index) bind(c, name="lua_settop")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
    end subroutine lua_settop

    function lua_tolstring(L, index, len) bind(c, name="lua_tolstring")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      integer(kind=c_size_t) :: len
      type(c_ptr) :: lua_tolstring
    end function lua_tolstring

    function lua_tonumberx(L, index, isnum) bind(c, name="lua_tonumberx")
      use, intrinsic :: iso_c_binding
      use lua_parameters, only: lua_num
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      integer(kind=c_int) :: isnum
      real(kind=lua_num) :: lua_tonumberx
    end function lua_tonumberx

    function lua_toboolean(L, index) bind(c, name="lua_toboolean")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      integer(kind=c_int) :: lua_toboolean
    end function lua_toboolean

    function lua_touserdata(L, index) bind(c, name="lua_touserdata")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      type(c_ptr) :: lua_touserdata
    end function lua_touserdata

    function lua_topointer(L, index) bind(c, name="lua_topointer")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      integer(kind=c_intptr_t) :: lua_topointer
    end function lua_topointer

    function lua_type(L, index) bind(c, name="lua_type")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: index
      integer(kind=c_int) :: lua_type
    end function lua_type

    subroutine lua_pushcclosure(L, c_fn, n) bind(c, name="lua_pushcclosure")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      type(c_funptr), value :: c_fn
      integer(c_int), value :: n
    end subroutine lua_pushcclosure

    subroutine lua_pushlightuserdata(L, ptr) bind(c, name="lua_pushlightuserdata")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      type(c_ptr), value :: ptr
    end subroutine lua_pushlightuserdata

    function lua_getmetatable(L, index) bind(c, name="lua_getmetatable")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(c_int), value :: index
      integer(c_int) :: lua_getmetatable
    end function lua_getmetatable

  end interface
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !



  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
  ! Lua auxiliary library interfaces
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
  interface

    subroutine luaL_openlibs(L) bind(c, name="luaL_openlibs")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
    end subroutine luaL_openlibs

    function luaL_newstate() bind(c, name="luaL_newstate")
      use, intrinsic :: iso_c_binding
      type(c_ptr) :: luaL_newstate
    end function luaL_newstate

    function luaL_loadfilex(L, filename, mode) bind(c, name="luaL_loadfilex")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      character(kind=c_char), dimension(*) :: filename
      character(kind=c_char), dimension(*) :: mode
      integer(kind=c_int) :: luaL_loadfilex
    end function luaL_loadfilex

    function luaL_loadbufferx(L, buff, sz, name, mode) bind(c, name="luaL_loadbufferx")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      character(kind=c_char), dimension(*) :: buff
      integer(kind=c_size_t), value :: sz
      character(kind=c_char), dimension(*) :: name
      character(kind=c_char), dimension(*) :: mode
      integer(kind=c_int) :: luaL_loadbufferx
    end function luaL_loadbufferx

    function luaL_loadstring(L, string) bind(c, name="luaL_loadstring")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      character(kind=c_char), dimension(*) :: string
      integer(kind=c_int) :: luaL_loadstring
    end function luaL_loadstring

    subroutine luaL_setmetatable(L, tname) bind(c, name="luaL_setmetatable")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      character(kind=c_char), dimension(*) :: tname
    end subroutine luaL_setmetatable

    function luaL_newmetatable(L, tname) bind(c, name="luaL_newmetatable")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      character(kind=c_char), dimension(*) :: tname
      integer(kind=c_int) :: luaL_newmetatable
    end function luaL_newmetatable

    function luaL_ref(L, t) bind(c, name="luaL_ref")
      use, intrinsic :: iso_c_binding
      type(c_ptr), value :: L
      integer(kind=c_int), value :: t
      integer(kind=c_int) :: luaL_ref
    end function luaL_ref

  end interface
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !

end module lua_fif
