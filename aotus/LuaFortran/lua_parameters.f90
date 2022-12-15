! Copyright (c) 2011-2012, 2016 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2012 Kannan Masilamani <k.masilamani@grs-sim.de>
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

!> This module provides some parameters defined in the
!! Lua header file that are needed in the wrapper
!! functions for the Fortran interface.
!!
!! Lua_int and Lua_num are system dependent, and you
!! might need to adapt them on your system.
!! The type constants have to be consistent with the
!! Lua header definition, and thus should be checked
!! after version upgrades of the Lua library.
module lua_parameters
  use, intrinsic :: iso_c_binding

  implicit none

  ! System dependent, might need to be adapted:
  integer, parameter :: lua_int = c_long
  integer, parameter :: lua_num = c_double

  ! Lua config constants (see luaconf.h)
  ! Attention: might need to be adapted!
  integer(kind=c_int), parameter :: LUAI_MAXSTACK = 1000000

  ! Lua constants (see lua.h)
  integer(kind=c_int), parameter :: LUA_TNONE          = -1
  integer(kind=c_int), parameter :: LUA_TNIL           =  0
  integer(kind=c_int), parameter :: LUA_TBOOLEAN       =  1
  integer(kind=c_int), parameter :: LUA_TLIGHTUSERDATA =  2
  integer(kind=c_int), parameter :: LUA_TNUMBER        =  3
  integer(kind=c_int), parameter :: LUA_TSTRING        =  4
  integer(kind=c_int), parameter :: LUA_TTABLE         =  5
  integer(kind=c_int), parameter :: LUA_TFUNCTION      =  6
  integer(kind=c_int), parameter :: LUA_TUSERDATA      =  7
  integer(kind=c_int), parameter :: LUA_TTHREAD        =  8

  integer(kind=c_int), parameter :: LUA_REGISTRYINDEX  = -LUAI_MAXSTACK - 1000

end module lua_parameters
