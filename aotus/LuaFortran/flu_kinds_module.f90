! Copyright (c) 2016, 2018 Harald Klimach <harald@klimachs.de>
!
! Parts of this file were written by Harald Klimach for University of
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

!> Global definitions of some handy kind declarations
!! with the help of the intrinsic selected_*_kind
!! functions.
module flu_kinds_module
  implicit none

  integer, parameter :: quad_k = selected_real_kind(33)
  integer, parameter :: double_k = selected_real_kind(15)
  integer, parameter :: single_k = selected_real_kind(6)
  integer, parameter :: int_k = selected_int_kind(6)
  integer, parameter :: long_k = selected_int_kind(15)

end module flu_kinds_module
