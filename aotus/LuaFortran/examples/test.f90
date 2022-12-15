! Copyright (c) 2011, 2016 Harald Klimach <harald@klimachs.de>
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

program test
  use flu_binding

  implicit none

  type(flu_State) :: L
  character(len=89) :: filename
  integer :: str_len
  integer :: widthtype, heighttype

  L = fluL_newstate()

  filename = "config.test"
  if (fluL_loadfile(L, filename).ne.0) then
     write(*,*) "cannot load configuration file: ", flu_tolstring(L,-1, str_len)
     STOP
  endif
  if (flu_pcall(L, 0, 0, 0).ne.0) then
     write(*,*) "cannot run configuration file: ", flu_tolstring(L,-1, str_len)
     STOP
  endif

  widthtype = flu_getglobal(L, "width")
  heighttype = flu_getglobal(L, "height")

  if (widthtype == FLU_TNUMBER) then
    write(*,*) 'width=', flu_tonumber(L, -2)
  end if
  if (heighttype == FLU_TNUMBER) then
    write(*,*) 'height=', flu_tonumber(L, -1)
  end if

  call flu_close(L)

end program test
