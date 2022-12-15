! Copyright (c) 2015 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Parts of this file were written by Peter Vitt for University of Siegen.
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

program aot_path_test

  use aot_path_module,  only: aot_path_type,     &
    &                         aot_init_path,     &
    &                         aot_path_addNode,  &
    &                         aot_path_delNode,  &
    &                         aot_path_toString, &
    &                         aot_path_dump

  implicit none

  type(aot_path_type) :: path
  character(len=256) :: buffer
  character(len=10) :: shortBuffer
  logical :: passed = .true.

  call aot_init_path(path)

  call aot_path_addNode( me = path, NodeType = 'table', key = 'level1' )
  call aot_path_addNode( me = path, NodeType = 'table', key = 'level2' )
  call aot_path_addNode( me = path, NodeType = 'table', key = 'level3' )

  ! Test whether the conversion works at all
  call aot_path_toString( path, buffer )
  if (buffer /= 'level1.level2.level3') then
    write(*,*) buffer
    write(*,*) 'Converting path into string failed'
    passed = .false.
  end if

  call aot_path_delNode( path )

  call aot_path_toString( path, buffer )
  if (buffer /= 'level1.level2') then
    write(*,*) 'Converting path into string failed after deleting last node'
    passed = .false.
  end if

  ! Does it return an empty buffer when the buffer is to short?
  shortBuffer = 'test'
  call aot_path_toString( path, shortBuffer )
  if (shortBuffer /= '') then
    write(*,*) 'Buffer not empty when result is too long.'
    passed = .false.
  end if

  if (passed) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

end program aot_path_test
