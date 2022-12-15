! Copyright (c) 2016 Harald Klimach <harald@klimachs.de>
! Copyright (c) 2016 Peter Vitt <peter.vitt2@uni-siegen.de>
!
! Parts of this file were written by Harald Klimach and Peter Vitt for
! University of Siegen.
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

!> A small program to test the usage of vectors.
program aot_vec_test
  use flu_binding, only: flu_State

  use aotus_module, only: open_config_file, close_config, aot_get_val, &
    &                     aot_top_get_val
  use aot_fun_module, only: aot_fun_type, aot_fun_top, aot_fun_do, &
    &                       aot_fun_put, aot_fun_id, aot_fun_open, &
    &                       aot_fun_close
  use aot_references_module, only: aot_reference_for, aot_reference_to_top
  use aot_table_module, only: aot_table_open, aot_table_close

  implicit none

  type(flu_State) :: conf
  integer, allocatable :: fibo(:)
  integer, allocatable :: vError(:)
  character(len=22), allocatable :: names(:)
  real, allocatable :: res(:)

  integer :: iError
  integer :: thandle
  character(len=80) :: ErrString
  logical :: passed

  passed = .true.

  call create_script('aot_vec_test_config.lua')
  write(*,*)
  write(*,*) 'Running aot_vec_test...'
  write(*,*) ' * open_config_file (aot_vec_test_config.lua)'
  call open_config_file(L = conf, filename = 'aot_vec_test_config.lua', &
    &                   ErrCode = iError, ErrString = ErrString)
  if (iError /= 0) then
    write(*,*) ' : unexpected FATAL Error occured !!!'
    write(*,*) ' : Could not open the config file aot_vec_test_config.lua:'
    write(*,*) trim(ErrString)
    STOP
  end if
  write(*,*) '  : success.'

  ! Testing for references
  write(*,*) ' * Reading a global list of numbers'
  call aot_get_val( val       = fibo,   &
    &               ErrCode   = vError, &
    &               maxlength = 10,     &
    &               L         = conf,   &
    &               key       = 'fibo', &
    &               default   = [1, 1]  )

  if (size(fibo) /= 8) then
    write(*,*) 'ERROR: got wrong number of elements. Should be 8,'
    write(*,*) '       but is: ', size(fibo)
    passed = .false.
  end if

  if (any(vError /= 0)) then
    write(*,*) 'ERROR: Something went wrong while reading the global list of'
    write(*,*) '       numbers!'
    passed = .false.
  end if

  deallocate(fibo)
  deallocate(vError)

  write(*,*) ' * Reading a global scalar string into an allocatable array.'
  call aot_get_val( val       = names,  &
    &               ErrCode   = vError, &
    &               maxlength = 10,     &
    &               L         = conf,   &
    &               key       = 'aname' )

  if (size(names) /= 1) then
    write(*,*) 'ERROR: got wrong number of elements. Should be 1,'
    write(*,*) '       but is: ', size(names)
    passed = .false.
  end if

  if (any(vError /= 0)) then
    write(*,*) 'ERROR: Something went wrong while reading the global scalar'
    write(*,*) '       string!'
    passed = .false.
  end if

  deallocate(names)
  deallocate(vError)

  write(*,*) ' * Looking for vectors inside other tables.'
  call aot_table_open(L = conf, key = 'stuff', thandle = thandle)
  if (thandle > 0) then
    call aot_get_val( val       = names,      &
      &               ErrCode   = vError,     &
      &               maxlength = 10,         &
      &               L         = conf,       &
      &               thandle   = thandle,    &
      &               key       = 'namelist'  )

    if (size(names) /= 3) then
      write(*,*) 'ERROR: got wrong number of elements. Should be 3,'
      write(*,*) '       but is: ', size(names)
      passed = .false.
    end if

    if (any(vError /= 0)) then
      write(*,*) 'ERROR: Something went wrong while reading the list of strings'
      write(*,*) '       inside the table "stuff"!'
      passed = .false.
    end if

    deallocate(names)
    deallocate(vError)

    call aot_get_val( val       = res,        &
      &               ErrCode   = vError,     &
      &               maxlength = 10,         &
      &               L         = conf,       &
      &               thandle   = thandle,    &
      &               key       = 'anumber'   )

    if (size(res) /= 1) then
      write(*,*) 'ERROR: got wrong number of elements. Should be 1,'
      write(*,*) '       but is: ', size(names)
      passed = .false.
    end if

    if (any(vError /= 0)) then
      write(*,*) 'ERROR: Something went wrong while reading the scalar number'
      write(*,*) '       inside the table "stuff"!'
      passed = .false.
    end if

  else

    write(*,*) 'FATAL Error occured: could not open table "stuff"!'
    passed = .false.

  end if

  call aot_table_close(L = conf, thandle = thandle)

  write(*,*) ' * close_conf'
  call close_config(conf)
  write(*,*) '  : success.'
  write(*,*) '... Done with aot_vec_test.'
  if (passed) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

  ! for the sake of the address sanitizer
  deallocate(res, vError)

contains

  subroutine create_script(filename)
    character(len=*) :: filename

    open(file=trim(filename), unit=22, action='write', status='replace')
    write(22,*) '-- test script for aotus_vec_test'
    write(22,*) 'fibo = {1, 1, 2, 3, 5, 8, 13, 21}'
    write(22,*) 'aname = "hello"'
    write(22,*) 'stuff = {'
    write(22,*) '  namelist = {"Angelo", "Bertie", "Caesar"},'
    write(22,*) '  anumber = 42.1'
    write(22,*) '}'
    close(22)
  end subroutine create_script

end program aot_vec_test
