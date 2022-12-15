! Copyright (c) 2012-2013, 2016 Harald Klimach <harald@klimachs.de>
!
! Parts of this file were written by Harald Klimach for German Research School
! of Simulation Sciences and University of Siegen.
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

program aotus_test
  use flu_binding, only: flu_State

  use flu_kinds_module, only: long_k
  use aotus_module, only: open_config_file, close_config, aot_get_val
  use aot_top_module, only: aoterr_Fatal, aoterr_NonExistent, aoterr_WrongType

  implicit none

  type(flu_State) :: conf
  integer :: iError
  integer :: glob_int
  integer(kind=long_k) :: glob_long
  real :: glob_real
  logical :: glob_log
  character(len=20) :: glob_string
  character(len=80) :: ErrString
  logical :: passed

  passed = .true.

  call create_script('aotus_test_config.lua')
  write(*,*)
  write(*,*) 'Running aotus_test...'
  write(*,*) ' * open_config_file (aotus_test_config.lua)'
  call open_config_file(L = conf, filename = 'aotus_test_config.lua', &
    &                   ErrCode = iError, ErrString = ErrString)
  if (iError /= 0) then
    write(*,*) ' : unexpected FATAL Error occured !!!'
    write(*,*) ' : Could not open the config file aotus_test_config.lua:'
    write(*,*) trim(ErrString)
    STOP
  end if
  write(*,*) '  : success.'

  ! Testing for global INTEGER
  write(*,*) ' * reading a global integer'
  call aot_get_val(L = conf, key = 'int_test', &
    &                 val = glob_int, ErrCode = iError)

  if (btest(iError, aoterr_Fatal)) then
    write(*,*) '  : unexpected FATAL Error occured !!!'
    if (btest(iError, aoterr_NonExistent)) &
      &   write(*,*) '  : Variable not existent!'
    if (btest(iError, aoterr_WrongType)) &
      &   write(*,*) '  : Variable has wrong type!'
    passed = .false.
  else
    if (glob_int == 5) then
      write(*,*) '  : success.'
    else
      write(*,*) '  : unexpected ERROR, value mismatch, got: ', glob_int
      write(*,*) '  :                             should be: ', 5
      passed = .false.
    end if
  end if
  ! -------------------------------- !


  ! Testing for global LONG
  write(*,*) ' * reading a global long'
  call aot_get_val(L = conf, key = 'long_test', &
    &                 val = glob_long, ErrCode = iError)

  if (btest(iError, aoterr_Fatal)) then
    write(*,*) '  : unexpected FATAL Error occured !!!'
    if (btest(iError, aoterr_NonExistent)) &
      &   write(*,*) '  : Variable not existent!'
    if (btest(iError, aoterr_WrongType)) &
      &   write(*,*) '  : Variable has wrong type!'
    passed = .false.
  else
    if (glob_long == 5000000000_long_k) then
      write(*,*) '  : success.'
    else
      write(*,*) '  : unexpected ERROR, value mismatch, got: ', glob_long
      write(*,*) '  :                             should be: ', &
        &                                         5000000000_long_k
      passed = .false.
    end if
  end if
  ! -------------------------------- !

  ! Testing for global REAL
  write(*,*) ' * reading a global real'
  call aot_get_val(L = conf, key = 'real_test', &
    &              val = glob_real, ErrCode = iError)

  if (btest(iError, aoterr_Fatal)) then
    write(*,*) '  : unexpected FATAL Error occured !!!'
    if (btest(iError, aoterr_NonExistent)) &
      &   write(*,*) '  : Variable not existent!'
    if (btest(iError, aoterr_WrongType)) &
      &   write(*,*) '  : Variable has wrong type!'
    passed = .false.
  else
    if ((glob_real > 0.5*(1.0-epsilon(glob_real))) &
      & .or. (glob_real < 0.5*(1.0+epsilon(glob_real)))) then
      write(*,*) '  : success.'
    else
      write(*,*) '  : unexpected ERROR, value mismatch, got: ', glob_real
      write(*,*) '  :                             should be: ', 0.5
      passed = .false.
    end if
  end if
  ! -------------------------------- !

  ! Testing for global LOGICAL
  glob_log = .false.
  write(*,*) ' * reading a global logical'
  call aot_get_val(L = conf, key = 'log_test', &
    &              val = glob_log, ErrCode = iError)

  if (btest(iError, aoterr_Fatal)) then
    write(*,*) '  : unexpected FATAL Error occured !!!'
    if (btest(iError, aoterr_NonExistent)) &
      &   write(*,*) '  : Variable not existent!'
    if (btest(iError, aoterr_WrongType)) &
      &   write(*,*) '  : Variable has wrong type!'
    passed = .false.
  else
    if (glob_log) then
      write(*,*) '  : success.'
    else
      write(*,*) '  : unexpected ERROR, value mismatch, got: ', glob_log
      write(*,*) '  :                             should be: ', .true.
      passed = .false.
    end if
  end if
  ! -------------------------------- !

  ! Testing for global STRING
  write(*,*) ' * reading a global string'
  call aot_get_val(L = conf, key = 'string_test', &
    &              val = glob_string, ErrCode = iError)

  if (btest(iError, aoterr_Fatal)) then
    write(*,*) '  : unexpected FATAL Error occured !!!'
    if (btest(iError, aoterr_NonExistent)) &
      &   write(*,*) '  : Variable not existent!'
    if (btest(iError, aoterr_WrongType)) &
      &   write(*,*) '  : Variable has wrong type!'
    passed = .false.
  else
    if (trim(glob_string) == 'last words') then
      write(*,*) '  : success.'
    else
      write(*,*) '  : unexpected ERROR, value mismatch, got: ', glob_string
      write(*,*) '  :                             should be: ', 'last words'
      passed = .false.
    end if
  end if
  ! -------------------------------- !

  ! Testing for global non-fatal NONEXISTENT
  write(*,*) ' * reading a nonexistent global integer with default'
  call aot_get_val(L = conf, key = 'nonexist', &
    &                 val = glob_int, ErrCode = iError, &
    &                 default = 1)

  if (btest(iError, aoterr_Fatal)) then
    write(*,*) '  : unexpected FATAL Error occured !!!'
    passed = .false.
  else
    if (btest(iError, aoterr_WrongType)) then
      write(*,*) '  : unexpected ERROR, found WrongType !!!'
      passed = .false.
    else
      if (btest(iError, aoterr_NonExistent)) then
        if (glob_int == 1) then
          write(*,*) '  : success.'
        else
          write(*,*) '  : unexpected ERROR, value mismatch, got: ', glob_int
          write(*,*) '  :                             should be: ', 1
          passed = .false.
        end if
      else
        write(*,*) '  : ERROR: missing aoterr_NonExistent !!!'
        write(*,*) '  :        should receive error when looking up "nonexist"'
        passed = .false.
      end if
    end if
  end if
  ! -------------------------------- !

  ! Testing for global fatal NONEXISTENT
  write(*,*) ' * reading a nonexistent global integer NO default'
  call aot_get_val(L = conf, key = 'nonexist', &
    &                 val = glob_int, ErrCode = iError)

  if (btest(iError, aoterr_Fatal)) then
    if (btest(iError, aoterr_NonExistent)) then
      write(*,*) '  : success.'
    else
      write(*,*) '  : ERROR no aoterr_NonExistent returned !!!'
      passed = .false.
    end if
  else
    write(*,*) '  : ERROR unexpectly no fatal Error code returned !!!'
    passed = .false.
  end if
  ! -------------------------------- !

  ! Testing for global wrong type
  write(*,*) ' * reading a global integer from variable with wrong type'
  call aot_get_val(L = conf, key = 'string_test', &
    &                 val = glob_int, ErrCode = iError)

  if (btest(iError, aoterr_Fatal)) then
    if (btest(iError, aoterr_WrongType)) then
      write(*,*) '  : success.'
    else
      write(*,*) '  : ERROR no aoterr_WrongType returned !!!'
      passed = .false.
    end if
  else
    write(*,*) '  : ERROR unexpectly no fatal Error code returned !!!'
    passed = .false.
  end if
  ! -------------------------------- !


  write(*,*) ' * close_conf'
  call close_config(conf)
  write(*,*) '  : success.'
  write(*,*) '... Done with aotus_test.'
  if (passed) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

contains

  subroutine create_script(filename)
    character(len=*) :: filename

    open(file=trim(filename), unit=22, action='write', status='replace')
    write(22,*) '-- test script for aotus_test'
    write(22,*) 'int_test = 5'
    write(22,*) 'long_test = 5000000000'
    write(22,*) 'real_test = 0.5'
    write(22,*) 'log_test = true'
    write(22,*) "string_test = 'last words'"
    close(22)
  end subroutine create_script

end program aotus_test
