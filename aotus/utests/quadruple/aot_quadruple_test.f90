! Copyright (c) 2012-2013 Harald Klimach <harald@klimachs.de>
!
! Parts of this file were written by Harald Klimach for German Research School
! of Simulation Sciences.
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

program quadruple_test
  use flu_binding, only: flu_State

  use aotus_module, only: open_config_file, close_config, aot_get_val
  use aot_top_module, only: aoterr_Fatal, aoterr_NonExistent, aoterr_WrongType
  use aot_table_module, only: aot_table_open, aot_table_close

  implicit none

  integer, parameter :: quad_k = selected_real_kind(33)

  type(flu_State) :: conf
  integer :: iError
  real(kind=quad_k) :: glob_quad
  real(kind=quad_k) :: tab_quad
  integer :: thandle
  character(len=80) :: ErrString
  logical :: passed

  passed = .true.

  call create_script('quad_test_config.lua')
  write(*,*)
  write(*,*) 'Running aotus_test...'
  write(*,*) ' * open_config_file (aotus_test_config.lua)'
  call open_config_file(L = conf, filename = 'quad_test_config.lua', &
    &                   ErrCode = iError, ErrString = ErrString)
  if (iError /= 0) then
    write(*,*) ' : unexpected FATAL Error occured !!!'
    write(*,*) ' : Could not open the config file quad_test_config.lua:'
    write(*,*) trim(ErrString)
    STOP
  end if
  write(*,*) '  : success.'

  ! Testing for global INTEGER
  write(*,*) ' * reading a global quadruple'
  call aot_get_val(L = conf, key = 'real_test', &
    &              val = glob_quad, ErrCode = iError)

  if (btest(iError, aoterr_Fatal)) then
    write(*,*) '  : unexpected FATAL Error occured !!!'
    if (btest(iError, aoterr_NonExistent)) &
      &   write(*,*) '  : Variable not existent!'
    if (btest(iError, aoterr_WrongType)) &
      &   write(*,*) '  : Variable has wrong type!'
    passed = .false.
  else
    if ((glob_quad > 0.5_quad_k*(1.0_quad_k-epsilon(glob_quad))) &
      & .and. (glob_quad < 0.5_quad_k*(1.0_quad_k+epsilon(glob_quad)))) then
      write(*,*) '  : success.'
    else
      write(*,*) '  : unexpected ERROR, value mismatch, got: ', glob_quad
      write(*,*) '  :                             should be: ', 0.5
      passed = .false.
    end if
  end if

  ! Testing for global Table
  write(*,*) ' * opening a global table'
  call aot_table_open(L = conf, thandle = thandle, key = 'tab')
  if (thandle == 0) then
    write(*,*) '  : unexpected FATAL Error occured !!!'
    write(*,*) '  : could not open global table primes.'
    passed = .false.
  end if

  call aot_get_val(L = conf, thandle = thandle, &
    &              key = 'real_in_tab', &
    &              val = tab_quad, ErrCode = iError)
  if (btest(iError, aoterr_Fatal)) then
    write(*,*) '  : unexpected FATAL Error occured !!!'
    passed = .false.
    if (btest(iError, aoterr_NonExistent)) &
      &   write(*,*) '  : Variable not existent!'
    if (btest(iError, aoterr_WrongType)) &
      &   write(*,*) '  : Variable has wrong type!'
  else
    if ((tab_quad < 2.0_quad_k*(1._quad_k-epsilon(tab_quad))) &
      & .or. (tab_quad > 2.0_quad_k*(1._quad_k+epsilon(tab_quad)))) then
      write(*,*) '  : unexpected ERROR, value mismatch, got: ', tab_quad
      write(*,*) '  :                             should be: ', 2.0
      passed = .false.
      iError = 42
    end if
  end if
  if (iError == 0) write(*,*) '  : success.'

  write(*,*) ' * Closing table'
  call aot_table_close(L = conf, thandle = thandle)
  write(*,*) '  : success.'

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
    write(22,*) 'tab = {real_in_tab = 2.0}'
    write(22,*) 'log_test = true'
    write(22,*) "string_test = 'last words'"
    close(22)
  end subroutine create_script

end program quadruple_test
