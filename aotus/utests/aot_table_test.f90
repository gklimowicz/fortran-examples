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

program aot_table_test
  use flu_binding, only: flu_State

  use aotus_module, only: open_config_file, close_config
  use aot_top_module, only: aoterr_Fatal, aoterr_NonExistent, aoterr_WrongType
  use aot_table_module, only: aot_table_open, aot_table_close, &
    &                         aot_table_length, aot_table_get_val

  implicit none

  integer, parameter :: primes(5) = [ 1, 2, 3, 5, 7 ]

  type(flu_State) :: conf
  integer :: globhandle
  integer :: tabint
  integer :: tablen
  integer :: iError
  integer :: i
  character(len=80) :: ErrString
  logical :: passed

  call create_script('aot_table_test_config.lua')

  write(*,*)
  write(*,*) 'Running aot_table_test...'

  passed = .true.
  call open_config_file(L = conf, filename = 'aot_table_test_config.lua', &
    &                   ErrCode = iError, ErrString = ErrString)
  if (iError /= 0) then
    write(*,*) 'Unexpected FATAL Error occured !!!'
    write(*,*) 'Could not open the config file aot_table_test_config.lua:'
    write(*,*) trim(ErrString)
    STOP
  end if

  ! Testing for global Table
  write(*,*) ' * opening a global table'
  call aot_table_open(L = conf, thandle = globhandle, key = 'primes')
  if (globhandle == 0) then
    write(*,*) '  : unexpected FATAL Error occured !!!'
    write(*,*) '  : could not open global table primes.'
    passed = .false.
  else
    write(*,*) '  : success.'
    write(*,*) ' * getting the length of the table'
    tablen = aot_table_length(L = conf, thandle = globhandle)
    if (tablen /= 5) then
      write(*,*) '  : unexpected FATAL Error occured !!!'
      write(*,*) '  : found a table length of ', tablen
      write(*,*) '  :    but should have been ', 5
      passed = .false.
    else
      write(*,*) '  : success.'
      write(*,*) ' * retrieving entries of table by position'
      do i=1,5
        call aot_table_get_val(L = conf, thandle = globhandle, &
          &                pos = i, &
          &                val = tabint, ErrCode = iError)
        if (btest(iError, aoterr_Fatal)) then
          write(*,*) '  : unexpected FATAL Error occured !!!'
          passed = .false.
          if (btest(iError, aoterr_NonExistent)) &
            &   write(*,*) '  : Variable not existent!'
          if (btest(iError, aoterr_WrongType)) &
            &   write(*,*) '  : Variable has wrong type!'
          exit
        else
          if (tabint /= primes(i)) then
            write(*,*) '  : unexpected ERROR, value mismatch, got: ', tabint
            write(*,*) '  :                             should be: ', primes(i)
            passed = .false.
            iError = 42
            exit
          end if
        end if
      end do
      if (iError == 0) write(*,*) '  : success.'
      write(*,*) ' * Attempting read of nonexistent entry'
      call aot_table_get_val(L = conf, thandle = globhandle, &
        &                pos = 10, &
        &                val = tabint, ErrCode = iError)
      if (btest(iError, aoterr_Fatal)) then
        write(*,*) '  : success.'
      else
        passed = .false.
        write(*,*) '  : ERROR, unexpected success in reading nonexistent entry'
        if (btest(iError, aoterr_NonExistent)) &
          &   write(*,*) '  : Variable not existent, but should be fatal!'
        if (btest(iError, aoterr_WrongType)) &
          &   write(*,*) '  : Variable has wrong type but should not exist!'
      end if
    end if
  end if

  write(*,*) ' * Closing table'
  call aot_table_close(L = conf, thandle = globhandle)
  write(*,*) '  : success.'

  call close_config(conf)
  write(*,*) '... Done with aot_table_test.'

  if (passed) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

contains

  subroutine create_script(filename)
    character(len=*) :: filename

    open(file=trim(filename), unit=22, action='write', status='replace')
    write(22,*) '-- test script for aot_table_test'
    write(22,*) 'primes = { 1, 2, 3, 5, 7 }'
    write(22,*) 'tabtab = { origin = {0.0, 0.0, 0.0} }'
    write(22,*) 'unnamed = { { kind = 1 } }'
    close(22)
  end subroutine create_script

end program aot_table_test
