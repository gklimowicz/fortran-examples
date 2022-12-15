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

!> A small program to test the usage of Lua references for functions.
program aot_ref_test
  use flu_binding, only: flu_State, flu_topointer, flu_pop

  use aotus_module, only: open_config_file, close_config, aot_get_val, &
    &                     aot_top_get_val
  use aot_fun_module, only: aot_fun_type, aot_fun_top, aot_fun_do, &
    &                       aot_fun_put, aot_fun_id, aot_fun_open, &
    &                       aot_fun_close
  use aot_references_module, only: aot_reference_for, aot_reference_to_top
  use aot_err_module, only: aoterr_Fatal, aoterr_NonExistent, aoterr_WrongType
  use aot_table_module, only: aot_table_open, aot_table_close

  implicit none

  type(flu_State) :: conf
  type(aot_fun_type) :: luafun
  integer :: iError
  integer :: globref, tabref
  integer :: thandle
  real :: res
  character(len=80) :: ErrString
  logical :: passed

  passed = .true.

  call create_script('aot_ref_test_config.lua')
  write(*,*)
  write(*,*) 'Running aot_ref_test...'
  write(*,*) ' * open_config_file (aot_ref_test_config.lua)'
  call open_config_file(L = conf, filename = 'aot_ref_test_config.lua', &
    &                   ErrCode = iError, ErrString = ErrString)
  if (iError /= 0) then
    write(*,*) ' : unexpected FATAL Error occured !!!'
    write(*,*) ' : Could not open the config file aot_ref_test_config.lua:'
    write(*,*) trim(ErrString)
    STOP
  end if
  write(*,*) '  : success.'

  ! Testing for references
  write(*,*) ' * reference from the global variable with the function'
  globref = aot_reference_for(L = conf, key = 'funref')
  write(*,*) '   : got the reference: ', globref
  call aot_reference_to_top(L = conf, ref = globref)
  write(*,*) '   : pointer: ', flu_topointer(conf, -1)
  call flu_pop(conf)

  write(*,*) ' * reference from inside a table'
  call aot_table_open(L = conf, key = 'test_table', thandle = thandle)

  if (thandle > 0) then
    tabref = aot_reference_for(L = conf, thandle = thandle, key = 'myfun')
    write(*,*) '   : got the reference: ', tabref
  else
    write(*,*) 'FATAL Error occured: could not open table test_table!'
    passed = .false.
  end if

  call aot_table_close(L = conf, thandle = thandle)

  ! Recalling the function and executing it
  res = 0.0
  write(*,*) ' * Recalling the function in the reference and executing it'
  call aot_fun_open(L = conf, fun=luafun, ref = tabref)
  write(*,*) '   : id: ', aot_fun_id(luafun)
  call aot_fun_put(L = conf, fun = luafun, arg = 2.0)
  call aot_fun_do(L = conf, fun = luafun, nresults = 1)
  call aot_top_get_val(L = conf, val = res, ErrCode = iError)

  if (btest(iError, aoterr_Fatal)) then
    write(*,*) '  : unexpected FATAL Error'
    write(*,*) '  : after executing the function occured !!!'
    if (btest(iError, aoterr_NonExistent)) &
      &   write(*,*) '  : Variable not existent!'
    if (btest(iError, aoterr_WrongType)) &
      &   write(*,*) '  : Variable has wrong type!'
    passed = .false.
  else
    if (nint(res) /= 4) then
      write(*,*) '  : Function evaluation returned wrong result, should be 4'
      write(*,*) '  : but is: ', res
      passed = .false.
    else
      write(*,*) '  : success'
    end if
  end if

  call aot_fun_close(L = conf, fun=luafun)

  write(*,*) ' * close_conf'
  call close_config(conf)
  write(*,*) '  : success.'
  write(*,*) '... Done with aot_ref_test.'
  if (passed) then
    write(*,*) 'PASSED'
  else
    write(*,*) 'FAILED'
  end if

contains

  subroutine create_script(filename)
    character(len=*) :: filename

    open(file=trim(filename), unit=22, action='write', status='replace')
    write(22,*) '-- test script for aot_ref_test'
    write(22,*) 'function squaring(x)'
    write(22,*) '  return x^2'
    write(22,*) 'end'
    write(22,*) "funref = squaring"
    write(22,*) "test_table = {"
    write(22,*) "  myfun = squaring"
    write(22,*) "}"
    close(22)
  end subroutine create_script

end program aot_ref_test
