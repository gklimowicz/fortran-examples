subroutine initModelEtests()
  use stop_model_mod, only: set_stop_model_ptr
  implicit none
  external stop_model_pfunit
  call set_stop_model_ptr(stop_model_pfunit)

end subroutine initModelEtests

subroutine stop_model_pfunit(message, retcode)
  use pfunit_mod
  implicit none
!@var message an error message (reason to stop)
  character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
  integer, intent(in) :: retcode
      
  call throw(message)
end subroutine stop_model_pfunit

