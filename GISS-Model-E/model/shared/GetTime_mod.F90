! Use F90 system_clock for portable accuracy
module GetTime_mod
  implicit none
  private
  public :: GetTime
contains

  subroutine GetTime(ctime, cmax)
    implicit none
    real*8, intent(out) :: ctime
    real*8, intent(out), optional :: cmax
    integer :: counter, count_rate, count_max
    call system_clock(counter,count_rate,count_max)
    if( present(cmax) ) cmax = count_max/real(count_rate,kind=8)
    ctime=counter/real(count_rate,kind=8)  ! convert to seconds
  end subroutine GetTime

end module GetTime_mod


