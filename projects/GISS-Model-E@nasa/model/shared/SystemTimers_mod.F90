! A simple wrapper module for Igor's low-overhead timers based upon
! sysusage.
module SystemTimers_mod
  implicit none
  private

  public :: initializeSysTimers
  public :: startSysTimer
  public :: stopSysTimer
  public :: printSysTimers

  integer, parameter :: SYS_INITIALIZE = 0
  integer, parameter :: SYS_START = 1
  integer, parameter :: SYS_STOP  = 2
  integer, parameter :: SYS_PRINT = 3

  integer, parameter :: MAX_SU = 3

contains

  subroutine initializeSysTimers()
    integer :: timerIdx

#ifdef USE_SYSUSAGE
    do timerIdx = 0, MAX_SU
      call sysusage(timerIdx,SYS_INITIALIZE)
    enddo
#endif
  end subroutine initializeSysTimers

  subroutine startSysTimer(timerIdx)
    integer, intent(in) :: timerIdx

#ifdef USE_SYSUSAGE
    call sysusage(timerIdx,SYS_START)
#endif
  end subroutine startSysTimer

  subroutine stopSysTimer(timerIdx)
    integer, intent(in) :: timerIdx

#ifdef USE_SYSUSAGE
    call sysusage(timerIdx,SYS_STOP)
#endif
  end subroutine stopSysTimer

  subroutine printSysTimers()

#ifdef USE_SYSUSAGE
    do timerIdx = 0, MAX_SU
      call sysusage(timerIdx,SYS_PRINT)
    end do
#endif
  end subroutine printSysTimers

end module SystemTimers_mod
