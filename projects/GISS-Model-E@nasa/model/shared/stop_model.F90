#include "rundeck_opts.h"

module stop_model_mod

  implicit none
  private
  public :: set_stop_model_ptr
  public :: stop_model_generic, stop_model_default, stop_model_segfault
#ifdef USE_FEXCEPTION
  public :: stop_model_fexception
#endif

  interface
    subroutine stop_model_cb(message, retcode)
      implicit none
!@var message an error message (reason to stop)
      character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
      integer, intent(in) :: retcode
    end subroutine stop_model_cb
  end interface
  
  procedure(stop_model_cb), pointer, protected :: stop_model_ptr => null()

contains

  subroutine stop_model_generic( message, retcode )
!@sum Aborts the execution of the program. Passes an error message and
!@+ a return code to the calling script. Should be used instead of STOP
!@var message an error message (reason to stop)
    character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
    integer, intent(in) :: retcode
  
    if (.not. associated(stop_model_ptr)) then
      call set_stop_model_ptr(stop_model_default)
    end if
    call stop_model_ptr(message, retcode)

  end subroutine stop_model_generic

  subroutine set_stop_model_ptr(proc_ptr)
    procedure(stop_model_cb) :: proc_ptr
    stop_model_ptr => proc_ptr
  end subroutine set_stop_model_ptr

  subroutine stop_model_default( message, retcode )
!@sum Aborts the execution of the program. Passes an error message and
!@+ a return code to the calling script. Should be used instead of STOP
    use Dictionary_mod
!@var message an error message (reason to stop)
    character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
    integer, intent(in) :: retcode
    integer, parameter :: iu_err = 9
    integer :: rank
    logical :: flag_mpi
#ifdef USE_MPI
    integer :: mpi_err
#  ifdef MPI_DEFS_HACK
#  include "mpi_defs.h"
#  endif
#include "mpif.h"
#endif

    rank =0
#ifdef USE_MPI
    call MPI_Initialized(flag_mpi, mpi_err)
    if (flag_mpi) call MPI_COMM_RANK(MPI_COMM_WORLD, rank, mpi_err)
#endif
  ! skip writing status file for retcode<0
    if ( retcode >= 0 ) call write_run_status( message, retcode )
    if (rank == 0) then
      write (6,'(//2(" ",132("*")/))')
      write (6,*) ' Program terminated due to the following reason:'
      write (6,*) ' >>  ', message, '  <<'
      write (6,'(/2(" ",132("*")/))')
    endif

    call sys_flush(6)

    if ( retcode > 13 ) then
      write (0,*) 'Model crashed due to ',message
#ifdef USE_MPI
      if (flag_mpi) call mpi_abort(MPI_COMM_WORLD, retcode, iu_err)
#endif
      call sys_abort
    else
#ifdef USE_MPI
      if (flag_mpi) call mpi_finalize(mpi_err)
#endif
      call exit_rc (0)
  endif
  end subroutine stop_model_default


  subroutine stop_model_segfault( message, retcode )
!@sum Aborts the execution of the program. Passes an error message and
!@+ a return code to the calling script. Should be used instead of STOP
    use Dictionary_mod
!@var message an error message (reason to stop)
    character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
    integer, intent(in) :: retcode
    integer :: rank
    integer, pointer :: crash_me

  ! skip writing status file for retcode<0
    if ( retcode >= 0 ) call write_run_status( message, retcode )
    if (rank == 0) then
      write (6,'(//2(" ",132("*")/))')
      write (6,*) ' Program terminated due to the following reason:'
      write (6,*) ' >>  ', message, '  <<'
      write (6,'(/2(" ",132("*")/))')
    endif

    call sys_flush(6)

    if ( retcode > 13 ) then
      write (0,*) 'Model crashed due to ',message
    endif

    ! Cause a segfault, which will produce a stack trace
    nullify(crash_me)
    crash_me = 17

  end subroutine stop_model_segfault


#ifdef USE_FEXCEPTION
  subroutine stop_model_fexception(message, retcode)
    use fexception_mod  
!@var message an error message (reason to stop)
    character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
    integer, intent(in) :: retcode
      
    call throw(message, retcode)
      
  end subroutine stop_model_fexception
#endif

  subroutine exit_rc (code)
!@sum  exit_rc stops the run and sets a return code
!@auth Reto A Ruedy
#if ( defined(COMPILER_NAG) )
    use f90_unix_proc
#endif
    integer, intent(IN) :: code !@var code return code set by user
#if defined(MACHINE_SGI) || defined(MACHINE_Linux) || defined(MACHINE_DEC) \
  || ( defined(MACHINE_MAC) && ! defined(COMPILER_XLF) )
      call exit(code) !!! should check if it works for Absoft and DEC
#elif defined( MACHINE_IBM ) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_XLF) )
      call exit_(code)
#else
  none of supported architectures was specified.
  This will crash the compiling process.
#endif

  end subroutine exit_rc

end module stop_model_mod

! =====================================================================

subroutine stop_model( message, retcode )
!@sum Aborts the execution of the program. Passes an error message and
!@+ a return code to the calling script. Should be used instead of STOP
!@+ Method for stopping model is delegated to stop_module_mod.
  use stop_model_mod, only: stop_model_generic
  implicit none
!@var message an error message (reason to stop)
  character*(*), intent (in) :: message
!@var retcode return code to be passed to the calling script
  integer, intent(in) :: retcode
  
  call stop_model_generic(message, retcode)

end subroutine stop_model

