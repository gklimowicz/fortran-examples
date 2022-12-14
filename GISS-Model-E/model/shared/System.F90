!@sum This file contains architecture specific code for SGI, IBM, Linux, DEC
  
! hack for PGI : use the same settings as G95 (except for iargc)
#ifdef COMPILER_PGI
#define COMPILER_G95
#endif

subroutine sys_flush (unit)
!@sum system call to flush corresponding I/O unit
!@auth I. Aleinov
#if ( defined(COMPILER_NAG) )
  use f90_unix_io
#endif
  implicit none
  integer, intent(IN) :: unit !@var unit
#if defined(MACHINE_SGI)
#if defined(COMPILER_G95)
  call flush(unit)
#else
  integer status
  call flush(unit,status)
#endif
#elif defined(MACHINE_Linux) || defined(MACHINE_DEC) \
  || ( defined(MACHINE_MAC) && ! defined(COMPILER_XLF) )
       call flush(unit) !!! should check if it works for Absoft and DEC
#elif defined( MACHINE_IBM ) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_XLF) )
  call flush_(unit)
#else
  none of supported architectures was specified.
  This will crash the compiling process.
#endif
  return
end subroutine sys_flush

subroutine sys_signal (sig, prog)
!@sum system call to "signal"
!@auth I. Aleinov
  implicit none
!@var unit signal number to catch
  integer, intent(IN) :: sig
!@var prog handler subroutine for given signal
  external prog
#if defined(MACHINE_SGI) \
  || ( defined(MACHINE_Linux) && ! defined(COMPILER_G95) && ! defined(COMPILER_NAG) ) \
       || defined(MACHINE_DEC) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_Intel8) ) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_ABSOFT) )
  call signal( sig, prog, -1 )
#elif defined( MACHINE_IBM ) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_XLF) ) \
  || ( defined(MACHINE_Linux) && defined(COMPILER_G95) ) \
  || ( defined(MACHINE_MAC) && defined(COMPILER_G95) )
  call signal( sig, prog )
#elif ( defined(COMPILER_NAG) )
  ! do nothing if "signal" is not supported by NAG
#else
  none of supported architectures was specified.
  This will crash the compiling process.
#endif
  return
end subroutine sys_signal


subroutine sys_abort
!@sum system call to "abort" (to dump core)
#if ( defined(COMPILER_NAG) )
  use f90_unix_proc
#endif
  call abort
end subroutine sys_abort


subroutine nextarg( arg, opt )
!@sum returns next argument on the command line
!@+  arg - returned argument, or returns "" if no more arguments
!@+  if opt==1 return arg only if it is an option (starts with -)
#if ( defined(COMPILER_NAG) )
  use f90_unix_env
#endif
  implicit none
  character(*), intent(out) :: arg
#if ((! defined(COMPILER_NAG) ) && (! defined(COMPILER_G95) )) || (defined COMPILER_PGI)
  integer, external :: iargc
#endif
  integer, intent(in) :: opt
  integer, save :: count = 1
  if ( count > iargc() ) then
    arg=""
    return
  endif
  call getarg( count, arg )
  !if ( present(opt) ) then
  if ( opt == 1 .and. arg(1:1) .ne. '-' ) then
    arg=""
    return
  endif
  !endif
  count = count + 1
  return
end subroutine nextarg
