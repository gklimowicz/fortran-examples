Subroutine assert_(line, fname, msg)
#ifdef USE_ESMF
  Use ESMF, only: ESMF_Finalize
#endif
  Implicit None
  Integer,          Intent(In) :: line
  Character(Len=*), Intent(In) :: fname
  Character(Len=*), Intent(In) :: msg
  integer :: status

  Write(*,*) 'Assertion failed at line',line,'in file',fname
  Write(*,*) msg
#ifdef USE_ESMF
  Call ESMF_Finalize(rc=status)
#else
#ifdef USE_MPI
  call MPI_Finalize(status)
#endif
#endif
  Stop

End Subroutine assert_

Subroutine warn_(line, fname, msg)
  Implicit None
  Integer,          Intent(In) :: line
  Character(Len=*), Intent(In) :: fname
  Character(Len=*), Intent(In) :: msg
  
  Write(*,*) 'Warning: assertion failed at line',line,'in file',fname
  Write(*,*) msg
  Write(*,*) 'Continuing anyway ...'
  Write(*,*) ' '

End Subroutine warn_
