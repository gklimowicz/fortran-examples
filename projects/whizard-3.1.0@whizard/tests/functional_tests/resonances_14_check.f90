!!! Auxiliary program
!!! Read the weight-list output from resonances_14
!!! and compare the weights from the original run and the rescan.

program main
  use iso_fortran_env, only: real64

  implicit none

  character(*), parameter :: infile1 = "resonances_14_a.weights.dat"
  character(*), parameter :: infile2 = "resonances_14_b.weights.dat"
  character(*), parameter :: outfile = "resonances_14_check.out"

  integer :: u_in1, u_in2, u_out
  integer :: iostat

  integer :: idum
  real(real64) :: wdum, sqme1, sqme2, ratio, tolerance
  logical :: ok

  open (newunit = u_in1, file = infile1, action = "read", status = "old")
  open (newunit = u_in2, file = infile2, action = "read", status = "old")
  open (newunit = u_out, file = outfile, action = "write", status = "replace")

  ok = .true.
  tolerance = 1e-5

  do
     read (u_in1, *, iostat=iostat) idum, wdum, sqme1
     read (u_in2, *, iostat=iostat) idum, wdum, sqme2
     if (iostat /= 0) exit
     ratio = sqme2 / sqme1
     write (u_out, "(F12.10)")  ratio
     if (abs (ratio - 1) > tolerance)  ok = .false.
  end do
     
  close (u_out)
  close (u_in2)
  close (u_in1)

  print *, "weight ratio deviation < tolerance:", ok

end program main
