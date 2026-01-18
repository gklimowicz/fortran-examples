!!! Auxiliary program
!!! Read the weight-list output from resonances_15
!!! and compare the weights from the original run and the rescan.

program main
  use iso_fortran_env, only: real64

  implicit none

  character(*), parameter :: infile1 = "resonances_15_a.weights.dat"
  character(*), parameter :: infile2 = "resonances_15_b.weights.dat"
  character(*), parameter :: outfile = "resonances_15_check.out"

  integer :: u_in1, u_in2, u_out
  integer :: iostat

  integer :: idum
  real(real64) :: wdum, sqme1, sqme2, sqme3, ratio21, ratio32, tol1, tol2, tol3
  logical :: ok1, ok2, ok3

  open (newunit = u_in1, file = infile1, action = "read", status = "old")
  open (newunit = u_in2, file = infile2, action = "read", status = "old")
  open (newunit = u_out, file = outfile, action = "write", status = "replace")

  ok1 = .true.
  tol1 = 1e-5

  ok2 = .true.
  tol2 = 0.20

  ok3 = .true.
  tol3 = 1e-8

  do
     read (u_in1, *, iostat=iostat) idum, wdum, sqme1
     read (u_in2, *, iostat=iostat) idum, wdum, sqme2
     read (u_in2, *, iostat=iostat) idum, wdum, sqme3
     if (iostat /= 0) exit
     ratio21 = sqme2 / sqme1
     ratio32 = sqme3 / sqme2
     write (u_out, "(F12.10,2x,F12.10)")  ratio21, ratio32
     if (abs (ratio21 - 1) > tol1)  ok1 = .false.
     if (abs (ratio21 - 1) > tol2)  ok2 = .false.
     if (abs (ratio32 - 1) > tol3)  ok3 = .false.
  end do
     
  close (u_out)
  close (u_in2)
  close (u_in1)

  print *, "deviation rescan vs simulation < tolerance:", ok1
  print *, "deviation rescan vs simulation < 20%:      ", ok2
  print *, "deviation rescan alt_setup     < tolerance:", ok3

end program main
