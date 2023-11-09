!!! Auxiliary program
!!! Read the histogram output from epa_2
!!! and count the number of events below and above Q threshold.

program count_evt
  implicit none

  character(*), parameter :: infile = "epa_2a.dat"
  character(*), parameter :: outfile = "epa_2_count.out"
  integer, parameter :: n_bin = 20
  character(100) :: buffer
  integer :: u_in, u_out
  real :: mid, val, err, exc
  integer :: i, n
  integer, dimension(n_bin) :: n_entry

  open (newunit = u_in, file = infile, action = "read", status = "old")

  n_entry = 0

  ! skip comment lines
  SKIP_COMMENT: do
     read (u_in, "(A)")  buffer
     if (buffer(1:1) /= "#")  exit SKIP_COMMENT
  end do SKIP_COMMENT
  
  ! read histogram entries
  READ_ENTRIES: do i = 1, n_bin
     if (buffer == "")  exit READ_ENTRIES
     read (buffer, *) mid, val, err, exc, n
     n_entry(i) = n
     read (u_in, "(A)")  buffer
  end do READ_ENTRIES

  close (u_in)
  
  open (newunit = u_out, file = outfile, action = "write", status = "replace")
  write (u_out, "(A)")  "Event breakdown:"
! Enable if Qmin is supported (currently not)
!   write (u_out, "(2x,A,1x,I3)")  "Q < Qmin:        ", sum (n_entry(1:1))
!   write (u_out, "(2x,A,1x,I3)")  "Q_min < Q < Qmax:", sum (n_entry(2:n_bin/2))
  write (u_out, "(2x,A,1x,I3)")  "Q < Qmax:", sum (n_entry(1:n_bin/2))
  write (u_out, "(2x,A,1x,I3)")  "Q > Qmax:", sum (n_entry(n_bin/2+1:n_bin))
  close (u_out)

end program count_evt
