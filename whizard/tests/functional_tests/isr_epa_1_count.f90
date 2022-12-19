!!! Auxiliary program
!!! Read the histogram output from isr_epa_1
!!! and count the number of events below and above Q threshold.

program count_evt
  implicit none

  character(*), parameter :: infile = "isr_epa_1a.dat"
  character(*), parameter :: outfile = "isr_epa_1_count.out"
  integer, parameter :: n_bin = 20
  character(100) :: buffer
  integer :: u_in, u_out
  real :: mid, val, err, exc
  integer :: i, n
  integer, dimension(n_bin) :: n_entry

  open (newunit = u_in, file = infile, action = "read", status = "old")

  n_entry = 0

  ! find next histogram
  FIND_HIST1: do
     read (u_in, "(A)")  buffer
     if (buffer(1:11) == "# Histogram")  exit FIND_HIST1
  end do FIND_HIST1

  ! skip comment lines
  SKIP_COMMENT1: do
     read (u_in, "(A)")  buffer
     if (buffer(1:1) /= "#")  exit SKIP_COMMENT1
  end do SKIP_COMMENT1
  
  ! read histogram entries (ISR)
  READ_ENTRIES1: do i = 1, n_bin
     if (buffer == "")  exit READ_ENTRIES1
     read (buffer, *) mid, val, err, exc, n
     n_entry(i) = n
     read (u_in, "(A)")  buffer
  end do READ_ENTRIES1
  
  open (newunit = u_out, file = outfile, action = "write", status = "replace")
  write (u_out, "(A)")  "Event breakdown (ISR):"
  write (u_out, "(2x,A,1x,I3)")  "Q < Qmax:", sum (n_entry(1:n_bin/2))
  write (u_out, "(2x,A,1x,I3)")  "Q > Qmax:", sum (n_entry(n_bin/2+1:n_bin))
 
  n_entry = 0

  ! find next histogram
  FIND_HIST2: do
     read (u_in, "(A)")  buffer
     if (buffer(1:11) == "# Histogram")  exit FIND_HIST2
  end do FIND_HIST2

  ! skip comment lines
  SKIP_COMMENT2: do
     read (u_in, "(A)")  buffer
     if (buffer(1:1) /= "#")  exit SKIP_COMMENT2
  end do SKIP_COMMENT2
  
  ! read histogram entries (ISR)
  READ_ENTRIES2: do i = 1, n_bin
     if (buffer == "")  exit READ_ENTRIES2
     read (buffer, *) mid, val, err, exc, n
     n_entry(i) = n
     read (u_in, "(A)")  buffer
  end do READ_ENTRIES2
  
  close (u_in)

  write (u_out, "(A)")  "Event breakdown (EPA):"
  write (u_out, "(2x,A,1x,I3)")  "Q < Qmax:", sum (n_entry(1:n_bin/4))
  write (u_out, "(2x,A,1x,I3)")  "Q > Qmax:", sum (n_entry(n_bin/4+1:n_bin))
  close (u_out)

end program count_evt
