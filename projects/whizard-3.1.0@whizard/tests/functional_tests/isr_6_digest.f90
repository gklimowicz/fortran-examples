program main
  use iso_fortran_env, only: default => real64
  implicit none

  character(*), parameter :: infile = "isr_6.dat"
  character(*), parameter :: outfile = "isr_6_digest.out"
  integer :: u_in, u_out
  integer :: iostat
  character(10) :: key
  character(19) :: val
  real(default) :: avg

  open (newunit = u_in, file = infile, action = "read", status = "old")
  open (newunit = u_out, file = outfile, action = "write", status = "replace")

  READ_LINES: do
     read (u_in, "(A10, 4x, A19)", iostat = iostat)  key, val
     if (iostat /= 0)  exit READ_LINES
     if (key(1:1) == "#")  cycle READ_LINES
     if (key == "average   ") then
        read (val, *)  avg
        write (u_out, "(F10.7)") abs (avg)
     end if
  end do READ_LINES

  close (u_in)
  close (u_out)

end program main
