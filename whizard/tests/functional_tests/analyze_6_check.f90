!!! Auxiliary program
!!! Read the logfile and the histogram output from analyze_6
!!! and check whether the relevant entry matches the expectation, within uncertainty

program main
  implicit none

  character(*), parameter :: infile = "analyze_6.dat"
  character(*), parameter :: logfile = "analyze_6.log"
  character(*), parameter :: outfile = "analyze_6_check.out"

  character(80) :: buffer
  real :: bin, val, f
  integer :: u_in, u_log, u_out, i, j, pos

  open (newunit = u_in, file = infile, action = "read", status = "old")
  open (newunit = u_log, file = logfile, action = "read", status = "old")
  open (newunit = u_out, file = outfile, action = "write", status = "replace")

  ! Five histograms (analyze_6_a .. e)
  do i = 1, 5
     ! Find and read the relevant histogram entry
     SKIP_TO_NEXT_HIST: do
        read (u_in, "(A)") buffer
        if (buffer(1:6) == "# Hist")  exit SKIP_TO_NEXT_HIST
     end do SKIP_TO_NEXT_HIST
     read (u_in, *)             ! extra header line
     read (u_in, *)  bin, val   ! bin midpoint, bin value, error

     ! Find and read the rescaling factor, if any
     SKIP_TO_NEXT_SIMULATION: do
        read (u_log, "(A)") buffer
        if (buffer(1:21) == "| Starting simulation")  exit SKIP_TO_NEXT_SIMULATION
     end do SKIP_TO_NEXT_SIMULATION
     f = 1
     LOOK_FOR_DROPPED: do j = 1, 10
        read (u_log, "(A)") buffer
        if (buffer(1:6) == "| Drop") then
           read (u_log, "(A)")  buffer
           pos = scan (buffer, "=")
           read (buffer(pos+1:), *)  f
           exit LOOK_FOR_DROPPED
        end if
     end do LOOK_FOR_DROPPED

     ! Print the rescaled value, rounded to one decimal digit.
     write (u_out, "(A,1x,ES7.1)") &
          "Normalized entry, should be 1.0E+03:", pacify (val * f)
  end do

  close (u_in)
  close (u_log)
  close (u_out)

contains

  ! 1,000 with error margin +- 10% 
  function pacify (val) result (r)
    real, intent(in) :: val
    real :: r

    r = nint (val / 200) * 200.

  end function pacify
  
end program main
