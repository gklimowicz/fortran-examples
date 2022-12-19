!!! Auxiliary program
!!! Read the event output (LHE format) from resonances_1
!!! and count the number of events for each resonance pattern.

program count_res
  implicit none

  character(*), parameter :: infile = "resonances_3_p.lhe"
  character(*), parameter :: outfile = "resonances_3_count.out"
  character(80) :: buffer
  integer :: u, iostat
  integer :: pdg, status
  integer :: ev_z
  integer :: n_evt
  integer :: tot_0, tot_1z

  n_evt = 0
  tot_0 = 0
  tot_1z = 0

  open (newunit = u, file = infile, action = "read", status = "old")
  SCAN_EVENTS: do

     ! Find <event> tag
     NEXT_EVENT: do
        read (u, "(A80)", iostat=iostat)  buffer
        if (iostat /= 0)  exit SCAN_EVENTS
        if (buffer(1:7) == "<event>")  exit NEXT_EVENT
     end do NEXT_EVENT

     ! reset counters for this event
     n_evt = n_evt + 1
     ev_z = 0

     ! skip event header line
     read (u, *)

     ! digest particle lines
     READ_PRT: do
        read (u, "(A80)")  buffer
        if (buffer(1:1) == "<")  exit READ_PRT
        read (buffer, *)  pdg, status
        select case (status)
        case (2)
           select case (pdg)
           case (23)
              ev_z = ev_z + 1
           end select
        end select
     end do READ_PRT

     ! analyze event
     select case (ev_z)
     case (1);  tot_1z = tot_1z + 1
     case default
        tot_0 = tot_0 + 1
     end select
  end do SCAN_EVENTS
  close (u)
  
  open (newunit = u, file = outfile, action = "write", status = "replace")
  write (u, "(A)")  "Event breakdown:"
  write (u, "(2x,A,1x,I3)")  "total:", n_evt
  write (u, "(2x,A,1x,I3)")  "Z:    ", tot_1z
  write (u, "(2x,A,1x,I3)")  "bkgd: ", tot_0
  close (u)

end program count_res
