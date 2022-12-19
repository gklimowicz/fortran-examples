!!! Auxiliary program
!!! Read the event output (LHE format) from resonances_4
!!! and count the number of events for each resonance pattern.

program count_res
  implicit none

  character(*), parameter :: infile = "resonances_4_p.lhe"
  character(*), parameter :: outfile = "resonances_4_count.out"
  character(80) :: buffer
  integer :: u, iostat
  integer :: pdg, status
  integer :: ev_w, ev_z, ev_u, ev_d
  integer :: n_evt
  integer :: tot_0, tot_1w, tot_2w, tot_1z, tot_2z
  integer :: tot_uuuu, tot_uudd, tot_dddd

  n_evt = 0
  tot_0 = 0
  tot_1w = 0
  tot_2w = 0
  tot_1z = 0
  tot_2z = 0
  tot_uuuu = 0
  tot_uudd = 0
  tot_dddd = 0

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
     ev_u = 0
     ev_d = 0
     ev_w = 0
     ev_z = 0

     ! skip event header line
     read (u, *)

     ! digest particle lines
     READ_PRT: do
        read (u, "(A80)")  buffer
        if (buffer(1:1) == "<")  exit READ_PRT
        read (buffer, *)  pdg, status
        select case (status)
        case (1)
           select case (pdg)
           case (1, -1)
              ev_d = ev_d + 1
           case (2, -2)
              ev_u = ev_u + 1
           end select
        case (2)
           select case (pdg)
           case (24, -24)
              ev_w = ev_w + 1
           case (23)
             ev_z = ev_z + 1
           end select
        end select
     end do READ_PRT

     ! analyze event
     select case (ev_u)
     case (0)
        select case (ev_d)
        case (4)
           tot_dddd = tot_dddd + 1
        end select
     case (2)
        select case (ev_d)
        case (2)
           tot_uudd = tot_uudd + 1
        end select
     case (4)
        select case (ev_d)
        case (0)
           tot_uuuu = tot_uuuu + 1
        end select
     end select

     select case (ev_w)
     case (1);  tot_1w = tot_1w + 1
     case (2);  tot_2w = tot_2w + 1
     case default
        select case (ev_z)
        case (1);  tot_1z = tot_1z + 1
        case (2);  tot_2z = tot_2z + 1
        case default
           tot_0 = tot_0 + 1
        end select
     end select
  end do SCAN_EVENTS
  close (u)
  
  open (newunit = u, file = outfile, action = "write", status = "replace")
  write (u, "(A)")  "Event breakdown:"
  write (u, "(2x,A,1x,I4)")  "total:", n_evt
  write (u, "(2x,A,1x,I4)")  "uuuu: ", tot_uuuu
  write (u, "(2x,A,1x,I4)")  "uudd: ", tot_uudd
  write (u, "(2x,A,1x,I4)")  "dddd: ", tot_dddd
  write (u, "(2x,A,1x,I4)")  "WW:   ", tot_2w
  write (u, "(2x,A,1x,I4)")  "W:    ", tot_1w
  write (u, "(2x,A,1x,I4)")  "ZZ:   ", tot_2z
  write (u, "(2x,A,1x,I4)")  "Z:    ", tot_1z
  write (u, "(2x,A,1x,I4)")  "bkgd: ", tot_0
  close (u)

end program count_res
