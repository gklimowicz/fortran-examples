!!! Auxiliary program
!!! Read the event output (LHE format) from resonances_9
!!! and count the number of events for each resonance pattern.

program count_res
  implicit none

  character(*), parameter :: infile = "resonances_9_a.lhe"
  character(*), parameter :: outfile = "resonances_9_count.out"
  character(80) :: buffer
  integer :: u, iostat
  integer :: pdg, status
  integer :: ev_w, ev_z, ev_c
  integer :: n_evt
  integer :: tot_uddu, tot_udcs
  integer :: tot_0_u, tot_1w_u, tot_2w_u, tot_1z_u, tot_2z_u
  integer :: tot_0_c, tot_1w_c, tot_2w_c, tot_1z_c, tot_2z_c

  n_evt = 0
  tot_uddu = 0
  tot_0_u = 0
  tot_1w_u = 0
  tot_2w_u = 0
  tot_1z_u = 0
  tot_2z_u = 0
  tot_udcs = 0
  tot_0_c = 0
  tot_1w_c = 0
  tot_2w_c = 0
  tot_1z_c = 0
  tot_2z_c = 0

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
     ev_w = 0
     ev_z = 0
     ev_c = 0

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
           case (4,-4)
              ev_c = ev_c + 1
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
     select case (ev_c)
     case (0)
        tot_uddu = tot_uddu + 1
        select case (ev_w)
        case (1);  tot_1w_u = tot_1w_u + 1
        case (2);  tot_2w_u = tot_2w_u + 1
        case default
           select case (ev_z)
           case (1);  tot_1z_u = tot_1z_u + 1
           case (2);  tot_2z_u = tot_2z_u + 1
           case default
              tot_0_u = tot_0_u + 1
           end select
        end select
     case (1);  tot_udcs = tot_udcs + 1
        select case (ev_w)
        case (1);  tot_1w_c = tot_1w_c + 1
        case (2);  tot_2w_c = tot_2w_c + 1
        case default
           select case (ev_z)
           case (1);  tot_1z_c = tot_1z_c + 1
           case (2);  tot_2z_c = tot_2z_c + 1
           case default
              tot_0_c = tot_0_c + 1
           end select
        end select
     end select

  end do SCAN_EVENTS
  close (u)
  
  open (newunit = u, file = outfile, action = "write", status = "replace")
  write (u, "(A)")  "Event breakdown:"
  write (u, "(2x,A,1x,I3)")  "total:", n_evt
  write (u, "(2x,A,1x,I3)")  "uddu: ", tot_uddu
  write (u, "(2x,A,1x,I3)")  "  WW: ", tot_2w_u
  write (u, "(2x,A,1x,I3)")  "  W:  ", tot_1w_u
  write (u, "(2x,A,1x,I3)")  "  ZZ: ", tot_2z_u
  write (u, "(2x,A,1x,I3)")  "  Z:  ", tot_1z_u
  write (u, "(2x,A,1x,I3)")  "  bkg:", tot_0_u
  write (u, "(2x,A,1x,I3)")  "udcs: ", tot_udcs
  write (u, "(2x,A,1x,I3)")  "  WW: ", tot_2w_c
  write (u, "(2x,A,1x,I3)")  "  W:  ", tot_1w_c
  write (u, "(2x,A,1x,I3)")  "  ZZ: ", tot_2z_c
  write (u, "(2x,A,1x,I3)")  "  Z:  ", tot_1z_c
  write (u, "(2x,A,1x,I3)")  "  bkg:", tot_0_c
  close (u)

end program count_res
