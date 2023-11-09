! Example Fortran code for calling WHIZARD via the Fortran API
! For detailed explanations, cf. the WHIZARD manual
!    
! ########################################################################
! #
! # Copyright (C) 1999-2020 by 
! #     Wolfgang Kilian <kilian@physik.uni-siegen.de>
! #     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
! #     Juergen Reuter <juergen.reuter@desy.de>
! #     with contributions from
! #     cf. main AUTHORS file
! #
! # WHIZARD is free software; you can redistribute it and/or modify it
! # under the terms of the GNU General Public License as published by 
! # the Free Software Foundation; either version 2, or (at your option)
! # any later version.
! #
! # WHIZARD is distributed in the hope that it will be useful, but
! # WITHOUT ANY WARRANTY; without even the implied warranty of
! # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! # GNU General Public License for more details.
! #
! # You should have received a copy of the GNU General Public License
! # along with this program; if not, write to the Free Software
! # Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
! #
! ########################################################################

program main

  ! WHIZARD API as a module
  use api

  ! Standard numeric types
  use iso_fortran_env, only: real64, int32

  implicit none

  ! WHIZARD and event-sample objects
  type(whizard_api_t) :: whizard
  type(simulation_api_t) :: sample

  ! Local variables
  real(real64) :: integral, error
  real(real64) :: sqme, weight
  integer(int32) :: idx
  integer(int32) :: i, it_begin, it_end
  
  ! Initialize WHIZARD, setting some global option
  call whizard%option ("model", "QED")
  call whizard%init ()

  ! Define a process, set some variables
  call whizard%command ("process mupair = e1, E1 => e2, E2")
  call whizard%set_var ("sqrts", 100._real64)
  call whizard%set_var ("seed", 0)
  
  ! Integrate and retrieve result
  call whizard%command ("integrate (mupair)")
  call whizard%get_integration_result ("mupair", integral, error)
  
  ! Print result
  print 1, "cross section =", integral / 1000, "pb"
  print 1, "error =", error / 1000, "pb"
1 format (2x,A,1x,F5.1,1x,A)
2 format (2x,A,1x,L1)

  ! Settings for event generation
  call whizard%set_var ("$sample", "mupair_events")
  call whizard%set_var ("n_events", 2)

  ! Create an event-sample object and generate events
  call whizard%new_sample ("mupair", sample)
  call sample%open (it_begin, it_end)
  do i = it_begin, it_end
     call sample%next_event ()
     call sample%get_event_index (idx)
     call sample%get_weight (weight)
     call sample%get_sqme (sqme)
     print "(A,I0)", "Event #", idx
     print 3, "sqme =", sqme
     print 3, "weight =", weight
3    format (2x,A,1x,ES10.3)
  end do

  ! Finalize the event-sample object
  call sample%close ()

  ! Finalize the WHIZARD object
  call whizard%final ()

end program main
