! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'

submodule (eio_data) eio_data_s

  use io_units
  use numeric_utils
  use diagnostics

  implicit none

contains

  module subroutine event_sample_data_init (data, n_proc, n_alt)
    class(event_sample_data_t), intent(out) :: data
    integer, intent(in) :: n_proc
    integer, intent(in), optional :: n_alt
    data%n_proc = n_proc
    allocate (data%proc_num_id (n_proc), source = 0)
    allocate (data%cross_section (n_proc), source = 0._default)
    allocate (data%error (n_proc), source = 0._default)
    if (present (n_alt)) then
       data%n_alt = n_alt
       allocate (data%md5sum_alt (n_alt))
       data%md5sum_alt = ""
    end if
  end subroutine event_sample_data_init

  module subroutine event_sample_data_write (data, unit)
    class(event_sample_data_t), intent(in) :: data
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Event sample properties:"
    write (u, "(3x,A,A,A)")  "MD5 sum (proc)   = '", data%md5sum_prc, "'"
    write (u, "(3x,A,A,A)")  "MD5 sum (config) = '", data%md5sum_cfg, "'"
    write (u, "(3x,A,L1)")  "unweighted       = ", data%unweighted
    write (u, "(3x,A,L1)")  "negative weights = ", data%negative_weights
    write (u, "(3x,A,A)")   "normalization    = ", &
         char (event_normalization_string (data%norm_mode))
    write (u, "(3x,A,I0)")  "number of beams  = ", data%n_beam
    write (u, "(5x,A,2(1x,I19))")  "PDG    = ", &
         data%pdg_beam(:data%n_beam)
    write (u, "(5x,A,2(1x,ES19.12))")  "Energy = ", &
         data%energy_beam(:data%n_beam)
    if (data%n_evt > 0) then
       write (u, "(3x,A,I0)")  "number of events = ", data%n_evt
    end if
    if (.not. vanishes (data%total_cross_section)) then
       write (u, "(3x,A,ES19.12)")  "total cross sec. = ", &
            data%total_cross_section
    end if
    write (u, "(3x,A,I0)")  "num of processes = ", data%n_proc
    do i = 1, data%n_proc
       write (u, "(3x,A,I0)")  "Process #", data%proc_num_id (i)
       select case (data%n_beam)
       case (1)
          write (u, "(5x,A,ES19.12)")  "Width = ", data%cross_section(i)
       case (2)
          write (u, "(5x,A,ES19.12)")  "CSec  = ", data%cross_section(i)
       end select
       write (u, "(5x,A,ES19.12)")  "Error = ", data%error(i)
    end do
    if (data%n_alt > 0) then
       write (u, "(3x,A,I0)")  "num of alt wgt   = ", data%n_alt
       do i = 1, data%n_alt
          write (u, "(5x,A,A,A,1x,I0)")  "MD5 sum (cfg)  = '", &
               data%md5sum_alt(i), "'", i
       end do
    end if
  end subroutine event_sample_data_write


end submodule eio_data_s

