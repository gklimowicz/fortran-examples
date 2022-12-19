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

submodule (pdf) pdf_s

  use io_units
  use system_dependencies, only: LHAPDF5_AVAILABLE, LHAPDF6_AVAILABLE
  use diagnostics

  implicit none

contains

  module subroutine pdf_data_init (pdf_data, pdf_data_in)
    class(pdf_data_t), intent(out) :: pdf_data
    type(pdf_data_t), target, intent(in) :: pdf_data_in
    pdf_data%xmin = pdf_data_in%xmin
    pdf_data%xmax = pdf_data_in%xmax
    pdf_data%qmin = pdf_data_in%qmin
    pdf_data%qmax = pdf_data_in%qmax
    pdf_data%set = pdf_data_in%set
    pdf_data%type = pdf_data_in%type
    if (pdf_data%type == STRF_LHAPDF6) then
       if (pdf_data_in%pdf%is_associated ()) then
          call lhapdf_copy_pointer (pdf_data_in%pdf, pdf_data%pdf)
       else
          call msg_bug ('pdf_data_init: pdf_data%pdf was not associated!')
       end if
    end if
  end subroutine pdf_data_init

  module subroutine pdf_data_write (pdf_data, unit)
    class(pdf_data_t), intent(in) :: pdf_data
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit); if (u < 0) return
    write (u, "(3x,A,I0)") "PDF set  = ", pdf_data%set
    write (u, "(3x,A,I0)") "PDF type = ", pdf_data%type
  end subroutine pdf_data_write

  module subroutine pdf_data_setup &
       (pdf_data, caller, beam_structure, lhapdf_member, set)
    class(pdf_data_t), intent(inout) :: pdf_data
    character(len=*), intent(in) :: caller
    type(beam_structure_t), intent(in) :: beam_structure
    integer, intent(in) :: lhapdf_member, set
    real(default) :: xmin, xmax, q2min, q2max
    pdf_data%set = set
    if (beam_structure%contains ("lhapdf")) then
       if (LHAPDF6_AVAILABLE) then
          pdf_data%type = STRF_LHAPDF6
       else if (LHAPDF5_AVAILABLE) then
          pdf_data%type = STRF_LHAPDF5
       end if
       write (msg_buffer, "(A,I0)")  caller &
            // ": interfacing LHAPDF set #", pdf_data%set
       call msg_message ()
    else if (beam_structure%contains ("pdf_builtin")) then
       pdf_data%type = STRF_PDF_BUILTIN
       write (msg_buffer, "(A,I0)")  caller &
            // ": interfacing PDF builtin set #", pdf_data%set
       call msg_message ()
    end if
    select case (pdf_data%type)
    case (STRF_LHAPDF6)
       pdf_data%xmin = pdf_data%pdf%getxmin ()
       pdf_data%xmax = pdf_data%pdf%getxmax ()
       pdf_data%qmin = sqrt(pdf_data%pdf%getq2min ())
       pdf_data%qmax = sqrt(pdf_data%pdf%getq2max ())
    case (STRF_LHAPDF5)
       call GetXminM (1, lhapdf_member, xmin)
       call GetXmaxM (1, lhapdf_member, xmax)
       call GetQ2minM (1, lhapdf_member, q2min)
       call GetQ2maxM (1, lhapdf_member, q2max)
       pdf_data%xmin = xmin
       pdf_data%xmax = xmax
       pdf_data%qmin = sqrt(q2min)
       pdf_data%qmax = sqrt(q2max)
    end select
  end subroutine pdf_data_setup

  module subroutine pdf_data_evolve (pdf_data, x, q_in, f)
    class(pdf_data_t), intent(inout) :: pdf_data
    real(double), intent(in) :: x, q_in
    real(double), dimension(-6:6), intent(out) :: f
    real(double) :: q
    select case (pdf_data%type)
    case (STRF_PDF_BUILTIN)
       call pdf_evolve_LHAPDF (pdf_data%set, x, q_in, f)
    case (STRF_LHAPDF6)
       q = min (pdf_data%qmax, q_in)
       q = max (pdf_data%qmin, q)
       call pdf_data%pdf%evolve_pdfm (x, q, f)
    case (STRF_LHAPDF5)
       q = min (pdf_data%qmax, q_in)
       q = max (pdf_data%qmin, q)
       call evolvePDFM (pdf_data%set, x, q, f)
    case default
       call msg_fatal ("PDF function: unknown PDF method.")
    end select
  end subroutine pdf_data_evolve

end submodule pdf_s

