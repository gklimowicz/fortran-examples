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

module phs_fks_ut
  use unit_tests
  use phs_fks_uti

  implicit none
  private

  public :: phs_fks_generator_test

contains

  subroutine phs_fks_generator_test (u, results)
    integer, intent(in) :: u
    type(test_results_t), intent(inout) :: results
    call test(phs_fks_generator_1, "phs_fks_generator_1", &
         "Test the generation of FKS phase spaces", u, results)
    call test(phs_fks_generator_2, "phs_fks_generator_2", &
         "Test the generation of an ISR FKS phase space", u, results)
    call test(phs_fks_generator_3, "phs_fks_generator_3", &
         "Test the generation of a real phase space for decays", &
         u, results)
    call test(phs_fks_generator_4, "phs_fks_generator_4", &
         "Test the generation of an FSR phase space with "&
         &"conserved invariant resonance masses", u, results)
    call test(phs_fks_generator_5, "phs_fks_generator_5", &
         "Test on-shell projection of a Born phase space and the generation"&
         &" of a real phase-space from that", u, results)
    call test(phs_fks_generator_6, "phs_fks_generator_6", &
         "Test the generation of a real phase space for 1 -> 3 decays", &
         u, results)
    call test(phs_fks_generator_7, "phs_fks_generator_7", &
         "Test the generation of an ISR FKS phase space for fixed beam energy", &
         u, results)
  end subroutine phs_fks_generator_test


end module phs_fks_ut
