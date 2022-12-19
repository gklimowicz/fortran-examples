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

module jets_uti

  use kinds, only: default
  use fastjet !NODEP!

  use jets

  implicit none
  private

  public :: jets_1

contains

  subroutine jets_1 (u)
    integer, intent(in) :: u

    type(pseudojet_t), dimension(:), allocatable :: prt, jets, constituents
    type(jet_definition_t) :: jet_def
    type(cluster_sequence_t) :: cs

    integer, parameter :: dp = default
    integer :: i, j

    write (u, "(A)")  "* Test output: jets_1"
    write (u, "(A)")  "*   Purpose: test basic FastJet functionality"
    write (u, "(A)")

    write (u, "(A)")  "* Print banner"
    call print_banner ()

    write (u, *)
    write (u, "(A)")  "* Prepare input particles"
    allocate (prt (3))
    call prt(1)%init ( 99._dp, 0.1_dp, 0._dp, 100._dp)
    call prt(2)%init (  4._dp,-0.1_dp, 0._dp,   5._dp)
    call prt(3)%init (-99._dp, 0._dp,  0._dp,  99._dp)

    write (u, *)
    write (u, "(A)")  "* Define jet algorithm"
    call jet_def%init (antikt_algorithm, 0.7_dp)

    write (u, *)
    write (u, "(A)")  "* Cluster particles according to jet algorithm"

    write (u, *)
    write (u, "(A,A)")  "Clustering with ", jet_def%description ()
    call cs%init (pseudojet_vector (prt), jet_def)

    write (u, *)
    write (u, "(A)")  "* Sort output jets"
    jets = sorted_by_pt (cs%inclusive_jets ())

    write (u, *)
    write (u, "(A)")  "* Print jet observables and constituents"
    write (u, *)
    write (u, "(4x,3(7x,A3))") "pt", "y", "phi"
    do i = 1, size (jets)
       write (u, "(A,1x,I0,A,3(1x,F9.5))") &
            "jet", i, ":", jets(i)%perp (), jets(i)%rap (), jets(i)%phi ()
       constituents = jets(i)%constituents ()
       do j = 1, size (constituents)
          write (u, "(4x,A,1x,I0,A,F9.5)") &
               "constituent", j, "'s pt:", constituents(j)%perp ()
       end do
       do j = 1, size (constituents)
          call constituents(j)%final ()
       end do
    end do

    write (u, *)
    write (u, "(A)")  "* Cleanup"

    do i = 1, size (prt)
       call prt(i)%final ()
    end do
    do i = 1, size (jets)
       call jets(i)%final ()
    end do
    call jet_def%final ()
    call cs%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: jets_1"

  end subroutine jets_1


end module jets_uti
