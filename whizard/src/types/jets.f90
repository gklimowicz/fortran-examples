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

module jets

  use fastjet !NODEP!

  implicit none
  private

  public :: fastjet_available
  public :: fastjet_init
  public :: jet_definition_t
  public :: pseudojet_t
  public :: pseudojet_vector_t
  public :: cluster_sequence_t
  public :: assignment (=)

  public :: kt_algorithm
  public :: cambridge_algorithm
  public :: antikt_algorithm
  public :: genkt_algorithm
  public :: cambridge_for_passive_algorithm
  public :: genkt_for_passive_algorithm
  public :: ee_kt_algorithm
  public :: ee_genkt_algorithm
  public :: plugin_algorithm
  public :: undefined_jet_algorithm


contains

  subroutine fastjet_init ()
    call print_banner ()
  end subroutine fastjet_init


end module jets
