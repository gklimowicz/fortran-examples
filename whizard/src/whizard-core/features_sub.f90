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

submodule (features) features_s

  use string_utils, only: lower_case
  use system_dependencies, only: WHIZARD_VERSION
  use kinds, only: default
  use system_dependencies, only: openmp_is_active
  use system_dependencies, only: GOSAM_AVAILABLE
  use system_dependencies, only: OPENLOOPS_AVAILABLE
  use system_dependencies, only: RECOLA_AVAILABLE
  use system_dependencies, only: LHAPDF5_AVAILABLE
  use system_dependencies, only: LHAPDF6_AVAILABLE
  use system_dependencies, only: HOPPET_AVAILABLE
  use jets, only: fastjet_available
  use system_dependencies, only: PYTHIA6_AVAILABLE
  use system_dependencies, only: PYTHIA8_AVAILABLE
  use hepmc_interface, only: hepmc_is_available
  use lcio_interface, only: lcio_is_available
  use system_dependencies, only: EVENT_ANALYSIS

  implicit none

contains

  module subroutine print_features ()
    print "(A)", "WHIZARD " // WHIZARD_VERSION
    print "(A)", "Build configuration:"
    call print_check ("precision")
    print "(A)", "Optional features available in this build:"
    call print_check ("OpenMP")
    call print_check ("GoSam")
    call print_check ("OpenLoops")
    call print_check ("Recola")
    call print_check ("LHAPDF")
    call print_check ("HOPPET")
    call print_check ("fastjet")
    call print_check ("Pythia6")
    call print_check ("Pythia8")
    call print_check ("StdHEP")
    call print_check ("HepMC")
    call print_check ("LCIO")
    call print_check ("MetaPost")
  end subroutine print_features

  subroutine check (feature, recognized, result, help)
    character(*), intent(in) :: feature
    logical, intent(out) :: recognized
    character(*), intent(out) :: result, help
    recognized = .true.
    result = "no"
    select case (lower_case (trim (feature)))
    case ("precision")
       write (result, "(I0)")  precision (1._default)
       help = "significant decimals of real/complex numbers"
    case ("openmp")
       if (openmp_is_active ()) then
          result = "yes"
       end if
       help = "OpenMP parallel execution"
    case ("gosam")
       if (GOSAM_AVAILABLE) then
          result = "yes"
       end if
       help = "external NLO matrix element provider"
    case ("openloops")
       if (OPENLOOPS_AVAILABLE) then
          result = "yes"
       end if
       help = "external NLO matrix element provider"
    case ("recola")
       if (RECOLA_AVAILABLE) then
          result = "yes"
       end if
       help = "external NLO matrix element provider"
    case ("lhapdf")
       if (LHAPDF5_AVAILABLE) then
          result = "v5"
       else if (LHAPDF6_AVAILABLE) then
          result = "v6"
       end if
       help = "PDF library"
    case ("hoppet")
       if (HOPPET_AVAILABLE) then
          result = "yes"
       end if
       help = "PDF evolution package"
    case ("fastjet")
       if (fastjet_available ()) then
          result = "yes"
       end if
       help = "jet-clustering package"
    case ("pythia6")
       if (PYTHIA6_AVAILABLE) then
          result = "yes"
       end if
       help = "direct access for shower/hadronization"
    case ("pythia8")
       if (PYTHIA8_AVAILABLE) then
          result = "yes"
       end if
       help = "direct access for shower/hadronization"
    case ("stdhep")
       result = "yes"
       help = "event I/O format"
    case ("hepmc")
       if (hepmc_is_available ()) then
          result = "yes"
       end if
       help = "event I/O format"
    case ("lcio")
       if (lcio_is_available ()) then
          result = "yes"
       end if
       help = "event I/O format"
    case ("metapost")
       result = EVENT_ANALYSIS
       help = "graphical event analysis via LaTeX/MetaPost"
    case default
       recognized = .false.
    end select
  end subroutine check

  subroutine print_check (feature)
    character(*), intent(in) :: feature
    character(16) :: f
    logical :: recognized
    character(10) :: result
    character(48) :: help
    call check (feature, recognized, result, help)
    if (.not. recognized) then
       result = "unknown"
       help = ""
    end if
    f = feature
    print "(2x,A,1x,A,'(',A,')')", f, result, trim (help)
  end subroutine print_check


end submodule features_s

