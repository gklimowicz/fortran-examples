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

module physics_defs

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use constants, only: one, two, three

  implicit none
  private

  real(default), parameter, public :: &
       conv = 0.38937966e12_default
  real(default), parameter, public :: &
       ns_per_mm = 1.e6_default / 299792458._default
  real(default), parameter, public :: &
       pb_per_fb = 1.e-3_default
  character(*), parameter, public :: &
       energy_unit = "GeV"
  character(*), parameter, public :: &
       cross_section_unit = "fb"
  real(default), parameter, public :: &
       NC = three, &
       CF = (NC**2 - one) / two / NC, &
       CA = NC, &
       TR = one / two
  real(default), public, parameter :: MZ_REF = 91.188_default
  real(default), public, parameter :: ME_REF = 0.000510998928_default
  real(default), public, parameter :: ALPHA_QCD_MZ_REF = 0.1178_default
  real(default), public, parameter :: ALPHA_QED_ME_REF = 0.0072973525693_default
  real(default), public, parameter :: LAMBDA_QCD_REF = 200.e-3_default
  integer, parameter, public :: UNDEFINED = 0

  integer, parameter, public :: DOWN_Q = 1
  integer, parameter, public :: UP_Q = 2
  integer, parameter, public :: STRANGE_Q = 3
  integer, parameter, public :: CHARM_Q = 4
  integer, parameter, public :: BOTTOM_Q = 5
  integer, parameter, public :: TOP_Q = 6
  integer, parameter, public :: ELECTRON = 11
  integer, parameter, public :: ELECTRON_NEUTRINO = 12
  integer, parameter, public :: MUON = 13
  integer, parameter, public :: MUON_NEUTRINO = 14
  integer, parameter, public :: TAU = 15
  integer, parameter, public :: TAU_NEUTRINO = 16

  integer, parameter, public :: GLUON = 21
  integer, parameter, public :: PHOTON = 22
  integer, parameter, public :: PHOTON_OFFSHELL = -2002
  integer, parameter, public :: PHOTON_ONSHELL = 2002
  integer, parameter, public :: Z_BOSON = 23
  integer, parameter, public :: W_BOSON = 24

  integer, parameter, public :: PION = 111
  integer, parameter, public :: PIPLUS = 211
  integer, parameter, public :: PIMINUS = - PIPLUS

  integer, parameter, public :: UD0 = 2101
  integer, parameter, public :: UD1 = 2103
  integer, parameter, public :: UU1 = 2203

  integer, parameter, public :: K0L = 130
  integer, parameter, public :: K0S = 310
  integer, parameter, public :: K0 = 311
  integer, parameter, public :: KPLUS = 321
  integer, parameter, public :: DPLUS = 411
  integer, parameter, public :: D0 = 421
  integer, parameter, public :: B0 = 511
  integer, parameter, public :: BPLUS = 521

  integer, parameter, public :: PROTON = 2212
  integer, parameter, public :: NEUTRON = 2112
  integer, parameter, public :: DELTAPLUSPLUS = 2224
  integer, parameter, public :: DELTAPLUS = 2214
  integer, parameter, public :: DELTA0 = 2114
  integer, parameter, public :: DELTAMINUS = 1114

  integer, parameter, public :: SIGMAPLUS = 3222
  integer, parameter, public :: SIGMA0 = 3212
  integer, parameter, public :: SIGMAMINUS = 3112

  integer, parameter, public :: SIGMACPLUSPLUS = 4222
  integer, parameter, public :: SIGMACPLUS = 4212
  integer, parameter, public :: SIGMAC0 = 4112

  integer, parameter, public :: SIGMAB0 = 5212
  integer, parameter, public :: SIGMABPLUS = 5222

  integer, parameter, public :: BEAM_REMNANT = 9999
  integer, parameter, public :: HADRON_REMNANT = 90
  integer, parameter, public :: HADRON_REMNANT_SINGLET = 91
  integer, parameter, public :: HADRON_REMNANT_TRIPLET = 92
  integer, parameter, public :: HADRON_REMNANT_OCTET = 93

  integer, parameter, public :: INTERNAL = 94
  integer, parameter, public :: INVALID = 97

  integer, parameter, public :: COMPOSITE = 99

  integer, parameter, public:: UNKNOWN = 0
  integer, parameter, public :: SCALAR = 1, SPINOR = 2, VECTOR = 3, &
                                VECTORSPINOR = 4, TENSOR = 5

  integer, parameter, public :: BORN = 0
  integer, parameter, public :: NLO_REAL = 1
  integer, parameter, public :: NLO_VIRTUAL = 2
  integer, parameter, public :: NLO_MISMATCH = 3
  integer, parameter, public :: NLO_DGLAP = 4
  integer, parameter, public :: NLO_SUBTRACTION = 5
  integer, parameter, public :: NLO_FULL = 6
  integer, parameter, public :: GKS = 7
  integer, parameter, public :: COMPONENT_UNDEFINED = 99

  integer, parameter, public :: n_beams_rescaled = 2

  integer, parameter, public :: EPDF_LL = 0
  integer, parameter, public :: EPDF_NLL = 1

  integer, parameter, public :: THR_POS_WP = 3
  integer, parameter, public :: THR_POS_WM = 4
  integer, parameter, public :: THR_POS_B = 5
  integer, parameter, public :: THR_POS_BBAR = 6
  integer, parameter, public :: THR_POS_GLUON = 7

  integer, parameter, public :: THR_EMITTER_OFFSET = 4

  integer, parameter, public :: NO_FACTORIZATION = 0
  integer, parameter, public :: FACTORIZATION_THRESHOLD = 1

  integer, dimension(2), parameter, public :: ass_quark = [5, 6]
  integer, dimension(2), parameter, public :: ass_boson = [3, 4]

  integer, parameter, public :: PROC_MODE_UNDEFINED = 0
  integer, parameter, public :: PROC_MODE_TT = 1
  integer, parameter, public :: PROC_MODE_WBWB = 2


  public :: component_status
  public :: is_nlo_component
  public :: is_subtraction_component
  public :: thr_leg

  interface component_status
     module procedure component_status_of_string
     module procedure component_status_to_string
  end interface

  interface
    elemental module function component_status_of_string (string) result (i)
      integer :: i
      type(string_t), intent(in) :: string
    end function component_status_of_string
    elemental module function component_status_to_string (i) result (string)
      type(string_t) :: string
      integer, intent(in) :: i
    end function component_status_to_string
    elemental module function is_nlo_component (comp) result (is_nlo)
      logical :: is_nlo
      integer, intent(in) :: comp
    end function is_nlo_component
    module function is_subtraction_component (emitter, nlo_type) result (is_subtraction)
      logical :: is_subtraction
      integer, intent(in) :: emitter, nlo_type
    end function is_subtraction_component
    module function thr_leg (emitter) result (leg)
      integer :: leg
      integer, intent(in) :: emitter
    end function thr_leg
  end interface

end module physics_defs
