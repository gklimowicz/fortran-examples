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

module sf_lhapdf_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use system_dependencies, only: LHAPDF5_AVAILABLE
  use system_dependencies, only: LHAPDF6_AVAILABLE
  use os_interface
  use physics_defs, only: PROTON
  use sm_qcd
  use lorentz
  use pdg_arrays
  use flavors
  use interactions, only: reset_interaction_counter
  use model_data
  use sf_base

  use sf_lhapdf

  implicit none
  private

  public :: sf_lhapdf_1
  public :: sf_lhapdf_2
  public :: sf_lhapdf_3

contains

  subroutine sf_lhapdf_1 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(pdg_array_t) :: pdg_in
    type(pdg_array_t), dimension(1) :: pdg_out
    integer, dimension(:), allocatable :: pdg1
    class(sf_data_t), allocatable :: data

    write (u, "(A)")  "* Test output: sf_lhapdf_1"
    write (u, "(A)")  "*   Purpose: initialize and display &
         &test structure function data"
    write (u, "(A)")

    write (u, "(A)")  "* Create empty data object"
    write (u, "(A)")

    call model%init_sm_test ()
    pdg_in = PROTON

    allocate (lhapdf_data_t :: data)
    call data%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize"
    write (u, "(A)")

    select type (data)
    type is (lhapdf_data_t)
       call data%init (model, pdg_in)
    end select

    call data%write (u)

    write (u, "(A)")

    write (u, "(1x,A)")  "Outgoing particle codes:"
    call data%get_pdg_out (pdg_out)
    pdg1 = pdg_out(1)
    write (u, "(2x,99(1x,I0))")  pdg1

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_lhapdf_1"

  end subroutine sf_lhapdf_1

  subroutine sf_lhapdf_2 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(flavor_t) :: flv
    type(pdg_array_t) :: pdg_in
    class(sf_data_t), allocatable, target :: data
    class(sf_int_t), allocatable :: sf_int
    type(vector4_t) :: k
    type(vector4_t), dimension(2) :: q
    real(default) :: E
    real(default), dimension(:), allocatable :: r, rb, x, xb
    real(default) :: f

    write (u, "(A)")  "* Test output: sf_lhapdf_2"
    write (u, "(A)")  "*   Purpose: initialize and fill &
         &test structure function object"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call model%init_sm_test ()
    call flv%init (PROTON, model)
    pdg_in = PROTON
    call lhapdf_global_reset ()

    call reset_interaction_counter ()

    allocate (lhapdf_data_t :: data)
    select type (data)
    type is (lhapdf_data_t)
       call data%init (model, pdg_in)
    end select

    write (u, "(A)")  "* Initialize structure-function object"
    write (u, "(A)")

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])

    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Initialize incoming momentum with E=500"
    write (u, "(A)")
    E = 500
    k = vector4_moving (E, sqrt (E**2 - flv%get_mass ()**2), 3)
    call vector4_write (k, u)
    call sf_int%seed_kinematics ([k])

    write (u, "(A)")
    write (u, "(A)")  "* Set kinematics for x=0.5"
    write (u, "(A)")

    allocate (r (data%get_n_par ()))
    allocate (rb(size (r)))
    allocate (x (size (r)))
    allocate (xb(size (r)))

    r = 0.5_default
    rb = 1 - r
    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%write (u)

    write (u, "(A)")
    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb
    write (u, "(A,9(1x,F10.7))")  "f =", f

    write (u, "(A)")
    write (u, "(A)")  "* Recover x from momenta"
    write (u, "(A)")

    q = sf_int%get_momenta (outgoing=.true.)
    call sf_int%final ()
    deallocate (sf_int)

    call data%allocate_sf_int (sf_int)
    call sf_int%init (data)
    call sf_int%set_beam_index ([1])

    call sf_int%seed_kinematics ([k])
    call sf_int%set_momenta (q, outgoing=.true.)
    call sf_int%recover_x (x, xb)

    write (u, "(A,9(1x,F10.7))")  "x =", x
    write (u, "(A,9(1x,F10.7))")  "xb=", xb

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate for Q = 100 GeV"
    write (u, "(A)")

    call sf_int%complete_kinematics (x, xb, f, r, rb, map=.false.)
    call sf_int%apply (scale = 100._default)
    call sf_int%write (u, testflag = .true.)


    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call sf_int%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_lhapdf_2"

  end subroutine sf_lhapdf_2

  subroutine sf_lhapdf_3 (u)
    integer, intent(in) :: u
    type(qcd_t) :: qcd
    type(string_t) :: name, path
    integer :: member

    write (u, "(A)")  "* Test output: sf_lhapdf_3"
    write (u, "(A)")  "*   Purpose: initialize and evaluate alpha_s"
    write (u, "(A)")

    write (u, "(A)")  "* Initialize configuration data"
    write (u, "(A)")

    call lhapdf_global_reset ()

    if (LHAPDF5_AVAILABLE) then
       name = "cteq6ll.LHpdf"
       member = 1
       path = ""
    else if (LHAPDF6_AVAILABLE) then
       name = "CT10"
       member = 1
       path = ""
    end if

    write (u, "(A)")  "* Initialize qcd object"
    write (u, "(A)")

    allocate (alpha_qcd_lhapdf_t :: qcd%alpha)
    select type (alpha => qcd%alpha)
    type is (alpha_qcd_lhapdf_t)
       call alpha%init (name, member, path)
    end select
    call qcd%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate for Q = 100"
    write (u, "(A)")

    write (u, "(1x,A,F8.5)")  "alpha = ", qcd%alpha%get (100._default)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: sf_lhapdf_3"

  end subroutine sf_lhapdf_3


end module sf_lhapdf_uti
