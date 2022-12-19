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

module resonances_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use format_defs, only: FMF_12
  use lorentz, only: vector4_t, vector4_at_rest
  use model_data, only: model_data_t
  use flavors, only: flavor_t

  use resonances, only: resonance_history_t

  use resonances

  implicit none
  private

  public :: resonances_1
  public :: resonances_2
  public :: resonances_3
  public :: resonances_4
  public :: resonances_5
  public :: resonances_6
  public :: resonances_7

contains

  subroutine resonances_1 (u)
    integer, intent(in) :: u
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(model_data_t), target :: model

    write (u, "(A)")  "* Test output: resonances_1"
    write (u, "(A)")  "*   Purpose: test resonance history setup"
    write (u, "(A)")

    write (u, "(A)")  "* Read model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Empty resonance history"
    write (u, "(A)")

    call res_history%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Add resonance"
    write (u, "(A)")

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)


    write (u, "(A)")
    write (u, "(A)")  "* Add another resonance"
    write (u, "(A)")

    call res_info%init (7, 23, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Remove resonance"
    write (u, "(A)")

    call res_history%remove_resonance (1)
    call res_history%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonances_1"

  end subroutine resonances_1

  subroutine resonances_2 (u)
    integer, intent(in) :: u
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(model_data_t), target :: model
    type(string_t) :: restrictions

    write (u, "(A)")  "* Test output: resonances_2"
    write (u, "(A)")  "*   Purpose: test OMega restrictions strings &
         &for resonance history"
    write (u, "(A)")

    write (u, "(A)")  "* Read model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Empty resonance history"
    write (u, "(A)")

    restrictions = res_history%as_omega_string (2)
    write (u, "(A,A,A)")  "restrictions = '", char (restrictions), "'"

    write (u, "(A)")
    write (u, "(A)")  "* Add resonance"
    write (u, "(A)")

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    restrictions = res_history%as_omega_string (2)
    write (u, "(A,A,A)")  "restrictions = '", char (restrictions), "'"

    write (u, "(A)")
    write (u, "(A)")  "* Add another resonance"
    write (u, "(A)")

    call res_info%init (7, 23, model, 5)
    call res_history%add_resonance (res_info)
    restrictions = res_history%as_omega_string (2)
    write (u, "(A,A,A)")  "restrictions = '", char (restrictions), "'"

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonances_2"

  end subroutine resonances_2

  subroutine resonances_3 (u)
    integer, intent(in) :: u
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(resonance_history_t), dimension(:), allocatable :: res_histories
    type(resonance_history_set_t) :: res_set
    type(model_data_t), target :: model
    integer :: i

    write (u, "(A)")  "* Test output: resonances_3"
    write (u, "(A)")  "*   Purpose: test resonance history set"
    write (u, "(A)")

    write (u, "(A)")  "* Read model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Initialize resonance history set"
    write (u, "(A)")

    call res_set%init (initial_size = 2)

    write (u, "(A)")  "* Add resonance histories, one at a time"
    write (u, "(A)")

    call res_history%write (u)
    call res_set%enter (res_history)
    call res_history%clear ()

    write (u, *)

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)
    call res_set%enter (res_history)
    call res_history%clear ()

    write (u, *)

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_info%init (7, 23, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)
    call res_set%enter (res_history)
    call res_history%clear ()

    write (u, *)

    call res_info%init (7, 23, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)
    call res_set%enter (res_history)
    call res_history%clear ()

    write (u, *)

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)
    call res_set%enter (res_history)
    call res_history%clear ()

    write (u, *)

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_info%init (7, 25, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_set%freeze ()

    write (u, "(A)")
    write (u, "(A)")  "* Result"
    write (u, "(A)")

    call res_set%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Queries"
    write (u, "(A)")

    write (u, "(A,1x,I0)")  "n_history =", res_set%get_n_history ()

    write (u, "(A)")
    write (u, "(A)")  "History #2:"

    res_history = res_set%get_history (2)
    call res_history%write (u, indent=1)
    call res_history%clear ()

    write (u, "(A)")
    write (u, "(A)")  "* Result in array form"

    call res_set%to_array (res_histories)
    do i = 1, size (res_histories)
       write (u, *)
       call res_histories(i)%write (u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Re-initialize resonance history set with filter n=2"
    write (u, "(A)")

    call res_set%init (n_filter = 2)

    write (u, "(A)")  "* Add resonance histories, one at a time"
    write (u, "(A)")

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)
    call res_set%enter (res_history)
    call res_history%clear ()

    write (u, *)

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_info%init (7, 23, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)
    call res_set%enter (res_history)
    call res_history%clear ()

    write (u, *)

    call res_info%init (7, 23, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)
    call res_set%enter (res_history)
    call res_history%clear ()

    write (u, *)

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_set%freeze ()

    write (u, "(A)")
    write (u, "(A)")  "* Result"
    write (u, "(A)")

    call res_set%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonances_3"

  end subroutine resonances_3

  subroutine resonances_4 (u)
    integer, intent(in) :: u
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(model_data_t), target :: model
    type(flavor_t) :: fw, fz
    real(default) :: mw, mz, ww, wz
    type(vector4_t), dimension(3) :: p
    real(default), dimension(2) :: dist
    real(default) :: gw, factor
    integer :: i

    write (u, "(A)")  "* Test output: resonances_4"
    write (u, "(A)")  "*   Purpose: test resonance history evaluation"
    write (u, "(A)")

    write (u, "(A)")  "* Read model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* W and Z parameters"
    write (u, "(A)")

    call fw%init (24, model)
    call fz%init (23, model)
    mw = fw%get_mass ()
    ww = fw%get_width ()
    mz = fz%get_mass ()
    wz = fz%get_width ()

    write (u, "(A,1x," // FMF_12 // ")")  "mW  =", mw
    write (u, "(A,1x," // FMF_12 // ")")  "wW  =", ww
    write (u, "(A,1x," // FMF_12 // ")")  "mZ  =", mz
    write (u, "(A,1x," // FMF_12 // ")")  "wZ  =", wz

    write (u, "(A)")
    write (u, "(A)")  "* Gaussian width parameter"
    write (u, "(A)")

    gw = 2
    write (u, "(A,1x," // FMF_12 // ")")  "gw  =", gw

    write (u, "(A)")
    write (u, "(A)")  "* Setup resonance histories"
    write (u, "(A)")

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)

    call res_info%init (7, 23, model, 5)
    call res_history%add_resonance (res_info)

    call res_history%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Setup zero momenta"
    write (u, "(A)")

    do i = 1, 3
       call p(i)%write (u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate distances from resonances"
    write (u, "(A)")

    call res_history%evaluate_distances (p, dist)
    write (u, "(A,1x," // FMF_12 // ")")  "distance (W) =", dist(1)
    write (u, "(A,1x," // FMF_12 // ")")  "m/w (W)      =", mw / ww
    write (u, "(A,1x," // FMF_12 // ")")  "distance (Z) =", dist(2)
    write (u, "(A,1x," // FMF_12 // ")")  "m/w (Z)      =", mz / wz

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate Gaussian turnoff factor"
    write (u, "(A)")

    factor = res_history%evaluate_gaussian (p, gw)
    write (u, "(A,1x," // FMF_12 // ")")  "gaussian fac =", factor

    write (u, "(A)")
    write (u, "(A)")  "* Set momenta on W peak"
    write (u, "(A)")

    p(1) = vector4_at_rest (mw/2)
    p(2) = vector4_at_rest (mw/2)
    do i = 1, 3
       call p(i)%write (u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate distances from resonances"
    write (u, "(A)")

    call res_history%evaluate_distances (p, dist)
    write (u, "(A,1x," // FMF_12 // ")")  "distance (W) =", dist(1)
    write (u, "(A,1x," // FMF_12 // ")")  "distance (Z) =", dist(2)
    write (u, "(A,1x," // FMF_12 // ")")  "expected     =", &
         abs (mz**2 - mw**2) / (mz*wz)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate Gaussian turnoff factor"
    write (u, "(A)")

    factor = res_history%evaluate_gaussian (p, gw)
    write (u, "(A,1x," // FMF_12 // ")")  "gaussian fac =", factor
    write (u, "(A,1x," // FMF_12 // ")")  "expected     =", &
         exp (- (abs (mz**2 - mw**2) / (mz*wz))**2 / (gw * wz)**2)

    write (u, "(A)")
    write (u, "(A)")  "* Set momenta on both peaks"
    write (u, "(A)")

    p(3) = vector4_at_rest (mz - mw)
    do i = 1, 3
       call p(i)%write (u)
    end do

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate distances from resonances"
    write (u, "(A)")

    call res_history%evaluate_distances (p, dist)
    write (u, "(A,1x," // FMF_12 // ")")  "distance (W) =", dist(1)
    write (u, "(A,1x," // FMF_12 // ")")  "distance (Z) =", dist(2)

    write (u, "(A)")
    write (u, "(A)")  "* Evaluate Gaussian turnoff factor"
    write (u, "(A)")

    factor = res_history%evaluate_gaussian (p, gw)
    write (u, "(A,1x," // FMF_12 // ")")  "gaussian fac =", factor

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonances_4"

  end subroutine resonances_4

  subroutine resonances_5 (u)
    integer, intent(in) :: u
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(resonance_history_set_t) :: res_set
    type(model_data_t), target :: model
    type(flavor_t) :: fw, fz
    real(default) :: mw, mz, ww, wz
    real(default) :: on_shell_limit
    integer, dimension(:), allocatable :: on_shell
    type(vector4_t), dimension(4) :: p

    write (u, "(A)")  "* Test output: resonances_5"
    write (u, "(A)")  "*   Purpose: resonance history on-shell test"
    write (u, "(A)")

    write (u, "(A)")  "* Read model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* W and Z parameters"
    write (u, "(A)")

    call fw%init (24, model)
    call fz%init (23, model)
    mw = fw%get_mass ()
    ww = fw%get_width ()
    mz = fz%get_mass ()
    wz = fz%get_width ()

    write (u, "(A,1x," // FMF_12 // ")")  "mW  =", mw
    write (u, "(A,1x," // FMF_12 // ")")  "wW  =", ww
    write (u, "(A,1x," // FMF_12 // ")")  "mZ  =", mz
    write (u, "(A,1x," // FMF_12 // ")")  "wZ  =", wz

    write (u, "(A)")
    write (u, "(A)")  "* On-shell parameter: distance as multiple of width"
    write (u, "(A)")

    on_shell_limit = 3
    write (u, "(A,1x," // FMF_12 // ")")  "on-shell limit =", on_shell_limit


    write (u, "(A)")
    write (u, "(A)")  "* Setup resonance history set"
    write (u, "(A)")

    call res_set%init ()

    call res_info%init (3, -24, model, 6)
    call res_history%add_resonance (res_info)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_info%init (12, 24, model, 6)
    call res_history%add_resonance (res_info)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_info%init (15, 23, model, 6)
    call res_history%add_resonance (res_info)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_info%init (3, -24, model, 6)
    call res_history%add_resonance (res_info)
    call res_info%init (15, 23, model, 6)
    call res_history%add_resonance (res_info)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_info%init (12, 24, model, 6)
    call res_history%add_resonance (res_info)
    call res_info%init (15, 23, model, 6)
    call res_history%add_resonance (res_info)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_set%freeze ()
    call res_set%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Setup zero momenta"
    write (u, "(A)")

    call write_momenta (p)

    call res_set%determine_on_shell_histories (p, on_shell_limit, on_shell)
    call write_on_shell_histories (on_shell)

    write (u, "(A)")
    write (u, "(A)")  "* Setup momenta near W- resonance (2 widths off)"
    write (u, "(A)")

    p(1) = vector4_at_rest (82.5_default)
    call write_momenta (p)

    call res_set%determine_on_shell_histories (p, on_shell_limit, on_shell)
    call write_on_shell_histories (on_shell)

    write (u, "(A)")
    write (u, "(A)")  "* Setup momenta near W- resonance (4 widths off)"
    write (u, "(A)")

    p(1) = vector4_at_rest (84.5_default)
    call write_momenta (p)

    call res_set%determine_on_shell_histories (p, on_shell_limit, on_shell)
    call write_on_shell_histories (on_shell)

    write (u, "(A)")
    write (u, "(A)")  "* Setup momenta near Z resonance"
    write (u, "(A)")

    p(1) = vector4_at_rest (45._default)
    p(3) = vector4_at_rest (45._default)
    call write_momenta (p)

    call res_set%determine_on_shell_histories (p, on_shell_limit, on_shell)
    call write_on_shell_histories (on_shell)

    write (u, "(A)")
    write (u, "(A)")  "* Setup momenta near W- and W+ resonances"
    write (u, "(A)")

    p(1) = vector4_at_rest (40._default)
    p(2) = vector4_at_rest (40._default)
    p(3) = vector4_at_rest (40._default)
    p(4) = vector4_at_rest (40._default)
    call write_momenta (p)

    call res_set%determine_on_shell_histories (p, on_shell_limit, on_shell)
    call write_on_shell_histories (on_shell)

    write (u, "(A)")
    write (u, "(A)")  "* Setup momenta near W- and Z resonances, &
         &shadowing single resonances"
    write (u, "(A)")

    p(1) = vector4_at_rest (40._default)
    p(2) = vector4_at_rest (40._default)
    p(3) = vector4_at_rest (10._default)
    p(4) = vector4_at_rest ( 0._default)
    call write_momenta (p)

    call res_set%determine_on_shell_histories (p, on_shell_limit, on_shell)
    call write_on_shell_histories (on_shell)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonances_5"

  contains

    subroutine write_momenta (p)
      type(vector4_t), dimension(:), intent(in) :: p
      integer :: i
      do i = 1, size (p)
         call p(i)%write (u)
      end do
    end subroutine write_momenta

    subroutine write_on_shell_histories (on_shell)
      integer, dimension(:), intent(in) :: on_shell
      integer :: i
      write (u, *)
      write (u, "(A)", advance="no")  "on-shell = ("
      do i = 1, size (on_shell)
         if (i > 1)  write (u, "(',')", advance="no")
         write (u, "(I0)", advance="no")  on_shell(i)
      end do
      write (u, "(')')")
    end subroutine write_on_shell_histories

  end subroutine resonances_5

  subroutine resonances_6 (u)
    integer, intent(in) :: u
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(resonance_tree_t) :: res_tree
    type(model_data_t), target :: model

    write (u, "(A)")  "* Test output: resonances_6"
    write (u, "(A)")  "*   Purpose: retrieve resonance histories as trees"
    write (u, "(A)")

    write (u, "(A)")  "* Read model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Empty resonance history"
    write (u, "(A)")

    call res_history%write (u)

    write (u, "(A)")
    call res_history%to_tree (res_tree)
    call res_tree%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Single resonance"
    write (u, "(A)")

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)

    write (u, "(A)")
    call res_history%to_tree (res_tree)
    call res_tree%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Nested resonances"
    write (u, "(A)")

    call res_info%init (7, 23, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%write (u)

    write (u, "(A)")
    call res_history%to_tree (res_tree)
    call res_tree%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Disjunct resonances"
    write (u, "(A)")

    call res_history%clear ()

    call res_info%init (5, 24, model, 7)
    call res_history%add_resonance (res_info)

    call res_info%init (7, 6, model, 7)
    call res_history%add_resonance (res_info)

    call res_info%init (80, -24, model, 7)
    call res_history%add_resonance (res_info)

    call res_info%init (112, -6, model, 7)
    call res_history%add_resonance (res_info)

    call res_history%write (u)

    write (u, "(A)")
    call res_history%to_tree (res_tree)
    call res_tree%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonances_6"

  end subroutine resonances_6

  subroutine resonances_7 (u)
    integer, intent(in) :: u
    type(resonance_info_t) :: res_info
    type(resonance_history_t) :: res_history
    type(resonance_tree_t) :: res_tree
    type(resonance_history_set_t) :: res_set
    type(model_data_t), target :: model
    type(flavor_t) :: flv

    write (u, "(A)")  "* Test output: resonances_7"
    write (u, "(A)")  "*   Purpose: test tree format"
    write (u, "(A)")

    write (u, "(A)")  "* Read model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Initialize, fill and freeze resonance history set"
    write (u, "(A)")

    call res_set%init (initial_size = 2)

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_history%clear ()

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_info%init (7, 23, model, 5)
    call res_history%add_resonance (res_info)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_info%init (7, 23, model, 5)
    call res_history%add_resonance (res_info)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_info%init (3, -24, model, 5)
    call res_history%add_resonance (res_info)
    call res_info%init (7, 25, model, 5)
    call res_history%add_resonance (res_info)
    call res_set%enter (res_history)
    call res_history%clear ()

    call res_set%freeze ()

    call res_set%write (u, show_trees = .true.)

    write (u, "(A)")
    write (u, "(A)")  "* Extract tree #1"
    write (u, "(A)")

    call res_set%get_tree (1, res_tree)
    call res_tree%write (u)

    write (u, *)
    write (u, "(1x,A,1x,I0)")  "n_resonances =", res_tree%get_n_resonances ()

    write (u, *)
    write (u, "(1x,A,1x)", advance="no")  "flv(r1) ="
    flv = res_tree%get_flv (1)
    call flv%write (u)
    write (u, *)
    write (u, "(1x,A,1x)", advance="no")  "flv(r2) ="
    flv = res_tree%get_flv (2)
    call flv%write (u)
    write (u, *)

    write (u, *)
    write (u, "(1x,A)")  "[offset = 2, 4]"
    write (u, "(1x,A,9(1x,I0))")  "children(r1) =", &
         res_tree%get_children(1, 2, 4)
    write (u, "(1x,A,9(1x,I0))")  "children(r2) =", &
         res_tree%get_children(2, 2, 4)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: resonances_7"

  end subroutine resonances_7


end module resonances_uti
