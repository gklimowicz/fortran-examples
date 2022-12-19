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

module auto_components_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use pdg_arrays
  use model_data
  use model_testbed, only: prepare_model, cleanup_model

  use auto_components

  implicit none
  private

  public :: auto_components_1
  public :: auto_components_2
  public :: auto_components_3

contains

  subroutine auto_components_1 (u)
    integer, intent(in) :: u
    class(model_data_t), pointer :: model
    type(field_data_t), pointer :: prt
    type(ds_table_t) :: ds_table
    type(split_constraints_t) :: constraints

    write (u, "(A)")  "* Test output: auto_components_1"
    write (u, "(A)")  "*   Purpose: determine Higgs decay table"
    write (u, *)

    write (u, "(A)")  "* Read Standard Model"

    model => null ()
    call prepare_model (model, var_str ("SM"))

    prt => model%get_field_ptr (25)

    write (u, *)
    write (u, "(A)")  "* Higgs decays n = 2"
    write (u, *)

    call constraints%init (2)
    call constraints%set (1, constrain_n_tot (2))
    call constraints%set (2, constrain_mass_sum (prt%get_mass ()))

    call ds_table%make (model, 25, constraints)
    call ds_table%write (u)
    call ds_table%final ()

    write (u, *)
    write (u, "(A)")  "* Higgs decays n = 3 (w/o radiative)"
    write (u, *)

    call constraints%init (3)
    call constraints%set (1, constrain_n_tot (3))
    call constraints%set (2, constrain_mass_sum (prt%get_mass ()))
    call constraints%set (3, constrain_radiation ())

    call ds_table%make (model, 25, constraints)
    call ds_table%write (u)
    call ds_table%final ()

    write (u, *)
    write (u, "(A)")  "* Higgs decays n = 3 (w/ radiative)"
    write (u, *)

    call constraints%init (2)
    call constraints%set (1, constrain_n_tot (3))
    call constraints%set (2, constrain_mass_sum (prt%get_mass ()))

    call ds_table%make (model, 25, constraints)
    call ds_table%write (u)
    call ds_table%final ()

    write (u, *)
    write (u, "(A)")  "* Cleanup"

    call cleanup_model (model)
    deallocate (model)

    write (u, *)
    write (u, "(A)")  "* Test output end: auto_components_1"

  end subroutine auto_components_1

  subroutine auto_components_2 (u)
    integer, intent(in) :: u
    class(model_data_t), pointer :: model
    type(pdg_list_t), dimension(:), allocatable :: pl, pl_zzh
    type(pdg_list_t) :: pl_match
    type(fs_table_t) :: fs_table
    type(split_constraints_t) :: constraints
    real(default) :: sqrts
    integer :: i

    write (u, "(A)")  "* Test output: auto_components_2"
    write (u, "(A)")  "*   Purpose: generate radiation (NLO)"
    write (u, *)


    write (u, "(A)")  "* Read Standard Model"

    model => null ()
    call prepare_model (model, var_str ("SM"))

    write (u, *)
    write (u, "(A)")  "* LO final state"
    write (u, *)

    allocate (pl (2))
    call pl(1)%init (2)
    call pl(1)%set (1, 1)
    call pl(1)%set (2, -1)
    call pl(2)%init (2)
    call pl(2)%set (1, 21)
    call pl(2)%set (2, 21)
    do i = 1, 2
       call pl(i)%write (u);  write (u, *)
    end do

    write (u, *)
    write (u, "(A)")  "* Initialize FS table"
    write (u, *)

    call constraints%init (1)
    call constraints%set (1, constrain_n_tot (3))

    call fs_table%init (model, pl, constraints)
    call fs_table%write (u)

    write (u, *)
    write (u, "(A)")  "* Generate NLO corrections, unconstrained"
    write (u, *)

    call fs_table%radiate (constraints)
    call fs_table%write (u)
    call fs_table%final ()

    write (u, *)
    write (u, "(A)")  "* Generate NLO corrections, &
         &complete but mass-constrained"
    write (u, *)

    sqrts = 50

    call constraints%init (2)
    call constraints%set (1, constrain_n_tot (3))
    call constraints%set (2, constrain_mass_sum (sqrts))

    call fs_table%init (model, pl, constraints)
    call fs_table%radiate (constraints)
    call fs_table%write (u)
    call fs_table%final ()

    write (u, *)
    write (u, "(A)")  "* Generate NLO corrections, restricted"
    write (u, *)

    call pl_match%init ([1, -1, 21])

    call constraints%init (2)
    call constraints%set (1, constrain_n_tot (3))
    call constraints%set (2, constrain_insert (pl_match))

    call fs_table%init (model, pl, constraints)
    call fs_table%radiate (constraints)
    call fs_table%write (u)
    call fs_table%final ()

    write (u, *)
    write (u, "(A)")  "* Generate NNLO corrections, restricted, with one loop"
    write (u, *)

    call pl_match%init ([1, -1, 21])

    call constraints%init (3)
    call constraints%set (1, constrain_n_tot (4))
    call constraints%set (2, constrain_n_loop (1))
    call constraints%set (3, constrain_insert (pl_match))

    call fs_table%init (model, pl, constraints)
    call fs_table%enable_loops ()
    call fs_table%radiate (constraints)
    call fs_table%write (u)
    call fs_table%final ()

    write (u, *)
    write (u, "(A)")  "* Generate NNLO corrections, restricted, with loops"
    write (u, *)

    call constraints%init (2)
    call constraints%set (1, constrain_n_tot (4))
    call constraints%set (2, constrain_insert (pl_match))

    call fs_table%init (model, pl, constraints)
    call fs_table%enable_loops ()
    call fs_table%radiate (constraints)
    call fs_table%write (u)
    call fs_table%final ()

    write (u, *)
    write (u, "(A)")  "* Generate NNLO corrections, restricted, to Z Z H, &
         &no loops"
    write (u, *)

    allocate (pl_zzh (1))
    call pl_zzh(1)%init (3)
    call pl_zzh(1)%set (1, 23)
    call pl_zzh(1)%set (2, 23)
    call pl_zzh(1)%set (3, 25)

    call constraints%init (3)
    call constraints%set (1, constrain_n_tot (5))
    call constraints%set (2, constrain_mass_sum (500._default))
    call constraints%set (3, constrain_require (pl_zzh(1)))

    call fs_table%init (model, pl_zzh, constraints)
    call fs_table%radiate (constraints)
    call fs_table%write (u)
    call fs_table%final ()

    call cleanup_model (model)
    deallocate (model)

    write (u, *)
    write (u, "(A)")  "* Test output end: auto_components_2"

  end subroutine auto_components_2

  subroutine auto_components_3 (u)
    integer, intent(in) :: u
    class(model_data_t), pointer :: model
    type(pdg_list_t), dimension(:), allocatable :: pl_in, pl_out
    type(pdg_list_t) :: pl_match, pl_beam
    type(if_table_t) :: if_table
    type(split_constraints_t) :: constraints
    real(default) :: sqrts
    integer :: i

    write (u, "(A)")  "* Test output: auto_components_3"
    write (u, "(A)")  "*   Purpose: generate radiation (NLO)"
    write (u, *)

    write (u, "(A)")  "* Read Standard Model"

    model => null ()
    call prepare_model (model, var_str ("SM"))

    write (u, *)
    write (u, "(A)")  "* LO initial state"
    write (u, *)

    allocate (pl_in (2))
    call pl_in(1)%init (2)
    call pl_in(1)%set (1, 1)
    call pl_in(1)%set (2, -1)
    call pl_in(2)%init (2)
    call pl_in(2)%set (1, -1)
    call pl_in(2)%set (2, 1)
    do i = 1, 2
       call pl_in(i)%write (u);  write (u, *)
    end do

    write (u, *)
    write (u, "(A)")  "* LO final state"
    write (u, *)

    allocate (pl_out (1))
    call pl_out(1)%init (1)
    call pl_out(1)%set (1, 23)
    call pl_out(1)%write (u);  write (u, *)

    write (u, *)
    write (u, "(A)")  "* Initialize FS table"
    write (u, *)

    call constraints%init (1)
    call constraints%set (1, constrain_n_tot (4))

    call if_table%init (model, pl_in, pl_out, constraints)
    call if_table%write (u)

    write (u, *)
    write (u, "(A)")  "* Generate NLO corrections, unconstrained"
    write (u, *)

    call if_table%radiate (constraints)
    call if_table%write (u)
    call if_table%final ()

    write (u, *)
    write (u, "(A)")  "* Generate NLO corrections, &
         &complete but mass-constrained"
    write (u, *)

    sqrts = 100
    call constraints%init (2)
    call constraints%set (1, constrain_n_tot (4))
    call constraints%set (2, constrain_mass_sum (sqrts))

    call if_table%init (model, pl_in, pl_out, constraints)
    call if_table%radiate (constraints)
    call if_table%write (u)
    call if_table%final ()

    write (u, *)
    write (u, "(A)")  "* Generate NLO corrections, &
         &mass-constrained, restricted beams"
    write (u, *)

    call pl_beam%init (3)
    call pl_beam%set (1, 1)
    call pl_beam%set (2, -1)
    call pl_beam%set (3, 21)

    call constraints%init (3)
    call constraints%set (1, constrain_n_tot (4))
    call constraints%set (2, constrain_in_state (pl_beam))
    call constraints%set (3, constrain_mass_sum (sqrts))

    call if_table%init (model, pl_in, pl_out, constraints)
    call if_table%radiate (constraints)
    call if_table%write (u)
    call if_table%final ()

    write (u, *)
    write (u, "(A)")  "* Generate NLO corrections, restricted"
    write (u, *)

    call pl_match%init ([1, -1, 21])

    call constraints%init (4)
    call constraints%set (1, constrain_n_tot (4))
    call constraints%set (2, constrain_in_state (pl_beam))
    call constraints%set (3, constrain_mass_sum (sqrts))
    call constraints%set (4, constrain_insert (pl_match))

    call if_table%init (model, pl_in, pl_out, constraints)
    call if_table%radiate (constraints)
    call if_table%write (u)
    call if_table%final ()

    write (u, *)
    write (u, "(A)")  "* Generate NNLO corrections, restricted, Z preserved, &
         &with loops"
    write (u, *)

    call constraints%init (5)
    call constraints%set (1, constrain_n_tot (5))
    call constraints%set (2, constrain_in_state (pl_beam))
    call constraints%set (3, constrain_mass_sum (sqrts))
    call constraints%set (4, constrain_insert (pl_match))
    call constraints%set (5, constrain_require (pl_out(1)))

    call if_table%init (model, pl_in, pl_out, constraints)
    call if_table%enable_loops ()
    call if_table%radiate (constraints)
    call if_table%write (u)
    call if_table%final ()

    call cleanup_model (model)
    deallocate (model)

    write (u, *)
    write (u, "(A)")  "* Test output end: auto_components_3"

  end subroutine auto_components_3


end module auto_components_uti
