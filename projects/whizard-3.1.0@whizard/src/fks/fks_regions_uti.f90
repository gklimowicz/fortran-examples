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

module fks_regions_uti

  use iso_varying_string, string_t => varying_string
  use format_utils, only: write_separator
  use os_interface
  use models

  use fks_regions

  implicit none
  private

  public :: fks_regions_1
  public :: fks_regions_2
  public :: fks_regions_3
  public :: fks_regions_4
  public :: fks_regions_5
  public :: fks_regions_6
  public :: fks_regions_7
  public :: fks_regions_8

contains

  subroutine fks_regions_1 (u)
    integer, intent(in) :: u
    type(flv_structure_t) :: flv_born, flv_real
    type(model_t), pointer :: test_model => null ()
    write (u, "(A)") "* Test output: fks_regions_1"
    write (u, "(A)") "* Purpose: Test utilities of flavor structure manipulation"
    write (u, "(A)")

    call create_test_model (var_str ("SM"), test_model)

    flv_born = [11, -11, 2, -2]
    flv_real = [11, -11, 2, -2, 21]
    flv_born%n_in = 2; flv_real%n_in = 2
    write (u, "(A)") "* Valid splittings of ee -> uu"
    write (u, "(A)") "Born Flavors: "
    call flv_born%write (u)
    write (u, "(A)") "Real Flavors: "
    call flv_real%write (u)
    write (u, "(A,L1)") "3, 4 (2, -2) : ", flv_real%valid_pair (3, 4, flv_born, test_model)
    write (u, "(A,L1)") "4, 3 (-2, 2) : ", flv_real%valid_pair (4, 3, flv_born, test_model)
    write (u, "(A,L1)") "3, 5 (2, 21) : ", flv_real%valid_pair (3, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 3 (21, 2) : ", flv_real%valid_pair (5, 3, flv_born, test_model)
    write (u, "(A,L1)") "4, 5 (-2, 21): ", flv_real%valid_pair (4, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 4 (21, -2): ", flv_real%valid_pair (5, 4, flv_born, test_model)
    call write_separator (u)

    call flv_born%final ()
    call flv_real%final ()

    flv_born = [2, -2, 11, -11]
    flv_real = [2, -2, 11, -11, 21]
    flv_born%n_in = 2; flv_real%n_in = 2
    write (u, "(A)") "* Valid splittings of uu -> ee"
    write (u, "(A)") "Born Flavors: "
    call flv_born%write (u)
    write (u, "(A)") "Real Flavors: "
    call flv_real%write (u)
    write (u, "(A,L1)") "1, 2 (2, -2) : " , flv_real%valid_pair (1, 2, flv_born, test_model)
    write (u, "(A,L1)") "2, 1 (-2, 2) : " , flv_real%valid_pair (2, 1, flv_born, test_model)
    write (u, "(A,L1)") "5, 2 (21, -2): " , flv_real%valid_pair (5, 2, flv_born, test_model)
    write (u, "(A,L1)") "2, 5 (-2, 21): " , flv_real%valid_pair (2, 5, flv_born, test_model)
    write (u, "(A,L1)") "1, 5 (21, 2) : " , flv_real%valid_pair (5, 1, flv_born, test_model)
    write (u, "(A,L1)") "5, 1 (2, 21) : " , flv_real%valid_pair (1, 5, flv_born, test_model)
    call flv_real%final ()
    flv_real = [21, -2, 11, -11, -2]
    flv_real%n_in = 2
    write (u, "(A)") "Real Flavors: "
    call flv_real%write (u)
    write (u, "(A,L1)") "1, 2 (21, -2): " , flv_real%valid_pair (1, 2, flv_born, test_model)
    write (u, "(A,L1)") "2, 1 (-2, 21): " , flv_real%valid_pair (2, 1, flv_born, test_model)
    write (u, "(A,L1)") "5, 2 (-2, -2): " , flv_real%valid_pair (5, 2, flv_born, test_model)
    write (u, "(A,L1)") "2, 5 (-2, -2): " , flv_real%valid_pair (2, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 1 (-2, 21): " , flv_real%valid_pair (5, 1, flv_born, test_model)
    write (u, "(A,L1)") "1, 5 (21, -2): " , flv_real%valid_pair (1, 5, flv_born, test_model)
    call flv_real%final ()
    flv_real = [2, 21, 11, -11, 2]
    flv_real%n_in = 2
    write (u, "(A)") "Real Flavors: "
    call flv_real%write (u)
    write (u, "(A,L1)") "1, 2 (2, 21) : " , flv_real%valid_pair (1, 2, flv_born, test_model)
    write (u, "(A,L1)") "2, 1 (21, 2) : " , flv_real%valid_pair (2, 1, flv_born, test_model)
    write (u, "(A,L1)") "5, 2 (2, 21) : " , flv_real%valid_pair (5, 2, flv_born, test_model)
    write (u, "(A,L1)") "2, 5 (21, 2) : " , flv_real%valid_pair (2, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 1 (2, 2)  : " , flv_real%valid_pair (5, 1, flv_born, test_model)
    write (u, "(A,L1)") "1, 5 (2, 2)  : " , flv_real%valid_pair (1, 5, flv_born, test_model)
    call write_separator (u)

    call flv_born%final ()
    call flv_real%final ()

    flv_born = [11, -11, 2, -2, 21]
    flv_real = [11, -11, 2, -2, 21, 21]
    flv_born%n_in = 2; flv_real%n_in = 2
    write (u, "(A)") "* Valid splittings of ee -> uug"
    write (u, "(A)") "Born Flavors: "
    call flv_born%write (u)
    write (u, "(A)") "Real Flavors: "
    call flv_real%write (u)
    write (u, "(A,L1)") "3, 4 (2, -2) : " , flv_real%valid_pair (3, 4, flv_born, test_model)
    write (u, "(A,L1)") "4, 3 (-2, 2) : " , flv_real%valid_pair (4, 3, flv_born, test_model)
    write (u, "(A,L1)") "3, 5 (2, 21) : " , flv_real%valid_pair (3, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 3 (21, 2) : " , flv_real%valid_pair (5, 3, flv_born, test_model)
    write (u, "(A,L1)") "4, 5 (-2, 21): " , flv_real%valid_pair (4, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 4 (21, -2): " , flv_real%valid_pair (5, 4, flv_born, test_model)
    write (u, "(A,L1)") "3, 6 (2, 21) : " , flv_real%valid_pair (3, 6, flv_born, test_model)
    write (u, "(A,L1)") "6, 3 (21, 2) : " , flv_real%valid_pair (6, 3, flv_born, test_model)
    write (u, "(A,L1)") "4, 6 (-2, 21): " , flv_real%valid_pair (4, 6, flv_born, test_model)
    write (u, "(A,L1)") "6, 4 (21, -2): " , flv_real%valid_pair (6, 4, flv_born, test_model)
    write (u, "(A,L1)") "5, 6 (21, 21): " , flv_real%valid_pair (5, 6, flv_born, test_model)
    write (u, "(A,L1)") "6, 5 (21, 21): " , flv_real%valid_pair (6, 5, flv_born, test_model)
    call flv_real%final ()
    flv_real = [11, -11, 2, -2, 1, -1]
    flv_real%n_in = 2
    write (u, "(A)") "Real Flavors (exemplary g -> dd splitting): "
    call flv_real%write (u)
    write (u, "(A,L1)") "3, 4 (2, -2) : " , flv_real%valid_pair (3, 4, flv_born, test_model)
    write (u, "(A,L1)") "4, 3 (-2, 2) : " , flv_real%valid_pair (4, 3, flv_born, test_model)
    write (u, "(A,L1)") "3, 5 (2, 1)  : " , flv_real%valid_pair (3, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 3 (1, 2)  : " , flv_real%valid_pair (5, 3, flv_born, test_model)
    write (u, "(A,L1)") "4, 5 (-2, 1) : " , flv_real%valid_pair (4, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 4 (1, -2) : " , flv_real%valid_pair (5, 4, flv_born, test_model)
    write (u, "(A,L1)") "3, 6 (2, -1) : " , flv_real%valid_pair (3, 6, flv_born, test_model)
    write (u, "(A,L1)") "6, 3 (-1, 2) : " , flv_real%valid_pair (6, 3, flv_born, test_model)
    write (u, "(A,L1)") "4, 6 (-2, -1): " , flv_real%valid_pair (4, 6, flv_born, test_model)
    write (u, "(A,L1)") "6, 4 (-1, -2): " , flv_real%valid_pair (6, 4, flv_born, test_model)
    write (u, "(A,L1)") "5, 6 (1, -1) : " , flv_real%valid_pair (5, 6, flv_born, test_model)
    write (u, "(A,L1)") "6, 5 (-1, 1) : " , flv_real%valid_pair (6, 5, flv_born, test_model)
    call write_separator (u)

    call flv_born%final ()
    call flv_real%final ()

    flv_born = [6, -5, 2, -1 ]
    flv_real = [6, -5, 2, -1, 21]
    flv_born%n_in = 1; flv_real%n_in = 1
    write (u, "(A)") "* Valid splittings of t -> b u d~"
    write (u, "(A)") "Born Flavors: "
    call flv_born%write (u)
    write (u, "(A)") "Real Flavors: "
    call flv_real%write (u)
    write (u, "(A,L1)") "1, 2 (6, -5) : " , flv_real%valid_pair (1, 2, flv_born, test_model)
    write (u, "(A,L1)") "1, 3 (6, 2)  : " , flv_real%valid_pair (1, 3, flv_born, test_model)
    write (u, "(A,L1)") "1, 4 (6, -1) : " , flv_real%valid_pair (1, 4, flv_born, test_model)
    write (u, "(A,L1)") "2, 1 (-5, 6) : " , flv_real%valid_pair (2, 1, flv_born, test_model)
    write (u, "(A,L1)") "3, 1 (2, 6)  : " , flv_real%valid_pair (3, 1, flv_born, test_model)
    write (u, "(A,L1)") "4, 1 (-1, 6) : " , flv_real%valid_pair (4, 1, flv_born, test_model)
    write (u, "(A,L1)") "2, 3 (-5, 2) : " , flv_real%valid_pair (2, 3, flv_born, test_model)
    write (u, "(A,L1)") "2, 4 (-5, -1): " , flv_real%valid_pair (2, 4, flv_born, test_model)
    write (u, "(A,L1)") "3, 2 (2, -5) : " , flv_real%valid_pair (3, 2, flv_born, test_model)
    write (u, "(A,L1)") "4, 2 (-1, -5): " , flv_real%valid_pair (4, 2, flv_born, test_model)
    write (u, "(A,L1)") "3, 4 (2, -1) : " , flv_real%valid_pair (3, 4, flv_born, test_model)
    write (u, "(A,L1)") "4, 3 (-1, 2) : " , flv_real%valid_pair (4, 3, flv_born, test_model)
    write (u, "(A,L1)") "1, 5 (6, 21) : " , flv_real%valid_pair (1, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 1 (21, 6) : " , flv_real%valid_pair (5, 1, flv_born, test_model)
    write (u, "(A,L1)") "2, 5 (-5, 21): " , flv_real%valid_pair (2, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 2 (21, 5) : " , flv_real%valid_pair (5, 2, flv_born, test_model)
    write (u, "(A,L1)") "3, 5 (2, 21) : " , flv_real%valid_pair (3, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 3 (21, 2) : " , flv_real%valid_pair (5, 3, flv_born, test_model)
    write (u, "(A,L1)") "4, 5 (-1, 21): " , flv_real%valid_pair (4, 5, flv_born, test_model)
    write (u, "(A,L1)") "5, 4 (21, -1): " , flv_real%valid_pair (5, 4, flv_born, test_model)

    call flv_born%final ()
    call flv_real%final ()

  end subroutine fks_regions_1

  subroutine fks_regions_2 (u)
    integer, intent(in) :: u
    integer :: n_flv_born, n_flv_real
    integer :: n_legs_born, n_legs_real
    integer :: n_in
    integer, dimension(:,:), allocatable :: flv_born, flv_real
    type(region_data_t) :: reg_data
    write (u, "(A)") "* Test output: fks_regions_2"
    write (u, "(A)") "* Create singular regions for processes with up to four singular regions"
    write (u, "(A)") "* ee -> qq with QCD corrections"
    write (u, "(A)")

    n_flv_born = 1; n_flv_real = 1
    n_legs_born = 4; n_legs_real = 5
    n_in = 2

    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))
    flv_born (:, 1) = [11, -11, 2, -2]
    flv_real (:, 1) = [11, -11, 2, -2, 21]
    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 2, 0)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()
    call write_separator (u)

    write (u, "(A)") "* ee -> qq with EW corrections"
    write (u, "(A)")

    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))
    flv_born (:, 1) = [11, -11, 2, -2]
    flv_real (:, 1) = [11, -11, 2, -2, 22]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("EW"), 2, 0)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()
    call write_separator (u)

    write (u, "(A)") "* ee -> tt"
    write (u, "(A)")
    write (u, "(A)") "* This process has four singular regions because they are not equivalent."

    n_flv_born = 1; n_flv_real = 1
    n_legs_born = 6; n_legs_real = 7
    n_in = 2

    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))
    flv_born (:, 1) = [11, -11, 6, -6, 6, -6]
    flv_real (:, 1) = [11, -11, 6, -6, 6, -6, 21]
    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 4, 0)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()
  end subroutine fks_regions_2

  subroutine fks_regions_3 (u)
    integer, intent(in) :: u
    integer :: n_flv_born, n_flv_real
    integer :: n_legs_born, n_legs_real
    integer :: n_in, i, j
    integer, dimension(:,:), allocatable :: flv_born, flv_real
    type(region_data_t) :: reg_data
    write (u, "(A)") "* Test output: fks_regions_3"
    write (u, "(A)") "* Create singular regions for processes with three final-state particles"

    write (u, "(A)") "* ee -> qqg"
    write (u, "(A)")
    n_flv_born = 1; n_flv_real = 2
    n_legs_born = 5; n_legs_real = 6
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))

    flv_born (:, 1) = [11, -11, 2, -2, 21]
    flv_real (:, 1) = [11, -11, 2, -2, 21, 21]
    flv_real (:, 2) = [11, -11, 2, -2, 1, -1]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 2, 1)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()
    call write_separator (u)

    write (u, "(A)") "* ee -> qqA"
    write (u, "(A)")
    n_flv_born = 1; n_flv_real = 2
    n_legs_born = 5; n_legs_real = 6
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))

    flv_born (:, 1) = [11, -11, 2, -2, 22]
    flv_real (:, 1) = [11, -11, 2, -2, 22, 22]
    flv_real (:, 2) = [11, -11, 2, -2, 11, -11]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("EW"), 3, 0)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()
    call write_separator (u)

    write (u, "(A)") "* ee -> jet jet jet"
    write (u, "(A)") "* with jet = u:U:d:D:s:S:c:C:b:B:gl"
    write (u, "(A)")
    n_flv_born = 5; n_flv_real = 22
    n_legs_born = 5; n_legs_real = 6
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))

    flv_born (:, 1)  = [11, -11, -4, 4, 21]
    flv_born (:, 2)  = [11, -11, -2, 2, 21]
    flv_born (:, 3)  = [11, -11, -5, 5, 21]
    flv_born (:, 4)  = [11, -11, -3, 3, 21]
    flv_born (:, 5)  = [11, -11, -1, 1, 21]
    flv_real (:, 1)  = [11, -11, -4, -4, 4, 4]
    flv_real (:, 2)  = [11, -11, -4, -2, 2, 4]
    flv_real (:, 3)  = [11, -11, -4, 4, 21, 21]
    flv_real (:, 4)  = [11, -11, -4, -5, 4, 5]
    flv_real (:, 5)  = [11, -11, -4, -3, 4, 3]
    flv_real (:, 6)  = [11, -11, -4, -1, 2, 3]
    flv_real (:, 7)  = [11, -11, -4, -1, 4, 1]
    flv_real (:, 8)  = [11, -11, -2, -2, 2, 2]
    flv_real (:, 9)  = [11, -11, -2, 2, 21, 21]
    flv_real (:, 10) = [11, -11, -2, -5, 2, 5]
    flv_real (:, 11) = [11, -11, -2, -3, 2, 3]
    flv_real (:, 12) = [11, -11, -2, -3, 4, 1]
    flv_real (:, 13) = [11, -11, -2, -1, 2, 1]
    flv_real (:, 14) = [11, -11, -5, -5, 5, 5]
    flv_real (:, 15) = [11, -11, -5, -3, 3, 5]
    flv_real (:, 16) = [11, -11, -5, -1, 1, 5]
    flv_real (:, 17) = [11, -11, -5, 5, 21, 21]
    flv_real (:, 18) = [11, -11, -3, -3, 3, 3]
    flv_real (:, 19) = [11, -11, -3, -1, 1, 3]
    flv_real (:, 20) = [11, -11, -3, 3, 21, 21]
    flv_real (:, 21) = [11, -11, -1, -1, 1, 1]
    flv_real (:, 22) = [11, -11, -1, 1, 21, 21]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 2, 1)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()
    call write_separator (u)

    write (u, "(A)") "* ee -> L L A"
    write (u, "(A)") "* with L = e2:E2:e3:E3"
    write (u, "(A)")
    n_flv_born = 2; n_flv_real = 6
    n_legs_born = 5; n_legs_real = 6
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))

    flv_born (:, 1) = [11, -11, -15, 15, 22]
    flv_born (:, 2) = [11, -11, -13, 13, 22]

    flv_real (:, 1) = [11, -11, -15, -15, 15, 15]
    flv_real (:, 2) = [11, -11, -15, -13, 13, 13]
    flv_real (:, 3) = [11, -11, -13, -15, 13, 15]
    flv_real (:, 4) = [11, -11, -15, 15, 22, 22]
    flv_real (:, 5) = [11, -11, -13, -13, 13, 13]
    flv_real (:, 6) = [11, -11, -13, 13, 22, 22]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("EW"), 3, 0)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()

  end subroutine fks_regions_3

  subroutine fks_regions_4 (u)
    integer, intent(in) :: u
    integer :: n_flv_born, n_flv_real
    integer :: n_legs_born, n_legs_real
    integer :: n_in
    integer, dimension(:,:), allocatable :: flv_born, flv_real
    type(region_data_t) :: reg_data
    write (u, "(A)") "* Test output: fks_regions_4"
    write (u, "(A)") "* Create singular regions for processes with four final-state particles"
    write (u, "(A)") "* ee -> 4 jet"
    write (u, "(A)") "* with jet = u:U:d:D:s:S:c:C:b:B:gl"
    write (u, "(A)")
    n_flv_born = 22; n_flv_real = 22
    n_legs_born = 6; n_legs_real = 7
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))

    flv_born (:, 1)   = [11, -11, -4, -4, 4, 4]
    flv_born (:, 2)   = [11, -11, -4, -2, 2, 4]
    flv_born (:, 3)   = [11, -11, -4, 4, 21, 21]
    flv_born (:, 4)   = [11, -11, -4, -5, 4, 5]
    flv_born (:, 5)   = [11, -11, -4, -3, 4, 3]
    flv_born (:, 6)   = [11, -11, -4, -1, 2, 3]
    flv_born (:, 7)   = [11, -11, -4, -1, 4, 1]
    flv_born (:, 8)   = [11, -11, -2, -2, 2, 2]
    flv_born (:, 9)   = [11, -11, -2, 2, 21, 21]
    flv_born (:, 10)  = [11, -11, -2, -5, 2, 5]
    flv_born (:, 11)  = [11, -11, -2, -3, 2, 3]
    flv_born (:, 12)  = [11, -11, -2, -3, 4, 1]
    flv_born (:, 13)  = [11, -11, -2, -1, 2, 1]
    flv_born (:, 14)  = [11, -11, -5, -5, 5, 5]
    flv_born (:, 15)  = [11, -11, -5, -3, 3, 5]
    flv_born (:, 16)  = [11, -11, -5, -1, 1, 5]
    flv_born (:, 17)  = [11, -11, -5, 5, 21, 21]
    flv_born (:, 18)  = [11, -11, -3, -3, 3, 3]
    flv_born (:, 19)  = [11, -11, -3, -1, 1, 3]
    flv_born (:, 20)  = [11, -11, -3, -3, 21, 21]
    flv_born (:, 21)  = [11, -11, -1, -1, 1, 1]
    flv_born (:, 22)  = [11, -11, -1, 1, 21, 21]
    flv_real (:, 1)   = [11, -11, -4, -4, 4, 4, 21]
    flv_real (:, 2)   = [11, -11, -4, -2, 2, 4, 21]
    flv_real (:, 3)   = [11, -11, -4, 4, 21, 21, 21]
    flv_real (:, 4)   = [11, -11, -4, -5, 4, 5, 21]
    flv_real (:, 5)   = [11, -11, -4, -3, 4, 3, 21]
    flv_real (:, 6)   = [11, -11, -4, -1, 2, 3, 21]
    flv_real (:, 7)   = [11, -11, -4, -1, 4, 1, 21]
    flv_real (:, 8)   = [11, -11, -2, -2, 2, 2, 21]
    flv_real (:, 9)   = [11, -11, -2, 2, 21, 21, 21]
    flv_real (:, 10)  = [11, -11, -2, -5, 2, 5, 21]
    flv_real (:, 11)  = [11, -11, -2, -3, 2, 3, 21]
    flv_real (:, 12)  = [11, -11, -2, -3, 4, 1, 21]
    flv_real (:, 13)  = [11, -11, -2, -1, 2, 1, 21]
    flv_real (:, 14)  = [11, -11, -5, -5, 5, 5, 21]
    flv_real (:, 15)  = [11, -11, -5, -3, 3, 5, 21]
    flv_real (:, 16)  = [11, -11, -5, -1, 1, 5, 21]
    flv_real (:, 17)  = [11, -11, -5, 5, 21, 21, 21]
    flv_real (:, 18)  = [11, -11, -3, -3, 3, 3, 21]
    flv_real (:, 19)  = [11, -11, -3, -1, 1, 3, 21]
    flv_real (:, 20)  = [11, -11, -3, 3, 21, 21, 21]
    flv_real (:, 21)  = [11, -11, -1, -1, 1, 1, 21]
    flv_real (:, 22)  = [11, -11, -1, 1, 21, 21, 21]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 2, 2)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()
    call write_separator (u)

    write (u, "(A)") "* ee -> bbmumu with QCD corrections"
    write (u, "(A)")
    n_flv_born = 1; n_flv_real = 1
    n_legs_born = 6; n_legs_real = 7
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))

    flv_born (:, 1)   = [11, -11, -5, 5, -13, 13]
    flv_real (:, 1)   = [11, -11, -5, 5, -13, 13, 21]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 4, 0)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()
    call write_separator (u)

    write (u, "(A)") "* ee -> bbmumu with EW corrections"
    write (u, "(A)")
    n_flv_born = 1; n_flv_real = 1
    n_legs_born = 6; n_legs_real = 7
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))

    flv_born (:, 1)   = [11, -11, -5, 5, -13, 13]
    flv_real (:, 1)   = [11, -11, -5, 5, -13, 13, 22]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 4, 0)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()

  end subroutine fks_regions_4

  subroutine fks_regions_5 (u)
    integer, intent(in) :: u
    integer :: n_flv_born, n_flv_real
    integer :: n_legs_born, n_legs_real
    integer :: n_in
    integer, dimension(:,:), allocatable :: flv_born, flv_real
    type(region_data_t) :: reg_data
    write (u, "(A)") "* Test output: fks_regions_5"
    write (u, "(A)") "* Create singular regions for processes with five final-state particles"
    write (u, "(A)") "* ee -> 5 jet"
    write (u, "(A)") "* with jet = u:U:d:D:s:S:c:C:b:B:gl"
    write (u, "(A)")
    n_flv_born = 22; n_flv_real = 67
    n_legs_born = 7; n_legs_real = 8
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))


    flv_born (:,1)  = [11,-11,-4,-4,4,4,21]
    flv_born (:,2)  = [11,-11,-4,-2,2,4,21]
    flv_born (:,3)  = [11,-11,-4,4,21,21,21]
    flv_born (:,4)  = [11,-11,-4,-5,4,5,21]
    flv_born (:,5)  = [11,-11,-4,-3,4,3,21]
    flv_born (:,6)  = [11,-11,-4,-1,2,3,21]
    flv_born (:,7)  = [11,-11,-4,-1,4,1,21]
    flv_born (:,8)  = [11,-11,-2,-2,2,2,21]
    flv_born (:,9)  = [11,-11,-2,2,21,21,21]
    flv_born (:,10) = [11,-11,-2,-5,2,5,21]
    flv_born (:,11) = [11,-11,-2,-3,2,3,21]
    flv_born (:,12) = [11,-11,-2,-3,4,1,21]
    flv_born (:,13) = [11,-11,-2,-1,2,1,21]
    flv_born (:,14) = [11,-11,-5,-5,5,5,21]
    flv_born (:,15) = [11,-11,-5,-3,3,5,21]
    flv_born (:,16) = [11,-11,-5,-1,1,5,21]
    flv_born (:,17) = [11,-11,-5,5,21,21,21]
    flv_born (:,18) = [11,-11,-3,-3,3,3,21]
    flv_born (:,19) = [11,-11,-3,-1,1,3,21]
    flv_born (:,20) = [11,-11,-3,3,21,21,21]
    flv_born (:,21) = [11,-11,-1,-1,1,1,21]
    flv_born (:,22) = [11,-11,-1,1,21,21,21]

    flv_real (:,1)  = [11,-11,-4,-4,-4,4,4,4]
    flv_real (:,2)  = [11,-11,-4,-4,-2,2,4,4]
    flv_real (:,3)  = [11,-11,-4,-4,4,4,21,21]
    flv_real (:,4)  = [11,-11,-4,-4,-5,4,4,5]
    flv_real (:,5)  = [11,-11,-4,-4,-3,4,4,3]
    flv_real (:,6)  = [11,-11,-4,-4,-1,2,4,3]
    flv_real (:,7)  = [11,-11,-4,-4,-1,4,4,1]
    flv_real (:,8)  = [11,-11,-4,-2,-2,2,2,4]
    flv_real (:,9)  = [11,-11,-4,-2,2,4,21,21]
    flv_real (:,10) = [11,-11,-4,-2,-5,2,4,5]
    flv_real (:,11) = [11,-11,-4,-2,-3,2,4,3]
    flv_real (:,12) = [11,-11,-4,-2,-3,4,4,1]
    flv_real (:,13) = [11,-11,-4,-2,-1,2,2,3]
    flv_real (:,14) = [11,-11,-4,-2,-1,2,4,1]
    flv_real (:,15) = [11,-11,-4,4,21,21,21,21]
    flv_real (:,16) = [11,-11,-4,-5,4,5,21,21]
    flv_real (:,17) = [11,-11,-4,-5,-5,4,5,5]
    flv_real (:,18) = [11,-11,-4,-5,-3,4,3,5]
    flv_real (:,19) = [11,-11,-4,-5,-1,2,3,5]
    flv_real (:,20) = [11,-11,-4,-5,-1,4,1,5]
    flv_real (:,21) = [11,-11,-4,-3,4,3,21,21]
    flv_real (:,22) = [11,-11,-4,-3,-3,4,3,3]
    flv_real (:,23) = [11,-11,-4,-3,-1,2,3,3]
    flv_real (:,24) = [11,-11,-4,-3,-1,4,1,3]
    flv_real (:,25) = [11,-11,-4,-1,2,3,21,21]
    flv_real (:,26) = [11,-11,-4,-1,4,1,21,21]
    flv_real (:,27) = [11,-11,-4,-1,-1,2,1,3]
    flv_real (:,28) = [11,-11,-4,-1,-1,4,1,1]
    flv_real (:,29) = [11,-11,-2,-2,-2,2,2,2]
    flv_real (:,30) = [11,-11,-2,-2,2,2,21,21]
    flv_real (:,31) = [11,-11,-2,-2,-5,2,2,5]
    flv_real (:,32) = [11,-11,-2,-2,-3,2,2,3]
    flv_real (:,33) = [11,-11,-2,-2,-3,2,4,1]
    flv_real (:,34) = [11,-11,-2,-2,-1,2,2,1]
    flv_real (:,35) = [11,-11,-2,2,21,21,21,21]
    flv_real (:,36) = [11,-11,-2,-5,2,5,21,21]
    flv_real (:,37) = [11,-11,-2,-5,-5,2,5,5]
    flv_real (:,38) = [11,-11,-2,-5,-3,2,3,5]
    flv_real (:,39) = [11,-11,-2,-5,-3,4,1,5]
    flv_real (:,40) = [11,-11,-2,-5,-1,2,1,5]
    flv_real (:,41) = [11,-11,-2,-3,2,3,21,21]
    flv_real (:,42) = [11,-11,-2,-3,4,1,21,21]
    flv_real (:,43) = [11,-11,-2,-3,-3,2,3,3]
    flv_real (:,44) = [11,-11,-2,-3,-3,4,1,3]
    flv_real (:,45) = [11,-11,-2,-3,-1,2,1,3]
    flv_real (:,46) = [11,-11,-2,-3,-1,4,1,1]
    flv_real (:,47) = [11,-11,-2,-1,2,1,21,21]
    flv_real (:,48) = [11,-11,-2,-1,-1,2,1,1]
    flv_real (:,49) = [11,-11,-5,-5,-5,5,5,5]
    flv_real (:,50) = [11,-11,-5,-5,-3,3,5,5]
    flv_real (:,51) = [11,-11,-5,-5,-1,1,5,5]
    flv_real (:,52) = [11,-11,-5,-5,5,5,21,21]
    flv_real (:,53) = [11,-11,-5,-3,-3,3,3,5]
    flv_real (:,54) = [11,-11,-5,-3,-1,1,3,5]
    flv_real (:,55) = [11,-11,-5,-3,3,5,21,21]
    flv_real (:,56) = [11,-11,-5,-1,-1,1,1,5]
    flv_real (:,57) = [11,-11,-5,-1,1,5,21,21]
    flv_real (:,58) = [11,-11,-5,5,21,21,21,21]
    flv_real (:,59) = [11,-11,-3,-3,-3,3,3,3]
    flv_real (:,60) = [11,-11,-3,-3,-1,1,3,3]
    flv_real (:,61) = [11,-11,-3,-3,3,3,21,21]
    flv_real (:,62) = [11,-11,-3,-1,-1,1,1,3]
    flv_real (:,63) = [11,-11,-3,-1,1,3,21,21]
    flv_real (:,64) = [11,-11,-3,3,21,21,21,21]
    flv_real (:,65) = [11,-11,-1,-1,-1,1,1,1]
    flv_real (:,66) = [11,-11,-1,-1,1,1,21,21]
    flv_real (:,67) = [11,-11,-1,1,21,21,21,21]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 2, 3)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()

  end subroutine fks_regions_5

  subroutine fks_regions_6 (u)
    integer, intent(in) :: u
    integer :: n_flv_born, n_flv_real
    integer :: n_legs_born, n_legs_real
    integer :: n_in
    integer, dimension(:,:), allocatable :: flv_born, flv_real
    type(region_data_t) :: reg_data
    integer :: i, j
    integer, dimension(10) :: flavors
    write (u, "(A)") "* Test output: fks_regions_6"
    write (u, "(A)") "* Create table of singular regions for Drell Yan"
    write (u, "(A)")

    n_flv_born = 10; n_flv_real = 30
    n_legs_born = 4; n_legs_real = 5
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))
    flavors = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
    do i = 1, n_flv_born
       flv_born (3:4, i) = [11, -11]
    end do
    do j = 1, n_flv_born
       flv_born (1, j) = flavors (j)
       flv_born (2, j) = -flavors (j)
    end do

    do i = 1, n_flv_real
       flv_real (3:4, i) = [11, -11]
    end do
    i = 1
    do j = 1, n_flv_real
       if (mod (j, 3) == 1) then
          flv_real (1, j) = flavors (i)
          flv_real (2, j) = -flavors (i)
          flv_real (5, j) = 21
       else if (mod (j, 3) == 2) then
          flv_real (1, j) = flavors (i)
          flv_real (2, j) = 21
          flv_real (5, j) = flavors (i)
       else
          flv_real (1, j) = 21
          flv_real (2, j) = -flavors (i)
          flv_real (5, j) = -flavors (i)
          i = i + 1
       end if
    end do

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 2, 0)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    call write_separator (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()

    write (u, "(A)") "* Create table of singular regions for hadronic top decay"
    write (u, "(A)")
    n_flv_born = 1; n_flv_real = 1
    n_legs_born = 4; n_legs_real = 5
    n_in = 1
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))

    flv_born (:, 1) = [6, -5, 2, -1]
    flv_real (:, 1) = [6, -5, 2, -1, 21]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 0, 2)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)

    call write_separator (u)

    deallocate (flv_born, flv_real)
    call reg_data%final ()

    write (u, "(A)") "* Create table of singular regions for dijet s sbar -> jet jet"
    write (u, "(A)") "* With jet = u:d:gl"
    write (u, "(A)")

    n_flv_born = 3; n_flv_real = 3
    n_legs_born = 4; n_legs_real = 5
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))
    do i = 1, n_flv_born
       flv_born (1:2, i) = [3, -3]
    end do
    flv_born (3, :) = [1, 2, 21]
    flv_born (4, :) = [-1, -2, 21]

    do i = 1, n_flv_real
       flv_real (1:2, i) = [3, -3]
    end do
    flv_real (3, :) = [1, 2, 21]
    flv_real (4, :) = [-1, -2, 21]
    flv_real (5, :) = [21, 21, 21]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 0, 2)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)
    call reg_data%final ()

  end subroutine fks_regions_6

  subroutine fks_regions_7 (u)
    integer, intent(in) :: u
    integer :: n_flv_born, n_flv_real
    integer :: n_legs_born, n_legs_real
    integer :: n_in
    integer, dimension(:,:), allocatable :: flv_born, flv_real
    type(region_data_t) :: reg_data
    write (u, "(A)") "* Test output: fks_regions_7"
    write (u, "(A)") "* Create table of singular regions for ee -> qq"
    write (u, "(A)")

    n_flv_born = 1; n_flv_real = 1
    n_legs_born = 4; n_legs_real = 5
    n_in = 2

    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))
    flv_born (:, 1) = [11, -11, 2, -2]
    flv_real (:, 1) = [11, -11, 2, -2, 21]
    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("QCD"), 2, 0)
    call reg_data%write_latex (u)
    call reg_data%final ()

  end subroutine fks_regions_7

  subroutine fks_regions_8 (u)
    integer, intent(in) :: u
    integer :: n_flv_born, n_flv_real
    integer :: n_legs_born, n_legs_real
    integer :: n_in
    integer, dimension(:,:), allocatable :: flv_born, flv_real
    type(region_data_t) :: reg_data
    integer :: i, j
    integer, dimension(10) :: flavors
    write (u, "(A)") "* Test output: fks_regions_8"
    write (u, "(A)") "* Create table of singular regions for ee -> ee"
    write (u, "(A)")

    n_flv_born = 1; n_flv_real = 3
    n_legs_born = 4; n_legs_real = 5
    n_in = 2
    allocate (flv_born (n_legs_born, n_flv_born))
    allocate (flv_real (n_legs_real, n_flv_real))
    flv_born (:, 1) = [11, -11, -11, 11]

    flv_real (:, 1) = [11, -11, -11, 11, 22]
    flv_real (:, 2) = [11, 22, -11, 11, 11]
    flv_real (:, 3) = [22, -11, 11, -11, -11]

    call setup_region_data_for_test (n_in, flv_born, flv_real, reg_data, var_str ("EW"), 2, 0)
    call reg_data%check_consistency (.false., u)
    call reg_data%write (u)
    call reg_data%final ()

  end subroutine fks_regions_8


end module fks_regions_uti
