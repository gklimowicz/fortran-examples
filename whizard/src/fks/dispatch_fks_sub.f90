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

submodule (dispatch_fks) dispatch_fks_s

  use string_utils, only: split_string

  implicit none

contains

  module subroutine dispatch_fks_setup (fks_template, var_list)
    type(fks_template_t), intent(inout) :: fks_template
    type(var_list_t), intent(in) :: var_list
    real(default) :: fks_dij_exp1, fks_dij_exp2
    type(string_t) :: fks_mapping_type
    logical :: subtraction_disabled
    type(string_t) :: exclude_from_resonance
    fks_dij_exp1 = &
         var_list%get_rval (var_str ("fks_dij_exp1"))
    fks_dij_exp2 = &
         var_list%get_rval (var_str ("fks_dij_exp2"))
    fks_mapping_type = &
         var_list%get_sval (var_str ("$fks_mapping_type"))
    subtraction_disabled = &
         var_list%get_lval (var_str ("?disable_subtraction"))
    exclude_from_resonance = &
         var_list%get_sval (var_str ("$resonances_exclude_particles"))
    if (exclude_from_resonance /= var_str ("default")) &
       call split_string (exclude_from_resonance, var_str (":"), &
       fks_template%excluded_resonances)
    call fks_template%set_parameters ( &
         exp1 = fks_dij_exp1, exp2 = fks_dij_exp2, &
         xi_min = var_list%get_rval (var_str ("fks_xi_min")), &
         y_max = var_list%get_rval (var_str ("fks_y_max")), &
         xi_cut = var_list%get_rval (var_str ("fks_xi_cut")), &
         delta_o = var_list%get_rval (var_str ("fks_delta_o")), &
         delta_i = var_list%get_rval (var_str ("fks_delta_i")))
    select case (char (fks_mapping_type))
    case ("default")
       call fks_template%set_mapping_type (FKS_DEFAULT)
    case ("resonances")
       call fks_template%set_mapping_type (FKS_RESONANCES)
    end select
    fks_template%subtraction_disabled = subtraction_disabled
    fks_template%n_f = var_list%get_ival (var_str ("alphas_nf"))
  end subroutine dispatch_fks_setup


end submodule dispatch_fks_s

