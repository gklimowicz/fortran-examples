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

submodule (dglap_remnant) dglap_remnant_s

  use numeric_utils
  use diagnostics
  use constants
  use physics_defs
  use pdg_arrays

  implicit none

contains

  module subroutine dglap_remnant_init &
       (dglap, settings, reg_data, isr_kinematics)
    class(dglap_remnant_t), intent(inout) :: dglap
    type(nlo_settings_t), intent(in), target :: settings
    type(region_data_t), intent(in), target :: reg_data
    integer :: n_flv_born
    type(isr_kinematics_t), intent(in), target :: isr_kinematics
    dglap%reg_data => reg_data
    n_flv_born = reg_data%get_n_flv_born ()
    allocate (dglap%sf_factors (reg_data%n_regions, 0:reg_data%n_in))
    dglap%sf_factors = zero
    dglap%settings => settings
    allocate (dglap%sqme_born(n_flv_born))
    dglap%sqme_born = zero
    allocate (dglap%sqme_color_c_extra (reg_data%n_legs_born, &
         reg_data%n_legs_born, reg_data%n_flv_born))
    dglap%sqme_color_c_extra = zero
    dglap%isr_kinematics => isr_kinematics
  end subroutine dglap_remnant_init

  module subroutine dglap_remnant_set_parameters (dglap, CA, CF, TR)
    class(dglap_remnant_t), intent(inout) :: dglap
    real(default), intent(in) :: CA, CF, TR
    dglap%CA = CA
    dglap%CF = CF
    dglap%TR = TR
  end subroutine dglap_remnant_set_parameters

  module subroutine dglap_remnant_evaluate &
       (dglap, alpha_coupling, separate_uborns, sqme_dglap)
    class(dglap_remnant_t), intent(inout) :: dglap
    real(default), dimension(2), intent(in) :: alpha_coupling
    logical, intent(in) :: separate_uborns
    real(default), intent(inout), dimension(:) :: sqme_dglap
    integer :: alr, emitter, i_corr
    real(default) :: sqme_alr
    logical, dimension(:,:,:), allocatable :: evaluated
    real(default) :: sb, fac_scale2
    sb = dglap%isr_kinematics%sqrts_born**2
    fac_scale2 = dglap%isr_kinematics%fac_scale**2
    allocate (evaluated(dglap%reg_data%get_n_flv_born (), &
         dglap%reg_data%get_n_flv_real (), dglap%reg_data%n_in))
    evaluated = .false.
    do alr = 1, dglap%reg_data%n_regions
       i_corr = 0
       if (dglap%reg_data%regions(alr)%nlo_correction_type == "QCD") then
          i_corr = 1
       else if (dglap%reg_data%regions(alr)%nlo_correction_type == "EW") then
          i_corr = 2
       end if
       if (allocated (dglap%settings%selected_alr)) then
          if (.not. any (dglap%settings%selected_alr == alr)) cycle
       end if
       sqme_alr = zero
       emitter = dglap%reg_data%regions(alr)%emitter
       if (emitter > dglap%reg_data%n_in .or. i_corr == 0) cycle
       associate (i_flv_born => dglap%reg_data%regions(alr)%uborn_index, &
               i_flv_real => dglap%reg_data%regions(alr)%real_index)
          if (emitter == 0) then
             do emitter = 1, 2
                if (evaluated(i_flv_born, i_flv_real, emitter)) cycle
                call evaluate_alr (alr, emitter, i_flv_born, &
                     i_flv_real, sqme_alr, evaluated)
             end do
          else if (emitter > 0) then
             if (evaluated(i_flv_born, i_flv_real, emitter)) cycle
             call evaluate_alr (alr, emitter, i_flv_born, &
                  i_flv_real, sqme_alr, evaluated)
          end if
          if (separate_uborns) then
             sqme_dglap(i_flv_born) = sqme_dglap(i_flv_born) &
                  + alpha_coupling (i_corr)/ twopi * sqme_alr
          else
             sqme_dglap(1) = sqme_dglap(1) &
                  + alpha_coupling (i_corr) / twopi * sqme_alr
          end if
       end associate
    end do

  contains
    function p_hat_gtogg (z)
      real(default) :: p_hat_gtogg
      real(default), intent(in) :: z
      real(default) :: onemz
      onemz = one - z

      p_hat_gtogg = two * dglap%CA * (z + onemz**2 / z + z * onemz**2)
    end function p_hat_gtogg

    function p_hat_gtoqq (z)
      real(default) :: p_hat_gtoqq
      real(default), intent(in) :: z
      real(default) :: onemz
      onemz = one - z

      p_hat_gtoqq = dglap%CF * onemz / z * (one + onemz**2)
    end function p_hat_gtoqq

    function p_hat_qtogq (z)
      real(default) :: p_hat_qtogq
      real(default), intent(in) :: z
      real(default) :: onemz
      onemz = one - z

      p_hat_qtogq = dglap%TR * (onemz - two * z * onemz**2)
    end function p_hat_qtogq

    function p_hat_qtoqg (z)
      real(default) :: p_hat_qtoqg
      real(default), intent(in) :: z
      p_hat_qtoqg = dglap%CF * (one + z**2)
    end function p_hat_qtoqg

    function p_derived_gtogg (z)
      real(default) :: p_derived_gtogg
      real(default), intent(in) :: z
      p_derived_gtogg = zero
    end function p_derived_gtogg

    function p_derived_gtoqq (z)
      real(default) :: p_derived_gtoqq
      real(default), intent(in) :: z
      p_derived_gtoqq = -dglap%CF * z
    end function p_derived_gtoqq

    function p_derived_qtogq (z)
      real(default) :: p_derived_qtogq
      real(default), intent(in) :: z
      real(default) :: onemz
      onemz = one - z

      p_derived_qtogq = -two * dglap%TR * z * onemz
    end function p_derived_qtogq

    function p_derived_qtoqg (z)
      real(default) :: p_derived_qtoqg
      real(default), intent(in) :: z
      real(default) :: onemz
      onemz = one - z

      p_derived_qtoqg = -dglap%CF * onemz
    end function p_derived_qtoqg

  subroutine evaluate_alr (alr, emitter, i_flv_born, i_flv_real, sqme_alr, evaluated)
    integer, intent(in) :: alr, emitter, i_flv_born, i_flv_real
    real(default), intent(inout) :: sqme_alr
    logical, intent(inout), dimension(:,:,:) :: evaluated
    real(default) :: z, jac
    real(default) :: factor, factor_soft, plus_dist_remnant
    real(default) :: xb, onemz
    real(default) :: sqme_scaled, sqme_born_dglap
    real(default) :: charge_rad2, charge_em2
    integer :: flv_em, flv_rad, N_col, i
    N_col = 1
    sqme_born_dglap = zero
    associate (template => dglap%settings%fks_template)
       z = dglap%isr_kinematics%z(emitter)
       flv_rad = dglap%reg_data%regions(alr)%flst_real%flst(dglap%reg_data%n_legs_real)
       flv_em = dglap%reg_data%regions(alr)%flst_real%flst(emitter)
       charge_rad2 = dglap%reg_data%regions(alr)%flst_real%charge(dglap%reg_data%n_legs_real)**2
       charge_em2 = dglap%reg_data%regions(alr)%flst_real%charge(emitter)**2
       if (dglap%reg_data%regions(alr)%nlo_correction_type == "QCD") then
          call dglap%set_parameters (CA = CA, CF = CF, TR = TR)
       else if (dglap%reg_data%regions(alr)%nlo_correction_type == "EW") then
          if (is_quark(flv_rad)) N_col = NC
          call dglap%set_parameters (CA = zero, CF = charge_em2, TR = N_col*charge_rad2)
       end if
       jac = dglap%isr_kinematics%jacobian(emitter)
       onemz = one - z
       factor = log (sb * template%delta_i / two / z / fac_scale2) / &
            onemz + two * log (onemz) / onemz
       factor_soft = log (sb * template%delta_i / two / fac_scale2) / &
            onemz + two * log (onemz) / onemz
       xb = dglap%isr_kinematics%x(emitter)
       plus_dist_remnant = log ((one - xb) / template%xi_cut) * log (sb * template%delta_i / &
            two / fac_scale2) + (log (one - xb)**2 - log (template%xi_cut)**2)
    end associate
    if (dglap%reg_data%nlo_correction_type == "EW" .and. &
       dglap%reg_data%regions(alr)%nlo_correction_type == "QCD" .and. &
       qcd_ew_interferences (dglap%reg_data%regions(alr)%flst_uborn%flst)) then
       do i = 1, size (dglap%reg_data%regions(alr)%flst_uborn%flst)
          if (is_quark (dglap%reg_data%regions(alr)%flst_uborn%flst (i))) then
             sqme_born_dglap = -dglap%sqme_color_c_extra (i, i, i_flv_born)/CF
             exit
          end if
       end do
    else
       sqme_born_dglap = dglap%sqme_born(i_flv_born)
    end if
    sqme_scaled = sqme_born_dglap * dglap%sf_factors(alr, emitter)
    if (is_massless_vector (flv_em) .and. is_massless_vector (flv_rad)) then
       sqme_alr = sqme_alr + p_hat_gtogg(z) * factor / z * sqme_scaled * jac &
            - p_hat_gtogg(one) * factor_soft * sqme_born_dglap * jac &
            + p_hat_gtogg(one) * plus_dist_remnant * sqme_born_dglap
    else if (is_fermion (flv_em) .and. is_massless_vector (flv_rad)) then
       sqme_alr = sqme_alr + p_hat_qtoqg(z) * factor / z * sqme_scaled * jac &
            - p_derived_qtoqg(z) / z * sqme_scaled * jac &
            - p_hat_qtoqg(one) * factor_soft * sqme_born_dglap * jac &
            + p_hat_qtoqg(one) * plus_dist_remnant * sqme_born_dglap
    else if (is_fermion (flv_em) .and. is_fermion (flv_rad)) then
       sqme_alr = sqme_alr + (p_hat_gtoqq(z) * factor - p_derived_gtoqq(z)) / z * jac * &
            sqme_scaled
    else if (is_massless_vector (flv_em) .and. is_fermion (flv_rad)) then
       sqme_alr = sqme_alr + (p_hat_qtogq(z) * factor - p_derived_qtogq(z)) / z * sqme_scaled * jac
    else
       sqme_alr = sqme_alr + zero
    end if
    evaluated(i_flv_born, i_flv_real, emitter) = .true.
  end subroutine evaluate_alr
  end subroutine dglap_remnant_evaluate

  module subroutine dglap_remnant_final (dglap)
    class(dglap_remnant_t), intent(inout) :: dglap
    if (associated (dglap%isr_kinematics)) nullify (dglap%isr_kinematics)
    if (associated (dglap%reg_data)) nullify (dglap%reg_data)
    if (associated (dglap%settings)) nullify (dglap%settings)
    if (allocated (dglap%sqme_born)) deallocate (dglap%sqme_born)
    if (allocated (dglap%sf_factors)) deallocate (dglap%sf_factors)
  end subroutine dglap_remnant_final


end submodule dglap_remnant_s

