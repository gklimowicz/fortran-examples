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

submodule (virtual) virtual_s

  use debug_master, only: debug_on
  use constants
  use numeric_utils
  use diagnostics
  use physics_defs
  use sm_physics
  use flavors
  use nlo_data, only: ASSOCIATED_LEG_PAIR, get_threshold_momenta

  implicit none

contains

  module subroutine virtual_init &
       (virt, flv_born, n_in, settings, model, has_pdfs)
    class(virtual_t), intent(inout) :: virt
    integer, intent(in), dimension(:,:) :: flv_born
    integer, intent(in) :: n_in
    type(nlo_settings_t), intent(in) :: settings
    class(model_data_t), intent(in) :: model
    logical, intent(in) :: has_pdfs
    integer :: i_flv, n_corr
    n_corr = 2
    virt%n_legs = size (flv_born, 1); virt%n_flv = size (flv_born, 2)
    virt%n_in = n_in
    allocate (virt%sqme_born (virt%n_flv))
    allocate (virt%sqme_virt_fin (virt%n_flv))
    allocate (virt%sqme_color_c (virt%n_legs, virt%n_legs, virt%n_flv))
    allocate (virt%sqme_charge_c (virt%n_legs, virt%n_legs, virt%n_flv))
    allocate (virt%gamma_0 (virt%n_legs, virt%n_flv, n_corr), &
       virt%gamma_p (virt%n_legs, virt%n_flv, n_corr), &
       virt%c_flv (virt%n_legs, virt%n_flv, n_corr))
    call virt%init_constants (flv_born, settings%fks_template%n_f, model)
    allocate (virt%n_is_neutrinos (virt%n_flv))
    virt%n_is_neutrinos = 0
    do i_flv = 1, virt%n_flv
       if (is_neutrino (flv_born(1, i_flv))) &
          virt%n_is_neutrinos(i_flv) = virt%n_is_neutrinos(i_flv) + 1
       if (is_neutrino (flv_born(2, i_flv))) &
          virt%n_is_neutrinos(i_flv) = virt%n_is_neutrinos(i_flv) + 1
    end do
    select case (char (settings%virtual_selection))
    case ("Full", "OLP", "Subtraction")
       virt%selection = settings%virtual_selection
    case default
       call msg_fatal ('Virtual selection: Possible values are "Full", "OLP" or "Subtraction')
    end select
    virt%settings = settings
    virt%has_pdfs = has_pdfs
  contains

    function is_neutrino (flv) result (neutrino)
      integer, intent(in) :: flv
      logical :: neutrino
      neutrino = (abs(flv) == 12 .or. abs(flv) == 14 .or. abs(flv) == 16)
    end function is_neutrino

  end subroutine virtual_init

  module subroutine virtual_init_constants (virt, flv_born, nf_input, model)
    class(virtual_t), intent(inout) :: virt
    integer, intent(in), dimension(:,:) :: flv_born
    integer, intent(in) :: nf_input
    type(string_t), dimension(2) :: corr_type
    class(model_data_t), intent(in) :: model
    type(field_data_t), pointer :: field
    integer :: i_part, i_flv, pdg, i_corr
    real(default) :: nf, CA_factor, TR_sum
    real(default), dimension(:,:), allocatable :: CF_factor, TR_factor
    type(flavor_t) :: flv
    allocate (CF_factor (size (flv_born, 1), size (flv_born, 2)), &
         TR_factor (size (flv_born, 1), size (flv_born, 2)))
    corr_type (1) = "QCD"; corr_type (2) = "EW"
    do i_corr = 1, 2
       TR_sum = 0
       if (i_corr == 1) then
          CA_factor = CA; CF_factor = CF; TR_factor = TR
          nf = real(nf_input, default)
          TR_sum = nf * TR
       else
          CA_factor = zero
          do i_flv = 1, size (flv_born, 2)
             do i_part = 1, size (flv_born, 1)
                call flv%init (flv_born(i_part, i_flv), model)
                CF_factor(i_part, i_flv) = (flv%get_charge ())**2
                TR_factor(i_part, i_flv) = (flv%get_charge ())**2
                if (is_quark (flv_born (i_part, i_flv))) &
                     TR_factor(i_part, i_flv) = NC * TR_factor(i_part, i_flv)
             end do
          end do
          do pdg = 1, nf_input
             field => model%get_field_ptr (pdg)
             TR_sum = TR_sum + NC*field%get_charge()**2
          end do
          do pdg = 11, 15, 2
             field => model%get_field_ptr (pdg)
             if (field%get_mass() > 0) exit
             TR_sum = TR_sum + field%get_charge()**2
          end do
       end if
       do i_flv = 1, size (flv_born, 2)
          do i_part = 1, size (flv_born, 1)
             if (is_massless_vectorboson (flv_born(i_part, i_flv), corr_type (i_corr))) then
                virt%gamma_0(i_part, i_flv, i_corr) = 11._default / 6._default * CA_factor &
                     - two / three * TR_sum
                virt%gamma_p(i_part, i_flv, i_corr) = (67._default / 9._default &
                     - two * pi**2 / three) * CA_factor &
                     - 23._default / 9._default * TR_sum
                virt%c_flv(i_part, i_flv, i_corr) = CA_factor
             else if (is_corresponding_fermion (flv_born(i_part, i_flv), corr_type (i_corr))) then
                virt%gamma_0(i_part, i_flv, i_corr) = 1.5_default * CF_factor(i_part, i_flv)
                virt%gamma_p(i_part, i_flv, i_corr) = (6.5_default - two * pi**2 / three) * CF_factor(i_part, i_flv)
                virt%c_flv(i_part, i_flv, i_corr) = CF_factor(i_part, i_flv)
             else if (is_massive_vectorboson (flv_born(i_part, i_flv), corr_type (i_corr))) then
                virt%gamma_0(i_part, i_flv, i_corr) = zero
                virt%gamma_p(i_part, i_flv, i_corr) = zero
                virt%c_flv(i_part, i_flv, i_corr) = CF_factor(i_part, i_flv)
             else
                virt%gamma_0(i_part, i_flv, i_corr) = zero
                virt%gamma_p(i_part, i_flv, i_corr) = zero
                virt%c_flv(i_part, i_flv, i_corr) = zero
             end if
          end do
       end do
    end do
  contains
    function is_massless_vectorboson (pdg_nr, nlo_corr_type)
      logical :: is_massless_vectorboson
      integer, intent(in) :: pdg_nr
      type(string_t), intent(in) :: nlo_corr_type
      is_massless_vectorboson = .false.
      if (nlo_corr_type == "QCD") then
         is_massless_vectorboson = is_gluon (pdg_nr)
      else if (nlo_corr_type == "EW") then
         is_massless_vectorboson = is_photon (pdg_nr)
      end if
    end function is_massless_vectorboson
    function is_corresponding_fermion (pdg_nr, nlo_corr_type)
      logical :: is_corresponding_fermion
      integer, intent(in) :: pdg_nr
      type(string_t), intent(in) :: nlo_corr_type
      is_corresponding_fermion = .false.
      if (nlo_corr_type == "QCD") then
         is_corresponding_fermion = is_quark (pdg_nr)
      else if (nlo_corr_type == "EW") then
         is_corresponding_fermion = is_fermion (pdg_nr)
      end if
    end function is_corresponding_fermion
    function is_massive_vectorboson (pdg_nr, nlo_corr_type)
      logical :: is_massive_vectorboson
      integer, intent(in) :: pdg_nr
      type(string_t), intent(in) :: nlo_corr_type
      is_massive_vectorboson = .false.
      if (nlo_corr_type == "EW") then
         is_massive_vectorboson = is_massive_vector (pdg_nr)
      end if
    end function is_massive_vectorboson
  end subroutine virtual_init_constants

  module subroutine virtual_set_ren_scale (virt, ren_scale)
    class(virtual_t), intent(inout) :: virt
    real(default), allocatable, intent(in) :: ren_scale
    if (allocated (ren_scale)) then
       if (allocated (virt%ren_scale)) then
          virt%ren_scale = ren_scale
       else
          allocate (virt%ren_scale, source=ren_scale)
       end if
    end if
  end subroutine virtual_set_ren_scale

  module subroutine virtual_set_fac_scale (virt, p, fac_scale)
    class(virtual_t), intent(inout) :: virt
    type(vector4_t), dimension(:), intent(in) :: p
    real(default), intent(in), optional :: fac_scale
    if (present (fac_scale)) then
       virt%fac_scale = fac_scale
    else
       virt%fac_scale = (p(1) + p(2))**1
    end if
  end subroutine virtual_set_fac_scale

  module subroutine virtual_set_ellis_sexton_scale (virt, Q)
    class(virtual_t), intent(inout) :: virt
    real(default), allocatable, intent(in) :: Q
    if (allocated (Q)) then
       if (allocated (virt%es_scale2)) then
          virt%es_scale2 = Q * Q
       else
          allocate (virt%es_scale2, source=Q*Q)
       end if
    end if
  end subroutine virtual_set_ellis_sexton_scale

  module subroutine virtual_evaluate (virt, reg_data, alpha_coupling, &
       p_born, separate_uborns, sqme_virt)
    class(virtual_t), intent(inout) :: virt
    type(region_data_t), intent(in) :: reg_data
    real(default), dimension(2), intent(in) :: alpha_coupling
    type(vector4_t), intent(in), dimension(:)  :: p_born
    logical, intent(in) :: separate_uborns
    real(default), dimension(:), intent(inout) :: sqme_virt
    integer, dimension(:), allocatable :: eqv_flv_index
    real(default), dimension(:), allocatable :: sqme_virt_arr
    real(default) :: s, s_o_Q2, es_scale2
    real(default), dimension(reg_data%n_flv_born) :: QB, BI
    integer :: i_flv, ii_flv, alr, i_corr
    logical, dimension(:), allocatable :: flv_evaluated
    integer, dimension(:), allocatable :: corr_index
    logical :: alr_qcd, alr_ew
    allocate (flv_evaluated(reg_data%n_flv_born))
    allocate (sqme_virt_arr(reg_data%n_flv_born))
    sqme_virt_arr = zero
    flv_evaluated = .false.
    if (virt%bad_point) return
    if (allocated (virt%es_scale2)) then
       es_scale2 = virt%es_scale2
    else
       if (allocated (virt%ren_scale)) then
          es_scale2 = virt%ren_scale**2
       else
          es_scale2 = virt%fac_scale**2
       end if
    end if
    if (debug2_active (D_VIRTUAL)) then
       print *, 'Compute virtual component using alpha = ', alpha_coupling
       print *, 'Virtual selection: ', char (virt%selection)
       print *, 'virt%es_scale2 =    ', es_scale2 !!! Debugging
    end if
    s = sum (p_born(1 : virt%n_in))**2
    if (virt%settings%factorization_mode == FACTORIZATION_THRESHOLD) &
         call set_s_for_threshold ()
    s_o_Q2 = s / es_scale2 * virt%settings%fks_template%xi_cut**2
    eqv_flv_index = reg_data%eqv_flv_index_born
    do i_flv = 1, reg_data%n_flv_born
       alr_qcd = .false.; alr_ew = .false.
       do alr = 1, reg_data%n_regions
          if (i_flv == reg_data%regions(alr)%uborn_index) then
             if (reg_data%regions(alr)%nlo_correction_type == "QCD") then
                alr_qcd = .true.
             else if (reg_data%regions(alr)%nlo_correction_type == "EW") then
                alr_ew = .true.
             end if
          end if
       end do
       if (alr_qcd .and. alr_ew) then
          allocate (corr_index (2))
          corr_index (1) = 1; corr_index (2) = 2
       else
          allocate (corr_index (1))
          corr_index (1) = 0
          if (alr_qcd) then
             corr_index (1) = 1
          else if (alr_ew) then
             corr_index (1) = 2
          end if
       end if
       if (.not. flv_evaluated(eqv_flv_index(i_flv)) .and. corr_index(1) > 0) then
          if (virt%selection == var_str ("Full") .or. virt%selection == var_str ("OLP")) then
             !!! A factor of alpha_coupling/twopi is assumed to be included in vfin
             sqme_virt_arr(i_flv) = sqme_virt_arr(i_flv) + virt%sqme_virt_fin(i_flv)
          end if
          do i_corr = 1, size (corr_index)
             QB = zero; BI = zero
             if (virt%selection == var_str ("Full") .or. &
                  virt%selection == var_str ("Subtraction")) then
                call virt%evaluate_initial_state (i_flv, corr_index (i_corr), reg_data, QB)
                call virt%compute_collinear_contribution &
                     (i_flv, corr_index (i_corr), p_born, sqrt(s), reg_data, QB)
                select case (virt%settings%factorization_mode)
                case (FACTORIZATION_THRESHOLD)
                   call virt%compute_eikonals_threshold (i_flv, p_born, s_o_Q2, QB, BI)
                case default
                   call virt%compute_massive_self_eikonals &
                        (i_flv, corr_index (i_corr), p_born, s_o_Q2, reg_data, QB)
                   call virt%compute_eikonals &
                        (i_flv, corr_index (i_corr), p_born, s_o_Q2, reg_data, BI)
                end select
                if (debug2_active (D_VIRTUAL)) then
                   if (corr_index (i_corr) == 1) then
                      print *, 'Correction type: QCD'
                   else
                      print *, 'Correction type: EW'
                   end if
                   print *, 'Evaluate i_flv: ', i_flv
                   print *, 'sqme_born: ', virt%sqme_born (i_flv)
                   print *, 'Q * sqme_born: ', alpha_coupling / twopi * QB(i_flv)
                   print *, 'BI: ', alpha_coupling / twopi * BI(i_flv)
                   print *, 'vfin: ', virt%sqme_virt_fin (i_flv)
                end if
                sqme_virt_arr(i_flv) = sqme_virt_arr(i_flv) &
                      + alpha_coupling (corr_index (i_corr))/ twopi * (QB(i_flv) + BI(i_flv))
             end if
          end do
          if (.not. (debug_active (D_VIRTUAL) .or. &
               debug2_active (D_VIRTUAL))) flv_evaluated(eqv_flv_index(i_flv)) = .true.
       else
          sqme_virt_arr(i_flv) = sqme_virt_arr(eqv_flv_index(i_flv))
       end if
       if (separate_uborns) then
          sqme_virt(i_flv) = sqme_virt(i_flv) + sqme_virt_arr(i_flv)
       else
          sqme_virt(1) = sqme_virt(1) + sqme_virt_arr(i_flv)
       end if
       deallocate (corr_index)
    end do
    if (debug2_active (D_VIRTUAL)) then
       call msg_debug2 (D_VIRTUAL, "virtual-subtracted matrix element(s): ")
       print *, sqme_virt
    end if
    do i_flv = 1, reg_data%n_flv_born
       if (virt%n_is_neutrinos(i_flv) > 0) &
            sqme_virt = sqme_virt * virt%n_is_neutrinos(i_flv) * two
    end do
  contains
    subroutine set_s_for_threshold ()
      use ttv_formfactors, only: m1s_to_mpole
      real(default) :: mtop2
      mtop2 = m1s_to_mpole (sqrt(s))**2
      if (s < four * mtop2) s = four * mtop2
    end subroutine set_s_for_threshold

  end subroutine virtual_evaluate

  module subroutine virtual_compute_eikonals (virtual, i_flv, i_corr, &
       p_born, s_o_Q2, reg_data, BI)
    class(virtual_t), intent(inout) :: virtual
    integer, intent(in) :: i_flv, i_corr
    type(vector4_t), intent(in), dimension(:)  :: p_born
    real(default), intent(in) :: s_o_Q2
    type(region_data_t), intent(in) :: reg_data
    real(default), intent(inout), dimension(:) :: BI
    integer :: i, j
    real(default) :: I_ij, BI_tmp
    BI_tmp = zero
    ! TODO vincent_r: Split the procedure into one computing QCD eikonals
    ! and one computing QED eikonals.
    ! TODO vincent_r: In the best case, remove the dependency on
    ! reg_data completely.
    associate (flst_born => reg_data%flv_born(i_flv))
       do i = 1, virtual%n_legs
          do j = 1, virtual%n_legs
             if (i /= j) then
                if (i_corr == 1) then
                   if (flst_born%colored(i) .and. flst_born%colored(j)) then
                      I_ij = compute_eikonal_factor &
                           (p_born, flst_born%massive, i, j, s_o_Q2)
                      BI_tmp = BI_tmp + &
                           virtual%sqme_color_c (i, j, i_flv) * I_ij
                      if (debug2_active (D_VIRTUAL)) &
                           print *, 'b_ij: ', i, j, &
                           virtual%sqme_color_c (i, j, i_flv), 'I_ij: ', I_ij
                   end if
                else if (i_corr == 2) then
                   if (flst_born%charge (i) /= 0 .and. flst_born%charge(j) /= 0) then
                      I_ij = compute_eikonal_factor (p_born, flst_born%massive, &
                           i, j, s_o_Q2)
                      BI_tmp = BI_tmp + virtual%sqme_charge_c (i, j, i_flv) * I_ij
                      if (debug2_active (D_VIRTUAL)) &
                           print *, 'b_ij: ', &
                           virtual%sqme_charge_c (i, j, i_flv), 'I_ij: ', I_ij
                   end if
                end if
             else if (debug2_active (D_VIRTUAL)) then
                if (i_corr == 1) then
                   print *, 'b_ij: ', i, j, &
                        virtual%sqme_color_c (i, j, i_flv), 'I_ij: ', I_ij
                else if (i_corr == 2) then
                   print *, 'b_ij: ', i, j, &
                        virtual%sqme_charge_c (i, j, i_flv), 'I_ij: ', I_ij
                end if
             end if
          end do
       end do
       if (virtual%settings%use_internal_color_correlations .or. i_corr == 2) &
            BI_tmp = BI_tmp * virtual%sqme_born (i_flv)
    end associate
    BI(i_flv) = BI(i_flv) + BI_tmp
  end subroutine virtual_compute_eikonals

  module subroutine virtual_compute_eikonals_threshold (virtual, i_flv, &
         p_born, s_o_Q2, QB, BI)
    class(virtual_t), intent(in) :: virtual
    integer, intent(in) :: i_flv
    type(vector4_t), intent(in), dimension(:) :: p_born
    real(default), intent(in) :: s_o_Q2
    real(default), intent(inout), dimension(:) :: QB
    real(default), intent(inout), dimension(:) :: BI
    type(vector4_t), dimension(4) :: p_thr
    integer :: leg
    BI = zero; p_thr = get_threshold_momenta (p_born)
    call compute_massive_self_eikonals (virtual%sqme_born(i_flv), QB(i_flv))
    do leg = 1, 2
       BI(i_flv) = BI(i_flv) + evaluate_leg_pair (ASSOCIATED_LEG_PAIR(leg), i_flv)
    end do
  contains
    subroutine compute_massive_self_eikonals (sqme_born, QB)
      real(default), intent(in) :: sqme_born
      real(default), intent(inout) :: QB
      integer :: i
      if (debug_on) call msg_debug2 (D_VIRTUAL, "compute_massive_self_eikonals")
      if (debug_on) call msg_debug2 (D_VIRTUAL, "s_o_Q2", s_o_Q2)
      if (debug_on) call msg_debug2 (D_VIRTUAL, "log (s_o_Q2)", log (s_o_Q2))
      do i = 1, 4
         QB = QB - (cf * (log (s_o_Q2) - 0.5_default * I_m_eps (p_thr(i)))) &
              * sqme_born
      end do
    end subroutine compute_massive_self_eikonals

    function evaluate_leg_pair (i_start, i_flv) result (b_ij_times_I)
      real(default) :: b_ij_times_I
      integer, intent(in) :: i_start, i_flv
      real(default) :: I_ij
      integer :: i, j
      b_ij_times_I = zero
      do i = i_start, i_start + 1
         do j = i_start, i_start + 1
            if (i /= j) then
               I_ij = compute_eikonal_factor &
                    (p_thr, [.true., .true., .true., .true.], i, j, s_o_Q2)
               b_ij_times_I = b_ij_times_I + &
                    virtual%sqme_color_c (i, j, i_flv) * I_ij
               if (debug2_active (D_VIRTUAL)) &
                  print *, 'b_ij: ', virtual%sqme_color_c (i, j, i_flv), 'I_ij: ', I_ij
            end if
         end do
      end do
      if (virtual%settings%use_internal_color_correlations) &
           b_ij_times_I = b_ij_times_I * virtual%sqme_born (i_flv)
      if (debug2_active (D_VIRTUAL)) then
         print *, 'internal color: ', virtual%settings%use_internal_color_correlations
         print *, 'b_ij_times_I =    ', b_ij_times_I
         print *, 'QB           =    ', QB
      end if
    end function evaluate_leg_pair
  end subroutine virtual_compute_eikonals_threshold

  module subroutine virtual_set_bad_point (virt, value)
    class(virtual_t), intent(inout) :: virt
    logical, intent(in) :: value
    virt%bad_point = value
  end subroutine virtual_set_bad_point

  module subroutine virtual_evaluate_initial_state &
       (virt, i_flv, i_corr, reg_data, QB)
    class(virtual_t), intent(inout) :: virt
    type(region_data_t), intent(in) :: reg_data
    integer, intent(in) :: i_flv, i_corr
    real(default), intent(inout), dimension(:) :: QB
    real(default) :: sqme_born_virt, es_scale2
    integer :: i
    if (allocated (virt%es_scale2)) then
       es_scale2 = virt%es_scale2
    else
       if (allocated (virt%ren_scale)) then
          es_scale2 = virt%ren_scale**2
       else
          es_scale2 = virt%fac_scale**2
       end if
    end if
    sqme_born_virt = zero
    if (reg_data%nlo_correction_type == "EW" .and. i_corr == 1 &
       .and. qcd_ew_interferences (reg_data%flv_born(i_flv)%flst)) then
       do i = 1, size (reg_data%flv_born(i_flv)%flst)
          if (is_quark (reg_data%flv_born(i_flv)%flst (i))) then
             sqme_born_virt = -virt%sqme_color_c (i, i, i_flv)/CF
             exit
          end if
       end do
    else
       sqme_born_virt = virt%sqme_born (i_flv)
    end if
    if (virt%n_in == 2) then
       do i = 1, virt%n_in
          QB(i_flv) = QB(i_flv) - (virt%gamma_0 (i, i_flv, i_corr) &
               + two * virt%c_flv(i, i_flv, i_corr) &
               * log (virt%settings%fks_template%xi_cut)) &
               * log(virt%fac_scale**2 / es_scale2) * sqme_born_virt
       end do
    end if
  end subroutine virtual_evaluate_initial_state

  module subroutine virtual_compute_collinear_contribution &
       (virt, i_flv, i_corr, p_born, sqrts, reg_data, QB)
    class(virtual_t), intent(inout) :: virt
    integer, intent(in) :: i_flv, i_corr
    type(vector4_t), dimension(:), intent(in) :: p_born
    real(default), intent(in) :: sqrts
    type(region_data_t), intent(in) :: reg_data
    real(default), intent(inout), dimension(:) :: QB
    real(default) :: s1, s2, s3, s4, s5
    real(default) :: sqme_born_virt
    integer :: alr, em, i
    real(default) :: E_em, xi_max, log_xi_max, E_tot2, es_scale2
    logical, dimension(virt%n_flv, virt%n_legs) :: evaluated
    integer :: i_contr
    type(vector4_t) :: k_res
    type(lorentz_transformation_t) :: L_to_resonance
    evaluated = .false.
    if (allocated (virt%es_scale2)) then
       es_scale2 = virt%es_scale2
    else
       if (allocated (virt%ren_scale)) then
          es_scale2 = virt%ren_scale**2
       else
          es_scale2 = virt%fac_scale**2
       end if
    end if
    sqme_born_virt = zero
    if (reg_data%nlo_correction_type == "EW" .and. i_corr == 1 &
       .and. qcd_ew_interferences (reg_data%flv_born(i_flv)%flst)) then
       do i = 1, size (reg_data%flv_born(i_flv)%flst)
          if (is_quark (reg_data%flv_born(i_flv)%flst (i))) then
             sqme_born_virt = -virt%sqme_color_c (i, i, i_flv)/CF
             exit
          end if
       end do
    else
       sqme_born_virt = virt%sqme_born (i_flv)
    end if
    do alr = 1, reg_data%n_regions
       if (i_flv /= reg_data%regions(alr)%uborn_index) cycle
       em = reg_data%regions(alr)%emitter
       if (em <= virt%n_in) cycle
       if (evaluated(i_flv, em)) cycle
       !!! Collinear terms only for massless particles
       if (reg_data%regions(alr)%flst_uborn%massive(em)) cycle
       E_em = p_born(em)%p(0)
       if (allocated (reg_data%alr_contributors)) then
          i_contr = reg_data%alr_to_i_contributor (alr)
          k_res = get_resonance_momentum (p_born, reg_data%alr_contributors(i_contr)%c)
          E_tot2 = k_res%p(0)**2
          L_to_resonance = inverse (boost (k_res, k_res**1))
          xi_max = two * space_part_norm (L_to_resonance * p_born(em)) / k_res%p(0)
       else
          E_tot2 = sqrts**2
          xi_max = two * E_em / sqrts
       end if
       log_xi_max = log (xi_max)
       associate (xi_cut => virt%settings%fks_template%xi_cut, delta_o => virt%settings%fks_template%delta_o)
         if (virt%settings%virtual_resonance_aware_collinear) then
            if (debug_active (D_VIRTUAL)) &
                 call msg_debug (D_VIRTUAL, "Using resonance-aware collinear subtraction")
            s1 = virt%gamma_p(em, i_flv, i_corr)
            s2 = two * (log (sqrts / (two * E_em)) + log_xi_max) * &
                 (log (sqrts / (two * E_em)) + log_xi_max + log (es_scale2 / sqrts**2)) &
                 * virt%c_flv(em, i_flv, i_corr)
            s3 = two * log_xi_max * &
                 (log_xi_max - log (es_scale2 / E_tot2)) * virt%c_flv(em, i_flv, i_corr)
            s4 = (log (es_scale2 / E_tot2) - two * log_xi_max) &
                 * virt%gamma_0(em, i_flv, i_corr)
            QB(i_flv) = QB(i_flv) + (s1 + s2 + s3 + s4) * sqme_born_virt
         else
            if (debug_active (D_VIRTUAL)) &
                 call msg_debug (D_VIRTUAL, "Using old-fashioned collinear subtraction")
            s1 = virt%gamma_p(em, i_flv, i_corr)
            s2 = log (delta_o * sqrts**2 / (two * es_scale2)) &
                 * virt%gamma_0(em,i_flv, i_corr)
            s3 = log (delta_o * sqrts**2 / (two * es_scale2)) * two &
                 * virt%c_flv(em,i_flv, i_corr) * log (two * E_em / (xi_cut * sqrts))
            ! s4 = two * virt%c_flv(em,i_flv, i_corr) * (log (two * E_em / sqrts)**2 - log (xi_cut)**2)
            s4 = two * virt%c_flv(em,i_flv, i_corr) * & ! a**2 - b**2 = (a - b) * (a + b), for better numerical performance
                 (log (two * E_em / sqrts) + log (xi_cut)) * (log (two * E_em / sqrts) - log (xi_cut))
            s5 = two * virt%gamma_0(em,i_flv, i_corr) * log (two * E_em / sqrts)
            QB(i_flv) = QB(i_flv) + (s1 - s2 + s3 + s4 - s5) * sqme_born_virt
         end if
       end associate
       evaluated(i_flv, em) = .true.
    end do
  end subroutine virtual_compute_collinear_contribution

  module subroutine virtual_compute_massive_self_eikonals (virt, &
       i_flv, i_corr, p_born, s_over_Q2, reg_data,  QB)
    class(virtual_t), intent(inout) :: virt
    integer, intent(in) :: i_flv, i_corr
    type(vector4_t), intent(in), dimension(:) :: p_born
    real(default), intent(in) :: s_over_Q2
    type(region_data_t), intent(in) :: reg_data
    real(default), intent(inout), dimension(:) :: QB
    real(default) :: sqme_born_virt
    integer :: i
    logical :: massive
    sqme_born_virt = zero
    if (reg_data%nlo_correction_type == "EW" .and. i_corr == 1 &
       .and. qcd_ew_interferences (reg_data%flv_born(i_flv)%flst)) then
       do i = 1, size (reg_data%flv_born(i_flv)%flst)
          if (is_quark (reg_data%flv_born(i_flv)%flst (i))) then
             sqme_born_virt = -virt%sqme_color_c (i, i, i_flv)/CF
             exit
          end if
       end do
    else
       sqme_born_virt = virt%sqme_born (i_flv)
    end if
    do i = 1, virt%n_legs
       massive = reg_data%flv_born(i_flv)%massive(i)
       if (massive) then
          if (virt%c_flv (i, i_flv, i_corr) /= 0) then
             QB(i_flv) = QB(i_flv) - (virt%c_flv (i, i_flv, i_corr) &
                  * (log (s_over_Q2) - 0.5_default * I_m_eps (p_born(i)))) &
                  * sqme_born_virt
          end if
       end if
    end do
  end subroutine virtual_compute_massive_self_eikonals

  function compute_eikonal_factor (p_born, massive, i, j, s_o_Q2) result (I_ij)
    real(default) :: I_ij
    type(vector4_t), intent(in), dimension(:) :: p_born
    logical, dimension(:), intent(in) :: massive
    integer, intent(in) :: i, j
    real(default), intent(in) :: s_o_Q2
    if (massive(i) .and. massive(j)) then
       I_ij = compute_Imm (p_born(i), p_born(j), s_o_Q2)
    else if (.not. massive(i) .and. massive(j)) then
       I_ij = compute_I0m (p_born(i), p_born(j), s_o_Q2)
    else if (massive(i) .and. .not. massive(j)) then
       I_ij = compute_I0m (p_born(j), p_born(i), s_o_Q2)
    else
       I_ij = compute_I00 (p_born(i), p_born(j), s_o_Q2)
    end if
  end function compute_eikonal_factor

  function compute_I00 (pi, pj, s_o_Q2) result (I)
    type(vector4_t), intent(in) :: pi, pj
    real(default), intent(in) :: s_o_Q2
    real(default) :: I
    real(default) :: Ei, Ej
    real(default) :: pij, Eij
    real(default) :: s1, s2, s3, s4, s5
    real(default) :: arglog
    real(default), parameter :: tiny_value = epsilon(1.0)
    s1 = 0; s2 = 0; s3 = 0; s4 = 0; s5 = 0
    Ei = pi%p(0); Ej = pj%p(0)
    pij = pi * pj; Eij = Ei * Ej
    s1 = 0.5_default * log(s_o_Q2)**2
    s2 = log(s_o_Q2) * log(pij / (two * Eij))
    s3 = Li2 (pij / (two * Eij))
    s4 = 0.5_default * log (pij / (two * Eij))**2
    arglog = one - pij / (two * Eij)
    if (arglog > tiny_value) then
      s5 = log(arglog) * log(pij / (two * Eij))
    else
      s5 = zero
    end if
    I = s1 + s2 - s3 + s4 - s5
  end function compute_I00

  function compute_I0m (ki, kj, s_o_Q2) result (I)
    type(vector4_t), intent(in) :: ki, kj
    real(default), intent(in) :: s_o_Q2
    real(default) :: I
    real(default) :: logsomu
    real(default) :: s1, s2, s3
    s1 = 0; s2 = 0; s3 = 0
    logsomu = log(s_o_Q2)
    s1 = 0.5 * (0.5 * logsomu**2 - pi**2 / 6)
    s2 = 0.5 * I_0m_0 (ki, kj) * logsomu
    s3 = 0.5 * I_0m_eps (ki, kj)
    I = s1 + s2 - s3
  end function compute_I0m

  function compute_Imm (pi, pj, s_o_Q2) result (I)
    type(vector4_t), intent(in) :: pi, pj
    real(default), intent(in) :: s_o_Q2
    real(default) :: I
    real(default) :: s1, s2
    s1 = 0.5 * log(s_o_Q2) * I_mm_0(pi, pj)
    s2 = 0.5 * I_mm_eps(pi, pj)
    I = s1 - s2
  end function compute_Imm

  function I_m_eps (p) result (I)
    type(vector4_t), intent(in) :: p
    real(default) :: I
    real(default) :: beta
    beta = space_part_norm (p)/p%p(0)
    if (beta < tiny_07) then
       I = four * (one + beta**2/3 + beta**4/5 + beta**6/7)
    else
       I = two * log((one + beta) / (one - beta)) / beta
    end if
  end function I_m_eps

  function I_0m_eps (p, k) result (I)
    type(vector4_t), intent(in) :: p, k
    real(default) :: I
    type(vector4_t) :: pp, kp
    real(default) :: beta

    pp = p / p%p(0); kp = k / k%p(0)

    beta = sqrt (one - kp*kp)
    I = -2*(log((one - beta) / (one + beta))**2/4 + log((pp*kp) / (one + beta))*log((pp*kp) / (one - beta)) &
        + Li2(one - (pp*kp) / (one + beta)) + Li2(one - (pp*kp) / (one - beta)))
  end function I_0m_eps

  function I_0m_0 (p, k) result (I)
    type(vector4_t), intent(in) :: p, k
    real(default) :: I
    type(vector4_t) :: pp, kp

    pp = p / p%p(0); kp = k / k%p(0)
    I = log((pp*kp)**2 / kp**2)
  end function I_0m_0

  function I_mm_eps (p1, p2) result (I)
    type(vector4_t), intent(in) :: p1, p2
    real(default) :: I
    type(vector4_t) :: q1, q2
    type(vector3_t) :: beta1, beta2
    real(default) :: a, b
    real(default) :: zp, zm, z1, z2, x1, x2
    real(default) :: zmb, z1b
    real(default) :: K1, K2, b1, b2
    real(default) :: nu, a_kl, m12, m22
    beta1 = space_part (p1) / energy(p1)
    beta2 = space_part (p2) / energy(p2)
    if (min (one - beta1**1, one - beta2**1) < tiny_07) then
       if (beta1**1 < beta2**1) then
          call switch_beta (beta1, beta2)
          q1 = p2
          q2 = p1
       else
          q1 = p1
          q2 = p2
       end if
       b1 = beta1**1
       b2 = beta2**1
       m12 = q1**2
       m22 = q2**2
       a_kl = ((q1*q2) + sqrt((q1*q2)**2 - q1**2*q2**2))/m12
       nu = (a_kl**2 * m12 - m22) / two / (a_kl * p1%p(0) - q2%p(0))
       K1 =   0.5_default * log ((one - b1) / (one + b1))**2 &
            + two * Li2 (one - (one - b1)*(a_kl*q1%p(0)/nu)) &
            + two * Li2 (one - (one + b1)*(a_kl*q1%p(0)/nu))
       K2 =   0.5_default * log((one - b2) / (one + b2))**2 &
            + two * Li2 (one - (one - b2)*(q2%p(0)/nu)) &
            + two * Li2 (one - (one + b2)*(q2%p(0)/nu))
       I = (K2 - K1) / sqrt(one - m12*m22/(q1*q2)**2)
    else
       a = beta1**2 + beta2**2 - 2 * beta1 * beta2
       b = beta1**2 * beta2**2 - (beta1 * beta2)**2
       if (beta1**1 > beta2**1) call switch_beta (beta1, beta2)
       b2 = beta2**1
       if (beta1 == vector3_null) then
          I = (-0.5 * log ((one - b2) / (one + b2))**2 - two * Li2 (-two * b2 / (one - b2))) &
               * one / sqrt (a - b)
          return
       end if
       x1 = beta1**2 - beta1 * beta2
       x2 = beta2**2 - beta1 * beta2
       zp = sqrt (a) + sqrt (a - b)
       zm = sqrt (a) - sqrt (a - b)
       zmb = one  / zp
       z1 = sqrt (x1**2 + b) - x1
       z2 = sqrt (x2**2 + b) + x2
       z1b = one / (sqrt (x1**2 + b) + x1)
       K1 = - 0.5_default * log (((z1b - zmb) * (zp - z1)) / ((zp + z1) * (z1b + zmb)))**2 &
            - two * Li2 ((two * zmb * (zp - z1)) / ((zp - zm) * (zmb + z1b))) &
            - two * Li2 ((-two * zp * (zm + z1)) / ((zp - zm) * (zp - z1)))
       K2 = - 0.5_default * log ((( z2 - zm) * (zp - z2)) / ((zp + z2) * (z2 + zm)))**2 &
            - two * Li2 ((two * zm * (zp - z2)) / ((zp - zm) * (zm + z2))) &
            - two * Li2 ((-two * zp * (zm + z2)) / ((zp - zm) * (zp - z2)))
       I = (K2 - K1) * (one - beta1 * beta2) / sqrt (a - b)
    end if
  contains
    subroutine switch_beta (beta1, beta2)
      type(vector3_t), intent(inout) :: beta1, beta2
      type(vector3_t) :: beta_tmp
      beta_tmp = beta1
      beta1 = beta2
      beta2 = beta_tmp
    end subroutine switch_beta
  end function I_mm_eps

  function I_mm_0 (k1, k2) result (I)
    type(vector4_t), intent(in) :: k1, k2
    real(default) :: I
    real(default) :: beta, kquotient
    kquotient = k1**2 * k2**2 / (k1 * k2)**2
    if (kquotient > tiny_13) then
       beta = sqrt (one - kquotient)
       I = log ((one + beta) / (one - beta)) / beta
    else
       beta = one - kquotient / two
       I = log (two * (one + beta) / kquotient) / beta
    end if
  end function I_mm_0

  module subroutine virtual_final (virtual)
    class(virtual_t), intent(inout) :: virtual
    if (allocated (virtual%gamma_0)) deallocate (virtual%gamma_0)
    if (allocated (virtual%gamma_p)) deallocate (virtual%gamma_p)
    if (allocated (virtual%c_flv)) deallocate (virtual%c_flv)
    if (allocated (virtual%n_is_neutrinos)) deallocate (virtual%n_is_neutrinos)
  end subroutine virtual_final


end submodule virtual_s

