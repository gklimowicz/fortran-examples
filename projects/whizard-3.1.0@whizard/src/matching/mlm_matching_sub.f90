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

submodule (mlm_matching) mlm_matching_s

  use debug_master, only: debug_on
  use io_units
  use format_utils, only: write_separator
  use diagnostics
  use file_utils
  use subevents, only: PRT_OUTGOING
  use shower_base
  use ktclus

  implicit none

contains

  module subroutine mlm_matching_settings_init (settings, var_list)
    class(mlm_matching_settings_t), intent(out) :: settings
    type(var_list_t), intent(in) :: var_list
    settings%mlm_Qcut_ME = &
         var_list%get_rval (var_str ("mlm_Qcut_ME"))
    settings%mlm_Qcut_PS = &
         var_list%get_rval (var_str ("mlm_Qcut_PS"))
    settings%mlm_ptmin = &
         var_list%get_rval (var_str ("mlm_ptmin"))
    settings%mlm_etamax = &
         var_list%get_rval (var_str ("mlm_etamax"))
    settings%mlm_Rmin = &
         var_list%get_rval (var_str ("mlm_Rmin"))
    settings%mlm_Emin = &
         var_list%get_rval (var_str ("mlm_Emin"))
    settings%mlm_nmaxMEjets = &
         var_list%get_ival (var_str ("mlm_nmaxMEjets"))

    settings%mlm_ETclusfactor = &
         var_list%get_rval (var_str ("mlm_ETclusfactor"))
    settings%mlm_ETclusminE = &
         var_list%get_rval (var_str ("mlm_ETclusminE"))
    settings%mlm_etaclusfactor = &
         var_list%get_rval (var_str ("mlm_etaclusfactor"))
    settings%mlm_Rclusfactor = &
         var_list%get_rval (var_str ("mlm_Rclusfactor"))
    settings%mlm_Eclusfactor = &
         var_list%get_rval (var_str ("mlm_Eclusfactor"))
  end subroutine mlm_matching_settings_init

  module subroutine mlm_matching_settings_write (settings, unit)
    class(mlm_matching_settings_t), intent(in) :: settings
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(3x,A,ES19.12)") &
         "mlm_Qcut_ME                  = ", settings%mlm_Qcut_ME
    write (u, "(3x,A,ES19.12)") &
         "mlm_Qcut_PS                  = ", settings%mlm_Qcut_PS
    write (u, "(3x,A,ES19.12)") &
         "mlm_ptmin                    = ", settings%mlm_ptmin
    write (u, "(3x,A,ES19.12)") &
         "mlm_etamax                   = ", settings%mlm_etamax
    write (u, "(3x,A,ES19.12)") &
         "mlm_Rmin                     = ", settings%mlm_Rmin
    write (u, "(3x,A,ES19.12)") &
         "mlm_Emin                     = ", settings%mlm_Emin
    write (u, "(3x,A,1x,I0)") &
         "mlm_nmaxMEjets               = ", settings%mlm_nmaxMEjets
    write (u, "(3x,A,ES19.12)") &
         "mlm_ETclusfactor  (D=0.2)    = ", settings%mlm_ETclusfactor
    write (u, "(3x,A,ES19.12)") &
         "mlm_ETclusminE    (D=5.0)    = ", settings%mlm_ETclusminE
    write (u, "(3x,A,ES19.12)") &
         "mlm_etaclusfactor (D=1.0)    = ", settings%mlm_etaClusfactor
    write (u, "(3x,A,ES19.12)") &
         "mlm_Rclusfactor   (D=1.0)    = ", settings%mlm_RClusfactor
    write (u, "(3x,A,ES19.12)") &
         "mlm_Eclusfactor   (D=1.0)    = ", settings%mlm_EClusfactor
  end subroutine mlm_matching_settings_write

  module subroutine mlm_matching_init (matching, var_list, process_name)
    class(mlm_matching_t), intent(out) :: matching
    type(var_list_t), intent(in) :: var_list
    type(string_t), intent(in) :: process_name
    if (debug_on) call msg_debug (D_MATCHING, "matching_init")
    call matching%settings%init (var_list)
    matching%process_name = process_name
  end subroutine mlm_matching_init

  module subroutine mlm_matching_write (matching, unit)
    class(mlm_matching_t), intent(in) :: matching
    integer, intent(in), optional :: unit
    integer :: i, u
    u = given_output_unit (unit);  if (u < 0)  return

    write (u, "(1x,A)") "MLM matching:"
    call matching%settings%write (u)
    write (u, "(3x,A)") "Momenta of ME partons:"
    if (allocated (matching%P_ME)) then
       do i = 1, size (matching%P_ME)
          write (u, "(4x)", advance = "no")
          call vector4_write (matching%P_ME(i), unit = u)
       end do
    else
       write (u, "(5x,A)")  "[empty]"
    end if
    call write_separator (u)
    write (u, "(3x,A)")  "Momenta of ME jets:"
    if (allocated (matching%JETS_ME)) then
       do i = 1, size (matching%JETS_ME)
          write (u, "(4x)", advance = "no")
          call vector4_write (matching%JETS_ME(i), unit = u)
       end do
    else
       write (u, "(5x,A)")  "[empty]"
    end if
    call write_separator (u)
    write(u, "(3x,A)")  "Momenta of shower partons:"
    if (allocated (matching%P_PS)) then
       do i = 1, size (matching%P_PS)
          write (u, "(4x)", advance = "no")
          call vector4_write (matching%P_PS(i), unit = u)
       end do
    else
       write (u, "(5x,A)")  "[empty]"
    end if
    call write_separator (u)
    write (u, "(3x,A)")  "Momenta of shower jets:"
    if (allocated (matching%JETS_PS)) then
       do i = 1, size (matching%JETS_PS)
          write (u, "(4x)", advance = "no")
          call vector4_write (matching%JETS_PS(i), unit = u)
       end do
    else
       write (u, "(5x,A)")  "[empty]"
    end if
    call write_separator (u)
  end subroutine mlm_matching_write

  module function mlm_matching_get_method (matching) result (method)
     type(string_t) :: method
     class(mlm_matching_t), intent(in) :: matching
     method = matching_method (MATCH_MLM)
  end function mlm_matching_get_method

  module subroutine mlm_matching_before_shower &
       (matching, particle_set, vetoed)
    class(mlm_matching_t), intent(inout) :: matching
    type(particle_set_t), intent(inout) :: particle_set
    logical, intent(out) :: vetoed
    vetoed = .false.
  end subroutine mlm_matching_before_shower

  module subroutine mlm_matching_after_shower (matching, particle_set, vetoed)
    class(mlm_matching_t), intent(inout) :: matching
    type(particle_set_t), intent(inout) :: particle_set
    logical, intent(out) :: vetoed
    if (debug_on) call msg_debug (D_MATCHING, "mlm_matching_after_shower")
    call matching%shower%get_final_colored_ME_momenta (matching%P_ME)
    call matching%fill_P_PS (particle_set)
    !!! MLM stage 3 -> reconstruct and possibly reject
    call matching%apply (vetoed)
    if (debug_active (D_MATCHING)) call matching%write ()
    if (allocated (matching%P_ME))  deallocate (matching%P_ME)
    if (allocated (matching%P_PS))  deallocate (matching%P_PS)
    if (allocated (matching%JETS_ME))  deallocate (matching%JETS_ME)
    if (allocated (matching%JETS_PS))  deallocate (matching%JETS_PS)
  end subroutine mlm_matching_after_shower

  module subroutine mlm_matching_fill_P_PS (matching, particle_set)
    class(mlm_matching_t), intent(inout) :: matching
    type(particle_set_t), intent(in) :: particle_set
    integer :: i, j, n_jets_PS
    integer, dimension(2) :: col
    type(particle_t) :: tempprt
    real(double) :: eta
    type(vector4_t) :: p_tmp

    !!! loop over particles and extract final colored ones with eta<etamax
    n_jets_PS = 0
    do i = 1, particle_set%get_n_tot ()
       if (signal_is_pending ()) return
       tempprt = particle_set%get_particle (i)
       if (tempprt%get_status () /= PRT_OUTGOING) cycle
       col = tempprt%get_color ()
       if (all (col == 0)) cycle
! TODO: (bcn 2015-04-28) where is the corresponding part for lepton colliders?
       if (matching%is_hadron_collision) then
          p_tmp = tempprt%get_momentum ()
          if (energy (p_tmp) - longitudinal_part (p_tmp) < 1.E-10_default .or. &
               energy (p_tmp) + longitudinal_part (p_tmp) < 1.E-10_default) then
             eta = pseudorapidity (p_tmp)
          else
             eta = rapidity (p_tmp)
          end if
          if (eta > matching%settings%mlm_etaClusfactor * &
               matching%settings%mlm_etamax)  then
             if (debug_active (D_MATCHING)) then
                call msg_debug (D_MATCHING, "Rejecting this particle")
                call tempprt%write ()
             end if
             cycle
          end if
       end if
       n_jets_PS = n_jets_PS + 1
    end do

    allocate (matching%P_PS(1:n_jets_PS))
    if (debug_on) call msg_debug (D_MATCHING, "n_jets_ps", n_jets_ps)

    j = 1
    do i = 1, particle_set%get_n_tot ()
       tempprt = particle_set%get_particle (i)
       if (tempprt%get_status () /= PRT_OUTGOING) cycle
       col = tempprt%get_color ()
       if (all(col == 0)) cycle
! TODO: (bcn 2015-04-28) where is the corresponding part for lepton colliders?
       if (matching%is_hadron_collision) then
          p_tmp = tempprt%get_momentum ()
          if (energy (p_tmp) - longitudinal_part (p_tmp) < 1.E-10_default .or. &
               energy (p_tmp) + longitudinal_part (p_tmp) < 1.E-10_default) then
             eta = pseudorapidity (p_tmp)
          else
             eta = rapidity (p_tmp)
          end if
          if (eta > matching%settings%mlm_etaClusfactor * &
               matching%settings%mlm_etamax) cycle
       end if
       matching%P_PS(j) = tempprt%get_momentum ()
       j = j + 1
    end do
  end subroutine mlm_matching_fill_P_PS

  module subroutine mlm_matching_apply (matching, vetoed)
    class(mlm_matching_t), intent(inout) :: matching
    logical, intent(out) :: vetoed
    integer :: i, j
    integer :: n_jets_ME, n_jets_PS, n_jets_PS_atycut
    real(double) :: ycut
    real(double), dimension(:, :), allocatable :: PP
    real(double), dimension(:), allocatable :: Y
    real(double), dimension(:,:), allocatable :: P_JETS
    real(double), dimension(:,:), allocatable :: P_ME
    integer, dimension(:), allocatable :: JET
    integer :: NJET, NSUB
    integer :: imode
!!! TODO: (bcn 2014-03-26) Why is ECUT hard coded to 1?
!!! It is the denominator of the KT measure. Candidate for removal
    real(double) :: ECUT = 1._double
    integer :: ip1,ip2

    ! KTCLUS COMMON BLOCK
    INTEGER NMAX,NUM,HIST
    PARAMETER (NMAX=512)
    DOUBLE PRECISION P,KT,KTP,KTS,ETOT,RSQ,KTLAST
    COMMON /KTCOMM/ETOT,RSQ,P(9,NMAX),KTP(NMAX,NMAX),KTS(NMAX), &
         KT(NMAX),KTLAST(NMAX),HIST(NMAX),NUM

    vetoed = .true.
    if (signal_is_pending ())  return

    if (allocated (matching%P_ME)) then
       ! print *, "number of partons after ME: ", size(matching%P_ME)
       n_jets_ME = size (matching%P_ME)
    else
       n_jets_ME = 0
    end if
    if (allocated (matching%p_PS)) then
       ! print *, "number of partons after PS: ", size(matching%p_PS)
       n_jets_PS = size (matching%p_PS)
    else
       n_jets_PS = 0
    end if

    if (n_jets_ME > 0) then
       ycut = (matching%settings%mlm_ptmin)**2
       allocate (PP(1:4, 1:N_jets_ME))
       do i = 1, n_jets_ME
          PP(1:3,i) = matching%p_ME(i)%p(1:3)
          PP(4,i) = matching%p_ME(i)%p(0)
       end do

       if (matching%is_hadron_collision) then
          imode = matching%settings%kt_imode_hadronic
       else
          imode = matching%settings%kt_imode_leptonic
       end if

       allocate (P_ME(1:4,1:n_jets_ME))
       allocate (JET(1:n_jets_ME))
       allocate (Y(1:n_jets_ME))

       if (signal_is_pending ())  return
       call KTCLUR (imode, PP, n_jets_ME, &
            dble (matching%settings%mlm_Rclusfactor * matching%settings%mlm_Rmin), ECUT, y, *999)
       call KTRECO (1, PP, n_jets_ME, ECUT, ycut, ycut, P_ME, JET, &
            NJET, NSUB, *999)

       n_jets_ME = NJET
       if (NJET > 0) then
          allocate (matching%JETS_ME (1:NJET))
          do i = 1, NJET
             matching%JETS_ME(i) = vector4_moving (REAL(P_ME(4,i), default), &
                  vector3_moving([REAL(P_ME(1,i), default), &
                  REAL(P_ME(2,i), default), REAL(P_ME(3,i), default)]))
          end do
       end if
       deallocate (P_ME)
       deallocate (JET)
       deallocate (Y)
       deallocate (PP)
    end if

    if (n_jets_PS > 0) then
       ycut = (matching%settings%mlm_ptmin + max (matching%settings%mlm_ETclusminE, &
            matching%settings%mlm_ETclusfactor * matching%settings%mlm_ptmin))**2
       allocate (PP(1:4, 1:n_jets_PS))
       do i = 1, n_jets_PS
          PP(1:3,i) = matching%p_PS(i)%p(1:3)
          PP(4,i) = matching%p_PS(i)%p(0)
       end do

       if (matching%is_hadron_collision) then
          imode = matching%settings%kt_imode_hadronic
       else
          imode = matching%settings%kt_imode_leptonic
       end if

       allocate (P_JETS(1:4,1:n_jets_PS))
       allocate (JET(1:n_jets_PS))
       allocate (Y(1:n_jets_PS))

       if (signal_is_pending ()) return
       call KTCLUR (imode, PP, n_jets_PS, &
            dble (matching%settings%mlm_Rclusfactor * matching%settings%mlm_Rmin), &
            ECUT, y, *999)
       call KTRECO (1, PP, n_jets_PS, ECUT, ycut, ycut, P_JETS, JET, &
            NJET, NSUB, *999)
       n_jets_PS_atycut = NJET
       if (n_jets_ME == matching%settings%mlm_nmaxMEjets .and. NJET > 0) then
          ! print *, " resetting ycut to ", Y(matching%settings%mlm_nmaxMEjets)
          ycut = y(matching%settings%mlm_nmaxMEjets)
          call KTRECO (1, PP, n_jets_PS, ECUT, ycut, ycut, P_JETS, JET, &
               NJET, NSUB, *999)
       end if

       ! !Sample of code for a FastJet interface
       ! palg = 1d0         ! 1.0d0 = kt, 0.0d0 = Cam/Aachen, -1.0d0 = anti-kt
       ! R = 0.7_double     ! radius parameter
       ! f = 0.75_double    ! overlap threshold
       ! !call fastjetppgenkt(PP,n,R,palg,P_JETS,NJET)   ! KT-Algorithm
       ! !call fastjetsiscone(PP,n,R,f,P_JETS,NJET)      ! SiSCone-Algorithm

       if (NJET > 0) then
          allocate (matching%JETS_PS(1:NJET))
          do i = 1, NJET
             matching%JETS_PS(i) = vector4_moving (REAL(P_JETS(4,i), default), &
                  vector3_moving([REAL(P_JETS(1,i), default), &
                  REAL(P_JETS(2,i), default), REAL(P_JETS(3,i), default)]))
          end do
       end if

       deallocate (P_JETS)
       deallocate (JET)
       deallocate (Y)
    else
       n_jets_PS_atycut = 0
    end if

    if (n_jets_PS_atycut < n_jets_ME) then
       ! print *, "DISCARDING: Not enough PS jets: ", n_jets_PS_atycut
       return
    end if
    if (n_jets_PS_atycut > n_jets_ME .and. n_jets_ME /= matching%settings%mlm_nmaxMEjets) then
       ! print *, "DISCARDING: Too many PS jets: ", n_jets_PS_atycut
       return
    end if

    if (allocated(matching%JETS_PS)) then
       ! print *, "number of jets after PS: ", size(matching%JETS_PS)
       n_jets_PS = size (matching%JETS_PS)
    else
       n_jets_PS = 0
    end if
    if (n_jets_ME > 0 .and. n_jets_PS > 0) then
       n_jets_PS = size (matching%JETS_PS)
       if (allocated (PP))  deallocate(PP)
       allocate (PP(1:4, 1:n_jets_PS + 1))
       do i = 1, n_jets_PS
          if (signal_is_pending ()) return
          PP(1:3,i) = matching%JETS_PS(i)%p(1:3)
          PP(4,i) = matching%JETS_PS(i)%p(0)
       end do
       if (allocated (Y))  deallocate(Y)
       allocate (Y(1:n_jets_PS + 1))
       y = zero
       do i = 1, n_jets_ME
          PP(1:3,n_jets_PS + 2 - i) = matching%JETS_ME(i)%p(1:3)
          PP(4,n_jets_PS + 2 - i) = matching%JETS_ME(i)%p(0)
          !!! This makes more sense than hardcoding
          ! call KTCLUS (4313, PP, (n_jets_PS + 2 - i), 1.0_double, Y, *999)
          call KTCLUR (imode, PP, (n_jets_PS + 2 - i), &
            dble (matching%settings%mlm_Rclusfactor * matching%settings%mlm_Rmin), &
            ECUT, y, *999)
          if (0.99 * y(n_jets_PS + 1 - (i - 1)).gt.ycut) then
             ! print *, "DISCARDING: Jet ", i, " not clusterd"
             return
          end if
          !!! search for and remove PS jet clustered with ME Jet
          ip1 = HIST(n_jets_PS + 2 - i) / NMAX
          ip2 = mod(hist(n_jets_PS + 2 - i), NMAX)
          if ((ip2 /= n_jets_PS + 2 - i) .or. (ip1 <= 0)) then
             ! print *, "DISCARDING: Jet ", i, " not clustered ", ip1, ip2, &
             !      hist(n_jets_PS + 2 - i)
             return
          else
             ! print *, "PARTON clustered", ip1, ip2, hist(n_jets_PS + 2 - i)
             PP(:,IP1) = zero
             do j = IP1, n_jets_PS - i
                PP(:, j) = PP(:,j + 1)
             end do
          end if
       end do
    end if

    vetoed = .false.
999 continue
  end subroutine mlm_matching_apply


end submodule mlm_matching_s

