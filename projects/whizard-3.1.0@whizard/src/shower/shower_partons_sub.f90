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

submodule (shower_partons) shower_partons_s

  use debug_master, only: debug_on
  use io_units
  use system_defs, only: TAB
  use diagnostics
  use physics_defs
  use sm_physics
  use colors
  use flavors
  use subevents

  implicit none

contains

  module function parton_to_particle &
       (parton, model, from_hard_int) result (particle)
    type(particle_t) :: particle
    class(parton_t), intent(in) :: parton
    class(model_data_t), pointer, intent(in) :: model
    logical, intent(in), optional :: from_hard_int
    integer :: col, anti_col
    call parton%to_color (col, anti_col, from_hard_int)
    call particle%init (parton%to_status (from_hard_int), parton%type, &
         model, col, anti_col, parton%momentum)
  end function parton_to_particle

  pure module function parton_of_particle (particle, nr) result (parton)
    type(parton_t) :: parton
    type(particle_t), intent(in) :: particle
    integer, intent(in) :: nr
    integer, dimension(2) :: col_array
    parton%nr = nr
    parton%momentum = particle%p
    parton%t = particle%p2
    parton%type = particle%flv%get_pdg ()
    col_array = particle%get_color ()
    parton%c1 = col_array (1)
    parton%c2 = col_array (2)
    parton%interactionnr = 1
    parton%mass2 = particle%flv%get_mass () ** 2
  end function parton_of_particle

  pure module function parton_to_status &
       (parton, from_hard_int) result (status)
    integer :: status
    class(parton_t), intent(in) :: parton
    logical, intent(in), optional :: from_hard_int
    logical :: fhi
    fhi = .false.; if (present (from_hard_int))  fhi = from_hard_int
    if (fhi .or. parton%is_colored ()) then
       if (associated (parton%initial) .and. .not. parton%belongstoFSR) then
          status = PRT_INCOMING
       else
          status = PRT_OUTGOING
       end if
    else
       status = PRT_BEAM_REMNANT
    end if
  end function parton_to_status

  pure module subroutine parton_to_color (parton, c1, c2, from_hard_int)
    class(parton_t), intent(in) :: parton
    integer, intent(out) :: c1, c2
    logical, intent(in), optional :: from_hard_int
    logical :: fhi
    fhi = .false.; if (present (from_hard_int))  fhi = from_hard_int
    c1 = 0
    c2 = 0
    if (parton%is_colored ()) then
       if (fhi) then
          if (parton%c1 /= 0) c1 = parton%c1
          if (parton%c2 /= 0) c2 = parton%c2
       else
          if (parton%c1 /= 0) c1 = 500 + parton%c1
          if (parton%c2 /= 0) c2 = 500 + parton%c2
       end if
    end if
  end subroutine parton_to_color

  module subroutine parton_copy (prt1, prt2)
    type(parton_t), intent(in) :: prt1
    type(parton_t), intent(out) :: prt2
    if (associated (prt1%settings))  prt2%settings => prt1%settings
    prt2%nr = prt1%nr
    prt2%type = prt1%type
    prt2%momentum = prt1%momentum
    prt2%t = prt1%t
    prt2%mass2 = prt1%mass2
    prt2%scale = prt1%scale
    prt2%z = prt1%z
    prt2%costheta = prt1%costheta
    prt2%x = prt1%x
    prt2%simulated = prt1%simulated
    prt2%belongstoFSR = prt1%belongstoFSR
    prt2%belongstointeraction = prt1%belongstointeraction
    prt2%interactionnr = prt1%interactionnr
    if (associated (prt1%parent))  prt2%parent  => prt1%parent
    if (associated (prt1%child1))  prt2%child1  => prt1%child1
    if (associated (prt1%child2))  prt2%child2  => prt1%child2
    if (associated (prt1%initial)) prt2%initial => prt1%initial
    prt2%c1 = prt1%c1
    prt2%c2 = prt1%c2
    prt2%aux_pt = prt1%aux_pt
  end subroutine parton_copy

  elemental module function parton_get_costheta (prt) result (costheta)
    class(parton_t), intent(in) :: prt
    real(default) :: costheta
    real(default) :: denom
    denom = two * prt%z * (one - prt%z) * prt%momentum%p(0)**2
    if (denom > eps0) then
       costheta = one - prt%t / denom
    else
       costheta = - one
    end if
  end function parton_get_costheta

  elemental module function parton_get_costheta_mass (prt) result (costheta)
    class(parton_t), intent(in) :: prt
    real(default) :: costheta
    real(default) :: sqrt12
    if (prt%is_branched ()) then
       if (prt%child1%simulated .and. &
           prt%child2%simulated) then
          sqrt12 = sqrt (max (zero, (prt%z)**2 * prt%momentum%p(0)**2 &
                                    - prt%child1%t)) * &
                   sqrt (max (zero, (one - prt%z)**2 * prt%momentum%p(0)**2 &
                                    - prt%child2%t))
          if (sqrt12 > eps0) then
             costheta = (prt%t - prt%child1%t - prt%child2%t - &
                  two * prt%z * (one - prt%z) * prt%momentum%p(0)**2) / &
                  (- two * sqrt12)
             return
          end if
       end if
    end if
    costheta = prt%get_costheta ()
  end function parton_get_costheta_mass

  elemental module function parton_get_costheta_motherfirst &
       (prt) result (costheta)
    class(parton_t), intent(in) :: prt
    real(default) :: costheta
    if (prt%is_branched ()) then
       if ((prt%child1%simulated .or. &
            prt%child1%is_final () .or. &
            prt%child1%is_branched ()) .and. &
            (prt%child2%simulated .or. &
             prt%child2%is_final () .or. &
             prt%child2%is_branched ())) then
          costheta = enclosed_angle_ct (prt%momentum, prt%child1%momentum)
          return
       end if
    end if
    costheta = - two
  end function parton_get_costheta_motherfirst

  pure function get_beta (t,E) result (beta)
    real(default), intent(in) :: t,E
    real(default) :: beta
    beta = sqrt (max (tiny_07, one - t /(E**2)))
  end function get_beta

  elemental module function parton_get_beta (prt) result (beta)
    class(parton_t), intent(in) :: prt
    real(default) :: beta
    beta = sqrt (max (tiny_07, one - prt%t / prt%momentum%p(0)**2))
  end function parton_get_beta

  module subroutine parton_write (prt, unit)
    class(parton_t), intent(in) :: prt
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit); if (u < 0) return

    write (u, "(1x,7A)") "Shower parton <nr>", TAB, "<type>", TAB // TAB, &
         "<parent>", TAB, "<mom(0:3)>"
    write (u, "(2x,I5,3A)", advance = "no")  prt%nr, TAB, TAB, TAB
    if (prt%is_final ()) then
       write (u, "(1x,I5,1x,A)", advance = "no") prt%type, TAB // TAB
    else
       write (u, "('[',I5,']',A)", advance = "no") prt%type, TAB // TAB
    end if
    if (associated (prt%parent)) then
       write (u, "(I5,A)", advance = "no") prt%parent%nr, TAB // TAB
    else
       write (u, "(5x,2A)", advance = "no") TAB, TAB
    end if
    write (u, "(4(ES12.5,A))") prt%momentum%p(0), TAB, &
                               prt%momentum%p(1), TAB, &
                               prt%momentum%p(2), TAB, &
                               prt%momentum%p(3)
    write (u, "(1x,9A)") "<p4square>", TAB // TAB, "<t>", TAB // TAB, &
         "<scale>", TAB // TAB, "<c1>", TAB, "<c2>", TAB, "<mass2>"
    write (u, "(1x,3(ES12.5,A))", advance = "no") &
         prt%momentum ** 2, TAB // TAB, prt%t, TAB, prt%scale, TAB, prt%mass2
    write (u, "(2(I4,A))") prt%c1, TAB, prt%c2, TAB
    if (prt%is_branched ()) then
       if (prt%belongstoFSR) then
          write (u, "(1x,9A)") "costheta(prt)", TAB, &
               "costheta_correct(prt)", TAB, &
               "prt%costheta", TAB, "prt%z", TAB, &
               "costheta_motherfirst(prt)"
          write (u, "(1X,5(ES12.5,A))") &
               prt%get_costheta (), TAB, &
               prt%get_costheta_mass (), TAB // TAB, &
               prt%costheta, TAB, prt%z, TAB, &
               prt%get_costheta_motherfirst (), TAB
       else
          write (u, "(1x,9A)") "prt%z", TAB, "prt%x", TAB, &
               "costheta_correct(prt)", TAB, &
               "prt%costheta", TAB, &
               "costheta_motherfirst(prt)"
          write (u, "(1X,5(ES12.5,A))") &
               prt%z, TAB, prt%x, TAB, &
               prt%get_costheta_mass (), TAB, &
               prt%costheta, TAB, &
               prt%get_costheta_motherfirst (), TAB
       end if
    else
       if (prt%belongstoFSR) then
          write (u, "(1X,A)") "not branched."
       else
          write (u, "(1X,A,ES12.5)") "not branched. x = ",  prt%x
       end if
    end if
    write (u, "(A)", advance = "no") " Parton"
    if (prt%belongstoFSR) then
       write (u, "(A)", advance = "no")  " is FSR,"
    else
       if (associated (prt%initial)) then
          write (u, "(A,I1)", advance = "no")  " from hadron,", prt%initial%nr
       else
          write (u, "(A)", advance = "no")  ""
       end if
    end if
    if (prt%is_final ()) then
       write (u, "(A)", advance = "no")  " is final,"
    else
       write (u, "(A)", advance = "no")  ""
    end if
    if (prt%simulated) then
       write (u, "(A)", advance = "no")  " is simulated,"
    else
       write (u, "(A)", advance = "no")  ""
    end if
    if (associated (prt%child1) .and. associated (prt%child2)) then
       write (u, "(A,2(I5),A)", advance = "no") &
            " has children: ", prt%child1%nr, prt%child2%nr, ","
    else if (associated (prt%child1)) then
       write (u, "(A,1(I5),A)", advance = "no") &
            " has one child: ", prt%child1%nr, ", "
    end if
    if (prt%belongstointeraction) then
       write (u, "(A,I2)") " belongs to interaction ", &
            prt%interactionnr
    else
       write (u, "(A,I2)") " does not belong to interaction ", &
            prt%interactionnr
    end if
    write (u,"(A)")  TAB
  end subroutine parton_write

  elemental module function parton_is_final (prt) result (is_final)
    class(parton_t), intent(in) :: prt
    logical :: is_final
    is_final = .false.
    if (prt%belongstoFSR) then
       is_final = .not. associated (prt%child1) .and. &
            (.not. prt%belongstointeraction .or. &
            (prt%belongstointeraction .and. prt%simulated))
    end if
  end function parton_is_final

  elemental module function parton_is_branched (prt) result (is_branched)
    class(parton_t), intent(in) :: prt
    logical :: is_branched
    is_branched = associated (prt%child1) .and. associated (prt%child2)
  end function parton_is_branched

  pure module subroutine parton_set_simulated (prt, sim)
    class(parton_t), intent(inout) :: prt
    logical, intent(in), optional :: sim
    if (present (sim)) then
       prt%simulated = sim
    else
       prt%simulated = .true.
    end if
  end subroutine parton_set_simulated

  module subroutine parton_set_parent (prt, parent)
    type(parton_t), intent(inout) :: prt
    type(parton_t), intent(in) , target :: parent
    prt%parent => parent
  end subroutine parton_set_parent

  module function parton_get_parent (prt) result (parent)
    type(parton_t), intent(in) :: prt
    type(parton_t), pointer :: parent
    parent => prt%parent
  end function parton_get_parent

  module subroutine parton_set_initial (prt, initial)
    type(parton_t), intent(inout) :: prt
    type(parton_t), intent(in) , target :: initial
    prt%initial => initial
  end subroutine parton_set_initial

  module function parton_get_initial (prt) result (initial)
    type(parton_t), intent(in) :: prt
    type(parton_t), pointer :: initial
    initial => prt%initial
  end function parton_get_initial

  module subroutine parton_set_child (prt, child, i)
    type(parton_t), intent(inout) :: prt
    type(parton_t), intent(in), target :: child
    integer, intent(in) ::  i
    if (i == 1) then
       prt%child1 => child
    else
       prt%child2 => child
    end if
  end subroutine parton_set_child

  module function parton_get_child (prt, i) result (child)
    type(parton_t), pointer :: child
    type(parton_t), intent(in) :: prt
    integer, intent(in) :: i
    child => null ()
    if (i == 1) then
       child => prt%child1
    else
       child => prt%child2
    end if
  end function parton_get_child

  elemental module function parton_is_quark (prt) result (is_quark)
    class(parton_t), intent(in) ::prt
    logical :: is_quark
    is_quark = abs (prt%type) <= 6 .and. prt%type /= 0
  end function parton_is_quark

  elemental module function parton_is_squark (prt) result (is_squark)
    class(parton_t), intent(in) ::prt
    logical :: is_squark
    is_squark = ((abs(prt%type) >= 1000001) .and. (abs(prt%type) <= 1000006)) &
             .or. ((abs(prt%type) >= 2000001) .and. (abs(prt%type) <= 2000006))
  end function parton_is_squark

  elemental module function parton_is_gluon (prt) result (is_gluon)
    class(parton_t), intent(in) :: prt
    logical :: is_gluon
    is_gluon = prt%type == GLUON .or. prt%type == 9
  end function parton_is_gluon

  elemental module function parton_is_gluino (prt) result (is_gluino)
    class(parton_t), intent(in) :: prt
    logical :: is_gluino
    is_gluino = prt%type == 1000021
  end function parton_is_gluino

  elemental module function parton_is_proton (prt) result (is_hadron)
    class(parton_t), intent(in) :: prt
    logical :: is_hadron
    is_hadron = abs (prt%type) == PROTON
  end function parton_is_proton

  pure module function parton_is_colored (parton) result (is_colored)
    logical :: is_colored
    class(parton_t), intent(in) :: parton
    is_colored = parton_is_quark (parton) .or. parton_is_gluon (parton)
  end function parton_is_colored

  elemental module function parton_mass (prt) result (mass)
    class(parton_t), intent(in) :: prt
    real(default) :: mass
    mass = mass_type (prt%type, prt%mass2)
  end function parton_mass

  elemental module function parton_mass_squared (prt) result (mass_squared)
    class(parton_t), intent(in) :: prt
    real(default) :: mass_squared
    mass_squared = mass_squared_type (prt%type, prt%mass2)
  end function parton_mass_squared

  pure module function parton_momentum_to_pythia6 (prt) result (p)
    real(double), dimension(1:5) :: p
    class(parton_t), intent(in) :: prt
    p = prt%momentum%to_pythia6 (prt%mass ())
  end function parton_momentum_to_pythia6

  module function P_prt_to_child1 (prt) result (retvalue)
    type(parton_t), intent(in) :: prt
    real(default) :: retvalue
    retvalue = zero
    if (prt%is_gluon ()) then
       if (prt%child1%is_quark ()) then
          retvalue = P_gqq (prt%z)
       else if (prt%child1%is_gluon ()) then
          retvalue = P_ggg (prt%z) + P_ggg (one - prt%z)
       end if
    else if (prt%is_quark ()) then
       if (prt%child1%is_quark ()) then
          retvalue = P_qqg (prt%z)
       else if (prt%child1%is_gluon ()) then
          retvalue = P_qqg (one - prt%z)
       end if
    end if
  end function P_prt_to_child1

  module function thetabar (prt, recoiler, isr_ang, E3out) result (retvalue)
    type(parton_t), intent(inout) :: prt
    type(parton_t), intent(in) :: recoiler
    real(default), intent(out), optional :: E3out
    logical, intent(in) :: isr_ang
    logical :: retvalue
    real(default) :: ctheta, cthetachild1
    real(default) p1, p4, p3, E3, shat

    shat = (prt%child1%momentum + recoiler%momentum)**2
    E3 = 0.5_default * (shat / prt%z -recoiler%t + prt%child1%t - &
         prt%child2%mass_squared ()) / sqrt(shat)
    if (present (E3out)) then
       E3out = E3
    end if
    !!! absolute values of momenta in a 3 -> 1 + 4 branching
    p3 = sqrt (E3**2 - prt%t)
    p1 = sqrt (prt%child1%momentum%p(0)**2 - prt%child1%t)
    p4 = sqrt (max (zero, (E3 - prt%child1%momentum%p(0))**2 &
                          - prt%child2%t))
    if (p3 > zero) then
       retvalue = ((p1 + p4 >= p3) .and. (p3 >= abs(p1 - p4)) )
       if (retvalue .and. isr_ang) then
          !!! check angular ordering
          if (associated (prt%child1)) then
             if (associated (prt%child1%child2)) then
                ctheta = (E3**2 - p1**2 - p4**2 + prt%t) / (two * p1 * p4)
                cthetachild1 = (prt%child1%momentum%p(0)**2 - &
                     space_part (prt%child1%child1%momentum)**2 &
                     - space_part (prt%child1%child2%momentum)**2 + prt%child1%t) &
                     / (two * space_part (prt%child1%child1%momentum)**1 * &
                              space_part (prt%child1%child2%momentum)**1)
                retvalue = (ctheta > cthetachild1)
             end if
          end if
       end if
    else
       retvalue = .false.
    end if
  end function thetabar

  recursive module subroutine parton_apply_costheta (prt, rng)
    type(parton_t), intent(inout) :: prt
    class(rng_t), intent(inout), allocatable :: rng
    if (debug2_active (D_SHOWER)) then
       print *, "D: parton_apply_costheta for parton " , prt%nr
       print *, 'prt%momentum%p =    ', prt%momentum%p
       if (debug_on) call msg_debug2 (D_SHOWER, "prt%type", prt%type)
    end if
    prt%z = 0.5_default * (one + prt%get_beta () * prt%costheta)
    if (associated (prt%child1) .and. associated (prt%child2)) then
       if (prt%child1%simulated .and. prt%child2%simulated) then
          prt%z = 0.5_default * (one + (prt%child1%t - prt%child2%t) / &
               prt%t + prt%get_beta () * prt%costheta * &
                sqrt((prt%t - prt%child1%t - prt%child2%t)**2 - &
                4 * prt%child1%t * prt%child2%t) / prt%t)
          if (prt%type /= INTERNAL) then
             prt%child1%momentum%p(0) = prt%z * prt%momentum%p(0)
             prt%child2%momentum%p(0) = (one - prt%z) * prt%momentum%p(0)
          end if
          call prt%generate_ps (rng)
          call parton_apply_costheta (prt%child1, rng)
          call parton_apply_costheta (prt%child2, rng)
       end if
    end if
  end subroutine parton_apply_costheta

  module subroutine parton_apply_lorentztrafo (prt, L)
    type(parton_t), intent(inout) :: prt
    type(lorentz_transformation_t), intent(in) :: L
    prt%momentum = L * prt%momentum
  end subroutine parton_apply_lorentztrafo

  recursive module subroutine parton_apply_lorentztrafo_recursive (prt, L)
    type(parton_t), intent(inout) :: prt
    type(lorentz_transformation_t) ,intent(in) :: L
    if (prt%type /= PROTON .and. prt%type /= BEAM_REMNANT) then
       !!! don't boost hadrons and beam-remnants
       call parton_apply_lorentztrafo (prt, L)
    end if
    if (associated (prt%child1) .and. associated (prt%child2)) then
       if ((space_part_norm (prt%child1%momentum) < eps0) .and. &
           (space_part_norm (prt%child2%momentum) < eps0) .and. &
           (.not. prt%child1%belongstointeraction) .and. &
           (.not. prt%child2%belongstointeraction)) then
          !!! don't boost unevolved timelike partons
       else
          call parton_apply_lorentztrafo_recursive (prt%child1, L)
          call parton_apply_lorentztrafo_recursive (prt%child2, L)
       end if
    else
       if (associated (prt%child1)) then
          call parton_apply_lorentztrafo_recursive (prt%child1, L)
       end if
       if (associated (prt%child2)) then
          call parton_apply_lorentztrafo_recursive (prt%child2, L)
       end if
    end if
  end subroutine parton_apply_lorentztrafo_recursive

  module subroutine parton_generate_ps (prt, rng)
    class(parton_t), intent(inout) :: prt
    class(rng_t), intent(inout), allocatable :: rng
    real(default), dimension(1:3, 1:3) :: directions
    integer i,j
    real(default) :: scproduct, pabs, p1abs, p2abs, x, ptabs, phi
    real(default), dimension(1:3) :: momentum
    type(vector3_t) :: pchild1_direction
    type(lorentz_transformation_t) :: L, rotation
    if (debug2_active (D_SHOWER)) &
         print *, "D: parton_generate_ps for parton " , prt%nr
    if (debug_active (D_SHOWER)) then
       if (.not. (associated (prt%child1) .and. associated (prt%child2))) then
          call msg_fatal ("no children for generate_ps")
       end if
    end if
    !!! test if parton is a virtual parton from the imagined parton shower history
    if (prt%type == INTERNAL) then
       L = inverse (boost (prt%momentum, sqrt(prt%t)))
       !!! boost to restframe of mother
       call parton_apply_lorentztrafo (prt, L)
       call parton_apply_lorentztrafo (prt%child1, L)
       call parton_apply_lorentztrafo (prt%child2, L)
       !!! Store child1's momenta
       pchild1_direction = direction (space_part (prt%child1%momentum))
       !!! Redistribute energy
       prt%child1%momentum%p(0) = (prt%momentum%p(0)**2 - &
            prt%child2%t + prt%child1%t) / (two * prt%momentum%p(0))
       prt%child2%momentum%p(0) = prt%momentum%p(0) - &
            prt%child1%momentum%p(0)

       ! rescale momenta and set momenta to be along z-axis
       prt%child1%momentum = vector4_moving (prt%child1%momentum%p(0), &
            vector3_canonical(3) * &
            sqrt(prt%child1%momentum%p(0)**2 - prt%child1%t))
       prt%child2%momentum = vector4_moving (prt%child2%momentum%p(0), &
            - vector3_canonical(3) * &
            sqrt(prt%child2%momentum%p(0)**2 - prt%child2%t))

       !!! rotate so that total momentum is along former total momentum
       rotation = rotation_to_2nd (space_part (prt%child1%momentum), &
            pchild1_direction)
       call parton_apply_lorentztrafo (prt%child1, rotation)
       call parton_apply_lorentztrafo (prt%child2, rotation)

       L = inverse (L)             !!! inverse of the boost to restframe of mother
       call parton_apply_lorentztrafo (prt, L)
       call parton_apply_lorentztrafo (prt%child1, L)
       call parton_apply_lorentztrafo (prt%child2, L)
    else
       !!! directions(1,:) -> direction of the parent parton
       if (space_part_norm (prt%momentum) < eps0) return
       directions(1,1:3) = prt%momentum%p(1:3) / space_part_norm (prt%momentum)
       !!! directions(2,:) and directions(3,:) -> two random directions
       !!!   perpendicular to the direction of the parent parton
        do j = 2, 3
           call rng%generate (directions(j,:))
        end do
       do i = 2, 3
          scproduct = zero
          do j = 1, i - 1
             scproduct = directions(i,1) * directions(j,1) + &
                  directions(i,2) * directions(j,2) + &
                  directions(i,3) * directions(j,3)
             directions(i,1) = directions(i,1) - directions(j,1) * scproduct
             directions(i,2) = directions(i,2) - directions(j,2) * scproduct
             directions(i,3) = directions(i,3) - directions(j,3) * scproduct
          end do
          scproduct = directions(i,1)**2 + directions(i,2)**2 + &
               directions(i,3)**2
          do j = 1, 3
             directions(i,j) = directions(i,j) / sqrt(scproduct)
          end do
       end do
       if ((directions(1,1) * (directions(2,2) * directions(3,3) - &
            directions(2,3) * directions(3,2)) + &
            directions(1,2) * (directions(2,3) * directions(3,1) - &
            directions(2,1) * directions(3,3)) + &
            directions(1,3) * (directions(2,1) * directions(3,2) - &
            directions(2,2) * directions(3,1))) < 0) then
          directions(3,:) = - directions(3,:)
       end if

       pabs = space_part_norm (prt%momentum)
       if ((prt%child1%momentum%p(0)**2 - prt%child1%t < 0) .or. &
           (prt%child2%momentum%p(0)**2 - prt%child2%t < 0)) then
          if (debug_on) call msg_debug(D_SHOWER, "generate_ps error at E^2 < t")
          return
       end if
       p1abs = sqrt (prt%child1%momentum%p(0)**2 - prt%child1%t)
       p2abs = sqrt (prt%child2%momentum%p(0)**2 - prt%child2%t)
       x = (pabs**2 + p1abs**2 - p2abs**2) / (two * pabs)
       if (pabs > p1abs + p2abs .or. &
            pabs < abs(p1abs - p2abs)) then
          if (debug_active (D_SHOWER)) then
             print *, "D: parton_generate_ps Dreiecksungleichung error &
                  &for parton ", prt%nr, " ", &
                  space_part_norm (prt%momentum), " ", p1abs, " ", p2abs
             call prt%write ()
             call prt%child1%write ()
             call prt%child2%write ()
          end if
          return
       end if
       !!! Due to numerical problems transverse momentum could be imaginary ->
       !!!     set transverse momentum to zero
       ptabs = sqrt (max (p1abs * p1abs - x * x, zero))
       call rng%generate (phi)
       phi = twopi * phi
       do i = 1, 3
          momentum(i) = x * directions(1,i) + ptabs * &
                (cos(phi) * directions(2,i) + sin(phi) * directions(3,i))
       end do
       prt%child1%momentum%p(1:3) = momentum(1:3)
       do i = 1, 3
          momentum(i) = (space_part_norm (prt%momentum) - x) * directions(1,i) - &
               ptabs * (cos(phi) * directions(2,i) + sin(phi) * directions(3,i))
       end do
       prt%child2%momentum%p(1:3) = momentum(1:3)
    end if
  end subroutine parton_generate_ps

  module subroutine parton_generate_ps_ini (prt, rng)
    class(parton_t), intent(inout) :: prt
    class(rng_t), intent(inout), allocatable :: rng
    real(default), dimension(1:3, 1:3) :: directions
    integer :: i,j
    real(default) :: scproduct, pabs, p1abs, p2abs, x, ptabs, phi
    real(default), dimension(1:3) :: momentum
    if (debug_active (D_SHOWER)) &
         print *, "D: parton_generate_ps_ini: for parton " , prt%nr
    if (debug_active (D_SHOWER)) then
       if (.not. (associated (prt%child1) .and. associated (prt%child2))) then
          call msg_fatal ("no children for generate_ps")
       end if
    end if

    if (.not. prt%is_proton()) then
       !!! generate ps for normal partons
       do i = 1, 3
          directions(1,i) = prt%child1%momentum%p(i) / &
               space_part_norm(prt%child1%momentum)
       end do
       do j = 2, 3
          call rng%generate (directions(j,:))
       end do
       do i = 2, 3
          scproduct = zero
          do j = 1, i - 1
             scproduct = directions(i,1) * directions(j,1) + &
                  directions(i,2) * directions(j,2) + &
                  directions(i,3) * directions(j,3)
             directions(i,1) = directions(i,1) - directions(j,1) * scproduct
             directions(i,2) = directions(i,2) - directions(j,2) * scproduct
             directions(i,3) = directions(i,3) - directions(j,3) * scproduct
          end do
          scproduct = directions(i,1)**2 + directions(i,2)**2 + &
               directions(i,3)**2
          do j = 1, 3
             directions(i,j) = directions(i,j) / sqrt(scproduct)
          end do
       end do
       if ((directions(1,1) * (directions(2,2) * directions(3,3) - &
            directions(2,3) * directions(3,2)) + &
            directions(1,2) * (directions(2,3) * directions(3,1) - &
            directions(2,1) * directions(3,3)) + &
            directions(1,3) * (directions(2,1) * directions(3,2) - &
            directions(2,2) * directions(3,1))) < 0) then
          directions(3,:) = - directions(3,:)
       end if

       pabs = space_part_norm (prt%child1%momentum)
       p1abs = sqrt (prt%momentum%p(0)**2 - prt%t)
       p2abs = sqrt (max(zero, prt%child2%momentum%p(0)**2 - &
            prt%child2%t))

       x = (pabs**2 + p1abs**2 - p2abs**2) / (two * pabs)
       if (debug_active (D_SHOWER)) then
          if (pabs > p1abs + p2abs .or. pabs < abs(p1abs - p2abs)) then
             print *, "error at generate_ps, Dreiecksungleichung for parton ", &
                  prt%nr, " ", pabs," ",p1abs," ",p2abs
             call prt%write ()
             call prt%child1%write ()
             call prt%child2%write ()
             call msg_fatal ("parton_generate_ps_ini: Dreiecksungleichung")
          end if
       end if
       if (debug_active (D_SHOWER)) print *, "D: parton_generate_ps_ini: x = ", x
       ptabs = sqrt (p1abs * p1abs - x**2)
       call rng%generate (phi)
       phi = twopi * phi
       do i = 1,3
          momentum(i) = x * directions(1,i) + ptabs * (cos(phi) * &
               directions(2,i) + sin(phi) * directions(3,i))
       end do
       prt%momentum%p(1:3) = momentum
       do i = 1, 3
          momentum(i) = (x - pabs) * directions(1,i) + ptabs * (cos(phi) * &
               directions(2,i) + sin(phi) * directions(3,i))
       end do
       prt%child2%momentum%p(1:3) = momentum(1:3)
    else
       !!! for first partons just set beam remnants momentum
       prt%child2%momentum = prt%momentum - prt%child1%momentum
    end if
  end subroutine parton_generate_ps_ini

  module subroutine parton_next_t_ana (prt, rng)
    class(parton_t), intent(inout) :: prt
    class(rng_t), intent(inout), allocatable :: rng
    integer :: gtoqq
    real(default) :: integral, random
    if (signal_is_pending ()) return
    if (debug_on) call msg_debug (D_SHOWER, "next_t_ana")
    ! check if branchings are possible at all
    if (min (prt%t, prt%momentum%p(0)**2) < &
         prt%mass_squared () + prt%settings%min_virtuality) then
       prt%t = prt%mass_squared ()
       call prt%set_simulated ()
       return
    end if
    integral = zero
    call rng%generate (random)
    do
       call parton_simulate_stept (prt, rng, integral, random, gtoqq, .false.)
       if (prt%simulated) then
          if (prt%is_gluon ()) then
             !!! Abusing the x-variable to store the information to which
             !!! quark flavor the gluon branches (if any)
             prt%x = one * gtoqq + 0.1_default
             !!! x = gtoqq + 0.1 -> int(x) will be the quark flavor or
             !!! zero for g -> gg
          end if
          exit
       end if
    end do
  end subroutine parton_next_t_ana

  function cmax (prt, tt) result (cmaxx)
    type(parton_t), intent(in) :: prt
    real(default), intent(in), optional :: tt
    real(default) :: t, cost, cmaxx, radicand
    t = prt%t;  if (present (tt))  t = tt
    if (associated (prt%parent)) then
       cost = prt%parent%get_costheta ()
       radicand = max(zero, one - &
            t / (prt%get_beta () * prt%momentum%p(0))**2 * &
            (one + cost) / (one - cost))
       if (debug_on) call msg_debug2 (D_SHOWER, "cmax: sqrt (radicand)", sqrt (radicand))
       cmaxx = min (0.99999_default, sqrt (radicand))
    else
       cmaxx = 0.99999_default
    end if
  end function cmax

  module subroutine parton_simulate_stept &
       (prt, rng, integral, random, gtoqq, lookatsister)
    type(parton_t), intent(inout) :: prt
    class(rng_t), intent(inout), allocatable :: rng
    real(default), intent(inout) :: integral
    real(default), intent(inout) :: random
    integer, intent(out) :: gtoqq
    logical, intent(in), optional :: lookatsister

    type(parton_t), pointer :: sister
    real(default) :: tstep, tmin, oldt
    real(default) :: c, cstep
    real(default), dimension(3) :: z, P
    real(default) :: to_integral
    real(default) :: a11,a12,a13,a21,a22,a23
    real(default) :: cmax_t
    real(default) :: temprand
    real(default), dimension(3) :: a, x

    ! higher values -> faster but coarser
    real(default), parameter :: tstepfactor = 0.02_default
    real(default), parameter :: tstepmin = 0.5_default
    real(default), parameter :: cstepfactor = 0.8_default
    real(default), parameter :: cstepmin = 0.03_default

    if (signal_is_pending ()) return
    if (debug_on) call msg_debug (D_SHOWER, "parton_simulate_stept")
    gtoqq = 111 ! illegal value
    call prt%set_simulated (.false.)

    sister => null()
    SET_SISTER: do
       if (present (lookatsister)) then
          if (.not. lookatsister) then
             exit SET_SISTER
          end if
       end if
       if (prt%nr == prt%parent%child1%nr) then
          sister => prt%parent%child2
       else
          sister => prt%parent%child1
       end if
       exit SET_SISTER
    end do SET_SISTER

    tmin = prt%settings%min_virtuality + prt%mass_squared ()
    if (prt%is_quark ()) then
       to_integral = three *pi * log(one / random)
    else if (prt%is_gluon ()) then
       to_integral = four *pi * log(one / random)
    else
       prt%t = prt%mass_squared ()
       call prt%set_simulated ()
       return
    end if

    if (associated (sister)) then
       if (sqrt(prt%t) > sqrt(prt%parent%t) - &
            sqrt(sister%mass_squared ())) then
          prt%t = (sqrt (prt%parent%t) - sqrt (sister%mass_squared ()))**2
       end if
    end if
    if (prt%t > prt%momentum%p(0)**2) then
        prt%t = prt%momentum%p(0)**2
    end if

    if (prt%t <= tmin) then
       prt%t = prt%mass_squared ()
       call prt%set_simulated ()
       return
    end if

    ! simulate the branchings between prt%t and prt%t - tstep
    tstep = max(tstepfactor * (prt%t - 0.9_default * tmin), tstepmin)
    cmax_t = cmax(prt)
    c = - cmax_t ! take highest t -> minimal constraint
    cstep = max(cstepfactor * (one - abs(c)), cstepmin)
    ! get values at border of "previous" bin -> to be used in first bin
    z(3) = 0.5_default + 0.5_default * get_beta (prt%t - &
         0.5_default * tstep, prt%momentum%p(0)) * c
    if (prt%is_gluon ()) then
       P(3) = P_ggg (z(3)) + P_gqq (z(3)) * number_of_flavors &
            (prt%t, prt%settings%max_n_flavors, prt%settings%min_virtuality)
    else
       P(3) = P_qqg (z(3))
    end if
    a(3) = D_alpha_s_fsr (z(3) * (one - z(3)) * prt%t, &
         prt%settings) * P(3) / (prt%t - 0.5_default * tstep)

    do while (c < cmax_t .and. (integral < to_integral))
       if (signal_is_pending ()) return
       cmax_t = cmax (prt)
       cstep = max (cstepfactor * (one - abs(c)**2), cstepmin)
       if (c + cstep > cmax_t) then
          cstep = cmax_t - c
       end if
       if (cstep < 1E-9_default) then
          !!! reject too small bins
          exit
       end if
       z(1) = z(3)
       z(2) = 0.5_default + 0.5_default * get_beta &
            (prt%t - 0.5_default * tstep, prt%momentum%p(0)) * &
            (c + 0.5_default * cstep)
       z(3) = 0.5_default + 0.5_default * get_beta &
            (prt%t - 0.5_default * tstep, prt%momentum%p(0)) * (c + cstep)
       P(1) = P(3)
       if (prt%is_gluon ()) then
          P(2) = P_ggg(z(2)) + P_gqq(z(2)) * number_of_flavors &
               (prt%t, prt%settings%max_n_flavors, prt%settings%min_virtuality)
          P(3) = P_ggg(z(3)) + P_gqq(z(3)) * number_of_flavors &
               (prt%t, prt%settings%max_n_flavors, prt%settings%min_virtuality)
       else
          P(2) = P_qqg(z(2))
          P(3) = P_qqg(z(3))
       end if
       ! get values at borders of the intgral and in the middle
       a(1) = a(3)
       a(2) = D_alpha_s_fsr (z(2) * (one - z(2)) * prt%t, &
            prt%settings) * P(2) / &
            (prt%t - 0.5_default * tstep)
       a(3) = D_alpha_s_fsr (z(3) * (one - z(3)) * prt%t, &
            prt%settings) * P(3) / &
            (prt%t - 0.5_default * tstep)

       !!! a little tricky:
       !!! fit x(1) + x(2)/(1 + c) + x(3)/(1 - c) to these values
       a11 = (one+c+0.5_default*cstep) * (one-c-0.5_default*cstep) - &
             (one-c) * (one+c+0.5_default*cstep)
       a12 = (one-c-0.5_default*cstep) - (one+c+0.5_default*cstep) * &
             (one-c) / (one+c)
       a13 = a(2) * (one+c+0.5_default*cstep) * (one-c-0.5_default*cstep) - &
             a(1) * (one-c) * (one+c+0.5_default*cstep)
       a21 = (one+c+cstep) * (one-c-cstep) - (one+c+cstep) * (one-c)
       a22 = (one-c-cstep) - (one+c+cstep) * (one-c) / (one+c)
       a23 = a(3) * (one+c+cstep) * (one-c-cstep) - &
            a(1) * (one-c) * (one+c+cstep)

       x(2) = (a23 - a21 * a13 / a11) / (a22 - a12 * a21 / a11)
       x(1) = (a13 - a12 * x(2)) / a11
       x(3) = a(1) * (one - c) - x(1) * (one - c) - x(2) * (one - c) / (one + c)

       integral = integral + tstep * (x(1) * cstep + x(2) * &
            log((one + c + cstep) / (one + c)) - x(3) * &
            log((one - c - cstep) / (one - c)))

       if (integral > to_integral) then
          oldt = prt%t
          call rng%generate (temprand)
          prt%t = prt%t - temprand * tstep
          call rng%generate (temprand)
          prt%costheta = c + (0.5_default - temprand) * cstep
          call prt%set_simulated ()

          if (prt%t < prt%settings%min_virtuality + prt%mass_squared ()) then
             prt%t = prt%mass_squared ()
          end if
          if (abs(prt%costheta) > cmax_t) then
             ! reject branching due to violation of costheta-limits
             call rng%generate (random)
             if (prt%is_quark ()) then
                to_integral = three * pi * log(one / random)
             else if (prt%is_gluon ()) then
                to_integral = four * pi * log(one / random)
             end if
             integral = zero
             prt%t = oldt
             call prt%set_simulated (.false.)
          end if
          if (prt%is_gluon ()) then
             ! decide between g->gg and g->qqbar splitting
             z(1) = 0.5_default + 0.5_default * prt%costheta
             call rng%generate (temprand)
             if (P_ggg(z(1)) > temprand * (P_ggg (z(1)) + P_gqq (z(1)) * &
                  number_of_flavors(prt%t, prt%settings%max_n_flavors, &
                  prt%settings%min_virtuality))) then
                gtoqq = 0
             else
                call rng%generate (temprand)
                gtoqq = 1 + int (temprand * number_of_flavors &
                     (prt%t, prt%settings%max_n_flavors, &
                      prt%settings%min_virtuality))
             end if
          end if
       else
          c = c + cstep
       end if
       cmax_t = cmax (prt)
    end do
    if (integral <= to_integral) then
       prt%t = prt%t - tstep
       if (prt%t < prt%settings%min_virtuality + prt%mass_squared ()) then
          prt%t = prt%mass_squared ()
          call prt%set_simulated ()
       end if
    end if
  end subroutine parton_simulate_stept

  module function maxzz (shat, s, maxz_isr, minenergy_timelike) result (maxz)
    real(default), intent(in) :: shat, s, minenergy_timelike, maxz_isr
    real(default) :: maxz
    maxz = min (maxz_isr, one - (two * minenergy_timelike * sqrt(shat)) / s)
  end function maxzz


end submodule shower_partons_s

