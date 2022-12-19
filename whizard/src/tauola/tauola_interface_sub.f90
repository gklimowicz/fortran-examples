!
!  WHIZARD Tauola interface
!  Adapted from ilc_tauola_mod.f90
!  for Whizard1 developed by Timothy Barklow (SLAC)
!
!    Akiya Miyamoto
!    Bug fixes: Mikael Berggren, Juergen Reuter
!

submodule (tauola_interface) tauola_interface_s

  use io_units
  use constants
  use format_utils, only: write_separator
  use debug_master, only: debug_on
  use diagnostics
  use hep_common
  use hepev4_aux

  implicit none

  integer :: jtau2, jorig, jforig
  integer :: nproducts

  integer, parameter :: n_pyjets_max = 4000

  integer, save :: nsub_call = 0

  double precision :: spin_dexay

  double precision, dimension(4) :: pol
  logical :: higgs_dec

  !!! Probability of tau- to be left-handed in Z decays
  double precision, parameter :: a_tau = 0.15
  double precision, parameter :: prob_tau_left_z = (a_tau+1.) / 2.

  type(pyjets_spin_t), dimension(n_pyjets_max), save :: pyjets_spin_data
  integer :: last_event_number = -999

contains

  subroutine fill_pyjets_spin_data
    integer :: ip
    integer :: hepeup_index
    integer :: iorig, iii, ihep_start, ipy_start, first_pdg
    integer :: idau1, idau2, n_doc_lines
    integer, dimension(200) :: mstp
    double precision, dimension(200) :: parp
    integer, dimension(200) :: msti
    double precision, dimension(200) :: pari
    common/pypars/mstp,parp,msti,pari
    save/pypars/
    integer :: n
    integer :: npad
    integer, dimension(4000,5) :: k
    double precision, dimension(4000,5) :: p
    double precision, dimension(4000,5) :: v
    common/pyjets/n,npad,k,p,v
    save/pyjets/
    !!! Set helicity information of document lines at the first call of
    !!! this event
    !!! MSTI(4) ; number of documentation lines
    !!! MSTI(5) ; number of events generated
    if (last_event_number == MSTI(5)) then
       return
    end if
    do ip = 1, n_pyjets_max
       pyjets_spin_data(ip)%index_to_hepeup = 0
       pyjets_spin_data(ip)%helicity = 0
    end do

    ip = 1
    hepeup_index = 0

    if (debug_active (D_TAUOLA)) then
       call msg_debug (D_TAUOLA, "TAUOLA interface: fill_pyjets_spin_data")
       print *,"nup in hepev4 = ", nup
       print *, "        index       spin" // &
            "                      pdg              mothers "
       do iii=1,nup
          print *, iii, spinup(iii), idup(iii),mothup(:,iii)
       end do
       if (debug2_active (D_TAUOLA)) then
          call pylist(2)
       end if
    end if

    do ip = 1, Nup
      if (mothup(1,ip) /= 0) exit
    end do
    ihep_start = ip
    first_pdg = idup(ip)
    do ip = 1, N
      if (k(ip,1) == 21) then
        if (k(ip,2) == first_pdg) exit
      end if
    end do
    ipy_start = ip
    if (debug_active (D_TAUOLA)) then
      print *, ' First hepev4 index to consider: ', ihep_start, &
               ' , which correspinds to  pyjets index : ', ipy_start
    end if


    hepeup_index = ihep_start
    n_doc_lines = count (k(1:N,1) == 21)


    if (signal_is_pending ()) return
    do ip=1,n_doc_lines
       if (k(ip,1) /= 21) exit
       pyjets_spin_data(ip)%pid     = k(ip,2)
       pyjets_spin_data(ip)%id_orig = k(ip,3)
       pyjets_spin_data(ip)%index_orig = ip
       pyjets_spin_data(ip)%n_daughter = 0
       iorig = k(ip,3)
       if ((ip >= ipy_start) .and. (.not. any(spinup(hepeup_index) == [0.0, 9.0]))) then
          pyjets_spin_data(ip)%index_to_hepeup = hepeup_index
          pyjets_spin_data(ip)%helicity = spinup(hepeup_index)
       else
          pyjets_spin_data(ip)%index_to_hepeup = -1
          if (ip >= ipy_start) then
            pyjets_spin_data(ip)%helicity = spinup(hepeup_index)
          else
            pyjets_spin_data(ip)%helicity = 0.
          end if
          if (iorig > 0) then
            pyjets_spin_data(iorig)%n_daughter = &
            pyjets_spin_data(iorig)%n_daughter + 1
            pyjets_spin_data(iorig)%index_daughter(pyjets_spin_data(iorig)%n_daughter)=ip
          end if
          if (debug_active (D_TAUOLA)) then
             if (abs(k(ip,2)) == 15) then
                print *, " No input spin information for tau at ip = ", ip, &
                     ". Will be set in fill_pyjets_spin_data"
             end if
          end if
       end if
       if (ip >= ipy_start) hepeup_index = hepeup_index + 1
       if (debug2_active (D_TAUOLA)) then
          call msg_debug2 (D_TAUOLA, "TAUOLA interface: fill_pyjets_spin_data")
          write (msg_buffer, "(A,I0,A,I0,A,I0,A,ES19.12,A,I0)") &
               "ip = ", ip, " iorig = ", iorig, " pid  ", k(ip,2), &
               " spin = ", pyjets_spin_data(ip)%helicity,&
               " hepeup index = ", pyjets_spin_data(ip)%index_to_hepeup
          call msg_message ()
       end if
       if (abs(k(ip,2)) == 15 .and. pyjets_spin_data(ip)%helicity == 9) then
          if (nsub_call .lt. 5) then
             write (msg_buffer, "(A)") &
                  "Subroutine fill_pyjets_spin_data: tau helicity information"
             call msg_message ()
             write (msg_buffer, "(A)") &
                  "is not set, though polarized tau decay was requested."
             call msg_message ()
             write (msg_buffer, "(A)") &
                  "Most likely, the SINDARIN file does not include polarized"
             call msg_message ()
             write (msg_buffer, "(A)") &
                  "for particles and/or not ?polarized_events=true"
             call msg_message ()
             write (msg_buffer, "(A,I0,A,I0,A,I0,A,I0,A,ES19.12)") &
                  "Number of calls:", nsub_call, " ip = ", ip, " iorig = ", &
                  iorig, " pid = ", k(ip,2), " spin = ", &
                  pyjets_spin_data(ip)%helicity
          end if
       end if
    end do

    if (signal_is_pending ()) return
    do ip = 1, n_doc_lines
       iorig = k(ip,3)
       if (pyjets_spin_data(ip)%n_daughter == 2) then
          !!! h0/H0/A0 -> tau tau
          if (pyjets_spin_data(ip)%pid == 25 .or. &
               pyjets_spin_data(ip)%pid == 35 .or. &
               pyjets_spin_data(ip)%pid == 36) then
             idau1 = pyjets_spin_data(ip)%index_daughter(1)
             idau2 = pyjets_spin_data(ip)%index_daughter(2)
             if (abs(pyjets_spin_data(idau1)%pid) == 15 .and. &
                  (pyjets_spin_data(idau1)%pid + &
                  pyjets_spin_data(idau2)%pid) == 0) then
                if (pyr(0) .lt. 0.5) then
                   pyjets_spin_data(idau1)%helicity = -1
                   pyjets_spin_data(idau2)%helicity = -1
                else
                   pyjets_spin_data(idau1)%helicity = +1
                   pyjets_spin_data(idau2)%helicity = +1
                end if
             end if

             !!! Z0 -> tau tau
          else if (pyjets_spin_data(ip)%pid == 23) then
             idau1 = pyjets_spin_data(ip)%index_daughter(1)
             idau2 = pyjets_spin_data(ip)%index_daughter(2)
             if (abs(pyjets_spin_data(idau1)%pid) == 15 .and. &
                  (pyjets_spin_data(idau1)%pid + &
                  pyjets_spin_data(idau2)%pid) == 0) then
                if (((pyr(0) - prob_tau_left_z) * &
                     pyjets_spin_data(idau1)%pid) .gt. 0.0) then
                   pyjets_spin_data(idau1)%helicity = +1
                   pyjets_spin_data(idau2)%helicity = -1
                else
                   pyjets_spin_data(idau1)%helicity = -1
                   pyjets_spin_data(idau2)%helicity = +1
                end if
             end if

             !!! W+(24)/H+(37) -> tau+(-15) and neu_tau
          else if (pyjets_spin_data(ip)%pid == 24 .or. &
               pyjets_spin_data(ip)%pid == 37) then
             idau1 = pyjets_spin_data(ip)%index_daughter(1)
             idau2 = pyjets_spin_data(ip)%index_daughter(2)
             if (pyjets_spin_data(idau1)%pid == -15) then
                pyjets_spin_data(idau1)%helicity = 1
             else if (pyjets_spin_data(idau2)%pid == -15) then
                pyjets_spin_data(idau2)%helicity = 1
             end if

             !!! W-(-24)/H-(-37) -> tau-(+15) and neu_tau_bar
          else if (pyjets_spin_data(ip)%pid == -24 .or. &
               pyjets_spin_data(ip)%pid == -37) then
             idau1 = pyjets_spin_data(ip)%index_daughter(1)
             idau2 = pyjets_spin_data(ip)%index_daughter(2)
             if (pyjets_spin_data(idau1)%pid == 15) then
                pyjets_spin_data(idau1)%helicity = -1
             else if (pyjets_spin_data(idau2)%pid == 15) then
                pyjets_spin_data(idau2)%helicity = -1
             end if
          end if
       end if
       if (debug2_active (D_TAUOLA)) then
          call msg_debug2 (D_TAUOLA, "TAUOLA interface: end of fill_pyjets_spin_data")
          write (msg_buffer, "(A,I0,A,I0,A,I0,A,ES19.12,A,I0)") &
               "ip = ", ip, " iorig = ", iorig, " pid  ", k(ip,2), &
               " spin = ", pyjets_spin_data(ip)%helicity,&
               " hepeup index = ", pyjets_spin_data(ip)%index_to_hepeup
          call msg_message ()
       end if
    end do
  end subroutine fill_pyjets_spin_data

! =====================================================================
!  Main interface to tauola.
!  Called by PYTAUD and calls TAUOLA
! =====================================================================
  module subroutine wo_tauola_pytaud (itau, iorig, kforig, ndecay)
    !!! Line number in /JETSET/ where the tau is stored
    integer, intent(in)  :: itau
    !!! Line number where the mother is stored. =0 if the mother is not stored
    integer, intent(in)  :: iorig
    !!! Flavour code of the mother. 0 unknown. H0(25), W+-(+-24),
    !!!   gamma*/Z=23, H+-(+-37)
    integer, intent(in)  :: kforig
    !!! Number of decay products to be given by user routine.
    integer, intent(out) :: ndecay
    double precision, dimension(5) :: p_dexay
    integer :: id_dexay
    integer :: ip
    integer :: n
    integer :: npad
    integer, dimension(4000,5) :: k
    double precision, dimension(4000,5) :: p
    double precision, dimension(4000,5) :: v
    common /pyjets/ n, npad, k, p, v
    save /pyjets/
    integer, dimension(200) :: mstp
    double precision, dimension(200) :: parp
    integer, dimension(200) :: msti
    double precision, dimension(200) :: pari
    common /pypars/ mstp, parp, msti, pari
    save /pypars/
    !!! TODO: (bcn 2016-03-11) this should only be
    !!!       called once per event (not per tau)
    integer :: idau1, idau2

    integer, allocatable :: these_taus(:), these_tau_daughters(:)
    logical, allocatable :: chain(:)
    integer :: i,ntaus, itp , offset

    higgs_dec = .false.
    if (kforig == 25 .or. kforig==35 .or. kforig== 36) then
      higgs_dec = .true.
    end if

    if (signal_is_pending ()) return
    !!! JRR: Tau decays are very sensitive to numerical noise, momenta
    !!!      should be, in principle, strictly on the z axis
    if (abs (p(itau,1)) < 1.d-13)  p(itau,1) = 0
    if (abs (p(itau,2)) < 1.d-13)  p(itau,2) = 0

    !!! MSTI(4): number of documentation lines
    !!! MSTI(5): number of events generated
    if (last_event_number .ne. MSTI(5)) then
       call fill_pyjets_spin_data
       last_event_number = MSTI(5)
       jtau2 = -1000
       nsub_call = nsub_call + 1
    end if

    if (debug_active (D_TAUOLA)) then
       write (msg_buffer, "(A)") "wo_tauola_pytaud was called."
       call msg_message ()
       write (msg_buffer, "(A,I0,A,I0,A,I0,A,I0)") &
            "ncall = ", nsub_call, " itau = ", itau, " iorig = ", iorig, &
            " kforig = ", kforig
       call msg_message ()
       if (debug2_active (D_TAUOLA)) then
          call pylist(2)
       end if
    end if

    jorig = iorig
    jforig = kforig

    offset=0 ; if (MSTP(125) > 1) offset = mstp(126)-count(k(1:mstp(126)-1,1)==21)
    these_taus=PACK([(i,i=1,n)],(k(1:n,2)==k(itau,2).and.k(1:n,1)/=21)) ; ntaus=size(these_taus) 
    these_tau_daughters=modulo(k(these_taus,4),10000)  ! 10000 = mstu(5) in /pydat1/

    if (debug2_active (D_TAUOLA)) then
       print *, ' offsets etc. : ',offset,  MSTP(125),  MSTP(126)
       print *, ' ntaus = ', ntaus
       print *, ' these_taus          = ', these_taus
       print *, ' these_tau_daughters = ', these_tau_daughters
       print *, ' these_taus parent on line(s)  = ', k(these_taus(1:ntaus),3)
       print *, ' helicty of parents : ', pyjets_spin_data(k(these_taus(1:ntaus),3))%helicity
    end if

    if (.not. all(these_tau_daughters(1:ntaus-1) == these_taus(2:ntaus))) then
      ! maybe multi-tau, sqeeze out the wrong chain(s)
       allocate (chain(ntaus))
       chain=.FALSE.
       itp = transfer (maxloc(these_taus,(these_taus==itau)),itp)

       if (debug2_active (D_TAUOLA)) then
          print *, ' tau and daughter lists NOT equal. Checking for multiple taus ... '
          print *, ' initial itp = ', itp
       end if

       chain(itp) = .true.
       do while (itp > 0)
          if (signal_is_pending ()) return
          if (size (these_taus) /= 0 .and. any (these_tau_daughters == these_taus(itp))) then
             itp=transfer(maxloc(these_taus,(these_tau_daughters==these_taus(itp))),itp)
          else
             itp = 0
          end if
          if (itp > 0) chain(itp) = .true.
       end do

       ntaus = count(chain)
       these_taus(1:ntaus) = pack (these_taus, chain)
       these_tau_daughters(1:ntaus) = pack (these_tau_daughters, chain)

       if (debug2_active (D_TAUOLA)) then
          print *, ' after multi-tau loop : '
          print *, ' ntaus = ', ntaus
          print *, ' chain = ', chain
          print *, ' these_taus          = ', these_taus
          print *, ' these_tau_daughters = ', these_tau_daughters
       end if

    end if
    if (all(these_tau_daughters(1:ntaus-1) == these_taus(2:ntaus)) .and. ntaus > 0) then
       pyjets_spin_data(these_taus(1:ntaus))%helicity = pyjets_spin_data(k(these_taus(1),3))%helicity
       if (offset /= 0) pyjets_spin_data(these_taus(1:ntaus)-offset)%helicity = &
            pyjets_spin_data(k(these_taus(1),3))%helicity
    else
       call msg_warning ("Tau lepton not found: ntaus is " // int2char (ntaus) // ".")
    end if

    if (signal_is_pending ()) return

    if (iorig /= 0 .and. higgs_dec .and. trans_spin) then

      idau1 = pyjets_spin_data(iorig)%index_daughter(1)
      idau2 = pyjets_spin_data(iorig)%index_daughter(2)
             !!! Parent Higgs is not in the documentation line (K(,1) != 21)
             !!! Get pointer to daughter directly from JETSET
      if (idau1 == 0) then
        idau1 = k(iorig,4)
        idau2 = k(iorig,5)
      end if
      if (idau1 .ne. itau) then
        write (msg_buffer, "(A,I0,A,I0,A)") &
                     "idau1 = ", idau1, "itau = ", itau, " are not equal."
        call msg_fatal ("wo_tauola_pytaud: " // &
                     "Something is wrong in parent-daughter relation.")
      end if
      jtau2 = idau2
             !!! Reset tau spin information because it is decided internally
      pyjets_spin_data(itau)%helicity = 0
      pyjets_spin_data(jtau2)%helicity = 0
    end if

    id_dexay = k(itau,2)
    p_dexay = p(itau,1:5)
    spin_dexay = pyjets_spin_data(itau)%helicity

    call do_dexay (itau, p_dexay, id_dexay, kforig)
    ndecay = nproducts
  end subroutine wo_tauola_pytaud

  subroutine do_dexay (itau, p_dexay, id_dexay, kforig)
    !!! Main routine to call Tauola. Three type of tau decay:
    !!! (A) Higgs to tau+tau-
    !!! (B) single tau+ decay
    !!! (C) single tau- decay, are treated separately
    integer, intent(in) :: itau, id_dexay
    double precision, dimension(5), intent(in) :: p_dexay
    integer, intent(in) :: kforig

    integer :: i, IFPHOT
    logical :: ifpseudo, is_swapped
    double precision, dimension(4) :: pol1, pol2
    integer :: im
    double precision :: rrr(1), wt
    double precision :: hh1(4), hh2(4)
    integer :: ion(3), np
    common /PHOACT/ IFPHOT
    integer :: n
    integer :: npad
    integer, dimension(4000,5) :: k
    double precision, dimension(4000,5) :: p
    double precision, dimension(4000,5) :: v
    common /pyjets/ n, npad, k, p, v
    save /pyjets/

    integer :: n1, n2

    is_swapped = .false.

    !!! For transverse spin of the Higgs, Higgs and the two taus
    !!!    have to be considered
    if (trans_spin .and. IFPHOT == 1) then
       n1 = 2
       n2 = 3
       nhep = 3
    else
       n1 = 1
       n2 = 2
       nhep = 2
    end if

    tauspin_pyjets(itau) = spin_dexay

    !!! Does SPINHIGGS in tauface_jetset.f
    ifpseudo = kforig == 36

    isthep(n1)   = 1
    idhep(n1)    = id_dexay
    jmohep(:,n1) = 0
    jdahep(:,n1) = 0
    phep(:,n1)   = p_dexay

    isthep(n2)   = 1
    idhep(n2)    = - id_dexay
    jmohep(:,n2) = 0
    jdahep(:,n2) = 0
    phep(1:3,n2) = - phep(1:3,n1)
    phep(4:5,n2) = phep(4:5,n1)

!!! NOTE (Akiya Miyamoto, 25-March-2016)
!!!  Higgs (h0/H0/A0) to tau+tau- decay is handled here
!!!  in order to implement a transverse spin correlation.
!!!  For this algorithm to work, photon emission from tau
!!!  before decay should be turned off. Since photon emission
!!!  from tau is handled by PYTHIA, photon emission from ALL tau
!!!  should be turned off.  It is done by setting MSTJ(39)=15.
!!!  Instead, PHOTOS is called after tau decay and generate
!!!  photons.

!!! ****************************************************************
!!! (A) Higgs to tau+ tau- decay .
!!! ****************************************************************
    if (higgs_dec .and. trans_spin) then
       if (idhep(2) .gt. 0) then
          idhep(3) = id_dexay
          idhep(2) = - id_dexay
          is_swapped = .true.
       end if

       phep(:,n1) = p_dexay
       phep(1:3,n2) = - phep(1:3,2)
       phep(4:5,n2) = phep(4:5,2)

       isthep(1) = 11
       idhep(1) = kforig
       jmohep(:,1) = 0
       jdahep(1,1) = n1
       jdahep(2,1) = n2
       phep(:,1) = phep(:,2) + phep(:,3)
       phep(5,1) = sqrt(phep(4,1)**2 - phep(1,1)**2 - phep(2,1)**2 - &
            phep(3,1)**2)
       jmohep(:,n1) = 1
       jmohep(:,n2) = 1

       p1=phep(1:4,np1)  ! tau+ momentum
       p2=phep(1:4,np2)  ! tau- momentum

       q1 = p1 + p2
       im = 1
    end if

!!!  tau+ momentum should have positive Pz
!!!  tau- momentum should have negative Pz

!!! ********************************************************
!!! (B) Single Tau+ decay
!!! ********************************************************
    if (.not. trans_spin) then
       check_tau_sign: if (idhep(n1) .lt. 0) then
          np1 = n1
          np2 = n2
          pol = 0.
          pol(3) = - spin_dexay
          p1 = phep(1:4,n1)
          p2 = phep(1:4,n2)
          q1 = p1 + p2
          if (debug_active (D_TAUOLA)) then
             call msg_message ("Tau+ decay with pol(3) = " // &
                  real2char (real (pol(3), kind=default)) // ".")
             write (*, "(A,4(1x,ES19.12))") "Antiparticle decay, q1 = ", q1
             write (*, "(A,4(1x,ES19.12))") "Antiparticle decay, p1 = ", p1
             write (*, "(A,4(1x,ES19.12))") "Antiparticle decay, p2 = ", p2
          end if
          if (debug_on) call msg_debug2 (D_TAUOLA, "TAUOLA is called here")
          call dexay (1,pol)
          if (IFPHOT == 1)  call photos (np1)
!!! ********************************************************
!!! (C) Single Tau- decay
!!! ********************************************************
       else check_tau_sign
          idhep(3) = id_dexay
          idhep(2) = - id_dexay
          np2 = n1
          np1 = n2
          pol = 0.
          pol(3) = spin_dexay
         !!! Akiya now has a relation with the negative spin_dexay
          ! pol(3) = - spin_dexay
          p2 = phep(1:4,n1)
          p1 = phep(1:4,n2)
          q1 = p1 + p2
          if (debug_active (D_TAUOLA)) then
             call msg_message ("Tau- decay with pol(3) = " // &
                  real2char (real (pol(3), kind=default)) // ".")
             write (*, "(A,4(1x,ES19.12))") "Antiparticle decay, q1 = ", q1
             write (*, "(A,4(1x,ES19.12))") "Antiparticle decay, p1 = ", p1
             write (*, "(A,4(1x,ES19.12))") "Antiparticle decay, p2 = ", p2
          end if
          if (debug_on) call msg_debug2 (D_TAUOLA, "TAUOLA is called here")
          call dexay (2,pol)
          if (IFPHOT == 1)  call photos (np2)
          is_swapped = .true.
       end if check_tau_sign
    end if

!!! TODO (Akiya Miyamoto, 25-march-2016)
!!!   In the following code, the tau helicity (polarization vector)
!!!   information is not stored in pyjets_spin_data(jtau)%helicity
!!!   and /HEPEV4/, because the tau polarization vector is determined
!!!   here in order to have a transverse spin correlation between
!!!   tau+ and tau-, but the decided polarization vectors are not
!!!   calculated here.  It would be possible to calculate them
!!!   from the polarimetric vectors, hh1 and hh2,  after
!!!   the end of the rejection loop.

    if (signal_is_pending ()) return

    if (trans_spin) then
       if (.not. tau_pol_vec) then
          pol1 = 0
          pol2 = 0
          if (pyr(0) .gt. 0.5) then
             pol1(3) = 1
             pol2(3) = -1
          else
             pol1(3) = -1
             pol2(3) = 1
          end if
          call dexay (1, pol1)
          call dexay (2, pol2)
       else
          !!! Decide polarimetric vector to have a spin correlation
          REJECTION: do
             call ranmar (rrr, 1)
             !!! tau+ decay
             call dekay (1, hh1)
             !!! tau- decay
             call dekay (2, hh2)
             wt = wthiggs (ifpseudo, hh1, hh2)
             if (rrr(1) .lt. wt)  exit REJECTION
          end do REJECTION
          ion = 0
          call dekay(11, hh1)
          call taupi0 (0, 1, ion)
          call dekay(12, hh2)
          call taupi0 (0, 2, ion)
       end if
       if (IFPHOT == 1)  call photos (im)
    end if

!!! **********************************************************
!!! Now copies /HEPEVT/ to /PYJETS/
!!! Higgs tau pair decay and single tau decay are treated
!!! separately.
!!! **********************************************************

    nproducts = 0
    np = nproducts

!!! =========================================================
!!! Higgs to tau pair decay case.
!!! =========================================================

    if (higgs_dec .and. trans_spin .and. jtau2 .gt. 0) then
       if (is_swapped) then
          !!! invert all momentum
          do i = n1, nhep
             phep(1:3,i) = - phep(1:3,i)
          end do
          do i = n2+1, nhep
             if (jmohep(1,i) == n1) then
                jmohep(1,i) = n2
                jmohep(2,i) = n2
             else if (jmohep(1,i) == n2) then
                jmohep(1,i) = n1
                jmohep(2,i) = n1
             end if
          end do
       end if

!!! Overwrite tau+ and tau- data in /PYJETS/, because tau+tau- momentum
!!! could have been changed due to photon emmision in Higgs --> tau+ tau-
!!! system. Their momentum should be boosted and rotate back to the lab frame
!!! in the calling routine, PYDCAY.

       if (is_swapped) then
          p(itau,:) = phep(:,3)
          k(itau,4) = jdahep(1,3) - n2 + n
          k(itau,5) = jdahep(2,3) - n2 + n
          p(jtau2,:) = phep(:,2)
          k(jtau2,4) = jdahep(1,2) - n2 + n
          k(jtau2,5) = jdahep(2,2) - n2 + n
       else
          p(itau,:) = phep(:,2)
          k(itau,4) = jdahep(1,2) - n2 + n
          k(itau,5) = jdahep(2,2) - n2 + n
          p(jtau2,:) = phep(:,3)
          k(jtau2,4) = jdahep(1,3) - n2 + n
          k(jtau2,5) = jdahep(2,3) - n2 + n
       end if
       k(itau,  1) = 11
       k(jtau2, 1) = 11
       k(itau,  3) = jorig
       k(jtau2, 3) = jorig

!!! TODO : Akiya Miyamoto, 12-April-2016
!!!   Reset daughter pointer of Higgs, because Higgs daughters
!!!   increase when photons are emitted. This may not work well if
!!!   additional particles exist after second tau.
       k(jorig,4) = itau
       k(jorig,5) = jtau2  ! jtau2 > jtau allways
       if (jdahep(2,1)-jdahep(1,1)+1 .gt. 2) then
          if (n .gt. jtau2) then
             write (msg_buffer, "(A)") &
                  "Tau decay routine do_dexay: necessary to update " // &
                  "index of Higgs daughter in order to include photons " // &
                  "produced by PHOTOS."
             call msg_message ()
             write (msg_buffer, "(A)") &
                  "Run continues without modifying the 2nd daughter pointer."
             call msg_message ()
          else
             k(jorig,5) = jdahep(2,1) - n2 + n
          end if
       end if

!!! Now, fill the information of tau daughters to /PYJETS/

       nproducts = 0
       loop_products_higgs: do i = n2+1, nhep
          nproducts = nproducts + 1
          p(n+nproducts,:) = phep(:,i)
          k(n+nproducts,2) = idhep(i)
          k(n+nproducts,3) = jmohep(1,i) - n2 + n
          if (isthep(i) == 1) then
             k(n+nproducts,1) = 1
             k(n+nproducts,4) = 0
             k(n+nproducts,5) = 0
          else
             k(n+nproducts,1) = 11
             k(n+nproducts,4) = jdahep(1,i) - n2 + n
             k(n+nproducts,5) = jdahep(2,i) - n2 + n
          end if
       end do loop_products_higgs

!!! ***************************************************************
!!! Single tau decay case.
!!! This case, parent tau daghter momentum is not over-wtitten
!!! ***************************************************************

    else
       loop_products_nohiggs: do i = n2+1, nhep
          nproducts = nproducts + 1
          p(n+nproducts,:) = phep(:,i)
          if (isthep(i) == 1) then
             k(n+nproducts,1) = 1
             k(n+nproducts,4) = 0
             k(n+nproducts,5) = 0
          else
             k(n+nproducts,1) = 11
             k(n+nproducts,4) = jdahep(1,i) - n2 + n
             k(n+nproducts,5) = jdahep(2,i) - n2 + n
          end if
          k(n+nproducts,2) = idhep(i)
          if (abs(idhep(jmohep(1,i))) .ne. 15) then
             k(n+nproducts,3) = jmohep(1,i) - n2 + n
          else
             k(n+nproducts,3) = itau
          end if
       end do loop_products_nohiggs
       k(itau,4) = jdahep(1,2) - n2 + n
       k(itau,5) = jdahep(2,2) - n2 + n
    end if

    if (debug2_active (D_TAUOLA)) then
      call msg_message ("TAUOLA interface: PYLIST at the end of do_dexay")
      n = n + nproducts
      call pylist(2)
      n = n - nproducts
    end if

  end subroutine do_dexay

  module subroutine taudec_settings_init (taudec_settings, var_list, model)
    class(taudec_settings_t), intent(out) :: taudec_settings
    type(var_list_t), intent(in) :: var_list
    class(model_data_t), intent(in) :: model
    type(field_data_t), pointer :: field
    taudec_settings%photos = &
         var_list%get_lval (var_str ("?ps_tauola_photos"))
    taudec_settings%transverse = &
         var_list%get_lval (var_str ("?ps_tauola_transverse"))
    taudec_settings%dec_rad_cor = &
         var_list%get_lval (var_str ("?ps_tauola_dec_rad_cor"))
    taudec_settings%dec_mode1 = &
         var_list%get_ival (var_str ("ps_tauola_dec_mode1"))
    taudec_settings%dec_mode2 = &
         var_list%get_ival (var_str ("ps_tauola_dec_mode2"))
    taudec_settings%mh = &
         var_list%get_rval (var_str ("ps_tauola_mh"))
    taudec_settings%mix_angle = &
         var_list%get_rval (var_str ("?ps_tauola_mix_angle"))
    taudec_settings%use_pol_vec = &
         var_list%get_lval (var_str ("?ps_tauola_pol_vector"))
    select case (char (model%get_name ()))
    case ("QCD", "Test")
       call msg_fatal ("taudec_settings_init: Model has no tau.")
    case default
       field => model%get_field_ptr (15)
       taudec_settings%mtau = field%get_mass ()
    end select
  end subroutine taudec_settings_init

  module subroutine taudec_settings_write (taudec_settings, unit)
    class(taudec_settings_t), intent(in) :: taudec_settings
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit); if (u<0)  return
    write (u, "(1x,A)")  "Tau decay settings:"
    call write_separator (u)
    write (u, "(3x,A,1x,L1)") &
         "ps_tauola_photos      = ", taudec_settings%photos
    write (u, "(3x,A,1x,L1)") &
         "ps_tauola_transverse  = ", taudec_settings%transverse
    write (u, "(3x,A,1x,L1)") &
         "ps_tauola_dec_rad_cor = ", taudec_settings%dec_rad_cor
    write (u, "(3x,A,1x,I2)") &
         "ps_tauola_dec_mode1   = ", taudec_settings%dec_mode1
    write (u, "(3x,A,1x,I2)") &
         "ps_tauola_dec_mode2   = ", taudec_settings%dec_mode2
    write (u, "(3x,A,1x,ES19.12)") &
         "ps_tauola_mh          = ", taudec_settings%mh
    write (u, "(3x,A,1x,ES19.12)") &
         "ps_tauola_mix_angle   = ", taudec_settings%mix_angle
    write (u, "(3x,A,1x,L1)") &
         "ps_tauola_use_pol_vec = ", taudec_settings%use_pol_vec
  end subroutine taudec_settings_write

  module function wo_tauola_get_helicity_mod (ip) result (the_helicity)
    integer, intent(in) :: ip
    integer :: the_helicity
    integer :: n, npad
    integer, dimension(4000,5) :: k
    double precision, dimension(4000,5) :: p
    double precision, dimension(4000,5) :: v
    common /pyjets/ n, npad, k, p, v
    save /pyjets/
    integer, dimension(200) :: mstu
    double precision, dimension(200) :: paru
    integer, dimension(200) :: mstj
    double precision, dimension(200) :: parj
    common /pydat1/ mstu, paru, mstj, parj
    save /pydat1/
    if (MSTJ(28) .NE. 2) then
       the_helicity=0
    else
      if (ip .le. 0 .or. ip .gt. n) then
        the_helicity = 0
      else
        the_helicity = int(pyjets_spin_data(ip)%helicity)
      end if
    end if
  end function wo_tauola_get_helicity_mod

  module subroutine wo_tauola_get_helicity (ip, the_helicity)
    integer, intent(in)  :: ip
    integer, intent(out) :: the_helicity
    the_helicity = wo_tauola_get_helicity_mod(ip)
    if (debug_active (D_TAUOLA)) then
       if (abs(the_helicity) .gt. 1) then
          write (msg_buffer, "(A,I0,A,I0,A)") &
               "Stored helicity information is wrong: ", the_helicity, &
               "for ip = ", ip, "."
          call msg_warning ()
       end if
    end if
  end subroutine wo_tauola_get_helicity

  module subroutine wo_tauola_init_call (taudec_settings)
    !!! Tauola initialization.
    !!! (default defined in rt_data)
    !!!  JAK1   ! (0) decay mode of first tau
    !!!  JAK2   ! (0) decay mode of second tau
    !!!  ITDKRC ! (1) switch on radiative corrections in decay
    !!!  IFPHOT ! (1) PHOTOS switch
    type(taudec_settings_t), intent(in) :: taudec_settings
    INTEGER JAK1, JAK2, JAKP, JAKM, KTOM
    COMMON /JAKI/ JAK1, JAK2, JAKP, JAKM, KTOM

    integer, dimension(200) :: MSTP
    double precision, dimension(200) :: PARP
    integer, dimension(200) :: MSTI
    double precision, dimension(200) :: PARI
    common /PYPARS/ MSTP, PARP, MSTI, PARI
    save /PYPARS/

    integer, dimension(200) :: MSTU
    double precision, dimension(200) :: PARU
    integer, dimension(200) :: MSTJ
    double precision, dimension(200) :: PARJ
    common /PYDAT1/ MSTU, PARU, MSTJ, PARJ
    save /PYDAT1/

    integer :: ITDKRC, IFPHOT
    double precision :: psi, betah
    double precision :: csc, ssc
    common /pseudocoup/ csc, ssc
    save /pseudocoup/

    integer, dimension(3) :: ion
    double precision, dimension(4) :: pol1x

    JAK1 = taudec_settings%dec_mode1
    JAK2 = taudec_settings%dec_mode2
    if (taudec_settings%dec_rad_cor) then
       ITDKRC = 1
    else
       ITDKRC = 0
    end if
    if (taudec_settings%photos) then
       IFPHOT = 1
    else
       IFPHOT = 0
    end if

    trans_spin = taudec_settings%transverse
    tau_pol_vec = taudec_settings%use_pol_vec

    psi = dble (taudec_settings%mix_angle * degree)
    betah = dble (sqrt (one - four * taudec_settings%mtau**2 / &
         taudec_settings%mh**2))
    csc = cos(psi) * betah
    ssc = sin(psi)

    if (signal_is_pending ()) return
    if (trans_spin) then
       if (mstj(39) .ne. 15) then
          call msg_warning ("wo_tauola_init_call: transverse spin " // &
               "correlation requested for H -> tau tau. Photon radiation " // &
               "from PYTHIA will be switched off.")
          mstj(39) = 15
       end if
    end if

    call phoini
    call inietc (JAK1, JAK2, ITDKRC, IFPHOT)
    call inimas
    call iniphx (0.01d0)
    call initdk
    ! !!! Deactivation of pi0 and eta decays: (1) means on, (0) off
    ion = 0
    call taupi0 (-1, 1, ion)
    call dekay (-1, pol1x)

    if (debug2_active (D_TAUOLA)) then
       call msg_debug2 (D_TAUOLA, "TAUOLA initialization")
       call taudec_settings%write ()
       call msg_debug2 (D_TAUOLA, " check if TAUOLA common block has been set")
       call msg_debug2 (D_TAUOLA, "Tau decay modes set")
       print *, " Tau decay modes: tau+(JAK1) = ", jak1, &
            " tau-(JAK2) = ", JAK2
       call msg_message ("   JAK =  0    : All decay mode")
       call msg_message ("   JAK =  1    : electron mode")
       call msg_message ("   JAK =  2    : muon mode")
       call msg_message ("   JAK =  3    : pion mode")
       call msg_message ("   JAK =  4    : rho mode")
       call msg_message ("   JAK =  5    : a1 mode")
       call msg_message ("   JAK =  6    : K mode")
       call msg_message ("   JAK =  7    : K* mode")
       call msg_message ("   JAK =  8-13 : n pion modes")
       call msg_message ("   JAK = 14-19 : K K pi and K pi pi modes")
       call msg_message ("   JAK = 20-21 : eta pi pi; gamma pi pi modes")
       call msg_debug2 (D_TAUOLA, "Radiative corrections in decay ON(1),Off(0)")
       print *, " ITDKRC = ", ITDKRC
       call msg_debug2 (D_TAUOLA, "PHOTOS switch: ON(1), OFF(0)")
       print *, " IFPHOT = ", IFPHOT
    end if
  end subroutine wo_tauola_init_call

end submodule tauola_interface_s

