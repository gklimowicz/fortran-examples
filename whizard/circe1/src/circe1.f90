! circe1.f90 -- canonical beam spectra for linear collider physics
! 
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     with contributions from
!     Christian Speckner <cnspeckn@googlemail.com>
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
! to the source 'circe1.nw'
module circe1
  use kinds

  implicit none
  private
    public :: circe
    public :: circes
    public :: circe1_params_t
    public :: circex
    public :: circel
    public :: circee
    public :: circeg
    public :: circgg
    public :: kirke
    public :: kirkee
    public :: kirkeg
    public :: kirkgg
    public :: girce
    public :: gircee
    public :: girceg
    public :: gircgg
    public :: girceb
    public :: circem

    public :: rng_type
  
    integer, parameter, public :: C1_ELECTRON = 11
    integer, parameter, public :: C1_POSITRON = -11
    integer, parameter, public :: C1_PHOTON = 22
    integer, parameter :: SBAND  = 1
    integer, parameter :: TESLA  = 2
    integer, parameter :: XBAND  = 3
    integer, parameter :: JLCNLC = 3
    integer, parameter :: SBNDEE = 4
    integer, parameter :: TESLEE = 5
    integer, parameter :: XBNDEE = 6
    integer, parameter :: NLCH   = 7
    integer, parameter :: ILC    = 8
    integer, parameter :: CLIC   = 9
    integer, parameter :: NACC   = 9
    integer :: e, r, ehi, elo

  integer, parameter, public :: MAGIC0 = 19040616  
  real(kind=double), parameter :: KIREPS = 1D-6

    type :: circe1_params_t
      real(kind=double) :: x1m = 0d0
      real(kind=double) :: x2m = 0d0
      real(kind=double) :: roots = 500D0
      real(kind=double) :: lumi
      real(kind=double) :: a1(0:7)
      real(kind=double) :: elect0, gamma0
      integer :: acc = TESLA
      integer :: ver = 0
      integer :: rev = 0
      integer :: chat = 1
      integer :: magic
    end type circe1_params_t

  type(circe1_params_t), public, save :: circe1_params

    type, abstract :: rng_type
     contains
       procedure(rng_generate), deferred :: generate
    end type rng_type
    

    abstract interface
      subroutine rng_proc (u)
        import :: double
        real(kind=double), intent(out) :: u
      end subroutine rng_proc
    end interface
    
    abstract interface
       subroutine rng_generate (rng_obj, u)
         import :: rng_type, double
         class(rng_type), intent(inout) :: rng_obj
         real(kind=double), intent(out) :: u
       end subroutine rng_generate
    end interface
    

contains

    function circe (x1, x2, p1, p2)
      real(kind=double) :: x1, x2
      integer :: p1, p2
      real(kind=double) :: circe
            if (circe1_params%magic .ne. MAGIC0) then
               call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
            endif
        circe = -1.0
        if (abs(p1) .eq. C1_ELECTRON) then
           if (abs(p2) .eq. C1_ELECTRON) then
              circe = circee (x1, x2)
           else if (p2 .eq. C1_PHOTON) then
              circe = circeg (x1, x2)
           end if
        else if (p1 .eq. C1_PHOTON) then
           if (abs(p2) .eq. C1_ELECTRON) then
              circe = circeg (x2, x1)
           else if (p2 .eq. C1_PHOTON) then
              circe = circgg (x1, x2)
           end if
        end if
    end function circe

    subroutine circes (xx1m, xx2m, xroots, xacc, xver, xrev, xchat)
      real(kind=double) :: xx1m, xx2m, xroots   
      integer :: xacc, xver, xrev, xchat
        character(len=60) :: msgbuf
        character(len=6), dimension(NACC) :: accnam
        integer :: ver34
      integer, parameter :: EINVAL = -2
      integer, parameter :: GEV090 = -1
      integer, parameter :: GEV170 = 0
      integer, parameter :: GEV350 = 1
      integer, parameter :: GEV500 = 2
      integer, parameter :: GEV800 = 3
      integer, parameter :: TEV1   = 4
      integer, parameter :: TEV16  = 5
      integer, parameter :: GEV250 = 6
      integer, parameter :: TEV12  = 7
      integer, parameter :: TEV15  = 8
      integer, parameter :: GEV200 = 9
      integer, parameter :: GEV230 = 10
      integer, parameter :: A1NEGY = 5
      integer, parameter :: A1NREV = 5
      integer :: i
      real(kind=double), dimension(A1NEGY,NACC,0:A1NREV), save :: xa1lum = 0
      real(kind=double), dimension(0:7,A1NEGY,NACC,0:A1NREV), save :: xa1 = 0
      integer, parameter :: A3NEGY = 5, A3NREV = 5
      real, dimension(A3NEGY,NACC,0:A3NREV), save :: xa3lum = -1
      real, dimension(0:7,A3NEGY,NACC,0:A3NREV), save :: xa3 = 0
      integer, parameter :: A5NEGY = 5, A5NREV = 1
      real, dimension(A5NEGY,NACC,0:A5NREV), save :: xa5lum
      real, dimension(0:7,A5NEGY,NACC,0:A5NREV), save :: xa5
      integer, parameter :: A6NEGY = 2, A6NREV = 1
      real, dimension(GEV090:A6NEGY,NACC,0:A6NREV), save :: xa6lum
      real, dimension(0:7,GEV090:A6NEGY,NACC,0:A6NREV), save :: xa6
      real(kind=double) :: eloval, ehival
      real(kind=double), parameter :: DELTAE = 0.5d0
      integer, parameter :: A7NEGY = TEV1, A7NREV = 1
      real, dimension(GEV090:A7NEGY,NACC,0:A7NREV), save :: xa7lum
      real, dimension(0:7,GEV090:A7NEGY,NACC,0:A7NREV), save :: xa7
      integer, parameter :: A8NEGY = TEV1, A8NREV = 1
      real, dimension(GEV090:A8NEGY,NACC,0:A8NREV), save :: xa8lum
      real, dimension(0:7,GEV090:A8NEGY,NACC,0:A8NREV), save :: xa8
      integer, parameter ::  A9NEGY = TEV15, A9NREV = 1
      real, dimension(GEV090:A9NEGY,NACC,0:A9NREV) :: xa9lum
      real, dimension(0:7,GEV090:A9NEGY,NACC,0:A9NREV) :: xa9
      integer, parameter ::  A10NEGY = GEV230, A10NREV = 1
      real, dimension(GEV090:A10NEGY,ILC:ILC,0:A10NREV) :: xa10lum
      real, dimension(0:7,GEV090:A10NEGY,ILC:ILC,0:A10NREV) :: xa10
        data accnam(SBAND)  /'SBAND'/
        data accnam(TESLA)  /'TESLA'/
        data accnam(JLCNLC) /'JLCNLC'/
        data accnam(SBNDEE) /'SBNDEE'/
        data accnam(TESLEE) /'TESLEE'/
        data accnam(XBNDEE) /'XBNDEE'/
        data accnam(NLCH) /'NLC H'/
        data accnam(ILC) /'ILC'/
        data accnam(CLIC) /'CLIC'/
      xa1lum(GEV500,SBAND,1) = 5.212299E+01 
      xa1(0:7,GEV500,SBAND,1) = (/ &
         .39192E+00, .66026E+00, .11828E+02,-.62543E+00, &
         .52292E+00,-.69245E+00, .14983E+02, .65421E+00 /)
      xa1lum(GEV500,TESLA,1) = 6.066178E+01
      xa1(0:7,GEV500,TESLA,1) = (/ &
         .30196E+00, .12249E+01, .21423E+02,-.57848E+00, &
         .68766E+00,-.69788E+00, .23121E+02, .78399E+00 /)
      xa1lum(GEV500,XBAND,1) = 5.884699E+01
      xa1(0:7,GEV500,XBAND,1) = (/ & 
         .48594E+00, .52435E+00, .83585E+01,-.61347E+00, &
         .30703E+00,-.68804E+00, .84109E+01, .44312E+00 /)
      xa1lum(TEV1,SBAND,1) = 1.534650E+02
      xa1(0:7,TEV1,SBAND,1) = (/ &
          .24399E+00, .87464E+00, .66751E+01,-.56808E+00, &
          .59295E+00,-.68921E+00, .94232E+01, .83351E+00 /)
      xa1lum(TEV1,TESLA,1) = 1.253381E+03
      xa1(0:7,TEV1,TESLA,1) = (/ &
          .39843E+00, .70097E+00, .11602E+02,-.61061E+00, &
          .40737E+00,-.69319E+00, .14800E+02, .51382E+00 /)
      xa1lum(TEV1,XBAND,1) = 1.901783E+02 
      xa1(0:7,TEV1,XBAND,1) = (/ &
          .32211E+00,   .61798E+00,   .28298E+01,  -.54644E+00, &
          .45674E+00,  -.67301E+00,   .41703E+01,   .74536E+00 /)
      xa1lum(GEV350,1:NACC,1) = NACC * (-1d0)
      xa1lum(GEV800,1:NACC,1) = NACC * (-1d0)
      xa1lum(GEV500,SBNDEE:NACC,1) = 4 * (-1d0)
      xa1lum(TEV1,SBNDEE:NACC,1) = 4 * (-1d0)
      xa1lum(TEV16,1:NACC,1) = 7 * (-1d0)
      xa1lum(GEV500,SBAND,2) = .31057E+02
      xa1(0:7,GEV500,SBAND,2) = (/ &
          .38504E+00, .79723E+00, .14191E+02,-.60456E+00, &
          .53411E+00,-.68873E+00, .15105E+02, .65151E+00 /)
      xa1lum(TEV1,SBAND,2) = .24297E+03
      xa1(0:7,TEV1,SBAND,2) = (/ & 
          .24374E+00, .89466E+00, .70242E+01,-.56754E+00, &
          .60910E+00,-.68682E+00, .96083E+01, .83985E+00 /)
      xa1lum(GEV350,TESLA,2) = .73369E+02
      xa1(0:7,GEV350,TESLA,2) = (/ &
          .36083E+00, .12819E+01, .37880E+02,-.59492E+00, &
          .69109E+00,-.69379E+00, .40061E+02, .65036E+00 /)
      xa1lum(GEV500,TESLA,2) = .10493E+03
      xa1(0:7,GEV500,TESLA,2) = (/ &
          .29569E+00, .11854E+01, .21282E+02,-.58553E+00, &
          .71341E+00,-.69279E+00, .24061E+02, .77709E+00 /)
      xa1lum(GEV800,TESLA,2) = .28010E+03
      xa1(0:7,GEV800,TESLA,2) = (/ &
          .22745E+00, .11265E+01, .10483E+02,-.55711E+00, &
          .69579E+00,-.69068E+00, .13093E+02, .89605E+00 /)
      xa1lum(TEV1,TESLA,2) = .10992E+03
      xa1(0:7,TEV1,TESLA,2) = (/ & 
          .40969E+00, .66105E+00, .11972E+02,-.62041E+00, &
          .40463E+00,-.69354E+00, .14669E+02, .51281E+00 /)
      xa1lum(GEV500,XBAND,2) = .35689E+02
      xa1(0:7,GEV500,XBAND,2) = (/ &
          .48960E+00, .46815E+00, .75249E+01,-.62769E+00, &
          .30341E+00,-.68754E+00, .85545E+01, .43453E+00 /)
      xa1lum(TEV1,XBAND,2) = .11724E+03
      xa1(0:7,TEV1,XBAND,2) =  (/ &
          .31939E+00, .62415E+00, .30763E+01,-.55314E+00, &
          .45634E+00,-.67089E+00, .41529E+01, .73807E+00 /)
      xa1lum(GEV350,SBAND,2) = -1d0
      xa1lum(GEV350,XBAND,2) = -1d0
      xa1lum(GEV800,SBAND,2) = -1d0
      xa1lum(GEV800,XBAND,2) = -1d0
      xa1lum(GEV350,SBNDEE:NACC,2) = 4 * (-1d0)
      xa1lum(GEV500,SBNDEE:NACC,2) = 4 * (-1d0)
      xa1lum(GEV800,SBNDEE:NACC,2) = 4 * (-1d0)
      xa1lum(TEV1,SBNDEE:NACC,2) = 4 * (-1d0)
      xa1lum(TEV16,1:NACC,2) = 7 * (-1d0)
      xa1lum(GEV500,SBAND,3) = .31469E+02
      xa1(0:7,GEV500,SBAND,3) = (/ & 
          .38299E+00, .72035E+00, .12618E+02,-.61611E+00, &
          .51971E+00,-.68960E+00, .15066E+02, .63784E+00 /)
      xa1lum(TEV1,SBAND,3) = .24566E+03
      xa1(0:7,TEV1,SBAND,3) = (/ & 
          .24013E+00, .95763E+00, .69085E+01,-.55151E+00, &
          .59497E+00,-.68622E+00, .94494E+01, .82158E+00 /)
      xa1lum(GEV350,TESLA,3) = .74700E+02
      xa1(0:7,GEV350,TESLA,3) = (/ &
          .34689E+00, .12484E+01, .33720E+02,-.59523E+00, &
          .66266E+00,-.69524E+00, .38488E+02, .63775E+00 /)
      xa1lum(GEV500,TESLA,3) = .10608E+03
      xa1(0:7,GEV500,TESLA,3) = (/ &
          .28282E+00, .11700E+01, .19258E+02,-.58390E+00, &
          .68777E+00,-.69402E+00, .23638E+02, .75929E+00 /)
      xa1lum(GEV800,TESLA,3) = .28911E+03
      xa1(0:7,GEV800,TESLA,3) = (/ &
          .21018E+00, .12039E+01, .96763E+01,-.54024E+00, &
          .67220E+00,-.69083E+00, .12733E+02, .87355E+00 /)
      xa1lum(TEV1,TESLA,3) = .10936E+03
      xa1(0:7,TEV1,TESLA,3) = (/ &
          .41040E+00, .68099E+00, .11610E+02,-.61237E+00, &
          .40155E+00,-.69073E+00, .14698E+02, .49989E+00 /)
      xa1lum(GEV500,XBAND,3) = .36145E+02
      xa1(0:7,GEV500,XBAND,3) = (/ &
          .51285E+00, .45812E+00, .75135E+01,-.62247E+00, &
          .30444E+00,-.68530E+00, .85519E+01, .43062E+00 /)
      xa1lum(TEV1,XBAND,3) = .11799E+03
      xa1(0:7,TEV1,XBAND,3) = (/ &
          .31241E+00, .61241E+00, .29938E+01,-.55848E+00, &
          .44801E+00,-.67116E+00, .41119E+01, .72753E+00 /)
      xa1lum(GEV350,SBAND,3) = -1d0
      xa1lum(GEV350,XBAND,3) = -1d0
      xa1lum(GEV800,SBAND,3) = -1d0
      xa1lum(GEV800,XBAND,3) = -1d0
      xa1lum(GEV350,SBNDEE:NACC,3) = 4 * (-1d0)
      xa1lum(GEV500,SBNDEE:NACC,3) = 4 * (-1d0)
      xa1lum(GEV800,SBNDEE:NACC,3) = 4 * (-1d0)
      xa1lum(TEV1,SBNDEE:NACC,3) = 4 * (-1d0)
      xa1lum(TEV16,1:NACC,3) = 7 * (-1d0)
      xa1lum(GEV500,SBAND,4) = .31528E+02
      xa1(0:7,GEV500,SBAND,4) = (/ & 
          .38169E+00, .73949E+00, .12543E+02,-.61112E+00, &
          .51256E+00,-.69009E+00, .14892E+02, .63314E+00 /)
      xa1lum(TEV1,SBAND,4) = .24613E+03
      xa1(0:7,TEV1,SBAND,4) = (/ &
          .24256E+00, .94117E+00, .66775E+01,-.55160E+00, &
          .57484E+00,-.68891E+00, .92271E+01, .81162E+00 /)
      xa1lum(GEV350,TESLA,4) = .74549E+02
      xa1(0:7,GEV350,TESLA,4) = (/ &
          .34120E+00, .12230E+01, .32932E+02,-.59850E+00, &
          .65947E+00,-.69574E+00, .38116E+02, .63879E+00 /)
      xa1lum(GEV500,TESLA,4) = .10668E+03
      xa1(0:7,GEV500,TESLA,4) = (/ &
          .28082E+00, .11074E+01, .18399E+02,-.59118E+00, &
          .68880E+00,-.69375E+00, .23463E+02, .76073E+00 /)
      xa1lum(GEV800,TESLA,4) = .29006E+03
      xa1(0:7,GEV800,TESLA,4) = (/ &
          .21272E+00, .11443E+01, .92564E+01,-.54657E+00, &
          .66799E+00,-.69137E+00, .12498E+02, .87571E+00 /)
      xa1lum(TEV1,TESLA,4) = .11009E+03
      xa1(0:7,TEV1,TESLA,4) = (/ &
          .41058E+00, .64745E+00, .11271E+02,-.61996E+00, &
          .39801E+00,-.69150E+00, .14560E+02, .49924E+00 /)
      xa1lum(GEV500,XBAND,4) = .36179E+02
      xa1(0:7,GEV500,XBAND,4) = (/ &
          .51155E+00, .43313E+00, .70446E+01,-.63003E+00, &
          .29449E+00,-.68747E+00, .83489E+01, .42458E+00 /)
      xa1lum(TEV1,XBAND,4) = .11748E+03
      xa1(0:7,TEV1,XBAND,4) = (/ &
          .32917E+00, .54322E+00, .28493E+01,-.57959E+00, &
          .39266E+00,-.68217E+00, .38475E+01, .68478E+00 /)
      xa1lum(GEV350,SBAND,4) = -1d0
      xa1lum(GEV350,XBAND,4) = -1d0
      xa1lum(GEV800,SBAND,4) = -1d0
      xa1lum(GEV800,XBAND,4) = -1d0
      xa1lum(GEV350,SBNDEE:NACC,4) = 4 * (-1d0)
      xa1lum(GEV500,SBNDEE:NACC,4) = 4 * (-1d0)
      xa1lum(GEV800,SBNDEE:NACC,4) = 4 * (-1d0)
      xa1lum(TEV1,SBNDEE:NACC,4) = 4 * (-1d0)
    xa1lum(TEV16,1:NACC,4) = 7 * (-1d0)
      xa1lum(GEV350,SBAND,5) = 0.21897E+02
      xa1(0:7,GEV350,SBAND,5) = (/ &
          0.57183E+00, 0.53877E+00, 0.19422E+02,-0.63064E+00, &
          0.49112E+00,-0.69109E+00, 0.24331E+02, 0.52718E+00 /)
      xa1lum(GEV500,SBAND,5) = 0.31383E+02
      xa1(0:7,GEV500,SBAND,5) = (/ &
          0.51882E+00, 0.49915E+00, 0.11153E+02,-0.63017E+00, &
          0.50217E+00,-0.69113E+00, 0.14935E+02, 0.62373E+00 /)
      xa1lum(GEV800,SBAND,5) = 0.95091E+02
      xa1(0:7,GEV800,SBAND,5) = (/ &
          0.47137E+00, 0.46150E+00, 0.56562E+01,-0.61758E+00, &
          0.46863E+00,-0.68897E+00, 0.85876E+01, 0.67577E+00 /)
      xa1lum(TEV1,SBAND,5) = 0.11900E+03
      xa1(0:7,TEV1,SBAND,5) = (/ &
          0.43956E+00, 0.45471E+00, 0.42170E+01,-0.61180E+00, &
          0.48711E+00,-0.68696E+00, 0.67145E+01, 0.74551E+00 /)
      xa1lum(TEV16,SBAND,5) = 0.11900E+03
      xa1(0:7,TEV16,SBAND,5) = (/ &
          0.43956E+00, 0.45471E+00, 0.42170E+01,-0.61180E+00, &
          0.48711E+00,-0.68696E+00, 0.67145E+01, 0.74551E+00 /)
      xa1lum(GEV350,TESLA,5) = 0.97452E+02
      xa1(0:7,GEV350,TESLA,5) = (/ &
          0.39071E+00, 0.84996E+00, 0.17614E+02,-0.60609E+00, &
          0.73920E+00,-0.69490E+00, 0.28940E+02, 0.77286E+00 /)
      xa1lum(GEV500,TESLA,5) = 0.10625E+03
      xa1(0:7,GEV500,TESLA,5) = (/ &
          0.42770E+00, 0.71457E+00, 0.15284E+02,-0.61664E+00, &
          0.68166E+00,-0.69208E+00, 0.24165E+02, 0.73806E+00 /)
      xa1lum(GEV800,TESLA,5) = 0.17086E+03
      xa1(0:7,GEV800,TESLA,5) = (/ &
          0.36025E+00, 0.69118E+00, 0.76221E+01,-0.59440E+00, &
          0.71269E+00,-0.69077E+00, 0.13117E+02, 0.91780E+00 /)
      xa1lum(TEV1,TESLA,5) = 0.21433E+03
      xa1(0:7,TEV1,TESLA,5) = (/ &
          0.33145E+00, 0.67075E+00, 0.55438E+01,-0.58468E+00, &
          0.72503E+00,-0.69084E+00, 0.99992E+01, 0.10112E+01 /)
      xa1lum(TEV16,TESLA,5) = 0.34086E+03
      xa1(0:7,TEV16,TESLA,5) = (/ &
          0.49058E+00, 0.42609E+00, 0.50550E+01,-0.61867E+00, &
          0.39225E+00,-0.68916E+00, 0.75514E+01, 0.58754E+00 /)
      xa1lum(GEV350,XBAND,5) = 0.31901E+02
      xa1(0:7,GEV350,XBAND,5) = (/ &
          0.65349E+00, 0.31752E+00, 0.94342E+01,-0.64291E+00, &
          0.30364E+00,-0.68989E+00, 0.11446E+02, 0.40486E+00 /)
      xa1lum(GEV500,XBAND,5) = 0.36386E+02
      xa1(0:7,GEV500,XBAND,5) = (/ &
          0.65132E+00, 0.28728E+00, 0.69853E+01,-0.64440E+00, &
          0.28736E+00,-0.68758E+00, 0.83227E+01, 0.41492E+00 /)
      xa1lum(GEV800,XBAND,5) = 0.10854E+03
      xa1(0:7,GEV800,XBAND,5) = (/ &
          0.49478E+00, 0.36221E+00, 0.30116E+01,-0.61548E+00, &
          0.39890E+00,-0.68418E+00, 0.45183E+01, 0.67243E+00 /)
      xa1lum(TEV1,XBAND,5) = 0.11899E+03
      xa1(0:7,TEV1,XBAND,5) = (/ &
          0.49992E+00, 0.34299E+00, 0.26184E+01,-0.61584E+00, &
          0.38450E+00,-0.68342E+00, 0.38589E+01, 0.67408E+00 /)
      xa1lum(TEV16,XBAND,5) = 0.13675E+03
      xa1(0:7,TEV16,XBAND,5) = (/ &
          0.50580E+00, 0.30760E+00, 0.18339E+01,-0.61421E+00, &
          0.35233E+00,-0.68315E+00, 0.26708E+01, 0.67918E+00 /)
      xa1lum(GEV500,SBNDEE,0) = .92914E+01
      xa1(0:7,GEV500,SBNDEE,0) = (/ &
          .34866E+00, .78710E+00, .10304E+02,-.59464E+00, &
          .40234E+00,-.69741E+00, .20645E+02, .47274E+00 /)
      xa1lum(TEV1,SBNDEE,0) = .45586E+02
      xa1(0:7,TEV1,SBNDEE,0) = (/ &
          .21084E+00, .99168E+00, .54407E+01,-.52851E+00, &
          .47493E+00,-.69595E+00, .12480E+02, .64027E+00 /)
      xa1lum(GEV350,TESLEE,0) = .15175E+02
      xa1(0:7,GEV350,TESLEE,0) = (/ &
          .33093E+00, .11137E+01, .25275E+02,-.59942E+00, &
          .49623E+00,-.70403E+00, .60188E+02, .44637E+00 /)
      xa1lum(GEV500,TESLEE,0) = .21622E+02
      xa1(0:7,GEV500,TESLEE,0) = (/ &
          .27175E+00, .10697E+01, .14858E+02,-.58418E+00, &
          .50824E+00,-.70387E+00, .36129E+02, .53002E+00 /)
      xa1lum(GEV800,TESLEE,0) = .43979E+02
      xa1(0:7,GEV800,TESLEE,0) = (/ &
          .22994E+00, .10129E+01, .81905E+01,-.55751E+00, &
          .46551E+00,-.70461E+00, .19394E+02, .58387E+00 /)
      xa1lum(TEV1,TESLEE,0) = .25465E+02
      xa1(0:7,TEV1,TESLEE,0) = (/ &
          .37294E+00, .67522E+00, .87504E+01,-.60576E+00, &
          .35095E+00,-.69821E+00, .18526E+02, .42784E+00 /)
      xa1lum(GEV500,XBNDEE,0) = .13970E+02
      xa1(0:7,GEV500,XBNDEE,0) = (/ &
          .47296E+00, .46800E+00, .58897E+01,-.61689E+00, &
          .27181E+00,-.68923E+00, .10087E+02, .37462E+00 /)
      xa1lum(TEV1,XBNDEE,0) = .41056E+02
      xa1(0:7,TEV1,XBNDEE,0) = (/ & 
          .27965E+00, .74816E+00, .27415E+01,-.50491E+00, &
          .38320E+00,-.67945E+00, .47506E+01, .62218E+00 /)
      xa1lum(GEV350,SBNDEE,0) = -1d0  
      xa1lum(GEV350,XBNDEE,0) = -1d0  
      xa1lum(GEV800,SBNDEE,0) = -1d0  
      xa1lum(GEV800,XBNDEE,0) = -1d0  
      xa1lum(GEV500,SBAND,0) = .31528E+02
      xa1(0:7,GEV500,SBAND,0) = (/ &
          .38169E+00, .73949E+00, .12543E+02,-.61112E+00, &
          .51256E+00,-.69009E+00, .14892E+02, .63314E+00 /)
      xa1lum(TEV1,SBAND,0) = .24613E+03
      xa1(0:7,TEV1,SBAND,0) = (/ &
          .24256E+00, .94117E+00, .66775E+01,-.55160E+00, &
          .57484E+00,-.68891E+00, .92271E+01, .81162E+00 /)
      xa1lum(GEV350,TESLA,0) = .74549E+02
      xa1(0:7,GEV350,TESLA,0) = (/ &
          .34120E+00, .12230E+01, .32932E+02,-.59850E+00, &
          .65947E+00,-.69574E+00, .38116E+02, .63879E+00 /)
      xa1lum(GEV500,TESLA,0) = .10668E+03
      xa1(0:7,GEV500,TESLA,0) = (/ &
          .28082E+00, .11074E+01, .18399E+02,-.59118E+00, &
          .68880E+00,-.69375E+00, .23463E+02, .76073E+00 /)
      xa1lum(GEV800,TESLA,0) = .29006E+03
      xa1(0:7,GEV800,TESLA,0) = (/ &
          .21272E+00, .11443E+01, .92564E+01,-.54657E+00, &
          .66799E+00,-.69137E+00, .12498E+02, .87571E+00 /)
      xa1lum(TEV1,TESLA,0) = .11009E+03
      xa1(0:7,TEV1,TESLA,0) = (/ &
          .41058E+00, .64745E+00, .11271E+02,-.61996E+00, &
          .39801E+00,-.69150E+00, .14560E+02, .49924E+00 /)
      xa1lum(GEV500,XBAND,0) = .36179E+02
      xa1(0:7,GEV500,XBAND,0) = (/ &
          .51155E+00, .43313E+00, .70446E+01,-.63003E+00, &
          .29449E+00,-.68747E+00, .83489E+01, .42458E+00 /)
      xa1lum(TEV1,XBAND,0) = .11748E+03
      xa1(0:7,TEV1,XBAND,0) = (/ &
          .32917E+00, .54322E+00, .28493E+01,-.57959E+00, &
          .39266E+00,-.68217E+00, .38475E+01, .68478E+00 /)
      xa1lum(GEV350,SBAND,0) = -1d0
      xa1lum(GEV350,XBAND,0) = -1d0
      xa1lum(GEV800,SBAND,0) = -1d0
      xa1lum(GEV800,XBAND,0) = -1d0
      xa3lum(GEV800,TESLA,3) = .17196E+03
      xa3(0:7,GEV800,TESLA,3) = (/ &
          .21633E+00, .11333E+01, .95928E+01,-.55095E+00, &
          .73044E+00,-.69101E+00, .12868E+02, .94737E+00 /)
      xa3lum(GEV800,TESLA, 4) = .16408E+03
      xa3(0:7,GEV800,TESLA, 4) = (/ &
          .41828E+00, .72418E+00, .14137E+02,-.61189E+00, &
          .36697E+00,-.69205E+00, .17713E+02, .43583E+00 /)
      xa3lum(GEV350,TESLA,5) = 0.66447E+02
      xa3(0:7,GEV350,TESLA,5) = (/ &
          0.69418E+00, 0.50553E+00, 0.48430E+02,-0.63911E+00, &
          0.34074E+00,-0.69533E+00, 0.55502E+02, 0.29397E+00 /)
      xa3lum(GEV500,TESLA,5) = 0.95241E+02
      xa3(0:7,GEV500,TESLA,5) = (/ &
          0.64882E+00, 0.45462E+00, 0.27103E+02,-0.64535E+00, &
          0.35101E+00,-0.69467E+00, 0.33658E+02, 0.35024E+00 /)
      xa3lum(GEV800,TESLA,5) = 0.16974E+03
      xa3(0:7,GEV800,TESLA,5) = (/ &
          0.58706E+00, 0.43771E+00, 0.13422E+02,-0.63804E+00, &
          0.35541E+00,-0.69467E+00, 0.17528E+02, 0.43051E+00 /)
      xa3lum(TEV1,TESLA,5) = 0.21222E+03
      xa3(0:7,TEV1,TESLA,5) = (/ &
          0.55525E+00, 0.42577E+00, 0.96341E+01,-0.63587E+00, &
          0.36448E+00,-0.69365E+00, 0.13161E+02, 0.47715E+00 /)
      xa3lum(TEV16,TESLA,5) = 0.34086E+03
      xa3(0:7,TEV16,TESLA,5) = (/ &
          0.49058E+00, 0.42609E+00, 0.50550E+01,-0.61867E+00, &
          0.39225E+00,-0.68916E+00, 0.75514E+01, 0.58754E+00 /)
      xa3lum(GEV350,TESLA,0) = 0.66447E+02
      xa3(0:7,GEV350,TESLA,0) = (/ &
          0.69418E+00, 0.50553E+00, 0.48430E+02,-0.63911E+00, &
          0.34074E+00,-0.69533E+00, 0.55502E+02, 0.29397E+00 /)
      xa3lum(GEV500,TESLA,0) = 0.95241E+02
      xa3(0:7,GEV500,TESLA,0) = (/ &
          0.64882E+00, 0.45462E+00, 0.27103E+02,-0.64535E+00, &
          0.35101E+00,-0.69467E+00, 0.33658E+02, 0.35024E+00 /)
      xa3lum(GEV800,TESLA,0) = 0.16974E+03
      xa3(0:7,GEV800,TESLA,0) = (/ &
          0.58706E+00, 0.43771E+00, 0.13422E+02,-0.63804E+00, &
          0.35541E+00,-0.69467E+00, 0.17528E+02, 0.43051E+00 /)
      xa3lum(TEV1,TESLA,0) = 0.21222E+03
      xa3(0:7,TEV1,TESLA,0) = (/ &
          0.55525E+00, 0.42577E+00, 0.96341E+01,-0.63587E+00, &
          0.36448E+00,-0.69365E+00, 0.13161E+02, 0.47715E+00 /)
      xa3lum(TEV16,TESLA,0) = 0.34086E+03
      xa3(0:7,TEV16,TESLA,0) = (/ &
          0.49058E+00, 0.42609E+00, 0.50550E+01,-0.61867E+00, &
          0.39225E+00,-0.68916E+00, 0.75514E+01, 0.58754E+00 /)
      xa5lum(GEV350,TESLA,1) = -1.0
      xa5lum(GEV500,TESLA,1) = 0.33980E+03
      xa5(0:7,GEV500,TESLA,1) = (/ &
         0.49808E+00, 0.54613E+00, 0.12287E+02,-0.62756E+00, &
         0.42817E+00,-0.69120E+00, 0.17067E+02, 0.51143E+00 /)
      xa5lum(GEV800,TESLA,1) = 0.35936E+03
      xa5(0:7,GEV800,TESLA,1) = (/ &
         0.58751E+00, 0.43128E+00, 0.13324E+02,-0.64006E+00, & 
         0.30682E+00,-0.69235E+00, 0.16815E+02, 0.37078E+00 /)
      xa5lum(TEV1, TESLA,1) = -1.0
      xa5lum(TEV16,TESLA,1) = -1.0
      xa5lum(GEV350,TESLA,0) = -1.0
      xa5lum(GEV500,TESLA,0) = 0.33980E+03
      xa5(0:7,GEV500,TESLA,0) = (/ &
         0.49808E+00, 0.54613E+00, 0.12287E+02,-0.62756E+00, &
         0.42817E+00,-0.69120E+00, 0.17067E+02, 0.51143E+00 /)
      xa5lum(GEV800,TESLA,0) = 0.35936E+03
      xa5(0:7,GEV800,TESLA,0) = (/ &
         0.58751E+00, 0.43128E+00, 0.13324E+02,-0.64006E+00, &
         0.30682E+00,-0.69235E+00, 0.16815E+02, 0.37078E+00 /)
      xa5lum(TEV1, TESLA,0) = -1.0
      xa5lum(TEV16,TESLA,0) = -1.0
      xa6lum(GEV090,TESLA,1) = 0.62408E+02
      xa6(0:7,GEV090,TESLA,1) = (/ &
         0.72637E+00, 0.75534E+00, 0.18180E+03,-0.63426E+00, &
         0.36829E+00,-0.69653E+00, 0.18908E+03, 0.22157E+00 /)
      xa6lum(GEV170,TESLA,1) = 0.11532E+02
      xa6(0:7,GEV170,TESLA,1) = (/ &
         0.65232E+00, 0.67249E+00, 0.66862E+02,-0.63315E+00, &
         0.38470E+00,-0.69477E+00, 0.75120E+02, 0.30162E+00 /)
      xa6lum(GEV350,TESLA,1) = 0.24641E+03
      xa6(0:7,GEV350,TESLA,1) = (/ &
         0.54610E+00, 0.59105E+00, 0.20297E+02,-0.62747E+00, &
         0.41588E+00,-0.69188E+00, 0.26345E+02, 0.43818E+00 /)
      xa6lum(GEV500,TESLA,1) = 0.30340E+03
      xa6(0:7,GEV500,TESLA,1) = (/ &
         0.52744E+00, 0.52573E+00, 0.13895E+02,-0.63145E+00, &
         0.40824E+00,-0.69150E+00, 0.18645E+02, 0.47585E+00 /)
      xa6lum(GEV090,TESLA,0) = 0.62408E+02
      xa6(0:7,GEV090,TESLA,0) = (/ &
         0.72637E+00, 0.75534E+00, 0.18180E+03,-0.63426E+00, &
         0.36829E+00,-0.69653E+00, 0.18908E+03, 0.22157E+00 /)
      xa6lum(GEV170,TESLA,0) = 0.11532E+02
      xa6(0:7,GEV170,TESLA,0) = (/ &
         0.65232E+00, 0.67249E+00, 0.66862E+02,-0.63315E+00, &
         0.38470E+00,-0.69477E+00, 0.75120E+02, 0.30162E+00 /)
      xa6lum(GEV350,TESLA,0) = 0.24641E+03
      xa6(0:7,GEV350,TESLA,0) = (/ &
         0.54610E+00, 0.59105E+00, 0.20297E+02,-0.62747E+00, &
         0.41588E+00,-0.69188E+00, 0.26345E+02, 0.43818E+00 /)
      xa6lum(GEV500,TESLA,0) = 0.30340E+03
      xa6(0:7,GEV500,TESLA,0) = (/ &
         0.52744E+00, 0.52573E+00, 0.13895E+02,-0.63145E+00, &
         0.40824E+00,-0.69150E+00, 0.18645E+02, 0.47585E+00 /)
      xa7lum(GEV090,TESLA,1) = 0.62408E+02
      xa7(0:7,GEV090,TESLA,1) = (/ &
         0.72637E+00, 0.75534E+00, 0.18180E+03,-0.63426E+00, &
         0.36829E+00,-0.69653E+00, 0.18908E+03, 0.22157E+00 /)
      xa7lum(GEV170,TESLA,1) = 0.11532E+02
      xa7(0:7,GEV170,TESLA,1) = (/ &
         0.65232E+00, 0.67249E+00, 0.66862E+02,-0.63315E+00, &
         0.38470E+00,-0.69477E+00, 0.75120E+02, 0.30162E+00 /)
      xa7lum(GEV350,TESLA,1) = 0.24641E+03
      xa7(0:7,GEV350,TESLA,1) = (/ &
         0.54610E+00, 0.59105E+00, 0.20297E+02,-0.62747E+00, &
         0.41588E+00,-0.69188E+00, 0.26345E+02, 0.43818E+00 /)
      xa7lum(GEV500,TESLA,1) = 0.34704E+03
      xa7(0:7,GEV500,TESLA,1) = (/ &
         0.51288E+00, 0.49025E+00, 0.99716E+01,-0.62850E+00, &
         0.41048E+00,-0.69065E+00, 0.13922E+02, 0.51902E+00 /)
      xa7lum(GEV800,TESLA,1) = 0.57719E+03
      xa7(0:7,GEV800,TESLA,1) = (/ &
         0.52490E+00, 0.42573E+00, 0.69069E+01,-0.62649E+00, &
         0.32380E+00,-0.68958E+00, 0.93819E+01, 0.45671E+00 /)
      xa7lum(TEV1,TESLA,1) = -1.0
      xa7lum(GEV090,JLCNLC,1) = -1.0
      xa7lum(GEV170,JLCNLC,1) = -1.0
      xa7lum(GEV350,JLCNLC,1) = -1.0
      xa7lum(GEV500,JLCNLC,1) = 0.63039E+02
      xa7(0:7,GEV500,JLCNLC,1) = (/ &
          0.58967E+00, 0.34035E+00, 0.63631E+01,-0.63683E+00, &
          0.33383E+00,-0.68803E+00, 0.81005E+01, 0.48702E+00 /)
      xa7lum(TEV1,JLCNLC,1) = 0.12812E+03
      xa7(0:7,TEV1,JLCNLC,1) = (/ &
          0.50222E+00, 0.33773E+00, 0.25681E+01,-0.61711E+00, &
          0.36826E+00,-0.68335E+00, 0.36746E+01, 0.65393E+00 /)
      xa7lum(GEV090,TESLA,0) = 0.62408E+02
      xa7(0:7,GEV090,TESLA,0) = (/ &
          0.72637E+00, 0.75534E+00, 0.18180E+03,-0.63426E+00, &
          0.36829E+00,-0.69653E+00, 0.18908E+03, 0.22157E+00 /)
      xa7lum(GEV170,TESLA,0) = 0.11532E+02
      xa7(0:7,GEV170,TESLA,0) = (/ &
          0.65232E+00, 0.67249E+00, 0.66862E+02,-0.63315E+00, &
          0.38470E+00,-0.69477E+00, 0.75120E+02, 0.30162E+00 /)
      xa7lum(GEV350,TESLA,0) = 0.24641E+03
      xa7(0:7,GEV350,TESLA,0) = (/ &
          0.54610E+00, 0.59105E+00, 0.20297E+02,-0.62747E+00, &
          0.41588E+00,-0.69188E+00, 0.26345E+02, 0.43818E+00 /)
      xa7lum(GEV500,TESLA,0) = 0.34704E+03
      xa7(0:7,GEV500,TESLA,0) = (/ &
          0.51288E+00, 0.49025E+00, 0.99716E+01,-0.62850E+00, &
          0.41048E+00,-0.69065E+00, 0.13922E+02, 0.51902E+00 /)
      xa7lum(GEV800,TESLA,0) = 0.57719E+03
      xa7(0:7,GEV800,TESLA,0) = (/ &
          0.52490E+00, 0.42573E+00, 0.69069E+01,-0.62649E+00, &
          0.32380E+00,-0.68958E+00, 0.93819E+01, 0.45671E+00 /)
      xa7lum(TEV1,TESLA,0) = -1.0
      xa7lum(GEV090,JLCNLC,0) = -1.0
      xa7lum(GEV170,JLCNLC,0) = -1.0
      xa7lum(GEV350,JLCNLC,0) = -1.0
      xa7lum(GEV500,JLCNLC,0) = 0.63039E+02
      xa7(0:7,GEV500,JLCNLC,0) = (/ &
          0.58967E+00, 0.34035E+00, 0.63631E+01,-0.63683E+00, &
          0.33383E+00,-0.68803E+00, 0.81005E+01, 0.48702E+00 /)
      xa7lum(TEV1,JLCNLC,0) = 0.12812E+03
      xa7(0:7,TEV1,JLCNLC,0) = (/ &
          0.50222E+00, 0.33773E+00, 0.25681E+01,-0.61711E+00, &
          0.36826E+00,-0.68335E+00, 0.36746E+01, 0.65393E+00 /)
      xa8lum(GEV090,TESLA,1) = -1.0
      xa8lum(GEV170,TESLA,1) = -1.0
      xa8lum(GEV350,TESLA,1) = -1.0
      xa8lum(GEV500,TESLA,1) = -1.0
      xa8lum(GEV800,TESLA,1) = -1.0
      xa8lum(TEV1,  TESLA,1) = -1.0
      xa8lum(GEV090,JLCNLC,1) = -1.0
      xa8lum(GEV170,JLCNLC,1) = -1.0
      xa8lum(GEV350,JLCNLC,1) = -1.0
      xa8lum(GEV500,JLCNLC,1) = 0.239924E+03
      xa8(0:7,GEV500,JLCNLC,1) = (/ &
         0.57025E+00, 0.34004E+00, 0.52864E+01,-0.63405E+00, &
         0.31627E+00,-0.68722E+00, 0.69629E+01, 0.47973E+00 /)
      xa8lum(TEV1,JLCNLC,1) = 0.40858E+03
      xa8(0:7,TEV1,JLCNLC,1) = (/ &
         0.52344E+00, 0.31536E+00, 0.25244E+01,-0.62215E+00, &
         0.31935E+00,-0.68424E+00, 0.35877E+01, 0.57315E+00 /)
      xa8lum(GEV090,TESLA,0) = -1.0
      xa8lum(GEV170,TESLA,0) = -1.0
      xa8lum(GEV350,TESLA,0) = -1.0
      xa8lum(GEV500,TESLA,0) = -1.0
      xa8lum(GEV800,TESLA,0) = -1.0
      xa8lum(TEV1,  TESLA,0) = -1.0
      xa8lum(GEV090,JLCNLC,0) = -1.0
      xa8lum(GEV170,JLCNLC,0) = -1.0
      xa8lum(GEV350,JLCNLC,0) = -1.0
      xa8lum(GEV500,JLCNLC,0) = 0.239924E+03
      xa8(0:7,GEV500,JLCNLC,0) = (/ &
         0.57025E+00, 0.34004E+00, 0.52864E+01,-0.63405E+00, &
         0.31627E+00,-0.68722E+00, 0.69629E+01, 0.47973E+00 /)
      xa8lum(TEV1,JLCNLC,0) = 0.40858E+03
      xa8(0:7,TEV1,JLCNLC,0) = (/ &
         0.52344E+00, 0.31536E+00, 0.25244E+01,-0.62215E+00, &
         0.31935E+00,-0.68424E+00, 0.35877E+01, 0.57315E+00 /)
      xa9lum(GEV090,TESLA,1) = -1.0
      xa9lum(GEV170,TESLA,1) = -1.0
      xa9lum(GEV350,TESLA,1) = -1.0
      xa9lum(GEV500,TESLA,1) = -1.0
      xa9lum(GEV800,TESLA,1) = -1.0
      xa9lum(TEV1,  TESLA,1) = -1.0
      xa9lum(TEV12, TESLA,1) = -1.0
      xa9lum(TEV15, TESLA,1) = -1.0
      xa9lum(TEV16, TESLA,1) = -1.0
      xa9lum(GEV090,JLCNLC,1) = -1.0
      xa9lum(GEV170,JLCNLC,1) = -1.0
      xa9lum(GEV250,JLCNLC,1) = 109.886976
      xa9(0:7,GEV250,JLCNLC,1) = (/ &
          0.65598E+00, 0.34993E+00, 0.13766E+02,-0.64698E+00, &
          0.29984E+00,-0.69053E+00, 0.16444E+02, 0.36060E+00 /)
      xa9lum(GEV350,JLCNLC,1) = -1.0
      xa9lum(GEV500,JLCNLC,1) = 220.806144
      xa9(0:7,GEV500,JLCNLC,1) = (/ &
          0.57022E+00, 0.33782E+00, 0.52811E+01,-0.63540E+00, &
          0.32035E+00,-0.68776E+00, 0.69552E+01, 0.48751E+00 /)
      xa9lum(GEV800,JLCNLC,1) = 304.63488
      xa9(0:7,GEV800,JLCNLC,1) = (/ &
          0.54839E+00, 0.31823E+00, 0.33071E+01,-0.62671E+00, &
          0.31655E+00,-0.68468E+00, 0.45325E+01, 0.53449E+00 /)
      xa9lum(TEV1, JLCNLC,1) = 319.95648
      xa9(0:7,TEV1, JLCNLC,1) = (/ &
          0.56047E+00, 0.29479E+00, 0.28820E+01,-0.62856E+00, & 
          0.29827E+00,-0.68423E+00, 0.39138E+01, 0.52297E+00 /)
      xa9lum(TEV12,JLCNLC,1) = 349.90848
      xa9(0:7,TEV12,JLCNLC,1) = (/ &
          0.56102E+00, 0.28503E+00, 0.24804E+01,-0.62563E+00, &
          0.29002E+00,-0.68376E+00, 0.33854E+01, 0.52736E+00 /)
      xa9lum(TEV15,JLCNLC,1) = 363.15648
      xa9(0:7,TEV15,JLCNLC,1) = (/ &
          0.57644E+00, 0.26570E+00, 0.22007E+01,-0.62566E+00, &
          0.27102E+00,-0.68283E+00, 0.29719E+01, 0.50764E+00 /)
      xa9lum(TEV16,JLCNLC,1) = -1.0
      xa9lum(GEV090,NLCH,1) = -1.0 
      xa9lum(GEV170,NLCH,1) = -1.0 
      xa9lum(GEV250,NLCH,1) = -1.0 
      xa9lum(GEV350,NLCH,1) = -1.0 
      xa9lum(GEV500,NLCH,1) = 371.4624
      xa9(0:7,GEV500,NLCH,1)= (/ &
          0.33933E+00, 0.55165E+00, 0.29138E+01,-0.57341E+00, &
          0.54323E+00,-0.68590E+00, 0.51786E+01, 0.88956E+00 /)
      xa9lum(GEV800,NLCH,1) = -1.0
      xa9lum(TEV1,NLCH,1) = 516.41856
      xa9(0:7,TEV1,NLCH,1)= (/ &
          0.35478E+00, 0.46474E+00, 0.17666E+01,-0.56949E+00, &
          0.49269E+00,-0.68384E+00, 0.31781E+01, 0.91121E+00 /)
      xa9lum(TEV12,NLCH,1) = -1.0
      xa9lum(TEV15,NLCH,1) = 575.06688
      xa9(0:7,TEV15,NLCH,1)= (/ &
          0.38183E+00, 0.40310E+00, 0.13704E+01,-0.57742E+00, &
          0.44548E+00,-0.68341E+00, 0.24956E+01, 0.87448E+00 /)
      xa9lum(TEV16,NLCH,  1) = -1.0
      xa9lum(GEV090,TESLA,0) = -1.0
      xa9lum(GEV170,TESLA,0) = -1.0
      xa9lum(GEV350,TESLA,0) = -1.0
      xa9lum(GEV500,TESLA,0) = -1.0
      xa9lum(GEV800,TESLA,0) = -1.0
      xa9lum(TEV1,  TESLA,0) = -1.0
      xa9lum(TEV12, TESLA,0) = -1.0
      xa9lum(TEV15, TESLA,0) = -1.0
      xa9lum(TEV16, TESLA,0) = -1.0
      xa9lum(GEV090,JLCNLC,0) = -1.0
      xa9lum(GEV170,JLCNLC,0) = -1.0
      xa9lum(GEV250,JLCNLC,0) = 109.886976
      xa9(0:7,GEV250,JLCNLC,0) = (/ &
          0.65598E+00, 0.34993E+00, 0.13766E+02,-0.64698E+00, &
          0.29984E+00,-0.69053E+00, 0.16444E+02, 0.36060E+00 /)
      xa9lum(GEV350,JLCNLC,0) = -1.0
      xa9lum(GEV500,JLCNLC,0) = 220.806144
      xa9(0:7,GEV500,JLCNLC,0) = (/ &
          0.57022E+00, 0.33782E+00, 0.52811E+01,-0.63540E+00, &
          0.32035E+00,-0.68776E+00, 0.69552E+01, 0.48751E+00 /)
      xa9lum(GEV800,JLCNLC,0) = 304.63488
      xa9(0:7,GEV800,JLCNLC,0) = (/ &
          0.54839E+00, 0.31823E+00, 0.33071E+01,-0.62671E+00, &
          0.31655E+00,-0.68468E+00, 0.45325E+01, 0.53449E+00 /)
      xa9lum(TEV1, JLCNLC,0) = 319.95648
      xa9(0:7,TEV1, JLCNLC,0) = (/ &
          0.56047E+00, 0.29479E+00, 0.28820E+01,-0.62856E+00, &
          0.29827E+00,-0.68423E+00, 0.39138E+01, 0.52297E+00 /)
      xa9lum(TEV12,JLCNLC,0) = 349.90848
      xa9(0:7,TEV12,JLCNLC,0) = (/ &
          0.56102E+00, 0.28503E+00, 0.24804E+01,-0.62563E+00, &
          0.29002E+00,-0.68376E+00, 0.33854E+01, 0.52736E+00 /)
      xa9lum(TEV15,JLCNLC,0) = 363.15648
      xa9(0:7,TEV15,JLCNLC,0) = (/ & 
          0.57644E+00, 0.26570E+00, 0.22007E+01,-0.62566E+00, &
          0.27102E+00,-0.68283E+00, 0.29719E+01, 0.50764E+00 /)
      xa9lum(TEV16,JLCNLC,0) = -1.0
      xa9lum(GEV090,NLCH,0) = -1.0
      xa9lum(GEV170,NLCH,0) = -1.0
      xa9lum(GEV250,NLCH,0) = -1.0
      xa9lum(GEV350,NLCH,0) = -1.0
      xa9lum(GEV500,NLCH,0) = 371.4624
      xa9(0:7,GEV500,NLCH,0) = (/ &
          0.33933E+00, 0.55165E+00, 0.29138E+01,-0.57341E+00, &
          0.54323E+00,-0.68590E+00, 0.51786E+01, 0.88956E+00 /)
      xa9lum(GEV800,NLCH,0) = -1.0
      xa9lum(TEV1,NLCH,0)   = 516.41856
      xa9(0:7,TEV1,NLCH,0) = (/ & 
          0.35478E+00, 0.46474E+00, 0.17666E+01,-0.56949E+00, &
          0.49269E+00,-0.68384E+00, 0.31781E+01, 0.91121E+00 /)
      xa9lum(TEV12,NLCH,0) = -1.0
      xa9lum(TEV15,NLCH,0) = 575.06688
      xa9(0:7,TEV15,NLCH,0) = (/ &
          0.38183E+00, 0.40310E+00, 0.13704E+01,-0.57742E+00, &
          0.44548E+00,-0.68341E+00, 0.24956E+01, 0.87448E+00 /)
      xa9lum(TEV16,NLCH,0) = -1.0
      xa10lum = -1
      xa10 = -1
      xa10lum(GEV200,ILC,1) =  56
      xa10(:,GEV200,ILC,1) = (/ &
           0.66253E+00,  0.51646E+00,  0.43632E+02, -0.64508E+00, &
           0.35915E+00, -0.69716E+00,  0.51645E+02,  0.32097E+00 /)
      xa10lum(GEV230,ILC,1) =  83
      xa10(:,GEV230,ILC,1) = (/ &
           0.62360E+00,  0.52780E+00,  0.31915E+02, -0.64171E+00, &
           0.38375E+00, -0.69529E+00,  0.39717E+02,  0.36597E+00 /)
      xa10lum(GEV250,ILC,1) =  97
      xa10(:,GEV250,ILC,1) = (/ &
           0.59996E+00,  0.52141E+00,  0.26647E+02, -0.64331E+00, &
           0.39186E+00, -0.69687E+00,  0.33764E+02,  0.39669E+00 /)
      xa10lum(GEV350,ILC,1) = 100
      xa10(:,GEV350,ILC,1) = (/ &
           0.58875E+00,  0.50027E+00,  0.18594E+02, -0.63380E+00, &
           0.38659E+00, -0.69239E+00,  0.23964E+02,  0.42049E+00 /)
      xa10lum(GEV500,ILC,1) = 180
      xa10(:,GEV500,ILC,1) = (/ &
           0.46755E+00,  0.51768E+00,  0.83463E+01, -0.62311E+00, &
           0.45704E+00, -0.69165E+00,  0.12372E+02,  0.60192E+00 /)
      xa10lum(:,:,0) = xa10lum(:,:,A10NREV)
      xa10(:,:,:,0) = xa10(:,:,:,A10NREV)
      if (circe1_params%magic .ne. 19040616) then
         circe1_params%magic = 19040616
               circe1_params%x1m = 0d0
               circe1_params%x2m = 0d0
               circe1_params%roots = 500D0
               circe1_params%acc = TESLA
               circe1_params%ver = 0
               circe1_params%rev = 0
               circe1_params%chat = 1
               if (xchat .ne. 0) then
                  call circem ('MESSAGE', 'starting up ...')
               endif
      end if
        if ((xchat .ge. 0) .and. (xchat .ne. circe1_params%chat)) then
           circe1_params%chat = xchat
           if (circe1_params%chat .ge. 1) then
              write (msgbuf, 1000) 'chat', circe1_params%chat
   1000       format ('updating `', A, ''' to ', I2)
              call circem ('MESSAGE', msgbuf)
           endif
        else
           if (circe1_params%chat .ge. 2) then
              write (msgbuf, 1100) 'chat', circe1_params%chat
   1100       format ('keeping `', A, ''' at ', I2)
              call circem ('MESSAGE', msgbuf)
           endif
        endif
        if ((xx1m .ge. 0d0) .and. (xx1m .ne. circe1_params%x1m)) then
           circe1_params%x1m = xx1m
           if (circe1_params%chat .ge. 1) then
              write (msgbuf, 1001) 'x1min', circe1_params%x1m
   1001       format ('updating `', A, ''' to ', E12.4)
              call circem ('MESSAGE', msgbuf)
           endif
        else
           if (circe1_params%chat .ge. 2) then
              write (msgbuf, 1101) 'x1min', circe1_params%x1m
   1101       format ('keeping `', A, ''' at ', E12.4)
              call circem ('MESSAGE', msgbuf)
           endif
        endif
        if ((xx2m .ge. 0d0) .and. (xx2m .ne. circe1_params%x2m)) then
           circe1_params%x2m = xx2m
           if (circe1_params%chat .ge. 1) then
              write (msgbuf, 1001) 'x2min', circe1_params%x2m
              call circem ('MESSAGE', msgbuf)
           endif
        else
           if (circe1_params%chat .ge. 2) then
              write (msgbuf, 1101) 'x2min', circe1_params%x2m
              call circem ('MESSAGE', msgbuf)
           endif
        endif
        if ((xroots .ge. 0d0) .and.(xroots .ne. circe1_params%roots)) then
           circe1_params%roots = xroots
           if (circe1_params%chat .ge. 1) then
              write (msgbuf, 1002) 'roots', circe1_params%roots
   1002       format ('updating `', A, ''' to ', F6.1)
              call circem ('MESSAGE', msgbuf)
           endif
        else
           if (circe1_params%chat .ge. 2) then
              write (msgbuf, 1102) 'roots', circe1_params%roots
   1102       format ('keeping `', A, ''' at ', F6.1)
              call circem ('MESSAGE', msgbuf)
           endif
        endif
        if ((xacc .ge. 0) .and.(xacc .ne. circe1_params%acc)) then
           if ((xacc .ge. 1) .and. (xacc .le. NACC)) then
              circe1_params%acc = xacc
              if (circe1_params%chat .ge. 1) then
                 write (msgbuf, 1003) 'acc', accnam(circe1_params%acc)
   1003          format ('updating `', A, ''' to ', A)
                 call circem ('MESSAGE', msgbuf)
              endif
           else
              write (msgbuf, 1203) xacc
   1203       format ('invalid `acc'': ', I8)
              call circem ('ERROR', msgbuf)
              write (msgbuf, 1103) 'acc', accnam(circe1_params%acc)
   1103       format ('keeping `', A, ''' at ', A)
              call circem ('MESSAGE', msgbuf)
           endif
        else
           if (circe1_params%chat .ge. 2) then
              write (msgbuf, 1003) 'acc', accnam(circe1_params%acc)
              call circem ('MESSAGE', msgbuf)
           endif
        endif
        if ((circe1_params%acc .eq. SBNDEE) .or. (circe1_params%acc .eq. TESLEE) &
           .or. (circe1_params%acc .eq. XBNDEE)) then
        call circem ('WARNING', '***********************************')
        call circem ('WARNING', '* The accelerator parameters have *')
        call circem ('WARNING', '* not been endorsed for use in    *')
        call circem ('WARNING', '* an e-e- collider yet!!!         *')
        call circem ('WARNING', '***********************************')
        endif
        if (xver .ge. 0) then
           circe1_params%ver = xver
           if (circe1_params%chat .ge. 1) then
              write (msgbuf, 1000) 'ver', circe1_params%ver
              call circem ('MESSAGE', msgbuf)
           endif
        else
           if (circe1_params%chat .ge. 2) then
              write (msgbuf, 1100) 'ver', circe1_params%ver
              call circem ('MESSAGE', msgbuf)
           endif
        endif
        if ((xrev .ge. 0) .and.(xrev .ne. circe1_params%rev)) then
           circe1_params%rev = xrev
           if (circe1_params%chat .ge. 1) then
              write (msgbuf, 1004) 'rev', circe1_params%rev
   1004       format ('updating `', A, ''' to ', I8)
              call circem ('MESSAGE', msgbuf)
           endif
        else
           if (circe1_params%chat .ge. 2) then
              write (msgbuf, 1104) 'rev', circe1_params%rev
   1104       format ('keeping `', A, ''' at ', I8)
              call circem ('MESSAGE', msgbuf)
           endif
        endif
        ver34 = 0
        if ((circe1_params%ver .eq. 1) .or. (circe1_params%ver .eq. 0)) then
        if (circe1_params%rev .eq. 0) then
           r = 0
        elseif (circe1_params%rev .ge. 19970417) then
           r = 5
        elseif (circe1_params%rev .ge. 19960902) then
           r = 4
        elseif (circe1_params%rev .ge. 19960729) then
           r = 3
        elseif (circe1_params%rev .ge. 19960711) then
           r = 2
        elseif (circe1_params%rev .ge. 19960401) then
           r = 1
        elseif (circe1_params%rev .lt. 19960401) then
           call circem ('ERROR', &
               'no revision of version 1 before 96/04/01 available')
           call circem ('MESSAGE', 'falling back to default')
           r = 1
        endif
        if (circe1_params%chat .ge. 2) then
           write (msgbuf, 2000) circe1_params%rev, r
  2000    format ('mapping date ', I8, ' to revision index ', I2)
           call circem ('MESSAGE', msgbuf)
        endif
        if (circe1_params%roots .eq. 350d0) then
           e = GEV350
        else if ((circe1_params%roots .ge. 340d0) .and. (circe1_params%roots .le. 370d0)) then
           write (msgbuf, 2001) circe1_params%roots, 350d0
           call circem ('MESSAGE', msgbuf)
           e = GEV350
        else if (circe1_params%roots .eq. 500d0) then
           e = GEV500
        else if ((circe1_params%roots .ge. 480d0) .and. (circe1_params%roots .le. 520d0)) then
           write (msgbuf, 2001) circe1_params%roots, 500d0
           call circem ('MESSAGE', msgbuf)
           e = GEV500
        else if (circe1_params%roots .eq. 800d0) then
           e = GEV800
        else if ((circe1_params%roots .ge. 750d0) .and. (circe1_params%roots .le. 850d0)) then
           write (msgbuf, 2001) circe1_params%roots, 800d0
           call circem ('MESSAGE', msgbuf)
           e = GEV800
        else if (circe1_params%roots .eq. 1000d0) then
           e = TEV1
        else if ((circe1_params%roots .ge. 900d0) .and. (circe1_params%roots .le. 1100d0)) then
           write (msgbuf, 2001) circe1_params%roots, 1000d0
           call circem ('MESSAGE', msgbuf)
           e = TEV1
        else if (circe1_params%roots .eq. 1600d0) then
           e = TEV16
        else if ((circe1_params%roots .ge. 1500d0) .and. (circe1_params%roots .le. 1700d0)) then
           write (msgbuf, 2001) circe1_params%roots, 1600d0
           call circem ('MESSAGE', msgbuf)
           e = TEV16
        else
           call circem ('ERROR', &
               'only ROOTS = 350, 500, 800, 1000 and 1600GeV available')
           call circem ('MESSAGE', 'falling back to 500GeV')
           e = GEV500
        endif
        if (xa1lum(e,circe1_params%acc,r) .lt. 0d0) then
           write (msgbuf, 2002) circe1_params%roots, accnam(circe1_params%acc), r
           call circem ('ERROR', msgbuf)
           call circem ('MESSAGE', 'falling back to 500GeV')
           e = GEV500
        end if
        if (circe1_params%chat .ge. 2) then
           if (e .ge. GEV090) then
              write (msgbuf, 2003) circe1_params%roots, e
              call circem ('MESSAGE', msgbuf)
           else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
              write (msgbuf, 2013) circe1_params%roots, elo, ehi
              call circem ('MESSAGE', msgbuf)
           end if
        endif
        circe1_params%lumi = xa1lum (e,circe1_params%acc,r)
        do i = 0, 7
           circe1_params%a1(i) = xa1(i,e,circe1_params%acc,r)
        end do
        else if ((circe1_params%ver .eq. 3) .or. (circe1_params%ver .eq. 4)) then
           ver34 = circe1_params%ver
           circe1_params%ver = 1
        if (circe1_params%rev .eq. 0) then
           r = 0
        elseif (circe1_params%rev .ge. 19970417) then
           r = 5
           if (ver34 .eq. 3) then
              call circem ('WARNING', 'version 3 retired after 97/04/17')
              call circem ('MESSAGE', 'falling back to version 4')
           end if
        else if (circe1_params%rev .ge. 19961022) then
           r = ver34
           if ((circe1_params%roots .ne. 800d0) .or. (circe1_params%acc .ne. TESLA)) then
              call circem ('ERROR', 'versions 3 and 4 before 97/04/17')
              call circem ('ERROR', 'apply to TESLA at 800 GeV only')
              call circem ('MESSAGE', 'falling back to TESLA at 800GeV')
              circe1_params%acc = TESLA
              e = GEV800
           end if
        else if (circe1_params%rev .lt. 19961022) then
           call circem ('ERROR', &
             'no revision of versions 3 and 4 available before 96/10/22')
           call circem ('MESSAGE', 'falling back to default')
           r = 5
        end if
        if (circe1_params%chat .ge. 2) then
           write (msgbuf, 2000) circe1_params%rev, r
           call circem ('MESSAGE', msgbuf)
        endif
        if (circe1_params%roots .eq. 350d0) then
           e = GEV350
        else if ((circe1_params%roots .ge. 340d0) .and. (circe1_params%roots .le. 370d0)) then
           write (msgbuf, 2001) circe1_params%roots, 350d0
           call circem ('MESSAGE', msgbuf)
           e = GEV350
        else if (circe1_params%roots .eq. 500d0) then
           e = GEV500
        else if ((circe1_params%roots .ge. 480d0) .and. (circe1_params%roots .le. 520d0)) then
           write (msgbuf, 2001) circe1_params%roots, 500d0
           call circem ('MESSAGE', msgbuf)
           e = GEV500
        else if (circe1_params%roots .eq. 800d0) then
           e = GEV800
        else if ((circe1_params%roots .ge. 750d0) .and. (circe1_params%roots .le. 850d0)) then
           write (msgbuf, 2001) circe1_params%roots, 800d0
           call circem ('MESSAGE', msgbuf)
           e = GEV800
        else if (circe1_params%roots .eq. 1000d0) then
           e = TEV1
        else if ((circe1_params%roots .ge. 900d0) .and. (circe1_params%roots .le. 1100d0)) then
           write (msgbuf, 2001) circe1_params%roots, 1000d0
           call circem ('MESSAGE', msgbuf)
           e = TEV1
        else if (circe1_params%roots .eq. 1600d0) then
           e = TEV16
        else if ((circe1_params%roots .ge. 1500d0) .and. (circe1_params%roots .le. 1700d0)) then
           write (msgbuf, 2001) circe1_params%roots, 1600d0
           call circem ('MESSAGE', msgbuf)
           e = TEV16
        else
           call circem ('ERROR', &
               'only ROOTS = 350, 500, 800, 1000 and 1600GeV available')
           call circem ('MESSAGE', 'falling back to 500GeV')
           e = GEV500
        endif
        if (xa3lum(e,circe1_params%acc,r) .lt. 0d0) then
           write (msgbuf, 2002) circe1_params%roots, accnam(circe1_params%acc), r
           call circem ('ERROR', msgbuf)
           call circem ('MESSAGE', 'falling back to 500GeV')
           e = GEV500
        endif
        if (circe1_params%chat .ge. 2) then
           if (e .ge. GEV090) then
              write (msgbuf, 2003) circe1_params%roots, e
              call circem ('MESSAGE', msgbuf)
           else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
              write (msgbuf, 2013) circe1_params%roots, elo, ehi
              call circem ('MESSAGE', msgbuf)
           end if
        endif
        circe1_params%lumi = xa3lum (e,circe1_params%acc,r)
        do i = 0, 7
           circe1_params%a1(i) = xa3(i,e,circe1_params%acc,r)
        end do
        else if (circe1_params%ver .eq. 5) then
           circe1_params%ver = 1
        if (circe1_params%rev .eq. 0) then
           r = 0
        elseif (circe1_params%rev .ge. 19980505) then
           r = 1
        elseif (circe1_params%rev .lt. 19980505) then
           call circem ('ERROR', &
            'no revision of version 5 available before 98/05/05')
           call circem ('MESSAGE', 'falling back to default')
           r = 1
        endif
        if (circe1_params%chat .ge. 2) then
           write (msgbuf, 2000) circe1_params%rev, r
           call circem ('MESSAGE', msgbuf)
        endif
        if (circe1_params%acc .ne. TESLA) then
           call circem ('ERROR', 'versions 5 applies to TESLA only')
           circe1_params%acc = TESLA
        end if
        if (circe1_params%roots .eq. 350d0) then
           e = GEV350
        else if ((circe1_params%roots .ge. 340d0) .and. (circe1_params%roots .le. 370d0)) then
           write (msgbuf, 2001) circe1_params%roots, 350d0
           call circem ('MESSAGE', msgbuf)
           e = GEV350
        else if (circe1_params%roots .eq. 500d0) then
           e = GEV500
        else if ((circe1_params%roots .ge. 480d0) .and. (circe1_params%roots .le. 520d0)) then
           write (msgbuf, 2001) circe1_params%roots, 500d0
           call circem ('MESSAGE', msgbuf)
           e = GEV500
        else if (circe1_params%roots .eq. 800d0) then
           e = GEV800
        else if ((circe1_params%roots .ge. 750d0) .and. (circe1_params%roots .le. 850d0)) then
           write (msgbuf, 2001) circe1_params%roots, 800d0
           call circem ('MESSAGE', msgbuf)
           e = GEV800
        else if (circe1_params%roots .eq. 1000d0) then
           e = TEV1
        else if ((circe1_params%roots .ge. 900d0) .and. (circe1_params%roots .le. 1100d0)) then
           write (msgbuf, 2001) circe1_params%roots, 1000d0
           call circem ('MESSAGE', msgbuf)
           e = TEV1
        else if (circe1_params%roots .eq. 1600d0) then
           e = TEV16
        else if ((circe1_params%roots .ge. 1500d0) .and. (circe1_params%roots .le. 1700d0)) then
           write (msgbuf, 2001) circe1_params%roots, 1600d0
           call circem ('MESSAGE', msgbuf)
           e = TEV16
        else
           call circem ('ERROR', &
               'only ROOTS = 350, 500, 800, 1000 and 1600GeV available')
           call circem ('MESSAGE', 'falling back to 500GeV')
           e = GEV500
        endif
        if (xa5lum(e,circe1_params%acc,r) .lt. 0d0) then
           write (msgbuf, 2002) circe1_params%roots, accnam(circe1_params%acc), r
           call circem ('ERROR', msgbuf)
           call circem ('MESSAGE', 'falling back to 500GeV')
           e = GEV500
        endif
        if (circe1_params%chat .ge. 2) then
           if (e .ge. GEV090) then
              write (msgbuf, 2003) circe1_params%roots, e
              call circem ('MESSAGE', msgbuf)
           else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
              write (msgbuf, 2013) circe1_params%roots, elo, ehi
              call circem ('MESSAGE', msgbuf)
           end if
        endif
        circe1_params%lumi = xa5lum (e,circe1_params%acc,r)
        do i = 0, 7
           circe1_params%a1(i) = xa5(i,e,circe1_params%acc,r)
        end do   
        else if (circe1_params%ver .eq. 6) then
           circe1_params%ver = 1
        if (circe1_params%rev .eq. 0) then
           r = 0
        else if (circe1_params%rev .ge. 19990415) then
           r = 1
        else if (circe1_params%rev .lt. 19990415) then
           call circem ('ERROR', &
             'no revision of version 6 available before 1999/04/15')
           call circem ('MESSAGE', 'falling back to default')
           r = 1
        end if
        if (circe1_params%chat .ge. 2) then
           write (msgbuf, 2000) circe1_params%rev, r
           call circem ('MESSAGE', msgbuf)
        endif
        if (circe1_params%acc .ne. TESLA) then
           call circem ('ERROR', 'versions 6 applies to TESLA only')
           circe1_params%acc = TESLA
        end if
        if (circe1_params%roots .eq.  90d0) then
           e = GEV090
        elseif ((circe1_params%roots .ge. 85d0) .and. (circe1_params%roots .le. 95d0)) then
           write (msgbuf, 2001) circe1_params%roots, 90d0
           call circem ('MESSAGE', msgbuf)
           e = GEV090
        elseif (circe1_params%roots .eq. 170d0) then
           e = GEV170
        elseif ((circe1_params%roots .ge. 160d0) .and. (circe1_params%roots .le. 180d0)) then
           write (msgbuf, 2001) circe1_params%roots, 170d0
           call circem ('MESSAGE', msgbuf)
           e = GEV170
        elseif (circe1_params%roots .eq. 350d0) then
           e = GEV350
        elseif ((circe1_params%roots .ge. 340d0) .and. (circe1_params%roots .le. 370d0)) then
           write (msgbuf, 2001) circe1_params%roots, 350d0
           call circem ('MESSAGE', msgbuf)
           e = GEV350
        elseif (circe1_params%roots .eq. 500d0) then
           e = GEV500
        elseif ((circe1_params%roots .ge. 480d0) .and. (circe1_params%roots .le. 520d0)) then
           write (msgbuf, 2001) circe1_params%roots, 500d0
           call circem ('MESSAGE', msgbuf)
           e = GEV500
        else
           call circem ('ERROR', &
               'only ROOTS = 90, 170, 350, and 500GeV available')
           call circem ('MESSAGE', 'falling back to 500GeV')
           e = GEV500
        endif
        if (xa6lum(e,circe1_params%acc,r) .lt. 0d0) then
           write (msgbuf, 2002) circe1_params%roots, accnam(circe1_params%acc), r
           call circem ('ERROR', msgbuf)
           call circem ('MESSAGE', 'falling back to 500GeV')
           e = GEV500
        endif
        if (circe1_params%chat .ge. 2) then
           if (e .ge. GEV090) then
              write (msgbuf, 2003) circe1_params%roots, e
              call circem ('MESSAGE', msgbuf)
           else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
              write (msgbuf, 2013) circe1_params%roots, elo, ehi
              call circem ('MESSAGE', msgbuf)
           end if
        endif
        circe1_params%lumi = xa6lum (e,circe1_params%acc,r)
        do i = 0, 7
           circe1_params%a1(i) = xa6(i,e,circe1_params%acc,r)
        end do
        else if (circe1_params%ver .eq. 7) then
           circe1_params%ver = 1
        if (circe1_params%rev .eq. 0) then
           r = 0
        elseif (circe1_params%rev .ge. 20000426) then
           r = 1
        elseif (circe1_params%rev .lt. 20000426) then
           call circem ('ERROR', &
            'no revision of version 7 available before 2000/04/26')
           call circem ('MESSAGE', 'falling back to default')
           r = 1
        endif
        if (circe1_params%chat .ge. 2) then
           write (msgbuf, 2000) circe1_params%rev, r
           call circem ('MESSAGE', msgbuf)
        endif
        if (circe1_params%acc .ne. TESLA .and. circe1_params%acc .ne. JLCNLC) then
           call circem ('ERROR', &
                       'version 7 applies to TESLA and JLCNLC only')
           call circem ('ERROR', 'falling back to TESLA')
           circe1_params%acc = TESLA
        end if
        e = GEV090 - 1
        elo = e
        ehi = e
        if (circe1_params%acc .eq. TESLA) then
           if (circe1_params%roots .lt.  90d0 - DELTAE) then
              write (msgbuf, 2004) circe1_params%roots, 90d0
              call circem ('MESSAGE', msgbuf)
              e = GEV090
           elseif (abs (circe1_params%roots-090d0) .le. DELTAE) then
              e = GEV090
           elseif (circe1_params%roots .lt. 170d0 - DELTAE) then
              write (msgbuf, 2005) circe1_params%roots, 170d0
              call circem ('MESSAGE', msgbuf)
              e = GEV170
           elseif (abs (circe1_params%roots-170d0) .le. DELTAE) then
              e = GEV170
           elseif (circe1_params%roots .lt. 350d0-DELTAE) then
              write (msgbuf, 2006) circe1_params%roots, 170d0, 350d0
              call circem ('MESSAGE', msgbuf)
              elo = GEV170
              ehi = GEV350
              eloval = 170d0
              ehival = 350d0
           elseif (abs (circe1_params%roots-350d0) .le. DELTAE) then
              e = GEV350
           elseif (circe1_params%roots .lt. 500d0 - DELTAE) then
              write (msgbuf, 2006) circe1_params%roots, 350d0, 500d0
              call circem ('MESSAGE', msgbuf)
              elo = GEV350
              ehi = GEV500
              eloval = 350d0
              ehival = 500d0
           elseif (abs (circe1_params%roots-500d0) .le. DELTAE) then
              e = GEV500
           elseif (circe1_params%roots .lt. 800d0 - DELTAE) then
              write (msgbuf, 2006) circe1_params%roots, 500d0, 800d0
              call circem ('MESSAGE', msgbuf)
              elo = GEV500
              ehi = GEV800
              eloval = 500d0
              ehival = 800d0
           elseif (abs (circe1_params%roots-800d0) .le. DELTAE) then
              e = GEV800
           else
              write (msgbuf, 2005) circe1_params%roots, 800d0
              call circem ('MESSAGE', msgbuf)
              e = GEV800
           endif
        elseif (circe1_params%acc .eq. XBAND) then
           if (circe1_params%roots .lt.  500d0 - DELTAE) then
              write (msgbuf, 2004) circe1_params%roots, 500d0
              call circem ('MESSAGE', msgbuf)
              e = GEV500
           elseif (abs (circe1_params%roots-500d0) .le. DELTAE) then
              e = GEV500
           elseif (circe1_params%roots .lt. 1000d0 - DELTAE) then
              write (msgbuf, 2006) circe1_params%roots, 500d0, 1000d0
              call circem ('MESSAGE', msgbuf)
              elo = GEV500
              ehi = TEV1
              eloval =  500d0
              ehival = 1000d0
           elseif (abs (circe1_params%roots-1000d0) .le. DELTAE) then
              e = TEV1
           else
              write (msgbuf, 2005) circe1_params%roots, 1000d0
              call circem ('MESSAGE', msgbuf)
              e = TEV1
           endif
        endif
        if (circe1_params%chat .ge. 2) then
           if (e .ge. GEV090) then
              write (msgbuf, 2003) circe1_params%roots, e
              call circem ('MESSAGE', msgbuf)
           else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
              write (msgbuf, 2013) circe1_params%roots, elo, ehi
              call circem ('MESSAGE', msgbuf)
           end if
        endif
        if (e .ge. GEV090) then
           circe1_params%lumi = xa7lum(e,circe1_params%acc,r)
           do i = 0, 7
              circe1_params%a1(i) = xa7(i,e,circe1_params%acc,r)
           end do
        else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
           circe1_params%lumi = ((circe1_params%roots-eloval)*xa7lum(ehi,circe1_params%acc,r) &
               + (ehival-circe1_params%roots)*xa7lum(elo,circe1_params%acc,r)) / (ehival - eloval)
           do i = 1, 6
              circe1_params%a1(i) = ((circe1_params%roots-eloval)*xa7(i,ehi,circe1_params%acc,r) &
                  + (ehival-circe1_params%roots)*xa7(i,elo,circe1_params%acc,r)) / (ehival - eloval)
           end do       
           circe1_params%a1(0) = 1d0 - circe1_params%a1(1) * beta(circe1_params%a1(2)+1d0,circe1_params%a1(3)+1d0)
           circe1_params%a1(7) = circe1_params%a1(4) * beta(circe1_params%a1(5)+1d0,circe1_params%a1(6)+1d0)
        endif
        else if (circe1_params%ver .eq. 8) then
           circe1_params%ver = 1
        if (circe1_params%rev .eq. 0) then
           r = 0
        elseif (circe1_params%rev .ge. 20010617) then
           r = 1
        elseif (circe1_params%rev .lt. 20010617) then
           call circem ('ERROR', &
            'no revision of version 8 available before 2001/06/17')
           call circem ('MESSAGE', 'falling back to default')
           r = 1
        endif
        if (circe1_params%chat .ge. 2) then
           write (msgbuf, 2000) circe1_params%rev, r
           call circem ('MESSAGE', msgbuf)
        endif
        if (circe1_params%acc .eq. NLCH) then
           circe1_params%acc = JLCNLC
        end if
        if (circe1_params%acc .ne. JLCNLC) then
           call circem ('ERROR', &
                       'version 8 applies to JLCNLC (NLC H) only')
           call circem ('ERROR', 'falling back to JLCNLC')
           circe1_params%acc = JLCNLC
        end if
        e = GEV090 - 1
        elo = e
        ehi = e
        if (circe1_params%acc .eq. TESLA) then
           if (circe1_params%roots .lt.  90d0 - DELTAE) then
              write (msgbuf, 2004) circe1_params%roots, 90d0
              call circem ('MESSAGE', msgbuf)
              e = GEV090
           elseif (abs (circe1_params%roots-090d0) .le. DELTAE) then
              e = GEV090
           elseif (circe1_params%roots .lt. 170d0 - DELTAE) then
              write (msgbuf, 2005) circe1_params%roots, 170d0
              call circem ('MESSAGE', msgbuf)
              e = GEV170
           elseif (abs (circe1_params%roots-170d0) .le. DELTAE) then
              e = GEV170
           elseif (circe1_params%roots .lt. 350d0-DELTAE) then
              write (msgbuf, 2006) circe1_params%roots, 170d0, 350d0
              call circem ('MESSAGE', msgbuf)
              elo = GEV170
              ehi = GEV350
              eloval = 170d0
              ehival = 350d0
           elseif (abs (circe1_params%roots-350d0) .le. DELTAE) then
              e = GEV350
           elseif (circe1_params%roots .lt. 500d0 - DELTAE) then
              write (msgbuf, 2006) circe1_params%roots, 350d0, 500d0
              call circem ('MESSAGE', msgbuf)
              elo = GEV350
              ehi = GEV500
              eloval = 350d0
              ehival = 500d0
           elseif (abs (circe1_params%roots-500d0) .le. DELTAE) then
              e = GEV500
           elseif (circe1_params%roots .lt. 800d0 - DELTAE) then
              write (msgbuf, 2006) circe1_params%roots, 500d0, 800d0
              call circem ('MESSAGE', msgbuf)
              elo = GEV500
              ehi = GEV800
              eloval = 500d0
              ehival = 800d0
           elseif (abs (circe1_params%roots-800d0) .le. DELTAE) then
              e = GEV800
           else
              write (msgbuf, 2005) circe1_params%roots, 800d0
              call circem ('MESSAGE', msgbuf)
              e = GEV800
           endif
        elseif (circe1_params%acc .eq. XBAND) then
           if (circe1_params%roots .lt.  500d0 - DELTAE) then
              write (msgbuf, 2004) circe1_params%roots, 500d0
              call circem ('MESSAGE', msgbuf)
              e = GEV500
           elseif (abs (circe1_params%roots-500d0) .le. DELTAE) then
              e = GEV500
           elseif (circe1_params%roots .lt. 1000d0 - DELTAE) then
              write (msgbuf, 2006) circe1_params%roots, 500d0, 1000d0
              call circem ('MESSAGE', msgbuf)
              elo = GEV500
              ehi = TEV1
              eloval =  500d0
              ehival = 1000d0
           elseif (abs (circe1_params%roots-1000d0) .le. DELTAE) then
              e = TEV1
           else
              write (msgbuf, 2005) circe1_params%roots, 1000d0
              call circem ('MESSAGE', msgbuf)
              e = TEV1
           endif
        endif
        if (circe1_params%chat .ge. 2) then
           if (e .ge. GEV090) then
              write (msgbuf, 2003) circe1_params%roots, e
              call circem ('MESSAGE', msgbuf)
           else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
              write (msgbuf, 2013) circe1_params%roots, elo, ehi
              call circem ('MESSAGE', msgbuf)
           end if
        endif
        if (e .ge. GEV090) then
           circe1_params%lumi = xa8lum(e,circe1_params%acc,r)
           do i = 0, 7
              circe1_params%a1(i) = xa8(i,e,circe1_params%acc,r)
           end do
        elseif (elo .ge. GEV090 .and. ehi .ge. GEV090) then
           circe1_params%lumi = ((circe1_params%roots-eloval)*xa8lum(ehi,circe1_params%acc,r) &
               + (ehival-circe1_params%roots)*xa8lum(elo,circe1_params%acc,r)) / (ehival - eloval)
           do i = 1, 6
              circe1_params%a1(i) = ((circe1_params%roots-eloval)*xa8(i,ehi,circe1_params%acc,r) &
                  + (ehival-circe1_params%roots)*xa8(i,elo,circe1_params%acc,r)) / (ehival - eloval)
           end do
           circe1_params%a1(0) = 1d0 - circe1_params%a1(1) * beta(circe1_params%a1(2)+1d0,circe1_params%a1(3)+1d0)
           circe1_params%a1(7) = circe1_params%a1(4) * beta(circe1_params%a1(5)+1d0,circe1_params%a1(6)+1d0)
        endif
        else if (circe1_params%ver .eq. 9) then
           circe1_params%ver = 1
        if (circe1_params%rev .eq. 0) then
           r = 0
        elseif (circe1_params%rev .ge. 20020328) then
           r = 1
        elseif (circe1_params%rev .lt. 20020328) then
           call circem ('ERROR', &
            'no revision of version 9 available before 2002/03/28')
           call circem ('MESSAGE', 'falling back to default')
           r = 1
        endif
              if (circe1_params%chat .ge. 2) then
                 write (msgbuf, 2000) circe1_params%rev, r
                 call circem ('MESSAGE', msgbuf)
              endif
        if (circe1_params%acc .ne. JLCNLC .and. circe1_params%acc .ne. NLCH) then
           call circem ('ERROR', &
                       'version 9 applies to JLCNLC and NLCH only')
           call circem ('ERROR', 'falling back to JLCNLC')
           circe1_params%acc = JLCNLC
        end if
        if (circe1_params%acc .eq. JLCNLC) then
                e = GEV090 - 1
                elo = e
                ehi = e
                if (circe1_params%roots .lt. 250d0 - DELTAE) then
                   write (msgbuf, 2004) circe1_params%roots, 250d0
                   call circem ('MESSAGE', msgbuf)
                   e = GEV250
                elseif (abs (circe1_params%roots-250d0) .le. DELTAE) then
                   e = GEV250
                elseif (circe1_params%roots .lt. 500d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 250d0, 500d0
                   call circem ('MESSAGE', msgbuf)
                   elo = GEV250
                   ehi = GEV500
                   eloval = 250d0
                   ehival = 500d0
                elseif (abs (circe1_params%roots-500d0) .le. DELTAE) then
                   e = GEV500
                elseif (circe1_params%roots .lt. 800d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 500d0, 800d0
                   call circem ('MESSAGE', msgbuf)
                   elo = GEV500
                   ehi = GEV800
                   eloval = 500d0
                   ehival = 800d0
                elseif (abs (circe1_params%roots-800d0) .le. DELTAE) then
                   e = GEV800
                elseif (circe1_params%roots .lt. 1000d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 800d0, 1000d0
                   call circem ('MESSAGE', msgbuf)
                   elo = GEV800
                   ehi = TEV1
                   eloval =  800d0
                   ehival = 1000d0
                elseif (abs (circe1_params%roots-1000d0) .le. DELTAE) then
                   e = TEV1
                elseif (circe1_params%roots .lt. 1200d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 1000d0, 1200d0
                   call circem ('MESSAGE', msgbuf)
                   elo = TEV1
                   ehi = TEV12
                   eloval = 1000d0
                   ehival = 1200d0
                elseif (abs (circe1_params%roots-1200d0) .le. DELTAE) then
                   e = TEV12
                elseif (circe1_params%roots .lt. 1500d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 1200d0, 1500d0
                   call circem ('MESSAGE', msgbuf)
                   elo = TEV12
                   ehi = TEV15
                   eloval = 1200d0
                   ehival = 1500d0
                elseif (abs (circe1_params%roots-1500d0) .le. DELTAE) then
                   e = TEV15
                else
                   write (msgbuf, 2005) circe1_params%roots, 1500d0
                   call circem ('MESSAGE', msgbuf)
                   e = TEV15
                endif
        else if (circe1_params%acc .eq. NLCH) then
                e = GEV090 - 1
                elo = e
                ehi = e
                if (circe1_params%roots .lt. 500d0 - DELTAE) then
                   write (msgbuf, 2004) circe1_params%roots, 500d0
                   call circem ('MESSAGE', msgbuf)
                   e = GEV500
                elseif (abs (circe1_params%roots-500d0) .le. DELTAE) then
                   e = GEV500
                elseif (circe1_params%roots .lt. 1000d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 500d0, 1000d0
                   call circem ('MESSAGE', msgbuf)
                   elo = GEV500
                   ehi = TEV1
                   eloval =  500d0
                   ehival = 1000d0
                elseif (abs (circe1_params%roots-1000d0) .le. DELTAE) then
                   e = TEV1
                elseif (circe1_params%roots .lt. 1500d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 1000d0, 1500d0
                   call circem ('MESSAGE', msgbuf)
                   elo = TEV1
                   ehi = TEV15
                   eloval = 1000d0
                   ehival = 1500d0
                elseif (abs (circe1_params%roots-1500d0) .le. DELTAE) then
                   e = TEV15
                else
                   write (msgbuf, 2005) circe1_params%roots, 1500d0
                   call circem ('MESSAGE', msgbuf)
                   e = TEV15
                endif
        end if
              if (circe1_params%chat .ge. 2) then
                 if (e .ge. GEV090) then
                    write (msgbuf, 2003) circe1_params%roots, e
                    call circem ('MESSAGE', msgbuf)
                 else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
                    write (msgbuf, 2013) circe1_params%roots, elo, ehi
                    call circem ('MESSAGE', msgbuf)
                 end if
              endif
        if (e .ge. GEV090) then
           circe1_params%lumi = xa9lum(e,circe1_params%acc,r)
           do i = 0, 7
              circe1_params%a1(i) = xa9(i,e,circe1_params%acc,r)
           end do   
        else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
           circe1_params%lumi = ((circe1_params%roots-eloval)*xa9lum(ehi,circe1_params%acc,r) &
               + (ehival-circe1_params%roots)*xa9lum(elo,circe1_params%acc,r)) / (ehival - eloval)
           do i = 1, 6
              circe1_params%a1(i) = ((circe1_params%roots-eloval)*xa9(i,ehi,circe1_params%acc,r) &
                  + (ehival-circe1_params%roots)*xa9(i,elo,circe1_params%acc,r)) / (ehival - eloval)
           end do       
           circe1_params%a1(0) = 1d0 - circe1_params%a1(1) * beta(circe1_params%a1(2)+1d0,circe1_params%a1(3)+1d0)
           circe1_params%a1(7) = circe1_params%a1(4) * beta(circe1_params%a1(5)+1d0,circe1_params%a1(6)+1d0)
        end if
        else if (circe1_params%ver .eq. 10) then
           circe1_params%ver = 1
        if (circe1_params%rev .eq. 0) then
           r = 0
        elseif (circe1_params%rev .ge. 20140305) then
           r = 1
        elseif (circe1_params%rev .lt. 20140305) then
           call circem ('ERROR', &
            'no revision of version 10 available before 2014/03/05')
           call circem ('MESSAGE', 'falling back to default')
           r = 1
        endif
              if (circe1_params%chat .ge. 2) then
                 write (msgbuf, 2000) circe1_params%rev, r
                 call circem ('MESSAGE', msgbuf)
              endif
        if (circe1_params%acc .ne. ILC) then
           call circem ('ERROR', 'version 10 applies to ILC only')
           call circem ('ERROR', 'falling back to ILC')
           circe1_params%acc = ILC
        end if
        if (circe1_params%acc .eq. ILC) then
                e = -EINVAL
                elo = -EINVAL
                ehi = -EINVAL
                if (circe1_params%roots .lt. 200d0 - DELTAE) then
                   write (msgbuf, 2004) circe1_params%roots, 200d0
                   call circem ('MESSAGE', msgbuf)
                   e = GEV200
                elseif (abs (circe1_params%roots-200d0) .le. DELTAE) then
                   e = GEV200
                elseif (circe1_params%roots .lt. 230d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 200d0, 230d0
                   call circem ('MESSAGE', msgbuf)
                   elo = GEV200
                   ehi = GEV230
                   eloval = 200d0
                   ehival = 230d0
                elseif (abs (circe1_params%roots-230d0) .le. DELTAE) then
                   e = GEV230
                elseif (circe1_params%roots .lt. 250d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 230d0, 250d0
                   call circem ('MESSAGE', msgbuf)
                   elo = GEV230
                   ehi = GEV250
                   eloval = 230d0
                   ehival = 250d0
                elseif (abs (circe1_params%roots-250d0) .le. DELTAE) then
                   e = GEV250
                elseif (circe1_params%roots .lt. 350d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 250d0, 350d0
                   call circem ('MESSAGE', msgbuf)
                   elo = GEV250
                   ehi = GEV350
                   eloval = 250d0
                   ehival = 350d0
                elseif (abs (circe1_params%roots-350d0) .le. DELTAE) then
                   e = GEV350
                elseif (circe1_params%roots .lt. 500d0 - DELTAE) then
                   write (msgbuf, 2006) circe1_params%roots, 350d0, 500d0
                   call circem ('MESSAGE', msgbuf)
                   elo = GEV350
                   ehi = GEV500
                   eloval = 350d0
                   ehival = 500d0
                elseif (abs (circe1_params%roots-500d0) .le. DELTAE) then
                   e = GEV500
                else
                   write (msgbuf, 2005) circe1_params%roots, 500d0
                   call circem ('MESSAGE', msgbuf)
                   e = GEV500
                endif
        end if
              if (circe1_params%chat .ge. 2) then
                 if (e .ge. GEV090) then
                    write (msgbuf, 2003) circe1_params%roots, e
                    call circem ('MESSAGE', msgbuf)
                 else if (elo .ge. GEV090 .and. ehi .ge. GEV090) then
                    write (msgbuf, 2013) circe1_params%roots, elo, ehi
                    call circem ('MESSAGE', msgbuf)
                 end if
              endif
        if (e .ne. EINVAL) then
           circe1_params%lumi = xa10lum(e,circe1_params%acc,r)
           do i = 0, 7
              circe1_params%a1(i) = xa10(i,e,circe1_params%acc,r)
           end do   
        else if (elo .ne. EINVAL .and. ehi .ne. EINVAL) then
           circe1_params%lumi = ((circe1_params%roots-eloval)*xa10lum(ehi,circe1_params%acc,r) &
               + (ehival-circe1_params%roots)*xa10lum(elo,circe1_params%acc,r)) / (ehival - eloval)
           do i = 1, 6
              circe1_params%a1(i) = ((circe1_params%roots-eloval)*xa10(i,ehi,circe1_params%acc,r) &
                  + (ehival-circe1_params%roots)*xa10(i,elo,circe1_params%acc,r)) / (ehival - eloval)
           end do       
           circe1_params%a1(0) = 1d0 - circe1_params%a1(1) * beta(circe1_params%a1(2)+1d0,circe1_params%a1(3)+1d0)
           circe1_params%a1(7) = circe1_params%a1(4) * beta(circe1_params%a1(5)+1d0,circe1_params%a1(6)+1d0)
        end if
        else if (circe1_params%ver .eq. 2) then
        call circem ('PANIC', '*********************************')
        call circem ('PANIC', '* version 2 has been retired,   *')
        call circem ('PANIC', '* please use version 1 instead! *')
        call circem ('PANIC', '*********************************')
        return
        else if (circe1_params%ver .gt. 10) then
           call circem ('PANIC', 'versions >10 not available yet')
           return
        else
           call circem ('PANIC', 'version must be positive')
           return
        end if
        circe1_params%elect0 = circe1_params%a1(0) + circe1_params%a1(1) * KIREPS**(circe1_params%a1(3)+1) / (circe1_params%a1(3)+1)
        circe1_params%elect0 = circe1_params%elect0 / KIREPS
        circe1_params%gamma0 = circe1_params%a1(4) * KIREPS**(circe1_params%a1(5)+1) / (circe1_params%a1(5)+1)
        circe1_params%gamma0 = circe1_params%gamma0 / KIREPS
   2001 format ('treating energy ', F6.1, 'GeV as ',  F6.1, 'GeV')
  2002 format ('energy ', F6.1, ' not available for ', A6,' in revison ', I2)
  2003 format ('mapping energy ', F6.1, ' to energy index ', I2)
  2013 format ('mapping energy ', F6.1, ' to energy indices ', I2, ' and ', I2)
  2004 format ('energy ', F6.1, 'GeV too low, using spectrum for ', F6.1, 'GeV')
  2005 format ('energy ', F6.1, 'GeV too high, using spectrum for ', F6.1, 'GeV')
  2006 format ('energy ', F6.1, 'GeV interpolated between ', F6.1, ' and ', F6.1, 'GeV')
    end subroutine circes
    subroutine circex (xx1m, xx2m, xroots, cacc, xver, xrev, xchat)
      real(kind=double) :: xx1m, xx2m, xroots   
      character(*) :: cacc
      integer :: xver, xrev, xchat
      integer :: xacc, i
    integer, parameter :: SBAND  = 1
    integer, parameter :: TESLA  = 2
    integer, parameter :: XBAND  = 3
    integer, parameter :: JLCNLC = 3
    integer, parameter :: SBNDEE = 4
    integer, parameter :: TESLEE = 5
    integer, parameter :: XBNDEE = 6
    integer, parameter :: NLCH   = 7
    integer, parameter :: ILC    = 8
    integer, parameter :: CLIC   = 9
    integer, parameter :: NACC   = 9
        character(len=6), dimension(NACC) :: accnam
        data accnam(SBAND)  /'SBAND'/
        data accnam(TESLA)  /'TESLA'/
        data accnam(JLCNLC) /'JLCNLC'/
        data accnam(SBNDEE) /'SBNDEE'/
        data accnam(TESLEE) /'TESLEE'/
        data accnam(XBNDEE) /'XBNDEE'/
        data accnam(NLCH) /'NLC H'/
        data accnam(ILC) /'ILC'/
        data accnam(CLIC) /'CLIC'/
      xacc = -1
      do i = 1, NACC
         if (trim (accnam(i)) == trim (cacc)) then
           xacc = i
         end if
      end do
      call circes (xx1m, xx2m, xroots, xacc, xver, xrev, xchat)
    end subroutine circex
    subroutine circel (l)
      real(kind=double), intent(out) :: l    
      l = circe1_params%lumi
    end subroutine circel

    function circee (x1, x2)
      real(kind=double) :: x1, x2
      real(kind=double) :: circee
      real(kind=double) :: d1, d2
        if (circe1_params%magic .ne. MAGIC0) then
           call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
        endif
      circee = -1.0
      if ((circe1_params%ver .eq. 1) .or. (circe1_params%ver .eq. 0)) then
        if (x1 .eq. 1d0) then
           d1 = circe1_params%a1(0)
        elseif (x1 .lt. 1d0 .and. x1 .gt. 0d0) then
           d1 = circe1_params%a1(1) * x1**circe1_params%a1(2) * (1d0 - x1)**circe1_params%a1(3)
        elseif (x1 .eq. -1d0) then
           d1 = 1d0 - circe1_params%a1(0)
        else
           d1 = 0d0
        endif
        if (x2 .eq. 1d0) then
           d2 = circe1_params%a1(0)
        elseif (x2 .lt. 1d0 .and. x2 .gt. 0d0) then
           d2 = circe1_params%a1(1) * x2**circe1_params%a1(2) * (1d0 - x2)**circe1_params%a1(3)
        elseif (x2 .eq. -1d0) then
           d2 = 1d0 - circe1_params%a1(0)
        else
           d2 = 0d0
        endif
        circee = d1 * d2
        else if (circe1_params%ver .eq. 2) then
        call circem ('PANIC', '*********************************')
        call circem ('PANIC', '* version 2 has been retired,   *')
        call circem ('PANIC', '* please use version 1 instead! *')
        call circem ('PANIC', '*********************************')
        return
        else if (circe1_params%ver .gt. 10) then
           call circem ('PANIC', 'versions >10 not available yet')
           return
        else
           call circem ('PANIC', 'version must be positive')
           return
        end if
    end function circee

    function circeg (x1, x2)
      real(kind=double) :: x1, x2
      real(kind=double) :: circeg
      real(kind=double) :: d1, d2
        if (circe1_params%magic .ne. MAGIC0) then
           call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
        endif
      circeg = -1.0
      if ((circe1_params%ver .eq. 1) .or. (circe1_params%ver .eq. 0)) then
      if (x1 .eq. 1d0) then
         d1 = circe1_params%a1(0)
      else if (x1 .lt. 1d0 .and. x1 .gt. 0d0) then
         d1 = circe1_params%a1(1) * x1**circe1_params%a1(2) * (1d0 - x1)**circe1_params%a1(3)
      else if (x1 .eq. -1d0) then
         d1 = 1d0 - circe1_params%a1(0)
      else
         d1 = 0d0
      end if
      if (x2 .lt. 1d0 .and. x2 .gt. 0d0) then
         d2 = circe1_params%a1(4) * x2**circe1_params%a1(5) * (1d0 - x2)**circe1_params%a1(6)
      else if (x2 .eq. -1d0) then
         d2 = circe1_params%a1(7)
      else
         d2 = 0d0
      end if
      circeg = d1 * d2
        else if (circe1_params%ver .eq. 2) then
        call circem ('PANIC', '*********************************')
        call circem ('PANIC', '* version 2 has been retired,   *')
        call circem ('PANIC', '* please use version 1 instead! *')
        call circem ('PANIC', '*********************************')
        return
        else if (circe1_params%ver .gt. 10) then
           call circem ('PANIC', 'versions >10 not available yet')
           return
        else
           call circem ('PANIC', 'version must be positive')
           return
        end if
    end function circeg

    function circgg (x1, x2)
      real(kind=double) :: x1, x2
      real(kind=double) :: circgg
      real(kind=double) :: d1, d2
        if (circe1_params%magic .ne. MAGIC0) then
           call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
        endif
      circgg = -1.0
      if ((circe1_params%ver .eq. 1) .or. (circe1_params%ver .eq. 0)) then
        if (x1 .lt. 1d0 .and. x1 .gt. 0d0) then
           d1 = circe1_params%a1(4) * x1**circe1_params%a1(5) * (1d0 - x1)**circe1_params%a1(6)
        elseif (x1 .eq. -1d0) then
           d1 = circe1_params%a1(7)
        else
           d1 = 0d0
        endif
        if (x2 .lt. 1d0 .and. x2 .gt. 0d0) then
           d2 = circe1_params%a1(4) * x2**circe1_params%a1(5) * (1d0 - x2)**circe1_params%a1(6)
        elseif (x2 .eq. -1d0) then
           d2 = circe1_params%a1(7)
        else
           d2 = 0d0
        endif
        circgg = d1 * d2
        else if (circe1_params%ver .eq. 2) then
        call circem ('PANIC', '*********************************')
        call circem ('PANIC', '* version 2 has been retired,   *')
        call circem ('PANIC', '* please use version 1 instead! *')
        call circem ('PANIC', '*********************************')
        return
        else if (circe1_params%ver .gt. 10) then
           call circem ('PANIC', 'versions >10 not available yet')
           return
        else
           call circem ('PANIC', 'version must be positive')
           return
        end if
    end function circgg

    function beta (a, b)      
      real(kind=double) :: a, b, beta
      beta = exp (dlogam(a) + dlogam(b) - dlogam(a+b))
    end function beta

  !!! CERNLIB C304

    function dlogam (x)
      real(kind=double) :: dlogam
      real(kind=double), dimension(7) :: p1, q1, p2, q2, p3, q3
      real(kind=double), dimension(5) :: c, xl
      real(kind=double) :: x, y, zero, one, two, half, ap, aq
      integer :: i
      data ZERO /0.0D0/, ONE /1.0D0/, TWO /2.0D0/, HALF /0.5D0/
      data XL /0.0D0,0.5D0,1.5D0,4.0D0,12.0D0/
      data p1 /+3.8428736567460D+0, +5.2706893753010D+1, &
               +5.5584045723515D+1, -2.1513513573726D+2, &
               -2.4587261722292D+2, -5.7500893603041D+1, &
               -2.3359098949513D+0/
      data q1 /+1.0000000000000D+0, +3.3733047907071D+1, &
               +1.9387784034377D+2, +3.0882954973424D+2, &
               +1.5006839064891D+2, +2.0106851344334D+1, &
               +4.5717420282503D-1/
      data p2 /+4.8740201396839D+0, +2.4884525168574D+2, &
               +2.1797366058896D+3, +3.7975124011525D+3, &
               -1.9778070769842D+3, -3.6929834005591D+3, &
               -5.6017773537804D+2/
      data q2 /+1.0000000000000D+0, +9.5099917418209D+1, &
               +1.5612045277929D+3, +7.2340087928948D+3, &
               +1.0459576594059D+4, +4.1699415153200D+3, &
               +2.7678583623804D+2/
      data p3 /-6.8806240094594D+3, -4.3069969819571D+5, &
               -4.7504594653440D+6, -2.9423445930322D+6, &
               +3.6321804931543D+7, -3.3567782814546D+6, &
               -2.4804369488286D+7/
      data q3 /+1.0000000000000D+0, -1.4216829839651D+3, &
               -1.5552890280854D+5, -3.4152517108011D+6, &
               -2.0969623255804D+7, -3.4544175093344D+7, &
               -9.1605582863713D+6/
      data c / 1.1224921356561D-1,  7.9591692961204D-2, &
              -1.7087794611020D-3,  9.1893853320467D-1, &
               1.3469905627879D+0/
      if (x .le. xl(1)) then
           print *, 'ERROR: DLOGAM non positive argument: ', X
           dlogam = zero
      end if
      if (x .le. xl(2)) then
        y = x + one
        ap = p1(1)
        aq = q1(1)
        do i = 2, 7
          ap = p1(i) + y * ap
          aq = q1(i) + y * aq
        end do
        y = - log(x) + x * ap / aq
      else if (x .le. xl(3)) then
        ap = p1(1)
        aq = q1(1)
        do i = 2, 7
           ap = p1(i) + x * ap
           aq = q1(i) + x * aq
        end do
        y = (x - one) * ap / aq
      else if (x .le. xl(4)) then
        ap = p2(1)
        aq = q2(1)
        do i = 2, 7
           ap = p2(i) + x * ap
           aq = q2(i) + x * aq
        end do     
        y = (x-two) * ap / aq
      else if (x .le. xl(5)) then
        ap = p3(1)
        aq = q3(1)
        do i = 2, 7
           ap = p3(i) + x * ap
           aq = q3(i) + x * aq
        end do 
        y = ap / aq
      else
       y = one / x**2
       y = (x-half) * log(x) - x + c(4) + &
           (c(1) + y * (c(2) + y * c(3))) / ((c(5) + y) * x)
      end if
      dlogam = y
    end function dlogam

    function kirke (x1, x2, p1, p2)
      real(kind=double) :: x1, x2
      real(kind=double) :: kirke
      integer :: p1, p2
            if (circe1_params%magic .ne. MAGIC0) then
               call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
            endif
      kirke = -1.0
      if (abs(p1) .eq. C1_ELECTRON) then
         if (abs(p2) .eq. C1_ELECTRON) then
            kirke = kirkee (x1, x2)
         else if (p2 .eq. C1_PHOTON) then
            kirke = kirkeg (x1, x2)
         end if
      else if (p1 .eq. C1_PHOTON) then
         if (abs(p2) .eq. C1_ELECTRON) then
            kirke = kirkeg (x2, x1)
         else if (p2 .eq. C1_PHOTON) then
            kirke = kirkgg (x1, x2)
         end if
      endif
    end function kirke

    function kirkee (x1, x2)
      real(kind=double) :: x1, x2
      real(kind=double) :: kirkee
      real(kind=double) :: d1, d2
            if (circe1_params%magic .ne. MAGIC0) then
               call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
            endif
      kirkee = -1.0
      if ((circe1_params%ver .eq. 1) .or. (circe1_params%ver .eq. 0)) then
               if (x1 .gt. 1d0) then
                  d1 = 0d0
               elseif (x1 .ge. (1d0 - KIREPS)) then
                  d1 = circe1_params%elect0
               elseif (x1 .ge. 0d0) then
                  d1 = circe1_params%a1(1) * x1**circe1_params%a1(2) * (1d0 - x1)**circe1_params%a1(3)
               else
                  d1 = 0d0
               endif
               if (x2 .gt. 1d0) then
                  d2 = 0d0
               elseif (x2 .ge. (1d0 - KIREPS)) then
                  d2 = circe1_params%elect0
               elseif (x2 .ge. 0d0) then
                  d2 = circe1_params%a1(1) * x2**circe1_params%a1(2) * (1d0 - x2)**circe1_params%a1(3)
               else
                  d2 = 0d0
               endif
               kirkee = d1 * d2         
            else if (circe1_params%ver .eq. 2) then
            call circem ('PANIC', '*********************************')
            call circem ('PANIC', '* version 2 has been retired,   *')
            call circem ('PANIC', '* please use version 1 instead! *')
            call circem ('PANIC', '*********************************')
            return
            else if (circe1_params%ver .gt. 10) then
               call circem ('PANIC', 'versions >10 not available yet')
               return
            else
               call circem ('PANIC', 'version must be positive')
               return
            end if
    end function kirkee

    function kirkeg (x1, x2)
      real(kind=double) :: x1, x2
      real(kind=double) :: kirkeg
      real(kind=double) :: d1, d2
            if (circe1_params%magic .ne. MAGIC0) then
               call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
            endif
      kirkeg = -1.0
      if ((circe1_params%ver .eq. 1) .or. (circe1_params%ver .eq. 0)) then
               if (x1 .gt. 1d0) then
                  d1 = 0d0
               elseif (x1 .ge. (1d0 - KIREPS)) then
                  d1 = circe1_params%elect0
               elseif (x1 .ge. 0d0) then
                  d1 = circe1_params%a1(1) * x1**circe1_params%a1(2) * (1d0 - x1)**circe1_params%a1(3)
               else
                  d1 = 0d0
               endif
               if (x2 .gt. 1d0) then
                  d2 = 0d0
               elseif (x2 .gt. KIREPS) then
                  d2 = circe1_params%a1(4) * x2**circe1_params%a1(5) * (1d0 - x2)**circe1_params%a1(6)
               elseif (x2 .ge. 0d0) then
                  d2 = circe1_params%gamma0
               else
                  d2 = 0d0
               endif
               kirkeg = d1 * d2         
            else if (circe1_params%ver .eq. 2) then
            call circem ('PANIC', '*********************************')
            call circem ('PANIC', '* version 2 has been retired,   *')
            call circem ('PANIC', '* please use version 1 instead! *')
            call circem ('PANIC', '*********************************')
            return
            else if (circe1_params%ver .gt. 10) then
               call circem ('PANIC', 'versions >10 not available yet')
               return
            else
               call circem ('PANIC', 'version must be positive')
               return
            end if
    end function kirkeg

    function kirkgg (x1, x2)
      real(kind=double) :: x1, x2
      real(kind=double) :: kirkgg
      real(kind=double) :: d1, d2
            if (circe1_params%magic .ne. MAGIC0) then
               call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
            endif
      kirkgg = -1.0
      if ((circe1_params%ver .eq. 1) .or. (circe1_params%ver .eq. 0)) then
               if (x1 .gt. 1d0) then
                  d1 = 0d0
               elseif (x1 .gt. KIREPS) then
                  d1 = circe1_params%a1(4) * x1**circe1_params%a1(5) * (1d0 - x1)**circe1_params%a1(6)
               elseif (x1 .ge. 0d0) then
                  d1 = circe1_params%gamma0
               else
                  d1 = 0d0
               endif
               if (x2 .gt. 1d0) then
                  d2 = 0d0
               elseif (x2 .gt. KIREPS) then
                  d2 = circe1_params%a1(4) * x2**circe1_params%a1(5) * (1d0 - x2)**circe1_params%a1(6)
               elseif (x2 .ge. 0d0) then
                  d2 = circe1_params%gamma0
               else
                  d2 = 0d0
               endif
               kirkgg = d1 * d2         
            else if (circe1_params%ver .eq. 2) then
            call circem ('PANIC', '*********************************')
            call circem ('PANIC', '* version 2 has been retired,   *')
            call circem ('PANIC', '* please use version 1 instead! *')
            call circem ('PANIC', '*********************************')
            return
            else if (circe1_params%ver .gt. 10) then
               call circem ('PANIC', 'versions >10 not available yet')
               return
            else
               call circem ('PANIC', 'version must be positive')
               return
            end if
    end function kirkgg

    subroutine rng_call (u, rng, rng_obj)
      real(kind=double), intent(out) :: u
      procedure(rng_proc), optional :: rng
      class(rng_type), intent(inout), optional :: rng_obj
      if (present (rng)) then
         call rng (u)
      else if (present (rng_obj)) then
         call rng_obj%generate (u)
      else
         call circem ('PANIC', &
              'generator requires either rng or rng_obj argument')
      end if
    end subroutine rng_call
    
    subroutine girce (x1, x2, p1, p2, rng, rng_obj)
      real(kind=double), intent(out) :: x1, x2
      integer :: p1, p2
      procedure(rng_proc), optional :: rng
      class(rng_type), intent(inout), optional :: rng_obj
      real(kind=double) :: u, w
        if (circe1_params%magic .ne. MAGIC0) then
           call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
        endif
      do
        w = 1d0 / (1d0 + circgg (-1d0, -1d0))
        call rng_call (u, rng, rng_obj)
        if (u*u .le. w) then
           p1 = C1_POSITRON
        else
           p1 = C1_PHOTON
        end if
        call rng_call (u, rng, rng_obj)
        if (u*u .le. w) then
           p2 = C1_ELECTRON
        else
           p2 = C1_PHOTON
        end if
        if (abs(p1) .eq. C1_ELECTRON) then
           if (abs(p2) .eq. C1_ELECTRON) then
              call gircee (x1, x2, rng, rng_obj)
           else if (p2 .eq. C1_PHOTON) then
              call girceg (x1, x2, rng, rng_obj)
           end if
        else if (p1 .eq. C1_PHOTON) then
           if (abs(p2) .eq. C1_ELECTRON) then
              call girceg (x2, x1, rng, rng_obj)
           else if (p2 .eq. C1_PHOTON) then
              call gircgg (x1, x2, rng, rng_obj)
           end if
        end if
        if ((x1 .ge. circe1_params%x1m) .and. (x2 .ge. circe1_params%x2m)) exit
      end do   
    end subroutine girce

    subroutine gircee (x1, x2, rng, rng_obj)
      real(kind=double), intent(out) :: x1, x2
      procedure(rng_proc), optional :: rng
      class(rng_type), intent(inout), optional :: rng_obj
      real(kind=double) :: u
            if (circe1_params%magic .ne. MAGIC0) then
               call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
            endif
      x1 = 1
      x2 = 1
      if ((circe1_params%ver .eq. 1) .or. (circe1_params%ver .eq. 0)) then
               call rng_call (u, rng, rng_obj)
               if (u .le. circe1_params%a1(0)) then
                  x1 = 1d0
               else
                  x1 = 1d0 - girceb (0d0, 1d0-circe1_params%x1m, &
                                     circe1_params%a1(3)+1d0, circe1_params%a1(2)+1d0, &
                                     rng, rng_obj)
               endif
               call rng_call (u, rng, rng_obj)
               if (u .le. circe1_params%a1(0)) then
                  x2 = 1d0
               else
                  x2 = 1d0 - girceb (0d0, 1d0-circe1_params%x2m, &
                                     circe1_params%a1(3)+1d0, circe1_params%a1(2)+1d0, &
                                     rng, rng_obj)
               endif
            else if (circe1_params%ver .eq. 2) then
            call circem ('PANIC', '*********************************')
            call circem ('PANIC', '* version 2 has been retired,   *')
            call circem ('PANIC', '* please use version 1 instead! *')
            call circem ('PANIC', '*********************************')
            return
            else if (circe1_params%ver .gt. 10) then
               call circem ('PANIC', 'versions >10 not available yet')
               return
            else
               call circem ('PANIC', 'version must be positive')
               return
            end if
    end subroutine gircee

    subroutine girceg (x1, x2, rng, rng_obj)
      real(kind=double), intent(out) :: x1, x2
      procedure(rng_proc), optional :: rng
      class(rng_type), intent(inout), optional :: rng_obj
      real(kind=double) :: u
            if (circe1_params%magic .ne. MAGIC0) then
               call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
            endif
      x1 = 1
      x2 = 1
      if ((circe1_params%ver .eq. 1) .or. (circe1_params%ver .eq. 0)) then
               call rng_call (u, rng, rng_obj)
               if (u .le. circe1_params%a1(0)) then
                  x1 = 1d0
               else
                  x1 = 1d0 - girceb (0d0, 1d0-circe1_params%x1m, &
                                     circe1_params%a1(3)+1d0, circe1_params%a1(2)+1d0, &
                                     rng, rng_obj)
               endif
               x2 = girceb (circe1_params%x2m, 1d0, &
                            circe1_params%a1(5)+1d0, circe1_params%a1(6)+1d0, &
                            rng, rng_obj)
            else if (circe1_params%ver .eq. 2) then
            call circem ('PANIC', '*********************************')
            call circem ('PANIC', '* version 2 has been retired,   *')
            call circem ('PANIC', '* please use version 1 instead! *')
            call circem ('PANIC', '*********************************')
            return
            else if (circe1_params%ver .gt. 10) then
               call circem ('PANIC', 'versions >10 not available yet')
               return
            else
               call circem ('PANIC', 'version must be positive')
               return
            end if
    end subroutine girceg

    subroutine gircgg (x1, x2, rng, rng_obj)
      real(kind=double), intent(out) :: x1, x2
      procedure(rng_proc), optional :: rng
      class(rng_type), intent(inout), optional :: rng_obj
            if (circe1_params%magic .ne. MAGIC0) then
               call circes (-1d0, -1d0, -1d0, -1, -1, -1, -1)
            endif
      x1 = 1
      x2 = 1
      if ((circe1_params%ver .eq. 1) .or. (circe1_params%ver .eq. 0)) then
               x1 = girceb (circe1_params%x1m, 1d0, &
                            circe1_params%a1(5)+1d0, circe1_params%a1(6)+1d0, &
                            rng, rng_obj)
               x2 = girceb (circe1_params%x2m, 1d0, &
                            circe1_params%a1(5)+1d0, circe1_params%a1(6)+1d0, &
                            rng, rng_obj)
            else if (circe1_params%ver .eq. 2) then
            call circem ('PANIC', '*********************************')
            call circem ('PANIC', '* version 2 has been retired,   *')
            call circem ('PANIC', '* please use version 1 instead! *')
            call circem ('PANIC', '*********************************')
            return
            else if (circe1_params%ver .gt. 10) then
               call circem ('PANIC', 'versions >10 not available yet')
               return
            else
               call circem ('PANIC', 'version must be positive')
               return
            end if
    end subroutine gircgg

    function girceb (xmin, xmax, a, b, rng, rng_obj)
      real(kind=double) :: xmin, xmax, a, b
      real(kind=double) :: girceb
      procedure(rng_proc), optional :: rng
      class(rng_type), intent(inout), optional :: rng_obj
      real(kind=double) :: t, p, u, umin, umax, x, w
            if ((a .ge. 1d0) .or. (b .le. 1d0)) then
               girceb = -1d0
               call circem ('ERROR', 'beta-distribution expects a<1<b')
               return
            end if
                  t = (1d0 - a) / (b + 1d0 - a)
            p = b*t / (b*t + a * (1d0 - t)**b)
            if (xmin .le. 0d0) then
               umin = 0d0
            elseif (xmin .lt. t) then
               umin = p * (xmin/t)**a
            elseif (xmin .eq. t) then
               umin = p
            elseif (xmin .lt. 1d0) then
               umin = 1d0 - (1d0 - p) * ((1d0 - xmin)/(1d0 - t))**b
            else
               umin = 1d0
            endif
            if (xmax .ge. 1d0) then
               umax = 1d0
            elseif (xmax .gt. t) then
               umax = 1d0 - (1d0 - p) * ((1d0 - xmax)/(1d0 - t))**b
            elseif (xmax .eq. t) then
               umax = p
            elseif (xmax .gt. 0d0) then
               umax = p * (xmax/t)**a
            else
               umax = 0d0
            endif
            if (umax .lt. umin) then
               girceb = -1d0
               return
            endif
      do 
               call rng_call (u, rng, rng_obj)
               u = umin + (umax - umin) * u
               if (u .le. p) then
                  x = t * (u/p)**(1d0/a)
                  w = (1d0 - x)**(b-1d0)
               else
                  x = 1d0 - (1d0 - t) * ((1d0 - u)/(1d0 - p))**(1d0/b)
                  w = (x/t)**(a-1d0)
               end if
         call rng_call (u, rng, rng_obj)
         if (w .gt. u) exit
      end do 
      girceb = x
    end function girceb

    subroutine circem (errlvl, errmsg)
      character(len=*) :: errlvl, errmsg
      integer, save :: errcnt = 0
      if (errlvl .eq. 'MESSAGE') then
         print *, 'circe1:message: ', errmsg
      else if (errlvl .eq. 'WARNING') then
         if (errcnt .lt. 100) then
            errcnt = errcnt + 1
            print *, 'circe1:warning: ', errmsg
         else if (errcnt .eq. 100) then
            errcnt = errcnt + 1
            print *, 'circe1:message: more than 100 messages'
            print *, 'circe1:message: turning warnings off'
         end if
      else if (errlvl .eq. 'ERROR') then
         if (errcnt .lt. 200) then
            errcnt = errcnt + 1
            print *, 'circe1:error:   ', errmsg
         else if (errcnt .eq. 200) then
            errcnt = errcnt + 1
            print *, 'circe1:message: more than 200 messages'
            print *, 'circe1:message: turning error messages off'
         endif
      else if (errlvl .eq. 'PANIC') then
         if (errcnt .lt. 300) then
            errcnt = errcnt + 1
            print *, 'circe1:panic:   ', errmsg
         else if (errcnt .eq. 300) then
            errcnt = errcnt + 1
            print *, 'circe1:message: more than 300 messages'
            print *, 'circe1:message: turning panic messages off'
         end if
      else
         print *, 'circe1:panic:    invalid error code ', errlvl
      end if
    end subroutine circem

end module circe1
