module const_proc
use hycom_dimen
implicit none

      integer :: indoi, indoj1, indoj2		&! Indo Throughflow
        	,idrk1, jdrk1, idrk2, jdrk2	&! Drake Passage
       		,iberi, jberi			&! Bering Strait
      	        ,ikuro1, ikuro2, jkuro		&! Kuroshio
             	,igulf1, igulf2, jgulf		&! Gulf Stream
             	,imed, jmed		 	 ! Med outflow

      real, parameter :: spval = -0.03125, onem=9806., onecm=98.06	&
          ,grvty=9.806, rho =1000., dz=1., spcifh=4185, flag=-999.	&
          ,epsil_p=0.01                             ! 1cm cutoff
      integer, parameter :: iorign=1, jorign=1,julian=365
      logical, parameter :: diag=.true., rhodot=.false.
      integer, dimension(12),parameter :: 				&
!        jdofm=(/30,60,90,120,150,180,210,240,270,300,330,360/)		&
         jdofm=(/31,59,90,120,151,181,212,243,273,304,334,365/)
      integer, dimension(12), parameter :: 				&
         monlg=(/31,28,31,30,31,30,31,31,30,31,30,31/)
      integer, parameter :: mo1=1,mo2=12
      character*3,dimension(13),parameter:: 				&
        amon=(/'JAN','FEB','MAR','APR','MAY','JUN',			&
               'JUL','AUG','SEP','OCT','NOV','DEC','ANN'/)
      real, dimension(33), parameter ::					&
       z33=(/0.,  10.,  20.,  30.,  50.,  75., 100., 125., 150., 200.,	&
           250., 300., 400., 500., 600., 700., 800., 900.,1000.,1100.,	&
          1200.,1300.,1400.,1500.,1750.,2000.,2500.,3000.,3500.,4000.,	&
          4500.,5000.,5500./)

!     integer, parameter:: lp=6,itest=301,jtest=181,iatest=181,jatest=31	! for 387
!     integer, parameter:: lp=6,itest=229,jtest=181,iatest=181,jatest=91
      integer, parameter:: lp=6,itest=300,jtest=150,iatest=181,jatest=91

      character(len=128) :: &
        path0='/discover/nobackup/projects/giss/prod_input_files/', & ! directory input files for modelE
        path1='./',     &! directory hycom input files (.out)
        path2='./'       ! directory output diagnostic files
      character(len=80) :: hycomtopo, latlonij, basinmask, flnmcoso, flnmo2a
      character(len=20) :: runid="Exxx"
      integer           :: ny1=1800, ny2=1800
      logical           :: monave_convert,solo_convert,timav,cnvert
!
contains

  subroutine const
      if (idm==387) then
        indoi=251; indoj1=117; indoj2=134            ! Indo Throughflow: 2 land points
        idrk1=323; jdrk1=291; idrk2=342; jdrk2=300   ! Drake Passage: 2 land points
        iberi=138; jberi=189                         ! Bering Strait: single point in the channel
        ikuro1=197; ikuro2=198; jkuro=129            ! Kuroshio
        igulf1=199; igulf2=200; jgulf=280            ! Gulf Stream
        imed=190; jmed=356                           ! Med outflow: single point in the channel
        hycomtopo='depth387x360_2009.4big'           ! same topo as in rundeck
        latlonij='latlon387x360.4bin'
        basinmask='ibasin387x360_topo2009_col3_oct2015.txt' ! same mask as in rundeck
        flnmcoso='cososino387x360.8bin'
        flnmo2a='wgt_o2a_h387x360topo2009_360x181_may2016.8bin'
      else if (idm==359) then
        indoi=233; indoj1=106; indoj2=135            ! Indo Throughflow: 2 land points
        idrk1=295; jdrk1=294; idrk2=312; jdrk2=302   ! Drake Passage: 2 land points
        iberi=138; jberi=189                         ! Bering Strait: single point in the channel
        ikuro1=197; ikuro2=198; jkuro=129            ! Kuroshio
        igulf1=199; igulf2=200; jgulf=280            ! Gulf Stream
        imed=191; jmed=355                           ! Med outflow: single point in the channel
        hycomtopo='depth359x360_oct2015.4big'        ! same topo as in rundeck
        latlonij='latlon359x360.4bin'
        basinmask='ibasin359x360_oct2015_col3.txt'   ! same mask as in rundeck
        flnmcoso='cososino359x360.8bin'
        flnmo2a='wgt_o2a_h359x360_360x181_may2016.8bin'
      end if
!
#if defined HYCOM_RES_387x360x26
#endif
#if defined HYCOM_RES_387x360x32
#endif
#if defined HYCOM_RES_359x360x26
#endif
#if defined HYCOM_RES_359x360x32
#endif
      return
      end subroutine const
end module const_proc
