      module const_proc
      implicit none
      real, parameter :: spval = -0.03125, onem=9806., onecm=98.06
     .    ,g=9.806, rho =1000., dz=1., spcifh=4185, flag=-999.
     .    ,epsil_p=0.01                             ! 1cm cutoff
      integer, parameter :: iorign=1, jorign=1,julian=365
      logical, parameter :: diag=.true., rhodot=.false.
      integer, dimension(12),parameter :: 
c    .   jdofm=(/30,60,90,120,150,180,210,240,270,300,330,360/)
     .   jdofm=(/31,59,90,120,151,181,212,243,273,304,334,365/)
      integer, dimension(12), parameter :: 
     .   monlg=(/31,28,31,30,31,30,31,31,30,31,30,31/)
      integer, parameter :: mo1=1,mo2=12
      character*3,dimension(13),parameter:: 
     .  amon=(/'JAN','FEB','MAR','APR','MAY','JUN',
     .         'JUL','AUG','SEP','OCT','NOV','DEC','ANN'/)
      real, dimension(33), parameter ::
     .   z33=(/0.,  10.,  20.,  30.,  50.,  75., 100., 125., 150., 200.,
     .       250., 300., 400., 500., 600., 700., 800., 900.,1000.,1100.,
     .      1200.,1300.,1400.,1500.,1750.,2000.,2500.,3000.,3500.,4000.,
     .      4500.,5000.,5500./)

      integer, parameter :: lp=6,itest=331,jtest=181
     .                     ,iatest=181,jatest=31
      character(len=128) :: 
     .  path0='/discover/nobackup/projects/giss/prod_input_files/' ! directory input files for modelE
     . ,path1='./'   ! directory hycom input files (.out)
     . ,path2='./'   ! directory output diagnostic files 
     . ,hycomtopo='depth387x360.4bin_1'
     . ,latlonij='latlon387x360.4bin'
     . ,basinmask='ibasin387x360.txt_1'
     . ,flnmcoso='cososino387x360.8bin'
     . ,flnmo2a='ssto2a_2_1deg.8bin'
c
      integer, parameter :: indoi=252,indoj1=118,indoj2=133 ! Indo Throughflow
     .  ,idrk1=323,jdrk1=292,idrk2=342,jdrk2=299            ! Drake Passage
     .  ,iberi=138, jberi=189                               ! Bering Strait
     .  ,ikuro1=197, ikuro2=198, jkuro=129                  ! Kuroshio
     .  ,igulf1=202, igulf2=202, jgulf=281                  ! Gulf Stream
     .  ,imed1=189, imed2=191, jmed=356                     ! Med outflow
c
      character(len=20) :: runid="Exxx"
      integer           :: ny1=1800, ny2=1800
      logical           :: monave_convert,solo_convert,timav,cnvert

      end module const_proc
