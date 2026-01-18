
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: readinput
! !INTERFACE:
subroutine readinput
! !USES:
use modmain
use moddftu
use modrdm
use modphonon
use modtest
use modrandom
use modpw
use modtddft
use modulr
use modvars
use modgw
use modbog
use modw90
use modmpi
use modomp
use modramdisk
! !DESCRIPTION:
!   Reads in the input parameters from the file {\tt elk.in}. Also sets default
!   values for the input parameters.
!
! !REVISION HISTORY:
!   Created September 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
logical lv
integer is,ia,ias,ios
integer i,j,k,l
real(8) sc,sc1,sc2,sc3
real(8) scx,scy,scz
real(8) scu,scu1,scu2,scu3
real(8) solscf,zn
real(8) axang(4),rot(3,3)
real(8) v(3),t1
character(256) block,symb,str

!------------------------!
!     default values     !
!------------------------!
ntasks=0
avec(:,:)=0.d0
avec(1,1)=1.d0
avec(2,2)=1.d0
avec(3,3)=1.d0
davec(:,:)=0.d0
sc=1.d0
sc1=1.d0
sc2=1.d0
sc3=1.d0
scx=1.d0
scy=1.d0
scz=1.d0
epslat=1.d-6
primcell=.false.
tshift=.true.
ngridk(:)=1
dngridk(:)=0
vkloff(:)=0.d0
autokpt=.false.
radkpt=40.d0
reducek=1
ngridq(:)=-1
reduceq=1
rgkmax=7.d0
drgkmax=0.d0
gmaxvr=12.d0
dgmaxvr=0.d0
lmaxapw=8
dlmaxapw=0
lmaxo=6
dlmaxo=0
lmaxi=1
fracinr=0.01d0
trhonorm=.true.
xctype(1)=3
xctype(2:3)=0
xctsp(1)=3
xctsp(2:3)=0
ktype(1)=51
ktype(2:3)=0
stype=3
swidth=0.001d0
autoswidth=.false.
mstar=10.d0
epsocc=1.d-8
epschg=1.d-3
nempty0=4.d0
dnempty0=0.d0
maxscl=200
mixtype=3
amixpm(1)=0.05d0
amixpm(2)=1.d0
! Broyden parameters recommended by M. Meinert
mixsdb=5
broydpm(1)=0.4d0
broydpm(2)=0.15d0
epspot=1.d-6
epsengy=1.d-4
epsforce=5.d-3
epsstress=2.d-3
molecule=.false.
nspecies=0
natoms(:)=0
atposl(:,:,:)=0.d0
datposl(:,:,:)=0.d0
atposc(:,:,:)=0.d0
bfcmt0(:,:,:)=0.d0
sppath=''
scrpath=''
nvp1d=2
if (allocated(vvlp1d)) deallocate(vvlp1d)
allocate(vvlp1d(3,nvp1d))
vvlp1d(:,1)=0.d0
vvlp1d(:,2)=1.d0
npp1d=200
vclp2d(:,:)=0.d0
vclp2d(1,1)=1.d0
vclp2d(2,2)=1.d0
np2d(:)=40
vclp3d(:,:)=0.d0
vclp3d(1,1)=1.d0
vclp3d(2,2)=1.d0
vclp3d(3,3)=1.d0
np3d(:)=20
nwplot=500
ngrkf=100
nswplot=1
wplot(1)=-0.5d0
wplot(2)=0.5d0
dosocc=.false.
dosmsum=.false.
dosssum=.false.
lmirep=.true.
spinpol=.false.
spinorb=.false.
socscf=1.d0
atpopt=1
maxatpstp=200
tau0atp=0.25d0
deltast=0.005d0
latvopt=0
maxlatvstp=30
tau0latv=0.2d0
lradstp=4
chgexs=0.d0
scissor=0.d0
noptcomp=1
optcomp(:,1)=1
intraband=.false.
evaltol=-1.d0
epsband=1.d-12
demaxbnd=2.5d0
autolinengy=.false.
dlefe=-0.1d0
deapwlo=0.05d0
bfieldc0(:)=0.d0
efieldc(:)=0.d0
afieldc(:)=0.d0
fsmtype=0
momfix(:)=0.d0
dmomfix(:)=0.d0
mommtfix(:,:,:)=1.d6
taufsm=0.01d0
rmtdelta=0.05d0
isgkmax=-1
symtype=1
deltaph=0.01d0
nphwrt=1
if (allocated(vqlwrt)) deallocate(vqlwrt)
allocate(vqlwrt(3,nphwrt))
vqlwrt(:,:)=0.d0
notelns=0
tforce=.false.
maxitoep=300
tauoep(1)=0.05d0
tauoep(2)=0.95d0
nkstlist=1
kstlist(:,1)=1
vklem(:)=0.d0
deltaem=0.025d0
ndspem=1
nosource=.false.
spinsprl=.false.
ssdph=.true.
vqlss(:)=0.d0
nwrite=0
dftu=0
inpdftu=1
ndftu=0
ujdu(:,:)=0.d0
fdu(:,:)=0.d0
edu(:,:)=0.d0
lamdu(:)=0.d0
udufix(:)=0.d0
dudufix(:)=0.d0
tmwrite=.false.
rdmxctype=2
rdmmaxscl=2
maxitn=200
maxitc=0
taurdmn=0.5d0
taurdmc=0.25d0
rdmalpha=0.656d0
rdmtemp=0.d0
reducebf=1.d0
ptnucl=.true.
tefvr=.true.
tefvit=.false.
minitefv=6
maxitefv=4
befvit=0.25d0
epsefvit=1.d-5
vecql(:)=0.d0
mustar=0.15d0
sqados(1:2)=0.d0
sqados(3)=1.d0
test=.false.
spincore=.false.
solscf=1.d0
emaxelnes=-1.2d0
wsfac(1)=-1.1d6; wsfac(2)=1.1d6
vhmat(:,:)=0.d0
vhmat(1,1)=1.d0
vhmat(2,2)=1.d0
vhmat(3,3)=1.d0
reduceh=.true.
hybrid0=.false.
hybridc=1.d0
ecvcut=-3.5d0
esccut=-0.4d0
gmaxrf=3.d0
emaxrf=1.d6
ntemp=40
nvbse0=2
ncbse0=3
nvxbse=0
ncxbse=0
bsefull=.false.
hxbse=.true.
hdbse=.true.
fxctype=-1
fxclrc(1)=0.d0
fxclrc(2)=0.d0
rndatposc=0.d0
rndbfcmt=0.d0
rndavec=0.d0
c_tb09=0.d0
tc_tb09=.false.
hmaxvr=20.d0
hkmax=12.d0
lorbcnd=.false.
lorbordc=3
nrmtscf=1.d0
dnrmtscf=0.d0
lmaxdos=3
epsdev=0.005d0
msmooth=0
npmae0=-1
wrtvars=.false.
ftmtype=0
ntmfix=0
tauftm=0.1d0
cmagz=.false.
axang(:)=0.d0
dncgga=1.d-8
tstime=1000.d0
dtimes=0.1d0
npulse=0
nramp=0
ntswrite(1)=500
ntswrite(2)=1
nxoapwlo=0
nxlo=0
tdrho1d=.false.
tdrho2d=.false.
tdrho3d=.false.
tdmag1d=.false.
tdmag2d=.false.
tdmag3d=.false.
tdjr1d=.false.
tdjr2d=.false.
tdjr3d=.false.
tddos=.false.
tdlsj=.false.
rndevt0=0.d0
sxcscf=1.d0
avecu(:,:)=0.d0
avecu(1,1)=1.d0
avecu(2,2)=1.d0
avecu(3,3)=1.d0
scu=1.d0
scu1=1.d0
scu2=1.d0
scu3=1.d0
q0cut=0.d0
rndbfcu=0.d0
efielduc(:)=0.d0
tplotq0=.true.
trdvclr=.false.
wmaxgw=-10.d0
twdiag=.false.
tsediag=.true.
actype=10
npole=3
nspade=100
tfav0=.true.
rmtall=-1.d0
maxthd=0
maxthd1=0
maxthdmkl=0
maxlvl=4
tdphi=0.d0
thetamld=45.d0*pi/180.d0
ntsbackup=0
! Wannier90 variables
seedname='wannier'
num_wann=0
num_bands=0
num_iter=400
trial_step=1.d-3
nxlwin=0
tbdip=.false.
tjr=.false.
tauefm=0.01d0
epsefm=1.d-6
ehfb=1.d0
t0gclq0=.false.
tafindt=.false.
afindpm(:)=0.d0
afindpm(2)=1.d0
nkspolar=4
ntsforce=10
wphcut=1.d-6
ephscf(1)=8.d0
ephscf(2)=0.02d0
anomalous=.false.
tephde=.false.
bdiag=.false.
ecutb=0.001d0
ediag=.false.
pwxpsn=2
ramdisk=.true.
wrtdsk=.true.
epsdmat=1.d-8
tm3old=.false.
batch=.false.

!--------------------------!
!     read from elk.in     !
!--------------------------!
open(50,file='elk.in',status='OLD',form='FORMATTED',iostat=ios)
if (ios.ne.0) then
  write(*,*)
  write(*,'("Error(readinput): error opening elk.in")')
  write(*,*)
  stop
end if
10 continue
read(50,*,end=30) block
! check for a comment
if ((scan(trim(block),'!').eq.1).or.(scan(trim(block),'#').eq.1)) goto 10
select case(trim(block))
case('tasks')
  do i=1,maxtasks
    read(50,'(A)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): no tasks to perform")')
        write(*,*)
        stop
      end if
      ntasks=i-1
      goto 10
    end if
    read(str,*,iostat=ios) tasks(i)
    if (ios.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading tasks")')
      write(*,'("(blank line required after tasks block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): too many tasks")')
  write(*,'("Adjust maxtasks in modmain and recompile code")')
  write(*,*)
  stop
case('species')
! generate a species file
  call genspecies(50)
case('fspecies')
! generate fractional species files
  do is=1,maxspecies
    read(50,'(A)',err=20) str
    if (trim(str).eq.'') goto 10
    read(str,*,iostat=ios) zn,symb
    if (ios.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading fractional species")')
      write(*,'("(blank line required after fspecies block)")')
      write(*,*)
      stop
    end if
    if (zn.gt.-1.d0+epsocc) then
      write(*,*)
      write(*,'("Error(readinput): fractional nuclear Z > -1 : ",G18.10)') zn
      write(*,*)
      stop
    end if
    call genfspecies(zn,symb)
  end do
  write(*,*)
  write(*,'("Error(readinput): too many fractional nucleus species")')
  write(*,*)
  stop
case('avec')
  read(50,'(A)',err=20) str
  read(str,*,err=20) avec(:,1)
  read(str,*,iostat=ios) avec(:,1),davec(:,1)
  read(50,'(A)',err=20) str
  read(str,*,err=20) avec(:,2)
  read(str,*,iostat=ios) avec(:,2),davec(:,2)
  read(50,'(A)',err=20) str
  read(str,*,err=20) avec(:,3)
  read(str,*,iostat=ios) avec(:,3),davec(:,3)
case('scale')
  read(50,*,err=20) sc
case('scale1')
  read(50,*,err=20) sc1
case('scale2')
  read(50,*,err=20) sc2
case('scale3')
  read(50,*,err=20) sc3
case('scalex')
  read(50,*,err=20) scx
case('scaley')
  read(50,*,err=20) scy
case('scalez')
  read(50,*,err=20) scz
case('epslat')
  read(50,*,err=20) epslat
  if (epslat.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epslat <= 0 : ",G18.10)') epslat
    write(*,*)
    stop
  end if
case('primcell')
  read(50,*,err=20) primcell
case('tshift')
  read(50,*,err=20) tshift
case('autokpt')
  read(50,*,err=20) autokpt
case('radkpt')
  read(50,*,err=20) radkpt
  if (radkpt.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): radkpt <= 0 : ",G18.10)') radkpt
    write(*,*)
    stop
  end if
case('ngridk')
  read(50,'(A)',err=20) str
  read(str,*,err=20) ngridk(:)
  read(str,*,iostat=ios) ngridk(:),dngridk(:)
  if (any(ngridk(:).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridk : ",3I8)') ngridk
    write(*,*)
    stop
  end if
  autokpt=.false.
case('vkloff')
  read(50,*,err=20) vkloff(:)
  if (any(vkloff(:).lt.0.d0).or.any(vkloff(:).ge.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): vkloff components should be in [0,1) : ",&
     &3G18.10)') vkloff
    write(*,*)
    stop
  end if
case('reducek')
  read(50,*,err=20) reducek
case('ngridq')
  read(50,*,err=20) ngridq(:)
  if (any(ngridq(:).le.0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid ngridq : ",3I8)') ngridq
    write(*,*)
    stop
  end if
case('reduceq')
  read(50,*,err=20) reduceq
case('rgkmax')
  read(50,'(A)',err=20) str
  read(str,*,err=20) rgkmax
  read(str,*,iostat=ios) rgkmax,drgkmax
  if (rgkmax.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rgkmax <= 0 : ",G18.10)') rgkmax
    write(*,*)
    stop
  end if
case('gmaxvr')
  read(50,'(A)',err=20) str
  read(str,*,err=20) gmaxvr
  read(str,*,iostat=ios) gmaxvr,dgmaxvr
case('lmaxapw')
  read(50,'(A)',err=20) str
  read(str,*,err=20) lmaxapw
  read(str,*,iostat=ios) lmaxapw,dlmaxapw
  if (lmaxapw.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxapw < 0 : ",I8)') lmaxapw
    write(*,*)
    stop
  end if
  if (lmaxapw.ge.maxlapw) then
    write(*,*)
    write(*,'("Error(readinput): lmaxapw too large : ",I8)') lmaxapw
    write(*,'("Adjust maxlapw in modmain and recompile code")')
    write(*,*)
    stop
  end if
case('lmaxo','lmaxvr')
  read(50,'(A)',err=20) str
  read(str,*,err=20) lmaxo
  read(str,*,iostat=ios) lmaxo,dlmaxo
  if (lmaxo.lt.3) then
    write(*,*)
    write(*,'("Error(readinput): lmaxo < 3 : ",I8)') lmaxo
    write(*,*)
    stop
  end if
case('lmaxi','lmaxinr')
  read(50,*,err=20) lmaxi
  if (lmaxi.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): lmaxi < 1 : ",I8)') lmaxi
    write(*,*)
    stop
  end if
case('lmaxmat')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''lmaxmat'' is no longer used")')
case('fracinr')
  read(50,*,err=20) fracinr
case('trhonorm')
  read(50,*,err=20) trhonorm
case('spinpol')
  read(50,*,err=20) spinpol
case('spinorb')
  read(50,*,err=20) spinorb
case('socscf')
  read(50,*,err=20) socscf
  if (socscf.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): socscf < 0 : ",G18.10)') socscf
    write(*,*)
    stop
  end if
case('xctype')
  read(50,'(A)',err=20) str
  str=trim(str)//' 0 0'
  read(str,*,err=20) xctype
case('xctsp')
  read(50,'(A)',err=20) str
  str=trim(str)//' 0 0'
  read(str,*,err=20) xctsp
case('ktype')
  read(50,'(A)',err=20) str
  str=trim(str)//' 0 0'
  read(str,*,err=20) ktype
  if (ktype(3).ne.0) then
    write(*,*)
    write(*,'("Error(readinput): ktype(3) should be zero : ",I8)') ktype(3)
    write(*,*)
    stop
  end if
case('stype')
  read(50,*,err=20) stype
case('swidth')
  read(50,*,err=20) swidth
  if (swidth.lt.1.d-9) then
    write(*,*)
    write(*,'("Error(readinput): swidth too small or negative : ",G18.10)') &
     swidth
    write(*,*)
    stop
  end if
case('autoswidth')
  read(50,*,err=20) autoswidth
case('mstar')
  read(50,*,err=20) mstar
  if (mstar.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): mstar <= 0 : ",G18.10)') mstar
    write(*,*)
    stop
  end if
case('epsocc')
  read(50,*,err=20) epsocc
  if (epsocc.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsocc <= 0 : ",G18.10)') epsocc
    write(*,*)
    stop
  end if
case('epschg')
  read(50,*,err=20) epschg
  if (epschg.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epschg <= 0 : ",G18.10)') epschg
    write(*,*)
    stop
  end if
case('nempty','nempty0')
  read(50,'(A)',err=20) str
  read(str,*,err=20) nempty0
  read(str,*,iostat=ios) nempty0,dnempty0
  if (nempty0.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): nempty <= 0 : ",G18.10)') nempty0
    write(*,*)
    stop
  end if
case('mixtype')
  read(50,*,err=20) mixtype
case('amixpm','beta0','betamax')
  if (trim(block).eq.'amixpm') then
    read(50,*,err=20) amixpm(:)
  else if (trim(block).eq.'beta0') then
    read(50,*,err=20) amixpm(1)
  else
    read(50,*,err=20) amixpm(2)
  end if
  if (amixpm(1).lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): beta0 [amixpm(1)] < 0 : ",G18.10)') amixpm(1)
    write(*,*)
    stop
  end if
  if ((amixpm(2).lt.0.d0).or.(amixpm(2).gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): betamax [amixpm(2)] not in [0,1] : ",G18.10)')&
     amixpm(2)
    write(*,*)
    stop
  end if
case('mixsdb')
  read(50,*,err=20) mixsdb
  if (mixsdb.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): mixsdb < 2 : ",I8)') mixsdb
    write(*,*)
    stop
  end if
case('broydpm')
  read(50,*,err=20) broydpm(:)
  if ((broydpm(1).lt.0.d0).or.(broydpm(1).gt.1.d0).or. &
      (broydpm(2).lt.0.d0).or.(broydpm(2).gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid Broyden mixing parameters : ",&
     &2G18.10)') broydpm
    write(*,*)
    stop
  end if
case('maxscl')
  read(50,*,err=20) maxscl
  if (maxscl.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): maxscl < 0 : ",I8)') maxscl
    write(*,*)
    stop
  end if
case('epspot')
  read(50,*,err=20) epspot
case('epsengy')
  read(50,*,err=20) epsengy
case('epsforce')
  read(50,*,err=20) epsforce
case('epsstress')
  read(50,*,err=20) epsstress
case('sppath')
  read(50,*,err=20) sppath
  sppath=adjustl(sppath)
case('scrpath')
  read(50,*,err=20) scrpath
case('molecule')
  read(50,*,err=20) molecule
case('atoms')
  read(50,*,err=20) nspecies
  if (nspecies.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nspecies <= 0 : ",I8)') nspecies
    write(*,*)
    stop
  end if
  if (nspecies.gt.maxspecies) then
    write(*,*)
    write(*,'("Error(readinput): nspecies too large : ",I8)') nspecies
    write(*,'("Adjust maxspecies in modmain and recompile code")')
    write(*,*)
    stop
  end if
  do is=1,nspecies
    read(50,*,err=20) spfname(is)
    spfname(is)=adjustl(spfname(is))
    read(50,*,err=20) natoms(is)
    if (natoms(is).le.0) then
      write(*,*)
      write(*,'("Error(readinput): natoms <= 0 : ",I8)') natoms(is)
      write(*,'(" for species ",I4)') is
      write(*,*)
      stop
    end if
    if (natoms(is).gt.maxatoms) then
      write(*,*)
      write(*,'("Error(readinput): natoms too large : ",I8)') natoms(is)
      write(*,'(" for species ",I4)') is
      write(*,'("Adjust maxatoms in modmain and recompile code")')
      write(*,*)
      stop
    end if
    do ia=1,natoms(is)
      read(50,'(A)',err=20) str
      read(str,*,err=20) atposl(:,ia,is)
      read(str,*,iostat=ios) atposl(:,ia,is),bfcmt0(:,ia,is),datposl(:,ia,is)
    end do
  end do
case('plot1d')
  read(50,*,err=20) nvp1d,npp1d
  if (nvp1d.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): nvp1d < 1 : ",I8)') nvp1d
    write(*,*)
    stop
  end if
  if (npp1d.lt.nvp1d) then
    write(*,*)
    write(*,'("Error(readinput): npp1d < nvp1d : ",2I8)') npp1d,nvp1d
    write(*,*)
    stop
  end if
  if (allocated(vvlp1d)) deallocate(vvlp1d)
  allocate(vvlp1d(3,nvp1d))
  do i=1,nvp1d
    read(50,*,err=20) vvlp1d(:,i)
  end do
case('plot2d')
  read(50,*,err=20) vclp2d(:,0)
  read(50,*,err=20) vclp2d(:,1)
  read(50,*,err=20) vclp2d(:,2)
  read(50,*,err=20) np2d(:)
  if ((np2d(1).lt.1).or.(np2d(2).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): np2d < 1 : ",2I8)') np2d
    write(*,*)
    stop
  end if
case('plot3d')
  read(50,*,err=20) vclp3d(:,0)
  read(50,*,err=20) vclp3d(:,1)
  read(50,*,err=20) vclp3d(:,2)
  read(50,*,err=20) vclp3d(:,3)
  read(50,*,err=20) np3d(:)
  if ((np3d(1).lt.1).or.(np3d(2).lt.1).or.(np3d(3).lt.1)) then
    write(*,*)
    write(*,'("Error(readinput): np3d < 1 : ",3I8)') np3d
    write(*,*)
    stop
  end if
case('wplot','dos')
  read(50,*,err=20) nwplot,ngrkf,nswplot
  if (nwplot.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): nwplot < 2 : ",I8)') nwplot
    write(*,*)
    stop
  end if
  if (ngrkf.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): ngrkf < 1 : ",I8)') ngrkf
    write(*,*)
    stop
  end if
  if (nswplot.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nswplot < 0 : ",I8)') nswplot
    write(*,*)
    stop
  end if
  read(50,*,err=20) wplot(:)
  if (wplot(1).gt.wplot(2)) then
    write(*,*)
    write(*,'("Error(readinput): wplot(1) > wplot(2) : ",2G18.10)') wplot
    write(*,*)
    stop
  end if
case('dosocc')
  read(50,*,err=20) dosocc
case('dosmsum')
  read(50,*,err=20) dosmsum
case('dosssum')
  read(50,*,err=20) dosssum
case('lmirep')
  read(50,*,err=20) lmirep
case('atpopt')
  read(50,*,err=20) atpopt
case('maxatpstp','maxatmstp')
  read(50,*,err=20) maxatpstp
  if (maxatpstp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): maxatpstp <= 0 : ",I8)') maxatpstp
    write(*,*)
    stop
  end if
case('tau0atp','tau0atm')
  read(50,*,err=20) tau0atp
case('deltast')
  read(50,*,err=20) deltast
  if (deltast.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): deltast <= 0 : ",G18.10)') deltast
    write(*,*)
    stop
  end if
case('latvopt')
  read(50,*,err=20) latvopt
case('maxlatvstp')
  read(50,*,err=20) maxlatvstp
  if (maxlatvstp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): maxlatvstp <= 0 : ",I8)') maxlatvstp
    write(*,*)
    stop
  end if
case('tau0latv')
  read(50,*,err=20) tau0latv
case('nstfsp')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''nstfsp'' is no longer used")')
case('lradstp')
  read(50,*,err=20) lradstp
  if (lradstp.le.0) then
    write(*,*)
    write(*,'("Error(readinput): lradstp <= 0 : ",I8)') lradstp
    write(*,*)
    stop
  end if
case('chgexs')
  read(50,*,err=20) chgexs
case('nprad')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''nprad'' is no longer used")')
case('scissor')
  read(50,*,err=20) scissor
case('optcomp')
  do i=1,27
    read(50,'(A)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty optical component list")')
        write(*,*)
        stop
      end if
      noptcomp=i-1
      goto 10
    end if
    str=trim(str)//' 1 1'
    read(str,*,iostat=ios) optcomp(:,i)
    if (ios.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading optical component list")')
      write(*,'("(blank line required after optcomp block)")')
      write(*,*)
      stop
    end if
    if (any(optcomp(:,i).lt.1).or.any(optcomp(:,i).gt.3)) then
      write(*,*)
      write(*,'("Error(readinput): invalid optcomp : ",3I8)') optcomp(:,i)
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): optical component list too long")')
  write(*,*)
  stop
case('intraband')
  read(50,*,err=20) intraband
case('evaltol')
  read(50,*,err=20) evaltol
case('deband')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''deband'' is no longer used")')
case('epsband')
  read(50,*,err=20) epsband
  if (epsband.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsband <= 0 : ",G18.10)') epsband
    write(*,*)
    stop
  end if
case('demaxbnd')
  read(50,*,err=20) demaxbnd
  if (demaxbnd.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): demaxbnd <= 0 : ",G18.10)') demaxbnd
    write(*,*)
    stop
  end if
case('autolinengy')
  read(50,*,err=20) autolinengy
case('dlefe')
  read(50,*,err=20) dlefe
case('deapwlo')
  read(50,*,err=20) deapwlo
  if (abs(deapwlo).lt.1.d-8) then
    write(*,*)
    write(*,'("Error(readinput): invalid deapwlo : ",G18.10)') deapwlo
    write(*,*)
    stop
  end if
case('bfieldc')
  read(50,*,err=20) bfieldc0(:)
case('efieldc')
  read(50,*,err=20) efieldc(:)
case('afieldc')
  read(50,*,err=20) afieldc(:)
case('fsmtype','fixspin')
  read(50,*,err=20) fsmtype
case('momfix')
  read(50,'(A)',err=20) str
  read(str,*,err=20) momfix(:)
  read(str,*,iostat=ios) momfix(:),dmomfix(:)
case('mommtfix')
  do ias=1,maxspecies*maxatoms
    read(50,'(A)',err=20) str
    if (trim(str).eq.'') goto 10
    read(str,*,iostat=ios) is,ia,mommtfix(:,ia,is)
    if (ios.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading muffin-tin fixed spin &
       &moments")')
      write(*,'("(blank line required after mommtfix block")')
      write(*,*)
      stop
    end if
  end do
case('taufsm')
  read(50,*,err=20) taufsm
  if (taufsm.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taufsm < 0 : ",G18.10)') taufsm
    write(*,*)
    stop
  end if
case('autormt')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''autormt'' is no longer used")')
case('rmtdelta')
  read(50,*,err=20) rmtdelta
  if (rmtdelta.lt.0.d0) then
    write(*,*)
    write(*,'("Warning(readinput): rmtdelta < 0 : ",G18.10)') rmtdelta
  end if
case('isgkmax')
  read(50,*,err=20) isgkmax
case('nosym')
  read(50,*,err=20) lv
  if (lv) symtype=0
case('symtype')
  read(50,*,err=20) symtype
  if ((symtype.lt.0).or.(symtype.gt.2)) then
    write(*,*)
    write(*,'("Error(readinput): symtype not defined : ",I8)') symtype
    write(*,*)
    stop
  end if
case('deltaph')
  read(50,*,err=20) deltaph
  if (deltaph.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): deltaph <= 0 : ",G18.10)') deltaph
    write(*,*)
    stop
  end if
case('phwrite')
  read(50,*,err=20) nphwrt
  if (nphwrt.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nphwrt <= 0 : ",I8)') nphwrt
    write(*,*)
    stop
  end if
  if (allocated(vqlwrt)) deallocate(vqlwrt)
  allocate(vqlwrt(3,nphwrt))
  do i=1,nphwrt
    read(50,*,err=20) vqlwrt(:,i)
  end do
case('notes')
  if (allocated(notes)) deallocate(notes)
  allocate(notes(0))
  notelns=0
  do
    read(50,'(A)') str
    if (trim(str).eq.'') goto 10
    notelns=notelns+1
    call addstr(notes)
  end do
case('tforce')
  read(50,*,err=20) tforce
case('tfibs')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''tfibs'' is no longer used")')
case('maxitoep')
  read(50,*,err=20) maxitoep
  if (maxitoep.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): maxitoep < 1 : ",I8)') maxitoep
    write(*,*)
    stop
  end if
case('tauoep')
  read(50,*,err=20) tauoep(:)
  if ((tauoep(1).lt.0.d0).or.(tauoep(2).lt.0.d0)) then
    write(*,*)
    write(*,'("Error(readinput): tauoep < 0 : ",3G18.10)') tauoep
    write(*,*)
    stop
  end if
case('kstlist')
  do i=1,maxkst
    read(50,'(A)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty k-point and state list")')
        write(*,*)
        stop
      end if
      nkstlist=i-1
      goto 10
    end if
    read(str,*,iostat=ios) kstlist(:,i)
    if (ios.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading k-point and state list")')
      write(*,'("(blank line required after kstlist block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): k-point and state list too long")')
  write(*,*)
  stop
case('vklem')
  read(50,*,err=20) vklem
case('deltaem')
  read(50,*,err=20) deltaem
case('ndspem')
  read(50,*,err=20) ndspem
  if ((ndspem.lt.1).or.(ndspem.gt.4)) then
    write(*,*)
    write(*,'("Error(readinput): ndspem out of range : ",I8)') ndspem
    write(*,*)
    stop
  end if
case('nosource')
  read(50,*,err=20) nosource
case('spinsprl')
  read(50,*,err=20) spinsprl
case('ssdph')
  read(50,*,err=20) ssdph
case('vqlss')
  read(50,*,err=20) vqlss
case('nwrite')
  read(50,*,err=20) nwrite
case('DFT+U','dft+u','lda+u')
  read(50,*,err=20) dftu,inpdftu
  do i=1,maxdftu
    read(50,'(A)',err=20) str
    if (trim(str).eq.'') then
      ndftu=i-1
      goto 10
    end if
    select case(inpdftu)
    case(1)
      read(str,*,iostat=ios) is,l,ujdu(1:2,i)
    case(2)
      read(str,*,iostat=ios) is,l,(fdu(k,i),k=0,2*l,2)
    case(3)
      read(str,*,iostat=ios) is,l,(edu(k,i),k=0,l)
    case(4)
      read(str,*,iostat=ios) is,l,lamdu(i)
    case(5)
      read(str,*,iostat=ios) is,l,udufix(i),dudufix(i)
      read(str,*,iostat=ios) is,l,udufix(i)
    case default
      write(*,*)
      write(*,'("Error(readinput): invalid inpdftu : ",I8)') inpdftu
      write(*,*)
      stop
    end select
    if (ios.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading DFT+U parameters")')
      write(*,'("(blank line required after dft+u block)")')
      write(*,*)
      stop
    end if
    if ((is.le.0).or.(is.ge.maxspecies)) then
      write(*,*)
      write(*,'("Error(readinput): invalid species number in dft+u block : ", &
       &I8)') is
      write(*,*)
      stop
    end if
    if (l.lt.0) then
      write(*,*)
      write(*,'("Error(readinput): l < 0 in dft+u block : ",I8)') l
      write(*,*)
      stop
    end if
    if (l.gt.lmaxdm) then
      write(*,*)
      write(*,'("Error(readinput): l > lmaxdm in dft+u block : ",2I8)') l,lmaxdm
      write(*,*)
      stop
    end if
! check for repeated entries
    do j=1,i-1
      if ((is.eq.isldu(1,j)).and.(l.eq.isldu(2,j))) then
        write(*,*)
        write(*,'("Error(readinput): repeated entry in DFT+U block")')
        write(*,*)
        stop
      end if
    end do
    isldu(1,i)=is
    isldu(2,i)=l
  end do
  write(*,*)
  write(*,'("Error(readinput): too many DFT+U entries")')
  write(*,'("Adjust maxdftu in modmain and recompile code")')
  write(*,*)
  stop
case('tmwrite','tmomlu')
  read(50,*,err=20) tmwrite
case('readadu','readalu')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''readadu'' is no longer used")')
case('rdmxctype')
  read(50,*,err=20) rdmxctype
case('rdmmaxscl')
  read(50,*,err=20) rdmmaxscl
  if (rdmmaxscl.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): rdmmaxscl < 0 : ",I8)') rdmmaxscl
    write(*,*)
  end if
case('maxitn')
  read(50,*,err=20) maxitn
case('maxitc')
  read(50,*,err=20) maxitc
case('taurdmn')
  read(50,*,err=20) taurdmn
  if (taurdmn.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taurdmn < 0 : ",G18.10)') taurdmn
    write(*,*)
    stop
  end if
case('taurdmc')
  read(50,*,err=20) taurdmc
  if (taurdmc.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): taurdmc < 0 : ",G18.10)') taurdmc
    write(*,*)
    stop
  end if
case('rdmalpha')
  read(50,*,err=20) rdmalpha
  if ((rdmalpha.le.0.d0).or.(rdmalpha.ge.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): rdmalpha not in (0,1) : ",G18.10)') rdmalpha
    write(*,*)
    stop
  end if
case('rdmtemp')
  read(50,*,err=20) rdmtemp
  if (rdmtemp.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): rdmtemp < 0 : ",G18.10)') rdmtemp
    write(*,*)
    stop
  end if
case('reducebf')
  read(50,*,err=20) reducebf
  if ((reducebf.lt.0.49d0).or.(reducebf.gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): reducebf not in [0.5,1] : ",G18.10)') reducebf
    write(*,*)
    stop
  end if
case('ptnucl')
  read(50,*,err=20) ptnucl
case('tefvr','tseqr')
  read(50,*,err=20) tefvr
case('tefvit','tseqit')
  read(50,*,err=20) tefvit
case('minitefv','minseqit')
  read(50,*,err=20) minitefv
  if (minitefv.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): minitefv < 1 : ",I8)') minitefv
    write(*,*)
    stop
  end if
case('maxitefv','maxseqit')
  read(50,*,err=20) maxitefv
  if (maxitefv.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): maxitefv < 1 : ",I8)') maxitefv
    write(*,*)
    stop
  end if
case('befvit','bseqit')
  read(50,*,err=20) befvit
  if (befvit.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): befvit <= 0 : ",G18.10)') befvit
    write(*,*)
    stop
  end if
case('epsefvit','epsseqit')
  read(50,*,err=20) epsefvit
  if (epsefvit.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsefvit < 0 : ",G18.10)') epsefvit
    write(*,*)
    stop
  end if
case('nseqit')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''nseqit'' is no longer used")')
case('tauseq')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''tauseq'' is no longer used")')
case('vecql')
  read(50,*,err=20) vecql(:)
case('mustar')
  read(50,*,err=20) mustar
case('sqados')
  read(50,*,err=20) sqados(:)
case('test')
  read(50,*,err=20) test
case('frozencr')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''frozencr'' is no longer used")')
case('spincore')
  read(50,*,err=20) spincore
case('solscf')
  read(50,*,err=20) solscf
  if (solscf.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): solscf < 0 : ",G18.10)') solscf
    write(*,*)
    stop
  end if
case('emaxelnes')
  read(50,*,err=20) emaxelnes
case('wsfac')
  read(50,*,err=20) wsfac(:)
case('vhmat')
  read(50,*,err=20) vhmat(1,:)
  read(50,*,err=20) vhmat(2,:)
  read(50,*,err=20) vhmat(3,:)
case('reduceh')
  read(50,*,err=20) reduceh
case('hybrid')
  read(50,*,err=20) hybrid0
case('hybridc','hybmix')
  read(50,*,err=20) hybridc
  if ((hybridc.lt.0.d0).or.(hybridc.gt.1.d0)) then
    write(*,*)
    write(*,'("Error(readinput): invalid hybridc : ",G18.10)') hybridc
    write(*,*)
    stop
  end if
case('ecvcut')
  read(50,*,err=20) ecvcut
case('esccut')
  read(50,*,err=20) esccut
case('nvbse')
  read(50,*,err=20) nvbse0
  if (nvbse0.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nvbse < 0 : ",I8)') nvbse0
    write(*,*)
    stop
  end if
case('ncbse')
  read(50,*,err=20) ncbse0
  if (ncbse0.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): ncbse < 0 : ",I8)') ncbse0
    write(*,*)
    stop
  end if
case('istxbse')
  do i=1,maxxbse
    read(50,'(A)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty BSE extra valence state list")')
        write(*,*)
        stop
      end if
      nvxbse=i-1
      goto 10
    end if
    read(str,*,iostat=ios) istxbse(i)
    if (ios.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading BSE valence state list")')
      write(*,'("(blank line required after istxbse block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): BSE extra valence state list too long")')
  write(*,*)
  stop
case('jstxbse')
  do i=1,maxxbse
    read(50,'(A)',err=20) str
    if (trim(str).eq.'') then
      if (i.eq.1) then
        write(*,*)
        write(*,'("Error(readinput): empty BSE extra conduction state list")')
        write(*,*)
        stop
      end if
      ncxbse=i-1
      goto 10
    end if
    read(str,*,iostat=ios) jstxbse(i)
    if (ios.ne.0) then
      write(*,*)
      write(*,'("Error(readinput): error reading BSE conduction state list")')
      write(*,'("(blank line required after jstxbse block)")')
      write(*,*)
      stop
    end if
  end do
  write(*,*)
  write(*,'("Error(readinput): BSE extra conduction state list too long")')
  write(*,*)
  stop
case('bsefull')
  read(50,*,err=20) bsefull
case('hxbse')
  read(50,*,err=20) hxbse
case('hdbse')
  read(50,*,err=20) hdbse
case('gmaxrf','gmaxrpa')
  read(50,*,err=20) gmaxrf
  if (gmaxrf.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): gmaxrf < 0 : ",G18.10)') gmaxrf
    write(*,*)
    stop
  end if
case('emaxrf')
  read(50,*,err=20) emaxrf
  if (emaxrf.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): emaxrf < 0 : ",G18.10)') emaxrf
    write(*,*)
    stop
  end if
case('fxctype')
  read(50,'(A)',err=20) str
  str=trim(str)//' 0 0'
  read(str,*,err=20) fxctype
case('fxclrc')
  read(50,'(A)',err=20) str
  str=trim(str)//' 0.0'
  read(str,*,err=20) fxclrc(:)
case('ntemp')
  read(50,*,err=20) ntemp
  if (ntemp.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): ntemp < 1 : ",I8)') ntemp
    write(*,*)
    stop
  end if
case('trimvg')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''trimvg'' is no longer used")')
  write(*,'(" setting msmooth=4 instead")')
  msmooth=4
case('rndstate','rndseed')
  read(50,*,err=20) rndstate(0)
  rndstate(0)=abs(rndstate(0))
case('rndatposc')
  read(50,*,err=20) rndatposc
case('rndbfcmt')
  read(50,*,err=20) rndbfcmt
case('rndavec')
  read(50,*,err=20) rndavec
case('c_tb09')
  read(50,*,err=20) c_tb09
! set flag to indicate Tran-Blaha constant has been read in
  tc_tb09=.true.
case('lowq','highq','vhighq','uhighq')
  read(50,*,err=20) lv
  if (lv) then
    if (trim(block).eq.'lowq') then
      rgkmax=6.5d0
      gmaxvr=10.d0
      lmaxapw=7
      lmaxo=5
      nxlo=2
      lorbcnd=.true.
      radkpt=25.d0
      autokpt=.true.
      vkloff(:)=0.5d0
      nempty0=4.d0
      epspot=1.d-5
      epsengy=5.d-4
      epsforce=1.d-2
      epsstress=3.d-3
      autolinengy=.true.
      gmaxrf=2.5d0
      lradstp=6
    else if (trim(block).eq.'highq') then
! parameter set for high quality calculation
      rgkmax=max(rgkmax,8.d0)
      gmaxvr=max(gmaxvr,16.d0)
      msmooth=max(msmooth,4)
      lmaxapw=max(lmaxapw,9)
      lmaxo=max(lmaxo,7)
      nrmtscf=max(nrmtscf,1.5d0)
      nxlo=max(nxlo,2)
      lorbcnd=.true.
      radkpt=max(radkpt,50.d0)
      autokpt=.true.
      vkloff(:)=0.d0
      nempty0=max(nempty0,10.d0)
      epspot=min(epspot,1.d-7)
      epsengy=min(epsengy,1.d-5)
      epsforce=min(epsforce,5.d-4)
      epsstress=min(epsstress,1.d-3)
      autolinengy=.true.
      gmaxrf=max(gmaxrf,4.d0)
    else if (trim(block).eq.'vhighq') then
! parameter set for very high quality calculation
      rgkmax=max(rgkmax,9.d0)
      gmaxvr=max(gmaxvr,18.d0)
      msmooth=max(msmooth,4)
      lmaxapw=max(lmaxapw,12)
      lmaxo=max(lmaxo,8)
      nrmtscf=max(nrmtscf,2.d0)
      nxlo=max(nxlo,3)
      lorbcnd=.true.
      radkpt=max(radkpt,90.d0)
      autokpt=.true.
      vkloff(:)=0.d0
      nempty0=max(nempty0,20.d0)
      epspot=min(epspot,1.d-7)
      epsengy=min(epsengy,1.d-6)
      epsforce=min(epsforce,2.d-4)
      epsstress=min(epsstress,5.d-4)
      autolinengy=.true.
      gmaxrf=max(gmaxrf,5.d0)
    else
! parameter set for ultra high quality calculation
      rgkmax=max(rgkmax,10.d0)
      gmaxvr=max(gmaxvr,20.d0)
      msmooth=max(msmooth,4)
      lmaxapw=max(lmaxapw,14)
      lmaxo=max(lmaxo,9)
      nrmtscf=max(nrmtscf,4.d0)
      nxlo=max(nxlo,3)
      lorbcnd=.true.
      radkpt=max(radkpt,120.d0)
      autokpt=.true.
      vkloff(:)=0.d0
      nempty0=max(nempty0,40.d0)
      epspot=min(epspot,1.d-7)
      epsengy=min(epsengy,1.d-6)
      epsforce=min(epsforce,1.d-4)
      epsstress=min(epsstress,2.d-4)
      autolinengy=.true.
      gmaxrf=max(gmaxrf,6.d0)
    end if
    if (mp_mpi) then
      write(*,*)
      write(*,'("Info(readinput): parameters set by ",A," option")') trim(block)
      write(*,'(" rgkmax : ",G18.10)') rgkmax
      write(*,'(" gmaxvr : ",G18.10)') gmaxvr
      write(*,'(" msmooth : ",I4)') msmooth
      write(*,'(" lmaxapw : ",I4)') lmaxapw
      write(*,'(" lmaxo : ",I4)') lmaxo
      write(*,'(" nrmtscf : ",G18.10)') nrmtscf
      write(*,'(" nxlo : ",I4)') nxlo
      write(*,'(" lorbcnd : ",L1)') lorbcnd
      write(*,'(" radkpt : ",G18.10)') radkpt
      write(*,'(" autokpt : ",L1)') autokpt
      write(*,'(" vkloff : ",3G18.10)') vkloff
      write(*,'(" nempty0 : ",G18.10)') nempty0
      write(*,'(" epspot : ",G18.10)') epspot
      write(*,'(" epsengy : ",G18.10)') epsengy
      write(*,'(" epsforce : ",G18.10)') epsforce
      write(*,'(" epsstress : ",G18.10)') epsstress
      write(*,'(" autolinengy : ",L1)') autolinengy
      write(*,'(" gmaxrf : ",G18.10)') gmaxrf
      if (trim(block).eq.'lowq') then
        write(*,'(" lradstp : ",I4)') lradstp
      end if
    end if
  end if
case('hmaxvr')
  read(50,*,err=20) hmaxvr
  if (hmaxvr.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): hmaxvr < 0 : ",G18.10)') hmaxvr
    write(*,*)
    stop
  end if
case('hkmax')
  read(50,*,err=20) hkmax
  if (hkmax.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): hkmax <= 0 : ",G18.10)') hkmax
    write(*,*)
    stop
  end if
case('lorbcnd')
  read(50,*,err=20) lorbcnd
case('lorbordc')
  read(50,*,err=20) lorbordc
  if (lorbordc.lt.2) then
    write(*,*)
    write(*,'("Error(readinput): lorbordc < 2 : ",I8)') lorbordc
    write(*,*)
    stop
  end if
  if (lorbordc.gt.maxlorbord) then
    write(*,*)
    write(*,'("Error(readinput): lorbordc too large : ",I8)') lorbordc
    write(*,'("Adjust maxlorbord in modmain and recompile code")')
    write(*,*)
    stop
  end if
case('nrmtscf')
  read(50,'(A)',err=20) str
  read(str,*,err=20) nrmtscf
  read(str,*,iostat=ios) nrmtscf,dnrmtscf
  if (nrmtscf.lt.0.5d0) then
    write(*,*)
    write(*,'("Error(readinput): nrmtscf < 0.5 : ",G18.10)') nrmtscf
    write(*,*)
    stop
  end if
case('lmaxdos')
  read(50,*,err=20) lmaxdos
  if (lmaxdos.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): lmaxdos < 0 : ",I8)') lmaxdos
    write(*,*)
    stop
  end if
case('epsdev')
  read(50,*,err=20) epsdev
  if (epsdev.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): epsdev <= 0 : ",G18.10)') epsdev
    write(*,*)
    stop
  end if
case('msmooth')
  read(50,*,err=20) msmooth
  if (msmooth.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): msmooth < 0 : ",I8)') msmooth
    write(*,*)
    stop
  end if
case('npmae')
  read(50,*,err=20) npmae0
case('wrtvars')
  read(50,*,err=20) wrtvars
case('ftmtype')
  read(50,*,err=20) ftmtype
case('tmomfix')
  write(*,*)
  write(*,'("Error(readinput): variable ''tmomfix'' is no longer used")')
  write(*,'(" use tm3fix instead")')
  write(*,*)
  stop
case('tm3fix')
  read(50,*,err=20) ntmfix
  if (ntmfix.le.0) then
    write(*,*)
    write(*,'("Error(readinput): ntmfix <= 0 : ",I8)') ntmfix
    write(*,*)
    stop
  end if
  if (allocated(itmfix)) deallocate(itmfix)
  allocate(itmfix(7,ntmfix))
  if (allocated(wkprfix)) deallocate(wkprfix)
  allocate(wkprfix(ntmfix))
  do i=1,ntmfix
    read(50,*,err=20) is,ia,l
    if ((is.le.0).or.(ia.le.0).or.(l.lt.0)) then
      write(*,*)
      write(*,'("Error(readinput): invalid is, ia or l in tm3fix block : ",&
       &4I8)') is,ia,l
      write(*,*)
      stop
    end if
    itmfix(1,i)=is
    itmfix(2,i)=ia
    itmfix(3,i)=l
! read k, p, r, t for the 3-index tensor
    read(50,*,err=20) itmfix(4:7,i)
! read 3-index tensor component with conventional normalisation
    read(50,*,err=20) wkprfix(i)
  end do
case('tauftm')
  read(50,*,err=20) tauftm
  if (tauftm.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): tauftm < 0 : ",G18.10)') tauftm
    write(*,*)
    stop
  end if
case('ftmstep')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''ftmstep'' is no longer used")')
case('cmagz','forcecmag')
  read(50,*,err=20) cmagz
case('rotavec')
  read(50,*,err=20) axang(:)
case('tstime')
  read(50,*,err=20) tstime
  if (tstime.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): tstime <= 0 : ",G18.10)') tstime
    write(*,*)
    stop
  end if
case('dtimes')
  read(50,*,err=20) dtimes
  if (dtimes.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): dtimes <= 0 : ",G18.10)') dtimes
    write(*,*)
    stop
  end if
case('pulse')
  read(50,*,err=20) npulse
  if (npulse.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): npulse < 1 : ",I8)') npulse
    write(*,*)
    stop
  end if
  if (allocated(pulse)) deallocate(pulse)
  allocate(pulse(12,npulse))
  do i=1,npulse
    read(50,'(A)',err=20) str
    str=trim(str)//' 1.0 0.0 0.0 0.0'
    read(str,*,err=20) pulse(:,i)
  end do
case('ramp')
  read(50,*,err=20) nramp
  if (nramp.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): nramp < 1 : ",I8)') nramp
    write(*,*)
    stop
  end if
  if (allocated(ramp)) deallocate(ramp)
  allocate(ramp(12,nramp))
  do i=1,nramp
    read(50,'(A)',err=20) str
    str=trim(str)//' 1.0 0.0 0.0 0.0'
    read(str,*,err=20) ramp(:,i)
  end do
case('ncgga')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''ncgga'' is no longer used")')
case('dncgga')
  read(50,*,err=20) dncgga
  if (dncgga.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): dncgga < 0 : ",G18.10)') dncgga
    write(*,*)
    stop
  end if
case('ntswrite')
  read(50,'(A)',err=20) str
  str=trim(str)//' 1'
  read(str,*,err=20) ntswrite(:)
case('nxoapwlo','nxapwlo')
  read(50,*,err=20) nxoapwlo
  if (nxoapwlo.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nxoapwlo < 0 : ",I8)') nxoapwlo
    write(*,*)
    stop
  end if
case('nxlo')
  read(50,*,err=20) nxlo
  if (nxlo.lt.0) then
    write(*,*)
    write(*,'("Error(readinput): nxlo < 0 : ",I8)') nxlo
    write(*,*)
    stop
  end if
case('tdrho1d')
  read(50,*,err=20) tdrho1d
case('tdrho2d')
  read(50,*,err=20) tdrho2d
case('tdrho3d')
  read(50,*,err=20) tdrho3d
case('tdmag1d')
  read(50,*,err=20) tdmag1d
case('tdmag2d')
  read(50,*,err=20) tdmag2d
case('tdmag3d')
  read(50,*,err=20) tdmag3d
case('tdjr1d','tdcd1d')
  read(50,*,err=20) tdjr1d
case('tdjr2d','tdcd2d')
  read(50,*,err=20) tdjr2d
case('tdjr3d','tdcd3d')
  read(50,*,err=20) tdjr3d
case('tddos')
  read(50,*,err=20) tddos
case('tdlsj')
  read(50,*,err=20) tdlsj
case('epseph')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''epseph'' is no longer used")')
case('rndevt0')
  read(50,*,err=20) rndevt0
case('sxcscf','ssxc','rstsf')
  read(50,*,err=20) sxcscf
case('tempk')
  read(50,*,err=20) tempk
  if (tempk.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): tempk <= 0 : ",G18.10)') tempk
    write(*,*)
    stop
  end if
! set Fermi-Dirac smearing
  stype=3
! set the smearing width
  swidth=kboltz*tempk
case('avecu')
  read(50,*,err=20) avecu(:,1)
  read(50,*,err=20) avecu(:,2)
  read(50,*,err=20) avecu(:,3)
case('scaleu')
  read(50,*,err=20) scu
case('scaleu1')
  read(50,*,err=20) scu1
case('scaleu2')
  read(50,*,err=20) scu2
case('scaleu3')
  read(50,*,err=20) scu3
case('q0cut')
  read(50,*,err=20) q0cut
  if (q0cut.lt.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): q0cut < 0 : ",G18.10)') q0cut
    write(*,*)
    stop
  end if
case('rndbfcu')
  read(50,*,err=20) rndbfcu
case('efielduc')
  read(50,*,err=20) efielduc
case('tplotq0')
  read(50,*,err=20) tplotq0
case('trdvclr')
  read(50,*,err=20) trdvclr
case('evtype')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''evtype'' is no longer used")')
case('wmaxgw')
  read(50,*,err=20) wmaxgw
case('twdiag')
  read(50,*,err=20) twdiag
case('tsediag')
  read(50,*,err=20) tsediag
case('actype')
  read(50,*,err=20) actype
case('npole')
  read(50,*,err=20) npole
  if (npole.le.0) then
    write(*,*)
    write(*,'("Error(readinput): npole <= 0 : ",I8)') npole
    write(*,*)
    stop
  end if
case('nspade')
  read(50,*,err=20) nspade
  if (nspade.le.0) then
    write(*,*)
    write(*,'("Error(readinput): nspade <= 0 : ",I8)') nspade
    write(*,*)
    stop
  end if
case('tfav0')
  read(50,*,err=20) tfav0
case('rmtall')
  read(50,*,err=20) rmtall
case('maxthd','omp_num_threads','OMP_NUM_THREADS')
  read(50,*,err=20) maxthd
case('maxthd1')
  read(50,*,err=20) maxthd1
case('maxthdmkl')
  read(50,*,err=20) maxthdmkl
case('maxlvl','omp_max_active_levels','OMP_MAX_ACTIVE_LEVELS')
  read(50,*,err=20) maxlvl
  if (maxlvl.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): maxlvl < 1 : ",I8)') maxlvl
    write(*,*)
    stop
  end if
case('stable')
  read(50,*,err=20) lv
  if (lv) then
    mixtype=3
    broydpm(1)=min(broydpm(1),0.04d0)
    broydpm(2)=min(broydpm(2),0.01d0)
    msmooth=max(msmooth,8)
    isgkmax=-4
    if (mp_mpi) then
      write(*,*)
      write(*,'("Info(readinput): parameters set by stable option")')
      write(*,'(" mixtype : ",I4)') mixtype
      write(*,'(" broydpm : ",2G18.10)') broydpm
      write(*,'(" msmooth : ",I4)') msmooth
      write(*,'(" isgkmax : ",I4)') isgkmax
    end if
  end if
case('metagga')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''metagga'' is no longer used")')
case('t0tdlr')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''t0tdlr'' is no longer used")')
case('tdphi')
  read(50,*,err=20) tdphi
! convert phase from degrees to radians
  tdphi=tdphi*pi/180.d0
case('thetamld')
  read(50,*,err=20) thetamld
! convert MLD angle from degrees to radians
  thetamld=thetamld*pi/180.d0
case('ntsbackup')
  read(50,*,err=20) ntsbackup
case('seedname')
  read(50,*,err=20) seedname
  seedname=adjustl(seedname)
case('num_wann')
  read(50,*,err=20) num_wann
case('idxw90','wann_bands')
  read(50,'(A)',err=20) str
  num_bands=1024
  if (allocated(idxw90)) deallocate(idxw90)
  allocate(idxw90(num_bands))
  call numlist(str,num_bands,idxw90)
case('num_iter')
  read(50,*,err=20) num_iter
case('trial_step')
  read(50,*,err=20) trial_step
case('xlwin','wannierExtra')
  if (allocated(xlwin)) deallocate(xlwin)
  allocate(xlwin(0))
  nxlwin=0
  do
    read(50,'(A)',err=20) str
    if (trim(str).eq.'') goto 10
    nxlwin=nxlwin+1
    call addstr(xlwin)
  end do
case('tbdip')
  read(50,*,err=20) tbdip
case('tjr','tcden')
  read(50,*,err=20) tjr
case('tauefm')
  read(50,*,err=20) tauefm
case('epsefm')
  read(50,*,err=20) epsefm
case('ehfb')
  read(50,*,err=20) ehfb
case('t0gclq0')
  read(50,*,err=20) t0gclq0
case('tafindt')
  read(50,*,err=20) tafindt
case('afindscf')
  read(50,*,err=20)
  write(*,*)
  write(*,'("Info(readinput): variable ''afindscf'' is no longer used")')
case('afindpm')
  read(50,*,err=20) afindpm(:)
case('nkspolar')
  read(50,*,err=20) nkspolar
  if (nkspolar.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): nkspolar < 1 : ",I8)') nkspolar
    write(*,*)
    stop
  end if
case('ntsforce')
  read(50,*,err=20) ntsforce
  if (ntsforce.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): ntsforce < 1 : ",I8)') ntsforce
    write(*,*)
    stop
  end if
case('wphcut')
  read(50,*,err=20) wphcut
  if (wphcut.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): wphcut <= 0 : ",G18.10)') wphcut
    write(*,*)
    stop
  end if
case('ephscf')
  read(50,*,err=20) ephscf(:)
case('anomalous')
  read(50,*,err=20) anomalous
case('tephde')
  read(50,*,err=20) tephde
case('bdiag')
  read(50,*,err=20) bdiag
case('ecutb')
  read(50,*,err=20) ecutb
  if (ecutb.le.0.d0) then
    write(*,*)
    write(*,'("Error(readinput): ecutb <= 0 : ",G18.10)') ecutb
    write(*,*)
    stop
  end if
case('ediag')
  read(50,*,err=20) ediag
case('pwxpsn')
  read(50,*,err=20) pwxpsn
  if (pwxpsn.lt.1) then
    write(*,*)
    write(*,'("Error(readinput): pwxpsn < 1 : ",I8)') pwxpsn
    write(*,*)
    stop
  end if
case('ramdisk')
  read(50,*,err=20) ramdisk
case('wrtdsk')
  read(50,*,err=20) wrtdsk
case('epsdmat')
  read(50,*,err=20) epsdmat
case('tm3old')
  read(50,*,err=20) tm3old
case('batch')
  read(50,*,err=20) batch
case('')
  goto 10
case default
  write(*,*)
  write(*,'("Error(readinput): invalid block name : ",A)') trim(block)
  write(*,*)
  stop
end select
goto 10
20 continue
write(*,*)
write(*,'("Error(readinput): error reading from elk.in")')
write(*,'("Problem occurred in ''",A,"'' block")') trim(block)
write(*,'("Check input convention in manual")')
write(*,*)
stop
30 continue
close(50)
! scale the speed of light
solsc=sol*solscf
! scale and rotate the lattice vectors (not referenced again in code)
avec(:,:)=sc*avec(:,:)
avec(:,1)=sc1*avec(:,1)
avec(:,2)=sc2*avec(:,2)
avec(:,3)=sc3*avec(:,3)
avec(1,:)=scx*avec(1,:)
avec(2,:)=scy*avec(2,:)
avec(3,:)=scz*avec(3,:)
t1=axang(4)
if (t1.ne.0.d0) then
  t1=t1*pi/180.d0
  call axangrot(axang(:),t1,rot)
  do i=1,3
    v(:)=avec(:,i)
    call r3mv(rot,v,avec(:,i))
  end do
end if
! randomise lattice vectors if required
if (rndavec.gt.0.d0) then
  do i=1,3
    do j=1,3
      t1=rndavec*(randomu()-0.5d0)
      avec(i,j)=avec(i,j)+t1
    end do
  end do
end if
! case of isolated molecule
if (molecule) then
! convert atomic positions from Cartesian to lattice coordinates
  call r3minv(avec,ainv)
  do is=1,nspecies
    do ia=1,natoms(is)
      call r3mv(ainv,atposl(:,ia,is),v)
      atposl(:,ia,is)=v(:)
    end do
  end do
end if
! randomise atomic positions if required
if (rndatposc.gt.0.d0) then
  call r3minv(avec,ainv)
  do is=1,nspecies
    do ia=1,natoms(is)
      call r3mv(avec,atposl(:,ia,is),v)
      do i=1,3
        t1=rndatposc*(randomu()-0.5d0)
        v(i)=v(i)+t1
      end do
      call r3mv(ainv,v,atposl(:,ia,is))
    end do
  end do
end if
! randomise the muffin-tin magnetic fields if required
if (rndbfcmt.gt.0.d0) then
  do is=1,nspecies
    do ia=1,natoms(is)
      do i=1,3
        t1=rndbfcmt*(randomu()-0.5d0)
        bfcmt0(i,ia,is)=bfcmt0(i,ia,is)+t1
      end do
    end do
  end do
end if
! set fxctype to fxctype if required
if (fxctype(1).eq.-1) fxctype(:)=xctype(:)
! find primitive cell if required
if (primcell) call findprimcell
! scale the ultracell vectors if required
avecu(:,1)=scu1*avecu(:,1)
avecu(:,2)=scu2*avecu(:,2)
avecu(:,3)=scu3*avecu(:,3)
avecu(:,:)=scu*avecu(:,:)
! read in atomic species data
call readspecies
return

contains

subroutine addstr(slist)
implicit none
! arguments
character(256), intent(inout), allocatable :: slist(:)
! local variables
integer n
! allocatable arrays
character(256), allocatable :: stmp(:)
n=size(slist)
allocate(stmp(n))
stmp(1:n)=slist(1:n)
deallocate(slist)
allocate(slist(n+1))
slist(1:n)=stmp(1:n)
slist(n+1)=str
deallocate(stmp)
end subroutine

end subroutine
!EOC

