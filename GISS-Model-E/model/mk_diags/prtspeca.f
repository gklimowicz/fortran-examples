      module dft_mod
!@sum dft_mod quick implemention of DFT for standing eddy KE diagnostic
      integer, private :: im
      real, private :: byim
      real*4, dimension(:,:), allocatable, private :: sini,cosi
      contains
      subroutine dft_init(im_in)
      implicit none
      integer :: im_in
      real :: twopibyN,rk
      real, allocatable :: twopiri(:)
      integer :: i,k
      im = im_in
      byim = 1./real(im,kind=4)
      twopibyN = 2.*acos(-1.)*byim
      allocate(sini(im,0:im/2),cosi(im,0:im/2))
      allocate(twopiri(im))
      do i=1,im
        twopiri(i) = twopibyN*real(i,kind=4)
      enddo
      do k=0,im/2
        rk = real(k,kind=4)
        sini(:,k) = sin(rk*twopiri)
        cosi(:,k) = cos(rk*twopiri)
      enddo
      end subroutine dft_init
      subroutine dfte(x,e)
      implicit none
      real*4 :: x(im),e(0:im/2)
      real*4 :: ck,sk
      integer :: k
      do k=0,im/2
        ck = sum(x*cosi(:,k))
        sk = sum(x*sini(:,k))
        e(k) = (ck*ck + sk*sk)*byim
      enddo
      e(0) = e(0)/2.
      e(im/2) = e(im/2)/2.
      end subroutine dfte
      end module dft_mod

      program prtspeca
! Print spectral analysis tables.  Note the speca array has an extra dimension
! relative to the instance in the GCM (for readability)
      use dft_mod
      implicit none
      include 'netcdf.inc'
      integer :: status,acc_fid,ijk_fid

      integer, parameter :: njspeca=4,nlspeca_max=4 ! params for now
      character*8 :: latitd(njspeca) = (/
     &     'SOUTHERN','NORTHERN',' EQUATOR','45 NORTH'/)
      character*16 :: sphere(nlspeca_max)=
     &     (/'TROPOSPHERE     ','LOW STRATOSPHERE',
     &       'MID STRATOSPHERE','UPP STRATOSPHERE'/)

      real*4, dimension(nlspeca_max) :: scalel=(/1.,1.,10.,10./)
      real*4, dimension(njspeca) :: scalej=(/1.,1.,10.,10./)

      integer :: im,jm,lm,nm,nlspeca,nspher,kspeca

      integer :: idacc(12)

      real*4, dimension(:,:,:,:), allocatable :: speca
      integer, parameter :: ktpe=8,nhemi=2 ! params for now
      real*4, dimension(ktpe,nhemi) :: atpe

      integer :: nargs
      integer :: jspeca,lspeca,i,j,l,m,n,jeq

      real*4 :: pi,radian,radius

      real*4, parameter :: grav = 9.80665

      character(len=132) :: xlabel
      character(len=100) :: fromto

      logical :: have_aijk
      integer, dimension(njspeca) :: jmin_speca,jmax_speca
      integer, dimension(nlspeca_max) :: lmin_speca,lmax_speca
      real*4 :: dpav,fac
      real*4, dimension(:), allocatable :: x,upow,vpow,dxyv,lats,dxyp
      real*4, dimension(:,:,:), allocatable :: u,v,dp,uvpow,seke

      character(len=80) :: accfile,ijkfile,run_name,time_period,fmtstr,
     &     fmtbase

      character(len=30) :: cunitjw

      nargs = iargc()
      if(nargs.ne.2) then
        write(6,*) 'usage: prtspeca run_name time_period'
        write(6,*) '  e.g. prtspeca Exyz JAN1901'
        write(6,*) '       prtspeca Exyz ANN1901-1910'
        stop
      endif

      call getarg(1,run_name)
      call getarg(2,time_period)

      accfile=trim(time_period)//'.acc'//trim(run_name)//'.nc'
      ijkfile=trim(time_period)//'.aijk'//trim(run_name)//'.nc'

      write(6,*) 'Searching for the following input files:'
      write(6,*) '   ',trim(accfile)
      write(6,*) '   ',trim(ijkfile),' (optional)'

c
c open input files
c
      call handle_err(nf_open(accfile,nf_nowrite,acc_fid),
     &     'opening '//trim(accfile))
      status = nf_open(ijkfile,nf_nowrite,ijk_fid)
      have_aijk = status==nf_noerr
      if(.not.have_aijk) write(6,*)
     &     'no aijk-file: skipping standing eddy calculations'

      call get_dimsize(acc_fid,'imlonh_plus_1',nm)
      call get_dimsize(acc_fid,'kspeca',kspeca)
      call get_dimsize(acc_fid,'nlspeca',nlspeca)

      call get_dimsize(acc_fid,'nspher',nspher)
      if(njspeca*nlspeca.ne.nspher) stop 'bad nspher'

      call get_var_int(acc_fid,'idacc',idacc)

      allocate(speca(nm,kspeca,njspeca,nlspeca))
      call get_var_real(acc_fid,'speca',speca)
      call get_var_real(acc_fid,'atpe',atpe)
      speca = speca/idacc(1)
      atpe = atpe/idacc(1)

      if(have_aijk) then
C****
C**** spectral decomposition of the KE of the time-mean flow
C****

        call get_dimsize(ijk_fid,'lon2',im)
        call get_dimsize(ijk_fid,'lat2',jm)
        call get_dimsize(ijk_fid,'plm',lm)

        if(nm .ne. 1+im/2) stop 'mismatched im,nm'

        allocate(x(im),upow(nm),vpow(nm),uvpow(nm,jm,lm))
        allocate(u(im,jm,lm),v(im,jm,lm),dp(im,jm,lm))
        allocate(seke(nm,njspeca,nlspeca))
        allocate(dxyv(jm),dxyp(jm),lats(jm))

        pi = acos(-1.)
        radian = pi/180.
        radius = 6371000.

        call get_var_real(ijk_fid,'lat2',lats(1))
        dxyp(1) = 2.*pi*radius*radius*(sin(radian*lats(2))+1.)
        do j=2,jm-1
          dxyp(j) = 2.*pi*radius*radius*
     &         (sin(radian*lats(j+1))-sin(radian*lats(j)))
        enddo
        dxyp(jm) = dxyp(1)
        do j=2,jm
          dxyv(j) = .5*(dxyp(j-1)+dxyp(j))/real(im,kind=4)
        enddo

        call get_var_real(ijk_fid,'dpb',dp)
        call get_var_real(ijk_fid,'ub',u)
        call get_var_real(ijk_fid,'vb',v)

        call get_var_int(acc_fid,'lmax_speca',lmax_speca)
        lmin_speca(1) = 1
        do lspeca=2,nlspeca
          lmin_speca(lspeca) = 1 + lmax_speca(lspeca-1)
        enddo

        jeq = 1 + jm/2

        jmin_speca(1) = 2
        jmax_speca(1) = jeq

        jmin_speca(2) = jeq
        jmax_speca(2) = jm

        jmin_speca(3) = jeq
        jmax_speca(3) = jmin_speca(3)

        jmin_speca(4) = 2.+.75*(jm-1.)
        jmax_speca(4) = jmin_speca(4)

        call dft_init(im)

        uvpow = 0.
        do l=1,lm
        do j=2,jm
          dpav = sum(dp(:,j,l))/real(im,kind=8)
          if(dpav.le.0.) cycle
          x(:) = u(:,j,l)*dp(:,j,l)
          call dfte(x,upow)
          x(:) = v(:,j,l)*dp(:,j,l)
          call dfte(x,vpow)
          uvpow(:,j,l) = (upow + vpow)/dpav
        enddo
        enddo

        seke = 0.
        do lspeca=1,nlspeca
        do l=lmin_speca(lspeca),lmax_speca(lspeca)
          do jspeca=1,njspeca
          do j=jmin_speca(jspeca),jmax_speca(jspeca)
            if(j.eq.jeq .and.
     &           jmin_speca(jspeca).ne.jmax_speca(jspeca)) then
              fac = .5*dxyv(j)
            else
              fac = 1.*dxyv(j)
            endif
            seke(:,jspeca,lspeca) = seke(:,jspeca,lspeca)
     &         + uvpow(:,j,l)*fac
          enddo
          enddo
        enddo
        enddo
        speca(:,1,:,:) = seke*100.e-17 / grav
      else
        speca(:,1,:,:) = 0.
      endif

c
c get run ID, time/date info, number of latitudes and levels
c
      xlabel=''; fromto=''
      status = nf_get_att_text(acc_fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(acc_fid,nf_global,'fromto',fromto)

      do lspeca=1,nlspeca
      do jspeca=1,njspeca
        speca(:,:,jspeca,lspeca) =
     &  speca(:,:,jspeca,lspeca)*scalej(jspeca)*scalel(lspeca)
      enddo
      enddo

      fmtbase='i5,i6,i8,4i6,i8,i6,i7,2i6,i7,i6,i7,2i6,i8,i6)'

      do jspeca=1,njspeca ! loop over latitude zones

      if(jspeca.ge.3) then ! eq and 45 n are infinitesimal areas.  change power of 10
        cunitjw = '10**16 JOULES AND 10**11 WATTS'
      else
        cunitjw = '10**17 JOULES AND 10**12 WATTS'
      endif

      write(6,'(a)') '1'//xlabel
      write(6,'(a)') '0**  Spectral Analysis **      '//
     &     fromto(1:57)//'       UNITS '//cunitjw

      do lspeca=1,nlspeca !one for each level (trp/lstr/mstr/ustr)
        if (jm.ge.25.and.lspeca.eq.2) write (6,'(a1)') '1'
        write (6,'(a)')
     &       '0                                                  '//
     &       latitd(jspeca)//' '//sphere(lspeca)
        write(6,'(a)')
     &     '             MEAN                   '//
     &     'DYNAMICS                         SOURCES'//
     &     '                FILTER        DAILY    PR SURF     LAST'
        write(6,'(a)')
     &     '   N    SKE   KE   APE    KADV  KCOR   P-K  KDYN  PDYN   '//
     &     'KCNDS PCNDS   PRAD KSURF PSURF   KFIL  PFIL   KGMP  PGMP'//
     &     '    KE      KE   APE'
C**** Write KE and APE for each wavenumber
        fmtstr = '(a4,i7,'//fmtbase
        write (6,fmtstr) '0  0',nint(speca(1,1:kspeca,jspeca,lspeca))
        write(6,'(/)')
        fmtstr = '(i4,i7,'//fmtbase
        do n=2,nm
          write (6,fmtstr) n-1,nint(speca(n,1:kspeca,jspeca,lspeca))
        enddo
        fmtstr = '(a5,i6,'//fmtbase
        write(6,fmtstr) ' EDDY',
     &       nint(sum(speca(2:nm,:,jspeca,lspeca),dim=1))
        write(6,fmtstr) '0TOTL',
     &       nint(sum(speca(1:nm,:,jspeca,lspeca),dim=1))
      enddo ! end loop over altitude zones

      if(jspeca.lt.3) then
C**** total potential energy for each hemisphere
        fmtstr='(/a4,i18,i32,i14,i7,i12,2i13,i20)'
        write(6,fmtstr) '0TPE',nint(atpe(:,jspeca))
      endif

      enddo ! end loop over latitude zones

      end program prtspeca
