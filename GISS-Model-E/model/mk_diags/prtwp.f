      program prtwp
! Print wave power tables from acc-file arrays.
! Usage: prtwp nwin accfile(s)
!   where nwin is the number of desired time windows and accfile(s) are
!   one or more monthly acc-files containing the GCM diagnostic array "wave".
!   It is the responsibility of the user to specify the list of acc-files
!   in time-increasing order.
!
! Code adapted from online routine DIAG7P.  Much of the hard-coding in
! DIAG7P was eliminated for conciseness/readability, but a great deal
! of inflexibility remains.
! Output should be identical to DIAG7P for single months, but when
! months are combined, it differs from the old-I/O pdE
! This version combines months as follows:
! (1) timeseries from consecutive months are concatenated
! (2) the resulting longer sequence is evenly split into nwin time windows
! (3) subroutine MEM (Maximum Entropy Method) is applied to each window
! (4) MEM outputs from each time window are averaged over windows
! 
! The anticipated use of multiple windows is to aggregate seasonal results
! over multiple years, in which case nwin would be the number of years to
! be averaged.  A multi-year timeseries containing all months can
! be processed with nwin=1.  It is desirable to choose nwin so that
! no time window contains a time-axis discontinuity.
!
! The following checks are currently applied to the list of input files to
! facilitate combining months and to avoid unintended results:
! (1) If multiple years are present, each year should have the same number
!     of months
! (2) mod(number of files, nwin) == mod(number of timesteps, nwin) == 0
! (3) Each acc-file should be for a single month, not a sum over months.
! (4) Month gaps are not allowed, e.g. no January and March without Febrary.
      implicit none
      include 'netcdf.inc'
      integer :: status,acc_fid

      integer :: ifile,nfiles,iwin,nwin
      character(len=132) :: xlabel
      character(len=100) :: fromto
      character(len=256) :: accfile
      integer, dimension(12) :: monacc_part
      integer, dimension(0:13) :: monacc

      real*4, dimension(:,:,:,:), allocatable :: wave
      real*4, dimension(:,:,:,:,:), allocatable :: wave_series
      integer, parameter :: mmax=12,nuamax=120,nubmax=15,nuxmax=41
      integer, dimension(nuxmax) :: nu1_,nu2_
      real*8, dimension(nuamax) :: power

      real*8, dimension(:,:,:), allocatable :: xpower
      real*8, dimension(:,:), allocatable :: fpe,pnu,var

      real*8 :: xpower_(nuxmax),fpe_(mmax+1),pnu_,var_

!@var comp_wave complex form of wave. correct arg. to subr. mem
      complex*16, dimension(:), allocatable :: comp_wave
      integer :: kq,n,nmax,nu1,nu2,nux,nuxo,kwp,nwav_dag,
     &     max12hr_sequ,itime,ntime,ntime_win,ic,ic0
      character(len=132), dimension(:), allocatable :: title
      real*8, dimension(:), allocatable :: scalet
      character(len=4) :: cnwin

      nfiles = iargc()-1
      if(nfiles.le.0) then
        write(6,*) 'usage: prtwp nwin accfile(s)'
        stop
      endif

      cnwin = ''
      call getarg(1,cnwin)
      read(cnwin,*) nwin

      if(mod(nfiles,nwin).ne.0) stop 'mod(nfiles,nwin) not zero'

c
c read the list of accfiles
c
      monacc = 0
      do ifile=1,nfiles

      call getarg(ifile+1,accfile)

      call handle_err(nf_open(accfile,nf_nowrite,acc_fid),
     &     'opening '//trim(accfile))

c
c get run ID info from the first file
c
      if(ifile.eq.1) then
        xlabel=''; fromto=''
        status = nf_get_att_text(acc_fid,nf_global,'xlabel',xlabel)
        status = nf_get_att_text(acc_fid,nf_global,'fromto',fromto)
      endif

      call get_dimsize(acc_fid,'kwp',kwp)
      if(kwp.ne.12) stop 'bad kwp: expecting 12'

      call get_dimsize(acc_fid,'max12hr_sequ',max12hr_sequ)
      call get_dimsize(acc_fid,'nwav_dag',nwav_dag)

      call get_var_int(acc_fid,'monacc',monacc_part)
      if(sum(monacc_part).gt.1)
     &     stop 'do not use summed accfiles with prtwp'
      monacc(1:12) = monacc(1:12) + monacc_part

      if(ifile.eq.1) then
        allocate(wave_series(2,max12hr_sequ,nwav_dag,kwp,nfiles))
        allocate(wave(2,max12hr_sequ,nwav_dag,kwp))
      endif

      call get_var_real(acc_fid,'wave',wave)
      wave_series(:,:,:,:,ifile) = wave

      status = nf_close(acc_fid)

      enddo ! end loop over files

      monacc(0) = monacc(12)
      monacc(13) = monacc(1)
! other sanity checks on input sequence
      if(minval(monacc,mask=monacc.gt.0).ne.
     &   maxval(monacc,mask=monacc.gt.0)) stop
     &     'unequal numbers of months across years'
      do n=1,12
        if(monacc(n).eq.0 .and. monacc(n-1).gt.0 .and.
     &       monacc(n+1).gt.0) stop 'monacc gap'
      enddo

! concatenate timesteps
      deallocate(wave)
      allocate(wave(2,nfiles*max12hr_sequ,nwav_dag,kwp))
      ntime = 0
      do ifile=1,nfiles
        do itime=1,max12hr_sequ
          if(all(wave_series(:,itime,:,:,ifile).eq.0.)) exit
          ntime = ntime + 1
          wave(:,ntime,:,:) = wave_series(:,itime,:,:,ifile)
        enddo
      enddo
      if(ntime.le.mmax) stop 'ntime < mmax'

      allocate(title(kwp),scalet(kwp))

      kq = 0

      kq = kq + 1
      scalet(kq) = 1.
      title(kq) =
     &     'WAVE POWER FOR U NEAR 850 MB AND EQUATOR (DAY*(m/s)^2)'

      kq = kq + 1
      scalet(kq) = 1.
      title(kq) =
     &     'WAVE POWER FOR V NEAR 850 MB AND EQUATOR (DAY*(m/s)^2)'

      kq = kq + 1
      scalet(kq) = .1d0
      title(kq) =
     &     'WAVE POWER FOR U NEAR 300 MB AND EQUATOR (10 DAY*(m/s)^2)'

      kq = kq + 1
      scalet(kq) = 1.
      title(kq) =
     &     'WAVE POWER FOR V NEAR 300 MB AND EQUATOR (DAY*(m/s)^2)'

      kq = kq + 1
      scalet(kq) = .1d0
      title(kq) =
     &     'WAVE POWER FOR U NEAR 50 MB AND EQUATOR (10 DAY*(m/s)^2)'

      kq = kq + 1
      scalet(kq) = 1d0
      title(kq) =
     &     'WAVE POWER FOR V NEAR 50 MB AND EQUATOR (DAY*(m/s)^2)'

      kq = kq + 1
      scalet(kq) = 1d-3
      title(kq) =
     &     'WAVE POWER FOR PHI AT 922 MB AND 50 DEG N. (10**3 DAY*m^2)'

      kq = kq + 1
      scalet(kq) = 1d-3
      title(kq) =
     &     'WAVE POWER FOR PHI AT 700 MB AND 50 DEG N. (10**3 DAY*m^2)'

      kq = kq + 1
      scalet(kq) = 1d-3
      title(kq) =
     &     'WAVE POWER FOR PHI AT 500 MB AND 50 DEG N. (10**3 DAY*m^2)'

      kq = kq + 1
      scalet(kq) = 1d-3
      title(kq) =
     &     'WAVE POWER FOR PHI AT 300 MB AND 50 DEG N. (10**3 DAY*m^2)'

      kq = kq + 1
      scalet(kq) = 1d-4
      title(kq) =
     &     'WAVE POWER FOR PHI AT 100 MB AND 50 DEG N. (10**4 DAY*m^2)'

      kq = kq + 1
      scalet(kq) = 1d-5
      title(kq) =
     &     'WAVE POWER FOR PHI AT 10 MB AND 50 DEG N. (10**5 DAY*m^2)'


      if(mod(ntime,nwin).ne.0) stop 'mod(ntime,nwin) not zero'

      ntime_win = ntime/nwin
      allocate(comp_wave(ntime_win))

      allocate(xpower(nuxmax,nwav_dag,kwp))
      allocate(fpe(mmax+1,kwp))
      allocate(pnu(nwav_dag,kwp))
      allocate(var(nwav_dag,kwp))

      ! set up wavenumber-bounds lists for subsequent averaging over wavenumbers
      nu2_(1:4) = (/27, 34, 38, 40 /)
      do nux=5,40
        nu2_(nux) = 36+nux
      enddo
      nu2_(nuxmax) = nuamax
      nu1_(1) = 2
      do nux=2,nuxmax
        nu1_(nux) = 1 + nu2_(nux-1)
      enddo

      ! loop over time windows and sum up contributions from each
      xpower = 0.
      fpe = 0.
      var = 0.
      pnu = 0.

      do iwin=1,nwin
      ic0 = (iwin-1)*ntime_win
      do kq=1,kwp
        do n=nwav_dag,1,-1
          do ic=1,ntime_win
            comp_wave(ic) =
     &           cmplx( wave(1,ic0+ic,n,kq) , wave(2,ic0+ic,n,kq) )
          enddo
          call mem(comp_wave,ntime_win,mmax,nuamax,nubmax,
     &         power,fpe_,var_,pnu_)
          do nux=1,nuxmax
            if(kq.le.6) then    ! Equator
              nuxo = nux
              nu1 = nu1_(nux)
              nu2 = nu2_(nux)
            else                ! 50 N
              nuxo = nuxmax-nux+1
              nu1 = nuamax-nu2_(nux)+2
              nu2 = nuamax-nu1_(nux)+2
            endif
            if(nux.eq.1 .or. nux.eq.nuxmax) then
              xpower_(nuxo) = (.5d0*power(1) + sum(power(nu1:nu2)))/
     &             (.5d0 + real(nu2-nu1+1,kind=8))
            else
              xpower_(nuxo) =
     &             sum(power(nu1:nu2))/real(nu2-nu1+1,kind=8)
            endif
          enddo

          xpower(:,n,kq) = xpower(:,n,kq) + xpower_
          var(n,kq) = var(n,kq) + var_
          pnu(n,kq) = pnu(n,kq) + pnu_

        enddo ! n
        fpe(:,kq) = fpe(:,kq) + fpe_
      enddo ! kq
      enddo ! iwin

      xpower = xpower/nwin
      fpe = fpe/nwin
      var = var/nwin
      pnu = pnu/nwin


      ! Print the averages over all time windows

      do kq=1,kwp
        if(mod(kq-1,3).eq.0) then
          write(6,'(a105)') xlabel(1:105)
c          write(6,'(a57)') fromto(1:57) ! todo: properly combine calendar info
        endif
        if(kq.le.6) then        ! Equator
          write(6,901) trim(title(kq))
        else                    ! 50 N
          write(6,911) trim(title(kq))
        endif
        do n=nwav_dag,1,-1
          write(6,902) n,
     &         int(scalet(kq)*xpower(:,n,kq)+.5d0),
     &         int(10.*scalet(kq)*var(n,kq)+.5d0),
     &         int(1000.*scalet(kq)*(var(n,kq)-pnu(n,kq))+.5d0)
        enddo ! n
        if(kq.gt.6) fpe(:,kq) = 1000.*scalet(kq)*fpe(:,kq) ! why not for eq?
        write(6,903) fpe(:,kq)
      enddo ! kq

C****
  901 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',
     *   35('---')/' N    -2      *-3   -3.3      -4       -5    -6   -7
     *.5  -10-12-15-20-30-60    60 30 20 15 12 10    7.5    6     5
     *   4*   VAR ERR'/'   --',40('---'))
  902 FORMAT (I2,41I3,I4,I4)
  903 FORMAT ('   --',40('---')/(1X,13F10.4))
  907 FORMAT ('1',A,I3,1X,A3,I5,' - ',I3,1X,A3,I5)
  911 FORMAT ('0',30X,A64,8X,'*1/60 (1/DAY)'/'   PERIOD EASTWARD--',
     *  35('---')/               ' N   *-4       -5    -6   -7.5  -10-12
     *-15-20-30-60    60 30 20 15 12 10    7.5    6     5        4
     * 3.3    3*       2    VAR ERR'/'   --',40('---'))

      end program prtwp


      SUBROUTINE MEM (SERIES,ITM,MMAX,NUAMAX,NUBMAX,POWER,FPE,VAR,PNU)
      IMPLICIT NONE
      DIMENSION C(1800),S(1800),A(12),AA(11),P(13)
      DIMENSION SERIES(ITM),POWER(NUAMAX),FPE(MMAX+1)
      REAL*8 ARG,PP,POWERX,P,C,S,POWER,FPE
      COMPLEX*16 CI,CSUM,A,AA,ANOM,ADEN
      COMPLEX*16 SERIES
      REAL*8 :: PNU,VAR
      COMPLEX*16, DIMENSION(:), ALLOCATABLE :: B1,B2
      INTEGER ::
     &     I,ITM,L,M,MMAX,MMAXP1,NU,NUA,
     &     NUAMAX,NUB,NUBMAX,NUMAX,NUTM
      REAL*8, PARAMETER :: PI=ACOS(-1D0)

      ALLOCATE(B1(ITM),B2(ITM))

      CI=CMPLX(0.D0,1.D0)
      MMAXP1=MMAX+1
C**COSINE AND SINE FUNCTION
      NUMAX=NUAMAX*NUBMAX
      DO 20 NU=1,NUMAX
      ARG=2.0*PI*FLOAT(NU)/FLOAT(NUMAX)
      C(NU)=DCOS(ARG)
   20 S(NU)=DSIN(ARG)
   50 PP=0.0
      DO 60 I=1,ITM
   60 PP=PP+SERIES(I)*CONJG(SERIES(I))
      P(1)=PP/FLOAT(ITM)
      VAR=P(1)
      M=1
      B1(1)=SERIES(1)
      B2(ITM-1)=SERIES(ITM)
      DO 70 I=2,ITM-1
      B1(I)=SERIES(I)
   70 B2(I-1)=SERIES(I)
      GO TO 80
  100 DO 110 I=1,M
  110 AA(I)=A(I)
      M=M+1
      DO 120 I=1,ITM-M
      B1(I)=B1(I)-CONJG(AA(M-1))*B2(I)
  120 B2(I)=B2(I+1)-AA(M-1)*B1(I+1)
   80 ANOM=CMPLX(0.D0,0.D0)
      ADEN=CMPLX(0.D0,0.D0)
      DO 90 I=1,ITM-M
      ANOM=ANOM+CONJG(B1(I))*B2(I)
   90 ADEN=ADEN+B1(I)*CONJG(B1(I))+B2(I)*CONJG(B2(I))
      A(M)=(ANOM+ANOM)/ADEN
      P(M+1)=P(M)*(1.0-CONJG(A(M))*A(M))
      IF (M.EQ.1) GO TO 100
  130 CONTINUE
      DO 140 I=1,M-1
  140 A(I)=AA(I)-A(M)*CONJG(AA(M-I))
      IF (M.LT.MMAX) GO TO 100
C**FINAL PREDICTION ERROR
      DO 150 M=1,MMAXP1
  150 FPE(M)=P(M)*FLOAT(ITM+M-1)/FLOAT(ITM-M+1)
      DO 180 NUA=1,NUAMAX
      POWERX=0.
C**FREQUENCY BAND AVERAGE
      DO 170 NUB=1,NUBMAX
      NU=NUB+NUA*NUBMAX+(NUMAX-3*NUBMAX-1)/2
      CSUM=1.
      DO 160 M=1,MMAX
      NUTM=MOD(NU*M-1,NUMAX)+1
  160 CSUM=CSUM-A(M)*(C(NUTM)-CI*S(NUTM))
  170 POWERX=POWERX+P(MMAXP1)/(CSUM*CONJG(CSUM))
      POWER(NUA)=.5*POWERX/FLOAT(NUBMAX)
  180 CONTINUE
      PNU=0.0
      DO 210 L=1,NUAMAX
  210 PNU=PNU+POWER(L)
      PNU=PNU/(.5*NUAMAX)
      RETURN
      END SUBROUTINE MEM
