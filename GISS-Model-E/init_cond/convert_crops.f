C**** convert_crops.f converts crops data from 720x360 to model grid
C****     and crop fraction per grid box to crop fraction/soil fraction
      parameter (ima=720,jma=360,offia=0.,divja=360.) !  input grid
!     parameter (imb=144,jmb=90 ,offib=0.,divjb=90. ) !  output grid
      parameter (skip=-9.999)

      character*80 :: title='  yr crops/soil (1) - Julia Pongratz et al.
     *, Hamburg 2007', titleT
      real*4 crops(ima,jma),wt(ima,jma)
!     real*4 cropso(imb,jmb),psoil(imb,jmb),fil(imb,jmb),ones(imb,jmb)
      real*4, allocatable :: cropso(:,:),psoil(:,:),fil(:,:),ones(:,:)
      character*150 line,filein,fileTopo,fileout

      if(iargc()<3) then
        write(0,*)'Usage: convert_crops IMxJM topo_file crop_file(s)'
        stop
      end if

!**** get model parameters (imb,jmb and related parameters)
      call getarg(1,titleT)
      n=index(titleT,'x') ; if (n<=0) n=index(titleT,'X')
      if (n<=0) then
         write(0,*) 'Unable to get model-dim from ',trim(titleT)
         stop
      end if
      read(titleT(1:n-1),*) imb ; read(titleT(n+1:n+5),*) jmb
      write(*,*) 'output grid:',imb,'x',jmb
      allocate (cropso(imb,jmb),psoil(imb,jmb))
      allocate (   fil(imb,jmb), ones(imb,jmb))
      if(jmb==90) divjb=90.
      if(jmb==46) divjb=45.
      offib=0.

!**** read in land fraction (not counting land ice)
      call getarg(2,fileTopo)
      open(3,file=fileTopo,form='unformatted',status='old')
      read (3) ; read(3) ; read (3) titleT,psoil ; close (3)

!**** prepare interpolation to model grid
      call hntrp0 (ima,jma,offia,divja,
     *             imb,jmb,offib,divjb, skip)

!**** output file(s)
   !! open(12,file='CROPS2007',form='unformatted')
      open(13,file='CROPS2007.ext',form='unformatted')

!**** Read in crops files (.5x.5 degree)
      do n=3,iargc()
        call getarg(n,filein)
        open(1,file=filein,form='formatted',status='old')

        crops=-10. !! may be used for testing completeness
        do
          read(1,'(a)',end=100) line
          call getinfo(line, lon, lat, frac , iyr)
          crops(lon,lat) = frac
        end do
  100   continue
        write(title(1:4),'(i4)') iyr
   !!   fileout='crops_720x360_yr_yyyy.pongratz.bin'
   !!   write(fileout(18:21),'(i4.4)') iyr
   !!   open(2,file=fileout,form='unformatted')
   !!   write(2) title,crops
   !!   close(2)
        close(1)

!**** create the wt-array
        do j=1,jma
        do i=1,ima
          wt(i,j)=1.
          if(crops(i,j) < 0.) wt(i,j)=0.
        end do
        end do

!**** interpolate to model grid
        call hntrpp(wt,crops,cropso)

!**** convert box fraction to soil fraction: crops->crops/soil
        do j=1,jmb
        do i=1,imb
          if(psoil(i,j).le.0..or.cropso(i,j) < 0.) then
           cropso(i,j)=skip
          else
           cropso(i,j)=min(1.,cropso(i,j)/psoil(i,j))
          end if
          ones(i,j)=1. ! prepare to extend data to all points
        end do
        end do

!**** Fill in missing data from nearest neighbors over SOIL only
        call fillin (fil, cropso,psoil, imb,jmb, skip, .true.)

!!*** write to disk data restricted to soil only (not done here)
   !!   write(12) title,fil

!**** Extend data to ocean and glaciers
        call fillin (cropso, fil,ones, imb,jmb, skip, .false.)
        write(13) title,cropso
      end do
      stop
      end

      subroutine getinfo(line, i ,j , f,iyr)
      character*150 line
      read(line(9:12),*) iyr

      nlat=index(line,'latitude[')+9
      n1=nlat+index(line(nlat:nlat+5),']')-2
      read (line(nlat:n1),*) j ; j=j+1       ! c->fortran

      nlon=index(line,'longitude[')+10
      n1=nlon+index(line(nlon:nlon+5),']')-2
      read (line(nlon:n1),*) i ; i=i+1       ! c->fortran

      nfrac=n1+10+index(line(n1+10:150),']')+1
      n1=nfrac+index(line(nfrac:150),' fraction')-2
      read (line(nfrac:n1),*) f
      return
      end subroutine getinfo

!!!!! include 'fillin.f'
      subroutine fillin (fil, dat,pe, im,jm, skip, vrb)
C**** nearest neighbor fill of 'skip's for all points where pe>0
      real fil(im,jm),dat(im,jm),pe(im,jm)
      logical vrb
      do j=1,jm
      do i=1,im
        fil(i,j)=dat(i,j)
        if(pe(i,j).gt.0..and.dat(i,j).eq.skip) then
          call getnn(im,jm,dat, skip, i,j, in,jn,nd)
          if(vrb) write(*,*) i,j,' filled from',in,jn,dat(in,jn),nd
          fil(i,j)=dat(in,jn)
        end if
      end do
      end do

      return
      end

      subroutine getnn(im,jm,dat, skip, i0,j0, in,jn,nd)
C**** dat(i0,j0)=skip; find (i,j) with dat(i,j).ne.skip 'closest'
C**** to point (i0,j0) - 'distance' = (j-j0)^2 + min(|i-i0|,im-|i-i0|)
C**** taking into account that lat is more important than lon: (in,jn)
      real dat(im,jm)

C****    the distance function
      ndist(ix,jx,iy,jy)=(jx-jy)**2+min(abs(ix-iy),im-abs(ix-iy))

      nd=jm*jm+im       ! larger than ndist for any pair of points
      do 20 id=1,im/2                   ! first try j=j0
      if(id.ge.nd) go to 40
      do 10 m=-1,1,2
      i=i0+m*id
      if(i.gt.im) i=i-im
      if(i.lt.1) i=i+im
      if (dat(i,j0).eq.skip) go to 10
        nd0=id  !  =ndist(i0,j0,i,j0)
        if (nd0.lt.nd) then
          in=i
          jn=j0
          nd=nd0
          go to 40
        end if
   10 continue
   20 continue

   40 do 100 jd=1,jm-1
      if (jd*jd.ge.nd) return
      do 60 m=-1,1,2
      j=j0+m*jd
      if(j.gt.jm) go to 60
      if(j.le.0) go to 60
      do 50 i=1,im
      if (dat(i,j).eq.skip) go to 50
        nd0=ndist(i0,j0,i,j)
        if (nd0.lt.nd) then
          in=i
          jn=j
          nd=nd0
        end if
   50 continue
   60 continue
  100 continue

      return
      end
!!!!! end of  fillin.f

      include 'HNTRPS.f'
