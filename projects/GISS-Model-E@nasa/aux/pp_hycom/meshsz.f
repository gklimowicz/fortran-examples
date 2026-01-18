      subroutine meshsz
c
      use const_proc, only: path0,latlonij,lp,iorign,jorign
      use hycom_dimen, only: idm,jdm,i,j
      use hycom_arrays, only: scp2,scux,scvx,scpx,scuy,scvy,scpy
     .    ,latij,lonij
c --- compute mesh size at u,v,p,q points
      implicit none
c
      real*4 lat4(idm,jdm,4),lon4(idm,jdm,4)
      real sphdis,hilat
      external sphdis
      integer iz,jz,ipol,jpol,ja,jb,n
c
c --- read in array of lat/lon values
      write (lp,'(a/a)') 'subr. meshsz -- about to open lat/lon file',
     .   latlonij
      open (unit=16,file=trim(path0)//latlonij,status='old',
     .      form='unformatted',action='read')
      read (16) iz,jz
      if (iz.ne.idm .or. jz.ne.jdm) then
        write (lp,'(2(a,2i5))') 'error - idm,jdm =',iz,jz,
     .   ' in lat/lon file should be',idm,jdm
        stop '(error)'
      end if
      rewind 16
      read (16) iz,jz,lat4,lon4
      close (unit=16)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- optional: locate north pole in grid
      hilat=-99.
      do 1 j=jz,1,-1
      do 1 i=1,iz
      if (lat4(i,j,4).gt.hilat) then
        hilat=lat4(i,j,4)
        ipol=i
        jpol=j
      end if
 1    continue
      if (hilat.gt.89.99) then
        write (*,'(a,2i5)') 'north pole at -q- point',ipol,jpol
      else
        write (*,*) 'unable to locate north pole'
        ipol=0
        jpol=0
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do n=1,4
        call extrct(lat4(1,1,n),iorign,jorign,latij(1,1,n))
        call extrct(lon4(1,1,n),iorign,jorign,lonij(1,1,n))
      end do
c
      do 56 j=1,jdm
      ja=mod(j-2+jdm,jdm)+1
      jb=mod(j      ,jdm)+1
      do 56 i=1,idm
      if (i.lt.idm) then
        scpx(i,j)=sphdis(latij(i  ,j,1),lonij(i  ,j,1),
     .                   latij(i+1,j,1),lonij(i+1,j,1))
        scvx(i,j)=sphdis(latij(i  ,j,4),lonij(i  ,j,4),
     .                   latij(i+1,j,4),lonij(i+1,j,4))
      else
        scpx(i,j)=scpx(i-1,j)
        scvx(i,j)=scvx(i-1,j)
      endif
c
      if (i.gt.1) then
        scux(i,j)=sphdis(latij(i  ,j,3),lonij(i  ,j,3),
     .                   latij(i-1,j,3),lonij(i-1,j,3))
      else
        scux(i,j)=scpx(i,j)
      endif
c
        scuy(i,j)=sphdis(latij(i,j  ,4),lonij(i,j  ,4),
     .                   latij(i,jb ,4),lonij(i,jb ,4))
        scpy(i,j)=sphdis(latij(i,j  ,2),lonij(i,j  ,2),
     .                   latij(i,jb ,2),lonij(i,jb ,2))
        scvy(i,j)=sphdis(latij(i,j  ,3),lonij(i,j  ,3),
     .                   latij(i,ja ,3),lonij(i,ja ,3))
      scp2(i,j)=scpx(i,j)*scpy(i,j)
 56   continue
c
      return
      end
c
      real function sphdis(x1,y1,x2,y2)
c --- dist.(m) between 2 points on sphere, lat/lon (x1,y1) and lat/lon (x2,y2)
      implicit none
      real x1,y1,x2,y2,ang,radius,radian
      data radius/6376.e3/,radian/57.2957795/
c
      ang=mod(y2-y1+540.,360.)-180.
      sphdis=radius*acos(min(1.,cosd(90.-x1)*cosd(90.-x2)
     .                         +sind(90.-x1)*sind(90.-x2)*cosd(ang)))
      if (sphdis.eq.0.) 
     .  sphdis=radius*sqrt((x2-x1)**2+(ang*cosd(.5*(x1+x2)))**2)/radian
c     if (sphdis.eq.0.) write (*,'(a,2f8.3,2x,2f8.3)')
c    .  'warning - zero distance between lat/lon points',x1,y1,x2,y2
      sphdis=max(sphdis,1.)
      return
      end
