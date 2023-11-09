#include "rundeck_opts.h"

!global ?
      subroutine geopar(iniOCEAN)
c
c --- set up model parameters related to geography
c
c --- hycom version 0.9 -- cyclic in j
css   USE GEOM, only : dxyp
c
      USE CONSTANT, only: radian  ! radian=pi/180.
      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT,broadcast
cddd      USE HYCOM_DIM_GLOB, only : ii,jj,kk,ii1,isp,ifp,ilp,ip,isq,ifq,ilq
cddd     &     ,isu,ifu,ilu,jsv,jfv,jlv,ntrcr,jsp,jfp,jlp,msk,iio,jjo
cddd     &     ,iia,jja,idm,jdm, iu,iv,iq
      USE HYCOM_DIM_GLOB
      USE HYCOM_SCALARS, only : pi,area,avgbot,huge,flnmlat,flnmdep
     &   ,flnmbas,ipacn,ipacs,jpac,iatln,iatls,jatl,beropn,ocnvol
     &   ,init_pr1d,zonarea
      USE HYCOM_ARRAYS_GLOB
      USE KPRF_ARRAYS
      USE HYCOM_CPLER
      USE HYCOM_DYNSI_CPLER
      use filemanager, only : findunit
      use hycom_dim, only : ogrid

      implicit none
      integer i,j,k,l,n,nn,ia,ib,ja,jb,jp,iu1,iu2,iu3
c
      logical, intent(in) :: iniOCEAN
      real realat,sphdis,sphrec,q,loncor(4),latcor(4)
      integer idim,jdim,length,iz,jz,nt
      character util(idm*jdm+14)*2,preambl(5)*79
      real*4 real4(idm,jdm),lat4(idm,jdm,4),lon4(idm,jdm,4)
c --- 'glufac' = regional viscosity enhancement factor
      real, parameter :: glufac=3., zero=0.
      logical:: succes
c
c --- read basin depth array
      if (AM_I_ROOT())
     .write (*,'(2a)') ' reading bathymetry file from ',flnmdep
      call findunit(iu1)
      open (unit=iu1,file=flnmdep,form='unformatted',status='old'
     .     ,convert='big_endian')
      read (iu1) iz,jz,real4
      close (unit=iu1)

      if (iz.ne.idm .or. jz.ne.jdm) then
        write (*,'(2(a,2i5))') 'depth file dimensions',iz,jz,
     .   '  should be',idm,jdm
        stop '(geopar)'
      end if
      do 9 j=1,jj
      do 9 i=1,ii
 9    depths(i,j)=real4(i,j)
c
c     do 7 j=1,jj
c     do 7 i=1,ii
c7    if (depths(i,j).gt.0.) depths(i,j)=max(botmin.,depths(i,j))
c
c --- reset the Denmark Strait - done in advance
c
c     depths(69,168)=798.
c     depths(69,169)=798.
c     depths(70,167)=798.
c     depths(70,168)=798.
c     depths(71,167)=798.
c
c --- reset the Iceland-Faeroes ridge - done in advance
c     do i=72,74
c     depths(i,175)=798.
c     end do
c
      !write(0,*) "ok ",__FILE__,__LINE__
      if (AM_I_ROOT()) then ! print only on root
        write (*,*) 'shown below: bottom depth'
        call zebra(depths,idm,ii1,jj)
      endif

      !write(0,*) "ok ",__FILE__,__LINE__
c
c --- determine do-loop limits for u,v,p,q points
      call bigrid(depths)

      !write(0,*) "ok ",__FILE__,__LINE__

!! copy hycom_dim arrays to global grid
      call gather_hycom_dim

      if (AM_I_ROOT()) then
ccc      do 3 i=1,ii1
ccc 3    write (*,'('' i='',i3,'' jfp,jlp='',7(1x,2i5))') i,
ccc     . (jfp(i,l),jlp(i,l),l=1,jsp(i))
ccc      do 5 j=1,jj
ccc 5    write (*,'('' j='',i3,'' ifp,ilp='',7(1x,2i5))') j,
ccc     . (ifp(j,l),ilp(j,l),l=1,isp(j))
c
c --- smooth bottom topography (optional)
ccc      call psmoo(depths)
c
c     call prtmsk(ip,depths,util1,idm,ii1,jj,0.,1.,
c    .     'bottom depth (m)')
c
      write (*,'(2a)') 'read lat/lon from ',flnmlat
      call findunit(iu2)
      open (iu2,file=flnmlat,form='unformatted',status='old')
      read (iu2) iz,jz
      if (iz.ne.idm .or. jz.ne.jdm) then
        write (*,'(2(a,2i5))') 'error - idm,jdm =',iz,jz,
     .   ' in lat/lon file should be',idm,jdm
        stop '(geopar)'
      end if
      rewind iu2
      read (iu2) iz,jz,latij,lonij
      close(iu2)
c
c     write (*,*) 'shown below: latitude of vorticity points'
c     call zebra(latij(1,1,4),idm,ii,jj)
c     write (*,*) 'shown below: longitude of vorticity points'
c     call zebra(lonij(1,1,4),idm,ii,jj)
c
c --- define coriolis parameter and grid size
      do 56 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 56 i=1,ii
c
      corio(i,j)=sin(latij(i,j,4)*radian)*4.*pi/86164.        !  86400 * 365 / 366
c
      scpy(i,j)=sphdis(latij(i,j ,2),lonij(i,j ,2),
     .                 latij(i,jb,2),lonij(i,jb,2))
c
      if (i.lt.ii) then
        scpx(i,j)=sphdis(latij(i  ,j,1),lonij(i  ,j,1),
     .                   latij(i+1,j,1),lonij(i+1,j,1))
        loncor(1:2) = (/ lonij(i:i+1:+1,j ,4) /) ! List the 4 lons/lats
        latcor(1:2) = (/ latij(i:i+1:+1,j ,4) /) ! of cell corners in
        loncor(3:4) = (/ lonij(i+1:i:-1,jb,4) /) ! counterclockwise order
        latcor(3:4) = (/ latij(i+1:i:-1,jb,4) /) ! for the area calculation
#ifdef CUBED_SPHERE
c Define cell areas using great-circle polygons.
c In future, will be done for all cases.
      n = 4
      do nn=1,4  ! check whether this cell is a triangle
        if(abs(loncor(nn)-loncor(1+mod(nn,4))).lt.1d-3 .and.
     &     abs(latcor(nn)-latcor(1+mod(nn,4))).lt.1d-3) then
          loncor = cshift(loncor,nn)
          latcor = cshift(latcor,nn)
          n = 3
          exit
        endif
      enddo
      call gc_polyarea(loncor,latcor,n,scp2(i,j))
#else
        scp2(i,j)=sphrec(latcor,loncor,succes)
#endif
        scp2i(i,j)=1./scp2(i,j)
      end if
c
      scuy(i,j)=sphdis(latij(i,j ,4),lonij(i,j ,4),
     .                 latij(i,jb,4),lonij(i,jb,4))
c
      if (i.gt.1) then
        scux(i,j)=sphdis(latij(i  ,j,3),lonij(i  ,j,3),
     .                   latij(i-1,j,3),lonij(i-1,j,3))
        scu2(i,j)=scux(i,j)*scuy(i,j)
        scuxi(i,j)=1./scux(i,j)
      end if
c
      scvy(i,j)=sphdis(latij(i,j ,3),lonij(i,j ,3),
     .                 latij(i,ja,3),lonij(i,ja,3))
c
      if (i.lt.ii) then
        scvx(i,j)=sphdis(latij(i  ,j,4),lonij(i  ,j,4),
     .                   latij(i+1,j,4),lonij(i+1,j,4))
        scv2(i,j)=scvx(i,j)*scvy(i,j)
        scvyi(i,j)=1./scvy(i,j)
      end if
c
      scqy(i,j)=sphdis(latij(i,j ,1),lonij(i,j ,1),
     .                 latij(i,ja,1),lonij(i,ja,1))
c
      if (i.gt.1) then
        scqx(i,j)=sphdis(latij(i  ,j,2),lonij(i  ,j,2),
     .                   latij(i-1,j,2),lonij(i-1,j,2))
        scq2(i,j)=scqx(i,j)*scqy(i,j)
        scq2i(i,j)=1./scq2(i,j)
      end if
c
 56   continue
c
      if (beropn .and. scu2(ipacs,jpac).ne.scu2(iatln,jatl))
     .  write(*,'(a,6f13.5)') ' chk WRONG scu2'
     . ,scu2(ipacs,jpac),scu2(iatln,jatl)
c
      if ( latij(ipacs,jpac,3).ne.latij(iatls,jatl,3)
     . .or.latij(ipacn,jpac,3).ne.latij(iatln,jatl,3)
     . .or.lonij(ipacs,jpac,3).ne.lonij(iatls,jatl,3)
     . .or.lonij(ipacn,jpac,3).ne.lonij(iatln,jatl,3))
     .  write(*,'(a,8f9.2)') ' chk WRONG lat/lon '
     .,latij(ipacs,jpac,3),latij(iatls,jatl,3)
     .,latij(ipacn,jpac,3),latij(iatln,jatl,3)
     .,lonij(ipacs,jpac,3),lonij(iatls,jatl,3)
     .,lonij(ipacn,jpac,3),lonij(iatln,jatl,3)
c
      write(*,'(a,3f8.2)') 'lat/lon/depth of Bering Strait:'
     . ,latij(ipacs,jpac,3),lonij(ipacs,jpac,3),depths(ipacs,jpac)
c     write (*,'('' shown below: coriolis parameter'')')
c     call zebra(corio,idm,ii,jj)
c     write (*,'('' shown below: grid cell size'')')
c     call zebra(scp2,idm,ii,jj)
c
      area=0.
      zonarea(:)=0.
      avgbot=0.
      ocnvol=0.
      zone(:,:,:)=0.
c
      do 57 j=1,jj
      do 57 l=1,isp(j)
      do 57 i=ifp(j,l),ilp(j,l)
      ocnvol=ocnvol+depths(i,j)*scp2(i,j)
      area=area+scp2(i,j)
      if (latij(i,j,3).gt.40.) then		! arctic zone is 1
        zone(i,j,1)=1.
        zonarea(1)=zonarea(1)+scp2(i,j)
      else if (latij(i,j,3).lt.-40.) then	! antarctic zone is 3
        zone(i,j,3)=1.
        zonarea(3)=zonarea(3)+scp2(i,j)
      else					! remaining area is 2
        zone(i,j,2)=1.
        zonarea(2)=zonarea(2)+scp2(i,j)
      end if
 57   continue
      avgbot=ocnvol/area
      write (*,100) avgbot,area*1.e-12,ocnvol*1.e-18
 100  format('mean basin depth(m), area(Mm^2) & volumn(Mm^3):',3f9.3)
      write (*,101) (zonarea(l)*1.e-12,l=1,3)
 101  format ('latitude zones (arctic,antarctic,other):',3f9.3)

      write(*,'(a,es22.14)') "Total surface area in hycom model (m2)",
     &      sum(scp2)
      write(*,'(a,es22.14)') "Ocean surface area in hycom model (m2)",
     &      sum(scp2, mask = depths > 0.)
      write(*,'(a,es12.3)') "(Ocean surface area)/"
     &      //"(Total surface area) %",
     &      100*sum(scp2, mask = depths > 0.)/sum(scp2)
c
c --- initialize some arrays
c
      call init_pr1d()  ! TNL : isobaric depth levels

      ! uncommented by IA
      !if (nstep0.eq.0) then
      if (iniOCEAN) then
      write (*,*) 'laying out arrays in memory ...'
      do 209 j=1,jj
      do 209 i=1,ii
      p(i,j,1)=huge
      if (ip(i,j).eq.1) p(i,j,1)=zero
      pu(i,j,1)=huge
      pv(i,j,1)=huge
      pbot(i,j)=huge
      ubavg(i,j,1)=huge
      ubavg(i,j,2)=huge
      ubavg(i,j,3)=huge
      vbavg(i,j,1)=huge
      vbavg(i,j,2)=huge
      vbavg(i,j,3)=huge
      utotm(i,j)=huge
      vtotm(i,j)=huge
      utotn(i,j)=huge
      vtotn(i,j)=huge
      uflux (i,j)=huge
      vflux (i,j)=huge
      uflux1(i,j)=huge
      vflux1(i,j)=huge
      uflux2(i,j)=huge
      vflux2(i,j)=huge
      uflux3(i,j)=huge
      vflux3(i,j)=huge
      uja(i,j)=huge
      ujb(i,j)=huge
      via(i,j)=huge
      vib(i,j)=huge
      pgfx(i,j)=huge
      pgfy(i,j)=huge
      depthu(i,j)=huge
      depthv(i,j)=huge
      tprime(i,j)=huge
c
      srfhgt(i,j)=huge
      dpmixl(i,j,:)= 1.0  ! TNL: avoid NaN on the first step
      oice(i,j)=huge
      taux(i,j)=huge
      tauy(i,j)=huge
      oflxa2o(i,j)=huge
      osalt(i,j)=huge
      oemnp(i,j)=huge
      ustar(i,j)=huge
      sswflx(i,j)=huge
c
      pbavav(i,j)=huge
      sfhtav(i,j)=huge
      dpmxav(i,j)=huge
      oiceav(i,j)=huge
      eminpav(i,j)=huge
      surflav(i,j)=huge
      salflav(i,j)=huge
      brineav(i,j)=huge
      tauxav(i,j)=huge
      tauyav(i,j)=huge
c
      do 209 k=1,kk
      u  (i,j,k   )=huge
      u  (i,j,k+kk)=huge
      v  (i,j,k   )=huge
      v  (i,j,k+kk)=huge
      uflx(i,j,k)=huge
      vflx(i,j,k)=huge
      ufxcum(i,j,k)=huge
      vfxcum(i,j,k)=huge
      dpinit(i,j,k)=huge
      dpold (i,j,k)=huge
      dp (i,j,k   )=huge
      dp (i,j,k+kk)=huge
      dpu(i,j,k   )=huge
      dpu(i,j,k+kk)=huge
      dpv(i,j,k   )=huge
      dpv(i,j,k+kk)=huge
      p (i,j,k+1)=huge
      pu(i,j,k+1)=huge
      pv(i,j,k+1)=huge
c
      th3d(i,j,k+kk)=huge
      th3d(i,j,k   )=huge
      thstar(i,j,k)=huge
      do nt=1,ntrcr
        tracer(i,j,k,nt)=zero
      end do
      uav(i,j,k)=huge
      vav(i,j,k)=huge
      dpuav(i,j,k)=huge
      dpvav(i,j,k)=huge
      dpav (i,j,k)=huge
      temav(i,j,k)=huge
      salav(i,j,k)=huge
      th3av(i,j,k)=huge
      uflxav(i,j,k)=huge
      vflxav(i,j,k)=huge
      ufxavp(i,j,k)=huge
      vfxavp(i,j,k)=huge
      diaflx(i,j,k)=huge
 209  continue
c
      do 210 j=1,jj
      ja=mod(j-2+jj,jj)+1
      do 210 l=1,isq(j)
      do 210 i=ifq(j,l),ilq(j,l)
      pbot(i  ,j  )=0.
      pbot(i-1,j  )=0.
      pbot(i  ,ja )=0.
      pbot(i-1,ja )=0.
      p(i  ,j  ,1)=0.
      p(i-1,j  ,1)=0.
      p(i  ,ja ,1)=0.
      p(i-1,ja ,1)=0.
      do 210 k=1,kk
      dp(i  ,j  ,k   )=0.
      dp(i  ,j  ,k+kk)=0.
      dp(i-1,j  ,k   )=0.
      dp(i-1,j  ,k+kk)=0.
      dp(i  ,ja ,k   )=0.
      dp(i  ,ja ,k+kk)=0.
      dp(i-1,ja ,k   )=0.
 210  dp(i-1,ja ,k+kk)=0.
c
c --- initialize  u,ubavg,utotm,uflx,uflux,uflux2/3,uja,ujb  at points
c --- located upstream and downstream (in i direction) of p points.
c --- initialize  depthu,dpu,utotn,pgfx  upstream and downstream of p points
c --- as well as at lateral neighbors of interior u points.
c
      do 156 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 156 l=1,isu(j)
      do 156 i=ifu(j,l),ilu(j,l)
      pu(i,j,1)=0.
c
      depthu(i,ja)=0.
      utotn (i,ja)=0.
      pgfx  (i,ja)=0.
c
      depthu(i,jb)=0.
      utotn (i,jb)=0.
      pgfx  (i,jb)=0.
c
      do 156 k=1,kk
      dpu(i,ja,k   )=0.
      dpu(i,ja,k+kk)=0.
c
      dpu(i,jb,k   )=0.
      dpu(i,jb,k+kk)=0.
 156  continue
c
      do 158 j=1,jj
      do 158 l=1,isp(j)
      do 158 i=ifp(j,l),ilp(j,l)+1
      depthu(i,j)=0.
      utotn (i,j)=0.
      pgfx  (i,j)=0.
      ubavg(i,j,1)=0.
      ubavg(i,j,2)=0.
      ubavg(i,j,3)=0.
      utotm (i,j)=0.
      uflux (i,j)=0.
      uflux2(i,j)=0.
      uflux3(i,j)=0.
      uja(i,j)=0.
      ujb(i,j)=0.
c
      do 158 k=1,kk
      dpu(i,j,k   )=0.
      dpu(i,j,k+kk)=0.
      uflx(i,j,k)=0.
      ufxcum(i,j,k)=0.
      u(i,j,k   )=0.
 158  u(i,j,k+kk)=0.
c
c --- initialize  v,vbavg,vtotm,vflx,vflux,vflux2/3,via,vib  at points
c --- located upstream and downstream (in j direction) of p points.
c --- initialize  depthv,dpv,vtotn,pgfy  upstream and downstream of p points
c --- as well as at lateral neighbors of interior v points.
c
      do 166 i=1,ii1
      ia=mod(i-2+ii,ii)+1
      ib=i+1
      do 166 l=1,jsv(i)
      do 166 j=jfv(i,l),jlv(i,l)
      pv(i,j,1)=0.
c
      depthv(ia,j)=0.
      vtotn (ia,j)=0.
      pgfy  (ia,j)=0.
c
      depthv(ib,j)=0.
      vtotn (ib,j)=0.
      pgfy  (ib,j)=0.
c
      do 166 k=1,kk
      dpv(ia,j,k   )=0.
      dpv(ia,j,k+kk)=0.
c
      dpv(ib,j,k   )=0.
      dpv(ib,j,k+kk)=0.
 166  continue
c
      do 168 i=1,ii1
      do 168 l=1,jsp(i)
      do 168 jp=jfp(i,l),jlp(i,l)+1
      j=mod(jp-1+jj,jj)+1
      depthv(i,j)=0.
      vtotn (i,j)=0.
      pgfy  (i,j)=0.
      vbavg(i,j,1)=0.
      vbavg(i,j,2)=0.
      vbavg(i,j,3)=0.
      vtotm (i,j)=0.
      vflux (i,j)=0.
      vflux2(i,j)=0.
      vflux3(i,j)=0.
      via(i,j)=0.
      vib(i,j)=0.
c
      do 168 k=1,kk
      dpv(i,j,k   )=0.
      dpv(i,j,k+kk)=0.
      vflx(i,j,k)=0.
      vfxcum(i,j,k)=0.
      v(i,j,k   )=0.
 168  v(i,j,k+kk)=0.
      write (*,*) '... array layout completed'
      ! uncommented by IA
      endif                    ! end of nstep=0
c
c --- set 'glue' to values > 1 in regions where extra viscosity is needed
c
      do 154 j=1,jj
      ja=mod(j-2+jj,jj)+1
      jb=mod(j     ,jj)+1
      do 154 i=1,ii
      ia=mod(i-2+ii,ii)+1
      ib=mod(i     ,ii)+1
      glue(i,j)=1.
c --- add glue in coastal areas
      if (depths(i,j).lt.200.
     .    .or. depths(ia,j).lt.200. .or. depths(i,ja).lt.200.
     .    .or. depths(ib,j).lt.200. .or. depths(i,jb).lt.200.)
     .    glue(i,j)=max(glue(i,j),glufac)
      if (depths(i,j).lt.100.
     .    .or. depths(ia,j).lt.100. .or. depths(i,ja).lt.100.
     .    .or. depths(ib,j).lt.100. .or. depths(i,jb).lt.100.)
     .    glue(i,j)=max(glue(i,j),2.*glufac)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- add glue to mediterranean:
#ifdef HYCOM2deg
        if (i.ge.  91 .and. i.le.  98 .and. j.le.  25)
     .        glue(i,j)=glufac
#endif
#ifdef HYCOM1degRefined
        if ((i.ge. 180 .and. i.le. 198 .and. j.le.  37)
     . .or. (i .ge.188 .and. i.le. 191 .and. j.ge. 356))
     .        glue(i,j)=glufac
#endif
#ifdef HYCOM1degUnrefined
        if ((i.ge. 180 .and. i.le. 198 .and. j.le.  37)
     . .or. (i .ge.188 .and. i.le. 191 .and. j.ge. 356))
     .        glue(i,j)=glufac
#endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 154  continue
c
c --- 1:9 represent NAT, SAT, NIN, SIN, NPA, SPA, ARC, SO, MED
      call findunit(iu3)
      open (iu3,file=flnmbas,form='formatted',status='old')
#ifdef HYCOM2deg
        do n=1,2	! reading in basinmask in 2 columns
        read(iu3,*)
        read(iu3,'(90i1)') ((msk(i,j),j=(n-1)*jj/2+1,n*jj/2),i=1,ii)
        enddo
#endif
#ifdef HYCOM1degRefined
!       do n=1,2	! read basinmask in 2 columns
!       read(iu3,*)
!       read(iu3,'(4x,180i1)') ((msk(i,j),j=(n-1)*jj/2+1,n*jj/2),i=1,ii)
!       enddo
        do n=1,3	! read basinmask in 3 columns
        read(iu3,*)
        read(iu3,'(4x,120i1)') ((msk(i,j),j=(n-1)*jj/3+1,n*jj/3),i=1,ii)
        enddo
#endif
#ifdef HYCOM1degUnrefined
        do n=1,3	! read basinmask in 3 columns
        read(iu3,*)
        read(iu3,'(4x,120i1)') ((msk(i,j),j=(n-1)*jj/3+1,n*jj/3),i=1,ii)
        enddo
#endif
      close(iu3)
c
      do i=1,ii
      do j=1,jj
      ijlist(i,j)=1000*i+j
      enddo
      enddo
c
      write(*,'(a,20i12)') 'ijlist ',((ijlist(i,j),i=30,32),j=4,5)
c
      wgtkap=0.
      do 159 j=1,jj
      do 159 l=1,isp(j)
      do 159 i=ifp(j,l),ilp(j,l)+1
c
c --- in indopacific, wgtkap varies between 2 in the south and 3 in the north
c --- in atlantic, wgtkap varies between 2 in the south and 1 in the north
c --- in mediterranean, wgtkap is set to 4
c
c --- linear variation between 30 S and 30 N
      q=min(1.,max(0.,(latij(i,j,3)+30.)/60.))
c
      wgtkap(i,j)=2.
      if (msk(i,j).eq.1.or.msk(i,j).eq.2.or.msk(i,j).eq.7) then ! Atl. & Arctic
        wgtkap(i,j)=2.*(1.-q)+1.*q
      elseif (msk(i,j).ge.3.and.msk(i,j).le.6) then ! Pacific & Indian
        wgtkap(i,j)=2.*(1.-q)+3.*q
      elseif (msk(i,j).eq.9) then       ! Med
        wgtkap(i,j)=4.
      endif
 159  continue
c     call prtmsk(ip,wgtkap,util1,idm,ii1,jj,0.,100.,
c    .     'wgtkap')
c
      endif ! AM_I_ROOT
c
      call cpl_wgt                      ! read in weights for coupler
      call init_hycom_dynsi_cpler
c
      call broadcast(ogrid,area)

      return
      end
c
      real function sphdis(lat1,lon1,lat2,lon2)
c
c --- great-circle distance between 2 points on sphere
c
      USE CONSTANT, only: rad_earth => radius,radian    ! radian = pi/180.
      implicit none
      real,intent(IN) :: lat1,lon1,lat2,lon2            ! in degrees
      real onex,oney,onez,twox,twoy,twoz,qq
      logical, parameter:: vrbos=.false.
c
c --- step 1: define vertical unit vectors at the 2 locations
      onex=cosd(lat1)*cosd(lon1)
      oney=cosd(lat1)*sind(lon1)
      onez=sind(lat1)
      qq=1./sqrt(onex**2+oney**2+onez**2)
      onex=onex*qq
      oney=oney*qq
      onez=onez*qq

      twox=cosd(lat2)*cosd(lon2)
      twoy=cosd(lat2)*sind(lon2)
      twoz=sind(lat2)
      qq=1./sqrt(twox**2+twoy**2+twoz**2)
      twox=twox*qq
      twoy=twoy*qq
      twoz=twoz*qq

c --- step 2: get the angle between the 2 vectors from their dot product
      sphdis=acos(max(-1.,min(1.,onex*twox+oney*twoy+onez*twoz)))
     .  *rad_earth
      sphdis=max(1.,sphdis)
      if (vrbos) print '(a,es11.3,2(3x,a,2f9.3))','sphdis (km):',
     .  sphdis*6371./rad_earth,'lat1/2:',lat1,lat2,'lon1/2:',lon1,lon2

      return
      end

      real function sphtri(rlat,rlon,succes)
c
c --- compute area of spherical triangle delineated by 3 lat/lon points.
c
      USE CONSTANT, only: rad_earth => radius,radian    ! radian = pi/180.
      USE HYCOM_SCALARS, only : pi
      implicit none
      real   ,intent(IN)  :: rlat(3),rlon(3)		! in degrees
      logical,intent(OUT) :: succes
      real vec(3,3),hor1(3),hor2(3),hlfway(3),ang(3),qq
      integer n,n1,n2
CTNL  logical, parameter:: vrbos=.true.
      logical, parameter:: vrbos=.false.
c
      if (vrbos) print '(/a,3f10.4)','(sphtri) corner lat:',rlat
      if (vrbos) print '( a,3f10.4)','(sphtri) corner lon:',rlon
      succes=.true.

c --- make sure corner points are not fused
      do n=1,3
       n1=mod(n,3)+1
       if ((rlat(n)-rlat(n1))**2+(rlon(n)-rlon(n1))**2.lt.1.e-13) then
        print '(a,3f10.4)','(sphtri) corner lat:',rlat
        print '(a,3f10.4)','(sphtri) corner lon:',rlon
        print '(a,2i3,2a)','corners',n,n1,' are are identical',
     .    ' -- set area to zero'
        sphtri=0.
        return
       end if
      end do

c --- step 1: define vertical unit vectors
      do n=1,3
       vec(1,n)=cosd(rlat(n))*cosd(rlon(n))
       vec(2,n)=cosd(rlat(n))*sind(rlon(n))
       vec(3,n)=sind(rlat(n))
       qq=1./sqrt(vec(1,n)**2+vec(2,n)**2+vec(3,n)**2)
       vec(1,n)=vec(1,n)*qq
       vec(2,n)=vec(2,n)*qq
       vec(3,n)=vec(3,n)*qq
      end do
      if (vrbos) print '(a,3f14.5/(22x,3f14.5))',
     .   '(sphtri) vert.vectors:',((vec(n1,n),n=1,3),n1=1,3)
c
c --- step 2: use double cross product ax(bxc)=b(a.c)-c(a.b) to onstruct
c --- horizontal vectors tangential to sphere & parallel to edges
      do n=1,3
       n1=mod(n,3)+1
       call vecprod(vec(1,n),vec(1,n1),hlfway)
       call vecprod(vec(1,n),hlfway,hor1)
       n2=mod(n+1,3)+1
       call vecprod(vec(1,n),vec(1,n2),hlfway)
       call vecprod(vec(1,n),hlfway,hor2)
       if (vrbos) print '(a,i3/(2es17.4))',
     .   '(sphtri) horz.vectors, node',n,(hor1(n1),hor2(n1),n1=1,3)
c
c --- step 3: use cosine form of dot product to get angle between edge pairs
       ang(n)=acos(max(-1.,min(1.,
     .   (hor1(1)*hor2(1)+hor1(2)*hor2(2)+hor1(3)*hor2(3))/
     .   sqrt((hor1(1)**2+hor1(2)**2+hor1(3)**2)
     .       *(hor2(1)**2+hor2(2)**2+hor2(3)**2)))))
      end do

c --- step 4: get area
      sphtri=(ang(1)+ang(2)+ang(3)-pi)*rad_earth**2
      if(sphtri.lt.0.) succes=.false.
      if (vrbos)
     . print '(a,es13.5,a,3f8.3)','(sphtri) cell area [km^2]:',
     .  sphtri*(6371./rad_earth)**2,'  angles:',ang(:)*radian
      return
      end function sphtri
c
c
      real function sphrec(rlat,rlon,succes)
c
c --- compute area of spherical quadrilateral delineated by 4 lat/lon points.
c --- points must be ordered in either direction around quadrilateral
c
      USE CONSTANT, only: rad_earth => radius,radian    ! radian = pi/180.
      USE HYCOM_SCALARS, only : pi
      implicit none
      real   ,intent(INOUT)  :: rlat(4),rlon(4)		! in degrees
      logical,intent(OUT)    :: succes
      real vec(3,4),hor1(3),hor2(3),hlfway(3),ang(4),qq,sphtri
      integer n,n1,n2
CTNL  logical, parameter:: vrbos=.true.
      logical, parameter:: vrbos=.false.
c
      if (vrbos) print '(a,4f10.4)','(sphrec) corner lat:',rlat
      if (vrbos) print '(a,4f10.4)','(sphrec) corner lon:',rlon
      succes=.true.

c --- make sure corner points are not fused
      do n=1,4
       n1=mod(n,4)+1
       if ((rlat(n)-rlat(n1))**2+(rlon(n)-rlon(n1))**2.lt.1.e-13) then
        if (vrbos) then
          print '(a,4f10.4)','(sphrec) corner lat:',rlat
          print '(a,4f10.4)','(sphrec) corner lon:',rlon
          print '(a,2i3,2a)','corners',n,n1,' are are identical',
     .    ' -- treat as triangle'
        end if
        do n2=min(n,n1),3
         rlat(n2)=rlat(n2+1)
         rlon(n2)=rlon(n2+1)
        end do
        sphrec=sphtri(rlat,rlon,succes)
        return
       end if
      end do

c --- step 1: define vertical unit vectors
      do n=1,4
       vec(1,n)=cosd(rlat(n))*cosd(rlon(n))
       vec(2,n)=cosd(rlat(n))*sind(rlon(n))
       vec(3,n)=sind(rlat(n))
       qq=1./sqrt(vec(1,n)**2+vec(2,n)**2+vec(3,n)**2)
       vec(1,n)=vec(1,n)*qq
       vec(2,n)=vec(2,n)*qq
       vec(3,n)=vec(3,n)*qq
      end do
!     if (vrbos) print '(a,4f14.5/(22x,4f14.5))',
!    .   '(sphrec) vert.vectors:',((vec(n1,n),n=1,4),n1=1,3)
c
c --- step 2: use double cross product ax(bxc)=b(a.c)-c(a.b) to onstruct
c --- horizontal vectors tangential to sphere & parallel to edges
      do n=1,4
       n1=mod(n,4)+1
       call vecprod(vec(1,n),vec(1,n1),hlfway)
       call vecprod(vec(1,n),hlfway,hor1)
       n2=mod(n+2,4)+1
       call vecprod(vec(1,n),vec(1,n2),hlfway)
       call vecprod(vec(1,n),hlfway,hor2)
!      if (vrbos) print '(a,i3/(2es17.4))',
!    .   '(sphrec) horz.vectors, node',n,(hor1(n1),hor2(n1),n1=1,3)
c
c --- step 3: use cosine form of dot product to get angle between edge pairs
       ang(n)=acos(max(-1.,min(1.,
     .   (hor1(1)*hor2(1)+hor1(2)*hor2(2)+hor1(3)*hor2(3))/
     .   sqrt((hor1(1)**2+hor1(2)**2+hor1(3)**2)
     .       *(hor2(1)**2+hor2(2)**2+hor2(3)**2)))))
       if (abs(ang(n)*radian-90.).gt.5.) succes=.false.
      end do

c --- step 4: get area
      sphrec=(ang(1)+ang(2)+ang(3)+ang(4)-2.*pi)*rad_earth**2
      if(sphrec.lt.0.) succes=.false.
!     if (vrbos)
!    .  print '(a,es13.5,a,4f8.3)','(sphrec) cell area [km^2]:',
!    .  sphrec*(6371./rad_earth)**2,'  angles:',ang(:)*radian
      return
      end function sphrec
c
c
      subroutine vecprod(vec1,vec2,vecout)
! --- cross product of 2 vectors
      implicit none
      real,intent(IN)  :: vec1(3),vec2(3)
      real,intent(OUT) :: vecout(3)
      vecout(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
      vecout(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
      vecout(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
      return
      end subroutine vecprod
      subroutine gc_polyarea(lon,lat,n,area)
!@sum gc_polyarea calculates the area of a polygon on a sphere
!@+   whose edges are great circles
!@+   M. Kelley
      USE CONSTANT, only: radius,twopi
      implicit none
      integer :: n
      real*8, dimension(n) :: lon,lat ! input: lon,lat in degrees
      real*8 :: area
      real*8, dimension(3) :: vi,vip1,vn,cri,crip1,crn
      integer :: i
      real*8 :: twopibyn,cc

      twopibyn = twopi/n
      vn = v3d(lon(n),lat(n))
      vi = v3d(lon(1),lat(1))
      crn = cross3d(vn,vi)
      cri = crn
      area = 0.
      do i=1,n-1
        vip1 = v3d(lon(i+1),lat(i+1))
        crip1 = cross3d(vi,vip1)
        cc = sum(vi*cross3d(cri,crip1))
        area = area + (twopibyn-atan2(cc,sum(cri*crip1)))
        vi = vip1
        cri = crip1
      enddo
      cc = sum(vi*cross3d(cri,crn))
      area = area + (twopibyn-atan2(cc,sum(cri*crn)))
      area = area*radius*radius
      return
      contains
      function v3d(lon,lat)
      USE CONSTANT, only: radian   ! radian = pi/180.
      real*8 :: lon,lat,v3d(3)
      real*8 :: rlon,rlat  ! convert to radians
      rlon=lon*radian
      rlat=lat*radian
      v3d(1:2) = cos(rlat)*(/cos(rlon),sin(rlon)/)
      v3d(3)   = sin(rlat)
      end function v3d
      function cross3d(v1,v2)
      real*8, dimension(3) :: v1,v2,cross3d
      cross3d = cshift(v1,1)*cshift(v2,-1)-cshift(v1,-1)*cshift(v2,1)
      end function cross3d
      end subroutine gc_polyarea
c
c> Revision history
c>
c> Mar. 2000 - conversion to SI units
c> May  2000 - changed loop 56 from i/j to j/i to improve memory layout
c> Oct. 2000 - added code to compute 'glue'
c> Apr. 2001 - eliminated stmt_funcs.h
c> Dec. 2001 - added clause to assure p=0 in single-width channels (loop 209)
