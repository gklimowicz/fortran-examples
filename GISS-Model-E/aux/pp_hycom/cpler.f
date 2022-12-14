      module hycom_o2a
      use hycom_dimen
      use const_proc
      implicit none
c
      private

      public o2a_wgt,o2a_sfc,o2a_2d,o2a_2dvec,o2a_3d,o2a_3dvec,kij

      public nwgto2a

      public wlisto2a, coso, sino, ilisto2a,  jlisto2a,  nlisto2a

      integer nwgto2a
      parameter (nwgto2a=16)      ! 1deg hycom => 360x180
c
      real*8 wlisto2a(iia,jja,nwgto2a), wlisto2a_f(iia,jja,nwgto2a)
     .    ,coso(idm,jdm),sino(idm,jdm)
      integer ilisto2a(iia,jja,nwgto2a),jlisto2a(iia,jja,nwgto2a)
     .                                 ,nlisto2a(iia,jja)
      integer kij(iia,jja)   ! last k of mass layer in 1x1 deg

      contains

      subroutine o2a_2d(dp,fldo,flda)
c --- mapping 2d from 'o' grid to 'a' grd at depth; 
c --- and wgt has to be modified by dp
c
c --- fldo/dp:  input field from ogcm grid / layer thickness
c     flda: output field onto agcm grid
c
      implicit none
      integer n,ia,ja,io,jo
      real flda(iia,jja),fldo(idm,jdm),dp(idm,jdm),wgt
c
      do 16 ja=1,jja
      do 16 ia=1,iia
      if (kij(ia,ja).eq.0) then  
        flda(ia,ja)=flag
      else
        flda(ia,ja)=0.
c
        wgt=0.
        do 17 n=1,nlisto2a(ia,ja)
        io=ilisto2a(ia,ja,n)
        jo=jlisto2a(ia,ja,n)
        if (diag.and.abs(ia-iatest).le.1.and.abs(ja-jatest).le.1) 
     .    write(*,'(a,4i4,2f7.1,f6.2)')
     .  ' ia,ja1<= ',ia,ja,io,jo,fldo(io,jo),dp(io,jo),wlisto2a(ia,ja,n)
        if (fldo(io,jo).gt.flag.and.dp(io,jo).ge.epsil_p) then
          flda(ia,ja)=flda(ia,ja)+fldo(io,jo)*wlisto2a(ia,ja,n)
          wgt=wgt+wlisto2a(ia,ja,n)
        endif
 17     continue
        if (wgt.gt.0) then
          flda(ia,ja)=flda(ia,ja)/wgt
        else
          flda(ia,ja)=flag
        endif
        if (diag.and.abs(ia-iatest).le.1.and.abs(ja-jatest).le.1) 
     .    write(*,'(a,2i4,2f7.1)')
     .  ' ia,ja2<= ',ia,ja,flda(ia,ja)
      endif
 16   continue
c
      return
      end subroutine o2a_2d
c
      subroutine o2a_sfc(fldo,flda)
c --- mapping surface fields from 'o' grid to 'a' grd
c
c --- fldo:  input field from ogcm grid
c     flda: output field onto agcm grid
c
      implicit none
      integer n,ia,ja
      real flda(iia,jja),fldo(idm,jdm)
c
      do 16 ja=1,jja
      do 16 ia=1,iia
      flda(ia,ja)=0.
c
      do 17 n=1,nlisto2a(ia,ja)
      flda(ia,ja)=flda(ia,ja)+fldo(ilisto2a(ia,ja,n),jlisto2a(ia,ja,n))
     .                                              *wlisto2a(ia,ja,n)
 17   continue
 16   continue
c
      return
      end subroutine o2a_sfc

      subroutine o2a_2dvec(dp,tauxo,tauyo,tauxa,tauya)
c --- mapping 2dvec from 'o' grid to 'a' grd at depth; wgt has to be modified
c --- tauxo/tauyo: input taux (S-ward)/tauy (E-ward) on ogcm grid (N/m*m)
c ---              c-grid
c --- tauxa/tauya:output taux (E-ward)/tauy (N-ward) on agcm grid (N/m*m)
c ---              center grid location
c    
      implicit none
      integer i,j,l,n,ia,ja,jb,io,jo
      real tauxa(iia,jja),tauya(iia,jja),tauxo(idm,jdm),tauyo(idm,jdm)
     .    ,nward(idm,jdm),eward(idm,jdm),dp(idm,jdm),tta,tto,sine,wgt
c
      nward=flag
      eward=flag
c --- rotate taux/tauy to n/e orientation at local p point on panam grid
      do 12 j=1,jdm
      jb=mod(j,jdm)+1
      do 12 i=1,idm-1
      if (tauxo(i  ,j).gt.flag.and.tauyo(i,j ).gt.flag .and.
     .    tauxo(i+1,j).gt.flag.and.tauyo(i,jb).gt.flag) then
      sine=sino(i,j)*sino(i,j)+coso(i,j)*coso(i,j)
      nward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*sino(i,j)
     .           -(tauxo(i,j)+tauxo(i+1,j))*coso(i,j))/(2.*sine)
      eward(i,j)=((tauyo(i,j)+tauyo(i ,jb))*coso(i,j)
     .           +(tauxo(i,j)+tauxo(i+1,j))*sino(i,j))/(2.*sine)
      endif
 12   continue
c
c --- mapping nward/eward from ogcm grid to agcm grid
c
      do 16 ja=1,jja
      do 16 ia=1,iia
      if (kij(ia,ja).eq.0) then
        tauxa(ia,ja)=flag
        tauya(ia,ja)=flag
      else
        tauxa(ia,ja)=0.
        tauya(ia,ja)=0.
        wgt=0.
c
        do 17 n=1,nlisto2a(ia,ja)
        io=ilisto2a(ia,ja,n)
        jo=jlisto2a(ia,ja,n)
        if (eward(io,jo).gt.flag.and.dp(io,jo).ge.epsil_p) then
          tauxa(ia,ja)=tauxa(ia,ja)+eward(io,jo)*wlisto2a(ia,ja,n)
          tauya(ia,ja)=tauya(ia,ja)+nward(io,jo)*wlisto2a(ia,ja,n)
          wgt=wgt+wlisto2a(ia,ja,n)
        endif
 17     continue
c
        if (wgt.gt.0.) then
          tauxa(ia,ja)=tauxa(ia,ja)/wgt
          tauya(ia,ja)=tauya(ia,ja)/wgt
        else
          tauxa(ia,ja)=flag
          tauya(ia,ja)=flag
        endif
      endif
 16   continue
c
      return
      end subroutine o2a_2dvec
c
      subroutine o2a_3d(p,xrho,xijz)
c
c --- transfrom scalar from rho to z space
c
      real,   intent (in)  :: xrho(idm,jdm,kdm),p(idm,jdm,kdm+1)
      real,   intent (out) :: xijz(iia,jja,k33)

      real :: puv(idm,jdm,2*kdm),coord(2*kdm),xyij(idm,jdm,2),
     .   ptop,pbot,pmid,pedg,dpth,dum,xrhoz(idm,jdm,k33)
     .  ,dp(idm,jdm),x2k(idm,jdm,2*kdm),pij(iia,jja)
      integer, dimension(6) :: lrfi,lrfo
c
c --- 'step' controls the appearance of velocity contours
c ---        step = 1.0 produces stairstep type contour lines
c ---                   with discontinuities at the layer interfaces
c ---        step = 0.0 produces gently curved contour lines
      real, parameter :: step=.3
c
      call o2a_sfc(p(1,1,kdm+1),pij)

      if (diag) then
        do k=1,kdm
        write(*,'(2i4,a,i5/(5f9.2))')
     .   itest,jtest, '  p    input grd k=',k,
     .  ((p(i,j,k+1),j=jtest-2,jtest+2),i=itest-2,itest+2)
        enddo
        do k=1,kdm
        write(*,'(2i4,a,i5/(5f9.2))')
     .   itest,jtest, '  xho   input grd k=',k,
     .   ((xrho(i,j,k),j=jtest-2,jtest+2),i=itest-2,itest+2)
        enddo
      endif

      lrfi=0
      lrfo=0
      lrfi(1)=idm
      lrfi(2)=jdm
c
      do 10 k=1,2*kdm
 10   coord(k)=float(k)
c
c --- define pressure at mass points (2 per layer) for use in grdtrns
c
      dpth=z33(2)
      do 26 j=1,jdm
      do 26 i=1,idm
      if(p(i,j,kdm+1).le.0.) goto 26

      xyij(i,j,1)=float(i)
      xyij(i,j,2)=float(j)
      do 25 k=1,kdm
      x2k(i,j,2*k-1)=xrho(i,j,k)
      x2k(i,j,2*k  )=xrho(i,j,k)
      ptop=p(i,j,k  )
      pbot=p(i,j,k+1)
      pmid=.51*ptop+.49*pbot
      pedg=min(ptop+dpth,pmid)
      puv(i,j,2*k-1)=(1.-step)*pmid+step*pedg
      pmid=.51*pbot+.49*ptop
      pedg=max(pbot-dpth,pmid)
      puv(i,j,2*k  )=(1.-step)*pmid+step*pedg
 25   continue
 
c --- put uppermost depth point at sea surface and lowest point at bottom
      puv(i,j,1)=0.
      puv(i,j,2*kdm)=p(i,j,kdm+1)

      if (diag.and.i.eq.itest.and.j.eq.jtest) then
        do k=1,2*kdm
        write(*,'(2i4,a,i3,5f9.2)')
     .   i,j, ' input to grd k=',k,puv(i,j,k),x2k(i,j,k)
        enddo
      endif
 26    continue

      xrhoz=flag
      call grdtrns(0,puv,x2k,dum,idm,jdm,2*kdm,coord,lrfi,xrhoz,xrhoz,
     .             idm,jdm,k33,z33,lrfo,xyij)
c
      if (diag) then
        do k=1,k33
        write(*,'(2i4,a,f8.1/(5f9.2))')
     .   itest,jtest,'  xrhoz after grdtrns z=',z33(k),
     .   ((xrhoz(i,j,k),j=jtest-2,jtest+2),i=itest-2,itest+2)
        enddo
      endif

c --- done with vertical interpolation; next: horizontal interpolation

      do 17 k=1,k33
c
      do j=1,jdm
      do i=1,idm
        if (k.ge.2) then
          dp(i,j)=min(z33(k),p(i,j,kdm+1))-min(z33(k-1),p(i,j,kdm+1))
        elseif (k.eq.1) then
          dp(i,j)=min(z33(2)-z33(1),p(i,j,kdm+1))
        endif
      enddo
      enddo

      call o2a_2d(dp,xrhoz(1,1,k),xijz(1,1,k))

      do i=1,iia
      do j=1,jja
      if (z33(k).gt.pij(i,j)) xijz(i,j,k)=flag
      enddo
      enddo

 17   continue

      if (diag) then
        do k=1,k33
        write(*,'(2i4,a,f8.1/(5f9.2))')
     .   iatest,jatest,'  xyz after o2a z=',z33(k),
     .   ((xijz(i,j,k),i=iatest-2,iatest+2),j=jatest+2,jatest-2,-1)
        enddo
      endif
c
      return
      end subroutine o2a_3d
c
      subroutine o2a_3dvec(p,u,v,uz,vz)
c
c --- transfrom vector from rho to z space
c
      real,intent (in) :: u(idm,jdm,kdm),v(idm,jdm,kdm),p(idm,jdm,kdm+1)
      real,intent (out):: uz(iia,jja,k33),vz(iia,jja,k33)

      real :: puv(idm,jdm,2*kdm),coord(2*kdm),xyij(idm,jdm,2),
     .   ptop,pbot,pmid,pedg,dpth,dum,urhoz(idm,jdm,k33)
     .  ,vrhoz(idm,jdm,k33),u2k(idm,jdm,2*kdm),v2k(idm,jdm,2*kdm)
     .  ,pij(iia,jja),dp(idm,jdm)
      integer, dimension(6) :: lrfi,lrfo
c
c --- 'step' controls the appearance of velocity contours
c ---        step = 1.0 produces stairstep type contour lines
c ---                   with discontinuities at the layer interfaces
c ---        step = 0.0 produces gently curved contour lines
      real, parameter :: step=.3
c
      call o2a_sfc(p(1,1,kdm+1),pij)

      lrfi=0
      lrfo=0
      lrfi(1)=idm
      lrfi(2)=jdm
c
      do 10 k=1,2*kdm
 10   coord(k)=float(k)
c
c --- define pressure at mass points (2 per layer) for use in grdtrns
c
      dpth=z33(2)
      do 26 j=1,jdm
      do 26 i=1,idm
      if(p(i,j,kdm+1).le.0.) goto 26

      xyij(i,j,1)=float(i)
      xyij(i,j,2)=float(j)

      do 25 k=1,kdm
      u2k(i,j,2*k-1)=u(i,j,k)
      u2k(i,j,2*k  )=u(i,j,k)
      v2k(i,j,2*k-1)=v(i,j,k)
      v2k(i,j,2*k  )=v(i,j,k)
      ptop=p(i,j,k  )
      pbot=p(i,j,k+1)
      pmid=.51*ptop+.49*pbot
      pedg=min(ptop+dpth,pmid)
      puv(i,j,2*k-1)=(1.-step)*pmid+step*pedg
      pmid=.51*pbot+.49*ptop
      pedg=max(pbot-dpth,pmid)
      puv(i,j,2*k  )=(1.-step)*pmid+step*pedg
 25   continue
 
c --- put uppermost depth point at sea surface and lowest point at bottom
      puv(i,j,1)=0.
      puv(i,j,2*kdm)=p(i,j,kdm+1)

      if (diag.and.i.eq.itest.and.j.eq.jtest) then
        do k=1,2*kdm
        write(*,'(2i4,a,i3,5f9.2)')
     .   i,j, ' input to grd k=',k,puv(i,j,k),u2k(i,j,k)
        enddo
      endif
 26   continue

      urhoz=flag
      vrhoz=flag
      call grdtrns(0,puv,u2k,dum,idm,jdm,2*kdm,coord,lrfi,urhoz,urhoz,
     .             idm,jdm,k33,z33,lrfo,xyij)
      call grdtrns(0,puv,v2k,dum,idm,jdm,2*kdm,coord,lrfi,vrhoz,vrhoz,
     .             idm,jdm,k33,z33,lrfo,xyij)
c
c --- done with vertical interpolation; next: horizontal interpolation

      do 17 k=1,k33
c
      do i=1,idm
      do j=1,jdm
        if (k.ge.2) then
          dp(i,j)=min(z33(k),p(i,j,kdm+1))-min(z33(k-1),p(i,j,kdm+1))
        elseif (k.eq.1) then
          dp(i,j)=min(z33(2)-z33(1),p(i,j,kdm+1))
        endif
      enddo
      enddo

      call o2a_2dvec(dp,urhoz(1,1,k),vrhoz(1,1,k),uz(1,1,k),vz(1,1,k))

      do i=1,iia
      do j=1,jja
      if (z33(k).gt.pij(i,j)) then
        uz(i,j,k)=flag
        vz(i,j,k)=flag
      endif
      enddo
      enddo

 17   continue

      if (diag) then
        do k=1,k33
        write(*,'(2i4,a,f8.1/(5f9.2))')
     .   iatest,jatest,'  uz after o2a z=',z33(k),
     .  ((uz(i,j,k),i=iatest-2,iatest+2),j=jatest+2,jatest-2,-1)
        enddo
        do k=1,k33
        write(*,'(2i4,a,f8.1/(5f9.2))')
     .   iatest,jatest,'  vz after o2a z=',z33(k),
     .  ((vz(i,j,k),i=iatest-2,iatest+2),j=jatest+2,jatest-2,-1)
        enddo
      endif
c
      return
      end subroutine o2a_3dvec

      subroutine o2a_wgt
      implicit none
      integer :: iz,jz
c
      integer, parameter :: nsize2=16848000   ! 1deg hycom => 360x180

c --- read in all weights
      if (iia*jja*((nwgto2a*2+1)*4+nwgto2a*8).ne.nsize2) then
        write(lp,'(a,i12,a,i12)') 'wrong size in cpler '
     .  ,iia*jja*((nwgto2a*2+1)*4+nwgto2a*8)
     .  ,' should be ',nsize2
        stop ' wrong size in cpler'
      endif
c
      open(24,file=trim(path0)//flnmo2a,form='unformatted',
     .     status='old',access='direct',recl=nsize2)
      read(24,rec=1) ilisto2a,jlisto2a,wlisto2a,nlisto2a
      close(24)
c
      open(28,file=trim(path0)//flnmcoso,form='unformatted',
     .        status='old')
      read(28) iz,jz,coso,sino
      close(28)
      if (iz.ne.idm .or. jz.ne.jdm) then
        write(lp,*) ' iz,jz=',iz,jz
        stop '(wrong iz/jz in cososino.8bin)'
      endif
c
      return
      end subroutine o2a_wgt

      end module hycom_o2a
