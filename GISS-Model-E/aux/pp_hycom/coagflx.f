      subroutine coagflx(uflx,vflx,ylo,southfl,xlo,eastfl)
c
c --- v e c t o r   field coagulation routine.
c --- sums up -uflx,vflx- over intervals between consecutive sign changes.
c --- ignores uflx/vflx values where iu(i,j)=0 or iv(i,j)=0 resp.
c
c --- nonzero values returned in southfl [eastfl] represent zonally
c --- [meridionally] integrated meridional [zonal] fluxes that are
c --- associated with locations given in ylo [xlo]
c 
      use hycom_dimen, only: idm,jdm,ifull,jfull,i,j,k,iu,iv
c
      implicit none
      real,intent(INOUT)  :: uflx(idm,jdm),vflx(idm,jdm)
      real,intent(OUT) :: ylo(idm,jdm),southfl(idm,jdm),
     .                    xlo(idm,jdm),eastfl (idm,jdm)
      integer ia,ib,ja,jb,jfirst,jp
      real wt,huge
      logical smooth,cyclic
      data smooth/.false./
      data huge/1.e33/
c
c --- cyclic domain is indicated by jdm=jfull+2
      cyclic=.false.
      if (jdm.eq.jfull+2) cyclic=.true.
c
c --- ------------
c --- zonal fluxes
c --- ------------
c
      do 22 j=1,min(jdm-1,jfull)
c     write (*,*) ' processing  j =',j
c
      if (smooth) then
c --- smooth in lateral direction
c
       do k=1,2			!  number of smoothing passes (must be even)
c
        do 12 i=1,idm-1
        ia=max(    1,i-1)
        ib=min(idm-1,i+1)
        eastfl(i,j)=huge
        if (iv(i,j).gt.0) then
         if (iv(ia,j).eq.0) ia=i
         if (iv(ib,j).eq.0) ib=i
         eastfl(i,j)=.25*(vflx(ia,j)+vflx(ib,j))+.5*vflx(i,j)
        endif
 12     continue
c
        do 13 i=1,idm-1
        ia=max(  1,i-1)
        ib=min(idm-1,i+1)
        if (iv(i,j).gt.0) then
         if (iv(ia,j).eq.0) ia=i
         if (iv(ib,j).eq.0) ib=i
         vflx(i,j)=.25*(eastfl(ia,j)+eastfl(ib,j))+.5*eastfl(i,j)
        else
         vflx(i,j)=0.
        endif
 13     continue
       end do				!  number of smoothing passes
      end if				!  smooth = .true.
c
      do 8 i=1,idm-1
 8    eastfl(i,j)=0.
c
      eastfl(1,j)=vflx(1,j)
      wt=vflx(1,j)
c
      do 23 i=2,idm-1
c
c --- check whether flux changes sign between i and i-1
      if (vflx(i,j)*vflx(i-1,j).gt.0.) then
        eastfl(i,j)=eastfl(i-1,j)+vflx(i,j)
        eastfl(i-1,j)=0.
        wt=wt+i*vflx(i,j)
        if (i.eq.idm-1) xlo(i,j)=wt/eastfl(i,j)
c
      else				!  sign change found
        xlo(i-1,j)=0.
        if (eastfl(i-1,j).ne.0.) xlo(i-1,j)=wt/eastfl(i-1,j)
        eastfl(i,j)=vflx(i,j)
        wt=vflx(i,j)*i
        if (i.eq.idm-1) xlo(i,j)=i
      endif
 23   continue
 22   continue
c
c --- -----------------
c --- meridional fluxes
c --- -----------------
c
      do 42 i=1,idm-1
cdiag write (*,*) ' processing  i =',i
c
      if (smooth) then
c --- smooth in lateral direction
c	
       do k=1,2                 !  number of smoothing passes (must be even)
        do 14 j=1,min(jdm-1,jfull)
        if (cyclic) then
         ja=mod(j-2+jfull,jfull)+1
         jb=mod(j        ,jfull)+1
        else
         ja=max(  1,j-1)
         jb=min(jdm-1,j+1)
        end if
        southfl(i,j)=huge
        if (iu(i,j).gt.0) then
         if (iu(i,ja).eq.0) ja=j
         if (iu(i,jb).eq.0) jb=j
         southfl(i,j)=.25*(uflx(i,ja)+uflx(i,jb))+.5*uflx(i,j)
        endif
 14     continue
c
        do 15 j=1,min(jdm-1,jfull)
        if (cyclic) then
         ja=mod(j-2+jfull,jfull)+1
         jb=mod(j        ,jfull)+1
        else
         ja=max(  1,j-1)
         jb=min(jdm-1,j+1)
        end if
        if (iu(i,j).gt.0) then
         if (iu(i,ja).eq.0) ja=j
         if (iu(i,jb).eq.0) jb=j
         uflx(i,j)=.25*(southfl(i,ja)+southfl(i,jb))+.5*southfl(i,j)
        else
         uflx(i,j)=0.
        endif
 15     continue
       end do				!  number of smoothing passes
      end if				!  smooth = .true.
c
      do 9 j=1,min(jdm-1,jfull)
 9    southfl(i,j)=0.
c
      if (cyclic) then
c --- find a land point, or a flux sign change location, to start from
        do 75 j=1,jfull
        if (iu(i,j).eq.0) go to 71
 75     continue
c
        do 76 j=1,jfull
        ja=mod(j-2+jfull,jfull)+1
        if (uflx(i,ja)*uflx(i,j).lt.0.) go to 71
 76     continue
        write (*,'(a,i5)') ' no starting point found in row',i
        go to 42
      else
        j=1
      end if
c
 71   jfirst=j
cdiag write (*,'(2(a,i5))') ' row i =',i,'  starts at j =',jfirst
      southfl(i,jfirst)=uflx(i,jfirst)
      wt=uflx(i,jfirst)
c
      do 43 jp=2,min(jdm-1,jfull)
      if (cyclic) then
        j=mod(jp-2+jfirst,jfull)+1
        ja=mod(j-2+jfull,jfull)+1
      else
        j=jp-1+jfirst
        ja=j-1
      end if
c
c --- check whether flux changes sign between j and j-1
      if (uflx(i,ja)*uflx(i,j).gt.0.) then 
        southfl(i,j)=southfl(i,ja)+uflx(i,j)
        southfl(i,ja)=0.
        wt=wt+jp*uflx(i,j)
        if (jp.eq.min(jdm-1,jfull)) ylo(i,j)=wt/southfl(i,j)
      else				!  sign change found
        southfl(i,j)=uflx(i,j)
        ylo(i,ja)=0.
        if (southfl(i,ja).ne.0.) ylo(i,ja)=wt/southfl(i,ja)
        wt=uflx(i,j)*jp
        if (jp.eq.min(jdm-1,jfull)) ylo(i,j)=jp
      endif
 43   continue
 42   continue
c
      return
      end subroutine coagflx
