      subroutine zebra(array,idim,ii,jj)
c
c --- find nice contour interval resulting in 7 to 10 contour lines and
c --- draw contours on line printer through the following set of grid points:
c
c         array( 1, 1) . . . . . . . . .  array( 1,jj)
c              .                               .
c              .                               .          (plot will appear
c              .                               .          on paper as shown,
c              .                               .          i down, j across,
c              .                               .          mimicking matrix
c              .                               .          layout)
c              .                               .
c         array(ii,jj) . . . . . . . . .  array(ii,jj)
c
c --- ii  may be smaller than  idim, the first (row) dimension of 'array'
c --- in the calling program. thus, plotting of array subsets is possible.
c
      implicit none
      integer,intent(IN) :: idim,ii,jj
      real,   intent(IN) :: array(idim,jj)
      integer i,j
      integer imn,imx,jmn,jmx,imnj(jj),imxj(jj),jmnj(jj),jmxj(jj)
      real sqrt2,contur,q,ratio,amn,amx,amnj(jj),amxj(jj)
      data sqrt2/1.414/
c
      write (*,'(a,3i6)') 'ZEBRA call with arguments',idim,ii,jj
c
      do 1 j=1,jj
      amxj(j)=-1.e33
      amnj(j)= 1.e33
      imnj(j)=-1
      imxj(j)=-1
      do 1 i=1,ii
      if (amxj(j).lt.array(i,j)) then
        amxj(j)=array(i,j)
        imxj(j)=i
        jmxj(j)=j
      end if
      if (amnj(j).gt.array(i,j)) then
        amnj(j)=array(i,j)
        imnj(j)=i
        jmnj(j)=j
      end if
 1    continue
c
      amx=-1.e33
      amn= 1.e33
      imn=-1
      imx=-1
      jmn=-1
      jmx=-1
      do 2 j=1,jj
      if (amx.lt.amxj(j)) then
        amx=amxj(j)
        imx=imxj(j)
        jmx=jmxj(j)
      end if
      if (amn.gt.amnj(j)) then
        amn=amnj(j)
        imn=imnj(j)
        jmn=jmnj(j)
      end if
 2    continue
c
      if (amx.gt.amn) go to 3
      write (*,100) array(1,1)
 100  format (//' field to be contoured is constant ...',1pe15.5/)
      return
c
 3    contur=(amx-amn)/6.			!  tunable
      q=10.**int(log10(contur))
      if (contur.lt.1.) q=q/10.
      ratio=contur/q
      if (ratio.gt.sqrt2*5.)  contur=q*10.
      if (ratio.le.sqrt2*5.)  contur=q*5.
      if (ratio.le.sqrt2*2.)  contur=q*2.
      if (ratio.le.sqrt2)     contur=q
      write (*,101) contur,amn,imn,jmn,amx,imx,jmx
 101  format ('contour interval in plot below is',es9.1,
     . 6x,'min =',es11.3'  at',2i5/48x,'max =',es11.3'  at',2i5)
      call digplt(array,idim,ii,jj,contur)
c
      return
      end
c
c
      subroutine digplt(array,idim,ii,jj,contur)
c
c --- simulate a contour line plot on the printer
c
      implicit none
      integer,intent(IN) :: idim,ii,jj
      real   ,intent(IN) :: array(idim,jj),contur
c
      character*1 digit(130),dig(20)
      data dig/'0',' ','1',' ','2',' ','3',' ','4',' ',
     .         '5',' ','6',' ','7',' ','8',' ','9',' '/
c
c     nchar = number of character increments in 'j' direction
c     ratio = character width / line spacing
c
      integer,parameter :: nchar=74
      real   ,parameter :: ratio=.58
      integer i,ia,j,ja,k,n
      real xinc,yinc,x,y,dx,dy,dxdy,value
c
      xinc=float(jj-1)/(float(nchar)*ratio)
      yinc=float(jj-1)/ float(nchar)
      k=float(nchar)*ratio*float(ii-1)/float(jj-1)+1.00001
      do 1 i=1,k
      x=1.+float(i-1)*xinc
      ia=min(ii-1,int(x))
      dx=x-float(ia)
      do 2 j=1,nchar+1
      y=1.+float(j-1)*yinc
      ja=min(jj-1,int(y))
      dy=y-float(ja)
      dxdy=dx*dy
      value=array(ia,ja)*(1.-dx-dy+dxdy)
     .     +array(ia+1,ja)*(dx-dxdy)
     .     +array(ia,ja+1)*(dy-dxdy)
     .     +array(ia+1,ja+1)*dxdy
      n=mod(mod(int(2.*value/contur+sign(.5,value)),20)+20,20)+1
 2    digit(j)=dig(n)
 1    write (*,100) 'i',' ',(digit(j),j=1,nchar+1),' ','i'
 100  format(1x,130a1)
      return
      end
