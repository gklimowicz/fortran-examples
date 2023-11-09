      subroutine findmx(mask,array,idm,ii,jj,name)
c
c --- find maximum and minimum in 'array'. only check points where mask > 0
c
      implicit none
      integer jmax
      parameter (jmax=2000)
      integer lp,i,j,idm,ii,jj,mask(idm,jj),ipos,jpos,ineg,jneg,
     .        jpoj(jmax),ipoj(jmax),jnej(jmax),inej(jmax)
      real array(idm,jj),difpos,difneg,difpoj(jmax),difnej(jmax),huge
      character name*(*)
      data huge/1.e33/
      common/linepr/lp
c
      if (jj.gt.jmax) stop '(error subr.findmx -- jmax < jj)'
c
      do j=1,jj
         difpos=-huge
         difneg= huge
         do i=1,ii
           if (mask(i,j).gt.0) then
             if (array(i,j).gt.difpos) then
               difpos=array(i,j)
               ipos=i
               jpos=j
             end if
             if (array(i,j).lt.difneg) then
               difneg=array(i,j)
               ineg=i
               jneg=j
             end if
           end if
         end do
c
         difpoj(j)=difpos
         difnej(j)=difneg
         ipoj(j)=ipos
         jpoj(j)=jpos
         inej(j)=ineg
         jnej(j)=jneg
      end do
c
      difpos=-huge
      difneg= huge
      ipos=-1
      jpos=-1
      ineg=-1
      jneg=-1
c
      do j=1,jj
        if (difpoj(j).gt.difpos) then
          difpos=difpoj(j)
          ipos=ipoj(j)
          jpos=jpoj(j)
        end if
        if (difnej(j).lt.difneg) then
          difneg=difnej(j)
          ineg=inej(j)
          jneg=jnej(j)
        end if
      end do
c
      write (lp,'(2a,1p,2(e11.2,2i5))') name,'  min,max =',
     .    difneg,ineg,jneg,difpos,ipos,jpos
c
      return
      end
