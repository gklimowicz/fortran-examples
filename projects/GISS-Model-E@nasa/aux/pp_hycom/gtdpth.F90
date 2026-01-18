      subroutine gtdpth(depths,im)
!
      use hycom_dimen, only: idm,jdm,i,j,n
      use const_proc, only: lp,path0,hycomtopo,basinmask,diag
      implicit none
      real, intent (out) :: depths(idm,jdm)
      integer, intent (out) :: im(idm,jdm)    ! basinmask
      integer, parameter :: n1=20,n2=21
!
      real*4 work(idm,jdm)
!
! --- acquire basin depths
      if (diag) write (lp,'(a,a)') 'subr. gtdpth - open depth file ',trim(hycomtopo)
!
      open (unit=n1,file=trim(path0)//trim(hycomtopo),status='old',form='unformatted',	&
          convert='big_endian')
      read (n1) i,j
      if (i.ne.idm .or.j.ne.jdm) then
        write (lp,'(2(a,2i5))') 'depth file dimensions',i,j,' should be',idm,jdm
        stop '(depth file dimensions)'
      end if
      rewind n1
      read (n1) i,j,work
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      close(unit=n1)
!
      depths=work
!
      if (diag) write (lp,'(a,a)') 'subr. gtdpth - open mask file ',trim(basinmask)
      open (unit=n2,file=trim(path0)//trim(basinmask),status='old',form='formatted')
      do n=1,3
        read(n2,*)
        read(n2,'(4x,120i1)')((im(i,j),j=(n-1)*jdm/3+1,n*jdm/3),i=1,idm)
      enddo
      close(n2)

      return
 6    print *,'error reading depth file'
      stop
      end
