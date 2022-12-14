      subroutine gtdpth(depths,im)
c
      use hycom_dimen, only: idm,jdm,i,j,n
      use const_proc, only: lp,path0,hycomtopo,basinmask,diag
      implicit none
      real, intent (out) :: depths(idm,jdm)
      integer, intent (out) :: im(idm,jdm)    ! basinmask
      integer, parameter :: n1=20,n2=21
c
      real*4 work(idm,jdm)
c
c --- acquire basin depths
      if (diag)
     .   write (lp,'(a,a)') 'subr. gtdpth - open depth file ',hycomtopo
c
      open (unit=n1,file=trim(path0)//hycomtopo,status='old'
     .    ,form='unformatted')
      read (n1) i,j
      if (i.ne.idm .or.j.ne.jdm) then
        write (lp,'(2(a,2i5))') 'depth file dimensions',i,j,
     .   ' should be',idm,jdm
        stop '(depth file dimensions)'
      end if
      rewind n1
      read (n1) i,j,work
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      close(unit=n1)
c
      depths=work
c
      if (diag)
     .   write (lp,'(a,a)') 'subr. gtdpth - open mask file ',basinmask
      open (unit=n2,file=trim(path0)//basinmask,status='old',
     .      form='formatted')
      do n=1,3
        read(n2,*)
        read(n2,'(4x,120i1)')((im(i,j),j=(n-1)*jdm/3+1,n*jdm/3),i=1,idm)
      enddo
      close(n2)

      return
 6    print *,'error reading depth file'
      stop
      end
