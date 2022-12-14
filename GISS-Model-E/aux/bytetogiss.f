c***  ifort -convert big_endian bytetogiss.f -o bytetogiss
      program atogiss
      integer, parameter :: rows=2160,cols=4320,im=720,jm=360
c     *     ,im=360,jm=180
c     *     ,im=180,jm=90
      integer :: icur(rows*cols),i
      real*4 :: arr(im,jm)
      character*80 :: title
      integer :: ratio,iuout

      ratio = cols/im
      write(*,*) "ratio=",ratio
      open(unit=51,file='input',
     &     form='formatted',status='old',action='read')
            

      do i=1,rows*cols
         read (51,'(I)') icur(i)
      enddo

      close (51)
      
      do i=1,im
         do j=1,jm
            arr(i,jm+1-j)=0.
            do k=1,ratio
               ires=(i-1)*ratio+k
               do m=1,ratio
                  jres=(j-1)*ratio+m
                  index=ires+(jres-1)*cols
                  arr(i,jm+1-j)=arr(i,jm+1-j)+icur(index)
               enddo
            enddo
            arr(i,jm+1-j)=arr(i,jm+1-j)/(real(ratio*ratio))
         enddo
      enddo
      
      iuout=100
      open(unit=iuout,file='output.giss',form='unformatted',
     &     status='UNKNOWN')
      title='fraction of soil texture 720x360'
      write(unit=iuout) title,arr
      close(iuout)

      end program atogiss
