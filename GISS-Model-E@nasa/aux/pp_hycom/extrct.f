      subroutine extrct(work,iorign,jorign,array)
c
c --- copy [work(iorign,jorign)...work(iorign+idm-1,jorign+jdm-1)] into 'array'
c --- this version   c y c l i c   in j
      use hycom_dimen, only : idm,jdm,i,j
      implicit none
c
      real*4, intent (in)  :: work(idm,jdm)	!input array is real*4
      real,   intent (out) :: array(idm,jdm)
      integer,intent(IN)   :: iorign,jorign
      integer :: jp
c
      do 1 j=1,jdm
      jp=mod(jorign+j-2,jdm)+1
      do 1 i=1,min(idm,idm-iorign+1)
 1    array(i,j)=work(iorign+i-1,jp)
c
      return
      end
