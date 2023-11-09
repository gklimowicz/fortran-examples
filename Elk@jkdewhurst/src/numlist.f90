
! Copyright (C) 2015 Manh Duc Le, 2017-18 Arsenii Gerasimov, Yaroslav Kvashnin,
! Lars Nordstrom and J. K. Dewhurst. This file is distributed under the terms of
! the GNU General Public License. See the file COPYING for license details.

subroutine numlist(str,n,list)
implicit none
! arguments
character(256), intent(in) :: str
integer, intent(inout) :: n
integer, intent(out) :: list(n)
! local variables
integer i0,i1,i,j,m,ios
! automatic arrays
integer l(n)
i=0
i0=1
do
  m=index(str(i0:),'-')
  if (m.eq.0) then
    i1=256
  else
    i1=i0+m-2
  end if
  l(:)=0
  read(str(i0:i1),*,iostat=ios) l
  if (i.gt.0) then
    do j=list(i)+1,l(1)-1
      if (i.eq.n) goto 10
      i=i+1
      list(i)=j
    end do
  end if
  do j=1,n
    if (l(j).eq.0) exit
    if (i.eq.n) goto 10
    i=i+1
    list(i)=l(j)
  end do
  if (m.eq.0) exit
  i0=i0+m
end do
10 continue
n=i
end subroutine

