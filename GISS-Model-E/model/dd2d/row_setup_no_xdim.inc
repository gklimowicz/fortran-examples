      type(dist_grid), intent(in) :: grid
      real*8 arr1d(:)
      integer :: jdim
      integer :: nj,hwidth,j1,j2,nl,nk,xdim,ndim
      ndim = size(shape(arr))
      nj = size(arr,jdim)
      hwidth = (nj-(grid%je-grid%js+1))/2
      j1 = 1+hwidth; j2 = nj-hwidth
      nl = 1
      nk = 1
      do xdim=1,jdim-1
        nl = nl*size(arr,xdim)
      enddo
      do xdim=jdim+1,ndim
        nk = nk*size(arr,xdim)
      enddo
