      integer :: fid
      integer, intent(in) :: npx,npy,ntiles
      character(len=*) :: varinfo
      integer :: nd,shp(7),rc,vid
      nd = size(shape(arr))
      shp(1:nd) = shape(arr)
      call define_var(fid,dtype,varinfo,nd,shp,
     &     npx,npy,ntiles,rc,vid)