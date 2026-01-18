!@sum prtostat extracts some key ocean statistics from
!@+   input files and prints them
!@auth M. Kelley
      program prtostat
      implicit none
      include 'netcdf.inc'
      character(len=80) :: run_name,time_period,jlfile,ijfile
      integer :: iargc,nargs
      integer :: status,varid,jl_fid,ij_fid
      character(len=100) :: fromto
      character(len=132) :: xlabel
      integer :: imo,jmo,lmo
      real*4, dimension(:), allocatable :: lons,lats,zoce
      real*4, dimension(:,:), allocatable ::
     &     sfij,sf_atl,sf_pac,sf_glo
      integer :: j24n,j38n,j48n,j52s,j72s,j54s,l900m,l3km
      integer :: i80w,i60w, i135e,i155e, i65w

      nargs = iargc()
      if(nargs.ne.2) then
        write(6,*) 'usage: prtostat run_name time_period'
        write(6,*) '  e.g. prtostat Exyz JAN1901'
        write(6,*) '       prtostat Exyz ANN1901-1910'
        stop
      endif

      call getarg(1,run_name)
      call getarg(2,time_period)

      jlfile=trim(time_period)//'.ojl'//trim(run_name)//'.nc'
      ijfile=trim(time_period)//'.oij'//trim(run_name)//'.nc'

      write(6,*) 'Searching for the following input files:'
      write(6,*) '   ',trim(jlfile)
      write(6,*) '   ',trim(ijfile)

c
c open input files
c
      call handle_err(nf_open(jlfile,nf_nowrite,jl_fid),
     &     'opening '//trim(jlfile))
      call handle_err(nf_open(ijfile,nf_nowrite,ij_fid),
     &     'opening '//trim(ijfile))

      write(6,*) 'Searching for the following vars/coords in '//
     &     trim(jlfile)//':'
      write(6,*) '   streamfunctions: sf_Atl sf_Pac sf_Glo'
      write(6,*) '   coord names: lato2,zoce (latitude,depth)'

      write(6,*) 'Searching for the following vars/coords in '//
     &     trim(ijfile)//':'
      write(6,*) '   streamfunction: osfij'
      write(6,*) '   coord names: lono2,lato2 (longitude,latitude)'


      call get_dimsize(ij_fid,'lono2',imo)
      call get_dimsize(ij_fid,'lato2',jmo)

      call get_dimsize(jl_fid,'lato2',jmo)
      call get_dimsize(jl_fid,'zoce',lmo)

      allocate(lons(imo),lats(jmo),zoce(lmo))

      call get_var_real(ij_fid,'lono2',lons)
      call get_var_real(ij_fid,'lato2',lats)
      call get_var_real(jl_fid,'zoce',zoce)

      allocate(sfij(imo,jmo),
     &     sf_atl(jmo,lmo),sf_pac(jmo,lmo),sf_glo(jmo,lmo))

      call get_var_real(ij_fid,'osfij',sfij)
      call get_var_real(jl_fid,'sf_Atl',sf_atl)
      call get_var_real(jl_fid,'sf_Pac',sf_pac)
      call get_var_real(jl_fid,'sf_Glo',sf_glo)

      call get_index(lats,jmo,+48.,j48n)
      call get_index(lats,jmo,-52.,j52s)
      call get_index(zoce,lmo,900.,l900m)
      call get_index(zoce,lmo,3000.,l3km)

      call get_index(lons,imo,-80.,i80w)
      call get_index(lons,imo,-60.,i60w)
      call get_index(lons,imo,135.,i135e)
      call get_index(lons,imo,155.,i155e)
      call get_index(lats,jmo,+24.,j24n)
      call get_index(lats,jmo,+38.,j38n)

      call get_index(lons,imo,-65.,i65w)
      call get_index(lats,jmo,-72.,j72s)
      call get_index(lats,jmo,-54.,j54s)

c
c get run ID, time/date info, etc.
c
      xlabel=''; fromto=''
      status = nf_get_att_text(jl_fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(jl_fid,nf_global,'fromto',fromto)

c
c print 
c
      write(6,*)
      write(6,'(a)') xlabel
      write(6,'(a)') fromto
      write(6,*)
      write(6,'(a)') " Key overturning stream function diags (Sv):"
      write(6,'(a46,f7.2)') 
     &     " North Atlantic overturning: 900m 48N:",sf_atl(j48n,l900m)
      write(6,'(a46,f7.2)') 
     &     " North Pacific overturning: 900m 48N:",sf_pac(j48n,l900m)
      write(6,'(a46,f7.2)')
     &     " Antarctic Bottom Water production: 3000m 52S:",
     &     sf_glo(j52s,l3km)
      write(6,*)
      write(6,'(a)') " Key horizontal stream function diags (Sv):"
      write(6,'(a,f10.3)') " Gulf Stream (MAX in (24-38N,60-80W)):",
     &      maxval(sfij(i80w:i60w,j24n:j38n))
     &     -minval(sfij(i80w:i60w,j24n:j38n))
      write(6,'(a,f10.3)') " Kuroshio  (MAX in (24-38N,135-150E)):",
     &      maxval(sfij(i135e:i155e,j24n:j38n))
     &     -minval(sfij(i135e:i155e,j24n:j38n))
      write(6,'(a,f10.3)') " ACC                 (Drakes Passage):",
     &      maxval(sfij(i65w,j72s:j54s))
     &     -minval(sfij(i65w,j72s:j54s))

c
c close input files
c
      status = nf_close(jl_fid)
      status = nf_close(ij_fid)

c
c free workspace
c
      deallocate(lons,lats,zoce)
      deallocate(sfij, sf_atl,sf_pac,sf_glo)

      end program prtostat

      subroutine get_index(xvals,n,x,j)
      implicit none
      integer :: n,j
      real*4 :: xvals(n),x
      integer :: i
      real*4 :: dx,dxmin
      dxmin = 1e30
      do i=1,n
        dx = abs(xvals(i)-x)
        if(dx.lt.dxmin) then
          dxmin = dx
          j = i
        endif
      enddo
      return
      end subroutine get_index
