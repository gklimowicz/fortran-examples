module cdf_io
use netcdf
use hycom_scalars, only : huge
contains

   subroutine out1cdf(ncid,kdm,array,time,shortname,longname,units)

! --- add  1 - d i m e n s i o n a l  array to a previously opened netcdf file

   implicit none
   integer,intent(IN) :: ncid			! id of opened file
   integer,intent(IN) :: kdm			! spatial dimensions
   real   ,intent(IN) :: array(kdm)
   real   ,intent(IN) :: time
   character(len=*),intent(IN) :: shortname,longname,units

   integer dimids(2),start(2),count(2),fldid,timeid,status,old
   character(len=4) name
   real*4 real4(kdm)

! --- define spatial and time dimensions if unknown at this point

   status=nf90_inq_dimid (ncid,'kdm',dimids(1))
   if (status.eq.nf90_noerr) then
     status=nf90_inquire_dimension (ncid,dimids(1),name,old)
     print 100,'previously defined: ','kdm',old
     if (old.ne.kdm) then
       print 104,trim(name),old,kdm
       stop '(error)'
     end if
   else
     print 100,'defining dimension: ','kdm',kdm
     call errhandl (nf90_def_dim (ncid,'kdm',kdm,dimids(1)))
   end if

   status=nf90_inq_dimid (ncid,'tdm',dimids(2))
   if (status.eq.nf90_noerr) then
     print 102,'previously defined: ','tdm'
   else
     print 100,'defining dimension: ','tdm',1
     call errhandl (nf90_def_dim (ncid,'tdm',1,dimids(2)))
   end if

! --- define 'time' variable if unknown at this point

   status=nf90_inq_varid (ncid,'time',timeid)
   if (status.eq.nf90_noerr) then
     print 102,'previously defined: ','time'
   else
     print 102,'defining variable: ','time'
     call errhandl (nf90_def_var (ncid,'time',nf90_float,dimids(2),timeid))
     print 101,'attributing: ','time','model days'
     call errhandl (nf90_put_att (ncid,timeid,'long_name','model days'))
   end if

! --- declare attributes for the array to be stored

   print 102,'defining variable: ',shortname
   call errhandl (nf90_def_var (ncid,shortname,nf90_float,dimids(1:1),fldid))
   print 101,'attributing: ',shortname,longname
   call errhandl (nf90_put_att (ncid,fldid,'long_name',longname))
   print 101,'attributing: units ',units
   call errhandl (nf90_put_att (ncid,fldid,'units',units))
   print *,'attributing: missing_value ',huge
   call errhandl (nf90_put_att (ncid,fldid,'missing_value',huge))

! --- now store the array

   call errhandl (nf90_enddef (ncid))		! switch to data mode
   print 103,'storing: ','time',time
   call errhandl (nf90_put_var (ncid,timeid,time))
   real4(:)=array(:)
   start = (/  1,  1/)
   count = (/kdm,  1/)
   print 101,'storing: ',shortname,longname
   call errhandl (nf90_put_var (ncid,fldid,real4,start=start,count=count))

    call errhandl (nf90_redef (ncid))		! back to define mode
100 format (a20,a10,'=',i7)
101 format (a20,a10,'=',2x,a)
102 format (a20,a10)
103 format (a20,a10,'=',f9.1)
104 format ('outcdf error: attempt to change ',a,'=',i5,' to',i5)

   return
   end subroutine out1cdf


   subroutine out2cdf(ncid,idm,jdm,array,time,shortname,longname,units)

! --- add  2 - d i m e n s i o n a l  array to a previously opened netcdf file

   implicit none
   integer,intent(IN) :: ncid			! id of opened file
   integer,intent(IN) :: idm,jdm		! spatial dimensions
   real   ,intent(IN) :: array(idm,jdm)
   real   ,intent(IN) :: time
   character(len=*),intent(IN) :: shortname,longname,units

   integer dimids(3),start(3),count(3),fldid,timeid,status,old
   character(len=4) name
   real*4 real4(idm,jdm)

! --- define spatial and time dimensions if unknown at this point

   status=nf90_inq_dimid (ncid,'idm',dimids(1))
   if (status.eq.nf90_noerr) then
     status=nf90_inquire_dimension (ncid,dimids(1),name,old)
     print 100,'previously defined: ','idm',old
     if (old.ne.idm) then
       print 104,trim(name),old,idm
       stop '(error)'
     end if
   else
     print 100,'defining dimension: ','idm',idm
     call errhandl (nf90_def_dim (ncid,'idm',idm,dimids(1)))
   end if

   status=nf90_inq_dimid (ncid,'jdm',dimids(2))
   if (status.eq.nf90_noerr) then
     status=nf90_inquire_dimension (ncid,dimids(2),name,old)
     print 100,'previously defined: ','jdm',old
     if (old.ne.jdm) then
       print 104,trim(name),old,jdm
       stop '(error)'
     end if
   else
     print 100,'defining dimension: ','jdm',jdm
     call errhandl (nf90_def_dim (ncid,'jdm',jdm,dimids(2)))
   end if

   status=nf90_inq_dimid (ncid,'tdm',dimids(3))
   if (status.eq.nf90_noerr) then
     print 102,'previously defined: ','tdm'
   else
     print 100,'defining dimension: ','tdm',1
     call errhandl (nf90_def_dim (ncid,'tdm',1,dimids(3)))
   end if

! --- define 'time' variable if unknown at this point

   status=nf90_inq_varid (ncid,'time',timeid)
   if (status.eq.nf90_noerr) then
     print 102,'previously defined: ','time'
   else
     print 102,'defining variable: ','time'
     call errhandl (nf90_def_var (ncid,'time',nf90_float,dimids(3),timeid))
     print 101,'attributing: ','time','model days'
     call errhandl (nf90_put_att (ncid,timeid,'long_name','model days'))
   end if

! --- declare attributes for the array to be stored

   print 102,'defining variable: ',shortname
   call errhandl (nf90_def_var (ncid,shortname,nf90_float,dimids(1:2),fldid))
   print 101,'attributing: ',shortname,longname
   call errhandl (nf90_put_att (ncid,fldid,'long_name',longname))
   print 101,'attributing: ','units',units
   call errhandl (nf90_put_att (ncid,fldid,'units',units))
   print *,'attributing: missing_value ',huge
   call errhandl (nf90_put_att (ncid,fldid,'missing_value',huge))

! --- now store the array

   call errhandl (nf90_enddef (ncid))		! switch to data mode
   print 103,'storing: ','time',time
   call errhandl (nf90_put_var (ncid,timeid,time))
   real4(:,:)=array(:,:)
   start = (/  1,  1,  1/)
   count = (/idm,jdm,  1/)
   print 101,'storing: ',shortname,longname
   call errhandl (nf90_put_var (ncid,fldid,real4,start=start,count=count))

    call errhandl (nf90_redef (ncid))		! back to define mode
100 format (a20,a10,'=',i7)
101 format (a20,a10,'=',2x,a)
102 format (a20,a10)
103 format (a20,a10,'=',f9.1)
104 format ('outcdf error: attempt to change ',a,'=',i5,' to',i5)

   return
   end subroutine out2cdf


   subroutine out3cdf(ncid,idm,jdm,kdm,array,time,shortname,longname,units)

! --- add  3 - d i m e n s i o n a l  array to a previously opened netcdf file

   implicit none
   integer,intent(IN) :: ncid			! id of opened file
   integer,intent(IN) :: idm,jdm,kdm		! spatial dimensions
   real   ,intent(IN) :: array(idm,jdm,kdm)
   real   ,intent(IN) :: time
   character(len=*),intent(IN) :: shortname,longname,units

   integer dimids(4),start(4),count(4),fldid,timeid,status,old
   character(len=4) name
   real*4 real4(idm,jdm,kdm)

! --- define spatial and time dimensions if unknown at this point

   status=nf90_inq_dimid (ncid,'idm',dimids(1))
   if (status.eq.nf90_noerr) then
     status=nf90_inquire_dimension (ncid,dimids(1),name,old)
     print 100,'previously defined: ','idm',old
     if (old.ne.idm) then
       print 104,trim(name),old,idm
       stop '(error)'
     end if
   else
     print 100,'defining dimension: ','idm',idm
     call errhandl (nf90_def_dim (ncid,'idm',idm,dimids(1)))
   end if

   status=nf90_inq_dimid (ncid,'jdm',dimids(2))
   if (status.eq.nf90_noerr) then
     status=nf90_inquire_dimension (ncid,dimids(2),name,old)
     print 100,'previously defined: ','jdm',old
     if (old.ne.jdm) then
       print 104,trim(name),old,jdm
       stop '(error)'
     end if
   else
     print 100,'defining dimension: ','jdm',jdm
     call errhandl (nf90_def_dim (ncid,'jdm',jdm,dimids(2)))
   end if

   status=nf90_inq_dimid (ncid,'kdm',dimids(3))
   if (status.eq.nf90_noerr) then
     status=nf90_inquire_dimension (ncid,dimids(3),name,old)
     print 100,'previously defined: ','kdm',old
     if (old.ne.kdm) then
       print 104,trim(name),old,kdm
       stop '(error)'
     end if
   else
     print 100,'defining dimension: ','kdm',kdm
     call errhandl (nf90_def_dim (ncid,'kdm',kdm,dimids(3)))
   end if

   status=nf90_inq_dimid (ncid,'tdm',dimids(4))
   if (status.eq.nf90_noerr) then
     print 102,'previously defined: ','tdm'
   else
     print 100,'defining dimension: ','tdm',1
     call errhandl (nf90_def_dim (ncid,'tdm',  1,dimids(4)))
   end if

! --- define 'time' variable if unknown at this point

   status=nf90_inq_varid (ncid,'time',timeid)
   if (status.eq.nf90_noerr) then
     print 102,'previously defined: ','time'
   else
     print 102,'defining variable: ','time'
     call errhandl (nf90_def_var (ncid,'time',nf90_float,dimids(4),timeid))
     print 101,'attributing: ','time','model days'
     call errhandl (nf90_put_att (ncid,timeid,'long_name','model days'))
   end if

! --- declare attributes for the array to be stored

   print 102,'defining variable: ',shortname
   call errhandl (nf90_def_var (ncid,shortname,nf90_float,dimids(1:3),fldid))
   print 101,'attributing: ',shortname,longname
   call errhandl (nf90_put_att (ncid,fldid,'long_name',longname))
   print 101,'attributing: ','units',units
   call errhandl (nf90_put_att (ncid,fldid,'units',units))
   print *,'attributing: missing_value ',huge
   call errhandl (nf90_put_att (ncid,fldid,'missing_value',huge))

! --- now store the array

   call errhandl (nf90_enddef (ncid))		! switch to data mode
   print 103,'storing: ','time',time
   call errhandl (nf90_put_var (ncid,timeid,time))
   real4(:,:,:)=array(:,:,:)
   start = (/  1,  1,  1,  1/)
   count = (/idm,jdm,kdm,  1/)
   print 101,'storing: ',shortname,longname
   call errhandl (nf90_put_var (ncid,fldid,real4,start=start,count=count))

    call errhandl (nf90_redef (ncid))		! back to define mode
100 format (a20,a10,'=',i7)
101 format (a20,a10,'=',2x,a)
102 format (a20,a10)
103 format (a20,a10,'=',f9.1)
104 format ('outcdf error: attempt to change ',a,'=',i5,' to',i5)

   return
   end subroutine out3cdf


   subroutine readcdf(ncid,varname,field,idm,jdm,start,count,check)
!
! --- extract a 2-D field (real*4) from a previously opened netcdf file
!
   implicit none
   integer      ,intent(IN)     :: ncid,idm,jdm,start(3),count(3)
   character*(*),intent(IN)     :: varname
   logical,optional,intent(IN)  :: check
   real*4       ,intent(OUT)    :: field(idm,jdm)
   integer :: iz,jz,tz,id_idm,id_jdm,id_tdm,id_fld

   if (present(check)) then
    if (check) then
! --- verify that dimensions of archived field agree with specified dimensions
     call errhandl (nf90_inq_dimid (ncid, 'idm', id_idm))
     call errhandl (nf90_inq_dimid (ncid, 'jdm', id_jdm))
     call errhandl (nf90_inq_dimid (ncid, 'tdm', id_tdm))

     call errhandl (nf90_inquire_dimension (ncid, id_idm, len=iz))
     call errhandl (nf90_inquire_dimension (ncid, id_jdm, len=jz))
     call errhandl (nf90_inquire_dimension (ncid, id_tdm, len=tz))

     if (iz.ne.idm .or. jz.ne.jdm) then
      print '(2(a,2i5))','(readcdf) array dimensions are',iz,jz,	&
        ', not the expected values',idm,jdm
      stop '(readcdf error)'
     end if

     if (start(3)+count(3)-1.gt.tz) then
      print '(2(a,i5))','(readcdf) requested record',			&
        start(3)+count(3)-1,'  does not exist. records in file:',tz
      stop '(readcdf error)'
     end if
    end if
   end if

   call errhandl (nf90_inq_varid (ncid, trim(varname), id_fld))
   call errhandl (nf90_get_var (ncid, id_fld, field, start=start, count=count))
   return
   end subroutine readcdf


   subroutine errhandl (ret)
   implicit none
   integer, intent(in) :: ret
   if (ret /= NF90_NOERR) then
     write(6,*) nf90_strerror (ret)
     stop 999
   end if
   return
   end subroutine errhandl

end module cdf_io
