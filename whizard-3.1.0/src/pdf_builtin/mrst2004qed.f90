! $Id: mrst2004qed.f90 7179 2015-08-27 12:32:33Z jr_reuter $

!
! ---- FORTRAN 90 ----
!
! This file is identical to the original mrst2004qed.f with minor
! modifications to account for the variable data file prefix and to
! sanitize the unit handling. In addition, the code has been wrapped in
! a module to avoid symbol collisions (in particular with LHApdf).

module mrst2004qed

contains

  subroutine mrstqed(x,q,mode,upv,dnv,usea,dsea,str,chm,bot, &
       glu,phot,prefix)
!***************************************************************C
!								C
!  This is a package for the new MRST 2004-QED parton           C
!  distributions.                                               C     
!  Reference: A.D. Martin, R.G. Roberts, W.J. Stirling and      C
!  R.S. Thorne, hep-ph/0411040                                  C
!                                                               C
!  There are 2 pdf sets:  proton distributions (mode = 1)       C
!                        neutron distributions (mode = 2)       C
!                                                               C
!  These sets read grids qed6-10gridp.dat and qed6-10gridn.dat  C
!  respectively.                                                C
!  Note that the neutron grid was replaced in January 2007 to   C
!  correct a bug in the original version that affected the      C
!  ubar and dbar distributions at high x.                       C
!								                C
!								C
!  Note the extra argument "phot" for the photon distribution   C
!  in the proton/neutron. As always, the quantity returned      C
!  is x (f,x,Q^2).                                              C
!								C
!   This subroutine uses an improved interpolation procedure    C 
!   for extracting values of the pdf's from the grid            C
!                                                               C
!         Comments to : W.J.Stirling@durham.ac.uk               C
!                                                               C
!***************************************************************C
    implicit real*8(a-h,o-z)
    character(*) prefix
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    q2=q*q
    if(q2.lt.qsqmin.or.q2.gt.qsqmax) print 99,q2
    if(x.lt.xmin.or.x.gt.xmax)       print 98,x
    if(mode.eq.1) then
       call mrst1(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu,phot,prefix)
    elseif(mode.eq.2) then
       call mrst2(x,q2,upv,dnv,usea,dsea,str,chm,bot,glu,phot,prefix)
    endif
99  format('  WARNING:  Q^2 VALUE IS OUT OF RANGE   ','q2= ',e10.5)
98  format('  WARNING:   X  VALUE IS OUT OF RANGE   ','x= ',e10.5)
    return
  end subroutine mrstqed

  subroutine mrst1(x,qsq,upv,dnv,usea,dsea,str,chm, &
       bot,glu,phot,prefix)
    use io_units
    implicit real*8(a-h,o-z)
    character(*) :: prefix
    parameter(nx=49,nq=37,np=9,nqc0=2,nqb0=11,nqc=35,nqb=26)
    real*8 :: f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq), &
         f6(nx,nq),f7(nx,nq),f8(nx,nq),f9(nx,nq),fc(nx,nqc), &
         fb(nx,nqb)
    real*8 :: qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4), &
         cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4), &
         cc9(nx,nq,4,4),ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
    real*8 :: xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
         1d-4,2d-4,4d-4,6d-4,8d-4, &
         1d-3,2d-3,4d-3,6d-3,8d-3, &
         1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
         .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
         .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
         .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
         .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data init/0/
    integer runit, estat
    save
    xsave=x
    q2save=qsq
    if(init.ne.0) goto 10
    runit = free_unit ()
    open(unit=runit,file=prefix // '/qed6-10gridp.dat', &
         status='old', iostat=estat)
    if (estat /= 0) then
       print *, "ERROR opening grid file " // prefix // &
            "/qed6-10gridp.dat in mrst1!!"
       stop
    endif
    do n=1,nx-1
       do m=1,nq
          read(runit,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m), &
               f5(n,m),f7(n,m),f6(n,m),f8(n,m),f9(n,m)
! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea 9=photon
       end do
    end do
    close (runit)
    do m=1,nq
       f1(nx,m)=0.d0
       f2(nx,m)=0.d0
       f3(nx,m)=0.d0
       f4(nx,m)=0.d0
       f5(nx,m)=0.d0
       f6(nx,m)=0.d0
       f7(nx,m)=0.d0
       f8(nx,m)=0.d0
       f9(nx,m)=0.d0
    end do
    do n=1,nx
       xxl(n)=dlog(xx(n))
    end do
    do m=1,nq
       qql(m)=dlog(qq(m))
    end do

    call jeppe1(nx,nq,xxl,qql,f1,cc1)
    call jeppe1(nx,nq,xxl,qql,f2,cc2)
    call jeppe1(nx,nq,xxl,qql,f3,cc3)
    call jeppe1(nx,nq,xxl,qql,f4,cc4)
    call jeppe1(nx,nq,xxl,qql,f6,cc6)
    call jeppe1(nx,nq,xxl,qql,f8,cc8)
    call jeppe1(nx,nq,xxl,qql,f9,cc9)

    emc2=2.045
    emb2=18.5

    do m=1,nqc
       qqlc(m)=qql(m+nqc0)
       do n=1,nx
          fc(n,m)=f5(n,m+nqc0)
       end do
    end do
    qqlc(1)=dlog(emc2)
    call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

    do m=1,nqb
       qqlb(m)=qql(m+nqb0)
       do n=1,nx
          fb(n,m)=f7(n,m+nqb0)
       end do
    end do
    qqlb(1)=dlog(emb2)
    call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)


    init=1
10  continue
      
    xlog=dlog(x)
    qsqlog=dlog(qsq)

    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc9,phot)

    chm=0.d0
    if(qsq.gt.emc2) then 
       call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
    endif

    bot=0.d0
    if(qsq.gt.emb2) then 
       call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
    endif

    x=xsave
    qsq=q2save
    return
50  format(9f10.5)
  end subroutine mrst1

  subroutine mrst2(x,qsq,upv,dnv,usea,dsea,str,chm, &
       bot,glu,phot,prefix)
    use io_units
    implicit real*8(a-h,o-z)
    character(*) prefix
    parameter(nx=49,nq=37,np=9,nqc0=2,nqb0=11,nqc=35,nqb=26)
    real*8 :: f1(nx,nq),f2(nx,nq),f3(nx,nq),f4(nx,nq),f5(nx,nq), &
         f6(nx,nq),f7(nx,nq),f8(nx,nq),f9(nx,nq),fc(nx,nqc), &
         fb(nx,nqb)
    real*8 :: qq(nq),xx(nx),cc1(nx,nq,4,4),cc2(nx,nq,4,4), &
         cc3(nx,nq,4,4),cc4(nx,nq,4,4),cc6(nx,nq,4,4),cc8(nx,nq,4,4), &
         cc9(nx,nq,4,4),ccc(nx,nqc,4,4),ccb(nx,nqb,4,4)
    real*8 :: xxl(nx),qql(nq),qqlc(nqc),qqlb(nqb)
    data xx/1d-5,2d-5,4d-5,6d-5,8d-5, &
         1d-4,2d-4,4d-4,6d-4,8d-4, &
         1d-3,2d-3,4d-3,6d-3,8d-3, &
         1d-2,1.4d-2,2d-2,3d-2,4d-2,6d-2,8d-2, &
         .1d0,.125d0,.15d0,.175d0,.2d0,.225d0,.25d0,.275d0, &
         .3d0,.325d0,.35d0,.375d0,.4d0,.425d0,.45d0,.475d0, &
         .5d0,.525d0,.55d0,.575d0,.6d0,.65d0,.7d0,.75d0, &
         .8d0,.9d0,1d0/
    data qq/1.25d0,1.5d0,2d0,2.5d0,3.2d0,4d0,5d0,6.4d0,8d0,1d1, &
         1.2d1,1.8d1,2.6d1,4d1,6.4d1,1d2, &
         1.6d2,2.4d2,4d2,6.4d2,1d3,1.8d3,3.2d3,5.6d3,1d4, &
         1.8d4,3.2d4,5.6d4,1d5,1.8d5,3.2d5,5.6d5,1d6, &
         1.8d6,3.2d6,5.6d6,1d7/
    data xmin,xmax,qsqmin,qsqmax/1d-5,1d0,1.25d0,1d7/
    data init/0/
    integer runit, estat
    save
    xsave=x
    q2save=qsq
    if(init.ne.0) goto 10
    runit = free_unit ()
    open(unit=runit,file=prefix // '/qed6-10gridn.dat', &
         status='old', iostat=estat)
    if (estat /= 0) then
       print *, "ERROR opening grid file " // prefix // &
            "/qed6-10gridn.dat in mrst1!!"
       stop
    endif
    do n=1,nx-1
       do m=1,nq
          read(runit,50)f1(n,m),f2(n,m),f3(n,m),f4(n,m), &
               f5(n,m),f7(n,m),f6(n,m),f8(n,m),f9(n,m)
! notation: 1=uval 2=val 3=glue 4=usea 5=chm 6=str 7=btm 8=dsea 9=photon
       end do
    end do          
    close (runit)
    do m=1,nq
       f1(nx,m)=0.d0
       f2(nx,m)=0.d0
       f3(nx,m)=0.d0
       f4(nx,m)=0.d0
       f5(nx,m)=0.d0
       f6(nx,m)=0.d0
       f7(nx,m)=0.d0
       f8(nx,m)=0.d0
       f9(nx,m)=0.d0
    end do
    do n=1,nx
       xxl(n)=dlog(xx(n))
    end do
    do m=1,nq
       qql(m)=dlog(qq(m))
    end do

    call jeppe1(nx,nq,xxl,qql,f1,cc1)
    call jeppe1(nx,nq,xxl,qql,f2,cc2)
    call jeppe1(nx,nq,xxl,qql,f3,cc3)
    call jeppe1(nx,nq,xxl,qql,f4,cc4)
    call jeppe1(nx,nq,xxl,qql,f6,cc6)
    call jeppe1(nx,nq,xxl,qql,f8,cc8)
    call jeppe1(nx,nq,xxl,qql,f9,cc9)

    emc2=2.045
    emb2=18.5

    do m=1,nqc
       qqlc(m)=qql(m+nqc0)
       do n=1,nx
          fc(n,m)=f5(n,m+nqc0)
       end do
    end do
    qqlc(1)=dlog(emc2)
    call jeppe1(nx,nqc,xxl,qqlc,fc,ccc)

    do m=1,nqb
       qqlb(m)=qql(m+nqb0)
       do n=1,nx
          fb(n,m)=f7(n,m+nqb0)
       end do
    end do
    qqlb(1)=dlog(emb2)
    call jeppe1(nx,nqb,xxl,qqlb,fb,ccb)

    init=1
10  continue
      
    xlog=dlog(x)
    qsqlog=dlog(qsq)

    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc1,upv)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc2,dnv)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc3,glu)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc4,usea)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc6,str)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc8,dsea)
    call jeppe2(xlog,qsqlog,nx,nq,xxl,qql,cc9,phot)

    chm=0.d0
    if(qsq.gt.emc2) then 
       call jeppe2(xlog,qsqlog,nx,nqc,xxl,qqlc,ccc,chm)
    endif

    bot=0.d0
    if(qsq.gt.emb2) then 
       call jeppe2(xlog,qsqlog,nx,nqb,xxl,qqlb,ccb,bot)
    end if

    x=xsave
    qsq=q2save
    return
50  format(9f10.5)
  end subroutine mrst2

  subroutine jeppe1(nx,my,xx,yy,ff,cc)
    implicit real*8(a-h,o-z)
    dimension xx(nx),yy(my),ff(nx,my),ff1(nx,my),ff2(nx,my), &
         ff12(nx,my),yy0(4),yy1(4),yy2(4),yy12(4),z(16),wt(16,16), &
         cl(16),cc(nx,my,4,4),iwt(16,16)

    data iwt/1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, &
         0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0, &
         -3,0,0,3,0,0,0,0,-2,0,0,-1,0,0,0,0, &
         2,0,0,-2,0,0,0,0,1,0,0,1,0,0,0,0, &
         0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0, &
         0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0, &
         0,0,0,0,-3,0,0,3,0,0,0,0,-2,0,0,-1, &
         0,0,0,0,2,0,0,-2,0,0,0,0,1,0,0,1, &
         -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0, &
         0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0, &
         9,-9,9,-9,6,3,-3,-6,6,-6,-3,3,4,2,1,2, &
         -6,6,-6,6,-4,-2,2,4,-3,3,3,-3,-2,-1,-1,-2, &
         2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0, &
         0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0, &
         -6,6,-6,6,-3,-3,3,3,-4,4,2,-2,-2,-2,-1,-1, &
         4,-4,4,-4,2,2,-2,-2,2,-2,-2,2,1,1,1,1/


    do m=1,my
       dx=xx(2)-xx(1)
       ff1(1,m)=(ff(2,m)-ff(1,m))/dx
       dx=xx(nx)-xx(nx-1)
       ff1(nx,m)=(ff(nx,m)-ff(nx-1,m))/dx
       do n=2,nx-1
          ff1(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff(n-1,m),ff(n,m), &
               ff(n+1,m))
       end do
    end do

    do n=1,nx
       dy=yy(2)-yy(1)
       ff2(n,1)=(ff(n,2)-ff(n,1))/dy
       dy=yy(my)-yy(my-1)
       ff2(n,my)=(ff(n,my)-ff(n,my-1))/dy
       do m=2,my-1
          ff2(n,m)=polderiv(yy(m-1),yy(m),yy(m+1),ff(n,m-1),ff(n,m), &
               ff(n,m+1))
       end do
    end do       

    do m=1,my
       dx=xx(2)-xx(1)
       ff12(1,m)=(ff2(2,m)-ff2(1,m))/dx
       dx=xx(nx)-xx(nx-1)
       ff12(nx,m)=(ff2(nx,m)-ff2(nx-1,m))/dx
       do n=2,nx-1
          ff12(n,m)=polderiv(xx(n-1),xx(n),xx(n+1),ff2(n-1,m),ff2(n,m), &
               ff2(n+1,m))
       end do
    end do

    do n=1,nx-1
       do m=1,my-1
          d1=xx(n+1)-xx(n)
          d2=yy(m+1)-yy(m)
          d1d2=d1*d2

          yy0(1)=ff(n,m)
          yy0(2)=ff(n+1,m)
          yy0(3)=ff(n+1,m+1)
          yy0(4)=ff(n,m+1)

          yy1(1)=ff1(n,m)
          yy1(2)=ff1(n+1,m)
          yy1(3)=ff1(n+1,m+1)
          yy1(4)=ff1(n,m+1)

          yy2(1)=ff2(n,m)
          yy2(2)=ff2(n+1,m)
          yy2(3)=ff2(n+1,m+1)
          yy2(4)=ff2(n,m+1)

          yy12(1)=ff12(n,m)
          yy12(2)=ff12(n+1,m)
          yy12(3)=ff12(n+1,m+1)
          yy12(4)=ff12(n,m+1)

          do k=1,4
             z(k)=yy0(k)
             z(k+4)=yy1(k)*d1
             z(k+8)=yy2(k)*d2
             z(k+12)=yy12(k)*d1d2
          end do
             
          do l=1,16
             xxd=0.
             do k=1,16
                xxd=xxd+iwt(k,l)*z(k)
             end do
             cl(l)=xxd
          end do
          l=0
          do k=1,4
             do j=1,4
                l=l+1
                cc(n,m,k,j)=cl(l)
             end do
          end do
       end do
    end do
    return
  end subroutine jeppe1

  subroutine jeppe2(x,y,nx,my,xx,yy,cc,z)
!--   G.W. 02/07/2007 Allow extrapolation to small x and large q.
    implicit real*8(a-h,o-z)
    dimension xx(nx),yy(my),cc(nx,my,4,4)      
      
    n=locx(xx,nx,x)
    m=locx(yy,my,y)
      
    if (n.gt.0.and.n.lt.nx.and.m.gt.0.and.m.lt.my) then
!--   Do usual interpolation.
       t=(x-xx(n))/(xx(n+1)-xx(n))
       u=(y-yy(m))/(yy(m+1)-yy(m))
       z=0.d0
       do l=4,1,-1
          z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u &
               +cc(n,m,l,2))*u+cc(n,m,l,1)
       end do
         
    else if (n.eq.0.and.m.gt.0.and.m.lt.my) then
!--   Extrapolate to small x.
       call jeppe3(xx(1),y,nx,my,xx,yy,cc,f0)
       call jeppe3(xx(2),y,nx,my,xx,yy,cc,f1)
       if (f0.gt.0.d0.and.f1.gt.0.d0) then
          z = exp(log(f0)+(log(f1)-log(f0))/(xx(2)-xx(1))*(x-xx(1)))
       else
          z = f0+(f1-f0)/(xx(2)-xx(1))*(x-xx(1))
       end if
         
    else if (n.gt.0.and.m.eq.my) then
!--   Extrapolate to large q.
       call jeppe3(x,yy(my),nx,my,xx,yy,cc,f0)
       call jeppe3(x,yy(my-1),nx,my,xx,yy,cc,f1)
       if (f0.gt.0.d0.and.f1.gt.0.d0) then
          z = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))* &
               (y-yy(my)))
       else
          z = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
       end if
         
    else if (n.eq.0.and.m.eq.my) then
!--   Extrapolate to small x AND large q.
       call jeppe3(xx(1),yy(my),nx,my,xx,yy,cc,f0)
       call jeppe3(xx(1),yy(my-1),nx,my,xx,yy,cc,f1)
       if (f0.gt.0.d0.and.f1.gt.0.d0) then
          z0 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))* &
               (y-yy(my)))
       else
          z0 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
       end if
       call jeppe3(xx(2),yy(my),nx,my,xx,yy,cc,f0)
       call jeppe3(xx(2),yy(my-1),nx,my,xx,yy,cc,f1)
       if (f0.gt.0.d0.and.f1.gt.0.d0) then
          z1 = exp(log(f0)+(log(f0)-log(f1))/(yy(my)-yy(my-1))* &
               (y-yy(my)))
       else
          z1 = f0+(f0-f1)/(yy(my)-yy(my-1))*(y-yy(my))
       end if
       if (z0.gt.0.d0.and.z1.gt.0.d0) then
          z = exp(log(z0)+(log(z1)-log(z0))/(xx(2)-xx(1))*(x-xx(1)))
       else
          z = z0+(z1-z0)/(xx(2)-xx(1))*(x-xx(1))
       end if

    else
!--   Set parton distribution to zero otherwise.
       z = 0.d0

    end if
      
    return
  end subroutine jeppe2

!--   G.W. 02/07/2007 Copy of the original jeppe2,
!--   only used for extrapolation.
  subroutine jeppe3(x,y,nx,my,xx,yy,cc,z)
    implicit real*8(a-h,o-z)
    dimension xx(nx),yy(my),cc(nx,my,4,4)      
    n=locx(xx,nx,x)
    m=locx(yy,my,y)
    t=(x-xx(n))/(xx(n+1)-xx(n))
    u=(y-yy(m))/(yy(m+1)-yy(m))
    z=0.d0
    do l=4,1,-1
       z=t*z+((cc(n,m,l,4)*u+cc(n,m,l,3))*u &
            +cc(n,m,l,2))*u+cc(n,m,l,1)
    end do
    return
  end subroutine jeppe3
      
  integer function locx(xx,nx,x)
    implicit real*8(a-h,o-z)
    dimension xx(nx)
    ! if(x.le.xx(1)) then
    if(x.eq.xx(1)) then ! G.W. 02/07/2007
       locx=1
       return
    endif
    ! if(x.ge.xx(nx)) then 
    if(x.eq.xx(nx)) then ! G.W. 02/07/2007
       locx=nx-1  
       return
    endif
    ju=nx+1
    jl=0
1   if((ju-jl).le.1) go to 2
    jm=(ju+jl)/2
    if(x.ge.xx(jm)) then
       jl=jm
    else
       ju=jm
    endif
    go to 1
2   locx=jl
    return
  end function locx

  real*8 function  polderiv(x1,x2,x3,y1,y2,y3)
    implicit real*8(a-h,o-z)
    polderiv=(x3*x3*(y1-y2)-2.0*x2*(x3*(y1-y2)+x1* &
         (y2-y3))+x2*x2*(y1-y3)+x1*x1*(y2-y3))/((x1-x2)*(x1-x3)*(x2-x3))
    return
  end function polderiv
      
end module mrst2004qed
