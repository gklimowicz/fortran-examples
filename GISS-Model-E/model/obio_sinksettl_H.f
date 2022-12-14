#include "rundeck_opts.h"
      subroutine obio_sinksettl(vrbos,kmax,errcon,i,j,
     &                    kdm,nstep,dtsrc,ddxypo)

      USE obio_dim
      USE obio_incom,only: wsdeth,mgchltouMC
      USE obio_com, only: P_tend,obio_deltat,D_tend,C_tend
     .                   ,obio_P,det,car
     .                   ,dp1d,wsdet,p1d,obio_ws
     .                   ,rhs,cexp,kzc
      use TimeConstants_mod, only: HOURS_PER_DAY, DAYS_PER_YEAR,
     &                             SECONDS_PER_HOUR



      implicit none

      integer, intent(in) ::kdm,nstep
      real, intent(in) :: dtsrc,ddxypo

      integer :: i,j,k,nt,kmax
      real    :: trnd
      logical :: vrbos,errcon


!---------------------------------------------------------------
! --- phyto sinking and detrital settling
!---------------------------------------------------------------

!originally: the sinking term is given in units (m/hr)*(mgr,chl/m3)
!in order to be converted into mgr,chl/m3/hr as the tendency
!terms are in the phytoplankton equations, 
!we need to multiply by dz of each layer:
!  dz(k  ) * P_tend(k  ) = dz(k  ) * P_tend(k  ) - trnd
!  dz(k+1) * P_tend(k+1) = dz(k+1) * P_tend(k+1) + trnd
!this way we ensure conservation of tracer after vertical adjustment
!the /hr factor is bcz the obio timestep is in hrs.
! now: all terms in /s, hence tendencies in /s. July 2016

      if (kmax.le.1) return

       !phyto sinking
       do nt=1,nchl

          do k=kmax+1,2,-1
            !obio_ws is in m/s   July 2016
            obio_ws(k,nt)=.5*(obio_ws(k,nt)+obio_ws(k-1,nt))
          end do
          obio_ws(     1,nt)=0.                !  no flux through sea surface

          !no flux through the sea floor
          !and convert to distance (m/timestep)
          do k=1,kmax
             obio_ws(k,nt)=min(obio_ws(k,nt),p1d(kmax+1)-p1d(k))
!    .                    * baclin/SECONDS_PER_HOUR
     .                    !* baclin        !July 2016
     &                     *dtsrc
          enddo

           do k=1,kmax
            rhs(k,nnut+nt,16)=obio_P(k,nnut+nt)
           enddo
           call advc1d(kmax,obio_P(1,nnut+nt),p1d,
     .                 obio_ws(1,nt),vrbos,errcon)
           if (errcon) then
             write(*,'(a,4i8)')
     .       'error in phyto sinking: nt=',nt,i,j,kmax
             do k=1,kdm
               write (*,'(i5,3e12.4)')
     .           k,dp1d(k),p1d(k),obio_ws(k,nt)
             enddo
           endif

           do k=1,kmax
            rhs(k,nnut+nt,16)=
     ,          (obio_P(k,nnut+nt)-rhs(k,nnut+nt,16))/obio_deltat
           enddo

       enddo   !nchl
       !detrital settling
       do nt=1,ndet
          do k=kmax+1,2,-1
            !detritus settling rates in m/hr
            wsdet(k,nt)=.5*(wsdet(k,nt)+wsdet(k-1,nt))
            !for the moment let wsdet be constant and not depending on T
            !when I tried to change this I would get crossover errors
            !for layers of near-zero thickness
            wsdet(k,nt)=wsdeth(nt)
          end do
          wsdet(1,nt)=0.                !  no flux through sea surface

          !no flux through the sea floor
          !and convert to distance
          do k=1,kmax
             wsdet(k,nt)=min(wsdet(k,nt),p1d(kmax+1)-p1d(k))
!    .                  *baclin/SECONDS_PER_HOUR
     .                  !*baclin
     &                   *dtsrc
          enddo
!need to change that (?) and create an array that will actually
!hold the excess stuff (sediment array) to be used in
!sequastration studies.
!(sediment=stuff(bottom layer)-stuff(layer above)

          do k=1,kmax
           !this is not really rhs, just a temp storage space
           rhs(k,nnut+nchl+nzoo+nt,16)=det(k,nt)
          enddo
          call advc1d(kmax,det(1,nt),p1d,
     .      wsdet(1,nt),vrbos,errcon)
          if (errcon) then
            write(*,'(a,4i8)')
     .     'error in detritus component: nt=',nt,i,j,kmax
             do k=1,kdm
               write (*,'(i5,3e12.4)')
     .           k,dp1d(k),p1d(k),wsdet(k,nt)
             enddo
          stop
          endif

          do k=1,kmax
           !this is really rhs
           rhs(k,nnut+nchl+nzoo+nt,16)=
     .            (det(k,nt)-rhs(k,nnut+nchl+nzoo+nt,16))/obio_deltat
          enddo

       end do  !ndet

      !diagnostic for carbon export at compensation depth
      cexp = 0.
      do  k=1,kzc    
        do nt=nnut+1,nnut+nchl
           cexp = cexp
     .        + obio_P(k,nt)*obio_ws(k,nt-nnut)   
     .        * mgchltouMC              
     .        * 12.d0              
     .        * SECONDS_PER_HOUR    !July 2016
     .        * HOURS_PER_DAY * DAYS_PER_YEAR         
     .        * 1.d-15            !mgm3 -> PgC/yr              
     &        * ddxypo
 
        enddo

      !term2: settling C detritus contribution
      !dont set cexp = 0 here, because adds to before
        nt= 1            !only the for carbon detritus
        cexp = cexp 
     .        + det(k,nt)*wsdet(k,nt)
     .        * SECONDS_PER_HOUR     !July 2016
     .        * HOURS_PER_DAY * DAYS_PER_YEAR
     .        * 1.d-15                 !ugC/l -> PgC/yr
     &        * ddxypo
   
      enddo


      end subroutine obio_sinksettl
