      subroutine inikpp
ccc   use mod_xc  ! HYCOM communication interface
c
c --- hycom version 2.1
!     USE HYCOM_DIM_GLOB, only : jj,isp,ifp,ilp
      USE HYCOM_SCALARS, only :epsil,onem
!     USE HYCOM_ARRAYS_GLOB
!     USE KPRF_ARRAYS
      implicit none
c
      integer i,j
      include 'kprf_scalars.h'
c
      integer    nzehat,nustar
      parameter (nzehat=890,nustar=192)
c
      real, dimension (0:nzehat+1,0:nustar+1) ::
     & wmt            ! momentum velocity scale table
     &,wst            ! scalar   velocity scale table
      common/kppltr/ wmt,wst
      save  /kppltr/
c
c -------------------------------------------------------------------
c --- initialize large, mc williams, doney kpp vertical mixing scheme
c -------------------------------------------------------------------
c
      real zehat,zeta,am,cm,c22,zetam,as,c33,zetas,usta,
     .     afourth,athird,ahalf
      parameter (afourth=.25,athird=1./3,ahalf=.5)
c
      data am,cm,c22,zetam/1.257,8.380,16.0,-0.2/
      data as,c33,zetas/-28.86,16.0,-1.0/
c
c --- 'vonk'      = von karman constant
c --- 'zmin,zmax' = zehat limits for velocity scale lookup table, m**3/s**3
c --- 'umin,umax' = ustar limits for velocity scale lookup table
c --- 'epsilon'   = vertical coordinate scale factor
c
      vonk   =  0.4
      zmin   = -0.4e-6
      zmax   =  0.0
      umin   =  0.0
      umax   =  0.16
      epsilon=  0.1
c
c --- construct the velocity-scale lookup tables
c
      deltaz = (zmax-zmin)/(nzehat+1)
      deltau = (umax-umin)/(nustar+1)
c
      do i=0,nzehat+1
        zehat=deltaz*i+zmin
        do j=0,nustar+1
          usta=deltau*j+umin
          zeta=zehat/(usta**3+epsil)
          if (zehat.ge.0.) then
            wmt(i,j)=vonk*usta/(1.+c11*zeta)
            wst(i,j)=wmt(i,j)
          else
            if (zeta.gt.zetam) then
              wmt(i,j)=vonk*usta*(1.-c22*zeta)**afourth
            else
              wmt(i,j)=vonk*(am*usta**3-cm*zehat)**athird
            endif
            if (zeta.gt.zetas) then
              wst(i,j)=vonk*usta*(1.-c33*zeta)**ahalf
            else
              wst(i,j)=vonk*(as*usta**3-cs*zehat)**athird
            endif
          endif
        enddo
      enddo
c
c --- set derived constants
      vtc=sqrt(.2/cs/epsilon)/vonk**2/ricr
      cg=cstar*vonk*(cs*vonk*epsilon)**athird
      dp0enh=2.0*dp00
c
      qdif0 =difm0 /difs0
      qdifiw=difmiw/difsiw
c
cdiag write(*,*) 'shown below: lookup table wst'
cdiag call zebra(wst(0,0),nzehat+2,nzehat+2,nustar+2)
cdiag write(*,*) 'shown below: lookup table wmt'
cdiag call zebra(wmt(0,0),nzehat+2,nzehat+2,nustar+2)
c
      call initurb
c
      return
      end
c
c
      subroutine initurb
c
      USE HYCOM_DIM, only : J_0,J_1,isp,ifp,ilp,jchunk
      USE HYCOM_SCALARS, only : onem, jerlv0
      USE HYCOM_ARRAYS_GLOB, only : depths
      USE HYCOM_ARRAYS,      only : latij
      USE KPRF_ARRAYS,   only : jerlov_loc, betard, betabl, redfac
      implicit none
      integer i,j,l
      include 'kprf_scalars.h'
c
      do 199 j=J_0,J_1
      do 199 l=1,isp(j)
      do 199 i=ifp(j,l),ilp(j,l)
c --- map shallow depths ('brown' water) to high jerlov numbers
      jerlov_loc(i,j)=6-max(1,min(5,int(depths(i,j)/15.0)))
      jerlov_loc(i,j)=max(jerlv0,jerlov_loc(i,j))
c --- reduce turbidity of subtropical oceans
      if (abs(abs(latij(i,j,3))-25.).lt.20.) jerlov_loc(i,j)=1
 199  continue
c
c --- red and blue light extinction coefficients (1/pressure units)
c --- for jerlov water types 1 to 5 - fraction of penetrating red light
      betard(1) = 1.0/( 0.35*onem)
      betard(2) = 1.0/( 0.6 *onem)
      betard(3) = 1.0/( 1.0 *onem)
      betard(4) = 1.0/( 1.5 *onem)
      betard(5) = 1.0/( 1.4 *onem)
      betabl(1) = 1.0/(23.0 *onem)
      betabl(2) = 1.0/(20.0 *onem)
      betabl(3) = 1.0/(17.0 *onem)
      betabl(4) = 1.0/(14.0 *onem)
      betabl(5) = 1.0/( 7.9 *onem)
      redfac(1) = 0.58
      redfac(2) = 0.62
      redfac(3) = 0.67
      redfac(4) = 0.77
      redfac(5) = 0.78
      print *,'turbidity parameters initialized'
      return
      end
c
c
c> Revision history:
c>
c> May  2001 - increased nustar and umax by a factor of 4
