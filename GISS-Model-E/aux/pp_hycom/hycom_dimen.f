      module hycom_dimen
      implicit none

      integer, public, parameter :: idm=387,jdm=360,kdm=26,ms=15
     .   ,iold=359, ntrcr = 1, equat=229, iia=360,jja=180,k33=33
     ,   ,ifull=idm, jfull=jdm
      integer, public, parameter :: I_0H=1, I_1H=idm, J_0H=1, J_1H=jdm
      integer :: i,j,k,l,n,jchunk
c
c --- information in common block gindex keeps do loops from running into land
      integer, allocatable :: ip(:,:),iu(:,:),iv(:,:),iq(:,:),
     .ifp(:,:),ilp(:,:),isp(:),jfp(:,:),jlp(:,:),jsp(:),
     .ifq(:,:),ilq(:,:),isq(:),jfq(:,:),jlq(:,:),jsq(:),
     .ifu(:,:),ilu(:,:),isu(:),jfu(:,:),jlu(:,:),jsu(:),
     .ifv(:,:),ilv(:,:),isv(:),jfv(:,:),jlv(:,:),jsv(:),
     .msk(:,:)

      contains
c
      subroutine alloc_hycom_dimen

      allocate(
     . ip(I_0H:I_1H,J_0H:J_1H),iu(I_0H:I_1H,J_0H:J_1H), 
     . iv(I_0H:I_1H,J_0H:J_1H),iq(I_0H:I_1H,J_0H:J_1H),
     . ifp(J_0H:J_1H,ms),ilp(J_0H:J_1H,ms),isp(J_0H:J_1H),
     . jfp(I_0H:I_1H,ms),jlp(I_0H:I_1H,ms),jsp(I_0H:I_1H),
     . ifq(J_0H:J_1H,ms),ilq(J_0H:J_1H,ms),isq(J_0H:J_1H), 
     . jfq(I_0H:I_1H,ms),jlq(I_0H:I_1H,ms),jsq(I_0H:I_1H),
     . ifu(J_0H:J_1H,ms),ilu(J_0H:J_1H,ms),isu(J_0H:J_1H),
     . jfu(I_0H:I_1H,ms),jlu(I_0H:I_1H,ms),jsu(I_0H:I_1H),
     . ifv(J_0H:J_1H,ms),ilv(J_0H:J_1H,ms),isv(J_0H:J_1H), 
     . jfv(I_0H:I_1H,ms),jlv(I_0H:I_1H,ms),jsv(I_0H:I_1H),
     . msk(I_0H:I_1H,J_0H:J_1H) )
      end subroutine alloc_hycom_dimen

      end module hycom_dimen
