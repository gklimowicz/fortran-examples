      program MPI_GEVP_InvIter_solve_qp
      USE MPI
      implicit real(16)(a-h,o-z), integer(4)(i-n)
 
! Solve the generalized symmetric indefinite eigenvalue problem
! H x = L S x using an MPI version of the inverse iteration 
! routine gevp_InvIter_qp.  
!
      INTERFACE
         SUBROUTINE MPI_GEVP_InvIter_qp(N, H, S, Etrial, x, Eritz,
     &                                  exitcode)
         integer(4), INTENT(IN)                :: N        !dimension of H, S, size of problem 
         real(16), INTENT(INOUT), DIMENSION(*) :: H        !H partition
         real(16), INTENT(IN), DIMENSION(*)    :: S        !S partition
         real(16), INTENT(IN)                  :: Etrial   !Estimate of desired eigenvalue
         real(16), INTENT(OUT), DIMENSION(*)   :: x        !converged vector
         real(16), INTENT(OUT)                 :: Eritz    !returned eigenvalue
         integer(4), INTENT(OUT)               :: exitcode ! 0 - convergence obtained
         END
      END INTERFACE
 
      real(16), ALLOCATABLE :: H(:), S(:), x(:)
      real(16) Etrial, Eritz
      real(4) :: time1, time2, time3
      character(100) :: Printstr, Spaces_qp
      character(72)  :: LINE
      character(12)  :: real_clock(3)
      integer(4) :: mytype, exitcode
      logical :: ROOT, Already_open
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
      call MPI_type_contiguous(16,MPI_BYTE,mytype,ierr)
      call MPI_type_commit(mytype,ierr)

      ROOT = myrank.eq.0  

***************************************************************************
***************************************************************************

      if(ROOT) THEN
         call DATE_AND_TIME(real_clock(1),real_clock(2))
         write(6,*)' '
         write(6,*)' DATE and TIME: ',real_clock(1:2)
      ENDIF  !ROOT

      one = 1
      open(10,form='formatted',status='old',file='aa')
! Each process reads the input files, only process 0 writes output
      kH = 7
      kS = 8
      open(unit=kH,form='unformatted')
      open(unit=kS,form='unformatted')

! Loopback starts here
10    rewind kH; rewind kS
      read (10,'(a)'), LINE
      read(LINE,*),N, Etrial
      if(N.eq.0) then
         close(10,status='keep')
         go to 1000
      endif
      IF(ROOT) print*
      IF(ROOT)THEN               !list output is by process 0 (ROOT)...
         write(6,'(a)')' Begin ROOT finding step:'
         write(6,'(a)')
     &     '  Step will be done in quadruple precision arithmetic.'
         write(6,'(a)')'  Inverse iteration will be used...'
         write(6,'(a,18x,f12.6)')'  Etrial =', Etrial
!    &      Etrial_str(1:LEN_TRIM(Etrial_str))
         write(6,'(a,i6)')'  Number terms in expansion =',N
         write(6,'(a,i3)')'  Number processors used =', nproc
      ENDIF   !on ROOT
! Now gather the columns of H and S...
! Preparing matrices for MPI_GEVP_InvIter_qp...each process needs to
! determine the size of vector H and vector S.
      call CPU_time(time1)
      LMAX = 0
      do j = myrank+1, N, nproc
         LMAX = LMAX + j
      enddo
! and allocate space.
      ALLOCATE(S(LMAX))
      ALLOCATE(H(LMAX))
      ALLOCATE(x(N))

      j = myrank+1
      kkk = 0
      do k = 1, N
         if(j.eq.k)then
            read(kH) H(kkk+1:kkk+k)
            read(kS) S(kkk+1:kkk+k)        
            kkk = kkk + k
            j = j+nproc
         else
            read(kH)            !skip over this record (column)
            read(kS)            ! "    "    "     "
         endif
      enddo   !on k
      if(.false.) print *,'i am',myrank,'No of terms',kkk
! Synchronize...
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call CPU_time(time2)
      IF(ROOT)print '(f8.3,a)', time2-time1,
     &           ' seconds required to build H and S partitions.'
!      IF(ROOT)print*,'Calling MPI_GEVP_InvIter_qp:'

! Now calculate the lowest eigenvalue...do not save wavefuction (x) for now.  
      call MPI_GEVP_InvIter_qp(N, H, S, Etrial, x, Eritz, exitcode)
      IF(ROOT)THEN
! Root reports the results...
         Printstr = Spaces_qp('(f44.34)',Eritz,11)
         write(6,'(a,i5,a,a)') '   E(', N, ') = ',
     &      Printstr(1:LEN_TRIM(Printstr))
         write(6,*)
         call CPU_time(time3)
         write(6,'(a,f8.3,a)')
     &   '  Time secular equation step:',time3-time2,' seconds.'
      ENDIF   !on ROOT
      DEALLOCATE (H, S, x)
      go to 10
! Close files kH (H-matrix) and kS (S-matrix)
 1000 inquire(kH, OPENED = Already_open)
      if(Already_open) close(kH)               !H matrix
      inquire(kS, OPENED = Already_open)
      if(Already_open) close(kS)               !S matrix
! And deallocate all working space...
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_FINALIZE(ierr)
      call CPU_time(time3)
      IF(ROOT)print '(f8.3,a)', time3-time1,
     &           ' seconds required to to solve the eigenvalue problem.'
* End of Main Program
      END PROGRAM MPI_GEVP_InvIter_solve_qp
