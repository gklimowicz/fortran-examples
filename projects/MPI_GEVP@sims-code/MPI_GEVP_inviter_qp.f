      SUBROUTINE MPI_GEVP_InvIter_qp(N, H, S, Etrial, x, Eritz,
     &                               exitcode)
! Solve for the eigenvalue 'Eritz' and associated eigenvector 'x' of the
! generalized eigenvalue problem  H*x = Eritz*S*x  using inverse 
! iteration and real(16) arithmetic.  Partitions H and S present as 
! vectors in packed upper triangle by columns round-robin format. The 
! number of processes (nproc) and a processes rank (myrank) are obtained
! by calls to MPI_COMM_SIZE and MPI_COMM_RANK, respectively.  N=Nsize is 
! the dimension of the problem.  'Etrial' is the initial guess for the 
! eigenvalue required by the inverse iteration method.  'b' is a working
! vector of length n.  H is replaced by the implicit form of H-inverse.  
      USE MPI
      implicit real(16)(a-h,o-z), integer(4)(i-n)
      integer(4), INTENT(IN)                :: N         !dimension of H, S, size of problem 
      real(16), INTENT(INOUT), DIMENSION(*) :: H         !H partition
      real(16), INTENT(IN), DIMENSION(*)    :: S         !S partition
      real(16), INTENT(IN)                  :: Etrial    !Estimate of desired eigenvalue
      real(16), INTENT(OUT), DIMENSION(*)   :: x         !converged vector
      real(16), INTENT(OUT)                 :: Eritz     !returned eigenvalue
      integer(4), INTENT(OUT)               :: exitcode  ! 0 - convergence obtained
                                                         ! 1 - convergence not obtained
      logical :: timing = .true.                         !set the timing flag
      integer(4) :: Num_Iter = 10                        !number of iterations
      real(16)   :: ConvEPS = 1.0d-27                    !Iteration convergence constant
      integer(4), ALLOCATABLE               :: L0tab(:,:)
      real(16), allocatable :: L0(:)
      real(16) x0(N), b(N), t(N), tt(N), trecv(N)
      real(8) tin, time1, timeiter, tout                 !timing stuff
      character(84) string84
      logical :: ROOT
      integer(4) istat(mpi_status_size)
      character(100) :: Printstr, Spaces_qp
!
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)
      call MPI_TYPE_CONTIGUOUS(16,MPI_BYTE,mytype,ierror)
      call MPI_TYPE_COMMIT(mytype, ierror)
      ROOT = myrank .eq. 0
      if(ROOT.and.(timing)) tin = MPI_WTIME()

      one = 1
      x0 = one                                           !initial vector guess
      exitcode = 0
      Nsize = N
      if(N.eq.0) then 
         exitcode = 1
         if(ROOT) then 
            print*,'MPI_GEVP_InvIter_qp call with zero sized arrays.' 
         endif
 2000    IF(ROOT) THEN
            print *,'Abnormal termination on processor ',myrank
         ENDIF
         call MPI_ABORT(MPI_COMM_WORLD,exitcode,ierror)
         STOP ' '
      endif

      ALLOCATE(L0tab(3,0:nproc-1))
! Each process determines length "nn" of H, S partitions (vectors) 
      nn=0
      do j=myrank+1,n,nproc         !RR column ordering
         nn=nn+j
      enddo   !on j
! and sets up partition H = (H - Etrial*S), i.e., entry H is lost.
      do k = 1,nn
         H(k) = H(k) - Etrial * S(k)    !H = H - Etrial * S
      enddo
! Memory for implicit inverse L0 is distributed and must be dynamically
! ALLOCATED. Each process therefore needs to determine the size of its
! L0 partition (L0ArrSize) as well as the first (L0Start) and last
! (L0End) 
! rows of L0 in the partition.
      call IILTable(Nsize-1, nproc, M, L0tab)
      L0size_max=0
      do np = 0,nproc-1
         L0size_max = max(L0size_max, L0tab(3,np))
      enddo   !on np
! Allocate the space for the L0 partitions...
      if(ROOT)then
         ALLOCATE( L0(L0size_max+1) )
      else
         ALLOCATE( L0(L0tab(3,myrank)+1) )
      endif

! Now solve for one root (close to Etrial)...control will be by process 0.
      IF(ROOT)THEN
         print*,'Solving for one root using MPI_GEVP_InvIter_qp:'
         print *,' '
      ENDIF  !on ROOT 

! Calculate b = S*x0...x0 is the initial x vector everyone has.  
      do i = 1,N
         b(i) = 0       !all processes do b_i = S_i * x
      enddo
      M = 0
      do jq=myrank+1, N, nproc
         b(jq) = b(jq)+DOT_PRODUCT(x0(1:jq),S(M+1:M+jq))
         do i = 1,jq-1
            b(i) = b(i)+x0(jq)*S(M+i)
         enddo       
      enddo   !on jq     
      if(ROOT)then      !ROOT gathers the b_i to form vector b
         do i = 1,N
            t(i)=b(i)  
         enddo
         do km=1,nproc-1  
            call MPI_RECV(trecv, N, mytype, km, km, 
     *                    MPI_COMM_WORLD, istat, ierr) 
            do i = 1,N
               t(i) = t(i) + trecv(i)
            enddo
         enddo   !on k
         do i = 1,N
            b(i) = t(i)     !b = S*x0 is done
         enddo
      else
         call MPI_SEND(b,N,mytype,0,myrank,MPI_COMM_WORLD,ierr)  
      endif
  
! Calculate the L partitions of the L*U = H factorization. 
      if(ROOT.and.timing)time1 = MPI_WTIME()   !Proc 0 times the LU step. 
      if(ROOT)print*,'LU step for problem size', Nsize
      call LU_QP(Nsize, H, L0, tt, L0tab(1,myrank), L0tab(2,myrank)) 
      if(ROOT.and.timing)print '(a,f10.3,a)',  
     *      '  Time for L*U step:', MPI_WTIME()-time1,' seconds'

      if(ROOT) print*
! Inverse Iteration (II) iterations start here...note that the FULL 
! x-vector will be returned only to the ROOT instance of Ltdm1L_qp.
      timeiter = MPI_WTIME() 

      if(ROOT) then
         print*,'iterate using Ltdm1L_qp', N, Nsize
         print*
      endif
      kk=1      !counts number of inverse iterations done...
  300 continue  
      call Ltdm1L_qp(N, Nsize, L0, tt, b, x, L0tab(1,myrank), 
     *            L0tab(2,myrank))
      time1 = MPI_WTIME()
      if(ROOT)xb = DOT_PRODUCT( x(1:N),b(1:N) )
! Now calculate a new b = S*x...first send x to all processes.
      call MPI_BCAST(x, N, mytype, 0, MPI_COMM_WORLD, ierr)  !Bcast x 
      do i = 1,N
         b(i) = 0       !all processes do b_i = S_i * x
      enddo
      M = 0
      do jq=myrank+1, N, nproc
         b(jq) = b(jq)+DOT_PRODUCT(x(1:jq),S(M+1:M+jq))
         do i = 1,jq-1
            b(i) = b(i)+x(jq)*S(M+i)
         enddo
         M = M+jq
      enddo   !on jq     
      if(ROOT)then      !ROOT gathers the b_i to form vector b
         t(1:N)=b(1:N)  
         do km=1,nproc-1  
            call MPI_RECV(trecv, N, mytype, km, km, 
     *                    MPI_COMM_WORLD, istat, ierr) 
            do i = 1,N
               t(i) = t(i) + trecv(i)
            enddo
         enddo   !on k
         do i = 1,N
            b(i) = t(i)     !b = S*x is done
         enddo
      else
         call MPI_SEND(b,N,mytype,0,myrank,MPI_COMM_WORLD,ierr)  
      endif
  
! Proc 0 checks if we are done...
      if(ROOT)then         
         xsx = DOT_PRODUCT( x(1:N), b(1:N) )
         Eritz = Etrial + xb/xsx   
         Printstr = Spaces_qp('(f44.34)',Eritz,11)
         write(*,'(i3,1x,i5,a,a54,a,f12.6)') kk,n,' Eritz = ',
     &          Printstr,' Etrial =',etrial
         call flush(6)
         r = 1/sqrt(xsx)
         do i = 1,N
            b(i) = r*b(i)       !scale the b-vector
            x(i) = r*x(i)       !scale the x-vector also
         enddo
      endif   !on ROOT  
! Inform the other processes of the value of Eritz...
      call MPI_BCAST(Eritz, 1, mytype, 0, MPI_COMM_WORLD, ierr)  
      if(kk.gt.1)then
         if(abs((Eritz - Eold)/Eritz).lt.ConvEPS)go to 200
      endif  
  311 Eold = Eritz  
      if(kk.eq.Num_Iter)then
         if(ROOT)Print*,'Num_Iter:',Num_Iter,' Too many iterations.'
         goto 200  
      endif
      kk=kk+1  
      go to 300               !loop back for another iteration
  200 tout = MPI_WTIME()

      if(ALLOCATED(L0)) DEALLOCATE(L0)
      if(ALLOCATED(L0tab)) DEALLOCATE(L0tab)
      IF(ROOT)THEN
         Printstr = Spaces_qp('(f44.34)',Eritz,11)
         print*
         write(*,'(1x,i5,a,a54,a,f16.6)') N,' Eritz=',
     *            Printstr
         print*
         if(timing)then 
            print '(a,f10.2,a)',
     *     ' Iteration solution time =', time1-timeiter, ' seconds'
            print '(a,f10.2,a)',
     *     ' Secular equation solution time =', tout-tin, ' seconds'
	    call flush(6)
         endif
         print*
      ENDIF      ! IF(ROOT) 
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      RETURN  

      CONTAINS

         SUBROUTINE IILTable( N2, nproc, M, L0tab)
! Procedure for calculationg array L0tab used in partitioning the L0 
! implicit inverse array (used by MPI_GEVP_InvIter_* routines).
         integer(4), INTENT(IN) :: N2, nproc
         integer(4), INTENT(OUT) :: M
         integer(4), INTENT(OUT) :: L0tab(3,0:*) 
         integer(8) :: Msize, N1x, N2x
         integer(4) :: np, k
         N1x = 1
         N2x = N2
         Msize = (N2x*(N2x+1))/2 - (N1x*(N1x-1))/2 
         M = Msize/nproc
         np = 0
         L0tab(1,0) = 1
         L0tab(3,0) = 0
         do k = 1, N2
            if((L0tab(3,np)+k).le.M .or. np.eq.(nproc-1))then
               L0tab(3,np) = L0tab(3,np)+k
               L0tab(2,np) = k
            else
               np=np+1
               L0tab(1,np) = k
               L0tab(2,np) = k
               L0tab(3,np) = k
            endif
         enddo   !on k
         N1x=0
         do k=0,nproc-1
            N1x=N1x+L0tab(3,k)
         enddo
         if(Msize.ne.N1x) STOP 'IILTable:  Msize fatal error'
         RETURN
         END SUBROUTINE IILTable

      END SUBROUTINE MPI_GEVP_InvIter_qp
!
      SUBROUTINE Ltdm1L_qp(Nt, nx, L, D, b, x, Lstart, Lend) 
! This program carries out the calculation x = L(transpose)*D**(-1)*L*b  
! where L in its implicit form is spread across the several processes
! and b (in P0, which will Bcast it to the remaining processes) is a 
! starting vector of size Nt.  The L partitions and diagonal matrix D
! are generated by routine LU_QP and are of size nx, where nx >= Nt.  
! D belongs to P_lastproc, where lastproc=nproc-1. 

! To calculate L*b, P0 applies all its L(Lstart:Lend) to b(1:Nt), 
! then passes the result to P1, which in turn applies all its L_ks and 
! passes the result on to Proc 2 and so on. P(last) then applies D^{-1}, 
! then begins the process of applying L(transpose).  Final x-vector 
! winds up in P0.
      USE MPI
      implicit real(16)(a-h,o-z), integer(4)(i-n)
      integer(4), INTENT(IN)     :: Nt         !length of truncated expansion
      integer(4), INTENT(IN)     :: nx         !length of original expansion
      real(16), INTENT(IN)    :: L(*)
      real(16), INTENT(IN)    :: D(*)
      real(16), INTENT(IN)    :: b(*)
      real(16), INTENT(OUT)   :: x(*)
      integer(4), INTENT(IN)     :: Lstart 
      integer(4), INTENT(INOUT)  :: Lend


      integer istat(MPI_STATUS_SIZE)
      logical :: ROOT
      integer :: kflag=0
      SAVE :: mytype, myrank, nproc 
!
      if(kflag.eq.0)then
         kflag = 1
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
         call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)
         call MPI_TYPE_CONTIGUOUS(16, MPI_BYTE, mytype, ierror)
         call MPI_TYPE_COMMIT(mytype, ierror)
      endif
! 
      ROOT = myrank.eq.0
      n = nx        
      if(ROOT)then 
         x(1:Nt) = b(1:Nt)      
      else
         call MPI_Recv( x, Nt, mytype, myrank-1, myrank-1,
     *                  MPI_COMM_WORLD, istat, ierr )  
      endif
! All processes attempt L_{myrank)*x...
      k=1
      k1=Lstart
      do i = Lstart, min(Nt-1,Lend)         
         x(i+1) = x(i+1) + DOT_PRODUCT( x(1:i),L(k:k1) )
         k=k1+1
         k1=k+i
      enddo   !on i
      k1=k-1   !k,k1 values to be used coming back 
      k=k-i+1
  100 continue 
! All processes send x to next process except the last process..
      if(myrank.ne.(nproc-1)) then
         call MPI_Send( x, Nt, mytype, myrank+1, myrank, 
     *                  MPI_COMM_WORLD, ierr )
      else
         ! Last process does x(1:Nt) = D^{-1}*x(1:Nt)
         do i=1,Nt
            x(i)=x(i)*D(i)
         enddo   !on i
      endif
           
! We have done x = D^{-1}*L*x and all processes are synchronized...
! Finally, a left multiply of x by L(transpose), a coding sequence 
! similar to that for doing L*x, but in reverse L_i order.
      if(myrank.ne.(nproc-1))then   
         call MPI_Recv( x, Nt, mytype, myrank+1, myrank+1, 
     *                 MPI_COMM_WORLD, istat, ierr )
      endif
! Process myrank applies its transposed L terms to x... 
      do i = min(Nt-1,Lend), Lstart, -1
         do j = 1,i
            x(j) = x(j) + x(i+1)*L(k1+j-i)
         enddo
         k1=k-1
         k=k-i+1
      enddo   !on i
      continue
      if(myrank.ne.0) then      
         call MPI_Send( x, Nt, mytype, myrank-1, myrank, 
     *                  MPI_COMM_WORLD, ierr )
      endif   
      RETURN
      END SUBROUTINE Ltdm1L_qp
!
      SUBROUTINE LU_QP( Nsize, H, L, D, Lstart, Lend) 
      USE MPI
      implicit real(16) (a-h,o-z), integer(4)(i-n)

! Calculate the L partitions and diagonal matrix D (stored as D^-1)
! in the LU factorization of H

      integer(4), INTENT(IN)  :: Nsize
      real(16), INTENT(INOUT) :: H(*)
      real(16), INTENT(OUT)   :: L(*)
      real(16), INTENT(OUT)   :: D(*)
      integer(4), INTENT(IN) :: Lstart
      integer(4), INTENT(IN) :: Lend

      integer istat(MPI_STATUS_SIZE)
      real(16) :: t(Nsize), x(Nsize)
      logical :: ROOT
      integer :: kflag=0
      SAVE :: mytype, myrank, nproc 
!
      if(kflag.eq.0)then
         kflag = 1
         call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
         call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierror)
         call MPI_TYPE_CONTIGUOUS(16, MPI_BYTE, mytype, ierror)
         call MPI_TYPE_COMMIT(mytype, ierror)
      endif

      i1 = 1        !index in array H of the 1st element in a col.
      jL = 1        !indexs the packing of rows of the L matrix
      ROOT = myrank.eq.0
      if(ROOT)then
         i1 = 2
         do i = 1,Nsize
            D(i) = 0
         enddo
         D(1) = 1/H(1)
      endif
      nn = 0
      do k = myrank+1,Nsize,nproc
         nn = nn + k
      enddo   !on k

      do jcol = 2,Nsize
         k = mod(jcol-1,nproc)   !process that owns column jcol
         if(ROOT)then
            if(k.eq.0)then
               do j=1,jcol-1           !column  belongs to ROOT
                  x(j) = -H(i1)*D(j)
                  D(jcol) = D(jcol)+x(j)*H(i1)
                  i1 = i1+1
               enddo   !on j
               D(jcol) = ( 1/(D(jcol)+H(i1)) )   !invert diag. element
               i1 = i1+1
            else
               call MPI_Recv(t, jcol, mytype, k, jcol,
     *                       MPI_COMM_WORLD, istat, ierr)
! Column belongs to Proc k but Proc 0 will build L_jcol
               do j = 1,jcol-1
                  x(j) = -t(j)*D(j)
                  D(jcol) = D(jcol)+t(j)*x(j)
               enddo   !on j
               D(jcol) = ( 1/(D(jcol)+t(jcol)) )
            endif 
         elseif(myrank.eq.k)then
            call MPI_Send(H(i1), jcol, mytype, 0, jcol, 
     *                    MPI_COMM_WORLD, ierr)
            i1 = i1+jcol  
         endif

! All processes get a copy of L_jcol...
         call MPI_Bcast(x, jcol-1, mytype, 0, MPI_COMM_WORLD, ierr)
! and all Procs apply L_jcol to their H partition...
         nn = jcol-1
         jj = 1
         do j1a = myrank+1,Nsize,nproc
            if(jcol.le.j1a)exit
            jj = jj+j1a
         enddo   !on j1a        
         do j1 = j1a, Nsize, nproc
            H(jj+nn) = H(jj+nn) + DOT_PRODUCT(x(1:nn),H(jj:jj+nn-1))
            jj = jj + j1
         enddo   !on j1

! Then save L_jcol (vector x) in array L if Lstart <= jcol <= Lend...
         if(Lstart.le.nn .and. nn.le.Lend)then
            L(jL:jL+nn-1) = x(1:nn)
            jL = jL+nn       !set jL for next L-row
         endif
      enddo   !on jcol

! For convenience we send D off to lastproc (nproc-1)...
      if(ROOT)then
         call MPI_send(D, Nsize, mytype, nproc-1, Nsize,
     *                 MPI_COMM_WORLD, ierr)
      elseif(myrank.eq.(nproc-1))then
         call MPI_Recv(D, Nsize, mytype, 0, Nsize, MPI_COMM_WORLD,
     *                 istat, ierr)
      endif
      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      RETURN
      END SUBROUTINE LU_QP

      character(100) FUNCTION Spaces_qp(fmt, e, n)
      character(*), INTENT(IN) :: fmt    !(fw.d) only
      real(16), INTENT(IN) :: e
      integer(4), INTENT(IN) :: n   !number spaces inserted
! Routine for inserting spaces after every fourth digit of the 
! fractional part. Works for the fw.d format descriptor only.
! Typical usage...
!     character(50) :: string
!          . . .
!     string = spaces_qp('(f32.23)', value, 5)
!     print '(f10.5,a37)', zeta, string   !37 = 32+5 spaces
! To change the number digits between spaces set "length" accordingly.
! Or make length an optional argument. 
      character(100) :: str
      integer(4), parameter :: length = 4   !blank spacing
      write(str,fmt) e
      i = INDEX(str,'.')
      k = LEN_TRIM(str)
      do kk=1,n
         kkk = (length+1)*kk + i
         str(kkk+1:k+1) = str(kkk:k)
         str(kkk:kkk) = ' '
      enddo   !on kk
      spaces_qp = str(1:k+n)
      RETURN
      END FUNCTION Spaces_qp
