! ==============================================================================
!  FILE / MODULE: inverse.f90 (INVERSE)
!
!  PURPOSE & CONTEXT
!    Advanced numerical solver module for the waveguide QED simulation. It 
!    implements specialized matrix inversion and linear system solvers required 
!    to evaluate the variational equations of motion (EOMs).
!
!  CORE RESPONSIBILITIES
!    1. Direct Solvers  : Implements `directInverse` to solve the full system 
!                         of variational equations via standard LU decomposition.
!    2. Iterative Solvers: Implements "SuperInverse" routines using a stabilized 
!                         BiCGStab(l) Krylov subspace method with LU 
!                         preconditioning for high-dimensional bases (np > 12).
!    3. Data Management : Includes professional sorting (MergeSort) and 
!                         diagnostic printing routines for complex tensors.
!
!  BUILD / DEPENDENCIES
!    * Requires BLAS/LAPACK (ZGETRF, ZTRSM, ZLASWP).
!
! ==============================================================================
!
MODULE INVERSE

  USE consts
  USE lapackmodule
  USE typedefs, only : cx => c_type, rl=> r_type
  IMPLICIT none

  INTERFACE printArray
    module procedure printMatrix
    module procedure printVector
!> -------------------------------------------------------------------------
    !> INTERFACE: printArray
    !> -------------------------------------------------------------------------
    !> Purpose / context:
    !>   Generic interface for diagnostic output. Allows the developer to call 
    !>   `printArray` for either 1D complex vectors or 2D complex matrices 
    !>   during debugging of the variational parameters.
    !> -------------------------------------------------------------------------
  END INTERFACE


  CONTAINS

!> -------------------------------------------------------------------------
  !> SUBROUTINE: directInverse
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   A standard direct solver for the multi-polaron variational system. 
  !>   It maps the 3D tensor representation of the overlaps `A(n,n,n)` into a 
  !>   flattened 2D matrix `A2D(n^2, n^2)` and solves the resulting linear 
  !>   system for the parameter derivatives using LU decomposition (`SolveEq_c`).
  !> Arguments:
  !>   - n   : Number of polarons in the basis.
  !>   - A   : Input 3D overlap tensor.
  !>   - rhs : Right-hand side driving vector.
  !>   - sol : Output solution vector (parameter derivatives).
  !>
  subroutine directInverse(n,A,sol,rhs)
    complex(cx),intent(in)     :: A(n,n,n)
    complex(cx),intent(in)     :: rhs(n**2)
    complex(cx),intent(in out) :: sol(n**2)
    complex(cx)           :: A2D(n**2, n**2)
    integer               :: n,nn, i,j,k  !-- n=block size, nn=array size


    !-- Build Matrix
    A2D=0._cx
    do i=1, n
      do j=1, n
        do k=1, n
          A2D((j-1)*n+i,(j-1)*n+k)=A(i,j,k)
        end do
        A2D((j-1)*n+i,(i-1)*n+j)=1.0_cx
      end do
    end do
    ! -- resolve
    sol=rhs
    Call solveEq_c(A2D, sol)
  end subroutine directInverse

!> -------------------------------------------------------------------------
  !> SUBROUTINE: superInverse_f / superInverse_h
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   High-performance iterative solver for the bosonic "f" or "h" manifolds. 
  !>   Designed for large basis sizes where direct O(N^6) inversion is costly. 
  !>   It utilizes an LU-preconditioned BiCGStab(2) algorithm to converge on 
  !>   the solution. If the iterative method fails to converge within tolerance, 
  !>   it automatically falls back to `directInverse` for stability.
  !> Arguments:
  !>   - n   : Number of polarons.
  !>   - A   : Overlap tensor.
  !>   - sol : Current solution / seed.
  !>   - rhs : Driving terms.
  !>
  subroutine superInverse_f(n,A,sol,rhs)
    complex(cx),intent(in)     :: A(n,n,n)
    complex(cx),intent(in)     :: rhs(n**2)
    complex(cx),intent(in out) :: sol(n**2)
    integer               :: n,nn   !-- n=block size, nn=array size
    integer               ::i,j,k
    !--Krylov parameters
    logical                      :: print_resid, precondActive, nonZeroSeed
    integer                      :: lKry, info, mxmatvec
    real(rl)                     :: toler
    complex(cx)                  :: work(n**2,9)
    complex(cx),allocatable,save :: seed(:)
    ! -- time prec parameters
    complex(cx)                    :: A2D(n**2,n**2)
    complex(cx), allocatable,save  :: Lo(:,:), Up(:,:)
    integer, allocatable,save      :: Perm(:)
    complex(cx)                    :: alpha
    !--ALLOCATION
    nn=n**2
    info=0
    nonZeroSeed=.true.
    ! -- initialise time preconditionner
    if ((ALLOCATED(Lo) .eqv. .true.) .and. (size(Lo,1) /= nn)) then
      deallocate(Lo)
      deallocate(Up)
      deallocate(Perm)
      deallocate(seed)
    end if
    if (ALLOCATED(Lo) .eqv. .false.) THEN
      !print*, 'begin allocation ...'
      allocate(Lo(nn,nn))
      allocate(Up(nn,nn))
      allocate(Perm(nn))
      allocate(seed(nn))
      Lo=0._cx
      Up=0._cx
      seed=0._cx
      nonZeroSeed=.false.
      do i=1, nn
        perm(i)=0
      end do
      A2D=0._cx
      do i=1, n
        do j=1, n
          do k=1, n
            A2D((j-1)*n+i,(j-1)*n+k)=A(i,j,k)
          end do
          A2D((j-1)*n+i,(i-1)*n+j)=1.0_cx
        end do
      end do
      call ZGETRF(nn,nn,A2D,nn,Perm,info)
      do i=1, nn
        do j=i+1, nn
          Lo(j,i)=A2D(j,i)
          Up(i,j)=A2D(i,j)
        end do
        Lo(i,i)=1.0_cx
        Up(i,i)=A2D(i,i)
      end do
      !print*, 'LU decomposed. backwrd sub ...'
      sol=rhs
      alpha=1.0_cx
      call ZLASWP(1,sol,nn,1,nn,Perm,1)
      call ZTRSM('l','l','n','u',nn,1,alpha,Lo,nn,sol,nn)
      call ZTRSM('l','u','n','n',nn,1,alpha,Up,nn,sol,nn)
      seed=sol
      !print*, 'PREC initialized.'
      return
    end if
    !--call krylovroutine
    precondActive=.true.
    print_resid=.false.
    sol=seed
    lKry=2
    mxmatvec=n**2
    toler=0.000000000001
    CALL ZBCG2(print_resid,lKry,nn,sol,nonZeroSeed,rhs,n,A,toler, &
               mxmatvec,work,precondActive,Lo,Up,Perm,info)
    seed=sol
    if (isnan(toler).or.(mxmatvec>=n**2-1)) Then
      call directInverse(n,A,sol,rhs)
      deallocate(Lo)
      deallocate(Up)
      deallocate(Perm)
      deallocate(seed)
    elseif (mxmatvec > n) then
      deallocate(Lo)
      deallocate(Up)
      deallocate(Perm)
      deallocate(seed)
    end if
  end subroutine

  subroutine superInverse_h(n,A,sol,rhs)
    complex(cx),intent(in)     :: A(n,n,n)
    complex(cx),intent(in)     :: rhs(n**2)
    complex(cx),intent(in out) :: sol(n**2)
    integer               :: n,nn   !-- n=block size, nn=array size
    integer               ::i,j,k
    !--Krylov parameters
    logical                      :: print_resid, precondActive, nonZeroSeed
    integer                      :: lKry, info, mxmatvec
    real(rl)                     :: toler
    complex(cx)                  :: work(n**2,9)
    complex(cx),allocatable,save :: seed(:)
    ! -- time prec parameters
    complex(cx)                    :: A2D(n**2,n**2)
    complex(cx), allocatable,save  :: Lo(:,:), Up(:,:)
    integer, allocatable,save      :: Perm(:)
    complex(cx)                    :: alpha
    !--ALLOCATION
    nn=n**2
    info=0
    nonZeroSeed=.true.
    ! -- initialise time preconditionner
    if ((ALLOCATED(Lo) .eqv. .true.) .and. (size(Lo,1) /= nn)) then
      deallocate(Lo)
      deallocate(Up)
      deallocate(Perm)
      deallocate(seed)
    end if
    if (ALLOCATED(Lo) .eqv. .false.) THEN
      allocate(Lo(nn,nn))
      allocate(Up(nn,nn))
      allocate(Perm(nn))
      allocate(seed(nn))
      Lo=0._cx
      Up=0._cx
      seed=0._cx
      nonZeroSeed=.false.
      do i=1, nn
        perm(i)=0
      end do
      A2D=0._cx
      do i=1, n
        do j=1, n
          do k=1, n
            A2D((j-1)*n+i,(j-1)*n+k)=A(i,j,k)
          end do
          A2D((j-1)*n+i,(i-1)*n+j)=1.0_cx
        end do
      end do
      call ZGETRF(nn,nn,A2D,nn,Perm,info)
      do i=1, nn
        do j=i+1, nn
          Lo(j,i)=A2D(j,i)
          Up(i,j)=A2D(i,j)
        end do
        Lo(i,i)=1.0_cx
        Up(i,i)=A2D(i,i)
      end do
      ! print*, 'LU decomposed. backwrd sub ...'
      sol=rhs
      alpha=1.0_cx
      call ZLASWP(1,sol,nn,1,nn,Perm,1)
      call ZTRSM('l','l','n','u',nn,1,alpha,Lo,nn,sol,nn)
      call ZTRSM('l','u','n','n',nn,1,alpha,Up,nn,sol,nn)
      seed=sol
      return
    end if
    !--call krylovroutine
    precondActive=.true.
    print_resid=.false.
    sol=seed
    lKry=2
    mxmatvec=n**2
    toler=0.000000000001
    CALL ZBCG2(print_resid,lKry,nn,sol,nonZeroSeed,rhs,n,A,toler, &
               mxmatvec,work,precondActive,Lo,Up,Perm,info)
    seed=sol
    if (isnan(toler).or.(mxmatvec>=n**2-1)) Then
      ! print*, ' !! --------- Krylov failure! lauch direct inverse.'
      call directInverse(n,A,sol,rhs)
      deallocate(Lo)
      deallocate(Up)
      deallocate(Perm)
      deallocate(seed)
    elseif (mxmatvec > n) then
      deallocate(Lo)
      deallocate(Up)
      deallocate(Perm)
      deallocate(seed)
    end if
  end subroutine

!> -------------------------------------------------------------------------
  !> SUBROUTINE: ZBCG2
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Core implementation of the Bi-Conjugate Gradient Stabilized method 
  !>   of order l=2. This routine minimizes the residual of the non-symmetric 
  !>   linear systems generated by the multi-polaron ansatz. Includes 
  !>   optimizations for reliable residual updates in finite precision.
  !>
  subroutine ZBCG2 (print_resid,l,n,x,nonzero_x,rhs,m,A,toler, &
                    mxmatvec,work,precondActive,Lo,Up,Perm,&
                    info)

    ! subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
    !                   mxmatvec,work,info)
    !
    ! Improved "vanilla" BiCGStab(2) iterative method
    !
    ! Copyright (c) 2001 by M.A.Botchev (http://www.math.utwente.nl/~botchev/),
    !                       University of Twente
    ! Permission to copy all or part of this work is granted,
    ! provided that the copies are not made or distributed
    ! for resale, and that the copyright notice and this
    ! notice are retained.
    !
    ! This is the "vanilla" version of BiCGstab(\ell) as described
    ! in PhD thesis of D.R.Fokkema, Chapter 3 (also available as
    ! Preprint 976, Dept. of Mathematics, Utrecht University, URL
    ! http://www.math.uu.nl/publications/).  It includes two enhancements
    ! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
    ! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence
    !    properties of BiCGstab methods in finite precision arithmetic",
    !    Numerical Algorithms, 10, 1995, pp.203-223
    ! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
    !    hybrid BiCG methods", Computing, 56, 1996, pp.141-163
    !
    ! {{ This code based on original work of D.R.Fokkema:
    !
    ! subroutine zbistbl v1.1 1998
    ! Copyright (c) 1995-1998 by D.R. Fokkema.
    ! Permission to copy all or part of this work is granted,
    ! provided that the copies are not made or distributed
    ! for resale, and that the copyright notice and this
    ! notice are retained.
    !
    ! }}
    !
    ! Your bug reports, comments, etc. are welcome:
    ! m.a.botchev@math.utwente.nl
    !
    ! ------------------------------
    ! Description of the parameters:
    ! ------------------------------
    !
    ! print_resid (input) LOGICAL. If print_resid=.true. the number of
    !            matrix-vector multiplications done so far and residual norm will
    !            be printed to the standard output each iteration
    !
    ! l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
    !            in this simple version it is required that l <= 2
    !            l=2 is often useful for systems with nonsymmetric matrices
    !
    ! n          (input) INTEGER size of the linear system to solve
    !
    ! x          (input/output) COMPLEX(cx) array dimension n
    !            initial guess on input, solution on output
    !
    ! rhs        (input) COMPLEX(cx) array dimension n
    !            the right-hand side (r.h.s.) vector
    !
    ! matvec     (input) EXTERNAL name of matrix vector subroutine
    !            to deliver y:=A*x by CALL matvec(n,x,y,A)
    !
    ! A          the matrix called by matvec
    !
    ! nonzero_x  (input) LOGICAL tells
    !            BiCGstab(\ell) if the initial guess x is zero or not.
    !            If nonzero_x is .FALSE., initial residual equals r.h.s. vector
    !            and one MATVEC call is avoided
    !
    ! toler      (input/output) DOUBLE PRECISION tolerance: the iterations are
    !            stopped as soon as || residual ||/|| initial residual|| <= toler,
    !            the norm is Euclidean.  On output, if info>=0, the value of
    !            toler is set to the actually achieved residual reduction
    !
    ! mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix
    !            vector multiplications allowed to be done.  On output:
    !            if info>=0, mxmatvec is set to the actual number of matrix
    !            vector multiplications done
    !
    ! work       (workspace) COMPLEX(cx) array of dimension (n,2*l+5)
    !
    ! info       (output) INTEGER.  info = 0 in case of succesful computations
    !            and
    !            info = -m (<0) - means paramater number m has an illegal value
    !            info = 1 - means no convergence achieved (stopping criterion
    !            is not fulfilled)
    !            info = 2 - means breakdown of the algorithm (taking a larger
    !            value of parameter l usually helps)
    !
    ! WARNING: If the iterations are ended normally (info=0 or info=1),
    ! the true residual norm is computed and returned as an output value
    ! of the parameter toler.  The true residual norm can be slightly larger
    ! than the projected residual norm used by the algorithm to stop the
    ! iterations.  It may thus happen that on output info=0 but the value
    ! of toler is (slightly) larger than tolerance prescribed on input.
    ! ----------------------------------------------------------
    implicit none

    ! Parameters:
    logical,   intent(in)   :: print_resid,nonzero_x
    integer,   intent(in)   :: l, n, m
    integer,   intent(inout):: mxmatvec
    integer,   intent(out)  :: info
    complex(cx),intent(inout):: x(n)
    complex(cx),intent(in)   :: rhs(n)
    complex(cx),intent(in)  :: A(m,m,m)
    complex(cx)             :: Lo(n,n), Up(n,n)
    integer                 :: Perm(n)
    real(rl),  intent(inout):: toler
    complex(cx),intent(out)  :: work(n,3+2*(l+1))
    logical                  :: precondActive

    ! Local variables:
    complex(cx) :: matrix_z(l+1,l+1),y0(l+1),yl(l+1),zy0(l+1),zyl(l+1)
    logical    :: rcmp, xpdt
    integer    :: i, j, k, nmatvec
    complex(cx) :: alpha, beta, omega, rho0, rho1, sigma
    complex(cx) :: varrho, hatgamma
    real(rl)     :: rnrm0, rnrm
    real(rl)     :: mxnrmx, mxnrmr
    complex(cx) :: kappa0, kappal

    ! Aliases for the parts of the work array:
    integer          :: rr, r, u, xp, bp

    ! Constants:
    real(rl),    parameter :: delta = 1d-2
    complex(cx),parameter :: zzero = (0d0,0d0), zone = (1d0,0d0)

    ! Functions:
    !real(8)     :: dnorm2_bcg
    !complex(cx) :: zdot_bcg

    info = 0

    if (l<1 .or. l>2) info = -2
    if (n<1) info = -3
    if (toler<=0d0) info = -9
    if (mxmatvec<0) info = -10

    rr = 1
    r = rr+1
    u = r+(l+1)
    xp = u+(l+1)
    bp = xp+1

    if (info/=0) return


    ! Initialize first residual




    if (nonzero_x) then
       call blockmatvec (m, x, work(1:n,r), A)
       work(1:n,r) = rhs - work(1:n,r)
       nmatvec = 1
    else
       work(1:n,r) = rhs
       nmatvec = 0
    end if
    call precond (n, work(1:n,r),precondActive,Lo,Up,Perm)


    ! Initialize iteration loop

    work(1:n,rr) = work(1:n,r)
    work(1:n,bp) = work(1:n,r)
    work(1:n,xp) = x
    x = zzero



    rnrm0 = dnorm2_bcg (n, work(1:n,r))
    rnrm = rnrm0

    mxnrmx = rnrm0
    mxnrmr = rnrm0
    rcmp = .false.
    xpdt = .false.

    alpha = zzero
    omega = zone
    sigma = zone
    rho0  = zone

    ! Iterate
    do while (rnrm > toler*rnrm0 .and. nmatvec < mxmatvec)

    ! =====================
    ! The BiCG part ---
    ! =====================

       rho0 = -omega*rho0

       do k=1,l
          rho1 = zdot_bcg (n, work(1:n,rr), work(1:n,r+k-1))
          if (rho0.eq.zzero) then
             info = 2
             toler = rnrm/rnrm0
             mxmatvec = nmatvec
             return
          endif
          beta = alpha*(rho1/rho0)
          rho0 = rho1
          do j=0,k-1
             work(1:n,u+j) = work(1:n,r+j) - beta*work(1:n,u+j)
          enddo
          call blockmatvec(m,work(1:n,u+k-1), work(1:n,u+k), A)
          call precond(n, work(1:n,u+k),precondActive,Lo,Up,Perm)
          nmatvec = nmatvec+1

          sigma = zdot_bcg (n, work(1:n,rr), work(1:n,u+k))
          if (sigma.eq.zzero) then
             info = 2
             toler = rnrm/rnrm0
             mxmatvec = nmatvec
             return
          endif
          alpha = rho1/sigma
          x(1:n) = alpha*work(1:n,u) + x(1:n)
          do j=0,k-1
             work(1:n,r+j) = -alpha*work(1:n,u+j+1) + work(1:n,r+j)
          enddo
          call blockmatvec (m, work(1:n,r+k-1), work(1:n,r+k), A)
          call precond (n, work(1:n,r+k),precondActive,Lo,Up,Perm)
          nmatvec = nmatvec+1
          rnrm = dnorm2_bcg (n, work(1,r))
          mxnrmx = max (mxnrmx, rnrm)
          mxnrmr = max (mxnrmr, rnrm)
       enddo

    ! ==================================
    ! The convex polynomial part ---
    ! ==================================

    !  --- Z = R'R

       do i=1,l+1
          do j=1,i
             matrix_z(i,j) = conjg(zdot_bcg( n, work(1:n,r+j-1), work(1:n,r+i-1) ))
          end do
       end do

    !  lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
       do j=2,l+1
          matrix_z(1:j-1,j) = conjg( matrix_z(j,1:j-1) )
       end do

    !  small vectors y0 and yl

    y0(1) = -zone
    y0(2) =      ( matrix_z(2,1) / matrix_z(2,2) )   ! works only for l=2
    y0(l+1) = zzero

    yl(1) = zzero
    yl(2) =      ( matrix_z(2,3) / matrix_z(2,2) )   ! works only for l=2
    yl(l+1) = -zone

    !  --- Convex combination

    ! compute Z*y0 and Z*yl
    zy0 = zzero
    zyl = zzero
    do j=1,l+1
       zy0 = zy0 + matrix_z(:,j)*y0(j)
       zyl = zyl + matrix_z(:,j)*yl(j)
    end do

    kappa0 = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )
    kappal = sqrt( abs(zdot_bcg(l+1,yl,zyl)) )

    varrho = zdot_bcg(l+1,yl,zy0) / (kappa0*kappal)

    hatgamma = varrho/abs(varrho) * max( abs(varrho),7d-1 )

    y0 = y0 - (hatgamma*kappa0/kappal)*yl


    !  --- Update

    omega = y0(l+1)

    do j=1,l
       work(1:n,u) = work(1:n,u) - y0(j+1)*work(1:n,u+j)
       x(1:n)      = x(1:n)      + y0(j+1)*work(1:n,r+j-1)
       work(1:n,r) = work(1:n,r) - y0(j+1)*work(1:n,r+j)
    enddo

    ! y0 has changed; compute Z*y0 once more
    zy0 = zzero
    do j=1,l+1
       zy0 = zy0 + matrix_z(:,j)*y0(j)
    end do

    rnrm = sqrt( abs(zdot_bcg(l+1,y0,zy0)) )

    ! ================================
    ! The reliable update part ---
    ! ================================

    mxnrmx = max (mxnrmx, rnrm)
    mxnrmr = max (mxnrmr, rnrm)
    xpdt =  (rnrm < delta*rnrm0  .and. rnrm0 < mxnrmx)
    rcmp = ((rnrm < delta*mxnrmr .and. rnrm0 < mxnrmr) .or. xpdt)
    if (rcmp) then
       call blockmatvec (m, x, work(1:n,r),A)
       call precond (n, work(1:n,r),precondActive,Lo,Up,Perm)
       nmatvec = nmatvec + 1
       work(1:n,r) =  work(1:n,bp) - work(1:n,r)
       mxnrmr = rnrm
       if (xpdt) then

          work(1:n,xp) = x(1:n) + work(1:n,xp)
          x = zzero
          work(1:n,bp) = work(1:n,r)

          mxnrmx = rnrm
       endif
    endif

    if (print_resid) print *,nmatvec,' ',rnrm

    enddo

    ! =========================
    ! End of iterations ---
    ! =========================

    x(1:n) = x(1:n) + work(1:n,xp)

    if (rnrm>toler*rnrm0) info = 1

    ! compute the true residual:

    ! --------------------- One matvec can be saved by commenting out this:
    ! call blockmatvec (m, x, work(1:n,r), A )
    ! work(1:n,r) = rhs(1:n) - work(1:n,r)
    ! call precond (n, work(1:n,r),precondActive,param,PREC,Lo,Up,Perm)
    ! rnrm = dnorm2_bcg(n,work(1:n,r))
    ! nmatvec = nmatvec+1
    ! --------------------- One matvec can be saved by commenting out this^
    !print *,nmatvec,' ',rnrm

    toler = rnrm
    mxmatvec = nmatvec
  end subroutine ZBCG2


!> -------------------------------------------------------------------------
  !> FUNCTION: zdot_bcg
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the complex inner product of two vectors, defined as the 
  !>   sum of `conjg(zx) * zy`. This is the fundamental operation for 
  !>   calculating projections and orthogonality within the Krylov subspace 
  !>   methods.
  !> Arguments:
  !>   - n      : Vector dimension.
  !>   - zx, zy : Complex vectors.
  !> Return:
  !>   - complex(cx) : The resulting complex scalar inner product.
  !>
  function zdot_bcg(n,zx,zy)

    ! complex inner product function

    implicit none
    integer,       intent(in):: n
    complex(cx),intent(in):: zx(n),zy(n)
    complex(cx)           :: zdot_bcg

    zdot_bcg = sum( conjg(zx) * zy )
  end function zdot_bcg


!> -------------------------------------------------------------------------
  !> FUNCTION: dnorm2_bcg
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the standard Euclidean L2 norm of a complex vector. Used 
  !>   internally by the BiCGStab(2) algorithm to monitor the residual 
  !>   magnitude and determine numerical convergence.
  !> Arguments:
  !>   - n  : Vector dimension.
  !>   - zx : The complex vector to be evaluated.
  !> Return:
  !>   - real(8) : The scalar L2 norm.
  !>
  function dnorm2_bcg(n,zx)

    ! l2 norm function

    implicit none
    integer,       intent(in):: n
    complex(cx),intent(in):: zx(n)
    !complex(cx),external  :: zdot_bcg
    real(8)           :: dnorm2_bcg

    dnorm2_bcg = sqrt( abs( zdot_bcg(n, zx, zx) ) )
  end function dnorm2_bcg

!> -------------------------------------------------------------------------
  !> SUBROUTINE: matvec
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   A generic wrapper for complex matrix-vector multiplication (y = A * x). 
  !>   Utilizes the Fortran intrinsic `MATMUL` for standard 2D dense matrices 
  !>   used within the iterative solver routines.
  !> Arguments:
  !>   - n : Dimension of the system.
  !>   - x : Input complex vector.
  !>   - y : Output complex vector.
  !>   - A : Input complex coefficient matrix.
  !>
  SUBROUTINE matvec(n, x, y, A)
    integer               :: n
    complex(cx)           :: x(:), y(:)
    complex(cx)           :: A(:,:)

    y=MATMUL(A,x)

    RETURN
  end subroutine

!> -------------------------------------------------------------------------
  !> SUBROUTINE: blockMatvec
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the specialized matrix-vector product required for the 
  !>   multi-polaron variational derivatives. It treats the flattened vector 
  !>   `x` as a block-structure and applies the 3D tensor `Alpha` to evaluate 
  !>   the coupling between different coherent state components across 
  !>   the polaron basis.
  !> Arguments:
  !>   - n     : Block size (number of polarons).
  !>   - x     : Input flattened complex vector of size n^2.
  !>   - y     : Output flattened complex vector of size n^2.
  !>   - Alpha : 3D tensor representing the polaron overlap couplings.
  !>
  subroutine blockMatvec(n, x, y, Alpha)
    integer                  :: n
    integer                  :: i,j

    complex(cx)              :: x(n**2)
    complex(cx)              :: y(n**2)
    complex(cx)              :: tempMat(n,n)
    complex(cx)              :: tempVec(n)
    complex(cx)              :: xBlock(n,n), yBlock(n,n)
    complex(cx)              :: Alpha(n,n,n)

    do i=1, n
      do j=1, n
        xBlock(i,j)=x((i-1)*n+j)
        yBlock(j,i)= xBlock(i,j)
      end do
    end do
    do i=1, n
      tempMat= Alpha(:,i,:)
      tempVec= xBlock(i,:)
      yBlock(i,:)= matmul(tempMat, tempVec)+yblock(i,:)
    end do
    do i=1, n
      do j=1, n
        y((i-1)*n+j)=yBlock(i,j)
      end do
    end do
  end subroutine blockMatvec


!> -------------------------------------------------------------------------
  !> SUBROUTINE: precond
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Applies an LU-based preconditioner to a vector `y` during the iterative 
  !>   BiCGStab solver. It performs a forward and backward substitution using 
  !>   the cached `Lo` and `Up` triangular matrices and the pivot permutation 
  !>   array `Perm`. This accelerates convergence by transforming the 
  !>   variational system into a better-conditioned identity-proximal form.
  !> Arguments:
  !>   - n      : Vector dimension.
  !>   - y      : The vector to be preconditioned (modified in-place).
  !>   - active : Boolean flag to toggle preconditioning.
  !>   - Lo, Up : Cached lower and upper triangular LU factors.
  !>   - Perm   : Cached pivot indices for row swapping.
  !>
  SUBROUTINE precond (n,y,active,Lo,Up,Perm)
    integer                  :: n
    complex(cx)              :: y(n)
    LOGICAL                  :: active
    ! time prec parameters
    complex(cx)              :: Lo(n,n), Up(n,n)
    integer                  :: Perm(n)
    ! intern variables
    complex(cx)   :: seed(n)
    complex(cx)   :: alpha

    if (active .eqv. .FALSE.) THEN
      return
    end if
    alpha=1.0_cx
    call ZLASWP(1,y,n,1,n,Perm,1)
    call ZTRSM('l','l','n','u',n,1,alpha,Lo,n,y,n)
    call ZTRSM('l','u','n','n',n,1,alpha,Up,n,y,n)
  end subroutine

  !-- sorting routines

!> -------------------------------------------------------------------------
  !> SUBROUTINE: Merge
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   The merging phase of the MergeSort algorithm. Combines two sorted 
  !>   sub-arrays (`A` and `B`) into a single sorted array (`C`). Crucially, 
  !>   it simultaneously reorders an associated complex array (`indx`), 
  !>   ensuring that the variational parameters remain linked to their 
  !>   corresponding sorted physical indices.
  !> Arguments:
  !>   - NA, NB, NC : Dimensions of sub-arrays and destination array.
  !>   - A, B       : Input sorted integer sub-arrays.
  !>   - C          : Output unified sorted integer array.
  !>   - indx, indxT: Complex index arrays associated with A and B.
  !>   - indxC      : Unified complex index array associated with C.
  !>
  subroutine Merge(A,NA,B,NB,C,NC,indx,indxT,indxC)

   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
   integer, intent(in)     :: B(NB)
   integer, intent(in out) :: C(NC)
   complex(cx), intent(in out) :: indx(NA)
   complex(cx), intent(in)     :: indxT(NB)
   complex(cx), intent(in out) :: indxC(NC)

   integer :: I,J,K

   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
         indxC(K) = indx(I)
         I = I+1
      else
         C(K) = B(J)
         indxC(K)=indxT(J)
         J = J+1
      endif
      K = K + 1
   enddo
   do while (I <= NA)
      C(K) = A(I)
      indxC(K) = indx(I)
      I = I + 1
      K = K + 1
   enddo
   return
 end subroutine merge

!> -------------------------------------------------------------------------
  !> SUBROUTINE: Sort / MergeSort
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Implements a stable recursive MergeSort algorithm. Used to organize 
  !>   variational parameters or grid indices based on magnitude. This is 
  !>   essential for pruning the basis when two polarons become nearly 
  !>   identical (high overlap).
  !>
 recursive subroutine MergeSort(A,N,T,indx,indxT)

   integer, intent(in) :: N
   integer, dimension(N), intent(in out) :: A
   integer, dimension((N+1)/2), intent (out) :: T
   complex(cx), dimension(N), intent (in out)  :: indx
   complex(cx), dimension((N+1)/2), intent (out) :: indxT

   integer :: NA,NB,V
   complex(cx) :: Vcx

   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then
         V = A(1)
         A(1) = A(2)
         A(2) = V
         Vcx=indx(1)
         indx(1)=indx(2)
         indx(2)=Vcx
      endif
      return
   endif
   NA=(N+1)/2
   NB=N-NA

   call MergeSort(A,NA,T, indx,indxT)
   call MergeSort(A(NA+1),NB,T,indx(NA+1),indxT)

   if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      indxT(1:NA)=indx(1:NA)
      call Merge(T,NA,A(NA+1),NB,A,N,indxT,indx(NA+1),indx)
   endif
   return
 end subroutine MergeSort

 subroutine Sort(n,A,indx)
   integer, intent (in)        :: n
   integer, intent(in out)     :: A(n)
   complex(cx), intent (in out):: indx(n)
   integer                     :: i
   integer, allocatable        :: T(:)
   complex(cx), allocatable    :: indxT(:)

   allocate(T((n+1)/2))
   allocate(indxT((n+1)/2))
   call MergeSort(A,n,T,indx,indxT)
 end subroutine Sort

  ! -- math routines

!> -------------------------------------------------------------------------
  !> FUNCTION: kroneckerDelta
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Mathematical utility evaluating the discrete Kronecker delta function 
  !>   δ_{a,b}. Returns 1 if indices match and 0 otherwise.
  !>
  FUNCTION kroneckerDelta(a,b)

 		integer,intent(in)   ::  a,b
 		real(rl)				   ::  kroneckerDelta

 		if (a==b) then
 		  kroneckerDelta = 1
 		else
 		  kroneckerDelta = 0
 		end if
  END FUNCTION

!> -------------------------------------------------------------------------
  !> SUBROUTINE: printMatrix / printVector
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Formatted diagnostic output for complex-valued arrays. It prints a 
  !>   readable grid of the matrix elements, representing unity as "1" and 
  !>   near-zero values as "." to allow for quick visual inspection of 
  !>   orthogonality and numerical stability.
  !>
  SUBROUTINE printMatrix(A)

    complex(cx)                      :: A(:,:)
    integer                          :: i,j
    print*, '______'
    print*, ' '
    do i=1, size(A,1)
      do j=1, size(A,2)
        if ((ABS(real(A(i,j))-1) .LT. 0.00001) .AND. (abs(aimag(A(i,j))) .LT. 0.00001)) THEN
          write(*,'(a7, a1, a6, a2)', advance='no') '       ', '1', '      ', ' |'
        elseif ((ABS(real(A(i,j))) .GT. 0.00001) .OR. (ABS(aimag(A(i,j))) .GT. 0.00001)) THEN
            write(*,'(f6.3,a2, f6.3, a2)', advance='no') real(A(i,j)),'+i', aimag(A(i,j)), ' |'
        else
          write(*,'(a7, a1, a6, a2)', advance='no') '       ', '.', '      ', ' |'
        end if
      end do
      write(*,*)
    end do
    print*, '______'
  END SUBROUTINE

  SUBROUTINE printVector(A)
    complex(cx)                      :: A(:)
    integer                          :: i,j
    print*, '______'
    print*, ' '
    do i=1, size(A,1)
        write(*,'(f9.6,a2, f9.6, a2)', advance='no') real(A(i)),'+i', aimag(A(i)), ' |'
    end do
    print*, ' '
    print*, '______'
  END SUBROUTINE

END MODULE INVERSE
