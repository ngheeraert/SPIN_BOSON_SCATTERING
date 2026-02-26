! ==============================================================================
!  FILE / MODULE: lapackmodule.f90 (lapackmodule)
!
!  PURPOSE & CONTEXT
!    Provides high-level Fortran wrappers for standard BLAS and LAPACK routines 
!    used throughout the waveguide QED simulation. This ensures the dense linear 
!    algebra operations required by the variational equations of motion (EOMs) 
!    are executed efficiently.
!
!  CORE RESPONSIBILITIES
!    1. Matrix Inversion : Wrappers for Cholesky (`ZPOTRF`), LU (`ZGETRF`), 
!                          and Bunch-Kaufman (`ZHETRF`) factorizations to invert 
!                          Hermitian, Positive-Definite, and General matrices.
!    2. Linear Solvers   : Direct solvers for AX = B systems (`DGESV`, `ZGESV`).
!    3. Multiplication   : Operator overloading for fast BLAS matrix 
!                          multiplication (`ZGEMM`).
!
! ==============================================================================
!
MODULE lapackmodule

  USE typedefs, only : cx => c_type, rl=> r_type
  USE consts

  implicit none

  INTERFACE OPERATOR(.matprod.)
    module procedure matmultiply_c
!> -------------------------------------------------------------------------
    !> INTERFACE: OPERATOR(.matprod.)
    !> -------------------------------------------------------------------------
    !> Purpose / context:
    !>   Overloads the custom `.matprod.` operator to seamlessly compute the 
    !>   product of two complex matrices using the highly optimized BLAS `zgemm` 
    !>   routine under the hood, rather than standard unoptimized Fortran `matmul`.
    !> -------------------------------------------------------------------------
  END INTERFACE OPERATOR(.matprod.)

CONTAINS

!> -------------------------------------------------------------------------
  !> SUBROUTINE: InvertHPD
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Inverts a Hermitian Positive Definite (HPD) matrix in-place. 
  !>   It utilizes the LAPACK Cholesky factorization routine (`ZPOTRF`) 
  !>   followed by the specific inversion routine (`ZPOTRI`). Since LAPACK 
  !>   only returns the upper triangular part, this routine explicitly reconstructs 
  !>   the full matrix by mirroring the upper half into the lower half.
  !> Arguments:
  !>   - A : The complex HPD matrix to be inverted (modified in-place to its inverse).
  !>
  SUBROUTINE InvertHPD(A)

    COMPLEX(8), intent(in out)                      ::  A(:,:)
	 INTEGER                	   	            ::  INFO, LDA,N,i,j

	 LDA = size(A,1)
	 N = size(A,2)
	 info=0

	 CALL ZPOTRF('U',N,A,LDA,INFO)  !-- Performs the Choelesky factorisation

	 if (info==0) then

		CALL ZPOTRI('U',N,A,LDA,INFO) !-- CAREFUL: returns only triangular part
		do i=1,N
		  do j=1,N
			 if (i>j) then
				a(i,j) = a(j,i)
			 end if
		  end do
		end do

		if (info /= 0) then
		  print*, "Failure in the inversion step, ZPOTRI"
		end if
	 else
	 	print*, "Failure in ZPOTRF, choelesky factorisation"
	 end if

  END SUBROUTINE InvertHPD

!> -------------------------------------------------------------------------
  !> SUBROUTINE: InvertGeneral
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Inverts a general dense complex matrix in-place. 
  !>   Computes the LU factorization using `ZGETRF` (with partial pivoting) 
  !>   and subsequently generates the full matrix inverse using `ZGETRI`.
  !> Arguments:
  !>   - A    : The general complex matrix to be inverted (modified in-place).
  !>   - info : Output flag from LAPACK (0 indicates successful inversion).
  !>
  SUBROUTINE InvertGeneral(A,info)
    COMPLEX(cx), intent(in out)                  ::  A(:,:)
    integer, intent(out)								 ::  info
	 INTEGER                	   		          ::  LDA,N,LWORK
	 INTEGER, dimension(size(A,1))		          ::  IPIV
	 COMPLEX(cx), allocatable 				          ::  WORK(:)

	 LDA = size(A,1)
	 N = size(A,1)
	 info=0
	 LWORK = N
	 allocate(WORK(LWORK))

	 CALL ZGETRF(N,N,A,LDA,IPIV,INFO)

	 if (info==0) then

		CALL ZGETRI(N,A,LDA,IPIV,WORK,-1,INFO)

		if (info /= 0) then
		  print*, "Failure in the inversion step, ZGETRI"
		  print*,"info=", info
		end if

	 else
		print*, "Failure in ZGETRF"
		print*,"info=", info
	 end if

  END SUBROUTINE InvertGeneral

!> -------------------------------------------------------------------------
  !> SUBROUTINE: InvertH
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Inverts a generic complex Hermitian matrix in-place. 
  !>   Uses the Bunch-Kaufman diagonal pivoting method via `ZHETRF`, followed 
  !>   by `ZHETRI`. Because LAPACK only overwrites the upper triangle, this 
  !>   routine explicitly reconstructs the lower half by taking the complex 
  !>   conjugate of the mirrored upper elements.
  !> Arguments:
  !>   - A    : The complex Hermitian matrix (modified in-place to its inverse).
  !>   - info : Output flag from LAPACK (0 indicates successful inversion).
  !>
  SUBROUTINE InvertH(A,info)
    COMPLEX(cx), intent(in out)                  ::  A(:,:)
    integer, intent(out)								 ::  info
	 INTEGER                	   		          ::  LDA,N,LWORK,i,j
	 INTEGER, dimension(size(A,1))		          ::  IPIV
	 COMPLEX(cx), allocatable 				          ::  WORK(:)

	 LDA = size(A,1)
	 N = size(A,2)
	 info=0
	 LWORK = N
	 allocate(WORK(LWORK))

	 CALL ZHETRF('U',N,A,LDA,IPIV,WORK,LWORK,INFO)  !-- Performs the Bunch-Kaufman factorisation

	 if (info==0) then

		CALL ZHETRI('U',N,A,LDA,IPIV,WORK,INFO) !-- CAREFUL: returns only triangular part
		do i=1,N
		  do j=1,N
			 if (i>j) then
				a(i,j) = conjg(a(j,i))
			 end if
		  end do
		end do

		if (info /= 0) then
		  print*, "Failure in the inversion step, ZHETRI"
		  print*,"info=", info
		end if

	 else
		print*, "Failure in ZHETRF, Bunch-Kaufman factorisation"
		print*,"info=", info
	 end if

  END SUBROUTINE InvertH

!> -------------------------------------------------------------------------
  !> SUBROUTINE: SolveEq_r
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Solves a system of real linear equations of the form A * X = B. 
  !>   Wraps the standard LAPACK double-precision general solver `DGESV`, 
  !>   which performs an LU decomposition with partial pivoting.
  !> Arguments:
  !>   - A : The real square coefficient matrix (overwritten by LU factors).
  !>   - B : On entry, the right-hand side vector; on exit, the solution vector X.
  !>
  SUBROUTINE SolveEq_r(A,B)
    REAL(rl), intent(in out)               ::  A(:,:)
	 REAL(rl), intent(in out)					 ::  B(:)
	 INTEGER                	   		    ::  INFO,LDA,LDB,N,NRHS
	 INTEGER, dimension(size(A,1))		    ::  IPIV   !-- pivot indices

	 NRHS = 1  						!-- number of right hand sides
	 N = size(A,1)					!-- the number of linear equations
	 LDA = size(A,1)				!-- the leading dimension of A multiply_c
	 LDB = size(B,1)				!-- the leading dimension of B
	 info=0							!-- 0 is successful

	 CALL DGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)  !-- Solve by performing the Bunch-Kaufman factorisation

	 if (info /= 0) then
	 	print*, "Failure in DGESV - solving real system of equations"
	 	print*,"INFO = ",info
	 end if

  END SUBROUTINE SolveEq_r

!> -------------------------------------------------------------------------
  !> SUBROUTINE: SolveEq_c
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Solves a system of complex linear equations of the form A * X = B. 
  !>   Wraps the LAPACK complex general solver `ZGESV`, performing an LU 
  !>   decomposition to evaluate the solution vector.
  !> Arguments:
  !>   - A : The complex square coefficient matrix (overwritten by LU factors).
  !>   - B : On entry, the right-hand side vector; on exit, the solution vector X.
  !>
  SUBROUTINE SolveEq_c(A,B)

    COMPLEX(cx), intent(in out)               ::  A(:,:)
	 COMPLEX(cx), intent(in out)					 ::  B(:)
	 INTEGER                	   		    	 ::  INFO,LDA,LDB,N,NRHS
	 INTEGER, dimension(size(A,1))		    	 ::  IPIV   !-- pivot indices

	 NRHS = 1  						!-- number of right hand sides
	 N = size(A,1)					!-- the number of linear equations
	 LDA = size(A,1)				!-- the leading dimension of A multiply_c
	 LDB = size(B,1)				!-- the leading dimension of B
	 info=0							!-- 0 is successful

	 CALL ZGESV(N,NRHS,A,LDA,IPIV,B,LDB,INFO)  !-- Solve by performing the Bunch-Kaufman factorisation

	 if (info /= 0) then
	 	print*, "Failure in DGESV - solving real system of equations"
	 	print*,"INFO = ",info
	 end if

 END SUBROUTINE

!> -------------------------------------------------------------------------
  !> FUNCTION: matmultiply_c
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Performs high-performance complex matrix multiplication (C = A * B). 
  !>   Directly wraps the level-3 BLAS routine `ZGEMM`, bypassing slow, unoptimized 
  !>   nested Fortran loops. This is the underlying procedure mapped to the 
  !>   custom `.matprod.` operator.
  !> Arguments:
  !>   - amat : The left-hand complex matrix.
  !>   - bmat : The right-hand complex matrix.
  !> Return:
  !>   - outmat : The resulting complex matrix product.
  !>
  FUNCTION matmultiply_c(amat,bmat) RESULT(outmat)

    complex(cx), dimension(:,:), intent(in)            :: amat
    complex(cx), dimension(:,:), intent(in)            :: bmat
    complex(cx), dimension(size(amat,1),size(bmat,2))  :: outmat

    call zgemm ('N','N',size(amat,1),size(bmat,2),size(amat,2),one,amat,&
         size(amat,1),bmat,size(bmat,1), zero,outmat,size(outmat,1))

  END FUNCTION matmultiply_c

END MODULE lapackmodule
