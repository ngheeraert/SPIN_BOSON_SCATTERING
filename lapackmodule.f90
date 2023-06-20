MODULE lapackmodule

  USE typedefs, only : cx => c_type, rl=> r_type
  USE consts

  implicit none

  INTERFACE OPERATOR(.matprod.)
    module procedure matmultiply_c
  END INTERFACE OPERATOR(.matprod.)

CONTAINS

  ! --> Inverse of a Hermitian Positive Definite matrix
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

  ! --> Solve General Complex Equations
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

  ! --> Matrix muliplication
  FUNCTION matmultiply_c(amat,bmat) RESULT(outmat)

    complex(cx), dimension(:,:), intent(in)            :: amat
    complex(cx), dimension(:,:), intent(in)            :: bmat
    complex(cx), dimension(size(amat,1),size(bmat,2))  :: outmat

    call zgemm ('N','N',size(amat,1),size(bmat,2),size(amat,2),one,amat,&
         size(amat,1),bmat,size(bmat,1), zero,outmat,size(outmat,1))

  END FUNCTION matmultiply_c

END MODULE lapackmodule
