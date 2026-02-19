! ==============================================================================
!  FILE: main.f90
!
!  PURPOSE
!    Commented / documented version of the original research code used to produce
!    the numerical data underlying the figures of:
!      N. Gheeraert et al., Phys. Rev. A 98, 043816 (2018) "Particle production in
!      ultrastrong-coupling waveguide QED".
!
!  OVERVIEW
!    This code implements a time-dependent variational simulation of the spin-boson
!    model (a two-level system coupled to a continuum of bosonic modes) using a
!    superposition of multimode coherent states (sometimes called the multi-polaron
!    or MCS ansatz). The main workflow is:
!      main.f90 -> output:printTrajectory_DL -> output:evolveState_DL -> RK4 time-step
!      with systm:CalcDerivatives computing the variational equations of motion.
!
!  BUILD / DEPENDENCIES
!    * Free-form Fortran 90/95 code.
!    * Requires BLAS/LAPACK (ZGESV, ZGETRF, ZTRSM, ZPOTRF, ...).
!    * Preprocessor symbol DP may be used to select double precision in typedefs.f90.
!
!  NOTE ON COMMENTING
!    Only comments have been added. Computational statements and control flow are
!    intentionally left unchanged to preserve bitwise-identical behavior as much as
!    possible (compiler-dependent).
!
!  Generated on: 2026-02-19
! ==============================================================================
!
PROGRAM main
!> -------------------------------------------------------------------------
!> PROGRAM: main
!> -------------------------------------------------------------------------
!> Purpose / context:
!>   Program `main`: main entry point.
!>   Parses parameters and runs a full simulation trajectory.
!> Arguments:
!>   (none)
!>

  USE systm
  USE output

  IMPLICIT NONE

  TYPE(param) 	 					::  sys,sys2
  type(state)						::  st, ini_st
  real(rl)							::  t1,t2

  !-- initialise the system parameters
  print*, "-- Starting parameter initialisation - MAIN"
  CALL getParameters(sys)
  print*, "-- Parameters initialised - MAIN"

  !== DOUBLE-LINE TRAJECTORY
  !IF ((sys%prep .ge. 50) .and. (sys%prep < 100)) then
		
	 print*, "-- Simulating a double-line - MAIN"
	 CALL CPU_TIME(t1)
	 CALL printTrajectory_DL(sys,st,ini_st)
	 CALL CPU_TIME(t2)


	 print*,"==================================================="
	 write(*,'(a15,f10.2)') "COMPUTING TIME", t1-t2
	 print*," "

	 !== REPLOTTING CERTAIN QUANTITIES

  IF ( (sys%prep > 2600) .and. ( sys%prep < 2700 ) ) then

	 sys2=sys
	 sys2%npini = sys%npini + sys%npadd
	 sys%prep = sys%prep - 2600
	 CALL allocate_state(sys2,st)
	 CALL print_param(sys)
	 CALL initialise_from_file_eo(sys,st)
	 CALL print_g2_0(sys,st,"  fst")
	 CALL print_g2(sys,-681._rl,0._rl,st)

  END IF

END PROGRAM main








