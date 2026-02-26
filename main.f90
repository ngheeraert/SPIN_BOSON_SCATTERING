! ==============================================================================
!  FILE / PROGRAM: main.f90 (main)
!
!  PURPOSE & CONTEXT
!    Main entry point for the numerical simulation of ultrastrong-coupling 
!    waveguide QED:
!      N. Gheeraert et al., "Particle production in ultrastrong-coupling 
!      waveguide QED", Phys. Rev. A 98, 043816 (2018).
!
!  PHYSICAL MODEL & METHODOLOGY
!    Implements a time-dependent variational simulation of the spin-boson model 
!    (a two-level system coupled to a one-dimensional continuum of bosonic modes). 
!    The system's wave-function is approximated using a superposition of multimode 
!    coherent states, commonly referred to as the multi-polaron or MCS ansatz.
!
!  CORE RESPONSIBILITIES & WORKFLOW
!    Acts as the top-level driver for the simulation suite. Its primary 
!    responsibilities include:
!      1. Initialization : Parses command-line arguments to populate system parameters.
!      2. Execution      : Triggers the main time-evolution trajectory (`printTrajectory_DL`).
!      3. Benchmarking   : Logs core CPU execution time.
!      4. Post-Analysis  : Contains conditional blocks (based on `prep` flags) 
!                          to reload states from file and evaluate specific 
!                          observables (like g2 correlation functions) without 
!                          rerunning the full time evolution.
!
!  EXECUTION HIERARCHY
!      main.f90 
!       └── output.f90:printTrajectory_DL 
!            └── output.f90:evolveState_DL 
!                 └── output.f90:Evolve_RK4 (integrator)
!                      └── systm.f90:CalcDerivatives (computes variational EOMs)
!
!  BUILD / DEPENDENCIES
!    * BLAS / LAPACK (Requires standard routines: ZGESV, ZGETRF, ZTRSM, ZPOTRF).
!
! ==============================================================================
!
PROGRAM main

  USE systm
  USE output

  IMPLICIT NONE

  TYPE(param) 	 					::  sys,sys2
  type(state)						::  st, ini_st
  real(rl)							::  t1,t2

! ----------------------------------------------------------------------------
  ! 1. INITIALIZATION
  ! ----------------------------------------------------------------------------
  ! Parse the command-line arguments to build the `sys` (param) object.
  ! This sets up the physical constants (alpha, delta), grid sizes, 
  ! and wavepacket characteristics before allocating the state arrays.
  print*, "-- Starting parameter initialisation - MAIN"
  CALL getParameters(sys)
  print*, "-- Parameters initialised - MAIN"

! ----------------------------------------------------------------------------
  ! 2. PRIMARY TIME EVOLUTION TRAJECTORY
  ! ----------------------------------------------------------------------------
  ! The main execution block. It passes the uninitialized states (`st`, `ini_st`) 
  ! to the output driver, which handles memory allocation, state preparation, 
  ! RK4 integration, and observable generation.

  !== DOUBLE-LINE TRAJECTORY
  !IF ((sys%prep .ge. 50) .and. (sys%prep < 100)) then
		
print*, "-- Simulating a double-line - MAIN"
! Start benchmarking timer
CALL CPU_TIME(t1)

! Execute the core simulation workflow
CALL printTrajectory_DL(sys,st,ini_st)

! End benchmarking timer and output duration
CALL CPU_TIME(t2)


print*,"==================================================="
write(*,'(a15,f10.2)') "COMPUTING TIME", t1-t2
print*," "

! ----------------------------------------------------------------------------
  ! 3. CONDITIONAL POST-PROCESSING / REPLOTTING
  ! ----------------------------------------------------------------------------
  ! If a specific `prep` flag block is triggered (2601 - 2699), the program skips 
  ! the time evolution. Instead, it reloads a previously finished state from disk 
  ! and calculates expensive spatial multi-photon correlations (g2 functions) 
  ! that might have been skipped during the initial run.

IF ( (sys%prep > 2600) .and. ( sys%prep < 2700 ) ) then

! Clone the system parameters and adjust the basis size to account 
     ! for polarons added dynamically during the original run
	 sys2=sys
	 sys2%npini = sys%npini + sys%npadd
	 sys%prep = sys%prep - 2600

! Allocate the state arrays for the target dimension
	 CALL allocate_state(sys2,st)
	 CALL print_param(sys)

! Reconstruct the multi-polaron state from the Even/Odd basis files
	 CALL initialise_from_file_eo(sys,st)
	 CALL print_g2_0(sys,st,"  fst")

! Compute zero-delay and spatially-scanned 2nd-order coherences (g2)
	 CALL print_g2(sys,-681._rl,0._rl,st)

 END IF

END PROGRAM main








