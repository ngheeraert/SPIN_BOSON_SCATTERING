! ==============================================================================
!  FILE: output.f90
!
!  PURPOSE & CONTEXT
!    Core physics engine and data structures for the numerical simulation of 
!    ultrastrong-coupling waveguide QED: 
!      N. Gheeraert et al., "Particle production in ultrastrong-coupling 
!      waveguide QED", Phys. Rev. A 98, 043816 (2018).
!
!  PHYSICAL MODEL & METHODOLOGY
!    Implements a time-dependent variational simulation of the spin-boson model 
!    (a two-level system coupled to a one-dimensional continuum of bosonic modes). 
!    The system's wave-function is approximated using a superposition of multimode 
!    coherent states, commonly referred to as the multi-polaron or MCS (Multiple 
!    Coherent States) ansatz.
!
!  CORE RESPONSIBILITIES & WORKFLOW
!    This module orchestrates the time evolution and data export. It drives 
!    the adaptive 4th-order Runge-Kutta (RK4) integration scheme, with dynamic 
!    basis enlargement.
!
!      1. Data Structures : Manages the fundamental types (`param`, `state`, `traj`).
!      2. Time Evolution  : Drives the adaptive RK4 loop and handles basis 
!                           dimension changes (adding/removing polarons).
!      3. Observables     : Computes macroscopic physical quantities (reflection, 
!                           transmission, photon densities, and g2 correlations).
!      4. Data Export     : Formats and writes diagnostic/physical data to disk.
!
!  EXECUTION HIERARCHY
!      main.f90 
!       └── output.f90:printTrajectory_DL 
!            └── output.f90:evolveState_DL 
!                 └── output.f90:Evolve_RK4 (integrator)
!                      └── systm.f90:CalcDerivatives (computes variational EOMs)
!
!  BUILD / DEPENDENCIES
!    * BLAS / LAPACK (Requires standard routines: ZGESV, ZGETRF, ZTRSM, ZPOTRF)
!
! ==============================================================================
!
MODULE output

  USE systm
  USE consts

  implicit none

CONTAINS

!> -------------------------------------------------------------------------
  !> SUBROUTINE: printTrajectory_DL
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Top-level driver for the time-dependent variational simulation.
  !>   Initializes the system state based on the preparation flag, handles the 
  !>   time-evolution loop, triggers basis expansion (polaron adding), and 
  !>   computes final physical observables (reflection, transmission, correlations).
  !> Arguments:
  !>   - sys  : Parameter structure defining the physical system and grid.
  !>   - st   : The dynamic variational state (current time).
  !>   - st_0 : The initial state at t=0, preserved for reference and observable calculations.
  !>
  SUBROUTINE printTrajectory_DL(sys,st,st_0) 

    type(param), intent(in out)			::  sys
	 type(state), intent(in out)			::  st,st_0
	 type(traj) 								::  tr
	 type(state) 								::  wp_st,wp_st_0, cl_st_0, st_ba
	 real(rl)									::  R,T
	 integer										::  adding
	 character(len=200)						::  k0char, name_check_sol
	 
	 !========================================================
	 !-- Allocation and initialisation
	 !=======================================================

	 IF_MAKEGIF: IF ( sys%makegif == 1 ) THEN
		CALL open_gif_files(sys)
	 ENDIF IF_MAKEGIF

	 CALL allocate_state(sys,st_0)
	 CALL allocate_state(sys,st)
	 CALL allocate_trajectory(sys,tr)
	 print*, "-- ARRAYS ALLOCATED"

	!-- FULL-CHAIN PREPARATIONS
	 !-- Select the initial physical configuration of the wavepacket and qubit
	 if ( sys%prep == -1 ) then
		 !-- Custom or pre-allocated state (do nothing)
	 else if ( sys%prep == 50 ) then
	 	!-- Gaussian coherent state wavepacket + qubit in ground state
		CALL initialise_cs_gs(sys,st_0,sys%k0,sys%x0,sys%sigma,sys%n_wp,cl_st_0)
		CALL remove_gs(sys,st_0,cl_st_0,wp_st_0)
		st = st_0
	 else if ( sys%prep == 51 ) then
	 	!-- Single photon Fock state + qubit in ground state
		CALL initialise_photon_gs(sys,st_0,sys%k0,sys%x0,sys%sigma,sys%n_wp,cl_st_0)
		CALL remove_gs(sys,st_0,cl_st_0,wp_st_0)
		st = st_0
	 else if ( sys%prep == 52 ) then
	 	!-- Single photon wavepacket (no explicit ground state separation)
		CALL initialise_photon(sys,st_0,sys%k0,sys%x0,sys%sigma,sys%n_wp,cl_st_0)
		st = st_0
	 else if ( sys%prep == 54 ) then
	 	!-- Bi-chromatic coherent state wavepacket + qubit in ground state
		CALL initialise_cs_2freqs_gs(sys,st_0,sys%k0,sys%k0_2,sys%x0,sys%sigma,sys%n_wp,cl_st_0)
		CALL remove_gs(sys,st_0,cl_st_0,wp_st_0)
		st = st_0
	 else if ( sys%prep == 60 ) then
	    !-- Prepare as prep=50 but subsequently overwrite with a state from file (Even/Odd basis)
		sys%prep=sys%prep-10
		CALL initialise_cs_gs(sys,st_0,sys%k0,sys%x0,sys%sigma,sys%n_wp,cl_st_0)
		CALL remove_gs(sys,st_0,cl_st_0,wp_st_0)
		CALL initialise_from_file_eo(sys,st)
	 else
	 	print*, "Error in the value of prep"
	 	stop
	 end if

	 print*,"-- state initialisation complete"

	 !========================================================
	 !-- Evolving the state until sys%tmax
	 !=======================================================

	 !-- Evolution including the adding of polarons
	 ADDING_IF: IF (sys%npadd .ne. 0) THEN
	 	adding = 1
	 ELSE
	 	adding=0
	 END IF ADDING_IF
	 st_ba=st

	 TIME_EVOLUTION: DO

		CALL evolveState_DL(sys,st,tr,sys%tmax,adding,st_ba) !-- The 1 specifies that routine should stop for adding pol

		if ( st%t >= sys%tmax ) exit time_evolution

		if ( adding==1 ) then

		  st_ba = st
		  CALL add_coherent_state(sys,st)
		  tr%np_ar(tr%i_time) = st%np
		  if ( sys%npini + sys%npadd <= st%np ) then
		  	 adding=0
		  end if

		end if

	 END DO TIME_EVOLUTION

	 CALL end_of_GIF(sys)

	 !========================================================
	 !-- Plotting the final state of the system
	 !========================================================

	 if ( st%t < sys%tmax-5._rl ) then
		sys%tmax = int(st%t)
	 end if
	 CALL print_evolution_data(sys,tr)
	 CALL print_param(sys)

	 !-- print out f(x) and f(k) at t=tmax
	 !CALL print_fxs(sys,st,"  fst")
	 CALL print_nk_eo(sys,wp_st_0,"  iwp")
	 CALL print_nx_eo(sys,st_0,"  ist")
	 CALL print_fks(sys,st,"  fst")
	 CALL print_ps(sys,st,"  fst")
	 CALL print_nx_eo(sys,st,"  fst")

	 wp_st = filtered_wp_eo(sys,st)

	 !CALL print_1ph_losses(sys,wp_st,wp_st_0)
	 CALL print_nk_eo(sys,wp_st,"  fwp")
	 CALL print_nphk_eo(sys,wp_st_0,"  fwp",3)
	 CALL print_nphk_eo(sys,wp_st,"  fwp",3)


	 R = reflection(sys,wp_st,wp_st_0)
	 T = transmission(sys,wp_st,wp_st_0)

	 write(*,'(a30,f12.9)') "reflection ", R
	 write(*,'(a30,f12.9)') "transmission ", T

	 write(k0char,'(f8.4)') sys%k0
	 open (unit=105,file="file_"//trim(adjustl(k0char))//".d",action="write",status="replace")
	 write(105,'(5f25.15)',advance='no') sys%k0, R, T
	 write(105,*)
	 close(105)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: evolveState_DL
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Inner time evolution loop executing the variational equations of motion.
  !>   It adaptively monitors the local energy/norm error and adjusts the Runge-Kutta 
  !>   time-step width (slowfactor). Can interrupt the loop to return control to 
  !>   the main driver for dynamically adding a coherent component to the basis.
  !> Arguments:
  !>   - sys   : Parameter structure.
  !>   - st    : Current state to be evolved.
  !>   - tr    : Trajectory object recording observables over time.
  !>   - tmax  : Maximum simulation time for this specific evolution block.
  !>   - add   : Flag indicating if the system should attempt to add polarons.
  !>   - st_ba : Backup state used for rollbacks when integration error spikes.
  !>
  SUBROUTINE evolveState_DL(sys,st,tr,tmax, add, st_ba) 

    type(param), intent(in)				::  sys
	 type(state), intent(in out)	  		::  st, st_ba
	 type(traj), intent(in out)	  		::  tr
    real(rl), intent(in)					::  tmax
	 integer,intent(in out)					::  add
	 type(state)								::  ost,oost,st_t3
	 integer  									::  i, counter,nb_clean_evolve
	 real(rl)									::  tini,dtprint,tprinted,slowfactor
	 real(rl)									::  err_instant, err_reference, last_E_check
	 logical										::  err_ref_defined, superinverseflag

	 !-- initialise the states corresponding to previous time-steps
	 ost = st
	 oost = st
		
	 tini = st%t

	 !-- error monitoring variables
	 err_instant = 0._rl
	 err_reference = 0._rl
	 err_ref_defined = .false. 
	 nb_clean_evolve = 0

	 !-- variables for printing interval
	 tprinted = 0._rl
	 counter=0

	 last_E_check = st%t
	 st_t3=st

	 slowfactor = 1._rl
	 dtprint = 20*sys%dt/dble(slowfactor)
	 superinverseflag = .true.

	 EVOLUTION: DO

		if ( tmax <= st%t ) exit evolution

		!if ( ( add == 1 ) &
		!		.and. ( ((st%t >= tini + sys%tref) .and. (st%np > sys%npini)) &
		!		.or.  ( ( err_ref_defined ) .and. (abs(err_instant-err_reference) >= sys%merr) .and. (st%np==sys%npini) ) ) ) then 

		!  	exit evolution

		!end if
		!print*, abs(err_instant-err_reference), err_instant,err_reference
		if ( st%t >= tini + sys%tref  &
				.and.  ( add == 1 ) &
				.and.  ( err_ref_defined ) &
				.and.  ( ( (abs(err_instant-err_reference) >= sys%merr/2).and.(st%np==sys%npini) ) &
							 .or. (abs(err_instant-err_reference) >= sys%merr) ) ) then

		  !print*, abs(err_instant-err_reference), err_instant, err_reference
		  	exit evolution

		end if


		IF ( st%t < tini+sys%tref ) THEN
		  dtprint = 1000._rl
		  superinverseflag = .false.
		  if (st%t < 0.0000001_rl+tini)  then
			 slowfactor = 1.0e7_rl
		  elseif ( (st%t < 0.000001_rl+tini) .and. (st%t > 0.0000001_rl+tini) ) then
			 slowfactor = 1.0e6_rl
		  elseif ( (st%t < 0.00001_rl+tini) .and. (st%t > 0.000001_rl+tini) ) then
			 slowfactor = 1.0e5_rl
		  elseif ( (st%t < 0.0001_rl+tini) .and. (st%t > 0.00001_rl+tini) ) then
			 slowfactor = 1.0e4_rl
		  elseif ( (st%t < 0.001_rl+tini) .and. (st%t > 0.0001_rl+tini) ) then
			 slowfactor = 1.0e3_rl
		  elseif ( (st%t < 0.01_rl+tini) .and. (st%t > 0.001_rl+tini) ) then
			 slowfactor = 1.0e2_rl
		  elseif ( (st%t < 0.05_rl+tini) .and. (st%t > 0.01_rl+tini) ) then
			 slowfactor = 50
		  elseif ( (st%t < 0.5_rl+tini) .and. (st%t > 0.05_rl+tini) ) then
			 slowfactor = 5
		  elseif ( (st%t < tini+sys%tref) .and. (st%t > 0.5_rl+tini) ) then
			 slowfactor = 1
		  end if
		ELSE IF ( st%np < sys%npini+sys%npadd ) THEN
		  slowfactor = 1
		  dtprint = 20*sys%dt
		  superinverseflag = .true.
		ELSE
		  slowfactor = 1
		  dtprint = 200*sys%dt/dble(slowfactor)
		  superinverseflag = .true.
		END IF


		IF_MAKEGIF: IF ( sys%makegif == 1 ) THEN
		  CALL print_gif_image(sys,st)
		ENDIF IF_MAKEGIF

		!-- PRINTOUT AND CALCULATE THE ERROR
		PRINT_IF: IF ( ( (st%t-tprinted > dtprint) .and. (err_ref_defined) )  &
							 .OR. (st%t<1e-13) &
							 .OR. ((st%t - tini > sys%tref ) .and. (.not. err_ref_defined)) ) THEN

		  !-- calculate 3 points with dt/10000, to obtain smooth curve for error
		  !nb_clean_evolve = 0
		  DO_ERR_PREP: DO i=1,4

			 CALL evolve_rk4( sys,st,ost,oost,1.0e5_rl,superInverseFlag )
			 nb_clean_evolve = nb_clean_evolve+1

		  END DO DO_ERR_PREP

		  err_instant = real( error(sys,oost,ost,st) )

		  tr%time_ar(tr%i_time) = st%t
		  tr%error_ar(tr%i_time) = err_instant
		  tr%energy_ar(tr%i_time) = energy(sys,st)
		  tr%norm_ar(tr%i_time) = norm(st)
		  tr%spinXYZ_ar(tr%i_time,1) = sigmaX(st)
		  tr%spinXYZ_ar(tr%i_time,2) = sigmaY(st)
		  tr%spinXYZ_ar(tr%i_time,3) = sigmaZ(st)
		  !do n=1,st%np
		  !  tr%ps_ar(tr%i_time,n) =  real( conjg(st%p(n))*st%p(n) ) 
		  !end do

		  tr%i_time = tr%i_time + 1
		  tprinted=st%t

		  if ((.not. err_ref_defined) .and. (st%t - tini > sys%tref )) then
			 err_ref_defined = .true.
			 err_reference = err_instant
			 tr%tref_counter=tr%tref_counter+1
			 tr%tref_ar(tr%tref_counter,1) = st%t
			 tr%tref_ar(tr%tref_counter,2) = err_reference
			 write( *,'(a30,f9.2)' ) "ERROR REFERENCE DEFINED AT t= ", st%t 
		  end if

		  counter=counter+1
		  if (counter>10)  then
			 write(*,*) st%t,'Err=',err_instant,superInverseFlag
			 counter=0
		  end if
		
		ENDIF PRINT_IF

		!========================================================
	 	!== Evolve the system by one time-step

		CALL evolve_RK4(sys,st,ost,oost,slowfactor,superInverseFlag)	!-- 1 means randomising of dt

		IF ( (sys%max_deltaE > 0) &
		  .and. (ost%np == st%np ) &
		  .and. (oost%np == st%np ) ) THEN
		  CALL checkTimestep( sys,st,ost,oost,slowfactor,sys%max_deltaE,tini,st_ba )
		  if (st%np < sys%npini + sys%npadd) then
		  	 add=1
		  end if
		END IF

		IF (st%t > last_E_check + 3._rl) THEN
		  CALL checkTimestep_t3( sys,st,ost,oost,st_t3,5.0e-7_rl )
		  st_t3=st
		  last_E_check = st%t
		END IF

	 END DO EVOLUTION

  END SUBROUTINE

	  !> -------------------------------------------------------------------------
	 !> SUBROUTINE: checkTimestep
	 !> -------------------------------------------------------------------------
	 !> Purpose / context:
	 !>   Local error validation routine. Compares the energy of the system across 
	 !>   successive steps. If the energy drift (deltaE) exceeds `errorlimit`, it 
	 !>   rewinds the state to a previous backup (`oost_tofix`) and locally refines 
	 !>   the RK4 integration by increasing `slowfactorFix` (reducing the time step).
	 !> Arguments:
	 !>   - sys        : Parameter structure.
	 !>   - st, ost, oost : Current, old, and older states (t, t-dt, t-2dt).
	 !>   - slowfactor : Base time-step reduction factor.
	 !>   - errorlimit : Maximum allowable energy drift.
	 !>   - tini       : Initial time reference for backtracking.
	 !>   - st_ba      : Absolute safe backup state in case local rewinding fails.
	 !>
	 SUBROUTINE checkTimestep(sys,st,ost,oost,slowfactor,errorlimit,tini,st_ba)

	   type(param), intent(in)			:: sys
	   type(state), intent(in out)   :: st, ost, oost, st_ba
	   real(rl), intent(in) 	      :: slowfactor
	   real(rl), intent(in out) 	   :: tini
	   integer, save       				:: cnt
	   real(rl), intent(in)          :: errorlimit
	   type(state)							:: st_tofix,ost_tofix,oost_tofix,best_st, best_ost
	   integer                       :: fixtry
	   real(rl)								:: best_deltaE, deltaE, ini_deltaE,slowfactorFix,odeltaE, final_time
	   character(200)						:: name_fixes

		st_tofix = st
		ost_tofix = ost
		oost_tofix = oost
	   deltaE = abs(energy(sys,st) - energy(sys,ost))
	   odeltaE = deltaE
	   ini_deltaE = deltaE
	   best_deltaE = deltaE
	   best_st = st
	   best_ost = ost
	   slowfactorFix = 1._rl
	   fixtry = 0

		FIXING: DO

	     if ( (deltaE < errorlimit) .or. (fixtry > 8) ) exit fixing

	     fixtry = fixtry + 1 	!-- keep track of the number of tries
	     st = oost_tofix				 	!-- rewind evolution by 2 steps
		  print*, "FIXING: t,deltaE,slowfactor,fixtry=",st%t,deltaE,slowfactor,fixtry

		  if (fixtry<5) then
			 slowfactorFix = slowfactor*(fixtry*5)  !-- first try to decrease dt
		  else
			 slowfactorFix = slowfactor*(fixtry*25)  !-- first try to decrease dt
		  end if

		  RE_EVOLVE: do

			 if ( st_tofix%t+0.1*(sys%dt/slowfactor) < st%t ) exit RE_EVOLVE
			 CALL evolve_RK4(sys,st,ost,oost,slowfactorFix,.false.)
	     
		  end do RE_EVOLVE

	     !-- recalculate the energy error
	     odeltaE=deltaE
	     deltaE = abs((energy(sys,st) - energy(sys,oost_tofix)))/sqrt(2.0_rl)

	     if (abs(deltaE) < abs(best_deltaE)) then
	   	 best_deltaE = deltaE
	   	 best_st = st
	   	 best_ost = ost
	     end if


		END DO FIXING

		IF ( (fixtry == 0) .and. (abs(deltaE) < errorlimit) ) then

		  !-- count the number of steps betweeen 2 fixings
		  cnt = cnt + 1

		ELSE IF ( (fixtry .ne. 0) .and. ( deltaE < errorlimit) ) THEN

	     write(*,*) "No fix for ",(cnt)," steps"
	     write(*,*) "SUCCESSFUL, fix attempts:"&
	   	 ,(fixtry)," time=",(st%t) ," New_deltaE=",(deltaE)," Ini_deltaE=", (ini_deltaE)
	     write(*,*) "EVOLVING FURTHER BY 4 TIME STEPS WITH slowfactor=", slowfactorfix
	     write(*,*)
	     cnt=0

		  final_time = st%t + 4*sys%dt/slowfactor
		  FURTHER_EVOLVE: do

			 if ( final_time < st%t ) exit FURTHER_EVOLVE
			 CALL evolve_RK4(sys,st,ost,oost,slowfactorFix,.false.)
	     
		  end do FURTHER_EVOLVE

		ELSE IF ( (fixtry > 8) .and. (abs(deltaE) > errorlimit) ) THEN

		  IF ( ( st_ba%np < st%np ) &
					 .and. ( sys%back_track_on==1 ) &
					 .and. ( st_tofix%t - tini < 8._rl ) ) THEN

		  	 st=st_ba
		  	 tini = st%t
			 print*," "
			 print*,"================================="
			 write(*,'(a23,I2,a16)') "BACK TRACKING BACK TO: ",st%np," COHERENT STATES"
			 print*,"================================="

		  ELSE

		    write(*,*) "FAILED FIXING -- ABORT"
		    write(*,*)
			 name_fixes=trim(adjustl(sys%file_path))//"/FAILED_"//trim(adjustl(parameterchar(sys)))//".d"
			 open (unit=21,file= name_fixes,action="write",status="replace")
			 write(21,*) "No fix for ",cnt," steps"
			 write(21,*) "Time=", st%t

			 close(21)

		    !== CLOSING THE PROGRAM
		    !=======================
		  	 stop
		    !=======================

		  END IF

		END IF


	 END SUBROUTINE
!> -------------------------------------------------------------------------
	 !> SUBROUTINE: checkTimestep_t3
	 !> -------------------------------------------------------------------------
	 !> Purpose / context:
	 !>   Macroscopic error validation routine. Tracks the cumulative energy drift 
	 !>   over a longer, predefined time horizon (default delta_t = 3.0). 
	 !>   Operates similarly to `checkTimestep` but mitigates slow, progressive 
	 !>   integration errors that evade step-by-step local detection.
	 !> Arguments:
	 !>   - sys        : Parameter structure.
	 !>   - st, ost, oost : Trajectory tracking states.
	 !>   - st_t3      : Reference state captured 3 time units prior.
	 !>   - errorlimit : Maximum allowable energy deviation over the interval.
	 !>
	 SUBROUTINE checkTimestep_t3(sys,st,ost,oost,st_t3,errorlimit)

	   type(param), intent(in)			:: sys
	   type(state), intent(in out)   :: st, ost, oost,st_t3
	   integer, save       				:: cnt_t3
	   real(rl), intent(in)          :: errorlimit
	   type(state)							:: st_tofix,ost_tofix,oost_tofix
	   integer                       :: fixtry
	   real(rl)								:: deltaE, ini_deltaE,slowfactorFix,odeltaE, final_time
	   character(200)						:: name_fixes

		st_tofix = st
		ost_tofix = ost
		oost_tofix = oost
	   deltaE = abs(energy(sys,st_t3) - energy(sys,st))
	   odeltaE = deltaE
	   ini_deltaE = deltaE
	   slowfactorFix = 1._rl
	   fixtry = 0

		FIXING: DO

	     if ( (deltaE < errorlimit) .or. (fixtry > 8) ) exit fixing

	     fixtry = fixtry + 1 	!-- keep track of the number of tries
	     st = st_t3				 	!-- rewind evolution by 2 steps
		  print*, "FIXING (OVER t=3): t,deltaE,slowfactor,fixtry=",st%t,deltaE,fixtry

		  slowfactorFix = (fixtry*5)  !-- first try to decrease dt

		  RE_EVOLVE: do

			 if ( st_tofix%t < st%t ) exit RE_EVOLVE
			 CALL evolve_RK4(sys,st,ost,oost,slowfactorFix,.false.)
	     
		  end do RE_EVOLVE

	     !-- recalculate the energy error
	     deltaE = abs((energy(sys,st) - energy(sys,st_t3)))

		END DO FIXING

		IF ( (fixtry == 0) .and. (abs(deltaE) < errorlimit) ) then

		  !-- count the number of steps betweeen 2 fixings
		  cnt_t3 = cnt_t3 + 1

		ELSE IF ( (fixtry .ne. 0) .and. ( deltaE < errorlimit) ) THEN

	     write(*,*) "No fix for ",(cnt_t3)," steps"
	     write(*,*) "SUCCESSFUL (OVER t=3), fix attempts:"&
	   	 ,(fixtry)," time=",(st%t) ," New_deltaE=",(deltaE)," Ini_deltaE=", (ini_deltaE)
	     write(*,*)
	     cnt_t3=0

		ELSE IF ( (fixtry > 8) .and. (abs(deltaE) > errorlimit) ) THEN

		  write(*,*) "FAILED FIXING (OVER t=3) -- ABORT"
		  write(*,*)
		  name_fixes=trim(adjustl(sys%file_path))//"/FAILED_"//trim(adjustl(parameterchar(sys)))//".d"
		  open (unit=21,file= name_fixes,action="write",status="replace")
		  write(21,*) "No fix for ",cnt_t3," steps"
		  write(21,*) "Time=", st%t

		  close(21)

		  !== CLOSING THE PROGRAM
		  !=======================
		  stop
		  !=======================

		END IF


	 END SUBROUTINE
  !> -------------------------------------------------------------------------
  !> SUBROUTINE: print_evolution_data
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Write an observable / diagnostic to a text file (in the `data/` folder).
  !>   The filename is generally parameter-tagged to support parameter sweeps.
  !> Arguments:
  !>   - sys
  !>   - tr
  !>
  SUBROUTINE print_evolution_data(sys,tr)

		type(param), intent(in)  		::  sys
		type(traj), intent(in)			::  tr
		character(len=200)		  		::  name_eren,name_spinXYZ,name_np,name_ps,name_trefs
		real(rl)								::  t
		integer								::  i,n

		name_spinXYZ=trim(adjustl(sys%file_path))//"/spinXYZ_"//trim(adjustl(parameterchar(sys)))//".d"
		name_np=trim(adjustl(sys%file_path))//"/np_"//trim(adjustl(parameterchar(sys)))//".d"
		name_ErEN=trim(adjustl(sys%file_path))//"/ErEN_"//trim(adjustl(parameterchar(sys)))//".d"
		name_trefs=trim(adjustl(sys%file_path))//"/trefs_"//trim(adjustl(parameterchar(sys)))//".d"
		name_ps=trim(adjustl(sys%file_path))//"/ps_"//trim(adjustl(parameterchar(sys)))//".d"
		!name_info=trim(adjustl(sys%file_path))//"/info_"//trim(adjustl(parameterchar(sys)))//".d"
		!name_deltaR=trim(adjustl(sys%file_path))//"/deltaR_"//trim(adjustl(parameterchar(sys)))//".d"

		open (unit=10,file=name_ErEN,action="write",status="replace")
		open (unit=11,file= name_spinXYZ,action="write",status="replace")
		open (unit=12,file= name_np,action="write",status="replace")
		!open (unit=13,file= name_ps,action="write",status="replace")
		open (unit=14,file= name_trefs,action="write",status="replace")

		do n=1, size(tr%tref_ar,1)
		  write(14,*) tr%tref_ar(n,1), tr%tref_ar(n,2)
		end do

		do i=1, tr%i_time	-1

		  t = tr%time_ar(i)
		  write(10,*) t, tr%error_ar(i), tr%energy_ar(i), tr%norm_ar(i)
		  write(11,'(f25.15,3f25.10)') t, tr%spinXYZ_ar(i,1), tr%spinXYZ_ar(i,2), tr%spinXYZ_ar(i,3)
		  write(12,*) t, tr%np_ar(i)

		  !write(13,'(f25.15)',advance='no') t
		  !do n=1,sys%npini + sys%npadd
		!	 write(13,'(f25.15)',advance='no') tr%ps_ar(i,n)
		 ! end do
		  !write(13,*)

		end do

		close(10)
		close(11)
		close(12)
		!close(13)
		close(14)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: Evolve_RK4
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Executes a single step of the standard 4th-order Runge-Kutta integration scheme.
  !>   Evaluates the variational derivatives via `calcderivatives` at four intermediate 
  !>   points to calculate the increment for state arrays `f, h, p, q`. Also applies 
  !>   the free-evolution phase terms `-Ic * sys%w * dt` explicitly to the odd mode 
  !>   arrays `fo, ho`.
  !> Arguments:
  !>   - sys              : Parameter structure.
  !>   - st               : Current state (updated in place to t+dt)
  !>   - ost              : Old state (assigned the incoming value of st).
  !>   - oost             : Older state (assigned the incoming value of ost).
  !>   - slowfactor       : Divisor applied to the base dt to increase resolution.
  !>   - superinverseflag : Boolean flag to toggle full/partial inverse calculations in derivatives.
  !>
  SUBROUTINE Evolve_RK4(sys,st,ost,oost,slowfactor,superinverseflag)

	 type(state),intent(in out)   ::  st, ost, oost
	 type(param),intent(in)   		::  sys
	 real(rl)							::  slowfactor
	 type(state)   					::  midSt
	 real(rl)							::  dt
	 integer								::  n
	 logical								:: superinverseflag

	 oost = ost
	 ost = st 
	 midst = st

	 dt = sys%dt/dble(slowfactor)

	 !== FIRST STEP OF RK4 ==========! 
	 CALL calcderivatives(sys,st,superinverseflag)

	 st%f = st%f + st%fdot*dt/6._rl
	 st%h = st%h + st%hdot*dt/6._rl
	 st%p = st%p + st%pdot*dt/6._rl
	 st%q = st%q + st%qdot*dt/6._rl
	 CALL update_sums(sys,st)

	 !== SECOND STEP OF RK4 ==========! 
	 midSt%f = ost%f + 0.5_rl*dt*st%fdot
	 midSt%h = ost%h + 0.5_rl*dt*st%hdot
	 midSt%p = ost%p + 0.5_rl*dt*st%pdot
	 midSt%q = ost%q + 0.5_rl*dt*st%qdot
	 CALL update_sums(sys,midst)

	 !-- for 1 mode only one cannot use the intrinsic fortran functions
	 !-- hence we use here the manual summations in the case of 1 mode
	 CALL calcderivatives(sys,midst,superinverseflag)

	 st%f = st%f + midst%fdot*dt/3._rl
	 st%h = st%h + midst%hdot*dt/3._rl
	 st%p = st%p + midst%pdot*dt/3._rl
	 st%q = st%q + midst%qdot*dt/3._rl
	 CALL update_sums(sys,st)

	 !== THIRD STEP OF RK4 ==========! 
	 midSt%f = ost%f + 0.5_rl*dt*midst%fdot 
	 midSt%h = ost%h + 0.5_rl*dt*midst%hdot 
	 midSt%p = ost%p + 0.5_rl*dt*midst%pdot 
	 midSt%q = ost%q + 0.5_rl*dt*midst%qdot 
	 CALL update_sums(sys,midst)

	 CALL calcderivatives(sys,midst,superinverseflag)

	 st%f = st%f + midst%fdot*dt/3._rl
	 st%h = st%h + midst%hdot*dt/3._rl
	 st%p = st%p + midst%pdot*dt/3._rl
	 st%q = st%q + midst%qdot*dt/3._rl
	 CALL update_sums(sys,st)

	 !== FOURTH STEP OF RK4 ==========! 
	 midSt%f = ost%f + dt*midst%fdot 
	 midSt%h = ost%h + dt*midst%hdot 
	 midSt%p = ost%p + dt*midst%pdot 
	 midSt%q = ost%q + dt*midst%qdot 
	 CALL update_sums(sys,midst)

	 CALL calcderivatives(sys,midst,superinverseflag)

	 st%f = st%f + midst%fdot*dt/6._rl
	 st%h = st%h + midst%hdot*dt/6._rl
	 st%p = st%p + midst%pdot*dt/6._rl
	 st%q = st%q + midst%qdot*dt/6._rl
	 CALL update_sums(sys,st)
	 do n=1, st%np
		st%fo(n,:) = st%fo(n,:) * exp( -Ic * sys%w(:) * dt )
		st%ho(n,:) = st%ho(n,:) * exp( -Ic * sys%w(:) * dt )
	 end do

	 st%t = st%t + dt

  END SUBROUTINE evolve_RK4

!> -------------------------------------------------------------------------
  !> SUBROUTINE: add_coherent_state
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Implements basis enlargement for the multi-polaron ansatz. Appends a 
  !>   new coherent component to the basis to capture scattering complexity. The new 
  !>   polaron is initialized with vacuum modes (`f` and `h` = 0) and predefined 
  !>   atomic displacement amplitudes `p0`, adjusting phase based on current parity.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state (reallocated to size np + 1).
  !>
  SUBROUTINE add_coherent_state(sys,st)

		type(param), intent(in)			:: sys
		type(state), intent(in out)   :: st
		type(state)						   :: st_add
		type(state)						   :: cl_st_0
		type(state)					 		:: ba_st
		integer								:: n
		real(rl)								::	E1,E2,SZ1,SZ2,SX1,SX2

		E1 = energy(sys,st)
		SZ1 = sigmaZ(st)
		SX1 = sigmaX(st)
		ba_st = st 

		!-- Re-allocate st with more polarons
		CALL allocate_state(sys,st,ba_st%np+1)
		CALL allocate_state(sys,st_add,1)
		CALL allocate_state(sys,cl_st_0,1)
		print*, "-- ARRAYS ALLOCATED"

		st%t = ba_st%t
		st_add%f(:,:) = 0._cx
		st_add%h(:,:) = 0._cx

		do n=1, ba_st%np
		  st%f(n,:) = ba_st%f(n,:)
		  st%h(n,:) = ba_st%h(n,:)
		  st%p(n) = ba_st%p(n)
		  st%q(n) = ba_st%q(n)
		end do
		st%f(ba_st%np+1,:) = st_add%f(1,:)
		st%h(ba_st%np+1,:) = st_add%h(1,:)
		st%p(ba_st%np+1) = sys%p0
		st%q(ba_st%np+1) = (-1)**(st%np) * sys%p0

		CALL update_sums(sys,st)
		CALL normalise(st)

		print*,"================================="
		write(*,'(a21,I5,a4,I3,a6,f9.2)') "POLARON ADDED: from ",st%np-st_add%np," to ",st%np, "at t=", st%t 
		E2 = energy(sys,st)
		SZ2 = sigmaZ(st)
		SX2 = sigmaX(st)
		print*, "delta( Sz ) = ",SZ2-SZ1
		print*, "delta( Sx ) = ",SX2-SX1
		print*, "delta( E ) = ",E2-E1
		print*,"================================="
		print*," "

	 END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: remove_coherent_state
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Reduces the computational overhead by pruning the variational basis. 
  !>   It searches for two highly overlapping coherent states using `find_biggest_Mnp`, 
  !>   merges their amplitudes (`p` and `q`), and collapses the state array dimension 
  !>   by removing the redundant components, resulting in state array `np - 1`.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state (reallocated and compressed to size np - 1).
  !>
  SUBROUTINE remove_coherent_state(sys,st)

		type(param), intent(in)			:: sys
		type(state), intent(in out)   :: st
		type(state)					 		:: br_st
		integer								:: n,p_f,n_f,p_h,n_h
		real(rl)							::	E1,E2,SZ1,SZ2,SX1,SX2,min_M_f,min_M_h

		E1 = energy(sys,st)
		SZ1 = sigmaZ(st)
		SX1 = sigmaX(st)
		br_st = st 

		CALL find_biggest_Mnp( st, min_M_f, n_f, p_f, min_M_h, n_h, p_h  )
		CALL allocate_state(sys,st,br_st%np-1)

		st%t = br_st%t

		do n=1, p_f-1
		  st%f(n,:) = br_st%f(n,:)
		  st%p(n) = br_st%p(n)
		end do
		do n=p_f+1, br_st%np
		  st%f(n-1,:) = br_st%f(n,:)
		  st%p(n-1) = br_st%p(n)
		end do
		st%p(n_f) = br_st%p(n_f) + br_st%p(p_f)
		do n=1, p_h-1
		  st%h(n,:) = br_st%h(n,:)
		  st%q(n) = br_st%q(n)
		end do
		do n=p_h+1,br_st%np
		  st%h(n-1,:) = br_st%h(n,:)
		  st%q(n-1) = br_st%q(n)
		end do
		st%q(n_h) = br_st%q(n_h) + br_st%q(p_h)

		CALL update_sums(sys,st)
		CALL normalise(st)

		print*,"================================="
		write(*,'(a23,I5,a4,I3,a6,f9.2)') "COHERENT STATE REMOVED: from ",br_st%np," to ",st%np, "at t=", st%t 
		write(*,'(a27,f25.15,f25.15)')  "	--BIGGEST Mnp_f, Mnp_h= ", min_M_f, min_M_h
		write(*,'(a30,I3,I3,a8,I3,I3)') "	--COHERENT STATE NUMBERS f:", n_f, p_f,"      h:",n_h,p_h 
		E2 = energy(sys,st)
		SZ2 = sigmaZ(st)
		SX2 = sigmaX(st)
		print*, "Delta( Sz ) = ",SZ2-SZ1
		print*, "Delta( Sx ) = ",SX2-SX1
		print*, "Delta( E ) = ",E2-E1
		print*,"================================="
		print*," "

		!p2min=1000._rl
		!q2min=1000._rl
		!n_p2min = 1
		!n_q2min = 1
		!do n=1, st%np
		!  if ( abs(st%p(n))**2 <p2min ) then
		!  	 p2min = abs(st%p(n))**2
		!  	 n_p2min = n
		!  end if
		!  if ( abs(st%q(n))**2 <q2min ) then
		!  	 q2min = abs(st%q(n))**2
		!  	 n_q2min = n
		!  end if
		!end do
		!if ( (p2min>900._rl) .or. (p2min>900._rl) ) then
		!  print*,"Error in removal of coherent state"
		!end if

	 END SUBROUTINE

!> -------------------------------------------------------------------------
	 !> SUBROUTINE: find_biggest_Mnp
	 !> -------------------------------------------------------------------------
	 !> Purpose / context:
	 !>   Utility for basis pruning. Iterates over the off-diagonal elements of the 
	 !>   bosonic overlap matrices (`ov_ff` and `ov_hh`) to identify the pair of polarons 
	 !>   (i, j) that share the highest quantum overlap, signifying redundancy.
	 !> Arguments:
	 !>   - st    : Current state containing the overlap matrices.
	 !>   - Mnp_f : Output scalar holding the maximum absolute overlap for the 'f' modes.
	 !>   - n_f, p_f : Output indices of the maximally overlapping pair in 'f'.
	 !>   - Mnp_h : Output scalar holding the maximum absolute overlap for the 'h' modes.
	 !>   - n_h, p_h : Output indices of the maximally overlapping pair in 'h'.
	 !>
	 SUBROUTINE find_biggest_Mnp(st,Mnp_f,n_f,p_f,Mnp_h,n_h,p_h)

		type(state), intent(in)    ::  st
		real(rl), intent(out)		::  Mnp_f, Mnp_h
		integer, intent(out)			::  n_f,n_h,p_f,p_h
		integer							::  i,j

		Mnp_f=-1
		Mnp_h=-1

		do i=1, st%np
		  do j=i+1, st%np
			 if ( abs(st%ov_ff(i,j)) > Mnp_f ) then
				Mnp_f = abs(st%ov_ff(i,j))
				n_f=min(i,j)
				p_f=max(i,j)
			 end if
			 if ( abs(st%ov_hh(i,j)) > Mnp_h ) then
				Mnp_h = abs(st%ov_hh(i,j))
				n_h=min(i,j)
				p_h=max(i,j)
			 end if
		  end do
		end do

	 END SUBROUTINE

!> -------------------------------------------------------------------------
  !> FUNCTION: filtered_wp_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the scattered outgoing wavepacket by spatially filtering out 
  !>   components inside the scattering region (`x < xmin`). Operates in the 
  !>   Even/Odd parity basis (EO). It transforms the spatial modes back into 
  !>   momentum space (k-space) and reconstructs the symmetric/antisymmetric 
  !>   modal arrays (`f, fo, h, ho`)[cite: 58, 59].
  !> Arguments:
  !>   - sys     : Parameter structure.
  !>   - st      : The full state containing all spatial components.
  !>   - xmin_in : Optional cutoff spatial coordinate (defaults to `sys%xmin`).
  !> Return:
  !>   - Type(state) : A normalized cloned state containing only the filtered wavepacket.
  !>
  FUNCTION filtered_wp_eo(sys,st,xmin_in) 

    type(param), intent(in)			  	::  sys
	 type(state), intent(in)				::  st
	 real(rl), intent(in), optional		::  xmin_in
	 type(state)								::  filtered_wp_eo
	 type(state)								::	 tmp_st
	 complex(cx), dimension( st%np, sys%nmode )	::  fnk, fnmk, hnk, hnmk
	 real(rl)									::  xmin

	 if (present(xmin_in)) then
		xmin = xmin_in	
	 else
	 	xmin = sys%xmin
	 end if

	 tmp_st = st
	 fnk = FT_X_to_K_EO(sys,f_nx_eo(sys,st),-1,xmin)
	 fnmk = FT_X_to_K_EO(sys,f_nx_eo(sys,st),1,xmin)
	 hnk = FT_X_to_K_EO(sys,h_nx_eo(sys,st),-1,xmin)
	 hnmk = FT_X_to_K_EO(sys,h_nx_eo(sys,st),1,xmin)

	 tmp_st%f = sqrt(0.5_rl) * ( fnk + fnmk ) 
	 tmp_st%fo = sqrt(0.5_rl) * ( fnk - fnmk ) 
	 tmp_st%h = sqrt(0.5_rl) * ( hnk + hnmk ) 
	 tmp_st%ho = sqrt(0.5_rl) * ( hnk - hnmk ) 
	 
	 CALL update_sums(sys,tmp_st)
	 CALL normalise(tmp_st)
	 
	 filtered_wp_eo = tmp_st


  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: left_going_st_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts only the left-propagating components of the photonic field.
  !>   It recombines the Even/Odd (EO) basis modal arrays (f, fo, h, ho) to 
  !>   filter out right-propagating modes, then reconstructs the EO arrays 
  !>   for the purely left-going state.
  !> Arguments:
  !>   - sys : Parameter structure defining the grid and modes.
  !>   - st  : The full input state from which to extract the left-propagating field.
  !> Return:
  !>   - Type(state) : A normalized cloned state containing only the left-going wavepacket.
  !>
  FUNCTION left_going_st_eo(sys,st) 

    type(param), intent(in)			  	::  sys
	 type(state), intent(in)				::  st
	 type(state)								::  left_going_st_eo
	 type(state)								::	 tmp_st
	 complex(cx), dimension( st%np, sys%nmode )	:: fnmk, hnmk

	 tmp_st = st
	 fnmk = sqrt(0.5_rl) * ( st%f - st%fo )  !-- left propagating field
	 hnmk = sqrt(0.5_rl) * ( st%h - st%ho )  !-- left propagating field

	 tmp_st%f = sqrt(0.5_rl)* fnmk
	 tmp_st%fo = - sqrt(0.5_rl)*fnmk
	 tmp_st%h = sqrt(0.5_rl)*hnmk
	 tmp_st%ho = -sqrt(0.5_rl)* hnmk
	 
	 CALL update_sums(sys,tmp_st)
	 CALL normalise(tmp_st)
	 
	 left_going_st_eo = tmp_st

  END FUNCTION

!> -------------------------------------------------------------------------
  !> SUBROUTINE: remove_gs_fromfile
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Subtracts the bound ground state (GS) from the current state to isolate 
  !>   the scattered wavepacket. Instead of passing a computed GS, this routine 
  !>   reads the GS modal arrays directly from a pre-generated data file and 
  !>   subtracts them from the state arrays `f` and `h`.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state, modified in-place to remove the GS contribution.
  !>
  SUBROUTINE remove_gs_fromfile(sys,st)

	 type(param),intent(in)								 ::  sys
	 type(state),intent(in out)						 ::  st
	 integer 												 ::  k,i
	 complex(cx), dimension(sys%npini,sys%nmode)  ::  gs_f, gs_h
	 real(rl) 												 ::  a,b
	 logical													 ::  ex

	 gs_f = 0._cx
	 gs_h = 0._cx

	 inquire( file=gs_filename(sys,st), exist = ex )
	 if ( .not. ex ) then
		print*, "error: Ground State file does not exist:"
		print*, trim(gs_filename(sys,st))
		print*
		stop
	 end if
	 open(unit=1000,file=gs_filename(sys,st),status="old",action="read")
	 read(1000,*)
	 read(1000,'(a1)', advance="no") b
	 read(1000,*)
	 do k=1,sys%nmode
		read(1000,'(f11.5)', advance="no") b
		do i=1,st%np
		  read(1000,'(f11.5)', advance="no") a
		  gs_f(i,k) = a
		  gs_h(i,k) = - a
		end do
		read(1000,*)
	 end do
	 close(1000)

	 st%f = st%f - gs_f
	 st%h = st%h - gs_h

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: remove_gs
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Subtracts an explicitly provided ground state (`gst`) from the current 
  !>   variational state (`st`). This is typically done prior to evaluating 
  !>   scattering observables (like reflection/transmission) to ensure the 
  !>   qubit's bound cloud is not mistakenly counted as free propagating photons.
  !> Arguments:
  !>   - sys    : Parameter structure.
  !>   - st     : Current total state.
  !>   - gst    : The bare ground state of the strongly coupled qubit.
  !>   - out_st : Output state containing only the separated, free-propagating wavepacket. 
  !>
  SUBROUTINE remove_gs(sys,st,gst,out_st)

	 type(param),intent(in)								 ::  sys
	 type(state),intent(in)								 ::  st,gst
	 type(state),intent(out)							 ::  out_st

	 out_st = st

	 out_st%f = st%f - gst%f
	 out_st%h = st%h - gst%h

	 CALL update_sums(sys,out_st)
	 CALL normalise(out_st)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> FUNCTION: transmission
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the energy transmission coefficient. It computes the total mean 
  !>   photon number for all positive momentum modes (k > 0) in the final state, 
  !>   and divides it by the total input energy (mean photon number of the initial 
  !>   state).
  !> Arguments:
  !>   - sys    : Parameter structure.
  !>   - st     : The final scattered state (ground state typically removed).
  !>   - ini_st : The initial input state before scattering.
  !> Return:
  !>   - real(rl) : The transmission probability (T).
  !>
  FUNCTION transmission(sys,st,ini_st)

	 type(param),intent(in) 					 		::  sys
	 type(state),intent(in)						 		::  st, ini_st
	 real(rl)										 		::  transmission
	 real(rl)										 		::  t_output,input
	 real(rl),dimension(-sys%nmode+1:sys%nmode) 	::  nk_ini,nk

	 nk_ini = 0._rl
	 nk = 0._rl

	 !nk_ini = n_up_k_eo(sys,upst_ini)
	 !nk = n_up_k_eo(sys,upst)
	 nk_ini = n_k_eo(sys,ini_st)
	 nk = n_k_eo(sys,st)

	 input = sum( nk_ini(1:sys%nmode) )
	 t_output = sum( nk(1:sys%nmode) )

	 transmission = t_output/input

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: reflection
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the energy reflection coefficient. It sums the mean photon 
  !>   number across all negative momentum modes (k < 0) in the final state, 
  !>   and normalizes it against the total mean photon number of the initial state.
  !> Arguments:
  !>   - sys    : Parameter structure.
  !>   - st     : The final scattered state (ground state typically removed).
  !>   - ini_st : The initial input state before scattering.
  !> Return:
  !>   - real(rl) : The reflection probability (R).
  !>
  FUNCTION reflection(sys,st,ini_st)

	 type(param),intent(in) 					 			::  sys
	 type(state),intent(in)						 			::  st, ini_st
	 real(rl)										 			::  reflection
	 real(rl)										 			::  r_output,input
	 real(rl), dimension(-sys%nmode+1:sys%nmode) 	::  nk_ini,nk,w
	 integer							 				 			::  k

	 nk_ini = 0._rl
	 nk = 0._rl

	 nk_ini = n_k_eo(sys,ini_st)
	 nk = n_k_eo(sys,st)

	 do k=-sys%nmode+1, 0
	 	w(k) = sys%w(abs(k-1))
	 end do
	 do k=1, sys%nmode
		w(k) = sys%w(k)
	 end do

	 input = sum( nk_ini(1:sys%nmode) )
	 r_output = sum( nk(-sys%nmode+1:0) )

	 reflection = r_output/input

  END FUNCTION

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_1ph_losses
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes and writes out the single-photon reflection, transmission, and 
  !>   loss spectra. It isolates the one-photon amplitudes (projected on both 
  !>   qubit up/down states) within the momentum spread (`k0 +/- sigma`) of the 
  !>   initial pulse, calculating the inelastic losses strictly in the 1-photon subspace. 
  !> Arguments:
  !>   - sys    : Parameter structure (contains file path definitions).
  !>   - st     : The final scattered state. 
  !>   - ini_st : The initial input wavepacket. 
  !>
  SUBROUTINE print_1ph_losses(sys,st,ini_st)

	 type(param),intent(in) 					 			::  sys
	 type(state),intent(in)						 			::  st, ini_st
	 real(rl)										 			::  k_min, k_max, loss_1ph,R_1ph,T_1ph
	 real(rl), dimension(-sys%nmode+1 : sys%nmode)	::  nk_1ph_ini, nk_1ph_end
	 integer							 				 			::  k_ind_max,k_ind_min,k
	 character(len=300)					   				::  name_RT

	 name_RT=trim(adjustl(sys%file_path))//"/RT_1mRmT_"//trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100,file= name_RT,action="write",status="replace")

	 nk_1ph_ini = ( abs( one_photon_k_amp_up(sys,ini_st) )**2 + abs( one_photon_k_amp_down(sys,ini_st) )**2 ) / sys%dk1
	 nk_1ph_end = ( abs( one_photon_k_amp_up(sys,st) )**2 + abs( one_photon_k_amp_down(sys,st) )**2 ) / sys%dk1

	 k_min = sys%k0 - sys%sigma
	 k_max = sys%k0 + sys%sigma
	 k_ind_max = int(k_max/sys%dk1+0.5)
	 k_ind_min = int(k_min/sys%dk1+0.5)

	 do k= k_ind_min, k_ind_max

		if ( k <1 ) then
		  print*, "ERROR in RT calculation"
		  exit
		end if
		R_1ph = nk_1ph_end(-k+1) / nk_1ph_ini(k)
		T_1ph = nk_1ph_end(k) / nk_1ph_ini(k)
		loss_1ph = 1 - ( nk_1ph_end(k) + nk_1ph_end(-k+1) )  / nk_1ph_ini(k)

		write(100,*) sys%w(k), R_1ph, T_1ph, loss_1ph

	 end do
	 close(100)

	 print*,"-- Reflection, transmission and losses printed to file"

  END SUBROUTINE


!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_fks
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Writes the momentum-space (k-space) wavepacket amplitudes (`fnk`, `hnk`) 
  !>   to a text file. It explicitly transforms the symmetric and antisymmetric 
  !>   Even/Odd modes back into the full negative and positive momentum grid 
  !>   before writing the real and imaginary components. 
  !> Arguments:
  !>   - sys   : Parameter structure.
  !>   - st    : State object containing the photonic modal arrays.
  !>   - label : String tag appended to the output filename (e.g., 'fst', 'iwp').
  !>
  SUBROUTINE print_fks(sys,st,label)

	 type(param), intent(in)											::  sys
	 type(state), intent(in)											::  st
	 character(len=5), intent(in)    								::  label
	 integer																	::  k,i
	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode) 	::  fnk, hnk
	 real(rl), dimension(-sys%nmode+1:sys%nmode)					::  w
	 character(len=200)													::  name_fks

	 name_fks=trim(adjustl(sys%file_path))//"/fks_"//trim(adjustl(label))//"_"//&
									 trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100,file= name_fks,action="write",status="replace")

	 do i = -sys%nmode+1, 0
		fnk(:,i) = sqrt(0.5_rl)*(st%f(:,-i+1) - st%fo(:,-i+1))
		hnk(:,i) = sqrt(0.5_rl)*(st%h(:,-i+1) - st%ho(:,-i+1))
		w(i) = -sys%w(abs(i)+1)
	 end do
	 fnk(:,1:sys%nmode) = sqrt(0.5_rl)*(st%f + st%fo)
	 hnk(:,1:sys%nmode) = sqrt(0.5_rl)*(st%h + st%ho)
	 w(1:sys%nmode) = sys%w

	 do k = - sys%nmode+1, sys%nmode
		write(100,'(f25.15)',advance='no') w(k)
		do i=1,st%np
		  write(100,'(2f25.15)',advance='no') real( fnk(i,k) ), real( hnk(i,k) )
		end do
		do i=1,st%np
		  write(100,'(2f25.15)',advance='no') aimag( fnk(i,k) ), aimag( hnk(i,k) )
		end do
		write(100,*)
	 end do
	 write(100,*)
	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_ps
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Writes the qubit/atomic variational parameters (`p` and `q`) to a text file. 
  !>   These complex amplitudes represent the weight of the qubit's up/down 
  !>   displacements for each individual polaron in the multi-polaron superposition. 
  !> Arguments:
  !>   - sys   : Parameter structure. 
  !>   - st    : State object containing the `p` and `q` scalar arrays.
  !>   - label : String tag appended to the output filename. 
  !>
  SUBROUTINE print_ps(sys, st, label)

	 type(param), intent(in)			::  sys
	 type(state), intent(in)			::  st
	 character(len=5), intent(in)    ::  label
	 integer									::  i
	 character(len=200)					::  name_ps

	 name_ps=trim(adjustl(sys%file_path))//"/ps_"//trim(adjustl(label))//"_"//&
									 trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100,file= name_ps,action="write",status="replace")

	 do i=1,st%np
		write(100,'(2f25.15)',advance='no') real( st%p(i) ), real( st%q(i) )
	 end do
	 do i=1,st%np
		write(100,'(2f25.15)',advance='no') aimag( st%p(i) ), aimag( st%q(i) ) 
	 end do
	 write(100,*)

	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_fxs
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Writes the spatial (x-space) photonic amplitudes (`fnx`, `hnx`) to a text file. 
  !>   It maps the momentum modes onto the defined spatial grid `sys%dx`, evaluating 
  !>   the real and imaginary spatial wavefunctions for all polarons. 
  !> Arguments:
  !>   - sys   : Parameter structure defining the spatial grid `dx`. 
  !>   - st    : State object to be evaluated.
  !>   - label : String tag appended to the output filename. 
  !>
  SUBROUTINE print_fxs(sys,st,label)

	 type(param), intent(in)		 ::  sys
	 type(state), intent(in)		 ::  st
	 character(len=5), intent(in)	 ::  label
	 integer								 ::  n,i,min_index,max_index
	 complex(cx), allocatable  	 ::  fnx(:,:), hnx(:,:)
	 character(len=200)				 ::  name_fxs

	 name_fxs=trim(adjustl(sys%file_path))//"/fxs_"//trim(adjustl(label))//"_"//&
									 trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100,file= name_fxs,action="write",status="replace")
	 
	 allocate(fnx(st%np,-sys%nmode+1:sys%nmode))
	 allocate(hnx(st%np,-sys%nmode+1:sys%nmode))
	 fnx = f_nx_eo(sys,st)
	 hnx = h_nx_eo(sys,st)
	 min_index = - sys%nmode+1
	 max_index = sys%nmode+1

	 do i= min_index, max_index

		write(100,'(f25.15)',advance='no') sys%dx*dble(i-0.5_rl)
		do n=1, st%np
		  write(100,'(2f25.15)',advance='no') real( fnx(n,i) ), real( hnx(n,i) )
		end do
		do n=1, st%np
		  write(100,'(2f25.15)',advance='no') aimag( fnx(n,i) ), aimag( hnx(n,i) )
		end do
		write(100,*)

	 end do
	 write(100,*)
	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_nk_EO
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes and outputs the mean photon number momentum spectrum (n_k). 
  !>   It writes the total photon number, as well as the spin-projected 
  !>   (qubit up and qubit down) photon numbers across both negative and 
  !>   positive frequency axes. 
  !> Arguments:
  !>   - sys   : Parameter structure. 
  !>   - st    : Current state object. 
  !>   - label : String tag appended to the output filename. 
  !>
  SUBROUTINE print_nk_EO(sys,st,label) 

	 type(param), intent(in)							  ::  sys
	 type(state), intent(in)							  ::  st
	 character(len=5), intent(in)	 					  ::  label
	 character(len=200)									  ::  xminchar,filename
	 integer													  ::  i, k
	 real(rl), dimension( -sys%nmode+1:sys%nmode ) ::  n_k, n_up_k, n_down_k , w
	 real(rl)												  ::  col2, col3, col4

	 write(xminchar, '(I5)') int(sys%xmin)
	 filename=trim(adjustl(sys%file_path))//"/nk_"//trim(adjustl(label))//"_"&
	 //trim(adjustl(parameterchar(sys)))//"_"&
	 //trim(adjustl(xminchar))//".d"
	 open (unit=105,file= filename,action="write",status="replace")

	 n_k = n_k_EO(sys,st)
	 n_up_k = n_up_k_EO(sys,st)
	 n_down_k = n_down_k_EO(sys,st)

	 do k=-sys%nmode+1, 0
	 	w(k) = -sys%w(abs(k-1))
	 end do
	 do k=1, sys%nmode
		w(k) = sys%w(k)
	 end do

	 !--careful with array ranges!
	 do i=-sys%nmode+1, sys%nmode
		col2=n_k(i)/sys%dk1
		col3=n_up_k(i)/sys%dk1
		col4=n_down_k(i)/sys%dk1
		write(105,'(f25.10,3f25.15)') w(i), col2,col3,col4
	 end do
	 close(105)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_nk_EO_k0
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes and outputs the momentum-space mean photon number (n_k), similar 
  !>   to `print_nk_EO`, but selectively zeroes out or flags modes outside a defined 
  !>   bandwidth (k0 +/- 2*sigma) around the central wavepacket momentum. 
  !>   Useful for isolating the signal within the primary excitation bandwidth.
  !> Arguments:
  !>   - sys   : Parameter structure defining the momentum grid (w, dk1) and k0.
  !>   - st    : State object containing the multi-polaron amplitudes.
  !>   - label : String tag appended to the output filename.
  !>
  SUBROUTINE print_nk_EO_k0(sys,st,label) 

	 type(param), intent(in)							  ::  sys
	 type(state), intent(in)							  ::  st
	 character(len=5), intent(in)	 					  ::  label
	 character(len=200)									  :: xminchar,filename
	 integer													  ::  i, k, k_ind_max_n,k_ind_max_p,k_ind_min_n,k_ind_min_p
	 real(rl), dimension( -sys%nmode+1:sys%nmode ) ::  n_k, n_up_k, n_down_k , w
	 real(rl)												  ::  col2, col3, col4, k_min, k_max

	 write(xminchar, '(I5)') int(sys%xmin)
	 filename=trim(adjustl(sys%file_path))//"/nk_k0_"//trim(adjustl(label))//"_"&
	 //trim(adjustl(parameterchar(sys)))//"_"&
	 //trim(adjustl(xminchar))//".d"
	 open (unit=105,file= filename,action="write",status="replace")

	 n_k = n_k_EO(sys,st)
	 n_up_k = n_up_k_EO(sys,st)
	 n_down_k = n_down_k_EO(sys,st)

	 do k=-sys%nmode+1, 0
	 	w(k) = -sys%w(abs(k-1))
	 end do
	 do k=1, sys%nmode
		w(k) = sys%w(k)
	 end do

	 k_min = sys%k0 - 2*sys%sigma
	 k_max = sys%k0 + 2*sys%sigma

	 k_ind_max_p = int(k_max/sys%dk1+0.5)
	 k_ind_min_p = int(k_min/sys%dk1+0.5)
	 k_ind_max_n = int(-k_min/sys%dk1+0.5)
	 k_ind_min_n = int(-k_max/sys%dk1+0.5)

	 !--careful with array ranges!
	 do i=-sys%nmode+1, k_ind_min_n-1
		col2=n_k(i)/sys%dk1
		col3=n_up_k(i)/sys%dk1
		col4=n_down_k(i)/sys%dk1
		write(105,'(f25.10,3f25.15,I4)') w(i), 0.0,0.0,0.0,0
	 end do
	 do i=k_ind_min_n, k_ind_max_n
		col2=n_k(i)/sys%dk1
		col3=n_up_k(i)/sys%dk1
		col4=n_down_k(i)/sys%dk1
		write(105,'(f25.10,3f25.15,I4)') w(i), col2,col3,col4,1
	 end do
	 do i=k_ind_max_n+1, k_ind_min_p-1
		col2=n_k(i)/sys%dk1
		col3=n_up_k(i)/sys%dk1
		col4=n_down_k(i)/sys%dk1
		write(105,'(f25.10,3f25.15,I4)') w(i), 0.0,0.0,0.0,0
	 end do
	 do i=k_ind_min_p, k_ind_max_p
		col2=n_k(i)/sys%dk1
		col3=n_up_k(i)/sys%dk1
		col4=n_down_k(i)/sys%dk1
		write(105,'(f25.10,3f25.15,I4)') w(i), col2,col3,col4,1
	 end do
	 do i=k_ind_max_p+1, sys%nmode
		col2=n_k(i)/sys%dk1
		col3=n_up_k(i)/sys%dk1
		col4=n_down_k(i)/sys%dk1
		write(105,'(f25.10,3f25.15,I4)') w(i), 0.0,0.0,0.0,0
	 end do
	 close(105)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_nx_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates and writes the spatial photon density distribution, n(x). 
  !>   It projects the symmetric (Even) and antisymmetric (Odd) momentum modes 
  !>   back into real space to compute the local photon number expectation value 
  !>   at each spatial grid point.
  !> Arguments:
  !>   - sys          : Parameter structure defining the spatial grid (dx).
  !>   - st           : State object evaluated at the current time.
  !>   - label        : String tag for the standard data filename or frame index.
  !>   - makefilename : Optional flag to toggle output to the dedicated directory.
  !>
  SUBROUTINE print_nx_eo(sys,st,label,makefilename) 

    type(param), intent(in)		 ::  sys
	 type(state), intent(in)		 ::  st
	 character(len=5), intent(in)	 ::  label
	 type(state)						 ::  upst
	 character(len=200)				 ::  filename
	 integer, intent(in), optional ::  makefilename
	 integer								 ::  n,m,i
	 complex(cx)						 ::  fnx( st%np,-sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  tmp

	 if ( .not. present(makefilename) ) then
		filename=trim(adjustl(sys%file_path))//"/nx_"//trim(adjustl(label))//"_"&
		  //trim(adjustl(parameterchar(sys)))//".d"
	 else
	 	filename = trim(adjustl(gif_dirname(sys)))//"/"//trim(adjustl(label))//".d"
	 end if
	 open (unit=105,file= filename,action="write",status="replace")

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)

	 fnx = f_nx_eo(sys,st)

	 tmp = 0._cx
	 do i=-sys%nmode+1,sys%nmode
		do n=1,st%np
		  do m=1,st%np

			 tmp = tmp + conjg(upst%p(n))*upst%p(m)*(conjg(fnx(n,i))*fnx(m,i))*ov( fnx(n,:) , fnx(m,:) )

		  end do
		end do
		write(105,'(f25.15,f25.15)') sys%dx*(dble(i)-0.5_rl), real(tmp)/sys%dx
		tmp=0._cx
	 end do
	 close(105)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_g2
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the spatially resolved second-order optical coherence function, 
  !>   g2(x1, x2). This measures two-photon spatial correlations (photon bunching 
  !>   or antibunching). Fixes coordinate x1 and scans x2 across a specified range, 
  !>   evaluating expectation values of creation/annihilation operator pairs.
  !> Arguments:
  !>   - sys : Parameter structure defining grid spacing.
  !>   - x1  : The fixed reference spatial coordinate.
  !>   - x2  : A reference coordinate used to define the scanning boundary (imin to imax).
  !>   - st  : Current state object containing the photonic fields.
  !>
  SUBROUTINE print_g2(sys,x1,x2,st)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 real(rl)									::  g2,g2_f, g2_h
	 real(rl),intent(in)						::  x1,x2
	 integer								   	::  i1,i2,i,m,n,imin,imax
	 complex(cx)								::  tmp_f,tmp_f2,tmp_f3,tmp_h,tmp_h2,tmp_h3
	 complex(cx)								::  fnx(sys%npini+sys%npadd,-sys%nmode+1:sys%nmode)
	 complex(cx)								::  hnx(sys%npini+sys%npadd,-sys%nmode+1:sys%nmode)
	 character(len=300)					   ::  name_g2, x1char

	 write(x1char,'(I10)') int(x1)
	 print*, "x1",x1char
	 name_g2=trim(adjustl(sys%file_path))//"/g2_new_"//trim(adjustl(x1char))//"_"//&
				  	 trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=200, file= name_g2, action="write",status="replace")

	 i1 =  int( x1 / sys%dx )!+ 0.5_rl)
	 i2 =  int( x2 / sys%dx )!+ 0.5_rl)
	 print*,"dx=", sys%dx
	 !imin = min(i1,i2)
	 !imax = max(i1,i2)
	 imin = i1 !- abs(i1-i2)
	 imax = i1 + abs(i1-i2)

	 fnx = f_nx_eo(sys,st)
	 hnx = h_nx_eo(sys,st)

	 do i=imin, imax

		tmp_f = 0._cx
		tmp_f2= 0._cx
		tmp_f3= 0._cx
		tmp_h = 0._cx
		tmp_h2= 0._cx
		tmp_h3= 0._cx
		do n=1,st%np
		  do m=1,st%np
			 tmp_f = tmp_f + conjg(st%p(m))*st%p(n)*conjg(fnx(m,i1))*conjg(fnx(m,i))*fnx(n,i)*fnx(n,i1)*st%ov_ff(m,n)
			 tmp_f2 = tmp_f2 + conjg(st%p(m))*st%p(n)*conjg(fnx(m,i1))*fnx(n,i1)*st%ov_ff(m,n)
			 tmp_f3 = tmp_f3 + conjg(st%p(m))*st%p(n)*conjg(fnx(m,i))*fnx(n,i)*st%ov_ff(m,n)

			 tmp_h = tmp_h + conjg(st%q(m))*st%q(n)*conjg(hnx(m,i1))*conjg(hnx(m,i))*hnx(n,i)*hnx(n,i1)*ov( hnx(m,:) , hnx(n,:))!*st%ov_hh(m,n)
			 tmp_h2 = tmp_h2 + conjg(st%q(m))*st%q(n)*conjg(hnx(m,i1))*hnx(n,i1)*ov( hnx(m,:) , hnx(n,:))!st%ov_hh(m,n)
			 tmp_h3 = tmp_h3 + conjg(st%q(m))*st%q(n)*conjg(hnx(m,i))*hnx(n,i)*ov( hnx(m,:) , hnx(n,:))!st%ov_hh(m,n)
		  end do
		end do

		g2 =  real( tmp_f+tmp_h ) / (real(tmp_f2 + tmp_h2)*real(tmp_f3 + tmp_h3)  )
		g2_f = real(tmp_f) / ( real(tmp_f2 + tmp_h2)*real(tmp_f3 + tmp_h3)  )
		g2_h = real(tmp_h) / ( real(tmp_f2 + tmp_h2)*real(tmp_f3 + tmp_h3) )

		write(200,'(4f25.15)') dble(i-i1)*sys%dx, g2, 2*g2_f, 2*g2_h
		print*,dble(i-i1)*sys%dx, g2, 2*g2_f, 2*g2_h

	 end do

	 close(200)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_g2_0
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the equal-position second-order coherence function, g2(x, x), 
  !>   across the entire spatial grid. It quantifies the probability of finding 
  !>   two photons at the exact same location, decomposed into contributions 
  !>   from the symmetric (f) and antisymmetric (h) fields.
  !> Arguments:
  !>   - sys   : Parameter structure.
  !>   - st    : State object to evaluate.
  !>   - label : String tag appended to the output filename.
  !>
  SUBROUTINE print_g2_0(sys,st, label)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 real(rl)									::  g2,g2_f, g2_h
	 character(len=5), intent(in)		   ::  label
	 integer								   	::  i,m,n
	 complex(cx)								::  tmp_f,tmp_f2,tmp_h,tmp_h2
	 complex(cx)								::  fnx(sys%npini+sys%npadd,-sys%nmode+1:sys%nmode)
	 complex(cx)								::  hnx(sys%npini+sys%npadd,-sys%nmode+1:sys%nmode)
	 character(len=300)					   ::  name_g2

	 name_g2=trim(adjustl(sys%file_path))//"/g2_0_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100, file= name_g2, action="write",status="replace")

	 fnx = f_nx_eo(sys,st)
	 hnx = h_nx_eo(sys,st)

	 do i= - sys%nmode+1, sys%nmode

		tmp_f = 0._cx
		tmp_f2= 0._cx
		tmp_h = 0._cx
		tmp_h2= 0._cx
		do n=1,st%np
		  do m=1,st%np
			 tmp_f = tmp_f + conjg(st%p(m))*st%p(n)*conjg(fnx(m,i))*conjg(fnx(m,i))*fnx(n,i)*fnx(n,i)*st%ov_ff(m,n)
			 tmp_f2 = tmp_f2 + conjg(st%p(m))*st%p(n)*conjg(fnx(m,i))*fnx(n,i)*st%ov_ff(m,n)

			 tmp_h = tmp_h + conjg(st%q(m))*st%q(n)*conjg(hnx(m,i))*conjg(hnx(m,i))*hnx(n,i)*hnx(n,i)*st%ov_hh(m,n)
			 tmp_h2 = tmp_h2 + conjg(st%q(m))*st%q(n)*conjg(hnx(m,i))*hnx(n,i)*st%ov_hh(m,n)
		  end do
		end do

		g2_f = real(tmp_f)/real(tmp_f2)**2
		g2_h = real(tmp_h)/real(tmp_h2)**2
		g2 = real(tmp_h+tmp_f)/real(tmp_h2+tmp_f2)**2

		write(100,'(4f25.15)') dble(i-0.5)*sys%dx, g2, g2_f, g2_h

	 end do

	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_g2_0fh
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   An alternative evaluation of the equal-position coherence g2(x, x) that 
  !>   computes the raw combined field intensity correlations without decomposing 
  !>   the normalization strictly by symmetric/antisymmetric terms in the output.
  !>   Can initialize an Even/Odd state from file if `st_in` is not provided.
  !> Arguments:
  !>   - sys   : Parameter structure.
  !>   - st_in : Optional input state object; if absent, state is loaded from file.
  !>   - label : String tag appended to the output filename.
  !>
  SUBROUTINE print_g2_0fh(sys,st_in, label)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in),optional	::  st_in
	 type(state)								::  st
	 character(len=5), intent(in)		   ::  label
	 integer								   	::  i,m,n
	 complex(cx)								::  tmp,tmp2
	 complex(cx)								::  fnx(sys%npini+sys%npadd,-sys%nmode+1:sys%nmode)
	 complex(cx)								::  hnx(sys%npini+sys%npadd,-sys%nmode+1:sys%nmode)
	 character(len=300)					   ::  name_g2

	 if (present(st_in)) then
		st = st_in	
	 else
		CALL allocate_state(sys,st,sys%npini+sys%npadd)
		CALL initialise_from_file_eo(sys,st)
	 end if

	 name_g2=trim(adjustl(sys%file_path))//"/g2_0fh_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100, file= name_g2, action="write",status="replace")

	 fnx = f_nx_eo(sys,st)
	 hnx = h_nx_eo(sys,st)

	 do i= - sys%nmode+1, sys%nmode

		tmp = 0._cx
		tmp2= 0._cx
		do n=1,st%np
		  do m=1,st%np
			 tmp = tmp + conjg(st%p(m))*st%p(n)*conjg(fnx(m,i))*conjg(fnx(m,i))*fnx(n,i)*fnx(n,i)*st%ov_ff(m,n) &
					 + conjg(st%q(m))*st%q(n)*conjg(hnx(m,i))*conjg(hnx(m,i))*hnx(n,i)*hnx(n,i)*st%ov_hh(m,n)
			 tmp2 = tmp2 + conjg(st%p(m))*st%p(n)*conjg(fnx(m,i))*fnx(n,i)*st%ov_ff(m,n) &
					 + conjg(st%q(m))*st%q(n)*conjg(hnx(m,i))*hnx(n,i)*st%ov_hh(m,n)
		  end do
		end do

		write(100,'(2f25.15)') dble(i)*sys%dx, real(tmp)/real(tmp2)**2

	 end do

	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_2ph_prob_kk_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the two-photon momentum probability distribution, P(k1, k2). 
  !>   Outputs the normalized squared amplitudes of the multi-polaron state 
  !>   projected onto the two-photon Fock basis, resolving inelastic scattering 
  !>   where a single high-energy photon splits into two lower-energy ones.
  !> Arguments:
  !>   - sys   : Parameter structure detailing the momentum grid ranges.
  !>   - st    : Current state object.
  !>   - label : String tag appended to the output filename.
  !>
  SUBROUTINE print_2ph_prob_kk_eo(sys,st,label)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 character(len=8), intent(in)		   ::  label
	 character(len=300)					   ::  name_numbst
	 integer							 			::  k1,k2
	 real(rl)									::  w( -sys%nmode+1 : sys%nmode )
	 complex(cx)								::  two_ph_amp_up( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )
	 complex(cx)								::  two_ph_amp_down( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )


	 name_numbst=trim(adjustl(sys%file_path))//"/2ph_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100,file= name_numbst,action="write",status="replace")
	 two_ph_amp_up = two_photon_kk_amp_up(sys,st)
	 two_ph_amp_down = two_photon_kk_amp_down(sys,st)


	 do k1=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)

		do k2=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)

		  write(100,*) w(k1), w(k2),&
		  	 ( abs( two_ph_amp_up(k1,k2) )**2 + abs( two_ph_amp_down(k1,k2) )**2 ) / (sys%dk1)**2

		end do
		write(100,* )

	 end do

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_2ph_prob_kk_eo_PY
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Identical physics context to `print_2ph_prob_kk_eo` (two-photon momentum 
  !>   spectra), but the output text file is explicitly formatted as a 2D matrix 
  !>   with axis headers. This matrix format is directly usable by Python's 
  !>   matplotlib or numpy (e.g., for `imshow` or `contourf` plotting).
  !> Arguments:
  !>   - sys   : Parameter structure.
  !>   - st    : Current state object.
  !>   - label : String tag appended to the output filename (suffixed with '_PY').
  !>
  SUBROUTINE print_2ph_prob_kk_eo_PY(sys,st,label)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 character(len=8), intent(in)		   ::  label
	 character(len=300)					   ::  name_numbst
	 integer							 			::  k1,k2,k
	 real(rl)									::  w( -sys%nmode+1 : sys%nmode )
	 complex(cx)								::  two_ph_amp_up( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )
	 complex(cx)								::  two_ph_amp_down( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )

	 name_numbst=trim(adjustl(sys%file_path))//"/2ph_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//"_PY.d"
	 open (unit=100,file= name_numbst,action="write",status="replace")
	 two_ph_amp_up = two_photon_kk_amp_up(sys,st)
	 two_ph_amp_down = two_photon_kk_amp_down(sys,st)

	 do k=-sys%nmode+1, 0
	 	w(k) = -sys%w(abs(k-1))
	 end do
	 do k=1, sys%nmode
		w(k) = sys%w(k)
	 end do

	 write(100,'(f25.15)',advance='no') 0.0_rl
	 do k1=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)
		  write(100,'(f25.15)',advance='no') w(k1)
	 end do
	 write(100,* )

	 do k2=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)

		write(100,'(f25.15)',advance='no') w(k2)
		do k1=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)
		  write(100,'(f25.15)', advance='no') (abs( two_ph_amp_up(k1,k2) )**2 +abs( two_ph_amp_down(k1,k2) )**2 / (sys%dk1))**2
		end do
		write(100,* )

	 end do
	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_2ph_prob_xx_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the two-photon spatial probability distribution, P(x1, x2). 
  !>   Transforms the multi-polaron representation into the spatial two-photon 
  !>   Fock basis to analyze relative coordinate clustering (e.g., bound photon 
  !>   molecules) in real space.
  !> Arguments:
  !>   - sys   : Parameter structure defining the spatial grid boundaries and dx.
  !>   - st    : Current state object.
  !>   - label : String tag appended to the output filename.
  !>
  SUBROUTINE print_2ph_prob_xx_eo(sys,st,label)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 character(len=5), intent(in)		   ::  label
	 type(state)								::  upst
	 character(len=300)					   ::  name_numbst
	 integer							 			::  i1,i2
	 complex(cx)								::  two_ph_amp( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)

	 name_numbst=trim(adjustl(sys%file_path))//"/2phxx_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100,file= name_numbst,action="write",status="replace")

	 two_ph_amp = two_photon_xx_amp_up(sys,upst)

	 do i1=-sys%nmode+1,sys%nmode

		do i2=-sys%nmode+1, sys%nmode

		  write(100,'(3f25.15)')  i1*sys%dx - 0.5_rl, i2*sys%dx - 0.5_rl, abs( two_ph_amp(i1,i2) )**2 / (sys%dx)**2

		end do
		write(100,* )

	 end do
	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_2ph_prob_xx_eo_PY
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Identical physics context to `print_2ph_prob_xx_eo` (two-photon spatial 
  !>   distribution), but reorganizes the output data into a grid-matrix format 
  !>   with coordinate headers for native ingestion into Python plotting routines.
  !> Arguments:
  !>   - sys   : Parameter structure.
  !>   - st    : Current state object.
  !>   - label : String tag appended to the output filename (suffixed with '_PY').
  !>
  SUBROUTINE print_2ph_prob_xx_eo_PY(sys,st,label)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 character(len=5), intent(in)		   ::  label
	 type(state)								::  upst
	 character(len=300)					   ::  name_numbst
	 integer							 			::  i1,i2
	 complex(cx)								::  two_ph_amp( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)

	 name_numbst=trim(adjustl(sys%file_path))//"/2phxx_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//"_PY.d"
	 open (unit=100,file= name_numbst,action="write",status="replace")

	 two_ph_amp = two_photon_xx_amp_up(sys,upst)

	 write(100,'(f25.15)',advance='no') 0.0_rl
	 do i1=-sys%nmode+1,sys%nmode
		  write(100,'(f25.15)',advance='no') i1*sys%dx - 0.5_rl
	 end do
	 write(100,* )

	 do i2=-sys%nmode+1, sys%nmode

		write(100,'(f25.15)',advance='no') i2*sys%dx - 0.5_rl
		do i1=-sys%nmode+1,sys%nmode
		  write(100,'(f25.15)', advance='no') abs( two_ph_amp(i1,i2) )**2 / (sys%dx)**2
		end do
		write(100,* )

	 end do
	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_3ph_prob_kk_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the three-photon momentum probability distribution. Because 
  !>   visualizing a 3D parameter space (k1, k2, k3) is complex, this routine 
  !>   fixes one of the photon momenta (k3, derived from `w_in`) and plots a 2D 
  !>   slice mapping the remaining coordinates k1 and k2. Includes combinatoric 
  !>   multipliers to account for particle permutations.
  !> Arguments:
  !>   - sys   : Parameter structure.
  !>   - st    : Current state object.
  !>   - label : String tag appended to the output filename.
  !>   - w_in  : The specified frequency (momentum) of the fixed third photon (k3).
  !>
  SUBROUTINE print_3ph_prob_kk_eo(sys,st,label,w_in)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 real(rl), intent(in)					::  w_in
	 character(len=5), intent(in)		   ::  label
	 type(state)								::  upst
	 character(len=300)					   ::  name_numbst,k3_char
	 integer							 			::  k1,k2,k3,k
	 real(rl)									::  w( -sys%nmode+1 : sys%nmode )
	 complex(cx)								::  three_ph_amp( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )
	 integer										::  factor

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)
	 factor=-1

	 k3 = int ( w_in / sys%dk1 )
	 write(k3_char,'(f8.4)') sys%w(k3)

	 name_numbst=trim(adjustl(sys%file_path))//"/3ph_"//trim(adjustl(k3_char))//"_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100,file= name_numbst,action="write",status="replace")

	 three_ph_amp = three_photon_kk_amp_up(sys,upst,k3)

	 do k=-sys%nmode+1, 0
	 	w(k) = -sys%w(abs(k-1))
	 end do
	 do k=1, sys%nmode
		w(k) = sys%w(k)
	 end do

	 do k1=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)

		do k2=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)

		  !-- adding a factor to account for the permutations for k3
		  if ( ( (k3==k2) .and. (k2 .ne. k1) ) &
		    .or. ( (k3==k1) .and. (k1 .ne. k2) ) ) then
		    factor = 2
		  else if  ( (k1 .ne. k2) .and. (k1 .ne. k3) ) then
		  	 factor = 3
		  else
		  	 print*, "ERROR IN print_3ph_prob_kk_eo"
		  	 stop
		  end if

		  write(100,*) w(k1), w(k2), factor * abs( three_ph_amp(k1,k2) )**2 / (sys%dk1)**2

		end do
		write(100,* )

	 end do
	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_3ph_prob_kk_eo_PY
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Identical physics context to `print_3ph_prob_kk_eo` (three-photon momentum 
  !>   spectra evaluated by fixing the third photon's momentum k3). This version 
  !>   formats the output text file as a 2D matrix with explicit axis headers 
  !>   to allow seamless direct ingestion by Python (e.g., matplotlib/numpy).
  !> Arguments:
  !>   - sys   : Parameter structure.
  !>   - st    : Current state object.
  !>   - label : String tag appended to the output filename (suffixed with '_PY').
  !>   - w_in  : The specified momentum/frequency of the fixed third photon (k3).
  !>
  SUBROUTINE print_3ph_prob_kk_eo_PY(sys,st,label,w_in)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 real(rl), intent(in)					::  w_in
	 character(len=5), intent(in)		   ::  label
	 type(state)								::  upst
	 character(len=300)					   ::  name_numbst,k3_char
	 integer							 			::  k1,k2,k3,k
	 real(rl)									::  w( -sys%nmode+1 : sys%nmode )
	 complex(cx)								::  three_ph_amp( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )
	 integer										::  factor

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)
	 factor=-1

	 k3 = int ( w_in / sys%dk1 )
	 write(k3_char,'(f8.4)') sys%w(k3)

	 name_numbst=trim(adjustl(sys%file_path))//"/3ph_"//trim(adjustl(k3_char))//"_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//"_PY.d"
	 open (unit=100,file= name_numbst,action="write",status="replace")

	 three_ph_amp = three_photon_kk_amp_up(sys,upst,k3)

	 do k=-sys%nmode+1, 0
	 	w(k) = -sys%w(abs(k-1))
	 end do
	 do k=1, sys%nmode
		w(k) = sys%w(k)
	 end do

	 write(100,'(f25.15)',advance='no') 0.0_rl
	 do k1=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)
		  write(100,'(f25.15)',advance='no') w(k1)
	 end do
	 write(100,* )

	 do k2=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)

		write(100,'(f25.15)',advance='no') w(k2)
		do k1=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)

		  !-- adding a factor to account for the permutations for k3
		  if ( ( (k3==k2) .and. (k2 .ne. k1) ) &
		    .or. ( (k3==k1) .and. (k1 .ne. k2) ) ) then
		    factor = 2
		  else if  ( (k1 .ne. k2) .and. (k1 .ne. k3) ) then
		  	 factor = 3
		  end if

		  write(100,'(f25.15)', advance='no') factor * abs( three_ph_amp(k1,k2) )**2 / (sys%dk1)**3

		end do
		write(100,* )

	 end do
	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_3ph_prob_xx_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the three-photon spatial probability distribution, P(x1, x2, x3). 
  !>   Fixes one spatial coordinate (x3) and computes the 2D spatial correlation 
  !>   matrix for the remaining two photons across the real-space grid. Includes 
  !>   combinatorial factors to account for bosonic coordinate permutations.
  !> Arguments:
  !>   - sys   : Parameter structure defining the spatial grid `dx`.
  !>   - st    : Current state object.
  !>   - label : String tag appended to the output filename.
  !>   - x_in  : The fixed spatial coordinate for the third photon (x3).
  !>
  SUBROUTINE print_3ph_prob_xx_eo(sys,st,label,x_in)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 real(rl), intent(in)					::  x_in
	 character(len=5), intent(in)		   ::  label
	 type(state)								::  upst
	 character(len=300)					   ::  name_numbst,i3_char
	 integer							 			::  i1,i2,i3
	 complex(cx)								::  three_ph_amp( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )
	 integer										::  factor

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)
	 factor=-1

	 print*,"1"
	 i3 = int ( x_in / sys%dx )
	 write(i3_char,'(I6)') int(x_in)

	 name_numbst=trim(adjustl(sys%file_path))//"/3phxx_"//trim(adjustl(i3_char))//"_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100,file= name_numbst,action="write",status="replace")

	 print*,"2"
	 three_ph_amp = three_photon_xx_amp_up(sys,upst,i3)


	 do i1=-sys%nmode+1,sys%nmode

		do i2=-sys%nmode+1, sys%nmode

		  !-- adding a factor to account for the permutations for i3
		  if ( ( (i3==i2) .and. (i2 .ne. i1) ) &
			 .or. ( (i3==i1) .and. (i1 .ne. i2) ) ) then
			 factor = 2
		  else if  ( (i1 .ne. i2) .and. (i1 .ne. i3) ) then
			 factor = 3
		  else
			 stop
			 print*, "ERROR in print_fks"
		  end if

		  write(100,'(3f25.15)') (i1-0.5)*sys%dx,&
							 (i2-0.5)*sys%dx, &
							 factor * abs( three_ph_amp(i1,i2) )**2 / (sys%dx)**3

		end do
		write(100,*)

	 end do
	 close(100)

	 print*,"3"




  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_3ph_prob_xx_eo_PY
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Identical physics context to `print_3ph_prob_xx_eo` (three-photon spatial 
  !>   distribution mapped over x1 and x2 for a fixed x3). It explicitly reformats 
  !>   the output into a 2D grid-matrix with spatial headers for Python integration.
  !> Arguments:
  !>   - sys   : Parameter structure.
  !>   - st    : Current state object.
  !>   - label : String tag appended to the output filename (suffixed with '_PY').
  !>   - x_in  : The fixed spatial coordinate for the third photon (x3).
  !>
  SUBROUTINE print_3ph_prob_xx_eo_PY(sys,st,label,x_in)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 real(rl), intent(in)					::  x_in
	 character(len=5), intent(in)		   ::  label
	 type(state)								::  upst
	 character(len=300)					   ::  name_numbst,i3_char
	 integer							 			::  i1,i2,i3
	 complex(cx)								::  three_ph_amp( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )
	 integer										::  factor

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)
	 factor=-1

	 print*,"1"
	 i3 = int ( x_in / sys%dx )
	 write(i3_char,'(I6)') int(x_in)

	 name_numbst=trim(adjustl(sys%file_path))//"/3phxx_"//trim(adjustl(i3_char))//"_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//"_PY.d"
	 open (unit=100,file= name_numbst,action="write",status="replace")

	 print*,"2"
	 three_ph_amp = three_photon_xx_amp_up(sys,upst,i3)


	 write(100,'(f25.15)',advance='no') 0.0_rl
	 do i1=-sys%nmode+1,sys%nmode
		  write(100,'(f25.15)',advance='no') (i1-0.5)*sys%dx 
	 end do
	 write(100,* )

	 print*,"3"
	 do i2=-sys%nmode+1, sys%nmode

		write(100,'(f25.15)',advance='no') (i2-0.5)*sys%dx
		do i1=-sys%nmode+1,sys%nmode

		  !-- adding a factor to account for the permutations for i3
		  if ( ( (i3==i2) .and. (i2 .ne. i1) ) &
			 .or. ( (i3==i1) .and. (i1 .ne. i2) ) ) then
			 factor = 2
		  else if  ( (i1 .ne. i2) .and. (i1 .ne. i3) ) then
			 factor = 3
		  else
			 stop
			 print*, "ERROR in print_3ph"
		  end if

		  write(100,'(f25.15)', advance='no') factor * abs( three_ph_amp(i1,i2) )**2 / (sys%dx)**3

		end do
		write(100,* )

	 end do
	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_4ph_prob_kk_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the four-photon momentum probability distribution. Maps a 2D 
  !>   slice of the 4D parameter space by fixing the momenta of two outgoing 
  !>   photons (k3 and k4) and calculating the normalized squared amplitudes 
  !>   for the remaining two variables (k1 and k2).
  !> Arguments:
  !>   - sys    : Parameter structure.
  !>   - st     : Current state object.
  !>   - label  : String tag appended to the output filename.
  !>   - w_in_1 : Frequency/momentum for the fixed third photon (k3).
  !>   - w_in_2 : Frequency/momentum for the fixed fourth photon (k4).
  !>
  SUBROUTINE print_4ph_prob_kk_eo(sys,st,label,w_in_1,w_in_2)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 real(rl), intent(in)					::  w_in_1, w_in_2
	 character(len=5), intent(in)		   ::  label
	 type(state)								::  upst
	 character(len=300)					   ::  name_numbst,k3_char, k4_char
	 integer							 			::  k1,k2,k3,k,k4
	 real(rl)									::  w( -sys%nmode+1 : sys%nmode )
	 complex(cx)								::  four_ph_amp( -sys%nmode+1 : sys%nmode , -sys%nmode+1 : sys%nmode )
	 integer										::  factor

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)
	 factor=-1

	 k3 = int ( w_in_1 / sys%dk1 )
	 k4 = int ( w_in_2 / sys%dk1 )
	 write(k3_char,'(f8.4)') sys%w(k3)
	 write(k4_char,'(f8.4)') sys%w(k4)

	 name_numbst=trim(adjustl(sys%file_path))//"/4ph"//&
		"_"//trim(adjustl(k3_char))//&
		"_"//trim(adjustl(k4_char))//&
		"_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100,file= name_numbst,action="write",status="replace")

	 four_ph_amp = four_photon_kk_amp_up(sys,upst,k3,k4)

	 do k=-sys%nmode+1, 0
	 	w(k) = -sys%w(abs(k-1))
	 end do
	 do k=1, sys%nmode
		w(k) = sys%w(k)
	 end do

	 do k1=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)

		do k2=-int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)

		  write(100,*) w(k1), w(k2), abs( four_ph_amp(k1,k2) )**2 / (sys%dk1)**2

		end do
		write(100,* )

	 end do
	 close(100)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_nphk_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the full momentum-space photonic decomposition of the scattered 
  !>   wavepacket. It recursively projects the state onto N-photon Fock sub-spaces 
  !>   (up to an optional `max_number`) and writes the individual 1-photon, 2-photon, 
  !>   etc. spectral contributions, alongside their sum to verify energy conservation.
  !> Arguments:
  !>   - sys        : Parameter structure.
  !>   - st         : Current state object.
  !>   - label      : String tag appended to the output filename.
  !>   - max_number : Optional upper limit for the N-photon subspace evaluation.
  !>
  SUBROUTINE print_nphk_eo(sys,st,label,max_number)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 character(len=5), intent(in)		   ::  label
	 integer, intent(in),optional			::  max_number
	 character(len=300)					   ::  name_numbst,xminchar
	 integer							 			::  k1,k
	 real(rl)									::  w( -sys%nmode+1 : sys%nmode )
	 real(rl), dimension(-sys%nmode+1 : sys%nmode)	::  tmp1, tmp2, tmp3, tmp4, tmp5
	 real(rl)									::  sum_12345

	 write(xminchar, '(I5)') int(sys%xmin)
	 name_numbst=trim(adjustl(sys%file_path))//"/nphk_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//"_"&
		//trim(adjustl(xminchar))//".d"
	 open (unit=100,file= name_numbst,action="write",status="replace")


	 tmp2(:) = 0._rl
	 tmp3(:) = 0._rl
	 tmp4(:) = 0._rl
	 tmp5(:) = 0._rl

	 print*,"-- Start of spectrum spectrum photonic decomposition"

	 tmp1 = ( abs( one_photon_k_amp_up(sys,st) )**2 + abs( one_photon_k_amp_down(sys,st) )**2 ) / sys%dk1
	 if ( present( max_number ) .and. ( max_number > 1) ) then
		tmp2 = nk_2_photon(sys,st)
	 end if
	 if ( present( max_number ) .and. ( max_number > 2) ) then
		tmp3 = nk_3_photon(sys,st)
	 end if

	 if (present( max_number ) .and. ( max_number > 3) .and. (sys%k0<0.21)) then

		!print*,"      4-photon states START"
	 	!call system('date' )
	 	!tmp4 = nk_4_photon_up(sys,upst)
	 	!print*,"      4-photon states DONE"
	 	!call system('date' )

	 	!print*,"-- Spectrum printed up to 4 photons"

	 end if

	 do k=-sys%nmode+1, 0
	 	w(k) = -sys%w(abs(k-1))
	 end do
	 do k=1, sys%nmode
		w(k) = sys%w(k)
	 end do

	 do k1= -int(sys%nmode/sys%wmax)+1, int(sys%nmode/sys%wmax)

	 	sum_12345 = tmp1(k1) + tmp2(k1) + tmp3(k1) + tmp4(k1) + tmp5(k1) 
		write(100,*) w(k1), tmp1(k1), tmp2(k1), tmp3(k1) , tmp4(k1), tmp5(k1), sum_12345
		if (k1==0) then
		  write(100,*) 0.0_rl, 0.0_rl,0.0_rl, 0.0_rl,0.0_rl,0.0_rl,0.0_rl
		end if

	 end do
	 close(100)

	 print*,"-- Spectrum photonic decomposition printed to file"
	 print*," "

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_nphx_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the real-space spatial photonic density decomposition. It maps 
  !>   the separate N-photon projected states onto the spatial grid `x` to show 
  !>   where single-photon vs. multi-photon bound packets are spatially localized.
  !>   (Note: Output flags suggest this routine is currently under development).
  !> Arguments:
  !>   - sys        : Parameter structure defining the spatial grid `dx`.
  !>   - st         : Current state object.
  !>   - label      : String tag appended to the output filename.
  !>   - max_number : Optional upper limit for the N-photon subspace evaluation.
  !>
  SUBROUTINE print_nphx_eo(sys,st,label,max_number)

	 type(param), intent(in)      		::  sys
	 type(state), intent(in)				::  st
	 character(len=5), intent(in)		   ::  label
	 integer, intent(in),optional			::  max_number
	 type(state)								::  upst
	 character(len=300)					   ::  name_numbst
	 integer							 			::  i1
	 complex(cx)								::  one_ph_amp( -sys%nmode+1 : sys%nmode )
	 real(rl), dimension(-sys%nmode+1 : sys%nmode)	::  tmp1, tmp2, tmp3, tmp4, tmp5
	 real(rl)									::  sum_12345

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)

	 name_numbst=trim(adjustl(sys%file_path))//"/nphx_"//trim(adjustl(label))//"_"&
		//trim(adjustl(parameterchar(sys)))//".d"
	 open (unit=100,file= name_numbst,action="write",status="replace")


	 tmp2(:) = 0._rl
	 tmp3(:) = 0._rl
	 tmp4(:) = 0._rl
	 tmp5(:) = 0._rl

	 print*,"-- ROUTINE REQUIRES UPDATE"

	 one_ph_amp = one_photon_x_amp_up(sys,upst)
	 tmp1 = abs( one_ph_amp )**2 / sys%dx
	 if ( present( max_number ) .and. ( max_number > 2) ) then
		tmp2 = nx_2_photon_up(sys,upst)
	 end if
	 if ( present( max_number ) .and. ( max_number > 3) ) then
		tmp3 = nx_3_photon_up(sys,upst)
	 end if

	 !if (present( max_number ) .and. ( max_number > 3) .and. (sys%k0<0.18)) then

	 !  print*,"      4-photon states START"
	 !  call system('date' )
	 !  tmp4 = nk_4_photon_up(sys,upst)
	 !  print*,"      4-photon states DONE"
	 !  call system('date' )

	 !  !k_del = int ( 0.08_rl/sys%dk1 )
	 !  !tmp4_del_del = nk_4_photon_up_k1_k2(sys,upst,k_del , k_del)

	 !  print*,"-- Spectrum printed up to 4 photons"

	 !end if

	 !if (present( max_number ) .and. ( max_number > 4) .and. (sys%k0<0.18) ) then

	 !  k_5_near = 0.06
	 !  n_k5 = int(k_5_near/sys%dk1) 
	 !  k_5 = sys%w(n_k5)
	 !  print*,"      5-photon states START"
	 !  tmp5(n_k5) = nk_5_photon_at_k_up(sys,st,n_k5)
	 !  print*,"      5-photon states DONE"

	 !end if

	 !do k=-sys%nmode+1, 0
	 !	w(k) = -sys%w(abs(k-1))
	 !end do
	 !do k=1, sys%nmode
	 !  w(k) = sys%w(k)
	 !end do

	 do i1= -sys%nmode+1, sys%nmode
	 	sum_12345 = tmp1(i1) + tmp2(i1) + tmp3(i1) + tmp4(i1) + tmp5(i1) 
		write(100,*) i1*sys%dx-0.5_rl, tmp1(i1), tmp2(i1), tmp3(i1)!, tmp4(i1), tmp5(i1), sum_12345

	 end do
	 close(100)

	 print*,"-- Spectrum photonic decomposition printed to file"
	 print*," "

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: destroy_photon
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Applies the real-space bosonic annihilation operator `a(x)` to the current 
  !>   variational state at a specific coordinate `x`. Modifies the qubit/atomic 
  !>   amplitudes `p(n)` across all coherent polarons based on their local spatial 
  !>   field amplitudes `fnx`.
  !> Arguments:
  !>   - sys : Parameter structure defining the grid.
  !>   - st  : State object evaluated and updated in-place.
  !>   - x   : The spatial coordinate targeted for photon annihilation.
  !>
  SUBROUTINE destroy_photon(sys,st,x)

	 type(param), intent(in)	   ::  sys
	 type(state), intent(in out)	::  st
	 type(state)					   ::  upst
	 real(rl), intent(in)		   ::  x
	 complex(cx)					   ::  fnx(sys%npini+sys%npadd,-sys%nmode+1:sys%nmode)
	 integer								::  n

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)

	 fnx = f_nx_eo(sys,st)

	 do n=1, st%np
		st%p(n) = st%p(n) * fnx( n , int(x/sys%dx) )
	 end do


  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: destroy_photon_2
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Applies a spatially smeared annihilation operator over a region defined by 
  !>   `delta_x` centered at coordinate `x`. After localized photon destruction, 
  !>   it immediately renormalizes the full multi-polaron state to restore valid 
  !>   quantum probability amplitudes.
  !> Arguments:
  !>   - sys     : Parameter structure.
  !>   - st      : State object evaluated and updated in-place.
  !>   - x       : Central spatial coordinate for photon annihilation.
  !>   - delta_x : Integration/smearing width over which the photon is destroyed.
  !>
  SUBROUTINE destroy_photon_2(sys,st,x,delta_x)

	 type(param), intent(in)	::  sys
	 type(state), intent(in out)	::  st
	 type(state)					::  upst
	 real(rl), intent(in)		::  delta_x,x
	 complex(cx)					::  fnx(sys%npini+sys%npadd,-sys%nmode+1:sys%nmode)
	 integer							::  imax, imin
	 real(rl)						::  norm_var

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)

	 fnx = f_nx_eo(sys,st)

	 imin = int( (x - 0.5*delta_x)/sys%dx )
	 imax = int( (x + 0.5*delta_x)/sys%dx )

	 st%p(:) = st%p(:) * fnx( : , int(x/sys%dx) )
	 norm_var= norm(st)
	 st%p(:) = st%p(:) / norm_var

  END SUBROUTINE
		  
!> -------------------------------------------------------------------------
  !> FUNCTION: gif_dirname
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Utility function that constructs and returns the specific directory path 
  !>   used to store the individual frames and gnuplot scripts generated during 
  !>   GIF animation routines.
  !> Arguments:
  !>   - sys : Parameter structure containing base file paths and parameter character tags.
  !> Return:
  !>   - character(len=200) : Formatted directory path string.
  !>
  FUNCTION gif_dirname(sys)

	 type(param), intent(in)			::   sys
	 character(len=200)				:: gif_dirname

	 gif_dirname=trim(adjustl(sys%file_path))//"/gif_"//trim(adjustl(parameterchar(sys)))

  END FUNCTION

!> -------------------------------------------------------------------------
  !> SUBROUTINE: open_gif_files
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Utility setup for real-time visual diagnostics. Initializes the animation 
  !>   directory, constructs the persistent Gnuplot script (`gif_fxs.gnu`) defining 
  !>   the plotting bounds and styles, and writes the initial bash script (`gif_fxs.sh`) 
  !>   for ImageMagick conversion.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>
  SUBROUTINE open_gif_files(sys)

	 type(param), intent(in)	::  sys
	 character(len=250)			::  dirname

	 dirname=trim(adjustl(gif_dirname(sys)))
	 call system('mkdir -p ./' // adjustl(trim(dirname)) )
	 open(unit=40001,file=adjustl(trim(dirname))//"/gif_fxs.sh",action="write",status="replace") 
	 open(unit=40000,file=adjustl(trim(dirname))//"/gif_fxs.gnu",action="write",status="replace") 
	 write(40001,*) 'convert -delay 200 -size 500x500 \'

	 write(40000, * ) "#!/opt/local/bin/gnuplot -persist"
	 write(40000, * ) "set term pngcairo enhanced"

	 write(40000, * ) "set style line 2 linecolor rgb '#2C6AAB' lw 3"
	 write(40000, *) "set xrange [-1350:1350]"
	 write(40000, *) "set yrange [-0.0003:0.0025]"
	 write(40000, *) "set style fill transparent solid 0.3"
	 write(40000, * ) "unset key"
	 write(40000, * ) 

	 !open(unit=20001,file=adjustl(trim(dirname))//"/gif.sh",action="write",status="replace") 
	 !open(unit=20000,file=adjustl(trim(dirname))//"/gif.gnu",action="write",status="replace") 
	 !write(20001,*) 'convert -delay 200 -size 500x500 \'
	 !open(unit=30001,file=adjustl(trim(dirname))//"/gif_wig.sh",action="write",status="replace") 
	 !open(unit=30000,file=adjustl(trim(dirname))//"/gif_wig.gnu",action="write",status="replace") 
	 !write(30001,*) 'convert -delay 200 -size 500x500 \'
	 !		  write(20000, * ) "#!/opt/local/bin/gnuplot -persist"
	 !		  write(20000, * ) "set term pngcairo enhanced"
	 !		  write(20000, * ) "set yrange [-"//trim(adjustl(alchar3))//":"//trim(adjustl(alchar3))//"]"
	 !		  write(20000, * ) "set xrange [-1:"//trim(adjustl(tmaxchar))//"]"
	 !		  write(20000, * ) 

	 !		  write(30000, * ) "#!/opt/local/bin/gnuplot -persist"
	 !		  write(30000, * ) "set term pngcairo enhanced"
	 !		  write(30000, * ) "unset key"
	 !		  write(30000, * ) "set pm3d map"
	 !		  write(30000, * ) "set size square"
	 !		  write(30000, * ) "set xrange [-2.5:2.5]"
	 !		  write(30000, * ) "set yrange [-2.5:2.5]"
	 !		  write(30000, * ) "set palette defined (-0.6 'black', -0.3 'blue', 0 '#fffaf0', 0.3 'red', 0.6 'yellow')"
	 !		  write(30000, * ) "set cbrange [-0.6:0.6]"
	 !		  write(30000, * ) 

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: end_of_GIF
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Finalizes the animation export process. Appends the ImageMagick `-loop 0` 
  !>   commands to the bash scripts to stitch the sequential PNG files into a 
  !>   continuous looping GIF, and closes out the active file handles.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>
  SUBROUTINE end_of_GIF(sys)
	 
	 type(param), intent(in)	::  sys

	 if ( sys%makegif == 1 ) then
		write(20001,*) "-loop 0 animation_fxs.gif"
		write(30001,*) "-loop 0 animation_wig.gif"
		write(40001,*) "-loop 0 animation_wig.gif"
		close(20000)
		close(20001)
		close(30000)
		close(30001)
		close(40000)
		close(40001)
	 end if

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_gif_image
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Core animation frame generator. Called periodically during the time-evolution 
  !>   loop (controlled by an adaptive threshold `t_gifplot`). Calculates the spatial 
  !>   wavepacket density via `print_nx_eo` and dynamically appends plot rendering 
  !>   commands to the active Gnuplot scripts for the current time step.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object representing the snapshot at time `t`.
  !>
  SUBROUTINE print_gif_image(sys,st)

	 type(param), intent(in)	::  sys
	 type(state), intent(in)	::  st
	 character(len=200)			::  tchar,tchar2,np2char,npchar,delaychar, &
												  tcharmx0
	 real(rl),save						::  t_gifplot

	 if ( st%t < 1.e-10_rl ) then
		t_gifplot = 0._rl
	 end if

	 if ( st%t > t_gifplot ) then

		 print*,"enter gif"
		 print*,"t=", st%t

		 write(tchar, '(f12.4)') st%t
		 write(tcharmx0, '(f12.4)') st%t-sys%x0
		 write(tchar2, '(I5)') int(st%t)
		 write(np2char, '(I10)') st%np*2
		 write(npchar, '(I10)') st%np
		 !dirname=trim(adjustl(sys%file_path))//"/gif_"//trim(adjustl(parameterchar(sys)))
		 !wig1_file=trim(adjustl(dirname))//"/t_wig1_"//trim(adjustl(tchar))//".d"
		 !wig2_file=trim(adjustl(dirname))//"/t_wig2_"//trim(adjustl(tchar))//".d"
		 !open(unit=30002,file=wig1_file,action="write",status="replace") 
		 !open(unit=30003,file=wig2_file,action="write",status="replace") 
		 !open(unit=20002,file=fxs_file,action="write",status="replace") 

		 !CALL printwigner(30002,sys,st,10,0._rl,200._rl)
		 !close(30002)
		 !CALL printwigner(30003,sys,st,10)
		 !close(30003)

		 ! do m=0,sys%nmode
		 !	write(20002,'(f8.2)',advance='no') sys%dx*dble(m)
		 !	do i=1, st%np
		 !	  write(20002,'(2f11.5)',advance='no') real( fn_x(sys,st,i,sys%dx*m) ), real( hn_x(sys,st,i,sys%dx*m) )
		 !	end do
		 !	write(20002,*)
		 ! end do
		 ! close(20002)

		 CALL print_nx_eo(sys,st,tchar2,1)

		 t_gifplot = t_gifplot + 10_rl
		! if (st%t < 1000._rl) then
		!	t_gifplot = t_gifplot + 100_rl
		! else if ( (st%t > 1000._rl) .and. (st%t < 2000._rl) ) then
		!	t_gifplot = t_gifplot + 10_rl
		! else if ( st%t > 2000._rl ) then
		!	t_gifplot = t_gifplot + 10_rl
		! end if
		 write(delaychar,'(I3)') 5

		 !write(20001,*) '-page +5+10 t_fxs_'//trim(adjustl(tchar))//'.png \'
		 !write(20001,*) '-delay '//trim(adjustl(delaychar))//' \'
		 !write(20000, *) "set output 't_fxs_"//trim(adjustl(tchar))//".png'"
		 !write(20000, *) "plot for [i=2:"//trim(adjustl(np2char))//":2]'t_fxs"//trim(adjustl(tchar))//".d' u 1:i w l"

		 !write(30001,*) '-page +5+10 t_wig_'//trim(adjustl(tchar))//'.png \'
		 !write(30001,*) '-delay '//trim(adjustl(delaychar))//' \'
		 !write(30000, *) "set output 't_wig_"//trim(adjustl(tchar))//".png'"
		 !write(30000, *) "set multiplot layout 1,2"
		 !write(30000, *) "set title 't="//trim(adjustl(tchar))//"'"
		 !write(30000, *) "splot 't_wig1_"//trim(adjustl(tchar))//".d' "
		 !write(30000, *) "splot 't_wig2_"//trim(adjustl(tchar))//".d' "
		 !write(30000, *) "unset multiplot"
		 write(40001,*) '-page +5+10 t_nx_'//trim(adjustl(tchar))//'.png \'
		 write(40001,*) '-delay '//trim(adjustl(delaychar))//' \'

		 write(40000, *) "set output 't_nx_"//trim(adjustl(tchar))//".png'"
		 write(40000, *) "set title 't="//trim(adjustl(tchar))//"'"
		 write(40000, *) "unset arrow"
		 write(40000, *) "set arrow from "//trim(adjustl(tcharmx0))//",-1 to "//trim(adjustl(tcharmx0))//",1 lc 'red'"
		 write(40000, *) "set key"
		 write(40000, *) "plot  '"//trim(adjustl(tchar2))//".d' w filledcurve ls 2"
		 !write(40000, *) "plot for [i=2:4:2] '"//trim(adjustl(tchar2))//".d' u 1:i w l title 'np='.i "

		  print*,"out gif"
	  end if

  END SUBROUTINE
  

END MODULE output
