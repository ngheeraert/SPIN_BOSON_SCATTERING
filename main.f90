PROGRAM main

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








