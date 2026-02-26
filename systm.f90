! ==============================================================================
!  FILE: systm.f90
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
!  CORE WORKFLOW
!    This module defines the primary data structures and contains the critical 
!    equations of motion (EOM) solver called by the integrator in `output.f90`. 
!    Its main responsibilities are:
!
!      1. Data Structures : Defines `param`, `state`, and `traj`.
!      2. Initialization  : Sets up ground states, single photons, and wavepackets.
!      3. EOM Evaluation  : `CalcDerivatives` computes the variational overlaps 
!                           and solves the dense linear system for parameter 
!                           time-derivatives.
!      4. Observables     : Evaluates energy, Pauli matrices, and multi-photon 
!                           correlations across momentum and spatial domains.
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

MODULE SYSTM 

  USE consts
  USE lapackmodule
  USE inverse
  USE typedefs, only : cx => c_type, rl => r_type

  IMPLICIT NONE

  TYPE param
  	 integer						 ::  fastcalcderiv, grid, back_track_on
	 integer    				 ::  npini
	 integer    				 ::  npadd
	 real(rl)    				 ::  p0 			!-- initial value for added polarons 

	 character(len=4)			 ::  file_path
	 real(rl)   				 ::  del    	!-- delta
	 integer    				 ::  nmode  	!-- number of modes
	 real(rl)					 ::  length    !-- length on which to plot the field distribution
	 real(rl)					 ::  dx     	!-- distance between 2 sites
	 real(rl)   				 ::  dt     	!-- time step
	 real(rl)   				 ::  tmax   	!-- final time
	 real(rl)   				 ::  tref   	!-- time to add polaron
	 real(rl)   				 ::  merr   	!-- time for adding additional polaron
	 real(rl)   				 ::  wc     	!-- frequency cutoff
	 real(rl)					 ::  wmax      !-- maximum frequency (cutoff at wc is exponential)
	 real(rl)					 ::  alpha     
	 real(rl), allocatable   ::  w(:)      !-- frequency spectrum
	 real(rl)					 ::  dk1
	 real(rl), allocatable   ::  g(:)      !-- coupling
	 integer						 ::  makegif
	 real(rl)		 			 ::  max_deltaE
	 real(rl)		 			 ::  xmin

	 integer						 ::  prep
	 real(rl)					 ::  k0
	 real(rl)					 ::  k0_2
	 real(rl)					 ::  x0
	 real(rl)					 ::  xg2
	 real(rl)					 ::  sigma
	 real(rl)					 ::  x1
	 real(rl)					 ::  n_wp
	 real(rl)					 ::  eta
  END TYPE param

  TYPE STATE
	 real(rl)					 	::  t 								!-- Time
	 integer    				 	::  np
	 complex(cx), allocatable  ::  f(:,:),h(:,:)   				!-- value of the fiel displacements and momenta
	 complex(cx), allocatable  ::  fo(:,:),ho(:,:)   				!-- value of the fiel displacements and momenta
	 complex(cx), allocatable  ::  p(:),q(:)   					!-- probability amplitudes of the polarons
	 complex(cx), allocatable  ::  fdot(:,:),hdot(:,:)   		!-- the time derivatives
	 complex(cx), allocatable  ::  pdot(:),qdot(:)   			!-- the time derivatives
	 !-- for storing the large sums over k
	 complex(cx), allocatable  ::  ov_ff(:,:), ov_hh(:,:)  	!-- matrices of overlaps
	 complex(cx), allocatable  ::  ov_fh(:,:), ov_hf(:,:)  	!-- matrices of overlaps
	 complex(cx), allocatable  ::  bigW_f(:,:), bigW_h(:,:)  !-- Ws (as in the notes)
	 complex(cx), allocatable  ::  bigL_f(:,:), bigL_h(:,:)  !-- Ls (as in the notes)
	 complex(cx), allocatable  ::  savepackedsol(:)
  END TYPE STATE

  TYPE TRAJ
	 integer							::  i_time,tref_counter
	 real(rl), allocatable  	::  time_ar(:), energy_ar(:), tref_ar(:,:), &
												norm_ar(:), error_ar(:), spinXYZ_ar(:,:)!, ps_ar(:,:) !,spinXYZ(:,:)
	 integer, allocatable  		::  np_ar(:)
  END TYPE TRAJ

  COMPLEX(cx), PARAMETER :: Ic = ( 0._rl , 1._rl )


CONTAINS

!> -------------------------------------------------------------------------
  !> SUBROUTINE: getParameters
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Parses command-line arguments to populate the primary `param` structure (`sys`).
  !>   Initializes the physical parameters of the spin-boson model (coupling alpha, 
  !>   detuning delta, frequency cutoff wc) and the incident wavepacket properties 
  !>   (k0, sigma, x0). Also defines the momentum grid `w(k)` and coupling `g(k)`.
  !> Arguments:
  !>   - sys : Parameter structure to be populated and allocated.
  !>
  SUBROUTINE getParameters(sys)

	 type(param)           			  ::  sys
	 character(len=100)    			  ::  buffer
	 integer 				  			  ::  i, nargs
	 character							  ::  path_input

	 !== these are the parameters provided in the command line 
	 nargs = iargc()
	 call get_command_argument(1, buffer)

	 sys%prep = 50
	 do i=1,nargs
		call get_command_argument(i, buffer)
		if(buffer=='-prep')then 
		  call get_command_argument(i+1, buffer)
		  read(buffer,*) sys%prep
		end if
	 end do

	 sys%npini= 1
	 sys%npadd= 0
	 sys%p0 = 0.000001_rl
	 sys%merr = 0.0000001_rl
	 sys%tref = 0.2_rl
	 sys%dt = 0.05_rl
	 sys%makegif = 0
	 sys%file_path = "data"
	 sys%xg2 = -600._rl

	 sys%max_deltaE= 1.0e-7_rl
	 sys%sigma = 0.005_rl
	 sys%x1 = -2250._rl
	 sys%alpha = 0.1_rl
	 sys%del = 0.1_rl
	 sys%nmode = 300
	 sys%tmax = 1300._rl
	 sys%wc = 1000._rl
	 sys%wmax = 1._rl
	 sys%k0=0.0783_rl  !-- resonance for delta=0.2, alpha=0.1 
	 sys%k0_2=0.0783_rl  !-- resonance for delta=0.2, alpha=0.1 
	 sys%x0=-500._rl
	 sys%n_wp=0.5_rl
	 sys%back_track_on=1


	 if(buffer == '--help' .or. buffer == '-h' .or. buffer=='-help')then
	 	print*, '========================================'
		print*, '-npini [2]','     	', 'Number of polarons'
		print*, '-nm [200]','      	', 'Number of modes'
		print*, '-del [0.1]','     	', 'Delta'
		print*, '-al [0.1]','      	', 'Alpha'
		print*, '-dt [0.1]','      	', 'Time step'
		print*, '-tmax [300]','    	', 'Final time'
		print*, '-wc [1000]','     	', 'Frequency cutoff'
		print*, '-wmax [1]','      	', 'Maximum frequency (smooth cutoff)'
		print*, ' '
	 	print*, '-- Polaron adding parameters ------------------------------'
		print*, '-npadd [0]','     	', 'Number of polarons to be added'
		print*, '-rtadd [1]','     	', 'Rate of adding'
		print*, '-addings [0]','     	', 'Add GS or Vaccum'
		print*, '-me ','  			', 'max error allowed before adding polaron'
		print*, '-ta2','    		', 'Add time [varies with alpha]'
		print*, '-ta2','    		', 'Add time for np>4'
		print*, '-ta3','    		', 'Add time [varies with alpha]'
		print*, '-tac1','   		', 'Change add time1 at [np pol]'
		print*, '-tac2','   		', 'Change add time2 at [np pol]'
		print*, ' '
	 	print*, '-- State preparation parameters ----------------------------'
		print*, '-A [0.01]','	      ', 'Initial amplitude of fks (if vaccum initialisation)'
		print*, '-sgn [1]','   	 	', 'Initial sign of q(:)'
		print*, '-p1i [1]','   		 ', 'Initial value for p(1), (before normalisation)'

		print*, '-prep [0]','		', '0: vaccum, 1:scat_cs, 2:scat_phton,3:get_chi'
		print*, '-getwig [0]','		', 'Outputs file with values fo the wigner function'
		print*, '-makewig [0]','	', 'Outputs file with values fo the wigner function'
	 	print*, '========================================'
		stop

	 elseif(nargs==0)then
		stop '# Use input parameters parameters'
	 else
		do i=1,nargs
		  call get_command_argument(i, buffer)
		  if(buffer=='-k0')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%k0
		  else if(buffer=='-bt')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%back_track_on
		  else if(buffer=='-k0_2')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%k0_2
		  else if(buffer=='-x0')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%x0
		  else if(buffer=='-sigma')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%sigma
		  else if(buffer=='-p0')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%p0
		  else if(buffer=='-n')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%n_wp
		  else if(buffer=='-mde')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%max_deltaE
		  else if(buffer=='-me')then
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%merr
		  elseif(buffer=='-npini')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%npini
		  elseif(buffer=='-npadd')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%npadd
		  else if(buffer=='-nm')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%nmode
		  else if(buffer=='-del')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%del
		  else if(buffer=='-al')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%alpha
		  else if(buffer=='-wc')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%wc
		  else if(buffer=='-wmax')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%wmax
		  else if(buffer=='-dt')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%dt
		  else if(buffer=='-tmax')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%tmax
		  else if(buffer=='-tref')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%tref
		  else if(buffer=='-xmin')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%xmin
		  else if(buffer=='-makegif')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%makegif
		  else if(buffer=='-x1')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%x1
		  else if(buffer=='-xg2')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) sys%xg2
		  else if(buffer=='-path')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) path_input
			 if ( path_input == "." ) then
			 	sys%file_path = "   ."
			 else if (path_input == "d") then
			 	sys%file_path = "data"
			 end if

		  end if
		
		end do

	 end if

	 allocate(sys%w(sys%nmode))
	 allocate(sys%g(sys%nmode))

	 sys%dk1 = sys%wmax / dble( sys%nmode )
	 if (sys%nmode == 1) then
		sys%w(1) = 1._rl
	 else if (sys%nmode == 2) then
		sys%w(1) = 0._rl
		sys%w(2) = sys%wmax
	 else if (sys%nmode == 5) then
		sys%w(1) = 0._rl
		sys%w(2) = 0._rl
		sys%w(3) = 0._rl
		sys%w(4) = 0._rl
		sys%w(5) = sys%wmax
	 else

		do i=1, sys%nmode
		  sys%w(i) = sys%dk1 * (i-0.5_rl)
		end do

	 end if

	 !-- define de coupling g(k)
	 sys%g(:) = sqrt( 2._rl*sys%alpha * sys%w(:) * sys%dk1 * exp(-sys%w(:)/sys%wc) ) 

	 !-- system variables that depend on a possible input
	 sys%length = pi/sys%dk1
	 sys%dx  = sys%length/sys%nmode

  END SUBROUTINE  

!> -------------------------------------------------------------------------
  !> SUBROUTINE: allocate_state
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Allocates memory for all arrays within a `state` object based on the 
  !>   number of coherent components (`np`) and the size of the momentum grid (`nmode`).
  !>   Initializes the multi-polaron amplitudes (p, q), modal arrays (f, h, fo, ho),
  !>   their time derivatives, and the extensive overlap/matrix arrays.
  !> Arguments:
  !>   - sys     : Parameter structure defining `nmode`.
  !>   - st      : The state object to be allocated and zeroed out.
  !>   - npvalue : Optional override for the number of polarons (defaults to sys%npini).
  !>
  SUBROUTINE allocate_state(sys,st,npvalue)

	 type(param), intent(in)      	::  sys
	 type(state), intent(out)  		::  st
	 type(state)							::  tmpst !-- if st has to be reallocated
	 integer, intent(in),optional		::  npvalue
	 integer									::  np

	 if (present(npvalue)) then
	 	tmpst%np = npvalue
	 else
		tmpst%np = sys%npini
	 end if
	 np = tmpst%np

	 !spinXYZ_dim = (sys%tmax/sys%dt)*2._rl
	 !allocate(tmpst%spinXYZ(spinXYZ_dim,7))

	 allocate(tmpst%f(np,sys%nmode))
	 allocate(tmpst%h(np,sys%nmode))
	 allocate(tmpst%fo(np,sys%nmode))
	 allocate(tmpst%ho(np,sys%nmode))
	 allocate(tmpst%p(np))
	 allocate(tmpst%q(np))
	 allocate(tmpst%fdot(np,sys%nmode))
	 allocate(tmpst%hdot(np,sys%nmode))
	 allocate(tmpst%pdot(np))
	 allocate(tmpst%qdot(np))
	 allocate(tmpst%ov_ff(np,np))
	 allocate(tmpst%ov_hh(np,np)) 
	 allocate(tmpst%ov_fh(np,np)) 
	 allocate(tmpst%ov_hf(np,np)) 
	 allocate(tmpst%bigW_f(np,np)) 
	 allocate(tmpst%bigW_h(np,np)) 
	 allocate(tmpst%bigL_f(np,np)) 
	 allocate(tmpst%bigL_h(np,np))

	 tmpst%t = 0._rl
	 !tmpst%spinXYZ(:,:) = 0._rl
	 tmpst%f(:,:) = 0._cx
	 tmpst%h(:,:) = 0._cx
	 tmpst%fo(:,:) = 0._cx
	 tmpst%ho(:,:) = 0._cx
	 tmpst%p(:) = 0._cx
	 tmpst%q(:) = 0._cx
	 tmpst%fdot(:,:) = 0.0_rl
	 tmpst%hdot(:,:) = 0.0_rl
	 tmpst%pdot(:) = 0.0_rl
	 tmpst%qdot(:) = 0.0_rl
	 tmpst%ov_ff(:,:) = 0._cx
	 tmpst%ov_hh(:,:) = 0._cx
	 tmpst%ov_fh(:,:) = 0._cx
	 tmpst%ov_hf(:,:) = 0._cx

	 st = tmpst

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: allocate_trajectory
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Allocates memory for tracking the system's time evolution. Dimensions 
  !>   the arrays inside the `traj` object based on the total simulation time 
  !>   (`tmax`) and integration step (`dt`), preparing arrays to log energy, 
  !>   norm, numerical error, basis size, and Pauli spin vectors.
  !> Arguments:
  !>   - sys : Parameter structure defining `tmax` and `dt`.
  !>   - tr  : The trajectory tracking object to be allocated.
  !>
  SUBROUTINE allocate_trajectory(sys,tr)

	 type(param), intent(in)      	::  sys
	 type(traj), intent(out)      	::  tr
	 integer									::  step_numb

	 print*, "-- BEGINNING TRANJECTORY ALLOCATION"

	 !spinXYZ_dim = (sys%tmax/sys%dt)*2._rl
	 !allocate(tmpst%spinXYZ(spinXYZ_dim,7))
	 step_numb = int( (sys%tmax / sys%dt)*1.1 ) 

	 allocate( tr%time_ar(step_numb) )
	 allocate( tr%error_ar(step_numb) )
	 allocate( tr%norm_ar(step_numb) )
	 allocate( tr%energy_ar(step_numb) )
	 allocate( tr%np_ar( step_numb ) )
	 allocate( tr%spinXYZ_ar(step_numb,3) )
	 allocate( tr%tref_ar(100,2) )
	 !allocate( tr%ps_ar(step_numb,sys%npini+sys%npadd) )

    tr%i_time=1
    tr%tref_counter=0
	 tr%time_ar = 0._rl
	 tr%norm_ar = 0._rl
	 tr%error_ar = 0._rl
	 tr%energy_ar = 0._rl
	 tr%np_ar = 0
	 tr%spinXYZ_ar = 0._rl
	 tr%tref_ar = 0._rl
	 !tr%ps_ar = 0._rl

	 print*, "-- TRAJECTORY ALLOCATED"

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: initialise_gs
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Initializes the system exclusively in its numerically pre-calculated 
  !>   ground state (GS) for the ultrastrong coupling regime. Reads the GS 
  !>   amplitudes and spatial profiles from an external data file generated 
  !>   by a separate minimization routine.
  !> Arguments:
  !>   - sys     : Parameter structure.
  !>   - st      : Target state to be overwritten with the GS data.
  !>   - time_in : The initial simulation start time.
  !>
  SUBROUTINE initialise_gs(sys,st,time_in)

	 type(param), intent(in)      	::  sys
	 type(state), intent(in out)  	::  st
	 real(rl), intent(in)				::  time_in
	 integer									::  i,k
	 real(rl)								::  a,b
	 logical									::  ex

	 st%t = time_in

	 inquire( file=gs_filename(sys,st), exist = ex )
	 if ( .not. ex ) then
		print*, "error: Ground State file does not exist:"
		print*, trim(gs_filename(sys,st))
		print*
		stop
	 end if
	 open(unit=1000,file=gs_filename(sys,st),status="old",action="read")
	 print*, "Initialising from: ", gs_filename(sys,st)

	 read(1000,*)
	 read(1000,'(a1)', advance="no") b !-- b is to pass the #
	 do i=1,st%np
		read(1000,'(f20.10)', advance="no") a
		st%p(i) = a
		st%q(i) = - a
	 end do
	 print*, " "
	 read(1000,*)


	 do k=1,sys%nmode
		read(1000,'(f20.10)', advance="no") b
		do i=1,st%np
		  read(1000,'(f20.10)', advance="no") a
		  st%f(i,k) = a
		  st%h(i,k) = - a
		end do
		read(1000,*)
	 end do
	 close(1000)

	 !-- updating the sums over k
	 CALL update_sums(sys,st)
	 CALL normalise(st)

  END SUBROUTINE 

!> -------------------------------------------------------------------------
  !> SUBROUTINE: initialise_cs_gs
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Prepares the primary scattering configuration: an incident Gaussian 
  !>   coherent state wavepacket scattering off the fully dressed ground state 
  !>   of the qubit. Superimposes the bare incoming photon amplitudes with the 
  !>   bound virtual photon cloud of the ground state.
  !> Arguments:
  !>   - sys          : Parameter structure.
  !>   - st           : Target state initialized with the combined configuration.
  !>   - k0           : Central momentum of the incident wavepacket.
  !>   - x0           : Initial spatial center of the wavepacket.
  !>   - sigma        : Spatial spread of the wavepacket.
  !>   - nb           : Mean photon number of the coherent state.
  !>   - ini_cloud_st : Output backup containing just the bare ground state cloud.
  !>   - t_in         : Optional initial time offset.
  !>
  SUBROUTINE initialise_cs_gs(sys,st, k0, x0, sigma, nb, ini_cloud_st, t_in)

	 type(param), intent(in)   						::  sys
	 type(state),intent(in out)						::  st, ini_cloud_st
	 real(rl), intent(in)      						::  x0, sigma, k0, nb
	 real(rl), intent(in), optional  				::  t_in   !-- time of free evolution
	 real(rl) 								 				::  t   !-- time of free evolution

	 integer 												:: k,i
	 complex(cx), dimension(size( st%f,1 ),sys%nmode) :: zke, gs_f, gs_h
	 complex(cx), dimension(size( st%f,1 )) 				:: gs_p, gs_q
	 complex(cx) 											:: phi
	 real(rl) 												:: a,b
	 logical													:: ex

	 !-- t is the time of free evolution
	 if (present(t_in)) then
	 	t=t_in
	 else
	 	t=0._rl
	 end if

	 st%f(:,:) = 0._rl
	 st%h(:,:) = 0._rl

	 gs_f = 0._rl
	 gs_h = 0._rl
	 gs_p = 0._rl
	 gs_q = 0._rl
	 zke(:,:) = 0._rl

	 do i=1,size( st%f,1 )
	   zke(i,:) = zke_f(sys,sigma, k0, x0, nb)! * exp( -Ic * sys%w(:) * t )
	   !-- the second part makes the wave-packet evolve freely for time t
	   !-- Note that zke_f(sys,sigma, k0, x0 + t, nb) != zke_f(sys,sigma, k0, x0, nb) * exp( -Ic * sys%w(:) * t )
	 end do

	 open(100,file='zke.d',action='write',status='replace')
	 do i=1, sys%nmode
		write(100,*) sys%w(i), real(zke(1,i)), aimag(zke(1,i))
	 end do
	 close(100)

	 inquire( file=gs_filename(sys,st), exist = ex )
	 if ( .not. ex ) then
		print*, "error: Ground State file does not exist:"
		print*, trim(gs_filename(sys,st))
		print*
		stop
	 end if
	 print*, "Initialising from: ", gs_filename(sys,st)
	 !==> create the gaussian
	 open(unit=1000,file=gs_filename(sys,st),status="old",action="read")
	 read(1000,*)
	 read(1000,'(a1)', advance="no") b
	 do i=1,st%np
		read(1000,'(f20.10)', advance="no") a
		gs_p(i) = a
		gs_q(i) = - a
	 end do
	 read(1000,*)
	 do k=1,sys%nmode
		read(1000,'(f20.10)', advance="no") b
		do i=1,st%np
		  read(1000,'(f20.10)', advance="no") a
		  gs_f(i,k) = a
		  gs_h(i,k) = - a
		end do
		read(1000,*)
	 end do
	 close(1000)

	 !-- transfer initial to initial_cloud state for later reference
	 st%p = gs_p
	 st%q = gs_q
	 st%f = gs_f
	 st%h = gs_h
	 ini_cloud_st = st
	 CALL update_sums(sys,ini_cloud_st) 
	 CALL normalise(ini_cloud_st)

	 do i=1,size( st%f,1 )
	   phi = 0._rl
	   phi = 0.5_rl*sum( zke(i,:)*conjg(gs_f(i,:)) - conjg(zke(i,:))*gs_f(i,:) )/dble(sys%nmode)
	   st%p(i) = exp(phi)*st%p(i)

	   phi = 0.5_rl*sum( zke(i,:)*conjg(gs_h(i,:)) - conjg(zke(i,:))*gs_h(i,:) )/dble(sys%nmode)
	   st%q(i) = exp(phi)*st%q(i)
	 end do

	 !--  zk^e = z_k + z_-k
	 !-- st%f(k) = fgs(k) + fe(k)
	 !-- st%fo(k) = fo(k)

	 st%f = gs_f + zke
	 st%h = gs_h + zke
	 do i=1,size( st%f,1 )
	   st%fo(i,:) = zko_f(sys, sigma, k0, x0, nb)!* exp( -Ic * sys%w(:) * t )
	   st%ho(i,:) = zko_f(sys, sigma, k0, x0, nb)!* exp( -Ic * sys%w(:) * t )
	 end do
	 CALL update_sums(sys,st)
	 CALL normalise(st)

	 print*, "=="
	 print*, "SYSTEM PREPARED IN THE GROUND STATE "
	 print*, "COHERENT STATE PREPAERD "
	 print*, "=="

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: initialise_photon_gs
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Prepares a single-photon Fock state scattering off the dressed qubit 
  !>   ground state. Within the multi-polaron ansatz, a single photon is 
  !>   constructed by creating an antisymmetric superposition of two opposing 
  !>   coherent states in the small-amplitude limit.
  !> Arguments:
  !>   - sys          : Parameter structure.
  !>   - st           : Target state initialized with the single photon + GS.
  !>   - k0           : Central momentum of the incident wavepacket.
  !>   - x0           : Initial spatial center of the wavepacket.
  !>   - sigma        : Spatial spread of the wavepacket.
  !>   - nb           : Amplitude scaling parameter for the coherent components.
  !>   - ini_cloud_st : Output backup containing just the bare ground state cloud.
  !>
  SUBROUTINE initialise_photon_gs(sys,st, k0, x0, sigma, nb, ini_cloud_st)

	 type(param), intent(in)   						::  sys
	 type(state),intent(in out)						::  st, ini_cloud_st
	 real(rl), intent(in)      						::  x0, sigma, k0, nb

	 type(state)											::  st2
	 integer 												:: k,i
	 complex(cx), dimension(size( st%f,1 ),sys%nmode) :: zke, gs_f, gs_h
	 complex(cx), dimension(size( st%f,1 )) 				:: gs_p, gs_q
	 complex(cx) 											:: phi
	 real(rl) 												:: a,b
	 logical													::  ex

	 st%f(:,:) = 0._rl
	 st%h(:,:) = 0._rl

	 gs_f = 0._rl
	 gs_h = 0._rl
	 gs_p = 0._rl
	 gs_q = 0._rl
	 zke(:,:) = 0._rl

	 do i=1,size( st%f,1 )/2
	   zke(i,:) = zke_f(sys,sigma, k0, x0, nb)
	 end do
	 do i=1+size( st%f,1 )/2,size( st%f,1 )
	   zke(i,:) = - zke_f(sys,sigma, k0, x0, nb)
	 end do

	 st2 = st
	 st2%np = st%np/2
	 !==> create the gaussian
	 inquire( file=gs_filename(sys,st2), exist = ex )
	 if ( .not. ex ) then
		print*, "error: Ground State file does not exist:"
		print*, trim(gs_filename(sys,st2))
		print*
		stop
	 end if
	 print*, "Initialising from: ", gs_filename(sys,st2)
	 open(unit=1000,file=gs_filename(sys,st2),status="old",action="read")
	 read(1000,*)
	 read(1000,'(a1)', advance="no") b
	 do i=1,st%np/2
		read(1000,'(f20.10)', advance="no") a
		gs_p(i) = a
		gs_q(i) = - a
		gs_p(i+size( st%f,1 )/2) = - a !--for the negatively disp. need p(3)=-p(1) to get photon
		gs_q(i+size( st%f,1 )/2) = + a
	 end do
	 read(1000,*)
	 do k=1,sys%nmode
		read(1000,'(f20.10)', advance="no") b
		do i=1,st%np/2
		  read(1000,'(f20.10)', advance="no") a
		  gs_f(i,k) = a
		  gs_h(i,k) = - a
		  gs_f(i+size( st%f,1 )/2,k) = a
		  gs_h(i+size( st%f,1 )/2,k) = - a
		end do
		read(1000,*)
	 end do
	 close(1000)

	 !-- transfer initial to initial_cloud state for later reference
	 st%p = gs_p
	 st%q = gs_q
	 st%f = gs_f
	 st%h = gs_h
	 ini_cloud_st = st

	 do i=1,size( st%f,1 )/2
	   phi = 0._rl
	   phi = 0.5_rl*sum( zke(i,:)*conjg(gs_f(i,:)) - conjg(zke(i,:))*gs_f(i,:) )/dble(sys%nmode)
	   st%p(i) = exp(phi)*st%p(i)

	   phi = 0.5_rl*sum( zke(i,:)*conjg(gs_h(i,:)) - conjg(zke(i,:))*gs_h(i,:) )/dble(sys%nmode)
	   st%q(i) = exp(phi)*st%q(i)
	 end do
	 do i=1+size( st%f,1 )/2,size( st%f,1 )
	   phi = 0._rl
	   phi = 0.5_rl*sum( -zke(i,:)*conjg(gs_f(i,:)) + conjg(zke(i,:))*gs_f(i,:) )/dble(sys%nmode)
	   st%p(i) = exp(phi)*st%p(i)

	   phi = 0.5_rl*sum( -zke(i,:)*conjg(gs_h(i,:)) + conjg(zke(i,:))*gs_h(i,:) )/dble(sys%nmode)
	   st%q(i) = exp(phi)*st%q(i)
	 end do

	 !--  zk^e = z_k + z_-k
	 !-- st%f(k) = fgs(k) + fe(k)
	 !-- st%fo(k) = fo(k)

	 st%f = gs_f + zke
	 st%h = gs_h + zke
	 do i=1,size( st%f,1 )/2
	   st%fo(i,:) = zko_f(sys, sigma, k0, x0, nb)
	   st%ho(i,:) = zko_f(sys, sigma, k0, x0, nb)
	   st%fo(i+size( st%f,1 )/2,:) = - zko_f(sys, sigma, k0, x0, nb)
	   st%ho(i+size( st%f,1 )/2,:) = - zko_f(sys, sigma, k0, x0, nb)
	 end do
	 CALL update_sums(sys,st)
	 CALL normalise(st)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: initialise_cs_2freqs_gs
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Prepares a bi-chromatic (two-color) coherent wavepacket incident on the 
  !>   qubit ground state. Sums two distinct Gaussian wavepackets in momentum 
  !>   space (centered at `k0_1` and `k0_2`) before superimposing them onto the 
  !>   ground state fields. Useful for analyzing frequency mixing processes.
  !> Arguments:
  !>   - sys          : Parameter structure.
  !>   - st           : Target state initialized with the combined configuration.
  !>   - k0_1, k0_2   : Central momenta of the two wavepacket components.
  !>   - x0, sigma    : Initial spatial center and spread.
  !>   - nb           : Mean photon number amplitude.
  !>   - ini_cloud_st : Output backup containing just the bare ground state cloud.
  !>   - t_in         : Optional initial time offset.
  !>
  SUBROUTINE initialise_cs_2freqs_gs(sys,st, k0_1, k0_2, x0, sigma, nb, ini_cloud_st, t_in)

	 type(param), intent(in)   						::  sys
	 type(state),intent(in out)						::  st, ini_cloud_st
	 real(rl), intent(in)      						::  x0, sigma, k0_1, k0_2, nb
	 real(rl), intent(in), optional  				::  t_in   !-- time of free evolution
	 real(rl) 								 				::  t   !-- time of free evolution

	 integer 												:: k,i
	 complex(cx), dimension(size( st%f,1 ),sys%nmode) :: zke, gs_f, gs_h
	 complex(cx), dimension(size( st%f,1 )) 				:: gs_p, gs_q
	 complex(cx) 											:: phi
	 real(rl) 												:: a,b,dk
	 logical													:: ex

	 !-- t is the time of free evolution
	 if (present(t_in)) then
	 	t=t_in
	 else
	 	t=0._rl
	 end if

	 st%f(:,:) = 0._rl
	 st%h(:,:) = 0._rl

	 gs_f = 0._rl
	 gs_h = 0._rl
	 gs_p = 0._rl
	 gs_q = 0._rl
	 zke(:,:) = 0._rl

	 do i=1,size( st%f,1 )
	   zke(i,:) = zke_f(sys,sigma, k0_1, x0, nb) + zke_f(sys,sigma, k0_2, x0, nb)! * exp( -Ic * sys%w(:) * t )
	   !-- the second part makes the wave-packet evolve freely for time t
	   !-- Note that zke_f(sys,sigma, k0_1, x0 + t, nb) != zke_f(sys,sigma, k0_1, x0, nb) * exp( -Ic * sys%w(:) * t )
	 end do

	 dk = sys%wmax/sys%nmode
	 inquire( file=gs_filename(sys,st), exist = ex )
	 if ( .not. ex ) then
		print*, "error: Ground State file does not exist:"
		print*, trim(gs_filename(sys,st))
		print*
		stop
	 end if
	 print*, "Initialising from: ", gs_filename(sys,st)
	 !==> create the gaussian
	 open(unit=1000,file=gs_filename(sys,st),status="old",action="read")
	 read(1000,*)
	 read(1000,'(a1)', advance="no") b
	 do i=1,st%np
		read(1000,'(f20.10)', advance="no") a
		gs_p(i) = a
		gs_q(i) = - a
	 end do
	 read(1000,*)
	 do k=1,sys%nmode
		read(1000,'(f20.10)', advance="no") b
		do i=1,st%np
		  read(1000,'(f20.10)', advance="no") a
		  gs_f(i,k) = a
		  gs_h(i,k) = - a
		end do
		read(1000,*)
	 end do
	 close(1000)

	 !-- transfer initial to initial_cloud state for later reference
	 st%p = gs_p
	 st%q = gs_q
	 st%f = gs_f
	 st%h = gs_h
	 ini_cloud_st = st
	 CALL update_sums(sys,ini_cloud_st) 
	 CALL normalise(ini_cloud_st)

	 do i=1,size( st%f,1 )
	   phi = 0._rl
	   phi = 0.5_rl*sum( zke(i,:)*conjg(gs_f(i,:)) - conjg(zke(i,:))*gs_f(i,:) )/dble(sys%nmode)
	   st%p(i) = exp(phi)*st%p(i)

	   phi = 0.5_rl*sum( zke(i,:)*conjg(gs_h(i,:)) - conjg(zke(i,:))*gs_h(i,:) )/dble(sys%nmode)
	   st%q(i) = exp(phi)*st%q(i)
	 end do

	 !--  zk^e = z_k + z_-k
	 !-- st%f(k) = fgs(k) + fe(k)
	 !-- st%fo(k) = fo(k)

	 st%f = gs_f + zke
	 st%h = gs_h + zke
	 do i=1,size( st%f,1 )
	   st%fo(i,:) = zko_f(sys, sigma, k0_1, x0, nb) + zko_f(sys, sigma, k0_2, x0, nb)!* exp( -Ic * sys%w(:) * t )
	   st%ho(i,:) = zko_f(sys, sigma, k0_1, x0, nb) + zko_f(sys, sigma, k0_2, x0, nb)!* exp( -Ic * sys%w(:) * t )
	 end do
	 CALL update_sums(sys,st)
	 CALL normalise(st)

	 print*, "=="
	 print*, "SYSTEM PREPARED IN THE GROUND STATE "
	 print*, "COHERENT STATE PREPAERD "
	 print*, "=="

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: initialise_photon
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Initializes a bare single-photon Fock state wavepacket without merging it 
  !>   with the explicitly loaded qubit ground state. Built via an antisymmetric 
  !>   superposition of two coherent states. Often used for reference baseline 
  !>   calculations or specific decoupled system tests.
  !> Arguments:
  !>   - sys          : Parameter structure.
  !>   - st           : Target state initialized with the single-photon fields.
  !>   - k0, x0       : Central momentum and initial spatial center.
  !>   - sigma, nb    : Spatial spread and amplitude scalar.
  !>   - ini_cloud_st : Output backup for compatibility.
  !>
  SUBROUTINE initialise_photon(sys,st, k0, x0, sigma, nb, ini_cloud_st)

	 type(param), intent(in)   						::  sys
	 type(state),intent(in out)						::  st, ini_cloud_st
	 real(rl), intent(in)      						::  x0, sigma, k0, nb
	 integer 												:: i
	 complex(cx), dimension(size( st%f,1 ),sys%nmode) :: zke, gs_f, gs_h
	 complex(cx), dimension(size( st%f,1 )) 				:: gs_p, gs_q

	 st%f(:,:) = 0._rl
	 st%h(:,:) = 0._rl

	 gs_f = 0._rl
	 gs_h = 0._rl
	 gs_p = 0._rl
	 gs_q = 0._rl
	 zke(:,:) = 0._rl

	 do i=1,size( st%f,1 )/2
	   zke(i,:) = zke_f(sys,sigma, k0, x0, nb)
	 end do
	 do i=1+size( st%f,1 )/2,size( st%f,1 )
	   zke(i,:) = - zke_f(sys,sigma, k0, x0, nb)
	 end do

	 !-- transfer initial to initial_cloud state for later reference
	 st%p(1) = 1._rl
	 st%p(2) = -1._rl
	 st%q(1) = 1._rl
	 st%q(2) = -1._rl
	 st%f = gs_f
	 st%h = gs_h
	 ini_cloud_st = st

	 st%f = gs_f + zke
	 st%h = gs_h + zke
	 do i=1,size( st%f,1 )/2
	   st%fo(i,:) = zko_f(sys, sigma, k0, x0, nb)
	   st%ho(i,:) = zko_f(sys, sigma, k0, x0, nb)
	   st%fo(i+size( st%f,1 )/2,:) = - zko_f(sys, sigma, k0, x0, nb)
	   st%ho(i+size( st%f,1 )/2,:) = - zko_f(sys, sigma, k0, x0, nb)
	 end do
	 CALL update_sums(sys,st)
	 CALL normalise(st)
	 print*,"here"

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: initialise_from_file_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Restores a previously saved simulation state from disk. Specifically 
  !>   designed to ingest data formatted in the Even/Odd (EO) parity momentum 
  !>   basis (`fks_fst` and `ps_fst` files). Reconstructs the symmetric and 
  !>   antisymmetric arrays (`f, h, fo, ho`) and normalizes the restored state.
  !> Arguments:
  !>   - sys : Parameter structure (dictates file path and parameter tag).
  !>   - st  : Target state to be populated with the loaded data.
  !>
  SUBROUTINE initialise_from_file_eo(sys,st)

    type(param), intent(in)			  	::  sys
	 type(state), intent(in out)	  		::  st
	 integer  									::  k,i,npol, io, items
	 character(len=200)						::  fks_file,ps_file
	 real(rl)									::  f_r,f_i,h_r,h_i,p_r,q_r,p_i,q_i,a
	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode)	::  fks, hks

	 print*, "Initialising from: ", parameterchar(sys)
	 fks_file=trim(adjustl(sys%file_path))//"/fks_fst_"//trim(adjustl(parameterchar(sys)))//".d"
	 ps_file=trim(adjustl(sys%file_path))//"/ps_fst_"//trim(adjustl(parameterchar(sys)))//".d"

	 !== reading the number of coherent states from the file
	 items=0
	 open (unit=101,file=fks_file,action="read",status="old")
	 DO
		READ(101,'(f25.15)',iostat=io,advance='no') f_r
		IF (io/=0) EXIT
		items = items + 1
	 END DO
	 items=items-1 ! accounting for the last line
	 CLOSE (101)
	 npol=(items/4)
	 print*,"-- Number of coherent states in file= ",npol

	 open (unit=101,file=ps_file,action="read",status="old")
	 open (unit=100,file=fks_file,action="read",status="old")
	 fks = 0._rl
	 hks = 0._rl
	 do  k=-sys%nmode+1,sys%nmode
		read(100,'(f25.15)',advance='no') a
	   do i=1,npol
	     read(100,'(2f25.15)',advance='no') f_r, h_r
	     fks(i,k) = f_r
	     hks(i,k) = h_r
	   end do
	   do i=1,npol
	     read(100,'(2f25.15)',advance='no') f_i, h_i
	     fks(i,k) = fks(i,k) + Ic*f_i
	     hks(i,k) = hks(i,k) + Ic*h_i
	   end do
	   read(100,*)
	 end do

	 do  k=1,sys%nmode
		st%f(:,k) = sqrt(0.5_rl)*( fks(:,k) + fks(:,-k+1) )
		st%h(:,k) = sqrt(0.5_rl)*( hks(:,k) + hks(:,-k+1) )
		st%fo(:,k) = sqrt(0.5_rl)*( fks(:,k) - fks(:,-k+1) )
		st%ho(:,k) = sqrt(0.5_rl)*( hks(:,k) - hks(:,-k+1) )
	 end do

	 do i=1,npol
		read(101,'(2f25.15)',advance='no') p_r, q_r
		st%p(i) = p_r
		st%q(i) = q_r
	 end do
	 do i=1,npol
		read(101,'(2f25.15)',advance='no') p_i, q_i
		st%p(i) = st%p(i) + Ic*p_i
		st%q(i) = st%q(i) + Ic*q_i
	 end do

	 st%t = sys%tmax

	 !-- updating the sums over k
	 CALL update_sums(sys,st)
	 CALL normalise(st)

	 close(100)
	 close(101)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: initialise_from_files_eo_2
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   An extended state-restore routine that reads both the Even/Odd state data 
  !>   and the exact system parameters from separate files specified via command 
  !>   line arguments (`-sys`, `-fks`, `-ps`). It dynamically adjusts basis sizes 
  !>   (`npadd`) before allocating memory and loading the state, allowing for 
  !>   flexible simulation continuation.
  !> Arguments:
  !>   - sys : Parameter structure (overwritten by the parsed parameter file).
  !>   - st  : Target state to be allocated and populated with the loaded data.
  !>
  SUBROUTINE initialise_from_files_eo_2(sys,st) 
    type(param), intent(in out)			::  sys
	 type(state), intent(in out)	  		::  st
	 integer  									::  k,i,npol, io, items
	 character(len=300)						::  fks_file,ps_file, sys_file, buffer, char1
	 character(len=10)						::  char2
	 real(rl)									::  f_r,f_i,h_r,h_i,p_r,q_r,p_i,q_i,a
	 integer										::  nargs
	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode)	::  fks, hks

	 nargs = iargc()
	 call get_command_argument(1, buffer)

	 if(nargs==0)then
		stop '# NO fILES'
	 else
		do i=1,nargs
		  call get_command_argument(i, buffer)
		  if(buffer=='-sys')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) char1
			 sys_file=trim(adjustl( char1 ))
		  else if(buffer=='-fks')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) char1
			 fks_file=trim(adjustl( char1 ))
		  else if(buffer=='-ps')then 
			 call get_command_argument(i+1, buffer)
			 read(buffer,*) char1
			 ps_file=trim(adjustl( char1 ))
		  end if
		end do
	 end if

	 open (unit=102,file=sys_file,action="read",status="old")

	 read(102,'(a10,I3)'    ) char2,  sys%npini
	 read(102,'(a10,I3)'    ) char2,  sys%npadd
	 read(102,'(a10,f25.15)') char2,  sys%tref
	 read(102,'(a10,f25.15)') char2,  sys%merr
	 read(102,'(a10,f25.15)') char2,  sys%p0
	 read(102,'(a10,f25.15)') char2,  sys%k0
	 read(102,'(a10,f25.15)') char2,  sys%x0
	 read(102,'(a10,f25.15)') char2,  sys%sigma
	 read(102,'(a10,f25.15)') char2,  sys%n_wp
	 read(102,'(a10,f25.15)') char2,  sys%alpha
	 read(102,'(a10,f25.15)') char2,  sys%del
	 read(102,'(a10,f25.15)') char2,  sys%dt
	 read(102,'(a10,f25.15)') char2,  sys%tmax
	 read(102,'(a10,I3)'    ) char2,  sys%prep
	 read(102,'(a10,I5)'    ) char2,  sys%nmode

	 print*, "Initialising from: ", parameterchar(sys)
	 
	 sys%prep= sys%prep-2700

	 !== reading the number of coherent states from the file
	 items=0
	 open (unit=101,file=fks_file,action="read",status="old")
	 DO
		READ(101,'(f25.15)',iostat=io,advance='no') f_r
		IF (io/=0) EXIT
		items = items + 1
	 END DO
	 items=items-1 ! accounting for the last line
	 CLOSE (101)
	 npol=(items/4)
	 print*,"-- Number of coherent states in file= ",npol

	 sys%npadd = npol - sys%npini
	 CALL allocate_state(sys,st)


	 open (unit=101,file=ps_file,action="read",status="old")
	 open (unit=100,file=fks_file,action="read",status="old")
	 print*,"here"
	 fks = 0._rl
	 hks = 0._rl
	 do  k=-sys%nmode+1,sys%nmode
		read(100,'(f25.15)',advance='no') a
	   do i=1,npol
	     read(100,'(2f25.15)',advance='no') f_r, h_r
	     fks(i,k) = f_r
	     hks(i,k) = h_r
	   end do
	   do i=1,npol
	     read(100,'(2f25.15)',advance='no') f_i, h_i
	     fks(i,k) = fks(i,k) + Ic*f_i
	     hks(i,k) = hks(i,k) + Ic*h_i
	   end do
	   read(100,*)
	 end do
	 print*,"here1"

	 do  k=1,sys%nmode
		st%f(:,k) = sqrt(0.5_rl)*( fks(:,k) + fks(:,-k+1) )
		st%h(:,k) = sqrt(0.5_rl)*( hks(:,k) + hks(:,-k+1) )
		st%fo(:,k) = sqrt(0.5_rl)*( fks(:,k) - fks(:,-k+1) )
		st%ho(:,k) = sqrt(0.5_rl)*( hks(:,k) - hks(:,-k+1) )
	 end do
	 print*,"here2"

	 do i=1,npol
		print*,"here2.2"
		read(101,'(2f25.15)',advance='no') p_r, q_r
		print*,"here2.3"
		print*,p_r
		print*,"here2.4"
		print*,st%p(i)
		print*,"here2.5"
		st%p(i) = p_r
		print*,"here2.6"
		st%q(i) = q_r
	 end do
	 print*,"here3"
	 do i=1,npol
		read(101,'(2f25.15)',advance='no') p_i, q_i
		st%p(i) = st%p(i) + Ic*p_i
		st%q(i) = st%q(i) + Ic*q_i
	 end do
	 print*,"here4"

	 st%t = sys%tmax

	 !-- updating the sums over k
	 CALL update_sums(sys,st)
	 CALL normalise(st)

	 close(100)
	 close(101)

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> FUNCTION: zke_f
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the symmetric (Even parity) momentum-space amplitude of an 
  !>   incident Gaussian wavepacket. It computes the superposition of right- 
  !>   and left-propagating components: z_k^e = z_k + z_{-k}. 
  !> Arguments:
  !>   - sys   : Parameter structure defining the momentum grid (`w` and `dk1`).
  !>   - sigma : Spatial spread of the wavepacket in real space.
  !>   - k0    : Central momentum of the wavepacket.
  !>   - x0    : Initial spatial center coordinate.
  !>   - nb    : Mean photon number amplitude.
  !> Return:
  !>   - complex(cx) : Array of size `sys%nmode` containing the symmetric amplitudes.
  !>
  FUNCTION zke_f(sys, sigma, k0 , x0, nb) result(zkout)

	 !--  this function calculates: zk^e = z_k + z_-k

	 type(param),intent(in) 	::  sys
	 real(rl), intent(in) 		::  sigma, k0, x0, nb
	 complex(cx)          		::  zkout(sys%nmode)

	 complex(cx)   ::  prefac
	 real(rl)      ::  k
	 integer       ::  i

	 prefac = 0._rl
	 prefac  = sqrt(nb)*((1._rl/(2._rl*pi*sigma**2))**0.25_rl)*exp(-Ic*0.5_rl*k0*x0)
	 zkout = 0._rl

	 do i=1, sys%nmode

		k = sys%w(i)
		zkout(i) =  sqrt(0.5_rl) *prefac* sqrt(sys%dk1) * ( &
								exp(-Ic*(k-k0)*x0) * exp(-(k-k0)**2 / (4._rl*sigma**2)) &
								+ exp(-Ic*(-k-k0)*x0) * exp(-(-k-k0)**2/(4._rl*sigma**2)) )
	 end do

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: zko_f
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the antisymmetric (Odd parity) momentum-space amplitude of an 
  !>   incident Gaussian wavepacket. It computes the difference between right- 
  !>   and left-propagating components: z_k^o = z_k - z_{-k}.
  !> Arguments:
  !>   - sys   : Parameter structure.
  !>   - sigma : Spatial spread of the wavepacket in real space.
  !>   - k0    : Central momentum of the wavepacket.
  !>   - x0    : Initial spatial center coordinate.
  !>   - nb    : Mean photon number amplitude.
  !> Return:
  !>   - complex(cx) : Array of size `sys%nmode` containing the antisymmetric amplitudes.
  !>
  FUNCTION zko_f(sys, sigma, k0 , x0, nb) result(zkout)

	 !--  this function calculates: zk^e = z_k + z_-k

	 type(param),intent(in) 	::  sys
	 real(rl), intent(in) 		::  sigma, k0, x0, nb
	 complex(cx)          		::  zkout(sys%nmode)

	 complex(cx)   ::  prefac
	 real(rl)      ::  k
	 integer       :: i

	 prefac = 0._rl
	 prefac  = sqrt(nb)*((1._rl/(2._rl*pi*sigma**2))**0.25_rl)*exp(-Ic*0.5_rl*k0*x0)
	 zkout = 0._rl

	 do i=1, sys%nmode

		k = sys%w(i)
		zkout(i) = sqrt(0.5_rl)*prefac* sqrt(sys%dk1) * ( &
		exp(-Ic*(k-k0)*x0) * exp(-(k-k0)**2/(4._rl*sigma**2)) &
		- exp(-Ic*(-k-k0)*x0) * exp(-(-k-k0)**2/(4._rl*sigma**2)) )
	 end do

  END FUNCTION

!> -------------------------------------------------------------------------
	 !> FUNCTION: parameterchar
	 !> -------------------------------------------------------------------------
	 !> Purpose / context:
	 !>   Utility function that constructs a unique string identifier based on the 
	 !>   current physical and numerical parameters (e.g., alpha, delta, nmode, dt). 
	 !>   This string is appended to output filenames to ensure data generated 
	 !>   during parameter sweeps is properly labeled and organized without overwriting.
	 !> Arguments:
	 !>   - sys : Parameter structure.
	 !> Return:
	 !>   - character(len=200) : The concatenated parameter string.
	 !>
	 FUNCTION parameterchar(sys)

		type(param), intent(in)		::   sys
		character(len=100)      		:: delchar,alChar,npiniChar,nmChar,npaddChar,&
														tmaxchar,dtchar,&
														merrchar,trefchar, &
														p0char, bt_char, mde_char, &
														wmaxchar,wcchar,k0char,nchar,&
														x0char,prepchar,sigchar
		character(len=200)				:: parameterChar, addchar, scatchar

		write(delchar, '(f6.1)') sys%del
		write(sigchar, '(f9.4)') sys%sigma
		write(alChar, '(f8.4)') sys%alpha
		write(npinichar, '(i2)') sys%npini
		write(npaddchar, '(i2)') sys%npadd
		write(p0char, '(f6.3)') sys%p0*1000000._rl
		write(merrchar, '(f6.3)') sys%merr*1000000._rl
		write(nmChar, '(I5)') sys%nmode
		write(tmaxchar, '(I10)') int(sys%tmax)
		write(mde_char, '(f6.3)') sys%max_deltaE*1000000._rl
		write(trefchar, '(f7.2)') sys%tref
		write(dtchar, '(f8.4)') sys%dt
		write(wmaxchar, '(I4)') int(sys%wmax)
		write(wcchar, '(I4)') int(sys%wc)
		write(k0char, '(f10.4)') sys%k0
		write(x0char, '(I8)') int(abs(sys%x0))
		write(nchar, '(f6.3)') sys%n_wp
		write(prepChar, '(I2)') sys%prep
		write(bt_char, '(I1)') sys%back_track_on

		addchar="_"
		scatchar="_"
		if (sys%npadd .ne. 0) then
		  addchar="_tr"//trim(adjustl(trefchar))//&
						 !"_"//trim(adjustl(tac1char))//&
						 "_me"//trim(adjustl(merrchar))//&
						 "_bt"//trim(adjustl(bt_char))//&
						 "_p"//trim(adjustl(p0char))
		end if
		if ( (sys%prep .ge. 50) .and. (sys%prep < 100) ) then
		  scatchar="_k"//trim(adjustl(k0char))//&
					 "_x"//trim(adjustl(x0char))//&
					 "_sig"//trim(adjustl(sigchar))//&
					 "_n"//trim(adjustl(nchar))
		end if

		parameterchar=trim(adjustl(nmchar))//"m"//&
					 "_np"//trim(adjustl(npinichar))//&
					 "_"//trim(adjustl(npaddchar))//&
					 "_al"//trim(adjustl(alchar))//&
					 "_del"//trim(adjustl(delchar))//&
					 "_dt"//trim(adjustl(dtchar))//&
					 trim(adjustl(addchar))//&
					 trim(adjustl(scatchar))//&
					 "tmax"//trim(adjustl(tmaxchar))//&
					 "_mde"//trim(adjustl(mde_char))//&
					 "_p"//trim(adjustl(prepchar))

	 END FUNCTION

!> -------------------------------------------------------------------------
  !> SUBROUTINE: print_param
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Writes out the complete set of simulation parameters to a dedicated text file 
  !>   (usually prefixed with `sys_`). This ensures absolute reproducibility of the 
  !>   data, allowing post-processing scripts to automatically parse the system 
  !>   configuration alongside the physical observables.
  !> Arguments:
  !>   - sys : Parameter structure to be exported.
  !>
  SUBROUTINE print_param(sys)
  
	 type(param), intent(in)	::  sys
	 character(len=200)			::  name_sys

    name_sys=trim(adjustl(sys%file_path))//"/sys_"//trim(adjustl(parameterchar(sys)))//".d"

	 open (unit=100,file= name_sys,action="write",status="replace")
 
	 write(100,'(a10,I3)'    ) "npini=    ",  sys%npini
	 write(100,'(a10,I3)'    ) "npadd=    ",  sys%npadd
	 write(100,'(a10,f25.15)') "tref=     ",  sys%tref
	 write(100,'(a10,f25.15)') "merr=     ",  sys%merr
	 write(100,'(a10,f25.15)') "p0=       ",  sys%p0
	 write(100,'(a10,f25.15)') "k0=       ",  sys%k0
	 write(100,'(a10,f25.15)') "x0=       ",  sys%x0
	 write(100,'(a10,f25.15)') "sigma=    ",  sys%sigma
	 write(100,'(a10,f25.15)') "n=		  ",  sys%n_wp
	 write(100,'(a10,f25.15)') "alpha=    ",  sys%alpha
	 write(100,'(a10,f25.15)') "delta=    ",  sys%del
	 write(100,'(a10,f25.15)') "dt=		  ",  sys%dt
	 write(100,'(a10,f25.15)') "tmax=     ",  sys%tmax
	 write(100,'(a10,I3)'    ) "prep=     ",  sys%prep
	 write(100,'(a10,I5)'    ) "nmodes=   ",  sys%nmode

	 close(100)


  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> FUNCTION: gs_filename
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Constructs the expected file path and name for the pre-calculated ground 
  !>   state data. The filename format rigidly depends on the coupling `alpha`, 
  !>   detuning `delta`, basis size `pols`, and cutoff constraints, ensuring the 
  !>   scatterer is initialized with a ground state matching the current parameters.
  !> Arguments:
  !>   - sys     : Parameter structure.
  !>   - st      : State structure determining the number of polarons (`np`).
  !>   - path_in : Optional override for the target directory (defaults to `gs_data`).
  !> Return:
  !>   - character(len=200) : The fully qualified ground state file path.
  !>
  FUNCTION gs_filename(sys,st,path_in)

	 type(param), intent(in)   		::  sys
	 type(state), intent(in)			::  st
	 character, intent(in),optional 	::  path_in
	 character(len=200)				 	::  gs_path
	 character(len=10)					::  alchar,delchar,npchar,nmchar,wmchar,wcchar
	 character(len=200)					::  gs_filename

	 gs_path = "gs_data"
	 if (present(path_in)) then
	 	gs_path = path_in
	 end if

	 write(alchar,'(f6.3)') sys%alpha
	 write(delchar,'(f5.2)') sys%del
	 write(npchar,'(I2)') st%np
	 write(nmchar,'(I7)') sys%nmode
	 write(wcchar,'(I4)') int(sys%wc+0.0001)		
	 write(wmchar,'(I4)') int(sys%wmax+0.0001)

	 gs_filename= trim(adjustl(gs_path))//"/FinalFksP_al"//trim(adjustl(alchar))//&
							 "_del"//trim(adjustl(delchar))//&
							 "_pols"//trim(adjustl(npchar))//"_nk"//trim(adjustl(nmchar))//&
							 "_wm"//trim(adjustl(wmchar))//"_wc"//trim(adjustl(wcchar))//".d"

  END FUNCTION

!> -------------------------------------------------------------------------
 !> SUBROUTINE: CalcDerivatives
 !> -------------------------------------------------------------------------
 !> Purpose / context:
 !>   The core variational equations of motion (EOM) solver. Computes the time 
 !>   derivatives of all variational parameters (f, h, p, q) by constructing the 
 !>   multi-polaron overlap matrices, building the right-hand side (RHS) driving 
 !>   terms, and solving the resulting dense linear system via LAPACK matrix inversion.
 !>   Supports both direct inverse and optimized "super inverse" methods.
 !> Arguments:
 !>   - sys              : Parameter structure.
 !>   - st               : State object (updated in place with fDot, hDot, pDot, qDot).
 !>   - superInverseFlag : Boolean to toggle the block-matrix optimized inverse routine.
 !>
 SUBROUTINE CalcDerivatives(sys,st,superInverseFlag)

  	type(param), intent(in)                         :: sys
  	type(state), intent(in out)                     :: st
    logical                                         :: superInverseFlag
  	complex(cx), dimension(st%np)                   :: bigP,bigQ
  	complex(cx), dimension(st%np, sys%nmode)        :: bigF, bigH
    complex(cx), dimension(st%np,st%np)             :: inv_ov_ff,inv_ov_hh
  	complex(cx), dimension(st%np,st%np)             :: a_f, a_h, b_f, b_h
  	complex(cx), dimension(st%np, sys%nmode)        :: rhsA_f, rhsA_h
  	complex(cx), dimension(st%np**2)                :: packed_dRHS_f, packed_dRHS_h, packedSol_f, packedSol_h, solsave_f, solsave_h
    complex(cx), dimension(st%np,st%np)             :: tempMatrix_f, tempMatrix_h, dRHS_f, dRHS_h
    complex(cx), dimension(st%np)                   :: tempMatrixTer_f, tempMatrixTer_h
    complex(cx), dimension(st%np, sys%nmode)        :: tempfDot, temphDot
    complex(cx), dimension(st%np,st%np,st%np)       :: alphaT_f,alphaT_h
    integer                                         :: info,i,j,k,n,m,l

    integer,save                                    :: counter
    real(rl)                                        :: startTime, endTime


    !print*, '-----------------------------------'
  	!print*, "fastCalcDerivatives used ! number : "
    ! print*, counter
    ! if (counter > 100000) then
    !   stop "counter>100000 !"
    ! end if

  	!-- initialisations
    !print*, 'initialisations'
    info=0
    do i=1,st%np
     bigP(i) = P_j(sys,st,i)
     bigQ(i) = Q_j(sys,st,i)
     bigF(i,:) =  F_j(sys,st,i)
     bigH(i,:) =  H_j(sys,st,i)
    end do

  	!-- invert overlap matrices
    !print*, 'invert overlap matrix'
  	inv_ov_ff=st%ov_ff
  	inv_ov_hh=st%ov_hh
  	CALL invertH(inv_ov_ff,info)
  	CALL invertH(inv_ov_hh,info)

  	!-- build b matrices
    !print*, 'build b matrices'
  	b_f=matmultiply_c(CONJG(st%f),TRANSPOSE(st%f))
  	b_h=matmultiply_c(CONJG(st%h),TRANSPOSE(st%h))

  	!-- build RHS
    !print*, 'build rhs'
    do k=1, sys%nmode
      do n=1, st%np
          rhsA_f(n,k)=sum(inv_ov_ff(n,:)*(bigF(:,k)-st%f(n,k)*bigP(:)))
          rhsA_h(n,k)=sum(inv_ov_hh(n,:)*(bigH(:,k)-st%h(n,k)*bigQ(:)))
      end do
    end do

    dRHS_f=matmultiply_c(CONJG(st%f),TRANSPOSE(rhsA_f))
    dRHS_f=st%ov_ff*dRHS_f
    dRHS_f=matmultiply_c(inv_ov_ff,dRHS_f)
    dRHS_h=matmultiply_c(CONJG(st%h),TRANSPOSE(rhsA_h))
    dRHS_h=st%ov_hh*dRHS_h
    dRHS_h=matmultiply_c(inv_ov_hh,dRHS_h)

    do i=1, st%np
      do n=1, st%np
        packed_dRHS_f((n-1)*st%np+i)=dRHS_f(i,n)
        packed_dRHS_h((n-1)*st%np+i)=dRHS_h(i,n)
      end do
    end do


    !-- build alphaTensor
    do i=1, st%np
      do n=1, st%np
        do m=1, st%np
            alphaT_f(i,n,m)=sum(inv_ov_ff(i,:)*st%ov_ff(:,n)*(b_f(:,m)-b_f(:,n)) )
            alphaT_h(i,n,m)=sum(inv_ov_hh(i,:)*st%ov_hh(:,n)*(b_h(:,m)-b_h(:,n)) )
        end do
      end do
    end do

    ! -- SuperInverse detection
    if ((superInverseFlag .eqv. .false.) &
         .OR. (st%np < 12) &
         .OR. (sys%fastCalcDeriv==1)) then
      CALL directInverse(st%np, alphaT_f, packedSol_f, packed_dRHS_f)
      CALL directInverse(st%np, alphaT_h, packedSol_h, packed_dRHS_h)
    else
      CALL superInverse_f(st%np, alphaT_f, packedSol_f, packed_dRHS_f)
      CALL superInverse_h(st%np, alphaT_h, packedSol_h, packed_dRHS_h)
    end if

    !-- system unpack
    ! print*, 'unpack'
    do i=1, st%np
      do n=1, st%np
        a_f(n,i)=packedSol_f((n-1)*st%np+i)
        a_h(n,i)=packedSol_h((n-1)*st%np+i)
      end do
    end do

    a_f=matmultiply_c(st%ov_ff, a_f)
    a_f=TRANSPOSE(a_f/st%ov_ff)  !! -- This TRANSPOSE is a mystery.

    a_h=matmultiply_c(st%ov_hh, a_h)
    a_h=TRANSPOSE(a_h/st%ov_hh)  !! -- This TRANSPOSE is a mystery.

   !-- fDot and hdot extraction
   ! print*, 'fDot and hdot'

   tempMatrix_f=matmultiply_c(inv_ov_ff,st%ov_ff*TRANSPOSE(a_f))
   do i=1, st%np
     tempfDot(i,:)=SUM(tempMatrix_f(i,:))*st%f(i,:)
   end do
   st%fDot= rhsA_f-matmultiply_c(tempMatrix_f,st%f)+tempfDot
   do i=1, st%np
     st%fDot(i,:)=st%fDot(i,:)/st%p(i)
   end do


   tempMatrix_h=matmultiply_c(inv_ov_hh,st%ov_hh*TRANSPOSE(a_h))
   do i=1, st%np
     temphDot(i,:)=SUM(tempMatrix_h(i,:))*st%h(i,:)
   end do
   st%hDot= rhsA_h-matmultiply_c(tempMatrix_h,st%h)+temphDot
   do i=1, st%np
     st%hDot(i,:)=st%hDot(i,:)/st%q(i)
   end do

   !-- evaluate pDot and qDot
   ! print*, 'pdot and qdot'
   tempMatrixTer_f= MATMUL(inv_ov_ff,bigP)
   st%pDot= 0.5_rl*( (/ (a_f(n,n), n=1, st%np) /) + st%p*(/ (conjg(a_f(m,m)), m=1, st%np) /) /CONJG(st%p) )
   st%pdot=st%pDot + tempMatrixTer_f
   st%pDot=st%pDot - SUM(tempMatrix_f, dim=2)

   tempMatrixTer_h= MATMUL(inv_ov_hh,bigQ)
   st%qDot= 0.5_rl*( (/ (a_h(m,m), m=1, st%np) /) + st%q*(/ (conjg(a_h(m,m)), m=1, st%np) /) /CONJG(st%q) )
   st%qdot=st%qDot + tempMatrixTer_h
   st%qDot=st%qDot - SUM(tempMatrix_h, dim=2)

 END SUBROUTINE

!> -------------------------------------------------------------------------
  !> FUNCTION: P_j
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Evaluates the gradient of the system energy with respect to the complex 
  !>   conjugate of the qubit's "up" amplitude, p_j^*. Used to construct the 
  !>   RHS driving terms for the variational equations of motion. Includes 
  !>   contributions from the free field energy and the spin-boson coupling.
  !> Arguments:
  !>   - sys : Parameter structure (provides `delta` and coupling definitions).
  !>   - st  : Current state object containing overlap matrices and amplitudes.
  !>   - m   : Index of the specific polaron component being evaluated.
  !> Return:
  !>   - complex(cx) : The derivative magnitude evaluated at polaron `m`.
  !>
  FUNCTION P_j(sys,st,m)

	 type(param),intent(in)         ::  sys
	 type(state),intent(in)         ::  st
	 complex(cx) 					     ::  P_j
	 integer,intent(in)             ::  m     !-- index of the p
	 integer 							  ::  j

	 P_j = 0._rl
	 do j=1, size(st%p,1)
		P_j = P_j  &
		  + 0.5_cx * sys%del * st%q(j) * st%ov_fh(m,j) & 
		  + st%p(j)*st%ov_ff(m,j) * ( st%bigW_f(m,j) - 0.5_cx * st%bigL_f(m,j) )
	 end do
	 P_j = -Ic*P_j
  
  END FUNCTION 

!> -------------------------------------------------------------------------
  !> FUNCTION: Q_j
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Evaluates the gradient of the system energy with respect to the complex 
  !>   conjugate of the qubit's "down" amplitude, q_j^*. Acts as the direct 
  !>   counterpart to `P_j` for the lower two-level state, aggregating modal 
  !>   overlaps and field expectation values.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !>   - m   : Index of the specific polaron component being evaluated.
  !> Return:
  !>   - complex(cx) : The derivative magnitude evaluated at polaron `m`.
  !>
  FUNCTION Q_j(sys,st,m)

	 type(param),intent(in)  		   ::  sys
	 type(state),intent(in)          ::  st
	 complex(cx) 				         ::  Q_j
	 integer,intent(in)              ::  m     !-- index of the q
	 integer 								::  j

	 Q_j = 0._rl
	 do j=1, size(st%q,1)
		Q_j = Q_j &
		  + 0.5_rl*sys%del * st%p(j) * st%ov_hf(m,j) & 
		  + st%q(j)*st%ov_hh(m,j) * ( st%bigW_h(m,j) + 0.5_rl * st%bigL_h(m,j) )
	 end do
	 Q_j = -Ic*Q_j

  END FUNCTION 

!> -------------------------------------------------------------------------
  !> FUNCTION: F_j
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Evaluates the functional derivative of the system energy with respect 
  !>   to the continuous photonic mode amplitudes associated with the qubit 
  !>   "up" state, f_j^*(k). Produces the local momentum-space driving forces 
  !>   across all field frequencies.
  !> Arguments:
  !>   - sys : Parameter structure (provides `g(k)` and `w(k)` arrays).
  !>   - st  : Current state object.
  !>   - j   : Index of the specific polaron component being evaluated.
  !> Return:
  !>   - complex(cx) : Array of size `sys%nmode` containing the derivative across the k-grid.
  !>
  FUNCTION F_j(sys,st,j)

	 type(param),intent(in)  		   			 ::  sys
	 type(state),intent(in)          			 ::  st
	 complex(cx) 				         			 ::  F_j(sys%nmode)
	 integer,intent(in)              			 ::  j     !-- index of the f
	 integer 											 ::  m

	 F_j  = 0._cx
	 do m=1, size(st%f,1)
		F_j(:) = F_j(:) &
		+ 0.5_rl*sys%del*st%q(m)*st%h(m,:)*st%ov_fh(j,m) &
		+ st%p(m)*( st%f(m,:)*st%bigW_f(j,m) + sys%w(:)*st%f(m,:) )*st%ov_ff(j,m) &
		- 0.5_rl*st%p(m)*st%ov_ff(j,m)*( st%f(m,:)*st%bigL_f(j,m) + sys%g(:) )
	 end do
	 F_j(:) = - Ic*F_j(:) 

	 END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: H_j
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Evaluates the functional derivative of the system energy with respect 
  !>   to the continuous photonic mode amplitudes associated with the qubit 
  !>   "down" state, h_j^*(k). Acts as the direct counterpart to `F_j` for 
  !>   the lower energy manifold of the spin-boson system.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !>   - j   : Index of the specific polaron component being evaluated.
  !> Return:
  !>   - complex(cx) : Array of size `sys%nmode` containing the derivative across the k-grid.
  !>
  FUNCTION H_j(sys,st,j)

	 type(param),intent(in)  		   			 ::  sys
	 type(state),intent(in)          			 ::  st
	 complex(cx) 				         			 ::  H_j(sys%nmode)
	 integer,intent(in)              			 ::  j     !-- index of the f
	 integer 											 ::  m

	 H_j  = 0._cx
	 do m=1, size(st%f,1)
		H_j(:) = H_j(:) &
		+ 0.5_rl*sys%del*st%p(m)*st%f(m,:)*st%ov_hf(j,m) &
		+ st%q(m)*( st%h(m,:)*st%bigW_h(j,m) + sys%w(:)*st%h(m,:) )*st%ov_hh(j,m) &
		+ 0.5_rl*st%q(m)*st%ov_hh(j,m)*( st%h(m,:)*st%bigL_h(j,m) + sys%g(:) )
	 end do
	 H_j(:) = - Ic*H_j(:) 

	 END FUNCTION
	 
!> -------------------------------------------------------------------------
  !> SUBROUTINE: update_sums
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   A core performance-critical utility that precomputes and caches the extensive 
  !>   matrix elements needed for the variational equations of motion (EOM). 
  !>   It calculates the bosonic overlap matrices (`ov_ff`, `ov_hh`, `ov_fh`, `ov_hf`) 
  !>   between all polarons, as well as the weighted frequency matrices (`bigW`) 
  !>   and coupling matrices (`bigL`) representing the free-field and interaction 
  !>   energies, respectively.
  !> Arguments:
  !>   - sys : Parameter structure defining the `w` (frequency) and `g` (coupling) grids.
  !>   - st  : Target state whose matrix elements are evaluated and updated in-place.
  !>
  SUBROUTINE update_sums(sys,st)
	 
	 type(param),intent(in)			::  sys
	 type(state),intent(in out)   ::  st
	 integer								::  i,j

	 do i=1, size(st%f,1)
		do j=1, size(st%f,1)

			 st%ov_ff(i,j) = ov(st%f(i,:),st%f(j,:))
			 st%ov_hh(i,j) = ov(st%h(i,:),st%h(j,:))
			 st%ov_fh(i,j) = ov(st%f(i,:),st%h(j,:))
			 st%ov_hf(i,j) = ov(st%h(i,:),st%f(j,:))

		  if (sys%nmode .ne. 1) then
			 st%bigW_f(i,j) = dot_product( sys%w , conjg(st%f(i,:)) * st%f(j,:) )
			 st%bigW_h(i,j) = dot_product( sys%w , conjg(st%h(i,:)) * st%h(j,:) )
			 st%bigL_f(i,j) = dot_product( sys%g , conjg(st%f(i,:)) + st%f(j,:) )
			 st%bigL_h(i,j) = dot_product( sys%g , conjg(st%h(i,:)) + st%h(j,:) )
		  else if (sys%nmode == 1) then
			 st%bigW_f(i,j) = sys%w(1)*conjg(st%f(i,1)) * st%f(j,1)
			 st%bigW_h(i,j) = sys%w(1)*conjg(st%h(i,1)) * st%h(j,1)
			 st%bigL_f(i,j) = sys%g(1)*( conjg(st%f(i,1)) + st%f(j,1) )
			 st%bigL_h(i,j) = sys%g(1)*( conjg(st%h(i,1)) + st%h(j,1) )
		  end if

		end do
	 end do

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> FUNCTION: Energy
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the exact expectation value of the spin-boson Hamiltonian 
  !>   for the current multi-polaron variational state. It sums three physical 
  !>   contributions: the qubit bare energy (proportional to `delta`), the 
  !>   free bosonic field energy (using `bigW`), and the dipole interaction 
  !>   energy between the qubit and the field (using `bigL`).
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The total system energy. Throws a warning if the imaginary part is non-zero.
  !>
  FUNCTION Energy(sys,st)
	 type(param),intent(in)         				::  sys
	 type(state),intent(in)         				::  st
	 complex(cx)						  				::  tmp
	 real(rl)							  				::  energy
	 integer 							  				::  i,j

	 tmp=0._cx
	 energy=0._rl
	 
	 do i=1,size(st%p,1)
		do j=1,size(st%p,1)

		  tmp = tmp &
				  + 0.5_rl * sys%del*( conjg(st%p(i))*st%q(j)*st%ov_fh(i,j) + conjg(st%q(j))*st%p(i)*st%ov_hf(j,i) ) &
				  + ( conjg(st%p(i))*st%p(j)*st%ov_ff(i,j) * st%bigW_f(i,j) + conjg(st%q(i))*st%q(j)*st%ov_hh(i,j) * st%bigW_h(i,j) ) &
				  - 0.5_rl * ( conjg(st%p(i))*st%p(j)*st%ov_ff(i,j) * st%bigL_f(i,j) - conjg(st%q(i))*st%q(j)*st%ov_hh(i,j) * st%bigL_h(i,j) )

		end do
	 end do

	 if (aimag(tmp) > 1e-10_rl) then
	 	print*, "Error: Energy is complex"
	 end if

	 energy = real(tmp)

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: norm
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the overall L2 norm of the multi-polaron state to verify total 
  !>   probability conservation. It evaluates the double sum over all coherent 
  !>   state cross-terms, factoring in both the qubit amplitudes (`p`, `q`) 
  !>   and the full bosonic inner products (Even/Odd parity spaces).
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The scalar norm of the variational wavefunction.
  !>
  FUNCTION norm(st)

	 real(rl) 			       		::  norm
	 complex(cx)						::  tmp
	 type(state),intent(in)  		::  st
	 integer						 		::  m,n

	 tmp = 0._cx
	 do m=1,size(st%p,1)
		do n=1,size(st%p,1)
		  tmp = tmp  &
		  + conjg(st%p(m))*st%p(n)*st%ov_ff(m,n) *ov( st%fo(m,:) , st%fo(n,:) ) &
		  + conjg(st%q(m))*st%q(n)*st%ov_hh(m,n)*ov( st%ho(m,:) , st%ho(n,:) )
		!  + conjg(st%p(m))*st%p(n)*ov( sqrt(0.5_rl)*(st%f(m,:)+st%fo(m,:)) , sqrt(0.5_rl)*(st%f(n,:)+st%fo(n,:)) )&
		!							 *ov( sqrt(0.5_rl)*(st%f(m,:)-st%fo(m,:)) , sqrt(0.5_rl)*(st%f(n,:)-st%fo(n,:)) ) &
		!  + conjg(st%q(m))*st%q(n)*ov( sqrt(0.5_rl)*(st%h(m,:)+st%ho(m,:)) , sqrt(0.5_rl)*(st%h(n,:)+st%ho(n,:)) )&
		!							 *ov( sqrt(0.5_rl)*(st%h(m,:)-st%ho(m,:)) , sqrt(0.5_rl)*(st%h(n,:)-st%ho(n,:)) )
		end do
	 end do

	 norm = sqrt(  real(tmp) )

  END FUNCTION

!> -------------------------------------------------------------------------
  !> SUBROUTINE: normalise
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Enforces probability conservation by normalizing the multi-polaron state. 
  !>   It calculates the total state norm via `norm(st)` and scales the qubit 
  !>   probability amplitudes (`p` and `q`) accordingly. The bosonic displacement 
  !>   amplitudes (`f`, `h`) are left untouched, as scaling them would change the 
  !>   physical photon number rather than the state normalization.
  !> Arguments:
  !>   - st : Target state to be normalized in-place.
  !>
  SUBROUTINE normalise(st)

	 type(state), intent(in out)  ::  st
	 real(rl)							:: normval

	 normval=norm(st)

	 st%p = st%p/normval
	 st%q = st%q/normval

  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> FUNCTION: error
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Evaluates the McLachlan/Dirac-Frenkel variational error metric. This 
  !>   quantifies how well the restricted multi-polaron ansatz is capturing the 
  !>   true Schrödinger time evolution. It computes the distance between the 
  !>   projected time derivatives and the exact Hamiltonian action on the state. 
  !>   Used adaptively to trigger the addition of new polarons (basis enlargement).
  !> Arguments:
  !>   - sys  : Parameter structure.
  !>   - oost : State at t - 2*dt.
  !>   - ost  : State at t - dt.
  !>   - st   : Current state at time t.
  !> Return:
  !>   - complex(cx) : The variational error metric.
  !>
 FUNCTION error(sys,oost,ost,st)
	 
	 type(param),intent(in)		::  sys
	 type(state),intent(in) 	::  oost,ost,st
	 complex(cx)					::  tmp1,tmp2, tmp3, tmp4
	 complex(cx)	   			::  error
    complex(cx), dimension(ost%np,sys%nmode)  ::  off, ohh, ofh, ohf,ofdd,ohdd
    complex(cx), dimension(ost%np)  ::  opdd,oqdd
	 integer					   	::  i,j

	 tmp1 = 0._cx
	 tmp2 = 0._cx
	 tmp3 = 0._cx
	 tmp4 = 0._cx
    error = 0._cx

	 ofdd(:,:) = (st%fdot(:,:) - oost%fdot(:,:))/(st%t-oost%t)
	 ohdd(:,:) = (st%hdot(:,:) - oost%hdot(:,:))/(st%t-oost%t)
	 opdd(:) = (st%pdot(:) - oost%pdot(:))/(st%t-oost%t)
	 oqdd(:) = (st%qdot(:) - oost%qdot(:))/(st%t-oost%t)

	 do i=1, size(st%f,1)
		do j=1, size(ost%f,1)
		  off(i,j) = dot_product( ost%fdot(i,:) , ost%f(j,:) )
		  ohh(i,j) = dot_product( ost%hdot(i,:) , ost%h(j,:) )
		  ofh(i,j) = dot_product( ost%fdot(i,:) , ost%h(j,:) )
		  ohf(i,j) = dot_product( ost%hdot(i,:) , ost%f(j,:) )
		end do
	 end do

	 tmp3 = tmp3 + 0.25_rl*sys%del**2 + 0.25_rl*sum(sys%g(:)*sys%g(:))

	 do i=1,ost%np
		do j=1,ost%np

		  tmp1 = tmp1 + ost%ov_ff(i,j)*( &
						+ conjg(ost%pdot(i))*ost%pdot(j) - 0.5_rl*conjg(ost%pdot(i))*ost%p(j)*( conjg(off(j,j)) + off(j,j) - 2_rl*conjg(off(j,i)) )&
								- 0.5_rl*conjg(ost%p(i))*ost%pdot(j)*( conjg(off(i,i)) + off(i,i) - 2_rl*off(i,j) ) &
						+ 0.25_rl*conjg(ost%p(i))*ost%p(j)*( (conjg(off(i,i))+off(i,i))*( conjg(off(j,j)) + off(j,j) - 2_rl*conjg(off(j,i))) &
												    - 2_rl*off(i,j)*( conjg(off(j,j)) + off(j,j) ) &
													 + 4_rl*( off(i,j)*conjg(off(j,i)) + sum(conjg(ost%fdot(i,:))*ost%fdot(j,:)) ) ) )
		  
		  tmp2 = tmp2 + ost%ov_ff(i,j)*ost%p(j)* ( &
								( conjg(ost%pdot(i)) - 0.5_rl*conjg(ost%p(i))*(conjg(off(i,i)) + off(i,i)) )*( ost%bigW_f(i,j) &
																										  - 0.5_rl*ost%bigL_f(i,j) ) &
								  + conjg(ost%p(i))*off(i,j)*( ost%bigW_f(i,j) - 0.5_rl*ost%bigL_f(i,j) ) &
								  + conjg(ost%p(i))*sum( sys%w(:)*conjg(ost%fdot(i,:))*ost%f(j,:) - 0.5_rl*conjg(ost%fdot(i,:))*sys%g(:)) ) &
							 + 0.5_rl*sys%del*ost%q(j)*ost%ov_fh(i,j)*( conjg(ost%pdot(i)) & 
										  - 0.5_rl*conjg(ost%p(i))*( conjg(off(i,i)) + off(i,i) - 2_rl*ofh(i,j) ) )


		  tmp3 = tmp3 + conjg(ost%p(i))*ost%p(j)*ost%ov_ff(i,j)*( sum( sys%w(:)*sys%w(:)*conjg(ost%f(i,:))*ost%f(j,:) ) &
								+ (ost%bigW_f(i,j))**2 &
								+ 0.25_rl*(ost%bigL_f(i,j))**2 &
								- 0.5_rl*sum( sys%g(:)*sys%w(:)*(conjg(ost%f(i,:)) + ost%f(j,:)) ) &
								- ost%bigL_f(i,j)*ost%bigW_f(i,j) ) &
						  + sys%del*conjg(ost%p(i))*ost%q(j)*ost%ov_fh(i,j)*sum( sys%w(:)*conjg(ost%f(i,:))*ost%h(j,:) )

		  tmp4 = tmp4 + conjg(ost%p(j))*ost%ov_ff(j,i)*( &
										+ opdd(i) &
										- ost%pdot(i)*( off(i,i) + conjg(off(i,i)) - 2_rl*conjg(off(i,j)) ) &
										- 0.5_rl*ost%p(i)*( sum( 2_rl*real(conjg(ofdd(i,:))*ost%f(i,:)) &
																				+ 2_rl*conjg(ost%fdot(i,:))*ost%fdot(i,:) &
																				- 2_rl*ofdd(i,:)*conjg(ost%f(j,:)) ) ) &
										+ 0.25_rl*ost%p(i)*( off(i,i) + conjg(off(i,i)) - 2_rl*conjg(off(i,j)) )**2 &
										)
								

		  !-- the q and h parts 

		  tmp1 = tmp1 + ost%ov_hh(i,j)*( &
						+ conjg(ost%qdot(i))*ost%qdot(j) - 0.5_rl*conjg(ost%qdot(i))*ost%q(j)*( conjg(ohh(j,j)) + ohh(j,j) - 2_rl*conjg(ohh(j,i)) )&
								- 0.5_rl*conjg(ost%q(i))*ost%qdot(j)*( conjg(ohh(i,i)) + ohh(i,i) - 2_rl*ohh(i,j) ) &
						+ 0.25_rl*conjg(ost%q(i))*ost%q(j)*( (conjg(ohh(i,i))+ohh(i,i))*( conjg(ohh(j,j)) + ohh(j,j) - 2_rl*conjg(ohh(j,i))) &
												    - 2_rl*ohh(i,j)*( conjg(ohh(j,j)) + ohh(j,j) ) &
													 + 4_rl*( ohh(i,j)*conjg(ohh(j,i)) + sum(conjg(ost%hdot(i,:))*ost%hdot(j,:)) ) ) )
													 

		  tmp2 = tmp2 + ost%ov_hh(i,j)*ost%q(j)* ( &
								( conjg(ost%qdot(i)) - 0.5_rl*conjg(ost%q(i))*(conjg(ohh(i,i)) + ohh(i,i)) )*( ost%bigW_h(i,j) &
																										- 0.5_rl*(-ost%bigL_h(i,j)) ) &
								  + conjg(ost%q(i))*ohh(i,j)*( ost%bigW_h(i,j) - 0.5_rl*(-ost%bigL_h(i,j)) ) &
								  + conjg(ost%q(i))*sum( sys%w(:)*conjg(ost%hdot(i,:))*ost%h(j,:) - 0.5_rl*conjg(ost%hdot(i,:))*(-sys%g(:)) ) ) &
							 + 0.5_rl*sys%del*ost%p(j)*ost%ov_hf(i,j)*( conjg(ost%qdot(i)) & 
										  - 0.5_rl*conjg(ost%q(i))*( conjg(ohh(i,i)) + ohh(i,i) - 2_rl*ohf(i,j) ) )

		  tmp3 = tmp3 + conjg(ost%q(i))*ost%q(j)*ost%ov_hh(i,j)*( sum( sys%w(:)*sys%w(:)*conjg(ost%h(i,:))*ost%h(j,:) ) &
								+ (ost%bigW_h(i,j))**2 &
								+ 0.25_rl*(-ost%bigL_h(i,j))**2 &
								- 0.5_rl*sum( (-sys%g(:))*sys%w(:)*(conjg(ost%h(i,:)) + ost%h(j,:)) ) &
								- (-ost%bigL_h(i,j))*ost%bigW_h(i,j) ) &
						  + sys%del*conjg(ost%q(i))*ost%p(j)*ost%ov_hf(i,j)*sum( sys%w(:)*conjg(ost%h(i,:))*ost%f(j,:) )

		  tmp4 = tmp4 + conjg(ost%q(j))*ost%ov_hh(j,i)*( &
										+ oqdd(i) &
										- ost%qdot(i)*( ohh(i,i) + conjg(ohh(i,i)) - 2_rl*conjg(ohh(i,j)) ) &
										- 0.5_rl*ost%q(i)*( sum( 2_rl*real(conjg(ohdd(i,:))*ost%h(i,:)) &
																				+ 2_rl*conjg(ost%hdot(i,:))*ost%hdot(i,:) &
																				- 2_rl*ohdd(i,:)*conjg(ost%h(j,:)) ) ) &
										+ 0.25_rl*ost%q(i)*( ohh(i,i) + conjg(ohh(i,i)) - 2_rl*conjg(ohh(i,j)) )**2 &
										)

		end do
	 end do

	 error =  -0.5_rl*real(tmp4) + 0.5_rl*tmp1 - 2._rl*aimag(tmp2) + tmp3

  END FUNCTION	

!> -------------------------------------------------------------------------
  !> FUNCTION: ov_states
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the full quantum mechanical overlap (inner product) between 
  !>   two entirely separate multi-polaron states, `<st1 | st2>`. It loops 
  !>   over all cross-combinations of their respective coherent bases, incorporating 
  !>   both the spin (`p`, `q`) and bosonic field (`f`, `h`) overlaps.
  !> Arguments:
  !>   - st1 : First (bra) state object.
  !>   - st2 : Second (ket) state object.
  !> Return:
  !>   - complex(cx) : The complex inner product.
  !>
  FUNCTION ov_states(st1,st2)

	 type(state), intent(in)  ::  st1,st2
	 complex(cx) 				  ::  ov_states 					
	 complex(cx)				  ::  tmp
	 integer						  ::  i,j

	 tmp = 0._cx

	 do i=1,st1%np
		do j=1,st1%np
		  tmp = tmp + conjg(st1%p(i))*st2%p(j)*ov( st1%f(i,:),st2%f(j,:) ) &
		  + conjg(st1%q(i))*st2%q(j)*ov( st1%h(i,:),st2%h(j,:) )
		end do
	 end do

	 ov_states = tmp

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: ov_scalar
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the fundamental bosonic overlap between two scalar coherent 
  !>   state amplitudes. Evaluates the standard coherent state inner product 
  !>   formula: exp(-0.5*|f1|^2 - 0.5*|f2|^2 + f1^* f2).
  !> Arguments:
  !>   - f1 : First scalar complex amplitude.
  !>   - f2 : Second scalar complex amplitude.
  !> Return:
  !>   - complex(cx) : The overlap scalar.
  !>
  FUNCTION ov_scalar(f1,f2)

        complex(cx), intent(in) :: f1, f2
        complex(cx)             :: ov_scalar
        complex(cx)				  :: tmp1, tmp2, tmp3

		  tmp1 = conjg(f1)*f1
		  tmp2 = conjg(f2)*f2
		  tmp3 = conjg(f1)*f2

        ov_scalar = exp( -0.5_rl*tmp1 - 0.5_rl*tmp2 + tmp3 ) 

  END FUNCTION	ov_scalar

!> -------------------------------------------------------------------------
  !> FUNCTION: ov
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the macroscopic bosonic overlap between two multimode continuous 
  !>   coherent states. Generalizes the scalar coherent overlap to vectors of 
  !>   modes by taking the dot product across the entire momentum grid `k`.
  !> Arguments:
  !>   - f1 : First complex mode array (bra).
  !>   - f2 : Second complex mode array (ket).
  !> Return:
  !>   - complex(cx) : The total multimode coherent overlap.
  !>
  FUNCTION ov(f1,f2)

        complex(cx), intent(in) :: f1( : ), f2( : )
        complex(cx)             :: ov
        != internal variables
        complex(cx)    :: tmp1, tmp2, tmp3

        != initialize
        tmp1 = 0._rl
        tmp2 = 0._rl
        tmp3 = 0._rl

		  if (size(f1,1) .ne. 1) then
			 tmp1 = dot_product(f1, f1)
			 tmp2 = dot_product(f2, f2)
			 tmp3 = dot_product(f1, f2)
		  else if (size(f1,1) == 1) then
			 tmp1 = conjg(f1(1))*f1(1)
			 tmp2 = conjg(f2(1))*f2(1)
			 tmp3 = conjg(f1(1))*f2(1)
		  end if

        ov = exp( -0.5_rl*tmp1 - 0.5_rl*tmp2 + tmp3 ) 

  END FUNCTION	ov

!> -------------------------------------------------------------------------
  !> FUNCTION: upProb
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the reduced density matrix probability of finding the qubit 
  !>   (two-level system) in the excited "up" state. It traces out the bosonic 
  !>   environment by summing the squared magnitudes of the `p` amplitudes, 
  !>   weighted by their associated bosonic overlaps `ov_ff`.
  !> Arguments:
  !>   - st : Current state object.
  !> Return:
  !>   - real(rl) : The spin-up probability.
  !>
  FUNCTION upProb(st)

	 real(rl) 			       		::  upProb
	 type(state),intent(in)  		::  st
	 integer						 		::  m,n
	 complex(cx)						::  tmp

	 tmp = 0._rl
	 do m=1,size(st%p,1)
		do n=1,size(st%p,1)
		  tmp = tmp + conjg(st%p(m))*st%p(n)*st%ov_ff(m,n)
		end do
	 end do
	 upprob=real(tmp)

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: downProb
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the reduced density matrix probability of finding the qubit 
  !>   in the ground "down" state. Counterpart to `upProb`, it traces out the 
  !>   environment by evaluating the `q` amplitudes and the `ov_hh` overlaps.
  !> Arguments:
  !>   - st : Current state object.
  !> Return:
  !>   - real(rl) : The spin-down probability.
  !>
  FUNCTION downProb(st)

	 real(rl) 			       		::  downProb
	 type(state),intent(in)  		::  st
	 integer						 		::  m,n
	 complex(cx)						::  tmp

	 tmp = 0._rl
	 do m=1,size(st%p,1)
		do n=1,size(st%p,1)
		  tmp  = tmp + conjg(st%q(m))*st%q(n)*st%ov_hh(m,n)
		end do
	 end do
	 downprob=real(tmp)

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: sigmaX
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the expectation value of the Pauli X operator, ⟨σ_X⟩, for the 
  !>   artificial atom (qubit). It traces out the bosonic bath by evaluating 
  !>   the off-diagonal coherence terms between the qubit "up" (p) and "down" (q) 
  !>   amplitudes, weighted by the cross-state bosonic overlaps `ov_fh` and `ov_hf`.
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The expectation value of σ_X.
  !>
  FUNCTION sigmaX(st)
  
	 type(state), intent(in)  :: st
	 real(rl) 					  :: sigmaX
	 complex(cx) 				  :: tmp
	 integer						  :: n,m
	 
	 tmp = 0_cx
	 sigmaX = 0_rl

	 do n=1,size(st%p,1)
		do m=1,size(st%p,1)
		  tmp = tmp + conjg(st%p(n))*st%q(m)*st%ov_fh(n,m) &
							   + conjg(st%q(n))*st%p(m)*st%ov_hf(n,m)
		end do	
	 end do
	 sigmaX = real(tmp)
  
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: sigmaZ
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the expectation value of the Pauli Z operator, ⟨σ_Z⟩, representing 
  !>   the population inversion of the two-level system. Evaluated simply as the 
  !>   difference between the spin-up and spin-down probabilities.
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The expectation value of σ_Z (ranges from -1 to 1).
  !>
  FUNCTION sigmaZ(st)
  
	 type(state), intent(in)  :: st
	 real(rl) 					  :: sigmaZ

	 sigmaZ = upProb(st) - downProb(st)
  
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: d_sigmaZ
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the numerical time derivative of the population inversion, 
  !>   d⟨σ_Z⟩/dt, using a simple finite-difference approach between the current 
  !>   and previous time steps. Useful for tracking emission rates or Rabi oscillations.
  !> Arguments:
  !>   - ost : Old state object (at time t - dt).
  !>   - st  : Current state object (at time t).
  !>   - dt  : The integration time step width.
  !> Return:
  !>   - real(rl) : The rate of change of the population inversion.
  !>
  FUNCTION d_sigmaZ(ost,st,dt)
  
	 type(state), intent(in)  :: ost, st
	 real(rl), intent(in)	  :: dt
	 real(rl) 					  :: d_sigmaZ

	 d_sigmaZ = ( sigmaZ(st) - sigmaZ(ost) ) / dt
  
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: sigmaY
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the expectation value of the Pauli Y operator, ⟨σ_Y⟩. Similar 
  !>   to σ_X, it traces over the bosonic environment, but applies the appropriate 
  !>   complex phase (i) to the off-diagonal off-diagonal coherences between 
  !>   the "up" and "down" qubit manifolds.
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : The expectation value of σ_Y.
  !>
  FUNCTION sigmaY(st)
  
	 type(state), intent(in)  :: st
	 real(rl) 					  :: sigmaY
	 real(rl) 					  :: tmp
	 integer						  :: n,m
	 
	 tmp = 0_cx
	 sigmaY = 0_rl

	 do n=1,size(st%p,1)
		do m=1,size(st%p,1)
		  tmp = tmp + 2._rl*real( Ic * conjg(st%q(n))*st%p(m)*st%ov_hf(n,m) )
		end do	
	 end do
	 sigmaY = tmp
  
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: FT_X_to_k_EO
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Performs a spatially windowed Fourier transform from real space (x) back 
  !>   into momentum space (k) within the Even/Odd (EO) parity basis. It uses 
  !>   smoothed arctangent envelope functions to explicitly filter out the 
  !>   bound state localized near the qubit (x=0), isolating the freely propagating 
  !>   scattered wavepacket.
  !> Arguments:
  !>   - sys    : Parameter structure.
  !>   - fnx    : The real-space spatial distribution array.
  !>   - sign_k : Phase sign convention for the transform (+1 or -1).
  !>   - xmin   : Optional spatial cutoff boundary to exclude the interaction region.
  !> Return:
  !>   - complex(cx) : The filtered momentum-space array.
  !>
  FUNCTION FT_X_to_k_EO(sys,fnx,sign_k,xmin)

	 type(param),intent(in)   													::  sys
	 complex(cx), intent(in)													::  fnx(:,:)!(:,-sys%nmode+1:sys%nmode)
	 integer, intent(in)															::  sign_k
	 real(rl),intent(in),optional												::  xmin
	 real(rl)																		::  delta, delta2, delta3
	 complex(cx), dimension(size(fnx,1),-sys%nmode+1:sys%nmode) 	::  new_fnx
	 complex(cx), dimension(size(fnx,1),sys%nmode)  					::	 FT_X_to_k_EO
	 real(rl), dimension(-sys%nmode+1:sys%nmode)		  				::  x_arr,atan_arr,atan_arr2,atan_arr3
	 integer																			::  n,i, imin, i_delta

	 if (present(xmin)) then
		imin = int(xmin/sys%dx) + 1
	 else
		imin = -sys%nmode +1
	 end if

	 delta = 0.010_rl
	 delta2 = 0.01_rl
	 !delta3 = 0.0001_rl
	 delta3 = 0.001_rl
	 i_delta = int(delta/sys%dx)

	 do i=-sys%nmode+1, sys%nmode
		x_arr(i) = (pi/sys%wmax) * dble(i-0.5)
		atan_arr(i) = 1 + datan((x_arr(i) - xmin)/delta)/pi - datan((x_arr(i)+xmin)/delta)/pi
		atan_arr2(i) = 1 + datan((x_arr(i) - xmin/2._rl)/delta2)/pi - datan((x_arr(i)+xmin/2._rl)/delta2)/pi
		atan_arr3(i) = 1 + datan((x_arr(i) - xmin/4._rl)/delta3)/pi - datan((x_arr(i)+xmin/4._rl)/delta3)/pi
	 end do

	 do n=1,size(fnx,1)
		new_fnx(n,:) = fnx(n,:)*atan_arr(:)*atan_arr(:)*atan_arr2(:)*atan_arr2(:)!*atan_arr3(:)
	 end do

	 do n=1,size(fnx,1)
		do i= 1, sys%nmode
		  FT_X_to_K_EO(n,i) = sqrt(0.5_rl/pi) * sqrt(sys%dk1) * sqrt(sys%dx) &
					 * SUM( new_fnx(n,:) * exp( sign_k * Ic*sys%w(i)*x_arr(:) ) )
		end do
	 end do

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: f_nx_EO
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Constructs the real-space bosonic spatial distribution `f(x)` associated 
  !>   with the qubit "up" state manifold. It transforms the symmetric (`f`) and 
  !>   antisymmetric (`fo`) momentum arrays into real-space coordinates using the 
  !>   dispersion relation embedded in the grid.
  !> Arguments:
  !>   - sys : Parameter structure defining the spatial grid `dx`.
  !>   - st  : Current state object containing the momentum arrays.
  !> Return:
  !>   - complex(cx) : The real-space field array for the "up" branch.
  !>
  FUNCTION f_nx_EO(sys,st)

	 type(param),intent(in)   	::  sys
	 type(state),intent(in)    ::  st
	 real(rl)				 	   ::  x
	 integer 						::  n,i
	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode)  ::  f_nx_eo

	 do n=1,st%np
		do i=-sys%nmode+1,sys%nmode
		  x = sys%dx * (i-0.5_rl)
		  f_nx_eo(n,i) = sqrt(0.5_rl/(2._rl*pi)) * sqrt(sys%dx) &
					 * SUM( sqrt(sys%dk1) * (st%f(n,:)+st%fo(n,:)) * exp( Ic*sys%w(:)*x ) &
								+ sqrt(sys%dk1) * (st%f(n,:)-st%fo(n,:)) * exp( -Ic*sys%w(:)*x ) )
		end do
	 end do

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: h_nx_EO
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Constructs the real-space bosonic spatial distribution `h(x)` associated 
  !>   with the qubit "down" state manifold. Counterpart to `f_nx_EO`, converting 
  !>   the `h` and `ho` arrays into their coordinate-space representations.
  !> Arguments:
  !>   - sys : Parameter structure defining the spatial grid `dx`.
  !>   - st  : Current state object containing the momentum arrays.
  !> Return:
  !>   - complex(cx) : The real-space field array for the "down" branch.
  !>
  FUNCTION h_nx_EO(sys,st)


	 type(param),intent(in)   	::  sys
	 type(state),intent(in)    ::  st
	 real(rl)				 	   ::  x
	 integer 						::  n,i
	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode)  ::  h_nx_eo


	 do n=1,st%np
		do i=-sys%nmode+1,sys%nmode
		  x = sys%dx*dble(i-0.5_rl)
		  h_nx_eo(n,i) = sqrt(0.5_rl/(2._rl*pi)) *  sqrt(sys%dx)  &
					 * SUM( sqrt(sys%dk1) * (st%h(n,:)+st%ho(n,:)) * exp( Ic*sys%w(:)*x ) &
								+ sqrt(sys%dk1) *(st%h(n,:)-st%ho(n,:)) * exp( -Ic*sys%w(:)*x ) )
		end do
	 end do

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: n_up_k_EO
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the momentum-space mean photon number spectrum strictly projected 
  !>   onto the qubit's excited "up" state subspace. It zeroes out the "down" 
  !>   amplitudes, normalizes the intermediate state, and evaluates the expectation 
  !>   value of the photon number operator a^†_k a_k.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : Array containing the spin-projected photon spectrum.
  !>
  FUNCTION n_up_k_EO(sys,st)

	 type(state), intent(in)   						::  st
	 type(param), intent(in)							::  sys
	 real(rl)												::  n_up_k_EO( -sys%nmode+1:sys%nmode )
	 type(state)											::  upst
	 complex(cx)											::  tmp
	 integer													::  n,m,j,i
	 complex(cx),dimension( st%np,-sys%nmode+1:sys%nmode ) ::  fk

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( upst%f + upst%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( upst%f(:,(abs(j)+1)) - upst%fo(:,(abs(j)+1)) )
	 end do

	 do i=-sys%nmode+1, sys%nmode
		tmp = 0._cx
		do n=1,st%np
		  do m=1,st%np

			 tmp = tmp + conjg(upst%p(n))*upst%p(m)*(conjg(fk(n,i))*fk(m,i))*ov( fk(n,:) , fk(m,:) )

		  end do
		end do
		n_up_k_EO(i) = real( tmp )
	 end do

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: n_down_k_EO
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the momentum-space mean photon number spectrum strictly projected 
  !>   onto the qubit's ground "down" state subspace. It zeroes out the "up" 
  !>   amplitudes, normalizes, and calculates the expected photon count across `k`.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : Array containing the spin-projected photon spectrum.
  !>
  FUNCTION n_down_k_EO(sys,st)

	 type(state), intent(in)   									 ::  st
	 type(param), intent(in)										 ::  sys
	 real(rl)															 ::  n_down_k_EO( -sys%nmode+1:sys%nmode )
	 type(state)														 ::  downst
	 complex(cx)														 ::  tmp
	 integer																 ::  n,m,j,i
	 complex(cx),dimension( st%np,-sys%nmode+1:sys%nmode ) ::  hk

	 downst = st
	 downst%p(:) = 0._cx
	 CALL normalise(downst)

	 hk(:,1:sys%nmode) = sqrt(0.5_rl)*( downst%h + downst%ho )
	 do j = -sys%nmode +1,0
		hk(:,j) = sqrt(0.5_rl)*( downst%h(:,(abs(j)+1)) - downst%ho(:,(abs(j)+1)) )
	 end do

	 do i=-sys%nmode+1,sys%nmode
		tmp = 0._cx
		do n=1,st%np
		  do m=1,st%np

			 tmp = tmp + conjg(downst%q(n))*downst%q(m)*(conjg(hk(n,i))*hk(m,i))*ov( hk(n,:) , hk(m,:) )

		  end do
		end do

		n_down_k_EO(i) = real( tmp )
	 end do

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: n_k_EO
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the total momentum-space mean photon number spectrum, n_k, for the 
  !>   full system. It recombines the symmetric and antisymmetric modes into standard 
  !>   plane-wave bases and calculates the total photon count by summing the 
  !>   contributions from both the "up" and "down" qubit branches.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : Array containing the total photon spectrum.
  !>
  FUNCTION n_k_EO(sys,st)

	 type(state), intent(in)   										::  st
	 type(param), intent(in)											::  sys
	 real(rl)														 		::  n_k_EO( -sys%nmode+1:sys%nmode )
	 complex(cx)													  		::  tmp1, tmp2
	 integer																	::  n,m,j,i
	 complex(cx),dimension( st%np,-sys%nmode+1:sys%nmode )	::  fk, hk

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 hk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%h + st%ho )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
		hk(:,j) = sqrt(0.5_rl)*( st%h(:,(abs(j)+1)) - st%ho(:,(abs(j)+1)) )
	 end do

	 do i=-sys%nmode+1, sys%nmode
		tmp1 = 0._cx
		tmp2 = 0._cx
		do n=1,st%np
		  do m=1,st%np

			 tmp1 = tmp1 + conjg(st%p(n))*st%p(m)*(conjg(fk(n,i))*fk(m,i))*ov( fk(n,:) , fk(m,:) )
			 tmp2 = tmp2 + conjg(st%q(n))*st%q(m)*(conjg(hk(n,i))*hk(m,i))*ov( hk(n,:) , hk(m,:) )

		  end do
		end do

		n_k_EO(i) = real(tmp1)+real(tmp2)
	 end do

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: n_x_EO
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the total spatial photon density profile, n(x). Transforms the 
  !>   Even/Odd momentum arrays into real space via `f_nx_eo` and `h_nx_eo`, 
  !>   and evaluates the expectation value of the local photon number operator 
  !>   across both the "up" and "down" qubit manifolds.
  !> Arguments:
  !>   - sys : Parameter structure defining the spatial grid `dx`.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array containing the total real-space photon density.
  !>
  FUNCTION n_x_EO(sys,st)

	 type(state), intent(in)   										::  st
	 type(param), intent(in)											::  sys
	 real(rl)														 		::  n_x_EO( -sys%nmode+1:sys%nmode )
	 complex(cx)													  		::  tmp
	 integer																	::  n,m,i
	 complex(cx),dimension( st%np,-sys%nmode+1:sys%nmode )	::  fnx, hnx

	 fnx = f_nx_eo(sys,st)
	 hnx = h_nx_eo(sys,st)

	 do i=-sys%nmode+1, sys%nmode
		tmp = 0._cx
		do n=1,st%np
		  do m=1,st%np

			 tmp = tmp + conjg(st%p(n))*st%p(m)*(conjg(fnx(n,i))*fnx(m,i))*ov( fnx(n,:) , fnx(m,:) ) &
					 + conjg(st%q(n))*st%q(m)*(conjg(hnx(n,i))*hnx(m,i))*ov( hnx(n,:) , hnx(m,:) )

		  end do
		end do

		n_x_EO(i) = real(tmp)
	 end do

  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: n_up_eo
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Calculates the total integrated, macroscopic number of photons associated 
  !>   exclusively with the qubit "up" state. Sums the total mode occupation 
  !>   across the entire spatial/momentum grid while projecting out the "down" 
  !>   state contributions.
  !> Arguments:
  !>   - st : Current state object to evaluate.
  !> Return:
  !>   - real(rl) : Total macroscopic photon count in the "up" subspace.
  !>
  FUNCTION n_up_eo(st)

	 type(state),intent(in)			::  st
	 real(rl)							::  n_up_eo
	 type(state)						::  upst
	 complex(cx)						::  sum_k, tmp
	 integer							 	::  n,m

	 upst = st
	 upst%q(:) = 0._cx
	 CALL normalise(upst)

	 tmp = 0._cx
	 do n=1, st%np
		do m=1,st%np

		  sum_k = 0.5_rl*sum( conjg(st%f(n,:)+st%fo(n,:)) * (st%f(m,:)+st%fo(m,:)) &
					 + conjg(st%f(n,:)-st%fo(n,:)) * (st%f(m,:)-st%fo(m,:)) )
		  tmp = tmp + conjg(upst%p(n))*upst%p(m)*sum_k*upst%ov_ff(n,m)

		end do
	 end do
	 
	 n_up_eo = real(tmp)

  END FUNCTION

  !======================================================
  !== PHOTON DECOMPOSITION FUNCITONS
  !======================================================
  
!> -------------------------------------------------------------------------
  !> FUNCTION: one_photon_k_amp_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the probability amplitude for the single-photon momentum state 
  !>   |1_k, up>. Projects the complex multi-polaron ansatz onto the single-photon 
  !>   Fock space for the excited qubit manifold. Essential for evaluating elastic 
  !>   scattering and single-photon transmission/reflection spectra.
  !> Arguments:
  !>   - sys : Parameter structure defining the momentum grid.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : Array containing the complex single-photon amplitudes across k.
  !>
  FUNCTION one_photon_k_amp_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  one_photon_k_amp_up( -sys%nmode+1 : sys%nmode)
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k

    one_photon_k_amp_up = 0._cx
    null_arr = 0._cx

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k = -sys%nmode+1, sys%nmode

		  one_photon_k_amp_up(k) = sum( st%p(:) * fk(:,k) * ov_0f(:) )

	 end do
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: one_photon_k_amp_down
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the probability amplitude for the single-photon momentum state 
  !>   |1_k, down>. Counterpart to `one_photon_k_amp_up`, mapping the single-photon 
  !>   Fock state projection for the qubit ground state manifold.
  !> Arguments:
  !>   - sys : Parameter structure defining the momentum grid.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : Array containing the complex single-photon amplitudes across k.
  !>
  FUNCTION one_photon_k_amp_down(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  one_photon_k_amp_down( -sys%nmode+1 : sys%nmode)
	 complex(cx) 						 ::  hk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0h( st%np )
	 integer								 ::  j,n,k

    one_photon_k_amp_down = 0._cx
    null_arr = 0._cx

	 hk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%h + st%ho )
	 do j = -sys%nmode +1,0
		hk(:,j) = sqrt(0.5_rl)*( st%h(:,(abs(j)+1)) - st%ho(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0h(n) = ov( null_arr , hk(n,:) )
	 end do

	 do k = -sys%nmode+1, sys%nmode

		  one_photon_k_amp_down(k) = sum( st%q(:) * hk(:,k) * ov_0h(:) )

	 end do
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: one_photon_x_amp_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the probability amplitude for the single-photon spatial state 
  !>   |1_x, up>. Transforms the multi-polaron projection into the real-space 
  !>   single-photon basis. Useful for tracking the shape and dispersion of the 
  !>   fundamental scattered wavepacket in time.
  !> Arguments:
  !>   - sys : Parameter structure defining the spatial grid `dx`.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : Array of complex single-photon spatial amplitudes.
  !>
  FUNCTION one_photon_x_amp_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  one_photon_x_amp_up( -sys%nmode+1 : sys%nmode)
	 complex(cx) 						 ::  fx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  n,i

    one_photon_x_amp_up = 0._cx
    null_arr = 0._cx

	 fx = f_nx_eo(sys,st)

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fx(n,:) )
	 end do

	 do i = -sys%nmode+1, sys%nmode

		  one_photon_x_amp_up(i) = sum( st%p(:) * fx(:,i) * ov_0f(:) )

	 end do
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: one_photon_x_amp_down
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the probability amplitude for the single-photon spatial state 
  !>   |1_x, down>. Counterpart to `one_photon_x_amp_up`, evaluating the spatial 
  !>   wavefunction of the single photon conditioned on the qubit being in its ground state.
  !> Arguments:
  !>   - sys : Parameter structure defining the spatial grid `dx`.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : Array of complex single-photon spatial amplitudes.
  !>
  FUNCTION one_photon_x_amp_down(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  one_photon_x_amp_down( -sys%nmode+1 : sys%nmode)
	 complex(cx) 						 ::  hx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0h( st%np )
	 integer								 ::  n,i

    one_photon_x_amp_down = 0._cx
    null_arr = 0._cx

	 hx = h_nx_eo(sys,st)

	 do n=1,st%np
		ov_0h(n) = ov( null_arr , hx(n,:) )
	 end do

	 do i = -sys%nmode+1, sys%nmode

		  one_photon_x_amp_down(i) = sum( st%q(:) * hx(:,i) * ov_0h(:) )

	 end do
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: two_photon_kk_amp_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the joint probability amplitude for the two-photon momentum state 
  !>   |1_k1, 1_k2, up>. Projects the multi-polaron state onto the two-photon 
  !>   Fock basis. Crucial for observing inelastic scattering phenomena where 
  !>   a single high-energy photon splits into two lower-energy photons.
  !> Arguments:
  !>   - sys : Parameter structure defining the momentum grid.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : 2D array matrix containing the complex two-photon amplitudes P(k1, k2).
  !>
  FUNCTION two_photon_kk_amp_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  two_photon_kk_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k1,k2

    two_photon_kk_amp_up = 0._cx
    null_arr = 0._cx


	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k1=-sys%nmode+1,sys%nmode

		do k2=-sys%nmode+1,sys%nmode

		  two_photon_kk_amp_up(k1,k2) = &
				  sum( st%p(:) * fk(:,k1) * fk(:,k2) * ov_0f(:) )/2._rl
		  

		end do

	 end do

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: two_photon_kk_amp_down
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the joint probability amplitude for the two-photon momentum state 
  !>   |1_k1, 1_k2, down>. Evaluates the inelastic frequency conversion spectra 
  !>   associated with the qubit returning to its ground state after scattering.
  !> Arguments:
  !>   - sys : Parameter structure defining the momentum grid.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : 2D array matrix containing the complex two-photon amplitudes P(k1, k2).
  !>
  FUNCTION two_photon_kk_amp_down(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  two_photon_kk_amp_down( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  hk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0h( st%np )
	 integer								 ::  j,n,k1,k2

    two_photon_kk_amp_down = 0._cx
    null_arr = 0._cx

	 hk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%h + st%ho )
	 do j = -sys%nmode +1,0
		hk(:,j) = sqrt(0.5_rl)*( st%h(:,(abs(j)+1)) - st%ho(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0h(n) = ov( null_arr , hk(n,:) )
	 end do

	 do k1=-sys%nmode+1,sys%nmode

		do k2=-sys%nmode+1,sys%nmode

		  two_photon_kk_amp_down(k1,k2) = &
				  sum( st%q(:) * hk(:,k1) * hk(:,k2) * ov_0h(:) )/2._rl
		  

		end do

	 end do

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: two_photon_xx_amp_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the joint probability amplitude for the two-photon spatial state 
  !>   |1_x1, 1_x2, up>. Maps the two-photon correlation into real space to 
  !>   investigate spatial bunching, antibunching, or the formation of bound 
  !>   photon molecules emerging from the ultrastrong scatterer.
  !> Arguments:
  !>   - sys : Parameter structure defining the spatial grid `dx`.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : 2D array matrix containing spatial two-photon amplitudes P(x1, x2).
  !>
  FUNCTION two_photon_xx_amp_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  two_photon_xx_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  n,i1,i2

    two_photon_xx_amp_up = 0._cx
    null_arr = 0._cx

	 fx = f_nx_eo(sys,st)

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fx(n,:) )
	 end do

	 do i1=-sys%nmode+1,sys%nmode

		do i2=-sys%nmode+1,sys%nmode

		  two_photon_xx_amp_up(i1,i2) = &
				  sum( st%p(:) * fx(:,i1) * fx(:,i2) * ov_0f(:) )/2._rl
		  

		end do

	 end do

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: two_photon_xx_amp_p
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Similar physical context to `two_photon_xx_amp_up` (spatial two-photon 
  !>   bunching amplitudes). This formulation isolates the evaluation explicitly 
  !>   to the coherent amplitude weights `p` without decomposing the full spin 
  !>   subspace, acting as a lightweight diagnostic for up-manifold correlations.
  !> Arguments:
  !>   - sys : Parameter structure defining the spatial grid `dx`.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : 2D array matrix of spatial two-photon amplitudes.
  !>
  FUNCTION two_photon_xx_amp_p(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  two_photon_xx_amp_p( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  n,i1,i2

    two_photon_xx_amp_p = 0._cx
    null_arr = 0._cx

	 fx = f_nx_eo(sys,st)

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fx(n,:) )
	 end do

	 do i1=-sys%nmode+1,sys%nmode

		do i2=-sys%nmode+1,sys%nmode

		  two_photon_xx_amp_p(i1,i2) = sum( st%p(:) * fx(:,i1) * fx(:,i2) * ov_0f(:) )/2._rl !+ sum( st%q(:) * hx(:,i1) * hx(:,i2) * ov_0h(:) )/2._rl

		end do

	 end do

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: two_photon_xx_amp_q
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the joint probability amplitude for the two-photon spatial state 
  !>   |1_x1, 1_x2, down>. Acts as the counterpart to `two_photon_xx_amp_p`, 
  !>   evaluating real-space two-photon clustering (bunching/antibunching) 
  !>   specifically for the down-manifold using the `q` amplitudes and `h` modes.
  !> Arguments:
  !>   - sys : Parameter structure defining the spatial grid `dx`.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : 2D array matrix of spatial two-photon amplitudes.
  !>
  FUNCTION two_photon_xx_amp_q(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  two_photon_xx_amp_q( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  hx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0h( st%np )
	 integer								 ::  n,i1,i2

    two_photon_xx_amp_q = 0._cx
    null_arr = 0._cx

	 hx = h_nx_eo(sys,st)

	 do n=1,st%np
		ov_0h(n) = ov( null_arr , hx(n,:) )
	 end do

	 do i1=-sys%nmode+1,sys%nmode

		do i2=-sys%nmode+1,sys%nmode

		  two_photon_xx_amp_q(i1,i2) = sum( st%q(:) * hx(:,i1) * hx(:,i2) * ov_0h(:) )/2._rl !+ sum( st%q(:) * hx(:,i1) * hx(:,i2) * ov_0h(:) )/2._rl

		end do

	 end do

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: three_photon_kk_amp_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the joint probability amplitude for a three-photon momentum state 
  !>   |1_k1, 1_k2, 1_k_in, up>. By fixing one of the photon momenta (`k_in`), 
  !>   it collapses the 3D tensor into a 2D matrix, allowing visualization of 
  !>   three-body scattering channels where the artificial atom is left excited.
  !> Arguments:
  !>   - sys  : Parameter structure.
  !>   - st   : Current state object.
  !>   - k_in : The specific momentum index fixed for the third photon.
  !> Return:
  !>   - complex(cx) : 2D array matrix of the constrained three-photon amplitudes.
  !>
  FUNCTION three_photon_kk_amp_up(sys,st,k_in)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 integer, intent(in)		 ::  k_in
	 complex(cx)					  	 ::  three_photon_kk_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k1,k2

    three_photon_kk_amp_up = 0._cx
    null_arr = 0._cx


	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k1=-sys%nmode+1,sys%nmode

		do k2=-sys%nmode+1,sys%nmode

		  three_photon_kk_amp_up(k1,k2) = &
				  sum( st%p(:) * fk(:,k1) * fk(:,k2) * fk(:,k_in) * ov_0f(:) )/6._rl

		end do

	 end do 
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: three_photon_xx_amp_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the joint probability amplitude for a three-photon spatial state 
  !>   |1_x1, 1_x2, 1_x_in, up>. Fixes one spatial coordinate (`i_in`) to 
  !>   evaluate the real-space two-body correlations of the remaining photons 
  !>   relative to the fixed third photon.
  !> Arguments:
  !>   - sys  : Parameter structure.
  !>   - st   : Current state object.
  !>   - i_in : The specific spatial grid index fixed for the third photon.
  !> Return:
  !>   - complex(cx) : 2D array matrix of constrained spatial amplitudes.
  !>
  FUNCTION three_photon_xx_amp_up(sys,st,i_in)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 integer, intent(in)		 		 ::  i_in
	 complex(cx)					  	 ::  three_photon_xx_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  n,i1,i2

    three_photon_xx_amp_up = 0._cx
    null_arr = 0._cx

	 fx = f_nx_eo(sys,st)

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fx(n,:) )
	 end do

	 do i1=-sys%nmode+1,sys%nmode

		do i2=-sys%nmode+1,sys%nmode

		  three_photon_xx_amp_up(i1,i2) = &
				  sum( st%p(:) * fx(:,i1) * fx(:,i2) * fx(:,i_in) * ov_0f(:) )/6._rl

		end do

	 end do

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: three_photon_xxx_amp_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Evaluates the full 3D spatial probability amplitude for the three-photon 
  !>   state |1_x1, 1_x2, 1_x3, up>. Resolves the complete spatial correlation 
  !>   tensor without constraints, incorporating the bosonic permutation 
  !>   combinatoric factor (1/6 = 1/3!).
  !> Arguments:
  !>   - sys : Parameter structure defining the spatial bounds.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : 3D array tensor containing the full spatial amplitudes.
  !>
  FUNCTION three_photon_xxx_amp_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  three_photon_xxx_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  n,i1,i2,i3

    three_photon_xxx_amp_up = 0._cx
    null_arr = 0._cx

	 fx = f_nx_eo(sys,st)

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fx(n,:) )
	 end do

	 do i1=-sys%nmode+1,sys%nmode

		do i2=-sys%nmode+1,sys%nmode

		  do i3=-sys%nmode+1,sys%nmode

			 three_photon_xxx_amp_up(i1,i2,i3) = &
				  sum( st%p(:) * fx(:,i1) * fx(:,i2) * fx(:,i3) * ov_0f(:) )/6._rl

		  end do

		end do

	 end do

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: three_photon_xxx_amp_p
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Evaluates the 3D spatial three-photon correlation strictly using the 
  !>   coherent weight amplitudes `p` and spatial fields `fx`. Acts as a 
  !>   lightweight alias/variant for the excited manifold projection, identical 
  !>   in physics context to `three_photon_xxx_amp_up`.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : 3D array tensor of spatial amplitudes.
  !>
  FUNCTION three_photon_xxx_amp_p(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  three_photon_xxx_amp_p( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  n,i1,i2,i3

    three_photon_xxx_amp_p = 0._cx
    null_arr = 0._cx

	 fx = f_nx_eo(sys,st)

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fx(n,:) )
	 end do

	 do i1=-sys%nmode+1,sys%nmode

		do i2=-sys%nmode+1,sys%nmode

		  do i3=-sys%nmode+1,sys%nmode

			 three_photon_xxx_amp_p(i1,i2,i3) = sum( st%p(:) * fx(:,i1) * fx(:,i2) * fx(:,i3) * ov_0f(:) )/6._rl

		  end do

		end do

	 end do

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> SUBROUTINE: sub_three_photon_xxx_amp_p
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   A memory-optimized block evaluation of the 3D three-photon spatial 
  !>   probability. To manage stack limits and optimize cache performance 
  !>   during the heavy O(N^3) tensor population, it breaks the coordinate 
  !>   space into distinct contiguous blocks (`ar1` and `ar2`) and populates 
  !>   them sequentially.
  !> Arguments:
  !>   - sys      : Parameter structure.
  !>   - st       : Current state object.
  !>   - ar1, ar2 : Target 3D output arrays representing spatial blocks.
  !>
  SUBROUTINE sub_three_photon_xxx_amp_p(sys,st,ar1,ar2)

	 type(param), intent(in)    	 		::  sys
	 type(state), intent(in)		 		::  st
	 complex(cx)								::  ar1( -sys%nmode+1:0 , -sys%nmode+1:0, -sys%nmode+1:0 )
	 complex(cx)								::  ar2( -sys%nmode+1:0 , -sys%nmode+1:0, 1:sys%nmode )
	 complex(cx) 						 		::  fx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 		::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 		::  ov_0f( st%np )
	 integer								 		::  n,i1,i2,i3

    null_arr = 0._cx

	 fx = f_nx_eo(sys,st)

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fx(n,:) )
	 end do


	 print*,"start"
	 do i1=-sys%nmode+1, 0
		do i2=-sys%nmode+1, 0
		  do i3=-sys%nmode+1, 0
			 ar1(i1,i2,i3) = sum( st%p(:) * fx(:,i1) * fx(:,i2) * fx(:,i3) * ov_0f(:) )/6._rl
		  end do
		end do
	 end do
	 print*,"1/2 done"
	 do i1=-sys%nmode+1, 0
		do i2=-sys%nmode+1, 0
		  do i3=1, sys%nmode
			 ar2(i1,i2,i3) = sum( st%p(:) * fx(:,i1) * fx(:,i2) * fx(:,i3) * ov_0f(:) )/6._rl
		  end do
		end do
	 end do
	 print*,"done"

	 
  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> SUBROUTINE: sub_three_photon_xxx_amp_p
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   A memory-optimized block evaluation of the 3D three-photon spatial 
  !>   probability. To manage stack limits and optimize cache performance 
  !>   during the heavy O(N^3) tensor population, it breaks the coordinate 
  !>   space into distinct contiguous blocks (`ar1` and `ar2`) and populates 
  !>   them sequentially.
  !> Arguments:
  !>   - sys      : Parameter structure.
  !>   - st       : Current state object.
  !>   - ar1, ar2 : Target 3D output arrays representing spatial blocks.
  !>
  SUBROUTINE sub_three_photon_xxx_amp_q(sys,st,ar1,ar2)

	 type(param), intent(in)    	 		::  sys
	 type(state), intent(in)		 		::  st
	 complex(cx)								::  ar1( -sys%nmode+1:0 , -sys%nmode+1:0, -sys%nmode+1:0 )
	 complex(cx)								::  ar2( -sys%nmode+1:0 , -sys%nmode+1:0, 1:sys%nmode )
	 complex(cx) 						 		::  hx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 		::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 		::  ov_0h( st%np )
	 integer								 		::  n,i1,i2,i3

    null_arr = 0._cx

	 hx = h_nx_eo(sys,st)

	 do n=1,st%np
		ov_0h(n) = ov( null_arr , hx(n,:) )
	 end do


	 do i1=-sys%nmode+1, 0
		do i2=-sys%nmode+1, 0
		  do i3=-sys%nmode+1, 0
			 ar1(i1,i2,i3) = sum( st%q(:) * hx(:,i1) * hx(:,i2) * hx(:,i3) * ov_0h(:) )/6._rl
		  end do
		end do
	 end do
	 do i1=-sys%nmode+1, 0
		do i2=-sys%nmode+1, 0
		  do i3=1, sys%nmode
			 ar2(i1,i2,i3) = sum( st%q(:) * hx(:,i1) * hx(:,i2) * hx(:,i3) * ov_0h(:) )/6._rl
		  end do
		end do
	 end do
	 
  END SUBROUTINE

!> -------------------------------------------------------------------------
  !> FUNCTION: three_photon_xxx_amp_q
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Evaluates the full 3D spatial probability amplitude for the three-photon 
  !>   state |1_x1, 1_x2, 1_x3, down>. To dramatically reduce O(N^3) runtime, 
  !>   it exploits bosonic permutation symmetry by calculating only unique index 
  !>   combinations (i1 <= i2 <= i3) and mirroring the symmetric results.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - complex(cx) : 3D array tensor containing the populated spatial amplitudes.
  !>
  FUNCTION three_photon_xxx_amp_q(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 complex(cx)					  	 ::  three_photon_xxx_amp_q( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  hx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0h( st%np )
	 integer								 ::  n,i1,i2,i3

    three_photon_xxx_amp_q = 0._cx
    null_arr = 0._cx

	 hx = h_nx_eo(sys,st)

	 do n=1,st%np
		ov_0h(n) = ov( null_arr , hx(n,:) )
	 end do

	 do i1=-sys%nmode+1,sys%nmode

		do i2=i1,sys%nmode

		  do i3=i2,sys%nmode

			 three_photon_xxx_amp_q(i1,i2,i3) = sum( st%q(:) * hx(:,i1) * hx(:,i2) * hx(:,i3) * ov_0h(:) )/6._rl
			 three_photon_xxx_amp_q(i1,i3,i2) = three_photon_xxx_amp_q(i1,i2,i3)

		  end do
		  three_photon_xxx_amp_q(i2,i1,:) = three_photon_xxx_amp_q(i1,i2,:)

		end do

	 end do

	 
  END FUNCTION
  
!> -------------------------------------------------------------------------
  !> FUNCTION: four_photon_kk_amp_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Extracts the joint probability amplitude for the four-photon momentum 
  !>   state |1_k1, 1_k2, 1_k_in_1, 1_k_in_2, up>. Fixes two of the photon 
  !>   momenta to render a 2D slice of the full 4D parameter space. Incorporates 
  !>   the permutation combinatoric factor (1/24 = 1/4!).
  !> Arguments:
  !>   - sys            : Parameter structure.
  !>   - st             : Current state object.
  !>   - k_in_1, k_in_2 : Fixed momentum grid indices for the third and fourth photons.
  !> Return:
  !>   - complex(cx) : 2D array matrix of the constrained four-photon amplitudes.
  !>
  FUNCTION four_photon_kk_amp_up(sys,st,k_in_1,k_in_2)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 integer, intent(in)		 		 ::  k_in_1, k_in_2
	 complex(cx)					  	 ::  four_photon_kk_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k1,k2

    four_photon_kk_amp_up = 0._cx
    null_arr = 0._cx

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k1=-sys%nmode+1,sys%nmode

		do k2=-sys%nmode+1,sys%nmode

		  four_photon_kk_amp_up(k1,k2) = &
				  sum( st%p(:) * fk(:,k1) * fk(:,k2) * fk(:,k_in_1) * fk(:,k_in_2) * ov_0f(:) )/24._rl

		end do

	 end do

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nk_2_photon
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the integrated momentum-space two-photon spectrum. Evaluates 
  !>   the normalized squared amplitudes of the multi-polaron state projected 
  !>   onto the two-photon Fock basis. It traces over both the up and down 
  !>   qubit manifolds to extract the total two-photon cross-section.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array representing the marginal two-photon momentum density.
  !>
  FUNCTION nk_2_photon(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 real(rl)							 ::  nk_2_photon( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  two_photon_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  two_photon_amp_down( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 real(rl)					  	 	 ::  two_photon_amp_norm2( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  hk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np ),ov_0h( st%np )
	 integer								 ::  j,n,k1,k2

	 nk_2_photon = 0._rl
    two_photon_amp_up = 0._cx
    two_photon_amp_down = 0._cx
    two_photon_amp_norm2 = 0._rl
    null_arr = 0._cx

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 hk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%h + st%ho )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
		hk(:,j) = sqrt(0.5_rl)*( st%h(:,(abs(j)+1)) - st%ho(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
		ov_0h(n) = ov( null_arr , hk(n,:) )
	 end do

	 do k1=-sys%nmode+1, sys%nmode

		do k2=-sys%nmode+1, sys%nmode

			 two_photon_amp_up(k1,k2) = &
				  sum( st%p(:) * fk(:,k1) * fk(:,k2) * ov_0f(:) ) / 2._rl
			 two_photon_amp_down(k1,k2) = &
				  sum( st%q(:) * hk(:,k1) * hk(:,k2) * ov_0h(:) ) / 2._rl

		end do

		two_photon_amp_norm2 = real( conjg(two_photon_amp_up)*two_photon_amp_up + conjg(two_photon_amp_down)*two_photon_amp_down )

		do k2=-sys%nmode+1, sys%nmode

			 nk_2_photon( k1 ) = nk_2_photon( k1 ) + two_photon_amp_norm2(k1,k2)

		end do

	 end do
		
	 nk_2_photon = 4 * nk_2_photon/sys%dk1
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nk_3_photon
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the integrated momentum-space three-photon spectrum. Reconstructs 
  !>   the density of three-photon emission events by squaring the three-photon 
  !>   amplitudes. It loops over a reduced index space and heavily relies on 
  !>   integration symmetry factors (1, 2) to correctly account for indistinguishable 
  !>   bosonic configurations while optimizing O(N^3) runtimes.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array representing the marginal three-photon momentum density.
  !>
  FUNCTION nk_3_photon(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 real(rl)							 ::  nk_3_photon( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  three_photon_amp_up( -int(sys%nmode/12.0)+1:int(sys%nmode/12.0) &
														, -int(sys%nmode/12.0)+1:int(sys%nmode/12.0) &
														, -int(sys%nmode/12.0)+1:int(sys%nmode/12.0) )
	 complex(cx)					  	 ::  three_photon_amp_down( -int(sys%nmode/12.0)+1:int(sys%nmode/12.0) &
														, -int(sys%nmode/12.0)+1:int(sys%nmode/12.0) &
														, -int(sys%nmode/12.0)+1:int(sys%nmode/12.0) )
	 real(rl)						  	 ::  three_photon_amp_norm2( -int(sys%nmode/12.0)+1:int(sys%nmode/12.0) &
														, -int(sys%nmode/12.0)+1:int(sys%nmode/12.0) &
														, -int(sys%nmode/12.0)+1:int(sys%nmode/12.0) )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  hk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np ), ov_0h( st%np )
	 integer								 ::  j,n,k1,k2,k3,factor

	 nk_3_photon = 0._rl
    three_photon_amp_up = 0._cx
    three_photon_amp_down = 0._cx
    three_photon_amp_norm2 = 0._rl
    null_arr = 0._cx
	 factor=2

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 hk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%h + st%ho )
	 do j = -sys%nmode +1,0
		hk(:,j) = sqrt(0.5_rl)*( st%h(:,(abs(j)+1)) - st%ho(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
		ov_0h(n) = ov( null_arr , hk(n,:) )
	 end do

	 do k1=-int(sys%nmode/12.0)+1,int(sys%nmode/12.0)

		do k2=-int(sys%nmode/12.0)+1,int(sys%nmode/12.0)

		  do k3=k2,int(sys%nmode/12.0)

			 three_photon_amp_up(k1,k2,k3) = &
				  sum( st%p(:) * fk(:,k1) * fk(:,k2) * fk(:,k3) * ov_0f(:) ) / 6._rl
			 three_photon_amp_down(k1,k2,k3) = &
				  sum( st%q(:) * hk(:,k1) * hk(:,k2) * hk(:,k3) * ov_0h(:) ) / 6._rl

		  end do

		end do

		three_photon_amp_norm2 = real( conjg(three_photon_amp_up)*three_photon_amp_up &
														+ conjg(three_photon_amp_down)*three_photon_amp_down )

		do k2=-int(sys%nmode/12.0)+1, int(sys%nmode/12.0)

		  do k3= k2, int(sys%nmode/12.0)

			 if ( k2 == k3 ) then
				factor = 1
			 elseif ( k2 .ne. k3 ) then
				factor = 2
			 end if

			 nk_3_photon( k1 ) = nk_3_photon( k1 ) + factor*three_photon_amp_norm2(k1,k2,k3)

		  end do

		end do

	 end do
		
	 nk_3_photon = 18*nk_3_photon/sys%dk1
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nk_2_photon_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the integrated momentum-space two-photon spectrum strictly 
  !>   projected onto the qubit's excited "up" manifold. Useful for correlating 
  !>   multiphoton scattering events specifically with the atomic state.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array representing the up-projected two-photon spectrum.
  !>
  FUNCTION nk_2_photon_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 real(rl)							 ::  nk_2_photon_up( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  two_photon_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 real(rl)					  	 	 ::  two_photon_amp_up_norm2( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k1,k2

	 nk_2_photon_up = 0._rl
    two_photon_amp_up = 0._cx
    two_photon_amp_up_norm2 = 0._rl
    null_arr = 0._cx

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k1=-sys%nmode+1, sys%nmode

		do k2=-sys%nmode+1, sys%nmode

			 two_photon_amp_up(k1,k2) = &
				  sum( st%p(:) * fk(:,k1) * fk(:,k2) * ov_0f(:) ) / 2._rl

		end do

		two_photon_amp_up_norm2 = real( conjg(two_photon_amp_up)*two_photon_amp_up )

		do k2=-sys%nmode+1, sys%nmode

			 nk_2_photon_up( k1 ) = nk_2_photon_up( k1 ) + two_photon_amp_up_norm2(k1,k2)

		end do

	 end do
		
	 nk_2_photon_up = 4 * nk_2_photon_up/sys%dk1
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nk_2_photon_up_RR_LL
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the two-photon momentum spectrum constrained to strictly 
  !>   co-propagating photon pairs (both Right-moving or both Left-moving) 
  !>   within the "up" manifold. The inner loops conditionally filter `k1` and 
  !>   `k2` indices to enforce same-sign momentum matching.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array representing the co-propagating two-photon density.
  !>
  FUNCTION nk_2_photon_up_RR_LL(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 real(rl)							 ::  nk_2_photon_up_RR_LL( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  two_photon_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 real(rl)					  	 	 ::  two_photon_amp_up_norm2( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k1,k2

	 nk_2_photon_up_RR_LL = 0._rl
    two_photon_amp_up = 0._cx
    two_photon_amp_up_norm2 = 0._rl
    null_arr = 0._cx

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k1=-sys%nmode+1, sys%nmode

		do k2=-sys%nmode+1, sys%nmode

			 two_photon_amp_up(k1,k2) = &
				  sum( st%p(:) * fk(:,k1) * fk(:,k2) * ov_0f(:) ) / 2._rl

		end do

		two_photon_amp_up_norm2 = real( conjg(two_photon_amp_up)*two_photon_amp_up )

		do k2=-sys%nmode+1, sys%nmode

		  if ( ( (k2>=0) .and. (k1>=0) ) .or. ( (k2<0) .and. (k1<0) ) ) then
			 nk_2_photon_up_RR_LL( k1 ) = nk_2_photon_up_RR_LL( k1 ) + two_photon_amp_up_norm2(k1,k2)
		  end if

		end do

	 end do
		
	 nk_2_photon_up_RR_LL = 4 * nk_2_photon_up_RR_LL/sys%dk1
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nk_2_photon_up_RL_RL
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the two-photon momentum spectrum constrained to counter-propagating 
  !>   photon pairs (one Right-moving, one Left-moving) within the "up" manifold. 
  !>   Filters indices to enforce opposite-sign momenta, characterizing bidirectional 
  !>   multiphoton emission events.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array representing the counter-propagating two-photon density.
  !>
  FUNCTION nk_2_photon_up_RL_RL(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 real(rl)							 ::  nk_2_photon_up_RL_RL( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  two_photon_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 real(rl)					  	 	 ::  two_photon_amp_up_norm2( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k1,k2

	 nk_2_photon_up_RL_RL = 0._rl
    two_photon_amp_up = 0._cx
    two_photon_amp_up_norm2 = 0._rl
    null_arr = 0._cx

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k1=-sys%nmode+1, sys%nmode

		do k2=-sys%nmode+1, sys%nmode

			 two_photon_amp_up(k1,k2) = &
				  sum( st%p(:) * fk(:,k1) * fk(:,k2) * ov_0f(:) ) / 2._rl

		end do

		two_photon_amp_up_norm2 = real( conjg(two_photon_amp_up)*two_photon_amp_up )

		do k2=-sys%nmode+1, sys%nmode

		  if ( ( (k2<=0) .and. (k1>=0) ) .or. ( (k2>0) .and. (k1<0) ) ) then
			 nk_2_photon_up_RL_RL( k1 ) = nk_2_photon_up_RL_RL( k1 ) + two_photon_amp_up_norm2(k1,k2)
		  end if

		end do

	 end do
		
	 nk_2_photon_up_RL_RL = 4 * nk_2_photon_up_RL_RL/sys%dk1
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nx_2_photon_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the spatial two-photon density profile for the "up" manifold. 
  !>   It projects the 2D spatial two-photon amplitudes and sums the squared 
  !>   modulus over the second coordinate, yielding the marginal probability 
  !>   distribution of finding at least one photon at position x1 given a 
  !>   two-photon event.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array representing the spatial density of two-photon events.
  !>
  FUNCTION nx_2_photon_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 real(rl)							 ::  nx_2_photon_up( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  two_photon_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 real(rl)					  	 	 ::  two_photon_amp_up_norm2( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  fx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  n,i1,i2

	 nx_2_photon_up = 0._rl
    two_photon_amp_up = 0._cx
    two_photon_amp_up_norm2 = 0._rl
    null_arr = 0._cx

	 fx = f_nx_eo(sys,st)

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fx(n,:) )
	 end do

	 do i1=-sys%nmode+1, sys%nmode

		do i2=-sys%nmode+1, sys%nmode

			 two_photon_amp_up(i1,i2) = &
				  sum( st%p(:) * fx(:,i1) * fx(:,i2) * ov_0f(:) ) / 2._rl

		end do

		two_photon_amp_up_norm2 = real( conjg(two_photon_amp_up)*two_photon_amp_up )

		do i2=-sys%nmode+1, sys%nmode

			 nx_2_photon_up( i1 ) = nx_2_photon_up( i1 ) + two_photon_amp_up_norm2(i1,i2)

		end do

	 end do
		
	 nx_2_photon_up = 4 * nx_2_photon_up/sys%dx
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nxx_2_photon_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the equal-coordinate spatial two-photon density profile 
  !>   (x1 = x2) for the "up" manifold. This directly evaluates the probability 
  !>   of finding two photons at the exact same spatial location, which serves 
  !>   as the unnormalized basis for the zero-delay bunching statistic, g2(0).
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array mapping the local two-photon overlap probability.
  !>
  FUNCTION nxx_2_photon_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 real(rl)							 ::  nxx_2_photon_up( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  two_photon_amp_up( -sys%nmode+1:sys%nmode , -sys%nmode+1:sys%nmode )
	 real(rl)					  	 	 ::  two_photon_amp_up_norm2
	 complex(cx) 						 ::  fx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  n,i1

	 nxx_2_photon_up = 0._rl
    two_photon_amp_up = 0._cx
    two_photon_amp_up_norm2 = 0._rl
    null_arr = 0._cx

	 fx = f_nx_eo(sys,st)

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fx(n,:) )
	 end do

	 do i1=-sys%nmode+1, sys%nmode


		two_photon_amp_up(i1,i1) = &
				  sum( st%p(:) * fx(:,i1) * fx(:,i1) * ov_0f(:) ) / 2._rl

		two_photon_amp_up_norm2 = real( conjg(two_photon_amp_up(i1,i1))*two_photon_amp_up(i1,i1) )

		nxx_2_photon_up( i1 ) = two_photon_amp_up_norm2

	 end do
		
	 nxx_2_photon_up = 4 * nxx_2_photon_up/sys%dx
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nk_3_photon_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the integrated momentum-space three-photon spectrum isolated to 
  !>   the qubit "up" state. Sums the squared three-photon projection amplitudes 
  !>   across the restricted frequency grid, carefully applying indistinguishability 
  !>   factors to preserve normalization without redundant coordinate swapping.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array representing the up-projected three-photon spectrum.
  !>
  FUNCTION nk_3_photon_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 real(rl)							 ::  nk_3_photon_up( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  three_photon_amp_up( -int(sys%nmode/6.0)+1:int(sys%nmode/6.0) &
														, -int(sys%nmode/6.0)+1:int(sys%nmode/6.0) &
														, -int(sys%nmode/6.0)+1:int(sys%nmode/6.0) )
	 real(rl)						  	 ::  three_photon_amp_up_norm2( -int(sys%nmode/6.0)+1:int(sys%nmode/6.0) &
														, -int(sys%nmode/6.0)+1:int(sys%nmode/6.0) &
														, -int(sys%nmode/6.0)+1:int(sys%nmode/6.0) )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k1,k2,k3,factor

	 nk_3_photon_up = 0._rl
    three_photon_amp_up = 0._cx
    three_photon_amp_up_norm2 = 0._rl
    null_arr = 0._cx
	 factor=2

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k1=-int(sys%nmode/6.0)+1,int(sys%nmode/6.0)

		do k2=-int(sys%nmode/6.0)+1,int(sys%nmode/6.0)

		  do k3=k2,int(sys%nmode/6.0)

			 three_photon_amp_up(k1,k2,k3) = &
				  sum( st%p(:) * fk(:,k1) * fk(:,k2) * fk(:,k3) * ov_0f(:) ) / 6._rl

		  end do

		end do

		three_photon_amp_up_norm2 = real( conjg(three_photon_amp_up)*three_photon_amp_up )

		do k2=-int(sys%nmode/6.0)+1, int(sys%nmode/6.0)

		  do k3= k2, int(sys%nmode/6.0)

			 if ( k2 == k3 ) then
				factor = 1
			 elseif ( k2 .ne. k3 ) then
				factor = 2
			 end if

			 nk_3_photon_up( k1 ) = nk_3_photon_up( k1 ) + factor*three_photon_amp_up_norm2(k1,k2,k3)

		  end do

		end do

	 end do
		
	 nk_3_photon_up = 18*nk_3_photon_up/sys%dk1
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nx_3_photon_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the spatial three-photon density profile for the "up" manifold. 
  !>   Evaluates the marginal spatial probability distribution for finding a 
  !>   photon at x1, conditioned on the combined spatial spread of the remaining 
  !>   two photons integrated across the entire 3D grid.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array mapping the spatial bounds of three-photon packets.
  !>
  FUNCTION nx_3_photon_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 real(rl)							 ::  nx_3_photon_up( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  three_photon_amp_up( -int(sys%nmode)+1:int(sys%nmode) &
														, -int(sys%nmode)+1:int(sys%nmode) &
														, -int(sys%nmode)+1:int(sys%nmode) )
	 real(rl)						  	 ::  three_photon_amp_up_norm2( -int(sys%nmode)+1:int(sys%nmode) &
														, -int(sys%nmode)+1:int(sys%nmode) &
														, -int(sys%nmode)+1:int(sys%nmode) )
	 complex(cx) 						 ::  fx( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  n,x1,x2,x3,factor

	 nx_3_photon_up = 0._rl
    three_photon_amp_up = 0._cx
    three_photon_amp_up_norm2 = 0._rl
    null_arr = 0._cx
	 factor=2

	 fx = f_nx_eo(sys,st)

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fx(n,:) )
	 end do

	 do x1=-int(sys%nmode)+1,int(sys%nmode)

		do x2=-int(sys%nmode)+1,int(sys%nmode)

		  do x3=x2,int(sys%nmode)

			 three_photon_amp_up(x1,x2,x3) = &
				  sum( st%p(:) * fx(:,x1) * fx(:,x2) * fx(:,x3) * ov_0f(:) ) / 6._rl

		  end do

		end do

		three_photon_amp_up_norm2 = real( conjg(three_photon_amp_up)*three_photon_amp_up )

		do x2=-int(sys%nmode)+1, int(sys%nmode)

		  do x3= x2, int(sys%nmode)

			 if ( x2 == x3 ) then
				factor = 1
			 elseif ( x2 .ne. x3 ) then
				factor = 2
			 end if

			 nx_3_photon_up( x1 ) = nx_3_photon_up( x1 ) + factor*three_photon_amp_up_norm2(x1,x2,x3)

		  end do

		end do

	 end do
		
	 nx_3_photon_up = 18*nx_3_photon_up/sys%dx
	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nk_4_photon_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes the integrated momentum-space four-photon spectrum for the "up" 
  !>   manifold. Evaluates the absolute squares of the complex 4D projection 
  !>   amplitudes. Relies on advanced nested permutation symmetry conditions 
  !>   (factors of 1, 3, 6) to make the massive O(N^4) integral computationally tractable.
  !> Arguments:
  !>   - sys : Parameter structure.
  !>   - st  : Current state object.
  !> Return:
  !>   - real(rl) : Array representing the up-projected four-photon momentum density.
  !>
  FUNCTION nk_4_photon_up(sys,st)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 real(rl)							 ::  nk_4_photon_up( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  four_photon_amp_up( -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) &
														, -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) &
														, -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) )
	 real(rl)						  	 ::  four_photon_amp_up_norm2( -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) &
														, -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) &
														, -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k1,k2,k3,k4,factor

	 nk_4_photon_up = 0._rl
    four_photon_amp_up = 0._cx
    four_photon_amp_up_norm2 = 0._rl
    null_arr = 0._cx
	 factor=2

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k1=-int(sys%nmode/10.0)+1,int(sys%nmode/10.0)


		do k2=-int(sys%nmode/10.0)+1,int(sys%nmode/10.0)

		  do k3=k2,int(sys%nmode/10.0)

			 do k4=k3,int(sys%nmode/10.0)

				four_photon_amp_up(k2,k3,k4) = &
				  sum( st%p(:) * fk(:,k1) * fk(:,k2) * fk(:,k3) * fk(:,k4) * ov_0f(:) ) / 24._rl

			 end do

		  end do

		end do

		four_photon_amp_up_norm2 = real( four_photon_amp_up*conjg(four_photon_amp_up) )

		do k2=-int(sys%nmode/10.0)+1,int(sys%nmode/10.0)

		  do k3=k2,int(sys%nmode/10.0)

			 do k4=k3,int(sys%nmode/10.0)

	   		if ( ( k2 .ne. k3 ) .and. (k2 .ne. k4) .and. ( k3 .ne. k4) ) then
	   			 factor = 6
	   		else if ( ( ( k2 == k3 ) .and. (k2 .ne. k4) .and. ( k3 .ne. k4) ) &
	   		         .or. ( ( k2 .ne. k3 ) .and. (k2 == k4) .and. ( k3 .ne. k4) ) &
	   		         .or. ( ( k2 .ne. k3 ) .and. (k2 .ne. k4) .and. ( k3 == k4) ) ) then
	   			 factor = 3
	   		else if ( ( k2 == k3 ) .and. (k2 == k4) .and. ( k3 == k4) )then
	   		  	 factor = 1
				end if

				nk_4_photon_up( k1 ) = nk_4_photon_up( k1 ) + factor*four_photon_amp_up_norm2(k2,k3,k4)

			 end do

		  end do

		end do
		
	 end do

	 nk_4_photon_up = 92*nk_4_photon_up/sys%dk1

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nk_4_photon_up_k1_k2
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes a selectively constrained four-photon momentum spectrum where 
  !>   two outgoing photon momenta are fixed (`k_in_1` and `k_in_2`). Used 
  !>   primarily to isolate and analyze specific four-body correlated scattering 
  !>   channels extracted from the overarching multi-polaron state. 
  !>   (Note: Current implementation enforces `k_in_1 == k_in_2`).
  !> Arguments:
  !>   - sys            : Parameter structure.
  !>   - st             : Current state object.
  !>   - k_in_1, k_in_2 : Fixed momenta for the reference correlated pair.
  !> Return:
  !>   - real(rl) : Array representing the restricted four-photon spectrum.
  !>
  FUNCTION nk_4_photon_up_k1_k2(sys,st,k_in_1,k_in_2) 
  	  
	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 integer, intent(in)		 		 ::  k_in_1, k_in_2
	 real(rl)							 ::  nk_4_photon_up_k1_k2( -sys%nmode+1:sys%nmode )
	 complex(cx)					  	 ::  four_photon_amp_up( -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) )
	 real(rl)						  	 ::  four_photon_amp_up_norm2( -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k1,k2,factor

	 !--  FOR NOW FUNCTION ONLY VALID IF K_IN_1 = K_IN_2
	 if ( k_in_1 .ne. k_in_2 ) then
		print*, "update function nk: nk_4_photon_up_k1_k2"
		stop
	 end if

	 nk_4_photon_up_k1_k2 = 0._rl
    four_photon_amp_up = 0._cx
    four_photon_amp_up_norm2 = 0._rl
    null_arr = 0._cx
	 factor=3

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k1=-int(sys%nmode/10.0)+1,int(sys%nmode/10.0)

		do k2=k1, int(sys%nmode/10.0)

		  four_photon_amp_up(k2) = &
			 sum( st%p(:) * fk(:,k1) * fk(:,k2) * fk(:,k_in_1) * fk(:,k_in_2) * ov_0f(:) ) / 24._rl

		end do

		four_photon_amp_up_norm2 = real( four_photon_amp_up*conjg(four_photon_amp_up) )

		do k2=k1,int(sys%nmode/10.0)

		  if ( ( k1 == k2 ) .and. ( k1 == k_in_1 ) .and. ( k1 == k_in_2 ) ) then
		  	 factor = 1
		  else
		  	 factor = 3
		  end if

			 nk_4_photon_up_k1_k2( k1 ) = nk_4_photon_up_k1_k2( k1 ) + factor*four_photon_amp_up_norm2( k2 )

		end do
		
	 end do

	 nk_4_photon_up_k1_k2 = 92*nk_4_photon_up_k1_k2/sys%dk1

	 
  END FUNCTION

!> -------------------------------------------------------------------------
  !> FUNCTION: nk_5_photon_at_k_up
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Computes a highly constrained slice of the five-photon momentum spectrum 
  !>   projected onto the "up" manifold. By fixing one of the five photon momenta 
  !>   at a specific grid index (`n_k_in`) and integrating over the remaining 
  !>   four variables, it isolates specific high-order multiphoton scattering 
  !>   resonances. Includes combinatoric permutations (1/120 = 1/5!) and 
  !>   symmetry factors (1, 3, 6) for indistinguishable bosons.
  !> Arguments:
  !>   - sys    : Parameter structure.
  !>   - st     : Current state object.
  !>   - n_k_in : The fixed momentum grid index for the reference photon.
  !> Return:
  !>   - real(rl) : A scalar representing the total 5-photon probability density 
  !>                associated with the emission of at least one photon at `n_k_in`.
  !>
  FUNCTION nk_5_photon_at_k_up(sys,st,n_k_in)

	 type(param), intent(in)    	 ::  sys
	 type(state), intent(in)		 ::  st
	 integer, intent(in)				 ::  n_k_in
	 real(rl)							 ::  nk_5_photon_at_k_up
	 complex(cx)					  	 ::  five_photon_amp_up( -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) &
														, -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) &
														, -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) )
	 real(rl)						  	 ::  five_photon_amp_up_norm2( -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) &
														, -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) &
														, -int(sys%nmode/10.0)+1:int(sys%nmode/10.0) )
	 complex(cx) 						 ::  fk( st%np, -sys%nmode+1:sys%nmode )
	 complex(cx) 						 ::  null_arr( -sys%nmode+1:sys%nmode )
	 complex(cx)						 ::  ov_0f( st%np )
	 integer								 ::  j,n,k1,k2,k3,k4,factor

	 nk_5_photon_at_k_up = 0._rl
    five_photon_amp_up = 0._cx
    five_photon_amp_up_norm2 = 0._rl
    null_arr = 0._cx
    factor=0

	 fk(:,1:sys%nmode) = sqrt(0.5_rl)*( st%f + st%fo )
	 do j = -sys%nmode +1,0
		fk(:,j) = sqrt(0.5_rl)*( st%f(:,(abs(j)+1)) - st%fo(:,(abs(j)+1)) )
	 end do

	 do n=1,st%np
		ov_0f(n) = ov( null_arr , fk(n,:) )
	 end do

	 do k1=-int(sys%nmode/10.0)+1,int(sys%nmode/10.0)


		do k2=-int(sys%nmode/10.0)+1,int(sys%nmode/10.0)

		  do k3=k2,int(sys%nmode/10.0)

			 do k4=k3,int(sys%nmode/10.0)

				five_photon_amp_up(k2,k3,k4) = &
				  sum( st%p(:) * fk(:,n_k_in) * fk(:,k1) * fk(:,k2) * fk(:,k3) * fk(:,k4) * ov_0f(:) ) / 120._rl

			 end do

		  end do

		end do

		five_photon_amp_up_norm2 = real( five_photon_amp_up*conjg(five_photon_amp_up) )

		do k2=-int(sys%nmode/10.0)+1,int(sys%nmode/10.0)

		  do k3=k2,int(sys%nmode/10.0)

			 do k4=k3,int(sys%nmode/10.0)

	   		if ( ( k2 .ne. k3 ) .and. (k2 .ne. k4) .and. ( k3 .ne. k4) ) then
	   			 factor = 6
	   		else if ( ( ( k2 == k3 ) .and. (k2 .ne. k4) .and. ( k3 .ne. k4) ) &
	   		         .or. ( ( k2 .ne. k3 ) .and. (k2 == k4) .and. ( k3 .ne. k4) ) &
	   		         .or. ( ( k2 .ne. k3 ) .and. (k2 .ne. k4) .and. ( k3 == k4) ) ) then
	   			 factor = 3
	   		else if ( ( k2 == k3 ) .and. (k2 == k4) .and. ( k3 == k4) )then
	   		  	 factor = 1
				end if

				nk_5_photon_at_k_up = nk_5_photon_at_k_up + factor*five_photon_amp_up_norm2(k2,k3,k4)

			 end do

		  end do

		end do
		
	 end do

	 nk_5_photon_at_k_up = 600*nk_5_photon_at_k_up/sys%dk1

	 
  END FUNCTION

  !======================================================
  !== MATH FUNCTIONS
  !======================================================

!> -------------------------------------------------------------------------
  !> FUNCTION: factorial
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Standard mathematical utility function used internally to compute the 
  !>   factorial of an integer (`n!`). Used primarily to evaluate the correct 
  !>   bosonic permutation and normalization factors required when projecting 
  !>   the multi-polaron ansatz onto higher-order N-photon Fock states.
  !> Arguments:
  !>   - n : The integer input.
  !> Return:
  !>   - integer : The calculated factorial, `n!`.
  !>
  FUNCTION factorial(n)

	 integer, intent(in)		::  n
	 integer						::  factorial
	 integer						::  i 

	 factorial=1
	 if (n > 1) then
		do i=1,n
		  factorial = factorial*i
		end do
	 end if

  END FUNCTION


END MODULE SYSTM 
