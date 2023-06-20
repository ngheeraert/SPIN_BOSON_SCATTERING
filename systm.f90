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

  !== initialize the parameters and get those in the command line
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

  !== Initialisation routines
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

  !=======================================
  !== Calculation of the derivatives
  !======================================

	 !-- Calculate the state derivatives form the st variables with KRYLOV
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

	 !-- Derivatives of E with respect to sepcified variable STARRED
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

  !=======================================
  !== State functions
  !======================================
  
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
  SUBROUTINE normalise(st)

	 type(state), intent(in out)  ::  st
	 real(rl)							:: normval

	 normval=norm(st)

	 st%p = st%p/normval
	 st%q = st%q/normval

  END SUBROUTINE
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

  !-- Calculate the overlap between two coherent states
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
  FUNCTION ov_scalar(f1,f2)

        complex(cx), intent(in) :: f1, f2
        complex(cx)             :: ov_scalar
        complex(cx)				  :: tmp1, tmp2, tmp3

		  tmp1 = conjg(f1)*f1
		  tmp2 = conjg(f2)*f2
		  tmp3 = conjg(f1)*f2

        ov_scalar = exp( -0.5_rl*tmp1 - 0.5_rl*tmp2 + tmp3 ) 

  END FUNCTION	ov_scalar
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

  !-- Probability of being in the up or down state
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

  !-- Expectation value of sigmaX
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
  FUNCTION sigmaZ(st)
  
	 type(state), intent(in)  :: st
	 real(rl) 					  :: sigmaZ

	 sigmaZ = upProb(st) - downProb(st)
  
  END FUNCTION
  FUNCTION d_sigmaZ(ost,st,dt)
  
	 type(state), intent(in)  :: ost, st
	 real(rl), intent(in)	  :: dt
	 real(rl) 					  :: d_sigmaZ

	 d_sigmaZ = ( sigmaZ(st) - sigmaZ(ost) ) / dt
  
  END FUNCTION
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

  !-- Fourrier routine EO
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

  !-- total number of photons
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
!
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
!
!  SUBROUTINE calcDerivatives(sys,st)
!
!	 type(param), intent(in)                         :: sys
!	 type(state), intent(in out)                     :: st
!	 complex(cx), dimension(st%np)                   :: bigP,bigQ
!	 complex(cx), dimension(st%np, sys%nmode)        :: bigF, bigH
!	 complex(cx), dimension(st%np,st%np)             :: a_f, a_h, b_f, b_h
!	 complex(cx), dimension(st%np, sys%nmode)        :: rhsA_f, rhsA_h
!	 complex(cx), dimension(st%np**2)                :: FullRHS_f, FullRHS_h, FullRHS_f_old, FullRHS_f_diff, FullRHS_f_new
!	 complex(cx), dimension(st%np**2,st%np**2)       :: packed_D_f, packed_D_h, packed_D_f_old
!	 complex(cx), dimension(st%np,st%np)             :: inv_ov_ff,inv_ov_hh, ov_ff_inv_noinv
!	 integer                                         :: info,i,j,k,n,m
!	 complex(cx), dimension(st%np,st%np)             :: tempMatrix_f, tempMatrix_h, tempMatrix_f_bis, tempMatrix_h_bis
!	 complex(cx), dimension(st%np)                   :: tempMatrixTer_f, tempMatrixTer_h
!	 complex(cx), dimension(st%np, sys%nmode)        :: tempfDot, temphDot
!	 integer, dimension(st%np,st%np)						 :: KD
!	 real(rl)													 :: tmp, norm2_check, norm2_inv_check
!	 complex(cx)												 :: tmp_cx
!
!	 KD=0
!	 do i=1, st%np
!		KD(i,i)=1
!	 end do
!	 b_f=0._rl
!	 b_h=0._rl
!	 rhsA_f=0._rl
!	 rhsA_h=0._rl
!	 info=0
!	 do i=1,st%np
!		bigP(i) = P_j(sys,st,i)
!		bigQ(i) = Q_j(sys,st,i)
!		bigF(i,:) =  F_j(sys,st,i)
!		bigH(i,:) =  H_j(sys,st,i)
!	 end do
!
!	 !-- invert overlap matrices
!	 !print*, 'invert overlap matrix'
!	 inv_ov_ff=st%ov_ff
!	 inv_ov_hh=st%ov_hh
!	 CALL invertH(inv_ov_ff,info)
!	 CALL invertH(inv_ov_hh,info)
!	 ov_ff_inv_noinv = inv_ov_ff .matprod. st%ov_ff
!
!	 norm2_inv_check = 0._rl
!	 do i=1,size(ov_ff_inv_noinv,1)
!	   do j=1,size(ov_ff_inv_noinv,1)
!		  if (i==j) then
!			 norm2_inv_check = norm2_inv_check + abs(ov_ff_inv_noinv(i,j)-1._rl)**2
!		  else 
!			 norm2_inv_check = norm2_inv_check + abs(ov_ff_inv_noinv(i,j))**2
!		  end if
!	   end do
!	 end do
!
!	 !-- build b matrices
!	 !print*, 'build b matrices'
!	 b_f=matmul(CONJG(st%f),TRANSPOSE(st%f))
!	 b_h=matmul(CONJG(st%h), TRANSPOSE(st%h))
!
!	 !-- build RHS
!	 ! print*, 'build rhs'
!	 do k=1, sys%nmode
!		do n=1, st%np
!		  rhsA_f(n,k)=sum(inv_ov_ff(n,:)*(bigF(:,k)-st%f(n,k)*bigP(:)))
!		  rhsA_h(n,k)=sum(inv_ov_hh(n,:)*(bigH(:,k)-st%h(n,k)*bigQ(:)))
!		end do
!	 end do
!
!	 tempMatrix_f_bis=matmul(CONJG(st%f),TRANSPOSE(rhsA_f))
!	 tempMatrix_h_bis=matmul(CONJG(st%h),TRANSPOSE(rhsA_h))
!
!	 do i=1, st%np
!		do n=1, st%np
!		  FullRHS_f((i-1)*st%np+n)=tempMatrix_f_bis(i,n)
!		  FullRHS_h((i-1)*st%np+n)=tempMatrix_h_bis(i,n)
!		end do
!	 end do
!
!	 !-- build Big D matrix
!	 !print*, 'build big d matrix'
!	 do m=1, st%np
!		do j=1, st%np
!		  do i=1, st%np
!			 do n=1, st%np
!				packed_D_f((i-1)*st%np+n,(m-1)*st%np+j) = KD(i,j)*KD(n,m)+ inv_ov_ff(n,j)*st%ov_ff(j,m)*(b_f(i,m)-b_f(i,n))
!				packed_D_h((i-1)*st%np+n,(m-1)*st%np+j) = KD(i,j)*KD(n,m)+ inv_ov_hh(n,j)*st%ov_hh(j,m)*(b_h(i,m)-b_h(i,n))
!			 end do
!		  end do
!		end do
!	 end do
!	 !-- system inversion
!	 !print*, 'system inversion'
!	 packed_D_f_old = packed_D_f
!	 FullRHS_f_old = FullRHS_f
!	 CALL solveEq_c(packed_D_f,FullRHS_f)
!	 CALL solveEq_c(packed_D_h,FullRHS_h)
!	 do j=1, st%np
!		do m=1, st%np
!		  a_f(m,j)=FullRHS_f((m-1)*st%np+j)
!		  a_h(m,j)=FullRHS_h((m-1)*st%np+j)
!		end do
!	 end do
!	 
!	 norm2_check = 0._rl
!	 do i=1,size(FullRHS_f,1)
!	   tmp_cx = 0._rl
!	   do j=1,size(FullRHS_f,1)
!	     tmp_cx = tmp_cx + packed_D_f_old(i,j)*FullRHS_f(j)
!	   end do
!		FullRHS_f_new(i) = tmp_cx
!	   norm2_check = norm2_check + abs(FullRHS_f_old(i) - FullRHS_f_new(i))**2
!    end do
!	 !FullRHS_f_diff = FullRHS_f_old - FullRHS_f_new
!	 write(200,*) st%t,norm2_inv_check, norm2_check
!	 !print*, st%t, real(abs(FullRHS_f_diff))
!	 !print*, st%t, real(abs(FullRHS_f_old))
!	 !print*,"------------------------"
!	 !print*,
!
!
!	 !-- fDot and hdot extraction
!	 !print*, 'fDot and hdot'
!
!	 tempMatrix_f=inv_ov_ff .matprod. st%ov_ff*TRANSPOSE(a_f)
!	 tempMatrix_h=inv_ov_hh .matprod. st%ov_hh*TRANSPOSE(a_h)
!	 do i=1, st%np
!		tempfDot(i,:)=SUM(tempMatrix_f(i,:))*st%f(i,:)
!		temphDot(i,:)=SUM(tempMatrix_h(i,:))*st%h(i,:)
!	 end do
!	 st%fDot= rhsA_f-matmul(tempMatrix_f,st%f)+tempfDot
!	 st%hDot= rhsA_h-matmul(tempMatrix_h,st%h)+temphDot
!	 do i=1, st%np
!		st%fDot(i,:)=st%fDot(i,:)/st%p(i)
!		st%hDot(i,:)=st%hDot(i,:)/st%q(i)
!	 end do
!
!
!	 !-- evaluate pDot and qDot
!	 !print*, 'evaluate pdot and qdot'
!	 tempMatrixTer_f= MATMUL(inv_ov_ff,bigP)
!	 st%pDot= 0.5_rl*( (/ (a_f(n,n), n=1, st%np) /) + st%p*(/ (conjg(a_f(m,m)), m=1, st%np) /) /CONJG(st%p) )
!	 st%pdot=st%pDot + tempMatrixTer_f - SUM(tempMatrix_f, dim=2)
!
!	 tempMatrixTer_h= MATMUL(inv_ov_hh,bigQ)
!	 st%qDot= 0.5_rl*( (/ (a_h(m,m), m=1, st%np) /) + st%q*(/ (conjg(a_h(m,m)), m=1, st%np) /) /CONJG(st%q) )
!	 st%qdot=st%qDot + tempMatrixTer_h - SUM(tempMatrix_h, dim=2)
!
!	 if (info==1) then
!		stop "fatal error in inversion"
!	 end if
!
!  END SUBROUTINE
!!  SUBROUTINE initialise_from_file_eo_13oct(sys,st)
!
!    type(param), intent(in)			  	::  sys
!	 type(state), intent(in out)	  		::  st
!	 integer  									::  k,i,j,m,npol, nlines, io, items
!	 character(len=200)						::  fks_file,ps_file
!	 real(rl)									::  f_r,f_i,h_r,h_i,p_r,q_r,p_i,q_i,a
!	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode)	::  fks, hks
!
!	 print*, "Initialising from: ", parameterchar_13oct(sys)
!	 fks_file=trim(adjustl(sys%file_path))//"/fks_fst_"//trim(adjustl(parameterchar_13oct(sys)))//".d"
!	 ps_file=trim(adjustl(sys%file_path))//"/ps_fst_"//trim(adjustl(parameterchar_13oct(sys)))//".d"
!
!	 !== reading the number of coherent states from the file
!	 items=0
!	 open (unit=101,file=fks_file,action="read",status="old")
!	 DO
!		READ(101,'(f25.15)',iostat=io,advance='no') f_r
!		IF (io/=0) EXIT
!		items = items + 1
!	 END DO
!	 items=items-1 ! accounting for the last line
!	 CLOSE (101)
!	 npol=(items/4)
!	 print*,"-- Number of coherent states in file= ",npol
!
!	 open (unit=101,file=ps_file,action="read",status="old")
!	 open (unit=100,file=fks_file,action="read",status="old")
!	 fks = 0._rl
!	 hks = 0._rl
!	 do  k=-sys%nmode+1,sys%nmode
!		read(100,'(f25.15)',advance='no') a
!	   do i=1,npol
!	     read(100,'(2f25.15)',advance='no') f_r, h_r
!	     fks(i,k) = f_r
!	     hks(i,k) = h_r
!	   end do
!	   do i=1,npol
!	     read(100,'(2f25.15)',advance='no') f_i, h_i
!	     fks(i,k) = fks(i,k) + Ic*f_i
!	     hks(i,k) = hks(i,k) + Ic*h_i
!	   end do
!	   read(100,*)
!	 end do
!
!	 do  k=1,sys%nmode
!		st%f(:,k) = sqrt(0.5_rl)*( fks(:,k) + fks(:,-k+1) )
!		st%h(:,k) = sqrt(0.5_rl)*( hks(:,k) + hks(:,-k+1) )
!		st%fo(:,k) = sqrt(0.5_rl)*( fks(:,k) - fks(:,-k+1) )
!		st%ho(:,k) = sqrt(0.5_rl)*( hks(:,k) - hks(:,-k+1) )
!	 end do
!
!	 do i=1,npol
!		read(101,'(2f25.15)',advance='no') p_r, q_r
!		st%p(i) = p_r
!		st%q(i) = q_r
!	 end do
!	 do i=1,npol
!		read(101,'(2f25.15)',advance='no') p_i, q_i
!		st%p(i) = st%p(i) + Ic*p_i
!		st%q(i) = st%q(i) + Ic*q_i
!	 end do
!
!	 st%t = sys%tmax
!
!	 !-- updating the sums over k
!	 CALL update_sums(sys,st)
!	 CALL normalise(st)
!
!	 close(100)
!	 close(101)
!
!  END SUBROUTINE
!  SUBROUTINE initialise_from_file_eo_1oct(sys,st)
!
!    type(param), intent(in)			  	::  sys
!	 type(state), intent(in out)	  		::  st
!	 integer  									::  k,i,j,m,npol, nlines, io, items
!	 character(len=200)						::  fks_file,ps_file
!	 real(rl)									::  f_r,f_i,h_r,h_i,p_r,q_r,p_i,q_i,a
!	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode)	::  fks, hks
!
!	 print*, "Initialising from: ", parameterchar_1oct(sys)
!	 fks_file=trim(adjustl(sys%file_path))//"/fks_fst_"//trim(adjustl(parameterchar_1oct(sys)))//".d"
!	 ps_file=trim(adjustl(sys%file_path))//"/ps_fst_"//trim(adjustl(parameterchar_1oct(sys)))//".d"
!
!	 !== reading the number of coherent states from the file
!	 items=0
!	 open (unit=101,file=fks_file,action="read",status="old")
!	 DO
!		READ(101,'(f25.15)',iostat=io,advance='no') f_r
!		IF (io/=0) EXIT
!		items = items + 1
!	 END DO
!	 items=items-1 ! accounting for the last line
!	 CLOSE (101)
!	 npol=(items/4)
!	 print*,"-- Number of coherent states in file= ",npol
!
!	 open (unit=101,file=ps_file,action="read",status="old")
!	 open (unit=100,file=fks_file,action="read",status="old")
!	 fks = 0._rl
!	 hks = 0._rl
!	 do  k=-sys%nmode+1,sys%nmode
!		read(100,'(f25.15)',advance='no') a
!	   do i=1,npol
!	     read(100,'(2f25.15)',advance='no') f_r, h_r
!	     fks(i,k) = f_r
!	     hks(i,k) = h_r
!	   end do
!	   do i=1,npol
!	     read(100,'(2f25.15)',advance='no') f_i, h_i
!	     fks(i,k) = fks(i,k) + Ic*f_i
!	     hks(i,k) = hks(i,k) + Ic*h_i
!	   end do
!	   read(100,*)
!	 end do
!
!	 do  k=1,sys%nmode
!		st%f(:,k) = sqrt(0.5_rl)*( fks(:,k) + fks(:,-k+1) )
!		st%h(:,k) = sqrt(0.5_rl)*( hks(:,k) + hks(:,-k+1) )
!		st%fo(:,k) = sqrt(0.5_rl)*( fks(:,k) - fks(:,-k+1) )
!		st%ho(:,k) = sqrt(0.5_rl)*( hks(:,k) - hks(:,-k+1) )
!	 end do
!
!	 do i=1,npol
!		read(101,'(2f25.15)',advance='no') p_r, q_r
!		st%p(i) = p_r
!		st%q(i) = q_r
!	 end do
!	 do i=1,npol
!		read(101,'(2f25.15)',advance='no') p_i, q_i
!		st%p(i) = st%p(i) + Ic*p_i
!		st%q(i) = st%q(i) + Ic*q_i
!	 end do
!
!	 st%t = sys%tmax
!
!	 !-- updating the sums over k
!	 CALL update_sums(sys,st)
!	 CALL normalise(st)
!
!	 close(100)
!	 close(101)
!
!  END SUBROUTINE
!  SUBROUTINE initialise_from_file_eo_28sept(sys,st)
!
!    type(param), intent(in)			  	::  sys
!	 type(state), intent(in out)	  		::  st
!	 integer  									::  k,i,j,m,npol, nlines, io, items
!	 character(len=200)						::  fks_file,ps_file
!	 real(rl)									::  f_r,f_i,h_r,h_i,p_r,q_r,p_i,q_i,a
!	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode)	::  fks, hks
!
!	 print*, "Initialising from: ", parameterchar_28sept(sys)
!	 fks_file=trim(adjustl(sys%file_path))//"/fks_fst_"//trim(adjustl(parameterchar_28sept(sys)))//".d"
!	 ps_file=trim(adjustl(sys%file_path))//"/ps_fst_"//trim(adjustl(parameterchar_28sept(sys)))//".d"
!
!	 !== reading the number of coherent states from the file
!	 items=0
!	 open (unit=101,file=fks_file,action="read",status="old")
!	 DO
!		READ(101,'(f25.15)',iostat=io,advance='no') f_r
!		IF (io/=0) EXIT
!		items = items + 1
!	 END DO
!	 items=items-1 ! accounting for the last line
!	 CLOSE (101)
!	 npol=(items/4)
!	 print*,"-- Number of coherent states in file= ",npol
!
!	 open (unit=101,file=ps_file,action="read",status="old")
!	 open (unit=100,file=fks_file,action="read",status="old")
!	 fks = 0._rl
!	 hks = 0._rl
!	 do  k=-sys%nmode+1,sys%nmode
!		read(100,'(f25.15)',advance='no') a
!	   do i=1,npol
!	     read(100,'(2f25.15)',advance='no') f_r, h_r
!	     fks(i,k) = f_r
!	     hks(i,k) = h_r
!	   end do
!	   do i=1,npol
!	     read(100,'(2f25.15)',advance='no') f_i, h_i
!	     fks(i,k) = fks(i,k) + Ic*f_i
!	     hks(i,k) = hks(i,k) + Ic*h_i
!	   end do
!	   read(100,*)
!	 end do
!
!	 do  k=1,sys%nmode
!		st%f(:,k) = sqrt(0.5_rl)*( fks(:,k) + fks(:,-k+1) )
!		st%h(:,k) = sqrt(0.5_rl)*( hks(:,k) + hks(:,-k+1) )
!		st%fo(:,k) = sqrt(0.5_rl)*( fks(:,k) - fks(:,-k+1) )
!		st%ho(:,k) = sqrt(0.5_rl)*( hks(:,k) - hks(:,-k+1) )
!	 end do
!
!	 do i=1,npol
!		read(101,'(2f25.15)',advance='no') p_r, q_r
!		st%p(i) = p_r
!		st%q(i) = q_r
!	 end do
!	 do i=1,npol
!		read(101,'(2f25.15)',advance='no') p_i, q_i
!		st%p(i) = st%p(i) + Ic*p_i
!		st%q(i) = st%q(i) + Ic*q_i
!	 end do
!
!	 st%t = sys%tmax
!
!	 !-- updating the sums over k
!	 CALL update_sums(sys,st)
!	 CALL normalise(st)
!
!	 close(100)
!	 close(101)
!
!  END SUBROUTINE
!  SUBROUTINE initialise_from_file_eo_2tref(sys,st)
!
!    type(param), intent(in)			  	::  sys
!	 type(state), intent(in out)	  		::  st
!	 integer  									::  k,i,j,m
!	 character(len=200)						::  fks_file,ps_file
!	 real(rl)									::  f_r,f_i,h_r,h_i,p_r,q_r,p_i,q_i,a
!	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode)	::  fks, hks
!
!	 print*, "Initialising from: ", parameterchar_2tref(sys)
!	 fks_file=trim(adjustl(sys%file_path))//"/fks_fst_"//trim(adjustl(parameterchar_2tref(sys)))//".d"
!	 ps_file=trim(adjustl(sys%file_path))//"/ps_fst_"//trim(adjustl(parameterchar_2tref(sys)))//".d"
!	 open (unit=101,file=ps_file,action="read",status="old")
!	 open (unit=100,file=fks_file,action="read",status="old")
!
!
!	 fks = 0._rl
!	 hks = 0._rl
!	 do  k=-sys%nmode+1,sys%nmode
!		read(100,'(f25.15)',advance='no') a
!	   do i=1,st%np
!	     read(100,'(2f25.15)',advance='no') f_r, h_r
!	     fks(i,k) = f_r
!	     hks(i,k) = h_r
!	   end do
!	   do i=1,st%np
!	     read(100,'(2f25.15)',advance='no') f_i, h_i
!	     fks(i,k) = fks(i,k) + Ic*f_i
!	     hks(i,k) = hks(i,k) + Ic*h_i
!	   end do
!	   read(100,*)
!	 end do
!
!	 do  k=1,sys%nmode
!		st%f(:,k) = sqrt(0.5_rl)*( fks(:,k) + fks(:,-k+1) )
!		st%h(:,k) = sqrt(0.5_rl)*( hks(:,k) + hks(:,-k+1) )
!		st%fo(:,k) = sqrt(0.5_rl)*( fks(:,k) - fks(:,-k+1) )
!		st%ho(:,k) = sqrt(0.5_rl)*( hks(:,k) - hks(:,-k+1) )
!	 end do
!
!	 do i=1,st%np
!		read(101,'(2f25.15)',advance='no') p_r, q_r
!		st%p(i) = p_r
!		st%q(i) = q_r
!	 end do
!	 do i=1,st%np
!		read(101,'(2f25.15)',advance='no') p_i, q_i
!		st%p(i) = st%p(i) + Ic*p_i
!		st%q(i) = st%q(i) + Ic*q_i
!	 end do
!
!	 st%t = sys%tmax
!
!	 !-- updating the sums over k
!	 CALL update_sums(sys,st)
!	 CALL normalise(st)
!
!	 close(100)
!	 close(101)
!
!  END SUBROUTINE
!  SUBROUTINE initialise_from_file_eo_3tref(sys,st)
!
!    type(param), intent(in)			  	::  sys
!	 type(state), intent(in out)	  		::  st
!	 integer  									::  k,i,j,m
!	 character(len=200)						::  fks_file,ps_file
!	 real(rl)									::  f_r,f_i,h_r,h_i,p_r,q_r,p_i,q_i,a
!	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode)	::  fks, hks
!
!	 print*, "Initialising from: ", parameterchar_3tref(sys)
!	 fks_file=trim(adjustl(sys%file_path))//"/fks_fst_"//trim(adjustl(parameterchar_3tref(sys)))//".d"
!	 ps_file=trim(adjustl(sys%file_path))//"/ps_fst_"//trim(adjustl(parameterchar_3tref(sys)))//".d"
!	 open (unit=101,file=ps_file,action="read",status="old")
!	 open (unit=100,file=fks_file,action="read",status="old")
!
!
!	 fks = 0._rl
!	 hks = 0._rl
!	 do  k=-sys%nmode+1,sys%nmode
!		read(100,'(f25.15)',advance='no') a
!	   do i=1,st%np
!	     read(100,'(2f25.15)',advance='no') f_r, h_r
!	     fks(i,k) = f_r
!	     hks(i,k) = h_r
!	   end do
!	   do i=1,st%np
!	     read(100,'(2f25.15)',advance='no') f_i, h_i
!	     fks(i,k) = fks(i,k) + Ic*f_i
!	     hks(i,k) = hks(i,k) + Ic*h_i
!	   end do
!	   read(100,*)
!	 end do
!
!	 do  k=1,sys%nmode
!		st%f(:,k) = sqrt(0.5_rl)*( fks(:,k) + fks(:,-k+1) )
!		st%h(:,k) = sqrt(0.5_rl)*( hks(:,k) + hks(:,-k+1) )
!		st%fo(:,k) = sqrt(0.5_rl)*( fks(:,k) - fks(:,-k+1) )
!		st%ho(:,k) = sqrt(0.5_rl)*( hks(:,k) - hks(:,-k+1) )
!	 end do
!
!	 do i=1,st%np
!		read(101,'(2f25.15)',advance='no') p_r, q_r
!		st%p(i) = p_r
!		st%q(i) = q_r
!	 end do
!	 do i=1,st%np
!		read(101,'(2f25.15)',advance='no') p_i, q_i
!		st%p(i) = st%p(i) + Ic*p_i
!		st%q(i) = st%q(i) + Ic*q_i
!	 end do
!
!	 st%t = sys%tmax
!
!	 !-- updating the sums over k
!	 CALL update_sums(sys,st)
!	 CALL normalise(st)
!
!	 close(100)
!	 close(101)
!
!  END SUBROUTINE
!	 FUNCTION parameterchar_13oct(sys)
!
!		!-- without rand_dt
!		type(param), intent(in)		::   sys
!		character(len=100)      		:: delchar,alChar,npiniChar,nmChar,npaddChar,&
!														tmaxchar,dtchar,aachar,p1ichar,sgnchar,&
!														fstepchar, merrchar,trefchar,trefchar2,trefchar3, &
!														p0char,wigxminchar,wigxmaxchar,tac1char,tac2char, &
!														wmaxchar,wcchar,alchar3,k0char,nchar,rtaddchar,&
!														x0char,prepchar,addstylechar,sigchar,dkratiochar,&
!														dkcchar,dkrchar,inter10char,nptrefchangechar
!		character(len=200)				:: parameterchar_13oct, addchar, scatchar, interchar
!
!		write(delchar, '(f6.1)') sys%del
!		write(sigchar, '(f9.4)') sys%sigma
!		write(alChar, '(f8.4)') sys%alpha
!		write(alChar3, '(f6.2)') sys%alpha*3.0
!		write(p1iChar, '(f6.2)') sys%p1i
!		write(npinichar, '(i2)') sys%npini
!		write(npaddchar, '(i2)') sys%npadd
!		write(rtaddchar, '(i2)') sys%rtadd
!		write(addstylechar, '(i2)') sys%adding_style
!		write(p0char, '(f6.3)') sys%p0*1000._rl
!		write(merrchar, '(f6.3)') sys%merr*1000000._rl
!		write(sgnchar, '(i2)') sys%sgn
!		write(nmChar, '(I5)') sys%nmode
!		write(tmaxchar, '(I10)') int(sys%tmax)
!		write(trefchar, '(f7.2)') sys%tref
!		write(trefchar2, '(f6.1)') sys%tref2
!		write(nptrefchangechar, '(I3)') sys%nptrefchange
!		write(wigxminchar, '(I5)') int(sys%wigxmin)
!		write(wigxmaxchar, '(I5)') int(sys%wigxmax)
!		write(dtchar, '(f6.4)') sys%dt
!		write(aachar, '(I10)') int(1._rl/sys%A)
!		write(wmaxchar, '(I4)') int(sys%wmax)
!		write(wcchar, '(I4)') int(sys%wc)
!		write(k0char, '(f10.4)') sys%k0
!		write(x0char, '(I8)') int(abs(sys%x0))
!		write(nchar, '(f6.2)') sys%n_wp
!		write(prepChar, '(I2)') sys%prep
!		write(dkrchar,'(I3)') sys%dk_ratio
!		write(dkcchar,'(f4.2)') sys%dk_change
!		write(inter10char,'(I1)') sys%inter
!
!		addchar="_"
!		scatchar="_"
!		interchar="_"
!		if (sys%inter == 1) then
!		  interchar="_dkc"//trim(adjustl(dkcchar))//&
!						"_dkr"//trim(adjustl(dkrchar))//"_"
!		end if
!		if (sys%npadd .ne. 0) then
!		  addchar="_tr"//trim(adjustl(trefchar))//&
!						 !"_"//trim(adjustl(tac1char))//&
!						 "_me"//trim(adjustl(merrchar))//&
!						 "_p"//trim(adjustl(p0char))//"_"
!		end if
!		if ( (sys%prep .ge. 50) .and. (sys%prep < 100) ) then
!		  scatchar="_k"//trim(adjustl(k0char))//&
!					 "_x"//trim(adjustl(x0char))//&
!					 "_sig"//trim(adjustl(sigchar))//&
!					 "_n"//trim(adjustl(nchar))
!		end if
!
!		parameterchar_13oct=trim(adjustl(nmchar))//"m"//&
!					 "_I"//trim(adjustl(inter10char))//&
!					 trim(adjustl(interchar))//&
!					 "np"//trim(adjustl(npinichar))//&
!					 "_"//trim(adjustl(npaddchar))//&
!					 "_al"//trim(adjustl(alchar))//&
!					 "_del"//trim(adjustl(delchar))//&
!					 "_dt"//trim(adjustl(dtchar))//&
!					 trim(adjustl(addchar))//&
!					 trim(adjustl(scatchar))//&
!					 "tmax"//trim(adjustl(tmaxchar))//&
!					 "_p"//trim(adjustl(prepchar))
!
!	 END FUNCTION
!	 FUNCTION parameterchar_1oct(sys)
!
!		type(param), intent(in)		::   sys
!		character(len=100)      		:: delchar,alChar,npiniChar,nmChar,npaddChar,&
!														tmaxchar,dtchar,aachar,p1ichar,sgnchar,&
!														fstepchar, merrchar,trefchar,trefchar2,trefchar3, &
!														p0char,wigxminchar,wigxmaxchar,tac1char,tac2char, &
!														wmaxchar,wcchar,alchar3,k0char,nchar,rtaddchar,&
!														x0char,prepchar,addstylechar,sigchar,dkratiochar,&
!														dkcchar,dkrchar,inter10char,nptrefchangechar
!		character(len=200)				:: parameterChar_1oct, addchar, scatchar, interchar
!
!		write(delchar, '(f6.1)') sys%del
!		write(sigchar, '(f8.3)') sys%sigma
!		write(alChar, '(f8.4)') sys%alpha
!		write(alChar3, '(f6.2)') sys%alpha*3.0
!		write(p1iChar, '(f6.2)') sys%p1i
!		write(npinichar, '(i2)') sys%npini
!		write(npaddchar, '(i2)') sys%npadd
!		write(p0char, '(f6.3)') sys%p0*1000._rl
!		write(merrchar, '(f6.3)') sys%merr*1000000._rl
!		write(nmChar, '(I5)') sys%nmode
!		write(tmaxchar, '(I10)') int(sys%tmax)
!		write(trefchar, '(f7.2)') sys%tref
!		write(wigxminchar, '(I5)') int(sys%wigxmin)
!		write(wigxmaxchar, '(I5)') int(sys%wigxmax)
!		write(dtchar, '(f6.4)') sys%dt
!		write(wmaxchar, '(I4)') int(sys%wmax)
!		write(wcchar, '(I4)') int(sys%wc)
!		write(k0char, '(f10.4)') sys%k0
!		write(x0char, '(I8)') int(abs(sys%x0))
!		write(nchar, '(f6.2)') sys%n_wp
!		write(prepChar, '(I2)') sys%prep
!
!		addchar="_"
!		scatchar="_"
!		interchar="_"
!		if (sys%inter == 1) then
!		  interchar="_dkc"//trim(adjustl(dkcchar))//&
!						"_dkr"//trim(adjustl(dkrchar))//"_"
!		end if
!		if (sys%npadd .ne. 0) then
!		  addchar="_tr"//trim(adjustl(trefchar))//&
!						 !"_"//trim(adjustl(tac1char))//&
!						 "_me"//trim(adjustl(merrchar))//&
!						 "_p"//trim(adjustl(p0char))//"_"
!		end if
!		if ( (sys%prep .ge. 50) .and. (sys%prep < 100) ) then
!		  scatchar="_k"//trim(adjustl(k0char))//&
!					 "_x"//trim(adjustl(x0char))//&
!					 "_sig"//trim(adjustl(sigchar))//&
!					 "_n"//trim(adjustl(nchar))
!		end if
!
!		parameterchar_1oct=trim(adjustl(nmchar))//"m"//&
!					 "_I"//trim(adjustl(inter10char))//&
!					 trim(adjustl(interchar))//&
!					 "np"//trim(adjustl(npinichar))//&
!					 "_"//trim(adjustl(npaddchar))//&
!					 "_al"//trim(adjustl(alchar))//&
!					 "_del"//trim(adjustl(delchar))//&
!					 "_dt"//trim(adjustl(dtchar))//&
!					 trim(adjustl(addchar))//&
!					 trim(adjustl(scatchar))//&
!					 "tmax"//trim(adjustl(tmaxchar))//&
!					 "_p"//trim(adjustl(prepchar))
!
!	 END FUNCTION
!	 FUNCTION parameterchar_28sept(sys)
!
!		type(param), intent(in)		::   sys
!		character(len=100)      		:: delchar,alChar,npiniChar,nmChar,npaddChar,&
!														tmaxchar,dtchar,aachar,p1ichar,sgnchar,&
!														fstepchar, merrchar,trefchar,trefchar2,trefchar3, &
!														p0char,wigxminchar,wigxmaxchar,tac1char,tac2char, &
!														wmaxchar,wcchar,alchar3,k0char,nchar,rtaddchar,&
!														x0char,prepchar,addstylechar,sigchar,dkratiochar,&
!														dkcchar,dkrchar,inter10char,nptrefchangechar
!		character(len=200)				:: parameterChar_28sept, addchar, scatchar, interchar
!
!		write(delchar, '(f6.1)') sys%del
!		write(sigchar, '(f8.3)') sys%sigma
!		write(alChar, '(f8.4)') sys%alpha
!		write(alChar3, '(f6.2)') sys%alpha*3.0
!		write(p1iChar, '(f6.2)') sys%p1i
!		write(npinichar, '(i2)') sys%npini
!		write(npaddchar, '(i2)') sys%npadd
!		write(rtaddchar, '(i2)') sys%rtadd
!		write(addstylechar, '(i2)') sys%adding_style
!		write(p0char, '(f6.3)') sys%p0*1000._rl
!		write(merrchar, '(f6.3)') sys%merr*1000000._rl
!		write(sgnchar, '(i2)') sys%sgn
!		write(nmChar, '(I5)') sys%nmode
!		write(tmaxchar, '(I10)') int(sys%tmax)
!		write(trefchar, '(f6.1)') sys%tref
!		write(trefchar2, '(f6.1)') sys%tref2
!		write(nptrefchangechar, '(I3)') sys%nptrefchange
!		write(wigxminchar, '(I5)') int(sys%wigxmin)
!		write(wigxmaxchar, '(I5)') int(sys%wigxmax)
!		write(dtchar, '(f6.4)') sys%dt
!		write(aachar, '(I10)') int(1._rl/sys%A)
!		write(wmaxchar, '(I4)') int(sys%wmax)
!		write(wcchar, '(I4)') int(sys%wc)
!		write(k0char, '(f10.4)') sys%k0
!		write(x0char, '(f5.1)') abs(sys%x0/1000._rl)
!		write(nchar, '(f6.2)') sys%n_wp
!		write(prepChar, '(I2)') sys%prep
!		write(dkrchar,'(I3)') sys%dk_ratio
!		write(dkcchar,'(f4.2)') sys%dk_change
!		write(inter10char,'(I1)') sys%inter
!
!		addchar="_"
!		scatchar="_"
!		interchar="_"
!		if (sys%inter == 1) then
!		  interchar="_dkc"//trim(adjustl(dkcchar))//&
!						"_dkr"//trim(adjustl(dkrchar))//"_"
!		end if
!		if (sys%npadd .ne. 0) then
!		  addchar="_tr"//trim(adjustl(trefchar))//&
!						 !"_"//trim(adjustl(tac1char))//&
!						 "_me"//trim(adjustl(merrchar))//&
!						 "_p"//trim(adjustl(p0char))//"_"
!		end if
!		if ( (sys%prep .ge. 50) .and. (sys%prep < 100) ) then
!		  scatchar="_k"//trim(adjustl(k0char))//&
!					 "_x"//trim(adjustl(x0char))//&
!					 "_sig"//trim(adjustl(sigchar))//&
!					 "_n"//trim(adjustl(nchar))
!		end if
!
!		parameterchar_28sept=trim(adjustl(nmchar))//"m"//&
!					 "_I"//trim(adjustl(inter10char))//&
!					 trim(adjustl(interchar))//&
!					 "np"//trim(adjustl(npinichar))//&
!					 "_"//trim(adjustl(npaddchar))//&
!					 "_al"//trim(adjustl(alchar))//&
!					 "_del"//trim(adjustl(delchar))//&
!					 "_dt"//trim(adjustl(dtchar))//&
!					 trim(adjustl(addchar))//&
!					 trim(adjustl(scatchar))//&
!					 "tmax"//trim(adjustl(tmaxchar))//&
!					 "_p"//trim(adjustl(prepchar))
!
!	 END FUNCTION
!	 FUNCTION parameterchar_3tref(sys)
!
!		type(param), intent(in)		::   sys
!		character(len=100)      		:: delchar,alChar,npiniChar,nmChar,npaddChar,&
!														tmaxchar,dtchar,aachar,p1ichar,sgnchar,&
!														fstepchar, merrchar,trefchar,trefchar2,trefchar3, &
!														p0char,wigxminchar,wigxmaxchar,tac1char,tac2char, &
!														wmaxchar,wcchar,alchar3,k0char,nchar,rtaddchar,&
!														x0char,prepchar,addstylechar,sigchar,dkratiochar,&
!														dkcchar,dkrchar,inter10char,nptrefchangechar,nptrefchangechar2
!		character(len=200)				:: parameterChar_3tref, addchar, scatchar, interchar
!
!		write(delchar, '(f6.1)') sys%del
!		write(sigchar, '(f8.4)') sys%sigma
!		write(alChar, '(f8.4)') sys%alpha
!		write(alChar3, '(f6.2)') sys%alpha*3.0
!		write(p1iChar, '(f6.2)') sys%p1i
!		write(npinichar, '(i2)') sys%npini
!		write(npaddchar, '(i2)') sys%npadd
!		write(rtaddchar, '(i2)') sys%rtadd
!		write(addstylechar, '(i2)') sys%adding_style
!		write(p0char, '(f6.3)') sys%p0*1000._rl
!		write(merrchar, '(f6.3)') sys%merr*1000000._rl
!		write(sgnchar, '(i2)') sys%sgn
!		write(nmChar, '(I5)') sys%nmode
!		write(tmaxchar, '(I10)') int(sys%tmax)
!		write(trefchar, '(f6.1)') sys%tref
!		write(trefchar2, '(f6.1)') sys%tref2
!		write(trefchar3, '(f6.1)') sys%tref3
!		write(nptrefchangechar, '(I3)') sys%nptrefchange
!		write(nptrefchangechar2, '(I3)') sys%nptrefchange2
!		write(wigxminchar, '(I5)') int(sys%wigxmin)
!		write(wigxmaxchar, '(I5)') int(sys%wigxmax)
!		write(dtchar, '(f6.4)') sys%dt
!		write(aachar, '(I10)') int(1._rl/sys%A)
!		write(wmaxchar, '(I4)') int(sys%wmax)
!		write(wcchar, '(I4)') int(sys%wc)
!		write(k0char, '(f10.4)') sys%k0
!		write(x0char, '(f5.1)') abs(sys%x0/1000._rl)
!		write(nchar, '(f6.2)') sys%n_wp
!		write(prepChar, '(I2)') sys%prep
!		write(dkrchar,'(I3)') sys%dk_ratio
!		write(dkcchar,'(f4.2)') sys%dk_change
!		write(inter10char,'(I1)') sys%inter
!
!		addchar="_"
!		scatchar="_"
!		interchar="_"
!		if (sys%inter == 1) then
!		  interchar="_dkc"//trim(adjustl(dkcchar))//&
!						"_dkr"//trim(adjustl(dkrchar))//"_"
!		end if
!		if (sys%npadd .ne. 0) then
!		  addchar="_tr"//trim(adjustl(trefchar))//&
!						 "_"//trim(adjustl(trefchar2))//&
!						 "_"//trim(adjustl(trefchar3))//&
!						 "_"//trim(adjustl(nptrefchangechar))//&
!						 "_"//trim(adjustl(nptrefchangechar2))//&
!						 !"_"//trim(adjustl(tac1char))//&
!						 "_me"//trim(adjustl(merrchar))//&
!						 "_p"//trim(adjustl(p0char))//"_"
!		end if
!		if ( (sys%prep .ge. 50) .and. (sys%prep < 100) ) then
!		  scatchar="_k"//trim(adjustl(k0char))//&
!					 "_x"//trim(adjustl(x0char))//&
!					 "_sig"//trim(adjustl(sigchar))//&
!					 "_n"//trim(adjustl(nchar))
!		end if
!
!		parameterchar_3tref=trim(adjustl(nmchar))//"m"//&
!					 "_I"//trim(adjustl(inter10char))//&
!					 trim(adjustl(interchar))//&
!					 "np"//trim(adjustl(npinichar))//&
!					 "_"//trim(adjustl(npaddchar))//&
!					 "_al"//trim(adjustl(alchar))//&
!					 "_del"//trim(adjustl(delchar))//&
!					 "_dt"//trim(adjustl(dtchar))//&
!					 trim(adjustl(addchar))//&
!					 trim(adjustl(scatchar))//&
!					 "tmax"//trim(adjustl(tmaxchar))//&
!					 "_p"//trim(adjustl(prepchar))
!
!	 END FUNCTION
!	 FUNCTION parameterchar_no_nmode(sys)
!
!		type(param), intent(in)		::   sys
!		character(len=100)      		:: delchar,alChar,npiniChar,nmChar,npaddChar,&
!														tmaxchar,dtchar,aachar,p1ichar,sgnchar,&
!														fstepchar, merrchar,trefchar,trefchar2,trefchar3, &
!														p0char,wigxminchar,wigxmaxchar,tac1char,tac2char, &
!														wmaxchar,wcchar,alchar3,k0char,nchar,rtaddchar,&
!														x0char,prepchar,addstylechar,sigchar,dkratiochar,&
!														dkcchar,dkrchar,inter10char,nptrefchangechar,nptrefchangechar2
!		character(len=200)				:: parameterChar_no_nmode, addchar, scatchar, interchar
!
!		write(delchar, '(f6.1)') sys%del
!		write(sigchar, '(f8.4)') sys%sigma
!		write(alChar, '(f8.4)') sys%alpha
!		write(alChar3, '(f6.2)') sys%alpha*3.0
!		write(p1iChar, '(f6.2)') sys%p1i
!		write(npinichar, '(i2)') sys%npini
!		write(npaddchar, '(i2)') sys%npadd
!		write(rtaddchar, '(i2)') sys%rtadd
!		write(addstylechar, '(i2)') sys%adding_style
!		write(p0char, '(f6.3)') sys%p0*1000._rl
!		write(merrchar, '(f6.3)') sys%merr*1000000._rl
!		write(sgnchar, '(i2)') sys%sgn
!		write(nmChar, '(I5)') sys%nmode
!		write(tmaxchar, '(I10)') int(sys%tmax)
!		write(trefchar, '(f6.1)') sys%tref
!		write(trefchar2, '(f6.1)') sys%tref2
!		write(trefchar3, '(f6.1)') sys%tref3
!		write(nptrefchangechar, '(I3)') sys%nptrefchange
!		write(nptrefchangechar2, '(I3)') sys%nptrefchange2
!		write(wigxminchar, '(I5)') int(sys%wigxmin)
!		write(wigxmaxchar, '(I5)') int(sys%wigxmax)
!		write(dtchar, '(f6.4)') sys%dt
!		write(aachar, '(I10)') int(1._rl/sys%A)
!		write(wmaxchar, '(I4)') int(sys%wmax)
!		write(wcchar, '(I4)') int(sys%wc)
!		write(k0char, '(f10.4)') sys%k0
!		write(x0char, '(f5.1)') abs(sys%x0/1000._rl)
!		write(nchar, '(f6.2)') sys%n_wp
!		write(prepChar, '(I2)') sys%prep
!		write(dkrchar,'(I3)') sys%dk_ratio
!		write(dkcchar,'(f4.2)') sys%dk_change
!		write(inter10char,'(I1)') sys%inter
!
!		addchar="_"
!		scatchar="_"
!		interchar="_"
!		if (sys%inter == 1) then
!		  interchar="_dkc"//trim(adjustl(dkcchar))//&
!						"_dkr"//trim(adjustl(dkrchar))//"_"
!		end if
!		if (sys%npadd .ne. 0) then
!		  addchar="_tr"//trim(adjustl(trefchar))//&
!						 "_"//trim(adjustl(trefchar2))//&
!						 "_"//trim(adjustl(trefchar3))//&
!						 "_"//trim(adjustl(nptrefchangechar))//&
!						 "_"//trim(adjustl(nptrefchangechar2))//&
!						 !"_"//trim(adjustl(tac1char))//&
!						 "_me"//trim(adjustl(merrchar))//&
!						 "_p"//trim(adjustl(p0char))//"_"
!		end if
!		if ( (sys%prep .ge. 50) .and. (sys%prep < 100) ) then
!		  scatchar="_k"//trim(adjustl(k0char))//&
!					 "_x"//trim(adjustl(x0char))//&
!					 "_sig"//trim(adjustl(sigchar))//&
!					 "_n"//trim(adjustl(nchar))
!		end if
!
!		parameterchar_no_nmode="_I"//trim(adjustl(inter10char))//&
!					 trim(adjustl(interchar))//&
!					 "np"//trim(adjustl(npinichar))//&
!					 "_"//trim(adjustl(npaddchar))//&
!					 "_al"//trim(adjustl(alchar))//&
!					 "_del"//trim(adjustl(delchar))//&
!					 "_dt"//trim(adjustl(dtchar))//&
!					 trim(adjustl(addchar))//&
!					 trim(adjustl(scatchar))//&
!					 "tmax"//trim(adjustl(tmaxchar))//&
!					 "_p"//trim(adjustl(prepchar))
!
!	 END FUNCTION
!	 FUNCTION parameterchar_2tref(sys)
!
!		type(param), intent(in)		::   sys
!		character(len=100)      		:: delchar,alChar,npiniChar,nmChar,npaddChar,&
!														tmaxchar,dtchar,aachar,p1ichar,sgnchar,&
!														fstepchar, merrchar,trefchar,trefchar2,trefchar3, &
!														p0char,wigxminchar,wigxmaxchar,tac1char,tac2char, &
!														wmaxchar,wcchar,alchar3,k0char,nchar,rtaddchar,&
!														x0char,prepchar,addstylechar,sigchar,dkratiochar,&
!														dkcchar,dkrchar, inter10char, nptrefchangechar
!		character(len=200)				:: parameterChar_2tref, addchar, scatchar, interchar
!
!		write(delchar, '(f6.1)') sys%del
!		write(sigchar, '(f8.3)') sys%sigma
!		write(alChar, '(f8.4)') sys%alpha
!		write(alChar3, '(f6.2)') sys%alpha*3.0
!		write(p1iChar, '(f6.2)') sys%p1i
!		write(npinichar, '(i2)') sys%npini
!		write(npaddchar, '(i2)') sys%npadd
!		write(rtaddchar, '(i2)') sys%rtadd
!		write(addstylechar, '(i2)') sys%adding_style
!		write(p0char, '(f10.3)') sys%p0*1000._rl
!		write(merrchar, '(f6.3)') sys%merr*1000000._rl
!		write(sgnchar, '(i2)') sys%sgn
!		write(nmChar, '(I5)') sys%nmode
!		write(tmaxchar, '(I10)') int(sys%tmax)
!		write(trefchar, '(f6.1)') sys%tref
!		write(nptrefchangechar, '(I3)') sys%nptrefchange
!		write(trefchar2, '(f6.1)') sys%tref2
!		write(wigxminchar, '(I5)') int(sys%wigxmin)
!		write(wigxmaxchar, '(I5)') int(sys%wigxmax)
!		write(dtchar, '(f6.4)') sys%dt
!		write(aachar, '(I10)') int(1._rl/sys%A)
!		write(wmaxchar, '(I4)') int(sys%wmax)
!		write(wcchar, '(I4)') int(sys%wc)
!		write(k0char, '(f10.4)') sys%k0
!		write(x0char, '(f5.1)') abs(sys%x0/1000._rl)
!		write(nchar, '(f6.2)') sys%n_wp
!		write(prepChar, '(I2)') sys%prep
!		write(dkrchar,'(I3)') sys%dk_ratio
!		write(dkcchar,'(f4.2)') sys%dk_change
!		write(inter10char,'(I1)') sys%inter
!
!		addchar="_"
!		scatchar="_"
!		interchar="_"
!		if (sys%inter == 1) then
!		  interchar="_dkc"//trim(adjustl(dkcchar))//&
!						"_dkr"//trim(adjustl(dkrchar))//"_"
!		end if
!		if (sys%npadd .ne. 0) then
!		  addchar="_tr"//trim(adjustl(trefchar))//&
!						 "_"//trim(adjustl(trefchar2))//&
!						 "_"//trim(adjustl(nptrefchangechar))//&
!						 "_me"//trim(adjustl(merrchar))//&
!						 "_p"//trim(adjustl(p0char))//"_"
!		end if
!		if ( (sys%prep .ge. 50) .and. (sys%prep < 100) ) then
!		  scatchar="_k"//trim(adjustl(k0char))//&
!					 "_x"//trim(adjustl(x0char))//&
!					 "_sig"//trim(adjustl(sigchar))//&
!					 "_n"//trim(adjustl(nchar))
!		end if
!
!		parameterchar_2tref=trim(adjustl(nmchar))//"m"//&
!					 "_I"//trim(adjustl(inter10char))//&
!					 trim(adjustl(interchar))//&
!					 "np"//trim(adjustl(npinichar))//&
!					 "_"//trim(adjustl(npaddchar))//&
!					 "_al"//trim(adjustl(alchar))//&
!					 "_del"//trim(adjustl(delchar))//&
!					 "_dt"//trim(adjustl(dtchar))//&
!					 trim(adjustl(addchar))//&
!					 trim(adjustl(scatchar))//&
!					 "tmax"//trim(adjustl(tmaxchar))//&
!					 "_p"//trim(adjustl(prepchar))
!
!	 END FUNCTION

  !SUBROUTINE calcDerivatives_slow(sys,st) 

  !  type(param), intent(in)                            		::  sys
  !  type(state), intent(in out)	                      		::  st
  !  complex(cx), dimension(st%np)							 		::  bigP, bigQ, v_p, v_q
  !  complex(cx), dimension(st%np,sys%nmode)				 		::  bigF, bigH, v_f, v_h
  !  complex(cx), dimension(st%np,st%np)              		::  M_p, M_q, M_pi, M_qi,N_p, N_q, N_pi, N_qi
  !  complex(cx), dimension(st%np,st%np,st%np)       		::  A_p, A_q
  !  complex(cx), dimension(st%np,st%np,st%np,sys%nmode)  ::  A_f, A_h
  !  real(rl), dimension(2*st%np**2,2*st%np**2)   	 		::  EqMat_p, EqMat_q
  !  real(rl), dimension(2*st%np**2)    		 			 		::  EqRhs_p, EqRhs_q, kappaVec_p, kappaVec_q
  !  real(rl), dimension(st%np,st%np)	   				 		::  Kappa_p_r,Kappa_q_r,kappa_p_c,kappa_q_c
  !  complex(cx), dimension(st%np,st%np)					 		::  Kappa_p,Kappa_q  !-- FINAL KAPPA MATRICES
  !  complex(cx),dimension(st%np,sys%nmode)  			 		::  f,h,fc,hc  		
  !  complex(cx),dimension(st%np)			   			 		::  p,q,pc,qc
  !  complex(cx),dimension(st%np)			   			 		::  der_E_pc_save,der_E_qc_save
  !  complex(cx),dimension(st%np,sys%nmode)			   	::  der_E_fc_save,der_E_hc_save
  !  integer												          		::  i,j,k,l,m,cj,ii,jj,np,info,s
  !  
  !  complex(cx), dimension(st%np) 				          		::  tmpV1_p,tmpV1_q,tmpV2_p,tmpV2_q
  !  complex(cx), dimension(st%np,st%np,st%np,sys%nmode) 	::  tmpA_f1,tmpA_h1,tmpA_f2,tmpA_h2
  !  complex(cx), dimension(st%np,sys%nmode) 				   ::  tmpV1_f,tmpV1_h
  !  complex(cx), dimension(st%np,st%np,st%np,st%np)		::  D_f,D_h
  !  complex(cx), dimension(st%np,st%np,st%np)				::  E_f,E_h
  !  complex(cx), dimension(st%np,st%np)						::  K_f,K_h


  !  !-- A few shortcuts
  !  f=st%f;h=st%h;p=st%p;q=st%q
  !  fc=conjg(f);hc=conjg(h);pc=conjg(p);qc=conjg(q)
  !  np = st%np

  !  !- set all variables to zero
  !  bigP = 0._cx;bigQ = 0._cx;bigF = 0._cx;bigH = 0._cx
  !  M_p=0._cx;M_q=0._cx;N_p=0._cx;N_q=0._cx
  !  M_pi=0._cx;M_qi=0._cx;N_pi=0._cx;N_qi=0._cx
  !  A_p=0._cx;A_q=0._cx;A_f=0._cx;A_h=0._cx
  !  eqMat_p=0._cx;eqMat_q=0._cx;eqrhs_p=0._cx;eqrhs_q=0._cx
  !	 eqRhs_p = 0._cx;eqRhs_q = 0._cx
  !  v_p=0._cx;v_q=0._cx;v_f=0._cx;v_h=0._cx
  !  Kappa_p_r=0._cx;Kappa_q_r=0._cx;kappa_p_c=0._cx;kappa_q_c=0._cx
  !  kappaVec_p=0._cx; kappaVec_q=0._cx;Kappa_p=0._cx;Kappa_q=0._cx
  !  der_E_fc_save(:,:) = 0._cx
  !  der_E_hc_save(:,:) = 0._cx
  !  der_E_pc_save(:) = 0._cx
  !  der_E_qc_save(:) = 0._cx
  !  D_f = 0._cx;D_h = 0._cx
  !  E_f = 0._cx;E_h = 0._cx
  !  K_f = 0._cx;K_h = 0._cx

  !  !==================================================
  !  !-- STEP 1: Computation of A_p,A_q and v_p,v_q
  !  !-- pDot(i) = sum_(j,k) [ A_p(i,j,k)*Kappa(k,j) + v_p(i) ] 
  !  !-- qDot(i) = sum_(j,k) [ A_q(i,j,k)*Kappa(k,j) + v_q(i) ] 

  !  do i=1,np
  ! 	  bigP(i) = P_j(sys,st,i)
  ! 	  bigQ(i) = Q_j(sys,st,i)
  ! 	  bigF(i,:) =  F_j(sys,st,i) 
  ! 	  bigH(i,:) =  H_j(sys,st,i) 
  !  end do

  !  !-- Matrix M: to be inverted
  !  M_p = st%ov_ff
  !  M_q = st%ov_hh

  !  !-- perform the inversion
  !  M_pi = M_p
  !  M_qi = M_q
  !  CALL invertH(M_pi,info)
  !  CALL invertH(M_qi,info)

  !  !-- Define v_p, v_q and A_p, A_q
  !  v_p = M_pi .matprod. bigP
  !  v_q = M_qi .matprod. bigQ

  !  do i=1,size(A_p,1)
  ! 	do j=1,size(A_p,2)
  ! 	  do k=1,size(A_p,3)
  ! 	    A_p(i,j,k) = 0.5_rl* M_pi(i,j) * st%ov_ff(j,k) * p(k)
  ! 	   A_q(i,j,k) = 0.5_rl* M_qi(i,j) * st%ov_hh(j,k) * q(k)
  ! 	  end do
  ! 	end do
  !  end do
  !  

  !  !==================================================
  !  !-- STEP 2: Computation of A_f,A_h and v_f,v_h
  !  !-- fDot(i) = sum_(j,k) [ A_f(i,j,k)*Kappa(k,j) ] + v_f(i)  
  !  !-- hDot(i) = sum_(j,k) [ A_h(i,j,k)*Kappa(k,j) ] + v_h(i)  

  !  do j=1,np
  ! 	do m=1,np
  ! 	  N_p(j,m) = p(m)*st%ov_ff(j,m)
  ! 	  N_q(j,m) = q(m)*st%ov_hh(j,m) 
  ! 	end do
  !  end do

  !  !-- Calculate Nij and is inverse
  !  N_pi = N_p
  !  N_qi = N_q
  !  CALL invertGeneral(N_pi,info)
  !  CALL invertGeneral(N_qi,info)
  !  

  !  !-- Computation of v_f and v_h (H_i in the notes)
  !  tmpV1_f=0._cx
  !  tmpV1_h=0._cx

  !  !-- First term of v_f (and v_h)
  !  tmpV1_f = bigF
  !  tmpV1_h = bigH
  !  !-- Second term of v_f
  !  do s=1,sys%nmode
  ! 	do j=1,np
  ! 		 tmpV1_f(j,s) = tmpV1_f(j,s) &
  ! 								 - sum( st%ov_ff(j,:)*v_p(:)*f(:,s) )
  ! 		 tmpV1_h(j,s) = tmpV1_h(j,s) &
  ! 								 - sum( st%ov_hh(j,:)*v_q(:)*h(:,s) )
  ! 	end do
  ! 	v_f(:,s) = N_pi .matprod. tmpV1_f(:,s)
  ! 	v_h(:,s) = N_qi .matprod. tmpV1_h(:,s)
  !  end do

  !  tmpA_f1=0._cx; tmpA_h1=0._cx
  !  tmpA_f2=0._cx; tmpA_h2=0._cx

  !  !-- Defining the Betref matrices, here designated by A_f and A_h

  !  do s=1,sys%nmode
  ! 	do k=1,size(A_f,2)
  ! 	  do j=1,size(A_f,3)

  ! 		 do i=1,size(A_f,1)
  ! 			tmpA_f1(i,j,k,s) = 0.5_rl* N_pi(i,j) * p(k)*f(k,s) * st%ov_ff(j,k)
  ! 			tmpA_h1(i,j,k,s) = 0.5_rl* N_qi(i,j) * q(k)*h(k,s) * st%ov_hh(j,k)
  ! 		 end do

  ! 		 do l=1,np
  ! 			  tmpA_f2(l,j,k,s) = tmpA_f2(l,j,k,s) & 
  ! 						 - sum( A_p(:,j,k)*f(:,s)*st%ov_ff(l,:) )
  ! 			  tmpA_h2(l,j,k,s) = tmpA_h2(l,j,k,s) &
  ! 					 - sum( A_q(:,j,k)*h(:,s)*st%ov_hh(l,:) )
  ! 		 end do

  ! 	  end do
  ! 	

  ! 	  A_f(:,:,k,s) = tmpA_f1(:,:,k,s) + ( N_pi .matprod. tmpA_f2(:,:,k,s) )
  ! 	  A_h(:,:,k,s) = tmpA_h1(:,:,k,s) + ( N_qi .matprod. tmpA_h2(:,:,k,s) )

  ! 	end do
  !  end do

  !  !==================================================
  !  !-- STEP 3: Solve the system of equations for Kappa
  !  !-- Define Q_ij, a 2*np**2 matrix: j in [1,np**2] correpsonds to Re(kappa)
  !  !--    while j in [1+np**2,2*np**2] corresponds to Im(Kappa)

  !  cj = st%np**2   !-- a shortcut to get to the imagainary part of kappaVec

  !  do i=1,st%np
  ! 	do j=1,st%np
  ! 	  K_f(i,j) = sum( conjg(f(i,:))*v_f(i,:) + f(i,:)*conjg(v_f(i,:)) - 2._rl*conjg(f(j,:))*v_f(i,:) )
  ! 	  K_h(i,j) = sum( conjg(h(i,:))*v_h(i,:) + h(i,:)*conjg(v_h(i,:)) - 2._rl*conjg(h(j,:))*v_h(i,:) )
  ! 	  do k=1,st%np
  ! 		 E_f(i,j,k) = sum( f(i,:)*conjg(A_f(i,j,k,:))  )
  ! 		 E_h(i,j,k) = sum( h(i,:)*conjg(A_h(i,j,k,:))  )
  ! 		 do m=1,st%np
  ! 			D_f(i,m,j,k) = sum( (conjg(f(i,:)) - 2._rl*conjg(f(m,:)))*A_f(i,j,k,:) )
  ! 			D_h(i,m,j,k) = sum( (conjg(h(i,:)) - 2._rl*conjg(h(m,:)))*A_h(i,j,k,:) )
  ! 		 end do
  ! 	  end do
  ! 	end do
  !  end do
  !  

  !  do i=1, np
  ! 	do m=1, np

  ! 	  ii = np*(i-1)+m

  ! 	  eqRhs_p(ii) 	  =  real( K_f(i,m) )
  ! 	  eqRhs_q(ii) 	  =  real( K_h(i,m) )
  ! 	  eqRhs_p(ii+cj)  =  aimag( K_f(i,m) )
  ! 	  eqRhs_q(ii+cj)  =  aimag( K_h(i,m) )

  ! 	  do j=1, np
  ! 		 do k=1, np

  ! 		   jj = np*(k-1) + j

  ! 				 eqMat_p(ii,jj) = kroneckerDelta(ii,jj) - real(  D_f(i,m,j,k) + E_f(i,j,k) )
  ! 				 eqMat_p(ii,jj+cj) = kroneckerDelta(ii,jj+cj) + aimag(  D_f(i,m,j,k) - E_f(i,j,k) )
  ! 				 eqMat_p(ii+cj,jj) = kroneckerDelta(ii+cj,jj) - aimag(  D_f(i,m,j,k) + E_f(i,j,k) )
  ! 				 eqMat_p(ii+cj,jj+cj) = kroneckerDelta(ii+cj,jj+cj) - real(  D_f(i,m,j,k) - E_f(i,j,k) )

  ! 				 eqMat_q(ii,jj) = kroneckerDelta(ii,jj) - real(  D_h(i,m,j,k) + E_h(i,j,k) )
  ! 				 eqMat_q(ii,jj+cj) = kroneckerDelta(ii,jj+cj) + aimag(  D_h(i,m,j,k) - E_h(i,j,k) )
  ! 				 eqMat_q(ii+cj,jj) = kroneckerDelta(ii+cj,jj) - aimag(  D_h(i,m,j,k) + E_h(i,j,k) )
  ! 				 eqMat_q(ii+cj,jj+cj) = kroneckerDelta(ii+cj,jj+cj) - real(  D_h(i,m,j,k) - E_h(i,j,k) )

  ! 		 end do
  ! 	  end do
  ! 	end do
  !  end do


  !  !-- Apply the Lapack algorithm to solve the real system of equations
  !  !-- the solution replaces the 2nd argument
  !  kappaVec_p = eqRhs_p
  !  CALL solveEq_r(eqMat_p,kappaVec_p)
  !  kappaVec_q = eqRhs_q
  !  CALL solveEq_r(eqMat_q,kappaVec_q)


  !  ! Convert the matrices to np*np matrix
  !  kappa_p_r = transpose(reshape(kappaVec_p(1:np**2),(/np,np/)) )
  !  kappa_p_c = transpose(reshape(kappaVec_p(1+np**2:2*np**2),(/np,np/)) )
  !  kappa_q_r = transpose(reshape(kappaVec_q(1:np**2),(/np,np/)) )
  !  kappa_q_c = transpose(reshape(kappaVec_q(1+np**2:2*np**2),(/np,np/)) )

  !  kappa_p = kappa_p_r + Ic * kappa_p_c
  !  kappa_q = kappa_q_r + Ic * kappa_q_c

  !  !-- Compute the pdots and qdots
  !  !-- pdot(i) = sum_j,k 0.5*M_pi(i,j)*ov(f(j),f(k))*p(k)*Kappa_p(k,j) + v_p(i)

  !  tmpV1_p=0._rl;tmpV2_p=0._rl
  !  tmpV1_q=0._rl;tmpV2_q=0._rl
  !  do j=1,np
  ! 	  tmpV1_p(j)=tmpV1_p(j) + 0.5_rl*sum( p(:)*st%ov_ff(j,:)*Kappa_p(:,j) )
  ! 	  tmpV1_q(j)=tmpV1_q(j) + 0.5_rl*sum( q(:)*st%ov_hh(j,:)*Kappa_q(:,j) )
  !  end do

  !  tmpV2_p = M_pi .matprod. tmpV1_p
  !  tmpV2_q = M_qi .matprod. tmpV1_q

  !  st%pdot = tmpV2_p + v_p
  !  st%qdot = tmpV2_q + v_q


  !  !-- Compute the fdots and hdots
  !  !-- fdot(i) = sum_j,k [ -N_pi(i,j)*p(k)*Kappa_p(k,j)*pc(j)*f(k)*ov(f(j),f(k)) 
  !  !								+ sum_l,m [ -N_pi(i,l)*A_p(m,j,k)*Kappa_p(k,j)*pc(l)f(m)*ov(f(l),f(m)) ] ] + v_f(i)
  !  tmpV1_f=0._cx
  !  tmpV1_h=0._cx

  !  do s=1, sys%nmode
  ! 	do i=1,np
  ! 	  do j=1,np
  ! 			tmpV1_f(i,s) = tmpV1_f(i,s) + sum( A_f(i,j,:,s) * Kappa_p(:,j) )
  ! 			tmpV1_h(i,s) = tmpV1_h(i,s) + sum( A_h(i,j,:,s) * Kappa_q(:,j) )
  ! 	  end do
  ! 	end do
  !  end do
  !  
  !  st%fdot = tmpV1_f + v_f
  !  st%hdot = tmpV1_h + v_h

  !END SUBROUTINE calcDerivatives_slow
  !SUBROUTINE calcDerivatives_symmetric(sys,st) 

  !  type(param), intent(in)                            		::  sys
  !  type(state), intent(in out)	                      		::  st
  !  complex(cx), dimension(st%np)							 		::  bigP, bigQ, v_p, v_q
  !  complex(cx), dimension(st%np,sys%nmode)				 		::  bigF, bigH, v_f, v_h
  !  complex(cx), dimension(st%np,st%np)              		::  M_p, M_q, M_pi, M_qi,N_p, N_q, N_pi, N_qi
  !  complex(cx), dimension(st%np,st%np,st%np)       		::  A_p, A_q
  !  complex(cx), dimension(st%np,st%np,st%np,sys%nmode)  ::  A_f, A_h
  !  real(rl), dimension(2*st%np**2,2*st%np**2)   	 		::  EqMat_p, EqMat_q
  !  real(rl), dimension(2*st%np**2)    		 			 		::  EqRhs_p, EqRhs_q, kappaVec_p, kappaVec_q
  !  real(rl), dimension(st%np,st%np)	   				 		::  Kappa_p_r,Kappa_q_r,kappa_p_c,kappa_q_c
  !  complex(cx), dimension(st%np,st%np)					 		::  Kappa_p,Kappa_q  !-- FINAL KAPPA MATRICES
  !  complex(cx),dimension(st%np,sys%nmode)  			 		::  f,h,fc,hc  		
  !  complex(cx),dimension(st%np)			   			 		::  p,q,pc,qc
  !  complex(cx),dimension(st%np)			   			 		::  der_E_pc_save,der_E_qc_save
  !  complex(cx),dimension(st%np,sys%nmode)			   	::  der_E_fc_save,der_E_hc_save
  !  integer												          		::  i,j,k,l,m,cj,ii,jj,np,info,s
  !  
  !  complex(cx), dimension(st%np) 				          		::  tmpV1_p,tmpV1_q,tmpV2_p,tmpV2_q
  !  complex(cx), dimension(st%np,st%np,st%np,sys%nmode) 	::  tmpA_f1,tmpA_h1,tmpA_f2,tmpA_h2
  !  complex(cx), dimension(st%np,sys%nmode) 				   ::  tmpV1_f,tmpV1_h
  !  complex(cx), dimension(st%np,st%np,st%np,st%np)		::  D_f,D_h
  !  complex(cx), dimension(st%np,st%np,st%np)				::  E_f,E_h
  !  complex(cx), dimension(st%np,st%np)						::  K_f,K_h


  !  !-- A few shortcuts
  !  f=st%f;h=st%h;p=st%p;q=st%q
  !  fc=conjg(f);hc=conjg(h);pc=conjg(p);qc=conjg(q)
  !  np = st%np

  !  !- set all variables to zero
  !  bigP = 0._cx;bigQ = 0._cx;bigF = 0._cx;bigH = 0._cx
  !  M_p=0._cx;M_q=0._cx;N_p=0._cx;N_q=0._cx
  !  M_pi=0._cx;M_qi=0._cx;N_pi=0._cx;N_qi=0._cx
  !  A_p=0._cx;A_q=0._cx;A_f=0._cx;A_h=0._cx
  !  eqMat_p=0._cx;eqMat_q=0._cx;eqrhs_p=0._cx;eqrhs_q=0._cx
  !	 eqRhs_p = 0._cx;eqRhs_q = 0._cx
  !  v_p=0._cx;v_q=0._cx;v_f=0._cx;v_h=0._cx
  !  Kappa_p_r=0._cx;Kappa_q_r=0._cx;kappa_p_c=0._cx;kappa_q_c=0._cx
  !  kappaVec_p=0._cx; kappaVec_q=0._cx;Kappa_p=0._cx;Kappa_q=0._cx
  !  der_E_fc_save(:,:) = 0._cx
  !  der_E_hc_save(:,:) = 0._cx
  !  der_E_pc_save(:) = 0._cx
  !  der_E_qc_save(:) = 0._cx
  !  D_f = 0._cx;D_h = 0._cx
  !  E_f = 0._cx;E_h = 0._cx
  !  K_f = 0._cx;K_h = 0._cx

  !  !==================================================
  !  !-- STEP 1: Computation of A_p,A_q and v_p,v_q
  !  !-- pDot(i) = sum_(j,k) [ A_p(i,j,k)*Kappa(k,j) + v_p(i) ] 
  !  !-- qDot(i) = sum_(j,k) [ A_q(i,j,k)*Kappa(k,j) + v_q(i) ] 

  !  do i=1,np
  ! 	  bigP(i) = P_j(sys,st,i)
  ! 	  bigQ(i) = Q_j(sys,st,i)
  ! 	  bigF(i,:) =  F_j(sys,st,i) 
  ! 	  bigH(i,:) =  H_j(sys,st,i) 
  !  end do

  !  !-- Matrix M: to be inverted
  !  M_p = st%ov_ff
  !  !M_q = st%ov_hh

  !  !-- perform the inversion
  !  M_pi = M_p
  !  !M_qi = M_q
  !  CALL invertH(M_pi,info)
  !  !CALL invertH(M_qi,info)

  !  !-- Define v_p, v_q and A_p, A_q
  !  v_p = M_pi .matprod. bigP
  !  !v_q = M_qi .matprod. bigQ

  !  do i=1,size(A_p,1)
  ! 	do j=1,size(A_p,2)
  ! 	  do k=1,size(A_p,3)
  ! 	    A_p(i,j,k) = 0.5_rl* M_pi(i,j) * st%ov_ff(j,k) * p(k)
  ! 		 !A_q(i,j,k) = 0.5_rl* M_qi(i,j) * st%ov_hh(j,k) * q(k)
  ! 	  end do
  ! 	end do
  !  end do
  !  

  !  !==================================================
  !  !-- STEP 2: Computation of A_f,A_h and v_f,v_h
  !  !-- fDot(i) = sum_(j,k) [ A_f(i,j,k)*Kappa(k,j) ] + v_f(i)  
  !  !-- hDot(i) = sum_(j,k) [ A_h(i,j,k)*Kappa(k,j) ] + v_h(i)  

  !  do j=1,np
  ! 	do m=1,np
  ! 	  N_p(j,m) = p(m)*st%ov_ff(j,m)
  ! 	  !N_q(j,m) = q(m)*st%ov_hh(j,m) 
  ! 	end do
  !  end do

  !  !-- Calculate Nij and is inverse
  !  N_pi = N_p
  !  !N_qi = N_q
  !  CALL invertGeneral(N_pi,info)
  !  !CALL invertGeneral(N_qi,info)
  !  

  !  !-- Computation of v_f and v_h (H_i in the notes)
  !  tmpV1_f=0._cx
  !  tmpV1_h=0._cx

  !  !-- First term of v_f (and v_h)
  !  tmpV1_f = bigF
  !  tmpV1_h = bigH
  !  !-- Second term of v_f
  !  do s=1,sys%nmode
  ! 	do j=1,np
  ! 		 tmpV1_f(j,s) = tmpV1_f(j,s) &
  ! 								 - sum( st%ov_ff(j,:)*v_p(:)*f(:,s) )
  ! 		! tmpV1_h(j,s) = tmpV1_h(j,s) &
  ! 		!						 - sum( st%ov_hh(j,:)*v_q(:)*h(:,s) )
  ! 	end do
  ! 	v_f(:,s) = N_pi .matprod. tmpV1_f(:,s)
  ! 	!v_h(:,s) = N_qi .matprod. tmpV1_h(:,s)
  !  end do

  !  tmpA_f1=0._cx; tmpA_h1=0._cx
  !  tmpA_f2=0._cx; tmpA_h2=0._cx

  !  !-- Defining the Betref matrices, here designated by A_f and A_h

  !  do s=1,sys%nmode
  ! 	do k=1,size(A_f,2)
  ! 	  do j=1,size(A_f,3)

  ! 		 do i=1,size(A_f,1)
  ! 			tmpA_f1(i,j,k,s) = 0.5_rl* N_pi(i,j) * p(k)*f(k,s) * st%ov_ff(j,k)
  ! 	!		tmpA_h1(i,j,k,s) = 0.5_rl* N_qi(i,j) * q(k)*h(k,s) * st%ov_hh(j,k)
  ! 		 end do

  ! 		 do l=1,np
  ! 			  tmpA_f2(l,j,k,s) = tmpA_f2(l,j,k,s) & 
  ! 						 - sum( A_p(:,j,k)*f(:,s)*st%ov_ff(l,:) )
  ! 	!		  tmpA_h2(l,j,k,s) = tmpA_h2(l,j,k,s) &
  ! 	!				 - sum( A_q(:,j,k)*h(:,s)*st%ov_hh(l,:) )
  ! 		 end do

  ! 	  end do
  ! 	

  ! 	  A_f(:,:,k,s) = tmpA_f1(:,:,k,s) + ( N_pi .matprod. tmpA_f2(:,:,k,s) )
  ! 	!  A_h(:,:,k,s) = tmpA_h1(:,:,k,s) + ( N_qi .matprod. tmpA_h2(:,:,k,s) )

  ! 	end do
  !  end do

  !  !==================================================
  !  !-- STEP 3: Solve the system of equations for Kappa
  !  !-- Define Q_ij, a 2*np**2 matrix: j in [1,np**2] correpsonds to Re(kappa)
  !  !--    while j in [1+np**2,2*np**2] corresponds to Im(Kappa)

  !  cj = st%np**2   !-- a shortcut to get to the imagainary part of kappaVec

  !  do i=1,st%np
  ! 	do j=1,st%np
  ! 	  K_f(i,j) = sum( conjg(f(i,:))*v_f(i,:) + f(i,:)*conjg(v_f(i,:)) - 2._rl*conjg(f(j,:))*v_f(i,:) )
  ! 	  !K_h(i,j) = sum( conjg(h(i,:))*v_h(i,:) + h(i,:)*conjg(v_h(i,:)) - 2._rl*conjg(h(j,:))*v_h(i,:) )
  ! 	  do k=1,st%np
  ! 		 E_f(i,j,k) = sum( f(i,:)*conjg(A_f(i,j,k,:))  )
  ! 		 !E_h(i,j,k) = sum( h(i,:)*conjg(A_h(i,j,k,:))  )
  ! 		 do m=1,st%np
  ! 			D_f(i,m,j,k) = sum( (conjg(f(i,:)) - 2._rl*conjg(f(m,:)))*A_f(i,j,k,:) )
  ! 			!D_h(i,m,j,k) = sum( (conjg(h(i,:)) - 2._rl*conjg(h(m,:)))*A_h(i,j,k,:) )
  ! 		 end do
  ! 	  end do
  ! 	end do
  !  end do
  !  

  !  do i=1, np
  ! 	do m=1, np

  ! 	  ii = np*(i-1)+m

  ! 	  eqRhs_p(ii) 	  =  real( K_f(i,m) )
  ! 	  !eqRhs_q(ii) 	  =  real( K_h(i,m) )
  ! 	  eqRhs_p(ii+cj)  =  aimag( K_f(i,m) )
  ! 	  !eqRhs_q(ii+cj)  =  aimag( K_h(i,m) )

  ! 	  do j=1, np
  ! 		 do k=1, np

  ! 		   jj = np*(k-1) + j

  ! 				 eqMat_p(ii,jj) = kroneckerDelta(ii,jj) - real(  D_f(i,m,j,k) + E_f(i,j,k) )
  ! 				 eqMat_p(ii,jj+cj) = kroneckerDelta(ii,jj+cj) + aimag(  D_f(i,m,j,k) - E_f(i,j,k) )
  ! 				 eqMat_p(ii+cj,jj) = kroneckerDelta(ii+cj,jj) - aimag(  D_f(i,m,j,k) + E_f(i,j,k) )
  ! 				 eqMat_p(ii+cj,jj+cj) = kroneckerDelta(ii+cj,jj+cj) - real(  D_f(i,m,j,k) - E_f(i,j,k) )

  ! 		!		 eqMat_q(ii,jj) = kroneckerDelta(ii,jj) - real(  D_h(i,m,j,k) + E_h(i,j,k) )
  ! 		!		 eqMat_q(ii,jj+cj) = kroneckerDelta(ii,jj+cj) + aimag(  D_h(i,m,j,k) - E_h(i,j,k) )
  ! 		!		 eqMat_q(ii+cj,jj) = kroneckerDelta(ii+cj,jj) - aimag(  D_h(i,m,j,k) + E_h(i,j,k) )
  ! 		!		 eqMat_q(ii+cj,jj+cj) = kroneckerDelta(ii+cj,jj+cj) - real(  D_h(i,m,j,k) - E_h(i,j,k) )

  ! 		 end do
  ! 	  end do
  ! 	end do
  !  end do


  !  !-- Apply the Lapack algorithm to solve the real system of equations
  !  !-- the solution replaces the 2nd argument
  !  kappaVec_p = eqRhs_p
  !  CALL solveEq_r(eqMat_p,kappaVec_p)
  !  !kappaVec_q = eqRhs_q
  !  !CALL solveEq_r(eqMat_q,kappaVec_q)


  !  ! Convert the matrices to np*np matrix
  !  kappa_p_r = transpose(reshape(kappaVec_p(1:np**2),(/np,np/)) )
  !  kappa_p_c = transpose(reshape(kappaVec_p(1+np**2:2*np**2),(/np,np/)) )
  !  !kappa_q_r = transpose(reshape(kappaVec_q(1:np**2),(/np,np/)) )
  !  !kappa_q_c = transpose(reshape(kappaVec_q(1+np**2:2*np**2),(/np,np/)) )

  !  kappa_p = kappa_p_r + Ic * kappa_p_c
  !  !kappa_q = kappa_q_r + Ic * kappa_q_c

  !  !-- Compute the pdots and qdots
  !  !-- pdot(i) = sum_j,k 0.5*M_pi(i,j)*ov(f(j),f(k))*p(k)*Kappa_p(k,j) + v_p(i)

  !  tmpV1_p=0._rl;tmpV2_p=0._rl
  !  !tmpV1_q=0_rl;tmpV2_q=0_rl
  !  do j=1,np
  ! 	  tmpV1_p(j)=tmpV1_p(j) + 0.5_rl*sum( p(:)*st%ov_ff(j,:)*Kappa_p(:,j) )
  ! !	  tmpV1_q(j)=tmpV1_q(j) + 0.5_rl*sum( q(:)*st%ov_hh(j,:)*Kappa_q(:,j) )
  !  end do

  !  tmpV2_p = M_pi .matprod. tmpV1_p
  !  !tmpV2_q = M_qi .matprod. tmpV1_q

  !  st%pdot = tmpV2_p + v_p
  !  !st%qdot = tmpV2_q + v_q
  !  st%qdot = st%pdot


  !  !-- Compute the fdots and hdots
  !  !-- fdot(i) = sum_j,k [ -N_pi(i,j)*p(k)*Kappa_p(k,j)*pc(j)*f(k)*ov(f(j),f(k)) 
  !  !								+ sum_l,m [ -N_pi(i,l)*A_p(m,j,k)*Kappa_p(k,j)*pc(l)f(m)*ov(f(l),f(m)) ] ] + v_f(i)
  !  tmpV1_f=0._cx
  !  !tmpV1_h=0_cx

  !  do s=1, sys%nmode
  ! 	do i=1,np
  ! 	  do j=1,np
  ! 			tmpV1_f(i,s) = tmpV1_f(i,s) + sum( A_f(i,j,:,s) * Kappa_p(:,j) )
  ! !			tmpV1_h(i,s) = tmpV1_h(i,s) + sum( A_h(i,j,:,s) * Kappa_q(:,j) )
  ! 	  end do
  ! 	end do
  !  end do
  !  
  !  st%fdot = tmpV1_f + v_f
  !  !st%hdot = tmpV1_h + v_h
  !  st%hdot = - st%fdot

  !END SUBROUTINE

  !FUNCTION f_nx_eo_lowk(sys,st)

  !  type(param),intent(in)   	::  sys
  !  type(state),intent(in)    ::  st
  !  real(rl)				 	   ::  x
  !  integer 						::  n,i,kmin,kmax
  !  complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode)  ::  f_nx_eo_lowk

  !  kmin = int( (-sys%nmode+1)/10._rl )
  !  kmax = int( sys%nmode/10._rl )

  !  do n=1,st%np
  ! 	do i=-sys%nmode+1,sys%nmode
  ! 	  x = sys%dx * (i-0.5_rl)
  ! 	  f_nx_eo_lowk(n,i) = sqrt(0.5_rl/(2._rl*pi)) * sqrt(sys%dx) &
  ! 				 * SUM( sqrt(sys%dk(1)) * (st%f(n,kmin:kmax)+st%fo(n,kmin:kmax)) * exp( Ic*sys%w(kmin:kmax)*x ) &
  ! 							+ sqrt(sys%dk(1)) * (st%f(n,kmin:kmax)-st%fo(n,kmin:kmax)) * exp( -Ic*sys%w(kmin:kmax)*x ) )
  ! 	end do
  !  end do

  !END FUNCTION
!  FUNCTION transmission_n(sys,st,initial_st)
!
!	 type(param),intent(in) 					 ::  sys
!	 type(state),intent(in)						 ::  st, initial_st
!	 type(state)									 ::  upst,upst_ini
!	 real(rl)										 ::  transmission_n
!	 real(rl)										 ::  t_output,input
!	 complex(cx),dimension(st%np,sys%nmode) ::  zpk,zpk_ini
!	 complex(cx)					  				 ::  sum_k, tmp
!	 integer							 				 ::  n,m
!
!	 upst = st
!	 upst_ini = initial_st !	 upst%q(:) = 0._cx
!	 upst_ini%q(:) = 0._cx
!	 CALL normalise(upst)
!	 CALL normalise(upst_ini)
!
!	 zpk_ini = sqrt(0.5_rl) * ( upst_ini%f + upst_ini%fo ) 
!	 zpk = sqrt(0.5_rl) * ( upst%f + upst%fo )
!
!	 tmp = 0._cx
!	 do n=1, upst_ini%np
!		do m=1, upst_ini%np
!		  sum_k = sum( conjg(zpk_ini(n,:))*zpk_ini(m,:) )
!		  tmp = tmp + sum_k*conjg(upst%p(n))*upst%p(m)*upst%ov_ff(n,m)
!		end do
!	 end do
!	 input = real(tmp)
!
!	 tmp = 0._cx
!	 do n=1, upst%np
!		do m=1,upst%np
!		  sum_k = sum( conjg(zpk(n,:))*zpk(m,:) ) 
!		  tmp = tmp + sum_k*conjg(upst%p(n))*upst%p(m)*upst%ov_ff(n,m)
!		end do
!	 end do
!	 t_output = real(tmp)
!
!	 transmission_n = t_output/input
!
!  END FUNCTION
!  FUNCTION reflection_n(sys,st,initial_st)
!
!	 type(param),intent(in) 					 ::  sys
!	 type(state),intent(in)						 ::  st, initial_st
!	 type(state)									 ::  upst,upst_ini
!	 real(rl)										 ::  reflection_n
!	 real(rl)										 ::  r_output,input
!	 complex(cx),dimension(st%np,sys%nmode) ::  zmk,zmk_ini
!	 complex(cx)					  				 ::  sum_k, tmp
!	 integer							 				 ::  n,m
!
!	 upst = st
!	 upst_ini = initial_st
!	 upst%q(:) = 0._cx
!	 upst_ini%q(:) = 0._cx
!	 CALL normalise(upst)
!	 CALL normalise(upst_ini)
!
!	 zmk_ini = sqrt(0.5_rl) * ( upst_ini%f + upst_ini%fo ) 
!	 zmk = sqrt(0.5_rl) * ( upst%f - upst%fo )
!
!	 tmp = 0._cx
!	 do n=1, upst_ini%np
!		do m=1, upst_ini%np
!		  sum_k = sum( conjg(zmk_ini(n,:))*zmk_ini(m,:) )
!		  tmp = tmp + sum_k*conjg(upst%p(n))*upst%p(m)*upst%ov_ff(n,m)
!		end do
!	 end do
!	 input = real(tmp)
!
!	 tmp = 0._cx
!	 do n=1, upst%np
!		do m=1,upst%np
!		  sum_k = sum( conjg(zmk(n,:))*zmk(m,:) ) 
!		  tmp = tmp + sum_k*conjg(upst%p(n))*upst%p(m)*upst%ov_ff(n,m)
!		end do
!	 end do
!	 r_output = real(tmp)
!
!	 reflection_n = r_output/input
!
!  END FUNCTION
!  SUBROUTINE reflection_n_routine(sys,st,initial_st)
!
!	 type(param),intent(in) 					 ::  sys
!	 type(state),intent(in)						 ::  st, initial_st
!	 type(state)									 ::  upst,upst_ini
!	 real(rl)										 ::  r_output,input
!	 complex(cx),dimension(st%np,sys%nmode) ::  zmk,zmk_ini
!	 complex(cx)					  				 ::  sum_k, tmp
!	 integer							 				 ::  n,m
!
!	 print*,"0"
!	 upst = st
!	 print*,"01"
!	 upst_ini = initial_st
!	 print*,"02"
!	 upst%q(:) = 0._cx
!	 print*,"03"
!	 upst_ini%q(:) = 0._cx
!	 print*,"1"
!	 CALL normalise(upst)
!	 CALL normalise(upst_ini)
!
!	 print*,"2"
!	 zmk_ini = sqrt(0.5_rl) * ( upst_ini%f + upst_ini%fo ) 
!	 zmk = sqrt(0.5_rl) * ( upst%f - upst%fo )
!
!	 print*,"3"
!	 tmp = 0._cx
!	 do n=1, upst_ini%np
!		do m=1, upst_ini%np
!		  sum_k = sum( conjg(zmk_ini(n,:))*zmk_ini(m,:) )
!		  tmp = tmp + sum_k*conjg(upst%p(n))*upst%p(m)*upst%ov_ff(n,m)
!		end do
!	 end do
!	 input = real(tmp)
!
!	 print*,"4"
!	 tmp = 0._cx
!	 do n=1, upst%np
!		do m=1,upst%np
!		  sum_k = sum( conjg(zmk(n,:))*zmk(m,:) ) 
!		  tmp = tmp + sum_k*conjg(upst%p(n))*upst%p(m)*upst%ov_ff(n,m)
!		end do
!	 end do
!	 r_output = real(tmp)
!	 print*,"5"
!
!	 !reflection_n = r_output/input
!
!	 END SUBROUTINE
!  FUNCTION fo_nx(sys,st)
!
!	 type(param),intent(in)   						::  sys
!	 type(state),intent(in)    					::  st
!	 real(rl)				 	   					::  x
!	 complex(cx), dimension(st%np,sys%nmode)  ::  fo_nx
!	 integer												::  n,i
!
!	 do n=1,st%np
!		do i=1,sys%nmode
!		  x = sys%dx * (i-1)
!		  fo_nx(n,i) = ( 2.0/sqrt(2._rl*pi) )*sqrt(sys%dx)*sqrt(sys%dk) &
!		  * SUM( st%fo(n,:) * cos( sys%w(:)*x )  )
!		end do
!		fo_nx(n,1) = fo_nx(n,1)/sqrt(2._rl)
!	 end do
!
!  END FUNCTION
!  FUNCTION ho_nx(sys,st)
!
!	 type(param),intent(in)   						::  sys
!	 type(state),intent(in)    					::  st
!	 real(rl)				 	   					::  x
!	 complex(cx), dimension(st%np,sys%nmode)  ::  ho_nx
!	 integer												::  n,i
!
!	 do n=1,st%np
!		do i=1,sys%nmode
!		  x = sys%dx * (i-1)
!		  ho_nx(n,i) = ( 2.0/sqrt(2._rl*pi) )*sqrt(sys%dx)*sqrt(sys%dk) &
!		  * SUM( st%ho(n,:) * cos( sys%w(:)*x )  )
!		end do
!		ho_nx(n,1) = ho_nx(n,1)/sqrt(2._rl)
!	 end do
!
!  END FUNCTION

!  FUNCTION fks(sys,st,xmin,xmax,xmin2,xmax2)
!
!	 type(param),intent(in)				::  sys
!	 type(state),intent(in)				::  st
!	 real(rl), intent(in)				::  xmin, xmax, xmin2, xmax2
!	 complex(cx), dimension(st%np,-sys%nmode+1:sys%nmode-1)  ::  fnx
!	 complex(cx), dimension(st%np,sys%nmode)  ::  fks
!	 real(rl)							  	::  xx(-sys%nmode+1:sys%nmode-1)
!	 integer									::  k,n,i
!
!	 fnx(:,:int(xmin/sys%dx)) = 0._rl
!	 fnx(:,int(xmax/sys%dx):int(xmin2/sys%dx)) = 0._rl
!	 fnx(:,int(xmax2/sys%dx):) = 0._rl
!
!	 do i=-sys%nmode+1,sys%nmode-1
!	 	xx(i) = sys%dx * i
!	 end do
!
!	 do n=1,st%np
!		do k=1,sys%nmode
!		  fks(n,k) = sys%dx*sqrt(1._rl/(2._rl*pi))*sum( fnx(n,:)*exp( -Ic*sys%w(k)*xx(:) ) )
!		end do
!	 end do
!	 
!
!  END FUNCTION
!
!!
!  FUNCTION zx(sys, sigma, k0 , x0, nb) result(zkout)
!
!!--  this function calculates: zx^e = z_x + z_-x
!
!	 type(param),intent(in) 	::  sys
!	 real(rl), intent(in) 		::  sigma, k0, x0, nb
!	 complex(cx)          		::  zkout(-sys%nmode+1:sys%nmode-1)
!
!	 complex(cx)   ::  prefac
!	 real(rl)      ::  x
!	 integer       ::  i
!
!	 prefac = 0._rl
!	 prefac  = sqrt(nb)*((2._rl*sigma**2/pi)**0.25_rl)*exp(+Ic*0.5_rl*k0*x0)
!	 zkout = 0._rl
!
!	 do i= -sys%nmode+1, sys%nmode-1
!
!		x = dble(i)*sys%dx
!		zkout(i) = prefac * ( exp(Ic*(x-x0)*k0) * exp( - sigma**2 * (x-x0)**2 ) )
!
!	 end do
!
!  END FUNCTION
!












!  FUNCTION aDag2_up(sys,st,x)
!
!	 type(param), intent(in)  :: sys
!	 type(state), intent(in)  ::  st
!	 real(rl), intent(in)	  ::  x
!	 complex(cx)				  ::  aDag2_up
!	 complex(cx)				  ::  tmp
!	 type(state)				  ::  upst
!	 integer 					  ::  n,m
!	 
!
!	 upst = st
!	 upst%q(:) = 0._cx
!	 CALL normalise(upst)
!
!	 tmp = 0._cx
!	 do n=1,size(st%p,1)
!		do m=1,size(st%p,1)
!		  tmp = tmp + conjg(upst%p(m))*upst%p(n) * conjg(fn_x(sys,upst,m,x))**2 * upst%ov_ff(m,n)
!		end do	
!	 end do
!
!	 aDag2_up = tmp





!  SUBROUTINE isolateGroundState(sys,st)
!
!	 type(param), intent(in)						::  sys
!	 type(state), intent(in out)					::  st
!	 type(state)										::  gs
!	 integer												::  i,l,j
!	 real(rl)											::  x,k
!	 real(rl)											::  lengthToKeep
!	 integer												::  sitesToKeep
!
!	 lengthToKeep = 3._rl
!	 sitesToKeep = int(lengthToKeep/sys%dx)+1
!
!	 gs=st
!	 gs%f(:,:) = 0._cx
!	 gs%h(:,:) = 0._cx
!
!	 do i=1, st%np
!		siteDo2: do l=1, sitesToKeep
!		  x = (sys%length/sys%nmode) * (l-1)
!		  if ( real(fn_x(sys,st,i,x )) < 0 ) then
!			 if (sitesToKeep > l-1 ) then
!				sitesToKeep = l-1
!			 end if
!			 exit siteDo2
!		  end if
!		end do siteDo2
!	 end do
!
!	 sitestokeep = int(sys%nmode*(0.25_rl))
!	 do i=1, st%np
!		do l=1, sitestokeep
!		  x = (sys%length/sys%nmode) * (l-1)
!		  gs%f(i,:) = gs%f(i,:) + fn_x(sys,st,i,x)*exp( - Ic*sys%w(:)*2_rl*PI*x ) * sqrt(1/(sys%nmode*sys%wmax))
!		  gs%h(i,:) = gs%h(i,:) + hn_x(sys,st,i,x)*exp( - Ic*sys%w(:)*2_rl*PI*x ) * sqrt(1/(sys%nmode*sys%wmax)) 
!		end do
!	 end do
!
!	 CALL update_sums(sys,gs)
!	 CALL normalise(gs)
!
!	 st = gs
!
!  END SUBROUTINE
  !-- VARIANTS FOR ROUTINES
  !-- a routine to caluclate the derivatives using manual summations (no SUM)
!
! SUBROUTINE calcDerivatives_manualsum(sys,st) 
!
!	 type(param), intent(in)                            		::  sys
!	 type(state), intent(in out)	                      		::  st
!	 complex(cx), dimension(st%np)							 		::  bigP, bigQ, v_p, v_q
!	 complex(cx), dimension(st%np,sys%nmode)				 		::  bigF, bigH, v_f, v_h
!	 complex(cx), dimension(st%np,st%np)              		::  M_p, M_q, M_pi, M_qi,N_p, N_q, N_pi, N_qi
!	 complex(cx), dimension(st%np,st%np,st%np)       		::  A_p, A_q
!	 complex(cx), dimension(st%np,st%np,st%np,sys%nmode)  ::  A_f, A_h
!	 real(rl), dimension(2*st%np**2,2*st%np**2)   	 		::  EqMat_p, EqMat_q
!	 real(rl), dimension(2*st%np**2)    		 			 		::  EqRhs_p, EqRhs_q, kappaVec_p, kappaVec_q
!	 real(rl), dimension(st%np,st%np)	   				 		::  Kappa_p_r,Kappa_q_r,kappa_p_c,kappa_q_c
!	 complex(cx), dimension(st%np,st%np)					 		::  Kappa_p,Kappa_q  !-- FINAL KAPPA MATRICES
!	 complex(cx),dimension(st%np,sys%nmode)  			 		::  f,h,fc,hc  		
!	 complex(cx),dimension(st%np)			   			 		::  p,q,pc,qc
!	 complex(cx),dimension(st%np)			   			 		::  der_E_pc_save,der_E_qc_save
!	 integer												          		::  i,j,k,l,m,cj,ii,jj,np,info,s
!	 
!	 complex(cx), dimension(st%np) 				          		::  tmpV1_p,tmpV1_q,tmpV2_p,tmpV2_q
!	 complex(cx), dimension(st%np,st%np,st%np,sys%nmode) 	::  tmpA_f1,tmpA_h1,tmpA_f2,tmpA_h2
!	 complex(cx), dimension(st%np,sys%nmode) 				   ::  tmpV1_f,tmpV1_h
!
!
!	 !-- A few shortcuts
!	 f=st%f;h=st%h;p=st%p;q=st%q
!	 fc=conjg(f);hc=conjg(h);pc=conjg(p);qc=conjg(q)
!	 np = st%np
!
!	 !- set all variables to zero
!	 bigP = 0;bigQ = 0;bigF = 0;bigH = 0
!	 M_p=0;M_q=0;N_p=0;N_q=0
!	 M_pi=0;M_qi=0;N_pi=0;N_qi=0
!	 A_p=0;A_q=0;A_f=0;A_h=0
!	 eqMat_p=0;eqMat_q=0;eqrhs_p=0;eqrhs_q=0
!	 v_p=0;v_q=0;v_f=0;v_h=0
!	 Kappa_p_r=0;Kappa_q_r=0;kappa_p_c=0;kappa_q_c=0
!	 kappaVec_p=0; kappaVec_q=0;Kappa_p=0;Kappa_q=0
!
!
!	 !==================================================
!	 !-- STEP 1: Computation of A_p,A_q and v_p,v_q
!	 !-- pDot(i) = sum_(j,k) [ A_p(i,j,k)*Kappa(k,j) + v_p(i) ] 
!	 !-- qDot(i) = sum_(j,k) [ A_q(i,j,k)*Kappa(k,j) + v_q(i) ] 
!
!	 do i=1,np
!		  der_E_pc_save(i) = der_E_pc(sys,st,i)
!		  der_E_qc_save(i) = der_E_qc(sys,st,i)
!		  bigP(i) = - Ic * der_E_pc_save(i)
!		  bigQ(i) = - Ic * der_E_qc_save(i)
!		  do s=1,sys%nmode
!			 bigF(i,s) = - Ic * ( der_E_fc(sys,st,i,s) &
!							 + 0.5_rl * f(i,s) * ( der_E_pc_save(i)*pc(i) + conjg(der_E_pc_save(i))*p(i) ) )
!			 bigH(i,s) = - Ic * ( der_E_hc(sys,st,i,s) &
!							 + 0.5_rl * h(i,s) * ( der_E_qc_save(i)*qc(i) + conjg(der_E_qc_save(i))*q(i) ) )
!		  end do
!	 end do
!
!	 !-- Matrix M: to be inverted
!	 M_p = st%ov_ff
!	 M_q = st%ov_hh
!
!	 !-- perform the inversion
!	 M_pi = M_p
!	 M_qi = M_q
!	 CALL invertH(M_pi,info)
!	 CALL invertH(M_qi,info)
!
!	 
!	 !-- Define v_p, v_q and A_p, A_q
!	 v_p = M_pi .matprod. bigP
!	 v_q = M_qi .matprod. bigQ
!
!	 do i=1,size(A_p,1)
!		do j=1,size(A_p,2)
!		  do k=1,size(A_p,3)
!		    A_p(i,j,k) = 0.5_rl* M_pi(i,j) * st%ov_ff(j,k) * p(k)
!		    A_q(i,j,k) = 0.5_rl* M_qi(i,j) * st%ov_hh(j,k) * q(k)
!		  end do
!		end do
!	 end do
!	 
!
!	 !==================================================
!	 !-- STEP 2: Computation of A_f,A_h and v_f,v_h
!	 !-- fDot(i) = sum_(j,k) [ A_f(i,j,k)*Kappa(k,j) ] + v_f(i)  
!	 !-- hDot(i) = sum_(j,k) [ A_h(i,j,k)*Kappa(k,j) ] + v_h(i)  
!
!	 do j=1,np
!		do m=1,np
!		  N_p(j,m) = pc(j)*p(m)*st%ov_ff(j,m)
!		  N_q(j,m) = qc(j)*q(m)*st%ov_hh(j,m) 
!		end do
!	 end do
!
!	 !-- Calculate Nij and is inverse
!	 N_pi = N_p
!	 N_qi = N_q
!	 CALL invertH(N_pi,info)
!	 CALL invertH(N_qi,info)
!	 
!
!	 !-- Computation of v_f and v_h (H_i in the notes)
!	 tmpV1_f=0_cx
!	 tmpV1_h=0_cx
!
!	 !-- First term of v_f (and v_h)
!	 tmpV1_f = bigF
!	 tmpV1_h = bigH
!	 !-- Second term of v_f
!	 do s=1,sys%nmode
!		do j=1,np
!		  do m=1,np
!			 tmpV1_f(j,s) = tmpV1_f(j,s) &
!									 - pc(j)*st%ov_ff(j,m)*v_p(m)*f(m,s)
!			 tmpV1_h(j,s) = tmpV1_h(j,s) &
!									 - qc(j)*st%ov_hh(j,m)*v_q(m)*h(m,s)
!		  end do
!		end do
!		v_f(:,s) = N_pi .matprod. tmpV1_f(:,s)
!		v_h(:,s) = N_qi .matprod. tmpV1_h(:,s)
!	 end do
!
!	 tmpA_f1=0_cx; tmpA_h1=0_cx
!	 tmpA_f2=0_cx; tmpA_h2=0_cx
!
!	 !-- Defining the Betref matrices, here designated by A_f and A_h
!
!	 do s=1,sys%nmode
!		do k=1,size(A_f,2)
!		  do j=1,size(A_f,3)
!
!			 do i=1,size(A_f,1)
!				tmpA_f1(i,j,k,s) = 0.5_rl* N_pi(i,j) * p(k)*pc(j)*f(k,s) * st%ov_ff(j,k)
!				tmpA_h1(i,j,k,s) = 0.5_rl* N_qi(i,j) * q(k)*qc(j)*h(k,s) * st%ov_hh(j,k)
!			 end do
!
!			 do l=1,np
!				do m=1,np
!				  tmpA_f2(l,j,k,s) = tmpA_f2(l,j,k,s) &
!											 - A_p(m,j,k)*f(m,s)*pc(l)*st%ov_ff(l,m)
!				  tmpA_h2(l,j,k,s) = tmpA_h2(l,j,k,s) &
!											 - A_q(m,j,k)*h(m,s)*qc(l)*st%ov_hh(l,m)
!				end do
!			 end do
!
!		  end do
!		
!
!		  A_f(:,:,k,s) = tmpA_f1(:,:,k,s) + ( N_pi .matprod. tmpA_f2(:,:,k,s) )
!		  A_h(:,:,k,s) = tmpA_h1(:,:,k,s) + ( N_qi .matprod. tmpA_h2(:,:,k,s) )
!
!		end do
!	 end do
!
!	 !==================================================
!	 !-- STEP 3: Solve the system of equations for Kappa
!	 !-- Define Q_ij, a 2*np**2 matrix: j in [1,np**2] correpsonds to Re(kappa)
!	 !--    while j in [1+np**2,2*np**2] corresponds to Im(Kappa)
!
!	 cj = st%np**2   !-- a shortcut to get to the imagainary part of kappaVec
!
!	 do i=1, np
!		do m=1, np
!
!		  ii = np*(i-1)+m
!
!		  do j=1, np
!			 do k=1, np
!
!			   jj = np*(k-1) + j
!
!
!				  eqMat_p(ii,jj) =  kroneckerDelta(ii,jj) 
!				  eqMat_p(ii,jj+cj) = kroneckerDelta(ii,jj+cj) 
!				  eqMat_p(ii+cj,jj) = kroneckerDelta(ii+cj,jj) 
!				  eqMat_p(ii+cj,jj+cj) = kroneckerDelta(ii+cj,jj+cj) 
!
!				  eqMat_q(ii,jj) = kroneckerDelta(ii,jj) 
!				  eqMat_q(ii,jj+cj) = kroneckerDelta(ii,jj+cj) 
!				  eqMat_q(ii+cj,jj) = kroneckerDelta(ii+cj,jj) 
!				  eqMat_q(ii+cj,jj+cj) = kroneckerDelta(ii+cj,jj+cj) 
!				  
!				  do s=1,sys%nmode    !-- sum over the modes
!
!					 eqMat_p(ii,jj) =  eqMat_p(ii,jj) &
!												- real( ( conjg(f(i,s))-2_rl*conjg(f(m,s)) ) * A_f(i,j,k,s) ) &
!												- real( f(i,s) * conjg(A_f(i,j,k,s)) )
!					 eqMat_p(ii,jj+cj) = eqMat_p(ii,jj+cj) &
!											  + aimag( ( conjg(f(i,s))-2_rl*conjg(f(m,s)) ) * A_f(i,j,k,s) ) &
!											  - aimag( f(i,s)*conjg(A_f(i,j,k,s)) )
!					 eqMat_p(ii+cj,jj) = eqMat_p(ii+cj,jj) &
!												- aimag( ( conjg(f(i,s))-2_rl*conjg(f(m,s)) ) * A_f(i,j,k,s) ) &
!												- aimag( f(i,s) * conjg(A_f(i,j,k,s)) )
!					 eqMat_p(ii+cj,jj+cj) = eqMat_p(ii+cj,jj+cj) &
!												- real( ( conjg(f(i,s))-2_rl*conjg(f(m,s)) ) * A_f(i,j,k,s) ) &
!												+ real( f(i,s) * conjg(A_f(i,j,k,s)) )
!
!					 eqMat_q(ii,jj) = eqMat_q(ii,jj) &
!												- real( ( conjg(h(i,s))-2_rl*conjg(h(m,s)) ) * A_h(i,j,k,s) ) &
!												- real( h(i,s)*conjg(A_h(i,j,k,s)) )
!					 eqMat_q(ii,jj+cj) = eqMat_q(ii,jj+cj) &
!											  + aimag( ( conjg(h(i,s))-2_rl*conjg(h(m,s)) ) * A_h(i,j,k,s) ) &
!											  - aimag( h(i,s)*conjg(A_h(i,j,k,s)) )
!					 eqMat_q(ii+cj,jj) = eqMat_q(ii+cj,jj) &
!												- aimag( ( conjg(h(i,s))-2_rl*conjg(h(m,s)) ) * A_h(i,j,k,s) ) &
!												- aimag( h(i,s) * conjg(A_h(i,j,k,s)) )
!					 eqMat_q(ii+cj,jj+cj) = eqMat_q(ii+cj,jj+cj) &
!												- real( ( conjg(h(i,s))-2_rl*conjg(h(m,s)) ) * A_h(i,j,k,s) ) &
!												+ real( h(i,s) * conjg(A_h(i,j,k,s)) )
!				  end do
!
!			 end do
!		  end do
!		end do
!	 end do
!
!
!  	 !-- define the RHS of the equations for kappa
!  	 eqRhs_p = 0_cx
!  	 eqRhs_q = 0_cx
!	 do i=1,np
!		do m=1,np
!		  ii = np*(i-1) + m
!
!		  do s=1,sys%nmode
!			 eqRhs_p(ii) 	  =  eqRhs_p(ii) +  real( conjg(f(i,s))*v_f(i,s) + f(i,s)*conjg(v_f(i,s)) - 2_rl*conjg(f(m,s))*v_f(i,s) )
!			 eqRhs_q(ii) 	  =  eqRhs_q(ii) +  real( conjg(h(i,s))*v_h(i,s) + h(i,s)*conjg(v_h(i,s)) - 2_rl*conjg(h(m,s))*v_h(i,s) )
!			 eqRhs_p(ii+cj)  =  eqRhs_p(ii+cj) + aimag( - 2_rl*conjg(f(m,s))*v_f(i,s) )
!			 eqRhs_q(ii+cj)  =  eqRhs_q(ii+cj) + aimag( - 2_rl*conjg(h(m,s))*v_h(i,s) )
!		  end do
!
!		end do
!	 end do
!
!	 !-- Apply the Lapack algorithm to solve the real system of equations
!	 !-- the solution replaces the 2nd argument
!	 kappaVec_p = eqRhs_p
!	 CALL solveEq_r(eqMat_p,kappaVec_p)
!	 kappaVec_q = eqRhs_q
!	 CALL solveEq_r(eqMat_q,kappaVec_q)
!
!
!	 ! Convert the matrices to np*np matrix
!	 kappa_p_r = transpose(reshape(kappaVec_p(1:np**2),(/np,np/)) )
!	 kappa_p_c = transpose(reshape(kappaVec_p(1+np**2:2*np**2),(/np,np/)) )
!	 kappa_q_r = transpose(reshape(kappaVec_q(1:np**2),(/np,np/)) )
!	 kappa_q_c = transpose(reshape(kappaVec_q(1+np**2:2*np**2),(/np,np/)) )
!
!	 kappa_p = kappa_p_r + Ic * kappa_p_c
!	 kappa_q = kappa_q_r + Ic * kappa_q_c
!
!	 !-- Compute the pdots and qdots
!	 !-- pdot(i) = sum_j,k 0.5*M_pi(i,j)*ov(f(j),f(k))*p(k)*Kappa_p(k,j) + v_p(i)
!
!	 tmpV1_p=0_rl;tmpV2_p=0_rl
!	 tmpV1_q=0_rl;tmpV2_q=0_rl
!	 do j=1,np
!		do k=1,np
!		  tmpV1_p(j)=tmpV1_p(j) &
!						+ 0.5_rl*p(k)*st%ov_ff(j,k)*Kappa_p(k,j)
!		  tmpV1_q(j)=tmpV1_q(j) &
!						+ 0.5_rl*q(k)*st%ov_hh(j,k)*Kappa_q(k,j)
!		end do
!	 end do
!
!	 tmpV2_p = M_pi .matprod. tmpV1_p
!	 tmpV2_q = M_qi .matprod. tmpV1_q
!
!	 st%pdot = tmpV2_p + v_p
!	 st%qdot = tmpV2_q + v_q
!
!
!	 !-- Compute the fdots and hdots
!	 !-- fdot(i) = sum_j,k [ -N_pi(i,j)*p(k)*Kappa_p(k,j)*pc(j)*f(k)*ov(f(j),f(k)) 
!	 !								+ sum_l,m [ -N_pi(i,l)*A_p(m,j,k)*Kappa_p(k,j)*pc(l)f(m)*ov(f(l),f(m)) ] ] + v_f(i)
!	 tmpV1_f=0_cx
!	 tmpV1_h=0_cx
!
!	 do s=1, sys%nmode
!		do i=1,np
!		  do j=1,np
!			 do k=1,np
!				tmpV1_f(i,s) = tmpV1_f(i,s) + A_f(i,j,k,s) * Kappa_p(k,j)
!				tmpV1_h(i,s) = tmpV1_h(i,s) + A_h(i,j,k,s) * Kappa_q(k,j)
!			 end do
!		  end do
!		end do
!	 end do
!	 
!	 st%fdot = tmpV1_f + v_f
!	 st%hdot = tmpV1_h + v_h
!
!  END SUBROUTINE calcDerivatives_manualsum
!
!  FUNCTION rehn_x(sys,st,n,x)
!
!	 type(param),intent(in)   	::  sys
!	 type(state),intent(in)    ::  st
!	 real(rl), intent(in) 	   ::  x
!	 integer, intent(in)       ::  n
!	 real(rl)	 					::  rehn_x
!	 complex(cx)					::  tmp
!	 integer							::  s
!
!	 rehn_x = 0_rl
!	 tmp = 0_cx
!	 do s=1, sys%nmode
!		tmp = tmp + real( st%h(n,s) )* exp( Ic*(sys%w(s)/c_light)*2_rl*PI*x ) &
!														  * sqrt(1._rl/sys%nmode)
!	 end do
!	 rehn_x = real(tmp)
!
!  END FUNCTION 
!  FUNCTION imfn_x(sys,st,n,x)
!
!	 type(param),intent(in)   	::  sys
!	 type(state),intent(in)    ::  st
!	 real(rl), intent(in) 	   ::  x
!	 integer, intent(in)       ::  n
!	 real(rl) 						::  imfn_x
!	 complex(cx)					::  tmp
!	 integer							::  s
!
!	 imfn_x = 0_rl
!	 tmp = 0_cx
!	 do s=1, sys%nmode
!		tmp = tmp + aimag( st%f(n,s) )* exp( Ic*(sys%w(s)/c_light)*2_rl*PI*x ) &
!														  * sqrt(1._rl/sys%nmode)
!	 end do
!	 imfn_x = real(tmp)
!
!  END FUNCTION 
!  FUNCTION imhn_x(sys,st,n,x)
!
!	 type(param),intent(in)   	::  sys
!	 type(state),intent(in)    ::  st
!	 real(rl), intent(in) 	   ::  x
!	 integer, intent(in)       ::  n
!	 real(rl) 						::  imhn_x
!	 complex(cx)					::  tmp
!	 integer							::  s
!
!	 imhn_x = 0_rl
!	 tmp = 0_cx
!	 do s=1, sys%nmode
!		tmp = tmp + aimag( st%h(n,s) )* exp( Ic*(sys%w(s)/c_light)*2_rl*PI*x ) &
!														  * sqrt(1._rl/sys%nmode)
!	 end do
!	 imhn_x = real(tmp)
!
!  END FUNCTION 
!
!  FUNCTION ref_x(sys,st,x)
!
!	 type(param),intent(in)   	::  sys
!	 type(state),intent(in)    ::  st
!	 real(rl), intent(in) 	   ::  x
!	 complex(cx) 					::  ref_x
!	 complex(cx)					::  tmp
!	 integer							::  n,m,s
!
!	 ref_x = 0_cx
!	 do s=1, sys%nmode
!		tmp = 0_cx
!		do n=1, st%np
!		  do m=1, st%np
!			 tmp = tmp + conjg(st%p(n))*st%p(m) * ( conjg(st%f(n,s)) + st%f(m,s) ) * st%ov_ff(n,m) / 2._rl
!		  end do
!		end do
!		tmp = tmp * exp( Ic*sys%w(s)*2_rl*PI*x ) * sqrt(sys%wmax/sys%nmode)
!		ref_x = ref_x + tmp
!	 end do
!
!  END FUNCTION 
!  FUNCTION imf_x(sys,st,x)
!
!	 type(param),intent(in)   	::  sys
!	 type(state),intent(in)    ::  st
!	 real(rl), intent(in) 	   ::  x
!	 complex(cx) 					::  imf_x
!	 complex(cx)					::  tmp
!	 integer							::  n,m,s
!
!	 imf_x = 0_cx
!	 do s=1, sys%nmode
!		tmp = 0_cx
!		do n=1, st%np
!		  do m=1, st%np
!			 tmp = tmp + conjg(st%p(n))*st%p(m) * ( st%f(m,s) - conjg(st%f(n,s)) ) * st%ov_ff(n,m) / (Ic*2._rl)
!		  end do
!		end do
!		tmp = tmp * exp( Ic*sys%w(s)*2_rl*PI*x ) * sqrt(sys%wmax/sys%nmode)
!		imf_x = imf_x + tmp
!	 end do
!
!  END FUNCTION
!  FUNCTION reh_x(sys,st,x)
!
!	 type(param),intent(in)   	::  sys
!	 type(state),intent(in)    ::  st
!	 real(rl), intent(in) 	   ::  x
!	 complex(cx) 					::  reh_x
!	 complex(cx)					::  tmp
!	 integer							::  n,m,s
!
!	 reh_x = 0_cx
!	 do s=1, sys%nmode
!	 	tmp = 0_cx
!		do n=1, st%np
!		  do m=1, st%np
!			 tmp = tmp + conjg(st%q(n))*st%q(m) * ( conjg(st%h(n,s)) + st%h(m,s) ) * st%ov_hh(n,m) / 2._rl
!		  end do
!		end do
!		tmp = tmp * exp( Ic*sys%w(s)*2_rl*PI*x ) * sqrt(sys%wmax/sys%nmode)
!		reh_x = reh_x + tmp
!	 end do
!
!  END FUNCTION 
!  FUNCTION imh_x(sys,st,x)
!
!	 type(param),intent(in)   	::  sys
!	 type(state),intent(in)    ::  st
!	 real(rl), intent(in) 	   ::  x
!	 complex(cx) 					::  imh_x
!	 complex(cx)					::  tmp
!	 integer							::  n,m,s
!
!	 imh_x = 0_cx
!	 do s=1, sys%nmode
!		tmp = 0_cx
!		do n=1, st%np
!		  do m=1, st%np
!			 tmp = tmp + conjg(st%q(n))*st%q(m) * ( st%h(m,s) - conjg(st%h(n,s)) ) * st%ov_hh(n,m) / (2._rl*Ic)
!		  end do
!		end do
!		tmp = tmp * exp( Ic*sys%w(s)*2_rl*PI*x ) * sqrt(sys%wmax/sys%nmode)
!		imh_x = imh_x + tmp
!	 end do
!
!  END FUNCTION
!
!  FUNCTION wigFH(x,p,f,h)
!  
!	 complex(cx), intent(in)  ::  f,h
!	 real(rl), intent(in)	  ::  x,p
!	 complex(cx)				  ::  wigFH
!
!	 wigFH = (2._rl/pi) * exp( - 0.5_rl*( 2._rl*p + Ic*(f-conjg(h) ) )**2  &
!										   - 0.5_rl*( 2._rl*x - (f+conjg(h) ) )**2	 &
!										   - 0.5_rl*( abs(f)**2 + abs(h)**2 - 2._rl*f*conjg(h) ) ) 
!
!  END FUNCTION
!  FUNCTION wigner_2pol(x,p,p1,p2,f1,f2)
!
!	 complex(cx), intent(in)  ::  p1,p2,f1,f2
!	 real(rl), intent(in)	  ::  x,p
!	 real(rl)	   			  ::  wigner_2pol
!	 real(rl)					  ::  normSq
!
!	 normSq = abs(p1)**2 + abs(p2)**2 &
!					 + 2._rl*real( p1*conjg(p2)*exp( - 0.5_rl*( abs(f1)**2 + abs(f2)**2 - 2._rl*conjg(f2)*f1 ) ) )
!
!	 wigner_2pol = real( (1._rl/normSq) * ( abs(p1)**2*wigFH(x,p,f1,f1) + abs(p2)**2*wigFH(x,p,f2,f2) &
!										  + p1*conjg(p2)*wigFH(x,p,f1,f2) + p2*conjg(p1)*wigFH(x,p,f2,f1) ) )
!
!
!  END FUNCTION
















