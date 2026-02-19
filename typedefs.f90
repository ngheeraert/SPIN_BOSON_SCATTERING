! ==============================================================================
!  FILE: typedefs.f90
!
!  PURPOSE
!    Original research code used to produce the numerical data underlying the figures of:
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
! ==============================================================================
!
  module typedefs
  !> -------------------------------------------------------------------------
  !> MODULE: typedefs
  !> -------------------------------------------------------------------------
  !> Purpose / context:
  !>   Module `typedefs`: central container for simulation types and routines.
  !>   This module defines the core data structures (param/state/traj) and the
  !>   variational equations of motion used during time evolution.
  !> Arguments:
  !>   (none)
  !>
  implicit none
  integer, parameter :: i1b=selected_int_kind(2)
  integer, parameter :: i2b=selected_int_kind(4)
  integer, parameter :: i4b=selected_int_kind(8)
  integer, parameter :: sp=kind(1.0)
  integer, parameter :: spc=kind((1.0,1.0))
  integer, parameter :: dp=selected_real_kind(15)
  integer, parameter :: qp=selected_real_kind(2*precision(1.0_dp))
  integer, parameter :: dpc=kind((1.0_dp,1.0_dp))
  integer, parameter :: qpc=kind((1.0_qp,1.0_qp))
#ifdef DP
  integer, parameter :: r_type = dp 
  integer, parameter :: q_type = qp
  integer, parameter :: c_type = dpc 
  integer, parameter :: qc_type = qpc 
#else
  integer, parameter :: r_type = sp 
  integer, parameter :: c_type = spc 
#endif
end module typedefs
