! ==============================================================================
!  FILE: consts.f90
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
!    * Requires BLAS/LAPACK (ZGESV, ZGETRF, ZTRSM, ZPOTRF, ...).
!
! ==============================================================================
!
module consts
!> -------------------------------------------------------------------------
!> MODULE: consts
!> -------------------------------------------------------------------------
!> Purpose / context:
!>   Module `consts`: central container for simulation types and routines.
!>   This module defines the core data structures (param/state/traj) and the
!>   variational equations of motion used during time evolution.
!> Arguments:
!>   (none)
!>

  use typedefs 
  implicit none

  real(r_type),   parameter :: c_light = 1_r_type
  real(r_type),   parameter :: pi=(4.0_r_type)*atan(1.0_r_type)
  real(r_type),   parameter :: one_r=1.0_r_type
  real(r_type),   parameter :: zero_r=0.0_r_type

  real(q_type),   parameter :: pi_q=(4.0_q_type)*atan(1.0_q_type)
  real(q_type),   parameter :: one_q=1.0_q_type
  real(q_type),   parameter :: zero_q=0.0_q_type

  complex(c_type),parameter :: im=(0.0d0,1.0d0)
  complex(c_type),parameter :: zero=(0.0d0,0.0d0)
  complex(c_type),parameter :: one=(1.0d0,0.0d0),two=(2.0d0,0.0d0)

contains

 !> -------------------------------------------------------------------------
 !> SUBROUTINE: alloc_check
 !> -------------------------------------------------------------------------
 !> Arguments:
 !>   - stat
 !>   - str
 !>
 subroutine alloc_check(stat,str)
    implicit none 
    integer           :: stat
    character(len=50) :: str

    if(stat.ne.0)then 
       print*, ': Allocation faliure of the matrix... quitting', &
            ' - ', str
       stop  
    end if
  end subroutine alloc_check

end module consts
