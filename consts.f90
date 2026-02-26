! ==============================================================================
!  FILE / MODULE: consts.f90 (consts)
!
!  PURPOSE & CONTEXT
!    Central container for physical constants, mathematical constants, and 
!    core utility routines for the waveguide QED simulation.
!
!  CORE RESPONSIBILITIES
!    1. Constants : Defines precision-independent mathematical and physical 
!                   constants (e.g., pi, the imaginary unit, zero/one aliases) 
!                   to ensure numerical consistency across the entire codebase.
!    2. Utilities : Provides defensive programming safety checks (e.g., memory 
!                   allocation validation) to handle runtime errors.
!
! ==============================================================================
!
module consts

  use typedefs 
  implicit none

! -- Real (r_type) Constants --
  real(r_type),   parameter :: c_light = 1_r_type
  real(r_type),   parameter :: pi=(4.0_r_type)*atan(1.0_r_type)
  real(r_type),   parameter :: one_r=1.0_r_type
  real(r_type),   parameter :: zero_r=0.0_r_type

! -- Quadruple/Extended Precision (q_type) Constants --
  real(q_type),   parameter :: pi_q=(4.0_q_type)*atan(1.0_q_type)
  real(q_type),   parameter :: one_q=1.0_q_type
  real(q_type),   parameter :: zero_q=0.0_q_type

! -- Complex (c_type) Constants --
  complex(c_type),parameter :: im=(0.0d0,1.0d0)
  complex(c_type),parameter :: zero=(0.0d0,0.0d0)
  complex(c_type),parameter :: one=(1.0d0,0.0d0),two=(2.0d0,0.0d0)

contains

!> -------------------------------------------------------------------------
 !> SUBROUTINE: alloc_check
 !> -------------------------------------------------------------------------
 !> Purpose / context:
 !>   Defensive programming utility for safe dynamic memory allocation. 
 !>   Evaluates the status flag returned by Fortran's `allocate` statements. 
 !>   If the operating system fails to grant memory (status != 0), it catches 
 !>   the error, prints a custom traceback string identifying the point of 
 !>   failure, and safely halts execution to prevent uncontrolled segmentation faults.
 !> Arguments:
 !>   - stat : The integer status flag returned by an `allocate(..., stat=...)` call.
 !>   - str  : A string tag identifying the array or module where allocation failed.
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
