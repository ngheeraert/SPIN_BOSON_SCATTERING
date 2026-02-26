! ==============================================================================
!  FILE / MODULE: typedefs.f90 (typedefs)
!
!  PURPOSE & CONTEXT
!    Central definition of numerical precision (kinds) for the waveguide QED 
!    simulation. This infrastructure ensures that the code can be compiled in 
!    either single or double precision without modifying the physics routines. 
!
!  CORE RESPONSIBILITIES
!    1. Kind Selection : Defines standard integer and real kind parameters 
!                        (single, double, and quadruple precision) using 
!                        portable Fortran intrinsics.
!    2. Type Aliasing  : Maps the generic aliases `r_type` (real) and `c_type` 
!                        (complex) to specific precisions via preprocessor 
!                        directives.
!
!  COMPILATION NOTE
!    - If the preprocessor symbol `DP` is defined during compilation (e.g., 
!      -DDP), the simulation runs in 64-bit double precision.
!    - Otherwise, it defaults to 32-bit single precision.
!
! ==============================================================================
!
MODULE typedefs

  implicit none

! -- Integer Kinds --
  integer, parameter :: i1b=selected_int_kind(2)
  integer, parameter :: i2b=selected_int_kind(4)
  integer, parameter :: i4b=selected_int_kind(8)

! -- Floating Point Kinds --
  integer, parameter :: sp=kind(1.0)			! Single precision
  integer, parameter :: spc=kind((1.0,1.0))		! Single complex
  integer, parameter :: dp=selected_real_kind(15)	! Double precision
  integer, parameter :: qp=selected_real_kind(2*precision(1.0_dp))	! Quadruple precision
  integer, parameter :: dpc=kind((1.0_dp,1.0_dp))	! Double complex
  integer, parameter :: qpc=kind((1.0_qp,1.0_qp))	! Quadruple complex
#ifdef DP
! -- Double Precision Configuration --
  integer, parameter :: r_type = dp 
  integer, parameter :: q_type = qp
  integer, parameter :: c_type = dpc 
  integer, parameter :: qc_type = qpc 
#else
! -- Single Precision Configuration (Lightweight) --
  integer, parameter :: r_type = sp 
  integer, parameter :: c_type = spc 
#endif

END module typedefs
