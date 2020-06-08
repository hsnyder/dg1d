module mod_gaussian_quadrature

	use iso_fortran_env, only: real32, real64
	implicit none

	private
	public :: gaussxw, shiftxw

	interface gaussxw
		module procedure :: gaussxw_fp32
		module procedure :: gaussxw_fp64
	end interface gaussxw

	interface shiftxw
		module procedure :: shiftxw_fp32
		module procedure :: shiftxw_fp64
	end interface shiftxw


contains


	subroutine gaussxw_fp32 (N, x, w)
#define TYPE real32
#define EPSVAL 1e-7_real32
#define ONE 1.0_real32
#include "impl_gaussxw.F90"
#undef ONE
#undef EPSVAL
#undef TYPE
	end subroutine 

	subroutine gaussxw_fp64 (N, x, w)
#define TYPE real64
#define EPSVAL 1e-15_real64
#define ONE 1.0_real64
#include "impl_gaussxw.F90"
#undef ONE
#undef EPSVAL
#undef TYPE
	end subroutine 

	subroutine shiftxw_fp32(x,w,a,b) 
		real(real32), intent(in out) :: x(:), w(:)
	       	real(real32), intent(in) :: a, b
		x = 0.5*(b-a)*x+0.5*(b+a)
		w = 0.5*(b-a)*w
	end subroutine

	subroutine shiftxw_fp64(x,w,a,b)
		real(real64), intent(in out) :: x(:), w(:)
	       	real(real64), intent(in) :: a, b
		x = 0.5*(b-a)*x+0.5*(b+a)
		w = 0.5*(b-a)*w
	end subroutine



end module
