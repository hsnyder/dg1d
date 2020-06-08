! Not real source file. Rather, this is included several times with different preprocessor defines
! When including this, define the following constants: 
! - TYPE: real32 or real64
! - EPSVAL: value of epsilon for Newton's method (use _kind)
! - ONE: value of 1 in the exact kind you're using for TYPE
		integer, intent(in) :: N 
		real(TYPE), intent(out) :: x(N), w(N) 
		real(TYPE), parameter :: pi = acos( -ONE )
		real(TYPE), parameter :: eps = EPSVAL

		real(TYPE) :: a(N), linsp_start, linsp_end
		integer :: k , i
		real(TYPE) :: delta 
		real(TYPE) :: dx(N), dp(N), p0(N), p1(N), ptmp(N) 

		delta = 1.0
		dp = 0
		p1 = 0
		p0 = 0
		ptmp = 0

		! linspace(3, 4*N-1, N) 
		linsp_start = 3.0
		linsp_end = 4.0*N-1
		a = [(linsp_start + (linsp_end-linsp_start) * (i-1) / (N-1), i=1, N)] 

		a = a / (4*N+2) 

		! Initial guess of Legendre polynomial roots 
		x = cos(pi*a + 1/(8*N*N*tan(a))) 
		
		! Find roots using Newton's method 
		do while (delta > eps) 
			p0 = [(1, i = 1, N)] 
			p1 = x 
			do k = 1, N-1 
				ptmp = p1 
				p1 = ((2*k+1)*x*p1-k*p0)/(k+1)
				p0 = ptmp 
			end do 
			dp = (N+1)*(p0-x*p1)/(1-x*x)
			dx = p1/dp 
			x = x - dx 
			delta = maxval(abs(dx)) 
		end do 
		! Calculate the weights 
		w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)
