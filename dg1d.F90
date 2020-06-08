!========================================================================================================== !  
! 	Arbitrary order discontinuous Galerkin solver for 1D advection equation for conserved quantity q 
!  
! 		d/dt q(x) + d/dx f(q(x)) = 0 
!  
!========================================================================================================== 
program dg1d 

	use iso_fortran_env, only: real64, int64 
	use mod_gaussian_quadrature 
	
	implicit none 
	
	!-------------------------------------------------------------------------------------------------- !
	!	PARAMETERS 
	!-------------------------------------------------------------------------------------------------- 
	integer, parameter :: Np = 12 ! Number of Gauss-Legendre points (so order is N-1) 
	integer, parameter :: Nel = 12 ! Number of elements in our domain 
	integer, parameter :: Nt = 2000 ! Number of time steps 
	integer, parameter :: write_interval = 20 ! write data every x timesteps 
	
	real(real64), parameter :: DomainX1 = 0.0, DomainX2 = 12.0 ! Total domain to be discretized
	real(real64), parameter :: c = 0.6 ! Externally-prescribed constant velocity
	real(real64), parameter :: dt = 0.01 ! Time differential


	!--------------------------------------------------------------------------------------------------
	!	SOLVER VARIABLES
	!--------------------------------------------------------------------------------------------------

	real(real64) :: elementPtsLR(2,Nel) ! L and R points that make up each element, i.e. the mesh
	real(real64) :: linspace(Nel+1) ! Unique element endpoints, temporary var used in discretisation

	real(real64) :: gx(Np), gw(Np)  ! Gaussian quadrature points and weights

	integer(int64) :: i, j, l, m    ! Iterators, no special meaning
	integer :: t ! current time
	real(real64) :: stiff(Np,Np)    ! Stiffness matrix


	! Quantities of interest
	real(real64) :: q(Np,Nel)      ! nodal values of unknown (q), same for the upwind elmnt
	real(real64) :: qprime(Np,Nel) ! Time derivative of q-vector

	! Time-stepping intermediates for Runge-Kutta 4th order
	real(real64), dimension(Np,Nel) :: k1, k2, k3, k4


	!--------------------------------------------------------------------------------------------------
	!	SOLVER SETUP
	!--------------------------------------------------------------------------------------------------

	! Create mesh (discretize domain into elements)
	linspace = [(DomainX1 + (DomainX2-DomainX1) * (i-1) / Nel, i=1, Nel+1)] 
	elementPtsLR(1,:) = linspace(1:Nel)
	elementPtsLR(2,:) = linspace(2:Nel+1)

	! Calculate gaussian quadrature points and weights
	call gaussxw(Np, gx, gw) 


	! Initial conditions
	initcond: block
		real(real64) :: x_L, x_R            ! x coord of left and right element endpoint
		real(real64) :: xr(Np), w(Np)	    ! Gaussian quad. pts., shifted to the domain of a given elmnt
		real(real64) :: twopi = 2*acos(-1.0_real64)

		do i = 1, Nel
			x_L = elementPtsLR(1,i)
			x_R = elementPtsLR(2,i)

			xr = gx
			w = gw
			call shiftxw(xr,w,x_L,x_R)
			q(:,i) = sin( xr * twopi/(DomainX2-DomainX1) )
		end do
	end block initcond

	! Populate stiffness matrix
	! (never changes)
	do i = 1, Np
		do j = 1, Np
			! integral of L_m, L'_l
			stiff(j,i) = sum([( gw(m) * Legendre(gx,i,gx(m)) * LegendrePrime(gx+2*epsilon(gx),j,gx(m)) , m=1, Np )])
		end do
	end do	


	!--------------------------------------------------------------------------------------------------
	!	SOLVER 
	!--------------------------------------------------------------------------------------------------


	solver: block

		real(real64) :: xreal(Np), garbage(Np), pltx(Np*Nel), plty(Np*Nel)
		character(len=1024) :: t_string
		integer :: fileunit

		do t = 0, Nt-1

			if((mod(t,write_interval) .eq. 0) .or. (t .eq. Nt-1)) then
				do l = 1, Nel
					xreal = gx
					garbage = gw
					call shiftxw(xreal, garbage, elementPtsLR(1,l), elementPtsLR(2,l))

					pltx(Np*(l-1)+1:Np*l) = xreal
					plty(Np*(l-1)+1:Np*l) = q(:,l)
				end do

				print *, t
				write (t_string,"(i6.6)") t

				open (newunit=fileunit, file="dg1d."//trim(t_string)//".dat")
				write (fileunit, *) pltx
				write (fileunit, *) plty
				flush (fileunit)
				close (fileunit)

			end if

			! Runge-kutta 4th order

			call calc_qprime(Nel, Np, elementPtsLR, gx, gw, stiff, q, qprime)	
			k1 = qprime * dt

			call calc_qprime(Nel, Np, elementPtsLR, gx, gw, stiff, q + 0.5*k1, qprime)	
			k2 = qprime * dt

			call calc_qprime(Nel, Np, elementPtsLR, gx, gw, stiff, q + 0.5*k2, qprime)	
			k3 = qprime * dt

			call calc_qprime(Nel, Np, elementPtsLR, gx, gw, stiff, q + k3, qprime)	
			k4 = qprime * dt

			q = q + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0_real64

		end do


	end block solver


	

contains

	subroutine calc_qprime(Nel, Np, elementPtsLR, gx, gw, stiff, q, qprime)
		integer, intent(in) :: Nel, Np
		real(real64), intent(in) :: elementPtsLR(2,Nel), gx(Np), gw(Np), stiff(Np,Np)
		real(real64), intent(in) :: q(Np,Nel)
		real(real64), intent(out) :: qprime(Np,Nel)

		integer(int64) :: i, j ! Iterators, no special meaning

		integer :: iPrev   ! Index of previous element
		real(real64) :: x_L, x_R            ! x coord of left and right element endpoint
		real(real64) :: qR(Np), qRprev(Np)  ! q-vector on the right endpt of this elmnt and upwind elmnt
		real(real64) :: flx(Np) ! flux-vector



		! This is probably inefficient, but I'm going to loop over elements and solve them one at a time. 
		do i = 1, Nel
			iPrev = i-1
			if (iPrev .eq. 0) iPrev = Nel

			x_L = elementPtsLR(1,i)
			x_R = elementPtsLR(2,i)

			! Populate flx
			! Upwind difference, technically only valid for c > 0 I think
			qR     = sum([( Legendre(gx,j,1.0_real64)  * q(j,i)      , j=1, Np )])
			qRprev = sum([( Legendre(gx,j,1.0_real64) * q(j,iPrev)  , j=1, Np )])

			flx = c* [(  qR(j) * Legendre(gx,j,1.0_real64) - qRprev(j) * Legendre(gx,j,-1.0_real64) ,j=1,Np)]

			! Calculate qprime
			qprime(:,i) = 2.0/(gw*(x_R-x_L)) * (c * matmul(stiff,q(:,i)) - flx) 
		end do
	end subroutine calc_qprime


	! value of the k-th Legendre polynomial on the [-1,1] domain at point x
	! xr and w are passed in as the roots for Gauss-Legendre quadrature
	real(real64) pure function Legendre(xr, k, x)
		integer(int64),intent(in) :: k
		real(real64),intent(in) :: xr(:),  x

		real(real64) :: L
		integer(int64) :: m
		L = 1.0

		! Newman, p. 165
		do m = 1,size(xr)
			if (m /= k) then
				L = L * (x - xr(m))/(xr(k) - xr(m))
			end if
		end do
		Legendre = L
	end function Legendre

	! value of the derivative of the k-th Legendre polynomial on the [-1,1] domain at point x
	! xr and w are passed in as the roots for Gauss-Legendre quadrature
	! Calculation takes advantage of log differentiation, f'/f = log(f)'
	real(real64) pure function LegendrePrime(xr,k,x)
		integer(int64),intent(in) :: k
		real(real64),intent(in) :: xr(:), x

		real(real64) :: L
		integer(int64) :: m
		L = 0.0

		do m = 1,size(xr)
			if (m /= k) then
				L = L + 1.0/(x-xr(m))
			end if
		end do

		LegendrePrime = L * Legendre(xr,k,x)
	end function

#if 0
	subroutine simple_plot(x, y)
		real(real64), asynchronous, intent(in) :: x(:)
		real(real64), asynchronous, intent(in) :: y(:)

		integer :: pyerr
		type(module_py) :: plt
		type(tuple) :: args    
		type(ndarray) :: x_arr, y_arr

		pyerr = import_py(plt, "matplotlib.pyplot")
		errchk

		pyerr = ndarray_create_nocopy(x_arr, x)
		errchk

		pyerr = ndarray_create_nocopy(y_arr, y)
		errchk

		pyerr = tuple_create(args, 3)
		errchk

		pyerr = args%setitem(0, x_arr)
		errchk
		pyerr = args%setitem(1, y_arr)
		errchk
		pyerr = args%setitem(2, "o")
		errchk

		pyerr = call_py_noret(plt, "plot", args)
		errchk
		pyerr = call_py_noret(plt, "grid")
		errchk
		pyerr = call_py_noret(plt, "show")
		errchk

		call x_arr%destroy
		call y_arr%destroy
		call args%destroy
		call plt%destroy

	end subroutine

	subroutine wrp_forpy_initialize
	integer :: pyerr
	pyerr = forpy_initialize()
	errchk
	end subroutine wrp_forpy_initialize

	subroutine wrp_forpy_finalize
	call forpy_finalize
	end subroutine wrp_forpy_finalize
#endif

end program dg1d
