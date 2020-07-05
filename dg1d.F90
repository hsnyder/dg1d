!========================================================================================================== !  
! 	Arbitrary order discontinuous Galerkin solver for 1D advection equation for conserved quantity q 
!  
! 		d/dt q(x) + d/dx f(q(x)) = 0 
!  
!========================================================================================================== 

program dg1d 

	use iso_fortran_env, only: real64, int64 
	use m_quadrature
	
	implicit none 
	
	!-------------------------------------------------------------------------------------------------- !
	!	PARAMETERS 
	!-------------------------------------------------------------------------------------------------- 
	integer, parameter :: Np = 10 ! Number of Gauss-Lagrange points (so order is N-1) 
	integer, parameter :: Nel = 5 ! Number of elements in our domain 
	integer, parameter :: Nt = 2000 ! Number of time steps to run for
	integer, parameter :: write_interval = 20 ! write data every x timesteps 
	
	real(real64), parameter :: DomainX1 = 0.0, DomainX2 = 19.0 ! Total domain to be discretized
	real(real64), parameter :: c = 1.2 ! Externally-prescribed constant velocity
	real(real64), parameter :: dt = 0.01 ! Time differential


	!--------------------------------------------------------------------------------------------------
	!	SOLVER VARIABLES
	!--------------------------------------------------------------------------------------------------

	real(real64) :: elementPtsLR(2,Nel) ! L and R points that make up each element, i.e. the mesh
	real(real64) :: linspace(Nel+1) ! Unique element endpoints, temporary var used in discretisation

	real(real64) :: gx(Np), gw(Np)  ! Gaussian quadrature points and weights

	integer :: i, j, l, m    ! Iterators, no special meaning
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
	call GLL_Points(gx) 
	call GLL_Weights(gx, gw)


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
			stiff(j,i) = sum([( gw(m) * Lagrange(gx,i,gx(m)) * DLagrange(gx,j,gx(m)+epsilon(gx)) , m=1, Np )])
		end do
	end do	


	!--------------------------------------------------------------------------------------------------
	!	SOLVER 
	!--------------------------------------------------------------------------------------------------


	solver: block

		real(real64) :: xreal(Np), garbage(Np), pltx(Np*Nel), plty(Np*Nel)

		do t = 0, Nt-1

			! If we need to output some data for plotting...
			if((mod(t,write_interval) .eq. 0) .or. (t .eq. Nt-1)) then
				do l = 1, Nel
					xreal = gx
					garbage = gw
					call shiftxw(xreal, garbage, elementPtsLR(1,l), elementPtsLR(2,l))

					pltx(Np*(l-1)+1:Np*l) = xreal
					plty(Np*(l-1)+1:Np*l) = q(:,l)
				end do

				print *, t
				call savedata(t,pltx,plty)

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
! GL version
!	pure function qprime_el(Nel, Np, gx, stiff, massinv, lineJac, c, q, i )
!		integer, intent(in) :: Nel, Np, i ! i is the ith element
!		real(real64), intent(in) :: gx(Np), stiff(Np,Np)
!		real(real64), intent(in) :: q(Np,Nel), massinv(Np,Np), lineJac(Np)
!		real(real64), intent(in) :: c
!		real(real64) :: qprime_el(Np)
!
!		real(real64) :: x_L, x_R ! "left" and "right" element boundaries, in [-1,1] (i.e. either -1 or 1)
!
!		integer :: j ! Iterators, no special meaning
!
!		integer :: iPrev   ! Index of previous element
!		real(real64) :: qR(Np), qRprev(Np)  ! q-vector on the right endpt of this elmnt and upwind elmnt
!		real(real64) :: flx(Np) ! flux-vector
!
!		! To accomodate our upwind difference scheme for the flux, we need to know whether c is positive 
!		if (c >= 0) then
!			! Apply periodic BC, set iPrev to be the upwind element
!			iPrev = i-1
!			if (iPrev == 0) iPrev = Nel
!			x_L = -1.0_real64
!			x_R = 1.0_real64
!		else 
!			! Apply periodic BC, set iPrev to be the upwind element
!			! if c is negative we swap xL and xR because we need to evaluate at the opposite element boundary
!			iPrev = i+1
!			if (iPrev > Nel) iPrev = 1
!			x_L = 1.0_real64
!			x_R = -1.0_real64
!		end if 
!
!
!		! Populate flx
!		! Upwind difference 
!		! Flux = loss - gain
!		qR     = sum([( Lagrange(gx,j,x_R)  * q(j,i)      , j=1, Np )])
!		qRprev = sum([( Lagrange(gx,j,x_R)  * q(j,iPrev)  , j=1, Np )])
!
!		flx = abs(c)* [(  qR(j) * Lagrange(gx,j,x_R) - qRprev(j) * Lagrange(gx,j,x_L) ,j=1,Np)]
!
!		! Calculate qprime
!		qprime_el = lineJac * matmul(massinv, (c * matmul(stiff,q(:,i)) - flx) )
!
!	end function


	pure function qprime_el(Nel, Np, gx, stiff, massinv, lineJac, cBkwd, cFwd, q , prevBdyVal, nextBdyVal, i )
		integer, intent(in) :: Nel, Np, i ! i is the ith element
		real(real64), intent(in) :: gx(Np), stiff(Np,Np)
		real(real64), intent(in) :: q(Np), massinv(Np,Np), lineJac(Np), prevBdyVal, nextBdyVal
		real(real64), intent(in) :: cFwd, cBkwd
		real(real64) :: qprime_el(Np)

		integer :: j ! Iterators, no special meaning

		real(real64) :: flx(Np) ! flux-vector

				! Populate flx
		! Upwind difference 
		! Flux = loss - gain
		
		flx = 0
		if (cFwd>=0) then ! losing at top end, gaining at lower end
			flx(1) = cBkwd*prevBdyVal
			flx(Np) = cFwd*q(Np)
		else ! losing at bottom end, gaining at top end
			flx(1) = cBkwd*q(1)
			flx(Np) = cFwd*nextBdyVal
		end if

		! Calculate qprime
		qprime_el = lineJac * matmul(massinv, (c * matmul(stiff,q) - flx) )

	end function




	subroutine calc_qprime(Nel, Np, elementPtsLR, gx, gw, stiff, q, qprime)
		integer, intent(in) :: Nel, Np
		real(real64), intent(in) :: elementPtsLR(2,Nel), gx(Np), gw(Np), stiff(Np,Np)
		real(real64), intent(in) :: q(Np,Nel)
		real(real64), intent(out) :: qprime(Np,Nel) 

		real(real64) :: delta_x, Jac(Np) ! size of element

		integer :: i ! Iterators, no special meaning
		integer :: iPrev, iNext

		real(real64) :: massinv1d(Np,Np)
		massinv1d = 0
		do i = 1, Np
			massinv1d(i,i) = 1.0/gw(i)
		end do



		! This is probably inefficient, but I'm going to loop over elements and solve them one at a time. 
		do i = 1, Nel
			delta_x = elementPtsLR(2,i) - elementPtsLR(1,i)
			Jac = [(2.0/(delta_x), i = 1, Np)]
			! To accomodate our upwind difference scheme for the flux, we need to know whether c is positive 
			iNext = i+1
			iPrev = i-1

			! Periodic BC
			if(iNext > Nel) iNext = 1
			if(iPrev < 1) iPrev = Nel

			! Calculate qprime
			qprime(:,i) = qprime_el(Nel, Np, gx, stiff, massinv1d, Jac, -c, c, q(:,i), q(Np,iPrev), q(1,iNext), i)
		end do
	end subroutine calc_qprime

	subroutine savedata(t,x,y)
		integer, intent(in) :: t
		real(real64), intent(in) :: x(:), y(:)

		character(len=1024) :: t_string
		integer :: fileunit

		write (t_string,"(i6.6)") t

		open (newunit=fileunit, file="output."//trim(t_string)//".dat")
		write (fileunit, *) x
		write (fileunit, *) y
		flush (fileunit)
		close (fileunit)

	end subroutine savedata

	subroutine shiftxw(x,w,a,b)
		real(real64), intent(in out) :: x(:), w(:)
	       	real(real64), intent(in) :: a, b
		x = 0.5*(b-a)*x+0.5*(b+a)
		w = 0.5*(b-a)*w
	end subroutine


end program dg1d
