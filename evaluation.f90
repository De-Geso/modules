! evaluation.f90 by Brad Friesen
! Stuff from Ch.5 on evaluation of functions in numerical recipes
! compile with: gfortran -o3 -fdefault-real-8 -c evaluation.f90
! updated: 2020-05-01

module evaluation
implicit none


contains

!=======================================================================
! NUMERICAL DERIVATIVES
!=======================================================================

function dfridr (func, x, h, err)
integer ntab
real dfridr, err, h, x, func, con, con2, big, safe
parameter (con=1.4, con2=con*con, big=1.0e30, ntab=10, safe=2.0)
external func
! Returns the derivative of a function func at a point x by Ridders'
! method  of polynomial extrapolation. The value h is input as an
! estimated initial stepsize; it need not be small, but rather should be
! an increment in x over which func changes substantially. An estimate
! of the error in the derivative is returned as err. Parameters stepsize
! is decreased by con at each iteration. Max size of tableau is set by
! ntab. Return when error is safe worse than the best so far.
integer i, j
real errt, fac, hh, a(ntab,ntab)

if (h .eq. 0.0) stop 'h must be nonzero in dfridr'
hh=h
a(1,1) = (func(x+hh)-func(x-hh))/(2.0*hh)
err = big

! Successive columns in the Neville tableau will go to smaller stepsizes
! and higher orders of extrapolation.
do i = 2, ntab
	hh = hh/con
	! Try new, smaller step size.
	a(1,i) = (func(x+hh)-func(x-hh))/(2.0*hh)
	fac = con2
	! Compute extrapolations of various orders, requiring no new
	! function evaluations.
	do j = 2, i
		a(j,i) = (a(j-1,i)*fac-a(j-1,i-1))/(fac-1.0)
		fac = con2*fac
		! The error strategy is to compare each new extrapolation to one
		! order lower, both at the present step size and the previous
		! one.
		errt = max(abs(a(j,i)-a(j-1,i)), abs(a(j,i)-a(j-1,i-1)))
		! If error is decreased, save the improved answer.
		if (errt .le. err) then
			err = errt
			dfridr = a(j,i)
		end if
	end do
	! If higher order is worse by a significant factor SAFE, then quit
	! early.
	if (abs(a(i,i)-a(i-1,i-1)) .ge. safe*err) return
end do

return
end function

!=======================================================================
! CHEBYSHEV APPROXIMATION
!=======================================================================

subroutine chebft (a, b, c, n, func)
integer n, nmax
real a, b, c(n), func
external func
parameter (nmax=50)
! Chebyshev fit: Given a function func, lower and upper limits of the
! interval [a,b], and a maximum degree n, this routine computes the n
! coefficients. This routine is to be used with moderately large n
! (e.g., 30 or 50), the array of c's subsequently to be truncated at the
! smaller value mm such that c_{m+1} and subsequent elements are
! negligible.
! Parameters: Maximum expected value of n, and pi.
integer j, k
real bma, bpa, fac, y, f(nmax)
double precision sum
bma = 0.5*(b-a)
bpa = 0.5*(b+a)

! We evaluate the function at the n points required by (5.8.7)
do k = 1, n
	y = cos(pi*(k-0.5d0)/n)
	f(k) = func(y*bma+bpa)
end do

fac = 2.0/n

do j = 1, n
	sum = 0.d0
	do k = 1, n
		sum = sum+f(k)*cos((pi*(j-1))*((k-0.5d0)/n))
	end do
	c(j) = fac*sum
end do

return
end subroutine

end program
