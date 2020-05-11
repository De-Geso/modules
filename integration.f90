! integration.f90 by Brad Friesen
! Stuff from Ch.4 on integration in numerical recipes
! More or less ordered in terms of power.
! compile with: gfortran -o3 -fdefault-real-8 -c integration.f90
! updated: 2020-05-01

module integration
use interpolation
implicit none


contains

!=======================================================================
! ELEMENTARY ALGORITHMS
!=======================================================================

subroutine trapzd (func, a, b, s, n)
integer n
real a, b, s, func
external func
! This routine computes the nth stage of refinement of an extended
! trapezoidal rule. fun is input as the name of the function to be
! integrated between the limits a and b, also input. When called with
! n=1, the routine returns as s the crudest estimate of the integral.
! subsequent calls with n=2,3,... (in that sequential order) will
! improve the accuracy of s by adding 2^(n-2) additional interior
! points. s should not be modified between sequential calls.
integer it, j
real del, sum, tnm, x

if (n .eq. 1) then
	s = 0.5*(b-a)*(func(a)+func(b))
else
	it = 2**(n-2)
	tnm = it
	del = (b-a)/tnm
	! Spacing of points to be added.
	x = a+0.5*del
	sum = 0.0
	do j=1, it
		sum = sum+func(x)
		x = x+del
	end do
	s = 0.5*(s+(b-a)*sum/tnm)
	! replaces s by its refined value
end if

return
end subroutine


subroutine qtrap (func, a, b, s)
integer jmax
real a, b, func, s, eps
external func
parameter (eps=1.0e-6, jmax=20)
! uses trapzd

! Returns as s the integral of the function func from a to b. The
! parameters eps can be set to the desired fractional accuracy and jmax
! so that 2^(jmax-1) is the maximum allowed number of steps. Integration
! is performed by the trapezoidal rule.
integer j
real olds

olds = 0.0		! initial value of olds is arbitrary.
do j = 1, jmax
	call trapzd(func, a, b, s, j)
	if (j .gt. 5) then		! avoid spurious early convergence.
		if (abs(s-olds) .lt. eps*abs(olds) .or. &
			(s .eq. 0.0 .and. olds .eq. 0.0)) return
	end if
	olds = s
end do

stop 'too many steps in qtrap'
end subroutine


subroutine qsimp (func, a, b, s)
integer jmax
real a, b, func, s, eps
external func
parameter (eps=1.0e-6, jmax=20)
!uses trapzd

! Returns as s the integral of the function func from a to b. The
! parameters eps can be set to the desired fractional accuracy and jmax
! so that 2^(jmax-1) is the maximum allowed number of steps. Integration
! is performed by Simpson's rule.
integer j
real os, ost, st

ost = 0.0
os = 0.0
do j = 1, jmax
	call trapzd(func, a, b, st, j)
	s = (4.0*st-ost)/3.0
	if (j .gt. 5) then		! avoid spurious early convergence.
		if (abs(s-os) .lt. eps*abs(os) .or. &
			(s .eq. 0.0 .and. os .eq. 0.0)) return
	end if
	os = s
	ost = st
end do

stop 'too many steps in qsimp'
end subroutine

!=======================================================================
! ROMBERG INTEGRATION
!=======================================================================

subroutine qromb (func, a, b, ss)
integer jmax, jmaxp, k, km
real a, b, func, ss, eps
external func
parameter (eps=1.0e-6, jmax=20, jmaxp=jmax+1, k=5, km=k-1)
! Uses polint, trapzd
! Returns as ss the integral of the function func from a to b.
! Integration is performed by Romberg's method of order 2k, where e.g.,
! k=2 is Simpson's rule.
! Parameters: eps is the fractional accuracy desired, as determined by
! the extrapolation error estimate; jmax limits the total number of
! steps; k is the number of points used in the extrapolation.
integer j
! These store the successive trapezoidal approximations and their
! respective step sizes.
real dss, h(jmaxp), s(jmaxp)

h(1) = 1.0
do j = 1, jmax
	call trapzd (func, a, b, s(j), j)
	if (j .ge. k) then
		call polint(h(j-km), s(j-km), k, 0.0, ss, dss)
		if (abs(dss) .le. eps*abs(ss)) return
	end if
	s(j+1) = s(j)
	h(j+1) = 0.25*h(j)
end do

stop 'too many steps in qromb'
end subroutine

!=======================================================================
! IMPROPER INTEGRALS
!=======================================================================
! Integrand goes to a finite limiting value at finite upper and lower
!	limits, but cannot be evaluated right on one of those limits.
! The upper limit is infinity, or the lower limit is -infinity.
! It has an integrable singularity at either limit (ex. x^-1/2 at x=0).
! It has an integrable singularity at a known place between the upper
!	and lower limits.

subroutine midpnt (func, a, b, s, n)
integer n
real a, b, s, func
external func
! This routine computes the nth stage of refinement of an extended
! midpoint rule. func is input as the name of the function to be
! integrated between limits a and b, also input. When called with n=1,
! the routine returns as s the crudest estimate of the integral.
! subsequent calls with n=2,3,... (in that sequential order) will
! improve the accuracy of s by adding (2/3)*3^(n-1) additional interior
! points. s should not be modified between sequential calls.
integer it, j
real ddel, del, sum, tnm, x
if (n .eq. 1) then
	s = (b-a)*func(0.5*(a+b))
else
	it = 3**(n-2)
	tnm = it
	del = (b-a)/(3.0*tnm)
	ddel = del+del
	! The added points alternate in spacing between del and ddel
	x = a+0.5*del
	sum = 0.0
	do j = 1, it
		sum = sum+func(x)
		x = x+ddel
		sum = sum+func(x)
		x = x+del
	end do
	s = (s+(b-a)*sum/tnm)/3.0
	! The new sum is combined with the old integral to give a refined
	! integral
end if

return
end subroutine

subroutine qromo (func, a, b, ss, choose)
integer jmax, jmaxp, k, km
real a, b, func, ss, eps
external func, choose
parameter (eps=1.e-6, jmax=14, jmaxp=jmax+1, k=5, km=k-1)
! USES polint
! Romberg integration on an open interval. Returns as ss the integral of
! the function func from a to b, using any specified integrating
! subroutine choose and Romberg's method. Normally choose will be an
! open formula, not evaluating the function at the endpoints. It is
! assumed that choose triples the number of steps on each call, and that
! its error series contains only even powers of the number of steps. The
! routines midpnt, midinf, midsql, midsqu, are possible choices for
! choose. The parameters have the same meaning as in qromb.
integer j
real dss, h(jmaxp), s(jmaxp)

h(1) = 1.0
do j = 1, jmax
	call choose(func, a, b, s(j), j)
	if (j .ge. k) then
		call polint(h(j-km), s(j-km), k, 0.0, ss, dss)
		if (abs(dss) .le. eps*abs(ss)) return
	end if
	s(j+1) = s(j)
	h(j+1) = h(j)/9.0
	! This is where the assumption of step tripling and an even error
	! series is used.
end do

stop 'too many steps in qromo'
end subroutine

subroutine midinf (funk, aa, bb, s, n)
integer n
real aa, bb, s, funk
external funk
! This integral is an exact replacement for midpnt, i.e. returns as s
! the nth stage of refinement of the integral of funk from aa to bb,
! except that the function is evaluated at evenly spaced points in 1/x
! rather than in x. This allows the upper limit bb to be as large and
! positive as the computer allows, or the lower limit aa to be as large
! and negative, but not both. aa dn bb must have the same sign.
integer it, j
real a, b, ddel, del, sum, tnm, func, x

func(x) = funk(1.0/x)/x**2
! This statement function effects the change of variable.
b = 1.0/aa
a = 1.0/bb
! These two statements change the limits of integration accordingly
! From this point on, the routine is exactly identical to midpnt.
if (n .eq. 1) then
	s = (b-a)*func(0.5*(a+b))
else
	it = 3**(n-2)
	tnm = it
	del = (b-a)/(3.0*tnm)
	ddel = del+del
	! The added points alternate in spacing between del and ddel
	x = a+0.5*del
	sum = 0.0
	do j = 1, it
		sum = sum+func(x)
		x = x+ddel
		sum = sum+func(x)
		x = x+del
	end do
	s = (s+(b-a)*sum/tnm)/3.0
	! The new sum is combined with the old integral to give a refined
	! integral
end if

return
end subroutine

end module
