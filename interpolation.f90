! interpolation.f90 by Brad Friesen
! Stuff from Ch.3 on interpolation in numerical recipes
! compile with: gfortran -o3 -fdefault-real-8 -c functions.f90
! updated: 2020-04-26

module interpolation
implicit none


contains


!=======================================================================
! 3.1 POLYNOMIAL INTERPOLATION AND EXTRAPOLATION
!=======================================================================
subroutine polint (xa, ya, n, x, y, dy)
integer n, nmax
real dy, x, y, xa(n), ya(n)
! Largest anticipated value of n
parameter (nmax=10)
! Given arrays xa and ya, each of length n, and given a value x, this
! routine returns a value y, and an error estimate dy. If p(x) is the
! polynomial of degree n-1 such that p(xa_i) = ya_i, i =1,...,n, then
! the returned value y=p(x).
integer i, m, ns
real den, dif, dift, ho, hp, w, c(nmax), d(nmax)

ns = 1
dif = abs(x-xa(1))

! Find the index, ns, of the closest table entry.
do i = 1, n
	dift = abs(x-xa(i))
	if (dift .lt. dif) then
		ns = i
		dif = dift
	end if
	c(i) = ya(i)
	d(i) = ya(i)
end do

y = ya(ns)
ns = ns-1

do m = 1, n-1
	do i = 1, n-m
		ho = xa(i)-x
		hp = xa(i+m)-x
		w = c(i+1)-d(i)
		den = ho-hp
		! This error occurs if two xa's are identical.
		if (den .eq. 0.0) write(*,*) 'failure in polint'

		den = w/den
		! Update c and d
		d(i) = hp*den
		c(i) = ho*den
	end do
	if (2*ns .lt. n-m) then
		dy = c(ns+1)
	else
		dy = d(ns)
		ns = ns-1
	end if
	y = y+dy
end do

return
end subroutine


!=======================================================================
! 3.2 RATIONAL FUNCTION INTERPOLATION AND EXTRAPOLATION
!=======================================================================


subroutine ratint (xa, ya, n, x, y, dy)
integer n, nmax
real dy, x, y, xa(n), ya(n), tiny
! Largest expected value of n, and a small number
parameter (nmax = 10, tiny = 1.e-25)
! Given arrays xa and ya, each of length n, and given a value of x, this
! routine returns a value of y and an accuracy estimate dy. The value
! returned is that of the diagonal rational function, evaluated at x,
! which passes through n points (xa_i, ya_i), i=1,...,n.
integer i, m, ns
real dd, h, hh, t, w, c(nmax), d(nmax)

ns = 1
hh = abs(x-xa(1))

do  i = 1, n
	h = abs(x-xa(i))
	if (h .eq. 0.0) then
		y = ya(i)
		dy = 0.0
		return
	else if (h .lt. hh) then
		ns = i
		hh = h
	end if
	c(i) = ya(i)
	d(i) = ya(i) + tiny
	! Tiny i s needed to prevent a zero-over-zero condition
end do

y = ya(ns)
ns = ns-1

do m = 1, n-1
	do i = 1, n-m
		w = c(i+1)-d(i)
		h = xa(i+m)-x
		! h will never be zero, since this was tested in the initial loop
		t = (xa(i)-x)*d(i)/h
		dd = t-c(i+1)

		! This error condition indicates that the interpolating function
		! has a pole at the requested value of x.
		if (dd .eq. 0.0) write(*,*) 'failure in ratint'

		dd = w/dd
		d(i) = c(i+1)*dd
		c(i) = t*dd
	end do
	if (2*ns .lt. n-m) then
		dy = c(ns+1)
	else
		dy = d(ns)
		ns = ns-1
	end if
	y = y+dy
end do
return
end subroutine


!=======================================================================
! 3.3 CUBIC SPLINE INTERPOLATION
!=======================================================================


subroutine spline (x, y, n, yp1, ypn, y2)
integer n, nmax
real yp1, ypn, x(n), y(n), y2(n)
parameter (nmax=100001)
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.,
! y_i = f(x_i), with x_1<x_2<x_N, and given values yp1 and ypn for the
! first derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2(1:n) of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x_i. If yp1 and/or ypn are equal to 1.e30 or larger,
! the routine is signaled to set the correspoinding boundary condition
! for a natural spline, with zero second derivative on that boundary.
! Parameter: nmax is the largest anticipated value of n.
integer i, k
real p, qn, sig, un, u(nmax)

! The lower boundary condition is set either to be
! "natural"
if (yp1 .gt. 0.99e30) then
	y2(1) = 0.
	u(1) = 0.
! or else to have a specified first derivative
else
	y2(1) = -0.5
	u(1) = (3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
end if
! This is the decomposition loop of the tridiagonal algorithm. y2 and u
! are used for temporary storage of the decomposed factors.
do i = 2,n
	sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
	p = sig*y2(i-1)+2.
	y2(i) = (sig-1.)/p
	u(i) = (6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
		/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
end do
! The upper boundary condition is set either to be
! "natural"
if (ypn .gt. 0.99e30) then
	qn = 0.
	un = 0.
! or else to have a specified first derivative
else
	qn = 0.5
	un = (3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
end if
y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.)
! This is the backsubstitution loop of the tridiagonal algorithm
do k = n-1, 1, -1
	y2(k) = y2(k)*y2(k+1)+u(k)
end do
return
end subroutine



subroutine splint (xa, ya, y2a, n, x, y)
integer n
real x, y, xa(n), y2a(n), ya(n)
! Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a
! function (with the xa_i's in order), and given the array y2a(1:n),
! which is the output from spline above, and given a value of x, this
! routine returns a cubic-spline interpolated value y.
integer k, khi, klo
real a, b, h

! We will find the right place in the table by means of bisection. This
! is optimal if sequential calls to this routine are at random values of
! x.
klo = 1
khi = n
1 if (khi-klo .gt. 1) then
	k = (khi+klo)/2
	if (xa(k) .gt. x) then
		khi = k
	else
		klo = k
	end if
	go to 1
end if
! klo and khi now bracket the input value of x.
h = xa(khi)-xa(klo)
! The xa's must be distinct.
if (h .eq. 0.) stop 'bad xa input in splint'
! Cubic spline polynomial is now evaluated
a = (xa(khi)-x)/h
b = (x-xa(klo))/h
y = a*ya(klo)+b*ya(khi) + &
	((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
return
end subroutine

end module
