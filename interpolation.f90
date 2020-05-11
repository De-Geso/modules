! interpolation.f90 by Brad Friesen
! Stuff from Ch.3 on interpolation in numerical recipes
! compile with: gfortran -o3 -fdefault-real-8 -c functions.f90
! updated: 2020-04-26

module interpolation
implicit none


contains


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

end module
