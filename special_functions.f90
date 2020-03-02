! special_functions.f90 by Brad Friesen
! Special functions not usually built in.(ex. fresnel integrals, gamma
! function, etc.)
! Compile with: gfortran -O3 -fdefault-real-8 -c special_functions.f90
! Updated: 2020-03-01

module special_functions
use parameters
implicit none

contains

subroutine fresnel (x, c, s)
! Computes the Fresnel integrals S(x) and C(x) for all real x.
! EPS is the relative error; MAXIT is the maximum number of iterations
! allowed; FPMIN is a number near the smallest representable
! floating-point number; XMIN is the dividing line between using the
! series and continued fraction; PI = :^) and PIBY2 = PI/2.
	integer MAXIT
	real x, c, s, EPS, FPMIN, PIBY2, XMIN
	parameter (EPS = 6.e-8, MAXIT = 100, FPMIN = 1.e-30, XMIN = 1.5, &
		PIBY2 = PI/2.0)
	integer k, n
	real a, absc, ax, fact, pix2, sign, sum, sumc, sums, term, test
	complex b, cc, d, h, del, cs
	logical odd
	
	! Statement function.
	absc(h) = abs(real(h))+abs(aimag(h))
	ax = abs(x)
	! Special case: avoid failure of convergence test due to underflow.
	if (ax .lt. sqrt(FPMIN)) then
		s = 0.
		c = ax
	! Evaluate both series simultaneously.
	else if (ax .le. XMIN) then
		sum = 0.
		sums = 0.
		sumc = ax
		sign = 1.
		fact = PIBY2*ax*ax
		odd = .true.
		term = ax
		n = 3
		do k = 1, MAXIT
			term = term*fact/k
			sum = sum+sign*term/n
			test = abs(sum)*EPS
			if (odd) then
				sign = -sign
				sums = sum
				sum = sumc
			else
				sumc = sum
				sum = sums
			end if
			if (term .lt. test) goto 1
			odd = .not. odd
			n = n+2
		end do
		write (stderr, *) 'series failed in fresnel'
1		s = sums
		c = sumc
	! Evaluate continued fraction by modified Lentz's method
	else
		pix2 = PI*ax*ax
		b = cmplx(1., -pix2)
		cc = 1./FPMIN
		d = 1./b
		h = d
		n = -1
		do k = 2, MAXIT
			n = n+2
			a = -n*(n+1)
			b = b+4
			! Denominators cannot be zero
			d = 1./(a*d+b)
			cc = b+a/cc
			del = cc*d
			h = h*del
			if (absc(del-1.) .lt. EPS) goto 2
		end do
		write (stderr, *) 'cf failed in fresnel'
2		h = h*cmplx(ax,-ax)
		cs = cmplx(.5, .5)*(1.-cmplx(cos(.5*pix2), sin(.5*pix2))*h)
		c = real(cs)
		s = aimag(cs)
	end if
	if (x .lt. 0.) then
		c = -c
		s = -s
	end if
	return
end subroutine fresnel

end module special_functions
