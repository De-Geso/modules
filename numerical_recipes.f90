! numerical_recipes.f90 by Brad Friesen
! Compile with: gfortran -O3 -fdefault-real-8 -c functions.f90
! Updated: 2019-11-29

module numerical_recipes
implicit none

contains

function ran2 (idum)
	integer (idum, im1, im2, imm1, ia1, ia2, iq1, iq2, ir1, ir2, ntab, 
		ndiv)
	real ran2, am, eps, rnmx
	parameter (im1=2147483563, im2=2147483399, am=1./im1, imm1=im1-1, & 
		ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211, &
		ir2=3791, ntab=32, ndiv=1+imm1/ntab, eps=1.2e-7, rnmx=1.-eps)
	! Long period (>2x10^18 random number generator of L'Ecuyer with
	! Bays-Durham shuffle and added safeguards. Returns a uniform
	! random deviate between 0.0 and 1.0 (exclusive of the endpoint
	! values). Call with idum a negative integer to initialize;
	! thereafter, do not alter idum between successive deviates in a
	! sequence. rnmx should approximate the largest floating value that
	! is less than 1.
	integer idum2, j, k, iv(ntab), iy
	save iv, iy, idum2
	data idum2/123456789/, iv/ntab*0/, iy/0/
	
	! initialize
	if (idum <= 0) then
	! be sure to prevent idum=0
		idum = max(-idum,1)
		idum2 = idum
		! load the shuffle table (after 8 warmups)
		do j = ntab+8, 1, -1
			k = idum/iq1
			idum = ia1*(idum-k*iq1)-k*ir1
			if (idum < 0) idum = idum+im1
			if (j <= ntab) iv(j) = idum
		end do
		iy = iv(1)
	end if
	! start here when not initializing
	k = idum/iq1
	! compute idum = mod(ia1*idum,im1) without overflows by Schrage's
	! method
	idum = ia1*(idum-k*iq1)-k*ir1
	if (idum < 0) idum = idum+im1
	k = idum2/iq2
	idum2 = ia2*(idum2-k*iq2)-k*ir2
	if (idum2 < 0) idum2 = idum2+im2
	j = 1+iy/ndiv
	iy = iv(j)-idum2
	iv(j)=idum
	if (iy < 1) iy = iy+imm1
	ran2 = min(am*iy,rnmx)
	return
end

end module numerical_recipes
