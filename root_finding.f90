! root_finding.f90 by Brad Friesen
! Root finding utilities
! Compile with: gfortran -O3 -fdefault-real-8 -c root_finding.f90
! Updated: 2020-03-01

module root_finding
implicit none

contains

!=======================================================================
! BRACKETING AND BISECTION
!=======================================================================

subroutine zbrak (fx, x1, x2, n, xb1, xb2, nb)
integer n, nb
real x1, x2, xb1(nb), xb2(nb), fx
external fx
! Given a function fx defined on the interval from x1-x2 subdivide the
! interval into n equally spaced segments, and search for zero crossings
! of the function. nb is input as the maximum number of roots sought,
! and is reset to the number of bracketing pairs xb1(1:nb), xb2(1:nb)
! that are found.
integer i, nbb
real  dx, fc, fp, x

nbb = 0
x = x1
! Determine the spacing appropriate to the mesh
dx = (x2-x1)/n
fp = fx(x)
! Loop over all intervals
do i = 1,n
	x = x+dx
	fc = fx(x)
	! If a sign change occurs then record values for the bounds
	if (fc*fp .le. 0.0) then
		nbb = nbb+1
		xb1(nbb) = x-dx
		xb2(nbb) = x
		write (*,*) nbb
		if (nbb .eq. nb) goto 10
		write (*,*) nbb
	end if
	fp = fc
end do
10 continue
write (*,*) 'test'
! nb = nbb
return
end subroutine



function rtbis (func, x1, x2, xacc)
integer jmax
real rtbis, x1, x2, xacc, func
external func
! Maximum allowed number of bisections
parameter (jmax=40)
! Using bisection, find the root of a function func known to lie between
! x1 and x2. The root, returned as rtbis, will be refined until accuracy
! is +- xacc.
integer j
real dx, f, fmid, xmid

fmid = func(x2)
f = func(x1)
if (f*fmid .ge. 0.) stop 'Root must be bracketed in rtbis'

! Orient the search so that f>0  lies at x+dx
if (f .lt. 0.) then
	rtbis = x1
	dx = x2-x1
else
	rtbis = x2
	dx = x1-x2
end if

do j = 1,jmax
	dx = dx*.5
	xmid = rtbis+dx
	fmid = func(xmid)
	if (fmid .le. 0.) rtbis = xmid
	if (abs(dx) .lt. xacc .or. fmid .eq. 0.) return
end do
stop 'too many bisections in rtbis'
end function

end module root_finding
