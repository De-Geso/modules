! linear_algebra.f90 by Brad Friesen
! Stuff from Ch.2 on solutions of linear algebraic equations
! compile with: gfortran -o3 -fdefault-real-8 -c linear_algebra.f90
! updated: 2020-05-23

module linear_algebra
implicit none


contains


!=======================================================================
! TRIDIAGONAL AND BAND DIAGONAL SYSTEMS OF EQUATIONS
!=======================================================================


subroutine tridag (a, b, c, r, u, n)
integer n, nmax
real a(n), b(n), c(n), r(n), u(n)
parameter (nmax=1000)
! solves for a vector u(1:n) of length n the tridiagonal linear set
! given by eq. (2.4.1). a(1:n), b(1:n), c(1:n), and r(1:n) are input
! vectors and are not modified. Parameter: nmax is the maximum expected
! value of n.
integer j
! one vector of workspace, gam is needed.
real bet, gam(nmax)
! if this happens then you should rewrite your equations as a set of
! order n-1, with u2 trivially eliminated.
if (b(1) .eq. 0.) stop 'tridag: rewrite equtions'
bet = b(1)
u(1) = r(1)/bet
do j = 2,n
	gam(j) = c(j-1)/bet
	bet = b(j) - a(j)*gam(j)
	if (bet .eq. 0.) stop 'tridag failed'
	u(j) = (r(j) - a(j)*u(j-1))/bet
end do
do j = n-1, 1, -1
	u(j) = u(j) - gam(j+1)*u(j+1)
end do
return
end subroutine


subroutine banmul (a, n, m1, m2, np, mp, x, b)
integer m1, m2, mp, n, np
real a(np,mp), b(n), x(n)
! Matrix multiply Ax=b where A is band diagonal with m1 rows below the
! diagonal and m2 rows above. The input vector x and output vector b are
! stored as x(1:n) and b(1:n) respectively. The array A(1:n,1:m1+m2+1)
! stores A as follows: The diagonal elements are in A(1:n,m1+1).
! Subdiagonal elements are in A(j:n,1:m1) with j>1 appropriate to the
! number of elements on each subdiagonal. Superdiagonal elements are in
! A(1:j,m1+2:m1+m2+1) with j<n appropriate to the number of elements on
! each superdiagonal.
integer i, j, k

do i = 1,n
	b(i) = 0.
	k = i-m1-1
	do j = max(1,1-k), min(m1+m2+2,n-k)
		b(i) = b(i) + a(i,j)*x(j+k)
	end do
end do
return
end subroutine


end module
