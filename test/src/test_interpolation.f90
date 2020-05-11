! test_interpolation.f90
! Brad Friesen
! Updated 2020-04-26
! Test interpolation
! Compile with make_interpolation.sh

program test_interpolation
use interpolation
implicit none

! keep n <= nmax in interpolation
integer, parameter :: n = 5
integer, parameter :: nx = 1000
character(len=*), parameter :: path = 'data/interpolation/'
integer, allocatable :: seed(:)
real x(n), y(n)
real xx, xstep
real ypol, dypol
real yrat, dyrat
integer i, a

call init()

open (20, file=path//'interpolation.dat', status='replace')
do i = 1, n
	write (20, *) x(i), y(i)
end do
close (20)

open (30, file=path//'polint.dat', status='replace')
open (40, file=path//'ratint.dat', status='replace')
do i = 1, nx
	xx = x(1) + xstep*(i-1.0)
	call polint(x, y, n, xx, ypol, dypol)
	call ratint(x, y, n, xx, yrat, dyrat)
	write (30, *) xx, ypol, dypol
	write (40, *) xx, yrat, dyrat
end do
close (30)
close (40)


contains


subroutine init()
call random_seed()
call random_seed(size = a)
allocate (seed(a))
call random_seed(get = seed(1:a))
open (10, file=path//'seed', status='replace')
write (10,*) seed
close (10)

do i = 1,n
	x(i) = (10.0)/(n-1.0) * (i-1.0)
	call random_number (y(i))
end do

xstep = (x(n)-x(1))/(nx-1.0)
end subroutine

end program
