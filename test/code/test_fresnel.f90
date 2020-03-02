! test_fresnel.f90
! Brad Friesen
! Updated 2020-03-01
! Test the fresnel integrals from special functions.
! Compile with make_fresnel.sh


program test_fresnel
use parameters
use special_functions
implicit none

real, parameter :: xmin = -5.0
real, parameter :: xmax = 5.0
integer, parameter :: n = 10000
real, parameter :: dx = (xmax-xmin)/(n-1.0)

real(8) x, c, s
integer i

open (10, file='output/fresnel.dat', status='replace')
do i = 1, n
	x = xmin + dx * (i-1.0)
	call fresnel(x, c, s)
	write (10,*) x, c, s, PI
end do
close (10)

end program test_fresnel
