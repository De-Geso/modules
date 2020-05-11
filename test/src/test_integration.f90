! test_integration.f90
! Brad Friesen
! Updated 2020-04-26
! Test integration module
! Compile with make_integration.sh

program test_integration
use integration
implicit none

integer, parameter :: n = 1000
real, parameter :: a = -2
real, parameter :: b = 2
character(len=*), parameter :: path = 'data/integration/'
real x, xstep, sqtrap, sqsimp, sqromb, smidpnt, s1midinf, s2midinf
integer i


call init()

call open_all()

do i = 1, n
	x = a + xstep*(i-1.0)
	call qtrap(f, a, x, sqtrap)
	call qsimp(f, a, x, sqsimp)
	call qromb(f, a, x, sqromb)
	call qromo(f, a, x, smidpnt, midpnt)
	s1midinf = 0.0
	s2midinf = 0.0
	if (x .lt. 0.0) then
		call qromo(f, -2.0, x, s1midinf, midinf)
	else
		call qromo(f, -2.0, -0.5, s1midinf, midinf)
		call qromo(f, -0.5, x, s2midinf, midpnt)
	end if
	call write_all()	
end do


contains


subroutine init()
xstep = (b-a)/(n-1.0)
end subroutine

pure function f(x)
! Function that we want to integrate
intent(in) x
real f, x

! Exact answer is 94/25 = 3.76 for limits a=-2, b=2.
f = 1.0/5.0*(1.0/100.0*(322.0+3.0*x*(98.0+x*(37.0+x)))-24.0*x/(1.0+x**2))
end function

subroutine open_all()
open (10, file=path//'function.dat', status='replace')
open (20, file=path//'qtrap.dat', status='replace')
open (30, file=path//'qsimp.dat', status='replace')
open (40, file=path//'qromb.dat', status='replace')
open (50, file=path//'midpnt.dat', status='replace')
open (60, file=path//'midinf.dat', status='replace')
end subroutine

subroutine write_all()
write (10, *) x, f(x)
write (20, *) x, sqtrap
write (30, *) x, sqsimp
write (40, *) x, sqromb
write (50, *) x, smidpnt
write (60, *) x, s1midinf+s2midinf
end subroutine
end program
