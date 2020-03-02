! parameters.f90 by Brad Friesen
! Constants that may be needed. (ex. pi, planck's constant, etc.)
! Compile with: gfortran -O3 -fdefault-real-8 -c parameters.f90
! Updated: 2020-03-01

module parameters
use, intrinsic :: iso_fortran_env, only: stdin=>input_unit, &
	stdout=>output_unit, stderr=>error_unit

real(16), parameter :: PI = 4.0*atan(1.0_16)

end module parameters
