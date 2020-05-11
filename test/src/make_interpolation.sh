# Brad Friesen
# Updated: 2020-04-26
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

echo 'Building interpolation...'
fbuild -Jbin/ -o bin/interpolation.o -c ../interpolation.f90
echo 'interpolation built'

echo 'Building test_interpolation...'
fbuild -Ibin/ -c src/test_interpolation.f90 -o bin/test_interpolation.o
fbuild -o bin/test_interpolation.out bin/interpolation.o bin/test_interpolation.o
echo 'test_interpolation built'
