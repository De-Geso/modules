# Brad Friesen
# Updated: 2020-04-27
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

echo 'Building interpolation...'
fbuild -Jbin/ -o bin/interpolation.o -c ../interpolation.f90
echo 'interpolation built'

echo 'Building integration...'
fbuild -Jbin/ -o bin/integration.o -c ../integration.f90
echo 'integration built'

echo 'Building test_integration...'
fbuild -Ibin/ -c src/test_integration.f90 -o bin/test_integration.o
fbuild -o bin/test_integration.out bin/interpolation.o bin/integration.o bin/test_integration.o
echo 'test_integration built'
