# Brad Friesen
# Updated: 2020-03-01
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

echo 'Building parameters...'
fbuild -Jbuild/ -o build/parameters.o -c ../parameters.f90
echo 'parameters built'

echo 'Building numerical_recipes...'
fbuild -Jbuild/ -o build/numerical_recipes.o -c ../numerical_recipes.f90
echo 'numerical_recipes built'

echo 'Building special_functions...'
fbuild -Jbuild/ -o build/special_functions.o -c ../special_functions.f90
echo 'special_functions built'

echo 'Building test_fresnel...'
gfortran -Ibuild/ -c code/test_fresnel.f90 -o build/test_fresnel.o
gfortran -o fresnel.out build/test_fresnel.o build/parameters.o build/special_functions.o
echo 'test_fresnel built'
