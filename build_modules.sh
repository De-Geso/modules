# author: Brad Friesen
#!/bin/bash

fbuild() {
	gfortran -O3 -fdefault-real-8 $@
}

program=$1

echo 'Building numerical_recipes...'
fbuild -c numerical_recipes.f90
echo 'numerical_recipes built'

echo 'Compiling' "$program"'...'
fbuild -c -o script.o ../$program
fbuild -o script.out numerical_recipes.o script.o

echo 'Running' "$program"'...'
./script.out
echo 'Removing intermediate files...'
#rm *.o *.out *.mod
