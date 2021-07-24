#! /bin/bash
#perl -ni -e 'print unless $. == 1' pre-program.txt
#echo '0' | cat - pre-program.txt > temp && mv temp pre-program.txt
clear
cd ../Output/
rm -rf Track_particle
rm -rf *.dat
mkdir Track_particle

cd -
cd ../Source

if [ -f "./dsmcStart" ]; 
then
rm -rf ./dsmcStart
echo "Removed ./dsmcStart"
fi

if [ -f "./dsmc" ]; 
then
rm -rf ./dsmc
echo "Removed ./dsmc"
fi

if [ -f "./a.out" ]; 
then
rm -rf ./a.out
echo "Removed ./a.out"
fi 

gcc -std=c99 -o dsmcStart dsmcStart.c -lm
./dsmcStart
#gcc -std=c99 -g DSMC.c dsmcFunctions.c -lm
#./a.out
gcc -std=c99 -o dsmc -g DSMC.c dsmcFunctions.c -lm
#icc -qopenmp  -DMKL_ILP64  -mkl=parallel  -L${MKLROOT}/lib/intel64 -liomp5 -lpthread  -ldl -o dsmc DSMC.c dsmcFunctions.c -lm
./dsmc
