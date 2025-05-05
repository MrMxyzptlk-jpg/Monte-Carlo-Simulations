gfortran ./modulos/precision.f90\
        ./modulos/constantes.f90\
        ./modulos/formats.f90\
        ./modulos/mzranmod.f90\
        ./modulos/subrutinas.f90\
        law_of_large_numbers.f90 -o run.exe \
        -ffpe-trap=invalid,overflow,zero -O -Wall -fcheck=all -g -fbacktrace
#./run.exe

#gnuplot plot.gp

echo