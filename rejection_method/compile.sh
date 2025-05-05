gfortran ./modulos/precision.f90\
        ./modulos/formats.f90\
        ./modulos/funciones.f90\
        ./modulos/mzranmod.f90\
        ./modulos/subrutinas.f90\
        rejection_method.f90 -o run.exe -ffpe-trap=invalid,overflow,zero -O -Wall -fcheck=all -g -fbacktrace
#./run.exe

echo
