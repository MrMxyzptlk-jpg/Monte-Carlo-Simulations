gfortran ./modulos/precision.f90\
        ./modulos/constantes.f90\
        ./modulos/formats.f90\
        ./modulos/mzranmod.f90\
        ./modulos/funciones.f90\
        ./modulos/subrutinas.f90\
        uniform_sampling_integration.f90 -o run.exe  \
        -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=2 \
        -ffpe-trap=invalid,overflow,zero -O -Wall -fcheck=all -g -fbacktrace
#./run.exe

echo
