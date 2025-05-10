gfortran ./modulos/precision.f90\
        ./modulos/formats.f90\
        ./modulos/mzranmod.f90\
        ./modulos/subrutinas.f90\
        ejemplo-Gillespie.f90 -o run.x -O3 \
        -march=native -ftree-vectorize -ftree-vectorizer-verbose=2 \
        -ffpe-trap=invalid,overflow,zero -O -Wall -fcheck=all -g -fbacktrace