gfortran modulos/constantes.f90 \
        modulos/formats.f90 \
        modulos/precision.f90 \
        modulos/materials.f90 \
        modulos/mzranmod.f90 \
        modulos/electron_module.f90 \
        electron_scatter.f90 -o run.exe \
        -O3 -march=native -ftree-vectorize -ftree-vectorizer-verbose=2 \
        -ffpe-trap=invalid,overflow,zero -O -Wall -fcheck=all -g -fbacktrace

#./run.exe

echo