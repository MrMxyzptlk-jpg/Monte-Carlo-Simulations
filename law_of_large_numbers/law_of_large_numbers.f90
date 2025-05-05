!***************************************************************
! Program: law_of_large_numbers.f90
! Purpose: calcular histogramas para la distribución de la suma de n términos aleatorios
!   (con n de 1 hasta n). Calcular los primeros 4 momentos de estas distribuciones y compararlas
!   con los de una distribución normal.
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program law_of_large_numbers
    use precision
    use subrutinas
    use formats
    use mzranmod
    implicit none

    integer(kind=int_large)                             :: i, j, k, N_variables, samples, unitnum
    real(kind=pr), dimension(:,:), allocatable          :: moments
    real(kind=pr), dimension(:), allocatable            :: random_vec, data_points
    character(len=:), allocatable                       :: filename, prefix, suffix, file_moments

!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /calculation/ N_variables, samples

    ! DEFAULT SETTINGS
    N_variables = 100
    samples = 100000

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=calculation)
    close(unitnum)
!###################################################################################################

    prefix = "./datos/variables_"
    suffix = ".out"
    file_moments = "./datos/moments.out"
!##############################################################################################################################

    call mzran_init()

    allocate(data_points(samples), moments(N_variables,4))
    do i = 1, N_variables
        allocate(random_vec(i))
        call create_file_name(prefix, i, suffix, filename)
        open(newunit=unitnum, file=filename,status='replace')
            do j = 1, samples
                random_vec = (/(rmzran(), k = 1, i)/)
                data_points(j) = sum(random_vec)/real(size(random_vec),pr)
                write(unitnum,format_style0) data_points(j)
            end do
        close(unitnum)
        deallocate(filename)

        moments(i,1) = sum(data_points)/real(samples,pr) ! Mean
        moments(i,2) = sum((data_points - moments(i,1))**2)/real(samples-1,pr) ! Variance
        moments(i,3) = sum((data_points - moments(i,1))**3)/real(samples-2,pr) / moments(i,2)**(3._pr/2._pr) ! Skewness
        moments(i,4) = sum((data_points - moments(i,1))**4)/real(samples-3,pr) / moments(i,2)**2 - 3._pr ! Excess Kurtosis

        deallocate(random_vec)
    end do

    open(newunit=unitnum, file=file_moments,status='replace')
        do i = 1, N_variables
            write(unitnum,format_style0) moments(i,1), moments(i,2), moments(i,3), moments(i,4)
        end do
    close(unitnum)

end program law_of_large_numbers