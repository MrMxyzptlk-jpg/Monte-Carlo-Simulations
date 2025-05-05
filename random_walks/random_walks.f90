!***************************************************************
! Program: random_walks.f90
! Purpose: random walk analysis
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program random_walks
    use precision
    use subrutinas
    use formats
    use mzranmod
    implicit none

    integer(kind=int_large)                 :: j, k, dim, walkers, steps, min_steps, max_steps, walkers_saved
    integer(kind=int_large)                 :: unit_endPosition, unit_moments, unit_trips
    real(kind=pr)                           :: mean, stddev
    real(kind=pr), allocatable              :: r_avg(:), r_var(:), r_mean_diff(:,:)
    integer(kind=int_large), allocatable    :: walkers_position(:,:)
    logical                                 :: save_walks
    character(len=:), allocatable           :: file_endPosition, file_moments, file_trips, prefix, suffix


!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /calculation/ dim, walkers, min_steps, max_steps, walkers_saved, save_walks

    ! DEFAULT SETTINGS
    dim = 3
    walkers = 1e5
    min_steps = 1
    max_steps = 100
    walkers_saved = 0
    save_walks = .false.

    ! Read from input file
    open(newunit=unit_endPosition, file="input.nml", status="old", action="read")
        read(unit_endPosition, nml=calculation)
    close(unit_endPosition)
!###################################################################################################

    prefix = "./datos/steps_"
    suffix = ".out"
    file_moments = "./datos/moments.out"

!##############################################################################################################################

    call mzran_init()

    allocate(walkers_position(walkers,dim), r_avg(dim),r_var(dim), r_mean_diff(walkers,dim))

    open(newunit=unit_moments, file=file_moments,status='replace')
        write(unit_moments, *) "## steps  | mean | stddev"
        do steps = min_steps, max_steps
            call create_file_name(prefix, steps, suffix, file_endPosition)
            walkers_position = 0

            if (save_walks) then
                call create_file_name("datos/trips_steps",steps,suffix,file_trips)
                open(newunit=unit_trips, file=file_trips,status='replace')
                    write(unit_trips, *) "## steps  | coordinates in canonical basis"
                    do j = 1, steps
                        call random_step(walkers_position)
                        do k = 1, walkers_saved
                            write(unit_trips, format_style1) steps, walkers_position(k,:)
                        end do
                    end do
                close(unit_trips)
                deallocate(file_trips)
            else
                do j = 1, steps
                    call random_step(walkers_position)
                end do
            end if

            open(newunit=unit_endPosition, file=file_endPosition,status='replace')
            do j = 1, walkers
                write(unit_endPosition,format_style1) walkers_position(j,:)
            end do
            close(unit_endPosition)
            deallocate(file_endPosition)

            r_avg = sum(real(walkers_position,pr),1)
            mean = sum(r_avg) / real(walkers, pr)
            r_mean_diff = real(walkers_position, pr) - mean
            r_var = sum(r_mean_diff*r_mean_diff,1) / real(walkers - 1, pr)
            stddev = sqrt(sum(r_var))
            write(unit_moments, *) steps, mean, stddev

        end do
    close(unit_moments)

end program random_walks