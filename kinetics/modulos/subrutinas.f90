MODULE subrutinas
    use mzranmod
    use precision

    integer                     :: num_species
    integer                     :: num_eqs
    integer                     :: max_iterations
    character(21), parameter    :: salida='datos/KMC-prueba1.dat'    ! output file
    real(pr)                    :: velocities_sum, volume, inv_vol, inv_vol2
    real(pr), allocatable       :: species_total(:)   ! species' net quantities
    real(pr), allocatable       :: velocity_cte(:)   ! velocity constants
    real(pr), allocatable       :: velocity(:)

contains

subroutine update_velocitySum()
    velocities_sum = sum(velocity)
    if (velocities_sum == 0) then
        print*, "Program stopping: no more possible progress. "
        STOP
    end if
end subroutine

subroutine Gillespie_init()

    velocity(1) = velocity_cte(1)*species_total(1)*inv_vol
    velocity(2) = velocity_cte(2)*species_total(2)*species_total(3)*inv_vol2

    call update_velocitySum()

end subroutine Gillespie_init

subroutine consecutive_init()

    velocity(1) = velocity_cte(1)*species_total(1)*inv_vol
    velocity(2) = velocity_cte(2)*species_total(2)*inv_vol

    call update_velocitySum()

end subroutine consecutive_init

subroutine preyPredator_init()

    velocity(1) = velocity_cte(1)*species_total(1)
    velocity(2) = velocity_cte(2)*species_total(2)*species_total(1)
    velocity(3) = velocity_cte(3)*species_total(2)

    call update_velocitySum()

end subroutine preyPredator_init

subroutine Gillespie_process(selected_process)
    integer, intent(in)       :: selected_process

    select case(selected_process)
        case(1)
            species_total(1) = species_total(1) - 1._pr
            species_total(2) = species_total(2) + 1._pr
        case(2)
            species_total(2) = species_total(2) - 1._pr
            species_total(3) = species_total(3) - 1._pr
            species_total(4) = species_total(4) + 1._pr
        case default
            write(*,*)'Wrong process'
    end select

    velocity(1) = velocity_cte(1)*species_total(1)*inv_vol
    velocity(2) = velocity_cte(2)*species_total(2)*species_total(3)*inv_vol2

    call update_velocitySum()

end subroutine Gillespie_process

subroutine consecutive_process(selected_process)
    integer, intent(in)       :: selected_process

    select case(selected_process)
        case(1)
            species_total(1) = species_total(1) - 1._pr
            species_total(2) = species_total(2) + 1._pr
        case(2)
            species_total(2) = species_total(2) - 1._pr
            species_total(3) = species_total(3) + 1._pr
        case default
            write(*,*)'Wrong process'
    end select

    velocity(1) = velocity_cte(1)*species_total(1)*inv_vol
    velocity(2) = velocity_cte(2)*species_total(2)*inv_vol

    call update_velocitySum()

end subroutine consecutive_process

subroutine preyPredator_process(selected_process)
    integer, intent(in)       :: selected_process

    select case(selected_process)
        case(1)
            species_total(1) = species_total(1) + 1._pr
        case(2)
            species_total(1) = species_total(1) - 1._pr
            species_total(2) = species_total(2) + 1._pr
        case(3)
            species_total(2) = species_total(2) - 1._pr
        case default
            write(*,*)'Wrong process'
    end select

    velocity(1) = velocity_cte(1)*species_total(1)
    velocity(2) = velocity_cte(2)*species_total(2)*species_total(1)
    velocity(3) = velocity_cte(3)*species_total(2)

    call update_velocitySum()

end subroutine preyPredator_process


subroutine select_process(selected_process)
    integer, intent(out)    :: selected_process
    real(pr)                :: c
    integer                 :: i
    real(pr)                :: sum

    c = rmzran()*velocities_sum

    sum = 0._pr

    ! This section could be done with SUM() and int() functions but is left like this for clarity
    do i = 1, num_eqs
        sum = sum + velocity(i)
        if (c < sum) then
            selected_process = i
            exit
        end if
    end do

end subroutine select_process

subroutine time(time_increment)
    real(pr), intent(out)    :: time_increment

    time_increment = -log(rmzran())/velocities_sum

end subroutine time


END MODULE subrutinas
