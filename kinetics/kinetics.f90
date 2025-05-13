PROGRAM KMC
    use formats
    use subrutinas
    use mzranmod

    implicit none

    real(pr)                :: t, tau       ! time, time increment
    integer                 :: i, unitnum   ! counter, unit
    integer                 :: selected_process, iostat           ! selected process
    character(len=15)       :: equations    ! system of equations

    abstract interface
        subroutine init()
        end subroutine init

        subroutine process(selected_process)
            integer, intent(in)    :: selected_process
        end subroutine process
    end interface
    procedure (init), pointer :: initialize => null()
    procedure (process), pointer :: doProcess => null()

!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /physical/ volume
    namelist /calculation/ max_iterations
    namelist /tasks/ equations
    namelist /initial_conditions/ species_total, velocity_cte

    ! DEFAULT SETTINGS
        ! Physical problems' characteristics
        volume = 10._pr

        !Calculation settings
        max_iterations = 100

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=physical)
        read(unitnum, nml=calculation)
        read(unitnum, nml=tasks)

        ! Select the system of equations to read the initial conditions of said system. Note that Fortran is NOT case-sensitive in case select for character variables.
        select case(equations)
            case("Gillespie")
                initialize => Gillespie_init
                doProcess  => Gillespie_process

                num_species = 4
                num_eqs = 2

            case("Consecutive")
                initialize => consecutive_init
                doProcess  => consecutive_process

                num_species = 3
                num_eqs = 2

            case("Prey-Predator")
                initialize => preyPredator_init
                doProcess  => preyPredator_process

                num_species = 2
                num_eqs = 3
        end select

        ! This allocation MUST be done in order to read the initial conditions data
        allocate(species_total(num_species), velocity_cte(num_eqs), velocity(num_eqs))

        iostat = 0
        read(unitnum, nml=initial_conditions, iostat=iostat)
        if (iostat /= 0) then
            print *, "Error reading initial_conditions namelist. iostat =", iostat
            stop 1
        end if
    close(unitnum)

!##################################################################################################
!      Initializations
!##################################################################################################

    call mzran_init()

    inv_vol = 1._pr/volume
    inv_vol2 = 1._pr/(volume*volume)

    t=0._pr

!##################################################################################################
!      Start of the calculations
!##################################################################################################


    call initialize()
    print '(a,*(F8.4))', "species_total = ", species_total
    print '(a,*(F8.4))', "velocity_cte  =", velocity_cte


    open(newunit=unitnum, file=salida)
        write(unitnum,*) "##    t     | species"

        do i = 1, max_iterations
            call select_process(selected_process)
            call doProcess(selected_process)
            call time(tau)
            t = t + tau
            write(unitnum, format_style1) t, species_total
        end do

    close(unitnum)

END PROGRAM KMC