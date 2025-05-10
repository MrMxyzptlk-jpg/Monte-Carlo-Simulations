PROGRAM KMC
    use formats
    use subrutinas
    use mzranmod

    implicit none

    real(pr)                :: r            ! suma de velocidades
    real(pr)                :: t, tau       ! tiempo, incremento de tiempo
    integer                 :: i, unitnum, equations   ! contador, unidad, systema de eqs.
    integer                 :: nu           ! proceso seleccionado

    abstract interface
        subroutine init(r)
            use precision
            real(pr), intent(out)    :: r        ! suma de velocidades
        end subroutine init
    end interface
    procedure (init), pointer :: initialize => null()

!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /physical/ volume
    namelist /calculation/ max_iterations
    namelist /tasks/ equations

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
    close(unitnum)

!##################################################################################################
!      Initializations
!##################################################################################################



inv_vol = 1._pr/volume
inv_vol2 = 1._pr/(volume*volume)


t=0._pr

select case(equations)
    case(1)
        initialize => inicio_test
!       sistema a estudiar:
!       m=2, n=4
!          x(1)     --> x(2)    (k1)
!       x(2) + x(3) --> x(4)    (k2)

    case(2)
        initialize => inicio_test2
end select

call initialize(r)

!##################################################################################################
!      Start of the calculations
!##################################################################################################


    open(newunit=unitnum,file=salida)
        write(unitnum,*) "##    t     | species"

        do i = 1, max_iterations
            call select_process(r, nu)
            call tiempo(r, tau)
            call ejecuta_process(r, nu)
            t=t+tau
            write(unitnum, format_style1) t, x
        enddo

    close(unitnum)

END PROGRAM KMC