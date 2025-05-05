!***************************************************************
! Program: uniform_sampling_integration.f90
! Purpose:
!
!
! Description:
!
! Input:
!
!
! Output:
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program uniform_sampling_integration
    use precision
    use formats
    use subrutinas
    use funciones
    use mzranmod
    implicit none

    real(kind=pr)             :: exponent, MC_integral, integral, lower_lim, upper_lim, stddev
    integer(kind=int_large)   :: i, samples, unitnum, max_sample_exponent
    logical                   :: rel_err = .false.

    abstract interface
        function funcion(x_x,constants)
            use precision
            implicit none
            real(kind=pr)                :: funcion
            real(kind=pr), intent(in)    :: x_x
            real(kind=pr), intent(in)    :: constants
        end function
    end interface

procedure (funcion), pointer :: function_pointer => null()

!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /calculation/ exponent, lower_lim, upper_lim, max_sample_exponent

    ! DEFAULT SETTINGS
        exponent = 3._pr
        lower_lim = 0._pr
        upper_lim = 1._pr
        max_sample_exponent = 1000

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=calculation)
    close(unitnum)

function_pointer => raise

!###################################################################################################
!   Start of the calculations
!###################################################################################################

    call mzran_init()
    integral = raise_integral(upper_lim, exponent) - raise_integral(lower_lim, exponent)
    if (integral /= 0._pr) rel_err = .true.

    open(newunit=unitnum, file="./datos/MC_integral_regular_sampling.out")
        if (rel_err) then
            write(unitnum,*) "## Samples | MC integral | Analytic integral | Abs error | stddev | Rel Error"
            do i = 1, max_sample_exponent
                samples = 10*i
                call MC_with_uncertainty(lower_lim, upper_lim, samples, function_pointer, exponent, MC_integral, stddev)
                write(unitnum,format_style1) samples, MC_integral, integral, abs(integral - MC_integral), stddev &
                , abs((integral - MC_integral)/integral)
            end do
        else
            write(unitnum,*) "## Samples | MC integral | Analytic integral | Abs error | stddev"
            do i = 1, max_sample_exponent
                samples = 10*i
                call MC_with_uncertainty(lower_lim, upper_lim, samples, function_pointer, exponent, MC_integral, stddev)
                write(unitnum,format_style1) samples, MC_integral, integral, abs(integral - MC_integral), stddev
            end do
        end if
    close(unitnum)

end program uniform_sampling_integration