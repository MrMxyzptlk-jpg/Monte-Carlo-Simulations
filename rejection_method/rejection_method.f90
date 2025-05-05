!***************************************************************
! Program: rejection_method.f90
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
!
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program rejection_method
    use precision
    use formats
    use subrutinas
    use funciones
    use mzranmod
    implicit none

    real(kind=pr)             :: constants(4), integral, x_lowerLim, x_upperLim, y_lowerLim, y_upperLim, efficiency
    integer(kind=int_large)   :: samples, unitnum

    abstract interface
        function funcion(x_x,constants)
            use precision
            implicit none
            real(kind=pr)                :: funcion
            real(kind=pr), intent(in)    :: x_x
            real(kind=pr), intent(in)    :: constants(:)
        end function
    end interface

procedure (funcion), pointer :: function_pointer => null()

!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /calculation/ x_lowerLim, x_upperLim, y_lowerLim, y_upperLim, samples


    ! DEFAULT SETTINGS
        x_lowerLim = 0._pr
        x_upperLim = 1._pr
        y_lowerLim = 0._pr
        y_upperLim = 1._pr
        samples = 1000

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=calculation)
    close(unitnum)

function_pointer => sinusoid

!###################################################################################################
!   Start of the calculations
!###################################################################################################

    constants = (/1._pr,1._pr,1._pr,1._pr/)
    call mzran_init()

    call MC_rejection(x_lowerLim, x_upperLim, y_lowerLim, y_upperLim, function_pointer, constants, samples, integral, efficiency)
    print*, "samples | MC integral | analytic integral | Abs error | Efficiency"
    print format_style1, samples, integral, sinusoid_integral(x_lowerLim, x_upperLim, constants)&
        , abs(sinusoid_integral(x_lowerLim, x_upperLim, constants) - integral ), efficiency

end program rejection_method