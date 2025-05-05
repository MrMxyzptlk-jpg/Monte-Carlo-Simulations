MODULE funciones
    USE precision
    implicit none

    contains

real(pr) function sinusoid(x, constants)
    real(kind=pr), intent(in)    :: x
    real(kind=pr), intent(in)    :: constants(:)

    sinusoid = constants(1)*sin(x*constants(2)) + constants(3)*cos(x*constants(4))

end function

real(pr) function sinusoid_integral(x_lowerLim, x_upperLim, constants)
    real(kind=pr), intent(in)    :: x_lowerLim, x_upperLim
    real(kind=pr), intent(in)    :: constants(:)

    sinusoid_integral = - constants(1)*cos(x_upperLim*constants(2))/constants(2) &
                        + constants(3)*sin(x_upperLim*constants(4))/constants(4) &
                        + constants(1)*cos(x_lowerLim*constants(2))/constants(2) &
                        - constants(3)*sin(x_lowerLim*constants(4))/constants(4)

end function

END module