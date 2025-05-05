MODULE funciones
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32, real64
    USE precision
    implicit none

    contains

    real(pr) function raise(x, exponent)
        real(kind=pr), intent(in)    :: x
        real(kind=pr), intent(in)    :: exponent

        raise = x**exponent

    end function

    real(pr) function raise_integral(x, exponent)
        real(kind=pr), intent(in)    :: x
        real(kind=pr), intent(in)    :: exponent

        raise_integral = x**(exponent+1._pr)/(exponent+1._pr)

    end function


END module