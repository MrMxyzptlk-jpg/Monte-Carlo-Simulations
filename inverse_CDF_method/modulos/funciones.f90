MODULE funciones
    USE precision
    implicit none
    private
    public :: PDF_function

    contains

    ! Point density function
    function PDF_function(x, exponent)
        real(kind=pr)                :: PDF_function
        real(kind=pr), intent(in)    :: x
        real(kind=pr), intent(in)    :: exponent

        PDF_function = (exponent + 1._pr)*x**exponent

    end function

    ! Cumulative density function
    function CDF_function(x, exponent)
        real(kind=pr)                :: CDF_function
        real(kind=pr), intent(in)    :: x
        real(kind=pr), intent(in)    :: exponent

        CDF_function = x**exponent

    end function


END module