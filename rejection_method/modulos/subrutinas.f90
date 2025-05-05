MODULE subrutinas
    USE precision
    use formats
    use mzranmod
    implicit none

    contains

subroutine MC_rejection(x_lowerLim, x_upperLim, y_lowerLim, y_upperLim, f, constants, N_samples, integral, efficiency)
    use precision
    real(kind=pr), intent(in)                    :: x_lowerLim, x_upperLim, y_lowerLim, y_upperLim, constants(:)
    integer(kind=int_large), intent(in)          :: N_samples
    real(kind=pr), intent(out)                   :: integral, efficiency

    real(kind=pr)                                :: x, y, fx
    integer(kind=int_large)                      :: i, accepted_positive, accepted_negative
    real(kind=pr)                                :: scale_x, scale_y
    integer(kind=int_small)                      :: unitnum

    interface
        function f(xx, constants)
            use precision
            real(kind=pr)                :: f
            real(kind=pr), intent(in)    :: xx
            real(kind=pr), intent(in)    :: constants(:)
        end function
    end interface

    scale_x = x_upperLim - x_lowerLim
    scale_y = y_upperLim - y_lowerLim
    accepted_positive = 0
    accepted_negative = 0
    open(newunit=unitnum, file="datos/MC_histogram.out", status="replace")
        write(unitnum,*) "## sign(y) | x"
        do i = 1, N_samples
            x = x_lowerLim + rmzran() * scale_x
            y = y_lowerLim + rmzran() * scale_y
            fx = f(x, constants)

            ! Count if point (x, y) is under the curve (considering sign of fx)
            if (fx >= 0._pr .and. y <= fx .and. y >= 0._pr) then
                accepted_positive = accepted_positive + 1
                write(unitnum,format_style1) 1, x
            else if (fx < 0._pr .and. y >= fx .and. y <= 0._pr) then
                accepted_negative = accepted_negative + 1
                write(unitnum,format_style1) -1, x
            end if
        end do
    close (unitnum)

    ! Estimate area
    efficiency = real(accepted_positive + accepted_negative, pr) / real(N_samples, pr)
    integral = real(accepted_positive - accepted_negative, pr) / real(N_samples, pr) * scale_x * scale_y

end subroutine MC_rejection


END module subrutinas