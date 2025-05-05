MODULE subrutinas
    USE precision
    use mzranmod
    implicit none

    contains

subroutine MC_with_uncertainty(lower_lim, upper_lim, N_samples, f, constants, estimate, stddev)
    integer, intent(in)          ::  N_samples
    real(pr), intent(in)         :: lower_lim, upper_lim
    real(pr), intent(out)        :: estimate, stddev
    real(pr), intent(out)        :: constants
    real(pr)                     :: fx, x
    real(pr)                     :: sum_w, sum_w2
    integer                      :: j

    interface
        real(pr) function f(x_x, constants)
            use precision
            real(pr), intent(in) :: x_x
            real(pr), intent(in) :: constants
        end function
    end interface


    sum_w = 0.0_pr
    sum_w2 = 0.0_pr

    do j = 1,  N_samples
        x = lower_lim + rmzran()*(upper_lim - lower_lim)
        fx = f(x, constants)
        sum_w  = sum_w  + fx
        sum_w2 = sum_w2 + fx*fx
    end do

    estimate = sum_w / real(N_samples, pr)
    stddev   = sqrt( (sum_w2 / real(N_samples, pr) - estimate**2) / real(N_samples, pr) )

end subroutine MC_with_uncertainty

END module subrutinas