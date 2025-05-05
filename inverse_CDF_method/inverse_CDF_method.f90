!***************************************************************
! Program: inverse_CDF_method.f90
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
! Room for improvement: function used could be optimized for the case x^3 by using f(x) = x*x*x
!
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program inverse_CDF_method
    use precision
    use formats
    use funciones
    use mzranmod
    implicit none

    real(kind=pr)           :: exponent, dx
    integer(kind=int_large) :: i, j, Integration_points, MC_samples
    integer                 :: unitnum
    real(pr), allocatable   :: x(:), pdf(:), cdf(:), inv_cdf_x(:), inv_cdf_u(:), samples(:), rnd_vals(:)

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
    namelist /calculation/ exponent, Integration_points, MC_samples

    ! DEFAULT SETTINGS
        exponent  = 3._pr
        Integration_points  = 1000
        MC_samples = 100000

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=calculation)
    close(unitnum)

function_pointer => PDF_function

!###################################################################################################
!   Start of the calculations
!###################################################################################################

    call mzran_init()
    allocate(x(Integration_points), pdf(Integration_points), cdf(Integration_points))
    dx = 1.0_pr / (Integration_points - 1)

    ! Build x and PDF
    do i = 1, Integration_points
        x(i) = (i - 1) * dx
        pdf(i) = function_pointer(x(i), exponent)
    end do

    ! Build CDF with trapezoidal rule
    cdf(1) = 0.0_pr
    do i = 2, Integration_points
        cdf(i) = cdf(i-1) + 0.5_pr * dx * (pdf(i-1) + pdf(i))
    end do

    ! Normalize (due to numerical integration error)
    cdf = cdf / cdf(Integration_points)

    ! Inverse CDF lookup table
    inv_cdf_u = cdf
    inv_cdf_x = x

    ! Sample from inverse CDF
    allocate(samples(MC_samples), rnd_vals(MC_samples))
    do i = 1, MC_samples
        rnd_vals(i) = rmzran()
        ! Find bracketing indices for interpolation
        do j = 2, Integration_points
            if (inv_cdf_u(j) >= rnd_vals(i)) exit
        end do
        ! Linear interpolation
        samples(i) = inv_cdf_x(j-1) + (rnd_vals(i) - inv_cdf_u(j-1)) * &
                    (inv_cdf_x(j) - inv_cdf_x(j-1)) / (inv_cdf_u(j) - inv_cdf_u(j-1))
    end do

    ! Output results
    open(newunit=unitnum, file="datos/samples.out")
        write(unitnum,*) "## x (uniform rng) | CDF(x)"
    do i = 1, MC_samples
        write(unitnum,format_style0) samples(i), rnd_vals(i)
    end do
    close(unitnum)

end program inverse_CDF_method