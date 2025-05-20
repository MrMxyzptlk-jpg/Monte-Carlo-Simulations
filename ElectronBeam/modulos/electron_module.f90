module electron_module
    use precision
    use formats
    use materials
    use mzranmod
    use constantes
    implicit none
    private
    public :: Electron, track_trajectory, init_angles, init_materials, init_energy_range, MC_step

    type :: Electron
        real(pr) :: x, y, z
        real(pr) :: x_old, y_old, z_old
        real(pr) :: v_x, v_y, v_z
        real(pr) :: r_x, r_y, r_z
        real(pr) :: energy
        real(pr) :: step_length
        contains
        procedure :: init => electron_init
        procedure :: move
        procedure :: do_step
        procedure :: lambda
        procedure :: sigma_e
        procedure :: alpha
        procedure :: step
        procedure :: cos_scat_angle
        procedure :: dE_dS
        procedure :: update
    end type Electron

    ! Bulk material constants
    real(pr)    :: Z, A, Rho

    ! Film cover constants
    real(pr)    :: Zf, Af, Rhof
    real(pr)    :: film_thickness

    ! Auxiliary constants
    real(pr)    :: Z67, Zf67, J, Jf

    ! Ray tilt
    real(pr)    :: cos_tilt(3)

    ! Energy range
    real(pr)    :: Emin, Emax



contains

subroutine init_energy_range(E_max, E_min)
    real(pr)    :: E_min, E_max
    Emax = E_max
    Emin = E_min
end subroutine init_energy_range

subroutine init_angles(ray_tilt)
    real(pr)    :: ray_tilt, comp_tilt
    comp_tilt = 90.0_pr - ray_tilt
    cos_tilt = [0.0_pr, cos(comp_tilt*pi/180.0_pr), cos(ray_tilt*pi/180.0_pr)]
end subroutine init_angles

subroutine init_materials(bulk_material, surface_material, covering, covering_thickness)
    character(len=*)    :: bulk_material, surface_material
    logical             :: covering
    real(pr)            :: covering_thickness
    type(Material)      :: mat

    mat = get_material(bulk_material)
    Z = mat%Z
    A = mat%A
    Rho = mat%rho
    Z67 = mat%Z**0.67_pr
    J = 1e-3_pr*(9.76_pr*Z + 58.5_pr*Z**(-0.19_pr))

    if (covering) then
        mat = get_material(surface_material)
        Zf = mat%Z
        Af = mat%A
        Rhof = mat%rho
        Zf67 = mat%Z**0.67_pr
        Jf = 1e-3_pr*(9.76_pr*mat%Z + 58.5_pr*mat%Z**(-0.19_pr))
        film_thickness = covering_thickness
    else
        mat = get_material(bulk_material)
        Zf = mat%Z
        Af = mat%A
        Rhof = mat%rho
        Zf67 = mat%Z**0.67_pr
        Jf = 1e-3_pr*(9.76_pr*mat%Z + 58.5_pr*mat%Z**(-0.19_pr))
        film_thickness = 0._pr
    end if

end subroutine init_materials

subroutine electron_init(particle, x0, y0, z0)
    class(Electron), intent(inout) :: particle
    real(pr), intent(in) :: x0, y0, z0

    particle%x = x0; particle%y = y0; particle%z = z0
    particle%v_x = cos_tilt(1); particle%v_y = cos_tilt(2); particle%v_z = cos_tilt(3)
    particle%energy = Emax
end subroutine electron_init

subroutine move(particle, Rnd_1, Rnd_2, Rnd_3)
    class(Electron), intent(inout) :: particle
    real(pr), intent(in) :: Rnd_1, Rnd_2, Rnd_3
    real(pr) :: theta, cos_theta, sin_theta, cos_phi, sin_phi
    real(pr) :: U1, U2, norm

    if (particle%z >= 0.0_pr) then
        particle%x_old = particle%x; particle%y_old = particle%y; particle%z_old = particle%z   ! Saving previous coordinates

        ! Calculating in the particle's coordinate system (velocity along z')
        particle%step_length = particle%step(Rnd_1)  ! Random step using exponential decay
        theta = 2._pr*pi*Rnd_2                       ! Random angle from uniform distribution
        cos_phi = particle%cos_scat_angle(Rnd_3)     ! Random azimuth angle from formula (see eq. 10 from Saenz. et al. 2016)

        sin_phi = sqrt(1.0_pr - cos_phi*cos_phi)    ! To avoid recalculation
        cos_theta = cos(theta)                      ! To avoid recalculation
        sin_theta = sin(theta)                      ! To avoid recalculation

        norm = sqrt(particle%v_x*particle%v_x + particle%v_z*particle%v_z)

        if (norm > tiny(1._pr)) then    ! Avoid unstable division
            ! Realative to y axis
            U1  = -particle%v_z * sin_phi / norm
            U2  = particle%v_x * sin_phi / norm

            ! Rotate back to original reference frame
            particle%r_x = particle%v_x*cos_phi + U1*cos_theta + particle%v_y*U2*sin_theta     ! New displacement component along x
            particle%r_y = particle%v_y*cos_phi + sin_theta*(particle%v_z*U1 - particle%v_x*U2) ! New displacement component along y
            particle%r_z = particle%v_z*cos_phi + U2*cos_theta - particle%v_y*U1*sin_theta     ! New displacement component along z

        else    ! Use alternate reference frame
            norm = sqrt(particle%v_y*particle%v_y + particle%v_z*particle%v_z)

            ! Realative to x axis
            U1  = -particle%v_z * sin_phi / norm
            U2  = particle%v_y * sin_phi / norm

            ! Rotate back to original reference frame
            particle%r_x = particle%v_x*cos_phi + sin_theta*(particle%v_z*U1 - particle%v_y*U2)
            particle%r_y = particle%v_y*cos_phi + U1*cos_theta + particle%v_x*U2*sin_theta
            particle%r_z = particle%v_z*cos_phi + U2*cos_theta - particle%v_x*U1*sin_theta

        end if

        particle%x = particle%x + particle%step_length*particle%r_x
        particle%y = particle%y + particle%step_length*particle%r_y
        particle%z = particle%z + particle%step_length*particle%r_z
    else
        particle%x = particle%x + (particle%x - particle%x_old)
        particle%y = particle%y + (particle%y - particle%y_old)
        particle%z = particle%z + (particle%z - particle%z_old)
    end if
end subroutine move

subroutine update(particle)
    class(Electron), intent(inout) :: particle
    if (particle%z < film_thickness .and. particle%z > 0.0_pr) then
        particle%energy = particle%energy - particle%step_length * Rhof * particle%dE_dS()
    else
        particle%energy = particle%energy - particle%step_length * Rho * particle%dE_dS()
    end if
    particle%v_x = particle%r_x; particle%v_y = particle%r_y; particle%v_z = particle%r_z
end subroutine update

subroutine do_step(particle, Rnd_1, Rnd_2, Rnd_3)
    class(Electron), intent(inout) :: particle
    real(pr), intent(in) :: Rnd_1, Rnd_2, Rnd_3
    call particle%move(Rnd_1, Rnd_2, Rnd_3)
    call particle%update()
end subroutine do_step

real(pr) function step(particle, Rnd_1)
    class(Electron), intent(in) :: particle
    real(pr), intent(in) :: Rnd_1
    step = -particle%lambda() * log(Rnd_1)
end function step

real(pr) function lambda(particle)
    class(Electron), intent(in) :: particle
    if (particle%z < film_thickness .and. particle%z > 0.0_pr) then
        lambda = Af / (NA * Rhof * particle%sigma_e())
    else
        lambda = A / (NA * Rho * particle%sigma_e())
    end if
end function lambda

real(pr) function sigma_e(particle)
    class(Electron), intent(in) :: particle
    real(pr) :: Zval, E
    real(pr) :: alpha_val

    if (particle%z < 2e-6_pr .and. particle%z > 0.0_pr) then
        Zval = Zf
        alpha_val = particle%alpha(Zf67)
    else
        Zval = Z
        alpha_val = particle%alpha(Z67)
    end if

    E = particle%energy
    sigma_e = 5.21e-21_pr * (Zval / E)**2 * (4._pr*pi / (alpha_val * (1._pr + alpha_val))) * ((E + 511.0_pr)/(E + 1024.0_pr))**2
end function sigma_e

real(pr) function alpha(particle, Zp67)
    class(Electron), intent(in) :: particle
    real(pr), intent(in) :: Zp67
    alpha = 3.4e-3_pr * Zp67 / particle%energy
end function alpha

function cos_scat_angle(particle, Rnd_3) result(cos_phi)
    class(Electron), intent(in) :: particle
    real(pr), intent(in) :: Rnd_3
    real(pr) :: cos_phi
    real(pr) :: alpha_val

    if (particle%z < film_thickness .and. particle%z > 0.0_pr) then
        alpha_val = particle%alpha(Zf67)
    else
        alpha_val = particle%alpha(Z67)
    end if

    cos_phi = 1.0_pr - (2.0_pr * alpha_val * Rnd_3) / (1.0_pr + alpha_val - Rnd_3)
end function cos_scat_angle

function dE_dS(particle) result(de)
    class(Electron), intent(in) :: particle
    real(pr) :: de, Jval, Zval, Aval
    if (particle%z < film_thickness .and. particle%z > 0.0_pr) then
        Zval = Zf; Aval = Af; Jval = Jf
    else
        Zval = Z; Aval = A; Jval = J
    end if
    de = 78500.0_pr * (Zval / (Aval * particle%energy)) * log(1.166_pr*(particle%energy + 0.85_pr*Jval)/Jval)
end function dE_dS

subroutine track_trajectory(unit, step, x, y, z, energy)
    implicit none
    integer, intent(in) :: unit, step
    real(pr), intent(in) :: x, y, z, energy
    write(unit, format_style_data) step, x, y, z, energy
end subroutine track_trajectory

subroutine MC_step(sample, i)
    type(Electron), intent(inout)   :: sample(:)
    integer, intent(in)             :: i
    real(pr)                        :: R1, R2, R3

    ! Walk randomly
    R1 = rmzran()
    R2 = rmzran()
    R3 = rmzran()
    call sample(i)%do_step(R1, R2, R3)

end subroutine MC_step

end module electron_module
