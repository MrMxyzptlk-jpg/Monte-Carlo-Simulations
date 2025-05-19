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
        real(pr) :: cx, cy, cz
        real(pr) :: ca, cb, cc
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
    particle%cx = cos_tilt(1); particle%cy = cos_tilt(2); particle%cz = cos_tilt(3)
    particle%energy = Emax
end subroutine electron_init

subroutine move(particle, RanA, RanB, RanC)
    class(Electron), intent(inout) :: particle
    real(pr), intent(in) :: RanA, RanB, RanC
    real(pr) :: theta, cos_theta, sin_theta, cos_phi, sin_phi
    real(pr) :: xzDirection_x, xzDirection_z, direction_z, direction_x, xz_norm

    if (particle%z >= 0.0_pr) then
        particle%x_old = particle%x; particle%y_old = particle%y; particle%z_old = particle%z   ! Saving previous coordinates

        ! Calculating in the particle's reference frame (velocity along z')
        particle%step_length = particle%step(RanA)  ! Random step using exponential decay
        theta = 2._pr*pi*RanB                       ! Random angle from uniform distribution
        cos_phi = particle%cos_scat_angle(RanC)     ! Random azimuth angle from formula (see eq. 10 from Saenz. et al. 2016)

        sin_phi = sqrt(1.0_pr - cos_phi*cos_phi)    ! To avoid recalculation
        cos_theta = cos(theta)                      ! To avoid recalculation
        sin_theta = sin(theta)                      ! To avoid recalculation

        xz_norm = sqrt(particle%cx*particle%cx + particle%cz*particle%cz)

        if (xz_norm > tiny(1._pr)) then                 ! Avoid unstable division
            xzDirection_z  = particle%cz / xz_norm      ! z direction cosine in the x-z plane
            xzDirection_x  = particle%cx / xz_norm      ! x direction cosine in the x-z plane
            direction_z    = xzDirection_z * sin_phi    ! z direction cosine
            direction_x    = -xzDirection_x * sin_phi   ! x direction cosine
        else    ! As a first approximation
            direction_z = 0._pr     ! z direction cosine
            direction_x = 0._pr     ! x direction cosine
        end if

        ! Rotate back to original reference frame
        particle%ca = particle%cx*cos_phi + direction_z*cos_theta + particle%cy*direction_x*sin_theta     ! New displacement component along x
        particle%cb = particle%cy*cos_phi + sin_theta*(particle%cz*direction_z - particle%cx*direction_x) ! New displacement component along y
        particle%cc = particle%cz*cos_phi + direction_x*cos_theta - particle%cy*direction_z*sin_theta     ! New displacement component along z

        particle%x = particle%x + particle%step_length*particle%ca
        particle%y = particle%y + particle%step_length*particle%cb
        particle%z = particle%z + particle%step_length*particle%cc
    else
        particle%x = particle%x + (particle%x - particle%x_old)
        particle%y = particle%y + (particle%y - particle%y_old)
        particle%z = particle%z + (particle%z - particle%z_old)
    end if
end subroutine move

subroutine do_step(particle, RanA, RanB, RanC)
    class(Electron), intent(inout) :: particle
    real(pr), intent(in) :: RanA, RanB, RanC
    call particle%move(RanA, RanB, RanC)
    call particle%update()
end subroutine do_step

function step(particle, RanA) result(s)
    class(Electron), intent(in) :: particle
    real(pr), intent(in) :: RanA
    real(pr) :: s
    s = -particle%lambda() * log(RanA)
end function step

function lambda(particle) result(lmbd)
    class(Electron), intent(in) :: particle
    real(pr) :: lmbd
    if (particle%z < film_thickness .and. particle%z > 0.0_pr) then
        lmbd = Af / (NA * Rhof * particle%sigma_e())
    else
        lmbd = A / (NA * Rho * particle%sigma_e())
    end if
end function lambda

function sigma_e(particle) result(sigma)
    class(Electron), intent(in) :: particle
    real(pr) :: sigma, Zval, E
    real(pr) :: alpha_val

    if (particle%z < 2e-6_pr .and. particle%z > 0.0_pr) then
        Zval = Zf
        alpha_val = particle%alpha(Zf67)
    else
        Zval = Z
        alpha_val = particle%alpha(Z67)
    end if

    E = particle%energy
    sigma = 5.21e-21_pr * (Zval / E)**2 * (4._pr*pi / (alpha_val * (1._pr + alpha_val))) * ((E + 511.0_pr)/(E + 1024.0_pr))**2
end function sigma_e

function alpha(particle, Zp67) result(a)
    class(Electron), intent(in) :: particle
    real(pr), intent(in) :: Zp67
    real(pr) :: a
    a = 3.4e-3_pr * Zp67 / particle%energy
end function alpha

function cos_scat_angle(particle, RanC) result(cos_phi)
    class(Electron), intent(in) :: particle
    real(pr), intent(in) :: RanC
    real(pr) :: cos_phi
    real(pr) :: alpha_val

    if (particle%z < film_thickness .and. particle%z > 0.0_pr) then
        alpha_val = particle%alpha(Zf67)
    else
        alpha_val = particle%alpha(Z67)
    end if

    cos_phi = 1.0_pr - (2.0_pr * alpha_val * RanC) / (1.0_pr + alpha_val - RanC)
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

subroutine update(particle)
    class(Electron), intent(inout) :: particle
    if (particle%z < film_thickness .and. particle%z > 0.0_pr) then
        particle%energy = particle%energy - particle%step_length * Rhof * particle%dE_dS()
    else
        particle%energy = particle%energy - particle%step_length * Rho * particle%dE_dS()
    end if
    particle%cx = particle%ca; particle%cy = particle%cb; particle%cz = particle%cc
end subroutine update

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
