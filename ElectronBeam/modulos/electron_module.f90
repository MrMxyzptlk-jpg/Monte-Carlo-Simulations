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

subroutine electron_init(this, x0, y0, z0)
    class(Electron), intent(inout) :: this
    real(pr), intent(in) :: x0, y0, z0

    this%x = x0; this%y = y0; this%z = z0
    this%cx = cos_tilt(1); this%cy = cos_tilt(2); this%cz = cos_tilt(3)
    this%energy = Emax
end subroutine electron_init

subroutine move(this, RanA, RanB, RanC)
    class(Electron), intent(inout) :: this
    real(pr), intent(in) :: RanA, RanB, RanC
    real(pr) :: psi, c_phi, AM, AN, V1, V2, V3, V4

    if (this%z >= 0.0_pr) then
        this%x_old = this%x; this%y_old = this%y; this%z_old = this%z
        this%step_length = this%step(RanA)
        psi = 2*pi*RanB
        c_phi = this%cos_scat_angle(RanC)
        AM = -(this%cx / this%cz)
        AN = 1.0_pr / sqrt(1.0_pr + AM*AM)
        V1 = AN * sqrt(1.0_pr - c_phi*c_phi)
        V2 = AN * AM * sqrt(1.0_pr - c_phi*c_phi)
        V3 = cos(psi)
        V4 = sin(psi)

        this%ca = this%cx*c_phi + V1*V3 + this%cy*V2*V4
        this%cb = this%cy*c_phi + V4*(this%cz*V1 - this%cx*V2)
        this%cc = this%cz*c_phi + V2*V3 - this%cy*V1*V4

        this%x = this%x + this%step_length*this%ca
        this%y = this%y + this%step_length*this%cb
        this%z = this%z + this%step_length*this%cc
    else
        this%x = this%x + (this%x - this%x_old)
        this%y = this%y + (this%y - this%y_old)
        this%z = this%z + (this%z - this%z_old)
    end if
end subroutine move

subroutine do_step(this, RanA, RanB, RanC)
    class(Electron), intent(inout) :: this
    real(pr), intent(in) :: RanA, RanB, RanC
    call this%move(RanA, RanB, RanC)
    call this%update()
end subroutine do_step

function step(this, RanA) result(s)
    class(Electron), intent(in) :: this
    real(pr), intent(in) :: RanA
    real(pr) :: s
    s = -this%lambda() * log(RanA)
end function step

function lambda(this) result(lmbd)
    class(Electron), intent(in) :: this
    real(pr) :: lmbd
    if (this%z < 2e-6_pr .and. this%z > 0.0_pr) then
        lmbd = Af / (NA * Rhof * this%sigma_e())
    else
        lmbd = A / (NA * Rho * this%sigma_e())
    end if
end function lambda

function sigma_e(this) result(sigma)
    class(Electron), intent(in) :: this
    real(pr) :: sigma, Zval, E
    real(pr) :: alpha_val

    if (this%z < 2e-6_pr .and. this%z > 0.0_pr) then
        Zval = Zf
        alpha_val = this%alpha(Zf67)
    else
        Zval = Z
        alpha_val = this%alpha(Z67)
    end if

    E = this%energy
    sigma = 5.21e-21_pr * (Zval / E)**2 * (4._pr*pi / (alpha_val * (1 + alpha_val))) * ((E + 511.0_pr)/(E + 1024.0_pr))**2
end function sigma_e

function alpha(this, Zp67) result(a)
    class(Electron), intent(in) :: this
    real(pr), intent(in) :: Zp67
    real(pr) :: a
    a = 3.4e-3_pr * Zp67 / this%energy
end function alpha

function cos_scat_angle(this, RanC) result(cos_phi)
    class(Electron), intent(in) :: this
    real(pr), intent(in) :: RanC
    real(pr) :: cos_phi
    real(pr) :: alpha_val

    if (this%z < film_thickness .and. this%z > 0.0_pr) then
        alpha_val = this%alpha(Zf67)
    else
        alpha_val = this%alpha(Z67)
    end if

    cos_phi = 1.0_pr - (2.0_pr * alpha_val * RanC) / (1.0_pr + alpha_val - RanC)
end function cos_scat_angle

function dE_dS(this) result(de)
    class(Electron), intent(in) :: this
    real(pr) :: de, Jval, Zval, Aval
    if (this%z < film_thickness .and. this%z > 0.0_pr) then
        Zval = Zf; Aval = Af; Jval = Jf
    else
        Zval = Z; Aval = A; Jval = J
    end if
    de = 78500.0_pr * (Zval / (Aval * this%energy)) * log(1.166_pr*(this%energy + 0.85_pr*Jval)/Jval)
end function dE_dS

subroutine update(this)
    class(Electron), intent(inout) :: this
    if (this%z < film_thickness .and. this%z > 0.0_pr) then
        this%energy = this%energy - this%step_length * Rhof * this%dE_dS()
    else
        this%energy = this%energy - this%step_length * Rho * this%dE_dS()
    end if
    this%cx = this%ca; this%cy = this%cb; this%cz = this%cc
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
