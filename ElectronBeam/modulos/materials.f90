module materials
    use precision
    implicit none

    type :: Material
        real(pr) :: Z      ! Atomic number
        real(pr) :: A      ! Atomic mass [g/mol]
        real(pr) :: rho    ! Density [g/cm^3]
        real(pr) :: J      ! Ionization potential [MeV]
    end type Material

    private
    public :: Material, get_material

contains

    function get_material(name) result(mat)
        character(len=*), intent(in) :: name
        type(Material) :: mat

        select case (trim(adjustl(name)))
        case ("PMMA")
            mat = Material(Z=4.38_pr, A=39.97_pr, rho=1.18_pr, J=68.5e-6_pr)
        case ("Au", "Gold")
            mat = Material(Z=79.0_pr, A=196.97_pr, rho=19.3_pr, J=790.0e-6_pr)
        case ("C")
            mat = Material(Z=6.0_pr, A=12.011_pr, rho=2.0_pr, J=81.0e-6_pr)
        case ("Fe")
            mat = Material(Z=26.0_pr, A=55.845_pr, rho=7.874_pr, J=286.0e-6_pr)
        case ("Na")
            mat = Material(Z=11.0_pr, A=23.0_pr, rho=0.970_pr, J=790.0e-6_pr)
        case ("Si")
            mat = Material(Z=14.0_pr, A=28.085_pr, rho=2.33_pr, J=790.0e-6_pr)
        case ("K")
            mat = Material(Z=19.0_pr, A=39.1_pr, rho=0.862_pr, J=790.0e-6_pr)
        case ("V")
            mat = Material(Z=23.0_pr, A=50.9415_pr, rho=6.11_pr, J=790.0e-6_pr)
        case ("Cu")
            mat = Material(Z=29.0_pr, A=63.546_pr, rho=7.6_pr, J=790.0e-6_pr)
        case ("As")
            mat = Material(Z=33.0_pr, A=74.92159_pr, rho=5.776_pr, J=790.0e-6_pr)
        case ("Nb")
            mat = Material(Z=41.0_pr, A=92.90637_pr, rho=8.57_pr, J=790.0e-6_pr)
        case ("Sn")
            mat = Material(Z=50.0_pr, A=118.710_pr, rho=17.287_pr, J=790.0e-6_pr)
        case ("Cs")
            mat = Material(Z=55.0_pr, A=132.9054_pr, rho=1.873_pr, J=790.0e-6_pr)
        case ("Eu")
            mat = Material(Z=63.0_pr, A=151.964_pr, rho=5.243_pr, J=790.0e-6_pr)
        case ("Yb")
            mat = Material(Z=70.0_pr, A=173.045_pr, rho=6.965_pr, J=790.0e-6_pr)
        case ("Hg")
            mat = Material(Z=80.0_pr, A=200.592_pr, rho=13.53_pr, J=790.0e-6_pr)
        case default
            print *, "Unknown material: ", trim(name)
            mat = Material(Z=0._pr, A=0._pr, rho=0._pr, J=0._pr)
        end select
    end function get_material

end module materials
