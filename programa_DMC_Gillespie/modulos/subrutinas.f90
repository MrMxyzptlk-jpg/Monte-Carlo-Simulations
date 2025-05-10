MODULE subrutinas
    use mzranmod
    use precision

    integer                     :: num_species=4              ! (numero de especies)
    integer                     :: num_eqs=2              ! (numero de ecuaciones)
    integer                     :: max_iterations
    real(pr)                    :: volume
    character(21),parameter     :: salida='datos/KMC-prueba1.dat'    ! archivo de salida
    real(pr)                    :: inv_vol, inv_vol2
    real(pr), allocatable       :: x(:)   ! especies
    real(pr), allocatable       :: k(:)   ! constantes de velocidad
    real(pr), allocatable       :: a(:)   ! velocidades

contains

subroutine inicio_test(r)
    real(pr), intent(out)    :: r        ! suma de velocidades

    num_eqs = 2
    num_species = 4
    allocate(x(num_species), a(num_eqs),k(num_eqs))


    x = 0._pr

    x(1) = 100._pr
    x(3) = 20._pr

    k = 1._pr

    a(1) = k(1)*x(1)*inv_vol
    a(2) = k(2)*x(2)*x(3)*inv_vol2

    r=0._pr

    r = r + sum(a)

    call mzran_init()

end subroutine inicio_test

subroutine inicio_test2(r)
    real(pr), intent(out)    :: r        ! suma de velocidades

    num_eqs = 2
    num_species = 4
    allocate(x(num_species), a(num_eqs),k(num_eqs))


    x = 0._pr

    x(1) = 100._pr
    x(3) = 20._pr

    k = 1._pr

    a(1) = k(1)*x(1)*inv_vol
    a(2) = k(2)*x(2)*x(3)*inv_vol2

    r=0._pr

    r = r + sum(a)

    call mzran_init()

end subroutine inicio_test2

subroutine ejecuta_process(r,nu)
    integer, intent(in)       :: nu         ! proceso selecionado
    real(pr),intent(inout)    :: r          ! suma de velocidades

    select case(nu)
        case(1)
            x(1)=x(1)-1._pr
            x(2)=x(2)+1._pr
        case(2)
            x(2)=x(2)-1._pr
            x(3)=x(3)-1._pr
            x(4)=x(4)+1._pr
        case default
            write(*,*)'proceso equivocado'
    end select

    a(1)=k(1)*x(1)*inv_vol
    a(2)=k(2)*x(2)*x(3)*inv_vol2

    r=0._pr

    r = r + sum(a)

end subroutine ejecuta_process


subroutine select_process(r,nu)
    real(pr),intent(in)     :: r        ! suma de velocidades
    integer,intent(out)     :: nu       ! proceso seleccionado
    real(pr)                :: c        ! numero aleatorio entre 0 y r
    integer                 :: i        ! contador
    real(pr)                :: suma

    c = rmzran()*r

    suma = 0._pr

    do i = 1, num_eqs
        suma=suma+a(i)

        if(c<suma)then
            nu=i
            exit
        endif

    enddo

end subroutine select_process

subroutine tiempo(r,tau)
    real(pr),intent(in)     :: r        ! suma de velocidades
    real(pr),intent(out)    :: tau      ! incremento de tiempo

    tau = -log(rmzran())/r

end subroutine    tiempo


END MODULE subrutinas
