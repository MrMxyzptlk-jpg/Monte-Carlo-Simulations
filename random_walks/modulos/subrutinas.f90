MODULE subrutinas
    USE precision
    use constantes
    use mzranmod
    implicit none

    contains

subroutine create_file_name(prefix, num, suffix, filename)
    character(len=8)                :: fmt  ! Format descriptor
    character(len=12)               :: x1   ! Temporary string for formatted real
    character(len=*), intent(in)    :: prefix, suffix
    character(len=:), allocatable   :: filename
    integer, intent(in)             :: num

    fmt = '(I8.8)'  ! adjust as needed
    write(x1, fmt) num  ! Convert to string

    ! Trim spaces in formatted number
    x1 = adjustl(trim(x1))

    ! Concatenate strings
    filename = prefix // trim(x1) // suffix

end subroutine create_file_name

subroutine random_step(random_array)
    integer(kind=int_large), allocatable    :: random_array(:,:)
    real(kind=pr)                           :: rnd_num
    integer                                 :: index, dim
    integer(kind=int_large)                 :: i

    dim = size(random_array,2)

    do i = 1, size(random_array,1)
        ! index = floor(rmzran()*real(dim)) + 1 ! Implementation if we consider 1 step ONLY
        do index = 1, dim
            rnd_num = rmzran()
            if (rnd_num<0.5) then
                random_array(i,index) = random_array(i,index) + 1
            elseif (rnd_num>=0.5) then
                random_array(i,index) = random_array(i,index) - 1
            end if
        end do
    end do

end subroutine random_step

END module subrutinas