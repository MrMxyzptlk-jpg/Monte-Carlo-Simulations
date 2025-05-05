MODULE subrutinas
    USE precision
    use constantes
    use mzranmod
    implicit none
    
    contains

    subroutine create_file_name(prefix, num, suffix, filename)
        character(len=8)                :: fmt  ! Format descriptor
        character(len=12)                :: x1   ! Temporary string for formatted real
        character(len=:), allocatable    :: prefix, suffix, filename
        integer, intent(in)        :: num  ! Input real number
    
        fmt = '(I3.0)'  ! Format integer with 0 decimal places (adjust as needed)
        write(x1, fmt) num  ! Convert real to string
    
        ! Trim spaces in formatted number
        x1 = adjustl(trim(x1))
    
        ! Concatenate strings
        filename = prefix // trim(x1) // suffix
    
    end subroutine create_file_name
    


END module subrutinas