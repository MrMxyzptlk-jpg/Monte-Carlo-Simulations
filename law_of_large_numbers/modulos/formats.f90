module formats
    implicit none

    character (len=13)   :: format_style0 = "(*(E14.7,3x))"             ! 0 integers, all floats
    character (len=23)   :: format_style1 = "(1(I10,3x),*(E14.7,3x))"   ! 1 integers, the rest floats

END MODULE formats
