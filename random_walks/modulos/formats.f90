module formats
    implicit none
    character (len=13)   :: format_style0 = "(*(E14.7,3x))"         ! 0 integers, all floats
    character (len=11)   :: format_style1 = "(*(I10,3x))"           ! all integers, 0 floats

END MODULE formats