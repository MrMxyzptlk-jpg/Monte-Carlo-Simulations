module formats
    implicit none
    character (len=43)   :: format_style_header  = '(a,2x,a,5x,a,6x,a,6x,a,6x,a,6x,a,6x,a,6x,a)'
    character (len=40)   :: format_style_data    = '(I5,2X,E13.6,2X,E13.6,2X,E13.6,2X,E13.6)'
END MODULE formats