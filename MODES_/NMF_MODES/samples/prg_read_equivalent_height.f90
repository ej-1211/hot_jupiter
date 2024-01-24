  program main
    implicit none
    integer, parameter :: num_vmode=72
    double precision :: eveh(num_vmode)

    open(1,file='~/NMF_DATABASE/vsf/equivalent_height_nume_FD_sigma.dat', &
         form='unformatted',access='direct',recl=8*num_vmode)

    read(1,rec=1) eveh
    write(*,*) eveh

  end program main
