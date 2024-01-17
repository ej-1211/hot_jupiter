  program main
    implicit none
    integer, parameter :: mp=73
    double precision :: vgrid(mp)
    double precision :: vweight(mp)

    open(1,file='~/NMF_DATABASE/vsf/ecmwf_sigma_levels_right.dat', &
         form='unformatted',access='sequential')

    read(1) vgrid
    write(*,*) vgrid

  end program main
