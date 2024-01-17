  program main
    implicit none
    integer :: i
    integer, parameter :: mp=46
    double precision :: grid(mp)
    double precision :: weight(mp)

    open(1,file='~/NMF_DATABASE/gauss/gauss46.dat', &
         form='unformatted',access='direct',recl=8*mp*2)

    read(1,rec=1) grid, weight
    do i = 1,mp
      write(*,'(i8,2f20.16)') i,grid(i), weight(i)
    end do 

  end program main
