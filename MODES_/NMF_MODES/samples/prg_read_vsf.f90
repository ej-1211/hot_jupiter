  program main
    implicit none
    integer :: i
    integer, parameter :: num_vmode=72
    integer, parameter :: mp=73
    double precision :: vsf(mp,num_vmode)

    open(1,file='~/NMF_DATABASE/vsf/vsf_nume_FD_sigma.data', &
         form='unformatted',access='direct',recl=8*num_vmode*mp)

    read(1,rec=1) vsf
    do i = 1,mp
      write(*,'(i8,3f15.8)')  i, vsf(i,1), vsf(i,2), vsf(i,3)
    end do

  end program main
