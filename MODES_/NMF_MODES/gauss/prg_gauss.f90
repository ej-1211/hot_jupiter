  program main
    use mod_gauss

    integer :: N
    real(8), allocatable :: x(:)
    real(8), allocatable :: w(:)
    character(128) :: gauss_fname
  
    namelist / gaussian / &
      N, gauss_fname
 
    open(1,file='gauss.cnf')
    read(1,gaussian)

    allocate( x(N) )
    allocate( w(N) )

    call gauss(N,x,w)

    do i = 1,N
      write(*,*) i, x(i), w(i)
    end do

    open(2,file=trim(gauss_fname),form='unformatted',access='direct',&
         recl=8*N*2)
    write(2,rec=1) x(:), w(:)
    close(2)

    write(*,*) sum(w)

  end program main
