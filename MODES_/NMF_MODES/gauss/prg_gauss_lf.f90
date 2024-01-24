  program main
    use mod_gauss_lf

    integer :: N
    real(16), allocatable :: x(:)
    real(16), allocatable :: w(:)
    real(8), allocatable :: xa(:)
    real(8), allocatable :: wa(:)
    character(128) :: gauss_fname
  
    namelist / gaussian / &
      N, gauss_fname
 
    open(1,file='gauss.cnf')
    read(1,gaussian)

    allocate( x(N) )
    allocate( w(N) )
    allocate( xa(N) )
    allocate( wa(N) )

    call gauss(N,x,w)

    do i = 1,N
    !do i = N,1,-1
      xa(i)=dble(x(i))
      wa(i)=dble(w(i))
      write(*,'(i5,2f25.20)') i, x(i), w(i)
    end do

    open(2,file=trim(gauss_fname),form='unformatted',access='direct',&
         recl=8*N*2)
    write(2,rec=1) xa,wa

  end program main
