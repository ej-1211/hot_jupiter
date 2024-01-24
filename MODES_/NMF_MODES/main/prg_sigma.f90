  program main
    use mod_sigma
    implicit none
    integer :: is
    integer :: nn
    integer :: maxl
    integer :: ierr
    double precision :: eps

    double precision, allocatable :: eastgs(:)
    double precision, allocatable :: westgs(:)
    double precision, allocatable :: rotats(:)
    double precision :: w(60000)

    is   = 1
    maxl = 50
    eps  = 20.0d0

    allocate( eastgs(0:maxl) )
    allocate( westgs(0:maxl) )
    allocate( rotats(0:maxl) )

    call sigma(is, maxl, nn, ierr, eps, eastgs, westgs, rotats, w)

  end program main
