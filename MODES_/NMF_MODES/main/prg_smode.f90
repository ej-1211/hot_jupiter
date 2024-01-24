  program main
    use mod_adm, only : &
        houghcalc_init,     &
        mk, maxl
    use mod_smode
    implicit none
    integer :: is
    integer :: nn
    integer :: ierr
    double precision :: eps
    double precision :: an(40)
    double precision :: h0(40, 18, 3)

    is   = 0
    eps  = 8.0d0

    call houghcalc_init

    call smtrx_smode(eps)
    call shige(eps, h0)

  end program main
