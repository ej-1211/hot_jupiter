module mod_uvzder 

  use mod_const, only : r8

  implicit none
  public :: uvzder
contains
  !----------------------------------------------------------------------------
  subroutine uvzder( is, maxl, l, iewr, n, eps, nt, phi, &
                     a, b, c, u, v, z,                   &
                     strmfn, velpot, dvrgnc, vrtcty,     &
                     grdntu, grdntv, grdntz, w       ) 
    implicit none
    integer,  intent(in)    :: is, maxl, l, iewr, n, nt
    real(r8), intent(in)    :: eps, phi(nt)
    real(r8), intent(inout) :: a(2*n), b(2*n), c(2*n)
    real(r8), intent(out)   :: u(nt), v(nt), z(nt)
    real(r8), intent(out)   :: strmfn(nt), velpot(nt), dvrgnc(nt)
    real(r8), intent(out)   :: vrtcty(nt), grdntu(nt), grdntv(nt)
    real(r8), intent(out)   :: grdntz(nt)
    real(r8), intent(inout) :: w(27*n+32)
    
    integer :: isym

    !write(*,*) '#############################################'
    !write(*,*) '####          subroutine uvzder          ####'
    !write(*,*) '#############################################'

    isym = mod(l,2)
    if(iewr==3) isym = 1 - isym

    call uvzdr1( is, n, isym, nt, eps, phi, a, b, c, u, v, z, &
                 strmfn, velpot, dvrgnc, vrtcty,              &
                 grdntu, grdntv, grdntz )
                 !grdntu, grdntv, grdntz, w, w(3*n+1:27*n+32)      )

  end subroutine uvzder
  !----------------------------------------------------------------------------
  subroutine uvzdr1(is, n, isym, nt, eps, phi, a, b, c, u, v, z, &
                 strmfn, velpot, dvrgnc, vrtcty,                 &
                 grdntu, grdntv, grdntz )
                 !grdntu, grdntv, grdntz, p, w      )
    use mod_const, only:  &
        pai, pis2 => hfpai

    implicit none
    integer,  intent(in)    :: is, n, isym, nt
    real(r8), intent(in)    :: eps, phi(nt)
    real(r8), intent(inout) :: a(2*n), b(2*n), c(2*n)
    real(r8), intent(out)   :: u(nt), v(nt), z(nt)
    real(r8), intent(out)   :: strmfn(nt), velpot(nt), dvrgnc(nt)
    real(r8), intent(out)   :: vrtcty(nt), grdntu(nt), grdntv(nt)
    real(r8), intent(out)   :: grdntz(nt)
    real(r8), allocatable   :: p(:)
    real(r8), allocatable   :: work(:)

    integer :: i
    integer :: j, j1
    integer :: npn, npnps

    real(r8), parameter :: conv=1.0d-17
    real(r8) :: amax, bmax, cmax
    real(r8) :: s, ts
    real(r8) :: c1, seps, th
    real(r8) :: ph, sinp, cosp
    real(r8) :: t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13
    real(r8) :: alp1, alp2, cfn, cpse, fn

    amax = 0.0_r8
    bmax = 0.0_r8
    cmax = 0.0_r8

    do i = 1, n
       amax = dmax1( amax, dabs(a(i)) )
       bmax = dmax1( bmax, dabs(b(i)) )
       cmax = dmax1( cmax, dabs(c(i)) )
    end do

    !write(*,*)  dabs(a(n)), dabs(a(n)) - conv * amax
    !write(*,*)  dabs(b(n)), dabs(b(n)) - conv * bmax
    !write(*,*)  dabs(c(n)), dabs(c(n)) - conv * cmax

    !if( ( dabs(a(n)) - conv * amax ) > 0.0_r8 .or. &
    !    ( dabs(b(n)) - conv * bmax ) > 0.0_r8     )  return
    !write(*,*) "#######"
    !write(*,*) abs(c(n)), amax ,abs(c(n)) - conv * amax
    !if( ( abs(c(n)) - conv * cmax ) < 0.0_r8 ) then
    !if( ( dabs(c(n)) - conv * cmax ) <= 0.0_r8 ) then

       npn = 2*n
       npnps = npn + is
       s     = dfloat(is)
       ts    = 2.0_r8 * s
       seps  = 1.0_r8 / dsqrt(eps)
       !write(*,*) 'seps',seps
       c1    = -0.5_r8 * seps
 
       allocate( p(npnps+1) )
       allocate( work(npnps*5) )

       call dlfnc( 0, is+1, npnps, th, p, work ) 


       do i = 1, nt
          ph   = phi(i)
          sinp = dsin(ph) 
          th   = pis2 - ph
          call dlfnc( 1, is+1, npnps, th, p, work )
             !write(*,*) "p"
          do j = 1, npn
             j1   = j + is
             p(j) = p(j1) 
             !write(*,*) p
          end do
          fn  = s - 2.0_r8
          t5  = 0.0_r8
          t6  = 0.0_r8
          t12 = 0.0_r8
 
          if(isym == 0) then
             do j = 1, n
                fn = fn + 2.0_r8
                alp1 = dsqrt( (fn + s + 1.0_r8) * (fn - s) )
                alp2 = dsqrt( (fn + s + 2.0_r8) * (fn - s + 1.0_r8) )
                t5   = t5  + alp1 * a(j) * p(2*j-1)
                t6   = t6  + alp2 * b(j) * p(2*j)
                t12  = t12 + alp1 * c(j) * p(2*j-1)
             end do
          else
             do j = 1, n
                fn = fn + 2.0_r8
                alp1 = dsqrt( (fn + s + 1.0_r8) * (fn - s) )
                alp2 = dsqrt( (fn + s + 2.0_r8) * (fn - s + 1.0_r8) )
                t5   = t5  + alp2 * a(j) * p(2*j)
                t6   = t6  + alp1 * b(j) * p(2*j-1)
                t12  = t12 + alp2 * c(j) * p(2*j)
             end do
          end if 
          grdntz(i) = t12

          if( is == 0 ) then
             !write(*,*) c1, t6, t5
             u(i) = c1 * ( t6 + t6 )
             v(i) = c1 * ( t5 + t5 )
          else
             u(i) = sinp * t5 + t6
             v(i) = sinp * t6 + t5
          end if

       end do

       call dlfnc( 0, is, npnps, th, p, work )

       do i = 1, nt
          ph   = phi(i)
          !sinp = dcos(ph)
          cosp = dcos(ph)
          th   = pis2 - ph
          call dlfnc( 1, is, npnps, th, p, work )
          do j = 1, npn
             j1   = j + is
             p(j) = p(j1)
          end do
          t3  = 0.0_r8
          t4  = 0.0_r8
          t7  = 0.0_r8
          t10 = 0.0_r8
          t11 = 0.0_r8
          fn  = s - 2.0_r8

          if( isym == 0 ) then
             do j = 1, n 
                t8 = a(j) * p(2*j-1)
                t9 = b(j) * p(2*j)
                t3 = t3 + t8
                t4 = t4 + t9
                fn = fn + 2.0_r8
                cfn = fn * ( fn + 1.0_r8 ) 
                t10 = t10 + cfn * t8
                t11 = t11 + ( cfn + fn + fn + 2.0_r8 ) * t9
                t7  = t7 + c(j) * p(2*j-1)
             end do
          else
             do j = 1, n 
                t8 = a(j) * p(2*j)
                t9 = b(j) * p(2*j-1)
                t3 = t3 + t8
                t4 = t4 + t9
                fn = fn + 2.0_r8
                cfn = fn * ( fn + 1.0_r8 ) 
                t10 = t10 + ( cfn + fn + fn + 2.0_r8 ) * t8
                t11 = t11 + cfn * t9
                t7  = t7 + c(j) * p(2*j)
             end do
          end if

          if( is == 0 ) then
             z(i) = t7
             velpot(i) = - t3
             strmfn(i) =   t4
             dvrgnc(i) =   t10
             vrtcty(i) = - t11
          else 
             u(i) = u(i) + ts * cosp * t3
             v(i) = v(i) + ts * cosp * t4
             z(i) = t7
             velpot(i) = - t3
             strmfn(i) =   t4
             dvrgnc(i) =   t10
             vrtcty(i) = - t11
          end if

       end do
 
       if( is == 0 ) then
          do i = 1, nt
             cpse = seps * dcos(phi(i))
             !write(*,*) 'cpse',cpse
             grdntu(i) = s * v(i) - cpse * vrtcty(i)
             grdntv(i) = s * u(i) + cpse * dvrgnc(i)
          end do
       else

          call dlfnc( 0, is-1, npnps, th, p, work )

          do i = 1, nt
             ph   = phi(i)
             sinp = dsin(ph)
             th   = pis2 - ph
             call dlfnc( 1, is-1, npnps, th, p, work )
             do j = 1, npn
                j1   = j + is
                p(j) = p(j1)
             end do
             t1   = 0.0_r8
             t2   = 0.0_r8
             t13  = 0.0_r8
             fn   = s - 2.0_r8
   
             if( isym == 0 ) then
                do j = 1, n
                   fn = fn + 2.0_r8
                   alp1 = dsqrt( (fn + s) * (fn - s + 1.0_r8) )
                   alp2 = dsqrt( (fn + s + 1.0_r8) * (fn - s + 2.0_r8) )
                   t1   = t1  + alp1 * a(j) * p(2*j-1)
                   t2   = t2  + alp2 * b(j) * p(2*j)
                   t13  = t13 + alp1 * c(j) * p(2*j-1)
                end do
             else
                do j = 1, n
                   fn = fn + 2.0_r8
                   alp1 = dsqrt( (fn + s) * (fn - s + 1.0_r8) )
                   alp2 = dsqrt( (fn + s + 1.0_r8) * (fn - s + 2.0_r8) )
                   t1   = t1  + alp2 * a(j) * p(2*j)
                   t2   = t2  + alp1 * b(j) * p(2*j-1)
                   t13  = t13 + alp2 * c(j) * p(2*j)
                end do
             end if
             u(i) = c1 * ( sinp * t1 - t2 + u(i) )
             v(i) = c1 * ( sinp * t2 - t1 + v(i) )
             grdntz(i) = 0.5_r8 * ( grdntz(i) - t13 )
          end do
          do i = 1, nt
             cpse = seps * dcos(phi(i))
             grdntu(i) = s * v(i) - cpse * vrtcty(i)
             grdntv(i) = s * u(i) + cpse * dvrgnc(i)
          end do
       end if
    !end if
  end subroutine uvzdr1
  !----------------------------------------------------------------------------
end module mod_uvzder 
