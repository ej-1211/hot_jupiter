module mod_abcoef

  use mod_const, only : r8

  implicit none
  private

  public :: abcoef
  private :: abcof1

contains
  !----------------------------------------------------------------------------
  subroutine abcoef(is, maxl, l, iewr, sig, eps, n, a, b, c, w)
    implicit none
    integer, intent(in)    :: is
    integer, intent(in)    :: maxl
    integer, intent(in)    :: l
    integer, intent(in)    :: iewr
    integer, intent(in)    :: n
    real(r8), intent(inout) :: a(n), b(n), c(n), w(5*n)
    real(r8), intent(in) :: sig, eps

    integer :: lp1, ir, ie, if, ig

    lp1 = l + 1
    ir  = n + 1
    ie  = ir + n
    if  = ie + n
    ig  = if + n
    call abcof1( is, maxl, lp1, iewr, sig, eps, a, b, c, n, w, &
                     w(ir:ir+n-1), w(ie:ie+n-1), w(if:if+n-1), w(ig:ig+n-1) )


  end subroutine abcoef
  !----------------------------------------------------------------------------
  subroutine abcof1( is, maxl, l, iewr, sig, eps, a, b, c, n, w, r, e, f, g)
    use mod_sigma, only: &
        kpqr
    implicit none
    integer, intent(in)    :: is
    integer, intent(in)    :: maxl
    integer, intent(in)    :: l
    integer, intent(in)    :: iewr
    integer, intent(in)    :: n
    real(r8), intent(inout) :: a(n), b(n), c(n), w(5*n)
    real(r8), intent(inout) :: r(n), e(n), f(n), g(n)
    real(r8), intent(in) :: sig, eps
  
    real(r8) :: qhl, fn, tn, ph, qh, qe
    real(r8) :: t1, t2
    real(r8) :: pivot, pvct
    real(r8) :: summ, cmax, absc, ss
    integer :: lmax, lmxm
    integer :: i1
    integer :: nt
    integer :: lc
   
    ! loop variable
    integer :: i , ii
    integer :: ib , if, ip2
    integer :: ns2
    integer :: isym
   
    ns2 = n / 2
    a(:) = 0.0_r8
    b(:) = 0.0_r8
    c(:) = 0.0_r8

    isym = mod(l-1,2)
    if( iewr == 3 ) isym = 1 - isym

    if( iewr == 3 .and. is == 0 ) then

       if( (l-2) < 0 ) then
         b(1) = 0.0_r8
         return
       elseif( (l-2) == 0 ) then
         b(1) = 1.0_r8 / dsqrt(2.0_r8 / eps + 1.0_r8 / 15.0_r8)
         c(2) = b(1) / dsqrt(15.0_r8)
         return
       end if

       if( isym == 0 ) then
          e(1) = 1.0_r8 / dsqrt(2.0_r8 / eps + 1.0_r8 / 15.0_r8) 
          qhl  = 1.0_r8 / dsqrt(15.0_r8)
          fn   = 1.0_r8
          tn   = 2.0_r8
          lmax = l / 2
       else
          e(1) = 1.0_r8
          qhl  = 0.0_r8
          fn   = 0.0_r8
          tn   = 0.0_r8
          lmax = (l+1) / 2
       end if

       do i = 2, lmax
          fn = fn + 2.0_r8
          tn = tn + 4.0_r8
          ph = ( fn + 1.0_r8 ) / dsqrt( (tn - 1.0_r8) * (tn + 1.0_r8) )
          qh = fn / dsqrt( (tn + 1.0_r8) * (tn + 3.0_r8) )
          qe = qhl * e(i-1)
          e(i) = 1.0_r8 / dsqrt( fn * ( fn + 1.0_r8 ) / eps + qh**2 + ph**2 &
                                * ( 1.0_r8 + qe ) * ( 1.0_r8 - qe ) )
          f(i) = - ph * qe * e(i-1)
          qhl = qh
       end do

       b(lmax) = e(lmax)
       do ii = 2, lmax
          i = lmax - ii + 1
          b(i) = f(i+1) * b(i+1)
       end do
       b(lmax+1) = 0.0_r8

       if( isym == 0 ) then
          fn   = 0.0_r8
          i1   = 1
       else
          fn   = -1.0_r8
          i1   = 0
       end if

       do i = 1, lmax
          fn = fn + 2.0_r8
          tn = fn * 2.0_r8
          i1 = i1 + 1
          t1 = ( fn - 1.0_r8 ) / dsqrt( (tn - 1.0_r8) * (tn + 1.0_r8) )
          t2 = ( fn + 2.0_r8 ) / dsqrt( (tn + 1.0_r8) * (tn + 3.0_r8) )
          c(i1) = t1 * b(i) + t2 * b(i+1)
       end do

       if( c(lmax) > 0 ) return

       lmxm = lmax + 1
       do i = 1, lmxm
          b(i) = -b(i)
          c(i) = -c(i)
       end do

     else
        if( is == 0 .and. l <= 1 ) then
          if(iewr==1) a(1) = 1.0_r8 
          if(iewr==2) c(1) = 1.0_r8 
          return
        end if

        call kpqr( n, is, eps, w, e, f, r )

        do i = 2, n
           a(i) = f(i-1)
           c(i-1) = e(i)
        end do
        a(1) = 0.0_r8
        c(n) = 0.0_r8

        if( isym == 0 ) then
           do i = 1, n, 2
              b(i) = w(i) - sig - r(i) / sig
              b(i+1) = w(i+1) - sig
           end do
        else 
           do i = 1, n, 2
              b(i) = w(i) - sig
              b(i+1) = w(i+1) - sig - r(i+1) / sig
           end do
        end if
        !if( isym == 0 ) then
        !   do i = 1, n, 2
        !      b(i) = w(i) - sig
        !      b(i+1) = w(i+1) - sig - r(i+1) / sig
        !   end do
        !else 
        !   do i = 1, n, 2
        !      b(i) = w(i) - sig - r(i) / sig
        !      b(i+1) = w(i+1) - sig
        !   end do
        !end if
        i = l
        nt = n - l
        call tripf( l, n, a, b, c, e, f, g )
        call tripb( nt, a(l+1:n), b(l+1:n), c(l+1:n), &
                    e(l+1:n), f(l+1:n), g(l+1:n) )
        !write(*,*) e,f,g

        !if( abs(e(i)) < abs(f(i+1)) )  then
        if( dabs(e(i)) >= dabs(f(i+1)) )  then
           pivot = e(i+1) - f(i+1) / e(i) * f(i)
           pvct  = dabs(e(i+1))
           e(i+1) = 1.0_r8
           e(i)  = -f(i) / e(i)
        else
           pivot = f(i) - e(i)/ f(i+1) * e(i+1)
           e(i)  = -e(i+1) / f(i+1)
           e(i+1) = 1.0_r8
           pvct  = dabs(f(i))
        end if
     
        do ii = 2, i
           ib = i - ii
           e(ib+1) = - (f(ib+1)*e(ib+2)+g(ib+1)*e(ib+3))/e(ib+1)
        end do

        ip2 = i + 2
        do if = ip2, n
           e(if)=-(f(if)*e(if-1)+g(if)*e(if-2))/e(if)
        end do

        a(:) = 0.0_r8 
        b(:) = 0.0_r8 
        c(:) = 0.0_r8 

        if( isym == 0 ) then
           do i = 1, ns2
              a(i) = e(2*i-1)
              b(i) = e(2*i)
              c(i) = r(2*i-1) / sig * e(2*i-1)
           end do
        else
           do i = 1, ns2
              a(i) = e(2*i)
              b(i) = e(2*i-1)
              c(i) = r(2*i) / sig * e(2*i)
           end do
        end if
 
        fn = dfloat(is-2)
        summ = 0.0_r8
        lc = 1
        cmax = 0.0_r8
        if(isym == 0) then
           do i = 1, ns2
              absc = dabs(c(i))
              if(absc > cmax) then
                 cmax = absc
                 lc = i
              end if
              fn = fn + 2.0_r8
              summ = summ + ( fn + 1.0_r8 ) * &
                            ( fn * a(i)**2 + ( fn + 2.0_r8 ) * b(i)**2 ) / &
                            eps + c(i)**2
           end do
        else
           do i = 1, ns2
              absc = dabs(c(i))
              if(absc > cmax) then
                 cmax = absc
                 lc = i
              end if
              fn = fn + 2.0_r8
              summ = summ + ( fn + 1.0_r8 ) * &
                            ( fn * b(i)**2 + ( fn + 2.0_r8 ) * a(i)**2 ) / &
                            eps + c(i)**2
           end do
        end if
        if(summ == 0.0_r8) return
     end if

     ss = dsign( 1.0_r8 / dsqrt(summ), c(lc) )

     do i = 1, ns2
        a(i) = ss * a(i)
        b(i) = ss * b(i)
        c(i) = ss * c(i)
        !write(*,*) i, ss, a(i), b(i), c(i)
        !write(*,*) i, ss, a(i), b(i), c(i)
     end do

  end subroutine abcof1
  !----------------------------------------------------------------------------
  subroutine tripf( l, n, a, b, c, e, f, g )
    implicit none
    integer, intent(in) :: l
    integer, intent(in) :: n
    real(r8), intent(inout) :: a(n), b(n), c(n)
    real(r8), intent(inout) :: e(n), f(n), g(n)

    integer :: j
    real(r8) :: ah, bh, factor

    e(1) = b(1)
    f(1) = c(1)
    g(1) = 0.0_r8

    if(l==1) return

    do j = 2, l
       if( dabs(e(j-1)) <= dabs( a(j) ) ) then
          ah = e(j-1)
          bh = f(j-1)
          e(j-1) = a(j)
          f(j-1) = b(j)
          g(j-1) = c(j)
          factor = -ah / e(j-1)
          e(j) = bh + factor * f(j-1)
          f(j) = factor * g(j-1)
          g(j) = 0.0_r8
       else
          e(j) = b(j) - a(j) / e(j-1) * f(j-1)
          f(j) = c(j)
          g(j) = 0.0_r8
       end if
    end do

  end subroutine tripf
  !----------------------------------------------------------------------------
  subroutine tripb( m, a, b, c, s, t, u )
    implicit none
    integer, intent(in) :: m
    real(r8), intent(inout) :: a(m), b(m), c(m)
    real(r8), intent(inout) :: s(m), t(m), u(m)

    integer :: j, jj
    real(r8) :: bh, ch, factor


    s(m) = b(m)
    t(m) = c(m)
    u(m) = 0.0_r8

    do jj = 2, m
       j = m - jj + 1
       if( dabs(s(j+1)) < dabs( c(j) ) ) then
          bh = t(j+1)
          ch = s(j+1)
          s(j+1) = c(j)
          t(j+1) = b(j)
          u(j+1) = a(j)
          factor = -ch / s(j+1)
          s(j) = bh + factor * t(j+1)
          t(j) = factor * u(j+1)
          u(j) = 0.0_r8
       else
          s(j) = b(j) - c(j) / s(j+1) * t(j+1)
          t(j) = a(j)
          u(j) = 0.0_r8
       end if
    end do
  end subroutine tripb
  !----------------------------------------------------------------------------
end module mod_abcoef
