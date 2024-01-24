module mod_sigma

  use mod_const, only : r8

  public  :: sigma
  private :: sigma1
  public  :: kpqr

contains
  !----------------------------------------------------------------------------
  subroutine sigma(is, maxl, nn, ierr, eps, eastgs, westgs, rotats, w)
        
    implicit none

    integer, intent(in)  :: is
    integer, intent(in)  :: maxl
    integer, intent(in)  :: nn
    integer, intent(inout) :: ierr
    integer :: ip, iq, ir, ia, ib, ic, ie, if, ig, ih
    !integer :: i
    integer :: n

    real(r8), intent(in)  :: eps
    real(r8), intent(out) :: eastgs(0:maxl)
    real(r8), intent(out) :: westgs(0:maxl)
    real(r8), intent(out) :: rotats(0:maxl)
    real(r8), intent(inout)  :: w(11*nn)

    write(*,*) '#############################################'
    write(*,*) '#### Computation of eigenfrequency start ####'
    write(*,*) '#############################################'

    !n = 3 * max0( 20,maxl,int(dsqrt(eps)) )
    n = nn
    write(*,*) 'n=',11*nn
    ip = n + 1
    iq = ip + n
    ir = iq + n
    ia = ir + n
    ib = ia + n
    ic = ib + n
    ie = ic + n
    if = ie + n
    ig = if + n
    ih = ig + n

    call sigma1(is, maxl, ierr, eps, eastgs, westgs, rotats, n, w, &
                w(ip:ip+n-1), w(iq:iq+n-1), w(ir:ir+n-1), &
                w(ia:ia+n-1), w(ib:ib+n-1), w(ic:ic+n-1), &
                w(ie:ie+n-1), w(if:if+n-1), w(ig:ig+n-1), w(ih:ih+n-1))

  end subroutine sigma
  !----------------------------------------------------------------------------
  subroutine sigma1(is, maxl, ierr, eps, eastgs, westgs, rotats, &
                    n, w, p, q, r, a, b, c, e, f, g, e2 )

    implicit none
    integer, intent(in)    :: is
    integer, intent(in)    :: maxl
    integer, intent(inout) :: n
    integer, intent(inout) :: ierr
    real(r8) :: p(n), q(n), r(n), a(n), b(n), c(n), e(n), f(n), g(n), e2(n) ! TJH FIXME intent
    real(r8), intent(in)  :: eps
    real(r8), intent(inout) :: eastgs(0:maxl)
    real(r8), intent(inout) :: westgs(0:maxl)
    real(r8), intent(inout) :: rotats(0:maxl)
    real(r8), intent(inout)  :: w(11*n)

    integer :: isym
    integer :: i1
    integer :: i, it
    real(r8) :: sig
    real(r8) :: deast
    real(r8) :: dwest
    real(r8) :: drotat

    real(r8), parameter :: error = 1.0d-12

    ierr = 0

    if( is.eq.0 ) then

      !call kpqr(n, is, eps, w, p, q, r) 
      call kpqr(n, is, eps, w(1:n), p, q, r) 
      call sigest(n, n, maxl, ierr, w, p, q, r, a, e, f, e2, &
                  eastgs, westgs, rotats ) 
      n = 2 * n / 3
      do i = 2, n
        a(i)   = q(i-1)
        c(i-1) = p(i)
      end do
      a(1) = 0.0_r8
      c(n) = 0.0_r8
      isym = 0

      do i = 2, maxl
        isym = 1 - isym
        sig  = eastgs(i)
        call sigcom(it, i, isym, sig, eastgs(i), n, w, r, a, b, c, g)
      end do

      do i = 2, maxl
        deast = eastgs(i) - eastgs(i-1)
        if(deast < 0 ) ierr = 8
        if( dabs(deast) < error * dabs( eastgs(i) ) ) ierr = 9
      end do

      eastgs(1) = 0.0_r8
      westgs(1) = 0.0_r8
      rotats(1) = 0.0_r8
      do i = 2, maxl
         westgs(i) = -eastgs(i)
         rotats(i) = 0.0_r8
      end do
    else

      !call kpqr(n, is, eps, w, p, q, r)
      call kpqr(n, is, eps, w(1:n), p, q, r)
      call sigest(n, n, maxl, ierr, w, p, q, r, a, e, f, e2, &
                  eastgs, westgs, rotats ) 
      n = 2 * n / 3
      do i = 2, n
        a(i)   = q(i-1)
        c(i-1) = p(i)
      end do
      a(1) = 0.0_r8
      c(n) = 0.0_r8
!
!     COMPUTE EASTWARD GRAVITY WAVE MODES
!
      isym = 1
      do i1 = 1, maxl
        isym = 1 - isym
        sig  = eastgs(i1)
        call sigcom(it, i1, isym, sig, eastgs(i1), n, w, r, a, b, c, g)
      end do

      do i = 2, maxl
        deast = eastgs(i) - eastgs(i-1)
        if( deast <= 0.0_r8 ) ierr = 1
        if( dabs(deast) < error * dabs(eastgs(i) ) ) ierr = 2
      end do
!
!     COMPUTE WESTWARD GRAVITY MODES
!
      isym = 1
      do i1 = 1, maxl
        isym = 1 - isym
        sig  = westgs(i1)
        call sigcom(it,  i1, isym, sig, westgs(i1), &                  
                    n, w, r, a, b, c, g)
      end do

      do i = 2, maxl
        dwest = westgs(i) - westgs(i-1)
        if( dwest >= 0.0_r8 ) ierr = 2
        if( dabs(dwest) < error * dabs(westgs(i) ) ) ierr = 4
      end do
!
!     COMPUTE ROTATIONAL MODES
!
      isym = 0
      do i1 = 1, maxl
        isym = 1 - isym
        sig  = rotats(i1)
        call sigcom(it,  i1, isym, sig, rotats(i1), &
                    n, w, r, a, b, c, g)
      end do

      do i = 2, maxl
        drotat = rotats(i) - rotats(i-1)
        if( drotat <= 0.0_r8 ) ierr = 5
        if( dabs(drotat) < error * dabs(rotats(i) ) ) ierr = 6
      end do
      if(rotats(1) <= westgs(1) ) ierr = 7

    end if

    write(*,*) 'Eigenfrequency of eastward gravity mode'
    write(*,'(5f16.10)') eastgs(1:maxl)
    write(*,*) 'Eigenfrequency of westward gravity mode'
    write(*,'(5f16.10)') westgs(1:maxl)
    write(*,*) 'Eigenfrequency of rotational mode'
    write(*,'(5f16.10)') rotats(1:maxl)

    !do i = 1, 11*n
    !   write(*,*) i, w(i)
    !end do

  end subroutine sigma1
  !----------------------------------------------------------------------------
  subroutine kpqr(n, is, eps, w, p, q, r)

    implicit none
    integer, intent(in)  :: n
    integer, intent(in)  :: is
    real(r8), intent(in)  :: eps

    real(r8), intent(out) :: p(n)
    real(r8), intent(out) :: q(n)
    real(r8), intent(out) :: r(n)
    !real(r8), intent(inout)  :: w(11*n)
    !real(r8), allocatable, intent(inout)  :: w(11*n)
    real(r8), intent(inout)  :: w(n)
    
    real(r8) :: fn, fnp1, tn, sqns, s
    real(r8) :: fnpsp1, fnmsp1, fnnp1, fnsnp1
    integer :: i

    if( is == 0 ) then
      w(1) = 0.0_r8
      p(1)  = 0.0_r8
      q(1)  = 0.0_r8
      r(1)  = 0.0_r8
      fn    = 0.0_r8
      fnp1  = 1.0_r8
      tn    = 0.0_r8
      sqns  = dsqrt(3.0_r8)
      do i = 2, n
        fn   = fn + 1.0_r8
        fnp1 = fnp1 + 1.0_r8
        tn   = tn + 2.0_r8
        w(i) = 0.0_r8
        p(i) = fnp1 / sqns
        sqns = dsqrt( (tn+1.0_r8)*(tn+3.0_r8) )
        q(i) = fn / sqns
        r(i) = -fn * fnp1 / eps
      end do
    else
      s      = dfloat(is)
      fn     = s - 1.0_r8       ! s-1
      fnp1   = fn + 1.0_r8      ! s
      fnpsp1 = s + s           ! s*2
      fnmsp1 = 0.0_r8
      tn     = fnpsp1 - 2.0_r8  ! 2*s - 2
      sqns   = dsqrt( fnmsp1 * fnpsp1 ) ! 0
      do i = 1, n
        fn     = fn + 1.0_r8       ! n + s - 1
        fnp1   = fnp1 + 1.0_r8     ! n + s
        fnpsp1 = fnpsp1 + 1.0_r8   ! n + 2*s
        fnmsp1 = fnmsp1 + 1.0_r8   ! n
        tn     = tn + 2.0_r8       ! 2*n + 2*s - 2
        fnnp1  = fn * fnp1        ! (n + s - 1) * (n + s)
        fnsnp1 = fn / fnp1        ! (n + s - 1) / (n + s)
        w(i)   = -s / fnnp1       ! -s / (n + s) / (n + s - 1)
        p(i)   = sqns / fnsnp1    ! 
        sqns   = dsqrt( fnmsp1 * fnpsp1 / ( (tn+1.0_r8)*(tn+3.0_r8) ) )
        q(i)   = fnsnp1 * sqns
        r(i)   = - fnnp1 / eps    ! - ( n + s - 1 ) * ( n + s ) / eps
      end do
    endif

  end subroutine kpqr
  !----------------------------------------------------------------------------
  subroutine sigest(n, ndim, maxl, ierr, w, p, q, r, a, d, e, e2, &
                    eastgs, westgs, rotats )
    implicit none
    integer, intent(in)    :: maxl
    integer, intent(in)    :: n
    integer, intent(in)    :: ndim
    integer, intent(inout) :: ierr
    !real(r8) :: p(n), q(n), r(n), a(n), d(n), e(n), e2(n)
    real(r8) :: p(n), q(n), r(n), a(n), b(n), c(n), d(n), e(n), e2(n)
    real(r8), intent(inout) :: eastgs(0:maxl)
    real(r8), intent(inout) :: westgs(0:maxl)
    real(r8), intent(inout) :: rotats(0:maxl)
    real(r8), intent(inout)  :: w(11*n)
    logical :: matz
    integer :: i1, i2, i3
    integer :: mls2, mdm, ns3

    integer :: i
    real(r8) :: tmp(n,3)
    real(r8) :: tmp1(3,n)
    !real(r8) :: z(ndim,n)
    real(r8) :: work(n)
    real(r8) :: work1(2*n-2)
    real(r8) :: qq(n,n)

    tmp(:,1) = a(:)
    tmp(:,2) = b(:)
    tmp(:,3) = c(:)

    matz=.false.

    call symtry(0, n, ndim, w, p, q, r, a, b, c)
    tmp1(1,:) = c(:)
    do i = 2,n
      tmp1(2,i-1) = b(i)
    end do
    do i = 3,n
      tmp1(3,i-2) = a(i)
    end do
    call dsbtrd('N','L',n, 2, tmp1, 3, d, e(1:n-1), qq, n, work, ierr)
    call dsteqr( 'N', n, d, e(1:n-1),qq ,n , work1, ierr) 

    mls2 = maxl / 2
    mdm  = maxl - mls2 - mls2
    ns3  = n / 3
    i1   = ns3 + 1
    i2   = ns3
    i3   = i2 + ns3

    do i = 1, mls2
      i1 = i1 - 1
      i2 = i2 + 1
      i3 = i3 + 1
      eastgs(2*i-1) = d(i3)
      westgs(2*i-1) = d(i1)
      rotats(2*i)   = d(i2)
    end do
    if(mdm==0) go to 1

    eastgs(maxl) = d(i3+1)
    westgs(maxl) = d(i1-1)

  1 call symtry(1, n, ndim, w, p, q, r, a, b, c)
    tmp1 = 0.0_r8
    tmp1(1,:) = c(:)
    do i = 2,n
      tmp1(2,i-1) = b(i)
    end do
    do i = 3,n
      tmp1(3,i-2) = a(i)
    end do
    call dsbtrd('N','L' ,n, 2, tmp1, 3, d, e(1:n-1), qq, n, work, ierr)
    call dsteqr( 'N', n, d, e(1:n-1),qq ,n , work1, ierr)
 
    i1 = ns3 + 1
    i2 = ns3
    i3 = i2  + ns3
    do i = 1, mls2
      i1 = i1 - 1
      i2 = i2 + 1
      i3 = i3 + 1
      eastgs(2*i)   = d(i3)
      westgs(2*i)   = d(i1)
      rotats(2*i-1) = d(i2)
    end do

    if(mdm /= 0) rotats(maxl) = d(i2+1)

  end subroutine sigest
  !----------------------------------------------------------------------------
  subroutine symtry( isym, n, ndim, w, p, q, r, a, b, c)
    implicit none
    integer, intent(in)    :: isym
    integer, intent(in)    :: n
    integer, intent(in)    :: ndim
    real(r8) :: p(n), q(n), r(n), a(n), b(n), c(n)
    real(r8), intent(inout)  :: w(11*n)

    integer :: j, i
    integer :: nm5

    a=0.0_r8
    b=0.0_r8
    c=0.0_r8

    j = -1

    nm5 = n-5

    do i = 1,nm5,3
      j = j + 2
      c(i)   = w(j)
      c(i+1) = w(j+1)
      b(i+1) = dsqrt( q(j)*p(j+1) )
      a(i+3) = dsqrt( q(j+1)*p(j+2) )
    end do

    j=j+2
    c(n-2) = w(j)
    c(n-1) = w(j+1)
    b(n-1) = dsqrt( q(j)*p(j+1) )

    j = -1
    if( isym.le.0 ) then
      do i = 1,n,3
        j=j+2
        a(i+2) = - dsqrt( -r(j) )
      end do
    else
      do i = 1,n,3
        j=j+2
        b(i+2) = - dsqrt( -r(j+1) )
      end do
    end if

  end subroutine symtry
  !----------------------------------------------------------------------------
  subroutine sigcom( it, lp1, isym, sig, sigc, n, w, r, a, b, c, g)
    implicit none
    integer, intent(inout) :: it
    integer, intent(in)    :: lp1
    integer, intent(in)    :: isym
    integer, intent(in)    :: n
    real(r8) :: r(n), a(n), b(n), c(n), g(n)
    real(r8), intent(inout)  :: w(11*n)
    real(r8) :: sh(3), dh(3)
    real(r8) :: sig, sigc, sigh, sigl, dsig, sigp
    real(r8) :: det, dd, detl
    real(r8) :: aq, bq, cq
    real(r8) :: r1, r2, r3
    real(r8) :: disc
    real(r8) :: error
    integer :: iflg

    error=1.0d-12
    sigh=sig
    it=0
    iflg=0
 10 det=detcom(lp1,n,isym,sig,w,r,a,b,c,g,error)
! 10 det=detcom(lp1,n,ierr,isym,sig,w,r,a,b,c,g,error)
    it=it+1
    if(it.gt.3) go to 8
    sh(it)=sig
    dh(it)=det
    if(it.eq.1) sig=.99999_r8*sigh
    if(it.eq.2) sig=1.00001_r8*sigh
    if(it.lt.3) go to 10
    detl=dh(1)
    sigl=sh(1)
    dsig=.00001_r8*sigh
    aq=0.5_r8*(dh(3)-2.0_r8*dh(1)+dh(2))/dsig**2
    bq=0.5_r8*(dh(3)-dh(2))/dsig
    cq=dh(1)
    r1=-cq/bq
    if(aq.eq.0.0_r8) go to 5
    disc=bq**2-4.0_r8*aq*cq
    if(disc.lt.0.0_r8) go to 5
    disc=dsqrt(disc)
  3 if(bq.gt.0.0_r8) go to 4
    r2=0.5_r8*(disc-bq)/aq
    r3=cq/aq/r2
    go to 6
  4 r2=-0.5_r8*(disc+bq)/aq
    r3=cq/aq/r2
  6 if(dabs(r3).lt.dabs(r2)) r2=r3
    if(dabs(r2).lt.dabs(r1)) r1=r2
  5 sig=sigl+r1
    sigh=sig
    go to 10
  8 dd=det-detl
    sigp=sig
    if(dd.eq.0.0_r8) go to 26
    sigp=sig-det*(sig-sigl)/dd
    if(iflg.ne.0) go to 26
    if(dabs(sigp-sig).lt.error*dabs(sig)) iflg=1
    sigl=sig
    sig=sigp
    detl=det
    go to 10
 26 sigc=sigp
!
  end subroutine sigcom
  !----------------------------------------------------------------------------
  real(r8) function detcom( lp1, n, isym, sig, &
                                    w, r, a, b, c, g, error )
    implicit none
    integer, intent(in) :: lp1, n
    real(r8), intent(in) :: r(n), a(n), c(n), g(n)
    real(r8), intent(inout) :: b(n)
    real(r8), intent(in) :: w(11*n)
    real(r8), intent(in) :: error, sig
    !integer :: ierr
    integer :: isym
    integer :: i, nt

    real(r8) :: e1, e2, f1, f2

    if( isym == 0 ) then

      do i = 1, n, 2
        b(i)   = w(i) - sig - r(i) / sig
        b(i+1) = w(i+1) - sig
      end do
      nt = n - lp1
      call tripvf(lp1, a, b, c, e1, f1)
      call tripvb(nt, a(lp1+1), b(lp1+1), c(lp1+1), e2, f2)
      detcom = e1 * e2 - f1 * f2

    else

      do i = 1, n, 2
        b(i) = w(i) - sig
        b(i+1) = w(i+1) - sig - r(i+1) / sig
      end do
      nt = n - lp1
      call tripvf(lp1, a, b, c, e1, f1)
      call tripvb(nt, a(lp1+1), b(lp1+1), c(lp1+1), e2, f2)
      detcom = e1 * e2 - f1 * f2

    end if

  end function detcom
  !----------------------------------------------------------------------------
  subroutine tripvf(n, a, b, c, e1, f1)
    implicit none
    integer, intent(in) :: n
    real(r8), intent(in)  :: a(n), b(n), c(n)
    real(r8), intent(out) :: e1, f1
    real(r8) :: g1, t1, ah, bh
    integer :: i

    e1 = b(1)
    f1 = c(1)
    g1 = 0.0_r8
    if( n == 1 ) return

    do i = 2, n
      if( dabs(e1) <= dabs(a(i) ) ) then
        ah = e1
        bh = f1
        e1 = a(i)
        f1 = b(i)
        g1 = c(i)
        t1 = - ah / e1
        e1 = bh + t1 * f1
        f1 = t1 * g1
        g1 = 0.0_r8
      else
        e1 = b(i) - a(i) / e1 * f1
        f1 = c(i) 
        g1 = 0.0_r8
      end if
    end do

  end subroutine tripvf
  !----------------------------------------------------------------------------
  subroutine tripvb(n, a, b, c, s, t)
    implicit none
    integer, intent(in) :: n
    real(r8), intent(in)  :: a(n), b(n), c(n)
    real(r8), intent(out) :: s, t
    real(r8) :: u, bh, ch, t1

    integer :: jj, j
 
    s = b(n)
    t = a(n)
    u = 0.0_r8

    if ( n == 1 ) return

    do jj = 2, n
      j = n - jj + 1
      if( dabs(s) < dabs( c(j) ) ) then
        bh = t 
        ch = s
         s = c(j)
         t = b(j)
         u = a(j)
        t1 = - ch / s
         s = bh + t1 * t
         t = t1 * u
      else
         s = b(j) - c(j) / s * t
         t = a(j)
         u = 0.0_r8
      end if
    end do

  end subroutine tripvb
  !----------------------------------------------------------------------------
end module mod_sigma
