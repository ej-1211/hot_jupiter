module mod_gauss

  use mod_const, only : r8

  implicit none

  contains


  subroutine gauss(N, x, w)

    integer,  intent(in)    :: N
    real(r8), intent(inout) :: x(N)
    real(r8), intent(inout) :: w(N)

    real(r8) :: a, b
    real(r8) :: lambak
    real(r8) :: pi
    integer  :: i, m, k

    pi = dacos(-1.0_r8)
    
    if(mod(N,2)==0) then
      m = N/2
    else
      m = (N-1)/2
    end if 

    do i = 1,m
      a = dsin(pi*dfloat(N-2*i)/dfloat(2*N)) 
      b = dsin(pi*dfloat(N+2-2*i)/dfloat(2*N)) 
      x(N+1-i) = NIBUN(N,a,b)
      x(i) = -x(N+1-i)
    end do

    if(mod(N,2)==1) x(m) = 0.0_r8

!$omp parallel
!$omp do
    do i = 1,N
      w(i) = 0.0_r8
      do k = 0,N
        lambak = 2.0_r8 / ( 2.0_r8*dfloat(k)+1.0_r8)
        w(i) = w(i) + P_k(k,x(i))**2 / lambak
      end do
      w(i) = 1.0_r8 / w(i)
    end do
!$omp end do
!$omp end parallel

  end subroutine gauss


  function nibun(N,a,b)
    integer,  intent(in)    :: N
    real(r8), intent(inout) :: a
    real(r8), intent(inout) :: b
    real(r8)                :: nibun

    real(r8), parameter :: delta = 1.0d-15
    real(r8), parameter :: epsi  = 1.0d-15

    integer  :: k
    real(r8) :: c

    do k = 1, 100000
      c = (a+b)/2.0_r8
      if( dabs(P_k(N,c)) < delta .or. dabs(a-b) < epsi ) then
        nibun = c
        return
      end if 
      if( P_k(N,c)*P_k(N,a) < 0.0_r8) then
        b=c
      else
        a=c
      end if
    end do

    write(*,*) 'Error in : nibun()'
    stop

  end function nibun


  function P_k(k,x)
    integer,  intent(in) :: k
    real(r8), intent(in) :: x
    real(r8)             :: P_k

    real(r8) :: p0, p1,p2
    integer  :: n
    
    p0=1.0_r8
    p1=x
    if(k==0) then
      P_k = p0
      return
    else
      !do n = 1, k
      do n = 1, k-1
        p2 = ( (2.0_r8*dfloat(n)+1.0_r8)*x*p1 - dfloat(n)*p0) / ( dfloat(n)+1.0_r8)
        p0 = p1
        p1 = p2
      end do
    end if

    P_k = p1

  end function P_k

end module mod_gauss
