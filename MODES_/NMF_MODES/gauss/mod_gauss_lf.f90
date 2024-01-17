module mod_gauss_lf
  implicit none
  contains
  !----------------------------------------------------------------------------
  subroutine gauss(N, x, w)
    implicit none
    integer, intent(in) :: N
    real(16), intent(inout) :: x(N)
    real(16), intent(inout) :: w(N)
    real(16) :: a, b
    real(16) :: lambak
    real(16) :: pi
    integer  :: i, m, k, kk

!    call omp_set_num_threads(12)

    pi = qacos(-1.0q0)
    
    if(mod(N,2)==0) then
      m = N/2
    else
      m = (N-1)/2
    end if 

!$omp parallel
!$omp do
    do i = 1,m
      write(*,*) i
      a = qsin(pi*qfloat(N-2*i)/qfloat(2*N)) 
      b = qsin(pi*qfloat(N+2-2*i)/qfloat(2*N)) 
      x(N+1-i) = NIBUN(N,a,b)
      x(i) = -x(N+1-i)
    end do
!$omp end do
!$omp end parallel

    if(mod(N,2)==1) x(m) = 0.0q0

!$omp parallel
!$omp do
    do i = 1,N
      w(i) = 0.0q0
      w(i) = 0.5q0
      do k = 1,N
        lambak = 2.0q0 / ( 2.0q0*qfloat(k)+1.0q0)
        w(i) = w(i) + P_k(k,x(i))**2 / lambak
      end do
      w(i) = 1.0q0 / w(i)
    end do
!$omp end do
!$omp end parallel

  end subroutine gauss
  !----------------------------------------------------------------------------
  function nibun(N,a,b)
    integer :: N
    integer :: k
    real(16) :: a, b, c
    real(16) :: delta, epsi
    real(16) :: nibun

    delta = 1.0q-15
    epsi  = 1.0q-15

    do k = 1, 100000
      c = (a+b)/2.0q0
      if( qabs(P_k(N,c)) < delta .or. qabs(a-b) < epsi ) then
        nibun = c
        return
      end if 
      if( P_k(N,c)*P_k(N,a) < 0.0q0) then
        b=c
      else
        a=c
      end if
    end do

    write(*,*) 'Error'
    stop

  end function nibun
  !----------------------------------------------------------------------------
  function P_k(k,x)
    integer :: n, k
    real(16) :: p0, p1,p2, x
    real(16) :: p_k
    
    p0=1.0q0
    p1=x
    do n = 1, k-1
      p2 = ( (2.0q0*qfloat(n)+1.0q0)*x*p1 - qfloat(n)*p0) / ( qfloat(n)+1.0q0)
      p0 = p1
      p1 = p2
    end do

    P_k = p1

  end function P_k
  !----------------------------------------------------------------------------
end module mod_gauss_lf
