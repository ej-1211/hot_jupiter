module mod_fft

  use mod_const, only : r8

  implicit none

  public :: fft_forward
  public :: fft_backward
  private :: fft_1d_forward
  private :: fft_1d_backward
  !----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------
  subroutine fft_forward
    use mod_adm
    implicit none
    integer :: iy, m

    do m = 1, num_vmode
      do iy = 1,ny
        call fft_1d_forward(nx, num_zw, vu(:,iy,m), wave(:,iy,m,1))
        call fft_1d_forward(nx, num_zw, vv(:,iy,m), wave(:,iy,m,2))
        call fft_1d_forward(nx, num_zw, vz(:,iy,m), wave(:,iy,m,3))
      end do
    end do

  end subroutine fft_forward
  !----------------------------------------------------------------------------
  subroutine fft_backward
    use mod_adm
    implicit none
    integer :: iy, m

    do m = 1, num_vmode
      do iy = 1,ny
        call fft_1d_backward(nx, num_zw, vu(:,iy,m), wave(:,iy,m,1))
        call fft_1d_backward(nx, num_zw, vv(:,iy,m), wave(:,iy,m,2))
        call fft_1d_backward(nx, num_zw, vz(:,iy,m), wave(:,iy,m,3))
      end do
    end do

  end subroutine fft_backward
  !----------------------------------------------------------------------------
  subroutine fft_1d_forward(nx, num_zw, ddd, www)
    implicit none
    integer,        intent(in)    :: nx
    integer,        intent(in)    :: num_zw
    real(r8),       intent(inout) :: ddd(nx)
    double complex, intent(inout) :: www(0:num_zw)

    double complex :: ddd1(nx)
    integer :: i
    integer :: lensav
    integer :: inc
    integer :: ier
    real(r8), allocatable :: wsave(:)
    real(r8) :: work(2*nx)

    inc = 1
    lensav=(nx+int(log(real(nx)))+4)*2
    allocate(wsave(lensav))

    ddd1(:)=dcmplx(ddd(:),0.0_r8)
    call cfft1i(nx, wsave, lensav, ier)        
    call cfft1f(nx, inc, ddd1, nx, wsave, lensav, work, 2*nx, ier)

    do i = 0,num_zw
      www(i)=ddd1(i+1)
    end do

    deallocate(wsave)

  end subroutine fft_1d_forward
  !----------------------------------------------------------------------------
  subroutine fft_1d_backward(nx, num_zw, ddd, www)
    implicit none
    integer,        intent(in)    :: nx
    integer,        intent(in)    :: num_zw
    real(r8),       intent(inout) :: ddd(nx)
    double complex, intent(inout) :: www(0:num_zw)

    double complex :: www1(nx)
    integer :: i
    integer :: lensav
    integer :: inc
    integer :: ier
    real(r8), allocatable :: wsave(:)
    real(r8) :: work(2*nx)

    inc = 1
    lensav=(nx+int(log(real(nx)))+4)*2
    allocate(wsave(lensav))
    www1=0.0_r8
    do i = 0,num_zw 
      www1(i+1)=www(i)
    end do
    do i = nx,nx-num_zw+1,-1
      www1(i)=dconjg(www1(nx-i+2))
    end do
    call cfft1i(nx, wsave, lensav, ier)        
    call cfft1b(nx, inc, www1, nx, wsave, lensav, work, 2*nx, ier)
    do i = 1,nx
      ddd(i)=dreal(www1(i))
    end do

    deallocate(wsave)

  end subroutine fft_1d_backward
  !----------------------------------------------------------------------------
end module mod_fft
