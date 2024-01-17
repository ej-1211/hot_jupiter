module mod_special_function
  use mod_const, only : r8
  implicit none
  private

  real(r8), allocatable, public, save :: p0(:,:)
  real(r8), allocatable, public, save :: p1(:,:)

  public :: legendre

contains
  !----------------------------------------------------------------------------
  subroutine legendre(n1, n2, gp, gw)
    use mod_const
    implicit none
    integer,  intent(in) :: n1
    integer,  intent(in) :: n2
    real(r8), intent(in) :: gp(n1)
    real(r8), intent(in) :: gw(n1)

    real(r8) :: colat(n1)
    real(r8) :: del, pb, sq
    real(r8), allocatable :: cp(:)
    integer :: i, n

    allocate(p0(n1+1,0:n2))
    allocate(p1(n1+1,0:n2))

    do i = 1, n1
      colat(i) = pai*(90.0_r8 - dasin(gp(i))*180.0_r8/pai) / 180.0_r8
    end do

    del=0.5_r8
    do n = 0, n2
      allocate( cp(n/2+1) )
      call dalfk(n,0,cp)
      do i = 1, n1
        call dlfpt(n,0,colat(i),cp,pb)
        p0(i,n)=pb*dsqrt( del*4.0_r8/(2.0_r8*dfloat(n)+1.0_r8) )
      end do
      deallocate( cp )
    end do
    p0(n1+1,:)=1.0_r8
     
    do n = 0,n2
      sq=0.0_r8
      do i = 1,n1
        sq=sq+p0(i,n)**2*gw(i)
      end do
      sq=dsqrt(sq)
      p0(:,n)=p0(:,n)/sq
    end do

    do n = 1,n2-1
      sq=dsqrt((2.0_r8*dfloat(n)+1.0_r8)/(2.0_r8*dfloat(n)+3.0_r8))
      do i = 1,n1
        p1(i,n)=(sq*p0(i,n+1)-gp(i)*p0(i,n))*dfloat(n+1) / &
                (gp(i)**2 - 1.0_r8)
      end do
    end do

  end subroutine legendre
  !----------------------------------------------------------------------------
end module mod_special_function
