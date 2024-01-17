module mod_gaussg

  use mod_const, only : r8

  implicit none
  private
  real(r8), allocatable, public, save :: gusl(:)
  real(r8), allocatable, public, save :: gusw(:)
  !--- co-latitude in radian
  real(r8), allocatable, public, save :: an(:)
  !--- co-latitude in degree
  real(r8), allocatable, public, save :: ang(:)
  !--- latitude in radian
  real(r8), allocatable, public, save :: phi(:)
  !--- latitude in degree
  real(r8), allocatable, public, save ::  deg(:)

!  character(len=128), public, save :: ygrid_fname = 'grid.dat'
  !----------------------------------------------------------------------------
  public  :: gauss_init
  public  :: gaussg
  private :: ordleg
  !----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------
  subroutine gauss_init
    use mod_adm, only : mk
        !houghcalccnf_fname
    implicit none
!    namelist / merid_grid / &
!               ygrid_fname

!    open(unit=2, file=trim(houghcalccnf_fname))
!    read(2,merid_grid)
!    write(*,merid_grid)
!    close(2)

    allocate( gusl(mk) )
    allocate( gusw(mk) )
    allocate(   an(mk) )
    allocate(  ang(mk) )
    allocate(  phi(mk) )
    allocate(  deg(mk) )


  end subroutine gauss_init
  !----------------------------------------------------------------------------
  subroutine gaussg(nzero,sia,rad)
    use mod_const, only : &
        pai, dtr, hfpai

    implicit none
    integer,  intent(in)  :: nzero
    real(r8), intent(out) :: sia(nzero)
    real(r8), intent(out) :: rad(nzero)

    integer  :: i
    integer  :: ir, irp, irm
    real(r8) :: xlim, fi, fi1
    real(r8) :: piov2,fn
    real(r8) :: dn, dn1, a, b
    real(r8) :: g, gm, gp, gt
    real(r8) :: zz, ftemp, gtemp

    write(*,*) '#############################################'
    write(*,*) '#### Computation of Gaussian point start ####'
    write(*,*) '#############################################'

    xlim=1.0d-15
    ir=nzero+nzero
    fi=ir
    fi1=fi+1.0_r8
    piov2=0.5_r8*pai
    fn=piov2/dfloat(nzero)

    do i=1,nzero
      gusw(i)=dfloat(i) - 0.5_r8
    end do

    do i=1,nzero
      gusl(i)=dsin(gusw(i)*fn+piov2)
    end do

    dn=fi/dsqrt(4.0_r8*fi*fi-1.0_r8)
    dn1=fi1/dsqrt(4.0_r8*fi1*fi1-1.0_r8)
    a=dn1*fi
    b=dn*fi1
    irp=ir+1
    irm=ir-1

    do i = 1, nzero
   5  call ordleg(  g, gusl(i),  ir )
      call ordleg( gm, gusl(i), irm )
      call ordleg( gp, gusl(i), irp )
      gt = ( gusl(i) * gusl(i) - 1.0_r8 ) / ( a * gp - b * gm )
      zz = gusl(i)
      ftemp = zz - g * gt
      gtemp = zz - ftemp
      gusl(i) = ftemp
      if( dabs(gtemp) .gt. xlim ) go to 5
    end do

    do i=1,nzero
      a = 2.0_r8 * ( 1.0_r8 - gusl(i) * gusl(i) )
      call ordleg(b,gusl(i),irm)
      b=b*b*fi*fi
      gusw(i) = a* ( fi - 0.5_r8 ) / b
      rad(i) = dacos(gusl(i))
      sia(i) = dsin(rad(i))
    end do

    write(*,'("No.",10X,"gusl",10X,"an",10X,"ang",10X,"phi",10X,"deg",7X,"weight")') 
    do i = 1,nzero
      an(i)  = dacos( gusl(i) )
      ang(i) = an(i) / dtr
      phi(i) = hfpai - an(i)
      deg(i) = phi(i) / dtr
      write(*,'(i4,6f12.8)') i, gusl(i), an(i), ang(i), &
                             phi(i), deg(i), gusw(i)
    end do


  end subroutine gaussg
  !----------------------------------------------------------------------------
  subroutine ordleg( sx, coa, ir )
    implicit none
    real(r8), intent(inout) :: sx
    real(r8), intent(inout) :: coa
    integer,  intent(in)    :: ir

    real(r8) :: delta
    real(r8) :: sqr2
    real(r8) :: theta
    real(r8) :: c1
    real(r8) :: fn, fn2, fn2sq
    real(r8) :: ang, s1, c4, a, b
    real(r8) :: fk

    integer :: n, n1, k, kk
    integer :: irpp, irppm

    irpp=ir+1
    irppm=irpp-1
    delta=dacos(coa)
    sqr2=dsqrt(2.0_r8)
    theta=delta
    c1=sqr2

    do n=1,irppm
      fn=dfloat(n)
      fn2=fn+fn
      fn2sq=fn2*fn2
      c1=c1*dsqrt(1.0_r8-1.0_r8/fn2sq)
    end do

    n   =  irppm
    ang =  fn * theta
    s1  =  0.0_r8
    c4  =  1.0_r8
    a   = -1.0_r8
    b   =  0.0_r8
    n1  =  n + 1

    do kk=1,n1,2
      k = kk-1
      if(k.eq.n) c4 = 0.5_r8 * c4
      s1  = s1 + c4 * dcos(ang)
      a   = a + 2.0_r8
      b   = b + 1.0_r8
      fk  = dfloat(k)
      ang = theta * ( fn - fk - 2.0_r8)
      c4  = ( a * (fn - b + 1.0_r8 ) / ( b* ( fn2 - a ))) * c4
    end do

    sx = s1 * c1

  end subroutine ordleg
  !----------------------------------------------------------------------------
end module mod_gaussg
