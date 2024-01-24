module mod_smode

  use mod_const, only : r8

  implicit none
  private
  real(r8), allocatable, save, private :: z(:,:)
  real(r8), allocatable, save, private :: sig(:)
  !----------------------------------------------------------------------------
  public :: shige_init
  public :: smtrx_smode
  public :: shige
  !----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------
  subroutine shige_init
    use mod_adm, only : &
        mk, maxl
    implicit none

    write(*,*) 'shige_init',mk,maxl
    allocate( z(-1:2*mk,-1:maxl-1) )
    allocate( sig(-1:maxl-1) )

  end subroutine shige_init
  !----------------------------------------------------------------------------
  subroutine smtrx_smode( eps ) 
    use mod_adm, only : &
        mk, maxl
    !use f95_lapack
    implicit none
    real(r8), intent(in) :: eps
    real(r8) :: s(2*mk+1,2*mk+1)
    real(r8) :: a(2*mk,2*mk)
    real(r8) :: a0(2*mk)
    real(r8) :: a1(2*mk-1)

    !real(r8) :: qq(2*mk+1,2*mk+1)
    !real(r8) :: work1(2*(2*mk+1)-2,2*(2*mk+1)-2)
    !real(r8) :: work2(2*(2*mk+1)-2)
    !real(r8) :: work3(3*(2*mk+1)-2)
    real(r8) :: alp
    real(r8) :: summ
    integer :: i, j
    integer :: n, nn, l
    !integer :: ierr
    integer :: info

    !real(r8) :: ab(2,2*mk+1)
    !real(r8) :: zz(2*mk+1,2*mk+1)
    !real(r8) :: rwork(22000)
    real(r8) :: work(22000)
    !real(r8) :: ss(2*mk+1,2*mk+1)
    real(r8) :: vr(2*mk+1,2*mk+1)
    real(r8) :: vl(2*mk+1,2*mk+1)
    real(r8) :: wr(2*mk+1)
    real(r8) :: wi(2*mk+1)
    !real(r8) :: aa(2*mk,2*mk)
    !real(r8) :: zzz(2*mk,2*mk)

    !real(r8) :: vl1(2*mk,2*mk)
    !real(r8) :: vr1(2*mk,2*mk)
    !real(r8) :: w1(2*mk)
    !real(r8) :: wi1(2*mk)

    real(r8) :: dn, en

    dn(n,alp)=(n-1)/alp/dsqrt(1.0_r8*(2*n-1)*(2*n+1))
    en(n,alp)=(n+2)/alp/dsqrt(1.0_r8*(2*n+1)*(2*n+3))

    z=0.0_r8
    s=0.0_r8

    alp = 1.0_r8 / dsqrt(eps)

    nn = 2 * mk + 1
    s(1,1) = -eps
    s(1,2) = 2.0_r8 * eps / dsqrt(3.0_r8)
    s(2,1) = -s(1,2)
    s(2,2) = 2.0_r8 + en(0,alp)**2 + dn(2,alp)**2

    do i = 3, nn
       n = 2*i-3
       s(i,i)   = en(n-1,alp)**2 + dn(n+1,alp)**2 + dble(n) * ( dble(n) + 1.0_r8 )
       s(i,i-1) = dn(n-1,alp)*en(n-1,alp)
       s(i-1,i) = s(i,i-1) 
    end do

    call dgeev('N','V',NN,S,NN,WR,WI,vl,NN,vr,NN,work,9900*2,info)

     do j = 1, nn
        summ = - ( vr(1,j) )**2
        do i = 2, nn
           summ = summ + ( vr(i,j) )**2
        end do
        summ = dsqrt( summ * wr(j) )
        vr(:,j) = vr(:,j) / summ
     end do

     do l = -1, maxl-1, 2
        j = (l+3) / 2
        sig(l) = -1.0_r8 / wr(j)
        do n = -1, 2*mk, 2
           i = (n+3) / 2
           z(n,l) = vr(i,j)
        end do
     end do

     a = 0.0_r8

     a0(1)  = 6.0_r8 + en(1,alp)**2 + dn(3,alp)**2
     a1(1)  = 0.0_r8

     do i = 2, 2*mk
        n = 2*i
        a0(i)   = en(n-1,alp)**2 + dn(n+1,alp)**2 + dble(n*(n+1))
        a1(i-1) = dn(n-1,alp) * en(n-1,alp)
     end do

     call dstev('V',NN-1,a0,a1,a,nn-1,work,info)

     do j = 1, 2*mk
        summ = 0.0_r8
        do i = 1, 2*mk
           summ = summ + a(i,j)**2
        end do
        summ = dsqrt( summ * a0(j) )
        a(:,j) = a(:,j) / summ
     end do
        
     do l = 2, maxl-1, 2
        j = l/ 2
        sig(l) = -1.0_r8 / a0(j)
        do n = 2, 2*mk, 2
           i = n / 2
           z(n,l) = a(i,j)
        end do
     end do

  end subroutine smtrx_smode
  !----------------------------------------------------------------------------
  !subroutine shige(eps, phi, h0)
  subroutine shige(eps, h0)
    use mod_adm, only : &
        mk, maxl
    use mod_gaussg
    use mod_const, only:  &
        pis2 => hfpai

    implicit none
    real(r8), intent(in) :: eps
!    real(r8), intent(in) :: phi(mk)
    real(r8), intent(inout) :: h0(2*mk,maxl,3)

    real(r8) ::  rn, dn, en
    real(r8) ::  alp

    real(r8) :: b(2*mk)
    real(r8) :: c(0:2*mk)
    real(r8) :: w0(maxl)
    !real(r8) :: sig(maxl)

    real(r8), allocatable :: cp(:)
    real(r8) :: phi2(2*mk)
    real(r8) :: p0(2*mk,-1:2*mk)
    real(r8) :: p1(2*mk,0:2*mk)
    real(r8) :: pb
    real(r8) :: y

    real(r8) :: cs
    real(r8) :: sn
    real(r8) :: sq
    real(r8) :: sinp
    real(r8) :: th
    

    integer :: n, l, ll
    integer :: i, iy
    integer :: isym

    real(r8) :: orth

    real(r8) :: yw(2*mk)
    real(r8) :: summ, sym
    !real(r8) :: oth, summ, sym
    real(r8) :: uu, vv, zz

    rn(n)=dsqrt(1.0_r8*n*(n+1))
    dn(n,alp)=(n-1)/alp/dsqrt(1.0_r8*(2*n-1)*(2*n+1))
    en(n,alp)=(n+2)/alp/dsqrt(1.0_r8*(2*n+1)*(2*n+3))
    alp=1.0_r8/dsqrt(eps)

    do i = 1,mk
       phi2(i) = -phi(i)
       yw(i)   = gusw(i)
       !write(*,*) i, phi2(i), yw(i)
    end do

    do i = mk+1, 2*mk
       phi2(i) =  phi(2*mk-i+1)
       yw(i)   = gusw(2*mk-i+1)
       !write(*,*) i, phi2(i), yw(i)
    end do

    do n = 0, 2*mk
       allocate(cp ( n/2 + 1) )
       cp = 0.0_r8
       orth= 0.0_r8
       call dalfk (n,0,cp)
          do i = 1,2*mk
             y=phi2(i)
             sinp = dsin(y)
             th   = pis2 - y
             call dlfpt (n,0,th,cp,pb)
             p0(i,n)=pb
             !p0(i,n)=pb*dsqrt( 0.5_r8*4.0_r8/(2.0_r8*dfloat(n)+1.0_r8) )
             orth=orth+p0(i,n)**2*yw(i)
          end do
         !write(*,*) 'orth', orth, dsqrt( 0.5_r8*4.0_r8/(2.0_r8*dfloat(n)+1.0_r8) )
       deallocate(cp)
    end do

    do n = 0,2*mk-1
       sq = dsqrt( (2.0_r8 * dfloat(n) + 1.0_r8 ) / (2.0_r8 * dfloat(n) + 3.0_r8 ) )
       do iy = 1, 2*mk
          !y=phi(i)
          !sinp = dsin(y)
          !th   = pis2 - y
          !sn = dsin(th)
          !cs = dcos(th)
          sn = dsin(phi2(iy))
          cs = dcos(phi2(iy))
          !cs = dsqrt( 1.0_r8 - sn**2 )
          p1(iy,n) = ( sn * p0(iy,n) - sq * p0(iy,n+1) ) * &
                       dfloat(n+1) / cs
          !if(iy.eq.mk.and.n.lt.3) !write(*,'(2i4,5f15.9)') n,iy,sq,sn,cs,p0(iy,n),p0(iy,n+1)
       end do
    end do

    !write(*,*) 'phi'
    !write(*,'(8f10.6)') phi 
    !write(*,*) 'p0'
    !write(*,'(8f10.6)') p0(:,1)
    !write(*,*) 'p1'
    !write(*,'(8f10.6)') p1(:,1)
    !write(*,*) 'sig in shige'
    !write(*,*) sig


    do l = -1, maxl-1
       ll = l + 1
       if( l == -1 ) ll = 1
       if( l /=  0 ) then
          w0(ll) = sig(l) 
          b(:)   = 0.0_r8
          c(:)   = 0.0_r8

          isym = mod(l+2,2)
          !write(*,*) 'isym=', isym, ll

          do n = 2-isym, 2*mk, 2
             b(n) = rn(n) * z(n,l)
          end do
          if(isym == 1) c(0) = z(-1,l) / alp - en(0,alp) * z(1,l)
          do n = 1+isym, 2*mk-1, 2
             c(n) = -dn(n,alp) * z(n-1,l) - en(n,alp) * z(n+1,l)
          end do

          do iy = 1, mk
             !sn = dsin(phi2(2*mk-iy+1))
             sn = dsin(phi2(iy))
             cs = dsqrt( 1.0_r8 - sn**2 )
             h0(iy,ll,1) = 0.0_r8
             h0(iy,ll,3) = 0.0_r8
             do n = 2-isym, 2*mk-1, 2
                h0(iy,ll,1) = h0(iy,ll,1) - b(n) / rn(n) * p1(iy,n)
                !h0(iy,ll,1) = h0(iy,ll,1) - b(n) / rn(n) * p1(2*mk-iy+1,n)
             end do
             do n = 1-isym, 2*mk, 2
                h0(iy,ll,3) = h0(iy,ll,3) - c(n) * p0(iy,n)
                !h0(iy,ll,3) = h0(iy,ll,3) - c(n) * p0(2*mk-iy+1,n)
             end do
             h0(iy,ll,2) = sig(l) * h0(iy,ll,1) / sn &
                         - alp * h0(iy,ll,3) / sn / cs 
          end do
       !end if

       if(isym == 1) sym =  1.0_r8
       if(isym == 0) sym = -1.0_r8
       do iy = 1, mk
          h0(2*mk-iy+1,ll,1) =  sym * h0(iy,ll,1)
          h0(2*mk-iy+1,ll,2) = -sym * h0(iy,ll,2)
          h0(2*mk-iy+1,ll,3) =  sym * h0(iy,ll,3)
       end do
       end if
    end do
    h0(:,:,2) = 0.0_r8
    do l = 1, maxl
      do ll = l,maxl
         summ = 0.0_r8
         do iy = 1, 2*mk
            uu = h0(iy,l,1)*h0(iy,ll,1)
            vv = h0(iy,l,2)*h0(iy,ll,2)
            zz = h0(iy,l,3)*h0(iy,ll,3)
            if(l.le.maxl) vv = 0.0_r8
            summ = summ + ( uu + vv + zz ) * yw(iy)
         end do
         !write(*,*) l, ll,  summ
      end do
   end do


  end subroutine shige
  !----------------------------------------------------------------------------
end module mod_smode
