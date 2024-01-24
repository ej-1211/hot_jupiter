module mod_numvsf_sigma
  implicit none
  !----------------------------------------------------------------------------
  public :: numvsf_sigma
  !----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------
  subroutine numvsf_sigma
    use mod_const
    use mod_adm
    implicit none
    integer :: k, kk, i
    real(r8), allocatable :: dsig(:)
    real(r8), allocatable :: c(:)
    real(r8), allocatable :: d(:)
    real(r8), allocatable :: rcon(:)
    real(r8), allocatable :: qm(:,:)
    real(r8), allocatable :: ort(:)
    real(r8) :: sufst
    real(r8) :: rsufc
    real(r8) :: acon

    real(r8), allocatable :: w(:)
    real(r8), allocatable :: wr(:)
    real(r8), allocatable :: wi(:)
    real(r8), allocatable :: vr(:,:)
    real(r8), allocatable :: vi(:,:)
    real(r8), allocatable :: work(:)
    real(r8), allocatable :: oth(:)
    real(r8) :: summ
    integer :: lwork
    integer :: info

    allocate(dsig(mp+1))
    allocate(c(mp))
    allocate(d(mp))
    allocate(rcon(mp))
    allocate(qm(mp,mp))
    allocate(ort(mp))

    lwork=4*mp
    allocate(w(mp))
    allocate(wr(mp))
    allocate(wi(mp))
    allocate(vr(mp,mp))
    allocate(vi(mp,mp))
    allocate(work(lwork))
    allocate(oth(mp))

    sufst=stab(1)
    acon =sufst/suft

    do k = 2,mp
      dsig(k)=vgrid(k-1)-vgrid(k)
    end do
    dsig(   1)=2.0_r8*(1.0_r8-vgrid(1))
    dsig(mp+1)=2.0_r8*vgrid(mp)

    do k = 1,mp
      d(k)=dsqrt(0.5_r8*(dsig(k)+dsig(k+1)))
    end do

    do k = 1,mp-1
      rcon(k)=g*hstd/(rd*stab(k+1))
    end do
    rsufc=g*hstd/(rd*sufst)
    
    do k = 2,mp
      c(k)=-0.5_r8*(vgrid(k)+vgrid(k-1))*rcon(k-1)/dsig(k)
    end do
    c(1)=-rsufc*acon/(1.0_r8+0.5_r8*acon*dsig(1))

    do k = 1,mp+1
      write(*,*) k, dsig(k)
    end do 
    do k = 1,mp
      write(*,'(i4,3f20.9)') k, c(k), d(k), rcon(k)
    end do 

    qm(:,:)=0.0_r8

    do k = 2,mp-1
      qm(k,k-1)=c(k)/(d(k-1)*d(k))
      qm(k,k  )=-(c(k)+c(K+1))/(d(k)*d(k))
      qm(k,k+1)=c(k+1)/(d(k)*d(k+1))
    end do

!    qm(1,1)=-(c(1)+c(2))/(d(1)*d(1)) 
!    qm(1,2)= c(2)/(d(1)*d(2)) 
    qm(1,1)=-(c(2))/(d(1)*d(1)) 
    qm(1,2)= c(2)/(d(1)*d(2)) 
    qm(mp,mp-1)=c(mp)/(d(mp-1)*d(mp))
    qm(mp,mp  )=-c(mp)/(d(mp)*d(mp))

    call dsyev('V','U',mp,qm,mp,w,work,lwork,info)
    write(*,*) hstd/w

    do k = 1,num_vmode
      evht(k)=hstd/w(k)
    end do
    
    do k = 1,mp
      summ=0.0_r8
      do kk = 1,mp
        summ = summ + qm(kk,k)**2
      end do
      summ = dsqrt(summ)
      qm(:,k) = qm(:,k)/summ
    end do

    do k = 1,mp
      oth(:)=0.0_r8
      do kk = 1,mp
        do i = 1,mp
          oth(kk)=oth(kk)+qm(i,k)*qm(i,kk)
        end do
      end do
      if(ocheck) then
        write(*,*) oth
      end if
    end do

    vsf(1:mp,1:num_vmode)=qm(1:mp,1:num_vmode)

  end subroutine numvsf_sigma
  !----------------------------------------------------------------------------
end module mod_numvsf_sigma
