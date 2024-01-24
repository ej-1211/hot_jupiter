module mod_interpolation

  use mod_const, only : r8

  implicit none

  public  :: linear_interpolation
  private :: find_closest_grid
  !----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------
  subroutine linear_interpolation(n, m, x1, y1, x2, y2)
    implicit none
    integer,  intent(in)  :: n     ! Num. of points for original data
    integer,  intent(in)  :: m     ! Num. of points which you want to interpolate
    real(r8), intent(in)  :: x1(n) ! value of grid for original data
    real(r8), intent(in)  :: y1(n) ! value of original data
    real(r8), intent(in)  :: x2(m) ! value of grid which you want to interpolate
    real(r8), intent(out) :: y2(m) ! value of interpolated data

    integer  :: iz, ibelow, iabove
    real(r8) :: pbelow, pbbelow, pabove, deltap

    do iz = 1,m
      call find_closest_grid(x2(iz), n, x1(:), ibelow, iabove)
      pbelow = dlog(x1(ibelow))
      if (iabove == 0) then
        pabove = dlog(0.10_r8)
      else
        pabove = dlog(x1(iabove))
      end if
      deltap = pabove - pbelow 

      if ( iabove /= 0 ) then
       if(ibelow==iabove) then
         y2(iz)=y1(ibelow)
       elseif(ibelow > iabove) then
         y2(iz)=( y1(ibelow)*(pabove-dlog(x2(iz))) +  &
                  y1(iabove)*(dlog(x2(iz))-pbelow) ) / deltap 
       end if
      elseif ( iabove == 0 ) then
       pbbelow = dlog(x1(ibelow+1))      
       y2(iz) = y1(ibelow)*(dlog(x2(iz))-pbbelow)*(dlog(x2(iz))-pabove)/ & 
                ((pbelow-pbbelow)*(pbelow-pabove)) &
              + y1(ibelow+1)*(dlog(x2(iz))-pbelow)*(dlog(x2(iz))-pabove)/ & 
                ((pbbelow-pabove)*(pbbelow-pbelow))
!       write(*,*) 'linear_interpolation: extrapolation above top level of input data'
      else
       write(*,*) 'linear_interpolation: level ',ibelow,' not under/same as level ',iabove
      end if 
    end do

  end subroutine linear_interpolation 
  !----------------------------------------------------------------------------
  subroutine find_closest_grid(x1,m,x2,ibelow,iabove)
    implicit none
    real(r8), intent(in)  :: x1
    integer,  intent(in)  :: m
    real(r8), intent(in)  :: x2(m)
    integer,  intent(out) :: ibelow, iabove

    integer :: ibe, iab, k
    real(r8) :: pabove, pbelow, pz, zlittle

    pz = 0.0_r8
    ibelow = m
    iabove = 0
    zlittle = 1.0d-15

! Find a nearest grid point below a given point

    pbelow = x2(m) - x1 
    ibe = m
    do k = m-1, 1, -1
      pz = x2(k) - x1
      if(pz>0.0_r8 .and. pz<pbelow) then
        pbelow = pz
        ibe    = k
      end if
    end do
    ibelow = ibe
    
! Find a nearest grid point above a given point
    pabove = x1 - x2(1) 
    iab    = 1
    if(pabove < 0.0_r8) then
!     new level is above the top level of original grid  
      pabove = 0.0_r8
      iab    = 0
    end if    
    do k = 2, m
      pz = x1 - x2(k) 
      if(pz>0.0_r8 .and. pz<pabove) then
!       new sigma level is above the current level 
        pabove = pz
        iab    = k
      end if
    end do
    iabove = iab

  end subroutine find_closest_grid
  !----------------------------------------------------------------------------
end module mod_interpolation
