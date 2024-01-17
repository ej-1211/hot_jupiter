module mod_hough_driver
  implicit none

  public :: hough_expansion
  public :: hough_inverse
  !----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------
  subroutine hough_expansion
    use mod_adm
    use mod_const
    use mod_interpolation
    integer :: n, m, l, iy
    real(r8) :: enrg(0:num_zw)
    real(r8) :: venrg(num_vmode)
    real(r8) :: cc
    double complex   :: wave2(my,3)

    w0(:,:,:) = 0.0_r8
    ww0(:,:,:) = 0.0_r8
    do m = 1, num_vmode 
      do n = 0, num_zw
        wave2(:,:)=wave(n,:,m,:)
        do l = 1, 3*maxl
          do iy = 1,my
            w0(l,m,n) = w0(l,m,n) + &
            ( wave2(iy,1) * hough(my-iy+1,l,1,m,n)       &
            + wave2(iy,2) * hough(my-iy+1,l,2,m,n) * img &
            + wave2(iy,3) * hough(my-iy+1,l,3,m,n) ) * ygrid_weight(iy)
            if ( savehoughwind ) then
              ww0(l,m,n) = ww0(l,m,n) + &
              ( wave2(iy,1) * hough(my-iy+1,l,1,m,n)       &
              + wave2(iy,2) * hough(my-iy+1,l,2,m,n) * img ) * ygrid_weight(iy)
            end if 
          end do
        end do
      end do
    end do

    enrg(:) = 0.0_r8
    venrg(:) = 0.0_r8

    do n = 0, num_zw
      cc=1.0_r8
      if(n==0) cc=0.5_r8
      do m = 1, num_vmode
        do l = 1, maxl*3
          enrg(n) = enrg(n) + cc * g* evht(m) * &
                     dreal( w0(l,m,n)*dconjg(w0(l,m,n)) )
        end do
      end do
    end do

    do n = 1, num_zw
      cc=1.0_r8
      if(n==0) cc=0.5_r8
      do m = 1, num_vmode
        do l = 1, maxl*3
          venrg(m) = venrg(m) + cc * g* evht(m) * &
                     dreal( w0(l,m,n)*dconjg(w0(l,m,n)) )
        end do
      end do
    end do
    
    write(*,*) 'Energy per zonal wavenumber: '
    do n = 0, num_zw
      write(*,*) n, enrg(n)
    end do

    write(*,*) 'Energy in vertical modes w/t k=0: '
    do m = 1, num_vmode
      write(*,*) m, venrg(m)
    end do

  end subroutine hough_expansion
  !----------------------------------------------------------------------------
  subroutine hough_inverse
    use mod_adm
    use mod_const
    integer :: n, m, l, iy

    wave(:,:,:,:)=0.0_r8

    do n = 0, num_zw
      do m = 1, num_vmode
        do l = shough, ehough
          do iy = 1, my
            wave(n,iy,m,1)=wave(n,iy,m,1)+w0(l,m,n)*hough(my-iy+1,l,1,m,n)
            wave(n,iy,m,2)=wave(n,iy,m,2)+w0(l,m,n)*hough(my-iy+1,l,2,m,n)*(-img)
            wave(n,iy,m,3)=wave(n,iy,m,3)+w0(l,m,n)*hough(my-iy+1,l,3,m,n)
          end do
        end do
      end do
    end do
  end subroutine hough_inverse
  !----------------------------------------------------------------------------
end module mod_hough_driver
