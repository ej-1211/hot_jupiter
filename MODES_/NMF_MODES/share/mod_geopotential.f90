module mod_geopotential

  implicit none
  
contains
  !----------------------------------------------------------------------------
  subroutine compute_geoheight
    use mod_adm, only : &
        nx, ny, nz,     &
        aa, bb,         &
        t_input,        &
        q_input,        &
        ppp,            &
        orog,           &
        z_input 
    use mod_const

! computes geopotential on hybrid model levels from surface pressure,  t_input. and spec. humidity 
! includes the use of the surface geopotential 

 implicit none

 real(r8), dimension(0:nz) :: phalf
 real(r8), dimension(nx,ny,0:nz) :: hhalf
 real(r8) :: pfull, pdel, dP, tv
 integer :: ix, iy, iz

 z_input(:,:,:) = 0.0_r8
 hhalf(:,:,:) = 0.0_r8
 phalf(:) = 0.0_r8

! Pressure at half levels and height fields at half and full model levels
!
  do ix = 1, nx
  do iy = 1, ny

! surface height (=orography) taken from the input surface geopotential field
   hhalf(ix,iy,nz) = orog(ix,iy)/g

   do iz = nz,1,-1

     phalf(iz) = aa(iz) + bb(iz) * dexp(ppp(ix,iy))
     phalf(iz-1) = aa(iz-1) + bb(iz-1) * dexp(ppp(ix,iy))
     pfull = (phalf(iz)+phalf(iz-1))*0.5_r8 ! mid-level pressure
     pdel = phalf(iz)-phalf(iz-1) ! layer thickness
     dP = 0.5_r8*pdel/pfull

! compute virtual  t_inputerature
     tv = (1.0_r8+0.61_r8* q_input(ix,iy,iz))* t_input(ix,iy,iz)

! mid-level (height) and half level (hhalf) geopotential 
     z_input(ix,iy,iz) = hhalf(ix,iy,iz) + tv*dP*Rd/g
     hhalf(ix,iy,iz-1) = hhalf(ix,iy,iz) + tv*2*dP*Rd/g

   end do
  end do
  end do
  end subroutine compute_geoheight
  !----------------------------------------------------------------------------
end module mod_geopotential
