module mod_geopotential_ecmwf

  implicit none
  
contains
  !----------------------------------------------------------------------------
  subroutine compute_geoheight_ecmwf
    use mod_adm, only : &
        nx, ny, nz,     &
        aa, bb,         &
        t_input,        &
        q_input,        &
        ppp,            &
        orog,           &
        z_input,        &
        output_3d_gmt, output_2d_gmt
    use mod_const
      
    implicit none
    real(r8), dimension(0:nz) :: phalf
    real(r8), dimension(nx,ny,0:nz) :: hhalf
    real(r8) ::  alpha, pfull, dP, dlogP, &
                        R, Cp, tpres
    integer :: ix, iy, iz

    z_input(:,:,:) = 0.0_r8
    hhalf(:,:,:) = 0.0_r8
    phalf(:) = 0.0_r8
 
!    call output_3d_gmt('InputT_','(512f12.2,1x)','200001010000000',t_input, nx, ny, nz)
!    call output_3d_gmt('InputQ_','(512f12.8,1x)','200001010000000',q_input, nx, ny, nz)
!     call output_2d_gmt('orog.dat','(512f12.4,1x)',orog, nx, ny)
!     call output_2d_gmt('p_surface.dat','(512f12.4,1x)',ppp, nx, ny)  

!
! Pressure at half levels and height fields at half and full model
! levels
!
    do ix = 1, nx
      do iy = 1, ny
        hhalf(ix,iy,nz) = orog(ix,iy)/g
        do iz = nz,1,-1
          phalf(iz) = aa(iz) + bb(iz) * dexp(ppp(ix,iy))
          phalf(iz-1) = aa(iz-1) + bb(iz-1) * dexp(ppp(ix,iy))
          pfull = (phalf(iz)+phalf(iz-1))*0.5_r8
          dP = phalf(iz)-phalf(iz-1)
          if(iz == 1 ) then
            tpres = 0.1_r8
!            hydr = 1.0_r8
            dlogP = dlog(phalf(iz)/tpres)
            alpha = dlog(2.0_r8) !hydr  
          else
            dlogP = dlog(phalf(iz)/phalf(iz-1))
            alpha = 1.0_r8 - phalf(iz-1)/dP*dlogP
          end if
          Cp = cpd*(1.0_r8-q_input(ix,iy,iz)) + Cvd*q_input(ix,iy,iz)
          R  =  Rd*(1.0_r8-q_input(ix,iy,iz)) +  Rv*q_input(ix,iy,iz)
          z_input(ix,iy,iz) = hhalf(ix,iy,iz) + R*t_input(ix,iy,iz)*alpha/g
          hhalf(ix,iy,iz-1) = hhalf(ix,iy,iz) + R*t_input(ix,iy,iz)*dlogP/g
         end do

      end do
    end do
 
  end subroutine compute_geoheight_ecmwf
  !----------------------------------------------------------------------------
end module mod_geopotential_ecmwf
