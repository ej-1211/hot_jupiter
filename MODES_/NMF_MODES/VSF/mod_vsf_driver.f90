module mod_vsf_driver

  use mod_const,        only : r8
  use mod_write_netcdf, only : write_netcdf_vsf

  implicit none
  public :: vsf_driver
  public :: vertical_expansion
  public :: vertical_inverse
  !----------------------------------------------------------------------------
contains
  !----------------------------------------------------------------------------
  subroutine vsf_driver
    use mod_adm
    use mod_numvsf_sigma
    use mod_special_function

    implicit none

    call numvsf_sigma

    if(trim(equiheight_fname)=='') then

    else
      call output_1d(equiheight_fname, evht(1:num_vmode), num_vmode, &
                     1, .true., 8)
    end if 

    if(trim(vsf_fname)=='') then

    else
      call output_2d(vsf_fname, vsf, mp, num_vmode, 1, .true., 8)
      call write_netcdf_vsf(vsf_fname, mp, num_vmode, vsf)
    end if

  end subroutine vsf_driver
  !----------------------------------------------------------------------------
  subroutine vertical_expansion
    use mod_adm
    implicit none
    integer :: ix, iy, ip, m
    real(r8) :: Vst(mp), Eq
    
    vu(:,:,:) = 0.0_r8
    vv(:,:,:) = 0.0_r8
    vz(:,:,:) = 0.0_r8
    
    Vst(:) = 0.0_r8
    Eq = 0.0_r8

    do m = 1,num_vmode
      Vst(1:mp) =  vtmp(1:mp,m)
      Eq = sqgh(m)
      do ip = 1,mp 
        do ix = 1,nx
          do iy = 1,ny
            vu(ix,iy,m)=vu(ix,iy,m)+uuu(ip,ix,iy)*Vst(ip)/Eq
            vv(ix,iy,m)=vv(ix,iy,m)+vvv(ip,ix,iy)*Vst(ip)/Eq
            vz(ix,iy,m)=vz(ix,iy,m)+zzz(ip,ix,iy)*Vst(ip)/evht(m)
           ! vu(ix,iy,m)=vu(ix,iy,m)+uuu(ip,ix,iy)*Vst(ip)
           ! vv(ix,iy,m)=vv(ix,iy,m)+vvv(ip,ix,iy)*Vst(ip)
           ! vz(ix,iy,m)=vz(ix,iy,m)+zzz(ip,ix,iy)*Vst(ip)
            !print *, (vz(ix, iy, m))
          end do
          !print *, vv(ix,:,m)
        end do
      end do
    end do

  end subroutine vertical_expansion
  
  !----------------------------------------------------------------------------
  subroutine vertical_inverse
    use mod_adm
    implicit none
    integer :: ix, iy, ip, m
    real(r8) :: Eq
    
    u_input(:,:,:) = 0.0_r8
    v_input(:,:,:) = 0.0_r8
    z_input(:,:,:) = 0.0_r8

    Eq = 0.0_r8
    
    do ip = 1, mp 
      do m = 1, num_vmode
        Eq = sqgh(m)
        do ix = 1,nx 
          do iy = 1,ny
            !u_input(ix,iy,ip)=u_input(ix,iy,ip)+vu(ix,iy,m)*vsf(mp-ip+1,m)
            !v_input(ix,iy,ip)=v_input(ix,iy,ip)+vv(ix,iy,m)*vsf(mp-ip+1,m)
            !z_input(ix,iy,ip)=z_input(ix,iy,ip)+vz(ix,iy,m)*vsf(mp-ip+1,m)
            u_input(ix,iy,ip)=u_input(ix,iy,ip)+vu(ix,iy,m)*vsf(mp-ip+1,m)*Eq
            v_input(ix,iy,ip)=v_input(ix,iy,ip)+vv(ix,iy,m)*vsf(mp-ip+1,m)*Eq
            z_input(ix,iy,ip)=z_input(ix,iy,ip)+vz(ix,iy,m)*vsf(mp-ip+1,m)*evht(m)
          end do
        end do
      end do
    end do


!    open(100,file='output.data',form='unformatted',access='direct', &
!         recl=4*nx*my*mp)
!    !write(100,rec=1) (((sngl(uuu(ip,ix,iy)),ix=1,nx),iy=1,my),ip=1,mp)
!    !write(100,rec=2) (((sngl(vvv(ip,ix,iy)),ix=1,nx),iy=1,my),ip=1,mp)
!    !write(100,rec=3) (((sngl(zzz(ip,ix,iy)),ix=1,nx),iy=1,my),ip=1,mp)
!    write(100,rec=1) sngl(u_input)
!    write(100,rec=2) sngl(v_input)
!    write(100,rec=3) sngl(z_input)
!    close(100)

  end subroutine vertical_inverse
  !----------------------------------------------------------------------------
end module mod_vsf_driver
