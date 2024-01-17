module mod_normal

  use mod_const, only : r8
  implicit none
  private
  public :: read_inputdata
  public :: NMF_output
  public :: NMF_outputw
  private :: global_mean
  public :: grid_data_order
  !
contains
  !----------------------------------------------------------------------------
  subroutine read_inputdata
    use mod_adm
    use mod_time
    use mod_const, only : g, rd
!!!    use mod_read_grib 
    use mod_read_netcdf
    use mod_interpolation
    use mod_geopotential_ecmwf
    use mod_geopotential
    implicit none
    integer :: iz, ix, iy, ip
    integer :: i 
    real(r8) :: glob(nz)
    real(r8) :: phalf(0:nz), pfull(nz)
    real(r8) :: vgrid_lnsigma(mp)
    character(256) :: tmp_fname
    real(r8) :: zzz_(mp, nx, ny)  
    
! IFs for the specified datatype
    if (trim(dataformat_input) == 'grib' ) then  
      write(*,*) "============= Input data file(s) =============="
      do i = 1, numoffile
        ifile_grib(i)=trim(ifile_head(i))//trim(cdate) 
!        write(*,*) trim(ifile_head(i)), trim(cdate)
        write(*,*) trim(ifile_grib(i))
      end do
      if (len(trim(ifile_orog)) /= 0)  then
         ifile_grib(numoffile+1) = trim(ifile_orog)
         write(*,*) "Orography file: ", ifile_grib(numoffile+1)
      end if
      write(*,*) "==============================================="       
      write(*,*) ""

      ! Reading GRIB file
      ! Output data are sorted NP to SP and bottom to top level upward
!!!      call read_grib_hybrid
    
      ! change  from the bottom to top to order from top to bottom for the input fields
!!!      call grid_data_order

!    call output_3d_3var_gmt('InputData_',aformat,cdate,u_input,v_input,t_input,nx,ny,nz)

    elseif (trim(dataformat_input) == 'binary' ) then
      write(*,*) '============ Input binary files ==============="'
!      call input_3d(trim(ifile_U), u_input, nx, ny, nz, 1, .true., 4)
!      call input_3d(trim(ifile_V), v_input, nx, ny, nz, 1, .true., 4)
!      call input_3d(trim(ifile_Z), z_input, nx, ny, nz, 1, .true., 4)
!      call input_3d(trim(ifile_T), t_input, nx, ny, nz, 1, .true., 4)
!      call input_2d(trim(ifile_P), ppp, nx, ny, 1, .true., 4)
          
    elseif (trim(dataformat_input) == 'netcdf' ) then
      write(*,*) ""
      write(*,*) "============ Input NetCDF files ==============="
      do i = 1, numoffile
        ifile_netcdf(i)=trim(ifile_head(i))//trim(cdate)//trim(".nc")
!        write(*,*) trim(ifile_head(i)), trim(cdate), trim(".nc") 
        write(*,*) trim(ifile_netcdf(i))
      end do
      if (len(trim(ifile_orog)) /= 0)  then
        ifile_netcdf(numoffile+1) = trim(ifile_orog)
        write(*,*) "Orography file: ", ifile_netcdf(numoffile+1)
      end if
      write(*,*) "==============================================="
      write(*,*) ""

    ! Reading NetCDF file
    ! Output data are sorted NP to SP and bottom to top level upward
      call read_netcdf_hybrid
    ! change from the bottom to top to order from top to bottom for the input fields 
!##########
! uncomment if JRA in netcdf is reading 
! call grid_data_order        !<< this part is already done in read_netcdf_hybrid!! >>
    else
      write(*,*) "Please specify a proper dataformat_input parameter in the namelist. Choose: 'grib' or 'netcdf'"
    end if 
        
! Compute geopotential height
    if (orig == 'ECMWF') then
      call compute_geoheight_ecmwf
    else
      call compute_geoheight
    end if
    
!   save input files on model levels
    if ( saveasci ) then
      call output_3d_3var_gmt(afname,aformat,cdate,u_input,v_input,z_input,nx,ny,nz)
    end if
          
!   if(z-Interpolate) then
    ! Interpolate from model levels to other coordinate
    if (zgrid_type == 'hybrid' ) then
     do iy = 1,ny
       do ix = 1,nx
          phalf(0) = aa(0) + bb(0) * dexp(ppp(ix,iy)) 
          do iz = 1,nz
           phalf(iz) = aa(iz) + bb(iz) * dexp(ppp(ix,iy))
           pfull(iz) = (phalf(iz-1)+phalf(iz))*0.5_r8 ! hybrid full levels pressure
           orig_lnsigma(iz)=pfull(iz) !sorted top to bottom
          end do
          do ip = 1,mp
            vgrid_lnsigma(ip)=vgrid(mp-ip+1)*dexp(ppp(ix,iy)) !now sorted top to bottom
          end do    
          call linear_interpolation(nz,mp,  &
                                    orig_lnsigma(1:nz), u_input(ix,iy,1:nz), &
                                    vgrid_lnsigma(1:mp), uuu(1:mp,ix,iy))
          call linear_interpolation(nz,mp,  &
                                    orig_lnsigma(1:nz), v_input(ix,iy,1:nz), &
                                    vgrid_lnsigma(1:mp), vvv(1:mp,ix,iy))
          call linear_interpolation(nz,mp,  &
                                    orig_lnsigma(1:nz), t_input(ix,iy,1:nz), &
                                    vgrid_lnsigma(1:mp), ttt(1:mp,ix,iy))
          call linear_interpolation(nz,mp,  &
                                    orig_lnsigma(1:nz), z_input(ix,iy,1:nz), &
                                    vgrid_lnsigma(1:mp), zzz(1:mp,ix,iy))
       end do
     end do
    
    else
     print *,'Other vertical interpolation than hybrid to sigma levels are not available'
    end if
!   end if (end of z_interpolate)


! compute the globally averaged mean level temperature
!
    call global_mean(nx,ny,mp,dcos(dasin(ygrid)),ttt,glob)

    write(*,*) 'Check after interpolation for Z'
    write(*,*) 'minimum, maximum and average T at every sigma model level'
    do ip = 1,mp
      write(*,'(i4,3F10.4)') ip, minval(ttt(ip,:,:)), maxval(ttt(ip,:,:)), &
                             glob(ip)
    end do
    write(*,*) 'maximum and minimum for surface pressure'
    write(*,*) minval(dexp(ppp)/100.0d0), maxval(dexp(ppp)/100.0d0)

! save the profile of the global mean T for later use during the inversion
      if(savemeant) then
          tmp_fname=trim(meant_fname)//trim(cdate)
          call output_1d(tmp_fname, glob, mp, 1, .true., 8)
      end if
        
!   Kasahara&Puri modified geopotential height
! 
      do ix = 1,nx
       do iy = 1,ny
        do ip = 1,mp 
         zzz(ip,ix,iy) = zzz(ip,ix,iy)+Rd*glob(ip)*ppp(ix,iy)/g
        end do
       end do 
      end do
    call read_3d_array(trim(source_fname), zzz_, mp, nx, ny)    
!  save output files on sigma levels
!    open(unit=101,file='Inputdata_sigma.dat')
!    do iz = 1, mp
!     do iy = 1, ny
!        write(101,'(96E14.4,1x)') uuu(iz,:,iy)
!    end do
!    end do
!    do iz = 1, mp
!     do iy = 1, ny
!        write(101,'(96E14.4,1x)') vvv(iz,:,iy)
!     end do
!    end do
!    do iz = 1, mp
!     do iy = 1, ny
!        write(101,'(96E14.4,1x)') zzz(iz,:,iy)
!     end do
!    end do
!    close(101)

  end subroutine read_inputdata
  
  !----------------------------------------------------------------------------
  subroutine NMF_output(itime)
    use mod_adm
    use mod_time, only: &
        cdate
    implicit none
    integer, intent(in) :: itime
    character(256) :: tmp_fname

    if(output_3DNMF) then
      if(separate_time) then
        tmp_fname=trim(coef3DNMF_fname)//trim(cdate)
        open(1,file=trim(tmp_fname),form='unformatted',access='direct', &
             recl=16*(num_zw+1)*maxl*3*num_vmode)
        write(1,rec=1) w0
        close(1)
      else
        open(1,file=trim(coef3DNMF_fname),form='unformatted',access='direct', &
             recl=16*(num_zw+1)*maxl*3*num_vmode)
        write(1,rec=itime) w0
        close(1)
      end if
    end if

  end subroutine NMF_output
  !----------------------------------------------------------------------------
  subroutine NMF_outputw(itime)
    use mod_adm

    use mod_time, only: &
        cdate
    implicit none
    integer, intent(in) :: itime
    character(256) :: tmp_fname

    if(output_3DNMF) then
      if(separate_time) then
        tmp_fname=trim(coef3DNMF_windfname)//trim(cdate)
        open(1,file=trim(tmp_fname),form='unformatted',access='direct', &
             recl=16*(num_zw+1)*maxl*3*num_vmode)
        write(1,rec=1) ww0
        close(1)
      else
        open(1,file=trim(coef3DNMF_windfname),form='unformatted',access='direct', &
             recl=16*(num_zw+1)*maxl*3*num_vmode)
        write(1,rec=itime) ww0
        close(1)
      end if
    end if

  end subroutine NMF_outputw
  !----------------------------------------------------------------------------
  subroutine global_mean(nx,ny,nz,ygrid_weight,ttt,glob)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: ny
    integer, intent(in) :: nz
    real(r8), intent(in) :: ygrid_weight(ny)
    real(r8), intent(inout) :: ttt(nz,nx,ny)
    real(r8), intent(inout) :: glob(nz)
    integer :: iy, ip

    glob(:)=0.0_r8
    do ip = 1,nz
      do iy = 1,ny
        glob(ip)=glob(ip)+sum(ttt(ip,:,iy))*ygrid_weight(iy)
      end do
    end do
    glob(:)=glob(:)/dfloat(nx)/sum(ygrid_weight)

  end subroutine global_mean
  !----------------------------------------------------------------------------
  subroutine grid_data_order 
  
! This subroutine reorders the input data so that the level 1 is the top level
! it also removes negative orography 
    
    use mod_adm
    implicit none
    integer ::  iz, iy, ix
    real(r8) :: x(nx,ny,nz)
    real(r8) :: flilla
    
    x(:,:,:) = 0.0_r8
    do iz = 1, nz
      x(:,:,iz) = u_input(:,:,nz-iz+1)
    end do
    u_input(:,:,:) = x(:,:,:)    
    x(:,:,:) = 0.0_r8
    do iz = 1, nz
      x(:,:,iz) = v_input(:,:,nz-iz+1)
    end do
    v_input(:,:,:) = x(:,:,:)
    x(:,:,:) = 0.0_r8
    do iz = 1, nz
      x(:,:,iz) = t_input(:,:,nz-iz+1)
    end do
    t_input(:,:,:) = x(:,:,:)
    x(:,:,:) = 0.0_r8
    do iz = 1, nz
      x(:,:,iz) = q_input(:,:,nz-iz+1)
    end do
    q_input(:,:,:) = x(:,:,:)      

    flilla = 1.0_r8-15
    do ix = 1, nx
     do iy = 1, ny
      if ( orog(ix,iy) < 0.0_r8 )  orog(ix,iy) = flilla
     end do
    end do
  
  end subroutine grid_data_order 
!----------------------------------------------------------------------------
  subroutine read_3d_array(filename, zzz, mp, nx, ny)
    implicit none
    character(len=*), intent(in) :: filename
    real(8), intent(out) :: zzz(mp, nx, ny)
    integer, intent(in) :: mp, nx, ny
    integer :: ip, ix, iy, iunit, ios
    character(len=100) :: line

    iunit = 10  ! Arbitrary unit number
    open(unit=iunit, file=filename, status='old', action='read', iostat=ios)

    if (ios /= 0) then
        print *, 'Error opening file.'
        return
    end if

    ! Read and ignore the header line (dimensions)
    read(iunit, '(A)', iostat=ios) line
    if (ios /= 0) then
        print *, 'Error reading file.'
        return
    end if

    do ip = 1, mp
        ! Read and ignore the slice header
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) then
            print *, 'Error reading file.'
            return
        end if

        do ix = 1, nx
            read(iunit, *, iostat=ios) (zzz(ip, ix, iy), iy = 1, ny)
            if (ios /= 0) then
                print *, 'Error reading file.'
                return
            end if
        end do
    end do

    close(iunit)
    print *, 'Contents of zzz array:'
    do ip = 1, mp
        print *, 'Slice:', ip
        do ix = 1, nx
            print *, (zzz(ip, ix, iy), iy = 1, ny)
        end do
    end do
  end subroutine read_3d_array
end module mod_normal
