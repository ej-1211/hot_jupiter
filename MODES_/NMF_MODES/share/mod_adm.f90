module mod_adm

  use mod_const, only : r8

  implicit none
  private

  !--- zonal wavenumber 
  integer, public :: ls
  !--- number of zonal wavenumber 
  integer, public :: num_zw
  !--- start zonal wavenumber
  integer, public, save :: szw
  !--- end zonal wavenumber
  integer, public, save :: ezw
  !--- start Hough mode number
  integer, public, save :: shough
  !--- end Hough mode number
  integer, public, save :: ehough
  !--- start vertical mode number
  integer, public, save :: svmode
  !--- end vertical mode number
  integer, public, save :: evmode
  !--- number of vertical modes
  integer, public, save :: num_vmode = 40
  !--- half number of grid to meridional direction for Hough functions
  integer, public :: mk
  !--- number of Hough modes in each type
  integer, public, save :: maxl
  !--- number of meridional grids for gaussian
  integer, public, save :: my = 128
  !--- number of vertical grids for gaussian
  integer, public, save :: mp = 62
  !--- number of zonal grids in original data
  integer, public, save :: nx = 256
  !--- number of meridional grids in original data
  integer, public, save :: ny = 128
  !--- number of vertical grids in original data
  integer, public, save :: nz = 91
  !--- number of data time
  integer, public, save :: nstep = 1
  !--- zonal wind (input)
  real(r8), allocatable, public, save :: u_input(:,:,:)
  !--- meridional wind (input)
  real(r8), allocatable, public, save :: v_input(:,:,:)
  !--- geopotential (input)
  real(r8), allocatable, public, save :: z_input(:,:,:)
  !--- temperature (input)
  real(r8), allocatable, public, save :: t_input(:,:,:)
  !--- specific humidity (input)
  real(r8), allocatable, public, save :: q_input(:,:,:)
  !--- Surface pressure (input)
  real(r8), allocatable, public, save :: ppp(:,:)
  !--- Orography (input)
  real(r8), allocatable, public, save :: orog(:,:)

  !--- zonal wind (input)
  real(r8), allocatable, public, save :: uuu(:,:,:)
  !--- meridional wind (input)
  real(r8), allocatable, public, save :: vvv(:,:,:)
  !--- geopotential (input)
  real(r8), allocatable, public, save :: zzz(:,:,:)
  !--- temperature (input)
  real(r8), allocatable, public, save :: ttt(:,:,:)

  !--- zonal wind after vertical expansion
  real(r8), allocatable, public, save :: vu(:,:,:)
  !--- meridional wind after vertical expansion
  real(r8), allocatable, public, save :: vv(:,:,:)
  !--- geopotential after vertical expansion
  real(r8), allocatable, public, save :: vz(:,:,:)

  !--- expansion coefficient of 3D normal mode
  double complex,   allocatable, public, save :: w0(:,:,:)
  !--- expansion coefficient of 3D normal mode, wind only part
  double complex,   allocatable, public, save :: ww0(:,:,:)
  !--- expansion coefficient of Fourier transform
  double complex,   allocatable, public, save :: wave(:,:,:,:)
  !--- hough function
  real(r8), allocatable, public, save :: hough(:,:,:,:,:)

  !--- equivalent height
  real(r8), allocatable, public, save :: evht(:)
  !--- gh
  real(r8), allocatable, public, save :: gh(:)
  !--- sqrt(gh)
  real(r8), allocatable, public, save :: sqgh(:)
  !--- vertical structure function
  real(r8), allocatable, public, save :: vsf(:,:)
  !--- vertical grid in gaussian point
  real(r8), allocatable, public, save :: vgrid(:)
  !--- gaussian weight
  real(r8), allocatable, public, save :: vgrid_weight(:)
  !--- vertical grid in sigma
  real(r8), allocatable, public, save :: vsigma(:)
  !--- vertical grid in Pressure
  real(r8), allocatable, public, save :: vpres(:)
  !--- vertical grid in Logarithm of sigma
  real(r8), allocatable, public, save :: vlnsig(:)
  !--- vertical grid in Logarithm of sigma for original data
  real(r8), allocatable, public, save :: orig_lnsigma(:)
  !--- vertical grid in Logarithm of sigma for interpolated data
  real(r8), allocatable, public, save :: model_lnsigma(:)
  !--- global mean temperature
  real(r8), allocatable, public, save :: temperature(:)
  !--- static stability
  real(r8), allocatable, public, save :: stab(:)
  !--- vertical integration
  real(r8), allocatable, public, save :: vtmp(:,:)
  !--- meridional grid in gaussian point
  real(r8), allocatable, public, save :: ygrid(:)
  !--- gaussian weight
  real(r8), allocatable, public, save :: ygrid_weight(:)
  !--- coefficient for model levels ( Hybrid sigma-p grid )
  real(r8), allocatable, public, save :: aa(:)
  real(r8), allocatable, public, save :: bb(:)
  !
  !--- longitude for input data in radius (for equally scaping grid)
  real(r8), allocatable, public, save :: lon(:)
  !--- latitude for input data in radius (for equally scaping grid)
  real(r8), allocatable, public, save :: lat(:)
  
  !real(r8), public, save :: suft = 281.1_r8
  real(r8), public, save :: suft = 2068.6_r8
  real(r8), public, save :: hstd = 8000.0_r8
  !
  ! These variables are used only when you choose analytical for 'solution_type'
  !--- Top of the atmosphere in sigma
  real(r8), public :: eps
  !--- Top of the atmosphere in Pressure
  real(r8), public :: top
  !--- Top of the atmosphere
  real(r8), public :: zt
  ! 
  !
  !--- Hough function for zonal wavenumber 0
  character(len=1),   public, save :: ks_mode = 'K'
                                            ! K = Kasahara mode (default)
                                            ! S = Shigehisa mode
  !--- Vertical coordinate for vertical expansion
  character(len=128), public, save :: zgrid   = 'Sigma'
!  !--- How to solve the vertical structure equation
  !
  ! These are the coordinates for input file

  !--- grid type for zonal direction for input
  character(len=128), public, save :: xgrid_type = 'equally'
  !--- grid type for meridional direction for input
  character(len=128), public, save :: ygrid_type = 'gaussian'
  !--- grid type for vertical direction for input
  character(len=128), public, save :: zgrid_type = 'hybrid'
  !--- Interpolation method for meridional direction
  character(len=128), public, save :: yinterpolation = ''
  !--- Interpolation method for vertical direction
  character(len=128), public, save :: zinterpolation = 'linear'
  !--- Filenames for I/O
  character(len=128), public, save :: normalcnf_fname  = './normal.cnf'
  character(len=128), public, save :: normalinvcnf_fname  = './normal_inverse.cnf'
  character(len=128), public, save :: vsfcalccnf_fname = './vsfcalc.cnf'
  character(len=128), public, save :: houghcalccnf_fname = './houghcalc.cnf'
  character(len=128), public, save :: model_level_fname = ''
  character(len=128), public, save :: equiheight_fname = ''
  character(len=128), public, save :: globalMeanTemp_fname = ''

  character(len=128), public, save :: ifile_U = ''
  character(len=128), public, save :: ifile_V = ''
  character(len=128), public, save :: ifile_Z = ''
  character(len=128), public, save :: ifile_T = ''
  character(len=128), public, save :: ifile_P = ''

  character(len=128), dimension(7), public, save :: ifile_head = ''
  character(len=128), dimension(7), public, save :: ifile_grib = ''
  character(len=128), dimension(7), public, save :: ifile_netcdf = ''
  character(len=128), public, save :: ifile_orog = ''

  character(len=128), dimension(7), public, save :: ifile_source = ''

  integer, public, save :: numoffile = 1
  ! dataformat_input can be 'grib', 'binary' and 'netcdf'
  character(len=128), public, save :: dataformat_input = ''
  ! data originating centre or model, can be 'ECMWF', 'SPEEDY'
  character(len=128), public, save :: orig = 'ECMWF'
  ! vertical coordinate of the data, can be 'hybrid', 'sigma', 'pressure' 
!  character(len=128), public, save :: vertco = 'hybrid'
  
  !--- Filename for coefficients of 3DNMF
  character(len=128), public, save :: coef3DNMF_fname = ''
  !--- Filename for coefficients of 3DNMF, only wind part
  character(len=128), public, save :: coef3DNMF_windfname = ''
  !--- Filename for inverse projection
  character(len=128), public, save :: inverse_fname = ''
  character(len=128), public, save :: ps_fname = ''
  character(len=128), public, save :: meant_fname = ''
  logical, public, save :: inv2hybrid  = .true.

  character(len=128), public, save :: source_fname = ''
  !--- save the surface pressure upon reading
  logical, public, save :: saveps = .false.
  !--- save the global temp profile on model levels 
  logical, public, save :: savemeant = .false.
  !--- save the Hough expansion coefficient only for the wind part 
  logical, public, save :: savehoughwind = .false.
    
!--- start and end indices of modes kept in inversion
!--- eig_s starting wavenumber to be filtered out
!--- eig_e ending wavenumber to be filtered out 
  integer, public, save :: vmode_s, vmode_e
  integer, public, save :: eig_n_s, eig_n_e
  integer, public, save :: wig_n_s, wig_n_e
  integer, public, save :: rot_n_s, rot_n_e
  integer, public, save :: kmode_s, kmode_e

  !--- save output of inversion as asci format data  
  logical, public, save :: saveasci = .false.
  !--- save output of inversion as asci format data  
  character(len=128), public, save :: afname = ''
  !--- save output of inversion as asci format data  
  character(len=128), public, save :: aformat = ''

  !--- save outputs in netCDF format  
  logical, public, save :: savenc = .false.
  !--- netCDF format: output name   
  character(128), public, save :: ncname = ''
        
  character(len=128), public, save :: stab_fname = ''
  character(len=128), public, save :: vsf_fname  = ''

  character(len=128), public, save :: hough_fname = ''
  character(len=128), public, save :: freq_fname  = ''
  character(len=128), public, save :: ygrid_fname = ''
  character(len=128), public, save :: vgrid_fname = ''
  character(len=128), public, save :: ofname_gmt  = ''
  character(len=128), public, save :: ofname_bin  = ''
  character(len=128), public, save :: ofname_netCDF = '' ! TJH FIXME
  character(len=128), public, save :: bin_combine = 'zonal'
  !
  logical, public, save :: SurfpresInPa    = .true.
  logical, public, save :: separate_time   = .true.
  logical, public, save :: x_interpolate   = .false.
  logical, public, save :: y_interpolate   = .false.
  logical, public, save :: z_interpolate   = .false.
  logical, public, save :: given_stability = .true.
  logical, public, save :: given_vsf       = .false.
  logical, public, save :: given_hough     = .false.
  logical, public, save :: output_gmt      = .false.
  logical, public, save :: output_global   = .true.
  logical, public, save :: output_3DNMF    = .true.
  logical, public, save :: output_inv      = .true.
  logical, public, save :: ocheck          = .true.
  !
  public :: normal_init
  public :: normal_inverse_init
  public :: houghcalc_init
  public :: vsfcalc_init
  public :: vsf_init
  public :: make_idstr
  public :: output_1d
  public :: output_2d
  public :: output_3d
  public :: input_1d
  public :: input_2d
  public :: input_3d
  public :: output_3d_gmt
  public :: output_3d_3var_gmt
  public :: output_2d_gmt  
  !
  !--- use in sub[MISC_make_idstr].
  logical, save :: NSTR_ZERO_START = .false.
  !
  !--- use in sub[MISC_make_idstr].
  integer, save :: NSTR_MAX_DIGIT  = 5
contains
  !----------------------------------------------------------------------------
  subroutine normal_init
!
! This subroutine prepares for computing 3DNMF projections
!
    use mod_const
    implicit none
    integer :: iy, m, i, iz
    character(len=128) :: fname
    real(4) :: a1, b1, z1, z2     
    namelist / normal_cnf / &
         nx, ny, nz, nstep, &
         coef3DNMF_fname,   &
         output_3DNMF,      &
         source_fname,      &
         saveps, savemeant, &
         ps_fname,          &
         meant_fname,       &
         saveasci,          &
         afname, aformat,   & 
         savenc, ncname,    &
         savehoughwind,     &
         coef3DNMF_windfname

    namelist / input_data /   &
           dataformat_input,  &
           orig,              &
           zgrid_type,        &
           ifile_orog,        &
           numoffile,         &
           ifile_head,        &
           separate_time
           
    namelist / hough_cnf /    &
         hough_fname, num_zw, &
         maxl, my

    namelist / meridional_grid / &
         ygrid_fname

    open(unit=1, file=trim(normalcnf_fname))
    read( 1,normal_cnf)
    write(*,normal_cnf)
    read( 1,input_data)
    write(*,input_data)
    read( 1,hough_cnf)
    write(*,hough_cnf)
    read( 1,meridional_grid)
    write(*,meridional_grid)
    close(1)

    allocate(lon(nx))
    allocate(lat(ny))
    allocate(ygrid(my))
    allocate(ygrid_weight(my))
    allocate(orig_lnsigma(nz))
!
! setting for vertical coordinate for input data
!
    allocate( aa(0:nz) )
    allocate( bb(0:nz) )
    open(unit=1,file=trim(model_level_fname))
    do iz = 0, nz
      read(1,*) i, a1, b1, z1, z2 
      aa(iz) = dble(a1)
      bb(iz) = dble(b1)
    end do
    write(*,*) 'Vertical coordinate coefficients:'
    do iz = 0, nz
      write(*,*) iz, aa(iz), bb(iz)
    end do

! setting the meridional coordinate for input data
!
    open(1,file=trim(ygrid_fname),form='unformatted',access='direct', &
         recl=8*my*2)
    read(1,rec=1) ygrid(:), ygrid_weight(:)
    do iy = 1,my
      write(*,'(i5,3f15.10)') iy, ygrid(iy), ygrid_weight(iy), &
           dcos(dasin(ygrid(iy)))/sum(dcos(dasin(ygrid(:))))
    end do
    write(*,*) sum(ygrid_weight)
    write(*,*) 'lat'
    write(*,*) dasin(ygrid(my:1:-1))*(-90.0_r8)*2.0_r8/pai

    allocate(w0(maxl*3,num_vmode,0:num_zw))
    allocate(ww0(maxl*3,num_vmode,0:num_zw))
    allocate(wave(0:num_zw,my,num_vmode,3))
    allocate(hough(my,maxl*3,3,num_vmode,0:num_zw))
    allocate(u_input(nx,ny,nz))
    allocate(v_input(nx,ny,nz))
    allocate(z_input(nx,ny,nz))
    allocate(t_input(nx,ny,nz))
    allocate(q_input(nx,ny,nz))
    allocate(ppp(nx,ny))
    allocate(orog(nx,ny))
    allocate(uuu(mp,nx,ny))
    allocate(vvv(mp,nx,ny))
    allocate(zzz(mp,nx,ny))
    allocate(ttt(mp,nx,ny))

    allocate(vu(nx,ny,num_vmode))
    allocate(vv(nx,ny,num_vmode))
    allocate(vz(nx,ny,num_vmode))

! reading direct access binary file
    do ls = 0, num_zw
      write(*,*) 'Reading meridional Hough structure functions for zwn ', ls
      call make_idstr( fname, trim(hough_fname), 'wn', ls)
      open(11,file=trim(fname), form='unformatted', &
           access='direct', recl=8*my*maxl*3*3, status='unknown')
      do m = 1,num_vmode
        read(11,rec=m) hough(:,:,:,m,ls)
      end do 
      close(11)
    end do    

! reading formatted ascii file
!    do ls = 0, num_zw
!        call make_idstr( fname, trim(hough_fname), 'wn', ls)
!        write(*,*) 'Reading Hough files: ', ls 
!        open(100+ls, file=trim(fname), form='formatted', action='read', status='old')
!      do m = 1,num_vmode
!       do i = 1, maxl * 3
!        do iy = 1, my
!          read(100+ls,'(3f20.12,1X)') hough(iy,i,1,m,ls),hough(iy,i,2,m,ls), &
!                                      hough(iy,i,3,m,ls)
!        end do
!       end do       
!      end do
!     close(100+ls)
!    end do
      
  end subroutine normal_init
  !----------------------------------------------------------------------------
  subroutine normal_inverse_init
    use mod_const
    implicit none
    integer :: iy, m, iz, i
    character(len=128) :: fname
    real(4) :: a1, b1, z1, z2   
    
    namelist / normal_cnf_inverse / &
         nx, ny, nz, nstep,         &
         coef3DNMF_fname,           &
         inverse_fname,             &
         inv2hybrid,                &
         ps_fname,                  &
         meant_fname,               &
         saveasci, afname, aformat, &
         savenc, ncname

    namelist / hough_cnf /    &
         hough_fname, num_zw, &
         maxl, my,            &
         shough, ehough

    namelist / meridional_grid / &
         ygrid_fname

    namelist / filter_cnf / &
         vmode_s, vmode_e,      &
         eig_n_s, eig_n_e,      &
         wig_n_s, wig_n_e,      &
         rot_n_s, rot_n_e,      &
         kmode_s, kmode_e
                  
    open(unit=1, file=trim(normalinvcnf_fname))
    read(1,normal_cnf_inverse)
    write(*,normal_cnf_inverse)
    read(1,hough_cnf)
    write(*,hough_cnf)
    read(1,meridional_grid)
    write(*,meridional_grid)
    read(1,filter_cnf)
    write(*,filter_cnf)
    close(1)
 
    allocate(lon(nx))
    allocate(lat(ny))
    allocate(ygrid(my))
    allocate(ygrid_weight(my))
    allocate(hough(my,maxl*3,3,num_vmode,0:num_zw))
    allocate(w0(maxl*3,num_vmode,0:num_zw))
    allocate(wave(0:num_zw,my,num_vmode,3))
    allocate(vu(nx,my,num_vmode))
    allocate(vv(nx,my,num_vmode))
    allocate(vz(nx,my,num_vmode))
    allocate(u_input(nx,my,mp))
    allocate(v_input(nx,my,mp))
    allocate(z_input(nx,my,mp))

! setting for vertical coordinate for input data
!
    allocate( aa(0:nz) )
    allocate( bb(0:nz) )
    open(unit=1,file=trim(model_level_fname))
    do iz = 0, nz
      read(1,*) i, a1, b1, z1, z2   
      aa(iz) = dble(a1) 
      bb(iz) = dble(b1) 
    end do
    write(*,*) 'Vertical coordinate coefficients from ',trim(model_level_fname),':'
    do iz = 0, nz
      write(*,*) iz, aa(iz), bb(iz)
    end do
!
    write(*,*) 'Gaussian grid coefficients from ',trim(ygrid_fname),':'
    open(1,file=trim(ygrid_fname),form='unformatted',access='direct', &
         recl=8*my*2)
    read(1,rec=1) ygrid(:), ygrid_weight(:)
    do iy = 1,my
      write(*,'(i5,3f15.10)') iy, ygrid(iy), ygrid_weight(iy), &
           dcos(dasin(ygrid(iy)))/sum(dcos(dasin(ygrid(:))))
    end do
    write(*,*) sum(ygrid_weight)

    do ls = 0, num_zw
      write(*,*) 'Reading meridional Hough structure functions for zwn ', ls
      call make_idstr( fname, trim(hough_fname), 'wn', ls)
      open(11,file=trim(fname), form='unformatted', &
           access='direct', recl=8*my*maxl*3*3, status='unknown')
      do m = 1,num_vmode
        read(11,rec=m) hough(:,:,:,m,ls)
      end do
      close(11)
    end do

  end subroutine normal_inverse_init
  !----------------------------------------------------------------------------
  subroutine houghcalc_init
!
! This subroutine prepares for computing Hough functions
!
    implicit none
    character(len=128) :: fname
    integer :: iy

    namelist / houghcalc_cnf / &
        maxl, my,              &
        szw, ezw,              &
        ks_mode,               &
        freq_fname,            &
        ocheck

    namelist / output /    &
        output_gmt,        &
        ofname_gmt,        &
        ofname_bin,        &
        ofname_netCDF,     &
        bin_combine

    namelist / vsf_cnf /   &
        equiheight_fname,  &
        num_vmode

    namelist / meridional_grid / &
         ygrid_fname

    open(unit=1, file=trim(houghcalccnf_fname))
    read(1,houghcalc_cnf)
    write(*,houghcalc_cnf)
    read(1,meridional_grid)
    write(*,meridional_grid)
    read(1,vsf_cnf)
    write(*,vsf_cnf)
    read(1,output)
    write(*,output)
    close(1)

    mk = my/2

    allocate( evht(num_vmode) )
    
    open(unit=1,file=trim(equiheight_fname),form='unformatted', &
         access='direct', recl=8*num_vmode)
    read(1,rec=1) evht(:)
    close(1)
    write(*,*) evht(:)

    allocate(ygrid(my))
    allocate(ygrid_weight(my))
   
    write(*,*) 'Gaussian grid: '
    open(1,file=trim(ygrid_fname),form='unformatted',access='direct', &
         recl=8*my*2)
    read(1,rec=1) ygrid(:), ygrid_weight(:)
    close(1)
    do iy = 1,my
      write(*,*) iy, ygrid(iy), ygrid_weight(iy)
    end do

    if(output_gmt) then
      if( trim(bin_combine) == 'zonal' ) then
       do iy = szw, ezw 
        call make_idstr( hough_fname, trim(ofname_gmt), 'wn', iy)
        write(*,*) iy, trim(hough_fname)
        open(100+iy, file=trim(hough_fname), status='unknown')
!        close(100+iy)
       end do
      end if
     end if
      
    if(output_gmt) then
      fname=trim(ofname_gmt)//'_U.dat'
      open(unit=101,file=trim(fname)) 
      fname=trim(ofname_gmt)//'_V.dat'
      open(unit=102,file=trim(fname)) 
      fname=trim(ofname_gmt)//'_Z.dat'
      open(unit=103,file=trim(fname)) 
    end if

  end subroutine houghcalc_init
  !----------------------------------------------------------------------------
  subroutine vsfcalc_init
!
! This subroutine prepares for computing Vertical structure functions.
!
    use mod_const, only : &
        ps
    implicit none

    namelist / vsfcalc_cnf / &
        equiheight_fname,    &
        vsf_fname,           &
        vgrid_fname,         &
        num_vmode,           &
        mp, zgrid,           &
        top,                 &
        suft, hstd,          &
        given_stability,     &
        globalMeanTemp_fname,&
        given_vsf,           &
        stab_fname,          &
        ocheck

    open(unit=1, file=trim(vsfcalccnf_fname))
    read(1,vsfcalc_cnf)
    write(*,vsfcalc_cnf)
    close(1)
!
    allocate( evht(num_vmode) )
    allocate(  vsf(mp,num_vmode) )
    allocate(  vgrid(mp) )
    allocate(  vgrid_weight(mp) )
    allocate(  vsigma(mp) )
    allocate(  vpres(mp) )
    allocate(  stab(mp) )
!
    !open(unit=1,file=trim(vgrid_fname), form='unformatted',&
    !       access='direct',recl=mp*8)
    !read(1,rec=1) vgrid(:)
    open(unit=1,file=trim(vgrid_fname), form='formatted',&
           status='old', action='read',access='sequential')
    read(1,*) vgrid(:)
    close(1)

    if(given_stability) then
      open(unit=1,file=trim(stab_fname), form='formatted',&
           status='old', action='read',access='sequential')
      read(1,*) stab(:)

      !open(unit=1,file=trim(stab_fname), form='unformatted',&
      !     access='direct',recl=mp*8)
      !read(1,rec=1) stab(:)
      close(1)      
    else
      allocate(temperature(mp))
      open(unit=1,file=trim(globalMeanTemp_fname),  form='unformatted',&
           access='sequential')
      read(1) temperature
      close(1)
    end if
    write(*,*) 'Vgrid as read from ',trim(vgrid_fname),': ', vgrid
    write(*,*) 'Stability: ', stab
    
  end subroutine vsfcalc_init

  !----------------------------------------------------------------------------

  subroutine vsf_init
  !
  ! This subroutine prepares for vertical expansion when you compute
  ! 3D NMF projections.
  ! 
    use mod_const, only : ps, g
    implicit none
    integer :: ip, m

    namelist / vsf_cnf /   &
        equiheight_fname,  &
        model_level_fname, &
        vsf_fname,         &
        vgrid_fname,       &
        num_vmode,         &
        mp, zgrid,         &
        top,               &
        given_vsf,         &
        stab_fname

    open(unit=1, file=trim(normalcnf_fname))
    read(1,vsf_cnf)
    write(*,vsf_cnf)
    close(1)

    allocate( evht(num_vmode) )
    allocate(   gh(num_vmode) )
    allocate( sqgh(num_vmode) )
    allocate(  vsf(mp,num_vmode) )
    allocate(  vtmp(mp,num_vmode) )
    allocate(  vgrid(mp) )
    allocate(  vgrid_weight(mp) )
    allocate(  model_lnsigma(mp) )
    allocate(  vsigma(mp) )
    allocate(  vpres(mp) )
    allocate(  stab(mp) )
!
    if(given_vsf) then ! read equivalent depths and vertical structure functions
      write(*,*) ' Reading equivalent Height'
      call input_1d(trim(equiheight_fname), evht, num_vmode, 1, .true., 8)
      write(*,*) evht
      write(*,*) ' Reading vertical structure functions'  
      call input_2d(trim(vsf_fname), vsf, mp, num_vmode, 1, .true., 8)  

    end if
    gh(:)=evht(:)*g
    sqgh(:)=dsqrt(gh(:))
!
    !open(unit=1,file=trim(vgrid_fname), form='unformatted',&
    !       access='direct',recl=mp*8)
    !read(1,rec=1) vgrid(:)
    !close(1)
    open(unit=1,file=trim(vgrid_fname), form='formatted',&
           status='old', action='read',access='sequential')
    read(1,*) vgrid(:)
    close(1)
    write(*,*) 'Vgrid sigma values from ',trim(vgrid_fname),': '
    write(*,*) vgrid
    ! vgrid is by default in the input file sorted from the bottom upward

    do ip = 2,mp
      vgrid_weight(ip)=vgrid(ip-1)-vgrid(ip)
    end do
    vgrid_weight(1)=2.0_r8*(1.0_r8-vgrid(1))
    model_lnsigma(:)=-dlog(vgrid(:))

    ! vertical structure functions saved for later use in the program vertical_expansion
    do m = 1,num_vmode
     do ip = 1,mp
      vtmp(ip,m) = vsf(mp-ip+1,m) !now sorted top to bottom
     end do
    end do

!    write(*,*) 'VSF: ip,vsf(ip,1),vsf(ip,2), vsf(ip,3)'
!    do ip = 1,mp
!      write(*,'(i4,3f15.4)') ip,vsf(mp-ip+1,1),vsf(mp-ip+1,2), vsf(mp-ip+1,3)
!    end do
!
  end subroutine vsf_init
  !----------------------------------------------------------------------------
  subroutine make_idstr( &
       str,                   & !--- OUT : file name ! 05/11/15 M.Satoh
       headstr,               & !--- IN : header name
       ext,                   & !--- IN : extention string
       rank )                   !--- IN : ID number
    !
    implicit none
    !
    character(len=*), intent(out) :: str     !--- strings 
    character(len=*), intent(in) :: headstr  !--- header strings
    character(len=*), intent(in) :: ext      !--- extention( eg. ***.dat ) 
    integer, intent(in) :: rank          !--- number(+1)
    !
    character(len=128) :: cnum
    character(len=NSTR_MAX_DIGIT) cnum1
    !
    if(NSTR_ZERO_START) then
       write(cnum,'(I128.128)') rank-1
    else
       write(cnum,'(I128.128)') rank
    end if
    cnum1(1:NSTR_MAX_DIGIT) = cnum(128-(NSTR_MAX_DIGIT-1):128)
    str=trim(headstr)//'.'//trim(ext)//cnum1
    !
  end subroutine make_idstr

  !-----------------------------------------------------------------------------

  subroutine output_1d(fname, var, nobs, recnum, direct_access, output_size)
    implicit none
    character(len=*), intent(in) :: fname
    real(r8),         intent(in) :: var(nobs)
    integer,          intent(in) :: nobs
    integer,          intent(in) :: recnum
    integer,          intent(in) :: output_size
    logical,          intent(in) :: direct_access

    if(direct_access) then
      open(10,file=trim(fname), form='unformatted', &
           access='direct',recl=output_size*nobs)
      write(10,rec=recnum) var(:)
      close(10)
    else
      open(10,file=trim(fname), form='unformatted', access='sequential')
      write(10) var(:)
      close(10)
    end if

  end subroutine output_1d

  !-----------------------------------------------------------------------------

  subroutine output_2d(fname, var, nobs1, nobs2, recnum, &
                       direct_access, output_size)
    implicit none
    character(len=*), intent(in) :: fname
    real(r8),         intent(in) :: var(nobs1,nobs2)
    integer,          intent(in) :: nobs1
    integer,          intent(in) :: nobs2
    integer,          intent(in) :: recnum
    integer,          intent(in) :: output_size
    logical,          intent(in) :: direct_access

    if(direct_access) then
      open(10,file=trim(fname), form='unformatted', &
           access='direct',recl=output_size*nobs1*nobs2)
      write(10,rec=recnum) var(:,:)
      close(10)
    else
      open(10,file=trim(fname), form='unformatted', access='sequential')
      write(10) var(:,:)
      close(10)
    end if

  end subroutine output_2d

  !-----------------------------------------------------------------------------

  subroutine output_3d(fname, var, nobs1, nobs2, nobs3, recnum, &
                       direct_access, output_size)
    implicit none
    character(len=*), intent(in) :: fname
    real(r8),         intent(in) :: var(nobs1,nobs2,nobs3)
    integer,          intent(in) :: nobs1
    integer,          intent(in) :: nobs2
    integer,          intent(in) :: nobs3
    integer,          intent(in) :: recnum
    integer,          intent(in) :: output_size
    logical,          intent(in) :: direct_access

    if(direct_access) then
      open(10,file=trim(fname), form='unformatted', &
           access='direct',recl=output_size*nobs1*nobs2*nobs3)
      write(10,rec=recnum) var(:,:,:)
      close(10)
    else
      open(10,file=trim(fname), form='unformatted', access='sequential')
      write(10) var(:,:,:)
      close(10)
    end if

  end subroutine output_3d

  !-----------------------------------------------------------------------------

  subroutine input_1d(fname, var, nobs, recnum, direct_access, input_size)
    implicit none
    character(len=*), intent(in)    :: fname
    real(r8),         intent(inout) :: var(nobs)
    integer,          intent(in)    :: nobs
    integer,          intent(in)    :: recnum
    integer,          intent(in)    :: input_size
    logical,          intent(in)    :: direct_access

    if(direct_access) then
      open(10,file=trim(fname), form='unformatted', &
           access='direct',recl=input_size*nobs)
      read(10,rec=recnum) var(:)
      close(10)
    else
      open(10,file=trim(fname), form='unformatted', access='sequential')
      read(10) var(:)
      close(10)
    end if

  end subroutine input_1d
  !-----------------------------------------------------------------------------
  subroutine input_2d(fname, var, nobs1, nobs2, &
                      recnum, direct_access, input_size)
    implicit none
    character(len=*), intent(in) :: fname
    real(r8), intent(inout) :: var(nobs1,nobs2)
    integer, intent(in) :: nobs1
    integer, intent(in) :: nobs2
    integer, intent(in) :: recnum
    integer, intent(in) :: input_size
    logical, intent(in) :: direct_access
    integer(8) :: nobs
    real(4) :: temp(nobs1, nobs2)      

    nobs=nobs1*nobs2

    if(direct_access) then
      if(input_size==8) then
        open(10,file=trim(fname), form='unformatted', &
             access='direct',recl=input_size*nobs)
        read(10,rec=recnum) var(:,:)
        close(10)
      else if(input_size==4) then
        open(10,file=trim(fname), form='unformatted', &
             access='direct',recl=input_size*nobs)
        read(10,rec=recnum) temp(:,:)
        var(:,:) = dble(temp(:,:))
        close(10)
      end if

    else
      open(10,file=trim(fname), form='unformatted', access='sequential')
      read(10) var(:,:)
      close(10)
    end if
      
  end subroutine input_2d
  !-----------------------------------------------------------------------------
  subroutine input_3d(fname, var, nobs1, nobs2, nobs3, &
                      recnum, direct_access, input_size)
    implicit none
    character(len=*), intent(in)    :: fname
    real(r8),         intent(inout) :: var(nobs1,nobs2,nobs3)
    integer,          intent(in)    :: nobs1
    integer,          intent(in)    :: nobs2
    integer,          intent(in)    :: nobs3
    integer,          intent(in)    :: recnum
    integer,          intent(in)    :: input_size
    logical,          intent(in)    :: direct_access

    integer(8) :: nobs
    real(4) :: temp(nobs1, nobs2, nobs3)      

    nobs=nobs1*nobs2*nobs3

    if(direct_access) then
      if(input_size==8) then
        open(10,file=trim(fname), form='unformatted', &
             access='direct',recl=input_size*nobs)
        read(10,rec=recnum) var(:,:,:)
        close(10)
      else if(input_size==4) then
        open(10,file=trim(fname), form='unformatted', &
             access='direct',recl=input_size*nobs)
        read(10,rec=recnum) temp(:,:,:)
        var(:,:,:) = dble(temp(:,:,:))
        close(10)
      end if
    else
      open(10,file=trim(fname), form='unformatted', access='sequential')
      read(10) var(:,:,:)
      close(10)
    end if

  end subroutine input_3d
  !----------------------------------------------------------------------------- 
    subroutine output_3d_gmt(afname, aformat, ccdate, var, nobs1, nobs2, nobs3)
    implicit none
    character(len=*),  intent(in)    :: afname
    character(len=*),  intent(in)    :: aformat
    character(len=15), intent(in)    :: ccdate
    real(r8),          intent(inout) :: var(nobs1,nobs2,nobs3)
    integer,           intent(in)    :: nobs1
    integer,           intent(in)    :: nobs2
    integer,           intent(in)    :: nobs3
    integer :: i2, i3, ierr
    character(len=256) :: tmp_fname
    
    tmp_fname=trim(afname)//trim(ccdate)
    open(unit=101,file=trim(tmp_fname),status='unknown', &
         action='write',iostat=ierr)
    if (ierr /= 0) then
      print *, 'ERROR: Opening gp 3d file ',trim(tmp_fname),' failed with err=', ierr
      stop
    end if   
    write(*,*) 'Writing gp file ', trim(tmp_fname)
    do i3 = 1, nobs3
     do i2 = 1, nobs2
        write(101,fmt=trim(aformat)) var(:,i2,i3)
     end do
    end do
    close(101)
    
  end subroutine output_3d_gmt
  !----------------------------------------------------------------------------- 
  subroutine output_3d_3var_gmt(afname, aformat, ccdate, &
             var1, var2, var3, nobs1, nobs2, nobs3)
    implicit none
    character(len=*), intent(in) :: afname
    character(len=*), intent(in) :: aformat
    character(len=15), intent(in) :: ccdate
    real(r8), intent(inout) :: var1(nobs1,nobs2,nobs3)
    real(r8), intent(inout) :: var2(nobs1,nobs2,nobs3)
    real(r8), intent(inout) :: var3(nobs1,nobs2,nobs3)
    integer, intent(in) :: nobs1
    integer, intent(in) :: nobs2
    integer, intent(in) :: nobs3
    integer :: i2,i3, ierr
    character(len=256) :: tmp_fname

    tmp_fname=trim(afname)//trim(ccdate)
    open(unit=101,file=trim(tmp_fname),status='unknown', &
         action='write',iostat=ierr)
    if (ierr /= 0) then
      print *, 'ERROR: Opening gp 3d file ',trim(tmp_fname),' failed with err=', ierr
      stop
    end if   
    write(*,*) 'Writing gp file ', trim(tmp_fname)
    do i3 = 1, nobs3
     do i2 = 1, nobs2
        write(101,fmt=trim(aformat)) var1(:,i2,i3)
     end do
    end do
    do i3 = 1, nobs3
     do i2 = 1, nobs2
        write(101,fmt=trim(aformat)) var2(:,i2,i3)
     end do
    end do
    do i3 = 1, nobs3
     do i2 = 1, nobs2
        write(101,fmt=trim(aformat)) var3(:,i2,i3)
     end do
    end do
    close(101)
    
  end subroutine output_3d_3var_gmt

 !-----------------------------------------------------------------------------

 subroutine output_2d_gmt(afname, aformat, var, nobs1, nobs2)
    implicit none
    character(len=*), intent(in)    :: afname
    character(len=*), intent(in)    :: aformat
    real(r8),         intent(inout) :: var(nobs1,nobs2)
    integer,          intent(in)    :: nobs1
    integer,          intent(in)    :: nobs2
    integer :: i2

    open(unit=102,file=trim(afname))
      
    do i2 = 1, nobs2
        write(102,fmt=trim(aformat)) var(:,i2)
    end do
    close(102)
    
  end subroutine output_2d_gmt
 !-----------------------------------------------------------------------------
       
end module mod_adm
