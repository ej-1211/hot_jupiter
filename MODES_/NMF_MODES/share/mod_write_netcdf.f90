module mod_write_netcdf

  use netcdf
  use mod_adm, only : my        ! meridional grid size
  use mod_adm, only : maxl      ! number of meridional modes
  use mod_adm, only : num_vmode ! number of vertical structure function modes
  use mod_adm, only : szw       ! first/starting zonal wavenumber
  use mod_adm, only : ezw       ! last/ending    zonal wavenumber
  use mod_adm, only : ygrid, ygrid_weight  ! gaussian grid/weights

  use mod_adm, only : nx, ny, nz, vsf, evht, stab, vgrid, vgrid_weight, &
                      vpres, vsigma, wig_n_s, wig_n_e, vmode_s, vmode_e, &
                      rot_n_s, rot_n_e, kmode_s, kmode_e, eig_n_s, eig_n_e, &
                      savehoughwind, num_zw, equiheight_fname, ofname_bin, &
                      given_stability, given_vsf, hstd, ocheck, stab_fname, top, &
                      vgrid_fname, vsf_fname, zgrid, make_idstr, freq_fname, &
                      ks_mode

  use mod_const
  use mod_time

  implicit none
  private

  public :: write_netcdf_inv
  public :: write_netcdf_proj
  public :: write_netcdf_Tmean
  public :: write_netcdf_lnps
  public :: write_netcdf_vsf
  public :: initialize_hough_output
  public :: write_netcdf_hough
  public :: finalize_hough_output
  public :: check

  ! nc_table will relate the filename to the netcdf id.
  ! since each netcdf file will be left open and written to many times,
  ! we need a convenient way to figure out which ncid to write to
  ! given each wavenumber.

  type hough_nc_table
     private
     character(len=512) :: filename
     integer :: ncid
  end type hough_nc_table

  type(hough_nc_table), allocatable, dimension(:) :: netcdf_table
  logical :: debug = .true.

  integer, parameter :: EWR = 3
  integer, parameter :: EIG_INDEX = 1
  integer, parameter :: WIG_INDEX = 2
  integer, parameter :: BAL_INDEX = 3

  integer, parameter :: UVZ = 3
  integer, parameter :: U_INDEX = 1
  integer, parameter :: V_INDEX = 2
  integer, parameter :: Z_INDEX = 3

  !----------------------------------------------------------------------------
  contains
  !----------------------------------------------------------------------------

  subroutine write_netcdf_inv(outfile, ccdate, var1, var2, var3, nobs1, nobs2, nobs3)

  character(len=*),  intent(in)    :: outfile
  character(len=15), intent(in)    :: ccdate
  integer,           intent(in)    :: nobs1
  integer,           intent(in)    :: nobs2
  integer,           intent(in)    :: nobs3
  real(r8),          intent(inout) :: var1(nobs1,nobs2,nobs3)
  real(r8),          intent(inout) :: var2(nobs1,nobs2,nobs3)
  real(r8),          intent(inout) :: var3(nobs1,nobs2,nobs3)

  ! Local variables
  integer, PARAMETER :: Nfields=3, Ndim=4, nobs4=1

  ! NetCDF - IDs and handlers
  integer                  :: ncid, I, J
  integer                  :: t_dimid, x_dimid, y_dimid, z_dimid
  integer                  :: t_varid, x_varid, y_varid, z_varid
  integer                  :: u_varid, v_varid, P_varid
  integer, dimension(Ndim) :: dimids4D

  ! NetCDF attribute and storage info
  character(len=128), dimension(3)       :: infname
  character(len=128), dimension(Nfields) :: vrname, lvrname, uvrname
  character(len=128), dimension(Ndim)    :: dmname, ldmname, udmname
  character(len=256) :: outnc, filtername

  ! TIME handling
  integer,dimension(8) :: ivals, tvals
  character(len=256)   :: curr_date, st_date

  ! Output fields (including UNLIMITED time dimension)
  real(r8) :: out_lon(nobs1), out_lat(nobs2), out_lev(nobs3)
  real(r8) :: nhours
  real(r8) :: var4D1(nobs1,nobs2,nobs3,nobs4)
  real(r8) :: var4D2(nobs1,nobs2,nobs3,nobs4)
  real(r8) :: var4D3(nobs1,nobs2,nobs3,nobs4)

  !---------------------------------------------------------------------------!
  !                               Preparations
  !---------------------------------------------------------------------------!

  ! mod_adm:ygrid(my)        are the gaussian points
  ! mod_adm:ygrid_weight(my) are the gaussian weights
  ! Define coodinate arrays (Gaussian grid + sigma levels)

  out_lon = (/ (I, I = 0, nx-1) /)*(360.0_r8/nx)

  !out_lat = (/ (J, J = ny-1, 0, -1) /)*(180.0_r8/ny) - 90.0_r8*(1 - 1.0_r8/ny)   ! NP to SP
  !out_lat = (/ (J, J = 0, ny-1) /)*(180./ny) - 90*(1 - 1./ny)       ! SP to NP

  ! this is the same formula used by mod_adm:normal_init() to calculate the gaussian lats
  out_lat = dasin(ygrid(my:1:-1))*(-90.0_r8)*2.0_r8/pai   ! [South-to-North]
  out_lat = dasin(ygrid(1:my))*(-90.0_r8)*2.0_r8/pai      ! [North-to-South]

  out_lev = vgrid(nz:1:-1)                                          ! Top to bottom

  ! Add UNLIMITED time dimension to variables
  ! > NP to SP <
  var4D1(:,:,:,1) = var1
  var4D2(:,:,:,1) = var2
  var4D3(:,:,:,1) = var3
  ! > SP to NP <
  !do J = 1, ny
  !  var4D1(:,J,:,1) = var1(:,ny-J+1,:)
  !  var4D2(:,J,:,1) = var2(:,ny-J+1,:)
  !  var4D3(:,J,:,1) = var3(:,ny-J+1,:)
  !end do

  ! Define time value
  ! > current date and time
  call date_and_time(VALUES=tvals)
  ! > save as current date in specified format
  write(curr_date, "(I2.2,A1,I2.2,A1,I4,A4,I2.2,A1,I2.2,A1,I2.2)") tvals(3), "-", tvals(2), "-", tvals(1), &
                  "   ", tvals(5), ":", tvals(6), ":", tvals(7)

  ! > Time value (in hours) since 1st dataset  (CHANGE LATER)
  ivals = (/ syear, smon, sday, 60, shour, smins, ssec, 0/)
  write(st_date, "(A1,I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)") " ", ivals(1), "-", ivals(2), "-", ivals(3), &
                  " ", ivals(5), ":", ivals(6), ":", ivals(7)
  nhours = (itime-1)*dt/3600

  ! netCDF filename specification
  outnc = trim(outfile)//trim(ccdate)//trim(".nc")

  ! Attribute names
  ! > variables
  vrname(1)  = "u"
  vrname(2)  = "v"
  vrname(3)  = "Z"
  lvrname(1) = "U component of wind"
  lvrname(2) = "V component of wind"
  lvrname(3) = "modified geopotential heigth, Z = P/g ... P = gh + RT_0 ln(p_s)"
  uvrname(1) = "m s**-1"
  uvrname(2) = "m s**-1"
  uvrname(3) = "gpm"

  ! > coordinates
  dmname(1)  = "lon"
  dmname(2)  = "lat"
  dmname(3)  = "lev"
  dmname(4)  = "time"
  ldmname(1) = "longitude"
  ldmname(2) = "latitude"
  ldmname(3) = "sigma"
  ldmname(4) = "time"
  udmname(1) = "degrees_east"
  udmname(2) = "degrees_north"
  udmname(3) = "~"
  udmname(4) = trim("hours since")//trim(st_date)

  ! > global
  if((eig_n_e .gt. eig_n_s) .or. (wig_n_e .gt. wig_n_s) .or. (rot_n_e .gt. rot_n_s)) then
    write(infname(1), "(A31,I3,A4,I3,A5,I3,A4,I3,A9,I3,A4,I3)") "Filtered meridional modes: EIG ", &
              eig_n_s," to ",eig_n_e, ", WIG ",wig_n_s," to ",wig_n_e, &
             " and BAL ",rot_n_s," to ",rot_n_e
  else
    infname(1) = "No filtering of meridional modes"
  end if

  if(vmode_e .gt. vmode_s) then
    write(infname(2), "(A24,I3,A4,I3)") "Filtered vertical modes ", vmode_s," to ",vmode_e
  else
    infname(2) = "No filtering of vertical modes"
  end if

  if(kmode_e .gt. kmode_s) then
    write(infname(3), "(A27,I3,A4,I3)") "Filtered zonal wavenumbers ", kmode_s," to ",kmode_e
  else
    infname(3) = "No filtering of zonal modes"
  end if

  filtername = "Modes filtered: "//NEW_LINE('\n')//trim(infname(1))//&
              NEW_LINE('\n')//trim(infname(2))//NEW_LINE('\n')//trim(infname(3))

  !---------------------------------------------------------------------------!
  !                            Write out to netCDF
  !---------------------------------------------------------------------------!

  ! <1> Open NetCDF file
  call check(nf90_create(outnc,NF90_CLOBBER, ncid))

  ! <2> Define dimensions
  call check(nf90_def_dim(ncid, dmname(1), nobs1, x_dimid))
  call check(nf90_def_dim(ncid, dmname(2), nobs2, y_dimid))
  call check(nf90_def_dim(ncid, dmname(3), nobs3, z_dimid))
  call check(nf90_def_dim(ncid, dmname(4), NF90_UNLIMITED, t_dimid))

  ! <3> Define coordinate variables
  call check(nf90_def_var(ncid, dmname(1), NF90_DOUBLE, x_dimid, x_varid))
  call check(nf90_def_var(ncid, dmname(2), NF90_DOUBLE, y_dimid, y_varid))
  call check(nf90_def_var(ncid, dmname(3), NF90_DOUBLE, z_dimid, z_varid))
  call check(nf90_def_var(ncid, dmname(4), NF90_DOUBLE, t_dimid, t_varid))
  ! coordinate array
  dimids4D = (/ x_dimid, y_dimid, z_dimid, t_dimid /)

  ! <4> Define variables
  call check(nf90_def_var(ncid, vrname(1), NF90_FLOAT, dimids4D, u_varid))
  call check(nf90_def_var(ncid, vrname(2), NF90_FLOAT, dimids4D, v_varid))
  call check(nf90_def_var(ncid, vrname(3), NF90_FLOAT, dimids4D, P_varid))

  ! <5> Assign attributes
  call check(nf90_put_att(ncid, x_varid, "long_name", ldmname(1)))
  call check(nf90_put_att(ncid, x_varid, "units",     udmname(1)))

  call check(nf90_put_att(ncid, y_varid, "long_name", ldmname(2)))
  call check(nf90_put_att(ncid, y_varid, "units",     udmname(2)))

  call check(nf90_put_att(ncid, z_varid, "long_name", ldmname(3)))
  call check(nf90_put_att(ncid, z_varid, "units",     udmname(3)))

  call check(nf90_put_att(ncid, t_varid, "long_name", ldmname(4)))
  call check(nf90_put_att(ncid, t_varid, "units",     udmname(4)))
  call check(nf90_put_att(ncid, t_varid, "calendar", "standard"))

  call check(nf90_put_att(ncid, u_varid, "long_name", lvrname(1)))
  call check(nf90_put_att(ncid, u_varid, "units",     uvrname(1)))
  call check(nf90_put_att(ncid, u_varid, "grid_type", "gaussian"))

  call check(nf90_put_att(ncid, v_varid, "long_name", lvrname(2)))
  call check(nf90_put_att(ncid, v_varid, "units",     uvrname(2)))
  call check(nf90_put_att(ncid, v_varid, "grid_type", "gaussian"))

  call check(nf90_put_att(ncid, P_varid, "long_name", lvrname(3)))
  call check(nf90_put_att(ncid, P_varid, "units",     uvrname(3)))
  call check(nf90_put_att(ncid, P_varid, "grid_type", "gaussian"))

  ! Global attributes
  call write_global_attributes(ncid, "inverse fields", curr_date, filtername)

  ! <6>  End define mode
  call check(nf90_enddef(ncid))

  ! <7> Write data
  call check(nf90_put_var(ncid, t_varid, nhours))
  call check(nf90_put_var(ncid, z_varid, out_lev))
  call check(nf90_put_var(ncid, y_varid, out_lat))
  call check(nf90_put_var(ncid, x_varid, out_lon))

  call check(nf90_put_var(ncid, u_varid, sngl(var4D1)))
  call check(nf90_put_var(ncid, v_varid, sngl(var4D2)))
  call check(nf90_put_var(ncid, P_varid, sngl(var4D3)))

  ! <8> close NetCDF file
  call check(nf90_close(ncid))

  ! If we got this far, everything worked as expected. Yipee!
  write(*,*) ""
  write(*,*) "*** SUCCESSFULLY written inversion fields to NetCDF file ***"
  write(*,*) "*** <",trim(outnc),"> ***"
  write(*,*) ""

  end subroutine write_netcdf_inv


!-------------------------------------------------------------------------------


  subroutine write_netcdf_proj(outfile, ccdate, var, nobs1, nobs2, nobs3)

  ! Subroutine given arguments TJH FIXME ... nobs123 should be before var in the call
  character(len=*),  intent(in)    :: outfile
  character(len=15), intent(in)    :: ccdate
  integer,           intent(in)    :: nobs1
  integer,           intent(in)    :: nobs2
  integer,           intent(in)    :: nobs3
  double complex,    intent(inout) :: var(nobs1*3,nobs2,nobs3)

  ! #
  integer, PARAMETER :: Nfields=3, Ndim=5
  integer, PARAMETER :: nobs4=2, nobs5=1

  ! NetCDF - IDs and handlers
  integer                  :: ncid, I, J, K
  integer                  :: t_dimid, k_dimid, n_dimid, m_dimid, c_dimid
  integer                  :: t_varid, k_varid, n_varid, m_varid
  integer                  :: EIG_varid, WIG_varid, BAL_varid
  integer, dimension(Ndim) :: dimids5D

  ! NetCDF attribute and storage info
  character(len=128), dimension(Nfields) :: vrname, lvrname, uvrname
  character(len=128), dimension(Ndim)    :: dmname, ldmname, udmname
  character(len=256)                     :: outnc, infoname

  ! TIME handling
  integer,dimension(8) :: ivals, tvals
  character(len=256)   :: curr_date, st_date

  ! Output fields (including UNLIMITED time dimension)
  real(r8) :: out_n(nobs1), out_m(nobs2), out_k(nobs3)
  real(r8) :: nhours
  real(r8) :: var5D1(nobs1,nobs2,nobs3,nobs4,nobs5)
  real(r8) :: var5D2(nobs1,nobs2,nobs3,nobs4,nobs5)
  real(r8) :: var5D3(nobs1,nobs2,nobs3,nobs4,nobs5)

  !---------------------------------------------------------------------------!
  !                               Preparations
  !---------------------------------------------------------------------------!

  ! Define coodinate arrays (Gaussian grid + sigma levels)
  out_k = (/ (I, I = 0, num_zw) /)
  out_n = (/ (J, J = 0, maxl-1) /)
  out_m = (/ (K, K = 1, num_vmode) /)

  ! Add UNLIMITED time dimension to variables
  ! Divide in REAL and COMPLEX part
  var5D1(:,:,:,1,1) = DREAL(var((0*nobs1+1):(1*nobs1),:,:))
  var5D1(:,:,:,2,1) = IMAG (var((0*nobs1+1):(1*nobs1),:,:))
  var5D2(:,:,:,1,1) = DREAL(var((1*nobs1+1):(2*nobs1),:,:))
  var5D2(:,:,:,2,1) = IMAG (var((1*nobs1+1):(2*nobs1),:,:))
  var5D3(:,:,:,1,1) = DREAL(var((2*nobs1+1):(3*nobs1),:,:))
  var5D3(:,:,:,2,1) = IMAG (var((2*nobs1+1):(3*nobs1),:,:))

  ! Define time value
  ! > current date and time
  call date_and_time(VALUES=tvals)
  ! > save as current date in specified format
  write(curr_date, "(I2.2,A1,I2.2,A1,I4,A4,I2.2,A1,I2.2,A1,I2.2)") tvals(3), "-", tvals(2), "-", tvals(1), &
                    "   ", tvals(5), ":", tvals(6), ":", tvals(7)

  ! > Time value (in hours) since 1st dataset  (CHANGE LATER)
  ivals = (/ syear, smon, sday, 60, shour, smins, ssec, 0/)
  write(st_date, "(A1,I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)") " ", ivals(1), "-", ivals(2), "-", ivals(3), &
                    " ", ivals(5), ":", ivals(6), ":", ivals(7)
  nhours = (itime-1)*dt/3600

  ! netCDF filename specification
  outnc=trim(outfile)//trim(ccdate)//trim(".nc")

  ! Attribute names
  ! > variables

  if ( savehoughwind ) then
    vrname(1) = "H_EIGw"
    vrname(2) = "H_WIGw"
    vrname(3) = "H_BALw"
    lvrname(1) = "Hough Coefficients for Eastward Inertio-Gravity modes - wind only"
    lvrname(2) = "Hough Coefficients for Westward Inertio-Gravity modes - wind only"
    lvrname(3) = "Hough Coefficients for balanced modes - wind only"
  else
    vrname(1) = "EIG"
    vrname(2) = "WIG"
    vrname(3) = "BAL"
    lvrname(1) = "Hough expansion coefficients for Eastward Inertio-Gravity modes"
    lvrname(2) = "Hough expansion coefficients for Westward Inertio-Gravity modes"
    lvrname(3) = "Hough expansion coefficients for balanced modes"
  end if
  uvrname(1) = "~"
  uvrname(2) = "~"
  uvrname(3) = "~"

  ! > coordinates
  dmname(1)  = "n"
  dmname(2)  = "m"
  dmname(3)  = "k"
  dmname(4)  = "Re+Im"
  dmname(5)  = "time"
  ldmname(1) = "meridional mode index"
  ldmname(2) = "vertical mode index"
  ldmname(3) = "zonal wavenumber"
  ldmname(4) = "Real + Imaginary part"
  ldmname(5) = "time"
  udmname(1) = "~"
  udmname(2) = "~"
  udmname(3) = "~"
  udmname(4) = "~"
  udmname(5) = trim("hours since")//trim(st_date)

  ! > global
  if ( savehoughwind ) then
    infoname = "The complex Hough expansion coefficients (wind only)"//&
               ", given separately for EIG, WIG and BAL modes"//&
               ", REAL and IMAG parts"
  else
    infoname = "The complex Hough expansion coefficients"//&
               ", given separately for EIG, WIG and BAL modes"//&
               ", REAL and IMAG parts"
  end if

  !---------------------------------------------------------------------------!
  !                            Write out to netCDF
  !---------------------------------------------------------------------------!

  ! <1> Open NetCDF file
  call check(nf90_create(outnc,NF90_CLOBBER, ncid))

  ! <2> Define dimensions
  call check(nf90_def_dim(ncid, dmname(1), nobs1, n_dimid))
  call check(nf90_def_dim(ncid, dmname(2), nobs2, m_dimid))
  call check(nf90_def_dim(ncid, dmname(3), nobs3, k_dimid))
  call check(nf90_def_dim(ncid, dmname(4), nobs4, c_dimid))
  call check(nf90_def_dim(ncid, dmname(5), NF90_UNLIMITED, t_dimid))

  ! <3> Define coordinate variables
  call check(nf90_def_var(ncid, dmname(1), NF90_DOUBLE, n_dimid, n_varid))
  call check(nf90_def_var(ncid, dmname(2), NF90_DOUBLE, m_dimid, m_varid))
  call check(nf90_def_var(ncid, dmname(3), NF90_DOUBLE, k_dimid, k_varid))
  call check(nf90_def_var(ncid, dmname(5), NF90_DOUBLE, t_dimid, t_varid))
  ! coordinate array
  dimids5D = (/ n_dimid, m_dimid, k_dimid, c_dimid, t_dimid /)

  ! <4> Define variables
  call check(nf90_def_var(ncid, vrname(1), NF90_DOUBLE, dimids5D, EIG_varid))
  call check(nf90_def_var(ncid, vrname(2), NF90_DOUBLE, dimids5D, WIG_varid))
  call check(nf90_def_var(ncid, vrname(3), NF90_DOUBLE, dimids5D, BAL_varid))

  ! <5> Assign attributes
  call check(nf90_put_att(ncid, n_varid, "long_name", ldmname(1)))
  call check(nf90_put_att(ncid, n_varid, "units",     udmname(1)))

  call check(nf90_put_att(ncid, m_varid, "long_name", ldmname(2)))
  call check(nf90_put_att(ncid, m_varid, "units",     udmname(2)))

  call check(nf90_put_att(ncid, k_varid, "long_name", ldmname(3)))
  call check(nf90_put_att(ncid, k_varid, "units",     udmname(3)))

  call check(nf90_put_att(ncid, t_varid, "long_name", ldmname(5)))
  call check(nf90_put_att(ncid, t_varid, "units",     udmname(5)))
  call check(nf90_put_att(ncid, t_varid, "calendar", "standard"))

  call check(nf90_put_att(ncid, EIG_varid, "long_name", lvrname(1)))
  call check(nf90_put_att(ncid, EIG_varid, "units",     uvrname(1)))
  call check(nf90_put_att(ncid, EIG_varid, "grid_type", "gaussian"))

  call check(nf90_put_att(ncid, WIG_varid, "long_name", lvrname(2)))
  call check(nf90_put_att(ncid, WIG_varid, "units",     uvrname(2)))
  call check(nf90_put_att(ncid, WIG_varid, "grid_type", "gaussian"))

  call check(nf90_put_att(ncid, BAL_varid, "long_name", lvrname(3)))
  call check(nf90_put_att(ncid, BAL_varid, "units",     uvrname(3)))
  call check(nf90_put_att(ncid, BAL_varid, "grid_type", "gaussian"))

  ! Global attributes
  if ( savehoughwind ) then
    call write_global_attributes(ncid, "Hough coefficients - wind only", curr_date, infoname)
  else
    call write_global_attributes(ncid, "Hough coefficients", curr_date, infoname)
  end if

  ! <6>  End define mode
  call check(nf90_enddef(ncid))

  ! <7> Write data
  call check(nf90_put_var(ncid, t_varid, nhours))
  call check(nf90_put_var(ncid, m_varid, out_m))
  call check(nf90_put_var(ncid, n_varid, out_n))
  call check(nf90_put_var(ncid, k_varid, out_k))

  call check(nf90_put_var(ncid, EIG_varid, var5D1))
  call check(nf90_put_var(ncid, WIG_varid, var5D2))
  call check(nf90_put_var(ncid, BAL_varid, var5D3))

  ! <8> close NetCDF file
  call check(nf90_close(ncid))
  ! If we got this far, everything worked as expected. Yipee!
  write(*,*) ""
  if ( savehoughwind ) then
    write(*,*) "*** SUCCESSFULLY written Hough coefficients (wind only) to NetCDF file ***"
  else
    write(*,*) "*** SUCCESSFULLY written Hough coefficients to NetCDF file ***"
  end if
  write(*,*) "*** <",trim(outnc),"> ***"
  write(*,*) ""
  end subroutine write_netcdf_proj


  !----------------------------------------------------------------------------


  subroutine write_netcdf_Tmean(outfile, var, nobs1)

  ! Subroutine given arguments
  character(len=*), intent(in) :: outfile
  integer,          intent(in) :: nobs1
  real(r8),         intent(in) :: var(nobs1)

  integer, PARAMETER :: Ndim=2, nobs2=1

  ! NetCDF - IDs and handlers
  integer :: ncid
  integer :: z_dimid, t_dimid
  integer :: z_varid, t_varid
  integer :: Tm_varid

  ! NetCDF attribute and storage info
  character(len=128)                  :: vrname, lvrname, uvrname
  character(len=128), dimension(Ndim) :: dmname, ldmname, udmname
  integer, dimension(Ndim)            :: dimids2D
  character(len=256)                  :: outnc

  ! TIME handling
  integer,dimension(8) :: ivals, tvals
  character(len=256)   :: curr_date, st_date

  ! Output fields
  real(r8) :: out_lev(nobs1)
  real(r8) :: nhours
  real(r8) :: out_var(nobs1,nobs2)

  !---------------------------------------------------------------------------!
  !                               Preparations
  !---------------------------------------------------------------------------!

  ! Define arrays (in sigma levels)
  out_lev = vgrid(nz:1:-1)                    ! Top to bottom
  out_var(:,1) = var

  ! Define time value
  ! > current date and time
  call date_and_time(VALUES=tvals)
  ! > save as current date in specified format
  write(curr_date, "(I2.2,A1,I2.2,A1,I4,A4,I2.2,A1,I2.2,A1,I2.2)") &
    tvals(3), "-", tvals(2), "-", tvals(1), "   ", tvals(5), ":", tvals(6), ":", tvals(7)

  ! > Time value (in hours) since 1st dataset  (CHANGE LATER)
  ivals = (/ syear, smon, sday, 60, shour, smins, ssec, 0/)
  write(st_date, "(A1,I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)") &
    " ", ivals(1), "-", ivals(2), "-", ivals(3), " ", ivals(5), ":", ivals(6), ":", ivals(7)
  nhours = (itime-1)*dt/3600

  ! netCDF filename specification
  outnc = trim(outfile)//trim(".nc")

  ! Attribute names
  ! > variables
  vrname  = "glob"
  lvrname = "Global mean temperature profile"
  uvrname = "K"

  ! > coordinates
  dmname(1)  = "lev"
  dmname(2)  = "time"
  ldmname(1) = "sigma"
  ldmname(2) = "time"
  udmname(1) = "~"
  udmname(2) = trim("hours since")//trim(st_date)

  !---------------------------------------------------------------------------!
  !                            Write out to netCDF
  !---------------------------------------------------------------------------!

  ! <1> Open NetCDF file
  call check(nf90_create(outnc,NF90_CLOBBER, ncid))

  ! <2> Define dimensions
  call check(nf90_def_dim(ncid, dmname(1), nobs1, z_dimid))
  call check(nf90_def_dim(ncid, dmname(2), NF90_UNLIMITED, t_dimid))

  ! <3> Define coordinate variables
  call check(nf90_def_var(ncid, dmname(1), NF90_DOUBLE, z_dimid, z_varid))
  call check(nf90_def_var(ncid, dmname(2), NF90_DOUBLE, t_dimid, t_varid))
  ! coordinate array
  dimids2D = (/ z_dimid, t_dimid /)

  ! <4> Define variables
  call check(nf90_def_var(ncid, vrname, NF90_FLOAT, dimids2D, Tm_varid))

  ! <5> Assign attributes
  call check(nf90_put_att(ncid, z_varid, "long_name", ldmname(1)))
  call check(nf90_put_att(ncid, z_varid, "units", udmname(1)))

  call check(nf90_put_att(ncid, t_varid, "long_name", ldmname(2)))
  call check(nf90_put_att(ncid, t_varid, "units", udmname(2)))
  call check(nf90_put_att(ncid, t_varid, "calendar", "standard"))

  call check(nf90_put_att(ncid, Tm_varid, "long_name", lvrname))
  call check(nf90_put_att(ncid, Tm_varid, "units", uvrname))

  ! Global attributes
  call write_global_attributes(ncid, "Global mean temperature profile", curr_date)

  ! <6>  End define mode
  call check(nf90_enddef(ncid))

  ! <7> Write data
  call check(nf90_put_var(ncid, z_varid, out_lev))
  call check(nf90_put_var(ncid, t_varid, nhours))

  call check(nf90_put_var(ncid, Tm_varid, sngl(out_var)))

  ! <8> close NetCDF file
  call check(nf90_close(ncid))

  ! If we got this far, everything worked as expected. Yipee!
  write(*,*) ""
  write(*,*) "*** SUCCESSFULLY written average T fields to NetCDF file ***"
  write(*,*) "*** <",trim(outnc),"> ***"
  write(*,*) ""

  end subroutine write_netcdf_Tmean


  !----------------------------------------------------------------------------


  subroutine write_netcdf_lnps(outfile, var, nobs1, nobs2)

  ! Subroutine given arguments
  character(*), intent(in) :: outfile
  integer,      intent(in) :: nobs1, nobs2
  real(r8),     intent(in) :: var(nobs1,nobs2)

  integer, PARAMETER :: Ndim=3, nobs3=1

  ! NetCDF - IDs and handlers
  integer :: ncid, I, J
  integer :: t_dimid, x_dimid, y_dimid
  integer :: t_varid, x_varid, y_varid
  integer :: lnps_varid

  ! NetCDF attribute and storage info
  character(len=128)                  :: vrname, lvrname, uvrname
  character(len=128), dimension(Ndim) :: dmname, ldmname, udmname
  integer, dimension(Ndim)            :: dimids3D
  character(len=256)                  :: outnc

  ! TIME handling
  integer,dimension(8) :: ivals, tvals
  character(len=256)   :: curr_date, st_date

  ! Output fields
  real(r8) :: out_lon(nobs1), out_lat(nobs2)
  real(r8) :: nhours
  real(r8) :: out_var(nobs1, nobs2, nobs3)

  !---------------------------------------------------------------------------!
  !                               Preparations
  !---------------------------------------------------------------------------!

  ! Define arrays
  out_lon = (/ (I, I = 0, nx-1) /)*(3600.0_r8/nx)
  out_lat = (/ (J, J = ny-1, 0, -1) /)*(180.0_r8/ny) - 90.0_r8*(1 - 1.0_r8/ny)   ! NP to SP
  !out_lat = (/ (J, J = 0, ny-1) /)*(180./ny) - 90*(1 - 1./ny)       ! SP to NP

  ! > NP to SP
  out_var(:,:,1) = var
  ! > SP to NP
  !do J = 1, ny
  !  out_var(:,J,1) = var(:,ny-J+1)
  !end do

  ! Define time value
  ! > current date and time
  call date_and_time(VALUES=tvals)
  ! > save as current date in specified format
  write(curr_date, "(I2.2,A1,I2.2,A1,I4,A4,I2.2,A1,I2.2,A1,I2.2)") tvals(3), "-", tvals(2), "-", tvals(1), &
                  "   ", tvals(5), ":", tvals(6), ":", tvals(7)

  ! > Time value (in hours) since 1st dataset  (CHANGE LATER)
  ivals = (/ syear, smon, sday, 60, shour, smins, ssec, 0/)
  write(st_date, "(A1,I4,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2,A1,I2.2)") " ", ivals(1), "-", ivals(2), "-", ivals(3), &
                  " ", ivals(5), ":", ivals(6), ":", ivals(7)
  nhours = (itime-1)*dt/3600

  ! netCDF filename specification
  outnc=trim(outfile)//trim(".nc")

  ! Attribute names
  ! > variables
  vrname = "lnsp"
  lvrname = "Logarithm of surface pressure"
  uvrname = "~"

  ! > coordinates
  dmname(1)  = "lon"
  dmname(2)  = "lat"
  dmname(3)  = "time"
  ldmname(1) = "longitude"
  ldmname(2) = "latitude"
  ldmname(3) = "time"
  udmname(1) = "degrees_east"
  udmname(2) = "degrees_north"
  udmname(3) = trim("hours since")//trim(st_date)

  !---------------------------------------------------------------------------!
  !                            Write out to netCDF
  !---------------------------------------------------------------------------!

  ! <1> Open NetCDF file
  call check(nf90_create(outnc,NF90_CLOBBER, ncid))

  ! <2> Define dimensions
  call check(nf90_def_dim(ncid, dmname(1), nobs1, x_dimid))
  call check(nf90_def_dim(ncid, dmname(2), nobs2, y_dimid))
  call check(nf90_def_dim(ncid, dmname(3), NF90_UNLIMITED, t_dimid))

  ! <3> Define coordinate variables
  call check(nf90_def_var(ncid, dmname(1), NF90_DOUBLE, x_dimid, x_varid))
  call check(nf90_def_var(ncid, dmname(2), NF90_DOUBLE, y_dimid, y_varid))
  call check(nf90_def_var(ncid, dmname(3), NF90_DOUBLE, t_dimid, t_varid))
  ! coordinate array
  dimids3D = (/ x_dimid, y_dimid, t_dimid /)

  ! <4> Define variables
  call check(nf90_def_var(ncid, vrname, NF90_FLOAT, dimids3D, lnps_varid))

  ! <5> Assign attributes
  call check(nf90_put_att(ncid, x_varid, "long_name", ldmname(1)))
  call check(nf90_put_att(ncid, x_varid, "units",     udmname(1)))

  call check(nf90_put_att(ncid, y_varid, "long_name", ldmname(2)))
  call check(nf90_put_att(ncid, y_varid, "units",     udmname(2)))

  call check(nf90_put_att(ncid, t_varid, "long_name", ldmname(3)))
  call check(nf90_put_att(ncid, t_varid, "units",     udmname(3)))
  call check(nf90_put_att(ncid, t_varid, "calendar", "standard"))

  call check(nf90_put_att(ncid, lnps_varid, "long_name", lvrname))
  call check(nf90_put_att(ncid, lnps_varid, "units",     uvrname))
  call check(nf90_put_att(ncid, lnps_varid, "grid_type", "gaussian"))

  ! Global attributes
  call write_global_attributes(ncid, "Surface log pressure map", curr_date)

  ! <6>  End define mode
  call check(nf90_enddef(ncid))

  ! <7> Write data
  call check(nf90_put_var(ncid, t_varid, nhours))
  call check(nf90_put_var(ncid, x_varid, out_lon))
  call check(nf90_put_var(ncid, y_varid, out_lat))

  call check(nf90_put_var(ncid, lnps_varid, sngl(out_var)))

  ! <8> close NetCDF file
  call check(nf90_close(ncid))

  ! If we got this far, everything worked as expected. Yipee!
  write(*,*) ""
  write(*,*) "*** SUCCESSFULLY written surface log. pressure field to NetCDF file ***"
  write(*,*) "*** <",trim(outnc),"> ***"
  write(*,*) ""

  end subroutine write_netcdf_lnps


  !----------------------------------------------------------------------------
  !> write_netcdf_vsf


  subroutine write_netcdf_vsf(outfile, mp, num_vmode, var)
  ! routine to write the vertical structure functions to a netCDF file.
  ! Some information is coming from the vsfcalc_cnf namelist read
  ! and stored in module mod_adm

  ! TJH : these are things in the namelist. They might be useful in
  ! the netCDF file (at some point). Because we use the whole mod_atm,
  ! we should be able to write them to the netCDF file.
  !
  ! stab_fname          input file with values of stability on sigma levels
  ! vgrid_fname         input file with value of full sigma levels
  ! vsf_fname           (binary) output file with vertical structure functions
  ! equiheight_fname    output  file with equivalent depths
  ! num_vmode           number of vertical modes for which structure functions
  !                     are stored in the output file
  ! mp                  number of grid points in the input vertical grid
  ! zgrid
  ! top
  ! suft                global surface temperature
  ! hstd                scale height
  ! globalMeanTemp_fname

  ! Subroutine given arguments
  character(*), intent(in) :: outfile
  integer,      intent(in) :: mp
  integer,      intent(in) :: num_vmode
  real(r8),     intent(in) :: var(mp, num_vmode)

! From mod_adm:vsfcalc_init() ?all? things related to vertical structure functions
!  allocate(vsf(mp,num_vmode))
!  allocate(  evht(num_vmode))
!  allocate(        vgrid(mp))
!  allocate( vgrid_weight(mp))
!  allocate(       vsigma(mp))
!  allocate(        vpres(mp))
!  allocate(         stab(mp))

  ! NetCDF - IDs and handlers
  integer :: ncid
  integer ::    mp_dimid,    mp_varid
  integer :: num_vmode_dimid, num_vmode_varid

  integer ::          vsf_varid
  integer ::         evht_varid
  integer ::        vgrid_varid
  integer :: vgrid_weight_varid
  integer ::       vsigma_varid
  integer ::        vpres_varid
  integer ::         stab_varid

  ! NetCDF attribute and storage info
  character(len=256) :: outnc

  integer :: istatus, i

  !---------------------------------------------------------------------------!
  !                               Preparations
  !---------------------------------------------------------------------------!

  ! netCDF filename specification
  outnc = trim(outfile)//trim(".nc")

  !---------------------------------------------------------------------------!
  !                            Write out to netCDF
  !---------------------------------------------------------------------------!

  ! <1> Open NetCDF file
  istatus = nf90_create(outnc, NF90_CLOBBER, ncid)
  call check(istatus, 'write_netcdf_vsf', 'create '//trim(outnc))

  ! <2> Define dimensions
  istatus = nf90_def_dim(ncid, 'mp', mp, mp_dimid)
  call check(istatus, 'write_netcdf_vsf', 'def_dim mp')

  istatus = nf90_def_dim(ncid, 'num_vmode', num_vmode, num_vmode_dimid)
  call check(istatus, 'write_netcdf_vsf', 'def_dim sigma_levels')

  ! <3> Define coordinate variables
  istatus = nf90_def_var(ncid, 'mp',    NF90_INT,    mp_dimid,    mp_varid)
  call check(istatus, 'write_netcdf_vsf', 'def_var mp')

  istatus = nf90_def_var(ncid, 'num_vmode', NF90_INT, num_vmode_dimid, num_vmode_varid)
  call check(istatus, 'write_netcdf_vsf', 'def_var num_vmode')

  ! <4> Define variables
  istatus = nf90_def_var(ncid, 'vsf', NF90_DOUBLE, (/ mp_dimid, num_vmode_dimid /), vsf_varid)
  call check(istatus, 'write_netcdf_vsf', 'def_var vsf')

  istatus = nf90_def_var(ncid, 'evht', NF90_DOUBLE, num_vmode_dimid, evht_varid)
  call check(istatus, 'write_netcdf_vsf', 'def_var evht')

  istatus = nf90_def_var(ncid, 'vgrid', NF90_DOUBLE, mp_dimid, vgrid_varid)
  call check(istatus, 'write_netcdf_vsf', 'def_var vgrid')

  istatus = nf90_def_var(ncid, 'vgrid_weight', NF90_DOUBLE, mp_dimid, vgrid_weight_varid)
  call check(istatus, 'write_netcdf_vsf', 'def_var vgrid_weight')

  istatus = nf90_def_var(ncid, 'vsigma', NF90_DOUBLE, mp_dimid, vsigma_varid)
  call check(istatus, 'write_netcdf_vsf', 'def_var vsigma')

  istatus = nf90_def_var(ncid, 'vpres', NF90_DOUBLE, mp_dimid, vpres_varid)
  call check(istatus, 'write_netcdf_vsf', 'def_var vpres')

  istatus = nf90_def_var(ncid, 'stab', NF90_DOUBLE, mp_dimid, stab_varid)
  call check(istatus, 'write_netcdf_vsf', 'def_var stab')

  ! <5> Assign attributes
  ! call check(nf90_put_att(ncid,    mp_varid, "long_name", 'number of grid points of the input vertical grid'))
  call check(nf90_put_att(ncid, num_vmode_varid, "long_name", 'number of vertical structure functions'))
  call check(istatus, 'write_netcdf_vsf', 'put_att num_vmode long_name')

  istatus = nf90_put_att(ncid, vsf_varid, "long_name", 'vertical structure functions')
  call check(istatus, 'write_netcdf_vsf', 'put_att vsf long_name')

  ! Global attributes
  call write_global_attributes(ncid, "Vertical Structure Functions")

  ! Write the namelist values as global attributes
  istatus = nf90_put_att(ncid, NF90_GLOBAL, "vsfcalc.cnf", "values follow")
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att vsfcalc.cnf')

  istatus = nf90_put_att(ncid, NF90_GLOBAL, "stab_fname", trim(stab_fname))
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att stab_fname')

  istatus = nf90_put_att(ncid, NF90_GLOBAL, "vgrid_fname", trim(vgrid_fname))
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att vgrid_fname')

  istatus = nf90_put_att(ncid, NF90_GLOBAL, "vsf_fname", trim(vsf_fname))
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att vsf_fname')

  istatus = nf90_put_att(ncid, NF90_GLOBAL, "equiheight_fname", trim(equiheight_fname))
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att equiheight_fname')

  istatus = nf90_put_att(ncid, NF90_GLOBAL, "num_vmode", num_vmode)
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att num_vmode')

  istatus = nf90_put_att(ncid, NF90_GLOBAL, "zgrid", trim(zgrid))
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att zgrid')

  istatus = nf90_put_att(ncid, NF90_GLOBAL, "mp", mp)
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att mp')

  istatus = nf90_put_att(ncid, NF90_GLOBAL, "top", top)
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att top')

  istatus = nf90_put_att(ncid, NF90_GLOBAL, "hstd", hstd)
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att hstd')

  if (given_stability) then
     istatus = nf90_put_att(ncid, NF90_GLOBAL, "given_stability", ".true.")
  else
     istatus = nf90_put_att(ncid, NF90_GLOBAL, "given_stability", ".false.")
  endif
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att given_stability')

  if (given_vsf) then
     istatus = nf90_put_att(ncid, NF90_GLOBAL, "given_vsf", ".true.")
  else
     istatus = nf90_put_att(ncid, NF90_GLOBAL, "given_vsf", ".false.")
  endif
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att given_vsf')

  if (ocheck) then
     istatus = nf90_put_att(ncid, NF90_GLOBAL, "ocheck", ".true.")
  else
     istatus = nf90_put_att(ncid, NF90_GLOBAL, "ocheck", ".false.")
  endif
  call check(istatus, 'write_netcdf_vsf', 'GLOBAL put_att ocheck')

  ! <6>  End define mode
  istatus = nf90_enddef(ncid)
  call check(istatus, 'write_netcdf_vsf', 'enddef <'//trim(outnc)//'>')

  ! <7> Write data
  istatus = nf90_put_var(ncid, mp_varid, (/ (i,i=1,mp) /))
  call check(istatus, 'write_netcdf_vsf', 'put_var mp')

  istatus = nf90_put_var(ncid, num_vmode_varid, (/ (i,i=1,num_vmode) /))
  call check(istatus, 'write_netcdf_vsf', 'put_var num_vmode')

  istatus = nf90_put_var(ncid, vsf_varid, var)
  call check(istatus, 'write_netcdf_vsf', 'put_var vsf')

  istatus = nf90_put_var(ncid, evht_varid, evht)
  call check(istatus, 'write_netcdf_evht', 'put_var evht')

  istatus = nf90_put_var(ncid, vgrid_varid, vgrid)
  call check(istatus, 'write_netcdf_vgrid', 'put_var vgrid')

  istatus = nf90_put_var(ncid, vgrid_weight_varid, vgrid_weight)
  call check(istatus, 'write_netcdf_vgrid_weight', 'put_var vgrid_weight')

  istatus = nf90_put_var(ncid, vsigma_varid, vsigma)
  call check(istatus, 'write_netcdf_vsigma', 'put_var vsigma')

  istatus = nf90_put_var(ncid, vpres_varid, vpres)
  call check(istatus, 'write_netcdf_vpres', 'put_var vpres')

  istatus = nf90_put_var(ncid, stab_varid, stab)
  call check(istatus, 'write_netcdf_stab', 'put_var stab')

  ! <8> close NetCDF file
  istatus = nf90_close(ncid)

  ! If we got this far, everything worked as expected.
  write(*,*) ""
  write(*,*) "*** SUCCESSFULLY written vertical structure functions to netCDF file ***"
  write(*,*) "*** <",trim(outnc),"> ***"
  write(*,*) ""

  end subroutine write_netcdf_vsf


  !----------------------------------------------------------------------------
  ! initialize_hough_output
  !> Create the table relating the wavenumber to the netcdf id and filename.
  !> Since each netcdf file will be left open and written to many times,
  !> we need a convenient way to figure out which ncid to write to
  !> given each wavenumber.

  subroutine initialize_hough_output()

     character(len=256) :: filebase
     character(len=256) :: filename

     integer :: zonal_wavenumber, i
     integer :: istatus
     integer :: ncid

     integer :: my_dimid, maxl_dimid, uvz_dimid, num_vmode_dimid
     integer :: my_varid, maxl_varid, uvz_varid, num_vmode_varid
     integer :: varid, y_varid

     real(r8) :: out_lat(my)

     if (allocated(netcdf_table)) return

     out_lat = dasin(ygrid(my:1:-1))*(-90.0_r8)*2.0_r8/pai   ! [South-to-North]

     allocate( netcdf_table(szw:ezw) )

     do zonal_wavenumber = szw,ezw

        call make_idstr( filebase, ofname_bin, 'wn', zonal_wavenumber)

        filename = trim(filebase)//'.nc'

        ! <1> Open NetCDF file
        istatus = nf90_create(filename, NF90_CLOBBER, ncid)
        call check(istatus, 'initialize_hough_output', 'create '//trim(filename))

        netcdf_table(zonal_wavenumber)%filename = filename
        netcdf_table(zonal_wavenumber)%ncid = ncid

        if (debug) write(*,*)'initialize_hough_output: ', &
                   trim(netcdf_table(zonal_wavenumber)%filename), &
                   ' is netCDF ID ', &
                   netcdf_table(zonal_wavenumber)%ncid

        ! <2> Define dimensions
        istatus = nf90_def_dim(ncid, 'num_vmode', num_vmode, num_vmode_dimid)
        call check(istatus, 'initialize_hough_output', 'def_dim num_vmode')

        istatus = nf90_def_dim(ncid, 'uvz', UVZ, uvz_dimid)
        call check(istatus, 'initialize_hough_output', 'def_dim uvz')

        istatus = nf90_def_dim(ncid, 'my', my, my_dimid)
        call check(istatus, 'initialize_hough_output', 'def_dim my')

        istatus = nf90_def_dim(ncid, 'maxl', maxl, maxl_dimid)
        call check(istatus, 'initialize_hough_output', 'def_dim maxl')

        ! <3> Define coordinate variables and attributes

        istatus = nf90_def_var(ncid, 'num_vmode', NF90_INT, num_vmode_dimid, num_vmode_varid)
        call check(istatus, 'initialize_hough_output', 'def_var num_vmode')
        istatus = nf90_put_att(ncid, num_vmode_varid, 'long_name', 'number of vertical modes of structure functions')
        call check(istatus, 'initialize_hough_output', 'put_att num_vmode long_name')

        ! TJH FIXME ... should use U_INDEX, V_INDEX, Z_INDEX to construct strings
        istatus = nf90_def_var(ncid, 'uvz', NF90_INT, uvz_dimid, uvz_varid)
        call check(istatus, 'initialize_hough_output', 'def_var uvz')
        istatus = nf90_put_att(ncid, uvz_varid, 'U', 'uvz = 1')
        call check(istatus, 'initialize_hough_output', 'put_att uvz u')
        istatus = nf90_put_att(ncid, uvz_varid, 'V', 'uvz = 2')
        call check(istatus, 'initialize_hough_output', 'put_att uvz v')
        istatus = nf90_put_att(ncid, uvz_varid, 'Z', 'uvz = 3')
        call check(istatus, 'initialize_hough_output', 'put_att uvz z')

        istatus = nf90_def_var(ncid, 'my', NF90_INT, my_dimid, my_varid)
        call check(istatus, 'initialize_hough_output', 'def_var my')
        istatus = nf90_put_att(ncid, my_varid, 'long_name', 'meridional grid size')
        call check(istatus, 'initialize_hough_output', 'put_att my long_name')

        istatus = nf90_def_var(ncid, 'maxl', NF90_INT, maxl_dimid, maxl_varid)
        call check(istatus, 'initialize_hough_output', 'def_var maxl')
        istatus = nf90_put_att(ncid, maxl_varid, 'long_name', 'number of meridional modes')
        call check(istatus, 'initialize_hough_output', 'put_att maxl long_name')

        ! <4> Define variables and attributes

        istatus = nf90_def_var(ncid, 'EIG', NF90_DOUBLE, &
                  (/ maxl_dimid, my_dimid, uvz_dimid, num_vmode_dimid /), varid)
        call check(istatus, 'initialize_hough_output', 'def_var EIG hough')
        istatus = nf90_put_att(ncid, varid, 'long_name', 'eastward inertio-gravity')
        call check(istatus, 'initialize_hough_output', 'put_att eastward inertio-gravity long_name')

        istatus = nf90_def_var(ncid, 'WIG', NF90_DOUBLE, &
                  (/ maxl_dimid, my_dimid, uvz_dimid, num_vmode_dimid /), varid)
        call check(istatus, 'initialize_hough_output', 'def_var WIG hough')
        istatus = nf90_put_att(ncid, varid, 'long_name', 'westward inertio-gravity')
        call check(istatus, 'initialize_hough_output', 'put_att westward inertio-gravity long_name')

        istatus = nf90_def_var(ncid, 'BAL', NF90_DOUBLE, &
                  (/ maxl_dimid, my_dimid, uvz_dimid, num_vmode_dimid /), varid)
        call check(istatus, 'initialize_hough_output', 'def_var BAL hough')
        istatus = nf90_put_att(ncid, varid, 'long_name', 'balanced')
        call check(istatus, 'initialize_hough_output', 'put_att balanced long_name')

        istatus = nf90_def_var(ncid, 'lat', NF90_DOUBLE, my_dimid, y_varid)
        call check(istatus, 'initialize_hough_output', 'def_var lat hough')
        istatus = nf90_put_att(ncid, y_varid, 'long_name', 'meridional latitudes')
        call check(istatus, 'initialize_hough_output', 'put_att lat long_name')

        ! Global attributes
        call write_global_attributes(ncid, 'Horizontal Structure Functions')

        ! Write houghcalc.cnf values as global attributes
        call check(nf90_put_att(ncid, NF90_GLOBAL, "houghcalc_cnf", "values follow"))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "szw", szw))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "ezw", ezw))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "maxl", maxl))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "my", my))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "ks_mode", trim(ks_mode)))
        call check(nf90_put_att(ncid, NF90_GLOBAL, "freq_fname", trim(freq_fname)))

        ! <6>  End define mode
        istatus = nf90_enddef(ncid)
        call check(istatus, 'initialize_hough_output', 'enddef <'//trim(filename)//'>')

        ! <7>  Fill the coordinate variables.

        istatus = nf90_put_var(ncid, my_varid, (/ (i,i=1,my) /))
        call check(istatus, 'initialize_hough_output', 'put_var my')

        istatus = nf90_put_var(ncid, maxl_varid, (/ (i,i=0,maxl-1) /))
        call check(istatus, 'initialize_hough_output', 'put_var maxl')

        istatus = nf90_put_var(ncid, uvz_varid, (/ (i,i=1,UVZ) /))
        call check(istatus, 'initialize_hough_output', 'put_var uvz')

        istatus = nf90_put_var(ncid, num_vmode_varid, (/ (i,i=1,num_vmode) /))
        call check(istatus, 'initialize_hough_output', 'put_var num_vmode')

        istatus = nf90_put_var(ncid, y_varid, out_lat )
        call check(istatus, 'initialize_hough_output', 'put_var num_vmode')

     enddo

  end subroutine initialize_hough_output


  !----------------------------------------------------------------------------
  !> write_netcdf_hough

  subroutine write_netcdf_hough(wavenumber, mode, datmat)

  ! zgrid
  ! top
  ! suft                global surface temperature
  ! hstd                scale height

  ! Subroutine given arguments
  integer,      intent(in) :: wavenumber
  integer,      intent(in) :: mode
  real(r8),     intent(in) :: datmat(3*maxl,my,UVZ)

  ! NetCDF - IDs and handlers
  integer :: ncid
  integer :: eig_varid
  integer :: wig_varid
  integer :: rot_varid

  ! NetCDF attribute and storage info
  character(len=256) :: outnc

  integer :: istatus, ind1, ind2

  real(r8) :: east(maxl, my, UVZ)
  real(r8) :: west(maxl, my, UVZ)
  real(r8) :: rota(maxl, my, UVZ)

  integer :: nccount(4)
  integer :: ncstart(4)

  !---------------------------------------------------------------------------!
  !                               Preparations
  !---------------------------------------------------------------------------!

  outnc = netcdf_table(wavenumber)%filename
  ncid  = netcdf_table(wavenumber)%ncid

  ! TJH FIXME ... should use [U,V,Z]_INDEX somehow

  ind1 = 1
  ind2 = maxl
  east = datmat(ind1:ind2,:,:)

  ind1 = ind2 + 1
  ind2 = ind1 + maxl -1
  west = datmat(ind1:ind2,:,:)

  ind1 = ind2 + 1
  ind2 = ind1 + maxl -1
  rota = datmat(ind1:ind2,:,:)

  !---------------------------------------------------------------------------!
  !                            Write out to netCDF
  !---------------------------------------------------------------------------!

   istatus = nf90_inq_varid(ncid, 'EIG', eig_varid)
   call check(istatus, 'write_netcdf_hough', 'inq_varid EIG')

   istatus = nf90_inq_varid(ncid, 'WIG', wig_varid)
   call check(istatus, 'write_netcdf_hough', 'inq_varid WIG')

   istatus = nf90_inq_varid(ncid, 'BAL', rot_varid)
   call check(istatus, 'write_netcdf_hough', 'inq_varid BAL')

   ncstart = (/    1,  1,   1, mode /)
   nccount = (/ maxl, my, UVZ,    1 /)

   ! <7> Write data
   istatus = nf90_put_var(ncid, eig_varid, east, start=ncstart, count=nccount)
   call check(istatus, 'write_netcdf_hough', 'put_var east')

   istatus = nf90_put_var(ncid, wig_varid, west, start=ncstart, count=nccount)
   call check(istatus, 'write_netcdf_hough', 'put_var west')

   istatus = nf90_put_var(ncid, rot_varid, rota, start=ncstart, count=nccount)
   call check(istatus, 'write_netcdf_hough', 'put_var rota')

  end subroutine write_netcdf_hough


  !----------------------------------------------------------------------------
  ! finalize_hough_output
  !> Create the table relating the wavenumber to the netcdf id and filename.
  !> Since each netcdf file will be left open and written to many times,
  !> we need a convenient way to figure out which ncid to write to
  !> given each wavenumber.

  subroutine finalize_hough_output()

     integer :: wavenumber, istatus

     do wavenumber = szw,ezw

        istatus = nf90_close(netcdf_table(wavenumber)%ncid)
        call check(istatus, 'finalize_hough_output', &
                   'close '//trim(netcdf_table(wavenumber)%filename))

        ! If we got this far, everything worked as expected.
        write(*,*) ''
        write(*,*) '*** SUCCESSFULLY written hough wavenumber data to netCDF file ***'
        write(*,*) '*** <',trim(netcdf_table(wavenumber)%filename),'> ***'
        write(*,*) ''

     enddo

     if (allocated(netcdf_table)) deallocate(netcdf_table)

  end subroutine finalize_hough_output


  !----------------------------------------------------------------------------

  subroutine write_global_attributes(ncid, titlestring, curr_date, infostring)
  integer,          intent(in) :: ncid
  character(len=*), intent(in) :: titlestring
  character(len=*), intent(in), optional :: curr_date
  character(len=*), intent(in), optional :: infostring

  character(len=256) :: mystring

  write(mystring,*)'MODES: ',trim(titlestring)

  call check(nf90_put_att(ncid, NF90_GLOBAL, "title", trim(mystring)))
  if (present(infostring)) &
  call check(nf90_put_att(ncid, NF90_GLOBAL, "info", infostring))
  if (present(curr_date)) &
  call check(nf90_put_att(ncid, NF90_GLOBAL, "creation_date", curr_date))
  call check(nf90_put_att(ncid, NF90_GLOBAL, "Conventions", "None"))
  call check(nf90_put_att(ncid, NF90_GLOBAL, "support", "http://meteo.fmf.uni-lj.si/MODES/"))
  call check(nf90_put_att(ncid, NF90_GLOBAL, "grant_id", &
                  "European Research Council, Grant Agreement no. 28015"))
  call check(nf90_put_att(ncid, NF90_GLOBAL, "project", "MODES"))
  call check(nf90_put_att(ncid, NF90_GLOBAL, "reference", &
            "Zagar et. al. (2015), Geosci. Model Dev., 8, 1169-1195"))
  call check(nf90_put_att(ncid, NF90_GLOBAL, "reference_doi", &
            "doi:10.5194/gmd-8-1169-2015"))

  end subroutine write_global_attributes


  !----------------------------------------------------------------------------


  subroutine check(istatus, routine, context)
  integer,          intent(in) :: istatus
  character(len=*), intent(in), OPTIONAL :: routine
  character(len=*), intent(in), OPTIONAL :: context

  character(len=512) :: string

  if(istatus == nf90_noerr) return

  !write(*,*)'Correct netCDF return code is ',nf90_noerr
  !write(*,*)'        netCDF return code is ',istatus

  if (present(routine) .and. present(context)) then
     write(string,'(A,'' : '',A,'' : '',A)') trim(routine), &
                trim(context), trim(nf90_strerror(istatus))
  elseif (present(routine)) then
     write(string,'(A,'' : '',A)') trim(routine), trim(nf90_strerror(istatus))
  elseif (present(context)) then
     write(string,'(A,'' : '',A)') trim(context), trim(nf90_strerror(istatus))
  endif

  print *, trim(string)
  stop "ERROR: Stopping."

  end subroutine check

end module mod_write_netcdf
