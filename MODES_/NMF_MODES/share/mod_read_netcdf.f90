module mod_read_netcdf

  implicit none

  public :: read_netcdf_hybrid

  !----------------------------------------------------------------------------
  contains
  !----------------------------------------------------------------------------

  subroutine read_netcdf_hybrid

  use netcdf
  use mod_adm
  use mod_const
  use mod_time
  use mod_write_netcdf, only : check, write_netcdf_lnps

  implicit none

  integer :: inum, idm, ivr, J
  integer :: debug = 0

  ! Input field requirements
  integer                          :: Nf=6
  character(len=128), dimension(4) :: Req_dmnms
  character(len=128), dimension(6) :: Req_vrnms
  character(len=128), dimension(6) :: Req_longvrnms
  logical, dimension(6)            :: inp_check(6) = .false.
  
  ! netCDF info storage
  integer                     :: nFiles, ncid, nDim, nVar, nAtt, unlimDimID 
  character(len=128), allocatable :: dimname(:)
  character(len=128), allocatable :: vrname(:)
  integer, allocatable        :: ldim(:)     ! # elements per dimension, 
  integer, allocatable        :: varNdim(:)  ! # dimensions per variable, 
  integer, allocatable        :: varNatt(:)  ! # attributes per variable
  integer, allocatable        :: vardimIDs(:,:)  ! All variables with their dimension IDs (e.g. lon(lon), u(x,y,z),lnsp(x,y),...)

  ! Input fields
  real(r8)          :: inp_lon(nx), inp_lat(ny), inp_lev(nz)
  real(r8), allocatable   :: inp_values(:,:,:,:)
  real(r8), allocatable   :: inp_values3D(:,:,:)
  real(r8), allocatable   :: inp_values2D(:,:)

  ! Input time
  integer               :: time_offset
  character(len=256)    :: time_units
  !character(len=256)        :: dataTime
  !character(len=256)        :: dataStep

  ! Variables to extract: long names
  Req_longvrnms(1) = "U component of wind"
  Req_longvrnms(2) = "V component of wind"
  Req_longvrnms(3) = "Temperature"
  Req_longvrnms(4) = "Specific humidity"
  Req_longvrnms(5) = "Logarithm of surface pressure"
  Req_longvrnms(6) = "Surface geopotential (orography)"

  ! Variable names as given in netCDF files per orig

  if (orig == 'MERRA') then
     Req_dmnms(1) = "lon"
     Req_dmnms(2) = "lat"
     Req_dmnms(3) = "levels"
     Req_dmnms(4) = "time"
     Req_vrnms(1) = "u"
     Req_vrnms(2) = "v"
     Req_vrnms(3) = "t"
     Req_vrnms(4) = "qv"
     Req_vrnms(5) = "ps"
     Req_vrnms(6) = "PHIS"
  end if

!  if (orig == 'ECMWF') then
!     Req_dmnms(1) = "lon"
!     Req_dmnms(2) = "lat"
!     Req_dmnms(3) = "lev"
!     Req_dmnms(4) = "time"
!     Req_vrnms(1) = "u"
!     Req_vrnms(2) = "v"
!     Req_vrnms(3) = "t"
!     Req_vrnms(4) = "q"
!     Req_vrnms(5) = "lnsp"
!     Req_vrnms(6) = "z"
!  end if

  if (orig == 'ECMWF') then
     Req_dmnms(1) = "lon"
     Req_dmnms(2) = "lat"
     Req_dmnms(3) = "lev"
     Req_dmnms(4) = "time"
     Req_vrnms(1) = "var131"
     Req_vrnms(2) = "var132"
     Req_vrnms(3) = "var130"
     Req_vrnms(4) = "var133"
     Req_vrnms(5) = "var152"
     Req_vrnms(6) = "var129"
  end if

  if (orig == 'CMIP5') then
     Req_dmnms(1) = "lon"
     Req_dmnms(2) = "lat"
     Req_dmnms(3) = "lev"
     Req_dmnms(4) = "time"
     Req_vrnms(1) = "ua"
     Req_vrnms(2) = "va"
     Req_vrnms(3) = "ta"
     Req_vrnms(4) = "hus"
     Req_vrnms(5) = "ps"
     Req_vrnms(6) = "orog"
  end if

  if (orig == 'CCSM') then
     Req_dmnms(1) = "lon"
     Req_dmnms(2) = "lat"
     Req_dmnms(3) = "lev"
     Req_dmnms(4) = "time"
     Req_vrnms(1) = "U"
     Req_vrnms(2) = "V"
     Req_vrnms(3) = "T"
     Req_vrnms(4) = "Q"
     Req_vrnms(5) = "PS"
     Req_vrnms(6) = "PHIS"
  end if

  ! We recommend that each variable carry a "units" attribute.
  !  character (len = *), parameter :: UNITS = "units"
  !  character (len = *), parameter :: PRES_UNITS = "hPa"
  !  character (len = *), parameter :: TEMP_UNITS = "celsius"
  !  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  !  character (len = *), parameter :: LON_UNITS = "degrees_east"

  if (len(trim(ifile_orog)) == 0) Nfiles = numoffile

  ! Include orography nc file in iteration process when specified by user
  if (len(trim(ifile_orog)) /= 0) Nfiles = numoffile+1

    do inum = 1, Nfiles 
      ! <1> Open netCDF file
      write(*,*) 'Reading netCDF file: <', trim(ifile_netcdf(inum)),'>'
      call check( nf90_open(ifile_netcdf(inum), nf90_nowrite, ncid) )

      if (debug > 0) &
         write(*,*) "1. netCDF file: ", trim(ifile_netcdf(inum)), &
                    " with ID ", ncid, " has been opened successfully."

      ! <2> Examine netCDF file (# dimensions, # variables, # attributes)
      call check( nf90_inquire(ncid, nDim, nVar, nAtt, unlimDimID) )

      !write(*,*) "2. # dimensions, # variables and # attributes"
      write(*,*) "Writing netCDF info in iterations= ", inum, nDim, nVar, nAtt
      !write(*,*) "   unlimited dimension ID = ", unlimDimID

      ! - allocate arrays
      allocate(dimname(nDim))
      allocate(ldim(nDim))
      allocate(vrname(nVar))
      allocate(varNdim(nVar))
      allocate(varNatt(nVar))
      allocate(vardimIDs(nVar,nf90_max_var_dims))

      ! - for dimensions
      !write(*,*) "3. Examine dimensions"
      do idm = 1, nDim
        call check( nf90_inquire_dimension(ncid, idm, dimname(idm), ldim(idm)) )
        ! write(*,*) idm, "Dimension name = ", trim(dimname(idm)), ".      Array length: ", ldim(idm)
      end do
      ! result: lon(256), lat(128), lev(60), nhym(60), nhyi(61), lev_2(1), time(1)

      ! - for variables
      !write(*,*) "4. Examine variables"
      do ivr = 1, nVar
        call check( nf90_inquire_variable(ncid, ivr, vrname(ivr), ndims=varNdim(ivr), natts=varNatt(ivr)) )
        ! write(*,*) "Variable name = ", trim(vrname(ivr)), ".       # dim: ", varNdim(ivr), " # att: ", varNatt(ivr)
        call check( nf90_inquire_variable(ncid, ivr, dimids=vardimIDs(ivr,:varNdim(ivr))) )
        ! write(*,*) "Variable dimension IDs: ", vardimIDs(ivr,:varNdim(ivr))
      end do

      ! <3> Overview of variables found in netCDF: 
      write(*,*) "========= netCDF file content overview =========="
      write(*,*) ">> Variables found <<"
      do ivr = 1, nVar
        write(*,*) trim(vrname(ivr)), "(", trim(dimname(vardimIDs(ivr,1))), &
          (",", trim(dimname(vardimIDs(ivr,J))), J=2,varNdim(ivr)), & 
                   ")       with dimensions: (", ldim(vardimIDs(ivr,1)), &
                   (",", ldim(vardimIDs(ivr,J)), J=2,varNdim(ivr)), ')' 
      end do
      write(*,*) ""

      ! <4> Extract grid data
      write(*,*) "============ Start extracting grid data ============"

      do ivr = 1, nVar
        ! - lon
        if (trim(vrname(ivr)) == trim(Req_dmnms(1))) then 
          call check( nf90_get_var(ncid, ivr, inp_lon) )
          write(*,*) "inp_lon:   from", inp_lon(1), "till", inp_lon(ldim(vardimIDs(ivr,1)))
        end if

        ! - lat
        if (trim(vrname(ivr)) == trim(Req_dmnms(2))) then 
          call check( nf90_get_var(ncid, ivr, inp_lat) )
          write(*,*) "inp_lat:   from", inp_lat(1), "till", inp_lat(ldim(vardimIDs(ivr,1)))
        end if

        ! - lev
        if (trim(vrname(ivr)) == trim(Req_dmnms(3))) then 
          call check( nf90_get_var(ncid, ivr, inp_lev) )
          write(*,*) "inp_lev:   from", inp_lev(1), "till", inp_lev(ldim(vardimIDs(ivr,1)))
        end if

        ! - time
        if (trim(vrname(ivr)) == trim(Req_dmnms(4))) then 
          call check( nf90_get_att(ncid, ivr, "units" , time_units) )
          call check( nf90_get_var(ncid, ivr, time_offset) )
          write(*,*) "time = ", trim(time_units), " + ", time_offset
        end if
      end do

   
      ! <5> Extract data per variable
      write(*,*) "============ Start extracting fields ============"
      do ivr = 1, nVar
        ! - u
        if (trim(vrname(ivr)) == trim(Req_vrnms(1))) then
          allocate( inp_values(ldim(vardimIDs(ivr,1)), &
                               ldim(vardimIDs(ivr,2)), &
                               ldim(vardimIDs(ivr,3)), &
                               ldim(vardimIDs(ivr,4))) )  
          call check( nf90_get_var(ncid, ivr, inp_values) )
          u_input = inp_values(:,:,:,1)
          inp_check(1) = .true.
          write(*,*) "Zonal wind field extracted"  
          deallocate(inp_values)
        end if

        ! - v
        if (trim(vrname(ivr)) == trim(Req_vrnms(2))) then
          allocate( inp_values(ldim(vardimIDs(ivr,1)), &
                               ldim(vardimIDs(ivr,2)), &
                               ldim(vardimIDs(ivr,3)), &
                               ldim(vardimIDs(ivr,4))) )  
          call check( nf90_get_var(ncid, ivr, inp_values) )
          v_input = inp_values(:,:,:,1)
          inp_check(2) = .true.
          write(*,*) "Meridional wind field extracted"  
          deallocate(inp_values)
        end if

        ! - t
        if (trim(vrname(ivr)) == trim(Req_vrnms(3))) then
          allocate( inp_values(ldim(vardimIDs(ivr,1)), &
                               ldim(vardimIDs(ivr,2)), &
                               ldim(vardimIDs(ivr,3)), &
                               ldim(vardimIDs(ivr,4))) )  
          call check( nf90_get_var(ncid, ivr, inp_values) )
          t_input = inp_values(:,:,:,1)
          inp_check(3) = .true.
          write(*,*) "Temperature field extracted" 
          deallocate(inp_values) 
        end if

        ! - q
        if (trim(vrname(ivr)) == trim(Req_vrnms(4))) then
          allocate( inp_values(ldim(vardimIDs(ivr,1)), &
                               ldim(vardimIDs(ivr,2)), &
                               ldim(vardimIDs(ivr,3)), &
                               ldim(vardimIDs(ivr,4))) )  
          call check( nf90_get_var(ncid, ivr, inp_values) )
          q_input = inp_values(:,:,:,1)
          inp_check(4) = .true.
          write(*,*) "Specific humidity field extracted" 
          deallocate(inp_values) 
        end if

        ! - lnsp (or sp)
        if (trim(vrname(ivr)) == trim(Req_vrnms(5))) then
          allocate( inp_values3D(ldim(vardimIDs(ivr,1)), &
                                 ldim(vardimIDs(ivr,2)), &
                                 ldim(vardimIDs(ivr,3))) )  
          call check( nf90_get_var(ncid, ivr, inp_values3D) )
          ppp = inp_values3D(:,:,1)
          inp_check(5) = .true.
          write(*,*) "(Log.) surface pressure field extracted"  
          deallocate(inp_values3D)
        end if

        ! - z
        if (trim(vrname(ivr)) == trim(Req_vrnms(6))) then
          allocate( inp_values2D(ldim(vardimIDs(ivr,1)), &
                                 ldim(vardimIDs(ivr,2))) )  
          call check( nf90_get_var(ncid, ivr, inp_values2D) )
          orog = inp_values2D(:,:)
          inp_check(6) = .true.
          write(*,*) "Orography (surface geopotential) extracted"  
          deallocate(inp_values2D)
        end if
      end do

      ! deallocate arrays for this netCDF file
      deallocate(dimname)
      deallocate(ldim)
      deallocate(vrname)
      deallocate(varNdim)
      deallocate(varNatt)
      deallocate(vardimIDs)

      ! close netCDF file
      call check( nf90_close(ncid) )
    end do

    ! <6> Check if all required variables (u,v,t,q,lnsp,z) are extracted from the user-given netCDF files.
    ! If orography file (var = z) is not found yet, then read it from 'ifile_orog' parameter in namelist 
    call check_fields(Nf, inp_check, Req_vrnms, Req_longvrnms)

    ! <7> Check grid data on their order per dimension to make sure:
    ! - Latitude: North pole to south pole
    ! - Altitude: Bottom to top >> also correct for negative orography here
    call check_lat( inp_lat )
    call check_lvl( inp_lev )
    ! - No negative orography, and check surface pressure field too (and save it)
    call check_surface

    ! If we got this far, everything worked as expected. Yipee! 
    write(*,*) "*** SUCCESSFULLY read and processed all netCDF files."
    write(*,*) ""
         
  !----------------------------------------------------------------------------
  contains
  !----------------------------------------------------------------------------

   subroutine check_fields(Nch, arrcheck, varnms, varlng)
   ! This subroutine checks if all fields are imported
   ! - and imports orography from separate file if necessary
   use mod_adm
   implicit none

   integer,            intent(in) :: Nch
   logical,            intent(in) :: arrcheck(Nch)
   character(len=4),   intent(in) :: varnms(Nch)
   character(len=128), intent(in) :: varlng(Nch)

   integer :: ic
  
   do ic=1, Nch

     if (.NOT. arrcheck(ic)) then
       write(*,*) "ERROR: Variable ", trim(varnms(ic)), " (", trim(varlng(ic)), & 
                    " ) not found in the provided netCDF files. Program aborted."
       if (trim(varnms(ic)) == 'z') then 
         write(*,*) "Please supply the orography file in the namelist either"
         write(*,*) "(i)  as a single netCDF file (without date-time string) through 'ifile_orog' or "
         write(*,*) "(ii) for each time step through 'ifile_head' parameter."
       end if
       stop
     end if

   end do
       
   end subroutine check_fields 

   subroutine check_lat(latin)
   ! This subroutine checks and reorders the input data such that it starts from the north pole to south pole
   use mod_adm
   implicit none

   real(r8), intent(inout) :: latin(ny)

   real(r8) :: templat(ny)
   logical  :: revlat
   real(r8) :: x(nx,ny,nz)
   real(r8) :: xs(nx,ny)
   
   revlat = .false. 

   ! lat array should be from NP down to SP
   if (latin(1) < latin(ny)) then
     templat = latin(ny:1:-1) 
     latin   = templat
     revlat  = .true.
     write(*,*) "NOTE: Input data will be latitudally reversed (NP to SP)"
   end if

   !Now, reorder the input data too
   if (revlat) then
     x(:,:,:) = 0.0_r8
     x(:,:,:) = u_input(:,ny:1:-1,:)
     u_input(:,:,:) = x(:,:,:)    

     x(:,:,:) = 0.0_r8
     x(:,:,:) = v_input(:,ny:1:-1,:)
     v_input(:,:,:) = x(:,:,:)

     x(:,:,:) = 0.0_r8
     x(:,:,:) = t_input(:,ny:1:-1,:)
     t_input(:,:,:) = x(:,:,:)

     x(:,:,:) = 0.0_r8
     x(:,:,:) = q_input(:,ny:1:-1,:)
     q_input(:,:,:) = x(:,:,:) 

     xs(:,:) = 0.0_r8
     xs(:,:) = ppp(:,ny:1:-1)
     ppp(:,:) = xs(:,:) 

     xs(:,:) = 0.0_r8
     xs(:,:) = orog(:,ny:1:-1)
     orog(:,:) = xs(:,:) 
   end if
   end subroutine check_lat  

   subroutine check_lvl(levin)
   ! This subroutine checks and reorders the input data such that the level 1 is the top level
   use mod_adm
   implicit none

   real(r8), intent(inout)  :: levin(nz)

   real(r8) :: templev(nz)
   logical  :: revlev
   real(r8) :: x(nx,ny,nz)
    
   revlev = .false. 

   if (levin(1) > levin(nz)) then
     templev = levin(nz:1:-1) 
     levin   = templev
     revlev  = .true.
     write(*,*) "NOTE: Input data will be altitudally reversed (TOP to BOTTOM)"
   end if

       !Now, reorder the input data too
   if (revlev) then
     x(:,:,:) = 0.0_r8
     x(:,:,:) = u_input(:,:,nz:1:-1)
     u_input(:,:,:) = x(:,:,:)    

     x(:,:,:) = 0.0_r8
     x(:,:,:) = v_input(:,:,nz:1:-1)
     v_input(:,:,:) = x(:,:,:)

     x(:,:,:) = 0.0_r8
     x(:,:,:) = t_input(:,:,nz:1:-1)
     t_input(:,:,:) = x(:,:,:)

     x(:,:,:) = 0.0_r8
     x(:,:,:) = q_input(:,:,nz:1:-1)
     q_input(:,:,:) = x(:,:,:) 
   end if     
   end subroutine check_lvl 

   subroutine check_surface
   ! This subroutine checks surface values given for pressure and orography
   use mod_adm
   implicit none

!   integer ::  iy, ix
   real(r8) :: flilla
   real(r8) :: xs(nx,ny)
   character(len=256) :: fname, tmp_fname

   ! Remove negative "orography" here 
   flilla = 1.0d0-15

   where (orog < 0.0_r8) orog = flilla

   ! Surface geopotential (in m^2s^-2) or surface orography (in m) given
   ! If given in m above sea level, then multiply with g=9.81
   if (maxval(maxval(orog,1)) < 10000.0) then
     write(*,*) "Max value of orography, likely given by user as 'height above sea'.", maxval(maxval(orog,1))
     write(*,*) "Too low value for surface geopotential >> multiply by 9.81 to get correct values. "
       xs(:,:) = orog(:,:)*g
     orog(:,:) = xs(:,:)
     write(*,*) "New max value of surface geopotential: ", maxval(maxval(orog,1))
   end if

   ! Surface pressure values given in logarithmic form or not?
   ! If given in log(sp) form then averaged field value wouldn't exceed 50
   ! If it does exceed, derive lnsp = dlog(sp)
   if (sum(ppp)/(nx*ny) > 50.0_r8) then
     write(*,*) "Average surface pressure: ", sum(ppp)/(nx*ny)
     write(*,*) "Too high values >> ln(Ps) values are expected."
      xs(:,:) = dlog(ppp(:,:))
     ppp(:,:) = xs(:,:)
     write(*,*) "New values of log(surface pressure): ", sum(ppp)/(nx*ny)
   end if

   ! save surface pressure for later use during the inversion
   if(saveps) then
     fname='Ps_'
     tmp_fname=trim(ps_fname)//trim(cdate)
     call output_2d(trim(tmp_fname), ppp, nx, ny, 1, .true., 8) 
     if(savenc) call write_netcdf_lnps(trim(tmp_fname), ppp, nx, ny)
   end if

   end subroutine check_surface 

  end subroutine read_netcdf_hybrid
end module mod_read_netcdf
