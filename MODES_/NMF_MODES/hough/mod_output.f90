module mod_output

  use mod_const, only : r8
  use mod_adm,   only : mk, maxl, szw, ezw, bin_combine, make_idstr, &
                        ofname_gmt, ofname_bin, num_vmode
 
  use mod_write_netcdf, only : initialize_hough_output, &
                               write_netcdf_hough, &
                               finalize_hough_output 

  implicit none
  private

  public :: output_for_gmt
  public :: output_for_binary
  public :: output_for_netcdf

contains

  !----------------------------------------------------------------------------

  subroutine output_for_gmt( ls, m, hftbl, iparj )

! Koji's need    use mod_gaussg, only : deg

    integer,  intent(in) :: ls
    integer,  intent(in) :: m
    real(r8), intent(in) :: hftbl(mk,3*maxl, 4)
    integer,  intent(in) :: iparj(   3*maxl   )

    integer  :: my, i, iy
    real(r8) :: sym
    character(len=128) :: hough_fname
    real(r8), allocatable :: h0(:,:,:)

    my = 2*mk
    allocate( h0( my, 3*maxl, 3 ) )

    call compute_h0(maxl, mk, iparj, hftbl, h0)

    if( trim(bin_combine) == 'zonal' ) then
!      call make_idstr( hough_fname, trim(ofname_gmt), 'wn', ls)
      write(*,*) 'Writing: ', ls
!      open(100+ls, file=trim(hough_fname), position='append', status='unknown')
      do i = 1, maxl * 3
      do iy = 1, my
        write(100+ls,'(3f20.12,1X)') h0(iy, i, 1),h0(iy, i, 2),h0(iy, i, 3)
      end do
      end do
!      close(100+ls)
    end if

    deallocate(h0)

! Koji's code
!    do i = 1, maxl*3
!      if( iparj(i) == 1 ) sym =  1.0_r8
!      if( iparj(i) == 2 ) sym = -1.0_r8
!      do iy = 1, mk
!        write(101,'(3i4,2f20.12)') ls, m, i, deg(iy), hftbl(iy, i, 1)
!        write(102,'(3i4,2f20.12)') ls, m, i, deg(iy), hftbl(iy, i, 2)
!        write(103,'(3i4,2f20.12)') ls, m, i, deg(iy), hftbl(iy, i, 3)
!      end do
!      do iy = mk, 1, -1
!        write(101,'(3i4,2f20.12)') ls, m, i, -deg(iy), -sym*hftbl(iy, i, 1)
!        write(102,'(3i4,2f20.12)') ls, m, i, -deg(iy),  sym*hftbl(iy, i, 2)
!        write(103,'(3i4,2f20.12)') ls, m, i, -deg(iy), -sym*hftbl(iy, i, 3)
!      end do
!    end do

  end subroutine output_for_gmt

  !----------------------------------------------------------------------------

  subroutine output_for_binary( ls, m, hftbl, iparj )

    integer,  intent(in) :: ls
    integer,  intent(in) :: m
    real(r8), intent(in) :: hftbl(mk,3*maxl, 4)
    integer,  intent(in) :: iparj(   3*maxl   )

    character(len=128) :: hough_fname
    real(r8), allocatable :: h0(:,:,:)
    integer :: my     ! size of grid in meridional direction

    ! mk is the number of meridional points of one hemisphere.
    ! magic 3 is U,V,ZG
    ! This routine gets called num_vertical_modes times.

    my = 2*mk
    allocate( h0( my, 3*maxl, 3 ) )

    call compute_h0(maxl, mk, iparj, hftbl, h0)

    if( trim(bin_combine) == 'zonal' ) then
      call make_idstr( hough_fname, trim(ofname_bin), 'wn', ls)
      write(*,*) ls, trim(hough_fname),my,maxl,m
      open(11, file=trim(hough_fname), form='unformatted', &
        access='direct', recl=8*my*maxl*3*3, status='unknown')
        !access='direct', recl=8*mk*2*maxl*3*3, status='unknown')
      write(11,rec=m) h0
      close(11)
    end if

    deallocate(h0)

    ! TJH FIXME what happens if bin_combine /= zonal

  end subroutine output_for_binary

  !----------------------------------------------------------------------------

  subroutine output_for_netcdf( zonal_wavenumber, mode_index, hftbl, iparj )

    integer,  intent(in) :: zonal_wavenumber ! [szw:ezw]
    integer,  intent(in) :: mode_index   ! index of current vertical mode
    real(r8), intent(in) :: hftbl(mk,3*maxl, 4)
    integer,  intent(in) :: iparj(   3*maxl   )

    ! maxl is the number of meridional modes
    !      [eig, wig, rot] == 3 modes

    character(len=128) :: hough_fname
    real(r8), allocatable :: h0(:,:,:), hnew(:,:,:)
    integer  :: i
    integer  :: meridional_grid_size, iy
    integer  :: num_meridional_modes
    integer  :: maxcount

    real(r8) :: sym

    integer, parameter :: UVZ = 3
    integer, parameter :: EWR = 3

    integer, save :: call_count = 0

    meridional_grid_size = 2*mk
    num_meridional_modes = EWR*maxl
    maxcount = (ezw - szw + 1) * num_vmode

    if ( call_count == 0 ) call initialize_hough_output()
    call_count = call_count + 1

    ! write(*,*)'TJH: szw, zonal_wavenum, ezw, vert_mode_index, call, max', szw, zonal_wavenumber, ezw, mode_index, call_count, maxcount

    allocate(   h0( meridional_grid_size, num_meridional_modes, UVZ ) )
    allocate( hnew( num_meridional_modes, meridional_grid_size, UVZ ) )  ! permuted dimensions

    call compute_h0(maxl, mk, iparj, hftbl, h0)

    call permute(h0, hnew)

    ! do the permutation

    call write_netcdf_hough(zonal_wavenumber, mode_index, hnew)

    deallocate(h0)
    deallocate(hnew)

    ! TJH FIXME what happens if bin_combine /= zonal

    if ( call_count == maxcount ) call finalize_hough_output()

  end subroutine output_for_netcdf

  !----------------------------------------------------------------------------

  subroutine compute_h0(maxl, mk, iparj, hftbl, h0)

    integer,  intent(in)    :: maxl
    integer,  intent(in)    :: mk
    integer,  intent(in)    :: iparj(maxl*3)
    real(r8), intent(in)    :: hftbl(mk  , maxl*3, 3)
    real(r8), intent(inout) ::    h0(mk*2, maxl*3, 3)

    integer  :: my, i, iy
    real(r8) :: sym

    integer, parameter :: EIG_INDEX = 1
    integer, parameter :: WIG_INDEX = 2
    integer, parameter :: ROT_INDEX = 3

    my = 2*mk
    do i = 1, maxl * 3
      if( iparj(i) == 1 ) sym =  1.0_r8
      if( iparj(i) == 2 ) sym = -1.0_r8
      do iy = 1, mk
        h0(     iy, i, EIG_INDEX) = -sym*hftbl(iy, i, EIG_INDEX)
        h0(     iy, i, WIG_INDEX) =  sym*hftbl(iy, i, WIG_INDEX)
        h0(     iy, i, ROT_INDEX) = -sym*hftbl(iy, i, ROT_INDEX)
        h0(my-iy+1, i, EIG_INDEX) =      hftbl(iy, i, EIG_INDEX)
        h0(my-iy+1, i, WIG_INDEX) =      hftbl(iy, i, WIG_INDEX)
        h0(my-iy+1, i, ROT_INDEX) =      hftbl(iy, i, ROT_INDEX)
      end do
    end do

  end subroutine compute_h0

  !----------------------------------------------------------------------------

  subroutine permute(x, y)

  !  allocate(   h0( meridional_grid_size, num_meridional_modes, UVZ ) )
  !  allocate( hnew( num_meridional_modes, meridional_grid_size, UVZ ) )  ! permuted dimensions
    real(r8), intent(in)  :: x(:,:,:)
    real(r8), intent(out) :: y(:,:,:)

    integer  :: i1, i2, i3

    do i3 = 1, size(x,3)
    do i2 = 1, size(x,2)
    do i1 = 1, size(x,1)
        y(i2, i1, i3) = x(i1, i2, i3)
    end do
    end do
    end do

  end subroutine permute

end module mod_output
