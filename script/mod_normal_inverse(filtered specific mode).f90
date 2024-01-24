module mod_normal_inverse
  implicit none
  private
  public :: NMFinv_output
  public :: read_3dNMFcoefdata
  public :: filter_3dNMFcoefdata
  public :: sigma2hybrid
  !
contains
  !----------------------------------------------------------------------------
  subroutine read_3DNMFcoefdata(itime)
    use mod_adm
    use mod_time, only: &
        cdate
    use mod_const, only :&
        g, rd
    implicit none
    integer :: itime
    character(len=256) :: tmp_fname

    if(separate_time) then
      tmp_fname=trim(coef3DNMF_fname)//trim(cdate)
      open(1,file=trim(tmp_fname),form='unformatted',access='direct', &
               recl=16*(num_zw+1)*maxl*3*num_vmode)
      read(1,rec=1)  w0
    else
      open(1,file=trim(coef3DNMF_fname),form='unformatted',access='direct', &
               recl=16*(num_zw+1)*maxl*3*num_vmode)
      read(1,rec=itime)  w0
    end if

  end subroutine read_3DNMFcoefdata
!----------------------------------------------------------------------------
  subroutine filter_3DNMFcoefdata
!
! all vertical modes specified are filtered within specified range 
! all meridional modes are filtered separately for every mode type
! for all vertical and zonal modes
! zonal modes are filtered for the whole specified range  
    use mod_adm
    use mod_const
    implicit none
    integer :: n, k, m
    double complex,   allocatable :: w0_buffer(:,:,:)
    allocate(w0_buffer(maxl*3,num_vmode,0:num_zw))
    ! Copy data from w0 to w0_buffer and set w0 to 0 
    do k = 0, num_zw
        do m = 1, num_vmode
            do n = 1, maxl*3
                w0_buffer(n, m, k) = w0(n, m, k)
                w0(n,m,k) = (0.0d0,0.0d0)
            end do
        end do
    end do
    
    ! save the user specific modes

        ! Kelvin
    w0(1,2,1) = w0_buffer(1,2,1)

    w0(maxl+1,:,0)=w0_buffer(maxl+1,:,0)


!    ! Top 100 EIG
!     w0(1,10,1)=w0_buffer(1,10,1)
!     w0(1,13,1)=w0_buffer(1,13,1)
!     w0(1,7,1)=w0_buffer(1,7,1)
!     w0(1,9,1)=w0_buffer(1,9,1)
!     w0(1,8,1)=w0_buffer(1,8,1)
!     w0(6,3,0)=w0_buffer(6,3,0)
!     w0(1,18,3)=w0_buffer(1,18,3)
!     w0(1,19,3)=w0_buffer(1,19,3)
!     w0(2,9,0)=w0_buffer(2,9,0)
!     w0(1,20,3)=w0_buffer(1,20,3)
!     w0(1,12,1)=w0_buffer(1,12,1)
!     w0(1,20,1)=w0_buffer(1,20,1)
!     w0(1,20,2)=w0_buffer(1,20,2)
!     w0(1,18,1)=w0_buffer(1,18,1)
!     w0(1,9,2)=w0_buffer(1,9,2)
!     w0(2,11,0)=w0_buffer(2,11,0)
!     w0(2,14,0)=w0_buffer(2,14,0)
!     w0(1,14,1)=w0_buffer(1,14,1)
!     w0(4,3,0)=w0_buffer(4,3,0)
!     w0(1,11,2)=w0_buffer(1,11,2)
!     w0(1,6,1)=w0_buffer(1,6,1)
!     w0(1,19,2)=w0_buffer(1,19,2)
!     w0(1,14,2)=w0_buffer(1,14,2)
!     w0(2,12,0)=w0_buffer(2,12,0)
!     w0(1,11,1)=w0_buffer(1,11,1)
!     w0(1,17,3)=w0_buffer(1,17,3)
!     w0(1,19,1)=w0_buffer(1,19,1)
!     w0(1,8,2)=w0_buffer(1,8,2)
!     w0(1,18,2)=w0_buffer(1,18,2)
!     w0(1,12,2)=w0_buffer(1,12,2)
!     w0(3,11,0)=w0_buffer(3,11,0)
!     w0(1,17,1)=w0_buffer(1,17,1)
!     w0(1,16,3)=w0_buffer(1,16,3)
!     w0(3,8,0)=w0_buffer(3,8,0)
!     w0(3,9,0)=w0_buffer(3,9,0)
!     w0(3,14,0)=w0_buffer(3,14,0)
!     w0(2,13,0)=w0_buffer(2,13,0)
!     w0(1,15,1)=w0_buffer(1,15,1)
!     w0(7,3,0)=w0_buffer(7,3,0)
!     w0(3,2,0)=w0_buffer(3,2,0)
!     w0(3,19,0)=w0_buffer(3,19,0)
!     w0(11,8,0)=w0_buffer(11,8,0)
!     w0(3,12,0)=w0_buffer(3,12,0)
!     w0(1,10,4)=w0_buffer(1,10,4)
!     w0(1,13,2)=w0_buffer(1,13,2)
!     w0(1,16,2)=w0_buffer(1,16,2)
!     w0(1,15,3)=w0_buffer(1,15,3)
!     w0(1,4,1)=w0_buffer(1,4,1)
!     w0(9,6,0)=w0_buffer(9,6,0)
!     w0(12,8,0)=w0_buffer(12,8,0)
!     w0(9,5,0)=w0_buffer(9,5,0)
!     w0(3,16,0)=w0_buffer(3,16,0)
!     w0(9,3,0)=w0_buffer(9,3,0)
!     w0(2,16,0)=w0_buffer(2,16,0)
!     w0(3,17,0)=w0_buffer(3,17,0)
!     w0(10,8,0)=w0_buffer(10,8,0)
!     w0(1,15,2)=w0_buffer(1,15,2)
!     w0(1,15,4)=w0_buffer(1,15,4)
!     w0(3,20,0)=w0_buffer(3,20,0)
!     w0(3,1,0)=w0_buffer(3,1,0)
!     w0(2,2,0)=w0_buffer(2,2,0)
!     w0(1,9,4)=w0_buffer(1,9,4)
!     w0(2,8,0)=w0_buffer(2,8,0)
!     w0(1,16,1)=w0_buffer(1,16,1)
!     w0(3,13,0)=w0_buffer(3,13,0)
!     w0(2,10,0)=w0_buffer(2,10,0)
!     w0(10,6,0)=w0_buffer(10,6,0)
!     w0(2,15,0)=w0_buffer(2,15,0)
!     w0(6,3,1)=w0_buffer(6,3,1)
!     w0(12,7,0)=w0_buffer(12,7,0)
!     w0(7,5,0)=w0_buffer(7,5,0)
!     w0(1,14,5)=w0_buffer(1,14,5)
!     w0(2,17,0)=w0_buffer(2,17,0)
!     w0(1,17,2)=w0_buffer(1,17,2)
!     w0(1,16,4)=w0_buffer(1,16,4)
!     w0(9,3,1)=w0_buffer(9,3,1)
!     w0(8,6,0)=w0_buffer(8,6,0)
!     w0(2,19,0)=w0_buffer(2,19,0)
!     w0(4,6,0)=w0_buffer(4,6,0)
!     w0(1,10,3)=w0_buffer(1,10,3)
!     w0(3,10,0)=w0_buffer(3,10,0)
!     w0(9,8,0)=w0_buffer(9,8,0)
!     w0(11,7,0)=w0_buffer(11,7,0)
!     w0(1,11,3)=w0_buffer(1,11,3)
!     w0(7,6,0)=w0_buffer(7,6,0)
!     w0(11,6,0)=w0_buffer(11,6,0)
!     w0(1,6,2)=w0_buffer(1,6,2)
!     w0(3,4,0)=w0_buffer(3,4,0)
!     w0(2,20,0)=w0_buffer(2,20,0)
!     w0(9,5,1)=w0_buffer(9,5,1)
!     w0(13,10,0)=w0_buffer(13,10,0)
!     w0(3,10,2)=w0_buffer(3,10,2)
!     w0(1,13,4)=w0_buffer(1,13,4)
!     w0(13,8,1)=w0_buffer(13,8,1)
!     w0(12,7,1)=w0_buffer(12,7,1)
!     w0(2,18,0)=w0_buffer(2,18,0)
!     w0(6,6,0)=w0_buffer(6,6,0)
!     w0(2,10,1)=w0_buffer(2,10,1)
!     w0(2,6,0)=w0_buffer(2,6,0)
!     w0(13,8,0)=w0_buffer(13,8,0)


!     w0(31,20,0)=w0_buffer(31,20,0)
!     w0(31,19,0)=w0_buffer(31,19,0)
!     w0(31,18,0)=w0_buffer(31,18,0)
!     w0(31,15,0)=w0_buffer(31,15,0)
!     w0(31,16,0)=w0_buffer(31,16,0)
!     w0(31,17,0)=w0_buffer(31,17,0)
!     w0(31,14,0)=w0_buffer(31,14,0)
!     w0(31,13,0)=w0_buffer(31,13,0)
!     w0(31,12,0)=w0_buffer(31,12,0)
!     w0(31,11,0)=w0_buffer(31,11,0)
!     w0(31,10,0)=w0_buffer(31,10,0)
!     w0(31,9,0)=w0_buffer(31,9,0)
!     w0(31,1,0)=w0_buffer(31,1,0)
!     w0(31,8,0)=w0_buffer(31,8,0)
!     w0(31,6,0)=w0_buffer(31,6,0)
!     w0(31,7,0)=w0_buffer(31,7,0)
!     w0(31,4,0)=w0_buffer(31,4,0)
!     w0(31,3,0)=w0_buffer(31,3,0)
!     w0(31,2,0)=w0_buffer(31,2,0)
!     w0(31,5,0)=w0_buffer(31,5,0)
!     w0(36,3,0)=w0_buffer(36,3,0)
!     w0(32,9,0)=w0_buffer(32,9,0)
!     w0(32,11,0)=w0_buffer(32,11,0)
!     w0(32,14,0)=w0_buffer(32,14,0)
!     w0(34,3,0)=w0_buffer(34,3,0)
!     w0(32,12,0)=w0_buffer(32,12,0)
!     w0(33,11,0)=w0_buffer(33,11,0)
!     w0(33,8,0)=w0_buffer(33,8,0)
!     w0(33,9,0)=w0_buffer(33,9,0)
!     w0(33,14,0)=w0_buffer(33,14,0)
!     w0(32,13,0)=w0_buffer(32,13,0)
!     w0(37,3,0)=w0_buffer(37,3,0)
!     w0(33,2,0)=w0_buffer(33,2,0)
!     w0(33,19,0)=w0_buffer(33,19,0)
!     w0(41,8,0)=w0_buffer(41,8,0)
!     w0(33,12,0)=w0_buffer(33,12,0)
!     w0(39,6,0)=w0_buffer(39,6,0)
!     w0(42,8,0)=w0_buffer(42,8,0)
!     w0(39,5,0)=w0_buffer(39,5,0)
!     w0(33,16,0)=w0_buffer(33,16,0)
!     w0(39,3,0)=w0_buffer(39,3,0)
!     w0(32,16,0)=w0_buffer(32,16,0)
!     w0(33,17,0)=w0_buffer(33,17,0)
!     w0(40,8,0)=w0_buffer(40,8,0)
!     w0(34,7,2)=w0_buffer(34,7,2)
!     w0(33,20,0)=w0_buffer(33,20,0)
!     w0(33,1,0)=w0_buffer(33,1,0)
!     w0(32,2,0)=w0_buffer(32,2,0)
!     w0(39,3,1)=w0_buffer(39,3,1)
!     w0(32,8,0)=w0_buffer(32,8,0)
!     w0(33,13,0)=w0_buffer(33,13,0)
!     w0(32,10,0)=w0_buffer(32,10,0)
!     w0(40,6,0)=w0_buffer(40,6,0)
!     w0(32,15,0)=w0_buffer(32,15,0)
!     w0(42,7,0)=w0_buffer(42,7,0)
!     w0(31,8,3)=w0_buffer(31,8,3)
!     w0(37,5,0)=w0_buffer(37,5,0)
!     w0(32,17,0)=w0_buffer(32,17,0)
!     w0(36,3,1)=w0_buffer(36,3,1)
!     w0(38,6,0)=w0_buffer(38,6,0)
!     w0(32,19,0)=w0_buffer(32,19,0)
!     w0(42,8,1)=w0_buffer(42,8,1)
!     w0(34,6,0)=w0_buffer(34,6,0)
!     w0(33,10,0)=w0_buffer(33,10,0)
!     w0(39,8,0)=w0_buffer(39,8,0)
!     w0(41,7,0)=w0_buffer(41,7,0)
!     w0(37,6,0)=w0_buffer(37,6,0)
!     w0(41,6,0)=w0_buffer(41,6,0)
!     w0(42,7,1)=w0_buffer(42,7,1)
!     w0(33,4,0)=w0_buffer(33,4,0)
!     w0(32,20,0)=w0_buffer(32,20,0)
!     w0(43,10,0)=w0_buffer(43,10,0)
!     w0(32,3,2)=w0_buffer(32,3,2)
!     w0(32,18,0)=w0_buffer(32,18,0)
!     w0(36,6,0)=w0_buffer(36,6,0)
!     w0(35,7,2)=w0_buffer(35,7,2)
!     w0(39,5,1)=w0_buffer(39,5,1)
!     w0(32,6,0)=w0_buffer(32,6,0)
!     w0(43,8,0)=w0_buffer(43,8,0)
!     w0(39,7,0)=w0_buffer(39,7,0)
!     w0(34,8,0)=w0_buffer(34,8,0)
!     w0(38,3,0)=w0_buffer(38,3,0)
!     w0(34,7,1)=w0_buffer(34,7,1)
!     w0(33,18,0)=w0_buffer(33,18,0)
!     w0(31,8,4)=w0_buffer(31,8,4)
!     w0(44,10,0)=w0_buffer(44,10,0)
!     w0(32,7,0)=w0_buffer(32,7,0)
!     w0(34,20,0)=w0_buffer(34,20,0)
!     w0(34,19,0)=w0_buffer(34,19,0)
!     w0(33,3,2)=w0_buffer(33,3,2)
!     w0(31,9,3)=w0_buffer(31,9,3)
!     w0(34,16,0)=w0_buffer(34,16,0)
!     w0(38,5,0)=w0_buffer(38,5,0)
!     w0(44,9,0)=w0_buffer(44,9,0)
!     w0(41,9,0)=w0_buffer(41,9,0)
!     w0(43,8,1)=w0_buffer(43,8,1)
!     w0(42,5,1)=w0_buffer(42,5,1)
!     w0(40,5,0)=w0_buffer(40,5,0)
!     w0(31,10,1)=w0_buffer(31,10,1)
!     w0(44,8,0)=w0_buffer(44,8,0)
!     w0(42,3,1)=w0_buffer(42,3,1)
!     w0(42,5,0)=w0_buffer(42,5,0)
!     w0(41,7,1)=w0_buffer(41,7,1)
!     w0(35,2,0)=w0_buffer(35,2,0)
!     w0(31,20,2)=w0_buffer(31,20,2)
!     w0(40,3,0)=w0_buffer(40,3,0)
!     w0(37,3,1)=w0_buffer(37,3,1)
!     w0(44,8,1)=w0_buffer(44,8,1)
!     w0(34,7,4)=w0_buffer(34,7,4)
!     w0(33,7,2)=w0_buffer(33,7,2)
!     w0(41,6,1)=w0_buffer(41,6,1)
!     w0(35,16,0)=w0_buffer(35,16,0)
!     w0(35,14,0)=w0_buffer(35,14,0)
!     w0(35,6,0)=w0_buffer(35,6,0)
!     w0(35,13,0)=w0_buffer(35,13,0)
!     w0(33,6,0)=w0_buffer(33,6,0)
!     w0(36,8,0)=w0_buffer(36,8,0)
!     w0(42,9,0)=w0_buffer(42,9,0)
!     w0(44,9,1)=w0_buffer(44,9,1)
!     w0(31,19,2)=w0_buffer(31,19,2)

!      w0(62,13,0)=w0_buffer(62,13,0)
!      w0(62,11,0)=w0_buffer(62,11,0)
!      w0(62,12,0)=w0_buffer(62,12,0)
!      w0(62,10,0)=w0_buffer(62,10,0)
!      w0(62,14,0)=w0_buffer(62,14,0)
!      w0(62,9,0)=w0_buffer(62,9,0)
!      w0(62,8,0)=w0_buffer(62,8,0)
!      w0(62,15,0)=w0_buffer(62,15,0)
!      w0(62,7,0)=w0_buffer(62,7,0)
!      w0(62,16,0)=w0_buffer(62,16,0)
!      w0(63,2,0)=w0_buffer(63,2,0)
!      w0(63,14,0)=w0_buffer(63,14,0)
!      w0(63,12,0)=w0_buffer(63,12,0)
!      w0(62,17,0)=w0_buffer(62,17,0)
!      w0(63,16,0)=w0_buffer(63,16,0)
!      w0(63,10,0)=w0_buffer(63,10,0)
!      w0(63,7,0)=w0_buffer(63,7,0)
!      w0(63,15,0)=w0_buffer(63,15,0)
!      w0(62,6,0)=w0_buffer(62,6,0)
!      w0(63,3,0)=w0_buffer(63,3,0)
!      w0(63,17,0)=w0_buffer(63,17,0)
!      w0(63,13,0)=w0_buffer(63,13,0)
!      w0(66,8,0)=w0_buffer(66,8,0)
!      w0(63,20,0)=w0_buffer(63,20,0)
!      w0(63,18,0)=w0_buffer(63,18,0)
!      w0(66,9,0)=w0_buffer(66,9,0)
!      w0(63,19,0)=w0_buffer(63,19,0)
!      w0(63,1,0)=w0_buffer(63,1,0)
!      w0(65,13,0)=w0_buffer(65,13,0)
!      w0(62,18,0)=w0_buffer(62,18,0)
!      w0(65,12,0)=w0_buffer(65,12,0)
!      w0(64,1,0)=w0_buffer(64,1,0)
!      w0(65,14,0)=w0_buffer(65,14,0)
!      w0(63,8,0)=w0_buffer(63,8,0)
!      w0(68,9,0)=w0_buffer(68,9,0)
!      w0(65,11,0)=w0_buffer(65,11,0)
!      w0(63,5,0)=w0_buffer(63,5,0)
!      w0(62,2,0)=w0_buffer(62,2,0)
!      w0(63,4,0)=w0_buffer(63,4,0)
!      w0(65,10,0)=w0_buffer(65,10,0)
!      w0(63,9,0)=w0_buffer(63,9,0)
!      w0(65,15,0)=w0_buffer(65,15,0)
!      w0(64,3,0)=w0_buffer(64,3,0)
!      w0(68,8,0)=w0_buffer(68,8,0)
!      w0(63,11,0)=w0_buffer(63,11,0)
!      w0(62,19,0)=w0_buffer(62,19,0)
!      w0(65,3,0)=w0_buffer(65,3,0)
!      w0(62,3,0)=w0_buffer(62,3,0)
!      w0(65,8,0)=w0_buffer(65,8,0)
!      w0(64,10,0)=w0_buffer(64,10,0)
!      w0(66,10,0)=w0_buffer(66,10,0)
!      w0(64,20,0)=w0_buffer(64,20,0)
!      w0(65,16,0)=w0_buffer(65,16,0)
!      w0(66,11,0)=w0_buffer(66,11,0)
!      w0(66,14,0)=w0_buffer(66,14,0)
!      w0(64,9,0)=w0_buffer(64,9,0)
!      w0(65,9,0)=w0_buffer(65,9,0)
!      w0(64,8,0)=w0_buffer(64,8,0)
!      w0(68,14,0)=w0_buffer(68,14,0)
!      w0(66,12,0)=w0_buffer(66,12,0)
!      w0(66,4,0)=w0_buffer(66,4,0)
!      w0(64,5,0)=w0_buffer(64,5,0)
!      w0(67,12,0)=w0_buffer(67,12,0)
!      w0(67,13,0)=w0_buffer(67,13,0)
!      w0(67,10,0)=w0_buffer(67,10,0)
!      w0(65,4,0)=w0_buffer(65,4,0)
!      w0(68,11,0)=w0_buffer(68,11,0)
!      w0(64,12,0)=w0_buffer(64,12,0)
!      w0(62,1,0)=w0_buffer(62,1,0)
!      w0(68,12,0)=w0_buffer(68,12,0)
!      w0(66,13,0)=w0_buffer(66,13,0)
!      w0(62,9,1)=w0_buffer(62,9,1)
!      w0(67,14,0)=w0_buffer(67,14,0)
!      w0(62,10,1)=w0_buffer(62,10,1)
!      w0(63,6,0)=w0_buffer(63,6,0)
!      w0(65,7,0)=w0_buffer(65,7,0)
!      w0(68,13,0)=w0_buffer(68,13,0)
!      w0(68,15,0)=w0_buffer(68,15,0)
!      w0(70,9,0)=w0_buffer(70,9,0)
!      w0(64,11,0)=w0_buffer(64,11,0)
!      w0(64,3,1)=w0_buffer(64,3,1)
!      w0(65,17,0)=w0_buffer(65,17,0)
!      w0(64,7,0)=w0_buffer(64,7,0)
!      w0(68,16,0)=w0_buffer(68,16,0)
!      w0(67,9,0)=w0_buffer(67,9,0)
!      w0(64,19,0)=w0_buffer(64,19,0)
!      w0(67,11,0)=w0_buffer(67,11,0)
!      w0(62,12,1)=w0_buffer(62,12,1)
!      w0(62,11,1)=w0_buffer(62,11,1)
!      w0(65,3,1)=w0_buffer(65,3,1)
!      w0(64,6,0)=w0_buffer(64,6,0)
!      w0(68,17,0)=w0_buffer(68,17,0)
!      w0(65,5,1)=w0_buffer(65,5,1)
!      w0(64,2,0)=w0_buffer(64,2,0)
!      w0(67,15,0)=w0_buffer(67,15,0)
!      w0(67,7,0)=w0_buffer(67,7,0)
!      w0(67,8,0)=w0_buffer(67,8,0)
!      w0(70,8,0)=w0_buffer(70,8,0)
!      w0(66,15,0)=w0_buffer(66,15,0)
!      w0(68,19,0)=w0_buffer(68,19,0)


    write(*,*) 'Filtering Done'




!! w0(maxl*3,num_vmode,0:num_zw)
!
!! filtering vertical modes
!    if(vmode_e .gt. vmode_s) then
!       write(*,*) 'Filtering vertical modes', vmode_s,' to ',vmode_e
!    do m = vmode_s, vmode_e
!        w0(:,m,:) = (0.0d0,0.0d0)
!    end do
!    else 
!       write(*,*) 'No filtering of vertical modes'
!    end if
!! ej-filter out all the unnecessary vertical mode
!!    do m = 16,20
!!       w0(:,m,:) = (0.0d0,0.0d0)
!!    end do
!! filtering meridional modes
!    if((eig_n_e .gt. eig_n_s) .or. (wig_n_e .gt. wig_n_s) .or. (rot_n_e .gt. rot_n_s)) then
!       write(*,*) 'Filtering meridional modes: EIG ', eig_n_s,' to ',eig_n_e, & 
!               ', WIG ',wig_n_s,' to ',wig_n_e, &
!               ' and ROT ',rot_n_s,' to ',rot_n_e
!    do n = eig_n_s, eig_n_e
!        w0(n,:,:) = (0.0d0,0.0d0)  
!    end do
!!    write(*,*) w0(0,:,:)
!!    w0(1,:,:) = (0.0d0,0.0d0)
!    do n = (maxl+wig_n_s), (maxl+wig_n_e)
!        w0(n,:,:) = (0.0d0,0.0d0)
!    end do
!!    w0(maxl+1,:,:) = (0.0d0,0.0d0)
!    do n = (2*maxl+rot_n_s), (2*maxl+rot_n_e)
!        w0(n,:,:) = (0.0d0,0.0d0)
!    end do  
!!    w0(2*maxl+1,:,:) = (0.0d0,0.0d0)
!!    w0(2*maxl+2,:,:) = (0.0d0,0.0d0) 
!!    w0(2*maxl+3,:,:) = (0.0d0,0.0d0)
!    else
!       write(*,*) 'No filtering of meridional modes'     
!    end if
!!! ej-filter out all the unnecessary meri mode
!!    do m = 1,5
!!       w0(m,:,:) = (0.0d0,0.0d0)
!!    end do
!
!!    write(*,*) w0(1,:,:)    
!
!! filtering zonal modes
!    if(kmode_e .gt. kmode_s) then
!       write(*,*) 'Filtering zonal wavenumbers ', kmode_s,' to ',kmode_e
!    do k = kmode_s, kmode_e
!        w0(:,:,k) = (0.0d0,0.0d0)
!    end do
!!    w0(:,:,0)=(0.0d0,0.0d0)
!!    w0(:,:,2)=(0.0d0,0.0d0)
!    else
!       write(*,*) 'No filtering of zonal modes'
!    end if
!!###### filtering constant field
!     w0(maxl+1,:,0)=(0.0d0,0.0d0)


  end subroutine filter_3DNMFcoefdata
  !----------------------------------------------------------------------------
  subroutine sigma2hybrid

  
  end subroutine sigma2hybrid
  !----------------------------------------------------------------------------    
  subroutine NMFinv_output(itime)
    use mod_adm
    use mod_time, only: &
        cdate
    implicit none
    integer, intent(in) :: itime
    character(len=256) :: tmp_fname

    if(output_inv) then
      if(separate_time) then
        tmp_fname=trim(inverse_fname)//trim(cdate)
        open(1,file=trim(tmp_fname),form='unformatted',access='direct', &
             recl=4*nx*my*mp)
        write(1,rec=1) sngl(u_input)
        write(1,rec=2) sngl(v_input)
        write(1,rec=3) sngl(z_input)
        close(1)
      else
        open(1,file=trim(inverse_fname),form='unformatted',access='direct', &
             recl=4*nx*my*mp)
        write(1,rec=3*(itime-1)+1) sngl(u_input)
        write(1,rec=3*(itime-1)+2) sngl(v_input)
        write(1,rec=3*(itime-1)+3) sngl(z_input)
        close(1)
      end if
    end if

  end subroutine NMFinv_output
  !----------------------------------------------------------------------------
end module mod_normal_inverse
