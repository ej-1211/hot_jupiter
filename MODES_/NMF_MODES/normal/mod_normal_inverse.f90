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
! w0(maxl*3,num_vmode,0:num_zw)

! filtering vertical modes
    if(vmode_e .gt. vmode_s) then
       write(*,*) 'Filtering vertical modes', vmode_s,' to ',vmode_e
    do m = vmode_s, vmode_e
        w0(:,m,:) = (0.0d0,0.0d0)
    end do
    else 
       write(*,*) 'No filtering of vertical modes'
    end if
    
! filtering meridional modes
    if((eig_n_e .gt. eig_n_s) .or. (wig_n_e .gt. wig_n_s) .or. (rot_n_e .gt. rot_n_s)) then
       write(*,*) 'Filtering meridional modes: EIG ', eig_n_s,' to ',eig_n_e, & 
               ', WIG ',wig_n_s,' to ',wig_n_e, &
               ' and ROT ',rot_n_s,' to ',rot_n_e
    do n = eig_n_s, eig_n_e
        w0(n,:,:) = (0.0d0,0.0d0)
    end do
    do n = (maxl+wig_n_s), (maxl+wig_n_e)
        w0(n,:,:) = (0.0d0,0.0d0)
    end do
    do n = (2*maxl+rot_n_s), (2*maxl+rot_n_e)
        w0(n,:,:) = (0.0d0,0.0d0)
    end do   
    else
       write(*,*) 'No filtering of meridional modes'     
    end if
    
! filtering zonal modes
    if(kmode_e .gt. kmode_s) then
       write(*,*) 'Filtering zonal wavenumbers ', kmode_s,' to ',kmode_e
    do k = kmode_s, kmode_e
        w0(:,:,k) = (0.0d0,0.0d0)
    end do
    else
       write(*,*) 'No filtering of zonal modes'
    end if
!###### filtering constant field
     w0(maxl+1,:,0)=(0.0d0,0.0d0)
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