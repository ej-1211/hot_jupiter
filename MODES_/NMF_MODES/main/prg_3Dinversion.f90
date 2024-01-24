  program inversion3D
    use mod_adm
    use mod_time
    use mod_normal_inverse
    use mod_const
    use mod_vsf_driver
    use mod_fft
    use mod_hough_driver
    use mod_interpolation
		use mod_write_netcdf   !<-- CHANGED in DEC 2014 by Marten
    
    implicit none
!    integer :: ip, ix, iy
 
    call vsf_init

    call normal_inverse_init

    call init_time

    do

      call read_3DNMFcoefdata(itime)
      
      call filter_3DNMFcoefdata
      
      call hough_inverse

      call fft_backward

      call vertical_inverse
      
!    print *, 'Contents of z_input array:'
!    do ip = 1, mp
!        print *, 'Slice:', ip
!        do ix = 1, nx
!            print *, (z_input(ip, ix, iy), iy = 1, ny)
!        end do
!    end do
     if(saveasci) then
     call output_3d_gmt(afname,aformat,cdate,u_input,nx,my,mp)  
      	call output_3d_3var_gmt(afname,aformat,cdate,u_input,v_input,z_input,nx,my,mp)
      end if
      
!<-- CHANGED in DEC 2014 by Marten
			if(savenc) then
				call write_netcdf_inv(ncname,cdate,u_input,v_input,z_input,nx,my,mp)
      !else
      	!call NMFinv_output(itime)
			end if
!---> END OF CHANGE
      
      call NMFinv_output(itime)
            
!      if(inv2hybrid) then
!        call sigma2hybrid_ecmwf(cdate)
!      call output_3d_gmt('KWh_uwind_',cdate,u_input,nx,my,mp) 
!      call output_3d_3var_gmt('InvFlow_hybrid',cdate,u_input,v_input,z_input,nx,my,mp) 
!      end if  
           
      if(iyear==eyear .and. imon==emon .and. iday==eday .and. ihour==ehour &
         .and. imins==emins .and. ilen==elen) exit

      call time_advance

    end do

  end program inversion3D
