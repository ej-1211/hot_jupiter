  program expansion3D
    use mod_adm
    use mod_time
    use mod_normal
    use mod_const
    use mod_vsf_driver
    use mod_fft
    use mod_hough_driver
		use mod_write_netcdf   !<-- CHANGED in DEC 2014 by Marten

    implicit none

    call vsf_init

    call normal_init

    call init_time

    do

      call read_inputdata

      call vertical_expansion
      
      call fft_forward

      call hough_expansion

!<-- CHANGED in DEC 2014 by Marten
			if(savenc) then
				call write_netcdf_proj(ncname,cdate,w0,maxl,num_vmode,num_zw+1)
      !else
      	!call NMFinv_output(itime)
			end if
!---> END OF CHANGE

      call NMF_output(itime)
  
      if ( savehoughwind ) then
       call NMF_outputw(itime)
			 if(savenc) then
				call write_netcdf_proj(coef3DNMF_windfname,cdate,ww0,maxl,num_vmode,num_zw+1)
      !else
      	!call NMF_outputw(itime)
			 end if
      end if 

      if(iyear==eyear .and. imon==emon .and. iday==eday .and. ihour==ehour &
         .and. imins==emins .and. ilen==elen) exit

      call time_advance

    end do

  end program expansion3D
