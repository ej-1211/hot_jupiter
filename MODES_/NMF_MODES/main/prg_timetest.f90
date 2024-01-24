  program timetest
    use mod_adm
    use mod_time

    implicit none
!    integer :: itime

    call init_time


    do

      if(iyear==eyear .and. imon==emon .and. iday==eday .and. ihour==ehour &
         .and. imins==emins .and. ilen==elen) exit

      call time_advance

    end do

  end program timetest
