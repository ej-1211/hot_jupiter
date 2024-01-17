  program main
    use mod_hort_struct
    use mod_adm
    use mod_smode
    implicit none

    call houghcalc_init

    if(trim(ks_mode)=='S') then
      call shige_init
    end if

    call hort_struct

  end program main
