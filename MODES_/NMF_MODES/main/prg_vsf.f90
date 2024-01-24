  program anl_vsf
    use mod_adm
    use mod_vsf_driver
    implicit none

    call vsfcalc_init

    call vsf_driver

  end program anl_vsf
