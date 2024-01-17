  program global_average
    integer :: nx=144
    integer :: ny= 73
    integer :: np= 17 
    real(4), allocatable :: idata(:,:,:)

    character(128) :: globavecnf_fname = 'globave.cnf'
    character(128) :: ygrid_fname = ''
    character(128) :: ifile_fname = ''
    character(128) :: ofile_fname = ''

    character(128) :: ygrid = 'latitude'

    namelist / globave / &
        nx, ny, np,      &
        ifile_fname,     &
        ofile_fname,     &
        ygrid,           &
        ygrid_fname     

  end program global_average
