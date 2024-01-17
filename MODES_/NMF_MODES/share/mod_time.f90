module mod_time
  implicit none
  
  integer, public, save :: itime
  !--- Begining time
  integer, public, save :: syear
  integer, public, save ::  smon
  integer, public, save ::  sday
  integer, public, save :: shour
  integer, public, save :: smins 
  integer, public, save ::  ssec
  integer, public, save ::  slen
  !--- End time
  integer, public, save :: eyear
  integer, public, save ::  emon
  integer, public, save ::  eday
  integer, public, save :: ehour
  integer, public, save :: emins 
  integer, public, save ::  esec
  integer, public, save ::  elen
  !--- loop index
  integer, public, save :: iyear
  integer, public, save ::  imon
  integer, public, save ::  iday
  integer, public, save :: ihour
  integer, public, save :: imins 
  integer, public, save ::  isec
  integer, public, save ::  ilen
  !--- date name of input file
  character(len=10), public, save :: datetype 
  ! 'yyyymmddhh' 
  ! 'forclength' 
 
!--- Current date
!!  character(len=12), public, save :: cdate
  character(len=15), public, save :: cdate
  !--- Time interval in seconds (s)
  integer, public, save :: dt

  integer, private, save :: monday ( 12,2 )
  data   monday /                            &
       31,28,31,30,31,30,31,31,30,31,30,31,  &
       31,29,31,30,31,30,31,31,30,31,30,31 /

  public :: init_time
  public :: time_advance

contains
  !-----------------------------------------------------------------------------
  subroutine init_time
    use mod_adm
    implicit none
    integer :: i

    namelist / time /  &
      datetype,        &
      syear, smon, sday, shour, smins, ssec, slen, &
      eyear, emon, eday, ehour, emins, esec, elen, &
      dt

    itime = 1    
    
    open(1,file=trim(normalcnf_fname))
    read(1,time)
    write(*,time)
    ! TJH FIXME ... should we close unit 1 at this point
  
    iyear = syear
    imon  = smon
    iday  = sday
    ihour = shour
    imins = smins
    isec  = ssec
    ilen  = slen

    if(trim(datetype)=='yyyymmddhh') then
      write(cdate,'(i4.4,4i2.2,i3.3)') iyear,imon,iday,ihour,imins,ilen
      do i = 1,15
        if(cdate(i:i)==' ') cdate(i:i)='0'
      end do 
      write(*,*) 'Initial cdate: ',cdate
      
    elseif(trim(datetype)=='forclength') then  
      write(cdate,'(i4.4,4i2.2,i3.3)') iyear,imon,iday,ihour,imins,ilen
      do i = 1,15
        if(cdate(i:i)==' ') cdate(i:i)='0'
      end do 
      write(*,*) 'Initial cdate: ',cdate    
    end if

  end subroutine init_time
  !-----------------------------------------------------------------------------
  subroutine time_advance
    implicit none
!    integer :: iyearin
!    integer :: imonin
    integer :: idayin
    integer :: ihourin
    integer :: iminsin
    integer :: isecin
    integer :: ilenin
    integer :: i

    ! TJH FIXME ...  iyearin, imonin are not used ...

    if(trim(datetype)=='yyyymmddhh') then
    
    idayin = int(dt/60/60/24)
    ihourin = int(dt/60/60) - idayin*24
    iminsin = int(dt/60) - (idayin*24 + ihourin)*60
    isecin = dt - int(dt/60) * 60
      
    iday=iday+idayin
    ihour=ihour+ihourin
    imins=imins+iminsin
    isec=isec+isecin

    if(isec>=60) then
      isec=isec-60
      imins=imins+1
    end if
    if(imins>=60) then
      imins=imins-60
      ihour=ihour+1
    end if
    if(ihour>=24) then
      ihour=ihour-24
      iday=iday+1
    end if
    if(mod(iyear,4)==0 .and. mod(iyear,100)/=0) then
      if(iday>monday(imon,2)) then
        imon=imon+1
        iday=1
      end if
    else
      if(iday>monday(imon,1)) then
        imon=imon+1
        iday=1
      end if
    end if
    if(imon>12) then
      imon=1
      iyear=iyear+1
    end if
    write(cdate,'(i4.4,4i2.2,i3.3)') iyear,imon,iday,ihour,imins,ilen
     do i = 1,15
      if(cdate(i:i)==' ') cdate(i:i)='0'
     end do
    write(*,*) 'Updated cdate: ',cdate

    elseif(trim(datetype)=='forclength') then

      ilenin=int(dt/60/60)
      ilen=ilen+ilenin 
      write(cdate,'(i4.4,4i2.2,i3.3)') iyear,imon,iday,ihour,imins,ilen
      do i = 1,15
        if(cdate(i:i)==' ') cdate(i:i)='0'
      end do
      write(*,*) 'Updated cdate: ',cdate
         
    end if
    
    itime=itime+1
         
  end subroutine time_advance
  !-----------------------------------------------------------------------------
end module mod_time
