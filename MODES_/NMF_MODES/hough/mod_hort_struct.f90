module mod_hort_struct

  use mod_const, only : r8

  implicit none
  private

  !--- non-dimensional paramter 
  real(r8),              public, save :: eps

  !real(r8), allocatable, public, save :: gusl(:)
  !real(r8), allocatable, public, save :: gusw(:)
  real(r8), allocatable, public, save ::  sia(:)
  real(r8), allocatable, public, save ::  rad(:)

  !--- co-latitude in radian
  real(r8), allocatable, public, save :: an(:)
  !--- co-latitude in degree
  real(r8), allocatable, public, save :: ang(:)
  !--- latitude in radian
  real(r8), allocatable, public, save :: phi(:)
  !--- latitude in degree
  real(r8), allocatable, public, save ::  deg(:)
  !--- number of points between the pole and equator
  integer, public, save :: klast

  !--- number of modes in each type
  integer, public, save :: maxl
  !--- eigenfrequency for eastward ig
  real(r8), allocatable, public, save :: eastgs(:)
  !--- eigenfrequency for westward ig
  real(r8), allocatable, public, save :: westgs(:)
  !--- eigenfrequency for balanced motion
  real(r8), allocatable, public, save :: rotats(:)
  !--- 
  real(r8), allocatable, public, save :: u(:)
  !--- 
  real(r8), allocatable, public, save :: v(:)
  !--- 
  real(r8), allocatable, public, save :: z(:)
  !--- 
  real(r8), allocatable, public, save :: h0(:,:,:)
  !--- 
  real(r8), allocatable, public, save :: strmfn(:)
  !--- 
  real(r8), allocatable, public, save :: velpot(:)
  !--- 
  real(r8), allocatable, public, save :: dvrgnc(:)
  !--- 
  real(r8), allocatable, public, save :: vrtcty(:)
  !--- 
  real(r8), allocatable, public, save :: grdntu(:)
  !--- 
  real(r8), allocatable, public, save :: grdntv(:)
  !--- 
  real(r8), allocatable, public, save :: grdntz(:)
  !--- 
  real(r8), allocatable, public, save :: hftbl(:,:,:)
  !--- 
  integer, allocatable, public, save :: iparj(:)
  !--- 
  real(r8), allocatable, public, save :: sigmb(:)
  !--- workspace at least 33*max0(20,maxl,int(dsqrt(eps)))
  real(r8), allocatable, public, save :: w(:)
  !--- for checking status
  integer :: ierr

!  character(len=128) :: ygrid_fname = './grid.dat'

  !----------------------------------------------------------------------------
  public  :: hort_struct
  
  !----------------------------------------------------------------------------
contains  
  !----------------------------------------------------------------------------
  subroutine hort_struct
    use mod_adm
    use mod_gaussg, only : &
        gaussg,            &
        gauss_init,        &
        gusl,              &
        gusw,              &
        an,                &
        ang,               &
        phi,               &
        deg
    use mod_sigma,  only : &
        sigma
    use mod_smode
    use mod_abcoef
    use mod_uvzder
    use mod_const
    use mod_output

    implicit none

    integer :: lp, lpst
    integer :: l, m
    integer :: i
    integer :: n
    integer :: iewr
    real(r8) :: t1
    real(r8) :: sig

    real(r8), allocatable :: a(:), b(:), c(:)

    call gauss_init

    allocate(  sia(mk) )
    allocate(  rad(mk) )

    allocate(      u(mk) )
    allocate(      v(mk) )
    allocate(      z(mk) )
    allocate( h0(2*mk,maxl,3) )
    allocate( strmfn(mk) )
    allocate( velpot(mk) )
    allocate( dvrgnc(mk) )
    allocate( vrtcty(mk) )
    allocate( grdntu(mk) )
    allocate( grdntv(mk) )
    allocate( grdntz(mk) )

    allocate( iparj(3*maxl) )
    allocate( sigmb(3*maxl) )

    allocate( hftbl(mk,3*maxl,4) )

    klast = mk
!
!   Computation of Gaussian latitude and weight
!
!    call gaussg( klast, sia, rad )

    !write(*,*) "Output grid information to => ", trim(ygrid_fname)
    !open(21,file=trim(ygrid_fname),form='unformatted',access='sequential')
    !write(21) gusl, gusw, phi, deg 
    !close(21)

    do i = 1,mk
      gusl(i) = ygrid(my-i+1)
      gusw(i) = ygrid_weight(my-i+1)
      rad(i)  = dacos(gusl(i))
      sia(i)  = dsin(rad(i))
      an(i)   = dacos(gusl(i))
      ang(i)  = an(i)/dtr
      phi(i)  = hfpai - an(i)
      deg(i)  = phi(i) / dtr
      write(*,*) i, gusl(i), gusw(i)
    end do

    allocate( eastgs(0:maxl) )
    allocate( westgs(0:maxl) )
    allocate( rotats(0:maxl) )

    do ls = szw, ezw        ! loop over zonal wavenumbers
      do m = 1, num_vmode   ! loop over vertical modes
        eps =  4.0_r8 * ( omega * ae ) ** 2 / ( g * evht(m) )
        write(*,'("              eps =",f15.5)') eps
        write(*,'("            omega =",f15.5)') omega
        write(*,'("               ae =",f15.5)') ae
        write(*,'("                g =",f15.5)') g
        write(*,'("          evht(m) =",f15.5)') evht(m)
        
        t1  = dsqrt(eps)
        !n   = max0( 20, 8*maxl )
        n   = max0( 20, maxl, int(dsqrt(eps)) )
        allocate( w(33*n) )
        allocate( a(2*n) )
        allocate( b(2*n) )
        allocate( c(2*n) )

        write(*,'(i4, "  Equivalent Height =",f15.5," (m)")') m, evht(m)
        write(*,'("              eps =",f15.5)') eps
        write(*,'("The square of eps =",f15.5)') t1
        write(*,'("Size of workspace =",i9)'  ) n
!
!   Computation of eigenfrequency
!
        call sigma( ls, maxl, 3*n, ierr, eps, eastgs, westgs, rotats, w )

        if( ierr /= 0 ) write(*,*) 'The error status', ierr
        open(31,file=trim(freq_fname),form='unformatted',access='sequential')
        write(31) eastgs(1:maxl), westgs(1:maxl), rotats(1:maxl)
        close(31)
!
!   Computation of Hough function
!
!
! For eastward inertio-gravity modes
!
        lpst = 0
        do lp = 1,maxl
          lpst = lpst + 1
          l    = lp - 1 
          iewr = 1
          sig  = eastgs(l+1)
          w    = 0.0_r8
          call abcoef( ls, maxl, l, iewr, sig, eps, 2*n, a, b, c, w(1:10*n) )
          call uvzder( ls, maxl, l, iewr, n, eps, klast, phi, a, b, c, u, v, z, &
                      strmfn, velpot, dvrgnc, vrtcty, grdntu, grdntv, grdntz, w )
          if(mod(l,2)==0) then
             iparj(lpst) = 2
          else
             iparj(lpst) = 1
          end if
          sigmb(lpst) = sig
          do i = 1, mk
             hftbl(i,lpst,1) = u(i)      
             hftbl(i,lpst,2) = v(i)      
             hftbl(i,lpst,3) = z(i)      
             hftbl(i,lpst,4) = grdntv(i)      
          end do
        end do
!
! For westward inertio-gravity modes
!
        do lp = 1,maxl
          lpst = lpst + 1
          l    = lp - 1
          iewr = 2
          sig   = westgs(l+1)
          call abcoef( ls, maxl, l, iewr, sig, eps, 2*n, a, b, c, w(1:10*n) )
          call uvzder( ls, maxl, l, iewr, n, eps, klast, phi, a, b, c, u, v, z, &
                       strmfn, velpot, dvrgnc, vrtcty, grdntu, grdntv, grdntz, w )
          if(mod(l,2)==0) then
            iparj(lpst) = 2
          else
            iparj(lpst) = 1
          end if
          sigmb(lpst) = sig
          do i = 1, mk
            hftbl(i,lpst,1) = u(i)
            hftbl(i,lpst,2) = v(i)
            hftbl(i,lpst,3) = z(i)
            hftbl(i,lpst,4) = grdntv(i)
          end do
        end do
!
! For rotational modes
!
        if(ks_mode == 'K' .or. ls /= 0 ) then
          do lp = 1,maxl
            lpst = lpst + 1
            l    = lp - 1
            iewr = 3
            sig  = rotats(l+1)
            call abcoef( ls, maxl, l, iewr, sig, eps, 2*n, a, b, c, w(1:10*n) )
            call uvzder( ls, maxl, l, iewr, n, eps, klast, phi, a, b, c, u, v, z, &
                         strmfn, velpot, dvrgnc, vrtcty, grdntu, grdntv, grdntz, w )
            if(mod(l,2)==0) then
              iparj(lpst) = 1
            else
              iparj(lpst) = 2
            end if
            sigmb(lpst) = sig
            do i = 1, mk
              hftbl(i,lpst,1) = u(i)
              hftbl(i,lpst,2) = v(i)
              hftbl(i,lpst,3) = z(i)
              hftbl(i,lpst,4) = grdntv(i)
            end do
          end do
        else if(ks_mode == 'S' .and. ls == 0) then
          call smtrx_smode(eps)
          call shige(eps, h0)
          !call shige(eps, phi, h0)
          do lp = 1,maxl
            lpst = lpst + 1
            l    = lp - 1
            if(mod(l,2)==0) then
              iparj(lpst) = 2
            else
              iparj(lpst) = 1
            end if
            !if(lp==1) iparj(lpst) = 1
            if(lp==1) iparj(lpst) = 2
            do i = 1, mk
              hftbl(i,lpst,1) = h0(2*mk-i+1,lp,1)
              hftbl(i,lpst,2) = h0(2*mk-i+1,lp,2)
              hftbl(i,lpst,3) = h0(2*mk-i+1,lp,3)
            end do
          end do
        end if

        if( ocheck ) then
          call cekoth( mk, maxl*3, mk*2, hftbl, klast, maxl*3, gusw, iparj)
        end if

        if( output_gmt ) then
          call output_for_gmt( ls, m, hftbl, iparj ) 
        end if

        call output_for_binary( ls, m, hftbl, iparj )

        call output_for_netcdf( ls, m, hftbl, iparj )

        !write(*,*) 'Structure of Hough functions'
        !do i = 1,mk
        !  write(*,'(i5,3F20.15)') i, (hftbl(i,maxl+1,l),l=1,3)
        !end do

        deallocate( w )
        deallocate( a )
        deallocate( b )
        deallocate( c )

      end do
    end do

  end subroutine hort_struct
  !----------------------------------------------------------------------------
  subroutine cekoth( mk, ld, mj, hftbl, klast, n, gw, ipar )
    use mod_adm, only : &
        ocheck
    implicit none
    integer, intent(in) :: mk, ld, mj, klast, n
    integer, intent(in) :: ipar(ld)
    real(r8), intent(in) :: hftbl(mk,ld,4)
    real(r8), intent(in) :: gw(mk)
    
    integer :: jend, jdp
    integer :: j, k, ic

    real(r8) :: gusw(mj)
    real(r8) :: ufi(mj), vfi(mj), zfi(mj) 
    real(r8) :: ufj(mj), vfj(mj), zfj(mj) 
    real(r8) :: usq, vsq, zsq, asq
    real(r8) :: orth(ld,ld)

    !write(*,*) '###################################################'
    !write(*,*) '#### Check the orthogonality of Hough Function ####'
    !write(*,*) '###################################################'

    orth(:,:) = 0.0_r8
    jend = 2 * klast
    jdp  = jend + 1
    do j = 1, klast
      gusw(j) = gw(j)
      gusw(jdp-j) = gw(j)
    end do

    !do j = 1, jend
    !   write(*,*) j, gusw(j)
    !end do

    do ic = 1, n

      if(ipar(ic)==1) then

        do j = 1, klast
          ufi(jdp-j) =  hftbl(j,ic,1)
          ufi(j)     = -hftbl(j,ic,1)
          vfi(jdp-j) =  hftbl(j,ic,2)
          vfi(j)     =  hftbl(j,ic,2)
          zfi(jdp-j) =  hftbl(j,ic,3)
          zfi(j)     = -hftbl(j,ic,3)
        end do

        do k = ic, n
          if(ipar(k)==1) then
            do j = 1, klast
              ufj(jdp-j) =  hftbl(j,k,1)
              ufj(j)     = -hftbl(j,k,1)
              vfj(jdp-j) =  hftbl(j,k,2)
              vfj(j)     =  hftbl(j,k,2)
              zfj(jdp-j) =  hftbl(j,k,3)
              zfj(j)     = -hftbl(j,k,3)
            end do
          else
            do j = 1, klast
              ufj(jdp-j) =  hftbl(j,k,1)
              ufj(j)     =  hftbl(j,k,1)
              vfj(jdp-j) =  hftbl(j,k,2)
              vfj(j)     = -hftbl(j,k,2)
              zfj(jdp-j) =  hftbl(j,k,3)
              zfj(j)     =  hftbl(j,k,3)
            end do
          end if

          usq = 0.0_r8
          vsq = 0.0_r8
          zsq = 0.0_r8
          do j = 1, jend
            usq = usq + gusw(j) * ufi(j) * ufj(j)
            vsq = vsq + gusw(j) * vfi(j) * vfj(j)
            zsq = zsq + gusw(j) * zfi(j) * zfj(j)
          end do
          asq = usq + vsq + zsq
          orth(k,ic) = asq
          !write(*,'(2i5,4e15.7)') ic, k, usq, vsq, zsq, asq
        end do

      else ! ipar

        do j = 1, klast
          ufi(jdp-j) =  hftbl(j,ic,1)
          ufi(j)     =  hftbl(j,ic,1)
          vfi(jdp-j) =  hftbl(j,ic,2)
          vfi(j)     = -hftbl(j,ic,2)
          zfi(jdp-j) =  hftbl(j,ic,3)
          zfi(j)     =  hftbl(j,ic,3)
        end do

        do k = ic, n
          if(ipar(k)==1) then
            do j = 1, klast
              ufj(jdp-j) =  hftbl(j,k,1)
              ufj(j)     = -hftbl(j,k,1)
              vfj(jdp-j) =  hftbl(j,k,2)
              vfj(j)     =  hftbl(j,k,2)
              zfj(jdp-j) =  hftbl(j,k,3)
              zfj(j)     = -hftbl(j,k,3)
            end do
          else
            do j = 1, klast
              ufj(jdp-j) =  hftbl(j,k,1)
              ufj(j)     =  hftbl(j,k,1)
              vfj(jdp-j) =  hftbl(j,k,2)
              vfj(j)     = -hftbl(j,k,2)
              zfj(jdp-j) =  hftbl(j,k,3)
              zfj(j)     =  hftbl(j,k,3)
            end do
          end if

          usq = 0.0_r8
          vsq = 0.0_r8
          zsq = 0.0_r8
          do j = 1, jend
            usq = usq + gusw(j) * ufi(j) * ufj(j)
            vsq = vsq + gusw(j) * vfi(j) * vfj(j)
            zsq = zsq + gusw(j) * zfi(j) * zfj(j)
          end do
          asq = usq + vsq + zsq
          orth(k,ic) = asq
          !write(*,'(2i5,4e15.7)') ic, k, usq, vsq, zsq, asq
        end do
      end if
    end do

    !if( ocheck ) then
    !  do ic = 1,n
    !    write(*,'(60f10.6)') (orth(k,ic),k=ic,n)
    !  end do
!
!      do ic = 1,n
!        write(*,*) ic, ipar(ic)
!      end do
!
!      do ic = 1,n
!         do k = ic, n
!            if(orth(k,ic).gt.1.0d-10) then
!               write(*,'(2i5,e20.10)') ic,k,orth(k,ic)
!            end if
!         end do
!      end do
    !end if

    do k = 1, n
      if( dabs( orth(k,k) - 1.0_r8 ) .gt. 1.0d-10 ) then
        write(*,*) k, orth(k,k)
      end if
    end do 

    do ic = 1, n
      do k = ic+1, n
        if( dabs( orth(k,ic) ) .gt. 1.0d-12 ) then
          write(*,*) ic, k, orth(k,ic)
        end if
      end do
    end do


  end subroutine cekoth
  !----------------------------------------------------------------------------
end module mod_hort_struct
