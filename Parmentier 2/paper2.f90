      PROGRAM     paper2
      IMPLICIT NONE
!!!!!!INPUT VARIABLES!!!!!!!!!!!!!!!!
      !Model parameters 
      REAL*8 :: Pmin,Pmax !minimal and maximal pressure of the model. 
      INTEGER, PARAMETER :: N=100 ! Number of levels in the model
      CHARACTER*150 :: OutputFile ! Name of the output file
      CHARACTER*150 :: OutputFileCSV ! Name of the output file
      !Planet parameters
      REAL*8 :: Tint !Temperature for the internal flux
      REAL*8 :: Teq0 !Temperature of the incoming flux
      REAL*8 :: mu !Cosine of the stellar angle mu=1/sqrt(3) for averaged profiles
      REAL*8 :: f !Parameter for averaged profiles f=0.5 for dayside average and f=0.25 for planet-average
      REAL*8 :: grav !Gravity of the planet (m/s**2)

      !Ej parameter
      REAL*8 :: Omega
      REAL*8 :: R_planet !planet radius
      REAL*8 :: M_planet !planet radius
      REAL*8 :: R
      REAL*8 :: cp
      REAL*8 :: notyet
      REAL*8 :: G
      REAL*8 :: M_star
      REAL*8 :: R_star
      REAL*8 :: T_star
      REAL*8 :: C1
      REAL*8 :: C3

      !For normalization





      !Options
      CHARACTER*10 :: ROSS ! Can be either 'USER' when Kappa is defined by the user or 'AUTO' to use the opacities of Freedman et al. 2008 as fitted by Valencia et al. 2013
      CHARACTER*10 :: COEFF ! Can be either 'USER' when the input coefficients for Gv1,Gv2,Betav,Beta,Gp or 'AUTO' to use the fit provided by Parmentier et al. 2014
      CHARACTER*10 :: COMP ! Can be 'SOLAR' for a solar-composition atmosphere or 'NOTIO' for an atmosphere without TiO/VO.
      CHARACTER*10 :: STAR ! Can only be 'SUN' for now. 
      CHARACTER*10 :: ALBEDO !If 'USER' uses the value of Ab, if 'AUTO', uses the fit of Parmentier at al. 2014.
      CHARACTER*10 :: CONV !If YES calculates a radiative/convective profile, if NO, calculates only the radiative solution

      !Atmosphere parameters
      REAL*8 :: Gv1 !Gamma_v1=Kv1/Kr is the ratio of the first visible opacity to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8 :: Gv2 !Gamma_v2=Kv2/Kr is the ratio of the second visible opacity to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8 :: Gv3 !Gamma_v2=Kv3/Kr is the ratio of the third visible opacity to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8 :: Gp  !Gamma_P=Kp/Kr ratio of the Planck to the Rosseland mean opacity see Parmentier et al. 2014
      REAL*8 :: Beta ! Relative spectral width of the first thermal band
      REAL*8 :: Ab ! Bond albedo of the planet


!!!!!!INPUTS/OUTPUTS!!!!!!!!!!
      REAL*8, DIMENSION(N) ::Kappa ! Rosseland mean opacity at each atmospheric level, either calculated or set by the user

!!!!!!!OUTPUTS!!!!!!!!!!!!!
      REAL*8 :: Teff ! Effective temperature of the model 
      REAL*8, DIMENSION(N) :: TAU ! Rosseland optical depth
      REAL*8, DIMENSION(N) :: P   ! Pressure
      REAL*8, DIMENSION(N) :: sigma   ! sigma=P/Pmax
      REAL*8, DIMENSION(N) :: rho
      REAL*8, DIMENSION(N) :: T   ! Temperature
      REAL*8, DIMENSION(N) :: gradad ! Adiabatic gradient
      REAL*8, DIMENSION(N) :: gradrad ! Radiative gradients
      REAL*8, DIMENSION(N) :: lpkpt
      REAL*8, DIMENSION(N) :: lpkpr
      REAL*8, DIMENSION(N) :: J1
      REAL*8, DIMENSION(N) :: J2
      REAL*8, DIMENSION(N) :: H1
      REAL*8, DIMENSION(N) :: H2
      REAL*8, DIMENSION(N) :: GAMMA_0
      REAL*8, DIMENSION(N) :: density
      REAL*8, DIMENSION(N) :: t_th



!      REAL*8, DIMENSION(N) :: pk
!      REAL*8, DIMENSION(N) :: pr
      REAL*8 :: PRC ! Pressure at the radiative/convective boudary
      REAL*8 :: Tmu
      REAL*8 :: Tirr
      REAL*8 :: Betav1
      REAL*8 :: Betav2
      REAL*8 :: Betav3
      REAL*8 :: Gamma1
      REAL*8 :: Gamma2
      REAL*8 :: Hv0

!!!!!!!INTERNAL VARIABLES
      INTEGER :: i,K
      CHARACTER*4 :: xtemp

       OutputFile='PTprofile.dat'
       !OutputFileCSV='PTprofile(20bar)(v0911)(normalized).csv'
       OutputFileCSV='PTprofile(20bar)(Hd209458b)(1107).csv'

!!!!!!Model parameters
      G = 6.67E-11
      Pmin=1E1 !Minimal pressure (Pa) surface pressure
      Pmax=2E6 ! Maximal pressure (Pa) 乘以10^-5＝bar 2e6
      R = 8.314/0.00236 !specific gas constant
      cp = R*3.5
      notyet = 123.456


!!!!!!Planet parameters
      Tint=278.4425!Internal temperature (K)
      Teq0=1448!!Equilibrium temperature for zero albedo (K) 1448
      f=0.25 !f=0.25 for a planet-average profile, f=0.5 for a dayside average profile
      mu=1/sqrt(3.0) ! Cosine of the irradiation angle. mu=1/sqrt(3) is the mean mu to use for a dayside or a planet average.
!      Teff0=1000 !Effective temperature for zero albedo
      R_planet = 15.6 * (6371000) !15.6 Earth Radius (m)
      M_planet = 232 * (5.9722E24) !kg
      grav=G*M_planet/(R_planet**2)! gravity (m/s**2) (general version)
!      grav = 25 !initial setting
      write(*,*)'g=',grav


!!!!!!Star parameters
      M_star = 1.069175195875 *(1988500E24) !kg
      R_star = 1.19997599424 * (695700E3) !m
      T_star = 6026.354945475


!!!!!Options
      ROSS="AUTO" !This can be AUTO or USER. If USER, the Rosseland mean opacity must be provided below for each level. If AUTO, the Rosseland opacities are calculated from the fit of Valencia et al. 2013
      COEFF="AUTO" !This can be AUTO or USER whether you want to use the fit of the coefficients provided in Parmentier et al. 2013 or set your own coefficients.
      STAR="SUN"
      COMP="SOLAR"    !This can be SOLAR for a solar-composition atmosphere or NOTIO for an atmosphere without TiO/VO
      ALBEDO="AUTO" !This can be AUTO or USER. If USER, the Bond albedo must be provided below. If AUTO, the Bond Albedo is calculated from the fit of Parmentier et al. 2014
      CONV="NO" ! If "YES", a convective zone is calculated at the bottom of the model using the convective gradient of Parmentier et al. 2014. If "NO" only the radiative solution is provided

!!!!!Input variables, used only if COEFF="USER"
      if(COEFF.eq."USER") then
      Gv1=10 
      Gv2=0.1
      Gv3=0.01 ! 0.5 is two bands of equal width for the absorption of the stellar flux
      Gp=10! Gp=10 Planck opacities are ten times higher than Rosseland ones.
      beta=0.5! The two thermal bands have the same width
      endif

!!!!!Input variables, used only if ALBEDO="USER"
      if(ALBEDO.eq."USER") then
      Ab=0.3!!Bond albedo set by the user if Albedo="USER".
      endif

!!!!!Input variables, used only if ROSS="USER"
      if(ROSS.eq."USER") then
      DO i=1,N
      Kappa(i)=1e-1 ! kg/m^2 Here I set the Rosseland mean opacity, if ROSS="CONSTANT". Kappa(i) can be any functional form of the pressure P(i)
      ENDDO
      endif

!!!!!!!!!!!!!!END OF INPUTS!!!!!!!!!!!!!!!!!!!!!

      !Set the pressure levels
      DO i=1,N
       P(i)=10**(log10(Pmin)+(log10(Pmax)-log10(Pmin))*                 &
     &1.0*(i-1)/(N-1))!The pressure varies from Pmin to Pmax in logspace with N points!     
       sigma(i) = P(i)/Pmax
      ENDDO

      !Calculate the Temperature 
      call tprofile(Tint,Teq0,mu,f,grav,Gv1,Gv2,Gv3,Gp,Beta,Kappa,      &
     &Ab,Teff,ROSS,COEFF,COMP,STAR,CONV,ALBEDO,TAU,P,T,N,PRC,gradad,    &
     &gradrad,Tmu,Betav1,Betav2,Betav3,Gamma1,Gamma2,J1,J2,H1,H2,Tirr,&
     &lpkpt,lpkpr,C1,C3,Hv0)

     ! Calculate Omega
     Omega = sqrt(G*M_star/(R_star*((T_star/Tmu)**2))**3)
     write(*,*) 'D=',R_planet/(R_star*((T_star/Tmu)**2))
     write(*,*) 'mstara=',Omega

      Do i=1,N
        rho(i) = P(i)/(R*T(i))
        GAMMA_0(i) = R*T(i)/(mu*cp*sigma(i))-(T(i)-T(i+1))/(sigma(i)-sigma(i+1))
        density(i) = P(i)/(R*T(i))
        t_th(i) = 1/((16*density(i)*Kappa(i)*5.670374419E-8*T(i)**4/P(i))*(0.4))
      ENDDO



!      DO i=1,N
!        pk(i)=log(kappa(i))-log(kappa(i+1))
!        pr(i)=log(rho(i))-log(rho(i+1))
!        pkpr(i) = pk(i)/pr(i)
!!        pkpr(i)=(log(kappa(i))-log(kappa(i+1)))/(log(rho(i))-log(rho(i+1)))
!      ENDDO
!      print *, log(kappa(110))-log(kappa(111)),(log(rho(110))-log(rho(111))),&
!      &(log(kappa(110))-log(kappa(111)))/(log(rho(110))-log(rho(111))),pk(110),pr(110),pkpr(110)

      !Output formats
 1400    format (f10.2,1X,f10.2,1X,f10.2,1X,f10.2,1X,f10.2,1X,f6.2,1X,  &
      &e9.3,1X,e9.3,1X,f5.4,1X,e9.3,1X,f5.4,1X,f5.4,1X                  & 
      &,A10,A10,A10,A10,A10,A10' Teff,Tint,Teq0,mu,f,grav,Gv1,Gv2,Gv3,Gp&  
      &,Beta,Ab,ROSS,COEFF,COMP,STAR,CONV,ALBEDO')

 1500   format(                                                         &
     &      'c',t1,'(1)N',t6,'(2)P(Pa)',t17,'(3)T(K)'                   &
     &      ,t27,'(4)Tau', t37,'(5)Kr(m^2/kg)', t50, '(6)Gradad'        &
     &      , t60, '(7)Gradrad', t72,'(8)PRC (Pa)')

!9900   format (1X,I3,",",E9.3,",",E9.3,",",F10.3,",",E9.3,",",E9.3,",",F8.5,"," &
!     &,F8.5,",",E9.3,",",E9.3,",",E9.3,",",F5.3,",",E9.3,",",F5.3,",",  &
!     &F6.3,",",F6.1,",",F6.1,",",F6.1,",",F5.3,",",F5.3,",",F8.3,",",   &
!     &A10,",",A10,",",A10,",",A10,",",A10,",",A10)

 9600   format (1X,I3,",",E9.3,",",E10.5,",",F10.3,",",E10.3,",",F6.3,",",F5.3,",", &
     & F10.5,",",F10.5,",",F10.5,",",F6.1,",",E9.3,",",F10.5,",",F10.5,",",F10.5,",",F10.5,&
     & ",",F10.5,",",E9.3,",",E9.3,",",E9.3,",",E10.5,",",E9.3,",",F10.5,",",F10.5,",",E20.10,&
     &",",E20.10,",",E20.10,",",E20.10,",",E20.10,",",E20.10,",",E20.10,",",E20.10)


 9700   format (1X,I3,",",E9.3,",",E9.3,",",F10.3,",",F10.3,",",F10.3,",",F5.3,",", &
     & F10.3,",",F10.5,",",F10.5,",",F10.5,",",E9.3,",",E10.5,",",F10.5,",",F10.5,",",F10.5,&
     & ",",F10.5,",",F10.5,",",E9.3,",",E9.3,",",E10.3,",",F10.3,",",E10.3,",",E10.3,",",E20.5,&
     &",",E20.10,",",E20.10,",",E20.10,",",E20.10,",",E20.10&
     &,",",F10.5,",",E20.10,",",E20.10,",",E20.10)

 9800   format (1X,I3,1X,E9.3,1X,F10.3,1X,E9.3,1X,E9.3,1X,E9.3,1X,      &
     &E9.3,1X,E9.3,1X,E9.3,1X,E9.3,1X,E9.3,1X,E9.3)
      write(xtemp,'(I4.4)') Int((1-Ab)**0.25*(4*f)**0.25*Teq0)
!       OutputFile='./profilesNoTiO/PTprofile-Tmu-'//xtemp//'.dat'



      !Open and write output file
       open(unit=40,form='formatted',status='UNKNOWN',file=OutputFileCSV&
     &)
!       write(40,'(A)')'N,sigma,P(Pa),T(K),Tau,Kappa(m^2/kg),gradad,gradrad,PRC&
!     &,Gv1,Gv2,Gv3,Gp,Beta,Ab,Teff(K),Tint(K),Teq0(K),mu,f,grav(m/s^2), &
!     &ROSS,COEFF,COMP,STAR,CONV,ALBEDO'

!     Desired
!        write(40,'(A)')'N,Omega(s^-1),R_p(m),R(J/(K*g)),cp(J/(mole*K)),g(m/s^2),beta,Betav1,Betav2,Betav3,T_irr(K),P_H(Pa),&
!        &Gv1,Gv2,Gv3,Gamma1,Gamma2,sigma,P(Pa),Kappa(m^2/kg),T(K),rho(g/m^3),pkpt,pkpr,&
!        &J1(W/m^2),J2(W/m^2),H1(W/m^2),H2(W/m^2),tau,mu,Gp,C1,C3'
!
!        DO i=1,N
!         !write(40,9900)i,sigma(i),P(i),T(i),Tau(i),Kappa(i),gradad(i),gradrad(i),PRC,Gv1,Gv2,Gv3,Gp,Beta,Ab,Teff,Tint,Teq0,mu,f,grav,ROSS,COEFF,COMP,STAR,CONV,ALBEDO
!         write(40,9700)i,Omega,R_planet,R,cp,grav,Beta,Betav1,Betav2,Betav3,Tirr, &
!         &Pmax,Gv1,Gv2,Gv3,Gamma1,Gamma2,sigma(i),P(i),Kappa(i),T(i),rho(i)&
!         &,lpkpt(i),lpkpr(i),J1(i),J2(i),H1(i),H2(i),TAU(i),mu,Gp,C1,C3
!        ENDDO

        !normalization
!         write(40,'(A)')'N,Omega(s^-1),R_p(m),R,cp,g,beta,Betav1,Betav2,Betav3,T_irr,P_H,&
!        &Tau,Gv1,Gv2,Gv3,Gamma1,Gamma2,sigma,Kappa,T,rho,pkpt,pkpr,&
!        &J1,J2,H1,H2,tau,mu,Gp,C1,C3,GAMMA_0'
!
!        DO i=1,N
!         write(40,9700)i,Omega,R_planet,R,cp/R,grav/(R_planet*Omega**2),Beta,Betav1,Betav2,Betav3,&
!         &Tirr/(Omega**2*R_planet**2/R), &
!         &Pmax/(Omega**2*R_planet/(Kappa(1))),TAU(i),Gv1,Gv2,Gv3,Gamma1,Gamma2,sigma(i),Kappa(i)/Kappa(1)&
!         &,T(i)/(Omega**2*R_planet**2/R),rho(i)/(1/(R_planet*Kappa(1)))&
!         &,lpkpt(i),lpkpr(i),J1(i)/(R_planet**2*Omega**3/Kappa(1)),J2(i)/(R_planet**2*Omega**3/Kappa(1))&
!         &,H1(i)/(R_planet**2*Omega**3/Kappa(1)),H2(i)/(R_planet**2*Omega**3/Kappa(1)),TAU(i),mu,Gp,C1,C3,&
!         &GAMMA_0(i)/(Omega**2*R_planet**2/R)
!        ENDDO
       !unnormalized()


        write(40,'(A)')'N,Gp,mu,R(J/(K*g)),tau,g(m/s^2),beta,Betav1,Betav2,Betav3,T_irr(K),P_H(Pa),&
        &Gv1,Gv2,Gv3,Gamma1,Gamma2,sigma,P(Pa),Kappa(m^2/kg),T(K),rho,pkpt,pkpr,&
        &J1(W/m^2),J2(W/m^2),H1(W/m^2),H2(W/m^2),C1,C3,GAMMA_0,t_th'
        DO i=1,N
         write(40,9600)i,Gp,mu,R,TAU(i),grav,Beta,Betav1,Betav2,Betav3,Tirr,&
         &Pmax,Gv1,Gv2,Gv3,Gamma1,Gamma2,sigma(i),P(i),Kappa(i)&
         &,T(i),rho(i)&
         &,lpkpt(i),lpkpr(i),J1(i),J2(i)&
         &,H1(i),H2(i),C1,C3,GAMMA_0(i),t_th(i)
        ENDDO


!        write(40,'(A)')'N,Omega(s^-1),R_p(m),R(J/(K*g)),cp(J/(mole*K)),g(m/s^2),beta,Betav1,Betav2,Betav3,T_irr(K),P_H(Pa),&
!        &Gv1,Gv2,Gv3,Gamma1,Gamma2,sigma,P(Pa),Kappa(m^2/kg),T(K),rho(g/m^3),pkpt,pkpr,&
!        &J1(W/m^2),J2(W/m^2),H1(W/m^2),H2(W/m^2),tau,mu,Gp,Hv0,Tint'
!
!        DO i=1,N
!         !write(40,9900)i,sigma(i),P(i),T(i),Tau(i),Kappa(i),gradad(i),gradrad(i),PRC,Gv1,Gv2,Gv3,Gp,Beta,Ab,Teff,Tint,Teq0,mu,f,grav,ROSS,COEFF,COMP,STAR,CONV,ALBEDO
!         write(40,9700)i,Omega,R_planet,R,cp,grav,Beta,Betav1,Betav2,Betav3,Tirr, &
!         &Pmax,Gv1,Gv2,Gv3,Gamma1,Gamma2,sigma(i),P(i),Kappa(i),T(i),rho(i)&
!         &,lpkpt(i),lpkpr(i),J1(i),J2(i),H1(i),H2(i),TAU(i),mu,Gp&
!         &,Hv0,Tint
!        ENDDO


       close(40)

       open(unit=40,form='formatted',status='UNKNOWN',file=OutputFile)
       write(40,1400) Teff,Tint,Teq0,mu,f,grav,Gv1,Gv2,Gv3,Gp,Beta,Ab,&
     &ROSS,COEFF,COMP,STAR,CONV,ALBEDO
       write(40,*),' '
       write(40,*),' '
       write(40,1500)
       DO i=1,N
         write(40,9800) i,P(i),T(i),Tau(i),Kappa(i),gradad(i),gradrad(i)&
     &,PRC
       ENDDO
      close(40)
!      ENDDO


!      !Open and write output file original
!       open(unit=40,form='formatted',status='UNKNOWN',file=OutputFileCSV&
!     &)
!       write(40,'(A)')'N,P(Pa),T(K),Tau,Kappa(m^2/kg),gradad,gradrad,PRC&
!     &,Gv1,Gv2,Gv3,Gp,Beta,Ab,Teff(K),Tint(K),Teq0(K),mu,f,grav(m/s^2), &
!     &ROSS,COEFF,COMP,STAR,CONV,ALBEDO'
!        DO i=1,N
!         write(40,9900)i,P(i),T(i),Tau(i),Kappa(i),gradad(i),gradrad(i),&
!     &PRC,Gv1,Gv2,Gv3,Gp,Beta,Ab,Teff,Tint,Teq0,mu,f,grav,ROSS,COEFF,   &
!     &COMP,STAR,CONV,ALBEDO
!        ENDDO
!       close(40)
!
!       open(unit=40,form='formatted',status='UNKNOWN',file=OutputFile)
!       write(40,1400) Teff,Tint,Teq0,mu,f,grav,Gv1,Gv2,Gv3,Gp,Beta,Ab,&
!     &ROSS,COEFF,COMP,STAR,CONV,ALBEDO
!       write(40,*),' '
!       write(40,*),' '
!       write(40,1500)
!       DO i=1,N
!         write(40,9800) i,P(i),T(i),Tau(i),Kappa(i),gradad(i),gradrad(i)&
!     &,PRC
!       ENDDO
!      close(40)
!!      ENDDO



!      !EJ Open and write output file
!       open(unit=40,form='formatted',status='UNKNOWN',file=OutputFileCSV&
!     &)
!       write(40,'(A)')'N,sigma,Omega_0,Omega,R_rho,R,c_rho,grav(m/s^2),Beta,Betav1,Betav2,Betav3,A1,A2,Teff,P_H,Gamma1,Gamma2,Gv1,Gv2,Gv3,Kappa_R,T,rho,pkpt,pkpr,J1,J2,H1,H2'
!        DO i=1,N
!         write(40,9900)i,P(i),T(i),Tau(i),Kappa(i),gradad(i),gradrad(i),&
!     &PRC,Gv1,Gv2,Gv3,Gp,Beta,Ab,Teff,Tint,Teq0,mu,f,grav,ROSS,COEFF,   &
!     &COMP,STAR,CONV,ALBEDO
!        ENDDO
!       close(40)
!
!       open(unit=40,form='formatted',status='UNKNOWN',file=OutputFile)
!       write(40,1400) Teff,Tint,Teq0,mu,f,grav,Gv1,Gv2,Gv3,Gp,Beta,Ab,&
!     &ROSS,COEFF,COMP,STAR,CONV,ALBEDO
!       write(40,*),' '
!       write(40,*),' '
!       write(40,1500)
!       DO i=1,N
!         write(40,9800) i,P(i),T(i),Tau(i),Kappa(i),gradad(i),gradrad(i)&
!     &,PRC
!       ENDDO
!      close(40)
!!      ENDDO

!
      write(*,*) Tirr,'Tirr'
      write(*,*) Tmu,'Tmu'
      write(*,*) (Omega**2*R_planet**2/R),'tcoefficient'
      write(*,*) 'Done, PTprofile is available in files: "',            &
     &TRIM(OutputFile),'" or "',TRIM(OutputFileCSV),'"'
      END PROGRAM paper2
