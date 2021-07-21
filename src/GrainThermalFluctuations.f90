! The subroutines for calculations of thermal evolution of solid particles
! upon absorption of UV photons of interstellar radiation field written for publication 
! Zhukovska, Henning and Dobbs, 2018, ApJ, 857, 94. Please cite the paper if you use the code.
! Main subroutine: computation of the cooling timescales for a grain
! with initial temperature T0 and radius agr 
! Steps:
! 1) calculate normal mode spectrum assuming Debye normal modes  
! 2) compute temperatures corresponding to the energy absorbed on a given energy value grid  
! 3) compute radiative cooling timescales for the given energy value grid
! Implemented materials: silicate or iron
! Supplementary subroutines include calculations of the maximum temperature to which a grain
! is heated, mean time between UV absorptions and mean energy of absorbed photons, cooling rates, etc.

!-----------------------------------------------------------------

      Module Constants
        use nrtype
        real, parameter :: clight=2.99792e10,me=9.14823e-28,mp=1.67962e-24,mh=1.66053e-24
        real, parameter :: boltz=1.38062e-16,h=6.62619e-27,guniv=6.67323e-8, boltzevK=8.6173303e-5, & !k_b in eV/K
                                Navogadro = 6.02214086e23
        real, parameter :: echarge=4.8e-10,acons=7.57e-15, stefan=5.67e-5
! astronomical constants (cgs)
        real, parameter :: Msol=1.989e33, lsol=3.826e33, Rsol=7.0e10
! conversion factors
        real, parameter :: cmau=1.49597e13,cmpc=3.0856e18,cmkpc=3.0856e21,cmkm=1.0e5, &
                           everg=1.60184e-12
! grain materials
        real(dp), parameter :: kmrn = 3.5d0, dsil = 4.d-3, dcar = 3.4d-3, & ! dust2h mass ratios for sil and car, Tagach+03
             rhosil = 3.5d0 , rhocar = 2.24d0 !grain material densities  g cm-3, WD2001
        real(dp), parameter :: Ac=12., rhoiron=7.87, Afe=56.         
      End Module Constants

!-----------------------------------------------------------------

      Module Data
        use nrtype
        integer, parameter :: nq=2000, nc=25, nTemp=1500
        real(dp) :: labs(nq), qabs(nq), Tsol(nc), Csol(nc), agr, int2
        integer, parameter ::nsolfine=1000000
        real(dp) :: T0, Ephoton, Tsolfine(nsolfine), csolfine(nsolfine), Temp(nTemp), En(nTemp)

        integer :: nm
        real(dp) :: Eupper
        real(dp), allocatable :: Enm(:)
      End Module Data     

!-----------------------------------------------------------------


  
!*************************************************************************!
! Compute Radiative cooling time in thermal descret approaximation eq 44 Draine&Li2001 (grain as a vybrational system)
! a - grain radius, Eu - vibration energy of the system, Tu - corresponding temperature 
      Function taucooltd(a, Eu, Tu)
        use nrtype
        use Data
        use Constants
        implicit none
        integer ::  i, ic, io, il
        integer, parameter :: nmax=10000
        real(dp), intent(in) :: a, Tu, Eu
        real(dp) :: taucooltd, E(nmax), arr(nmax), int1, lam, Qabsinterpol
        real(dp), external :: TrapInt, Interpol1D

        call FillArrayEquidlog(nmax, 1.d0*h*clight, Eu*h*clight, E)

        do i=1,nmax
           lam = h*clight/E(i)
           if(lam.gt.labs(nq))then
             Qabsinterpol = Qabs(nq) * (labs(nq)/lam)**2
             else
                if(lam.le.labs(1))then
                    Qabsinterpol=Qabs(1)
                   else
                    Qabsinterpol = Interpol1D(nq, labs, Qabs, lam)
                endif   
           endif
           arr(i) = E(i)**3 * Qabsinterpol/(exp(E(i)/boltz/Tu)  - 1.d0) 
        enddo
        int1 = TrapInt(nmax, E, arr)
        taucooltd = 1./( 1./(Eu*h*clight)  * int1*  pi*a**2 * 8*pi/h**3/clight**2 )
                
      End Function taucooltd
!*************************************************************************!
! T is derived from E-T relation for bulk iron
      Subroutine tauRadCoolIron
        use nrtype
        use Constants
        use Data
        implicit none
        integer, parameter :: nmax=500
        integer :: i, nclsize
        real(dp) :: Egrid(nmax), Tgrid(nmax), Tu, TempCl, Tgrainmax, &
                   taucooltd, Interpol1D
        character(200) :: fname

       write(*,'(/, A)') '--- Radiative cooling time calculations for iron grains ---' 
        call ReadDataIron(agr)
        call FillArrayEquidLog(nmax, 30.d0, 2d5, Egrid)
        do i=1,nmax
           Tgrid(i) = TempCl(agr, Egrid(i))
        enddo   

        write(*,'(A, g11.2)') 'Min cooling time (for 91 nm photons) [s]: ', &
             taucooltd(agr,  109690.8d0, Interpol1D(nmax, Egrid, Tgrid, 109690.8d0)) 
        
         call OpenOutFile('iro', fname)

        do i=1,nmax
           write(3,*) Egrid(i), taucooltd(agr, Egrid(i), Tgrid(i))
        enddo   
        close(3)
      end Subroutine tauRadCoolIron

!------------------------------------------------------------
! See description for tauRadCoolsil_harmonosc
      Subroutine tauRadCoolIron_harmonosc
        use nrtype
        use Constants
        use Data
        implicit none
        integer, parameter :: nmax=150
        integer :: i, nclsize, Na
        real(dp) :: Eugrid(nmax), Tugrid(nmax), taucooltd, interpol1d 
        character(200) :: fname

       write(*,'(/,A)') '--- Radiative cooling time calculations for iron grains ---' 
       call ReadDataIron(agr)
       
        Na= nclsize(agr); nm=3*Na-6;
        allocate(Enm(Nm))
        call NormalModesIron(Na, nm, Enm)
        call TemperatureLevelu(Na, nmax, Eugrid, Tugrid)
        deallocate(Enm)
        
        call OpenOutFile('iro', fname)
               
        do i=1,nmax
           write(44,*) Eugrid(i), Tugrid(i)
           write(3,*) Eugrid(i), taucooltd(agr, Eugrid(i), Tugrid(i))
        enddo   
        write(44,*); write(44,*) 
        close(3)

       write(*,'(/,A, A)') 'Cooling times are successfully written to: ', fname(1:len_trim(fname))         
        write(*,'(A, g11.2)') 'Min cooling time (for 91 nm photons) [s]: ', &
             taucooltd(agr,  109690.8d0, Interpol1D(nmax,  Eugrid, Tugrid, 109690.8d0)) 

      End Subroutine tauRadCoolIron_harmonosc      

!-----------------------------------------------------------------
! E1 normal mode spectrum in cm-1 (energies of vibrational states)
      Subroutine NormalModesIron(Na, nm, e1)
        use nrtype
        use Constants
        implicit none
       integer :: i,j,nm, na
       real(dp) :: td1, betan, delt, n
       real(dp) :: e1(Nm)
 
       Td1=415.d0
       n=3.d0
       do j=1,nm
          delt=0.5d0
          if(j.eq.3.or.j.eq.2)delt=1d0
          e1(j) = boltz*Td1/h/clight * ( (1-betan(n,Nm))/Nm * (j-delt) + betan(n,Nm) )**(1./n)
       enddo

     End Subroutine NormalModesIron
!-----------------------------------------------------------------
! Radiative cooling time calculations for silicate grains as a function of absorbed photon energy E
! given on an arbitrary energy grid Eugrid. Corresponding grain temperatures Tugrid are derived
! using expression for the expected energy of harmonic oscilators
! Output: Eugrid, taucool are written in a file
     Subroutine tauRadCoolsil_harmonosc
        use nrtype
        use Constants
        use data, only: Eupper, Enm, nm, agr
        implicit none
        integer, parameter :: nupper=250
        integer :: i,j, nm1,nm2, na, nclsize
        real(dp) ::  Eugrid(nupper), Tugrid(nupper), taucooltd, Interpol1D
        character(2) :: stragr
        character(200) :: fname

       Na= nclsize(agr); Nm = 3*Na-6
       write(*,'(/,A)') '--- Radiative cooling time calculations for silicate grains ---'
       
       call DraineQabsSil(agr)
       
       allocate(Enm(Nm))
       call NormalModesSilicates(Na, nm, Enm)      
       call TemperatureLevelu(Na, nupper, Eugrid, Tugrid)
       deallocate(Enm)
 
       call OpenOutFile('sil', fname)
       
       do i=1,nupper
          write(3, '(1p, 2e11.3)' )  Eugrid(i), taucooltd(agr,  Eugrid(i), Tugrid(i))
       enddo
       close(3)

       write(*,'(/,2A)') 'Cooling times are successfully written to: ', fname(1:len_trim(fname))         
       write(*,'(A, g11.2)') 'Min cooling time (for 91 nm photons) [s]: ', &
            taucooltd(agr,  109690.8d0, Interpol1D(nupper, Eugrid, tugrid, 109690.8d0)) 
        
      End Subroutine tauRadCoolsil_harmonosc

!-----------------------------------------------------------------
! Calculates temperatures of the vibrational system with number of atoms Na by solving
! Eq. 25 from Zhukovska+2018 for the expectation value of vibration energy with Brent method.
! Requires pre-computed normal mode spectrum Enm
! T is set to constant for the first 20 modes 
! Output: arbitrary energy grid Eugrid and corresponding grid of values Tugrid 
      Subroutine TemperatureLevelu(Na, nupper, Eugrid, Tugrid)
       use nrtype
        use Constants
        use data, only: Eupper, Enm, nm, agr
        implicit none
	INTERFACE
		FUNCTION ExpectedEnergyDiff(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: ExpectedEnergyDiff
		END FUNCTION ExpectedEnergyDiff
        END INTERFACE
        integer, intent(in) :: Na, nupper
        real(dp), intent(out) ::  Eugrid(nupper), Tugrid(nupper)
        real(dp) ::  tol, zbrent_mod
        integer :: i
        logical :: interval
        
        print'(/,A, i4)', 'Calculations of excitation temperatures T on an arbitrary energy grid E for Na=', Na
        tol=1d-7

        do i=1,20
            Eugrid(i) = Enm(i)
            Tugrid(i) = Enm(1)*h*clight/boltz/log(2.)
         enddo

        call FillArrayEquidlog(nupper-20, Enm(20)*1.01, 2d5, Eugrid(21:nupper))
! compute Tu for a grain at the vibrational state Eupper
        do i=21, nupper
           Eupper = Eugrid(i)
           Tugrid(i) = zbrent_mod(ExpectedEnergyDiff, 1.d0, 1d4, tol, interval) 
           if(interval) then
              print'("i, E_i,  T_i" , i5, 1pe12.2, 0p, f10.2)', i, Eugrid(i), tugrid(i)
           else
              print'(A)', 'Temperature not found. Check initial grain temperature.'
           endif   
        enddo      
      End Subroutine TemperatureLevelu

!-----------------------------------------------------------------
! ExpectedEnergyDiff = ExpectationEnergy - Eu 
! Benchmarked: the same results as Fig. 8 for silicates
      Function ExpectedEnergyDiff(x) 
        use data, only: Enm, nm, Eupper
        use Constants
        implicit none
        integer :: j
        real(dp) ::  Eexp, ExpectedEnergyDiff
        real(dp), intent(in) ::  x

       Eexp = 0.
         do j=1,nm
           Eexp = Eexp + Enm(j)*h*clight/( exp(Enm(j)*h*clight/x/boltz) -1. )
         enddo
         ExpectedEnergyDiff = (Eupper*h*clight - Eexp)/(Eupper*h*clight+Eexp)
      end function ExpectedEnergyDiff
     
!-----------------------------------------------------------------
! E1 - normal mode spectrum in cm-1 (energies of vibrational states)
! n-dimensional Debye mode spectrum  
      Subroutine NormalModesSilicates(Na, nm, e1)
        use nrtype
        use Constants
        implicit none
        integer, parameter :: nmax= 140000 
        integer :: i,j, nm1,nm2,nm, na
       real(dp) :: e1(nmax),e2(nmax), td1, td2, betan, delt, n
       if(Na*3.gt.nmax)then
          print*,'increase the size of Enm in NormalModes to ', Na*3; stop
       endif   
       Td1=500.
       nm1=2*na-4
       nm=nm1
       n=2.
       do j=1,nm
          delt=0.5
          if(j.eq.3.or.j.eq.2)delt=1
          e1(j) = boltz*Td1/h/clight * ( (1-betan(n,Nm))/Nm * (j-delt) + betan(n,Nm) )**(1./n)
       enddo
       Td2=1500.    
       nm2=na-2
       nm=nm2
       n=3.
       do j=1,nm
          delt=0.5
          if(j.eq.3.or.j.eq.2)delt=1
          e2(j) = boltz*Td2/h/clight * ( (1-betan(n,Nm))/Nm * (j-delt) + betan(n,Nm) )**(1./n)
          e1(j+nm1)=e2(j)
       enddo      

       nm=nm1+nm2
       call SortArrayShell(nm, e1(1:nm))
       return
       do j=1,nm
          write(8,*) j, e1(j)
       enddo   
       write(8,*)
       write(8,*) 
     End Subroutine NormalModesSilicates

     Function Betan(n,Nm)
        use nrtype
        implicit none
       integer :: Nm
       real(dp) :: betan,n
       betan = (Nm**(1.-n/3.) -1 )/( 2. *Nm -1.)
     end Function Betan
  
!------------------------------------------------------------
! a grain radius in cm       
       Function TempCl(a, E)
         use nrtype
         use Data
        implicit none
        integer :: ncl, nclsize, io,il
        real(dp) :: TempCl, a, E, Interpol1D

        ncl = nclsize(a)
        call findinarray(en,E/(ncl-2d0),ntemp,io,il)
        TempCl = Interpol1D(ntemp, En, Temp, E/(ncl-2d0))
       End Function TempCl
 
!------------------------------------------------------------
! agr in cm       
       Function nclsize(agr)
        use nrtype
        implicit none
        integer :: nclsize
        real(dp) :: agr
         nclsize = 352*(agr*1e7)**3
       End Function nclsize
!------------------------------------------------------------

! Calculates the cooling rate of a grain with the initialised size agr    
! Output written in '../out/Cool-grain-*' file
!-----------------------------------------------------------------     
      Subroutine IronCoolingRate (Tmax)
        use nrtype
        use Constants
        use Data
       implicit none
       integer :: i,IO, IL, ic
       real(dp), intent(in) :: Tmax 
       real(dp) :: arr(nq), t, dt, dto, dtn, DELT, & ! all in time units, seconds
                   yo, yn, dydtn, dydto, &
                   dTempoverdt, B_l, TrapInt, JisrfM94, Tgrainmax
       character(2) :: stragr
       character(3) :: mater
       character(4) :: Tmaxstr
       character(200) :: fnameout

!       mater='gra'
        mater='iro'
        if(mater.eq.'iro')then
           call ReadDataIron(agr)
          else
           call ReadDataGra
        endif

        if(agr.lt.10e-7)then
          write(stragr,'(a1,i1)') '0',int((agr+1e-15)*1e7)
         else
          write(stragr,'(i2)') int((agr+1e-15)*1e7)
        endif
       
        write(Tmaxstr, '(i4)' ) int(Tmax)
        Tmaxstr = adjustl(Tmaxstr)

        fnameout = '../out/Cool-grain-'//mater//stragr//'nm-Tmax'//Tmaxstr(1:len_trim(Tmaxstr))//'K.out'
        open(3, file = fnameout(1:len_trim(fnameout)), action='write')
        write(3,'(A)') "# Cooling rate evolution for an iron grain with size "//stragr//"nm"
        write(3,'(4A11)') "#    time  ", "timestep", "T[K]  ", "dTempoverdt"
        
        print'(/,A, i4, A)', "--- Cooling rate evolution for an iron grain with size "//stragr//"nm heated to T_max=", int(Tmax), " K ---"

!     Initial temperature spike, has to be <340K for iron & <100K for gra
        yo = Tmax 
        call FindInArray(labs, 950d-7,nq,ic,il)!Draine&Andersen assume lam_min=1000 mkm for continuous heating between events
        arr = 0d0
        do i=ic,nq
           arr(i) = qabs(i) *  (4d0*pi* ( B_l(2.73d0,labs(i)) + 1d-14*B_l (7500d0,labs(i)) + 1d-13*B_l(4000d0,labs(i)) + &
              4d-14*B_l(3000d0,labs(i)) ) ) !== JisrfM94(labs(i))*clight/labs(i)**2
        enddo
!     Heating by ISRF     
        int2 = TrapInt(nq, labs, arr)
        
!     Initial time
        t=0D0
!     Intial time step
        dt=yo/abs(dTempoverdt(yo))/100.
        dto=dt
!     Limitation of time step width
        DELT=1d0
      
        call ADAMBASH(T,DT,DTo,DTn,yo,yn,DYDTN, DYDTo, DELT)
        yo=yn
        print'(A, 1p,5e11.3)', 'Initial values for integration time, timestep, T, dTdt:', t, dt, yn, DYDTN
        DYDTo = DYDTn
        t = t + dt


! integration cycle:        
10      call ADAMBASH(T,DT,DTo,DTn,yo,yn,DYDTN, DYDTo, DELT)

        IF(DTN.GT.DELT) then
         DT=DELT
          else
         DT=DTN
      EndIF
      yo = yn
      DYDTo = DYDTN
      t = t + dt
      ! print'(1p,4e11.3)', t, dt, yn, DYDTN
      write(3,'(1p,4e11.3)')  t, dt, yn, dTempoverdt(yn)
      
      if(t.lt.15*3600..and.yo.gt.20.)goto 10
      
      close(3)
      
      print'(A, 1p,5e11.3)', 'Final values from integration time, timestep, T, dTdt:', t, dt, yn, dTempoverdt(yn)      
      print'(A,A)', 'Output is successfully written to ', fnameout(1:len_trim(fnameout))
      
    end Subroutine IronCoolingRate

!------------------------------------------------------------
! Derivative used to integrate the temperature evolution of stochastically heated grains
      Function dTempoverdt(T)
        use nrtype
        use Data
        integer :: i, ic, io, il
        real(dp) :: int1, T, dTempoverdt, arr(nq)
        real(dp), external :: TrapInt, B_l, Interpol1D, JisrfM94, B_nu

        do i=1,nq
          arr(i) = B_l(T, labs(i))*Qabs(i)
        enddo
        int1 = 4d0*pi*TrapInt(nq, labs, arr)

        dTempoverdt = 3./4./agr/Interpol1D(nc,Tsol,csol,T) * (int2 - int1)
        
      End Function dTempoverdt
!------------------------------------------------------------

!************************************************************
            
      Subroutine ReadDataIron(agr0)
        use nrtype
        use Data
        use Constants 
        implicit none
        integer :: i
        real(dp) :: tmp, natoms, agr0
        open(1, file='../data/Qabs-iron-Fischera04.dat', action='read')
        open(2, file='../data/Heat-Capac-iron-Dasai86.dat', action='read')
        open(3, file='../data/E-T_HD2017.dat', action='read')
        
        read(1,*);  read(2,*);  read(3,*)
        do i=1,nq
           if(agr0.le.2d-7)then
              read(1,*) labs(i), qabs(i)
           else
              if(agr0.ge.3d-7.and.agr0.lt.10d-7)then
              read(1,*) labs(i), tmp, qabs(i)
                else
                   if(agr0.eq.10d-7)then
                      read(1,*) labs(i), tmp, tmp, qabs(i)
                    else
                    print*, 'no Qabs data for iron grain with a=',agr0; stop
                 endif
              endif
           endif           
        enddo
        do i=1,nc
           read(2,*) Tsol(i), Csol(i)
        enddo
        do i=1,nTemp
           read(3,*) Temp(i), En(i)
        enddo
        close(1); close(2); close(3)
        En = En/clight/h
        labs = labs*1e-4 !mkm2cm
      End Subroutine ReadDataIron

!------------------------------------------------------------

      Subroutine ReadDataGra
        use Data
        use Constants 
        implicit none
        integer :: i
        real(16) :: tmp, natoms
        open(1, file='Qabs-gra-Draine.dat', action='read')
        open(2, file='Heat-Capac-gra-Horn07.dat', action='read')
        
        read(1,*)

        do i=1,nq
           if(agr.eq.1e-7)then
              read(1,*)labs(i), qabs(i)
           else
             read(1,*)labs(i), tmp, qabs(i)
           endif  
        enddo

        read(2,*)
        do i=1,nc
           read(2,*)Tsol(i), Csol(i)
        enddo
        close(1); close(2)

        labs = labs*1e-4 !mkm2cm
        
      End Subroutine ReadDataGra
!------------------------------------------------------------
! NOT USED 
      Subroutine ReadDataSil
        use Data
        use Constants 
        implicit none
        integer :: i
        real(16) :: tmp, natoms
        open(1, file='Qabs-sil-Draine.dat', action='read')       
        read(1,*)
        do i=1,nq
           if(agr.eq.1e-7)then
              read(1,*) labs(i), qabs(i)
           else
             read(1,*)labs(i), tmp, qabs(i)
           endif  
        enddo
         close(1); 
        labs = labs*1e-4 !mkm2cm
      End Subroutine ReadDataSil
!------------------------------------------------------------
!
! PLANCK FUNCTION B_lambda 
!    erg/cm/cm2/s?
! lam in cm
!------------------------------------------------------------
    Function B_l(T, lam)
      USE nrtype; 
      implicit none
      real(DP), parameter ::  clight = 2.99792d10, h = 6.62619d-27, stefan = 5.67d-5, kB = 1.38d-16
      real(DP), intent(in) :: T, lam
      real(DP) :: B_l, x

      x = h*clight/(lam*kB*T)
      B_l = 2d0*h*clight**2/lam**5/(dexp(x)-1d0) 
    End function B_l
!------------------------------------------------------------
! PLANCK FUNCTION
! B_nu (= B_lambda * c/nu^2)    B_lam = Bnu * nu^2/c= Bnu/lam^2*c
! lam in mkm
!------------------------------------------------------------
      Function B_nu(T, lam)
         USE nrtype; 
      implicit none
      real(DP) :: B_nu, T, lam, x, nu
      real(DP), parameter ::  clight = 2.99792d10, h = 6.62619d-27, &
                              stefan = 5.67d-5, kB = 1.38d-16
      nu = clight/(lam )
      x = h*nu/(kB*T)
      if(x.gt.300d0) then ! Wien law
        B_nu = 2d0*h*nu**3/clight**2 * dexp(-x)
      else
        B_nu = 2d0*h*nu**3/clight**2/(dexp(x)-1d0)
      endif
    End Function B_nu
!------------------------------------------------------------
! Mean intensity of interstellar radiation field
! Mathis et al 1994
!    erg/Hz/cm2/s     J_nu/c == u_nu
!    [lam]=cm
!------------------------------------------------------------
     Function JisrfM94(lam)
        USE nrtype; 
	implicit none
        real(dp), intent(in) :: lam
        real(DP), parameter ::  clight = 2.99792d10
        real(dp) :: JisrfM94, B_nu, B_l
! Mathis et al 1994
! Stellar componet is dominated by BB at  T=4000 K with corresponding max (for B_nu) lam_max = 1.2 mkm  
!        Jisrf = 4d0*pi(0d0 + 0d0*B_l(2.73d0,lam) + 1d-14*B_nu (7500d0,lam) + 1d-13*B_nu(4000d0,lam) + &
!             4d-14*B_nu(3000d0,lam) ) 

        JisrfM94 = 4d0*pi* ( B_nu(2.73d0,lam) + 1d-14*B_nu (7500d0,lam) + 1d-13*B_nu(4000d0,lam) + &
              4d-14*B_nu(3000d0,lam) ) 
      End Function JisrfM94

!------------------------------------------------------------

    Function natoms(agr,rho,A)
      use nrtype
        use Constants 
         implicit none
         real(dp) :: natoms,agr,rho,A
         natoms = 4./3.*pi*rho*agr**3/A/mp
      End Function natoms

!------------------------------------

      Function Tgrainmax(lam, a0, T00)
        use nrtype
        use Data
        use Constants
        implicit none
        real(dp),intent(in) :: lam, t00, a0
        real(dp) :: Tgrainmax, zbrent_mod
        real(dp) :: tol,Trapint
        logical :: interval
        integer :: i,io0,il0,io,il
	INTERFACE
		FUNCTION fdiff(x)
		USE nrtype
		IMPLICIT NONE
		REAL(DP), INTENT(IN) :: x
		REAL(DP) :: fdiff
		END FUNCTION fdiff
        END INTERFACE

        Ephoton = h*clight/lam
        tol = 1d-10
        T0=T00
        call FillArrayEquid(nsolfine,Tsol(1), Tsol(nc), Tsolfine)
        call NewGridArrayExtra(nc,nsolfine,Tsol,Tsolfine,Csol,Csolfine)    

        Csolfine = Csolfine *4./3.*pi*a0**3
       
        Tgrainmax = zbrent_mod(fdiff,T0,Tsol(nc),tol,interval)

      End Function Tgrainmax
!------------------------------------------------------------
     
      Function fdiff(x)
        use nrtype
        use Data
        implicit none
        real(dp), intent(in) :: x
        real(dp) :: fdiff, TrapInt
        integer :: i,io0,il0, io,il

        call findinarray(tsolfine,T0,nsolfine,io0,il0)
        call findinarray(tsolfine,x,nsolfine,io,il)
        
	TrapInt=0d0
	do i = io0, il
	  TrapInt = TrapInt + (Csolfine(i)+Csolfine(i-1))*0.5d0 * (tsolfine(i)-tsolfine(i-1))
        enddo

        fdiff = (Ephoton - TrapInt)/(Ephoton+TrapInt)
       
      End Function fdiff

!***********************************************************************C
      Subroutine TgrainmaxFileOut
        use nrtype
        use Constants
        use Data
        implicit none
        integer, parameter :: nmax=100
       integer :: i, nclsize
       real(dp) :: Egrid(nmax), taucooltd, Tu, TempCl, Tgrainmax
       character(2) :: stragr

        call ReadDataIron(agr)
        call FillArrayEquidLog(nmax, 30.d0, 1d5, Egrid)
        
        do i=1,nmax
           write(*,*)'i=',i
           !initial Temperatures of grains are taken from Fig.7 Fischera
!           write(44,*)   Egrid(i), Tgrainmax(1d0/Egrid(i),1d-7,67d0), &
!               Tgrainmax(1d0/Egrid(i),3d-7,70d0), Tgrainmax(1d0/Egrid(i),10d-7,75d0)
           write(44,*)   Egrid(i), Tgrainmax(1d0/Egrid(i),1d-7,8d0), &
               Tgrainmax(1d0/Egrid(i),3d-7,40d0), Tgrainmax(1d0/Egrid(i),10d-7,75d0)
        enddo
      End Subroutine TgrainmaxFileOut
!***********************************************************************C
!      Stepper for Differential equations. Method: Adams-Bashforth      C
!***********************************************************************C
     SUBROUTINE ADAMBASH(T,DT,DTo,DTn,yo,yn,DYDTN, DYDTo, DELT)
       use nrtype
       implicit none
       integer :: n
       real(dp):: x, dTempoverdt, &
                  T,DT,DTo,DTn,yo,yn,DYDTN, DYDTo, DELT

        DYDTN = dTempoverdt(Yo)
        
        YN = YO + (DYDTN + 0.5D0*(DYDTN - DYDTO)*DT/DTO) * DT

!     New timestep width
      DTN=1D30

      IF(YN.GT.1D-15) THEN
          X=0.02D0*DABS(YN)/(DABS(DYDTN)+1D-30)
        ELSE
          X=1D30
        END IF
        IF(X.LT.DTN) THEN
          DTN=X
        END IF

      Dto=DT
!     New timestep width
      IF(DTN.GT.1.2D0*DT) DT=1.2D0*DT
      IF(DT.GT.DELT) DT=DELT
      IF(DTN.LT.DT) DT=DTN
      END SUBROUTINE ADAMBASH
!*************************************************************************!

!*************************************************************************!

      Subroutine DraineQabsSil(agr0)
        use nrtype
        use data, only: labs, qabs
        implicit none
        integer, parameter :: n=241, nn=2000
        real(dp) :: l1(n), q1(n),l2(n), q2(n), ln(nn), q1n(nn), q2n(nn), amin, amax, agrgrid(81),agr0
        integer :: i,IO, IL, nsize, isize
        
        open(1, file='../data/Sil_81', action='read')
        open(2, file='../data/Qabs-sil-Draine.dat')

        call ShiftReadFile (1,3)
        read(1,*)nsize, amin, amax
        call ShiftReadFile (1,4)
        call FillArrayEquidLog(nsize, amin*1d-4, amax*1d-4, agrgrid)
        call FindInArray(agrgrid, agr0, nsize, IO, IL)
        if(agrgrid(io)-agr0.lt.agr0-agrgrid(il))then
           isize=io
        else
           isize=il
        endif   
!        print'(1p4e11.2)',agrgrid(io),agr0,agrgrid(il), agrgrid(isize)

        if(isize.gt.1) call ShiftReadFile (1,(n+3)*(isize-1))
        
        do i=1,n
           read(1,*)l1(i), q1(i)
        enddo
        call  Sort2ArrayShell(n, l1, q1)
        call FillArrayEquidLog(nn,1d-3, 1d3, ln)      
        ln=log10(ln); l1=log10(l1); q1=log10(q1); 
        call NewGridArrayExtra(n,nn,l1,ln,q1,q1n)      
        ln=10**ln; q1n=10**q1n; l1=10**l1; q1=10**q1
        
        write(2,*)'# lamda[mkm], Qabs'
        write(2,*)'# a(cm)=',agrgrid(isize)
        q1n(1) = q1(1)
        do i=1,nn
           write(2,*) ln(i), q1n(i)
        enddo
        close(1)
        close(2)
        labs=ln
        qabs=q1n
        labs = labs*1e-4 !mkm2cm
      End Subroutine DraineQabsSil
!*************************************************************************!
! Time between photon absorptions of UV photons
      Function tauabs(a)
        use nrtype
        use Constants, only: h
        use data, only: labs, qabs, nq
        implicit none
        integer, parameter :: n=100000
        real(dp) :: tauabs, a, tmp(nq), JisrfM94, TrapInt, lam(n), tmpn(n)
        integer :: i

        do i=1,nq
           tmp(i)=Qabs(i)*JisrfM94(labs(i))/labs(i)
           if(labs(i).gt.400d-7.or.labs(i).lt.91.16d-7) tmp(i)=0.d0
        enddo
          
        call FillArrayEquid(n,labs(1), labs(nq), lam)
        call NewGridArrayExtra(nq,n,labs,lam,tmp,tmpn)
        
!        tauabs = 1d0/TrapInt(nq,labs,tmp)/pi/a**2*h
        tauabs = 1d0/TrapInt(n,lam,tmpn)/pi/a**2*h
      End Function tauabs
!*************************************************************************!
      Function meanabsenergy(a)
        use nrtype
        use Constants, only: h, clight,everg
        use data, only: labs, qabs, nq
        implicit none
        real(dp) :: tauabs, a, tmp(nq), JisrfM94, TrapInt, meanabsenergy
        integer :: i

        do i=1,nq
           tmp(i)=Qabs(i)*JisrfM94(labs(i)) /labs(i)**2
           if(labs(i).gt.400d-7.or.labs(i).lt.91.16d-7) tmp(i)=0.d0
        enddo
        meanabsenergy = tauabs(a) * TrapInt(nq,labs,tmp)*pi*a**2*clight
      End Function meanabsenergy
!*************************************************************************!
      Subroutine OpenOutFile(nam, fname)
        use Data, only: agr
        character(2) :: stragr
        character(3) :: nam
        character(200), intent(out) :: fname
       
        if(agr.lt.1d-7)then
          write(stragr,'(a1,i1)') '0',int((agr+1e-15)*1e8)
         else
          write(stragr,'(i2)') int((agr+1e-15)*1e8)
       endif
       fname = '../out/tauradcool-'//nam//stragr//'A.out'
        open(3, file=fname(1:len_trim(fname)), action='write')
        write(3, '(2A11)') "#   Eu     ", "tau_cool[s]"

      End Subroutine OpenOutFile
              
