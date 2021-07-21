! A sample program to use subroutines for thermal evolution of iron and silicate grains
     Program Main
       use nrtype
       use Data
       use Constants
       implicit none
       real(dp) :: tauabs, meanabsenergy, Tgrainmax, Tmax, t, cc = 1.d0*3600.d0*24.d0
       integer :: nclsize

       t0 = 30d0
       print'("Computations for initial grain T= ",i4,"K")', int(t0)
      
       print'("tau_abs - mean time between absorption of two UV interstellar photons")'
       print'("E_abs - mean energy of absorbed UV photons")'
       print'("lam_abs - wavelength corresponding to the mean energy")'

! grain radius in cm
       agr = 1.d-7 

       print'("Computations for grain size", f4.1,"nm, number of atoms Na ", i6)', agr*1e7, nclsize(agr)
       print'("          tau_abs[days] lam_abs[nm] E_abs[eV]")'

       call DraineQabsSil(agr)
       print'(a,1p3e11.2)', 'silicates:', tauabs(agr)/cc,  clight*h/meanabsenergy(agr)*1d7,  meanabsenergy(agr)/everg
       call ReadDataIron(agr)
       print'(a,1p3e11.2)', 'iron     :', tauabs(agr)/cc,  clight*h/meanabsenergy(agr)*1d7,  meanabsenergy(agr)/everg

       t = Tgrainmax(91.16d-7,agr, t0)
       print'("Max grain temperature (for 90nm photons) T=",i4,"K")', int(t)
       
       call tauRadCoolSil_harmonosc
       call tauRadCoolIron_harmonosc

! Calculate cooling rates for a grain heated to the temperature Tmax
       Tmax = 1d2
       call IronCoolingRate(Tmax)

      End Program Main
