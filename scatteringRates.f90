subroutine scatteringRates
    use globalVariables
    use scatteringVariables
    implicit none

    real :: g3dAcoustic
    real :: Ld
    real :: gamma

    ! GaAs
    effm = 0.067*m0 ![kg] for gamma, L, X respectively
    rho = 5.36e3 ![kg/m^3]
    vs = 5.24e3 ![m/s]

    ! Polar Optical Phonon Scattering
    w0 = E0/hbar ! [1/s]
    N0 = (exp(E0/(kb*T))-1)**(-1)

    allocate(Energy(nE), k(nE))

    ! Initialize Energy
    do i = 1, nE
        Energy(i) = (2.0/nE)*i
        k(i) = sqrt(2.0*effm*Energy(i)/(hbar*hbarJ))
    enddo

    allocate(GammaMAcoustic(nE),GammaMIonImp(nE), GammaMPop(nE), GammaTot(nE))

    open(unit=10, file='Data/Energy', status="unknown")
    open(unit=11, file='Data/GammaMAcoustic', status="unknown")
    open(unit=12, file='Data/GammaMPop', status="unknown")
    open(unit=13, file='Data/GammaMIonImp', status="unknown")
    open(unit=14, file='Data/GammaTot', status="unknown")

    g3dFac = 1.0/(2.0*(pi**2.0))*(2.0*effm/(hbar**2.0))**(1.5)
    MAcFac = pi/((hbar)**(0.5)*(hbarJ)**(0.5))*(Dac**2.0)*kb*T/(2.0*rho*(vs**2.0))
    MPopFac = ((e**2.0)*w0*(epr0/eprInf-1))/(4.0*pi*epr0*ep0*((hbar)**(0.25)*(hbarJ)**(0.25)))
    MIonFac = (e**2.0)/(eprInf**2*ep0**2) * &
    (hbar/hbarJ)**1.5*(NI*e**2.0) / &
    (16.0*sqrt(2.0*effm)*pi)
    print *, g3dFac
    print *, MAcFac
    print *, MPopFac
    print *, MIonFac

    do i = 1, nE

    ! Density of States
        g3dAcoustic = g3dFac*sqrt(Energy(i)-Ec)

    ! Acoustic Phonon Scattering
        GammaMAcoustic(i) = MAcFac*g3dAcoustic

    ! Polar Optical Phonon Scattering
        GammaMPop(i) = MPopFac * (1.0/(hbarJ*sqrt(Energy(i)/effm*2.0))) * &
        (N0*sqrt(E0/Energy(i) + 1.0) + &
        (N0+1.0)*sqrt(-E0/Energy(i) + 1.0) - &
        E0*N0/Energy(i)*asinh(sqrt(Energy(i)/E0)) + &
        E0*(N0+1.0)/Energy(i)*asinh(sqrt(Energy(i)/E0 - 1.0)))

    ! Ionized Impurity Scattering
        Ld=sqrt(ep0*eprInf*kbJ*T/(e**2.0*NI)) ! [m]
        gamma=sqrt(8.0*effm/hbarJ*((Energy(i)*(Ld**2.0))/(hbar)));
        GammaMIonImp(i)= MIonFac / ((log(1.0+gamma**2.0)-(gamma**2.0)/(1.0+gamma**2.0))*Energy(i)**(-1.5))
    ! Total Scattering
        GammaTot(i) = GammaMAcoustic(i) + GammaMPop(i) + GammaMIonImp(i)
    
    ! Write Scattering Rates 
        write(10, *) Energy(i)
        write(11, *) GammaMAcoustic(i)
        write(12, *) GammaMPop(i)
        write(13, *) GammaMIonImp(i)
        write(14, *) GammaTot(i)

    enddo

    close(10)
    close(11)
    close(12)
    close(13)
    close(14)
        
    ! Numerical Integration
    
    !ExpTauEAcoustic = 0;
   ! ExpEAcoustic = 0;   
    
    !ExpTauEIonImp = 0;
    !ExpEIonImp = 0; 
    
    !ExpTauEPop = 0;
    !ExpEPop = 0; 
    
    !ExpTauETot = 0;
    !ExpETot = 0;
    
    !for i=1:nE
    !    f0=exp(-Energy(i)/(kb*T));
    !    g3d = sqrt(Energy(i));
        
    !    ExpTauEAcoustic = ExpTauEAcoustic + 1/GammaMAcoustic(i)*Energy(i)*f0*g3d;
     !   ExpEAcoustic    = ExpEAcoustic + Energy(i)*f0*g3d;
        
    !    ExpTauEIonImp = ExpTauEIonImp + 1/GammaMIonImp(i)*Energy(i)*f0*g3d;
    !    ExpEIonImp    = ExpEIonImp + Energy(i)*f0*g3d;
    !    
    !    ExpTauEPop = ExpTauEPop + 1/GammaMPop(i)*Energy(i)*f0*g3d;
    !    ExpEPop    = ExpEPop + Energy(i)*f0*g3d;
        
    !    ExpTauETot = ExpTauETot + 1/GammaTot(i)*Energy(i)*f0*g3d;
    !    ExpETot    = ExpETot + Energy(i)*f0*g3d;
        
    !end
    
    ! Find Tau
    !AvgTauAcoustic  = ExpTauEAcoustic/ExpEAcoustic;
    !AvgTauIonImp  = ExpTauEIonImp/ExpEIonImp;
    !AvgTauPop  = ExpTauEPop/ExpEPop;
    !AvgTauTot  = ExpTauETot/ExpETot;
    
    ! Calculate Mobilities
    !MobilityAcoustic = 100^2*e*AvgTauAcoustic/effm; % cm^2/(Vs)
    !MobilityIonImp = 100^2*e*AvgTauIonImp/effm; % cm^2/(Vs)
    !MobilityPop = 100^2*e*AvgTauPop/effm; % cm^2/(Vs)
    
    ! Calculate Total Mobility
    !MobilityTot = 100^2*e*AvgTauTot/effm; % cm^2/(Vs)
    
    ! Matthiessen's Rule
    !Mobility = (1/MobilityAcoustic + 1/MobilityIonImp + 1/MobilityPop).^(-1); % cm^2/(Vs)
    
    !disp(['The acoustic mobility is ' num2str(MobilityAcoustic)])
    !disp(['The ionized impurity mobility is ' num2str(MobilityIonImp)])
    !disp(['The POP mobility is ' num2str(MobilityPop)])
    !disp(['The mobility from Matthiessens Rule is ' num2str(Mobility)])
    !disp(['The mobility from total momentum relaxation rate is ' num2str(MobilityTot)])
    
    ! Conductivities
    !conductivity = NI*e*Mobility;
    !conductivityTot = NI*e*MobilityTot;
    
    !disp(['Conductivity is ' num2str(conductivity)])
    !disp(['Conductivity is ' num2str(conductivityTot)])


end subroutine scatteringRates