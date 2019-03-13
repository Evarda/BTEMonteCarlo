! Calculates Scattering Rates for Each Scattering Mechanism
subroutine scatteringRates
    use GaAsConstants
    use scatteringVariables
    implicit none

    ! GaAs for [Gamma, L, X] Valleys
    real :: rho
    real :: vs
    real :: epr0 = 12.90
    real :: eprInf = 10.92

    ! Useful Values
    real :: g3dAcoustic
    real :: Ld
    real :: gamma

    ! Factors
    real :: g3dFac
    real :: AcFac
    real :: PopFac
    real :: IonFac
    real, dimension(8) :: g3DIVFac

    ! Acoustic Phonon Scattering
    real, dimension(3) :: Dac = (/ 7.01, 9.2, 9.0 /) ! [eV] for gamma, L, X respectively
    
    ! Ionized Impurity Scattering
    real :: NI = 1e23 ! [1/m^3] (Equiv to 10^17 1/cm^3)
    real :: Z = 1

    ! Polar Optical Phonon Scattering
    real :: E0 = 3.536e-2 ![eV]
    real :: w0
    real :: N0

    ! Intervalley Scattering
    real, dimension(8) :: Div   ! [eV/m]
    real, dimension(8) :: Eiv     ! [eV]
    real, dimension(8) :: deltaE                                              ! [eV]
    real, dimension(8) :: effmIV                                              ! [kg]
    real, dimension(8) :: ValleyN                                             ! Number of Valleys
    real :: Niv
    real, dimension(8) :: g3dIVAbs
    real, dimension(8) :: g3dIVEmi

    ! Ordering in Previous Code
    ! GL, GX, LL, LX, XX, (Book Order) 
    ! LG, XG, LL, XL, XX  (Flipped Book Order)

    ! Book: 10, 10, 10, 5, 7 (*10^8 eV/cm)
    real :: DivGL = 10.0e10 ! [eV/m] (1)
    real :: DivGX = 10.0e10 ! [eV/m] (2)
    real :: DivLG = 10.0e10 ! [eV/m] (1 flip)
    real :: DivLL = 10.0e10 ! [eV/m] (3)
    real :: DivLX =  5.0e10 ! [eV/m] (4)
    real :: DivXG = 10.0e10 ! [eV/m] (2 flip)
    real :: DivXL =  5.0e10 ! [eV/m] (4 flip)
    real :: DivXX =  7.0e10 ! [eV/m] (5)

    ! Book: 0.0278, 0.0299, 0.0290, 0.0293, 0.0299
    real :: EivGL = 2.78e-2 ! [eV] (1)
    real :: EivGX = 2.99e-2 ! [eV] (2)
    real :: EivLG = 2.78e-2 ! [eV] (1 flip)
    real :: EivLL = 2.90e-2 ! [eV] (3)
    real :: EivLX = 2.93e-2 ! [eV] (4)
    real :: EivXG = 2.99e-2 ! [eV] (2 flip)
    real :: EivXL = 2.93e-2 ! [eV] (4 flip)
    real :: EivXX = 2.99e-2 ! [eV] (5)

    ! Book: 0.29, 0.49, 0, 0.19, 0
    real :: deltaEGL =  0.29 ! [eV] (1)
    real :: deltaEGX =  0.48 ! [eV] (2)
    real :: deltaELG = -0.29 ! [eV] (1 flip)
    real :: deltaELL =  0.0  ! [eV] (3)
    real :: deltaELX =  0.19 ! [eV] (4)
    real :: deltaEXG = -0.48 ! [eV] (2 flip)
    real :: deltaEXL = -0.19 ! [eV] (4 flip)
    real :: deltaEXX =  0.0  ! [eV] (5)

    ! Counter
    integer :: i        ! Energy
    integer :: valley   ! Valley
    integer :: ivstep   ! Intervalley Pair

    ! GaAs Constants
    effm = (/ 0.067, 0.22, 0.58 /)*m0
    rho = 5.36e3 ![kg/m^3]
    vs = 5.24e3 ![m/s]

    ! Ionized Impurity Scattering
    Ld=sqrt(ep0*eprInf*kbJ*T/(q**2.0*NI)) ! [m]

    ! Polar Optical Phonon Scattering
    w0 = E0/hbar ! [1/s]
    N0 = (exp(E0/(kb*T))-1)**(-1)

    ! Ordering in Previous Code
    ! GL, GX, LL, LX, XX, 
    ! LG, XG, LL, XL, XX

    ! Ordering
    ! GL GX LG LL LX XG XL XX
    ! G has 1 equiv valley, L has 4, X has 3
    Div = (/ DivGL, DivGX, DivLG, DivLL, DivLX, DivXG, DivXL, DivXX /)
    Eiv = (/ EivGL, EivGX, EivLG, EivLL, EivLX, EivXG, EivXL, EivXX /)
    deltaE = (/ deltaEGL, deltaEGX, deltaELG, deltaELL, deltaELX, deltaEXG, deltaEXL, deltaEXX /)
    ! Final Valleys
    ! L X G LL X G L XX
    effmIV = (/ effm(2), effm(3), effm(1), effm(2), effm(3), effm(1), effm(2), effm(3)/)
    ValleyN = (/ 4, 3, 1, 3, 3, 1, 4, 2 /)

    print*, w0, N0

    call allocateEnergy()

    ! Initialize Energy and Wavevector

    open(unit=10, file='Data/Energy', status="unknown")

    do i = 1, nE
        Energy(i) = (2.0/nE)*i
        write(10, *) Energy(i)
        do valley = 1, 3
            k(valley, i) = sqrt(2.0*effm(valley)/hbarJ*Energy(i)/hbar) 
        enddo
    enddo

    close(10)

    allocate(GammaAcousticAbs(3,nE))
    allocate(GammaAcousticEmi(3,nE))
    allocate(GammaIonImp(3,nE))
    allocate(GammaPopAbs(3,nE))
    allocate(GammaPopEmi(3,nE))
    allocate(GammaIVAbsCalc(nE, 8))
    allocate(GammaIVEmiCalc(nE, 8))
    allocate(GammaIVAbs(3, 3, nE))
    allocate(GammaIVEmi(3, 3, nE))
    
    
    ! Write GammaAcoustic
    open(unit=21, file='Data/ScatRates/gamma/GammaAcousticAbs', status="unknown")
    open(unit=22, file='Data/ScatRates/L/GammaAcousticAbs',     status="unknown")
    open(unit=23, file='Data/ScatRates/X/GammaAcousticAbs',     status="unknown")
    open(unit=24, file='Data/ScatRates/gamma/GammaAcousticEmi', status="unknown")
    open(unit=25, file='Data/ScatRates/L/GammaAcousticEmi',     status="unknown")
    open(unit=26, file='Data/ScatRates/X/GammaAcousticEmi',     status="unknown")

    ! Write GammaIonImp
    open(unit=31, file='Data/ScatRates/gamma/GammaIonImp', status="unknown")
    open(unit=32, file='Data/ScatRates/L/GammaIonImp',     status="unknown")
    open(unit=33, file='Data/ScatRates/X/GammaIonImp',     status="unknown")

    ! Write GammaPop
    open(unit=41, file='Data/ScatRates/gamma/GammaPopAbs', status="unknown")
    open(unit=42, file='Data/ScatRates/L/GammaPopAbs',     status="unknown")
    open(unit=43, file='Data/ScatRates/X/GammaPopAbs',     status="unknown")
    open(unit=44, file='Data/ScatRates/gamma/GammaPopEmi', status="unknown")
    open(unit=45, file='Data/ScatRates/L/GammaPopEmi',     status="unknown")
    open(unit=46, file='Data/ScatRates/X/GammaPopEmi',     status="unknown")

    ! Write GammaIV
    open(unit=51, file='Data/ScatRates/gamma/GammaIVLAbs', status="unknown")
    open(unit=52, file='Data/ScatRates/gamma/GammaIVXAbs', status="unknown")
    open(unit=53, file='Data/ScatRates/L/GammaIVGAbs', status="unknown")
    open(unit=54, file='Data/ScatRates/L/GammaIVLAbs', status="unknown")
    open(unit=55, file='Data/ScatRates/L/GammaIVXAbs', status="unknown")
    open(unit=56, file='Data/ScatRates/X/GammaIVGAbs', status="unknown")
    open(unit=57, file='Data/ScatRates/X/GammaIVLAbs', status="unknown")
    open(unit=58, file='Data/ScatRates/X/GammaIVXAbs', status="unknown")

    open(unit=61, file='Data/ScatRates/gamma/GammaIVLEmi', status="unknown")
    open(unit=62, file='Data/ScatRates/gamma/GammaIVXEmi', status="unknown")
    open(unit=63, file='Data/ScatRates/L/GammaIVGEmi', status="unknown")
    open(unit=64, file='Data/ScatRates/L/GammaIVLEmi', status="unknown")
    open(unit=65, file='Data/ScatRates/L/GammaIVXEmi', status="unknown")
    open(unit=66, file='Data/ScatRates/X/GammaIVGEmi', status="unknown")
    open(unit=67, file='Data/ScatRates/X/GammaIVLEmi', status="unknown")
    open(unit=68, file='Data/ScatRates/X/GammaIVXEmi', status="unknown")

    do valley = 1, 3
    g3dFac = 1.0/(2.0*(pi**2.0))*(2.0*effm(valley)/(hbar**2.0))**(1.5)
    AcFac = pi/(2.0*(hbar)**(0.5)*(hbarJ)**(0.5))*(Dac(valley)**2.0)*kb*T/(rho*(vs**2.0))
    IonFac = (Z**2.0)/(pi*eprInf**2.0)*(q**2.0/hbarJ)*(NI*q)*(effm(valley)/hbarJ)*((Ld**4.0)/hbarJ)*(q/ep0)*(1/ep0)
    PopFac = sqrt(2.0)/(8.0*pi)*(1.0/eprInf-1.0/epr0)*(q/sqrt(hbarJ))*(q/ep0)*(w0*sqrt(hbar))*(sqrt(effm(valley))/hbarJ)

    !sqrt(2.0)/(pi**2.0)*effmIV(n,m)**(0.5)/hbar*effmIV(n,m)/hbar/(hbarJ**(0.5))/(hbar**(0.5))
    
    do i = 1, nE

    ! Density of States
        g3dAcoustic = g3dFac*sqrt(Energy(i)-Ec)

    ! Acoustic Phonon Scattering (Reduced by AcFac)
        GammaAcousticAbs(valley, i) = g3dAcoustic
        GammaAcousticEmi(valley, i) = g3dAcoustic

    ! Ionized Impurity Scattering
        gamma=sqrt(8.0*effm(valley)/hbarJ*((Energy(i)*(Ld**2.0))/hbar))
        GammaIonImp(valley, i) = k(valley,i)/(4*(k(valley,i)**2)*(Ld**2)+1)

    ! Polar Optical Phonon Scattering
        GammaPopAbs(valley, i) = 1 / sqrt(Energy(i)) * N0 * &
            log(abs((1.0+sqrt(1.0+(E0/Energy(i))))/(-1.0+sqrt(1.0+(E0/Energy(i))))))
        
        !print *, (Energy(i).gt.E0)
        if (Energy(i).gt.E0) then
            GammaPopEmi(valley, i) = 1 / sqrt(Energy(i)) * (N0+1) * &
                log(abs((1.0+sqrt(1.0-(E0/Energy(i))))/(1.0-sqrt(1.0-(E0/Energy(i))))))
        else
            GammaPopEmi(valley, i) = 0
        endif
    
    ! Intervalley Scattering
        do ivstep = 1, 8
            ! Density of States Factor
            g3dIVFac(ivstep) = sqrt(2.0)/(pi**2.0)*(effmIV(ivstep)**(0.5))/hbar*effmIV(ivstep)/hbar/(hbarJ**(0.5))/(hbar**(0.5))

            ! Density of States without Factor
            g3dIVAbs(ivstep) = sqrt(Energy(i)-deltaE(ivstep)+Eiv(ivstep)-Ec)
            g3dIVEmi(ivstep) = sqrt(Energy(i)-deltaE(ivstep)-Eiv(ivstep)-Ec)

            ! Bose-Einstein Distribution
            Niv=1/(exp(Eiv(ivstep)/(kb*T))-1.0)

            ! Intervalley Scattering Calculation
            GammaIVAbsCalc(i, ivstep) = pi*Div(ivstep)*hbar*Div(ivstep)/(2*rho*Eiv(ivstep))*ValleyN(ivstep) *  Niv
            GammaIVEmiCalc(i, ivstep) = pi*Div(ivstep)*hbar*Div(ivstep)/(2*rho*Eiv(ivstep))*ValleyN(ivstep) * (Niv+1)

        enddo

        ! Reformat (initial valley, final valley, Energy)
        ! From: GL GX LG LL LX XG XL XX
            GammaIVAbs(1, 2, i) = g3dIVAbs(1)*GammaIVAbsCalc(i, 1)*g3dIVFac(1)
            GammaIVAbs(1, 3, i) = g3dIVAbs(2)*GammaIVAbsCalc(i, 2)*g3dIVFac(2)
            GammaIVAbs(2, 1, i) = g3dIVAbs(3)*GammaIVAbsCalc(i, 3)*g3dIVFac(3)
            GammaIVAbs(2, 2, i) = g3dIVAbs(4)*GammaIVAbsCalc(i, 4)*g3dIVFac(4)
            GammaIVAbs(2, 3, i) = g3dIVAbs(5)*GammaIVAbsCalc(i, 5)*g3dIVFac(5)
            GammaIVAbs(3, 1, i) = g3dIVAbs(6)*GammaIVAbsCalc(i, 6)*g3dIVFac(6)
            GammaIVAbs(3, 2, i) = g3dIVAbs(7)*GammaIVAbsCalc(i, 7)*g3dIVFac(7)
            GammaIVAbs(3, 3, i) = g3dIVAbs(8)*GammaIVAbsCalc(i, 8)*g3dIVFac(8)

            GammaIVEmi(1, 2, i) = g3dIVEmi(1)*GammaIVEmiCalc(i, 1)*g3dIVFac(1)
            GammaIVEmi(1, 3, i) = g3dIVEmi(2)*GammaIVEmiCalc(i, 2)*g3dIVFac(2)
            GammaIVEmi(2, 1, i) = g3dIVEmi(3)*GammaIVEmiCalc(i, 3)*g3dIVFac(3)
            GammaIVEmi(2, 2, i) = g3dIVEmi(4)*GammaIVEmiCalc(i, 4)*g3dIVFac(4)
            GammaIVEmi(2, 3, i) = g3dIVEmi(5)*GammaIVEmiCalc(i, 5)*g3dIVFac(5)
            GammaIVEmi(3, 1, i) = g3dIVEmi(6)*GammaIVEmiCalc(i, 6)*g3dIVFac(6)
            GammaIVEmi(3, 2, i) = g3dIVEmi(7)*GammaIVEmiCalc(i, 7)*g3dIVFac(7)
            GammaIVEmi(3, 3, i) = g3dIVEmi(8)*GammaIVEmiCalc(i, 8)*g3dIVFac(8)

        ! Write GammaAcoustic
        if (valley.eq.1) then
            write(21, *) AcFac*GammaAcousticAbs(valley, i)
            write(24, *) AcFac*GammaAcousticEmi(valley, i)
        else if (valley.eq.2) then
            write(22, *) AcFac*GammaAcousticAbs(valley, i)
            write(25, *) AcFac*GammaAcousticEmi(valley, i)
        else if (valley.eq.3) then
            write(23, *) AcFac*GammaAcousticAbs(valley, i)
            write(26, *) AcFac*GammaAcousticEmi(valley, i)
        else
            print *, "Error: Unknown Valley Encountered (Acoustic)"
        end if

        ! Write GammaIonImp
        if (valley.eq.1) then
            write(31, *) IonFac*GammaIonImp(valley, i)
        else if (valley.eq.2) then
            write(32, *) IonFac*GammaIonImp(valley, i)
        else if (valley.eq.3) then
            write(33, *) IonFac*GammaIonImp(valley, i)
        else
            print *, "Error: Unknown Valley Encountered (Ion Imp)"
        end if

        ! Write GammaPop
        if (valley.eq.1) then
            write(41, *) PopFac*GammaPopAbs(valley, i)
            write(44, *) PopFac*GammaPopEmi(valley, i)
        else if (valley.eq.2) then
            write(42, *) PopFac*GammaPopAbs(valley, i)
            write(45, *) PopFac*GammaPopEmi(valley, i)
        else if (valley.eq.3) then
            write(43, *) PopFac*GammaPopAbs(valley, i)
            write(46, *) PopFac*GammaPopEmi(valley, i)
        else
            print *, "Error: Unknown Valley Encountered (Pop)"
        end if


        ! Write GammaIV
        if (valley.eq.1) then
        write(51, *) GammaIVAbs(1, 2, i)
        write(52, *) GammaIVAbs(1, 3, i)
        write(61, *) GammaIVEmi(1, 2, i)
        write(62, *) GammaIVEmi(1, 3, i)
        else if (valley.eq.2) then
        write(53, *) GammaIVAbs(2, 1, i)
        write(54, *) GammaIVAbs(2, 2, i)
        write(55, *) GammaIVAbs(2, 3, i)
        write(63, *) GammaIVEmi(2, 1, i)
        write(64, *) GammaIVEmi(2, 2, i)
        write(65, *) GammaIVEmi(2, 3, i)
        else if (valley.eq.3) then
        write(56, *) GammaIVAbs(3, 1, i)
        write(57, *) GammaIVAbs(3, 2, i)
        write(58, *) GammaIVAbs(3, 3, i)
        write(66, *) GammaIVEmi(3, 1, i)
        write(67, *) GammaIVEmi(3, 2, i)
        write(68, *) GammaIVEmi(3, 3, i)
        else
            print *, "Error: Unknown Valley Encountered (IV)"
        end if

    enddo
    enddo

    close(21)
    close(22)
    close(23)
    close(24)
    close(25)
    close(26)

    close(31)
    close(32)
    close(33)

    close(41)
    close(42)
    close(43)
    close(44)
    close(45)
    close(46)

    close(51)
    close(52)
    close(53)
    close(54)
    close(55)
    close(56)
    close(57)
    close(58)

    close(61)
    close(62)
    close(63)
    close(64)
    close(65)
    close(66)
    close(67)
    close(68)

end subroutine scatteringRates