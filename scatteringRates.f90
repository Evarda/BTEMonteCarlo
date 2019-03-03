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

    ! Acoustic Phonon Scattering
    real, dimension(3) :: Dac = (/ 7.01, 9.2, 9.0 /) ! [eV] for gamma, L, X respectively
    
    ! Ionized Impurity Scattering
    real :: NI = 1e23 ! [1/m^3] (Equiv to 10^17 1/cm^3)
    real :: Z = 1

    ! Polar Optical Phonon Scattering
    real :: E0 = 3.536e-2 ![eV]
    real :: w0
    real :: N0

    ! Counter
    integer :: i
    integer :: valley

    ! GaAs Constants
    effm = (/ 0.067, 0.022, 0.058 /)*m0
    rho = 5.36e3 ![kg/m^3]
    vs = 5.24e3 ![m/s]

    ! Ionized Impurity Scattering
    Ld=sqrt(ep0*eprInf*kbJ*T/(e**2.0*NI)) ! [m]

    ! Polar Optical Phonon Scattering
    w0 = E0/hbar ! [1/s]
    N0 = (exp(E0/(kb*T))-1)**(-1)

    print*, w0, N0

    allocate(Energy(nE), k(3, nE))

    ! Initialize Energy and Wavevector
    do i = 1, nE
        Energy(i) = (2.0/nE)*i
        do valley = 1, 3
            k(valley, i) = sqrt(2.0*effm(valley)/hbarJ*Energy(i)/hbar) 
        enddo
    enddo

    allocate(GammaAcousticAbs(3,nE), GammaAcousticEmi(3,nE), &
             GammaIonImp(3,nE), &
             GammaPopAbs(3,nE), GammaPopEmi(3,nE), &
             )
    
    ! Write Energy
    open(unit=10, file='Data/Energy', status="unknown")
    
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

    do valley = 1, 3
    g3dFac = 1.0/(2.0*(pi**2.0))*(2.0*effm(valley)/(hbar**2.0))**(1.5)
    AcFac = pi/(2.0*(hbar)**(0.5)*(hbarJ)**(0.5))*(Dac(valley)**2.0)*kb*T/(rho*(vs**2.0))
    IonFac = (Z**2.0)/(pi*eprInf**2.0)*(e**2.0/hbarJ)*(NI*e)*(effm(valley)/hbarJ)*((Ld**4.0)/hbarJ)*(e/ep0)*(1/ep0)
    PopFac = sqrt(2.0)/(8.0*pi)*(1.0/eprInf-1.0/epr0)*(e/sqrt(hbarJ))*(e/ep0)*(w0*sqrt(hbar))*(sqrt(effm(valley))/hbarJ)
    
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


    
        ! Write Energy 
        write(10, *) Energy(i)

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

    enddo
    enddo

    close(10)

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

end subroutine scatteringRates