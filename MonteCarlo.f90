program MonteCarlo
    use GaAsConstants
    implicit none

    real :: timeLeft

    call random_seed

    ! Calculate Scattering Rates
    call scatteringRates()
    ! Make Scattering Table
    call makeScatTable()

    ! Allocate Arrays for eN0 Particles
    allocate(eValley(eN0))
    allocate(eEnergy(eN0))
    allocate(eMomentum(eN0,3))
    allocate(eMomentumMag(eN0))
    allocate(ePhi(eN0))
    allocate(eTheta(eN0))
    allocate(eTff(eN0))


    ! Problem 1 Open
    print *, 'Problem 1 Opened'
    open(unit=11, file='Data/Problem1/Energy', status="unknown")
    open(unit=12, file='Data/Problem1/Momentum', status="unknown")

    do nEfield = 1, 1
        print *, 'Electric field =', Efield(nEfield)
        ! Initialize Particles (t = 1)
        do particle = 1, eN0
            print *, 'particle', particle
            ! Initialize Random Numbers
            call random_number(rEnergy)
            call random_number(rtheta)
            call random_number(rphi)
            call random_number(rteff)
                print *, 'Random Numbers =', rEnergy, rtheta, rphi, rteff
            ! Initialize Particles in Gamma Valley
                eValley(particle) = 1
                Valleyindex = eValley(particle)
            ! Initialize Randomized Particle Values
                eEnergy(particle) = -3.0/2.0*kb*T*log(rEnergy) ! eV
                ePhi(particle) = 2.0*pi*rphi
                eTheta(particle) = acos(1.0-2.0*rtheta)
            ! Calculate Momentum
                print *, sqrt(2.0*effm(Valleyindex)*eEnergy(particle))
                eMomentumMag(particle) = sqrt(2.0*effm(Valleyindex)*eEnergy(particle)) ! CHECK UNITS
                eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
                eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
                eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))
            ! Initialize Free Flight Time
                eTff(particle) = -1.0/Gamma0(Valleyindex)*log(rteff)
            ! Write Problem 1
                write(11, *) eEnergy(particle)
                write(12, *) eMomentum(particle,3)

        enddo

        close(11)
        close(12)


        ! Begin Time Stepping Loop (t > 1)
        do timestep = 1, maxtimestep
            do particle = 1, eN0
                if (eTff(particle)>=dt) then
                    ! Drift pz
                    eMomentum(particle,3) = eMomentum(particle,3) + (-q)*Efield(nEfield)*dt
                    eMomentumMag(particle) = sqrt(eMomentum(particle,1)**2.0+eMomentum(particle,2)**2.0+eMomentum(particle,3)**2.0)
                    eEnergy(particle) = (eMomentumMag(particle)**2.0)/(2.0*effm(Valleyindex))
                    eTff(particle) = eTff(particle) - dt
                else if (eTff(particle)<dt) then
                99  Valleyindex = eValley(particle)
                    timeLeft = eTff(particle)-dt
                    call random_number(rScat)
                    call random_number(rteff)
                    ! Drift to get pz and calculate E
                    eMomentum(particle,3) = eMomentum(particle,3) + (-q)*Efield(nEfield)*eTff(particle)
                    eMomentumMag(particle) = sqrt(eMomentum(particle,1)**2.0+eMomentum(particle,2)**2.0+eMomentum(particle,3)**2.0)
                    eEnergy(particle) = (eMomentumMag(particle)**2.0)/(2.0*effm(Valleyindex))
                    ! Choose Scattering Mechanism, update E, p, theta, phi
                    call chooseScatMech()
                    ! Update tff
                    eTff(particle) = -1.0/Gamma0(Valleyindex)*log(rteff)
                    ! Check if time less than the rest, if not, loop
                    if (eTff(particle)<=timeLeft) then
                        goto 99
                    else
                        ! Drift rest of way
                        eMomentum(particle,3) = eMomentum(particle,3) + (-q)*Efield(nEfield)*dt
                    endif
                else
                    print *, "Error: Particle Time Error in Time Stepping Loop"
                endif
            enddo
        enddo

        ! Write Values

    enddo
    
end program MonteCarlo