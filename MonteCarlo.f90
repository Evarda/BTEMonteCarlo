program MonteCarlo
    use GaAsConstants
    implicit none

    ! Timestepping Loop
    real :: timeLeft
    real :: timeLastScat

    ! Time Averages
    real :: vxSum = 0
    real :: vySum = 0
    real :: vzSum = 0
    real :: KEGSum = 0
    real :: KELSum = 0
    real :: KEXSum = 0
    real :: KESum = 0
    real :: ValleyPopG = 0
    real :: ValleyPopL = 0
    real :: ValleyPopX = 0
    real :: EAVG = 0
    real :: PAVG = 0
    real :: PZAVG = 0


    ! Call Random Seed
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
    open(unit=11, file='Data/Problem1/Energy', status="unknown")
    open(unit=12, file='Data/Problem1/Momentum', status="unknown")

    ! Problem 2 Open
    open(unit=13, file='Data/Problem2/Efield', status="unknown")
    open(unit=14, file='Data/Problem2/vx', status="unknown")
    open(unit=15, file='Data/Problem2/vy', status="unknown")
    open(unit=16, file='Data/Problem2/vz', status="unknown")
    open(unit=17, file='Data/Problem2/KEG', status="unknown")
    open(unit=18, file='Data/Problem2/KEL', status="unknown")
    open(unit=19, file='Data/Problem2/KEX', status="unknown")
    open(unit=20, file='Data/Problem2/KEavg', status="unknown")
    open(unit=21, file='Data/Problem2/ValleyPopG', status="unknown")
    open(unit=22, file='Data/Problem2/ValleyPopL', status="unknown")
    open(unit=23, file='Data/Problem2/ValleyPopX', status="unknown")
    open(unit=24, file='Data/Problem2/timestep', status="unknown")
    open(unit=25, file='Data/Problem2/time', status="unknown")
    open(unit=26, file='Data/Problem2/Estep', status="unknown")
    open(unit=27, file='Data/Problem2/EAVG', status="unknown")
    open(unit=28, file='Data/Problem2/PAVG', status="unknown")
    open(unit=29, file='Data/Problem2/PZAVG', status="unknown")

    do nEfield = 1, numE
        print *, 'Electric field =', Efield(nEfield)

        ! Clear Arrays
        vxSum = 0
        vySum = 0
        vzSum = 0
        KEGSum = 0
        KELSum = 0
        KEXSum = 0
        KESum = 0
        ValleyPopG = 0
        ValleyPopL = 0
        ValleyPopX = 0

        EAVG = 0
        PAVG = 0
        PZAVG = 0

        
        ! Initialize Particles (t = 1)
        do particle = 1, eN0
            !print *, 'particle', particle
            ! Initialize Random Numbers
            call random_number(rEnergy)
            call random_number(rtheta)
            call random_number(rphi)
            call random_number(rteff)
            ! Initialize Particles in Gamma Valley
                eValley(particle) = 1
            ! Initialize Randomized Particle Values
                eEnergy(particle) = -3.0/2.0*kb*T*log(rEnergy) ! eV
                ePhi(particle) = 2.0*pi*rphi
                eTheta(particle) = acos(1.0-2.0*rtheta)
            ! Calculate Momentum
                eMomentumMag(particle) = sqrt(2.0*effm(eValley(particle))*eEnergy(particle))*sqrt(q)
                eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
                eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
                eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))
            ! Initialize Free Flight Time
                eTff(particle) = -1.0/Gamma0(eValley(particle))*log(rteff)
            if (nEfield.eq.1) then

            ! Write Problem 1
                write(11, *) eEnergy(particle)
                write(12, *) eMomentum(particle,3)
            endif

            ! Time Sums
            vxSum = vxSum + eMomentum(particle,1)/effm(eValley(particle))
            vySum = vySum + eMomentum(particle,2)/effm(eValley(particle))
            vzSum = vzSum + eMomentum(particle,3)/effm(eValley(particle))
            
            KESum = KESum + eEnergy(particle)

            if (Valleyindex.eq.1) then
                KEGSum = KEGSum+eEnergy(particle)
                ValleyPopG = ValleyPopG + 1
            elseif (Valleyindex.eq.2) then
                KELSum = KELSum+eEnergy(particle)
                ValleyPopL = ValleyPopL + 1
            elseif (Valleyindex.eq.3) then
                KEXSum = KEXSum+eEnergy(particle)
                ValleyPopX = ValleyPopX + 1
            endif

            EAVG = EAVG + eEnergy(particle)
            PAVG = PAVG + eMomentumMag(particle)
            PZAVG = PZAVG + eMomentum(particle,3)
            
        enddo

        ! Write Problem 2 (for 1st timestep)
        write(13, *) Efield(nEfield)
        write(14, *) vxSum / eN0
        write(15, *) vySum / eN0
        write(16, *) vzSum / eN0
        write(17, *) KEGSum / ValleyPopG
        write(18, *) KELSum / ValleyPopL
        write(19, *) KEXSum / ValleyPopX
        write(20, *) KESum / eN0
        write(21, *) ValleyPopG
        write(22, *) ValleyPopL
        write(23, *) ValleyPopX
        write(24, *) 1
        write(25, *) 0.0
        write(26, *) nEfield
        write(27, *) EAVG  / eN0
        write(28, *) PAVG  / eN0
        write(29, *) PZAVG / eN0


        close(11)
        close(12)


        ! Begin Time Stepping Loop (t > 1)
        do timestep = 1, maxtimestep

            ! Clear Arrays
            vxSum = 0
            vySum = 0
            vzSum = 0
            KEGSum = 0
            KELSum = 0
            KEXSum = 0
            KESum = 0
            ValleyPopG = 0
            ValleyPopL = 0
            ValleyPopX = 0

            EAVG = 0
            PAVG = 0
            PZAVG = 0


            do particle = 1, eN0
                ! If Free Flight Time is bigger than the step size dt
                if (eTff(particle)>=dt) then
                    ! Drift pz by dt
                    eMomentum(particle,3) = eMomentum(particle,3) + (-q)*Efield(nEfield)*dt
                    ! Update p
                    eMomentumMag(particle) = sqrt((eMomentum(particle,1)/q)**2.0 + &
                                                  (eMomentum(particle,2)/q)**2.0 + &
                                                  (eMomentum(particle,3)/q)**2.0)*q
                    ! Update Energy
                    eEnergy(particle) = ((eMomentumMag(particle)/q)**2.0)/(2.0*effm(eValley(particle)))*q
                    ! Update Free Flight Time
                    eTff(particle) = eTff(particle) - dt

                ! If Free Flight Time is smaller than the step size dt
                else if (eTff(particle)<dt) then
                    ! Record time left in time step after scattering event
                    timeLeft = dt-eTff(particle)
                    
                    ! Drift pz by eTff
                99  eMomentum(particle,3) = eMomentum(particle,3) + (-q)*Efield(nEfield)*eTff(particle)
                    ! Update p
                    eMomentumMag(particle) = sqrt((eMomentum(particle,1)/q)**2.0 + &
                                                  (eMomentum(particle,2)/q)**2.0 + &
                                                  (eMomentum(particle,3)/q)**2.0)*q
                    ! Update Energy
                    eEnergy(particle) = ((eMomentumMag(particle)/q)**2.0)/(2.0*effm(eValley(particle)))*q
                    
                    ! Scatter (new theta and phi, E and p calculated)
                    call random_number(rScat)
                    call chooseScatMech()
                    
                    ! New Free Flight Time
                    call random_number(rteff)
                    eTff(particle) = -1.0/Gamma0(eValley(particle))*log(rteff)
                    
                    ! Update time left in time step after scattering event
                    timeLastScat = timeLeft
                    timeLeft = timeLeft-eTff(particle)
                    
                    ! Check if time left is positive, if so, scattering event happens again in this timestep
                    if (0<timeLeft) then
                        ! Repeat drift and scattering with remaining time
                        goto 99
                    else
                        if (timeLastScat<0) then
                            print *, "Error: Scattering time is negative!"
                        endif
                        ! Drift pz by the remaining time after the last scattering in timestep
                        eMomentum(particle,3) = eMomentum(particle,3) + (-q)*Efield(nEfield)*timeLastScat
                        ! Update p
                        eMomentumMag(particle) = sqrt((eMomentum(particle,1)/q)**2.0 + &
                                                      (eMomentum(particle,2)/q)**2.0 + &
                                                      (eMomentum(particle,3)/q)**2.0)*q
                        ! Update Energy
                        eEnergy(particle) = ((eMomentumMag(particle)/q)**2.0)/(2.0*effm(eValley(particle)))*q
                    endif
                else
                    print *, "Error: Particle Time Error in Time Stepping Loop"
                endif

                ! Time Sums
                vxSum = vxSum + eMomentum(particle,1)/effm(eValley(particle))
                vySum = vySum + eMomentum(particle,2)/effm(eValley(particle))
                vzSum = vzSum + eMomentum(particle,3)/effm(eValley(particle))
                
                KESum = KESum + eEnergy(particle)

                if (eValley(particle).eq.1) then
                    KEGSum = KEGSum+eEnergy(particle)
                    ValleyPopG = ValleyPopG + 1
                elseif (eValley(particle).eq.2) then
                    KELSum = KELSum+eEnergy(particle)
                    ValleyPopL = ValleyPopL + 1
                elseif (eValley(particle).eq.3) then
                    KEXSum = KEXSum+eEnergy(particle)
                    ValleyPopX = ValleyPopX + 1
                endif

                EAVG = EAVG + eEnergy(particle)
                PAVG = PAVG + eMomentumMag(particle)
                PZAVG = PZAVG + eMomentum(particle,3)

            enddo

            ! Write Problem 2 (for remaining timesteps)
            write(13, *) Efield(nEfield)
            write(14, *) vxSum / eN0
            write(15, *) vySum / eN0
            write(16, *) vzSum / eN0
            write(17, *) KEGSum / ValleyPopG
            write(18, *) KELSum / ValleyPopL
            write(19, *) KEXSum / ValleyPopX
            write(20, *) KESum / eN0
            write(21, *) ValleyPopG
            write(22, *) ValleyPopL
            write(23, *) ValleyPopX
            write(24, *) timestep
            write(25, *) timestep*dt
            write(26, *) nEfield

            write(27, *) EAVG  / eN0
            write(28, *) PAVG  / eN0
            write(29, *) PZAVG / eN0

        enddo

    enddo

    close(13)
    close(14)
    close(15)
    close(16)
    close(17)
    close(18)
    close(19)
    close(20)
    close(21)
    close(22)
    close(23)
    close(24)
    close(25)
    close(26)

    close(27)
    close(28)
    close(29)
    
end program MonteCarlo