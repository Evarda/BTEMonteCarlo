program MonteCarlo
    use GaAsConstants
    implicit none

    ! Timestepping Loop
    real :: timeLeft

    ! Time Averages
    real :: timeavgvx = 0
    real :: timeavgvy = 0
    real :: timeavgvz = 0
    real :: timeavgKEG = 0
    real :: timeavgKEL = 0
    real :: timeavgKEX = 0
    real :: timeavgKEavg = 0
    real :: timeValleyPopG = 0
    real :: timeValleyPopL = 0
    real :: timeValleyPopX = 0


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

    do nEfield = 1, numE
        print *, 'Electric field =', Efield(nEfield)

        ! Clear Arrays
        timeavgvx = 0
        timeavgvy = 0
        timeavgvz = 0
        timeavgKEG = 0
        timeavgKEL = 0
        timeavgKEX = 0
        timeavgKEavg = 0
        timeValleyPopG = 0
        timeValleyPopL = 0
        timeValleyPopX = 0

        
        ! Initialize Particles (t = 1)
        do particle = 1, eN0
            print *, 'particle', particle
            ! Initialize Random Numbers
            call random_number(rEnergy)
            call random_number(rtheta)
            call random_number(rphi)
            call random_number(rteff)
            ! Initialize Particles in Gamma Valley
                eValley(particle) = 1
                Valleyindex = eValley(particle)
            ! Initialize Randomized Particle Values
                eEnergy(particle) = -3.0/2.0*kb*T*log(rEnergy) ! eV
                ePhi(particle) = 2.0*pi*rphi
                eTheta(particle) = acos(1.0-2.0*rtheta)
            ! Calculate Momentum
                eMomentumMag(particle) = sqrt(2.0*effm(Valleyindex)*eEnergy(particle))*sqrt(q)
                eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
                eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
                eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))
            ! Initialize Free Flight Time
                eTff(particle) = -1.0/Gamma0(Valleyindex)*log(rteff)
            if (nEfield.eq.1) then
            ! Write Problem 1
                write(11, *) eEnergy(particle)
                write(12, *) eMomentum(particle,3)
            endif

            ! Time Sums
            timeavgvx = timeavgvx + eMomentum(particle,1)/effm(Valleyindex)
            timeavgvy = timeavgvy + eMomentum(particle,2)/effm(Valleyindex)
            timeavgvz = timeavgvz + eMomentum(particle,3)/effm(Valleyindex)
            
            if (Valleyindex.eq.1) then
                timeavgKEG = timeavgKEG+eEnergy(particle)
                timeValleyPopG = timeValleyPopG + 1
            elseif (Valleyindex.eq.2) then
                timeavgKEL = timeavgKEG+eEnergy(particle)
                timeValleyPopL = timeValleyPopL + 1
            elseif (Valleyindex.eq.3) then
                timeavgKEX = timeavgKEG+eEnergy(particle)
                timeValleyPopX = timeValleyPopX + 1
            endif

            timeavgKEavg = timeavgKEavg + eEnergy(particle)
        enddo

        ! Time Averages
        timeavgvx = timeavgvx / eN0
        timeavgvy = timeavgvy / eN0
        timeavgvz = timeavgvz / eN0
        timeavgKEG = timeavgKEG / timeValleyPopG
        timeavgKEL = timeavgKEL / timeValleyPopL
        timeavgKEX = timeavgKEX / timeValleyPopX
        timeavgKEavg = timeavgKEavg / eN0

        write(13, *) Efield(nEfield)
        write(14, *) timeavgvx
        write(15, *) timeavgvy
        write(16, *) timeavgvz
        write(17, *) timeavgKEG
        write(18, *) timeavgKEL
        write(19, *) timeavgKEX
        write(20, *) timeavgKEavg
        write(21, *) timeValleyPopG
        write(22, *) timeValleyPopL
        write(23, *) timeValleyPopX


    close(11)
    close(12)


        ! Begin Time Stepping Loop (t > 1)
        do timestep = 1, maxtimestep

            ! Clear Arrays
            timeavgvx = 0
            timeavgvy = 0
            timeavgvz = 0
            timeavgKEG = 0
            timeavgKEL = 0
            timeavgKEX = 0
            timeavgKEavg = 0
            timeValleyPopG = 0
            timeValleyPopL = 0
            timeValleyPopX = 0

            do particle = 1, eN0
                if (eTff(particle)>=dt) then
                    ! Drift pz
                    eMomentum(particle,3) = eMomentum(particle,3) + (-q)*Efield(nEfield)*dt
                    eMomentumMag(particle) = sqrt((eMomentum(particle,1)/q)**2.0 + &
                                                  (eMomentum(particle,2)/q)**2.0 + &
                                                  (eMomentum(particle,3)/q)**2.0)*q
                    eEnergy(particle) = ((eMomentumMag(particle)/q)**2.0)/(2.0*effm(Valleyindex))*q
                    eTff(particle) = eTff(particle) - dt
                else if (eTff(particle)<dt) then
                99  Valleyindex = eValley(particle)
                    timeLeft = eTff(particle)-dt
                    call random_number(rScat)
                    call random_number(rteff)
                    ! Drift to get pz and calculate E
                    eMomentum(particle,3) = eMomentum(particle,3) + (-q)*Efield(nEfield)*eTff(particle)
                    eMomentumMag(particle) = sqrt((eMomentum(particle,1)/q)**2.0 + &
                                                  (eMomentum(particle,2)/q)**2.0 + &
                                                  (eMomentum(particle,3)/q)**2.0)*q
                    eEnergy(particle) = ((eMomentumMag(particle)/q)**2.0)/(2.0*effm(Valleyindex))*q
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

                print *, eMomentum(particle,1), eMomentum(particle,2), eMomentum(particle,3), &
                eMomentumMag(particle), eEnergy(particle)


                ! Time Sums
                timeavgvx = timeavgvx + eMomentum(particle,1)/effm(Valleyindex)
                timeavgvy = timeavgvy + eMomentum(particle,2)/effm(Valleyindex)
                timeavgvz = timeavgvz + eMomentum(particle,3)/effm(Valleyindex)
                
                if (Valleyindex.eq.1) then
                    timeavgKEG = timeavgKEG+eEnergy(particle)
                    timeValleyPopG = timeValleyPopG + 1
                elseif (Valleyindex.eq.2) then
                    timeavgKEL = timeavgKEG+eEnergy(particle)
                    timeValleyPopL = timeValleyPopL + 1
                elseif (Valleyindex.eq.3) then
                    timeavgKEX = timeavgKEG+eEnergy(particle)
                    timeValleyPopX = timeValleyPopX + 1
                endif

                timeavgKEavg = timeavgKEavg + eEnergy(particle)

            enddo

            ! Time Averages
            timeavgvx = timeavgvx / eN0
            timeavgvy = timeavgvy / eN0
            timeavgvz = timeavgvz / eN0
            timeavgKEG = timeavgKEG / timeValleyPopG
            timeavgKEL = timeavgKEL / timeValleyPopL
            timeavgKEX = timeavgKEX / timeValleyPopX
            timeavgKEavg = timeavgKEavg / eN0
    
            write(13, *) Efield(nEfield)
            write(14, *) timeavgvx
            write(15, *) timeavgvy
            write(16, *) timeavgvz
            write(17, *) timeavgKEG
            write(18, *) timeavgKEL
            write(19, *) timeavgKEX
            write(20, *) timeavgKEavg
            write(21, *) timeValleyPopG
            write(22, *) timeValleyPopL
            write(23, *) timeValleyPopX

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
    
end program MonteCarlo