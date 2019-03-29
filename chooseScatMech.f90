subroutine chooseScatMech
    use GaAsConstants
    use scatteringVariables
    implicit none

    ! Counter
    integer :: index
    integer :: ii

    ! Values
    real :: f

    ! Scattering Energies
    real :: EPOP = 3.536e-2 ![eV]

    ! Intervalley Rates
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

    real, dimension(3, 3) :: Eiv
    real, dimension(3, 3) :: deltaE

    real :: pxnew
    real :: pynew
    real :: pznew
    real :: pxold
    real :: pyold
    real :: pzold
    real :: pold
    real :: pnew

    real :: Eold
    
    Eiv(1, 1) = 0
    Eiv(1, 2) = EivGL
    Eiv(1, 3) = EivGX
    Eiv(2, 1) = EivLG
    Eiv(2, 2) = EivLL
    Eiv(2, 3) = EivLX
    Eiv(3, 1) = EivXG
    Eiv(3, 2) = EivXL
    Eiv(3, 3) = EivXX

    deltaE(1,1) = 0
    deltaE(1,2) = deltaEGL
    deltaE(1,3) = deltaEGX
    deltaE(2,1) = deltaELG
    deltaE(2,2) = deltaELL
    deltaE(2,3) = deltaELX
    deltaE(3,1) = deltaEXG
    deltaE(3,2) = deltaEXL
    deltaE(3,3) = deltaEXX
    


    ! Find Energy Index
    index = 0
    do ii = 1, nE
        if (Energy(ii)<=eEnergy(particle)) then
            index = index+1
        elseif (Energy(ii)>eEnergy(particle)) then
            exit
        endif
    enddo
    if (index>nE) then
        index = nE
        print *, "Index > nE"
    endif
    if(index.eq.0) then
        index = 1
    endif

    ! Choose Gamma


    ! Update From Scattering Mechanism
        ! Update E
        ! Update theta
        ! Update phi
        ! Update p

    if (rScat<=ScatteringTable(eValley(particle), index, 1)) then
        mechanism = 1
        ! AcousticAbs
            ! Generate New Theta, Phi
            call random_number(rtheta)
            call random_number(rphi)
            ePhi(particle) = 2.0*pi*rphi
            eTheta(particle) = acos(1.0-2.0*rtheta)
            ! Calculate New Components of Momentum
            eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
            eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
            eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))
        
    elseif (rScat<=ScatteringTable(eValley(particle), index, 2)) then
        mechanism = 2
        ! AcoutsticEmi
            ! Generate New Theta, Phi
            call random_number(rtheta)
            call random_number(rphi)
            ePhi(particle) = 2.0*pi*rphi
            eTheta(particle) = acos(1.0-2.0*rtheta)
            ! Calculate New Components of Momentum
            eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
            eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
            eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))

    elseif (rScat<=ScatteringTable(eValley(particle), index, 3)) then
        mechanism = 3
        ! GammaPopAbs
            ! Calculate Energy, Momentum
            Eold = eEnergy(particle)
            eEnergy(particle) = eEnergy(particle)+EPOP
            ! Generate New Theta, Phi
            call random_number(rtheta)
            call random_number(rphi)
            ePhi(particle) = 2.0*pi*rphi
            f = 2.0*sqrt(Eold*eEnergy(particle))/((sqrt(Eold)-sqrt(eEnergy(particle)))**2)
            eTheta(particle) = acos((1.0+f-(1.0+2.0*f)**rtheta)/f)
            ! Calculate Momentum
            eMomentumMag(particle) = sqrt(2.0*effm(eValley(particle))/q*eEnergy(particle))
            pxold = eMomentum(particle,1)
            pyold = eMomentum(particle,2)
            pzold = eMomentum(particle,3)
            pxnew = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
            pynew = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
            pznew = eMomentumMag(particle)*cos(eTheta(particle))
            pold = sqrt(pxold**2+pyold**2+pzold**2)
            pnew = eMomentumMag(particle)
            eMomentum(particle,1) = pxnew*(pyold/sqrt(pxold**2.0+pyold**2.0)) + &
                                    pynew*(pxold*pzold/pold/sqrt(pxold**2.0+pyold**2.0)) + &
                                    pznew*(pxold/pold)
            eMomentum(particle,2) = pxnew*(-pxold/sqrt(pxold**2.0+pyold**2.0)) + &
                                    pynew*(pyold*pzold/pold/sqrt(pxold**2.0+pyold**2.0)) + &
                                    pznew*(pyold/pold)
            eMomentum(particle,3) = pynew*(-sqrt(pxold**2.0+pyold**2.0)/pold) + &
                                    pznew*(pzold/pold)
            
            !eMomentum(particle,1) = pxnew*(pzold*pxold/sqrt(pxold**2.0+pyold**2.0)/pold) + &
            !                        pynew*(-pyold/sqrt(pxold**2.0+pyold**2.0)) + &
            !                        pznew*(pxold/pold)
            !eMomentum(particle,2) = pxnew*(pyold*pzold/sqrt(pxold**2.0+pyold**2.0)/pold) + &
            !                        pynew*(pxold/sqrt(pxold**2.0+pyold**2.0)) + &
            !                        pznew*(pyold/pold)
            !eMomentum(particle,3) = pxnew*(-sqrt(pxold**2.0+pyold**2.0)/pold) + &
            !                        pznew*(pzold/pold)


    elseif (rScat<=ScatteringTable(eValley(particle), index, 4)) then
        mechanism = 4
        ! GammaPopEmi
            ! Calculate Energy
            Eold = eEnergy(particle)
            eEnergy(particle) = eEnergy(particle)-EPOP
            ! Generate New Theta, Phi
            call random_number(rtheta)
            call random_number(rphi)
            ePhi(particle) = 2.0*pi*rphi
            f = 2.0*sqrt(Eold*eEnergy(particle))/((sqrt(Eold)-sqrt(eEnergy(particle)))**2)
            !print *, "Theta", Eold, EPOP
            eTheta(particle) = acos((1.0+f-((1.0+2.0*f)**rtheta))/f)
            

            ! Calculate Momentum
            eMomentumMag(particle) = sqrt(2.0*effm(eValley(particle))/q*eEnergy(particle))
            pxold = eMomentum(particle,1)
            pyold = eMomentum(particle,2)
            pzold = eMomentum(particle,3)
            pxnew = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
            pynew = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
            pznew = eMomentumMag(particle)*cos(eTheta(particle))
            pold = sqrt(pxold**2+pyold**2+pzold**2)
            pnew = eMomentumMag(particle)
            eMomentum(particle,1) = pxnew*(pyold/sqrt(pxold**2.0+pyold**2.0)) + &
                                    pynew*(pxold*pzold/pold/sqrt(pxold**2.0+pyold**2.0)) + &
                                    pznew*(pxold/pold)
            eMomentum(particle,2) = pxnew*(-pxold/sqrt(pxold**2.0+pyold**2.0)) + &
                                    pynew*(pyold*pzold/pold/sqrt(pxold**2.0+pyold**2.0)) + &
                                    pznew*(pyold/pold)
            eMomentum(particle,3) = pynew*(-sqrt(pxold**2.0+pyold**2.0)/pold) + &
                                    pznew*(pzold/pold)


            !eMomentum(particle,1) = pxnew*(pzold*pxold/sqrt(pxold**2.0+pyold**2.0)/pold) + &
            !                        pynew*(-pyold/sqrt(pxold**2.0+pyold**2.0)) + &
            !                        pznew*(pxold/pold)
            !eMomentum(particle,2) = pxnew*(pyold*pzold/sqrt(pxold**2.0+pyold**2.0)/pold) + &
            !                        pynew*(pxold/sqrt(pxold**2.0+pyold**2.0)) + &
            !                        pznew*(pyold/pold)
            !eMomentum(particle,3) = pxnew*(-sqrt(pxold**2.0+pyold**2.0)/pold) + &
            !                        pznew*(pzold/pold)

    elseif (rScat<=ScatteringTable(eValley(particle), index, 5)) then
        mechanism = 5
        ! GammaIVAbs(to 2) Energy(i)-deltaE(ivstep)+Eiv(ivstep)

        ! Calculate Energy
        eEnergy(particle) = eEnergy(particle)-deltaE(eValley(particle), 2)+Eiv(eValley(particle), 2)
        ! Change Valley
        eValley(particle) = 2

        ! Generate New Theta, Phi
        call random_number(rtheta)
        call random_number(rphi)
        ePhi(particle) = 2.0*pi*rphi
        eTheta(particle) = acos(1.0-2.0*rtheta)
        ! Calculate New Components of Momentum
        eMomentumMag(particle) = sqrt(2.0*effm(eValley(particle))/q*eEnergy(particle))
        eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))

    elseif (rScat<=ScatteringTable(eValley(particle), index, 6)) then
        mechanism = 6
        ! GammaIVEmi(to 2)

        ! Calculate Energy
        eEnergy(particle) = eEnergy(particle)-deltaE(eValley(particle), 2)-Eiv(eValley(particle), 2)
        !print *, "Intervalley Scattering", eEnergy(particle)
        ! Change Valley
        eValley(particle) = 2

        ! Generate New Theta, Phi
        call random_number(rtheta)
        call random_number(rphi)
        ePhi(particle) = 2.0*pi*rphi
        eTheta(particle) = acos(1.0-2.0*rtheta)
        ! Calculate New Components of Momentum
        eMomentumMag(particle) = sqrt(2.0*effm(eValley(particle))/q*eEnergy(particle))
        eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))

        

    elseif (rScat<=ScatteringTable(eValley(particle), index, 7)) then
        mechanism = 7
        ! GammaIVAbs(to 3)

        ! Calculate Energy
        eEnergy(particle) = eEnergy(particle)-deltaE(eValley(particle), 3)+Eiv(eValley(particle), 3)
        !print *, "Intervalley Scattering", eEnergy(particle)
        ! Change Valley
        eValley(particle) = 3

        ! Generate New Theta, Phi
        call random_number(rtheta)
        call random_number(rphi)
        ePhi(particle) = 2.0*pi*rphi
        eTheta(particle) = acos(1.0-2.0*rtheta)
        ! Calculate New Components of Momentum
        eMomentumMag(particle) = sqrt(2.0*effm(eValley(particle))/q*eEnergy(particle))
        eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))

        
        
    elseif (rScat<=ScatteringTable(eValley(particle), index, 8)) then
        mechanism = 8
        ! GammaIVEmi(to 3)
        ! Calculate Energy

        eEnergy(particle) = eEnergy(particle)-deltaE(eValley(particle), 3)-Eiv(eValley(particle), 3)
        !print *, "Intervalley Scattering", eEnergy(particle)
        ! Change Valley
        eValley(particle) = 3

        ! Generate New Theta, Phi
        call random_number(rtheta)
        call random_number(rphi)
        ePhi(particle) = 2.0*pi*rphi
        eTheta(particle) = acos(1.0-2.0*rtheta)
        ! Calculate New Components of Momentum
        eMomentumMag(particle) = sqrt(2.0*effm(eValley(particle))/q*eEnergy(particle))
        eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))



    elseif (rScat<=ScatteringTable(eValley(particle), index, 9)) then
        mechanism = 9
        ! GammaIVAbs(to 1) NOT FOR 1
        if (Valleyindex.eq.1) then
            print *, "ERROR: Gamma-Gamma Scattering DOES NOT EXIST"
            call abort
        endif

        ! Calculate Energy
        eEnergy(particle) = eEnergy(particle)-deltaE(eValley(particle), 1)+Eiv(eValley(particle), 1)
        !print *, "Intervalley Scattering", eEnergy(particle)
        ! Change Valley
        eValley(particle) = 1

        ! Generate New Theta, Phi
        call random_number(rtheta)
        call random_number(rphi)
        ePhi(particle) = 2.0*pi*rphi
        eTheta(particle) = acos(1.0-2.0*rtheta)
        ! Calculate New Components of Momentum
        eMomentumMag(particle) = sqrt(2.0*effm(eValley(particle))/q*eEnergy(particle))
        eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))

        

    elseif (rScat<=ScatteringTable(eValley(particle), index, 10)) then
        mechanism = 10
        ! GammaIVEmi(to 1) NOT FOR 1
        if (eValley(particle).eq.1) then
            print *, "ERROR: Gamma-Gamma Scattering DOES NOT EXIST"
            call abort
        endif

        ! Calculate Energy
        eEnergy(particle) = eEnergy(particle)-deltaE(eValley(particle), 1)-Eiv(eValley(particle), 1)
        !print *, "Intervalley Scattering", eEnergy(particle)
        ! Change Valley
        eValley(particle) = 1


        ! Generate New Theta, Phi
        call random_number(rtheta)
        call random_number(rphi)
        ePhi(particle) = 2.0*pi*rphi
        eTheta(particle) = acos(1.0-2.0*rtheta)
        ! Calculate New Components of Momentum
        eMomentumMag(particle) = sqrt(2.0*effm(eValley(particle))/q*eEnergy(particle))
        eMomentum(particle,1) = eMomentumMag(particle)*cos(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,2) = eMomentumMag(particle)*sin(ePhi(particle))*sin(eTheta(particle))
        eMomentum(particle,3) = eMomentumMag(particle)*cos(eTheta(particle))

    else
        ! Self Scattering
        mechanism = 11
    endif



end subroutine chooseScatMech