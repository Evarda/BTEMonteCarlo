! Populates ScatteringTable in GaAsConstants with Scattering Rates from scatteringVariables
subroutine makeScatTable
    use GaAsConstants
    use scatteringVariables
    implicit none

    ! Counter
    integer :: ii        ! Energy
    integer :: valley   ! Valley

    allocate(ScatteringTable(3, nE, 10))

    ! Order
    !GammaAcoustic Abs
    !GammaAcoustic Emi
    !GammaPop Abs
    !GammaPop Emi
    !GammaIV Abs to 2
    !GammaIV Emi to 2
    !GammaIV Abs to 3
    !GammaIV Emi to 3
    !GammaIV Emi to 1
    !GammaIV Emi to 1


   ! ScatteringTable(Initial Valley, Energy, Rate #)
    do ii = 1, nE
        do valley = 1, 3
        ScatteringTable(valley, ii, 1) = GammaAcousticAbs(valley, ii)
        ScatteringTable(valley, ii, 2) = ScatteringTable(valley, ii, 1) + GammaAcousticEmi(valley, ii)
        ScatteringTable(valley, ii, 3) = ScatteringTable(valley, ii, 2) + GammaPopAbs(valley, ii)
        ScatteringTable(valley, ii, 4) = ScatteringTable(valley, ii, 3) + GammaPopEmi(valley, ii)
        ScatteringTable(valley, ii, 5) = ScatteringTable(valley, ii, 4) + GammaIVAbs(valley, 2, ii)
        ScatteringTable(valley, ii, 6) = ScatteringTable(valley, ii, 5) + GammaIVEmi(valley, 2, ii)
        ScatteringTable(valley, ii, 7) = ScatteringTable(valley, ii, 6) + GammaIVAbs(valley, 3, ii)
        ScatteringTable(valley, ii, 8) = ScatteringTable(valley, ii, 7) + GammaIVEmi(valley, 3, ii)
        enddo
        ScatteringTable(1, ii,  9) = ScatteringTable(1, ii, 8)
        ScatteringTable(1, ii, 10) = ScatteringTable(1, ii, 9)
        do valley = 2, 3
        ScatteringTable(valley, ii,  9) = ScatteringTable(valley, ii, 8) + GammaIVEmi(valley, 1, ii)
        ScatteringTable(valley, ii, 10) = ScatteringTable(valley, ii, 9) + GammaIVEmi(valley, 1, ii)
        enddo
    enddo

    do valley = 1, 3
        Gamma0(valley)=maxval(ScatteringTable(valley, :, 10))
        print *, 'Gamma0 = ', Gamma0(valley)
    enddo

    do ii = 1, nE
        do valley = 1, 3
        ScatteringTable(valley, ii, 1) = ScatteringTable(valley, ii, 1)/Gamma0(valley)
        ScatteringTable(valley, ii, 2) = ScatteringTable(valley, ii, 2)/Gamma0(valley)
        ScatteringTable(valley, ii, 3) = ScatteringTable(valley, ii, 3)/Gamma0(valley)
        ScatteringTable(valley, ii, 4) = ScatteringTable(valley, ii, 4)/Gamma0(valley)
        ScatteringTable(valley, ii, 5) = ScatteringTable(valley, ii, 5)/Gamma0(valley)
        ScatteringTable(valley, ii, 6) = ScatteringTable(valley, ii, 6)/Gamma0(valley)
        ScatteringTable(valley, ii, 7) = ScatteringTable(valley, ii, 7)/Gamma0(valley)
        ScatteringTable(valley, ii, 8) = ScatteringTable(valley, ii, 8)/Gamma0(valley)
        enddo
        ScatteringTable(1, ii,  9) = 0
        ScatteringTable(1, ii, 10) = 0
        do valley = 2, 3
        ScatteringTable(valley, ii,  9) = ScatteringTable(valley, ii, 9)/Gamma0(valley)
        ScatteringTable(valley, ii, 10) = ScatteringTable(valley, ii, 10)/Gamma0(valley)
        enddo
    enddo

    ! Write Scattering Table
    open(unit=70, file='Data/ScatTable/ScatTable_G0', status="unknown")
    open(unit=71, file='Data/ScatTable/ScatTable_G1', status="unknown")
    open(unit=72, file='Data/ScatTable/ScatTable_G2', status="unknown")
    open(unit=73, file='Data/ScatTable/ScatTable_G3', status="unknown")
    open(unit=74, file='Data/ScatTable/ScatTable_G4', status="unknown")
    open(unit=75, file='Data/ScatTable/ScatTable_G5', status="unknown")
    open(unit=76, file='Data/ScatTable/ScatTable_G6', status="unknown")
    open(unit=77, file='Data/ScatTable/ScatTable_G7', status="unknown")
    open(unit=78, file='Data/ScatTable/ScatTable_G8', status="unknown")
    open(unit=79, file='Data/ScatTable/ScatTable_G9', status="unknown")
    
    open(unit=80, file='Data/ScatTable/ScatTable_L0', status="unknown")
    open(unit=81, file='Data/ScatTable/ScatTable_L1', status="unknown")
    open(unit=82, file='Data/ScatTable/ScatTable_L2', status="unknown")
    open(unit=83, file='Data/ScatTable/ScatTable_L3', status="unknown")
    open(unit=84, file='Data/ScatTable/ScatTable_L4', status="unknown")
    open(unit=85, file='Data/ScatTable/ScatTable_L5', status="unknown")
    open(unit=86, file='Data/ScatTable/ScatTable_L6', status="unknown")
    open(unit=87, file='Data/ScatTable/ScatTable_L7', status="unknown")
    open(unit=88, file='Data/ScatTable/ScatTable_L8', status="unknown")
    open(unit=89, file='Data/ScatTable/ScatTable_L9', status="unknown")

    open(unit=90, file='Data/ScatTable/ScatTable_X0', status="unknown")
    open(unit=91, file='Data/ScatTable/ScatTable_X1', status="unknown")
    open(unit=92, file='Data/ScatTable/ScatTable_X2', status="unknown")
    open(unit=93, file='Data/ScatTable/ScatTable_X3', status="unknown")
    open(unit=94, file='Data/ScatTable/ScatTable_X4', status="unknown")
    open(unit=95, file='Data/ScatTable/ScatTable_X5', status="unknown")
    open(unit=96, file='Data/ScatTable/ScatTable_X6', status="unknown")
    open(unit=97, file='Data/ScatTable/ScatTable_X7', status="unknown")
    open(unit=98, file='Data/ScatTable/ScatTable_X8', status="unknown")
    open(unit=99, file='Data/ScatTable/ScatTable_X9', status="unknown")

    do ii = 1, nE
        write (70, *) ScatteringTable(1, ii,  1)
        write (71, *) ScatteringTable(1, ii,  2)
        write (72, *) ScatteringTable(1, ii,  3)
        write (73, *) ScatteringTable(1, ii,  4)
        write (74, *) ScatteringTable(1, ii,  5)
        write (75, *) ScatteringTable(1, ii,  6)
        write (76, *) ScatteringTable(1, ii,  7)
        write (77, *) ScatteringTable(1, ii,  8)
        write (78, *) ScatteringTable(1, ii,  9)
        write (79, *) ScatteringTable(1, ii, 10)

        write (80, *) ScatteringTable(2, ii,  1)
        write (81, *) ScatteringTable(2, ii,  2)
        write (82, *) ScatteringTable(2, ii,  3)
        write (83, *) ScatteringTable(2, ii,  4)
        write (84, *) ScatteringTable(2, ii,  5)
        write (85, *) ScatteringTable(2, ii,  6)
        write (86, *) ScatteringTable(2, ii,  7)
        write (87, *) ScatteringTable(2, ii,  8)
        write (88, *) ScatteringTable(2, ii,  9)
        write (89, *) ScatteringTable(2, ii, 10)

        write (90, *) ScatteringTable(3, ii,  1)
        write (91, *) ScatteringTable(3, ii,  2)
        write (92, *) ScatteringTable(3, ii,  3)
        write (93, *) ScatteringTable(3, ii,  4)
        write (94, *) ScatteringTable(3, ii,  5)
        write (95, *) ScatteringTable(3, ii,  6)
        write (96, *) ScatteringTable(3, ii,  7)
        write (97, *) ScatteringTable(3, ii,  8)
        write (98, *) ScatteringTable(3, ii,  9)
        write (99, *) ScatteringTable(3, ii, 10)
    enddo

    close(70)
    close(71)
    close(72)
    close(73)
    close(74)
    close(75)
    close(76)
    close(77)
    close(78)
    close(79)

    close(80)
    close(81)
    close(82)
    close(83)
    close(84)
    close(85)
    close(86)
    close(87)
    close(88)
    close(89)

    close(90)
    close(91)
    close(92)
    close(93)
    close(94)
    close(95)
    close(96)
    close(97)
    close(98)
    close(99)

end subroutine makeScatTable