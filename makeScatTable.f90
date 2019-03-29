! Populates ScatteringTable in GaAsConstants with Scattering Rates from scatteringVariables
subroutine makeScatTable
    use GaAsConstants
    use scatteringVariables
    implicit none

    ! Counter
    integer :: ii        ! Energy
    integer :: valley   ! Valley

    allocate(ScatteringTable(3, nE, 10))

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
    open(unit=70, file='Data/ScatTable/ScatTable_0', status="unknown")
    open(unit=71, file='Data/ScatTable/ScatTable_1', status="unknown")
    open(unit=72, file='Data/ScatTable/ScatTable_2', status="unknown")
    open(unit=73, file='Data/ScatTable/ScatTable_3', status="unknown")
    open(unit=74, file='Data/ScatTable/ScatTable_4', status="unknown")
    open(unit=75, file='Data/ScatTable/ScatTable_5', status="unknown")
    open(unit=76, file='Data/ScatTable/ScatTable_6', status="unknown")
    open(unit=77, file='Data/ScatTable/ScatTable_7', status="unknown")
    open(unit=78, file='Data/ScatTable/ScatTable_8', status="unknown")
    open(unit=79, file='Data/ScatTable/ScatTable_9', status="unknown")
    do ii = 1, nE
        write (70, *) ScatteringTable(2, ii,  1)
        write (71, *) ScatteringTable(2, ii,  2)
        write (72, *) ScatteringTable(2, ii,  3)
        write (73, *) ScatteringTable(2, ii,  4)
        write (74, *) ScatteringTable(2, ii,  5)
        write (75, *) ScatteringTable(2, ii,  6)
        write (76, *) ScatteringTable(2, ii,  7)
        write (77, *) ScatteringTable(2, ii,  8)
        write (78, *) ScatteringTable(2, ii,  9)
        write (79, *) ScatteringTable(2, ii, 10)
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

end subroutine makeScatTable