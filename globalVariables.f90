module globalVariables
    implicit none

    ! Constants
    real :: pi = 3.14159265359 ! Pi

    ! Scattering Rates
    real, dimension(:,:), allocatable :: GammaAcoustic
    real, dimension(:,:), allocatable :: GammaMIonImp
    real, dimension(:,:), allocatable :: GammaMPop
    real, dimension(:,:), allocatable :: GammaTot
    
end module globalVariables