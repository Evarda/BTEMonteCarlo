module globalVariables
    implicit none

    ! Constants
    real :: pi = 3.14159265359 ! Pi

    ! Scattering Rates
    real, dimension(:,:), allocatable :: GammaAcousticAbs
    real, dimension(:,:), allocatable :: GammaAcousticEmi
    real, dimension(:,:), allocatable :: GammaMIonImp
    real, dimension(:,:), allocatable :: GammaMPop
    real, dimension(:,:), allocatable :: GammaTot
    
end module globalVariables