! Holds Scattering Rates for Each Scattering Mechanism
module scatteringVariables
    implicit none

    ! Scattering Rates
    real, dimension(:,:), allocatable :: GammaAcousticAbs
    real, dimension(:,:), allocatable :: GammaAcousticEmi
    real, dimension(:,:), allocatable :: GammaIonImp
    real, dimension(:,:), allocatable :: GammaPopAbs
    real, dimension(:,:), allocatable :: GammaPopEmi

end module scatteringVariables