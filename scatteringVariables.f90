! Holds Scattering Rates for Each Scattering Mechanism
module scatteringVariables
    implicit none

    ! Scattering Rates
    real, dimension(:,:), allocatable :: GammaAcousticAbs
    real, dimension(:,:), allocatable :: GammaAcousticEmi
    real, dimension(:,:), allocatable :: GammaIonImp
    real, dimension(:,:), allocatable :: GammaPopAbs
    real, dimension(:,:), allocatable :: GammaPopEmi
    real, dimension(:,:), allocatable :: GammaIVAbsCalc
    real, dimension(:,:), allocatable :: GammaIVEmiCalc
    real, dimension(:,:,:), allocatable :: GammaIVAbs
    real, dimension(:,:,:), allocatable :: GammaIVEmi

end module scatteringVariables