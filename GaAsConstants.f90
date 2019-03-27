! Generic Values used in the Simulation
module GaAsConstants
    implicit none

    ! Constants
    real :: pi = 3.14159265359 ! Pi

    real, dimension(3) :: effm       ![kg] for gamma, L, X respectively

    ! Constants
    real :: hbar  = 6.582119e-16     ! [eV*s]
    real :: hbarJ = 1.0545718e-34    ! [J*s]
    real :: kb    = 8.6173303e-5     ! [eV/K]
    real :: kbJ   = 1.38064852e-23   ! [J/K]
    real :: T     = 300              ! [K]
    real :: ep0   = 8.854187817e-12  ! [F/m]
    real :: q     = 1.6021766208e-19 ! [C]
    real :: m0    = 9.10938356e-31   ! [kg]
    real :: Ec    = 0                ! Energy of the Conduction Band [eV]

    ! Number of Energy Steps
    integer :: nE=1000
    real, dimension(:), allocatable :: Energy
    real, dimension(:, :), allocatable :: k

    ! Scattering Table (Valley, Scattering Mechanism, Energy)
    real, dimension(:,:,:), allocatable :: ScatteringTable
    real, dimension(3) :: Gamma0

    ! Random Numbers
    real :: rEnergy
    real :: rtheta
    real :: rphi
    real :: rteff
    real :: rScat

    ! Particle Values
    integer :: eN0 = 1000
    integer, dimension(:), allocatable :: eValley
    integer :: Valleyindex
    real, dimension(:), allocatable :: eEnergy
    real, dimension(:,:), allocatable :: eMomentum
    real, dimension(:), allocatable :: eMomentumMag
    real, dimension(:), allocatable :: ePhi
    real, dimension(:), allocatable :: eTheta
    real, dimension(:), allocatable :: eTff

    integer :: numE = 6
    real, dimension(6) :: Efield = (/ 0.5e5, 1.0e5, 2.0e5, 5.0e5, 8.0e5, 10.0e5 /)

    ! Counter
    integer :: particle
    integer :: timestep
    integer :: nEfield
    integer :: maxtimestep = 10

    real :: dt = 1e-15



    contains

    subroutine allocateEnergy()
        allocate(Energy(nE))
        allocate(k(3, nE))
    end subroutine allocateEnergy
    
end module GaAsConstants