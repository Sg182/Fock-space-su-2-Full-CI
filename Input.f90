module InputParams
  use Precision
  implicit none

  !-------------------------
  ! This Module controls the Hamiltonian and Geometry
  !-------------------------
  
  !========================================================================
  ! Model meanings:                                                       !
  !   Dim=1: 1=XXZ, 2=J1J2XXZ, 3=XXZ_P+JW_OBC (P+JW Hamiltonian)      !
  !   Dim=2: 1=XXZ, 2=J1J2XXZ                                             !
  !=======================================================================!
  integer :: Model = 3

  integer :: Dim = 1             ! 1 => 1D chain, 2 => 2D square
  integer :: NLevels = 6          ! used when Dim=1
  integer :: Nx = 4, Ny = 4       ! used when Dim=2   (NSites = Nx*Ny)

  ! Boundary conditions
  logical :: Periodic  = .False.  ! for Dim=1
  logical :: PBCx      = .True.  ! Modify for Dim=2
  logical :: PBCy      = .True.  ! Modify for Dim=2

  ! Hamiltonian parameters
  real(kind=pr) :: Delta = 0.80_pr
  real(kind=pr) :: J2    = 0.0_pr

  ! =====Optional model selector (if you want a clean switch in Main)======
  !  Dim=1 (1: 1D XXZ, 2: 1D J1J2XXZ), Dim=2(1: 2D XXZ, 2:2D J1J2)
   

  integer       :: LanczosMaxIt = 200
   real(kind=pr) :: LanczosTol   = 1.0e-12_pr

contains

  !=========================================================
  ! Return total number of sites/pair-levels used by Fock code
  !=========================================================
  integer function NSites()
    implicit none
    if (Dim == 1) then
      NSites = NLevels
    else if (Dim == 2) then
      NSites = Nx*Ny
    else
      stop "InputParams: Dim must be 1 or 2"
    end if
  end function NSites

  subroutine PrintInput()
    implicit none
    write(*,'(A,I0)')     "Dim       = ", Dim
    if (Dim == 1) then
      write(*,'(A,I0)')   "NLevels   = ", NLevels
      write(*,'(A,L1)')   "Periodic  = ", Periodic
    else
      write(*,'(A,I0)')   "Nx        = ", Nx
      write(*,'(A,I0)')   "Ny        = ", Ny
      write(*,'(A,I0)')   "NSites    = ", Nx*Ny
      write(*,'(A,L1)')   "PBCx      = ", PBCx
      write(*,'(A,L1)')   "PBCy      = ", PBCy
    end if
    write(*,'(A,I0)')     "Model     = ", Model
    write(*,'(A,F8.4)')   "Delta     = ", Delta
    write(*,'(A,F8.4)')   "J2        = ", J2
  end subroutine PrintInput

end module InputParams
