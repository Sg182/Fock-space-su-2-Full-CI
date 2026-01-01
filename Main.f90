program Fock_FCI
  use Precision
  use Constants
  use InputParams, only: PrintInput, NSites, Dim, Model, Delta, J2, Periodic, &
                       PBCx, PBCy, Nx, Ny, LanczosMaxIt, LanczosTol

  !use Operators
  use Integrals

  ! NEW:
  use BuildHamMatVec          ! provides ApplyHamiltonian(NS, v, w)
  use LanczosSolver          ! provides LanczosGroundState(...)

  implicit none

  integer :: NDet
  integer :: NS
  real(kind=pr) :: E0

  call PrintInput()

  ! total number of pair-levels/sites used by the Fock space
  NS = NSites()

  NDet = ishft(1, NS)
  write(*,'(A,I0)') "NDet      = ", NDet

  !call SetUpPairOperators(NS)
  call SetUpIntegrals(NS)

  !============================================
  ! Build integrals depending on geometry/model
  !============================================
  select case (Dim)

  case (1)
    if (Model == 1) then
      call DoIntegralsXXZ_1D(NS, Delta, Periodic)
    else if (Model == 2) then
      call DoIntegralsJ1J2XXZ_1D(NS, Delta, J2, Periodic)
    else
      stop "Main: unknown Model (use 1=XXZ, 2=J1J2XXZ)"
    end if

  case (2)
    if (Model == 1) then
      call DoIntegralsXXZ_2D(Nx, Ny, Delta, PBCx, PBCy)
    else if (Model == 2) then
      call DoIntegralsJ1J2XXZ_2D(Nx, Ny, Delta, J2, PBCx, PBCy)
    else
      stop "Main: unknown Model (use 1=XXZ, 2=J1J2XXZ)"
    end if

  case default
    stop "Main: Dim must be 1 or 2"

  end select

  !---------------------------------------------
  ! Lanczos ground state (NO dense H)
  !---------------------------------------------
  call LanczosGroundState(NS, NDet, ApplyHamiltonian, LanczosMaxIt, LanczosTol, E0)

  write(*,'(A,F22.16)') "Ground-state energy = ", E0

  call ShutDownIntegrals()
  !call ShutDownPairOperators()

end program Fock_FCI
