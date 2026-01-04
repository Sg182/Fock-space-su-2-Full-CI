program Fock_FCI
  use Precision
  use Constants
  use InputParams, only: PrintInput, NSites, Dim, Model, Delta, J2, Periodic, &
                         PBCx, PBCy, Nx, Ny, LanczosMaxIt, LanczosTol

  use Integrals
  use BuildHamMatVec
  use CIUtils
  use LanczosSolver

  implicit none

  integer :: NS, iters
  integer :: NDetSect
  
  real(kind=pr) :: E_Sz1, resid
  integer, allocatable :: basis(:)
  real(kind=pr), allocatable :: C0(:)
  integer :: p
  real(kind=pr) :: Szp
  real(kind=pr) :: sumSz

  call PrintInput()

  NS = NSites()    ! <<< MUST set NS first

  call SetUpIntegrals(NS)

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

  !======== Set Sz sector (Sz_tot = +1) ==================
  call SetSzSector(NS, 1)
  NDetSect = GetSzSectorDim()
  write(*,'(A,I0)') "NDet(Sz=1) = ", NDetSect
  !=======================================================

  allocate(C0(NDetSect))   ! <<< sector size

  call LanczosGroundState(NS, NDetSect, ApplyHamiltonian_SzSector, &
                          LanczosMaxIt, LanczosTol, E_Sz1, iters, resid, C0)

  write(*,'(A,F22.16)') "Ground-state energy (Sz=1) = ", E_Sz1
  write(*,'(A,I0)')     "Lanczos iterations  = ", iters

  call GetSzSectorBasis(basis)
  call WriteTopCICoeffs_WithBasis("CI_top.dat", NS, C0, basis,30)

  ! Sz expectation values: MUST be sector-aware
  do p = 0, NS-1
    Szp = Expect_Sz_Sector(NS, C0, p)
    write(*,'(A,I0,A,F12.8)') "p=", p, "  <Sz_p>=", Szp
  end do

  
  sumSz = 0.0_pr
  do p = 0, NS-1
    sumSz = sumSz + Expect_Sz_Sector(NS, C0, p)
  end do
  write(*,'(A,F12.8)') "<Sz_tot> = ", sumSz


  deallocate(C0)
  call ClearSzSector()

  call ShutDownIntegrals()
end program Fock_FCI
