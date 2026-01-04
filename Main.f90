!=========================================================================================================
  !USE THIS FILE TO RUN THE CALCULATION IN THE TOTAL FOCK SPACE (OR INCLUDE ALL Sz SECTORS IN THE SPACE)

program Fock_FCI
  use Precision
  use Constants
  use InputParams, only: PrintInput, NSites, Dim, Model, Delta, J2, Periodic, &
                       PBCx, PBCy, Nx, Ny, LanczosMaxIt, LanczosTol

  use Integrals
  use BuildHamMatVec  
  use CIUtils 
       ! provides ApplyHamiltonian(NS, v, w)
  use LanczosSolver           ! provides LanczosGroundState(...)

  implicit none

  integer :: NDet
  integer :: NS, iters
  real(kind=pr) :: E0, resid
  real(kind=pr), allocatable :: C0(:)
   
  integer :: p
  real(kind=pr) :: Szp

  call PrintInput()

  NS = NSites()
  NDet = ishft(1, NS)
  write(*,'(A,I0)') "NDet      = ", NDet

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

  !---------------------------------------------
  ! Lanczos ground state + CI coefficients
  !---------------------------------------------
  allocate(C0(NDet))

  call LanczosGroundState(NS, NDet, ApplyHamiltonian, LanczosMaxIt, LanczosTol, &
                          E0, iters, resid, C0)

  write(*,'(A,F22.16)') "Ground-state energy = ", E0
  write(*,'(A,I0)')     "Lanczos iterations  = ", iters
  !write(*,'(A,ES12.4)') "Ritz diff (proxy)   = ", resid

  ! Optional: dump coefficients to a file
  call WriteTopCICoeffs("CI_top.dat", NS, C0, 30)
   
!==============================================================
               !PRINTS OCCUPATION NUMBER Sz

!do p = 0, NS-1
!  Szp = Expect_Sz(C0, p)
!  write(*,'(A,I0,A,F12.8)') "p=", p, "  <Sz_p>=", Szp
!end do
!===============================================================


  deallocate(C0)

  call ShutDownIntegrals()
end program Fock_FCI

 