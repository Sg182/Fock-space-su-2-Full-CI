!=========================================================================================================
  !USE THIS FILE TO RUN THE CALCULATION IN THE TOTAL FOCK SPACE (OR INCLUDE ALL Sz SECTORS IN THE SPACE)

program Fock_FCI
  use Precision
  use Constants
  use InputParams, only: PrintInput, NSites, Dim, Model, Delta, J2, Periodic, &
                       PBCx, PBCy, Nx, Ny, LanczosMaxIt, LanczosTol, isCorr

  use Integrals
  use BuildHamMatVec  !Provides ApplyHamiltonian
  use CIUtils 
  use LanczosSolver           ! provides LanczosGroundState(...)

  implicit none

  integer :: NDet
  integer :: NS, iters
  real(kind=pr) :: E0, resid
  real(kind=pr), allocatable :: C0(:)
  !real(kind=pr), allocatable :: SpSq(:,:)
  !integer :: q
  !real(kind=pr), allocatable :: Corr(:)
  integer :: iPBC = 0
  Real(Kind=pr), Allocatable :: SpSq(:,:)
  Real(Kind=pr), Allocatable :: Corr1D(:)
  Real(Kind=pr), Allocatable :: Corr2D(:,:)
  Integer :: p, dx, dy

  real(kind=pr) :: normsq
  logical :: is_norm
  real(kind=pr) :: tol
   
   
  real(kind=pr) :: Szp
  

  if (Periodic) then
    iPBC = 1
  endif

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
    else if (Model == 3) then
    ! New Hamiltonian: no Integrals build needed (unless you want H000 shift)
    ! Optionally: Periodic must be false for your OBC version
      if (Periodic) stop "Model=3: Periodic not implemented"
    else if (Model == 4) Then
      If (Periodic) stop "Model=4: Periodic not implemented"

    else
      stop "Main: unknown Model (use 1=XXZ, 2=J1J2XXZ or 3)"
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

  !--------------------------------------------------------
  ! Evaluaion of the Lanczos ground state + CI coefficients
  !--------------------------------------------------------
  allocate(C0(NDet))

  if (Model == 3) then
    call SetDelta_PJW(Delta)
    call LanczosGroundState(NS, NDet, ApplyHamiltonian_XXZPJW_OBC, LanczosMaxIt, LanczosTol, &
                          E0, iters, resid, C0)
  else if (Model == 4) then
    call SetDelta_PJW(Delta)
    call LanczosGroundState(NS, NDet, ApplyHamiltonian_XXZJW_P_OBC, LanczosMaxIt, LanczosTol, &
                          E0, iters, resid, C0)
  
  
  else
    call LanczosGroundState(NS, NDet, ApplyHamiltonian, LanczosMaxIt, LanczosTol, &
                          E0, iters, resid, C0)
  end if


  !call LanczosGroundState(NS, NDet, ApplyHamiltonian, LanczosMaxIt, LanczosTol, &
  !                        E0, iters, resid, C0)

  call CheckNormalization_FCI(C0, normsq, is_norm)
  write(*,*)
  write(*,*) "========================================"
  write(*,*) "Normalization check of FCI wavefunction"
  write(*,*) "========================================"
  write(*,' (A,F20.12)') "<Psi|Psi> = ", normsq
  write(*,*)
  write(*,'(A,F22.16)') "Ground-state energy = ", E0
  write(*,'(A,I0)')     "Lanczos iterations  = ", iters
  !write(*,'(A,ES12.4)') "Ritz diff (proxy)   = ", resid

  ! Optional: dump coefficients to a file
  
  call WriteTopCICoeffs("CI_top.dat", NS, C0, 30)
!!===============PRINTS CORRELATION FUNCTIONS F(r)===========================
!===============PRINTS CORRELATION FUNCTIONS F(r)===========================
if (isCorr /= 0) then

  allocate(SpSq(NS,NS))
  call Build_SpSq_Matrix_Full(NS, C0, NS, SpSq)

  write(*,*) "iPBC = ", iPBC

  if (Dim == 1) then

    allocate(Corr1D(NS))

    call CalcCorr1D_FromSpSq(SpSq, NS, iPBC, Corr1D)

    write(*,*) "========================================"
    write(*,*) "1D averaged correlation function Corr(r)"
    write(*,*) "========================================"
    write(*,*) "    r         Corr(r)"
    do p = 1, NS
      write(*,'(I6,2X,F20.12)') p-1, Corr1D(p)
    end do

    open(unit=11, file="Corr1D.dat", status="replace")
    do p = 1, NS
      write(11,'(I6,2X,F20.12)') p-1, Corr1D(p)
    end do
    close(11)

    deallocate(Corr1D)

  else

    allocate(Corr2D(Nx,Ny))

    call CalcCorr2D_FromSpSq(SpSq, Nx, Ny, iPBC, 1.0e-10_pr, Corr2D)

    write(*,*) "=================================================="
    write(*,*) "2D averaged correlation function Corr(dx,dy)"
    write(*,*) "=================================================="
    write(*,*) "   dx    dy         Corr(dx,dy)"
    do dx = 1, Nx
      do dy = 1, Ny
        write(*,'(I6,2X,I6,2X,F20.12)') dx-1, dy-1, Corr2D(dx,dy)
      end do
    end do

    open(unit=11, file="Corr2D.dat", status="replace")
    do dx = 1, Nx
      do dy = 1, Ny
        write(11,'(I6,2X,I6,2X,F20.12)') dx-1, dy-1, Corr2D(dx,dy)
      end do
      write(11,*)
    end do
    close(11)

    deallocate(Corr2D)

  end if

  deallocate(SpSq)

end if
!==============================================================
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

 