module Integrals
  use Precision
  use Constants
  implicit none

  real(kind=pr) :: H000
  real(kind=pr), allocatable :: H010(:)
  real(kind=pr), allocatable :: H020(:,:)
  real(kind=pr), allocatable :: H101(:,:)

contains

  subroutine SetUpIntegrals(NLevels)
    integer, intent(in) :: NLevels
    integer :: IAlloc
    allocate(H010(NLevels), H020(NLevels,NLevels), H101(NLevels,NLevels), stat=IAlloc)
    if (IAlloc /= 0) stop "SetUpIntegrals: allocation failed"
    H000 = Zero; H010 = Zero; H020 = Zero; H101 = Zero
  end subroutine SetUpIntegrals

  subroutine ShutDownIntegrals()
    integer :: IAlloc
    if (allocated(H010)) deallocate(H010, stat=IAlloc)
    if (allocated(H020)) deallocate(H020, stat=IAlloc)
    if (allocated(H101)) deallocate(H101, stat=IAlloc)
  end subroutine ShutDownIntegrals

  subroutine DoIntegralsXXZ_1D(NLevels, Delta, Periodic)
    integer, intent(in) :: NLevels
    real(kind=pr), intent(in) :: Delta
    logical, intent(in) :: Periodic
    integer :: i, j

    H000 = Zero; H010 = Zero; H020 = Zero; H101 = Zero

    do i = 1, NLevels-1
      j = i + 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    end do

    if (Periodic) then
      i = NLevels; j = 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    end if
  end subroutine DoIntegralsXXZ_1D

  subroutine DoIntegralsJ1J2XXZ_1D(NLevels, Delta, J2, Periodic)
    integer, intent(in) :: NLevels
    real(kind=pr), intent(in) :: Delta, J2
    logical, intent(in) :: Periodic
    integer :: i, j

    H000 = Zero; H010 = Zero; H020 = Zero; H101 = Zero

    ! 1NN
    do i = 1, NLevels-1
      j = i + 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    end do
    if (Periodic) then
      i = NLevels; j = 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    end if

    ! 2NN
    do i = 1, NLevels-2
      j = i + 2
      H101(i,j) = H101(i,j) + Half*J2
      H101(j,i) = H101(j,i) + Half*J2
      H020(i,j) = H020(i,j) + Delta*J2
      H020(j,i) = H020(j,i) + Delta*J2
    end do
    if (Periodic) then
      i = NLevels;   j = 2
      H101(i,j) = H101(i,j) + Half*J2
      H101(j,i) = H101(j,i) + Half*J2
      H020(i,j) = H020(i,j) + Delta*J2
      H020(j,i) = H020(j,i) + Delta*J2

      i = NLevels-1; j = 1
      H101(i,j) = H101(i,j) + Half*J2
      H101(j,i) = H101(j,i) + Half*J2
      H020(i,j) = H020(i,j) + Delta*J2
      H020(j,i) = H020(j,i) + Delta*J2
    end if
  end subroutine DoIntegralsJ1J2XXZ_1D

end module Integrals
