module BuildHamDense
  use Precision
  use Constants
  use Operators
  use Integrals
  implicit none
contains

  subroutine BuildHamiltonianDense(NLevels, HMat)
    integer, intent(in) :: NLevels
    real(kind=pr), intent(out) :: HMat(:,:)
    integer :: NDet, p, q, IAlloc
    real(kind=pr), allocatable :: Mat1(:,:), Mat2(:,:)

    NDet = ishft(1, NLevels)
    if (size(HMat,1) /= NDet .or. size(HMat,2) /= NDet) stop "BuildHamiltonianDense: wrong size"
    HMat = Zero

    allocate(Mat1(NDet,NDet), Mat2(NDet,NDet), stat=IAlloc)
    if (IAlloc /= 0) stop "BuildHamiltonianDense: allocation failed"

    ! Constant term
    if (abs(H000) > 0.0_pr) then
      do p = 1, NDet
        HMat(p,p) = HMat(p,p) + H000
      end do
    end if

    ! Sum_p H010(p) Sz(p)
    do p = 1, NLevels
      if (abs(H010(p)) > 0.0_pr) HMat = HMat + H010(p) * Sz(:,:,p)
    end do

    ! Sum_{p>q} H020(p,q) Sz(p)Sz(q) (works since we stored both [p,q] and [q,p])
    do p = 1, NLevels
      Mat1 = Sz(:,:,p)
      do q = 1, p-1
        if (abs(H020(p,q)) > 0.0_pr) then
          Mat2 = Sz(:,:,q)
          HMat = HMat + H020(p,q) * matmul(Mat1, Mat2)
        end if
      end do
    end do

    ! Sum_{p,q} H101(p,q) S+(p)S-(q)
    do p = 1, NLevels
      Mat1 = Splus(:,:,p)
      do q = 1, NLevels
        if (abs(H101(p,q)) > 0.0_pr) then
          Mat2 = Sminus(:,:,q)
          HMat = HMat + H101(p,q) * matmul(Mat1, Mat2)
        end if
      end do
    end do

    deallocate(Mat1, Mat2, stat=IAlloc)
    if (IAlloc /= 0) stop "BuildHamiltonianDense: deallocation failed"
  end subroutine BuildHamiltonianDense

end module BuildHamDense
