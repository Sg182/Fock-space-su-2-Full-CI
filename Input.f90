module InputParams
  use Precision
  implicit none
  integer :: NLevels = 6
  real(kind=pr) :: Delta = -0.60_pr
  real(kind=pr) :: J2    = 0.0_pr
  logical :: Periodic    = .FALSE.
contains
  subroutine PrintInput()
    implicit none
    write(*,'(A,I0)')     "NLevels   = ", NLevels
    write(*,'(A,F8.4)')   "Delta     = ", Delta
    write(*,'(A,F6.4)')   "J2        = ", J2
    write(*,'(A,L1)')     "Periodic  = ", Periodic
  end subroutine PrintInput
end module InputParams
