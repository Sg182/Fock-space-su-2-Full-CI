!To execute the program , run the makefile by the 'make' command, and 
!running the executable 

program Fock_FCI
  use Precision
  use Constants
  use InputParams
  use Operators
  use Integrals
  use BuildHamDense
  use LapackIface
  implicit none

  integer :: NDet, IAlloc, info, lwork
  real(kind=pr), allocatable :: H(:,:), evals(:), work(:)

  call PrintInput()

  NDet = ishft(1, NLevels)
  write(*,'(A,I0)') "NDet      = ", NDet

  call SetUpPairOperators(NLevels)

  call SetUpIntegrals(NLevels)
  if (abs(J2) < 1.0e-14_pr) then
    call DoIntegralsXXZ_1D(NLevels, Delta, Periodic)
  else
    call DoIntegralsJ1J2XXZ_1D(NLevels, Delta, J2, Periodic)
  end if

  allocate(H(NDet,NDet), evals(NDet), stat=IAlloc)
  if (IAlloc /= 0) stop "Allocation failed for H/evals"

  call BuildHamiltonianDense(NLevels, H)

  lwork = max(1, 3*NDet - 1)
  allocate(work(lwork), stat=IAlloc)
  if (IAlloc /= 0) stop "Allocation failed for work"

  call dsyev('V','U', NDet, H, NDet, evals, work, lwork, info)
  if (info /= 0) then
    write(*,*) "DSYEV failed, info=", info
    stop
  end if

  write(*,'(A,F20.14)') "Ground-state energy = ", evals(1)

  deallocate(work, H, evals, stat=IAlloc)
  call ShutDownIntegrals()
  call ShutDownPairOperators()

end program Fock_FCI
