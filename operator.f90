program Build_Pdag
  implicit none

  integer, parameter :: NLevels = 2   !Number of sites
  integer :: NDet, p
  INTEGER :: Ialloc
  integer, allocatable :: Pdag(:,:)

  NDet = ishft(1, NLevels)     ! 2^NLevels
  allocate(Pdag(NDet, NDet),Stat=IAlloc)

  If(IAlloc /= 0) then
    print*, "Could not allocate in Pdag matrix"
      stop 1
  EndIf
  Pdag = 0

  ! CHOOSE which level you mean by Pdag_p.
  ! IMPORTANT: p is 0-based level index:
  ! p=0 acts on the RIGHTMOST bit in the printed string
  ! p=NLevels-1 acts on the LEFTMOST bit
  p = 1

  call BuildPdagMatrix(NLevels, p, Pdag)

  print *, "==========Mapping for Pdag_p with p=", p, " (0-based) ===================="
  call PrintPdagMapping(NLevels, p)

  print *
  print *, "========== Dense matrix Pdag(:, :) in basis state=0..2^NLevels-1 ===================="
  call PrintIntMat(Pdag, NDet)

  deallocate(Pdag)

contains

  subroutine ApplyPdag(state_in, p, state_out, ok)
    implicit none
    integer, intent(in)  :: state_in, p
    integer, intent(out) :: state_out
    logical, intent(out) :: ok

    if (iand(state_in, ishft(1, p)) /= 0) then
      ok = .false.
      state_out = state_in
    else
      ok = .true.
      state_out = ior(state_in, ishft(1, p))
    end if
  end subroutine ApplyPdag


  subroutine BuildPdagMatrix(NLevels,p, Pdag)
    implicit none
    integer, intent(in) :: NLevels, p
    integer, intent(out) :: Pdag(:,:)
    integer :: NDet, mu, nu, outState
    logical :: ok

    NDet = ishft(1, NLevels)
    if (size(Pdag,1) /= NDet .or. size(Pdag,2) /= NDet) then
      stop "BuildPdagMatrix: wrong matrix size"
    end if

    Pdag = 0

    ! Column mu corresponds to |mu-1>, row nu corresponds to |nu-1>
    do mu = 1, NDet
      call ApplyPdag(mu-1, p, outState, ok)
      if (ok) then
        nu = outState + 1
        Pdag(nu, mu) = 1
      end if
    end do
  end subroutine BuildPdagMatrix


  subroutine PrintStateBits(state, NLevels)
    implicit none
    integer, intent(in) :: state, NLevels
    integer :: q

    ! Print as |b_{NLevels-1} ... b_1 b_0>
    do q = NLevels-1, 0, -1
      if (iand(state, ishft(1, q)) /= 0) then
        write(*,'(A)', advance='no') '1'
      else
        write(*,'(A)', advance='no') '0'
      end if
    end do
  end subroutine PrintStateBits


  subroutine PrintPdagMapping(NLevels, p)
    implicit none
    integer, intent(in) :: NLevels, p
    integer :: state, outState, NDet
    logical :: ok

    NDet = ishft(1, NLevels)

    do state = 0, NDet-1
      write(*,'(A)', advance='no') '|'
      call PrintStateBits(state, NLevels)
      write(*,'(A)', advance='no') '>  ->  '

      call ApplyPdag(state, p, outState, ok)
      if (.not. ok) then
        write(*,'(A)') '0'
      else
        write(*,'(A)', advance='no') '|'
        call PrintStateBits(outState, NLevels)
        write(*,'(A)') '>'
      end if
    end do
  end subroutine PrintPdagMapping


  subroutine PrintIntMat(A, n)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: A(n,n)
    integer :: i, j

    do i = 1, n
      write(*,'(*(I2,1X))') (A(i,j), j=1,n)
    end do
  end subroutine PrintIntMat

end program Build_Pdag
