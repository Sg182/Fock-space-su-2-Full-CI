program Operators
  implicit none

  integer, parameter :: NLevels = 2
  integer :: NDet, p
  integer, allocatable :: Pdag(:,:), Pmat(:,:), Nmat(:,:)

  NDet = ishft(1, NLevels)     ! 2^NLevels

  allocate(Pdag(NDet,NDet), Pmat(NDet,NDet), Nmat(NDet,NDet))
  Pdag = 0
  Pmat = 0
  Nmat = 0

  ! Choose level (0-based):
  ! p=0 is the rightmost bit in |b_{NLevels-1}...b_0>
  p = 1

  call BuildPdagMatrix(NLevels, p, Pdag)
  call BuildPMatrix   (NLevels, p, Pmat)
  call BuildNMatrix   (NLevels, p, Nmat)

  print *,"====== Mapping for Pdag_p with p=", p, " (0-based) ======"
  call PrintOpMapping(NLevels, p, "Pdag")

  print *
  print *, "==== Mapping for P_p with p=", p, " (0-based) ===="
  call PrintOpMapping(NLevels, p, "P")

  print *
  print *, "==== Mapping for N_p with p=", p, " (0-based) ===="
  call PrintOpMapping(NLevels, p, "N")

  print *
  print *, "---- Dense matrix Pdag ----"
  call PrintIntMat(Pdag, NDet)

  print *
  print *, "---- Dense matrix P ----"
  call PrintIntMat(Pmat, NDet)

  print *
  print *, "---- Dense matrix N ----"
  call PrintIntMat(Nmat, NDet)

  deallocate(Pdag, Pmat, Nmat)

contains

  !========================
  ! Operator actions
  !========================
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

  subroutine ApplyP(state_in, p, state_out, ok)
    implicit none
    integer, intent(in)  :: state_in, p
    integer, intent(out) :: state_out
    logical, intent(out) :: ok

    if (iand(state_in, ishft(1, p)) == 0) then
      ok = .false.
      state_out = state_in
    else
      ok = .true.
      state_out = iand(state_in, not(ishft(1, p)))
    end if
  end subroutine ApplyP

  subroutine ApplyN(state_in, p, value)
    implicit none
    integer, intent(in)  :: state_in, p
    integer, intent(out) :: value
    ! Pair-basis number operator: N_p = 2 if occupied, else 0
    if (iand(state_in, ishft(1, p)) /= 0) then
      value = 2
    else
      value = 0
    end if
  end subroutine ApplyN

  !========================
  ! Matrix builders
  !========================
  subroutine BuildPdagMatrix(NLevels, p, Op)
    implicit none
    integer, intent(in) :: NLevels, p
    integer, intent(out) :: Op(:,:)
    integer :: NDet, mu, nu, outState
    logical :: ok

    NDet = ishft(1, NLevels)
    if (size(Op,1) /= NDet .or. size(Op,2) /= NDet) stop "BuildPdagMatrix: wrong size"
    Op = 0

    do mu = 1, NDet
      call ApplyPdag(mu-1, p, outState, ok)
      if (ok) then
        nu = outState + 1
        Op(nu, mu) = 1
      end if
    end do
  end subroutine BuildPdagMatrix

  subroutine BuildPMatrix(NLevels, p, Op)
    implicit none
    integer, intent(in) :: NLevels, p
    integer, intent(out) :: Op(:,:)
    integer :: NDet, mu, nu, outState
    logical :: ok

    NDet = ishft(1, NLevels)
    if (size(Op,1) /= NDet .or. size(Op,2) /= NDet) stop "BuildPMatrix: wrong size"
    Op = 0

    do mu = 1, NDet
      call ApplyP(mu-1, p, outState, ok)
      if (ok) then
        nu = outState + 1
        Op(nu, mu) = 1
      end if
    end do
  end subroutine BuildPMatrix

  subroutine BuildNMatrix(NLevels, p, Op)
    implicit none
    integer, intent(in) :: NLevels, p
    integer, intent(out) :: Op(:,:)
    integer :: NDet, mu, val

    NDet = ishft(1, NLevels)
    if (size(Op,1) /= NDet .or. size(Op,2) /= NDet) stop "BuildNMatrix: wrong size"
    Op = 0

    do mu = 1, NDet
      call ApplyN(mu-1, p, val)
      Op(mu, mu) = val
    end do
  end subroutine BuildNMatrix

  !========================
  ! Printing helpers
  !========================
  subroutine PrintStateBits(state, NLevels)
    implicit none
    integer, intent(in) :: state, NLevels
    integer :: q
    do q = NLevels-1, 0, -1
      if (iand(state, ishft(1, q)) /= 0) then
        write(*,'(A)', advance='no') '1'
      else
        write(*,'(A)', advance='no') '0'
      end if
    end do
  end subroutine PrintStateBits

  subroutine PrintOpMapping(NLevels, p, which)
    implicit none
    integer, intent(in) :: NLevels, p
    character(*), intent(in) :: which
    integer :: state, outState, NDet, val
    logical :: ok

    NDet = ishft(1, NLevels)

    do state = 0, NDet-1
      write(*,'(A)', advance='no') '|'
      call PrintStateBits(state, NLevels)
      write(*,'(A)', advance='no') '>  ->  '

      select case (trim(which))
      case ("Pdag")
        call ApplyPdag(state, p, outState, ok)
        if (.not. ok) then
          write(*,'(A)') '0'
        else
          write(*,'(A)', advance='no') '|'
          call PrintStateBits(outState, NLevels)
          write(*,'(A)') '>'
        end if

      case ("P")
        call ApplyP(state, p, outState, ok)
        if (.not. ok) then
          write(*,'(A)') '0'
        else
          write(*,'(A)', advance='no') '|'
          call PrintStateBits(outState, NLevels)
          write(*,'(A)') '>'
        end if

      case ("N")
        call ApplyN(state, p, val)
        write(*,'(I0)') val

      case default
        stop "PrintOpMapping: unknown operator"
      end select
    end do
  end subroutine PrintOpMapping

  subroutine PrintIntMat(A, n)
    implicit none
    integer, intent(in) :: n
    integer, intent(in) :: A(n,n)
    integer :: i, j
    do i = 1, n
      write(*,'(*(I2,1X))') (A(i,j), j=1,n)
    end do
  end subroutine PrintIntMat

end program Operators
