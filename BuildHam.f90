module BuildHamMatVec

!Uses Sparse diagonalization method Lanczos
  
  use Precision
  use Constants
  use Integrals
  implicit none
  private
  public :: ApplyHamiltonian

contains

  pure real(kind=pr) function SzVal(state, bit) result(sz)
    ! state is 0-based integer determinant bitstring
    ! bit is 0-based bit index
    integer, intent(in) :: state, bit
    if (iand(state, ishft(1, bit)) /= 0) then
      sz = Half
    else
      sz = -Half
    end if
  end function SzVal


  pure subroutine ApplySpSm(state_in, pbit, qbit, state_out, ok)
    ! Applies S^+_p S^-_q to |state_in>
    ! pbit,qbit are 0-based bit indices
    integer, intent(in)  :: state_in, pbit, qbit
    integer, intent(out) :: state_out
    logical, intent(out) :: ok
    integer :: s1

    ! S^-_q requires q occupied
    if (iand(state_in, ishft(1, qbit)) == 0) then
      ok = .false.
      state_out = state_in
      return
    end if

    ! First apply S^-_q : clear qbit
    s1 = iand(state_in, not(ishft(1, qbit)))

    if (pbit == qbit) then
      ! Then S^+_p restores it (always allowed because we just cleared it)
      ok = .true.
      state_out = ior(s1, ishft(1, pbit))   ! equals original state
      return
    end if

    ! S^+_p requires p empty AFTER the lowering
    if (iand(s1, ishft(1, pbit)) /= 0) then
      ok = .false.
      state_out = state_in
      return
    end if

    ok = .true.
    state_out = ior(s1, ishft(1, pbit))
  end subroutine ApplySpSm


  subroutine ApplyHamiltonian(NS, v, w)
    ! w = H v  in the full Fock space (dimension 2^NS)
    integer, intent(in) :: NS
    real(kind=pr), intent(in)  :: v(:)
    real(kind=pr), intent(out) :: w(:)

    integer :: NDet
    integer :: mu, outState
    integer :: p, q
    integer :: pbit, qbit
    real(kind=pr) :: diag, szp, szq, hij
    logical :: ok

    NDet = ishft(1, NS)
    if (size(v) /= NDet .or. size(w) /= NDet) stop "ApplyHamiltonian: wrong vector size"

    w = Zero

    do mu = 0, NDet-1

      !---------------------------
      ! Diagonal part
      !---------------------------
      diag = H000

      ! Sum_p H010(p) Sz(p)
      do p = 1, NS
        if (H010(p) /= Zero) then
          szp = SzVal(mu, p-1)
          diag = diag + H010(p) * szp
        end if
      end do

      ! Sum_{p>q} H020(p,q) Sz(p)Sz(q)
      do p = 2, NS
        szp = SzVal(mu, p-1)
        do q = 1, p-1
          if (H020(p,q) /= Zero) then
            szq = SzVal(mu, q-1)
            diag = diag + H020(p,q) * (szp * szq)
          end if
        end do
      end do

      ! apply diagonal
      w(mu+1) = w(mu+1) + diag * v(mu+1)

      !---------------------------
      ! Off-diagonal: Sum_{p,q} H101(p,q) S+(p) S-(q)
      !---------------------------
      do p = 1, NS
        pbit = p - 1
        do q = 1, NS
          hij = H101(p,q)
          if (hij == Zero) cycle
          qbit = q - 1

          call ApplySpSm(mu, pbit, qbit, outState, ok)
          if (ok) then
            w(outState+1) = w(outState+1) + hij * v(mu+1)
          end if
        end do
      end do

    end do

  end subroutine ApplyHamiltonian

end module BuildHamMatVec
