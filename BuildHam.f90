module BuildHamMatVec

!Uses Sparse diagonalization method Lanczos
  
  use Precision
  use Constants
  use Integrals
  implicit none
  private
  public :: ApplyHamiltonian, SetSzSector, ClearSzSector, ApplyHamiltonian_SzSector, GetSzSectorDim,Expect_Sz_Sector, GetSzSectorBasis

  integer, allocatable :: basis_g(:)
  integer, allocatable :: idx_g(:)     ! idx_g(0:2^NS-1) -> sector index (1..Nsect), 0 if not in sector
  integer :: NS_g = -1


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

    subroutine SetSzSector(NS, SzTot)
    ! Build basis for a fixed total Sz sector.
    ! Convention: bit=1 => Sz=+1/2, bit=0 => Sz=-1/2
    ! Then SzTot = Nup - NS/2  =>  Nup = SzTot + NS/2
    integer, intent(in) :: NS, SzTot
    integer :: NFull, mu, nsect, c, NupTarget

    if (mod(NS,2) /= 0) then
      stop "SetSzSector: NS must be even for integer SzTot sectors"
    end if

    NupTarget = (NS/2) + SzTot
    if (NupTarget < 0 .or. NupTarget > NS) stop "SetSzSector: invalid SzTot for given NS"

    NFull = ishft(1, NS)

    if (allocated(basis_g)) deallocate(basis_g)
    if (allocated(idx_g))   deallocate(idx_g)

    ! Count
    nsect = 0
    do mu = 0, NFull-1
      if (popcnt(mu) == NupTarget) nsect = nsect + 1
    end do

    allocate(basis_g(nsect))
    allocate(idx_g(0:NFull-1))
    idx_g = 0
    NS_g = NS

    c = 0
    do mu = 0, NFull-1
      if (popcnt(mu) == NupTarget) then
        c = c + 1
        basis_g(c) = mu
        idx_g(mu) = c
      end if
    end do
  end subroutine SetSzSector


  subroutine ClearSzSector()
    if (allocated(basis_g)) deallocate(basis_g)
    if (allocated(idx_g))   deallocate(idx_g)
    NS_g = -1
  end subroutine ClearSzSector


  integer function GetSzSectorDim() result(nsect)
    if (.not. allocated(basis_g)) then
      nsect = 0
    else
      nsect = size(basis_g)
    end if
  end function GetSzSectorDim


  subroutine ApplyHamiltonian_SzSector(NS, v, w)
    ! w = H v restricted to the CURRENT Sz sector defined by SetSzSector
    integer, intent(in) :: NS
    real(kind=pr), intent(in)  :: v(:)
    real(kind=pr), intent(out) :: w(:)

    integer :: NFull, Nsect
    integer :: i, j, mu, outState
    integer :: p, q, pbit, qbit
    real(kind=pr) :: diag, szp, szq, hij
    logical :: ok

    if (NS_g /= NS) stop "ApplyHamiltonian_SzSector: sector not set for this NS"
    if (.not. allocated(basis_g)) stop "ApplyHamiltonian_SzSector: call SetSzSector first"
    if (.not. allocated(idx_g))   stop "ApplyHamiltonian_SzSector: call SetSzSector first"

    NFull  = ishft(1, NS)
    Nsect  = size(basis_g)

    if (size(v) /= Nsect .or. size(w) /= Nsect) stop "ApplyHamiltonian_SzSector: wrong vector size"
    if (lbound(idx_g,1) /= 0 .or. ubound(idx_g,1) /= NFull-1) stop "ApplyHamiltonian_SzSector: idx bounds wrong"

    w = Zero

    do i = 1, Nsect
      mu = basis_g(i)

      !---------------------------
      ! Diagonal part
      !---------------------------
      diag = H000

      do p = 1, NS
        if (H010(p) /= Zero) then
          szp = SzVal(mu, p-1)
          diag = diag + H010(p) * szp
        end if
      end do

      do p = 2, NS
        szp = SzVal(mu, p-1)
        do q = 1, p-1
          if (H020(p,q) /= Zero) then
            szq = SzVal(mu, q-1)
            diag = diag + H020(p,q) * (szp * szq)
          end if
        end do
      end do

      w(i) = w(i) + diag * v(i)

      !---------------------------
      ! Off-diagonal
      !---------------------------
      do p = 1, NS
        pbit = p - 1
        do q = 1, NS
          hij = H101(p,q)
          if (hij == Zero) cycle
          qbit = q - 1

          call ApplySpSm(mu, pbit, qbit, outState, ok)
          if (ok) then
            j = idx_g(outState)
            if (j > 0) w(j) = w(j) + hij * v(i)
            ! j should always be > 0 for Sz-conserving H, but keep as safety.
          end if
        end do
      end do

    end do

  end subroutine ApplyHamiltonian_SzSector



  !====================================================================

   !EXPECTATION OF Sz for DIFFERENT Sz sectors like Sz=1,-1,.......

  real(kind=pr) function Expect_Sz_Sector(NS, C, pbit) result(val)
  integer, intent(in) :: NS, pbit
  real(kind=pr), intent(in) :: C(:)

  integer :: i, mu, Nsect
  real(kind=pr) :: w

  if (NS_g /= NS) stop "Expect_Sz_Sector: sector not set for this NS"
  if (.not. allocated(basis_g)) stop "Expect_Sz_Sector: call SetSzSector first"

  Nsect = size(basis_g)
  if (size(C) /= Nsect) stop "Expect_Sz_Sector: C has wrong length"

  val = 0.0_pr
  do i = 1, Nsect
    mu = basis_g(i)
    w  = C(i) * C(i)          ! |C_i|^2 (real)
    val = val + w * SzVal(mu, pbit)
  end do
end function Expect_Sz_Sector
!==========================================================================

subroutine GetSzSectorBasis(basis_out)
  implicit none
  integer, allocatable, intent(out) :: basis_out(:)

  if (.not. allocated(basis_g)) stop "GetSzSectorBasis: sector not set"
  allocate(basis_out(size(basis_g)))
  basis_out = basis_g
end subroutine GetSzSectorBasis


end module BuildHamMatVec
