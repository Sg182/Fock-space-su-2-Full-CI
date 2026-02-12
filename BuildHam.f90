module BuildHamMatVec

!Uses Sparse diagonalization method Lanczos
  
  use Precision
  use Constants
  use Integrals
  implicit none
  private
  public :: ApplyHamiltonian,ApplyHamiltonian_XXZPJW_OBC, SetSzSector, ClearSzSector, ApplyHamiltonian_SzSector
  public :: ApplyHamiltonian_XXZJW_P_OBC, GetSzSectorDim,Expect_Sz_Sector, GetSzSectorBasis , SetDelta_PJW

  integer, allocatable :: basis_g(:)
  integer, allocatable :: idx_g(:)     ! idx_g(0:2^NS-1) -> sector index (1..Nsect), 0 if not in sector
  integer :: NS_g = -1

  real(kind=pr) :: Delta_PJW = 0.0_pr
contains

  !=============================================================
  ! NEW: setter for Delta (call from Main before Lanczos)
  !=============================================================
  subroutine SetDelta_PJW(Delta)
    real(kind=pr), intent(in) :: Delta
    Delta_PJW = Delta
  end subroutine SetDelta_PJW

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

  pure real(kind=pr) function SzVal_LGhost(state, bit) result(sz)
    integer, intent(in) :: state, bit
    if (bit < 0) then
      sz = Half          ! Sz_0 = 1/2
    else
      sz = SzVal(state, bit)
    end if
  end function SzVal_LGhost



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
    ! BUILD BASIS FOR A FIXED Sz SECTOR.
    ! CONVENTION: BIT=1 => Sz=+1/2 (SPIN-UP), BIT=0 => Sz=-1/2 (SPIN-DOWN)
    
    ! Then SzTot = Nup - NS/2  =>  Nup = SzTot + NS/2
    integer, intent(in) :: NS, SzTot    !NS IS THE NUMBER OF SITES
    integer :: NFull, mu, nsect, c, NupTarget

    if (mod(NS,2) /= 0) then
      stop "SetSzSector: NS must be even for integer SzTot sectors"
    end if

    NupTarget = (NS/2) + SzTot   !NupTarget IS THE NUMBER OF UP-SPINS
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

  
  subroutine ApplyHamiltonian_XXZJW_P_OBC(NS, v, w)
  !THIS IS FOR 1D OBC ONLY! JW FOLLOWED BY PARITY
    integer, intent(in) :: NS
    real(kind=pr), intent(in)  :: v(:)
    real(kind=pr), intent(out) :: w(:)

    integer :: NDet
    integer :: mu, outState
    integer :: p
    integer :: pbit, leftbit, rightbit
    real(kind=pr) :: diag, szL, szR, szLR, cflip

    NDet = ishft(1, NS)
    if (size(v) /= NDet .or. size(w) /= NDet) stop "ApplyHamiltonian_XXZMAP_OBC: wrong vector size"

    w = Zero

    do mu = 0, NDet-1

      diag = Zero

    ! OBC: p = 1..NS-1
      do p = 1, NS-1
        pbit     = p - 1      ! site p
        leftbit  = p - 2      ! site p-1  (p=1 -> bit=-1 -> ghost Sz0)
        rightbit = p          ! site p+1  (max p=NS-1 -> bit=NS-1 OK)

        szL  = SzVal_LGhost(mu, leftbit) !Sz_{p-1}
        szR  = SzVal(mu, rightbit)       !Sz_{p+1}
        szLR = szL * szR

      ! Diagonal: + Delta * Sz_{p-1} Sz_{p+1}
        diag = diag + Delta_PJW * szLR

      ! Off-diagonal flip at site p:
      ! (1/2)Sx_p  ->  +1/4
      ! -2 Sz_{p-1} Sx_p Sz_{p+1} ->  -(szL*szR)
        cflip = 0.25_pr - szLR

        outState = ieor(mu, ishft(1, pbit))
        w(outState+1) = w(outState+1) + cflip * v(mu+1)
      end do

      w(mu+1) = w(mu+1) + diag * v(mu+1)

    end do
  end subroutine ApplyHamiltonian_XXZJW_P_OBC


  subroutine ApplyHamiltonian_XXZPJW_OBC(NS, v, w)
  ! w = H v in the full Fock space (dimension 2^NS)
  ! Hamiltonian (OBC, your convention):
  !   H = sum_{p=0..NS-2} [ Sx_p Sx_{p+2} + 2 Sx_p Sz_{p+1} Sx_{p+2} - (Delta/2) Sz_{p+1} ]
  ! with ghost Sx_{p+2}=1/2 when p=NS-2 (i.e., p+2 = NS out of range).
  integer, intent(in) :: NS
  real(kind=pr), intent(in)  :: v(:)
  real(kind=pr), intent(out) :: w(:)

  integer :: NDet
  integer :: mu, outState
  integer :: p
  integer :: pbit, midbit, p2bit
  integer :: mask
  real(kind=pr) :: diag, szmid, coeff

  NDet = ishft(1, NS)
  if (size(v) /= NDet .or. size(w) /= NDet) stop "ApplyHamiltonian_XXZPJW_OBC: wrong vector size"

  w = Zero

  do mu = 0, NDet-1

    !---------------------------
    ! Diagonal part:  -(Delta/2) * sum_{p=0..NS-2} Sz_{p+1}
    !---------------------------
    diag = 0.0
    do p = 1, NS-1
      midbit = p              ! (p+1) site has bit index p
      szmid  = SzVal(mu, midbit)
      diag   = diag - Half * Delta_PJW * szmid
    end do
    w(mu+1) = w(mu+1) + diag * v(mu+1)

    !---------------------------
    ! Off-diagonal:
    !   sum_{p=0..NS-2} [ Sx_p Sx_{p+2} + 2 Sx_p Sz_{p+1} Sx_{p+2} ]
    !
    ! In Sz basis: Sx flips a spin with matrix element 1/2.
    ! So:
    !   <out| Sx_p Sx_{p+2} |mu> = 1/4 (two flips)
    !   <out| 2 Sx_p Sz_{p+1} Sx_{p+2} |mu> = (1/2)*Sz_{p+1} (two flips)
    ! Total coeff = 1/4 + (1/2)*Sz_{p+1}
    !
    ! Boundary p=NS-2: take Sx_{p+2}=1/2 (ghost) => only flip p,
    ! and the same coeff formula still works with "one-flip" outState.
    !---------------------------
    do p = 1, NS-1
      pbit   = p - 1
      midbit = p
      szmid  = SzVal(mu, midbit)
      coeff  = 0.25_pr + Half * szmid

      if (p == NS-1) then
        ! boundary: flip only site p (since p+2 is ghost)
        outState = ieor(mu, ishft(1, pbit))
      else
        ! bulk: flip sites p and p+2
        p2bit = p + 1
        mask  = ior(ishft(1, pbit), ishft(1, p2bit))
        outState = ieor(mu, mask)
      end if

      w(outState+1) = w(outState+1) + coeff * v(mu+1)
    end do

  end do

end subroutine ApplyHamiltonian_XXZPJW_OBC




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
