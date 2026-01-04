!======================================================================
!  Lanczos.f90
!
!  Purpose:
!  --------
!  Implements a matrix-free Lanczos algorithm to compute the
!  lowest (ground-state) eigenvalue of a large, sparse,
!  real-symmetric Hamiltonian.
!
!  This module is designed for Full Configuration Interaction (FCI)
!  calculations in SU(2) pairing / spin models, where the Hilbert space
!  dimension grows as 2^NSites and explicit dense diagonalization is
!  impossible for NSites ≳ 16–18.
!
!  The Hamiltonian is never stored explicitly. Instead, it is accessed
!  through a user-supplied matrix–vector product:
!
!        w = H * v
!
!  supplied via a procedure argument.
!
!
!  Mathematical background:
!  ------------------------
!  Given a symmetric matrix H and an initial normalized vector q₁,
!  the Lanczos algorithm builds an orthonormal Krylov basis:
!
!      K_m = span{ q₁, H q₁, H² q₁, ..., H^{m-1} q₁ }
!
!  and reduces H to a symmetric tridiagonal matrix Tₘ:
!
!          | α₁  β₁           |
!          | β₁  α₂  β₂       |
!      T = |      β₂  α₃  β₃  |
!          |           ...    |
!
!  whose eigenvalues (Ritz values) converge rapidly to the extremal
!  eigenvalues of H. The smallest Ritz value approximates the
!  ground-state energy.
!
!
!  Algorithmic details:
!  --------------------
!  - Classical (three-term) Lanczos recursion
!  - No full reorthogonalization (sufficient for ground state)
!  - Convergence monitored via change in lowest Ritz value
!  - Tridiagonal eigenproblem solved using LAPACK DSTEV
!
!
!  Usage:
!  ------
!  The main entry point is:
!
!    call LanczosGroundState(NS, NDet, MatVec, maxit, tol, E0)
!
!  where:
!    NS      = number of pair-levels / sites
!    NDet    = dimension of Fock space (= 2^NS)
!    MatVec  = user-supplied subroutine computing w = H*v
!    maxit   = maximum Lanczos iterations
!    tol     = convergence tolerance on Ritz value
!    E0      = returned ground-state energy
!
!
!  Design philosophy:
!  ------------------
!  - Completely Hamiltonian-agnostic
!  - Works for 1D / 2D XXZ, J1–J2, or general SU(2) pairing models
!  - Compatible with matrix-free operator application
!  - Suitable for exact diagonalization benchmarks and scaling studies
!
!
!  Notes:
!  ------
!  - The initial vector is random with its mean removed to avoid
!    trivial overlap with symmetry-protected states.
!  - Breakdown (β → 0) signals Krylov space closure.
!  - For very large systems, consider adding reorthogonalization.
!
!======================================================================
!=============================================================
!=============================================================
module LanczosSolver
  use Precision
  implicit none
  private
  public :: LanczosGroundState

  ! Abstract interface for H*v
  abstract interface
    subroutine MatVecProc(NS, v, w)
      import :: pr
      integer, intent(in) :: NS
      real(kind=pr), intent(in)  :: v(:)
      real(kind=pr), intent(out) :: w(:)
    end subroutine MatVecProc
  end interface

contains

  subroutine LanczosGroundState(NS, NDet, MatVec, maxit, tol, E0, iters, resid, C0)
    implicit none
    integer, intent(in) :: NS, NDet, maxit
    real(kind=pr), intent(in) :: tol
    procedure(MatVecProc) :: MatVec
    real(kind=pr), intent(out) :: E0
    integer, intent(out), optional :: iters
    real(kind=pr), intent(out), optional :: resid
    real(kind=pr), intent(out), optional :: C0(:)   ! <-- CI coefficients (length NDet)

    integer :: k, m, info
    real(kind=pr) :: alpha, beta, beta_prev
    real(kind=pr) :: Eold, Enew, diff, normq
    real(kind=pr), allocatable :: q(:), qprev(:), w(:), q1(:)
    real(kind=pr), allocatable :: d(:), e(:), s(:)
    logical :: need_vec

    if (NDet <= 0) stop "LanczosGroundState: NDet must be positive"
    if (maxit < 2) stop "LanczosGroundState: maxit must be >= 2"

    need_vec = present(C0)
    if (need_vec) then
      if (size(C0) /= NDet) stop "LanczosGroundState: C0 has wrong length"
    end if

    allocate(q(NDet), qprev(NDet), w(NDet))
    allocate(d(maxit), e(maxit))

    !-----------------------
    ! Initial vector q1
    !-----------------------
    call random_seed()
    call random_number(q)
    q = q - sum(q)/real(NDet,kind=pr)

    normq = sqrt(dot_product(q,q))
    if (normq == 0.0_pr) then
      q = 0.0_pr
      q(1) = 1.0_pr
      normq = 1.0_pr
    end if
    q = q / normq

    allocate(q1(NDet))
    q1 = q

    qprev = 0.0_pr
    beta_prev = 0.0_pr

    Eold = huge(1.0_pr)
    Enew = huge(1.0_pr)
    diff = huge(1.0_pr)
    m = 0

    !=========================================================
    ! PASS 1: Lanczos to get d(1:m), e(1:m-1), and ground E0
    !=========================================================
    do k = 1, maxit

      call MatVec(NS, q, w)

      alpha = dot_product(q, w)
      w = w - alpha*q - beta_prev*qprev

      beta = sqrt(dot_product(w, w))

      d(k) = alpha
      if (k < maxit) e(k) = beta

      m = k
      Enew = LowestEigen_Tridiag(d, e, m, info)
      if (info /= 0) then
        write(*,*) "Lanczos: tridiagonal eigensolve failed at k=", k, " info=", info
        stop
      end if

      diff = abs(Enew - Eold)
      if (k > 2) then
        if (diff < tol) exit
      end if
      Eold = Enew

      if (beta < 1.0e-30_pr) exit

      qprev = q
      q = w / beta
      beta_prev = beta

    end do

    E0 = Enew
    if (present(iters)) iters = m
    if (present(resid)) resid = diff

    !=========================================================
    ! If CI vector requested: compute Ritz vector s (T_m eigvec)
    ! then PASS 2: reconstruct C0 = sum_k s_k q_k
    !=========================================================
    if (need_vec) then
      allocate(s(m))
      call LowestEigenVec_Tridiag(d, e, m, Enew, s, info)
      if (info /= 0) then
        write(*,*) "Lanczos: tridiagonal eigvec solve failed, info=", info
        stop
      end if

      C0 = 0.0_pr

      ! Re-run Lanczos recurrence, but use stored d/e (no convergence checks)
      q = q1
      qprev = 0.0_pr
      beta_prev = 0.0_pr

      do k = 1, m
        C0 = C0 + s(k) * q

        call MatVec(NS, q, w)
        alpha = d(k)
        w = w - alpha*q - beta_prev*qprev

        if (k < m) then
          beta = e(k)
          if (beta < 1.0e-30_pr) exit
          qprev = q
          q = w / beta
          beta_prev = beta
        end if
      end do

      ! Normalize CI vector (good practice)
      normq = sqrt(dot_product(C0, C0))
      if (normq > 0.0_pr) C0 = C0 / normq

      deallocate(s)
    end if

    deallocate(q, qprev, w, q1, d, e)
  end subroutine LanczosGroundState


  function LowestEigen_Tridiag(d_in, e_in, m, info) result(Emin)
    implicit none
    integer, intent(in) :: m
    real(kind=pr), intent(in) :: d_in(:), e_in(:)
    integer, intent(out) :: info
    real(kind=pr) :: Emin

    real(kind=pr), allocatable :: d(:), e(:), work(:)
    integer :: lwork
    real(kind=pr) :: Ztmp(1,1)

    interface
      subroutine dstev(jobz, n, d, e, z, ldz, work, info)
        character(len=1), intent(in) :: jobz
        integer, intent(in) :: n, ldz
        double precision, intent(inout) :: d(*), e(*)
        double precision, intent(out)   :: z(ldz,*)
        double precision, intent(out)   :: work(*)
        integer, intent(out) :: info
      end subroutine dstev
    end interface

    if (m < 1) then
      info = -1
      Emin = huge(1.0_pr)
      return
    end if

    allocate(d(m))
    d = d_in(1:m)

    if (m == 1) then
      info = 0
      Emin = d(1)
      deallocate(d)
      return
    end if

    allocate(e(m-1))
    e = e_in(1:m-1)

    lwork = max(1, 2*m - 2)
    allocate(work(lwork))

    Ztmp = 0.0_pr
    call dstev('N', m, d, e, Ztmp, 1, work, info)

    if (info == 0) then
      Emin = d(1)
    else
      Emin = huge(1.0_pr)
    end if

    deallocate(d, e, work)
  end function LowestEigen_Tridiag


  subroutine LowestEigenVec_Tridiag(d_in, e_in, m, Emin, s, info)
    implicit none
    integer, intent(in) :: m
    real(kind=pr), intent(in) :: d_in(:), e_in(:)
    real(kind=pr), intent(out) :: Emin
    real(kind=pr), intent(out) :: s(:)     ! length m
    integer, intent(out) :: info

    real(kind=pr), allocatable :: d(:), e(:), work(:), Z(:,:)
    integer :: lwork

    interface
      subroutine dstev(jobz, n, d, e, z, ldz, work, info)
        character(len=1), intent(in) :: jobz
        integer, intent(in) :: n, ldz
        double precision, intent(inout) :: d(*), e(*)
        double precision, intent(out)   :: z(ldz,*)
        double precision, intent(out)   :: work(*)
        integer, intent(out) :: info
      end subroutine dstev
    end interface

    if (m < 1) then
      info = -1
      Emin = huge(1.0_pr)
      s = 0.0_pr
      return
    end if
    if (size(s) /= m) stop "LowestEigenVec_Tridiag: s has wrong length"

    allocate(d(m)); d = d_in(1:m)

    if (m == 1) then
      info = 0
      Emin = d(1)
      s(1) = 1.0_pr
      deallocate(d)
      return
    end if

    allocate(e(m-1)); e = e_in(1:m-1)

    allocate(Z(m,m))
    lwork = max(1, 2*m - 2)
    allocate(work(lwork))

    call dstev('V', m, d, e, Z, m, work, info)

    if (info == 0) then
      Emin = d(1)
      s(1:m) = Z(1:m,1)   ! eigenvector for smallest eigenvalue
    else
      Emin = huge(1.0_pr)
      s = 0.0_pr
    end if

    deallocate(d, e, Z, work)
  end subroutine LowestEigenVec_Tridiag



   
end module LanczosSolver
