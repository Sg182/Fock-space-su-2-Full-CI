!============================================================
! Dense FCI in SU(2) pair Fock space (N_p in {0,2})
! Build XXZ / J1-J2(XXZ) 1D Hamiltonian and diagonalize (DSYEV)
!============================================================

Module Precision
  Implicit None
  Integer, Parameter :: pr = selected_real_kind(15, 300)  ! ~double
End Module Precision

Module Constants
  Use Precision
  Implicit None
  Real(Kind=pr), Parameter :: Zero=0.0_pr, One=1.0_pr
  Real(Kind=pr), Parameter :: Half=0.5_pr
End Module Constants


Module OperatorsPair
  Use Precision
  Use Constants
  Implicit None

  ! Dense operator matrices in pair basis
  Real(Kind=pr), Allocatable :: Pdag(:,:,:), Pmat(:,:,:), Nmat(:,:,:)
  Real(Kind=pr), Allocatable :: Sz(:,:,:), Splus(:,:,:), Sminus(:,:,:)

Contains

  Pure Integer Function GetBit(bits, p) Result(b)
    Implicit None
    Integer, Intent(In) :: bits, p  ! p is 0-based
    b = IAND(ISHFT(bits, -p), 1)
  End Function GetBit

  Pure Integer Function SetBit(bits, p) Result(newbits)
    Implicit None
    Integer, Intent(In) :: bits, p
    newbits = IOR(bits, ISHFT(1, p))
  End Function SetBit

  Pure Integer Function ClearBit(bits, p) Result(newbits)
    Implicit None
    Integer, Intent(In) :: bits, p
    newbits = IAND(bits, NOT(ISHFT(1, p)))
  End Function ClearBit

  Pure Integer Function ApplyPdagBits(bits_in, p) Result(bits_out)
    ! returns -1 if annihilated (occupied already), else new bits
    Implicit None
    Integer, Intent(In) :: bits_in, p
    If (GetBit(bits_in, p) == 1) Then
      bits_out = -1
    Else
      bits_out = SetBit(bits_in, p)
    End If
  End Function ApplyPdagBits

  Pure Integer Function ApplyPBits(bits_in, p) Result(bits_out)
    ! returns -1 if annihilated (empty), else new bits
    Implicit None
    Integer, Intent(In) :: bits_in, p
    If (GetBit(bits_in, p) == 0) Then
      bits_out = -1
    Else
      bits_out = ClearBit(bits_in, p)
    End If
  End Function ApplyPBits

  Pure Integer Function ApplyNval(bits_in, p) Result(nval)
    ! N_p in {0,2}
    Implicit None
    Integer, Intent(In) :: bits_in, p
    If (GetBit(bits_in, p) == 1) Then
      nval = 2
    Else
      nval = 0
    End If
  End Function ApplyNval


  Subroutine SetUpPairOperators(NLevels)
    Implicit None
    Integer, Intent(In) :: NLevels
    Integer :: NDet, p, mu, nu_bits, mu_bits
    Integer :: IAlloc
    Real(Kind=pr), Allocatable :: Iden(:,:)

    NDet = ISHFT(1, NLevels)

    Allocate(Pdag(NDet,NDet,NLevels), Pmat(NDet,NDet,NLevels), Nmat(NDet,NDet,NLevels), &
             Sz(NDet,NDet,NLevels), Splus(NDet,NDet,NLevels), Sminus(NDet,NDet,NLevels), &
             Iden(NDet,NDet), Stat=IAlloc)
    If (IAlloc /= 0) Stop "SetUpPairOperators: allocation failed"

    Pdag = Zero
    Pmat = Zero
    Nmat = Zero
    Sz   = Zero
    Splus = Zero
    Sminus = Zero

    Iden = Zero
    Do mu = 1, NDet
      Iden(mu,mu) = One
    End Do

    ! Build Pdag, P, N (dense) for each level p (0-based in bit operations)
    Do p = 1, NLevels
      Do mu = 1, NDet
        mu_bits = mu - 1

        nu_bits = ApplyPdagBits(mu_bits, p-1)
        If (nu_bits >= 0) Then
          Pdag(nu_bits+1, mu, p) = One
        End If

        nu_bits = ApplyPBits(mu_bits, p-1)
        If (nu_bits >= 0) Then
          Pmat(nu_bits+1, mu, p) = One
        End If

        Nmat(mu, mu, p) = Real(ApplyNval(mu_bits, p-1), Kind=pr)  ! 0 or 2
      End Do

      ! Spin operators in this pair basis:
      ! S+ = Pdag, S- = P
      Splus(:,:,p)  = Pdag(:,:,p)
      Sminus(:,:,p) = Pmat(:,:,p)

      ! Sz = (N - 1)/2 = N/2 - 1/2
      Sz(:,:,p) = Half * Nmat(:,:,p) - Half * Iden
    End Do

    DeAllocate(Iden, Stat=IAlloc)
    If (IAlloc /= 0) Stop "SetUpPairOperators: deallocation failed"
  End Subroutine SetUpPairOperators


  Subroutine ShutDownPairOperators()
    Implicit None
    Integer :: IAlloc
    If (Allocated(Pdag))   DeAllocate(Pdag, Stat=IAlloc)
    If (Allocated(Pmat))   DeAllocate(Pmat, Stat=IAlloc)
    If (Allocated(Nmat))   DeAllocate(Nmat, Stat=IAlloc)
    If (Allocated(Sz))     DeAllocate(Sz, Stat=IAlloc)
    If (Allocated(Splus))  DeAllocate(Splus, Stat=IAlloc)
    If (Allocated(Sminus)) DeAllocate(Sminus, Stat=IAlloc)
  End Subroutine ShutDownPairOperators

End Module OperatorsPair


Module IntegralsXXZJ1J2_1D
  Use Precision
  Use Constants
  Implicit None

  Real(Kind=pr) :: H000
  Real(Kind=pr), Allocatable :: H010(:)      ! coefficient of Sz(p)
  Real(Kind=pr), Allocatable :: H020(:,:)    ! coefficient of Sz(p)Sz(q)
  Real(Kind=pr), Allocatable :: H101(:,:)    ! coefficient of S+(p)S-(q)

Contains

  Subroutine SetUpIntegrals(NLevels)
    Implicit None
    Integer, Intent(In) :: NLevels
    Integer :: IAlloc
    Allocate(H010(NLevels), H020(NLevels,NLevels), H101(NLevels,NLevels), Stat=IAlloc)
    If (IAlloc /= 0) Stop "SetUpIntegrals: allocation failed"
    H000 = Zero
    H010 = Zero
    H020 = Zero
    H101 = Zero
  End Subroutine SetUpIntegrals

  Subroutine ShutDownIntegrals()
    Implicit None
    Integer :: IAlloc
    If (Allocated(H010)) DeAllocate(H010, Stat=IAlloc)
    If (Allocated(H020)) DeAllocate(H020, Stat=IAlloc)
    If (Allocated(H101)) DeAllocate(H101, Stat=IAlloc)
  End Subroutine ShutDownIntegrals

  Subroutine DoIntegralsXXZ_1D(NLevels, Delta, Periodic)
    ! Hamiltonian (spin form):
    ! H = Sum_<i,j> [ 1/2 (S+_i S-_j + S-_i S+_j) + Delta Sz_i Sz_j ]
    Implicit None
    Integer, Intent(In) :: NLevels
    Real(Kind=pr), Intent(In) :: Delta
    Logical, Intent(In) :: Periodic
    Integer :: i, j, Nbonds

    H000 = Zero
    H010 = Zero
    H020 = Zero
    H101 = Zero

    ! bonds: (i,i+1) for i=1..NLevels-1, plus (N,1) if periodic
    Nbonds = NLevels - 1
    If (Periodic) Nbonds = Nbonds + 1

    ! constant from SzSz is automatically included if you build Sz matrices,
    ! so here H000 is just "extra constants" (none needed).
    ! (Your requested "+1/4 factor" is already inside SzSz since Sz=(N-1)/2.)

    Do i = 1, NLevels-1
      j = i + 1
      ! XY part: 1/2(S+_i S-_j + S+_j S-_i)
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half

      ! Ising part: Delta * Sz_i Sz_j
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    End Do

    If (Periodic) Then
      i = NLevels
      j = 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    End If

  End Subroutine DoIntegralsXXZ_1D

  Subroutine DoIntegralsJ1J2XXZ_1D(NLevels, Delta, J2, Periodic)
    ! J1-J2 XXZ:
    ! H = Sum_<i,j> [ 1/2(S+_i S-_j + S-_i S+_j) + Delta Sz_i Sz_j ]
    !   + Sum_<<i,k>> [ J2 * ( 1/2(S+_i S-_k + S-_i S+_k) + Delta Sz_i Sz_k ) ]
    ! Here J1 = 1, next-nearest have factor J2.
    Implicit None
    Integer, Intent(In) :: NLevels
    Real(Kind=pr), Intent(In) :: Delta, J2
    Logical, Intent(In) :: Periodic
    Integer :: i, j

    H000 = Zero
    H010 = Zero
    H020 = Zero
    H101 = Zero

    ! 1NN bonds (i,i+1)
    Do i = 1, NLevels-1
      j = i + 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    End Do
    If (Periodic) Then
      i = NLevels; j = 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    End If

    ! 2NN bonds (i,i+2) with factor J2
    Do i = 1, NLevels-2
      j = i + 2
      H101(i,j) = H101(i,j) + Half*J2
      H101(j,i) = H101(j,i) + Half*J2
      H020(i,j) = H020(i,j) + Delta*J2
      H020(j,i) = H020(j,i) + Delta*J2
    End Do
    If (Periodic) Then
      ! wrap-around 2NN: (N,2) and (N-1,1)
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
    End If
  End Subroutine DoIntegralsJ1J2XXZ_1D

End Module IntegralsXXZJ1J2_1D


Module BuildHamDense
  Use Precision
  Use Constants
  Use OperatorsPair
  Use IntegralsXXZJ1J2_1D
  Implicit None
Contains

  Subroutine BuildHamiltonianDense(NLevels, HMat)
    Implicit None
    Integer, Intent(In) :: NLevels
    Real(Kind=pr), Intent(Out) :: HMat(:,:)
    Integer :: NDet, p, q, IAlloc
    Real(Kind=pr), Allocatable :: Mat1(:,:), Mat2(:,:)

    NDet = ISHFT(1, NLevels)
    If (Size(HMat,1) /= NDet .or. Size(HMat,2) /= NDet) Stop "BuildHamiltonianDense: wrong HMat size"

    HMat = Zero

    Allocate(Mat1(NDet,NDet), Mat2(NDet,NDet), Stat=IAlloc)
    If (IAlloc /= 0) Stop "BuildHamiltonianDense: allocation failed"

    ! Constant
    If (Abs(H000) > 0.0_pr) Then
      Do p = 1, NDet
        HMat(p,p) = HMat(p,p) + H000
      End Do
    End If

    ! One-body Sz terms
    Do p = 1, NLevels
      If (Abs(H010(p)) > 0.0_pr) HMat = HMat + H010(p) * Sz(:,:,p)
    End Do

    ! Two-body Sz Sz terms (use p>q style with symmetry)
    Do p = 1, NLevels
      Mat1 = Sz(:,:,p)
      Do q = 1, p-1
        If (Abs(H020(p,q)) > 0.0_pr) Then
          Mat2 = Sz(:,:,q)
          ! Since H020(p,q) and H020(q,p) are both filled, we can just use one
          ! and do MatMul once; safest is:
          HMat = HMat + H020(p,q) * MatMul(Mat1, Mat2)
        End If
      End Do
    End Do

    ! XY terms: H101(p,q) S+(p) S-(q)
    Do p = 1, NLevels
      Mat1 = Splus(:,:,p)
      Do q = 1, NLevels
        If (Abs(H101(p,q)) > 0.0_pr) Then
          Mat2 = Sminus(:,:,q)
          HMat = HMat + H101(p,q) * MatMul(Mat1, Mat2)
        End If
      End Do
    End Do

    DeAllocate(Mat1, Mat2, Stat=IAlloc)
    If (IAlloc /= 0) Stop "BuildHamiltonianDense: deallocation failed"
  End Subroutine BuildHamiltonianDense

End Module BuildHamDense


Program PairFCI_Dense_XXZ_J1J2
  Use Precision
  Use Constants
  Use OperatorsPair
  Use IntegralsXXZJ1J2_1D
  Use BuildHamDense
  Implicit None

  Integer :: NLevels, NDet, IAlloc, info, lwork
  Real(Kind=pr) :: Delta, J2
  Logical :: Periodic
  Real(Kind=pr), Allocatable :: H(:,:), evals(:), work(:)
  Character(len=1) :: jobz, uplo

  External :: dsyev

  !-----------------------
  ! User parameters
  !-----------------------
  NLevels  = 6              ! <-- increase carefully (dense!)
  Delta    = -0.40_pr
  J2       = 0.0_pr          ! 0 => pure XXZ, else J1-J2 XXZ
  Periodic = .false.

  NDet = ISHFT(1, NLevels)
  Write(*,'(A,F10.4)') "Delta   = ", Delta
  Write(*,'(A,I0)') "NLevels = ", NLevels
  Write(*,'(A,I0)') "NDet    = ", NDet

  ! Build operators
  Call SetUpPairOperators(NLevels)

  ! Build integrals (spin-form). Constant 1/4 in SzSz is included automatically via Sz=(N-1)/2.
  Call SetUpIntegrals(NLevels)
  If (Abs(J2) < 1.0e-14_pr) Then
    Call DoIntegralsXXZ_1D(NLevels, Delta, Periodic)
  Else
    Call DoIntegralsJ1J2XXZ_1D(NLevels, Delta, J2, Periodic)
  End If

  ! Allocate and build Hamiltonian
  Allocate(H(NDet,NDet), evals(NDet), Stat=IAlloc)
  If (IAlloc /= 0) Stop "Allocation failed for H/evals"

  Call BuildHamiltonianDense(NLevels, H)

  ! Diagonalize (real symmetric)
  jobz = 'V'
  uplo = 'U'
  lwork = Max(1, 3*NDet - 1)
  Allocate(work(lwork), Stat=IAlloc)
  If (IAlloc /= 0) Stop "Allocation failed for work"

  Call dsyev(jobz, uplo, NDet, H, NDet, evals, work, lwork, info)
  If (info /= 0) Then
    Write(*,*) "DSYEV failed, info=", info
    Stop
  End If

  Write(*,'(A,F20.14)') "Ground-state energy = ", evals(1)

  ! Cleanup
  DeAllocate(work, H, evals, Stat=IAlloc)
  Call ShutDownIntegrals()
  Call ShutDownPairOperators()

End Program PairFCI_Dense_XXZ_J1J2


!--------------------------
! LAPACK interface
!--------------------------
 
