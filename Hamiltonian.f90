Module Hamiltonian
  Use Precision
  Use Constants
  Implicit None

  Real(Kind=pr) :: H000
  Real(Kind=pr), Allocatable :: H010(:)       ! N_p
  Real(Kind=pr), Allocatable :: H101(:,:)     ! Pdag_p P_q
  Real(Kind=pr), Allocatable :: H020(:,:)     ! N_p N_q

Contains

  Subroutine AllocIntegrals(N)
    Implicit None
    Integer, Intent(In) :: N
    Integer :: IAlloc
    Allocate(H010(N), H101(N,N), H020(N,N), Stat=IAlloc)
    If (IAlloc /= 0) Stop "AllocIntegrals: allocation failed"
  End Subroutine AllocIntegrals

  Subroutine ZeroIntegrals()
    Implicit None
    H000 = Zero
    H010 = Zero
    H101 = Zero
    H020 = Zero
  End Subroutine ZeroIntegrals

  Subroutine FreeIntegrals()
    Implicit None
    Integer :: IAlloc
    If (Allocated(H010)) DeAllocate(H010, Stat=IAlloc)
    If (Allocated(H101)) DeAllocate(H101, Stat=IAlloc)
    If (Allocated(H020)) DeAllocate(H020, Stat=IAlloc)
  End Subroutine FreeIntegrals

  !===========================================================
  ! Add a single bond contribution:
  !   Hij = 1/2 (Pdag_i P_j + Pdag_j P_i) + Delta * Sz_i Sz_j
  ! with Sz=(N-1)/2 => Delta/4*(N_i N_j - N_i - N_j + 1)
  !
  ! fac counts multiplicity (usually 1).
  !===========================================================
  Subroutine AddBond_XXZ(i, j, fac, delta)
    Implicit None
    Integer, Intent(In) :: i, j, fac
    Real(Kind=pr), Intent(In) :: delta
    Real(Kind=pr) :: f

    f = Real(fac,Kind=pr)

    ! flip-flop
    H101(i,j) = H101(i,j) + f*F12
    H101(j,i) = H101(j,i) + f*F12

    ! SzSz expanded:
    ! + delta/4 * N_i N_j
    ! - delta/4 * N_i
    ! - delta/4 * N_j
    ! + delta/4 * 1  (THIS is the +1/4 you want included)
    H020(i,j) = H020(i,j) + f*delta*F14
    H020(j,i) = H020(j,i) + f*delta*F14
    H010(i)   = H010(i)   - f*delta*F14
    H010(j)   = H010(j)   - f*delta*F14
    H000      = H000      + f*delta*F14
  End Subroutine AddBond_XXZ

  !===========================================================
  ! Add a single Heisenberg (XXX) bond:
  !   Hij = 1/2 (Pdag_i P_j + Pdag_j P_i) + Sz_i Sz_j
  ! which is AddBond_XXZ with delta=1.
  !===========================================================
  Subroutine AddBond_XXX(i, j, fac, Jc)
    Implicit None
    Integer, Intent(In) :: i, j, fac
    Real(Kind=pr), Intent(In) :: Jc
    Real(Kind=pr) :: f

    f = Real(fac,Kind=pr)

    H101(i,j) = H101(i,j) + f*Jc*F12
    H101(j,i) = H101(j,i) + f*Jc*F12

    H020(i,j) = H020(i,j) + f*Jc*F14
    H020(j,i) = H020(j,i) + f*Jc*F14
    H010(i)   = H010(i)   - f*Jc*F14
    H010(j)   = H010(j)   - f*Jc*F14
    H000      = H000      + f*Jc*F14
  End Subroutine AddBond_XXX

  !===========================================================
  ! 1D XXZ on a chain (PBC/OBC)
  !===========================================================
  Subroutine Build_XXZ_1D(L, delta, periodic)
    Implicit None
    Integer, Intent(In) :: L
    Real(Kind=pr), Intent(In) :: delta
    Logical, Intent(In) :: periodic
    Integer :: i, j

    Call ZeroIntegrals()

    Do i = 1, L-1
      j = i + 1
      Call AddBond_XXZ(i, j, 1, delta)
    End Do
    If (periodic .And. L > 2) Then
      Call AddBond_XXZ(1, L, 1, delta)
    End If
  End Subroutine Build_XXZ_1D

  !===========================================================
  ! 1D J1-J2 XXZ:
  ! bonds at distance 1 with weight 1
  ! bonds at distance 2 with weight j2
  !===========================================================
  Subroutine Build_J1J2XXZ_1D(L, delta, j2, periodic)
    Implicit None
    Integer, Intent(In) :: L
    Real(Kind=pr), Intent(In) :: delta, j2
    Logical, Intent(In) :: periodic
    Integer :: i, j

    Call ZeroIntegrals()

    ! J1 bonds
    Do i = 1, L-1
      j = i + 1
      Call AddBond_XXZ(i, j, 1, delta)
    End Do
    If (periodic .And. L > 2) Then
      Call AddBond_XXZ(1, L, 1, delta)
    End If

    ! J2 bonds (scaled by j2)
    Do i = 1, L-2
      j = i + 2
      Call AddBond_XXZ_scaled(i, j, 1, delta, j2)
    End Do
    If (periodic .And. L > 3) Then
      Call AddBond_XXZ_scaled(1, L-1, 1, delta, j2)
      Call AddBond_XXZ_scaled(2, L,   1, delta, j2)
    End If

  Contains
    Subroutine AddBond_XXZ_scaled(i, j, fac, delta, scale)
      Implicit None
      Integer, Intent(In) :: i, j, fac
      Real(Kind=pr), Intent(In) :: delta, scale
      Integer :: k
      ! just reuse AddBond_XXZ by multiplying "fac" into scale safely
      ! simplest: apply manually:
      Real(Kind=pr) :: f
      f = Real(fac,Kind=pr) * scale

      H101(i,j) = H101(i,j) + f*F12
      H101(j,i) = H101(j,i) + f*F12

      H020(i,j) = H020(i,j) + f*delta*F14
      H020(j,i) = H020(j,i) + f*delta*F14
      H010(i)   = H010(i)   - f*delta*F14
      H010(j)   = H010(j)   - f*delta*F14
      H000      = H000      + f*delta*F14
    End Subroutine AddBond_XXZ_scaled
  End Subroutine Build_J1J2XXZ_1D


  !===========================================================
  ! 2D square XXZ with optional PBC
  ! Build bonds once: only +x and +y directions to avoid double counting.
  !===========================================================
  Subroutine Build_XXZ_Square(nx, ny, delta, periodic)
    Implicit None
    Integer, Intent(In) :: nx, ny
    Real(Kind=pr), Intent(In) :: delta
    Logical, Intent(In) :: periodic
    Integer :: ix, iy, i, jx, jy, j, nmo

    nmo = nx*ny
    Call ZeroIntegrals()

    Do ix = 0, nx-1
      Do iy = 0, ny-1
        i = ix*ny + iy + 1

        ! +x bond
        jx = ix + 1
        jy = iy
        If (periodic) Then
          jx = Mod(jx, nx)
        Else
          If (jx >= nx) GoTo 10
        End If
        j = jx*ny + jy + 1
        Call AddBond_XXZ(i, j, 1, delta)
10      Continue

        ! +y bond
        jx = ix
        jy = iy + 1
        If (periodic) Then
          jy = Mod(jy, ny)
        Else
          If (jy >= ny) Cycle
        End If
        j = jx*ny + jy + 1
        Call AddBond_XXZ(i, j, 1, delta)

      End Do
    End Do
  End Subroutine Build_XXZ_Square


  !===========================================================
  ! 2D square J1-J2 XXX (no delta): J1 on NN, J2 on diagonal NNN
  ! If you want XXZ anisotropy here too, tell me and I'll add delta.
  !===========================================================
  Subroutine Build_J1J2_Square(nx, ny, j2, periodic)
    Implicit None
    Integer, Intent(In) :: nx, ny
    Real(Kind=pr), Intent(In) :: j2
    Logical, Intent(In) :: periodic
    Integer :: ix, iy, i, j, nmo
    Integer :: jx, jy

    nmo = nx*ny
    Call ZeroIntegrals()

    Do ix = 0, nx-1
      Do iy = 0, ny-1
        i = ix*ny + iy + 1

        ! NN: +x and +y, weight J1=1
        ! +x
        jx = ix + 1; jy = iy
        If (periodic) Then
          jx = Mod(jx, nx)
        Else
          If (jx >= nx) GoTo 20
        End If
        j = jx*ny + jy + 1
        Call AddBond_XXX(i, j, 1, One)
20      Continue
        ! +y
        jx = ix; jy = iy + 1
        If (periodic) Then
          jy = Mod(jy, ny)
        Else
          If (jy >= ny) GoTo 30
        End If
        j = jx*ny + jy + 1
        Call AddBond_XXX(i, j, 1, One)
30      Continue

        ! Diagonal NNN: (+1,+1) and (+1,-1) to count each diagonal bond once
        ! (+1,+1)
        jx = ix + 1; jy = iy + 1
        If (periodic) Then
          jx = Mod(jx, nx); jy = Mod(jy, ny)
          j = jx*ny + jy + 1
          Call AddBond_XXX(i, j, 1, j2)
        Else
          If (jx < nx .And. jy < ny) Then
            j = jx*ny + jy + 1
            Call AddBond_XXX(i, j, 1, j2)
          End If
        End If

        ! (+1,-1)
        jx = ix + 1; jy = iy - 1
        If (periodic) Then
          jx = Mod(jx, nx)
          jy = Mod(jy, ny); If (jy < 0) jy = jy + ny
          j = jx*ny + jy + 1
          Call AddBond_XXX(i, j, 1, j2)
        Else
          If (jx < nx .And. jy >= 0) Then
            j = jx*ny + jy + 1
            Call AddBond_XXX(i, j, 1, j2)
          End If
        End If

      End Do
    End Do
  End Subroutine Build_J1J2_Square

End Module Hamiltonian
