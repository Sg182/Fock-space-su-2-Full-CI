module Integrals
  use Precision
  use Constants
  implicit none

  real(kind=pr) :: H000
  real(kind=pr), allocatable :: H010(:)
  real(kind=pr), allocatable :: H020(:,:)
  real(kind=pr), allocatable :: H101(:,:)

contains

  subroutine SetUpIntegrals(NLevels)
    integer, intent(in) :: NLevels
    integer :: IAlloc
    allocate(H010(NLevels), H020(NLevels,NLevels), H101(NLevels,NLevels), stat=IAlloc)
    if (IAlloc /= 0) stop "SetUpIntegrals: allocation failed"
    H000 = Zero; H010 = Zero; H020 = Zero; H101 = Zero
  end subroutine SetUpIntegrals

  subroutine ShutDownIntegrals()
    integer :: IAlloc
    if (allocated(H010)) deallocate(H010, stat=IAlloc)
    if (allocated(H020)) deallocate(H020, stat=IAlloc)
    if (allocated(H101)) deallocate(H101, stat=IAlloc)
  end subroutine ShutDownIntegrals

  subroutine DoIntegralsXXZ_1D(NLevels, Delta, Periodic)
    integer, intent(in) :: NLevels
    real(kind=pr), intent(in) :: Delta
    logical, intent(in) :: Periodic
    integer :: i, j

    H000 = Zero; H010 = Zero; H020 = Zero; H101 = Zero

    do i = 1, NLevels-1
      j = i + 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    end do

    if (Periodic) then
      i = NLevels; j = 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    end if
  end subroutine DoIntegralsXXZ_1D

  subroutine DoIntegralsJ1J2XXZ_1D(NLevels, Delta, J2, Periodic)
    integer, intent(in) :: NLevels
    real(kind=pr), intent(in) :: Delta, J2
    logical, intent(in) :: Periodic
    integer :: i, j

    H000 = Zero; H010 = Zero; H020 = Zero; H101 = Zero

    ! 1NN
    do i = 1, NLevels-1
      j = i + 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    end do
    if (Periodic) then
      i = NLevels; j = 1
      H101(i,j) = H101(i,j) + Half
      H101(j,i) = H101(j,i) + Half
      H020(i,j) = H020(i,j) + Delta
      H020(j,i) = H020(j,i) + Delta
    end if

    ! 2NN
    do i = 1, NLevels-2
      j = i + 2
      H101(i,j) = H101(i,j) + Half*J2
      H101(j,i) = H101(j,i) + Half*J2
      H020(i,j) = H020(i,j) + Delta*J2
      H020(j,i) = H020(j,i) + Delta*J2
    end do
    if (Periodic) then
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
    end if
  end subroutine DoIntegralsJ1J2XXZ_1D


  !===========================================================
  ! 2D square XXZ (pair language)
  ! Sites are indexed by: i = ix*ny + iy + 1   with ix=0..nx-1, iy=0..ny-1
  !
  ! Adds NN bonds once: (+x) and (+y) only.
  ! Convention matches your 1D: H101 += 1/2, H020 += Delta
  ! No constant (+1/4) and no one-body shifts here.
  !===========================================================
  subroutine DoIntegralsXXZ_2D(nx, ny, Delta, PBCx, PBCy)
    integer, intent(in) :: nx, ny
    real(kind=pr), intent(in) :: Delta
    logical, intent(in) :: PBCx, PBCy
    integer :: ix, iy
    integer :: i, j
    integer :: jx, jy
    integer :: nSites

    nSites = nx*ny
    H000 = Zero; H010 = Zero; H020 = Zero; H101 = Zero

    do ix = 0, nx-1
      do iy = 0, ny-1
        i = ix*ny + iy + 1

        ! +x neighbor
        jx = ix + 1
        jy = iy
        if (jx >= nx) then
          if (PBCx) then
            jx = 0
          else
            goto 10
          end if
        end if
        j = jx*ny + jy + 1
        H101(i,j) = H101(i,j) + Half
        H101(j,i) = H101(j,i) + Half
        H020(i,j) = H020(i,j) + Delta
        H020(j,i) = H020(j,i) + Delta
10      continue

        ! +y neighbor
        jx = ix
        jy = iy + 1
        if (jy >= ny) then
          if (PBCy) then
            jy = 0
          else
            cycle
          end if
        end if
        j = jx*ny + jy + 1
        H101(i,j) = H101(i,j) + Half
        H101(j,i) = H101(j,i) + Half
        H020(i,j) = H020(i,j) + Delta
        H020(j,i) = H020(j,i) + Delta

      end do
    end do

  end subroutine DoIntegralsXXZ_2D


  !===========================================================
  ! 2D square J1-J2 XXZ:
  !   J1 on NN bonds (+x,+y) with weight 1
  !   J2 on diagonal NNN bonds (+1,+1) and (+1,-1) with weight J2
  ! Anisotropy Delta enters the same way as your 1D:
  !   H101 += (1/2)*Jbond,   H020 += Delta*Jbond
  !===========================================================
  subroutine DoIntegralsJ1J2XXZ_2D(nx, ny, Delta, J2, PBCx, PBCy)
    integer, intent(in) :: nx, ny
    real(kind=pr), intent(in) :: Delta, J2
    logical, intent(in) :: PBCx, PBCy
    integer :: ix, iy
    integer :: i, j
    integer :: jx, jy

    H000 = Zero; H010 = Zero; H020 = Zero; H101 = Zero

    do ix = 0, nx-1
      do iy = 0, ny-1
        i = ix*ny + iy + 1

        !-------------------------
        ! J1: NN (+x)
        !-------------------------
        jx = ix + 1
        jy = iy
        if (jx >= nx) then
          if (PBCx) then
            jx = 0
          else
            goto 20
          end if
        end if
        j = jx*ny + jy + 1
        H101(i,j) = H101(i,j) + Half
        H101(j,i) = H101(j,i) + Half
        H020(i,j) = H020(i,j) + Delta
        H020(j,i) = H020(j,i) + Delta
20      continue

        !-------------------------
        ! J1: NN (+y)
        !-------------------------
        jx = ix
        jy = iy + 1
        if (jy >= ny) then
          if (PBCy) then
            jy = 0
          else
            goto 30
          end if
        end if
        j = jx*ny + jy + 1
        H101(i,j) = H101(i,j) + Half
        H101(j,i) = H101(j,i) + Half
        H020(i,j) = H020(i,j) + Delta
        H020(j,i) = H020(j,i) + Delta
30      continue

        !-------------------------
        ! J2: diagonal (+1,+1)
        !-------------------------
        jx = ix + 1
        jy = iy + 1
        if (jx >= nx) then
          if (PBCx) then
            jx = 0
          else
            goto 40
          end if
        end if
        if (jy >= ny) then
          if (PBCy) then
            jy = 0
          else
            goto 40
          end if
        end if
        j = jx*ny + jy + 1
        H101(i,j) = H101(i,j) + Half*J2
        H101(j,i) = H101(j,i) + Half*J2
        H020(i,j) = H020(i,j) + Delta*J2
        H020(j,i) = H020(j,i) + Delta*J2
40      continue

        !-------------------------
        ! J2: diagonal (+1,-1)
        !-------------------------
        jx = ix + 1
        jy = iy - 1
        if (jx >= nx) then
          if (PBCx) then
            jx = 0
          else
            cycle
          end if
        end if
        if (jy < 0) then
          if (PBCy) then
            jy = ny - 1
          else
            cycle
          end if
        end if
        j = jx*ny + jy + 1
        H101(i,j) = H101(i,j) + Half*J2
        H101(j,i) = H101(j,i) + Half*J2
        H020(i,j) = H020(i,j) + Delta*J2
        H020(j,i) = H020(j,i) + Delta*J2

      end do
    end do

  end subroutine DoIntegralsJ1J2XXZ_2D

end module Integrals
