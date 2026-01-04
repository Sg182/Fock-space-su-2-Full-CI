module CIUtils
  use Precision
  implicit none
  private
  public :: WriteTopCICoeffs

contains

  subroutine WriteTopCICoeffs(fname, NS, C, nprint)
    character(len=*), intent(in) :: fname
    integer, intent(in) :: NS, nprint
    real(kind=pr), intent(in) :: C(:)

    integer :: NDet, i, j, k, idx, u, nout
    real(kind=pr), allocatable :: absC(:)
    integer, allocatable :: order(:)

    NDet = size(C)
    nout = min(nprint, NDet)

    allocate(absC(NDet), order(NDet))
    do i = 1, NDet
      absC(i) = abs(C(i))
      order(i) = i
    end do

    ! Partial selection sort: only top nout
    do k = 1, nout
      idx = k
      do j = k+1, NDet
        if (absC(order(j)) > absC(order(idx))) idx = j
      end do
      call swap_int(order(k), order(idx))
    end do

    open(newunit=u, file=fname, status="replace", action="write", form="formatted")
    write(u,'(A)') "# rank   mu        bitstring            C(mu)                |C(mu)|"
    do k = 1, nout
      i = order(k)
      write(u,'(I5,2X,I10,2X,A,2X,ES22.14,2X,ES12.4)') &
           k, i-1, Bits(i-1, NS), C(i), absC(i)
    end do
    close(u)

    deallocate(absC, order)
  end subroutine WriteTopCICoeffs


  subroutine swap_int(a, b)
    integer, intent(inout) :: a, b
    integer :: t
    t = a
    a = b
    b = t
  end subroutine swap_int


  function Bits(mu, NS) result(str)
    integer, intent(in) :: mu, NS
    character(len=:), allocatable :: str
    integer :: p

    allocate(character(len=NS) :: str)

    ! left = most-significant bit
    do p = 0, NS-1
      if (btest(mu,p)) then
        str(NS-p:NS-p) = '1'
      else
        str(NS-p:NS-p) = '0'
      end if
    end do
  end function Bits

end module CIUtils
