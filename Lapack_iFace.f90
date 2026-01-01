module LapackIface
  use Precision
  implicit none
  interface
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import pr
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      real(kind=pr), intent(inout) :: a(lda,*)
      real(kind=pr), intent(out) :: w(*)
      real(kind=pr), intent(inout) :: work(*)
      integer, intent(out) :: info
    end subroutine dsyev
  end interface
end module LapackIface
