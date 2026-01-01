module Operators
  use Precision
  use Constants
  implicit none

  real(kind=pr), allocatable :: Pdag(:,:,:), Pmat(:,:,:), Nmat(:,:,:)
  real(kind=pr), allocatable :: Sz(:,:,:), Splus(:,:,:), Sminus(:,:,:)

contains

  pure integer function GetBit(bits, p) result(b)
    integer, intent(in) :: bits, p
    b = iand(ishft(bits, -p), 1)
  end function GetBit

  pure integer function SetBit(bits, p) result(newbits)
    integer, intent(in) :: bits, p
    newbits = ior(bits, ishft(1, p))
  end function SetBit

  pure integer function ClearBit(bits, p) result(newbits)
    integer, intent(in) :: bits, p
    newbits = iand(bits, not(ishft(1, p)))
  end function ClearBit

  pure integer function ApplyPdagBits(bits_in, p) result(bits_out)
    integer, intent(in) :: bits_in, p
    if (GetBit(bits_in, p) == 1) then
      bits_out = -1
    else
      bits_out = SetBit(bits_in, p)
    end if
  end function ApplyPdagBits

  pure integer function ApplyPBits(bits_in, p) result(bits_out)
    integer, intent(in) :: bits_in, p
    if (GetBit(bits_in, p) == 0) then
      bits_out = -1
    else
      bits_out = ClearBit(bits_in, p)
    end if
  end function ApplyPBits

  pure integer function ApplyNval(bits_in, p) result(nval)
    integer, intent(in) :: bits_in, p
    if (GetBit(bits_in, p) == 1) then
      nval = 2
    else
      nval = 0
    end if
  end function ApplyNval

  subroutine SetUpPairOperators(NLevels)
    integer, intent(in) :: NLevels
    integer :: NDet, p, mu, mu_bits, nu_bits, IAlloc
    real(kind=pr), allocatable :: Iden(:,:)

    NDet = ishft(1, NLevels)

    allocate(Pdag(NDet,NDet,NLevels), Pmat(NDet,NDet,NLevels), Nmat(NDet,NDet,NLevels), &
             Sz(NDet,NDet,NLevels), Splus(NDet,NDet,NLevels), Sminus(NDet,NDet,NLevels), &
             Iden(NDet,NDet), stat=IAlloc)
    if (IAlloc /= 0) stop "SetUpPairOperators: allocation failed"

    Pdag = Zero; Pmat = Zero; Nmat = Zero
    Sz = Zero; Splus = Zero; Sminus = Zero

    Iden = Zero
    do mu = 1, NDet
      Iden(mu,mu) = One
    end do

    do p = 1, NLevels
      do mu = 1, NDet
        mu_bits = mu - 1

        nu_bits = ApplyPdagBits(mu_bits, p-1)
        if (nu_bits >= 0) Pdag(nu_bits+1, mu, p) = One

        nu_bits = ApplyPBits(mu_bits, p-1)
        if (nu_bits >= 0) Pmat(nu_bits+1, mu, p) = One

        Nmat(mu, mu, p) = real(ApplyNval(mu_bits, p-1), kind=pr)  ! 0 or 2
      end do

      Splus(:,:,p)  = Pdag(:,:,p)
      Sminus(:,:,p) = Pmat(:,:,p)
      Sz(:,:,p)     = Half*Nmat(:,:,p) - Half*Iden
    end do

    deallocate(Iden, stat=IAlloc)
    if (IAlloc /= 0) stop "SetUpPairOperators: deallocation failed"
  end subroutine SetUpPairOperators

  subroutine ShutDownPairOperators()
    integer :: IAlloc
    if (allocated(Pdag))   deallocate(Pdag, stat=IAlloc)
    if (allocated(Pmat))   deallocate(Pmat, stat=IAlloc)
    if (allocated(Nmat))   deallocate(Nmat, stat=IAlloc)
    if (allocated(Sz))     deallocate(Sz, stat=IAlloc)
    if (allocated(Splus))  deallocate(Splus, stat=IAlloc)
    if (allocated(Sminus)) deallocate(Sminus, stat=IAlloc)
  end subroutine ShutDownPairOperators

end module Operators
