module Chebyshev

  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic

  implicit none
  private
  public ChebMat, ChebCalc

  real(kind = c_double), parameter :: ZERO = 0._c_double
  real(kind = c_double), parameter :: ONE = 1._c_double

contains

  pure elemental function chebPoly(x, j) result(cP)

  real(kind = c_double), intent(in)      :: x
  integer(kind = c_int), intent(in)      :: j
  real(kind = c_double)                  :: cP

    if (x < -ONE) then
      cP = -ONE ** j * cosh(j * acosh(-x))
    else if (x <= ONE) then
      cP = cos(j * acos(x))
    else
      cP = cosh(j * acosh(x))
    end if

  end function ChebPoly

  pure subroutine ChebMat(x, m, n, mat) bind(C, name = "chebM_f_")

  integer(kind = c_int), intent(in), value                 :: m, n
  real(kind = c_double), intent(in), dimension(m)          :: x
  real(kind = c_double), intent(out), dimension(m, n)      :: mat
  integer(kind = c_int)                                    :: i

    do i = 0, n
      mat(:, i) = chebPoly(x, i)
    end do

  end subroutine ChebMat

  subroutine ChebCalc(x, m, a, n, ret) bind(C, name = "chebC_f_")

  external dgemv

  integer(kind = c_int), intent(in), value                 :: m, n
  real(kind = c_double), intent(in), dimension(m)          :: x
  real(kind = c_double), intent(in), dimension(n)          :: a
  real(kind = c_double), intent(out), dimension(m)         :: ret
  real(kind = c_double), dimension(m, n)                   :: mat
  integer(kind = c_int)                                    :: i

    do i = 0, n
      mat(:, i) = chebPoly(x, i)
    end do

    call dgemv('N', m, n, ONE, mat, m, a, 1, ZERO, ret, 1)

  end subroutine ChebCalc

end module Chebyshev
