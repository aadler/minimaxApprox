module Chebyshev

  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic

  implicit none
  private

  real(kind = c_double), parameter :: ZERO = 0._c_double
  real(kind = c_double), parameter :: HALF = 0.5_c_double
  real(kind = c_double), parameter :: ONE = 1._c_double

contains

  pure elemental function chebPoly(x, j) result(cP)

  real(kind = c_double), intent(in)      :: x
  integer(kind = c_int), intent(in)      :: j
  real(kind = c_double)                  :: cP, jj

    jj = REAL(j, c_double)

    if (x < -ONE) then
      cP = -ONE ** jj * cosh(jj * acosh(-x))
    else if (x <= ONE) then
      cP = cos(jj * acos(x))
    else
      cP = cosh(jj * acosh(x))
    end if

  end function ChebPoly

end module Chebyshev
