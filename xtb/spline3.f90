! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

subroutine spline_cubic_set ( n, t, y, ibcbeg, ybcbeg, ibcend, ybcend, ypp )

!*****************************************************************************80
!
!! SPLINE_CUBIC_SET computes the second derivatives of a piecewise cubic spline.
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to
!    determine the second derivative data, passing in the data to be
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output,
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) )
!             + C(IVAL) * ( T - T(IVAL) )^2
!             + D(IVAL) * ( T - T(IVAL) )^3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )^2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1))
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))^2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      =
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL)
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!
!    Boundary conditions must be applied at the first and last knots.
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data points; N must be
!    at least 2.
!
!    Input, real ( kind = 8 ) T(N), the points where data is specified.
!    The values should be distinct, and increasing.
!
!    Input, real ( kind = 8 ) Y(N), the data values to be interpolated.
!
!    Input, integer ( kind = 4 ) IBCBEG, the left boundary condition flag:
!    0: the spline should be a quadratic over the first interval;
!    1: the first derivative at the left endpoint should be YBCBEG;
!    2: the second derivative at the left endpoint should be YBCBEG;
!    3: Not-a-knot: the third derivative is continuous at T(2).
!
!    Input, real ( kind = 8 ) YBCBEG, the left boundary value, if needed.
!
!    Input, integer ( kind = 4 ) IBCEND, the right boundary condition flag:
!    0: the spline should be a quadratic over the last interval;
!    1: the first derivative at the right endpoint should be YBCEND;
!    2: the second derivative at the right endpoint should be YBCEND;
!    3: Not-a-knot: the third derivative is continuous at T(N-1).
!
!    Input, real ( kind = 8 ) YBCEND, the right boundary value, if needed.
!
!    Output, real ( kind = 8 ) YPP(N), the second derivatives of
!    the cubic spline.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)
  real ( kind = 8 ) a4(n)
  real ( kind = 8 ) a5(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibcbeg
  integer ( kind = 4 ) ibcend
  integer ( kind = 4 ) info
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ybcbeg
  real ( kind = 8 ) ybcend
  real ( kind = 8 ) ypp(n)
!
!  Check.
!
  if ( n <= 1 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The number of knots must be at least 2.'
    write ( *, '(a,i8)' ) '  The input value of N = ', n
    stop 1
  end if

  do i = 1, n - 1
    if ( t(i+1) <= t(i) ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
      write ( *, '(a)' ) '  The knots must be strictly increasing, but'
      write ( *, '(a,i8,a,g14.6)' ) '  T(',  i,') = ', t(i)
      write ( *, '(a,i8,a,g14.6)' ) '  T(',i+1,') = ', t(i+1)
      stop 1
    end if
  end do
!
!  Zero out the matrix.
!
  a1(1:n) = 0.0D+00
  a2(1:n) = 0.0D+00
  a3(1:n) = 0.0D+00
  a4(1:n) = 0.0D+00
  a5(1:n) = 0.0D+00
!
!  Set the first equation.
!
  if ( ibcbeg == 0 ) then
    b(1) = 0.0D+00
    a3(1) =  1.0D+00
    a4(1) = -1.0D+00
  else if ( ibcbeg == 1 ) then
    b(1) = ( y(2) - y(1) ) / ( t(2) - t(1) ) - ybcbeg
    a3(1) = ( t(2) - t(1) ) / 3.0D+00
    a4(1) = ( t(2) - t(1) ) / 6.0D+00
  else if ( ibcbeg == 2 ) then
    b(1) = ybcbeg
    a3(1) = 1.0D+00
    a4(1) = 0.0D+00
  else if ( ibcbeg == 3 ) then
    b(1) = 0.0D+00
    a3(1) = - ( t(3) - t(2) )
    a4(1) =   ( t(3)        - t(1) )
    a5(1) = - (        t(2) - t(1) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCBEG must be 0, 1, 2 or 3.'
    write ( *, '(a,i8)' ) '  The input value is IBCBEG = ', ibcbeg
    stop 1
  end if
!
!  Set the intermediate equations.
!
  do i = 2, n - 1
    b(i) = ( y(i+1) - y(i) ) / ( t(i+1) - t(i) ) &
         - ( y(i) - y(i-1) ) / ( t(i) - t(i-1) )
    a2(i) = ( t(i+1) - t(i)   ) / 6.0D+00
    a3(i) = ( t(i+1) - t(i-1) ) / 3.0D+00
    a4(i) = ( t(i)   - t(i-1) ) / 6.0D+00
  end do
!
!  Set the last equation.
!
  if ( ibcend == 0 ) then
    b(n) = 0.0D+00
    a2(n) = -1.0D+00
    a3(n) = 1.0D+00
  else if ( ibcend == 1 ) then
    b(n) = ybcend - ( y(n) - y(n-1) ) / ( t(n) - t(n-1) )
    a2(n) = ( t(n) - t(n-1) ) / 6.0D+00
    a3(n) = ( t(n) - t(n-1) ) / 3.0D+00
  else if ( ibcend == 2 ) then
    b(n) = ybcend
    a2(n) = 0.0D+00
    a3(n) = 1.0D+00
  else if ( ibcend == 3 ) then
    b(n) = 0.0D+00
    a1(n) = - ( t(n) - t(n-1) )
    a2(n) =   ( t(n)          - t(n-2) )
    a3(n) = - (        t(n-1) - t(n-2) )
  else
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, '(a)' ) '  The boundary flag IBCEND must be 0, 1, 2 or 3.'
    write ( *, '(a,i8)' ) '  The input value is IBCEND = ', ibcend
    stop 1
  end if
!
!  Special case:
!    N = 2, IBCBEG = IBCEND = 0.
!
  if ( n == 2 .and. ibcbeg == 0 .and. ibcend == 0 ) then

    ypp(1) = 0.0D+00
    ypp(2) = 0.0D+00
!
!  Solve the linear system.
!
  else

    call penta ( n, a1, a2, a3, a4, a5, b, ypp )

  end if

  return
end

subroutine penta ( n, a1, a2, a3, a4, a5, b, x )

!*****************************************************************************80
!
!! PENTA solves a pentadiagonal system of linear equations.
!
!  Discussion:
!
!    The matrix A is pentadiagonal.  It is entirely zero, except for
!    the main diagaonal, and the two immediate sub- and super-diagonals.
!
!    The entries of Row I are stored as:
!
!      A(I,I-2) -> A1(I)
!      A(I,I-1) -> A2(I)
!      A(I,I)   -> A3(I)
!      A(I,I+1) -> A4(I)
!      A(I,I-2) -> A5(I)
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 June 2013
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Cheney, Kincaid,
!    Numerical Mathematics and Computing,
!    1985, pages 233-236.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the order of the matrix.
!
!    Input, real ( kind = 8 ) A1(N), A2(N), A3(N), A4(N), A5(N), the nonzero
!    elements of the matrix.  Note that the data in A2, A3 and A4
!    is overwritten by this routine during the solution process.
!
!    Input, real ( kind = 8 ) B(N), the right hand side of the linear system.
!
!    Output, real ( kind = 8 ) X(N), the solution of the linear system.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) a1(n)
  real ( kind = 8 ) a2(n)
  real ( kind = 8 ) a3(n)
  real ( kind = 8 ) a4(n)
  real ( kind = 8 ) a5(n)
  real ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xmult

  do i = 2, n - 1
    xmult = a2(i) / a3(i-1)
    a3(i) = a3(i) - xmult * a4(i-1)
    a4(i) = a4(i) - xmult * a5(i-1)
    b(i) = b(i) - xmult * b(i-1)
    xmult = a1(i+1) / a3(i-1)
    a2(i+1) = a2(i+1) - xmult * a4(i-1)
    a3(i+1) = a3(i+1) - xmult * a5(i-1)
    b(i+1) = b(i+1) - xmult * b(i-1)
  end do

  xmult = a2(n) / a3(n-1)
  a3(n) = a3(n) - xmult * a4(n-1)
  x(n) = ( b(n) - xmult * b(n-1) ) / a3(n)
  x(n-1) = ( b(n-1) - a4(n-1) * x(n) ) / a3(n-1)
  do i = n - 2, 1, -1
    x(i) = ( b(i) - a4(i) * x(i+1) - a5(i) * x(i+2) ) / a3(i)
  end do

  return
end

subroutine spline_cubic_val ( n, t, y, ypp, tval, yval, ypval, yppval )

!*****************************************************************************80
!
!! SPLINE_CUBIC_VAL evaluates a piecewise cubic spline at a point.
!
!  Discussion:
!
!    SPLINE_CUBIC_SET must have already been called to define the
!    values of YPP.
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A
!             + B * ( T - T(IVAL) )
!             + C * ( T - T(IVAL) )^2
!             + D * ( T - T(IVAL) )^3
!
!    Here:
!      A = Y(IVAL)
!      B = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C = YPP(IVAL) / 2
!      D = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    20 November 2000
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Carl deBoor,
!    A Practical Guide to Splines,
!    Springer, 2001,
!    ISBN: 0387953663.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number of data values.
!
!    Input, real ( kind = 8 ) T(N), the knot values.
!
!    Input, real ( kind = 8 ) Y(N), the data values at the knots.
!
!    Input, real ( kind = 8 ) YPP(N), the second derivatives of the
!    spline at the knots.
!
!    Input, real ( kind = 8 ) TVAL, a point, typically between T(1) and
!    T(N), at which the spline is to be evalulated.  If TVAL lies outside
!    this range, extrapolation is used.
!
!    Output, real ( kind = 8 ) YVAL, YPVAL, YPPVAL, the value of the spline, and
!    its first two derivatives at TVAL.
!
  implicit none

  integer ( kind = 4 ) n

  real ( kind = 8 ) dt
  real ( kind = 8 ) h
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) t(n)
  real ( kind = 8 ) tval
  real ( kind = 8 ) y(n)
  real ( kind = 8 ) ypp(n)
  real ( kind = 8 ) yppval
  real ( kind = 8 ) ypval
  real ( kind = 8 ) yval
!
!  Determine the interval [T(LEFT), T(RIGHT)] that contains TVAL.
!  Values below T(1) or above T(N) use extrapolation.
!
  call r8vec_bracket ( n, t, tval, left, right )
!
!  Evaluate the polynomial.
!
  dt = tval - t(left)
  h = t(right) - t(left)

  yval = y(left) &
       + dt * ( ( y(right) - y(left) ) / h &
              - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
       + dt * ( 0.5D+00 * ypp(left) &
       + dt * ( ( ypp(right) - ypp(left) ) / ( 6.0D+00 * h ) ) ) )

  ypval = ( y(right) - y(left) ) / h &
       - ( ypp(right) / 6.0D+00 + ypp(left) / 3.0D+00 ) * h &
       + dt * ( ypp(left) &
       + dt * ( 0.5D+00 * ( ypp(right) - ypp(left) ) / h ) )

  yppval = ypp(left) + dt * ( ypp(right) - ypp(left) ) / h

  return
end

subroutine r8vec_bracket ( n, x, xval, left, right )

!*****************************************************************************80
!
!! R8VEC_BRACKET searches a sorted R8VEC for successive brackets of a value.
!
!  Discussion:
!
!    An R8VEC is an array of double precision real values.
!
!    If the values in the vector are thought of as defining intervals
!    on the real line, then this routine searches for the interval
!    nearest to or containing the given value.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 April 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, length of input array.
!
!    Input, real ( kind = 8 ) X(N), an array sorted into ascending order.
!
!    Input, real ( kind = 8 ) XVAL, a value to be bracketed.
!
!    Output, integer ( kind = 4 ) LEFT, RIGHT, the results of the search.
!    Either:
!      XVAL < X(1), when LEFT = 1, RIGHT = 2;
!      X(N) < XVAL, when LEFT = N-1, RIGHT = N;
!    or
!      X(LEFT) <= XVAL <= X(RIGHT).
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  integer ( kind = 4 ) left
  integer ( kind = 4 ) right
  real ( kind = 8 ) x(n)
  real ( kind = 8 ) xval

  do i = 2, n - 1

    if ( xval < x(i) ) then
      left = i - 1
      right = i
      return
    end if

   end do

  left = n - 1
  right = n

  return
end
