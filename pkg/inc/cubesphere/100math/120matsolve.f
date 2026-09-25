





      SUBROUTINE MATSOLVE ( nSys, M, B, indx_B, vv_B)

! ----------------------------------------------------------------------

! Uses "ludcmp, lubksb" from Numerical Recipes to solve, "[M][X] = [B]"
! "M ( nSys, nSys), B ( nSys)". The output [X] is returned in "B" itself.
! The matrix "M" is destroyed in the process. "M" must be square matrix.

! ----------------------------------------------------------------------
!$acc routine seq

      IMPLICIT NONE

      INTEGER nSys, ifail, i, j, indx_B ( nSys)

      REAL dd, small

      REAL, DIMENSION ( nSys, nSys) :: M

      REAL, DIMENSION ( nSys) :: B, vv_B

! ----------------------------------------------------------------------
!$acc routine( ludcmp ) seq
!$acc routine( lubksb ) seq
! ----------------------------------------------------------------------

      small = 1.0e-30

      ifail = 0

      CALL ludcmp ( M, nSys, nSys, indx_B, dd, small, ifail, vv_B)

      CALL lubksb ( M, nSys, nSys, indx_B, B )

! ----------------------------------------------------------------------

      END SUBROUTINE MATSOLVE
