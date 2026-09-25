






      SUBROUTINE MATINV ( nSys, M, Minv, indx_B, vv_B)

! ----------------------------------------------------------------------

! Uses "ludcmp, lubksb" from Numerical Recipes to invert "M ( nSys, nSys)"
! and return it as "Minv ( nSys, nSys)". The matrix "M" is destroyed in the process.
! "M" must be square matrix.

! ----------------------------------------------------------------------
!$acc routine seq

      IMPLICIT NONE

      INTEGER nSys, ifail, i, j, indx_B ( nSys)

      REAL dd, small

      REAL, DIMENSION ( nSys, nSys) :: M, Minv

      REAL, DIMENSION ( nSys) :: vv_B

! ----------------------------------------------------------------------
!$acc routine( ludcmp ) seq
!$acc routine( lubksb ) seq
! ----------------------------------------------------------------------

      small = 1.0e-30

      ifail = 0

      CALL ludcmp ( M, nSys, nSys, indx_B, dd, small, ifail, vv_B)

      Minv = 0.0

      DO j = 1, nSys
      Minv ( j, j) = 1.0
      END DO

      DO j = 1, nSys
        CALL lubksb ( M, nSys, nSys, indx_B, Minv ( 1, j) )
      END DO

! ----------------------------------------------------------------------

      END SUBROUTINE MATINV

