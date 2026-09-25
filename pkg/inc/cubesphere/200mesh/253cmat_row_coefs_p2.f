





      SUBROUTINE CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                               CMat_Row_B)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Third order case. Evaluates the 30 coefficients that are needed for enforcing
! Almost Divergence Preserving (ADP) reconstruction of the magnetic fields.

! "x, y, z" position of the facial quadrature point relative to the zone's
! centroid, with scaling coefficient applied.

! "nx, ny, nz" normal vector components at that facial quadrature point.

! "Coef_0 ( 10)" The 10 third order coefficients at that zone's center.

! "CMat_Row_B ( 30)" the coefficients in front of the 30 unknowns.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL x, y, z, nx, ny, nz, Coef_0 ( 10), CMat_Row_B ( 30)

! ----------------------------------------------------------------------

      CMat_Row_B ( 1) = nx

      CMat_Row_B ( 2) = nx * x

      CMat_Row_B ( 3) = nx * y

      CMat_Row_B ( 4) = nx * z

      CMat_Row_B ( 5) = nx * ( x**2 - Coef_0 ( 5))

      CMat_Row_B ( 6) = nx * ( y**2 - Coef_0 ( 6))

      CMat_Row_B ( 7) = nx * ( z**2 - Coef_0 ( 7))

      CMat_Row_B ( 8) = nx * ( x * y - Coef_0 ( 8))

      CMat_Row_B ( 9) = nx * ( y * z - Coef_0 ( 9))

      CMat_Row_B ( 10) = nx * ( x * z - Coef_0 ( 10))

      CMat_Row_B ( 11) = ny

      CMat_Row_B ( 12) = ny * x

      CMat_Row_B ( 13) = ny * y

      CMat_Row_B ( 14) = ny * z

      CMat_Row_B ( 15) = ny * ( x**2 - Coef_0 ( 5))

      CMat_Row_B ( 16) = ny * ( y**2 - Coef_0 ( 6))

      CMat_Row_B ( 17) = ny * ( z**2 - Coef_0 ( 7))

      CMat_Row_B ( 18) = ny * ( x * y - Coef_0 ( 8))

      CMat_Row_B ( 19) = ny * ( y * z - Coef_0 ( 9))

      CMat_Row_B ( 20) = ny * ( x * z - Coef_0 ( 10))

      CMat_Row_B ( 21) = nz

      CMat_Row_B ( 22) = nz * x

      CMat_Row_B ( 23) = nz * y

      CMat_Row_B ( 24) = nz * z

      CMat_Row_B ( 25) = nz * ( x**2 - Coef_0 ( 5))

      CMat_Row_B ( 26) = nz * ( y**2 - Coef_0 ( 6))

      CMat_Row_B ( 27) = nz * ( z**2 - Coef_0 ( 7))

      CMat_Row_B ( 28) = nz * ( x * y - Coef_0 ( 8))

      CMat_Row_B ( 29) = nz * ( y * z - Coef_0 ( 9))

      CMat_Row_B ( 30) = nz * ( x * z - Coef_0 ( 10))

! ----------------------------------------------------------------------

      END SUBROUTINE CMat_Row_Coefs_P2



