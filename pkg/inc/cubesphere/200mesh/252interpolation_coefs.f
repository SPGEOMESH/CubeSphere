



      SUBROUTINE INTERPOLATION_COEFS ( iorder, n_cc_modes, x, y, z,
     1                                 Coef_0, Amat_Row10 )

! ----------------------------------------------------------------------

! Written by Dinshaw Balsara

! This subroutine is based on the WENO-AO reconstruction for unstructured meshes in
! Balsara et al. (2020) JCP. We come in with all the coefficients for
! the Taylor basis for zone "0" and just that zone.
! Reconstruction is to be carried out in zone "0".

! "x, y, z" position of zone "j" relative to zone "0"; input. These are pre-scaled with "ell_0".
! "Coef_0" is the Taylor basis coefficients for that zone; input.

! "Amat_Row10" Coefficients of the unknowns from. These will be
! inverted to get the reconstruction coefficients from a particular stencil; output.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER iorder, n_cc_modes

      REAL x, y, z

      REAL, DIMENSION ( n_cc_modes) :: Coef_0

      REAL Amat_Row10 ( n_cc_modes)

! ----------------------------------------------------------------------
! The 3 unknowns will be "ux, uy, uz" in that order.

      IF ( iorder .EQ. 2) THEN

      Amat_Row10 ( 1) = 1.0

      Amat_Row10 ( 2) = x

      Amat_Row10 ( 3) = y

      Amat_Row10 ( 4) = z

      END IF  ! End of "IF ( iorder .EQ. 2)".

! ----------------------------------------------------------------------
! The 9 unknowns will be "ux, uy, uz, uxx, uyy, uzz, uxy, uyz, uxz" in that order.

      IF ( iorder .EQ. 3) THEN

      Amat_Row10 ( 1) = 1.0

      Amat_Row10 ( 2) = x

      Amat_Row10 ( 3) = y

      Amat_Row10 ( 4) = z

      Amat_Row10 ( 5) = x**2 - Coef_0 ( 5)

      Amat_Row10 ( 6) = y**2 - Coef_0 ( 6)

      Amat_Row10 ( 7) = z**2 - Coef_0 ( 7)

      Amat_Row10 ( 8) = x * y - Coef_0 ( 8)

      Amat_Row10 ( 9) = y * z - Coef_0 ( 9)

      Amat_Row10 ( 10) = x * z - Coef_0 ( 10)

      END IF  ! End of "IF ( iorder .EQ. 3)".

! ----------------------------------------------------------------------
! The 19 unknowns will be "ux, uy, uz, uxx, uyy, uzz, uxy, uyz, uxz,
! uxxx, uyyy, uzzz, uxxy, uxxz, uxyy, uyyz, uxzz, uyzz, uxyz" in that order.

      IF ( iorder .EQ. 4) THEN

      Amat_Row10 ( 1) = 1.0

      Amat_Row10 ( 2) = x

      Amat_Row10 ( 3) = y

      Amat_Row10 ( 4) = z

      Amat_Row10 ( 5) = x**2 - Coef_0 ( 5)

      Amat_Row10 ( 6) = y**2 - Coef_0 ( 6)

      Amat_Row10 ( 7) = z**2 - Coef_0 ( 7)

      Amat_Row10 ( 8) = x * y - Coef_0 ( 8)

      Amat_Row10 ( 9) = y * z - Coef_0 ( 9)

      Amat_Row10 ( 10) = x * z - Coef_0 ( 10)

      Amat_Row10 ( 11) = x**3 - Coef_0 ( 11)

      Amat_Row10 ( 12) = y**3 - Coef_0 ( 12)

      Amat_Row10 ( 13) = z**3 - Coef_0 ( 13)

      Amat_Row10 ( 14) = x**2 * y - Coef_0 ( 14)

      Amat_Row10 ( 15) = x**2 * z - Coef_0 ( 15)

      Amat_Row10 ( 16) = x * y**2 - Coef_0 ( 16)

      Amat_Row10 ( 17) = y**2 * z - Coef_0 ( 17)

      Amat_Row10 ( 18) = x * z**2 - Coef_0 ( 18)

      Amat_Row10 ( 19) = y * z**2 - Coef_0 ( 19)

      Amat_Row10 ( 20) = x * y * z - Coef_0 ( 20)

      END IF  ! End of "IF ( iorder .EQ. 4)"

! ----------------------------------------------------------------------

      END SUBROUTINE INTERPOLATION_COEFS

