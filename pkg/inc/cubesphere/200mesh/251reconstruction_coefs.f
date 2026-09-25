

      SUBROUTINE RECONSTRUCTION_COEFS ( iorder, n_cc_modes, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

! ----------------------------------------------------------------------

! Written by Dinshaw Balsara

! This subroutine is based on the WENO-AO for unstructured meshes in
! Balsara et al. (2020) JCP. We come in with all the coefficients for
! the Taylor basis for zone "0" and zone "j". The target zone where the
! reconstruction is to be carried out is zone "0". Zone "j" enables us to
! fulfil eqns. (2.15) or (3.9) at third and fourth orders.

! "x, y, z" position of zone "j" relative to zone "0"; input. These are pre-scaled with "ell_0".
! "ell_0, ell_j" are the scaled lengths for the two zones, input.
! "Coef_0, Coef_j" are the Taylor basis coefficients for those two zones; input.

! "Amat_Row" Coefficients of the unknowns from (2.15) or (3.9). These will be
! inverted to get the reconstruction coefficients from a particular stencil; output.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER iorder, n_cc_modes

      REAL x, y, z, ell_0, ell_j, ell_j_by_ell_0

      REAL, DIMENSION ( n_cc_modes) :: Coef_0, Coef_j

      REAL Amat_Row ( n_cc_modes - 1)

! ----------------------------------------------------------------------

      ell_j_by_ell_0 = ell_j / ell_0

! ----------------------------------------------------------------------
! The 3 unknowns will be "ux, uy, uz" in that order.

      IF ( iorder .EQ. 2) THEN

      Amat_Row ( 1) = x

      Amat_Row ( 2) = y

      Amat_Row ( 3) = z

      END IF  ! End of "IF ( iorder .EQ. 2)".

! ----------------------------------------------------------------------
! The 9 unknowns will be "ux, uy, uz, uxx, uyy, uzz, uxy, uyz, uxz" in that order.

      IF ( iorder .EQ. 3) THEN

      Amat_Row ( 1) = x

      Amat_Row ( 2) = y

      Amat_Row ( 3) = z

      Amat_Row ( 4) = x**2 - Coef_0 ( 5)
     1              + Coef_j ( 5) * ell_j_by_ell_0**2

      Amat_Row ( 5) = y**2 - Coef_0 ( 6)
     1              + Coef_j ( 6) * ell_j_by_ell_0**2

      Amat_Row ( 6) = z**2 - Coef_0 ( 7)
     1              + Coef_j ( 7) * ell_j_by_ell_0**2

      Amat_Row ( 7) = x * y - Coef_0 ( 8)
     1              + Coef_j ( 8) * ell_j_by_ell_0**2

      Amat_Row ( 8) = y * z - Coef_0 ( 9)
     1              + Coef_j ( 9) * ell_j_by_ell_0**2

      Amat_Row ( 9) = x * z - Coef_0 ( 10)
     1              + Coef_j ( 10) * ell_j_by_ell_0**2

      END IF  ! End of "IF ( iorder .EQ. 3)".

! ----------------------------------------------------------------------
! The 19 unknowns will be "ux, uy, uz, uxx, uyy, uzz, uxy, uyz, uxz,
! uxxx, uyyy, uzzz, uxxy, uxxz, uxyy, uyyz, uxzz, uyzz, uxyz" in that order.

      IF ( iorder .EQ. 4) THEN

      Amat_Row ( 1) = x

      Amat_Row ( 2) = y

      Amat_Row ( 3) = z

      Amat_Row ( 4) = x**2 - Coef_0 ( 5)
     1              + Coef_j ( 5) * ell_j_by_ell_0**2

      Amat_Row ( 5) = y**2 - Coef_0 ( 6)
     1              + Coef_j ( 6) * ell_j_by_ell_0**2

      Amat_Row ( 6) = z**2 - Coef_0 ( 7)
     1              + Coef_j ( 7) * ell_j_by_ell_0**2

      Amat_Row ( 7) = x * y - Coef_0 ( 8)
     1              + Coef_j ( 8) * ell_j_by_ell_0**2

      Amat_Row ( 8) = y * z - Coef_0 ( 9)
     1              + Coef_j ( 9) * ell_j_by_ell_0**2

      Amat_Row ( 9) = x * z - Coef_0 ( 10)
     1              + Coef_j ( 10) * ell_j_by_ell_0**2

      Amat_Row ( 10) = x**3 - Coef_0 ( 11)
     1               + 3.0 * Coef_j ( 5) * x * ell_j_by_ell_0**2
     1               + Coef_j ( 11) * ell_j_by_ell_0**3

      Amat_Row ( 11) = y**3 - Coef_0 ( 12)
     1               + 3.0 * Coef_j ( 6) * y * ell_j_by_ell_0**2
     1               + Coef_j ( 12) * ell_j_by_ell_0**3

      Amat_Row ( 12) = z**3 - Coef_0 ( 13)
     1               + 3.0 * Coef_j ( 7) * z * ell_j_by_ell_0**2
     1               + Coef_j ( 13) * ell_j_by_ell_0**3

      Amat_Row ( 13) = x**2 * y - Coef_0 ( 14)
     1               + Coef_j ( 5) * y * ell_j_by_ell_0**2
     1               + 2.0 * Coef_j ( 8) * x * ell_j_by_ell_0**2
     1               + Coef_j ( 14) * ell_j_by_ell_0**3

      Amat_Row ( 14) = x**2 * z - Coef_0 ( 15)
     1               + Coef_j ( 5) * z * ell_j_by_ell_0**2
     1               + 2.0 * Coef_j ( 10) * x * ell_j_by_ell_0**2
     1               + Coef_j ( 15) * ell_j_by_ell_0**3

      Amat_Row ( 15) = x * y**2 - Coef_0 ( 16)
     1               + Coef_j ( 6) * x * ell_j_by_ell_0**2
     1               + 2.0 * Coef_j ( 8) * y * ell_j_by_ell_0**2
     1               + Coef_j ( 16) * ell_j_by_ell_0**3

      Amat_Row ( 16) = y**2 * z - Coef_0 ( 17)
     1               + Coef_j ( 6) * z * ell_j_by_ell_0**2
     1               + 2.0 * Coef_j ( 9) * y * ell_j_by_ell_0**2
     1               + Coef_j ( 17) * ell_j_by_ell_0**3

      Amat_Row ( 17) = x * z**2 - Coef_0 ( 18)
     1               + Coef_j ( 7) * x * ell_j_by_ell_0**2
     1               + 2.0 * Coef_j ( 10) * z * ell_j_by_ell_0**2
     1               + Coef_j ( 18) * ell_j_by_ell_0**3

      Amat_Row ( 18) = y * z**2 - Coef_0 ( 19)
     1               + Coef_j ( 7) * y * ell_j_by_ell_0**2
     1               + 2.0 * Coef_j ( 9) * z * ell_j_by_ell_0**2
     1               + Coef_j ( 19) * ell_j_by_ell_0**3

      Amat_Row ( 19) = x * y * z - Coef_0 ( 20)
     1               + Coef_j ( 8) * z * ell_j_by_ell_0**2
     1               + Coef_j ( 10) * y * ell_j_by_ell_0**2
     1               + Coef_j ( 9) * x * ell_j_by_ell_0**2
     1               + Coef_j ( 20) * ell_j_by_ell_0**3

      END IF  ! End of "IF ( iorder .EQ. 4)"

! ----------------------------------------------------------------------

      END SUBROUTINE RECONSTRUCTION_COEFS
