





      SUBROUTINE BUILD_MODE_INDICES

     1         ( n_cc_modes,

     1           ii_mode, jj_mode, kk_mode)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Builds the indices that enable us to zip through the zone-centered and
! face-centered reconstructions.

! ----------------------------------------------------------------------

      INTEGER n_cc_modes,
     1        ii, jj, kk, mm, ll

      INTEGER, DIMENSION ( n_cc_modes) ::
     1         ii_mode, jj_mode, kk_mode

! ----------------------------------------------------------------------

! Build the powers for the polynomials so that we can zip through them
! automatically during stencil construction.

! At second order, we minimally want the constant term and the first moments.
! We will need the quadratic moments for the magnetic field.

! 1 -- x^0 * y^0 * z^0 ; 2 -- x^1 * y^0 * z^0 ; 3 -- x^0 * y^1 * z^0 ;
! 4 -- x^0 * y^0 * z^1 ;

! ----------------------------------------

      ii_mode ( 01) = 0
      jj_mode ( 01) = 0
      kk_mode ( 01) = 0

      ii_mode ( 02) = 1
      jj_mode ( 02) = 0
      kk_mode ( 02) = 0

      ii_mode ( 03) = 0
      jj_mode ( 03) = 1
      kk_mode ( 03) = 0

      ii_mode ( 04) = 0
      jj_mode ( 04) = 0
      kk_mode ( 04) = 1

! ----------------------------------------

! At third order, we also want the quadratic moments.

! 5 -- x^2 * y^0 * z^0 ; 6 -- x^0 * y^2 * z^0 ; 7 -- x^0 * y^0 * z^2 ;
! 8 -- x^1 * y^1 * z^0 ; 9 -- x^0 * y^1 * z^1 ; 10 -- x^1 * y^0 * z^1 ;

! ----------------------------------------

#if ( NUMORDER >= 3)

      ii_mode ( 05) = 2
      jj_mode ( 05) = 0
      kk_mode ( 05) = 0

      ii_mode ( 06) = 0
      jj_mode ( 06) = 2
      kk_mode ( 06) = 0

      ii_mode ( 07) = 0
      jj_mode ( 07) = 0
      kk_mode ( 07) = 2

      ii_mode ( 08) = 1
      jj_mode ( 08) = 1
      kk_mode ( 08) = 0

      ii_mode ( 09) = 0
      jj_mode ( 09) = 1
      kk_mode ( 09) = 1

      ii_mode ( 10) = 1
      jj_mode ( 10) = 0
      kk_mode ( 10) = 1

! ----------------------------------------

#endif

! ----------------------------------------

#if ( NUMORDER >= 4)

! At fourth order, we will also need all the cubic moments.

! 11 -- x^3 * y^0 * z^0 ; 12 -- x^0 * y^3 * z^0 ; 13 -- x^0 * y^0 * z^3 ;
! 14 -- x^2 * y^1 * z^0 ; 15 -- x^2 * y^0 * z^1 ; 16 -- x^1 * y^2 * z^0 ;
! 17 -- x^0 * y^2 * z^1 ; 18 -- x^1 * y^0 * z^2 ; 19 -- x^0 * y^1 * z^2 ;
! 20 -- x^1 * y^1 * z^1 ;

! ----------------------------------------

      ii_mode ( 11) = 3
      jj_mode ( 11) = 0
      kk_mode ( 11) = 0

      ii_mode ( 12) = 0
      jj_mode ( 12) = 3
      kk_mode ( 12) = 0

      ii_mode ( 13) = 0
      jj_mode ( 13) = 0
      kk_mode ( 13) = 3

      ii_mode ( 14) = 2
      jj_mode ( 14) = 1
      kk_mode ( 14) = 0

      ii_mode ( 15) = 2
      jj_mode ( 15) = 0
      kk_mode ( 15) = 1

      ii_mode ( 16) = 1
      jj_mode ( 16) = 2
      kk_mode ( 16) = 0

      ii_mode ( 17) = 0
      jj_mode ( 17) = 2
      kk_mode ( 17) = 1

      ii_mode ( 18) = 1
      jj_mode ( 18) = 0
      kk_mode ( 18) = 2

      ii_mode ( 19) = 0
      jj_mode ( 19) = 1
      kk_mode ( 19) = 2

      ii_mode ( 20) = 1
      jj_mode ( 20) = 1
      kk_mode ( 20) = 1

! ----------------------------------------

#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END SUBROUTINE BUILD_MODE_INDICES

