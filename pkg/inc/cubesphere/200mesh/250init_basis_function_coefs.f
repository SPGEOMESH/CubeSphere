



      SUBROUTINE INIT_BASIS_FUNCTION_COEFS

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         imask_active_patch,

     1         v_quadpts_patch, v_quadwt_patch, v_volume_patch,
     1         v_centroid_patch, v_scale_length,

     1         Coef_Basis_Moments)

! ----------------------------------------------------------------------

! Written by Dinshaw Balsara

! The purpose of this routine is to build "Coef_Basis_Moments".
! These are the coefficients of the volumetric moments. They are useful
! for ensuring that all basis functions, except the zeroth order one,
! integrates to zero.

! At second order, the moments are:-

! 1 -- x^0 * y^0 * z^0 ; 2 -- x^1 * y^0 * z^0 ; 3 -- x^0 * y^1 * z^0 ;
! 4 -- x^0 * y^0 * z^1 ;

! At third order, the additional moments are:-

! 5 -- x^2 * y^0 * z^0 ; 6 -- x^0 * y^2 * z^0 ; 7 -- x^0 * y^0 * z^2 ;
! 8 -- x^1 * y^1 * z^0 ; 9 -- x^0 * y^1 * z^1 ; 10 -- x^1 * y^0 * z^1 ;

! At fourth order, the additional moments are:-

! 11 -- x^3 * y^0 * z^0 ; 12 -- x^0 * y^3 * z^0 ; 13 -- x^0 * y^0 * z^3 ;
! 14 -- x^2 * y^1 * z^0 ; 15 -- x^2 * y^0 * z^1 ; 16 -- x^1 * y^2 * z^0 ;
! 17 -- x^0 * y^2 * z^1 ; 18 -- x^1 * y^0 * z^2 ; 19 -- x^0 * y^1 * z^2 ;
! 20 -- x^1 * y^1 * z^1 ;

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_modes

      INTEGER, DIMENSION ( n_cc_modes) ::
     1         ii_mode, jj_mode, kk_mode

! ------------------------------

! Input variable, mask which is zone-centered and 1 if zone is active, 0 if inactive.
      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

! Input variable, gives us the 3-component position vector for each of the
! "NUMORDER * NUMORDER * NUMORDER" number of quadrature points in the volume.
      REAL v_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER * NUMORDER)

! Input variable, gives us the scalar weight  for each of the
! "NUMORDER * NUMORDER * NUMORDER" number of quadrature points in the volume.
      REAL v_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER * NUMORDER)

! Input variable, gives us the volume of the zone.
      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! Input variable, gives us the centroid of the zone.
      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! Input variable, gives us a scale length for each zone.
      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ------------------------------

! Output variable, these are the moments up to desired order. They are
! evaluated using the coordinates relative to the centroid of each zone.
! They are also scaled by the scale length of each zone.
      REAL Coef_Basis_Moments ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes)

! ----------------------------------------------------------------------
! Variables internal to this subroutine defined here.

      INTEGER ixx, iyy, izz, imomnt, ii, jj, kk

      REAL tempab, aa

      REAL, DIMENSION ( NUMORDER * NUMORDER * NUMORDER) ::
     1     wt, x, y, z, tempaa

! ----------------------------------------------------------------------
! Zero out "Coef_Basis_Moments". Then accumulate all the moments.

! First get coordinates of the quadrature points relative to the centroid and
! normalize them with the scale length of the zone.

! Then evaluate the weights and the moments and dot product them. That gives
! us the contribution to this moment from all the volumetric quadrature points
! in a vectorizing strategy.

! Because of the prior division by the scale length of the zone, these coefficients
! will be normalized by the length of the zone.

      DO imomnt = 1, n_cc_modes

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&                                             Coef_Basis_Moments)
!$acc loop gang vector collapse(3) independent private( ixx, iyy, izz)  &
      DO izz = 1 - NGhost, nr_zones + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
          DO ixx = 1 - NGhost, ntp_zones + NGhost

          Coef_Basis_Moments ( ixx, iyy, izz, imomnt) = 0.0

          END DO
        END DO
      END DO
!$acc end parallel
      END DO

! ------------------------------------------------------------

      DO imomnt = 2, n_cc_modes

      ii = ii_mode ( imomnt)
      jj = jj_mode ( imomnt)
      kk = kk_mode ( imomnt)

! ------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, x, y, z,                     &
!$OMP&    wt, tempaa, tempab, aa)                                       &
!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum,       &
!$OMP&    n_cc_modes,                                                   &
!$OMP&    imomnt, ii, jj, kk,                                           &

!$OMP&    ii_mode, jj_mode, kk_mode,                                    &

!$OMP&    imask_active_patch,                                           &
!$OMP&    v_quadpts_patch, v_quadwt_patch, v_volume_patch,              &
!$OMP&    v_centroid_patch, v_scale_length, Coef_Basis_Moments)         &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if( do_on_gpu) vector_length(64) present (               &
!$acc&    ii_mode, jj_mode, kk_mode,                                    &

!$acc&    imask_active_patch,                                           &
!$acc&    v_quadpts_patch, v_quadwt_patch, v_volume_patch,              &
!$acc&    v_centroid_patch, v_scale_length, Coef_Basis_Moments)         &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, x, y, z, wt, tempaa, tempab, aa)


      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

! --------------------------------------------------

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! ------------------------------

            aa = v_scale_length ( ixx, iyy, izz)

! "SIGN" is not used here, since v_scale_length is always +ve value.
            aa = 1.0 / AMAX1 ( smallnum, ABS ( aa) )

            wt ( :) = v_quadwt_patch ( ixx, iyy, izz, :)

            x ( :) = ( v_quadpts_patch ( ixx, iyy, izz, 1, :)
     1               - v_centroid_patch ( ixx, iyy, izz, 1) ) * aa

            y ( :) = ( v_quadpts_patch ( ixx, iyy, izz, 2, :)
     1               - v_centroid_patch ( ixx, iyy, izz, 2) ) * aa

            z ( :) = ( v_quadpts_patch ( ixx, iyy, izz, 3, :)
     1               - v_centroid_patch ( ixx, iyy, izz, 3) ) * aa

            tempaa ( :) = x ( :)**ii * y ( :)**jj * z ( :)**kk

#if ( NUMORDER >= 2)
             tempab = wt( 1) * tempaa ( 1) + wt( 2) * tempaa ( 2)
     1              + wt( 3) * tempaa ( 3) + wt( 4) * tempaa ( 4)
     1              + wt( 5) * tempaa ( 5) + wt( 6) * tempaa ( 6)
     1              + wt( 7) * tempaa ( 7) + wt( 8) * tempaa ( 8)
#endif
#if( NUMORDER >= 3)
     1              + wt(  9) * tempaa (  9) + wt( 10) * tempaa ( 10)
     1              + wt( 11) * tempaa ( 11) + wt( 12) * tempaa ( 12)
     1              + wt( 13) * tempaa ( 13) + wt( 14) * tempaa ( 14)
     1              + wt( 15) * tempaa ( 15) + wt( 16) * tempaa ( 16)
     1              + wt( 17) * tempaa ( 17) + wt( 18) * tempaa ( 18)
     1              + wt( 19) * tempaa ( 19) + wt( 20) * tempaa ( 20)
     1              + wt( 21) * tempaa ( 21) + wt( 22) * tempaa ( 22)
     1              + wt( 23) * tempaa ( 23) + wt( 24) * tempaa ( 24)
     1              + wt( 25) * tempaa ( 25) + wt( 26) * tempaa ( 26)
     1              + wt( 27) * tempaa ( 27)
#endif
#if( NUMORDER >= 4)
     1              + wt( 28) * tempaa ( 28) + wt( 29) * tempaa ( 29)
     1              + wt( 30) * tempaa ( 30) + wt( 31) * tempaa ( 31)
     1              + wt( 32) * tempaa ( 32) + wt( 33) * tempaa ( 33)
     1              + wt( 34) * tempaa ( 34) + wt( 35) * tempaa ( 35)
     1              + wt( 36) * tempaa ( 36) + wt( 37) * tempaa ( 37)
     1              + wt( 38) * tempaa ( 38) + wt( 39) * tempaa ( 39)
     1              + wt( 40) * tempaa ( 40) + wt( 41) * tempaa ( 41)
     1              + wt( 42) * tempaa ( 42) + wt( 43) * tempaa ( 43)
     1              + wt( 44) * tempaa ( 44) + wt( 45) * tempaa ( 45)
     1              + wt( 46) * tempaa ( 46) + wt( 47) * tempaa ( 47)
     1              + wt( 48) * tempaa ( 48) + wt( 49) * tempaa ( 49)
     1              + wt( 50) * tempaa ( 50) + wt( 51) * tempaa ( 51)
     1              + wt( 52) * tempaa ( 52) + wt( 53) * tempaa ( 53)
     1              + wt( 54) * tempaa ( 54) + wt( 55) * tempaa ( 55)
     1              + wt( 56) * tempaa ( 56) + wt( 57) * tempaa ( 57)
     1              + wt( 58) * tempaa ( 58) + wt( 59) * tempaa ( 59)
     1              + wt( 60) * tempaa ( 60) + wt( 61) * tempaa ( 61)
     1              + wt( 62) * tempaa ( 62) + wt( 63) * tempaa ( 63)
     1              + wt( 64) * tempaa ( 64)
#endif

            Coef_Basis_Moments ( ixx, iyy, izz, imomnt) =
     1      Coef_Basis_Moments ( ixx, iyy, izz, imomnt) + tempab

! ------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

! --------------------------------------------------

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------

      END DO  ! End of "DO imomnt = 1, n_cc_modes" loop.

! ----------------------------------------------------------------------

      END SUBROUTINE INIT_BASIS_FUNCTION_COEFS
