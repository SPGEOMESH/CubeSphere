




      SUBROUTINE BUILD_NON_STIFF_SOURCES

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,

     1           iproblem, time, dtmove,

     1           vert_coords_patch, imask_active_patch,

     1           v_scale_length,
     1           v_volume_patch, v_centroid_patch,

     1           U, V, Source, dU_dt)

! ----------------------------------------------------------------------
! Written By: Sethupathy S.

! Builds the non-stiff source terms, that will used in the RK-update step.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        iproblem

      REAL time, dtmove

! -----

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

! -----

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! -----

      REAL U ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL V ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost,
     1                                 n_cc_components) :: Source

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost,
     1                                 n_cc_components) :: dU_dt

! ----------------------------------------------------------------------
! Local variables.

      INTEGER ixx, iyy, izz, iexpand

      REAL rho_cen, obl0, rad, ob_rad3, GMbr3, aa,

     1     x_cen, vlx_cen, om_x, om_cross_r_x, centri_x, coriolis_x,
     1     y_cen, vly_cen, om_y, om_cross_r_y, centri_y, coriolis_y,
     1     z_cen, vlz_cen, om_z, om_cross_r_z, centri_z, coriolis_z,

     1     S_cen ( n_cc_components)

! ----------------------------------------------------------------------

!$acc kernels
      Source ( :, :, :, :) = 0.0
!$acc end kernels

      CALL MAKE_EXTERNAL_SOURCE

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,

     1           iproblem, time, dtmove,

     1           vert_coords_patch, imask_active_patch,

     1           v_scale_length,
     1           v_volume_patch, v_centroid_patch,

     1           U, V, Source)

! ----------------------------------------------------------------------

      iexpand = 1  ! Expand outwards by one zone; could be two if we wish.

!$acc data copyin ( iexpand)

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz,                              &
!$OMP&    rho_cen, obl0, rad, ob_rad3, GMbr3, aa,                       &
!$OMP&    x_cen, vlx_cen, om_x, om_cross_r_x, centri_x, coriolis_x,     &
!$OMP&    y_cen, vly_cen, om_y, om_cross_r_y, centri_y, coriolis_y,     &
!$OMP&    z_cen, vlz_cen, om_z, om_cross_r_z, centri_z, coriolis_z,     &
!$OMP&    S_cen)                                                        &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost, iexpand,                &
!$OMP&    smallnum, newtonsc, Mstar, omega_passive,                     &
!$OMP&    vert_coords_patch, imask_active_patch, v_scale_length,        &
!$OMP&    v_volume_patch, v_centroid_patch, U, V, Source, dU_dt)        &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 3)

!$acc parallel vector_length(64) present( iexpand,                      &
!$acc&    smallnum, newtonsc, Mstar, omega_passive,                     &
!$acc&    vert_coords_patch, imask_active_patch, v_scale_length,        &
!$acc&    v_volume_patch, v_centroid_patch, U, V, Source, dU_dt)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz,                                      &
!$acc&    rho_cen, obl0, rad, ob_rad3, GMbr3, aa,                       &
!$acc&    x_cen, vlx_cen, om_x, om_cross_r_x, centri_x, coriolis_x,     &
!$acc&    y_cen, vly_cen, om_y, om_cross_r_y, centri_y, coriolis_y,     &
!$acc&    z_cen, vlz_cen, om_z, om_cross_r_z, centri_z, coriolis_z,     &
!$acc&    S_cen)

      DO izz = 1 - iexpand,  nr_zones + iexpand
      DO iyy = 1 - iexpand, ntp_zones + iexpand
      DO ixx = 1 - iexpand, ntp_zones + iexpand

      IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! ----------------------------------------------------------------------

      x_cen = v_centroid_patch ( ixx, iyy, izz, 1)
      y_cen = v_centroid_patch ( ixx, iyy, izz, 2)
      z_cen = v_centroid_patch ( ixx, iyy, izz, 3)

      rho_cen = V ( ixx, iyy, izz, 1, 1)

      vlx_cen = V ( ixx, iyy, izz, 2, 1)
      vly_cen = V ( ixx, iyy, izz, 3, 1)
      vlz_cen = V ( ixx, iyy, izz, 4, 1)

! ----------------------------------------

      obl0 = 1.0 / MAX ( v_scale_length ( ixx, iyy, izz), smallnum)

! Use the above time-centered primitive variables to re-evaluate the
! source terms "S_cen". This time, the evaluation is at a time of t=1/2.

      S_cen = 0.0

! Gravitation force of the star.

      rad = SQRT ( x_cen**2 + y_cen**2 + z_cen**2)

      ob_rad3 = 1.0 / rad**3
      GMbr3 = newtonsc * Mstar * ob_rad3

! "omega_passive" is the rotational velocity of the frame, and not of the star.

      om_x = omega_passive ( 1)
      om_y = omega_passive ( 2)
      om_z = omega_passive ( 3)

      aa = SQRT ( om_x**2 + om_y**2 + om_z**2)

! --------------------

      IF ( (Mstar > smallnum) .OR. (aa > smallnum) ) THEN

! ----------
! Centrifugal force components.
      om_cross_r_x = om_y * z_cen - om_z * y_cen
      om_cross_r_y = om_z * x_cen - om_x * z_cen
      om_cross_r_z = om_x * y_cen - om_y * x_cen

      centri_x = - rho_cen * (om_y * om_cross_r_z - om_z * om_cross_r_y)
      centri_y = - rho_cen * (om_z * om_cross_r_x - om_x * om_cross_r_z)
      centri_z = - rho_cen * (om_x * om_cross_r_y - om_y * om_cross_r_x)

! ----------
! Coriolis force components.
      coriolis_x = - 2.0 * rho_cen * ( om_y * vlz_cen - om_z * vly_cen)
      coriolis_y = - 2.0 * rho_cen * ( om_z * vlx_cen - om_x * vlz_cen)
      coriolis_z = - 2.0 * rho_cen * ( om_x * vly_cen - om_y * vlx_cen)

! ----------

      S_cen ( 2) = - GMbr3 * x_cen * rho_cen + centri_x + coriolis_x
      S_cen ( 3) = - GMbr3 * y_cen * rho_cen + centri_y + coriolis_y
      S_cen ( 4) = - GMbr3 * z_cen * rho_cen + centri_z + coriolis_z

#if ( CHOOSE_ADIABATIC == 1)
      S_cen ( 5) = vlx_cen * S_cen ( 2) + vly_cen * S_cen ( 3)
     1           + vlz_cen * S_cen ( 4)
#else
      S_cen ( 5) = 0.0 ! Energy source is zero since it is isothermal.
#endif

! ----------

      END IF ! End of "IF ( (Mstar > smallnum) .OR. (aa > smallnum) )".

! ----------------------------------------

      Source ( ixx, iyy, izz, :) = Source ( ixx, iyy, izz, :)
     1                           + S_cen ( :)

! --------------------
! Add explicit non-stiff source terms to "dU_dt", in addition to the
! existing fluctuation terms.

      dU_dt ( ixx, iyy, izz, :) = dU_dt  ( ixx, iyy, izz, :)
     1                          + Source ( ixx, iyy, izz, :)

! ----------------------------------------------------------------------

      END IF ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN".

      END DO ! End of "DO ixx = 1 - iexpand, ntp_zones + iexpand".
      END DO ! End of "DO iyy = 1 - iexpand, ntp_zones + iexpand".
      END DO ! End of "DO izz = 1 - iexpand,  nr_zones + iexpand".

!$acc end parallel
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! The number of "end data" statements must match the "data copyin" statements.

!$acc end data

      END SUBROUTINE BUILD_NON_STIFF_SOURCES

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
#if ( NUMORDER == 2)








