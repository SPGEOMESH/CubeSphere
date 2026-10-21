






      SUBROUTINE INIT_PATCH_CUBE_SPHERE

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           iproblem, time, dtcur,
     1           rmin, rmax,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_quadpts_patch, v_quadwt_patch,

     1           U, V, entropygr, Bxf, Byf, Bzf,
     1           B0field,

     1           gamma_eos_gr,

     1           V_min_local, V_max_local)

! ----------------------------------------------------------------------

! Written by : Sethupathy S., Saurav Samantaray

! Initialize the conserved variable "U" using the volume integration of the
! pointwise values obtained from the "PROBLEM_GENERATOR_U_A" subroutine
! at each "x, y, z" locations.

! Then, initialize the "Bxf, Byf, Bzf" with their respective face centered
! values, by doing a second order accurate evaluation of the closed loop
! line integration of the vector potential - A.

! ----------------------------------------------------------------------

#if ( CHOOSE_APPL_HELIO == 1)
      USE HELIOSPHERE_MOD
#endif

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes

      INTEGER iproblem

      REAL time, dtcur, rmin, rmax

! ----------------------------------------------------------------------
! Arrays describing the geometry of the mesh.

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ----------------------------------------------------------------------
! Arrays used for curl-A calculation.

! ------------------------------
! FACIAL VARIABLES (set in COORDS_TO_WT_QUADPTS_FACES_ISO2):-

      REAL xf_quadpts_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL xf_quadwt_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NAREA_QUAD_PTS)

      REAL xf_quadnorml_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL xf_area_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost),
     1     xf_centroid_patch ( 0 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

! -----

      REAL yf_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL yf_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NAREA_QUAD_PTS)

      REAL yf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL yf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost),
     1     yf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

! -----

      REAL zf_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL zf_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   NAREA_QUAD_PTS)

      REAL zf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL zf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost),
     1     zf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost, 3)

! ----------------------------------------------------------------------
! Arrays used for volume integration.

      REAL v_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER * NUMORDER)

      REAL v_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER * NUMORDER)

! ----------------------------------------------------------------------
! Arrays representing the MHD variables in volume and on faces.

      REAL U ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL V ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components)

      REAL entropygr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   4)

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           n_fc_modes) ::
     1   Bxf, Byf, Bzf

      REAL B0field ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost) ::
     1      gamma_eos_gr

      REAL, DIMENSION ( n_cc_components + 1) :: V_min_local, V_max_local

! ----------------------------------------------------------------------
! Variables decalred below are local to this subroutine.

      INTEGER ixx, iyy, izz, ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,
     1        i, ii

! Used in Problem Generator part.
      REAL x, y, z, zone_size, xcen_sector, ycen_sector, zcen_sector,
     1     xvert1_sector, yvert1_sector, zvert1_sector,
     1     U_pt ( n_cc_components), A_pt ( 3)

! Local real variables used in the facial magnetic field assignement.
      REAL v1 ( 3), v2 ( 3), v3 ( 3), v4 ( 3),
     1     v1_mag,  v2_mag,  v3_mag,  v4_mag,

     1      v1_dot_v2, v3_dot_v4, v1_dot_v3, v2_dot_v4, angle, sum_A_dl

      REAL, DIMENSION ( 3) ::
     1      e1_cen, e2_cen, e3_cen, e4_cen,
     1      e1_dl,  e2_dl,  e3_dl,  e4_dl, uvec

! Variables V-min and V-max and dtcur
      REAL a11, a12, a13, a21, ob_8pi,
     1     tempaa, tempab, tempac, tempad, vx1, vy1, vz1, dr,

     1     rhomin, rhomax, vlxmin, vlxmax, vlymin, vlymax,
     1     vlzmin, vlzmax, prsmin, prsmax, bfxmin, bfxmax,
     1     bfymin, bfymax, bfzmin, bfzmax, div_min, div_max,
     1     flat_min, flat_max

! ----------------------------------------------------------------------
! Variables required by the ProblemGen_include.

      REAL :: rad, costheta, sintheta, cosphi, sinphi, R,
     1        rho, prs, v0, vr, vtheta, vphi, vx, vy, vz, bfx, bfy, bfz,
     1        rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0,
     1        rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,
     1        vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,
     1        radius, distance, bmag,
     1        m1, m2, m3, phi, theta

      REAL :: d0, u0, p0, gammaa, d, uu, p, a, b,
     1        gamma1, dcur, dprev, c, f1, f2, f, df,
     1        kappa, sq4pi

      REAL :: alpha, beta, tau, r0, omega, mach

      REAL :: Ar, Br, Br0, Bphi, Bphi1, Bphi2, eps,
     1        Aphi, B0, band_width, cosomegat, sinomegat,
     1        cosphi_mag, sinphi_mag, costheta_mag, sintheta_mag,
     1        costilt, sintilt, fast_slow, tilt,
     1        theta_mag, phi_mag, trans_width, x_mag, y_mag, z_mag

! ----------------------------------------------------------------------

!$acc routine( GET_MAGNETIC_FIELD ) seq
!$acc routine( PROBLEM_GENERATOR_U_A ) seq

! ----------------------------------------------------------------------
! initialize ratios.

      ob_8pi = 1.0 / ( 8.0 * pi)

!$acc update device ( time)

! ----------------------------------------------------------------------
! STEP 1. Initializing patch/sector related variables, this is similar
! to triangular Geomesh code.

! Needed in the Rotor problem for constructing the taper.
      zone_size = ABS( SQRT( vert_coords_patch ( 1, 1, nr_zones, 1)**2
     1                     + vert_coords_patch ( 1, 1, nr_zones, 2)**2
     1                     + vert_coords_patch ( 1, 1, nr_zones, 3)**2)
     1                - SQRT( vert_coords_patch ( 1, 1, 1, 1)**2
     1                      + vert_coords_patch ( 1, 1, 1, 2)**2
     1                      + vert_coords_patch ( 1, 1, 1, 3)**2) )
     1          / nr_zones

       xvert1_sector = vert_coords_patch ( 1, 1, 1, 1)
       yvert1_sector = vert_coords_patch ( 1, 1, 1, 2)
       zvert1_sector = vert_coords_patch ( 1, 1, 1, 3)

! Assuming a near rectangular patch for the test problems. Initializing
! the approximate x, y and z center of the patch.
      xcen_sector = 0.125
     1      * ( vert_coords_patch ( 1, 1, 1, 1)
     1        + vert_coords_patch ( ntp_zones, 1, 1, 1)
     1        + vert_coords_patch ( 1, ntp_zones, 1, 1)
     1        + vert_coords_patch ( ntp_zones, ntp_zones, 1, 1)
     1        + vert_coords_patch ( 1, 1, nr_zones, 1)
     1        + vert_coords_patch ( ntp_zones, 1, nr_zones, 1)
     1        + vert_coords_patch ( 1, ntp_zones, nr_zones, 1)
     1        + vert_coords_patch ( ntp_zones, ntp_zones, nr_zones, 1))

      ycen_sector = 0.125
     1      * ( vert_coords_patch ( 1, 1, 1, 2)
     1        + vert_coords_patch ( ntp_zones, 1, 1, 2)
     1        + vert_coords_patch ( 1, ntp_zones, 1, 2)
     1        + vert_coords_patch ( ntp_zones, ntp_zones, 1, 2)
     1        + vert_coords_patch ( 1, 1, nr_zones, 2)
     1        + vert_coords_patch ( ntp_zones, 1, nr_zones, 2)
     1        + vert_coords_patch ( 1, ntp_zones, nr_zones, 2)
     1        + vert_coords_patch ( ntp_zones, ntp_zones, nr_zones, 2))

      zcen_sector = 0.125
     1      * ( vert_coords_patch ( 1, 1, 1, 3)
     1        + vert_coords_patch ( ntp_zones, 1, 1, 3)
     1        + vert_coords_patch ( 1, ntp_zones, 1, 3)
     1        + vert_coords_patch ( ntp_zones, ntp_zones, 1, 3)
     1        + vert_coords_patch ( 1, 1, nr_zones, 3)
     1        + vert_coords_patch ( ntp_zones, 1, nr_zones, 3)
     1        + vert_coords_patch ( 1, ntp_zones, nr_zones, 3)
     1        + vert_coords_patch ( ntp_zones, ntp_zones, nr_zones, 3))

!$acc data copyin ( ob_8pi, zone_size, xvert1_sector, yvert1_sector,    &
!$acc&           zvert1_sector, xcen_sector, ycen_sector, zcen_sector)

! ----------------------------------------------------------------------
! STEP 2. Initializing the zone-centered values of "U", by integrating the
! pointwise values over the volume.

! ------------------------------------------------------------
! zero-out the "U" to facilitate the accumulation of weight*value.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost, U, B0field )            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present ( U, B0field )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)
      DO izz = 1 - NGhost, nr_zones + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
          DO ixx = 1 - NGhost, ntp_zones + NGhost
            U ( ixx, iyy, izz, :, :) = 0.0
            B0field ( ixx, iyy, izz, :) = 0.0
          END DO
        END DO
      END DO
!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Volume integration of point-wise "U".

      DO i = 1, NUMORDER * NUMORDER * NUMORDER

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz,                              &
!$OMP&    x, y, z, U_pt, A_pt,                                          &
!$OMP&    rad, costheta, sintheta, cosphi, sinphi, R,                   &
!$OMP&    rho, prs, v0, vr, vtheta, vphi, vx, vy, vz, bfx, bfy, bfz,    &
!$OMP&    rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0,               &
!$OMP&    rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,               &
!$OMP&    vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,               &
!$OMP&    radius, distance, bmag, m1, m2, m3, phi, theta)

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components, i,     &
!$OMP&    pi, gamma, rhofloor, prsfloor, ob_8pi, entropygr,             &
!$OMP&    isospeed, gamma_euler, barotropic_stiffness,                  &
!$OMP&    time, iproblem, rmin, rmax,                                   &
!$OMP&    zone_size, xcen_sector, ycen_sector, zcen_sector,             &
!$OMP&    xvert1_sector, yvert1_sector, zvert1_sector, B0field,         &
!$OMP&    imask_active_patch, v_quadpts_patch, v_quadwt_patch, U )      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    pi, gamma, rhofloor, prsfloor, ob_8pi, entropygr,             &
!$acc&    isospeed, gamma_euler, barotropic_stiffness,                  &
!$acc&    time, iproblem, rmin, rmax,                                   &
!$acc&    zone_size, xcen_sector, ycen_sector, zcen_sector,             &
!$acc&    xvert1_sector, yvert1_sector, zvert1_sector, B0field,         &
!$acc&    imask_active_patch, v_quadpts_patch, v_quadwt_patch, U )      &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz,                                      &
!$acc&    x, y, z, U_pt, A_pt,                                          &
!$acc&    rad, costheta, sintheta, cosphi, sinphi, R,                   &
!$acc&    rho, prs, v0, vr, vtheta, vphi, vx, vy, vz, bfx, bfy, bfz,    &
!$acc&    rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0,               &
!$acc&    rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,               &
!$acc&    vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,               &
!$acc&    radius, distance, bmag, m1, m2, m3, phi, theta)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! ----------------------------------------

      x = v_quadpts_patch ( ixx, iyy, izz, 1, i)
      y = v_quadpts_patch ( ixx, iyy, izz, 2, i)
      z = v_quadpts_patch ( ixx, iyy, izz, 3, i)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      U ( ixx, iyy, izz, :, 1) = U ( ixx, iyy, izz, :, 1)
     1             + U_pt ( :) * v_quadwt_patch ( ixx, iyy, izz, i)

! ----------------------------------------
! Initialize background magnetic field. Used only for imaging.

#if( USE_B0B1_FORMULATION == 1)

      CALL GET_MAGNETIC_FIELD ( 1, 0.0, x, y, z, bfx, bfy, bfz )

      B0field ( ixx, iyy, izz, 1) = B0field ( ixx, iyy, izz, 1)
     1               + bfx * v_quadwt_patch ( ixx, iyy, izz, i)

      B0field ( ixx, iyy, izz, 2) = B0field ( ixx, iyy, izz, 2)
     1               + bfy * v_quadwt_patch ( ixx, iyy, izz, i)

      B0field ( ixx, iyy, izz, 3) = B0field ( ixx, iyy, izz, 3)
     1               + bfz * v_quadwt_patch ( ixx, iyy, izz, i)

#endif

! ----------------------------------------
! Initialize "entropygr".

#if ( PRESSURE_POSITIVITY >= 1)

          rho = MAX ( U ( ixx, iyy, izz, 1, 1), rhofloor )

          v0  = rho / ( rho**2 + rhofloor**2)

          vx  = U ( ixx, iyy, izz, 2, 1) * v0
          vy  = U ( ixx, iyy, izz, 3, 1) * v0
          vz  = U ( ixx, iyy, izz, 4, 1) * v0

          bfx = U ( ixx, iyy, izz, 6, 1)
          bfy = U ( ixx, iyy, izz, 7, 1)
          bfz = U ( ixx, iyy, izz, 8, 1)

#if ( CHOOSE_ADIABATIC == 1)
          prs = ( gamma - 1.0) * ( U ( ixx, iyy, izz, 5, 1)
     1        - 0.5 * rho * ( vx**2 + vy**2 + vz**2)
     1        - ( bfx**2 + bfy**2 + bfz**2) * ob_8pi )
#else
          prs = rho * isospeed**2
     1        + barotropic_stiffness * rho**gamma_euler
#endif

          prs = MAX ( prs, prsfloor)

          entropygr ( ixx, iyy, izz, 1) = prs / rho**( gamma - 1.0)

#endif

! ----------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) THEN".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO i = 1, NUMORDER * NUMORDER" loop.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 3.a) Initializing the x-face magnetic field.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxp1, iyyp1, izzp1, ii,     &
!$OMP&     x, y, z, U_pt, A_pt,                                         &

!$OMP&     v1, v2, v3, v4, v1_mag,  v2_mag,  v3_mag,  v4_mag,           &
!$OMP&     v1_dot_v2, v3_dot_v4, v1_dot_v3, v2_dot_v4, angle, sum_A_dl, &
!$OMP&     e1_cen, e2_cen, e3_cen, e4_cen,                              &
!$OMP&     e1_dl,  e2_dl,  e3_dl,  e4_dl, uvec,                         &
!$OMP&    rad, costheta, sintheta, cosphi, sinphi, R,                   &
!$OMP&    rho, prs, v0, vr, vtheta, vphi, vx, vy, vz, bfx, bfy, bfz,    &
!$OMP&    rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0,               &
!$OMP&    rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,               &
!$OMP&    vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,               &
!$OMP&    radius, distance, bmag, m1, m2, m3, phi, theta)

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,        &
!$OMP&    pi, gamma,                                                    &
!$OMP&    time, iproblem, rmin, rmax,                                   &
!$OMP&    zone_size, xcen_sector, ycen_sector, zcen_sector,             &
!$OMP&    xvert1_sector, yvert1_sector, zvert1_sector,                  &
!$OMP&    vert_coords_patch, imask_active_patch,                        &
!$OMP&    Bxf, xf_area_patch)                                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    vert_coords_patch, imask_active_patch,                        &
!$acc&    Bxf, xf_area_patch)                                           &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxp1, iyyp1, izzp1, ii,             &
!$acc&     x, y, z, U_pt, A_pt,                                         &

!$acc&     v1, v2, v3, v4, v1_mag,  v2_mag,  v3_mag,  v4_mag,           &
!$acc&     v1_dot_v2, v3_dot_v4, v1_dot_v3, v2_dot_v4, angle, sum_A_dl, &
!$acc&     e1_cen, e2_cen, e3_cen, e4_cen,                              &
!$acc&     e1_dl,  e2_dl,  e3_dl,  e4_dl, uvec,                         &
!$acc&    rad, costheta, sintheta, cosphi, sinphi, R,                   &
!$acc&    rho, prs, v0, vr, vtheta, vphi, vx, vy, vz, bfx, bfy, bfz,    &
!$acc&    rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0,               &
!$acc&    rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,               &
!$acc&    vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,               &
!$acc&    radius, distance, bmag, m1, m2, m3, phi, theta)


      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 0 - NGhost, ntp_zones + NGhost

          ii = MAX0 ( ixx, 1 - NGhost)
          ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)

! active face is the one which has at least one active-zone on the either sides.
          IF ( ( imask_active_patch ( ii, iyy, izz) .EQ. 1) .OR.
     1         ( imask_active_patch ( ixxp1, iyy, izz) .EQ. 1) ) THEN

! --------------------------------------------------

      v1 ( 1:3) = vert_coords_patch ( ixx, iyy-1, izz-1, :)
      v2 ( 1:3) = vert_coords_patch ( ixx, iyy,   izz-1, :)
      v3 ( 1:3) = vert_coords_patch ( ixx, iyy-1, izz,   :)
      v4 ( 1:3) = vert_coords_patch ( ixx, iyy,   izz,   :)

      v1_mag = SQRT( v1 ( 1)**2 + v1 ( 2)**2 + v1 ( 3)**2)
      v2_mag = SQRT( v2 ( 1)**2 + v2 ( 2)**2 + v2 ( 3)**2)
      v3_mag = SQRT( v3 ( 1)**2 + v3 ( 2)**2 + v3 ( 3)**2)
      v4_mag = SQRT( v4 ( 1)**2 + v4 ( 2)**2 + v4 ( 3)**2)

      v1_dot_v2 = ( v1 ( 1) * v2 ( 1) + v1 ( 2) * v2 ( 2)
     1            + v1 ( 3) * v2 ( 3))

      v3_dot_v4 = ( v3 ( 1) * v4 ( 1) + v3 ( 2) * v4 ( 2)
     1            + v3 ( 3) * v4 ( 3))

! ------------------------------
! edges and vertices numbering for the x-face, edges are oriented as an
! anti-clockwise circulation.
!
!           e3
!   v3 -----------------  v4
!      \               /
!    e2 \             /  e1         z ^
!        \           /                |
!         \         /                 |
!     v1   --------- v2               --------> y
!              e4

! ----------
! edge 1 located in between v2 and v4.

      e1_cen ( :) = ( v2 ( :) + v4 ( :)) * 0.5

      e1_dl ( :)   = v4 ( :) - v2 ( :)     ! v2 --> v4

! ----------
! edge 2 located in between v3 and v1.

      e2_cen ( :) = ( v1 ( :) + v3 ( :)) * 0.5

      e2_dl ( :)   = v1 ( :) - v3 ( :)     ! v3 --> v1

! ----------
! edge 3 located in between v4 and v3.

      e3_cen ( :) = ( v3 ( :) + v4 ( :)) * 0.5
      e3_cen ( :) = v3_mag * e3_cen ( :) / SQRT ( e3_cen ( 1)**2
     1                         + e3_cen ( 2)**2 + e3_cen ( 3)**2)

      uvec ( :) = v3 ( :) - v4 ( :)      ! v4 --> v3

      uvec ( :) = uvec ( :) / SQRT ( uvec ( 1)**2  + uvec ( 2)**2
     1                                             + uvec ( 3)**2)

      angle = ACOS ( v3_dot_v4 / ( v3_mag * v4_mag))

      e3_dl ( :) = v3_mag * angle * uvec( :)  ! arc length = r*theta * unitvector

! ----------
! edge 4 located in between v1 and v2.

      e4_cen ( :) = ( v1 ( :) + v2 ( :)) * 0.5
      e4_cen ( :) = v1_mag * e4_cen ( :) / SQRT ( e4_cen ( 1)**2
     1                         + e4_cen ( 2)**2 + e4_cen ( 3)**2)

      uvec ( :) = v2 ( :) - v1 ( :)      ! v1 --> v2

      uvec ( :) = uvec ( :) / SQRT ( uvec ( 1)**2  + uvec ( 2)**2
     1                                             + uvec ( 3)**2)

      angle = ACOS ( v1_dot_v2 / ( v1_mag * v2_mag))

      e4_dl ( :) = v1_mag * angle * uvec( :)  ! arc length = r*theta * unitvector


! ------------------------------
! Now for the edges 1 - 4, the edge centers "e?_cen" and their
! vector length's "e?_dl" are available. Use these to obtain the
! pointwise vector potential "A_pt ( 3)" at each edge center
! and perform a 2nd order accurate line integration by doing the
! sum of "A_pt .dot. e?_dl".

      sum_A_dl = 0.0

! ----------
! A.dl along edge - 1.

      x = e1_cen ( 1) ;  y = e1_cen ( 2) ;  z = e1_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e1_dl ( 1)
     1                      + A_pt ( 2) * e1_dl ( 2)
     1                      + A_pt ( 3) * e1_dl ( 3))

! ----------
! A.dl along edge - 2.

      x = e2_cen ( 1) ;  y = e2_cen ( 2) ;  z = e2_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e2_dl ( 1)
     1                      + A_pt ( 2) * e2_dl ( 2)
     1                      + A_pt ( 3) * e2_dl ( 3))

! ----------
! A.dl along edge - 3.

      x = e3_cen ( 1) ;  y = e3_cen ( 2) ;  z = e3_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e3_dl ( 1)
     1                      + A_pt ( 2) * e3_dl ( 2)
     1                      + A_pt ( 3) * e3_dl ( 3))

! ----------
! A.dl along edge - 4.

      x = e4_cen ( 1) ;  y = e4_cen ( 2) ;  z = e4_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e4_dl ( 1)
     1                      + A_pt ( 2) * e4_dl ( 2)
     1                      + A_pt ( 3) * e4_dl ( 3))

! ------------------------------

      Bxf ( ixx, iyy, izz, 1)= sum_A_dl / xf_area_patch ( ixx, iyy, izz)

! --------------------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ii, iyy, izz) .EQ. 1) .OR. ...".

          END DO  ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! STEP 3.b) Initializing the y-face magnetic field.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxp1, iyyp1, izzp1, ii,     &
!$OMP&     x, y, z, U_pt, A_pt,                                         &

!$OMP&     v1, v2, v3, v4, v1_mag,  v2_mag,  v3_mag,  v4_mag,           &
!$OMP&     v1_dot_v2, v3_dot_v4, v1_dot_v3, v2_dot_v4, angle, sum_A_dl, &
!$OMP&     e1_cen, e2_cen, e3_cen, e4_cen,                              &
!$OMP&     e1_dl,  e2_dl,  e3_dl,  e4_dl, uvec,                         &
!$OMP&    rad, costheta, sintheta, cosphi, sinphi, R,                   &
!$OMP&    rho, prs, v0, vr, vtheta, vphi, vx, vy, vz, bfx, bfy, bfz,    &
!$OMP&    rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0,               &
!$OMP&    rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,               &
!$OMP&    vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,               &
!$OMP&    radius, distance, bmag, m1, m2, m3, phi, theta)

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,        &
!$OMP&    pi, gamma,                                                    &
!$OMP&    time, iproblem, rmin, rmax,                                   &
!$OMP&    zone_size, xcen_sector, ycen_sector, zcen_sector,             &
!$OMP&    xvert1_sector, yvert1_sector, zvert1_sector,                  &
!$OMP&    vert_coords_patch, imask_active_patch,                        &
!$OMP&    Byf, yf_area_patch)                                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    vert_coords_patch, imask_active_patch,                        &
!$acc&    Byf, yf_area_patch)                                           &
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxp1, iyyp1, izzp1, ii,             &
!$acc&     x, y, z, U_pt, A_pt,                                         &

!$acc&     v1, v2, v3, v4, v1_mag,  v2_mag,  v3_mag,  v4_mag,           &
!$acc&     v1_dot_v2, v3_dot_v4, v1_dot_v3, v2_dot_v4, angle, sum_A_dl, &
!$acc&     e1_cen, e2_cen, e3_cen, e4_cen,                              &
!$acc&     e1_dl,  e2_dl,  e3_dl,  e4_dl, uvec,                         &
!$acc&    rad, costheta, sintheta, cosphi, sinphi, R,                   &
!$acc&    rho, prs, v0, vr, vtheta, vphi, vx, vy, vz, bfx, bfy, bfz,    &
!$acc&    rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0,               &
!$acc&    rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,               &
!$acc&    vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,               &
!$acc&    radius, distance, bmag, m1, m2, m3, phi, theta)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 0 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          ii = MAX0 ( iyy, 1 - NGhost)
          iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)

! active face is the one which has at least one active-zone on the either sides.
          IF ( ( imask_active_patch ( ixx, ii, izz) .EQ. 1) .OR.
     1         ( imask_active_patch ( ixx, iyyp1, izz) .EQ. 1) ) THEN

! --------------------------------------------------

      v1 ( 1:3) = vert_coords_patch ( ixx-1, iyy, izz-1, :)
      v2 ( 1:3) = vert_coords_patch ( ixx-1, iyy, izz,   :)
      v3 ( 1:3) = vert_coords_patch ( ixx,   iyy, izz-1, :)
      v4 ( 1:3) = vert_coords_patch ( ixx,   iyy, izz,   :)

      v1_mag = SQRT( v1 ( 1)**2 + v1 ( 2)**2 + v1 ( 3)**2)
      v2_mag = SQRT( v2 ( 1)**2 + v2 ( 2)**2 + v2 ( 3)**2)
      v3_mag = SQRT( v3 ( 1)**2 + v3 ( 2)**2 + v3 ( 3)**2)
      v4_mag = SQRT( v4 ( 1)**2 + v4 ( 2)**2 + v4 ( 3)**2)

      v1_dot_v3 = ( v1 ( 1) * v3 ( 1) + v1 ( 2) * v3 ( 2)
     1            + v1 ( 3) * v3 ( 3))

      v2_dot_v4 = ( v2 ( 1) * v4 ( 1) + v2 ( 2) * v4 ( 2)
     1            + v2 ( 3) * v4 ( 3))

! ------------------------------
! edges and vertices numbering for the y-face, edges are oriented as an
! anti-clockwise circulation.
!
!               e3  ____ v4
!       v3 ___------    |           x ^
!         |             |             |
!      e2 |             | e1          |
!         |             |             |
!         |___          |             ---------> z
!       v1    ------____|
!               e4       v2
!

! ----------
! edge 1 located in between v2 and v4.

      e1_cen ( :) = ( v4 ( :) + v2 ( :)) * 0.5
      e1_cen ( :) = v2_mag * e1_cen ( :) / SQRT ( e1_cen ( 1)**2
     1                         + e1_cen ( 2)**2 + e1_cen ( 3)**2)

      uvec ( :) = v4 ( :) - v2 ( :)      ! v2 --> v4

      uvec ( :) = uvec ( :) / SQRT ( uvec ( 1)**2  + uvec ( 2)**2
     1                                             + uvec ( 3)**2)

      angle = ACOS ( v2_dot_v4 / ( v4_mag * v2_mag))

      e1_dl ( :) = v2_mag * angle * uvec( :)  ! arc length = r*theta * unitvector

! ----------
! edge 2 located in between v3 and v1.

      e2_cen ( :) = ( v3 ( :) + v1 ( :)) * 0.5
      e2_cen ( :) = v1_mag * e2_cen ( :) / SQRT ( e2_cen ( 1)**2
     1                         + e2_cen ( 2)**2 + e2_cen ( 3)**2)

      uvec ( :) = v1 ( :) - v3 ( :)      ! v3 --> v1

      uvec ( :) = uvec ( :) / SQRT ( uvec ( 1)**2  + uvec ( 2)**2
     1                                             + uvec ( 3)**2)

      angle = ACOS ( v1_dot_v3 / ( v3_mag * v1_mag))

      e2_dl ( :) = v1_mag * angle * uvec( :)  ! arc length = r*theta * unitvector

! ----------
! edge 3 located in between v4 and v3.

      e3_cen ( :) = ( v4 ( :) + v3 ( :)) * 0.5

      e3_dl ( :)   = v3 ( :) - v4 ( :)     ! v4 --> v3

! ----------
! edge 4 located in between v1 and v2.

      e4_cen ( :) = ( v1 ( :) + v2 ( :)) * 0.5

      e4_dl ( :)   = v2 ( :) - v1 ( :)     ! v1 --> v2


! ------------------------------
! Now for the edges 1 - 4, the edge centers "e?_cen" and their
! vector length's "e?_dl" are available. Use these to obtain the
! pointwise vector potential "A_pt ( 3)" at each edge center
! and perform a 2nd order accurate line integration by doing the
! sum of "A_pt .dot. e?_dl".

      sum_A_dl = 0.0

! ----------
! A.dl along edge - 1.

      x = e1_cen ( 1) ;  y = e1_cen ( 2) ;  z = e1_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e1_dl ( 1)
     1                      + A_pt ( 2) * e1_dl ( 2)
     1                      + A_pt ( 3) * e1_dl ( 3))

! ----------
! A.dl along edge - 2.

      x = e2_cen ( 1) ;  y = e2_cen ( 2) ;  z = e2_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e2_dl ( 1)
     1                      + A_pt ( 2) * e2_dl ( 2)
     1                      + A_pt ( 3) * e2_dl ( 3))

! ----------
! A.dl along edge - 3.

      x = e3_cen ( 1) ;  y = e3_cen ( 2) ;  z = e3_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e3_dl ( 1)
     1                      + A_pt ( 2) * e3_dl ( 2)
     1                      + A_pt ( 3) * e3_dl ( 3))

! ----------
! A.dl along edge - 4.

      x = e4_cen ( 1) ;  y = e4_cen ( 2) ;  z = e4_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e4_dl ( 1)
     1                      + A_pt ( 2) * e4_dl ( 2)
     1                      + A_pt ( 3) * e4_dl ( 3))

! ------------------------------

      Byf ( ixx, iyy, izz, 1)= sum_A_dl / yf_area_patch ( ixx, iyy, izz)

! --------------------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ii, iyy, izz) .EQ. 1) .OR. ...".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! STEP 3.c) Initializing the z-face magnetic field.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxp1, iyyp1, izzp1, ii,     &
!$OMP&     x, y, z, U_pt, A_pt,                                         &

!$OMP&     v1, v2, v3, v4, v1_mag,  v2_mag,  v3_mag,  v4_mag,           &
!$OMP&     v1_dot_v2, v3_dot_v4, v1_dot_v3, v2_dot_v4, angle, sum_A_dl, &
!$OMP&     e1_cen, e2_cen, e3_cen, e4_cen,                              &
!$OMP&     e1_dl,  e2_dl,  e3_dl,  e4_dl, uvec,                         &
!$OMP&    rad, costheta, sintheta, cosphi, sinphi, R,                   &
!$OMP&    rho, prs, v0, vr, vtheta, vphi, vx, vy, vz, bfx, bfy, bfz,    &
!$OMP&    rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0,               &
!$OMP&    rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,               &
!$OMP&    vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,               &
!$OMP&    radius, distance, bmag, m1, m2, m3, phi, theta)

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,        &
!$OMP&    pi, gamma,                                                    &
!$OMP&    time, iproblem, rmin, rmax,                                   &
!$OMP&    zone_size, xcen_sector, ycen_sector, zcen_sector,             &
!$OMP&    xvert1_sector, yvert1_sector, zvert1_sector,                  &
!$OMP&    vert_coords_patch, imask_active_patch,                        &
!$OMP&    Bzf, zf_area_patch)                                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    vert_coords_patch, imask_active_patch,                        &
!$acc&    Bzf, zf_area_patch)                                           &
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxp1, iyyp1, izzp1, ii,             &
!$acc&     x, y, z, U_pt, A_pt,                                         &

!$acc&     v1, v2, v3, v4, v1_mag,  v2_mag,  v3_mag,  v4_mag,           &
!$acc&     v1_dot_v2, v3_dot_v4, v1_dot_v3, v2_dot_v4, angle, sum_A_dl, &
!$acc&     e1_cen, e2_cen, e3_cen, e4_cen,                              &
!$acc&     e1_dl,  e2_dl,  e3_dl,  e4_dl, uvec,                         &
!$acc&    rad, costheta, sintheta, cosphi, sinphi, R,                   &
!$acc&    rho, prs, v0, vr, vtheta, vphi, vx, vy, vz, bfx, bfy, bfz,    &
!$acc&    rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0,               &
!$acc&    rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,               &
!$acc&    vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,               &
!$acc&    radius, distance, bmag, m1, m2, m3, phi, theta)

      DO izz = 0 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

! active face is the one which has at least one active-zone on the either sides.
! zones are masked along the "x" and "y" directions only and this is
! independent of the radial direction - "izz".

          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN

! --------------------------------------------------

      v1 ( 1:3) = vert_coords_patch ( ixx-1, iyy-1, izz, :)
      v2 ( 1:3) = vert_coords_patch ( ixx,   iyy-1, izz, :)
      v3 ( 1:3) = vert_coords_patch ( ixx-1, iyy,   izz, :)
      v4 ( 1:3) = vert_coords_patch ( ixx,   iyy,   izz, :)

      v1_mag = SQRT( v1 ( 1)**2 + v1 ( 2)**2 + v1 ( 3)**2)
      v2_mag = SQRT( v2 ( 1)**2 + v2 ( 2)**2 + v2 ( 3)**2)
      v3_mag = SQRT( v3 ( 1)**2 + v3 ( 2)**2 + v3 ( 3)**2)
      v4_mag = SQRT( v4 ( 1)**2 + v4 ( 2)**2 + v4 ( 3)**2)


      v1_dot_v2 = ( v1 ( 1) * v2 ( 1) + v1 ( 2) * v2 ( 2)
     1            + v1 ( 3) * v2 ( 3))

      v3_dot_v4 = ( v3 ( 1) * v4 ( 1) + v3 ( 2) * v4 ( 2)
     1            + v3 ( 3) * v4 ( 3))

      v1_dot_v3 = ( v1 ( 1) * v3 ( 1) + v1 ( 2) * v3 ( 2)
     1            + v1 ( 3) * v3 ( 3))

      v2_dot_v4 = ( v2 ( 1) * v4 ( 1) + v2 ( 2) * v4 ( 2)
     1            + v2 ( 3) * v4 ( 3))

! ------------------------------
! edges and vertices numbering for the z-face, edges are oriented as an
! anti-clockwise circulation.
!
!               e3
!       v3 -------------- v4
!         |              |           y ^
!         |              |             |
!      e2 |              | e1          |
!         |              |             |
!         |              |             ---------> x
!       v1 -------------- v2
!                 e4

! ----------
! edge 1 located in between v2 and v4.

      e1_cen ( :) = ( v4 ( :) + v2 ( :)) * 0.5
      e1_cen ( :) = v2_mag * e1_cen ( :) / SQRT ( e1_cen ( 1)**2
     1                         + e1_cen ( 2)**2 + e1_cen ( 3)**2)

      uvec ( :) = v4 ( :) - v2 ( :)      ! v2 --> v4

      uvec ( :) = uvec ( :) / SQRT ( uvec ( 1)**2  + uvec ( 2)**2
     1                                             + uvec ( 3)**2)

      angle = ACOS ( v2_dot_v4 / ( v4_mag * v2_mag))

      e1_dl ( :) = v2_mag * angle * uvec( :)  ! arc length = r*theta * unitvector

! ----------
! edge 2 located in between v3 and v1.

      e2_cen ( :) = ( v3 ( :) + v1 ( :)) * 0.5
      e2_cen ( :) = v1_mag * e2_cen ( :) / SQRT ( e2_cen ( 1)**2
     1                         + e2_cen ( 2)**2 + e2_cen ( 3)**2)

      uvec ( :) = v1 ( :) - v3 ( :)      ! v3 --> v1

      uvec ( :) = uvec ( :) / SQRT ( uvec ( 1)**2  + uvec ( 2)**2
     1                                             + uvec ( 3)**2)

      angle = ACOS ( v1_dot_v3 / ( v3_mag * v1_mag))

      e2_dl ( :) = v1_mag * angle * uvec( :)  ! arc length = r*theta * unitvector

! ----------
! edge 3 located in between v4 and v3.

      e3_cen ( :) = ( v3 ( :) + v4 ( :)) * 0.5
      e3_cen ( :) = v3_mag * e3_cen ( :) / SQRT ( e3_cen ( 1)**2
     1                         + e3_cen ( 2)**2 + e3_cen ( 3)**2)

      uvec ( :) = v3 ( :) - v4 ( :)      ! v4 --> v3

      uvec ( :) = uvec ( :) / SQRT ( uvec ( 1)**2  + uvec ( 2)**2
     1                                             + uvec ( 3)**2)

      angle = ACOS ( v3_dot_v4 / ( v3_mag * v4_mag))

      e3_dl ( :) = v3_mag * angle * uvec( :)  ! arc length = r*theta * unitvector

! ----------
! edge 4 located in between v1 and v2.

      e4_cen ( :) = ( v1 ( :) + v2 ( :)) * 0.5
      e4_cen ( :) = v1_mag * e4_cen ( :) / SQRT ( e4_cen ( 1)**2
     1                         + e4_cen ( 2)**2 + e4_cen ( 3)**2)

      uvec ( :) = v2 ( :) - v1 ( :)      ! v1 --> v2

      uvec ( :) = uvec ( :) / SQRT ( uvec ( 1)**2  + uvec ( 2)**2
     1                                             + uvec ( 3)**2)

      angle = ACOS ( v1_dot_v2 / ( v1_mag * v2_mag))

      e4_dl ( :) = v1_mag * angle * uvec( :)  ! arc length = r*theta * unitvector

! ------------------------------
! Now for the edges 1 - 4, the edge centers "e?_cen" and their
! vector length's "e?_dl" are available. Use these to obtain the
! pointwise vector potential "A_pt ( 3)" at each edge center
! and perform a 2nd order accurate line integration by doing the
! sum of "A_pt .dot. e?_dl".

      sum_A_dl = 0.0

! ----------
! A.dl along edge - 1.

      x = e1_cen ( 1) ;  y = e1_cen ( 2) ;  z = e1_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e1_dl ( 1)
     1                      + A_pt ( 2) * e1_dl ( 2)
     1                      + A_pt ( 3) * e1_dl ( 3))

! ----------
! A.dl along edge - 2.

      x = e2_cen ( 1) ;  y = e2_cen ( 2) ;  z = e2_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e2_dl ( 1)
     1                      + A_pt ( 2) * e2_dl ( 2)
     1                      + A_pt ( 3) * e2_dl ( 3))

! ----------
! A.dl along edge - 3.

      x = e3_cen ( 1) ;  y = e3_cen ( 2) ;  z = e3_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e3_dl ( 1)
     1                      + A_pt ( 2) * e3_dl ( 2)
     1                      + A_pt ( 3) * e3_dl ( 3))

! ----------
! A.dl along edge - 4.

      x = e4_cen ( 1) ;  y = e4_cen ( 2) ;  z = e4_cen ( 3)

       CALL PROBLEM_GENERATOR_U_A

     1     ( n_cc_components, x, y, z, time, iproblem, rmin, rmax,
     1       zone_size, xcen_sector, ycen_sector, zcen_sector,
     1       xvert1_sector, yvert1_sector, zvert1_sector,

     1       U_pt, A_pt )

      sum_A_dl = sum_A_dl + ( A_pt ( 1) * e4_dl ( 1)
     1                      + A_pt ( 2) * e4_dl ( 2)
     1                      + A_pt ( 3) * e4_dl ( 3))

! ------------------------------

      Bzf ( ixx, iyy, izz, 1)= sum_A_dl / zf_area_patch ( ixx, iyy, izz)

! --------------------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 0 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 4) Add the radial magnetic fields to the faces, directly.
! Since it is not possible to represent them using vector potential.

#if( 1 == 1)
! Initializing the z-face magnetic field, radial component.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxp1, iyyp1, izzp1, ii,     &
!$OMP&     x, y, z, bfx, bfy, bfz, bmag, v0, vx, vy, vz, R)             &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,        &
!$OMP&    pi, gamma,                                                    &
!$OMP&    time, iproblem, rmin, rmax,                                   &
!$OMP&    imask_active_patch, Bzf, zf_quadpts_patch, zf_quadwt_patch,   &
!$OMP&    zf_quadnorml_patch, zf_centroid_patch)                        &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    imask_active_patch, Bzf, zf_quadpts_patch, zf_quadwt_patch,   &
!$acc&    zf_quadnorml_patch, zf_centroid_patch)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxp1, iyyp1, izzp1, ii,             &
!$acc&     x, y, z, bfx, bfy, bfz, bmag, v0, vx, vy, vz, R)

      DO izz = 0 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

! active face is the one which has at least one active-zone on the either sides.
! zones are masked along the "x" and "y" directions only and this is
! independent of the radial direction - "izz".

          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN

! --------------------------------------------------

          bmag = 0.0

          R = SQRT ( zf_centroid_patch ( ixx, iyy, izz, 1)**2
     1             + zf_centroid_patch ( ixx, iyy, izz, 2)**2
     1             + zf_centroid_patch ( ixx, iyy, izz, 3)**2)


! -------------------
          IF ( iproblem == 1) THEN
! -------------------

            bmag = SQRT ( 4.0 * pi) / R**2

! -------------------
          ELSE IF ( iproblem == 10) THEN
! -------------------

            bmag = SIN ( x) * 0.75 * SQRT ( 4.0 * pi) * ( 4.5 / R)**2

! -------------------
          ELSE IF ( iproblem == 14) THEN
! -------------------

#if ( CHOOSE_APPL_HELIO == 1)
! By choosing the sign of "bmag" same or different, we can set up monopole or
! split monopole.

            IF ( zf_centroid_patch ( ixx, iyy, izz, 3) .GE. 0.0) THEN

              bmag = SQRT ( 2.0) * B0_sw * ( r_1**2 / R**2)

            ELSE

              bmag = SQRT ( 2.0) * B0_sw * ( r_1**2 / R**2)
#if (SPLIT_MONOPOLE_HELIO == 1)
              bmag = - bmag
#endif

            END IF
#endif
! End of "( CHOOSE_APPL_HELIO == 1)".

! -------------------
          END IF  ! End of "IF ( iproblem == 1) THEN".
! -------------------

          Bzf ( ixx, iyy, izz, 1) = Bzf ( ixx, iyy, izz, 1) + bmag

! --------------------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) ...".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 0 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 5) We now evaluate the local min and max of the MHD variables on this patch.
! Also keep track of the local divergence as it evolves in the course of a
! simulation.
! ----------------------------------------------------------------------

! These are used for obtaining global min and max fo all the variables which
! can be written out at each timestep. Bypassed for GPU.

      V_min_local ( :) = 1.0e30
      V_max_local ( :) = - 1.0e30

! These are used for obtaining global min and max fo all the variables which
! can be written out at each timestep. These variables are initialized
! on the CPU. Before reduction, OpenACC will automatically copy these
! initial values to the GPU.

      div_min = 1.0e30
      div_max = - 1.0e30

      flat_min = 1.0e30
      flat_max = - 1.0e30

      rhomin = 1.0e30
      rhomax = - 1.0e30

      vlxmin = 1.0e30
      vlxmax = - 1.0e30
      vlymin = 1.0e30
      vlymax = - 1.0e30
      vlzmin = 1.0e30
      vlzmax = - 1.0e30

      prsmin = 1.0e30
      prsmax = - 1.0e30

      bfxmin = 1.0e30
      bfxmax = - 1.0e30
      bfymin = 1.0e30
      bfymax = - 1.0e30
      bfzmin = 1.0e30
      bfzmax = - 1.0e30

! ------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    a11, a12, a13, a21)                                           &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,&
!$OMP&    rhofloor, prsfloor, ob_8pi,                                   &
!$OMP&    isospeed, gamma_euler, barotropic_stiffness,                  &
!$OMP&    imask_active_patch,                                           &
!$OMP&    xf_area_patch, yf_area_patch, zf_area_patch,                  &
!$OMP&    U, V, gamma_eos_gr, Bxf, Byf, Bzf)                            &

!$OMP&    REDUCTION ( MIN: V_min_local) REDUCTION ( MAX: V_max_local)   &
!$OMP&           REDUCTION ( MIN: div_min) REDUCTION ( MAX: div_max)    &
!$OMP&           REDUCTION ( MIN: flat_min) REDUCTION ( MAX: flat_max)  &
!$OMP&           REDUCTION ( MIN: rhomin) REDUCTION ( MAX: rhomax)      &
!$OMP&           REDUCTION ( MIN: vlxmin) REDUCTION ( MAX: vlxmax)      &
!$OMP&           REDUCTION ( MIN: vlymin) REDUCTION ( MAX: vlymax)      &
!$OMP&           REDUCTION ( MIN: vlzmin) REDUCTION ( MAX: vlzmax)      &
!$OMP&           REDUCTION ( MIN: prsmin) REDUCTION ( MAX: prsmax)      &
!$OMP&           REDUCTION ( MIN: bfxmin) REDUCTION ( MAX: bfxmax)      &
!$OMP&           REDUCTION ( MIN: bfymin) REDUCTION ( MAX: bfymax)      &
!$OMP&           REDUCTION ( MIN: bfzmin) REDUCTION ( MAX: bfzmax)      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    rhofloor, prsfloor, ob_8pi,                                   &
!$acc&    isospeed, gamma_euler, barotropic_stiffness,                  &
!$acc&    imask_active_patch,                                           &
!$acc&    xf_area_patch, yf_area_patch, zf_area_patch,                  &
!$acc&    U, V, gamma_eos_gr, Bxf, Byf, Bzf,                            &

!$acc&    V_min_local, V_max_local)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    a11, a12, a13, a21)                                           &

#if ( 0)
!$acc&    REDUCTION ( MIN: V_min_local) REDUCTION ( MAX: V_max_local)  ! ACC compilation error.
#endif
!$acc&           REDUCTION ( MIN: div_min) REDUCTION ( MAX: div_max)    &
!$acc&           REDUCTION ( MIN: flat_min) REDUCTION ( MAX: flat_max)  &
!$acc&           REDUCTION ( MIN: rhomin) REDUCTION ( MAX: rhomax)      &
!$acc&           REDUCTION ( MIN: vlxmin) REDUCTION ( MAX: vlxmax)      &
!$acc&           REDUCTION ( MIN: vlymin) REDUCTION ( MAX: vlymax)      &
!$acc&           REDUCTION ( MIN: vlzmin) REDUCTION ( MAX: vlzmax)      &
!$acc&           REDUCTION ( MIN: prsmin) REDUCTION ( MAX: prsmax)      &
!$acc&           REDUCTION ( MIN: bfxmin) REDUCTION ( MAX: bfxmax)      &
!$acc&           REDUCTION ( MIN: bfymin) REDUCTION ( MAX: bfymax)      &
!$acc&           REDUCTION ( MIN: bfzmin) REDUCTION ( MAX: bfzmax)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! ----------------------------------------

          V ( ixx, iyy, izz, 1) = AMAX1 ( U ( ixx, iyy, izz, 1, 1),
     1                                    rhofloor)
          U ( ixx, iyy, izz, 1, 1) = V ( ixx, iyy, izz, 1)

! Extract velocity vector field.

          a11 = V ( ixx, iyy, izz, 1) / ( V ( ixx, iyy, izz, 1)**2
     1                                  + rhofloor**2)

          V ( ixx, iyy, izz, 2) = a11 * U ( ixx, iyy, izz, 2, 1)
          V ( ixx, iyy, izz, 3) = a11 * U ( ixx, iyy, izz, 3, 1)
          V ( ixx, iyy, izz, 4) = a11 * U ( ixx, iyy, izz, 4, 1)

! Extract magnetic vector field.

#if ( 1 == 0)
! This might be more stable.

          V ( ixx, iyy, izz, 6) = B1_zc ( ixx, iyy, izz, 1)
          V ( ixx, iyy, izz, 7) = B2_zc ( ixx, iyy, izz, 1)
          V ( ixx, iyy, izz, 8) = B3_zc ( ixx, iyy, izz, 1)

#else
! We do update the zone centered magnetic fields with high accuracy. They
! may not be div-free. However, they can be used to extract a highly accurate
! pressure.

          V ( ixx, iyy, izz, 6) = U ( ixx, iyy, izz, 6, 1)
          V ( ixx, iyy, izz, 7) = U ( ixx, iyy, izz, 7, 1)
          V ( ixx, iyy, izz, 8) = U ( ixx, iyy, izz, 8, 1)

#endif

! This expression for the zone-averaged kinetic energy remains the same in this loop.

          a11 = 0.5 * V ( ixx, iyy, izz, 1)
     1        * ( V ( ixx, iyy, izz, 2)**2 + V ( ixx, iyy, izz, 3)**2
     1          + V ( ixx, iyy, izz, 4)**2 )

! Even the lowest order magnetic energy varies quadratically within the zone. So,
! this expression for the magnetic energy changes at different locations in the zone.

          a12 = ( V ( ixx, iyy, izz, 6)**2
     1          + V ( ixx, iyy, izz, 7)**2
     1          + V ( ixx, iyy, izz, 8)**2) * ob_8pi

          a13 = ( U ( ixx, iyy, izz, 5, 1) - a11 - a12 )
     1               * ( gamma_eos_gr ( ixx, iyy, izz) - 1.0)

#if ( CHOOSE_ADIABATIC == 1)
! Using adiabatic EOS.
          V ( ixx, iyy, izz, 5) = AMAX1 ( prsfloor, a13 )
#else
! Using isothermal/barotropic/polytropic EOS.
          V ( ixx, iyy, izz, 5) = isospeed**2 * V ( ixx, iyy, izz, 1)
     1    + barotropic_stiffness * V ( ixx, iyy, izz, 1)**gamma_euler
#endif


! ------------------------------
#if ( NFLUID >= 1)
          V ( ixx, iyy, izz, NVAR + 1: ) =
     1    U ( ixx, iyy, izz, NVAR + 1: , 1) / V ( ixx, iyy, izz, 1)
#endif


! ------------------------------
! Take min, max of all the fluid variables.

        rhomin = AMIN1 ( rhomin, V ( ixx, iyy, izz, 1) )
        vlxmin = AMIN1 ( vlxmin, V ( ixx, iyy, izz, 2) )
        vlymin = AMIN1 ( vlymin, V ( ixx, iyy, izz, 3) )
        vlzmin = AMIN1 ( vlzmin, V ( ixx, iyy, izz, 4) )
        prsmin = AMIN1 ( prsmin, V ( ixx, iyy, izz, 5) )

        bfxmin = AMIN1 ( bfxmin, V ( ixx, iyy, izz, 6) )
        bfymin = AMIN1 ( bfymin, V ( ixx, iyy, izz, 7) )
        bfzmin = AMIN1 ( bfzmin, V ( ixx, iyy, izz, 8) )

        rhomax = AMAX1 ( rhomax, V ( ixx, iyy, izz, 1) )
        vlxmax = AMAX1 ( vlxmax, V ( ixx, iyy, izz, 2) )
        vlymax = AMAX1 ( vlymax, V ( ixx, iyy, izz, 3) )
        vlzmax = AMAX1 ( vlzmax, V ( ixx, iyy, izz, 4) )
        prsmax = AMAX1 ( prsmax, V ( ixx, iyy, izz, 5) )

        bfxmax = AMAX1 ( bfxmax, V ( ixx, iyy, izz, 6) )
        bfymax = AMAX1 ( bfymax, V ( ixx, iyy, izz, 7) )
        bfzmax = AMAX1 ( bfzmax, V ( ixx, iyy, izz, 8) )

! ------------------------------
! Now evaluate the divergence here and take its min, max. Only do this
! for the interior zones, because divergence is not preserved in
! non-physical zones.

        IF ( ( 1 .LE. ixx) .AND. ( ixx .LE. ntp_zones) .AND.
     1       ( 1 .LE. iyy) .AND. ( iyy .LE. ntp_zones) .AND.
     1       ( 1 .LE. izz) .AND. ( izz .LE. nr_zones) ) THEN

        a21 = xf_area_patch ( ixx, iyy, izz)
     1             * Bxf ( ixx, iyy, izz, 1)
     1             + yf_area_patch ( ixx, iyy, izz)
     1             * Byf ( ixx, iyy, izz, 1)
     1             + zf_area_patch ( ixx, iyy, izz)
     1             * Bzf ( ixx, iyy, izz, 1)

     1             - xf_area_patch ( ixxm1, iyy, izz)
     1             * Bxf ( ixxm1, iyy, izz, 1)
     1             - yf_area_patch ( ixx, iyym1, izz)
     1             * Byf ( ixx, iyym1, izz, 1)
     1             - zf_area_patch ( ixx, iyy, izzm1)
     1             * Bzf ( ixx, iyy, izzm1, 1)

        div_min = AMIN1 ( div_min, a21 )
        div_max = AMAX1 ( div_max, a21 )

        END IF  ! End of divergence evaluation.

! ----------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! The REDUCTION clause in OpenACC automatically reduces and puts
! the values on the CPU. Therefore, the subsequent initialization will
! be done on the CPU.

      V_min_local ( 1) = rhomin
      V_min_local ( 2) = vlxmin
      V_min_local ( 3) = vlymin
      V_min_local ( 4) = vlzmin
      V_min_local ( 5) = prsmin

      V_min_local ( 6) = bfxmin
      V_min_local ( 7) = bfymin
      V_min_local ( 8) = bfzmin

      V_min_local ( 9) = div_min

      V_max_local ( 1) = rhomax
      V_max_local ( 2) = vlxmax
      V_max_local ( 3) = vlymax
      V_max_local ( 4) = vlzmax
      V_max_local ( 5) = prsmax

      V_max_local ( 6) = bfxmax
      V_max_local ( 7) = bfymax
      V_max_local ( 8) = bfzmax

      V_max_local ( 9) = div_max

! ----------------------------------------------------------------------
! STEP 5) Evaluate "dtcur" for this simulation. If we ever do the ensuing
! loop nest as a REDUCTION in OpenACC, then OpenACC will take the time
! from CPU, reduce it and send the reduced value back to GPU.

      dtcur = 1.0e32 ! Initialized on CPU

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, tempab, tempac,              &
!$OMP&    rho1, vlx1, vly1, vlz1, prs1, bfx1, bfy1, bfz1 )              &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    pi, gamma, cfl_coef, V, v_scale_length,                       &
!$OMP&    imask_active_patch )                                          &

!$OMP&            REDUCTION ( MIN: dtcur)                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel if ( do_on_gpu) vector_length(64) present (              &
!$acc&    pi, gamma, cfl_coef, V, v_scale_length,                       &
!$acc&    imask_active_patch )                                          &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, tempab, tempac,                      &
!$acc&    rho1, vlx1, vly1, vlz1, prs1, bfx1, bfy1, bfz1 )              &

!$acc&    REDUCTION ( MIN: dtcur )

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

            rho1 = V ( ixx, iyy, izz, 1)
            vlx1 = V ( ixx, iyy, izz, 2)
            vly1 = V ( ixx, iyy, izz, 3)
            vlz1 = V ( ixx, iyy, izz, 4)
            prs1 = V ( ixx, iyy, izz, 5)
            bfx1 = V ( ixx, iyy, izz, 6)
            bfy1 = V ( ixx, iyy, izz, 7)
            bfz1 = V ( ixx, iyy, izz, 8)

#if( USE_B0B1_FORMULATION == 1)
            bfx1 = bfx1 + B0field ( ixx, iyy, izz, 1)
            bfy1 = bfy1 + B0field ( ixx, iyy, izz, 2)
            bfz1 = bfz1 + B0field ( ixx, iyy, izz, 3)
#endif

            tempab = SQRT ( vlx1**2 + vly1**2 + vlz1**2)
            tempac = gamma * prs1
     1             + ( bfx1**2 + bfy1**2 + bfz1**2) / ( 4.0 * pi)
            tempac = SQRT ( tempac / rho1) + tempab

            tempab = cfl_coef * v_scale_length ( ixx, iyy, izz) / tempac

            dtcur = AMIN1 ( dtcur, tempab)

          END IF ! End of IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.
!$acc end parallel
!$OMP END PARALLEL DO

! Because "dtcur" is a REDUCTION variable in OpenACC, it will be
! automatically copied over to the CPU. No need to do anything further.

!$acc end data


      END SUBROUTINE INIT_PATCH_CUBE_SPHERE


