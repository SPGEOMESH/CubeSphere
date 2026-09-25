



      SUBROUTINE MAKE_EDGE_EFIELDS_P1

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1         ii_mode, jj_mode, kk_mode,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch,

     1           v_scale_length,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc, Ex_zc, Ey_zc, Ez_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr,
     1         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,

     1         S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr,
     1         Efx, Efy, Efz)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in all relevant geometric information, as described in the
! "LIMIT_ALL_SLOPES_P?" routines.

! Takes in the reconstructed fluid and magnetic field variables in
! "U, entropygr, B1_zc, B2_zc, B3_zc".

! Returns electric fields "Efx, Efy, Efz".
! Electric fields are already zeroed in "MAKE_FACIAL_FLUXES"
! They are set in "MAKE_FACIAL_FLUXES" if "USE_FACIAL_EFIELDS == 1";
! in which case we return.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes,
     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

! ------------------------------
! INDICES OF MOMENTS THAT HELP IN THE RECONSTRUCTION (set in BUILD_MODE_INDICES):-

      INTEGER, DIMENSION ( n_cc_modes) ::
     1         ii_mode, jj_mode, kk_mode

! ------------------------------
! VARIABLES THAT GENERATE THE MESH (set in EXTRACT_VERT_COORDS_PATCH):-

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

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

! ------------------------------
! Integration points for magnetic field in the multid riemann solver.

      REAL xb_quadpts_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

      REAL xb_quadwt_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

      REAL xb_quadnorml_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ----------

      REAL yb_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

      REAL yb_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

      REAL yb_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ----------

      REAL zb_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

      REAL zb_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

      REAL zb_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

! ------------------------------

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ------------------------------
! VOLUMETRIC VARIABLES (set in COORDS_TO_WT_QUADPTS_VOLUME_ISO2):-

      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------
! MOMENTS FOR BASIS FUNCTIONS (set in INIT_BASIS_FUNCTION_COEFS):-

      REAL Coef_Basis_Moments ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes)

! ------------------------------
! VARIABLES THAT FACILITATE STENCIL OPERATIONS (set in STENCIL_INDEXING_TRICORNERS):-

      INTEGER istencil_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                         1 - NGhost: ntp_zones + NGhost, 2, 12)

! ------------------------------
! FLUID VARIABLES:-

      REAL U ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL entropygr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   4)

      REAL V ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           n_fc_modes) ::
     1   Bxf, Byf, Bzf

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes) ::
     1   B1_zc, B2_zc, B3_zc

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost, n_cc_modes) ::
     1   Ex_zc, Ey_zc, Ez_zc

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost) ::
     1      gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1      divvelgr, msonicgr, flattengr

! ------------------------------
! Variables used for div-free reconstruction of magnetic fields:-

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           NBFLD_QUAD_PTS) ::
     1   Bxf_quadpts, Byf_quadpts, Bzf_quadpts

! ------------------------------
! FLUX AND ELECTRIC FIELD VARIABLES:-

! Store extremal speeds from the 1D Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost) ::
     1   S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr


! Electric fields at the edges of the mesh.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      Efx, Efy, Efz

! ----------------------------------------------------------------------

! Variables past this point are internal to the subroutine.

! LD, RD, RU, LU should be in a counter-clockwise sequence relative to
! our choice of local z-direction.

      INTEGER ixx, iyy, izz, ixxLD, iyyLD, izzLD, ixxRD, iyyRD, izzRD,
     1        ixxRU, iyyRU, izzRU, ixxLU, iyyLU, izzLU,
     1        bflag,
     1        ixx1, iyy1, ixx2, iyy2, ixx3, iyy3, j,
     1        itri_corner

      REAL xLD, yLD, zLD, xRD, yRD, zRD, xRU, yRU, zRU, xLU, yLU, zLU,
     1     oblLD, oblRD, oblRU, oblLU, xec, yec, zec, tempaa, wt

! ----------------------------------------------------------------------

! These variables go through the 2D Riemann solver interface. They can be declared
! SHARED and then the "TwoD_RS_Ez_4State, TwoD_RS_Ez_3State" set of multi-dimensional
! Riemann solvers can be used.

      REAL n1 ( 3), n2 ( 3), n3 ( 3), B0 ( 3), tau ( 3), neta ( 3),
     1     time_use, Bx0, By0, Bz0,

     1     tau_1 ( 3), tau_2 ( 3), tau_3 ( 3),
     1     neta_1 ( 3), neta_2 ( 3), neta_3( 3)

      REAL S_R, S_L, S_U, S_D, cos_t, sin_t,
     1     Bneta_U, Bneta_D, By_R, By_L,
     1     BSSx, BSSy, ESSz,

     1     Ex_LD, Ey_LD, Ez_LD, Ex_RD, Ey_RD, Ez_RD,
     1     Ex_LU, Ey_LU, Ez_LU, Ex_RU, Ey_RU, Ez_RU,

     1     obl_1, obl_2, obl_3,
     1     x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3,

     1     Ex_1, Ey_1, Ez_1, Ex_2, Ey_2, Ez_2, Ex_3, Ey_3, Ez_3,

     1     S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,
     1     Bneta_1, Bneta_2, Bneta_3

      REAL ob3

! ----------------------------------------------------------------------

!$acc routine( Two_D_RS_Ez_4State_HLL ) seq
!$acc routine( Two_D_RS_Ez_3State_HLL ) seq
!$acc routine( GET_MAGNETIC_FIELD ) seq

! ----------------------------------------------------------------------
! Zero out "Efx, Efy, Efz" that will be used as accumulation variables in this
! subroutine.

!$acc kernels
      Efx ( :, :, :) = 0.0
      Efy ( :, :, :) = 0.0
      Efz ( :, :, :) = 0.0
!$acc end kernels

! ------------------------------

      ob3 = 1.0 / 3.0

      IF ( USE_FACIAL_EFIELDS .EQ. 1) RETURN

! ----------------------------------------------------------------------
! If we get here, zero out all the electric fields and obtain them via
! the 2D Riemann solver.

#if ( USE_FACIAL_EFIELDS != 1)
! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    Efx, Efy, Efz)                                                &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present( Efx, Efy, Efz)
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          Efx ( ixx, iyy, izz) = 0.0
          Efy ( ixx, iyy, izz) = 0.0
          Efz ( ixx, iyy, izz) = 0.0

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
#endif

! End of zeroing out all the electric fields because we will obtain them via
! the 2D Riemann solver.
! ----------------------------------------------------------------------
! Move above initialization to GPU.

!$acc data copyin ( ob3 )


! ----------------------------------------------------------------------
! STEP 1) Build "Efx" in the positive x-direction. This electric field is
! zone-centered in the x-direction and edge-centered in the y- and z-directions.
! The global positive x-direction becomes the local z-direction for the
! 2D Riemann solver.

! At third and fourth order, loop over 3 quadrature points in the edge of interest.
      DO j = 1, 3

#if (NUMORDER == 2)
! At second order, we need only one nodal point.
      IF ( j .GT. 1) EXIT
#endif

      wt = 1.0
#if (NUMORDER >= 3)
      wt = 1.0 / 6.0  ! These are the weights for simpson integration.
      IF ( j .EQ. 1) wt = 4.0 / 6.0
#endif

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxLD, iyyLD, izzLD,         &
!$OMP&    ixxRD, iyyRD, izzRD, ixxRU, iyyRU, izzRU, ixxLU, iyyLU, izzLU,&

!$OMP&     xLD, yLD, zLD, xRD, yRD, zRD, xRU, yRU, zRU, xLU, yLU, zLU,  &
!$OMP&     oblLD, oblRD, oblRU, oblLU, xec, yec, zec, tempaa,           &
!$OMP&     bflag, time_use, Bx0, By0, Bz0,                              &
!$OMP&     Ex_LD, Ey_LD, Ez_LD, Ex_RD, Ey_RD, Ez_RD,                    &
!$OMP&     Ex_LU, Ey_LU, Ez_LU, Ex_RU, Ey_RU, Ez_RU,                    &

!$OMP&     n1, n2, n3, B0, S_R, S_L, S_U, S_D, cos_t, sin_t,            &
!$OMP&     Bneta_U, Bneta_D, By_R, By_L,                                &
!$OMP&     BSSx, BSSy, ESSz)                                            &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum, j, wt,&
!$OMP&    vert_coords_patch, xf_quadnorml_patch, v_scale_length,        &
!$OMP&    v_centroid_patch, Coef_Basis_Moments,                         &

!$OMP&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$OMP&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$OMP&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$OMP&    Efx )                                                         &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 3)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum, j, wt,                                              &
!$acc&    vert_coords_patch, yf_quadnorml_patch, v_scale_length,        &
!$acc&    v_centroid_patch, Coef_Basis_Moments,                         &

!$acc&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$acc&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$acc&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$acc&    Efx )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxLD, iyyLD, izzLD,                 &
!$acc&    ixxRD, iyyRD, izzRD, ixxRU, iyyRU, izzRU, ixxLU, iyyLU, izzLU,&

!$acc&     xLD, yLD, zLD, xRD, yRD, zRD, xRU, yRU, zRU, xLU, yLU, zLU,  &
!$acc&     oblLD, oblRD, oblRU, oblLU, xec, yec, zec, tempaa,           &
!$acc&     bflag, time_use, Bx0, By0, Bz0,                              &
!$acc&     Ex_LD, Ey_LD, Ez_LD, Ex_RD, Ey_RD, Ez_RD,                    &
!$acc&     Ex_LU, Ey_LU, Ez_LU, Ex_RU, Ey_RU, Ez_RU,                    &

!$acc&     n1, n2, n3, B0, S_R, S_L, S_U, S_D, cos_t, sin_t,            &
!$acc&     Bneta_U, Bneta_D, By_R, By_L,                                &
!$acc&     BSSx, BSSy, ESSz)

      DO izz = 0, nr_zones

        DO iyy = 0, ntp_zones

          DO ixx = 1, ntp_zones

! ------------------------------------------------------------

! Fill in "ixxLD, iyyLD, izzLD, ixxRD, iyyRD, izzRD,
! ixxRU, iyyRU, izzRU, ixxLU, iyyLU, izzLU" in counter-clockwise sense
! relative to the positive x-axis.

          ixxLD = ixx
          iyyLD = iyy
          izzLD = izz

          ixxRD = ixx
          iyyRD = iyy + 1
          izzRD = izz

          ixxRU = ixx
          iyyRU = iyy + 1
          izzRU = izz + 1

          ixxLU = ixx
          iyyLU = iyy
          izzLU = izz + 1

! --------------------------------------------------
! Find the edge-center "xec, yec, zec" by bisecting the two x-vertices.
! Then scale it up to be of the correct radius.
! These are needed in general if we are using the reconstruction.

          IF ( j .EQ. 1) THEN  ! First quadrature point is the center point.

            xec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 1)  ! Bisect the x-vertices.
     1                  + vert_coords_patch ( ixx - 1, iyy, izz, 1) )
            yec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 2)
     1                  + vert_coords_patch ( ixx - 1, iyy, izz, 2) )
            zec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 3)
     1                  + vert_coords_patch ( ixx - 1, iyy, izz, 3) )

            tempaa = SQRT ( vert_coords_patch ( ixx, iyy, izz, 1)**2
     1                    + vert_coords_patch ( ixx, iyy, izz, 2)**2
     1                    + vert_coords_patch ( ixx, iyy, izz, 3)**2)
     1             / SQRT ( xec**2 + yec**2 + zec**2)

            xec = xec * tempaa  ! Rescale to have correct radius.
            yec = yec * tempaa
            zec = zec * tempaa

          ELSE IF ( j .EQ. 2) THEN  ! Second quadrature point is one of the end points.

            xec = vert_coords_patch ( ixx, iyy, izz, 1)
            yec = vert_coords_patch ( ixx, iyy, izz, 2)
            zec = vert_coords_patch ( ixx, iyy, izz, 3)

          ELSE IF ( j .EQ. 3) THEN  ! Third quadrature point is the other end point.

            xec = vert_coords_patch ( ixx - 1, iyy, izz, 1)
            yec = vert_coords_patch ( ixx - 1, iyy, izz, 2)
            zec = vert_coords_patch ( ixx - 1, iyy, izz, 3)

          END IF

! ------------------------------
! Feed ( xec, yec, zec) to the tiny subroutine that evaluates "B0".
! Default is zero.

          B0 = 0.0
          Bx0 = 0.0
          By0 = 0.0
          Bz0 = 0.0

#if ( USE_B0B1_FORMULATION == 1)

          bflag = 1
          time_use = 0

          CALL GET_MAGNETIC_FIELD ( bflag, time_use,  xec, yec, zec,
     1                              Bx0, By0, Bz0)

          B0 ( 1) = Bx0
          B0 ( 2) = By0
          B0 ( 3) = Bz0

#endif

! ------------------------------
! Build "oblLD, oblRD, oblRU, oblLU" so that we can interpolate in the local
! coordinates of each of the four participating zones.
! These are needed in general if we are using the reconstruction.

          oblLD = 1.0 / AMAX1 ( v_scale_length ( ixxLD, iyyLD, izzLD),
     1                          smallnum)

          oblRD = 1.0 / AMAX1 ( v_scale_length ( ixxRD, iyyRD, izzRD),
     1                          smallnum)

          oblRU = 1.0 / AMAX1 ( v_scale_length ( ixxRU, iyyRU, izzRU),
     1                          smallnum)

          oblLU = 1.0 / AMAX1 ( v_scale_length ( ixxLU, iyyLU, izzLU),
     1                          smallnum)

! ------------------------------
! Build scaled "xLD, yLD, zLD, xRD, yRD, zRD, xRU, yRU, zRU, xLU, yLU, zLU".
! These are needed in general if we are using the reconstruction.

! ----------

          xLD = ( xec - v_centroid_patch ( ixxLD, iyyLD, izzLD, 1) )
     1        * oblLD
          yLD = ( yec - v_centroid_patch ( ixxLD, iyyLD, izzLD, 2) )
     1        * oblLD
          zLD = ( zec - v_centroid_patch ( ixxLD, iyyLD, izzLD, 3) )
     1        * oblLD

! ----------

          xRD = ( xec - v_centroid_patch ( ixxRD, iyyRD, izzRD, 1) )
     1        * oblRD
          yRD = ( yec - v_centroid_patch ( ixxRD, iyyRD, izzRD, 2) )
     1        * oblRD
          zRD = ( zec - v_centroid_patch ( ixxRD, iyyRD, izzRD, 3) )
     1        * oblRD

! ----------

          xRU = ( xec - v_centroid_patch ( ixxRU, iyyRU, izzRU, 1) )
     1        * oblRU
          yRU = ( yec - v_centroid_patch ( ixxRU, iyyRU, izzRU, 2) )
     1        * oblRU
          zRU = ( zec - v_centroid_patch ( ixxRU, iyyRU, izzRU, 3) )
     1        * oblRU

! ----------

          xLU = ( xec - v_centroid_patch ( ixxLU, iyyLU, izzLU, 1) )
     1        * oblLU
          yLU = ( yec - v_centroid_patch ( ixxLU, iyyLU, izzLU, 2) )
     1        * oblLU
          zLU = ( zec - v_centroid_patch ( ixxLU, iyyLU, izzLU, 3) )
     1        * oblLU

! ------------------------------
! Build local normals so that globally positive x-direction is indeed the
! locally positive z-direction. Store these in "n1, n2". "n3" not needed.
! "n2" is a vector in the radial direction; becomes the local y-direction.
! "n3" is a vector along the x-direction of the mesh; which becomes the local z-direction.
! "n1" is a vector along the y-direction of the mesh; which becomes the local x-direction.

          tempaa = 1.0 / SQRT ( xec**2 + yec**2 + zec**2)

! Here "n2" is in the global r-direction of the mesh, which is our local y-direction.
          n2 ( 1) = xec * tempaa   ! This is in the global r-direction.
          n2 ( 2) = yec * tempaa
          n2 ( 3) = zec * tempaa

! Here "n3" is in the x-direction of the mesh, which is our local z-direction.
! But as the quadrature points change "n2" will change slightly; so we make sure
! that "n3" is always a unit vector that is perpendicular to "n2".

          n3 ( 1) = vert_coords_patch ( ixx, iyy, izz, 1)
     1            - vert_coords_patch ( ixx - 1, iyy, izz, 1)
          n3 ( 2) = vert_coords_patch ( ixx, iyy, izz, 2)
     1            - vert_coords_patch ( ixx - 1, iyy, izz, 2)
          n3 ( 3) = vert_coords_patch ( ixx, iyy, izz, 3)
     1            - vert_coords_patch ( ixx - 1, iyy, izz, 3)

          tempaa = n3 ( 1) * n2 ( 1) + n3 ( 2) * n2 ( 2)
     1           + n3 ( 3) * n2 ( 3)

          n3 ( 1) = n3 ( 1) - tempaa * n2 ( 1)
          n3 ( 2) = n3 ( 2) - tempaa * n2 ( 2)
          n3 ( 3) = n3 ( 3) - tempaa * n2 ( 3)

          tempaa = 1.0 / SQRT ( n3 ( 1)**2 + n3 ( 2)**2 + n3 ( 3)**2)

          n3 ( :) = n3 ( :) * tempaa

! "n1 = n2 X n3". Here "n1" is in the y-direction of the mesh, which is our local x-direction.
          n1 ( 1) = n2 ( 2) * n3 ( 3) - n2 ( 3) * n3 ( 2)
          n1 ( 2) = n2 ( 3) * n3 ( 1) - n2 ( 1) * n3 ( 3)
          n1 ( 3) = n2 ( 1) * n3 ( 2) - n2 ( 2) * n3 ( 1)

! When the local x-axis, "n3" is mapped to the z-direction, we always have
! an orthogonal intersection of the axes.
          cos_t = 0.0  ;  sin_t = 1.0

! ------------------------------
! Now that we have the properly scaled relative coordinates from each of the
! four particpating zones, we can use these to obtain the four input electric field components
! "Ex,y,z_RU, Ex,y,z_LU, Ex,y,z_LD, Ex,y,z_RD" that go into the 2D Riemann solver.

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ == 0)
! --------------------
! Here the building of the electric fields requires that they have been reconstructed.

          Ex_LD = Ex_zc ( ixxLD, iyyLD, izzLD, 1)
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 2) * xLD
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 3) * yLD
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 4) * zLD
#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 5)
     1          * (xLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 5))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 6)
     1          * (yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 6))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 7)
     1          * (zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 7))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 8)
     1          * (xLD * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 8))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 9)
     1          * (yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 9))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 10)
     1          * (xLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 11)
     1          * (xLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 11))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 12)
     1          * (yLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 12))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 13)
     1          * (zLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 13))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 14)
     1          * (xLD**2 * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 14))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 15)
     1          * (xLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 15))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 16)
     1          * (xLD * yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 16))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 17)
     1          * (yLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 17))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 18)
     1          * (xLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 18))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 19)
     1          * (yLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 19))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 20)
     1          * (xLD * yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 20))
#endif

          Ey_LD = Ey_zc ( ixxLD, iyyLD, izzLD, 1)
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 2) * xLD
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 3) * yLD
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 4) * zLD
#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 5)
     1          * (xLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 5))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 6)
     1          * (yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 6))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 7)
     1          * (zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 7))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 8)
     1          * (xLD * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 8))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 9)
     1          * (yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 9))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 10)
     1          * (xLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 11)
     1          * (xLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 11))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 12)
     1          * (yLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 12))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 13)
     1          * (zLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 13))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 14)
     1          * (xLD**2 * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 14))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 15)
     1          * (xLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 15))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 16)
     1          * (xLD * yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 16))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 17)
     1          * (yLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 17))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 18)
     1          * (xLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 18))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 19)
     1          * (yLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 19))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 20)
     1          * (xLD * yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 20))
#endif

          Ez_LD = Ez_zc ( ixxLD, iyyLD, izzLD, 1)
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 2) * xLD
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 3) * yLD
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 4) * zLD

#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 5)
     1          * (xLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 5))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 6)
     1          * (yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 6))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 7)
     1          * (zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 7))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 8)
     1          * (xLD * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 8))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 9)
     1          * (yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 9))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 10)
     1          * (xLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 11)
     1          * (xLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 11))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 12)
     1          * (yLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 12))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 13)
     1          * (zLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 13))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 14)
     1          * (xLD**2 * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 14))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 15)
     1          * (xLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 15))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 16)
     1          * (xLD * yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 16))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 17)
     1          * (yLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 17))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 18)
     1          * (xLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 18))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 19)
     1          * (yLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 19))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 20)
     1          * (xLD * yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 20))
#endif

! ----------

          Ex_RD = Ex_zc ( ixxRD, iyyRD, izzRD, 1)
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 2) * xRD
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 3) * yRD
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 4) * zRD

#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 5)
     1          * (xRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 5))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 6)
     1          * (yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 6))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 7)
     1          * (zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 7))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 8)
     1          * (xRD * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 8))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 9)
     1          * (yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 9))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 10)
     1          * (xRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 11)
     1          * (xRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 11))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 12)
     1          * (yRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 12))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 13)
     1          * (zRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 13))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 14)
     1          * (xRD**2 * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 14))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 15)
     1          * (xRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 15))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 16)
     1          * (xRD * yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 16))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 17)
     1          * (yRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 17))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 18)
     1          * (xRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 18))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 19)
     1          * (yRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 19))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 20)
     1          * (xRD * yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 20))
#endif

          Ey_RD = Ey_zc ( ixxRD, iyyRD, izzRD, 1)
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 2) * xRD
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 3) * yRD
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 4) * zRD

#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 5)
     1          * (xRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 5))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 6)
     1          * (yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 6))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 7)
     1          * (zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 7))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 8)
     1          * (xRD * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 8))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 9)
     1          * (yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 9))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 10)
     1          * (xRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 11)
     1          * (xRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 11))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 12)
     1          * (yRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 12))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 13)
     1          * (zRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 13))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 14)
     1          * (xRD**2 * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 14))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 15)
     1          * (xRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 15))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 16)
     1          * (xRD * yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 16))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 17)
     1          * (yRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 17))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 18)
     1          * (xRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 18))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 19)
     1          * (yRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 19))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 20)
     1          * (xRD * yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 20))
#endif

          Ez_RD = Ez_zc ( ixxRD, iyyRD, izzRD, 1)
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 2) * xRD
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 3) * yRD
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 4) * zRD

#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 5)
     1          * (xRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 5))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 6)
     1          * (yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 6))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 7)
     1          * (zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 7))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 8)
     1          * (xRD * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 8))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 9)
     1          * (yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 9))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 10)
     1          * (xRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 11)
     1          * (xRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 11))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 12)
     1          * (yRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 12))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 13)
     1          * (zRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 13))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 14)
     1          * (xRD**2 * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 14))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 15)
     1          * (xRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 15))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 16)
     1          * (xRD * yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 16))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 17)
     1          * (yRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 17))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 18)
     1          * (xRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 18))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 19)
     1          * (yRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 19))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 20)
     1          * (xRD * yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 20))
#endif

! ----------

          Ex_RU = Ex_zc ( ixxRU, iyyRU, izzRU, 1)
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 2) * xRU
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 3) * yRU
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 4) * zRU

#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 5)
     1          * (xRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 5))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 6)
     1          * (yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 6))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 7)
     1          * (zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 7))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 8)
     1          * (xRU * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 8))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 9)
     1          * (yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 9))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 10)
     1          * (xRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 11)
     1          * (xRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 11))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 12)
     1          * (yRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 12))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 13)
     1          * (zRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 13))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 14)
     1          * (xRU**2 * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 14))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 15)
     1          * (xRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 15))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 16)
     1          * (xRU * yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 16))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 17)
     1          * (yRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 17))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 18)
     1          * (xRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 18))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 19)
     1          * (yRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 19))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 20)
     1          * (xRU * yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 20))
#endif

          Ey_RU = Ey_zc ( ixxRU, iyyRU, izzRU, 1)
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 2) * xRU
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 3) * yRU
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 4) * zRU

#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 5)
     1          * (xRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 5))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 6)
     1          * (yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 6))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 7)
     1          * (zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 7))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 8)
     1          * (xRU * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 8))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 9)
     1          * (yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 9))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 10)
     1          * (xRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 11)
     1          * (xRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 11))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 12)
     1          * (yRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 12))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 13)
     1          * (zRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 13))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 14)
     1          * (xRU**2 * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 14))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 15)
     1          * (xRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 15))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 16)
     1          * (xRU * yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 16))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 17)
     1          * (yRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 17))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 18)
     1          * (xRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 18))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 19)
     1          * (yRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 19))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 20)
     1          * (xRU * yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 20))
#endif

          Ez_RU = Ez_zc ( ixxRU, iyyRU, izzRU, 1)
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 2) * xRU
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 3) * yRU
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 4) * zRU

#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 5)
     1          * (xRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 5))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 6)
     1          * (yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 6))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 7)
     1          * (zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 7))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 8)
     1          * (xRU * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 8))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 9)
     1          * (yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 9))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 10)
     1          * (xRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 11)
     1          * (xRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 11))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 12)
     1          * (yRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 12))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 13)
     1          * (zRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 13))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 14)
     1          * (xRU**2 * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 14))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 15)
     1          * (xRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 15))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 16)
     1          * (xRU * yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 16))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 17)
     1          * (yRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 17))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 18)
     1          * (xRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 18))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 19)
     1          * (yRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 19))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 20)
     1          * (xRU * yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 20))
#endif

! ----------

          Ex_LU = Ex_zc ( ixxLU, iyyLU, izzLU, 1)
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 2) * xLU
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 3) * yLU
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 4) * zLU

#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 5)
     1          * (xLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 5))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 6)
     1          * (yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 6))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 7)
     1          * (zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 7))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 8)
     1          * (xLU * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 8))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 9)
     1          * (yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 9))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 10)
     1          * (xLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 11)
     1          * (xLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 11))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 12)
     1          * (yLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 12))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 13)
     1          * (zLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 13))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 14)
     1          * (xLU**2 * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 14))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 15)
     1          * (xLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 15))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 16)
     1          * (xLU * yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 16))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 17)
     1          * (yLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 17))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 18)
     1          * (xLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 18))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 19)
     1          * (yLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 19))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 20)
     1          * (xLU * yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 20))
#endif

          Ey_LU = Ey_zc ( ixxLU, iyyLU, izzLU, 1)
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 2) * xLU
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 3) * yLU
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 4) * zLU

#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 5)
     1          * (xLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 5))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 6)
     1          * (yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 6))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 7)
     1          * (zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 7))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 8)
     1          * (xLU * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 8))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 9)
     1          * (yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 9))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 10)
     1          * (xLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 11)
     1          * (xLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 11))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 12)
     1          * (yLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 12))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 13)
     1          * (zLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 13))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 14)
     1          * (xLU**2 * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 14))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 15)
     1          * (xLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 15))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 16)
     1          * (xLU * yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 16))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 17)
     1          * (yLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 17))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 18)
     1          * (xLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 18))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 19)
     1          * (yLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 19))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 20)
     1          * (xLU * yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 20))
#endif

          Ez_LU = Ez_zc ( ixxLU, iyyLU, izzLU, 1)
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 2) * xLU
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 3) * yLU
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 4) * zLU

#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 5)
     1          * (xLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 5))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 6)
     1          * (yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 6))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 7)
     1          * (zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 7))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 8)
     1          * (xLU * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 8))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 9)
     1          * (yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 9))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 10)
     1          * (xLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 11)
     1          * (xLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 11))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 12)
     1          * (yLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 12))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 13)
     1          * (zLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 13))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 14)
     1          * (xLU**2 * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 14))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 15)
     1          * (xLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 15))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 16)
     1          * (xLU * yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 16))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 17)
     1          * (yLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 17))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 18)
     1          * (xLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 18))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 19)
     1          * (yLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 19))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 20)
     1          * (xLU * yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 20))
#endif

! --------------------
#else
! --------------------
! Choose arithmetically averaged electric fields for he centered part of the
! 2D LLF Riemann solver. Electric field reconstruction is not needed.

          Ex_LD = ( Ex_zc ( ixxLD, iyyLD, izzLD, 1)
     1            + Ex_zc ( ixxRD, iyyRD, izzRD, 1)
     1            + Ex_zc ( ixxRU, iyyRU, izzRU, 1)
     1            + Ex_zc ( ixxLU, iyyLU, izzLU, 1) ) * 0.25
          Ex_RD = Ex_LD  ;  Ex_RU = Ex_LD  ;  Ex_LU = Ex_LD

          Ey_LD = ( Ey_zc ( ixxLD, iyyLD, izzLD, 1)
     1            + Ey_zc ( ixxRD, iyyRD, izzRD, 1)
     1            + Ey_zc ( ixxRU, iyyRU, izzRU, 1)
     1            + Ey_zc ( ixxLU, iyyLU, izzLU, 1) ) * 0.25
          Ey_RD = Ey_LD  ;  Ey_RU = Ey_LD  ;  Ey_LU = Ey_LD

          Ez_LD = ( Ez_zc ( ixxLD, iyyLD, izzLD, 1)
     1            + Ez_zc ( ixxRD, iyyRD, izzRD, 1)
     1            + Ez_zc ( ixxRU, iyyRU, izzRU, 1)
     1            + Ez_zc ( ixxLU, iyyLU, izzLU, 1) ) * 0.25
          Ez_RD = Ez_LD  ;  Ez_RU = Ez_LD  ;  Ez_LU = Ez_LD

! --------------------
#endif
! --------------------

! ------------------------------
! Build the electric fields in the local z-direction.

          tempaa = Ex_LD * n3 ( 1) + Ey_LD * n3 ( 2) + Ez_LD * n3 ( 3)
          Ez_LD = tempaa

          tempaa = Ex_RD * n3 ( 1) + Ey_RD * n3 ( 2) + Ez_RD * n3 ( 3)
          Ez_RD = tempaa

          tempaa = Ex_RU * n3 ( 1) + Ey_RU * n3 ( 2) + Ez_RU * n3 ( 3)
          Ez_RU = tempaa

          tempaa = Ex_LU * n3 ( 1) + Ey_LU * n3 ( 2) + Ez_LU * n3 ( 3)
          Ez_LU = tempaa

! ------------------------------
! Build "By_R, By_L" from the RD and LD radial faces respectively.
! These have to include the "Bx0, By0, Bz0" contributions in the "n2" direction.
! They also have to include the time-rate of change so that they can be centered
! in time. This time rate of change is obtained from half and half from
! both sides of the zone boundary.

         By_R = Bzf_quadpts ( ixxRD, iyyRD, izzRD, 5)
     1        + Bx0 * n2 ( 1) + By0 * n2 ( 2) + Bz0 * n2 ( 3)

         By_L = Bzf_quadpts ( ixxLD, iyyLD, izzLD, 4)
     1        + Bx0 * n2 ( 1) + By0 * n2 ( 2) + Bz0 * n2 ( 3)

! ------------------------------
! Build "Bneta_U, Bneta_D" from the LU and LD faces that are in the y-direction
! of the mesh.
! These have to include the "Bx0, By0, Bz0" contributions in the "n1" direction.
! They also have to include the time-rate of change so that they can be centered
! in time. This time rate of change is obtained from half and half from
! both sides of the zone boundary.

      Bneta_U = Byf_quadpts ( ixxLU, iyyLU, izzLU, 3)
     1        + Bx0 * n1 ( 1) + By0 * n1 ( 2) + Bz0 * n1 ( 3)

      Bneta_D = Byf_quadpts ( ixxLD, iyyLD, izzLD, 2)
     1        + Bx0 * n1 ( 1) + By0 * n1 ( 2) + Bz0 * n1 ( 3)

! ------------------------------
! Evaluate "S_R, S_L, S_U, S_D". The Up and Down waves from the y-faces of the mesh
! give us the "S_R, S_L". The The Top and Bottom waves from the z-faces of the mesh
! give us "S_U, S_D".

          S_R = AMAX1 ( S_U_gr ( ixxLU, iyyLU, izzLU),
     1                  S_U_gr ( ixxLD, iyyLD, izzLD) )

          S_L = AMIN1 ( S_D_gr ( ixxLU, iyyLU, izzLU),
     1                  S_D_gr ( ixxLD, iyyLD, izzLD) )

          S_U = AMAX1 ( S_T_gr ( ixxRD, iyyRD, izzRD),
     1                  S_T_gr ( ixxLD, iyyLD, izzLD) )

          S_D = AMIN1 ( S_B_gr ( ixxRD, iyyRD, izzRD),
     1                  S_B_gr ( ixxLD, iyyLD, izzLD) )

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ != 0)
! --------------------
! Make wave speeds fully subsonic for the 2D RS.

          tempaa = AMAX1 ( ABS ( S_R), ABS ( S_L),
     1                     ABS ( S_U), ABS ( S_D) )

          S_R = tempaa  ;  S_L = - tempaa
          S_U = tempaa  ;  S_D = - tempaa

! --------------------
#endif
! --------------------

! ------------------------------------------------------------
! All x-directional calls to the 2D Riemann solver will have four
! input states.

          CALL Two_D_RS_Ez_4State_HLL

     1         ( S_R, S_L, S_U, S_D, cos_t, sin_t,

     1           Ez_RU, Ez_LU, Ez_LD, Ez_RD,

     1           Bneta_U, Bneta_D, By_R, By_L,

     1           BSSx, BSSy, ESSz)

! ------------------------------------------------------------
! Hand back the global electric field "E_z" from the 2D Riemann solver.
! This is a weighted sum over quadrature points.

          Efx ( ixx, iyy, izz) = Efx ( ixx, iyy, izz) + wt * ESSz

! ------------------------------------------------------------

          END DO  ! End of "DO ixx = 1, ntp_zones" loop.

        END DO  ! End of "DO iyy = 0, ntp_zones" loop.

      END DO  ! End of "DO izz = 0, nr_zones" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, 3" loop over quadrature points.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 2) Build "Efy" in the positive y-direction. This electric field is
! zone-centered in the y-direction and edge-centered in the z- and x-directions.
! The global positive y-direction becomes the local z-direction for the
! 2D Riemann solver.

! At third and fourth order, loop over 3 quadrature points in the edge of interest.
      DO j = 1, 3

#if (NUMORDER == 2)
! At second order, we need only one nodal point.
      IF ( j .GT. 1) EXIT
#endif

      wt = 1.0
#if (NUMORDER >= 3)
      wt = 1.0 / 6.0  ! These are the weights for simpson integration.
      IF ( j .EQ. 1) wt = 4.0 / 6.0
#endif

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxLD, iyyLD, izzLD,         &
!$OMP&    ixxRD, iyyRD, izzRD, ixxRU, iyyRU, izzRU, ixxLU, iyyLU, izzLU,&

!$OMP&     xLD, yLD, zLD, xRD, yRD, zRD, xRU, yRU, zRU, xLU, yLU, zLU,  &
!$OMP&     oblLD, oblRD, oblRU, oblLU, xec, yec, zec, tempaa,           &
!$OMP&     bflag, time_use, Bx0, By0, Bz0,                              &
!$OMP&     Ex_LD, Ey_LD, Ez_LD, Ex_RD, Ey_RD, Ez_RD,                    &
!$OMP&     Ex_LU, Ey_LU, Ez_LU, Ex_RU, Ey_RU, Ez_RU,                    &

!$OMP&     n1, n2, n3, B0, S_R, S_L, S_U, S_D, cos_t, sin_t,            &
!$OMP&     Bneta_U, Bneta_D, By_R, By_L,                                &
!$OMP&     BSSx, BSSy, ESSz)                                            &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum, j, wt,&
!$OMP&    vert_coords_patch, yf_quadnorml_patch, v_scale_length,        &
!$OMP&    v_centroid_patch, Coef_Basis_Moments,                         &

!$OMP&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$OMP&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$OMP&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$OMP&    Efy )                                                         &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 3)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum, j, wt,                                              &
!$acc&    vert_coords_patch, yf_quadnorml_patch, v_scale_length,        &
!$acc&    v_centroid_patch, Coef_Basis_Moments,                         &

!$acc&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$acc&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$acc&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$acc&    Efy )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxLD, iyyLD, izzLD,                 &
!$acc&    ixxRD, iyyRD, izzRD, ixxRU, iyyRU, izzRU, ixxLU, iyyLU, izzLU,&

!$acc&     xLD, yLD, zLD, xRD, yRD, zRD, xRU, yRU, zRU, xLU, yLU, zLU,  &
!$acc&     oblLD, oblRD, oblRU, oblLU, xec, yec, zec, tempaa,           &
!$acc&     bflag, time_use, Bx0, By0, Bz0,                              &
!$acc&     Ex_LD, Ey_LD, Ez_LD, Ex_RD, Ey_RD, Ez_RD,                    &
!$acc&     Ex_LU, Ey_LU, Ez_LU, Ex_RU, Ey_RU, Ez_RU,                    &

!$acc&     n1, n2, n3, B0, S_R, S_L, S_U, S_D, cos_t, sin_t,            &
!$acc&     Bneta_U, Bneta_D, By_R, By_L,                                &
!$acc&     BSSx, BSSy, ESSz)

      DO izz = 0, nr_zones

        DO ixx = 0, ntp_zones

          DO iyy = 1, ntp_zones

! ------------------------------------------------------------

! Fill in "ixxLD, iyyLD, izzLD, ixxRD, iyyRD, izzRD,
! ixxRU, iyyRU, izzRU, ixxLU, iyyLU, izzLU" in counter-clockwise sense
! relative to the positive y-axis.

          ixxLD = ixx
          iyyLD = iyy
          izzLD = izz

          ixxRD = ixx
          iyyRD = iyy
          izzRD = izz + 1

          ixxLU = ixx + 1
          iyyLU = iyy
          izzLU = izz

          ixxRU = ixx + 1
          iyyRU = iyy
          izzRU = izz + 1

! --------------------------------------------------
! Find the edge-center "xec, yec, zec" by bisecting the two y-vertices.
! Then scale it up to be of the correct radius.
! These are needed in general if we are using the reconstruction.

          IF ( j .EQ. 1) THEN  ! First quadrature point is the center point.

            xec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 1)  ! Bisect the y-vertices.
     1                  + vert_coords_patch ( ixx, iyy - 1, izz, 1) )
            yec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 2)
     1                  + vert_coords_patch ( ixx, iyy - 1, izz, 2) )
            zec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 3)
     1                  + vert_coords_patch ( ixx, iyy - 1, izz, 3) )

            tempaa = SQRT ( vert_coords_patch ( ixx, iyy, izz, 1)**2
     1                    + vert_coords_patch ( ixx, iyy, izz, 2)**2
     1                    + vert_coords_patch ( ixx, iyy, izz, 3)**2)
     1             / SQRT ( xec**2 + yec**2 + zec**2)

            xec = xec * tempaa  ! Rescale to have correct radius.
            yec = yec * tempaa
            zec = zec * tempaa

          ELSE IF ( j .EQ. 2) THEN  ! Second quadrature point is one of the end points.

            xec = vert_coords_patch ( ixx, iyy, izz, 1)
            yec = vert_coords_patch ( ixx, iyy, izz, 2)
            zec = vert_coords_patch ( ixx, iyy, izz, 3)

          ELSE IF ( j .EQ. 3) THEN  ! Third quadrature point is the other end point.

            xec = vert_coords_patch ( ixx, iyy - 1, izz, 1)
            yec = vert_coords_patch ( ixx, iyy - 1, izz, 2)
            zec = vert_coords_patch ( ixx, iyy - 1, izz, 3)

          END IF

! ------------------------------
! Feed ( xec, yec, zec) to the tiny subroutine that evaluates "B0".
! Default is zero.

          B0 = 0.0
          Bx0 = 0.0
          By0 = 0.0
          Bz0 = 0.0

#if ( USE_B0B1_FORMULATION == 1)

          bflag = 1
          time_use = 0

          CALL GET_MAGNETIC_FIELD ( bflag, time_use,  xec, yec, zec,
     1                              Bx0, By0, Bz0)

          B0 ( 1) = Bx0
          B0 ( 2) = By0
          B0 ( 3) = Bz0

#endif

! ------------------------------
! Build "oblLD, oblRD, oblRU, oblLU" so that we can interpolate in the local
! coordinates of each of the four participating zones.
! These are needed in general if we are using the reconstruction.

          oblLD = 1.0 / AMAX1 ( v_scale_length ( ixxLD, iyyLD, izzLD),
     1                          smallnum)

          oblRD = 1.0 / AMAX1 ( v_scale_length ( ixxRD, iyyRD, izzRD),
     1                          smallnum)

          oblRU = 1.0 / AMAX1 ( v_scale_length ( ixxRU, iyyRU, izzRU),
     1                          smallnum)

          oblLU = 1.0 / AMAX1 ( v_scale_length ( ixxLU, iyyLU, izzLU),
     1                          smallnum)

! ------------------------------
! Build scaled "xLD, yLD, zLD, xRD, yRD, zRD, xRU, yRU, zRU, xLU, yLU, zLU".
! These are needed in general if we are using the reconstruction.

! ----------

          xLD = ( xec - v_centroid_patch ( ixxLD, iyyLD, izzLD, 1) )
     1        * oblLD
          yLD = ( yec - v_centroid_patch ( ixxLD, iyyLD, izzLD, 2) )
     1        * oblLD
          zLD = ( zec - v_centroid_patch ( ixxLD, iyyLD, izzLD, 3) )
     1        * oblLD

! ----------

          xRD = ( xec - v_centroid_patch ( ixxRD, iyyRD, izzRD, 1) )
     1        * oblRD
          yRD = ( yec - v_centroid_patch ( ixxRD, iyyRD, izzRD, 2) )
     1        * oblRD
          zRD = ( zec - v_centroid_patch ( ixxRD, iyyRD, izzRD, 3) )
     1        * oblRD

! ----------

          xRU = ( xec - v_centroid_patch ( ixxRU, iyyRU, izzRU, 1) )
     1        * oblRU
          yRU = ( yec - v_centroid_patch ( ixxRU, iyyRU, izzRU, 2) )
     1        * oblRU
          zRU = ( zec - v_centroid_patch ( ixxRU, iyyRU, izzRU, 3) )
     1        * oblRU

! ----------

          xLU = ( xec - v_centroid_patch ( ixxLU, iyyLU, izzLU, 1) )
     1        * oblLU
          yLU = ( yec - v_centroid_patch ( ixxLU, iyyLU, izzLU, 2) )
     1        * oblLU
          zLU = ( zec - v_centroid_patch ( ixxLU, iyyLU, izzLU, 3) )
     1        * oblLU

! ------------------------------
! Build local normals so that globally positive y-direction is indeed the
! locally positive z-direction.
! "n1" is a vector in the radial direction; which becomes the local x-direction.
! "n2" is a vector along the x-direction of the mesh; which becomes the local y-direction.
! "n3" is a vector along the y-direction of the mesh; which becomes the local z-direction.

          tempaa = 1.0 / SQRT ( xec**2 + yec**2 + zec**2)

          n1 ( 1) = xec * tempaa
          n1 ( 2) = yec * tempaa
          n1 ( 3) = zec * tempaa

          n3 ( 1) = vert_coords_patch ( ixx, iyy, izz, 1)
     1            - vert_coords_patch ( ixx, iyy - 1, izz, 1)
          n3 ( 2) = vert_coords_patch ( ixx, iyy, izz, 2)
     1            - vert_coords_patch ( ixx, iyy - 1, izz, 2)
          n3 ( 3) = vert_coords_patch ( ixx, iyy, izz, 3)
     1            - vert_coords_patch ( ixx, iyy - 1, izz, 3)

          tempaa = n3 ( 1) * n1 ( 1) + n3 ( 2) * n1 ( 2)
     1           + n3 ( 3) * n1 ( 3)

          n3 ( 1) = n3 ( 1) - tempaa * n1 ( 1)
          n3 ( 2) = n3 ( 2) - tempaa * n1 ( 2)
          n3 ( 3) = n3 ( 3) - tempaa * n1 ( 3)

          tempaa = 1.0 / SQRT ( n3 ( 1)**2 + n3 ( 2)**2 + n3 ( 3)**2)

          n3 ( :) = n3 ( :) * tempaa

! "n2 = n3 X n1".
          n2 ( 1) = n3 ( 2) * n1 ( 3) - n3 ( 3) * n1 ( 2)
          n2 ( 2) = n3 ( 3) * n1 ( 1) - n3 ( 1) * n1 ( 3)
          n2 ( 3) = n3 ( 1) * n1 ( 2) - n3 ( 2) * n1 ( 1)

! The local x and y-axes are orthogonal.
          cos_t = 0.0  ;  sin_t = 1.0

! ------------------------------
! Now that we have the properly scaled relative coordinates from each of the
! four particpating zones, we can use these to obtain the four input electric field components.
! "Ex,y,z_RU, Ex,y,z_LU, Ex,y,z_LD, Ex,y,z_RD" that go into the 2D Riemann solver.

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ == 0)
! --------------------
! Here the building of the electric fields requires that they have been reconstructed.

          Ex_LD = Ex_zc ( ixxLD, iyyLD, izzLD, 1)
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 2) * xLD
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 3) * yLD
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 4) * zLD

#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 5)
     1          * (xLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 5))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 6)
     1          * (yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 6))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 7)
     1          * (zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 7))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 8)
     1          * (xLD * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 8))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 9)
     1          * (yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 9))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 10)
     1          * (xLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 11)
     1          * (xLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 11))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 12)
     1          * (yLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 12))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 13)
     1          * (zLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 13))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 14)
     1          * (xLD**2 * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 14))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 15)
     1          * (xLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 15))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 16)
     1          * (xLD * yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 16))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 17)
     1          * (yLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 17))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 18)
     1          * (xLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 18))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 19)
     1          * (yLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 19))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 20)
     1          * (xLD * yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 20))

#endif

          Ey_LD = Ey_zc ( ixxLD, iyyLD, izzLD, 1)
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 2) * xLD
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 3) * yLD
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 4) * zLD
#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 5)
     1          * (xLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 5))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 6)
     1          * (yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 6))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 7)
     1          * (zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 7))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 8)
     1          * (xLD * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 8))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 9)
     1          * (yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 9))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 10)
     1          * (xLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 11)
     1          * (xLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 11))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 12)
     1          * (yLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 12))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 13)
     1          * (zLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 13))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 14)
     1          * (xLD**2 * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 14))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 15)
     1          * (xLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 15))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 16)
     1          * (xLD * yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 16))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 17)
     1          * (yLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 17))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 18)
     1          * (xLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 18))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 19)
     1          * (yLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 19))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 20)
     1          * (xLD * yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 20))

#endif

          Ez_LD = Ez_zc ( ixxLD, iyyLD, izzLD, 1)
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 2) * xLD
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 3) * yLD
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 4) * zLD
#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 5)
     1          * (xLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 5))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 6)
     1          * (yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 6))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 7)
     1          * (zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 7))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 8)
     1          * (xLD * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 8))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 9)
     1          * (yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 9))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 10)
     1          * (xLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 11)
     1          * (xLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 11))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 12)
     1          * (yLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 12))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 13)
     1          * (zLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 13))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 14)
     1          * (xLD**2 * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 14))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 15)
     1          * (xLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 15))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 16)
     1          * (xLD * yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 16))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 17)
     1          * (yLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 17))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 18)
     1          * (xLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 18))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 19)
     1          * (yLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 19))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 20)
     1          * (xLD * yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 20))

#endif

! ----------

          Ex_RD = Ex_zc ( ixxRD, iyyRD, izzRD, 1)
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 2) * xRD
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 3) * yRD
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 4) * zRD
#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 5)
     1          * (xRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 5))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 6)
     1          * (yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 6))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 7)
     1          * (zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 7))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 8)
     1          * (xRD * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 8))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 9)
     1          * (yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 9))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 10)
     1          * (xRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 11)
     1          * (xRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 11))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 12)
     1          * (yRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 12))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 13)
     1          * (zRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 13))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 14)
     1          * (xRD**2 * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 14))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 15)
     1          * (xRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 15))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 16)
     1          * (xRD * yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 16))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 17)
     1          * (yRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 17))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 18)
     1          * (xRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 18))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 19)
     1          * (yRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 19))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 20)
     1          * (xRD * yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 20))

#endif

          Ey_RD = Ey_zc ( ixxRD, iyyRD, izzRD, 1)
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 2) * xRD
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 3) * yRD
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 4) * zRD
#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 5)
     1          * (xRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 5))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 6)
     1          * (yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 6))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 7)
     1          * (zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 7))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 8)
     1          * (xRD * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 8))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 9)
     1          * (yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 9))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 10)
     1          * (xRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 11)
     1          * (xRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 11))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 12)
     1          * (yRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 12))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 13)
     1          * (zRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 13))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 14)
     1          * (xRD**2 * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 14))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 15)
     1          * (xRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 15))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 16)
     1          * (xRD * yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 16))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 17)
     1          * (yRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 17))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 18)
     1          * (xRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 18))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 19)
     1          * (yRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 19))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 20)
     1          * (xRD * yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 20))

#endif

          Ez_RD = Ez_zc ( ixxRD, iyyRD, izzRD, 1)
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 2) * xRD
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 3) * yRD
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 4) * zRD
#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 5)
     1          * (xRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 5))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 6)
     1          * (yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 6))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 7)
     1          * (zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 7))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 8)
     1          * (xRD * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 8))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 9)
     1          * (yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 9))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 10)
     1          * (xRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 11)
     1          * (xRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 11))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 12)
     1          * (yRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 12))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 13)
     1          * (zRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 13))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 14)
     1          * (xRD**2 * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 14))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 15)
     1          * (xRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 15))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 16)
     1          * (xRD * yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 16))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 17)
     1          * (yRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 17))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 18)
     1          * (xRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 18))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 19)
     1          * (yRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 19))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 20)
     1          * (xRD * yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 20))

#endif

! ----------

          Ex_RU = Ex_zc ( ixxRU, iyyRU, izzRU, 1)
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 2) * xRU
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 3) * yRU
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 4) * zRU
#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 5)
     1          * (xRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 5))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 6)
     1          * (yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 6))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 7)
     1          * (zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 7))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 8)
     1          * (xRU * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 8))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 9)
     1          * (yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 9))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 10)
     1          * (xRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 11)
     1          * (xRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 11))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 12)
     1          * (yRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 12))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 13)
     1          * (zRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 13))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 14)
     1          * (xRU**2 * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 14))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 15)
     1          * (xRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 15))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 16)
     1          * (xRU * yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 16))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 17)
     1          * (yRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 17))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 18)
     1          * (xRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 18))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 19)
     1          * (yRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 19))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 20)
     1          * (xRU * yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 20))

#endif

          Ey_RU = Ey_zc ( ixxRU, iyyRU, izzRU, 1)
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 2) * xRU
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 3) * yRU
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 4) * zRU
#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 5)
     1          * (xRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 5))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 6)
     1          * (yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 6))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 7)
     1          * (zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 7))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 8)
     1          * (xRU * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 8))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 9)
     1          * (yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 9))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 10)
     1          * (xRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 11)
     1          * (xRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 11))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 12)
     1          * (yRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 12))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 13)
     1          * (zRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 13))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 14)
     1          * (xRU**2 * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 14))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 15)
     1          * (xRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 15))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 16)
     1          * (xRU * yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 16))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 17)
     1          * (yRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 17))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 18)
     1          * (xRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 18))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 19)
     1          * (yRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 19))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 20)
     1          * (xRU * yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 20))

#endif

          Ez_RU = Ez_zc ( ixxRU, iyyRU, izzRU, 1)
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 2) * xRU
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 3) * yRU
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 4) * zRU
#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 5)
     1          * (xRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 5))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 6)
     1          * (yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 6))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 7)
     1          * (zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 7))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 8)
     1          * (xRU * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 8))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 9)
     1          * (yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 9))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 10)
     1          * (xRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 11)
     1          * (xRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 11))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 12)
     1          * (yRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 12))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 13)
     1          * (zRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 13))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 14)
     1          * (xRU**2 * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 14))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 15)
     1          * (xRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 15))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 16)
     1          * (xRU * yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 16))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 17)
     1          * (yRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 17))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 18)
     1          * (xRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 18))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 19)
     1          * (yRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 19))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 20)
     1          * (xRU * yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 20))

#endif

! ----------

          Ex_LU = Ex_zc ( ixxLU, iyyLU, izzLU, 1)
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 2) * xLU
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 3) * yLU
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 4) * zLU
#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 5)
     1          * (xLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 5))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 6)
     1          * (yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 6))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 7)
     1          * (zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 7))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 8)
     1          * (xLU * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 8))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 9)
     1          * (yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 9))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 10)
     1          * (xLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 11)
     1          * (xLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 11))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 12)
     1          * (yLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 12))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 13)
     1          * (zLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 13))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 14)
     1          * (xLU**2 * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 14))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 15)
     1          * (xLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 15))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 16)
     1          * (xLU * yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 16))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 17)
     1          * (yLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 17))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 18)
     1          * (xLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 18))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 19)
     1          * (yLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 19))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 20)
     1          * (xLU * yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 20))

#endif

          Ey_LU = Ey_zc ( ixxLU, iyyLU, izzLU, 1)
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 2) * xLU
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 3) * yLU
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 4) * zLU
#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 5)
     1          * (xLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 5))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 6)
     1          * (yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 6))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 7)
     1          * (zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 7))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 8)
     1          * (xLU * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 8))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 9)
     1          * (yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 9))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 10)
     1          * (xLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 11)
     1          * (xLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 11))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 12)
     1          * (yLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 12))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 13)
     1          * (zLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 13))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 14)
     1          * (xLU**2 * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 14))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 15)
     1          * (xLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 15))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 16)
     1          * (xLU * yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 16))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 17)
     1          * (yLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 17))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 18)
     1          * (xLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 18))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 19)
     1          * (yLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 19))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 20)
     1          * (xLU * yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 20))

#endif

          Ez_LU = Ez_zc ( ixxLU, iyyLU, izzLU, 1)
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 2) * xLU
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 3) * yLU
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 4) * zLU
#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 5)
     1          * (xLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 5))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 6)
     1          * (yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 6))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 7)
     1          * (zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 7))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 8)
     1          * (xLU * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 8))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 9)
     1          * (yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 9))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 10)
     1          * (xLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 11)
     1          * (xLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 11))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 12)
     1          * (yLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 12))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 13)
     1          * (zLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 13))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 14)
     1          * (xLU**2 * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 14))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 15)
     1          * (xLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 15))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 16)
     1          * (xLU * yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 16))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 17)
     1          * (yLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 17))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 18)
     1          * (xLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 18))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 19)
     1          * (yLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 19))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 20)
     1          * (xLU * yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 20))

#endif

! --------------------
#else
! --------------------
! Choose arithmetically averaged electric fields for he centered part of the
! 2D LLF Riemann solver. Electric field reconstruction is not needed.

          Ex_LD = ( Ex_zc ( ixxLD, iyyLD, izzLD, 1)
     1            + Ex_zc ( ixxRD, iyyRD, izzRD, 1)
     1            + Ex_zc ( ixxRU, iyyRU, izzRU, 1)
     1            + Ex_zc ( ixxLU, iyyLU, izzLU, 1) ) * 0.25
          Ex_RD = Ex_LD  ;  Ex_RU = Ex_LD  ;  Ex_LU = Ex_LD

          Ey_LD = ( Ey_zc ( ixxLD, iyyLD, izzLD, 1)
     1            + Ey_zc ( ixxRD, iyyRD, izzRD, 1)
     1            + Ey_zc ( ixxRU, iyyRU, izzRU, 1)
     1            + Ey_zc ( ixxLU, iyyLU, izzLU, 1) ) * 0.25
          Ey_RD = Ey_LD  ;  Ey_RU = Ey_LD  ;  Ey_LU = Ey_LD

          Ez_LD = ( Ez_zc ( ixxLD, iyyLD, izzLD, 1)
     1            + Ez_zc ( ixxRD, iyyRD, izzRD, 1)
     1            + Ez_zc ( ixxRU, iyyRU, izzRU, 1)
     1            + Ez_zc ( ixxLU, iyyLU, izzLU, 1) ) * 0.25
          Ez_RD = Ez_LD  ;  Ez_RU = Ez_LD  ;  Ez_LU = Ez_LD

! --------------------
#endif
! --------------------

! ------------------------------
! Build the electric fields in the local z-direction.

          tempaa = Ex_LD * n3 ( 1) + Ey_LD * n3 ( 2) + Ez_LD * n3 ( 3)
          Ez_LD = tempaa

          tempaa = Ex_RD * n3 ( 1) + Ey_RD * n3 ( 2) + Ez_RD * n3 ( 3)
          Ez_RD = tempaa

          tempaa = Ex_RU * n3 ( 1) + Ey_RU * n3 ( 2) + Ez_RU * n3 ( 3)
          Ez_RU = tempaa

          tempaa = Ex_LU * n3 ( 1) + Ey_LU * n3 ( 2) + Ez_LU * n3 ( 3)
          Ez_LU = tempaa

! ------------------------------
! Build "By_R, By_L" from the RD and LD radial faces respectively.
! These have to include the "Bx0, By0, Bz0" contributions in the "n2" direction.

         By_R = Bxf_quadpts ( ixxRD, iyyRD, izzRD, 5)
     1        + Bx0 * n2 ( 1) + By0 * n2 ( 2) + Bz0 * n2 ( 3)

         By_L = Bxf_quadpts ( ixxLD, iyyLD, izzLD, 4)
     1        + Bx0 * n2 ( 1) + By0 * n2 ( 2) + Bz0 * n2 ( 3)

! ------------------------------
! Build "Bneta_U, Bneta_D" from RU and RD faces that are in the negative x-direction
! of the mesh.
! These have to include the "Bx0, By0, Bz0" contributions in the "n1" direction.

          Bneta_U = Bzf_quadpts ( ixxLU, iyyLU, izzLU, 3)
     1            + Bx0 * n1 ( 1) + By0 * n1 ( 2) + Bz0 * n1 ( 3)

          Bneta_D = Bzf_quadpts ( ixxLD, iyyLD, izzLD, 2)
     1            + Bx0 * n1 ( 1) + By0 * n1 ( 2) + Bz0 * n1 ( 3)

! ------------------------------
! Evaluate "S_R, S_L, S_U, S_D". The minus of Left and Right waves from the x-faces
! of the mesh give us "S_R, S_L". The Top and Bottom waves from the z-faces of the mesh
! give us "S_U, S_D".

          S_R = AMAX1 ( S_T_gr ( ixxLU, iyyLU, izzLU),
     1                  S_T_gr ( ixxLD, iyyLD, izzLD) )

          S_L = AMIN1 ( S_B_gr ( ixxLU, iyyLU, izzLU),
     1                  S_B_gr ( ixxLD, iyyLD, izzLD) )

          S_U = AMAX1 ( S_R_gr ( ixxLD, iyyLD, izzLD),
     1                  S_R_gr ( ixxRD, iyyRD, izzRD) )

          S_D = AMIN1 ( S_L_gr ( ixxLD, iyyLD, izzLD),
     1                  S_L_gr ( ixxRD, iyyRD, izzRD) )

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ != 0)
! --------------------
! Make wave speeds fully subsonic for the 2D RS.

          tempaa = AMAX1 ( ABS ( S_R), ABS ( S_L),
     1                     ABS ( S_U), ABS ( S_D) )

          S_R = tempaa  ;  S_L = - tempaa
          S_U = tempaa  ;  S_D = - tempaa

! --------------------
#endif
! --------------------

! ------------------------------------------------------------
! All y-directional calls to the 2D Riemann solver will have four
! input states.

          CALL Two_D_RS_Ez_4State_HLL

     1         ( S_R, S_L, S_U, S_D, cos_t, sin_t,

     1           Ez_RU, Ez_LU, Ez_LD, Ez_RD,

     1           Bneta_U, Bneta_D, By_R, By_L,

     1           BSSx, BSSy, ESSz)

! ------------------------------------------------------------
! Hand back the global electric field "E_z" from the 2D Riemann solver.

          Efy ( ixx, iyy, izz) = Efy ( ixx, iyy, izz) + wt * ESSz

! ------------------------------------------------------------

          END DO  ! End of "DO iyy = 1, ntp_zones" loop.

        END DO  ! End of "DO ixx = 0, ntp_zones" loop.

      END DO  ! End of "DO izz = 0, nr_zones" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, 3" loop over quadrature points.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 3) Build "Efz" in the positive z-direction. This electric field is
! zone-centered in the z-direction and edge-centered in the x- and y-directions.
! The global positive z-direction becomes the local z-direction for the
! 2D Riemann solver.

! Here we only treat quad-corners. The tri-corners are treated separately.

! At third and fourth order, loop over 3 quadrature points in the edge of interest.
      DO j = 1, 3

#if (NUMORDER == 2)
! At second order, we need only one nodal point.
      IF ( j .GT. 1) EXIT
#endif

      wt = 1.0
#if (NUMORDER >= 3)
      wt = 1.0 / 6.0  ! These are the weights for simpson integration.
      IF ( j .EQ. 1) wt = 4.0 / 6.0
#endif

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxLD, iyyLD, izzLD,         &
!$OMP&    ixxRD, iyyRD, izzRD, ixxRU, iyyRU, izzRU, ixxLU, iyyLU, izzLU,&

!$OMP&     xLD, yLD, zLD, xRD, yRD, zRD, xRU, yRU, zRU, xLU, yLU, zLU,  &
!$OMP&     oblLD, oblRD, oblRU, oblLU, xec, yec, zec, tempaa,           &
!$OMP&     bflag, time_use, Bx0, By0, Bz0, tau, neta,                   &
!$OMP&     Ex_LD, Ey_LD, Ez_LD, Ex_RD, Ey_RD, Ez_RD,                    &
!$OMP&     Ex_LU, Ey_LU, Ez_LU, Ex_RU, Ey_RU, Ez_RU,                    &

!$OMP&     n1, n2, n3, B0, S_R, S_L, S_U, S_D, cos_t, sin_t,            &
!$OMP&     Bneta_U, Bneta_D, By_R, By_L,                                &
!$OMP&     BSSx, BSSy, ESSz)                                            &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum, j, wt,&
!$OMP&    vert_coords_patch, yf_quadnorml_patch, v_scale_length,        &
!$OMP&    v_centroid_patch, Coef_Basis_Moments, imask_active_patch,     &

!$OMP&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$OMP&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$OMP&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$OMP&    Efz )                                                         &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 3)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum, j, wt,                                              &
!$acc&    vert_coords_patch, yf_quadnorml_patch, v_scale_length,        &
!$acc&    v_centroid_patch, Coef_Basis_Moments, imask_active_patch,     &

!$acc&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$acc&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$acc&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$acc&    Efz)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxLD, iyyLD, izzLD,                 &
!$acc&    ixxRD, iyyRD, izzRD, ixxRU, iyyRU, izzRU, ixxLU, iyyLU, izzLU,&

!$acc&     xLD, yLD, zLD, xRD, yRD, zRD, xRU, yRU, zRU, xLU, yLU, zLU,  &
!$acc&     oblLD, oblRD, oblRU, oblLU, xec, yec, zec, tempaa,           &
!$acc&     bflag, time_use, Bx0, By0, Bz0, tau, neta,                   &
!$acc&     Ex_LD, Ey_LD, Ez_LD, Ex_RD, Ey_RD, Ez_RD,                    &
!$acc&     Ex_LU, Ey_LU, Ez_LU, Ex_RU, Ey_RU, Ez_RU,                    &

!$acc&     n1, n2, n3, B0, S_R, S_L, S_U, S_D, cos_t, sin_t,            &
!$acc&     Bneta_U, Bneta_D, By_R, By_L,                                &
!$acc&     BSSx, BSSy, ESSz)                                            &

      DO ixx = 0, ntp_zones

        DO iyy = 0, ntp_zones

          DO izz = 1, nr_zones

! ------------------------------------------------------------

! Fill in "ixxLD, iyyLD, izzLD, ixxRD, iyyRD, izzRD,
! ixxRU, iyyRU, izzRU, ixxLU, iyyLU, izzLU" in counter-clockwise sense
! relative to the z-axis.

          ixxLD = ixx
          iyyLD = iyy
          izzLD = izz

          ixxRD = ixx + 1
          iyyRD = iyy
          izzRD = izz

          ixxRU = ixx + 1
          iyyRU = iyy + 1
          izzRU = izz

          ixxLU = ixx
          iyyLU = iyy + 1
          izzLU = izz

! Prevent processing of tri-corners.

          IF ( ( imask_active_patch ( ixxLD, iyyLD, izzLD) .EQ. 0) .OR.
     1         ( imask_active_patch ( ixxRD, iyyRD, izzRD) .EQ. 0) .OR.
     1         ( imask_active_patch ( ixxRU, iyyRU, izzRU) .EQ. 0) .OR.
     1         ( imask_active_patch ( ixxLU, iyyLU, izzLU) .EQ. 0)) THEN

          CYCLE

          END IF

! --------------------------------------------------
! Find the edge-center "xec, yec, zec" by bisecting the two z-vertices.
! No need to rescale by radius since this is in the radial direction by definition.

          IF ( j .EQ. 1) THEN  ! First quadrature point is the center point.

            xec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 1)  ! Bisect the z-vertices.
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 1) )
            yec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 2)
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 2) )
            zec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 3)
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 3) )

          ELSE IF ( j .EQ. 2) THEN  ! Second quadrature point is one of the end points.

            xec = vert_coords_patch ( ixx, iyy, izz, 1)
            yec = vert_coords_patch ( ixx, iyy, izz, 2)
            zec = vert_coords_patch ( ixx, iyy, izz, 3)

          ELSE IF ( j .EQ. 3) THEN  ! Third quadrature point is the other end point.

            xec = vert_coords_patch ( ixx, iyy, izz - 1, 1)
            yec = vert_coords_patch ( ixx, iyy, izz - 1, 2)
            zec = vert_coords_patch ( ixx, iyy, izz - 1, 3)

          END IF

! ------------------------------
! Feed ( xec, yec, zec) to the tiny subroutine that evaluates "B0".
! Default is zero.

          B0 = 0.0
          Bx0 = 0.0
          By0 = 0.0
          Bz0 = 0.0

#if ( USE_B0B1_FORMULATION == 1)

          bflag = 1
          time_use = 0

          CALL GET_MAGNETIC_FIELD ( bflag, time_use,  xec, yec, zec,
     1                              Bx0, By0, Bz0)

          B0 ( 1) = Bx0
          B0 ( 2) = By0
          B0 ( 3) = Bz0

#endif

! ------------------------------
! Build "oblLD, oblRD, oblRU, oblLU" so that we can interpolate in the local
! coordinates of each of the four participating zones.
! These are needed in general if we are using the reconstruction.

          oblLD = 1.0 / AMAX1 ( v_scale_length ( ixxLD, iyyLD, izzLD),
     1                          smallnum)

          oblRD = 1.0 / AMAX1 ( v_scale_length ( ixxRD, iyyRD, izzRD),
     1                          smallnum)

          oblRU = 1.0 / AMAX1 ( v_scale_length ( ixxRU, iyyRU, izzRU),
     1                          smallnum)

          oblLU = 1.0 / AMAX1 ( v_scale_length ( ixxLU, iyyLU, izzLU),
     1                          smallnum)

! ------------------------------
! Build scaled "xLD, yLD, zLD, xRD, yRD, zRD, xRU, yRU, zRU, xLU, yLU, zLU".
! These are needed in general if we are using the reconstruction.

! ----------

          xLD = ( xec - v_centroid_patch ( ixxLD, iyyLD, izzLD, 1) )
     1        * oblLD
          yLD = ( yec - v_centroid_patch ( ixxLD, iyyLD, izzLD, 2) )
     1        * oblLD
          zLD = ( zec - v_centroid_patch ( ixxLD, iyyLD, izzLD, 3) )
     1        * oblLD

! ----------

          xRD = ( xec - v_centroid_patch ( ixxRD, iyyRD, izzRD, 1) )
     1        * oblRD
          yRD = ( yec - v_centroid_patch ( ixxRD, iyyRD, izzRD, 2) )
     1        * oblRD
          zRD = ( zec - v_centroid_patch ( ixxRD, iyyRD, izzRD, 3) )
     1        * oblRD

! ----------

          xRU = ( xec - v_centroid_patch ( ixxRU, iyyRU, izzRU, 1) )
     1        * oblRU
          yRU = ( yec - v_centroid_patch ( ixxRU, iyyRU, izzRU, 2) )
     1        * oblRU
          zRU = ( zec - v_centroid_patch ( ixxRU, iyyRU, izzRU, 3) )
     1        * oblRU

! ----------

          xLU = ( xec - v_centroid_patch ( ixxLU, iyyLU, izzLU, 1) )
     1        * oblLU
          yLU = ( yec - v_centroid_patch ( ixxLU, iyyLU, izzLU, 2) )
     1        * oblLU
          zLU = ( zec - v_centroid_patch ( ixxLU, iyyLU, izzLU, 3) )
     1        * oblLU

! ------------------------------
! Build local normals so that globally positive x-direction is indeed the
! locally positive x-direction. This is "n1". The radial direction is
! the locally positive z-direction. This is "n3". "n2 = n3 X n1" is
! the locally positive y-direction. But it may not coincide with the
! actual y-axis of the mesh. "tau" is aligned with the actual y-axis of
! the mesh, and "neta = tau X n3" is perpendicular to it.

          tempaa = 1.0 / SQRT ( xec**2 + yec**2 + zec**2)

          n3 ( 1) = xec * tempaa
          n3 ( 2) = yec * tempaa
          n3 ( 3) = zec * tempaa

! Obtain "n1" in the xy-tangent space of the sphere in the x-coordinate direction.
! Our construction is such that "n1" will always be perpendicular to "n3".

          n1 ( 1) = vert_coords_patch ( ixx + 1, iyy, izz, 1)
     1            - vert_coords_patch ( ixx - 1, iyy, izz, 1)
          n1 ( 2) = vert_coords_patch ( ixx + 1, iyy, izz, 2)
     1            - vert_coords_patch ( ixx - 1, iyy, izz, 2)
          n1 ( 3) = vert_coords_patch ( ixx + 1, iyy, izz, 3)
     1            - vert_coords_patch ( ixx - 1, iyy, izz, 3)

          tempaa = n1 ( 1) * n3 ( 1) + n1 ( 2) * n3 ( 2)
     1           + n1 ( 3) * n3 ( 3)

          n1 ( 1) = n1 ( 1) - tempaa * n3 ( 1)
          n1 ( 2) = n1 ( 2) - tempaa * n3 ( 2)
          n1 ( 3) = n1 ( 3) - tempaa * n3 ( 3)

          tempaa = 1.0 / SQRT ( n1 ( 1)**2 + n1 ( 2)**2 + n1 ( 3)**2)

          n1 ( :) = n1 ( :) * tempaa

! Obtain "n2 = n3 X n1" in the xy-tangent space of the sphere. This usually
! will not be the y-direction on the mesh because the xy-faces can be diamond shaped.

          n2 ( 1) = n3 ( 2) * n1 ( 3) - n3 ( 3) * n1 ( 2)
          n2 ( 2) = n3 ( 3) * n1 ( 1) - n3 ( 1) * n1 ( 3)
          n2 ( 3) = n3 ( 1) * n1 ( 2) - n3 ( 2) * n1 ( 1)

! ----------
! Obtain "tau" in the xy-tangent space of the sphere in the y-coordinate direction.
! Our construction is such that "tau" will always be perpendicular to "n3".

          tau ( 1) = vert_coords_patch ( ixx, iyy + 1, izz, 1)
     1            - vert_coords_patch ( ixx, iyy - 1, izz, 1)
          tau ( 2) = vert_coords_patch ( ixx, iyy + 1, izz, 2)
     1            - vert_coords_patch ( ixx, iyy - 1, izz, 2)
          tau ( 3) = vert_coords_patch ( ixx, iyy + 1, izz, 3)
     1            - vert_coords_patch ( ixx, iyy - 1, izz, 3)

          tempaa = tau ( 1) * n3 ( 1) + tau ( 2) * n3 ( 2)
     1           + tau ( 3) * n3 ( 3)

          tau = tau - tempaa * n3

          tempaa = 1.0 / SQRT ( tau ( 1)**2 + tau ( 2)**2 + tau ( 3)**2)

          tau ( :) = tau ( :) * tempaa

! Now we do "neta = tau X n3".

          neta ( 1) = tau ( 2) * n3 ( 3) - tau ( 3) * n3 ( 2)
          neta ( 2) = tau ( 3) * n3 ( 1) - tau ( 1) * n3 ( 3)
          neta ( 3) = tau ( 1) * n3 ( 2) - tau ( 2) * n3 ( 1)

! ----------
! Find the cosine and sine of the angle between the x-direction of the mesh
! and the y-direction of the mesh.
          cos_t = tau ( 1) * n1 ( 1) + tau ( 2) * n1 ( 2)
     1          + tau ( 3) * n1 ( 3)

          sin_t = SQRT ( 1.0 - cos_t**2)

! --------------------------------------------------
! Now that we have the properly scaled relative coordinates from each of the
! four particpating zones, we can use these to obtain the four input electric field components
! "Ex,y,z_RU, Ex,y,z_LU, Ex,y,z_LD, Ex,y,z_RD" that go into the 2D Riemann solver
! with 4 states.

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ == 0)
! --------------------
! Here the building of the electric fields requires that they have been reconstructed.

          Ex_LD = Ex_zc ( ixxLD, iyyLD, izzLD, 1)
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 2) * xLD
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 3) * yLD
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 4) * zLD
#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 5)
     1          * (xLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 5))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 6)
     1          * (yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 6))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 7)
     1          * (zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 7))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 8)
     1          * (xLD * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 8))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 9)
     1          * (yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 9))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 10)
     1          * (xLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 11)
     1          * (xLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 11))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 12)
     1          * (yLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 12))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 13)
     1          * (zLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 13))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 14)
     1          * (xLD**2 * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 14))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 15)
     1          * (xLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 15))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 16)
     1          * (xLD * yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 16))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 17)
     1          * (yLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 17))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 18)
     1          * (xLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 18))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 19)
     1          * (yLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 19))
     1          + Ex_zc ( ixxLD, iyyLD, izzLD, 20)
     1          * (xLD * yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 20))

#endif

          Ey_LD = Ey_zc ( ixxLD, iyyLD, izzLD, 1)
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 2) * xLD
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 3) * yLD
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 4) * zLD
#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 5)
     1          * (xLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 5))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 6)
     1          * (yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 6))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 7)
     1          * (zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 7))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 8)
     1          * (xLD * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 8))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 9)
     1          * (yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 9))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 10)
     1          * (xLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 11)
     1          * (xLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 11))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 12)
     1          * (yLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 12))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 13)
     1          * (zLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 13))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 14)
     1          * (xLD**2 * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 14))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 15)
     1          * (xLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 15))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 16)
     1          * (xLD * yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 16))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 17)
     1          * (yLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 17))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 18)
     1          * (xLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 18))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 19)
     1          * (yLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 19))
     1          + Ey_zc ( ixxLD, iyyLD, izzLD, 20)
     1          * (xLD * yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 20))

#endif

          Ez_LD = Ez_zc ( ixxLD, iyyLD, izzLD, 1)
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 2) * xLD
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 3) * yLD
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 4) * zLD
#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 5)
     1          * (xLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 5))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 6)
     1          * (yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 6))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 7)
     1          * (zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 7))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 8)
     1          * (xLD * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 8))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 9)
     1          * (yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 9))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 10)
     1          * (xLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 11)
     1          * (xLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 11))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 12)
     1          * (yLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 12))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 13)
     1          * (zLD**3
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 13))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 14)
     1          * (xLD**2 * yLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 14))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 15)
     1          * (xLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 15))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 16)
     1          * (xLD * yLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 16))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 17)
     1          * (yLD**2 * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 17))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 18)
     1          * (xLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 18))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 19)
     1          * (yLD * zLD**2
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 19))
     1          + Ez_zc ( ixxLD, iyyLD, izzLD, 20)
     1          * (xLD * yLD * zLD
     1          - Coef_Basis_Moments ( ixxLD, iyyLD, izzLD, 20))

#endif

! ----------

          Ex_RD = Ex_zc ( ixxRD, iyyRD, izzRD, 1)
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 2) * xRD
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 3) * yRD
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 4) * zRD
#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 5)
     1          * (xRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 5))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 6)
     1          * (yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 6))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 7)
     1          * (zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 7))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 8)
     1          * (xRD * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 8))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 9)
     1          * (yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 9))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 10)
     1          * (xRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 11)
     1          * (xRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 11))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 12)
     1          * (yRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 12))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 13)
     1          * (zRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 13))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 14)
     1          * (xRD**2 * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 14))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 15)
     1          * (xRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 15))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 16)
     1          * (xRD * yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 16))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 17)
     1          * (yRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 17))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 18)
     1          * (xRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 18))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 19)
     1          * (yRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 19))
     1          + Ex_zc ( ixxRD, iyyRD, izzRD, 20)
     1          * (xRD * yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 20))

#endif

          Ey_RD = Ey_zc ( ixxRD, iyyRD, izzRD, 1)
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 2) * xRD
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 3) * yRD
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 4) * zRD
#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 5)
     1          * (xRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 5))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 6)
     1          * (yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 6))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 7)
     1          * (zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 7))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 8)
     1          * (xRD * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 8))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 9)
     1          * (yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 9))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 10)
     1          * (xRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 11)
     1          * (xRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 11))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 12)
     1          * (yRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 12))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 13)
     1          * (zRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 13))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 14)
     1          * (xRD**2 * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 14))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 15)
     1          * (xRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 15))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 16)
     1          * (xRD * yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 16))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 17)
     1          * (yRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 17))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 18)
     1          * (xRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 18))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 19)
     1          * (yRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 19))
     1          + Ey_zc ( ixxRD, iyyRD, izzRD, 20)
     1          * (xRD * yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 20))

#endif

          Ez_RD = Ez_zc ( ixxRD, iyyRD, izzRD, 1)
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 2) * xRD
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 3) * yRD
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 4) * zRD
#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 5)
     1          * (xRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 5))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 6)
     1          * (yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 6))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 7)
     1          * (zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 7))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 8)
     1          * (xRD * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 8))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 9)
     1          * (yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 9))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 10)
     1          * (xRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 11)
     1          * (xRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 11))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 12)
     1          * (yRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 12))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 13)
     1          * (zRD**3
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 13))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 14)
     1          * (xRD**2 * yRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 14))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 15)
     1          * (xRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 15))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 16)
     1          * (xRD * yRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 16))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 17)
     1          * (yRD**2 * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 17))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 18)
     1          * (xRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 18))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 19)
     1          * (yRD * zRD**2
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 19))
     1          + Ez_zc ( ixxRD, iyyRD, izzRD, 20)
     1          * (xRD * yRD * zRD
     1          - Coef_Basis_Moments ( ixxRD, iyyRD, izzRD, 20))

#endif

! ----------

          Ex_RU = Ex_zc ( ixxRU, iyyRU, izzRU, 1)
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 2) * xRU
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 3) * yRU
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 4) * zRU
#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 5)
     1          * (xRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 5))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 6)
     1          * (yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 6))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 7)
     1          * (zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 7))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 8)
     1          * (xRU * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 8))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 9)
     1          * (yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 9))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 10)
     1          * (xRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 11)
     1          * (xRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 11))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 12)
     1          * (yRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 12))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 13)
     1          * (zRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 13))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 14)
     1          * (xRU**2 * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 14))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 15)
     1          * (xRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 15))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 16)
     1          * (xRU * yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 16))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 17)
     1          * (yRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 17))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 18)
     1          * (xRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 18))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 19)
     1          * (yRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 19))
     1          + Ex_zc ( ixxRU, iyyRU, izzRU, 20)
     1          * (xRU * yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 20))

#endif

          Ey_RU = Ey_zc ( ixxRU, iyyRU, izzRU, 1)
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 2) * xRU
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 3) * yRU
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 4) * zRU
#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 5)
     1          * (xRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 5))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 6)
     1          * (yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 6))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 7)
     1          * (zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 7))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 8)
     1          * (xRU * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 8))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 9)
     1          * (yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 9))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 10)
     1          * (xRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 11)
     1          * (xRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 11))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 12)
     1          * (yRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 12))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 13)
     1          * (zRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 13))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 14)
     1          * (xRU**2 * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 14))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 15)
     1          * (xRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 15))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 16)
     1          * (xRU * yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 16))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 17)
     1          * (yRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 17))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 18)
     1          * (xRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 18))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 19)
     1          * (yRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 19))
     1          + Ey_zc ( ixxRU, iyyRU, izzRU, 20)
     1          * (xRU * yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 20))

#endif

          Ez_RU = Ez_zc ( ixxRU, iyyRU, izzRU, 1)
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 2) * xRU
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 3) * yRU
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 4) * zRU
#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 5)
     1          * (xRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 5))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 6)
     1          * (yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 6))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 7)
     1          * (zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 7))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 8)
     1          * (xRU * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 8))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 9)
     1          * (yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 9))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 10)
     1          * (xRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 11)
     1          * (xRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 11))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 12)
     1          * (yRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 12))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 13)
     1          * (zRU**3
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 13))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 14)
     1          * (xRU**2 * yRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 14))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 15)
     1          * (xRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 15))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 16)
     1          * (xRU * yRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 16))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 17)
     1          * (yRU**2 * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 17))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 18)
     1          * (xRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 18))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 19)
     1          * (yRU * zRU**2
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 19))
     1          + Ez_zc ( ixxRU, iyyRU, izzRU, 20)
     1          * (xRU * yRU * zRU
     1          - Coef_Basis_Moments ( ixxRU, iyyRU, izzRU, 20))

#endif

! ----------

          Ex_LU = Ex_zc ( ixxLU, iyyLU, izzLU, 1)
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 2) * xLU
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 3) * yLU
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 4) * zLU
#if (NUMORDER >= 3)
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 5)
     1          * (xLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 5))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 6)
     1          * (yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 6))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 7)
     1          * (zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 7))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 8)
     1          * (xLU * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 8))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 9)
     1          * (yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 9))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 10)
     1          * (xLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 10))

#elif (NUMORDER >= 4)
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 11)
     1          * (xLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 11))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 12)
     1          * (yLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 12))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 13)
     1          * (zLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 13))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 14)
     1          * (xLU**2 * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 14))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 15)
     1          * (xLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 15))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 16)
     1          * (xLU * yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 16))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 17)
     1          * (yLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 17))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 18)
     1          * (xLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 18))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 19)
     1          * (yLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 19))
     1          + Ex_zc ( ixxLU, iyyLU, izzLU, 20)
     1          * (xLU * yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 20))

#endif

          Ey_LU = Ey_zc ( ixxLU, iyyLU, izzLU, 1)
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 2) * xLU
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 3) * yLU
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 4) * zLU
#if (NUMORDER >= 3)
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 5)
     1          * (xLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 5))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 6)
     1          * (yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 6))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 7)
     1          * (zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 7))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 8)
     1          * (xLU * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 8))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 9)
     1          * (yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 9))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 10)
     1          * (xLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 10))

#elif (NUMORDER >= 4)
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 11)
     1          * (xLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 11))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 12)
     1          * (yLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 12))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 13)
     1          * (zLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 13))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 14)
     1          * (xLU**2 * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 14))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 15)
     1          * (xLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 15))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 16)
     1          * (xLU * yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 16))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 17)
     1          * (yLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 17))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 18)
     1          * (xLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 18))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 19)
     1          * (yLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 19))
     1          + Ey_zc ( ixxLU, iyyLU, izzLU, 20)
     1          * (xLU * yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 20))

#endif

          Ez_LU = Ez_zc ( ixxLU, iyyLU, izzLU, 1)
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 2) * xLU
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 3) * yLU
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 4) * zLU
#if (NUMORDER >= 3)
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 5)
     1          * (xLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 5))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 6)
     1          * (yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 6))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 7)
     1          * (zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 7))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 8)
     1          * (xLU * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 8))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 9)
     1          * (yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 9))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 10)
     1          * (xLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 10))

#elif (NUMORDER >= 4)
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 11)
     1          * (xLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 11))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 12)
     1          * (yLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 12))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 13)
     1          * (zLU**3
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 13))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 14)
     1          * (xLU**2 * yLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 14))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 15)
     1          * (xLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 15))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 16)
     1          * (xLU * yLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 16))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 17)
     1          * (yLU**2 * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 17))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 18)
     1          * (xLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 18))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 19)
     1          * (yLU * zLU**2
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 19))
     1          + Ez_zc ( ixxLU, iyyLU, izzLU, 20)
     1          * (xLU * yLU * zLU
     1          - Coef_Basis_Moments ( ixxLU, iyyLU, izzLU, 20))

#endif

! --------------------
#else
! --------------------
! Choose arithmetically averaged electric fields for he centered part of the
! 2D LLF Riemann solver. Electric field reconstruction is not needed.

          Ex_LD = ( Ex_zc ( ixxLD, iyyLD, izzLD, 1)
     1            + Ex_zc ( ixxRD, iyyRD, izzRD, 1)
     1            + Ex_zc ( ixxRU, iyyRU, izzRU, 1)
     1            + Ex_zc ( ixxLU, iyyLU, izzLU, 1) ) * 0.25
          Ex_RD = Ex_LD  ;  Ex_RU = Ex_LD  ;  Ex_LU = Ex_LD

          Ey_LD = ( Ey_zc ( ixxLD, iyyLD, izzLD, 1)
     1            + Ey_zc ( ixxRD, iyyRD, izzRD, 1)
     1            + Ey_zc ( ixxRU, iyyRU, izzRU, 1)
     1            + Ey_zc ( ixxLU, iyyLU, izzLU, 1) ) * 0.25
          Ey_RD = Ey_LD  ;  Ey_RU = Ey_LD  ;  Ey_LU = Ey_LD

          Ez_LD = ( Ez_zc ( ixxLD, iyyLD, izzLD, 1)
     1            + Ez_zc ( ixxRD, iyyRD, izzRD, 1)
     1            + Ez_zc ( ixxRU, iyyRU, izzRU, 1)
     1            + Ez_zc ( ixxLU, iyyLU, izzLU, 1) ) * 0.25
          Ez_RD = Ez_LD  ;  Ez_RU = Ez_LD  ;  Ez_LU = Ez_LD

! --------------------
#endif
! --------------------

! ------------------------------
! Build the electric fields in the local z-direction.

          tempaa = Ex_LD * n3 ( 1) + Ey_LD * n3 ( 2) + Ez_LD * n3 ( 3)
          Ez_LD = tempaa

          tempaa = Ex_RD * n3 ( 1) + Ey_RD * n3 ( 2) + Ez_RD * n3 ( 3)
          Ez_RD = tempaa

          tempaa = Ex_RU * n3 ( 1) + Ey_RU * n3 ( 2) + Ez_RU * n3 ( 3)
          Ez_RU = tempaa

          tempaa = Ex_LU * n3 ( 1) + Ey_LU * n3 ( 2) + Ez_LU * n3 ( 3)
          Ez_LU = tempaa

! ------------------------------
! Build "By_R, By_L" from the RD and LD y-faces respectively.
! These have to include the "Bx0, By0, Bz0" contributions in the "n2" direction.

          By_R = Byf_quadpts ( ixxRD, iyyRD, izzRD, 5)
     1        + Bx0 * n2 ( 1) + By0 * n2 ( 2) + Bz0 * n2 ( 3)

          By_L = Byf_quadpts ( ixxLD, iyyLD, izzLD, 4)
     1        + Bx0 * n2 ( 1) + By0 * n2 ( 2) + Bz0 * n2 ( 3)

! ------------------------------
! Build "Bneta_U, Bneta_D" from the LU and LD faces that are in the x-direction
! of the mesh.
! These have to include the "Bx0, By0, Bz0" contributions in the "neta" direction.

          Bneta_U = Bxf_quadpts ( ixxLU, iyyLU, izzLU, 3)
     1            + Bx0 * neta ( 1) + By0 * neta ( 2) + Bz0 * neta ( 3)

          Bneta_D = Bxf_quadpts ( ixxLD, iyyLD, izzLD, 2)
     1            + Bx0 * neta ( 1) + By0 * neta ( 2) + Bz0 * neta ( 3)

! ------------------------------
! Evaluate "S_R, S_L, S_U, S_D". The Right and Left waves from the x-faces of the mesh
! give us the "S_R, S_L". The Up and Down Left waves from the y-faces of the mesh
! give us "S_U, S_D".

          S_R = AMAX1 ( S_R_gr ( ixxLU, iyyLU, izzLU),
     1                  S_R_gr ( ixxLD, iyyLD, izzLD) )

          S_L = AMIN1 ( S_L_gr ( ixxLU, iyyLU, izzLU),
     1                  S_L_gr ( ixxLD, iyyLD, izzLD) )

          S_U = AMAX1 ( S_U_gr ( ixxRD, iyyRD, izzRD),
     1                  S_U_gr ( ixxLD, iyyLD, izzLD) )

          S_D = AMIN1 ( S_D_gr ( ixxRD, iyyRD, izzRD),
     1                  S_D_gr ( ixxLD, iyyLD, izzLD) )

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ != 0)
! --------------------
! Make wave speeds fully subsonic for the 2D RS.

          tempaa = AMAX1 ( ABS ( S_R), ABS ( S_L),
     1                     ABS ( S_U), ABS ( S_D) )

          S_R = tempaa  ;  S_L = - tempaa
          S_U = tempaa  ;  S_D = - tempaa

! --------------------
#endif
! --------------------

! ------------------------------------------------------------
! We are only making z-directional calls to the 2D Riemann solver that has four
! input states.

          CALL Two_D_RS_Ez_4State_HLL

     1         ( S_R, S_L, S_U, S_D, cos_t, sin_t,

     1           Ez_RU, Ez_LU, Ez_LD, Ez_RD,

     1           Bneta_U, Bneta_D, By_R, By_L,

     1           BSSx, BSSy, ESSz)

! ------------------------------------------------------------
! Hand back the global electric field "ESSz" from the 2D Riemann solver.

          Efz ( ixx, iyy, izz) = Efz ( ixx, iyy, izz) + wt * ESSz

! ------------------------------------------------------------

          END DO  ! End of "DO izz = 1, nr_zones" loop.

        END DO  ! End of "DO iyy = 0, ntp_zones" loop.

      END DO  ! End of "DO ixx = 0, ntp_zones" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, 3" loop over quadrature points.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 4) Now we have between 0 to four 1D loops in the r-direction which are only
! invoked if we are at any of the four tri-corners.

! ----------------------------------------------------------------------
! STEP 4.a) The Left-Down corner of this mesh is a tri-corner.

      IF ( itri_patch_LD .EQ. 1) THEN

      ixx = 0  ! These are the indices of the z-edge along which we want the e-field.
      iyy = 0

! At third and fourth order, loop over 3 quadrature points in the edge of interest.
      DO j = 1, 3

#if (NUMORDER == 2)
! At second order, we need only one nodal point.
      IF ( j .GT. 1) EXIT
#endif

      wt = 1.0
#if (NUMORDER >= 3)
      wt = 1.0 / 6.0  ! These are the weights for simpson integration.
      IF ( j .EQ. 1) wt = 4.0 / 6.0
#endif

!$OMP PARALLEL DO PRIVATE (
!$OMP&    izz, ixx1, iyy1, ixx2, iyy2, ixx3, iyy3,                      &
!$OMP&    xec, yec, zec, tempaa, obl_1, obl_2, obl_3,                   &
!$OMP&    n1 , n2 , n3, B0, time_use, Bx0, By0, Bz0,                    &
!$OMP&    tau_1, neta_1, tau_2, neta_2, tau_3, neta_3,                  &

!$OMP&    x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3,                  &
!$OMP&    Ex_1, Ey_1, Ez_1, Ex_2, Ey_2, Ez_2, Ex_3, Ey_3, Ez_3,         &
!$OMP&    S_R1, S_L1, S_R2, S_L2, S_R3, S_L3, Bneta_1, Bneta_2, Bneta_3,&
!$OMP&    BSSx, BSSy, ESSz)                                             &

!$OMP&             SHARED ( nr_zones, ixx, iyy, smallnum, ob3, j, wt,   &
!$OMP&    vert_coords_patch, v_scale_length, v_centroid_patch,          &
!$OMP&    Coef_Basis_Moments,                                           &

!$OMP&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$OMP&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$OMP&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$OMP&    Efz)                                                          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum, ob3, j, wt,                                         &
!$acc&    vert_coords_patch, v_scale_length, v_centroid_patch,          &
!$acc&    Coef_Basis_Moments,                                           &

!$acc&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$acc&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$acc&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$acc&    Efz)

!$acc loop vector independent                                           &
!$acc&    private ( izz, ixx1, iyy1, ixx2, iyy2, ixx3, iyy3,            &
!$acc&    xec, yec, zec, tempaa, obl_1, obl_2, obl_3,                   &
!$acc&    n1 , n2 , n3, B0, time_use, Bx0, By0, Bz0,                    &
!$acc&    tau_1, neta_1, tau_2, neta_2, tau_3, neta_3,                  &

!$acc&    x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3,                  &
!$acc&    Ex_1, Ey_1, Ez_1, Ex_2, Ey_2, Ez_2, Ex_3, Ey_3, Ez_3,         &
!$acc&    S_R1, S_L1, S_R2, S_L2, S_R3, S_L3, Bneta_1, Bneta_2, Bneta_3,&
!$acc&    BSSx, BSSy, ESSz)


        DO izz = 1, nr_zones

! ------------------------------------------------------------

! Fill in "ixx1, iyy1, ixx2, iyy2, ixx3, iyy3" in counter-clockwise sense.

        ixx1 = ixx + 1  ;  iyy1 = iyy + 1
        ixx2 = ixx      ;  iyy2 = iyy + 1
        ixx3 = ixx + 1  ;  iyy3 = iyy

! --------------------------------------------------
! Find the edge-center "xec, yec, zec" by bisecting the two z-vertices.
! No need to rescale by radius since this is in the radial direction by definition.

          IF ( j .EQ. 1) THEN  ! First quadrature point is the center point.

            xec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 1)  ! Bisect the z-vertices.
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 1) )
            yec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 2)
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 2) )
            zec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 3)
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 3) )

          ELSE IF ( j .EQ. 2) THEN  ! Second quadrature point is one of the end points.

            xec = vert_coords_patch ( ixx, iyy, izz, 1)
            yec = vert_coords_patch ( ixx, iyy, izz, 2)
            zec = vert_coords_patch ( ixx, iyy, izz, 3)

          ELSE IF ( j .EQ. 3) THEN  ! Third quadrature point is the other end point.

            xec = vert_coords_patch ( ixx, iyy, izz - 1, 1)
            yec = vert_coords_patch ( ixx, iyy, izz - 1, 2)
            zec = vert_coords_patch ( ixx, iyy, izz - 1, 3)

          END IF

! ------------------------------
! Feed ( xec, yec, zec) to the tiny subroutine that evaluates "B0".
! Default is zero.

          B0 = 0.0
          Bx0 = 0.0
          By0 = 0.0
          Bz0 = 0.0

#if ( USE_B0B1_FORMULATION == 1)

          bflag = 1
          time_use = 0

          CALL GET_MAGNETIC_FIELD ( bflag, time_use,  xec, yec, zec,
     1                              Bx0, By0, Bz0)

          B0 ( 1) = Bx0
          B0 ( 2) = By0
          B0 ( 3) = Bz0

#endif

! ------------------------------

        obl_1 = 1.0 / AMAX1 ( v_scale_length ( ixx1, iyy1, izz),
     1                        smallnum)

        obl_2 = 1.0 / AMAX1 ( v_scale_length ( ixx2, iyy2, izz),
     1                        smallnum)

        obl_3 = 1.0 / AMAX1 ( v_scale_length ( ixx3, iyy3, izz),
     1                        smallnum)

! ------------------------------
! Build scaled "x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3".
! These are needed in general if we are using the reconstruction. But can
! be bypassed at second order.

        x_1 = ( xec - v_centroid_patch ( ixx1, iyy1, izz, 1) )
     1        * obl_1

        y_1 = ( yec - v_centroid_patch ( ixx1, iyy1, izz, 2) )
     1        * obl_1

        z_1 = ( zec - v_centroid_patch ( ixx1, iyy1, izz, 3) )
     1        * obl_1

! ----------

        x_2 = ( xec - v_centroid_patch ( ixx2, iyy2, izz, 1) )
     1        * obl_2

        y_2 = ( yec - v_centroid_patch ( ixx2, iyy2, izz, 2) )
     1        * obl_2

        z_2 = ( zec - v_centroid_patch ( ixx2, iyy2, izz, 3) )
     1        * obl_2

! ----------

        x_3 = ( xec - v_centroid_patch ( ixx3, iyy3, izz, 1) )
     1        * obl_3

        y_3 = ( yec - v_centroid_patch ( ixx3, iyy3, izz, 2) )
     1        * obl_3

        z_3 = ( zec - v_centroid_patch ( ixx3, iyy3, izz, 3) )
     1        * obl_3

! ------------------------------
! Since this is a tri-corner, we have to build:- "n3" along radial direction.
! "tau_1" along x-axis of mesh and "neta_1 = tau_1 X n3" perpendicular to it by using "n3".
! "tau_2" along y-axis of mesh and "neta_2 = tau_2 X n3" perpendicular to it by using "n3".
! "tau_3" along -ve x-axis of mesh and "neta_3 = tau_3 X n3" perpendicular to it by using "n3".

          tempaa = 1.0 / SQRT ( xec**2 + yec**2 + zec**2)

! Here "n3" is in the global r-direction of the mesh, which is our local z-direction.
          n3 ( 1) = xec * tempaa   ! This is in the global r-direction.
          n3 ( 2) = yec * tempaa
          n3 ( 3) = zec * tempaa

! ----------

          tau_1 ( 1) = vert_coords_patch ( ixx + 1, iyy, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_1 ( 2) = vert_coords_patch ( ixx + 1, iyy, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_1 ( 3) = vert_coords_patch ( ixx + 1, iyy, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_1 ( 1) * n3 ( 1) + tau_1 ( 2) * n3 ( 2)
     1           + tau_1 ( 3) * n3 ( 3)

          tau_1 = tau_1 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_1 ( 1)**2 + tau_1 ( 2)**2
     1                        + tau_1 ( 3)**2)

          tau_1 ( :) = tau_1 ( :) * tempaa

          neta_1 ( 1) = tau_1 ( 2) * n3 ( 3) - tau_1 ( 3) * n3 ( 2)
          neta_1 ( 2) = tau_1 ( 3) * n3 ( 1) - tau_1 ( 1) * n3 ( 3)
          neta_1 ( 3) = tau_1 ( 1) * n3 ( 2) - tau_1 ( 2) * n3 ( 1)

! ----------

          tau_2 ( 1) = vert_coords_patch ( ixx, iyy + 1, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_2 ( 2) = vert_coords_patch ( ixx, iyy + 1, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_2 ( 3) = vert_coords_patch ( ixx, iyy + 1, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_2 ( 1) * n3 ( 1) + tau_2 ( 2) * n3 ( 2)
     1           + tau_2 ( 3) * n3 ( 3)

          tau_2 = tau_2 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_2 ( 1)**2 + tau_2 ( 2)**2
     1                        + tau_2 ( 3)**2)

          tau_2 ( :) = tau_2 ( :) * tempaa

          neta_2 ( 1) = tau_2 ( 2) * n3 ( 3) - tau_2 ( 3) * n3 ( 2)
          neta_2 ( 2) = tau_2 ( 3) * n3 ( 1) - tau_2 ( 1) * n3 ( 3)
          neta_2 ( 3) = tau_2 ( 1) * n3 ( 2) - tau_2 ( 2) * n3 ( 1)

! ----------

          tau_3 ( 1) = vert_coords_patch ( ixx - 1, iyy, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_3 ( 2) = vert_coords_patch ( ixx - 1, iyy, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_3 ( 3) = vert_coords_patch ( ixx - 1, iyy, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_3 ( 1) * n3 ( 1) + tau_3 ( 2) * n3 ( 2)
     1           + tau_3 ( 3) * n3 ( 3)

          tau_3 = tau_3 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_3 ( 1)**2 + tau_3 ( 2)**2
     1                        + tau_3 ( 3)**2)

          tau_3 ( :) = tau_3 ( :) * tempaa

          neta_3 ( 1) = tau_3 ( 2) * n3 ( 3) - tau_3 ( 3) * n3 ( 2)
          neta_3 ( 2) = tau_3 ( 3) * n3 ( 1) - tau_3 ( 1) * n3 ( 3)
          neta_3 ( 3) = tau_3 ( 1) * n3 ( 2) - tau_3 ( 2) * n3 ( 1)

! --------------------------------------------------
! Now that we have the properly scaled relative coordinates from each of the
! three particpating zones, we can use these to obtain the three input electric field
! components "Ez_1, Ez_2, Ez_3" that go into the 2D Riemann solver
! with 3 states.

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ == 0)
! --------------------
! Here the building of the electric fields requires that they have been reconstructed.

          Ex_1 = Ex_zc ( ixx1, iyy1, izz, 1)
     1         + Ex_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ex_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ex_zc ( ixx1, iyy1, izz, 4) * z_1

#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ex_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ex_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ex_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ex_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ex_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ex_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ex_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ex_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ex_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ex_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ex_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ex_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ex_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ex_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

          Ey_1 = Ey_zc ( ixx1, iyy1, izz, 1)
     1         + Ey_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ey_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ey_zc ( ixx1, iyy1, izz, 4) * z_1

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ey_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ey_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ey_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ey_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ey_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ey_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ey_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ey_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ey_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ey_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ey_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ey_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ey_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ey_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

          Ez_1 = Ez_zc ( ixx1, iyy1, izz, 1)
     1         + Ez_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ez_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ez_zc ( ixx1, iyy1, izz, 4) * z_1

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ez_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ez_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ez_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ez_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ez_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ez_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ez_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ez_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ez_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ez_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ez_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ez_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ez_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ez_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

! ----------

          Ex_2 = Ex_zc ( ixx2, iyy2, izz, 1)
     1         + Ex_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ex_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ex_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ex_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ex_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ex_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ex_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ex_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ex_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ex_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ex_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ex_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ex_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ex_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ex_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ex_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ex_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

          Ey_2 = Ey_zc ( ixx2, iyy2, izz, 1)
     1         + Ey_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ey_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ey_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ey_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ey_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ey_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ey_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ey_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ey_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ey_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ey_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ey_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ey_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ey_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ey_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ey_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ey_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

          Ez_2 = Ez_zc ( ixx2, iyy2, izz, 1)
     1         + Ez_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ez_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ez_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ez_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ez_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ez_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ez_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ez_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ez_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ez_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ez_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ez_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ez_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ez_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ez_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ez_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ez_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

! ----------

          Ex_3 = Ex_zc ( ixx3, iyy3, izz, 1)
     1         + Ex_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ex_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ex_zc ( ixx3, iyy3, izz, 4) * z_3
#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ex_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ex_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ex_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ex_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ex_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ex_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ex_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ex_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ex_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ex_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ex_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ex_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ex_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ex_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

          Ey_3 = Ey_zc ( ixx3, iyy3, izz, 1)
     1         + Ey_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ey_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ey_zc ( ixx3, iyy3, izz, 4) * z_3

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ey_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ey_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ey_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ey_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ey_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ey_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ey_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ey_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ey_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ey_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ey_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ey_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ey_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ey_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

          Ez_3 = Ez_zc ( ixx3, iyy3, izz, 1)
     1         + Ez_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ez_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ez_zc ( ixx3, iyy3, izz, 4) * z_3
#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ez_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ez_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ez_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ez_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ez_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ez_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ez_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ez_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ez_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ez_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ez_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ez_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ez_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ez_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

! --------------------
#else
! --------------------
! Choose arithmetically averaged electric fields for the centered part of the
! 2D LLF Riemann solver. Electric field reconstruction is not needed.

          Ex_1 = ( Ex_zc ( ixx1, iyy1, izz, 1)
     1           + Ex_zc ( ixx2, iyy2, izz, 1)
     1           + Ex_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ex_2 = Ex_1  ;  Ex_3 = Ex_1

          Ey_1 = ( Ey_zc ( ixx1, iyy1, izz, 1)
     1           + Ey_zc ( ixx2, iyy2, izz, 1)
     1           + Ey_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ey_2 = Ey_1  ;  Ey_3 = Ey_1

          Ez_1 = ( Ez_zc ( ixx1, iyy1, izz, 1)
     1           + Ez_zc ( ixx2, iyy2, izz, 1)
     1           + Ez_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ez_2 = Ez_1  ;  Ez_3 = Ez_1

! --------------------
#endif
! --------------------

! ------------------------------
! Build the electric fields in the local z-direction; which is the global r-direction.

          tempaa = Ex_1 * n3 ( 1) + Ey_1 * n3 ( 2) + Ez_1 * n3 ( 3)
          Ez_1 = tempaa

          tempaa = Ex_2 * n3 ( 1) + Ey_2 * n3 ( 2) + Ez_2 * n3 ( 3)
          Ez_2 = tempaa

          tempaa = Ex_3 * n3 ( 1) + Ey_3 * n3 ( 2) + Ez_3 * n3 ( 3)
          Ez_3 = tempaa

! ------------------------------
! Build "Bneta_1, Bneta_2, Bneta_3" They may not be oriented in the positive directions.
! "Bneta_1" is along "neta_1"; "Bneta_2" is along "neta_2"; "Bneta_3" is along "neta_3".

          Bneta_1 = - Byf_quadpts ( ixx3, iyy3, izz, 5)
     1  + Bx0 * neta_1 ( 1) + By0 * neta_1 ( 2) + Bz0 * neta_1 ( 3)

          Bneta_2 = Bxf_quadpts ( ixx2, iyy2, izz, 3)
     1  + Bx0 * neta_2 ( 1) + By0 * neta_2 ( 2) + Bz0 * neta_2 ( 3)

          Bneta_3 = Byf_quadpts ( ixx2, iyy2 - 1, izz, 4)
     1  + Bx0 * neta_3 ( 1) + By0 * neta_3 ( 2) + Bz0 * neta_3 ( 3)

! ------------------------------
! Evaluate "S_R1, S_L1, S_R2, S_L2, S_R3, S_L3". Watch out for the orientatations
! because in each of the three zone faces we take the positive direction to be
! along "neta_1, neta_2, neta_3".

          S_R1 = - S_D_gr ( ixx3, iyy3, izz)
          S_L1 = - S_U_gr ( ixx3, iyy3, izz)

          S_R2 = S_R_gr ( ixx2, iyy2, izz)
          S_L2 = S_L_gr ( ixx2, iyy2, izz)

          S_R3 = S_U_gr ( ixx2, iyy2 - 1, izz)
          S_L3 = S_D_gr ( ixx2, iyy2 - 1, izz)

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ != 0)
! --------------------
! Make wave speeds fully subsonic for the 2D RS.

          tempaa = AMAX1 ( ABS ( S_R1), ABS ( S_L1),
     1                     ABS ( S_R2), ABS ( S_L2),
     1                     ABS ( S_R3), ABS ( S_L3))

          S_R1 = tempaa  ;  S_L1 = - tempaa
          S_R2 = tempaa  ;  S_L2 = - tempaa
          S_R3 = tempaa  ;  S_L3 = - tempaa

! --------------------
#endif
! --------------------

! ------------------------------------------------------------
! We are only making z-directional calls to the 2D Riemann solver that has three
! input states.

          CALL Two_D_RS_Ez_3State_HLL

     1         ( S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,

     1           Ez_1, Ez_2, Ez_3,

     1           Bneta_1, Bneta_2, Bneta_3,

     1           BSSx, BSSy, ESSz)

! ------------------------------------------------------------
! Hand back the global electric field "ESSz" from the 2D Riemann solver.

          Efz ( ixx, iyy, izz) = Efz ( ixx, iyy, izz) + wt * ESSz

! ------------------------------------------------------------

        END DO  ! End of "DO izz = 1, nr_zones".

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, 3" loop over quadrature points.

      END IF  ! End of "IF ( itri_patch_LD .EQ. 1) THEN".

! ----------------------------------------------------------------------
! STEP 4.b) The Right-Down corner of this mesh is a tri-corner.

      IF ( itri_patch_RD .EQ. 1) THEN

      ixx = ntp_zones  ! These are the indices of the z-edge along which we want the e-field.
      iyy = 0

! At third and fourth order, loop over 3 quadrature points in the edge of interest.
      DO j = 1, 3

#if (NUMORDER == 2)
! At second order, we need only one nodal point.
      IF ( j .GT. 1) EXIT
#endif

      wt = 1.0
#if (NUMORDER >= 3)
      wt = 1.0 / 6.0  ! These are the weights for simpson integration.
      IF ( j .EQ. 1) wt = 4.0 / 6.0
#endif


!$OMP PARALLEL DO PRIVATE (
!$OMP&    izz, ixx1, iyy1, ixx2, iyy2, ixx3, iyy3,                      &
!$OMP&    xec, yec, zec, tempaa, obl_1, obl_2, obl_3,                   &
!$OMP&    n1 , n2 , n3, B0, time_use, Bx0, By0, Bz0,                    &
!$OMP&    tau_1, neta_1, tau_2, neta_2, tau_3, neta_3,                  &

!$OMP&    x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3,                  &
!$OMP&    Ex_1, Ey_1, Ez_1, Ex_2, Ey_2, Ez_2, Ex_3, Ey_3, Ez_3,         &
!$OMP&    S_R1, S_L1, S_R2, S_L2, S_R3, S_L3, Bneta_1, Bneta_2, Bneta_3,&
!$OMP&    BSSx, BSSy, ESSz)                                             &

!$OMP&             SHARED ( nr_zones, ixx, iyy, smallnum, ob3, j, wt,   &
!$OMP&    vert_coords_patch, v_scale_length, v_centroid_patch,          &
!$OMP&    Coef_Basis_Moments,                                           &

!$OMP&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$OMP&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$OMP&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$OMP&    Efz)                                                          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum, ob3, j, wt,                                         &
!$acc&    vert_coords_patch, v_scale_length, v_centroid_patch,          &
!$acc&    Coef_Basis_Moments,                                           &

!$acc&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$acc&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$acc&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$acc&    Efz)

!$acc loop vector independent                                           &
!$acc&    private ( izz, ixx1, iyy1, ixx2, iyy2, ixx3, iyy3,            &
!$acc&    xec, yec, zec, tempaa, obl_1, obl_2, obl_3,                   &
!$acc&    n1 , n2 , n3, B0, time_use, Bx0, By0, Bz0,                    &
!$acc&    tau_1, neta_1, tau_2, neta_2, tau_3, neta_3,                  &

!$acc&    x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3,                  &
!$acc&    Ex_1, Ey_1, Ez_1, Ex_2, Ey_2, Ez_2, Ex_3, Ey_3, Ez_3,         &
!$acc&    S_R1, S_L1, S_R2, S_L2, S_R3, S_L3, Bneta_1, Bneta_2, Bneta_3,&
!$acc&    BSSx, BSSy, ESSz)

        DO izz = 1, nr_zones

! ------------------------------------------------------------

! Fill in "ixx1, iyy1, ixx2, iyy2, ixx3, iyy3" in counter-clockwise sense.

        ixx1 = ixx + 1  ;  iyy1 = iyy + 1
        ixx2 = ixx      ;  iyy2 = iyy + 1
        ixx3 = ixx      ;  iyy3 = iyy

! --------------------------------------------------
! Find the edge-center "xec, yec, zec" by bisecting the two z-vertices.
! No need to rescale by radius since this is in the radial direction by definition.

          IF ( j .EQ. 1) THEN  ! First quadrature point is the center point.

            xec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 1)  ! Bisect the z-vertices.
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 1) )
            yec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 2)
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 2) )
            zec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 3)
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 3) )

          ELSE IF ( j .EQ. 2) THEN  ! Second quadrature point is one of the end points.

            xec = vert_coords_patch ( ixx, iyy, izz, 1)
            yec = vert_coords_patch ( ixx, iyy, izz, 2)
            zec = vert_coords_patch ( ixx, iyy, izz, 3)

          ELSE IF ( j .EQ. 3) THEN  ! Third quadrature point is the other end point.

            xec = vert_coords_patch ( ixx, iyy, izz - 1, 1)
            yec = vert_coords_patch ( ixx, iyy, izz - 1, 2)
            zec = vert_coords_patch ( ixx, iyy, izz - 1, 3)

          END IF

! ------------------------------
! Feed ( xec, yec, zec) to the tiny subroutine that evaluates "B0".
! Default is zero.

          B0 = 0.0
          Bx0 = 0.0
          By0 = 0.0
          Bz0 = 0.0

#if ( USE_B0B1_FORMULATION == 1)

          bflag = 1
          time_use = 0

          CALL GET_MAGNETIC_FIELD ( bflag, time_use,  xec, yec, zec,
     1                              Bx0, By0, Bz0)

          B0 ( 1) = Bx0
          B0 ( 2) = By0
          B0 ( 3) = Bz0

#endif

! ------------------------------

        obl_1 = 1.0 / AMAX1 ( v_scale_length ( ixx1, iyy1, izz),
     1                        smallnum)

        obl_2 = 1.0 / AMAX1 ( v_scale_length ( ixx2, iyy2, izz),
     1                        smallnum)

        obl_3 = 1.0 / AMAX1 ( v_scale_length ( ixx3, iyy3, izz),
     1                        smallnum)

! ------------------------------
! Build scaled "x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3".
! These are needed in general if we are using the reconstruction. But can
! be bypassed at second order.

        x_1 = ( xec - v_centroid_patch ( ixx1, iyy1, izz, 1) )
     1        * obl_1

        y_1 = ( yec - v_centroid_patch ( ixx1, iyy1, izz, 2) )
     1        * obl_1

        z_1 = ( zec - v_centroid_patch ( ixx1, iyy1, izz, 3) )
     1        * obl_1

! ----------

        x_2 = ( xec - v_centroid_patch ( ixx2, iyy2, izz, 1) )
     1        * obl_2

        y_2 = ( yec - v_centroid_patch ( ixx2, iyy2, izz, 2) )
     1        * obl_2

        z_2 = ( zec - v_centroid_patch ( ixx2, iyy2, izz, 3) )
     1        * obl_2

! ----------

        x_3 = ( xec - v_centroid_patch ( ixx3, iyy3, izz, 1) )
     1        * obl_3

        y_3 = ( yec - v_centroid_patch ( ixx3, iyy3, izz, 2) )
     1        * obl_3

        z_3 = ( zec - v_centroid_patch ( ixx3, iyy3, izz, 3) )
     1        * obl_3

! ------------------------------
! Since this is a tri-corner, we have to build:- "n3" along radial direction.
! "tau_1" along x-axis of mesh and "neta_1 = tau_1 X n3" perpendicular to it by using "n3".
! "tau_2" along y-axis of mesh and "neta_2 = tau_2 X n3" perpendicular to it by using "n3".
! "tau_3" along -ve x-axis of mesh and "neta_3 = tau_3 X n3" perpendicular to it by using "n3".

          tempaa = 1.0 / SQRT ( xec**2 + yec**2 + zec**2)

! Here "n3" is in the global r-direction of the mesh, which is our local z-direction.
          n3 ( 1) = xec * tempaa   ! This is in the global r-direction.
          n3 ( 2) = yec * tempaa
          n3 ( 3) = zec * tempaa

! ----------

          tau_1 ( 1) = vert_coords_patch ( ixx + 1, iyy, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_1 ( 2) = vert_coords_patch ( ixx + 1, iyy, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_1 ( 3) = vert_coords_patch ( ixx + 1, iyy, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_1 ( 1) * n3 ( 1) + tau_1 ( 2) * n3 ( 2)
     1           + tau_1 ( 3) * n3 ( 3)

          tau_1 = tau_1 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_1 ( 1)**2 + tau_1 ( 2)**2
     1                        + tau_1 ( 3)**2)

          tau_1 ( :) = tau_1 ( :) * tempaa

          neta_1 ( 1) = tau_1 ( 2) * n3 ( 3) - tau_1 ( 3) * n3 ( 2)
          neta_1 ( 2) = tau_1 ( 3) * n3 ( 1) - tau_1 ( 1) * n3 ( 3)
          neta_1 ( 3) = tau_1 ( 1) * n3 ( 2) - tau_1 ( 2) * n3 ( 1)

! ----------

          tau_2 ( 1) = vert_coords_patch ( ixx, iyy + 1, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_2 ( 2) = vert_coords_patch ( ixx, iyy + 1, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_2 ( 3) = vert_coords_patch ( ixx, iyy + 1, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_2 ( 1) * n3 ( 1) + tau_2 ( 2) * n3 ( 2)
     1           + tau_2 ( 3) * n3 ( 3)

          tau_2 = tau_2 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_2 ( 1)**2 + tau_2 ( 2)**2
     1                        + tau_2 ( 3)**2)

          tau_2 ( :) = tau_2 ( :) * tempaa

          neta_2 ( 1) = tau_2 ( 2) * n3 ( 3) - tau_2 ( 3) * n3 ( 2)
          neta_2 ( 2) = tau_2 ( 3) * n3 ( 1) - tau_2 ( 1) * n3 ( 3)
          neta_2 ( 3) = tau_2 ( 1) * n3 ( 2) - tau_2 ( 2) * n3 ( 1)

! ----------

          tau_3 ( 1) = vert_coords_patch ( ixx - 1, iyy, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_3 ( 2) = vert_coords_patch ( ixx - 1, iyy, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_3 ( 3) = vert_coords_patch ( ixx - 1, iyy, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_3 ( 1) * n3 ( 1) + tau_3 ( 2) * n3 ( 2)
     1           + tau_3 ( 3) * n3 ( 3)

          tau_3 = tau_3 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_3 ( 1)**2 + tau_3 ( 2)**2
     1                        + tau_3 ( 3)**2)

          tau_3 ( :) = tau_3 ( :) * tempaa

          neta_3 ( 1) = tau_3 ( 2) * n3 ( 3) - tau_3 ( 3) * n3 ( 2)
          neta_3 ( 2) = tau_3 ( 3) * n3 ( 1) - tau_3 ( 1) * n3 ( 3)
          neta_3 ( 3) = tau_3 ( 1) * n3 ( 2) - tau_3 ( 2) * n3 ( 1)

! --------------------------------------------------
! Now that we have the properly scaled relative coordinates from each of the
! three particpating zones, we can use these to obtain the three input electric field
! components "Ez_1, Ez_2, Ez_3" that go into the 2D Riemann solver
! with 3 states.

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ == 0)
! --------------------
! Here the building of the electric fields requires that they have been reconstructed.

          Ex_1 = Ex_zc ( ixx1, iyy1, izz, 1)
     1         + Ex_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ex_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ex_zc ( ixx1, iyy1, izz, 4) * z_1
#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ex_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ex_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ex_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ex_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ex_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ex_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ex_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ex_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ex_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ex_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ex_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ex_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ex_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ex_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

          Ey_1 = Ey_zc ( ixx1, iyy1, izz, 1)
     1         + Ey_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ey_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ey_zc ( ixx1, iyy1, izz, 4) * z_1

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ey_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ey_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ey_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ey_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ey_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ey_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ey_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ey_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ey_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ey_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ey_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ey_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ey_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ey_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

          Ez_1 = Ez_zc ( ixx1, iyy1, izz, 1)
     1         + Ez_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ez_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ez_zc ( ixx1, iyy1, izz, 4) * z_1

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ez_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ez_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ez_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ez_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ez_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ez_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ez_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ez_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ez_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ez_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ez_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ez_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ez_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ez_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

! ----------

          Ex_2 = Ex_zc ( ixx2, iyy2, izz, 1)
     1         + Ex_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ex_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ex_zc ( ixx2, iyy2, izz, 4) * z_2
#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ex_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ex_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ex_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ex_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ex_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ex_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ex_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ex_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ex_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ex_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ex_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ex_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ex_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ex_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

          Ey_2 = Ey_zc ( ixx2, iyy2, izz, 1)
     1         + Ey_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ey_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ey_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ey_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ey_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ey_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ey_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ey_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ey_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ey_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ey_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ey_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ey_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ey_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ey_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ey_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ey_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

          Ez_2 = Ez_zc ( ixx2, iyy2, izz, 1)
     1         + Ez_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ez_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ez_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ez_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ez_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ez_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ez_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ez_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ez_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ez_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ez_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ez_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ez_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ez_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ez_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ez_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ez_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

! ----------

          Ex_3 = Ex_zc ( ixx3, iyy3, izz, 1)
     1         + Ex_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ex_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ex_zc ( ixx3, iyy3, izz, 4) * z_3
#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ex_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ex_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ex_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ex_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ex_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ex_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ex_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ex_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ex_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ex_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ex_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ex_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ex_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ex_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

          Ey_3 = Ey_zc ( ixx3, iyy3, izz, 1)
     1         + Ey_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ey_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ey_zc ( ixx3, iyy3, izz, 4) * z_3

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ey_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ey_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ey_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ey_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ey_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ey_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ey_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ey_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ey_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ey_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ey_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ey_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ey_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ey_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

          Ez_3 = Ez_zc ( ixx3, iyy3, izz, 1)
     1         + Ez_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ez_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ez_zc ( ixx3, iyy3, izz, 4) * z_3

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ez_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ez_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ez_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ez_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ez_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))
#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ez_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ez_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ez_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ez_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ez_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ez_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ez_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ez_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ez_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

! --------------------
#else
! --------------------
! Choose arithmetically averaged electric fields for he centered part of the
! 2D LLF Riemann solver. Electric field reconstruction is not needed.

          Ex_1 = ( Ex_zc ( ixx1, iyy1, izz, 1)
     1           + Ex_zc ( ixx2, iyy2, izz, 1)
     1           + Ex_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ex_2 = Ex_1  ;  Ex_3 = Ex_1

          Ey_1 = ( Ey_zc ( ixx1, iyy1, izz, 1)
     1           + Ey_zc ( ixx2, iyy2, izz, 1)
     1           + Ey_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ey_2 = Ey_1  ;  Ey_3 = Ey_1

          Ez_1 = ( Ez_zc ( ixx1, iyy1, izz, 1)
     1           + Ez_zc ( ixx2, iyy2, izz, 1)
     1           + Ez_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ez_2 = Ez_1  ;  Ez_3 = Ez_1

! --------------------
#endif
! --------------------

! ------------------------------
! Build the electric fields in the local z-direction; which is the global r-direction.

          tempaa = Ex_1 * n3 ( 1) + Ey_1 * n3 ( 2) + Ez_1 * n3 ( 3)
          Ez_1 = tempaa

          tempaa = Ex_2 * n3 ( 1) + Ey_2 * n3 ( 2) + Ez_2 * n3 ( 3)
          Ez_2 = tempaa

          tempaa = Ex_3 * n3 ( 1) + Ey_3 * n3 ( 2) + Ez_3 * n3 ( 3)
          Ez_3 = tempaa

! ------------------------------
! Build "Bneta_1, Bneta_2, Bneta_3" They may not be oriented in the positive directions.
! "Bneta_1" is along "neta_1"; "Bneta_2" is along "neta_2"; "Bneta_3" is along "neta_3".

          Bneta_1 = - Byf_quadpts ( ixx1, iyy1 - 1, izz, 5)
     1  + Bx0 * neta_1 ( 1) + By0 * neta_1 ( 2) + Bz0 * neta_1 ( 3)

          Bneta_2 = Bxf_quadpts ( ixx2, iyy2, izz, 3)
     1  + Bx0 * neta_2 ( 1) + By0 * neta_2 ( 2) + Bz0 * neta_2 ( 3)

          Bneta_3 = Byf_quadpts ( ixx3, iyy3, izz, 4)
     1  + Bx0 * neta_3 ( 1) + By0 * neta_3 ( 2) + Bz0 * neta_3 ( 3)

! ------------------------------

          S_R1 = - S_L_gr ( ixx3, iyy3, izz)
          S_L1 = - S_R_gr ( ixx3, iyy3, izz)

          S_R2 = S_R_gr ( ixx2, iyy2, izz)
          S_L2 = S_L_gr ( ixx2, iyy2, izz)

          S_R3 = S_U_gr ( ixx3, iyy3, izz)
          S_L3 = S_D_gr ( ixx3, iyy3, izz)

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ != 0)
! --------------------
! Make wave speeds fully subsonic for the 2D RS.

          tempaa = AMAX1 ( ABS ( S_R1), ABS ( S_L1),
     1                     ABS ( S_R2), ABS ( S_L2),
     1                     ABS ( S_R3), ABS ( S_L3))

          S_R1 = tempaa  ;  S_L1 = - tempaa
          S_R2 = tempaa  ;  S_L2 = - tempaa
          S_R3 = tempaa  ;  S_L3 = - tempaa

! --------------------
#endif
! --------------------

! ------------------------------------------------------------
! We are only making z-directional calls to the 2D Riemann solver that has three
! input states.

          CALL Two_D_RS_Ez_3State_HLL

     1         ( S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,

     1           Ez_1, Ez_2, Ez_3,

     1           Bneta_1, Bneta_2, Bneta_3,

     1           BSSx, BSSy, ESSz)

! ------------------------------------------------------------
! Hand back the global electric field "ESSz" from the 2D Riemann solver.

          Efz ( ixx, iyy, izz) = Efz ( ixx, iyy, izz) + wt * ESSz

! ------------------------------------------------------------

        END DO  ! End of "DO izz = 1, nr_zones".

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, 3" loop over quadrature points.

      END IF  ! End of "IF ( itri_patch_RD .EQ. 1) THEN".

! ----------------------------------------------------------------------
! STEP 4.c) The Right-Up corner of this mesh is a tri-corner.

      IF ( itri_patch_RU .EQ. 1) THEN

      ixx = ntp_zones  ! These are the indices of the z-edge along which we want the e-field.
      iyy = ntp_zones

! At third and fourth order, loop over 3 quadrature points in the edge of interest.
      DO j = 1, 3

#if (NUMORDER == 2)
! At second order, we need only one nodal point.
      IF ( j .GT. 1) EXIT
#endif

      wt = 1.0
#if (NUMORDER >= 3)
      wt = 1.0 / 6.0  ! These are the weights for simpson integration.
      IF ( j .EQ. 1) wt = 4.0 / 6.0
#endif

!$OMP PARALLEL DO PRIVATE (
!$OMP&    izz, ixx1, iyy1, ixx2, iyy2, ixx3, iyy3,                      &
!$OMP&    xec, yec, zec, tempaa, obl_1, obl_2, obl_3,                   &
!$OMP&    n1 , n2 , n3, B0, time_use, Bx0, By0, Bz0,                    &
!$OMP&    tau_1, neta_1, tau_2, neta_2, tau_3, neta_3,                  &

!$OMP&    x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3,                  &
!$OMP&    Ex_1, Ey_1, Ez_1, Ex_2, Ey_2, Ez_2, Ex_3, Ey_3, Ez_3,         &
!$OMP&    S_R1, S_L1, S_R2, S_L2, S_R3, S_L3, Bneta_1, Bneta_2, Bneta_3,&
!$OMP&    BSSx, BSSy, ESSz)                                             &

!$OMP&             SHARED ( nr_zones, ixx, iyy, smallnum, ob3, j, wt,   &
!$OMP&    vert_coords_patch, v_scale_length, v_centroid_patch,          &
!$OMP&    Coef_Basis_Moments,                                           &

!$OMP&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$OMP&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$OMP&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$OMP&    Efz)                                                          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum, ob3, j, wt,                                         &
!$acc&    vert_coords_patch, v_scale_length, v_centroid_patch,          &
!$acc&    Coef_Basis_Moments,                                           &

!$acc&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$acc&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$acc&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$acc&    Efz)

!$acc loop vector independent                                           &
!$acc&    private ( izz, ixx1, iyy1, ixx2, iyy2, ixx3, iyy3,            &
!$acc&    xec, yec, zec, tempaa, obl_1, obl_2, obl_3,                   &
!$acc&    n1 , n2 , n3, B0, time_use, Bx0, By0, Bz0,                    &
!$acc&    tau_1, neta_1, tau_2, neta_2, tau_3, neta_3,                  &

!$acc&    x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3,                  &
!$acc&    Ex_1, Ey_1, Ez_1, Ex_2, Ey_2, Ez_2, Ex_3, Ey_3, Ez_3,         &
!$acc&    S_R1, S_L1, S_R2, S_L2, S_R3, S_L3, Bneta_1, Bneta_2, Bneta_3,&
!$acc&    BSSx, BSSy, ESSz)

        DO izz = 1, nr_zones

! ------------------------------------------------------------

! Fill in "ixx1, iyy1, ixx2, iyy2, ixx3, iyy3" in counter-clockwise sense.

        ixx1 = ixx      ;  iyy1 = iyy + 1
        ixx2 = ixx      ;  iyy2 = iyy
        ixx3 = ixx + 1  ;  iyy3 = iyy

! --------------------------------------------------
! Find the edge-center "xec, yec, zec" by bisecting the two z-vertices.
! No need to rescale by radius since this is in the radial direction by definition.

          IF ( j .EQ. 1) THEN  ! First quadrature point is the center point.

            xec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 1)  ! Bisect the z-vertices.
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 1) )
            yec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 2)
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 2) )
            zec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 3)
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 3) )

          ELSE IF ( j .EQ. 2) THEN  ! Second quadrature point is one of the end points.

            xec = vert_coords_patch ( ixx, iyy, izz, 1)
            yec = vert_coords_patch ( ixx, iyy, izz, 2)
            zec = vert_coords_patch ( ixx, iyy, izz, 3)

          ELSE IF ( j .EQ. 3) THEN  ! Third quadrature point is the other end point.

            xec = vert_coords_patch ( ixx, iyy, izz - 1, 1)
            yec = vert_coords_patch ( ixx, iyy, izz - 1, 2)
            zec = vert_coords_patch ( ixx, iyy, izz - 1, 3)

          END IF

! ------------------------------
! Feed ( xec, yec, zec) to the tiny subroutine that evaluates "B0".
! Default is zero.

          B0 = 0.0
          Bx0 = 0.0
          By0 = 0.0
          Bz0 = 0.0

#if ( USE_B0B1_FORMULATION == 1)

          bflag = 1
          time_use = 0

          CALL GET_MAGNETIC_FIELD ( bflag, time_use,  xec, yec, zec,
     1                              Bx0, By0, Bz0)

          B0 ( 1) = Bx0
          B0 ( 2) = By0
          B0 ( 3) = Bz0

#endif

! ------------------------------

        obl_1 = 1.0 / AMAX1 ( v_scale_length ( ixx1, iyy1, izz),
     1                        smallnum)

        obl_2 = 1.0 / AMAX1 ( v_scale_length ( ixx2, iyy2, izz),
     1                        smallnum)

        obl_3 = 1.0 / AMAX1 ( v_scale_length ( ixx3, iyy3, izz),
     1                        smallnum)

! ------------------------------
! Build scaled "x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3".
! These are needed in general if we are using the reconstruction. But can
! be bypassed at second order.

        x_1 = ( xec - v_centroid_patch ( ixx1, iyy1, izz, 1) )
     1        * obl_1

        y_1 = ( yec - v_centroid_patch ( ixx1, iyy1, izz, 2) )
     1        * obl_1

        z_1 = ( zec - v_centroid_patch ( ixx1, iyy1, izz, 3) )
     1        * obl_1

! ----------

        x_2 = ( xec - v_centroid_patch ( ixx2, iyy2, izz, 1) )
     1        * obl_2

        y_2 = ( yec - v_centroid_patch ( ixx2, iyy2, izz, 2) )
     1        * obl_2

        z_2 = ( zec - v_centroid_patch ( ixx2, iyy2, izz, 3) )
     1        * obl_2

! ----------

        x_3 = ( xec - v_centroid_patch ( ixx3, iyy3, izz, 1) )
     1        * obl_3

        y_3 = ( yec - v_centroid_patch ( ixx3, iyy3, izz, 2) )
     1        * obl_3

        z_3 = ( zec - v_centroid_patch ( ixx3, iyy3, izz, 3) )
     1        * obl_3

! ------------------------------
! Since this is a tri-corner, we have to build:- "n3" along radial direction.
! "tau_1" along y-axis of mesh and "neta_1 = tau_1 X n3" perpendicular to it by using "n3".
! "tau_2" along -ve x-axis of mesh and "neta_2 = tau_2 X n3" perpendicular to it by using "n3".
! "tau_3" along -ve y-axis of mesh and "neta_3 = tau_3 X n3" perpendicular to it by using "n3".

          tempaa = 1.0 / SQRT ( xec**2 + yec**2 + zec**2)

! Here "n3" is in the global r-direction of the mesh, which is our local z-direction.
          n3 ( 1) = xec * tempaa   ! This is in the global r-direction.
          n3 ( 2) = yec * tempaa
          n3 ( 3) = zec * tempaa

! ----------

          tau_1 ( 1) = vert_coords_patch ( ixx, iyy + 1, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_1 ( 2) = vert_coords_patch ( ixx, iyy + 1, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_1 ( 3) = vert_coords_patch ( ixx, iyy + 1, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_1 ( 1) * n3 ( 1) + tau_1 ( 2) * n3 ( 2)
     1           + tau_1 ( 3) * n3 ( 3)

          tau_1 = tau_1 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_1 ( 1)**2 + tau_1 ( 2)**2
     1                        + tau_1 ( 3)**2)

          tau_1 ( :) = tau_1 ( :) * tempaa

          neta_1 ( 1) = tau_1 ( 2) * n3 ( 3) - tau_1 ( 3) * n3 ( 2)
          neta_1 ( 2) = tau_1 ( 3) * n3 ( 1) - tau_1 ( 1) * n3 ( 3)
          neta_1 ( 3) = tau_1 ( 1) * n3 ( 2) - tau_1 ( 2) * n3 ( 1)

! ----------

          tau_2 ( 1) = vert_coords_patch ( ixx - 1, iyy, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_2 ( 2) = vert_coords_patch ( ixx - 1, iyy, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_2 ( 3) = vert_coords_patch ( ixx - 1, iyy, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_2 ( 1) * n3 ( 1) + tau_2 ( 2) * n3 ( 2)
     1           + tau_2 ( 3) * n3 ( 3)

          tau_2 = tau_2 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_2 ( 1)**2 + tau_2 ( 2)**2
     1                        + tau_2 ( 3)**2)

          tau_2 ( :) = tau_2 ( :) * tempaa

          neta_2 ( 1) = tau_2 ( 2) * n3 ( 3) - tau_2 ( 3) * n3 ( 2)
          neta_2 ( 2) = tau_2 ( 3) * n3 ( 1) - tau_2 ( 1) * n3 ( 3)
          neta_2 ( 3) = tau_2 ( 1) * n3 ( 2) - tau_2 ( 2) * n3 ( 1)

! ----------

          tau_3 ( 1) = vert_coords_patch ( ixx, iyy - 1, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_3 ( 2) = vert_coords_patch ( ixx, iyy - 1, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_3 ( 3) = vert_coords_patch ( ixx, iyy - 1, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_3 ( 1) * n3 ( 1) + tau_3 ( 2) * n3 ( 2)
     1           + tau_3 ( 3) * n3 ( 3)

          tau_3 = tau_3 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_3 ( 1)**2 + tau_3 ( 2)**2
     1                        + tau_3 ( 3)**2)

          tau_3 ( :) = tau_3 ( :) * tempaa

          neta_3 ( 1) = tau_3 ( 2) * n3 ( 3) - tau_3 ( 3) * n3 ( 2)
          neta_3 ( 2) = tau_3 ( 3) * n3 ( 1) - tau_3 ( 1) * n3 ( 3)
          neta_3 ( 3) = tau_3 ( 1) * n3 ( 2) - tau_3 ( 2) * n3 ( 1)

! --------------------------------------------------
! Now that we have the properly scaled relative coordinates from each of the
! three particpating zones, we can use these to obtain the three input electric field
! components "Ez_1, Ez_2, Ez_3" that go into the 2D Riemann solver
! with 3 states.

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ == 0)
! --------------------
! Here the building of the electric fields requires that they have been reconstructed.

          Ex_1 = Ex_zc ( ixx1, iyy1, izz, 1)
     1         + Ex_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ex_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ex_zc ( ixx1, iyy1, izz, 4) * z_1
#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ex_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ex_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ex_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ex_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ex_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ex_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ex_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ex_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ex_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ex_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ex_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ex_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ex_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ex_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

          Ey_1 = Ey_zc ( ixx1, iyy1, izz, 1)
     1         + Ey_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ey_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ey_zc ( ixx1, iyy1, izz, 4) * z_1

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ey_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ey_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ey_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ey_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ey_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ey_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ey_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ey_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ey_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ey_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ey_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ey_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ey_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ey_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

          Ez_1 = Ez_zc ( ixx1, iyy1, izz, 1)
     1         + Ez_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ez_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ez_zc ( ixx1, iyy1, izz, 4) * z_1

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ez_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ez_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ez_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ez_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ez_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ez_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ez_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ez_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ez_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ez_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ez_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ez_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ez_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ez_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

! ----------

          Ex_2 = Ex_zc ( ixx2, iyy2, izz, 1)
     1         + Ex_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ex_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ex_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ex_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ex_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ex_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ex_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ex_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ex_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ex_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ex_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ex_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ex_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ex_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ex_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ex_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ex_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

          Ey_2 = Ey_zc ( ixx2, iyy2, izz, 1)
     1         + Ey_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ey_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ey_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ey_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ey_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ey_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ey_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ey_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ey_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ey_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ey_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ey_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ey_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ey_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ey_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ey_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ey_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

          Ez_2 = Ez_zc ( ixx2, iyy2, izz, 1)
     1         + Ez_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ez_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ez_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ez_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ez_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ez_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ez_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ez_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ez_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ez_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ez_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ez_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ez_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ez_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ez_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ez_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ez_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

! ----------

          Ex_3 = Ex_zc ( ixx3, iyy3, izz, 1)
     1         + Ex_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ex_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ex_zc ( ixx3, iyy3, izz, 4) * z_3

#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ex_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ex_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ex_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ex_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ex_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ex_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ex_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ex_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ex_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ex_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ex_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ex_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ex_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ex_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

          Ey_3 = Ey_zc ( ixx3, iyy3, izz, 1)
     1         + Ey_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ey_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ey_zc ( ixx3, iyy3, izz, 4) * z_3

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ey_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ey_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ey_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ey_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ey_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ey_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ey_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ey_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ey_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ey_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ey_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ey_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ey_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ey_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

          Ez_3 = Ez_zc ( ixx3, iyy3, izz, 1)
     1         + Ez_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ez_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ez_zc ( ixx3, iyy3, izz, 4) * z_3

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ez_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ez_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ez_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ez_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ez_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ez_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ez_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ez_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ez_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ez_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ez_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ez_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ez_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ez_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

! --------------------
#else
! --------------------
! Choose arithmetically averaged electric fields for he centered part of the
! 2D LLF Riemann solver. Electric field reconstruction is not needed.

          Ex_1 = ( Ex_zc ( ixx1, iyy1, izz, 1)
     1           + Ex_zc ( ixx2, iyy2, izz, 1)
     1           + Ex_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ex_2 = Ex_1  ;  Ex_3 = Ex_1

          Ey_1 = ( Ey_zc ( ixx1, iyy1, izz, 1)
     1           + Ey_zc ( ixx2, iyy2, izz, 1)
     1           + Ey_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ey_2 = Ey_1  ;  Ey_3 = Ey_1

          Ez_1 = ( Ez_zc ( ixx1, iyy1, izz, 1)
     1           + Ez_zc ( ixx2, iyy2, izz, 1)
     1           + Ez_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ez_2 = Ez_1  ;  Ez_3 = Ez_1

! --------------------
#endif
! --------------------

! ------------------------------
! Build the electric fields in the local z-direction; which is the global r-direction.

          tempaa = Ex_1 * n3 ( 1) + Ey_1 * n3 ( 2) + Ez_1 * n3 ( 3)
          Ez_1 = tempaa

          tempaa = Ex_2 * n3 ( 1) + Ey_2 * n3 ( 2) + Ez_2 * n3 ( 3)
          Ez_2 = tempaa

          tempaa = Ex_3 * n3 ( 1) + Ey_3 * n3 ( 2) + Ez_3 * n3 ( 3)
          Ez_3 = tempaa

! ------------------------------
! Build "Bneta_1, Bneta_2, Bneta_3" They may not be oriented in the positive directions.
! "Bneta_1" is along "neta_1"; "Bneta_2" is along "neta_2"; "Bneta_3" is along "neta_3".

          Bneta_1 = Bxf_quadpts ( ixx1, iyy1, izz, 3)
     1  + Bx0 * neta_1 ( 1) + By0 * neta_1 ( 2) + Bz0 * neta_1 ( 3)

          Bneta_2 = Byf_quadpts ( ixx2, iyy2, izz, 4)
     1  + Bx0 * neta_2 ( 1) + By0 * neta_2 ( 2) + Bz0 * neta_2 ( 3)

          Bneta_3 = - Bxf_quadpts ( ixx2, iyy2, izz, 2)
     1  + Bx0 * neta_3 ( 1) + By0 * neta_3 ( 2) + Bz0 * neta_3 ( 3)

! ------------------------------
! Evaluate "S_R1, S_L1, S_R2, S_L2, S_R3, S_L3". Watch out for the orientatations
! because in each of the three zone faces we take the positive direction to be
! along "neta_1, neta_2, neta_3".

          S_R1 = S_R_gr ( ixx1, iyy1, izz)
          S_L1 = S_L_gr ( ixx1, iyy1, izz)

          S_R2 = S_U_gr ( ixx2, iyy2, izz)
          S_L2 = S_D_gr ( ixx2, iyy2, izz)

          S_R3 = - S_L_gr ( ixx2, iyy2, izz)
          S_L3 = - S_R_gr ( ixx2, iyy2, izz)

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ != 0)
! --------------------
! Make wave speeds fully subsonic for the 2D RS.

          tempaa = AMAX1 ( ABS ( S_R1), ABS ( S_L1),
     1                     ABS ( S_R2), ABS ( S_L2),
     1                     ABS ( S_R3), ABS ( S_L3))

          S_R1 = tempaa  ;  S_L1 = - tempaa
          S_R2 = tempaa  ;  S_L2 = - tempaa
          S_R3 = tempaa  ;  S_L3 = - tempaa

! --------------------
#endif
! --------------------

! ------------------------------------------------------------
! We are only making z-directional calls to the 2D Riemann solver that has three
! input states.

          CALL Two_D_RS_Ez_3State_HLL

     1         ( S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,

     1           Ez_1, Ez_2, Ez_3,

     1           Bneta_1, Bneta_2, Bneta_3,

     1           BSSx, BSSy, ESSz)

! ------------------------------------------------------------
! Hand back the global electric field "ESSz" from the 2D Riemann solver.

          Efz ( ixx, iyy, izz) = Efz ( ixx, iyy, izz) + wt * ESSz

! ------------------------------------------------------------

        END DO  ! End of "DO izz = 1, nr_zones".

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, 3" loop over quadrature points.

      END IF  ! End of "IF ( itri_patch_RU .EQ. 1) THEN".

! ----------------------------------------------------------------------
! STEP 4.d) The Left-Up corner of this mesh is a tri-corner.

      IF ( itri_patch_LU .EQ. 1) THEN

      ixx = 0  ! These are the indices of the z-edge along which we want the e-field.
      iyy = ntp_zones

! At third and fourth order, loop over 3 quadrature points in the edge of interest.
      DO j = 1, 3

#if (NUMORDER == 2)
! At second order, we need only one nodal point.
      IF ( j .GT. 1) EXIT
#endif

      wt = 1.0
#if (NUMORDER >= 3)
      wt = 1.0 / 6.0  ! These are the weights for simpson integration.
      IF ( j .EQ. 1) wt = 4.0 / 6.0
#endif

!$OMP PARALLEL DO PRIVATE (
!$OMP&    izz, ixx1, iyy1, ixx2, iyy2, ixx3, iyy3,                      &
!$OMP&    xec, yec, zec, tempaa, obl_1, obl_2, obl_3,                   &
!$OMP&    n1 , n2 , n3, B0, time_use, Bx0, By0, Bz0,                    &
!$OMP&    tau_1, neta_1, tau_2, neta_2, tau_3, neta_3,                  &

!$OMP&    x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3,                  &
!$OMP&    Ex_1, Ey_1, Ez_1, Ex_2, Ey_2, Ez_2, Ex_3, Ey_3, Ez_3,         &
!$OMP&    S_R1, S_L1, S_R2, S_L2, S_R3, S_L3, Bneta_1, Bneta_2, Bneta_3,&
!$OMP&    BSSx, BSSy, ESSz)                                             &

!$OMP&             SHARED ( nr_zones, ixx, iyy, smallnum, ob3, j, wt,   &
!$OMP&    vert_coords_patch, v_scale_length, v_centroid_patch,          &
!$OMP&    Coef_Basis_Moments,                                           &

!$OMP&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$OMP&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$OMP&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$OMP&    Efz)                                                          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum, ob3, j, wt,                                         &
!$acc&    vert_coords_patch, v_scale_length, v_centroid_patch,          &
!$acc&    Coef_Basis_Moments,                                           &

!$acc&    U, Ex_zc, Ey_zc, Ez_zc,                                       &
!$acc&    S_R_gr, S_L_gr, S_U_gr, S_D_gr, S_T_gr, S_B_gr,               &
!$acc&    Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                        &
!$acc&    Efz)

!$acc loop vector independent                                           &
!$acc&    private ( izz, ixx1, iyy1, ixx2, iyy2, ixx3, iyy3,            &
!$acc&    xec, yec, zec, tempaa, obl_1, obl_2, obl_3,                   &
!$acc&    n1 , n2 , n3, B0, time_use, Bx0, By0, Bz0,                    &
!$acc&    tau_1, neta_1, tau_2, neta_2, tau_3, neta_3,                  &

!$acc&    x_1, x_2, x_3, y_1, y_2, y_3, z_1, z_2, z_3,                  &
!$acc&    Ex_1, Ey_1, Ez_1, Ex_2, Ey_2, Ez_2, Ex_3, Ey_3, Ez_3,         &
!$acc&    S_R1, S_L1, S_R2, S_L2, S_R3, S_L3, Bneta_1, Bneta_2, Bneta_3,&
!$acc&    BSSx, BSSy, ESSz)

        DO izz = 1, nr_zones

! ------------------------------------------------------------

! Fill in "ixx1, iyy1, ixx2, iyy2, ixx3, iyy3" in counter-clockwise sense.

        ixx1 = ixx + 1  ;  iyy1 = iyy + 1
        ixx2 = ixx      ;  iyy2 = iyy
        ixx3 = ixx + 1  ;  iyy3 = iyy

! --------------------------------------------------
! Find the edge-center "xec, yec, zec" by bisecting the two z-vertices.
! No need to rescale by radius since this is in the radial direction by definition.

          IF ( j .EQ. 1) THEN  ! First quadrature point is the center point.

            xec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 1)  ! Bisect the z-vertices.
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 1) )
            yec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 2)
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 2) )
            zec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 3)
     1                  + vert_coords_patch ( ixx, iyy, izz - 1, 3) )

          ELSE IF ( j .EQ. 2) THEN  ! Second quadrature point is one of the end points.

            xec = vert_coords_patch ( ixx, iyy, izz, 1)
            yec = vert_coords_patch ( ixx, iyy, izz, 2)
            zec = vert_coords_patch ( ixx, iyy, izz, 3)

          ELSE IF ( j .EQ. 3) THEN  ! Third quadrature point is the other end point.

            xec = vert_coords_patch ( ixx, iyy, izz - 1, 1)
            yec = vert_coords_patch ( ixx, iyy, izz - 1, 2)
            zec = vert_coords_patch ( ixx, iyy, izz - 1, 3)

          END IF

! ------------------------------
! Feed ( xec, yec, zec) to the tiny subroutine that evaluates "B0".
! Default is zero.

          B0 = 0.0
          Bx0 = 0.0
          By0 = 0.0
          Bz0 = 0.0

#if ( USE_B0B1_FORMULATION == 1)

          bflag = 1
          time_use = 0

          CALL GET_MAGNETIC_FIELD ( bflag, time_use,  xec, yec, zec,
     1                              Bx0, By0, Bz0)

          B0 ( 1) = Bx0
          B0 ( 2) = By0
          B0 ( 3) = Bz0

#endif

! ------------------------------

        obl_1 = 1.0 / AMAX1 ( v_scale_length ( ixx1, iyy1, izz),
     1                        smallnum)

        obl_2 = 1.0 / AMAX1 ( v_scale_length ( ixx2, iyy2, izz),
     1                        smallnum)

        obl_3 = 1.0 / AMAX1 ( v_scale_length ( ixx3, iyy3, izz),
     1                        smallnum)

! ------------------------------
! Build scaled "x_1, y_1, z_1, x_2, y_2, z_2, x_3, y_3, z_3".
! These are needed in general if we are using the reconstruction. But can
! be bypassed at second order.

        x_1 = ( xec - v_centroid_patch ( ixx1, iyy1, izz, 1) )
     1        * obl_1

        y_1 = ( yec - v_centroid_patch ( ixx1, iyy1, izz, 2) )
     1        * obl_1

        z_1 = ( zec - v_centroid_patch ( ixx1, iyy1, izz, 3) )
     1        * obl_1

! ----------

        x_2 = ( xec - v_centroid_patch ( ixx2, iyy2, izz, 1) )
     1        * obl_2

        y_2 = ( yec - v_centroid_patch ( ixx2, iyy2, izz, 2) )
     1        * obl_2

        z_2 = ( zec - v_centroid_patch ( ixx2, iyy2, izz, 3) )
     1        * obl_2

! ----------

        x_3 = ( xec - v_centroid_patch ( ixx3, iyy3, izz, 1) )
     1        * obl_3

        y_3 = ( yec - v_centroid_patch ( ixx3, iyy3, izz, 2) )
     1        * obl_3

        z_3 = ( zec - v_centroid_patch ( ixx3, iyy3, izz, 3) )
     1        * obl_3

! ------------------------------
! Since this is a tri-corner, we have to build:- "n3" along radial direction.
! "tau_1" along x-axis of mesh and "neta_1 = tau_1 X n3" perpendicular to it by using "n3".
! "tau_2" along y-axis of mesh and "neta_2 = tau_2 X n3" perpendicular to it by using "n3".
! "tau_3" along -ve y-axis of mesh and "neta_3 = tau_3 X n3" perpendicular to it by using "n3".

          tempaa = 1.0 / SQRT ( xec**2 + yec**2 + zec**2)

! Here "n3" is in the global r-direction of the mesh, which is our local z-direction.
          n3 ( 1) = xec * tempaa   ! This is in the global r-direction.
          n3 ( 2) = yec * tempaa
          n3 ( 3) = zec * tempaa

! ----------

          tau_1 ( 1) = vert_coords_patch ( ixx + 1, iyy, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_1 ( 2) = vert_coords_patch ( ixx + 1, iyy, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_1 ( 3) = vert_coords_patch ( ixx + 1, iyy, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_1 ( 1) * n3 ( 1) + tau_1 ( 2) * n3 ( 2)
     1           + tau_1 ( 3) * n3 ( 3)

          tau_1 = tau_1 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_1 ( 1)**2 + tau_1 ( 2)**2
     1                        + tau_1 ( 3)**2)

          tau_1 ( :) = tau_1 ( :) * tempaa

          neta_1 ( 1) = tau_1 ( 2) * n3 ( 3) - tau_1 ( 3) * n3 ( 2)
          neta_1 ( 2) = tau_1 ( 3) * n3 ( 1) - tau_1 ( 1) * n3 ( 3)
          neta_1 ( 3) = tau_1 ( 1) * n3 ( 2) - tau_1 ( 2) * n3 ( 1)

! ----------

          tau_2 ( 1) = vert_coords_patch ( ixx, iyy + 1, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_2 ( 2) = vert_coords_patch ( ixx, iyy + 1, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_2 ( 3) = vert_coords_patch ( ixx, iyy + 1, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_2 ( 1) * n3 ( 1) + tau_2 ( 2) * n3 ( 2)
     1           + tau_2 ( 3) * n3 ( 3)

          tau_2 = tau_2 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_2 ( 1)**2 + tau_2 ( 2)**2
     1                        + tau_2 ( 3)**2)

          tau_2 ( :) = tau_2 ( :) * tempaa

          neta_2 ( 1) = tau_2 ( 2) * n3 ( 3) - tau_2 ( 3) * n3 ( 2)
          neta_2 ( 2) = tau_2 ( 3) * n3 ( 1) - tau_2 ( 1) * n3 ( 3)
          neta_2 ( 3) = tau_2 ( 1) * n3 ( 2) - tau_2 ( 2) * n3 ( 1)

! ----------

          tau_3 ( 1) = vert_coords_patch ( ixx, iyy - 1, izz, 1)
     1               - vert_coords_patch ( ixx, iyy, izz, 1)
          tau_3 ( 2) = vert_coords_patch ( ixx, iyy - 1, izz, 2)
     1               - vert_coords_patch ( ixx, iyy, izz, 2)
          tau_3 ( 3) = vert_coords_patch ( ixx, iyy - 1, izz, 3)
     1               - vert_coords_patch ( ixx, iyy, izz, 3)

          tempaa = tau_3 ( 1) * n3 ( 1) + tau_3 ( 2) * n3 ( 2)
     1           + tau_3 ( 3) * n3 ( 3)

          tau_3 = tau_3 - tempaa * n3

          tempaa = 1.0 / SQRT ( tau_3 ( 1)**2 + tau_3 ( 2)**2
     1                        + tau_3 ( 3)**2)

          tau_3 ( :) = tau_3 ( :) * tempaa

          neta_3 ( 1) = tau_3 ( 2) * n3 ( 3) - tau_3 ( 3) * n3 ( 2)
          neta_3 ( 2) = tau_3 ( 3) * n3 ( 1) - tau_3 ( 1) * n3 ( 3)
          neta_3 ( 3) = tau_3 ( 1) * n3 ( 2) - tau_3 ( 2) * n3 ( 1)

! --------------------------------------------------
! Now that we have the properly scaled relative coordinates from each of the
! three particpating zones, we can use these to obtain the three input electric field
! components "Ez_1, Ez_2, Ez_3" that go into the 2D Riemann solver
! with 3 states.

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ == 0)
! --------------------
! Here the building of the electric fields requires that they have been reconstructed.

          Ex_1 = Ex_zc ( ixx1, iyy1, izz, 1)
     1         + Ex_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ex_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ex_zc ( ixx1, iyy1, izz, 4) * z_1

#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ex_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ex_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ex_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ex_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ex_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ex_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ex_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ex_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ex_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ex_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ex_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ex_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ex_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ex_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

          Ey_1 = Ey_zc ( ixx1, iyy1, izz, 1)
     1         + Ey_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ey_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ey_zc ( ixx1, iyy1, izz, 4) * z_1

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ey_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ey_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ey_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ey_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ey_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ey_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ey_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ey_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ey_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ey_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ey_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ey_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ey_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ey_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

          Ez_1 = Ez_zc ( ixx1, iyy1, izz, 1)
     1         + Ez_zc ( ixx1, iyy1, izz, 2) * x_1
     1         + Ez_zc ( ixx1, iyy1, izz, 3) * y_1
     1         + Ez_zc ( ixx1, iyy1, izz, 4) * z_1

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx1, iyy1, izz, 5)
     1         * (x_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 5))
     1         + Ez_zc ( ixx1, iyy1, izz, 6)
     1         * (y_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 6))
     1         + Ez_zc ( ixx1, iyy1, izz, 7)
     1         * (z_1**2 - Coef_Basis_Moments ( ixx1, iyy1, izz, 7))
     1         + Ez_zc ( ixx1, iyy1, izz, 8)
     1         * (x_1 * y_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 8))
     1         + Ez_zc ( ixx1, iyy1, izz, 9)
     1         * (y_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 9))
     1         + Ez_zc ( ixx1, iyy1, izz, 10)
     1         * (x_1 * z_1 - Coef_Basis_Moments ( ixx1, iyy1, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx1, iyy1, izz, 11)
     1         * (x_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 11))
     1         + Ez_zc ( ixx1, iyy1, izz, 12)
     1         * (y_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 12))
     1         + Ez_zc ( ixx1, iyy1, izz, 13)
     1         * (z_1**3 - Coef_Basis_Moments ( ixx1, iyy1, izz, 13))
     1         + Ez_zc ( ixx1, iyy1, izz, 14)
     1         * (x_1**2 * y_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 14))
     1         + Ez_zc ( ixx1, iyy1, izz, 15)
     1         * (x_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 15))
     1         + Ez_zc ( ixx1, iyy1, izz, 16)
     1         * (x_1 * y_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 16))
     1         + Ez_zc ( ixx1, iyy1, izz, 17)
     1         * (y_1**2 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 17))
     1         + Ez_zc ( ixx1, iyy1, izz, 18)
     1         * (x_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 18))
     1         + Ez_zc ( ixx1, iyy1, izz, 19)
     1         * (y_1 * z_1**2
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 19))
     1         + Ez_zc ( ixx1, iyy1, izz, 20)
     1         * (x_1 * y_1 * z_1
     1         - Coef_Basis_Moments ( ixx1, iyy1, izz, 20))

#endif

! ----------

          Ex_2 = Ex_zc ( ixx2, iyy2, izz, 1)
     1         + Ex_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ex_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ex_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ex_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ex_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ex_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ex_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ex_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ex_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ex_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ex_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ex_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ex_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ex_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ex_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ex_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ex_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

          Ey_2 = Ey_zc ( ixx2, iyy2, izz, 1)
     1         + Ey_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ey_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ey_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ey_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ey_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ey_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ey_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ey_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ey_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ey_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ey_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ey_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ey_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ey_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ey_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ey_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ey_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

          Ez_2 = Ez_zc ( ixx2, iyy2, izz, 1)
     1         + Ez_zc ( ixx2, iyy2, izz, 2) * x_2
     1         + Ez_zc ( ixx2, iyy2, izz, 3) * y_2
     1         + Ez_zc ( ixx2, iyy2, izz, 4) * z_2

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx2, iyy2, izz, 5)
     1         * (x_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 5))
     1         + Ez_zc ( ixx2, iyy2, izz, 6)
     1         * (y_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 6))
     1         + Ez_zc ( ixx2, iyy2, izz, 7)
     1         * (z_2**2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 7))
     1         + Ez_zc ( ixx2, iyy2, izz, 8)
     1         * (x_2 * y_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 8))
     1         + Ez_zc ( ixx2, iyy2, izz, 9)
     1         * (y_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 9))
     1         + Ez_zc ( ixx2, iyy2, izz, 10)
     1         * (x_2 * z_2 - Coef_Basis_Moments ( ixx2, iyy2, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx2, iyy2, izz, 11)
     1         * (x_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 11))
     1         + Ez_zc ( ixx2, iyy2, izz, 12)
     1         * (y_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 12))
     1         + Ez_zc ( ixx2, iyy2, izz, 13)
     1         * (z_2**3 - Coef_Basis_Moments ( ixx2, iyy2, izz, 13))
     1         + Ez_zc ( ixx2, iyy2, izz, 14)
     1         * (x_2**2 * y_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 14))
     1         + Ez_zc ( ixx2, iyy2, izz, 15)
     1         * (x_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 15))
     1         + Ez_zc ( ixx2, iyy2, izz, 16)
     1         * (x_2 * y_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 16))
     1         + Ez_zc ( ixx2, iyy2, izz, 17)
     1         * (y_2**2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 17))
     1         + Ez_zc ( ixx2, iyy2, izz, 18)
     1         * (x_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 18))
     1         + Ez_zc ( ixx2, iyy2, izz, 19)
     1         * (y_2 * z_2**2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 19))
     1         + Ez_zc ( ixx2, iyy2, izz, 20)
     1         * (x_2 * y_2 * z_2
     1         - Coef_Basis_Moments ( ixx2, iyy2, izz, 20))

#endif

! ----------

          Ex_3 = Ex_zc ( ixx3, iyy3, izz, 1)
     1         + Ex_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ex_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ex_zc ( ixx3, iyy3, izz, 4) * z_3

#if (NUMORDER >= 3)
     1         + Ex_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ex_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ex_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ex_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ex_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ex_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ex_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ex_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ex_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ex_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ex_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ex_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ex_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ex_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ex_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ex_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

          Ey_3 = Ey_zc ( ixx3, iyy3, izz, 1)
     1         + Ey_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ey_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ey_zc ( ixx3, iyy3, izz, 4) * z_3
#if (NUMORDER >= 3)
     1         + Ey_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ey_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ey_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ey_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ey_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ey_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ey_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ey_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ey_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ey_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ey_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ey_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ey_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ey_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ey_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ey_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

          Ez_3 = Ez_zc ( ixx3, iyy3, izz, 1)
     1         + Ez_zc ( ixx3, iyy3, izz, 2) * x_3
     1         + Ez_zc ( ixx3, iyy3, izz, 3) * y_3
     1         + Ez_zc ( ixx3, iyy3, izz, 4) * z_3

#if (NUMORDER >= 3)
     1         + Ez_zc ( ixx3, iyy3, izz, 5)
     1         * (x_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 5))
     1         + Ez_zc ( ixx3, iyy3, izz, 6)
     1         * (y_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 6))
     1         + Ez_zc ( ixx3, iyy3, izz, 7)
     1         * (z_3**2 - Coef_Basis_Moments ( ixx3, iyy3, izz, 7))
     1         + Ez_zc ( ixx3, iyy3, izz, 8)
     1         * (x_3 * y_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 8))
     1         + Ez_zc ( ixx3, iyy3, izz, 9)
     1         * (y_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 9))
     1         + Ez_zc ( ixx3, iyy3, izz, 10)
     1         * (x_3 * z_3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 10))

#elif (NUMORDER >= 4)
     1         + Ez_zc ( ixx3, iyy3, izz, 11)
     1         * (x_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 11))
     1         + Ez_zc ( ixx3, iyy3, izz, 12)
     1         * (y_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 12))
     1         + Ez_zc ( ixx3, iyy3, izz, 13)
     1         * (z_3**3 - Coef_Basis_Moments ( ixx3, iyy3, izz, 13))
     1         + Ez_zc ( ixx3, iyy3, izz, 14)
     1         * (x_3**2 * y_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 14))
     1         + Ez_zc ( ixx3, iyy3, izz, 15)
     1         * (x_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 15))
     1         + Ez_zc ( ixx3, iyy3, izz, 16)
     1         * (x_3 * y_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 16))
     1         + Ez_zc ( ixx3, iyy3, izz, 17)
     1         * (y_3**2 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 17))
     1         + Ez_zc ( ixx3, iyy3, izz, 18)
     1         * (x_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 18))
     1         + Ez_zc ( ixx3, iyy3, izz, 19)
     1         * (y_3 * z_3**2
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 19))
     1         + Ez_zc ( ixx3, iyy3, izz, 20)
     1         * (x_3 * y_3 * z_3
     1         - Coef_Basis_Moments ( ixx3, iyy3, izz, 20))

#endif

! --------------------
#else
! --------------------
! Choose arithmetically averaged electric fields for he centered part of the
! 2D LLF Riemann solver. Electric field reconstruction is not needed.

          Ex_1 = ( Ex_zc ( ixx1, iyy1, izz, 1)
     1           + Ex_zc ( ixx2, iyy2, izz, 1)
     1           + Ex_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ex_2 = Ex_1  ;  Ex_3 = Ex_1

          Ey_1 = ( Ey_zc ( ixx1, iyy1, izz, 1)
     1           + Ey_zc ( ixx2, iyy2, izz, 1)
     1           + Ey_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ey_2 = Ey_1  ;  Ey_3 = Ey_1

          Ez_1 = ( Ez_zc ( ixx1, iyy1, izz, 1)
     1           + Ez_zc ( ixx2, iyy2, izz, 1)
     1           + Ez_zc ( ixx3, iyy3, izz, 1) ) * ob3
          Ez_2 = Ez_1  ;  Ez_3 = Ez_1

! --------------------
#endif
! --------------------

! ------------------------------
! Build the electric fields in the local z-direction; which is the global r-direction.

          tempaa = Ex_1 * n3 ( 1) + Ey_1 * n3 ( 2) + Ez_1 * n3 ( 3)
          Ez_1 = tempaa

          tempaa = Ex_2 * n3 ( 1) + Ey_2 * n3 ( 2) + Ez_2 * n3 ( 3)
          Ez_2 = tempaa

          tempaa = Ex_3 * n3 ( 1) + Ey_3 * n3 ( 2) + Ez_3 * n3 ( 3)
          Ez_3 = tempaa

! ------------------------------
! Build "Bneta_1, Bneta_2, Bneta_3" They may not be oriented in the positive directions.
! "Bneta_1" is along "neta_1"; "Bneta_2" is along "neta_2"; "Bneta_3" is along "neta_3".

          Bneta_1 = - Byf_quadpts ( ixx3, iyy3, izz, 5)
     1  + Bx0 * neta_1 ( 1) + By0 * neta_1 ( 2) + Bz0 * neta_1 ( 3)

          Bneta_2 = Byf_quadpts ( ixx2, iyy2, izz, 4)
     1  + Bx0 * neta_2 ( 1) + By0 * neta_2 ( 2) + Bz0 * neta_2 ( 3)

          Bneta_3 = - Bxf_quadpts ( ixx2, iyy2, izz, 2)
     1  + Bx0 * neta_3 ( 1) + By0 * neta_3 ( 2) + Bz0 * neta_3 ( 3)

! ------------------------------
! Evaluate "S_R1, S_L1, S_R2, S_L2, S_R3, S_L3". Watch out for the orientatations
! because in each of the three zone faces we take the positive direction to be
! along "neta_1, neta_2, neta_3".

          S_R1 = - S_D_gr ( ixx3, iyy3, izz)
          S_L1 = - S_U_gr ( ixx3, iyy3, izz)

          S_R2 = S_U_gr ( ixx2, iyy2, izz)
          S_L2 = S_D_gr ( ixx2, iyy2, izz)

          S_R3 = - S_L_gr ( ixx2, iyy2, izz)
          S_L3 = - S_R_gr ( ixx2, iyy2, izz)

! --------------------
#if (USE_2D_LLF_RS_FOR_EZ != 0)
! --------------------
! Make wave speeds fully subsonic for the 2D RS.

          tempaa = AMAX1 ( ABS ( S_R1), ABS ( S_L1),
     1                     ABS ( S_R2), ABS ( S_L2),
     1                     ABS ( S_R3), ABS ( S_L3))

          S_R1 = tempaa  ;  S_L1 = - tempaa
          S_R2 = tempaa  ;  S_L2 = - tempaa
          S_R3 = tempaa  ;  S_L3 = - tempaa

! --------------------
#endif
! --------------------

! ------------------------------------------------------------
! We are only making z-directional calls to the 2D Riemann solver that has three
! input states.

          CALL Two_D_RS_Ez_3State_HLL

     1         ( S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,

     1           Ez_1, Ez_2, Ez_3,

     1           Bneta_1, Bneta_2, Bneta_3,

     1           BSSx, BSSy, ESSz)

! ------------------------------------------------------------
! Hand back the global electric field "ESSz" from the 2D Riemann solver.

          Efz ( ixx, iyy, izz) = Efz ( ixx, iyy, izz) + wt * ESSz

! ------------------------------------------------------------

        END DO  ! End of "DO izz = 1, nr_zones".

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, 3" loop over quadrature points.

      END IF  ! End of "IF ( itri_patch_LU .EQ. 1) THEN".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! End data for constants initialization.
!$acc end data

! ----------------------------------------------------------------------

      END SUBROUTINE MAKE_EDGE_EFIELDS_P1







