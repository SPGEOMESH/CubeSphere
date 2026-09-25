



      SUBROUTINE MAKE_FACIAL_FLUXES

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         ibdy_outer,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr,

     1         S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr,

     1         Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1         Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi, dU_dt,

     1         EntFlux_x, EntFlux_y, EntFlux_z,
     1         Efx, Efy, Efz)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes in all relevant geometric information, as described in the
! "LIMIT_ALL_SLOPES_P?" routines.

! Takes in the reconstructed fluid and magnetic field variables in
! "U, entropygr, B1_zc, B2_zc, B3_zc".

! Returns all fluid fluctuations and entropy fluxes.

! If "USE_FACIAL_EFIELDS == 1" then it also builds electric fields along
! the edges by averaging, "Efx, Efy, Efz".

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes,
     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      ibdy_outer

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

! -----

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
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost) ::
     1      gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1      divvelgr, msonicgr, flattengr

! ------------------------------
! FLUX AND ELECTRIC FIELD VARIABLES:-

! Store extremal speeds from the 1D Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost) ::
     1   S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr

! Fluid fluxes, fluctuations and resolved states from the 1D Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components) ::
     1      UStar_x, UStar_y, UStar_z

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                                          n_cc_components) ::
     1      Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1      Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost,
     1                                 n_cc_components) :: dU_dt

! Entropy flux from the 1D Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      EntFlux_x, EntFlux_y, EntFlux_z

! Electric fields at the edges of the mesh.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1       Efx, Efy, Efz

! ------------------------------
! This is just to prevent accumulation variables from being doubly added to.
! Needed when accumulating to the electric field.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1       Efx_m1,   Efy_m1,   Efz_m1

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost, 3) ::
     1       Efx_uvec, Efy_uvec, Efz_uvec

! counter for Efx, Efy, Efz accumulations from Fstar.
      INTEGER, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1       i_Efx, i_Efy, i_Efz

! ----------------------------------------------------------------------

! Variables past this point are internal to the subroutine.

      INTEGER ixx, iyy, izz, izzm1, izzp1,
     1        ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd, bflag

      REAL x, y, z, x2, y2, z2, obl0, obl02, time_use,
     1     tempaa, tempab, tempac, v1 ( 3), v2 ( 3),
     1     Efld ( 3), uvec ( 3)

      REAL Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz

! ----------------------------------------------------------------------

! These variables go through the Riemann solver interface. They can be declared
! SHARED and then the dime-a-dozen 1D HLL Riemann solver can be inlined.

      REAL fourpi, ob_4pi, ob_8pi, flat

      REAL S_L, S_R, ent_L, ent_R, ent_F_Star

      REAL, DIMENSION ( NVAR + NFLUID) ::
     1      U_L, U_R, U_Star, F_Star, V_L, V_R, F_L, F_R, F_LL, F_RR,

     1      Fluct_L, Fluct_R, Fluct_Pl, Fluct_Mi,
     1      V_Star

! ----------

      INTEGER iexpand, j

      REAL prs_L, prs_L0, del_prs_L, prs_R, prs_R0, del_prs_R,
     1     temp_01, temp_02, ob_rho,
     1     P_tot, Bf_2, Vel_2, V_dot_B, B0_dot_B,

     1     rhobar, vlxbar, vlybar, vlzbar, prsbar,
     1     bfxre, bfybar, bfzbar,

     1     va_x, ob_4pi_rho, c_s, v_a, v_f, wt

      REAL, DIMENSION ( 3) :: B0

! ----------

      REAL, DIMENSION ( n_cc_components) :: V_zc,
     1        dV_zc_dx, dV_zc_dy, dV_zc_dz

      REAL B0x, B0y, B0z,
     1     dB0x_dx, dB0x_dy, dB0x_dz, dB0y_dx, dB0y_dy, dB0y_dz,
     1     dB0z_dx, dB0z_dy, dB0z_dz,
     1     nc_product ( 2)

! ----------------------------------------------------------------------
! Index varibles, which helps in reconstruction using "Coef_Basis".
! It is mainly useful for the reconstruction of B1, B2, B3 where the
! sequence changes for higher order modes.

      REAL Coef_Basis ( n_cc_modes)

      INTEGER, PARAMETER ::
! 0, 1st, 2nd order modes.
     1      iC0  = 1,
     1      iCx  = 2,  iCy  = 3,  iCz  = 4

      INTEGER iOriginal_RS

      REAL c_sL, v_anL, v_aL, d_L2, v_fnL,
     1     c_sR, v_anR, v_aR, d_R2, v_fnR,
     1     v_atL, v_snL, alpha, beta, theta, c_ss,
     1     ent_temp, prs_avg

      INTEGER iprs_flag

! ----------------------------------------------------------------------

!$acc routine( OneD_RS_MHD_Ptwise ) seq
!$acc routine( GET_MAGNETIC_FIELD ) seq
!$acc routine( GET_MAGNETIC_FIELD_Gradients ) seq
!$acc routine( MAKE_NC_CHARMAT_Prod_Prim ) seq

! ----------------------------------------------------------------------

      fourpi = 4.0 * pi

      ob_8pi = 1.0 / ( 8.0 * pi)
      ob_4pi = 1.0 / ( 4.0 * pi)

! Minimum expansion is by 1 so that we have S_R, S_L etc. at zone boundaries
! for at least one layer of ghost zones. These are needed for electric
! field evaluation.

      iexpand = 1

! ------------------------------
! Zero out all variables that will be used as accumulation variables in this
! subroutine.

!$acc kernels

      dU_dt ( :, :, :, :) = 0.0

      Fluct_x_Pl ( :, :, :, :) = 0.0
      Fluct_x_Mi ( :, :, :, :) = 0.0

      Fluct_y_Pl ( :, :, :, :) = 0.0
      Fluct_y_Mi ( :, :, :, :) = 0.0

      Fluct_z_Pl ( :, :, :, :) = 0.0
      Fluct_z_Mi ( :, :, :, :) = 0.0

!$acc end kernels


! ------------------------------

! Move the initialized constants over to the GPU.

!$acc data copyin ( fourpi, ob_8pi, ob_4pi, iexpand,                    &
!$acc&              Efx_m1, Efy_m1, Efz_m1, i_Efx, i_Efy, i_Efz,        &
!$acc&              Efx_uvec, Efy_uvec, Efz_uvec,                       &
!$acc&              UStar_x, UStar_y, UStar_z)

!$acc update device ( itri_patch_LD, itri_patch_RD, itri_patch_RU,      &
!$acc&                itri_patch_LU, ibdy_outer

! ----------------------------------------------------------------------
! Zero out all the variables that will be accumulation variables in this
! subroutine. On a GPU this should be strongly parallelized.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, v1, v2, uvec)                &
!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    smallnum,                                                     &
!$OMP&    vert_coords_patch, Efx, Efy, Efz, i_Efx, i_Efy, i_Efz,        &
!$OMP&    Efx_uvec, Efy_uvec, Efz_uvec, dU_dt)                          &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum,                                                     &
!$acc&    vert_coords_patch, Efx, Efy, Efz, i_Efx, i_Efy, i_Efz,        &
!$acc&    Efx_uvec, Efy_uvec, Efz_uvec, dU_dt)
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, v1, v2, uvec)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

! ------------------------------
! Initialize the variables to zero.

          Efx ( ixx, iyy, izz) = 0.0
          Efy ( ixx, iyy, izz) = 0.0
          Efz ( ixx, iyy, izz) = 0.0

          dU_dt ( ixx, iyy, izz, :) = 0.0

! Initialize the counter for accumalation of "Ef_" to zero.
          i_Efx ( ixx, iyy, izz) = 0
          i_Efy ( ixx, iyy, izz) = 0
          i_Efz ( ixx, iyy, izz) = 0

! ------------------------------
! Compute and store the unit vectors of Efx, Efy, Efz

! --------------------
! Efx-unit vectors.

      v1 ( 1:3) = vert_coords_patch ( ixx-1, iyy, izz, :)
      v2 ( 1:3) = vert_coords_patch ( ixx,   iyy, izz, :)

      uvec ( :) = v2 ( :) - v1 ( :)      ! v1 --> v2

      Efx_uvec ( ixx, iyy, izz, :) = uvec ( :)
     1                 / MAX ( SQRT ( uvec ( 1)**2 + uvec ( 2)**2
     1                              + uvec ( 3)**2), smallnum)

! --------------------
! Efy-unit vectors.

      v1 ( 1:3) = vert_coords_patch ( ixx, iyy-1, izz, :)
      v2 ( 1:3) = vert_coords_patch ( ixx, iyy,   izz, :)

      uvec ( :) = v2 ( :) - v1 ( :)      ! v1 --> v2

      Efy_uvec ( ixx, iyy, izz, :) = uvec ( :)
     1                 / MAX ( SQRT ( uvec ( 1)**2 + uvec ( 2)**2
     1                              + uvec ( 3)**2), smallnum)

! --------------------
! Efz-unit vectors.

      v1 ( 1:3) = vert_coords_patch ( ixx, iyy, izz-1, :)
      v2 ( 1:3) = vert_coords_patch ( ixx, iyy, izz,   :)

      uvec ( :) = v2 ( :) - v1 ( :)      ! v1 --> v2

      Efz_uvec ( ixx, iyy, izz, :) = uvec ( :)
     1                 / MAX ( SQRT ( uvec ( 1)**2 + uvec ( 2)**2
     1                              + uvec ( 3)**2), smallnum)

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 1) Access the x-faces. From either side of the x-face we use
! "U ( ixxl, iyyl, izz, :, 1:4)" and "U ( ixxr, iyyr, izz, :, 1:4)" to obtain
! "U_L ( :), U_R ( :)" which are fed into the 1D Riemann solver.

! The higher order modes from "B1_zc, B2_zc, B3_zc" on either side of the
! x-face are also used for the magnetic field.

! Since this is only second order, we obtain "U_L ( :), U_R ( :)"
! at the centroid of the x-face.

! If "PRESSURE_POSITIVITY >= 1" we also use "entropygr" to obtain
! "ent_L, ent_R"
! ----------------------------------------------------------------------

!$acc kernels
      Efx_m1 ( :, :, :) = 0.0
      Efy_m1 ( :, :, :) = 0.0
      Efz_m1 ( :, :, :) = 0.0
!$acc end kernels

! Short loop over all quadrature points. Do not parallelize this outer loop.
! "j = 1" is the only central point and it is the only quadrature point needed at
! second order.
      DO j = 1, NAREA_QUAD_PTS

#if (NUMORDER == 2)
      IF ( j .GT. 1) EXIT
#endif


!$acc update device ( iexpand, j)

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxl, iyyl, ixxr, iyyr,      &
!$OMP&    x, y, z, x2, y2, z2, obl0, obl02, v1, Efld,                   &
!$OMP&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, tempaa, tempab, tempac,   &
!$OMP&    bflag, time_use, B0, Fluct_Pl, Fluct_Mi, Coef_Basis,          &

!$OMP&    S_L, S_R, ent_L, ent_R, ent_F_Star, Fluct_L, Fluct_R, flat,   &
!$OMP&    U_L, U_R, U_Star, F_Star, V_L, V_R, F_L, F_R, F_LL, F_RR,     &
!$OMP&    prs_L, prs_L0, del_prs_L, prs_R, prs_R0, del_prs_R,           &
!$OMP&    temp_01, temp_02, ob_rho, P_tot, Bf_2, Vel_2, V_dot_B,        &
!$OMP&    rhobar, vlxbar, vlybar, vlzbar, prsbar, bfxre, bfybar, bfzbar,&
!$OMP&    va_x, ob_4pi_rho, c_s, v_a, v_f, wt )

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    fourpi, ob_4pi, ob_8pi, smallnum, iexpand, j,                 &
!$OMP&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$OMP&    molewt, gasconst, rhofloor, prsfloor, tempfloor,              &
!$OMP&    imask_active_patch, istencil_patch, v_scale_length,           &
!$OMP&    xf_centroid_patch, xf_quadnorml_patch, Coef_Basis_Moments,    &
!$OMP&    xf_quadpts_patch, xf_quadwt_patch,                            &
!$OMP&    v_centroid_patch, xf_area_patch, v_volume_patch,              &
!$OMP&    U, V, entropygr, B1_zc, B2_zc, B3_zc, flattengr,              &
!$OMP&    UStar_x, EntFlux_x, S_L_gr, S_R_gr,                           &
!$OMP&    Fluct_x_Pl, Fluct_x_Mi, dU_dt,                                &
!$OMP&    magB_small_cgl, eps_bdry_cgl,                                 &
!$OMP&    Efy, Efz, Efy_m1, Efz_m1, i_Efy, i_Efz, Efy_uvec, Efz_uvec,   &
!$OMP&    gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr)

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    fourpi, ob_4pi, ob_8pi, smallnum, iexpand, j,                 &
!$acc&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$acc&    molewt, gasconst, rhofloor, prsfloor, tempfloor,              &
!$acc&    imask_active_patch, istencil_patch, v_scale_length,           &
!$acc&    xf_centroid_patch, xf_quadnorml_patch, Coef_Basis_Moments,    &
!$acc&    xf_quadpts_patch, xf_quadwt_patch,                            &
!$acc&    v_centroid_patch, xf_area_patch, v_volume_patch,              &
!$acc&    U, V, entropygr, B1_zc, B2_zc, B3_zc, flattengr,              &
!$acc&    UStar_x, EntFlux_x, S_L_gr, S_R_gr,                           &
!$acc&    Fluct_x_Pl, Fluct_x_Mi, dU_dt,                                &
!$acc&    magB_small_cgl, eps_bdry_cgl,                                 &
!$acc&    Efy, Efz, Efy_m1, Efz_m1, i_Efy, i_Efz, Efy_uvec, Efz_uvec,   &
!$acc&    gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxl, iyyl, ixxr, iyyr,              &
!$acc&    x, y, z, x2, y2, z2, obl0, obl02, v1, Efld,                   &
!$acc&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, tempaa, tempab, tempac,   &
!$acc&    bflag, time_use, B0, Fluct_Pl, Fluct_Mi, Coef_Basis,          &

!$acc&    S_L, S_R, ent_L, ent_R, ent_F_Star, Fluct_L, Fluct_R, flat,   &
!$acc&    U_L, U_R, U_Star, F_Star, V_L, V_R, F_L, F_R, F_LL, F_RR,     &
!$acc&    prs_L, prs_L0, del_prs_L, prs_R, prs_R0, del_prs_R,           &
!$acc&    temp_01, temp_02, ob_rho, P_tot, Bf_2, Vel_2, V_dot_B,        &
!$acc&    rhobar, vlxbar, vlybar, vlzbar, prsbar, bfxre, bfybar, bfzbar,&
!$acc&    va_x, ob_4pi_rho, c_s, v_a, v_f, wt )

      DO izz = 1 - iexpand, nr_zones + iexpand

        DO iyy = 1 - iexpand, ntp_zones + iexpand

          DO ixx = 0 - iexpand, ntp_zones + iexpand

! ------------------------------------------------------------

          IF ( ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) .AND.
     1         ( imask_active_patch ( ixx + 1, iyy, izz) .EQ. 1) ) THEN

            ixxl = ixx
            iyyl = iyy

            ixxr = ixx + 1
            iyyr = iyy

          ELSE IF ( ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) .AND.
     1         ( imask_active_patch ( ixx + 1, iyy, izz) .EQ. 1) ) THEN

! Physical zone to the Rt of x-face.
            ixxr = ixx + 1
            iyyr = iyy

! Physical zone to the Lf of x-face.
            ixxl = istencil_patch ( ixxr, iyyr, 1, 2)
            iyyl = istencil_patch ( ixxr, iyyr, 2, 2)

          ELSE IF ( ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) .AND.
     1         ( imask_active_patch ( ixx + 1, iyy, izz) .EQ. 0) ) THEN

! Physical zone to the Lf of x-face.
            ixxl = ixx
            iyyl = iyy

! Physical zone to the Rt of x-face.
            ixxr = istencil_patch ( ixxl, iyyl, 1, 1)
            iyyr = istencil_patch ( ixxl, iyyl, 2, 1)

          ELSE

            V_L ( :) = 0.0
            V_R ( :) = 0.0

            ent_L = prsfloor / rhofloor**(gamma - 1.0)
            ent_R = ent_L

            CYCLE  ! Go to the next iteration of the "DO ixx =" loop.

          END IF

! --------------------------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixxl, iyyl, izz), smallnum)  ! Scaling for Lf zone.

          obl02 = 1.0
     1  / AMAX1 ( v_scale_length ( ixxr, iyyr, izz), smallnum)  ! Scaling for Rt zone.


          Coef_Basis ( :)  = Coef_Basis_Moments ( ixx, iyy, izz, :)

! ------------------------------
! Obtain "Ax, Ay, Az" the normal to this face. This is the local x-direction.
! The entire x-face has the same normal. "Bx, By, Bz" is in the radial direction.
! This is also the local y-direction. The local z-direction is given by the
! vector cross product "C = A X B".

! In this loop, we rotate vectors into the local xyz-directions. Once we
! solve the Riemann problem, we will rotate flux components back into the global
! xyz-directions. For this reason, we save "Ax, Ay, Az" etc. as persistent 1d arrays.


          Ax = xf_quadnorml_patch ( ixx, iyy, izz, 1)  ! Unit vector in Local x-direction.
          Ay = xf_quadnorml_patch ( ixx, iyy, izz, 2)
          Az = xf_quadnorml_patch ( ixx, iyy, izz, 3)

          Cx = xf_centroid_patch ( ixx, iyy, izz, 1)
          Cy = xf_centroid_patch ( ixx, iyy, izz, 2)
          Cz = xf_centroid_patch ( ixx, iyy, izz, 3)

          tempaa = AMAX1 ( SQRT ( Cx**2 + Cy**2 + Cz**2), smallnum)
          tempaa = 1.0 / tempaa

          Cx = Cx * tempaa  ! Unit vector in Local z-direction.
          Cy = Cy * tempaa
          Cz = Cz * tempaa


          Bx = Cy * Az - Cz * Ay  ! Unit vector in Local y-direction.
          By = Cz * Ax - Cx * Az
          Bz = Cx * Ay - Cy * Ax

#if( 0 == 1)
! Validate C is unit vector. In turn validates A, B are perfectly normal.
          IF ( ABS ( SQRT ( Bx**2 + By**2 + Bz**2)
     1             - 1.0 ) .GE. smallnum ) THEN
            WRITE ( 6, *)"Stopping at ",__LINE__, "in ",__FILE__
            CALL FLUSH ( 6)
            STOP
          END IF
#endif

! ------------------------------
! Build scaled relative coordinates from Lf zone centroid and Rt zone centroid to
! facial centroid. Only one Riemann solve at second order. (BTW, at higher orders
! the facial centroid would be replaced by looping over facial quadrature points.)

          x = ( xf_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixxl, iyyl, izz, 1) ) * obl0

          y = ( xf_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixxl, iyyl, izz, 2) ) * obl0

          z = ( xf_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixxl, iyyl, izz, 3) ) * obl0

! ----------

          x2 = ( xf_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixxr, iyyr, izz, 1) ) * obl02

          y2 = ( xf_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixxr, iyyr, izz, 2) ) * obl02

          z2 = ( xf_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixxr, iyyr, izz, 3) ) * obl02

! ----------
! Obtain the weight for this quadrature point. Weight is 1 at second order.

          wt = 1.0
#if (NUMORDER >= 3)
          wt = xf_quadwt_patch ( ixx, iyy, izz, j)
#endif

! ------------------------------
! Now use the above relative coordinates as well as
! "U ( ixxl, iyyl, izz, :, 1:4)" and "U ( ixxr, iyyr, izz, :, 1:4)"
! to build "U_L ( :), U_R ( :)" with second order accuracy.
! Note the use of "x, y, z" from the left and "x2, y2, z2" from the right.

          V_L ( :) = V ( ixxl, iyyl, izz, :, 1)
     1             + V ( ixxl, iyyl, izz, :, 2) * x
     1             + V ( ixxl, iyyl, izz, :, 3) * y
     1             + V ( ixxl, iyyl, izz, :, 4) * z
#if (NUMORDER >= 3)
     1             + V ( ixxl, iyyl, izz, :, 5)
     1             * (x**2 - Coef_Basis_Moments ( ixxl, iyyl, izz, 5))
     1             + V ( ixxl, iyyl, izz, :, 6)
     1             * (y**2 - Coef_Basis_Moments ( ixxl, iyyl, izz, 6))
     1             + V ( ixxl, iyyl, izz, :, 7)
     1             * (z**2 - Coef_Basis_Moments ( ixxl, iyyl, izz, 7))
     1             + V ( ixxl, iyyl, izz, :, 8)
     1             * (x * y - Coef_Basis_Moments ( ixxl, iyyl, izz, 8))
     1             + V ( ixxl, iyyl, izz, :, 9)
     1             * (y * z - Coef_Basis_Moments ( ixxl, iyyl, izz, 9))
     1             + V ( ixxl, iyyl, izz, :, 10)
     1             * (x * z - Coef_Basis_Moments ( ixxl, iyyl, izz, 10))
#elif (NUMORDER >= 4)
     1             + V ( ixxl, iyyl, izz, :, 11)
     1             * (x**3 - Coef_Basis_Moments ( ixxl, iyyl, izz, 11))
     1             + V ( ixxl, iyyl, izz, :, 12)
     1             * (y**3 - Coef_Basis_Moments ( ixxl, iyyl, izz, 12))
     1             + V ( ixxl, iyyl, izz, :, 13)
     1             * (z**3 - Coef_Basis_Moments ( ixxl, iyyl, izz, 13))
     1             + V ( ixxl, iyyl, izz, :, 14)
     1             * (x**2 * y
     1             - Coef_Basis_Moments ( ixxl, iyyl, izz, 14))
     1             + V ( ixxl, iyyl, izz, :, 15)
     1             * (x**2 * z
     1             - Coef_Basis_Moments ( ixxl, iyyl, izz, 15))
     1             + V ( ixxl, iyyl, izz, :, 16)
     1             * (x * y**2
     1             - Coef_Basis_Moments ( ixxl, iyyl, izz, 16))
     1             + V ( ixxl, iyyl, izz, :, 17)
     1             * (y**2 * z
     1             - Coef_Basis_Moments ( ixxl, iyyl, izz, 17))
     1             + V ( ixxl, iyyl, izz, :, 18)
     1             * (x * z**2
     1             - Coef_Basis_Moments ( ixxl, iyyl, izz, 18))
     1             + V ( ixxl, iyyl, izz, :, 19)
     1             * (y * z**2
     1             - Coef_Basis_Moments ( ixxl, iyyl, izz, 19))
     1             + V ( ixxl, iyyl, izz, :, 20)
     1             * (x * y * z
     1             - Coef_Basis_Moments ( ixxl, iyyl, izz, 20))
#endif

          V_R ( :) = V ( ixxr, iyyr, izz, :, 1)
     1             + V ( ixxr, iyyr, izz, :, 2) * x2
     1             + V ( ixxr, iyyr, izz, :, 3) * y2
     1             + V ( ixxr, iyyr, izz, :, 4) * z2
#if (NUMORDER >= 3)
     1             + V ( ixxr, iyyr, izz, :, 5)
     1             * (x2**2 - Coef_Basis_Moments ( ixxr, iyyr, izz, 5))
     1             + V ( ixxr, iyyr, izz, :, 6)
     1             * (y2**2 - Coef_Basis_Moments ( ixxr, iyyr, izz, 6))
     1             + V ( ixxr, iyyr, izz, :, 7)
     1             * (z2**2 - Coef_Basis_Moments ( ixxr, iyyr, izz, 7))
     1             + V ( ixxr, iyyr, izz, :, 8)
     1             * (x2 * y2
     1             - Coef_Basis_Moments ( ixxr, iyyr, izz, 8))
     1             + V ( ixxr, iyyr, izz, :, 9)
     1             * (y2 * z2
     1             - Coef_Basis_Moments ( ixxr, iyyr, izz, 9))
     1             + V ( ixxr, iyyr, izz, :, 10)
     1             * (x2 * z2
     1             - Coef_Basis_Moments ( ixxr, iyyr, izz, 10))
#elif (NUMORDER >= 4)
     1             + V ( ixxr, iyyr, izz, :, 11)
     1             * (x2**3 - Coef_Basis_Moments ( ixxr, iyyr, izz, 11))
     1             + V ( ixxr, iyyr, izz, :, 12)
     1             * (y2**3 - Coef_Basis_Moments ( ixxr, iyyr, izz, 12))
     1             + V ( ixxr, iyyr, izz, :, 13)
     1             * (z2**3 - Coef_Basis_Moments ( ixxr, iyyr, izz, 13))
     1             + V ( ixxr, iyyr, izz, :, 14)
     1             * (x2**2 * y2
     1             - Coef_Basis_Moments ( ixxr, iyyr, izz, 14))
     1             + V ( ixxr, iyyr, izz, :, 15)
     1             * (x2**2 * z2
     1             - Coef_Basis_Moments ( ixxr, iyyr, izz, 15))
     1             + V ( ixxr, iyyr, izz, :, 16)
     1             * (x2 * y2**2
     1             - Coef_Basis_Moments ( ixxr, iyyr, izz, 16))
     1             + V ( ixxr, iyyr, izz, :, 17)
     1             * (y2**2 * z2
     1             - Coef_Basis_Moments ( ixxr, iyyr, izz, 17))
     1             + V ( ixxr, iyyr, izz, :, 18)
     1             * (x2 * z2**2
     1             - Coef_Basis_Moments ( ixxr, iyyr, izz, 18))
     1             + V ( ixxr, iyyr, izz, :, 19)
     1             * (y2 * z2**2
     1             - Coef_Basis_Moments ( ixxr, iyyr, izz, 19))
     1             + V ( ixxr, iyyr, izz, :, 20)
     1             * (x2 * y2 * z2
     1             - Coef_Basis_Moments ( ixxr, iyyr, izz, 20))
#endif

! ------------------------------
! Now use the above relative coordinates as well as
! "entropygr ( ixxl, iyyl, izz, 1:4)" and "entropygr ( ixxr, iyyr, izz, 1:4)"
! to build "ent_L, ent_R" with second order accuracy.

#if ( PRESSURE_POSITIVITY >= 1)

          ent_L = entropygr ( ixxl, iyyl, izz, 1)
     1                + entropygr ( ixxl, iyyl, izz, 2) * x
     1                + entropygr ( ixxl, iyyl, izz, 3) * y
     1                + entropygr ( ixxl, iyyl, izz, 4) * z

          ent_R = entropygr ( ixxr, iyyr, izz, 1)
     1                + entropygr ( ixxr, iyyr, izz, 2) * x2
     1                + entropygr ( ixxr, iyyr, izz, 3) * y2
     1                + entropygr ( ixxr, iyyr, izz, 4) * z2

#else
          ent_L = prsfloor / rhofloor**(gamma - 1.0)
          ent_R = ent_L
#endif

! ------------------------------
! Notice that the 1D Riemann solver always assumes that the waves expand
! in the local x-direction. Therefore, we have to rotate the vectorial entities
! like velocities and magnetic fields into the local coordinate system.

          tempaa = Ax * V_L ( 2) + Ay * V_L ( 3) + Az * V_L ( 4)
          tempab = Bx * V_L ( 2) + By * V_L ( 3) + Bz * V_L ( 4)
          tempac = Cx * V_L ( 2) + Cy * V_L ( 3) + Cz * V_L ( 4)

          V_L ( 2) = tempaa
          V_L ( 3) = tempab
          V_L ( 4) = tempac

          tempaa = Ax * V_L ( 6) + Ay * V_L ( 7) + Az * V_L ( 8)
          tempab = Bx * V_L ( 6) + By * V_L ( 7) + Bz * V_L ( 8)
          tempac = Cx * V_L ( 6) + Cy * V_L ( 7) + Cz * V_L ( 8)

          V_L ( 6) = tempaa
          V_L ( 7) = tempab
          V_L ( 8) = tempac

! ----------

          tempaa = Ax * V_R ( 2) + Ay * V_R ( 3) + Az * V_R ( 4)
          tempab = Bx * V_R ( 2) + By * V_R ( 3) + Bz * V_R ( 4)
          tempac = Cx * V_R ( 2) + Cy * V_R ( 3) + Cz * V_R ( 4)

          V_R ( 2) = tempaa
          V_R ( 3) = tempab
          V_R ( 4) = tempac

          tempaa = Ax * V_R ( 6) + Ay * V_R ( 7) + Az * V_R ( 8)
          tempab = Bx * V_R ( 6) + By * V_R ( 7) + Bz * V_R ( 8)
          tempac = Cx * V_R ( 6) + Cy * V_R ( 7) + Cz * V_R ( 8)

          V_R ( 6) = tempaa
          V_R ( 7) = tempab
          V_R ( 8) = tempac

! --------------------------------------------------
! Obtain the background magnetic field B0 and rotate it.

          B0 ( :) = 0.0

#if ( USE_B0B1_FORMULATION == 1)

          x = xf_centroid_patch ( ixx, iyy, izz, 1)
          y = xf_centroid_patch ( ixx, iyy, izz, 2)
          z = xf_centroid_patch ( ixx, iyy, izz, 3)

          bflag = 1
          time_use = 0

! Calling acc-seq function to get the B0, using temp. variables x2, y2, z2.
          CALL GET_MAGNETIC_FIELD ( bflag, time_use,  x, y, z,
     1                              x2, y2, z2 )

          B0 ( 1) = Ax * x2 + Ay * y2 + Az * z2
          B0 ( 2) = Bx * x2 + By * y2 + Bz * z2
          B0 ( 3) = Cx * x2 + Cy * y2 + Cz * z2

#endif

! --------------------------------------------------

          flat = AMAX1 ( flattengr ( ixxl, iyyl, izz),
     1                   flattengr ( ixxr, iyyr, izz) )

! --------------------------------------------------
! Now call the 1D Riemann solver. This RS operates in the local x-direction.
! It takes in "V_L, V_R" and "ent_L, ent_R" and returns speeds "S_L, S_R".
! Also returns resolved state "U_Star", resolved flux "F_Star", entropy flux "ent_F_Star".
! Also returns fluctuation "Fluct_L, Fluct_R" and "U_L, U_R, F_L, F_R".

          CALL OneD_RS_MHD_Ptwise

     1         ( gamma, gamma_euler, isospeed, barotropic_stiffness,
     1           molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1           magB_small_cgl, eps_bdry_cgl,

     1           fourpi, ob_4pi, ob_8pi,
     1           smallnum, B0, flat,

     1           V_L, V_R, ent_L, ent_R,
     1           S_L, S_R, U_Star, F_Star, ent_F_Star,
     1           Fluct_L, Fluct_R,
     1           U_L, U_R, F_L, F_R)

! ------------------------------------------------------------
! Hand back the resolved states, fluxes and extremal wave speeds from the
! 1D Riemann solver.

          IF ( ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) .AND.
     1         ( imask_active_patch ( ixx + 1, iyy, izz) .EQ. 0) ) THEN
! If we get here, we are at an unphysical boundary, so set zero values and move on.

            UStar_x ( ixx, iyy, izz, :) = 0.0
            EntFlux_x ( ixx, iyy, izz) = 0.0

            CYCLE  ! Go to the next iteration of the "DO ixx =" loop.

          END IF

! ------------------------------
! If we get here, we are at a physical boundary, so assign the good values
! from the 1D dime-a-dozen Riemann solver.

! ----------
! If "USE_FACIAL_EFIELDS == 1" then "Ez = -F7" and "Ey = F8" in the local
! coordinates. Therefore, before we rotate back, we assign "Efy, Efz".

#if ( USE_FACIAL_EFIELDS == 1)
          IF ( j .EQ. 1) THEN

          v1 ( 1)   = 0.0
          v1 ( 2)   =   F_Star ( 8)
          v1 ( 3)   = - F_Star ( 7)

! Accumulate the upper edge.
          Efy ( ixx, iyy, izz) = Efy ( ixx, iyy, izz) + v1 ( 2)
          Efz ( ixx, iyy, izz) = Efz ( ixx, iyy, izz) + v1 ( 3)

! Count the number of accumulations.
          i_Efz ( ixx, iyy, izz) = i_Efz ( ixx, iyy, izz) + 1
          i_Efy ( ixx, iyy, izz) = i_Efy ( ixx, iyy, izz) + 1

! Assign directly for the lower edge.
          Efy_m1 ( ixx, iyy, izz - 1) = v1 ( 2)
          Efz_m1 ( ixx, iyy - 1, izz) = v1 ( 3)

          END IF
#endif

! ----------
! Save the extremal speeds. They will be needed in the 2D RS for electric field evaluation.
! This is only done for the face-centered quadrature point, because it has the best speeds.

          IF ( j .EQ. 1) THEN

          S_L_gr ( ixx, iyy, izz) = S_L

          S_R_gr ( ixx, iyy, izz) = S_R

          END IF

! ----------
! For CGL we have to use fluctuation form. For that, we return "Fluct_L, Fluct_R".
! These will also have to undergo the vectorial components in the same way.

! ----------
! Fluctuations in + direction; rotated into global frame.

          Fluct_Pl ( 1) = Fluct_R ( 1)                        ! Density

          Fluct_Pl ( 2) = Ax * Fluct_R ( 2)                   ! x-Momentum
     1               + Bx * Fluct_R ( 3) + Cx * Fluct_R ( 4)
          Fluct_Pl ( 3) = Ay * Fluct_R ( 2)                   ! y-Momentum
     1               + By * Fluct_R ( 3) + Cy * Fluct_R ( 4)
          Fluct_Pl ( 4) = Az * Fluct_R ( 2)                   ! z-Momentum
     1               + Bz * Fluct_R ( 3) + Cz * Fluct_R ( 4)

          Fluct_Pl ( 5) = Fluct_R ( 5)                        ! Energy

          Fluct_Pl ( 6) = Ax * Fluct_R ( 6)                   ! x-Magnetic field
     1               + Bx * Fluct_R ( 7) + Cx * Fluct_R ( 8)
          Fluct_Pl ( 7) = Ay * Fluct_R ( 6)                   ! y-Magnetic field
     1               + By * Fluct_R ( 7) + Cy * Fluct_R ( 8)
          Fluct_Pl ( 8) = Az * Fluct_R ( 6)                   ! z-Magnetic field
     1               + Bz * Fluct_R ( 7) + Cz * Fluct_R ( 8)

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3)
          Fluct_Pl ( 9: ) = Fluct_R ( 9: )  ! pressures and tracer fields
#endif

#if ( IPDE_SYSTEM == 1 && NFLUID >= 1)
          Fluct_Pl ( NVAR + 1: NVAR + NFLUID) =  ! Species
     1                Fluct_R ( NVAR + 1: NVAR + NFLUID)
#endif

! -----

          Fluct_x_Pl ( ixx, iyy, izz, :) =Fluct_x_Pl ( ixx, iyy, izz, :)
     1    + wt * Fluct_Pl ( :)

! ----------
! Fluctuations in - direction; rotated into global frame.

          Fluct_Mi ( 1) = Fluct_L ( 1)                        ! Density

          Fluct_Mi ( 2) = Ax * Fluct_L ( 2)                   ! x-Momentum
     1               + Bx * Fluct_L ( 3) + Cx * Fluct_L ( 4)
          Fluct_Mi ( 3) = Ay * Fluct_L ( 2)                   ! y-Momentum
     1               + By * Fluct_L ( 3) + Cy * Fluct_L ( 4)
          Fluct_Mi ( 4) = Az * Fluct_L ( 2)                   ! z-Momentum
     1               + Bz * Fluct_L ( 3) + Cz * Fluct_L ( 4)

          Fluct_Mi ( 5) = Fluct_L ( 5)                        ! Energy

          Fluct_Mi ( 6) = Ax * Fluct_L ( 6)                   ! x-Magnetic field
     1               + Bx * Fluct_L ( 7) + Cx * Fluct_L ( 8)
          Fluct_Mi ( 7) = Ay * Fluct_L ( 6)                   ! y-Magnetic field
     1               + By * Fluct_L ( 7) + Cy * Fluct_L ( 8)
          Fluct_Mi ( 8) = Az * Fluct_L ( 6)                   ! z-Magnetic field
     1               + Bz * Fluct_L ( 7) + Cz * Fluct_L ( 8)

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3)
          Fluct_Mi ( 9: ) = Fluct_L ( 9: )  ! pressures and tracer fields
#endif

#if ( IPDE_SYSTEM == 1 && NFLUID >= 1)
          Fluct_Mi ( NVAR + 1: NVAR + NFLUID) =  ! Species
     1                Fluct_L ( NVAR + 1: NVAR + NFLUID)
#endif

! -----

          Fluct_x_Mi ( ixx, iyy, izz, :) =Fluct_x_Mi ( ixx, iyy, izz, :)
     1    + wt * Fluct_Mi ( :)

! ----------
! Rotate F_L to F_LL; rotate F_R to F_RR and add in the internal terms to the update.
! This is only done for the momentum fluxes and the magnetic field fluxes.

          F_LL = F_L

          F_LL ( 2) = Ax * F_L ( 2) + Bx * F_L ( 3) + Cx * F_L ( 4)
          F_LL ( 3) = Ay * F_L ( 2) + By * F_L ( 3) + Cy * F_L ( 4)
          F_LL ( 4) = Az * F_L ( 2) + Bz * F_L ( 3) + Cz * F_L ( 4)

          F_LL ( 6) = Ax * F_L ( 6) + Bx * F_L ( 7) + Cx * F_L ( 8)
          F_LL ( 7) = Ay * F_L ( 6) + By * F_L ( 7) + Cy * F_L ( 8)
          F_LL ( 8) = Az * F_L ( 6) + Bz * F_L ( 7) + Cz * F_L ( 8)

! -----

          F_RR = F_R

          F_RR ( 2) = Ax * F_R ( 2) + Bx * F_R ( 3) + Cx * F_R ( 4)
          F_RR ( 3) = Ay * F_R ( 2) + By * F_R ( 3) + Cy * F_R ( 4)
          F_RR ( 4) = Az * F_R ( 2) + Bz * F_R ( 3) + Cz * F_R ( 4)

          F_RR ( 6) = Ax * F_R ( 6) + Bx * F_R ( 7) + Cx * F_R ( 8)
          F_RR ( 7) = Ay * F_R ( 6) + By * F_R ( 7) + Cy * F_R ( 8)
          F_RR ( 8) = Az * F_R ( 6) + Bz * F_R ( 7) + Cz * F_R ( 8)

! -----
! Now add in the contributions from the internal fluxes to the update terms in
! the zone "( ixxl, iyyl, izz, :)" and "( ixxr, iyyr, izz, :)".

          dU_dt ( ixxl, iyyl, izz, :) = dU_dt ( ixxl, iyyl, izz, :)
     1  - wt * ( Fluct_Mi ( :) + F_LL ( :) )
     1       * xf_area_patch ( ixx, iyy, izz)
     1       / v_volume_patch ( ixxl, iyyl, izz)

! -----

          dU_dt ( ixxr, iyyr, izz, :) = dU_dt ( ixxr, iyyr, izz, :)
     1  - wt * ( Fluct_Pl ( :) - F_RR ( :) )
     1       * xf_area_patch ( ixx, iyy, izz)
     1       / v_volume_patch ( ixxr, iyyr, izz)

! ----------
! Likewise, rotate the state components back into the global xyz-directions.

          IF ( j .EQ. 1) THEN

          UStar_x ( ixx, iyy, izz, 1) = U_Star ( 1)  ! Density


          UStar_x ( ixx, iyy, izz, 2) = Ax * U_Star ( 2)  ! Momentum density
     1             + Bx * U_Star ( 3) + Cx * U_Star ( 4)

          UStar_x ( ixx, iyy, izz, 3) = Ay * U_Star ( 2)
     1             + By * U_Star ( 3) + Cy * U_Star ( 4)

          UStar_x ( ixx, iyy, izz, 4) = Az * U_Star ( 2)
     1             + Bz * U_Star ( 3) + Cz * U_Star ( 4)


          UStar_x ( ixx, iyy, izz, 5) = U_Star ( 5)  ! Energy density


          UStar_x ( ixx, iyy, izz, 6) = Ax * U_Star ( 6)  ! Magnetic field
     1             + Bx * U_Star ( 7) + Cx * U_Star ( 8)

          UStar_x ( ixx, iyy, izz, 7) = Ay * U_Star ( 6)
     1             + By * U_Star ( 7) + Cy * U_Star ( 8)

          UStar_x ( ixx, iyy, izz, 8) = Az * U_Star ( 6)
     1             + Bz * U_Star ( 7) + Cz * U_Star ( 8)

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3)
          UStar_x ( ixx, iyy, izz, 9: ) = U_Star ( 9: )  ! pressures and tracer fields
#endif

#if ( IPDE_SYSTEM == 1 && NFLUID >= 1)
          UStar_x ( ixx, iyy, izz, NVAR + 1: NVAR + NFLUID) =  ! Species density
     1                    U_Star ( NVAR + 1: NVAR + NFLUID)
#endif

! ----------
! Take care of entropy flux.

          EntFlux_x ( ixx, iyy, izz) =  ent_F_Star  ! Entropy density flux

          END IF  ! End of "If ( j .EQ. 1) THEN".

! --------------------------------------------------

          END DO  ! End of "DO ixx = 0, ntp_zones" loop.

        END DO  ! End of "DO iyy = 1, ntp_zones" loop.

      END DO  ! End of "DO izz = 1, nr_zones" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, NAREA_QUAD_PTS" loop over all facial quadrature points.

! ----------------------------------------------------------------------

! Make sure that the accumulation sum on "Efy, Efz" does not overwrite itself
! in a parallel setting.

#if ( USE_FACIAL_EFIELDS == 1)

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    iexpand,                                                      &
!$OMP&    Efy, Efz, Efy_m1, Efz_m1, i_Efy, i_Efz )                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present (                              &
!$acc&    iexpand,                                                      &
!$acc&    Efy, Efz, Efy_m1, Efz_m1, i_Efy, i_Efz )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)


      DO izz = 1, nr_zones

        DO iyy = 1, ntp_zones

          DO ixx = 0, ntp_zones

! --------------------------------------------------

          Efz ( ixx, iyy - 1, izz) = Efz ( ixx, iyy - 1, izz)
     1                             + Efz_m1 ( ixx, iyy - 1, izz)

          Efy ( ixx, iyy, izz - 1) = Efy ( ixx, iyy, izz - 1)
     1                             + Efy_m1 ( ixx, iyy, izz - 1)


          i_Efz ( ixx, iyy - 1, izz) = i_Efz ( ixx, iyy - 1, izz) + 1

          i_Efy ( ixx, iyy, izz - 1) = i_Efy ( ixx, iyy, izz - 1) + 1

! --------------------------------------------------

          END DO  ! End of "DO ixx = 0, ntp_zones" loop.

        END DO  ! End of "DO iyy = 1, ntp_zones" loop.

      END DO  ! End of "DO izz = 1, nr_zones" loop.

!$acc end parallel
!$OMP END PARALLEL DO

#endif
! End of "( USE_FACIAL_EFIELDS == 1)".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 2) Access the y-faces. From either side of the y-face we use
! "U ( ixxd, iyyd, izz, 6:8, 1:4)" and "U ( ixxu, iyyu, izz, 6:8, 1:4)" to obtain
! "U_L ( :), U_R ( :)" which are fed into the 1D Riemann solver.

! The higher order modes from "B1_zc, B2_zc, B3_zc" on either side of the
! y-face are also used for the magnetic field.

! Since this is only second order, we obtain "U_L ( :), U_R ( :)"
! at the centroid of the y-face.

! If "PRESSURE_POSITIVITY >= 1" we also use "entropygr" to obtain
! "ent_L, ent_R"
! ----------------------------------------------------------------------

!$acc kernels
      Efx_m1 ( :, :, :) = 0.0
      Efy_m1 ( :, :, :) = 0.0
      Efz_m1 ( :, :, :) = 0.0
!$acc end kernels

! Short loop over all quadrature points. Do not parallelize this outer loop.
! "j = 1" is the only central point and it is the only quadrature point needed at
! second order.
      DO j = 1, NAREA_QUAD_PTS

#if (NUMORDER == 2)
      IF ( j .GT. 1) EXIT
#endif


!$acc update device ( iexpand, j)

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxd, iyyd, ixxu, iyyu,      &
!$OMP&    x, y, z, x2, y2, z2, obl0, obl02, v1, Efld,                   &
!$OMP&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, tempaa, tempab, tempac,   &
!$OMP&    bflag, time_use, B0, Fluct_Pl, Fluct_Mi, Coef_Basis,          &

!$OMP&    S_L, S_R, ent_L, ent_R, ent_F_Star, Fluct_L, Fluct_R, flat,   &
!$OMP&    U_L, U_R, U_Star, F_Star, V_L, V_R, F_L, F_R, F_LL, F_RR,     &
!$OMP&    prs_L, prs_L0, del_prs_L, prs_R, prs_R0, del_prs_R,           &
!$OMP&    temp_01, temp_02, ob_rho, P_tot, Bf_2, Vel_2, V_dot_B,        &
!$OMP&    rhobar, vlxbar, vlybar, vlzbar, prsbar, bfxre, bfybar, bfzbar,&
!$OMP&    va_x, ob_4pi_rho, c_s, v_a, v_f, wt )

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    fourpi, ob_4pi, ob_8pi, smallnum, iexpand, j,                 &
!$OMP&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$OMP&    molewt, gasconst, rhofloor, prsfloor, tempfloor,              &
!$OMP&    imask_active_patch, istencil_patch, v_scale_length,           &
!$OMP&    yf_centroid_patch, yf_quadnorml_patch, Coef_Basis_Moments,    &
!$OMP&    yf_quadpts_patch, yf_quadwt_patch,                            &
!$OMP&    v_centroid_patch, yf_area_patch, v_volume_patch,              &
!$OMP&    U, V, entropygr, B1_zc, B2_zc, B3_zc, flattengr,              &
!$OMP&    UStar_y, EntFlux_y, S_D_gr, S_U_gr,                           &
!$OMP&    Fluct_y_Pl, Fluct_y_Mi, dU_dt,                                &
!$OMP&    magB_small_cgl, eps_bdry_cgl,                                 &
!$OMP&    Efx, Efz, Efx_m1, Efz_m1, i_Efx, i_Efz, Efx_uvec, Efz_uvec,   &
!$OMP&    gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr )           &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    fourpi, ob_4pi, ob_8pi, smallnum, iexpand, j,                 &
!$acc&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$acc&    molewt, gasconst, rhofloor, prsfloor, tempfloor,              &
!$acc&    imask_active_patch, istencil_patch, v_scale_length,           &
!$acc&    yf_centroid_patch, yf_quadnorml_patch, Coef_Basis_Moments,    &
!$acc&    yf_quadpts_patch, yf_quadwt_patch,                            &
!$acc&    v_centroid_patch, yf_area_patch, v_volume_patch,              &
!$acc&    U, V, entropygr, B1_zc, B2_zc, B3_zc, flattengr,              &
!$acc&    UStar_y, EntFlux_y, S_D_gr, S_U_gr,                           &
!$acc&    Fluct_y_Pl, Fluct_y_Mi, dU_dt,                                &
!$acc&    magB_small_cgl, eps_bdry_cgl,                                 &
!&acc&    Efx, Efz, Efx_m1, Efz_m1, i_Efx, i_Efz, Efx_uvec, Efz_uvec,   &
!$acc&    gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxd, iyyd, ixxu, iyyu,              &
!$acc&    x, y, z, x2, y2, z2, obl0, obl02, v1, Efld,                   &
!$acc&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, tempaa, tempab, tempac,   &
!$acc&    bflag, time_use, B0, Fluct_Pl, Fluct_Mi, Coef_Basis,          &

!$acc&    S_L, S_R, ent_L, ent_R, ent_F_Star, Fluct_L, Fluct_R, flat,   &
!$acc&    U_L, U_R, U_Star, F_Star, V_L, V_R, F_L, F_R, F_LL, F_RR,     &
!$acc&    prs_L, prs_L0, del_prs_L, prs_R, prs_R0, del_prs_R,           &
!$acc&    temp_01, temp_02, ob_rho, P_tot, Bf_2, Vel_2, V_dot_B,        &
!$acc&    rhobar, vlxbar, vlybar, vlzbar, prsbar, bfxre, bfybar, bfzbar,&
!$acc&    va_x, ob_4pi_rho, c_s, v_a, v_f, wt )

      DO izz = 1 - iexpand, nr_zones + iexpand

        DO iyy = 0 - iexpand, ntp_zones + iexpand

          DO ixx = 1 - iexpand, ntp_zones + iexpand

! ------------------------------------------------------------

          IF ( ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) .AND.
     1         ( imask_active_patch ( ixx, iyy + 1, izz) .EQ. 1) ) THEN

            ixxd = ixx
            iyyd = iyy

            ixxu = ixx
            iyyu = iyy + 1

          ELSE IF ( ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) .AND.
     1         ( imask_active_patch ( ixx, iyy + 1, izz) .EQ. 1) ) THEN

! Physical zone to the Up of y-face
            ixxu = ixx
            iyyu = iyy + 1

! Physical zone to the Dn of y-face
            ixxd = istencil_patch ( ixxu, iyyu, 1, 4)
            iyyd = istencil_patch ( ixxu, iyyu, 2, 4)

          ELSE IF ( ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) .AND.
     1         ( imask_active_patch ( ixx, iyy + 1, izz) .EQ. 0) ) THEN

! Physical zone to the Dn of y-face
            ixxd = ixx
            iyyd = iyy

! Physical zone to the Up of y-face
            ixxu = istencil_patch ( ixxd, iyyd, 1, 3)
            iyyu = istencil_patch ( ixxd, iyyd, 2, 3)

          ELSE

            V_L ( :) = 0.0
            V_R ( :) = 0.0

            ent_L = prsfloor / rhofloor**(gamma - 1.0)
            ent_R = ent_L

            CYCLE  ! Go to the next iteration of the "DO ixx =" loop.

          END IF

! --------------------------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixxd, iyyd, izz), smallnum)  ! Scaling for Dn zone.

          obl02 = 1.0
     1  / AMAX1 ( v_scale_length ( ixxu, iyyu, izz), smallnum)  ! Scaling for Up zone.


          Coef_Basis ( :)  = Coef_Basis_Moments ( ixx, iyy, izz, :)

! ------------------------------
! Obtain "Ax, Ay, Az" the normal to this face. This is the local x-direction.
! The entire y-face has the same normal. "Bx, By, Bz" is in the radial direction.
! This is also the local y-direction. The local z-direction is given by the
! vector cross product "C = A X B".

! In this loop, we rotate vectors into the local xyz-directions. Once we
! solve the Riemann problem, we will rotate flux components back into the global
! xyz-directions. For this reason, we save "Ax, Ay, Az" etc. as persistent 1d arrays.


          Ax = yf_quadnorml_patch ( ixx, iyy, izz, 1)  ! Unit vector in Local x-direction.
          Ay = yf_quadnorml_patch ( ixx, iyy, izz, 2)
          Az = yf_quadnorml_patch ( ixx, iyy, izz, 3)


          Bx = yf_centroid_patch ( ixx, iyy, izz, 1)
          By = yf_centroid_patch ( ixx, iyy, izz, 2)
          Bz = yf_centroid_patch ( ixx, iyy, izz, 3)

          tempaa = AMAX1 ( SQRT ( Bx**2 + By**2
     1                          + Bz**2), smallnum)
          tempaa = 1.0 / tempaa

          Bx = Bx * tempaa  ! Unit vector in Local y-direction.
          By = By * tempaa
          Bz = Bz * tempaa


          Cx = Ay * Bz - Az * By  ! Unit vector in Local z-direction.
          Cy = Az * Bx - Ax * Bz
          Cz = Ax * By - Ay * Bx


#if( 0 == 1)
! Validate C is unit vector. In turn validates A, B are perfectly normal.
          IF ( ABS ( SQRT ( Cx**2 + Cy**2 + Cz**2)
     1             - 1.0 ) .GE. smallnum ) THEN
            WRITE ( 6, *)"Stopping at ",__LINE__, "in ",__FILE__
            CALL FLUSH ( 6)
            STOP
          END IF
#endif
! ------------------------------
! Build scaled relative coordinates from Dn zone centroid and Up zone centroid to
! facial centroid. Only one Riemann solve at second order. (BTW, at higher orders
! the facial centroid would be replaced by looping over facial quadrature points.)

          x = ( yf_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixxd, iyyd, izz, 1) ) * obl0

          y = ( yf_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixxd, iyyd, izz, 2) ) * obl0

          z = ( yf_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixxd, iyyd, izz, 3) ) * obl0

! ----------

          x2 = ( yf_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixxu, iyyu, izz, 1) ) * obl02

          y2 = ( yf_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixxu, iyyu, izz, 2) ) * obl02

          z2 = ( yf_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixxu, iyyu, izz, 3) ) * obl02

! ----------
! Obtain the weight for this quadrature point. Weight is 1 at second order.

          wt = 1.0
#if (NUMORDER >= 3)
          wt = yf_quadwt_patch ( ixx, iyy, izz, j)
#endif

! ------------------------------
! Now use the above relative coordinates as well as
! "U ( ixxd, iyyd, izz, :, 1:4)" and "U ( ixxu, iyyu, izz, :, 1:4)"
! to build "U_L ( :), U_R ( :)" with second order accuracy.
! Note the use of "x, y, z" from down and "x2, y2, z2" from up.

          V_L ( :) = V ( ixxd, iyyd, izz, :, 1)
     1             + V ( ixxd, iyyd, izz, :, 2) * x
     1             + V ( ixxd, iyyd, izz, :, 3) * y
     1             + V ( ixxd, iyyd, izz, :, 4) * z
#if (NUMORDER >= 3)
     1             + V ( ixxd, iyyd, izz, :, 5)
     1             * (x**2 - Coef_Basis_Moments ( ixxd, iyyd, izz, 5))
     1             + V ( ixxd, iyyd, izz, :, 6)
     1             * (y**2 - Coef_Basis_Moments ( ixxd, iyyd, izz, 6))
     1             + V ( ixxd, iyyd, izz, :, 7)
     1             * (z**2 - Coef_Basis_Moments ( ixxd, iyyd, izz, 7))
     1             + V ( ixxd, iyyd, izz, :, 8)
     1             * (x * y - Coef_Basis_Moments ( ixxd, iyyd, izz, 8))
     1             + V ( ixxd, iyyd, izz, :, 9)
     1             * (y * z - Coef_Basis_Moments ( ixxd, iyyd, izz, 9))
     1             + V ( ixxd, iyyd, izz, :, 10)
     1             * (x * z - Coef_Basis_Moments ( ixxd, iyyd, izz, 10))

#elif (NUMORDER >= 4)
     1             + V ( ixxd, iyyd, izz, :, 11)
     1             * (x**3 - Coef_Basis_Moments ( ixxd, iyyd, izz, 11))
     1             + V ( ixxd, iyyd, izz, :, 12)
     1             * (y**3 - Coef_Basis_Moments ( ixxd, iyyd, izz, 12))
     1             + V ( ixxd, iyyd, izz, :, 13)
     1             * (z**3 - Coef_Basis_Moments ( ixxd, iyyd, izz, 13))
     1             + V ( ixxd, iyyd, izz, :, 14)
     1             * (x**2 * y
     1             - Coef_Basis_Moments ( ixxd, iyyd, izz, 14))
     1             + V ( ixxd, iyyd, izz, :, 15)
     1             * (x**2 * z
     1             - Coef_Basis_Moments ( ixxd, iyyd, izz, 15))
     1             + V ( ixxd, iyyd, izz, :, 16)
     1             * (x * y**2
     1             - Coef_Basis_Moments ( ixxd, iyyd, izz, 16))
     1             + V ( ixxd, iyyd, izz, :, 17)
     1             * (y**2 * z
     1             - Coef_Basis_Moments ( ixxd, iyyd, izz, 17))
     1             + V ( ixxd, iyyd, izz, :, 18)
     1             * (x * z**2
     1             - Coef_Basis_Moments ( ixxd, iyyd, izz, 18))
     1             + V ( ixxd, iyyd, izz, :, 19)
     1             * (y * z**2
     1             - Coef_Basis_Moments ( ixxd, iyyd, izz, 19))
     1             + V ( ixxd, iyyd, izz, :, 20)
     1             * (x * y * z
     1             - Coef_Basis_Moments ( ixxd, iyyd, izz, 20))
#endif

          V_R ( :) = V ( ixxu, iyyu, izz, :, 1)
     1             + V ( ixxu, iyyu, izz, :, 2) * x2
     1             + V ( ixxu, iyyu, izz, :, 3) * y2
     1             + V ( ixxu, iyyu, izz, :, 4) * z2
#if (NUMORDER >= 3)
     1             + V ( ixxu, iyyu, izz, :, 5)
     1             * (x2**2 - Coef_Basis_Moments ( ixxu, iyyu, izz, 5))
     1             + V ( ixxu, iyyu, izz, :, 6)
     1             * (y2**2 - Coef_Basis_Moments ( ixxu, iyyu, izz, 6))
     1             + V ( ixxu, iyyu, izz, :, 7)
     1             * (z2**2 - Coef_Basis_Moments ( ixxu, iyyu, izz, 7))
     1             + V ( ixxu, iyyu, izz, :, 8)
     1             * (x2 * y2
     1             - Coef_Basis_Moments ( ixxu, iyyu, izz, 8))
     1             + V ( ixxu, iyyu, izz, :, 9)
     1             * (y2 * z2
     1             - Coef_Basis_Moments ( ixxu, iyyu, izz, 9))
     1             + V ( ixxu, iyyu, izz, :, 10)
     1             * (x2 * z2
     1             - Coef_Basis_Moments ( ixxu, iyyu, izz, 10))
#elif (NUMORDER >= 4)
     1             + V ( ixxu, iyyu, izz, :, 11)
     1             * (x2**3 - Coef_Basis_Moments ( ixxu, iyyu, izz, 11))
     1             + V ( ixxu, iyyu, izz, :, 12)
     1             * (y2**3 - Coef_Basis_Moments ( ixxu, iyyu, izz, 12))
     1             + V ( ixxu, iyyu, izz, :, 13)
     1             * (z2**3 - Coef_Basis_Moments ( ixxu, iyyu, izz, 13))
     1             + V ( ixxu, iyyu, izz, :, 14)
     1             * (x2**2 * y2
     1             - Coef_Basis_Moments ( ixxu, iyyu, izz, 14))
     1             + V ( ixxu, iyyu, izz, :, 15)
     1             * (x2**2 * z2
     1             - Coef_Basis_Moments ( ixxu, iyyu, izz, 15))
     1             + V ( ixxu, iyyu, izz, :, 16)
     1             * (x2 * y2**2
     1             - Coef_Basis_Moments ( ixxu, iyyu, izz, 16))
     1             + V ( ixxu, iyyu, izz, :, 17)
     1             * (y2**2 * z2
     1             - Coef_Basis_Moments ( ixxu, iyyu, izz, 17))
     1             + V ( ixxu, iyyu, izz, :, 18)
     1             * (x2 * z2**2
     1             - Coef_Basis_Moments ( ixxu, iyyu, izz, 18))
     1             + V ( ixxu, iyyu, izz, :, 19)
     1             * (y2 * z2**2
     1             - Coef_Basis_Moments ( ixxu, iyyu, izz, 19))
     1             + V ( ixxu, iyyu, izz, :, 20)
     1             * (x2 * y2 * z2
     1             - Coef_Basis_Moments ( ixxu, iyyu, izz, 20))
#endif

! ------------------------------
! Now use the above relative coordinates as well as
! "entropygr ( ixxd, iyyd, izz, 1:4)" and "entropygr ( ixxu, iyyu, izz, 1:4)"
! to build "ent_L, ent_R" with second order accuracy.

#if ( PRESSURE_POSITIVITY >= 1)

          ent_L = entropygr ( ixxd, iyyd, izz, 1)
     1          + entropygr ( ixxd, iyyd, izz, 2) * x
     1          + entropygr ( ixxd, iyyd, izz, 3) * y
     1          + entropygr ( ixxd, iyyd, izz, 4) * z

          ent_R = entropygr ( ixxu, iyyu, izz, 1)
     1          + entropygr ( ixxu, iyyu, izz, 2) * x2
     1          + entropygr ( ixxu, iyyu, izz, 3) * y2
     1          + entropygr ( ixxu, iyyu, izz, 4) * z2

#else
          ent_L = prsfloor / rhofloor**(gamma - 1.0)
          ent_R = ent_L
#endif

! ------------------------------
! Notice that the 1D Riemann solver always assumes that the waves expand
! in the local x-direction. Therefore, we have to rotate the vectorial entities
! like velocities and magnetic fields into the local coordinate system.

          tempaa = Ax * V_L ( 2) + Ay * V_L ( 3) + Az * V_L ( 4)
          tempab = Bx * V_L ( 2) + By * V_L ( 3) + Bz * V_L ( 4)
          tempac = Cx * V_L ( 2) + Cy * V_L ( 3) + Cz * V_L ( 4)

          V_L ( 2) = tempaa
          V_L ( 3) = tempab
          V_L ( 4) = tempac

          tempaa = Ax * V_L ( 6) + Ay * V_L ( 7) + Az * V_L ( 8)
          tempab = Bx * V_L ( 6) + By * V_L ( 7) + Bz * V_L ( 8)
          tempac = Cx * V_L ( 6) + Cy * V_L ( 7) + Cz * V_L ( 8)

          V_L ( 6) = tempaa
          V_L ( 7) = tempab
          V_L ( 8) = tempac

! ----------

          tempaa = Ax * V_R ( 2) + Ay * V_R ( 3) + Az * V_R ( 4)
          tempab = Bx * V_R ( 2) + By * V_R ( 3) + Bz * V_R ( 4)
          tempac = Cx * V_R ( 2) + Cy * V_R ( 3) + Cz * V_R ( 4)

          V_R ( 2) = tempaa
          V_R ( 3) = tempab
          V_R ( 4) = tempac

          tempaa = Ax * V_R ( 6) + Ay * V_R ( 7) + Az * V_R ( 8)
          tempab = Bx * V_R ( 6) + By * V_R ( 7) + Bz * V_R ( 8)
          tempac = Cx * V_R ( 6) + Cy * V_R ( 7) + Cz * V_R ( 8)

          V_R ( 6) = tempaa
          V_R ( 7) = tempab
          V_R ( 8) = tempac

! --------------------------------------------------
! Obtain the background magnetic field B0 and rotate it.

          B0 ( :) = 0.0

#if ( USE_B0B1_FORMULATION == 1)

          x = yf_centroid_patch ( ixx, iyy, izz, 1)
          y = yf_centroid_patch ( ixx, iyy, izz, 2)
          z = yf_centroid_patch ( ixx, iyy, izz, 3)

          bflag = 1
          time_use = 0

! Calling acc-seq function to get the B0, using temp. variables x2, y2, z2.
          CALL GET_MAGNETIC_FIELD ( bflag, time_use,  x, y, z,
     1                              x2, y2, z2 )

          B0 ( 1) = Ax * x2 + Ay * y2 + Az * z2
          B0 ( 2) = Bx * x2 + By * y2 + Bz * z2
          B0 ( 3) = Cx * x2 + Cy * y2 + Cz * z2

#endif

! ------------------------------------------------------------

          flat = AMAX1 ( flattengr ( ixxd, iyyd, izz),
     1                   flattengr ( ixxu, iyyu, izz) )

! ------------------------------------------------------------
! Now call the 1D Riemann solver. This RS operates in the local x-direction.
! It takes in "V_L, V_R" and "ent_L, ent_R" and returns speeds "S_L, S_R".
! Also returns resolved state "U_Star", resolved flux "F_Star", entropy flux "ent_F_Star".
! Also returns fluctuation "Fluct_L, Fluct_R" and "U_L, U_R, F_L, F_R".

          CALL OneD_RS_MHD_Ptwise

     1         ( gamma, gamma_euler, isospeed, barotropic_stiffness,
     1           molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1           magB_small_cgl, eps_bdry_cgl,

     1           fourpi, ob_4pi, ob_8pi,
     1           smallnum, B0, flat,

     1           V_L, V_R, ent_L, ent_R,
     1           S_L, S_R, U_Star, F_Star, ent_F_Star,
     1           Fluct_L, Fluct_R,
     1           U_L, U_R, F_L, F_R)

! ------------------------------------------------------------
! Hand back the resolved states, fluxes and extremal wave speeds from the
! 1D Riemann solver.

          IF ( ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) .AND.
     1         ( imask_active_patch ( ixx, iyy + 1, izz) .EQ. 0) ) THEN
! If we get here, we are at an unphysical boundary, so set zero values and move on.

            UStar_y ( ixx, iyy, izz, :) = 0.0
            EntFlux_y ( ixx, iyy, izz) = 0.0

            CYCLE  ! Go to the next iteration of the "DO ixx =" loop.

          END IF

! ------------------------------
! If we get here, we are at a physical boundary, so assign the good values
! from the 1D dime-a-dozen Riemann solver.

! ----------
! If "USE_FACIAL_EFIELDS == 1" then "Ez = -F7" and "Ey = F8" in the local
! coordinates. Therefore, before we rotate back, we assign "Efx, Efz".

#if ( USE_FACIAL_EFIELDS == 1)
          IF ( j .EQ. 1) THEN

          v1 ( 1)   = 0.0
          v1 ( 2)   =   F_Star ( 8)
          v1 ( 3)   = - F_Star ( 7)

! Accumulate the upper edge.
          Efz ( ixx, iyy, izz) = Efz ( ixx, iyy, izz) + v1 ( 2)
          Efx ( ixx, iyy, izz) = Efx ( ixx, iyy, izz) + v1 ( 3)

! Count the number of accumulations.
          i_Efx ( ixx, iyy, izz) = i_Efx ( ixx, iyy, izz) + 1
          i_Efz ( ixx, iyy, izz) = i_Efz ( ixx, iyy, izz) + 1

! Assign directly for the lower edge.
          Efz_m1 ( ixx - 1, iyy, izz) = v1 ( 2)
          Efx_m1 ( ixx, iyy, izz - 1) = v1 ( 3)

          END IF
#endif

! ----------
! Save the extremal speeds. They will be needed in the 2D RS for electric field evaluation.
! This is only done for the face-centered quadrature point, because it has the best speeds.

          IF ( j .EQ. 1) THEN

          S_D_gr ( ixx, iyy, izz) = S_L

          S_U_gr ( ixx, iyy, izz) = S_R

          END IF

! ----------
! For CGL we have to use fluctuation form. For that, we return "Fluct_L, Fluct_R".
! These will also have to undergo the vectorial components in the same way.

! ----------
! Fluctuations in + direction; rotated into global frame.

          Fluct_Pl ( 1) = Fluct_R ( 1)                        ! Density

          Fluct_Pl ( 2) = Ax * Fluct_R ( 2)                   ! x-Momentum
     1               + Bx * Fluct_R ( 3) + Cx * Fluct_R ( 4)
          Fluct_Pl ( 3) = Ay * Fluct_R ( 2)                   ! y-Momentum
     1               + By * Fluct_R ( 3) + Cy * Fluct_R ( 4)
          Fluct_Pl ( 4) = Az * Fluct_R ( 2)                   ! z-Momentum
     1               + Bz * Fluct_R ( 3) + Cz * Fluct_R ( 4)

          Fluct_Pl ( 5) = Fluct_R ( 5)                        ! Energy

          Fluct_Pl ( 6) = Ax * Fluct_R ( 6)                   ! x-Magnetic field
     1               + Bx * Fluct_R ( 7) + Cx * Fluct_R ( 8)
          Fluct_Pl ( 7) = Ay * Fluct_R ( 6)                   ! y-Magnetic field
     1               + By * Fluct_R ( 7) + Cy * Fluct_R ( 8)
          Fluct_Pl ( 8) = Az * Fluct_R ( 6)                   ! z-Magnetic field
     1               + Bz * Fluct_R ( 7) + Cz * Fluct_R ( 8)

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3 )
          Fluct_Pl ( 9: ) = Fluct_R ( 9: )  ! pressures and tracer fields
#endif

#if ( IPDE_SYSTEM == 1 && NFLUID >= 1)
          Fluct_Pl ( NVAR + 1: NVAR + NFLUID) =  ! Species
     1     Fluct_R ( NVAR + 1: NVAR + NFLUID)
#endif

! -----

          Fluct_y_Pl ( ixx, iyy, izz, :) =Fluct_y_Pl ( ixx, iyy, izz, :)
     1    + wt * Fluct_Pl ( :)

! ----------
! Fluctuations in - direction; rotated into global frame.

          Fluct_Mi ( 1) = Fluct_L ( 1)                        ! Density

          Fluct_Mi ( 2) = Ax * Fluct_L ( 2)                   ! x-Momentum
     1               + Bx * Fluct_L ( 3) + Cx * Fluct_L ( 4)
          Fluct_Mi ( 3) = Ay * Fluct_L ( 2)                   ! y-Momentum
     1               + By * Fluct_L ( 3) + Cy * Fluct_L ( 4)
          Fluct_Mi ( 4) = Az * Fluct_L ( 2)                   ! z-Momentum
     1               + Bz * Fluct_L ( 3) + Cz * Fluct_L ( 4)

          Fluct_Mi ( 5) = Fluct_L ( 5)                        ! Energy

          Fluct_Mi ( 6) = Ax * Fluct_L ( 6)                   ! x-Magnetic field
     1               + Bx * Fluct_L ( 7) + Cx * Fluct_L ( 8)
          Fluct_Mi ( 7) = Ay * Fluct_L ( 6)                   ! y-Magnetic field
     1               + By * Fluct_L ( 7) + Cy * Fluct_L ( 8)
          Fluct_Mi ( 8) = Az * Fluct_L ( 6)                   ! z-Magnetic field
     1               + Bz * Fluct_L ( 7) + Cz * Fluct_L ( 8)

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3 )
          Fluct_Mi ( 9: ) = Fluct_L ( 9: )  ! pressures and tracer fields
#endif

#if ( IPDE_SYSTEM == 1 && NFLUID >= 1)
          Fluct_Mi ( NVAR + 1: NVAR + NFLUID) =  ! Species
     1     Fluct_L ( NVAR + 1: NVAR + NFLUID)
#endif

! -----

          Fluct_y_Mi ( ixx, iyy, izz, :) =Fluct_y_Mi ( ixx, iyy, izz, :)
     1    + wt * Fluct_Mi ( :)

! ----------
! Rotate F_L to F_LL; rotate F_R to F_RR and add in the internal terms to the update.
! This is only done for the momentum fluxes and the magnetic field fluxes.

          F_LL = F_L

          F_LL ( 2) = Ax * F_L ( 2) + Bx * F_L ( 3) + Cx * F_L ( 4)
          F_LL ( 3) = Ay * F_L ( 2) + By * F_L ( 3) + Cy * F_L ( 4)
          F_LL ( 4) = Az * F_L ( 2) + Bz * F_L ( 3) + Cz * F_L ( 4)

          F_LL ( 6) = Ax * F_L ( 6) + Bx * F_L ( 7) + Cx * F_L ( 8)
          F_LL ( 7) = Ay * F_L ( 6) + By * F_L ( 7) + Cy * F_L ( 8)
          F_LL ( 8) = Az * F_L ( 6) + Bz * F_L ( 7) + Cz * F_L ( 8)

! -----

          F_RR = F_R

          F_RR ( 2) = Ax * F_R ( 2) + Bx * F_R ( 3) + Cx * F_R ( 4)
          F_RR ( 3) = Ay * F_R ( 2) + By * F_R ( 3) + Cy * F_R ( 4)
          F_RR ( 4) = Az * F_R ( 2) + Bz * F_R ( 3) + Cz * F_R ( 4)

          F_RR ( 6) = Ax * F_R ( 6) + Bx * F_R ( 7) + Cx * F_R ( 8)
          F_RR ( 7) = Ay * F_R ( 6) + By * F_R ( 7) + Cy * F_R ( 8)
          F_RR ( 8) = Az * F_R ( 6) + Bz * F_R ( 7) + Cz * F_R ( 8)

! -----
! Now add in the contributions from the internal fluxes to the update terms in
! the zone "( ixxd, iyyd, izz, :)" and "( ixxu, iyyu, izz, :)".

          dU_dt ( ixxd, iyyd, izz, :) = dU_dt ( ixxd, iyyd, izz, :)
     1  - wt * ( Fluct_Mi ( :) + F_LL ( :) )
     1             * yf_area_patch ( ixx, iyy, izz)
     1             / v_volume_patch ( ixxd, iyyd, izz)

! -----

          dU_dt ( ixxu, iyyu, izz, :) = dU_dt ( ixxu, iyyu, izz, :)
     1  - wt * ( Fluct_Pl ( :) - F_RR ( :) )
     1             * yf_area_patch ( ixx, iyy, izz)
     1             / v_volume_patch ( ixxu, iyyu, izz)

! ----------
! Likewise, rotate the state components back into the global xyz-directions.

          IF ( j .EQ. 1) THEN

          UStar_y ( ixx, iyy, izz, 1) = U_Star ( 1)  ! Density


          UStar_y ( ixx, iyy, izz, 2) = Ax * U_Star ( 2)  ! Momentum density
     1             + Bx * U_Star ( 3) + Cx * U_Star ( 4)

          UStar_y ( ixx, iyy, izz, 3) = Ay * U_Star ( 2)
     1             + By * U_Star ( 3) + Cy * U_Star ( 4)

          UStar_y ( ixx, iyy, izz, 4) = Az * U_Star ( 2)
     1             + Bz * U_Star ( 3) + Cz * U_Star ( 4)


          UStar_y ( ixx, iyy, izz, 5) = U_Star ( 5)  ! Energy density


          UStar_y ( ixx, iyy, izz, 6) = Ax * U_Star ( 6)  ! Magnetic field
     1             + Bx * U_Star ( 7) + Cx * U_Star ( 8)

          UStar_y ( ixx, iyy, izz, 7) = Ay * U_Star ( 6)
     1             + By * U_Star ( 7) + Cy * U_Star ( 8)

          UStar_y ( ixx, iyy, izz, 8) = Az * U_Star ( 6)
     1             + Bz * U_Star ( 7) + Cz * U_Star ( 8)

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3 )
          UStar_y ( ixx, iyy, izz, 9: ) = U_Star ( 9: )  ! pressures and tracer fields
#endif

#if ( IPDE_SYSTEM == 1 && NFLUID >= 1)
          UStar_y ( ixx, iyy, izz, NVAR + 1: NVAR + NFLUID) =  ! Species density
     1                    U_Star ( NVAR + 1: NVAR + NFLUID)
#endif

! ----------
! Take care of entropy flux.

          EntFlux_y ( ixx, iyy, izz) =  ent_F_Star  ! Entropy density flux

          END IF  ! End of "If ( j .EQ. 1) THEN".

! ------------------------------------------------------------

          END DO  ! End of "DO ixx = 1, ntp_zones" loop.

        END DO  ! End of "DO iyy = 0, ntp_zones" loop.

      END DO  ! End of "DO izz = 1, nr_zones" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, NAREA_QUAD_PTS" loop over all facial quadrature points.

! ----------------------------------------------------------------------
! Make sure that the accumulation sum on "Efx, Efz" does not overwrite itself
! in a parallel setting.

#if ( USE_FACIAL_EFIELDS == 1)

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    iexpand,                                                      &
!$OMP&    Efx, Efz, Efx_m1, Efz_m1, i_Efx, i_Efz )                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present (                              &
!$acc&    iexpand, Efx, Efz, Efx_m1, Efz_m1, i_Efx, i_Efz )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 1, nr_zones

        DO iyy = 0, ntp_zones

          DO ixx = 1, ntp_zones

! --------------------------------------------------

          Efz ( ixx - 1, iyy, izz) = Efz ( ixx - 1, iyy, izz)
     1                             + Efz_m1 ( ixx - 1, iyy, izz)

          Efx ( ixx, iyy, izz - 1) = Efx ( ixx, iyy, izz - 1)
     1                             + Efx_m1 ( ixx, iyy, izz - 1)


          i_Efz ( ixx - 1, iyy, izz) = i_Efz ( ixx - 1, iyy, izz) + 1

          i_Efx ( ixx, iyy, izz - 1) = i_Efx ( ixx, iyy, izz - 1) + 1

! --------------------------------------------------

          END DO  ! End of "DO ixx = 1, ntp_zones" loop.

        END DO  ! End of "DO iyy = 0, ntp_zones" loop.

      END DO  ! End of "DO izz = 1, nr_zones" loop.

!$acc end parallel
!$OMP END PARALLEL DO

#endif
! End of "( USE_FACIAL_EFIELDS == 1)".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 3) Access the z-faces. From either side of the z-face we use
! "U ( ixx, iyy, izz, :, 1:4)" and "U ( ixx, iyy, izzp1, :, 1:4)" to obtain
! "U_L ( :), U_R ( :)" which are fed into the 1D Riemann solver.

! The higher order modes from "B1_zc, B2_zc, B3_zc" on either side of the
! z-face are also used for the magnetic field.

! Since this is only second order, we obtain "U_L ( :), U_R ( :)"
! at the centroid of the z-face.

! If "PRESSURE_POSITIVITY >= 1" we also use "entropygr" to obtain
! "ent_L, ent_R"
! ----------------------------------------------------------------------

!$acc kernels
      Efx_m1 ( :, :, :) = 0.0
      Efy_m1 ( :, :, :) = 0.0
      Efz_m1 ( :, :, :) = 0.0
!$acc end kernels

! Short loop over all quadrature points. Do not parallelize this outer loop.
! "j = 1" is the only central point and it is the only quadrature point needed at
! second order.
      DO j = 1, NAREA_QUAD_PTS

#if (NUMORDER == 2)
      IF ( j .GT. 1) EXIT
#endif


!$acc update device ( iexpand, j)

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izzp1,                       &
!$OMP&    x, y, z, x2, y2, z2, obl0, obl02, v1, Efld,                   &
!$OMP&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz,                           &
!$OMP&    tempaa, tempab, tempac, Fluct_Pl, Fluct_Mi, Coef_Basis,       &
!$OMP&    bflag, time_use, B0,                                          &

!$OMP&    S_L, S_R, ent_L, ent_R, ent_F_Star, Fluct_L, Fluct_R, flat,   &
!$OMP&    U_L, U_R, U_Star, F_Star, V_L, V_R, F_L, F_R, F_LL, F_RR,     &
!$OMP&    prs_L, prs_L0, del_prs_L, prs_R, prs_R0, del_prs_R,           &
!$OMP&    temp_01, temp_02, ob_rho, P_tot, Bf_2, Vel_2, V_dot_B,        &
!$OMP&    rhobar, vlxbar, vlybar, vlzbar, prsbar, bfxre, bfybar, bfzbar,&
!$OMP&    prs_avg, c_sL, v_anL, v_aL, d_L2, v_fnL, beta, theta,         &
!$OMP&    c_sR, v_anR, v_aR, d_R2, v_fnR, v_atL, v_snL, alpha,          &
!$OMP&    iOriginal_RS,                                                 &
!$OMP&    va_x, ob_4pi_rho, c_s, v_a, v_f, V_Star, c_ss, wt )

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    fourpi, ob_4pi, ob_8pi, smallnum, iexpand, j,                 &
!$OMP&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$OMP&    molewt, gasconst, rhofloor, prsfloor, tempfloor,              &
!$OMP&    imask_active_patch, istencil_patch, v_scale_length,           &
!$OMP&    zf_centroid_patch, zf_quadnorml_patch, vert_coords_patch,     &
!$OMP&    zf_quadpts_patch, zf_quadwt_patch,                            &
!$OMP&    v_centroid_patch, Coef_Basis_Moments, zf_area_patch,          &
!$OMP&    U, V, entropygr, B1_zc, B2_zc, B3_zc, v_volume_patch,         &
!$OMP&    UStar_z, EntFlux_z, S_B_gr, S_T_gr, flattengr,                &
!$OMP&    Fluct_z_Pl, Fluct_z_Mi, dU_dt,                                &
!$OMP&    magB_small_cgl, eps_bdry_cgl,                                 &
!$OMP&    Efx, Efy, Efx_m1, Efy_m1, i_Efx, i_Efy, Efx_uvec, Efy_uvec,   &
!$OMP&    ibdy_outer,                                                   &
!$OMP&    gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr )           &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    fourpi, ob_4pi, ob_8pi, smallnum, iexpand, j,                 &
!$acc&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$acc&    molewt, gasconst, rhofloor, prsfloor, tempfloor,              &
!$acc&    imask_active_patch, istencil_patch, v_scale_length,           &
!$acc&    zf_centroid_patch, zf_quadnorml_patch, vert_coords_patch,     &
!$acc&    zf_quadpts_patch, zf_quadwt_patch,                            &
!$acc&    v_centroid_patch, Coef_Basis_Moments,                         &
!$acc&    U, entropygr, B1_zc, B2_zc, B3_zc, flattengr,                 &
!$acc&    UStar_z, EntFlux_z, S_B_gr, S_T_gr,                           &
!$acc&    Fluct_z_Pl, Fluct_z_Mi, dU_dt,                                &
!$acc&    Efx, Efy, Efx_m1, Efy_m1, i_Efx, i_Efy, Efx_uvec, Efy_uvec,   &
!$acc&    magB_small_cgl, eps_bdry_cgl,                                 &
!$acc&    gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr )           &

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, izzp1,                               &
!$acc&    x, y, z, x2, y2, z2, obl0, obl02, v1, Efld,                   &
!$acc&    Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz,                           &
!$acc&    tempaa, tempab, tempac, Fluct_Pl, Fluct_Mi, Coef_Basis,       &
!$acc&    bflag, time_use, B0,                                          &

!$acc&    S_L, S_R, ent_L, ent_R, ent_F_Star, Fluct_L, Fluct_R, flat,   &
!$acc&    U_L, U_R, U_Star, F_Star, V_L, V_R, F_L, F_R, F_LL, F_RR,     &
!$acc&    prs_L, prs_L0, del_prs_L, prs_R, prs_R0, del_prs_R,           &
!$acc&    temp_01, temp_02, ob_rho, P_tot, Bf_2, Vel_2, V_dot_B,        &
!$acc&    rhobar, vlxbar, vlybar, vlzbar, prsbar, bfxre, bfybar, bfzbar,&
!$acc&    va_x, ob_4pi_rho, c_s, v_a, v_f, V_Star, c_ss, wt )

      DO izz = 0 - iexpand, nr_zones + iexpand

        DO iyy = 1 - iexpand, ntp_zones + iexpand

          DO ixx = 1 - iexpand, ntp_zones + iexpand

! ------------------------------------------------------------

          izzp1 = izz + 1

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) THEN

            V_L ( :) = 0.0
            V_R ( :) = 0.0

            ent_L = prsfloor / rhofloor**(gamma - 1.0)
            ent_R = ent_L

            CYCLE  ! Go to the next iteration of the "DO ixx =" loop.

          END IF

! --------------------------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for "izz" zone.

          obl02 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izzp1), smallnum)  ! Scaling for "izzp1" zone.


          Coef_Basis ( :)  = Coef_Basis_Moments ( ixx, iyy, izz, :)

! ------------------------------
! Obtain "Ax, Ay, Az" the normal to this radial face. This is the local x-direction.
! Different quadrature points in this radial face will have different normal
! directions, but we only pick the centroid at second order.
! "Bx, By, Bz" is in the global x direction, which is the local y-direction.
! The local z-direction is given by the vector cross product "C = A X B".

! Seen from the outside of the sphere, the zones in the surface look like
! parallelograms. Therefore "ty, tz" gives us the components of a unit vector in
! the local tangent plane that is parallel to the global y-direction of the mesh.


          Ax = zf_quadnorml_patch ( ixx, iyy, izz, 1, j)  ! Each node has a different normal.
          Ay = zf_quadnorml_patch ( ixx, iyy, izz, 2, j)
          Az = zf_quadnorml_patch ( ixx, iyy, izz, 3, j)

          tempaa = AMAX1 ( SQRT ( Ax**2 + Ay**2
     1                          + Az**2), smallnum)
          tempaa = 1.0 / tempaa

          Ax = Ax * tempaa  ! Unit vector in Local x-direction.
          Ay = Ay * tempaa
          Az = Az * tempaa


          Bx = vert_coords_patch ( ixx, iyy, izz, 1)
     1       + vert_coords_patch ( ixx, iyy - 1, izz, 1)
     1       - vert_coords_patch ( ixx - 1, iyy, izz, 1)
     1       - vert_coords_patch ( ixx - 1, iyy - 1, izz, 1)

          By = vert_coords_patch ( ixx, iyy, izz, 2)
     1       + vert_coords_patch ( ixx, iyy - 1, izz, 2)
     1       - vert_coords_patch ( ixx - 1, iyy, izz, 2)
     1       - vert_coords_patch ( ixx - 1, iyy - 1, izz, 2)

          Bz = vert_coords_patch ( ixx, iyy, izz, 3)
     1       + vert_coords_patch ( ixx, iyy - 1, izz, 3)
     1       - vert_coords_patch ( ixx - 1, iyy, izz, 3)
     1       - vert_coords_patch ( ixx - 1, iyy - 1, izz, 3)

          tempaa = Bx * Ax + By * Ay + Bz * Az

          Bx = Bx - tempaa * Ax  ! This ensures that B vector is perpendicular to A vector.
          By = By - tempaa * Ay
          Bz = Bz - tempaa * Az

          tempaa = AMAX1 ( SQRT ( Bx**2 + By**2 + Bz**2), smallnum)
          tempaa = 1.0 / tempaa

          Bx = Bx * tempaa  ! Unit vector in Local y-direction.
          By = By * tempaa
          Bz = Bz * tempaa

! ----------

          Cx = Ay * Bz - Az * By  ! Unit vector in Local z-direction.
          Cy = Az * Bx - Ax * Bz
          Cz = Ax * By - Ay * Bx

          tempaa = AMAX1 ( SQRT ( Cx**2 + Cy**2 + Cz**2), smallnum)
          tempaa = 1.0 / tempaa

          Cx = Cx * tempaa  ! Unit vector in Local z-direction.
          Cy = Cy * tempaa
          Cz = Cz * tempaa

          Bx = Cy * Az - Cz * Ay  ! unit-vec in Local y-dir.
          By = Cz * Ax - Cx * Az
          Bz = Cx * Ay - Cy * Ax

! ----------
#if( 0 == 1)
! Validate C is unit vector. In turn validates A, B are perfectly normal.
          IF ( ABS ( SQRT ( Bx**2 + By**2 + Bz**2)
     1             - 1.0 ) .GE. smallnum ) THEN
            WRITE ( 6, *)"Stopping at ",__LINE__, "in ",__FILE__
            CALL FLUSH ( 6)
            STOP
          END IF
#endif

! ------------------------------
! Build scaled relative coordinates from "izz" zone centroid and "izzp1" zone centroid to
! facial centroid. Only one Riemann solve at second order. (BTW, at higher orders
! the facial centroid would be replaced by looping over facial quadrature points.)

          x = ( zf_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zf_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zf_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! ----------

          x2 = ( zf_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixx, iyy, izzp1, 1) ) * obl02

          y2 = ( zf_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixx, iyy, izzp1, 2) ) * obl02

          z2 = ( zf_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixx, iyy, izzp1, 3) ) * obl02

! ----------
! Obtain the weight for this quadrature point. Weight is 1 at second order.

          wt = 1.0
#if (NUMORDER >= 3)
          wt = zf_quadwt_patch ( ixx, iyy, izz, j)
#endif

! ------------------------------
! Now use the above relative coordinates as well as
! "U ( ixx, iyy, izz, :, 1:4)" and "U ( ixx, iyy, izzp1, :, 1:4)"
! to build "U_L ( :), U_R ( :)" with second order accuracy.
! Note the use of "x, y, z" from bottom z-direction and
! "x2, y2, z2" from top z-direction.

          V_L ( :) = V ( ixx, iyy, izz, :, 1)
     1             + V ( ixx, iyy, izz, :, 2) * x
     1             + V ( ixx, iyy, izz, :, 3) * y
     1             + V ( ixx, iyy, izz, :, 4) * z
#if (NUMORDER >= 3)
     1             + V ( ixx, iyy, izz, :, 5)
     1             * (x**2 - Coef_Basis_Moments ( ixx, iyy, izz, 5))
     1             + V ( ixx, iyy, izz, :, 6)
     1             * (y**2 - Coef_Basis_Moments ( ixx, iyy, izz, 6))
     1             + V ( ixx, iyy, izz, :, 7)
     1             * (z**2 - Coef_Basis_Moments ( ixx, iyy, izz, 7))
     1             + V ( ixx, iyy, izz, :, 8)
     1             * (x * y - Coef_Basis_Moments ( ixx, iyy, izz, 8))
     1             + V ( ixx, iyy, izz, :, 9)
     1             * (y * z - Coef_Basis_Moments ( ixx, iyy, izz, 9))
     1             + V ( ixx, iyy, izz, :, 10)
     1             * (x * z - Coef_Basis_Moments ( ixx, iyy, izz, 10))
#elif (NUMORDER >= 4)
     1             + V ( ixx, iyy, izz, :, 11)
     1             * (x**3 - Coef_Basis_Moments ( ixx, iyy, izz, 11))
     1             + V ( ixx, iyy, izz, :, 12)
     1             * (y**3 - Coef_Basis_Moments ( ixx, iyy, izz, 12))
     1             + V ( ixx, iyy, izz, :, 13)
     1             * (z**3 - Coef_Basis_Moments ( ixx, iyy, izz, 13))
     1             + V ( ixx, iyy, izz, :, 14)
     1             * (x**2 * y
     1             - Coef_Basis_Moments ( ixx, iyy, izz, 14))
     1             + V ( ixx, iyy, izz, :, 15)
     1             * (x**2 * z
     1             - Coef_Basis_Moments ( ixx, iyy, izz, 15))
     1             + V ( ixx, iyy, izz, :, 16)
     1             * (x * y**2
     1             - Coef_Basis_Moments ( ixx, iyy, izz, 16))
     1             + V ( ixx, iyy, izz, :, 17)
     1             * (y**2 * z
     1             - Coef_Basis_Moments ( ixx, iyy, izz, 17))
     1             + V ( ixx, iyy, izz, :, 18)
     1             * (x * z**2
     1             - Coef_Basis_Moments ( ixx, iyy, izz, 18))
     1             + V ( ixx, iyy, izz, :, 19)
     1             * (y * z**2
     1             - Coef_Basis_Moments ( ixx, iyy, izz, 19))
     1             + V ( ixx, iyy, izz, :, 20)
     1             * (x * y * z
     1             - Coef_Basis_Moments ( ixx, iyy, izz, 20))
#endif

          V_R ( :) = V ( ixx, iyy, izzp1, :, 1)
     1             + V ( ixx, iyy, izzp1, :, 2) * x2
     1             + V ( ixx, iyy, izzp1, :, 3) * y2
     1             + V ( ixx, iyy, izzp1, :, 4) * z2
#if (NUMORDER >= 3)
     1             + V ( ixx, iyy, izzp1, :, 5)
     1             * (x2**2 - Coef_Basis_Moments ( ixx, iyy, izzp1, 5))
     1             + V ( ixx, iyy, izzp1, :, 6)
     1             * (y2**2 - Coef_Basis_Moments ( ixx, iyy, izzp1, 6))
     1             + V ( ixx, iyy, izzp1, :, 7)
     1             * (z2**2 - Coef_Basis_Moments ( ixx, iyy, izzp1, 7))
     1             + V ( ixx, iyy, izzp1, :, 8)
     1             * (x2 * y2
     1             - Coef_Basis_Moments ( ixx, iyy, izzp1, 8))
     1             + V ( ixx, iyy, izzp1, :, 9)
     1             * (y2 * z2
     1             - Coef_Basis_Moments ( ixx, iyy, izzp1, 9))
     1             + V ( ixx, iyy, izzp1, :, 10)
     1             * (x2 * z2
     1             - Coef_Basis_Moments ( ixx, iyy, izzp1, 10))
#elif (NUMORDER >= 4)
     1             + V ( ixx, iyy, izzp1, :, 11)
     1             * (x2**3 - Coef_Basis_Moments ( ixx, iyy, izzp1, 11))
     1             + V ( ixx, iyy, izzp1, :, 12)
     1             * (y2**3 - Coef_Basis_Moments ( ixx, iyy, izzp1, 12))
     1             + V ( ixx, iyy, izzp1, :, 13)
     1             * (z2**3 - Coef_Basis_Moments ( ixx, iyy, izzp1, 13))
     1             + V ( ixx, iyy, izzp1, :, 14)
     1             * (x2**2 * y2
     1             - Coef_Basis_Moments ( ixx, iyy, izzp1, 14))
     1             + V ( ixx, iyy, izzp1, :, 15)
     1             * (x2**2 * z2
     1             - Coef_Basis_Moments ( ixx, iyy, izzp1, 15))
     1             + V ( ixx, iyy, izzp1, :, 16)
     1             * (x2 * y2**2
     1             - Coef_Basis_Moments ( ixx, iyy, izzp1, 16))
     1             + V ( ixx, iyy, izzp1, :, 17)
     1             * (y2**2 * z2
     1             - Coef_Basis_Moments ( ixx, iyy, izzp1, 17))
     1             + V ( ixx, iyy, izzp1, :, 18)
     1             * (x2 * z2**2
     1             - Coef_Basis_Moments ( ixx, iyy, izzp1, 18))
     1             + V ( ixx, iyy, izzp1, :, 19)
     1             * (y2 * z2**2
     1             - Coef_Basis_Moments ( ixx, iyy, izzp1, 19))
     1             + V ( ixx, iyy, izzp1, :, 20)
     1             * (x2 * y2 * z2
     1             - Coef_Basis_Moments ( ixx, iyy, izzp1, 20))
#endif

! ------------------------------
! Now use the above relative coordinates as well as
! "entropygr ( ixx, iyy, izz, 1:4)" and "entropygr ( ixx, iyy, izzp1, 1:4)"
! to build "ent_L, ent_R" with second order accuracy.

#if ( PRESSURE_POSITIVITY >= 1)

          ent_L = entropygr ( ixx, iyy, izz, 1)
     1          + entropygr ( ixx, iyy, izz, 2) * x
     1          + entropygr ( ixx, iyy, izz, 3) * y
     1          + entropygr ( ixx, iyy, izz, 4) * z

          ent_R = entropygr ( ixx, iyy, izzp1, 1)
     1          + entropygr ( ixx, iyy, izzp1, 2) * x2
     1          + entropygr ( ixx, iyy, izzp1, 3) * y2
     1          + entropygr ( ixx, iyy, izzp1, 4) * z2

#else
          ent_L = prsfloor / rhofloor**(gamma - 1.0)
          ent_R = ent_L
#endif

! ------------------------------
! Notice that the 1D Riemann solver always assumes that the waves expand
! in the local x-direction. Therefore, we have to rotate the vectorial entities
! like velocities and magnetic fields into the local coordinate system.

          tempaa = Ax * V_L ( 2) + Ay * V_L ( 3) + Az * V_L ( 4)
          tempab = Bx * V_L ( 2) + By * V_L ( 3) + Bz * V_L ( 4)
          tempac = Cx * V_L ( 2) + Cy * V_L ( 3) + Cz * V_L ( 4)

          V_L ( 2) = tempaa
          V_L ( 3) = tempab
          V_L ( 4) = tempac

          tempaa = Ax * V_L ( 6) + Ay * V_L ( 7) + Az * V_L ( 8)
          tempab = Bx * V_L ( 6) + By * V_L ( 7) + Bz * V_L ( 8)
          tempac = Cx * V_L ( 6) + Cy * V_L ( 7) + Cz * V_L ( 8)

          V_L ( 6) = tempaa
          V_L ( 7) = tempab
          V_L ( 8) = tempac

! ----------

          tempaa = Ax * V_R ( 2) + Ay * V_R ( 3) + Az * V_R ( 4)
          tempab = Bx * V_R ( 2) + By * V_R ( 3) + Bz * V_R ( 4)
          tempac = Cx * V_R ( 2) + Cy * V_R ( 3) + Cz * V_R ( 4)

          V_R ( 2) = tempaa
          V_R ( 3) = tempab
          V_R ( 4) = tempac

          tempaa = Ax * V_R ( 6) + Ay * V_R ( 7) + Az * V_R ( 8)
          tempab = Bx * V_R ( 6) + By * V_R ( 7) + Bz * V_R ( 8)
          tempac = Cx * V_R ( 6) + Cy * V_R ( 7) + Cz * V_R ( 8)

          V_R ( 6) = tempaa
          V_R ( 7) = tempab
          V_R ( 8) = tempac

! --------------------------------------------------
! Obtain the background magnetic field B0 and rotate it.

          B0 ( :) = 0.0

#if ( USE_B0B1_FORMULATION == 1)

          x = zf_centroid_patch ( ixx, iyy, izz, 1)
          y = zf_centroid_patch ( ixx, iyy, izz, 2)
          z = zf_centroid_patch ( ixx, iyy, izz, 3)

          bflag = 1
          time_use = 0

! Calling acc-seq function to get the B0, using temp. variables x2, y2, z2.
          CALL GET_MAGNETIC_FIELD ( bflag, time_use,  x, y, z,
     1                              x2, y2, z2 )

          B0 ( 1) = Ax * x2 + Ay * y2 + Az * z2
          B0 ( 2) = Bx * x2 + By * y2 + Bz * z2
          B0 ( 3) = Cx * x2 + Cy * y2 + Cz * z2

#endif

! ------------------------------------------------------------
! Here, we apply non-reflective Riemann solver at an outflow boundary.
! Currently, this is done only for the MHD system. Later, we will extend it
! the CGL case.

! We set "iOriginal_RS" to 1 if we want to use Original RS. This is the default!

      iOriginal_RS = 1

! FOR NOW, SUCH A CRUDE CODE HACK IS ABSOLUTELY UNACCEPTABLE AND IS BEING REMOVED.
#if ( 1 == 1)

#if ( CHOOSE_APPL_HELIO == 1)

! --------------------------------------------------

#if ( IPDE_SYSTEM == 1)

      IF ( ( ibdy_outer == 1) .AND.
     1     ( izz .EQ. nr_zones) ) THEN

! -------

           c_sL = SQRT (gamma * V_L ( 5) / V_L (1 ) )

           v_anL = V_L ( 6) / SQRT ( fourpi * V_L (1 ) )

           v_aL = SQRT ( ( V_L (6)**2 + V_L (7)**2 + V_L (8)**2)
     1                 / ( fourpi * V_L (1 ) ) )

           d_L2 = SQRT ( ( v_aL**2 + c_sL**2 )**2
     1                 - 4.0*v_anL**2 * c_sL**2 )

           v_fnL = SQRT (0.5 * (v_aL**2 + c_sL**2 + d_L2))

! --------

           c_sR = SQRT (gamma * V_R ( 5) / V_R (1 ) )

           v_anR = V_R ( 6) / SQRT ( fourpi * V_R (1 ) )

           v_aR = SQRT ( ( V_R (6)**2 + V_R (7)**2 + V_R (8)**2)
     1                  / ( fourpi * V_R (1 ) ) )

           d_R2 = SQRT ( ( v_aR**2 + c_sR**2 )**2
     1                 - 4.0*v_anR**2 * c_sR**2 )

           v_fnR = SQRT (0.5 * (v_aR**2 + c_sR**2 + d_R2))


! -----------
        ! ---------------------------------
        ! 1. If an inflow, then skip this part
        ! ---------------------------------
        IF ( V_L ( 2) .LT. 0.0) THEN

           ! Do nothing

        ! ---------------------------------
        ! 2. Superfast outflow in the L-state
        ! ---------------------------------
        ELSE IF ( V_L (2) .GT. v_fnL ) THEN

           ! Do nothing

        ! ---------------------------------
        ! 3. Subfast outflow in both L and R states
        ! ---------------------------------
        ELSE IF ( V_R ( 2) .LT. v_fnR) THEN

           ! Do nothing

        ELSE
        ! Now, we modify the U_star state as given in the notes.

! -----------

        iOriginal_RS = 0

! -----------

        v_atL = SQRT ( V_L(7)**2 +  V_L(8)**2 ) / SQRT (fourpi*V_L(1))

        v_snL = SQRT (0.5 * (v_aL**2 + c_sL**2 - d_L2))

        alpha = c_sL*v_snL*SIGN(1.0, V_L(6))
     1        / (SQRT(fourpi)*V_L(1)**(2.0/3.0) *(c_sL**2-v_snL**2) )

        beta = c_sL**2/(V_L(1)*(c_sL**2 - v_snL**2))

        theta = (1.0/(2.0*d_L2))*( -v_fnL * v_aL**2
     1        + v_fnL * v_atL**2 * (c_sL**2 + v_snL**2)
     1          /(c_sL**2 - v_snL**2)
     1        + (gamma - 1.0)*v_fnL*c_sL**2
     1        - (gamma - 2.0)*c_sL*v_aL*v_snL)

! ------

        V_Star ( 1) = V_L(1) * (1.0 + (V_L(2) - v_fnL)/(theta+v_fnL))

        V_Star ( 2) = (theta*V_L(2) + v_fnL**2)/(theta+v_fnL)
        V_Star ( 3) = V_L(3) + alpha*V_L(7) * (V_Star ( 1) - V_L(1))
        V_Star ( 4) = V_L(4) + alpha*V_L(8) * (V_Star ( 1) - V_L(1))

        c_ss = c_sL * (1.0 + (gamma-1.0)*(V_Star(1)-V_L(1))
     1                     / (2.0*V_L(1)) )
        V_Star ( 5) = V_Star ( 1) * c_ss**2 / gamma

        V_Star ( 6) = V_L(6)
        V_Star ( 7) = V_L(7) * (1.0+beta*(V_Star(1)-V_L(1)))
        V_Star ( 8) = V_L(8) * (1.0+beta*(V_Star(1)-V_L(1)))

! ------

        flat = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                 flattengr ( ixx, iyy, izzp1) )

! ------

          CALL OneD_RS_MHD_Ptwise

     1         ( gamma, gamma_euler, isospeed, barotropic_stiffness,
     1           molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1           magB_small_cgl, eps_bdry_cgl,

     1           fourpi, ob_4pi, ob_8pi,
     1           smallnum, B0, flat,

     1           V_Star, V_Star, ent_L, ent_L,
     1           S_L, S_R, U_Star, F_Star, ent_F_Star,
     1           Fluct_L, Fluct_R,
     1           U_L, U_R, F_L, F_R)

! ------

        END IF ! End of "IF ( V_L ( 2) .LT. 0.0) THEN"

      END IF  ! End of "IF ( ( ibdy_outer == 1) .AND. ..."

#endif
! End of "#if ( IPDE_SYSTEM == 1)"

! --------------------------------------------------

#if ( IPDE_SYSTEM == 2)
! CGL system

      IF ( ( ibdy_outer == 1) .AND.
     1     ( izz .EQ. nr_zones) ) THEN

! -------

           prs_avg = ( V_L ( 5) + 2.0 * V_L ( 9) ) / 3.0

           c_sL = SQRT (gamma * prs_avg / V_L (1 ) )

           v_anL = V_L ( 6) / SQRT ( fourpi * V_L (1 ) )

           v_aL = SQRT ( ( V_L (6)**2 + V_L (7)**2 + V_L (8)**2)
     1                 / ( fourpi * V_L (1 ) ) )

           d_L2 = SQRT ( ( v_aL**2 + c_sL**2 )**2
     1                 - 4.0*v_anL**2 * c_sL**2 )

           v_fnL = SQRT (0.5 * (v_aL**2 + c_sL**2 + d_L2))

! --------

           prs_avg = ( V_R ( 5) + 2.0 * V_R ( 9) ) / 3.0

           c_sR = SQRT (gamma * prs_avg / V_R (1 ) )

           v_anR = V_R ( 6) / SQRT ( fourpi * V_R (1 ) )

           v_aR = SQRT ( ( V_R (6)**2 + V_R (7)**2 + V_R (8)**2)
     1                  / ( fourpi * V_R (1 ) ) )

           d_R2 = SQRT ( ( v_aR**2 + c_sR**2 )**2
     1                 - 4.0*v_anR**2 * c_sR**2 )

           v_fnR = SQRT (0.5 * (v_aR**2 + c_sR**2 + d_R2))

! -----------
        ! ---------------------------------
        ! 1. If an inflow, then skip this part
        ! ---------------------------------
        IF ( V_L ( 2) .LT. 0.0) THEN

           ! Do nothing

        ! ---------------------------------
        ! 2. Superfast outflow in the L-state
        ! ---------------------------------
        ELSE IF ( V_L (2) .GT. v_fnL ) THEN

           ! Do nothing

        ! ---------------------------------
        ! 3. Subfast outflow in both L and R states
        ! ---------------------------------
        ELSE IF ( V_R ( 2) .LT. v_fnR) THEN

           ! Do nothing

        ELSE
        ! Now, we modify the U_star state as given in the notes.

! -----------

        iOriginal_RS = 0

! -----------

        v_atL = SQRT ( V_L(7)**2 +  V_L(8)**2 ) / SQRT (fourpi*V_L(1))

        v_snL = SQRT (0.5 * (v_aL**2 + c_sL**2 - d_L2))

        alpha = c_sL*v_snL*SIGN(1.0, V_L(6))
     1        / (SQRT(fourpi)*V_L(1)**(2.0/3.0) *(c_sL**2-v_snL**2) )

        beta = c_sL**2/(V_L(1)*(c_sL**2 - v_snL**2))

        theta = (1.0/(2.0*d_L2))*( -v_fnL * v_aL**2
     1        + v_fnL * v_atL**2 * (c_sL**2 + v_snL**2)
     1          /(c_sL**2 - v_snL**2)
     1        + (gamma - 1.0)*v_fnL*c_sL**2
     1        - (gamma - 2.0)*c_sL*v_aL*v_snL)

! ------

        V_Star ( 1) = V_L(1) * (1.0 + (V_L(2) - v_fnL)/(theta+v_fnL))

        V_Star ( 2) = (theta*V_L(2) + v_fnL**2)/(theta+v_fnL)
        V_Star ( 3) = V_L(3) + alpha*V_L(7) * (V_Star ( 1) - V_L(1))
        V_Star ( 4) = V_L(4) + alpha*V_L(8) * (V_Star ( 1) - V_L(1))

        c_ss = c_sL * (1.0 + (gamma-1.0)*(V_Star(1)-V_L(1))
     1                     / (2.0*V_L(1)) )
        V_Star ( 5) = V_Star ( 1) * c_ss**2 / gamma
        V_Star ( 9) = V_Star ( 5)

        V_Star ( 6) = V_L(6)
        V_Star ( 7) = V_L(7) * (1.0+beta*(V_Star(1)-V_L(1)))
        V_Star ( 8) = V_L(8) * (1.0+beta*(V_Star(1)-V_L(1)))

        V_Star ( 10) = V_L ( 10)
        V_Star ( 11) = V_L ( 11)

! ------

        flat = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                 flattengr ( ixx, iyy, izzp1) )

! ------

          CALL OneD_RS_MHD_Ptwise

     1         ( gamma, gamma_euler, isospeed, barotropic_stiffness,
     1           molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1           magB_small_cgl, eps_bdry_cgl,

     1           fourpi, ob_4pi, ob_8pi,
     1           smallnum, B0, flat,

     1           V_Star, V_Star, ent_L, ent_L,
     1           S_L, S_R, U_Star, F_Star, ent_F_Star,
     1           Fluct_L, Fluct_R,
     1           U_L, U_R, F_L, F_R)

! ------

        END IF ! End of "IF ( V_L ( 2) .LT. 0.0) THEN"

      END IF  ! End of "IF ( ( ibdy_outer == 1) .AND. ..."

#endif
! End of "#if ( IPDE_SYSTEM == 2)"

#endif
! End of "#if ( CHOOSE_APPL_HELIO == 1)"

#endif
! FOR NOW, SUCH A CRUDE CODE HACK IS ABSOLUTELY UNACCEPTABLE AND IS BEING REMOVED.

! ------------------------------------------------------------

        flat = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                 flattengr ( ixx, iyy, izzp1) )

! ------------------------------------------------------------
! Now call the 1D Riemann solver. This RS operates in the local x-direction.
! It takes in "V_L, V_R" and "ent_L, ent_R" and returns speeds "S_L, S_R".
! Also returns resolved state "U_Star", resolved flux "F_Star", entropy flux "ent_F_Star".
! Also returns fluctuation "Fluct_L, Fluct_R" and "U_L, U_R, F_L, F_R".

      IF ( iOriginal_RS .EQ. 1) THEN

          CALL OneD_RS_MHD_Ptwise

     1         ( gamma, gamma_euler, isospeed, barotropic_stiffness,
     1           molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1           magB_small_cgl, eps_bdry_cgl,

     1           fourpi, ob_4pi, ob_8pi,
     1           smallnum, B0, flat,

     1           V_L, V_R, ent_L, ent_R,
     1           S_L, S_R, U_Star, F_Star, ent_F_Star,
     1           Fluct_L, Fluct_R,
     1           U_L, U_R, F_L, F_R)

      END IF

! ------------------------------------------------------------
! Hand back the resolved states, fluxes and extremal wave speeds from the
! 1D Riemann solver.

!          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) THEN
! If we get here, we are at an unphysical boundary, so set zero values and move on.

!            UStar_z ( ixx, iyy, izz, :) = 0.0
!            EntFlux_z ( ixx, iyy, izz) = 0.0

!            CYCLE  ! Go to the next iteration of the "DO ixx =" loop.

!          END IF

! ------------------------------
! If we get here, we are at a physical boundary, so assign the good values
! from the 1D dime-a-dozen Riemann solver.

! ----------
! If "USE_FACIAL_EFIELDS == 1" then "Ez = -F7" and "Ey = F8" in the local
! coordinates. Therefore, before we rotate back, we assign "Efx, Efy".

#if ( USE_FACIAL_EFIELDS == 1)
          IF ( j .EQ. 1) THEN

          v1 ( 1)   = 0.0
          v1 ( 2)   =   F_Star ( 8)
          v1 ( 3)   = - F_Star ( 7)

! Accumulate the upper edge.
          Efx ( ixx, iyy, izz) = Efx ( ixx, iyy, izz) + v1 ( 2)
          Efy ( ixx, iyy, izz) = Efy ( ixx, iyy, izz) + v1 ( 3)

! Count the number of accumulations.
          i_Efy ( ixx, iyy, izz) = i_Efy ( ixx, iyy, izz) + 1
          i_Efx ( ixx, iyy, izz) = i_Efx ( ixx, iyy, izz) + 1

! Assign directly for the lower edge.
          Efx_m1 ( ixx, iyy - 1, izz) = v1 ( 2)
          Efy_m1 ( ixx - 1, iyy, izz) = v1 ( 3)

          END IF
#endif

! ----------
! Save the extremal speeds. They will be needed in the 2D RS for electric field evaluation.
! This is only done for the face-centered quadrature point, because it has the best speeds.

          IF ( j .EQ. 1) THEN

          S_B_gr ( ixx, iyy, izz) = S_L

          S_T_gr ( ixx, iyy, izz) = S_R

          END IF

! ----------
! For CGL we have to use fluctuation form. For that, we return "Fluct_L, Fluct_R".
! These will also have to undergo the vectorial components in the same way.

! ----------
! Fluctuations in + direction; rotated into global frame.

          Fluct_Pl ( 1) = Fluct_R ( 1)                        ! Density

          Fluct_Pl ( 2) = Ax * Fluct_R ( 2)                   ! x-Momentum
     1               + Bx * Fluct_R ( 3) + Cx * Fluct_R ( 4)
          Fluct_Pl ( 3) = Ay * Fluct_R ( 2)                   ! y-Momentum
     1               + By * Fluct_R ( 3) + Cy * Fluct_R ( 4)
          Fluct_Pl ( 4) = Az * Fluct_R ( 2)                   ! z-Momentum
     1               + Bz * Fluct_R ( 3) + Cz * Fluct_R ( 4)

          Fluct_Pl ( 5) = Fluct_R ( 5)                        ! Energy

          Fluct_Pl ( 6) = Ax * Fluct_R ( 6)                   ! x-Magnetic field
     1               + Bx * Fluct_R ( 7) + Cx * Fluct_R ( 8)
          Fluct_Pl ( 7) = Ay * Fluct_R ( 6)                   ! y-Magnetic field
     1               + By * Fluct_R ( 7) + Cy * Fluct_R ( 8)
          Fluct_Pl ( 8) = Az * Fluct_R ( 6)                   ! z-Magnetic field
     1               + Bz * Fluct_R ( 7) + Cz * Fluct_R ( 8)

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3 )
          Fluct_Pl ( 9: ) = Fluct_R ( 9: )  ! pressures and tracer fields
#endif

#if ( IPDE_SYSTEM == 1 && NFLUID >= 1)
          Fluct_Pl ( NVAR + 1: NVAR + NFLUID) =  ! Species
     1     Fluct_R ( NVAR + 1: NVAR + NFLUID)
#endif

! -----

          Fluct_z_Pl ( ixx, iyy, izz, :) =Fluct_z_Pl ( ixx, iyy, izz, :)
     1    + wt * Fluct_Pl ( :)

! ----------
! Fluctuations in - direction; rotated into global frame.

          Fluct_Mi ( 1) = Fluct_L ( 1)                        ! Density

          Fluct_Mi ( 2) = Ax * Fluct_L ( 2)                   ! x-Momentum
     1               + Bx * Fluct_L ( 3) + Cx * Fluct_L ( 4)
          Fluct_Mi ( 3) = Ay * Fluct_L ( 2)                   ! y-Momentum
     1               + By * Fluct_L ( 3) + Cy * Fluct_L ( 4)
          Fluct_Mi ( 4) = Az * Fluct_L ( 2)                   ! z-Momentum
     1               + Bz * Fluct_L ( 3) + Cz * Fluct_L ( 4)

          Fluct_Mi ( 5) = Fluct_L ( 5)                        ! Energy

          Fluct_Mi ( 6) = Ax * Fluct_L ( 6)                   ! x-Magnetic field
     1               + Bx * Fluct_L ( 7) + Cx * Fluct_L ( 8)
          Fluct_Mi ( 7) = Ay * Fluct_L ( 6)                   ! y-Magnetic field
     1               + By * Fluct_L ( 7) + Cy * Fluct_L ( 8)
          Fluct_Mi ( 8) = Az * Fluct_L ( 6)                   ! z-Magnetic field
     1               + Bz * Fluct_L ( 7) + Cz * Fluct_L ( 8)

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3 )
          Fluct_Mi ( 9: ) = Fluct_L ( 9: )  ! pressures and tracer fields
#endif

#if ( IPDE_SYSTEM == 1 && NFLUID >= 1)
          Fluct_Mi ( NVAR + 1: NVAR + NFLUID) =  ! Species
     1     Fluct_L ( NVAR + 1: NVAR + NFLUID)
#endif

! -----

          Fluct_z_Mi ( ixx, iyy, izz, :) =Fluct_z_Mi ( ixx, iyy, izz, :)
     1    + wt * Fluct_Mi ( :)

! ----------
! Rotate F_L to F_LL; rotate F_R to F_RR and add in the internal terms to the update

          F_LL = F_L

          F_LL ( 2) = Ax * F_L ( 2) + Bx * F_L ( 3) + Cx * F_L ( 4)
          F_LL ( 3) = Ay * F_L ( 2) + By * F_L ( 3) + Cy * F_L ( 4)
          F_LL ( 4) = Az * F_L ( 2) + Bz * F_L ( 3) + Cz * F_L ( 4)

          F_LL ( 6) = Ax * F_L ( 6) + Bx * F_L ( 7) + Cx * F_L ( 8)
          F_LL ( 7) = Ay * F_L ( 6) + By * F_L ( 7) + Cy * F_L ( 8)
          F_LL ( 8) = Az * F_L ( 6) + Bz * F_L ( 7) + Cz * F_L ( 8)

! -----

          F_RR = F_R

          F_RR ( 2) = Ax * F_R ( 2) + Bx * F_R ( 3) + Cx * F_R ( 4)
          F_RR ( 3) = Ay * F_R ( 2) + By * F_R ( 3) + Cy * F_R ( 4)
          F_RR ( 4) = Az * F_R ( 2) + Bz * F_R ( 3) + Cz * F_R ( 4)

          F_RR ( 6) = Ax * F_R ( 6) + Bx * F_R ( 7) + Cx * F_R ( 8)
          F_RR ( 7) = Ay * F_R ( 6) + By * F_R ( 7) + Cy * F_R ( 8)
          F_RR ( 8) = Az * F_R ( 6) + Bz * F_R ( 7) + Cz * F_R ( 8)

! -----
! Now add in the contributions from the internal fluxes to the update terms in
! the zone "( ixx, iyy, izz, :)" and "( ixxp1, iyy, izz, :)".

          dU_dt ( ixx, iyy, izz, :) = dU_dt ( ixx, iyy, izz, :)
     1  - wt * ( Fluct_Mi ( :) + F_LL ( :) )
     1             * zf_area_patch ( ixx, iyy, izz)
     1             / v_volume_patch ( ixx, iyy, izz)

! -----

          dU_dt ( ixx, iyy, izzp1, :) = dU_dt ( ixx, iyy, izzp1, :)
     1  - wt * ( Fluct_Pl ( :) - F_RR ( :) )
     1             * zf_area_patch ( ixx, iyy, izz)
     1             / v_volume_patch ( ixx, iyy, izzp1)

! ----------
! Likewise, rotate the state components back into the global xyz-directions.

          IF ( j .EQ. 1) THEN

          UStar_z ( ixx, iyy, izz, 1) = U_Star ( 1)  ! Density


          UStar_z ( ixx, iyy, izz, 2) = Ax * U_Star ( 2)  ! Momentum density
     1             + Bx * U_Star ( 3) + Cx * U_Star ( 4)

          UStar_z ( ixx, iyy, izz, 3) = Ay * U_Star ( 2)
     1             + By * U_Star ( 3) + Cy * U_Star ( 4)

          UStar_z ( ixx, iyy, izz, 4) = Az * U_Star ( 2)
     1             + Bz * U_Star ( 3) + Cz * U_Star ( 4)

          UStar_z ( ixx, iyy, izz, 5) = U_Star ( 5)  ! Energy density

          UStar_z ( ixx, iyy, izz, 6) = Ax * U_Star ( 6)  ! Magnetic field
     1             + Bx * U_Star ( 7) + Cx * U_Star ( 8)

          UStar_z ( ixx, iyy, izz, 7) = Ay * U_Star ( 6)
     1             + By * U_Star ( 7) + Cy * U_Star ( 8)

          UStar_z ( ixx, iyy, izz, 8) = Az * U_Star ( 6)
     1             + Bz * U_Star ( 7) + Cz * U_Star ( 8)

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3 )
          UStar_z ( ixx, iyy, izz, 9: ) = U_Star ( 9: )  ! pressures and tracer fields
#endif

#if ( IPDE_SYSTEM == 1 && NFLUID >= 1)
          UStar_z ( ixx, iyy, izz, NVAR + 1: NVAR + NFLUID) =  ! Species density
     1                    U_Star ( NVAR + 1: NVAR + NFLUID)
#endif

! ----------
! Take care of entropy flux.

          EntFlux_z ( ixx, iyy, izz) =  ent_F_Star  ! Entropy density flux

          END IF  ! End of "If ( j .EQ. 1) THEN".

! ------------------------------------------------------------

          END DO  ! End of "DO ixx = 1, ntp_zones" loop.

        END DO  ! End of "DO iyy = 1, ntp_zones" loop.

      END DO  ! End of "DO izz = 0, nr_zones" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, NAREA_QUAD_PTS" loop over all facial quadrature points.

! ----------------------------------------------------------------------

! Make sure that the accumulation sum on "Efx, Efy" does not overwrite itself
! in a parallel setting.

#if ( USE_FACIAL_EFIELDS == 1)

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    iexpand,                                                      &
!$OMP&    Efx, Efy, Efx_m1, Efy_m1, i_Efx, i_Efy )                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present (                              &
!$acc&    iexpand, Efx, Efy, Efx_m1, Efy_m1, i_efx, i_Efy )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 0, nr_zones

        DO iyy = 1, ntp_zones

          DO ixx = 1, ntp_zones

! --------------------------------------------------

          Efx ( ixx, iyy - 1, izz) = Efx ( ixx, iyy - 1, izz)
     1                             + Efx_m1 ( ixx, iyy - 1, izz)

          Efy ( ixx - 1, iyy, izz) = Efy ( ixx - 1, iyy, izz)
     1                             + Efy_m1 ( ixx - 1, iyy, izz)


          i_Efx ( ixx, iyy - 1, izz) = i_Efx ( ixx, iyy - 1, izz) + 1

          i_Efy ( ixx - 1, iyy, izz) = i_Efy ( ixx - 1, iyy, izz) + 1

! --------------------------------------------------

          END DO  ! End of "DO ixx = 1, ntp_zones" loop.

        END DO  ! End of "DO iyy = 1, ntp_zones" loop.

      END DO  ! End of "DO izz = 0, nr_zones" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Divide the counter by the accumulation times. Since the accumulation
! has 0.25* multiplied there will be a multiplication by 4 is needed.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, smallnum,       &
!$OMP&    iexpand, Efx, Efy, Efz, i_Efx, i_Efy, i_Efz )                 &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present (                              &
!$acc&     iexpand, Efx, Efy, Efz, i_Efx, i_Efy, i_Efz )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 0, nr_zones

        DO iyy = 0, ntp_zones

          DO ixx = 0, ntp_zones

! --------------------------------------------------

          IF ( ixx .GE. 1)
     1    Efx ( ixx, iyy, izz) = Efx ( ixx, iyy, izz)
     1            / MAX( REAL( i_Efx ( ixx, iyy, izz)), smallnum)

          IF ( iyy .GE. 1)
     1    Efy ( ixx, iyy, izz) = Efy ( ixx, iyy, izz)
     1            / MAX( REAL( i_Efy ( ixx, iyy, izz)), smallnum)

          IF ( izz .GE. 1)
     1    Efz ( ixx, iyy, izz) = Efz ( ixx, iyy, izz)
     1            / MAX( REAL( i_Efz ( ixx, iyy, izz)), smallnum)

! --------------------------------------------------

          END DO  ! End of "DO ixx = 1 - iexpand, ntp_zones + iexpand" loop.

        END DO  ! End of "DO iyy = 1 - iexpand, ntp_zones + iexpand" loop.

      END DO  ! End of "DO izz = 1 - iexpand, nr_zones + iexpand" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! End of #if ( USE_FACIAL_EFIELDS == 1)
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 4) We are using a fluctuation form. In the previous steps, the contributions
! from the internal flux portions have already been added in to "dU_dt".
! However, now we have to add in the non-conservative parts of the characteristic
! matrices multiplied by the gradients.

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3 )

      iexpand = 2
!$acc update device( iexpand )

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz,                              &
!$OMP&    obl0, V_zc, dV_zc_dx, dV_zc_dy, dV_zc_dz,                     &
!$OMP&    bflag, time_use,  x, y, z, B0x, B0y, B0z,                     &
!$OMP&    dB0x_dx, dB0x_dy, dB0x_dz, dB0y_dx, dB0y_dy, dB0y_dz,         &
!$OMP&    dB0z_dx, dB0z_dy, dB0z_dz, nc_product)

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    rhofloor, prsfloor, smallnum, gamma, fourpi, ob_4pi, ob_8pi,  &
!$OMP&    isospeed, barotropic_stiffness, gamma_euler, iexpand,         &
!$OMP&    imask_active_patch, v_scale_length, V, v_centroid_patch,      &
!$OMP&    dU_dt)
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    rhofloor, prsfloor, smallnum, gamma, fourpi, ob_4pi, ob_8pi,  &
!$acc&    isospeed, barotropic_stiffness, gamma_euler, iexpand,         &
!$acc&    imask_active_patch, v_scale_length, V, v_centroid_patch,      &
!$acc&    dU_dt)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz,                                      &
!$acc&    obl0, V_zc, dV_zc_dx, dV_zc_dy, dV_zc_dz,                     &
!$acc&    bflag, time_use,  x, y, z, B0x, B0y, B0z,                     &
!$acc&    dB0x_dx, dB0x_dy, dB0x_dz, dB0y_dx, dB0y_dy, dB0y_dz,         &
!$acc&    dB0z_dx, dB0z_dy, dB0z_dz, nc_product)

      DO izz = 1 - iexpand, nr_zones + iexpand

        DO iyy = 1 - iexpand, ntp_zones + iexpand

          DO ixx = 1 - iexpand, ntp_zones + iexpand

! ------------------------------------------------------------
! The evaluation of "dV_zc_dx ( :), dV_zc_dy ( :), dV_zc_dz ( :)" is still only
! second order since it is evaluated only at zone-center. For higher order, this
! should be a weighted sum over the volumetric quadrature points.

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) CYCLE

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

          V_zc ( :) = V ( ixx, iyy, izz, :, 1)

          dV_zc_dx ( :) = V ( ixx, iyy, izz, :, 2) * obl0
          dV_zc_dy ( :) = V ( ixx, iyy, izz, :, 3) * obl0
          dV_zc_dz ( :) = V ( ixx, iyy, izz, :, 4) * obl0

! --------------------------------------------------
! Obtain the background magnetic field B0.

          B0x = 0.0;  B0y = 0.0;  B0z = 0.0;

#if ( USE_B0B1_FORMULATION == 1)

          x = v_centroid_patch ( ixx, iyy, izz, 1)
          y = v_centroid_patch ( ixx, iyy, izz, 2)
          z = v_centroid_patch ( ixx, iyy, izz, 3)

          bflag = 1
          time_use = 0

          CALL GET_MAGNETIC_FIELD_Gradients ( bflag, time_use,  x, y, z,
     1                              B0x, B0y, B0z,
     1                              dB0x_dx, dB0x_dy, dB0x_dz,
     1                              dB0y_dx, dB0y_dy, dB0y_dz,
     1                              dB0z_dx, dB0z_dy, dB0z_dz )

          V_zc ( 6) = V_zc ( 6) + B0x
          V_zc ( 7) = V_zc ( 7) + B0y
          V_zc ( 8) = V_zc ( 8) + B0z

          dV_zc_dx ( 6) = dV_zc_dx ( 6) + dB0x_dx
          dV_zc_dy ( 6) = dV_zc_dy ( 6) + dB0x_dy
          dV_zc_dz ( 6) = dV_zc_dz ( 6) + dB0x_dz

          dV_zc_dx ( 7) = dV_zc_dx ( 7) + dB0y_dx
          dV_zc_dy ( 7) = dV_zc_dy ( 7) + dB0y_dy
          dV_zc_dz ( 7) = dV_zc_dz ( 7) + dB0y_dz

          dV_zc_dx ( 8) = dV_zc_dx ( 8) + dB0z_dx
          dV_zc_dy ( 8) = dV_zc_dy ( 8) + dB0z_dy
          dV_zc_dz ( 8) = dV_zc_dz ( 8) + dB0z_dz

#endif

          CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_zc, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)


           dU_dt ( ixx, iyy, izz, 9) = dU_dt ( ixx, iyy, izz, 9)
     1                               - nc_product ( 1)

#if ( IPDE_SYSTEM == 3)
           dU_dt ( ixx, iyy, izz, 10) = dU_dt ( ixx, iyy, izz, 10)
     1                                - nc_product ( 2)
#endif

! ------------------------------------------------------------

          END DO  ! End of "DO ixx = 1, ntp_zones" loop.

        END DO  ! End of "DO iyy = 1, ntp_zones" loop.

      END DO  ! End of "DO izz = 0, nr_zones" loop.

!$acc end parallel
!$OMP END PARALLEL DO

#endif
! End of "( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3 )".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! The number of "end data" statements must match the "data copyin" statements.

!$acc end data

! ----------------------------------------------------------------------

      END SUBROUTINE MAKE_FACIAL_FLUXES







