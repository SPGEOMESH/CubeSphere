



#if ( NUMORDER == 3)



      SUBROUTINE LIMIT_ALL_SLOPES_P12

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           iz_patch, nr_patches,

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
     1         V, B1_zc, B2_zc, B3_zc, Ex_zc, Ey_zc, Ez_zc, Wzc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr, workgr,
     1         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,

     1         V_min_local, V_max_local)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! First zero's out the higher modes. Builds an approximate zone-centered
! magnetic field if appropriate. Also builds a flattener.

! Does the zone-centered limiting at second order; P1 polynomial reconstruction.
! Once that is done, it also builds the facial moments and does the
! face-centered limiting of the magnetic fields. We also build the
! divergence-free reconstruction for the magnetic fields.

! This subroutine is hardcoded for each order for the sake of efficiency.

! Uses 9 smaller piecewise linear stencils for the sake of efficiency. However,
! they are non-linearly hybridized with a large third order stencil so that the
! effective accuracy is third order.

! ------------------------------
! INDICES OF MOMENTS THAT HELP IN THE RECONSTRUCTION (set in BUILD_MODE_INDICES):-

! "ii_mode, jj_mode, kk_mode" Give us the powers of the Taylor series for
! volumetric reconstruction.

! ------------------------------
! CUBE SPHERE MESH VARIABLES:-

! ----------
! VARIABLES THAT GENERATE THE MESH (set in EXTRACT_VERT_COORDS_PATCH):-

! "vert_coords_patch" coordinates of the vertices of each patch. From this,
! all the other patch-based geometry information is generated via quadratic
! isoparametric mapping.

! "imask_active_patch" is 1 for good zones, 0 for bad ones.

! ----------
! FACIAL VARIABLES (set in COORDS_TO_WT_QUADPTS_FACES_ISO2):-

! "xf_quadpts_patch, yf_quadpts_patch, zf_quadpts_patch" locations of
! the quadrature points in the faces.

! "xf_quadwt_patch, yf_quadwt_patch, zf_quadwt_patch" weights associated
! with the quadrature points in the faces.

! "xf_quadnorml_patch, yf_quadnorml_patch, zf_quadnorml_patch" unit normal
! vector to that face. The z-face is curved and has a normal at each
! quadrature point.

! "xf_area_patch, yf_area_patch, zf_area_patch" the areas of the faces.

! "xf_centroid_patch, yf_centroid_patch, zf_centroid_patch" coordinates
! of the centroids of the faces.

! "v_scale_length" a scale length that tells us how big each zone is.

! ------------------------------
! MOMENTS FOR BASIS FUNCTIONS (set in INIT_BASIS_FUNCTION_COEFS):-

! "Coef_Basis_Moments" are the coefficients of the volumetric moments. This
! choice of Taylor series basis is useful for ensuring that all basis
! functions, except the zeroth order one, integrates to zero.

! ------------------------------
! VARIABLES THAT FACILITATE STENCIL OPERATIONS (set in STENCIL_INDEXING_TRICORNERS):-

! "itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU" (output) labels the
! four corners in the xy-plane. 1 if this is a tri-corner patch; 0 if this
! is a normal quad-corner patch.

! "istencil_patch ( :, :, 2, 12)" (output)
! "istencil_patch ( :, :, 2, 1:4)" gives us Right, Left, Up, Down zones in xy-plane.
! "istencil_patch ( :, :, 2, 5:8)" gives us the zones that are two ones over in the
! Right, Left, Up, Down directions.
! "istencil_patch ( :, :, 2, 9:12)" gives us the RU, LU, LD, RD zones.

! ------------------------------
! FLUID VARIABLES:-

! "U" comes in with just the zone-centered conserved variables and leaves with
! their moments included. The velocities are defined in a global Cartesian
! space, as are the magnetic fields.

! "V" is the vector of primitive variables that are built from the conserved
! variables. Useful for building flattener etc. variables.

! "Bxf, Byf, Bzf" are magnetic field components in the x,y,z-faces of the mesh.
! Includes moments in the face, that are filled in via reconstruction.

! "B1_zc, B2_zc, B3_zc" are the zone-centered, divergence-free reconstructed
! magnetic fields. They are defined in a global Cartesian space.

! "V_min_local, V_max_local" For this patch we evaluate min and max of local
! flow variables. Also keep track of divergence; which is why we have one
! more component.

! ------------------------------
! MODES OF FLUID VARIABLES:-

! U ( x, y, z) = u0 + ux * x + uy * y + uz * z

! Bxf ( x, y, z) = a0 + ax * x + ay * y + az * z

! Byf ( x, y, z) = b0 + bx * x + by * y + bz * z

! Bzf ( x, y, z) = c0 + cx * x + cy * y + cz * z

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes,

     1        iz_patch, nr_patches,

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
     1   Ex_zc, Ey_zc, Ez_zc, Wzc

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost) ::
     1      gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1      divvelgr, msonicgr, flattengr, workgr

      REAL, DIMENSION ( n_cc_components + 1) :: V_min_local, V_max_local

! ------------------------------
! Variables used for div-free reconstruction of magnetic fields:-

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           NBFLD_QUAD_PTS) ::
     1   Bxf_quadpts, Byf_quadpts, Bzf_quadpts

! ----------------------------------------------------------------------

! Variables past this point are internal to the subroutine.

      INTEGER ixx, iyy, izz, ixxm1, iyym1, izzm1, ixxp1, iyyp1, izzp1,
     1        ixxm2, iyym2, izzm2, ixxp2, iyyp2, izzp2,

     1        pow, icomp, irow, irr,

     1        ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,
     1        ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,
     1        ixxrr, iyyrr, ixxll, iyyll, ixxuu, iyyuu, ixxdd, iyydd,

     1        i, j, nn, ii, iorder

      REAL ob_8pi, ob_4pi, central_weight,
     1     prs_ratio_min, prs_ratio_max, tinynum, ratio_pres_pos,

     1     rhomin, rhomax, vlxmin, vlxmax, vlymin, vlymax,
     1     vlzmin, vlzmax, prsmin, prsmax, bfxmin, bfxmax,
     1     bfymin, bfymax, bfzmin, bfzmax, div_min, div_max,
     1     flat_min, flat_max, psi,

     1     pparlmin, pparlmax, pperpmin, pperpmax,
     1     pelectmin, pelectmax

#if( NFLUID <= 6)
      REAL sp1min, sp1max, sp2min, sp2max, sp3min, sp3max,
     1     sp4min, sp4max, sp5min, sp5max, sp6min, sp6max
#endif

! ----------------------------------------------------------------------
! Variables used for in-place, third order, WENO reconstruction.
! Entropy is still reconstructed in second order fashion, so we still need
! those above-mentioned second order variables too.

      REAL x, y, z, x2, y2, z2, obl0, obl02, ell_0, ell_j, nx, ny, nz,
     1     bx, by, bz, bx2, by2, bz2, bdiff,
     1     a11, a12, a13, a21, a22, a23, a31, a32, a33, obDet,
     1     b11, b12, b13, b21, b22, b23, b31, b32, b33,

     1     ISlo, IShi, Wtlo, Wthi,

     1     betalo, gammahi, gammalo, gammaratio, obgammahi, tau,

     1     ob12, ob180, ob144, E_pt,

     1     B0 ( 3)

! These matrices are needed for piecewise linear FV-WENO reconstruction of
! zone-centered densities. The linear variation is used in MRWENO style to
! eventually get piecewise quadratic reconstruction.

      REAL, DIMENSION ( 3, 8) :: R38

      REAL, DIMENSION ( 3, 9) :: X39

      REAL, DIMENSION ( 3, 3, 8) :: AInv338

      REAL, DIMENSION ( 9) :: IS9, Wt9

! These matrices are needed for piecewise quadratic FV-WENO reconstruction of
! the zone-centered densities.

      INTEGER indx_A9 ( 9), n_cc_10

      REAL, DIMENSION ( n_cc_modes) :: Coef_0, Coef_j

      REAL Amat99 ( 9, 9), Amat99_Inv ( 9, 9), vv_A9 ( 9),
     1     Amat_Row ( 9),

     1     Amat189 ( 18, 9), Amat189T ( 9, 18),
     1     R9 ( 9), X9 ( 9), Y9 ( 9), R18 ( 18), Phi9 ( 9)

! These matrices are needed for piecewise quadratic FD-WENO interpolation of
! the zone-centered electric fields.

      INTEGER indx_A10 ( 10)

      REAL Amat1010 ( 10, 10), Amat1010_Inv ( 10, 10), vv_A10 ( 10),
     1     Amat_Row10 ( 10),

     1     Amat1910 ( 19, 10), Amat1910T ( 10, 19),
     1     X10 ( 10), Y10 ( 10), R19 ( 19), Phi10 ( 10)

! ----------------------------------------------------------------------

      INTEGER iter, ipass, inode, bflag, iprs_flag

      REAL, DIMENSION ( n_cc_components) :: V_min, V_max, V_diff

      REAL v_cen ( 3), ent_min, ent_max, ent_diff, ff, omff, expandfac,
     1     aa, bb, cc, gamma_eos_pt,
     1     vx, vy, vz,

     1     Upt ( n_cc_components),
     1     Vpt ( n_cc_components), ent,

     1     Coef_Basis ( n_cc_modes),

     1     vert ( 3, 27), xyz_node ( 3, 27),
     1     U_node ( n_cc_components, 27),
     1     V_node ( n_cc_components, 27),
     1     ent_node ( 27),

     1     B0_node ( 3, 27), B0_x, B0_y, B0_z, time_use

! ----------------------------------------------------------------------
! Variables used for linear algebra that is needed in the
! zone-centered assignment of magnetic fields from the facial values.
! The first time we do this is based on LSQ and uses just minimal 2nd order.

      REAL nv ( 3),
     1     Am63 ( 6, 3),  Rm6 ( 6), Am33 ( 3, 3),
     1     ATm36 ( 3, 6), Rm3 ( 3), det, R1, R2, R3

      REAL A ( 5, 5), AInv ( 5, 5), R ( 5), Amat_5 ( 5, 5), R5 ( 5)

      INTEGER indx_B ( 5)
      REAL vv_B ( 5)

! ----------
! The full reconstruction that matches the facial values is a CLSQ procedure
! with "NConstraint = 6, NUnknown = 30, NLSQ = 64".

       INTEGER indx_CLSQ_B ( 30 + 6)

       REAL CMat_Row_B ( 30), wt

       REAL CMat_B ( 6, 30), DVec_B ( 6),
     1      AMat_B ( 64, 30), BVec_B ( 64),

     1      ATrans_B ( 30, 64), ATransA_B ( 30, 30),
     1      ATransB_B ( 30), CTrans_B ( 30, 6),

     1      CLSQMat_B ( 30 + 6, 30 + 6),
     1      CLSQVec_B ( 30 + 6),

     1      vv_CLSQ_B ( 30 + 6)

! ----------------------------------------------------------------------
! Index varibles, which helps in reconstruction using "Coef_Basis".
! It is mainly useful for the reconstruction of B1, B2, B3 where the
! sequence changes for higher order modes.

      INTEGER, PARAMETER ::
! 0, 1st, 2nd, 3rd order modes.
     1      iC0  = 1,
     1      iCx  = 2,  iCy  = 3,  iCz  = 4,
     1      iCxx = 5,  iCyy = 6,  iCzz = 7,
     1      iCxy = 8,  iCyz = 9,  iCxz = 10


! ----------------------------------------------------------------------

!$acc routine( MATSOLVE ) seq
!$acc routine( GET_MAGNETIC_FIELD ) seq
!$acc routine( CONS_TO_PRIM ) seq

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 1) Initialize various coefficients.

      ob_8pi = 1.0 / ( 8.0 * pi)
      ob_4pi = 1.0 / ( 4.0 * pi)

      ob12 = 1.0 / 12.0
      ob180 = 1.0 / 180.0
      ob144 = 1.0 / 144.0

! Weight of central zone relative to the weight of the other zones (unity).

      central_weight = 100.0

! The power to which the smoothness indicators are raised in the WENO
! reconstruction.

      pow = 4

! Set the minimum and maximum pressure ratios. These govern the settings
! for the scalar version of the flattener. I.e., there is a version of the
! flattener that does not rely on the divergence of the velocity; but only
! relies on the pressure ratio.

      prs_ratio_min = 1.4
      prs_ratio_max = 4.0

! Tiny number to prevent divide by zero specifically in the matrix inversion.

      tinynum = 1.0e-32

! ------------------------------
! Set the factor by which we want to expand the range of each variable.
! Setting this to 0 gives us the TVD range in multi-dimensions.
      expandfac = 0.2

! ------------------------------

! Move above initialization to GPU.

!$acc data copyin ( ob_8pi, ob_4pi, central_weight, pow, expandfac,     &
!$acc&              prs_ratio_min, prs_ratio_max, tinynum,              &
!$acc&              ratio_pres_pos, icomp, psi)

!$acc update device ( itri_patch_LD, itri_patch_RD, itri_patch_RU,      &
!$acc&                itri_patch_LU)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 2) For this very simple situation, we can initialize certain
! obvious scalars. Also zero out reconstructed variables. Make everything
! safe if we have a tri-corner.
! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    gamma, molewt, rhofloor, prsfloor, smallnum,                  &
!$OMP&    imask_active_patch, gamma_eos_gr, gamma_soundspeed_gr,        &
!$OMP&    mean_mol_wt_gr, msonicgr, U, V, entropygr, v_scale_length)    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)


!$acc parallel vector_length(64) present(                               &
!$acc&    gamma, molewt, rhofloor, prsfloor, smallnum,                  &
!$acc&    imask_active_patch, gamma_eos_gr, gamma_soundspeed_gr,        &
!$acc&    mean_mol_wt_gr, msonicgr, U, V, entropygr, v_scale_length)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)


      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          gamma_eos_gr ( ixx, iyy, izz) = gamma
          gamma_soundspeed_gr ( ixx, iyy, izz) = gamma
          mean_mol_wt_gr ( ixx, iyy, izz) = molewt

          msonicgr ( ixx, iyy, izz) = smallnum

          U ( ixx, iyy, izz, :, 2: ) = 0.0  ! Zero out any higher order zonal modes.
          entropygr ( ixx, iyy, izz, 2: ) = 0.0

          V ( ixx, iyy, izz, :, :) = 0.0

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) THEN
! If we get here, it means that this zone is in a tri-corner. Make it safe
! for nominal computation, i.e. prevent divide by zero.

          U ( ixx, iyy, izz, :, 1) = 0.0
          U ( ixx, iyy, izz, 1, 1) = rhofloor
          U ( ixx, iyy, izz, 5, 1) = prsfloor / ( gamma - 1.0)
#if ( IPDE_SYSTEM == 2)
          U ( ixx, iyy, izz, 9, 1) = prsfloor / ( gamma - 1.0)
#elif ( IPDE_SYSTEM == 3)
          U ( ixx, iyy, izz, 9: 10, 1) = prsfloor / ( gamma - 1.0)
#endif

          v_scale_length ( ixx, iyy, izz) = smallnum

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) THEN".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------
! Zero out higher order modes of "Bxf". If we are in a tri-corner, we zero
! out the entire "Bxf".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ii, ixxp1)                   &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    imask_active_patch, Bxf)                                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)


!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, Bxf)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ii, ixxp1)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 0 - NGhost, ntp_zones + NGhost

          ii = MAX0 ( ixx, 1 - NGhost)
          ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)

!          Bxf ( ixx, iyy, izz, 2:) = 0.0  ! Zero out any higher order facial modes.

          IF ( ( imask_active_patch ( ii, iyy, izz) .EQ. 0) .AND.
     1         ( imask_active_patch ( ixxp1, iyy, izz) .EQ. 0) )
     1    Bxf ( ixx, iyy, izz, 1) = 0.0

          END DO  ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------
! Zero out higher order modes of "Byf". If we are in a tri-corner, we zero
! out the entire "Byf".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ii, iyyp1)                   &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    imask_active_patch, Byf)                                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)


!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, Byf)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ii, iyyp1)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 0 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          ii = MAX0 ( iyy, 1 - NGhost)
          iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)

!          Byf ( ixx, iyy, izz, 2:) = 0.0  ! Zero out any higher order facial modes.

          IF ( ( imask_active_patch ( ixx, ii, izz) .EQ. 0) .AND.
     1         ( imask_active_patch ( ixx, iyyp1, izz) .EQ. 0) )
     1    Byf ( ixx, iyy, izz, 1) = 0.0

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------
! Zero out higher order modes of "Bzf". If we are in a tri-corner, we zero
! out the entire "Bzf".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    imask_active_patch, Bzf)                                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)


!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, Bzf)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 0 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

!          Bzf ( ixx, iyy, izz, 2:) = 0.0  ! Zero out any higher order facial modes.

! zones are masked along the "x" and "y" directions only and this is
! independent of the radial direction - "izz".
          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 0)
     1    Bzf ( ixx, iyy, izz, 1) = 0.0

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 0 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 3) If any of "itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU"
! are active, we have to synchronize the facial magnetic fields at the matching
! faces.
! ----------------------------------------------------------------------
! At this corner, we loop along the x-direction. The y-magnetic field components at
! the bottom of the zone should be synchronized with the field from the
! Down zone. Specifically, to the lower x-magnetic field components in the Down zone.
! The two fields will have opposite signs, so that has to be accounted for in
! the averaging that is needed for synchronization.

! This step is very intricate and easy to get wrong. So we don't reorder any
! loops here.

      IF ( itri_patch_LD .EQ. 1) THEN

      iyy = 1

!$OMP PARALLEL DO PRIVATE ( ixx, izz, ixxd, iyyd, a11)                  &

!$OMP&            SHARED ( iyy, nr_zones, ntp_zones, NGhost,            &
!$OMP&    imask_active_patch, istencil_patch, Bxf, Byf)                 &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, istencil_patch, Bxf, Byf)

!$acc loop gang independent                                             &
!$acc&    private ( ixx, izz, ixxd, iyyd, a11)

      DO ixx = 1 - NGhost, 0

      ixxd = istencil_patch ( ixx, iyy, 1, 4)
      iyyd = istencil_patch ( ixx, iyy, 2, 4)

!$acc loop vector
        DO izz = 1 - NGhost, nr_zones + NGhost

        a11 = 0.5 * ( Byf ( ixx, iyy - 1, izz, 1)
     1              - Bxf ( ixxd - 1, iyyd, izz, 1) )

        Byf ( ixx, iyy - 1, izz, 1) = a11

        Bxf ( ixxd - 1, iyyd, izz, 1) = - a11

        END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

      END DO  ! End of "DO ixx = 1 - NGhost, 0".

!$acc end parallel
!$OMP END PARALLEL DO

      END IF  ! End of "IF ( itri_patch_LD .EQ. 1) THEN".

! ------------------------------------------------------------
! At this corner, we loop along the x-direction. The y-magnetic field components at
! the bottom of the zone should be synchronized with the field from the
! Down zone. Specifically, to the upper x-magnetic field components in the Down zone.
! The two fields will have same signs, so that has to be accounted for in
! the averaging that is needed for synchronization.

      IF ( itri_patch_RD .EQ. 1) THEN

      iyy = 1

!$OMP PARALLEL DO PRIVATE ( ixx, izz, ixxd, iyyd, a11)                  &

!$OMP&            SHARED ( iyy, nr_zones, ntp_zones, NGhost,            &
!$OMP&    imask_active_patch, istencil_patch, Bxf, Byf)                 &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, istencil_patch, Bxf, Byf)

!$acc loop gang independent                                             &
!$acc&    private ( ixx, izz, ixxd, iyyd, a11)

      DO ixx = ntp_zones + 1, ntp_zones + NGhost

      ixxd = istencil_patch ( ixx, iyy, 1, 4)
      iyyd = istencil_patch ( ixx, iyy, 2, 4)

!$acc loop vector
        DO izz = 1 - NGhost, nr_zones + NGhost

        a11 = 0.5 * ( Byf ( ixx, iyy - 1, izz, 1)
     1              + Bxf ( ixxd, iyyd, izz, 1) )

        Byf ( ixx, iyy - 1, izz, 1) = a11

        Bxf ( ixxd, iyyd, izz, 1) = a11

        END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

      END DO  ! End of "DO ixx = ntp_zones + 1, ntp_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END IF  ! End of "IF ( itri_patch_RD .EQ. 1) THEN".

! ------------------------------------------------------------
! At this corner, we loop along the x-direction. The y-magnetic field components at
! the top of the zone should be synchronized with the field from the
! Up zone. Specifically, to the upper x-magnetic field components in the Up zone.
! The two fields will have opposite signs, so that has to be accounted for in
! the averaging that is needed for synchronization.

      IF ( itri_patch_RU .EQ. 1) THEN

      iyy = ntp_zones

!$OMP PARALLEL DO PRIVATE ( ixx, izz, ixxu, iyyu, a11)                  &

!$OMP&            SHARED ( iyy, nr_zones, ntp_zones, NGhost,            &
!$OMP&    imask_active_patch, istencil_patch, Bxf, Byf)                 &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, istencil_patch, Bxf, Byf)

!$acc loop gang independent                                             &
!$acc&    private ( ixx, izz, ixxu, iyyu, a11)

      DO ixx = ntp_zones + 1, ntp_zones + NGhost

      ixxu = istencil_patch ( ixx, iyy, 1, 3)
      iyyu = istencil_patch ( ixx, iyy, 2, 3)

!$acc loop vector
        DO izz = 1 - NGhost, nr_zones + NGhost

        a11 = 0.5 * ( Byf ( ixx, iyy, izz, 1)
     1              - Bxf ( ixxu, iyyu, izz, 1) )

        Byf ( ixx, iyy, izz, 1) = a11

        Bxf ( ixxu, iyyu, izz, 1) = - a11

        END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

      END DO  ! End of "DO ixx = ntp_zones + 1, ntp_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END IF  ! End of "IF ( itri_patch_RU .EQ. 1) THEN".

! ------------------------------------------------------------
! At this corner, we loop along the x-direction. The y-magnetic field components at
! the top of the zone should be synchronized with the field from the
! Up zone. Specifically, to the lower x-magnetic field components in the Up zone.
! The two fields will have same signs, so that has to be accounted for in
! the averaging that is needed for synchronization.

      IF ( itri_patch_LU .EQ. 1) THEN

      iyy = ntp_zones

!$OMP PARALLEL DO PRIVATE ( ixx, izz, ixxu, iyyu, a11)                  &

!$OMP&            SHARED ( iyy, nr_zones, ntp_zones, NGhost,            &
!$OMP&    imask_active_patch, istencil_patch, Bxf, Byf)                 &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, istencil_patch, Bxf, Byf)

!$acc loop gang independent                                             &
!$acc&    private ( ixx, izz, ixxu, iyyu, a11)

      DO ixx = 1 - NGhost, 0

      ixxu = istencil_patch ( ixx, iyy, 1, 3)
      iyyu = istencil_patch ( ixx, iyy, 2, 3)

!$acc loop vector
        DO izz = 1 - NGhost, nr_zones + NGhost

        a11 = 0.5 * ( Byf ( ixx, iyy, izz, 1)
     1              + Bxf ( ixxu - 1, iyyu, izz, 1) )

        Byf ( ixx, iyy, izz, 1) = a11

        Bxf ( ixxu - 1, iyyu, izz, 1) = a11

        END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

      END DO  ! End of "DO ixx = 1 - NGhost, 0" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END IF  ! End of "IF ( itri_patch_LU .EQ. 1) THEN".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
#if ( 0 == 1)
! This entire step is deprecated because the zone-centered magnetic field is
! set from its facial values in STEP 8.d.
! ----------------------------------------------------------------------
! STEP 4.a) Build estimate for the zone-centered magnetic field from the
! facial magnetic field values using the least squares minimization procedure.

! Since a zone-centered magnetic field "U ( ixx, iyy, izz, 6:8, 1)" is already
! available, for mild problems, this step may not be truly necessary.

! ----------------------------------------------------------------------
#if( 1 == 1)
! This is the better choice. It is exclusively based on least squares.
! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    ii, obl0, nv,                                                 &
!$OMP&    a11, a12, a13, a21, a22, a23, a31, a32, a33, obDet,           &
!$OMP&    b11, b12, b13, b21, b22, b23, b31, b32, b33,                  &
!$OMP&    det, Bx, By, Bz, R1, R2, R3,                                  &
!$OMP&    Am63, Rm6, ATm36, Am33, Rm3)                                  &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    ob_8pi, ob_4pi, smallnum, tinynum,                            &
!$OMP&    imask_active_patch, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc, U, V, &
!$OMP&    v_scale_length,                                               &
!$OMP&    xf_quadnorml_patch, yf_quadnorml_patch, zf_centroid_patch)    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    ob_8pi, ob_4pi, smallnum, tinynum,                            &
!$acc&    imask_active_patch, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc, U, V, &
!$acc&    v_scale_length,                                               &
!$acc&    xf_quadnorml_patch, yf_quadnorml_patch, zf_centroid_patch)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    ii, obl0, nv,                                                 &
!$acc&    a11, a12, a13, a21, a22, a23, a31, a32, a33, obDet,           &
!$acc&    b11, b12, b13, b21, b22, b23, b31, b32, b33,                  &
!$acc&    det, Bx, By, Bz, R1, R2, R3,                                  &
!$acc&    Am63, Rm6, ATm36, Am33, Rm3)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN

! ----------------------------------------------------------------------

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

          ii = ixx - (1 - NGhost) + 1

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

! ------------------------------
! Unit normals for the 6-x,y,z faces.
! x- and y- face unit normals are available from "x,yf_quadnorml_patch" array.

! z-face unit normal is calculated from radial unit vector.

          nv ( 1) = zf_centroid_patch ( ixx, iyy, izz, 1)
          nv ( 2) = zf_centroid_patch ( ixx, iyy, izz, 2)
          nv ( 3) = zf_centroid_patch ( ixx, iyy, izz, 3)

          a11 = 1.0 / AMAX1( smallnum, SQRT ( nv ( 1)**2 + nv ( 2)**2
     1                                                   + nv ( 3)**2))

          nv ( :) = nv ( :) * a11

! ----------

          Am63 ( 1, :) = xf_quadnorml_patch ( ixx  , iyy, izz, :)  ! x+ face
          Am63 ( 2, :) = xf_quadnorml_patch ( ixxm1, iyy, izz, :)  ! x- face

          Am63 ( 3, :) = yf_quadnorml_patch ( ixx, iyy  , izz, :)  ! y+ face
          Am63 ( 4, :) = yf_quadnorml_patch ( ixx, iyym1, izz, :)  ! y- face

          Am63 ( 5, :) = nv ( :)    ! z+ face
          Am63 ( 6, :) = nv ( :)    ! z- face

! ------------------------------
! RHS array "R6" consisting 6-facial magnetic field.

          Rm6 ( 1) = Bxf ( ixx  , iyy, izz, 1)
          Rm6 ( 2) = Bxf ( ixxm1, iyy, izz, 1)

          Rm6 ( 3) = Byf ( ixx, iyy  , izz, 1)
          Rm6 ( 4) = Byf ( ixx, iyym1, izz, 1)

          Rm6 ( 5) = Bzf ( ixx, iyy, izz  , 1)
          Rm6 ( 6) = Bzf ( ixx, iyy, izzm1, 1)


! ------------------------------
! Build final matrix and solve.

          ATm36 = TRANSPOSE ( Am63)

          Am33 = MATMUL ( ATm36, Am63)
          Rm3  = MATMUL ( ATm36, Rm6)

          a11 = Am33 ( 1, 1);   a12 = Am33 ( 1, 2);   a13 = Am33 ( 1, 3)
          a21 = Am33 ( 2, 1);   a22 = Am33 ( 2, 2);   a23 = Am33 ( 2, 3)
          a31 = Am33 ( 3, 1);   a32 = Am33 ( 3, 2);   a33 = Am33 ( 3, 3)

          R1 = Rm3 ( 1)
          R2 = Rm3 ( 2)
          R3 = Rm3 ( 3)

! ----------
! Obtain inverse of A33.

          det = a11 * ( a22 * a33 - a32 * a23)
     1        - a12 * ( a21 * a33 - a31 * a23)
     1        + a13 * ( a21 * a32 - a31 * a22)

          obdet = SIGN ( 1.0, det) / AMAX1 ( tinynum, ABS ( det) )

          b11 = ( a22 * a33 - a23 * a32) * obdet
          b12 = ( a13 * a32 - a12 * a33) * obdet
          b13 = ( a12 * a23 - a22 * a13) * obdet

          b21 = ( a23 * a31 - a21 * a33) * obdet
          b22 = ( a11 * a33 - a13 * a31) * obdet
          b23 = ( a13 * a21 - a11 * a23) * obdet

          b31 = ( a21 * a32 - a22 * a31) * obdet
          b32 = ( a12 * a31 - a11 * a32) * obdet
          b33 = ( a11 * a22 - a12 * a21) * obdet

! ------------------------------

          Bx = U ( ixx, iyy, izz, 6, 1)
          By = U ( ixx, iyy, izz, 7, 1)
          Bz = U ( ixx, iyy, izz, 8, 1)

          U ( ixx, iyy, izz, 6, 1) = b11 * R1 + b12 * R2 + b13 * R3
          U ( ixx, iyy, izz, 7, 1) = b21 * R1 + b22 * R2 + b23 * R3
          U ( ixx, iyy, izz, 8, 1) = b31 * R1 + b32 * R2 + b33 * R3

! ------------------------------
! Assign new magentic field from "Bn_zc"; this is also 2nd order accurate
! and consistent. But it violates energy conservation and should only be
! used in strong magnetic field situations where the pressure floor can be hit.

#if ( 1 == 1)
          U ( ixx, iyy, izz, 5, 1) = U ( ixx, iyy, izz, 5, 1)

     1        - ( Bx**2 + By**2 + Bz**2) * ob_8pi

     1        + ( U ( ixx, iyy, izz, 6, 1)**2
     1          + U ( ixx, iyy, izz, 7, 1)**2
     1          + U ( ixx, iyy, izz, 8, 1)**2) * ob_8pi
#endif

! ------------------------------

          B1_zc ( ixx, iyy, izz, 1) = U ( ixx, iyy, izz, 6, 1)
          B2_zc ( ixx, iyy, izz, 1) = U ( ixx, iyy, izz, 7, 1)
          B3_zc ( ixx, iyy, izz, 1) = U ( ixx, iyy, izz, 8, 1)

          V ( ixx, iyy, izz, 6, 1) = U ( ixx, iyy, izz, 6, 1)
          V ( ixx, iyy, izz, 7, 1) = U ( ixx, iyy, izz, 7, 1)
          V ( ixx, iyy, izz, 8, 1) = U ( ixx, iyy, izz, 8, 1)

! ----------------------------------------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
#else
! Also ok, but found to be not equally stable.
! ----------------------------------------------------------------------
! STEP 4.b) Build estimate for the zone-centered magnetic field that is at
! least second order accurate. This is done by requiring that a constrained
! div-free field matches the mean values at 5 out of the 6 faces of each zone.
! ------------------------------------------------------------

! We will be solving for a lot of 5X5 linear systems here.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, i, j, ii, nn,                &
!$OMP&    x, y, z, obl0, nx, ny, nz, a11, ixxm1, iyym1, izzm1, nt, nb,  &
!$OMP&    A, R, indx_B, vv_B, Amat_5, R5, vx, vy, vz, Bx, By, Bz)       &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    smallnum, tinynum, ob_8pi, rhofloor,                          &
!$OMP&    imask_active_patch, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc, U, V, &
!$OMP&    Ex_zc, Ey_zc, Ez_zc,                                          &
!$OMP&    xf_quadnorml_patch, yf_quadnorml_patch,                       &
!$OMP&    xf_centroid_patch, yf_centroid_patch, zf_centroid_patch,      &
!$OMP&    v_centroid_patch, v_scale_length )                            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum, tinynum, ob_8pi, rhofloor,                          &
!$acc&    imask_active_patch, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc, U, V, &
!$acc&    Ex_zc, Ey_zc, Ez_zc,                                          &
!$acc&    xf_quadnorml_patch, yf_quadnorml_patch,                       &
!$acc&    xf_centroid_patch, yf_centroid_patch, zf_centroid_patch,      &
!$acc&    v_centroid_patch, v_scale_length )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, i, j, ii, nn,                        &
!$acc&    x, y, z, obl0, nx, ny, nz, a11, ixxm1, iyym1, izzm1, nt, nb,  &
!$acc&    A, R, indx_B, vv_B, Amat_5, R5, vx, vy, vz, Bx, By, Bz)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

          ii = ixx - (1 - NGhost) + 1
! ------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

! ----------
! Face 1, top x-face with facial centroid "xf_centroid_patch ( ixx, iyy, izz, :)".

! Build scaled relative coordinates from zone centroid to facial centroid.

          x = ( xf_centroid_patch ( ixx, iyy, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xf_centroid_patch ( ixx, iyy, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xf_centroid_patch ( ixx, iyy, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Obtain the normal to that face.

          nx = xf_quadnorml_patch ( ixx, iyy, izz, 1)
          ny = xf_quadnorml_patch ( ixx, iyy, izz, 2)
          nz = xf_quadnorml_patch ( ixx, iyy, izz, 3)

! Assemble matrix elements for row 1 of the matrix.

          A ( 1, 1) = nx
          A ( 1, 2) = ny
          A ( 1, 3) = ny * y - nx * x
          A ( 1, 4) = nz
          A ( 1, 5) = nz * z - nx * x

! Assemble RHS for row 1

          R ( 1) = Bxf ( ixx, iyy, izz, 1)

! ----------
! Face 2, bottom x-face with facial centroid "xf_centroid_patch ( ixxm1, iyy, izz, :)".

! Build scaled relative coordinates from zone centroid to facial centroid.

          x = ( xf_centroid_patch ( ixxm1, iyy, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xf_centroid_patch ( ixxm1, iyy, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xf_centroid_patch ( ixxm1, iyy, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Obtain the normal to that face.

          nx = xf_quadnorml_patch ( ixxm1, iyy, izz, 1)
          ny = xf_quadnorml_patch ( ixxm1, iyy, izz, 2)
          nz = xf_quadnorml_patch ( ixxm1, iyy, izz, 3)

! Assemble matrix elements for row 2 of the matrix.

          A ( 2, 1) = nx
          A ( 2, 2) = ny
          A ( 2, 3) = ny * y - nx * x
          A ( 2, 4) = nz
          A ( 2, 5) = nz * z - nx * x

! Assemble RHS for row 2

          R ( 2) = Bxf ( ixxm1, iyy, izz, 1)

! ----------
! Face 3, top y-face with facial centroid "yf_centroid_patch ( ixx, iyy, izz, :)".

! Build scaled relative coordinates from zone centroid to facial centroid.

          x = ( yf_centroid_patch ( ixx, iyy, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yf_centroid_patch ( ixx, iyy, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yf_centroid_patch ( ixx, iyy, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Obtain the normal to that face.

          nx = yf_quadnorml_patch ( ixx, iyy, izz, 1)
          ny = yf_quadnorml_patch ( ixx, iyy, izz, 2)
          nz = yf_quadnorml_patch ( ixx, iyy, izz, 3)

! Assemble matrix elements for row 3 of the matrix.

          A ( 3, 1) = nx
          A ( 3, 2) = ny
          A ( 3, 3) = ny * y - nx * x
          A ( 3, 4) = nz
          A ( 3, 5) = nz * z - nx * x

! Assemble RHS for row 3

          R ( 3) = Byf ( ixx, iyy, izz, 1)

! ----------
! Face 4, bottom y-face with facial centroid "yf_centroid_patch ( ixx, iyym1, izz, :)".

! Build scaled relative coordinates from zone centroid to facial centroid.

          x = ( yf_centroid_patch ( ixx, iyym1, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yf_centroid_patch ( ixx, iyym1, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yf_centroid_patch ( ixx, iyym1, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Obtain the normal to that face.

          nx = yf_quadnorml_patch ( ixx, iyym1, izz, 1)
          ny = yf_quadnorml_patch ( ixx, iyym1, izz, 2)
          nz = yf_quadnorml_patch ( ixx, iyym1, izz, 3)

! Assemble matrix elements for row 4 of the matrix.

          A ( 4, 1) = nx
          A ( 4, 2) = ny
          A ( 4, 3) = ny * y - nx * x
          A ( 4, 4) = nz
          A ( 4, 5) = nz * z - nx * x

! Assemble RHS for row 4

          R ( 4) = Byf ( ixx, iyym1, izz, 1)

! ----------
! Face 5, top z-face with facial centroid "zf_centroid_patch ( ixx, iyy, izz, :)".

! Build scaled relative coordinates from zone centroid to facial centroid.

          x = ( zf_centroid_patch ( ixx, iyy, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zf_centroid_patch ( ixx, iyy, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zf_centroid_patch ( ixx, iyy, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Obtain the normal to that face. Basically, it is a unit normal in the radial
! direction for that curved face.

          nx = zf_centroid_patch ( ixx, iyy, izz, 1)
          ny = zf_centroid_patch ( ixx, iyy, izz, 2)
          nz = zf_centroid_patch ( ixx, iyy, izz, 3)

          a11 = 1.0 / AMAX1 ( smallnum, SQRT ( nx**2 + ny**2 + nz**2 ) )

          nx = nx * a11
          ny = ny * a11
          nz = nz * a11

! Assemble matrix elements for row 5 of the matrix.

          A ( 5, 1) = nx
          A ( 5, 2) = ny
          A ( 5, 3) = ny * y - nx * x
          A ( 5, 4) = nz
          A ( 5, 5) = nz * z - nx * x

! Assemble RHS for row 5

          R ( 5) = Bzf ( ixx, iyy, izz, 1)

! ------------------------------
! Now that "A" and "R" have been assembled, we can inline the linear system
! solve here. The linear solve is "A . X = R". "X" will be returned in "R".
! The explicit subroutine call has also been documented here.

! Even in tri-corners this logic gives us computable values.
          B1_zc ( ixx, iyy, izz, :) = 0.0
          B2_zc ( ixx, iyy, izz, :) = 0.0
          B3_zc ( ixx, iyy, izz, :) = 0.0

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

            Amat_5 ( 1:5, 1:5) = A ( 1:5, 1:5)

            R5 ( 1:5) = R ( 1:5)

            CALL MATSOLVE ( 5, Amat_5, R5, indx_B, vv_B)

            R ( 1:5) = R5 ( 1:5)

! ------------------------------
! In this, we scatter back the results from the linear system solve.
! We are only interested in the mean values, not the linear variation of
! Bx in the x-direction; By in the y-direction or Bz in the z-direction.

! ------------------------------
! We have made a quick and easy zone-centered B-field value which is
! at least second order accurate. We reset zone-centered magnetic field
! values here, if desired.

! Magnetic fields at zone centers evaluated consistent with facial
! values and div-free condition.

          B1_zc ( ixx, iyy, izz, 1) = R ( 1)
          B2_zc ( ixx, iyy, izz, 1) = R ( 2)
          B3_zc ( ixx, iyy, izz, 1) = R ( 4)

          B1_zc ( ixx, iyy, izz, 2) = - R ( 3) - R ( 5)
          B2_zc ( ixx, iyy, izz, 3) = R ( 3)
          B3_zc ( ixx, iyy, izz, 4) = R ( 5)

! ------------------------------
! Assign new magentic field from "Bn_zc"; this is also 2nd order accurate
! and consistent. But it violates energy conservation and should only be
! used in strong magnetic field situations where the pressure floor can be hit.

#if ( 1 == 1)
          U ( ixx, iyy, izz, 5, 1) = U ( ixx, iyy, izz, 5, 1)

     1        + ( B1_zc ( ixx, iyy, izz, 1)**2
     1          + B2_zc ( ixx, iyy, izz, 1)**2
     1          + B3_zc ( ixx, iyy, izz, 1)**2) * ob_8pi

     1        - ( U ( ixx, iyy, izz, 6, 1)**2
     1          + U ( ixx, iyy, izz, 7, 1)**2
     1          + U ( ixx, iyy, izz, 8, 1)**2) * ob_8pi
#endif

! ----------

          U ( ixx, iyy, izz, 6, 1) = B1_zc ( ixx, iyy, izz, 1)
          U ( ixx, iyy, izz, 7, 1) = B2_zc ( ixx, iyy, izz, 1)
          U ( ixx, iyy, izz, 8, 1) = B3_zc ( ixx, iyy, izz, 1)

          V ( ixx, iyy, izz, 6, 1) = B1_zc ( ixx, iyy, izz, 1)
          V ( ixx, iyy, izz, 7, 1) = B2_zc ( ixx, iyy, izz, 1)
          V ( ixx, iyy, izz, 8, 1) = B3_zc ( ixx, iyy, izz, 1)

! ------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

! --------------------------------------------------

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------
#endif
! End of choosing one or the other way of building zone-centered magnetic field.
! ----------------------------------------------------------------------
! This entire step is deprecated because the zone-centered magnetic field is
! set from its facial values in STEP 8.d.
#endif
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 5) From the conserved variables, we obtain the zone-centered primitive
! variables "V" that are at least second order accurate.
! Also build an estimate of the magnetosonic speed "msonicgr".
! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1, B0,     &
!$OMP&    a11, a12, a13, a21, x, y, z, obl0, iprs_flag, ent, Upt, Vpt)  &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,&
!$OMP&    rhofloor, prsfloor, tempfloor, smallnum, ob_4pi, ob_8pi,      &
!$OMP&    gasconst, molewt, gamma,                                      &
!$OMP&    isospeed, gamma_euler, barotropic_stiffness,                  &
!$OMP&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$OMP&    magB_small_cgl, eps_bdry_cgl,                                 &
!$OMP&    vert_coords_patch, U, V, B1_zc, B2_zc, B3_zc,                 &
!$OMP&    gamma_eos_gr, gamma_soundspeed_gr, entropygr, msonicgr)       &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    rhofloor, prsfloor, tempfloor, smallnum, ob_4pi, ob_8pi,      &
!$acc&    gasconst, molewt, gamma,                                      &
!$acc&    isospeed, gamma_euler, barotropic_stiffness,                  &
!$acc&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$acc&    magB_small_cgl, eps_bdry_cgl,                                 &
!$acc&    vert_coords_patch, U, V, B1_zc, B2_zc, B3_zc,                 &
!$acc&    gamma_eos_gr, gamma_soundspeed_gr, entropygr, msonicgr)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1, B0,             &
!$acc&    a11, a12, a13, a21, x, y, z, obl0, iprs_flag, ent, Upt, Vpt)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

! ------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

! ------------------------------

        B0 ( :) = 0.0

#if( USE_B0B1_FORMULATION == 1)
        x = v_centroid_patch ( ixx, iyy, izz, 1)
        y = v_centroid_patch ( ixx, iyy, izz, 2)
        z = v_centroid_patch ( ixx, iyy, izz, 3)

        bflag = 1         ! 1 == magnetic field;  0 == vector potential.
        time_use = 0.0    ! useful for the case of active rotation.

        CALL GET_MAGNETIC_FIELD ( bflag, time_use, x, y, z,
     1                              a11, a12, a13)

        B0 ( 1) = a11
        B0 ( 2) = a12
        B0 ( 3) = a13
#endif

! ------------------------------
! The call to "CONS_TO_PRIM" will also bring all variables within their physical
! ranges. As a result, after this call we feel free to reset "U, V, ent".

! If "PRESSURE_POSITIVITY > 0", "CONS_TO_PRIM" can even use "ent" to reset "U, V"
! in extreme situations.

          Upt ( :) = U ( ixx, iyy, izz, :, 1)

          ent = entropygr ( ixx, iyy, izz, 1)

          CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, Upt, Vpt)

          V ( ixx, iyy, izz, :, 1) = Vpt ( :)
          U ( ixx, iyy, izz, :, 1) = Upt ( :)
          entropygr ( ixx, iyy, izz, 1) = ent

! --------------------
! Take care of higher modes of the entropy variable. Obtain species fractions.
! Evaluate the magnetosonic speed.

          entropygr ( ixx, iyy, izz, 2:) = 0.0  ! Zero out gradients of entropy.

#if ( NFLUID >= 1)
          V ( ixx, iyy, izz, NVAR + 1: , 1) =
     1    U ( ixx, iyy, izz, NVAR + 1: , 1) / V ( ixx, iyy, izz, 1, 1)
#endif
#if ( NVAR > 10)
          V ( ixx, iyy, izz, 10 + 1: , 1) =
     1    U ( ixx, iyy, izz, 10 + 1: , 1) / V ( ixx, iyy, izz, 1, 1)
#endif

          a12 = ( ( V ( ixx, iyy, izz, 6, 1) + B0 ( 1))**2
     1          + ( V ( ixx, iyy, izz, 7, 1) + B0 ( 2))**2
     1          + ( V ( ixx, iyy, izz, 8, 1) + B0 ( 3))**2) * ob_4pi

! --------------------

#if (IPDE_SYSTEM == 1)
#if ( CHOOSE_ADIABATIC == 1)
! Using adiabatic EOS.

          msonicgr ( ixx, iyy, izz) =
     1    gamma_soundspeed_gr ( ixx, iyy, izz)* V ( ixx, iyy, izz, 5, 1)
     1  + a12

          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz)
     1                                     / V ( ixx, iyy, izz, 1, 1) )

#else
! Using isothermal EOS.

          msonicgr ( ixx, iyy, izz) = isospeed**2
     1     + gamma_euler * barotropic_stiffness
     1                   * V ( ixx, iyy, izz, 1, 1)**(gamma_euler - 1.0)
     1     + a12 / V ( ixx, iyy, izz, 1, 1)

          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz))

#endif
#endif
! End of "IPDE_SYSTEM == 1"

! --------------------

#if (IPDE_SYSTEM == 2)
          msonicgr ( ixx, iyy, izz) = 1.67 * ( V ( ixx, iyy, izz, 5, 1)
     1                          + 2.0 * V ( ixx, iyy, izz, 9, 1)) / 3.0
     1  + a12

          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz)
     1                                     / V ( ixx, iyy, izz, 1, 1) )
#endif

! --------------------

#if (IPDE_SYSTEM == 3)
          msonicgr ( ixx, iyy, izz) = 1.67 * ( V ( ixx, iyy, izz, 5, 1)
     1                          + 2.0 * V ( ixx, iyy, izz, 9, 1)) / 3.0
     1                              + 1.67 * V ( ixx, iyy, izz, 10, 1)
     1  + a12

          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz)
     1                                     / V ( ixx, iyy, izz, 1, 1) )
#endif

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------






! ----------------------------------------------------------------------
! Remove this for normal code; only needed for heliospheric wind code.
#if (1 == 1)
! ----------------------------------------------------------------------

#if ( CHOOSE_APPL_HELIO == 1)

      IF ( iz_patch .EQ. 1 ) THEN

! To prevent pressure in the first radial zone from sinking for the
! heliosphere problem, we provide a floor.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost,                &
!$OMP&    V, U, gamma, imask_active_patch, ob_8pi)                      &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)
!$acc parallel vector_length(64) present (                              &
!$acc&    V, U, gamma, imask_active_patch, ob_8pi)
!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

      izz = 1

      IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! --------------------------------------------------
! This value of "3.8914E-5" is actually "prs_sw" evaluated at the inner boundary
! of the heliosphere boundary. For CGL, we will have to fix two such pressures.

! ------------------------------
#if (IPDE_SYSTEM == 1)
! ------------------------------

! For MHD, we have a single pressure in the 5th component of "V" and a single
! total Energy in "U".

!        IF ( V ( ixx, iyy, izz, 5, 1) .LT. 3.8914E-5) THEN

        U ( ixx, iyy, izz, 5, 1) = U ( ixx, iyy, izz, 5, 1)
     1   + (3.8914E-5 - V ( ixx, iyy, izz, 5, 1)) / (gamma - 1.0)

        V ( ixx, iyy, izz, 5, 1) = 3.8914E-5

!        END IF

! ------------------------------
#elif (IPDE_SYSTEM == 2)
! ------------------------------

! For CGL we have "V (5) = P_parallel ; V (9) = P_perpendicular" and
! "U (5) = Total energy density ; U (9) = P_parallel - P_perpendicular".
! In the setup we want "P_parallel == P_perpendicular" exactly.

        V ( ixx, iyy, izz, 5, 1) = 3.8914E-5

        V ( ixx, iyy, izz, 9, 1) = 3.8914E-5

        U ( ixx, iyy, izz, 5, 1) = 1.5 * 3.8914E-5

     1+ 0.5 * V ( ixx, iyy, izz, 1, 1) * ( V ( ixx, iyy, izz, 2, 1)**2
     1  + V ( ixx, iyy, izz, 3, 1)**2 + V ( ixx, iyy, izz, 4, 1)**2 )

     1 + ( V ( ixx, iyy, izz, 6, 1)**2 + V ( ixx, iyy, izz, 7, 1)**2
     1   + V ( ixx, iyy, izz, 8, 1)**2 ) * ob_8pi

        U ( ixx, iyy, izz, 9, 1) = 0.0  ! Pressures are isotropic.

#if ( NVAR >= 10)
        V ( ixx, iyy, izz, 10, 1) = 1.0
#endif

#if ( NVAR >= 11)
        V ( ixx, iyy, izz, 11: , 1) = 0.0
#endif

! ------------------------------
#elif (IPDE_SYSTEM == 3)
! ------------------------------

! For CGL with electron pressure we have "V (5) = P_parallel ; V (9) = P_perpendicular",
! "V (10) = Electron_Pressure" and "U (5) = Total energy density",
! "U (9) = P_parallel - P_perpendicular", "U (10) = Density * Electron_Pressure"

        V ( ixx, iyy, izz, 5, 1) = 3.8914E-5

        V ( ixx, iyy, izz, 9, 1) = 3.8914E-5

        V ( ixx, iyy, izz, 10, 1) = 3.8914E-5

#if ( NVAR >= 11)
        V ( ixx, iyy, izz, 11, 1) = 1.0
#endif

#if ( NVAR >= 12)
        V ( ixx, iyy, izz, 12: , 1) = 0.0
#endif


        U ( ixx, iyy, izz, 5, 1) = 1.5 * 3.8914E-5 + 1.5 * 3.8914E-5

     1+ 0.5 * V ( ixx, iyy, izz, 1, 1) * ( V ( ixx, iyy, izz, 2, 1)**2
     1  + V ( ixx, iyy, izz, 3, 1)**2 + V ( ixx, iyy, izz, 4, 1)**2 )

     1 + ( V ( ixx, iyy, izz, 6, 1)**2 + V ( ixx, iyy, izz, 7, 1)**2
     1   + V ( ixx, iyy, izz, 8, 1)**2 ) * ob_8pi

        U ( ixx, iyy, izz, 9, 1) = 0.0  ! Ion Pressures are isotropic at start.

        U ( ixx, iyy, izz, 10, 1) = V ( ixx, iyy, izz, 1, 1) * 3.8914E-5

#if ( NVAR >= 11)
        U ( ixx, iyy, izz, 11, 1) = V ( ixx, iyy, izz, 1, 1) * 1.0
#endif

#if ( NVAR >= 12)
        U ( ixx, iyy, izz, 12: , 1) = 0.0
#endif

! ------------------------------
#endif
! End of "IPDE_SYSTEM == 1/2/3".

! --------------------------------------------------

      END IF    ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN".

      END DO    ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO    ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
!$acc end parallel
!$OMP END PARALLEL DO

      END IF    ! End of "IF ( iz_patch .EQ. 1 ) THEN".

#endif
! End of "CHOOSE_APPL_HELIO == 1".

! ----------------------------------------------------------------------
#endif
! Remove this for normal code; only needed for heliospheric wind code.
! ----------------------------------------------------------------------






! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 6) Now build the "divvelgr" and then use it to build "flattengr".
! This is done in a sequence of several loops.
! ----------------------------------------------------------------------

!$acc kernels
      flattengr ( :, :, :) = 0.0
!$acc end kernels

! "flattengr" applies a flattener in regions of strong shock.
! ("flattengr" is 0.0 near smooth flow and 1.0 in regions with v. strong
!   shocks and rarefactions.)

! ------------------------------------------------------------
! STEP 6.a) Loop over zones that are within three zone widths of the boundary of the
! dynamically active part of the mesh. This is true when "NGhost == 4". When that
! is not the case, we focus on "NGhost - 1" zone widths of the boundary.
! ------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    ixxp1, iyyp1, izzp1, a11, a12, a13, x, y, z, obl0, nx, ny, nz,&
!$OMP&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd)               &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, prsfloor,       &
!$OMP&    smallnum,                                                     &
!$OMP&    ob_8pi, ob_4pi, prs_ratio_max, prs_ratio_min, flatten_coef,   &
!$OMP&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$OMP&    xf_quadnorml_patch, yf_quadnorml_patch, zf_centroid_patch,    &
!$OMP&    istencil_patch,                                               &
!$OMP&    V, divvelgr, workgr, msonicgr, flattengr)                     &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    prsfloor, smallnum,                                           &
!$acc&    ob_8pi, ob_4pi, prs_ratio_max, prs_ratio_min, flatten_coef,   &
!$acc&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$acc&    xf_quadnorml_patch, yf_quadnorml_patch, zf_centroid_patch,    &
!$acc&    istencil_patch,                                               &
!$acc&    V, divvelgr, workgr, msonicgr, flattengr)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    ixxp1, iyyp1, izzp1, a11, a12, a13, x, y, z, obl0, nx, ny, nz,&
!$acc&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd)

      DO izz = 2 - NGhost, nr_zones + NGhost - 1

        DO iyy = 2 - NGhost, ntp_zones + NGhost - 1

          DO ixx = 2 - NGhost, ntp_zones + NGhost - 1

          izzm1 = izz - 1
          izzp1 = izz + 1

          iyym1 = iyy - 1
          iyyp1 = iyy + 1

          ixxm1 = ixx - 1
          ixxp1 = ixx + 1

          ixxr = istencil_patch ( ixx, iyy, 1, 1)  ! For indexing zone one Right.
          iyyr = istencil_patch ( ixx, iyy, 2, 1)

          ixxl = istencil_patch ( ixx, iyy, 1, 2)  ! For indexing zone one Left.
          iyyl = istencil_patch ( ixx, iyy, 2, 2)

          ixxu = istencil_patch ( ixx, iyy, 1, 3)  ! For indexing zone one Up.
          iyyu = istencil_patch ( ixx, iyy, 2, 3)

          ixxd = istencil_patch ( ixx, iyy, 1, 4)  ! For indexing zone one Down.
          iyyd = istencil_patch ( ixx, iyy, 2, 4)

! ----------------------------------------

          obl0 = 1.0
     1    / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

          divvelgr ( ixx, iyy, izz) = 0.0

! Find the minimal magnetosonic speed in a halo of zones around the current one.

          workgr ( ixx, iyy, izz) = AMIN1 ( msonicgr ( ixx, iyy, izz),
     1    msonicgr ( ixxr, iyyr, izz), msonicgr ( ixxl, iyyl, izz),
     1    msonicgr ( ixxu, iyyu, izz), msonicgr ( ixxd, iyyd, izz),
     1    msonicgr ( ixx, iyy, izzp1), msonicgr ( ixx, iyy, izzm1) )

! ------------------------------

! ----------
! Face 1, top x-face with normal "xf_quadnorml_patch ( ixx, iyy, izz, :)".
! Physical zone to the Right "( ixxr, iyyr, izz)".

! Build scaled relative coordinates from zone centroid to zone centroid.
! Then scale that with the scaling length for these two zones.

          x = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          a11 = 0.5 / AMAX1 ( SQRT ( x**2 + y**2 + z**2 ), smallnum)

! Obtain the normal to that face.

          nx = xf_quadnorml_patch ( ixx, iyy, izz, 1)
          ny = xf_quadnorml_patch ( ixx, iyy, izz, 2)
          nz = xf_quadnorml_patch ( ixx, iyy, izz, 3)

! Project the velocity along the normal and scale to the size of this zone.

          divvelgr ( ixx, iyy, izz) = divvelgr ( ixx, iyy, izz)

     1    + nx * a11 * ( V ( ixxr, iyyr, izz, 2, 1)
     1                 - V ( ixx, iyy, izz, 2, 1) )

     1    + ny * a11 * ( V ( ixxr, iyyr, izz, 3, 1)
     1                 - V ( ixx, iyy, izz, 3, 1) )

     1    + nz * a11 * ( V ( ixxr, iyyr, izz, 4, 1)
     1                 - V ( ixx, iyy, izz, 4, 1) )

! ----------
! Face 2, bottom x-face with normal "xf_quadnorml_patch ( ixxm1, iyy, izz, :)".
! Physical zone to the Left "( ixxl, iyyl, izz)".

! Build scaled relative coordinates from zone centroid to zone centroid.
! Then scale that with the scaling length for these two zones.

          x = ( v_centroid_patch ( ixx, iyy, izz, 1)
     1      - v_centroid_patch ( ixxl, iyyl, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixx, iyy, izz, 2)
     1      - v_centroid_patch ( ixxl, iyyl, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixx, iyy, izz, 3)
     1      - v_centroid_patch ( ixxl, iyyl, izz, 3) ) * obl0

          a11 = 0.5 / AMAX1 ( SQRT ( x**2 + y**2 + z**2 ), smallnum)

! Obtain the normal to that face.

          nx = xf_quadnorml_patch ( ixxm1, iyy, izz, 1)
          ny = xf_quadnorml_patch ( ixxm1, iyy, izz, 2)
          nz = xf_quadnorml_patch ( ixxm1, iyy, izz, 3)

! Project the velocity along the normal and scale to the size of this zone.

          divvelgr ( ixx, iyy, izz) = divvelgr ( ixx, iyy, izz)

     1    + nx * a11 * ( V ( ixx, iyy, izz, 2, 1)
     1                 - V ( ixxl, iyyl, izz, 2, 1) )

     1    + ny * a11 * ( V ( ixx, iyy, izz, 3, 1)
     1                 - V ( ixxl, iyyl, izz, 3, 1) )

     1    + nz * a11 * ( V ( ixx, iyy, izz, 4, 1)
     1                 - V ( ixxl, iyyl, izz, 4, 1) )

! ----------
! Face 3, top y-face with normal "yf_quadnorml_patch ( ixx, iyy, izz, :)".
! Physical zone to the Up "( ixxu, iyyu, izz)".

! Build scaled relative coordinates from zone centroid to zone centroid.
! Then scale that with the scaling length for these two zones.

          x = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          a11 = 0.5 / AMAX1 ( SQRT ( x**2 + y**2 + z**2 ), smallnum)

! Obtain the normal to that face.

          nx = yf_quadnorml_patch ( ixx, iyy, izz, 1)
          ny = yf_quadnorml_patch ( ixx, iyy, izz, 2)
          nz = yf_quadnorml_patch ( ixx, iyy, izz, 3)

! Project the velocity along the normal and scale to the size of this zone.

          divvelgr ( ixx, iyy, izz) = divvelgr ( ixx, iyy, izz)

     1    + nx * a11 * ( V ( ixxu, iyyu, izz, 2, 1)
     1                 - V ( ixx, iyy, izz, 2, 1) )

     1    + ny * a11 * ( V ( ixxu, iyyu, izz, 3, 1)
     1                 - V ( ixx, iyy, izz, 3, 1) )

     1    + nz * a11 * ( V ( ixxu, iyyu, izz, 4, 1)
     1                 - V ( ixx, iyy, izz, 4, 1) )

! ----------
! Face 4, bottom y-face with normal "yf_quadnorml_patch ( ixx, iyym1, izz, :)".
! Physical zone to the Down "( ixxd, iyyd, izz)".

! Build scaled relative coordinates from zone centroid to zone centroid.
! Then scale that with the scaling length for these two zones.

          x = ( v_centroid_patch ( ixx, iyy, izz, 1)
     1      - v_centroid_patch ( ixxd, iyyd, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixx, iyy, izz, 2)
     1      - v_centroid_patch ( ixxd, iyyd, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixx, iyy, izz, 3)
     1      - v_centroid_patch ( ixxd, iyyd, izz, 3) ) * obl0

          a11 = 0.5 / AMAX1 ( SQRT ( x**2 + y**2 + z**2 ), smallnum)

! Obtain the normal to that face.

          nx = yf_quadnorml_patch ( ixx, iyym1, izz, 1)
          ny = yf_quadnorml_patch ( ixx, iyym1, izz, 2)
          nz = yf_quadnorml_patch ( ixx, iyym1, izz, 3)

! Project the velocity along the normal and scale to the size of this zone.

          divvelgr ( ixx, iyy, izz) = divvelgr ( ixx, iyy, izz)

     1    + nx * a11 * ( V ( ixx, iyy, izz, 2, 1)
     1                 - V ( ixxd, iyyd, izz, 2, 1) )

     1    + ny * a11 * ( V ( ixx, iyy, izz, 3, 1)
     1                 - V ( ixxd, iyyd, izz, 3, 1) )

     1    + nz * a11 * ( V ( ixx, iyy, izz, 4, 1)
     1                 - V ( ixxd, iyyd, izz, 4, 1) )

! ----------
! Face 5, top z-face with normal in the direction "zf_centroid_patch ( ixx, iyy, izz, :)".
! Physical zone in the top radial direction "( ixx, iyy, izzp1)".

! Top and bottom z-faces will have normal in the same (radial) direction.

! Build scaled relative coordinates from zone centroid to zone centroid.
! Then scale that with the scaling length for these two zones.

          x = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          a11 = 0.5 / AMAX1 ( SQRT ( x**2 + y**2 + z**2 ), smallnum)

! Obtain the normal to that face. Basically, it is a unit normal in the radial
! direction for that curved face.

          nx = zf_centroid_patch ( ixx, iyy, izz, 1)
          ny = zf_centroid_patch ( ixx, iyy, izz, 2)
          nz = zf_centroid_patch ( ixx, iyy, izz, 3)

          a12 = 1.0 / AMAX1 ( SQRT ( nx**2 + ny**2 + nz**2 ), smallnum)

          nx = nx * a12
          ny = ny * a12
          nz = nz * a12

! Project the velocity along the normal and scale to the size of this zone.

          divvelgr ( ixx, iyy, izz) = divvelgr ( ixx, iyy, izz)

     1    + nx * a11 * ( V ( ixx, iyy, izzp1, 2, 1)
     1                 - V ( ixx, iyy, izz, 2, 1) )

     1    + ny * a11 * ( V ( ixx, iyy, izzp1, 3, 1)
     1                 - V ( ixx, iyy, izz, 3, 1) )

     1    + nz * a11 * ( V ( ixx, iyy, izzp1, 4, 1)
     1                 - V ( ixx, iyy, izz, 4, 1) )

! ----------
! Face 6, bottom z-face with normal in the direction "zf_centroid_patch ( ixx, iyy, izzm1, :)".
! Physical zone in the bottom radial direction "( ixx, iyy, izzm1)".

! Normal is not recomputed because it is the same as above.

! Build scaled relative coordinates from zone centroid to zone centroid.
! Then scale that with the scaling length for these two zones.

          x = ( v_centroid_patch ( ixx, iyy, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izzm1, 1) ) * obl0

          y = ( v_centroid_patch ( ixx, iyy, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izzm1, 2) ) * obl0

          z = ( v_centroid_patch ( ixx, iyy, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izzm1, 3) ) * obl0

          a11 = 0.5 / AMAX1 ( SQRT ( x**2 + y**2 + z**2 ), smallnum)

! Project the velocity along the normal and scale to the size of this zone.

          divvelgr ( ixx, iyy, izz) = divvelgr ( ixx, iyy, izz)

     1    + nx * a11 * ( V ( ixx, iyy, izz, 2, 1)
     1                 - V ( ixx, iyy, izzm1, 2, 1) )

     1    + ny * a11 * ( V ( ixx, iyy, izz, 3, 1)
     1                 - V ( ixx, iyy, izzm1, 3, 1) )

     1    + nz * a11 * ( V ( ixx, iyy, izz, 4, 1)
     1                 - V ( ixx, iyy, izzm1, 4, 1) )

! ------------------------------
! Now that we have built "divvelgr", use it to build "flattengr".
! This version of the flattener is based on scaling the divergence of the valocity.

          a13 = ABS ( divvelgr ( ixx, iyy, izz) )
     1     / ( flatten_coef * workgr ( ixx, iyy, izz) )

          a13 = a13 - 1.0

          IF ( a13 .LE. 0.0) THEN
            flattengr ( ixx, iyy, izz) = 0.0
          ELSE
            flattengr ( ixx, iyy, izz) = AMIN1 ( 1.0, a13 )
          END IF

          IF ( imask_active_patch ( ixx, iyy, izz) .NE. 1)
     1    flattengr ( ixx, iyy, izz) = 0.0

! ----------------------------------------

          END DO  ! End of "DO ixx = 2 - NGhost, ntp_zones + NGhost - 1" loop.

        END DO  ! End of "DO iyy = 2 - NGhost, ntp_zones + NGhost - 1" loop.

      END DO  ! End of "DO izz = 2 - NGhost, nr_zones + NGhost - 1" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

! STEP 6.b) By this point, a basic amount of flattener is set up. So we can
! even skip the rest of the flattener routine if we wish.

! Here we loop over zones that are within two zone widths of the boundary of the
! dynamically active part of the mesh.

! If neighboring zone has pressure that is larger than current zone's pressure
! by a factor of "prs_flatten_coef" then the current zone should take on the
! larger value of the flattener.
! ------------------------------------------------------------

!$acc kernels
      workgr ( :, :, :) = flattengr ( :, :, :)
!$acc end kernels

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    ixxp1, iyyp1, izzp1,                                          &
!$OMP&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd)               &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    prs_flatten_coef, imask_active_patch,                         &
!$OMP&    istencil_patch,                                               &
!$OMP&    V, workgr, flattengr)                                         &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    prs_flatten_coef, imask_active_patch,                         &
!$acc&    istencil_patch,                                               &
!$acc&    V, workgr, flattengr)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    ixxp1, iyyp1, izzp1,                                          &
!$acc&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd)

      DO izz = 2 - NGhost, nr_zones + NGhost - 1

        DO iyy = 2 - NGhost, ntp_zones + NGhost - 1

          DO ixx = 2 - NGhost, ntp_zones + NGhost - 1

          izzm1 = izz - 1
          izzp1 = izz + 1

          iyym1 = iyy - 1
          iyyp1 = iyy + 1

          ixxm1 = ixx - 1
          ixxp1 = ixx + 1

          ixxr = istencil_patch ( ixx, iyy, 1, 1)  ! For indexing zone one Right.
          iyyr = istencil_patch ( ixx, iyy, 2, 1)

          ixxl = istencil_patch ( ixx, iyy, 1, 2)  ! For indexing zone one Left.
          iyyl = istencil_patch ( ixx, iyy, 2, 2)

          ixxu = istencil_patch ( ixx, iyy, 1, 3)  ! For indexing zone one Up.
          iyyu = istencil_patch ( ixx, iyy, 2, 3)

          ixxd = istencil_patch ( ixx, iyy, 1, 4)  ! For indexing zone one Down.
          iyyd = istencil_patch ( ixx, iyy, 2, 4)

! ----------------------------------------

! Physical zone to the Right "( ixxr, iyyr, izz)".

          IF ( V ( ixxr, iyyr, izz, 5, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 5, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxr, iyyr, izz) )

! Physical zone to the Left "( ixxl, iyyl, izz)".

          IF ( V ( ixxl, iyyl, izz, 5, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 5, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxl, iyyl, izz) )

! Physical zone to the Up "( ixxu, iyyu, izz)".

          IF ( V ( ixxu, iyyu, izz, 5, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 5, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxu, iyyu, izz) )

! Physical zone to the Down "( ixxd, iyyd, izz)".

          IF ( V ( ixxd, iyyd, izz, 5, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 5, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxd, iyyd, izz) )

! Physical zone in the top radial direction "( ixx, iyy, izzp1)".

          IF ( V ( ixx, iyy, izzp1, 5, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 5, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixx, iyy, izzp1) )

! Physical zone in the bottom radial direction "( ixx, iyy, izzm1)".

          IF ( V ( ixx, iyy, izzm1, 5, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 5, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixx, iyy, izzm1) )

! ----------------------------------------
! For CGL, the above equation set will only have flattened based on parallel pressure.
! Now do the same based on perpendicular pressure, which is the 9th component of the
! primitive variables.

#if ( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3)

! Physical zone to the Right "( ixxr, iyyr, izz)".

          IF ( V ( ixxr, iyyr, izz, 9, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 9, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxr, iyyr, izz) )

! Physical zone to the Left "( ixxl, iyyl, izz)".

          IF ( V ( ixxl, iyyl, izz, 9, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 9, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxl, iyyl, izz) )

! Physical zone to the Up "( ixxu, iyyu, izz)".

          IF ( V ( ixxu, iyyu, izz, 9, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 9, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxu, iyyu, izz) )

! Physical zone to the Down "( ixxd, iyyd, izz)".

          IF ( V ( ixxd, iyyd, izz, 9, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 9, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxd, iyyd, izz) )

! Physical zone in the top radial direction "( ixx, iyy, izzp1)".

          IF ( V ( ixx, iyy, izzp1, 9, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 9, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixx, iyy, izzp1) )

! Physical zone in the bottom radial direction "( ixx, iyy, izzm1)".

          IF ( V ( ixx, iyy, izzm1, 9, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 9, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixx, iyy, izzm1) )

#endif
! End of "IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3".

! ----------------------------------------
! For CGL with electron pressure, we can also flatten based on electron pressure,
! which is the 10th component of the primitive variables.

#if ( IPDE_SYSTEM == 3)

! Physical zone to the Right "( ixxr, iyyr, izz)".

          IF ( V ( ixxr, iyyr, izz, 10, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 10, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxr, iyyr, izz) )

! Physical zone to the Left "( ixxl, iyyl, izz)".

          IF ( V ( ixxl, iyyl, izz, 10, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 10, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxl, iyyl, izz) )

! Physical zone to the Up "( ixxu, iyyu, izz)".

          IF ( V ( ixxu, iyyu, izz, 10, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 10, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxu, iyyu, izz) )

! Physical zone to the Down "( ixxd, iyyd, izz)".

          IF ( V ( ixxd, iyyd, izz, 10, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 10, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxd, iyyd, izz) )

! Physical zone in the top radial direction "( ixx, iyy, izzp1)".

          IF ( V ( ixx, iyy, izzp1, 10, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 10, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixx, iyy, izzp1) )

! Physical zone in the bottom radial direction "( ixx, iyy, izzm1)".

          IF ( V ( ixx, iyy, izzm1, 10, 1) .GT.
     1         prs_flatten_coef * V ( ixx, iyy, izz, 10, 1) )
     1 flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixx, iyy, izzm1) )

#endif
! End of "IPDE_SYSTEM == 3".

! ----------------------------------------

          END DO  ! End of "DO ixx = 2 - NGhost, ntp_zones + NGhost - 1" loop.

        END DO  ! End of "DO iyy = 2 - NGhost, ntp_zones + NGhost - 1" loop.

      END DO  ! End of "DO izz = 2 - NGhost, nr_zones + NGhost - 1" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

! STEP 6.c) If we wish, we can upwind the flattener. We do that here.

! Here we loop over zones that are within two zone widths of the boundary of the
! dynamically active part of the mesh.

! ------------------------------------------------------------

!$acc kernels
      workgr ( :, :, :) = flattengr ( :, :, :)
!$acc end kernels

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    ixxp1, iyyp1, izzp1,                                          &
!$OMP&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd)               &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    prs_flatten_coef, imask_active_patch,                         &
!$OMP&    istencil_patch,                                               &
!$OMP&    V, workgr, flattengr)                                         &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    prs_flatten_coef, imask_active_patch,                         &
!$acc&    istencil_patch,                                               &
!$acc&    V, workgr, flattengr)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    ixxp1, iyyp1, izzp1,                                          &
!$acc&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd)

      DO izz = 2 - NGhost, nr_zones + NGhost - 1

        DO iyy = 2 - NGhost, ntp_zones + NGhost - 1

          DO ixx = 2 - NGhost, ntp_zones + NGhost - 1

          izzm1 = izz - 1
          izzp1 = izz + 1

          iyym1 = iyy - 1
          iyyp1 = iyy + 1

          ixxm1 = ixx - 1
          ixxp1 = ixx + 1

          ixxr = istencil_patch ( ixx, iyy, 1, 1)  ! For indexing zone one Right.
          iyyr = istencil_patch ( ixx, iyy, 2, 1)

          ixxl = istencil_patch ( ixx, iyy, 1, 2)  ! For indexing zone one Left.
          iyyl = istencil_patch ( ixx, iyy, 2, 2)

          ixxu = istencil_patch ( ixx, iyy, 1, 3)  ! For indexing zone one Up.
          iyyu = istencil_patch ( ixx, iyy, 2, 3)

          ixxd = istencil_patch ( ixx, iyy, 1, 4)  ! For indexing zone one Down.
          iyyd = istencil_patch ( ixx, iyy, 2, 4)

! ----------------------------------------

! Physical zone to the Right "( ixxr, iyyr, izz)".

       flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxr, iyyr, izz) )

! Physical zone to the Left "( ixxl, iyyl, izz)".

       flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxl, iyyl, izz) )

! Physical zone to the Up "( ixxu, iyyu, izz)".

       flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxu, iyyu, izz) )

! Physical zone to the Down "( ixxd, iyyd, izz)".

       flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixxd, iyyd, izz) )

! Physical zone in the top radial direction "( ixx, iyy, izzp1)".

       flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixx, iyy, izzp1) )

! Physical zone in the bottom radial direction "( ixx, iyy, izzm1)".

       flattengr ( ixx, iyy, izz) = AMAX1 ( flattengr ( ixx, iyy, izz),
     1                                      workgr ( ixx, iyy, izzm1) )

! ----------------------------------------

          END DO  ! End of "DO ixx = 2 - NGhost, ntp_zones + NGhost - 1" loop.

        END DO  ! End of "DO iyy = 2 - NGhost, ntp_zones + NGhost - 1" loop.

      END DO  ! End of "DO izz = 2 - NGhost, nr_zones + NGhost - 1" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 7) We limit the zone-centered variables by using WENO-AO limiting.
! The idea here is that we have nine smaller piecewise linear stencils that
! will eventually be non-linearly hybridized with a piecewise quadratic
! central stencil. This makes the process very inexpensive.

! ----------------------------------------------------------------------
! STEP 7.a) WENO limiting of primitive variables "V" and "entropygr".
! This is a second order, but truly multid, WENO
! reconstruction. We have been able to simplify the process to the point where
! we only have to invert eight small 3X3 matrices for the eight one-sided
! stencils. The ninth, centered stencil is obtained by arithmetic averaging
! the linear moments from the previous eight stencils and is given a much higher
! weight. We only do the reconstruction for zones that are two zones out from
! the physical zones.

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1, icomp,  &
!$OMP&    ixxp1, iyyp1, izzp1, obl0, i, j, nn,                          &
!$OMP&    x, y, z,                                                      &
!$OMP&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,               &
!$OMP&    ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,       &
!$OMP&    a11, a12, a13, a21, a22, a23, a31, a32, a33, obDet,           &
!$OMP&    R38, X39, AInv338, IS9, Wt9)                                  &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,&
!$OMP&    pow, central_weight, tinynum, smallnum,                       &
!$OMP&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$OMP&    vert_coords_patch, Coef_Basis_Moments,                        &
!$OMP&    istencil_patch, flattengr, V, entropygr)                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    pow, central_weight, tinynum, smallnum,                       &
!$acc&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$acc&    vert_coords_patch, Coef_Basis_Moments,                        &
!$acc&    istencil_patch, flattengr, V, entropygr)

!$acc loop gang collapse(2) independent                                 &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1, icomp,          &
!$acc&    ixxp1, iyyp1, izzp1, obl0, i, j, nn,                          &
!$acc&    x, y, z,                                                      &
!$acc&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,               &
!$acc&    ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,       &
!$acc&    a11, a12, a13, a21, a22, a23, a31, a32, a33, obDet,           &
!$acc&    AInv338, R38, X39, IS9, Wt9)

! ----------------------------------------------------------------------

      DO izz = - 2, nr_zones + 3

        DO iyy = - 2, ntp_zones + 3

          DO ixx = - 2, ntp_zones + 3

          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN

! ----------------------------------------------------------------------

      izzm1 = izz - 1
      izzp1 = izz + 1

      iyym1 = iyy - 1
      iyyp1 = iyy + 1

      ixxm1 = ixx - 1
      ixxp1 = ixx + 1

      ixxr = istencil_patch ( ixx, iyy, 1, 1)  ! For indexing zone one Right.
      iyyr = istencil_patch ( ixx, iyy, 2, 1)

      ixxl = istencil_patch ( ixx, iyy, 1, 2)  ! For indexing zone one Left.
      iyyl = istencil_patch ( ixx, iyy, 2, 2)

      ixxu = istencil_patch ( ixx, iyy, 1, 3)  ! For indexing zone one Up.
      iyyu = istencil_patch ( ixx, iyy, 2, 3)

      ixxd = istencil_patch ( ixx, iyy, 1, 4)  ! For indexing zone one Down.
      iyyd = istencil_patch ( ixx, iyy, 2, 4)

      ixxru = istencil_patch ( ixx, iyy, 1, 9)  ! For indexing zone that is Right-Up.
      iyyru = istencil_patch ( ixx, iyy, 2, 9)

      ixxlu = istencil_patch ( ixx, iyy, 1, 10)  ! For indexing zone that is Left-Up.
      iyylu = istencil_patch ( ixx, iyy, 2, 10)

      ixxld = istencil_patch ( ixx, iyy, 1, 11)  ! For indexing zone that is Left-Down.
      iyyld = istencil_patch ( ixx, iyy, 2, 11)

      ixxrd = istencil_patch ( ixx, iyy, 1, 12)  ! For indexing zone that is Right-Down.
      iyyrd = istencil_patch ( ixx, iyy, 2, 12)


! ----------------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

! --------------------------------------------------
! We start on the expensive part associated with building "AInv338".
! --------------------------------------------------

! ------------------------------
! Stencil 1 in direction "( Rt, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

          a11 = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a12 = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a13 = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a21 = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a22 = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a23 = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 1 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 1  ! Stencil 1

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 2 in direction "( Rt, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

!          a11 = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a12 = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a13 = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


!          a21 = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a22 = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a23 = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzm1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzm1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzm1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 2 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 2  ! Stencil 2

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 3 in direction "( Lf, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

          a11 = ( v_centroid_patch ( ixxl, iyyl, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a12 = ( v_centroid_patch ( ixxl, iyyl, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a13 = ( v_centroid_patch ( ixxl, iyyl, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a21 = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a22 = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a23 = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 3 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 3  ! Stencil 3

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil in direction 4 "( Lf, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

!          a11 = ( v_centroid_patch ( ixxl, iyyl, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a12 = ( v_centroid_patch ( ixxl, iyyl, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a13 = ( v_centroid_patch ( ixxl, iyyl, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


!          a21 = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a22 = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a23 = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzm1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzm1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzm1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 4 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 4  ! Stencil 4

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 5 in direction "( Rt, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

          a11 = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a12 = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a13 = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a21 = ( v_centroid_patch ( ixxd, iyyd, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a22 = ( v_centroid_patch ( ixxd, iyyd, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a23 = ( v_centroid_patch ( ixxd, iyyd, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 5 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 5  ! Stencil 5

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 6 in direction "( Rt, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

!          a11 = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a12 = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a13 = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


!          a21 = ( v_centroid_patch ( ixxd, iyyd, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a22 = ( v_centroid_patch ( ixxd, iyyd, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a23 = ( v_centroid_patch ( ixxd, iyyd, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzm1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzm1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzm1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 6 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 6  ! Stencil 6

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 7 in direction "( Lf, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

          a11 = ( v_centroid_patch ( ixxl, iyyl, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a12 = ( v_centroid_patch ( ixxl, iyyl, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a13 = ( v_centroid_patch ( ixxl, iyyl, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a21 = ( v_centroid_patch ( ixxd, iyyd, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a22 = ( v_centroid_patch ( ixxd, iyyd, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a23 = ( v_centroid_patch ( ixxd, iyyd, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 7 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 7  ! Stencil 7

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 8 in direction "( Lf, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

!          a11 = ( v_centroid_patch ( ixxl, iyyl, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a12 = ( v_centroid_patch ( ixxl, iyyl, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a13 = ( v_centroid_patch ( ixxl, iyyl, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


!          a21 = ( v_centroid_patch ( ixxd, iyyd, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a22 = ( v_centroid_patch ( ixxd, iyyd, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a23 = ( v_centroid_patch ( ixxd, iyyd, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzm1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzm1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzm1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 8 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 8  ! Stencil 8

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! --------------------------------------------------
! By this point, the expensive part associated with building "AInv338" is done.
! --------------------------------------------------

      DO icomp = 1, n_cc_components  ! Loop over the components of "V".

! --------------------------------------------------

! ------------------------------
! Stencil 1 in direction "( Rt, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 1

          R38 ( 1, j) = V ( ixxr, iyyr, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 2, j) = V ( ixxu, iyyu, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 3, j) = V ( ixx, iyy, izzp1, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 2 in direction "( Rt, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 2

          R38 ( 1, j) = V ( ixxr, iyyr, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 2, j) = V ( ixxu, iyyu, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 3, j) = V ( ixx, iyy, izzm1, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 3 in direction "( Lf, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 3

          R38 ( 1, j) = V ( ixxl, iyyl, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 2, j) = V ( ixxu, iyyu, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 3, j) = V ( ixx, iyy, izzp1, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil in direction 4 "( Lf, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 4

          R38 ( 1, j) = V ( ixxl, iyyl, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 2, j) = V ( ixxu, iyyu, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 3, j) = V ( ixx, iyy, izzm1, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 5 in direction "( Rt, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 5

          R38 ( 1, j) = V ( ixxr, iyyr, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 2, j) = V ( ixxd, iyyd, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 3, j) = V ( ixx, iyy, izzp1, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 6 in direction "( Rt, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 6

          R38 ( 1, j) = V ( ixxr, iyyr, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 2, j) = V ( ixxd, iyyd, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 3, j) = V ( ixx, iyy, izzm1, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 7 in direction "( Lf, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 7

          R38 ( 1, j) = V ( ixxl, iyyl, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 2, j) = V ( ixxd, iyyd, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 3, j) = V ( ixx, iyy, izzp1, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 8 in direction "( Lf, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 8

          R38 ( 1, j) = V ( ixxl, iyyl, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 2, j) = V ( ixxd, iyyd, izz, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

          R38 ( 3, j) = V ( ixx, iyy, izzm1, icomp, 1)
     1                     - V ( ixx, iyy, izz, icomp, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! By this point, all three gradients from the first 8 stencils have been
! evaluated and stored in "X39 ( 1:3, 1:8)". Now obtain the
! centered, ninth stencil's gradients by arithmetic averaging.
! ------------------------------

          X39 ( :, 9) = 0.125 * ( X39 ( :, 1)
     1  + X39 ( :, 2) + X39 ( :, 3) + X39 ( :, 4)
     1  + X39 ( :, 5) + X39 ( :, 6) + X39 ( :, 7)
     1  + X39 ( :, 8) )

! ------------------------------
! Build all 9 smoothness indicators for all 9 stencils.
! ------------------------------

          IS9 ( 1: 9) = X39 ( 1, 1: 9)**2
     1  + X39 ( 2, 1: 9)**2 + X39 ( 3, 1: 9)**2

! ------------------------------
! Initialize 9 linear weights and use them to obtain non-linear weights
! for all 9 stencils. Note that here the ninth stencil is the central one.
! ------------------------------

          Wt9 ( 1: 8) = 1.0
          Wt9 ( 9) = central_weight

          Wt9 ( 1: 9) = Wt9 ( 1: 9) / ( IS9 ( 1: 9)
     1                                          + smallnum)**pow

          a11 = 1.0 / ( Wt9 ( 1) + Wt9 ( 2)
     1      + Wt9 ( 3) + Wt9 ( 4) + Wt9 ( 5)
     1      + Wt9 ( 6) + Wt9 ( 7) + Wt9 ( 8)
     1      + Wt9 ( 9) )

          Wt9 ( 1: 9) = Wt9 ( 1: 9) * a11

! ------------------------------
! Make all the three non-linearly hybridized gradients of
! "W ( ixx, iyy, izz, 2: 4)" for this entropy variable.
! ------------------------------

          V ( ixx, iyy, izz, icomp, 2) =
     1        Wt9 ( 1) * X39 ( 1, 1)
     1      + Wt9 ( 2) * X39 ( 1, 2)
     1      + Wt9 ( 3) * X39 ( 1, 3)
     1      + Wt9 ( 4) * X39 ( 1, 4)
     1      + Wt9 ( 5) * X39 ( 1, 5)
     1      + Wt9 ( 6) * X39 ( 1, 6)
     1      + Wt9 ( 7) * X39 ( 1, 7)
     1      + Wt9 ( 8) * X39 ( 1, 8)
     1      + Wt9 ( 9) * X39 ( 1, 9)

          V ( ixx, iyy, izz, icomp, 3) =
     1        Wt9 ( 1) * X39 ( 2, 1)
     1      + Wt9 ( 2) * X39 ( 2, 2)
     1      + Wt9 ( 3) * X39 ( 2, 3)
     1      + Wt9 ( 4) * X39 ( 2, 4)
     1      + Wt9 ( 5) * X39 ( 2, 5)
     1      + Wt9 ( 6) * X39 ( 2, 6)
     1      + Wt9 ( 7) * X39 ( 2, 7)
     1      + Wt9 ( 8) * X39 ( 2, 8)
     1      + Wt9 ( 9) * X39 ( 2, 9)

          V ( ixx, iyy, izz, icomp, 4) =
     1        Wt9 ( 1) * X39 ( 3, 1)
     1      + Wt9 ( 2) * X39 ( 3, 2)
     1      + Wt9 ( 3) * X39 ( 3, 3)
     1      + Wt9 ( 4) * X39 ( 3, 4)
     1      + Wt9 ( 5) * X39 ( 3, 5)
     1      + Wt9 ( 6) * X39 ( 3, 6)
     1      + Wt9 ( 7) * X39 ( 3, 7)
     1      + Wt9 ( 8) * X39 ( 3, 8)
     1      + Wt9 ( 9) * X39 ( 3, 9)

! ------------------------------
! "flattengr" will be applied later; so it is commented out here.
! ------------------------------

!          V ( ixx, iyy, izz, icomp, 2) =
!     1    V ( ixx, iyy, izz, icomp, 2)
!     1  * ( 1.0 - flattengr ( ixx, iyy, izz) )

!          V ( ixx, iyy, izz, icomp, 3) =
!     1    V ( ixx, iyy, izz, icomp, 3)
!     1  * ( 1.0 - flattengr ( ixx, iyy, izz) )

!          V ( ixx, iyy, izz, icomp, 4) =
!     1    V ( ixx, iyy, izz, icomp, 4)
!     1  * ( 1.0 - flattengr ( ixx, iyy, izz) )

! --------------------------------------------------

          END DO    ! End of "DO icomp = 1, n_cc_components" loop.

! --------------------------------------------------
! Start limiting the entropy variable "entropygr".
#if (PRESSURE_POSITIVITY > 0)
! --------------------------------------------------

! ------------------------------
! Stencil 1 in direction "( Rt, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 1

          R38 ( 1, j) = entropygr ( ixxr, iyyr, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 2, j) = entropygr ( ixxu, iyyu, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 3, j) = entropygr ( ixx, iyy, izzp1, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 2 in direction "( Rt, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 2

          R38 ( 1, j) = entropygr ( ixxr, iyyr, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 2, j) = entropygr ( ixxu, iyyu, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 3, j) = entropygr ( ixx, iyy, izzm1, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 3 in direction "( Lf, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 3

          R38 ( 1, j) = entropygr ( ixxl, iyyl, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 2, j) = entropygr ( ixxu, iyyu, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 3, j) = entropygr ( ixx, iyy, izzp1, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil in direction 4 "( Lf, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 4

          R38 ( 1, j) = entropygr ( ixxl, iyyl, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 2, j) = entropygr ( ixxu, iyyu, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 3, j) = entropygr ( ixx, iyy, izzm1, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 5 in direction "( Rt, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 5

          R38 ( 1, j) = entropygr ( ixxr, iyyr, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 2, j) = entropygr ( ixxd, iyyd, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 3, j) = entropygr ( ixx, iyy, izzp1, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 6 in direction "( Rt, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 6

          R38 ( 1, j) = entropygr ( ixxr, iyyr, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 2, j) = entropygr ( ixxd, iyyd, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 3, j) = entropygr ( ixx, iyy, izzm1, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 7 in direction "( Lf, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 7

          R38 ( 1, j) = entropygr ( ixxl, iyyl, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 2, j) = entropygr ( ixxd, iyyd, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 3, j) = entropygr ( ixx, iyy, izzp1, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 8 in direction "( Lf, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 8

          R38 ( 1, j) = entropygr ( ixxl, iyyl, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 2, j) = entropygr ( ixxd, iyyd, izz, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

          R38 ( 3, j) = entropygr ( ixx, iyy, izzm1, 1)
     1                     - entropygr ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! By this point, all three gradients from the first 8 stencils have been
! evaluated and stored in "X39 ( 1:3, 1:8)". Now obtain the
! centered, ninth stencil's gradients by arithmetic averaging.
! ------------------------------

          X39 ( :, 9) = 0.125 * ( X39 ( :, 1)
     1  + X39 ( :, 2) + X39 ( :, 3) + X39 ( :, 4)
     1  + X39 ( :, 5) + X39 ( :, 6) + X39 ( :, 7)
     1  + X39 ( :, 8) )

! ------------------------------
! Build all 9 smoothness indicators for all 9 stencils.
! ------------------------------

          IS9 ( 1: 9) = X39 ( 1, 1: 9)**2
     1  + X39 ( 2, 1: 9)**2 + X39 ( 3, 1: 9)**2

! ------------------------------
! Initialize 9 linear weights and use them to obtain non-linear weights
! for all 9 stencils. Note that here the ninth stencil is the central one.
! ------------------------------

          Wt9 ( 1: 8) = 1.0
          Wt9 ( 9) = central_weight

          Wt9 ( 1: 9) = Wt9 ( 1: 9) / ( IS9 ( 1: 9)
     1                                          + smallnum)**pow

          a11 = 1.0 / ( Wt9 ( 1) + Wt9 ( 2)
     1      + Wt9 ( 3) + Wt9 ( 4) + Wt9 ( 5)
     1      + Wt9 ( 6) + Wt9 ( 7) + Wt9 ( 8)
     1      + Wt9 ( 9) )

          Wt9 ( 1: 9) = Wt9 ( 1: 9) * a11

! ------------------------------
! Make all the three non-linearly hybridized gradients of
! "entropygr ( ixx, iyy, izz, 2: 4)" for this entropy variable.
! ------------------------------

          entropygr ( ixx, iyy, izz, 2) =
     1        Wt9 ( 1) * X39 ( 1, 1)
     1      + Wt9 ( 2) * X39 ( 1, 2)
     1      + Wt9 ( 3) * X39 ( 1, 3)
     1      + Wt9 ( 4) * X39 ( 1, 4)
     1      + Wt9 ( 5) * X39 ( 1, 5)
     1      + Wt9 ( 6) * X39 ( 1, 6)
     1      + Wt9 ( 7) * X39 ( 1, 7)
     1      + Wt9 ( 8) * X39 ( 1, 8)
     1      + Wt9 ( 9) * X39 ( 1, 9)

          entropygr ( ixx, iyy, izz, 3) =
     1        Wt9 ( 1) * X39 ( 2, 1)
     1      + Wt9 ( 2) * X39 ( 2, 2)
     1      + Wt9 ( 3) * X39 ( 2, 3)
     1      + Wt9 ( 4) * X39 ( 2, 4)
     1      + Wt9 ( 5) * X39 ( 2, 5)
     1      + Wt9 ( 6) * X39 ( 2, 6)
     1      + Wt9 ( 7) * X39 ( 2, 7)
     1      + Wt9 ( 8) * X39 ( 2, 8)
     1      + Wt9 ( 9) * X39 ( 2, 9)

          entropygr ( ixx, iyy, izz, 4) =
     1        Wt9 ( 1) * X39 ( 3, 1)
     1      + Wt9 ( 2) * X39 ( 3, 2)
     1      + Wt9 ( 3) * X39 ( 3, 3)
     1      + Wt9 ( 4) * X39 ( 3, 4)
     1      + Wt9 ( 5) * X39 ( 3, 5)
     1      + Wt9 ( 6) * X39 ( 3, 6)
     1      + Wt9 ( 7) * X39 ( 3, 7)
     1      + Wt9 ( 8) * X39 ( 3, 8)
     1      + Wt9 ( 9) * X39 ( 3, 9)

! ------------------------------
! Reduce the gradients if flattener is active.
! ("flattengr" is 0.0 near smooth flow and 1.0 in regions with v. strong
!   shocks and rarefactions.)
! Since "entropygr" only has three moments, it is ok to apply "flattengr" here.
! ------------------------------

          entropygr ( ixx, iyy, izz, 2) =
     1    entropygr ( ixx, iyy, izz, 2)
     1  * ( 1.0 - flattengr ( ixx, iyy, izz) )

          entropygr ( ixx, iyy, izz, 3) =
     1    entropygr ( ixx, iyy, izz, 3)
     1  * ( 1.0 - flattengr ( ixx, iyy, izz) )

          entropygr ( ixx, iyy, izz, 4) =
     1    entropygr ( ixx, iyy, izz, 4)
     1  * ( 1.0 - flattengr ( ixx, iyy, izz) )

! --------------------------------------------------
#endif
! End of "PRESSURE_POSITIVITY > 0".
! End of limiting the entropy variable "entropygr".
! ----------------------------------------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = - 2, ntp_zones + 3" loop.

        END DO  ! End of "DO iyy = - 2, ntp_zones + 3" loop.

      END DO  ! End of "DO izz = - 2, nr_zones + 3" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 7.b) In STEP 7.a, we have already obtained the piecewise linear
! variation of "V". Here we non-linearly hybridize it with just one
! large central stencil that is piecewise quadratic. This saves a lot of
! computational cost with only minimal degradation in accuracy.

! ----------------------------------------------------------------------

      iorder = 3  ! At third order, all stencils are third order.

      n_cc_10 = 10  ! These are the total number of modes that we will use at 3rd order.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1, icomp,  &
!$OMP&    ixxp1, iyyp1, izzp1, obl0, i, j, nn,                          &
!$OMP&    ixxm2, ixxp2, iyym2, iyyp2, izzm2, izzp2,                     &
!$OMP&    x, y, z, ell_0, ell_j, Coef_0, Coef_j, indx_A9, Amat99,       &
!$OMP&    Amat99_Inv, vv_A9, Amat_Row, Amat189, Amat189T,               &
!$OMP&    X9, Y9, R18, ISlo, IShi, Wtlo, Wthi, aa,                      &
!$OMP&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,               &
!$OMP&    ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,       &
!$OMP&    ixxrr, iyyrr, ixxll, iyyll, ixxuu, iyyuu, ixxdd, iyydd,       &
!$OMP&    Phi9, betalo, gammahi, gammalo, gammaratio, obgammahi, tau)   &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,&
!$OMP&    pow, central_weight, tinynum, smallnum, n_cc_10, n_cc_modes,  &
!$OMP&    ob12, ob180, ob144, iorder,                                   &
!$OMP&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$OMP&    vert_coords_patch, Coef_Basis_Moments,                        &
!$OMP&    istencil_patch, flattengr, V, entropygr)                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    pow, central_weight, tinynum, smallnum, n_cc_10, n_cc_modes,  &
!$acc&    ob12, ob180, ob144, iorder,                                   &
!$acc&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$acc&    vert_coords_patch, Coef_Basis_Moments,                        &
!$acc&    istencil_patch, flattengr, V, entropygr)

!$acc loop gang collapse(2) independent                                 &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1, icomp,          &
!$acc&    ixxp1, iyyp1, izzp1, obl0, i, j, nn,                          &
!$acc&    ixxm2, ixxp2, iyym2, iyyp2, izzm2, izzp2,                     &
!$acc&    x, y, z, ell_0, ell_j, Coef_0, Coef_j, indx_A9, Amat99,       &
!$acc&    Amat99_Inv, vv_A9, Amat_Row, Amat189, Amat189T,               &
!$acc&    X9, Y9, R18, ISlo, IShi, Wtlo, Wthi, aa,                      &
!$acc&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,               &
!$acc&    ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,       &
!$acc&    ixxrr, iyyrr, ixxll, iyyll, ixxuu, iyyuu, ixxdd, iyydd,       &
!$acc&    Phi9, betalo, gammahi, gammalo, gammaratio, obgammahi, tau)   &

! ----------------------------------------------------------------------

!      DO izz = - 2, nr_zones + 3

!        DO iyy = - 2, ntp_zones + 3

!          DO ixx = - 2, ntp_zones + 3

      DO izz = - 1, nr_zones + 2

        DO iyy = - 1, ntp_zones + 2

          DO ixx = - 1, ntp_zones + 2

          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN

! ----------------------------------------------------------------------

      izzm1 = izz - 1
      izzp1 = izz + 1
      izzm2 = izz - 2
      izzp2 = izz + 2

      iyym1 = iyy - 1
      iyyp1 = iyy + 1
      iyym2 = iyy - 2
      iyyp2 = iyy + 2

      ixxm1 = ixx - 1
      ixxp1 = ixx + 1
      ixxm2 = ixx - 2
      ixxp2 = ixx + 2

      ixxr = istencil_patch ( ixx, iyy, 1, 1)   ! For indexing zone one Right.
      iyyr = istencil_patch ( ixx, iyy, 2, 1)

      ixxl = istencil_patch ( ixx, iyy, 1, 2)   ! For indexing zone one Left.
      iyyl = istencil_patch ( ixx, iyy, 2, 2)

      ixxu = istencil_patch ( ixx, iyy, 1, 3)   ! For indexing zone one Up.
      iyyu = istencil_patch ( ixx, iyy, 2, 3)

      ixxd = istencil_patch ( ixx, iyy, 1, 4)   ! For indexing zone one Down.
      iyyd = istencil_patch ( ixx, iyy, 2, 4)

      ixxrr = istencil_patch ( ixx, iyy, 1, 5)  ! For indexing zone two Right.
      iyyrr = istencil_patch ( ixx, iyy, 2, 5)

      ixxll = istencil_patch ( ixx, iyy, 1, 6)  ! For indexing zone two Left.
      iyyll = istencil_patch ( ixx, iyy, 2, 6)

      ixxuu = istencil_patch ( ixx, iyy, 1, 7)  ! For indexing zone two Up.
      iyyuu = istencil_patch ( ixx, iyy, 2, 7)

      ixxdd = istencil_patch ( ixx, iyy, 1, 8)  ! For indexing zone two Down.
      iyydd = istencil_patch ( ixx, iyy, 2, 8)

      ixxru = istencil_patch ( ixx, iyy, 1, 9)  ! For indexing zone that is Right-Up.
      iyyru = istencil_patch ( ixx, iyy, 2, 9)

      ixxlu = istencil_patch ( ixx, iyy, 1, 10)  ! For indexing zone that is Left-Up.
      iyylu = istencil_patch ( ixx, iyy, 2, 10)

      ixxld = istencil_patch ( ixx, iyy, 1, 11)  ! For indexing zone that is Left-Down.
      iyyld = istencil_patch ( ixx, iyy, 2, 11)

      ixxrd = istencil_patch ( ixx, iyy, 1, 12)  ! For indexing zone that is Right-Down.
      iyyrd = istencil_patch ( ixx, iyy, 2, 12)


! ----------------------------------------

          ell_0 = AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)

          obl0 = 1.0 / ell_0  ! Scaling for this zone.

          Coef_0 ( 1: n_cc_10) = Coef_Basis_Moments ( ixx, iyy, izz,
     1                              1: n_cc_10)

! --------------------------------------------------
! We start on the expensive part associated with building "Amat99_Inv".
! --------------------------------------------------

! ------------------------------
! Stencil 9 is centered in all 3 directions and it uses the eighteen zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixxd, iyyd, izz),
!    ( ixx, iyy, izzp1), ( ixx, iyy, izzm1),
!    ( ixxru, iyyru, izz), ( ixxlu, iyylu, izz), ( ixxld, iyyld, izz), ( ixxrd, iyyrd, izz),
!    ( ixxr, iyyr, izzp1), ( ixxl, iyyl, izzp1), ( ixxr, iyyr, izzm1), ( ixxl, iyyl, izzm1),
!    ( ixxu, iyyu, izzp1), ( ixxd, iyyd, izzp1), ( ixxu, iyyu, izzm1), ( ixxd, iyyd, izzm1) }".
! ------------------------------

! Build 9 sets of scaled relative coordinates from 9 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 9X9 matrix for Stencil 9.

! ----------
! Stencil 9, Row 1 from "( ixxr, iyyr, izz)".

          x = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxr, iyyr, izz), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxr, iyyr, izz, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 1, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 2 from "( ixxl, iyyl, izz)".

          x = ( v_centroid_patch ( ixxl, iyyl, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxl, iyyl, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxl, iyyl, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxl, iyyl, izz), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxl, iyyl, izz, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 2, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 3 from "( ixxu, iyyu, izz)".

          x = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxu, iyyu, izz), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxu, iyyu, izz, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 3, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 4 from "( ixxd, iyyd, izz)".

          x = ( v_centroid_patch ( ixxd, iyyd, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxd, iyyd, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxd, iyyd, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxd, iyyd, izz), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxd, iyyd, izz, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 4, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 5 from "( ixx, iyy, izzp1)".

          x = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixx, iyy, izzp1), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixx, iyy, izzp1, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 5, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 6 from "( ixx, iyy, izzm1)".

          x = ( v_centroid_patch ( ixx, iyy, izzm1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixx, iyy, izzm1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixx, iyy, izzm1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixx, iyy, izzm1), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixx, iyy, izzm1, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 6, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 7 from "( ixxru, iyyru, izz)".

          x = ( v_centroid_patch ( ixxru, iyyru, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxru, iyyru, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxru, iyyru, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxru, iyyru, izz), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxru, iyyru, izz, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 7, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 8 from "( ixxlu, iyylu, izz)".

          x = ( v_centroid_patch ( ixxlu, iyylu, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxlu, iyylu, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxlu, iyylu, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxlu, iyylu, izz), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxlu, iyylu, izz, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 8, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 9 from "( ixxld, iyyld, izz)".

          x = ( v_centroid_patch ( ixxld, iyyld, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxld, iyyld, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxld, iyyld, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxld, iyyld, izz), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxld, iyyld, izz, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 9, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 10 from "( ixxrd, iyyrd, izz)".

          x = ( v_centroid_patch ( ixxrd, iyyrd, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxrd, iyyrd, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxrd, iyyrd, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxrd, iyyrd, izz), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxrd, iyyrd, izz, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 10, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 11 from "( ixxr, iyyr, izzp1)".

          x = ( v_centroid_patch ( ixxr, iyyr, izzp1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxr, iyyr, izzp1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxr, iyyr, izzp1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxr, iyyr, izzp1), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxr, iyyr, izzp1, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 11, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 12 from "( ixxl, iyyl, izzp1)".

          x = ( v_centroid_patch ( ixxl, iyyl, izzp1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxl, iyyl, izzp1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxl, iyyl, izzp1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxl, iyyl, izzp1), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxl, iyyl, izzp1, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 12, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 13 from "( ixxr, iyyr, izzm1)".

          x = ( v_centroid_patch ( ixxr, iyyr, izzm1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxr, iyyr, izzm1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxr, iyyr, izzm1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxr, iyyr, izzm1), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxr, iyyr, izzm1, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 13, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 14 from "( ixxl, iyyl, izzm1)".

          x = ( v_centroid_patch ( ixxl, iyyl, izzm1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxl, iyyl, izzm1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxl, iyyl, izzm1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxl, iyyl, izzm1), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxl, iyyl, izzm1, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 14, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 15 from "( ixxu, iyyu, izzp1)".

          x = ( v_centroid_patch ( ixxu, iyyu, izzp1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxu, iyyu, izzp1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxu, iyyu, izzp1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxu, iyyu, izzp1), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxu, iyyu, izzp1, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 15, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 16 from "( ixxd, iyyd, izzp1)".

          x = ( v_centroid_patch ( ixxd, iyyd, izzp1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxd, iyyd, izzp1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxd, iyyd, izzp1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxd, iyyd, izzp1), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxd, iyyd, izzp1, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 16, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 17 from "( ixxu, iyyu, izzm1)".

          x = ( v_centroid_patch ( ixxu, iyyu, izzm1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxu, iyyu, izzm1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxu, iyyu, izzm1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxu, iyyu, izzm1), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxu, iyyu, izzm1, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 17, :) = Amat_Row ( :)

! ----------
! Stencil 9, Row 18 from "( ixxd, iyyd, izzm1)".

          x = ( v_centroid_patch ( ixxd, iyyd, izzm1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxd, iyyd, izzm1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxd, iyyd, izzm1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          ell_j = AMAX1 ( v_scale_length ( ixxd, iyyd, izzm1), smallnum)

          Coef_j ( 1: n_cc_10) =
     1    Coef_Basis_Moments ( ixxd, iyyd, izzm1, 1: n_cc_10)

          CALL RECONSTRUCTION_COEFS ( iorder, n_cc_10, x, y, z,
     1           ell_0, ell_j, Coef_0, Coef_j, Amat_Row )

          Amat189 ( 18, :) = Amat_Row ( :)

! ----------
! Stencil 9, First build "Amat189T" from "Amat189". Use it to build "Amat99".
! Then invert "Amat99" and store it in "Amat99_Inv ( :, :)".
! Since this is a least squares process, we will also need to save the transpose
! "Amat189T" when we build the coefficients.

          Amat189T = TRANSPOSE ( Amat189)

          Amat99 = MATMUL ( Amat189T, Amat189)

          CALL MATINV ( 9, Amat99, Amat99_Inv, indx_A9, vv_A9)

! --------------------------------------------------
! By this point, the expensive part associated with building "Amat99_Inv" is done.
! This array contains the geometry information that is needed for obtaining the
! modes. It is evaluated only once for all the components.
! --------------------------------------------------

      DO icomp = 1, n_cc_components  ! Loop over the components of "V".

! --------------------------------------------------

! ------------------------------
! Stencil 9 is centered in all 3 directions and it uses the eighteen zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixxd, iyyd, izz),
!    ( ixx, iyy, izzp1), ( ixx, iyy, izzm1),
!    ( ixxru, iyyru, izz), ( ixxlu, iyylu, izz), ( ixxld, iyyld, izz), ( ixxrd, iyyrd, izz),
!    ( ixxr, iyyr, izzp1), ( ixxl, iyyl, izzp1), ( ixxr, iyyr, izzm1), ( ixxl, iyyl, izzm1),
!    ( ixxu, iyyu, izzp1), ( ixxd, iyyd, izzp1), ( ixxu, iyyu, izzm1), ( ixxd, iyyd, izzm1) }".
! ------------------------------

      R18 ( 01) = V ( ixxr, iyyr, izz, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 02) = V ( ixxl, iyyl, izz, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 03) = V ( ixxu, iyyu, izz, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 04) = V ( ixxd, iyyd, izz, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 05) = V ( ixx, iyy, izzp1, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 06) = V ( ixx, iyy, izzm1, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 07) = V ( ixxru, iyyru, izz, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 08) = V ( ixxlu, iyylu, izz, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 09) = V ( ixxld, iyyld, izz, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 10) = V ( ixxrd, iyyrd, izz, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 11) = V ( ixxr, iyyr, izzp1, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 12) = V ( ixxl, iyyl, izzp1, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 13) = V ( ixxr, iyyr, izzm1, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 14) = V ( ixxl, iyyl, izzm1, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 15) = V ( ixxu, iyyu, izzp1, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 16) = V ( ixxd, iyyd, izzp1, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 17) = V ( ixxu, iyyu, izzm1, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      R18 ( 18) = V ( ixxd, iyyd, izzm1, icomp, 1)
     1          - V ( ixx, iyy, izz, icomp, 1)

      Y9 = MATMUL ( Amat189T, R18)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X9 ( :) = MATMUL ( Amat99_Inv ( :, :), Y9)
#else
!$acc loop seq
          DO i = 1, 9
            x = 0.0
            DO nn = 1, 9
              x = x + Amat99_Inv ( i, nn) * Y9 ( nn)
            END DO
            X9 ( i) = x
          END DO
#endif

! ------------------------------
! By this point, all nine modes for the large central stencil have been
! evaluated and stored in "X9 ( 1:9)".
! ------------------------------
! Build the smoothness indicator for the large stencil and compare with the smoothness
! indicator for the small piecewise linear stencil, which has already been computed.

! The "Coef_0" variables have been scaled to a unit cube. Different moments
! in the Taylor basis will integrate to different values over such a cube.
! The skewness of the cube sphere zones has not been accounted for here.
! ------------------------------

! Coefficient by which the smoothness indicator for low polynomial is multiplied.
          betalo = flattengr ( ixx, iyy, izz)
     1           + (1.0 - flattengr ( ixx, iyy, izz)) * 2.0

! Build MRWENO linear weights.
          gammahi = flattengr ( ixx, iyy, izz) * 0.6
     1            + (1.0 - flattengr ( ixx, iyy, izz)) * 0.95  ! Linear weight for high polynomial

          gammalo = 1.0 - gammahi  ! Linear weight for low polynomial

          obgammahi = 1.0 / gammahi
          gammaratio = gammalo / gammahi


! Smoothness indicator for low polynomial.
          ISlo = betalo * ( V ( ixx, iyy, izz, icomp, 2)**2
     1                    + V ( ixx, iyy, izz, icomp, 3)**2
     1                    + V ( ixx, iyy, izz, icomp, 4)**2 ) * ob12

! Construct high polynomial for MRWENO by suitably weighted subtracting off
! of the low polynomial.
          Phi9 ( 1:3) = X9 ( 1:3) * obgammahi
     1                - V ( ixx, iyy, izz, icomp, 2:4) * gammaratio
          Phi9 ( 4:9) = X9 ( 4:9) * obgammahi

! Smoothness indicator for high polynomial, with coefficients "Phi9 ( 1:9)".
          IShi = ( Phi9 ( 1)**2 + Phi9 ( 2)**2 + Phi9 ( 3)**2) * ob12

     1  + ( Phi9 ( 4)**2) * ( ( Coef_0 ( 5) - ob12)**2 + ob180)
     1  + ( Phi9 ( 5)**2) * ( ( Coef_0 ( 6) - ob12)**2 + ob180)
     1  + ( Phi9 ( 6)**2) * ( ( Coef_0 ( 7) - ob12)**2 + ob180)

     1  + ( Phi9 ( 7)**2) * ( Coef_0 ( 8)**2 + ob144)
     1  + ( Phi9 ( 8)**2) * ( Coef_0 ( 9)**2 + ob144)
     1  + ( Phi9 ( 9)**2) * ( Coef_0 ( 10)**2 + ob144)


! Construct non-linear weights using the linear weights and the smoothness indicators.

          tau = ABS ( IShi - ISlo)

          Wtlo = gammalo * (1.0 + tau / (ISlo + smallnum) )
          Wthi = gammahi * (1.0 + tau / (IShi + smallnum) )

          aa = 1.0 / ( Wtlo + Wthi)

          Wtlo = Wtlo * aa
          Wthi = Wthi * aa

! Non-linearly hybridized combination of low and high order polynomials.
! Note that the low order polynomial has only linear slopes.

          V ( ixx, iyy, izz, icomp, 2) = Wthi * Phi9 (1)
     1                            + Wtlo * V ( ixx, iyy, izz, icomp, 2)

          V ( ixx, iyy, izz, icomp, 3) = Wthi * Phi9 (2)
     1                            + Wtlo * V ( ixx, iyy, izz, icomp, 3)

          V ( ixx, iyy, izz, icomp, 4) = Wthi * Phi9 (3)
     1                            + Wtlo * V ( ixx, iyy, izz, icomp, 4)

          V ( ixx, iyy, izz, icomp, 5: 10) = Wthi * Phi9 ( 4: 9)

! ------------------------------
! Reduce the gradients if flattener is active.
! ("flattengr" is 0.0 near smooth flow and 1.0 in regions with v. strong
!   shocks and rarefactions.)
! ------------------------------

          x = 1.0 - flattengr ( ixx, iyy, izz)

          V ( ixx, iyy, izz, icomp, 2: n_cc_modes) =
     1    V ( ixx, iyy, izz, icomp, 2: n_cc_modes) * x

! --------------------------------------------------

          END DO    ! End of "DO icomp = 1, n_cc_components" loop.

! ----------------------------------------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = - 1, ntp_zones + 2" loop.

        END DO  ! End of "DO iyy = - 1, ntp_zones + 2" loop.

      END DO  ! End of "DO izz = - 1, nr_zones + 2" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 7.c) Enforce positivity of the reconstructed solution.

#if ( PRESSURE_POSITIVITY >= 1)

      CALL ENFORCE_PRESSURE_POSITIVITY_PRIM

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_components, n_cc_modes,
     1         n_fc_modes,

     1         vert_coords_patch, imask_active_patch,

     1         xf_centroid_patch, yf_centroid_patch, zf_centroid_patch,
     1         v_scale_length, v_centroid_patch,

     1         istencil_patch,

     1         V, entropygr, B1_zc, B2_zc, B3_zc, flattengr,
     1         Coef_Basis_Moments)

#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------








! ----------------------------------------------------------------------
! Remove this for normal code; only needed for heliospheric wind code.
#if (1 == 1)
! ----------------------------------------------------------------------
#if ( CHOOSE_APPL_HELIO == 1)
! ----------------------------------------------------------------------

! For Wind type problems, we have a density and pressure gradient in the
! density and pressure profiles. As a result, any time-variation in the
! solution at the first inner boundary of ghost zones can throw off the
! calculation.

! For that reason:-

! Set slopes of inner boundary ghost zones to zero, for the Heliosphere
! simulations. This is to prevent the pressure from touching floor at the
! inner boundary. For CGL, remember to zero out gradients in both pressures.
! ------------------------------------------------------------

      IF ( iz_patch .EQ. 1 ) THEN         ! For the inner most patch.

! ------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost,                &
!$OMP&    V, B1_zc, B2_zc, B3_zc, imask_active_patch)                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)
!$acc parallel vector_length(64) present (                              &
!$acc&    V, B1_zc, B2_zc, B3_zc, imask_active_patch)
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 1 - NGhost, 1
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

      IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! --------------------------------------------------

        V ( ixx, iyy, izz, :, 2: ) = 0.0

        B1_zc ( ixx, iyy, izz, 2:) = 0.0
        B2_zc ( ixx, iyy, izz, 2:) = 0.0
        B3_zc ( ixx, iyy, izz, 2:) = 0.0

! --------------------------------------------------

      END IF    ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN".

      END DO    ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO    ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO    ! End of "DO izz = 1 - NGhost, 1".

!$acc end parallel
!$OMP END PARALLEL DO
! ------------------------------------------------------------

      END IF    ! End of "IF ( iz_patch .EQ. 1 ) THEN".

! ------------------------------------------------------------
#endif
! End of "CHOOSE_APPL_HELIO == 1".
! ----------------------------------------------------------------------
#endif
! Remove this for normal code; only needed for heliospheric wind code.
! ----------------------------------------------------------------------










! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 8) Now we focus on facial reconstruction of the magnetic fields at
! select nodes where we will invoke the 2D RS.

! Therefore, the very expensive div-free volumetric reconstruction is not needed.

! To begin with, zero out all temporary accumulation variables.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 8.a) Access the x-faces. From either side of the x-face we use
! "V ( ixxr, iyyr, izz, 6:8, 1:10)" and "V ( ixxl, iyyl, izz, 6:8, 1:10)"
! to fill in "Bxf_quadpts ( ixx, iyy, izz, 1:5)". This is done in a way
! that ensures that the weighted average of "Bxf_quadpts ( ixx, iyy, izz, 1:5)"
! is equal to the x-facial mean field "Bxf ( ixx, iyy, izz, 1)".
! ----------------------------------------------------------------------
! First we find all the Bx components at all the quadrature points of each
! active x-face by using the reconstructed "U" from either side of that face.

      n_cc_10 = 10

! Short loop over all quadrature points. Do not parallelize this outer loop.
      DO j = 1, NBFLD_QUAD_PTS

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxl, iyyl, ixxr, iyyr,      &
!$OMP&    x, y, z, x2, y2, z2, bx, by, bz, bx2, by2, bz2,               &
!$OMP&    bdiff, nx, ny, nz, obl0, obl02, Coef_0, Coef_j)               &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, j, smallnum,    &
!$OMP&    n_cc_10,                                                      &
!$OMP&    imask_active_patch, istencil_patch, v_scale_length,           &
!$OMP&    xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,        &
!$OMP&    B1_zc, B2_zc, B3_zc,                                          &
!$OMP&    v_centroid_patch, U, V, Coef_Basis_Moments, Bxf, Bxf_quadpts) &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present( smallnum,                     &
!$acc&    n_cc_10,                                                      &
!$acc&    imask_active_patch, istencil_patch, v_scale_length,           &
!$acc&    xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,        &
!$acc&    B1_zc, B2_zc, B3_zc,                                          &
!$acc&    v_centroid_patch, U, V, Coef_Basis_Moments, Bxf, Bxf_quadpts)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxl, iyyl, ixxr, iyyr,              &
!$acc&    x, y, z, x2, y2, z2, bx, by, bz, bx2, by2, bz2,               &
!$acc&    bdiff, nx, ny, nz, obl0, obl02, Coef_0, Coef_j)

      DO izz = - 2, nr_zones + 3

        DO iyy = - 2, ntp_zones + 3

          DO ixx = - 3, ntp_zones + 3

! --------------------------------------------------

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

            Bxf_quadpts ( ixx, iyy, izz, j) = 0.0

            CYCLE  ! Go to the next iteration of the "DO ixx =" loop.

          END IF

! --------------------------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixxl, iyyl, izz), smallnum)  ! Scaling for Lf zone.

          obl02 = 1.0
     1  / AMAX1 ( v_scale_length ( ixxr, iyyr, izz), smallnum)  ! Scaling for Rt zone.

          Coef_0 ( 1: n_cc_10) = Coef_Basis_Moments ( ixxl, iyyl, izz,
     1                           1: n_cc_10)

          Coef_j ( 1: n_cc_10) = Coef_Basis_Moments ( ixxr, iyyr, izz,
     1                           1: n_cc_10)

! ------------------------------
! Obtain "nx, ny, nz" the normal to this face. The entire x-face has
! the same normal.

          nx = xb_quadnorml_patch ( ixx, iyy, izz, 1)
          ny = xb_quadnorml_patch ( ixx, iyy, izz, 2)
          nz = xb_quadnorml_patch ( ixx, iyy, izz, 3)

! ------------------------------
! Build scaled relative coordinates from Lf zone centroid and Rt zone centroid to
! "j" facial quadrature point.

          x = ( xb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixxl, iyyl, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixxl, iyyl, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixxl, iyyl, izz, 3) ) * obl0

! ----------

          x2 = ( xb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixxr, iyyr, izz, 1) ) * obl02

          y2 = ( xb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixxr, iyyr, izz, 2) ) * obl02

          z2 = ( xb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixxr, iyyr, izz, 3) ) * obl02

! ------------------------------
! Use the above scaled relative coordinates to obtain "bx, by, bz" from
! the Lf zone and "bx2, by2, bz2" from the Rt zone.

          bx = V ( ixxl, iyyl, izz, 6, 1)
     1       + V ( ixxl, iyyl, izz, 6, 2) * x
     1       + V ( ixxl, iyyl, izz, 6, 3) * y
     1       + V ( ixxl, iyyl, izz, 6, 4) * z
     1       + V ( ixxl, iyyl, izz, 6, 5) * ( x**2 - Coef_0 ( 5))
     1       + V ( ixxl, iyyl, izz, 6, 6) * ( y**2 - Coef_0 ( 6))
     1       + V ( ixxl, iyyl, izz, 6, 7) * ( z**2 - Coef_0 ( 7))
     1       + V ( ixxl, iyyl, izz, 6, 8) * ( x * y - Coef_0 ( 8))
     1       + V ( ixxl, iyyl, izz, 6, 9) * ( y * z - Coef_0 ( 9))
     1       + V ( ixxl, iyyl, izz, 6, 10) * ( x * z - Coef_0 ( 10))

          by = V ( ixxl, iyyl, izz, 7, 1)
     1       + V ( ixxl, iyyl, izz, 7, 2) * x
     1       + V ( ixxl, iyyl, izz, 7, 3) * y
     1       + V ( ixxl, iyyl, izz, 7, 4) * z
     1       + V ( ixxl, iyyl, izz, 7, 5) * ( x**2 - Coef_0 ( 5))
     1       + V ( ixxl, iyyl, izz, 7, 6) * ( y**2 - Coef_0 ( 6))
     1       + V ( ixxl, iyyl, izz, 7, 7) * ( z**2 - Coef_0 ( 7))
     1       + V ( ixxl, iyyl, izz, 7, 8) * ( x * y - Coef_0 ( 8))
     1       + V ( ixxl, iyyl, izz, 7, 9) * ( y * z - Coef_0 ( 9))
     1       + V ( ixxl, iyyl, izz, 7, 10) * ( x * z - Coef_0 ( 10))

          bz = V ( ixxl, iyyl, izz, 8, 1)
     1       + V ( ixxl, iyyl, izz, 8, 2) * x
     1       + V ( ixxl, iyyl, izz, 8, 3) * y
     1       + V ( ixxl, iyyl, izz, 8, 4) * z
     1       + V ( ixxl, iyyl, izz, 8, 5) * ( x**2 - Coef_0 ( 5))
     1       + V ( ixxl, iyyl, izz, 8, 6) * ( y**2 - Coef_0 ( 6))
     1       + V ( ixxl, iyyl, izz, 8, 7) * ( z**2 - Coef_0 ( 7))
     1       + V ( ixxl, iyyl, izz, 8, 8) * ( x * y - Coef_0 ( 8))
     1       + V ( ixxl, iyyl, izz, 8, 9) * ( y * z - Coef_0 ( 9))
     1       + V ( ixxl, iyyl, izz, 8, 10) * ( x * z - Coef_0 ( 10))

! ----------

          bx2 = V ( ixxr, iyyr, izz, 6, 1)
     1        + V ( ixxr, iyyr, izz, 6, 2) * x2
     1        + V ( ixxr, iyyr, izz, 6, 3) * y2
     1        + V ( ixxr, iyyr, izz, 6, 4) * z2
     1        + V ( ixxr, iyyr, izz, 6, 5) * ( x2**2 - Coef_j ( 5))
     1        + V ( ixxr, iyyr, izz, 6, 6) * ( y2**2 - Coef_j ( 6))
     1        + V ( ixxr, iyyr, izz, 6, 7) * ( z2**2 - Coef_j ( 7))
     1        + V ( ixxr, iyyr, izz, 6, 8) * ( x2 * y2 - Coef_j ( 8))
     1        + V ( ixxr, iyyr, izz, 6, 9) * ( y2 * z2 - Coef_j ( 9))
     1        + V ( ixxr, iyyr, izz, 6, 10) * ( x2 * z2 - Coef_j ( 10))

          by2 = V ( ixxr, iyyr, izz, 7, 1)
     1        + V ( ixxr, iyyr, izz, 7, 2) * x2
     1        + V ( ixxr, iyyr, izz, 7, 3) * y2
     1        + V ( ixxr, iyyr, izz, 7, 4) * z2
     1        + V ( ixxr, iyyr, izz, 7, 5) * ( x2**2 - Coef_j ( 5))
     1        + V ( ixxr, iyyr, izz, 7, 6) * ( y2**2 - Coef_j ( 6))
     1        + V ( ixxr, iyyr, izz, 7, 7) * ( z2**2 - Coef_j ( 7))
     1        + V ( ixxr, iyyr, izz, 7, 8) * ( x2 * y2 - Coef_j ( 8))
     1        + V ( ixxr, iyyr, izz, 7, 9) * ( y2 * z2 - Coef_j ( 9))
     1        + V ( ixxr, iyyr, izz, 7, 10) * ( x2 * z2 - Coef_j ( 10))

          bz2 = V ( ixxr, iyyr, izz, 8, 1)
     1        + V ( ixxr, iyyr, izz, 8, 2) * x2
     1        + V ( ixxr, iyyr, izz, 8, 3) * y2
     1        + V ( ixxr, iyyr, izz, 8, 4) * z2
     1        + V ( ixxr, iyyr, izz, 8, 5) * ( x2**2 - Coef_j ( 5))
     1        + V ( ixxr, iyyr, izz, 8, 6) * ( y2**2 - Coef_j ( 6))
     1        + V ( ixxr, iyyr, izz, 8, 7) * ( z2**2 - Coef_j ( 7))
     1        + V ( ixxr, iyyr, izz, 8, 8) * ( x2 * y2 - Coef_j ( 8))
     1        + V ( ixxr, iyyr, izz, 8, 9) * ( y2 * z2 - Coef_j ( 9))
     1        + V ( ixxr, iyyr, izz, 8, 10) * ( x2 * z2 - Coef_j ( 10))

! ------------------------------
! Using "nx, ny, nz" and "bx, by, bz" we can find the normal component of
! the magnetic field from the Lf.
! Using "nx, ny, nz" and "bx2, by2, bz2" we can find the normal component of
! the magnetic field from the Rt.
! These can be averaged together and stored at the appropriate quadrature point
! "j" in "Bxf_quadpts ( ixx, iyy, izz, j)".

          Bxf_quadpts ( ixx, iyy, izz, j) = 0.5 * (
     1    nx * ( bx + bx2 ) + ny * ( by + by2 ) + nz * ( bz + bz2 ) )

! --------------------------------------------------

          END DO  ! End of "DO ixx = - 3, ntp_zones + 3" loop.

        END DO  ! End of "DO iyy = - 2, ntp_zones + 2" loop.

      END DO  ! End of "DO izz = - 2, nr_zones + 3" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, NBFLD_QUAD_PTS" loop.

! ------------------------------------------------------------
! Now that all the Bx components have been found at all the quadrature points
! in this row of x-faces, we make sure that the weighted average of
! "Bxf_quadpts ( ixx, iyy, izz, j)" is equal to "Bxf ( ixx, iyy, izz, 1)".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, bdiff)                       &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    xb_quadwt_patch, Bxf, Bxf_quadpts)                            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    xb_quadwt_patch, Bxf, Bxf_quadpts)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, bdiff)

      DO izz = - 2, nr_zones + 3

        DO iyy = - 2, ntp_zones + 3

          DO ixx = - 3, ntp_zones + 3

! --------------------------------------------------

#if (0 == 1)
! This will just turn all magnetic fields to first order; just for debugging purposes.
!          Bxf_quadpts ( ixx, iyy, izz, 1) = Bxf ( ixx, iyy, izz, 1)
!          Bxf_quadpts ( ixx, iyy, izz, 2) = Bxf ( ixx, iyy, izz, 1)
!          Bxf_quadpts ( ixx, iyy, izz, 3) = Bxf ( ixx, iyy, izz, 1)
!          Bxf_quadpts ( ixx, iyy, izz, 4) = Bxf ( ixx, iyy, izz, 1)
!          Bxf_quadpts ( ixx, iyy, izz, 5) = Bxf ( ixx, iyy, izz, 1)
#endif

          bdiff = Bxf ( ixx, iyy, izz, 1)

     1          - xb_quadwt_patch ( ixx, iyy, izz, 1)
     1          * Bxf_quadpts ( ixx, iyy, izz, 1)

     1          - xb_quadwt_patch ( ixx, iyy, izz, 2)
     1          * Bxf_quadpts ( ixx, iyy, izz, 2)

     1          - xb_quadwt_patch ( ixx, iyy, izz, 3)
     1          * Bxf_quadpts ( ixx, iyy, izz, 3)

     1          - xb_quadwt_patch ( ixx, iyy, izz, 4)
     1          * Bxf_quadpts ( ixx, iyy, izz, 4)

     1          - xb_quadwt_patch ( ixx, iyy, izz, 5)
     1          * Bxf_quadpts ( ixx, iyy, izz, 5)

          Bxf_quadpts ( ixx, iyy, izz, 1) =
     1    Bxf_quadpts ( ixx, iyy, izz, 1) + bdiff

          Bxf_quadpts ( ixx, iyy, izz, 2) =
     1    Bxf_quadpts ( ixx, iyy, izz, 2) + bdiff

          Bxf_quadpts ( ixx, iyy, izz, 3) =
     1    Bxf_quadpts ( ixx, iyy, izz, 3) + bdiff

          Bxf_quadpts ( ixx, iyy, izz, 4) =
     1    Bxf_quadpts ( ixx, iyy, izz, 4) + bdiff

          Bxf_quadpts ( ixx, iyy, izz, 5) =
     1    Bxf_quadpts ( ixx, iyy, izz, 5) + bdiff

! --------------------------------------------------

          END DO  ! End of "DO ixx = - 3, ntp_zones + 3" loop.

        END DO  ! End of "DO iyy = - 2, ntp_zones + 3" loop.

      END DO  ! End of "DO izz = - 2, nr_zones + 3" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 8.b) Access the y-faces. From either side of the y-face we use
! "V ( ixxu, iyyu, izz, 6:8, 1:10)" and "V ( ixxd, iyyd, izz, 6:8, 1:10)"
! to fill in "Byf_quadpts ( ixx, iyy, izz, 1:5)". This is done in a way
! that ensures that the weighted average of "Byf_quadpts ( ixx, iyy, izz, 1:5)"
! is equal to the y-facial mean field "Byf ( ixx, iyy, izz, 1)".
! ----------------------------------------------------------------------
! First we find all the By components at all the quadrature points of each
! active y-face by using the reconstructed "U" from either side of that face.

      n_cc_10 = 10

! Short loop over all quadrature points. Do not parallelize this outer loop.
      DO j = 1, NBFLD_QUAD_PTS

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxd, iyyd, ixxu, iyyu,      &
!$OMP&    x, y, z, x2, y2, z2, bx, by, bz, bx2, by2, bz2,               &
!$OMP&    bdiff, nx, ny, nz, obl0, obl02, Coef_0, Coef_j)               &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, j, smallnum,    &
!$OMP&    n_cc_10,                                                      &
!$OMP&    imask_active_patch, istencil_patch, v_scale_length,           &
!$OMP&    yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,        &
!$OMP&    B1_zc, B2_zc, B3_zc,                                          &
!$OMP&    v_centroid_patch, U, V, Coef_Basis_Moments, Byf, Byf_quadpts) &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present( smallnum,                     &
!$acc&    n_cc_10,                                                      &
!$acc&    imask_active_patch, istencil_patch, v_scale_length,           &
!$acc&    yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,        &
!$acc&    B1_zc, B2_zc, B3_zc,                                          &
!$acc&    v_centroid_patch, U, V, Coef_Basis_Moments, Byf, Byf_quadpts)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxd, iyyd, ixxu, iyyu,              &
!$acc&    x, y, z, x2, y2, z2, bx, by, bz, bx2, by2, bz2,               &
!$acc&    bdiff, nx, ny, nz, obl0, obl02, Coef_0, Coef_j)

      DO izz = - 2, nr_zones + 3

        DO iyy = - 3, ntp_zones + 3

          DO ixx = - 2, ntp_zones + 3

! --------------------------------------------------

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

            Byf_quadpts ( ixx, iyy, izz, j) = 0.0

            CYCLE  ! Go to the next iteration of the "DO ixx =" loop.

          END IF

! --------------------------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixxd, iyyd, izz), smallnum)  ! Scaling for Dn zone.

          obl02 = 1.0
     1  / AMAX1 ( v_scale_length ( ixxu, iyyu, izz), smallnum)  ! Scaling for Up zone.

          Coef_0 ( 1: n_cc_10) = Coef_Basis_Moments ( ixxd, iyyd, izz,
     1                              1: n_cc_10)

          Coef_j ( 1: n_cc_10) = Coef_Basis_Moments ( ixxu, iyyu, izz,
     1                              1: n_cc_10)

! ------------------------------
! Obtain "nx, ny, nz" the normal to this face. The entire y-face has
! the same normal.

          nx = yb_quadnorml_patch ( ixx, iyy, izz, 1)
          ny = yb_quadnorml_patch ( ixx, iyy, izz, 2)
          nz = yb_quadnorml_patch ( ixx, iyy, izz, 3)

! ------------------------------
! Build scaled relative coordinates from Dn zone centroid and Up zone centroid to
! "j" facial quadrature point.

          x = ( yb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixxd, iyyd, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixxd, iyyd, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixxd, iyyd, izz, 3) ) * obl0

! ----------

          x2 = ( yb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixxu, iyyu, izz, 1) ) * obl02

          y2 = ( yb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixxu, iyyu, izz, 2) ) * obl02

          z2 = ( yb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixxu, iyyu, izz, 3) ) * obl02

! ------------------------------
! Use the above scaled relative coordinates to obtain "bx, by, bz" from
! the Dn zone and "bx2, by2, bz2" from the Up zone.

          bx = V ( ixxd, iyyd, izz, 6, 1)
     1       + V ( ixxd, iyyd, izz, 6, 2) * x
     1       + V ( ixxd, iyyd, izz, 6, 3) * y
     1       + V ( ixxd, iyyd, izz, 6, 4) * z
     1       + V ( ixxd, iyyd, izz, 6, 5) * ( x**2 - Coef_0 ( 5))
     1       + V ( ixxd, iyyd, izz, 6, 6) * ( y**2 - Coef_0 ( 6))
     1       + V ( ixxd, iyyd, izz, 6, 7) * ( z**2 - Coef_0 ( 7))
     1       + V ( ixxd, iyyd, izz, 6, 8) * ( x * y - Coef_0 ( 8))
     1       + V ( ixxd, iyyd, izz, 6, 9) * ( y * z - Coef_0 ( 9))
     1       + V ( ixxd, iyyd, izz, 6, 10) * ( x * z - Coef_0 ( 10))

          by = V ( ixxd, iyyd, izz, 7, 1)
     1       + V ( ixxd, iyyd, izz, 7, 2) * x
     1       + V ( ixxd, iyyd, izz, 7, 3) * y
     1       + V ( ixxd, iyyd, izz, 7, 4) * z
     1       + V ( ixxd, iyyd, izz, 7, 5) * ( x**2 - Coef_0 ( 5))
     1       + V ( ixxd, iyyd, izz, 7, 6) * ( y**2 - Coef_0 ( 6))
     1       + V ( ixxd, iyyd, izz, 7, 7) * ( z**2 - Coef_0 ( 7))
     1       + V ( ixxd, iyyd, izz, 7, 8) * ( x * y - Coef_0 ( 8))
     1       + V ( ixxd, iyyd, izz, 7, 9) * ( y * z - Coef_0 ( 9))
     1       + V ( ixxd, iyyd, izz, 7, 10) * ( x * z - Coef_0 ( 10))

          bz = V ( ixxd, iyyd, izz, 8, 1)
     1       + V ( ixxd, iyyd, izz, 8, 2) * x
     1       + V ( ixxd, iyyd, izz, 8, 3) * y
     1       + V ( ixxd, iyyd, izz, 8, 4) * z
     1       + V ( ixxd, iyyd, izz, 8, 5) * ( x**2 - Coef_0 ( 5))
     1       + V ( ixxd, iyyd, izz, 8, 6) * ( y**2 - Coef_0 ( 6))
     1       + V ( ixxd, iyyd, izz, 8, 7) * ( z**2 - Coef_0 ( 7))
     1       + V ( ixxd, iyyd, izz, 8, 8) * ( x * y - Coef_0 ( 8))
     1       + V ( ixxd, iyyd, izz, 8, 9) * ( y * z - Coef_0 ( 9))
     1       + V ( ixxd, iyyd, izz, 8, 10) * ( x * z - Coef_0 ( 10))

! ----------

          bx2 = V ( ixxu, iyyu, izz, 6, 1)
     1        + V ( ixxu, iyyu, izz, 6, 2) * x2
     1        + V ( ixxu, iyyu, izz, 6, 3) * y2
     1        + V ( ixxu, iyyu, izz, 6, 4) * z2
     1        + V ( ixxu, iyyu, izz, 6, 5) * ( x2**2 - Coef_j ( 5))
     1        + V ( ixxu, iyyu, izz, 6, 6) * ( y2**2 - Coef_j ( 6))
     1        + V ( ixxu, iyyu, izz, 6, 7) * ( z2**2 - Coef_j ( 7))
     1        + V ( ixxu, iyyu, izz, 6, 8) * ( x2 * y2 - Coef_j ( 8))
     1        + V ( ixxu, iyyu, izz, 6, 9) * ( y2 * z2 - Coef_j ( 9))
     1        + V ( ixxu, iyyu, izz, 6, 10) * ( x2 * z2 - Coef_j ( 10))

          by2 = V ( ixxu, iyyu, izz, 7, 1)
     1        + V ( ixxu, iyyu, izz, 7, 2) * x2
     1        + V ( ixxu, iyyu, izz, 7, 3) * y2
     1        + V ( ixxu, iyyu, izz, 7, 4) * z2
     1        + V ( ixxu, iyyu, izz, 7, 5) * ( x2**2 - Coef_j ( 5))
     1        + V ( ixxu, iyyu, izz, 7, 6) * ( y2**2 - Coef_j ( 6))
     1        + V ( ixxu, iyyu, izz, 7, 7) * ( z2**2 - Coef_j ( 7))
     1        + V ( ixxu, iyyu, izz, 7, 8) * ( x2 * y2 - Coef_j ( 8))
     1        + V ( ixxu, iyyu, izz, 7, 9) * ( y2 * z2 - Coef_j ( 9))
     1        + V ( ixxu, iyyu, izz, 7, 10) * ( x2 * z2 - Coef_j ( 10))

          bz2 = V ( ixxu, iyyu, izz, 8, 1)
     1        + V ( ixxu, iyyu, izz, 8, 2) * x2
     1        + V ( ixxu, iyyu, izz, 8, 3) * y2
     1        + V ( ixxu, iyyu, izz, 8, 4) * z2
     1        + V ( ixxu, iyyu, izz, 8, 5) * ( x2**2 - Coef_j ( 5))
     1        + V ( ixxu, iyyu, izz, 8, 6) * ( y2**2 - Coef_j ( 6))
     1        + V ( ixxu, iyyu, izz, 8, 7) * ( z2**2 - Coef_j ( 7))
     1        + V ( ixxu, iyyu, izz, 8, 8) * ( x2 * y2 - Coef_j ( 8))
     1        + V ( ixxu, iyyu, izz, 8, 9) * ( y2 * z2 - Coef_j ( 9))
     1        + V ( ixxu, iyyu, izz, 8, 10) * ( x2 * z2 - Coef_j ( 10))

! ------------------------------
! Using "nx, ny, nz" and "bx, by, bz" we can find the normal component of
! the magnetic field from the Dn.
! Using "nx, ny, nz" and "bx2, by2, bz2" we can find the normal component of
! the magnetic field from the Up.
! These can be averaged together and stored at the appropriate quadrature point
! "j" in "Byf_quadpts ( ixx, iyy, izz, j)".

          Byf_quadpts ( ixx, iyy, izz, j) = 0.5 * (
     1    nx * ( bx + bx2 ) + ny * ( by + by2 ) + nz * ( bz + bz2 ) )

! --------------------------------------------------

          END DO  ! End of "DO ixx = - 2, ntp_zones + 3" loop.

        END DO  ! End of "DO iyy = - 3, ntp_zones + 3" loop.

      END DO  ! End of "DO izz = - 2, nr_zones + 3" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, NBFLD_QUAD_PTS" loop.

! ------------------------------------------------------------
! Now that all the By components have been found at all the quadrature points
! in this row of y-faces, we make sure that the weighted average of
! "Byf_quadpts ( ixx, iyy, izz, j)" is equal to "Byf ( ixx, iyy, izz, 1)".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, bdiff)                       &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    yb_quadwt_patch, Byf, Byf_quadpts)                            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&     yb_quadwt_patch, Byf, Byf_quadpts)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, bdiff)

      DO izz = - 2, nr_zones + 3

        DO iyy = - 3, ntp_zones + 3

          DO ixx = - 2, ntp_zones + 3

! --------------------------------------------------

#if (0 == 1)
! This will just turn all magnetic fields to first order; just for debugging purposes.
!          Byf_quadpts ( ixx, iyy, izz, 1) = Byf ( ixx, iyy, izz, 1)
!          Byf_quadpts ( ixx, iyy, izz, 2) = Byf ( ixx, iyy, izz, 1)
!          Byf_quadpts ( ixx, iyy, izz, 3) = Byf ( ixx, iyy, izz, 1)
!          Byf_quadpts ( ixx, iyy, izz, 4) = Byf ( ixx, iyy, izz, 1)
!          Byf_quadpts ( ixx, iyy, izz, 5) = Byf ( ixx, iyy, izz, 1)
#endif

          bdiff = Byf ( ixx, iyy, izz, 1)

     1          - yb_quadwt_patch ( ixx, iyy, izz, 1)
     1          * Byf_quadpts ( ixx, iyy, izz, 1)

     1          - yb_quadwt_patch ( ixx, iyy, izz, 2)
     1          * Byf_quadpts ( ixx, iyy, izz, 2)

     1          - yb_quadwt_patch ( ixx, iyy, izz, 3)
     1          * Byf_quadpts ( ixx, iyy, izz, 3)

     1          - yb_quadwt_patch ( ixx, iyy, izz, 4)
     1          * Byf_quadpts ( ixx, iyy, izz, 4)

     1          - yb_quadwt_patch ( ixx, iyy, izz, 5)
     1          * Byf_quadpts ( ixx, iyy, izz, 5)

          Byf_quadpts ( ixx, iyy, izz, 1) =
     1    Byf_quadpts ( ixx, iyy, izz, 1) + bdiff

          Byf_quadpts ( ixx, iyy, izz, 2) =
     1    Byf_quadpts ( ixx, iyy, izz, 2) + bdiff

          Byf_quadpts ( ixx, iyy, izz, 3) =
     1    Byf_quadpts ( ixx, iyy, izz, 3) + bdiff

          Byf_quadpts ( ixx, iyy, izz, 4) =
     1    Byf_quadpts ( ixx, iyy, izz, 4) + bdiff

          Byf_quadpts ( ixx, iyy, izz, 5) =
     1    Byf_quadpts ( ixx, iyy, izz, 5) + bdiff

! --------------------------------------------------

          END DO  ! End of "DO ixx = - 2, ntp_zones + 3" loop.

        END DO  ! End of "DO iyy = - 3, ntp_zones + 3" loop.

      END DO  ! End of "DO izz = - 2, nr_zones + 3" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 8.c) Access the z-faces. From either side of the z-face we use
! "V ( ixx, iyy, izz, 6:8, 1:10)" and "V ( ixx, iyy, izzp1, 6:8, 1:10)"
! to fill in "Bzf_quadpts ( ixx, iyy, izz, 1:5)". This is done in a way
! that ensures that the weighted average of "Bzf_quadpts ( ixx, iyy, izz, 1:5)"
! is equal to the z-facial mean field "Bzf ( ixx, iyy, izz, 1)".
! ----------------------------------------------------------------------
! First we find all the Bz components at all the quadrature points of each
! active z-face by using the reconstructed "U" from either side of that face.

      n_cc_10 = 10

! Short loop over all quadrature points. Do not parallelize this outer loop.
      DO j = 1, NBFLD_QUAD_PTS

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izzp1,                       &
!$OMP&    x, y, z, x2, y2, z2, bx, by, bz, bx2, by2, bz2,               &
!$OMP&    bdiff, nx, ny, nz, obl0, obl02, Coef_0, Coef_j)               &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, j, smallnum,    &
!$OMP&    n_cc_10,                                                      &
!$OMP&    imask_active_patch, istencil_patch, v_scale_length,           &
!$OMP&    zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch,        &
!$OMP&    B1_zc, B2_zc, B3_zc,                                          &
!$OMP&    v_centroid_patch, U, V, Coef_Basis_Moments, Bzf, Bzf_quadpts) &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present( smallnum,                     &
!$acc&    n_cc_10,                                                      &
!$acc&    imask_active_patch, istencil_patch, v_scale_length,           &
!$acc&    zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch,        &
!$acc&    B1_zc, B2_zc, B3_zc,                                          &
!$acc&    v_centroid_patch, U, V, Coef_Basis_Moments, Bzf, Bzf_quadpts)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, izzp1,                               &
!$acc&    x, y, z, x2, y2, z2, bx, by, bz, bx2, by2, bz2,               &
!$acc&    bdiff, nx, ny, nz, obl0, obl02, Coef_0, Coef_j)

      DO izz = - 3, nr_zones + 3

        DO iyy = - 2, ntp_zones + 3

          DO ixx = - 2, ntp_zones + 3

          izzp1 = izz + 1

! --------------------------------------------------

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) THEN

            Bzf_quadpts ( ixx, iyy, izz, j) = 0.0

            CYCLE  ! Go to the next iteration of the "DO ixx =" loop.

          END IF

! --------------------------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for "izz" zone.

          obl02 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izzp1), smallnum)  ! Scaling for "izzp1" zone.

          Coef_0 ( 1: n_cc_10) = Coef_Basis_Moments ( ixx, iyy, izz,
     1                              1: n_cc_10)

          Coef_j ( 1: n_cc_10) = Coef_Basis_Moments ( ixx, iyy, izzp1,
     1                              1: n_cc_10)

! ------------------------------
! Obtain "nx, ny, nz" the normal to this face. Different quadrature points
! on the same z-face will have slightly differently pointing normals.
! Notice the dimensioning "zb_quadnorml_patch ( ixx, iyy, izz, 1:3, j)"; i.e.
! it depends on the quadrature point "j".

          nx = zb_quadnorml_patch ( ixx, iyy, izz, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izz, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izz, 3, j)

! ------------------------------
! Build scaled relative coordinates from "izz" zone centroid and "izzp1" zone
! centroid to "j" facial quadrature point.

          x = ( zb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! ----------

          x2 = ( zb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1      - v_centroid_patch ( ixx, iyy, izzp1, 1) ) * obl02

          y2 = ( zb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1      - v_centroid_patch ( ixx, iyy, izzp1, 2) ) * obl02

          z2 = ( zb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1      - v_centroid_patch ( ixx, iyy, izzp1, 3) ) * obl02

! ------------------------------
! Use the above scaled relative coordinates to obtain "bx, by, bz" from
! the "izz" zone and "bx2, by2, bz2" from the "izzp1" zone.

          bx = V ( ixx, iyy, izz, 6, 1)
     1       + V ( ixx, iyy, izz, 6, 2) * x
     1       + V ( ixx, iyy, izz, 6, 3) * y
     1       + V ( ixx, iyy, izz, 6, 4) * z
     1       + V ( ixx, iyy, izz, 6, 5) * ( x**2 - Coef_0 ( 5))
     1       + V ( ixx, iyy, izz, 6, 6) * ( y**2 - Coef_0 ( 6))
     1       + V ( ixx, iyy, izz, 6, 7) * ( z**2 - Coef_0 ( 7))
     1       + V ( ixx, iyy, izz, 6, 8) * ( x * y - Coef_0 ( 8))
     1       + V ( ixx, iyy, izz, 6, 9) * ( y * z - Coef_0 ( 9))
     1       + V ( ixx, iyy, izz, 6, 10) * ( x * z - Coef_0 ( 10))

          by = V ( ixx, iyy, izz, 7, 1)
     1       + V ( ixx, iyy, izz, 7, 2) * x
     1       + V ( ixx, iyy, izz, 7, 3) * y
     1       + V ( ixx, iyy, izz, 7, 4) * z
     1       + V ( ixx, iyy, izz, 7, 5) * ( x**2 - Coef_0 ( 5))
     1       + V ( ixx, iyy, izz, 7, 6) * ( y**2 - Coef_0 ( 6))
     1       + V ( ixx, iyy, izz, 7, 7) * ( z**2 - Coef_0 ( 7))
     1       + V ( ixx, iyy, izz, 7, 8) * ( x * y - Coef_0 ( 8))
     1       + V ( ixx, iyy, izz, 7, 9) * ( y * z - Coef_0 ( 9))
     1       + V ( ixx, iyy, izz, 7, 10) * ( x * z - Coef_0 ( 10))

          bz = V ( ixx, iyy, izz, 8, 1)
     1       + V ( ixx, iyy, izz, 8, 2) * x
     1       + V ( ixx, iyy, izz, 8, 3) * y
     1       + V ( ixx, iyy, izz, 8, 4) * z
     1       + V ( ixx, iyy, izz, 8, 5) * ( x**2 - Coef_0 ( 5))
     1       + V ( ixx, iyy, izz, 8, 6) * ( y**2 - Coef_0 ( 6))
     1       + V ( ixx, iyy, izz, 8, 7) * ( z**2 - Coef_0 ( 7))
     1       + V ( ixx, iyy, izz, 8, 8) * ( x * y - Coef_0 ( 8))
     1       + V ( ixx, iyy, izz, 8, 9) * ( y * z - Coef_0 ( 9))
     1       + V ( ixx, iyy, izz, 8, 10) * ( x * z - Coef_0 ( 10))

! ----------

          bx2 = V ( ixx, iyy, izzp1, 6, 1)
     1        + V ( ixx, iyy, izzp1, 6, 2) * x2
     1        + V ( ixx, iyy, izzp1, 6, 3) * y2
     1        + V ( ixx, iyy, izzp1, 6, 4) * z2
     1        + V ( ixx, iyy, izzp1, 6, 5) * ( x2**2 - Coef_j ( 5))
     1        + V ( ixx, iyy, izzp1, 6, 6) * ( y2**2 - Coef_j ( 6))
     1        + V ( ixx, iyy, izzp1, 6, 7) * ( z2**2 - Coef_j ( 7))
     1        + V ( ixx, iyy, izzp1, 6, 8) * ( x2 * y2 - Coef_j ( 8))
     1        + V ( ixx, iyy, izzp1, 6, 9) * ( y2 * z2 - Coef_j ( 9))
     1        + V ( ixx, iyy, izzp1, 6, 10) * ( x2 * z2 - Coef_j ( 10))

          by2 = V ( ixx, iyy, izzp1, 7, 1)
     1        + V ( ixx, iyy, izzp1, 7, 2) * x2
     1        + V ( ixx, iyy, izzp1, 7, 3) * y2
     1        + V ( ixx, iyy, izzp1, 7, 4) * z2
     1        + V ( ixx, iyy, izzp1, 7, 5) * ( x2**2 - Coef_j ( 5))
     1        + V ( ixx, iyy, izzp1, 7, 6) * ( y2**2 - Coef_j ( 6))
     1        + V ( ixx, iyy, izzp1, 7, 7) * ( z2**2 - Coef_j ( 7))
     1        + V ( ixx, iyy, izzp1, 7, 8) * ( x2 * y2 - Coef_j ( 8))
     1        + V ( ixx, iyy, izzp1, 7, 9) * ( y2 * z2 - Coef_j ( 9))
     1        + V ( ixx, iyy, izzp1, 7, 10) * ( x2 * z2 - Coef_j ( 10))

          bz2 = V ( ixx, iyy, izzp1, 8, 1)
     1        + V ( ixx, iyy, izzp1, 8, 2) * x2
     1        + V ( ixx, iyy, izzp1, 8, 3) * y2
     1        + V ( ixx, iyy, izzp1, 8, 4) * z2
     1        + V ( ixx, iyy, izzp1, 8, 5) * ( x2**2 - Coef_j ( 5))
     1        + V ( ixx, iyy, izzp1, 8, 6) * ( y2**2 - Coef_j ( 6))
     1        + V ( ixx, iyy, izzp1, 8, 7) * ( z2**2 - Coef_j ( 7))
     1        + V ( ixx, iyy, izzp1, 8, 8) * ( x2 * y2 - Coef_j ( 8))
     1        + V ( ixx, iyy, izzp1, 8, 9) * ( y2 * z2 - Coef_j ( 9))
     1        + V ( ixx, iyy, izzp1, 8, 10) * ( x2 * z2 - Coef_j ( 10))

! ------------------------------
! Using "nx, ny, nz" and "bx, by, bz" we can find the normal component of
! the magnetic field from the "izz" zone.
! Using "nx, ny, nz" and "bx2, by2, bz2" we can find the normal component of
! the magnetic field from the "izzp1" zone.
! These can be averaged together and stored at the appropriate quadrature point
! "j" in "Bzf_quadpts ( ixx, iyy, izz, j)".

          Bzf_quadpts ( ixx, iyy, izz, j) = 0.5 * (
     1    nx * ( bx + bx2 ) + ny * ( by + by2 ) + nz * ( bz + bz2 ) )

! --------------------------------------------------

          END DO  ! End of "DO ixx = - 2, ntp_zones + 3" loop.

        END DO  ! End of "DO iyy = - 2, ntp_zones + 3" loop.

      END DO  ! End of "DO izz = - 3, nr_zones + 3" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END DO  ! End of "DO j = 1, NBFLD_QUAD_PTS" loop.

! ------------------------------------------------------------
! Now that all the Bz components have been found at all the quadrature points
! in this row of z-faces, we make sure that the weighted average of
! "Bzf_quadpts ( ixx, iyy, izz, j)" is equal to "Bzf ( ixx, iyy, izz, 1)".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, bdiff)                       &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    zb_quadwt_patch, Bzf, Bzf_quadpts)                            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    zb_quadwt_patch, Bzf, Bzf_quadpts)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, bdiff)

      DO izz = - 3, nr_zones + 3

        DO iyy = - 2, ntp_zones + 3

          DO ixx = - 2, ntp_zones + 3

! --------------------------------------------------

#if (0 == 1)
! This will just turn all magnetic fields to first order; just for debugging purposes.
!          Bzf_quadpts ( ixx, iyy, izz, 1) = Bzf ( ixx, iyy, izz, 1)
!          Bzf_quadpts ( ixx, iyy, izz, 2) = Bzf ( ixx, iyy, izz, 1)
!          Bzf_quadpts ( ixx, iyy, izz, 3) = Bzf ( ixx, iyy, izz, 1)
!          Bzf_quadpts ( ixx, iyy, izz, 4) = Bzf ( ixx, iyy, izz, 1)
!          Bzf_quadpts ( ixx, iyy, izz, 5) = Bzf ( ixx, iyy, izz, 1)
#endif

          bdiff = Bzf ( ixx, iyy, izz, 1)

     1          - zb_quadwt_patch ( ixx, iyy, izz, 1)
     1          * Bzf_quadpts ( ixx, iyy, izz, 1)

     1          - zb_quadwt_patch ( ixx, iyy, izz, 2)
     1          * Bzf_quadpts ( ixx, iyy, izz, 2)

     1          - zb_quadwt_patch ( ixx, iyy, izz, 3)
     1          * Bzf_quadpts ( ixx, iyy, izz, 3)

     1          - zb_quadwt_patch ( ixx, iyy, izz, 4)
     1          * Bzf_quadpts ( ixx, iyy, izz, 4)

     1          - zb_quadwt_patch ( ixx, iyy, izz, 5)
     1          * Bzf_quadpts ( ixx, iyy, izz, 5)

          Bzf_quadpts ( ixx, iyy, izz, 1) =
     1    Bzf_quadpts ( ixx, iyy, izz, 1) + bdiff

          Bzf_quadpts ( ixx, iyy, izz, 2) =
     1    Bzf_quadpts ( ixx, iyy, izz, 2) + bdiff

          Bzf_quadpts ( ixx, iyy, izz, 3) =
     1    Bzf_quadpts ( ixx, iyy, izz, 3) + bdiff

          Bzf_quadpts ( ixx, iyy, izz, 4) =
     1    Bzf_quadpts ( ixx, iyy, izz, 4) + bdiff

          Bzf_quadpts ( ixx, iyy, izz, 5) =
     1    Bzf_quadpts ( ixx, iyy, izz, 5) + bdiff

! --------------------------------------------------

          END DO  ! End of "DO ixx = - 2, ntp_zones + 3" loop.

        END DO  ! End of "DO iyy = - 2, ntp_zones + 3" loop.

      END DO  ! End of "DO izz = - 3, nr_zones + 3" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 8.d) Now that we have "Bxf_quadpts, Byf_quadpts, Bzf_quadpts"
! at the six faces of each zone, we can make an Almost Div-Preserving reconstruction.
! of the magnetic field using a CLSQ formulation. This means that the ADP
! vector fields match the six facial averages exactly (the constraints)
! while matching the FV-reconstructed moments almost exactly (the LSQ part).

      n_cc_10 = 10

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, x, y, z, nx, ny, nz, obl0, j,&
!$OMP&   ixxm1, iyym1, izzm1, CMat_B, AMat_B, DVec_B, BVec_B, Coef_0,wt,&
!$OMP&   aa, bb, ATrans_B, ATransA_B, ATransB_B, CTrans_B,              &
!$OMP&    CLSQMat_B, CLSQVec_B, indx_CLSQ_B, vv_CLSQ_B, CMat_Row_B)     &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    smallnum, pi, ob_8pi, n_cc_10,                                &
!$OMP&    xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,        &
!$OMP&    xf_centroid_patch, Bxf, Bxf_quadpts, B1_zc,                   &
!$OMP&    yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,        &
!$OMP&    yf_centroid_patch, Byf, Byf_quadpts, B2_zc,                   &
!$OMP&    zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch,        &
!$OMP&    zf_centroid_patch, Bzf, Bzf_quadpts, B3_zc,                   &
!$OMP&    imask_active_patch, istencil_patch, v_scale_length,           &
!$OMP&    v_centroid_patch, U, V, Coef_Basis_Moments)                   &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 3)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum, pi, ob_8pi, n_cc_10,                                &
!$acc&    xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,        &
!$acc&    xf_centroid_patch, Bxf, Bxf_quadpts, B1_zc,                   &
!$acc&    yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,        &
!$acc&    yf_centroid_patch, Byf, Byf_quadpts, B2_zc,                   &
!$acc&    zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch,        &
!$acc&    zf_centroid_patch, Bzf, Bzf_quadpts, B3_zc,                   &
!$acc&    imask_active_patch, istencil_patch, v_scale_length,           &
!$acc&    v_centroid_patch, U, V, Coef_Basis_Moments)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, x, y, z, nx, ny, nz, obl0, j, wt,    &
!$acc&    ixxm1, iyym1, izzm1, CMat_B, AMat_B, DVec_B, BVec_B, Coef_0,  &
!$acc&    aa, bb, ATrans_B, ATransA_B, ATransB_B, CTrans_B,             &
!$acc&    CLSQMat_B, CLSQVec_B, indx_CLSQ_B, vv_CLSQ_B, CMat_Row_B)

      DO izz = - 2, nr_zones + 3

        DO iyy = - 2, ntp_zones + 3

          DO ixx = - 2, ntp_zones + 3

          ixxm1 = ixx - 1
          iyym1 = iyy - 1
          izzm1 = izz - 1

! --------------------------------------------------

! If the "imask_active_patch" is zero, this is not a valid zone, so we skip it.
          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) CYCLE

          CMat_B = 0.0 ; AMat_B = 0.0 ; DVec_B = 0.0 ; BVec_B = 0.0

          Coef_0 ( 1: n_cc_10) = Coef_Basis_Moments ( ixx, iyy, izz,
     1                             1: n_cc_10)

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

! --------------------------------------------------
! Obtain "nx, ny, nz" the normal to this right x-face. The entire x-face has
! the same normal.

          nx = xb_quadnorml_patch ( ixx, iyy, izz, 1)
          ny = xb_quadnorml_patch ( ixx, iyy, izz, 2)
          nz = xb_quadnorml_patch ( ixx, iyy, izz, 3)

! ----------

! Build scaled relative coordinates from zone centroid to all five facial quadrature points.

          j = 1

          wt = xb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( xb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 1, :) = CMat_B ( 1, :) + wt * CMat_Row_B ( :)
          DVec_B ( 1) = DVec_B ( 1) + wt
     1                              * Bxf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 1, :) = CMat_Row_B ( :)
          BVec_B ( 1) = Bxf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 2

          wt = xb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( xb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 1, :) = CMat_B ( 1, :) + wt * CMat_Row_B ( :)
          DVec_B ( 1) = DVec_B ( 1) + wt
     1                              * Bxf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 2, :) = CMat_Row_B ( :)
          BVec_B ( 2) = Bxf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 3

          wt = xb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( xb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 1, :) = CMat_B ( 1, :) + wt * CMat_Row_B ( :)
          DVec_B ( 1) = DVec_B ( 1) + wt
     1                              * Bxf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 3, :) = CMat_Row_B ( :)
          BVec_B ( 3) = Bxf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 4

          wt = xb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( xb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 1, :) = CMat_B ( 1, :) + wt * CMat_Row_B ( :)
          DVec_B ( 1) = DVec_B ( 1) + wt
     1                              * Bxf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 4, :) = CMat_Row_B ( :)
          BVec_B ( 4) = Bxf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 5

          wt = xb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( xb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 1, :) = CMat_B ( 1, :) + wt * CMat_Row_B ( :)
          DVec_B ( 1) = DVec_B ( 1) + wt
     1                              * Bxf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 5, :) = CMat_Row_B ( :)
          BVec_B ( 5) = Bxf_quadpts ( ixx, iyy, izz, j)

! ------------------------------
! Obtain "nx, ny, nz" the normal to this left x-face. The entire x-face has
! the same normal.

          nx = xb_quadnorml_patch ( ixxm1, iyy, izz, 1)
          ny = xb_quadnorml_patch ( ixxm1, iyy, izz, 2)
          nz = xb_quadnorml_patch ( ixxm1, iyy, izz, 3)

! ----------

! Build scaled relative coordinates from zone centroid to all five facial quadrature points.

          j = 1

          wt = xb_quadwt_patch ( ixxm1, iyy, izz, j)

          x = ( xb_quadpts_patch ( ixxm1, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixxm1, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixxm1, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 2, :) = CMat_B ( 2, :) + wt * CMat_Row_B ( :)
          DVec_B ( 2) = DVec_B ( 2) + wt
     1                              * Bxf_quadpts ( ixxm1, iyy, izz, j)

          AMat_B ( 6, :) = CMat_Row_B ( :)
          BVec_B ( 6) = Bxf_quadpts ( ixxm1, iyy, izz, j)

! ----------

          j = 2

          wt = xb_quadwt_patch ( ixxm1, iyy, izz, j)

          x = ( xb_quadpts_patch ( ixxm1, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixxm1, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixxm1, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 2, :) = CMat_B ( 2, :) + wt * CMat_Row_B ( :)
          DVec_B ( 2) = DVec_B ( 2) + wt
     1                              * Bxf_quadpts ( ixxm1, iyy, izz, j)

          AMat_B ( 7, :) = CMat_Row_B ( :)
          BVec_B ( 7) = Bxf_quadpts ( ixxm1, iyy, izz, j)

! ----------

          j = 3

          wt = xb_quadwt_patch ( ixxm1, iyy, izz, j)

          x = ( xb_quadpts_patch ( ixxm1, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixxm1, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixxm1, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 2, :) = CMat_B ( 2, :) + wt * CMat_Row_B ( :)
          DVec_B ( 2) = DVec_B ( 2) + wt
     1                              * Bxf_quadpts ( ixxm1, iyy, izz, j)

          AMat_B ( 8, :) = CMat_Row_B ( :)
          BVec_B ( 8) = Bxf_quadpts ( ixxm1, iyy, izz, j)

! ----------

          j = 4

          wt = xb_quadwt_patch ( ixxm1, iyy, izz, j)

          x = ( xb_quadpts_patch ( ixxm1, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixxm1, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixxm1, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 2, :) = CMat_B ( 2, :) + wt * CMat_Row_B ( :)
          DVec_B ( 2) = DVec_B ( 2) + wt
     1                              * Bxf_quadpts ( ixxm1, iyy, izz, j)

          AMat_B ( 9, :) = CMat_Row_B ( :)
          BVec_B ( 9) = Bxf_quadpts ( ixxm1, iyy, izz, j)

! ----------

          j = 5

          wt = xb_quadwt_patch ( ixxm1, iyy, izz, j)

          x = ( xb_quadpts_patch ( ixxm1, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xb_quadpts_patch ( ixxm1, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xb_quadpts_patch ( ixxm1, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 2, :) = CMat_B ( 2, :) + wt * CMat_Row_B ( :)
          DVec_B ( 2) = DVec_B ( 2) + wt
     1                              * Bxf_quadpts ( ixxm1, iyy, izz, j)

          AMat_B ( 10, :) = CMat_Row_B ( :)
          BVec_B ( 10) = Bxf_quadpts ( ixxm1, iyy, izz, j)

! --------------------------------------------------
! Obtain "nx, ny, nz" the normal to this up y-face. The entire y-face has
! the same normal.

          nx = yb_quadnorml_patch ( ixx, iyy, izz, 1)
          ny = yb_quadnorml_patch ( ixx, iyy, izz, 2)
          nz = yb_quadnorml_patch ( ixx, iyy, izz, 3)

! ----------

! Build scaled relative coordinates from zone centroid to all five facial quadrature points.

          j = 1

          wt = yb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( yb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 3, :) = CMat_B ( 3, :) + wt * CMat_Row_B ( :)
          DVec_B ( 3) = DVec_B ( 3) + wt
     1                              * Byf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 11, :) = CMat_Row_B ( :)
          BVec_B ( 11) = Byf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 2

          wt = yb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( yb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 3, :) = CMat_B ( 3, :) + wt * CMat_Row_B ( :)
          DVec_B ( 3) = DVec_B ( 3) + wt
     1                              * Byf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 12, :) = CMat_Row_B ( :)
          BVec_B ( 12) = Byf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 3

          wt = yb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( yb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 3, :) = CMat_B ( 3, :) + wt * CMat_Row_B ( :)
          DVec_B ( 3) = DVec_B ( 3) + wt
     1                              * Byf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 13, :) = CMat_Row_B ( :)
          BVec_B ( 13) = Byf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 4

          wt = yb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( yb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 3, :) = CMat_B ( 3, :) + wt * CMat_Row_B ( :)
          DVec_B ( 3) = DVec_B ( 3) + wt
     1                              * Byf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 14, :) = CMat_Row_B ( :)
          BVec_B ( 14) = Byf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 5

          wt = yb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( yb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 3, :) = CMat_B ( 3, :) + wt * CMat_Row_B ( :)
          DVec_B ( 3) = DVec_B ( 3) + wt
     1                              * Byf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 15, :) = CMat_Row_B ( :)
          BVec_B ( 15) = Byf_quadpts ( ixx, iyy, izz, j)

! --------------------------------------------------
! Obtain "nx, ny, nz" the normal to this down y-face. The entire y-face has
! the same normal.

          nx = yb_quadnorml_patch ( ixx, iyym1, izz, 1)
          ny = yb_quadnorml_patch ( ixx, iyym1, izz, 2)
          nz = yb_quadnorml_patch ( ixx, iyym1, izz, 3)

! ----------

! Build scaled relative coordinates from zone centroid to all five facial quadrature points.

          j = 1

          wt = yb_quadwt_patch ( ixx, iyym1, izz, j)

          x = ( yb_quadpts_patch ( ixx, iyym1, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyym1, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyym1, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 4, :) = CMat_B ( 4, :) + wt * CMat_Row_B ( :)
          DVec_B ( 4) = DVec_B ( 4) + wt
     1                              * Byf_quadpts ( ixx, iyym1, izz, j)

          AMat_B ( 16, :) = CMat_Row_B ( :)
          BVec_B ( 16) = Byf_quadpts ( ixx, iyym1, izz, j)

! ----------

          j = 2

          wt = yb_quadwt_patch ( ixx, iyym1, izz, j)

          x = ( yb_quadpts_patch ( ixx, iyym1, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyym1, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyym1, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 4, :) = CMat_B ( 4, :) + wt * CMat_Row_B ( :)
          DVec_B ( 4) = DVec_B ( 4) + wt
     1                              * Byf_quadpts ( ixx, iyym1, izz, j)

          AMat_B ( 17, :) = CMat_Row_B ( :)
          BVec_B ( 17) = Byf_quadpts ( ixx, iyym1, izz, j)

! ----------

          j = 3

          wt = yb_quadwt_patch ( ixx, iyym1, izz, j)

          x = ( yb_quadpts_patch ( ixx, iyym1, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyym1, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyym1, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 4, :) = CMat_B ( 4, :) + wt * CMat_Row_B ( :)
          DVec_B ( 4) = DVec_B ( 4) + wt
     1                              * Byf_quadpts ( ixx, iyym1, izz, j)

          AMat_B ( 18, :) = CMat_Row_B ( :)
          BVec_B ( 18) = Byf_quadpts ( ixx, iyym1, izz, j)

! ----------

          j = 4

          wt = yb_quadwt_patch ( ixx, iyym1, izz, j)

          x = ( yb_quadpts_patch ( ixx, iyym1, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyym1, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyym1, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 4, :) = CMat_B ( 4, :) + wt * CMat_Row_B ( :)
          DVec_B ( 4) = DVec_B ( 4) + wt
     1                              * Byf_quadpts ( ixx, iyym1, izz, j)

          AMat_B ( 19, :) = CMat_Row_B ( :)
          BVec_B ( 19) = Byf_quadpts ( ixx, iyym1, izz, j)

! ----------

          j = 5

          wt = yb_quadwt_patch ( ixx, iyym1, izz, j)

          x = ( yb_quadpts_patch ( ixx, iyym1, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yb_quadpts_patch ( ixx, iyym1, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yb_quadpts_patch ( ixx, iyym1, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 4, :) = CMat_B ( 4, :) + wt * CMat_Row_B ( :)
          DVec_B ( 4) = DVec_B ( 4) + wt
     1                              * Byf_quadpts ( ixx, iyym1, izz, j)

          AMat_B ( 20, :) = CMat_Row_B ( :)
          BVec_B ( 20) = Byf_quadpts ( ixx, iyym1, izz, j)

! --------------------------------------------------
! Obtain "nx, ny, nz" the normal to this top z-face. different z-quadrature
! points have different normals. The first normal is the central one.

! ----------

! Build scaled relative coordinates from zone centroid to all five facial quadrature points.
! Different quadrature points also have different normals.

          j = 1

          nx = zb_quadnorml_patch ( ixx, iyy, izz, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izz, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izz, 3, j)

          wt = zb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( zb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 5, :) = CMat_B ( 5, :) + wt * CMat_Row_B ( :)
          DVec_B ( 5) = DVec_B ( 5) + wt
     1                              * Bzf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 21, :) = CMat_Row_B ( :)
          BVec_B ( 21) = Bzf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 2

          nx = zb_quadnorml_patch ( ixx, iyy, izz, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izz, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izz, 3, j)

          wt = zb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( zb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 5, :) = CMat_B ( 5, :) + wt * CMat_Row_B ( :)
          DVec_B ( 5) = DVec_B ( 5) + wt
     1                              * Bzf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 22, :) = CMat_Row_B ( :)
          BVec_B ( 22) = Bzf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 3

          nx = zb_quadnorml_patch ( ixx, iyy, izz, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izz, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izz, 3, j)

          wt = zb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( zb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 5, :) = CMat_B ( 5, :) + wt * CMat_Row_B ( :)
          DVec_B ( 5) = DVec_B ( 5) + wt
     1                              * Bzf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 23, :) = CMat_Row_B ( :)
          BVec_B ( 23) = Bzf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 4

          nx = zb_quadnorml_patch ( ixx, iyy, izz, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izz, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izz, 3, j)

          wt = zb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( zb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 5, :) = CMat_B ( 5, :) + wt * CMat_Row_B ( :)
          DVec_B ( 5) = DVec_B ( 5) + wt
     1                              * Bzf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 24, :) = CMat_Row_B ( :)
          BVec_B ( 24) = Bzf_quadpts ( ixx, iyy, izz, j)

! ----------

          j = 5

          nx = zb_quadnorml_patch ( ixx, iyy, izz, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izz, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izz, 3, j)

          wt = zb_quadwt_patch ( ixx, iyy, izz, j)

          x = ( zb_quadpts_patch ( ixx, iyy, izz, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izz, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izz, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 5, :) = CMat_B ( 5, :) + wt * CMat_Row_B ( :)
          DVec_B ( 5) = DVec_B ( 5) + wt
     1                              * Bzf_quadpts ( ixx, iyy, izz, j)

          AMat_B ( 25, :) = CMat_Row_B ( :)
          BVec_B ( 25) = Bzf_quadpts ( ixx, iyy, izz, j)

! --------------------------------------------------
! Obtain "nx, ny, nz" the normal to this bottom z-face. different z-quadrature
! points have different normals. The first normal is the central one.

! ----------

! Build scaled relative coordinates from zone centroid to all five facial quadrature points.
! Different quadrature points also have different normals.

          j = 1

          nx = zb_quadnorml_patch ( ixx, iyy, izzm1, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izzm1, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izzm1, 3, j)

          wt = zb_quadwt_patch ( ixx, iyy, izzm1, j)

          x = ( zb_quadpts_patch ( ixx, iyy, izzm1, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izzm1, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izzm1, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 6, :) = CMat_B ( 6, :) + wt * CMat_Row_B ( :)
          DVec_B ( 6) = DVec_B ( 6) + wt
     1                              * Bzf_quadpts ( ixx, iyy, izzm1, j)

          AMat_B ( 26, :) = CMat_Row_B ( :)
          BVec_B ( 26) = Bzf_quadpts ( ixx, iyy, izzm1, j)

! ----------

          j = 2

          nx = zb_quadnorml_patch ( ixx, iyy, izzm1, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izzm1, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izzm1, 3, j)

          wt = zb_quadwt_patch ( ixx, iyy, izzm1, j)

          x = ( zb_quadpts_patch ( ixx, iyy, izzm1, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izzm1, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izzm1, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 6, :) = CMat_B ( 6, :) + wt * CMat_Row_B ( :)
          DVec_B ( 6) = DVec_B ( 6) + wt
     1                              * Bzf_quadpts ( ixx, iyy, izzm1, j)

          AMat_B ( 27, :) = CMat_Row_B ( :)
          BVec_B ( 27) = Bzf_quadpts ( ixx, iyy, izzm1, j)

! ----------

          j = 3

          nx = zb_quadnorml_patch ( ixx, iyy, izzm1, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izzm1, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izzm1, 3, j)

          wt = zb_quadwt_patch ( ixx, iyy, izzm1, j)

          x = ( zb_quadpts_patch ( ixx, iyy, izzm1, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izzm1, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izzm1, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 6, :) = CMat_B ( 6, :) + wt * CMat_Row_B ( :)
          DVec_B ( 6) = DVec_B ( 6) + wt
     1                              * Bzf_quadpts ( ixx, iyy, izzm1, j)

          AMat_B ( 28, :) = CMat_Row_B ( :)
          BVec_B ( 28) = Bzf_quadpts ( ixx, iyy, izzm1, j)

! ----------

          j = 4

          nx = zb_quadnorml_patch ( ixx, iyy, izzm1, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izzm1, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izzm1, 3, j)

          wt = zb_quadwt_patch ( ixx, iyy, izzm1, j)

          x = ( zb_quadpts_patch ( ixx, iyy, izzm1, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izzm1, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izzm1, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 6, :) = CMat_B ( 6, :) + wt * CMat_Row_B ( :)
          DVec_B ( 6) = DVec_B ( 6) + wt
     1                              * Bzf_quadpts ( ixx, iyy, izzm1, j)

          AMat_B ( 29, :) = CMat_Row_B ( :)
          BVec_B ( 29) = Bzf_quadpts ( ixx, iyy, izzm1, j)

! ----------

          j = 5

          nx = zb_quadnorml_patch ( ixx, iyy, izzm1, 1, j)
          ny = zb_quadnorml_patch ( ixx, iyy, izzm1, 2, j)
          nz = zb_quadnorml_patch ( ixx, iyy, izzm1, 3, j)

          wt = zb_quadwt_patch ( ixx, iyy, izzm1, j)

          x = ( zb_quadpts_patch ( ixx, iyy, izzm1, 1, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zb_quadpts_patch ( ixx, iyy, izzm1, 2, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zb_quadpts_patch ( ixx, iyy, izzm1, 3, j)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL CMat_Row_Coefs_P2 ( x, y, z, nx, ny, nz, Coef_0,
     1                             CMat_Row_B)

          CMat_B ( 6, :) = CMat_B ( 6, :) + wt * CMat_Row_B ( :)
          DVec_B ( 6) = DVec_B ( 6) + wt
     1                              * Bzf_quadpts ( ixx, iyy, izzm1, j)

          AMat_B ( 30, :) = CMat_Row_B ( :)
          BVec_B ( 30) = Bzf_quadpts ( ixx, iyy, izzm1, j)

! --------------------------------------------------
! Now in LSQ sense match the zone-centered modes of "V" that correspond to bfx.

          AMat_B ( 31, 1) = 1.0
          BVec_B ( 31) = V ( ixx, iyy, izz, 6, 1)

          AMat_B ( 32, 2) = 1.0
          BVec_B ( 32) = V ( ixx, iyy, izz, 6, 2)

          AMat_B ( 33, 3) = 1.0
          BVec_B ( 33) = V ( ixx, iyy, izz, 6, 3)

          AMat_B ( 34, 4) = 1.0
          BVec_B ( 34) = V ( ixx, iyy, izz, 6, 4)

          AMat_B ( 35, 5) = 1.0
          BVec_B ( 35) = V ( ixx, iyy, izz, 6, 5)

          AMat_B ( 36, 6) = 1.0
          BVec_B ( 36) = V ( ixx, iyy, izz, 6, 6)

          AMat_B ( 37, 7) = 1.0
          BVec_B ( 37) = V ( ixx, iyy, izz, 6, 7)

          AMat_B ( 38, 8) = 1.0
          BVec_B ( 38) = V ( ixx, iyy, izz, 6, 8)

          AMat_B ( 39, 9) = 1.0
          BVec_B ( 39) = V ( ixx, iyy, izz, 6, 9)

          AMat_B ( 40, 10) = 1.0
          BVec_B ( 40) = V ( ixx, iyy, izz, 6, 10)

! --------------------------------------------------
! Now in LSQ sense match the zone-centered modes of "V" that correspond to bfy.

          AMat_B ( 41, 11) = 1.0
          BVec_B ( 41) = V ( ixx, iyy, izz, 7, 1)

          AMat_B ( 42, 12) = 1.0
          BVec_B ( 42) = V ( ixx, iyy, izz, 7, 2)

          AMat_B ( 43, 13) = 1.0
          BVec_B ( 43) = V ( ixx, iyy, izz, 7, 3)

          AMat_B ( 44, 14) = 1.0
          BVec_B ( 44) = V ( ixx, iyy, izz, 7, 4)

          AMat_B ( 45, 15) = 1.0
          BVec_B ( 45) = V ( ixx, iyy, izz, 7, 5)

          AMat_B ( 46, 16) = 1.0
          BVec_B ( 46) = V ( ixx, iyy, izz, 7, 6)

          AMat_B ( 47, 17) = 1.0
          BVec_B ( 47) = V ( ixx, iyy, izz, 7, 7)

          AMat_B ( 48, 18) = 1.0
          BVec_B ( 48) = V ( ixx, iyy, izz, 7, 8)

          AMat_B ( 49, 19) = 1.0
          BVec_B ( 49) = V ( ixx, iyy, izz, 7, 9)

          AMat_B ( 50, 20) = 1.0
          BVec_B ( 50) = V ( ixx, iyy, izz, 7, 10)

! --------------------------------------------------
! Now in LSQ sense match the zone-centered modes of "V" that correspond to bfz.

          AMat_B ( 51, 21) = 1.0
          BVec_B ( 51) = V ( ixx, iyy, izz, 8, 1)

          AMat_B ( 52, 22) = 1.0
          BVec_B ( 52) = V ( ixx, iyy, izz, 8, 2)

          AMat_B ( 53, 23) = 1.0
          BVec_B ( 53) = V ( ixx, iyy, izz, 8, 3)

          AMat_B ( 54, 24) = 1.0
          BVec_B ( 54) = V ( ixx, iyy, izz, 8, 4)

          AMat_B ( 55, 25) = 1.0
          BVec_B ( 55) = V ( ixx, iyy, izz, 8, 5)

          AMat_B ( 56, 26) = 1.0
          BVec_B ( 56) = V ( ixx, iyy, izz, 8, 6)

          AMat_B ( 57, 27) = 1.0
          BVec_B ( 57) = V ( ixx, iyy, izz, 8, 7)

          AMat_B ( 58, 28) = 1.0
          BVec_B ( 58) = V ( ixx, iyy, izz, 8, 8)

          AMat_B ( 59, 29) = 1.0
          BVec_B ( 59) = V ( ixx, iyy, izz, 8, 9)

          AMat_B ( 60, 30) = 1.0
          BVec_B ( 60) = V ( ixx, iyy, izz, 8, 10)

! --------------------------------------------------
! Now in LSQ sense match the divergence-free constraint.

          AMat_B ( 61, 2) = 1.0  ;  AMat_B ( 61, 13) = 1.0
          AMat_B ( 61, 24) = 1.0
          BVec_B ( 61) = 0.0

          AMat_B ( 62, 5) = 2.0  ;  AMat_B ( 62, 18) = 1.0
          AMat_B ( 62, 30) = 1.0
          BVec_B ( 62) = 0.0

          AMat_B ( 63, 8) = 1.0  ;  AMat_B ( 63, 16) = 2.0
          AMat_B ( 63, 29) = 1.0
          BVec_B ( 63) = 0.0

          AMat_B ( 64, 10) = 1.0  ;  AMat_B ( 64, 19) = 1.0
          AMat_B ( 64, 27) = 2.0
          BVec_B ( 64) = 0.0

! ----------------------------------------
! Solve a CLSQ system with "NConstraint = 6, NUnknown = 30, NLSQ = 64".

! Build transpose of "AMat_B", "ATrans_B ( 30, 64)". From that, build
! "ATransA_B ( 30, 30)" and "ATransB_B ( 30)". Also build
! "CTrans_B ( 30, 6)".

          ATrans_B = TRANSPOSE ( AMat_B)

          ATransA_B = MATMUL ( ATrans_B, AMat_B)

          ATransB_B = MATMUL ( ATrans_B, BVec_B)

          CTrans_B = TRANSPOSE ( CMat_B)

! ----------------------------------------
! Assemble the matrix "CLSQMat_B ( 30 + 6, 30 + 6)"
! and its right hand side "CLSQVec_B ( 30 + 6)".

          CLSQMat_B ( 1: 30, 1: 30) = ATransA_B ( 1: 30, 1: 30)

          CLSQMat_B ( 1: 30, 30 + 1: 30 + 6) = CTrans_B ( 1: 30, 1: 6)

          CLSQMat_B ( 30 + 1: 30 + 6, 1: 30) = CMat_B ( 1: 6, 1: 30)

          CLSQMat_B ( 30 + 1: 30 + 6, 30 + 1: 30 + 6) = 0.0

! ----------

          CLSQVec_B ( 1: 30) = ATransB_B ( 1: 30)

          CLSQVec_B ( 30 + 1: 30 + 6) = DVec_B ( 1: 6)

! ----------------------------------------
! Solve the system "CLSQMat_B . (x, z)^T = CLSQVec_B". Extract from that the vector of
! unknowns.

          CALL MATSOLVE ( 30 + 6, CLSQMat_B, CLSQVec_B,
     1                    indx_CLSQ_B, vv_CLSQ_B)

! ----------------------------------------
! Now hand back the solution to the zone-centered reconstruction of the magnetic
! field that is closer to being divergence-free. The conserved variable is also
! reinitialized.

          aa = ( U ( ixx, iyy, izz, 6, 1)**2
     1         + U ( ixx, iyy, izz, 7, 1)**2
     1         + U ( ixx, iyy, izz, 8, 1)**2) * ob_8pi

          B1_zc ( ixx, iyy, izz, 1: 10) = CLSQVec_B ( 1: 10)
          V ( ixx, iyy, izz, 6, 1: 10) = CLSQVec_B ( 1: 10)
          U ( ixx, iyy, izz, 6, 1) = V ( ixx, iyy, izz, 6, 1)

          B2_zc ( ixx, iyy, izz, 1: 10) = CLSQVec_B ( 11: 20)
          V ( ixx, iyy, izz, 7, 1: 10) = CLSQVec_B ( 11: 20)
          U ( ixx, iyy, izz, 7, 1) = V ( ixx, iyy, izz, 7, 1)

          B3_zc ( ixx, iyy, izz, 1: 10) = CLSQVec_B ( 21: 30)
          V ( ixx, iyy, izz, 8, 1: 10) = CLSQVec_B ( 21: 30)
          U ( ixx, iyy, izz, 8, 1) = V ( ixx, iyy, izz, 8, 1)

#if ( 1 == 1)
          bb = ( U ( ixx, iyy, izz, 6, 1)**2
     1         + U ( ixx, iyy, izz, 7, 1)**2
     1         + U ( ixx, iyy, izz, 8, 1)**2) * ob_8pi

          U ( ixx, iyy, izz, 5, 1) = U ( ixx, iyy, izz, 5, 1) - aa + bb
#endif

! --------------------------------------------------

          END DO  ! End of "DO ixx = - 2, ntp_zones + 3" loop.

        END DO  ! End of "DO iyy = - 2, ntp_zones + 3" loop.

      END DO  ! End of "DO izz = - 2, nr_zones + 3" loop.

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 9) We limit the zone-centered Electric fields by using WENO limiting.

! ----------------------------------------------------------------------
! STEP 9.a) Build "Ex_zc, Ey_zc, Ez_zc".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, B0, x, y, z, bflag,          &
!$OMP&    time_use, a11, a12, a13, vx, vy, vz, Bx, By, Bz)              &
!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    B1_zc, B2_zc, B3_zc, Coef_Basis_Moments,                      &
!$OMP&    V, Ex_zc, Ey_zc, Ez_zc, v_centroid_patch)                     &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 3)

!$acc parallel vector_length(64) present( V, Ex_zc, Ey_zc, Ez_zc,       &
!$acc&    B1_zc, B2_zc, B3_zc, Coef_Basis_Moments,                      &
!$acc&    v_centroid_patch)
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, B0, x, y, z, bflag, time_use,        &
!$acc&    a11, a12, a13, vx, vy, vz, Bx, By, Bz)

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

! ------------------------------

        B0 ( :) = 0.0

#if( USE_B0B1_FORMULATION == 1)
        x = v_centroid_patch ( ixx, iyy, izz, 1)
        y = v_centroid_patch ( ixx, iyy, izz, 2)
        z = v_centroid_patch ( ixx, iyy, izz, 3)

        bflag = 1         ! 1 == magnetic field;  0 == vector potential.
        time_use = 0.0    ! useful for the case of active rotation.

        CALL GET_MAGNETIC_FIELD ( bflag, time_use, x, y, z,
     1                              a11, a12, a13)

        B0 ( 1) = a11
        B0 ( 2) = a12
        B0 ( 3) = a13
#endif

! ------------------------------

! --------------------
! Make new Bx, By, Bz with third order reconstruction.

        Bx = B0 ( 1) + V ( ixx, iyy, izz, 6, 1)
     1 - V ( ixx, iyy, izz, 6, 5)*Coef_Basis_Moments ( ixx, iyy, izz, 5)
     1 - V ( ixx, iyy, izz, 6, 6)*Coef_Basis_Moments ( ixx, iyy, izz, 6)
     1 - V ( ixx, iyy, izz, 6, 7)*Coef_Basis_Moments ( ixx, iyy, izz, 7)
     1 - V ( ixx, iyy, izz, 6, 8)*Coef_Basis_Moments ( ixx, iyy, izz, 8)
     1 - V ( ixx, iyy, izz, 6, 9)*Coef_Basis_Moments ( ixx, iyy, izz, 9)
     1 - V ( ixx, iyy, izz, 6,10)*Coef_Basis_Moments ( ixx, iyy, izz,10)

        By = B0 ( 2) + V ( ixx, iyy, izz, 7, 1)
     1 - V ( ixx, iyy, izz, 7, 5)*Coef_Basis_Moments ( ixx, iyy, izz, 5)
     1 - V ( ixx, iyy, izz, 7, 6)*Coef_Basis_Moments ( ixx, iyy, izz, 6)
     1 - V ( ixx, iyy, izz, 7, 7)*Coef_Basis_Moments ( ixx, iyy, izz, 7)
     1 - V ( ixx, iyy, izz, 7, 8)*Coef_Basis_Moments ( ixx, iyy, izz, 8)
     1 - V ( ixx, iyy, izz, 7, 9)*Coef_Basis_Moments ( ixx, iyy, izz, 9)
     1 - V ( ixx, iyy, izz, 7,10)*Coef_Basis_Moments ( ixx, iyy, izz,10)

        Bz = B0 ( 3) + V ( ixx, iyy, izz, 8, 1)
     1 - V ( ixx, iyy, izz, 8, 5)*Coef_Basis_Moments ( ixx, iyy, izz, 5)
     1 - V ( ixx, iyy, izz, 8, 6)*Coef_Basis_Moments ( ixx, iyy, izz, 6)
     1 - V ( ixx, iyy, izz, 8, 7)*Coef_Basis_Moments ( ixx, iyy, izz, 7)
     1 - V ( ixx, iyy, izz, 8, 8)*Coef_Basis_Moments ( ixx, iyy, izz, 8)
     1 - V ( ixx, iyy, izz, 8, 9)*Coef_Basis_Moments ( ixx, iyy, izz, 9)
     1 - V ( ixx, iyy, izz, 8,10)*Coef_Basis_Moments ( ixx, iyy, izz,10)

! --------------------
! Make new vx, vy, vz with third order reconstruction.

        vx = V ( ixx, iyy, izz, 2, 1)
     1 - V ( ixx, iyy, izz, 2, 5)*Coef_Basis_Moments ( ixx, iyy, izz, 5)
     1 - V ( ixx, iyy, izz, 2, 6)*Coef_Basis_Moments ( ixx, iyy, izz, 6)
     1 - V ( ixx, iyy, izz, 2, 7)*Coef_Basis_Moments ( ixx, iyy, izz, 7)
     1 - V ( ixx, iyy, izz, 2, 8)*Coef_Basis_Moments ( ixx, iyy, izz, 8)
     1 - V ( ixx, iyy, izz, 2, 9)*Coef_Basis_Moments ( ixx, iyy, izz, 9)
     1 - V ( ixx, iyy, izz, 2,10)*Coef_Basis_Moments ( ixx, iyy, izz,10)

        vy = V ( ixx, iyy, izz, 3, 1)
     1 - V ( ixx, iyy, izz, 3, 5)*Coef_Basis_Moments ( ixx, iyy, izz, 5)
     1 - V ( ixx, iyy, izz, 3, 6)*Coef_Basis_Moments ( ixx, iyy, izz, 6)
     1 - V ( ixx, iyy, izz, 3, 7)*Coef_Basis_Moments ( ixx, iyy, izz, 7)
     1 - V ( ixx, iyy, izz, 3, 8)*Coef_Basis_Moments ( ixx, iyy, izz, 8)
     1 - V ( ixx, iyy, izz, 3, 9)*Coef_Basis_Moments ( ixx, iyy, izz, 9)
     1 - V ( ixx, iyy, izz, 3,10)*Coef_Basis_Moments ( ixx, iyy, izz,10)

        vz = V ( ixx, iyy, izz, 4, 1)
     1 - V ( ixx, iyy, izz, 4, 5)*Coef_Basis_Moments ( ixx, iyy, izz, 5)
     1 - V ( ixx, iyy, izz, 4, 6)*Coef_Basis_Moments ( ixx, iyy, izz, 6)
     1 - V ( ixx, iyy, izz, 4, 7)*Coef_Basis_Moments ( ixx, iyy, izz, 7)
     1 - V ( ixx, iyy, izz, 4, 8)*Coef_Basis_Moments ( ixx, iyy, izz, 8)
     1 - V ( ixx, iyy, izz, 4, 9)*Coef_Basis_Moments ( ixx, iyy, izz, 9)
     1 - V ( ixx, iyy, izz, 4,10)*Coef_Basis_Moments ( ixx, iyy, izz,10)

! --------------------

      Ex_zc ( ixx, iyy, izz, 1) = - ( vy * Bz - vz * By)
      Ey_zc ( ixx, iyy, izz, 1) = - ( vz * Bx - vx * Bz)
      Ez_zc ( ixx, iyy, izz, 1) = - ( vx * By - vy * Bx)

! ------------------------------

      END DO
      END DO
      END DO

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! STEP 9.b) WENO limiting of "Ex_zc, Ey_zc, Ez_zc".
! This is a second order, but truly multid, WENO
! reconstruction. We have been able to simplify the process to the point where
! we only have to invert eight small 3X3 matrices for the eight one-sided
! stencils. The ninth, centered stencil is obtained by arithmetic averaging
! the linear moments from the previous eight stencils and is given a much higher
! weight. We only do the reconstruction for zones that are two zones out from
! the physical zones.

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1, icomp,  &
!$OMP&    ixxp1, iyyp1, izzp1, obl0, i, j, nn,                          &
!$OMP&    x, y, z,                                                      &
!$OMP&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,               &
!$OMP&    ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,       &
!$OMP&    a11, a12, a13, a21, a22, a23, a31, a32, a33, obDet,           &
!$OMP&    R38, X39, AInv338, IS9, Wt9)                                  &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,&
!$OMP&    pow, central_weight, tinynum, smallnum,                       &
!$OMP&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$OMP&    vert_coords_patch, Coef_Basis_Moments,                        &
!$OMP&    istencil_patch, flattengr, Ex_zc, Ey_zc, Ez_zc)               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    pow, central_weight, tinynum, smallnum,                       &
!$acc&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$acc&    vert_coords_patch, Coef_Basis_Moments,                        &
!$acc&    istencil_patch, flattengr, Ex_zc, Ey_zc, Ez_zc)

!$acc loop gang collapse(2) independent                                 &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1, icomp,          &
!$acc&    ixxp1, iyyp1, izzp1, obl0, i, j, nn,                          &
!$acc&    x, y, z,                                                      &
!$acc&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,               &
!$acc&    ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,       &
!$acc&    a11, a12, a13, a21, a22, a23, a31, a32, a33, obDet,           &
!$acc&    AInv338, R38, X39, IS9, Wt9)

! ----------------------------------------------------------------------

      DO izz = - 2, nr_zones + 3

        DO iyy = - 2, ntp_zones + 3

          DO ixx = - 2, ntp_zones + 3

          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN

! ----------------------------------------------------------------------

      izzm1 = izz - 1
      izzp1 = izz + 1

      iyym1 = iyy - 1
      iyyp1 = iyy + 1

      ixxm1 = ixx - 1
      ixxp1 = ixx + 1

      ixxr = istencil_patch ( ixx, iyy, 1, 1)  ! For indexing zone one Right.
      iyyr = istencil_patch ( ixx, iyy, 2, 1)

      ixxl = istencil_patch ( ixx, iyy, 1, 2)  ! For indexing zone one Left.
      iyyl = istencil_patch ( ixx, iyy, 2, 2)

      ixxu = istencil_patch ( ixx, iyy, 1, 3)  ! For indexing zone one Up.
      iyyu = istencil_patch ( ixx, iyy, 2, 3)

      ixxd = istencil_patch ( ixx, iyy, 1, 4)  ! For indexing zone one Down.
      iyyd = istencil_patch ( ixx, iyy, 2, 4)

      ixxru = istencil_patch ( ixx, iyy, 1, 9)  ! For indexing zone that is Right-Up.
      iyyru = istencil_patch ( ixx, iyy, 2, 9)

      ixxlu = istencil_patch ( ixx, iyy, 1, 10)  ! For indexing zone that is Left-Up.
      iyylu = istencil_patch ( ixx, iyy, 2, 10)

      ixxld = istencil_patch ( ixx, iyy, 1, 11)  ! For indexing zone that is Left-Down.
      iyyld = istencil_patch ( ixx, iyy, 2, 11)

      ixxrd = istencil_patch ( ixx, iyy, 1, 12)  ! For indexing zone that is Right-Down.
      iyyrd = istencil_patch ( ixx, iyy, 2, 12)


! ----------------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

! --------------------------------------------------
! We start on the expensive part associated with building "AInv338".
! --------------------------------------------------

! ------------------------------
! Stencil 1 in direction "( Rt, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

          a11 = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a12 = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a13 = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a21 = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a22 = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a23 = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 1 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 1  ! Stencil 1

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 2 in direction "( Rt, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

!          a11 = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a12 = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a13 = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


!          a21 = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a22 = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a23 = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzm1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzm1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzm1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 2 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 2  ! Stencil 2

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 3 in direction "( Lf, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

          a11 = ( v_centroid_patch ( ixxl, iyyl, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a12 = ( v_centroid_patch ( ixxl, iyyl, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a13 = ( v_centroid_patch ( ixxl, iyyl, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a21 = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a22 = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a23 = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 3 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 3  ! Stencil 3

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil in direction 4 "( Lf, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

!          a11 = ( v_centroid_patch ( ixxl, iyyl, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a12 = ( v_centroid_patch ( ixxl, iyyl, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a13 = ( v_centroid_patch ( ixxl, iyyl, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


!          a21 = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a22 = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a23 = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzm1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzm1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzm1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 4 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 4  ! Stencil 4

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 5 in direction "( Rt, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

          a11 = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a12 = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a13 = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a21 = ( v_centroid_patch ( ixxd, iyyd, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a22 = ( v_centroid_patch ( ixxd, iyyd, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a23 = ( v_centroid_patch ( ixxd, iyyd, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 5 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 5  ! Stencil 5

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 6 in direction "( Rt, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

!          a11 = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a12 = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a13 = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


!          a21 = ( v_centroid_patch ( ixxd, iyyd, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a22 = ( v_centroid_patch ( ixxd, iyyd, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a23 = ( v_centroid_patch ( ixxd, iyyd, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzm1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzm1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzm1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 6 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 6  ! Stencil 6

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 7 in direction "( Lf, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

          a11 = ( v_centroid_patch ( ixxl, iyyl, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a12 = ( v_centroid_patch ( ixxl, iyyl, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a13 = ( v_centroid_patch ( ixxl, iyyl, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a21 = ( v_centroid_patch ( ixxd, iyyd, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a22 = ( v_centroid_patch ( ixxd, iyyd, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a23 = ( v_centroid_patch ( ixxd, iyyd, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 7 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 7  ! Stencil 7

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! ------------------------------
! Stencil 8 in direction "( Lf, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

! Build 3 sets of scaled relative coordinates from 3 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 3X3 matrix.

!          a11 = ( v_centroid_patch ( ixxl, iyyl, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a12 = ( v_centroid_patch ( ixxl, iyyl, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a13 = ( v_centroid_patch ( ixxl, iyyl, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


!          a21 = ( v_centroid_patch ( ixxd, iyyd, izz, 1)
!     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

!          a22 = ( v_centroid_patch ( ixxd, iyyd, izz, 2)
!     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

!          a23 = ( v_centroid_patch ( ixxd, iyyd, izz, 3)
!     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0


          a31 = ( v_centroid_patch ( ixx, iyy, izzm1, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          a32 = ( v_centroid_patch ( ixx, iyy, izzm1, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          a33 = ( v_centroid_patch ( ixx, iyy, izzm1, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Invert the 3X3 matrix by initially building "obDet"; then store the matrix in the
! j = 8 matrix location in "AInv338".

          obDet = (-a13 * a22 * a31 + a12 * a23 * a31 + a13 * a21 * a32
     1            - a11 * a23 * a32 - a12 * a21 * a33 + a11 * a22 * a33)

          obDet = SIGN ( 1.0, obDet) / AMAX1 ( tinynum, ABS ( obDet) )

          j = 8  ! Stencil 8

          AInv338 ( 1, 1, j) = ( - a23 * a32 + a22 *a33 ) * obDet

          AInv338 ( 1, 2, j) = ( a13 * a32 - a12 * a33 ) * obDet

          AInv338 ( 1, 3, j) = ( - a13 * a22 + a12 * a23 ) * obDet


          AInv338 ( 2, 1, j) = ( a23 * a31 - a21 * a33 ) * obDet

          AInv338 ( 2, 2, j) = ( - a13 * a31 + a11 * a33 ) * obDet

          AInv338 ( 2, 3, j) = ( a13 * a21 - a11 * a23 ) * obDet


          AInv338 ( 3, 1, j) = ( - a22 * a31 + a21 * a32 ) * obDet

          AInv338 ( 3, 2, j) = ( a12 * a31 - a11 * a32 ) * obDet

          AInv338 ( 3, 3, j) = ( - a12 * a21 + a11 * a22 ) * obDet

! --------------------------------------------------
! By this point, the expensive part associated with building "AInv338" is done.
! --------------------------------------------------
! Interpolation of "Ex_zc" done using the matrices that have been constructed above.
! --------------------------------------------------
! ------------------------------
! Stencil 1 in direction "( Rt, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 1

          R38 ( 1, j) = Ex_zc ( ixxr, iyyr, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ex_zc ( ixxu, iyyu, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ex_zc ( ixx, iyy, izzp1, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 2 in direction "( Rt, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 2

          R38 ( 1, j) = Ex_zc ( ixxr, iyyr, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ex_zc ( ixxu, iyyu, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ex_zc ( ixx, iyy, izzm1, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 3 in direction "( Lf, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 3

          R38 ( 1, j) = Ex_zc ( ixxl, iyyl, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ex_zc ( ixxu, iyyu, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ex_zc ( ixx, iyy, izzp1, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil in direction 4 "( Lf, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 4

          R38 ( 1, j) = Ex_zc ( ixxl, iyyl, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ex_zc ( ixxu, iyyu, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ex_zc ( ixx, iyy, izzm1, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 5 in direction "( Rt, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 5

          R38 ( 1, j) = Ex_zc ( ixxr, iyyr, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ex_zc ( ixxd, iyyd, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ex_zc ( ixx, iyy, izzp1, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 6 in direction "( Rt, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 6

          R38 ( 1, j) = Ex_zc ( ixxr, iyyr, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ex_zc ( ixxd, iyyd, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ex_zc ( ixx, iyy, izzm1, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 7 in direction "( Lf, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 7

          R38 ( 1, j) = Ex_zc ( ixxl, iyyl, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ex_zc ( ixxd, iyyd, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ex_zc ( ixx, iyy, izzp1, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 8 in direction "( Lf, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 8

          R38 ( 1, j) = Ex_zc ( ixxl, iyyl, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ex_zc ( ixxd, iyyd, izz, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ex_zc ( ixx, iyy, izzm1, 1)
     1                     - Ex_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! By this point, all three gradients from the first 8 stencils have been
! evaluated and stored in "X39 ( 1:3, 1:8)". Now obtain the
! centered, ninth stencil's gradients by arithmetic averaging.
! ------------------------------

          X39 ( :, 9) = 0.125 * ( X39 ( :, 1)
     1  + X39 ( :, 2) + X39 ( :, 3) + X39 ( :, 4)
     1  + X39 ( :, 5) + X39 ( :, 6) + X39 ( :, 7)
     1  + X39 ( :, 8) )

! ------------------------------
! Build all 9 smoothness indicators for all 9 stencils.
! ------------------------------

          IS9 ( 1: 9) = X39 ( 1, 1: 9)**2
     1  + X39 ( 2, 1: 9)**2 + X39 ( 3, 1: 9)**2

! ------------------------------
! Initialize 9 linear weights and use them to obtain non-linear weights
! for all 9 stencils. Note that here the ninth stencil is the central one.
! ------------------------------

          Wt9 ( 1: 8) = 1.0
          Wt9 ( 9) = central_weight

          Wt9 ( 1: 9) = Wt9 ( 1: 9) / ( IS9 ( 1: 9)
     1                                          + smallnum)**pow

          a11 = 1.0 / ( Wt9 ( 1) + Wt9 ( 2)
     1      + Wt9 ( 3) + Wt9 ( 4) + Wt9 ( 5)
     1      + Wt9 ( 6) + Wt9 ( 7) + Wt9 ( 8)
     1      + Wt9 ( 9) )

          Wt9 ( 1: 9) = Wt9 ( 1: 9) * a11

! ------------------------------
! Make all the three non-linearly hybridized gradients of
! "Ex_zc ( ixx, iyy, izz, 2: 4)" for this "Ex_zc" variable.
! ------------------------------

          Ex_zc ( ixx, iyy, izz, 2) =
     1        Wt9 ( 1) * X39 ( 1, 1)
     1      + Wt9 ( 2) * X39 ( 1, 2)
     1      + Wt9 ( 3) * X39 ( 1, 3)
     1      + Wt9 ( 4) * X39 ( 1, 4)
     1      + Wt9 ( 5) * X39 ( 1, 5)
     1      + Wt9 ( 6) * X39 ( 1, 6)
     1      + Wt9 ( 7) * X39 ( 1, 7)
     1      + Wt9 ( 8) * X39 ( 1, 8)
     1      + Wt9 ( 9) * X39 ( 1, 9)

          Ex_zc ( ixx, iyy, izz, 3) =
     1        Wt9 ( 1) * X39 ( 2, 1)
     1      + Wt9 ( 2) * X39 ( 2, 2)
     1      + Wt9 ( 3) * X39 ( 2, 3)
     1      + Wt9 ( 4) * X39 ( 2, 4)
     1      + Wt9 ( 5) * X39 ( 2, 5)
     1      + Wt9 ( 6) * X39 ( 2, 6)
     1      + Wt9 ( 7) * X39 ( 2, 7)
     1      + Wt9 ( 8) * X39 ( 2, 8)
     1      + Wt9 ( 9) * X39 ( 2, 9)

          Ex_zc ( ixx, iyy, izz, 4) =
     1        Wt9 ( 1) * X39 ( 3, 1)
     1      + Wt9 ( 2) * X39 ( 3, 2)
     1      + Wt9 ( 3) * X39 ( 3, 3)
     1      + Wt9 ( 4) * X39 ( 3, 4)
     1      + Wt9 ( 5) * X39 ( 3, 5)
     1      + Wt9 ( 6) * X39 ( 3, 6)
     1      + Wt9 ( 7) * X39 ( 3, 7)
     1      + Wt9 ( 8) * X39 ( 3, 8)
     1      + Wt9 ( 9) * X39 ( 3, 9)

! --------------------------------------------------
! Interpolation of "Ey_zc" done using the matrices that have been constructed above.
! --------------------------------------------------
! ------------------------------
! Stencil 1 in direction "( Rt, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 1

          R38 ( 1, j) = Ey_zc ( ixxr, iyyr, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ey_zc ( ixxu, iyyu, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ey_zc ( ixx, iyy, izzp1, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 2 in direction "( Rt, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 2

          R38 ( 1, j) = Ey_zc ( ixxr, iyyr, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ey_zc ( ixxu, iyyu, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ey_zc ( ixx, iyy, izzm1, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 3 in direction "( Lf, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 3

          R38 ( 1, j) = Ey_zc ( ixxl, iyyl, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ey_zc ( ixxu, iyyu, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ey_zc ( ixx, iyy, izzp1, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil in direction 4 "( Lf, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 4

          R38 ( 1, j) = Ey_zc ( ixxl, iyyl, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ey_zc ( ixxu, iyyu, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ey_zc ( ixx, iyy, izzm1, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 5 in direction "( Rt, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 5

          R38 ( 1, j) = Ey_zc ( ixxr, iyyr, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ey_zc ( ixxd, iyyd, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ey_zc ( ixx, iyy, izzp1, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 6 in direction "( Rt, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 6

          R38 ( 1, j) = Ey_zc ( ixxr, iyyr, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ey_zc ( ixxd, iyyd, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ey_zc ( ixx, iyy, izzm1, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 7 in direction "( Lf, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 7

          R38 ( 1, j) = Ey_zc ( ixxl, iyyl, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ey_zc ( ixxd, iyyd, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ey_zc ( ixx, iyy, izzp1, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 8 in direction "( Lf, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 8

          R38 ( 1, j) = Ey_zc ( ixxl, iyyl, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ey_zc ( ixxd, iyyd, izz, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ey_zc ( ixx, iyy, izzm1, 1)
     1                     - Ey_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! By this point, all three gradients from the first 8 stencils have been
! evaluated and stored in "X39 ( 1:3, 1:8)". Now obtain the
! centered, ninth stencil's gradients by arithmetic averaging.
! ------------------------------

          X39 ( :, 9) = 0.125 * ( X39 ( :, 1)
     1  + X39 ( :, 2) + X39 ( :, 3) + X39 ( :, 4)
     1  + X39 ( :, 5) + X39 ( :, 6) + X39 ( :, 7)
     1  + X39 ( :, 8) )

! ------------------------------
! Build all 9 smoothness indicators for all 9 stencils.
! ------------------------------

          IS9 ( 1: 9) = X39 ( 1, 1: 9)**2
     1  + X39 ( 2, 1: 9)**2 + X39 ( 3, 1: 9)**2

! ------------------------------
! Initialize 9 linear weights and use them to obtain non-linear weights
! for all 9 stencils. Note that here the ninth stencil is the central one.
! ------------------------------

          Wt9 ( 1: 8) = 1.0
          Wt9 ( 9) = central_weight

          Wt9 ( 1: 9) = Wt9 ( 1: 9) / ( IS9 ( 1: 9)
     1                                          + smallnum)**pow

          a11 = 1.0 / ( Wt9 ( 1) + Wt9 ( 2)
     1      + Wt9 ( 3) + Wt9 ( 4) + Wt9 ( 5)
     1      + Wt9 ( 6) + Wt9 ( 7) + Wt9 ( 8)
     1      + Wt9 ( 9) )

          Wt9 ( 1: 9) = Wt9 ( 1: 9) * a11

! ------------------------------
! Make all the three non-linearly hybridized gradients of
! "Ey_zc ( ixx, iyy, izz, 2: 4)" for this "Ey_zc" variable.
! ------------------------------

          Ey_zc ( ixx, iyy, izz, 2) =
     1        Wt9 ( 1) * X39 ( 1, 1)
     1      + Wt9 ( 2) * X39 ( 1, 2)
     1      + Wt9 ( 3) * X39 ( 1, 3)
     1      + Wt9 ( 4) * X39 ( 1, 4)
     1      + Wt9 ( 5) * X39 ( 1, 5)
     1      + Wt9 ( 6) * X39 ( 1, 6)
     1      + Wt9 ( 7) * X39 ( 1, 7)
     1      + Wt9 ( 8) * X39 ( 1, 8)
     1      + Wt9 ( 9) * X39 ( 1, 9)

          Ey_zc ( ixx, iyy, izz, 3) =
     1        Wt9 ( 1) * X39 ( 2, 1)
     1      + Wt9 ( 2) * X39 ( 2, 2)
     1      + Wt9 ( 3) * X39 ( 2, 3)
     1      + Wt9 ( 4) * X39 ( 2, 4)
     1      + Wt9 ( 5) * X39 ( 2, 5)
     1      + Wt9 ( 6) * X39 ( 2, 6)
     1      + Wt9 ( 7) * X39 ( 2, 7)
     1      + Wt9 ( 8) * X39 ( 2, 8)
     1      + Wt9 ( 9) * X39 ( 2, 9)

          Ey_zc ( ixx, iyy, izz, 4) =
     1        Wt9 ( 1) * X39 ( 3, 1)
     1      + Wt9 ( 2) * X39 ( 3, 2)
     1      + Wt9 ( 3) * X39 ( 3, 3)
     1      + Wt9 ( 4) * X39 ( 3, 4)
     1      + Wt9 ( 5) * X39 ( 3, 5)
     1      + Wt9 ( 6) * X39 ( 3, 6)
     1      + Wt9 ( 7) * X39 ( 3, 7)
     1      + Wt9 ( 8) * X39 ( 3, 8)
     1      + Wt9 ( 9) * X39 ( 3, 9)

! --------------------------------------------------
! Interpolation of "Ez_zc" done using the matrices that have been constructed above.
! --------------------------------------------------
! ------------------------------
! Stencil 1 in direction "( Rt, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 1

          R38 ( 1, j) = Ez_zc ( ixxr, iyyr, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ez_zc ( ixxu, iyyu, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ez_zc ( ixx, iyy, izzp1, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 2 in direction "( Rt, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 2

          R38 ( 1, j) = Ez_zc ( ixxr, iyyr, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ez_zc ( ixxu, iyyu, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ez_zc ( ixx, iyy, izzm1, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 3 in direction "( Lf, Up, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 3

          R38 ( 1, j) = Ez_zc ( ixxl, iyyl, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ez_zc ( ixxu, iyyu, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ez_zc ( ixx, iyy, izzp1, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil in direction 4 "( Lf, Up, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 4

          R38 ( 1, j) = Ez_zc ( ixxl, iyyl, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ez_zc ( ixxu, iyyu, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ez_zc ( ixx, iyy, izzm1, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 5 in direction "( Rt, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 5

          R38 ( 1, j) = Ez_zc ( ixxr, iyyr, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ez_zc ( ixxd, iyyd, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ez_zc ( ixx, iyy, izzp1, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 6 in direction "( Rt, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxr, iyyr, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 6

          R38 ( 1, j) = Ez_zc ( ixxr, iyyr, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ez_zc ( ixxd, iyyd, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ez_zc ( ixx, iyy, izzm1, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 7 in direction "( Lf, Dn, izzp1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzp1) }".
! ------------------------------

          j = 7

          R38 ( 1, j) = Ez_zc ( ixxl, iyyl, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ez_zc ( ixxd, iyyd, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ez_zc ( ixx, iyy, izzp1, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! Stencil 8 in direction "( Lf, Dn, izzm1)" which uses the three zone stencil
! "{ ( ixxl, iyyl, izz), ( ixxd, iyyd, izz), ( ixx, iyy, izzm1) }".
! ------------------------------

          j = 8

          R38 ( 1, j) = Ez_zc ( ixxl, iyyl, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 2, j) = Ez_zc ( ixxd, iyyd, izz, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

          R38 ( 3, j) = Ez_zc ( ixx, iyy, izzm1, 1)
     1                     - Ez_zc ( ixx, iyy, izz, 1)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X39 ( :, j) = MATMUL ( AInv338 ( :, :, j),
     1                                R38 ( :, j) )
#else
!$acc loop seq
          DO i = 1, 3
            x = 0.0
            DO nn = 1, 3
              x = x + AInv338 ( i, nn, j) * R38 ( nn, j)
            END DO
            X39 ( i, j) = x
          END DO
#endif

! ------------------------------
! By this point, all three gradients from the first 8 stencils have been
! evaluated and stored in "X39 ( 1:3, 1:8)". Now obtain the
! centered, ninth stencil's gradients by arithmetic averaging.
! ------------------------------

          X39 ( :, 9) = 0.125 * ( X39 ( :, 1)
     1  + X39 ( :, 2) + X39 ( :, 3) + X39 ( :, 4)
     1  + X39 ( :, 5) + X39 ( :, 6) + X39 ( :, 7)
     1  + X39 ( :, 8) )

! ------------------------------
! Build all 9 smoothness indicators for all 9 stencils.
! ------------------------------

          IS9 ( 1: 9) = X39 ( 1, 1: 9)**2
     1  + X39 ( 2, 1: 9)**2 + X39 ( 3, 1: 9)**2

! ------------------------------
! Initialize 9 linear weights and use them to obtain non-linear weights
! for all 9 stencils. Note that here the ninth stencil is the central one.
! ------------------------------

          Wt9 ( 1: 8) = 1.0
          Wt9 ( 9) = central_weight

          Wt9 ( 1: 9) = Wt9 ( 1: 9) / ( IS9 ( 1: 9)
     1                                          + smallnum)**pow

          a11 = 1.0 / ( Wt9 ( 1) + Wt9 ( 2)
     1      + Wt9 ( 3) + Wt9 ( 4) + Wt9 ( 5)
     1      + Wt9 ( 6) + Wt9 ( 7) + Wt9 ( 8)
     1      + Wt9 ( 9) )

          Wt9 ( 1: 9) = Wt9 ( 1: 9) * a11

! ------------------------------
! Make all the three non-linearly hybridized gradients of
! "Ez_zc ( ixx, iyy, izz, 2: 4)" for this "Ez_zc" variable.
! ------------------------------

          Ez_zc ( ixx, iyy, izz, 2) =
     1        Wt9 ( 1) * X39 ( 1, 1)
     1      + Wt9 ( 2) * X39 ( 1, 2)
     1      + Wt9 ( 3) * X39 ( 1, 3)
     1      + Wt9 ( 4) * X39 ( 1, 4)
     1      + Wt9 ( 5) * X39 ( 1, 5)
     1      + Wt9 ( 6) * X39 ( 1, 6)
     1      + Wt9 ( 7) * X39 ( 1, 7)
     1      + Wt9 ( 8) * X39 ( 1, 8)
     1      + Wt9 ( 9) * X39 ( 1, 9)

          Ez_zc ( ixx, iyy, izz, 3) =
     1        Wt9 ( 1) * X39 ( 2, 1)
     1      + Wt9 ( 2) * X39 ( 2, 2)
     1      + Wt9 ( 3) * X39 ( 2, 3)
     1      + Wt9 ( 4) * X39 ( 2, 4)
     1      + Wt9 ( 5) * X39 ( 2, 5)
     1      + Wt9 ( 6) * X39 ( 2, 6)
     1      + Wt9 ( 7) * X39 ( 2, 7)
     1      + Wt9 ( 8) * X39 ( 2, 8)
     1      + Wt9 ( 9) * X39 ( 2, 9)

          Ez_zc ( ixx, iyy, izz, 4) =
     1        Wt9 ( 1) * X39 ( 3, 1)
     1      + Wt9 ( 2) * X39 ( 3, 2)
     1      + Wt9 ( 3) * X39 ( 3, 3)
     1      + Wt9 ( 4) * X39 ( 3, 4)
     1      + Wt9 ( 5) * X39 ( 3, 5)
     1      + Wt9 ( 6) * X39 ( 3, 6)
     1      + Wt9 ( 7) * X39 ( 3, 7)
     1      + Wt9 ( 8) * X39 ( 3, 8)
     1      + Wt9 ( 9) * X39 ( 3, 9)

! --------------------------------------------------
! ----------------------------------------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = - 2, ntp_zones + 3" loop.

        END DO  ! End of "DO iyy = - 2, ntp_zones + 3" loop.

      END DO  ! End of "DO izz = - 2, nr_zones + 3" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 9.c) In STEP 9.b, we have already obtained the piecewise linear
! variation of "Ex_zc, Ey_zc, Ez_zc". Here we non-linearly hybridize it with just one
! large central stencil that is piecewise quadratic. This saves a lot of
! computational cost with only minimal degradation in accuracy.

! ----------------------------------------------------------------------

      iorder = 3  ! At third order, all stencils are third order.

      n_cc_10 = 10  ! These are the total number of modes that we will use at 3rd order.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1, icomp,  &
!$OMP&    ixxp1, iyyp1, izzp1, obl0, i, j, nn,                          &
!$OMP&    ixxm2, ixxp2, iyym2, iyyp2, izzm2, izzp2,                     &
!$OMP&    x, y, z, ell_0, ell_j, Coef_0, Coef_j, indx_A10, Amat1010,    &
!$OMP&    Amat1010_Inv, vv_A10, Amat_Row10, Amat1910, Amat1910T,        &
!$OMP&    X10, Y10, R19, ISlo, IShi, Wtlo, Wthi, aa, E_pt,              &
!$OMP&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,               &
!$OMP&    ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,       &
!$OMP&    ixxrr, iyyrr, ixxll, iyyll, ixxuu, iyyuu, ixxdd, iyydd,       &
!$OMP&    Phi10, betalo, gammahi, gammalo, gammaratio, obgammahi, tau)  &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,&
!$OMP&    pow, central_weight, tinynum, smallnum, n_cc_10, n_cc_modes,  &
!$OMP&    ob12, ob180, ob144, iorder,                                   &
!$OMP&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$OMP&    vert_coords_patch, Coef_Basis_Moments,                        &
!$OMP&    istencil_patch, flattengr, Ex_zc, Ey_zc, Ez_zc)               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    pow, central_weight, tinynum, smallnum, n_cc_10, n_cc_modes,  &
!$acc&    ob12, ob180, ob144, iorder,                                   &
!$acc&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$acc&    vert_coords_patch, Coef_Basis_Moments,                        &
!$acc&    istencil_patch, flattengr, Ex_zc, Ey_zc, Ez_zc)

!$acc loop gang collapse(2) independent                                 &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1, icomp,          &
!$acc&    ixxp1, iyyp1, izzp1, obl0, i, j, nn,                          &
!$acc&    ixxm2, ixxp2, iyym2, iyyp2, izzm2, izzp2,                     &
!$acc&    x, y, z, ell_0, ell_j, Coef_0, Coef_j, indx_A10, Amat1010,    &
!$acc&    Amat1010_Inv, vv_A10, Amat_Row10, Amat1910, Amat1910T,        &
!$acc&    X10, Y10, R19, ISlo, IShi, Wtlo, Wthi, aa, E_pt,              &
!$acc&    ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,               &
!$acc&    ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,       &
!$acc&    ixxrr, iyyrr, ixxll, iyyll, ixxuu, iyyuu, ixxdd, iyydd,       &
!$acc&    Phi10, betalo, gammahi, gammalo, gammaratio, obgammahi, tau)  &

! ----------------------------------------------------------------------

!      DO izz = - 2, nr_zones + 3

!        DO iyy = - 2, ntp_zones + 3

!          DO ixx = - 2, ntp_zones + 3

      DO izz = - 1, nr_zones + 2

        DO iyy = - 1, ntp_zones + 2

          DO ixx = - 1, ntp_zones + 2

          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN

! ----------------------------------------------------------------------

      izzm1 = izz - 1
      izzp1 = izz + 1
      izzm2 = izz - 2
      izzp2 = izz + 2

      iyym1 = iyy - 1
      iyyp1 = iyy + 1
      iyym2 = iyy - 2
      iyyp2 = iyy + 2

      ixxm1 = ixx - 1
      ixxp1 = ixx + 1
      ixxm2 = ixx - 2
      ixxp2 = ixx + 2

      ixxr = istencil_patch ( ixx, iyy, 1, 1)   ! For indexing zone one Right.
      iyyr = istencil_patch ( ixx, iyy, 2, 1)

      ixxl = istencil_patch ( ixx, iyy, 1, 2)   ! For indexing zone one Left.
      iyyl = istencil_patch ( ixx, iyy, 2, 2)

      ixxu = istencil_patch ( ixx, iyy, 1, 3)   ! For indexing zone one Up.
      iyyu = istencil_patch ( ixx, iyy, 2, 3)

      ixxd = istencil_patch ( ixx, iyy, 1, 4)   ! For indexing zone one Down.
      iyyd = istencil_patch ( ixx, iyy, 2, 4)

      ixxrr = istencil_patch ( ixx, iyy, 1, 5)  ! For indexing zone two Right.
      iyyrr = istencil_patch ( ixx, iyy, 2, 5)

      ixxll = istencil_patch ( ixx, iyy, 1, 6)  ! For indexing zone two Left.
      iyyll = istencil_patch ( ixx, iyy, 2, 6)

      ixxuu = istencil_patch ( ixx, iyy, 1, 7)  ! For indexing zone two Up.
      iyyuu = istencil_patch ( ixx, iyy, 2, 7)

      ixxdd = istencil_patch ( ixx, iyy, 1, 8)  ! For indexing zone two Down.
      iyydd = istencil_patch ( ixx, iyy, 2, 8)

      ixxru = istencil_patch ( ixx, iyy, 1, 9)  ! For indexing zone that is Right-Up.
      iyyru = istencil_patch ( ixx, iyy, 2, 9)

      ixxlu = istencil_patch ( ixx, iyy, 1, 10)  ! For indexing zone that is Left-Up.
      iyylu = istencil_patch ( ixx, iyy, 2, 10)

      ixxld = istencil_patch ( ixx, iyy, 1, 11)  ! For indexing zone that is Left-Down.
      iyyld = istencil_patch ( ixx, iyy, 2, 11)

      ixxrd = istencil_patch ( ixx, iyy, 1, 12)  ! For indexing zone that is Right-Down.
      iyyrd = istencil_patch ( ixx, iyy, 2, 12)


! ----------------------------------------

          ell_0 = AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)

          obl0 = 1.0 / ell_0  ! Scaling for this zone.

          Coef_0 ( 1: n_cc_10) = Coef_Basis_Moments ( ixx, iyy, izz,
     1                              1: n_cc_10)

! --------------------------------------------------
! We start on the expensive part associated with building "Amat99_Inv".
! --------------------------------------------------

! ------------------------------
! Stencil 9 is centered in all 3 directions and it uses the nineteen zone stencil
! "{ ( ixx, iyy, izz), ( ixxr, iyyr, izz), ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixxd, iyyd, izz),
!    ( ixx, iyy, izzp1), ( ixx, iyy, izzm1),
!    ( ixxru, iyyru, izz), ( ixxlu, iyylu, izz), ( ixxld, iyyld, izz), ( ixxrd, iyyrd, izz),
!    ( ixxr, iyyr, izzp1), ( ixxl, iyyl, izzp1), ( ixxr, iyyr, izzm1), ( ixxl, iyyl, izzm1),
!    ( ixxu, iyyu, izzp1), ( ixxd, iyyd, izzp1), ( ixxu, iyyu, izzm1), ( ixxd, iyyd, izzm1) }".
! ------------------------------

! Build 9 sets of scaled relative coordinates from 9 neighboring zone centroids
! to the central zone's centroid.
! Then scale that with the scaling length for these zones to build 9X9 matrix for Stencil 9.

! ----------
! Stencil 9, Row 1 from "( ixx, iyy, izz)".

          x = 0.0  ;  y = 0.0  ;  z = 0.0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 1, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 2 from "( ixxr, iyyr, izz)".

          x = ( v_centroid_patch ( ixxr, iyyr, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxr, iyyr, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxr, iyyr, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 2, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 3 from "( ixxl, iyyl, izz)".

          x = ( v_centroid_patch ( ixxl, iyyl, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxl, iyyl, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxl, iyyl, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 3, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 4 from "( ixxu, iyyu, izz)".

          x = ( v_centroid_patch ( ixxu, iyyu, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxu, iyyu, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxu, iyyu, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 4, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 5 from "( ixxd, iyyd, izz)".

          x = ( v_centroid_patch ( ixxd, iyyd, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxd, iyyd, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxd, iyyd, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 5, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 6 from "( ixx, iyy, izzp1)".

          x = ( v_centroid_patch ( ixx, iyy, izzp1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixx, iyy, izzp1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixx, iyy, izzp1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 6, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 7 from "( ixx, iyy, izzm1)".

          x = ( v_centroid_patch ( ixx, iyy, izzm1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixx, iyy, izzm1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixx, iyy, izzm1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 7, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 8 from "( ixxru, iyyru, izz)".

          x = ( v_centroid_patch ( ixxru, iyyru, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxru, iyyru, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxru, iyyru, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 8, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 9 from "( ixxlu, iyylu, izz)".

          x = ( v_centroid_patch ( ixxlu, iyylu, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxlu, iyylu, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxlu, iyylu, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 9, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 10 from "( ixxld, iyyld, izz)".

          x = ( v_centroid_patch ( ixxld, iyyld, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxld, iyyld, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxld, iyyld, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 10, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 11 from "( ixxrd, iyyrd, izz)".

          x = ( v_centroid_patch ( ixxrd, iyyrd, izz, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxrd, iyyrd, izz, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxrd, iyyrd, izz, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 11, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 12 from "( ixxr, iyyr, izzp1)".

          x = ( v_centroid_patch ( ixxr, iyyr, izzp1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxr, iyyr, izzp1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxr, iyyr, izzp1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 12, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 13 from "( ixxl, iyyl, izzp1)".

          x = ( v_centroid_patch ( ixxl, iyyl, izzp1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxl, iyyl, izzp1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxl, iyyl, izzp1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 13, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 14 from "( ixxr, iyyr, izzm1)".

          x = ( v_centroid_patch ( ixxr, iyyr, izzm1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxr, iyyr, izzm1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxr, iyyr, izzm1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 14, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 15 from "( ixxl, iyyl, izzm1)".

          x = ( v_centroid_patch ( ixxl, iyyl, izzm1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxl, iyyl, izzm1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxl, iyyl, izzm1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 15, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 16 from "( ixxu, iyyu, izzp1)".

          x = ( v_centroid_patch ( ixxu, iyyu, izzp1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxu, iyyu, izzp1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxu, iyyu, izzp1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 16, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 17 from "( ixxd, iyyd, izzp1)".

          x = ( v_centroid_patch ( ixxd, iyyd, izzp1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxd, iyyd, izzp1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxd, iyyd, izzp1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 17, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 18 from "( ixxu, iyyu, izzm1)".

          x = ( v_centroid_patch ( ixxu, iyyu, izzm1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxu, iyyu, izzm1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxu, iyyu, izzm1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 18, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, Row 19 from "( ixxd, iyyd, izzm1)".

          x = ( v_centroid_patch ( ixxd, iyyd, izzm1, 1)
     1        - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( v_centroid_patch ( ixxd, iyyd, izzm1, 2)
     1        - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( v_centroid_patch ( ixxd, iyyd, izzm1, 3)
     1        - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

          CALL INTERPOLATION_COEFS ( iorder, n_cc_10, x, y, z,
     1                               Coef_0, Amat_Row10 )

          Amat1910 ( 19, :) = Amat_Row10 ( :)

! ----------
! Stencil 9, First build "Amat1910T" from "Amat1910". Use it to build "Amat1010".
! Then invert "Amat1010" and store it in "Amat1010_Inv_9 ( :, :, 9)".
! Since this is a least squares process, we will also need to save the transpose
! "Amat1910T" when we build the coefficients.

          Amat1910T = TRANSPOSE ( Amat1910)

          Amat1010 = MATMUL ( Amat1910T, Amat1910)

          CALL MATINV ( 10, Amat1010, Amat1010_Inv, indx_A10, vv_A10)


! --------------------------------------------------
! By this point, the expensive part associated with building "Amat1010_Inv" is done.
! This array contains the geometry information that is needed for obtaining the
! modes. It is evaluated only once for all the components.
! --------------------------------------------------
! --------------------------------------------------
! Now we focus on obtaining the modes for "Ex_zc" within the zone of interest.
! --------------------------------------------------
! Stencil 9 is centered in all 3 directions and it uses the nineteen zone stencil
! "{ ( ixx, iyy, izz), ( ixxr, iyyr, izz), ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixxd, iyyd, izz),
!    ( ixx, iyy, izzp1), ( ixx, iyy, izzm1),
!    ( ixxru, iyyru, izz), ( ixxlu, iyylu, izz), ( ixxld, iyyld, izz), ( ixxrd, iyyrd, izz),
!    ( ixxr, iyyr, izzp1), ( ixxl, iyyl, izzp1), ( ixxr, iyyr, izzm1), ( ixxl, iyyl, izzm1),
!    ( ixxu, iyyu, izzp1), ( ixxd, iyyd, izzp1), ( ixxu, iyyu, izzm1), ( ixxd, iyyd, izzm1) }".
! ------------------------------

      E_pt = Ex_zc ( ixx, iyy, izz, 1)

      R19 ( 01) = Ex_zc ( ixx, iyy, izz, 1)

      R19 ( 02) = Ex_zc ( ixxr, iyyr, izz, 1)

      R19 ( 03) = Ex_zc ( ixxl, iyyl, izz, 1)

      R19 ( 04) = Ex_zc ( ixxu, iyyu, izz, 1)

      R19 ( 05) = Ex_zc ( ixxd, iyyd, izz, 1)

      R19 ( 06) = Ex_zc ( ixx, iyy, izzp1, 1)

      R19 ( 07) = Ex_zc ( ixx, iyy, izzm1, 1)

      R19 ( 08) = Ex_zc ( ixxru, iyyru, izz, 1)

      R19 ( 09) = Ex_zc ( ixxlu, iyylu, izz, 1)

      R19 ( 10) = Ex_zc ( ixxld, iyyld, izz, 1)

      R19 ( 11) = Ex_zc ( ixxrd, iyyrd, izz, 1)

      R19 ( 12) = Ex_zc ( ixxr, iyyr, izzp1, 1)

      R19 ( 13) = Ex_zc ( ixxl, iyyl, izzp1, 1)

      R19 ( 14) = Ex_zc ( ixxr, iyyr, izzm1, 1)

      R19 ( 15) = Ex_zc ( ixxl, iyyl, izzm1, 1)

      R19 ( 16) = Ex_zc ( ixxu, iyyu, izzp1, 1)

      R19 ( 17) = Ex_zc ( ixxd, iyyd, izzp1, 1)

      R19 ( 18) = Ex_zc ( ixxu, iyyu, izzm1, 1)

      R19 ( 19) = Ex_zc ( ixxd, iyyd, izzm1, 1)

      Y10 = MATMUL ( Amat1910T, R19)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X10 = MATMUL ( Amat1010_Inv, Y10)
#else
!$acc loop seq
          DO i = 1, 10
            x = 0.0
            DO nn = 1, 10
              x = x + Amat1010_Inv ( i, nn) * Y10 ( nn)
            END DO
            X10 ( i) = x
          END DO
#endif

! ------------------------------
! By this point, all nine modes for the large central stencil have been
! evaluated and stored in "X9 ( 1:9)".
! ------------------------------
! Build the smoothness indicator for the large stencil and compare with the smoothness
! indicator for the small piecewise linear stencil, which has already been computed.

! The "Coef_0" variables have been scaled to a unit cube. Different moments
! in the Taylor basis will integrate to different values over such a cube.
! The skewness of the cube sphere zones has not been accounted for here.
! ------------------------------

! Coefficient by which the smoothness indicator for low polynomial is multiplied.
          betalo = flattengr ( ixx, iyy, izz)
     1           + (1.0 - flattengr ( ixx, iyy, izz)) * 2.0

! Build MRWENO linear weights.
          gammahi = flattengr ( ixx, iyy, izz) * 0.6
     1            + (1.0 - flattengr ( ixx, iyy, izz)) * 0.95  ! Linear weight for high polynomial

          gammalo = 1.0 - gammahi  ! Linear weight for low polynomial

          obgammahi = 1.0 / gammahi
          gammaratio = gammalo / gammahi


! Smoothness indicator for low polynomial.
          ISlo = betalo * ( Ex_zc ( ixx, iyy, izz, 2)**2
     1                    + Ex_zc ( ixx, iyy, izz, 3)**2
     1                    + Ex_zc ( ixx, iyy, izz, 4)**2 ) * ob12

! Construct high polynomial for MRWENO by suitably weighted subtracting off
! of the low polynomial.
          Phi10 ( 1:4) = X10 ( 1:4) * obgammahi
     1                - Ex_zc ( ixx, iyy, izz, 1:4) * gammaratio
          Phi10 ( 5:10) = X10 ( 5:10) * obgammahi

! Smoothness indicator for high polynomial, with coefficients "Phi10 ( 1:10)".
          IShi = ( Phi10 ( 2)**2 + Phi10 ( 3)**2 + Phi10 ( 4)**2) * ob12

     1  + ( Phi10 ( 5)**2) * ( ( Coef_0 ( 5) - ob12)**2 + ob180)
     1  + ( Phi10 ( 6)**2) * ( ( Coef_0 ( 6) - ob12)**2 + ob180)
     1  + ( Phi10 ( 7)**2) * ( ( Coef_0 ( 7) - ob12)**2 + ob180)

     1  + ( Phi10 ( 8)**2) * ( Coef_0 ( 8)**2 + ob144)
     1  + ( Phi10 ( 9)**2) * ( Coef_0 ( 9)**2 + ob144)
     1  + ( Phi10 ( 10)**2) * ( Coef_0 ( 10)**2 + ob144)


! Construct non-linear weights using the linear weights and the smoothness indicators.

          tau = ABS ( IShi - ISlo)

          Wtlo = gammalo * (1.0 + tau / (ISlo + smallnum) )
          Wthi = gammahi * (1.0 + tau / (IShi + smallnum) )

          aa = 1.0 / ( Wtlo + Wthi)

          Wtlo = Wtlo * aa
          Wthi = Wthi * aa

! Non-linearly hybridized combination of low and high order polynomials.
! Note that the low order polynomial has only linear slopes.

          Ex_zc ( ixx, iyy, izz, 1:4) = Wthi * Phi10 ( 1:4)
     1                            + Wtlo * Ex_zc ( ixx, iyy, izz, 1:4)

          Ex_zc ( ixx, iyy, izz, 5:10) = Wthi * Phi10 ( 5:10)

! ------------------------------
! Reduce the gradients if flattener is active. For interpolation we can only
! choose between first order and high order.
! ("flattengr" is 0.0 near smooth flow and 1.0 in regions with v. strong
!   shocks and rarefactions.)
! ------------------------------

          x = 1.0 - flattengr ( ixx, iyy, izz)

          Ex_zc ( ixx, iyy, izz, 1) = E_pt * flattengr ( ixx, iyy, izz)
     1                              + Ex_zc ( ixx, iyy, izz, 1) * x

          Ex_zc ( ixx, iyy, izz, 2:10) =
     1                              Ex_zc ( ixx, iyy, izz, 2:10) * x

! --------------------------------------------------
! --------------------------------------------------
! Now we focus on obtaining the modes for "Ey_zc" within the zone of interest.
! --------------------------------------------------
! Stencil 9 is centered in all 3 directions and it uses the nineteen zone stencil
! "{ ( ixx, iyy, izz), ( ixxr, iyyr, izz), ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixxd, iyyd, izz),
!    ( ixx, iyy, izzp1), ( ixx, iyy, izzm1),
!    ( ixxru, iyyru, izz), ( ixxlu, iyylu, izz), ( ixxld, iyyld, izz), ( ixxrd, iyyrd, izz),
!    ( ixxr, iyyr, izzp1), ( ixxl, iyyl, izzp1), ( ixxr, iyyr, izzm1), ( ixxl, iyyl, izzm1),
!    ( ixxu, iyyu, izzp1), ( ixxd, iyyd, izzp1), ( ixxu, iyyu, izzm1), ( ixxd, iyyd, izzm1) }".
! ------------------------------

      E_pt = Ey_zc ( ixx, iyy, izz, 1)

      R19 ( 01) = Ey_zc ( ixx, iyy, izz, 1)

      R19 ( 02) = Ey_zc ( ixxr, iyyr, izz, 1)

      R19 ( 03) = Ey_zc ( ixxl, iyyl, izz, 1)

      R19 ( 04) = Ey_zc ( ixxu, iyyu, izz, 1)

      R19 ( 05) = Ey_zc ( ixxd, iyyd, izz, 1)

      R19 ( 06) = Ey_zc ( ixx, iyy, izzp1, 1)

      R19 ( 07) = Ey_zc ( ixx, iyy, izzm1, 1)

      R19 ( 08) = Ey_zc ( ixxru, iyyru, izz, 1)

      R19 ( 09) = Ey_zc ( ixxlu, iyylu, izz, 1)

      R19 ( 10) = Ey_zc ( ixxld, iyyld, izz, 1)

      R19 ( 11) = Ey_zc ( ixxrd, iyyrd, izz, 1)

      R19 ( 12) = Ey_zc ( ixxr, iyyr, izzp1, 1)

      R19 ( 13) = Ey_zc ( ixxl, iyyl, izzp1, 1)

      R19 ( 14) = Ey_zc ( ixxr, iyyr, izzm1, 1)

      R19 ( 15) = Ey_zc ( ixxl, iyyl, izzm1, 1)

      R19 ( 16) = Ey_zc ( ixxu, iyyu, izzp1, 1)

      R19 ( 17) = Ey_zc ( ixxd, iyyd, izzp1, 1)

      R19 ( 18) = Ey_zc ( ixxu, iyyu, izzm1, 1)

      R19 ( 19) = Ey_zc ( ixxd, iyyd, izzm1, 1)

      Y10 = MATMUL ( Amat1910T, R19)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X10 = MATMUL ( Amat1010_Inv, Y10)
#else
!$acc loop seq
          DO i = 1, 10
            x = 0.0
            DO nn = 1, 10
              x = x + Amat1010_Inv ( i, nn) * Y10 ( nn)
            END DO
            X10 ( i) = x
          END DO
#endif

! ------------------------------
! By this point, all nine modes for the large central stencil have been
! evaluated and stored in "X9 ( 1:9)".
! ------------------------------
! Build the smoothness indicator for the large stencil and compare with the smoothness
! indicator for the small piecewise linear stencil, which has already been computed.

! The "Coef_0" variables have been scaled to a unit cube. Different moments
! in the Taylor basis will integrate to different values over such a cube.
! The skewness of the cube sphere zones has not been accounted for here.
! ------------------------------

! Coefficient by which the smoothness indicator for low polynomial is multiplied.
          betalo = flattengr ( ixx, iyy, izz)
     1           + (1.0 - flattengr ( ixx, iyy, izz)) * 2.0

! Build MRWENO linear weights.
          gammahi = flattengr ( ixx, iyy, izz) * 0.6
     1            + (1.0 - flattengr ( ixx, iyy, izz)) * 0.95  ! Linear weight for high polynomial

          gammalo = 1.0 - gammahi  ! Linear weight for low polynomial

          obgammahi = 1.0 / gammahi
          gammaratio = gammalo / gammahi


! Smoothness indicator for low polynomial.
          ISlo = betalo * ( Ey_zc ( ixx, iyy, izz, 2)**2
     1                    + Ey_zc ( ixx, iyy, izz, 3)**2
     1                    + Ey_zc ( ixx, iyy, izz, 4)**2 ) * ob12

! Construct high polynomial for MRWENO by suitably weighted subtracting off
! of the low polynomial.
          Phi10 ( 1:4) = X10 ( 1:4) * obgammahi
     1                - Ey_zc ( ixx, iyy, izz, 1:4) * gammaratio
          Phi10 ( 5:10) = X10 ( 5:10) * obgammahi

! Smoothness indicator for high polynomial, with coefficients "Phi10 ( 1:10)".
          IShi = ( Phi10 ( 2)**2 + Phi10 ( 3)**2 + Phi10 ( 4)**2) * ob12

     1  + ( Phi10 ( 5)**2) * ( ( Coef_0 ( 5) - ob12)**2 + ob180)
     1  + ( Phi10 ( 6)**2) * ( ( Coef_0 ( 6) - ob12)**2 + ob180)
     1  + ( Phi10 ( 7)**2) * ( ( Coef_0 ( 7) - ob12)**2 + ob180)

     1  + ( Phi10 ( 8)**2) * ( Coef_0 ( 8)**2 + ob144)
     1  + ( Phi10 ( 9)**2) * ( Coef_0 ( 9)**2 + ob144)
     1  + ( Phi10 ( 10)**2) * ( Coef_0 ( 10)**2 + ob144)


! Construct non-linear weights using the linear weights and the smoothness indicators.

          tau = ABS ( IShi - ISlo)

          Wtlo = gammalo * (1.0 + tau / (ISlo + smallnum) )
          Wthi = gammahi * (1.0 + tau / (IShi + smallnum) )

          aa = 1.0 / ( Wtlo + Wthi)

          Wtlo = Wtlo * aa
          Wthi = Wthi * aa

! Non-linearly hybridized combination of low and high order polynomials.
! Note that the low order polynomial has only linear slopes.

          Ey_zc ( ixx, iyy, izz, 1:4) = Wthi * Phi10 ( 1:4)
     1                            + Wtlo * Ey_zc ( ixx, iyy, izz, 1:4)

          Ey_zc ( ixx, iyy, izz, 5:10) = Wthi * Phi10 ( 5:10)

! ------------------------------
! Reduce the gradients if flattener is active. For interpolation we can only
! choose between first order and high order.
! ("flattengr" is 0.0 near smooth flow and 1.0 in regions with v. strong
!   shocks and rarefactions.)
! ------------------------------

          x = 1.0 - flattengr ( ixx, iyy, izz)

          Ey_zc ( ixx, iyy, izz, 1) = E_pt * flattengr ( ixx, iyy, izz)
     1                              + Ey_zc ( ixx, iyy, izz, 1) * x

          Ey_zc ( ixx, iyy, izz, 2:10) =
     1                              Ey_zc ( ixx, iyy, izz, 2:10) * x

! --------------------------------------------------
! --------------------------------------------------
! Now we focus on obtaining the modes for "Ez_zc" within the zone of interest.
! --------------------------------------------------
! Stencil 9 is centered in all 3 directions and it uses the nineteen zone stencil
! "{ ( ixx, iyy, izz), ( ixxr, iyyr, izz), ( ixxl, iyyl, izz), ( ixxu, iyyu, izz), ( ixxd, iyyd, izz),
!    ( ixx, iyy, izzp1), ( ixx, iyy, izzm1),
!    ( ixxru, iyyru, izz), ( ixxlu, iyylu, izz), ( ixxld, iyyld, izz), ( ixxrd, iyyrd, izz),
!    ( ixxr, iyyr, izzp1), ( ixxl, iyyl, izzp1), ( ixxr, iyyr, izzm1), ( ixxl, iyyl, izzm1),
!    ( ixxu, iyyu, izzp1), ( ixxd, iyyd, izzp1), ( ixxu, iyyu, izzm1), ( ixxd, iyyd, izzm1) }".
! ------------------------------

      E_pt = Ez_zc ( ixx, iyy, izz, 1)

      R19 ( 01) = Ez_zc ( ixx, iyy, izz, 1)

      R19 ( 02) = Ez_zc ( ixxr, iyyr, izz, 1)

      R19 ( 03) = Ez_zc ( ixxl, iyyl, izz, 1)

      R19 ( 04) = Ez_zc ( ixxu, iyyu, izz, 1)

      R19 ( 05) = Ez_zc ( ixxd, iyyd, izz, 1)

      R19 ( 06) = Ez_zc ( ixx, iyy, izzp1, 1)

      R19 ( 07) = Ez_zc ( ixx, iyy, izzm1, 1)

      R19 ( 08) = Ez_zc ( ixxru, iyyru, izz, 1)

      R19 ( 09) = Ez_zc ( ixxlu, iyylu, izz, 1)

      R19 ( 10) = Ez_zc ( ixxld, iyyld, izz, 1)

      R19 ( 11) = Ez_zc ( ixxrd, iyyrd, izz, 1)

      R19 ( 12) = Ez_zc ( ixxr, iyyr, izzp1, 1)

      R19 ( 13) = Ez_zc ( ixxl, iyyl, izzp1, 1)

      R19 ( 14) = Ez_zc ( ixxr, iyyr, izzm1, 1)

      R19 ( 15) = Ez_zc ( ixxl, iyyl, izzm1, 1)

      R19 ( 16) = Ez_zc ( ixxu, iyyu, izzp1, 1)

      R19 ( 17) = Ez_zc ( ixxd, iyyd, izzp1, 1)

      R19 ( 18) = Ez_zc ( ixxu, iyyu, izzm1, 1)

      R19 ( 19) = Ez_zc ( ixxd, iyyd, izzm1, 1)

      Y10 = MATMUL ( Amat1910T, R19)

#if( DO_LOOP_FOR_MATMUL == 0 )
          X10 = MATMUL ( Amat1010_Inv, Y10)
#else
!$acc loop seq
          DO i = 1, 10
            x = 0.0
            DO nn = 1, 10
              x = x + Amat1010_Inv ( i, nn) * Y10 ( nn)
            END DO
            X10 ( i) = x
          END DO
#endif

! ------------------------------
! By this point, all nine modes for the large central stencil have been
! evaluated and stored in "X9 ( 1:9)".
! ------------------------------
! Build the smoothness indicator for the large stencil and compare with the smoothness
! indicator for the small piecewise linear stencil, which has already been computed.

! The "Coef_0" variables have been scaled to a unit cube. Different moments
! in the Taylor basis will integrate to different values over such a cube.
! The skewness of the cube sphere zones has not been accounted for here.
! ------------------------------

! Coefficient by which the smoothness indicator for low polynomial is multiplied.
          betalo = flattengr ( ixx, iyy, izz)
     1           + (1.0 - flattengr ( ixx, iyy, izz)) * 2.0

! Build MRWENO linear weights.
          gammahi = flattengr ( ixx, iyy, izz) * 0.6
     1            + (1.0 - flattengr ( ixx, iyy, izz)) * 0.95  ! Linear weight for high polynomial

          gammalo = 1.0 - gammahi  ! Linear weight for low polynomial

          obgammahi = 1.0 / gammahi
          gammaratio = gammalo / gammahi


! Smoothness indicator for low polynomial.
          ISlo = betalo * ( Ez_zc ( ixx, iyy, izz, 2)**2
     1                    + Ez_zc ( ixx, iyy, izz, 3)**2
     1                    + Ez_zc ( ixx, iyy, izz, 4)**2 ) * ob12

! Construct high polynomial for MRWENO by suitably weighted subtracting off
! of the low polynomial.
          Phi10 ( 1:4) = X10 ( 1:4) * obgammahi
     1                - Ez_zc ( ixx, iyy, izz, 1:4) * gammaratio
          Phi10 ( 5:10) = X10 ( 5:10) * obgammahi

! Smoothness indicator for high polynomial, with coefficients "Phi10 ( 1:10)".
          IShi = ( Phi10 ( 2)**2 + Phi10 ( 3)**2 + Phi10 ( 4)**2) * ob12

     1  + ( Phi10 ( 5)**2) * ( ( Coef_0 ( 5) - ob12)**2 + ob180)
     1  + ( Phi10 ( 6)**2) * ( ( Coef_0 ( 6) - ob12)**2 + ob180)
     1  + ( Phi10 ( 7)**2) * ( ( Coef_0 ( 7) - ob12)**2 + ob180)

     1  + ( Phi10 ( 8)**2) * ( Coef_0 ( 8)**2 + ob144)
     1  + ( Phi10 ( 9)**2) * ( Coef_0 ( 9)**2 + ob144)
     1  + ( Phi10 ( 10)**2) * ( Coef_0 ( 10)**2 + ob144)


! Construct non-linear weights using the linear weights and the smoothness indicators.

          tau = ABS ( IShi - ISlo)

          Wtlo = gammalo * (1.0 + tau / (ISlo + smallnum) )
          Wthi = gammahi * (1.0 + tau / (IShi + smallnum) )

          aa = 1.0 / ( Wtlo + Wthi)

          Wtlo = Wtlo * aa
          Wthi = Wthi * aa

! Non-linearly hybridized combination of low and high order polynomials.
! Note that the low order polynomial has only linear slopes.

          Ez_zc ( ixx, iyy, izz, 1:4) = Wthi * Phi10 ( 1:4)
     1                            + Wtlo * Ez_zc ( ixx, iyy, izz, 1:4)

          Ez_zc ( ixx, iyy, izz, 5:10) = Wthi * Phi10 ( 5:10)

! ------------------------------
! Reduce the gradients if flattener is active. For interpolation we can only
! choose between first order and high order.
! ("flattengr" is 0.0 near smooth flow and 1.0 in regions with v. strong
!   shocks and rarefactions.)
! ------------------------------

          x = 1.0 - flattengr ( ixx, iyy, izz)

          Ez_zc ( ixx, iyy, izz, 1) = E_pt * flattengr ( ixx, iyy, izz)
     1                              + Ez_zc ( ixx, iyy, izz, 1) * x

          Ez_zc ( ixx, iyy, izz, 2:10) =
     1                              Ez_zc ( ixx, iyy, izz, 2:10) * x

! --------------------------------------------------
! ----------------------------------------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = - 1, ntp_zones + 2" loop.

        END DO  ! End of "DO iyy = - 1, ntp_zones + 2" loop.

      END DO  ! End of "DO izz = - 1, nr_zones + 2" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 10) We now evaluate the local min and max of the MHD variables on this patch.
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

      pparlmin = 1.0e30
      pparlmax = - 1.0e30

      pperpmin = 1.0e30
      pperpmax = - 1.0e30

      pelectmin = 1.0e30
      pelectmax = - 1.0e30

#if( NFLUID <= 6)
      sp1min = 1.0e30
      sp1max = - 1.0e30

      sp2min = 1.0e30
      sp2max = - 1.0e30

      sp3min = 1.0e30
      sp3max = - 1.0e30

      sp4min = 1.0e30
      sp4max = - 1.0e30

      sp5min = 1.0e30
      sp5max = - 1.0e30

      sp6min = 1.0e30
      sp6max = - 1.0e30
#endif

!$acc kernels
      workgr ( :, :, :) = 0.0
!$acc end kernels

! ------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    a11, a12, a13, a21, x, y, z, bflag, time_use,                 &
!$OMP&    B0_x, B0_y, B0_z)                                             &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,&
!$OMP&    rhofloor, prsfloor, ob_8pi, ob_4pi,                           &
!$OMP&    isospeed, gamma_euler, barotropic_stiffness,                  &
!$OMP&    gamma_soundspeed_gr, msonicgr,                                &
!$OMP&    imask_active_patch, v_centroid_patch,                         &
!$OMP&    xf_area_patch, yf_area_patch, zf_area_patch,                  &
!$OMP&    U, V, gamma_eos_gr, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc,       &
!$OMP&    flattengr, workgr)                                            &

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
!$OMP&           REDUCTION ( MIN: sp1min) REDUCTION ( MAX: sp1max)      &
!$OMP&           REDUCTION ( MIN: sp2min) REDUCTION ( MAX: sp2max)      &
!$OMP&           REDUCTION ( MIN: sp3min) REDUCTION ( MAX: sp3max)      &
!$OMP&           REDUCTION ( MIN: sp4min) REDUCTION ( MAX: sp4max)      &
!$OMP&           REDUCTION ( MIN: sp5min) REDUCTION ( MAX: sp5max)      &
!$OMP&           REDUCTION ( MIN: sp6min) REDUCTION ( MAX: sp6max)      &
!$OMP&       REDUCTION ( MIN: pparlmin) REDUCTION ( MAX: pparlmax)      &
!$OMP&       REDUCTION ( MIN: pperpmin) REDUCTION ( MAX: pperpmax)      &
!$OMP&     REDUCTION ( MIN: pelectmin) REDUCTION ( MAX: pelectmax)      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    rhofloor, prsfloor, ob_8pi, ob_4pi,                           &
!$acc&    isospeed, gamma_euler, barotropic_stiffness,                  &
!$acc&    gamma_soundspeed_gr, msonicgr,                                &
!$acc&    imask_active_patch, v_centroid_patch,                         &
!$acc&    xf_area_patch, yf_area_patch, zf_area_patch,                  &
!$acc&    U, V, gamma_eos_gr, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc,       &
!$acc&    flattengr, workgr,                                            &

!$acc&    V_min_local, V_max_local)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    a11, a12, a13, a21, bflag, time_use, x, y, z,                 &
!$acc&    B0_x, B0_y, B0_z)                                             &

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
!$acc&           REDUCTION ( MIN: bfzmin) REDUCTION ( MAX: bfzmax)      &
!$acc&           REDUCTION ( MIN: sp1min) REDUCTION ( MAX: sp1max)      &
!$acc&           REDUCTION ( MIN: sp2min) REDUCTION ( MAX: sp2max)      &
!$acc&           REDUCTION ( MIN: sp3min) REDUCTION ( MAX: sp3max)      &
!$acc&           REDUCTION ( MIN: sp4min) REDUCTION ( MAX: sp4max)      &
!$acc&           REDUCTION ( MIN: sp5min) REDUCTION ( MAX: sp5max)      &
!$acc&           REDUCTION ( MIN: sp6min) REDUCTION ( MAX: sp6max)      &
!$acc&       REDUCTION ( MIN: pparlmin) REDUCTION ( MAX: pparlmax)      &
!$acc&       REDUCTION ( MIN: pperpmin) REDUCTION ( MAX: pperpmax)      &
!$acc&     REDUCTION ( MIN: pelectmin) REDUCTION ( MAX: pelectmax)

      DO izz = 1, nr_zones

        DO iyy = 1, ntp_zones

          DO ixx = 1, ntp_zones

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! ----------------------------------------
! Take min, max of all the fluid variables.

        rhomin = AMIN1 ( rhomin, V ( ixx, iyy, izz, 1, 1) )
        vlxmin = AMIN1 ( vlxmin, V ( ixx, iyy, izz, 2, 1) )
        vlymin = AMIN1 ( vlymin, V ( ixx, iyy, izz, 3, 1) )
        vlzmin = AMIN1 ( vlzmin, V ( ixx, iyy, izz, 4, 1) )

        bfxmin = AMIN1 ( bfxmin, V ( ixx, iyy, izz, 6, 1) )
        bfymin = AMIN1 ( bfymin, V ( ixx, iyy, izz, 7, 1) )
        bfzmin = AMIN1 ( bfzmin, V ( ixx, iyy, izz, 8, 1) )

        flat_min = AMIN1 ( flat_min, flattengr ( ixx, iyy, izz) )

        rhomax = AMAX1 ( rhomax, V ( ixx, iyy, izz, 1, 1) )
        vlxmax = AMAX1 ( vlxmax, V ( ixx, iyy, izz, 2, 1) )
        vlymax = AMAX1 ( vlymax, V ( ixx, iyy, izz, 3, 1) )
        vlzmax = AMAX1 ( vlzmax, V ( ixx, iyy, izz, 4, 1) )

        bfxmax = AMAX1 ( bfxmax, V ( ixx, iyy, izz, 6, 1) )
        bfymax = AMAX1 ( bfymax, V ( ixx, iyy, izz, 7, 1) )
        bfzmax = AMAX1 ( bfzmax, V ( ixx, iyy, izz, 8, 1) )

        flat_max = AMAX1 ( flat_max, flattengr ( ixx, iyy, izz) )

! ----------

#if( IPDE_SYSTEM == 1)

! ----------

        prsmin = AMIN1 ( prsmin, V ( ixx, iyy, izz, 5, 1) )
        prsmax = AMAX1 ( prsmax, V ( ixx, iyy, izz, 5, 1) )

! ----------

#elif( IPDE_SYSTEM == 2)

! ----------

        pparlmin = AMIN1 ( pparlmin, V ( ixx, iyy, izz, 5, 1) )
        pparlmax = AMAX1 ( pparlmax, V ( ixx, iyy, izz, 5, 1) )

        pperpmin = AMIN1 ( pperpmin, V ( ixx, iyy, izz, 9, 1) )
        pperpmax = AMAX1 ( pperpmax, V ( ixx, iyy, izz, 9, 1) )

! ----------

#elif( IPDE_SYSTEM == 3)

! ----------

        pparlmin = AMIN1 ( pparlmin, V ( ixx, iyy, izz, 5, 1) )
        pparlmax = AMAX1 ( pparlmax, V ( ixx, iyy, izz, 5, 1) )

        pperpmin = AMIN1 ( pperpmin, V ( ixx, iyy, izz, 9, 1) )
        pperpmax = AMAX1 ( pperpmax, V ( ixx, iyy, izz, 9, 1) )

        pelectmin = AMIN1 ( pelectmin, V ( ixx, iyy, izz, 10, 1) )
        pelectmax = AMAX1 ( pelectmax, V ( ixx, iyy, izz, 10, 1) )

! ----------

#endif

! ------------------------------

#if ( IPDE_SYSTEM == 1)

! ----------

#if( NFLUID >= 1)
        sp1min = AMIN1 ( sp1min, V ( ixx, iyy, izz, 8 + 1, 1) )
        sp1max = AMAX1 ( sp1max, V ( ixx, iyy, izz, 8 + 1, 1) )
#endif

#if( NFLUID >= 2)
        sp2min = AMIN1 ( sp2min, V ( ixx, iyy, izz, 8 + 2, 1) )
        sp2max = AMAX1 ( sp2max, V ( ixx, iyy, izz, 8 + 2, 1) )
#endif

#if( NFLUID >= 3)
        sp3min = AMIN1 ( sp3min, V ( ixx, iyy, izz, 8 + 3, 1) )
        sp3max = AMAX1 ( sp3max, V ( ixx, iyy, izz, 8 + 3, 1) )
#endif

#if( NFLUID >= 4)
        sp4min = AMIN1 ( sp4min, V ( ixx, iyy, izz, 8 + 4, 1) )
        sp4max = AMAX1 ( sp4max, V ( ixx, iyy, izz, 8 + 4, 1) )
#endif

#if( NFLUID >= 5)
        sp5min = AMIN1 ( sp5min, V ( ixx, iyy, izz, 8 + 5, 1) )
        sp5max = AMAX1 ( sp5max, V ( ixx, iyy, izz, 8 + 5, 1) )
#endif

#if( NFLUID >= 6)
        sp6min = AMIN1 ( sp6min, V ( ixx, iyy, izz, 8 + 6, 1) )
        sp6max = AMAX1 ( sp6max, V ( ixx, iyy, izz, 8 + 6, 1) )
#endif

! ----------

#endif
! End of "IPDE_SYSTEM == 1"

! ------------------------------

#if ( IPDE_SYSTEM == 2)

! ----------

#if( NVAR >= 10 )
        sp1min = AMIN1 ( sp1min, V ( ixx, iyy, izz, 9 + 1, 1) )
        sp1max = AMAX1 ( sp1max, V ( ixx, iyy, izz, 9 + 1, 1) )
#endif

#if( NVAR >= 11 )
        sp2min = AMIN1 ( sp2min, V ( ixx, iyy, izz, 9 + 2, 1) )
        sp2max = AMAX1 ( sp2max, V ( ixx, iyy, izz, 9 + 2, 1) )
#endif

! ----------

#endif
! End of "IPDE_SYSTEM == 2"

! ------------------------------

#if ( IPDE_SYSTEM == 3)

! ----------

#if( NVAR >= 11 )
        sp1min = AMIN1 ( sp1min, V ( ixx, iyy, izz, 10 + 1, 1) )
        sp1max = AMAX1 ( sp1max, V ( ixx, iyy, izz, 10 + 1, 1) )
#endif

#if( NVAR >= 12 )
        sp2min = AMIN1 ( sp2min, V ( ixx, iyy, izz, 10 + 2, 1) )
        sp2max = AMAX1 ( sp2max, V ( ixx, iyy, izz, 10 + 2, 1) )
#endif

#if( NVAR >= 13 )
        sp3min = AMIN1 ( sp3min, V ( ixx, iyy, izz, 10 + 3, 1) )
        sp3max = AMAX1 ( sp3max, V ( ixx, iyy, izz, 10 + 3, 1) )
#endif

#if( NVAR >= 14 )
        sp4min = AMIN1 ( sp4min, V ( ixx, iyy, izz, 10 + 4, 1) )
        sp4max = AMAX1 ( sp4max, V ( ixx, iyy, izz, 10 + 4, 1) )
#endif

#if( NVAR >= 15 )
        sp5min = AMIN1 ( sp5min, V ( ixx, iyy, izz, 10 + 5, 1) )
        sp5max = AMAX1 ( sp5max, V ( ixx, iyy, izz, 10 + 5, 1) )
#endif

#if( NVAR >= 16 )
        sp6min = AMIN1 ( sp6min, V ( ixx, iyy, izz, 10 + 6, 1) )
        sp6max = AMAX1 ( sp6max, V ( ixx, iyy, izz, 10 + 6, 1) )
#endif

! ----------

#endif
! End of "IPDE_SYSTEM == 3".

! ----------

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

!        a21 = B1_zc ( ixx, iyy, izz, 2)
!     1      + B2_zc ( ixx, iyy, izz, 3)
!     1      + B3_zc ( ixx, iyy, izz, 4)

! ----------

        workgr ( ixx, iyy, izz) = a21

! ------------------------------
! Normalize divergence with the magnitude of magnetic field.

        B0_x = 0.0;  B0_y = 0.0;  B0_z = 0.0

#if( USE_B0B1_FORMULATION == 1)
        x = v_centroid_patch ( ixx, iyy, izz, 1)
        y = v_centroid_patch ( ixx, iyy, izz, 2)
        z = v_centroid_patch ( ixx, iyy, izz, 3)

        bflag = 1         ! 1 == magnetic field;  0 == vector potential.
        time_use = 0.0    ! useful for the case of active rotation.

        CALL GET_MAGNETIC_FIELD ( bflag, time_use, x, y, z,
     1                              B0_x, B0_y, B0_z)
#endif

        a11 = ( V ( ixx, iyy, izz, 6, 1) + B0_x)**2
     1      + ( V ( ixx, iyy, izz, 7, 1) + B0_y)**2
     1      + ( V ( ixx, iyy, izz, 8, 1) + B0_z)**2

! --------------------
#if ( CHOOSE_ADIABATIC == 1)
! Using adiabatic EOS.

          msonicgr ( ixx, iyy, izz) =
     1  gamma_soundspeed_gr ( ixx, iyy, izz) * V ( ixx, iyy, izz, 5, 1)
     1  + ob_4pi * a11

          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz)
     1                                     / V ( ixx, iyy, izz, 1, 1) )

#else
! Using isothermal EOS.

          msonicgr ( ixx, iyy, izz) = isospeed**2
     1     + gamma_euler * barotropic_stiffness
     1                   * V ( ixx, iyy, izz, 1, 1)**(gamma_euler - 1.0)
     1     + ob_4pi * a11 / V ( ixx, iyy, izz, 1, 1)

          msonicgr ( ixx, iyy, izz) = SQRT ( msonicgr ( ixx, iyy, izz))

#endif

! --------------------

        a21 = a21 / MAX ( SQRT ( Bxf ( ixx, iyy, izz, 1)**2
     1                         + Byf ( ixx, iyy, izz, 1)**2
     1                         + Bzf ( ixx, iyy, izz, 1)**2),

     1                    SQRT ( Bxf ( ixxm1, iyy, izz, 1)**2
     1                         + Byf ( ixx, iyym1, izz, 1)**2
     1                         + Bzf ( ixx, iyy, izzm1, 1)**2),

     1             msonicgr ( ixx, iyy, izz)
     1             * SQRT ( 4.0 * 3.14 * V ( ixx, iyy, izz, 1, 1)))


!       WRITE(6,*)"ixx, iyy, izz, div = ",ixx, iyy, izz, a21

        div_min = AMIN1 ( div_min, a21 )
        div_max = AMAX1 ( div_max, a21 )

        END IF  ! End of divergence evaluation.

! ----------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN".

! ----------------------------------------

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

      V_min_local ( 6) = bfxmin
      V_min_local ( 7) = bfymin
      V_min_local ( 8) = bfzmin

      V_min_local ( NVAR+1) = div_min

      V_max_local ( 1) = rhomax
      V_max_local ( 2) = vlxmax
      V_max_local ( 3) = vlymax
      V_max_local ( 4) = vlzmax

      V_max_local ( 6) = bfxmax
      V_max_local ( 7) = bfymax
      V_max_local ( 8) = bfzmax

      V_max_local ( NVAR+1) = div_max

! ----------

#if( IPDE_SYSTEM == 1)

! ----------

      V_min_local ( 5) = prsmin
      V_max_local ( 5) = prsmax

! ----------

#elif( IPDE_SYSTEM == 2)

! ----------

      V_min_local ( 5) = pparlmin
      V_max_local ( 5) = pparlmax

      V_min_local ( 9) = pperpmin
      V_max_local ( 9) = pperpmax

! ----------

#elif( IPDE_SYSTEM == 3)

! ----------

      V_min_local ( 5) = pparlmin
      V_max_local ( 5) = pparlmax

      V_min_local ( 9) = pperpmin
      V_max_local ( 9) = pperpmax

      V_min_local ( 10) = pelectmin
      V_max_local ( 10) = pelectmax

! ----------

#endif

! ------------------------------

#if ( IPDE_SYSTEM == 1)

! ----------

#if( NFLUID >= 1)
      V_min_local ( 8 + 1) = sp1min
      V_max_local ( 8 + 1) = sp1max
#endif

#if( NFLUID >= 2)
      V_min_local ( 8 + 2) = sp2min
      V_max_local ( 8 + 2) = sp2max
#endif

#if( NFLUID >= 3)
      V_min_local ( 8 + 3) = sp3min
      V_max_local ( 8 + 3) = sp3max
#endif

#if( NFLUID >= 4)
      V_min_local ( 8 + 4) = sp4min
      V_max_local ( 8 + 4) = sp4max
#endif

#if( NFLUID >= 5)
      V_min_local ( 8 + 5) = sp5min
      V_max_local ( 8 + 5) = sp5max
#endif

#if( NFLUID >= 6)
      V_min_local ( 8 + 6) = sp6min
      V_max_local ( 8 + 6) = sp6max
#endif

! ----------

#endif
! End of "IPDE_SYSTEM == 1".

! ------------------------------

#if ( IPDE_SYSTEM == 2)

! ----------

#if( NVAR >= 10 )
      V_min_local ( 9 + 1) = sp1min
      V_max_local ( 9 + 1) = sp1max
#endif

#if( NVAR >= 11 )
      V_min_local ( 9 + 2) = sp2min
      V_max_local ( 9 + 2) = sp2max
#endif

! ----------

#endif
! End of "IPDE_SYSTEM == 2".

! ------------------------------

#if ( IPDE_SYSTEM == 3)

! ----------

#if( NVAR >= 11 )
      V_min_local ( 10 + 1) = sp1min
      V_max_local ( 10 + 1) = sp1max
#endif

#if( NVAR >= 12 )
      V_min_local ( 10 + 2) = sp2min
      V_max_local ( 10 + 2) = sp2max
#endif

#if( NVAR >= 13 )
      V_min_local ( 10 + 3) = sp3min
      V_max_local ( 10 + 3) = sp3max
#endif

#if( NVAR >= 14 )
      V_min_local ( 10 + 4) = sp4min
      V_max_local ( 10 + 4) = sp4max
#endif

#if( NVAR >= 15 )
      V_min_local ( 10 + 5) = sp5min
      V_max_local ( 10 + 5) = sp5max
#endif

#if( NVAR >= 16 )
      V_min_local ( 10 + 6) = sp6min
      V_max_local ( 10 + 6) = sp6max
#endif

! ----------

#endif
! End of "IPDE_SYSTEM == 3"
! ------------------------------
! ----------------------------------------------------------------------

! The number of "end data" statements must match the "data copyin" statements.

!$acc end data

! ----------------------------------------------------------------------

      END SUBROUTINE LIMIT_ALL_SLOPES_P12

#endif
! End of "if ( NUMORDER == 3)".






