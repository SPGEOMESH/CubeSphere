! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Focus is on second, third and fourth orders.
#define NUMORDER 4

! Don't need Characteristic Reconstruction.
#define CHARACTERISTIC_RECONSTRUCTION 0

      PROGRAM TEST_HiOrder_DivPreserving_Vect_Field

! This program texts methods for making a divergence-preserving vector field
! without analytical manipulation, i.e. purely by using the CLSQ procedure.

! It would be useful on general meshes, which don't have a Cartesian structure.
! However, we prototype it here on a Cartesian mesh to show that the idea works.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        ixx, iyy, izz, ii,

     1        n_phys_consts,

     1        n_cc_modes, n_fc_nodes, n_cc_components,

     1        NConstraint, NLSQ, NUnknown,

     1        n_fc_modes, n_cc_efld_modes,
     1        num_level, num_grid, num_pe


      PARAMETER ( ix1 = 16, iy1 = 16, iz1 = 16,
     1            ioffx = 4, ioffy = 4, ioffz = 4,

     1            n_phys_consts = 20,

     1            n_fc_modes = 1, n_cc_efld_modes = 1,
     1            num_level = 1, num_grid = 1, num_pe = 1,

#if ( NUMORDER == 2)
     1            n_cc_modes = 4, n_fc_nodes = 5,
     1            NConstraint = 1, NLSQ = 18,
#elif ( NUMORDER == 3)
     1            n_cc_modes = 10, n_fc_nodes = 5,
     1            NConstraint = 4, NLSQ = 60,
#elif ( NUMORDER == 4)
     1            n_cc_modes = 20, n_fc_nodes = 5,
     1            NConstraint = 10, NLSQ = 90,
#endif

     1            NUnknown = 3 * n_cc_modes,  ! The unknowns are always the modes for 3 fields.
     1            n_cc_components = 3 )  ! Because we only have 3 components in a vector field.

      REAL pi, xmin, xmax, ymin, ymax, zmin, zmax, smallnum, A, B, C,
     1     ob3, ob6, tempaa, tempab, tempac, tempad,
     1     tempae, tempaf, tempag, tempah, tempai, tempaj,
     1     nx, ny, nz,

     1     L1_Bfx, Linf_Bfx, L1_Bfy, Linf_Bfy, L1_Bfz, Linf_Bfz,
     1     L1_Div, Linf_Div, L1_u_Div, Linf_u_Div,
     1     L1_Div_lo, Linf_Div_lo

      REAL, DIMENSION ( n_cc_modes) :: u_pt_1, u_pt_2

      REAL phys_consts ( n_phys_consts)

      PARAMETER ( A = 1.0, B = 2.5, C = 1.5)

      REAL coefs ( n_cc_modes), wt_fc ( n_fc_nodes)

! The volumetrically averaged components of the vector field are given to us,
! say from a previous timestep or via initialization.
      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) ::
     1      u

! Flattener variable and work array for evaluating it.
      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) :: flattengr

! There are times when the flow is very troubled when it does not make sense to try
! for very high order. "lo_order_zone_avg = 1" in a zone allows us to drop the order
! in that zone.
      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) :: lo_order_zone_avg

! These are the facial components of the vector field that are also available to us,
! say from a previous timestep or via initialization.
      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz) ::
     1      bfx_fc_gr, bfy_fc_gr, bfz_fc_gr

! We want the volumetric reconstruction of the magnetic field so that it is almost
! divergence-preserving. The quadrature-based integral over the faces should be preserved.
      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_modes) ::
     1      bfx_zc_gr, bfy_zc_gr, bfz_zc_gr

! These nodal values of the magnetic field are obtained using "u" on either side of
! a zone face; with the exception that the quadrature averages to
! "bfx_zc_gr, bfy_zc_gr, bfz_zc_gr".
      REAL, DIMENSION ( 0 - ioffx : ix1 + ioffx,
     1                  0 - ioffy : iy1 + ioffy,
     1                  0 - ioffz : iz1 + ioffz, n_fc_nodes) ::
     1      bfx_fc_nodes, bfy_fc_nodes, bfz_fc_nodes

! ----------------------------------------

! Arrays that specify boundary conditions.

      INTEGER, DIMENSION ( 2,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         x_indx_limits, bcarr_x

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  2,
     1                  1 - ioffz : iz1 + ioffz) ::
     1         y_indx_limits, bcarr_y

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  2) ::
     1         z_indx_limits, bcarr_z

! ----------------------------------------

! Arrays that specify geometry.

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)


! ----------------------------------------------------------------------
! For CLSQ, we have to specify the numbers "NConstraint, NLSQ, NUnknown" for the number
! of constraints, number of LSQ equations and number of unknowns.

      REAL CMat_Row ( NUnknown)

      REAL CMat ( NConstraint, NUnknown), DVec ( NConstraint),
     1     AMat ( NLSQ, NUnknown), BVec ( NLSQ),

     1     ATrans ( NUnknown, NLSQ), ATransA ( NUnknown, NUnknown),
     1     ATransB ( NUnknown), CTrans ( NUnknown, NConstraint),

     1     CLSQMat ( NUnknown + NConstraint, NUnknown + NConstraint),
     1     CLSQVec ( NUnknown + NConstraint),

     1     vv_CLSQ ( NUnknown + NConstraint)

      INTEGER indx_CLSQ ( NUnknown + NConstraint)

! ----------------------------------------------------------------------
! Set up some div-free vector field and its analytical facial and volumetric integrals.

      REAL x, y, z, x1, x2, y1, y2, z1, z2, Bfx_pt, Bfy_pt, Bfz_pt,
     1     Bfx_VolAvg, Bfy_VolAvg, Bfz_VolAvg,
     1     Bfx_XArea, Bfy_YArea, Bfz_ZArea

      Bfx_pt ( x, y, z) = A * Sin ( z) + C * Cos ( y)

      Bfy_pt ( x, y, z) = B * Sin ( x) + A * Cos ( z)

      Bfz_pt ( x, y, z) = C * Sin ( y) + B * Cos ( x)


      Bfx_VolAvg ( x1, x2, y1, y2, z1, z2) =
     1 ( - A * (Cos (z2) - Cos (z1) )*(x2-x1)*(y2-y1)
     1   + C * (Sin (y2) - Sin (y1) )*(x2-x1)*(z2-z1) )
     1       / ((x2-x1)*(y2-y1)*(z2-z1))

      Bfy_VolAvg ( x1, x2, y1, y2, z1, z2) =
     1 ( - B * (Cos (x2) - Cos (x1))*(y2-y1)*(z2-z1)
     1   + A * (Sin (z2) - Sin (z1))*(x2-x1)*(y2-y1) )
     1       / ((x2-x1)*(y2-y1)*(z2-z1))

      Bfz_VolAvg ( x1, x2, y1, y2, z1, z2) =
     1 ( - C * (Cos (y2) - Cos (y1))*(x2-x1)*(z2-z1)
     1   + B * (Sin (x2) - Sin (x1))*(y2-y1)*(z2-z1) )
     1       / ((x2-x1)*(y2-y1)*(z2-z1))


      Bfx_XArea ( x, y1, y2, z1, z2) = 
     1 ( - A * ( Cos ( z2) - Cos ( z1))*(y2-y1)
     1   + C * ( Sin ( y2) - Sin ( y1))*(z2-z1) )
     1       / ((y2-y1)*(z2-z1))

      Bfy_YArea ( x1, x2, y, z1, z2) = 
     1 ( - B * (Cos (x2) - Cos (x1))*(z2-z1)
     1   + A * (Sin (z2) - Sin (z1))*(x2-x1) )
     1       / ((x2-x1)*(z2-z1))

      Bfz_ZArea ( x1, x2, y1, y2, z) =
     1 ( - C * (Cos (y2) - Cos (y1))*(x2-x1)
     1   + B * (Sin (x2) - Sin (x1))*(y2-y1) )
     1       / ((x2-x1)*(y2-y1))

! ----------------------------------------------------------------------

      smallnum = 1.0e-12  ;  phys_consts ( 3) = smallnum

      ob3 = 1.0 / 3.0  ;  ob6 = 1.0 / 6.0

      lo_order_zone_avg = 0  ! Here we don't want any lowering of the zone average values.
      flattengr = 0.0  ! We don't want any flattening in the reconstruction.

      pi = 4.0 * ATAN ( 1.0)

      xmin = 0.0
      xmax = 2.0 * pi

      ymin = 0.0
      ymax = 2.0 * pi

      zmin = 0.0
      zmax = 2.0 * pi


! Set up the zone-averaged integrals of the modes within a zone.
! This is done so that the zone-averaged higher moments average to zero.
      coefs = 0.0
#if ( NUMORDER >= 3)
      coefs ( 5:7) = 1.0 / 12.0
#endif

#if ( NUMORDER == 1)
      wt_fc ( 1) = 1.0
#else
      wt_fc ( 1) = 1.0 / 3.0  ! Central node of 5 point Simpson quadrature.
      wt_fc ( 2:) = 1.0 / 6.0  ! This is Simpson quadrature over a 2D rectangle using 5 nodes.
#endif

! ----------------------------------------------------------------------
! Set up mesh and boundary conditions.

! ----------

      dx_save ( 1) = ( xmax - xmin) / ix1
      xb_save ( 0) = xmin

      DO ixx = 1, ix1 + ioffx
      dx_save ( ixx) = dx_save ( 1)
      xb_save ( ixx) = xb_save ( ixx - 1) + dx_save ( ixx)
      xc_save ( ixx) = 0.5 * ( xb_save ( ixx) + xb_save ( ixx - 1))
      END DO

      DO ixx = 0, 1 - ioffx, -1
      dx_save ( ixx) = dx_save ( 1)
      xb_save ( ixx - 1) = xb_save ( ixx) - dx_save ( ixx)
      xc_save ( ixx) = 0.5 * ( xb_save ( ixx) + xb_save ( ixx - 1))
      END DO

! ----------

      dy_save ( 1) = ( ymax - ymin) / iy1
      yb_save ( 0) = ymin

      DO iyy = 1, iy1 + ioffy
      dy_save ( iyy) = dy_save ( 1)
      yb_save ( iyy) = yb_save ( iyy - 1) + dy_save ( iyy)
      yc_save ( iyy) = 0.5 * ( yb_save ( iyy) + yb_save ( iyy - 1))
      END DO

      DO iyy = 0, 1 - ioffy, -1
      dy_save ( iyy) = dy_save ( 1)
      yb_save ( iyy - 1) = yb_save ( iyy) - dy_save ( iyy)
      yc_save ( iyy) = 0.5 * ( yb_save ( iyy) + yb_save ( iyy - 1))
      END DO

! ----------

      dz_save ( 1) = ( zmax - zmin) / iz1
      zb_save ( 0) = zmin

      DO izz = 1, iz1 + ioffz
      dz_save ( izz) = dz_save ( 1)
      zb_save ( izz) = zb_save ( izz - 1) + dz_save ( izz)
      zc_save ( izz) = 0.5 * ( zb_save ( izz) + zb_save ( izz - 1))
      END DO

      DO izz = 0, 1 - ioffz, -1
      dz_save ( izz) = dz_save ( 1)
      zb_save ( izz - 1) = zb_save ( izz) - dz_save ( izz)
      zc_save ( izz) = 0.5 * ( zb_save ( izz) + zb_save ( izz - 1))
      END DO

! ------------------------------

! Bottom and top yz planes.

      DO iyy = 1 - ioffy, iy1 + ioffy

        DO izz = 1 - ioffz, iz1 + ioffz

        x_indx_limits ( 1, iyy, izz) = 1
        x_indx_limits ( 2, iyy, izz) = ix1

          bcarr_x ( 1, iyy, izz) = 4

          bcarr_x ( 2, iyy, izz) = 4

        END DO ! End DO izz = loop.

      END DO ! End DO iyy = loop.

! ------------------------------

! Bottom and top xz planes.

      DO ixx = 1 - ioffx, ix1 + ioffx

        DO izz = 1 - ioffz, iz1 + ioffz

        y_indx_limits ( ixx, 1, izz) = 1
        y_indx_limits ( ixx, 2, izz) = iy1

          bcarr_y ( ixx, 1, izz) = 4

          bcarr_y ( ixx, 2, izz) = 4

        END DO ! End DO izz = loop.

      END DO ! End DO ixx = loop.

! ------------------------------

! Bottom and top xy planes.

      DO ixx = 1 - ioffx, ix1 + ioffx

        DO iyy = 1 - ioffy, iy1 + ioffy

        z_indx_limits ( ixx, iyy, 1) = 1
        z_indx_limits ( ixx, iyy, 2) = iz1

          bcarr_z ( ixx, iyy, 1) = 4

          bcarr_z ( ixx, iyy, 2) = 4

        END DO ! End DO iyy = loop.

      END DO ! End DO ixx = loop.

! ----------------------------------------------------------------------
! Initialize the volume-averaged magnetic fields in the three components of "u"
! as well as the facially-averaged magnetic fields "bfx_fc_gr, bfy_fc_gr, bfz_fc_gr".

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

          x1 = xb_save ( ixx - 1)  ;  x2 = xb_save ( ixx)
          y1 = yb_save ( iyy - 1)  ;  y2 = yb_save ( iyy)
          z1 = zb_save ( izz - 1)  ;  z2 = zb_save ( izz)

          u ( ixx, iyy, izz, 1, 1)= Bfx_VolAvg ( x1, x2, y1, y2, z1, z2)

          u ( ixx, iyy, izz, 2, 1)= Bfy_VolAvg ( x1, x2, y1, y2, z1, z2)

          u ( ixx, iyy, izz, 3, 1)= Bfz_VolAvg ( x1, x2, y1, y2, z1, z2)

          END DO

        END DO

      END DO

! ----------

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 0 - ioffx, ix1 + ioffx

          x = xb_save ( ixx)
          y1 = yb_save ( iyy - 1)  ;  y2 = yb_save ( iyy)
          z1 = zb_save ( izz - 1)  ;  z2 = zb_save ( izz)

          bfx_fc_gr ( ixx, iyy, izz) = Bfx_XArea ( x, y1, y2, z1, z2)

          END DO

        END DO

      END DO

! ----------

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 0 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

          y = yb_save ( iyy)
          x1 = xb_save ( ixx - 1)  ;  x2 = xb_save ( ixx)
          z1 = zb_save ( izz - 1)  ;  z2 = zb_save ( izz)

          bfy_fc_gr ( ixx, iyy, izz) = Bfy_YArea ( x1, x2, y, z1, z2)

          END DO

        END DO

      END DO

! ----------

      DO izz = 0 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

          z = zb_save ( izz)
          x1 = xb_save ( ixx - 1)  ;  x2 = xb_save ( ixx)
          y1 = yb_save ( iyy - 1)  ;  y2 = yb_save ( iyy)

          bfz_fc_gr ( ixx, iyy, izz) = Bfz_ZArea ( x1, x2, y1, y2, z)

          END DO

        END DO

      END DO

! ----------------------------------------------------------------------
! STEP I) Reconstruct all the three zone-centered variables using a standard
! finite-volume reconstruction.

#if (NUMORDER == 2)
      CALL LIMIT_RECONST_ALL_ZC_P1
#elif (NUMORDER == 3)
      CALL LIMIT_RECONST_ALL_ZC_P2
#elif (NUMORDER == 4)
      CALL LIMIT_RECONST_ALL_ZC_P3
#endif

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_modes, n_fc_modes, n_cc_efld_modes, n_cc_components,
     1        num_level, num_grid, num_pe,
     1        n_phys_consts, phys_consts,

     1        u, flattengr, lo_order_zone_avg,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------
! STEP II) Now that "u (:,:,;,1,:)" is loaded with all the volumetric modes for "Bfx",
! We can evaluate the nodal values for "bfx_fc_nodes (:,:,:,:)" at the x-faces.
! We then want to make sure that the integral over those nodes matches the facially
! averaged value of "bfx_fc_gr (:,:,:)".

      ii = 1  ! Accessing "u (:,:,;,1,:)".

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 0 - ioffx, ix1 + ioffx

! ----------

          u_pt_1 ( :) = u ( ixx, iyy, izz, ii, :)

          u_pt_2 ( :) = u ( ixx + 1, iyy, izz, ii, :)

! ----------
! Central node evaluated from both sides of an x-face.

          CALL v_polynomial_modes ( 0.5, 0.0, 0.0,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( - 0.5, 0.0, 0.0,
     1                              n_cc_modes, u_pt_2, tempab)

          bfx_fc_nodes ( ixx, iyy, izz, 1) = 0.5 * ( tempaa + tempab)

! ----------
! y-up node evaluated from both sides of an x-face.

          CALL v_polynomial_modes ( 0.5, 0.5, 0.0,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( - 0.5, 0.5, 0.0,
     1                              n_cc_modes, u_pt_2, tempab)

          bfx_fc_nodes ( ixx, iyy, izz, 2) = 0.5 * ( tempaa + tempab)

! ----------
! y-down node evaluated from both sides of an x-face.

          CALL v_polynomial_modes ( 0.5, - 0.5, 0.0,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( - 0.5, - 0.5, 0.0,
     1                              n_cc_modes, u_pt_2, tempab)

          bfx_fc_nodes ( ixx, iyy, izz, 3) = 0.5 * ( tempaa + tempab)

! ----------
! z-top node evaluated from both sides of an x-face.

          CALL v_polynomial_modes ( 0.5, 0.0, 0.5,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( - 0.5, 0.0, 0.5,
     1                              n_cc_modes, u_pt_2, tempab)

          bfx_fc_nodes ( ixx, iyy, izz, 4) = 0.5 * ( tempaa + tempab)

! ----------
! z-bottom node evaluated from both sides of an x-face.

          CALL v_polynomial_modes ( 0.5, 0.0, - 0.5,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( - 0.5, 0.0, - 0.5,
     1                              n_cc_modes, u_pt_2, tempab)

          bfx_fc_nodes ( ixx, iyy, izz, 5) = 0.5 * ( tempaa + tempab)

! ----------
! Now normalize so that the integral over these nodes adds up to "bfx_fc_gr (:,:,:)".

          tempaa = ob3 * bfx_fc_nodes ( ixx, iyy, izz, 1)
     1           + ob6 * ( bfx_fc_nodes ( ixx, iyy, izz, 2)
     1                   + bfx_fc_nodes ( ixx, iyy, izz, 3)
     1                   + bfx_fc_nodes ( ixx, iyy, izz, 4)
     1                   + bfx_fc_nodes ( ixx, iyy, izz, 5) )

          tempab = bfx_fc_gr ( ixx, iyy, izz) - tempaa

          bfx_fc_nodes ( ixx, iyy, izz, 1) =
     1    bfx_fc_nodes ( ixx, iyy, izz, 1) + tempab

          bfx_fc_nodes ( ixx, iyy, izz, 2) =
     1    bfx_fc_nodes ( ixx, iyy, izz, 2) + tempab

          bfx_fc_nodes ( ixx, iyy, izz, 3) =
     1    bfx_fc_nodes ( ixx, iyy, izz, 3) + tempab

          bfx_fc_nodes ( ixx, iyy, izz, 4) =
     1    bfx_fc_nodes ( ixx, iyy, izz, 4) + tempab

          bfx_fc_nodes ( ixx, iyy, izz, 5) =
     1    bfx_fc_nodes ( ixx, iyy, izz, 5) + tempab

! ----------

          END DO

        END DO

      END DO

! ----------------------------------------------------------------------
! Now that "u (:,:,;,2,:)" is loaded with all the volumetric modes for "Bfy",
! We can evaluate the nodal values for "bfy_fc_nodes (:,:,:,:)" at the y-faces. We then
! want to make sure that the integral over those nodes matches the facially
! averaged value of "bfy_fc_gr (:,:,:)".

      ii = 2  ! Accessing "u (:,:,;,2,:)".

      DO izz = 1 - ioffz, iz1 + ioffz

        DO ixx = 1 - ioffx, ix1 + ioffx

          DO iyy = 0 - ioffy, iy1 + ioffy

! ----------

          u_pt_1 ( :) = u ( ixx, iyy, izz, ii, :)

          u_pt_2 ( :) = u ( ixx, iyy + 1, izz, ii, :)

! ----------
! Central node evaluated from both sides of a y-face.

          CALL v_polynomial_modes ( 0.0, 0.5, 0.0,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( 0.0, - 0.5, 0.0,
     1                              n_cc_modes, u_pt_2, tempab)

          bfy_fc_nodes ( ixx, iyy, izz, 1) = 0.5 * ( tempaa + tempab)

! ----------
! x-right node evaluated from both sides of a y-face.

          CALL v_polynomial_modes ( 0.5, 0.5, 0.0,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( 0.5, - 0.5, 0.0,
     1                              n_cc_modes, u_pt_2, tempab)

          bfy_fc_nodes ( ixx, iyy, izz, 2) = 0.5 * ( tempaa + tempab)

! ----------
! x-left node evaluated from both sides of a y-face.

          CALL v_polynomial_modes ( - 0.5, 0.5, 0.0,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( - 0.5, - 0.5, 0.0,
     1                              n_cc_modes, u_pt_2, tempab)

          bfy_fc_nodes ( ixx, iyy, izz, 3) = 0.5 * ( tempaa + tempab)

! ----------
! z-top node evaluated from both sides of a y-face.

          CALL v_polynomial_modes ( 0.0, 0.5, 0.5,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( 0.0, - 0.5, 0.5,
     1                              n_cc_modes, u_pt_2, tempab)

          bfy_fc_nodes ( ixx, iyy, izz, 4) = 0.5 * ( tempaa + tempab)

! ----------
! z-bottom node evaluated from both sides of a y-face.

          CALL v_polynomial_modes ( 0.0, 0.5, - 0.5,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( 0.0, - 0.5, - 0.5,
     1                              n_cc_modes, u_pt_2, tempab)

          bfy_fc_nodes ( ixx, iyy, izz, 5) = 0.5 * ( tempaa + tempab)

! ----------
! Now normalize so that the integral over these nodes adds up to "bfy_fc_gr (:,:,:)".

          tempaa = ob3 * bfy_fc_nodes ( ixx, iyy, izz, 1)
     1           + ob6 * ( bfy_fc_nodes ( ixx, iyy, izz, 2)
     1                   + bfy_fc_nodes ( ixx, iyy, izz, 3)
     1                   + bfy_fc_nodes ( ixx, iyy, izz, 4)
     1                   + bfy_fc_nodes ( ixx, iyy, izz, 5) )

          tempab = bfy_fc_gr ( ixx, iyy, izz) - tempaa

          bfy_fc_nodes ( ixx, iyy, izz, 1) =
     1    bfy_fc_nodes ( ixx, iyy, izz, 1) + tempab

          bfy_fc_nodes ( ixx, iyy, izz, 2) =
     1    bfy_fc_nodes ( ixx, iyy, izz, 2) + tempab

          bfy_fc_nodes ( ixx, iyy, izz, 3) =
     1    bfy_fc_nodes ( ixx, iyy, izz, 3) + tempab

          bfy_fc_nodes ( ixx, iyy, izz, 4) =
     1    bfy_fc_nodes ( ixx, iyy, izz, 4) + tempab

          bfy_fc_nodes ( ixx, iyy, izz, 5) =
     1    bfy_fc_nodes ( ixx, iyy, izz, 5) + tempab

! ----------

          END DO

        END DO

      END DO

! ----------------------------------------------------------------------
! Now that "u (:,:,;,3,:)" is loaded with all the volumetric modes for "Bfz",
! We can evaluate the nodal values for "bfz_fc_nodes (:,:,:,:)" at the z-faces. We then
! want to make sure that the integral over those nodes matches the facially
! averaged value of "bfz_fc_gr (:,:,:)".

      ii = 3  ! Accessing "u (:,:,;,3,:)".

      DO ixx = 1 - ioffx, ix1 + ioffx

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO izz = 0 - ioffz, iz1 + ioffz

! ----------

          u_pt_1 ( :) = u ( ixx, iyy, izz, ii, :)

          u_pt_2 ( :) = u ( ixx, iyy, izz + 1, ii, :)

! ----------
! Central node evaluated from both sides of a z-face.

          CALL v_polynomial_modes ( 0.0, 0.0, 0.5,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( 0.0, 0.0, - 0.5,
     1                              n_cc_modes, u_pt_2, tempab)

          bfz_fc_nodes ( ixx, iyy, izz, 1) = 0.5 * ( tempaa + tempab)

! ----------
! x-right node evaluated from both sides of a z-face.

          CALL v_polynomial_modes ( 0.5, 0.0, 0.5,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( 0.5, 0.0, - 0.5,
     1                              n_cc_modes, u_pt_2, tempab)

          bfz_fc_nodes ( ixx, iyy, izz, 2) = 0.5 * ( tempaa + tempab)

! ----------
! x-left node evaluated from both sides of a z-face.

          CALL v_polynomial_modes ( - 0.5, 0.0, 0.5,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( - 0.5, 0.0, - 0.5,
     1                              n_cc_modes, u_pt_2, tempab)

          bfz_fc_nodes ( ixx, iyy, izz, 3) = 0.5 * ( tempaa + tempab)

! ----------
! y-up node evaluated from both sides of a z-face.

          CALL v_polynomial_modes ( 0.0, 0.5, 0.5,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( 0.0, 0.5, - 0.5,
     1                              n_cc_modes, u_pt_2, tempab)

          bfz_fc_nodes ( ixx, iyy, izz, 4) = 0.5 * ( tempaa + tempab)

! ----------
! y-down node evaluated from both sides of a z-face.

          CALL v_polynomial_modes ( 0.0, - 0.5, 0.5,
     1                              n_cc_modes, u_pt_1, tempaa)

          CALL v_polynomial_modes ( 0.0, - 0.5, - 0.5,
     1                              n_cc_modes, u_pt_2, tempab)

          bfz_fc_nodes ( ixx, iyy, izz, 5) = 0.5 * ( tempaa + tempab)

! ----------
! Now normalize so that the integral over these nodes adds up to "bfz_fc_gr (:,:,:)".

          tempaa = ob3 * bfz_fc_nodes ( ixx, iyy, izz, 1)
     1           + ob6 * ( bfz_fc_nodes ( ixx, iyy, izz, 2)
     1                   + bfz_fc_nodes ( ixx, iyy, izz, 3)
     1                   + bfz_fc_nodes ( ixx, iyy, izz, 4)
     1                   + bfz_fc_nodes ( ixx, iyy, izz, 5) )

          tempab = bfz_fc_gr ( ixx, iyy, izz) - tempaa

          bfz_fc_nodes ( ixx, iyy, izz, 1) =
     1    bfz_fc_nodes ( ixx, iyy, izz, 1) + tempab

          bfz_fc_nodes ( ixx, iyy, izz, 2) =
     1    bfz_fc_nodes ( ixx, iyy, izz, 2) + tempab

          bfz_fc_nodes ( ixx, iyy, izz, 3) =
     1    bfz_fc_nodes ( ixx, iyy, izz, 3) + tempab

          bfz_fc_nodes ( ixx, iyy, izz, 4) =
     1    bfz_fc_nodes ( ixx, iyy, izz, 4) + tempab

          bfz_fc_nodes ( ixx, iyy, izz, 5) =
     1    bfz_fc_nodes ( ixx, iyy, izz, 5) + tempab

! ----------

          END DO

        END DO

      END DO

! ----------------------------------------------------------------------
! STEP III) Use the CLSQ procedure to obtain the almost div-preserving modes.

! ----------------------------------------------------------------------
#if (NUMORDER == 2)
! STEP III.a) At second order, we use the 6 face-centered values at the 6 faces
! as constraints.
! Because there is only linear variation in the faces at second order, we choose
! two more (non-orthogonal) nodal values at each face as least squares input.
! "NConstraint = 1, NLSQ = 18, NUnknown = 12".

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------------------------------------

          CMat = 0.0  ;  AMat = 0.0  ;  DVec = 0.0  ;  BVec = 0.0

! ----------------------------------------
! Load up the matrix of constraints "CMat ( NConstraint, NUnknown)" and its
! right hand side "DVec ( NConstraint)".

      CMat ( 1, 2) = 1.0  ;  CMat ( 1, 7) = 1.0  ;  CMat ( 1, 12) = 1.0
      DVec ( 1) = 0.0

! ----------------------------------------
! Load up the matrix of LSQ equations "AMat ( NLSQ, NUnknown)" and their
! right hand side "BVec ( NLSQ)".

! --------------------

      x = 0.5 ; y = 0.0 ; z = 0.0  ! Central node at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 1, :) = CMat_Row ( :)
      BVec ( 1) = bfx_fc_nodes ( ixx, iyy, izz, 1)

! ----------

      x = - 0.5 ; y = 0.0 ; z = 0.0  ! Central node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 2, :) = CMat_Row ( :)
      BVec ( 2) = bfx_fc_nodes ( ixx - 1, iyy, izz, 1)

! --------------------

      x = 0.0 ; y = 0.5 ; z = 0.0  ! Central node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 3, :) = CMat_Row ( :)
      BVec ( 3) = bfy_fc_nodes ( ixx, iyy, izz, 1)

! ----------

      x = 0.0 ; y = - 0.5 ; z = 0.0  ! Central node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 4, :) = CMat_Row ( :)
      BVec ( 4) = bfy_fc_nodes ( ixx, iyy - 1, izz, 1)

! --------------------

      x = 0.0 ; y = 0.0 ; z = 0.5  ! Central node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 5, :) = CMat_Row ( :)
      BVec ( 5) = bfz_fc_nodes ( ixx, iyy, izz, 1)

! ----------

      x = 0.0 ; y = 0.0 ; z = - 0.5  ! Central node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 6, :) = CMat_Row ( :)
      BVec ( 6) = bfz_fc_nodes ( ixx, iyy, izz - 1, 1)

! --------------------
! Now in LSQ sense match the zone-centered modes of "u" that correspond to bfx.

      AMat ( 7, 1) = 1.0
      BVec ( 7) = u ( ixx, iyy, izz, 1, 1)

      AMat ( 8, 2) = 1.0
      BVec ( 8) = u ( ixx, iyy, izz, 1, 2)

      AMat ( 9, 3) = 1.0
      BVec ( 9) = u ( ixx, iyy, izz, 1, 3)

      AMat ( 10, 4) = 1.0
      BVec ( 10) = u ( ixx, iyy, izz, 1, 4)

! --------------------
! Now in LSQ sense match the zone-centered modes of "u" that correspond to bfy.

      AMat ( 11, 5) = 1.0
      BVec ( 11) = u ( ixx, iyy, izz, 2, 1)

      AMat ( 12, 6) = 1.0
      BVec ( 12) = u ( ixx, iyy, izz, 2, 2)

      AMat ( 13, 7) = 1.0
      BVec ( 13) = u ( ixx, iyy, izz, 2, 3)

      AMat ( 14, 8) = 1.0
      BVec ( 14) = u ( ixx, iyy, izz, 2, 4)

! --------------------
! Now in LSQ sense match the zone-centered modes of "u" that correspond to bfz.

      AMat ( 15, 9) = 1.0
      BVec ( 15) = u ( ixx, iyy, izz, 3, 1)

      AMat ( 16, 10) = 1.0
      BVec ( 16) = u ( ixx, iyy, izz, 3, 2)

      AMat ( 17, 11) = 1.0
      BVec ( 17) = u ( ixx, iyy, izz, 3, 3)

      AMat ( 18, 12) = 1.0
      BVec ( 18) = u ( ixx, iyy, izz, 3, 4)

! ----------------------------------------
! Build transpose of "AMat", "ATrans ( NUnknown, NLSQ)". From that, build
! "ATransA ( NUnknown, NUnknown)" and "ATransB ( NUnknown)". Also build
! "CTrans ( NUnknown, NConstraint)".

      ATrans = TRANSPOSE ( AMat)

      ATransA = MATMUL ( ATrans, AMat)

      ATransB = MATMUL ( ATrans, BVec)

      CTrans = TRANSPOSE ( CMat)

! ----------------------------------------
! Assemble the matrix "CLSQMat ( NUnknown + NConstraint, NUnknown + NConstraint)"
! and its right hand side "CLSQVec ( NUnknown + NConstraint)".

      CLSQMat ( 1: NUnknown, 1: NUnknown) =
     1                            ATransA ( 1: NUnknown, 1: NUnknown)

      CLSQMat ( 1: NUnknown, NUnknown + 1: NUnknown + NConstraint) =
     1                            CTrans ( 1: NUnknown, 1: NConstraint)

      CLSQMat ( NUnknown + 1: NUnknown + NConstraint, 1: NUnknown) =
     1                            CMat ( 1: NConstraint, 1: NUnknown)

      CLSQMat ( NUnknown + 1: NUnknown + NConstraint, 
     1          NUnknown + 1: NUnknown + NConstraint) = 0.0

! ----------

      CLSQVec ( 1: NUnknown) = ATransB ( 1: NUnknown)

      CLSQVec ( NUnknown + 1: NUnknown + NConstraint) = 
     1                             DVec ( 1: NConstraint)

! ----------------------------------------
! Solve the system "CLSQMat . (x, z)^T = CLSQVec". Extract from that the vector of
! unknowns.

      CALL MATSOLVE ( NUnknown + NConstraint, CLSQMat, CLSQVec,
     1                indx_CLSQ, vv_CLSQ)

! ----------------------------------------
! Now hand back the solution to the zone-centered reconstruction of the magnetic
! field that is closer to being divergence-free.

      bfx_zc_gr ( ixx, iyy, izz, 1: 4) = CLSQVec ( 1: 4)

      bfy_zc_gr ( ixx, iyy, izz, 1: 4) = CLSQVec ( 5: 8)

      bfz_zc_gr ( ixx, iyy, izz, 1: 4) = CLSQVec ( 9: 12)

! ----------------------------------------

          END DO

        END DO

      END DO

#endif
! End of "(NUMORDER == 2)".
! ----------------------------------------------------------------------
#if (NUMORDER == 3)
! STEP III.b) At third order, we use the 6 face-integrated values at the 6 faces
! as constraints. We also impose the linear variations of the divergence.
! We then want all the modes to be as close to the 30 volumetrically reconstructed
! modes in a least squares sense.
! "NConstraint = 4, NLSQ = 60, NUnknown = 30".

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------------------------------------

          CMat = 0.0  ;  AMat = 0.0  ;  DVec = 0.0  ;  BVec = 0.0

! ----------------------------------------
! Load up the matrix of constraints "CMat ( NConstraint, NUnknown)" and its
! right hand side "DVec ( NConstraint)".

! --------------------
! Constant divergence-free constraint.

      CMat ( 1, 2) = 1.0 ; CMat ( 1, 13) = 1.0 ; CMat ( 1, 24) = 1.0
      DVec ( 1) = 0.0

! --------------------
! Linear in x value of divergence-free constraint.

      CMat ( 2, 5) = 2.0 ; CMat ( 2, 18) = 1.0 ; CMat ( 2, 30) = 1.0
      DVec ( 2) = 0.0

! ----------
! Linear in y value of divergence-free constraint.

      CMat ( 3, 8) = 1.0 ; CMat ( 3, 16) = 2.0 ; CMat ( 3, 29) = 1.0
      DVec ( 3) = 0.0

! ----------
! Linear in z value of divergence-free constraint.

      CMat ( 4, 10) = 1.0 ; CMat ( 4, 19) = 1.0 ; CMat ( 4, 27) = 2.0
      DVec ( 4) = 0.0

! ----------------------------------------
! Load up the matrix of LSQ equations "AMat ( NLSQ, NUnknown)" and their
! right hand side "BVec ( NLSQ)".

! --------------------
! We require matching 5 nodal values at the right x-face in LSQ sense.

! --------------------

      x = 0.5 ; y = 0.0 ; z = 0.0  ! Central node at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 1, :) = CMat_Row ( :)
      BVec ( 1) = bfx_fc_nodes ( ixx, iyy, izz, 1)

! ----------

      x = 0.5 ; y = 0.5 ; z = 0.0  ! y-up node at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 2, :) = CMat_Row ( :)
      BVec ( 2) = bfx_fc_nodes ( ixx, iyy, izz, 2)

! ----------

      x = 0.5 ; y = - 0.5 ; z = 0.0  ! y-down node at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 3, :) = CMat_Row ( :)
      BVec ( 3) = bfx_fc_nodes ( ixx, iyy, izz, 3)

! ----------

      x = 0.5 ; y = 0.0 ; z = 0.5  ! z-top node at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 4, :) = CMat_Row ( :)
      BVec ( 4) = bfx_fc_nodes ( ixx, iyy, izz, 4)

! ----------

      x = 0.5 ; y = 0.0 ; z = - 0.5  ! z-bottom at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 5, :) = CMat_Row ( :)
      BVec ( 5) = bfx_fc_nodes ( ixx, iyy, izz, 5)

! --------------------
! We require matching 5 nodal values at the left x-face in LSQ sense.

      x = - 0.5 ; y = 0.0 ; z = 0.0  ! Central node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 6, :) = CMat_Row ( :)
      BVec ( 6) = bfx_fc_nodes ( ixx - 1, iyy, izz, 1)

! ----------

      x = - 0.5 ; y = 0.5 ; z = 0.0  ! y-up node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 7, :) = CMat_Row ( :)
      BVec ( 7) = bfx_fc_nodes ( ixx - 1, iyy, izz, 2)

! ----------

      x = - 0.5 ; y = - 0.5 ; z = 0.0  ! y-down node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 8, :) = CMat_Row ( :)
      BVec ( 8) = bfx_fc_nodes ( ixx - 1, iyy, izz, 3)

! ----------

      x = - 0.5 ; y = 0.0 ; z = 0.5  ! z-top node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 9, :) = CMat_Row ( :)
      BVec ( 9) = bfx_fc_nodes ( ixx - 1, iyy, izz, 4)

! ----------

      x = - 0.5 ; y = 0.0 ; z = - 0.5  ! z-bottom node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 10, :) = CMat_Row ( :)
      BVec ( 10) = bfx_fc_nodes ( ixx - 1, iyy, izz, 5)

! --------------------
! We require matching 5 nodal values at the up y-face in LSQ sense.

      x = 0.0 ; y = 0.5 ; z = 0.0  ! Central node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 11, :) = CMat_Row ( :)
      BVec ( 11) = bfy_fc_nodes ( ixx, iyy, izz, 1)

! ----------

      x = 0.5 ; y = 0.5 ; z = 0.0  ! x-right node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 12, :) = CMat_Row ( :)
      BVec ( 12) = bfy_fc_nodes ( ixx, iyy, izz, 2)

! ----------

      x = - 0.5 ; y = 0.5 ; z = 0.0  ! x-left node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 13, :) = CMat_Row ( :)
      BVec ( 13) = bfy_fc_nodes ( ixx, iyy, izz, 3)

! ----------

      x = 0.0 ; y = 0.5 ; z = 0.5  ! z-top node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 14, :) = CMat_Row ( :)
      BVec ( 14) = bfy_fc_nodes ( ixx, iyy, izz, 4)

! ----------

      x = 0.0 ; y = 0.5 ; z = - 0.5  ! z-bottom node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 15, :) = CMat_Row ( :)
      BVec ( 15) = bfy_fc_nodes ( ixx, iyy, izz, 5)

! --------------------
! We require matching 5 nodal values at the down y-face in LSQ sense.

      x = 0.0 ; y = - 0.5 ; z = 0.0  ! Central node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)
      
      AMat ( 16, :) = CMat_Row ( :)
      BVec ( 16) = bfy_fc_nodes ( ixx, iyy - 1, izz, 1)

! ----------

      x = 0.5 ; y = - 0.5 ; z = 0.0  ! x-right node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 17, :) = CMat_Row ( :)
      BVec ( 17) = bfy_fc_nodes ( ixx, iyy - 1, izz, 2)

! ----------

      x = - 0.5 ; y = - 0.5 ; z = 0.0  ! x-left node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 18, :) = CMat_Row ( :)
      BVec ( 18) = bfy_fc_nodes ( ixx, iyy - 1, izz, 3)

! ----------

      x = 0.0 ; y = - 0.5 ; z = 0.5  ! z-top node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 19, :) = CMat_Row ( :)
      BVec ( 19) = bfy_fc_nodes ( ixx, iyy - 1, izz, 4)

! ----------

      x = 0.0 ; y = - 0.5 ; z = - 0.5  ! z-bottom node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 20, :) = CMat_Row ( :)
      BVec ( 20) = bfy_fc_nodes ( ixx, iyy - 1, izz, 5)

! --------------------
! We require matching 5 nodal values at the top z-face in LSQ sense.

      x = 0.0 ; y = 0.0 ; z = 0.5  ! Central node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 21, :) = CMat_Row ( :)
      BVec ( 21) = bfz_fc_nodes ( ixx, iyy, izz, 1)

! ----------
      
      x = 0.5 ; y = 0.0 ; z = 0.5  ! x-right node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 22, :) = CMat_Row ( :)
      BVec ( 22) = bfz_fc_nodes ( ixx, iyy, izz, 2)

! ----------

      x = - 0.5 ; y = 0.0 ; z = 0.5  ! x-left node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 23, :) = CMat_Row ( :)
      BVec ( 23) = bfz_fc_nodes ( ixx, iyy, izz, 3)

! ----------

      x = 0.0 ; y = 0.5 ; z = 0.5  ! y-up node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 24, :) = CMat_Row ( :)
      BVec ( 24) = bfz_fc_nodes ( ixx, iyy, izz, 4)

! ----------

      x = 0.0 ; y = - 0.5 ; z = 0.5  ! y-down node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 25, :) = CMat_Row ( :)
      BVec ( 25) = bfz_fc_nodes ( ixx, iyy, izz, 5)

! --------------------
! We require matching 5 nodal values at the bottom z-face in LSQ sense.

      x = 0.0 ; y = 0.0 ; z = - 0.5  ! Central node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 26, :) = CMat_Row ( :)
      BVec ( 26) = bfz_fc_nodes ( ixx, iyy, izz - 1, 1)

! ----------

      x = 0.5 ; y = 0.0 ; z = - 0.5  ! x-right node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 27, :) = CMat_Row ( :)
      BVec ( 27) = bfz_fc_nodes ( ixx, iyy, izz - 1, 2)

! ----------

      x = - 0.5 ; y = 0.0 ; z = - 0.5  ! x-left node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 28, :) = CMat_Row ( :)
      BVec ( 28) = bfz_fc_nodes ( ixx, iyy, izz - 1, 3)

! ----------

      x = 0.0 ; y = 0.5 ; z = - 0.5  ! y-up node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 29, :) = CMat_Row ( :)
      BVec ( 29) = bfz_fc_nodes ( ixx, iyy, izz - 1, 4)

! ----------

      x = 0.0 ; y = - 0.5 ; z = - 0.5  ! y-down node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 30, :) = CMat_Row ( :)
      BVec ( 30) = bfz_fc_nodes ( ixx, iyy, izz - 1, 5)

! ----------
! Now in LSQ sense match the zone-centered modes of "u" that correspond to bfx.

      AMat ( 31, 1) = 1.0
      BVec ( 31) = u ( ixx, iyy, izz, 1, 1)

      AMat ( 32, 2) = 1.0
      BVec ( 32) = u ( ixx, iyy, izz, 1, 2)

      AMat ( 33, 3) = 1.0
      BVec ( 33) = u ( ixx, iyy, izz, 1, 3)

      AMat ( 34, 4) = 1.0
      BVec ( 34) = u ( ixx, iyy, izz, 1, 4)

      AMat ( 35, 5) = 1.0
      BVec ( 35) = u ( ixx, iyy, izz, 1, 5)

      AMat ( 36, 6) = 1.0
      BVec ( 36) = u ( ixx, iyy, izz, 1, 6)

      AMat ( 37, 7) = 1.0
      BVec ( 37) = u ( ixx, iyy, izz, 1, 7)

      AMat ( 38, 8) = 1.0
      BVec ( 38) = u ( ixx, iyy, izz, 1, 8)

      AMat ( 39, 9) = 1.0
      BVec ( 39) = u ( ixx, iyy, izz, 1, 9)

      AMat ( 40, 10) = 1.0
      BVec ( 40) = u ( ixx, iyy, izz, 1, 10)

! ----------
! Now in LSQ sense match the zone-centered modes of "u" that correspond to bfy.

      AMat ( 41, 11) = 1.0
      BVec ( 41) = u ( ixx, iyy, izz, 2, 1)

      AMat ( 42, 12) = 1.0
      BVec ( 42) = u ( ixx, iyy, izz, 2, 2)

      AMat ( 43, 13) = 1.0
      BVec ( 43) = u ( ixx, iyy, izz, 2, 3)

      AMat ( 44, 14) = 1.0
      BVec ( 44) = u ( ixx, iyy, izz, 2, 4)

      AMat ( 45, 15) = 1.0
      BVec ( 45) = u ( ixx, iyy, izz, 2, 5)

      AMat ( 46, 16) = 1.0
      BVec ( 46) = u ( ixx, iyy, izz, 2, 6)

      AMat ( 47, 17) = 1.0
      BVec ( 47) = u ( ixx, iyy, izz, 2, 7)

      AMat ( 48, 18) = 1.0
      BVec ( 48) = u ( ixx, iyy, izz, 2, 8)

      AMat ( 49, 19) = 1.0
      BVec ( 49) = u ( ixx, iyy, izz, 2, 9)

      AMat ( 50, 20) = 1.0
      BVec ( 50) = u ( ixx, iyy, izz, 2, 10)


! ----------
! Now in LSQ sense match the zone-centered modes of "u" that correspond to bfz.

      AMat ( 51, 21) = 1.0
      BVec ( 51) = u ( ixx, iyy, izz, 3, 1)

      AMat ( 52, 22) = 1.0
      BVec ( 52) = u ( ixx, iyy, izz, 3, 2)

      AMat ( 53, 23) = 1.0
      BVec ( 53) = u ( ixx, iyy, izz, 3, 3)

      AMat ( 54, 24) = 1.0
      BVec ( 54) = u ( ixx, iyy, izz, 3, 4)

      AMat ( 55, 25) = 1.0
      BVec ( 55) = u ( ixx, iyy, izz, 3, 5)

      AMat ( 56, 26) = 1.0
      BVec ( 56) = u ( ixx, iyy, izz, 3, 6)

      AMat ( 57, 27) = 1.0
      BVec ( 57) = u ( ixx, iyy, izz, 3, 7)

      AMat ( 58, 28) = 1.0
      BVec ( 58) = u ( ixx, iyy, izz, 3, 8)

      AMat ( 59, 29) = 1.0
      BVec ( 59) = u ( ixx, iyy, izz, 3, 9)

      AMat ( 60, 30) = 1.0
      BVec ( 60) = u ( ixx, iyy, izz, 3, 10)

! ----------------------------------------
! Build transpose of "AMat", "ATrans ( NUnknown, NLSQ)". From that, build
! "ATransA ( NUnknown, NUnknown)" and "ATransB ( NUnknown)". Also build
! "CTrans ( NUnknown, NConstraint)".

      ATrans = TRANSPOSE ( AMat)

      ATransA = MATMUL ( ATrans, AMat)

      ATransB = MATMUL ( ATrans, BVec)

      CTrans = TRANSPOSE ( CMat)

! ----------------------------------------
! Assemble the matrix "CLSQMat ( NUnknown + NConstraint, NUnknown + NConstraint)"
! and its right hand side "CLSQVec ( NUnknown + NConstraint)".

      CLSQMat ( 1: NUnknown, 1: NUnknown) =
     1                            ATransA ( 1: NUnknown, 1: NUnknown)

      CLSQMat ( 1: NUnknown, NUnknown + 1: NUnknown + NConstraint) =
     1                            CTrans ( 1: NUnknown, 1: NConstraint)

      CLSQMat ( NUnknown + 1: NUnknown + NConstraint, 1: NUnknown) =
     1                            CMat ( 1: NConstraint, 1: NUnknown)

      CLSQMat ( NUnknown + 1: NUnknown + NConstraint,
     1          NUnknown + 1: NUnknown + NConstraint) = 0.0

! ----------

      CLSQVec ( 1: NUnknown) = ATransB ( 1: NUnknown)

      CLSQVec ( NUnknown + 1: NUnknown + NConstraint) =
     1                             DVec ( 1: NConstraint)

! ----------------------------------------
! Solve the system "CLSQMat . (x, z)^T = CLSQVec". Extract from that the vector of
! unknowns.

      CALL MATSOLVE ( NUnknown + NConstraint, CLSQMat, CLSQVec,
     1                indx_CLSQ, vv_CLSQ)

! ----------------------------------------
! Now hand back the solution to the zone-centered reconstruction of the magnetic
! field that is closer to being divergence-free.

      bfx_zc_gr ( ixx, iyy, izz, 1: 10) = CLSQVec ( 1: 10)

      bfy_zc_gr ( ixx, iyy, izz, 1: 10) = CLSQVec ( 11: 20)

      bfz_zc_gr ( ixx, iyy, izz, 1: 10) = CLSQVec ( 21: 30)

! ----------------------------------------

          END DO

        END DO

      END DO

#endif
! End of "(NUMORDER == 3)".
! ----------------------------------------------------------------------
#if (NUMORDER == 4)
! STEP III.c) At fourth order, we use the 6 face-centered values at the 6 faces
! as constraints. We also impose the linear and quadratic variations of the divergence.
! We then want all the modes to be as close to the 30 volumetrically reconstructed
! modes in a least squares sense.
! "NConstraint = 10, NLSQ = 90, NUnknown = 60".

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! ----------------------------------------

          CMat = 0.0  ;  AMat = 0.0  ;  DVec = 0.0  ;  BVec = 0.0

! ----------------------------------------
! Load up the matrix of constraints "CMat ( NConstraint, NUnknown)" and its
! right hand side "DVec ( NConstraint)".

! --------------------
! Constant divergence-free constraint.

      CMat ( 1, 2) = 1.0 ; CMat ( 1, 23) = 1.0 ; CMat ( 1, 44) = 1.0
      CMat ( 1, 11) = - 3.0 / 20.0 ; CMat ( 1, 32) = - 3.0 / 20.0 
      CMat ( 1, 53) = - 3.0 / 20.0
      DVec ( 1) = 0.0

! --------------------
! Linear in x value of divergence-free constraint. "2 axx + bxy + cxz == 0"

      CMat ( 2, 5) = 2.0 ; CMat ( 2, 28) = 1.0 ; CMat ( 2, 50) = 1.0
      DVec ( 3) = 0.0

! ----------
! Linear in y value of divergence-free constraint. "axy + 2 byy + cyz == 0"

      CMat ( 3, 8) = 1.0 ; CMat ( 3, 26) = 2.0 ; CMat ( 3, 49) = 1.0
      DVec ( 3) = 0.0

! ----------
! Linear in z value of divergence-free constraint. "axz + byz + 2 czz == 0"

      CMat ( 4, 10) = 1.0 ; CMat ( 4, 29) = 1.0 ; CMat ( 4, 47) = 2.0
      DVec ( 4) = 0.0

! --------------------
! Quadratic in x^2 value of divergence-free constraint. "3 axxx + bxxy + cxxz == 0"

      CMat ( 5, 11) = 3.0 ; CMat ( 5, 34) = 1.0 ; CMat ( 5, 55) = 1.0
      DVec ( 5) = 0.0

! --------------------
! Quadratic in y^2 value of divergence-free constraint. "axyy + 3 byyy + cyyz == 0"

      CMat ( 6, 16) = 1.0 ; CMat ( 6, 32) = 3.0 ; CMat ( 6, 57) = 1.0
      DVec ( 6) = 0.0

! --------------------
! Quadratic in z^2 value of divergence-free constraint. "axzz + byzz + 3 czzz == 0"

      CMat ( 7, 18) = 1.0 ; CMat ( 7, 39) = 1.0 ; CMat ( 7, 53) = 3.0
      DVec ( 7) = 0.0

! --------------------
! Quadratic in x * y value of divergence-free constraint. "2 axxy + 2 bxyy + cxyz == 0"

      CMat ( 8, 14) = 2.0 ; CMat ( 8, 36) = 2.0 ; CMat ( 8, 60) = 1.0
      DVec ( 8) = 0.0

! --------------------
! Quadratic in y * z value of divergence-free constraint. "axyz + 2 byyz + 2 cyzz == 0"

      CMat ( 9, 20) = 1.0 ; CMat ( 9, 37) = 2.0 ; CMat ( 9, 59) = 2.0
      DVec ( 9) = 0.0

! --------------------
! Quadratic in x * z value of divergence-free constraint. "2 axxz + bxyz + 2 cxzz == 0"

      CMat ( 10, 15) = 2.0 ; CMat ( 10, 40) = 1.0 ; CMat ( 10, 58) = 2.0
      DVec ( 10) = 0.0

! ----------------------------------------
! Load up the matrix of LSQ equations "AMat ( NLSQ, NUnknown)" and their
! right hand side "BVec ( NLSQ)".

! --------------------
! We require matching 5 nodal values at the right x-face in LSQ sense.

! --------------------

      x = 0.5 ; y = 0.0 ; z = 0.0  ! Central node at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 1, :) = CMat_Row ( :)
      BVec ( 1) = bfx_fc_nodes ( ixx, iyy, izz, 1)

! ----------

      x = 0.5 ; y = 0.5 ; z = 0.0  ! y-up node at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 2, :) = CMat_Row ( :)
      BVec ( 2) = bfx_fc_nodes ( ixx, iyy, izz, 2)

! ----------

      x = 0.5 ; y = - 0.5 ; z = 0.0  ! y-down node at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 3, :) = CMat_Row ( :)
      BVec ( 3) = bfx_fc_nodes ( ixx, iyy, izz, 3)

! ----------

      x = 0.5 ; y = 0.0 ; z = 0.5  ! z-top node at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 4, :) = CMat_Row ( :)
      BVec ( 4) = bfx_fc_nodes ( ixx, iyy, izz, 4)

! ----------

      x = 0.5 ; y = 0.0 ; z = - 0.5  ! z-bottom at right x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 5, :) = CMat_Row ( :)
      BVec ( 5) = bfx_fc_nodes ( ixx, iyy, izz, 5)

! --------------------
! We require matching 5 nodal values at the left x-face in LSQ sense.

      x = - 0.5 ; y = 0.0 ; z = 0.0  ! Central node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 6, :) = CMat_Row ( :)
      BVec ( 6) = bfx_fc_nodes ( ixx - 1, iyy, izz, 1)

! ----------

      x = - 0.5 ; y = 0.5 ; z = 0.0  ! y-up node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 7, :) = CMat_Row ( :)
      BVec ( 7) = bfx_fc_nodes ( ixx - 1, iyy, izz, 2)

! ----------

      x = - 0.5 ; y = - 0.5 ; z = 0.0  ! y-down node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 8, :) = CMat_Row ( :)
      BVec ( 8) = bfx_fc_nodes ( ixx - 1, iyy, izz, 3)

! ----------

      x = - 0.5 ; y = 0.0 ; z = 0.5  ! z-top node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 9, :) = CMat_Row ( :)
      BVec ( 9) = bfx_fc_nodes ( ixx - 1, iyy, izz, 4)

! ----------

      x = - 0.5 ; y = 0.0 ; z = - 0.5  ! z-bottom node at left x-face.
      nx = 1.0 ; ny = 0.0 ; nz = 0.0  ! Normal to the x-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 10, :) = CMat_Row ( :)
      BVec ( 10) = bfx_fc_nodes ( ixx - 1, iyy, izz, 5)

! --------------------
! We require matching 5 nodal values at the up y-face in LSQ sense.

      x = 0.0 ; y = 0.5 ; z = 0.0  ! Central node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 11, :) = CMat_Row ( :)
      BVec ( 11) = bfy_fc_nodes ( ixx, iyy, izz, 1)

! ----------

      x = 0.5 ; y = 0.5 ; z = 0.0  ! x-right node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 12, :) = CMat_Row ( :)
      BVec ( 12) = bfy_fc_nodes ( ixx, iyy, izz, 2)

! ----------

      x = - 0.5 ; y = 0.5 ; z = 0.0  ! x-left node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 13, :) = CMat_Row ( :)
      BVec ( 13) = bfy_fc_nodes ( ixx, iyy, izz, 3)

! ----------

      x = 0.0 ; y = 0.5 ; z = 0.5  ! z-top node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 14, :) = CMat_Row ( :)
      BVec ( 14) = bfy_fc_nodes ( ixx, iyy, izz, 4)

! ----------

      x = 0.0 ; y = 0.5 ; z = - 0.5  ! z-bottom node at up y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 15, :) = CMat_Row ( :)
      BVec ( 15) = bfy_fc_nodes ( ixx, iyy, izz, 5)


! --------------------
! We require matching 5 nodal values at the down y-face in LSQ sense.

      x = 0.0 ; y = - 0.5 ; z = 0.0  ! Central node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 16, :) = CMat_Row ( :)
      BVec ( 16) = bfy_fc_nodes ( ixx, iyy - 1, izz, 1)

! ----------

      x = 0.5 ; y = - 0.5 ; z = 0.0  ! x-right node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 17, :) = CMat_Row ( :)
      BVec ( 17) = bfy_fc_nodes ( ixx, iyy - 1, izz, 2)

! ----------

      x = - 0.5 ; y = - 0.5 ; z = 0.0  ! x-left node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 18, :) = CMat_Row ( :)
      BVec ( 18) = bfy_fc_nodes ( ixx, iyy - 1, izz, 3)

! ----------

      x = 0.0 ; y = - 0.5 ; z = 0.5  ! z-top node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 19, :) = CMat_Row ( :)
      BVec ( 19) = bfy_fc_nodes ( ixx, iyy - 1, izz, 4)

! ----------

      x = 0.0 ; y = - 0.5 ; z = - 0.5  ! z-bottom node at down y-face.
      nx = 0.0 ; ny = 1.0 ; nz = 0.0  ! Normal to the y-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 20, :) = CMat_Row ( :)
      BVec ( 20) = bfy_fc_nodes ( ixx, iyy - 1, izz, 5)

! --------------------
! We require matching 5 nodal values at the top z-face in LSQ sense.

      x = 0.0 ; y = 0.0 ; z = 0.5  ! Central node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 21, :) = CMat_Row ( :)
      BVec ( 21) = bfz_fc_nodes ( ixx, iyy, izz, 1)

! ----------

      x = 0.5 ; y = 0.0 ; z = 0.5  ! x-right node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 22, :) = CMat_Row ( :)
      BVec ( 22) = bfz_fc_nodes ( ixx, iyy, izz, 2)

! ----------

      x = - 0.5 ; y = 0.0 ; z = 0.5  ! x-left node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 23, :) = CMat_Row ( :)
      BVec ( 23) = bfz_fc_nodes ( ixx, iyy, izz, 3)

! ----------

      x = 0.0 ; y = 0.5 ; z = 0.5  ! y-up node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 24, :) = CMat_Row ( :)
      BVec ( 24) = bfz_fc_nodes ( ixx, iyy, izz, 4)

! ----------

      x = 0.0 ; y = - 0.5 ; z = 0.5  ! y-down node at top z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 25, :) = CMat_Row ( :)
      BVec ( 25) = bfz_fc_nodes ( ixx, iyy, izz, 5)

! --------------------
! We require matching 5 nodal values at the bottom z-face in LSQ sense.

      x = 0.0 ; y = 0.0 ; z = - 0.5  ! Central node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 26, :) = CMat_Row ( :)
      BVec ( 26) = bfz_fc_nodes ( ixx, iyy, izz - 1, 1)

! ----------

      x = 0.5 ; y = 0.0 ; z = - 0.5  ! x-right node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 27, :) = CMat_Row ( :)
      BVec ( 27) = bfz_fc_nodes ( ixx, iyy, izz - 1, 2)

! ----------

      x = - 0.5 ; y = 0.0 ; z = - 0.5  ! x-left node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 28, :) = CMat_Row ( :)
      BVec ( 28) = bfz_fc_nodes ( ixx, iyy, izz - 1, 3)

! ----------

      x = 0.0 ; y = 0.5 ; z = - 0.5  ! y-up node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 29, :) = CMat_Row ( :)
      BVec ( 29) = bfz_fc_nodes ( ixx, iyy, izz - 1, 4)

! ----------

      x = 0.0 ; y = - 0.5 ; z = - 0.5  ! y-down node at bottom z-face.
      nx = 0.0 ; ny = 0.0 ; nz = 1.0  ! Normal to the z-face.

      CALL CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz, CMat_Row)

      AMat ( 30, :) = CMat_Row ( :)
      BVec ( 30) = bfz_fc_nodes ( ixx, iyy, izz - 1, 5)

! ----------
! Now in LSQ sense match the zone-centered modes of "u" that correspond to bfx.

      AMat ( 31, 1) = 1.0
      BVec ( 31) = u ( ixx, iyy, izz, 1, 1)

      AMat ( 32, 2) = 1.0
      BVec ( 32) = u ( ixx, iyy, izz, 1, 2)

      AMat ( 33, 3) = 1.0
      BVec ( 33) = u ( ixx, iyy, izz, 1, 3)

      AMat ( 34, 4) = 1.0
      BVec ( 34) = u ( ixx, iyy, izz, 1, 4)

      AMat ( 35, 5) = 1.0
      BVec ( 35) = u ( ixx, iyy, izz, 1, 5)

      AMat ( 36, 6) = 1.0
      BVec ( 36) = u ( ixx, iyy, izz, 1, 6)

      AMat ( 37, 7) = 1.0
      BVec ( 37) = u ( ixx, iyy, izz, 1, 7)

      AMat ( 38, 8) = 1.0
      BVec ( 38) = u ( ixx, iyy, izz, 1, 8)

      AMat ( 39, 9) = 1.0
      BVec ( 39) = u ( ixx, iyy, izz, 1, 9)

      AMat ( 40, 10) = 1.0
      BVec ( 40) = u ( ixx, iyy, izz, 1, 10)

      AMat ( 41, 11) = 1.0
      BVec ( 41) = u ( ixx, iyy, izz, 1, 11)

      AMat ( 42, 12) = 1.0
      BVec ( 42) = u ( ixx, iyy, izz, 1, 12)

      AMat ( 43, 13) = 1.0
      BVec ( 43) = u ( ixx, iyy, izz, 1, 13)

      AMat ( 44, 14) = 1.0
      BVec ( 44) = u ( ixx, iyy, izz, 1, 14)

      AMat ( 45, 15) = 1.0
      BVec ( 45) = u ( ixx, iyy, izz, 1, 15)

      AMat ( 46, 16) = 1.0
      BVec ( 46) = u ( ixx, iyy, izz, 1, 16)

      AMat ( 47, 17) = 1.0
      BVec ( 47) = u ( ixx, iyy, izz, 1, 17)

      AMat ( 48, 18) = 1.0
      BVec ( 48) = u ( ixx, iyy, izz, 1, 18)

      AMat ( 49, 19) = 1.0
      BVec ( 49) = u ( ixx, iyy, izz, 1, 19)

      AMat ( 50, 20) = 1.0
      BVec ( 50) = u ( ixx, iyy, izz, 1, 20)


! ----------
! Now in LSQ sense match the zone-centered modes of "u" that correspond to bfy.

      AMat ( 51, 21) = 1.0
      BVec ( 51) = u ( ixx, iyy, izz, 2, 1)

      AMat ( 52, 22) = 1.0
      BVec ( 52) = u ( ixx, iyy, izz, 2, 2)

      AMat ( 53, 23) = 1.0
      BVec ( 53) = u ( ixx, iyy, izz, 2, 3)

      AMat ( 54, 24) = 1.0
      BVec ( 54) = u ( ixx, iyy, izz, 2, 4)

      AMat ( 55, 25) = 1.0
      BVec ( 55) = u ( ixx, iyy, izz, 2, 5)

      AMat ( 56, 26) = 1.0
      BVec ( 56) = u ( ixx, iyy, izz, 2, 6)

      AMat ( 57, 27) = 1.0
      BVec ( 57) = u ( ixx, iyy, izz, 2, 7)

      AMat ( 58, 28) = 1.0
      BVec ( 58) = u ( ixx, iyy, izz, 2, 8)

      AMat ( 59, 29) = 1.0
      BVec ( 59) = u ( ixx, iyy, izz, 2, 9)

      AMat ( 60, 30) = 1.0
      BVec ( 60) = u ( ixx, iyy, izz, 2, 10)

      AMat ( 61, 31) = 1.0
      BVec ( 61) = u ( ixx, iyy, izz, 2, 11)

      AMat ( 62, 32) = 1.0
      BVec ( 62) = u ( ixx, iyy, izz, 2, 12)

      AMat ( 63, 33) = 1.0
      BVec ( 63) = u ( ixx, iyy, izz, 2, 13)

      AMat ( 64, 34) = 1.0
      BVec ( 64) = u ( ixx, iyy, izz, 2, 14)

      AMat ( 65, 35) = 1.0
      BVec ( 65) = u ( ixx, iyy, izz, 2, 15)

      AMat ( 66, 36) = 1.0
      BVec ( 66) = u ( ixx, iyy, izz, 2, 16)

      AMat ( 67, 37) = 1.0
      BVec ( 67) = u ( ixx, iyy, izz, 2, 17)

      AMat ( 68, 38) = 1.0
      BVec ( 68) = u ( ixx, iyy, izz, 2, 18)

      AMat ( 69, 39) = 1.0
      BVec ( 69) = u ( ixx, iyy, izz, 2, 19)

      AMat ( 70, 40) = 1.0
      BVec ( 70) = u ( ixx, iyy, izz, 2, 20)

    
! ----------
! Now in LSQ sense match the zone-centered modes of "u" that correspond to bfz.

      AMat ( 71, 41) = 1.0
      BVec ( 71) = u ( ixx, iyy, izz, 3, 1)

      AMat ( 72, 42) = 1.0
      BVec ( 72) = u ( ixx, iyy, izz, 3, 2)

      AMat ( 73, 43) = 1.0
      BVec ( 73) = u ( ixx, iyy, izz, 3, 3)

      AMat ( 74, 44) = 1.0
      BVec ( 74) = u ( ixx, iyy, izz, 3, 4)

      AMat ( 75, 45) = 1.0
      BVec ( 75) = u ( ixx, iyy, izz, 3, 5)

      AMat ( 76, 46) = 1.0
      BVec ( 76) = u ( ixx, iyy, izz, 3, 6)

      AMat ( 77, 47) = 1.0
      BVec ( 77) = u ( ixx, iyy, izz, 3, 7)

      AMat ( 78, 48) = 1.0
      BVec ( 78) = u ( ixx, iyy, izz, 3, 8)

      AMat ( 79, 49) = 1.0
      BVec ( 79) = u ( ixx, iyy, izz, 3, 9)

      AMat ( 80, 50) = 1.0
      BVec ( 80) = u ( ixx, iyy, izz, 3, 10)

      AMat ( 81, 51) = 1.0
      BVec ( 81) = u ( ixx, iyy, izz, 3, 11)

      AMat ( 82, 52) = 1.0
      BVec ( 82) = u ( ixx, iyy, izz, 3, 12)

      AMat ( 83, 53) = 1.0
      BVec ( 83) = u ( ixx, iyy, izz, 3, 13)

      AMat ( 84, 54) = 1.0
      BVec ( 84) = u ( ixx, iyy, izz, 3, 14)

      AMat ( 85, 55) = 1.0
      BVec ( 85) = u ( ixx, iyy, izz, 3, 15)

      AMat ( 86, 56) = 1.0
      BVec ( 86) = u ( ixx, iyy, izz, 3, 16)

      AMat ( 87, 57) = 1.0
      BVec ( 87) = u ( ixx, iyy, izz, 3, 17)

      AMat ( 88, 58) = 1.0
      BVec ( 88) = u ( ixx, iyy, izz, 3, 18)

      AMat ( 89, 59) = 1.0
      BVec ( 89) = u ( ixx, iyy, izz, 3, 19)

      AMat ( 90, 60) = 1.0
      BVec ( 90) = u ( ixx, iyy, izz, 3, 20)

! ----------------------------------------
! Build transpose of "AMat", "ATrans ( NUnknown, NLSQ)". From that, build
! "ATransA ( NUnknown, NUnknown)" and "ATransB ( NUnknown)". Also build
! "CTrans ( NUnknown, NConstraint)".

      ATrans = TRANSPOSE ( AMat)

      ATransA = MATMUL ( ATrans, AMat)

      ATransB = MATMUL ( ATrans, BVec)

      CTrans = TRANSPOSE ( CMat)

! ----------------------------------------
! Assemble the matrix "CLSQMat ( NUnknown + NConstraint, NUnknown + NConstraint)"
! and its right hand side "CLSQVec ( NUnknown + NConstraint)".

      CLSQMat ( 1: NUnknown, 1: NUnknown) =
     1                            ATransA ( 1: NUnknown, 1: NUnknown)

      CLSQMat ( 1: NUnknown, NUnknown + 1: NUnknown + NConstraint) =
     1                            CTrans ( 1: NUnknown, 1: NConstraint)

      CLSQMat ( NUnknown + 1: NUnknown + NConstraint, 1: NUnknown) =
     1                            CMat ( 1: NConstraint, 1: NUnknown)

      CLSQMat ( NUnknown + 1: NUnknown + NConstraint,
     1          NUnknown + 1: NUnknown + NConstraint) = 0.0

! ----------

      CLSQVec ( 1: NUnknown) = ATransB ( 1: NUnknown)

      CLSQVec ( NUnknown + 1: NUnknown + NConstraint) =
     1                             DVec ( 1: NConstraint)

! ----------------------------------------
! Solve the system "CLSQMat . (x, z)^T = CLSQVec". Extract from that the vector of
! unknowns.

      CALL MATSOLVE ( NUnknown + NConstraint, CLSQMat, CLSQVec,
     1                indx_CLSQ, vv_CLSQ)

! ----------------------------------------
! Now hand back the solution to the zone-centered reconstruction of the magnetic
! field that is closer to being divergence-free.

      bfx_zc_gr ( ixx, iyy, izz, 1: 20) = CLSQVec ( 1: 20)

      bfy_zc_gr ( ixx, iyy, izz, 1: 20) = CLSQVec ( 21: 40)

      bfz_zc_gr ( ixx, iyy, izz, 1: 20) = CLSQVec ( 41: 60)

! ----------------------------------------

          END DO

        END DO

      END DO

#endif
! End of "(NUMORDER == 4)".
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP IV) Now make the diagnostics for evaluating order property and the div-free
! aspect. Here the order property is evaluated at the zone centers.

      L1_Bfx = 0.0 ; Linf_Bfx = 0.0 ; L1_Bfy = 0.0 ; Linf_Bfy = 0.0
      L1_Bfz = 0.0 ; Linf_Bfz = 0.0

      L1_Div_lo = 0.0 ; Linf_Div_lo = 0.0

      L1_Div = 0.0 ; Linf_Div = 0.0 ; L1_u_Div = 0.0 ; Linf_u_Div = 0.0

! ----------------------------------------------------------------------

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------------------------------------

          x = 0.0  ;  y = 0.0  ;  z = 0.0

! ----------------------------------------

          u_pt_1 ( :) = bfx_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

          tempab = Bfx_pt ( xc_save ( ixx), yc_save ( iyy),
     1                                      zc_save ( izz))

          L1_Bfx = L1_Bfx + ABS ( tempab - tempaa)

          Linf_Bfx = AMAX1 ( Linf_Bfx, ABS ( tempab - tempaa) )

! ----------------------------------------

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

          tempab = Bfy_pt ( xc_save ( ixx), yc_save ( iyy), 
     1                                      zc_save ( izz))

          L1_Bfy = L1_Bfy + ABS ( tempab - tempaa)

          Linf_Bfy = AMAX1 ( Linf_Bfy, ABS ( tempab - tempaa) )

! ----------------------------------------

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

          tempab = Bfz_pt ( xc_save ( ixx), yc_save ( iyy), 
     1                                      zc_save ( izz))

          L1_Bfz = L1_Bfz + ABS ( tempab - tempaa)

          Linf_Bfz = AMAX1 ( Linf_Bfz, ABS ( tempab - tempaa) )

! ----------------------------------------
#if (NUMORDER == 2)

          tempaa = bfx_zc_gr ( ixx, iyy, izz, 2)
     1           + bfy_zc_gr ( ixx, iyy, izz, 3)
     1           + bfz_zc_gr ( ixx, iyy, izz, 4)

          L1_Div_lo = L1_Div_lo + ABS ( tempaa)
          Linf_Div_lo = AMAX1 ( Linf_Div_lo, ABS ( tempaa) )

          tempaa = u ( ixx, iyy, izz, 1, 2)
     1           + u ( ixx, iyy, izz, 2, 3)
     1           + u ( ixx, iyy, izz, 3, 4)

          L1_u_Div = L1_u_Div + ABS ( tempaa)
          Linf_u_Div = AMAX1 ( Linf_u_Div, ABS ( tempaa) )

#endif
! End of "(NUMORDER == 2)".
! ----------------------------------------
#if (NUMORDER == 3)

          x = 0.0  ;  y = 0.0  ;  z = 0.0

          tempaa = bfx_zc_gr ( ixx, iyy, izz, 2)
     1           + bfy_zc_gr ( ixx, iyy, izz, 3)
     1           + bfz_zc_gr ( ixx, iyy, izz, 4)

          tempab = 2.0 * bfx_zc_gr ( ixx, iyy, izz, 5)
     1           +       bfy_zc_gr ( ixx, iyy, izz, 8)
     1           +       bfz_zc_gr ( ixx, iyy, izz, 10)

          tempac =       bfx_zc_gr ( ixx, iyy, izz, 8)
     1           + 2.0 * bfy_zc_gr ( ixx, iyy, izz, 6)
     1           +       bfz_zc_gr ( ixx, iyy, izz, 9)

          tempad =       bfx_zc_gr ( ixx, iyy, izz, 10)
     1           +       bfy_zc_gr ( ixx, iyy, izz, 9)
     1           + 2.0 * bfz_zc_gr ( ixx, iyy, izz, 7)

          L1_Div_lo = L1_Div_lo + ABS ( tempaa)
          Linf_Div_lo = AMAX1 ( Linf_Div_lo, ABS ( tempaa) )

          L1_Div = L1_Div + ABS ( tempab)
     1                    + ABS ( tempac) + ABS ( tempad)
          Linf_Div = AMAX1 ( Linf_Div, ABS ( tempab),
     1                                 ABS ( tempac), ABS ( tempad) )

          tempaa = u ( ixx, iyy, izz, 1, 2)
     1           + u ( ixx, iyy, izz, 2, 3)
     1           + u ( ixx, iyy, izz, 3, 4)

          tempab = 2.0 * u ( ixx, iyy, izz, 1, 5)
     1           +       u ( ixx, iyy, izz, 2, 8)
     1           +       u ( ixx, iyy, izz, 3, 10)

          tempac =       u ( ixx, iyy, izz, 1, 8)
     1           + 2.0 * u ( ixx, iyy, izz, 2, 6)
     1           +       u ( ixx, iyy, izz, 3, 9)

          tempad =       u ( ixx, iyy, izz, 1, 10)
     1           +       u ( ixx, iyy, izz, 2, 9)
     1           + 2.0 * u ( ixx, iyy, izz, 3, 7)

          L1_u_Div = L1_u_Div + ABS ( tempab)
     1                    + ABS ( tempac) + ABS ( tempad)
          Linf_u_Div = AMAX1 ( Linf_u_Div, ABS ( tempab),
     1                    ABS ( tempac), ABS ( tempad) )

#endif
! End of "(NUMORDER == 3)".
! ----------------------------------------
#if (NUMORDER == 4)
    
         x = 0.0  ;  y = 0.0  ;  z = 0.0

          tempaa = bfx_zc_gr ( ixx, iyy, izz, 2)
     1           + bfy_zc_gr ( ixx, iyy, izz, 3)
     1           + bfz_zc_gr ( ixx, iyy, izz, 4)
     1           - ( 3.0 / 20.0) * bfx_zc_gr ( ixx, iyy, izz, 11)
     1           - ( 3.0 / 20.0) * bfy_zc_gr ( ixx, iyy, izz, 12)
     1           - ( 3.0 / 20.0) * bfz_zc_gr ( ixx, iyy, izz, 13)

          tempab = 2.0 * bfx_zc_gr ( ixx, iyy, izz, 5)
     1           +       bfy_zc_gr ( ixx, iyy, izz, 8)
     1           +       bfz_zc_gr ( ixx, iyy, izz, 10)

          tempac =       bfx_zc_gr ( ixx, iyy, izz, 8)
     1           + 2.0 * bfy_zc_gr ( ixx, iyy, izz, 6)
     1           +       bfz_zc_gr ( ixx, iyy, izz, 9)

          tempad =       bfx_zc_gr ( ixx, iyy, izz, 10)
     1           +       bfy_zc_gr ( ixx, iyy, izz, 9)
     1           + 2.0 * bfz_zc_gr ( ixx, iyy, izz, 7)

          tempae = 3.0 * bfx_zc_gr ( ixx, iyy, izz, 11)
     1           +       bfy_zc_gr ( ixx, iyy, izz, 14)
     1           +       bfz_zc_gr ( ixx, iyy, izz, 15)

          tempaf =       bfx_zc_gr ( ixx, iyy, izz, 16)
     1           + 3.0 * bfy_zc_gr ( ixx, iyy, izz, 12)
     1           +       bfz_zc_gr ( ixx, iyy, izz, 17)

          tempag =       bfx_zc_gr ( ixx, iyy, izz, 18)
     1           +       bfy_zc_gr ( ixx, iyy, izz, 19)
     1           + 3.0 * bfz_zc_gr ( ixx, iyy, izz, 13)


          tempah = 2.0 * bfx_zc_gr ( ixx, iyy, izz, 14)
     1           + 2.0 * bfy_zc_gr ( ixx, iyy, izz, 16)
     1           +       bfz_zc_gr ( ixx, iyy, izz, 20)

          tempai =       bfx_zc_gr ( ixx, iyy, izz, 20)
     1           + 2.0 * bfy_zc_gr ( ixx, iyy, izz, 17)
     1           + 2.0 * bfz_zc_gr ( ixx, iyy, izz, 19)

          tempaj = 2.0 * bfx_zc_gr ( ixx, iyy, izz, 15)
     1           +       bfy_zc_gr ( ixx, iyy, izz, 20)
     1           + 2.0 * bfz_zc_gr ( ixx, iyy, izz, 18)

          L1_Div_lo = L1_Div_lo + ABS ( tempaa)
          Linf_Div_lo = AMAX1 ( Linf_Div_lo, ABS ( tempaa) )

          L1_Div = L1_Div + ABS ( tempab) + ABS ( tempac)
     1                    + ABS ( tempad) + ABS ( tempae)
     1                    + ABS ( tempaf) + ABS ( tempag)
     1                    + ABS ( tempah) + ABS ( tempai)
     1                    + ABS ( tempaj)

          Linf_Div = AMAX1 ( Linf_Div, ABS ( tempab), ABS ( tempac),
     1                                 ABS ( tempad), ABS ( tempae), 
     1                                 ABS ( tempaf), ABS ( tempag),
     1                                 ABS ( tempah), ABS ( tempai),
     1                                 ABS ( tempaj) )

  
          tempaa = u ( ixx, iyy, izz, 1, 2)
     1           + u ( ixx, iyy, izz, 2, 3)
     1           + u ( ixx, iyy, izz, 3, 4)
     1           - ( 3.0 / 20.0) * u ( ixx, iyy, izz, 1, 11)
     1           - ( 3.0 / 20.0) * u ( ixx, iyy, izz, 2, 12)
     1           - ( 3.0 / 20.0) * u ( ixx, iyy, izz, 3, 13)


          tempab = 2.0 * u ( ixx, iyy, izz, 1, 5)
     1           +       u ( ixx, iyy, izz, 2, 8)
     1           +       u ( ixx, iyy, izz, 3, 10)

          tempac =       u ( ixx, iyy, izz, 1, 8)
     1           + 2.0 * u ( ixx, iyy, izz, 2, 6)
     1           +       u ( ixx, iyy, izz, 3, 9)

          tempad =       u ( ixx, iyy, izz, 1, 10)
     1           +       u ( ixx, iyy, izz, 2, 9)
     1           + 2.0 * u ( ixx, iyy, izz, 3, 7)

          tempae = 3.0 * u ( ixx, iyy, izz, 1, 11)
     1           +       u ( ixx, iyy, izz, 2, 14)
     1           +       u ( ixx, iyy, izz, 3, 15)

          tempaf =       u ( ixx, iyy, izz, 1, 16)
     1           + 3.0 * u ( ixx, iyy, izz, 2, 12)
     1           +       u ( ixx, iyy, izz, 3, 17)

          tempag =       u ( ixx, iyy, izz, 1, 18)
     1           +       u ( ixx, iyy, izz, 2, 19)
     1           + 3.0 * u ( ixx, iyy, izz, 3, 13)


          tempah = 2.0 * u ( ixx, iyy, izz, 1, 14)
     1           + 2.0 * u ( ixx, iyy, izz, 2, 16)
     1           +       u ( ixx, iyy, izz, 3, 20)

          tempai =       u ( ixx, iyy, izz, 1, 20)
     1           + 2.0 * u ( ixx, iyy, izz, 2, 17)
     1           + 2.0 * u ( ixx, iyy, izz, 3, 19)

          tempaj = 2.0 * u ( ixx, iyy, izz, 1, 15)
     1           +       u ( ixx, iyy, izz, 2, 20)
     1           + 2.0 * u ( ixx, iyy, izz, 3, 18)

          L1_u_Div = L1_u_Div + ABS ( tempab) + ABS ( tempac)
     1                        + ABS ( tempad) + ABS ( tempae)
     1                        + ABS ( tempaf) + ABS ( tempag)
     1                        + ABS ( tempah) + ABS ( tempai)
     1                        + ABS ( tempaj)     

          Linf_u_Div = AMAX1 ( Linf_u_Div, ABS ( tempab), ABS ( tempac),
     1                                     ABS ( tempad), ABS ( tempae),
     1                                     ABS ( tempaf), ABS ( tempag),
     1                                     ABS ( tempah), ABS ( tempai),
     1                                     ABS ( tempaj) )



#endif
! End of "(NUMORDER == 4)".
! ----------------------------------------

          END DO

        END DO

      END DO

! ----------------------------------------------------------------------

      L1_Bfx = L1_Bfx / REAL ( ix1 * iy1 * iz1)
      L1_Bfy = L1_Bfy / REAL ( ix1 * iy1 * iz1)
      L1_Bfz = L1_Bfz / REAL ( ix1 * iy1 * iz1)

      L1_Div = L1_Div / REAL ( ix1 * iy1 * iz1)
      L1_Div_lo = L1_Div_lo / REAL ( ix1 * iy1 * iz1)
      L1_u_Div = L1_u_Div / REAL ( ix1 * iy1 * iz1)

      WRITE ( 6, *)" "
      WRITE ( 6, *)" L1_Bfx, Linf_Bfx, L1_Bfy, Linf_Bfy = ",
     1               L1_Bfx, Linf_Bfx, L1_Bfy, Linf_Bfy
      WRITE ( 6, *)" "
      WRITE ( 6, *)" L1_Bfz, Linf_Bfz, L1_Div, Linf_Div = ",
     1               L1_Bfz, Linf_Bfz, L1_Div, Linf_Div
      WRITE ( 6, *)" "
      WRITE ( 6, *)" L1_Div_lo, Linf_Div_lo = ",
     1               L1_Div_lo, Linf_Div_lo
      WRITE ( 6, *)" "
      WRITE ( 6, *)" L1_u_Div, Linf_u_Div = ", L1_u_Div, Linf_u_Div

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP V) Now make the diagnostics for evaluating order property.
! Here the order property is evaluated at the nodal points from either
! of the two sides of the faces.

      L1_Bfx = 0.0 ; Linf_Bfx = 0.0 ; L1_Bfy = 0.0 ; Linf_Bfy = 0.0
      L1_Bfz = 0.0 ; Linf_Bfz = 0.0

! ----------------------------------------

      DO izz = 1, iz1

        DO iyy = 1, iy1

          DO ixx = 0, ix1

! ----------------------------------------

          x = 0.5  ;  y = 0.0  ;  z = 0.0

          u_pt_1 ( :) = bfx_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = -0.5  ;  y = 0.0  ;  z = 0.0

          u_pt_1 ( :) = bfx_zc_gr ( ixx + 1, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xb_save ( ixx)
          y = yc_save ( iyy)
          z = zc_save ( izz)

          tempac = Bfx_pt ( x, y, z)

! ----------

          L1_Bfx = L1_Bfx + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfx = AMAX1 ( Linf_Bfx, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = 0.5  ;  y = 0.5  ;  z = 0.0

          u_pt_1 ( :) = bfx_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = -0.5  ;  y = 0.5  ;  z = 0.0

          u_pt_1 ( :) = bfx_zc_gr ( ixx + 1, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xb_save ( ixx)
          y = yc_save ( iyy) + 0.5 * dy_save ( iyy)
          z = zc_save ( izz)

          tempac = Bfx_pt ( x, y, z)

! ----------

          L1_Bfx = L1_Bfx + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfx = AMAX1 ( Linf_Bfx, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = 0.5  ;  y = - 0.5  ;  z = 0.0

          u_pt_1 ( :) = bfx_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = -0.5  ;  y = - 0.5  ;  z = 0.0

          u_pt_1 ( :) = bfx_zc_gr ( ixx + 1, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xb_save ( ixx)
          y = yc_save ( iyy) - 0.5 * dy_save ( iyy)
          z = zc_save ( izz)

          tempac = Bfx_pt ( x, y, z)

! ----------

          L1_Bfx = L1_Bfx + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfx = AMAX1 ( Linf_Bfx, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = 0.5  ;  y = 0.0  ;  z = 0.5

          u_pt_1 ( :) = bfx_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = -0.5  ;  y = 0.0  ;  z = 0.5

          u_pt_1 ( :) = bfx_zc_gr ( ixx + 1, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xb_save ( ixx)
          y = yc_save ( iyy)
          z = zc_save ( izz) + 0.5 * dz_save ( izz)

          tempac = Bfx_pt ( x, y, z)

! ----------

          L1_Bfx = L1_Bfx + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfx = AMAX1 ( Linf_Bfx, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = 0.5  ;  y = 0.0  ;  z = - 0.5

          u_pt_1 ( :) = bfx_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,   
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = -0.5  ;  y = 0.0  ;  z = - 0.5

          u_pt_1 ( :) = bfx_zc_gr ( ixx + 1, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xb_save ( ixx)
          y = yc_save ( iyy)
          z = zc_save ( izz) - 0.5 * dz_save ( izz)

          tempac = Bfx_pt ( x, y, z)

! ----------

          L1_Bfx = L1_Bfx + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac) 

          Linf_Bfx = AMAX1 ( Linf_Bfx, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          END DO

        END DO

      END DO

! ----------

      L1_Bfx = L1_Bfx / ( 10.0 * (ix1 + 1.0) * iy1 * iz1)
 
      WRITE ( 6, *)" At x-faces:"
      WRITE ( 6, *)" L1_Bfx, Linf_Bfx = ",
     1               L1_Bfx, Linf_Bfx
      WRITE ( 6, *)" "

! --------------------------------------------------

      DO izz = 1, iz1

        DO iyy = 0, iy1

          DO ixx = 1, ix1

! ----------------------------------------

          x = 0.0  ;  y = 0.5  ;  z = 0.0

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = 0.0  ;  y = - 0.5  ;  z = 0.0

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy + 1, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xc_save ( ixx)
          y = yb_save ( iyy)
          z = zc_save ( izz)

          tempac = Bfy_pt ( x, y, z)

! ----------

          L1_Bfy = L1_Bfy + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfy = AMAX1 ( Linf_Bfy, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = 0.5  ;  y = 0.5  ;  z = 0.0

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = 0.5  ;  y = - 0.5  ;  z = 0.0

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy + 1, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xc_save ( ixx) + 0.5 * dx_save ( ixx)
          y = yb_save ( iyy)
          z = zc_save ( izz)

          tempac = Bfy_pt ( x, y, z)

! ----------

          L1_Bfy = L1_Bfy + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfy = AMAX1 ( Linf_Bfy, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = - 0.5  ;  y = 0.5  ;  z = 0.0

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = - 0.5  ;  y = - 0.5  ;  z = 0.0

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy + 1, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xc_save ( ixx) - 0.5 * dx_save ( ixx)
          y = yb_save ( iyy)
          z = zc_save ( izz)

          tempac = Bfy_pt ( x, y, z)

! ----------

          L1_Bfy = L1_Bfy + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfy = AMAX1 ( Linf_Bfy, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = 0.0  ;  y = 0.5  ;  z = 0.5

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = 0.0  ;  y = - 0.5  ;  z = 0.5

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy + 1, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xc_save ( ixx)
          y = yb_save ( iyy)
          z = zc_save ( izz) + 0.5 * dz_save ( izz)

          tempac = Bfy_pt ( x, y, z)

! ----------

          L1_Bfy = L1_Bfy + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfy = AMAX1 ( Linf_Bfy, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = 0.0  ;  y = 0.5  ;  z = - 0.5

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = 0.0  ;  y = - 0.5  ;  z = - 0.5

          u_pt_1 ( :) = bfy_zc_gr ( ixx, iyy + 1, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xc_save ( ixx)
          y = yb_save ( iyy)
          z = zc_save ( izz) - 0.5 * dz_save ( izz)

          tempac = Bfy_pt ( x, y, z)

! ----------

          L1_Bfy = L1_Bfy + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfy = AMAX1 ( Linf_Bfy, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          END DO

        END DO

      END DO

! ----------

      L1_Bfy = L1_Bfy / ( 10.0 * ix1 * (iy1 + 1.0) * iz1)

      WRITE ( 6, *)" At y-faces:"
      WRITE ( 6, *)" L1_Bfy, Linf_Bfy = ",
     1               L1_Bfy, Linf_Bfy
      WRITE ( 6, *)" "

! --------------------------------------------------

      DO izz = 0, iz1

        DO iyy = 1, iy1

          DO ixx = 1, ix1

! ----------------------------------------

          x = 0.0  ;  y = 0.0  ;  z = 0.5

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = 0.0  ;  y = 0.0  ;  z = - 0.5

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz + 1, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xc_save ( ixx)
          y = yc_save ( iyy)
          z = zb_save ( izz)

          tempac = Bfz_pt ( x, y, z)

! ----------

          L1_Bfz = L1_Bfz + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfz = AMAX1 ( Linf_Bfz, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = 0.5  ;  y = 0.0  ;  z = 0.5

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = 0.5  ;  y = 0.0  ;  z = - 0.5

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz + 1, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xc_save ( ixx) + 0.5 * dx_save ( ixx)
          y = yc_save ( iyy)
          z = zb_save ( izz)

          tempac = Bfz_pt ( x, y, z)

! ----------

          L1_Bfz = L1_Bfz + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfz = AMAX1 ( Linf_Bfz, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = - 0.5  ;  y = 0.0  ;  z = 0.5

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = - 0.5  ;  y = 0.0  ;  z = - 0.5

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz + 1, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xc_save ( ixx) - 0.5 * dx_save ( ixx)
          y = yc_save ( iyy)
          z = zb_save ( izz)

          tempac = Bfz_pt ( x, y, z)

! ----------

          L1_Bfz = L1_Bfz + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfz = AMAX1 ( Linf_Bfz, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = 0.0  ;  y = 0.5  ;  z = 0.5

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = 0.0  ;  y = 0.5  ;  z = - 0.5

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz + 1, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xc_save ( ixx)
          y = yc_save ( iyy) + 0.5 * dy_save ( iyy)
          z = zb_save ( izz)

          tempac = Bfz_pt ( x, y, z)

! ----------

          L1_Bfz = L1_Bfz + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfz = AMAX1 ( Linf_Bfz, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          x = 0.0  ;  y = - 0.5  ;  z = 0.5

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempaa)

! ----------

          x = 0.0  ;  y = - 0.5  ;  z = - 0.5

          u_pt_1 ( :) = bfz_zc_gr ( ixx, iyy, izz + 1, :)

          CALL v_polynomial_modes ( x, y, z,
     1                              n_cc_modes, u_pt_1, tempab)

! ----------

          x = xc_save ( ixx)
          y = yc_save ( iyy) - 0.5 * dy_save ( iyy)
          z = zb_save ( izz)

          tempac = Bfz_pt ( x, y, z)

! ----------

          L1_Bfz = L1_Bfz + ABS ( tempaa - tempac)
     1                    + ABS ( tempab - tempac)

          Linf_Bfz = AMAX1 ( Linf_Bfz, ABS ( tempaa - tempac),
     1                                 ABS ( tempab - tempac))

! ----------------------------------------

          END DO

        END DO

      END DO

! ----------

      L1_Bfz = L1_Bfz / ( 10.0 * ix1 * iy1 * (iz1 + 1.0))

      WRITE ( 6, *)" At z-faces:"
      WRITE ( 6, *)" L1_Bfz, Linf_Bfz = ",
     1               L1_Bfz, Linf_Bfz
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END PROGRAM TEST_HiOrder_DivPreserving_Vect_Field

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if (NUMORDER == 2)

      SUBROUTINE LIMIT_RECONST_ALL_ZC_P1

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_modes, n_fc_modes, n_cc_efld_modes, n_cc_components,
     1        num_level, num_grid, num_pe,
     1        n_phys_consts, phys_consts,

     1        u, flattengr, lo_order_zone_avg,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Uses "efxbdy, efybdy, efzbdy to build "sourcegr". By construction, this
! "sourcegr" already has the factor of "4 * pi" built in.

! The higher moments of "sourcegr" are then built using third order
! WENO reconstruction.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_modes, n_fc_modes, n_cc_efld_modes, n_cc_components,
     1        num_level, num_grid, num_pe, n_phys_consts

      REAL phys_consts ( n_phys_consts)

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) :: lo_order_zone_avg

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::

     1      flattengr

! ----------

! Arrays that specify geometry.

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

      INTEGER ixx, iyy, izz, ixx0, iyy0, izz0,
     1        ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,
     1        ixxp2, iyyp2, izzp2, ixxm2, iyym2, izzm2,
     1        ixxp3, iyyp3, izzp3, ixxm3, iyym3, izzm3,
     1        iprim, ifail, idirection, ii

      REAL c, r13_3, central_zone_wt, MC_beta, x, y, smallnum, mc_coef

      REAL, DIMENSION ( n_cc_components) :: a, b,
     1     u_ip2, u_ip1, u_i, u_im1, u_im2,

     1     R1_x, R2_x, R3_x, R1_xx, R2_xx, R3_xx,

     1     R1_xy, R2_xy, R3_xy, R4_xy,

     1     Is1, Is2, Is3, Is4,
     1     Wt1, Wt2, Wt3, Wt4,
     1     u_pt, v_pt, Eigenval_pt

      REAL, DIMENSION ( n_cc_components, n_cc_components) ::
     1      Lf_Eigenvec_pt, Rt_Eigenvec_pt

! MC_beta limiter : mc_coef = 1.0 ==> MinMod; mc_coef = 2.0 ==> van Leer MC
      MC_beta ( x, y, mc_coef) = AMIN1 ( 0.5 * ( ABS ( x + y)),
     1                         mc_coef * ABS ( x), mc_coef * ABS ( y) )
     1                    * ( SIGN ( 0.5, x) + SIGN ( 0.5, y))

! ----------------------------------------------------------------------

      smallnum = phys_consts ( 3)

      mc_coef = 1.5

      r13_3 = 13.0 / 3.0

      central_zone_wt = 100.0

      iprim = 1

! ------------------------------
!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iprim, ifail, idirection, ii,&
!$OMP&     ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,                    &
!$OMP&     ixxp2, iyyp2, izzp2, ixxm2, iyym2, izzm2,                    &

!$OMP&     Eigenval_pt, Lf_Eigenvec_pt, Rt_Eigenvec_pt, u_pt, v_pt,     &
!$OMP&     a, b, c, u_ip2, u_ip1, u_i, u_im1, u_im2,                    &

!$OMP&     R1_x, R2_x, R3_x, R1_xx, R2_xx, R3_xx,                       &

!$OMP&     Is1, Is2, Is3, Wt1, Wt2, Wt3)                                &

!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&     n_cc_components, central_zone_wt, smallnum, r13_3,           &
!$OMP&     n_phys_consts, phys_consts,                                  &
!$OMP&     u, flattengr, mc_coef, lo_order_zone_avg, n_cc_modes)        &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE(3)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

! --------------------

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)
          izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
          izzm2 = MAX0 ( izz - 2, 1 - ioffz)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 1 - ioffy)

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
          ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
          ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)

          c = 1.0

! ----------------------------------------------------------------------
! r=3 WENO limiter at third order starts here.
! ----------------------------------------------------------------------
! Starting Multidimentional characteristic reconstruction.

! ----------
! x moment.

          u_ip2( :) = u ( ixxp2, iyy, izz, :, 1)
          u_ip1( :) = u ( ixxp1, iyy, izz, :, 1)
          u_i  ( :) = u ( ixx  , iyy, izz, :, 1)
          u_im1( :) = u ( ixxm1, iyy, izz, :, 1)
          u_im2( :) = u ( ixxm2, iyy, izz, :, 1)

          R1_x = - 2.0 * u_im1 + 0.5 * u_im2 + 1.5 * u_i
          R1_xx = 0.5 * ( u_im2 - 2.0 * u_im1 + u_i )
          Is1 = R1_x**2 + r13_3 * R1_xx**2

          R2_x = 0.5 * ( u_ip1 - u_im1 )
          R2_xx = 0.5 * ( u_im1 - 2.0 * u_i + u_ip1 )
          Is2 = R2_x**2 + r13_3 * R2_xx**2

          R3_x = - 1.5 * u_i + 2.0 * u_ip1 - 0.5 * u_ip2
          R3_xx = 0.5 * ( u_i - 2.0 * u_ip1 + u_ip2 )
          Is3 = R3_x**2 + r13_3 * R3_xx**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          b ( :) = ( Wt1(:) * R1_x(:) + Wt2(:) * R2_x(:)
     1             + Wt3(:) * R3_x(:)) * c

          u ( ixx, iyy, izz, :, 2) = b ( :)

! ----------
! y moment.

          u_ip2 ( :) = u ( ixx, iyyp2, izz, :, 1)
          u_ip1 ( :) = u ( ixx, iyyp1, izz, :, 1)
          u_i   ( :) = u ( ixx, iyy  , izz, :, 1)
          u_im1 ( :) = u ( ixx, iyym1, izz, :, 1)
          u_im2 ( :) = u ( ixx, iyym2, izz, :, 1)

          R1_x = - 2.0 * u_im1 + 0.5 * u_im2 + 1.5 * u_i
          R1_xx = 0.5 * ( u_im2 - 2.0 * u_im1 + u_i )
          Is1 = R1_x**2 + r13_3 * R1_xx**2

          R2_x = 0.5 * ( u_ip1 - u_im1 )
          R2_xx = 0.5 * ( u_im1 - 2.0 * u_i + u_ip1 )
          Is2 = R2_x**2 + r13_3 * R2_xx**2

          R3_x = - 1.5 * u_i + 2.0 * u_ip1 - 0.5 * u_ip2
          R3_xx = 0.5 * ( u_i - 2.0 * u_ip1 + u_ip2 )
          Is3 = R3_x**2 + r13_3 * R3_xx**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          b ( :) = ( Wt1(:) * R1_x(:) + Wt2(:) * R2_x(:)
     1             + Wt3(:) * R3_x(:)) * c

          u ( ixx, iyy, izz, :, 3) = b ( :)

! --------------------
! z moment.

          u_ip2 ( :) = u ( ixx, iyy, izzp2, :, 1)
          u_ip1 ( :) = u ( ixx, iyy, izzp1, :, 1)
          u_i   ( :) = u ( ixx, iyy, izz  , :, 1)
          u_im1 ( :) = u ( ixx, iyy, izzm1, :, 1)
          u_im2 ( :) = u ( ixx, iyy, izzm2, :, 1)


          R1_x = - 2.0 * u_im1 + 0.5 * u_im2 + 1.5 * u_i
          R1_xx = 0.5 * ( u_im2 - 2.0 * u_im1 + u_i )
          Is1 = R1_x**2 + r13_3 * R1_xx**2

          R2_x = 0.5 * ( u_ip1 - u_im1 )
          R2_xx = 0.5 * ( u_im1 - 2.0 * u_i + u_ip1 )
          Is2 = R2_x**2 + r13_3 * R2_xx**2

          R3_x = - 1.5 * u_i + 2.0 * u_ip1 - 0.5 * u_ip2
          R3_xx = 0.5 * ( u_i - 2.0 * u_ip1 + u_ip2 )
          Is3 = R3_x**2 + r13_3 * R3_xx**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          b ( :) = ( Wt1(:) * R1_x(:) + Wt2(:) * R2_x(:)
     1             + Wt3(:) * R3_x(:)) * c

          u ( ixx, iyy, izz, :, 4) = b ( :)

! ----------------------------------------------------------------------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END SUBROUTINE LIMIT_RECONST_ALL_ZC_P1

#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if (NUMORDER == 3)

      SUBROUTINE LIMIT_RECONST_ALL_ZC_P2

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_modes, n_fc_modes, n_cc_efld_modes, n_cc_components,
     1        num_level, num_grid, num_pe,
     1        n_phys_consts, phys_consts,

     1        u, flattengr, lo_order_zone_avg,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Uses "efxbdy, efybdy, efzbdy to build "sourcegr". By construction, this
! "sourcegr" already has the factor of "4 * pi" built in.

! The higher moments of "sourcegr" are then built using third order
! WENO reconstruction.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_modes, n_fc_modes, n_cc_efld_modes, n_cc_components,
     1        num_level, num_grid, num_pe, n_phys_consts

      REAL phys_consts ( n_phys_consts)

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) :: lo_order_zone_avg

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::

     1      flattengr

! ----------

! Arrays that specify geometry.

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

      INTEGER ixx, iyy, izz, ixx0, iyy0, izz0,
     1        ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,
     1        ixxp2, iyyp2, izzp2, ixxm2, iyym2, izzm2,
     1        iprim, ifail, idirection, ii

      REAL c, r13_3, central_zone_wt, MC_beta, x, y, smallnum, mc_coef

      REAL, DIMENSION ( n_cc_components) :: a, b,
     1     u_ip2, u_ip1, u_i, u_im1, u_im2,

     1     R1_x, R2_x, R3_x, R1_xx, R2_xx, R3_xx,

     1     R1_xy, R2_xy, R3_xy, R4_xy,

     1     Is1, Is2, Is3, Is4,
     1     Wt1, Wt2, Wt3, Wt4,
     1     u_pt, v_pt, Eigenval_pt

      REAL, DIMENSION ( n_cc_components, n_cc_components) ::
     1      Lf_Eigenvec_pt, Rt_Eigenvec_pt

! MC_beta limiter : mc_coef = 1.0 ==> MinMod; mc_coef = 2.0 ==> van Leer MC
      MC_beta ( x, y, mc_coef) = AMIN1 ( 0.5 * ( ABS ( x + y)),
     1                         mc_coef * ABS ( x), mc_coef * ABS ( y) )
     1                    * ( SIGN ( 0.5, x) + SIGN ( 0.5, y))

! ----------------------------------------------------------------------

      smallnum = phys_consts ( 3)

      mc_coef = 1.5

      r13_3 = 13.0 / 3.0

      central_zone_wt = 100.0

      iprim = 1

! ----------------------------------------------------------------------
! ------------------------------
!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iprim, ifail, idirection, ii,&
!$OMP&     ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,                    &
!$OMP&     ixxp2, iyyp2, izzp2, ixxm2, iyym2, izzm2,                    &

!$OMP&     Eigenval_pt, Lf_Eigenvec_pt, Rt_Eigenvec_pt, u_pt, v_pt,     &
!$OMP&     a, b, c, u_ip2, u_ip1, u_i, u_im1, u_im2,                    &

!$OMP&     R1_x, R2_x, R3_x, R1_xx, R2_xx, R3_xx,                       &

!$OMP&     R1_xy, R2_xy, R3_xy, R4_xy,                                  &

!$OMP&     Is1, Is2, Is3, Is4, Wt1, Wt2, Wt3, Wt4)                      &

!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&     n_cc_components, central_zone_wt, smallnum, r13_3,           &
!$OMP&     n_phys_consts, phys_consts,                                  &
!$OMP&     u, flattengr, mc_coef, lo_order_zone_avg, n_cc_modes)        &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE(3)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)
          izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
          izzm2 = MAX0 ( izz - 2, 1 - ioffz)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 1 - ioffy)

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
          ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
          ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)

! --------------------

          c = 1.0

! ----------------------------------------------------------------------
! WENO LIMITING:
! -------------------------------
!  x, xx moments.

          u_ip2( :) = u ( ixxp2, iyy, izz, :, 1)
          u_ip1( :) = u ( ixxp1, iyy, izz, :, 1)
          u_i  ( :) = u ( ixx  , iyy, izz, :, 1)
          u_im1( :) = u ( ixxm1, iyy, izz, :, 1)
          u_im2( :) = u ( ixxm2, iyy, izz, :, 1)

          R1_x = - 2.0 * u_im1 + 0.5 * u_im2 + 1.5 * u_i
          R1_xx = 0.5 * ( u_im2 - 2.0 * u_im1 + u_i )
          Is1 = R1_x**2 + r13_3 * R1_xx**2

          R2_x = 0.5 * ( u_ip1 - u_im1 )
          R2_xx = 0.5 * ( u_im1 - 2.0 * u_i + u_ip1 )
          Is2 = R2_x**2 + r13_3 * R2_xx**2

          R3_x = - 1.5 * u_i + 2.0 * u_ip1 - 0.5 * u_ip2
          R3_xx = 0.5 * ( u_i - 2.0 * u_ip1 + u_ip2 )
          Is3 = R3_x**2 + r13_3 * R3_xx**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          b ( :) = ( Wt1(:) * R1_x(:) + Wt2(:) * R2_x(:)
     1             + Wt3(:) * R3_x(:)) * c

          u ( ixx, iyy, izz, :, 2) = b ( :)

          b ( :) = ( Wt1(:) * R1_xx(:) + Wt2(:) * R2_xx(:)
     1             + Wt3(:) * R3_xx(:)) * c

          u ( ixx, iyy, izz, :, 5) = b ( :)

! -------------------------------
! y, yy moments.

          u_ip2 ( :) = u ( ixx, iyyp2, izz, :, 1)
          u_ip1 ( :) = u ( ixx, iyyp1, izz, :, 1)
          u_i   ( :) = u ( ixx, iyy  , izz, :, 1)
          u_im1 ( :) = u ( ixx, iyym1, izz, :, 1)
          u_im2 ( :) = u ( ixx, iyym2, izz, :, 1)

          R1_x = - 2.0 * u_im1 + 0.5 * u_im2 + 1.5 * u_i
          R1_xx = 0.5 * ( u_im2 - 2.0 * u_im1 + u_i )
          Is1 = R1_x**2 + r13_3 * R1_xx**2

          R2_x = 0.5 * ( u_ip1 - u_im1 )
          R2_xx = 0.5 * ( u_im1 - 2.0 * u_i + u_ip1 )
          Is2 = R2_x**2 + r13_3 * R2_xx**2

          R3_x = - 1.5 * u_i + 2.0 * u_ip1 - 0.5 * u_ip2
          R3_xx = 0.5 * ( u_i - 2.0 * u_ip1 + u_ip2 )
          Is3 = R3_x**2 + r13_3 * R3_xx**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          b ( :) = ( Wt1(:) * R1_x(:) + Wt2(:) * R2_x(:)
     1             + Wt3(:) * R3_x(:)) * c

          u ( ixx, iyy, izz, :, 3) = b ( :)

! ------

          b ( :) = ( Wt1(:) * R1_xx(:) + Wt2(:) * R2_xx(:)
     1             + Wt3(:) * R3_xx(:)) * c

          u ( ixx, iyy, izz, :, 6) = b ( :)


! -------------------------------
! z, zz moments.

          u_ip2 ( :) = u ( ixx, iyy, izzp2, :, 1)
          u_ip1 ( :) = u ( ixx, iyy, izzp1, :, 1)
          u_i   ( :) = u ( ixx, iyy, izz  , :, 1)
          u_im1 ( :) = u ( ixx, iyy, izzm1, :, 1)
          u_im2 ( :) = u ( ixx, iyy, izzm2, :, 1)

          R1_x = - 2.0 * u_im1 + 0.5 * u_im2 + 1.5 * u_i
          R1_xx = 0.5 * ( u_im2 - 2.0 * u_im1 + u_i )
          Is1 = R1_x**2 + r13_3 * R1_xx**2

          R2_x = 0.5 * ( u_ip1 - u_im1 )
          R2_xx = 0.5 * ( u_im1 - 2.0 * u_i + u_ip1 )
          Is2 = R2_x**2 + r13_3 * R2_xx**2

          R3_x = - 1.5 * u_i + 2.0 * u_ip1 - 0.5 * u_ip2
          R3_xx = 0.5 * ( u_i - 2.0 * u_ip1 + u_ip2 )
          Is3 = R3_x**2 + r13_3 * R3_xx**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = 1.0 / ( Is3 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3

          b ( :) = ( Wt1(:) * R1_x(:) + Wt2(:) * R2_x(:)
     1             + Wt3(:) * R3_x(:)) * c

          u ( ixx, iyy, izz, :, 4) = b ( :)

! ------

          b ( :) = ( Wt1(:) * R1_xx(:) + Wt2(:) * R2_xx(:)
     1             + Wt3(:) * R3_xx(:)) * c

          u ( ixx, iyy, izz, :, 7) = b ( :)

! ----------
! xy moment.

          a (  :) = 4.0 * ( u ( ixx, iyy, izz,  :, 5)**2
     1              + u ( ixx, iyy, izz,  :, 6)**2 ) + smallnum

! -----

          R1_xy (  :) = u ( ixxp1, iyyp1, izz,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 2)
     1          - u ( ixx, iyy, izz,  :, 3)
     1          - u ( ixx, iyy, izz,  :, 5)
     1          - u ( ixx, iyy, izz,  :, 6)

          Is1 (  :) = R1_xy (  :)**2 + a (  :)

! -----

          R2_xy (  :) = - u ( ixxp1, iyym1, izz,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 2)
     1          - u ( ixx, iyy, izz,  :, 3)
     1          + u ( ixx, iyy, izz,  :, 5)
     1          + u ( ixx, iyy, izz,  :, 6)

          Is2 (  :) = R2_xy (  :)**2 + a (  :)

! -----

          R3_xy (  :) = - u ( ixxm1, iyyp1, izz,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 2)
     1          + u ( ixx, iyy, izz,  :, 3)
     1          + u ( ixx, iyy, izz,  :, 5)
     1          + u ( ixx, iyy, izz,  :, 6)

          Is3 (  :) = R3_xy (  :)**2 + a (  :)

! -----

          R4_xy (  :) = u ( ixxm1, iyym1, izz,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 2)
     1          + u ( ixx, iyy, izz,  :, 3)
     1          - u ( ixx, iyy, izz,  :, 5)
     1          - u ( ixx, iyy, izz,  :, 6)

          Is4 (  :) = R4_xy (  :)**2 + a (  :)

! -----
            Wt1 (  :) = 1.0 / Is1 (  :)**4
            Wt2 (  :) = 1.0 / Is2 (  :)**4
            Wt3 (  :) = 1.0 / Is3 (  :)**4
            Wt4 (  :) = 1.0 / Is4 (  :)**4

            a (  :) = 1.0 / ( Wt1 (  :) + Wt2 (  :)
     1                  + Wt3 (  :) + Wt4 (  :))

            Wt1 (  :) = Wt1 (  :) * a (  :)
            Wt2 (  :) = Wt2 (  :) * a (  :)
            Wt3 (  :) = Wt3 (  :) * a (  :)
            Wt4 (  :) = Wt4 (  :) * a (  :)

! -----

            u ( ixx, iyy, izz,  :, 8) =
     1       ( Wt1 (  :) * R1_xy (  :)
     1       + Wt2 (  :) * R2_xy (  :)
     1       + Wt3 (  :) * R3_xy (  :)
     1       + Wt4 (  :) * R4_xy (  :)) * c

! ----------
! xz moment.

          a (  :) = 4.0 * ( u ( ixx, iyy, izz,  :, 5)**2
     1              + u ( ixx, iyy, izz,  :, 7)**2 ) + smallnum

! -----

          R1_xy (  :) = u ( ixxp1, iyy, izzp1,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 2)
     1          - u ( ixx, iyy, izz,  :, 4)
     1          - u ( ixx, iyy, izz,  :, 5)
     1          - u ( ixx, iyy, izz,  :, 7)

          Is1 (  :) = R1_xy (  :)**2 + a (  :)

! -----

          R2_xy (  :) = - u ( ixxp1, iyy, izzm1,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 2)
     1          - u ( ixx, iyy, izz,  :, 4)
     1          + u ( ixx, iyy, izz,  :, 5)
     1          + u ( ixx, iyy, izz,  :, 7)

          Is2 (  :) = R2_xy (  :)**2 + a (  :)

! -----

          R3_xy (  :) = - u ( ixxm1, iyy, izzp1,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 2)
     1          + u ( ixx, iyy, izz,  :, 4)
     1          + u ( ixx, iyy, izz,  :, 5)
     1          + u ( ixx, iyy, izz,  :, 7)

          Is3 (  :) = R3_xy (  :)**2 + a (  :)

! -----

          R4_xy (  :) = u ( ixxm1, iyy, izzm1,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 2)
     1          + u ( ixx, iyy, izz,  :, 4)
     1          - u ( ixx, iyy, izz,  :, 5)
     1          - u ( ixx, iyy, izz,  :, 7)

          Is4 (  :) = R4_xy (  :)**2 + a (  :)

! -----

            Wt1 (  :) = 1.0 / Is1 (  :)**4
            Wt2 (  :) = 1.0 / Is2 (  :)**4
            Wt3 (  :) = 1.0 / Is3 (  :)**4
            Wt4 (  :) = 1.0 / Is4 (  :)**4

            a (  :) = 1.0 / ( Wt1 (  :) + Wt2 (  :)
     1                  + Wt3 (  :) + Wt4 (  :))

            Wt1 (  :) = Wt1 (  :) * a (  :)
            Wt2 (  :) = Wt2 (  :) * a (  :)
            Wt3 (  :) = Wt3 (  :) * a (  :)
            Wt4 (  :) = Wt4 (  :) * a (  :)

! -----

            u ( ixx, iyy, izz,  :, 10) =
     1       ( Wt1 (  :) * R1_xy (  :)
     1       + Wt2 (  :) * R2_xy (  :)
     1       + Wt3 (  :) * R3_xy (  :)
     1       + Wt4 (  :) * R4_xy (  :)) * c

! ----------
! yz moment.

          a (  :) = 4.0 * ( u ( ixx, iyy, izz,  :, 6)**2
     1              + u ( ixx, iyy, izz,  :, 7)**2 ) + smallnum

! -----


          R1_xy (  :) = u ( ixx, iyyp1, izzp1,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 3)
     1          - u ( ixx, iyy, izz,  :, 4)
     1          - u ( ixx, iyy, izz,  :, 6)
     1          - u ( ixx, iyy, izz,  :, 7)

          Is1 (  :) = R1_xy (  :)**2 + a (  :)

! -----

          R2_xy (  :) = - u ( ixx, iyyp1, izzm1,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 3)
     1          - u ( ixx, iyy, izz,  :, 4)
     1          + u ( ixx, iyy, izz,  :, 6)
     1          + u ( ixx, iyy, izz,  :, 7)

          Is2 (  :) = R2_xy (  :)**2 + a (  :)

! -----

          R3_xy (  :) = - u ( ixx, iyym1, izzp1,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 3)
     1          + u ( ixx, iyy, izz,  :, 4)
     1          + u ( ixx, iyy, izz,  :, 6)
     1          + u ( ixx, iyy, izz,  :, 7)

          Is3 (  :) = R3_xy (  :)**2 + a (  :)

! -----

          R4_xy (  :) = u ( ixx, iyym1, izzm1,  :, 1)
     1          - u ( ixx, iyy, izz,  :, 1)
     1          + u ( ixx, iyy, izz,  :, 3)
     1          + u ( ixx, iyy, izz,  :, 4)
     1          - u ( ixx, iyy, izz,  :, 6)
     1          - u ( ixx, iyy, izz,  :, 7)

          Is4 (  :) = R4_xy (  :)**2 + a (  :)

! -----

            Wt1 (  :) = 1.0 / Is1 (  :)**4
            Wt2 (  :) = 1.0 / Is2 (  :)**4
            Wt3 (  :) = 1.0 / Is3 (  :)**4
            Wt4 (  :) = 1.0 / Is4 (  :)**4

            a (  :) = 1.0 / ( Wt1 (  :) + Wt2 (  :)
     1                  + Wt3 (  :) + Wt4 (  :))

            Wt1 (  :) = Wt1 (  :) * a (  :)
            Wt2 (  :) = Wt2 (  :) * a (  :)
            Wt3 (  :) = Wt3 (  :) * a (  :)
            Wt4 (  :) = Wt4 (  :) * a (  :)

! -----

            u ( ixx, iyy, izz,  :, 9) =
     1       ( Wt1 (  :) * R1_xy (  :)
     1       + Wt2 (  :) * R2_xy (  :)
     1       + Wt3 (  :) * R3_xy (  :)
     1       + Wt4 (  :) * R4_xy (  :)) * c

! ----------------------------------------------------------------------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END SUBROUTINE LIMIT_RECONST_ALL_ZC_P2

#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if ( NUMORDER == 4)

      SUBROUTINE LIMIT_RECONST_ALL_ZC_P3

     1      ( ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_modes, n_fc_modes, n_cc_efld_modes, n_cc_components,
     1        num_level, num_grid, num_pe,
     1        n_phys_consts, phys_consts,

     1        u, flattengr, lo_order_zone_avg,

     1        xb_save, dx_save, xc_save, yb_save, dy_save, yc_save,
     1        zb_save, dz_save, zc_save)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Uses "efxbdy, efybdy, efzbdy to build "sourcegr". By construction, this
! "sourcegr" already has the factor of "4 * pi" built in.

! The higher moments of "sourcegr" are then built using fourth order
! WENO reconstruction.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ix1, iy1, iz1, ioffx, ioffy, ioffz,
     1        n_cc_modes, n_fc_modes, n_cc_efld_modes, n_cc_components,
     1        num_level, num_grid, num_pe, n_phys_consts

      REAL phys_consts ( n_phys_consts)

      INTEGER, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) :: lo_order_zone_avg

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz,
     1                  n_cc_components, n_cc_modes) :: u

      REAL, DIMENSION ( 1 - ioffx : ix1 + ioffx,
     1                  1 - ioffy : iy1 + ioffy,
     1                  1 - ioffz : iz1 + ioffz) ::

     1      flattengr

! ----------

! Arrays that specify geometry.

      REAL xb_save ( 0 - ioffx : ix1 + ioffx)
      REAL dx_save ( 1 - ioffx : ix1 + ioffx)
      REAL xc_save ( 1 - ioffx : ix1 + ioffx)

      REAL yb_save ( 0 - ioffy : iy1 + ioffy)
      REAL dy_save ( 1 - ioffy : iy1 + ioffy)
      REAL yc_save ( 1 - ioffy : iy1 + ioffy)

      REAL zb_save ( 0 - ioffz : iz1 + ioffz)
      REAL dz_save ( 1 - ioffz : iz1 + ioffz)
      REAL zc_save ( 1 - ioffz : iz1 + ioffz)

! ----------------------------------------------------------------------

      INTEGER ixx, iyy, izz, ixx0, iyy0, izz0,
     1        ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,
     1        ixxp2, iyyp2, izzp2, ixxm2, iyym2, izzm2,
     1        ixxp3, iyyp3, izzp3, ixxm3, iyym3, izzm3,
     1        iprim, ifail, idirection, ii

      REAL c, r13_3, r781_20, central_zone_wt, MC_beta, x, y,
     1     smallnum, mc_coef

      REAL, DIMENSION ( n_cc_components) :: 

     1     a, b, u_ip3, u_ip2, u_ip1, u_i, u_im1, u_im2, u_im3,

     1     R1_x, R2_x, R3_x, R4_x, R1_xx, R2_xx, R3_xx, R4_xx,
     1     R1_xxx, R2_xxx, R3_xxx, R4_xxx,

     1     R1_xy, R2_xy, R3_xy, R4_xy, R5_xy,
     1     R1_xxy, R2_xxy, R3_xxy, R4_xxy, R5_xxy,
     1     R1_xyy, R2_xyy, R3_xyy, R4_xyy, R5_xyy,

     1     Is1, Is2, Is3, Is4, Is5, Is6, Is7, Is8,
     1     Wt1, Wt2, Wt3, Wt4, Wt5, Wt6, Wt7, Wt8,
     1  R1_xyz, R2_xyz, R3_xyz, R4_xyz, R5_xyz, R6_xyz, R7_xyz, R8_xyz,
     1     u_pt, v_pt, Eigenval_pt

      REAL, DIMENSION ( n_cc_components, n_cc_components) ::
     1      Lf_Eigenvec_pt, Rt_Eigenvec_pt

! MC_beta limiter : mc_coef = 1.0 ==> MinMod; mc_coef = 2.0 ==> van Leer MC
      MC_beta ( x, y, mc_coef) = AMIN1 ( 0.5 * ( ABS ( x + y)),
     1                         mc_coef * ABS ( x), mc_coef * ABS ( y) )
     1                    * ( SIGN ( 0.5, x) + SIGN ( 0.5, y))

! ----------------------------------------------------------------------

      smallnum = phys_consts ( 3)

      mc_coef = 1.5

      r13_3 = 13.0 / 3.0

      r781_20 = 781.0 / 20.0

      central_zone_wt = 100.0

      iprim = 1

! ----------------------------------------------------------------------
! Now make the modes (with flattening) for the zone-centered "sourcegr".

! ------------------------------
!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iprim, ifail, idirection, ii,&
!$OMP&     ixxp1, iyyp1, izzp1, ixxm1, iyym1, izzm1,                    &
!$OMP&     ixxp2, iyyp2, izzp2, ixxm2, iyym2, izzm2,                    &
!$OMP&     ixxp3, iyyp3, izzp3, ixxm3, iyym3, izzm3,                    &

!$OMP&     Eigenval_pt, Lf_Eigenvec_pt, Rt_Eigenvec_pt, u_pt, v_pt,     &
!$OMP&     a, b, c, u_ip3, u_ip2, u_ip1, u_i, u_im1, u_im2, u_im3,      &

!$OMP&     R1_x, R2_x, R3_x, R4_x, R1_xx, R2_xx, R3_xx, R4_xx,          &
!$OMP&     R1_xxx, R2_xxx, R3_xxx, R4_xxx,                              &

!$OMP&     R1_xy, R2_xy, R3_xy, R4_xy, R5_xy,                           &
!$OMP&     R1_xxy, R2_xxy, R3_xxy, R4_xxy, R5_xxy,                      &
!$OMP&     R1_xyy, R2_xyy, R3_xyy, R4_xyy, R5_xyy,                      &

!$OMP&   R1_xyz, R2_xyz, R3_xyz, R4_xyz, R5_xyz, R6_xyz, R7_xyz, R8_xyz,&

!$OMP&     Is1, Is2, Is3, Is4, Is5, Is6, Is7, Is8,                      &
!$OMP&     Wt1, Wt2, Wt3, Wt4, Wt5, Wt6, Wt7, Wt8)                      &

!$OMP&            SHARED ( ix1, ioffx, iy1, ioffy, iz1, ioffz,          &
!$OMP&     n_cc_components, central_zone_wt, smallnum, r13_3, r781_20,  &
!$OMP&     n_phys_consts, phys_consts,                                  &
!$OMP&     u, flattengr, mc_coef, lo_order_zone_avg, n_cc_modes)        &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)
#if ( TRUE_2D == 1)
!$OMP&            COLLAPSE(3)
#endif

      DO izz = 1 - ioffz, iz1 + ioffz

        DO iyy = 1 - ioffy, iy1 + ioffy

          DO ixx = 1 - ioffx, ix1 + ioffx

          izzp1 = MIN0 ( izz + 1, iz1 + ioffz)
          izzm1 = MAX0 ( izz - 1, 1 - ioffz)
          izzp2 = MIN0 ( izz + 2, iz1 + ioffz)
          izzm2 = MAX0 ( izz - 2, 1 - ioffz)
          izzp3 = MIN0 ( izz + 3, iz1 + ioffz)
          izzm3 = MAX0 ( izz - 3, 1 - ioffz)

          iyyp1 = MIN0 ( iyy + 1, iy1 + ioffy)
          iyym1 = MAX0 ( iyy - 1, 1 - ioffy)
          iyyp2 = MIN0 ( iyy + 2, iy1 + ioffy)
          iyym2 = MAX0 ( iyy - 2, 1 - ioffy)
          iyyp3 = MIN0 ( iyy + 3, iy1 + ioffy)
          iyym3 = MAX0 ( iyy - 3, 1 - ioffy)

          ixxp1 = MIN0 ( ixx + 1, ix1 + ioffx)
          ixxm1 = MAX0 ( ixx - 1, 1 - ioffx)
          ixxp2 = MIN0 ( ixx + 2, ix1 + ioffx)
          ixxm2 = MAX0 ( ixx - 2, 1 - ioffx)
          ixxp3 = MIN0 ( ixx + 3, ix1 + ioffx)
          ixxm3 = MAX0 ( ixx - 3, 1 - ioffx)

! ----------

          c = 1.0 - flattengr ( ixx, iyy, izz)

! ------------------------------
! Preparing for the characteristic reconstruction.

#if ( CHARACTERISTIC_RECONSTRUCTION == 1)

          u_pt ( :) = u ( ixx, iyy, izz, :, 1)

! Reconstruct in primitives

          v_pt = u_pt

#endif
! End of "#if ( CHARACTERISTIC_RECONSTRUCTION == 1)"

! ----------------------------------------------------------------------
! TVD Limiting if lo_order_zone_avg ( ixx, iyy, izz) == 1
! ----------------------------------------------------------------------

      IF ( lo_order_zone_avg ( ixx, iyy, izz) .EQ. 1) THEN

! ------------------------------

          u ( ixx, iyy, izz, :, 2:n_cc_modes) = 0.0

! ------------------------------
! Starting Multidimentional characteristic reconstruction.

! ------------------------------
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)

          idirection = 1

          CALL PRIM_TO_EIGENSYS
     1           ( idirection, iprim, n_phys_consts, phys_consts,
     1             v_pt,

     1             Eigenval_pt,
     1             Lf_Eigenvec_pt, Rt_Eigenvec_pt)

#endif
! End of "#if ( CHARACTERISTIC_RECONSTRUCTION == 1)"
! ───────────────────────────────────────────────────────
! x-moment

          u_ip1 ( :) = u ( ixxp1, iyy, izz, :, 1)
          u_i   ( :) = u ( ixx  , iyy, izz, :, 1)
          u_im1 ( :) = u ( ixxm1, iyy, izz, :, 1)

#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          u_ip1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip1 ( :) )
          u_i   ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_i   ( :) )
          u_im1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im1 ( :) )
#endif

          DO ii = 1, n_cc_components
            b ( ii) = MC_beta ( u_ip1 ( ii) - u_i   ( ii),
     1                          u_i   ( ii) - u_im1 ( ii), mc_coef)
          END DO

#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif

          u ( ixx, iyy, izz, :, 2) = b ( :)


! ───────────────────────────────────────────────────────
! y moment

! -------------------------------
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
! -------------------------------

          idirection = 2

          CALL PRIM_TO_EIGENSYS
     1           ( idirection, iprim, n_phys_consts, phys_consts,
     1             v_pt,

     1             Eigenval_pt,
     1             Lf_Eigenvec_pt, Rt_Eigenvec_pt)

#endif
! End of "#if ( CHARACTERISTIC_RECONSTRUCTION == 1)"
! ----------

          u_ip1 ( :) = u ( ixx, iyyp1, izz, :, 1)
          u_i   ( :) = u ( ixx, iyy  , izz, :, 1)
          u_im1 ( :) = u ( ixx, iyym1, izz, :, 1)

#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          u_ip1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip1 ( :) )
          u_i   ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_i   ( :) )
          u_im1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im1 ( :) )
#endif

          DO ii = 1, n_cc_components
            b ( ii) = MC_beta ( u_ip1 ( ii) - u_i   ( ii),
     1                          u_i   ( ii) - u_im1 ( ii), mc_coef)
          END DO

#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif

          u ( ixx, iyy, izz, :, 3) = b ( :)


! ───────────────────────────────────────────────────────
! z moment

! -------------------------------
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
! -------------------------------

          idirection = 3

          CALL PRIM_TO_EIGENSYS
     1           ( idirection, iprim, n_phys_consts, phys_consts,
     1             v_pt,

     1             Eigenval_pt,
     1             Lf_Eigenvec_pt, Rt_Eigenvec_pt)

#endif
! End of "#if ( CHARACTERISTIC_RECONSTRUCTION == 1)"
! ----------

          u_ip1 ( :) = u ( ixx, iyy, izzp1, :, 1)
          u_i   ( :) = u ( ixx, iyy, izz  , :, 1)
          u_im1 ( :) = u ( ixx, iyy, izzm1, :, 1)

#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          u_ip1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip1 ( :) )
          u_i   ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_i   ( :) )
          u_im1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im1 ( :) )
#endif

          DO ii = 1, n_cc_components
            b ( ii) = MC_beta ( u_ip1 ( ii) - u_i   ( ii),
     1                          u_i   ( ii) - u_im1 ( ii), mc_coef)
          END DO

#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif

          u ( ixx, iyy, izz, :, 4) = b ( :)

! ───────────────────────────────────────────────────────

      ELSE

! ----------------------------------------------------------------------
! r=4 WENO limiter at fourth order starts here.
! ----------------------------------------------------------------------

! WENO LIMITING:
! ----------

! -------------------------------
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
! -------------------------------

          idirection = 1

          CALL PRIM_TO_EIGENSYS
     1           ( idirection, iprim, n_phys_consts, phys_consts,
     1             v_pt,

     1             Eigenval_pt,
     1             Lf_Eigenvec_pt, Rt_Eigenvec_pt)

#endif
! End of "#if ( CHARACTERISTIC_RECONSTRUCTION == 1)"
! ----------
! x, xx, xxx moments.

          u_ip3 ( :) = u ( ixxp3, iyy, izz, :, 1)
          u_ip2 ( :) = u ( ixxp2, iyy, izz, :, 1)
          u_ip1 ( :) = u ( ixxp1, iyy, izz, :, 1)
          u_i   ( :) = u ( ixx,   iyy, izz, :, 1)
          u_im1 ( :) = u ( ixxm1, iyy, izz, :, 1)
          u_im2 ( :) = u ( ixxm2, iyy, izz, :, 1)
          u_im3 ( :) = u ( ixxm3, iyy, izz, :, 1)

#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          u_ip3 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip3 ( :))
          u_ip2 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip2 ( :))
          u_ip1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip1 ( :))
          u_i   ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_i   ( :))
          u_im1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im1 ( :))
          u_im2 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im2 ( :))
          u_im3 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im3 ( :))
#endif

          R1_x = ( - 177.0 * u_im1 + 87.0 * u_im2 - 19.0 * u_im3
     1           + 109.0 * u_i) / 60.0
          R1_xx = - 2.5 * u_im1 + 2.0 * u_im2 - 0.5 * u_im3 + u_i
          R1_xxx = ( - 3.0 * u_im1 + 3.0 * u_im2 - u_im3 + u_i) / 6.0
          Is1 = ( R1_x + 0.1 * R1_xxx)**2 + r13_3 * R1_xx**2
     1        + r781_20 * R1_xxx**2

          R2_x = ( - 63.0 * u_im1 + 11.0 * u_im2 + 33.0 * u_i
     1           + 19.0 * u_ip1) / 60.0
          R2_xx = 0.5 * u_im1 - u_i + 0.5 * u_ip1
          R2_xxx = ( 3.0 * u_im1 - u_im2 - 3.0 * u_i + u_ip1) / 6.0
          Is2 = ( R2_x + 0.1 * R2_xxx)**2 + r13_3 * R2_xx**2
     1        + r781_20 * R2_xxx**2

          R3_x = ( - 19.0 * u_im1 - 33.0 * u_i + 63.0 * u_ip1
     1           - 11.0 * u_ip2) / 60.0
          R3_xx = 0.5 * u_im1 - u_i + 0.5 * u_ip1
          R3_xxx = ( - u_im1 + 3.0 * u_i - 3.0 * u_ip1 + u_ip2) / 6.0
          Is3 = ( R3_x + 0.1 * R3_xxx)**2 + r13_3 * R3_xx**2
     1        + r781_20 * R3_xxx**2

          R4_x = ( - 109.0 * u_i + 177.0 * u_ip1 - 87.0 * u_ip2
     1           + 19.0 * u_ip3) / 60.0
          R4_xx = u_i - 2.5 * u_ip1 + 2.0 * u_ip2 - 0.5 * u_ip3
          R4_xxx = ( - u_i + 3.0 * u_ip1 - 3.0 * u_ip2 + u_ip3) / 6.0
          Is4 = ( R4_x + 0.1 * R4_xxx)**2 + r13_3 * R4_xx**2
     1        + r781_20 * R4_xxx**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = central_zone_wt / ( Is3 + smallnum)**4
          Wt4 = 1.0 / ( Is4 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3
          Wt4 = a * Wt4

! ------------

          b ( :) = ( Wt1(:) * R1_x(:) + Wt2(:) * R2_x(:)
     1             + Wt3(:) * R3_x(:) + Wt4(:) * R4_x(:) ) * c
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif
          u ( ixx, iyy, izz, :, 2) = b ( :)

! ------

          b ( :) = ( Wt1(:) * R1_xx(:) + Wt2(:) * R2_xx(:)
     1             + Wt3(:) * R3_xx(:) + Wt4(:) * R4_xx(:) ) * c
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif
          u ( ixx, iyy, izz, :, 5) = b ( :)

! ------

          b ( :) = ( Wt1(:) * R1_xxx(:) + Wt2(:) * R2_xxx(:)
     1             + Wt3(:) * R3_xxx(:) + Wt4(:) * R4_xxx(:) ) * c
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif
          u ( ixx, iyy, izz, :, 11) = b ( :)

! ----------------------------------------------------------------------

! -------------------------------
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
! -------------------------------

          idirection = 2

          CALL PRIM_TO_EIGENSYS
     1           ( idirection, iprim, n_phys_consts, phys_consts,
     1             v_pt,

     1             Eigenval_pt,
     1             Lf_Eigenvec_pt, Rt_Eigenvec_pt)

#endif
! End of "#if ( CHARACTERISTIC_RECONSTRUCTION == 1)"
! --------------------
! y, yy, yyy moments.


          u_ip3 ( :) = u ( ixx, iyyp3, izz, :, 1)
          u_ip2 ( :) = u ( ixx, iyyp2, izz, :, 1)
          u_ip1 ( :) = u ( ixx, iyyp1, izz, :, 1)
          u_i   ( :) = u ( ixx, iyy  , izz, :, 1)
          u_im1 ( :) = u ( ixx, iyym1, izz, :, 1)
          u_im2 ( :) = u ( ixx, iyym2, izz, :, 1)
          u_im3 ( :) = u ( ixx, iyym3, izz, :, 1)

#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          u_ip3 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip3 ( :))
          u_ip2 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip2 ( :))
          u_ip1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip1 ( :))
          u_i   ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_i   ( :))
          u_im1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im1 ( :))
          u_im2 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im2 ( :))
          u_im3 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im3 ( :))
#endif


          R1_x = ( - 177.0 * u_im1 + 87.0 * u_im2 - 19.0 * u_im3
     1           + 109.0 * u_i) / 60.0
          R1_xx = - 2.5 * u_im1 + 2.0 * u_im2 - 0.5 * u_im3 + u_i
          R1_xxx = ( - 3.0 * u_im1 + 3.0 * u_im2 - u_im3 + u_i) / 6.0
          Is1 = ( R1_x + 0.1 * R1_xxx)**2 + r13_3 * R1_xx**2
     1        + r781_20 * R1_xxx**2

          R2_x = ( - 63.0 * u_im1 + 11.0 * u_im2 + 33.0 * u_i
     1           + 19.0 * u_ip1) / 60.0
          R2_xx = 0.5 * u_im1 - u_i + 0.5 * u_ip1
          R2_xxx = ( 3.0 * u_im1 - u_im2 - 3.0 * u_i + u_ip1) / 6.0
          Is2 = ( R2_x + 0.1 * R2_xxx)**2 + r13_3 * R2_xx**2
     1        + r781_20 * R2_xxx**2

          R3_x = ( - 19.0 * u_im1 - 33.0 * u_i + 63.0 * u_ip1
     1           - 11.0 * u_ip2) / 60.0
          R3_xx = 0.5 * u_im1 - u_i + 0.5 * u_ip1
          R3_xxx = ( - u_im1 + 3.0 * u_i - 3.0 * u_ip1 + u_ip2) / 6.0
          Is3 = ( R3_x + 0.1 * R3_xxx)**2 + r13_3 * R3_xx**2
     1        + r781_20 * R3_xxx**2

          R4_x = ( - 109.0 * u_i + 177.0 * u_ip1 - 87.0 * u_ip2
     1           + 19.0 * u_ip3) / 60.0
          R4_xx = u_i - 2.5 * u_ip1 + 2.0 * u_ip2 - 0.5 * u_ip3
          R4_xxx = ( - u_i + 3.0 * u_ip1 - 3.0 * u_ip2 + u_ip3) / 6.0
          Is4 = ( R4_x + 0.1 * R4_xxx)**2 + r13_3 * R4_xx**2
     1        + r781_20 * R4_xxx**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = central_zone_wt / ( Is3 + smallnum)**4
          Wt4 = 1.0 / ( Is4 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3
          Wt4 = a * Wt4

! ------------

          b ( :) = ( Wt1(:) * R1_x(:) + Wt2(:) * R2_x(:)
     1             + Wt3(:) * R3_x(:) + Wt4(:) * R4_x(:) ) * c
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif
          u ( ixx, iyy, izz, :, 3) = b ( :)

! ------

          b ( :) = ( Wt1(:) * R1_xx(:) + Wt2(:) * R2_xx(:)
     1             + Wt3(:) * R3_xx(:) + Wt4(:) * R4_xx(:) ) * c
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif
          u ( ixx, iyy, izz, :, 6) = b ( :)

! ------

          b ( :) = ( Wt1(:) * R1_xxx(:) + Wt2(:) * R2_xxx(:)
     1             + Wt3(:) * R3_xxx(:) + Wt4(:) * R4_xxx(:) ) * c
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif
          u ( ixx, iyy, izz, :, 12) = b ( :)

! ----------------------------------------------------------------------

! -------------------------------
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
! -------------------------------

          idirection = 3

          CALL PRIM_TO_EIGENSYS
     1           ( idirection, iprim, n_phys_consts, phys_consts,
     1             v_pt,

     1             Eigenval_pt,
     1             Lf_Eigenvec_pt, Rt_Eigenvec_pt)

#endif
! End of "#if ( CHARACTERISTIC_RECONSTRUCTION == 1)"
! --------------------
! z, zz, zzz moments.


          u_ip3 ( :) = u ( ixx, iyy, izzp3, :, 1)
          u_ip2 ( :) = u ( ixx, iyy, izzp2, :, 1)
          u_ip1 ( :) = u ( ixx, iyy, izzp1, :, 1)
          u_i   ( :) = u ( ixx, iyy, izz  , :, 1)
          u_im1 ( :) = u ( ixx, iyy, izzm1, :, 1)
          u_im2 ( :) = u ( ixx, iyy, izzm2, :, 1)
          u_im3 ( :) = u ( ixx, iyy, izzm3, :, 1)

#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          u_ip3 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip3 ( :))
          u_ip2 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip2 ( :))
          u_ip1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_ip1 ( :))
          u_i   ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_i   ( :))
          u_im1 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im1 ( :))
          u_im2 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im2 ( :))
          u_im3 ( :) = MATMUL ( Lf_Eigenvec_pt(:,:), u_im3 ( :))
#endif

          R1_x = ( - 177.0 * u_im1 + 87.0 * u_im2 - 19.0 * u_im3
     1           + 109.0 * u_i) / 60.0
          R1_xx = - 2.5 * u_im1 + 2.0 * u_im2 - 0.5 * u_im3 + u_i
          R1_xxx = ( - 3.0 * u_im1 + 3.0 * u_im2 - u_im3 + u_i) / 6.0
          Is1 = ( R1_x + 0.1 * R1_xxx)**2 + r13_3 * R1_xx**2
     1        + r781_20 * R1_xxx**2

          R2_x = ( - 63.0 * u_im1 + 11.0 * u_im2 + 33.0 * u_i
     1           + 19.0 * u_ip1) / 60.0
          R2_xx = 0.5 * u_im1 - u_i + 0.5 * u_ip1
          R2_xxx = ( 3.0 * u_im1 - u_im2 - 3.0 * u_i + u_ip1) / 6.0
          Is2 = ( R2_x + 0.1 * R2_xxx)**2 + r13_3 * R2_xx**2
     1        + r781_20 * R2_xxx**2

          R3_x = ( - 19.0 * u_im1 - 33.0 * u_i + 63.0 * u_ip1
     1           - 11.0 * u_ip2) / 60.0
          R3_xx = 0.5 * u_im1 - u_i + 0.5 * u_ip1
          R3_xxx = ( - u_im1 + 3.0 * u_i - 3.0 * u_ip1 + u_ip2) / 6.0
          Is3 = ( R3_x + 0.1 * R3_xxx)**2 + r13_3 * R3_xx**2
     1        + r781_20 * R3_xxx**2

          R4_x = ( - 109.0 * u_i + 177.0 * u_ip1 - 87.0 * u_ip2
     1           + 19.0 * u_ip3) / 60.0
          R4_xx = u_i - 2.5 * u_ip1 + 2.0 * u_ip2 - 0.5 * u_ip3
          R4_xxx = ( - u_i + 3.0 * u_ip1 - 3.0 * u_ip2 + u_ip3) / 6.0
          Is4 = ( R4_x + 0.1 * R4_xxx)**2 + r13_3 * R4_xx**2
     1        + r781_20 * R4_xxx**2

          Wt1 = 1.0 / ( Is1 + smallnum)**4
          Wt2 = central_zone_wt / ( Is2 + smallnum)**4
          Wt3 = central_zone_wt / ( Is3 + smallnum)**4
          Wt4 = 1.0 / ( Is4 + smallnum)**4

          a = 1.0 / ( Wt1 + Wt2 + Wt3 + Wt4)
          Wt1 = a * Wt1
          Wt2 = a * Wt2
          Wt3 = a * Wt3
          Wt4 = a * Wt4

! ------------

          b ( :) = ( Wt1(:) * R1_x(:) + Wt2(:) * R2_x(:)
     1             + Wt3(:) * R3_x(:) + Wt4(:) * R4_x(:) ) * c
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif
          u ( ixx, iyy, izz, :, 4) = b ( :)

! ------

          b ( :) = ( Wt1(:) * R1_xx(:) + Wt2(:) * R2_xx(:)
     1             + Wt3(:) * R3_xx(:) + Wt4(:) * R4_xx(:) ) * c
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif
          u ( ixx, iyy, izz, :, 7) = b ( :)

! ------

          b ( :) = ( Wt1(:) * R1_xxx(:) + Wt2(:) * R2_xxx(:)
     1             + Wt3(:) * R3_xxx(:) + Wt4(:) * R4_xxx(:) ) * c
#if ( CHARACTERISTIC_RECONSTRUCTION == 1)
          b ( :) = MATMUL ( Rt_Eigenvec_pt( :, :), b ( :))
#endif
          u ( ixx, iyy, izz, :, 13) = b ( :)

! ----------------------------------------------------------------------
! xy, xxy, xyy moments.

          a ( :) = 39.0 * ( u ( ixx, iyy, izz, :, 11)**2
     1                 + u ( ixx, iyy, izz, :, 12)**2 )
     1        + 4.0 * ( u ( ixx, iyy, izz, :, 5)**2
     1                + u ( ixx, iyy, izz, :, 6)**2 )

! -----

          R1_xy ( :) = ( 60.0 * u ( ixxp1, iyyp1, izz, :, 1)
     1            - 10.0 * u ( ixxp1, iyyp2, izz, :, 1)
     1            - 10.0 * u ( ixxp2, iyyp1, izz, :, 1)
     1            - 40.0 * u ( ixx, iyy, izz, :, 1)
     1            - 30.0 * u ( ixx, iyy, izz, :, 2)
     1            - 30.0 * u ( ixx, iyy, izz, :, 3)
     1            - 10.0 * u ( ixx, iyy, izz, :, 5)
     1            - 10.0 * u ( ixx, iyy, izz, :, 6)
     1            + 27.0 * u ( ixx, iyy, izz, :, 11)
     1            + 27.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0

          R1_xxy ( :) = ( - 20.0 * u ( ixxp1, iyyp1, izz, :, 1)
     1               + 10.0 * u ( ixxp2, iyyp1, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 3)
     1               - 20.0 * u ( ixx, iyy, izz, :, 5)
     1               + 10.0 * u ( ixx, iyy, izz, :, 6)
     1               - 60.0 * u ( ixx, iyy, izz, :, 11)
     1               + 11.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0

          R1_xyy ( :) = ( - 20.0 * u ( ixxp1, iyyp1, izz, :, 1)
     1               + 10.0 * u ( ixxp1, iyyp2, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 2)
     1               + 10.0 * u ( ixx, iyy, izz, :, 5)
     1               - 20.0 * u ( ixx, iyy, izz, :, 6)
     1               + 11.0 * u ( ixx, iyy, izz, :, 11)
     1               - 60.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0

          Is1 ( :) = a ( :) + R1_xy ( :)**2 + 14.0 * ( R1_xxy ( :)**2
     1               + R1_xyy ( :)**2 ) / 3.0

! -----

          R2_xy ( :) =  ( - 60.0 * u ( ixxm1, iyyp1, izz, :, 1)
     1               + 10.0 * u ( ixxm1, iyyp2, izz, :, 1)
     1               + 10.0 * u ( ixxm2, iyyp1, izz, :, 1)
     1               + 40.0 * u ( ixx, iyy, izz, :, 1)
     1               - 30.0 * u ( ixx, iyy, izz, :, 2)
     1               + 30.0 * u ( ixx, iyy, izz, :, 3)
     1               + 10.0 * u ( ixx, iyy, izz, :, 5)
     1               + 10.0 * u ( ixx, iyy, izz, :, 6)
     1               + 27.0 * u ( ixx, iyy, izz, :, 11)
     1               - 27.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0

          R2_xxy ( :) = ( - 20.0 * u ( ixxm1, iyyp1, izz, :, 1)
     1               + 10.0 * u ( ixxm2, iyyp1, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 3)
     1               - 20.0 * u ( ixx, iyy, izz, :, 5)
     1               + 10.0 * u ( ixx, iyy, izz, :, 6)
     1               + 60.0 * u ( ixx, iyy, izz, :, 11)
     1               + 11.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0

          R2_xyy ( :) = ( 20.0 * u ( ixxm1, iyyp1, izz, :, 1)
     1             - 10.0 * u ( ixxm1, iyyp2, izz, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 2)
     1             - 10.0 * u ( ixx, iyy, izz, :, 5)
     1             + 20.0 * u ( ixx, iyy, izz, :, 6)
     1             + 11.0 * u ( ixx, iyy, izz, :, 11)
     1             + 60.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0

          Is2 ( :) = a ( :) + R2_xy ( :)**2 + 14.0 * ( R2_xxy ( :)**2
     1             + R2_xyy ( :)**2 ) / 3.0

! -----

          R3_xy ( :) = ( - 60.0 * u ( ixxp1, iyym1, izz, :, 1)
     1              + 10.0 * u ( ixxp1, iyym2, izz, :, 1)
     1              + 10.0 * u ( ixxp2, iyym1, izz, :, 1)
     1              + 40.0 * u ( ixx, iyy, izz, :, 1)
     1              + 30.0 * u ( ixx, iyy, izz, :, 2)
     1              - 30.0 * u ( ixx, iyy, izz, :, 3)
     1              + 10.0 * u ( ixx, iyy, izz, :, 5)
     1              + 10.0 * u ( ixx, iyy, izz, :, 6)
     1              - 27.0 * u ( ixx, iyy, izz, :, 11)
     1              + 27.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0


          R3_xxy ( :) = ( 20.0 * u ( ixxp1, iyym1, izz, :, 1)
     1             - 10.0 * u ( ixxp2, iyym1, izz, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 3)
     1             + 20.0 * u ( ixx, iyy, izz, :, 5)
     1             - 10.0 * u ( ixx, iyy, izz, :, 6)
     1             + 60.0 * u ( ixx, iyy, izz, :, 11)
     1             + 11.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0


          R3_xyy ( :) = ( - 20.0 * u ( ixxp1, iyym1, izz, :, 1)
     1               + 10.0 * u ( ixxp1, iyym2, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 2)
     1               + 10.0 * u ( ixx, iyy, izz, :, 5)
     1               - 20.0 * u ( ixx, iyy, izz, :, 6)
     1               + 11.0 * u ( ixx, iyy, izz, :, 11)
     1               + 60.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0


          Is3 ( :) = a ( :) + R3_xy ( :)**2 + 14.0 * ( R3_xxy ( :)**2
     1               + R3_xyy ( :)**2 ) / 3.0

! -----

          R4_xy ( :) = ( 60.0 * u ( ixxm1, iyym1, izz, :, 1)
     1            - 10.0 * u ( ixxm2, iyym1, izz, :, 1)
     1            - 10.0 * u ( ixxm1, iyym2, izz, :, 1)
     1            - 40.0 * u ( ixx, iyy, izz, :, 1)
     1            + 30.0 * u ( ixx, iyy, izz, :, 2)
     1            + 30.0 * u ( ixx, iyy, izz, :, 3)
     1            - 10.0 * u ( ixx, iyy, izz, :, 5)
     1            - 10.0 * u ( ixx, iyy, izz, :, 6)
     1            - 27.0 * u ( ixx, iyy, izz, :, 11)
     1            - 27.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0

          R4_xxy ( :) = ( 20.0 * u ( ixxm1, iyym1, izz, :, 1)
     1             - 10.0 * u ( ixxm2, iyym1, izz, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 3)
     1             + 20.0 * u ( ixx, iyy, izz, :, 5)
     1             - 10.0 * u ( ixx, iyy, izz, :, 6)
     1             - 60.0 * u ( ixx, iyy, izz, :, 11)
     1             + 11.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0


          R4_xyy ( :) = ( 20.0 * u ( ixxm1, iyym1, izz, :, 1)
     1             - 10.0 * u ( ixxm1, iyym2, izz, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 2)
     1             - 10.0 * u ( ixx, iyy, izz, :, 5)
     1             + 20.0 * u ( ixx, iyy, izz, :, 6)
     1             + 11.0 * u ( ixx, iyy, izz, :, 11)
     1             - 60.0 * u ( ixx, iyy, izz, :, 12) ) / 20.0

          Is4 ( :) = a ( :) + R4_xy ( :)**2 + 14.0 * ( R4_xxy ( :)**2
     1             + R4_xyy ( :)**2 ) / 3.0

! -----

          R5_xy ( :) = ( u ( ixxm1, iyym1, izz, :, 1)
     1            - u ( ixxm1, iyyp1, izz, :, 1)
     1            - u ( ixxp1, iyym1, izz, :, 1)
     1            + u ( ixxp1, iyyp1, izz, :, 1) ) / 4.0

          R5_xxy ( :) = ( - 5.0 * u ( ixxm1, iyym1, izz, :, 1)
     1               + 5.0 * u ( ixxm1, iyyp1, izz, :, 1)
     1               - 5.0 * u ( ixxp1, iyym1, izz, :, 1)
     1               + 5.0 * u ( ixxp1, iyyp1, izz, :, 1)
     1               - 22.0 * u ( ixx, iyy, izz, :, 12)
     1               - 20.0 * u ( ixx, iyy, izz, :, 3) ) / 20.0

          R5_xyy ( :) = ( - 5.0 * u ( ixxm1, iyym1, izz, :, 1)
     1               - 5.0 * u ( ixxm1, iyyp1, izz, :, 1)
     1               + 5.0 * u ( ixxp1, iyym1, izz, :, 1)
     1               + 5.0 * u ( ixxp1, iyyp1, izz, :, 1)
     1               - 22.0 * u ( ixx, iyy, izz, :, 11)
     1               - 20.0 * u ( ixx, iyy, izz, :, 2) ) / 20.0

          Is5 ( :) = a ( :) + R5_xy ( :)**2 + 14.0 * ( R5_xxy ( :)**2
     1             + R5_xyy ( :)**2 ) / 3.0

! -----

          Wt1 ( :) = 1.0 / ( Is1 ( :) + smallnum)**4
          Wt2 ( :) = 1.0 / ( Is2 ( :) + smallnum)**4
          Wt3 ( :) = 1.0 / ( Is3 ( :) + smallnum)**4
          Wt4 ( :) = 1.0 / ( Is4 ( :) + smallnum)**4
          Wt5 ( :) = central_zone_wt / ( Is5 ( :) + smallnum)**4

          b ( :) = 1.0 / ( Wt1(:) + Wt2(:) + Wt3(:) + Wt4(:) + Wt5(:))
          Wt1 ( :)= b ( :) * Wt1 ( :)
          Wt2 ( :)= b ( :) * Wt2 ( :)
          Wt3 ( :)= b ( :) * Wt3 ( :)
          Wt4 ( :)= b ( :) * Wt4 ( :)
          Wt5 ( :)= b ( :) * Wt5 ( :)

! -----

            u ( ixx, iyy, izz, :, 8) = ( Wt1( :) * R1_xy( :)
     1      + Wt2( :) * R2_xy( :) + Wt3( :) * R3_xy( :)
     1      + Wt4( :) * R4_xy( :) + Wt5( :) * R5_xy( :) ) * c

            u ( ixx, iyy, izz, :, 14)=( Wt1( :) * R1_xxy( :)
     1      + Wt2( :) * R2_xxy( :) + Wt3( :) * R3_xxy( :)
     1      + Wt4( :) * R4_xxy( :) + Wt5( :) * R5_xxy( :) ) * c

            u ( ixx, iyy, izz, :, 16)=( Wt1( :) * R1_xyy( :)
     1      + Wt2( :) * R2_xyy( :) + Wt3( :) * R3_xyy( :)
     1      + Wt4( :) * R4_xyy( :) + Wt5( :) * R5_xyy( :) ) * c

! ----------
! xz, xxz, xzz moments.

          a ( :) = 39.0 * ( u ( ixx, iyy, izz, :, 11)**2
     1                 + u ( ixx, iyy, izz, :, 13)**2 )
     1        + 4.0 * ( u ( ixx, iyy, izz, :, 5)**2
     1                + u ( ixx, iyy, izz, :, 7)**2 )

! -----

          R1_xy ( :) = ( 60.0 * u ( ixxp1, iyy, izzp1, :, 1)
     1            - 10.0 * u ( ixxp1, iyy, izzp2, :, 1)
     1            - 10.0 * u ( ixxp2, iyy, izzp1, :, 1)
     1            - 40.0 * u ( ixx, iyy, izz, :, 1)
     1            - 30.0 * u ( ixx, iyy, izz, :, 2)
     1            - 30.0 * u ( ixx, iyy, izz, :, 4)
     1            - 10.0 * u ( ixx, iyy, izz, :, 5)
     1            - 10.0 * u ( ixx, iyy, izz, :, 7)
     1            + 27.0 * u ( ixx, iyy, izz, :, 11)
     1            + 27.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          R1_xxy ( :) = ( - 20.0 * u ( ixxp1, iyy, izzp1, :, 1)
     1               + 10.0 * u ( ixxp2, iyy, izzp1, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 4)
     1               - 20.0 * u ( ixx, iyy, izz, :, 5)
     1               + 10.0 * u ( ixx, iyy, izz, :, 7)
     1               - 60.0 * u ( ixx, iyy, izz, :, 11)
     1               + 11.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          R1_xyy ( :) = ( - 20.0 * u ( ixxp1, iyy, izzp1, :, 1)
     1               + 10.0 * u ( ixxp1, iyy, izzp2, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 2)
     1               + 10.0 * u ( ixx, iyy, izz, :, 5)
     1               - 20.0 * u ( ixx, iyy, izz, :, 7)
     1               + 11.0 * u ( ixx, iyy, izz, :, 11)
     1               - 60.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          Is1 ( :) = a ( :) + R1_xy ( :)**2 + 14.0 * ( R1_xxy ( :)**2
     1             + R1_xyy ( :)**2 ) / 3.0

! -----

          R2_xy ( :) =  ( - 60.0 * u ( ixxm1, iyy, izzp1, :, 1)
     1               + 10.0 * u ( ixxm1, iyy, izzp2, :, 1)
     1               + 10.0 * u ( ixxm2, iyy, izzp1, :, 1)
     1               + 40.0 * u ( ixx, iyy, izz, :, 1)
     1               - 30.0 * u ( ixx, iyy, izz, :, 2)
     1               + 30.0 * u ( ixx, iyy, izz, :, 4)
     1               + 10.0 * u ( ixx, iyy, izz, :, 5)
     1               + 10.0 * u ( ixx, iyy, izz, :, 7)
     1               + 27.0 * u ( ixx, iyy, izz, :, 11)
     1               - 27.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          R2_xxy ( :) = ( - 20.0 * u ( ixxm1, iyy, izzp1, :, 1)
     1               + 10.0 * u ( ixxm2, iyy, izzp1, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 4)
     1               - 20.0 * u ( ixx, iyy, izz, :, 5)
     1               + 10.0 * u ( ixx, iyy, izz, :, 7)
     1               + 60.0 * u ( ixx, iyy, izz, :, 11)
     1               + 11.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          R2_xyy ( :) = ( 20.0 * u ( ixxm1, iyy, izzp1, :, 1)
     1             - 10.0 * u ( ixxm1, iyy, izzp2, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 2)
     1             - 10.0 * u ( ixx, iyy, izz, :, 5)
     1             + 20.0 * u ( ixx, iyy, izz, :, 7)
     1             + 11.0 * u ( ixx, iyy, izz, :, 11)
     1             + 60.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          Is2 ( :) = a ( :) + R2_xy ( :)**2 + 14.0 * ( R2_xxy ( :)**2
     1             + R2_xyy ( :)**2 ) / 3.0

! -----

          R3_xy ( :) = ( - 60.0 * u ( ixxp1, iyy, izzm1, :, 1)
     1              + 10.0 * u ( ixxp1, iyy, izzm2, :, 1)
     1              + 10.0 * u ( ixxp2, iyy, izzm1, :, 1)
     1              + 40.0 * u ( ixx, iyy, izz, :, 1)
     1              + 30.0 * u ( ixx, iyy, izz, :, 2)
     1              - 30.0 * u ( ixx, iyy, izz, :, 4)
     1              + 10.0 * u ( ixx, iyy, izz, :, 5)
     1              + 10.0 * u ( ixx, iyy, izz, :, 7)
     1              - 27.0 * u ( ixx, iyy, izz, :, 11)
     1              + 27.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0


          R3_xxy ( :) = ( 20.0 * u ( ixxp1, iyy, izzm1, :, 1)
     1             - 10.0 * u ( ixxp2, iyy, izzm1, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 4)
     1             + 20.0 * u ( ixx, iyy, izz, :, 5)
     1             - 10.0 * u ( ixx, iyy, izz, :, 7)
     1             + 60.0 * u ( ixx, iyy, izz, :, 11)
     1             + 11.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0


          R3_xyy ( :) = ( - 20.0 * u ( ixxp1, iyy, izzm1, :, 1)
     1               + 10.0 * u ( ixxp1, iyy, izzm2, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 2)
     1               + 10.0 * u ( ixx, iyy, izz, :, 5)
     1               - 20.0 * u ( ixx, iyy, izz, :, 7)
     1               + 11.0 * u ( ixx, iyy, izz, :, 11)
     1               + 60.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0


          Is3 ( :) = a ( :) + R3_xy ( :)**2 + 14.0 * ( R3_xxy ( :)**2
     1             + R3_xyy ( :)**2 ) / 3.0

! -----

          R4_xy ( :) = ( 60.0 * u ( ixxm1, iyy, izzm1, :, 1)
     1            - 10.0 * u ( ixxm2, iyy, izzm1, :, 1)
     1            - 10.0 * u ( ixxm1, iyy, izzm2, :, 1)
     1            - 40.0 * u ( ixx, iyy, izz, :, 1)
     1            + 30.0 * u ( ixx, iyy, izz, :, 2)
     1            + 30.0 * u ( ixx, iyy, izz, :, 4)
     1            - 10.0 * u ( ixx, iyy, izz, :, 5)
     1            - 10.0 * u ( ixx, iyy, izz, :, 7)
     1            - 27.0 * u ( ixx, iyy, izz, :, 11)
     1            - 27.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          R4_xxy ( :) = ( 20.0 * u ( ixxm1, iyy, izzm1, :, 1)
     1             - 10.0 * u ( ixxm2, iyy, izzm1, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 4)
     1             + 20.0 * u ( ixx, iyy, izz, :, 5)
     1             - 10.0 * u ( ixx, iyy, izz, :, 7)
     1             - 60.0 * u ( ixx, iyy, izz, :, 11)
     1             + 11.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0


          R4_xyy ( :) = ( 20.0 * u ( ixxm1, iyy, izzm1, :, 1)
     1             - 10.0 * u ( ixxm1, iyy, izzm2, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 2)
     1             - 10.0 * u ( ixx, iyy, izz, :, 5)
     1             + 20.0 * u ( ixx, iyy, izz, :, 7)
     1             + 11.0 * u ( ixx, iyy, izz, :, 11)
     1             - 60.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          Is4 ( :) = a ( :) + R4_xy ( :)**2 + 14.0 * ( R4_xxy ( :)**2
     1             + R4_xyy ( :)**2 ) / 3.0

! -----

          R5_xy ( :) = ( u ( ixxm1, iyy, izzm1, :, 1)
     1            - u ( ixxm1, iyy, izzp1, :, 1)
     1            - u ( ixxp1, iyy, izzm1, :, 1)
     1            + u ( ixxp1, iyy, izzp1, :, 1) ) / 4.0

          R5_xxy ( :) = ( - 5.0 * u ( ixxm1, iyy, izzm1, :, 1)
     1               + 5.0 * u ( ixxm1, iyy, izzp1, :, 1)
     1               - 5.0 * u ( ixxp1, iyy, izzm1, :, 1)
     1               + 5.0 * u ( ixxp1, iyy, izzp1, :, 1)
     1               - 22.0 * u ( ixx, iyy, izz, :, 13)
     1               - 20.0 * u ( ixx, iyy, izz, :, 4) ) / 20.0

          R5_xyy ( :) = ( - 5.0 * u ( ixxm1, iyy, izzm1, :, 1)
     1               - 5.0 * u ( ixxm1, iyy, izzp1, :, 1)
     1               + 5.0 * u ( ixxp1, iyy, izzm1, :, 1)
     1               + 5.0 * u ( ixxp1, iyy, izzp1, :, 1)
     1               - 22.0 * u ( ixx, iyy, izz, :, 11)
     1               - 20.0 * u ( ixx, iyy, izz, :, 2) ) / 20.0

          Is5 ( :) = a ( :) + R5_xy ( :)**2 + 14.0 * ( R5_xxy ( :)**2
     1             + R5_xyy ( :)**2 ) / 3.0

! -----

          Wt1 ( :) = 1.0 / ( Is1 ( :) + smallnum)**4
          Wt2 ( :) = 1.0 / ( Is2 ( :) + smallnum)**4
          Wt3 ( :) = 1.0 / ( Is3 ( :) + smallnum)**4
          Wt4 ( :) = 1.0 / ( Is4 ( :) + smallnum)**4
          Wt5 ( :) = central_zone_wt / ( Is5 ( :) + smallnum)**4

          b ( :) = 1.0 / ( Wt1(:) + Wt2(:) + Wt3(:) + Wt4(:) + Wt5(:))
          Wt1 ( :)= b ( :) * Wt1 ( :)
          Wt2 ( :)= b ( :) * Wt2 ( :)
          Wt3 ( :)= b ( :) * Wt3 ( :)
          Wt4 ( :)= b ( :) * Wt4 ( :)
          Wt5 ( :)= b ( :) * Wt5 ( :)

! -----

            u ( ixx, iyy, izz, :, 10) = ( Wt1(:) * R1_xy(:)
     1      + Wt2(:) * R2_xy(:) + Wt3(:) * R3_xy(:) + Wt4(:) * R4_xy(:)
     1      + Wt5(:) * R5_xy(:) ) * c

            u ( ixx, iyy, izz, :, 15)=( Wt1(:) * R1_xxy(:)
     1      + Wt2(:) * R2_xxy(:) + Wt3(:) * R3_xxy(:)
     1      + Wt4(:) * R4_xxy(:) + Wt5(:) * R5_xxy(:) ) * c

            u ( ixx, iyy, izz, :, 18)=( Wt1(:) * R1_xyy(:)
     1      + Wt2(:) * R2_xyy(:) + Wt3(:) * R3_xyy(:)
     1      + Wt4(:) * R4_xyy(:) + Wt5(:) * R5_xyy(:) ) * c

! ----------
! yz, yyz, yzz moments.

          a(:) = 39.0 * ( u ( ixx, iyy, izz, :, 12)**2
     1                 + u ( ixx, iyy, izz, :, 13)**2 )
     1        + 4.0 * ( u ( ixx, iyy, izz, :, 6)**2
     1                + u ( ixx, iyy, izz, :, 7)**2 )

! -----

          R1_xy(:) = ( 60.0 * u ( ixx, iyyp1, izzp1, :, 1)
     1            - 10.0 * u ( ixx, iyyp1, izzp2, :, 1)
     1            - 10.0 * u ( ixx, iyyp2, izzp1, :, 1)
     1            - 40.0 * u ( ixx, iyy, izz, :, 1)
     1            - 30.0 * u ( ixx, iyy, izz, :, 3)
     1            - 30.0 * u ( ixx, iyy, izz, :, 4)
     1            - 10.0 * u ( ixx, iyy, izz, :, 6)
     1            - 10.0 * u ( ixx, iyy, izz, :, 7)
     1            + 27.0 * u ( ixx, iyy, izz, :, 12)
     1            + 27.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          R1_xxy(:) = ( - 20.0 * u ( ixx, iyyp1, izzp1, :, 1)
     1               + 10.0 * u ( ixx, iyyp2, izzp1, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 4)
     1               - 20.0 * u ( ixx, iyy, izz, :, 6)
     1               + 10.0 * u ( ixx, iyy, izz, :, 7)
     1               - 60.0 * u ( ixx, iyy, izz, :, 12)
     1               + 11.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          R1_xyy(:) = ( - 20.0 * u ( ixx, iyyp1, izzp1, :, 1)
     1               + 10.0 * u ( ixx, iyyp1, izzp2, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 3)
     1               + 10.0 * u ( ixx, iyy, izz, :, 6)
     1               - 20.0 * u ( ixx, iyy, izz, :, 7)
     1               + 11.0 * u ( ixx, iyy, izz, :, 12)
     1               - 60.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          Is1(:) = a(:) + R1_xy(:)**2 + 14.0 * ( R1_xxy(:)**2
     1           + R1_xyy(:)**2 ) / 3.0

! -----

          R2_xy(:) =  ( - 60.0 * u ( ixx, iyym1, izzp1, :, 1)
     1               + 10.0 * u ( ixx, iyym1, izzp2, :, 1)
     1               + 10.0 * u ( ixx, iyym2, izzp1, :, 1)
     1               + 40.0 * u ( ixx, iyy, izz, :, 1)
     1               - 30.0 * u ( ixx, iyy, izz, :, 3)
     1               + 30.0 * u ( ixx, iyy, izz, :, 4)
     1               + 10.0 * u ( ixx, iyy, izz, :, 6)
     1               + 10.0 * u ( ixx, iyy, izz, :, 7)
     1               + 27.0 * u ( ixx, iyy, izz, :, 12)
     1               - 27.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          R2_xxy(:) = ( - 20.0 * u ( ixx, iyym1, izzp1, :, 1)
     1               + 10.0 * u ( ixx, iyym2, izzp1, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 4)
     1               - 20.0 * u ( ixx, iyy, izz, :, 6)
     1               + 10.0 * u ( ixx, iyy, izz, :, 7)
     1               + 60.0 * u ( ixx, iyy, izz, :, 12)
     1               + 11.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          R2_xyy(:) = ( 20.0 * u ( ixx, iyym1, izzp1, :, 1)
     1             - 10.0 * u ( ixx, iyym1, izzp2, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 3)
     1             - 10.0 * u ( ixx, iyy, izz, :, 6)
     1             + 20.0 * u ( ixx, iyy, izz, :, 7)
     1             + 11.0 * u ( ixx, iyy, izz, :, 12)
     1             + 60.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          Is2(:) = a (:)+ R2_xy(:)**2 + 14.0 * ( R2_xxy(:)**2
     1           + R2_xyy(:)**2 ) / 3.0

! -----

          R3_xy(:) = ( - 60.0 * u ( ixx, iyyp1, izzm1, :, 1)
     1              + 10.0 * u ( ixx, iyyp1, izzm2, :, 1)
     1              + 10.0 * u ( ixx, iyyp2, izzm1, :, 1)
     1              + 40.0 * u ( ixx, iyy, izz, :, 1)
     1              + 30.0 * u ( ixx, iyy, izz, :, 3)
     1              - 30.0 * u ( ixx, iyy, izz, :, 4)
     1              + 10.0 * u ( ixx, iyy, izz, :, 6)
     1              + 10.0 * u ( ixx, iyy, izz, :, 7)
     1              - 27.0 * u ( ixx, iyy, izz, :, 12)
     1              + 27.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0


          R3_xxy(:) = ( 20.0 * u ( ixx, iyyp1, izzm1, :, 1)
     1             - 10.0 * u ( ixx, iyyp2, izzm1, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 4)
     1             + 20.0 * u ( ixx, iyy, izz, :, 6)
     1             - 10.0 * u ( ixx, iyy, izz, :, 7)
     1             + 60.0 * u ( ixx, iyy, izz, :, 12)
     1             + 11.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0


          R3_xyy(:) = ( - 20.0 * u ( ixx, iyyp1, izzm1, :, 1)
     1               + 10.0 * u ( ixx, iyyp1, izzm2, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 1)
     1               + 10.0 * u ( ixx, iyy, izz, :, 3)
     1               + 10.0 * u ( ixx, iyy, izz, :, 6)
     1               - 20.0 * u ( ixx, iyy, izz, :, 7)
     1               + 11.0 * u ( ixx, iyy, izz, :, 12)
     1               + 60.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0


          Is3(:) = a(:) + R3_xy(:)**2 + 14.0 * ( R3_xxy(:)**2
     1           + R3_xyy(:)**2 ) / 3.0

! -----

          R4_xy(:) = ( 60.0 * u ( ixx, iyym1, izzm1, :, 1)
     1            - 10.0 * u ( ixx, iyym2, izzm1, :, 1)
     1            - 10.0 * u ( ixx, iyym1, izzm2, :, 1)
     1            - 40.0 * u ( ixx, iyy, izz, :, 1)
     1            + 30.0 * u ( ixx, iyy, izz, :, 3)
     1            + 30.0 * u ( ixx, iyy, izz, :, 4)
     1            - 10.0 * u ( ixx, iyy, izz, :, 6)
     1            - 10.0 * u ( ixx, iyy, izz, :, 7)
     1            - 27.0 * u ( ixx, iyy, izz, :, 12)
     1            - 27.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          R4_xxy(:) = ( 20.0 * u ( ixx, iyym1, izzm1, :, 1)
     1             - 10.0 * u ( ixx, iyym2, izzm1, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 4)
     1             + 20.0 * u ( ixx, iyy, izz, :, 6)
     1             - 10.0 * u ( ixx, iyy, izz, :, 7)
     1             - 60.0 * u ( ixx, iyy, izz, :, 12)
     1             + 11.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0


          R4_xyy(:) = ( 20.0 * u ( ixx, iyym1, izzm1, :, 1)
     1             - 10.0 * u ( ixx, iyym1, izzm2, :, 1)
     1             - 10.0 * u ( ixx, iyy, izz, :, 1)
     1             + 10.0 * u ( ixx, iyy, izz, :, 3)
     1             - 10.0 * u ( ixx, iyy, izz, :, 6)
     1             + 20.0 * u ( ixx, iyy, izz, :, 7)
     1             + 11.0 * u ( ixx, iyy, izz, :, 12)
     1             - 60.0 * u ( ixx, iyy, izz, :, 13) ) / 20.0

          Is4(:) = a(:) + R4_xy(:)**2 + 14.0 * ( R4_xxy(:)**2
     1           + R4_xyy(:)**2 ) / 3.0

! -----

          R5_xy(:) = ( u ( ixx, iyym1, izzm1, :, 1)
     1            - u ( ixx, iyym1, izzp1, :, 1)
     1            - u ( ixx, iyyp1, izzm1, :, 1)
     1            + u ( ixx, iyyp1, izzp1, :, 1) ) / 4.0

          R5_xxy(:) = ( - 5.0 * u ( ixx, iyym1, izzm1, :, 1)
     1               + 5.0 * u ( ixx, iyym1, izzp1, :, 1)
     1               - 5.0 * u ( ixx, iyyp1, izzm1, :, 1)
     1               + 5.0 * u ( ixx, iyyp1, izzp1, :, 1)
     1               - 22.0 * u ( ixx, iyy, izz, :, 13)
     1               - 20.0 * u ( ixx, iyy, izz, :, 4) ) / 20.0

          R5_xyy(:) = ( - 5.0 * u ( ixx, iyym1, izzm1, :, 1)
     1               - 5.0 * u ( ixx, iyym1, izzp1, :, 1)
     1               + 5.0 * u ( ixx, iyyp1, izzm1, :, 1)
     1               + 5.0 * u ( ixx, iyyp1, izzp1, :, 1)
     1               - 22.0 * u ( ixx, iyy, izz, :, 12)
     1               - 20.0 * u ( ixx, iyy, izz, :, 3) ) / 20.0

          Is5(:) = a(:) + R5_xy(:)**2 + 14.0 * ( R5_xxy(:)**2
     1           + R5_xyy(:)**2 ) / 3.0

! -----

          Wt1 ( :) = 1.0 / ( Is1 ( :) + smallnum)**4
          Wt2 ( :) = 1.0 / ( Is2 ( :) + smallnum)**4
          Wt3 ( :) = 1.0 / ( Is3 ( :) + smallnum)**4
          Wt4 ( :) = 1.0 / ( Is4 ( :) + smallnum)**4
          Wt5 ( :) = central_zone_wt / ( Is5 ( :) + smallnum)**4

          b ( :) = 1.0 / ( Wt1(:) + Wt2(:) + Wt3(:) + Wt4(:) + Wt5(:))
          Wt1 ( :)= b ( :) * Wt1 ( :)
          Wt2 ( :)= b ( :) * Wt2 ( :)
          Wt3 ( :)= b ( :) * Wt3 ( :)
          Wt4 ( :)= b ( :) * Wt4 ( :)
          Wt5 ( :)= b ( :) * Wt5 ( :)

! -----

            u ( ixx, iyy, izz, :, 9) = ( Wt1(:) * R1_xy(:)
     1      + Wt2(:) * R2_xy(:) + Wt3(:) * R3_xy(:)
     1      + Wt4(:) * R4_xy(:) + Wt5(:) * R5_xy(:) ) * c

            u ( ixx, iyy, izz, :, 17)=( Wt1(:) * R1_xxy(:)
     1      + Wt2(:) * R2_xxy(:) + Wt3(:) * R3_xxy(:)
     1      + Wt4(:) * R4_xxy(:) + Wt5(:) * R5_xxy(:) ) * c

            u ( ixx, iyy, izz, :, 19)=( Wt1(:) * R1_xyy(:)
     1      + Wt2(:) * R2_xyy(:) + Wt3(:) * R3_xyy(:)
     1      + Wt4(:) * R4_xyy(:) + Wt5(:) * R5_xyy(:) ) * c

! ----------
! xyz moment.

            a(:) = 36.0 *( u ( ixx, iyy, izz, :, 11)**2
     1                + u ( ixx, iyy, izz, :, 12)**2
     1                + u ( ixx, iyy, izz, :, 13)**2)
     1        + 4.0 * ( u ( ixx, iyy, izz, :, 14)**2
     1                + u ( ixx, iyy, izz, :, 15)**2
     1                + u ( ixx, iyy, izz, :, 16)**2
     1                + u ( ixx, iyy, izz, :, 17)**2
     1                + u ( ixx, iyy, izz, :, 18)**2
     1                + u ( ixx, iyy, izz, :, 19)**2)
     1        + smallnum

! -----

            R1_xyz(:) = u ( ixxp1, iyyp1, izzp1, :, 1)
     1         - 1.1 * ( u ( ixx, iyy, izz, :, 11)
     1                 + u ( ixx, iyy, izz, :, 12)
     1                 + u ( ixx, iyy, izz, :, 13))
     1         - u ( ixx, iyy, izz, :, 1)
     1         - u ( ixx, iyy, izz, :, 2)
     1         - u ( ixx, iyy, izz, :, 3)
     1         - u ( ixx, iyy, izz, :, 4)
     1         - u ( ixx, iyy, izz, :, 5)
     1         - u ( ixx, iyy, izz, :, 6)
     1         - u ( ixx, iyy, izz, :, 7)
     1         - u ( ixx, iyy, izz, :, 8)
     1         - u ( ixx, iyy, izz, :, 9)
     1         - u ( ixx, iyy, izz, :, 10)
     1         - u ( ixx, iyy, izz, :, 14)
     1         - u ( ixx, iyy, izz, :, 15)
     1         - u ( ixx, iyy, izz, :, 16)
     1         - u ( ixx, iyy, izz, :, 17)
     1         - u ( ixx, iyy, izz, :, 18)
     1         - u ( ixx, iyy, izz, :, 19)

           Is1(:) = R1_xyz(:)**2 + a(:)

! -----

            R2_xyz(:) = - u ( ixxm1, iyyp1, izzp1, :, 1)
     1         + 1.1 * (-u ( ixx, iyy, izz, :, 11)
     1                 + u ( ixx, iyy, izz, :, 12)
     1                 + u ( ixx, iyy, izz, :, 13))
     1         + u ( ixx, iyy, izz, :, 1)
     1         - u ( ixx, iyy, izz, :, 2)
     1         + u ( ixx, iyy, izz, :, 3)
     1         + u ( ixx, iyy, izz, :, 4)
     1         + u ( ixx, iyy, izz, :, 5)
     1         + u ( ixx, iyy, izz, :, 6)
     1         + u ( ixx, iyy, izz, :, 7)
     1         - u ( ixx, iyy, izz, :, 8)
     1         + u ( ixx, iyy, izz, :, 9)
     1         - u ( ixx, iyy, izz, :, 10)
     1         + u ( ixx, iyy, izz, :, 14)
     1         + u ( ixx, iyy, izz, :, 15)
     1         - u ( ixx, iyy, izz, :, 16)
     1         + u ( ixx, iyy, izz, :, 17)
     1         - u ( ixx, iyy, izz, :, 18)
     1         + u ( ixx, iyy, izz, :, 19)

           Is2(:) = R2_xyz(:)**2 + a(:)

! -----

            R3_xyz(:) = - u ( ixxp1, iyym1, izzp1, :, 1)
     1         + 1.1 * ( u ( ixx, iyy, izz, :, 11)
     1                 - u ( ixx, iyy, izz, :, 12)
     1                 + u ( ixx, iyy, izz, :, 13))
     1         + u ( ixx, iyy, izz, :, 1)
     1         + u ( ixx, iyy, izz, :, 2)
     1         - u ( ixx, iyy, izz, :, 3)
     1         + u ( ixx, iyy, izz, :, 4)
     1         + u ( ixx, iyy, izz, :, 5)
     1         + u ( ixx, iyy, izz, :, 6)
     1         + u ( ixx, iyy, izz, :, 7)
     1         - u ( ixx, iyy, izz, :, 8)
     1         - u ( ixx, iyy, izz, :, 9)
     1         + u ( ixx, iyy, izz, :, 10)
     1         - u ( ixx, iyy, izz, :, 14)
     1         + u ( ixx, iyy, izz, :, 15)
     1         + u ( ixx, iyy, izz, :, 16)
     1         + u ( ixx, iyy, izz, :, 17)
     1         + u ( ixx, iyy, izz, :, 18)
     1         - u ( ixx, iyy, izz, :, 19)

           Is3(:) = R3_xyz(:)**2 + a(:)

! -----

            R4_xyz(:) = - u ( ixxp1, iyyp1, izzm1, :, 1)
     1         + 1.1 * ( u ( ixx, iyy, izz, :, 11)
     1                 + u ( ixx, iyy, izz, :, 12)
     1                 - u ( ixx, iyy, izz, :, 13))
     1         + u ( ixx, iyy, izz, :, 1)
     1         + u ( ixx, iyy, izz, :, 2)
     1         + u ( ixx, iyy, izz, :, 3)
     1         - u ( ixx, iyy, izz, :, 4)
     1         + u ( ixx, iyy, izz, :, 5)
     1         + u ( ixx, iyy, izz, :, 6)
     1         + u ( ixx, iyy, izz, :, 7)
     1         + u ( ixx, iyy, izz, :, 8)
     1         - u ( ixx, iyy, izz, :, 9)
     1         - u ( ixx, iyy, izz, :, 10)
     1         + u ( ixx, iyy, izz, :, 14)
     1         - u ( ixx, iyy, izz, :, 15)
     1         + u ( ixx, iyy, izz, :, 16)
     1         - u ( ixx, iyy, izz, :, 17)
     1         + u ( ixx, iyy, izz, :, 18)
     1         + u ( ixx, iyy, izz, :, 19)

           Is4(:) = R4_xyz(:)**2 + a(:)

! -----

            R5_xyz(:) = u ( ixxm1, iyym1, izzp1, :, 1)
     1         + 1.1 * ( u ( ixx, iyy, izz, :, 11)
     1                 + u ( ixx, iyy, izz, :, 12)
     1                 - u ( ixx, iyy, izz, :, 13))
     1         - u ( ixx, iyy, izz, :, 1)
     1         + u ( ixx, iyy, izz, :, 2)
     1         + u ( ixx, iyy, izz, :, 3)
     1         - u ( ixx, iyy, izz, :, 4)
     1         - u ( ixx, iyy, izz, :, 5)
     1         - u ( ixx, iyy, izz, :, 6)
     1         - u ( ixx, iyy, izz, :, 7)
     1         - u ( ixx, iyy, izz, :, 8)
     1         + u ( ixx, iyy, izz, :, 9)
     1         + u ( ixx, iyy, izz, :, 10)
     1         + u ( ixx, iyy, izz, :, 14)
     1         - u ( ixx, iyy, izz, :, 15)
     1         + u ( ixx, iyy, izz, :, 16)
     1         - u ( ixx, iyy, izz, :, 17)
     1         + u ( ixx, iyy, izz, :, 18)
     1         + u ( ixx, iyy, izz, :, 19)

           Is5(:) = R5_xyz(:)**2 + a(:)

! -----

            R6_xyz(:) = u ( ixxm1, iyyp1, izzm1, :, 1)
     1         + 1.1 * ( u ( ixx, iyy, izz, :, 11)
     1                 - u ( ixx, iyy, izz, :, 12)
     1                 + u ( ixx, iyy, izz, :, 13))
     1         - u ( ixx, iyy, izz, :, 1)
     1         + u ( ixx, iyy, izz, :, 2)
     1         - u ( ixx, iyy, izz, :, 3)
     1         + u ( ixx, iyy, izz, :, 4)
     1         - u ( ixx, iyy, izz, :, 5)
     1         - u ( ixx, iyy, izz, :, 6)
     1         - u ( ixx, iyy, izz, :, 7)
     1         + u ( ixx, iyy, izz, :, 8)
     1         + u ( ixx, iyy, izz, :, 9)
     1         - u ( ixx, iyy, izz, :, 10)
     1         - u ( ixx, iyy, izz, :, 14)
     1         + u ( ixx, iyy, izz, :, 15)
     1         + u ( ixx, iyy, izz, :, 16)
     1         + u ( ixx, iyy, izz, :, 17)
     1         + u ( ixx, iyy, izz, :, 18)
     1         - u ( ixx, iyy, izz, :, 19)

           Is6(:) = R6_xyz(:)**2 + a(:)

! -----

            R7_xyz(:) = u ( ixxp1, iyym1, izzm1, :, 1)
     1         + 1.1 * (-u ( ixx, iyy, izz, :, 11)
     1                 + u ( ixx, iyy, izz, :, 12)
     1                 + u ( ixx, iyy, izz, :, 13))
     1         - u ( ixx, iyy, izz, :, 1)
     1         - u ( ixx, iyy, izz, :, 2)
     1         + u ( ixx, iyy, izz, :, 3)
     1         + u ( ixx, iyy, izz, :, 4)
     1         - u ( ixx, iyy, izz, :, 5)
     1         - u ( ixx, iyy, izz, :, 6)
     1         - u ( ixx, iyy, izz, :, 7)
     1         + u ( ixx, iyy, izz, :, 8)
     1         - u ( ixx, iyy, izz, :, 9)
     1         + u ( ixx, iyy, izz, :, 10)
     1         + u ( ixx, iyy, izz, :, 14)
     1         + u ( ixx, iyy, izz, :, 15)
     1         - u ( ixx, iyy, izz, :, 16)
     1         + u ( ixx, iyy, izz, :, 17)
     1         - u ( ixx, iyy, izz, :, 18)
     1         + u ( ixx, iyy, izz, :, 19)

           Is7(:) = R7_xyz(:)**2 + a(:)

! -----

            R8_xyz(:) = - u ( ixxm1, iyym1, izzm1, :, 1)
     1         - 1.1 * ( u ( ixx, iyy, izz, :, 11)
     1                 + u ( ixx, iyy, izz, :, 12)
     1                 + u ( ixx, iyy, izz, :, 13))
     1         + u ( ixx, iyy, izz, :, 1)
     1         - u ( ixx, iyy, izz, :, 2)
     1         - u ( ixx, iyy, izz, :, 3)
     1         - u ( ixx, iyy, izz, :, 4)
     1         + u ( ixx, iyy, izz, :, 5)
     1         + u ( ixx, iyy, izz, :, 6)
     1         + u ( ixx, iyy, izz, :, 7)
     1         + u ( ixx, iyy, izz, :, 8)
     1         + u ( ixx, iyy, izz, :, 9)
     1         + u ( ixx, iyy, izz, :, 10)
     1         - u ( ixx, iyy, izz, :, 14)
     1         - u ( ixx, iyy, izz, :, 15)
     1         - u ( ixx, iyy, izz, :, 16)
     1         - u ( ixx, iyy, izz, :, 17)
     1         - u ( ixx, iyy, izz, :, 18)
     1         - u ( ixx, iyy, izz, :, 19)

           Is8(:) = R8_xyz(:)**2 + a(:)

! -----

            Wt1(:) = 1.0 / Is1(:)**4
            Wt2(:) = 1.0 / Is2(:)**4
            Wt3(:) = 1.0 / Is3(:)**4
            Wt4(:) = 1.0 / Is4(:)**4
            Wt5(:) = 1.0 / Is5(:)**4
            Wt6(:) = 1.0 / Is6(:)**4
            Wt7(:) = 1.0 / Is7(:)**4
            Wt8(:) = 1.0 / Is8(:)**4

            b(:) = 1.0 / ( Wt1(:) + Wt2(:) + Wt3(:) + Wt4(:)
     1                   + Wt5(:) + Wt6(:) + Wt7(:) + Wt8(:) )

            Wt1(:) = Wt1(:) * b(:)
            Wt2(:) = Wt2(:) * b(:)
            Wt3(:) = Wt3(:) * b(:)
            Wt4(:) = Wt4(:) * b(:)
            Wt5(:) = Wt5(:) * b(:)
            Wt6(:) = Wt6(:) * b(:)
            Wt7(:) = Wt7(:) * b(:)
            Wt8(:) = Wt8(:) * b(:)

            u ( ixx, iyy, izz, :, 20) =
     1        ( Wt1(:) * R1_xyz(:) + Wt2(:) * R2_xyz(:)
     1        + Wt3(:) * R3_xyz(:) + Wt4(:) * R4_xyz(:)
     1        + Wt5(:) * R5_xyz(:) + Wt6(:) * R6_xyz(:)
     1        + Wt7(:) * R7_xyz(:) + Wt8(:) * R8_xyz(:)) * c

! --------------------

          END IF   ! End of "IF ( lo_order_zone_avg .EQ. 1) THEN"

! ----------------------------------------------------------------------

          END DO ! End of DO ixx = loop.

        END DO ! End of DO iyy = loop.

      END DO ! End of DO izz = loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END SUBROUTINE LIMIT_RECONST_ALL_ZC_P3

#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE v_polynomial_modes (
      
     1           x, y, z, n_cc_modes, u_coef,
     1           v_polynomial )
     
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER n_cc_modes
      
      REAL a, p1, p2, p3, p4, p5, x, y, z


      REAL, DIMENSION ( n_cc_modes) :: u_coef
      REAL :: v_polynomial

! ----------------------------------------------------------------------
! Variables past this point are internal to this subroutine.
   
      p1(a ) = a
      p2(a ) = a**2 - (1.0 / 12.0)
      p3(a ) = a**3 - (3.0 / 20.0) * a
      p4(a ) = a**4 - (3.0 / 14.0) * a**2 + (3.0 / 560.0)
      p5(a ) = a**5 - (5.0 / 18.0) * a**3 + (5.0 / 336.0) * a
      
! ------------------------------

#if ( NUMORDER >= 2)

      v_polynomial = 
     1   u_coef ( 01) + u_coef ( 02) * p1(x) + 
     1   u_coef ( 03) * p1(y) + u_coef ( 04) * p1(z)

#endif
! End of "NUMORDER >= 2"
! ------------------------------
#if ( NUMORDER >= 3)

      v_polynomial = v_polynomial +
     
     1   u_coef ( 05) * p2(x) + u_coef ( 06) * p2(y) +
     1   u_coef ( 07) * p2(z) + 
     1   u_coef ( 08) * p1(x) * p1(y) + 
     1   u_coef ( 09) * p1(y) * p1(z) + 
     1   u_coef ( 10) * p1(x) * p1(z)

#endif
! End of "NUMORDER >= 3"
! ------------------------------
#if ( NUMORDER >= 4)

      v_polynomial = v_polynomial +
     
     1   u_coef ( 11) * p3(x) + u_coef ( 12) * p3(y) +
     1   u_coef ( 13) * p3(z) + 
     1   u_coef ( 14) * p2(x) * p1(y) + 
     1   u_coef ( 15) * p2(x) * p1(z) +
     1   u_coef ( 16) * p1(x) * p2(y) + 
     1   u_coef ( 17) * p2(y) * p1(z) +
     1   u_coef ( 18) * p1(x) * p2(z) +
     1   u_coef ( 19) * p1(y) * p2(z) + 
     1   u_coef ( 20) * p1(x) * p1(y) * p1(z)

#endif
! End of "NUMORDER >= 4"
! ------------------------------
#if ( NUMORDER >= 5)
      
      v_polynomial = v_polynomial + 
        
     1   u_coef ( 21) * p4(x) + u_coef ( 22) * p4(y) + 
     1   u_coef ( 23) * p4(z) +
     1   u_coef ( 24) * p3(x) * p1(y) + 
     1   u_coef ( 25) * p3(x) * p1(z) +
     1   u_coef ( 26) * p1(x) * p3(y) + 
     1   u_coef ( 27) * p3(y) * p1(z) +
     1   u_coef ( 28) * p1(x) * p3(z) + 
     1   u_coef ( 29) * p1(y) * p3(z) +
     1   u_coef ( 30) * p2(x) * p2(y) + 
     1   u_coef ( 31) * p2(y) * p2(z) + 
     1   u_coef ( 32) * p2(x) * p2(z) +
     1   u_coef ( 33) * p2(x) * p1(y) * p1(z) + 
     1   u_coef ( 34) * p1(x) * p2(y) * p1(z) + 
     1   u_coef ( 35) * p1(x) * p1(y) * p2(z)
     
#endif
! End of "NUMORDER >= 5"
! ------------------------------
#if ( NUMORDER >= 6)

      v_polynomial = v_polynomial + 
     
     1   u_coef ( 36) * p5(x) + u_coef ( 37) * p5(y) +
     1   u_coef ( 38) * p5(z) + 
     1   u_coef ( 39) * p4(x) * p1(y) + u_coef ( 40) * p4(x) * p1(z) +
     1   u_coef ( 41) * p1(x) * p4(y) + u_coef ( 42) * p4(y) * p1(z) +
     1   u_coef ( 43) * p1(x) * p4(z) + u_coef ( 44) * p1(y) * p4(z) +
     1   u_coef ( 45) * p3(x) * p2(y) + u_coef ( 46) * p3(x) * p2(z) +
     1   u_coef ( 47) * p2(x) * p3(y) + u_coef ( 48) * p3(y) * p2(z) +
     1   u_coef ( 49) * p2(x) * p3(z) + u_coef ( 50) * p2(y) * p3(z) +
     1   u_coef ( 51) * p3(x) * p1(y) * p1(z) + 
     1   u_coef ( 52) * p1(x) * p3(y) * p1(z) + 
     1   u_coef ( 53) * p1(x) * p1(y) * p3(z) + 
     1   u_coef ( 54) * p2(x) * p2(y) * p1(z) + 
     1   u_coef ( 55) * p1(x) * p2(y) * p2(z) + 
     1   u_coef ( 56) * p2(x) * p1(y) * p2(z)

#endif
! End of "NUMORDER >= 6"

! ----------------------------------------------------------------------

      END SUBROUTINE v_polynomial_modes

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE CMat_Row_Coefs ( NUnknown, x, y, z, nx, ny, nz,
     1                                      CMat_Row)

! ----------------------------------------------------------------------
! Written by : Dinshaw Balsara

! Takes in the position "x, y, z" at a face as well as the normal to it,
! "nx, ny, nz". Builds the coefficients for that row. Typically, "NUnknown = 3 * n_cc_modes"
! Because we have three components for the magnetic field.
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NUnknown

      REAL x, y, z, nx, ny, nz, CMat_Row ( NUnknown)

! ----------------------------------------------------------------------
#if ( NUMORDER == 2)

      CMat_Row ( 1) = nx ; CMat_Row ( 2) = nx * x
      CMat_Row ( 3) = nx * y ; CMat_Row ( 4) = nx * z

      CMat_Row ( 5) = ny ; CMat_Row ( 6) = ny * x
      CMat_Row ( 7) = ny * y ; CMat_Row ( 8) = ny * z

      CMat_Row ( 9) = nz ; CMat_Row ( 10) = nz * x
      CMat_Row ( 11) = nz * y ; CMat_Row ( 12) = nz * z

#endif
! End of "( NUMORDER == 2)".
! ----------------------------------------------------------------------
#if ( NUMORDER == 3)
    
      CMat_Row ( 1) = nx ; CMat_Row ( 2) = nx * x
      CMat_Row ( 3) = nx * y ; CMat_Row ( 4) = nx * z
      CMat_Row ( 5) = nx * ( x**2 - 1.0 / 12.0) 
      CMat_Row ( 6) = nx * ( y**2 - 1.0 / 12.0)
      CMat_Row ( 7) = nx * ( z**2 - 1.0 / 12.0) 
      CMat_Row ( 8) = nx * x * y
      CMat_Row ( 9) = nx * y * z 
      CMat_Row ( 10) = nx * x * z

      CMat_Row ( 11) = ny ; CMat_Row ( 12) = ny * x
      CMat_Row ( 13) = ny * y ; CMat_Row ( 14) = ny * z
      CMat_Row ( 15) = ny * ( x**2 - 1.0 / 12.0)  
      CMat_Row ( 16) = ny * ( y**2 - 1.0 / 12.0)
      CMat_Row ( 17) = ny * ( z**2 - 1.0 / 12.0)
      CMat_Row ( 18) = ny * x * y
      CMat_Row ( 19) = ny * y * z  
      CMat_Row ( 20) = ny * x * z

      CMat_Row ( 21) = nz ; CMat_Row ( 22) = nz * x
      CMat_Row ( 23) = nz * y ; CMat_Row ( 24) = nz * z
      CMat_Row ( 25) = nz * ( x**2 - 1.0 / 12.0)
      CMat_Row ( 26) = nz * ( y**2 - 1.0 / 12.0)
      CMat_Row ( 27) = nz * ( z**2 - 1.0 / 12.0)
      CMat_Row ( 28) = nz * x * y
      CMat_Row ( 29) = nz * y * z
      CMat_Row ( 30) = nz * x * z

#endif
! End of "( NUMORDER == 3)".
! ----------------------------------------------------------------------
#if ( NUMORDER == 4)
      
      CMat_Row ( 1) = nx ; CMat_Row ( 2) = nx * x
      CMat_Row ( 3) = nx * y ; CMat_Row ( 4) = nx * z
      CMat_Row ( 5) = nx * ( x**2 - 1.0 / 12.0)
      CMat_Row ( 6) = nx * ( y**2 - 1.0 / 12.0)
      CMat_Row ( 7) = nx * ( z**2 - 1.0 / 12.0)
      CMat_Row ( 8) = nx * x * y
      CMat_Row ( 9) = nx * y * z
      CMat_Row ( 10) = nx * x * z
      CMat_Row ( 11) = nx * ( x**3 - (3.0 / 20.0) * x)
      CMat_Row ( 12) = nx * ( y**3 - (3.0 / 20.0) * y)
      CMat_Row ( 13) = nx * ( z**3 - (3.0 / 20.0) * z)
      CMat_Row ( 14) = nx * x**2 * y
      CMat_Row ( 15) = nx * x**2 * z
      CMat_Row ( 16) = nx * x * y**2
      CMat_Row ( 17) = nx * y**2 * z
      CMat_Row ( 18) = nx * x * z**2
      CMat_Row ( 19) = nx * y * z**2
      CMat_Row ( 20) = nx * x * y * z

      CMat_Row ( 21) = ny ; CMat_Row ( 22) = ny * x
      CMat_Row ( 23) = ny * y ; CMat_Row ( 24) = ny * z
      CMat_Row ( 25) = ny * ( x**2 - 1.0 / 12.0)
      CMat_Row ( 26) = ny * ( y**2 - 1.0 / 12.0)
      CMat_Row ( 27) = ny * ( z**2 - 1.0 / 12.0)
      CMat_Row ( 28) = ny * x * y
      CMat_Row ( 29) = ny * y * z
      CMat_Row ( 30) = ny * x * z
      CMat_Row ( 31) = ny * ( x**3 - (3.0 / 20.0) * x)
      CMat_Row ( 32) = ny * ( y**3 - (3.0 / 20.0) * y)
      CMat_Row ( 33) = ny * ( z**3 - (3.0 / 20.0) * z)
      CMat_Row ( 34) = ny * x**2 * y
      CMat_Row ( 35) = ny * x**2 * z
      CMat_Row ( 36) = ny * x * y**2
      CMat_Row ( 37) = ny * y**2 * z
      CMat_Row ( 38) = ny * x * z**2
      CMat_Row ( 39) = ny * y * z**2
      CMat_Row ( 40) = ny * x * y * z

      CMat_Row ( 41) = nz ; CMat_Row ( 42) = nz * x
      CMat_Row ( 43) = nz * y ; CMat_Row ( 44) = nz * z
      CMat_Row ( 45) = nz * ( x**2 - 1.0 / 12.0)
      CMat_Row ( 46) = nz * ( y**2 - 1.0 / 12.0)
      CMat_Row ( 47) = nz * ( z**2 - 1.0 / 12.0)
      CMat_Row ( 48) = nz * x * y
      CMat_Row ( 49) = nz * y * z
      CMat_Row ( 50) = nz * x * z
      CMat_Row ( 51) = nz * ( x**3 - (3.0 / 20.0) * x)
      CMat_Row ( 52) = nz * ( y**3 - (3.0 / 20.0) * y)
      CMat_Row ( 53) = nz * ( z**3 - (3.0 / 20.0) * z)
      CMat_Row ( 54) = nz * x**2 * y
      CMat_Row ( 55) = nz * x**2 * z
      CMat_Row ( 56) = nz * x * y**2
      CMat_Row ( 57) = nz * y**2 * z
      CMat_Row ( 58) = nz * x * z**2
      CMat_Row ( 59) = nz * y * z**2
      CMat_Row ( 60) = nz * x * y * z
     
#endif
! End of "( NUMORDER == 4)".
! ----------------------------------------------------------------------

      END SUBROUTINE CMat_Row_Coefs

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

