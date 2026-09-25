!  cp BuildGeomesh_c3_src a.F; gfortran -o x a.F -fdefault-real-8; ./x;
!  vi test_quadpoints.out



#include <cubesphere/200mesh/201mesh_sector_general.f>
#include <cubesphere/200mesh/202generate_global_radial_mesh.f>
#include <cubesphere/200mesh/203fill_ghost2d_global.f>
#include <cubesphere/200mesh/204build_mode_indices.f>
#include <cubesphere/200mesh/210extract_vert_coords_patch.f>
#include <cubesphere/200mesh/211stencil_indexing_tricorners.f>
#include <cubesphere/200mesh/220coords_to_wt_quadpts_faces_iso2.f>
#include <cubesphere/200mesh/230coords_to_wt_quadpts_bfld_iso2.f>




      PROGRAM Test_QuadPoints

! ----------------------------------------------------------------------

! Tests "Mesh_Sector" for its ability to set up the mesh for a single sector.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ntp_global, nr_global, NGhost_global,
     1        ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1        isector, ix_patch, iy_patch, iz_patch,
     1        ilevel, ilevel_max, iscale,

     1        n_cc_modes,
     1        n_cc_components, n_fc_modes,

     1        itri_patch_LD, itri_patch_RD, itri_patch_RU,
     1        itri_patch_LU,

     1        istep, ntstep

      PARAMETER ( ntp_zones = 1, nr_zones = 1,
     1            ntp_patches = 1, nr_patches = 1, NGhost = 0,

     1            ntstep = 1,

     1            ilevel = 2, ilevel_max = 3,
     1            iscale = 2**(ilevel_max - ilevel),

     1            ntp_global = ntp_zones * ntp_patches * iscale,
     1            nr_global = nr_zones * nr_patches * iscale,
     1            NGhost_global = NGhost * iscale,

     1            n_cc_components = 8,
     1            n_fc_modes = 1,

#if ( NUMORDER == 2)

     1            n_cc_modes = 4

#elif ( NUMORDER == 3)

     1            n_cc_modes = 10

#elif ( NUMORDER == 4)

     1            n_cc_modes = 20

#endif

     1            )


      REAL rmin, rmax, pi, smallnum, radial_mesh_ratio

      REAL sector_vertices ( 4, 3, 6)

      INTEGER ixx, iyy, izz, i

! ----------------------------------------------------------------------

      REAL coordinates_gr_global
     1     ( 0 - NGhost_global: ntp_global + NGhost_global,
     1       0 - NGhost_global: ntp_global + NGhost_global, 3, 6)

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                          1 - NGhost: ntp_zones + NGhost,
     1                          1 - NGhost: nr_zones + NGhost)

      REAL rb_save_global ( 0 - NGhost_global:
     1                      nr_global + NGhost_global)

      REAL, DIMENSION ( 1  - NGhost_global: nr_global + NGhost_global)::
     1     dr_save_global, rc_save_global

! ----------------------------------------------------------------------

      INTEGER, DIMENSION ( n_cc_modes) ::
     1         ii_mode, jj_mode, kk_mode

      REAL Coef_Basis_Moments ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes)

! ------------------------------

      INTEGER istencil_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                         1 - NGhost: ntp_zones + NGhost, 2, 12)

! ----------------------------------------------------------------------

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

! ----------------------------------------------------------------------

      INTEGER ixx1, iyy1, ixx2, iyy2, ixx01, iyy01, ixx02, iyy02,
     1        itri_corner

      REAL timing_info1, timing_info2, time_vec ( ntstep),
     1     n1 ( 3), n2 ( 3), n3 ( 3), tempaa, angle, xec, yec, zec,
     1     n1_mag, n2_mag

! ----------------------------------------------------------------------

      OPEN ( UNIT = 6, FILE = "test_quadpoints.out")

      isector = 1     ! local and global x,y,z directions coincide

      smallnum = 1.0e-12
      pi = 4.0 * ATAN ( 1.0)

      ix_patch = 1
      iy_patch = 1
      iz_patch = 1

      isector = 1

      rmin = 10.0
      rmax = 11.0

      radial_mesh_ratio = 1.0

! -----------

      CALL GENERATE_GLOBAL_RADIAL_LOG_MESH

     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global)

! ----------------------------------------------------------------------
! Build the powers for the polynomials so that we can zip through them
! automatically during stencil construction.

      CALL BUILD_MODE_INDICES

     1         ( n_cc_modes,

     1           ii_mode, jj_mode, kk_mode)

! ----------------------------------------------------------------------
! Build "sector_vertices" for all 6 sectors. Also build a global coordinate
! system on the UNIT sphere and store it in 2D mesh "coordinates_gr_global".

      coordinates_gr_global = 0.0

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL Mesh_Sector_General ( i, 1,
     1                   1, 1, ntp_global, NGhost_global,
     1sector_vertices ( :, :, i), coordinates_gr_global ( :, :, :, i) )

      END DO

! ----------------------------------------------------------------------
! Extract a patch and save it in "vert_coords_patch". Also make a mask for the
! active zones, "imask_active_patch", set to 1 for active zones. Make sure that
! "isector, ix_patch, iy_patch, iz_patch" and "ilevel, ilevel_max" have
! been set correctly.

      CALL EXTRACT_VERT_COORDS_PATCH

     1       ( ntp_global, nr_global, NGhost_global,
     1         ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         isector, ix_patch, iy_patch, iz_patch,
     1         ilevel, ilevel_max,

     1         coordinates_gr_global, rb_save_global,

     1         vert_coords_patch, imask_active_patch)

! ----------------------------------------------------------------------

! At the tri-corners of a cube sphere mesh, we cannot index with "ixx + 1"
! etc. This subroutine works out the stencils for the zones that are one or two
! layers from a tri-corner. It maps the stencil back into the physical part of
! the mesh.

      CALL STENCIL_INDEXING_TRICORNERS

     1       ( ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         ix_patch, iy_patch,

     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,

     1    imask_active_patch, istencil_patch)

! ------------------------------

      CALL COORDS_TO_WT_QUADPTS_FACES_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length)

! Build the quadrature points, weights and normals for the integration
! of magnetic field in the Multi-D Riemann solver.

      CALL COORDS_TO_WT_QUADPTS_BFLD_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)


! --------------------

      WRITE ( 6, *)"vert_coords_patch: "

      WRITE ( 6, '(A26, I1, A4)', advance="no")"(0, 0, 0) = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') vert_coords_patch ( 0,0,0,1),
     1                                    vert_coords_patch ( 0,0,0,2),
     1                                    vert_coords_patch ( 0,0,0,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"(1, 0, 0) = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') vert_coords_patch ( 1,0,0,1),
     1                                    vert_coords_patch ( 1,0,0,2),
     1                                    vert_coords_patch ( 1,0,0,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"(0, 1, 0) = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') vert_coords_patch ( 0,1,0,1),
     1                                    vert_coords_patch ( 0,1,0,2),
     1                                    vert_coords_patch ( 0,1,0,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"(1, 1, 0) = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') vert_coords_patch ( 1,1,0,1),
     1                                    vert_coords_patch ( 1,1,0,2),
     1                                    vert_coords_patch ( 1,1,0,3)

! --------------------

      WRITE ( 6, '(A26, I1, A4)', advance="no")"(0, 0, 1) = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') vert_coords_patch ( 0,0,1,1),
     1                                    vert_coords_patch ( 0,0,1,2),
     1                                    vert_coords_patch ( 0,0,1,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"(1, 0, 1) = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') vert_coords_patch ( 1,0,1,1),
     1                                    vert_coords_patch ( 1,0,1,2),
     1                                    vert_coords_patch ( 1,0,1,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"(0, 1, 1) = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') vert_coords_patch ( 0,1,1,1),
     1                                    vert_coords_patch ( 0,1,1,2),
     1                                    vert_coords_patch ( 0,1,1,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"(1, 1, 1) = "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') vert_coords_patch ( 1,1,1,1),
     1                                    vert_coords_patch ( 1,1,1,2),
     1                                    vert_coords_patch ( 1,1,1,3)

! --------------------
#if( NUMORDER == 2)

      WRITE ( 6, *)" "
      WRITE ( 6, *)"xb_quadpts_patch, face 0 (-x), Points:  ",
     1             "1 at +y,  2 at -y,  3 at +z,  4 at -z "

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 1, y+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') xb_quadpts_patch ( 0,1,1,1,1),
     1                                    xb_quadpts_patch ( 0,1,1,2,1),
     1                                    xb_quadpts_patch ( 0,1,1,3,1)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 2, y- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') xb_quadpts_patch ( 0,1,1,1,2),
     1                                    xb_quadpts_patch ( 0,1,1,2,2),
     1                                    xb_quadpts_patch ( 0,1,1,3,2)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 3, z+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') xb_quadpts_patch ( 0,1,1,1,3),
     1                                    xb_quadpts_patch ( 0,1,1,2,3),
     1                                    xb_quadpts_patch ( 0,1,1,3,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 4, z- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') xb_quadpts_patch ( 0,1,1,1,4),
     1                                    xb_quadpts_patch ( 0,1,1,2,4),
     1                                    xb_quadpts_patch ( 0,1,1,3,4)


      WRITE ( 6, *)" "
      WRITE ( 6, *)"xb_quadpts_patch, face 1 (+x), Points:  ",
     1             "1 at +y,  2 at -y,  3 at +z,  4 at -z "

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 1, y+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') xb_quadpts_patch ( 1,1,1,1,1),
     1                                    xb_quadpts_patch ( 1,1,1,2,1),
     1                                    xb_quadpts_patch ( 1,1,1,3,1)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 2, y- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') xb_quadpts_patch ( 1,1,1,1,2),
     1                                    xb_quadpts_patch ( 1,1,1,2,2),
     1                                    xb_quadpts_patch ( 1,1,1,3,2)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 3, z+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') xb_quadpts_patch ( 1,1,1,1,3),
     1                                    xb_quadpts_patch ( 1,1,1,2,3),
     1                                    xb_quadpts_patch ( 1,1,1,3,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 4, z- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') xb_quadpts_patch ( 1,1,1,1,4),
     1                                    xb_quadpts_patch ( 1,1,1,2,4),
     1                                    xb_quadpts_patch ( 1,1,1,3,4)


! --------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"yb_quadpts_patch, face 0 (-y), Points:  ",
     1             "1 at +z,  2 at -z,  3 at +x,  4 at -x "

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 1, y+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') yb_quadpts_patch ( 1,0,1,1,1),
     1                                    yb_quadpts_patch ( 1,0,1,2,1),
     1                                    yb_quadpts_patch ( 1,0,1,3,1)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 2, y- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') yb_quadpts_patch ( 1,0,1,1,2),
     1                                    yb_quadpts_patch ( 1,0,1,2,2),
     1                                    yb_quadpts_patch ( 1,0,1,3,2)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 3, z+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') yb_quadpts_patch ( 1,0,1,1,3),
     1                                    yb_quadpts_patch ( 1,0,1,2,3),
     1                                    yb_quadpts_patch ( 1,0,1,3,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 4, z- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') yb_quadpts_patch ( 1,0,1,1,4),
     1                                    yb_quadpts_patch ( 1,0,1,2,4),
     1                                    yb_quadpts_patch ( 1,0,1,3,4)


      WRITE ( 6, *)" "
      WRITE ( 6, *)"yb_quadpts_patch, face 1 (+y), Points:  ",
     1             "1 at +z,  2 at -z,  3 at +x,  4 at -x "

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 1, y+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') yb_quadpts_patch ( 1,1,1,1,1),
     1                                    yb_quadpts_patch ( 1,1,1,2,1),
     1                                    yb_quadpts_patch ( 1,1,1,3,1)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 2, y- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') yb_quadpts_patch ( 1,1,1,1,2),
     1                                    yb_quadpts_patch ( 1,1,1,2,2),
     1                                    yb_quadpts_patch ( 1,1,1,3,2)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 3, z+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') yb_quadpts_patch ( 1,1,1,1,3),
     1                                    yb_quadpts_patch ( 1,1,1,2,3),
     1                                    yb_quadpts_patch ( 1,1,1,3,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 4, z- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') yb_quadpts_patch ( 1,1,1,1,4),
     1                                    yb_quadpts_patch ( 1,1,1,2,4),
     1                                    yb_quadpts_patch ( 1,1,1,3,4)


! --------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"zb_quadpts_patch, face 0 (-z), Points:  ",
     1             "1 at +x,  2 at -x,  3 at +y,  4 at -y "

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 1, x+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') zb_quadpts_patch ( 1,1,0,1,1),
     1                                    zb_quadpts_patch ( 1,1,0,2,1),
     1                                    zb_quadpts_patch ( 1,1,0,3,1)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 2, x- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') zb_quadpts_patch ( 1,1,0,1,2),
     1                                    zb_quadpts_patch ( 1,1,0,2,2),
     1                                    zb_quadpts_patch ( 1,1,0,3,2)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 3, y+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') zb_quadpts_patch ( 1,1,0,1,3),
     1                                    zb_quadpts_patch ( 1,1,0,2,3),
     1                                    zb_quadpts_patch ( 1,1,0,3,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 4, y- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') zb_quadpts_patch ( 1,1,0,1,4),
     1                                    zb_quadpts_patch ( 1,1,0,2,4),
     1                                    zb_quadpts_patch ( 1,1,0,3,4)


      WRITE ( 6, *)" "
      WRITE ( 6, *)"zb_quadpts_patch, face 1 (+z), Points:  ",
     1             "1 at +x,  2 at -x,  3 at +y,  4 at -y "

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 1, x+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') zb_quadpts_patch ( 1,1,1,1,1),
     1                                    zb_quadpts_patch ( 1,1,1,2,1),
     1                                    zb_quadpts_patch ( 1,1,1,3,1)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 2, x- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') zb_quadpts_patch ( 1,1,1,1,2),
     1                                    zb_quadpts_patch ( 1,1,1,2,2),
     1                                    zb_quadpts_patch ( 1,1,1,3,2)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 3, y+ : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') zb_quadpts_patch ( 1,1,1,1,3),
     1                                    zb_quadpts_patch ( 1,1,1,2,3),
     1                                    zb_quadpts_patch ( 1,1,1,3,3)

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Pt 4, y- : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') zb_quadpts_patch ( 1,1,1,1,4),
     1                                    zb_quadpts_patch ( 1,1,1,2,4),
     1                                    zb_quadpts_patch ( 1,1,1,3,4)

#endif


! --------------------
#if( NUMORDER == 3)

      WRITE ( 6, *)" "
      WRITE ( 6, *)"xb_quadpts_patch, face 0 (-x), Points in sequence:",
     1             " -y -> +y and then -z to +z"

      DO i = 1, NBFLD_QUAD_PTS

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Point : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') xb_quadpts_patch ( 0,1,1,1,i),
     1                                    xb_quadpts_patch ( 0,1,1,2,i),
     1                                    xb_quadpts_patch ( 0,1,1,3,i)

      END DO

! ----------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"xb_quadpts_patch, face 1 (+x), Points in sequence:",
     1             " -y -> +y and then -z to +z"

      DO i = 1, NBFLD_QUAD_PTS

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Point : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') xb_quadpts_patch ( 1,1,1,1,i),
     1                                    xb_quadpts_patch ( 1,1,1,2,i),
     1                                    xb_quadpts_patch ( 1,1,1,3,i)

      END DO

! --------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"yb_quadpts_patch, face 0 (-y), Points in sequence:",
     1             " -z -> +z and then -x to +x"

      DO i = 1, NBFLD_QUAD_PTS

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Point : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') yb_quadpts_patch ( 1,0,1,1,i),
     1                                    yb_quadpts_patch ( 1,0,1,2,i),
     1                                    yb_quadpts_patch ( 1,0,1,3,i)

      END DO


      WRITE ( 6, *)" "
      WRITE ( 6, *)"yb_quadpts_patch, face 1 (+y), Points in sequence:",
     1             " -z -> +z and then -x to +x"

      DO i = 1, NBFLD_QUAD_PTS

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Point : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') yb_quadpts_patch ( 1,1,1,1,i),
     1                                    yb_quadpts_patch ( 1,1,1,2,i),
     1                                    yb_quadpts_patch ( 1,1,1,3,i)

      END DO

! --------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"zb_quadpts_patch, face 0 (-z), Points in sequence:",
     1             " -x -> +x and then -y to +y"

      DO i = 1, NBFLD_QUAD_PTS

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Point : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') zb_quadpts_patch ( 1,1,0,1,i),
     1                                    zb_quadpts_patch ( 1,1,0,2,i),
     1                                    zb_quadpts_patch ( 1,1,0,3,i)

      END DO


      WRITE ( 6, *)" "
      WRITE ( 6, *)"zb_quadpts_patch, face 1 (+z), Points in sequence:",
     1             " -x -> +x and then -y to +y"

      DO i = 1, NBFLD_QUAD_PTS

      WRITE ( 6, '(A26, I1, A4)', advance="no")"Point : "
      WRITE ( 6, '(f10.4, f10.4, f10.4)') zb_quadpts_patch ( 1,1,1,1,i),
     1                                    zb_quadpts_patch ( 1,1,1,2,i),
     1                                    zb_quadpts_patch ( 1,1,1,3,i)

      END DO

#endif

! ----------------------------------------------------------------------

      END PROGRAM Test_QuadPoints

