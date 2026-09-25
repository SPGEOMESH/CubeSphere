




      PROGRAM TEST_CUBE_SPHERE_Mesh

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER ntp_global, nr_global, NGhost_global,
     1        ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1        isector, ix_patch, iy_patch, iz_patch,
     1        ilevel, ilevel_max, iscale,

     1        n_cc_modes,

     1        itri_patch_LD, itri_patch_RD, itri_patch_RU,
     1        itri_patch_LU,

     1        istep, ntstep

      PARAMETER ( ntp_zones = 10, nr_zones = 10,
     1            ntp_patches = 1, nr_patches = 1, NGhost = 4,

     1            ntstep = 1,

     1            ilevel = 2, ilevel_max = 3,
     1            iscale = 2**(ilevel_max - ilevel),

     1            ntp_global = ntp_zones * ntp_patches * iscale,
     1            nr_global = nr_zones * nr_patches * iscale,
     1            NGhost_global = NGhost * iscale,
#if ( NUMORDER == 2)

     1            n_cc_modes = 4

#elif ( NUMORDER == 3)

     1            n_cc_modes = 10

#elif ( NUMORDER == 4)

     1            n_cc_modes = 20

#endif

     1            )

      REAL rmin, rmax

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

! This is just a tester program. Here we must put "ntp_patches = 1" and
! "nr_patches = 1" and "NGhost = 0".
! ----------------------------------------------------------------------

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

! ----------

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

! ----------

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

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ------------------------------

      REAL v_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER * NUMORDER)

      REAL v_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER * NUMORDER)

      REAL v_quadtangent_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, 3, NUMORDER * NUMORDER * NUMORDER)

      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------

      INTEGER, DIMENSION ( n_cc_modes) ::
     1         ii_mode, jj_mode, kk_mode

      REAL Coef_Basis_Moments ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes)

! ------------------------------

      INTEGER istencil_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                         1 - NGhost: ntp_zones + NGhost, 2, 12)

! ----------------------------------------------------------------------

      INTEGER ixx1, iyy1, ixx2, iyy2, ixx01, iyy01, ixx02, iyy02,
     1        itri_corner

      REAL timing_info1, timing_info2, time_vec ( ntstep),
     1     n1 ( 3), n2 ( 3), n3 ( 3), tempaa, angle, xec, yec, zec,
     1     n1_mag, n2_mag

! ----------------------------------------------------------------------
! Create arrays that only need to reside on the GPU. These arrays will
! never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

!$acc data create(xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,&
!$acc&           xf_area_patch, xf_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch, &
!$acc&           yf_area_patch, yf_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch, &
!$acc&           zf_area_patch, zf_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           v_scale_length,                                        &
!$acc&                                                                  &
!$acc&           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,  &
!$acc&           v_volume_patch, v_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           Coef_Basis_Moments, istencil_patch, time_vec)

! ----------------------------------------------------------------------
! Build a 1D mesh in the radial direction.

      OPEN ( UNIT = 6, FILE = "test_cube_sphere_mesh.out")

      smallnum = 1.0e-12
      pi = 4.0 * ATAN ( 1.0)

      ix_patch = 1
      iy_patch = 1
      iz_patch = 1

      isector = 1

      rmin = 1.0
      rmax = 10.0

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
! Previous call to "Mesh_Sector_General" will only fill the interior of
! each sector. The ghost zone information will have to be imported from
! the neighboring sectors.

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL FILL_GHOST2D_GLOBAL

     1           ( i, ntp_global, NGhost_global,
     1             coordinates_gr_global)

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

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
#define VERIFY_ANGLE_AT_TRICORNER 1

! ----------------------------------------------------------------------
#if( VERIFY_ANGLE_AT_TRICORNER == 1)
! ------------------------------------------------------------

          ixx = 0
          iyy = 0 !ntp_zones
          izz = nr_zones/2

! Find the edge-center "xec, yec, zec" by bisecting the two z-vertices.
! No need to rescale by radius since this is in the radial direction by definition.

          xec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 1)  ! Bisect the z-vertices.
     1                + vert_coords_patch ( ixx, iyy, izz - 1, 1) )
          yec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 2)
     1                + vert_coords_patch ( ixx, iyy, izz - 1, 2) )
          zec = 0.5 * ( vert_coords_patch ( ixx, iyy, izz, 3)
     1                + vert_coords_patch ( ixx, iyy, izz - 1, 3) )

! ------------------------------------------------------------

          tempaa = 1.0 / SQRT ( xec**2 + yec**2 + zec**2)

          n3 ( 1) = xec * tempaa
          n3 ( 2) = yec * tempaa
          n3 ( 3) = zec * tempaa

! ------------------------------------------------------------

          WRITE(6,*)" "


! ------------------------------------------------------------
! Left-Down Tri-corner.

          IF ( ( itri_patch_LD .EQ. 1) .AND. ( ixx .EQ. 0) .AND.
     1         ( iyy .EQ. 0) ) THEN
! --------------------------------------------------
! "RU --> U1; LU --> U2; RD --> U3" to be used in the 2D Riemann solver
! with 3 states.

          WRITE(6,*)"itri_patch_LD .EQ. 1"

!          U_1 = U_RU
!          U_2 = U_LU
!          U_3 = U_RD

! ------------------------------

! set the tri-corner flag to "1", to call the 3-State RS.
          itri_corner = 1

! "n1" goes from (0, 0) to (1, 0). "n2" goes from (0, 0) to (0, 1)

          ixx01 = 0 ;  iyy01 = 0 ;  ixx1 = 1 ;  iyy1 = 0  ! Defines gradient for "n1".

          ixx02 = 0 ;  iyy02 = 0 ;  ixx2 = 0 ;  iyy2 = 1  ! Defines gradient for "n2".

! --------------------------------------------------
! Right-Down Tri-corner.

          ELSE IF ( ( itri_patch_RD .EQ. 1) .AND. ( ixx .EQ. ntp_zones)
     1         .AND. ( iyy .EQ. 0) ) THEN
! --------------------------------------------------
! "RU --> U1; LU --> U2; LD --> U3" to be used in the 2D Riemann solver
! with 3 states.

          WRITE(6,*)"itri_patch_RD .EQ. 1"

!          U_1 = U_RU
!          U_2 = U_LU
!          U_3 = U_LD

! ------------------------------

! set the tri-corner flag to "1", to call the 3-State RS.
          itri_corner = 1

! "n1" goes from (ntp_zones, 0) to (ntp_zones + 1, 0).
! "n2" goes from (ntp_zones, 0) to (ntp_zones, 1).

          ixx01 = ntp_zones ;      iyy01 = 0
          ixx1  = ntp_zones + 1 ;  iyy1  = 0        ! Defines gradient for "n1".

          ixx02 = ntp_zones ;      iyy02 = 0
          ixx2  = ntp_zones ;      iyy2  = 1    ! Defines gradient for "n2".

! --------------------------------------------------
! Right-Up Tri-corner.

          ELSE IF ( ( itri_patch_RU .EQ. 1) .AND. ( ixx .EQ. ntp_zones)
     1         .AND. ( iyy .EQ. ntp_zones) ) THEN
! --------------------------------------------------
! "LU --> U1; LD --> U2; RD --> U3" to be used in the 2D Riemann solver
! with 3 states.

          WRITE(6,*)"itri_patch_RU .EQ. 1"

!          U_1 = U_LU
!          U_2 = U_LD
!          U_3 = U_RD

! ------------------------------

! set the tri-corner flag to "1", to call the 3-State RS.
          itri_corner = 1

! "n1" goes from (ntp_zones, ntp_zones) to (ntp_zones, ntp_zones + 1).
! "n2" goes from (ntp_zones, ntp_zones) to (ntp_zones - 1, ntp_zones).

          ixx01 = ntp_zones ;      iyy01 = ntp_zones
          ixx1  = ntp_zones ;      iyy1  = ntp_zones + 1  ! Defines gradient for "n1".

          ixx02 = ntp_zones ;      iyy02 = ntp_zones
          ixx2  = ntp_zones - 1 ;  iyy2  = ntp_zones      ! Defines gradient for "n2".

! --------------------------------------------------
! Left-Up Tri-corner.

          ELSE IF ( ( itri_patch_LU .EQ. 1) .AND. ( ixx .EQ. 0)
     1         .AND. ( iyy .EQ. ntp_zones) ) THEN
! --------------------------------------------------
! "RU --> U1; LD --> U2; RD --> U3" to be used in the 2D Riemann solver
! with 3 states.

          WRITE(6,*)"itri_patch_LU .EQ. 1"

!          U_1 = U_RU
!          U_2 = U_LD
!          U_3 = U_RD

! ------------------------------

! set the tri-corner flag to "1", to call the 3-State RS.
          itri_corner = 1

! We want all gradients for calculating unit vectors to be taken from the interior.
! "n1" goes from (0, ntp_zones) to (1, ntp_zones).
! "n2" goes from (0, ntp_zones) to (0, ntp_zones + 1).

          ixx01 = 0 ;    iyy01 = ntp_zones
          ixx1  = 1 ;    iyy1  = ntp_zones        ! Defines gradient for "n1".

          ixx02 = 0 ;    iyy02 = ntp_zones
          ixx2  = 0 ;    iyy2  = ntp_zones + 1    ! Defines gradient for "n2".

! ------------------------------

          ELSE

! ------------------------------
! Normal Quad-corner; there are 4 input states for the 2D Riemann solver.
! The pre-computed "U_RU, U_LU, U_LD, U_RD" will be used in the 2D Riemann solver
! with 4 states.

          WRITE(6,*)"Normal Quad-corner"

! ------------------------------

! set the tri-corner flag to "0", to call the 4-State RS.
          itri_corner = 0

! We want all gradients for calculating unit vectors to be centered around the corner.
! "n1" goes from (ixx - 1, iyy) to (ixx + 1, iyy).
! "n2" goes from (ixx, iyy - 1) to (ixx, iyy + 1).

          ixx01 = ixx - 1 ;  iyy01 = iyy ;  ixx1 = ixx + 1 ;  iyy1 = iyy  ! Defines gradient for "n1".
          ixx02 = ixx ;  iyy02 = iyy - 1 ;  ixx2 = ixx ;  iyy2 = iyy + 1  ! Defines gradient for "n2".

! ------------------------------

          END IF  ! End of "IF ( ( itri_patch_LD .EQ. 1) .AND.  ... ".

! --------------------------------------------------
! Build local normals so that globally positive z-direction is indeed the
! locally positive z-direction. Store these in "n1, n2". "n3" not needed.
! "n3" is a unit vector in the radial direction and it is easily built.
! "n1" is a unit vector in the local x-direction of the mesh.
! "n2" is a unit vector in the local y-direction of the mesh.
! (The mesh can be skewed in the xy-tangent space. x- and y-directions
!  need not be orthogonal.)

          tempaa = 1.0 / SQRT ( xec**2 + yec**2 + zec**2)

          n3 ( 1) = xec * tempaa
          n3 ( 2) = yec * tempaa
          n3 ( 3) = zec * tempaa

! Obtain "n1" in the xy-tangent space of the sphere in the x-coordinate
! direction.

          n1 ( 1) = vert_coords_patch ( ixx1,  iyy1,  izz, 1)
     1            - vert_coords_patch ( ixx01, iyy01, izz, 1)
          n1 ( 2) = vert_coords_patch ( ixx1,  iyy1,  izz, 2)
     1            - vert_coords_patch ( ixx01, iyy01, izz, 2)
          n1 ( 3) = vert_coords_patch ( ixx1,  iyy1,  izz, 3)
     1            - vert_coords_patch ( ixx01, iyy01, izz, 3)

          tempaa = n1 ( 1) * n3 ( 1) + n1 ( 2) * n3 ( 2)
     1           + n1 ( 3) * n3 ( 3)

          n1 = n1 - tempaa * n3

          tempaa = 1.0 / SQRT ( n1 ( 1)**2 + n1 ( 2)**2 + n1 ( 3)**2)

          n1 ( :) = n1 ( :) * tempaa

! Obtain "n2" in the xy-tangent space of the sphere in the y-coordinate
! direction.

          n2 ( 1) = vert_coords_patch ( ixx2,  iyy2,  izz, 1)
     1            - vert_coords_patch ( ixx02, iyy02, izz, 1)
          n2 ( 2) = vert_coords_patch ( ixx2,  iyy2,  izz, 2)
     1            - vert_coords_patch ( ixx02, iyy02, izz, 2)
          n2 ( 3) = vert_coords_patch ( ixx2,  iyy2,  izz, 3)
     1            - vert_coords_patch ( ixx02, iyy02, izz, 3)

          tempaa = n2 ( 1) * n3 ( 1) + n2 ( 2) * n3 ( 2)
     1           + n2 ( 3) * n3 ( 3)

          n2 = n2 - tempaa * n3

          tempaa = 1.0 / SQRT ( n2 ( 1)**2 + n2 ( 2)**2 + n2 ( 3)**2)

          n2 ( :) = n2 ( :) * tempaa
! --------------------------------------------------
! --------------------------------------------------

          n1_mag = SQRT ( n1 ( 1)**2 + n1 ( 2)**2 + n1 ( 3)**2)
          n2_mag = SQRT ( n2 ( 1)**2 + n2 ( 2)**2 + n2 ( 3)**2)

          angle = ACOS ( ( n1 ( 1) * n2 ( 1)
     1                   + n1 ( 2) * n2 ( 2)
     1                   + n1 ( 3) * n2 ( 3)) / ( n1_mag*n2_mag) )

          WRITE ( 6, *)" "
          WRITE ( 6, *)"n1 = ", n1
          WRITE ( 6, *)"n2 = ", n2
          WRITE ( 6, *)"n1.n2 = ",n1(1)*n2(1) + n1(2)*n2(2) +n1(3)*n2(3)
          WRITE ( 6, *)"n1_mag = ", n1_mag
          WRITE ( 6, *)"n2_mag = ", n2_mag
          WRITE ( 6, *)" "
          WRITE ( 6, *)"Angle = ", angle*180.0/pi
          WRITE ( 6, *)" "

      STOP

! ------------------------------------------------------------
#endif
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! We have built the two arrays "vert_coords_patch, imask_active_patch" for
! the specific patch given by "ix_patch, iy_patch, iz_patch". From
! these two arrays, all the other geometric information is built on the GPU.
! The geometry-building, and stencil-building, routine that are exercised below
! have been highly vectorized and parallelized for the GPU.

! In addition, move the polynomial modes to the GPU.

! ----------------------------------------------------------------------
!$acc enter data copyin(                                                &
!$acc&           vert_coords_patch, imask_active_patch,                 &

!$acc&    ii_mode, jj_mode, kk_mode)

! ----------------------------------------------------------------------

      DO istep = 1, ntstep

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      CALL CPU_TIME ( timing_info1)
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

! ----------------------------------------------------------------------

! Start with "vert_coords_patch" and build all the area-related quadrature points,
! weights and normals and areas for all the areas that bound each zone. Also build
! a length "v_scale_length" by which each zone can be scaled.

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

! ----------------------------------------------------------------------

! Start with "vert_coords_patch" and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

! ----------------------------------------------------------------------

! Now build the moments for that patch. These are useful for building the
! Taylor series basis functions into which the reconstruction problem is projected.

      CALL INIT_BASIS_FUNCTION_COEFS

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         imask_active_patch,

     1         v_quadpts_patch, v_quadwt_patch, v_volume_patch,
     1         v_centroid_patch, v_scale_length,

     1         Coef_Basis_Moments)

! ----------------------------------------------------------------------
      CALL CPU_TIME ( timing_info2)
! ----------------------------------------------------------------------

!$acc kernels
      time_vec ( istep) = timing_info2 - timing_info1
!$acc end kernels

! ----------------------------------------------------------------------

      END DO ! End of "istep = 1, ntstep"

! ----------------------------------------------------------------------

!$acc update host (                                                     &
!$acc&           vert_coords_patch, imask_active_patch,                 &
!$acc&                                                                  &
!$acc&           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch, &
!$acc&           xf_area_patch, xf_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch, &
!$acc&           yf_area_patch, yf_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch, &
!$acc&           zf_area_patch, zf_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           v_scale_length,                                        &
!$acc&                                                                  &
!$acc&           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,  &
!$acc&           v_volume_patch, v_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           Coef_Basis_Moments, istencil_patch, time_vec)

! ----------------------------------------------------------------------

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed the TEST_CUBE_SPHERE_Mesh : "
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "

      DO istep = 1, ntstep

      WRITE ( 6, *)"istep, time : ", istep, time_vec ( istep)
      WRITE ( 6, *)"--------------------------------------------------"

      END DO ! End of "istep = 1, ntstep"

! ----------------------------------------------------------------------
! We created a lot of temporary arrays that are used for processing a timestep
! on the GPU. Since we are at the end of the timestep, we can get rid of
! these arrays on the GPU.

!$acc end data

! ----------------------------------------------------------------------
! Diagnostics from facial geometry construction.
! ----------------------------------------------------------------------




! ----------------------------------------------------------------------
! Diagnostics from volumetric geometry construction.
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Diagnostics about the moments of each zone.
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Diagnostics about the stencil operations.
! ----------------------------------------------------------------------
#if( DEBUG_GEOM == 1)

      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1   "itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU = "
      WRITE ( 6, *)
     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU
      WRITE ( 6, *)" "

! ---------------------------------------------
! Print, Down - along x-axis, iyy = 1.

      iyy = 1

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Down - along x-axis, iyy = 1"

      DO ixx = 1 - NGhost, ntp_zones + NGhost

      WRITE ( 6, *)"istencil_patch( ixx, iyy, 1:2, 4)",
     1              istencil_patch( ixx, iyy, 1, 4),
     1              istencil_patch( ixx, iyy, 2, 4)

      END DO

! ---------------------------------------------
! Print, Up - along x-axis, iyy = ntp_zones.

      iyy = ntp_zones

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Up - along x-axis, iyy = ntp_zones"

      DO ixx = 1 - NGhost, ntp_zones + NGhost

      WRITE ( 6, *)"istencil_patch( ixx, iyy, 1:2, 3) = ",
     1              istencil_patch( ixx, iyy, 1, 3),
     1              istencil_patch( ixx, iyy, 2, 3)

      END DO

! ---------------------------------------------
! Print, Left - along y-axis, ixx = 1.

      ixx = 1

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Left - along y-axis, ixx = 1"

      DO iyy = 1 - NGhost, ntp_zones + NGhost

      WRITE ( 6, *)"istencil_patch( ixx, iyy, 1:2, 2)",
     1              istencil_patch( ixx, iyy, 1, 2),
     1              istencil_patch( ixx, iyy, 2, 2)

      END DO

! ---------------------------------------------
! Print, Right - along y-axis, ixx = ntp_zones.

      ixx = ntp_zones

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Right - along y-axis, ixx = ntp_zones"

      DO iyy = 1 - NGhost, ntp_zones + NGhost

      WRITE ( 6, *)"istencil_patch( ixx, iyy, 1:2, 1)",
     1              istencil_patch( ixx, iyy, 1, 1),
     1              istencil_patch( ixx, iyy, 2, 1)

      END DO

      WRITE ( 6, *)" "

#endif
! ----------------------------------------------------------------------

      END PROGRAM TEST_CUBE_SPHERE_Mesh
