






      PROGRAM TEST_FACIAL_VOL_QUADPTS

! ----------------------------------------------------------------------

! This program tests for a single zone. "vert_coords_patch" for a spherical
! zone are handed in. It shows the locations of the quadrature points and weights
! at the faces of that zone.

! It also cross-compares a heavy-weight set of volumetric quadrature points
! (from COORDS_TO_WT_QUADPTS_VOLUME_ISO2) with NUMORDER^3 quadrature points with
! a light-weight set of "NQUICK_VOL_QUAD_PTS = 6" volumetric quadrature points
! (from QUICK_QUADPTS_VOLUME_ISO2).

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nr_zones, ntp_zones, NGhost, ixx, iyy, izz, j

      REAL radius, tempaa

      PARAMETER ( nr_zones = 1, ntp_zones = 1, NGhost = 0)

! Input variable, gives us the vertices of a cube sphere mesh.
      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

! Input variable, mask which is zone-centered and 1 if zone is active, 0 if inactive.
      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

! ----------------------------------------------------------------------
! Output variable, gives us the 3-component position vector for each of the
! "NQUICK_VOL_QUAD_PTS" number of quick quadrature points in the volume.
      REAL vq_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NQUICK_VOL_QUAD_PTS)

! Output variable, gives us the scalar weight  for each of the
! "NQUICK_VOL_QUAD_PTS" number of quadrature points in the volume.
      REAL vq_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NQUICK_VOL_QUAD_PTS)

! Local, internal copy of "v_volume_patch". Can be compared to the one
! obtained from tensor product quadrature points to gauge the accuracy of our
! choice of quick quadrature points.
      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ----------------------------------------------------------------------
! Output variable, gives us the 3-component position vector for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the x-face of a zone.
      REAL xb_quadpts_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

! Output variable, gives us the scalar weight  for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the x-face of a zone.
      REAL xb_quadwt_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

! Output variable, gives us the unit normal in the x-face. Only one is needed
! because this face is flat.
      REAL xb_quadnorml_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------

! Output variable, gives us the 3-component position vector for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the y-face of a zone.
      REAL yb_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

! Output variable, gives us the scalar weight  for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the y-face of a zone.
      REAL yb_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

! Output variable, gives us the unit normal in the y-face. Only one is needed
! because this face is flat.
      REAL yb_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------

! Output variable, gives us the 3-component position vector for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the z-face of a zone.
      REAL zb_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

! Output variable, gives us the scalar weight  for each of the
! "NBFLD_QUAD_PTS" number of quadrature points in the z-face of a zone.
      REAL zb_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   NBFLD_QUAD_PTS)

! Output variable, gives us the unit normal in the z-face. Because this face is
! curved, each quadrature point will have a unit normal that points in a slightly
! different direction.
      REAL zb_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NBFLD_QUAD_PTS)

! ------------------------------

! Output variable, gives us the 3-component position vector for each of the
! "NUMORDER * NUMORDER * NUMORDER" number of quadrature points in the volume.
      REAL v_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER * NUMORDER)

! Output variable, gives us the scalar weight  for each of the
! "NUMORDER * NUMORDER * NUMORDER" number of quadrature points in the volume.
      REAL v_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER * NUMORDER)

! Output variable, gives us the derivative of the position vector relative to
! reference coordinates for each of the "NUMORDER * NUMORDER * NUMORDER" number of
! quadrature points in the volume. These are needed in the ADER update step.
      REAL v_quadtangent_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, 3, NUMORDER * NUMORDER * NUMORDER)

! Output variable, gives us the centroid of the zone.
      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ----------------------------------------------------------------------

      OPEN ( unit = 6, file = "COORDS_TO_WT.out")

      radius = 10.0

      vert_coords_patch ( 0, 0, 0, 1) = -0.5
      vert_coords_patch ( 0, 0, 0, 2) = -0.5
      vert_coords_patch ( 0, 0, 0, 3) = -0.5 + radius

      vert_coords_patch ( 1, 0, 0, 1) = 0.5
      vert_coords_patch ( 1, 0, 0, 2) = -0.5
      vert_coords_patch ( 1, 0, 0, 3) = -0.5 + radius

      vert_coords_patch ( 0, 1, 0, 1) = -0.5
      vert_coords_patch ( 0, 1, 0, 2) = 0.5
      vert_coords_patch ( 0, 1, 0, 3) = -0.5 + radius

      vert_coords_patch ( 1, 1, 0, 1) = 0.5
      vert_coords_patch ( 1, 1, 0, 2) = 0.5
      vert_coords_patch ( 1, 1, 0, 3) = -0.5 + radius


      vert_coords_patch ( 0, 0, 1, 1) = -0.5
      vert_coords_patch ( 0, 0, 1, 2) = -0.5
      vert_coords_patch ( 0, 0, 1, 3) = 0.5 + radius

      vert_coords_patch ( 1, 0, 1, 1) = 0.5
      vert_coords_patch ( 1, 0, 1, 2) = -0.5
      vert_coords_patch ( 1, 0, 1, 3) = 0.5 + radius

      vert_coords_patch ( 0, 1, 1, 1) = -0.5
      vert_coords_patch ( 0, 1, 1, 2) = 0.5
      vert_coords_patch ( 0, 1, 1, 3) = 0.5 + radius

      vert_coords_patch ( 1, 1, 1, 1) = 0.5
      vert_coords_patch ( 1, 1, 1, 2) = 0.5
      vert_coords_patch ( 1, 1, 1, 3) = 0.5 + radius

      imask_active_patch = 1

      CALL COORDS_TO_WT_QUADPTS_BFLD_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

! ----------------------------------------------------------------------
      WRITE ( 6, *)"-------------------------------------------------"

      ixx = 1
      iyy = 1
      izz = 1

      tempaa = 0.0

      DO j = 1, NBFLD_QUAD_PTS
      WRITE ( 6, *)"--------- x-face ", ixx, iyy, izz, j," -----------"
      tempaa = tempaa + xb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"xb_quadwt_patch ( ixx, iyy, izz, j) = ",
     1              xb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"xb_quadpts_patch ( ixx, iyy, izz, :, j) = ",
     1              xb_quadpts_patch ( ixx, iyy, izz, :, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"xb_quadnorml_patch ( ixx, iyy, izz, :) = ",
     1              xb_quadnorml_patch ( ixx, iyy, izz, :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)SQRT ( xb_quadnorml_patch ( ixx, iyy, izz, 1)**2
     1                  + xb_quadnorml_patch ( ixx, iyy, izz, 2)**2
     1                  + xb_quadnorml_patch ( ixx, iyy, izz, 3)**2)
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"sum of weights = ", tempaa

! ----------
      WRITE ( 6, *)"-------------------------------------------------"

      ixx = 0
      iyy = 1
      izz = 1

      tempaa = 0.0

      DO j = 1, NBFLD_QUAD_PTS
      WRITE ( 6, *)"--------- x-face ", ixx, iyy, izz, j," -----------"
      tempaa = tempaa + xb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"xb_quadwt_patch ( ixx, iyy, izz, j) = ",
     1              xb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"xb_quadpts_patch ( ixx, iyy, izz, :, j) = ",
     1              xb_quadpts_patch ( ixx, iyy, izz, :, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"xb_quadnorml_patch ( ixx, iyy, izz, :) = ",
     1              xb_quadnorml_patch ( ixx, iyy, izz, :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)SQRT ( xb_quadnorml_patch ( ixx, iyy, izz, 1)**2
     1                  + xb_quadnorml_patch ( ixx, iyy, izz, 2)**2
     1                  + xb_quadnorml_patch ( ixx, iyy, izz, 3)**2)
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"sum of weights = ", tempaa

! ----------
      WRITE ( 6, *)"-------------------------------------------------"

      ixx = 1
      iyy = 1
      izz = 1

      tempaa = 0.0

      DO j = 1, NBFLD_QUAD_PTS
      WRITE ( 6, *)"--------- y-face ", ixx, iyy, izz, j," -----------"
      tempaa = tempaa + yb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"yb_quadwt_patch ( ixx, iyy, izz, j) = ",
     1              yb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"yb_quadpts_patch ( ixx, iyy, izz, :, j) = ",
     1              yb_quadpts_patch ( ixx, iyy, izz, :, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"yb_quadnorml_patch ( ixx, iyy, izz, :) = ",
     1              yb_quadnorml_patch ( ixx, iyy, izz, :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)SQRT ( yb_quadnorml_patch ( ixx, iyy, izz, 1)**2
     1                  + yb_quadnorml_patch ( ixx, iyy, izz, 2)**2
     1                  + yb_quadnorml_patch ( ixx, iyy, izz, 3)**2)
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"sum of weights = ", tempaa

! ----------
      WRITE ( 6, *)"-------------------------------------------------"

      ixx = 1
      iyy = 0
      izz = 1

      tempaa = 0.0

      DO j = 1, NBFLD_QUAD_PTS
      WRITE ( 6, *)"--------- y-face ", ixx, iyy, izz, j," -----------"
      tempaa = tempaa + yb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"yb_quadwt_patch ( ixx, iyy, izz, j) = ",
     1              yb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"yb_quadpts_patch ( ixx, iyy, izz, :, j) = ",
     1              yb_quadpts_patch ( ixx, iyy, izz, :, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"yb_quadnorml_patch ( ixx, iyy, izz, :) = ",
     1              yb_quadnorml_patch ( ixx, iyy, izz, :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)SQRT ( yb_quadnorml_patch ( ixx, iyy, izz, 1)**2
     1                  + yb_quadnorml_patch ( ixx, iyy, izz, 2)**2
     1                  + yb_quadnorml_patch ( ixx, iyy, izz, 3)**2)
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"sum of weights = ", tempaa

! ----------
      WRITE ( 6, *)"-------------------------------------------------"

      ixx = 1
      iyy = 1
      izz = 1

      tempaa = 0.0

      DO j = 1, NBFLD_QUAD_PTS
      WRITE ( 6, *)"--------- z-face ", ixx, iyy, izz, j," -----------"
      tempaa = tempaa + zb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"zb_quadwt_patch ( ixx, iyy, izz, j) = ",
     1              zb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"zb_quadpts_patch ( ixx, iyy, izz, :, j) = ",
     1              zb_quadpts_patch ( ixx, iyy, izz, :, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"zb_quadnorml_patch ( ixx, iyy, izz, :, j) = ",
     1              zb_quadnorml_patch ( ixx, iyy, izz, :, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)SQRT ( zb_quadnorml_patch ( ixx, iyy, izz, 1,j)**2
     1                  + zb_quadnorml_patch ( ixx, iyy, izz, 2,j)**2
     1                  + zb_quadnorml_patch ( ixx, iyy, izz, 3,j)**2)
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"sum of weights = ", tempaa

! ----------
      WRITE ( 6, *)"-------------------------------------------------"

      ixx = 1
      iyy = 1
      izz = 0

      tempaa = 0.0

      DO j = 1, NBFLD_QUAD_PTS
      WRITE ( 6, *)"--------- z-face ", ixx, iyy, izz, j," -----------"
      tempaa = tempaa + zb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"zb_quadwt_patch ( ixx, iyy, izz, j) = ",
     1              zb_quadwt_patch ( ixx, iyy, izz, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"zb_quadpts_patch ( ixx, iyy, izz, :, j) = ",
     1              zb_quadpts_patch ( ixx, iyy, izz, :, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"zb_quadnorml_patch ( ixx, iyy, izz, :, j) = ",
     1              zb_quadnorml_patch ( ixx, iyy, izz, :, j)
      WRITE ( 6, *)" "
      WRITE ( 6, *)SQRT ( zb_quadnorml_patch ( ixx, iyy, izz, 1,j)**2
     1                  + zb_quadnorml_patch ( ixx, iyy, izz, 2,j)**2
     1                  + zb_quadnorml_patch ( ixx, iyy, izz, 3,j)**2)
      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"sum of weights = ", tempaa

! ----------------------------------------------------------------------

! The real test is that from inside this subroutine, and before the last
! loop nest, we are also evaluating "v_volume_patch ( ixx, iyy, izz)".
! This happens in STEP 2.b). we want to write out that value and compare
! it with the value that we get from "COORDS_TO_WT_QUADPTS_VOLUME_ISO2".

      WRITE ( 6, *)"From QUICK_QUADPTS_VOLUME_ISO2 we evaluate"
      WRITE ( 6, *)"v_volume_patch ( ixx, iyy, izz) after STEP 2.b = "

      CALL QUICK_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           vq_quadpts_patch, vq_quadwt_patch)

! ------------------------------

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

      WRITE ( 6, *)"From COORDS_TO_WT_QUADPTS_VOLUME_ISO2 we evaluate"
      WRITE ( 6, *)"v_volume_patch ( ixx, iyy, izz) = ",
     1              v_volume_patch ( ixx, iyy, izz)


! ----------------------------------------------------------------------

      END PROGRAM TEST_FACIAL_VOL_QUADPTS

