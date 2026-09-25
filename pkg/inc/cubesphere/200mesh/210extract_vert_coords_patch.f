





      SUBROUTINE EXTRACT_VERT_COORDS_PATCH

     1       ( ntp_global, nr_global, NGhost_global,
     1         ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         isector, ix_patch, iy_patch, iz_patch,
     1         ilevel, ilevel_max,

     1         coordinates_gr_global, rb_save_global,

     1         vert_coords_patch, imask_active_patch)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara and Sethu Subramanyan

! "ntp_global, nr_global, NGhost_global" are the number of zones on each sector
! in the theta-phi, radial and ghost directions.

! "ntp_zones, nr_zones, NGhost" are the number of zones in the theta-phi direction,
! the radial direction and the ghost zones on each patch. Any patch in the entire
! AMR hierarchy has this fixed number of zones and ghost zones.
! This is baked into the data structures for the whole code.

! "ntp_patches, nr_patches" are the total number of patches in the theta-phi and
! radial directions IF that level were indeed completely full. Other than the
! base level, the AMR levels do not have to be completely full.

! "ix_patch, iy_patch, iz_patch" is the patch at level "ilevel" that we indeed want.
! Again, this is based on expecting a full level whether there are physical levels
! there or not.

! Maximum levels assuming a refinement ratio of 2.

! "coordinates_gr_global" has ghost zones filled in and has the coordinates on the
! unit sphere. I.e., it is 2D in theta-phi space. It is the global mesh (in the
! theta-phi direction) that we would have obtained on the finest level, if it were
! fully filled in. The finest level is "ilevel_max".

! "rb_save_global" has ghost zones filled in and has the radial coordinates.
! It is the global mesh (in radial direction) that we would have obtained
! on the finest level, if it were fully filled in. The finest level is "ilevel_max".

! "vert_coords_patch" is the local patch whose vertices we fill in.
! This is on level "ilevel", so we will have to skip count on the global mesh. todo skip count ?

! "imask_active_patch" is 1 for active zones, 0 for inactive ones
! that may lie close to the corners. todo that may lie close to the corners ?

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER ntp_global, nr_global, NGhost_global,
     1        ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,
     1        isector, ix_patch, iy_patch, iz_patch,
     1        ilevel, ilevel_max

      REAL coordinates_gr_global
     1     ( 0 - NGhost_global: ntp_global + NGhost_global,
     1       0 - NGhost_global: ntp_global + NGhost_global, 3, 6)

      REAL rb_save_global ( 0 - NGhost_global:
     1                      nr_global + NGhost_global)

! Output variable, extracted vertex coordinates for this patch.
      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

! Output variable, mask which is zone-centered and 1 if zone is active, 0 if inactive.
      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

! ----------------------------------------------------------------------
! Variables below this line are internal to the subroutine.

      INTEGER ixx, iyy, izz, ixx_you, iyy_you, izz_you, ixxm1, iyym1,
     1        iscale

      REAL tempaa, tempab, tempac, tempad

      INTEGER imask_1_layer ( 1 - NGhost: ntp_zones + NGhost,
     1                        1 - NGhost: ntp_zones + NGhost)

      REAL vert_coords_1_layer ( 0 - NGhost: ntp_zones + NGhost,
     1                           0 - NGhost: ntp_zones + NGhost, 3)

! ----------------------------------------------------------------------
! STEP 1) Build "iscale" and check for consistencies.

      iscale = 2**(ilevel_max - ilevel)

! ----------
! "ntp_zones * ntp_patches" gives number of zones in theta-phi direction on
! level "ilevel". Therefore, "ntp_zones * ntp_patches * iscale" should give
! number of zones on the finest level in theta-phi direction.

      IF ( ntp_zones * ntp_patches * iscale .NE. ntp_global) THEN

      WRITE ( 6, *)"Stopping in EXTRACT_VERT_COORDS_PATCH because : "

      WRITE ( 6, *)"ntp_zones * ntp_patches * iscale .NE. ntp_global "
      WRITE ( 6, *)" "
      WRITE ( 6, *)"ntp_zones, ntp_patches, ntp_global = ",
     1              ntp_zones, ntp_patches, ntp_global
      WRITE ( 6, *)"iscale, ilevel, ilevel_max = ",
     1              iscale, ilevel, ilevel_max

      STOP
      END IF

! ----------
! "nr_zones * nr_patches" gives number of zones in radial direction on
! level "ilevel". Therefore, "nr_zones * nr_patches * iscale" should give
! number of zones on the finest level in radial direction.

      IF ( nr_zones * nr_patches * iscale .NE. nr_global) THEN

      WRITE ( 6, *)"Stopping in EXTRACT_VERT_COORDS_PATCH because : "

      WRITE ( 6, *)"nr_zones * nr_patches * iscale .NE. nr_global "
      WRITE ( 6, *)" "
      WRITE ( 6, *)"nr_zones, nr_patches, nr_global = ",
     1              nr_zones, nr_patches, nr_global
      WRITE ( 6, *)"iscale, ilevel, ilevel_max = ",
     1              iscale, ilevel, ilevel_max

      STOP
      END IF

! ----------
! Because the global mesh at finest level has to be conformant, the number
! of ghost zones should also scale.

      IF ( NGhost * iscale .NE. NGhost_global) THEN

      WRITE ( 6, *)"Stopping in EXTRACT_VERT_COORDS_PATCH because : "

      WRITE ( 6, *)"NGhost * iscale .NE. NGhost_global"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"NGhost, NGhost_global = ", NGhost, NGhost_global
      WRITE ( 6, *)"iscale, ilevel, ilevel_max = ",
     1              iscale, ilevel, ilevel_max

      STOP
      END IF

! ----------
! "ix_patch, iy_patch, iz_patch" should be within range; "NGhost"
! should not be too large large.

      IF ( ( ix_patch .LT. 1) .OR. ( ix_patch .GT. ntp_patches)) THEN

      WRITE ( 6, *)"Stopping in EXTRACT_VERT_COORDS_PATCH because : "

      WRITE ( 6, *)"ix_patch out of bounds ", ix_patch, 1, ntp_patches

      STOP
      END IF

! ----------

      IF ( ( iy_patch .LT. 1) .OR. ( iy_patch .GT. ntp_patches)) THEN

      WRITE ( 6, *)"Stopping in EXTRACT_VERT_COORDS_PATCH because : "

      WRITE ( 6, *)"iy_patch out of bounds ", iy_patch, 1, ntp_patches

      STOP
      END IF

! ----------

      IF ( ( iz_patch .LT. 1) .OR. ( iz_patch .GT.  nr_patches)) THEN

      WRITE ( 6, *)"Stopping in EXTRACT_VERT_COORDS_PATCH because : "

      WRITE ( 6, *)"iz_patch out of bounds ", iz_patch, 1, nr_patches

      STOP
      END IF

! ----------

      IF ( ( NGhost .GT. ntp_zones) .OR. ( NGhost .GT. nr_zones)) THEN

      WRITE ( 6, *)"Stopping in EXTRACT_VERT_COORDS_PATCH because : "

      WRITE ( 6, *)"NGhost, ntp_zones, nr_zones = ",
     1              NGhost, ntp_zones, nr_zones

      STOP
      END IF

! ----------------------------------------------------------------------
! STEP 2) Fill in one layer of "vert_coords_1_layer" using
! "coordinates_gr_global". This is still only one layer and it still lives
! on the unit sphere.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, ixx_you, iyy_you)                 &

!$OMP&            SHARED ( ntp_zones, NGhost, ix_patch, iy_patch,       &
!$OMP&    iscale, isector,                                              &
!$OMP&    vert_coords_1_layer, coordinates_gr_global)                   &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

      DO iyy = 0 - NGhost, ntp_zones + NGhost

      iyy_you = ( ntp_zones * ( iy_patch - 1) + iyy ) * iscale

        DO ixx = 0 - NGhost, ntp_zones + NGhost

        ixx_you = ( ntp_zones * ( ix_patch - 1) + ixx ) * iscale

        vert_coords_1_layer ( ixx, iyy, :) =
     1  coordinates_gr_global ( ixx_you, iyy_you, :, isector)

        END DO  ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost" loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! STEP 3) Examine "vert_coords_1_layer" zones for bad values and identify
! bad zones. Store that information in "imask_1_layer".
! Once this is done, we store that info. in "imask_active_patch".

      imask_1_layer = 1

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, ixxm1, iyym1,                     &
!$OMP&    tempaa, tempab, tempac, tempad)                               &

!$OMP&            SHARED ( ntp_zones, NGhost, smallnum,                 &
!$OMP&    vert_coords_1_layer, imask_1_layer)                           &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)


      DO iyy = 1 - NGhost, ntp_zones + NGhost
      iyym1 = iyy - 1

        DO ixx = 1 - NGhost, ntp_zones + NGhost
        ixxm1 = ixx - 1

        tempaa = SQRT ( vert_coords_1_layer ( ixx, iyy, 1)**2
     1                + vert_coords_1_layer ( ixx, iyy, 2)**2
     1                + vert_coords_1_layer ( ixx, iyy, 3)**2 )

        tempab = SQRT ( vert_coords_1_layer ( ixxm1, iyy, 1)**2
     1                + vert_coords_1_layer ( ixxm1, iyy, 2)**2
     1                + vert_coords_1_layer ( ixxm1, iyy, 3)**2 )

        tempac = SQRT ( vert_coords_1_layer ( ixx, iyym1, 1)**2
     1                + vert_coords_1_layer ( ixx, iyym1, 2)**2
     1                + vert_coords_1_layer ( ixx, iyym1, 3)**2 )

        tempad = SQRT ( vert_coords_1_layer ( ixxm1, iyym1, 1)**2
     1                + vert_coords_1_layer ( ixxm1, iyym1, 2)**2
     1                + vert_coords_1_layer ( ixxm1, iyym1, 3)**2 )

        IF ( ( tempaa .LT. smallnum) .OR. ( tempab .LT. smallnum) .OR.
     1       ( tempac .LT. smallnum) .OR. ( tempad .LT. smallnum)) THEN
           imask_1_layer ( ixx, iyy) = 0
!          WRITE(6, *) "imask_1_layer(ixx, iyy) = 0; ixx iyy", ixx, iyy
        END IF
        END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

!$OMP END PARALLEL DO

      DO izz = 1 - NGhost, nr_zones + NGhost
      imask_active_patch ( :, :, izz) = imask_1_layer ( :, :)
      END DO

! ----------------------------------------------------------------------
! STEP 4) Using "vert_coords_1_layer" and "rb_save_global" we now obtain
! "vert_coords_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izz_you)                     &

!$OMP&            SHARED ( ntp_zones, nr_zones, NGhost, iz_patch,       &
!$OMP&    iscale,                                                       &
!$OMP&    rb_save_global, vert_coords_1_layer, vert_coords_patch)       &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

      DO izz = 0 - NGhost, nr_zones + NGhost

        DO iyy = 0 - NGhost, ntp_zones + NGhost

        izz_you = ( nr_zones * ( iz_patch - 1) + izz ) * iscale

          DO ixx = 0 - NGhost, ntp_zones + NGhost

          vert_coords_patch ( ixx, iyy, izz, :) =
     1    rb_save_global ( izz_you) * vert_coords_1_layer ( ixx, iyy, :)

          END DO  ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 0 - NGhost, nr_zones + NGhost" loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if ( DEBUG_GEOM == 1)

      DO izz = 0 - NGhost, nr_zones + NGhost

      WRITE ( 6, *)"--------------------"

        DO iyy = 0 - NGhost, ntp_zones + NGhost

        WRITE ( 6, *)"----------"

          DO ixx = 0 - NGhost, ntp_zones + NGhost

          ixx_you = ( ntp_zones * ( ix_patch - 1) + ixx ) * iscale
          iyy_you = ( ntp_zones * ( iy_patch - 1) + iyy ) * iscale
          izz_you = ( nr_zones * ( iz_patch - 1) + izz ) * iscale

          WRITE ( 6, *)"Local ixx, iyy, izz = ", ixx, iyy, izz
          WRITE ( 6, *)"Global finest ixx_you, iyy_you, izz_you = ",
     1                  ixx_you, iyy_you, izz_you
          WRITE ( 6, *)"Coordinates : ",
     1                  vert_coords_patch ( ixx, iyy, izz, :)
          WRITE ( 6, *)" "

          END DO  ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 0 - NGhost, nr_zones + NGhost" loop.

#endif

! ----------------------------------------------------------------------

      END SUBROUTINE EXTRACT_VERT_COORDS_PATCH







