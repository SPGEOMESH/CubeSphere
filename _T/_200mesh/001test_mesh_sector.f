




#include <cubesphere/200mesh/201mesh_sector_general.f>



      PROGRAM TEST_Mesh_Sector

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Tests "Mesh_Sector" for its ability to set up the mesh for a single sector.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER isector, ntp_patches, ix_patch, iy_patch, ntp_zones,
     1        NGhost,
     1        i, j, im1, jm1, ii, jj

      PARAMETER ( ntp_zones = 6, NGhost = 0)

      REAL vertices ( 4, 3),
     1     coordinates_gr ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost, 3)

      REAL coplanarity, global_coplanarity, t,
     1     v1 (3), v2 (3), v12 (3), v12_prev (3)

! ----------------------------------------------------------------------

      global_coplanarity = 0.0

!     OPEN ( UNIT = 6, FILE = "test_mesh_sector.out")

! "isector = 1" is the easiest to use for debug purposes because x- and y-
! directions line up with the x- and y-coordinate system. Since all sectors
! use the exact same logic, verifying it for Sector 1 is a good enough test.

      isector = 1

      ntp_patches = 2

      WRITE ( 6, *)"isector, ntp_patches = ", isector, ntp_patches
      WRITE ( 6, *)"ntp_zones = ", ntp_zones
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------


      DO iy_patch = 1, ntp_patches

      DO ix_patch = 1, ntp_patches

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Starting Work on :- "
      WRITE ( 6, *)" "
      WRITE ( 6, *)"ix_patch, iy_patch = ", ix_patch, iy_patch
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

      CALL Mesh_Sector_General ( isector, ntp_patches,
     1                         ix_patch, iy_patch, ntp_zones, NGhost,
     1                         vertices, coordinates_gr)

! ----------------------------------------

      DO i = 1, 4

      WRITE ( 6, *)" "
      WRITE ( 6, *)"vertex i = ", i
      WRITE ( 6, *) vertices ( i, :)

      END DO

      WRITE ( 6, *)" "
      WRITE ( 6, *)"----------------------------------------"
      WRITE ( 6, *)" "

! ----------------------------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Writing out coordinates_gr for this patch"
      WRITE ( 6, *)" "

      DO j = 0, ntp_zones

      WRITE ( 6, *)"--------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"local y-index, j = ", j
      WRITE ( 6, *)" "

! ----------------------------------------

        DO i = 0, ntp_zones

        WRITE ( 6, *)"local xy-indices, i, j = ", i, j
        WRITE ( 6, *)"coordinates_gr = ", coordinates_gr ( i, j, :)
        WRITE ( 6, *)" "

        END DO  ! End of "DO i = 0, ntp_zones" loop.


        WRITE ( 6, *)"Checking for coplanarity in y-index j = ", j
        coplanarity = 0.0

        DO i = 1, ntp_zones
        im1 = i - 1

        v1 ( :) = coordinates_gr ( im1, j, :)
        v2 ( :) = coordinates_gr ( i, j, :)

        v12 (1) = v1 (2) * v2 (3) - v1 (3) * v2 (2)
        v12 (2) = v1 (3) * v2 (1) - v1 (1) * v2 (3)
        v12 (3) = v1 (1) * v2 (2) - v1 (2) * v2 (1)

        t = 1.0 / SQRT ( v12 (1)**2 + v12 (2)**2 + v12 (3)**2 )
        v12 = v12 * t

          IF ( i .GT. 1) coplanarity = AMAX1 ( coplanarity,
     1    SQRT ( (v12 (1)-v12_prev (1))**2 + (v12 (2)-v12_prev (2))**2
     1         + (v12 (3)-v12_prev (3))**2 ) )

        v12_prev = v12

        END DO  ! End of "DO i = 1, ntp_zones" loop


        global_coplanarity = AMAX1 ( global_coplanarity, coplanarity)

        WRITE ( 6, *)"Checked coplanarity factor = ", coplanarity
        WRITE ( 6, *)" "
        WRITE ( 6, *)"--------------------"

! ----------------------------------------

      END DO  ! End of "DO j = 0, ntp_zones" loop.

! ----------------------------------------

      DO i = 0, ntp_zones

      WRITE ( 6, *)"--------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Flipping dominant direction: local x-index, i = ", i
      WRITE ( 6, *)" "

! ----------------------------------------

        WRITE ( 6, *)"Checking for coplanarity in x-index, i = ", i
        coplanarity = 0.0

        DO j = 1, ntp_zones
        jm1 = j - 1

        v1 ( :) = coordinates_gr ( i, jm1, :)
        v2 ( :) = coordinates_gr ( i, j, :)

        v12 (1) = v1 (2) * v2 (3) - v1 (3) * v2 (2)
        v12 (2) = v1 (3) * v2 (1) - v1 (1) * v2 (3)
        v12 (3) = v1 (1) * v2 (2) - v1 (2) * v2 (1)

        t = 1.0 / SQRT ( v12 (1)**2 + v12 (2)**2 + v12 (3)**2 )
        v12 = v12 * t

          IF ( j .GT. 1) coplanarity = AMAX1 ( coplanarity,
     1    SQRT ( (v12 (1)-v12_prev (1))**2 + (v12 (2)-v12_prev (2))**2
     1         + (v12 (3)-v12_prev (3))**2 ) )

        v12_prev = v12

        END DO  ! End of "DO j = 1, ntp_zones" loop.


        global_coplanarity = AMAX1 ( global_coplanarity, coplanarity)

        WRITE ( 6, *)"Checked coplanarity factor = ", coplanarity
        WRITE ( 6, *)" "
        WRITE ( 6, *)"--------------------"

! ----------------------------------------

      END DO  ! End of "DO i = 0, ntp_zones" loop.

! ----------------------------------------

      WRITE ( 6, *)"----------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Finished Working on :- "
      WRITE ( 6, *)" "
      WRITE ( 6, *)"ix_patch, iy_patch = ", ix_patch, iy_patch
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

      END DO  ! End of "DO ix_patch = 1, ntp_patches" loop.

      END DO  ! End of "DO iy_patch = 1, ntp_patches" loop.

      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"global_coplanarity = ", global_coplanarity
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Getting almost zero on this number is a good thing."
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------
! We may also want to compare with a patched mesh with a global mesh on a single sector.
! This will be very slow because of multiple calls

      WRITE ( 6, *)"----------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Compare with Global Sectorial coverage with ",
     1             ntp_patches * ntp_zones, " zones in each direction. "

      DO jj = 0, ntp_patches * ntp_zones

      j = MOD ( jj, ntp_zones)
      iy_patch = ( jj / ntp_zones) + 1

      IF ( jj .EQ. ntp_patches * ntp_zones) THEN
      j = ntp_zones
      iy_patch = iy_patch - 1
      END IF

      WRITE ( 6, *)"--------------------"
      WRITE ( 6, *)"local iy_patch, j = ", iy_patch, j
      WRITE ( 6, *)"global y-index, jj = ", jj
      WRITE ( 6, *)" "

        DO ii = 0, ntp_patches * ntp_zones

        i = MOD ( ii, ntp_zones)
        ix_patch = ( ii / ntp_zones) + 1

        IF ( ii .EQ. ntp_patches * ntp_zones) THEN
        i = ntp_zones
        ix_patch = ix_patch - 1
        END IF

! ----------------------------------------

        CALL Mesh_Sector_General ( isector, ntp_patches,
     1                         ix_patch, iy_patch, ntp_zones, NGhost,
     1                         vertices, coordinates_gr)

! ----------------------------------------

        WRITE ( 6, *)"local ix_patch, i = ", ix_patch, i
        WRITE ( 6, *)"global xy-indices, ii, jj = ", ii, jj
        WRITE ( 6, *)"coordinates_gr = ", coordinates_gr ( i, j, :)
        WRITE ( 6, *)" "

! ----------------------------------------

        END DO  ! End of "DO i = 0, ntp_patches * ntp_zones" loop.

      END DO  ! End of "DO j = 0, ntp_patches * ntp_zones" loop.

! ----------------------------------------------------------------------

      END PROGRAM TEST_Mesh_Sector
