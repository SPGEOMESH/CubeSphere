! PROGRAM BUILD_2D_SECTOR
! compile with BuildGeomesh_c3_src and its dependencies

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      PROGRAM BUILD_2D_SECTOR

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER isector, ntp_patches, ix_patch, iy_patch, ntp_zones,
     1        NGhost, i, j

      PARAMETER ( NGhost = 4, ntp_zones = 64)

      REAL vertices ( 4, 3),
     1     coordinates_gr ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost, 3)

! ----------------------------------------------------------------------

! Treat entire sector as one patch with "ntp_zones" in each of the two directions.
! Vertices go from "0 : ntp_zones". Ghosts are not yet initialized.

      OPEN ( unit = 6, file = "sectors.out")

      ntp_patches = 1  ! Treat whole sector as one patch.
      ix_patch = 1  ;  iy_patch = 1  ! Treat whole sector as one patch.

      isector = 1  ! Choose a sector number, 1 to 6.

      WRITE ( 6, *)" "
      WRITE ( 6, *)"isector = ", isector
      WRITE ( 6, *)" "

      CALL Mesh_Sector_General ( isector, ntp_patches,
     1                         ix_patch, iy_patch, ntp_zones, NGhost,
     1                         vertices, coordinates_gr)

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Vertex 1 has coordinates : ", vertices ( 1, :)
      WRITE ( 6, *)"Vertex 2 has coordinates : ", vertices ( 2, :)
      WRITE ( 6, *)"Vertex 3 has coordinates : ", vertices ( 3, :)
      WRITE ( 6, *)"Vertex 4 has coordinates : ", vertices ( 4, :)
      WRITE ( 6, *)" "

      DO j = 0, ntp_zones

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Writing row j = ", j
      WRITE ( 6, *)" "

        DO i = 0, ntp_zones

        WRITE ( 6, *)"Coordinates for vertex i, j = ", i, j
        WRITE ( 6, *) coordinates_gr ( i, j, :)

        END DO  ! End of "DO i = 0, ntp_zones"

      END DO  ! End of "DO j = 0, ntp_zones"

      WRITE ( 6, *)" "
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------

      STOP

      END PROGRAM BUILD_2D_SECTOR

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
