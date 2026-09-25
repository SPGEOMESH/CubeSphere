






! helper, only used by Mesh_sector_general
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE SWIVEL_N ( a1, a2, ntp_zones, arr_in_patch)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Swivels from "a1" to "a2" subdividing that arc into "ntp_zones" number of arclets.
! "p1, p2" contain the end coordinates of arc that is to be subdivided.
! The swiveling action goes from "a1" to "a2".

! "arr_in_patch ( 0: ntp_zones, 3)" contains the subdivision of that arc.

! Vectors "a1, a2" are expected to be unit vectors with an acute angle
! between them.

! It is implicitly assumed that there is an acute angle between "a1" and "a2".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ntp_zones, i

      REAL a1 ( 3), a2 ( 3), arr_in_patch ( 0: ntp_zones, 3),
     1     x ( 3), y ( 3), dot, theta, d_theta, t

! ----------------------------------------------------------------------

      dot = a1 ( 1) * a2 ( 1) + a1 ( 2) * a2 ( 2) + a1 ( 3) * a2 ( 3)
      theta = ACOS ( dot)  ! Angle between "a1" and "a2" in radians.

      d_theta = theta / REAL ( ntp_zones)

! Make sure that the end points are absolutely identical to the input vectors.
! This will help in ensuring that abutting sectors match up exactly.
      arr_in_patch ( 0, :) = a1 ( :)
      arr_in_patch ( ntp_zones, :) = a2 ( :)

! Build two orthogonal vectors "x" and "y" with "x" coinciding with "a1".
      x = a1
      y = a2 - dot * x
      t = SQRT ( y ( 1)**2 + y ( 2)**2 + y ( 3)**2 )
      y = y / t

      DO i = 1, ntp_zones - 1
        t = d_theta * i
        arr_in_patch ( i, :) = COS ( t) * x ( :) + SIN ( t) * y ( :)
      END DO

! ----------------------------------------------------------------------

      END SUBROUTINE SWIVEL_N












! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE Mesh_Sector_General ( isector, ntp_patches,
     1                         ix_patch, iy_patch, ntp_zones, NGhost,
     1                         vertices, coordinates_gr)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Divides cubed sphere into six sectors. "isector = 1, 6" labels each of
! those six sectors. The center is at the origin and the sectors are symmetrically
! located about the origin. The vertices of the sectors coincide with the
! vertices of a cube that is inscribed in the unit sphere.

! Realistically, it is best to use this subroutine for meshing each entire sector.
! In other words, call this subroutine with "isector = 1, ..., 6" and
! "ntp_patches = ix_patch = iy_patch = 1" and "ntp_zones = ntp_global". We then
! load up "coordinates_gr_global" with the global mesh at the finest desired resolution.

! This does not fill the ghost zones in "coordinates_gr_global". That is subsequently
! done in "FILL_GHOST2D_GLOBAL".

! "ntp_patches" is the number of patches along each direction of one of those
! sectors. These are patches in the theta-phi direction, i.e. perpendicular
! to the radial direction.

! "ix_patch, iy_patch" tells us which patch we want in the chosen sector.

! "ntp_zones" gives us the number of zones in that block.

! "vertices ( 4, 3)" : gives us the coordinates (on the unit sphere) for
! the four end points that make up a patch. The coordinates are specified as a
! unit vector in Cartesian geometry.

! "coordinates_gr ( 0: ntp_zones, 0: ntp_zones, 3)" Vertices of the zones on that
! chosen patch. These are coordinates on the unit sphere given in Cartesian geometry.

! All the swiveling action takes place about the origin.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER isector, ntp_patches, ix_patch, iy_patch, ntp_zones,
     1        NGhost

      REAL vertices ( 4, 3),
     1     coordinates_gr ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost, 3)

! ----------------------------------------------------------------------
! Variables below this point are internal to this subroutine.

      INTEGER i, j

      REAL ob_sq3, v_endpts ( 4, 3), a1 ( 3), a2 ( 3),
     1     Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz, t

      REAL, DIMENSION ( 0: ntp_patches, 3) ::
     1                arr_in_sector_1, arr_in_sector_2,
     1                arr_in_sector_3, arr_in_sector_4

      REAL, DIMENSION ( 0: ntp_zones, 3) ::
     1                arr_in_patch_1, arr_in_patch_2,
     1                arr_in_patch_3, arr_in_patch_4

! ----------------------------------------------------------------------

      ob_sq3 = 1.0 / SQRT ( 3.0)

      IF ( ( isector .LT. 1) .OR. ( isector .GT. 6) ) THEN
      WRITE ( 6, *)
     1 "Stopping in Mesh_Sector_General, isector = ", isector
      STOP
      END IF

      IF ( ( ix_patch .LT. 1) .OR. ( ix_patch .GT. ntp_patches) ) THEN
      WRITE ( 6, *)
     1 "Stopping in Mesh_Sector_General, ix_patch, ntp_patches = ",
     1             ix_patch, ntp_patches
      STOP
      END IF

      IF ( ( iy_patch .LT. 1) .OR. ( iy_patch .GT. ntp_patches) ) THEN
      WRITE ( 6, *)
     1 "Stopping in Mesh_Sector_General, iy_patch, ntp_patches = ",
     1             iy_patch, ntp_patches
      STOP
      END IF


!     WRITE ( 6, *) "NGhost ntp_zones", NGhost, ntp_zones


! ----------------------------------------------------------------------
! STEP 1) Build the end points of the chosen sector :- "v_endpts ( 1, :)"
! effectively gives us the origin. "v_endpts ( 2, :) - v_endpts ( 1, :)"
! forms the effective x-axis. "v_endpts ( 3, :) - v_endpts ( 1, :)" forms the
! effective y-axis. "v_endpts ( 4, :)" gives the remaining vertex.
! This is a right-handed basis.

! ----------------------------------------

      IF ( isector .EQ. 1) THEN

! Sector 1 is based on (z^hat) = (x^hat) X (y^hat).

! Local x-direction is global x-direction; local y-direction is global y-direction.
! Neighbors in local x-direction :- 6, 3; Neighbors in local y-direction :- 5, 2.

! The vertices are (-1,-1,1), (1,-1,1), (-1,1,1) and (1,1,1).

! ----------------------------------------

        v_endpts ( 1, 1) = - ob_sq3 ; v_endpts ( 1, 2) = - ob_sq3
        v_endpts ( 1, 3) = ob_sq3

        v_endpts ( 2, 1) = ob_sq3 ; v_endpts ( 2, 2) = - ob_sq3
        v_endpts ( 2, 3) = ob_sq3

        v_endpts ( 3, 1) = - ob_sq3 ; v_endpts ( 3, 2) = ob_sq3
        v_endpts ( 3, 3) = ob_sq3

        v_endpts ( 4, 1) = ob_sq3 ; v_endpts ( 4, 2) = ob_sq3
        v_endpts ( 4, 3) = ob_sq3

! ----------------------------------------

      ELSE IF ( isector .EQ. 2) THEN

! Sector 2 is based on (y^hat) = (z^hat) X (x^hat).

! Local x-direction is global z-direction; local y-direction is global x-direction.
! Neighbors in local x-direction :- 4, 1; Neighbors in local y-direction :- 6, 3.

! The vertices are (-1,1,-1), (-1,1,1), (1,1,-1) and (1,1,1).

! ----------------------------------------

        v_endpts ( 1, 1) = - ob_sq3 ; v_endpts ( 1, 2) = ob_sq3
        v_endpts ( 1, 3) = - ob_sq3

        v_endpts ( 2, 1) = - ob_sq3 ; v_endpts ( 2, 2) = ob_sq3
        v_endpts ( 2, 3) = ob_sq3

        v_endpts ( 3, 1) = ob_sq3 ; v_endpts ( 3, 2) = ob_sq3
        v_endpts ( 3, 3) = - ob_sq3

        v_endpts ( 4, 1) = ob_sq3 ; v_endpts ( 4, 2) = ob_sq3
        v_endpts ( 4, 3) = ob_sq3

! ----------------------------------------

      ELSE IF ( isector .EQ. 3) THEN

! Sector 3 is based on (x^hat) = (y^hat) X (z^hat).

! Local x-direction is global y-direction; local y-direction is global z-direction.
! Neighbors in local x-direction :- 5, 2; Neighbors in local y-direction :- 4, 1.

! The vertices are (1,-1,-1), (1,1,-1), (1,-1,1) and (1,1,1).

! ----------------------------------------

        v_endpts ( 1, 1) = ob_sq3 ; v_endpts ( 1, 2) = - ob_sq3
        v_endpts ( 1, 3) = - ob_sq3

        v_endpts ( 2, 1) = ob_sq3 ; v_endpts ( 2, 2) = ob_sq3
        v_endpts ( 2, 3) = - ob_sq3

        v_endpts ( 3, 1) = ob_sq3 ; v_endpts ( 3, 2) = - ob_sq3
        v_endpts ( 3, 3) = ob_sq3

        v_endpts ( 4, 1) = ob_sq3 ; v_endpts ( 4, 2) = ob_sq3
        v_endpts ( 4, 3) = ob_sq3

! ----------------------------------------

      ELSE IF ( isector .EQ. 4) THEN

! Sector 4 is based on (-z^hat) = (-y^hat) X (-x^hat).

! Local x-direction is global -y-direction; local y-direction is global -x-direction.
! Neighbors in local x-direction :- 2, 5; Neighbors in local y-direction :- 3, 6.

! The vertices are (1,1,-1), (1,-1,-1), (-1,1,-1) and (-1,-1,-1).

! ----------------------------------------

        v_endpts ( 1, 1) = ob_sq3 ; v_endpts ( 1, 2) = ob_sq3
        v_endpts ( 1, 3) = - ob_sq3

        v_endpts ( 2, 1) = ob_sq3 ; v_endpts ( 2, 2) = - ob_sq3
        v_endpts ( 2, 3) = - ob_sq3

        v_endpts ( 3, 1) = - ob_sq3 ; v_endpts ( 3, 2) = ob_sq3
        v_endpts ( 3, 3) = - ob_sq3

        v_endpts ( 4, 1) = - ob_sq3 ; v_endpts ( 4, 2) = - ob_sq3
        v_endpts ( 4, 3) = - ob_sq3

! ----------------------------------------

      ELSE IF ( isector .EQ. 5) THEN

! Sector 5 is based on (-y^hat) = (-x^hat) X (-z^hat).

! Local x-direction is global -x-direction; local y-direction is global -z-direction.
! Neighbors in local x-direction :- 3, 6; Neighbors in local y-direction :- 1, 4.

! The vertices are (1,-1,1), (-1,-1,1), (1,-1,-1) and (-1,-1,-1).

! ----------------------------------------

        v_endpts ( 1, 1) = ob_sq3 ; v_endpts ( 1, 2) = - ob_sq3
        v_endpts ( 1, 3) = ob_sq3

        v_endpts ( 2, 1) = - ob_sq3 ; v_endpts ( 2, 2) = - ob_sq3
        v_endpts ( 2, 3) = ob_sq3

        v_endpts ( 3, 1) = ob_sq3 ; v_endpts ( 3, 2) = - ob_sq3
        v_endpts ( 3, 3) = - ob_sq3

        v_endpts ( 4, 1) = - ob_sq3 ; v_endpts ( 4, 2) = - ob_sq3
        v_endpts ( 4, 3) = - ob_sq3

! ----------------------------------------

      ELSE IF ( isector .EQ. 6) THEN

! Sector 6 is based on (-x^hat) = (-z^hat) X (-y^hat).

! Local x-direction is global -z-direction; local y-direction is global -y-direction.
! Neighbors in local x-direction :- 1, 4; Neighbors in local y-direction :- 2, 5.

! The vertices are (-1,1,1), (-1,1,-1), (-1,-1,1) and (-1,-1,-1).

! ----------------------------------------

        v_endpts ( 1, 1) = - ob_sq3 ; v_endpts ( 1, 2) = ob_sq3
        v_endpts ( 1, 3) = ob_sq3

        v_endpts ( 2, 1) = - ob_sq3 ; v_endpts ( 2, 2) = ob_sq3
        v_endpts ( 2, 3) = - ob_sq3

        v_endpts ( 3, 1) = - ob_sq3 ; v_endpts ( 3, 2) = - ob_sq3
        v_endpts ( 3, 3) = ob_sq3

        v_endpts ( 4, 1) = - ob_sq3 ; v_endpts ( 4, 2) = - ob_sq3
        v_endpts ( 4, 3) = - ob_sq3

! ----------------------------------------

      END IF

! ----------------------------------------------------------------------
! STEP 2) Now that "v_endpts ( 1:4, 1:3)" is built for this sector, we can build
! the endpoints of the present patch on the surface of the desired sector of the
! unit sphere. In other words, we fill in "vertices ( 1:4, 1:3)" which is the actual
! thing that is returned through the subroutine boundary.

! --------------------------------------------------
! If "ntp_patches == 1", then there is no need for the following code.

      IF ( ntp_patches .EQ. 1) THEN

        vertices ( 1:4, 1:3) = v_endpts ( 1:4, 1:3)
        GO TO 1000

      END IF

! --------------------------------------------------
! Swivel from "v_endpts ( 1, :)" to "v_endpts ( 2, :)" to build
! "arr_in_sector_1 ( 0: ntp_patches, 3)"

      a1 ( :) = v_endpts ( 1, :)
      a2 ( :) = v_endpts ( 2, :)
      CALL SWIVEL_N ( a1, a2, ntp_patches, arr_in_sector_1)

! ----------------------------------------
! Swivel from "v_endpts ( 3, :)" to "v_endpts ( 4, :)" to build
! "arr_in_sector_2 ( 0: ntp_patches, 3)"

      a1 ( :) = v_endpts ( 3, :)
      a2 ( :) = v_endpts ( 4, :)
      CALL SWIVEL_N ( a1, a2, ntp_patches, arr_in_sector_2)

! ----------------------------------------
! Swivel from "v_endpts ( 1, :)" to "v_endpts ( 3, :)" to build
! "arr_in_sector_3 ( 0: ntp_patches, 3)"

      a1 ( :) = v_endpts ( 1, :)
      a2 ( :) = v_endpts ( 3, :)
      CALL SWIVEL_N ( a1, a2, ntp_patches, arr_in_sector_3)

! ----------------------------------------
! Swivel from "v_endpts ( 2, :)" to "v_endpts ( 4, :)" to build
! "arr_in_sector_4 ( 0: ntp_patches, 3)"

      a1 ( :) = v_endpts ( 2, :)
      a2 ( :) = v_endpts ( 4, :)
      CALL SWIVEL_N ( a1, a2, ntp_patches, arr_in_sector_4)

! --------------------------------------------------
! Now use "arr_in_sector_1" to "arr_in_sector_4" at the boundary of the patch
! that we are considering to form the normals to the two planes. Then take their
! cross product to obtain "vertices ( 1:4, 1:3)".

! --------------------------------------------------
! Focus on "vertices ( 1, 1:3)".

      i = ix_patch - 1
      j = iy_patch - 1

! ------------------------------

! Build vector that is orthogonal to plane formed by "arr_in_sector_1 ( i, :)"
! and "arr_in_sector_2 ( i, :)".

        Ax = arr_in_sector_1 ( i, 2) * arr_in_sector_2 ( i, 3)
     1     - arr_in_sector_1 ( i, 3) * arr_in_sector_2 ( i, 2)

        Ay = arr_in_sector_1 ( i, 3) * arr_in_sector_2 ( i, 1)
     1     - arr_in_sector_1 ( i, 1) * arr_in_sector_2 ( i, 3)

        Az = arr_in_sector_1 ( i, 1) * arr_in_sector_2 ( i, 2)
     1     - arr_in_sector_1 ( i, 2) * arr_in_sector_2 ( i, 1)

! ------------------------------

! Build vector that is orthogonal to plane formed by "arr_in_sector_3 ( j, :)"
! and "arr_in_sector_4 ( j, :)".

        Bx = arr_in_sector_3 ( j, 2) * arr_in_sector_4 ( j, 3)
     1     - arr_in_sector_3 ( j, 3) * arr_in_sector_4 ( j, 2)

        By = arr_in_sector_3 ( j, 3) * arr_in_sector_4 ( j, 1)
     1     - arr_in_sector_3 ( j, 1) * arr_in_sector_4 ( j, 3)

        Bz = arr_in_sector_3 ( j, 1) * arr_in_sector_4 ( j, 2)
     1     - arr_in_sector_3 ( j, 2) * arr_in_sector_4 ( j, 1)

! ------------------------------

! Now take the cross product of the two vectors, flip for sector number if needed,
! and normalize the final vector. This gives us "coordinates_gr ( i, j, :)".

        Cx = By * Az - Bz * Ay
        Cy = Bz * Ax - Bx * Az
        Cz = Bx * Ay - By * Ax

        IF ( ( isector .EQ. 1) .AND. ( Cz .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 4) .AND. ( Cz .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 2) .AND. ( Cy .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 5) .AND. ( Cy .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 3) .AND. ( Cx .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 6) .AND. ( Cx .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

! ------------------------------

        t = 1.0 / SQRT ( Cx**2 + Cy**2 + Cz**2)
        Cx = t * Cx
        Cy = t * Cy
        Cz = t * Cz

        vertices ( 1, 1) = Cx
        vertices ( 1, 2) = Cy
        vertices ( 1, 3) = Cz

! --------------------------------------------------
! Focus on "vertices ( 2, 1:3)".

      i = ix_patch
      j = iy_patch - 1

! ------------------------------

! Build vector that is orthogonal to plane formed by "arr_in_sector_1 ( i, :)"
! and "arr_in_sector_2 ( i, :)".

        Ax = arr_in_sector_1 ( i, 2) * arr_in_sector_2 ( i, 3)
     1     - arr_in_sector_1 ( i, 3) * arr_in_sector_2 ( i, 2)

        Ay = arr_in_sector_1 ( i, 3) * arr_in_sector_2 ( i, 1)
     1     - arr_in_sector_1 ( i, 1) * arr_in_sector_2 ( i, 3)

        Az = arr_in_sector_1 ( i, 1) * arr_in_sector_2 ( i, 2)
     1     - arr_in_sector_1 ( i, 2) * arr_in_sector_2 ( i, 1)

! ------------------------------

! Build vector that is orthogonal to plane formed by "arr_in_sector_3 ( j, :)"
! and "arr_in_sector_4 ( j, :)".

        Bx = arr_in_sector_3 ( j, 2) * arr_in_sector_4 ( j, 3)
     1     - arr_in_sector_3 ( j, 3) * arr_in_sector_4 ( j, 2)

        By = arr_in_sector_3 ( j, 3) * arr_in_sector_4 ( j, 1)
     1     - arr_in_sector_3 ( j, 1) * arr_in_sector_4 ( j, 3)

        Bz = arr_in_sector_3 ( j, 1) * arr_in_sector_4 ( j, 2)
     1     - arr_in_sector_3 ( j, 2) * arr_in_sector_4 ( j, 1)

! ------------------------------

! Now take the cross product of the two vectors, flip for sector number if needed,
! and normalize the final vector. This gives us "coordinates_gr ( i, j, :)".

        Cx = By * Az - Bz * Ay
        Cy = Bz * Ax - Bx * Az
        Cz = Bx * Ay - By * Ax

        IF ( ( isector .EQ. 1) .AND. ( Cz .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 4) .AND. ( Cz .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 2) .AND. ( Cy .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 5) .AND. ( Cy .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 3) .AND. ( Cx .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 6) .AND. ( Cx .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

! ------------------------------

        t = 1.0 / SQRT ( Cx**2 + Cy**2 + Cz**2)
        Cx = t * Cx
        Cy = t * Cy
        Cz = t * Cz

        vertices ( 2, 1) = Cx
        vertices ( 2, 2) = Cy
        vertices ( 2, 3) = Cz

! --------------------------------------------------
! Focus on "vertices ( 3, 1:3)".

      i = ix_patch - 1
      j = iy_patch

! ------------------------------

! Build vector that is orthogonal to plane formed by "arr_in_sector_1 ( i, :)"
! and "arr_in_sector_2 ( i, :)".

        Ax = arr_in_sector_1 ( i, 2) * arr_in_sector_2 ( i, 3)
     1     - arr_in_sector_1 ( i, 3) * arr_in_sector_2 ( i, 2)

        Ay = arr_in_sector_1 ( i, 3) * arr_in_sector_2 ( i, 1)
     1     - arr_in_sector_1 ( i, 1) * arr_in_sector_2 ( i, 3)

        Az = arr_in_sector_1 ( i, 1) * arr_in_sector_2 ( i, 2)
     1     - arr_in_sector_1 ( i, 2) * arr_in_sector_2 ( i, 1)

! ------------------------------

! Build vector that is orthogonal to plane formed by "arr_in_sector_3 ( j, :)"
! and "arr_in_sector_4 ( j, :)".

        Bx = arr_in_sector_3 ( j, 2) * arr_in_sector_4 ( j, 3)
     1     - arr_in_sector_3 ( j, 3) * arr_in_sector_4 ( j, 2)

        By = arr_in_sector_3 ( j, 3) * arr_in_sector_4 ( j, 1)
     1     - arr_in_sector_3 ( j, 1) * arr_in_sector_4 ( j, 3)

        Bz = arr_in_sector_3 ( j, 1) * arr_in_sector_4 ( j, 2)
     1     - arr_in_sector_3 ( j, 2) * arr_in_sector_4 ( j, 1)

! ------------------------------

! Now take the cross product of the two vectors, flip for sector number if needed,
! and normalize the final vector. This gives us "coordinates_gr ( i, j, :)".

        Cx = By * Az - Bz * Ay
        Cy = Bz * Ax - Bx * Az
        Cz = Bx * Ay - By * Ax

        IF ( ( isector .EQ. 1) .AND. ( Cz .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 4) .AND. ( Cz .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 2) .AND. ( Cy .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 5) .AND. ( Cy .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 3) .AND. ( Cx .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 6) .AND. ( Cx .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

! ------------------------------

        t = 1.0 / SQRT ( Cx**2 + Cy**2 + Cz**2)
        Cx = t * Cx
        Cy = t * Cy
        Cz = t * Cz

        vertices ( 3, 1) = Cx
        vertices ( 3, 2) = Cy
        vertices ( 3, 3) = Cz

! --------------------------------------------------
! Focus on "vertices ( 4, 1:3)".

      i = ix_patch
      j = iy_patch

! ------------------------------

! Build vector that is orthogonal to plane formed by "arr_in_sector_1 ( i, :)"
! and "arr_in_sector_2 ( i, :)".

        Ax = arr_in_sector_1 ( i, 2) * arr_in_sector_2 ( i, 3)
     1     - arr_in_sector_1 ( i, 3) * arr_in_sector_2 ( i, 2)

        Ay = arr_in_sector_1 ( i, 3) * arr_in_sector_2 ( i, 1)
     1     - arr_in_sector_1 ( i, 1) * arr_in_sector_2 ( i, 3)

        Az = arr_in_sector_1 ( i, 1) * arr_in_sector_2 ( i, 2)
     1     - arr_in_sector_1 ( i, 2) * arr_in_sector_2 ( i, 1)

! ------------------------------

! Build vector that is orthogonal to plane formed by "arr_in_sector_3 ( j, :)"
! and "arr_in_sector_4 ( j, :)".

        Bx = arr_in_sector_3 ( j, 2) * arr_in_sector_4 ( j, 3)
     1     - arr_in_sector_3 ( j, 3) * arr_in_sector_4 ( j, 2)

        By = arr_in_sector_3 ( j, 3) * arr_in_sector_4 ( j, 1)
     1     - arr_in_sector_3 ( j, 1) * arr_in_sector_4 ( j, 3)

        Bz = arr_in_sector_3 ( j, 1) * arr_in_sector_4 ( j, 2)
     1     - arr_in_sector_3 ( j, 2) * arr_in_sector_4 ( j, 1)

! ------------------------------

! Now take the cross product of the two vectors, flip for sector number if needed,
! and normalize the final vector. This gives us "coordinates_gr ( i, j, :)".

        Cx = By * Az - Bz * Ay
        Cy = Bz * Ax - Bx * Az
        Cz = Bx * Ay - By * Ax

        IF ( ( isector .EQ. 1) .AND. ( Cz .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 4) .AND. ( Cz .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 2) .AND. ( Cy .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 5) .AND. ( Cy .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 3) .AND. ( Cx .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 6) .AND. ( Cx .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

! ------------------------------

        t = 1.0 / SQRT ( Cx**2 + Cy**2 + Cz**2)
        Cx = t * Cx
        Cy = t * Cy
        Cz = t * Cz

        vertices ( 4, 1) = Cx
        vertices ( 4, 2) = Cy
        vertices ( 4, 3) = Cz

! --------------------------------------------------
! Get here if "ntp_patches == 1".

1000  CONTINUE

! ----------------------------------------------------------------------
! STEP 3) Now that "vertices ( 1:4, 1:3)" is built for this patch, we can build
! the entire 2D mesh on the surface of the desired sector of the unit sphere.

! --------------------------------------------------
! Swivel from "vertices ( 1, :)" to "vertices ( 2, :)" to build
! "arr_in_patch_1 ( 0: ntp_zones, 3)"

      a1 ( :) = vertices ( 1, :)
      a2 ( :) = vertices ( 2, :)
      CALL SWIVEL_N ( a1, a2, ntp_zones, arr_in_patch_1)

! ----------------------------------------
! Swivel from "vertices ( 3, :)" to "vertices ( 4, :)" to build
! "arr_in_patch_2 ( 0: ntp_zones, 3)"

      a1 ( :) = vertices ( 3, :)
      a2 ( :) = vertices ( 4, :)
      CALL SWIVEL_N ( a1, a2, ntp_zones, arr_in_patch_2)

! ----------------------------------------
! Swivel from "vertices ( 1, :)" to "vertices ( 3, :)" to build
! "arr_in_patch_3 ( 0: ntp_zones, 3)"

      a1 ( :) = vertices ( 1, :)
      a2 ( :) = vertices ( 3, :)
      CALL SWIVEL_N ( a1, a2, ntp_zones, arr_in_patch_3)

! ----------------------------------------
! Swivel from "vertices ( 2, :)" to "vertices ( 4, :)" to build
! "arr_in_patch_4 ( 0: ntp_zones, 3)"

      a1 ( :) = vertices ( 2, :)
      a2 ( :) = vertices ( 4, :)
      CALL SWIVEL_N ( a1, a2, ntp_zones, arr_in_patch_4)

! --------------------------------------------------
! Now use "arr_in_patch_1" to "arr_in_patch_4" at the boundary of the patch
! that we are considering to form the normals to the two planes. Then take their
! cross product to obtain "coordinates_gr ( i, j, :)".

      DO j = 0, ntp_zones

        DO i = 0, ntp_zones

! ------------------------------

! Build vector that is orthogonal to plane formed by "arr_in_patch_1 ( i, :)"
! and "arr_in_patch_2 ( i, :)".

        Ax = arr_in_patch_1 ( i, 2) * arr_in_patch_2 ( i, 3)
     1     - arr_in_patch_1 ( i, 3) * arr_in_patch_2 ( i, 2)

        Ay = arr_in_patch_1 ( i, 3) * arr_in_patch_2 ( i, 1)
     1     - arr_in_patch_1 ( i, 1) * arr_in_patch_2 ( i, 3)

        Az = arr_in_patch_1 ( i, 1) * arr_in_patch_2 ( i, 2)
     1     - arr_in_patch_1 ( i, 2) * arr_in_patch_2 ( i, 1)

! ------------------------------

! Build vector that is orthogonal to plane formed by "arr_in_patch_3 ( j, :)"
! and "arr_in_patch_4 ( j, :)".

        Bx = arr_in_patch_3 ( j, 2) * arr_in_patch_4 ( j, 3)
     1     - arr_in_patch_3 ( j, 3) * arr_in_patch_4 ( j, 2)

        By = arr_in_patch_3 ( j, 3) * arr_in_patch_4 ( j, 1)
     1     - arr_in_patch_3 ( j, 1) * arr_in_patch_4 ( j, 3)

        Bz = arr_in_patch_3 ( j, 1) * arr_in_patch_4 ( j, 2)
     1     - arr_in_patch_3 ( j, 2) * arr_in_patch_4 ( j, 1)

! ------------------------------

! Now take the cross product of the two vectors, flip for sector number if needed,
! and normalize the final vector. This gives us "coordinates_gr ( i, j, :)".

        Cx = By * Az - Bz * Ay
        Cy = Bz * Ax - Bx * Az
        Cz = Bx * Ay - By * Ax

        IF ( ( isector .EQ. 1) .AND. ( Cz .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 4) .AND. ( Cz .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 2) .AND. ( Cy .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 5) .AND. ( Cy .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 3) .AND. ( Cx .LT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

        IF ( ( isector .EQ. 6) .AND. ( Cx .GT. 0.0) ) THEN
          Cx = - Cx   ;   Cy = - Cy   ;   Cz = - Cz
        END IF

! ------------------------------

        t = 1.0 / SQRT ( Cx**2 + Cy**2 + Cz**2)
        Cx = t * Cx
        Cy = t * Cy
        Cz = t * Cz

        coordinates_gr ( i, j, 1) = Cx
        coordinates_gr ( i, j, 2) = Cy
        coordinates_gr ( i, j, 3) = Cz

! ------------------------------

        END DO  ! End of "DO i = 0, ntp_zones" loop.

      END DO  ! End of "DO j = 0, ntp_zones" loop.

! ----------------------------------------------------------------------

      END SUBROUTINE Mesh_Sector_General