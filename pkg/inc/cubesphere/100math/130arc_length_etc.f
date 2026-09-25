





! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      REAL FUNCTION Scalar_Pot ( v )

      IMPLICIT NONE

      REAL :: v ( 3)

      Scalar_Pot = v ( 1) * v ( 2) * v ( 3)

      RETURN

      END FUNCTION

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      REAL FUNCTION arc_length ( v1, v2 )

      IMPLICIT NONE

      REAL :: v1 ( 3), v2 ( 3)

      REAL :: v1_dot_v2, v1_mag, v2_mag, angle

      v1_mag = SQRT( v1 ( 1)**2 + v1 ( 2)**2 + v1 ( 3)**2)
      v2_mag = SQRT( v2 ( 1)**2 + v2 ( 2)**2 + v2 ( 3)**2)

      v1_dot_v2 = ( v1 ( 1) * v2 ( 1) + v1 ( 2) * v2 ( 2)
     1            + v1 ( 3) * v2 ( 3))

      angle = ACOS ( v1_dot_v2 / ( v1_mag * v2_mag))

      arc_length = v1_mag * angle      ! arc length = r * theta

      RETURN

      END FUNCTION

