

! for silly code linking tests



! Add x and y
      REAL FUNCTION ADD ( x, y)

      IMPLICIT NONE

      REAL x, y

      ADD = x + y

      RETURN
      END


! Multiply x and y
      REAL FUNCTION MUL ( x, y)

      IMPLICIT NONE

      REAL x, y

      MUL = x * y

      RETURN
      END

! Subtract x and y (as subroutine)
      SUBROUTINE SUB ( x, y, z)

      IMPLICIT NONE

      REAL x, y, z

      z = x - y

      END

