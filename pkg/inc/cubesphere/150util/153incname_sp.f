
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INCNAME_SP ( filname, inc)

! ----------------------------------------------------------------------

      IMPLICIT NONE

      EXTERNAL MYSTOP

      CHARACTER*8 filname

      CHARACTER*1 unitchar ( 8)

      INTEGER unitint ( 8), zeroint, value, inc, i, j

! ----------------------------------------------------------------------

      unitchar ( 1) = filname ( 1:1)
      unitchar ( 2) = filname ( 2:2)
      unitchar ( 3) = filname ( 3:3)
      unitchar ( 4) = filname ( 4:4)
      unitchar ( 5) = filname ( 5:5)
      unitchar ( 6) = filname ( 6:6)
      unitchar ( 7) = filname ( 7:7)
      unitchar ( 8) = filname ( 8:8)

      zeroint = ichar ( '0')

      unitint ( 4) = ichar ( unitchar ( 4) ) - zeroint
      unitint ( 3) = ichar ( unitchar ( 3) ) - zeroint

      value = 10 * unitint ( 3) + unitint ( 4)

      value = value + inc

      unitint ( 3) = value / 10

      value = value - 10 * unitint ( 3)

      unitint ( 4) = value

      unitint ( 3) = unitint ( 3) + zeroint
      unitint ( 4) = unitint ( 4) + zeroint

      unitchar ( 3) = char ( unitint ( 3) )
      unitchar ( 4) = char ( unitint ( 4) )

      filname ( 3:3) = unitchar ( 3)
      filname ( 4:4) = unitchar ( 4)

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INCNAME_SP

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
