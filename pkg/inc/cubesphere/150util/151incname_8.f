

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE INCNAME_8 ( filname, inc)

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

      unitint ( 8) = ichar ( unitchar ( 8) ) - zeroint
      unitint ( 7) = ichar ( unitchar ( 7) ) - zeroint
      unitint ( 6) = ichar ( unitchar ( 6) ) - zeroint
      unitint ( 5) = ichar ( unitchar ( 5) ) - zeroint

      value = 1000 * unitint ( 5) + 100 * unitint ( 6)
     1              + 10 * unitint ( 7) + unitint ( 8)

      value = value + inc


      unitint ( 5) = value / 1000

      value = value - 1000 * unitint ( 5)

      unitint ( 6) = value / 100

      value = value - 100 * unitint ( 6)

      unitint ( 7) = value / 10

      value = value - 10 * unitint ( 7)

      unitint ( 8) = value

      unitint ( 5) = unitint ( 5) + zeroint
      unitint ( 6) = unitint ( 6) + zeroint
      unitint ( 7) = unitint ( 7) + zeroint
      unitint ( 8) = unitint ( 8) + zeroint

      unitchar ( 5) = char ( unitint ( 5) )
      unitchar ( 6) = char ( unitint ( 6) )
      unitchar ( 7) = char ( unitint ( 7) )
      unitchar ( 8) = char ( unitint ( 8) )

      filname ( 5:5) = unitchar ( 5)
      filname ( 6:6) = unitchar ( 6)
      filname ( 7:7) = unitchar ( 7)
      filname ( 8:8) = unitchar ( 8)

! ----------------------------------------------------------------------

      RETURN

      END SUBROUTINE INCNAME_8