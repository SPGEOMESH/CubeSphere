



      SUBROUTINE PRINT_VECT ( istart, iend, vector, print_fmt )

! ----------------------------------------------------------------------
! Written by: Sethupathy S.
!
! Prints a 1D array in a nice formatted way. Calling this function
! in DO-loop can print a 2D-array/Matrix in a nice easy looking format.
!
! Example:
!       DO iyy = ntp_zones + NGhost, 1 - NGhost, -1
!
!        Bxf_1D  ( 0 - NGhost: ntp_zones + NGhost) =
!     1  Bxf_amr ( 0 - NGhost: ntp_zones + NGhost, iyy, izz, ii)
!
!        CALL PRINT_VECT ( 0 - NGhost, ntp_zones + NGhost, Bxf_1D, "R" )
!
!      END DO
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER istart, iend

      REAL vector ( istart: iend)

      CHARACTER(LEN = 1) :: print_fmt

! ----------------------------------------------------------------------

      INTEGER i
      INTEGER :: vec ( istart: iend)

! ----------
      IF ( print_fmt == "R" ) THEN          ! Print as a REAL number
! ----------

      DO i = istart, iend
        WRITE ( 6, '(f4.1, 2x)', advance="no")vector ( i)
      END DO

! ----------
      ELSEIF ( print_fmt == "I" ) THEN      ! Print as an Integer
! ----------

      DO i = istart, iend
        vec ( i) = NINT ( vector ( i))
        WRITE ( 6, '(i4.1, 2x)', advance="no") vec ( i)
      END DO

! ----------
      ELSE
! ----------

      WRITE ( 6, *)"PRINT_VECT ( ): Error in specifying printing format"

! ----------
      ENDIF
! ----------

      WRITE(6, *)" "

! ----------------------------------------------------------------------

      END SUBROUTINE PRINT_VECT

