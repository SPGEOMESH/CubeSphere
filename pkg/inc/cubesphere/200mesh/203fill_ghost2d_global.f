





      SUBROUTINE FILL_GHOST2D_GLOBAL

     1           ( isector, ntp_global, NGhost_global,
     1             coordinates_gr_global)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! The interiors of "coordinates_gr_global" can be filled in by using
! "Mesh_Sector_General" for each of the sectors. But, for a cube sphere,
! the ghost zones on each sector need information from the neighboring sectors.
! This subroutine pulls in such information and initializes it correctly.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER isector, ntp_global, NGhost_global

      REAL coordinates_gr_global
     1     ( 0 - NGhost_global: ntp_global + NGhost_global,
     1       0 - NGhost_global: ntp_global + NGhost_global, 3, 6)

! ----------------------------------------------------------------------
! Variables internal to this subroutine.

      INTEGER ixx, iyy, ixx_you, iyy_you,
     1        isrt, islf, isup, isdn

! ----------------------------------------------------------------------

      IF ( NGhost_global .GE. ntp_global) THEN

      WRITE ( 6, *)"Stopping in FILL_GHOST2D_GLOBAL NGhost_global = ",
     1              NGhost_global

      STOP
      END IF

! ----------------------------------------------------------------------

      IF ( isector .EQ. 1) THEN

      isrt = 3
      islf = 6
      isup = 2
      isdn = 5

      ELSE IF ( isector .EQ. 2) THEN

      isrt = 1
      islf = 4
      isup = 3
      isdn = 6

      ELSE IF ( isector .EQ. 3) THEN

      isrt = 2
      islf = 5
      isup = 1
      isdn = 4

      ELSE IF ( isector .EQ. 4) THEN

      isrt = 5
      islf = 2
      isup = 6
      isdn = 3

      ELSE IF ( isector .EQ. 5) THEN

      isrt = 6
      islf = 3
      isup = 4
      isdn = 1

      ELSE IF ( isector .EQ. 6) THEN

      isrt = 4
      islf = 1
      isup = 5
      isdn = 2

      ELSE

      WRITE ( 6, *)"Stopping in FILL_GHOST2D_GLOBAL, isector = ",
     1             isector
      STOP

      END IF

! ----------------------------------------------------------------------
! Fill up ghost zones on "isector" using sector "isrt" to the right of it.

      DO ixx = ntp_global + 1, ntp_global + NGhost_global
      iyy_you = ntp_global - 1 - ( ixx - ntp_global - 1)

        DO iyy = 0, ntp_global
        ixx_you = iyy

        coordinates_gr_global ( ixx, iyy, :, isector) =
     1  coordinates_gr_global ( ixx_you, iyy_you, :, isrt)

        END DO  ! End of "DO iyy = 0 - NGhost_global, ntp_global + NGhost_global" loop.

      END DO  ! End of "DO ixx = ntp_global + 1, ntp_global + NGhost_global" loop.

! ----------------------------------------------------------------------
! Fill up ghost zones on "isector" using sector "islf" to the left of it.

      DO ixx = -1, 0 - NGhost_global, -1
      ixx_you = 1 - ( ixx + 1)

        DO iyy = 0, ntp_global
        iyy_you = ntp_global - iyy

        coordinates_gr_global ( ixx, iyy, :, isector) =
     1  coordinates_gr_global ( ixx_you, iyy_you, :, islf)

        END DO  ! End of "DO iyy = 0 - NGhost_global, ntp_global + NGhost_global" loop.

      END DO  ! End of "DO ixx = -1, 0 - NGhost_global, -1" loop.

! ----------------------------------------------------------------------
! Fill up ghost zones on "isector" using sector "isup" above it.

      DO iyy = ntp_global + 1, ntp_global + NGhost_global
      ixx_you = ntp_global - 1 - ( iyy - ntp_global - 1)

        DO ixx = 0, ntp_global
        iyy_you = ixx

       coordinates_gr_global ( ixx, iyy, :, isector) =
     1  coordinates_gr_global ( ixx_you, iyy_you, :, isup)

        END DO  ! End of "DO ixx = 0 - NGhost_global, ntp_global + NGhost_global" loop.

      END DO  ! End of "DO iyy = ntp_global + 1, ntp_global + NGhost_global" loop.

! ----------------------------------------------------------------------
! Fill up ghost zones on "isector" using sector "isdn" below it.

      DO iyy = -1, 0 - NGhost_global, -1
      iyy_you = 1 - ( iyy + 1)

        DO ixx = 0, ntp_global
        ixx_you = ntp_global - ixx

        coordinates_gr_global ( ixx, iyy, :, isector) =
     1  coordinates_gr_global ( ixx_you, iyy_you, :, isdn)

        END DO  ! End of "DO ixx = 0 - NGhost_global, ntp_global + NGhost_global" loop.

      END DO  ! End of "DO iyy = -1, 0 - NGhost_global, -1" loop.

! ----------------------------------------------------------------------

      END SUBROUTINE FILL_GHOST2D_GLOBAL



