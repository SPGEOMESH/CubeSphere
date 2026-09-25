







      SUBROUTINE PRINT_VARIABLES_ONE_PATCH_PER_SECTOR ( isector )

! ----------------------------------------------------------------------
! Written By: Sethupathy S.

! Prints the values of zones assuming 1 patch per sector.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER isector

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ipatch, ii, i, j, ixx, iyy, izz, icomp

       REAL U_1D ( 1 - NGhost: ntp_zones + NGhost),
     1      Bxf_1D ( 0 - NGhost: ntp_zones + NGhost),
     1      Byf_1D ( 0 - NGhost: ntp_zones + NGhost),
     1      Bzf_1D ( 0 - NGhost: ntp_zones + NGhost)


! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------

      IF ( me == isector) THEN

! ----------------------------------------------------------------------

      IF ( npatch_per_proc .NE. 1 ) THEN
        WRITE ( 6, *)"Stopping at PRINT_VARIABLES_ONE_PATCH_PER_SECTOR"
        WRITE ( 6, *)"Expecting npatch_per_proc == 1"
        WRITE ( 6, *)"npatch_per_proc = ",npatch_per_proc
        STOP
      END IF

      IF ( nimages .NE. 6 ) THEN
        WRITE ( 6, *)"Stopping at PRINT_VARIABLES_ONE_PATCH_PER_SECTOR"
        WRITE ( 6, *)"Expecting nimages == 6"
        WRITE ( 6, *)"nimages = ",nimages
        STOP
      END IF

      ii = 1

! ----------------------------------------------------------------------

      IF ( isector == iactive_FromProc ( ii, 1)) THEN

! ------------------------------------------------------------

      izz = 1
      icomp = 1

! ----------

      WRITE ( 6, *)"U_amr"

      DO iyy = ntp_zones + NGhost, 1 - NGhost, -1

        U_1D ( :) = U_amr ( 1 - NGhost: ntp_zones + NGhost,
     1                      iyy, izz, icomp, ii)

        CALL PRINT_VECT ( 1 - NGhost, ntp_zones + NGhost,
     1             U_1D ( 1 - NGhost: ntp_zones + NGhost), "I")

      END DO

      WRITE ( 6, *)"-----------------------------------------------"

! ----------

      WRITE ( 6, *)"Bxf_amr"

      DO iyy = ntp_zones + NGhost, 1 - NGhost, -1

        Bxf_1D  ( 0 - NGhost: ntp_zones + NGhost) =
     1  Bxf_amr ( 0 - NGhost: ntp_zones + NGhost, iyy, izz, ii)

        CALL PRINT_VECT ( 0 - NGhost, ntp_zones + NGhost, Bxf_1D, "R" )

      END DO

      WRITE ( 6, *)"-----------------------------------------------"

! ----------

      WRITE ( 6, *)"Byf_amr"

      DO iyy = ntp_zones + NGhost, 0 - NGhost, -1

        Byf_1D  ( 1 - NGhost: ntp_zones + NGhost) =
     1  Byf_amr ( 1 - NGhost: ntp_zones + NGhost, iyy, izz, ii)

        CALL PRINT_VECT ( 1 - NGhost, ntp_zones + NGhost,
     1           Byf_1D ( 1 - NGhost: ntp_zones + NGhost), "R" )

      END DO

      WRITE ( 6, *)"-----------------------------------------------"

! ----------

      WRITE ( 6, *)"Bzf_amr"

      DO iyy = ntp_zones + NGhost, 1 - NGhost, -1

        Bzf_1D  ( 1 - NGhost: ntp_zones + NGhost) =
     1  Bzf_amr ( 1 - NGhost: ntp_zones + NGhost, iyy, izz, ii)

        CALL PRINT_VECT ( 1 - NGhost, ntp_zones + NGhost,
     1           Bzf_1D ( 1 - NGhost: ntp_zones + NGhost), "I")

      END DO

      WRITE ( 6, *)"-----------------------------------------------"

! ------------------------------------------------------------

      END IF ! End of "IF ( isector == iactive_FromProc ( ii, 1)) THEN".

! ----------------------------------------------------------------------

      END IF ! End of "IF ( me == isector ) THEN".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE PRINT_VARIABLES_ONE_PATCH_PER_SECTOR


