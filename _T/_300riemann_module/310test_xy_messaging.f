







      SUBROUTINE TEST_XY_MESSAGING

! ----------------------------------------------------------------------
! Written By: Sethupathy S.
!
! This test subroutine so far prints-out the Fluid variables
! assuming 1 patch/sector. It can be modified to adopt for different cases
! as well.
!
! STEPS:
! ------
! - Initialize the active zones and faces in "U_amr, entropygr_amr,
!   Bxf_amr, Byf_amr, Bzf_amr".
! - Prints out the values before messaging
! - Do messaging
! - Prints out the values after messaging
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc

      INTEGER ixx, iyy, i, ii

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Fill in the values in U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      U_amr ( :, :, :, :, :) = 0.0
      entropygr_amr ( :, :, :, :) = 0.0

      Bxf_amr ( :, :, :, :) = 0.0
      Byf_amr ( :, :, :, :) = 0.0
      Bzf_amr ( :, :, :, :) = 0.0

      DO iyy = 0, ntp_zones
        DO ixx = 0, ntp_zones

          IF( ixx > 0 .AND. iyy > 0 ) THEN
            U_amr ( ixx, iyy, :,   :, :) = iyy + 1
            entropygr_amr ( ixx, iyy, :,  :) = iyy + 1

            Bzf_amr ( ixx, iyy, :,  :) = iyy + 1
          END IF

          IF( iyy > 0 ) THEN
            Bxf_amr ( ixx, iyy, :,  :) = iyy + 1
          END IF

          IF( ixx > 0 ) THEN
            Byf_amr ( ixx, iyy, :,  :) = iyy + 1
          END IF

        END DO
      END DO

! ----------------------------------------------------------------------

      i = 1
      CALL PRINT_VARIABLES_ONE_PATCH_PER_SECTOR ( i)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CALL WINCREATE_MESSAGING ()

! ----------------------------------------------------------------------
! Enable/disable the function calls as needed for testing

      CALL LOAD_LOCAL_BUFFERS_UXY ()

      CALL FETCH_BUFFERS_FILL_IN_UX ()

!      CALL LOAD_LOCAL_BUFFERS_UXY ()

!      CALL FETCH_BUFFERS_FILL_IN_UY ()

!      CALL LOAD_LOCAL_BUFFERS_UXY ()

!      CALL FETCH_BUFFERS_FILL_IN_UY ()

! ----------------------------------------------------------------------

      CALL WINFREE_MESSAGING ()

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CALL PRINT_VARIABLES_ONE_PATCH_PER_SECTOR ( i)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE TEST_XY_MESSAGING

