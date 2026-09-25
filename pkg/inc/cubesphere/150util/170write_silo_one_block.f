






      SUBROUTINE WRITE_SILO_ONE_BLOCK ( filenum, time, dtcur )

! ----------------------------------------------------------------------
! Written By: Sethupathy S.

! Tests the silo mesh write out for one block.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER filenum

      REAL time, dtcur

! ----------------------------------------------------------------------

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ii

      INTEGER block_index              ! global number for the patch.

      INTEGER itri_patch ( 4)          ! extracted from itri_patch_amr

! ----------------------------------------------------------------------

! Synchronize on entry.
!     CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------

! Global patch location, counted first along the processor direction,
      block_index = me + ( ii - 1) * nimages

      itri_patch ( :) = itri_patch_amr ( :, ii)

      vert_coords_patch( :, :, :, :) = vert_coords_amr ( :, :, :, :, ii)

      imask_active_patch ( :, :, :) = imask_active_amr ( :, :, :, ii)

! Call the tester subroutine for only one patch.
      IF ( ii == 1 ) THEN
#if 0
      CALL WRITE_ONE_BLOCK ( block_index,
     1                       nr_zones,  NGhost,
     1                       ntp_zones, NGhost,
     1                       itri_patch,
     1                       vert_coords_patch )
#endif
      END IF

      END DO     ! End of "DO ii = 1, npatch_per_proc".

      END SUBROUTINE WRITE_SILO_ONE_BLOCK



