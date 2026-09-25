! contained in RIEMANN_3D_CS_AMR_SUPPORT




      SUBROUTINE BUILD_BASE_PATCHES ()

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Builds "iactive_FromSector, isibling_number_t, isibling_sector_t, isibling_proc_t"
! for base level patches.
! Every processor will build this data locally so that it knows how to access
! data that might be on other sectors.

! When setting the number of processors and "Max_patches_per_proc", pay
! special attention to the logical flow of this subroutine.
! The check should be :-
! "6 * Max_base_patches_per_sector == nimages * npatch_per_proc".

! For load balance, "nimages" must be a multiple of 6; the number of sectors.

! So, once you know the total number of MPI processors "nimages",
! make sure that "Max_patches_per_proc .GE. npatch_per_proc".

! Once the above information is built, we set up "iactive_FromProc".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me,
     1        ixx, iyy, izz, isector, ii, jj, npatch_per_proc,
     1        ixxm1, ixxp1, iyym1, iyyp1, izzm1, izzp1,
     1        isector_rt, isector_lf, isector_up, isector_dn,
     1        ixxn, iyyn, iin, ipatch

      INTEGER neighbor_sector_list ( 6, 4)

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

! ----------

! We will have a total of "6 * Max_base_patches_per_sector" patches on the base mesh.
! If that does not divide evenly into "nimages", we have a load imbalance
! on the base level mesh. This should be treated as unacceptable.

! Please note that this check also forces "nimages" to be a multiple of six!

      IF ( ( MOD ( 6 * Max_base_patches_per_sector, nimages) .NE. 0)
     1   .AND. ( me .EQ. 1) ) THEN
      WRITE ( 6, *)"Stopping in BUILD_BASE_PATCHES : patches unbalanced"
      WRITE ( 6, *)"6 * Max_base_patches_per_sector, nimages = ",
     1             6 * Max_base_patches_per_sector, nimages
      CALL MPI_ABORT( MPI_COMM_WORLD, ERRORCODE, ierr)
      END IF

! ----------

! Find number of patches we can put per processor. It is bounded by
! "Max_patches_per_proc" so check for that.

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

      IF ( ( npatch_per_proc .GT. Max_patches_per_proc) .AND.
     1     ( me .EQ. 1)) THEN
      WRITE ( 6, *)"Stopping in BUILD_BASE_PATCHES : not enough patches"
      WRITE ( 6, *)"npatch_per_proc, Max_patches_per_proc = ",
     1              npatch_per_proc, Max_patches_per_proc
      CALL MPI_ABORT( MPI_COMM_WORLD, ERRORCODE, ierr)
      END IF

! ----------------------------------------------------------------------
! Set up "neighbor_sector_list ( 6, 4)". First index is parent sector.
! Second index is 1 (right), 2 (left), 3 (up), 4 (down) sectors relative
! to parent sector. z-direction does not need this mapping.

       neighbor_sector_list ( 1, 1) = 3
       neighbor_sector_list ( 1, 2) = 6
       neighbor_sector_list ( 1, 3) = 2
       neighbor_sector_list ( 1, 4) = 5

       neighbor_sector_list ( 2, 1) = 1
       neighbor_sector_list ( 2, 2) = 4
       neighbor_sector_list ( 2, 3) = 3
       neighbor_sector_list ( 2, 4) = 6

       neighbor_sector_list ( 3, 1) = 2
       neighbor_sector_list ( 3, 2) = 5
       neighbor_sector_list ( 3, 3) = 1
       neighbor_sector_list ( 3, 4) = 4

       neighbor_sector_list ( 4, 1) = 5
       neighbor_sector_list ( 4, 2) = 2
       neighbor_sector_list ( 4, 3) = 6
       neighbor_sector_list ( 4, 4) = 3

       neighbor_sector_list ( 5, 1) = 6
       neighbor_sector_list ( 5, 2) = 3
       neighbor_sector_list ( 5, 3) = 4
       neighbor_sector_list ( 5, 4) = 1

       neighbor_sector_list ( 6, 1) = 4
       neighbor_sector_list ( 6, 2) = 1
       neighbor_sector_list ( 6, 3) = 5
       neighbor_sector_list ( 6, 4) = 2

! ----------------------------------------------------------------------

! The first "Max_base_patches_per_sector" on each sector will be active and their
! sector number (first index) will point to themselves.

      DO isector = 1, 6

      iactive_FromSector ( isector, 1: Max_base_patches_per_sector, 1) =
     1                     isector

      END DO

! ----------

! The first "Max_base_patches_per_sector" on each sector will have their patch number
! (second index) be their ordinal number within the sector.

! The processor on which the patch lands is given by the third index.
! Each sector has "Max_Base_patches_per_sector" number of patches. But the patches
! in a single sector can be distributed across multiple processors.
! The data will be contiguous in order to facilitate efficient messaging.

      ii = 0
      jj = Max_Base_patches_per_sector

      DO izz = 1, nr_patches  ! z-directional patches are slowest moving.
        DO iyy = 1, ntp_patches  ! y-directional patches are intermediate.
          DO ixx = 1, ntp_patches  ! x-directional patches are contiguous.

          ii = ii + 1

          iactive_FromSector ( :, ii, 2) = ii

          END DO  ! End of "DO ixx = 1, ntp_patches" loop.
        END DO  ! End of "DO iyy = 1, ntp_patches" loop.
      END DO  ! End of "DO izz = 1, nr_patches" loop.

! ----------

      IF ( ( ii .NE. Max_base_patches_per_sector) .AND.
     1     ( me .EQ. 1)) THEN
      WRITE ( 6, *)"Stopping in BUILD_BASE_PATCHES :"
      WRITE ( 6, *)"ii .NE. Max_base_patches_per_sector = ",
     1              ii, Max_base_patches_per_sector
      CALL MPI_ABORT( MPI_COMM_WORLD, ERRORCODE, ierr)
      END IF

! ----------

! The processor on which the patch lands is given by the third index.
! Each sector has "Max_Base_patches_per_sector" number of patches. But the patches
! in a single sector can be distributed across multiple processors.
! The data will be contiguous in order to facilitate efficient messaging.

      DO isector = 1, 6

      ii = 0
      jj = Max_Base_patches_per_sector

      DO izz = 1, nr_patches  ! z-directional patches are slowest moving.
        DO iyy = 1, ntp_patches  ! y-directional patches are intermediate.
          DO ixx = 1, ntp_patches  ! x-directional patches are contiguous.

          ii = ii + 1

          iactive_FromSector ( isector, ii, 3) =
     1    1 + (ii + (isector - 1) * jj - 1) / npatch_per_proc

          END DO  ! End of "DO ixx = 1, ntp_patches" loop.
        END DO  ! End of "DO iyy = 1, ntp_patches" loop.
      END DO  ! End of "DO izz = 1, nr_patches" loop.

      END DO  ! End of "DO isector = 1, 6" loop.

! ----------------------------------------------------------------------
! Now start building "isibling_number_t, isibling_sector_t, isibling_proc_t".
! There are 7 neighbors because the first neighbor is to the self.

! First we fill in the self identification.

      isibling_sector_t ( :, 1: Max_base_patches_per_sector, 1) =
     1       iactive_FromSector ( :, 1: Max_base_patches_per_sector, 1)

      isibling_number_t ( :, 1: Max_base_patches_per_sector, 1) =
     1       iactive_FromSector ( :, 1: Max_base_patches_per_sector, 2)

      isibling_proc_t ( :, 1: Max_base_patches_per_sector, 1) =
     1       iactive_FromSector ( :, 1: Max_base_patches_per_sector, 3)

! ----------
! Loop over sectors and find their neighbors.

      DO isector = 1, 6

      ii = 0
      jj = Max_Base_patches_per_sector

      DO izz = 1, nr_patches  ! z-directional patches are slowest moving.
      izzp1 = izz + 1
      izzm1 = izz - 1

        DO iyy = 1, ntp_patches  ! y-directional patches are intermediate.
        iyyp1 = iyy + 1
        iyym1 = iyy - 1

          DO ixx = 1, ntp_patches  ! x-directional patches are contiguous.
          ixxp1 = ixx + 1
          ixxm1 = ixx - 1

          ii = ii + 1

! ----------------------------------------
! Find right neighbor; could be in-sector (easy) or out-of-sector (trickier).

          IF ( ixxp1 .LE. ntp_patches) THEN
! Right neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 2) = isector

            isibling_number_t ( isector, ii, 2) = ii + 1

            isibling_proc_t ( isector, ii, 2) =
     1      1 + (ii + (isector - 1) * jj) / npatch_per_proc

! ----------

          ELSE
! Right neighbor is out-of-sector.

! ----------

! This is the right neighboring sector.
            isector_rt = neighbor_sector_list ( isector, 1)

! On right neighboring sector, x and y get flipped; see figure in notes.
! Find "iin" as the right sibling's patch number on sector "isector_rt".
            ixxn = iyy
            iyyn = ntp_patches
            iin = ixxn + (iyyn - 1) * ntp_patches
     1          + (izz - 1) * ntp_patches**2

            isibling_sector_t ( isector, ii, 2) = isector_rt

            isibling_number_t ( isector, ii, 2) =
     1      iactive_FromSector ( isector_rt, iin, 2)

            isibling_proc_t ( isector, ii, 2) =
     1      iactive_FromSector ( isector_rt, iin, 3)

! ----------

          END IF  ! End of "IF ( ixxp1 .LE. ntp_patches) THEN".

! ----------------------------------------
! Find left neighbor; could be in-sector (easy) or out-of-sector (trickier).

          IF ( ixxm1 .GE. 1) THEN
! Left neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 3) = isector

            isibling_number_t ( isector, ii, 3) = ii - 1

            isibling_proc_t ( isector, ii, 3) =
     1      1 + (ii + (isector - 1) * jj - 2) / npatch_per_proc

! ----------

          ELSE
! Left neighbor is out-of-sector.

! ----------

! This is the left neighboring sector.
            isector_lf = neighbor_sector_list ( isector, 2)

! On left neighboring sector, x and y are not flipped but y-sequence
! is reversed; see figure in notes.
! Find "iin" as the left sibling's patch number on sector "isector_lf".
            ixxn = 1
            iyyn = ntp_patches - iyy + 1
            iin = ixxn + (iyyn - 1) * ntp_patches
     1          + (izz - 1) * ntp_patches**2

            isibling_sector_t ( isector, ii, 3) = isector_lf

            isibling_number_t ( isector, ii, 3) =
     1      iactive_FromSector ( isector_lf, iin, 2)

            isibling_proc_t ( isector, ii, 3) =
     1      iactive_FromSector ( isector_lf, iin, 3)

! ----------

          END IF  ! End of "IF ( ixxm1 .GE. 1) THEN".

! ----------------------------------------
! Find upper neighbor; could be in-sector (easy) or out-of-sector (trickier).

          IF ( iyyp1 .LE. ntp_patches) THEN
! Upper neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 4) = isector

            isibling_number_t ( isector, ii, 4) = ii + ntp_patches

            isibling_proc_t ( isector, ii, 4) =
     1      1 + (ii + ntp_patches + (isector - 1) * jj - 1)
     1          / npatch_per_proc

! ----------

          ELSE
! Upper neighbor is out-of-sector.

! ----------

! This is the upper neighboring sector.
            isector_up = neighbor_sector_list ( isector, 3)

! On the upper neighboring sector, x and y get flipped; see figure in notes.
! Find "iin" as the upper sibling's patch number on sector "isector_up".
            ixxn = ntp_patches
            iyyn = ixx
            iin = ixxn + (iyyn - 1) * ntp_patches
     1          + (izz - 1) * ntp_patches**2

            isibling_sector_t ( isector, ii, 4) = isector_up

            isibling_number_t ( isector, ii, 4) =
     1      iactive_FromSector ( isector_up, iin, 2)

            isibling_proc_t ( isector, ii, 4) =
     1      iactive_FromSector ( isector_up, iin, 3)

! ----------

          END IF  ! End of "IF ( iyyp1 .LE. ntp_patches) THEN".

! ----------------------------------------
! Find down neighbor; could be in-sector (easy) or out-of-sector (trickier).

          IF ( iyym1 .GE. 1) THEN
! Down neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 5) = isector

            isibling_number_t ( isector, ii, 5) = ii - ntp_patches

            isibling_proc_t ( isector, ii, 5) =
     1      1 + (ii - ntp_patches + (isector - 1) * jj - 1)
     1          / npatch_per_proc

! ----------

          ELSE
! Down neighbor is out-of-sector.

! ----------

! This is the down neighboring sector.
            isector_dn = neighbor_sector_list ( isector, 4)

! On the down neighboring sector, x and y are not flipped but x-sequence
! is reversed; see figure in notes.
! Find "iin" as the down sibling's patch number on sector "isector_dn".
            ixxn = ntp_patches - ixx + 1
            iyyn = 1
            iin = ixxn + (iyyn - 1) * ntp_patches
     1          + (izz - 1) * ntp_patches**2

            isibling_sector_t ( isector, ii, 5) = isector_dn

            isibling_number_t ( isector, ii, 5) =
     1      iactive_FromSector ( isector_dn, iin, 2)

            isibling_proc_t ( isector, ii, 5) =
     1      iactive_FromSector ( isector_dn, iin, 3)

! ----------

          END IF  ! End of "IF ( iyym1 .GE. 1) THEN".

! ----------------------------------------
! Find top neighbor; could be in-sector (easy) or non-existent.

          IF ( izzp1 .LE. nr_patches) THEN
! Top neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 6) = isector

            isibling_number_t ( isector, ii, 6) = ii + ntp_patches**2

            isibling_proc_t ( isector, ii, 6) =
     1      1 + (ii + ntp_patches**2 + (isector - 1) * jj - 1)
     1          / npatch_per_proc

! ----------

          ELSE
! Top neighbor is non-existent (will be filled in with boundary condition routine).

! ----------

            isibling_sector_t ( isector, ii, 6) = -1

            isibling_number_t ( isector, ii, 6) = -1

            isibling_proc_t ( isector, ii, 6) = -1

! ----------

          END IF  ! End of "IF ( izzp1 .LE. nr_patches) THEN".

! ----------------------------------------
! Find bottom neighbor; could be in-sector (easy) or non-existent.

          IF ( izzm1 .GE. 1) THEN
! Bottom neighbor is in-sector.

! ----------

            isibling_sector_t ( isector, ii, 7) = isector

            isibling_number_t ( isector, ii, 7) = ii - ntp_patches**2

            isibling_proc_t ( isector, ii, 7) =
     1      1 + (ii - ntp_patches**2 + (isector - 1) * jj - 1)
     1          / npatch_per_proc

! ----------

          ELSE
! Bottom neighbor is non-existent (will be filled in with boundary condition routine).

! ----------

            isibling_sector_t ( isector, ii, 7) = -1

            isibling_number_t ( isector, ii, 7) = -1

            isibling_proc_t ( isector, ii, 7) = -1

! ----------

          END IF  ! End of "IF ( izzm1 .GE. 1) THEN".

! ----------------------------------------

          END DO  ! End of "DO ixx = 1, ntp_patches" loop.

        END DO  ! End of "DO iyy = 1, ntp_patches" loop.

      END DO  ! End of "DO izz = 1, nr_patches" loop.

      END DO  ! End of "DO isector = 1, 6" loop.

! ----------------------------------------------------------------------
! Now we set up "iactive_FromProc". Here we are looping over the list of
! patches that reside on each processor. These could be on any sector and
! have any logical patch number. At least for the base mesh, this array
! gives us the (sector#, logical patch#).

      iactive_FromProc ( :, :) = -1

      DO ii = 1, npatch_per_proc

      jj = ii + (me - 1) * npatch_per_proc

      isector = 1 + ( ( jj - 1) / Max_base_patches_per_sector)
      ipatch = jj - (isector - 1) * Max_base_patches_per_sector

      iactive_FromProc ( ii, 1) = isector  ! Which sector is this patch in?
      iactive_FromProc ( ii, 2) = ipatch  ! logical patch number in "isector".

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE BUILD_BASE_PATCHES
