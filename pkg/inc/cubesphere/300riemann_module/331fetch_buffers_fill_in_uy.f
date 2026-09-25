! contained in RIEMANN_3D_CS_AMR_SUPPORT

















      SUBROUTINE FETCH_BUFFERS_FILL_IN_UY ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the xy-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! If the neighbouring patches lie within the sector:
! For the patch to its up,   it does one-sided get of "U_snd_d" from that patch.
! For the patch to its down, it does one-sided get of "U_snd_u" from that patch.

! If the neighbouring patches lie outside of the sector:
! For the patch to its up,   it does one-sided get of "U_snd_r" from that patch.
! For the patch to its down, it does one-sided get of "U_snd_d" from that patch.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, ixc, iyc, izc,
     1        ixn, iyn, izn, ixxp1, iyyp1, iyy0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_yface ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost)

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Evaluate number of real variables in the "U_recv_, U_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 ) * ( NGhost + 1 )
     1                * ( n_cc_components + 4 )

! ----------------------------------------------------------------------

! STEP 1.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the up sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! up sibling patch.

! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_d ( :, :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "U_recev_u"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_r ( :, :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "U_recev_l"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "U_recev_l"
! is unpacked into "U_recev_u" with the flipping.

! STEP 1.d) Now that the data in either of the two above cases is in
! "U_recev_u we can unpack it into "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! up sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access up sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for up sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 4)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 4)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 4)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc
! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which By locations will undergo averaging.
! At all y-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 0 - NGhost, ntp_zones + NGhost

        DO ixx = 1 - NGhost, ntp_zones + NGhost
        iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)
        iyy0 = MAX0 ( iyy, 1 - NGhost)

        imask_yface ( ixx, iyy) = imask_active_amr( ixx, iyy0, izz, ii)
     1                          * imask_active_amr( ixx, iyyp1, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_d ( :, :, :, :, yr_patch_in_proc)", and save in "U_recev_u".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "U_snd_d ( :, :, :, :, yr_patch_in_proc)" to "U_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_u, num_of_reals_xy, MPI_REAL8,    ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,  ! remote data on yr.
     1               MPI_REAL8, Data_U_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_r ( :, :, :, :, yr_patch_in_proc)". This is saved in "U_recev_l".
! Index flipping is done so that "U_recev_l" is unpacked into "U_recev_u"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "U_snd_r ( :, :, :, :, yr_patch_in_proc)" to "U_recev_l".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_XR_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_l, num_of_reals_xy, MPI_REAL8,    ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,  ! remote data on yr.
     1               MPI_REAL8, Data_U_XR_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_XR_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "U_recev_l" to "U_recev_u".
! Watch out, because "U_recev_l" does not have data in normal ordering.

! ----------
! Zone-centered data "U_amr, entropygr_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyc = 1 + NGhost        ! begin index + end index.

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.
        ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

        U_recev_u ( ixx, iyy, izz, 1: n_cc_components + 1) =
     1  U_recev_l ( ixn, iyn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-face-centered data "Bxf_amr".
!  my_sector x-face == (+ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        iyc = 1 + NGhost        ! begin index + end index.

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.
        ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

        U_recev_u ( ixx, iyy, izz, n_cc_components + 2) =
     1  U_recev_l ( ixn, iyn, izn, n_cc_components + 3)

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-face-centered data "Byf_amr".
!  my_sector y-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

          iyc = 0 + NGhost        ! begin index + end index.

          izn = izz
          iyn = ixx               ! yr_sector - y  is  my_sector - x.
          ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

          U_recev_u ( ixx, iyy, izz, n_cc_components + 3) =
     1  - U_recev_l ( ixn, iyn, izn, n_cc_components + 2)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-face-centered data "Bzf_amr".
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyc = 1 + NGhost        ! begin index + end index.

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.
        ixn = iyc - iyy         ! yr_sector upper-x  is  my_sector lower-y.

        U_recev_u ( ixx, iyy, izz, n_cc_components + 4) =
     1  U_recev_l ( ixn, iyn, izn, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 1.d) Unpack data from "U_recev_u" to
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 + ntp_zones, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

! ----------

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

        U_amr     ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  U_recev_u ( ixx, iyc, izz, 1: n_cc_components)

        entropygr_amr ( ixx, iyy, izz, ii) =
     1  U_recev_u     ( ixx, iyc, izz, n_cc_components + 1)

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Bxf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 + ntp_zones, ntp_zones + NGhost
      DO ixx = 0 - NGhost, ntp_zones  + NGhost

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

        Bxf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_u ( ixx, iyc, izz, n_cc_components + 2)

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Byf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 + ntp_zones, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

! ----------
! Exclude any copying or averaging at y- faces abutting tri-corners.
        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN

        IF ( ( iyy .EQ. ntp_zones) .AND.    ! Active Coinciding face.
     1       ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Byf_amr   ( ixx, iyy, izz, ii) =
     1  ( Byf_amr   ( ixx, iyy, izz, ii)
     1  + U_recev_u ( ixx, iyc, izz, n_cc_components + 3)) * 0.5

        ELSE                                ! Ghost face.

          Byf_amr   ( ixx, iyy, izz, ii) =
     1    U_recev_u ( ixx, iyc, izz, n_cc_components + 3)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Bxf_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 + ntp_zones, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        iyc = iyy - ntp_zones     ! location in "U_recev_u".

        Bzf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_u ( ixx, iyc, izz, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 + ntp_zones, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc" loop.

! ----------------------------------------------------------------------

! STEP 2.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the down sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! down sibling patch.

! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_u ( :, :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "U_recev_d"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_d ( :, :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "U_recev_u"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "U_recev_u"
! is unpacked into "U_recev_d" with the flipping.

! STEP 2.d) Now that the data in either of the two above cases is in
! "U_recev_d", we can unpack it into "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! down sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector)" will allow us to access down  sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for down sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 5)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 5)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 5)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------
! Identify Tri-corners for this patch.

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

! ------------------------------
! Also make a 2D mask which tells us which By locations will undergo averaging.
! At all y-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 0 - NGhost, ntp_zones + NGhost

        DO ixx = 1 - NGhost, ntp_zones + NGhost
        iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)
        iyy0 = MAX0 ( iyy, 1 - NGhost)

        imask_yface ( ixx, iyy) = imask_active_amr( ixx, iyy0, izz, ii)
     1                          * imask_active_amr( ixx, iyyp1, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "U_snd_u ( :, :, :, :, yr_patch_in_proc)", and save in "U_recev_d".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "U_snd_u ( :, :, :, :, yr_patch_in_proc)" to "U_recev_d".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_YU_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_d, num_of_reals_xy, MPI_REAL8,     ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,   ! remote data on yr.
     1               MPI_REAL8, Data_U_YU_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_YU_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "U_snd_d ( :, :, :, :, yr_patch_in_proc)". This is saved in "U_recev_u".
! Index flipping is done so that "U_recev_u" is unpacked into "U_recev_d"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "U_snd_d ( :, :, :, :, yr_patch_in_proc)" to "U_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_U_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( U_recev_u, num_of_reals_xy, MPI_REAL8,     ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,   ! remote data on yr.
     1               MPI_REAL8, Data_U_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_U_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "U_recev_u" to "U_recev_d".
! Watch out, because "U_recev_u" does not have data in normal ordering.

! ----------
! Zone-centered data "U_amr, entropygr_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        ixc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.
        iyc = 1 + NGhost                       ! begin index + end index.

        izn = izz
        ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
        iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

        U_recev_d ( ixx, iyy, izz, 1: n_cc_components + 1) =
     1  U_recev_u ( ixn, iyn, izn, 1: n_cc_components + 1)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-face-centered data "Bxf_amr".
!  my_sector x-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

         ixc = 0 - NGhost + ntp_zones + NGhost  ! begin index + end index.
         iyc = 1 + NGhost                       ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
         iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

         U_recev_d ( ixx, iyy, izz, n_cc_components + 2) =
     1 - U_recev_u ( ixn, iyn, izn, n_cc_components + 2)

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-face-centered data "Byf_amr".
!  my_sector y-face == (-ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

         ixc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.
         iyc = 0 + NGhost                       ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
         iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

         U_recev_d ( ixx, iyy, izz, n_cc_components + 3) =
     1 - U_recev_u ( ixn, iyn, izn, n_cc_components + 3)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0, NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-face-centered data "Bzf_amr".
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1, NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        ixc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.
        iyc = 1 + NGhost                       ! begin index + end index.

        izn = izz
        ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.
        iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

        U_recev_d ( ixx, iyy, izz, n_cc_components + 4) =
     1  U_recev_u ( ixn, iyn, izn, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1, NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 2.d) Unpack data from "U_recev_d" to
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, 0
      DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------

        iyc = iyy + NGhost        ! location in "U_recev_d".

        U_amr     ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  U_recev_d ( ixx, iyc, izz, 1: n_cc_components)

        entropygr_amr ( ixx, iyy, izz, ii) =
     1  U_recev_d     ( ixx, iyc, izz, n_cc_components + 1)

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, 0".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_d" to "Bxf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, 0
      DO ixx = 0 - NGhost, ntp_zones  + NGhost

        iyc = iyy + NGhost        ! location in "U_recev_d".

        Bxf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_d ( ixx, iyc, izz, n_cc_components + 2)

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, 0".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Byf_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, 0
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        iyc = iyy + NGhost        ! location in "U_recev_d".

! ----------
! Exclude any copying or averaging at y- faces abutting tri-corners.
        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN

        IF ( ( iyy .EQ. 0) .AND.         ! Active Coinciding face.
     1       ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          Byf_amr   ( ixx, iyy, izz, ii) =
     1  ( Byf_amr   ( ixx, iyy, izz, ii)
     1  + U_recev_d ( ixx, iyc, izz, n_cc_components + 3)) * 0.5

        ELSE                             ! Ghost face.

          Byf_amr   ( ixx, iyy, izz, ii) =
     1    U_recev_d ( ixx, iyc, izz, n_cc_components + 3)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 0 - NGhost, 0".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "U_recev_u" to "Bxf_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, 0
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        iyc = iyy + NGhost        ! location in "U_recev_d".

        Bzf_amr   ( ixx, iyy, izz, ii) =
     1  U_recev_d ( ixx, iyc, izz, n_cc_components + 4)

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, 0".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc" loop.

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_UY

