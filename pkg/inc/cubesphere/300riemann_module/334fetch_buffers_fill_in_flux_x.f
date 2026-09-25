! contained in RIEMANN_3D_CS_AMR_SUPPORT











      SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_X ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the xy-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! If the neighbouring patches lie within the sector:
! For the patch to its right, it does one-sided get of "Flux_snd_l" from that patch.
! For the patch to its left,  it does one-sided get of "Flux_snd_r" from that patch.

! If the neighbouring patches lie outside of the sector:
! For the patch to its right, it does one-sided get of "Flux_snd_u" from that patch.
! For the patch to its left,  it does one-sided get of "Flux_snd_l" from that patch.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, iyc,
     1        ixn, iyn, izn, ixxp1, iyyp1, ixx0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_xface ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost)

! Integer constants, that specify the beginning and end locations of
! fluctuation terms in the Flux-communication arrays.

      INTEGER, PARAMETER :: n_cc_comp = NVAR + NFLUID,

     1                      ip1 = 3 + 1,
     1                      ip2 = 3 + n_cc_comp,

     1                      im1 = 3 + n_cc_comp + 1,
     1                      im2 = 3 + n_cc_comp + n_cc_comp

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Evaluate number of real variables in the "Flux_recv_, Flux_snd_".
! In this cubesphere code, the mesh is symmetric w.r.t. "x" and "y" axis.
! Therefore, the message data size is same for "x" and "y" axes and it is
! different for "z" axis.

      num_of_reals_xy = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( nr_zones  + 2 * NGhost + 1 )
     1                * ( 3 + 2*n_cc_components )

      num_of_reals_z  = ( ntp_zones + 2 * NGhost + 1 )
     1                * ( ntp_zones + 2 * NGhost + 1 )
     1                * ( 3 + 2*n_cc_components )

! ----------------------------------------------------------------------

! STEP 1.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the right sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! right sibling patch.

! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_l ( :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "Flux_recev_r"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_u ( :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "Flux_recev_d"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "Flux_recev_d"
! is unpacked into "Flux_recev_r" with the flipping.

! STEP 1.d) Now that the data in either of the two above cases is in
! "Flux_recev_r", we can unpack it into "x-Fluctuations, EntFlux_x_amr,
! Efy_amr, Efz_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! right sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access right sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for right sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 2)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 2)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 2)

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
! Also make a 2D mask which tells us which Bx locations will undergo averaging.
! At all x-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 1 - NGhost, ntp_zones + NGhost

        DO ixx = 0 - NGhost, ntp_zones + NGhost
        ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)
        ixx0 = MAX0 ( ixx, 1 - NGhost)

        imask_xface ( ixx, iyy) = imask_active_amr( ixx0, iyy, izz, ii)
     1                          * imask_active_amr( ixxp1, iyy, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_l ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_r".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_l ( :, :, :, yr_patch_in_proc)" to "Flux_recev_r".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_XL_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_r, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_XL_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_XL_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_u ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_d".
! Index flipping is done so that "Flux_recev_d" is unpacked into "Flux_recev_r"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "Flux_snd_u ( :, :, :, yr_patch_in_proc)" to "Flux_recev_d".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YU_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_d, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YU_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YU_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "Flux_recev_d" to "Flux_recev_r".
! Watch out, because "Flux_recev_d" does not have data in normal ordering.

!  my_sector x-face == (-ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

         izn = izz
         ixn = iyy              ! yr_sector - x  is  my_sector - y.

         Flux_recev_r ( iyy, izz, 1: 1) =
     1 - Flux_recev_d ( ixn, izn, 1: 1)

         Flux_recev_r ( iyy, izz, ip1: ip2) =
     1   Flux_recev_d ( ixn, izn, im1: im2)       ! Fluct_x_Pl = + Fluct_y_Mi

         Flux_recev_r ( iyy, izz, im1: im2) =
     1   Flux_recev_d ( ixn, izn, ip1: ip2)       ! Fluct_x_Mi = + Fluct_y_Pl

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-edge-centered data "Efy_amr".
!  my_sector y-edge == (+ve) yr_sector x-edge.
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

        izn = izz
        ixn = iyy               ! yr_sector - x  is  my_sector - y.

        Flux_recev_r ( iyy, izz, 2) =    ! YZ: Y at +2
     1  Flux_recev_d ( ixn, izn, 3)      ! ZX: X at +3

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost

        izn = izz
        ixn = iyy               ! yr_sector - x  is  my_sector - y.

        Flux_recev_r ( iyy, izz, 3) =    ! YZ: Z at +3
     1  Flux_recev_d ( ixn, izn, 2)      ! ZX: Z at +2

      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 1.d) Unpack data from "Flux_recev_r" to
! "x-Fluctuations, EntFlux_x_amr, Efy_amr, Efz_amr".

! --------------------
! Unpack data from "Flux_recev_r" to "x-Fluctuations, EntFlux_x_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

        ixx = ntp_zones

! ----------

        IF ( imask_xface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at x- faces abutting tri-corners.

        IF ( ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.  ! Active Coinciding face.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          EntFlux_x_amr( ixx, iyy, izz, ii) =
     1  ( EntFlux_x_amr( ixx, iyy, izz, ii)
     1  + Flux_recev_r (      iyy, izz, 1)) * 0.5

          Fluct_x_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_x_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_r (      iyy, izz, ip1: ip2)) * 0.5

          Fluct_x_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_x_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_r (      iyy, izz, im1: im2)) * 0.5

        ELSE

          EntFlux_x_amr( ixx, iyy, izz, ii) =
     1    Flux_recev_r (      iyy, izz, 1)

          Fluct_x_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_r (      iyy, izz, ip1: ip2)

          Fluct_x_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_r (      iyy, izz, im1: im2)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_r" to "Efy_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

        ixx = ntp_zones

! ----------

        IF ( imask_xface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at x- faces abutting tri-corners.

        IF ( ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( izz .GE. 0) .AND. ( izz .LE. nr_zones )) THEN

          Efy_amr ( ixx, iyy, izz, ii) =
     1  ( Efy_amr ( ixx, iyy, izz, ii)
     1  + Flux_recev_r ( iyy, izz, 2)) * 0.5

        ELSE

          Efy_amr ( ixx, iyy, izz, ii) =
     1    Flux_recev_r ( iyy, izz, 2)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_r" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost

        ixx = ntp_zones

! ----------

        IF ( ( iyy .GE. 0) .AND. ( iyy .LE. ntp_zones) .AND.   ! Active Coinciding edge.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          IF ( ( ( iyy .EQ. 0 ) .AND. ( itri_patch_RD .EQ. 1)) .OR.
     1         ( ( iyy .EQ. ntp_zones) .AND.
     1                                ( itri_patch_RU .EQ. 1)) ) THEN

            Efz_amr ( ixx, iyy, izz, ii) =   ! accumulate tri-corner Efz; * 1/3 later.
     1      Efz_amr ( ixx, iyy, izz, ii)
     1    + Flux_recev_r ( iyy, izz, 3)

          ELSE

            Efz_amr ( ixx, iyy, izz, ii) =
     1    ( Efz_amr ( ixx, iyy, izz, ii)
     1    + Flux_recev_r ( iyy, izz, 3)) * 0.5

          END IF

        ELSE

          Efz_amr ( ixx, iyy, izz, ii) =
     1    Flux_recev_r ( iyy, izz, 3)

        END IF

! ----------

      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

! STEP 2.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the left sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! left sibling patch.

! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r ( :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "Flux_recev_l"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_l ( :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "Flux_recev_r"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "Flux_recev_r"
! is unpacked into "Flux_recev_l" with the flipping.

! STEP 2.d) Now that the data in either of the two above cases is in
! "Flux_recev_l", we can unpack it into "x--Fluctuations, EntFlux_x_amr,
! Efy_amr, Efz_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! left sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access left sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for left sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 3)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 3)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 3)

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
! Also make a 2D mask which tells us which Bx locations will undergo averaging.
! At all x-faces that do not abut the tri-corner, we can do averaging.

      izz = 1

      DO iyy = 1 - NGhost, ntp_zones + NGhost

        DO ixx = 0 - NGhost, ntp_zones + NGhost
        ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)
        ixx0 = MAX0 ( ixx, 1 - NGhost)

        imask_xface ( ixx, iyy) = imask_active_amr( ixx0, iyy, izz, ii)
     1                          * imask_active_amr( ixxp1, iyy, izz, ii)

        END DO

      END DO

! ------------------------------------------------------------
! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_l".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_r ( :, :, :, yr_patch_in_proc)" to "Flux_recev_l".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_XR_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_l, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_XR_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_XR_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_l ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_r".
! Index flipping is done so that "Flux_recev_r" is unpacked into "Flux_recev_l"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "Flux_snd_l ( :, :, :, yr_patch_in_proc)" to "Flux_recev_r".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_XL_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_r, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_XL_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_XL_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "Flux_recev_r" to "Flux_recev_l".
! Watch out, because "Flux_recev_r" does not have data in normal ordering.

!  my_sector x-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

         iyc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         izn = izz
         iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

         Flux_recev_l ( iyy, izz, 1: 1) =
     1 - Flux_recev_r ( iyn, izn, 1: 1)

         Flux_recev_l ( iyy, izz, ip1: ip2) =
     1   Flux_recev_r ( iyn, izn, im1: im2)       ! Fluct_x_Pl = + Fluct_x_Mi

         Flux_recev_l ( iyy, izz, im1: im2) =
     1   Flux_recev_r ( iyn, izn, ip1: ip2)       ! Fluct_x_Mi = + Fluct_x_Pl

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! y-edge-centered data "Efy_amr".
!  my_sector y-edge == (-ve) yr_sector y-edge.
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

         iyc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.
         izn = izz

         Flux_recev_l ( iyy, izz, 2) =    ! YZ: Y at +2
     1 - Flux_recev_r ( iyn, izn, 2)      ! YZ: Y at +2

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost

        iyc = 0 - NGhost + ntp_zones + NGhost  ! begin index + end index.

        izn = izz
        iyn = iyc - iyy     ! yr_sector lower-y  is  my_sector  upper-y.

        Flux_recev_l ( iyy, izz, 3) =    ! YZ: Z at +3
     1  Flux_recev_r ( iyn, izn, 3)      ! YZ: Z at +3

      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 2.d) Unpack data from "Flux_recev_l" to
! "x--Fluctuations, EntFlux_x_amr, Efy_amr, Efz_amr".

! --------------------
! Unpack data from "Flux_recev_l" to "x--Fluctuations, EntFlux_x_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

        ixx = 0

! ----------

        IF ( imask_xface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at x- faces abutting tri-corners.

        IF ( ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.  ! Active Coinciding face.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          EntFlux_x_amr( ixx, iyy, izz, ii) =
     1  ( EntFlux_x_amr( ixx, iyy, izz, ii)
     1  + Flux_recev_l (      iyy, izz, 1)) * 0.5

          Fluct_x_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_x_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_l (      iyy, izz, ip1: ip2)) * 0.5

          Fluct_x_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_x_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_l (      iyy, izz, im1: im2)) * 0.5

        ELSE

          EntFlux_x_amr( ixx, iyy, izz, ii) =
     1    Flux_recev_l (      iyy, izz, 1)

          Fluct_x_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_l (      iyy, izz, ip1: ip2)

          Fluct_x_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_l (      iyy, izz, im1: im2)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_l" to "Efy_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost

        ixx = 0

! ----------

        IF ( imask_xface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at x- faces abutting tri-corners.

        IF ( ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( izz .GE. 0) .AND. ( izz .LE. nr_zones )) THEN

          Efy_amr ( ixx, iyy, izz, ii) =
     1  ( Efy_amr ( ixx, iyy, izz, ii)
     1  + Flux_recev_l ( iyy, izz, 2)) * 0.5

        ELSE

          Efy_amr ( ixx, iyy, izz, ii) =
     1    Flux_recev_l ( iyy, izz, 2)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_l" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 0 - NGhost, ntp_zones + NGhost

        ixx = 0

! ----------

        IF ( ( iyy .GE. 0) .AND. ( iyy .LE. ntp_zones) .AND. ! Active Coinciding edge.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          IF ( ( ( iyy .EQ. 0 ) .AND. ( itri_patch_LD .EQ. 1)) .OR.
     1         ( ( iyy .EQ. ntp_zones) .AND.
     1                                ( itri_patch_LU .EQ. 1)) ) THEN

            Efz_amr ( ixx, iyy, izz, ii) = ! accumulate tri-corner Efz; * 1/3 later.
     1      Efz_amr ( ixx, iyy, izz, ii)
     1    + Flux_recev_l ( iyy, izz, 3)

          ELSE

            Efz_amr ( ixx, iyy, izz, ii) =
     1    ( Efz_amr ( ixx, iyy, izz, ii)
     1    + Flux_recev_l ( iyy, izz, 3)) * 0.5

          END IF

        ELSE

          Efz_amr ( ixx, iyy, izz, ii) =
     1    Flux_recev_l ( iyy, izz, 3)

        END IF

! ----------

      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_X
