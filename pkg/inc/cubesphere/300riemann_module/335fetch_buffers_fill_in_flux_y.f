! contained in RIEMANN_3D_CS_AMR_SUPPORT








      SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_Y ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the xy-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! If the neighbouring patches lie within the sector:
! For the patch to its up,   it does one-sided get of "Flux_snd_d" from that patch.
! For the patch to its down, it does one-sided get of "Flux_snd_u" from that patch.

! If the neighbouring patches lie outside of the sector:
! For the patch to its up,   it does one-sided get of "Flux_snd_r" from that patch.
! For the patch to its down, it does one-sided get of "Flux_snd_d" from that patch.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, ixc,
     1        ixn, iyn, izn, ixxp1, iyyp1, iyy0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_yface ( 0 - NGhost: ntp_zones + NGhost,
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
! is to be able to extract the up sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! up sibling patch.

! STEP 1.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d ( :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "Flux_recev_u"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r ( :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "Flux_recev_l"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "Flux_recev_l"
! is unpacked into "Flux_recev_u" with the flipping.

! STEP 1.d) Now that the data in either of the two above cases is in
! "Flux_recev_u", we can unpack it into "y-Fluctuations_amr, EntFlux_y_amr,
! Efz_amr, Efx_amr".

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
! of "Flux_snd_d ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_u".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_d ( :, :, :, yr_patch_in_proc)" to "Flux_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_u, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_l".
! Index flipping is done so that "Flux_recev_l" is unpacked into "Flux_recev_u"
! with the flipping.

      ELSE

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
! Unpack data from "Flux_recev_l" to "Flux_recev_u".
! Watch out, because "Flux_recev_l" does not have data in normal ordering.

!  my_sector y-face == (-ve) yr_sector x-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

         izn = izz
         iyn = ixx               ! yr_sector - y  is  my_sector - x.

         Flux_recev_u ( ixx, izz, 1: 1) =
     1 - Flux_recev_l ( iyn, izn, 1: 1)

         Flux_recev_u ( ixx, izz, ip1: ip2) =
     1   Flux_recev_l ( iyn, izn, im1: im2)        ! Fluct_y_Pl = + Fluct_x_Mi

         Flux_recev_u ( ixx, izz, im1: im2) =
     1   Flux_recev_l ( iyn, izn, ip1: ip2)        ! Fluct_y_Mi = + Fluct_x_Pl

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.

        Flux_recev_u ( ixx, izz, 2) =    ! ZX: Z at +2
     1  Flux_recev_l ( iyn, izn, 3)      ! YZ: Z at +3

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-edge-centered data "Efx_amr".
!  my_sector x-edge == (+ve) yr_sector y-edge.
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        izn = izz
        iyn = ixx               ! yr_sector - y  is  my_sector - x.

        Flux_recev_u ( ixx, izz, 3) =    ! ZX: X at +3
     1  Flux_recev_l ( iyn, izn, 2)      ! YZ: Y at +2

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 1.d) Unpack data from "Flux_recev_u" to
! "y-Fluctuations, EntFlux_y_amr, Efz_amr, Efx_amr".

! --------------------
! Unpack data from "Flux_recev_u" to "y-Fluctuations, EntFlux_y_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyy = ntp_zones

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at y- faces abutting tri-corners.

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding face.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          EntFlux_y_amr( ixx, iyy, izz, ii) =
     1  ( EntFlux_y_amr( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, 1)) * 0.5

          Fluct_y_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_y_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_u ( ixx,      izz, ip1: ip2)) * 0.5

          Fluct_y_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_y_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_u ( ixx,      izz, im1: im2)) * 0.5

        ELSE

          EntFlux_y_amr( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, 1)

          Fluct_y_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_u ( ixx,      izz, ip1: ip2)

          Fluct_y_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_u ( ixx,      izz, im1: im2)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_u" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        iyy = ntp_zones

! ----------

        IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) .AND.   ! Active Coinciding edge.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

! -----
! accumulate tri-corner Efz; * 1/3 later.

          IF ( ( ( ixx .EQ. 0) .AND. ( itri_patch_LU .EQ. 1)) .OR.
     1         ( ( ixx .EQ. ntp_zones) .AND.
     1                               ( itri_patch_RU .EQ. 1))) THEN

            Efz_amr      ( ixx, iyy, izz, ii) =
     1      Efz_amr      ( ixx, iyy, izz, ii)
     1    + Flux_recev_u ( ixx,      izz, 2)

! -----
! Exclude quad corner on "Efz" averaging for y-flux; this is handled separately.

          ELSE IF ( ( ixx .NE. 0) .AND. ( ixx .NE. ntp_zones)) THEN

            Efz_amr      ( ixx, iyy, izz, ii) =
     1    ( Efz_amr      ( ixx, iyy, izz, ii)
     1    + Flux_recev_u ( ixx,      izz, 2)) * 0.5

! -----

          END IF

! ----------

        ELSE

! ----------

          Efz_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, 2)

! ----------

        END IF  ! IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) ...

! ----------

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_u" to "Efx_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyy = ntp_zones

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at y- faces abutting tri-corners.

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.   ! Active Coinciding edge.
     1       ( izz .GE. 0) .AND. ( izz .LE. nr_zones )) THEN

          Efx_amr      ( ixx, iyy, izz, ii) =
     1  ( Efx_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, 3)) * 0.5

        ELSE

          Efx_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, 3)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

! STEP 2.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the down sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! down sibling patch.

! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_u ( :, :, :, yr_patch_in_proc)", just
! as we would do in a Cartesian mesh code. This is saved in "Flux_recev_d"
! on the processor that owns "(my_sector, my_patch_in_proc)".

! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d ( :, :, :, yr_patch_in_proc)". This is different from what
! we have in a Cartesian mesh code. This is saved in "Flux_recev_u"
! on the processor that owns "(my_sector, my_patch_in_proc)". However, note
! that massive flipping of indices has to be done. Therefore, "Flux_recev_u"
! is unpacked into "Flux_recev_d" with the flipping.

! STEP 2.d) Now that the data in either of the two above cases is in
! "Flux_recev_d", we can unpack it into "y-Fluctuations, EntFlux_y_amr,
! Efz_amr, Efx_amr".

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

! "(my_sector, my_patch_in_sector" will allow us to access down sibling
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
! of "Flux_snd_u ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_d".

      IF ( my_sector .EQ. yr_sector) THEN

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
! STEP 2.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_u".
! Index flipping is done so that "U_recev_u" is unpacked into "Flux_recev_d"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "Flux_snd_d ( :, :, :, yr_patch_in_proc)" to "Flux_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_u, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! Unpack data from "Flux_recev_u" to "Flux_recev_d".
! Watch out, because "Flux_recev_u" does not have data in normal ordering.

!  my_sector y-face == (-ve) yr_sector y-face.
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

         ixc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.

         Flux_recev_d ( ixx, izz, 1: 1) =
     1 - Flux_recev_u ( ixn, izn, 1: 1)

         Flux_recev_d ( ixx, izz, ip1: ip2) =
     1   Flux_recev_u ( ixn, izn, im1: im2)        ! Fluct_y_Pl = + Fluct_y_Mi

         Flux_recev_d ( ixx, izz, im1: im2) =
     1   Flux_recev_u ( ixn, izn, ip1: ip2)        ! Fluct_y_Mi = + Fluct_y_Pl

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! z-edge-centered data "Efz_amr".
      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

         ixc = 0 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.

         Flux_recev_d ( ixx, izz, 2) =    ! ZX: Z at +2
     1   Flux_recev_u ( ixn, izn, 2)      ! ZX: Z at +2

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ----------
! x-edge-centered data "Efx_amr".
!  my_sector x-edge == (-ve) yr_sector x-edge.
      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

         ixc = 1 - NGhost + ntp_zones + NGhost  ! begin index + end index.

         izn = izz
         ixn = ixc - ixx     ! yr_sector lower-x  is  my_sector  upper-x.

         Flux_recev_d ( ixx, izz, 3) =    ! ZX: X at +3
     1 - Flux_recev_u ( ixn, izn, 3)      ! ZX: X at +3

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 2.d) Unpack data from "Flux_recev_d" to
! "y-Fluctuations, EntFlux_y_amr, Efz_amr, Efx_amr".

! --------------------
! Unpack data from "Flux_recev_d" to "y-Fluctuations, EntFlux_y_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyy = 0

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at y- faces abutting tri-corners.

        IF ( ( ixx .GE. 1) .AND. ( iyy .LE. ntp_zones) .AND.   ! Active Coinciding face.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

          EntFlux_y_amr( ixx, iyy, izz, ii) =
     1  ( EntFlux_y_amr( ixx, iyy, izz, ii)
     1  + Flux_recev_d ( ixx,      izz, 1)) * 0.5

          Fluct_y_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_y_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_d ( ixx,      izz, ip1: ip2)) * 0.5

          Fluct_y_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_y_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_d ( ixx,      izz, im1: im2)) * 0.5

        ELSE

          EntFlux_y_amr( ixx, iyy, izz, ii) =
     1    Flux_recev_d ( ixx,      izz, 1)

          Fluct_y_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_d ( ixx,      izz, ip1: ip2)

          Fluct_y_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_d ( ixx,      izz, im1: im2)

        END IF

        END IF  ! Exclude any copying or averaging at y- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_d" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones + NGhost

        iyy = 0

! ----------

        IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) .AND.   ! Active Coinciding edge.
     1       ( izz .GE. 1) .AND. ( izz .LE. nr_zones )) THEN

! -----
! accumulate tri-corner Efz; * 1/3 later.

          IF ( ( ( ixx .EQ. 0) .AND. ( itri_patch_LD .EQ. 1)) .OR.
     1         ( ( ixx .EQ. ntp_zones) .AND.
     1                               ( itri_patch_RD .EQ. 1))) THEN

            Efz_amr      ( ixx, iyy, izz, ii) =
     1      Efz_amr      ( ixx, iyy, izz, ii)
     1    + Flux_recev_d ( ixx,      izz, 2)

! -----
! Exclude quad corner on "Efz" averaging for y-flux; this is handled separately.

          ELSE IF ( ( ixx .NE. 0) .AND. ( ixx .NE. ntp_zones)) THEN

            Efz_amr      ( ixx, iyy, izz, ii) =
     1    ( Efz_amr      ( ixx, iyy, izz, ii)
     1    + Flux_recev_d ( ixx,      izz, 2)) * 0.5

! -----

          END IF

! ----------

        ELSE

! ----------

          Efz_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_d ( ixx,      izz, 2)

! ----------

        END IF  ! IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) ...

! ----------

      END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_d" to "Efx_amr".

      DO izz = 0 - NGhost, nr_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

        iyy = 0

! ----------

        IF ( imask_yface ( ixx, iyy) .EQ. 1) THEN  ! Exclude any copying or averaging at x- faces abutting tri-corners.

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( izz .GE. 0) .AND. ( izz .LE. nr_zones )) THEN

          Efx_amr      ( ixx, iyy, izz, ii) =
     1  ( Efx_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_d ( ixx,      izz, 3)) * 0.5

        ELSE

          Efx_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_d ( ixx,      izz, 3)

        END IF

        END IF  ! Exclude any copying or averaging at x- faces abutting tri-corners.

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! The tri-corner "Efz" are accumulated. Now correctly average the "Efz"
! by dividing them by 3.

      DO ii = 1, npatch_per_proc

! ----------------------------------------

        itri_patch_LD = itri_patch_amr ( 1, ii)
        itri_patch_RD = itri_patch_amr ( 2, ii)
        itri_patch_RU = itri_patch_amr ( 3, ii)
        itri_patch_LU = itri_patch_amr ( 4, ii)

! ----------------------------------------

        DO izz = 1 - NGhost, nr_zones  + NGhost

! --------------------

          IF ( itri_patch_LD .EQ. 1)
     1    Efz_amr ( 0, 0, izz, ii) =
     1    Efz_amr ( 0, 0, izz, ii) / 3.0

          IF ( itri_patch_LU .EQ. 1)
     1    Efz_amr ( 0, ntp_zones, izz, ii) =
     1    Efz_amr ( 0, ntp_zones, izz, ii) / 3.0

          IF ( itri_patch_RD .EQ. 1)
     1    Efz_amr ( ntp_zones, 0, izz, ii) =
     1    Efz_amr ( ntp_zones, 0, izz, ii) / 3.0

          IF ( itri_patch_RU .EQ. 1)
     1    Efz_amr ( ntp_zones, ntp_zones, izz, ii) =
     1    Efz_amr ( ntp_zones, ntp_zones, izz, ii) / 3.0

! --------------------

        END DO  ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost"

! ----------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_Y
