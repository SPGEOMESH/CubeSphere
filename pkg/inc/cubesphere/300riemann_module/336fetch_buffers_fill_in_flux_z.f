! contained in RIEMANN_3D_CS_AMR_SUPPORT














      SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_Z ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S. and Dinshaw Balsara

! This subroutine is called after the buffers are loaded in the z-direction.
! Each patch should first know its own "( ix_patch, iy_patch, iz_patch )".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,
     1        ixx, iyy, izz, ixc, iyc, izc,
     1        ixn, iyn, izn

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank,
     1        ix_patch, iy_patch, iz_patch

      INTEGER num_of_reals_xy, num_of_reals_z

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
! is to be able to extract the top sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! top sibling patch.

! STEP 1.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch of
! "Flux_snd_b ( :, :, :, yr_patch_in_proc)", just as we would do
! in a Cartesian mesh code. This is saved in "Flux_recev_t" on the
! processor that owns "(my_sector, my_patch_in_proc)".

! STEP 1.c) Now that the data is in "Flux_recev_t" we can unpack it into
! "z-Fluctuations, EntFlux_z_amr, Efx_amr, Efy_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 1.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! top sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access top sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for top sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 6)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 6)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 6)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------------------------------------
! STEP 1.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch
! of "Flux_snd_b ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_t".

      IF ( yr_sector .GE. 1 ) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_b ( :, :, :, yr_patch_in_proc)" to "Flux_recev_t".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_ZB_Msg_handle, ierr)

      window_offset = num_of_reals_z * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_t, num_of_reals_z, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_z,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_ZB_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_ZB_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) Unpack data from "Flux_recev_t" to
! "z-Fluctuations, EntFlux_z_amr, Efx_amr, Efy_amr".

      DO iyy = 1 - NGhost, ntp_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        izz = nr_zones

! ----------

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding face.
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones)) THEN

          EntFlux_z_amr ( ixx, iyy, izz, ii) =
     1  ( EntFlux_z_amr ( ixx, iyy, izz, ii)
     1  + Flux_recev_t  ( ixx, iyy,      1)) * 0.5

          Fluct_z_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_z_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_t ( ixx, iyy,      ip1: ip2)) * 0.5

          Fluct_z_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_z_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_t ( ixx, iyy,      im1: im2)) * 0.5

        ELSE                              ! Ghost face.

          EntFlux_z_amr ( ixx, iyy, izz, ii) =
     1    Flux_recev_t  ( ixx, iyy,      1)

          Fluct_z_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_t ( ixx, iyy,      ip1: ip2)

          Fluct_z_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_t ( ixx, iyy,      im1: im2)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_t" to "Efx_amr".

      DO iyy = 0 - NGhost, ntp_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        izz = nr_zones

! ----------

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( iyy .GE. 0) .AND. ( iyy .LE. ntp_zones)) THEN

          Efx_amr      ( ixx, iyy, izz, ii) =
     1  ( Efx_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_t ( ixx, iyy, 2)) * 0.5

        ELSE                              ! Ghost face.

          Efx_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_t ( ixx, iyy, 2)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_t" to "Efy_amr".

      DO iyy = 1 - NGhost, ntp_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones  + NGhost

        izz = nr_zones

! ----------

        IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones)) THEN

          Efy_amr      ( ixx, iyy, izz, ii) =
     1  ( Efy_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_t ( ixx, iyy, 3)) * 0.5

        ELSE                              ! Ghost face.

          Efy_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_t ( ixx, iyy, 3)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "IF ( yr_sector .GE. 1 )"

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

! STEP 2.a) We loop over all the patches on the current processor "me".
! We extract "(my_sector, my_patch_in_proc)" for that patch. Our job
! is to be able to extract the bottom sibling of this patch in a way where
! we can do messaging.

! This means that we want "(yr_sector, yr_patch_in_proc, yr_proc)" for that
! bottom sibling patch.

! STEP 2.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch of
! "Flux_snd_t ( :, :, :, yr_patch_in_proc)", just as we would do
! in a Cartesian mesh code. This is saved in "Flux_recev_b" on the
! processor that owns "(my_sector, my_patch_in_proc)".

! STEP 2.c) Now that the data is in "Flux_recev_b" we can unpack it into
! "z-Fluctuations, EntFlux_z_amr, Efx_amr, Efy_amr".

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
! Use it to obtain "(yr_sector, yr_patch_in_proc, yr_proc)" for its
! top sibling patch.

      my_sector        = iactive_FromProc ( ii, 1)
      my_patch_in_proc = ii

! Gives us global patch number in a sequential space.

      jj = my_patch_in_proc + (me - 1) * npatch_per_proc

! "(my_sector, my_patch_in_sector" will allow us to access top sibling
! through "isibling_sector_t, isibling_number_t, isibling_proc_t".

      my_patch_in_sector = jj - (my_sector - 1)
     1                   * Max_base_patches_per_sector

! "(yr_sector, yr_patch_in_sector, yr_proc)" are easily obtained for top sibling.

      yr_sector = isibling_sector_t ( my_sector, my_patch_in_sector, 7)
      yr_patch_in_sector =
     1            isibling_number_t ( my_sector, my_patch_in_sector, 7)
      yr_proc =   isibling_proc_t   ( my_sector, my_patch_in_sector, 7)

! From "yr_patch_in_sector" we want "yr_patch_in_proc" in order to access
! sending array later on. Basically reverse the above logic.

! Obtain global patch number for "yr_patch_in_sector".
      jj = yr_patch_in_sector + (yr_sector - 1)
     1                   * Max_base_patches_per_sector

      yr_patch_in_proc = jj - (yr_proc - 1) * npatch_per_proc

! ------------------------------------------------------------
! STEP 1.b) "IF ( yr_sector .GE. 1 )" then we do a one-sided fetch
! of "Flux_snd_t ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_b".

      IF ( yr_sector .GE. 1 ) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_t ( :, :, :, yr_patch_in_proc)" to "Flux_recev_b".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_ZT_Msg_handle, ierr)

      window_offset = num_of_reals_z * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_b, num_of_reals_z, MPI_REAL8,    ! local data on my.
     1               yr_rank, window_offset, num_of_reals_z,     ! remote data on yr.
     1               MPI_REAL8, Data_Flux_ZT_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_ZT_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 1.c) Unpack data from "Flux_recev_b" to
! "z-Fluctuations, EntFlux_z_amr, Efx_amr, Efy_amr".

      DO iyy = 1 - NGhost, ntp_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        izz = 0

! ----------

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding face.
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones)) THEN

          EntFlux_z_amr ( ixx, iyy, izz, ii) =
     1  ( EntFlux_z_amr ( ixx, iyy, izz, ii)
     1  + Flux_recev_b  ( ixx, iyy,      1)) * 0.5

          Fluct_z_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_z_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_b ( ixx, iyy,      ip1: ip2)) * 0.5

          Fluct_z_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  ( Fluct_z_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1    + Flux_recev_b ( ixx, iyy,      im1: im2)) * 0.5

        ELSE                              ! Ghost face.

          EntFlux_z_amr ( ixx, iyy, izz, ii) =
     1    Flux_recev_t  ( ixx, iyy,      1)

          Fluct_z_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_b ( ixx, iyy,      ip1: ip2)

          Fluct_z_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1      Flux_recev_b ( ixx, iyy,      im1: im2)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_b" to "Efx_amr".

      DO iyy = 0 - NGhost, ntp_zones  + NGhost
      DO ixx = 1 - NGhost, ntp_zones  + NGhost

        izz = 0

! ----------

        IF ( ( ixx .GE. 1) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( iyy .GE. 0) .AND. ( iyy .LE. ntp_zones)) THEN

          Efx_amr      ( ixx, iyy, izz, ii) =
     1  ( Efx_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_b ( ixx, iyy, 2)) * 0.5

        ELSE                              ! Ghost face.

          Efx_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_b ( ixx, iyy, 2)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones  + NGhost".

! --------------------
! Unpack data from "Flux_recev_b" to "Efy_amr".

      DO iyy = 1 - NGhost, ntp_zones  + NGhost
      DO ixx = 0 - NGhost, ntp_zones  + NGhost

        izz = 0

! ----------

        IF ( ( ixx .GE. 0) .AND. ( ixx .LE. ntp_zones) .AND.  ! Active Coinciding edge.
     1       ( iyy .GE. 1) .AND. ( iyy .LE. ntp_zones)) THEN

          Efy_amr      ( ixx, iyy, izz, ii) =
     1  ( Efy_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_b ( ixx, iyy, 3)) * 0.5

        ELSE                              ! Ghost face.

          Efy_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_b ( ixx, iyy, 3)

        END IF

! ----------

      END DO ! End of "DO ixx = 1 - NGhost, ntp_zones  + NGhost".
      END DO ! End of "DO iyy = 0 - NGhost, ntp_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "IF ( yr_sector .GE. 1 )"

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

      END SUBROUTINE FETCH_BUFFERS_FILL_IN_Flux_Z
