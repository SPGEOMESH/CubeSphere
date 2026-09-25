! contained in RIEMANN_3D_CS_AMR_SUPPORT














      SUBROUTINE SYNCHRONIZE_EFZ_QUAD_CORNER ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after the synchronization of fluxes and electric
! fields at the coinciding faces. The "Efz" at the quad corner needs a
! further synchrnization so as to correctly average the "Efz" edge at the
! quad-corner.
!
! Consider the following quad-corner. The edge "Efz" is perpendicular to the screen.
!
!                            3  |  4
!                               |
!                        -------|-------
!                               |
!                            1  |  2
!
!
! The previous "FETCH..Flux_X ( )" subroutine synchronize at the x-face.
! And it assigned the following new electric fields at the quad-corner:
!
!            Efz_new_1 = ( Efz_1 + Efz_2 ) / 2
!            Efz_new_2 = ( Efz_1 + Efz_2 ) / 2
!
!            Efz_new_3 = ( Efz_3 + Efz_4 ) / 2
!            Efz_new_4 = ( Efz_3 + Efz_4 ) / 2
!
! Now the quad-corner electric fields are,
!
!        ( Efz_3 + Efz _4 ) / 2 |  ( Efz_3 + Efz _4 ) / 2
!                               |
!                        -------|-------
!                               |
!        ( Efz_1 + Efz _2 ) / 2 |  ( Efz_1 + Efz _2 ) / 2
!
!
! Note that the quad-corner Efz is not changed in the "FETCH..Flux_Y ( )" subroutine.
! This subroutine is called next, to handle the same.
!
! In this subroutine, the quad-corner Efz is synchronized in two distinct steps:
!     i) LOAD  - Efz_new values.
!    ii) FETCH - synchronize along the y-face.
!
! After the synchronization, the quad-corner Efz in all four shared patches,
! (1), (2), (3), (4) will have the value of:
!
!     "( Efz_1 + Efz _2 + Efz_3 + Efz_4 ) / 4"
!
! Thus completing the synchronization at the quad-corner Efz.
!
! --------------------
! The above mentioned technique works well in all the circumstances,
! except the for the "y-top" face when it interfaces with the next sector.
! This because of how the 'xy-axes' are defined in each sector.
!
! Hence, the Efz-quad corners which lie in the interfacing "y-top" face,
! are assigned separately; in step 5.
!
! ----------------------------------------------------------------------
!
! STEP 1 & 2: sync. "y-top"    quad corner Efz (y-top touching same sector patch)
! STEP 3    : sync. "y-bottom" quad corner Efz
! STEP 4 & 5: sync. "y-top"    quad corner Efz (y-top touching difft. sector patch)
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        ixs, iys, izs

      INTEGER ixc, iyc, ixn, iyn, izn, ixxp1, iyyp1, ixx0

      INTEGER my_sector, my_patch_in_proc, my_patch_in_sector,
     1        yr_sector, yr_patch_in_sector, yr_patch_in_proc, yr_proc,
     1        ii, jj, yr_rank

      INTEGER num_of_reals_xy, num_of_reals_z

      INTEGER itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

      INTEGER imask_xface ( 0 - NGhost: ntp_zones + NGhost,
     1                      0 - NGhost: ntp_zones + NGhost)

! ----------------------------------------------------------------------
! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

      Flux_snd_u_amr ( :, :, :, :) = 0.0
      Flux_snd_d_amr ( :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! STEP 1) Loop over all the patches in this processor "me".
! Store up and down y-data in "Flux_snd_u_amr, Flux_snd_d_amr", only the "Efz"
! at the quad-corner.

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Storing edge centered variable "Efz_amr" data in "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr, Efz_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

        DO izz = 1 - NGhost, nr_zones  + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.

! ----------

          ixx = 0

          Flux_snd_u_amr ( ixx,      izz, 2, ii) =
     1           Efz_amr ( ixx, iyu, izz, ii)

          Flux_snd_d_amr ( ixx,      izz, 2, ii) =
     1           Efz_amr ( ixx, iyd, izz, ii)

! ----------

          ixx = ntp_zones

          Flux_snd_u_amr ( ixx,      izz, 2, ii) =
     1           Efz_amr ( ixx, iyu, izz, ii)

          Flux_snd_d_amr ( ixx,      izz, 2, ii) =
     1           Efz_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO
! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Send complete; now synchronize and the recieve the Efz.

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
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
! STEP 2. Receive 'y-top' - "Flux_recev_u" if it is in the same sector.
! when "Flux_recev_u" lies in different sector, it has to be done
! separately at last; and it is done in step 5.

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 2.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
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
! STEP 2.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d_amr ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_u".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_d_amr ( :, :, :, yr_patch_in_proc)" to "Flux_recev_u".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YD_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_u, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YD_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YD_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 2.c) Unpack data from "Flux_recev_u" to "Efz_amr" at quad corner.

      DO izz = 1 - NGhost, nr_zones  + NGhost

! ------------------------------

        iyy = ntp_zones

! ----------
! Average "Efz" at the Left-Upper Quad corner.

        IF ( itri_patch_LU .NE. 1) THEN

          ixx = 0

          Efz_amr      ( ixx, iyy, izz, ii) =
     1  ( Efz_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, 2)) * 0.5

        ENDIF

! ----------
! Average "Efz" at the Right-Upper Quad corner.

        IF ( itri_patch_RU .NE. 1) THEN

          ixx = ntp_zones

          Efz_amr      ( ixx, iyy, izz, ii) =
     1  ( Efz_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_u ( ixx,      izz, 2)) * 0.5

        ENDIF

! ------------------------------

      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 3) Receive 'y-bottom' from "Flux_recev_d". Unlike the 'y-top',
! here the values are received from the neighbouring sectors as well.
! Averaging of Efz from the neighbouring sector works for the 'y-bottom'
! Efz; due to how the xy-axes of neighbouring sectors defined.

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 3.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
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

! ------------------------------------------------------------
! STEP 3.b) "If ( my_sector .EQ. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_u_amr ( :, :, :, yr_patch_in_proc)", and save in "Flux_recev_d".

      IF ( my_sector .EQ. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_u_amr ( :, :, :, yr_patch_in_proc)" to "Flux_recev_d".

      yr_rank = yr_proc - 1     ! To get to C indexing in MPI.

      CALL MPI_WIN_LOCK ( MPI_LOCK_SHARED, yr_rank, 0,
     1                    Data_Flux_YU_Msg_handle, ierr)

      window_offset = num_of_reals_xy * ( yr_patch_in_proc - 1)

      CALL MPI_GET ( Flux_recev_d, num_of_reals_xy, MPI_REAL8,   ! local data on my.
     1               yr_rank, window_offset, num_of_reals_xy,    ! remote data on yr.
     1               MPI_REAL8, Data_Flux_YU_Msg_handle, ierr)

      CALL MPI_WIN_UNLOCK ( yr_rank, Data_Flux_YU_Msg_handle, ierr)

! ------------------------------------------------------------
! STEP 3.c) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_d_amr ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_u".
! Index flipping is done so that "Flux_recev_u" is unpacked into "Flux_recev_d"
! with the flipping.

      ELSE

! ------------------------------------------------------------
! Fetch data from "Flux_snd_d_amr ( :, :, :, yr_patch_in_proc)" to "Flux_recev_u".

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

! Only the "Efz" is unpacked here.

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

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .EQ. yr_sector)".

! ------------------------------------------------------------
! STEP 3.d) Unpack data from "Flux_recev_d" to "Efz_amr".

      DO izz = 1 - NGhost, nr_zones  + NGhost

! ------------------------------

      iyy = 0

! ----------
! Average "Efz" at the Left-Down Quad corner.

        IF ( itri_patch_LD .NE. 1) THEN

          ixx = 0

          Efz_amr      ( ixx, iyy, izz, ii) =
     1  ( Efz_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_d ( ixx,      izz, 2)) * 0.5

        END IF

! ----------
! Average "Efz" at the Right-Down Quad corner.

        IF ( itri_patch_RD .NE. 1) THEN

          ixx = ntp_zones

          Efz_amr      ( ixx, iyy, izz, ii) =
     1  ( Efz_amr      ( ixx, iyy, izz, ii)
     1  + Flux_recev_d ( ixx,      izz, 2)) * 0.5

        END IF

! ------------------------------

      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 4) Till here, the quad-corner "Efz" at the bottom-face is complete;
! and the top-face is also completed for same-sector faces.
! Now synchronize the top-face Efz; when it is shared with a different sector.

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

      Flux_snd_r_amr ( :, :, :, :) = 0.0
      Flux_snd_l_amr ( :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! First store right and left x-data in "Flux_snd_r_amr, Flux_snd_l_amr".

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "Flux_snd_r_amr, Flux_snd_l_amr".
! Storing edge centered variable "Efz_amr" data in "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr, Efz_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

        DO izz = 1 - NGhost, nr_zones  + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.

! ----------

          iyy = 0

          Flux_snd_r_amr (      iyy, izz, 3, ii) =
     1           Efz_amr ( ixr, iyy, izz, ii)

          Flux_snd_l_amr (      iyy, izz, 3, ii) =
     1           Efz_amr ( ixl, iyy, izz, ii)

! ----------

          iyy = ntp_zones

          Flux_snd_r_amr (      iyy, izz, 3, ii) =
     1           Efz_amr ( ixr, iyy, izz, ii)

          Flux_snd_l_amr (      iyy, izz, 3, ii) =
     1           Efz_amr ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Send complete; now synchronize and the recieve the Efz.

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! STEP 5. Receive 'y-top' - "Flux_recev_u" if it is coming from different
! sector. Note that the same sector synchronization for y-top is completed
! in step 2.

! ----------------------------------------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------------------------------------
! STEP 5.a) Extract "(my_sector, my_patch_in_proc)" for that patch.
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

! ------------------------------------------------------------
! STEP 5.b) "If ( my_sector .NE. yr_sector)" then we do a one-sided fetch
! of "Flux_snd_r_amr ( :, :, :, :, yr_patch_in_proc)". This is saved in "Flux_recev_l".
! Index flipping is done so that "Flux_recev_l" is unpacked into "Flux_recev_u"
! with the flipping.

      IF ( my_sector .NE. yr_sector) THEN

! ------------------------------------------------------------
! Fetch data from "Flux_snd_r_amr ( :, :, :, yr_patch_in_proc)" to "Flux_recev_l".

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
! Only the "Efz" data.
! Watch out, because "Flux_recev_l" does not have data in normal ordering.

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

! ------------------------------------------------------------
! STEP 5.c) Unpack data from "Flux_recev_u" to "Efz_amr" at quad corner.

! NOTE: It is important to note here that the "Flux_recev_u" has
! correctly averged "Efz" from all 4-contributions at the quadcorner.
! Hence, the values are --assigned-- here. Not averaged.

! ------------------------------------------------------------

      DO izz = 1 - NGhost, nr_zones  + NGhost

! ------------------------------

        iyy = ntp_zones

! ----------
! Assign the average "Efz" at the Left-Upper Quad corner.

        IF ( itri_patch_LU .NE. 1) THEN

          ixx = 0

          Efz_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, 2)

        ENDIF

! ----------
! Assign the average "Efz" at the Right-Upper Quad corner.

        IF ( itri_patch_RU .NE. 1) THEN

          ixx = ntp_zones

          Efz_amr      ( ixx, iyy, izz, ii) =
     1    Flux_recev_u ( ixx,      izz, 2)

        ENDIF

! ------------------------------

      END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

! ------------------------------------------------------------

      END IF  ! End of "If ( my_sector .NE. yr_sector)".

! ------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE SYNCHRONIZE_EFZ_QUAD_CORNER
