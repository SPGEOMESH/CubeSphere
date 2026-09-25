! contained in RIEMANN_3D_CS_AMR_SUPPORT









      SUBROUTINE LOAD_LOCAL_BUFFERS_UZ ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after each timestep iteration, so that the
! final, updated data is in "U_amr, entropygr_amr, Bxf_amr, Byf_amr,
! Bzf_amr". The date is stored in the big container arrays with this
! ordering of variables. This data is loaded into "U_snd_t_amr, U_snd_b_amr".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        izt, izb, izs

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

      U_snd_t_amr ( :, :, :, :, :) = 0.0
      U_snd_b_amr ( :, :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! Store top and bottom z-data in "U_snd_t_amr, U_snd_b_amr".

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "U_snd_t_amr, U_snd_b_amr".  First for the zone centered vairables,
! "U_amr, entropygr_amr" then for the face centered variables,
! "Bxf_amr, Byf_amr, Bzf_amr".

! ------------------------------------------------------------
! Storing zone centered variables "U_amr, entropygr_amr" data in
! "U_snd_t_amr, U_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_t_amr, U_snd_b_amr, U_amr,             &
!$OMP&    entropygr_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, 0
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izs = izz + NGhost     ! location in "U_snd_"

          izt = izz + nr_zones   ! location of top slice
          izb = izz + NGhost     ! location of bottom slice


          U_snd_t_amr ( ixx, iyy, izs, 1: n_cc_components, ii) =
     1    U_amr   ( ixx, iyy, izt, 1: n_cc_components, ii)

          U_snd_b_amr ( ixx, iyy, izs, 1: n_cc_components, ii) =
     1    U_amr   ( ixx, iyy, izb, 1: n_cc_components, ii)


          U_snd_t_amr       ( ixx, iyy, izs, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixx, iyy, izt, ii)

          U_snd_b_amr       ( ixx, iyy, izs, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones  + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, 0".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Bxf_amr" data in "U_snd_t_amr, U_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_t_amr, U_snd_b_amr, Bxf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, 0
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izs = izz + NGhost     ! location in "U_snd_"

          izt = izz + nr_zones   ! location of top slice
          izb = izz + NGhost     ! location of bottom slice


          U_snd_t_amr ( ixx, iyy, izs, n_cc_components + 2, ii) =
     1    Bxf_amr ( ixx, iyy, izt, ii)

          U_snd_b_amr ( ixx, iyy, izs, n_cc_components + 2, ii) =
     1    Bxf_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, 0".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Byf_amr" data in "U_snd_t_amr, U_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_t_amr, U_snd_b_amr, Byf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, 0
        DO iyy = 0 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izs = izz + NGhost     ! location in "U_snd_"

          izt = izz + nr_zones   ! location of top slice
          izb = izz + NGhost     ! location of bottom slice


          U_snd_t_amr ( ixx, iyy, izs, n_cc_components + 3, ii) =
     1    Byf_amr ( ixx, iyy, izt, ii)

          U_snd_b_amr ( ixx, iyy, izs, n_cc_components + 3, ii) =
     1    Byf_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, 0".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Bzf_amr" data in "U_snd_t_amr, U_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_t_amr, U_snd_b_amr, Bzf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, 0
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izs = izz + NGhost     ! location in "U_snd_"

          izt = izz + nr_zones   ! location of top slice
          izb = izz + NGhost     ! location of bottom slice


          U_snd_t_amr ( ixx, iyy, izs, n_cc_components + 4, ii) =
     1    Bzf_amr ( ixx, iyy, izt, ii)

          U_snd_b_amr ( ixx, iyy, izs, n_cc_components + 4, ii) =
     1    Bzf_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, 0".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE LOAD_LOCAL_BUFFERS_UZ

