! contained in RIEMANN_3D_CS_AMR_SUPPORT







      SUBROUTINE LOAD_LOCAL_BUFFERS_UXY ()

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after each timestep iteration, so that the
! final, updated data is in "U_amr, entropygr_amr, Bxf_amr, Byf_amr,
! Bzf_amr". The date is stored in the big container arrays with this
! ordering of variables. This data is loaded into
! "U_snd_r_amr, U_snd_l_amr, U_snd_u_amr, U_snd_d_amr".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        ixs, iys, izs

! ----------------------------------------------------------------------

! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

      U_snd_r_amr ( :, :, :, :, :) = 0.0
      U_snd_l_amr ( :, :, :, :, :) = 0.0
      U_snd_u_amr ( :, :, :, :, :) = 0.0
      U_snd_d_amr ( :, :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! First store right and left x-data in "U_snd_r_amr, U_snd_l_amr".
! Then store up and down y-data in "U_snd_u_amr, U_snd_d_amr".

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "U_snd_r_amr, U_snd_l_amr". First for the zone centered vairables,
! "U_amr, entropygr_amr" then for the face centered variables,
! "Bxf_amr, Byf_amr, Bzf_amr".

! ------------------------------------------------------------
! Storing zone centered variables "U_amr, entropygr_amr" data in
! "U_snd_r_amr, U_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_r_amr, U_snd_l_amr, U_amr,             &
!$OMP&    entropygr_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, 0

! ----------------------------------------

          ixs = ixx + NGhost     ! location in "U_snd_"

          ixr = ixx + ntp_zones  ! location of right slice
          ixl = ixx + NGhost     ! location of left  slice


          U_snd_r_amr ( ixs, iyy, izz, 1: n_cc_components, ii) =
     1    U_amr       ( ixr, iyy, izz, 1: n_cc_components, ii)

          U_snd_l_amr ( ixs, iyy, izz, 1: n_cc_components, ii) =
     1    U_amr       ( ixl, iyy, izz, 1: n_cc_components, ii)


          U_snd_r_amr   ( ixs, iyy, izz, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixr, iyy, izz, ii)

          U_snd_l_amr   ( ixs, iyy, izz, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, 0".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variable "Bxf_amr" data in "U_snd_r_amr, U_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_r_amr, U_snd_l_amr, Bxf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 0 - NGhost, 0

! ----------------------------------------

          ixs = ixx + NGhost     ! location in "U_snd_"

          ixr = ixx + ntp_zones  ! location of right slice
          ixl = ixx + NGhost     ! location of left  slice


          U_snd_r_amr ( ixs, iyy, izz, n_cc_components + 2, ii) =
     1    Bxf_amr     ( ixr, iyy, izz, ii)

          U_snd_l_amr ( ixs, iyy, izz, n_cc_components + 2, ii) =
     1    Bxf_amr     ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, 0".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Byf_amr" data in "U_snd_r_amr, U_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_r_amr, U_snd_l_amr, Byf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 0 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, 0

! ----------------------------------------

          ixs = ixx + NGhost     ! location in "U_snd_"

          ixr = ixx + ntp_zones  ! location of right slice
          ixl = ixx + NGhost     ! location of left  slice


          U_snd_r_amr ( ixs, iyy, izz, n_cc_components + 3, ii) =
     1    Byf_amr     ( ixr, iyy, izz, ii)

          U_snd_l_amr ( ixs, iyy, izz, n_cc_components + 3, ii) =
     1    Byf_amr     ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, 0".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Bzf_amr" data in "U_snd_r_amr, U_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_r_amr, U_snd_l_amr, Bzf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, 0

! ----------------------------------------

          ixs = ixx + NGhost     ! location in "U_snd_"

          ixr = ixx + ntp_zones  ! location of right slice
          ixl = ixx + NGhost     ! location of left  slice


          U_snd_r_amr ( ixs, iyy, izz, n_cc_components + 4, ii) =
     1    Bzf_amr     ( ixr, iyy, izz, ii)

          U_snd_l_amr ( ixs, iyy, izz, n_cc_components + 4, ii) =
     1    Bzf_amr     ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, 0".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Loading "U_snd_u_amr, U_snd_d_amr". First for the zone centered vairables,
! "U_amr, entropygr_amr" then for the face centered variables,
! "Bxf_amr, Byf_amr, Bzf_amr".

! ------------------------------------------------------------
! Storing zone centered variables "U_amr, entropygr_amr" data in
! "U_snd_u_amr, U_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iys, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_u_amr, U_snd_d_amr, U_amr,             &
!$OMP&    entropygr_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, 0
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iys = iyy + NGhost     ! location in "U_snd_"

          iyu = iyy + ntp_zones  ! location of up slice
          iyd = iyy + NGhost     ! location of down slice


          U_snd_u_amr ( ixx, iys, izz, 1: n_cc_components, ii) =
     1    U_amr   ( ixx, iyu, izz, 1: n_cc_components, ii)

          U_snd_d_amr ( ixx, iys, izz, 1: n_cc_components, ii) =
     1    U_amr   ( ixx, iyd, izz, 1: n_cc_components, ii)


          U_snd_u_amr       ( ixx, iys, izz, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixx, iyu, izz, ii)

          U_snd_d_amr       ( ixx, iys, izz, n_cc_components + 1, ii) =
     1    entropygr_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, 0".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Bxf_amr" data in "U_snd_u_amr, U_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iys, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_u_amr, U_snd_d_amr, Bxf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, 0
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iys = iyy + NGhost     ! location in "U_snd_"

          iyu = iyy + ntp_zones  ! location of up slice
          iyd = iyy + NGhost     ! location of down slice


          U_snd_u_amr ( ixx, iys, izz, n_cc_components + 2, ii) =
     1    Bxf_amr ( ixx, iyu, izz, ii)

          U_snd_d_amr ( ixx, iys, izz, n_cc_components + 2, ii) =
     1    Bxf_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, 0".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Byf_amr" data in "U_snd_u_amr, U_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iys, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_u_amr, U_snd_d_amr, Byf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 0 - NGhost, 0
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iys = iyy + NGhost     ! location in "U_snd_"

          iyu = iyy + ntp_zones  ! location of up slice
          iyd = iyy + NGhost     ! location of down slice

! --------------------

          U_snd_u_amr ( ixx, iys, izz, n_cc_components + 3, ii) =
     1    Byf_amr ( ixx, iyu, izz, ii)

          U_snd_d_amr ( ixx, iys, izz, n_cc_components + 3, ii) =
     1    Byf_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, 0".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing face centered variables "Bzf_amr" data in "U_snd_u_amr, U_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iys, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, U_snd_u_amr, U_snd_d_amr, Bzf_amr )          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, 0
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iys = iyy + NGhost     ! location in "U_snd_"

          iyu = iyy + ntp_zones  ! location of up slice
          iyd = iyy + NGhost     ! location of down slice


          U_snd_u_amr ( ixx, iys, izz, n_cc_components + 4, ii) =
     1    Bzf_amr ( ixx, iyu, izz, ii)

          U_snd_d_amr ( ixx, iys, izz, n_cc_components + 4, ii) =
     1    Bzf_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, 0".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc"

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE LOAD_LOCAL_BUFFERS_UXY
