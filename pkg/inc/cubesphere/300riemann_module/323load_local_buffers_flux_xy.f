! contained in RIEMANN_3D_CS_AMR_SUPPORT















      SUBROUTINE LOAD_LOCAL_BUFFERS_Flux_XY ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after the build of fluxes and electric fields,
! so that the fluxes and electric fields are synchronized at the coinciding
! faces.  The date is stored in the big container arrays,
! "EntFlux_x_amr, EntFlux_y_amr, Efx_amr, Efy_amr, Efz_amr".
! This data is loaded into
! "Flux_snd_r_amr, Flux_snd_l_amr, Flux_snd_u_amr, Flux_snd_d_amr".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        ixs, iys, izs

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
! Initialize the send arrays with zero.

      Flux_snd_r_amr ( :, :, :, :) = 0.0
      Flux_snd_l_amr ( :, :, :, :) = 0.0
      Flux_snd_u_amr ( :, :, :, :) = 0.0
      Flux_snd_d_amr ( :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! First store right and left x-data in "Flux_snd_r_amr, Flux_snd_l_amr".
! Then store up and down y-data in "Flux_snd_u_amr, Flux_snd_d_amr".

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "Flux_snd_r_amr, Flux_snd_l_amr". First for the face centered fluxes,
! "x-Fluctuations, EntFlux_x_amr" then for the edge centered variables,
! "Efy_amr, Efz_amr".

! ------------------------------------------------------------
! Storing face centered fluxes "EntFlux_x_amr" data in
! "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr,              &
!$OMP&    EntFlux_x_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice
          ixr = ntp_zones      ! location of right slice

! ----------
! ip1 == 3 + 1
! ip2 == 3 + n_cc_components

          Flux_snd_r_amr (      iyy, izz, ip1: ip2,           ii) =
     1    Fluct_x_Pl_amr ( ixr, iyy, izz, 1: n_cc_components, ii)

          Flux_snd_l_amr (      iyy, izz, ip1: ip2,           ii) =
     1    Fluct_x_Pl_amr ( ixl, iyy, izz, 1: n_cc_components, ii)

! ----------
! im1 == 3 + n_cc_components + 1
! im2 == 3 + n_cc_components + n_cc_components

          Flux_snd_r_amr (      iyy, izz, im1: im2,           ii) =
     1    Fluct_x_Mi_amr ( ixr, iyy, izz, 1: n_cc_components, ii)

          Flux_snd_l_amr (      iyy, izz, im1: im2,           ii) =
     1    Fluct_x_Mi_amr ( ixl, iyy, izz, 1: n_cc_components, ii)

! ----------

          Flux_snd_r_amr (      iyy, izz, 1, ii) =
     1     EntFlux_x_amr ( ixr, iyy, izz, ii)

          Flux_snd_l_amr (      iyy, izz, 1, ii) =
     1     EntFlux_x_amr ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efy_amr" data in "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr, Efy_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.

          Flux_snd_r_amr (      iyy, izz, 2, ii) =
     1           Efy_amr ( ixr, iyy, izz, ii)

          Flux_snd_l_amr (      iyy, izz, 2, ii) =
     1           Efy_amr ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efz_amr" data in "Flux_snd_r_amr, Flux_snd_l_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_amr, Flux_snd_l_amr, Efz_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.

          Flux_snd_r_amr (      iyy, izz, 3, ii) =
     1           Efz_amr ( ixr, iyy, izz, ii)

          Flux_snd_l_amr (      iyy, izz, 3, ii) =
     1           Efz_amr ( ixl, iyy, izz, ii)

! ----------------------------------------

        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Loading "Flux_snd_u_amr, Flux_snd_d_amr". First for the face centered fluxes,
! "EntFlux_y_amr" then for the edge centered variables,
! "Efz_amr, Efx_amr".

! ------------------------------------------------------------
! Storing face centered fluxes "EntFlux_y_amr" data in
! "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr,              &
!$OMP&    EntFlux_y_amr )                                               &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice
          iyu = ntp_zones      ! location of up   slice

! ----------
! ip1 == 3 + 1
! ip2 == 3 + n_cc_components

          Flux_snd_u_amr ( ixx,      izz, ip1: ip2,           ii) =
     1    Fluct_y_Pl_amr ( ixx, iyu, izz, 1: n_cc_components, ii)

          Flux_snd_d_amr ( ixx,      izz, ip1: ip2,           ii) =
     1    Fluct_y_Pl_amr ( ixx, iyd, izz, 1: n_cc_components, ii)

! ----------
! im1 == 3 + n_cc_components + 1
! im2 == 3 + n_cc_components + n_cc_components

          Flux_snd_u_amr ( ixx,      izz, im1: im2,           ii) =
     1    Fluct_y_Mi_amr ( ixx, iyu, izz, 1: n_cc_components, ii)

          Flux_snd_d_amr ( ixx,      izz, im1: im2,           ii) =
     1    Fluct_y_Mi_amr ( ixx, iyd, izz, 1: n_cc_components, ii)

! ----------

          Flux_snd_u_amr ( ixx,      izz, 1, ii) =
     1     EntFlux_y_amr ( ixx, iyu, izz, ii)

          Flux_snd_d_amr ( ixx,      izz, 1, ii) =
     1     EntFlux_y_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efz_amr" data in "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr, Efz_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.

          Flux_snd_u_amr ( ixx,      izz, 2, ii) =
     1           Efz_amr ( ixx, iyu, izz, ii)

          Flux_snd_d_amr ( ixx,      izz, 2, ii) =
     1           Efz_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efx_amr" data in "Flux_snd_u_amr, Flux_snd_d_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_amr, Flux_snd_d_amr, Efx_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.

          Flux_snd_u_amr ( ixx,      izz, 3, ii) =
     1           Efx_amr ( ixx, iyu, izz, ii)

          Flux_snd_d_amr ( ixx,      izz, 3, ii) =
     1           Efx_amr ( ixx, iyd, izz, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE LOAD_LOCAL_BUFFERS_Flux_XY



