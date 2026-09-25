! contained in RIEMANN_3D_CS_AMR_SUPPORT













      SUBROUTINE LOAD_LOCAL_BUFFERS_Flux_Z ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after the build of fluxes and electric fields,
! so that the fluxes and electric fields are synchronized at the coinciding
! faces.  The date is stored in the big container arrays, "z-Fluctuations,
! EntFlux_z_amr, Efx_amr, Efy_amr". This data is loaded into
! "Flux_snd_t_amr, Flux_snd_b_amr".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        izt, izb, izs

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

      Flux_snd_t_amr ( :, :, :, :) = 0.0
      Flux_snd_b_amr ( :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! Loop over all the patches in this processor "me".
! Store top and bottom z--data in "Flux_snd_t_amr, Flux_snd_b_amr".


      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! Loading "Flux_snd_t_amr, Flux_snd_b_amr". First for the face centered fluxes,
! "z-Fluctuations, EntFlux_z_amr" then for the edge centered variables,
! "Efx_amr, Efy_amr".

! ------------------------------------------------------------
! Storing face centered fluxes "z-Fluctuations, EntFlux_z_amr" data in
! "Flux_snd_t_amr, Flux_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_amr, Flux_snd_b_amr,              &
!$OMP&    EntFlux_z_amr, Fluct_z_Pl_amr, Fluct_z_Mi_amr)                &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

! ----------
! ip1 == 3 + 1
! ip2 == 3 + n_cc_components

          Flux_snd_t_amr ( ixx, iyy,      ip1: ip2,           ii) =
     1    Fluct_z_Pl_amr ( ixx, iyy, izt, 1: n_cc_components, ii)

          Flux_snd_b_amr ( ixx, iyy,      ip1: ip2,           ii) =
     1    Fluct_z_Pl_amr ( ixx, iyy, izb, 1: n_cc_components, ii)

! ----------
! im1 == 3 + n_cc_components + 1
! im2 == 3 + n_cc_components + n_cc_components

          Flux_snd_t_amr ( ixx, iyy,      im1: im2,           ii) =
     1    Fluct_z_Mi_amr ( ixx, iyy, izt, 1: n_cc_components, ii)

          Flux_snd_b_amr ( ixx, iyy,      im1: im2,           ii) =
     1    Fluct_z_Mi_amr ( ixx, iyy, izb, 1: n_cc_components, ii)

! ----------

          Flux_snd_t_amr ( ixx, iyy,      1, ii) =
     1     EntFlux_z_amr ( ixx, iyy, izt, ii)

          Flux_snd_b_amr ( ixx, iyy,      1, ii) =
     1     EntFlux_z_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efx_amr" data in "Flux_snd_t_amr, Flux_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_amr, Flux_snd_b_amr, Efx_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO iyy = 0 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

          Flux_snd_t_amr ( ixx, iyy,      2, ii) =
     1           Efx_amr ( ixx, iyy, izt, ii)

          Flux_snd_b_amr ( ixx, iyy,      2, ii) =
     1           Efx_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".

!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efy_amr" data in "Flux_snd_t_amr, Flux_snd_b_amr".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_amr, Flux_snd_b_amr, Efy_amr )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

          Flux_snd_t_amr ( ixx, iyy,      3, ii) =
     1           Efy_amr ( ixx, iyy, izt, ii)

          Flux_snd_b_amr ( ixx, iyy,      3, ii) =
     1           Efy_amr ( ixx, iyy, izb, ii)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

! ----------------------------------------------------------------------

      END SUBROUTINE LOAD_LOCAL_BUFFERS_Flux_Z


