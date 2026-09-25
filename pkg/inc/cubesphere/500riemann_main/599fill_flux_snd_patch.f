









      SUBROUTINE FILL_FLUX_SND_PATCH

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           EntFlux_x, EntFlux_y, EntFlux_z,
     1           Efx, Efy, Efz,

     1           Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1           Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi,

     1           Flux_snd_r_patch, Flux_snd_l_patch, Flux_snd_u_patch,
     1           Flux_snd_d_patch, Flux_snd_t_patch, Flux_snd_b_patch)

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called after building the fluxes and electric fields,

! In this subroutine, only the boundary data that is stored in the
! "EntFlux_xyz, Efxyz, Fluct_xyz_Pl, Fluct_xyz_Mi" arrays are loaded into the
! "Flux_snd_r_patch", "Flux_snd_l_patch", "Flux_snd_u_patch",
! "Flux_snd_d_patch", "Flux_snd_t_patch", "Flux_snd_b_patch" arrays.

! After that the "Flux_snd_X_patch" arrays will be moved to the CPU and
! then they will fill their respective "_amr" arrays on the CPU.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes

! ----------

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost)::
     1      EntFlux_x, EntFlux_y, EntFlux_z

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost):: Efx, Efy, Efz

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                                          n_cc_components) ::
     1      Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1      Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi

! ----------

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  3 + 2*n_cc_components) ::
     1   Flux_snd_r_patch, Flux_snd_l_patch

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                  3 + 2*n_cc_components) ::
     1   Flux_snd_u_patch, Flux_snd_d_patch

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  3 + 2*n_cc_components) ::
     1   Flux_snd_t_patch, Flux_snd_b_patch

! ----------------------------------------------------------------------
! Variables that are declared below are internal to this subroutine.

      INTEGER ixx, iyy, izz, icomp,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        izt, izb,

     1        ixs, iys, izs


! Integer constants, that specify the beginning and end locations of
! fluctuation terms in the Flux-communication arrays.

      INTEGER, PARAMETER :: n_cc_comp = NVAR + NFLUID,

     1                      ip1 = 3 + 1,
     1                      ip2 = 3 + n_cc_comp,

     1                      im1 = 3 + n_cc_comp + 1,
     1                      im2 = 3 + n_cc_comp + n_cc_comp

! ----------------------------------------------------------------------
! Initialize the send arrays with zero.

!$acc kernels
      Flux_snd_r_patch ( :, :, :) = 0.0
      Flux_snd_l_patch ( :, :, :) = 0.0
      Flux_snd_u_patch ( :, :, :) = 0.0
      Flux_snd_d_patch ( :, :, :) = 0.0
      Flux_snd_t_patch ( :, :, :) = 0.0
      Flux_snd_b_patch ( :, :, :) = 0.0
!$acc end kernels

! ----------------------------------------------------------------------
! First store right and left x-data in "Flux_snd_r_patch, Flux_snd_l_patch".
! Then store up and down y-data in "Flux_snd_u_patch, Flux_snd_d_patch".

! ----------------------------------------------------------------------
! Loading "Flux_snd_r_patch, Flux_snd_l_patch". First for the face centered fluxes,
! "x-Fluctuations, EntFlux_x" then for the edge centered variables,
! "Efy, Efz".

! ------------------------------------------------------------
! Storing face centered fluxes "x-Fluctuations, EntFlux_x" data in
! "Flux_snd_r_patch, Flux_snd_l_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl, i1, i2 )      &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_patch, Flux_snd_l_patch,          &
!$OMP&    EntFlux_x, Fluct_x_Pl, Fluct_x_Mi)                            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_r_patch, Flux_snd_l_patch, EntFlux_x,                &
!$acc&    Fluct_x_Pl, Fluct_x_Mi)

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, ixr, ixl, i1, i2 )

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice
          ixr = ntp_zones      ! location of right slice

! ----------
! ip1 == 3 + 1
! ip2 == 3 + n_cc_components

          Flux_snd_r_patch (      iyy, izz, ip1: ip2) =
     1          Fluct_x_Pl ( ixr, iyy, izz,   1: n_cc_components)

          Flux_snd_l_patch (      iyy, izz, ip1: ip2) =
     1          Fluct_x_Pl ( ixl, iyy, izz,   1: n_cc_components)

! ----------
! im1 == 3 + n_cc_components + 1
! im2 == 3 + n_cc_components + n_cc_components

          Flux_snd_r_patch (      iyy, izz, im1: im2) =
     1          Fluct_x_Mi ( ixr, iyy, izz,   1: n_cc_components)

          Flux_snd_l_patch (      iyy, izz, im1: im2) =
     1          Fluct_x_Mi ( ixl, iyy, izz,   1: n_cc_components)

! ----------

          Flux_snd_r_patch (      iyy, izz, 1) =
     1           EntFlux_x ( ixr, iyy, izz)

          Flux_snd_l_patch (      iyy, izz, 1) =
     1           EntFlux_x ( ixl, iyy, izz)

! ----------------------------------------

        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efy" data in "Flux_snd_r_patch, Flux_snd_l_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_patch, Flux_snd_l_patch, Efy )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_r_patch, Flux_snd_l_patch, Efy )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, ixr, ixl )

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO iyy = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.

          Flux_snd_r_patch (      iyy, izz, 2) =
     1                 Efy ( ixr, iyy, izz)

          Flux_snd_l_patch (      iyy, izz, 2) =
     1                 Efy ( ixl, iyy, izz)

! ----------------------------------------

        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efz" data in "Flux_snd_r_patch, Flux_snd_l_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, ixr, ixl )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_r_patch, Flux_snd_l_patch, Efz )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_r_patch, Flux_snd_l_patch, Efz )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, ixr, ixl )

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO iyy = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          ixl = 0              ! location of left  slice.
          ixr = ntp_zones      ! location of right slice.

          Flux_snd_r_patch (      iyy, izz, 3) =
     1                 Efz ( ixr, iyy, izz)

          Flux_snd_l_patch (      iyy, izz, 3) =
     1                 Efz ( ixl, iyy, izz)

! ----------------------------------------

        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Loading "Flux_snd_u_patch, Flux_snd_d_patch". First for the face centered fluxes,
! "y-Fluctuations, EntFlux_y" then for the edge centered variables,
! "Efz, Efx".

! ------------------------------------------------------------
! Storing face centered fluxes "y-Fluctuations, EntFlux_y" data in
! "Flux_snd_u_patch, Flux_snd_d_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_patch, Flux_snd_d_patch,          &
!$OMP&    EntFlux_y, Fluct_y_Pl, Fluct_y_Mi)                            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_u_patch, Flux_snd_d_patch, EntFlux_y,                &
!$acc&    Fluct_y_Pl, Fluct_y_Mi)

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, iyu, iyd )

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice
          iyu = ntp_zones      ! location of up   slice

! ----------
! ip1 == 3 + 1
! ip2 == 3 + n_cc_components

          Flux_snd_u_patch ( ixx,      izz, ip1: ip2) =
     1          Fluct_y_Pl ( ixx, iyu, izz,   1: n_cc_components)

          Flux_snd_d_patch ( ixx,      izz, ip1: ip2) =
     1          Fluct_y_Pl ( ixx, iyd, izz,   1: n_cc_components)

! ----------
! im1 == 3 + n_cc_components + 1
! im2 == 3 + n_cc_components + n_cc_components

          Flux_snd_u_patch ( ixx,      izz, im1: im2) =
     1          Fluct_y_Mi ( ixx, iyu, izz,   1: n_cc_components)

          Flux_snd_d_patch ( ixx,      izz, im1: im2) =
     1          Fluct_y_Mi ( ixx, iyd, izz,   1: n_cc_components)

! ----------

          Flux_snd_u_patch ( ixx,      izz, 1) =
     1           EntFlux_y ( ixx, iyu, izz)

          Flux_snd_d_patch ( ixx,      izz, 1) =
     1           EntFlux_y ( ixx, iyd, izz)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efz" data in "Flux_snd_u_patch, Flux_snd_d_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_patch, Flux_snd_d_patch, Efz )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_u_patch, Flux_snd_d_patch, Efz )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, iyu, iyd )

        DO izz = 1 - NGhost, nr_zones  + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.

          Flux_snd_u_patch ( ixx,      izz, 2) =
     1                 Efz ( ixx, iyu, izz)

          Flux_snd_d_patch ( ixx,      izz, 2) =
     1                 Efz ( ixx, iyd, izz)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 1 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efx" data in "Flux_snd_u_patch, Flux_snd_d_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixs, iyu, iyd )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_u_patch, Flux_snd_d_patch, Efx )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_u_patch, Flux_snd_d_patch, Efx )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, ixs, iyu, iyd )

        DO izz = 0 - NGhost, nr_zones  + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          iyd = 0              ! location of down slice.
          iyu = ntp_zones      ! location of up   slice.

          Flux_snd_u_patch ( ixx,      izz, 3) =
     1                 Efx ( ixx, iyu, izz)

          Flux_snd_d_patch ( ixx,      izz, 3) =
     1                 Efx ( ixx, iyd, izz)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO izz = 0 - NGhost, nr_zones  + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Store top and bottom z--data in "Flux_snd_t_patch, Flux_snd_b_patch".

! ----------------------------------------------------------------------
! Loading "Flux_snd_t_patch, Flux_snd_b_patch". First for the face centered fluxes,
! "z-Fluctuations, EntFlux_z" then for the edge centered variables,
! "Efx, Efy".

! ------------------------------------------------------------
! Storing face centered fluxes "z-Fluctuations, EntFlux_z" data in
! "Flux_snd_t_patch, Flux_snd_b_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_patch, Flux_snd_b_patch,          &
!$OMP&    EntFlux_z, Fluct_z_Pl, Fluct_z_Mi)                            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_t_patch, Flux_snd_b_patch, EntFlux_z,                &
!$acc&    Fluct_z_Pl, Fluct_z_Mi)

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, izs, izt, izb )

        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

! ----------
! ip1 == 3 + 1
! ip2 == 3 + n_cc_components

          Flux_snd_t_patch ( ixx, iyy,      ip1: ip2) =
     1          Fluct_z_Pl ( ixx, iyy, izt,   1: n_cc_components)

          Flux_snd_b_patch ( ixx, iyy,      ip1: ip2) =
     1          Fluct_z_Pl ( ixx, iyy, izb,   1: n_cc_components)

! ----------
! im1 == 3 + n_cc_components + 1
! im2 == 3 + n_cc_components + n_cc_components

          Flux_snd_t_patch ( ixx, iyy,      im1: im2) =
     1          Fluct_z_Mi ( ixx, iyy, izt,   1: n_cc_components)

          Flux_snd_b_patch ( ixx, iyy,      im1: im2) =
     1          Fluct_z_Mi ( ixx, iyy, izb,   1: n_cc_components)

! ----------

          Flux_snd_t_patch ( ixx, iyy,      1) =
     1           EntFlux_z ( ixx, iyy, izt)

          Flux_snd_b_patch ( ixx, iyy,      1) =
     1           EntFlux_z ( ixx, iyy, izb)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efx" data in "Flux_snd_t_patch, Flux_snd_b_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_patch, Flux_snd_b_patch, Efx )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_t_patch, Flux_snd_b_patch, Efx )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, izs, izt, izb )

        DO iyy = 0 - NGhost, ntp_zones + NGhost
        DO ixx = 1 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

          Flux_snd_t_patch ( ixx, iyy,      2) =
     1                 Efx ( ixx, iyy, izt)

          Flux_snd_b_patch ( ixx, iyy,      2) =
     1                 Efx ( ixx, iyy, izb)

! ----------------------------------------

        END DO ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------------------------------------
! Storing edge centered variable "Efy" data in "Flux_snd_t_patch, Flux_snd_b_patch".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, izs, izt, izb )              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components, Flux_snd_t_patch, Flux_snd_b_patch, Efy )    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    Flux_snd_t_patch, Flux_snd_b_patch, Efy )

!$acc loop gang vector collapse(2) independent                          &
!$acc&    private ( ixx, iyy, izz, izs, izt, izb )

        DO iyy = 1 - NGhost, ntp_zones + NGhost
        DO ixx = 0 - NGhost, ntp_zones + NGhost

! ----------------------------------------

          izb = 0               ! location of bottom slice.
          izt = nr_zones        ! location of top slice.

          Flux_snd_t_patch ( ixx, iyy,      3) =
     1                 Efy ( ixx, iyy, izt)

          Flux_snd_b_patch ( ixx, iyy,      3) =
     1                 Efy ( ixx, iyy, izb)

! ----------------------------------------

        END DO ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost".
        END DO ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------

      END SUBROUTINE FILL_FLUX_SND_PATCH


