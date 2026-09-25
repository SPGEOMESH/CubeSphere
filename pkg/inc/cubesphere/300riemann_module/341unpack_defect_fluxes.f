! contained in RIEMANN_3D_CS_AMR_SUPPORT










      SUBROUTINE UNPACK_DEFECT_FLUXES ( )

! ----------------------------------------------------------------------

! Written By: Sethupathy S.

! This subroutine is called before the synchronization of fluxes and
! electric fields at the coinciding faces and edges, respectively.
! For that, the amr fluxes and electric fields at the coinciding faces and edges
! are populated in this subroutine from the "Flux_snd_X_amr" arrays.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, npatch_per_proc,

     1        ixx, iyy, izz, icomp, ii,

     1        ixr, iyr, izr, ixl, iyl, izl,
     1        ixu, iyu, izu, ixd, iyd, izd,
     1        izt, izb,

     1        ixs, iys, izs

! ----------------------------------------------------------------------
! Synchronize on entry.

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Assign Defect fluxes and electric fields at the coinciding faces and edges.

      DO ii = 1, npatch_per_proc

! --------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz )                             &

!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost, ii,            &
!$OMP&    n_cc_components,                                              &
!$OMP&    Defect_Efx_amr, Efx_amr,                                      &
!$OMP&    Defect_Efy_amr, Efy_amr,                                      &
!$OMP&    Defect_Efz_amr, Efz_amr,                                      &
!$OMP&    Defect_EntFlux_x_amr, EntFlux_x_amr,                          &
!$OMP&    Defect_EntFlux_y_amr, EntFlux_y_amr,                          &
!$OMP&    Defect_EntFlux_z_amr, EntFlux_z_amr)                          &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

        DO izz = 0,  nr_zones
        DO iyy = 0, ntp_zones
        DO ixx = 0, ntp_zones

! --------------------
! x-face.

        IF ( ( ixx .EQ. 0) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_EntFlux_x_amr ( ixx, iyy, izz, ii) =
     1  Defect_EntFlux_x_amr ( ixx, iyy, izz, ii)
     1       + EntFlux_x_amr ( ixx, iyy, izz, ii)

        Defect_Fluct_x_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  Defect_Fluct_x_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1       + Fluct_x_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

        Defect_Fluct_x_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  Defect_Fluct_x_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1       + Fluct_x_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

! ----------

        END IF   ! End of "x-face".

! --------------------
! y-face.

        IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 0) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_EntFlux_y_amr ( ixx, iyy, izz, ii) =
     1  Defect_EntFlux_y_amr ( ixx, iyy, izz, ii)
     1       + EntFlux_y_amr ( ixx, iyy, izz, ii)

        Defect_Fluct_y_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  Defect_Fluct_y_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1       + Fluct_y_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

        Defect_Fluct_y_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  Defect_Fluct_y_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1       + Fluct_y_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

! ----------

        END IF   ! End of "y-face".

! --------------------
! z-face.

        IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 0) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_EntFlux_z_amr ( ixx, iyy, izz, ii) =
     1  Defect_EntFlux_z_amr ( ixx, iyy, izz, ii)
     1       + EntFlux_z_amr ( ixx, iyy, izz, ii)

        Defect_Fluct_z_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  Defect_Fluct_z_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1       + Fluct_z_Pl_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

        Defect_Fluct_z_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii) =
     1  Defect_Fluct_z_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)
     1       + Fluct_z_Mi_amr ( ixx, iyy, izz, 1: n_cc_components, ii)

! ----------

        END IF   ! End of "z-face".

! --------------------
! x-edge.

        IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 0) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 0) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Efx_amr ( ixx, iyy, izz, ii) =
     1  Defect_Efx_amr ( ixx, iyy, izz, ii)
     1       + Efx_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "x-edge".

! --------------------
! y-edge.

        IF ( ( ixx .EQ. 0) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 0) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Efy_amr ( ixx, iyy, izz, ii) =
     1  Defect_Efy_amr ( ixx, iyy, izz, ii)
     1       + Efy_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "y-edge".

! --------------------
! z-edge.

        IF ( ( ixx .EQ. 0) .OR. ( ixx .EQ. ntp_zones) .OR.
     1       ( iyy .EQ. 0) .OR. ( iyy .EQ. ntp_zones) .OR.
     1       ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ----------

        Defect_Efz_amr ( ixx, iyy, izz, ii) =
     1  Defect_Efz_amr ( ixx, iyy, izz, ii)
     1       + Efz_amr ( ixx, iyy, izz, ii)

! ----------

        END IF   ! End of "z-edge".

! --------------------

        END DO ! End of "DO ixx = 0, ntp_zones".
        END DO ! End of "DO iyy = 0, ntp_zones".
        END DO ! End of "DO izz = 0,  nr_zones".

!$OMP END PARALLEL DO

! --------------------------------------------------

      END DO ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)
! Synchronize before exit.

      END SUBROUTINE UNPACK_DEFECT_FLUXES




