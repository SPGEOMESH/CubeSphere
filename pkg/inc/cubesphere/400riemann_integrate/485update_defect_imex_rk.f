


#if ( NUMORDER == 2)

      SUBROUTINE UPDATE_DEFECT_IMEX_RK

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_components, n_cc_modes,
     1         n_fc_modes, time, dtmove,

     1         irk,

     1         vert_coords_patch, imask_active_patch,

     1         xf_quadnorml_patch, yf_quadnorml_patch,
     1         zf_quadnorml_patch,
     1         xf_area_patch, yf_area_patch, zf_area_patch,
     1         v_scale_length, v_volume_patch, v_centroid_patch,

     1         U, entropygr, Bxf, Byf, Bzf, dU_dt,

     1         Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1         Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi,

     1         EntFlux_x, EntFlux_y, EntFlux_z,
     1         Efx, Efy, Efz, flag_entprs)

! ----------------------------------------------------------------------
! Written by: Dinshaw Balsara and Sethupathy S.

! This is a defect-corrector for the RK step subroutine.
! It directly operates on the conserved variable "U".
!
! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes

      REAL time, dtmove

      INTEGER irk

! ------------------------------
! VARIABLES THAT GENERATE THE MESH (set in EXTRACT_VERT_COORDS_PATCH):-

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

! ------------------------------
! FACIAL VARIABLES (set in COORDS_TO_WT_QUADPTS_FACES_ISO2):-

      REAL xf_quadnorml_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL yf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL zf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL xf_area_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL yf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL zf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ------------------------------
! VOLUMETRIC VARIABLES (set in COORDS_TO_WT_QUADPTS_VOLUME_ISO2):-

      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------
! FLUID VARIABLES:-

      REAL U ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL entropygr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   4)

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           n_fc_modes) ::
     1   Bxf, Byf, Bzf

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: nr_zones + NGhost,
     1                                     n_cc_components) ::
     1   dU_dt

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost: ntp_zones + NGhost,
     1                  0 - NGhost:  nr_zones + NGhost,
     1                                          n_cc_components) ::
     1      Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1      Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi

! ------------------------------
! FLUX AND ELECTRIC FIELD VARIABLES:-

! Entropy flux from the 1D Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      EntFlux_x, EntFlux_y, EntFlux_z

! Electric fields at the edges of the mesh.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      Efx, Efy, Efz

! Extra variable for silo output.

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost) :: flag_entprs

! ----------------------------------------------------------------------
! Local variables.

      INTEGER ixx, iyy, izz, ixxm1, iyym1, izzm1, j, iexpand,
     1        iprs_flag, bflag, irk_use

      REAL rin, rout, vx1, vy1, vz1, vx2, vy2, vz2,
     1     theta, lin, lout, dr, lxp, lxm, lyp, lym,
     1     tempaa, temp,
     1     rho, prs, vx, vy, vz, bx, by, bz, ob_8pi,

     1      x, y, z, a, b, c, time_use, B0 ( 3),
     1      ent, Upt ( n_cc_components), Vpt ( n_cc_components),
     1      s9, tau_2, G,

! For Bzf calculation.
     1     vrx1, vry1, vrz1, vrx2, vry2, vrz2,
     1     vlx1, vly1, vlz1, vlx2, vly2, vlz2,
     1     vux1, vuy1, vuz1, vux2, vuy2, vuz2,
     1     vdx1, vdy1, vdz1, vdx2, vdy2, vdz2

! ----------------------------------------------------------------------

      tau_2 = tau_phys_cgl

      irk_use = irk

      G = 1.0 - 1.0 / SQRT ( 2.0)

      ob_8pi = 1.0 / ( 8.0 * pi)

! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1)         &
!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,&
!$OMP&    irk_use, tau_2, G, dtmove, U, dU_dt,                          &
!$OMP&    xf_area_patch, yf_area_patch, zf_area_patch, v_volume_patch,  &
!$OMP&    Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,                           &
!$OMP&    Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi)                           &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

      DO izz = 1 - NGhost, nr_zones  + NGhost
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

      IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1     ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1     ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ------------------------------------------------------------
! In the absence of source terms, Defect correction is similar for both rk-stages.

      IF ( ( irk_use .EQ. 1) .OR. ( irk_use .EQ. 2)) THEN

! --------------------------------------------------
! RK-step 1:  ( Defect correction, no source terms)
! ----------
!
!   Un+1 = Un + dt*dU_dt
!
!
! RK-step 2:  ( Defect correction, no source terms)
! ----------
!
!   Un+1 = 0.5*(Un + U2) + 0.5*dt*dU_dt
!
!   0.5*(Un + U2) is pre-calculated and stored in "U" itself outside this subroutine.
!
!   "dtmove" is assigned to be "0.5 * dtcur".
!
!-------------------------------

!$acc loop seq
      DO j = 1, n_cc_components

! --------------------

      ixxm1 = ixx - 1
      iyym1 = iyy - 1
      izzm1 = izz - 1

      dU_dt ( ixx, iyy, izz, j) =
     1                        ( - Fluct_x_Mi ( ixx,   iyy, izz, j)
     1                       * xf_area_patch ( ixx,   iyy, izz)
     1                          - Fluct_x_Pl ( ixxm1, iyy, izz, j)
     1                       * xf_area_patch ( ixxm1, iyy, izz)

     1                          - Fluct_y_Mi ( ixx, iyy,   izz, j)
     1                       * yf_area_patch ( ixx, iyy,   izz)
     1                          - Fluct_y_Pl ( ixx, iyym1, izz, j)
     1                       * yf_area_patch ( ixx, iyym1, izz)

     1                          - Fluct_z_Mi ( ixx, iyy, izz,   j)
     1                       * zf_area_patch ( ixx, iyy, izz)
     1                          - Fluct_z_Pl ( ixx, iyy, izzm1, j)
     1                       * zf_area_patch ( ixx, iyy, izzm1))

     1     / v_volume_patch ( ixx, iyy, izz)


      U ( ixx, iyy, izz, j, 1) = U ( ixx, iyy, izz, j, 1)

     1      + dtmove * dU_dt ( ixx, iyy, izz, j)

! --------------------

      END DO   ! End of "DO j = 1, n_cc_components" loop.

! ----------------------------------------

      END IF  ! End of "IF ( ( irk_use .EQ. 1) .OR. ( irk_use .EQ. 2)) THEN".

! ------------------------------------------------------------

      END IF  ! End of "IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.".

      END DO  ! End of "DO ixx = 1, ntp_zones".
      END DO  ! End of "DO iyy = 1, ntp_zones".
      END DO  ! End of "DO izz = 1, nr_zones".

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! Also update "entropygr" variable with its corresponding defect correction
! flux.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    bflag, time_use, x, y, z, a, b, c, B0,                        &
!$OMP&    iprs_flag, ent, Upt, Vpt)                                     &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, dtmove, ob_8pi, &
!$OMP&    rhofloor, prsfloor, gamma, U, tempfloor, gasconst, molewt,    &
!$OMP&    isospeed, gamma_euler, n_cc_components, barotropic_stiffness, &
!$OMP&    xf_area_patch, yf_area_patch, zf_area_patch, v_volume_patch,  &
!$OMP&    v_centroid_patch, eps_bdry_cgl, magb_small_cgl,               &
!$OMP&    entropygr, EntFlux_x, EntFlux_y, EntFlux_z, flag_entprs)      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

      DO izz = 1, nr_zones

        DO iyy = 1, ntp_zones

          DO ixx = 1, ntp_zones

          IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1         ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1         ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ------------------------------

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

! ------------------------------

          entropygr ( ixx, iyy, izz, 1) = entropygr ( ixx, iyy, izz, 1)

     1 - dtmove * ( xf_area_patch ( ixx, iyy, izz)
     1            * EntFlux_x ( ixx, iyy, izz)
     1            - xf_area_patch ( ixxm1, iyy, izz)
     1            * EntFlux_x ( ixxm1, iyy, izz)

     1            + yf_area_patch ( ixx, iyy, izz)
     1            * EntFlux_y ( ixx, iyy, izz)
     1            - yf_area_patch ( ixx, iyym1, izz)
     1            * EntFlux_y ( ixx, iyym1, izz)

     1            + zf_area_patch ( ixx, iyy, izz)
     1            * EntFlux_z ( ixx, iyy, izz)
     1            - zf_area_patch ( ixx, iyy, izzm1)
     1            * EntFlux_z ( ixx, iyy, izzm1) )

     1 / v_volume_patch ( ixx, iyy, izz)

! ------------------------------
! Compute Primitive variables for "dtcur" calculations and enforce
! appropriate pressure positivity if necessary.

          Upt ( :) = U ( ixx, iyy, izz, :, 1)

          ent = entropygr ( ixx, iyy, izz, 1)

          B0 ( :) = 0.0

#if( USE_B0B1_FORMULATION == 1)
          x = v_centroid_patch ( ixx, iyy, izz, 1)
          y = v_centroid_patch ( ixx, iyy, izz, 2)
          z = v_centroid_patch ( ixx, iyy, izz, 3)

          bflag = 1         ! 1 == magnetic field;  0 == vector potential.
          time_use = 0.0    ! useful for the case of active rotation.

          CALL GET_MAGNETIC_FIELD ( bflag, time_use, x, y, z, a, b, c)

          B0 ( 1) = a
          B0 ( 2) = b
          B0 ( 3) = c
#endif

          CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, Upt, Vpt)

          U ( ixx, iyy, izz, :, 1) = Upt ( :)


! Check if the pressure is assigned from tempfloor or prsfloor or entropy.
! also check if it is not previously flagged in the update step.

          IF ( flag_entprs ( ixx, iyy, izz) .EQ. 0)
     1    flag_entprs ( ixx, iyy, izz) = iprs_flag   ! 1.0 - assigned; 0.0 - orig.

! ------------------------------

          END IF

          END DO  ! End of "DO ixx = 1, ntp_zones" loop.

        END DO  ! End of "DO iyy = 1, ntp_zones" loop.

      END DO  ! End of "DO izz = 1, nr_zones

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! STEP 2) Update "Bxf" at the patch boundary using the defect correction
! field "Efy, Efz". This is done by detecting the boundary faces that
! have the "ixx, iyy, izz" coinciding with "0" or "ntp_zones or nr_zones".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    rin, rout, vx1, vy1, vz1, vx2, vy2, vz2, theta, lin, lout, dr,&
!$OMP&    tempaa)                                                       &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, dtmove, iexpand,&
!$OMP&    vert_coords_patch, xf_area_patch,                             &
!$OMP&    Bxf, Efy, Efz)                                                &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

      DO izz = 1, nr_zones

        DO iyy = 1, ntp_zones

          DO ixx = 0, ntp_zones

          IF ( ( ixx .EQ. 0) .OR. ( ixx .EQ. ntp_zones) .OR.
     1         ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1         ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ------------------------------

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

          rout = SQRT ( vert_coords_patch ( ixx, iyy, izz, 1)**2  ! Radius of outer arc.
     1                + vert_coords_patch ( ixx, iyy, izz, 2)**2
     1                + vert_coords_patch ( ixx, iyy, izz, 3)**2 )

          rin = SQRT ( vert_coords_patch ( ixx, iyy, izzm1, 1)**2  ! Radius of inner arc.
     1               + vert_coords_patch ( ixx, iyy, izzm1, 2)**2
     1               + vert_coords_patch ( ixx, iyy, izzm1, 3)**2 )

          dr = rout - rin  ! Length of x-face in the radial direction.

          tempaa = 1.0 / rout

          vx1 = vert_coords_patch ( ixx, iyy, izz, 1) * tempaa  ! Unit vector in x-face.
          vy1 = vert_coords_patch ( ixx, iyy, izz, 2) * tempaa
          vz1 = vert_coords_patch ( ixx, iyy, izz, 3) * tempaa

          vx2 = vert_coords_patch ( ixx, iyym1, izz, 1) * tempaa  ! Unit vector in x-face.
          vy2 = vert_coords_patch ( ixx, iyym1, izz, 2) * tempaa
          vz2 = vert_coords_patch ( ixx, iyym1, izz, 3) * tempaa

          theta = ACOS ( vx1 * vx2 + vy1 * vy2 + vz1 * vz2 )  ! Dot product; find angle.

          lin = rin * theta  ! Length of inner arc.
          lout = rout * theta  ! Length of outer arc.

! ------------------------------

          Bxf ( ixx, iyy, izz, 1) = Bxf ( ixx, iyy, izz, 1)

     1  - dtmove * ( dr * Efz ( ixx, iyy, izz)
     1             - dr * Efz ( ixx, iyym1, izz)
     1             - lout * Efy ( ixx, iyy, izz)
     1             + lin * Efy ( ixx, iyy, izzm1) )

     1  / xf_area_patch ( ixx, iyy, izz)

! ------------------------------

          END IF

          END DO  ! End "DO ixx = 0, ntp_zones" loop.

        END DO  ! End "DO iyy = 1, ntp_zones" loop.

      END DO  ! End "DO izz = 1, nr_zones" loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! STEP 3) Update "Byf" at the patch boundary using the defect correction
! field "Efx, Efz". This is done by detecting the boundary faces that
! have the "ixx, iyy, izz" coinciding with "0" or "ntp_zones or nr_zones".


!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    rin, rout, vx1, vy1, vz1, vx2, vy2, vz2, theta, lin, lout, dr,&
!$OMP&    tempaa)                                                       &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, dtmove, iexpand,&
!$OMP&    vert_coords_patch, yf_area_patch,                             &
!$OMP&    Byf, Efx, Efz)                                                &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

      DO izz = 1, nr_zones

        DO iyy = 0, ntp_zones

          DO ixx = 1, ntp_zones

          IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1         ( iyy .EQ. 0) .OR. ( iyy .EQ. ntp_zones) .OR.
     1         ( izz .EQ. 1) .OR. ( izz .EQ.  nr_zones) ) THEN

! ------------------------------

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

          rout = SQRT ( vert_coords_patch ( ixx, iyy, izz, 1)**2  ! Radius of outer arc.
     1                + vert_coords_patch ( ixx, iyy, izz, 2)**2
     1                + vert_coords_patch ( ixx, iyy, izz, 3)**2 )

          rin = SQRT ( vert_coords_patch ( ixx, iyy, izzm1, 1)**2  ! Radius of inner arc.
     1               + vert_coords_patch ( ixx, iyy, izzm1, 2)**2
     1               + vert_coords_patch ( ixx, iyy, izzm1, 3)**2 )

          dr = rout - rin  ! Length of y-face in the radial direction.

          tempaa = 1.0 / rout

          vx1 = vert_coords_patch ( ixx, iyy, izz, 1) * tempaa  ! Unit vector in y-face.
          vy1 = vert_coords_patch ( ixx, iyy, izz, 2) * tempaa
          vz1 = vert_coords_patch ( ixx, iyy, izz, 3) * tempaa

          vx2 = vert_coords_patch ( ixxm1, iyy, izz, 1) * tempaa  ! Unit vector in y-face.
          vy2 = vert_coords_patch ( ixxm1, iyy, izz, 2) * tempaa
          vz2 = vert_coords_patch ( ixxm1, iyy, izz, 3) * tempaa

          theta = ACOS ( vx1 * vx2 + vy1 * vy2 + vz1 * vz2 )  ! Dot product; find angle.

          lin = rin * theta  ! Length of inner arc.
          lout = rout * theta  ! Length of outer arc.

! ------------------------------

          Byf ( ixx, iyy, izz, 1) = Byf ( ixx, iyy, izz, 1)

     1  - dtmove * ( - dr * Efz ( ixx, iyy, izz)
     1               + dr * Efz ( ixxm1, iyy, izz)
     1               + lout * Efx ( ixx, iyy, izz)
     1               - lin * Efx ( ixx, iyy, izzm1) )

     1  / yf_area_patch ( ixx, iyy, izz)

! ------------------------------

          END IF

          END DO  ! End "DO ixx = 1, ntp_zones" loop.

        END DO  ! End "DO iyy = 0, ntp_zones" loop.

      END DO  ! End "DO izz = 1, nr_zones" loop.

!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! STEP 3) Update "Bzf" at the patch boundary using the defect correction
! field "Efx, Efy". This is done by detecting the boundary faces that
! have the "ixx, iyy, izz" coinciding with "0" or "ntp_zones or nr_zones".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    rin, rout, lin, lout, dr,                                     &
!$OMP&    lxp, lxm, lyp, lym, tempaa,                                   &
!$OMP&    vrx1, vry1, vrz1, vrx2, vry2, vrz2,                           &
!$OMP&    vlx1, vly1, vlz1, vlx2, vly2, vlz2,                           &
!$OMP&    vux1, vuy1, vuz1, vux2, vuy2, vuz2,                           &
!$OMP&    vdx1, vdy1, vdz1, vdx2, vdy2, vdz2, theta)                    &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, dtmove, iexpand,&
!$OMP&    vert_coords_patch, zf_area_patch,                             &
!$OMP&    Bzf, Efx, Efy)                                                &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

      DO izz = 0, nr_zones

        DO iyy = 1, ntp_zones

          DO ixx = 1, ntp_zones

          IF ( ( ixx .EQ. 1) .OR. ( ixx .EQ. ntp_zones) .OR.
     1         ( iyy .EQ. 1) .OR. ( iyy .EQ. ntp_zones) .OR.
     1         ( izz .EQ. 0) .OR. ( izz .EQ.  nr_zones) ) THEN

! ------------------------------

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1


          rout = SQRT ( vert_coords_patch ( ixx, iyy, izz, 1)**2  ! Radius of outer arc.
     1                + vert_coords_patch ( ixx, iyy, izz, 2)**2
     1                + vert_coords_patch ( ixx, iyy, izz, 3)**2 )

          tempaa = 1.0 / rout

! ----------

! Making lxp :-
          vrx1 = vert_coords_patch ( ixx, iyy, izz, 1) * tempaa  ! Unit vector in +, +.
          vry1 = vert_coords_patch ( ixx, iyy, izz, 2) * tempaa
          vrz1 = vert_coords_patch ( ixx, iyy, izz, 3) * tempaa

          vrx2 = vert_coords_patch ( ixxm1, iyy, izz, 1) * tempaa  ! Unit vector in -, +.
          vry2 = vert_coords_patch ( ixxm1, iyy, izz, 2) * tempaa
          vrz2 = vert_coords_patch ( ixxm1, iyy, izz, 3) * tempaa

          theta = ACOS ( vrx1 * vrx2 + vry1 * vry2 + vrz1 * vrz2 )  ! Dot product; find angle.

          lxp = rout * theta

! ----------

! Making lxm :-
          vlx1 = vert_coords_patch ( ixx, iyym1, izz, 1) * tempaa  ! Unit vector in +, -.
          vly1 = vert_coords_patch ( ixx, iyym1, izz, 2) * tempaa
          vlz1 = vert_coords_patch ( ixx, iyym1, izz, 3) * tempaa

          vlx2 = vert_coords_patch ( ixxm1, iyym1, izz, 1) * tempaa  ! Unit vector in -, -.
          vly2 = vert_coords_patch ( ixxm1, iyym1, izz, 2) * tempaa
          vlz2 = vert_coords_patch ( ixxm1, iyym1, izz, 3) * tempaa

          theta = ACOS ( vlx1 * vlx2 + vly1 * vly2 + vlz1 * vlz2 )  ! Dot product; find angle.
          lxm = rout * theta

! ----------

! Making lyp :-
          vux1 = vert_coords_patch ( ixx, iyy, izz, 1) * tempaa  ! Unit vector in +, +.
          vuy1 = vert_coords_patch ( ixx, iyy, izz, 2) * tempaa
          vuz1 = vert_coords_patch ( ixx, iyy, izz, 3) * tempaa

          vux2 = vert_coords_patch ( ixx, iyym1, izz, 1) * tempaa  ! Unit vector in +, -.
          vuy2 = vert_coords_patch ( ixx, iyym1, izz, 2) * tempaa
          vuz2 = vert_coords_patch ( ixx, iyym1, izz, 3) * tempaa

          theta = ACOS ( vux1 * vux2 + vuy1 * vuy2 + vuz1 * vuz2 )  ! Dot product; find angle.

          lyp = rout * theta

! ----------

! Making lym :-
          vdx1 = vert_coords_patch ( ixxm1, iyy, izz, 1) * tempaa  ! Unit vector in -, +.
          vdy1 = vert_coords_patch ( ixxm1, iyy, izz, 2) * tempaa
          vdz1 = vert_coords_patch ( ixxm1, iyy, izz, 3) * tempaa

          vdx2 = vert_coords_patch ( ixxm1, iyym1, izz, 1) * tempaa  ! Unit vector in -, -.
          vdy2 = vert_coords_patch ( ixxm1, iyym1, izz, 2) * tempaa
          vdz2 = vert_coords_patch ( ixxm1, iyym1, izz, 3) * tempaa

          theta = ACOS ( vdx1 * vdx2 + vdy1 * vdy2 + vdz1 * vdz2 )  ! Dot product; find angle.

          lym = rout * theta

! ------------------------------

          Bzf ( ixx, iyy, izz, 1) = Bzf ( ixx, iyy, izz, 1)

     1  - dtmove * ( lyp * Efy ( ixx, iyy, izz)
     1             - lym * Efy ( ixxm1, iyy, izz)
     1             - lxp * Efx ( ixx, iyy, izz)
     1             + lxm * Efx ( ixx, iyym1, izz) )

     1  / zf_area_patch ( ixx, iyy, izz)

! ------------------------------

          END IF

          END DO  ! End "DO ixx = 1, ntp_zones" loop.

        END DO  ! End "DO iyy = 1, ntp_zones" loop.

      END DO  ! End "DO izz = 0, nr_zones" loop.

!$OMP END PARALLEL DO
! ----------------------------------------------------------------------

      END SUBROUTINE UPDATE_DEFECT_IMEX_RK

#endif
! End of "if ( NUMORDER == 2)".




