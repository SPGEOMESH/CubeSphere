


      SUBROUTINE CLEANUP_AFTER_UPDATE

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           iz_patch, nr_patches,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadnorml_patch, xf_centroid_patch,
     1           yf_quadnorml_patch, yf_centroid_patch,
     1           zf_quadnorml_patch, zf_centroid_patch,

     1           v_scale_length, v_centroid_patch,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! After the call to "UPDATE_IMEX_RK_STEP_P1" we may, in rare cases, have some zones
! with non-positive pressures. We bring them into the physical domain with this step.

! ----------------------------------------------------------------------

! ------------------------------
! CUBE SPHERE MESH VARIABLES:-

! ----------
! VARIABLES THAT GENERATE THE MESH (set in EXTRACT_VERT_COORDS_PATCH):-

! "vert_coords_patch" coordinates of the vertices of each patch. From this,
! all the other patch-based geometry information is generated via quadratic
! isoparametric mapping.

! "imask_active_patch" is 1 for good zones, 0 for bad ones.

! ----------
! FACIAL VARIABLES (set in COORDS_TO_WT_QUADPTS_FACES_ISO2):-

! "xf_quadnorml_patch, yf_quadnorml_patch, zf_quadnorml_patch" unit normal
! vector to that face. The z-face is curved and has a normal at each
! quadrature point.

! "xf_centroid_patch, yf_centroid_patch, zf_centroid_patch" coordinates
! of the centroids of the faces.

! "v_scale_length" a scale length that tells us how big each zone is.

! ------------------------------
! VARIABLES THAT FACILITATE STENCIL OPERATIONS (set in STENCIL_INDEXING_TRICORNERS):-

! "itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU" (output) labels the
! four corners in the xy-plane. 1 if this is a tri-corner patch; 0 if this
! is a normal quad-corner patch.

! "istencil_patch ( :, :, 2, 12)" (output)
! "istencil_patch ( :, :, 2, 1:4)" gives us Right, Left, Up, Down zones in xy-plane.
! "istencil_patch ( :, :, 2, 5:8)" gives us the zones that are two ones over in the
! Right, Left, Up, Down directions.
! "istencil_patch ( :, :, 2, 9:12)" gives us the RU, LU, LD, RD zones.

! ------------------------------
! FLUID VARIABLES:-

! "U" comes in with just the zone-centered conserved variables and leaves with
! their moments included. The velocities are defined in a global Cartesian
! space, as are the magnetic fields.

! "V" is the vector of primitive variables that are built from the conserved
! variables. Useful for building flattener etc. variables.

! "Bxf, Byf, Bzf" are magnetic field components in the x,y,z-faces of the mesh.
! Includes moments in the face, that are filled in via reconstruction.

! "B1_zc, B2_zc, B3_zc" are the zone-centered, divergence-free reconstructed
! magnetic fields. They are defined in a global Cartesian space.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes,

     1        iz_patch, nr_patches,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU

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

      REAL xf_centroid_patch ( 0 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

! -----

      REAL yf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL yf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

! -----

      REAL zf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL zf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost, 3)

! ------------------------------

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ------------------------------
! VOLUMETRIC VARIABLES (set in COORDS_TO_WT_QUADPTS_VOLUME_ISO2):-

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------
! VARIABLES THAT FACILITATE STENCIL OPERATIONS (set in STENCIL_INDEXING_TRICORNERS):-

      INTEGER istencil_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                         1 - NGhost: ntp_zones + NGhost, 2, 12)

! ------------------------------
! FLUID VARIABLES:-

      REAL U ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL entropygr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   4)

      REAL V ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           n_fc_modes) ::
     1   Bxf, Byf, Bzf

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes) ::
     1   B1_zc, B2_zc, B3_zc

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost) ::
     1      gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr

! ----------------------------------------------------------------------

! Variables past this point are internal to the subroutine.

      INTEGER ixx, iyy, izz, ixxm1, iyym1, izzm1,
     1        ixxp1, iyyp1, izzp1, pow, icomp, irow, irr,

     1        ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,
     1        ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,

     1        i, j, nn, ii

      REAL ob_8pi, ob_4pi, central_weight,
     1     prs_ratio_min, prs_ratio_max, tinynum, ratio_pres_pos

! ----------------------------------------------------------------------
! Variables used for in-place, second order, WENO reconstruction.

      REAL x, y, z, x2, y2, z2, obl0, obl02, nx, ny, nz,
     1     bx, by, bz, bx2, by2, bz2, bdiff,
     1     a11, a12, a13, a21, a22, a23, a31, a32, a33, obDet,
     1     b11, b12, b13, b21, b22, b23, b31, b32, b33

! ----------------------------------------------------------------------

      INTEGER iter, ipass, inode, bflag, iprs_flag

      REAL, DIMENSION ( n_cc_components) :: V_min, V_max, V_diff

      REAL v_cen ( 3), ent_min, ent_max, ent_diff, ff, omff, expandfac,
     1     aa, bb, cc, gamma_eos_pt,
     1     vx, vy, vz, B0 ( 3), l0,

     1     Upt ( n_cc_components),
     1     Vpt ( n_cc_components), ent,

     1     time_use

! ----------------------------------------------------------------------
! Variables used for linear algebra that is needed in the
! zone-centered assignment of magnetic fields from the facial values.

      REAL nv ( 3),
     1     Am63 ( 6, 3),  Rm6 ( 6), Am33 ( 3, 3),
     1     ATm36 ( 3, 6), Rm3 ( 3), det, R1, R2, R3

      REAL A ( 5, 5), AInv ( 5, 5), R ( 5), Amat_5 ( 5, 5), R5 ( 5)

      INTEGER indx_B ( 5)
      REAL vv_B ( 5)

! ----------------------------------------------------------------------
! Index varibles, which helps in reconstruction using "Coef_Basis".
! It is mainly useful for the reconstruction of B1, B2, B3 where the
! sequence changes for higher order modes.

      INTEGER, PARAMETER ::
! 0, 1st, 2nd order modes.
     1      iC0  = 1,
     1      iCx  = 2,  iCy  = 3,  iCz  = 4


! ----------------------------------------------------------------------

!$acc routine( MATSOLVE ) seq
!$acc routine( GET_MAGNETIC_FIELD ) seq
!$acc routine( CONS_TO_PRIM ) seq

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 1) Initialize various coefficients.

      ob_8pi = 1.0 / ( 8.0 * pi)
      ob_4pi = 1.0 / ( 4.0 * pi)

! Weight of central zone relative to the weight of the other zones (unity).

      central_weight = 100.0

! The power to which the smoothness indicators are raised in the WENO
! reconstruction.

      pow = 4

! Set the minimum and maximum pressure ratios. These govern the settings
! for the scalar version of the flattener. I.e., there is a version of the
! flattener that does not rely on the divergence of the velocity; but only
! relies on the pressure ratio.

      prs_ratio_min = 1.4
      prs_ratio_max = 4.0

! Tiny number to prevent divide by zero specifically in the matrix inversion.

      tinynum = 1.0e-32

! ------------------------------
! Set the factor by which we want to expand the range of each variable.
! Setting this to 0 gives us the TVD range in multi-dimensions.
      expandfac = 0.2

! ------------------------------

! Move above initialization to GPU.

!$acc data copyin ( ob_8pi, ob_4pi, central_weight, pow, expandfac,     &
!$acc&              prs_ratio_min, prs_ratio_max, tinynum,              &
!$acc&              ratio_pres_pos, icomp)

!$acc update device ( itri_patch_LD, itri_patch_RD, itri_patch_RU,      &
!$acc&                itri_patch_LU)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 2) For this very simple situation, we can initialize certain
! obvious scalars. Also zero out reconstructed variables. Make everything
! safe if we have a tri-corner.
! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    gamma, molewt, rhofloor, prsfloor, smallnum,                  &
!$OMP&    imask_active_patch, gamma_eos_gr, gamma_soundspeed_gr,        &
!$OMP&    mean_mol_wt_gr, U, V, entropygr, v_scale_length)              &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)


!$acc parallel vector_length(64) present(                               &
!$acc&    gamma, molewt, rhofloor, prsfloor, smallnum,                  &
!$acc&    imask_active_patch, gamma_eos_gr, gamma_soundspeed_gr,        &
!$acc&    mean_mol_wt_gr, U, V, entropygr, v_scale_length)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)


      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          gamma_eos_gr ( ixx, iyy, izz) = gamma
          gamma_soundspeed_gr ( ixx, iyy, izz) = gamma
          mean_mol_wt_gr ( ixx, iyy, izz) = molewt

          U ( ixx, iyy, izz, :, 2: ) = 0.0  ! Zero out any higher order zonal modes.
          entropygr ( ixx, iyy, izz, 2: ) = 0.0

          V ( ixx, iyy, izz, :, :) = 0.0

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) THEN
! If we get here, it means that this zone is in a tri-corner. Make it safe
! for nominal computation, i.e. prevent divide by zero.

          U ( ixx, iyy, izz, :, 1) = 0.0
          U ( ixx, iyy, izz, 1, 1) = rhofloor
          U ( ixx, iyy, izz, 5, 1) = prsfloor / ( gamma - 1.0)

          v_scale_length ( ixx, iyy, izz) = smallnum

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 0) THEN".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------
! Zero out higher order modes of "Bxf". If we are in a tri-corner, we zero
! out the entire "Bxf".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ii, ixxp1)                   &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    imask_active_patch, Bxf)                                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)


!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, Bxf)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ii, ixxp1)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 0 - NGhost, ntp_zones + NGhost

          ii = MAX0 ( ixx, 1 - NGhost)
          ixxp1 = MIN0 ( ixx + 1, ntp_zones + NGhost)

!          Bxf ( ixx, iyy, izz, 2:) = 0.0  ! Zero out any higher order facial modes.

          IF ( ( imask_active_patch ( ii, iyy, izz) .EQ. 0) .AND.
     1         ( imask_active_patch ( ixxp1, iyy, izz) .EQ. 0) )
     1    Bxf ( ixx, iyy, izz, 1) = 0.0

          END DO  ! End of "DO ixx = 0 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------
! Zero out higher order modes of "Byf". If we are in a tri-corner, we zero
! out the entire "Byf".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ii, iyyp1)                   &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    imask_active_patch, Byf)                                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)


!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, Byf)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ii, iyyp1)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 0 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          ii = MAX0 ( iyy, 1 - NGhost)
          iyyp1 = MIN0 ( iyy + 1, ntp_zones + NGhost)

!          Byf ( ixx, iyy, izz, 2:) = 0.0  ! Zero out any higher order facial modes.

          IF ( ( imask_active_patch ( ixx, ii, izz) .EQ. 0) .AND.
     1         ( imask_active_patch ( ixx, iyyp1, izz) .EQ. 0) )
     1    Byf ( ixx, iyy, izz, 1) = 0.0

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 0 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ------------------------------
! Zero out higher order modes of "Bzf". If we are in a tri-corner, we zero
! out the entire "Bzf".

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &

!$OMP&    SHARED ( nr_zones, ntp_zones, NGhost,                         &
!$OMP&    imask_active_patch, Bzf)                                      &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)


!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, Bzf)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 0 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

!          Bzf ( ixx, iyy, izz, 2:) = 0.0  ! Zero out any higher order facial modes.

! zones are masked along the "x" and "y" directions only and this is
! independent of the radial direction - "izz".
          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 0)
     1    Bzf ( ixx, iyy, izz, 1) = 0.0

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 0 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 3) If any of "itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU"
! are active, we have to synchronize the facial magnetic fields at the matching
! faces.
! ----------------------------------------------------------------------
! At this corner, we loop along the x-direction. The y-magnetic field components at
! the bottom of the zone should be synchronized with the field from the
! Down zone. Specifically, to the lower x-magnetic field components in the Down zone.
! The two fields will have opposite signs, so that has to be accounted for in
! the averaging that is needed for synchronization.

! This step is very intricate and easy to get wrong. So we don't reorder any
! loops here.

      IF ( itri_patch_LD .EQ. 1) THEN

      iyy = 1

!$OMP PARALLEL DO PRIVATE ( ixx, izz, ixxd, iyyd, a11)                  &

!$OMP&            SHARED ( iyy, nr_zones, ntp_zones, NGhost,            &
!$OMP&    imask_active_patch, istencil_patch, Bxf, Byf)                 &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, istencil_patch, Bxf, Byf)

!$acc loop gang independent                                             &
!$acc&    private ( ixx, izz, ixxd, iyyd, a11)

      DO ixx = 1 - NGhost, 0

      ixxd = istencil_patch ( ixx, iyy, 1, 4)
      iyyd = istencil_patch ( ixx, iyy, 2, 4)

!$acc loop vector
        DO izz = 1 - NGhost, nr_zones + NGhost

        a11 = 0.5 * ( Byf ( ixx, iyy - 1, izz, 1)
     1              - Bxf ( ixxd - 1, iyyd, izz, 1) )

        Byf ( ixx, iyy - 1, izz, 1) = a11

        Bxf ( ixxd - 1, iyyd, izz, 1) = - a11

        END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

      END DO  ! End of "DO ixx = 1 - NGhost, 0".

!$acc end parallel
!$OMP END PARALLEL DO

      END IF  ! End of "IF ( itri_patch_LD .EQ. 1) THEN".

! ------------------------------------------------------------
! At this corner, we loop along the x-direction. The y-magnetic field components at
! the bottom of the zone should be synchronized with the field from the
! Down zone. Specifically, to the upper x-magnetic field components in the Down zone.
! The two fields will have same signs, so that has to be accounted for in
! the averaging that is needed for synchronization.

      IF ( itri_patch_RD .EQ. 1) THEN

      iyy = 1

!$OMP PARALLEL DO PRIVATE ( ixx, izz, ixxd, iyyd, a11)                  &

!$OMP&            SHARED ( iyy, nr_zones, ntp_zones, NGhost,            &
!$OMP&    imask_active_patch, istencil_patch, Bxf, Byf)                 &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, istencil_patch, Bxf, Byf)

!$acc loop gang independent                                             &
!$acc&    private ( ixx, izz, ixxd, iyyd, a11)

      DO ixx = ntp_zones + 1, ntp_zones + NGhost

      ixxd = istencil_patch ( ixx, iyy, 1, 4)
      iyyd = istencil_patch ( ixx, iyy, 2, 4)

!$acc loop vector
        DO izz = 1 - NGhost, nr_zones + NGhost

        a11 = 0.5 * ( Byf ( ixx, iyy - 1, izz, 1)
     1              + Bxf ( ixxd, iyyd, izz, 1) )

        Byf ( ixx, iyy - 1, izz, 1) = a11

        Bxf ( ixxd, iyyd, izz, 1) = a11

        END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

      END DO  ! End of "DO ixx = ntp_zones + 1, ntp_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END IF  ! End of "IF ( itri_patch_RD .EQ. 1) THEN".

! ------------------------------------------------------------
! At this corner, we loop along the x-direction. The y-magnetic field components at
! the top of the zone should be synchronized with the field from the
! Up zone. Specifically, to the upper x-magnetic field components in the Up zone.
! The two fields will have opposite signs, so that has to be accounted for in
! the averaging that is needed for synchronization.

      IF ( itri_patch_RU .EQ. 1) THEN

      iyy = ntp_zones

!$OMP PARALLEL DO PRIVATE ( ixx, izz, ixxu, iyyu, a11)                  &

!$OMP&            SHARED ( iyy, nr_zones, ntp_zones, NGhost,            &
!$OMP&    imask_active_patch, istencil_patch, Bxf, Byf)                 &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, istencil_patch, Bxf, Byf)

!$acc loop gang independent                                             &
!$acc&    private ( ixx, izz, ixxu, iyyu, a11)

      DO ixx = ntp_zones + 1, ntp_zones + NGhost

      ixxu = istencil_patch ( ixx, iyy, 1, 3)
      iyyu = istencil_patch ( ixx, iyy, 2, 3)

!$acc loop vector
        DO izz = 1 - NGhost, nr_zones + NGhost

        a11 = 0.5 * ( Byf ( ixx, iyy, izz, 1)
     1              - Bxf ( ixxu, iyyu, izz, 1) )

        Byf ( ixx, iyy, izz, 1) = a11

        Bxf ( ixxu, iyyu, izz, 1) = - a11

        END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

      END DO  ! End of "DO ixx = ntp_zones + 1, ntp_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END IF  ! End of "IF ( itri_patch_RU .EQ. 1) THEN".

! ------------------------------------------------------------
! At this corner, we loop along the x-direction. The y-magnetic field components at
! the top of the zone should be synchronized with the field from the
! Up zone. Specifically, to the lower x-magnetic field components in the Up zone.
! The two fields will have same signs, so that has to be accounted for in
! the averaging that is needed for synchronization.

      IF ( itri_patch_LU .EQ. 1) THEN

      iyy = ntp_zones

!$OMP PARALLEL DO PRIVATE ( ixx, izz, ixxu, iyyu, a11)                  &

!$OMP&            SHARED ( iyy, nr_zones, ntp_zones, NGhost,            &
!$OMP&    imask_active_patch, istencil_patch, Bxf, Byf)                 &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE)

!$acc parallel vector_length(64) present(                               &
!$acc&    imask_active_patch, istencil_patch, Bxf, Byf)

!$acc loop gang independent                                             &
!$acc&    private ( ixx, izz, ixxu, iyyu, a11)

      DO ixx = 1 - NGhost, 0

      ixxu = istencil_patch ( ixx, iyy, 1, 3)
      iyyu = istencil_patch ( ixx, iyy, 2, 3)

!$acc loop vector
        DO izz = 1 - NGhost, nr_zones + NGhost

        a11 = 0.5 * ( Byf ( ixx, iyy, izz, 1)
     1              + Bxf ( ixxu - 1, iyyu, izz, 1) )

        Byf ( ixx, iyy, izz, 1) = a11

        Bxf ( ixxu - 1, iyyu, izz, 1) = a11

        END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

      END DO  ! End of "DO ixx = 1 - NGhost, 0" loop.

!$acc end parallel
!$OMP END PARALLEL DO

      END IF  ! End of "IF ( itri_patch_LU .EQ. 1) THEN".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 4.a) Build estimate for the zone-centered magnetic field from the
! facial magnetic field values using the least squares minimization procedure.

! Since a zone-centered magnetic field "U ( ixx, iyy, izz, 6:8, 1)" is already
! available, for mild problems, this step may not be truly necessary.

! ----------------------------------------------------------------------
#if( 1 == 1)
! This is the better choice. It is exclusively based on least squares.
! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1,         &
!$OMP&    ii, obl0, nv,                                                 &
!$OMP&    a11, a12, a13, a21, a22, a23, a31, a32, a33, obDet,           &
!$OMP&    b11, b12, b13, b21, b22, b23, b31, b32, b33,                  &
!$OMP&    det, Bx, By, Bz, R1, R2, R3,                                  &
!$OMP&    Am63, Rm6, ATm36, Am33, Rm3)                                  &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    ob_8pi, ob_4pi, smallnum, tinynum,                            &
!$OMP&    imask_active_patch, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc, U, V, &
!$OMP&    v_scale_length,                                               &
!$OMP&    xf_quadnorml_patch, yf_quadnorml_patch, zf_centroid_patch)    &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    ob_8pi, ob_4pi, smallnum, tinynum,                            &
!$acc&    imask_active_patch, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc, U, V, &
!$acc&    v_scale_length,                                               &
!$acc&    xf_quadnorml_patch, yf_quadnorml_patch, zf_centroid_patch)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1,                 &
!$acc&    ii, obl0, nv,                                                 &
!$acc&    a11, a12, a13, a21, a22, a23, a31, a32, a33, obDet,           &
!$acc&    b11, b12, b13, b21, b22, b23, b31, b32, b33,                  &
!$acc&    det, Bx, By, Bz, R1, R2, R3,                                  &
!$acc&    Am63, Rm6, ATm36, Am33, Rm3)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN

! ----------------------------------------------------------------------

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

          ii = ixx - (1 - NGhost) + 1

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

! ------------------------------
! Unit normals for the 6-x,y,z faces.
! x- and y- face unit normals are available from "x,yf_quadnorml_patch" array.

! z-face unit normal is calculated from radial unit vector.

          nv ( 1) = zf_centroid_patch ( ixx, iyy, izz, 1)
          nv ( 2) = zf_centroid_patch ( ixx, iyy, izz, 2)
          nv ( 3) = zf_centroid_patch ( ixx, iyy, izz, 3)

          a11 = 1.0 / AMAX1( smallnum, SQRT ( nv ( 1)**2 + nv ( 2)**2
     1                                                   + nv ( 3)**2))

          nv ( :) = nv ( :) * a11

! ----------

          Am63 ( 1, :) = xf_quadnorml_patch ( ixx  , iyy, izz, :)  ! x+ face
          Am63 ( 2, :) = xf_quadnorml_patch ( ixxm1, iyy, izz, :)  ! x- face

          Am63 ( 3, :) = yf_quadnorml_patch ( ixx, iyy  , izz, :)  ! y+ face
          Am63 ( 4, :) = yf_quadnorml_patch ( ixx, iyym1, izz, :)  ! y- face

          Am63 ( 5, :) = nv ( :)    ! z+ face
          Am63 ( 6, :) = nv ( :)    ! z- face

! ------------------------------
! RHS array "R6" consisting 6-facial magnetic field.

          Rm6 ( 1) = Bxf ( ixx  , iyy, izz, 1)
          Rm6 ( 2) = Bxf ( ixxm1, iyy, izz, 1)

          Rm6 ( 3) = Byf ( ixx, iyy  , izz, 1)
          Rm6 ( 4) = Byf ( ixx, iyym1, izz, 1)

          Rm6 ( 5) = Bzf ( ixx, iyy, izz  , 1)
          Rm6 ( 6) = Bzf ( ixx, iyy, izzm1, 1)


! ------------------------------
! Build final matrix and solve.

          ATm36 = TRANSPOSE ( Am63)

          Am33 = MATMUL ( ATm36, Am63)
          Rm3  = MATMUL ( ATm36, Rm6)

          a11 = Am33 ( 1, 1);   a12 = Am33 ( 1, 2);   a13 = Am33 ( 1, 3)
          a21 = Am33 ( 2, 1);   a22 = Am33 ( 2, 2);   a23 = Am33 ( 2, 3)
          a31 = Am33 ( 3, 1);   a32 = Am33 ( 3, 2);   a33 = Am33 ( 3, 3)

          R1 = Rm3 ( 1)
          R2 = Rm3 ( 2)
          R3 = Rm3 ( 3)

! ----------
! Obtain inverse of A33.

          det = a11 * ( a22 * a33 - a32 * a23)
     1        - a12 * ( a21 * a33 - a31 * a23)
     1        + a13 * ( a21 * a32 - a31 * a22)

          obdet = SIGN ( 1.0, det) / AMAX1 ( tinynum, ABS ( det) )

          b11 = ( a22 * a33 - a23 * a32) * obdet
          b12 = ( a13 * a32 - a12 * a33) * obdet
          b13 = ( a12 * a23 - a22 * a13) * obdet

          b21 = ( a23 * a31 - a21 * a33) * obdet
          b22 = ( a11 * a33 - a13 * a31) * obdet
          b23 = ( a13 * a21 - a11 * a23) * obdet

          b31 = ( a21 * a32 - a22 * a31) * obdet
          b32 = ( a12 * a31 - a11 * a32) * obdet
          b33 = ( a11 * a22 - a12 * a21) * obdet

! ------------------------------

          Bx = U ( ixx, iyy, izz, 6, 1)
          By = U ( ixx, iyy, izz, 7, 1)
          Bz = U ( ixx, iyy, izz, 8, 1)

          U ( ixx, iyy, izz, 6, 1) = b11 * R1 + b12 * R2 + b13 * R3
          U ( ixx, iyy, izz, 7, 1) = b21 * R1 + b22 * R2 + b23 * R3
          U ( ixx, iyy, izz, 8, 1) = b31 * R1 + b32 * R2 + b33 * R3

! ------------------------------
! Assign new magentic field from "Bn_zc"; this is also 2nd order accurate
! and consistent. But it violates energy conservation and should only be
! used in strong magnetic field situations where the pressure floor can be hit.

#if ( 1 == 1)
          U ( ixx, iyy, izz, 5, 1) = U ( ixx, iyy, izz, 5, 1)

     1        - ( Bx**2 + By**2 + Bz**2) * ob_8pi

     1        + ( U ( ixx, iyy, izz, 6, 1)**2
     1          + U ( ixx, iyy, izz, 7, 1)**2
     1          + U ( ixx, iyy, izz, 8, 1)**2) * ob_8pi
#endif

! ------------------------------

          B1_zc ( ixx, iyy, izz, 1) = U ( ixx, iyy, izz, 6, 1)
          B2_zc ( ixx, iyy, izz, 1) = U ( ixx, iyy, izz, 7, 1)
          B3_zc ( ixx, iyy, izz, 1) = U ( ixx, iyy, izz, 8, 1)

          V ( ixx, iyy, izz, 6, 1) = U ( ixx, iyy, izz, 6, 1)
          V ( ixx, iyy, izz, 7, 1) = U ( ixx, iyy, izz, 7, 1)
          V ( ixx, iyy, izz, 8, 1) = U ( ixx, iyy, izz, 8, 1)

! ----------------------------------------------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
#else
! Also ok, but found to be not equally stable.
! ----------------------------------------------------------------------
! STEP 4.b) Build estimate for the zone-centered magnetic field that is at
! least second order accurate. This is done by requiring that a constrained
! div-free field matches the mean values at 5 out of the 6 faces of each zone.
! ------------------------------------------------------------

! We will be solving for a lot of 5X5 linear systems here.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, i, j, ii, nn,                &
!$OMP&    x, y, z, obl0, nx, ny, nz, a11, ixxm1, iyym1, izzm1, nt, nb,  &
!$OMP&    A, R, indx_B, vv_B, Amat_5, R5, vx, vy, vz, Bx, By, Bz)       &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost,                 &
!$OMP&    smallnum, tinynum, ob_8pi, rhofloor,                          &
!$OMP&    imask_active_patch, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc, U, V, &
!$OMP&    xf_quadnorml_patch, yf_quadnorml_patch,                       &
!$OMP&    xf_centroid_patch, yf_centroid_patch, zf_centroid_patch,      &
!$OMP&    v_centroid_patch, v_scale_length )                            &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    smallnum, tinynum, ob_8pi, rhofloor,                          &
!$acc&    imask_active_patch, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc, U, V, &
!$acc&    xf_quadnorml_patch, yf_quadnorml_patch,                       &
!$acc&    xf_centroid_patch, yf_centroid_patch, zf_centroid_patch,      &
!$acc&    v_centroid_patch, v_scale_length )

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, i, j, ii, nn,                        &
!$acc&    x, y, z, obl0, nx, ny, nz, a11, ixxm1, iyym1, izzm1, nt, nb,  &
!$acc&    A, R, indx_B, vv_B, Amat_5, R5, vx, vy, vz, Bx, By, Bz)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

          ii = ixx - (1 - NGhost) + 1
! ------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

! ----------
! Face 1, top x-face with facial centroid "xf_centroid_patch ( ixx, iyy, izz, :)".

! Build scaled relative coordinates from zone centroid to facial centroid.

          x = ( xf_centroid_patch ( ixx, iyy, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xf_centroid_patch ( ixx, iyy, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xf_centroid_patch ( ixx, iyy, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Obtain the normal to that face.

          nx = xf_quadnorml_patch ( ixx, iyy, izz, 1)
          ny = xf_quadnorml_patch ( ixx, iyy, izz, 2)
          nz = xf_quadnorml_patch ( ixx, iyy, izz, 3)

! Assemble matrix elements for row 1 of the matrix.

          A ( 1, 1) = nx
          A ( 1, 2) = ny
          A ( 1, 3) = ny * y - nx * x
          A ( 1, 4) = nz
          A ( 1, 5) = nz * z - nx * x

! Assemble RHS for row 1

          R ( 1) = Bxf ( ixx, iyy, izz, 1)

! ----------
! Face 2, bottom x-face with facial centroid "xf_centroid_patch ( ixxm1, iyy, izz, :)".

! Build scaled relative coordinates from zone centroid to facial centroid.

          x = ( xf_centroid_patch ( ixxm1, iyy, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( xf_centroid_patch ( ixxm1, iyy, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( xf_centroid_patch ( ixxm1, iyy, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Obtain the normal to that face.

          nx = xf_quadnorml_patch ( ixxm1, iyy, izz, 1)
          ny = xf_quadnorml_patch ( ixxm1, iyy, izz, 2)
          nz = xf_quadnorml_patch ( ixxm1, iyy, izz, 3)

! Assemble matrix elements for row 2 of the matrix.

          A ( 2, 1) = nx
          A ( 2, 2) = ny
          A ( 2, 3) = ny * y - nx * x
          A ( 2, 4) = nz
          A ( 2, 5) = nz * z - nx * x

! Assemble RHS for row 2

          R ( 2) = Bxf ( ixxm1, iyy, izz, 1)

! ----------
! Face 3, top y-face with facial centroid "yf_centroid_patch ( ixx, iyy, izz, :)".

! Build scaled relative coordinates from zone centroid to facial centroid.

          x = ( yf_centroid_patch ( ixx, iyy, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yf_centroid_patch ( ixx, iyy, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yf_centroid_patch ( ixx, iyy, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Obtain the normal to that face.

          nx = yf_quadnorml_patch ( ixx, iyy, izz, 1)
          ny = yf_quadnorml_patch ( ixx, iyy, izz, 2)
          nz = yf_quadnorml_patch ( ixx, iyy, izz, 3)

! Assemble matrix elements for row 3 of the matrix.

          A ( 3, 1) = nx
          A ( 3, 2) = ny
          A ( 3, 3) = ny * y - nx * x
          A ( 3, 4) = nz
          A ( 3, 5) = nz * z - nx * x

! Assemble RHS for row 3

          R ( 3) = Byf ( ixx, iyy, izz, 1)

! ----------
! Face 4, bottom y-face with facial centroid "yf_centroid_patch ( ixx, iyym1, izz, :)".

! Build scaled relative coordinates from zone centroid to facial centroid.

          x = ( yf_centroid_patch ( ixx, iyym1, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( yf_centroid_patch ( ixx, iyym1, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( yf_centroid_patch ( ixx, iyym1, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Obtain the normal to that face.

          nx = yf_quadnorml_patch ( ixx, iyym1, izz, 1)
          ny = yf_quadnorml_patch ( ixx, iyym1, izz, 2)
          nz = yf_quadnorml_patch ( ixx, iyym1, izz, 3)

! Assemble matrix elements for row 4 of the matrix.

          A ( 4, 1) = nx
          A ( 4, 2) = ny
          A ( 4, 3) = ny * y - nx * x
          A ( 4, 4) = nz
          A ( 4, 5) = nz * z - nx * x

! Assemble RHS for row 4

          R ( 4) = Byf ( ixx, iyym1, izz, 1)

! ----------
! Face 5, top z-face with facial centroid "zf_centroid_patch ( ixx, iyy, izz, :)".

! Build scaled relative coordinates from zone centroid to facial centroid.

          x = ( zf_centroid_patch ( ixx, iyy, izz, 1)
     1      - v_centroid_patch ( ixx, iyy, izz, 1) ) * obl0

          y = ( zf_centroid_patch ( ixx, iyy, izz, 2)
     1      - v_centroid_patch ( ixx, iyy, izz, 2) ) * obl0

          z = ( zf_centroid_patch ( ixx, iyy, izz, 3)
     1      - v_centroid_patch ( ixx, iyy, izz, 3) ) * obl0

! Obtain the normal to that face. Basically, it is a unit normal in the radial
! direction for that curved face.

          nx = zf_centroid_patch ( ixx, iyy, izz, 1)
          ny = zf_centroid_patch ( ixx, iyy, izz, 2)
          nz = zf_centroid_patch ( ixx, iyy, izz, 3)

          a11 = 1.0 / AMAX1 ( smallnum, SQRT ( nx**2 + ny**2 + nz**2 ) )

          nx = nx * a11
          ny = ny * a11
          nz = nz * a11

! Assemble matrix elements for row 5 of the matrix.

          A ( 5, 1) = nx
          A ( 5, 2) = ny
          A ( 5, 3) = ny * y - nx * x
          A ( 5, 4) = nz
          A ( 5, 5) = nz * z - nx * x

! Assemble RHS for row 5

          R ( 5) = Bzf ( ixx, iyy, izz, 1)

! ------------------------------
! Now that "A" and "R" have been assembled, we can inline the linear system
! solve here. The linear solve is "A . X = R". "X" will be returned in "R".
! The explicit subroutine call has also been documented here.

! Even in tri-corners this logic gives us computable values.
          B1_zc ( ixx, iyy, izz, :) = 0.0
          B2_zc ( ixx, iyy, izz, :) = 0.0
          B3_zc ( ixx, iyy, izz, :) = 0.0

          IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

            Amat_5 ( 1:5, 1:5) = A ( 1:5, 1:5)

            R5 ( 1:5) = R ( 1:5)

            CALL MATSOLVE ( 5, Amat_5, R5, indx_B, vv_B)

            R ( 1:5) = R5 ( 1:5)

! ------------------------------
! In this, we scatter back the results from the linear system solve.
! We are only interested in the mean values, not the linear variation of
! Bx in the x-direction; By in the y-direction or Bz in the z-direction.

! ------------------------------
! We have made a quick and easy zone-centered B-field value which is
! at least second order accurate. We reset zone-centered magnetic field
! values here, if desired.

! Magnetic fields at zone centers evaluated consistent with facial
! values and div-free condition.

          B1_zc ( ixx, iyy, izz, 1) = R ( 1)
          B2_zc ( ixx, iyy, izz, 1) = R ( 2)
          B3_zc ( ixx, iyy, izz, 1) = R ( 4)

          B1_zc ( ixx, iyy, izz, 2) = - R ( 3) - R ( 5)
          B2_zc ( ixx, iyy, izz, 3) = R ( 3)
          B3_zc ( ixx, iyy, izz, 4) = R ( 5)

! ------------------------------
! Assign new magentic field from "Bn_zc"; this is also 2nd order accurate
! and consistent. But it violates energy conservation and should only be
! used in strong magnetic field situations where the pressure floor can be hit.

#if ( 1 == 1)
          U ( ixx, iyy, izz, 5, 1) = U ( ixx, iyy, izz, 5, 1)

     1        + ( B1_zc ( ixx, iyy, izz, 1)**2
     1          + B2_zc ( ixx, iyy, izz, 1)**2
     1          + B3_zc ( ixx, iyy, izz, 1)**2) * ob_8pi

     1        - ( U ( ixx, iyy, izz, 6, 1)**2
     1          + U ( ixx, iyy, izz, 7, 1)**2
     1          + U ( ixx, iyy, izz, 8, 1)**2) * ob_8pi
#endif

! ----------

          U ( ixx, iyy, izz, 6, 1) = B1_zc ( ixx, iyy, izz, 1)
          U ( ixx, iyy, izz, 7, 1) = B2_zc ( ixx, iyy, izz, 1)
          U ( ixx, iyy, izz, 8, 1) = B3_zc ( ixx, iyy, izz, 1)

          V ( ixx, iyy, izz, 6, 1) = B1_zc ( ixx, iyy, izz, 1)
          V ( ixx, iyy, izz, 7, 1) = B2_zc ( ixx, iyy, izz, 1)
          V ( ixx, iyy, izz, 8, 1) = B3_zc ( ixx, iyy, izz, 1)

! ------------------------------

          END IF  ! End of "IF ( imask_active_patch ( ixx, iyy, 1) .EQ. 1) THEN".

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

! --------------------------------------------------

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------
#endif
! End of choosing one or the other way of building zone-centered magnetic field.
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! STEP 5) From the conserved variables, we obtain the zone-centered primitive
! variables "V" that are at least second order accurate.
! ----------------------------------------------------------------------

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, iyym1, izzm1, B0,     &
!$OMP&    bflag, time_use,                                              &
!$OMP&    a11, a12, a13, a21, x, y, z, obl0, iprs_flag, ent, Upt, Vpt)  &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,&
!$OMP&    rhofloor, prsfloor, tempfloor, smallnum, ob_4pi, ob_8pi,      &
!$OMP&    gasconst, molewt, gamma, magB_small_cgl, eps_bdry_cgl,        &
!$OMP&    isospeed, gamma_euler, barotropic_stiffness,                  &
!$OMP&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$OMP&    vert_coords_patch, U, V, B1_zc, B2_zc, B3_zc,                 &
!$OMP&    gamma_eos_gr, gamma_soundspeed_gr, entropygr)                 &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

!$acc parallel vector_length(64) present(                               &
!$acc&    rhofloor, prsfloor, tempfloor, smallnum, ob_4pi, ob_8pi,      &
!$acc&    gasconst, molewt, gamma, magB_small_cgl, eps_bdry_cgl,        &
!$acc&    isospeed, gamma_euler, barotropic_stiffness,                  &
!$acc&    imask_active_patch, v_centroid_patch, v_scale_length,         &
!$acc&    vert_coords_patch, U, V, B1_zc, B2_zc, B3_zc,                 &
!$acc&    gamma_eos_gr, gamma_soundspeed_gr, entropygr)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, iyym1, izzm1, B0,             &
!$acc&    bflag, time_use,                                              &
!$acc&    a11, a12, a13, a21, x, y, z, obl0, iprs_flag, ent, Upt, Vpt)

      DO izz = 1 - NGhost, nr_zones + NGhost

        DO iyy = 1 - NGhost, ntp_zones + NGhost

          DO ixx = 1 - NGhost, ntp_zones + NGhost

          izzm1 = izz - 1
          iyym1 = iyy - 1
          ixxm1 = ixx - 1

! ------------------------------

          obl0 = 1.0
     1  / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

! ------------------------------

        B0 ( :) = 0.0

#if( USE_B0B1_FORMULATION == 1)
        x = v_centroid_patch ( ixx, iyy, izz, 1)
        y = v_centroid_patch ( ixx, iyy, izz, 2)
        z = v_centroid_patch ( ixx, iyy, izz, 3)

        bflag = 1         ! 1 == magnetic field;  0 == vector potential.
        time_use = 0.0    ! useful for the case of active rotation.

        CALL GET_MAGNETIC_FIELD ( bflag, time_use, x, y, z,
     1                              a11, a12, a13)

        B0 ( 1) = a11
        B0 ( 2) = a12
        B0 ( 3) = a13
#endif

! ------------------------------
! The call to "CONS_TO_PRIM" will also bring all variables within their physical
! ranges. As a result, after this call we feel free to reset "U, V, ent".

! If "PRESSURE_POSITIVITY > 0", "CONS_TO_PRIM" can even use "ent" to reset "U, V"
! in extreme situations.

          Upt ( :) = U ( ixx, iyy, izz, :, 1)

          ent = entropygr ( ixx, iyy, izz, 1)

          CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, Upt, Vpt)

          V ( ixx, iyy, izz, :, 1) = Vpt ( :)
          U ( ixx, iyy, izz, :, 1) = Upt ( :)
          entropygr ( ixx, iyy, izz, 1) = ent

! ------------------------------

          END DO  ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost" loop.

        END DO  ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost" loop.

      END DO  ! End of "DO izz = 1 - NGhost, nr_zones + NGhost" loop.

!$acc end parallel
!$OMP END PARALLEL DO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! The number of "end data" statements must match the "data copyin" statements.

!$acc end data

! ----------------------------------------------------------------------

      END SUBROUTINE CLEANUP_AFTER_UPDATE

