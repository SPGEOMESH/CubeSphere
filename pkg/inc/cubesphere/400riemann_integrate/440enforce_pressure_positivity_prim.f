




      SUBROUTINE ENFORCE_PRESSURE_POSITIVITY_PRIM

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_components, n_cc_modes,
     1         n_fc_modes,

     1         vert_coords_patch, imask_active_patch,

     1         xf_centroid_patch, yf_centroid_patch, zf_centroid_patch,
     1         v_scale_length, v_centroid_patch,

     1         istencil_patch,

     1         V, entropygr, B1_zc, B2_zc, B3_zc, flattengr,
     1         Coef_Basis_Moments)

! ----------------------------------------------------------------------
! Enforce positivity of the reconstructed solution.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes,
     1        n_fc_modes

! ------------------------------
! Variables that generate the mesh

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: ntp_zones + NGhost,
     1                             1 - NGhost: nr_zones + NGhost)

! ------------------------------
! Facial Variables

      REAL xf_centroid_patch ( 0 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

      REAL yf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

      REAL zf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost, 3)

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------

      INTEGER istencil_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                         1 - NGhost: ntp_zones + NGhost, 2, 12)

! ------------------------------
! Fluid and Magnetic Variables:-

      REAL V ( 1 - NGhost: ntp_zones + NGhost,
     1         1 - NGhost: ntp_zones + NGhost,
     1         1 - NGhost:  nr_zones + NGhost,
     1                      n_cc_components, n_cc_modes)

      REAL entropygr ( 1 - NGhost: ntp_zones + NGhost,
     1                 1 - NGhost: ntp_zones + NGhost,
     1                 1 - NGhost:  nr_zones + NGhost, 4)

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost,
     1      n_cc_modes) :: B1_zc, B2_zc, B3_zc

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost) :: flattengr

! ------------------------------
! MOMENTS FOR BASIS FUNCTIONS (set in INIT_BASIS_FUNCTION_COEFS):-

      REAL Coef_Basis_Moments ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes)

! ----------------------------------------------------------------------
! Local variables.

      INTEGER ixx, iyy, izz, ixxm1, iyym1, izzm1,
     1        ixxp1, iyyp1, izzp1, pow, icomp, irow, irr,

     1        ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd,
     1        ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,

     1        i, j, nn, ii, n_iter_pres_pos


      INTEGER iter, ipass, ipass_cgl, inode, bflag, iprs_flag

      REAL, DIMENSION ( n_cc_components) :: V_min, V_max, V_diff, Upt

      REAL, DIMENSION ( n_cc_modes) :: Coef_0

      REAL v_cen ( 3), ent_min, ent_max, ent_diff, ff, omff, expandfac,
     1     ratio_pres_pos,
     1     aa, bb, cc, gamma_eos_pt,
     1     vx, vy, vz, x, y, z, obl0,

     1     V_pt ( n_cc_components, n_cc_modes),
     1     V_pt_sav ( n_cc_components, n_cc_modes),

     1     B0_pt ( 3, n_cc_modes), B0 ( 3),

     1     ent_pt ( 4),

     1     vert ( 3, 27), xyz_node ( 3, 27),
     1     V_node ( n_cc_components, 27),
     1     ent_node ( 27),

     1     B0_node ( 3, 27), B0_x, B0_y, B0_z, time_use,
     1     p_par, p_perp, Bsq


! ----------------------------------------------------------------------
! Ratio by which the variation within a zone is reduced because of pressure
! positivity enforcement. If positivity cannot be achieved at 20% of the original
! variation, we just give up and fully flatten the zone. Increasing "n_iter_pres_pos"
! means that we reduce the variation much more slowly in a larger number of steps.

      PARAMETER (
     1      n_iter_pres_pos = 15,
     1      ratio_pres_pos  = ( 0.10)**( 1.0 / REAL(n_iter_pres_pos)),

! Set the factor by which we want to expand the range of each variable.
! Setting this to 0 gives us the TVD range in multi-dimensions.
     1      expandfac = 0.2)

!$acc routine( PASSED_CGL_PHYSICAL ) seq
!$acc routine( GET_MAGNETIC_FIELD ) seq
! ----------------------------------------------------------------------

! When "PRESSURE_POSITIVITY == 0" we don't enforce positivity at all.
! When "PRESSURE_POSITIVITY == 1" we only enforce positivity on density and pressure.
! When "PRESSURE_POSITIVITY == 2" we require all primitives to stay within the
!                                 range "V_min ( :)" to "V_max ( :)".

! ----------------------------------------------------------------------
#if ( PRESSURE_POSITIVITY >= 1)
! ----------------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, ixxm1, ixxp1, iyym1, iyyp1,  &
!$OMP&    izzm1, izzp1, ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd, &
!$OMP&    ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,       &

!$OMP&    obl0, v_cen, V_min, V_max, ent_min, ent_max, V_diff, ent_diff,&
!$OMP&    ff, omff, V_pt, V_pt_sav, ent_pt, vert, xyz_node,             &
!$OMP&    iter, inode, x, y, z, ent_node, V_node, iprs_flag,            &
!$OMP&    p_par, p_perp, ipass_cgl, Bsq,                                &
!$OMP&    aa, bb, cc, ipass, ii, bflag, time_use, B0, Upt, Coef_0)      &

!$OMP&            SHARED ( nr_zones, ntp_zones, NGhost, n_cc_components,&
!$OMP&    n_cc_modes,                                                   &
!$OMP&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$OMP&    molewt, gasconst, rhofloor, prsfloor, tempfloor, smallnum,    &
!$OMP&    istencil_patch, v_scale_length, v_centroid_patch,             &
!$OMP&    V, entropygr, flattengr, B1_zc, B2_zc, B3_zc,                 &
!$OMP&    vert_coords_patch, xf_centroid_patch,                         &
!$OMP&    yf_centroid_patch, zf_centroid_patch, Coef_Basis_Moments)     &

!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 3)

!$acc parallel vector_length(64) present (                              &
!$acc&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$acc&    molewt, gasconst, rhofloor, prsfloor, tempfloor, smallnum,    &
!$acc&    istencil_patch, v_scale_length, v_centroid_patch,             &
!$acc&    V, entropygr, flattengr, B1_zc, B2_zc, B3_zc,                 &
!$acc&    vert_coords_patch, xf_centroid_patch,                         &
!$acc&    yf_centroid_patch, zf_centroid_patch, Coef_Basis_Moments)

!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz, ixxm1, ixxp1, iyym1, iyyp1,          &
!$acc&    izzm1, izzp1, ixxr, iyyr, ixxl, iyyl, ixxu, iyyu, ixxd, iyyd, &
!$acc&    ixxru, iyyru, ixxlu, iyylu, ixxld, iyyld, ixxrd, iyyrd,       &

!$acc&    obl0, v_cen, V_min, V_max, ent_min, ent_max, V_diff, ent_diff,&
!$acc&    ff, omff, V_pt, V_pt_sav, ent_pt, vert, xyz_node,             &
!$acc&    iter, inode, x, y, z, ent_node, V_node, iprs_flag,            &
!$acc&    p_par, p_perp, ipass_cgl, Bsq,                                &
!$acc&    aa, bb, cc, ipass, ii, bflag, time_use, B0, Upt, Coef_0)

! ----------------------------------------------------------------------

      DO izz = - 2, nr_zones + 3

        DO iyy = - 2, ntp_zones + 3

          DO ixx = - 2, ntp_zones + 3

! ------------------------------------------------------------

      izzm1 = izz - 1
      izzp1 = izz + 1

      iyym1 = iyy - 1
      iyyp1 = iyy + 1

      ixxm1 = ixx - 1
      ixxp1 = ixx + 1

      ixxr = istencil_patch ( ixx, iyy, 1, 1)  ! For indexing zone one Right.
      iyyr = istencil_patch ( ixx, iyy, 2, 1)

      ixxl = istencil_patch ( ixx, iyy, 1, 2)  ! For indexing zone one Left.
      iyyl = istencil_patch ( ixx, iyy, 2, 2)

      ixxu = istencil_patch ( ixx, iyy, 1, 3)  ! For indexing zone one Up.
      iyyu = istencil_patch ( ixx, iyy, 2, 3)

      ixxd = istencil_patch ( ixx, iyy, 1, 4)  ! For indexing zone one Down.
      iyyd = istencil_patch ( ixx, iyy, 2, 4)

      ixxru = istencil_patch ( ixx, iyy, 1, 9)  ! For indexing zone that is Right-Up.
      iyyru = istencil_patch ( ixx, iyy, 2, 9)

      ixxlu = istencil_patch ( ixx, iyy, 1, 10)  ! For indexing zone that is Left-Up.
      iyylu = istencil_patch ( ixx, iyy, 2, 10)

      ixxld = istencil_patch ( ixx, iyy, 1, 11)  ! For indexing zone that is Left-Down.
      iyyld = istencil_patch ( ixx, iyy, 2, 11)

      ixxrd = istencil_patch ( ixx, iyy, 1, 12)  ! For indexing zone that is Right-Down.
      iyyrd = istencil_patch ( ixx, iyy, 2, 12)

! ----------------------------------------

      obl0 = 1.0
     1    / AMAX1 ( v_scale_length ( ixx, iyy, izz), smallnum)  ! Scaling for this zone.

      v_cen ( :) = v_centroid_patch ( ixx, iyy, izz, :)

      Coef_0 ( :) = Coef_Basis_Moments ( ixx, iyy, izz, :)

! ----------------------------------------

        B0 ( :) = 0.0

#if( USE_B0B1_FORMULATION == 1)
        x = v_cen ( 1)
        y = v_cen ( 2)
        z = v_cen ( 3)

        bflag = 1         ! 1 == magnetic field;  0 == vector potential.
        time_use = 0.0    ! useful for the case of active rotation.

        CALL GET_MAGNETIC_FIELD ( bflag, time_use, x, y, z, aa, bb, cc)

        B0 ( 1) = aa
        B0 ( 2) = bb
        B0 ( 3) = cc
#endif

! ----------------------------------------
! Find the min and max of "V"  and "entropy" active zones around (or)
! active von Neumann neighbor. First we pick up all the zones that we can
! even pick up from a tri-corner. This effectively gives us all 26 zones
! around a given zone.

      V_min ( :) = AMIN1 (
     1  V ( ixx,  iyy,  izz,   :, 1),
     1  V ( ixx,  iyy,  izzp1, :, 1),
     1  V ( ixx,  iyy,  izzm1, :, 1),

     1  V ( ixxr,  iyyr,  izz,   :, 1),  V ( ixxl,  iyyl,  izz,   :, 1),
     1  V ( ixxr,  iyyr,  izzp1, :, 1),  V ( ixxl,  iyyl,  izzp1, :, 1),
     1  V ( ixxr,  iyyr,  izzm1, :, 1),  V ( ixxl,  iyyl,  izzm1, :, 1),

     1  V ( ixxu,  iyyu,  izz,   :, 1),  V ( ixxd,  iyyd,  izz,   :, 1),
     1  V ( ixxu,  iyyu,  izzp1, :, 1),  V ( ixxd,  iyyd,  izzp1, :, 1),
     1  V ( ixxu,  iyyu,  izzm1, :, 1),  V ( ixxd,  iyyd,  izzm1, :, 1),

     1  V ( ixxru, iyyru, izz,   :, 1),  V ( ixxlu, iyylu, izz,   :, 1),
     1  V ( ixxru, iyyru, izzp1, :, 1),  V ( ixxlu, iyylu, izzp1, :, 1),
     1  V ( ixxru, iyyru, izzm1, :, 1),  V ( ixxlu, iyylu, izzm1, :, 1),

     1  V ( ixxld, iyyld, izz,   :, 1),  V ( ixxrd, iyyrd, izz,   :, 1),
     1  V ( ixxld, iyyld, izzp1, :, 1),  V ( ixxrd, iyyrd, izzp1, :, 1),
     1  V ( ixxld, iyyld, izzm1, :, 1),  V ( ixxrd, iyyrd, izzm1, :, 1))


      V_max ( :) = AMAX1 (
     1  V ( ixx,  iyy,  izz,   :, 1),
     1  V ( ixx,  iyy,  izzp1, :, 1),
     1  V ( ixx,  iyy,  izzm1, :, 1),

     1  V ( ixxr,  iyyr,  izz,   :, 1),  V ( ixxl,  iyyl,  izz,   :, 1),
     1  V ( ixxr,  iyyr,  izzp1, :, 1),  V ( ixxl,  iyyl,  izzp1, :, 1),
     1  V ( ixxr,  iyyr,  izzm1, :, 1),  V ( ixxl,  iyyl,  izzm1, :, 1),

     1  V ( ixxu,  iyyu,  izz,   :, 1),  V ( ixxd,  iyyd,  izz,   :, 1),
     1  V ( ixxu,  iyyu,  izzp1, :, 1),  V ( ixxd,  iyyd,  izzp1, :, 1),
     1  V ( ixxu,  iyyu,  izzm1, :, 1),  V ( ixxd,  iyyd,  izzm1, :, 1),

     1  V ( ixxru, iyyru, izz,   :, 1),  V ( ixxlu, iyylu, izz,   :, 1),
     1  V ( ixxru, iyyru, izzp1, :, 1),  V ( ixxlu, iyylu, izzp1, :, 1),
     1  V ( ixxru, iyyru, izzm1, :, 1),  V ( ixxlu, iyylu, izzm1, :, 1),

     1  V ( ixxld, iyyld, izz,   :, 1),  V ( ixxrd, iyyrd, izz,   :, 1),
     1  V ( ixxld, iyyld, izzp1, :, 1),  V ( ixxrd, iyyrd, izzp1, :, 1),
     1  V ( ixxld, iyyld, izzm1, :, 1),  V ( ixxrd, iyyrd, izzm1, :, 1))


      ent_min = AMIN1 (
     1   entropygr( ixx,  iyy,  izz,   1),
     1   entropygr( ixx,  iyy,  izzp1, 1),
     1   entropygr( ixx,  iyy,  izzm1, 1),

     1   entropygr( ixxr,  iyyr,  izz,   1),
     1   entropygr( ixxr,  iyyr,  izzp1, 1),
     1   entropygr( ixxr,  iyyr,  izzm1, 1),

     1   entropygr( ixxl,  iyyl,  izz,   1),
     1   entropygr( ixxl,  iyyl,  izzp1, 1),
     1   entropygr( ixxl,  iyyl,  izzm1, 1),

     1   entropygr( ixxu,  iyyu,  izz,   1),
     1   entropygr( ixxu,  iyyu,  izzp1, 1),
     1   entropygr( ixxu,  iyyu,  izzm1, 1),

     1   entropygr( ixxd,  iyyd,  izz,   1),
     1   entropygr( ixxd,  iyyd,  izzp1, 1),
     1   entropygr( ixxd,  iyyd,  izzm1, 1),

     1   entropygr( ixxru, iyyru, izz,   1),
     1   entropygr( ixxru, iyyru, izzp1, 1),
     1   entropygr( ixxru, iyyru, izzm1, 1),

     1   entropygr( ixxlu, iyylu, izz,   1),
     1   entropygr( ixxlu, iyylu, izzp1, 1),
     1   entropygr( ixxlu, iyylu, izzm1, 1),

     1   entropygr( ixxld, iyyld, izz,   1),
     1   entropygr( ixxld, iyyld, izzp1, 1),
     1   entropygr( ixxld, iyyld, izzm1, 1),

     1   entropygr( ixxrd, iyyrd, izz,   1),
     1   entropygr( ixxrd, iyyrd, izzp1, 1),
     1   entropygr( ixxrd, iyyrd, izzm1, 1))


      ent_max = AMAX1 (
     1   entropygr( ixx,  iyy,  izz,   1),
     1   entropygr( ixx,  iyy,  izzp1, 1),
     1   entropygr( ixx,  iyy,  izzm1, 1),

     1   entropygr( ixxr,  iyyr,  izz,   1),
     1   entropygr( ixxr,  iyyr,  izzp1, 1),
     1   entropygr( ixxr,  iyyr,  izzm1, 1),

     1   entropygr( ixxl,  iyyl,  izz,   1),
     1   entropygr( ixxl,  iyyl,  izzp1, 1),
     1   entropygr( ixxl,  iyyl,  izzm1, 1),

     1   entropygr( ixxu,  iyyu,  izz,   1),
     1   entropygr( ixxu,  iyyu,  izzp1, 1),
     1   entropygr( ixxu,  iyyu,  izzm1, 1),

     1   entropygr( ixxd,  iyyd,  izz,   1),
     1   entropygr( ixxd,  iyyd,  izzp1, 1),
     1   entropygr( ixxd,  iyyd,  izzm1, 1),

     1   entropygr( ixxru, iyyru, izz,   1),
     1   entropygr( ixxru, iyyru, izzp1, 1),
     1   entropygr( ixxru, iyyru, izzm1, 1),

     1   entropygr( ixxlu, iyylu, izz,   1),
     1   entropygr( ixxlu, iyylu, izzp1, 1),
     1   entropygr( ixxlu, iyylu, izzm1, 1),

     1   entropygr( ixxld, iyyld, izz,   1),
     1   entropygr( ixxld, iyyld, izzp1, 1),
     1   entropygr( ixxld, iyyld, izzm1, 1),

     1   entropygr( ixxrd, iyyrd, izz,   1),
     1   entropygr( ixxrd, iyyrd, izzp1, 1),
     1   entropygr( ixxrd, iyyrd, izzm1, 1))


! -------------------------------
! Now extend this range if "flattengr" allows it.

      V_diff ( :) = V_max ( :) - V_min ( :)

      ent_diff = ent_max - ent_min

      ff = AMAX1 ( AMIN1 ( (
     1           flattengr ( ixx, iyy, izz) - 0.1) / 0.4, 1.0), 0.0)

      omff = 1.0 - ff

! -------------------------------
! expansion of the range of variation that is allowed to the pressure positivity.
! For variables that are guaranteed to be positive, their "V_min, V_max" can be set
! by a multiplicative expansion. For variables that are not guaranteed positive, the
! min max should be set with "V_diff".

      V_min ( 1) = AMAX1 ( V_min ( 1) * ( 1.0 - expandfac * omff),
     1                     rhofloor)
      V_max ( 1) = AMAX1 ( V_max ( 1) * ( 1.0 + expandfac * omff),
     1                     rhofloor)

      V_min ( 5) = AMAX1 ( V_min ( 5) * ( 1.0 - expandfac * omff),
     1                     prsfloor)
      V_max ( 5) = AMAX1 ( V_max ( 5) * ( 1.0 + expandfac * omff),
     1                     prsfloor)

#if( IPDE_SYSTEM == 2)

      V_min ( 9) = AMAX1 ( V_min ( 9) * ( 1.0 - expandfac * omff),
     1                     prsfloor)
      V_max ( 9) = AMAX1 ( V_max ( 9) * ( 1.0 + expandfac * omff),
     1                     prsfloor)

#endif

#if( IPDE_SYSTEM == 3)

      V_min ( 9) = AMAX1 ( V_min ( 9) * ( 1.0 - expandfac * omff),
     1                     prsfloor)
      V_max ( 9) = AMAX1 ( V_max ( 9) * ( 1.0 + expandfac * omff),
     1                     prsfloor)

      V_min ( 10) = AMAX1 ( V_min ( 10) * ( 1.0 - expandfac * omff),
     1                     prsfloor)
      V_max ( 10) = AMAX1 ( V_max ( 10) * ( 1.0 + expandfac * omff),
     1                     prsfloor)

#endif

! ----------

      V_min ( 2:4) = V_min ( 2:4) - 0.5*expandfac * omff * V_diff ( 2:4)
      V_max ( 2:4) = V_max ( 2:4) + 0.5*expandfac * omff * V_diff ( 2:4)

      V_min ( 6:8) = V_min ( 6:8) - 0.5*expandfac * omff * V_diff ( 6:8)
      V_max ( 6:8) = V_max ( 6:8) + 0.5*expandfac * omff * V_diff ( 6:8)

! ----------
! For species fractions in MHD, or tracers in CGL with electron pressure, we have
! to take the min and max of the "mass fractions".

#if ( IPDE_SYSTEM == 1 && NFLUID >= 1)
      V_min ( 8 + 1: NVAR + NFLUID) = V_min ( 8 + 1: NVAR + NFLUID)
     1                              * (1.0 - expandfac*omff)
      V_max ( 8 + 1: NVAR + NFLUID) = V_max ( 8 + 1: NVAR + NFLUID)
     1                              * (1.0 + expandfac*omff)
#endif

#if ( IPDE_SYSTEM == 2 && NVAR >= 10)
      V_min ( 9 + 1: ) = V_min ( 9 + 1: ) * (1.0 - expandfac*omff)
      V_max ( 9 + 1: ) = V_max ( 9 + 1: ) * (1.0 + expandfac*omff)
#endif

#if ( IPDE_SYSTEM == 3 && NVAR >= 11)
      V_min ( 10 + 1: ) = V_min ( 10 + 1: ) * (1.0 - expandfac*omff)
      V_max ( 10 + 1: ) = V_max ( 10 + 1: ) * (1.0 + expandfac*omff)
#endif

! ----------

      ent_min = ent_min * ( 1.0 - expandfac * omff)
      ent_max = ent_max * ( 1.0 + expandfac * omff)

! --------------------
! Conserved variable and its higher moments; extracted just for this zone.
        V_pt ( :, :) = V ( ixx, iyy, izz, :, :)

        V_pt_sav ( :, :) = V_pt ( :, :)  ! Save it for use with magnetic field limiting.

! Entropy variable and its higher moments; extracted just for this zone.
        ent_pt ( :) = entropygr ( ixx, iyy, izz, :)

! ------------------------------------------------------------
! Pressure Positivity Enforcement Begins Here.

! ------------------------------
! 8-vertices

      vert ( :, 1) = vert_coords_patch ( ixxm1, iyym1, izzm1, :)
      vert ( :, 2) = vert_coords_patch ( ixx,   iyym1, izzm1, :)
      vert ( :, 3) = vert_coords_patch ( ixxm1, iyy,   izzm1, :)
      vert ( :, 4) = vert_coords_patch ( ixx,   iyy,   izzm1, :)
      vert ( :, 5) = vert_coords_patch ( ixxm1, iyym1, izz,   :)
      vert ( :, 6) = vert_coords_patch ( ixx,   iyym1, izz,   :)
      vert ( :, 7) = vert_coords_patch ( ixxm1, iyy,   izz,   :)
      vert ( :, 8) = vert_coords_patch ( ixx,   iyy,   izz,   :)

! ------------------------------
! 6-face centers

      vert ( :,  9) = xf_centroid_patch ( ixx,   iyy,   izz,   :)
      vert ( :, 10) = xf_centroid_patch ( ixxm1, iyy,   izz,   :)

      vert ( :, 11) = yf_centroid_patch ( ixx,   iyy,   izz,   :)
      vert ( :, 12) = yf_centroid_patch ( ixx,   iyym1, izz,   :)

      vert ( :, 13) = zf_centroid_patch ( ixx,   iyy,   izz,   :)
      vert ( :, 14) = zf_centroid_patch ( ixx,   iyy,   izzm1, :)

! ------------------------------
! 12 - edge centers; curvature of z-face edges not taken into account

      vert ( :, 15) = ( vert ( :, 1) + vert ( :, 2)) * 0.5
      vert ( :, 16) = ( vert ( :, 3) + vert ( :, 4)) * 0.5
      vert ( :, 17) = ( vert ( :, 5) + vert ( :, 6)) * 0.5
      vert ( :, 18) = ( vert ( :, 7) + vert ( :, 8)) * 0.5

      vert ( :, 19) = ( vert ( :, 1) + vert ( :, 3)) * 0.5
      vert ( :, 20) = ( vert ( :, 2) + vert ( :, 4)) * 0.5
      vert ( :, 21) = ( vert ( :, 5) + vert ( :, 7)) * 0.5
      vert ( :, 22) = ( vert ( :, 6) + vert ( :, 8)) * 0.5

      vert ( :, 23) = ( vert ( :, 1) + vert ( :, 5)) * 0.5
      vert ( :, 24) = ( vert ( :, 2) + vert ( :, 6)) * 0.5
      vert ( :, 25) = ( vert ( :, 3) + vert ( :, 7)) * 0.5
      vert ( :, 26) = ( vert ( :, 4) + vert ( :, 8)) * 0.5

! ------------------------------
! 1 - zone center.

      vert ( :, 27) = v_cen ( :)

! ------------------------------
! Compute the node locations on the local coordinates.

!$acc loop seq
      DO i = 1, 27
        xyz_node ( :, i) = ( vert ( :, i) - v_cen ( :)) * obl0
      END DO

! ------------------------------------------------------------
! ------------------------------------------------------------
! STEP 8.a) Start the iteration loop. In a very low cost way, this iteration
! makes sure that the primitive values within this zone are within the range
! of [V_min ( :), V_max ( :)].

!$acc loop seq
      DO iter = 1, n_iter_pres_pos + 1

! ------------------------------------------------------------
! On Fourth iteration, if we are still in this loop, we zero out the modes
! and exit. Magnetic fields will have a very restricted linear variation
! consistent with first order, div-free reconstruction of vector fields.

! If all the variation has been brought to its lowest value, then in STEP 5
! of "LIMIT_ALL_SLOPES" we have already ensured that this will have positive
! pressure.

      IF ( iter == n_iter_pres_pos + 1) THEN

        V_pt ( 1: n_cc_components, 2: n_cc_modes) = 0.0

        ent_pt ( :) = 0.0

        EXIT

      END IF

! ------------------------------------------------------------

        DO inode = 1, 27

! --------------------------------------------------
! Build the "U, entropy" at the xyz_nodes.

        x = xyz_node ( 1, inode)
        y = xyz_node ( 2, inode)
        z = xyz_node ( 3, inode)

        V_node ( :, inode) = V_pt ( :, 1)
     1                     + V_pt ( :, 2) * x
     1                     + V_pt ( :, 3) * y
     1                     + V_pt ( :, 4) * z
#if ( NUMORDER >= 3)
     1                     + V_pt ( :, 5) * ( x**2 - Coef_0 ( 5))
     1                     + V_pt ( :, 6) * ( y**2 - Coef_0 ( 6))
     1                     + V_pt ( :, 7) * ( z**2 - Coef_0 ( 7))
     1                     + V_pt ( :, 8) * ( x * y - Coef_0 ( 8))
     1                     + V_pt ( :, 9) * ( y * z - Coef_0 ( 9))
     1                     + V_pt ( :, 10) * ( x * z - Coef_0 ( 10))
#elif ( NUMORDER >= 4)
     1                     + V_pt ( :, 11) * ( x**3 - Coef_0 ( 11))
     1                     + V_pt ( :, 12) * ( y**3 - Coef_0 ( 12))
     1                     + V_pt ( :, 13) * ( z**3 - Coef_0 ( 13))
     1                     + V_pt ( :, 14) * ( x**2 * y - Coef_0 ( 14))
     1                     + V_pt ( :, 15) * ( x**2 * z - Coef_0 ( 15))
     1                     + V_pt ( :, 16) * ( x * y**2 - Coef_0 ( 16))
     1                     + V_pt ( :, 17) * ( y**2 * z - Coef_0 ( 17))
     1                     + V_pt ( :, 18) * ( x * z**2 - Coef_0 ( 18))
     1                     + V_pt ( :, 19) * ( y * z**2 - Coef_0 ( 19))
     1                     + V_pt ( :, 20) * ( x * y * z - Coef_0 ( 20))
#endif

        ent_node ( inode) = ent_pt ( 1)
     1                    + ent_pt ( 2) * x
     1                    + ent_pt ( 3) * y
     1                    + ent_pt ( 4) * z

! --------------------------------------------------

        END DO ! End of "DO inode = 1, 27".

! ------------------------------------------------------------

        ipass = 1

        DO inode = 1, 27

! --------------------------------------------------
! Check only for positivity of density or pressure. These are components 1 and 5.
! This is a less stringent test; based only on the scalar density and pressure.

         ii = 1
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0

         ii = 5
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0


! We also include the entropy variable.

        IF ( ( ent_node ( inode) < ent_min ) .OR.
     1       ( ent_node ( inode) > ent_max ) ) ipass = 0

! --------------------------------------------------

#if( IPDE_SYSTEM == 1)

! ------------------------------

        DO ii = NVAR + 1, NVAR + NFLUID

          IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1         ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0

        END DO  ! End of "DO ii = NVAR + 1, NVAR + NFLUID".

! ------------------------------

#endif
! End of "IPDE_SYSTEM == 1".

! --------------------------------------------------

#if( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3)

! ------------------------------
! Make sure that parallel and perpendicular pressures in the ions are within bounds.
! this is done for both types of CGL.

         ii = 9
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0


! Check whether "p_par" and "p_perp" are bounded by the mirror and firehose limits.

          p_par  = V_node ( 5, inode)
          p_perp = V_node ( 9, inode)

          Bsq = ( V_node ( 6, inode) + B0 ( 1))**2
     1        + ( V_node ( 7, inode) + B0 ( 2))**2
     1        + ( V_node ( 8, inode) + B0 ( 3))**2

         CALL PASSED_CGL_PHYSICAL ( p_par, p_perp, Bsq, ipass_cgl)

         IF ( ipass_cgl == 0) ipass = 0

! ------------------------------

#endif
! End of "IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3".

! --------------------------------------------------
#if ( IPDE_SYSTEM == 2)
! ------------------------------
! Make sure that species fractions are within bounds.

#if ( NVAR >= 10)
         ii = 10
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0
#endif

#if ( NVAR >= 11)
         ii = 11
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0
#endif

! ------------------------------
#endif
! End of "IPDE_SYSTEM == 2".
! --------------------------------------------------
#if ( IPDE_SYSTEM == 3)
! ------------------------------
! Make sure that electron pressure is within bounds. Also make sure species
! fractions are within bounds.

         ii = 10
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0

#if ( NVAR >= 11)
         ii = 11
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0
#endif

#if ( NVAR >= 12)
         ii = 12
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0
#endif

#if ( NVAR >= 13)
         ii = 13
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0
#endif

#if ( NVAR >= 14)
         ii = 14
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0
#endif

#if ( NVAR >= 15)
         ii = 15
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0
#endif

#if ( NVAR >= 16)
         ii = 16
         IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1        ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0
#endif

! ------------------------------
#endif
! End of "IPDE_SYSTEM == 3".
! --------------------------------------------------
#if( PRESSURE_POSITIVITY == 2 )
! Check for all variables to lie within range.
! This is a more stringent test; based on all variables.

        DO ii = 1, NVAR + NFLUID

          IF ( ( V_node ( ii, inode) < V_min ( ii) ) .OR.
     1         ( V_node ( ii, inode) > V_max ( ii) ) ) ipass = 0

        END DO  ! End of "DO ii = 1, NVAR + NFLUID".

#endif

! --------------------------------------------------

        END DO ! End of "DO inode = 1, 27".

! ------------------------------------------------------------
! All primitive variables are within their acceptable ranges, so we can
! move to next zone by exiting the "DO iter = 1, n_iter_pres_pos + 1" loop.

! For a zone that it a good zone ( i.e., with its primitive values lying
! within [V_min ( :), V_max ( :)] ), this pressure positivity step should
! cost almost nothing. In a good simulation, this adds almost no extra cost.

      IF ( ipass == 1) EXIT

! --------------------------------------------------
! At least one of the primitive variables is not within the acceptable range,
! so we reduce the variation in the modes.

! For fluid variables, this just entails reducing all the higher order modes.

      IF ( ipass == 0) THEN

      V_pt ( 1:NVAR + NFLUID, 2: n_cc_modes) =
     1     ratio_pres_pos * V_pt ( 1:NVAR + NFLUID, 2: n_cc_modes)

      ent_pt ( 2: 4) = ratio_pres_pos * ent_pt ( 2: 4)

      END IF

! ------------------------------------------------------------

      END DO ! End of "DO iter = 1, n_iter_pres_pos + 1".

! ------------------------------------------------------------
! Assign back the modes.

        V ( ixx, iyy, izz, :, :) = V_pt ( :, :)

        entropygr ( ixx, iyy, izz, :) = ent_pt ( :)

        B1_zc ( ixx, iyy, izz, :) = V ( ixx, iyy, izz, 6, :)
        B2_zc ( ixx, iyy, izz, :) = V ( ixx, iyy, izz, 7, :)
        B3_zc ( ixx, iyy, izz, :) = V ( ixx, iyy, izz, 8, :)

! ------------------------------------------------------------

          END DO  ! End of "DO ixx = - 2, ntp_zones + 3" loop.

        END DO  ! End of "DO iyy = - 2, ntp_zones + 3" loop.

      END DO  ! End of "DO izz = - 2, nr_zones + 3" loop.

!$acc end parallel
!$OMP END PARALLEL DO
! ----------------------------------------------------------------------
#endif
! ----------------------------------------------------------------------

      END SUBROUTINE ENFORCE_PRESSURE_POSITIVITY_PRIM



