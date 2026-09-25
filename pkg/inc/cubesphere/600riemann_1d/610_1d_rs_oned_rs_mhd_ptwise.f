! ----------------------------------------------------------------------
! This is new Oned_RS_that may or may not work.
! ----------------------------------------------------------------------









      SUBROUTINE OneD_RS_MHD_Ptwise

     1         ( gamma, gamma_euler, isospeed, barotropic_stiffness,
     1           molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1           magB_small_cgl, eps_bdry_cgl,

     1           fourpi, ob_4pi, ob_8pi,
     1           smallnum, B0, flat,

     1           V_L, V_R, ent_L, ent_R,
     1           S_L, S_R, U_Star, F_Star, ent_F_Star,
     1           Fluct_L, Fluct_R,
     1           U_L, U_R, F_L, F_R)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! 1D Riemann solver applied by convention in the x-direction. Since variables
! are rotated into the local variables at the face, we only use x-direction.
! Simplest HLL RS is used. A single, "gamma" law gas is assumed.

! "V_L, V_R" Left and Right states in primitive variables, fed in as input.

! "ent_L, ent_R" left and right entropy densities, inputs only used if
!       "PRESSURE_POSITIVITY >= 1".

! "S_L, S_R" extremal speeds, returned as outputs.

! "U_Star, F_Star" Resolved state and flux; these are the primary outputs.
! "ent_F_Star" flux of entropy density; output only provided if
!       "PRESSURE_POSITIVITY >= 1".

! "flat" ranges from 0.0 (can put in intermediate wave) to 1.0 (don't put in intermediate wave).

! "Fluct_L, Fluct_R" left and right-going fluctuations returned as outputs.
! These are most important because CGL has non-conservative products, which are
! contained in these variables.

! "U_L, U_R" conserved variables for the left and right states; auxiliary outputs.
! "F_L, F_R" fluxes associated with the left and right states; auxiliary outputs.

! ----------------------------------------------------------------------
!$acc routine seq

      IMPLICIT NONE

      REAL gamma, gamma_euler, isospeed, barotropic_stiffness,
     1     molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1     magB_small_cgl, eps_bdry_cgl, flat,
     1     fourpi, ob_4pi, ob_8pi, smallnum, tempaa, tempab, wt

      REAL, DIMENSION ( 3) :: B0

      REAL :: S_L, S_R,
     1      ent_L, ent_R, ent_F_Star

      REAL, DIMENSION ( NVAR + NFLUID) ::

     1      V_L, V_R, U_Star, F_Star, U_L, U_R, F_L, F_R,
     1      Fluct_L, Fluct_R

! ----------------------------------------------------------------------
! Variables past this point are internal to the subroutine.

      INTEGER idirection, nc_iter, i_nc, iprs_flag, iprimitive, i

      REAL S_min, S_max,
     1     temp_01, temp_02, temp_03, temp_04, bfxre,
     1     nc_product ( 2) , ent_U_Star, prs_from_ent, prs_from_eng,
     1     phi

      REAL, DIMENSION ( NVAR + NFLUID) :: V, V_use, V_Star,

     1    V_zc, dV_zc_dx, dV_zc_dy, dV_zc_dz

      REAL, DIMENSION ( NVAR + NFLUID) ::
     1     eigenval, Eigenwt, Antidiff

      REAL, DIMENSION ( NVAR + NFLUID, NVAR + NFLUID) ::
     1     lf_eigenvec, rt_eigenvec, delta_x

! ----------------------------------------------------------------------
! STEP 1) Impose continuity of the longitudinal magnetic field. Then
! obtain conserved variables "U_L, U_R" from primitive variables "V_L, V_R"
! and also entropy variables "ent_L, ent_R" for positivity.

      bfxre = 0.5 * ( V_L ( 6) + V_R ( 6) )
      V_L ( 6) = bfxre
      V_R ( 6) = bfxre

      nc_iter = 3  ! Number of iterations for the non-conservative products.

      idirection = 1  ! By default, choose x-direction.

! ----------

      CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent_L, B0, U_L, V_L)

! ----------

      CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent_R, B0, U_R, V_R)

! ------------------------------
! STEP 2) Obtain "F_L, F_R" for the fluxes.

! ----------

      CALL CONS_PRIM_TO_FLUX_3D ( idirection, gamma, B0, V_L, F_L)

! ----------

      CALL CONS_PRIM_TO_FLUX_3D ( idirection, gamma, B0, V_R, F_R)

! ------------------------------
! STEP 3) Obtain extremal wave speeds "S_R, S_L".

      S_R = - 1.0e30
      S_L = 1.0e30

! ----------

      CALL PRIM_TO_WAVSPD ( idirection, gamma, B0, V_L, S_min, S_max)

      S_L = AMIN1 ( S_L, S_min)  ;  S_R = AMAX1 ( S_R, S_max)

! ----------

      CALL PRIM_TO_WAVSPD ( idirection, gamma, B0, V_R, S_min, S_max)

      S_L = AMIN1 ( S_L, S_min)  ;  S_R = AMAX1 ( S_R, S_max)

! ----------

      V ( :) = 0.5 * ( V_L ( :) + V_R ( :) )

      CALL PRIM_TO_WAVSPD ( idirection, gamma, B0, V, S_min, S_max)

      S_L = AMIN1 ( S_L, S_min)  ;  S_R = AMAX1 ( S_R, S_max)

! Build the eigenvectors but do this only if this is called for.

#if ( COMPLETE_EIGENSYSTEM >= 1 && IPDE_SYSTEM == 1)

      iprimitive = 0

      CALL EIGENVEC_PTWISE_LinDeg

     1             ( iprimitive, idirection,

     1               gamma, smallnum, prsfloor, rhofloor, fourpi,

     1               V,

     1               eigenval,
     1               lf_eigenvec, rt_eigenvec )

#endif

! ------------------------------
! STEP 4) Make an initial "U_Star" that integrates from "U_L" to "U_R".
! We will iterate on this in the next step if we have a non-conservative product.

! If we have "PRESSURE_POSITIVITY > 0" then we also build a guaranteed positive
! "prs_from_ent" to guide "U_Star" (and fluxes and fluctuations resulting from it)
! towards a pressure positive state.

! ------------------------

      temp_03 = 1.0 / ( S_R - S_L)
      temp_01 = S_R * temp_03
      temp_02 = S_L * temp_03
      temp_04 = S_R * S_L * temp_03

      U_Star = temp_01 * U_R - temp_02 * U_L - temp_03 * ( F_R - F_L)

! ------------------------
! We also build "ent_U_Star" with the hope that it has a guaranteed pressure.
! We can use that pressure to correct for MHD and CGL in case "U_Star ( 5)" gives
! a negative pressure.

      ent_U_Star = temp_01 * ent_R - temp_02 * ent_L
     1           - temp_03 * ( ent_R * V_R ( 2) - ent_L * V_L ( 2) )

#if ( IPDE_SYSTEM == 1 && PRESSURE_POSITIVITY > 0)

      prs_from_ent = ent_U_Star * U_Star ( 1)**(gamma - 1.0)

      prs_from_eng = ( gamma - 1.0) * ( U_Star ( 5)
     1 - 0.5 * ( U_Star ( 2)**2 + U_Star ( 3)**2 + U_Star ( 4)**2)
     1 / U_Star ( 1)
     1 - ob_8pi * ( U_Star ( 6)**2 + U_Star ( 7)**2 + U_Star ( 8)**2) )

      IF ( prs_from_eng .LE. prsfloor) THEN

      U_Star ( 5) = U_Star ( 5) + ( prs_from_ent - prs_from_eng)
     1                          / ( gamma - 1.0)

      END IF


#elif ( IPDE_SYSTEM == 2 && PRESSURE_POSITIVITY > 0)

      prs_from_ent = ( ent_U_Star * U_Star ( 1)**4)**(1.0/3.0)

#elif ( IPDE_SYSTEM == 3 && PRESSURE_POSITIVITY > 0)

      prs_from_ent = ( ent_U_Star * U_Star ( 1)**4)**(1.0/3.0)

#endif

! ------------------------

! This part does not change in the fixed point iteration; so save it in "V".
! This "V" should not be touched during the course of the iteration.
      V = U_Star

! ------------------------
#if ( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3 )
! ------------------------
! Add in some contribution from the non-conservative products, needed for CGL.
! We use a 4-point gauss-lobatto quadrature rule for the integration.
! Only gradients in x-direction are needed.

      dV_zc_dx ( :) = V_R ( :) - V_L ( :)
      dV_zc_dy = 0.0  ;  dV_zc_dz = 0.0

! ----------
! First quadrature point using "V_L" and "V_R".

      tempab = 0.0
      tempaa = 1.0 - tempab
      wt = 1.0 / 12.0

      V_use = tempaa * V_L + tempab * V_R

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ----------
! Second quadrature point using "V_L" and "V_R".

      tempab = 0.5 * ( 1.0 - 1.0 / SQRT ( 5.0))
      tempaa = 1.0 - tempab
      wt = 5.0 / 12.0

      V_use = tempaa * V_L + tempab * V_R

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ----------
! Third quadrature point using "V_L" and "V_R".

      tempab = 0.5 * ( 1.0 + 1.0 / SQRT ( 5.0))
      tempaa = 1.0 - tempab
      wt = 5.0 / 12.0

      V_use = tempaa * V_L + tempab * V_R

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ----------
! Fourth quadrature point using "V_L" and "V_R".

      tempab = 1.0
      tempaa = 1.0 - tempab
      wt = 1.0 / 12.0

      V_use = tempaa * V_L + tempab * V_R

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ------------------------
! If we are doing CGL, we have no guarantee that "U_Star" is in the physical domain.
! Therefore, we push it back into the physical domain if it has gone out of the
! physical domain.

#if ( PRESSURE_POSITIVITY > 0)

      prs_from_eng = (2.0/3.0) * ( U_Star ( 5)
     1 - 0.5 * ( U_Star ( 2)**2 + U_Star ( 3)**2 + U_Star ( 4)**2)
     1 / U_Star ( 1)
     1 - ob_8pi * ( U_Star ( 6)**2 + U_Star ( 7)**2 + U_Star ( 8)**2) )

      IF ( prs_from_eng .LE. prsfloor) THEN

      U_Star ( 5) = U_Star ( 5) + 1.5 * ( prs_from_ent - prs_from_eng)

      END IF

#endif

! ------------------------
#endif
! End of "IPDE_SYSTEM == 2".
! ------------------------

! ------------------------------
! STEP 5) Iterate on "U_Star" with the non-conservative product.
! These are fixed point iterations. Only needed for CGL.

! ------------------------
#if ( IPDE_SYSTEM == 2 || IPDE_SYSTEM == 3 )
! ------------------------

      DO i_nc = 1, nc_iter

! Obtain "V_Star" from "U_Star".

      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, 0.5 * (ent_L + ent_R), B0, U_Star, V_Star)

! ----------
! Build gradients from "V_L" to "V_Star".

      dV_zc_dx ( :) = V_Star ( :) - V_L ( :)
      dV_zc_dy = 0.0  ;  dV_zc_dz = 0.0

! ----------
! First quadrature point using "V_L" and "V_Star".

      tempab = 0.0
      tempaa = 1.0 - tempab
      wt = 1.0 / 12.0

      V_use = tempaa * V_L + tempab * V_Star

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = V ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ----------
! Second quadrature point using "V_L" and "V_Star".

      tempab = 0.5 * ( 1.0 - 1.0 / SQRT ( 5.0))
      tempaa = 1.0 - tempab
      wt = 5.0 / 12.0

      V_use = tempaa * V_L + tempab * V_Star

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ----------
! Third quadrature point using "V_L" and "V_Star".

      tempab = 0.5 * ( 1.0 + 1.0 / SQRT ( 5.0))
      tempaa = 1.0 - tempab
      wt = 5.0 / 12.0

      V_use = tempaa * V_L + tempab * V_Star

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ----------
! Fourth quadrature point using "V_L" and "V_Star".

      tempab = 1.0
      tempaa = 1.0 - tempab
      wt = 1.0 / 12.0

      V_use = tempaa * V_L + tempab * V_Star

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ----------
! Build gradients from "V_Star" to "V_R".

      dV_zc_dx ( :) = V_R ( :) - V_Star ( :)
      dV_zc_dy = 0.0  ;  dV_zc_dz = 0.0

! ----------
! First quadrature point using "V_Star" and "V_R".

      tempab = 0.0
      tempaa = 1.0 - tempab
      wt = 1.0 / 12.0

      V_use = tempaa * V_Star + tempab * V_R

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ----------
! Second quadrature point using "V_Star" and "V_R".

      tempab = 0.5 * ( 1.0 - 1.0 / SQRT ( 5.0))
      tempaa = 1.0 - tempab
      wt = 5.0 / 12.0

      V_use = tempaa * V_Star + tempab * V_R

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ----------
! Third quadrature point using "V_Star" and "V_R".

      tempab = 0.5 * ( 1.0 + 1.0 / SQRT ( 5.0))
      tempaa = 1.0 - tempab
      wt = 5.0 / 12.0

      V_use = tempaa * V_Star + tempab * V_R

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ----------
! Fourth quadrature point using "V_Star" and "V_R".

      tempab = 1.0
      tempaa = 1.0 - tempab
      wt = 1.0 / 12.0

      V_use = tempaa * V_Star + tempab * V_R

      V_use ( 6) = V_use ( 6) + B0 ( 1)
      V_use ( 7) = V_use ( 7) + B0 ( 2)
      V_use ( 8) = V_use ( 8) + B0 ( 3)

      CALL MAKE_NC_CHARMAT_Prod_Prim
     1          (rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,

     1           V_use, dV_zc_dx, dV_zc_dy, dV_zc_dz, nc_product)

      U_Star ( 9) = U_Star ( 9) - temp_03 * nc_product ( 1) * wt
#if ( IPDE_SYSTEM == 3)
      U_Star ( 10) = U_Star ( 10) - temp_03 * nc_product ( 2) * wt
#endif

! ------------------------
! If we are doing CGL, we have no guarantee that "U_Star" is in the physical domain.
! Therefore, we push it back into the physical domain if it has gone out of the
! physical domain.

#if ( PRESSURE_POSITIVITY > 0)

      prs_from_eng = (2.0/3.0) * ( U_Star ( 5)
     1 - 0.5 * ( U_Star ( 2)**2 + U_Star ( 3)**2 + U_Star ( 4)**2)
     1 / U_Star ( 1)
     1 - ob_8pi * ( U_Star ( 6)**2 + U_Star ( 7)**2 + U_Star ( 8)**2) )

      IF ( prs_from_eng .LE. prsfloor) THEN

      U_Star ( 5) = U_Star ( 5) + 1.5 * ( prs_from_ent - prs_from_eng)

      END IF

#endif

! ------------------------

      END DO  ! End of "DO i_nc = 1, nc_iter".

! ------------------------
#endif
! End of "IPDE_SYSTEM == 2".

! ------------------------------
! STEP 6) Form the HLL flux for the conserved variables and for the entropy.
! Also build the fluctuations here.
! The logic provided here already enables us to pick out the supersonic limits.

! The "U_Star" that is used for building fluctuations and fluxes will have
! positive pressure if the pressure positivity has been turned on.

! --------------------------------------------------

      IF ( S_L >= 0.0) THEN

! --------------------------------------------------

        U_Star ( :) = U_L ( :)

        F_Star ( :) = F_L ( :)
        ent_F_Star = ent_L * V_L ( 2)

        Fluct_R = S_L * ( U_Star - U_L) + S_R * ( U_R - U_Star)
        Fluct_L = 0.0

! --------------------------------------------------

      ELSE IF ( S_R <= 0.0) THEN

! --------------------------------------------------

        U_Star ( :) = U_R ( :)

        F_Star ( :) = F_R ( :)
        ent_F_Star = ent_R * V_R ( 2)

        Fluct_L = S_L * ( U_Star - U_L) + S_R * ( U_R - U_Star)
        Fluct_R = 0.0

! --------------------------------------------------

      ELSE

! --------------------------------------------------

! The converged values of "U_Star, Fluct_L, Fluct_R" from above have already
! been iterated to their converged values. This form of "F_Star" is guaranteed
! to have positive pressure if we had non-zero pressure positivity.

        temp_01 = 1.0 / ( S_R - S_L)
        temp_02 = S_R * S_L * temp_01

        F_Star ( :) = 0.5 * ( F_L ( :) + F_R ( :))
     1              + 0.5 * S_L * ( U_Star ( :) - U_L ( :))
     1              - 0.5 * S_R * ( U_R ( :) - U_Star ( :))

        Fluct_L = S_L * ( U_Star - U_L)

        Fluct_R = S_R * ( U_R - U_Star)

        ent_F_Star = ( S_R * ent_L * V_L ( 2)
     1                      - S_L * ent_R * V_R ( 2))
     1                    * temp_01
     1  + ( ent_R - ent_L ) * temp_02

! ------------------------------
#if ( COMPLETE_EIGENSYSTEM >= 1 && IPDE_SYSTEM == 1)
! If the eigensystem is available then we can lower the dissipation in the Riemann solver.
! It is most beneficial to do this for the linearly degenerate waves.

! Decide what fraction of HLLI contribution we want in the RS.
      phi = 1.0 - AMIN1 ( AMAX1 ( flat, 0.0), 1.0)

! Include contributions from "S_L, S_R" for the anti-diffusive term in here.
      phi = phi * S_R * S_L / ( S_R - S_L)

! ----------

! The below evaluation of "delta_x" is safe and shuts down any waves
! with speeds that might slightly exceed the extremal wavespeeds.

      delta_x = 0.0

      DO i = 1, NVAR

#if (IPDE_SYSTEM == 1)
        IF ( ( i .EQ. 1) .OR. ( i .EQ. 3) .OR.
     1       ( i .EQ. 5) .OR. ( i .EQ. 7)) CYCLE  ! Only use linearly degenerate waves.
#endif

      IF ( ( S_L .LE. Eigenval ( i)) .AND.
     1     ( Eigenval ( i) .LE. S_R) ) THEN

        delta_x ( i, i) = 1.0
     1  - AMIN1 ( Eigenval ( i), 0.0) / AMIN1 ( S_L, - smallnum)
     1  - AMAX1 ( Eigenval ( i), 0.0) / AMAX1 ( S_R, smallnum)

      ELSE

        delta_x ( i, i) = 0.0

      END IF

      END DO

! ----------

! Once all those eigenweights are evaluated, we can obtain "Antidiff".
! Effect of flattener and extremal wave speeds is also put in here.

      Eigenwt = MATMUL ( lf_eigenvec, U_R - U_L)

      Antidiff = MATMUL ( rt_eigenvec,
     1                   MATMUL ( delta_x, Eigenwt) )

      Antidiff = phi * Antidiff

      Fluct_L = Fluct_L - Antidiff
      Fluct_R = Fluct_R + Antidiff

      F_Star = F_Star - Antidiff

#endif
! ------------------------------

      END IF

      END SUBROUTINE OneD_RS_MHD_Ptwise





