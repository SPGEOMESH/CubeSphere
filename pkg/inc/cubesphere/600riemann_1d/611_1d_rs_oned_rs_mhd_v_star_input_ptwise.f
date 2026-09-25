! ----------------------------------------------------------------------
! This is new Oned_RS_that may or may not work.
! ----------------------------------------------------------------------



#if ( CHOOSE_APPL_HELIO == 1)
      SUBROUTINE OneD_RS_MHD_V_Star_input_Ptwise

     1         ( gamma, gamma_euler, isospeed, barotropic_stiffness,
     1           molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1           magB_small_cgl, eps_bdry_cgl,

     1           fourpi, ob_4pi, ob_8pi,
     1           smallnum, B0,

     1           V_L, V_R, ent_L, ent_R,
     1           S_L, S_R,
     1           V_Star, U_Star, F_Star, ent_F_Star,
     1           Fluct_L, Fluct_R,
     1           U_L, U_R, F_L, F_R)

! ----------------------------------------------------------------------

! "V_Star" is fed in as input.

! 1D Riemann solver applied by convention in the x-direction. Since variables
! are rotated into the local variables at the face, we only use x-direction.
! Simplest HLL RS is used. A single, "gamma" law gas is assumed.

! "V_L, V_R" Left and Right states in primitive variables, fed in as input.

! "ent_L, ent_R" left and right entropy densities, inputs only used if
!       "PRESSURE_POSITIVITY >= 1".

! "S_L, S_R" extremal speeds, returned as outputs.

! "F_Star" Resolved flux; this is the primary output.
! "ent_F_Star" flux of entropy density; output only provided if
!       "PRESSURE_POSITIVITY >= 1".

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
     1     magB_small_cgl, eps_bdry_cgl,
     1     fourpi, ob_4pi, ob_8pi, smallnum, tempaa, tempab, wt

      REAL, DIMENSION ( 3) :: B0

      REAL :: S_L, S_R,
     1      ent_L, ent_R, ent_F_Star

      REAL, DIMENSION ( NVAR + NFLUID) ::

     1      V_L, V_R, U_Star, F_Star, U_L, U_R, F_L, F_R,
     1      Fluct_L, Fluct_R, V_Star

! ----------------------------------------------------------------------
! Variables past this point are internal to the subroutine.

      INTEGER idirection, nc_iter, i_nc, iprs_flag

      REAL S_min, S_max,
     1     temp_01, temp_02, temp_03, temp_04, bfxre,
     1     nc_product , ent_U_Star, prs_from_ent, prs_from_eng,
     1     ent_temp

      REAL, DIMENSION ( NVAR + NFLUID) :: V, V_use,

     1    V_zc, dV_zc_dx, dV_zc_dy, dV_zc_dz

! ----------------------------------------------------------------------
! STEP 1) Impose continuity of the longitudinal magnetic field. Then
! obtain conserved variables "U_L, U_R" from primitive variables "V_L, V_R"
! and also entropy variables "ent_L, ent_R" for positivity.

      bfxre = 0.5 * ( V_L ( 6) + V_R ( 6) )
      V_L ( 6) = bfxre
      V_R ( 6) = bfxre

      idirection = 1  ! By default, choose x-direction.

! ------
        ! Obtain "U_Star" from "V_Star"
        CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent_temp, B0, U_Star, V_Star)

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

! ------------------------------
! We also build "ent_U_Star" with the hope that it has a guaranteed pressure.
! We can use that pressure to correct for MHD and CGL in case "U_Star ( 5)" gives
! a negative pressure.

#if ( IPDE_SYSTEM == 1 )

      ent_U_Star =  V_Star ( 5) / V_Star ( 1)**(gamma - 1.0)

#endif

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

#endif
! ------------------------------

      CALL CONS_PRIM_TO_FLUX_3D ( idirection, gamma, B0, V_Star, F_Star)

      Fluct_L = F_Star - F_L

      Fluct_R = F_R - F_Star

      ent_F_Star = V_Star ( 2) * ent_U_Star

! ----------------------------------------------------------------------

      END SUBROUTINE OneD_RS_MHD_V_Star_input_Ptwise
#endif
! End of "CHOOSE_APPL_HELIO == 1"

