! ----------------------------------------------------------------------
! This is Oned_RS_DimeADozen
! ----------------------------------------------------------------------





      SUBROUTINE OneD_RS_MHD_DAD ( NumRS, NumRS_Max,

     1           rhofloor, prsfloor, smallnum, gamma,
     1           fourpi, ob_4pi, ob_8pi,

     1           isospeed, barotropic_stiffness, gamma_euler,
     1           B0,

     1           U_L, U_R, ent_L, ent_R,
     1           S_L, S_R, U_Star, F_Star, ent_F_Star,
     1           V_L, V_R, F_L, F_R)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! Takes a dime-a-dozen approach to 1D Riemann solvers. Simplest HLL
! RS is used. A single, "gamma" law gas is assumed.

! "NumRS" is the actual number of 1D Riemann problems to be solved.
! "NumRS_Max" is the maximum number of 1D Riemann problems that can be fed in.

! "U_L, U_R" Left and Right states, fed in as input.

! "ent_L, ent_R" left and right entropy densities, inputs only used if
!       "PRESSURE_POSITIVITY >= 1".

! "S_L, S_R" extremal speeds, returned as outputs.

! "U_Star, F_Star" Resolved state and flux; these are the primary outputs.
! "ent_F_Star" flux of entropy density; output only provided if
!       "PRESSURE_POSITIVITY >= 1".

! "V_L, V_R" primitive variables for the left and right states; auxiliary outputs.
! "F_L, F_R" fluxes associated with the left and right states; auxiliary outputs.

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER NumRS, NumRS_Max

      REAL rhofloor, prsfloor, smallnum, gamma,
     1     fourpi, ob_4pi, ob_8pi,

     1     isospeed, barotropic_stiffness, gamma_euler

      REAL, DIMENSION ( NumRS_Max, 3) :: B0

      REAL, DIMENSION ( NumRS_Max) :: S_L, S_R,
     1      ent_L, ent_R, ent_F_Star

      REAL, DIMENSION ( NumRS_Max, NFIELD + 1 + NFLUID) ::

     1      U_L, U_R, U_Star, F_Star, V_L, V_R, F_L, F_R

! ----------------------------------------------------------------------
! Variables past this point are internal to the subroutine.

      INTEGER ii

      REAL, DIMENSION ( NumRS_Max) ::
     1     prs_L, prs_L0, del_prs_L, prs_R, prs_R0, del_prs_R,

     1     temp_01, temp_02, ob_rho,
     1     P_tot, Bf_2, Vel_2, V_dot_B, B0_dot_B,

     1     rhobar, vlxbar, vlybar, vlzbar, prsbar,
     1     bfxre, bfybar, bfzbar,

     1     va_x, ob_4pi_rho, c_s, v_a, v_f

! ----------------------------------------------------------------------

!$acc loop vector
      DO ii = 1, NumRS

! ------------------------------
! STEP 1) Impose continuity of the longitudinal magnetic field. Then
! obtain density and the best possible positive pressure we can obtain.

      bfxre ( ii) = 0.5 * ( U_L ( ii, 6) + U_R ( ii, 6) )
      U_L ( ii, 6) = bfxre ( ii)
      U_R ( ii, 6) = bfxre ( ii)

! ----------

      U_L ( ii, 1) = AMAX1 ( U_L ( ii, 1), rhofloor)

      prs_L ( ii) = (gamma - 1.0) * ( U_L ( ii, 5)
     1 - 0.5 * ( U_L ( ii, 2)**2 + U_L ( ii, 3)**2 + U_L ( ii, 4)**2)
     1       / U_L ( ii, 1)
     1 - ob_8pi  * ( U_L ( ii, 6)**2 + U_L ( ii, 7)**2
     1             + U_L ( ii, 8)**2) )

      prs_L0 ( ii) = prs_L ( ii)  ! Save pressure before positivity fix.

      del_prs_L ( ii) = 0.0  ! Increase that we may get from entropy density.

#if ( PRESSURE_POSITIVITY >= 1)

      del_prs_L ( ii) =
     1            ( 0.5 + SIGN ( 0.5, prsfloor - prs_L ( ii)) )
     1  * ( AMAX1 ( ent_L ( ii) * U_L ( ii, 1)**( gamma - 1.0),
     1              prsfloor) - prs_L ( ii) )

#endif

      prs_L ( ii) = AMAX1 ( prs_L ( ii) + del_prs_L ( ii), prsfloor)

      U_L ( ii, 5) = U_L ( ii, 5) + ( prs_L ( ii) - prs_L0 ( ii) )
     1                            / ( gamma - 1.0)

! ----------

      U_R ( ii, 1) = AMAX1 ( U_R ( ii, 1), rhofloor)

      prs_R ( ii) = (gamma - 1.0) * ( U_R ( ii, 5)
     1 - 0.5 * ( U_R ( ii, 2)**2 + U_R ( ii, 3)**2 + U_R ( ii, 4)**2)
     1       / U_R ( ii, 1)
     1 - ob_8pi  * ( U_R ( ii, 6)**2 + U_R ( ii, 7)**2
     1             + U_R ( ii, 8)**2) )

      prs_R0 ( ii) = prs_R ( ii)  ! Save pressure before positivity fix.

      del_prs_R ( ii) = 0.0  ! Increase that we may get from entropy density.

#if ( PRESSURE_POSITIVITY >= 1)

      del_prs_R ( ii) =
     1            ( 0.5 + SIGN ( 0.5, prsfloor - prs_R ( ii)) )
     1  * ( AMAX1 ( ent_R ( ii) * U_R ( ii, 1)**( gamma - 1.0),
     1              prsfloor) - prs_R ( ii) )

#endif

      prs_R ( ii) = AMAX1 ( prs_R ( ii) + del_prs_R ( ii), prsfloor)

      U_R ( ii, 5) = U_R ( ii, 5) + ( prs_R ( ii) - prs_R0 ( ii) )
     1                            / ( gamma - 1.0)

! ------------------------------
! STEP 2) Obtain Primitive variables for states "U_R, U_L" and their fluxes.

! ----------
! Focus on state "U_R". Obtain corresponding primitives "V_R" and x-flux "F_R".

      U_R ( ii, 1) = AMAX1 ( U_R ( ii, 1), rhofloor)

      ob_rho ( ii) = 1.0 / U_R ( ii, 1)

      V_R ( ii, 1) = U_R ( ii, 1)                    ! density
      V_R ( ii, 2: 4) = U_R ( ii, 2: 4) * ob_rho ( ii)  ! x,y,z-velocities
      V_R ( ii, 6: 8) = U_R ( ii, 6: 8)              ! x,y,z-B fields
#if ( NFLUID >= 1)
      V_R ( ii, 9: ) = U_R ( ii, 9: ) * ob_rho ( ii)    ! species fraction
#endif

      Vel_2 ( ii) = V_R ( ii, 2)**2 + V_R ( ii, 3)**2
     1            + V_R ( ii, 4)**2

      Bf_2 ( ii) = V_R ( ii, 6)**2 + V_R ( ii, 7)**2
     1           + V_R ( ii, 8)**2

      V_dot_B ( ii) = V_R ( ii, 2) * V_R ( ii, 6)
     1              + V_R ( ii, 3) * V_R ( ii, 7)
     1              + V_R ( ii, 4) * V_R ( ii, 8)

      B0_dot_B ( ii) = B0 ( ii, 1) * V_R ( ii, 6)
     1               + B0 ( ii, 2) * V_R ( ii, 7)
     1               + B0 ( ii, 3) * V_R ( ii, 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_R ( ii, 5) = ( gamma - 1.0) * ( U_R ( ii, 5)
     1             - 0.5 * V_R ( ii, 1) * Vel_2 ( ii)
     1             - ob_8pi * Bf_2 ( ii) )
#else
      V_R ( ii, 5) = isospeed**2 * U_R ( ii, 1)
     1             + barotropic_stiffness * U_R ( ii, 1)**gamma_euler
#endif

      V_R ( ii, 5) = AMAX1 ( V_R ( ii, 5), prsfloor)  ! gas pressure

! -----

      P_tot ( ii) = V_R ( ii, 5) + ob_8pi * Bf_2 ( ii)
     1                           + ob_4pi * B0_dot_B ( ii)

! -----

      F_R ( ii, :) = V_R ( ii, 2) * U_R ( ii, :)


      F_R ( ii, 2) = F_R ( ii, 2) + P_tot ( ii)
     1          - ob_4pi * V_R ( ii, 6) * ( V_R ( ii, 6) + B0 ( ii, 1))
     1          - ob_4pi * V_R ( ii, 6) * B0 ( ii, 1)

      F_R ( ii, 3) = F_R ( ii, 3)
     1          - ob_4pi * V_R ( ii, 6) * V_R ( ii, 7)
     1          - ( B0 ( ii, 1) * V_R ( ii, 7)
     1            + B0 ( ii, 2) * V_R ( ii, 6)) * ob_4pi

      F_R ( ii, 4) = F_R ( ii, 4)
     1          - ob_4pi * V_R ( ii, 6) * V_R ( ii, 8)
     1          - ( B0 ( ii, 1) * V_R ( ii, 8)
     1            + B0 ( ii, 3) * V_R ( ii, 6)) * ob_4pi

#if ( CHOOSE_ADIABATIC == 1)
      F_R ( ii, 5) = F_R ( ii, 5)
     1          + P_tot ( ii) * V_R ( ii, 2)
     1          - ob_4pi * ( V_R ( ii, 6) + B0 ( ii, 1))* V_dot_B ( ii)
#else
      F_R ( ii, 5) = 0.0   ! isothermal has no energy flux.
#endif

      F_R ( ii, 6) = 0.0

      F_R ( ii, 7) = F_R ( ii, 7) - V_R ( ii, 6) * V_R ( ii, 3)
     1        + B0 ( ii, 2) * V_R ( ii, 2) - B0 ( ii, 1) * V_R ( ii, 3)

      F_R ( ii, 8) = F_R ( ii, 8) - V_R ( ii, 6) * V_R ( ii, 4)
     1        + B0 ( ii, 3) * V_R ( ii, 2) - B0 ( ii, 1) * V_R ( ii, 4)

! ----------
! Focus on state "U_L". Obtain corresponding primitives "V_L" and x-flux "F_L".

      U_L ( ii, 1) = AMAX1 ( U_L ( ii, 1), rhofloor)

      ob_rho ( ii) = 1.0 / U_L ( ii, 1)

      V_L ( ii, 1) = U_L ( ii, 1)                    ! density
      V_L ( ii, 2: 4) = U_L ( ii, 2: 4) * ob_rho ( ii)  ! x,y,z-velocities
      V_L ( ii, 6: 8) = U_L ( ii, 6: 8)              ! x,y,z-B fields
#if ( NFLUID >= 1)
      V_L ( ii, 9: ) = U_L ( ii, 9: ) * ob_rho ( ii)    ! species fraction
#endif

      Vel_2 ( ii) = V_L ( ii, 2)**2 + V_L ( ii, 3)**2
     1            + V_L ( ii, 4)**2

      Bf_2 ( ii) = V_L ( ii, 6)**2 + V_L ( ii, 7)**2
     1           + V_L ( ii, 8)**2

      V_dot_B ( ii) = V_L ( ii, 2) * V_L ( ii, 6)
     1              + V_L ( ii, 3) * V_L ( ii, 7)
     1              + V_L ( ii, 4) * V_L ( ii, 8)

      B0_dot_B ( ii) = B0 ( ii, 1) * V_L ( ii, 6)
     1               + B0 ( ii, 2) * V_L ( ii, 7)
     1               + B0 ( ii, 3) * V_L ( ii, 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_L ( ii, 5) = ( gamma - 1.0) * ( U_L ( ii, 5)
     1            - 0.5 * V_L ( ii, 1) * Vel_2 ( ii)
     1            - ob_8pi * Bf_2 ( ii) )
#else
      V_L ( ii, 5) = isospeed**2 * U_L ( ii, 1)
     1             + barotropic_stiffness * U_L ( ii, 1)**gamma_euler
#endif

      V_L ( ii, 5) = AMAX1 ( V_L ( ii, 5), prsfloor)  ! gas pressure

! -----

      P_tot ( ii) = V_L ( ii, 5) + ob_8pi * Bf_2 ( ii)
     1                           + ob_4pi * B0_dot_B ( ii)

! -----

      F_L ( ii, :) = V_L ( ii, 2) * U_L ( ii, :)


      F_L ( ii, 2) = F_L ( ii, 2) + P_tot ( ii)
     1          - ob_4pi * V_L ( ii, 6) * ( V_L ( ii, 6) + B0 ( ii, 1))
     1          - ob_4pi * V_L ( ii, 6) * B0 ( ii, 1)

      F_L ( ii, 3) = F_L ( ii, 3)
     1          - ob_4pi * V_L ( ii, 6) * V_L ( ii, 7)
     1          - ( B0 ( ii, 1) * V_L ( ii, 7)
     1            + B0 ( ii, 2) * V_L ( ii, 6)) * ob_4pi

      F_L ( ii, 4) = F_L ( ii, 4)
     1          - ob_4pi * V_L ( ii, 6) * V_L ( ii, 8)
     1          - ( B0 ( ii, 1) * V_L ( ii, 8)
     1            + B0 ( ii, 3) * V_L ( ii, 6)) * ob_4pi


#if ( CHOOSE_ADIABATIC == 1)
      F_L ( ii, 5) = F_L ( ii, 5)
     1          + P_tot ( ii) * V_L ( ii, 2)
     1          - ob_4pi * ( V_L ( ii, 6) + B0 ( ii, 1)) * V_dot_B ( ii)
#else
      F_L ( ii, 5) = 0.0  ! isothermal has no energy flux.
#endif


      F_L ( ii, 6) = 0.0

      F_L ( ii, 7) = F_L ( ii, 7) - V_L ( ii, 6) * V_L ( ii, 3)
     1        + B0 ( ii, 2) * V_L ( ii, 2) - B0 ( ii, 1) * V_L ( ii, 3)

      F_L ( ii, 8) = F_L ( ii, 8) - V_L ( ii, 6) * V_L ( ii, 4)
     1        + B0 ( ii, 3) * V_L ( ii, 2) - B0 ( ii, 1) * V_L ( ii, 4)

! ------------------------------
! STEP 3) Obtain extremal wave speeds "S_R, S_L".

      S_R ( ii) = - 1.0e30
      S_L ( ii) = 1.0e30

! ----------
! Work with the pair of states "U_L, U_R" to build "S_R ( ii), S_L ( ii)"
! from the mean state. This is an inexpensive substitute for a Roe average.

      bfxre ( ii) = 0.5 * ( V_L ( ii, 6) + V_R ( ii, 6) ) + B0 ( ii, 1)

      rhobar ( ii) = 0.5 * ( V_L ( ii, 1) + V_R ( ii, 1) )

      vlxbar ( ii) = 0.5 * ( V_L ( ii, 2) + V_R ( ii, 2) )
      vlybar ( ii) = 0.5 * ( V_L ( ii, 3) + V_R ( ii, 3) )
      vlzbar ( ii) = 0.5 * ( V_L ( ii, 4) + V_R ( ii, 4) )

      prsbar ( ii) = 0.5 * ( V_L ( ii, 5) + V_R ( ii, 5) )

      bfybar ( ii) = 0.5 * ( V_L ( ii, 7) + V_R ( ii, 7) ) + B0 ( ii, 2)
      bfzbar ( ii) = 0.5 * ( V_L ( ii, 8) + V_R ( ii, 8) ) + B0 ( ii, 3)

#if ( CHOOSE_ADIABATIC != 1)
      prsbar ( ii) = isospeed**2 * rhobar ( ii)
     1       + barotropic_stiffness * rhobar ( ii)**gamma_euler
#endif

! ----------

      ob_4pi_rho ( ii) = 1.0 / ( fourpi * rhobar ( ii) )
      va_x ( ii) = bfxre ( ii)**2 * ob_4pi_rho ( ii)
      va_x ( ii) = SQRT ( va_x ( ii) )

#if ( CHOOSE_ADIABATIC == 1)
      c_s ( ii) = SQRT ( gamma * prsbar ( ii) / rhobar ( ii) )
#else
      c_s ( ii) = SQRT ( isospeed**2
     1          + gamma_euler * barotropic_stiffness
     1                        * rhobar ( ii)**( gamma_euler - 1.0) )
#endif

      v_a ( ii) = ( bfxre ( ii)**2 + bfybar ( ii)**2
     1              + bfzbar ( ii)**2) * ob_4pi_rho ( ii)
      v_a ( ii) = SQRT ( v_a ( ii) )

      temp_01 ( ii) = v_a ( ii)**2 + c_s ( ii)**2
      temp_02 ( ii) = temp_01 ( ii)**2
     1                - 4.0 * va_x ( ii)**2 * c_s ( ii)**2
      temp_02 ( ii) = AMAX1 ( temp_02 ( ii), 0.0)
      temp_02 ( ii) = SQRT ( temp_02 ( ii))

      v_f ( ii) = 0.5 * ( temp_01 ( ii) + temp_02 ( ii))
      v_f ( ii) = SQRT ( v_f ( ii) )

! ----------
! Force a subsonic wave model. That will help with the easy implementation
! of a multidimensional Riemann solver later on.

! This was needed at one time. But it is not needed any more. Besides, it actually
! helps if an HLL RS can reach the supersonic limits.

      S_L ( ii) = AMIN1 ( S_L ( ii), vlxbar ( ii) - v_f ( ii)
!     1                 , - 0.001 * v_f ( ii)  ! This line forces a subsonic wave model.
     1                  )

      S_R ( ii) = AMAX1 ( S_R ( ii), vlxbar ( ii) + v_f ( ii)
!     1                 , 0.001 * v_f ( ii)  ! This line forces a subsonic wave model.
     1                  )

! ----------
! Work with the state "U_L" to expand "S_R ( ii), S_L ( ii)".

      bfxre ( ii) = V_L ( ii, 6) + B0 ( ii, 1)

      rhobar ( ii) = V_L ( ii, 1)

      vlxbar ( ii) = V_L ( ii, 2)
      vlybar ( ii) = V_L ( ii, 3)
      vlzbar ( ii) = V_L ( ii, 4)

      prsbar ( ii) = V_L ( ii, 5)

      bfybar ( ii) = V_L ( ii, 7) + B0 ( ii, 2)
      bfzbar ( ii) = V_L ( ii, 8) + B0 ( ii, 3)

! ----------

      ob_4pi_rho ( ii) = 1.0 / ( fourpi * rhobar ( ii) )
      va_x ( ii) = bfxre ( ii)**2 * ob_4pi_rho ( ii)
      va_x ( ii) = SQRT ( va_x ( ii) )

#if ( CHOOSE_ADIABATIC == 1)
      c_s ( ii) = SQRT ( gamma * prsbar ( ii) / rhobar ( ii) )
#else
      c_s ( ii) = SQRT ( isospeed**2
     1          + gamma_euler * barotropic_stiffness
     1                        * rhobar ( ii)**( gamma_euler - 1.0) )
#endif

      v_a ( ii) = ( bfxre ( ii)**2 + bfybar ( ii)**2
     1              + bfzbar ( ii)**2) * ob_4pi_rho ( ii)
      v_a ( ii) = SQRT ( v_a ( ii) )

      temp_01 ( ii) = v_a ( ii)**2 + c_s ( ii)**2
      temp_02 ( ii) = temp_01 ( ii)**2
     1                - 4.0 * va_x ( ii)**2 * c_s ( ii)**2
      temp_02 ( ii) = AMAX1 ( temp_02 ( ii), 0.0)
      temp_02 ( ii) = SQRT ( temp_02 ( ii))

      v_f ( ii) = 0.5 * ( temp_01 ( ii) + temp_02 ( ii))
      v_f ( ii) = SQRT ( v_f ( ii) )

! ----------

      S_L ( ii) = AMIN1 ( S_L ( ii), vlxbar ( ii) - v_f ( ii))
      S_R ( ii) = AMAX1 ( S_R ( ii), vlxbar ( ii) + v_f ( ii))

! ----------
! Work with the state "U_R" to expand "S_R ( ii), S_R ( ii)".

      bfxre ( ii) = V_R ( ii, 6) + B0 ( ii, 1)

      rhobar ( ii) = V_R ( ii, 1)

      vlxbar ( ii) = V_R ( ii, 2)
      vlybar ( ii) = V_R ( ii, 3)
      vlzbar ( ii) = V_R ( ii, 4)

      prsbar ( ii) = V_R ( ii, 5)

      bfybar ( ii) = V_R ( ii, 7) + B0 ( ii, 2)
      bfzbar ( ii) = V_R ( ii, 8) + B0 ( ii, 3)

! ----------

      ob_4pi_rho ( ii) = 1.0 / ( fourpi * rhobar ( ii) )
      va_x ( ii) = bfxre ( ii)**2 * ob_4pi_rho ( ii)
      va_x ( ii) = SQRT ( va_x ( ii) )

#if ( CHOOSE_ADIABATIC == 1)
      c_s ( ii) = SQRT ( gamma * prsbar ( ii) / rhobar ( ii) )
#else
      c_s ( ii) = SQRT ( isospeed**2
     1          + gamma_euler * barotropic_stiffness
     1                        * rhobar ( ii)**( gamma_euler - 1.0) )
#endif

      v_a ( ii) = ( bfxre ( ii)**2 + bfybar ( ii)**2
     1              + bfzbar ( ii)**2) * ob_4pi_rho ( ii)
      v_a ( ii) = SQRT ( v_a ( ii) )

      temp_01 ( ii) = v_a ( ii)**2 + c_s ( ii)**2
      temp_02 ( ii) = temp_01 ( ii)**2
     1                - 4.0 * va_x ( ii)**2 * c_s ( ii)**2
      temp_02 ( ii) = AMAX1 ( temp_02 ( ii), 0.0)
      temp_02 ( ii) = SQRT ( temp_02 ( ii))

      v_f ( ii) = 0.5 * ( temp_01 ( ii) + temp_02 ( ii))
      v_f ( ii) = SQRT ( v_f ( ii) )

! ----------

      S_L ( ii) = AMIN1 ( S_L ( ii), vlxbar ( ii) - v_f ( ii))
      S_R ( ii) = AMAX1 ( S_R ( ii), vlxbar ( ii) + v_f ( ii))

! ------------------------------
! STEP 4) Form the HLL flux for the conserved variables and for the entropy.
! The logic provided here already enables us to pick out the supersonic limits.

      IF ( S_L ( ii) >= 0.0) THEN

        U_Star ( ii, :) = U_L ( ii, :)
        F_Star ( ii, :) = F_L ( ii, :)
        ent_F_Star ( ii) = ent_L ( ii) * V_L ( ii, 2)

      ELSE IF ( S_R ( ii) <= 0.0) THEN

        U_Star ( ii, :) = U_R ( ii, :)
        F_Star ( ii, :) = F_R ( ii, :)
        ent_F_Star ( ii) = ent_R ( ii) * V_R ( ii, 2)

      ELSE

        temp_01 ( ii) = 1.0 / ( S_R ( ii) - S_L ( ii))
        temp_02 ( ii) = S_R ( ii) * S_L ( ii) * temp_01 ( ii)

        U_Star ( ii, :) =
     1  ( S_R ( ii) * U_R ( ii, :) - S_L ( ii) * U_L ( ii, :)
     1  - F_R ( ii, :) + F_L ( ii, :) ) * temp_01 ( ii)

        F_Star ( ii, :) =
     1  ( S_R ( ii) * F_L ( ii, :)
     1  - S_L ( ii) * F_R ( ii, :) ) * temp_01 ( ii)
     1  + ( U_R ( ii, :) - U_L ( ii, :) ) * temp_02 ( ii)

        ent_F_Star ( ii) = ( S_R ( ii) * ent_L ( ii) * V_L ( ii, 2)
     1                      - S_L ( ii) * ent_R ( ii) * V_R ( ii, 2))
     1                    * temp_01 ( ii)
     1  + ( ent_R ( ii) - ent_L ( ii) ) * temp_02 ( ii)

      END IF

! ------------------------------

      END DO  ! End of "DO ii = 1, NumRS" loop.

      END SUBROUTINE OneD_RS_MHD_DAD



