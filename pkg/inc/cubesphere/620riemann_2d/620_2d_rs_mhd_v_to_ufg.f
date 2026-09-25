







      SUBROUTINE MHD_V_TO_UFG ( n, gamma, B0, V, U, F, G )

! ----------------------------------------------------------------------
! Converts primitive variables "V" to conserved variables "U".
!$acc routine seq

      IMPLICIT NONE

      INTEGER n

      REAL gamma, B0 ( 3)

      REAL, DIMENSION ( n) :: V, U, F, G

      REAL v2, B2, V_dot_B, B0_dot_B, P_tot, ob_4pi, ob_8pi, pi

! ----------------------------------------------------------------------

      pi = 4.0 * ATAN ( 1.0)

      ob_4pi = 1.0 / ( 4.0 * pi)
      ob_8pi = 0.5 * ob_4pi

! ----------------------------------------------------------------------

      U ( 1) = V ( 1)

      U ( 2) = V ( 1) * V ( 2)
      U ( 3) = V ( 1) * V ( 3)
      U ( 4) = V ( 1) * V ( 4)

#if ( CHOOSE_ADIABATIC == 1)
      U ( 5) = V ( 5) / ( gamma - 1.0)
     1       + 0.5 * V ( 1) * ( V ( 2)**2 + V ( 3)**2 + V ( 4)**2)
     1             + ob_8pi * ( V ( 6)**2 + V ( 7)**2 + V ( 8)**2)
#else
      U ( 5) = 0.0
#endif

      U ( 6) = V ( 6)
      U ( 7) = V ( 7)
      U ( 8) = V ( 8)

! ----------------------------------------
! Building "F, G".

      v2 = V ( 2)**2 + V ( 3)**2 + V ( 4)**2
      B2 = V ( 6)**2 + V ( 7)**2 + V ( 8)**2

      V_dot_B  = V ( 2) * V ( 6)
     1         + V ( 3) * V ( 7)
     1         + V ( 4) * V ( 8)

      B0_dot_B  = B0 ( 1) * V ( 6)
     1          + B0 ( 2) * V ( 7)
     1          + B0 ( 3) * V ( 8)

      P_tot  = V ( 5) + ob_8pi * B2
     1                + ob_4pi * B0_dot_B

! ----------
! Building "F".

      F ( :) = V ( 2) * U ( :)


      F ( 2) = F ( 2) + P_tot
     1          - ob_4pi * V ( 6) * ( V ( 6) + B0 ( 1))
     1          - ob_4pi * V ( 6) * B0 ( 1)

      F ( 3) = F ( 3)
     1          - ob_4pi * V ( 6) * V ( 7)
     1          - ob_4pi * ( B0 ( 1) * V ( 7)
     1                     + B0 ( 2) * V ( 6))

      F ( 4) = F ( 4)
     1          - ob_4pi * V ( 6) * V ( 8)
     1          - ob_4pi * ( B0 ( 1) * V ( 8)
     1                     + B0 ( 3) * V ( 6))

#if ( CHOOSE_ADIABATIC == 1)
      F ( 5) = F ( 5)
     1          + P_tot  * V ( 2)
     1          - ob_4pi * ( V ( 6) + B0 ( 1))* V_dot_B
#else
      F ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      F ( 6) = 0.0

      F ( 7) = F ( 7) - V ( 6) * V ( 3)
     1       + B0 ( 2) * V ( 2) - B0 ( 1) * V ( 3)

      F ( 8) = F ( 8) - V ( 6) * V ( 4)
     1       + B0 ( 3) * V ( 2) - B0 ( 1) * V ( 4)

! ----------
! Building "G".

      G ( :) = V ( 3) * U ( :)

      G ( 2) = G ( 2)
     1          - ob_4pi * V ( 6) * V ( 7)
     1          - ob_4pi * ( B0 ( 1) * V ( 7)
     1                     + B0 ( 2) * V ( 6))

      G ( 3) = G ( 3) + P_tot
     1           - ob_4pi * V ( 7) * ( V ( 7) + B0 ( 2))
     1           - ob_4pi * V ( 7) * B0 ( 2)

      G ( 4) = G ( 4)
     1          - ob_4pi * V ( 7) * V ( 8)
     1          - ob_4pi * ( B0 ( 2) * V ( 8)
     1                     + B0 ( 3) * V ( 7))

#if ( CHOOSE_ADIABATIC == 1)
      G ( 5) = G ( 5)
     1          + P_tot  * V ( 3)
     1          - ob_4pi * ( V ( 7) + B0 ( 2))* V_dot_B
#else
      G ( 5) = 0.0             ! isothermal has no energy flux.
#endif

      G ( 6) = G ( 6) - V ( 2) * V ( 7)
     1       + V ( 3) * B0 ( 1) - V ( 2) * B0 ( 2)

      G ( 7) = 0.0

      G ( 8) = G ( 8) - V ( 4) * V ( 7)
     1       + V ( 3) * B0 ( 3) - V ( 4) * B0 ( 2)

! ----------------------------------------------------------------------

      END SUBROUTINE MHD_V_TO_UFG



