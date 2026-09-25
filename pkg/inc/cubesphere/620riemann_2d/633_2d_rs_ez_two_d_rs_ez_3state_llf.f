











      SUBROUTINE Two_D_RS_Ez_3State_LLF

     1         ( S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,

     1           Ez_1, Ez_2, Ez_3,

     1           Bneta_1, Bneta_2, Bneta_3,

     1           BSSx, BSSy, ESSz)

! ----------------------------------------------------------------------

! This is the simplest LLF Multidimensional Riemann Solver that is designed to handle
! the induction equation and equations like it on a tri-corner mesh with 3 input states.
! LLF only has a subsonic limit; but it has simple, closed-form expressions.
! While most Riemann solvers operate on a whole system, this one
! operates on just the induction equation:-

!   d/dt ( Bx, By, Bz) + d/dx ( 0, -Ez, Ey) + d/dy ( Ez, 0, -Ex) = 0

! We can think of the tri-corner mesh as a 3 state mesh with one axis along
! the x-axis; the other at 120 degrees to it; and the third at 240 degrees to it.

! The extremal speeds "S_R1, S_L1, S_R2, S_L2, S_R3, S_L3" in both directions
! are supplied externally. They can be obtained from suitable 1D Riemann solvers
! applied priorly to the three faces of the mesh.

! Likewise, the z-component of the electric fields "Ez_1, Ez_2, Ez_3"
! are evaluated at zone centers and interpolated (not reconstructed) at the edge centers.

! The dissipation is proportional to

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,

     1     Ez_1, Ez_2, Ez_3,

     1     Bneta_1, Bneta_2, Bneta_3,

     1     BSSx, BSSy, ESSz

! ----------------------------------------------------------------------
! Variables past this point are internal to this subroutine.

      REAL S, obS, Sq3, obSq3, ob3

! ----------------------------------------------------------------------

      S = AMAX1 ( ABS ( S_R1), ABS ( S_L1), ABS ( S_R2), ABS ( S_L2),
     1            ABS ( S_R3), ABS ( S_L3) )

      obS = 1.0 / S

      Sq3 = SQRT ( 3.0)
      obSq3 = 1.0 / Sq3
      ob3 = 1.0 / 3.0

! ----------------------------------------------------------------------

      BSSx = obSq3 * ( Bneta_2 - Bneta_3) + 0.5 * obS * ( Ez_3 - Ez_1)

      BSSy = ob3 * ( Bneta_2 + Bneta_3 - 2.0 * Bneta_1)
     1     + ( 0.5 * obSq3 * obS) * ( Ez_1 - 2.0 * Ez_2 + Ez_3)

      ESSz = ob3 * ( Ez_1 + Ez_2 + Ez_3)
     1     - ( 4.0 * Sq3 * S / 9.0) * ( Bneta_1 + Bneta_2 + Bneta_3)

! ----------------------------------------------------------------------

      END SUBROUTINE Two_D_RS_Ez_3State_LLF



