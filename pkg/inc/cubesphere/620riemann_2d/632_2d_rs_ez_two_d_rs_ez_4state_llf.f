








      SUBROUTINE Two_D_RS_Ez_4State_LLF

     1         ( S_R, S_L, S_U, S_D, cos_t, sin_t,

     1           Ez_RU, Ez_LU, Ez_LD, Ez_RD,

     1           Bneta_U, Bneta_D, By_R, By_L,

     1           BSSx, BSSy, ESSz)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! This is the simplest LLF Multidimensional Riemann Solver that is designed to handle
! the induction equation and equations like it on a skew Cartesian mesh with 4 input states.
! LLF only has a subsonic limit; but it has simple, closed-form expressions.
! While most Riemann solvers operate on a whole system, this one
! operates on just the induction equation:-

!   d/dt ( Bx, By, Bz) + d/dx ( 0, -Ez, Ey) + d/dy ( Ez, 0, -Ex) = 0

! The extremal speeds "S_R, S_L, S_U, S_D" in both directions are supplied externally.
! They can be obtained from suitable 1D Riemann solvers applied priorly
! in the x- and y-directions. Inside this subroutine, we take a local max of all 4 speeds.

! The x-axis is assumed given. "cos_t, sin_t" give the cos and sine of the angle
! between the x-axis and the skewed y-axis. "sin_t" can never be zero.

! Likewise, the z-component of the electric fields "Ez_RU, Ez_LU, Ez_LD, Ez_RD"
! are evaluated at zone centers and interpolated (not reconstructed) to the edge centers.

! The dissipation will just be proportional to the circulation of the magnetic field
! around that vertex.

! The output is the resolved x- and y-components of the magnetic field "BSSx, BSSy";
! and most importantly the z-component of the electric field "ESSz".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL S_R, S_L, S_U, S_D, cos_t, sin_t,

     1     Ez_RU, Ez_LU, Ez_LD, Ez_RD,

     1     Bneta_U, Bneta_D, By_R, By_L,

     1     BSSx, BSSy, ESSz

! ----------------------------------------------------------------------
! Variables past this point are internal to this subroutine.

      REAL csc_t, cot_t, S, obS

! ----------------------------------------------------------------------

        S = AMAX1 ( ABS ( S_R), ABS ( S_L), ABS ( S_U), ABS ( S_D) )
        obS = 1.0 / S

        csc_t = 1.0 / sin_t
        cot_t = cos_t / sin_t

! ----------------------------------------------------------------------

        BSSx = ( 0.5 * csc_t)
     1       * ( cos_t * ( By_R + By_L) + Bneta_U + Bneta_D)
     1       + ( 0.25 * obS * csc_t)
     1       * ( sin_t * ( Ez_LD - Ez_LU + Ez_RD - Ez_RU)
     1         - sin_t * cos_t * ( Ez_LD + Ez_LU - Ez_RD - Ez_RU) )

        BSSy = 0.5 * ( By_R + By_L)
     1       - ( 0.25 * obS * sin_t) * ( Ez_LD + Ez_LU - Ez_RD - Ez_RU)

        ESSz = 0.25 * ( Ez_LD + Ez_LU + Ez_RD + Ez_RU)
     1       + ( 0.5 * S * csc_t) * ( Bneta_D - Bneta_U + By_R - By_L)

! ----------------------------------------------------------------------

      END SUBROUTINE Two_D_RS_Ez_4State_LLF



