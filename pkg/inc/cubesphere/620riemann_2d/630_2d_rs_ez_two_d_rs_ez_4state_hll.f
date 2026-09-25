








      SUBROUTINE Two_D_RS_Ez_4State_HLL

     1         ( S_R, S_L, S_U, S_D, cos_t, sin_t,

     1           Ez_RU, Ez_LU, Ez_LD, Ez_RD,

     1           Bneta_U, Bneta_D, By_R, By_L,

     1           BSSx, BSSy, ESSz)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! This is the simplest HLL Multidimensional Riemann Solver that is designed to handle
! the induction equation and equations like it on a skew Cartesian mesh with 4 input states.
! Such equations require the invocation of a 2D Riemann solver and they
! greatly benefit if the Riemann solver has a well-defined supersonic limit.
! While most Riemann solvers operate on a whole system, this one
! operates on just the induction equation:-

!   d/dt ( Bx, By, Bz) + d/dx ( 0, -Ez, Ey) + d/dy ( Ez, 0, -Ex) = 0

! The extremal speeds "S_R, S_L, S_U, S_D" in both directions are supplied externally.
! They can be obtained from suitable 1D Riemann solvers applied priorly
! in the x- and y-directions.

! The x-axis is assumed given. "cos_t, sin_t" give the cos and sine of the angle
! between the x-axis and the skewed y-axis. "sin_t" can never be zero.

! Likewise, the z-component of the electric fields "Ez_RU, Ez_LU, Ez_LD, Ez_RD"
! are evaluated at zone centers and interpolated (not reconstructed) to the edge centers.

! The dissipation is proportional to "Bneta_U - Bneta_D" and "By_R - By_L". Therefore,
! higher order interpolation is applied in the faces that abut the edge and the
! magnetic fields "Bneta_U, Bneta_D, By_R, By_L" are obtained directly from those interpolations.

! The output is the resolved x- and y-components of the magnetic field "BSSx, BSSy";
! and most importantly the z-component of the electric field "ESSz".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL S_R, S_L, S_U, S_D, cos_t, sin_t,

     1     Ez_RU, Ez_LU, Ez_LD, Ez_RD,

     1     Bneta_U, Bneta_D, By_R, By_L,

     1     BSSx, BSSy, ESSz,

     1     obS_x, obS_y, BSneta, BStau, ESSz_1, ESSz_2, x_c, y_c,
     1     csc_t, cot_t, S, obS,
     1     x1, y1, x2, y2, x3, y3, x4, y4,
     1     x12, y12, x13, y13, x43, y43, x42, y42,

     1     Bx_RU, By_RU, Bneta_RU, Btau_RU,
     1     Bx_LU, By_LU, Bneta_LU, Btau_LU,
     1     Bx_LD, By_LD, Bneta_LD, Btau_LD,
     1     Bx_RD, By_RD, Bneta_RD, Btau_RD,

     1     BSx_R, BSy_R, ESz_R, BSx_L, BSy_L, ESz_L,
     1     BSx_U, BSy_U, ESz_U, BSx_D, BSy_D, ESz_D

! ----------------------------------------------------------------------

        obS_x = 1.0 / (S_R - S_L)
        obS_y = 1.0 / (S_U - S_D)

        csc_t = 1.0 / sin_t
        cot_t = cos_t / sin_t

        S = AMAX1 ( ABS ( S_R), ABS ( S_L), ABS ( S_U), ABS ( S_D) )
        obS = 1.0 / S

! ----------------------------------------------------------------------
! Set up the states in all of the four quadrants.

        Bx_RU = (Bneta_U + By_R * cos_t) / sin_t  ;  By_RU = By_R
        Bx_LU = (Bneta_U + By_L * cos_t) / sin_t  ;  By_LU = By_L
        Bx_RD = (Bneta_D + By_R * cos_t) / sin_t  ;  By_RD = By_R
        Bx_LD = (Bneta_D + By_L * cos_t) / sin_t  ;  By_LD = By_L

        Btau_RU = Bx_RU * cos_t + By_RU * sin_t  ;  Bneta_RU = Bneta_U
        Btau_LU = Bx_LU * cos_t + By_LU * sin_t  ;  Bneta_LU = Bneta_U
        Btau_RD = Bx_RD * cos_t + By_RD * sin_t  ;  Bneta_RD = Bneta_D
        Btau_LD = Bx_LD * cos_t + By_LD * sin_t  ;  Bneta_LD = Bneta_D

! ----------------------------------------------------------------------

      IF ( ( S_L .GE. 0.0) .AND. ( S_D .GE. 0.0)) THEN
! Left-Down state is supersonic in both directions.

! ------------------------------

        BSSx = (Bneta_D + cos_t * By_L) / sin_t
        BSSy = By_L
        ESSz = Ez_LD

        RETURN

! ------------------------------

      ELSE IF ( ( S_R .LE. 0.0) .AND. ( S_D .GE. 0.0)) THEN
! Right-Down state is supersonic in both directions.

! ------------------------------

        BSSx = (Bneta_D + cos_t * By_R) / sin_t
        BSSy = By_R
        ESSz = Ez_RD

        RETURN

! ------------------------------

      ELSE IF ( ( S_R .LE. 0.0) .AND. ( S_U .LE. 0.0)) THEN
! Right-Up state is supersonic in both directions.

! ------------------------------

        BSSx = (Bneta_U + cos_t * By_R) / sin_t
        BSSy = By_R
        ESSz = Ez_RU

        RETURN

! ------------------------------

      ELSE IF ( ( S_L .GE. 0.0) .AND. ( S_U .LE. 0.0)) THEN
! Left-Up state is supersonic in both directions.

! ------------------------------

        BSSx = (Bneta_U + cos_t * By_L) / sin_t
        BSSy = By_L
        ESSz = Ez_LU

        RETURN

! ------------------------------

      ELSE IF ( S_L .GE. 0.0) THEN
! Left-resolved state is supersonic in x-direction, subsonic in y-direction.

! ------------------------------

        BSx_L = (S_U * Bx_LU - S_D * Bx_LD) * obS_y
     1        - (Ez_LU - Ez_LD) * obS_y

        BSy_L = By_L

        ESz_L = (S_U * Ez_LD - S_D * Ez_LU) * obS_y
     1        + S_U * S_D * (Bx_LU - Bx_LD) * obS_y

        BSSx = BSx_L  ;  BSSy = BSy_L  ;  ESSz = ESz_L

        RETURN

! ------------------------------

      ELSE IF ( S_R .LE. 0.0) THEN
! Right-resolved state is supersonic in x-direction, subsonic in y-direction.

! ------------------------------

        BSx_R = (S_U * Bx_RU - S_D * Bx_RD) * obS_y
     1        - (Ez_RU - Ez_RD) * obS_y

        BSy_R = By_R

        ESz_R = (S_U * Ez_RD - S_D * Ez_RU) * obS_y
     1        + S_U * S_D * (Bx_RU - Bx_RD) * obS_y

        BSSx = BSx_R  ;  BSSy = BSy_R  ;  ESSz = ESz_R

        RETURN

! ------------------------------

      ELSE IF ( S_D .GE. 0.0) THEN
! Down-resolved state is supersonic in y-direction, subsonic in x-direction.

! ------------------------------

        BSneta = Bneta_D

        BStau = (S_R * Btau_RD - S_L * Btau_LD) * obS_x
     1        + (Ez_RD - Ez_LD) * obS_x

        BSx_D = BSneta * sin_t + BStau * cos_t

        BSy_D = - BSneta * cos_t + BStau * sin_t

        ESz_D = (S_R * Ez_LD - S_L * Ez_RD) * obS_x
     1        - S_R * S_L * (Btau_RD - Btau_LD) * obS_x

        BSSx = BSx_D  ;  BSSy = BSy_D  ;  ESSz = ESz_D

        RETURN

! ------------------------------

      ELSE IF ( S_U .LE. 0.0) THEN
! Up-resolved state is supersonic in y-direction, subsonic in x-direction.

! ------------------------------

        BSneta = Bneta_U

        BStau = (S_R * Btau_RU - S_L * Btau_LU) * obS_x
     1        + (Ez_RU - Ez_LU) * obS_x

        BSx_U = BSneta * sin_t + BStau * cos_t

        BSy_U = - BSneta * cos_t + BStau * sin_t

        ESz_U = (S_R * Ez_LU - S_L * Ez_RU) * obS_x
     1        - S_R * S_L * (Btau_RU - Btau_LU) * obS_x

        BSSx = BSx_U  ;  BSSy = BSy_U  ;  ESSz = ESz_U

        RETURN

! ------------------------------

      ELSE
! Fully subsonic state. This is the true 2D Riemann solver.

! ------------------------------

! ----------

! Left-resolved state.

        BSx_L = (S_U * Bx_LU - S_D * Bx_LD) * obS_y
     1        - (Ez_LU - Ez_LD) * obS_y

        BSy_L = By_L

        ESz_L = (S_U * Ez_LD - S_D * Ez_LU) * obS_y
     1        + S_U * S_D * (Bx_LU - Bx_LD) * obS_y

! For quick comparison/debugging, we also provide the LLF states; but we comment them out.
!        BSx_L = 0.5 * (Bx_LU + Bx_LD) - 0.5 * obS * (Ez_LU - Ez_LD)
!        BSy_L = By_L
!        ESz_L = 0.5 * (Ez_LD + Ez_LU) - 0.5 * S * (Bx_LU - Bx_LD)

! ----------

! Right-resolved state.

        BSx_R = (S_U * Bx_RU - S_D * Bx_RD) * obS_y
     1        - (Ez_RU - Ez_RD) * obS_y

        BSy_R = By_R

        ESz_R = (S_U * Ez_RD - S_D * Ez_RU) * obS_y
     1        + S_U * S_D * (Bx_RU - Bx_RD) * obS_y

! For quick comparison/debugging, we also provide the LLF states; but we comment them out.
!        BSx_R = 0.5 * (Bx_RU + Bx_RD) - 0.5 * obS * (Ez_RU - Ez_RD)
!        BSy_R = By_R
!        ESz_R = 0.5 * (Ez_RD + Ez_RU) - 0.5 * S * (Bx_RU - Bx_RD)

! ----------

! Down-resolved state.

        BSneta = Bneta_D

        BStau = (S_R * Btau_RD - S_L * Btau_LD) * obS_x
     1        + (Ez_RD - Ez_LD) * obS_x

        BSx_D = BSneta * sin_t + BStau * cos_t

        BSy_D = - BSneta * cos_t + BStau * sin_t

        ESz_D = (S_R * Ez_LD - S_L * Ez_RD) * obS_x
     1        - S_R * S_L * (Btau_RD - Btau_LD) * obS_x

! For quick comparison/debugging, we also provide the LLF states; but we comment them out.
!        BSneta = Bneta_D
!        BStau = 0.5 * ( Btau_RD + Btau_LD) + 0.5 * obS * (Ez_RD - Ez_LD)
!        ESz_D = 0.5 * ( Ez_LD + Ez_RD) + 0.5 * S * (Btau_RD - Btau_LD)
!        BSx_D = BSneta * sin_t + BStau * cos_t
!        BSy_D = - BSneta * cos_t + BStau * sin_t

! ----------

! Up-resolved state.

        BSneta = Bneta_U

        BStau = (S_R * Btau_RU - S_L * Btau_LU) * obS_x
     1        + (Ez_RU - Ez_LU) * obS_x

        BSx_U = BSneta * sin_t + BStau * cos_t

        BSy_U = - BSneta * cos_t + BStau * sin_t

        ESz_U = (S_R * Ez_LU - S_L * Ez_RU) * obS_x
     1        - S_R * S_L * (Btau_RU - Btau_LU) * obS_x

! For quick comparison/debugging, we also provide the LLF states; but we comment them out.
!        BSneta = Bneta_U
!        BStau = 0.5 * ( Btau_RU + Btau_LU) + 0.5 * obS * (Ez_RU - Ez_LU)
!        ESz_U = 0.5 * ( Ez_LU + Ez_RU) + 0.5 * S * (Btau_RU - Btau_LU)
!        BSx_U = BSneta * sin_t + BStau * cos_t
!        BSy_U = - BSneta * cos_t + BStau * sin_t

! ----------
! Build "x1, y1", "x2, y2", "x3, y3", "x4, y4" and "xc, yc". Again needed only for debugging.

!        x1 = (S_U * cot_t + S_R * csc_t)  ;  y1 = S_U  ! RU vertex of wave model.
!        x2 = (S_D * cot_t + S_R * csc_t)  ;  y2 = S_D  ! RD vertex of wave model.
!        x3 = (S_U * cot_t + S_L * csc_t)  ;  y3 = S_U  ! LU vertex of wave model.
!        x4 = (S_D * cot_t + S_L * csc_t)  ;  y4 = S_D  ! LD vertex of wave model.

!        xc = 0.5 * (x1 + x4)  ;  yc = 0.5 * (y1 + y4)

! Build centroids of the boundaries, relative to (xc, yc). Again needed only for debugging.

!        x12 = 0.5 * (x1 + x2) - xc  ;  y12 = 0.5 * (y1 + y2) - yc  ! Center of R boundary.
!        x13 = 0.5 * (x1 + x3) - xc  ;  y13 = 0.5 * (y1 + y3) - yc  ! Center of U boundary.
!        x43 = 0.5 * (x4 + x3) - xc  ;  y43 = 0.5 * (y4 + y3) - yc  ! Center of L boundary.
!        x42 = 0.5 * (x4 + x2) - xc  ;  y42 = 0.5 * (y4 + y2) - yc  ! Center of D boundary.

! Now use the above four solutions to the 1D Riemann problems to build the fully
! subsonic states and z-component of electric field.

        x_c = 0.5 * ( (S_U + S_D) * cos_t + (S_R + S_L)) / sin_t

        y_c = 0.5 * (S_U + S_D)

        BSSx = - 0.5 * ( ( - cos_t * ESz_R - S_R * BSx_R) * obS_x
     1                 + ( ESz_U - S_U * BSx_U) * obS_y
     1                 + ( cos_t * ESz_L + S_L * BSx_L) * obS_x
     1                 + ( - ESz_D + S_D * BSx_D) * obS_y )

        BSSy = - 0.5 * ( ( - sin_t * ESz_R - S_R * BSy_R) * obS_x
     1                 + ( - S_U * BSy_U) * obS_y
     1                 + ( sin_t * ESz_L + S_L * BSy_L) * obS_x
     1                 + ( S_D * BSy_D) * obS_y )

        ESSz_1 = - x_c * BSSy
     1         - 0.5 * ( csc_t * ( - sin_t * ESz_R - S_R * BSy_R)
     1                 + cot_t * ( - S_U * BSy_U)
     1                 - csc_t * ( sin_t * ESz_L + S_L * BSy_L)
     1                 - cot_t * ( S_D * BSy_D) )

        ESSz_2 = y_c * BSSx + 0.5 * ( ( ESz_U - S_U * BSx_U)
     1                              - ( - ESz_D + S_D * BSx_D) )

        ESSz = 0.5 * ( ESSz_1 + ESSz_2)

        RETURN

      END IF


      END SUBROUTINE Two_D_RS_Ez_4State_HLL





