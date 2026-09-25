







      SUBROUTINE Two_D_RS_Ez_3State_HLL

     1         ( S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,

     1           Ez_1, Ez_2, Ez_3,

     1           Bneta_1, Bneta_2, Bneta_3,

     1           BSSx, BSSy, ESSz)

! ----------------------------------------------------------------------

! This is the simplest HLL Multidimensional Riemann Solver that is designed to handle
! the induction equation and equations like it on a tri-corner mesh with 3 input states.
! Such equations require the invocation of a 2D Riemann solver and they
! greatly benefit if the Riemann solver has a well-defined supersonic limit.
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

! The dissipation is proportional to "Btau1_U3 - Btau1_U1" and
! "Btau2_U1 - Btau2_U2" and "Btau3_U2 - Btau3_U3".

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL S_R1, S_L1, S_R2, S_L2, S_R3, S_L3,

     1     Ez_1, Ez_2, Ez_3,

     1     Bneta_1, Bneta_2, Bneta_3,

     1     BSSx, BSSy, ESSz

! ----------------------------------------------------------------------
! local variables.

      REAL Btau1_U3, Btau1_U1, Btau2_U1, Btau2_U2,
     1     Btau3_U2, Btau3_U3, sq3_2, obsq3,
     1     obS_1, obS_2, obS_3,

     1     Bx1, Bx2, Bx3, By1, By2, By3,

     1     xc, yc, Ax, Ay, Bx, By, Cx, Cy,

     1     l1_x, l2_x, l3_x, l1_y, l2_y, l3_y,
     1     l1, l2, l3, Area, s,

     1     X13_x, X21_x, X32_x,
     1     X13_y, X21_y, X32_y,

     1     S13, S21, S32, xb13, xb21, xb32, yb13, yb21, yb32


      REAL BSneta_13, BSneta_21, BSneta_32,
     1     BStau_13, BStau_21, BStau_32, ESz_13, ESz_21, ESz_32,

     1     BSx_13, BSx_21, BSx_32, BSy_13, BSy_21, BSy_32,
     1     ESSz_1, ESSz_2

! ----------------------------------------------------------------------

      sq3_2 = SQRT ( 3.0) / 2.0
      obsq3 = 1.0 / SQRT ( 3.0)

      obS_1 = 1.0 / ( S_R1 - S_L1)
      obS_2 = 1.0 / ( S_R2 - S_L2)
      obS_3 = 1.0 / ( S_R3 - S_L3)

! --------------------
! vertices of triangular SI state.

      Ax =  ( S_L1 + 2.0 * S_R2) * obsq3 ;     Ay = - S_L1
      Bx =  ( S_L2 - S_R3) * obsq3 ;     By = S_L2 + S_R3
      Cx = -( S_R1 + 2.0 * S_L3) * obsq3 ;     Cy = - S_R1

      xc = ( Ax + Bx + Cx) / 3.0
      yc = ( Ay + By + Cy) / 3.0

! --------------------
! lengths of triangular SI state.

      l1_x = Ax - Cx ;   l2_x = Bx - Ax ;   l3_x = Cx - Bx
      l1_y = Ay - Cy ;   l2_y = By - Ay ;   l3_y = Cy - By

      l1 = SQRT ( l1_x**2 + l1_y**2)
      l2 = SQRT ( l2_x**2 + l2_y**2)
      l3 = SQRT ( l3_x**2 + l3_y**2)

! --------------------
! unit normals from the triangle. It is calculated by rotating the
! unit vectors of the edges of the triangle. The edge vectors are rotated
! by -90 degrees.

      X13_x = l1_y / l1 ;     X13_y = - l1_x / l1
      X21_x = l2_y / l2 ;     X21_y = - l2_x / l2
      X32_x = l3_y / l3 ;     X32_y = - l3_x / l3

! ----------
! Speeds of the triangular boundaries.

      S13 = X13_x * Ax + X13_y * Ay
      S21 = X21_x * Bx + X21_y * By
      S32 = X32_x * Cx + X32_y * Cy

! ----------
! midpoints of vectors.

      xb13 = 0.5 * ( Ax + Cx) ;     yb13 = 0.5 * ( Ay + Cy)
      xb21 = 0.5 * ( Bx + Ax) ;     yb21 = 0.5 * ( By + Ay)
      xb32 = 0.5 * ( Cx + Bx) ;     yb32 = 0.5 * ( Cy + By)


! --------------------
! Area of the Triangle.

      s = ( l1 + l2 + l3) / 2.0

      Area = SQRT ( s * ( s - l1) * ( s - l2) * ( s - l3))

! ----------------------------------------------------------------------
! xy - components of magnetic field.

      Bx1 =  ( 2.0 * Bneta_2 + Bneta_1) * obsq3
      Bx2 =  ( Bneta_2 - Bneta_3) * obsq3
      Bx3 = -( 2.0 * Bneta_3 + Bneta_1) * obsq3

      By1 = - Bneta_1
      By2 = Bneta_2 + Bneta_3
      By3 = - Bneta_1

! --------------------
! tau component of the magnetic field.

      Btau1_U1 = Bx1
      Btau1_U3 = Bx3

      Btau2_U1 = - 0.5 * Bx1 + sq3_2 * By1
      Btau2_U2 = - 0.5 * Bx2 + sq3_2 * By2

      Btau3_U2 = - 0.5 * Bx2 - sq3_2 * By2
      Btau3_U3 = - 0.5 * Bx3 - sq3_2 * By3

! ----------------------------------------------------------------------
! Building the magnetic field components to compute, 1D resolved quantities
! at the zone boundary - 1.

      BSneta_13 = Bneta_1          ! normal magnetic field

      BStau_13 = ( S_R1 * Btau1_U3
     1           - S_L1 * Btau1_U1 + Ez_3 - Ez_1) * obS_1  ! tangential magnetic field

      ESz_13 = ( S_R1 * Ez_1 - S_L1 * Ez_3
     1         - S_R1 * S_L1 * (Btau1_U3 - Btau1_U1)) * obS_1   ! 1D resolved e-field

      BSx_13 = BStau_13                    ! 1D resolved magnetic field-x
      BSy_13 = - BSneta_13                 ! 1D resolved magnetic field-y


! ----------------------------------------
! Building the magnetic field components to compute, 1D resolved quantities
! at the zone boundary - 2.

      BSneta_21 = Bneta_2          ! normal magnetic field

      BStau_21 = ( S_R2 * Btau2_U1
     1           - S_L2 * Btau2_U2 + Ez_1 - Ez_2) * obS_2  ! tangential magnetic field

      ESz_21 = ( S_R2 * Ez_2 - S_L2 * Ez_1
     1         - S_R2 * S_L2 * (Btau2_U1 - Btau2_U2)) * obS_2  ! 1D resolved e-field

      BSx_21 = - 0.5 * BStau_21 + sq3_2 * BSneta_21  ! 1D resolved magnetic field-x
      BSy_21 = sq3_2 * BStau_21 + 0.5 * BSneta_21    ! 1D resolved magnetic field-y


! ----------------------------------------
! Building the magnetic field components to compute, 1D resolved quantities
! at the zone boundary - 3.

      BSneta_32 = Bneta_3          ! normal magnetic field

      BStau_32 = ( S_R3 * Btau3_U2
     1           - S_L3 * Btau3_U3 + Ez_2 - Ez_3) * obS_3  ! tangential magnetic field

      ESz_32 = ( S_R3 * Ez_3 - S_L3 * Ez_2
     1         - S_R3 * S_L3 * ( Btau3_U2 - Btau3_U3)) * obS_3  ! 1D resolved e-field

      BSx_32 = - 0.5 * BStau_32 - sq3_2 * BSneta_32   ! 1D resolved magnetic field-x
      BSy_32 = - sq3_2 * BStau_32 + 0.5 * BSneta_32   ! 1D resolved magnetic field-y


! ----------------------------------------------------------------------

#if( TEST_3STATE_CODE == 1)
      WRITE ( 6, '(A26, I1, A4)',advance="no")"BSx_13, BSy_13, ESz_13 ="
      WRITE ( 6, '(f10.4, f10.4, f10.4)')      BSx_13, BSy_13, ESz_13

      WRITE ( 6, '(A26, I1, A4)',advance="no")"BSx_21, BSy_21, ESz_21 ="
      WRITE ( 6, '(f10.4, f10.4, f10.4)')      BSx_21, BSy_21, ESz_21

      WRITE ( 6, '(A26, I1, A4)',advance="no")"BSx_32, BSy_32, ESz_32 ="
      WRITE ( 6, '(f10.4, f10.4, f10.4)')      BSx_32, BSy_32, ESz_32
#endif

! ------------------------------------------------------------

      IF ( ( S_L1 .GE. 0.0) .AND. ( S_R2 .LE. 0.0)) THEN

! ------------------------------------------------------------

        BSSx = Bx1
        BSSy = By1
        ESSz = Ez_1

        RETURN

! ------------------------------------------------------------

      ELSE IF ( ( S_L2 .GE. 0.0) .AND. ( S_R3 .LE. 0.0)) THEN

! ------------------------------------------------------------

        BSSx = Bx2
        BSSy = By2
        ESSz = Ez_2

        RETURN

! ------------------------------------------------------------

      ELSE IF ( ( S_L3 .GE. 0.0) .AND. ( S_R1 .LE. 0.0)) THEN

! ------------------------------------------------------------

        BSSx = Bx3
        BSSy = By3
        ESSz = Ez_3

        RETURN

! ------------------------------------------------------------

      ELSE IF ( ( S_L1 * ( S_R1 + 2.0 * S_L3)
     1          + S_R1 * ( S_L1 + 2.0 * S_R2)) .LE. 0) THEN

! ------------------------------------------------------------

        BSSx = BSx_13
        BSSy = BSy_13
        ESSz = ESz_13

        RETURN

! ------------------------------------------------------------

      ELSE IF ( ( ( S_L2 + S_R3) * ( S_L1 + 2.0 * S_R2)
     1          + ( S_L2 - S_R3) * S_L1) .LE. 0) THEN

! ------------------------------------------------------------

        BSSx = BSx_21
        BSSy = BSy_21
        ESSz = ESz_21

        RETURN

! ------------------------------------------------------------

      ELSE IF ( ( ( S_L2 + S_R3) * ( S_R1 + 2.0 * S_L3)
     1          - ( S_L2 - S_R3) * S_R1) .LE. 0) THEN

! ------------------------------------------------------------

        BSSx = BSx_32
        BSSy = BSy_32
        ESSz = ESz_32

        RETURN

! ------------------------------------------------------------

      ELSE

! ------------------------------------------------------------

        BSSx = ( ( X13_y * ESz_13 - S13 * BSx_13) * l1
     1         + ( X21_y * ESz_21 - S21 * BSx_21) * l2
     1         + ( X32_y * ESz_32 - S32 * BSx_32) * l3) / ( -2.0*Area)

        BSSy = ( ( - X13_x * ESz_13 - S13 * BSy_13) * l1
     1         + ( - X21_x * ESz_21 - S21 * BSy_21) * l2
     1         + ( - X32_x * ESz_32 - S32 * BSy_32) * l3) / ( -2.0*Area)


        ESSz_1 = - xc * BSSy
     1     - ( ( xb13 - xc) * ( - X13_x * ESz_13 - S13 * BSy_13) * l1
     1       + ( xb21 - xc) * ( - X21_x * ESz_21 - S21 * BSy_21) * l2
     1       + ( xb32 - xc) * ( - X32_x * ESz_32 - S32 * BSy_32) * l3)
     1     / Area

        ESSz_2 = yc * BSSx
     1     + ( ( yb13 - yc) * ( X13_y * ESz_13 - S13 * BSx_13) * l1
     1       + ( yb21 - yc) * ( X21_y * ESz_21 - S21 * BSx_21) * l2
     1       + ( yb32 - yc) * ( X32_y * ESz_32 - S32 * BSx_32) * l3)
     1     / Area


        ESSz = 0.5 * ( ESSz_1 + ESSz_2)

        RETURN

! ------------------------------------------------------------

      END IF

! ----------------------------------------------------------------------

      END SUBROUTINE Two_D_RS_Ez_3State_HLL





