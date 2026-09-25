











      SUBROUTINE TwoD_RS_4State ( n_cc_components,

     1           pi, rhofloor, prsfloor, smallnum, gamma,
     1           isospeed, barotropic_stiffness, gamma_euler,

     1           n1, n2, n3, B0,

     1           U_RU, U_LU, U_LD, U_RD,
     1           U_Dbl_Star, F_Dbl_Star, G_Dbl_Star, E_z )

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara and Sethu Subramanian

! Takes in two vectors "n1, n2" that don't have to be mutually orthogonal.
! Makes a primed space that is locally Cartesian with the z-axis perpendicular
! to the plane in which the 2DRS is solved.

! The four input statess "U_RU, U_LU, U_LD, U_RD" have their interior components
! rotated into that locally Cartesian space. Only then can we begin a solution
! of the 2D Riemann problem.

! "B0" is the mean, externally-specified, multipolar, magnetic field.

! "U_Dbl_Star, F_Dbl_Star, G_Dbl_Star" are the resolved state and fluxes in
! the x- and y-directions from the multidimensional Riemann solver.
! They are back-rotated to the Global space before return.

! "E_z" is the electric field along whatever is the z-axis of this problem.
! It is very useful for obtaining the electric field along the edges.

! ----------------------------------------------------------------------
!$acc routine seq

      IMPLICIT NONE

      INTEGER n_cc_components

      REAL pi, rhofloor, prsfloor, smallnum, gamma,
     1     isospeed, barotropic_stiffness, gamma_euler

      REAL n1 ( 3), n2 ( 3), n3 ( 3), B0 ( 3), E_z

      REAL, DIMENSION ( NVAR + NFLUID) :: U_RU, U_LU, U_LD, U_RD,
     1      U_Dbl_Star, F_Dbl_Star, G_Dbl_Star

! ----------------------------------------------------------------------

! Variables past this point are internal to the subroutine.

      REAL fourpi, ob_4pi, ob_8pi,
     1     a, b, c, phi, c_phi, s_phi, n12 ( 3),
     1     np1 ( 3), np2 ( 3), np3 ( 3),

     1     v2, b2, v_dot_b, b0_dot_b, P_tot, ob_rho,
     1     rho, prs, vx, vy, vz, bx, by, bz, vphi, bphi, eng, B0phi,

     1     cs2, vel_long, ob_4prho, va2_long, va2_tot, vf,

     1     S_L, S_R, S_U, S_D, S1, S2,

     1     x1, y1, x2, y2, x3, y3, x4, y4, xc, yc,
     1     ell12, ell23, ell34, ell41, x12, x23, x34, x41,
     1     y12, y23, y34, y41, Area

      REAL, DIMENSION ( NVAR + NFLUID) :: V_RU, V_LU, V_LD, V_RD,
     1     F_RU, F_LU, F_LD, F_RD, G_RU, G_LU, G_LD, G_RD,

     1     F_RU_tau2, F_LU_tau2, F_LD_tau2, F_RD_tau2,

     1     U1, U2, U3, U4, F1, F2, F3, F4,

     1     U_Star_R, V_Star_R, F_Star_R, G_Star_R,
     1     U_Star_L, V_Star_L, F_Star_L, G_Star_L,
     1     U_Star_U, V_Star_U, F_Star_U, G_Star_U,
     1     U_Star_D, V_Star_D, F_Star_D, G_Star_D

! ----------------------------------------------------------------------
! STEP 1) Obtain the constants.

      fourpi = 4.0 * pi
      ob_4pi = 1.0 / fourpi
      ob_8pi = 0.5 * ob_4pi

! ----------------------------------------------------------------------
! STEP 2) Construct unit vectors "n1 (3), n12 (3), n3 (3)". Also construct the
! unit vectors that give us the reverse transform via "np1 (3), np2 (3), np3 (3)".

! "n1" and "n2", as they come in, are in a global basis. Besides, they may
! not be orthogonal.

      a = 1.0 / SQRT ( n1 ( 1)**2 + n1 ( 2)**2 + n1 ( 3)**2 )
      n1 ( :) = n1 ( :) * a

      a = 1.0 / SQRT ( n2 ( 1)**2 + n2 ( 2)**2 + n2 ( 3)**2 )
      n2 ( :) = n2 ( :) * a

      c_phi = n1 ( 1) * n2 ( 1) + n1 ( 2) * n2 ( 2) + n1 ( 3) * n2 ( 3)
      phi = ACOS ( c_phi)
      s_phi = SQRT ( 1.0 - c_phi**2)

      n3 ( 1) = n1 ( 2) * n2 ( 3) - n2 ( 2) * n1 ( 3)
      n3 ( 2) = n1 ( 3) * n2 ( 1) - n2 ( 3) * n1 ( 1)
      n3 ( 3) = n1 ( 1) * n2 ( 2) - n2 ( 1) * n1 ( 2)

      a = 1.0 / SQRT ( n3 ( 1)**2 + n3 ( 2)**2 + n3 ( 3)**2 )
      n3 ( :) = n3 ( :) * a

! By this point, we have "n1, n2, n3" expressed in terms of a global basis.
! "n1, n2" may not be orthogonal. We construct "n12" which is perpendicular
! to "n1, n3" and in the direction of "n2".

      n12 ( :) = ( n2 ( :) - c_phi * n1 ( :) ) / s_phi

! By this point, "n1 (3), n12 (3), n3 (3)" are orthonormal unit vectors
! which give us a local "xprime, yprime, zprime" frame in which the 2D Riemann
! problem can be solved in the "xprime, yprime" plane. They are expressed in terms
! of a global basis set, so they can be used to transform from global "x, y, z"
! frame to a rotated "xprime, yprime, zprime" frame.

! We also wish to have the reverse transformation. The inverse of an orthonormal
! matrix is indeed its transverse. So "np1 (3), np2 (3), np3 (3)" give us that
! reverse transformation.

      np1 ( 1) = n1 ( 1) ; np1 ( 2) = n12 ( 1) ; np1 ( 3) = n3 ( 1)
      np2 ( 1) = n1 ( 2) ; np2 ( 2) = n12 ( 2) ; np2 ( 3) = n3 ( 2)
      np3 ( 1) = n1 ( 3) ; np3 ( 2) = n12 ( 3) ; np3 ( 3) = n3 ( 3)

! By this point "np1 (3), np2 (3), np3 (3)" give us the ability to go from
! the local "xprime, yprime, zprime" frame to the global "x, y, z" frame.

! ----------------------------------------------------------------------
! STEP 3) Now make internal in-place rotations to the "2,3,4" and "6,7,8"
! components of the four input states "U_RU, U_LU, U_LD, U_RD".
! This makes all the velocities and magnetic fields defined in the
! "xprime, yprime, zprime" frame.

! Also rotate "B0" into the "xprime, yprime, zprime" frame.

! ----------
! rotating background magnetic field "B0" to the local frame.
      a = B0 ( 1)
      b = B0 ( 2)
      c = B0 ( 3)

      B0 ( 1) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      B0 ( 2) = a * n12( 1) + b * n12( 2) + c * n12( 3)
      B0 ( 3) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! ----------
! rotating momentum components in "U_RU" to the local frame.
      a = U_RU ( 2)
      b = U_RU ( 3)
      c = U_RU ( 4)

      U_RU ( 2) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_RU ( 3) = a * n12( 1) + b * n12( 2) + c * n12( 3)
      U_RU ( 4) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! rotating magnetic field components in "U_RU" to the local frame.
      a = U_RU ( 6)
      b = U_RU ( 7)
      c = U_RU ( 8)

      U_RU ( 6) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_RU ( 7) = a * n12( 1) + b * n12( 2) + c * n12( 3)
      U_RU ( 8) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! ----------
! rotating momentum components in "U_LU" to the local frame.
      a = U_LU ( 2)
      b = U_LU ( 3)
      c = U_LU ( 4)

      U_LU ( 2) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_LU ( 3) = a * n12( 1) + b * n12( 2) + c * n12( 3)
      U_LU ( 4) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! rotating magnetic field components in "U_LU" to the local frame.
      a = U_LU ( 6)
      b = U_LU ( 7)
      c = U_LU ( 8)

      U_LU ( 6) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_LU ( 7) = a * n12( 1) + b * n12( 2) + c * n12( 3)
      U_LU ( 8) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! ----------
! rotating momentum components in "U_LD" to the local frame.
      a = U_LD ( 2)
      b = U_LD ( 3)
      c = U_LD ( 4)

      U_LD ( 2) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_LD ( 3) = a * n12( 1) + b * n12( 2) + c * n12( 3)
      U_LD ( 4) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! rotating magnetic field components in "U_LD" to the local frame.
      a = U_LD ( 6)
      b = U_LD ( 7)
      c = U_LD ( 8)

      U_LD ( 6) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_LD ( 7) = a * n12( 1) + b * n12( 2) + c * n12( 3)
      U_LD ( 8) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! ----------
! rotating momentum components in "U_RD" to the local frame.
      a = U_RD ( 2)
      b = U_RD ( 3)
      c = U_RD ( 4)

      U_RD ( 2) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_RD ( 3) = a * n12( 1) + b * n12( 2) + c * n12( 3)
      U_RD ( 4) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! rotating magnetic field components in "U_RD" to the local frame.
      a = U_RD ( 6)
      b = U_RD ( 7)
      c = U_RD ( 8)

      U_RD ( 6) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_RD ( 7) = a * n12( 1) + b * n12( 2) + c * n12( 3)
      U_RD ( 8) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! ----------------------------------------------------------------------
! STEP 4) From "U_RU" etc. build "V_RU, F_RU, G_RU". We also build
! "F_RU_tau2 ( :) = s_phi * F_RU ( :) - c_phi * G_RU ( :)", so that we can apply the
! 1D RS perpendicular to the "n2" direction.

! The fluxes are defined in the "xprime, yprime, zprime" frame.

! ----------------------------------------
! Building "V_RU".

      U_RU ( 1) = AMAX1 ( U_RU ( 1), rhofloor)

      ob_rho  = U_RU ( 1) / ( U_RU ( 1)**2 + rhofloor**2)

      V_RU ( 1)    = U_RU ( 1)               ! density
      V_RU ( 2: 4) = U_RU ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_RU ( 6: 8) = U_RU ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_RU ( NVAR + 1: ) = U_RU ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_RU ( 2)**2 + V_RU ( 3)**2 + V_RU ( 4)**2
      B2 = V_RU ( 6)**2 + V_RU ( 7)**2 + V_RU ( 8)**2

      V_dot_B  = V_RU ( 2) * V_RU ( 6)
     1         + V_RU ( 3) * V_RU ( 7)
     1         + V_RU ( 4) * V_RU ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_RU ( 5) = ( gamma - 1.0) * ( U_RU ( 5)
     1             - 0.5 * V_RU ( 1) * v2  - ob_8pi * B2  )

      V_RU ( 5) = AMAX1 ( V_RU ( 5), prsfloor)  ! gas pressure

      U_RU ( 5) = V_RU ( 5) / ( gamma - 1.0)
     1          + 0.5 * V_RU ( 1) * v2 + ob_8pi * B2
#else
      V_RU ( 5) = isospeed**2 * U_RU ( 1)
     1          + barotropic_stiffness * U_RU ( 1)**gamma_euler
#endif

! --------------------
! Building "F_RU, G_RU".

      B0_dot_B  = B0 ( 1) * V_RU ( 6)
     1          + B0 ( 2) * V_RU ( 7)
     1          + B0 ( 3) * V_RU ( 8)

      P_tot  = V_RU ( 5) + ob_8pi * B2
     1                   + ob_4pi * B0_dot_B

! ----------
! Building "F_RU".

      F_RU ( :) = V_RU ( 2) * U_RU ( :)


      F_RU ( 2) = F_RU ( 2) + P_tot
     1          - ob_4pi * V_RU ( 6) * ( V_RU ( 6) + B0 ( 1))
     1          - ob_4pi * V_RU ( 6) * B0 ( 1)

      F_RU ( 3) = F_RU ( 3)
     1          - ob_4pi * V_RU ( 6) * V_RU ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_RU ( 7)
     1                     + B0 ( 2) * V_RU ( 6))

      F_RU ( 4) = F_RU ( 4)
     1          - ob_4pi * V_RU ( 6) * V_RU ( 8)
     1          - ob_4pi * ( B0 ( 1) * V_RU ( 8)
     1                     + B0 ( 3) * V_RU ( 6))

#if ( CHOOSE_ADIABATIC == 1)
      F_RU ( 5) = F_RU ( 5)
     1          + P_tot  * V_RU ( 2)
     1          - ob_4pi * ( V_RU ( 6) + B0 ( 1))* V_dot_B
#else
      F_RU ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      F_RU ( 6) = 0.0

      F_RU ( 7) = F_RU ( 7) - V_RU ( 6) * V_RU ( 3)
     1          + B0 ( 2) * V_RU ( 2) - B0 ( 1) * V_RU ( 3)

      F_RU ( 8) = F_RU ( 8) - V_RU ( 6) * V_RU ( 4)
     1          + B0 ( 3) * V_RU ( 2) - B0 ( 1) * V_RU ( 4)

! ----------
! Building "G_RU".

      G_RU ( :) = V_RU ( 3) * U_RU ( :)

      G_RU ( 2) = G_RU ( 2)
     1          - ob_4pi * V_RU ( 6) * V_RU ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_RU ( 7)
     1                     + B0 ( 2) * V_RU ( 6))

      G_RU ( 3) = G_RU ( 3) + P_tot
     1           - ob_4pi * V_RU ( 7) * ( V_RU ( 7) + B0 ( 2))
     1           - ob_4pi * V_RU ( 7) * B0 ( 2)

      G_RU ( 4) = G_RU ( 4)
     1          - ob_4pi * V_RU ( 7) * V_RU ( 8)
     1          - ob_4pi * ( B0 ( 2) * V_RU ( 8)
     1                     + B0 ( 3) * V_RU ( 7))

#if ( CHOOSE_ADIABATIC == 1)
      G_RU ( 5) = G_RU ( 5)
     1          + P_tot  * V_RU ( 3)
     1          - ob_4pi * ( V_RU ( 7) + B0 ( 2))* V_dot_B
#else
      G_RU ( 5) = 0.0             ! isothermal has no energy flux.
#endif

      G_RU ( 6) = G_RU ( 6) - V_RU ( 2) * V_RU ( 7)
     1          + V_RU ( 3) * B0 ( 1) - V_RU ( 2) * B0 ( 2)

      G_RU ( 7) = 0.0

      G_RU ( 8) = G_RU ( 8) - V_RU ( 4) * V_RU ( 7)
     1          + V_RU ( 3) * B0 ( 3) - V_RU ( 4) * B0 ( 2)

! ----------

      F_RU_tau2 ( :) = s_phi * F_RU ( :) - c_phi * G_RU ( :)

! ----------------------------------------
! Building "V_LU".

      U_LU ( 1) = AMAX1 ( U_LU ( 1), rhofloor)

      ob_rho  = U_LU ( 1) / ( U_LU ( 1)**2 + rhofloor**2)

      V_LU ( 1)    = U_LU ( 1)               ! density
      V_LU ( 2: 4) = U_LU ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_LU ( 6: 8) = U_LU ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_LU ( NVAR + 1: ) = U_LU ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_LU ( 2)**2 + V_LU ( 3)**2 + V_LU ( 4)**2
      B2 = V_LU ( 6)**2 + V_LU ( 7)**2 + V_LU ( 8)**2

      V_dot_B  = V_LU ( 2) * V_LU ( 6)
     1         + V_LU ( 3) * V_LU ( 7)
     1         + V_LU ( 4) * V_LU ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_LU ( 5) = ( gamma - 1.0) * ( U_LU ( 5)
     1             - 0.5 * V_LU ( 1) * v2  - ob_8pi * B2  )

      V_LU ( 5) = AMAX1 ( V_LU ( 5), prsfloor)  ! gas pressure

      U_LU ( 5) = V_LU ( 5) / ( gamma - 1.0)
     1          + 0.5 * V_LU ( 1) * v2 + ob_8pi * B2
#else
      V_LU ( 5) = isospeed**2 * U_LU ( 1)
     1          + barotropic_stiffness * U_LU ( 1)**gamma_euler
#endif

! --------------------
! Building "F_LU, G_LU".

      B0_dot_B  = B0 ( 1) * V_LU ( 6)
     1          + B0 ( 2) * V_LU ( 7)
     1          + B0 ( 3) * V_LU ( 8)

      P_tot  = V_LU ( 5) + ob_8pi * B2
     1                   + ob_4pi * B0_dot_B

! ----------
! Building "F_LU".

      F_LU ( :) = V_LU ( 2) * U_LU ( :)


      F_LU ( 2) = F_LU ( 2) + P_tot
     1          - ob_4pi * V_LU ( 6) * ( V_LU ( 6) + B0 ( 1))
     1          - ob_4pi * V_LU ( 6) * B0 ( 1)

      F_LU ( 3) = F_LU ( 3)
     1          - ob_4pi * V_LU ( 6) * V_LU ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_LU ( 7)
     1                     + B0 ( 2) * V_LU ( 6))

      F_LU ( 4) = F_LU ( 4)
     1          - ob_4pi * V_LU ( 6) * V_LU ( 8)
     1          - ob_4pi * ( B0 ( 1) * V_LU ( 8)
     1                     + B0 ( 3) * V_LU ( 6))

#if ( CHOOSE_ADIABATIC == 1)
      F_LU ( 5) = F_LU ( 5)
     1          + P_tot  * V_LU ( 2)
     1          - ob_4pi * ( V_LU ( 6) + B0 ( 1))* V_dot_B
#else
      F_LU ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      F_LU ( 6) = 0.0

      F_LU ( 7) = F_LU ( 7) - V_LU ( 6) * V_LU ( 3)
     1          + B0 ( 2) * V_LU ( 2) - B0 ( 1) * V_LU ( 3)

      F_LU ( 8) = F_LU ( 8) - V_LU ( 6) * V_LU ( 4)
     1          + B0 ( 3) * V_LU ( 2) - B0 ( 1) * V_LU ( 4)

! ----------
! Building "G_LU".

      G_LU ( :) = V_LU ( 3) * U_LU ( :)

      G_LU ( 2) = G_LU ( 2)
     1          - ob_4pi * V_LU ( 6) * V_LU ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_LU ( 7)
     1                     + B0 ( 2) * V_LU ( 6))

      G_LU ( 3) = G_LU ( 3) + P_tot
     1           - ob_4pi * V_LU ( 7) * ( V_LU ( 7) + B0 ( 2))
     1           - ob_4pi * V_LU ( 7) * B0 ( 2)

      G_LU ( 4) = G_LU ( 4)
     1          - ob_4pi * V_LU ( 7) * V_LU ( 8)
     1          - ob_4pi * ( B0 ( 2) * V_LU ( 8)
     1                     + B0 ( 3) * V_LU ( 7))

#if ( CHOOSE_ADIABATIC == 1)
      G_LU ( 5) = G_LU ( 5)
     1          + P_tot  * V_LU ( 3)
     1          - ob_4pi * ( V_LU ( 7) + B0 ( 2))* V_dot_B
#else
      G_LU ( 5) = 0.0             ! isothermal has no energy flux.
#endif

      G_LU ( 6) = G_LU ( 6) - V_LU ( 2) * V_LU ( 7)
     1          + V_LU ( 3) * B0 ( 1) - V_LU ( 2) * B0 ( 2)

      G_LU ( 7) = 0.0

      G_LU ( 8) = G_LU ( 8) - V_LU ( 4) * V_LU ( 7)
     1          + V_LU ( 3) * B0 ( 3) - V_LU ( 4) * B0 ( 2)

! ----------

      F_LU_tau2 ( :) = s_phi * F_LU ( :) - c_phi * G_LU ( :)

! ----------------------------------------
! Building "V_LD".

      U_LD ( 1) = AMAX1 ( U_LD ( 1), rhofloor)

      ob_rho  = U_LD ( 1) / ( U_LD ( 1)**2 + rhofloor**2)

      V_LD ( 1)    = U_LD ( 1)               ! density
      V_LD ( 2: 4) = U_LD ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_LD ( 6: 8) = U_LD ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_LD ( NVAR + 1: ) = U_LD ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_LD ( 2)**2 + V_LD ( 3)**2 + V_LD ( 4)**2
      B2 = V_LD ( 6)**2 + V_LD ( 7)**2 + V_LD ( 8)**2

      V_dot_B  = V_LD ( 2) * V_LD ( 6)
     1         + V_LD ( 3) * V_LD ( 7)
     1         + V_LD ( 4) * V_LD ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_LD ( 5) = ( gamma - 1.0) * ( U_LD ( 5)
     1             - 0.5 * V_LD ( 1) * v2  - ob_8pi * B2  )

      V_LD ( 5) = AMAX1 ( V_LD ( 5), prsfloor)  ! gas pressure

      U_LD ( 5) = V_LD ( 5) / ( gamma - 1.0)
     1          + 0.5 * V_LD ( 1) * v2 + ob_8pi * B2
#else
      V_LD ( 5) = isospeed**2 * U_LD ( 1)
     1          + barotropic_stiffness * U_LD ( 1)**gamma_euler
#endif

! --------------------
! Building "F_LD, G_LD".

      B0_dot_B  = B0 ( 1) * V_LD ( 6)
     1          + B0 ( 2) * V_LD ( 7)
     1          + B0 ( 3) * V_LD ( 8)

      P_tot  = V_LD ( 5) + ob_8pi * B2
     1                   + ob_4pi * B0_dot_B

! ----------
! Building "F_LD".

      F_LD ( :) = V_LD ( 2) * U_LD ( :)


      F_LD ( 2) = F_LD ( 2) + P_tot
     1          - ob_4pi * V_LD ( 6) * ( V_LD ( 6) + B0 ( 1))
     1          - ob_4pi * V_LD ( 6) * B0 ( 1)

      F_LD ( 3) = F_LD ( 3)
     1          - ob_4pi * V_LD ( 6) * V_LD ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_LD ( 7)
     1                     + B0 ( 2) * V_LD ( 6))

      F_LD ( 4) = F_LD ( 4)
     1          - ob_4pi * V_LD ( 6) * V_LD ( 8)
     1          - ob_4pi * ( B0 ( 1) * V_LD ( 8)
     1                     + B0 ( 3) * V_LD ( 6))

#if ( CHOOSE_ADIABATIC == 1)
      F_LD ( 5) = F_LD ( 5)
     1          + P_tot  * V_LD ( 2)
     1          - ob_4pi * ( V_LD ( 6) + B0 ( 1))* V_dot_B
#else
      F_LD ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      F_LD ( 6) = 0.0

      F_LD ( 7) = F_LD ( 7) - V_LD ( 6) * V_LD ( 3)
     1          + B0 ( 2) * V_LD ( 2) - B0 ( 1) * V_LD ( 3)

      F_LD ( 8) = F_LD ( 8) - V_LD ( 6) * V_LD ( 4)
     1          + B0 ( 3) * V_LD ( 2) - B0 ( 1) * V_LD ( 4)

! ----------
! Building "G_LD".

      G_LD ( :) = V_LD ( 3) * U_LD ( :)

      G_LD ( 2) = G_LD ( 2)
     1          - ob_4pi * V_LD ( 6) * V_LD ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_LD ( 7)
     1                     + B0 ( 2) * V_LD ( 6))

      G_LD ( 3) = G_LD ( 3) + P_tot
     1           - ob_4pi * V_LD ( 7) * ( V_LD ( 7) + B0 ( 2))
     1           - ob_4pi * V_LD ( 7) * B0 ( 2)

      G_LD ( 4) = G_LD ( 4)
     1          - ob_4pi * V_LD ( 7) * V_LD ( 8)
     1          - ob_4pi * ( B0 ( 2) * V_LD ( 8)
     1                     + B0 ( 3) * V_LD ( 7))

#if ( CHOOSE_ADIABATIC == 1)
      G_LD ( 5) = G_LD ( 5)
     1          + P_tot  * V_LD ( 3)
     1          - ob_4pi * ( V_LD ( 7) + B0 ( 2))* V_dot_B
#else
      G_LD ( 5) = 0.0             ! isothermal has no energy flux.
#endif

      G_LD ( 6) = G_LD ( 6) - V_LD ( 2) * V_LD ( 7)
     1          + V_LD ( 3) * B0 ( 1) - V_LD ( 2) * B0 ( 2)

      G_LD ( 7) = 0.0

      G_LD ( 8) = G_LD ( 8) - V_LD ( 4) * V_LD ( 7)
     1          + V_LD ( 3) * B0 ( 3) - V_LD ( 4) * B0 ( 2)

! ----------

      F_LD_tau2 ( :) = s_phi * F_LD ( :) - c_phi * G_LD ( :)

! ----------------------------------------
! Building "V_RD".

      U_RD ( 1) = AMAX1 ( U_RD ( 1), rhofloor)

      ob_rho  = U_RD ( 1) / ( U_RD ( 1)**2 + rhofloor**2)

      V_RD ( 1)    = U_RD ( 1)               ! density
      V_RD ( 2: 4) = U_RD ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_RD ( 6: 8) = U_RD ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_RD ( NVAR + 1: ) = U_RD ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_RD ( 2)**2 + V_RD ( 3)**2 + V_RD ( 4)**2
      B2 = V_RD ( 6)**2 + V_RD ( 7)**2 + V_RD ( 8)**2

      V_dot_B  = V_RD ( 2) * V_RD ( 6)
     1         + V_RD ( 3) * V_RD ( 7)
     1         + V_RD ( 4) * V_RD ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_RD ( 5) = ( gamma - 1.0) * ( U_RD ( 5)
     1             - 0.5 * V_RD ( 1) * v2  - ob_8pi * B2  )

      V_RD ( 5) = AMAX1 ( V_RD ( 5), prsfloor)  ! gas pressure

      U_RD ( 5) = V_RD ( 5) / ( gamma - 1.0)
     1          + 0.5 * V_RD ( 1) * v2 + ob_8pi * B2
#else
      V_RD ( 5) = isospeed**2 * U_RD ( 1)
     1          + barotropic_stiffness * U_RD ( 1)**gamma_euler
#endif

! --------------------
! Building "F_RD, G_RD".

      B0_dot_B  = B0 ( 1) * V_RD ( 6)
     1          + B0 ( 2) * V_RD ( 7)
     1          + B0 ( 3) * V_RD ( 8)

      P_tot  = V_RD ( 5) + ob_8pi * B2
     1                   + ob_4pi * B0_dot_B

! ----------
! Building "F_RD".

      F_RD ( :) = V_RD ( 2) * U_RD ( :)


      F_RD ( 2) = F_RD ( 2) + P_tot
     1          - ob_4pi * V_RD ( 6) * ( V_RD ( 6) + B0 ( 1))
     1          - ob_4pi * V_RD ( 6) * B0 ( 1)

      F_RD ( 3) = F_RD ( 3)
     1          - ob_4pi * V_RD ( 6) * V_RD ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_RD ( 7)
     1                     + B0 ( 2) * V_RD ( 6))

      F_RD ( 4) = F_RD ( 4)
     1          - ob_4pi * V_RD ( 6) * V_RD ( 8)
     1          - ob_4pi * ( B0 ( 1) * V_RD ( 8)
     1                     + B0 ( 3) * V_RD ( 6))

#if ( CHOOSE_ADIABATIC == 1)
      F_RD ( 5) = F_RD ( 5)
     1          + P_tot  * V_RD ( 2)
     1          - ob_4pi * ( V_RD ( 6) + B0 ( 1))* V_dot_B
#else
      F_RD ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      F_RD ( 6) = 0.0

      F_RD ( 7) = F_RD ( 7) - V_RD ( 6) * V_RD ( 3)
     1          + B0 ( 2) * V_RD ( 2) - B0 ( 1) * V_RD ( 3)

      F_RD ( 8) = F_RD ( 8) - V_RD ( 6) * V_RD ( 4)
     1          + B0 ( 3) * V_RD ( 2) - B0 ( 1) * V_RD ( 4)

! ----------
! Building "G_RD".

      G_RD ( :) = V_RD ( 3) * U_RD ( :)

      G_RD ( 2) = G_RD ( 2)
     1          - ob_4pi * V_RD ( 6) * V_RD ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_RD ( 7)
     1                     + B0 ( 2) * V_RD ( 6))

      G_RD ( 3) = G_RD ( 3) + P_tot
     1           - ob_4pi * V_RD ( 7) * ( V_RD ( 7) + B0 ( 2))
     1           - ob_4pi * V_RD ( 7) * B0 ( 2)

      G_RD ( 4) = G_RD ( 4)
     1          - ob_4pi * V_RD ( 7) * V_RD ( 8)
     1          - ob_4pi * ( B0 ( 2) * V_RD ( 8)
     1                     + B0 ( 3) * V_RD ( 7))

#if ( CHOOSE_ADIABATIC == 1)
      G_RD ( 5) = G_RD ( 5)
     1          + P_tot  * V_RD ( 3)
     1          - ob_4pi * ( V_RD ( 7) + B0 ( 2))* V_dot_B
#else
      G_RD ( 5) = 0.0             ! isothermal has no energy flux.
#endif

      G_RD ( 6) = G_RD ( 6) - V_RD ( 2) * V_RD ( 7)
     1          + V_RD ( 3) * B0 ( 1) - V_RD ( 2) * B0 ( 2)

      G_RD ( 7) = 0.0

      G_RD ( 8) = G_RD ( 8) - V_RD ( 4) * V_RD ( 7)
     1          + V_RD ( 3) * B0 ( 3) - V_RD ( 4) * B0 ( 2)

! ----------

      F_RD_tau2 ( :) = s_phi * F_RD ( :) - c_phi * G_RD ( :)

! ----------------------------------------------------------------------
! STEP 5) Find "S_U, S_D" which are the extremal speeds perpendicular to the
! "n1" axis. These are useful for the Riemann fan that is perpendicular to "n1".
! Also find "S_R, S_L" which are the extremal speeds perpendicular to the
! "n2" axis. These are useful for the Riemann fan that is perpendicular to "n2".

! ------------------------------

      S_R = -1.0e30
      S_L = 1.0e30

      S_U = -1.0e30
      S_D = 1.0e30

! ------------------------------
! Extremize "S_U, S_D" between "U_RU" and "U_RD" in three steps.

! ----------
! Extremize "S_U, S_D" using "U_RU".

      rho = V_RU ( 1)

      vx = V_RU ( 2)
      vy = V_RU ( 3)
      vz = V_RU ( 4)

      prs = V_RU ( 5)

      bx = V_RU ( 6) + B0 ( 1)
      by = V_RU ( 7) + B0 ( 2)
      bz = V_RU ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = by**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_U = AMAX1 ( S_U, vel_long + vf )
      S_D = AMIN1 ( S_D, vel_long - vf )

! ----------
! Extremize "S_U, S_D" using "U_RD".

      rho = V_RD ( 1)

      vx = V_RD ( 2)
      vy = V_RD ( 3)
      vz = V_RD ( 4)

      prs = V_RD ( 5)

      bx = V_RD ( 6) + B0 ( 1)
      by = V_RD ( 7) + B0 ( 2)
      bz = V_RD ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = by**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_U = AMAX1 ( S_U, vel_long + vf )
      S_D = AMIN1 ( S_D, vel_long - vf )

! ----------
! Extremize "S_U, S_D" using "0.5 * ( U_RU + U_RD)"

      rho = 0.5 * ( V_RU ( 1) + V_RD ( 1))

      vx = 0.5 * ( V_RU ( 2) + V_RD ( 2))
      vy = 0.5 * ( V_RU ( 3) + V_RD ( 3))
      vz = 0.5 * ( V_RU ( 4) + V_RD ( 4))

      prs = 0.5 * ( V_RU ( 5) + V_RD ( 5))

      bx = 0.5 * ( V_RU ( 6) + V_RD ( 6)) + B0 ( 1)
      by = 0.5 * ( V_RU ( 7) + V_RD ( 7)) + B0 ( 2)
      bz = 0.5 * ( V_RU ( 8) + V_RD ( 8)) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = by**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_U = AMAX1 ( S_U, vel_long + vf )
      S_D = AMIN1 ( S_D, vel_long - vf )

! ------------------------------
! Extremize "S_U, S_D" between "U_LU" and "U_LD" in three steps.

! ----------
! Extremize "S_U, S_D" using "U_LU".

      rho = V_LU ( 1)

      vx = V_LU ( 2)
      vy = V_LU ( 3)
      vz = V_LU ( 4)

      prs = V_LU ( 5)

      bx = V_LU ( 6) + B0 ( 1)
      by = V_LU ( 7) + B0 ( 2)
      bz = V_LU ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = by**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_U = AMAX1 ( S_U, vel_long + vf )
      S_D = AMIN1 ( S_D, vel_long - vf )

! ----------
! Extremize "S_U, S_D" using "U_LD".

      rho = V_LD ( 1)

      vx = V_LD ( 2)
      vy = V_LD ( 3)
      vz = V_LD ( 4)

      prs = V_LD ( 5)

      bx = V_LD ( 6) + B0 ( 1)
      by = V_LD ( 7) + B0 ( 2)
      bz = V_LD ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = by**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_U = AMAX1 ( S_U, vel_long + vf )
      S_D = AMIN1 ( S_D, vel_long - vf )

! ----------
! Extremize "S_U, S_D" using "0.5 * ( U_LU + U_LD)"

      rho = 0.5 * ( V_LU ( 1) + V_LD ( 1))

      vx = 0.5 * ( V_LU ( 2) + V_LD ( 2))
      vy = 0.5 * ( V_LU ( 3) + V_LD ( 3))
      vz = 0.5 * ( V_LU ( 4) + V_LD ( 4))

      prs = 0.5 * ( V_LU ( 5) + V_LD ( 5))

      bx = 0.5 * ( V_LU ( 6) + V_LD ( 6)) + B0 ( 1)
      by = 0.5 * ( V_LU ( 7) + V_LD ( 7)) + B0 ( 2)
      bz = 0.5 * ( V_LU ( 8) + V_LD ( 8)) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = by**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_U = AMAX1 ( S_U, vel_long + vf )
      S_D = AMIN1 ( S_D, vel_long - vf )

! ------------------------------
! Extremize "S_R, S_L" between "U_RU" and "U_LU" in three steps.

! ----------
! Extremize "S_R, S_L" using "U_RU".

      rho = V_RU ( 1)

      vx = V_RU ( 2)
      vy = V_RU ( 3)
      vz = V_RU ( 4)

      prs = V_RU ( 5)

      bx = V_RU ( 6) + B0 ( 1)
      by = V_RU ( 7) + B0 ( 2)
      bz = V_RU ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = s_phi * vx - c_phi * vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( s_phi * bx - c_phi * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R = AMAX1 ( S_R, vel_long + vf )
      S_L = AMIN1 ( S_L, vel_long - vf )

! ----------
! Extremize "S_R, S_L" using "U_LU".

      rho = V_LU ( 1)

      vx = V_LU ( 2)
      vy = V_LU ( 3)
      vz = V_LU ( 4)

      prs = V_LU ( 5)

      bx = V_LU ( 6) + B0 ( 1)
      by = V_LU ( 7) + B0 ( 2)
      bz = V_LU ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = s_phi * vx - c_phi * vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( s_phi * bx - c_phi * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R = AMAX1 ( S_R, vel_long + vf )
      S_L = AMIN1 ( S_L, vel_long - vf )

! ----------
! Extremize "S_R, S_L" using "0.5 * ( U_RU + U_LU)".

      rho = 0.5 * ( V_RU ( 1) + V_LU ( 1))

      vx = 0.5 * ( V_RU ( 2) + V_LU ( 2))
      vy = 0.5 * ( V_RU ( 3) + V_LU ( 3))
      vz = 0.5 * ( V_RU ( 4) + V_LU ( 4))

      prs = 0.5 * ( V_RU ( 5) + V_LU ( 5))

      bx = 0.5 * ( V_RU ( 6) + V_LU ( 6)) + B0 ( 1)
      by = 0.5 * ( V_RU ( 7) + V_LU ( 7)) + B0 ( 2)
      bz = 0.5 * ( V_RU ( 8) + V_LU ( 8)) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = s_phi * vx - c_phi * vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( s_phi * bx - c_phi * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R = AMAX1 ( S_R, vel_long + vf )
      S_L = AMIN1 ( S_L, vel_long - vf )

! ------------------------------
! Extremize "S_R, S_L" between "U_RD" and "U_LD" in three steps.

! ----------
! Extremize "S_R, S_L" using "U_RD".

      rho = V_RD ( 1)

      vx = V_RD ( 2)
      vy = V_RD ( 3)
      vz = V_RD ( 4)

      prs = V_RD ( 5)

      bx = V_RD ( 6) + B0 ( 1)
      by = V_RD ( 7) + B0 ( 2)
      bz = V_RD ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = s_phi * vx - c_phi * vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( s_phi * bx - c_phi * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R = AMAX1 ( S_R, vel_long + vf )
      S_L = AMIN1 ( S_L, vel_long - vf )

! ----------
! Extremize "S_R, S_L" using "U_LD".

      rho = V_LD ( 1)

      vx = V_LD ( 2)
      vy = V_LD ( 3)
      vz = V_LD ( 4)

      prs = V_LD ( 5)

      bx = V_LD ( 6) + B0 ( 1)
      by = V_LD ( 7) + B0 ( 2)
      bz = V_LD ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = s_phi * vx - c_phi * vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( s_phi * bx - c_phi * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R = AMAX1 ( S_R, vel_long + vf )
      S_L = AMIN1 ( S_L, vel_long - vf )

! ----------
! Extremize "S_R, S_L" using "0.5 * ( U_RD + U_LD)".

      rho = 0.5 * ( V_RD ( 1) + V_LD ( 1))

      vx = 0.5 * ( V_RD ( 2) + V_LD ( 2))
      vy = 0.5 * ( V_RD ( 3) + V_LD ( 3))
      vz = 0.5 * ( V_RD ( 4) + V_LD ( 4))

      prs = 0.5 * ( V_RD ( 5) + V_LD ( 5))

      bx = 0.5 * ( V_RD ( 6) + V_LD ( 6)) + B0 ( 1)
      by = 0.5 * ( V_RD ( 7) + V_LD ( 7)) + B0 ( 2)
      bz = 0.5 * ( V_RD ( 8) + V_LD ( 8)) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = s_phi * vx - c_phi * vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( s_phi * bx - c_phi * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R = AMAX1 ( S_R, vel_long + vf )
      S_L = AMIN1 ( S_L, vel_long - vf )

! ----------------------------------------------------------------------
! STEP 6) Obtain the resolved states and longitudinal fluxes for the four
! 1D Riemann problems that bound the strongly-interacting state.

! ------------------------------------------------------------
! 1D RP between "U_RU and U_RD" in the +ve y-direction.

      S1 = S_U  ;  S2 = S_D

      U1 = U_RU  ;  U2 = U_RD  ;  F1 = G_RU  ;  F2 = G_RD

      c = 1.0 / ( S1 - S2)
      a = S1 * c
      b = S2 * c

      U_Star_R = a * U1 - b * U2 - c * ( F1 - F2)

      c = S1 * S2 * c

      G_Star_R = a * F2 - b * F1 + c * ( U1 - U2)

! ------------------------------------------------------------
! 1D RP between "U_LU and U_LD" in the +ve y-direction.

      S1 = S_U  ;  S2 = S_D

      U1 = U_LU  ;  U2 = U_LD  ;  F1 = G_LU  ;  F2 = G_LD

      c = 1.0 / ( S1 - S2)
      a = S1 * c
      b = S2 * c

      U_Star_L = a * U1 - b * U2 - c * ( F1 - F2)

      c = S1 * S2 * c

      G_Star_L = a * F2 - b * F1 + c * ( U1 - U2)

! ------------------------------------------------------------
! 1D RP between "U_RU and U_LU" in the tau2-direction.

      S1 = S_R  ;  S2 = S_L

      U1 = U_RU  ;  U2 = U_LU  ;  F1 = F_RU_tau2  ; F2 = F_LU_tau2

      c = 1.0 / ( S1 - S2)
      a = S1 * c
      b = S2 * c

      U_Star_U = a * U1 - b * U2 - c * ( F1 - F2)

      c = S1 * S2 * c

      F_Star_U = a * F2 - b * F1 + c * ( U1 - U2)  ! Watch out; this is in tau2-direction.

! ------------------------------------------------------------
! 1D RP between "U_RD and U_LD" in the tau2-direction.

      S1 = S_R  ;  S2 = S_L

      U1 = U_RD  ;  U2 = U_LD  ;  F1 = F_RD_tau2  ;  F2 = F_LD_tau2

      c = 1.0 / ( S1 - S2)
      a = S1 * c
      b = S2 * c

      U_Star_D = a * U1 - b * U2 - c * ( F1 - F2)

      c = S1 * S2 * c

      F_Star_D = a * F2 - b * F1 + c * ( U1 - U2)  ! Watch out; this is in tau2-direction.

! ----------------------------------------------------------------------
! STEP 7) The RS has given us four resolved states and four longitudinal fluxes.
! We still need to construct the primitive variables from those states and use
! them to create the transverse fluxes.

! ------------------------------------------------------------
! "U_Star_R --> V_Star_R --> F_Star_R" gives us the flux in the +ve x-direction.

! --------------------
! Building "V_Star_R".

      U_Star_R ( 1) = AMAX1 ( U_Star_R ( 1), rhofloor)

      ob_rho  = U_Star_R ( 1) / ( U_Star_R ( 1)**2 + rhofloor**2)

      V_Star_R ( 1)    = U_Star_R ( 1)               ! density
      V_Star_R ( 2: 4) = U_Star_R ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_Star_R ( 6: 8) = U_Star_R ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_Star_R ( NVAR + 1: ) = U_Star_R ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_Star_R ( 2)**2 + V_Star_R ( 3)**2 + V_Star_R ( 4)**2
      B2 = V_Star_R ( 6)**2 + V_Star_R ( 7)**2 + V_Star_R ( 8)**2

      V_dot_B  = V_Star_R ( 2) * V_Star_R ( 6)
     1         + V_Star_R ( 3) * V_Star_R ( 7)
     1         + V_Star_R ( 4) * V_Star_R ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_Star_R ( 5) = ( gamma - 1.0) * ( U_Star_R ( 5)
     1             - 0.5 * V_Star_R ( 1) * v2  - ob_8pi * B2  )

      V_Star_R ( 5) = AMAX1 ( V_Star_R ( 5), prsfloor)  ! gas pressure

      U_Star_R ( 5) = V_Star_R ( 5) / ( gamma - 1.0)
     1        + 0.5 * V_Star_R ( 1) * v2 + ob_8pi * B2
#else
      V_Star_R ( 5) = isospeed**2 * U_Star_R ( 1)
     1          + barotropic_stiffness * U_Star_R ( 1)**gamma_euler
#endif


! --------------------
! Building "F_Star_R".

      B0_dot_B  = B0 ( 1) * V_Star_R ( 6)
     1          + B0 ( 2) * V_Star_R ( 7)
     1          + B0 ( 3) * V_Star_R ( 8)

      P_tot  = V_Star_R ( 5) + ob_8pi * B2
     1                       + ob_4pi * B0_dot_B

! ----------

      F_Star_R ( :) = V_Star_R ( 2) * U_Star_R ( :)


      F_Star_R ( 2) = F_Star_R ( 2) + P_tot
     1          - ob_4pi * V_Star_R ( 6) * ( V_Star_R ( 6) + B0 ( 1))
     1          - ob_4pi * V_Star_R ( 6) * B0 ( 1)

      F_Star_R ( 3) = F_Star_R ( 3)
     1          - ob_4pi * V_Star_R ( 6) * V_Star_R ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_Star_R ( 7)
     1                     + B0 ( 2) * V_Star_R ( 6))

      F_Star_R ( 4) = F_Star_R ( 4)
     1          - ob_4pi * V_Star_R ( 6) * V_Star_R ( 8)
     1          - ob_4pi * ( B0 ( 1) * V_Star_R ( 8)
     1                     + B0 ( 3) * V_Star_R ( 6))

#if ( CHOOSE_ADIABATIC == 1)
      F_Star_R ( 5) = F_Star_R ( 5)
     1          + P_tot  * V_Star_R ( 2)
     1          - ob_4pi * ( V_Star_R ( 6) + B0 ( 1))* V_dot_B
#else
      F_Star_R ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      F_Star_R ( 6) = 0.0

      F_Star_R ( 7) = F_Star_R ( 7) - V_Star_R ( 6) * V_Star_R ( 3)
     1          + B0 ( 2) * V_Star_R ( 2) - B0 ( 1) * V_Star_R ( 3)

      F_Star_R ( 8) = F_Star_R ( 8) - V_Star_R ( 6) * V_Star_R ( 4)
     1          + B0 ( 3) * V_Star_R ( 2) - B0 ( 1) * V_Star_R ( 4)

! ------------------------------------------------------------
! "U_Star_L --> V_Star_L --> F_Star_L" gives us the flux in the +ve x-direction.

! --------------------
! Building "V_Star_L".

      U_Star_L ( 1) = AMAX1 ( U_Star_L ( 1), rhofloor)

      ob_rho  = U_Star_L ( 1) / ( U_Star_L ( 1)**2 + rhofloor**2)

      V_Star_L ( 1)    = U_Star_L ( 1)               ! density
      V_Star_L ( 2: 4) = U_Star_L ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_Star_L ( 6: 8) = U_Star_L ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_Star_L ( NVAR + 1: ) = U_Star_L ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_Star_L ( 2)**2 + V_Star_L ( 3)**2 + V_Star_L ( 4)**2
      B2 = V_Star_L ( 6)**2 + V_Star_L ( 7)**2 + V_Star_L ( 8)**2

      V_dot_B  = V_Star_L ( 2) * V_Star_L ( 6)
     1         + V_Star_L ( 3) * V_Star_L ( 7)
     1         + V_Star_L ( 4) * V_Star_L ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_Star_L ( 5) = ( gamma - 1.0) * ( U_Star_L ( 5)
     1             - 0.5 * V_Star_L ( 1) * v2  - ob_8pi * B2  )

      V_Star_L ( 5) = AMAX1 ( V_Star_L ( 5), prsfloor)  ! gas pressure

      U_Star_L ( 5) = V_Star_L ( 5) / ( gamma - 1.0)
     1        + 0.5 * V_Star_L ( 1) * v2 + ob_8pi * B2
#else
      V_Star_L ( 5) = isospeed**2 * U_Star_L ( 1)
     1          + barotropic_stiffness * U_Star_L ( 1)**gamma_euler
#endif


! --------------------
! Building "F_Star_L".

      B0_dot_B  = B0 ( 1) * V_Star_L ( 6)
     1          + B0 ( 2) * V_Star_L ( 7)
     1          + B0 ( 3) * V_Star_L ( 8)

      P_tot  = V_Star_L ( 5) + ob_8pi * B2
     1                       + ob_4pi * B0_dot_B

! ----------

      F_Star_L ( :) = V_Star_L ( 2) * U_Star_L ( :)


      F_Star_L ( 2) = F_Star_L ( 2) + P_tot
     1          - ob_4pi * V_Star_L ( 6) * ( V_Star_L ( 6) + B0 ( 1))
     1          - ob_4pi * V_Star_L ( 6) * B0 ( 1)

      F_Star_L ( 3) = F_Star_L ( 3)
     1          - ob_4pi * V_Star_L ( 6) * V_Star_L ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_Star_L ( 7)
     1                     + B0 ( 2) * V_Star_L ( 6))

      F_Star_L ( 4) = F_Star_L ( 4)
     1          - ob_4pi * V_Star_L ( 6) * V_Star_L ( 8)
     1          - ob_4pi * ( B0 ( 1) * V_Star_L ( 8)
     1                     + B0 ( 3) * V_Star_L ( 6))

#if ( CHOOSE_ADIABATIC == 1)
      F_Star_L ( 5) = F_Star_L ( 5)
     1          + P_tot  * V_Star_L ( 2)
     1          - ob_4pi * ( V_Star_L ( 6) + B0 ( 1))* V_dot_B
#else
      F_Star_L ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      F_Star_L ( 6) = 0.0

      F_Star_L ( 7) = F_Star_L ( 7) - V_Star_L ( 6) * V_Star_L ( 3)
     1          + B0 ( 2) * V_Star_L ( 2) - B0 ( 1) * V_Star_L ( 3)

      F_Star_L ( 8) = F_Star_L ( 8) - V_Star_L ( 6) * V_Star_L ( 4)
     1          + B0 ( 3) * V_Star_L ( 2) - B0 ( 1) * V_Star_L ( 4)

! ------------------------------------------------------------
! "U_Star_U --> V_Star_U --> G_Star_U" gives us the flux in the neta2-direction.

! --------------------
! Building "V_Star_U".

      U_Star_U ( 1) = AMAX1 ( U_Star_U ( 1), rhofloor)

      ob_rho  = U_Star_U ( 1) / ( U_Star_U ( 1)**2 + rhofloor**2)

      V_Star_U ( 1)    = U_Star_U ( 1)               ! density
      V_Star_U ( 2: 4) = U_Star_U ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_Star_U ( 6: 8) = U_Star_U ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_Star_U ( NVAR + 1: ) = U_Star_U ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_Star_U ( 2)**2 + V_Star_U ( 3)**2 + V_Star_U ( 4)**2
      B2 = V_Star_U ( 6)**2 + V_Star_U ( 7)**2 + V_Star_U ( 8)**2

      V_dot_B  = V_Star_U ( 2) * V_Star_U ( 6)
     1         + V_Star_U ( 3) * V_Star_U ( 7)
     1         + V_Star_U ( 4) * V_Star_U ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_Star_U ( 5) = ( gamma - 1.0) * ( U_Star_U ( 5)
     1             - 0.5 * V_Star_U ( 1) * v2  - ob_8pi * B2  )

      V_Star_U ( 5) = AMAX1 ( V_Star_U ( 5), prsfloor)  ! gas pressure

      U_Star_U ( 5) = V_Star_U ( 5) / ( gamma - 1.0)
     1        + 0.5 * V_Star_U ( 1) * v2 + ob_8pi * B2
#else
      V_Star_U ( 5) = isospeed**2 * U_Star_U ( 1)
     1          + barotropic_stiffness * U_Star_U ( 1)**gamma_euler
#endif

! ----------
! "vphi, bphi" are the velocity and magnetic field in the neta2-direction.

      eng = U_star_U ( 5)

      rho = V_Star_U ( 1)
      prs = V_Star_U ( 5)

      vx = V_Star_U ( 2) ;   vy = V_Star_U ( 3) ;   vz = V_Star_U ( 4)
      bx = V_Star_U ( 6) ;   by = V_Star_U ( 7) ;   bz = V_Star_U ( 8)

      vphi = c_phi * vx  + s_phi* vy
      bphi = c_phi * bx  + s_phi* by

      B0phi = c_phi * B0 ( 1) + s_phi* B0 ( 2)

      B0_dot_B  = B0 ( 1) * V_Star_U ( 6)
     1          + B0 ( 2) * V_Star_U ( 7)
     1          + B0 ( 3) * V_Star_U ( 8)

      P_tot  = V_Star_U ( 5) + ob_8pi * B2
     1                       + ob_4pi * B0_dot_B

! ----------
! Building "G_Star_U".

      G_Star_U ( 1) = rho * vphi

      G_Star_U ( 2) = rho * vphi * vx + P_tot * c_phi
     1                - ob_4pi * ( bx * (B0phi + bphi) + B0 ( 1) * bphi)

      G_Star_U ( 3) = rho * vphi * vy + P_tot * s_phi
     1                - ob_4pi * ( by * (B0phi + bphi) + B0 ( 2) * bphi)

      G_Star_U ( 4) = rho * vphi * vz
     1                - ob_4pi * ( bz * (B0phi + bphi) + B0 ( 3) * bphi)

#if ( CHOOSE_ADIABATIC == 1)
      G_Star_U ( 5) = eng * vphi + P_tot * vphi
     1              - V_dot_B * (B0phi + bphi) * ob_4pi
#else
      G_Star_U ( 5) = 0.0
#endif

      G_Star_U ( 6) = vphi * ( bx + B0 ( 1)) - vx * (B0phi + bphi)
      G_Star_U ( 7) = vphi * ( by + B0 ( 2)) - vy * (B0phi + bphi)
      G_Star_U ( 8) = vphi * ( bz + B0 ( 3)) - vz * (B0phi + bphi)

! ----------
! Now rotate "F_Star_U, G_Star_U" back to the Cartesian space in which the problem
! is being solved.

        F1 = c_phi * G_Star_U + s_phi * F_Star_U
        F2 = s_phi * G_Star_U - c_phi * F_Star_U

        F_Star_U = F1
        G_Star_U = F2

! ------------------------------------------------------------
! "U_Star_D --> V_Star_D --> G_Star_D" gives us the flux in the neta2-direction.

! --------------------
! Building "V_Star_D".

      U_Star_D ( 1) = AMAX1 ( U_Star_D ( 1), rhofloor)

      ob_rho  = U_Star_D ( 1) / ( U_Star_D ( 1)**2 + rhofloor**2)

      V_Star_D ( 1)    = U_Star_D ( 1)               ! density
      V_Star_D ( 2: 4) = U_Star_D ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_Star_D ( 6: 8) = U_Star_D ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_Star_D ( NVAR + 1: ) = U_Star_D ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_Star_D ( 2)**2 + V_Star_D ( 3)**2 + V_Star_D ( 4)**2
      B2 = V_Star_D ( 6)**2 + V_Star_D ( 7)**2 + V_Star_D ( 8)**2

      V_dot_B  = V_Star_D ( 2) * V_Star_D ( 6)
     1         + V_Star_D ( 3) * V_Star_D ( 7)
     1         + V_Star_D ( 4) * V_Star_D ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_Star_D ( 5) = ( gamma - 1.0) * ( U_Star_D ( 5)
     1              - 0.5 * V_Star_D ( 1) * v2  - ob_8pi * B2  )

      V_Star_D ( 5) = AMAX1 ( V_Star_D ( 5), prsfloor)  ! gas pressure

      U_Star_D ( 5) = V_Star_D ( 5) / ( gamma - 1.0)
     1        + 0.5 * V_Star_D ( 1) * v2 + ob_8pi * B2
#else
      V_Star_D ( 5) = isospeed**2 * U_Star_D ( 1)
     1          + barotropic_stiffness * U_Star_D ( 1)**gamma_euler
#endif

! ----------
! "vphi, bphi" are the velocity and magnetic field in the neta2-direction.

      eng = U_star_D ( 5)

      rho = V_Star_D ( 1)
      prs = V_Star_D ( 5)

      vx = V_Star_D ( 2) ;   vy = V_Star_D ( 3) ;   vz = V_Star_D ( 4)
      bx = V_Star_D ( 6) ;   by = V_Star_D ( 7) ;   bz = V_Star_D ( 8)

      vphi = c_phi * vx  + s_phi* vy
      bphi = c_phi * bx  + s_phi* by

      B0phi = c_phi * B0 ( 1) + s_phi* B0 ( 2)

      B0_dot_B  = B0 ( 1) * V_Star_D ( 6)
     1          + B0 ( 2) * V_Star_D ( 7)
     1          + B0 ( 3) * V_Star_D ( 8)

      P_tot  = V_Star_D ( 5) + ob_8pi * B2
     1                       + ob_4pi * B0_dot_B

! ----------
! Building "G_Star_D".

      G_Star_D ( 1) = rho * vphi

      G_Star_D ( 2) = rho * vphi * vx + P_tot * c_phi
     1                - ob_4pi * ( bx * (B0phi + bphi) + B0 ( 1) * bphi)

      G_Star_D ( 3) = rho * vphi * vy + P_tot * s_phi
     1                - ob_4pi * ( by * (B0phi + bphi) + B0 ( 2) * bphi)

      G_Star_D ( 4) = rho * vphi * vz
     1                - ob_4pi * ( bz * (B0phi + bphi) + B0 ( 3) * bphi)

#if ( CHOOSE_ADIABATIC == 1)
      G_Star_D ( 5) = eng * vphi + P_tot * vphi
     1              - V_dot_B * (B0phi + bphi) * ob_4pi
#else
      G_Star_D ( 5) = 0.0
#endif

      G_Star_D ( 6) = vphi * ( bx + B0 ( 1)) - vx * (B0phi + bphi)
      G_Star_D ( 7) = vphi * ( by + B0 ( 2)) - vy * (B0phi + bphi)
      G_Star_D ( 8) = vphi * ( bz + B0 ( 3)) - vz * (B0phi + bphi)

! ----------
! Now rotate "F_Star_D, G_Star_D" back to the Cartesian space in which the problem
! is being solved.

        F1 = c_phi * G_Star_D + s_phi * F_Star_D
        F2 = s_phi * G_Star_D - c_phi * F_Star_D

        F_Star_D = F1
        G_Star_D = F2

! ----------------------------------------------------------------------
! STEP 8) Using the multidimensional wave model, we can identify situations
! where the resolved state is supersonic. If indeed it is supersonic, there is
! no need to construct the resolved state in the multidimensional Riemann solver.

! We still need to go through STEP 11 for a back rotation before exit.

      IF ( ( S_R .LE. 0.0) .AND. ( S_U .LE. 0.0)) THEN

        U_Dbl_Star = U_RU  ;  F_Dbl_Star = F_RU  ;  G_Dbl_Star = G_RU
        GO TO 11

      ELSE IF ( ( S_L .GE. 0.0) .AND. ( S_U .LE. 0.0)) THEN

        U_Dbl_Star = U_LU  ;  F_Dbl_Star = F_LU  ;  G_Dbl_Star = G_LU
        GO TO 11

      ELSE IF ( ( S_L .GE. 0.0) .AND. ( S_D .GE. 0.0)) THEN

        U_Dbl_Star = U_LD  ;  F_Dbl_Star = F_LD  ;  G_Dbl_Star = G_LD
        GO TO 11

      ELSE IF ( ( S_R .LE. 0.0) .AND. ( S_D .GE. 0.0)) THEN

        U_Dbl_Star = U_RD  ;  F_Dbl_Star = F_RD  ;  G_Dbl_Star = G_RD
        GO TO 11

      ELSE IF ( ( S_R .LE. 0.0) .AND. ( S_U .GE. 0.0) .AND.
     1                                ( S_D .LE. 0.0)) THEN

        U_Dbl_Star = U_Star_R
        F_Dbl_Star = F_Star_R  ;  G_Dbl_Star = G_Star_R
        GO TO 11

      ELSE IF ( ( S_L .GE. 0.0) .AND. ( S_U .GE. 0.0) .AND.
     1                                ( S_D .LE. 0.0)) THEN

        U_Dbl_Star = U_Star_L
        F_Dbl_Star = F_Star_L  ;  G_Dbl_Star = G_Star_L
        GO TO 11

      ELSE IF ( ( S_R .GE. 0.0) .AND. ( S_L .LE. 0.0) .AND.
     1                                ( S_D .GE. 0.0)) THEN

        U_Dbl_Star = U_Star_D
        F_Dbl_Star = F_Star_D   ! c_phi * G_Star_D + s_phi * F_Star_D
        G_Dbl_Star = G_Star_D   ! s_phi * G_Star_D - c_phi * F_Star_D
        GO TO 11

      ELSE IF ( ( S_R .GE. 0.0) .AND. ( S_L .LE. 0.0) .AND.
     1                                ( S_U .LE. 0.0)) THEN

        U_Dbl_Star = U_Star_U
        F_Dbl_Star = F_Star_U   ! c_phi * G_Star_U + s_phi * F_Star_U
        G_Dbl_Star = G_Star_U   ! s_phi * G_Star_U - c_phi * F_Star_U
        GO TO 11

      END IF

! ----------------------------------------------------------------------
! STEP 9) Now we identify the limits of the multidimensional wave model.
! (x1, y1), (x2, y2), (x3, y3) and (x4, y4) are labeled counterclockwise.

! We also obtain the lengths of the segments that make up the boundary of
! the multidimensional wave model, and the centers of those segments.

! ----------
! Obtain vertices of multidimensional wave model.

      x1 = ( S_R + S_U * c_phi) / s_phi  ;  y1 = S_U

      x2 = ( S_L + S_U * c_phi) / s_phi  ;  y2 = S_U

      x3 = ( S_L + S_D * c_phi) / s_phi  ;  y3 = S_D

      x4 = ( S_R + S_D * c_phi) / s_phi  ;  y4 = S_D

! ----------
! Obtain centroid of multidimensional wave model; which is shaped like a
! parallelogram here.

      xc = 0.25 * ( x1 + x2 + x3 + x4)
      yc = 0.25 * ( y1 + y2 + y3 + y4)

! ----------
! Obtain side lenghts of multidimensional wave model; and the centers
! of those length segments.

      ell12 = SQRT ( (x2 - x1)**2 + (y2 - y1)**2)
      x12 = 0.5 * ( x1 + x2)  ;  y12 = 0.5 * ( y1 + y2)

      ell23 = SQRT ( (x3 - x2)**2 + (y3 - y2)**2)
      x23 = 0.5 * ( x2 + x3)  ;  y23 = 0.5 * ( y2 + y3)

      ell34 = SQRT ( (x4 - x3)**2 + (y4 - y3)**2)
      x34 = 0.5 * ( x3 + x4)  ;  y34 = 0.5 * ( y3 + y4)

      ell41 = SQRT ( (x1 - x4)**2 + (y1 - y4)**2)
      x41 = 0.5 * ( x4 + x1)  ;  y41 = 0.5 * ( y4 + y1)

! ----------
! Obtain area of multidimensional wave model.

      Area = ABS ( x4 - x3) * ABS ( y1 - y4)  ! Area of parallelogram wave model.

! ----------------------------------------------------------------------
! STEP 10) If we get here, the resolved state is subsonic. In that case, we build
! it here along with the associated fluxes.

! We still need to go through STEP 11 for a back rotation before exit.

! ------------------------------------------------------------
! Build states and fluxes that lie just outside the multidimensional wave model.
! The fluxes are outward-going fluxes at those boundaries.

      U1 = U_Star_R
      F1 = s_phi * F_Star_R - c_phi * G_Star_R  ! Flux in tau2-direction.

      U2 = U_Star_U
      F2 = G_Star_U  ! Flux in y-direction.

      U3 = U_Star_L
      F3 = - s_phi * F_Star_L + c_phi * G_Star_L  ! Flux in - tau2-direction.

      U4 = U_Star_D
      F4 = - G_Star_D  ! Flux in - y-direction.

! ------------------------------------------------------------
! Build "U_Dbl_Star".

      a = - 0.5 / Area

      U_Dbl_Star = a * ( ( F1 - S_R * U1) * ell41
     1                 + ( F2 - S_U * U2) * ell12
     1                 + ( F3 + S_L * U3) * ell23
     1                 + ( F4 + S_D * U4) * ell34 )

! ------------------------------------------------------------
! Build "F_Dbl_Star".

      a = 1.0 / Area

      F_Dbl_Star = xc * U_Dbl_Star
     1           + a * ( ( F1 - S_R * U1) * ell41 * ( x41 - xc)
     1                 + ( F2 - S_U * U2) * ell12 * ( x12 - xc)
     1                 + ( F3 + S_L * U3) * ell23 * ( x23 - xc)
     1                 + ( F4 + S_D * U4) * ell34 * ( x34 - xc) )

! ------------------------------------------------------------
! Build "G_Dbl_Star".

      a = 1.0 / Area

      G_Dbl_Star = yc * U_Dbl_Star
     1           + a * ( ( F1 - S_R * U1) * ell41 * ( y41 - yc)
     1                 + ( F2 - S_U * U2) * ell12 * ( y12 - yc)
     1                 + ( F3 + S_L * U3) * ell23 * ( y23 - yc)
     1                 + ( F4 + S_D * U4) * ell34 * ( y34 - yc) )

! ----------------------------------------------------------------------
! STEP 11) Back rotate and exit subroutine. We can come to this point with a 1D
! Riemann solver solution from STEP 8, or with a 2D Riemann solver solution
! from STEP 10. Either way, the inner components ("2,3,4" and "6,7,8") of
! "U_Dbl_Star, F_Dbl_Star, G_Dbl_Star" have to be rotated back into the
! global space. In other words, we undo STEP 3.

! In STEP 2 we have obtained "np1 (3), np2 (3), np3 (3)" for the
! reverse transformation.

 11   CONTINUE

! ----------
! Build "E_z" before enforcing any back rotation.

! The fluxes (without back rotation) are in the "xprime, yprime" frame; so it
! is easy to read off "E_z".

      E_z = 0.5 * ( G_Dbl_Star ( 6) - F_Dbl_Star ( 7))

! ----------
! rotating momentum components in "U_Dbl_Star" to the global frame.

      a = U_Dbl_Star ( 2)
      b = U_Dbl_Star ( 3)
      c = U_Dbl_Star ( 4)

      U_Dbl_Star ( 2) = a * np1 ( 1) + b * np1 ( 2) + c * np1 ( 3)
      U_Dbl_Star ( 3) = a * np2 ( 1) + b * np2 ( 2) + c * np2 ( 3)
      U_Dbl_Star ( 4) = a * np3 ( 1) + b * np3 ( 2) + c * np3 ( 3)

! rotating magnetic field components in "U_Dbl_Star" to the global frame.
      a = U_Dbl_Star ( 6)
      b = U_Dbl_Star ( 7)
      c = U_Dbl_Star ( 8)

      U_Dbl_Star ( 6) = a * np1 ( 1) + b * np1 ( 2) + c * np1 ( 3)
      U_Dbl_Star ( 7) = a * np2 ( 1) + b * np2 ( 2) + c * np2 ( 3)
      U_Dbl_Star ( 8) = a * np3 ( 1) + b * np3 ( 2) + c * np3 ( 3)

! ----------
! rotating momentum components in "F_Dbl_Star" to the global frame.

      a = F_Dbl_Star ( 2)
      b = F_Dbl_Star ( 3)
      c = F_Dbl_Star ( 4)

      F_Dbl_Star ( 2) = a * np1 ( 1) + b * np1 ( 2) + c * np1 ( 3)
      F_Dbl_Star ( 3) = a * np2 ( 1) + b * np2 ( 2) + c * np2 ( 3)
      F_Dbl_Star ( 4) = a * np3 ( 1) + b * np3 ( 2) + c * np3 ( 3)

! rotating magnetic field components in "F_Dbl_Star" to the global frame.
      a = F_Dbl_Star ( 6)
      b = F_Dbl_Star ( 7)
      c = F_Dbl_Star ( 8)

      F_Dbl_Star ( 6) = a * np1 ( 1) + b * np1 ( 2) + c * np1 ( 3)
      F_Dbl_Star ( 7) = a * np2 ( 1) + b * np2 ( 2) + c * np2 ( 3)
      F_Dbl_Star ( 8) = a * np3 ( 1) + b * np3 ( 2) + c * np3 ( 3)

! ----------
! rotating momentum components in "G_Dbl_Star" to the global frame.

      a = G_Dbl_Star ( 2)
      b = G_Dbl_Star ( 3)
      c = G_Dbl_Star ( 4)

      G_Dbl_Star ( 2) = a * np1 ( 1) + b * np1 ( 2) + c * np1 ( 3)
      G_Dbl_Star ( 3) = a * np2 ( 1) + b * np2 ( 2) + c * np2 ( 3)
      G_Dbl_Star ( 4) = a * np3 ( 1) + b * np3 ( 2) + c * np3 ( 3)

! rotating magnetic field components in "G_Dbl_Star" to the global frame.
      a = G_Dbl_Star ( 6)
      b = G_Dbl_Star ( 7)
      c = G_Dbl_Star ( 8)

      G_Dbl_Star ( 6) = a * np1 ( 1) + b * np1 ( 2) + c * np1 ( 3)
      G_Dbl_Star ( 7) = a * np2 ( 1) + b * np2 ( 2) + c * np2 ( 3)
      G_Dbl_Star ( 8) = a * np3 ( 1) + b * np3 ( 2) + c * np3 ( 3)

! ----------------------------------------------------------------------

      END SUBROUTINE TwoD_RS_4State
