












      SUBROUTINE TwoD_RS_3State ( n_cc_components,

     1           pi, rhofloor, prsfloor, smallnum, gamma,
     1           isospeed, barotropic_stiffness, gamma_euler,

     1           n1, n2, n3, B0,

     1           U_1, U_2, U_3,
     1           U_Dbl_Star, F_Dbl_Star, G_Dbl_Star, E_z )

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara and Sethu Subramanian

! Takes in two vectors "n1, n2" that don't have to be mutually orthogonal.
! Makes a primed space that is locally Cartesian with the z-axis perpendicular
! to the plane in which the 2DRS is solved. Let the third vector in the xy-plane
! be called "n12".

! The three input states "U_1, U_2, U_3" have their interior components
! rotated into that locally Cartesian space. Only then can we begin a solution
! of the 2D Riemann problem.
! Note that "U_1, U_2, U_3" have a counterclockwise and cyclic arrangement.
! "U_1" is the state between "n1, n2"; "U_2" is the state between "n2, n12";
! "U_3" is the state between "n12, n1".

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

      REAL, DIMENSION ( NVAR + NFLUID) :: U_1, U_2, U_3,
     1      U_Dbl_Star, F_Dbl_Star, G_Dbl_Star

! ----------------------------------------------------------------------

! Variables past this point are internal to the subroutine.

      REAL fourpi, ob_4pi, ob_8pi,
     1     a, b, c, phi, s_phi, c_phi, s_30, c_30,
     1     ny (3), n12 (3), np1 (3), np2 (3), np3 (3),

     1     v2, b2, v_dot_b, b0_dot_b, P_tot, ob_rho,
     1     rho, prs, vx, vy, vz, bx, by, bz, vphi, bphi, eng, B0phi,

     1     cs2, vel_long, ob_4prho, va2_long, va2_tot, vf,

     1     S_L1, S_R1, S_L2, S_R2, S_L12, S_R12, S1, S2,

     1     n_x, n_y, vn, Bn, B0n, x1, x2, x3, y1, y2, y3,

     1     ell31, n31_x, n31_y, x31, y31, S_31,
     1     ell12, n12_x, n12_y, x12, y12, S_12,
     1     ell23, n23_x, n23_y, x23, y23, S_23, xc, yc, Area


      REAL, DIMENSION ( NVAR + NFLUID) :: V_1, V_2, V_3,
     1     F_1, F_2, F_3, G_1, G_2, G_3,
     1     U1, U2, U3, F1, F2, F3, F, G,
     1     U_Star_13, F_Star_13, V_Star_13, G_Star_13,
     1     U_Star_21, F_Star_21, V_Star_21, G_Star_21,
     1     U_Star_32, F_Star_32, V_Star_32, G_Star_32

! ----------------------------------------------------------------------
! STEP 1) Obtain the constants.

      fourpi = 4.0 * pi
      ob_4pi = 1.0 / fourpi
      ob_8pi = 0.5 * ob_4pi

      s_30 = 0.5
      c_30 = SQRT ( 3.0) / 2.0

! ----------------------------------------------------------------------
! STEP 2) Construct unit vectors "n1 (3), ny (3), n3 (3)". Also construct the
! unit vectors that give us the reverse transform via "np1 (3), np2 (3), np3 (3)".

! "n1 (3)" will become our x-direction in the xy-plane. We can always define
! "ny (3)" to be the y-direction in this plane.

! Note that now "n1 (3), n2 (3), n12 (3)" form a cyclic set with 120 degrees
! between them. The notation here is different from the notation in the 4-state
! multid RS.

! "n1" and "n2", as they come in, are in a global basis. Besides, they may
! not be orthogonal.

      a = 1.0 / SQRT ( n1 ( 1)**2 + n1 ( 2)**2 + n1 ( 3)**2 )
      n1 ( :) = n1 ( :) * a

      a = 1.0 / SQRT ( n2 ( 1)**2 + n2 ( 2)**2 + n2 ( 3)**2 )
      n2 ( :) = n2 ( :) * a

      n3 ( 1) = n1 ( 2) * n2 ( 3) - n2 ( 2) * n1 ( 3)
      n3 ( 2) = n1 ( 3) * n2 ( 1) - n2 ( 3) * n1 ( 1)
      n3 ( 3) = n1 ( 1) * n2 ( 2) - n2 ( 1) * n1 ( 2)

      a = 1.0 / SQRT ( n3 ( 1)**2 + n3 ( 2)**2 + n3 ( 3)**2 )
      n3 ( :) = n3 ( :) * a

! By this point, we have "n1, n2, n3" expressed in terms of a global basis.
! "n1, n2" may not be orthogonal. We construct "ny" which is perpendicular
! to "n1, n3" and in the direction of "n2".

      c_phi = n1 ( 1) * n2 ( 1) + n1 ( 2) * n2 ( 2) + n1 ( 3) * n2 ( 3)
      phi = ACOS ( c_phi)
      s_phi = SQRT ( 1.0 - c_phi**2)

      ny ( :) = ( n2 ( :) - c_phi * n1 ( :) ) / s_phi

! Since "n1 (3), n2 (3), n12 (3)" form a cyclic set with 120 degrees
! between them, we can obtain "n12 (3)" by a vector construction.

      n12 ( :) = n1 ( :) + n2 ( :)
      a = 1.0 / SQRT ( n12 ( 1)**2 + n12 ( 2)**2 + n12 ( 3)**2 )
      n12 ( :) = - n12 ( :) * a

! By this point, "n1 (3), ny (3), n3 (3)" are orthonormal unit vectors
! which give us a local "xprime, yprime, zprime" frame in which the 2D Riemann
! problem can be solved in the "xprime, yprime" plane. They are expressed in terms
! of a global basis set, so they can be used to transform from global "x, y, z"
! frame to a rotated "xprime, yprime, zprime" frame.

! We also wish to have the reverse transformation. The inverse of an orthonormal
! matrix is indeed its transverse. So "np1 (3), np2 (3), np3 (3)" give us that
! reverse transformation.

      np1 ( 1) = n1 ( 1) ; np1 ( 2) = ny ( 1) ; np1 ( 3) = n3 ( 1)
      np2 ( 1) = n1 ( 2) ; np2 ( 2) = ny ( 2) ; np2 ( 3) = n3 ( 2)
      np3 ( 1) = n1 ( 3) ; np3 ( 2) = ny ( 3) ; np3 ( 3) = n3 ( 3)

! By this point "np1 (3), np2 (3), np3 (3)" give us the ability to go from
! the local "xprime, yprime, zprime" frame to the global "x, y, z" frame.

! ----------------------------------------------------------------------
! STEP 3) Now make internal in-place rotations to the "2,3,4" and "6,7,8"
! components of the three input states "U_1, U_2, U_3".
! This makes all the velocities and magnetic fields defined in the
! "xprime, yprime, zprime" frame.

! Also rotate "B0" into the "xprime, yprime, zprime" frame.

! ----------
! rotating background magnetic field "B0" to the local frame.
      a = B0 ( 1)
      b = B0 ( 2)
      c = B0 ( 3)

      B0 ( 1) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      B0 ( 2) = a * ny ( 1) + b * ny ( 2) + c * ny ( 3)
      B0 ( 3) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! ----------
! rotating momentum components in "U_1" to the local frame.
      a = U_1 ( 2)
      b = U_1 ( 3)
      c = U_1 ( 4)

      U_1 ( 2) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_1 ( 3) = a * ny ( 1) + b * ny ( 2) + c * ny ( 3)
      U_1 ( 4) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! rotating magnetic field components in "U_1" to the local frame.
      a = U_1 ( 6)
      b = U_1 ( 7)
      c = U_1 ( 8)

      U_1 ( 6) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_1 ( 7) = a * ny ( 1) + b * ny ( 2) + c * ny ( 3)
      U_1 ( 8) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! ----------
! rotating momentum components in "U_2" to the local frame.
      a = U_2 ( 2)
      b = U_2 ( 3)
      c = U_2 ( 4)

      U_2 ( 2) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_2 ( 3) = a * ny ( 1) + b * ny ( 2) + c * ny ( 3)
      U_2 ( 4) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! rotating magnetic field components in "U_2" to the local frame.
      a = U_2 ( 6)
      b = U_2 ( 7)
      c = U_2 ( 8)

      U_2 ( 6) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_2 ( 7) = a * ny ( 1) + b * ny ( 2) + c * ny ( 3)
      U_2 ( 8) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! ----------
! rotating momentum components in "U_3" to the local frame.
      a = U_3 ( 2)
      b = U_3 ( 3)
      c = U_3 ( 4)

      U_3 ( 2) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_3 ( 3) = a * ny ( 1) + b * ny ( 2) + c * ny ( 3)
      U_3 ( 4) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! rotating magnetic field components in "U_3" to the local frame.
      a = U_3 ( 6)
      b = U_3 ( 7)
      c = U_3 ( 8)

      U_3 ( 6) = a * n1 ( 1) + b * n1 ( 2) + c * n1 ( 3)
      U_3 ( 7) = a * ny ( 1) + b * ny ( 2) + c * ny ( 3)
      U_3 ( 8) = a * n3 ( 1) + b * n3 ( 2) + c * n3 ( 3)

! ----------------------------------------------------------------------
! STEP 4) From "U_1" etc. build "V_1, F_1, G_1".

! The fluxes are defined in the "xprime, yprime, zprime" frame.

! ----------------------------------------
! Building "V_1".

      U_1 ( 1) = AMAX1 ( U_1 ( 1), rhofloor)

      ob_rho  = U_1 ( 1) / ( U_1 ( 1)**2 + rhofloor**2)

      V_1 ( 1)    = U_1 ( 1)               ! density
      V_1 ( 2: 4) = U_1 ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_1 ( 6: 8) = U_1 ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_1 ( NVAR + 1: ) = U_1 ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_1 ( 2)**2 + V_1 ( 3)**2 + V_1 ( 4)**2
      B2 = V_1 ( 6)**2 + V_1 ( 7)**2 + V_1 ( 8)**2

      V_dot_B  = V_1 ( 2) * V_1 ( 6)
     1         + V_1 ( 3) * V_1 ( 7)
     1         + V_1 ( 4) * V_1 ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_1 ( 5) = ( gamma - 1.0) * ( U_1 ( 5)
     1             - 0.5 * V_1 ( 1) * v2  - ob_8pi * B2  )

      V_1 ( 5) = AMAX1 ( V_1 ( 5), prsfloor)  ! gas pressure

      U_1 ( 5) = V_1 ( 5) / ( gamma - 1.0)
     1         + 0.5 * V_1 ( 1) * v2 + ob_8pi * B2
#else
      V_1 ( 5) = isospeed**2 * U_1 ( 1)
     1          + barotropic_stiffness * U_1 ( 1)**gamma_euler
#endif

! --------------------
! Building "F_1, G_1".

      B0_dot_B  = B0 ( 1) * V_1 ( 6)
     1          + B0 ( 2) * V_1 ( 7)
     1          + B0 ( 3) * V_1 ( 8)

      P_tot  = V_1 ( 5) + ob_8pi * B2
     1                   + ob_4pi * B0_dot_B

! ----------
! Building "F_1".

      F_1 ( :) = V_1 ( 2) * U_1 ( :)


      F_1 ( 2) = F_1 ( 2) + P_tot
     1          - ob_4pi * V_1 ( 6) * ( V_1 ( 6) + B0 ( 1))
     1          - ob_4pi * V_1 ( 6) * B0 ( 1)

      F_1 ( 3) = F_1 ( 3)
     1          - ob_4pi * V_1 ( 6) * V_1 ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_1 ( 7)
     1                     + B0 ( 2) * V_1 ( 6))

      F_1 ( 4) = F_1 ( 4)
     1          - ob_4pi * V_1 ( 6) * V_1 ( 8)
     1          - ob_4pi * ( B0 ( 1) * V_1 ( 8)
     1                     + B0 ( 3) * V_1 ( 6))

#if ( CHOOSE_ADIABATIC == 1)
      F_1 ( 5) = F_1 ( 5)
     1          + P_tot  * V_1 ( 2)
     1          - ob_4pi * ( V_1 ( 6) + B0 ( 1))* V_dot_B
#else
      F_1 ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      F_1 ( 6) = 0.0

      F_1 ( 7) = F_1 ( 7) - V_1 ( 6) * V_1 ( 3)
     1          + B0 ( 2) * V_1 ( 2) - B0 ( 1) * V_1 ( 3)

      F_1 ( 8) = F_1 ( 8) - V_1 ( 6) * V_1 ( 4)
     1          + B0 ( 3) * V_1 ( 2) - B0 ( 1) * V_1 ( 4)

! ----------
! Building "G_1".

      G_1 ( :) = V_1 ( 3) * U_1 ( :)

      G_1 ( 2) = G_1 ( 2)
     1          - ob_4pi * V_1 ( 6) * V_1 ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_1 ( 7)
     1                     + B0 ( 2) * V_1 ( 6))

      G_1 ( 3) = G_1 ( 3) + P_tot
     1           - ob_4pi * V_1 ( 7) * ( V_1 ( 7) + B0 ( 2))
     1           - ob_4pi * V_1 ( 7) * B0 ( 2)

      G_1 ( 4) = G_1 ( 4)
     1          - ob_4pi * V_1 ( 7) * V_1 ( 8)
     1          - ob_4pi * ( B0 ( 2) * V_1 ( 8)
     1                     + B0 ( 3) * V_1 ( 7))

#if ( CHOOSE_ADIABATIC == 1)
      G_1 ( 5) = G_1 ( 5)
     1          + P_tot  * V_1 ( 3)
     1          - ob_4pi * ( V_1 ( 7) + B0 ( 2))* V_dot_B
#else
      G_1 ( 5) = 0.0             ! isothermal has no energy flux.
#endif

      G_1 ( 6) = G_1 ( 6) - V_1 ( 2) * V_1 ( 7)
     1          + V_1 ( 3) * B0 ( 1) - V_1 ( 2) * B0 ( 2)

      G_1 ( 7) = 0.0

      G_1 ( 8) = G_1 ( 8) - V_1 ( 4) * V_1 ( 7)
     1          + V_1 ( 3) * B0 ( 3) - V_1 ( 4) * B0 ( 2)

! ----------------------------------------
! Building "V_2".

      U_2 ( 1) = AMAX1 ( U_2 ( 1), rhofloor)

      ob_rho  = U_2 ( 1) / ( U_2 ( 1)**2 + rhofloor**2)

      V_2 ( 1)    = U_2 ( 1)               ! density
      V_2 ( 2: 4) = U_2 ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_2 ( 6: 8) = U_2 ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_2 ( NVAR + 1: ) = U_2 ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_2 ( 2)**2 + V_2 ( 3)**2 + V_2 ( 4)**2
      B2 = V_2 ( 6)**2 + V_2 ( 7)**2 + V_2 ( 8)**2

      V_dot_B  = V_2 ( 2) * V_2 ( 6)
     1         + V_2 ( 3) * V_2 ( 7)
     1         + V_2 ( 4) * V_2 ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_2 ( 5) = ( gamma - 1.0) * ( U_2 ( 5)
     1             - 0.5 * V_2 ( 1) * v2  - ob_8pi * B2  )

      V_2 ( 5) = AMAX1 ( V_2 ( 5), prsfloor)  ! gas pressure

      U_2 ( 5) = V_2 ( 5) / ( gamma - 1.0)
     1         + 0.5 * V_2 ( 1) * v2 + ob_8pi * B2
#else
      V_2 ( 5) = isospeed**2 * U_2 ( 1)
     1          + barotropic_stiffness * U_2 ( 1)**gamma_euler
#endif

! --------------------
! Building "F_2, G_2".

      B0_dot_B  = B0 ( 1) * V_2 ( 6)
     1          + B0 ( 2) * V_2 ( 7)
     1          + B0 ( 3) * V_2 ( 8)

      P_tot  = V_2 ( 5) + ob_8pi * B2
     1                   + ob_4pi * B0_dot_B

! ----------
! Building "F_2".

      F_2 ( :) = V_2 ( 2) * U_2 ( :)


      F_2 ( 2) = F_2 ( 2) + P_tot
     1          - ob_4pi * V_2 ( 6) * ( V_2 ( 6) + B0 ( 1))
     1          - ob_4pi * V_2 ( 6) * B0 ( 1)

      F_2 ( 3) = F_2 ( 3)
     1          - ob_4pi * V_2 ( 6) * V_2 ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_2 ( 7)
     1                     + B0 ( 2) * V_2 ( 6))

      F_2 ( 4) = F_2 ( 4)
     1          - ob_4pi * V_2 ( 6) * V_2 ( 8)
     1          - ob_4pi * ( B0 ( 1) * V_2 ( 8)
     1                     + B0 ( 3) * V_2 ( 6))

#if ( CHOOSE_ADIABATIC == 1)
      F_2 ( 5) = F_2 ( 5)
     1          + P_tot  * V_2 ( 2)
     1          - ob_4pi * ( V_2 ( 6) + B0 ( 1))* V_dot_B
#else
      F_2 ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      F_2 ( 6) = 0.0

      F_2 ( 7) = F_2 ( 7) - V_2 ( 6) * V_2 ( 3)
     1          + B0 ( 2) * V_2 ( 2) - B0 ( 1) * V_2 ( 3)

      F_2 ( 8) = F_2 ( 8) - V_2 ( 6) * V_2 ( 4)
     1          + B0 ( 3) * V_2 ( 2) - B0 ( 1) * V_2 ( 4)

! ----------
! Building "G_2".

      G_2 ( :) = V_2 ( 3) * U_2 ( :)

      G_2 ( 2) = G_2 ( 2)
     1          - ob_4pi * V_2 ( 6) * V_2 ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_2 ( 7)
     1                     + B0 ( 2) * V_2 ( 6))

      G_2 ( 3) = G_2 ( 3) + P_tot
     1           - ob_4pi * V_2 ( 7) * ( V_2 ( 7) + B0 ( 2))
     1           - ob_4pi * V_2 ( 7) * B0 ( 2)

      G_2 ( 4) = G_2 ( 4)
     1          - ob_4pi * V_2 ( 7) * V_2 ( 8)
     1          - ob_4pi * ( B0 ( 2) * V_2 ( 8)
     1                     + B0 ( 3) * V_2 ( 7))

#if ( CHOOSE_ADIABATIC == 1)
      G_2 ( 5) = G_2 ( 5)
     1          + P_tot  * V_2 ( 3)
     1          - ob_4pi * ( V_2 ( 7) + B0 ( 2))* V_dot_B
#else
      G_2 ( 5) = 0.0             ! isothermal has no energy flux.
#endif

      G_2 ( 6) = G_2 ( 6) - V_2 ( 2) * V_2 ( 7)
     1          + V_2 ( 3) * B0 ( 1) - V_2 ( 2) * B0 ( 2)

      G_2 ( 7) = 0.0

      G_2 ( 8) = G_2 ( 8) - V_2 ( 4) * V_2 ( 7)
     1          + V_2 ( 3) * B0 ( 3) - V_2 ( 4) * B0 ( 2)


! ----------------------------------------
! Building "V_3".

      U_3 ( 1) = AMAX1 ( U_3 ( 1), rhofloor)

      ob_rho  = U_3 ( 1) / ( U_3 ( 1)**2 + rhofloor**2)

      V_3 ( 1)    = U_3 ( 1)               ! density
      V_3 ( 2: 4) = U_3 ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_3 ( 6: 8) = U_3 ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_3 ( NVAR + 1: ) = U_3 ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_3 ( 2)**2 + V_3 ( 3)**2 + V_3 ( 4)**2
      B2 = V_3 ( 6)**2 + V_3 ( 7)**2 + V_3 ( 8)**2

      V_dot_B  = V_3 ( 2) * V_3 ( 6)
     1         + V_3 ( 3) * V_3 ( 7)
     1         + V_3 ( 4) * V_3 ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_3 ( 5) = ( gamma - 1.0) * ( U_3 ( 5)
     1             - 0.5 * V_3 ( 1) * v2  - ob_8pi * B2  )

      V_3 ( 5) = AMAX1 ( V_3 ( 5), prsfloor)  ! gas pressure

      U_3 ( 5) = V_3 ( 5) / ( gamma - 1.0)
     1         + 0.5 * V_3 ( 1) * v2 + ob_8pi * B2
#else
      V_3 ( 5) = isospeed**2 * U_3 ( 1)
     1          + barotropic_stiffness * U_3 ( 1)**gamma_euler
#endif

! --------------------
! Building "F_3, G_3".

      B0_dot_B  = B0 ( 1) * V_3 ( 6)
     1          + B0 ( 2) * V_3 ( 7)
     1          + B0 ( 3) * V_3 ( 8)

      P_tot  = V_3 ( 5) + ob_8pi * B2
     1                   + ob_4pi * B0_dot_B

! ----------
! Building "F_3".

      F_3 ( :) = V_3 ( 2) * U_3 ( :)


      F_3 ( 2) = F_3 ( 2) + P_tot
     1          - ob_4pi * V_3 ( 6) * ( V_3 ( 6) + B0 ( 1))
     1          - ob_4pi * V_3 ( 6) * B0 ( 1)

      F_3 ( 3) = F_3 ( 3)
     1          - ob_4pi * V_3 ( 6) * V_3 ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_3 ( 7)
     1                     + B0 ( 2) * V_3 ( 6))

      F_3 ( 4) = F_3 ( 4)
     1          - ob_4pi * V_3 ( 6) * V_3 ( 8)
     1          - ob_4pi * ( B0 ( 1) * V_3 ( 8)
     1                     + B0 ( 3) * V_3 ( 6))

#if ( CHOOSE_ADIABATIC == 1)
      F_3 ( 5) = F_3 ( 5)
     1          + P_tot  * V_3 ( 2)
     1          - ob_4pi * ( V_3 ( 6) + B0 ( 1))* V_dot_B
#else
      F_3 ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      F_3 ( 6) = 0.0

      F_3 ( 7) = F_3 ( 7) - V_3 ( 6) * V_3 ( 3)
     1          + B0 ( 2) * V_3 ( 2) - B0 ( 1) * V_3 ( 3)

      F_3 ( 8) = F_3 ( 8) - V_3 ( 6) * V_3 ( 4)
     1          + B0 ( 3) * V_3 ( 2) - B0 ( 1) * V_3 ( 4)

! ----------
! Building "G_3".

      G_3 ( :) = V_3 ( 3) * U_3 ( :)

      G_3 ( 2) = G_3 ( 2)
     1          - ob_4pi * V_3 ( 6) * V_3 ( 7)
     1          - ob_4pi * ( B0 ( 1) * V_3 ( 7)
     1                     + B0 ( 2) * V_3 ( 6))

      G_3 ( 3) = G_3 ( 3) + P_tot
     1           - ob_4pi * V_3 ( 7) * ( V_3 ( 7) + B0 ( 2))
     1           - ob_4pi * V_3 ( 7) * B0 ( 2)

      G_3 ( 4) = G_3 ( 4)
     1          - ob_4pi * V_3 ( 7) * V_3 ( 8)
     1          - ob_4pi * ( B0 ( 2) * V_3 ( 8)
     1                     + B0 ( 3) * V_3 ( 7))

#if ( CHOOSE_ADIABATIC == 1)
      G_3 ( 5) = G_3 ( 5)
     1          + P_tot  * V_3 ( 3)
     1          - ob_4pi * ( V_3 ( 7) + B0 ( 2))* V_dot_B
#else
      G_3 ( 5) = 0.0             ! isothermal has no energy flux.
#endif

      G_3 ( 6) = G_3 ( 6) - V_3 ( 2) * V_3 ( 7)
     1          + V_3 ( 3) * B0 ( 1) - V_3 ( 2) * B0 ( 2)

      G_3 ( 7) = 0.0

      G_3 ( 8) = G_3 ( 8) - V_3 ( 4) * V_3 ( 7)
     1          + V_3 ( 3) * B0 ( 3) - V_3 ( 4) * B0 ( 2)

! ----------------------------------------------------------------------
! STEP 5) Build the extremal signal speeds for the 1D Riemann solvers.
! These are the speeds that go counterclockwise (left) and clockwise (right)
! relative to the unit vectors "n1 (3), n2 (3), n12 (3)".

! Relative to "n1 (3)" we have "S_L1, S_R1"; Relative to "n2 (3)" we have "S_L2, S_R2";
! Relative to "n12 (3)" we have "S_L12, S_R12".

! If flow is subsonic, we expect to have S_R1 < 0 < S_L1.

! ------------------------------
! Extremize "S_L1, S_R1" using "V_1 and V_3" in three steps.
! "V_1" serves as the left state and "V_3" as the right state in the RP.
! Longitudinal direction is in the ( 0, - 1) direction.

      S_L1 = 1.0e30
      S_R1 = - 1.0e30

! ----------
! Extremize "S_L1, S_R1" using "V_1".

      rho = V_1 ( 1)

      vx = V_1 ( 2)
      vy = V_1 ( 3)
      vz = V_1 ( 4)

      prs = V_1 ( 5)

      bx = V_1 ( 6) + B0 ( 1)
      by = V_1 ( 7) + B0 ( 2)
      bz = V_1 ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = - vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( - by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R1 = AMAX1 ( S_R1, vel_long + vf )
      S_L1 = AMIN1 ( S_L1, vel_long - vf )

! ----------
! Extremize "S_L1, S_R1" using "V_3".

      rho = V_3 ( 1)

      vx = V_3 ( 2)
      vy = V_3 ( 3)
      vz = V_3 ( 4)

      prs = V_3 ( 5)

      bx = V_3 ( 6) + B0 ( 1)
      by = V_3 ( 7) + B0 ( 2)
      bz = V_3 ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = - vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( - by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R1 = AMAX1 ( S_R1, vel_long + vf )
      S_L1 = AMIN1 ( S_L1, vel_long - vf )

! ----------
! Extremize "S_L1, S_R1" using "0.5 * ( V_1 + V_3)".

      rho = 0.5 * ( V_1 ( 1) + V_3 ( 1))

      vx = 0.5 * ( V_1 ( 2) + V_3 ( 2))
      vy = 0.5 * ( V_1 ( 3) + V_3 ( 3))
      vz = 0.5 * ( V_1 ( 4) + V_3 ( 4))

      prs = 0.5 * ( V_1 ( 5) + V_3 ( 5))

      bx = 0.5 * ( V_1 ( 6) + V_3 ( 6)) + B0 ( 1)
      by = 0.5 * ( V_1 ( 7) + V_3 ( 7)) + B0 ( 2)
      bz = 0.5 * ( V_1 ( 8) + V_3 ( 8)) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = - vy

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( - by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R1 = AMAX1 ( S_R1, vel_long + vf )
      S_L1 = AMIN1 ( S_L1, vel_long - vf )

! ------------------------------
! Extremize "S_L2, S_R2" using "V_2 and V_1" in three steps.
! "V_2" serves as the left state and "V_1" as the right state in the RP.
! Longitudinal direction is in the ( c_30, s_30) direction.

      S_L2 = 1.0e30
      S_R2 = - 1.0e30

! ----------
! Extremize "S_L2, S_R2" using "V_2".

      rho = V_2 ( 1)

      vx = V_2 ( 2)
      vy = V_2 ( 3)
      vz = V_2 ( 4)

      prs = V_2 ( 5)

      bx = V_2 ( 6) + B0 ( 1)
      by = V_2 ( 7) + B0 ( 2)
      bz = V_2 ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = ( c_30 * vx + s_30 * vy)

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( c_30 * bx + s_30 * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R2 = AMAX1 ( S_R2, vel_long + vf )
      S_L2 = AMIN1 ( S_L2, vel_long - vf )

! ----------
! Extremize "S_L2, S_R2" using "V_1".

      rho = V_1 ( 1)

      vx = V_1 ( 2)
      vy = V_1 ( 3)
      vz = V_1 ( 4)

      prs = V_1 ( 5)

      bx = V_1 ( 6) + B0 ( 1)
      by = V_1 ( 7) + B0 ( 2)
      bz = V_1 ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = ( c_30 * vx + s_30 * vy)

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( c_30 * bx + s_30 * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R2 = AMAX1 ( S_R2, vel_long + vf )
      S_L2 = AMIN1 ( S_L2, vel_long - vf )

! ----------
! Extremize "S_L2, S_R2" using "0.5 * ( V_2 + V_1)".

      rho = 0.5 * ( V_2 ( 1) + V_1 ( 1))

      vx = 0.5 * ( V_2 ( 2) + V_1 ( 2))
      vy = 0.5 * ( V_2 ( 3) + V_1 ( 3))
      vz = 0.5 * ( V_2 ( 4) + V_1 ( 4))

      prs = 0.5 * ( V_2 ( 5) + V_1 ( 5))

      bx = 0.5 * ( V_2 ( 6) + V_1 ( 6)) + B0 ( 1)
      by = 0.5 * ( V_2 ( 7) + V_1 ( 7)) + B0 ( 2)
      bz = 0.5 * ( V_2 ( 8) + V_1 ( 8)) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = ( c_30 * vx + s_30 * vy)

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( c_30 * bx + s_30 * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R2 = AMAX1 ( S_R2, vel_long + vf )
      S_L2 = AMIN1 ( S_L2, vel_long - vf )

! ------------------------------
! Extremize "S_L12, S_R12" using "V_3 and V_2" in three steps.
! "V_3" serves as the left state and "V_2" as the right state in the RP.
! Longitudinal direction is in the ( - c_30, s_30) direction.

      S_L12 = 1.0e30
      S_R12 = - 1.0e30

! ----------
! Extremize "S_L12, S_R12" using "V_3".

      rho = V_3 ( 1)

      vx = V_3 ( 2)
      vy = V_3 ( 3)
      vz = V_3 ( 4)

      prs = V_3 ( 5)

      bx = V_3 ( 6) + B0 ( 1)
      by = V_3 ( 7) + B0 ( 2)
      bz = V_3 ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = ( - c_30 * vx + s_30 * vy)

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( - c_30 * bx + s_30 * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R12 = AMAX1 ( S_R12, vel_long + vf )
      S_L12 = AMIN1 ( S_L12, vel_long - vf )

! ----------
! Extremize "S_L12, S_R12" using "V_2".

      rho = V_2 ( 1)

      vx = V_2 ( 2)
      vy = V_2 ( 3)
      vz = V_2 ( 4)

      prs = V_2 ( 5)

      bx = V_2 ( 6) + B0 ( 1)
      by = V_2 ( 7) + B0 ( 2)
      bz = V_2 ( 8) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = ( - c_30 * vx + s_30 * vy)

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( - c_30 * bx + s_30 * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R12 = AMAX1 ( S_R12, vel_long + vf )
      S_L12 = AMIN1 ( S_L12, vel_long - vf )

! ----------
! Extremize "S_L12, S_R12" using "0.5 * ( V_3 + V_2)".

      rho = 0.5 * ( V_3 ( 1) + V_2 ( 1))

      vx = 0.5 * ( V_3 ( 2) + V_2 ( 2))
      vy = 0.5 * ( V_3 ( 3) + V_2 ( 3))
      vz = 0.5 * ( V_3 ( 4) + V_2 ( 4))

      prs = 0.5 * ( V_3 ( 5) + V_2 ( 5))

      bx = 0.5 * ( V_3 ( 6) + V_2 ( 6)) + B0 ( 1)
      by = 0.5 * ( V_3 ( 7) + V_2 ( 7)) + B0 ( 2)
      bz = 0.5 * ( V_3 ( 8) + V_2 ( 8)) + B0 ( 3)

#if ( CHOOSE_ADIABATIC == 1)
      cs2 = gamma * prs / rho
#else
      cs2 = SQRT ( isospeed**2
     1           + gamma_euler * barotropic_stiffness
     1                         * rho**( gamma_euler - 1.0) )
#endif

      vel_long = ( - c_30 * vx + s_30 * vy)

      ob_4prho = 1.0 / ( fourpi * rho)

      va2_long = ( - c_30 * bx + s_30 * by)**2 * ob_4prho
      va2_tot = ( bx**2 + by**2 + bz**2) * ob_4prho

      b = - ( cs2 + va2_tot)
      c = cs2 * va2_long

      vf = 0.5 * ( - b + SQRT ( AMAX1 (b**2 - 4.0 * c, 0.0) ))
      vf = SQRT ( AMAX1 ( vf, 0.0) )

      S_R12 = AMAX1 ( S_R12, vel_long + vf )
      S_L12 = AMIN1 ( S_L12, vel_long - vf )

! ----------------------------------------------------------------------
! STEP 6) Obtain the resolved states and longitudinal fluxes for the three
! 1D Riemann problems that bound the strongly-interacting state.

! ------------------------------------------------------------
! 1D RP between "U_1" and "U_3" which have wavespeeds "S_L1, S_R1".
! "U_1" serves as the left state and "U_3" as the right state in the RP.
! Longitudinal direction is in the ( 0, - 1) direction.

      S1 = S_R1  ;  S2 = S_L1

      U1 = U_3  ;  U2 = U_1  ;  F1 = - G_3  ;  F2 = - G_1

      c = 1.0 / ( S1 - S2)
      a = S1 * c
      b = S2 * c

      U_Star_13 = a * U1 - b * U2 - c * ( F1 - F2)

      c = S1 * S2 * c

      F_Star_13 = a * F2 - b * F1 + c * ( U1 - U2)

! ------------------------------------------------------------
! 1D RP between "U_2" and "U_1" which have wavespeeds "S_L2, S_R2".
! "U_2" serves as the left state and "U_1" as the right state in the RP.
! Longitudinal direction is in the ( c_30, s_30) direction.

      S1 = S_R2  ;  S2 = S_L2

      U1 = U_1  ;  U2 = U_2
      F1 = c_30 * F_1 + s_30 * G_1  ;  F2 = c_30 * F_2 + s_30 * G_2

      c = 1.0 / ( S1 - S2)
      a = S1 * c
      b = S2 * c

      U_Star_21 = a * U1 - b * U2 - c * ( F1 - F2)

      c = S1 * S2 * c

      F_Star_21 = a * F2 - b * F1 + c * ( U1 - U2)

! ------------------------------------------------------------
! 1D RP between "U_3" and "U_2" which have wavespeeds "S_L12, S_R12".
! "U_3" serves as the left state and "U_2" as the right state in the RP.
! Longitudinal direction is in the ( - c_30, s_30) direction.

      S1 = S_R12  ;  S2 = S_L12

      U1 = U_2  ;  U2 = U_3
      F1 = - c_30 * F_2 + s_30 * G_2  ;  F2 = - c_30 * F_3 + s_30 * G_3

      c = 1.0 / ( S1 - S2)
      a = S1 * c
      b = S2 * c

      U_Star_32 = a * U1 - b * U2 - c * ( F1 - F2)

      c = S1 * S2 * c

      F_Star_32 = a * F2 - b * F1 + c * ( U1 - U2)

! ----------------------------------------------------------------------
! STEP 7) The RS has given us three resolved states and three longitudinal fluxes.
! We still need to construct the primitive variables from those states and use
! them to create the transverse fluxes.

! ------------------------------------------------------------
! U_Star_13 --> V_Star_13 --> G_Star_13

! --------------------
! Building "V_Star_13".

      U_Star_13 ( 1) = AMAX1 ( U_Star_13 ( 1), rhofloor)

      ob_rho  = U_Star_13 ( 1) / ( U_Star_13 ( 1)**2 + rhofloor**2)

      V_Star_13 ( 1)    = U_Star_13 ( 1)               ! density
      V_Star_13 ( 2: 4) = U_Star_13 ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_Star_13 ( 6: 8) = U_Star_13 ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_Star_13 ( NVAR + 1: ) = U_Star_13 ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_Star_13 ( 2)**2 + V_Star_13 ( 3)**2 + V_Star_13 ( 4)**2
      B2 = V_Star_13 ( 6)**2 + V_Star_13 ( 7)**2 + V_Star_13 ( 8)**2

      V_dot_B = V_Star_13 ( 2) * V_Star_13 ( 6)
     1        + V_Star_13 ( 3) * V_Star_13 ( 7)
     1        + V_Star_13 ( 4) * V_Star_13 ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_Star_13 ( 5) = ( gamma - 1.0) * ( U_Star_13 ( 5)
     1             - 0.5 * V_Star_13 ( 1) * v2  - ob_8pi * B2  )

      V_Star_13 ( 5) = AMAX1 ( V_Star_13 ( 5), prsfloor)  ! gas pressure

      U_Star_13 ( 5) = V_Star_13 ( 5) / ( gamma - 1.0)
     1         + 0.5 * V_Star_13 ( 1) * v2 + ob_8pi * B2
#else
      V_Star_13 ( 5) = isospeed**2 * U_Star_13 ( 1)
     1          + barotropic_stiffness * U_Star_13 ( 1)**gamma_euler
#endif

! --------------------
! Building "G_Star_13" which is a flux built out of "U_Star_13, V_Star_13" in
! the ( 1, 0) direction.

      rho = V_Star_13 ( 1)

      vx  = V_Star_13 ( 2)
      vy  = V_Star_13 ( 3)
      vz  = V_Star_13 ( 4)

      prs = V_Star_13 ( 5)

      Bx  = V_Star_13 ( 6)
      By  = V_Star_13 ( 7)
      Bz  = V_Star_13 ( 8)

      B0_dot_B  = B0 ( 1) * Bx + B0 ( 2) * By + B0 ( 3) * Bz

      P_tot  = prs + ob_8pi * B2 + ob_4pi * B0_dot_B

! ----------
! construct the velocity, magnetic field vectors in the (n_x, n_y)
! direction.

      n_x = 1.0  ;  n_y = 0.0

      vn  = n_x * vx + n_y * vy
      Bn  = n_x * Bx + n_y * By

      B0n = n_x * B0 ( 1) + n_y * B0 ( 2)

! ----------
! Building "G_Star_13".

      G_Star_13 ( :) = vn * U_Star_13 ( :)

      G_Star_13 ( 2) = G_Star_13 ( 2)
     1               + n_x * P_tot
     1               - ob_4pi * ( Bx * Bn + Bx * B0n + B0 ( 1) * Bn)

      G_Star_13 ( 3) = G_Star_13 ( 3)
     1               + n_y * P_tot
     1               - ob_4pi * ( By * Bn + By * B0n + B0 ( 2) * Bn)

      G_Star_13 ( 4) = G_Star_13 ( 4)
     1               - ob_4pi * ( Bz * Bn + Bz * B0n + B0 ( 3) * Bn)

#if ( CHOOSE_ADIABATIC == 1)
      G_Star_13 ( 5) = G_Star_13 ( 5)
     1               + P_tot  * vn
     1               - ob_4pi * ( Bn + B0n) * V_dot_B
#else
      G_Star_13 ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      G_Star_13 ( 6) = vn * ( B0 ( 1) + Bx) - vx * ( B0n + Bn)

      G_Star_13 ( 7) = vn * ( B0 ( 2) + By) - vy * ( B0n + Bn)

      G_Star_13 ( 8) = vn * ( B0 ( 3) + Bz) - vz * ( B0n + Bn)

! --------------------
! "F_Star_13" is in the direction "(0, -1)"; "G_Star_13"
! is in the direction "(1, 0)". We want the components in a global frame.

      F = F_Star_13
      G = G_Star_13

      F_Star_13 = G
      G_Star_13 = - F

! ------------------------------------------------------------
! U_Star_21 --> V_Star_21 --> G_Star_21

! --------------------
! Building "V_Star_21".

      U_Star_21 ( 1) = AMAX1 ( U_Star_21 ( 1), rhofloor)

      ob_rho  = U_Star_21 ( 1) / ( U_Star_21 ( 1)**2 + rhofloor**2)

      V_Star_21 ( 1)    = U_Star_21 ( 1)               ! density
      V_Star_21 ( 2: 4) = U_Star_21 ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_Star_21 ( 6: 8) = U_Star_21 ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_Star_21 ( NVAR + 1: ) = U_Star_21 ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_Star_21 ( 2)**2 + V_Star_21 ( 3)**2 + V_Star_21 ( 4)**2
      B2 = V_Star_21 ( 6)**2 + V_Star_21 ( 7)**2 + V_Star_21 ( 8)**2

      V_dot_B = V_Star_21 ( 2) * V_Star_21 ( 6)
     1        + V_Star_21 ( 3) * V_Star_21 ( 7)
     1        + V_Star_21 ( 4) * V_Star_21 ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_Star_21 ( 5) = ( gamma - 1.0) * ( U_Star_21 ( 5)
     1             - 0.5 * V_Star_21 ( 1) * v2  - ob_8pi * B2  )

      V_Star_21 ( 5) = AMAX1 ( V_Star_21 ( 5), prsfloor)  ! gas pressure

      U_Star_21 ( 5) = V_Star_21 ( 5) / ( gamma - 1.0)
     1          + 0.5 * V_Star_21 ( 1) * v2 + ob_8pi * B2
#else
      V_Star_21 ( 5) = isospeed**2 * U_Star_21 ( 1)
     1          + barotropic_stiffness * U_Star_21 ( 1)**gamma_euler
#endif

! --------------------
! Building "G_Star_21" which is a flux built out of "U_Star_21, V_Star_21" in
! the ( - s_30, c_30) direction.

      rho = V_Star_21 ( 1)

      vx  = V_Star_21 ( 2)
      vy  = V_Star_21 ( 3)
      vz  = V_Star_21 ( 4)

      prs = V_Star_21 ( 5)

      Bx  = V_Star_21 ( 6)
      By  = V_Star_21 ( 7)
      Bz  = V_Star_21 ( 8)

      B0_dot_B  = B0 ( 1) * Bx + B0 ( 2) * By + B0 ( 3) * Bz

      P_tot  = prs + ob_8pi * B2 + ob_4pi * B0_dot_B

! ----------
! construct the velocity, magnetic field vectors in the (n_x, n_y)
! direction.

      n_x = - s_30  ;  n_y = c_30

      vn  = n_x * vx + n_y * vy
      Bn  = n_x * Bx + n_y * By

      B0n = n_x * B0 ( 1) + n_y * B0 ( 2)

! ----------
! Building "G_Star_21".

      G_Star_21 ( :) = vn * U_Star_21 ( :)

      G_Star_21 ( 2) = G_Star_21 ( 2)
     1               + n_x * P_tot
     1               - ob_4pi * ( Bx * Bn + Bx * B0n + B0 ( 1) * Bn)

      G_Star_21 ( 3) = G_Star_21 ( 3)
     1               + n_y * P_tot
     1               - ob_4pi * ( By * Bn + By * B0n + B0 ( 2) * Bn)

      G_Star_21 ( 4) = G_Star_21 ( 4)
     1               - ob_4pi * ( Bz * Bn + Bz * B0n + B0 ( 3) * Bn)

#if ( CHOOSE_ADIABATIC == 1)
      G_Star_21 ( 5) = G_Star_21 ( 5)
     1               + P_tot  * vn
     1               - ob_4pi * ( Bn + B0n) * V_dot_B
#else
      G_Star_21 ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      G_Star_21 ( 6) = vn * ( B0 ( 1) + Bx) - vx * ( B0n + Bn)

      G_Star_21 ( 7) = vn * ( B0 ( 2) + By) - vy * ( B0n + Bn)

      G_Star_21 ( 8) = vn * ( B0 ( 3) + Bz) - vz * ( B0n + Bn)

! --------------------
! "F_Star_21" is in the direction "(c_30, s_30)"; "G_Star_21"
! is in the direction "(- s_30, c_30)". We want the components in a global frame.

      F = F_Star_21
      G = G_Star_21

      F_Star_21 = F * c_30 - G * s_30
      G_Star_21 = F * s_30 + G * c_30

! ------------------------------------------------------------
! U_Star_32 --> V_Star_32 --> G_Star_32

! --------------------
! Building "V_Star_32".

      U_Star_32 ( 1) = AMAX1 ( U_Star_32 ( 1), rhofloor)

      ob_rho  = U_Star_32 ( 1) / ( U_Star_32 ( 1)**2 + rhofloor**2)

      V_Star_32 ( 1)    = U_Star_32 ( 1)               ! density
      V_Star_32 ( 2: 4) = U_Star_32 ( 2: 4) * ob_rho   ! x,y,z-velocities
      V_Star_32 ( 6: 8) = U_Star_32 ( 6: 8)            ! x,y,z-B fields

#if ( NFLUID >= 1)
      V_Star_32 ( NVAR + 1: ) = U_Star_32 ( NVAR + 1: ) * ob_rho     ! species fraction
#endif

      v2 = V_Star_32 ( 2)**2 + V_Star_32 ( 3)**2 + V_Star_32 ( 4)**2
      B2 = V_Star_32 ( 6)**2 + V_Star_32 ( 7)**2 + V_Star_32 ( 8)**2

      V_dot_B = V_Star_32 ( 2) * V_Star_32 ( 6)
     1        + V_Star_32 ( 3) * V_Star_32 ( 7)
     1        + V_Star_32 ( 4) * V_Star_32 ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      V_Star_32 ( 5) = ( gamma - 1.0) * ( U_Star_32 ( 5)
     1             - 0.5 * V_Star_32 ( 1) * v2  - ob_8pi * B2  )

      V_Star_32 ( 5) = AMAX1 ( V_Star_32 ( 5), prsfloor)  ! gas pressure

      U_Star_32 ( 5) = V_Star_32 ( 5) / ( gamma - 1.0)
     1          + 0.5 * V_Star_32 ( 1) * v2 + ob_8pi * B2
#else
      V_Star_32 ( 5) = isospeed**2 * U_Star_32 ( 1)
     1          + barotropic_stiffness * U_Star_32 ( 1)**gamma_euler
#endif

! --------------------
! Building "G_Star_32" which is a flux built out of "U_Star_32, V_Star_32" in
! the ( - s_30, - c_30) direction.

      rho = V_Star_32 ( 1)

      vx  = V_Star_32 ( 2)
      vy  = V_Star_32 ( 3)
      vz  = V_Star_32 ( 4)

      prs = V_Star_32 ( 5)

      Bx  = V_Star_32 ( 6)
      By  = V_Star_32 ( 7)
      Bz  = V_Star_32 ( 8)

      B0_dot_B  = B0 ( 1) * Bx + B0 ( 2) * By + B0 ( 3) * Bz

      P_tot  = prs + ob_8pi * B2 + ob_4pi * B0_dot_B

! ----------
! construct the velocity, magnetic field vectors in the (n_x, n_y)
! direction.

      n_x = - s_30  ;  n_y = - c_30

      vn  = n_x * vx + n_y * vy
      Bn  = n_x * Bx + n_y * By

      B0n = n_x * B0 ( 1) + n_y * B0 ( 2)

! ----------
! Building "G_Star_32".

      G_Star_32 ( :) = vn * U_Star_32 ( :)

      G_Star_32 ( 2) = G_Star_32 ( 2)
     1               + n_x * P_tot
     1               - ob_4pi * ( Bx * Bn + Bx * B0n + B0 ( 1) * Bn)

      G_Star_32 ( 3) = G_Star_32 ( 3)
     1               + n_y * P_tot
     1               - ob_4pi * ( By * Bn + By * B0n + B0 ( 2) * Bn)

      G_Star_32 ( 4) = G_Star_32 ( 4)
     1               - ob_4pi * ( Bz * Bn + Bz * B0n + B0 ( 3) * Bn)

#if ( CHOOSE_ADIABATIC == 1)
      G_Star_32 ( 5) = G_Star_32 ( 5)
     1               + P_tot  * vn
     1               - ob_4pi * ( Bn + B0n) * V_dot_B
#else
      G_Star_32 ( 5) = 0.0   ! isothermal has no energy flux.
#endif

      G_Star_32 ( 6) = vn * ( B0 ( 1) + Bx) - vx * ( B0n + Bn)

      G_Star_32 ( 7) = vn * ( B0 ( 2) + By) - vy * ( B0n + Bn)

      G_Star_32 ( 8) = vn * ( B0 ( 3) + Bz) - vz * ( B0n + Bn)

! --------------------
! "F_Star_32" is in the direction "(-c_30, s_30)"; "G_Star_32"
! is in the direction "(- s_30, - c_30)". We want the components in a global frame.

      F = F_Star_32
      G = G_Star_32

      F_Star_32 = - F * c_30 - G * s_30
      G_Star_32 = F * s_30 - G * c_30

! ----------------------------------------------------------------------
! STEP 8) Now we build the three intersection points (x1, y1), (x2, y2)
! and (x3, y3) from the three 1D Riemann problems. We also check that they
! satisfy a counter-clockwise orientation.

! --------------------

      x1 = ( S_R2 + S_L1 * s_30) / c_30

      y1 = - S_L1

! --------------------

      x2 = ( S_L2 - S_R12) / ( 2.0 * c_30)

      y2 = ( S_L2 + S_R12) / ( 2.0 * s_30)

! --------------------

      x3 = ( - S_L12 - S_R1 * s_30) / c_30

      y3 = - S_R1

! --------------------
! Check for counter-clockwise orientation. Put stopping condition if orientation
! is not right.

      vx = ( x2 - x1) * ( y3 - y1) - ( y2 - y1) * ( x3 - x1)
      vy = ( x3 - x2) * ( y1 - y2) - ( y3 - y2) * ( x1 - x2)
      vz = ( x1 - x3) * ( y2 - y3) - ( y1 - y3) * ( x2 - x3)

! The orientation is right, if the "vx, vy, vz" are greater than zero.
! If not stop the execution.

      IF ( ( vx .LT. 0.0) .OR. ( vy .LT. 0.0) .OR. ( vz .LT. 0.0)) THEN

        WRITE ( 6, *)" "
        WRITE ( 6, *)"Stopping in TwoD_RS_3State ( )"
        WRITE ( 6, *)"The triangle formed by the intersection points of"
        WRITE ( 6, *)"S_L1, S_R1; S_L2, S_R2; S_L12, S_R12; does not"
        WRITE ( 6, *)"have the correct orientation"
        WRITE ( 6, *)" "

        STOP

      END IF

! ----------------------------------------------------------------------
! STEP 9) Now we go through a search for 6 conditions that identify any
! of the 6 supersonic states. Before exit, we will have to undo the frame
! rotation that was done initially.

      IF ( ( S_L1 .GE. 0.0) .AND. ( S_R2 .LE. 0.0) ) THEN

        U_Dbl_Star = U_1  ;  F_Dbl_Star = F_1  ;  G_Dbl_Star = G_1
        GO TO 11

      ELSE IF ( ( S_L2 .GE. 0.0) .AND. ( S_R12 .LE. 0.0) ) THEN

        U_Dbl_Star = U_2  ;  F_Dbl_Star = F_2  ;  G_Dbl_Star = G_2
        GO TO 11

      ELSE IF ( ( S_L12 .GE. 0.0) .AND. ( S_R1 .LE. 0.0) ) THEN

        U_Dbl_Star = U_3  ;  F_Dbl_Star = F_3  ;  G_Dbl_Star = G_3
        GO TO 11

      ELSE IF ( ( x3 * (y1 - y3) - y3 * (x1 - x3) .LE. 0.0) .AND.
     1          ( S_L1 .LE. 0.0) .AND. ( S_R1 .GE. 0.0) ) THEN

        U_Dbl_Star = U_Star_13
        F_Dbl_Star = F_Star_13  ;  G_Dbl_Star = G_Star_13
        GO TO 11

      ELSE IF ( ( x1 * (y2 - y1) - y1 * (x2 - x1) .LE. 0.0) .AND.
     1          ( S_L2 .LE. 0.0) .AND. ( S_R2 .GE. 0.0) ) THEN

        U_Dbl_Star = U_Star_21
        F_Dbl_Star = F_Star_21  ;  G_Dbl_Star = G_Star_21
        GO TO 11

      ELSE IF ( ( x2 * (y3 - y2) - y2 * (x3 - x2) .LE. 0.0) .AND.
     1          ( S_L12 .LE. 0.0) .AND. ( S_R12 .GE. 0.0) ) THEN

        U_Dbl_Star = U_Star_32
        F_Dbl_Star = F_Star_32  ;  G_Dbl_Star = G_Star_32
        GO TO 11

      END IF

! ----------------------------------------------------------------------
! STEP 10) Now build "ell31, n31_x, n31_y" for segment joining (x3, y3) to (x1, y1).
! Also build (x31, y31) as the midpoint for segment joining (x3, y3) to (x1, y1).
! "S_31" is the perpendicular speed with which this segment propagates outward from the SI state.

! Do same for "ell12, n12_x, n12_y" for segment joining (x1, y1) to (x2, y2).
! Also build (x12, y12) as the midpoint for segment joining (x1, y1) to (x2, y2).
! "S_12" is the perpendicular speed with which this segment propagates outward from the SI state.

! Do same for "ell23, n23_x, n23_y" for segment joining (x2, y2) to (x3, y3).
! Also build (x23, y23) as the midpoint for segment joining (x2, y2) to (x3, y3).
! "S_23" is the perpendicular speed with which this segment propagates outward from the SI state.

! Also build (xc, yc) the centroid of this triangle made by the SI state.
! Also build "Area" of this triangle using Heron's formula.

! --------------------

      ell31 = SQRT ( (x3 - x1)**2 + (y3 - y1)**2 )

      n31_x = (y1 - y3) / ell31  ;  n31_y = (x3 - x1) / ell31

      x31 = 0.5 * (x3 + x1)  ;  y31 = 0.5 * (y3 + y1)

      S_31 = n31_x * x3 + n31_y * y3

! --------------------

      ell12 = SQRT ( (x1 - x2)**2 + (y1 - y2)**2)

      n12_x = (y2 - y1) / ell12  ;  n12_y = (x1 - x2) / ell12

      x12 = 0.5 * (x1 + x2)  ;  y12 = 0.5 * (y1 + y2)

      S_12 = n12_x * x1 + n12_y * y1

! --------------------

      ell23 = SQRT ( (x2 - x3)**2 + (y2 - y3)**2 )

      n23_x = (y3 - y2) / ell23  ;  n23_y = (x2 - x3) / ell23

      x23 = 0.5 * (x2 + x3)  ;  y23 = 0.5 * (y2 + y3)

      S_23 = n23_x * x2 + n23_y * y2

! --------------------

      xc = (x1 + x2 + x3) / 3.0  ;  yc = (y1 + y2 + y3) / 3.0

      a = 0.5 * ( ell31 + ell12 + ell23)
      Area = SQRT ( a * (a - ell31) * (a - ell12) * (a - ell23) )

! ----------------------------------------------------------------------
! STEP 11) If we get here, the resolved state is subsonic. In that case, we build
! it here along with the associated fluxes.

! We still need to go through STEP 12 for a back rotation before exit.

! ------------------------------------------------------------
! Build states and fluxes that lie just outside the multidimensional wave model.
! The fluxes are outward-going fluxes at those boundaries.

      U1 = U_Star_13
      F1 = n31_x * F_Star_13 + n31_y * G_Star_13  ! Flux in (n31_x, n31_y) direction.

      U2 = U_Star_21
      F2 = n12_x * F_Star_21 + n12_y * G_Star_21  ! Flux in (n12_x, n12_y) direction.

      U3 = U_Star_32
      F3 = n23_x * F_Star_32 + n23_y * G_Star_32  ! Flux in (n23_x, n23_y) direction.

! ------------------------------------------------------------
! Build "U_Dbl_Star".

      a = - 0.5 / Area

      U_Dbl_Star = a * ( ( F1 - S_31 * U1) * ell31
     1                 + ( F2 - S_12 * U2) * ell12
     1                 + ( F3 - S_23 * U3) * ell23 )

! ------------------------------------------------------------
! Build "F_Dbl_Star".

      a = 1.0 / Area

      F_Dbl_Star = xc * U_Dbl_Star
     1           + a * ( ( F1 - S_31 * U1) * ell31 * ( x31 - xc)
     1                 + ( F2 - S_12 * U2) * ell12 * ( x12 - xc)
     1                 + ( F3 - S_23 * U3) * ell23 * ( x23 - xc) )

! ------------------------------------------------------------
! Build "G_Dbl_Star".

      a = 1.0 / Area

      G_Dbl_Star = yc * U_Dbl_Star
     1           + a * ( ( F1 - S_31 * U1) * ell31 * ( y31 - yc)
     1                 + ( F2 - S_12 * U2) * ell12 * ( y12 - yc)
     1                 + ( F3 - S_23 * U3) * ell23 * ( y23 - yc) )

! ----------------------------------------------------------------------
! STEP 12) Undo the initial rotation. Back rotate and exit subroutine.
! We can come to this point with a 1D
! Riemann solver solution from STEP 8, or with a 2D Riemann solver solution
! from STEP 11. Either way, the inner components ("2,3,4" and "6,7,8") of
! "U_Dbl_Star, F_Dbl_Star, G_Dbl_Star" have to be rotated back into the
! global space. In other words, we undo STEP 3.

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


      END SUBROUTINE TwoD_RS_3State



