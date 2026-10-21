! ----------------------------------------------------------------------
! This is new PDE_3D_src that works.
! ----------------------------------------------------------------------

! Must set "NVAR, IPDE_SYSTEM" for this to work.








      SUBROUTINE CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, Upt, Vpt)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! In 3D, it converts the conserved variables to primitive variables for
! the PDE system given by "IPDE_SYSTEM". The physical constants are
! stored in "phys_consts ( n_phys_consts)".

! "IPDE_SYSTEM == 1" for MHD
! "IPDE_SYSTEM == 2" for CGL System
! "IPDE_SYSTEM == 3" for CGL System with p_elec

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER iprs_flag

      REAL Upt ( NVAR + NFLUID), Vpt ( NVAR + NFLUID),
     1     gamma, gamma_euler, isospeed, barotropic_stiffness,
     1     molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1     ent, magB_small_cgl, eps_bdry_cgl, p_elec

      REAL, DIMENSION ( 3) :: B0

! ----------------------------------------------------------------------

! Variables past this point are internal to this subroutine.

      INTEGER i

      REAL rho, prs, vlx, vly, vlz,
     1     a, b, c, V2, B2,

     1     bfx, bfy, bfz, Bsq, temp,
     1     p_par_minus_p_perp, ener, p_perp, p_par, p_elect,

     1     updated_p_par, updated_p_perp, p_bar

      REAL, PARAMETER :: pi = 4.0 * ATAN ( 1.0),
     1                   smallnum = 1.0e-12,
     1                   ob8pi = 1.0 / ( 8.0 * pi)

! --------------------------------------------------

      iprs_flag = 0

! ----------------------------------------------------------------------
! Single fluid MHD.
#if ( IPDE_SYSTEM == 1)
! ----------------------------------------------------------------------

      Upt ( 1) = AMAX1 ( Upt ( 1), rhofloor)
      Vpt ( 1) = Upt ( 1)

      a = Vpt ( 1) / ( Vpt ( 1)**2 + rhofloor**2)

      Vpt ( 2) = Upt ( 2) * a
      Vpt ( 3) = Upt ( 3) * a
      Vpt ( 4) = Upt ( 4) * a

      Vpt ( 6) = Upt ( 6)
      Vpt ( 7) = Upt ( 7)
      Vpt ( 8) = Upt ( 8)

#if ( NFLUID >= 1)
      Vpt ( NVAR + 1: ) = Upt ( NVAR + 1: ) * a      ! species fraction
#endif


      V2 = Vpt ( 2)**2 + Vpt ( 3)**2 + Vpt ( 4)**2
      B2 = Vpt ( 6)**2 + Vpt ( 7)**2 + Vpt ( 8)**2

#if ( CHOOSE_ADIABATIC == 1)
      Vpt ( 5) = ( gamma - 1.0)
     1         * ( Upt ( 5) - 0.5 * Vpt ( 1) * V2 - ob8pi * B2)
#else
      Vpt ( 5) = isospeed**2 * Upt ( 1)
     1         + barotropic_stiffness * Upt ( 1)**gamma_euler
#endif


! ----------------------------------------
! Enforce temperature floor if we have a cooling function. If we don't have
! a cooling function, enforce pressure positivity, using entropy.

#if( ( CHOOSE_ADIABATIC == 1) && ( USE_COOLING_FUNCTION == 1))

      rho = Vpt ( 1)
      prs = Vpt ( 5)
      temp = molewt * prs / ( rho * gasconst)  ! Gas temperature.

      IF ( temp .LE. tempfloor) THEN
        Vpt ( 5) = rho * gasconst * tempfloor / molewt
        iprs_flag = 1
      END IF

#elif( PRESSURE_POSITIVITY >= 1)

      rho  = Vpt ( 1)
      IF ( Vpt ( 5) <= prsfloor) THEN
        Vpt ( 5) = AMAX1 ( prsfloor, ent * rho**( gamma - 1.0))
        iprs_flag = 1
      END IF

#else

      IF ( Vpt ( 5) <= prsfloor) THEN
        Vpt ( 5) = AMAX1 ( prsfloor, Vpt ( 5))
        iprs_flag = 1
      END IF

#endif

! ----------
! Since the "Vpt ( 5)" can be quite different from its original value;
! we have to reset the entropy variable "ent".

      ent = Vpt ( 5) / Vpt ( 1)**(gamma - 1.0)

! ----------
! Update Upt.

#if ( CHOOSE_ADIABATIC == 1)
      Upt ( 5) = Vpt ( 5) / ( gamma - 1.0) + 0.5 *  Vpt ( 1) * V2
     1         + ob8pi * B2
#else
      Upt ( 5) = 0.0
#endif

! ----------------------------------------------------------------------
#endif
! End of "IPDE_SYSTEM == 1"
! ----------------------------------------------------------------------
! Single fluid CGL.
#if ( IPDE_SYSTEM == 2)
! ----------------------------------------------------------------------


      Upt ( 1) = AMAX1 ( Upt ( 1), rhofloor)
      Vpt ( 1) = Upt ( 1)

      a = Vpt ( 1) / ( Vpt ( 1)**2 + rhofloor**2)

      Vpt ( 2) = Upt ( 2) * a
      Vpt ( 3) = Upt ( 3) * a
      Vpt ( 4) = Upt ( 4) * a

      Vpt ( 6) = Upt ( 6)
      Vpt ( 7) = Upt ( 7)
      Vpt ( 8) = Upt ( 8)

#if ( NVAR >= 10)

      Vpt ( 10: ) = AMAX1 ( 0.0, AMIN1 ( 1.0,
     1                           Upt ( 10: ) * a ) )  ! tracers must be between 0 and 1.

      b = 0.0
      DO i = 10, NVAR
      b = b + Vpt ( i)
      END DO

      b = 1.0 / b
      Vpt ( 10: ) = Vpt ( 10: ) * b

      Upt ( 10: ) = Vpt ( 1) * Vpt ( 10: )

#endif

      V2 = Vpt ( 2)**2 + Vpt ( 3)**2 + Vpt ( 4)**2
      B2 = Vpt ( 6)**2 + Vpt ( 7)**2 + Vpt ( 8)**2

! ----------
! Obtain P∥ using, energy U(5) and the pressure difference U(9)=P∥-P⊥.

      p_bar = ( 2.0 / 3.0) * ( Upt ( 5) - 0.5*Vpt (1)*V2 - ob8pi * B2)

      p_perp = p_bar - Upt ( 9) / 3.0
      p_par = 3.0 * p_bar - 2.0 * p_perp

#if( PRESSURE_POSITIVITY >= 1)
      IF ( (p_par <= prsfloor) .OR. (p_perp <= prsfloor) ) THEN

        p_bar = (ent * Vpt ( 1)**4)**(1.0/3.0)
        p_perp = p_bar  ;  p_par = p_bar

        iprs_flag = 1

      END IF
#else
      IF ( (p_par <= prsfloor) .OR. (p_perp <= prsfloor) ) THEN

        p_par = AMAX1 ( p_par, prsfloor)
        p_perp = AMAX1 ( p_perp, prsfloor)
        p_bar = ( p_par + 2.0 * p_perp) / 3.0

        iprs_flag = 1

      END IF
#endif


! ----------------------------------------
! If p_par and p_perp are not in the hyperbolic limit of CGL, we force them
! to remain in the physical hypebolic limit.

      Bsq = ( Vpt ( 6) + B0 ( 1))**2
     1    + ( Vpt ( 7) + B0 ( 2))**2
     1    + ( Vpt ( 8) + B0 ( 3))**2

      CALL FORCE_Pperp_Ppar_To_HYPERBOLIC_DOMAIN
     1     ( pi, magB_small_cgl, eps_bdry_cgl, prsfloor,
     1       Bsq, p_par, p_perp,
     1       updated_p_par, updated_p_perp )

      Vpt ( 5) = updated_p_par
      Vpt ( 9) = updated_p_perp

! ----------
! Since the "updated_p_par, updated_p_perp" can be quite different from "p_par, p_perp";
! we have to reset the entropy variable "ent".

      ent = updated_p_par * updated_p_perp**2 / Vpt ( 1)**4

! ----------
! Update Upt

      Upt ( 5) = 0.5 * Vpt ( 5) + Vpt ( 9)
     1   + 0.5 *  Vpt ( 1) * ( Vpt ( 2)**2 + Vpt ( 3)**2 + Vpt ( 4)**2)
     1             + ob8pi * ( Vpt ( 6)**2 + Vpt ( 7)**2 + Vpt ( 8)**2)


      Upt ( 9) = Vpt ( 5) - Vpt ( 9)


! ----------------------------------------------------------------------
#endif
! End of "IPDE_SYSTEM == 2".
! ----------------------------------------------------------------------
! Two fluids, ions + electron pressure, with CGL.
#if ( IPDE_SYSTEM == 3)
! ----------------------------------------------------------------------

      Upt ( 1) = AMAX1 ( Upt ( 1), rhofloor)
      Vpt ( 1) = Upt ( 1)

      a = Vpt ( 1) / ( Vpt ( 1)**2 + rhofloor**2)

      Vpt ( 2) = Upt ( 2) * a
      Vpt ( 3) = Upt ( 3) * a
      Vpt ( 4) = Upt ( 4) * a

      Vpt ( 6) = Upt ( 6)
      Vpt ( 7) = Upt ( 7)
      Vpt ( 8) = Upt ( 8)

      Vpt (10) = Upt (10) * a
      p_elect = Vpt (10)

#if ( NVAR >= 11)

      Vpt ( 11: ) = AMAX1 ( 0.0, AMIN1 ( 1.0,
     1                           Upt ( 11: ) * a ) )  ! tracers must be between 0 and 1.

      b = 0.0
      DO i = 11, NVAR
      b = b + Vpt ( i)
      END DO

      b = 1.0 / b
      Vpt ( 11: ) = Vpt ( 11: ) * b

      Upt ( 11: ) = Vpt ( 1) * Vpt ( 11: )

#endif


      V2 = Vpt ( 2)**2 + Vpt ( 3)**2 + Vpt ( 4)**2
      B2 = Vpt ( 6)**2 + Vpt ( 7)**2 + Vpt ( 8)**2

! ----------
! Obtain P∥ using, energy U(5) and the pressure difference U(9)=P∥-P⊥.

      p_bar = Upt ( 5) - 0.5*Vpt (1)*V2 - ob8pi * B2
      p_bar = ( 2.0 / 3.0) * p_bar - p_elect

      p_perp = p_bar - Upt ( 9) / 3.0
      p_par = 3.0 * p_bar - 2.0 * p_perp

#if( PRESSURE_POSITIVITY >= 1)
      IF ( (p_par <= prsfloor) .OR. (p_perp <= prsfloor) .OR.
     1     (p_elect <= prsfloor) ) THEN

        p_bar = (ent * Vpt ( 1)**4)**(1.0/3.0)
        p_perp = p_bar  ;  p_par = p_bar

        p_elect = AMAX1 ( p_elect, prsfloor)

        iprs_flag = 1

      END IF
#else
      IF ( (p_par <= prsfloor) .OR. (p_perp <= prsfloor) .OR.
     1     (p_elect <= prsfloor) ) THEN

        p_par = AMAX1 ( p_par, prsfloor)
        p_perp = AMAX1 ( p_perp, prsfloor)
        p_bar = ( p_par + 2.0 * p_perp) / 3.0

        p_elect = AMAX1 ( p_elect, prsfloor)

        iprs_flag = 1

      END IF
#endif

! ----------------------------------------
! If p_par and p_perp are not in the hyperbolic limit of CGL, we force them
! to remain in the physical hypebolic limit.

      Bsq = ( Vpt ( 6) + B0 ( 1))**2
     1    + ( Vpt ( 7) + B0 ( 2))**2
     1    + ( Vpt ( 8) + B0 ( 3))**2

      CALL FORCE_Pperp_Ppar_To_HYPERBOLIC_DOMAIN
     1     ( pi, magB_small_cgl, eps_bdry_cgl, prsfloor,
     1       Bsq, p_par, p_perp,
     1       updated_p_par, updated_p_perp )

      Vpt ( 5) = updated_p_par
      Vpt ( 9) = updated_p_perp
      Vpt ( 10) = p_elect

! ----------
! Since the "updated_p_par, updated_p_perp" can be quite different from "p_par, p_perp";
! we have to reset the entropy variable "ent".

      ent = updated_p_par * updated_p_perp**2 / Vpt ( 1)**4

! ----------
! Update Upt

      Upt ( 5) = 0.5 * Vpt ( 5) + Vpt ( 9) + 1.5 * Vpt ( 10)
     1   + 0.5 *  Vpt ( 1) * ( Vpt ( 2)**2 + Vpt ( 3)**2 + Vpt ( 4)**2)
     1             + ob8pi * ( Vpt ( 6)**2 + Vpt ( 7)**2 + Vpt ( 8)**2)


      Upt ( 9) = Vpt ( 5) - Vpt ( 9)

      Upt ( 10) = Vpt ( 1) * Vpt ( 10)

! ----------------------------------------------------------------------
#endif
! End of "IPDE_SYSTEM == 3".
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END SUBROUTINE CONS_TO_PRIM




