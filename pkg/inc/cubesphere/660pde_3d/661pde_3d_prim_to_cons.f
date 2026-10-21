! ----------------------------------------------------------------------
! This is new PDE_3D_src that works.
! ----------------------------------------------------------------------

! Must set "NVAR, IPDE_SYSTEM" for this to work.








      SUBROUTINE PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, Upt, Vpt)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara, Deepak Bhoriya

! In 3D, it converts the primitive variables to Conservative variables for
! the PDE system given by "IPDE_SYSTEM".

! "gamma, Vpt ( NVAR + NFLUID) and B0 ( 3)" are inputs,
! "Upt ( NVAR + NFLUID)" is output.

! "IPDE_SYSTEM == 1" for MHD
! "IPDE_SYSTEM == 2" for CGL System
! "IPDE_SYSTEM == 3" for CGL System with p_elec

! ----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER iprs_flag

      REAL Upt ( NVAR + NFLUID), Vpt ( NVAR + NFLUID),
     1     gamma, gamma_euler, isospeed, barotropic_stiffness,
     1     molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1     ent, magB_small_cgl, eps_bdry_cgl

      REAL, DIMENSION ( 3) :: B0

! ----------------------------------------------------------------------

! Variables past this point are internal to this subroutine.

      REAL rho, prs, vlx, vly, vlz,
     1     a, b, c, V2, B2,

     1     tempaa, tempab, tempac, bfx, bfy, bfz, temp,
     1     Bsq, p_perp, p_par,
     1     updated_p_par, updated_p_perp, p_elec

      REAL, PARAMETER :: pi = 4.0 * ATAN ( 1.0),
     1                   smallnum = 1.0e-12,
     1                   ob8pi = 1.0 / ( 8.0 * pi)

! --------------------------------------------------

      iprs_flag = 0

! ----------------------------------------------------------------------
#if ( IPDE_SYSTEM == 1)
! --------------------------------------------------
! Enforce density floor, pressure floor and temperature floor, if necessary.

      Vpt ( 1) = AMAX1 ( Vpt ( 1), rhofloor)

      IF ( Vpt ( 5) <= prsfloor)   iprs_flag = 1

#if( ( CHOOSE_ADIABATIC == 1) && ( USE_COOLING_FUNCTION == 1))

      rho  = Vpt ( 1)
      prs  = Vpt ( 5)
      temp = molewt * prs / ( rho * gasconst)  ! Gas temperature.

      IF ( temp .LE. tempfloor) THEN
        Vpt ( 5) = rho * gasconst * tempfloor / molewt    ! temperature floor.
        iprs_flag = 1
      END IF

#elif( PRESSURE_POSITIVITY >= 1)

      rho  = Vpt ( 1)
      IF ( Vpt ( 5) <= prsfloor)
     1  Vpt ( 5) = AMAX1 ( prsfloor, ent * rho**( gamma - 1.0))  ! pressure from entropy

#else

      IF ( Vpt ( 5) <= prsfloor)
     1   Vpt ( 5) = AMAX1 ( prsfloor, Vpt ( 5))    ! pressure floor.

#endif

! ------------------------------

      Upt ( 1) = Vpt ( 1)

      Upt ( 2) = Vpt ( 1) * Vpt ( 2)
      Upt ( 3) = Vpt ( 1) * Vpt ( 3)
      Upt ( 4) = Vpt ( 1) * Vpt ( 4)

      Upt ( 6) = Vpt ( 6)
      Upt ( 7) = Vpt ( 7)
      Upt ( 8) = Vpt ( 8)

#if ( CHOOSE_ADIABATIC == 1)
      Upt ( 5) = Vpt ( 5) / ( gamma - 1.0)
     1   + 0.5 *  Vpt ( 1) * ( Vpt ( 2)**2 + Vpt ( 3)**2 + Vpt ( 4)**2)
     1             + ob8pi * ( Vpt ( 6)**2 + Vpt ( 7)**2 + Vpt ( 8)**2)
#else
      Upt ( 5) = 0.0
#endif


#if ( NFLUID >= 1)
      Upt ( NVAR + 1: ) = Vpt ( NVAR + 1: ) * Vpt ( 1)      ! species fraction
#endif

! --------------------------------------------------
#endif
! ----------------------------------------------------------------------
#if ( IPDE_SYSTEM == 2)
! --------------------------------------------------
! Enforce density floor and pressure floor, if necessary.

      Vpt ( 1) = AMAX1 ( Vpt ( 1), rhofloor)

! If p_par and p_perp are not in the hyperbolic limit of CGL, we force them
! to remain in the physical hypebolic limit.

      p_par  = Vpt ( 5)
      p_perp = Vpt ( 9)

      Bsq = ( Vpt ( 6) + B0 ( 1))**2
     1    + ( Vpt ( 7) + B0 ( 2))**2
     1    + ( Vpt ( 8) + B0 ( 3))**2

      CALL FORCE_Pperp_Ppar_To_HYPERBOLIC_DOMAIN
     1     ( pi, magB_small_cgl, eps_bdry_cgl, prsfloor,
     1       Bsq, p_par, p_perp,
     1       updated_p_par, updated_p_perp )

      Vpt ( 5) = updated_p_par
      Vpt ( 9) = updated_p_perp

! --------------------

      Upt ( 1) = Vpt ( 1)

      Upt ( 2) = Vpt ( 1) * Vpt ( 2)
      Upt ( 3) = Vpt ( 1) * Vpt ( 3)
      Upt ( 4) = Vpt ( 1) * Vpt ( 4)

      Upt ( 6) = Vpt ( 6)
      Upt ( 7) = Vpt ( 7)
      Upt ( 8) = Vpt ( 8)

      Upt ( 5) = 0.5 * Vpt ( 5) + Vpt ( 9)
     1   + 0.5 *  Vpt ( 1) * ( Vpt ( 2)**2 + Vpt ( 3)**2 + Vpt ( 4)**2)
     1             + ob8pi * ( Vpt ( 6)**2 + Vpt ( 7)**2 + Vpt ( 8)**2)


      Upt ( 9) = Vpt ( 5) - Vpt ( 9)

#if ( NVAR >= 10)
      Upt ( 10: ) = Vpt ( 10: ) * Vpt ( 1)      ! species fraction
#endif

! --------------------------------------------------
#endif
! End of "IPDE_SYSTEM == 2"

! ----------------------------------------------------------------------
#if ( IPDE_SYSTEM == 3)
! --------------------------------------------------
! Enforce density floor and pressure floor, if necessary.

      Vpt ( 1) = AMAX1 ( Vpt ( 1), rhofloor)

! If p_par and p_perp are not in the hyperbolic limit of CGL, we force them
! to remain in the physical hypebolic limit.

      p_par  = Vpt ( 5)
      p_perp = Vpt ( 9)
      p_elec = Vpt (10)

      Bsq = ( Vpt ( 6) + B0 ( 1))**2
     1    + ( Vpt ( 7) + B0 ( 2))**2
     1    + ( Vpt ( 8) + B0 ( 3))**2

      CALL FORCE_Pperp_Ppar_To_HYPERBOLIC_DOMAIN
     1     ( pi, magB_small_cgl, eps_bdry_cgl, prsfloor,
     1       Bsq, p_par, p_perp,
     1       updated_p_par, updated_p_perp )

      Vpt ( 5) = updated_p_par
      Vpt ( 9) = updated_p_perp

! --------------------

      Upt ( 1) = Vpt ( 1)

      Upt ( 2) = Vpt ( 1) * Vpt ( 2)
      Upt ( 3) = Vpt ( 1) * Vpt ( 3)
      Upt ( 4) = Vpt ( 1) * Vpt ( 4)

      Upt ( 6) = Vpt ( 6)
      Upt ( 7) = Vpt ( 7)
      Upt ( 8) = Vpt ( 8)

      Upt ( 5) = 0.5 * Vpt ( 5) + Vpt ( 9) + 1.5 * Vpt ( 10)
     1   + 0.5 *  Vpt ( 1) * ( Vpt ( 2)**2 + Vpt ( 3)**2 + Vpt ( 4)**2)
     1             + ob8pi * ( Vpt ( 6)**2 + Vpt ( 7)**2 + Vpt ( 8)**2)


      Upt ( 9) = Vpt ( 5) - Vpt ( 9)

      Upt (10) = Vpt ( 1) * p_elec

#if ( NVAR >= 11)
      Upt ( 11: ) = Vpt ( 11: ) * Vpt ( 1)      ! species fraction
#endif

! --------------------------------------------------
#endif
! End of "IPDE_SYSTEM == 3"
! ----------------------------------------------------------------------

      END SUBROUTINE PRIM_TO_CONS
