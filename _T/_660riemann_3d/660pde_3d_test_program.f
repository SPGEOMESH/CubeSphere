








#define CHECK_CON_PRIM 1
#define CHECK_EIGENVALUES 1
#define CHECK_FLUXES 1
#define CHECK_NCMAT 0

! Compile using:
!     cp PDE_3D_src PDE_3D_src.F
!     gfortran -fdefault-integer-8 -fdefault-real-8 -O2 PDE_3D_src.F

      PROGRAM TEST_PROGRAM


! ----------------------------------------------------------------------

      IMPLICIT NONE


      INTEGER idirection, iprim, n_phys_consts, iprs_flag

      REAL phys_consts ( 20), S_min, S_max, pi, gamma, B0 ( 3),
     1       gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       ent

      REAL, DIMENSION ( NVAR + NFLUID) :: Upt, Vpt, Fpt, g_pt, EIG

      REAL, DIMENSION ( NVAR + NFLUID, NVAR + NFLUID) ::
     1                                          A_pt, B_pt, RT, LF

! ----------------------------------------------------------------------

      INTEGER i, j, k

      CHARACTER(LEN=1), DIMENSION(3) :: dir

      dir = (/'x', 'y', 'z'/)


      pi = 4.0 * ATAN ( 1.0)

      gamma = 5.0 / 3.0 ! gamma
      B0 = 0.0

      gamma_euler = 0.0
      isospeed = 0.0
      barotropic_stiffness = 0.0

      molewt = 0.0
      gasconst = 0.0
      rhofloor = 1.0e-12
      prsfloor = 1.0e-12
      tempfloor = 1.0e-12

      iprs_flag = 1
      ent = 0.0

! ----------------------------------------------------------------------
! Check NVAR for MHD, CGL

#if  ( ( IPDE_SYSTEM == 1) && ( NVAR == 8))
#elif( ( IPDE_SYSTEM == 2) && ( NVAR == 9))
#else
      WRITE ( 6, *)"Stopping in Program: IPDE_SYSTEM, NVAR = ",
     1              IPDE_SYSTEM, NVAR
      CALL FLUSH ( 6)
      STOP
#endif

! ----------------------------------------------------------------------
#if ( IPDE_SYSTEM == 1)


      Vpt (1) = 0.33333

      Vpt (2) = 1.3
      Vpt (3) = 0.3
      Vpt (4) = 0.7

      Vpt (5) = 1.3

      Vpt (6) = 0.88
      Vpt (7) = 0.3001
      Vpt (8) = 0.987

#endif

! ----------------------------------------------------------------------
#if ( IPDE_SYSTEM == 2)

      Vpt (1) = 0.33333

      Vpt (2) = 1.3
      Vpt (3) = 0.3
      Vpt (4) = 0.7

      Vpt (5) = 0.8

      Vpt (6) = 0.88
      Vpt (7) = 0.3001
      Vpt (8) = 0.987

      Vpt (9) = 0.84

#endif
! ----------------------------------------------------------------------

      WRITE ( 6, *) "starting with prim = ... "
      WRITE ( 6, *)"[ ", Vpt , "] ;"
      WRITE ( 6, *)" "

#if ( CHECK_CON_PRIM == 1)
      CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, Upt, Vpt)
      WRITE ( 6, *) "con = ... "
      WRITE ( 6, *)"[ ", Upt , "] ;"
      WRITE ( 6, *)" "
      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, Upt, Fpt)
      WRITE ( 6, *) "prim = ... "
      WRITE ( 6, *)"[ ", Fpt , "] ;"
      WRITE ( 6, *)" "

      WRITE ( 6, *)"! =============="
      WRITE ( 6, *)"! Error from Vpt --> Upt --> Vpt :: ",
     1   SUM(ABS(Vpt-Fpt))
      WRITE ( 6, *)"! =============="

#endif
! ----------------------------------------------------------------------
#if ( CHECK_EIGENVALUES == 1)

      CALL PRIM_TO_WAVSPD ( 1,
     1         gamma, B0, Vpt, S_min, S_max)
      WRITE ( 6, *) "X-direction: Smin, Smax = ... "
      WRITE ( 6, *)"[ ",S_min, S_max, "] ;"
      WRITE ( 6, *)" "

      CALL PRIM_TO_WAVSPD ( 2,
     1         gamma, B0, Vpt, S_min, S_max)
      WRITE ( 6, *) "Y-direction: Smin, Smax = ... "
      WRITE ( 6, *)"[ ",S_min, S_max, "] ;"
      WRITE ( 6, *)" "

      CALL PRIM_TO_WAVSPD ( 3,
     1         gamma, B0, Vpt, S_min, S_max)
      WRITE ( 6, *) "Z-direction: Smin, Smax = ... "
      WRITE ( 6, *)"[ ",S_min, S_max, "] ;"
      WRITE ( 6, *)" "
#endif

! ----------------------------------------------------------------------
#if ( CHECK_FLUXES == 1)

      CALL CONS_PRIM_TO_FLUX_3D
     1           ( 1, gamma, B0, Vpt, Fpt)
      WRITE ( 6, *) "Fpt = ... "
      WRITE ( 6, *)"[ ",Fpt, "] ;"
      WRITE ( 6, *)" "

      CALL CONS_PRIM_TO_FLUX_3D
     1           ( 2, gamma, B0, Vpt, Fpt)
      WRITE ( 6, *) "g_pt = ... "
      WRITE ( 6, *)"[ ",Fpt, "] ;"
      WRITE ( 6, *)" "

      CALL CONS_PRIM_TO_FLUX_3D
     1           ( 3, gamma, B0, Vpt, Fpt)
      WRITE ( 6, *) "h_pt = ... "
      WRITE ( 6, *)"[ ",Fpt, "] ;"
      WRITE ( 6, *)" "

#endif
! ----------------------------------------------------------------------
#if ( CHECK_NCMAT == 1 )

      CALL PRIM_TO_CHARMAT_NC_3D_X
     1           ( n_phys_consts, phys_consts,
     1             Vpt, A_pt)
      WRITE ( 6, *)"ANC = ... "
      WRITE ( 6, *) "[ ..."
      DO k = 1, NVAR
         WRITE ( 6, *) A_pt(k,:), "; ..."
      END DO
      WRITE ( 6, *) "] ;"
      WRITE ( 6, *)" "

      CALL PRIM_TO_CHARMAT_NC_3D_Y
     1           ( n_phys_consts, phys_consts,
     1             Vpt, B_pt)
      WRITE ( 6, *)"BNC = ... "
      WRITE ( 6, *) "[ ..."
      DO k = 1, NVAR
         WRITE ( 6, *) B_pt(k,:), "; ..."
      END DO
      WRITE ( 6, *) "] ;"
      WRITE ( 6, *)" "

      CALL PRIM_TO_CHARMAT_NC_3D_Z
     1           ( n_phys_consts, phys_consts,
     1             Vpt, B_pt)
      WRITE ( 6, *)"CNC = ... "
      WRITE ( 6, *) "[ ..."
      DO k = 1, NVAR
         WRITE ( 6, *) B_pt(k,:), "; ..."
      END DO
      WRITE ( 6, *) "] ;"
      WRITE ( 6, *)" "

#endif
! ----------------------------------------------------------------------

      END PROGRAM TEST_PROGRAM

