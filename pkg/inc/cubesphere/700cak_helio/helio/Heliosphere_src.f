




#define DEBUG_HELIO_INIT 0
     
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      MODULE HELIOSPHERE_MOD

! ----------------------------------------------------------------------

      IMPLICIT NONE
      REAL, PARAMETER :: M_PI = 3.141592653589793

      INTEGER, PARAMETER :: num_silo_images = 300      ! number of silo image files

      REAL, PARAMETER :: R_inner =   30.0,             ! inner radius in au
     1                   R_outer = 1000.0,             ! outer radius in au
     1                   stop_time = 365.2425 * 150.0   ! x days

! ----------------------------------------------------------------------

      REAL, PARAMETER :: Msun = 1.98847542e+33,  Gcgs = 6.6743e-08,
     1                   Lsun = 3.828e+33,       ccgs = 2.99792458e+10,
     1                   amu  = 1.66053904e-24,  kB   = 1.38064852e-16,

     1      one_au     = 1.49598E+13,                ! 1 AU (cm)
     1      one_day    = 24.0 * 60.0 * 60.0,         ! 1 day in (s) 
     1      one_year   = 365.2425 * one_day,         ! 1 year in (s)

! ----------

     1      mu          = 0.5,               ! mean molecular weight
     1      Rgas        = kB / amu ,               ! gas constant is chosen to be unity
     1      solar_T     = 1.0E6,             ! Temperature of solar wind gas at 30au in (K)

     1      solar_Mdot  = 2.5e-14 * Msun / one_year,          ! in cgs
     1      solar_Vsnd  = SQRT( kB * solar_T / (mu * amu) ),  ! Sound speed

     1      solar_uc    = 3.5 * solar_Vsnd,   ! Mach Number * sound speed at 'R_inner'

     1      solar_rho   = solar_Mdot 
     1          / ( 4.0 * M_PI * ((R_inner * one_au)**2) * solar_uc),

     1      solar_prs   =  solar_rho * Rgas * solar_T /  mu ,    ! Prs at inner boundary.

     1      solar_Omega = 0.0       ! Angular velocity of Sun's rotation.

      REAL, PARAMETER :: 
     1     unit_dist     = one_au,
     1     unit_time     = one_day,
     1     unit_mass     = solar_rho * one_au**3

! ------------------------------

      REAL, PARAMETER :: code_density = solar_rho
     1                        / ( unit_mass / unit_dist**3),
     1                   code_Omega   = solar_Omega * unit_time,

     1                   code_prs     = solar_prs
     1                   / ( unit_mass / (unit_dist * unit_time**2) ),

     1                  r_lism = 150.0 *
     1                            one_au / unit_dist,
     1                  r_1 = one_au / unit_dist,
     1                  B0_sw  =  50.0e-6 / SQRT ( unit_mass 
     1                  / (unit_dist * unit_time**2) )

! ----------------------------------------------------------------------

      END MODULE HELIOSPHERE_MOD

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE HELIO_INITIALIZE 

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,  

     1                    U_pt, A_pt )

! ----------------------------------------------------------------------
! Written By: Dinshaw Balsara, Vladimir Florinski, Sethu and Swati Sharma

! When CGL is used, "V_pt ( 5)" is pressure in the direction parallel to the field.
! "V_pt ( 9)" is the pressure perpendicular to the field.

! When dealing with conserved variables, "U_pt ( 5)" is the total energy.
! "U_pt ( 9)" is "P_parallel - P_perpendicular".

! IMPORTANT CHANGES:-

! 1) This problem instantaneously generates a Mach 8+ shock at the interface
! between the high speed, outgoing, solar wind and the LISM. This produces
! very rapid and extreme rise in pressure at the interface, which is bad for
! the timestep. For that reason, start with an initial CFL of 0.001 before
! "INIT_PATCH_CUBE_SPHERE".
! The line that needs to be changed is:- "cfl_coef = AMIN1 ( 0.0001, a)".

! 2) Furthermore, in "EVALUATE_TIMESTEP" we add the lines:-
! #if( CHOOSE_APPL_HELIO == 1)
!       dtcurp1 = 1.01 * dtcur
! #else
!       dtcurp1 = 1.1 * dtcur
! #endif
! This ensures that we don't have a sudden increase in the timestep.

! 3) Set inner boundary to to 17 so that the initial inner ghost zone values are
! untouched in the course of time-evolution. Set outer boundary to just reimpose
! it at every timestep so that any time-dependent initialization of magnetic fields
! at the outer boundary is smoothly convected into the interior. To achieve this,
! before "CALL PAD_BOUNDARY_R" set 
! "IF ( iz_patch .EQ. 1         )  bcarr_r ( :, :, 1) = 17" for inner radial boundary.
! also for the outer boundary, before "CALL ADD_BOUNDARY" set
! "IF ( iz_patch .EQ. nr_patches)  ibdy_outer = 1" for outer radial boundary.

! 4) The problem also requires a steady inflow at "rmin". For that reason,
! some extra lines have been added in the limiter "LIMIT_ALL_SLOPES_P1"
! so that the first ghost zone inside of "rmin" will have its slope zeroed out.
! See the part of the code with ifdef "#if ( CHOOSE_APPL_HELIO == 1)" right
! above STEP 8 in "LIMIT_ALL_SLOPES_P1".

! This important step is before STEP 6 in "LIMIT_ALL_SLOPES_P1" we put the
! following code that provides a floor to the pressure in the first layer of
! active zones at the inner radial boundary. This is done as follows:-
! ----------------------------------------------------------------------
! Remove this
#if (0 == 1)
! ----------------------------------------------------------------------

#if ( CHOOSE_APPL_HELIO == 1)

      IF ( iz_patch .EQ. 1 ) THEN

! To prevent pressure in the first radial zone from sinking for the     
! heliosphere problem, we provide a floor.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost,                &
!$OMP&    V, U, gamma, imask_active_patch)                              &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)   
!$acc parallel vector_length(64) present (                              &
!$acc&    V, U, gamma, imask_active_patch)
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

      izz = 1

      IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! --------------------------------------------------
! This value of "3.8914E-5" is actually "prs_sw" evaluated at the inner boundary
! of the heliosphere boundary. For CGL, we will have to fix two such pressures.

#if (IPDE_SYSTEM == 1)

! For MHD, we have a single pressure in the 5th component of "V" and a single
! total Energy in "U".

!        IF ( V ( ixx, iyy, izz, 5, 1) .LT. 3.8914E-5) THEN

        U ( ixx, iyy, izz, 5, 1) = U ( ixx, iyy, izz, 5, 1)
     1   + (3.8914E-5 - V ( ixx, iyy, izz, 5, 1)) / (gamma - 1.0)

        V ( ixx, iyy, izz, 5, 1) = 3.8914E-5

!        END IF

#elif (IPDE_SYSTEM == 2)

! For CGL we have "V (5) = P_parallel ; V (9) = P_perpendicular" and
! "U (5) = Total energy density ; U (9) = P_parallel - P_perpendicular".
! In the setup we want "P_parallel == P_perpendicular" exactly.

        V ( ixx, iyy, izz, 5, 1) = 3.8914E-5

        V ( ixx, iyy, izz, 9, 1) = 3.8914E-5

        U ( ixx, iyy, izz, 5, 1) = 1.5 * 3.8914E-5
     1+ 0.5 * V ( ixx, iyy, izz, 1, 1) * ( V ( ixx, iyy, izz, 2, 1)**2
     1  + V ( ixx, iyy, izz, 3, 1)**2 + V ( ixx, iyy, izz, 4, 1)**2 )
     1 + ( V ( ixx, iyy, izz, 6, 1)**2 + V ( ixx, iyy, izz, 7, 1)**2
     1   + V ( ixx, iyy, izz, 8, 1)**2 ) * ob_8pi

        U ( ixx, iyy, izz, 9, 1) = 0.0  ! Pressures are isotropic.

#endif

! --------------------------------------------------

      END IF    ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN".

      END DO    ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".     
      END DO    ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".     
!$acc end parallel
!$OMP END PARALLEL DO

      END IF    ! End of "IF ( iz_patch .EQ. 1 ) THEN".

#endif
! ----------------------------------------------------------------------
#endif  
! Remove this
! ----------------------------------------------------------------------

! 5) The other important step is before STEP 8 in "LIMIT_ALL_SLOPES_P1"
! we zero out the slopes in the ghost zone and the first layer of dynamically
! active zones at the innermost radial boundary. This is so that the wind
! conditions can be held stationary.
! ----------------------------------------------------------------------
! Remove this
#if (0 == 1)
! ----------------------------------------------------------------------
#if ( CHOOSE_APPL_HELIO == 1)
! ----------------------------------------------------------------------

! For Wind type problems, we have a density and pressure gradient in the
! density and pressure profiles. As a result, any time-variation in the
! solution at the first inner boundary of ghost zones can throw off the
! calculation.

! For that reason:-

! Set slopes of inner boundary ghost zones to zero, for the Heliosphere
! simulations. This is to prevent the pressure from touching floor at the
! inner boundary. For CGL, remember to zero out gradients in both pressures.
! ------------------------------------------------------------

      IF ( iz_patch .EQ. 1 ) THEN         ! For the inner most patch.

! ------------------------------------------------------------
!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz)                              &
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost,                &
!$OMP&    V, B1_zc, B2_zc, B3_zc, imask_active_patch)                   &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)
!$acc parallel vector_length(64) present (                              &
!$acc&    V, B1_zc, B2_zc, B3_zc, imask_active_patch)
!$acc loop gang vector collapse(3) independent                          &
!$acc&    private ( ixx, iyy, izz)

      DO izz = 1 - NGhost, 1
      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

      IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN

! --------------------------------------------------

        V ( ixx, iyy, izz, :, 2: ) = 0.0

        B1_zc ( ixx, iyy, izz, 2:) = 0.0
        B2_zc ( ixx, iyy, izz, 2:) = 0.0
        B3_zc ( ixx, iyy, izz, 2:) = 0.0

! --------------------------------------------------

      END IF    ! End of "IF ( imask_active_patch ( ixx, iyy, izz) .EQ. 1) THEN".

      END DO    ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO    ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
      END DO    ! End of "DO izz = 1 - NGhost, 1".

!$acc end parallel
!$OMP END PARALLEL DO
! ------------------------------------------------------------

      END IF    ! End of "IF ( iz_patch .EQ. 1 ) THEN".

! ------------------------------------------------------------
#endif
! End of "CHOOSE_APPL_HELIO == 1".
! ----------------------------------------------------------------------
#endif  
! Remove this
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
!$acc routine seq

      USE HELIOSPHERE_MOD

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      REAL :: x, y, z, time, rmin, rmax, zone_size,
     1        xcen_sector, ycen_sector, zcen_sector,
     1        xvert1_sector, yvert1_sector, zvert1_sector

      INTEGER :: n_cc_components, iproblem

      REAL, DIMENSION ( n_cc_components) :: U_pt, V_pt

      REAL :: A_pt ( 3)

! ----------------------------------------------------------------------
! local variables.

      INTEGER  bflag, iprs_flag

      REAL  rh, r, time_use, T_lism,

     1     rho_lism, prs_lism, Ax, Ay, Az, Bx, By, Bz, ux_lism, uy_lism,

     1     uz_lism, ent, uc_sw, rho_sw, prs_sw,B0 (3), a,

     1     ratio, dr_lism

      REAL, DIMENSION ( n_cc_components) :: V_parker, V_lism

! ----------------------------------------------------------------------
! Parker Solar wind problem.


      rh = SQRT ( x**2 + y**2)
      r  = SQRT ( x**2 + y**2 + z**2)

! ----------------------------------------------------------------------
! Setup the Parker spiral solar wind "sw" parameters.
! Physical density and pressure of SW at 1 AU.

      rho_sw =  (amu * 5.0) * (r_1/r)**2 
      prs_sw = (5.0*amu)* Rgas *(rho_sw/(amu* 5.0))**gamma  
     1                                        * solar_T / mu 

! ----------

      rho_sw = rho_sw  / (unit_mass / unit_dist**3)           ! In code units.
      prs_sw = prs_sw 
     1       / (unit_mass / (unit_dist * unit_time**2) )      ! In code units.

      uc_sw  = solar_uc / ( unit_dist / unit_time)            ! In code units.

! ------------------------------
! Initialize fluid parameters for the parker spiral wind.

      V_parker ( 1) = rho_sw
      V_parker ( 5) = prs_sw

#if ( IPDE_SYSTEM == 2)
      V_parker ( 9) = prs_sw
#endif

      V_parker ( 2) = uc_sw * x / r
      V_parker ( 3) = uc_sw * y / r
      V_parker ( 4) = uc_Sw * z / r

! ----------------------------------------------------------------------
! Setup the LISM inflow fluid parameters.
 
      rho_lism =  amu * 0.11                     ! Physical density of LISM.
      prs_lism = Rgas * rho_lism * 7000.0 / mu   ! Physical pressure of LISM.

      rho_lism = rho_lism / (unit_mass / unit_dist**3)        ! in code units.
      prs_lism = ( prs_lism)
     1           / (unit_mass / (unit_dist * unit_time**2) )  ! in code units


      ux_lism  = -26 * 1000.0 * 100.0 / ( unit_dist / unit_time)
      uy_lism  = 0.0
      uz_lism  = 0.0

! ----------
! Build a global unit vector that gives us the direction in which LISM velocity
! is oriented.

      a = 1.0 / SQRT ( ux_lism**2 + uy_lism**2 + uz_lism**2)

      inflow_uvec ( 1) = ux_lism * a
      inflow_uvec ( 2) = uy_lism * a
      inflow_uvec ( 3) = uz_lism * a

! ------------------------------
! Initialize fluid parameters for the LISM.

      V_lism ( 1) = rho_lism
      V_lism ( 5) = prs_lism

#if ( IPDE_SYSTEM == 2)
      V_lism ( 9) = prs_lism
#endif
! ----------

      V_lism ( 2) = (1.0 + r_lism**3 / (2.0 * r**3)) * ux_lism 
     1            - (3.0 * r_lism**3 / (2.0 * r**5))
     1            * (ux_lism * x + uy_lism * y + uz_lism * z ) * x

      V_lism ( 3) = (1.0 + r_lism**3 / (2.0 * r**3)) * uy_lism
     1            - (3.0 * r_lism**3 / (2.0 * r**5))
     1            * (ux_lism * x + uy_lism * y + uz_lism * z ) * y

      V_lism ( 4) = (1.0 + r_lism**3 / (2.0 * r**3)) * uz_lism
     1            - (3.0 * r_lism**3 / (2.0 * r**5))
     1            * (ux_lism * x + uy_lism * y + uz_lism * z ) * z

! ----------------------------------------------------------------------
! Setup the problem with a smooth TANH transition between Parker wind and LISM.

      dr_lism = 0.1 * r_lism

      ratio = ( r - r_lism) / dr_lism

      a = 0.5 * ( TANH ( ratio) + 1.0)

      IF ( ratio .GT.  10.0) a = 1.0
      IF ( ratio .LT. -10.0) a = 0.0

      V_pt ( 1:5) = a * V_lism ( 1:5) + ( 1.0 - a) * V_parker ( 1:5)

#if ( IPDE_SYSTEM == 2)
      V_pt ( 9)   = a * V_lism ( 9)   + ( 1.0 - a) * V_parker ( 9)
#endif

! ----------------------------------------------------------------------
! Obtain magnetic field and vector potential

      time_use = 0.0

! ----------

      bflag = 1      ! Magnetic field

      CALL GET_MAGNETIC_FIELD_HELIO ( bflag, time_use, x, y, z, 
     1                                                 Bx, By, Bz)

      V_pt ( 6) = Bx
      V_pt ( 7) = By
      V_pt ( 8) = Bz

! ----------

      bflag = 0      ! Vector potential

      CALL GET_MAGNETIC_FIELD_HELIO ( bflag, time_use, x, y, z, 
     1                                                 Ax, Ay, Az)

      A_pt ( 1) = Ax
      A_pt ( 2) = Ay
      A_pt ( 3) = Az

! ----------------------------------------------------------------------
! Obtain the conserved variables from the primitive variables.

        B0 ( :) = 0.0         ! Background magnetic field is zero.

! For MHD simulations, if necessary, the pressure is reset using the entropy.
        ent = V_pt ( 5) / V_pt ( 1)**( gamma - 1.0) 

        CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl, 
     1       iprs_flag, ent, B0, U_pt, V_pt)

! ----------------------------------------------------------------------

      END SUBROUTINE HELIO_INITIALIZE

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE GET_MAGNETIC_FIELD_HELIO ( bflag, time,
     1                                      x, y, z, ABx, ABy, ABz )

! ----------------------------------------------------------------------
! Written By: Dinshaw Balsara, Vladimir Florinski, Sethu and Swati Sharma

! This subroutine contains the magnetic field, that is necessary for
! Heliosphere simulations.

! ----------------------------------------------------------------------
!$acc routine seq

      USE HELIOSPHERE_MOD

      IMPLICIT NONE

      INTEGER bflag

      REAL time, x, y, z, ABx, ABy, ABz, Bx, By, Bz, nx, ny, nz

      REAL  Bxr1, Byr1, Bzr1, Bxr, Byr, Bzr, alpha1 , alpha2,
     1                     wx, wy, wz

      REAL  Dx, Dy, Dz, Cx, Cy, Cz

! ----------------------------------------------------------------------
! local variables.

      INTEGER i, j, k, ii

      REAL r, t, p, rh, b, tht, phi,
     1     Sr, St, Sp, xx

! ----------------------------------------------------------------------

      rh = SQRT ( x**2 + y**2)
      r  = SQRT ( x**2 + y**2 + z**2)

      phi = ATAN2 (  y, x)
      tht = ATAN2 ( rh, z)

! ----------------------------------------------------------------------
! Parker Spiral.


! Solar angluar rotation unit vector w

      wx = 0.0
      wy = 0.0
      wz = 1.0

      CALL CROSS_PRODUCT (  wx, wy, wz, x, y, z,
     1                                       Cx, Cy, Cz )


      CALL CROSS_PRODUCT (  Cx, Cy, Cz, x, y, z,
     1                                       Dx, Dy, Dz )


! ----------------------------------------------------------------------

      IF ( bflag == 1) THEN ! Compute Magnetic field

! ----------------------------------------------------------------------
! Parker spiral

      ABx =SQRT ( 2.0) * B0_sw * ( r_1 / r**2)* Cx
      ABy= SQRT ( 2.0) * B0_sw * ( r_1 / r**2)* Cy
      ABz= SQRT ( 2.0) * B0_sw * ( r_1 / r**2)* Cz

! ----------------------------------------------------------------------

      ELSE     ! Compute Magnetic Vector Potential

! ----------------------------------------------------------------------
! Parker spiral

      ABx= SQRT ( 2.0) * B0_sw * ( r_1 / r) * x +
     1           SQRT ( 2.0) * (B0_sw / 3.0) * ( r_1 / r**2)* Dx

      ABy= SQRT ( 2.0) * B0_sw * ( r_1 / r) * y +
     1           SQRT ( 2.0) * B0_sw * ( r_1 / r**2)* Dy


       ABz = SQRT ( 2.0) * B0_sw * ( r_1 / r) * z +
     1           SQRT ( 2.0) * B0_sw * ( r_1 / r**2)* Dz

! ----------------------------------------------------------------------

      END IF ! End of "IF ( bflag == 1) THEN".

! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! VLISM boundary conditions.

      IF ( r > r_lism) THEN


      Bx =  3.0e-6/ SQRT ( unit_mass / (unit_dist * unit_time**2) )  
      By =  0.0/ SQRT ( unit_mass / (unit_dist * unit_time**2) )
      Bz =  0.0 / SQRT ( unit_mass / (unit_dist * unit_time**2) )

      nx = 0.0 
      ny = -1.0
      nz = 0.0


      alpha1=140.0 * M_PI /180.0

        CALL ROTATE_VECTOR ( bflag, nx, ny, nz, Bx, By, Bz, alpha1,
     1                                       Bxr1, Byr1, Bzr1 )
      
      nx = 1.0
      ny = 0.0
      nz = 0.0

      alpha2=45.0 * M_PI /180.0


       CALL ROTATE_VECTOR ( bflag, nx, ny, nz, Bxr1, Byr1, Bzr1, alpha2,
     1                                       Bxr, Byr, Bzr )



        IF ( bflag == 1) THEN   ! Magnetic field

        ABx = (1.0 + r_lism**3 / (2.0 * r**3)) * Bxr
     1            - (3.0 * r_lism**3 / (2.0 * r**5))
     1            * (Bxr * x + Byr * y + Bzr * z ) * x

        ABy = (1.0 + r_lism**3 / (2.0 * r**3)) * Byr
     1            - (3.0 * r_lism**3 / (2.0 * r**5))
     1            * (Bxr * x + Byr * y + Bzr * z ) * y

        ABz= (1.0 + r_lism**3 / (2.0 * r**3)) * Bzr
     1            - (3.0 * r_lism**3 / (2.0 * r**5))
     1            * (Bxr * x + Byr * y + Bzr * z ) * z

! ----------------------------------------------------------------------

      ELSE     ! Compute Magnetic Vector Potential

! ----------------------------------------------------------------------

      CALL CROSS_PRODUCT (  Bxr, Byr, Bzr, x, y, z,
     1                                       Dx, Dy, Dz )


          ABx = (0.5) * (1 - (r_lism/r)**3) * Dx
          ABy = (0.5) * (1 - (r_lism/r)**3) * Dy
          ABz = (0.5) * (1 - (r_lism/r)**3) * Dz

         END IF

        END IF

! ----------------------------------------------------------------------

      END SUBROUTINE GET_MAGNETIC_FIELD_HELIO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE GET_HELIO_PARAMETERS (
     1     iproblem, nsilo_images,
     1     rmin, rmax, gamma, gasconst, molewt, timestop,
     1     omega_passive, isospeed, tempfloor,
     1     tau_phys_cgl, f_lowbd_cgl,
     1     neta_f_cgl, kappa_f_cgl, magB_small_cgl, eps_bdry_cgl) 

! ----------------------------------------------------------------------
! Written By: Dinshaw Balsara, Vladimir Florinski, Sethu and Swati Sharma

! subroutines returns the values of constants for initialization.

! ----------------------------------------------------------------------
!$acc routine seq

      USE HELIOSPHERE_MOD

      IMPLICIT NONE

      INTEGER iproblem, nsilo_images

      REAL rmin, rmax, gamma, gasconst, molewt, timestop,
     1     omega_passive ( 3), isospeed, tempfloor,
     1     tau_phys_cgl, f_lowbd_cgl,
     1     neta_f_cgl, kappa_f_cgl, magB_small_cgl, eps_bdry_cgl

! ----------------------------------------------------------------------
      iproblem = 14   ! Heliosphere problem is entered as 14th problem in Geomesh code.
                                                                         
      nsilo_images = num_silo_images                                     
                                                                         
      rmin = R_inner                                                     
      rmax = R_outer

! ---------------------------------------------------------------------- 
                                                                         
      gamma    = 1.67                                            
      gasconst = 1.0                                                     
      molewt   = mu                                                      
                                                                         
      timestop = stop_time

#if( IPDE_SYSTEM == 2)
      tau_phys_cgl = 30.0 * one_year / unit_time  ! CGL relaxation time - 30 years

      f_lowbd_cgl  = 1.0e-5

      neta_f_cgl     = 4.0        ! neta for f_beta
      kappa_f_cgl    = 4.0        ! kappa for f_beta
      magB_small_cgl = 1.0e-10    ! set p_par = p_perp if magB < magB_small

      eps_bdry_cgl   = 0.025    ! p_par and p_perp should be epsilon distance away 
                                ! from mirror and firehose limits. Take 0.01 or 0.025.
#endif

! ------------------------------                                         
                                                                         
      omega_passive ( :) = 0.0                                           
      omega_passive ( 3) = code_Omega         ! Angular velocity              
                                                                         
      isospeed = solar_Vsnd / ( unit_dist / unit_time)

! ----------------------------------------------------------------------

      END SUBROUTINE GET_HELIO_PARAMETERS

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE ROTATE_VECTOR ( bflag, nx, ny, nz, Bx, By, Bz, alpha,
     1                                       Bxr, Byr, Bzr )

! ----------------------------------------------------------------------
! Written By: Dinshaw Balsara, Vladimir Florinski, Sethu and Swati Sharma

! This subroutine rotates the vector by an angle alpha

! ----------------------------------------------------------------------
!$acc routine seq

      USE HELIOSPHERE_MOD

      IMPLICIT NONE

      INTEGER bflag

      REAL time, nx, ny, nz, Bxr, Byr, Bzr, Bx, By, Bz, alpha

      REAL Cx , Cy , Cz

! ----------------------------------------------------------------------
! local variables.

      REAL  dot_prod 

! ----------------------------------------------------------------------
! cross product n cross B = C

      CALL CROSS_PRODUCT (  nx, ny, nz, Bx, By, Bz,
     1                                       Cx, Cy, Cz )

!dot product n.B
      dot_prod= nx * Bx + ny* By + nz * Bz

!Rotated vector

      Bxr= COS(alpha) * Bx + (1- COS(alpha))* dot_prod * nx
     1                       + SIN(alpha) * Cx

      Byr= COS(alpha) * By + (1- COS(alpha))* dot_prod * ny
     1                       + SIN(alpha) * Cy

      Bzr= COS(alpha) * Bz + (1- COS(alpha))* dot_prod * nz
     1                       + SIN(alpha) * Cz
  
! ----------------------------------------------------------------------

      END SUBROUTINE ROTATE_VECTOR

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE CROSS_PRODUCT (  nx, ny, nz, Vx, Vy, Vz,
     1                                        Cx, Cy, Cz )

! ----------------------------------------------------------------------
! Written By: Dinshaw Balsara, Vladimir Florinski, Sethu and Swati Sharma

! This subroutine calculates the cross product of two vectors

! ----------------------------------------------------------------------
!$acc routine seq


      IMPLICIT NONE

      REAL  nx, ny, nz, Vx, Vy, Vz, Cx, Cy, Cz

! ----------------------------------------------------------------------
!cross product n cross B
      Cx= ny * Vz - nz * Vy
      Cy= nz * Vx - nx * Vz
      Cz= nx * Vy - ny * Vx


! ----------------------------------------------------------------------
      END SUBROUTINE CROSS_PRODUCT

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      SUBROUTINE FLOOR_HELIO_PRS_1_LAYER

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,

     1           iz_patch,

     1         U, V, entropygr)

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara

! If "iz_patch == 1" then this subroutine goes to the layer of zones with
! "izz = 1" and provides a minimum pressure floor to all the variables.

! ----------------------------------------------------------------------

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER nr_zones, ntp_zones, NGhost,
     1        n_cc_components, n_cc_modes, iz_patch

! ------------------------------
! FLUID VARIABLES:-

      REAL U ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

      REAL entropygr ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   4)

      REAL V ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components, n_cc_modes)

! ----------------------------------------------------------------------

! Variables past this point are internal to the subroutine.

      INTEGER ixx, iyy, izz, iprs_flag

      REAL Upt ( n_cc_components),
     1     Vpt ( n_cc_components), ent, B0 ( 3)

      REAL, PARAMETER :: ob_8pi = 1.0 / ( 8.0 * 4.0 * ATAN ( 1.0)) 

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#if ( CHOOSE_APPL_HELIO == 1 && (IPDE_SYSTEM == 1 || IPDE_SYSTEM == 2) )

      IF ( iz_patch .EQ. 1 ) THEN

! To prevent pressure in the first radial zone from sinking for the
! heliosphere problem, we provide a floor.

!$OMP PARALLEL DO PRIVATE ( ixx, iyy, izz, iprs_flag, ent, B0, Upt, Vpt)&
!$OMP&             SHARED ( nr_zones, ntp_zones, NGhost,                &
!$OMP&    gamma, gamma_euler, isospeed, barotropic_stiffness,           &
!$OMP&    molewt, gasconst, rhofloor, prsfloor, tempfloor,              &
!$OMP&    magB_small_cgl, eps_bdry_cgl,                                 &
!$OMP&    V, U, entropygr)                                              &
!$OMP&            SCHEDULE ( static, 1) DEFAULT ( NONE) COLLAPSE ( 2)

      DO iyy = 1 - NGhost, ntp_zones + NGhost
      DO ixx = 1 - NGhost, ntp_zones + NGhost

      izz = 1

! --------------------------------------------------

          B0 ( :) = 0.0

          Upt ( :) = U ( ixx, iyy, izz, :, 1)

          ent = entropygr ( ixx, iyy, izz, 1)

          CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, Upt, Vpt)

          V ( ixx, iyy, izz, :, 1) = Vpt ( :)
          U ( ixx, iyy, izz, :, 1) = Upt ( :)

! --------------------------------------------------
! This value of "3.8914E-5" is actually "prs_sw" evaluated at the inner boundary
! of the heliosphere boundary. For CGL, we will have to fix two such pressures.

#if (IPDE_SYSTEM == 1)

! For MHD, we have a single pressure in the 5th component of "V" and a single
! total Energy in "U".

!        IF ( V ( ixx, iyy, izz, 5, 1) .LT. 3.8914E-5) THEN

        U ( ixx, iyy, izz, 5, 1) = U ( ixx, iyy, izz, 5, 1)
     1   + (3.8914E-5 - V ( ixx, iyy, izz, 5, 1)) / (gamma - 1.0)

        V ( ixx, iyy, izz, 5, 1) = 3.8914E-5

!        END IF

#elif (IPDE_SYSTEM == 2)

! For CGL we have "V (5) = P_parallel ; V (9) = P_perpendicular" and    
! "U (5) = Total energy density ; U (9) = P_parallel - P_perpendicular".
! In the setup we want "P_parallel == P_perpendicular" exactly.

        V ( ixx, iyy, izz, 5, 1) = 3.8914E-5

        V ( ixx, iyy, izz, 9, 1) = 3.8914E-5

        U ( ixx, iyy, izz, 5, 1) = 1.5 * 3.8914E-5
     1+ 0.5 * V ( ixx, iyy, izz, 1, 1) * ( V ( ixx, iyy, izz, 2, 1)**2  
     1  + V ( ixx, iyy, izz, 3, 1)**2 + V ( ixx, iyy, izz, 4, 1)**2 )   
     1 + ( V ( ixx, iyy, izz, 6, 1)**2 + V ( ixx, iyy, izz, 7, 1)**2
     1   + V ( ixx, iyy, izz, 8, 1)**2 ) * ob_8pi

        U ( ixx, iyy, izz, 9, 1) = 0.0  ! Pressures are isotropic.

#endif

! --------------------------------------------------

      END DO    ! End of "DO ixx = 1 - NGhost, ntp_zones + NGhost".
      END DO    ! End of "DO iyy = 1 - NGhost, ntp_zones + NGhost".
!$OMP END PARALLEL DO

      END IF    ! End of "IF ( iz_patch .EQ. 1 ) THEN".

#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END SUBROUTINE FLOOR_HELIO_PRS_1_LAYER

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

#define HELIO_PRINT_VALUES 0

! ---------------------------------------------------------------------- 
! rm -f x a.F b.F ; cp Heliosphere_src a.F ; cp PDE_3D_src b.F ; gfortran -fdefault-real-8 -o x a.F b.F ; ./x

#if( HELIO_PRINT_VALUES == 1 )
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      PROGRAM TEST_HELIO

! ----------------------------------------------------------------------

      USE HELIOSPHERE_MOD
      IMPLICIT NONE
                                                                         
#include <cubesphere/000com/riemann.com>

      REAL :: x, y, z, time, rmin, rmax, zone_size,
     1        xcen_sector, ycen_sector, zcen_sector,
     1        xvert1_sector, yvert1_sector, zvert1_sector,

     1        ent, B0 ( 3)


      REAL ::  nx, ny, nz, Bx, By, Bz, alpha2,
     1          ratio_v, beta, b_tot, v_tot, timestop,
     1          prs_sw, beta1_au, c_s , c_s1, ratio_1au


      INTEGER :: n_cc_components, iproblem,nsilo_images,
     1            iprs_flag, bflag

      REAL, DIMENSION ( NVAR) :: U_pt, V_pt

      REAL :: A_pt ( 3), r

      n_cc_components = NVAR

      CALL GET_HELIO_PARAMETERS (
     1     iproblem, nsilo_images,
     1     rmin, rmax, gamma, gasconst, molewt, timestop,
     1     omega_passive, isospeed, tempfloor,
     1     tau_phys_cgl, f_lowbd_cgl,
     1     neta_f_cgl, kappa_f_cgl, magB_small_cgl, eps_bdry_cgl)


! ----------------------------------------------------------------------

      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" mu            = ", mu
      WRITE(*,fmt='(A,ES15.9E2)')" Rgas          = ", Rgas
      WRITE(*,fmt='(A,ES15.9E2)')" solar_Mdot    = ", solar_Mdot
      WRITE(*,fmt='(A,ES15.9E2)')" solar_Vsnd    = ", solar_Vsnd
      WRITE(*,fmt='(A,ES15.9E2)')" solar_uc      = ", solar_uc
      WRITE(*,fmt='(A,ES15.9E2)')" solar_rho     = ", solar_rho
      WRITE(*,fmt='(A,ES15.9E2)')" solar_prs     = ", solar_prs
      WRITE(*,fmt='(A,ES15.9E2)')" solar_Omega   = ", solar_Omega

      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" unit_dist     = ", unit_dist
      WRITE(*,fmt='(A,ES15.9E2)')" unit_time     = ", unit_time
      WRITE(*,fmt='(A,ES15.9E2)')" unit_mass     = ", unit_mass

      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" code_density  = ", code_density
      WRITE(*,fmt='(A,ES15.9E2)')" code_prs      = ", code_prs
      WRITE(*,fmt='(A,ES15.9E2)')" code_Omega    = ", code_Omega
      WRITE(*,fmt='(A,ES15.9E2)')" r_lism        = ", r_lism
      WRITE(*,fmt='(A,ES15.9E2)')" B0_sw        = ",B0_sw

! --------------------
     
      x = 1* R_inner     !/ SQRT(3.0);
      y = 0.0 * R_inner / SQRT(3.0);
      z = 0.0 * R_inner / SQRT(3.0);
     

      r = SQRT( x**2 + y**2 + z**2 )

      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,*)" At 1.0 R_inner:"
      WRITE(*,*)" ==========="

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )
 
      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

        v_tot= SQRT(V_pt(2)**2 + V_pt(3)**2+ V_pt(4)**2 )
        
        b_tot= SQRT(V_pt(6)**2 + V_pt(7)**2+ V_pt(8)**2 )
        
        c_s = SQRT(gamma* V_pt(5)/ V_pt(1))

        ratio_v= v_tot / c_s

        beta= V_pt(5) * 8 * M_PI / (b_tot**2 )

        prs_sw  = (5.0*amu)* Rgas * solar_T / mu                 ! Physical pressure of SW at 1 AU

       prs_sw =  prs_sw
     1         / (unit_mass / (unit_dist * unit_time**2) )           ! Pressure of SW in code units.

      beta1_au= prs_sw* 8*M_PI/ B0_sw**2
      c_s1= SQRT(gamma* (5.0*amu)* Rgas* solar_T / amu)
      ratio_1au = solar_uc/c_s1

      WRITE(*,*)" "
      WRITE(*,*)" x, y, z, r = ", x, y, z, r
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" rho = ", V_pt ( 1)
      WRITE(*,fmt='(A,ES15.9E2)')"|vx| = ", ABS(V_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|vy| = ", ABS(V_pt ( 3))
      WRITE(*,fmt='(A,ES15.9E2)')"|vz| = ", ABS(V_pt ( 4))
      WRITE(*,fmt='(A,ES15.9E2)')" prs = ", V_pt ( 5)
      WRITE(*,fmt='(A,ES15.9E2)')"|Bx| = ", ABS(V_pt ( 6))
      WRITE(*,fmt='(A,ES15.9E2)')"|By| = ", ABS(V_pt ( 7))
      WRITE(*,fmt='(A,ES15.9E2)')"|Bz| = ", ABS(V_pt ( 8))
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')"|Ax| = ", ABS(A_pt ( 1))
      WRITE(*,fmt='(A,ES15.9E2)')"|Ay| = ", ABS(A_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|Az| = ", ABS(A_pt ( 3))
      WRITE(*,*)"v_tot = ",v_tot
      WRITE(*,*)"ratio_v = ", ratio_v
      WRITE(*,*)"beta = ",beta




! --------------------

      x = 2* R_inner     !/ SQRT(3.0);
      y = 0.0 * R_inner / SQRT(3.0);
      z = 0.0 * R_inner / SQRT(3.0);


      r = SQRT( x**2 + y**2 + z**2 )

      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,*)" At 2.0 R_inner:"
      WRITE(*,*)" ==========="

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )

      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

        v_tot= SQRT(V_pt(2)**2 +V_pt(3)**2+ V_pt(4)**2 )

        b_tot= SQRT(V_pt(6)**2 +V_pt(7)**2+ V_pt(8)**2 )

        c_s = SQRT(gamma* V_pt(5)/ V_pt(1))

        ratio_v= v_tot / c_s

        beta= V_pt(5) * 8 * M_PI / (b_tot**2 )



      WRITE(*,*)" "
      WRITE(*,*)" x, y, z, r = ", x, y, z, r
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" rho = ", V_pt ( 1)
      WRITE(*,fmt='(A,ES15.9E2)')"|vx| = ", ABS(V_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|vy| = ", ABS(V_pt ( 3))
      WRITE(*,fmt='(A,ES15.9E2)')"|vz| = ", ABS(V_pt ( 4))
      WRITE(*,fmt='(A,ES15.9E2)')" prs = ", V_pt ( 5)
      WRITE(*,fmt='(A,ES15.9E2)')"|Bx| = ", ABS(V_pt ( 6))
      WRITE(*,fmt='(A,ES15.9E2)')"|By| = ", ABS(V_pt ( 7))
      WRITE(*,fmt='(A,ES15.9E2)')"|Bz| = ", ABS(V_pt ( 8))
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')"|Ax| = ", ABS(A_pt ( 1))
      WRITE(*,fmt='(A,ES15.9E2)')"|Ay| = ", ABS(A_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|Az| = ", ABS(A_pt ( 3))
      WRITE(*,*)"v_tot = ",v_tot
      WRITE(*,*)"ratio_v = ", ratio_v
      WRITE(*,*)"beta = ",beta


! --------------------

      x =3.0 * R_inner       !/ SQRT(3.0);
      y = 0.0 * R_inner / SQRT(3.0);
      z = 0.0 * R_inner / SQRT(3.0);

      r = SQRT( x**2 + y**2 + z**2 )
      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,*)" At 3.0 R_inner:"
      WRITE(*,*)" ============="

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )

      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

        v_tot= SQRT(V_pt(2)**2 +V_pt(3)**2+ V_pt(4)**2 )

        b_tot= SQRT(V_pt(6)**2 +V_pt(7)**2+ V_pt(8)**2 )

        c_s = SQRT(gamma* V_pt(5)/ V_pt(1))

        ratio_v= v_tot / c_s

        beta= V_pt(5) * 8 * M_PI / (b_tot**2 )


      WRITE(*,*)" "
      WRITE(*,*)" x, y, z, r = ", x, y, z, r
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" rho = ", V_pt ( 1)
      WRITE(*,fmt='(A,ES15.9E2)')"|vx| = ", ABS(V_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|vy| = ", ABS(V_pt ( 3))
      WRITE(*,fmt='(A,ES15.9E2)')"|vz| = ", ABS(V_pt ( 4))
      WRITE(*,fmt='(A,ES15.9E2)')" prs = ", V_pt ( 5)
      WRITE(*,fmt='(A,ES15.9E2)')"|Bx| = ", ABS(V_pt ( 6))
      WRITE(*,fmt='(A,ES15.9E2)')"|By| = ", ABS(V_pt ( 7))
      WRITE(*,fmt='(A,ES15.9E2)')"|Bz| = ", ABS(V_pt ( 8))
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')"|Ax| = ", ABS(A_pt ( 1))
      WRITE(*,fmt='(A,ES15.9E2)')"|Ay| = ", ABS(A_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|Az| = ", ABS(A_pt ( 3))
      WRITE(*,*)"v_tot = ",v_tot
      WRITE(*,*)"ratio_v = ", ratio_v
      WRITE(*,*)"beta = ",beta

! --------------------

      x =4.0 * R_inner       !/ SQRT(3.0);
      y = 0.0 * R_inner / SQRT(3.0);
      z = 0.0 * R_inner / SQRT(3.0);

      r = SQRT( x**2 + y**2 + z**2 )
      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,*)" At 4.0 R_inner:"
      WRITE(*,*)" ============="

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )

      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

        v_tot= SQRT(V_pt(2)**2 +V_pt(3)**2+ V_pt(4)**2 )

        b_tot= SQRT(V_pt(6)**2 +V_pt(7)**2+ V_pt(8)**2 )

        c_s = SQRT(gamma* V_pt(5)/ V_pt(1))

        ratio_v= v_tot / c_s

        beta= V_pt(5) * 8 * M_PI / (b_tot**2 )


      WRITE(*,*)" "
      WRITE(*,*)" x, y, z, r = ", x, y, z, r
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" rho = ", V_pt ( 1)
      WRITE(*,fmt='(A,ES15.9E2)')"|vx| = ", ABS(V_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|vy| = ", ABS(V_pt ( 3))
      WRITE(*,fmt='(A,ES15.9E2)')"|vz| = ", ABS(V_pt ( 4))
      WRITE(*,fmt='(A,ES15.9E2)')" prs = ", V_pt ( 5)
      WRITE(*,fmt='(A,ES15.9E2)')"|Bx| = ", ABS(V_pt ( 6))
      WRITE(*,fmt='(A,ES15.9E2)')"|By| = ", ABS(V_pt ( 7))
      WRITE(*,fmt='(A,ES15.9E2)')"|Bz| = ", ABS(V_pt ( 8))
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')"|Ax| = ", ABS(A_pt ( 1))
      WRITE(*,fmt='(A,ES15.9E2)')"|Ay| = ", ABS(A_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|Az| = ", ABS(A_pt ( 3))
      WRITE(*,*)"v_tot = ",v_tot
      WRITE(*,*)"ratio_v = ", ratio_v
      WRITE(*,*)"beta = ",beta


! --------------------

      x =5 * R_inner ! / SQRT(3.0);
      y = 0.0 * R_inner / SQRT(3.0);
      z = 0.0 * R_inner / SQRT(3.0);

      r = SQRT( x**2 + y**2 + z**2 )
      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,*)" At 5 R_inner:"
      WRITE(*,*)" ============="

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )

      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

        v_tot= SQRT(V_pt(2)**2 +V_pt(3)**2+ V_pt(4)**2 )

        b_tot= SQRT(V_pt(6)**2 +V_pt(7)**2+ V_pt(8)**2 )

        c_s = SQRT(gamma* V_pt(5)/ V_pt(1))

        ratio_v= v_tot / c_s


        beta= V_pt(5) * 8 * M_PI / (b_tot**2 )


      WRITE(*,*)" "
      WRITE(*,*)" x, y, z, r = ", x, y, z, r
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" rho = ", V_pt ( 1)
      WRITE(*,fmt='(A,ES15.9E2)')"|vx| = ", ABS(V_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|vy| = ", ABS(V_pt ( 3))
      WRITE(*,fmt='(A,ES15.9E2)')"|vz| = ", ABS(V_pt ( 4))
      WRITE(*,fmt='(A,ES15.9E2)')" prs = ", V_pt ( 5)
      WRITE(*,fmt='(A,ES15.9E2)')"|Bx| = ", ABS(V_pt ( 6))
      WRITE(*,fmt='(A,ES15.9E2)')"|By| = ", ABS(V_pt ( 7))
      WRITE(*,fmt='(A,ES15.9E2)')"|Bz| = ", ABS(V_pt ( 8))
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')"|Ax| = ", ABS(A_pt ( 1))
      WRITE(*,fmt='(A,ES15.9E2)')"|Ay| = ", ABS(A_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|Az| = ", ABS(A_pt ( 3))
      WRITE(*,*)"v_tot = ",v_tot
      WRITE(*,*)"ratio_v = ", ratio_v
      WRITE(*,*)"beta = ",beta


! --------------------
      x =6.0 * R_inner   ! / SQRT(3.0);
      y = 0.0 * R_inner / SQRT(3.0);
      z = 0.0 * R_inner / SQRT(3.0);

      r = SQRT( x**2 + y**2 + z**2 )
      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,*)" At 6 R_inner:"
      WRITE(*,*)" ============="

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )

      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

        v_tot= SQRT(V_pt(2)**2 +V_pt(3)**2+ V_pt(4)**2 )

        b_tot= SQRT(V_pt(6)**2 +V_pt(7)**2+ V_pt(8)**2 )

        c_s = SQRT(gamma* V_pt(5)/ V_pt(1))

        ratio_v= v_tot / c_s

        beta= V_pt(5) * 8 * M_PI / (b_tot**2 )


      WRITE(*,*)" "
      WRITE(*,*)" x, y, z, r = ", x, y, z, r
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" rho = ", V_pt ( 1)
      WRITE(*,fmt='(A,ES15.9E2)')"|vx| = ", ABS(V_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|vy| = ", ABS(V_pt ( 3))
      WRITE(*,fmt='(A,ES15.9E2)')"|vz| = ", ABS(V_pt ( 4))
      WRITE(*,fmt='(A,ES15.9E2)')" prs = ", V_pt ( 5)
      WRITE(*,fmt='(A,ES15.9E2)')"|Bx| = ", ABS(V_pt ( 6))
      WRITE(*,fmt='(A,ES15.9E2)')"|By| = ", ABS(V_pt ( 7))
      WRITE(*,fmt='(A,ES15.9E2)')"|Bz| = ", ABS(V_pt ( 8))
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')"|Ax| = ", ABS(A_pt ( 1))
      WRITE(*,fmt='(A,ES15.9E2)')"|Ay| = ", ABS(A_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|Az| = ", ABS(A_pt ( 3))
      WRITE(*,*)"v_tot = ",v_tot
      WRITE(*,*)"ratio_v = ", ratio_v
      WRITE(*,*)"beta = ",beta

! --------------------

      x =15.0 * R_inner !/ SQRT(3.0);
      y = 0.0 * R_inner / SQRT(3.0);
      z = 0.0 * R_inner / SQRT(3.0);

      r = SQRT( x**2 + y**2 + z**2 )
      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,*)" At 15 R_inner:"
      WRITE(*,*)" ============="

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )

      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

        v_tot= SQRT(V_pt(2)**2 +V_pt(3)**2+ V_pt(4)**2 )

        b_tot= SQRT(V_pt(6)**2 +V_pt(7)**2+ V_pt(8)**2 )

         c_s = SQRT(gamma* V_pt(5)/ V_pt(1))

        ratio_v= v_tot / c_s
        beta= V_pt(5) * 8 * M_PI / (b_tot**2 )

      WRITE(*,*)" "
      WRITE(*,*)" x, y, z, r = ", x, y, z, r
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" rho = ", V_pt ( 1)
      WRITE(*,fmt='(A,ES15.9E2)')"|vx| = ", ABS(V_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|vy| = ", ABS(V_pt ( 3))
      WRITE(*,fmt='(A,ES15.9E2)')"|vz| = ", ABS(V_pt ( 4))
      WRITE(*,fmt='(A,ES15.9E2)')" prs = ", V_pt ( 5)
      WRITE(*,fmt='(A,ES15.9E2)')"|Bx| = ", ABS(V_pt ( 6))
      WRITE(*,fmt='(A,ES15.9E2)')"|By| = ", ABS(V_pt ( 7))
      WRITE(*,fmt='(A,ES15.9E2)')"|Bz| = ", ABS(V_pt ( 8))
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')"|Ax| = ", ABS(A_pt ( 1))
      WRITE(*,fmt='(A,ES15.9E2)')"|Ay| = ", ABS(A_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|Az| = ", ABS(A_pt ( 3))
      WRITE(*,*)"v_tot = ",v_tot
      WRITE(*,*)"ratio_v = ", ratio_v
      WRITE(*,*)"beta = ",beta

! --------------------

      x =25.0 * R_inner !/ SQRT(3.0);
      y = 0.0 * R_inner / SQRT(3.0);
      z = 0.0 * R_inner / SQRT(3.0);

      r = SQRT( x**2 + y**2 + z**2 )
      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,*)" At 25 R_inner:"
      WRITE(*,*)" ============="

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )

      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

        v_tot= SQRT(V_pt(2)**2 +V_pt(3)**2+ V_pt(4)**2 )

        b_tot= SQRT(V_pt(6)**2 +V_pt(7)**2+ V_pt(8)**2 )

        c_s = SQRT(gamma* V_pt(5)/ V_pt(1))

        ratio_v= v_tot / c_s
         beta= V_pt(5) * 8 * M_PI / (b_tot**2 )


      WRITE(*,*)" "
      WRITE(*,*)" x, y, z, r = ", x, y, z, r
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" rho = ", V_pt ( 1)
      WRITE(*,fmt='(A,ES15.9E2)')"|vx| = ", ABS(V_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|vy| = ", ABS(V_pt ( 3))
      WRITE(*,fmt='(A,ES15.9E2)')"|vz| = ", ABS(V_pt ( 4))
      WRITE(*,fmt='(A,ES15.9E2)')" prs = ", V_pt ( 5)
      WRITE(*,fmt='(A,ES15.9E2)')"|Bx| = ", ABS(V_pt ( 6))
      WRITE(*,fmt='(A,ES15.9E2)')"|By| = ", ABS(V_pt ( 7))
      WRITE(*,fmt='(A,ES15.9E2)')"|Bz| = ", ABS(V_pt ( 8))
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')"|Ax| = ", ABS(A_pt ( 1))
      WRITE(*,fmt='(A,ES15.9E2)')"|Ay| = ", ABS(A_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|Az| = ", ABS(A_pt ( 3))
      WRITE(*,*)"v_tot = ",v_tot
      WRITE(*,*)"ratio_v = ", ratio_v
      WRITE(*,*)"beta = ",beta

! --------------------

      x =30* R_inner   !/ SQRT(3.0);
      y = 0.0 * R_inner / SQRT(3.0);
      z = 0.0 * R_inner / SQRT(3.0);

      r = SQRT( x**2 + y**2 + z**2 )
      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,*)" At  30 R_inner:"
      WRITE(*,*)" ============="

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )

      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

        v_tot= SQRT(V_pt(2)**2 +V_pt(3)**2+ V_pt(4)**2 )

        b_tot= SQRT(V_pt(6)**2 +V_pt(7)**2+ V_pt(8)**2 )

        c_s = SQRT(gamma* V_pt(5)/ V_pt(1))

        ratio_v= v_tot / c_s

        beta= V_pt(5) * 8 * M_PI / (b_tot**2 )


      WRITE(*,*)" "
      WRITE(*,*)" x, y, z, r = ", x, y, z, r
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" rho = ", V_pt ( 1)
      WRITE(*,fmt='(A,ES15.9E2)')"|vx| = ", ABS(V_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|vy| = ", ABS(V_pt ( 3))
      WRITE(*,fmt='(A,ES15.9E2)')"|vz| = ", ABS(V_pt ( 4))
      WRITE(*,fmt='(A,ES15.9E2)')" prs = ", V_pt ( 5)
      WRITE(*,fmt='(A,ES15.9E2)')"|Bx| = ", ABS(V_pt ( 6))
      WRITE(*,fmt='(A,ES15.9E2)')"|By| = ", ABS(V_pt ( 7))
      WRITE(*,fmt='(A,ES15.9E2)')"|Bz| = ", ABS(V_pt ( 8))
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')"|Ax| = ", ABS(A_pt ( 1))
      WRITE(*,fmt='(A,ES15.9E2)')"|Ay| = ", ABS(A_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|Az| = ", ABS(A_pt ( 3))
      WRITE(*,*)"v_tot = ",v_tot
      WRITE(*,*)"ratio_v = ", ratio_v
      WRITE(*,*)"beta = ",beta


! --------------------


      x =1 * R_outer - 1!/ SQRT(3.0);
      y = 0.0 * R_inner / SQRT(3.0);
      z = 0.0 * R_inner / SQRT(3.0);

      r = SQRT( x**2 + y**2 + z**2 )
      WRITE(*,*)" "
      WRITE(*,*)" "
      WRITE(*,*)" At R_outer-1:"
      WRITE(*,*)" ============="

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )

      CALL CONS_TO_PRIM (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

        v_tot= SQRT(V_pt(2)**2 +V_pt(3)**2+ V_pt(4)**2 )

        b_tot= SQRT(V_pt(6)**2 +V_pt(7)**2+ V_pt(8)**2 )
        c_s = SQRT(gamma* V_pt(5)/ V_pt(1))

        ratio_v= v_tot / c_s


        beta= V_pt(5) * 8 * M_PI / (b_tot**2 )


      WRITE(*,*)" "
      WRITE(*,*)" x, y, z, r = ", x, y, z, r
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')" rho = ", V_pt ( 1)
      WRITE(*,fmt='(A,ES15.9E2)')"|vx| = ", ABS(V_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|vy| = ", ABS(V_pt ( 3))
      WRITE(*,fmt='(A,ES15.9E2)')"|vz| = ", ABS(V_pt ( 4))
      WRITE(*,fmt='(A,ES15.9E2)')" prs = ", V_pt ( 5)
      WRITE(*,fmt='(A,ES15.9E2)')"|Bx| = ", ABS(V_pt ( 6))
      WRITE(*,fmt='(A,ES15.9E2)')"|By| = ", ABS(V_pt ( 7))
      WRITE(*,fmt='(A,ES15.9E2)')"|Bz| = ", ABS(V_pt ( 8))
      WRITE(*,*)" "
      WRITE(*,fmt='(A,ES15.9E2)')"|Ax| = ", ABS(A_pt ( 1))
      WRITE(*,fmt='(A,ES15.9E2)')"|Ay| = ", ABS(A_pt ( 2))
      WRITE(*,fmt='(A,ES15.9E2)')"|Az| = ", ABS(A_pt ( 3))
      WRITE(*,*)"v_tot = ",v_tot
      WRITE(*,*)"ratio_v = ", ratio_v
      WRITE(*,*)"beta = ",beta

! ----------------------------------------------------------------------

      END PROGRAM TEST_HELIO

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
#endif
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
