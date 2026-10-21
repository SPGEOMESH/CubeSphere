








      SUBROUTINE PROBLEM_GENERATOR_U_A

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )


      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      REAL :: x, y, z, time, rmin, rmax, zone_size,
     1        xcen_sector, ycen_sector, zcen_sector,
     1        xvert1_sector, yvert1_sector, zvert1_sector

      INTEGER :: n_cc_components, iproblem

      REAL, DIMENSION ( n_cc_components) :: U_pt, V_pt

      REAL :: A_pt ( 3)

! ----------------------------------------------------------------------
! Variables past here are local to this subroutine.

      REAL :: rad, costheta, sintheta, cosphi, sinphi, R,
     1        rho, prs, v0, vr, vtheta, vphi, vx, vy, vz, bfx, bfy, bfz,
     1        rho0, prs0, vlx0, vly0, vlz0, bfx0, bfy0, bfz0,
     1        rho1, rho2, prs1, prs2, vlx1, vlx2, vly1, vly2,
     1        vlz1, vlz2, bfx1, bfx2, bfy1, bfy2, bfz1, bfz2,
     1        radius, distance, bmag,

     1        m1, m2, m3, phi, theta

      REAL :: d0, u0, p0, gammaa, d, uu, p, a, b,
     1        gamma1, dcur, dprev, c, f1, f2, f, df,
     1        kappa, mu, sq4pi

      REAL :: alpha, beta, tau, r0, omega, mach

      REAL :: Ar, Br, Br0, Bphi, Bphi1, Bphi2, eps, tempaa,
     1        Aphi, B_0, band_width, cosomegat, sinomegat,
     1        cosphi_mag, sinphi_mag, costheta_mag, sintheta_mag,
     1        costilt, sintilt, fast_slow, tilt,
     1        theta_mag, phi_mag, trans_width, x_mag, y_mag, z_mag,
     1        time_use

      REAL :: V_CAK ( 3), cak_rhoc, xs, ys, zs, Rs, ent, B0 ( 3)

      INTEGER iprs_flag, bflag

! ----------------------------------------------------------------------

!$acc routine( GET_MAGNETIC_FIELD ) seq
!$acc routine( CAK_INIT_1D) seq
!$acc routine( PRIM_TO_CONS) seq
!$acc routine( HELIO_INITIALIZE) seq

! ----------------------------------------------------------------------

      R = SQRT ( x**2 + y**2 + z**2)

      theta = ATAN2 ( SQRT ( x**2 + y**2), z)
      IF ( theta < 0.0) theta = theta + 2.0 * pi

      phi = ATAN2 ( y, x)
      IF ( phi < 0.0) phi = phi + 2.0 * pi

      costheta = COS ( theta)
      sintheta = SIN ( theta)

      cosphi = COS ( phi)
      sinphi = SIN ( phi)

! --------------------

      B0 ( :) = 0.0

#if( USE_B0B1_FORMULATION == 1)

      bflag = 1
      time_use = time

      CALL GET_MAGNETIC_FIELD ( bflag, time_use, x, y, z, xs, ys, zs)

      B0 ( 1) = xs
      B0 ( 2) = ys
      B0 ( 3) = zs

#endif

! ----------------------------------------------------------------------
      IF ( iproblem == 0) THEN
! ----------------------------------------------------------------------
! Testing the MPI messaging

#define TESTING_OF_MPI_MESSAGING 1
#define TESTING_OF_CODE_T1 0

! ----------
#if( TESTING_OF_CODE_T1 == 1)
! ----------

!        IF ( R < rmin ) THEN
!          x = x * rmin / R
!          y = y * rmin / R
!          z = z * rmin / R
!          R = rmin
!        END IF

! fluid variables.
        rho0   = 1.0
        prs0   = 1.0 / gamma

        rho    = rho0 / R**2
        prs    = prs0 * ( ( rho / rho0)**gamma)

        vr     = 1.0
        vtheta = 0.0
        vphi   = 0.0     ! omega_active ( 3) * rmin * sintheta

        vx     = vr * sintheta * cosphi + vtheta * costheta * cosphi
     1         - vphi * sinphi
        vy     = vr * sintheta * sinphi + vtheta * costheta * sinphi
     1         + vphi * cosphi
        vz     = vr * costheta - vtheta * sintheta

! init the dipolar m-field.
        m1 = dipole_B ( 1)
        m2 = dipole_B ( 2)
        m3 = dipole_B ( 3)

#if( USE_B0B1_FORMULATION == 1)
        m1 = 0.0
        m2 = 0.0
        m3 = 0.0
#endif

        bfx  = - ( m1/R**3) + (3*x*(m1*x + m2*y + m3*z)) / R**5
        bfy  = - ( m2/R**3) + (3*y*(m1*x + m2*y + m3*z)) / R**5
        bfz  = - ( m3/R**3) + (3*z*(m1*x + m2*y + m3*z)) / R**5

! Initializa vector potential Ax, Ay, Az
        A_pt ( 1) = (m2*z - m3*y) / R**3
        A_pt ( 2) = (m3*x - m1*z) / R**3
        A_pt ( 3) = (m1*y - m2*x) / R**3

! --------------------

        V_pt ( 1) = rho

        V_pt ( 2) = vx
        V_pt ( 3) = vy
        V_pt ( 4) = vz

        V_pt ( 5) = prs

        V_pt ( 6) = bfx
        V_pt ( 7) = bfy
        V_pt ( 8) = bfz

        ent = V_pt ( 5) / V_pt ( 1)**( gamma - 1.0)

#if( IPDE_SYSTEM == 2 )
        V_pt ( 9) = V_pt ( 5)
#endif

#if( IPDE_SYSTEM == 3 )
        V_pt ( 9) = V_pt ( 5)
        V_pt ( 10) = V_pt ( 5)
#endif

        CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

! ----------
#elif( TESTING_OF_MPI_MESSAGING == 1)
! ----------
! Initializa vector potential Ax, Ay, Az
        A_pt ( 1) = y * z
        A_pt ( 2) = z * x * 2.0
        A_pt ( 3) = x * y * 3.0

        U_pt ( 1) = ABS( x )
        U_pt ( 2) = x
        U_pt ( 3) = y
        U_pt ( 4) = z
        U_pt ( 5) = ABS( z )

! The above choice of vector potential "A", would yield the
! following:
!   Bx  =  dAz/dy - dAy/dz  =  3x - 2x  =   x
!   By  =  dAx/dz - dAz/dx  =   y - 3y  = -2y
!   Bz  =  dAy/dx - dAx/dy  =  2z -  z  =   z

        U_pt ( 6) =  1.0 * x
        U_pt ( 7) = -2.0 * y
        U_pt ( 8) =  1.0 * z

#if( IPDE_SYSTEM == 2 )
        U_pt ( 9) = 0.0
#endif

#if( IPDE_SYSTEM == 10 )
        U_pt ( 9) = 0.0
        U_pt ( 10) = 0.0
#endif

! ----------
#endif
! ----------
#undef TESTING_OF_MPI_MESSAGING
#undef TESTING_OF_CODE_T1
! ----------------------------------------------------------------------
      ELSE IF ( iproblem == 1) THEN
! ----------------------------------------------------------------------

! Ref: Section 4.2.1 of L. Ivan et al., JCP 2013, 255, pp 205 - 227
! Eq. 12 is used here.
! sect_div = 0; face_division = 3; NRadial = 8
! gamma = 7.0 / 5.0
! rmin = 2.0; rmax = 3.5
! iproblem = 1
! timestop = 0.2
! set "bcarr_r" to 17 and call "ADD_BOUNDARY" after calling "PAD_BOUNDARY_R".

! switch on the accuracy calculation parts.
! switch on  directive "if ( 1 == 1)" after limit call for accuarcy
! calculation. And the same for printing the L1, Linf errors, located above
! "IF ( time .GE. timestop) EXIT" line.


! ----------------------------------------------------------------------
      ELSE IF ( iproblem == 2) THEN
! ----------------------------------------------------------------------

! ----------
! Set up MHD blast problem.
! ntp_global = ntp_zones * ntp_patches = 180 or higher
!  nr_global =  nr_zones *  nr_patches = 180
! rmin = 7.0; rmax = 10.0
! iproblem = 2
! timestop = 0.13
! "bfx1, bfy1, bfz1" below should be: 10.0 / SQRT ( 3.0)

! ----------
! 2nd way: Set up MHD blast problem.
! ntp_global = ntp_zones * ntp_patches = 180 or higher
!  nr_global =  nr_zones *  nr_patches = 180
! rmin = 7.0; rmax = 16.0
! iproblem = 2
! timestop = 0.06

! "bfx1, bfy1, bfz1" below should be: 100.0 / SQRT ( 3.0)

! Outflow boundary conditions at PAD_BOUNDARY_R. Assigned as,
!      bcarr_r ( :, :, 1)  = 2
!      bcarr_r ( :, :, 2)  = 2

! Inner Inflow boundary conditions must not be enforced at ADD_BOUNDARY.
!      ibdy_inner = 0
!      ibdy_outer = 0

! ----------

        rho1 = 1.0
        rho2 = 1.0

        prs1 = 1000.0
        prs2 = 0.1

        vlx1 = 0.0
        vlx2 = 0.0

        vly1 = 0.0
        vly2 = 0.0

        vlz1 = 0.0
        vlz2 = 0.0

        bfx1 = 100.0 / SQRT ( 3.0)
        bfx2 = bfx1

        bfy1 = 100.0 / SQRT ( 3.0)
        bfy2 = bfy1

        bfz1 = 100.0 / SQRT ( 3.0)
        bfz2 = bfz1

        radius = 1.0

!        distance = SQRT( ( x - xcen_sector)**2 + ( y - ycen_sector)**2
!     1                 + ( z - zcen_sector)**2 )

        distance = SQRT ( ( x - (rmax+rmin)*0.5)**2 + y**2 + z**2)

! Enable the following line to make this as a 2D problem, along the xy-plane.
!        bfz1 = 0.0; bfz2 = 0.0;  distance = SQRT ( x**2 + y**2 )

        IF ( distance <= radius ) THEN

        V_pt ( 1) = rho1
        V_pt ( 2) = vlx1
        V_pt ( 3) = vly1
        V_pt ( 4) = vlz1

        V_pt ( 5) = prs1

        V_pt ( 6) = bfx1
        V_pt ( 7) = bfy1
        V_pt ( 8) = bfz1

        A_pt ( 1) = bfy1 * z
        A_pt ( 2) = bfz1 * x
        A_pt ( 3) = bfx1 * y

        ELSE

        V_pt ( 1) = rho2

        V_pt ( 2) = vlx2
        V_pt ( 3) = vly2
        V_pt ( 4) = vlz2

        V_pt ( 5) = prs2

        V_pt ( 6) = bfx2
        V_pt ( 7) = bfy2
        V_pt ( 8) = bfz2

        A_pt ( 1) = bfy2 * z
        A_pt ( 2) = bfz2 * x
        A_pt ( 3) = bfx2 * y

        END IF

! ----------------------------------------

#if( IPDE_SYSTEM == 2 )
        V_pt ( 9) = V_pt ( 5)
#endif

#if( IPDE_SYSTEM == 3 )
        V_pt ( 9) = V_pt ( 5)
        V_pt ( 10) = V_pt ( 5)
#endif

        ent = V_pt ( 5) / V_pt ( 1)**( gamma - 1.0)

        CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

! ----------------------------------------------------------------------
      ELSE IF ( iproblem == 7) THEN
! ----------------------------------------------------------------------
! Parker Solar wind problem.
! ntp_global = ntp_zones * ntp_patches = 20
!  nr_global =  nr_zones *  nr_patches = 30   ! And double with each resolution increase.
! rmin = 2.0; rmax = 3.5
! iproblem = 7
! timestop = 0.25
! set "bcarr_r" to 17 and call "ADD_BOUNDARY" after calling "PAD_BOUNDARY_R".
! switch on the accuracy calculation parts.

        r = SQRT ( x**2 + y**2 + z**2)
        d0 = 1.0
        u0 = 1.0
        p0 = 1.0

        gammaa = gamma
        gamma1 = gammaa - 1.0

! Solve for density: a(𝝆₀/𝝆)²+b(𝝆/𝝆₀)^𝜸₁+c=0
        dcur = 0.5 * d0 / r / r
        dprev = 0.0
        a = u0 * u0 / (2.0 * r * r * r * r)
        b = gammaa / gamma1 * p0 / d0
        c = -u0 * u0 / 2.0 - gammaa / gamma1 * p0 / d0

! Newton's method
        do while (abs((dprev - dcur) / dcur) > smallnum)
        f1 = a * (d0 / dcur)**2
        f2 = b * (dcur / d0)**gamma1
        f = f1 + f2 + c
        df = (-2.0 * f1 + gamma1 * f2) / dcur
        dprev = dcur
        dcur = dprev - f / df
        if(dcur <= 0.0) dcur = 0.01
        end do

! Density and velocity from mass conservation
        d  = dcur
        uu = d0 * u0 / d / r / r

! Adiabatic law for pressure
        p = p0 * (d / d0)**gammaa

! ----------

        rho = d
        vx = uu * x / r
        vy = uu * y / r
        vz = uu * z / r
        prs = p

        bfx = 0.0
        bfy = 0.0
        bfz = 0.0

! ----------

        V_pt ( 1) = rho

        V_pt ( 2) = vx
        V_pt ( 3) = vy
        V_pt ( 4) = vz

        V_pt ( 5) = prs

        V_pt ( 6) = bfx
        V_pt ( 7) = bfy
        V_pt ( 8) = bfz

        A_pt ( 1) = 0.0
        A_pt ( 2) = 0.0
        A_pt ( 3) = 0.0

! ----------------------------------------

#if( IPDE_SYSTEM == 2 )
        V_pt ( 9) = V_pt ( 5)
#endif

#if( IPDE_SYSTEM == 3 )
        V_pt ( 9) = V_pt ( 5)
        V_pt ( 10) = V_pt ( 5)
#endif

        ent = V_pt ( 5) / V_pt ( 1)**( gamma - 1.0)

        CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)


! ----------------------------------------------------------------------
      ELSE IF ( iproblem == 9) THEN
! ----------------------------------------------------------------------
! SOD test
! Lower Resolution:
!     ntp_global = ntp_zones * ntp_patches = 40
!      nr_global =  nr_zones *  nr_patches = 120
! Good  Resolution:
!     ntp_global = ntp_zones * ntp_patches = 120
!      nr_global =  nr_zones *  nr_patches = 240
! rmin = 1.5; rmax = 2.5
! iproblem = 9
! timestop = 0.22
! Inner and outer boundary conditions :- "bcarr_r ( :, :, 1) = 2, bcarr_r ( :, :, 2) = 2".
! with "ibdy_inner = 0, ibdy_outer = 0".

        rad = SQRT ( x**2 + y**2 + z**2)

        rho1 = 1.0
        prs1 = 1.0
        vlx1 = 0.0
        vly1 = 0.0
        vlz1 = 0.0

        rho2 = 0.125
        prs2 = 0.1
        vlx2 = 0.0
        vly2 = 0.0
        vlz2 = 0.0

#if ( CHOOSE_ADIABATIC == 1)
! keep the previously assigned pressure.

#else

        prs1 = isospeed**2 * rho1
     1      + barotropic_stiffness * rho1**gamma_euler

        prs2 = isospeed**2 * rho2
     1      + barotropic_stiffness * rho2**gamma_euler

#endif

        IF ( rad <= 0.5 * ( rmin + rmax) ) THEN

          rho = rho1
          vx = vlx1
          vy = vly1
          vz = vlz1
          prs = prs1

        ELSE

          rho = rho2
          vx = vlx2
          vy = vly2
          vz = vlz2
          prs = prs2

        END IF

        bfx = 0.0
        bfy = 0.0
        bfz = 0.0

! ----------

        V_pt ( 1) = rho

        V_pt ( 2) = vx
        V_pt ( 3) = vy
        V_pt ( 4) = vz

        V_pt ( 5) = prs

        V_pt ( 6) = bfx
        V_pt ( 7) = bfy
        V_pt ( 8) = bfz

        A_pt ( 1) = 0.0
        A_pt ( 2) = 0.0
        A_pt ( 3) = 0.0

! ----------------------------------------

#if( IPDE_SYSTEM == 2 )
        V_pt ( 9) = V_pt ( 5)
#endif

#if( IPDE_SYSTEM == 3 )
        V_pt ( 9) = V_pt ( 5)
        V_pt ( 10) = V_pt ( 5)
#endif

        ent = V_pt ( 5) / V_pt ( 1)**( gamma - 1.0)

        CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

! ----------------------------------------------------------------------
      ELSE IF ( iproblem == 10) THEN
! ----------------------------------------------------------------------
! Brio-Wu test
! Lower Resolution:
!     ntp_global = ntp_zones * ntp_patches = 120
!      nr_global =  nr_zones *  nr_patches = 150
! Good  Resolution:
!     ntp_global = ntp_zones * ntp_patches = 240
!      nr_global =  nr_zones *  nr_patches = 240
! rmin = 20.0; rmax = 21.0
! iproblem = 10
! timestop = 0.08
!
! 2nd way:
! ========
! ntp_zones = 50, nr_zones = 60,
! ntp_patches = 2, nr_patches = 4,
! rmin = 4.0, rmax = 5.0
!
! gamma = 2.0, gamma_euler = 2.0,
! molewt = 1.4, gasconst = 1.0,
! iproblem = 10
! timestop = 0.1
!
! IF ( iz_patch .EQ. 1         )  bcarr_r ( :, :, 1) = 17     ! Inner boundary.
! IF ( iz_patch .EQ. nr_patches)  bcarr_r ( :, :, 2) = 17     ! Outer boundary.
!
! IF ( iz_patch .EQ. 1         )  ibdy_inner = 1   ! 1 == Impose Add_boundary
! IF ( iz_patch .EQ. nr_patches)  ibdy_outer = 1   ! 1 == Impose Add_boundary
!
! --------------------

        rad = SQRT ( x**2 + y**2 + z**2)

        theta = ATAN2 ( SQRT ( x**2 + y**2), z)
        IF ( theta < 0.0) theta = theta + 2.0 * pi

        phi = ATAN2 ( y, x)
        IF ( phi < 0.0) phi = phi + 2.0 * pi

        costheta = COS ( theta)
        sintheta = SIN ( theta)

        cosphi = COS ( phi)
        sinphi = SIN ( phi)

        rho1 = 1.0
        prs1 = 1.0
        vlx1 = 0.0
        vly1 = 0.0
        vlz1 = 0.0

        rho2 = 0.125
        prs2 = 0.1
        vlx2 = 0.0
        vly2 = 0.0
        vlz2 = 0.0

        Br0 = 0.75 * SQRT ( 4.0 * pi)

        Br = Br0 * sintheta * ( 4.5 / rad)**2

        Bphi1 =  1.0 * SQRT ( 4.0 * pi)
        Bphi2 = -1.0 * SQRT ( 4.0 * pi)


        tempaa = 0.5 * ( rmin + rmax)

        IF ( (NINT(rmin) == 4) .AND. (NINT(rmax) == 6) ) tempaa = 4.5


        eps = 0.0025

        Bphi =
     1   Bphi1 * 0.5*sintheta * ( - TANH (( rad - tempaa) / eps) + 1.0)
     1 + Bphi2 * 0.5*sintheta * (   TANH (( rad - tempaa) / eps) + 1.0)

        Ar =
     1   0.5 * costheta * Bphi1 * rad
     1                        * ( - TANH (( rad - tempaa) / eps) + 1.0)
     1 + 0.5 * costheta * Bphi2 * rad
     1                        * (   TANH (( rad - tempaa) / eps) + 1.0)


        IF ( rad <= tempaa ) THEN

          rho = rho1
          vx  = vlx1
          vy  = vly1
          vz  = vlz1
          prs = prs1

        ELSE

          rho = rho2
          vx  = vlx2
          vy  = vly2
          vz  = vlz2
          prs = prs2

        END IF

        bfx = ( Br * x / rad) - ( Bphi * y / SQRT ( x**2 + y**2))
        bfy = ( Br * y / rad) + ( Bphi * x / SQRT ( x**2 + y**2))
        bfz = Br * z / rad

        A_pt ( 1) = Ar * x / rad
        A_pt ( 2) = Ar * y / rad
        A_pt ( 3) = Ar * z / rad


! ----------
#if ( CHOOSE_ADIABATIC == 1)
! keep the previously assigned pressure.
#else
        prs = isospeed**2 * rho
     1      + barotropic_stiffness * rho**gamma_euler
#endif

        V_pt ( 1) = rho

        V_pt ( 2) = vx
        V_pt ( 3) = vy
        V_pt ( 4) = vz

        V_pt ( 5) = prs

        V_pt ( 6) = bfx
        V_pt ( 7) = bfy
        V_pt ( 8) = bfz

! ----------------------------------------

#if( IPDE_SYSTEM == 2 )
        V_pt ( 9) = V_pt ( 5)
#endif

#if( IPDE_SYSTEM == 3 )
        V_pt ( 9) = V_pt ( 5)
        V_pt ( 10) = V_pt ( 5)
#endif

        ent = V_pt ( 5) / V_pt ( 1)**( gamma - 1.0)

        CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

! ----------------------------------------------------------------------
      ELSE IF(iproblem == 13) THEN
! ----------------------------------------------------------------------
! VF: Rotating tilted split monopole with a band of slow solar
! wind near the magnetic equator. This configuration can produce
! corotating interation regions (CIRs).
!
! This problem is setup by Prof. Florinski. It is self-contained with
! the rotating split monopole. So, the following two 'omega_' variables
! of Geomesh-Code should be set to zero in "Riemann_Module_Main_CS_src".
!
! omega_passive ( :) = 0.0
! omega_active ( :) = 0.0
!
! rmin = 1.0
! rmax = 10.0
! iproblem = 13
! timestop = 1.0
!
! ntp_global = ntp_zones * ntp_patches = 100
!  nr_global =  nr_zones *  nr_patches = 100
!
! VF: Rotating tilted split monopole with a band of slow solar
! wind near the magnetic equator. This configuration can produce
! corotating interation regions (CIRs).

        r0 = rmin

! In the slow wind speed and number density = 1.
! The other parameters can be derived from:
        fast_slow = 2.0
        mach = 4.3
        beta = 1.4

! The rotation axis is assumed in the z direction. The magnetic
! axis is tilted by "tilt", initially in the -y direction.
! An equatorial band of width 2*band_width in the magnetic frame is
! occupied by slow wind and the polar regions by the fast wind
        tilt = 20.0 * pi / 180.0
        band_width = 30.0 * pi / 180.0
        omega = 1.0
        trans_width = 1.0 * pi / 180.0

        costilt = cos(tilt)
        sintilt = sin(tilt)
        cosomegat = cos(omega * time)
        sinomegat = sin(omega * time)

        rad = sqrt(x**2 + y**2 + z**2)

! Inertial frame
        theta = acos(z / rad)
        phi = atan2(y, x)
        IF(phi < 0.0) phi = phi + 2.0 * pi

        costheta = cos(theta)
        sintheta = sin(theta)

        cosphi = cos(phi)
        sinphi = sin(phi)

! Magnetic frame
        x_mag = x * cosomegat + y * sinomegat

        y_mag =-x * costilt * sinomegat
     1        + y * costilt * cosomegat
     1        + z * sintilt

        z_mag = x * sintilt * sinomegat
     1        - y * sintilt * cosomegat
     1        + z * costilt

        theta_mag = acos(z_mag / rad)
        phi_mag = atan2(y_mag, x_mag)
        IF(phi_mag < 0.0) phi_mag = phi_mag + 2.0 * pi

        costheta_mag = cos(theta_mag)
        sintheta_mag = sin(theta_mag)

        cosphi_mag = cos(phi_mag)
        sinphi_mag = sin(phi_mag)

! Velocity
        IF(ABS(0.5 * pi - theta_mag) < band_width) THEN
           vr = 1.0
        ELSE
           vr = fast_slow
        ENDIF

! Mass flux and pressure are equal in slow and fast wind
        rho0 = 1.0 / vr
        prs0 = 1.0 / (gamma * mach**2)

! Azimuthal speed - use corotating speed inside the star!
        vphi = omega * r0 * sintheta
        if(rad < r0) vphi = vphi / r0 * rad

! Flow properties in the inertial frame
        rho = rho0 * (r0 / rad)**2
        vx = vr * sintheta * cosphi - vphi * sinphi
        vy = vr * sintheta * sinphi + vphi * cosphi
        vz = vr * costheta
        prs = prs0 * ((rho / rho0)**gamma)

! Magnetic field is a split monopole. Its vector potential is
! singularity free if we take care of the sine function in the
! denominator.
        B_0 = sqrt(8.0 * pi * prs0 / beta)

        Br = B_0 * (r0 / rad)**2 * sign(1.0, 0.5 * pi - theta_mag)
     1       * (1.0 - exp(-abs(0.5 * pi - theta_mag)
     1         / trans_width) / max(sintheta_mag, smallnum))

        Aphi = B_0 * r0**2 / rad * (1.0 - abs(costheta_mag)
     1         - trans_width * exp(-abs(0.5 * pi - theta_mag)
     1         / trans_width)) / max(sintheta_mag, smallnum)

        bfx = Br * x / rad
        bfy = Br * y / rad
        bfz = Br * z / rad

! Conservative variables
        V_pt(1) = rho

        V_pt(2) = vx
        V_pt(3) = vy
        V_pt(4) = vz

        V_pt(5) = prs

        V_pt(6) = bfx
        V_pt(7) = bfy
        V_pt(8) = bfz

        A_pt(1) = Aphi * (-cosomegat * sinphi_mag
     1                       - costilt * sinomegat * cosphi_mag)
        A_pt(2) = Aphi * (-sinomegat * sinphi_mag
     1                       + costilt * cosomegat * cosphi_mag)
        A_pt(3) = Aphi * sintilt * cosphi_mag

! ----------------------------------------

#if( IPDE_SYSTEM == 2 )
        V_pt ( 9) = V_pt ( 5)
#endif

#if( IPDE_SYSTEM == 3 )
        V_pt ( 9) = V_pt ( 5)
        V_pt ( 10) = V_pt ( 5)
#endif

        ent = V_pt ( 5) / V_pt ( 1)**( gamma - 1.0)

        CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)


! ----------------------------------------------------------------------
! iproblem 14: Heliosphere simulation (Swati Sharma and V. Florinski)

#if( CHOOSE_APPL_HELIO == 1)
! ----------------------------------------------------------------------
      ELSE IF ( iproblem == 14) THEN
! ----------------------------------------------------------------------

      CALL HELIO_INITIALIZE

     1                  ( n_cc_components,
     1                    x, y, z, time, iproblem,
     1                    rmin, rmax, zone_size,
     1                    xcen_sector, ycen_sector, zcen_sector,
     1                    xvert1_sector, yvert1_sector, zvert1_sector,

     1                    U_pt, A_pt )

#endif

! ----------------------------------------------------------------------
      ELSE IF ( iproblem == 16) THEN
! ----------------------------------------------------------------------
#if( CHOOSE_APPL_CAK == 1)
! --------------------------------------------------
! Stellar wind - CAK

        xs = x
        ys = y
        zs = z
        Rs = R

#define CAK_PREV 1

#if( CAK_PREV == 1)
        IF ( R < rmin ) THEN
          xs = ( x / R) * rmin
          ys = ( y / R) * rmin
          zs = ( z / R) * rmin
          Rs = rmin
        END IF
#endif

        CALL CAK_INIT_1D ( xs, ys, zs, V_CAK, cak_rhoc)

        rho = V_CAK ( 1)
        vx  = V_CAK ( 2) * xs / Rs
        vy  = V_CAK ( 2) * ys / Rs
        vz  = V_CAK ( 2) * zs / Rs
        prs = V_CAK ( 3)

#if( CAK_PREV == 1)
        IF ( R < rmin ) THEN
          rho = cak_rhoc
        END IF
#endif

! --------------------

! init the dipolar m-field.
        m1 = dipole_B ( 1)
        m2 = dipole_B ( 2)
        m3 = dipole_B ( 3)

#if( USE_B0B1_FORMULATION == 1)
        m1 = 0.0
        m2 = 0.0
        m3 = 0.0
#endif

!       x = xs; y = ys; z = zs; R = Rs;

        bfx  = - ( m1/R**3) + (3*x*(m1*x + m2*y + m3*z)) / R**5
        bfy  = - ( m2/R**3) + (3*y*(m1*x + m2*y + m3*z)) / R**5
        bfz  = - ( m3/R**3) + (3*z*(m1*x + m2*y + m3*z)) / R**5

! Initializa vector potential Ax, Ay, Az
        A_pt ( 1) = (m2*z - m3*y) / R**3
        A_pt ( 2) = (m3*x - m1*z) / R**3
        A_pt ( 3) = (m1*y - m2*x) / R**3

! --------------------

        V_pt ( 1) = rho

        V_pt ( 2) = vx
        V_pt ( 3) = vy
        V_pt ( 4) = vz

        V_pt ( 5) = prs

        V_pt ( 6) = bfx
        V_pt ( 7) = bfy
        V_pt ( 8) = bfz

! ----------------------------------------

#if( IPDE_SYSTEM == 2 )
        V_pt ( 9) = V_pt ( 5)
#endif

#if( IPDE_SYSTEM == 3 )
        V_pt ( 9) = V_pt ( 5)
        V_pt ( 10) = V_pt ( 5)
#endif

        ent = V_pt ( 5) / V_pt ( 1)**( gamma - 1.0)

        CALL PRIM_TO_CONS (
     1       gamma, gamma_euler, isospeed, barotropic_stiffness,
     1       molewt, gasconst, rhofloor, prsfloor, tempfloor,
     1       magB_small_cgl, eps_bdry_cgl,
     1       iprs_flag, ent, B0, U_pt, V_pt)

! --------------------------------------------------
#endif
! if( CHOOSE_APPL_CAK == 1)
! ----------------------------------------------------------------------
      END IF  ! End of "IF ( iproblem == 1)".
! ----------------------------------------------------------------------

      END SUBROUTINE PROBLEM_GENERATOR_U_A



