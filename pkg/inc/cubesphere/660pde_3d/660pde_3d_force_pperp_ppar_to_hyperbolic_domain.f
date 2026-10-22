! ----------------------------------------------------------------------
! This is new PDE_3D_src that works.
! ----------------------------------------------------------------------

! Must set "NVAR, IPDE_SYSTEM" for this to work.




      SUBROUTINE FORCE_Pperp_Ppar_To_HYPERBOLIC_DOMAIN
     1           ( pi, magB_small_cgl, eps_bdry_cgl, prsfloor,
     1             Bsq, p_par, p_perp,
     1             updated_p_par, updated_p_perp )

! ----------------------------------------------------------------------

! Written by : Dinshaw Balsara, Deepak Bhoriya

! This subroutine forces to p_perp and p_par to remain in the physical domain.
! Only applicable to CGL system.
! Inputs: pi, Bsq, p_par, p_perp    Outputs: updated_p_par, updated_p_perp

! ----------------------------------------------------------------------

      IMPLICIT NONE

      REAL pi, Bsq, p_par, p_perp, updated_p_par, updated_p_perp

! ----------------------------------------------------------------------

! Variables past this point are internal to this subroutine.

      REAL beta_par, beta_perp, beta_bar, eps, tempaa, a, b, c,
     1     magB_small_cgl, eps_bdry_cgl, prsfloor, p_bar,
     1     beta_par_firehose, beta_par_mirror, beta_perp_temp,
     1     beta_par_firehose_L, beta_par_firehose_R,
     1     beta_perp_mirror,
     1     beta_par_mirror_orig, beta_perp_mirror_orig,
     1     beta_perp_mirror_orig_m_eps, eps_fix

! ----------------------------------------------------------------------

      p_bar = (p_par + 2.0 * p_perp) / 3.0

      eps = eps_bdry_cgl

      IF ( p_bar .LE. prsfloor) THEN
! Take lots of corrective action if we get here.

      p_par = prsfloor
      p_perp = prsfloor
      p_bar = prsfloor

      END IF

! Return with Average values if Magnetic Field is 0.

      IF ( SQRT(Bsq) .LT. magB_small_cgl) THEN

         updated_p_par = p_bar
         updated_p_perp = p_bar   ! Updated isotropic & floored pressures returned.
         RETURN

      END IF

! ------

      beta_par = 8.0 * pi * p_par /  ( Bsq + 1.0e-32)
      beta_perp = 8.0 * pi * p_perp / ( Bsq + 1.0e-32)

      beta_bar = ( beta_par + 2.0 * beta_perp) / 3.0

! ------
! Here we get "beta_par_mirror_orig" obtained from the quadratic without shifting by epsilon.
      beta_par_mirror_orig = 2.0 * (beta_bar + 1.0/3.0)
     1               - SQRT ( beta_bar**2 + 8.0*beta_bar/3.0 + 4.0/9.0)

      beta_perp_mirror_orig = 0.5 * ( 3.0 * beta_bar
     1                              - beta_par_mirror_orig)

      beta_perp_mirror_orig_m_eps = beta_perp_mirror_orig - eps

! ------
! Here we get "beta_par_mirror" obtained from the quadratic which is shifted by epsilon.

      eps_fix = ((-6.0 + 4.0*SQRT(6.0) - 3.0*eps)*eps)
     1        / (2.0 + 2.0*(-3.0 + SQRT(6.0))*eps)

      a = (3.0 - 2.0*eps_fix)
      b = (1.0 - eps_fix)*(2.0 - 3.0*beta_bar)
      c = -3.0*(1.0 - eps_fix)*beta_bar

      beta_perp_temp = ( - b + SQRT (b**2 - 4.0*a*c) ) / (2.0*a)

      ! ------

      beta_perp_mirror =  AMAX1 ( beta_perp_temp,
     1                            beta_perp_mirror_orig_m_eps)

      beta_par_mirror = 3.0 * beta_bar - 2.0 * beta_perp_mirror

! ------
! Here we get "beta_par_firehose"

      beta_par_firehose_L = 3.0 * beta_bar / (1.0 + eps)
      beta_par_firehose_R = beta_bar + (2.0/3.0) * (2.0 - eps)

      beta_par_firehose = AMIN1 ( beta_par_firehose_L,
     1                            beta_par_firehose_R)

! ------
! Use "beta_par_mirror" and "beta_par_firehose" to bound the actual "beta_par".

      beta_par = AMAX1 ( AMIN1 ( beta_par, beta_par_firehose),
     1                           beta_par_mirror)

      beta_perp = 1.5 * ( beta_bar - beta_par / 3.0)

! ------
! Updating the p_par and p_perp in case they have changed.

      updated_p_par  = beta_par * Bsq / ( 8.0 * pi )
      updated_p_perp = beta_perp * Bsq / ( 8.0 * pi )

! ----------------------------------------------------------------------

      END SUBROUTINE FORCE_Pperp_Ppar_To_HYPERBOLIC_DOMAIN


