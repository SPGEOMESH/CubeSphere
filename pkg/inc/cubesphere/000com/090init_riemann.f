





      SUBROUTINE INIT_RIEMANN ()

#include <cubesphere/000com/riemann.com>


      smallnum = 1.0e-12
      pi = 4.0 * ATAN ( 1.0)

      rhofloor = 1.0e-12   ! Density floor.
      prsfloor = 1.0e-9    ! Pressure floor.


      gamma = 1.667     ! Polytropic index for the gas.
      gasconst = 1.0  ! Set up scaled gas constant.
      molewt = 1.4    ! Set up mole wt. of the gas.

      gamma_euler = 1.4           ! Polytropic index for the gas.
      gasconst_euler = 1.0        ! Set up scaled gas constant.
      molewt_euler = 1.4          ! Set up mole wt. of the gas.

      isospeed = 1.0              ! Isothermal sound speed; only for isomhd problems.
      barotropic_stiffness = 0.0  ! Stiffens the isothermal EOS if non-zero.

      tempfloor = 1.0e-5          ! used only when the cooling fn is called.

      source_tau = 1.0e-10        ! CGL relaxation time


#if( IPDE_SYSTEM >= 2)
      gamma = 5.0 / 3.0
      gamma_euler = 5.0 / 3.0
#endif

! Note: The following numbers to be changed for suitable problems.
! If "Mstar" is non-zero then gravity will be switched on in the ADER subroutine.
! If "omega_passive" is non-zero then non-inertial terms will be switched
! on in the ADER subroutine.

      newtonsc = 6.674E-08        ! Newton's constant; in cgs unit
      Mstar = 0.0                 ! Mass of Sun; in grams  =  1.989E33
      omega_passive ( :) = 0.0    ! Rotation axis of star.
      omega_passive ( 3) = 0.0

      omega_active ( :) = 0.0

! ----------
! For multipole simulations.

      iflag_npoles ( 1:4) = 0
      iflag_npoles ( 1)   = 1     ! Enable only dipole by default.

      multipole_Bmag ( 1:4) = 0.0
      multipole_Bmag ( 1)   = 0.0

! ----------
! Define the dipole and qudrapole magnetic field strengths.

      a = pi / 6.0

! monopole.
      iflag_monopole  = 0    ! 0 - no monopole field, 1 - monpole, 2 - split-monopole.
      monopole_u0     = 0.0  ! radial velocity, used in parker spiral solution
      monopole_B ( :) = 0.0

! dipole and quadrupole.
      dipole_B ( 1) = 0.0 * SIN ( a)
      dipole_B ( 2) = 0.0
      dipole_B ( 3) = 0.0 * COS ( a)

      qpole_B ( 1:3, 1:3) = 0.0


! Courant number
      cfl_coef = 0.2

! If the divergence of the velocity at a zone center exceeds
! this coef times the sound speed then we have to flatten the slopes.
      flatten_coef = 0.4

! This is the factor by which the pressure on one or the other side
! of a zone has to be larger. If this condition is met, we look to
! that side to increase the flattener in the zone being considered.
      prs_flatten_coef = 1.5

      fracdiff = 0.1   ! Set fractional change needed for entropy fix.(deprecated)

! If the divergence of the velocity at a zone boundary exceeds this coef times
! the sound speed then we have to blend in some fraction of HLL RS.
      divratio_einfeldt_rs = 0.4

! Set a threshold value that triggers the Lapidus viscosity as well as
! its value. Set "lapidus_coef" between [0.0, 0.2]. Default is 0.0.
      lapidus_threshold = 0.2
      lapidus_coef = 0.05



! Using 100% multid flux can sometimes mean we take flux contributions from
! a stencil that does not depend on the zones that contribute to that face.
! For that reason, we blend multid RS with 1d RS associated with that face.
! Blending in more multidimensionality also seems to help out with pressure
! positivity without any artificial fixes.
! This value should range between [0.4, 0.6] for good performance.

! "blend_multid_flux = 1.0 / 6.0" means we use Simpson rule with 1/6 of the
! flux from each edge.
! "blend_multid_flux = 0.0" means we use full flux from suitable 1d RS.
! The electric field is always fully evaluated by the multid. RS.

      blend_multid_flux = 0.1666

! In the fully subsonic case (just) the multi-d HLLC RS can oversteepen
! (almost or fully) stationary contacts. To prevent that we blend in a small
! fraction of the multi-d HLL RS. Ranging between [0.1, 0.2]

! While the flattener also plays a role; "blend_mdhll_to_mdhllc = 1.0"
! converts HLLC into a pure HLL. "blend_mdhll_to_mdhllc = 0.0" does the
! reverse in smooth flow.

      blend_mdhll_to_mdhllc = 0.0


      ! todo initialize other quantities as needed


      END SUBROUTINE INIT_RIEMANN


