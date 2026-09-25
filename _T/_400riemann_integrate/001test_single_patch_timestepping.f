

#define DEBUG_SILO 0











      PROGRAM TEST_SINGLE_PATCH_TIMESTEPPING


      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

      INTEGER ntp_global, nr_global, NGhost_global,
     1        ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1        nRK_stage, irk, iproblem, ndumpstep, nimgstep,

     1        isector, ix_patch, iy_patch, iz_patch,
     1        ilevel, ilevel_max, iscale,

     1        n_cc_modes,
     1        n_cc_components, n_fc_modes,

     1        itri_patch_LD, itri_patch_RD, itri_patch_RU,
     1        itri_patch_LU,

     1        istep, ntstep

      PARAMETER ( ntp_zones = 80, nr_zones = 80,
     1            ntp_patches = 1, nr_patches = 1, NGhost = 4,

     1            ilevel = 1, ilevel_max = 1,
     1            iscale = 2**(ilevel_max - ilevel + 1) / 2,

     1            ntp_global = ntp_zones * ntp_patches * iscale,
     1            nr_global = nr_zones * nr_patches * iscale,
     1            NGhost_global = NGhost * iscale,

     1            n_cc_components = NVAR + NFLUID,
     1            n_fc_modes = 1,

#if ( NUMORDER == 2)

     1            n_cc_modes = 4,

#elif ( NUMORDER == 3)

     1            n_cc_modes = 10,

#elif ( NUMORDER == 4)

     1            n_cc_modes = 20,

#endif

     1            )

      REAL rmin, rmax, time, dtcur, dtcurp1, dtmove, a,
     1     timeimage, dtimage, timestop,

     1     err_U ( n_cc_components), err_Bxf, err_Byf, err_Bzf

!$acc declare create ( time, dtcur, dtcurp1, dtmove )

      REAL sector_vertices ( 4, 3, 6)

      INTEGER ixx, iyy, izz, i,
     1        ibdy_inner, ibdy_outer

! ----------------------------------------------------------------------

      REAL coordinates_gr_global
     1     ( 0 - NGhost_global: ntp_global + NGhost_global,
     1       0 - NGhost_global: ntp_global + NGhost_global, 3, 6)

      REAL vert_coords_patch
     1 ( 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost, 3)

      INTEGER imask_active_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                          1 - NGhost: ntp_zones + NGhost,
     1                          1 - NGhost: nr_zones + NGhost)

      REAL rb_save_global ( 0 - NGhost_global:
     1                      nr_global + NGhost_global)

      REAL, DIMENSION ( 1  - NGhost_global: nr_global + NGhost_global)::
     1     dr_save_global, rc_save_global

! ----------------------------------------------------------------------

      REAL xf_quadpts_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL xf_quadwt_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NAREA_QUAD_PTS)

      REAL xf_quadnorml_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL xf_area_patch ( 0 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost),
     1     xf_centroid_patch ( 0 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

! ----------

      REAL yf_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL yf_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NAREA_QUAD_PTS)

      REAL yf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL yf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost),
     1     yf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost, 3)

! ----------

      REAL zf_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL zf_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   NAREA_QUAD_PTS)

      REAL zf_quadnorml_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   3, NAREA_QUAD_PTS)

      REAL zf_area_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost),
     1     zf_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost, 3)

! ------------------------------

      REAL v_scale_length ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

! ------------------------------

      REAL v_quadpts_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, NUMORDER * NUMORDER * NUMORDER)

      REAL v_quadwt_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   NUMORDER * NUMORDER * NUMORDER)

      REAL v_quadtangent_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3, 3, NUMORDER * NUMORDER * NUMORDER)

      REAL v_volume_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL v_centroid_patch ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

! ------------------------------

      INTEGER, DIMENSION ( n_cc_modes) ::
     1         ii_mode, jj_mode, kk_mode

      REAL Coef_Basis_Moments ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_modes)

! ------------------------------

      INTEGER istencil_patch ( 1 - NGhost: ntp_zones + NGhost,
     1                         1 - NGhost: ntp_zones + NGhost, 2, 12)

! ----------------------------------------------------------------------
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

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           n_fc_modes) ::
     1   Bxf, Byf, Bzf

      REAL B0field ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   3)

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_bfld_modes) ::
     1   B1_zc, B2_zc, B3_zc

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost, n_cc_modes) ::
     1   Ex_zc, Ey_zc, Ez_zc

      REAL Source ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components)

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: nr_zones + NGhost) ::
     1      gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1      divvelgr, msonicgr, flattengr, workgr

      REAL, DIMENSION ( n_cc_components + 1) :: V_min_local, V_max_local

! ------------------------------
! Variables used for div-free reconstruction of magnetic fields:-

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1           NBFLD_QUAD_PTS) ::
     1   Bxf_quadpts, Byf_quadpts, Bzf_quadpts

! ------------------------------
! VARIABLES NEEDED FOR 2 STAGE RUNGE-KUTTA TIMESTEPPING:-

      REAL U0 ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost,
     1   n_cc_components)

      REAL entropygr0 ( 1 - NGhost: ntp_zones + NGhost,
     1   1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost) ::
     1   Bxf0, Byf0, Bzf0

! ------------------------------
! FLUX AND ELECTRIC FIELD VARIABLES:-

! Extremal wave speeds in the plus and minus directions at each zone boundary.
! These are obtained from the Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      Sxpf, Sxmf, Sypf, Symf, Szpf, Szmf

! Store extremal speeds from the 1D Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: nr_zones + NGhost) ::
     1   S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr

! Fluid fluxes and resolved states from the 1D Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1   0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost,
     1   n_cc_components) ::
     1      UStar_x, UStar_y, UStar_z

! Entropy flux from the 1D Riemann solver.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      EntFlux_x, EntFlux_y, EntFlux_z

! Electric fields at the edges of the mesh.

      REAL, DIMENSION ( 0 - NGhost: ntp_zones + NGhost,
     1 0 - NGhost: ntp_zones + NGhost, 0 - NGhost: nr_zones + NGhost)::
     1      Efx, Efy, Efz

! Extra variables for silo output; entropy, tempfloor flags.

      REAL, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost: ntp_zones + NGhost,
     1                  1 - NGhost:  nr_zones + NGhost) ::
     1   flag_entprs, flag_Tfloor

! ----------------------------------------------------------------------
! File names used for imaging.

      CHARACTER*8 rhofil, vlxfil, vlyfil, vlzfil, prsfil,
     1            bfxfil, bfyfil, bfzfil, dumpfil,
     1            speciefil

! Real*4 variable used for imaging variables.

      REAL*4, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1 1 - NGhost: ntp_zones + NGhost, 1 - NGhost: nr_zones + NGhost)::
     1   r4gr

! Enforcing radial boundary conditions. Used for PAD_BOUNDARY subroutine.

      INTEGER, DIMENSION ( 1 - NGhost: ntp_zones + NGhost,
     1                     1 - NGhost: ntp_zones + NGhost, 2) :: bcarr_r

! ----------------------------------------------------------------------

      REAL timing_info1, timing_info2

#if ( DEBUG_SILO == 1)
      INTEGER itri_patch ( 4)
      INTEGER top_or_bottom ( 2)   ! Indicate top block or bottom block or both.

#if( USE_B0B1_FORMULATION == 1 )
      CHARACTER*SILO_VARLEN silo_allvars ( NVAR + NFLUID + 3)
#else
      CHARACTER*SILO_VARLEN silo_allvars ( NVAR + NFLUID)
#endif

#endif

! ----------------------------------------------------------------------
! Create arrays that only need to reside on the GPU. These arrays will
! never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

! Also create fluid variable arrays that are needed by the algorithm and only need
! to reside on the GPU. These arrays are not directly used for the update
! and will never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

!$acc data create(xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,&
!$acc&           xf_area_patch, xf_centroid_patch,                      &

!$acc&           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch, &
!$acc&           yf_area_patch, yf_centroid_patch,                      &

!$acc&           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch, &
!$acc&           zf_area_patch, zf_centroid_patch,                      &

!$acc&           v_scale_length,                                        &

!$acc&           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,  &
!$acc&           v_volume_patch, v_centroid_patch,                      &

!$acc&           Coef_Basis_Moments,                                    &
!$acc&           itri_patch_LD, itri_patch_RD, itri_patch_RU,           &
!$acc&           itri_patch_LU, istencil_patch,                         &

!$acc&         U, entropygr, Bxf, Byf, Bzf,                             &
!$acc&         V, B1_zc, B2_zc, B3_zc, Source,                          &

!$acc&         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,       &
!$acc&         divvelgr, msonicgr, flattengr, workgr,                   &
!$acc&         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                   &

!$acc&         V_min_local, V_max_local,                                &

!$acc&         Sxpf, Sxmf, Sypf, Symf, Szpf, Szmf,                      &
!$acc&         UStar_x, UStar_y, UStar_z,                               &
!$acc&         EntFlux_x, EntFlux_y, EntFlux_z,                         &
!$acc&         Efx, Efy, Efz, flag_entprs, flag_Tfloor,                 &

!$acc&     U0, entropygr0, Bxf0, Byf0, Bzf0, bcarr_r,                   &
!$acc&     ibdy_inner, ibdy_outer, iproblem, rmin, rmax )

! ----------------------------------------------------------------------
! Build a 1D mesh in the radial direction.

      OPEN ( UNIT = 6, FILE = "riemann_cs.out")

      rhofloor = 1.0e-12   ! Density floor.
      prsfloor = 1.0e-8    ! Pressure floor.
      smallnum = 1.0e-12   ! Effective small number used in code.

      pi = 4.0 * ATAN ( 1.0)

      ix_patch = 1
      iy_patch = 1
      iz_patch = 1

      isector = 3

      rmin = 7.0  ! 14.0 - 1.5  ! for 5x5 ntp_patches and midpatch
      rmax = 10.0 ! 14.0 + 1.5

      gamma = 1.4     ! Polytropic index for the gas.
      gasconst = 1.0  ! Set up scaled gas constant.
      molewt = 1.4    ! Set up mole wt. of the gas.

      gamma_euler = 1.4     ! Polytropic index for the gas.
      gasconst_euler = 1.0  ! Set up scaled gas constant.
      molewt_euler = 1.4    ! Set up mole wt. of the gas.

      time = 0.0  ! Initial time is always zero.

! Decide which physical problem you want to initialize with "iproblem".
! iproblem = 1 for method of manufactored solution; iproblem = 2 is blast;
! iproblem = 3 is hydro vortex; iproblem = 4 is mhd vortex;
! iproblem = 5 is accuracy test; iproblem = 6 is rotor test;
! iproblem = 7 is Parker wind test; iproblem = 8 hydrodynamic manufactored solution;
! iproblem = 9 is SOD test; iproblem = 10 is Brio-Wu test.

      iproblem = 2

!$acc update device( iproblem, rmin, rmax )

! ----------------------------------------------------------------------

! These parameters control the numerics.

#if ( NUMORDER == 2)
      nRK_stage = 2  ! Two stage RK at second order.
#elif ( NUMORDER == 3)
      nRK_stage = 3  ! Three stage RK at third order.
#elif ( NUMORDER == 4)
      nRK_stage = 5  ! Four stage RK at fourth order.
#endif

! Courant number
      cfl_coef = 0.25

! If the divergence of the velocity at a zone center exceeds
! this coef times the sound speed then we have to flatten the slopes.
      flatten_coef = 0.4

! This is the factor by which the pressure on one or the other side
! of a zone has to be larger. If this condition is met, we look to
! that side to increase the flattener in the zone being considered.
      prs_flatten_coef = 1.5

      fracdiff = 0.1  !Set fractional change needed for entropy fix.(deprecated)

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

! ----------------------------------------------------------------------

! These parameters control the checkpoint dumps and imaging.

! Start with initial time set to 0.Restart will refresh this from "riemann.in".
      time = 0.0

! Time at which we want to output the next image.Set this yourself for restart.
      timeimage = 0.0

! Number of timesteps you want to take. Code stops after these many steps.
      ntstep = 25 ! 500000

! Write out full restart dumps every "ndumpstep" steps. By default, all the
! variables are saved in 64 bit unformatted files called "dump0001" etc.
      ndumpstep = 500

! Write out images at every "nimgstep" timesteps.(not used)
      nimgstep = 10

! Write out images every "dtimage" amount of time. By default, each variable
! is saved in its own file name in 32 bit formatted form, like "rhoa0001" etc.
      dtimage = 0.13 / 10.0

! Time at which you want the simulation to stop.
      timestop = 0.13

! ----------------------------------------------------------------------
! By this time we have updated a lot of variables on the CPU. Let us
! migrate that data to the GPU. We still want to continue working
! on CPU, so we set "do_on_gpu = .FALSE.".

      do_on_gpu = .FALSE.

!$acc update device  (/RIEMANN_COM/)

! ----------------------------------------------------------------------

! Set default multispecies values, just to make sure they are not NANS.

      gamma_specie ( :) = gamma
      molewt_specie ( :) = molewt
      heat_of_formation ( :) = 0.0
      cp_specie ( :) = gasconst * gamma_specie ( :)
     1            / ( ( gamma_specie ( :) - 1.0) * molewt_specie ( :) )
      cv_specie ( :) = cp_specie ( :) / gamma_specie ( :)

! ----------

      gamma_specie_euler ( :) = gamma_euler
      molewt_specie_euler ( :) = molewt_euler
      heat_of_formation_euler ( :) = 0.0
      cp_specie_euler ( :) = gasconst_euler * gamma_specie_euler ( :)
     1 / ( ( gamma_specie_euler ( :) - 1.0) * molewt_specie_euler ( :))
      cv_specie_euler ( :) = cp_specie_euler ( :)
     1                     / gamma_specie_euler ( :)

! ----------
! Make a default initialization of big fluid arrays. This is still done
! on CPU.

      gamma_eos_gr ( :, :, :) = gamma
      gamma_soundspeed_gr ( :, :, :) = gamma
      mean_mol_wt_gr ( :, :, :) = molewt

! ----------------------------------------------------------------------

! Set the variable names used in the silo file.

#if ( DEBUG_SILO == 1)
      silo_allvars ( 1) = "Density "
      silo_allvars ( 2) = "XMomentum "
      silo_allvars ( 3) = "YMomentum "
      silo_allvars ( 4) = "ZMomentum "
      silo_allvars ( 5) = "Energy "
      silo_allvars ( 6) = "XMagneticField "
      silo_allvars ( 7) = "YMagneticField "
      silo_allvars ( 8) = "ZMagneticField "

#if( USE_B0B1_FORMULATION == 1)
      silo_allvars ( 9) = "XTotalMagnField "
      silo_allvars (10) = "YTotalMagnField "
      silo_allvars (11) = "ZTotalMagnField "
#endif
#endif

! Set the variable names used in the silo file. Not done yet.

! Set file names for the image files. Note that this code just dumps
! out zone-centered variables for imaging. Can't restart from these files.

      rhofil = "rhoa0001"
      vlxfil = "vlxa0001"
      vlyfil = "vlya0001"
      vlzfil = "vlza0001"
      prsfil = "prsa0001"

      bfxfil = "bfxa0001"
      bfyfil = "bfya0001"
      bfzfil = "bfza0001"

      speciefil = "sp010001"

! Set file name for the first dump file.

      dumpfil = "dump0001"

! ----------------------------------------------------------------------

      IF ( ix_patch < 1 ) ix_patch = 1
      IF ( iy_patch < 1 ) iy_patch = 1
      IF ( iz_patch < 1 ) iz_patch = 1

      IF ( ix_patch > ntp_patches ) ix_patch = ntp_patches
      IF ( iy_patch > ntp_patches ) iy_patch = ntp_patches
      IF ( iz_patch >  nr_patches ) iz_patch =  nr_patches


! ----------------------------------------------------------------------

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Starting problem with the following parameters"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Chose physical problem given by iproblem = ",
     1              iproblem
      WRITE ( 6, *)" "


      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1"Starting timestep; Documenting constants in the module"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Order of accuracy = ", NUMORDER
      WRITE ( 6, *)"Number of zones theta-phi, Number of ghost cells= ",
     1             ntp_zones, NGhost
      WRITE ( 6, *)"Number of zones in radial direction = ", nr_zones
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Picking one of the patches from the sector"
      WRITE ( 6, *)"isector = ", isector
      WRITE ( 6, *)"ix_patch, iy_patch, iz_patch = ",
     1              ix_patch, iy_patch, iz_patch
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)
     1       "n_cc_components, n_cc_modes = ",
     1        n_cc_components, n_cc_modes
      WRITE ( 6, *)"n_fc_modes = ", n_fc_modes
      WRITE ( 6, *)" "

      WRITE ( 6, *)" "
      WRITE ( 6, *)"pi, rhofloor, prsfloor, smallnum = ",
     1              pi, rhofloor, prsfloor, smallnum
      WRITE ( 6, *)"gamma, gasconst, molewt = ",gamma, gasconst, molewt
      WRITE ( 6, *)" "
      WRITE ( 6, *)"cfl_coef = ", cfl_coef
      WRITE ( 6, *)" "

      WRITE ( 6, *)"flatten_coef, divratio_einfeldt_rs = ",
     1              flatten_coef, divratio_einfeldt_rs
      WRITE ( 6, *)" "
      WRITE ( 6, *)"lapidus_threshold, lapidus_coef = ",
     1              lapidus_threshold, lapidus_coef
      WRITE ( 6, *)" "
      WRITE ( 6, *)"blend_multid_flux, blend_mdhll_to_mdhllc = ",
     1              blend_multid_flux, blend_mdhll_to_mdhllc

      WRITE ( 6, *)"ntstep, ndumpstep, nimgstep, time, timestop = ",
     1              ntstep, ndumpstep, nimgstep, time, timestop
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Number of RK stages, nRK_stage = ", nRK_stage
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Build a 1D mesh in the radial direction. Always done on CPU.

      CALL GENERATE_GLOBAL_RADIAL_LOG_MESH

     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global)

! ----------------------------------------------------------------------
! Build the powers for the polynomials so that we can zip through them
! automatically during stencil construction. Always done on CPU.

      CALL BUILD_MODE_INDICES

     1         ( n_cc_modes,

     1           ii_mode, jj_mode, kk_mode)

! ----------------------------------------------------------------------
! Build "sector_vertices" for all 6 sectors. Also build a global coordinate
! system on the UNIT sphere and store it in 2D mesh "coordinates_gr_global".
! Always done on CPU.

      coordinates_gr_global = 0.0

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL Mesh_Sector_General ( i, 1,
     1                   1, 1, ntp_global, NGhost_global,
     1sector_vertices ( :, :, i), coordinates_gr_global ( :, :, :, i) )

      END DO

! ----------------------------------------------------------------------
! Previous call to "Mesh_Sector_General" will only fill the interior of
! each sector. The ghost zone information will have to be imported from
! the neighboring sectors. Always done on CPU.

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL FILL_GHOST2D_GLOBAL

     1           ( i, ntp_global, NGhost_global,
     1             coordinates_gr_global)

      END DO

! ----------------------------------------------------------------------
! Extract a patch and save it in "vert_coords_patch". Also make a mask for the
! active zones, "imask_active_patch", set to 1 for active zones. Make sure that
! "isector, ix_patch, iy_patch, iz_patch" and "ilevel, ilevel_max" have
! been set correctly. Always done on CPU.

      CALL EXTRACT_VERT_COORDS_PATCH

     1       ( ntp_global, nr_global, NGhost_global,
     1         ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         isector, ix_patch, iy_patch, iz_patch,
     1         ilevel, ilevel_max,

     1         coordinates_gr_global, rb_save_global,

     1         vert_coords_patch, imask_active_patch)

      WRITE ( 6, *)"-------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Writing out 8 endpoints of this patch."
      WRITE ( 6, *)" "

      ixx = 1  ;  iyy = 1  ;  izz = 1
      WRITE ( 6, *)"ixx, iyy, izz = ", ixx, iyy, izz
      WRITE ( 6, *)"vert_coords_patch = ",
     1              vert_coords_patch ( ixx, iyy, izz, 1),
     1              vert_coords_patch ( ixx, iyy, izz, 2),
     1              vert_coords_patch ( ixx, iyy, izz, 3)
      WRITE ( 6, *)" "

      ixx = ntp_zones  ;  iyy = 1  ;  izz = 1
      WRITE ( 6, *)"ixx, iyy, izz = ", ixx, iyy, izz
      WRITE ( 6, *)"vert_coords_patch = ",
     1              vert_coords_patch ( ixx, iyy, izz, 1),
     1              vert_coords_patch ( ixx, iyy, izz, 2),
     1              vert_coords_patch ( ixx, iyy, izz, 3)
      WRITE ( 6, *)" "

      ixx = 1  ;  iyy = ntp_zones  ;  izz = 1
      WRITE ( 6, *)"ixx, iyy, izz = ", ixx, iyy, izz
      WRITE ( 6, *)"vert_coords_patch = ",
     1              vert_coords_patch ( ixx, iyy, izz, 1),
     1              vert_coords_patch ( ixx, iyy, izz, 2),
     1              vert_coords_patch ( ixx, iyy, izz, 3)
      WRITE ( 6, *)" "

      ixx = ntp_zones  ;  iyy = ntp_zones  ;  izz = 1
      WRITE ( 6, *)"ixx, iyy, izz = ", ixx, iyy, izz
      WRITE ( 6, *)"vert_coords_patch = ",
     1              vert_coords_patch ( ixx, iyy, izz, 1),
     1              vert_coords_patch ( ixx, iyy, izz, 2),
     1              vert_coords_patch ( ixx, iyy, izz, 3)
      WRITE ( 6, *)" "

      ixx = 1  ;  iyy = 1  ;  izz = nr_zones
      WRITE ( 6, *)"ixx, iyy, izz = ", ixx, iyy, izz
      WRITE ( 6, *)"vert_coords_patch = ",
     1              vert_coords_patch ( ixx, iyy, izz, 1),
     1              vert_coords_patch ( ixx, iyy, izz, 2),
     1              vert_coords_patch ( ixx, iyy, izz, 3)
      WRITE ( 6, *)" "

      ixx = ntp_zones  ;  iyy = 1  ;  izz = nr_zones
      WRITE ( 6, *)"ixx, iyy, izz = ", ixx, iyy, izz
      WRITE ( 6, *)"vert_coords_patch = ",
     1              vert_coords_patch ( ixx, iyy, izz, 1),
     1              vert_coords_patch ( ixx, iyy, izz, 2),
     1              vert_coords_patch ( ixx, iyy, izz, 3)
      WRITE ( 6, *)" "

      ixx = 1  ;  iyy = ntp_zones  ;  izz = nr_zones
      WRITE ( 6, *)"ixx, iyy, izz = ", ixx, iyy, izz
      WRITE ( 6, *)"vert_coords_patch = ",
     1              vert_coords_patch ( ixx, iyy, izz, 1),
     1              vert_coords_patch ( ixx, iyy, izz, 2),
     1              vert_coords_patch ( ixx, iyy, izz, 3)
      WRITE ( 6, *)" "

      ixx = ntp_zones  ;  iyy = ntp_zones  ;  izz = nr_zones
      WRITE ( 6, *)"ixx, iyy, izz = ", ixx, iyy, izz
      WRITE ( 6, *)"vert_coords_patch = ",
     1              vert_coords_patch ( ixx, iyy, izz, 1),
     1              vert_coords_patch ( ixx, iyy, izz, 2),
     1              vert_coords_patch ( ixx, iyy, izz, 3)
      WRITE ( 6, *)" "

      WRITE ( 6, *)"-------------------------------------------------"

! ----------------------------------------------------------------------

! We have built the two arrays "vert_coords_patch, imask_active_patch" for
! the specific patch given by "ix_patch, iy_patch, iz_patch". From
! these two arrays, all the other geometric information is built on the GPU.
! The geometry-building, and stencil-building, routine that are exercised below
! have been highly vectorized and parallelized for the GPU.

! In addition, we also copy in the polynomial modes to the GPU.

! ------------------------------
! Copy in geometry variables from which we can build the entire geometry.

!$acc enter data copyin(                                                &
!$acc&           vert_coords_patch, imask_active_patch,                 &

!$acc&    ii_mode, jj_mode, kk_mode,                                    &
!$acc&    ii_mode_Bx, jj_mode_Bx, kk_mode_Bx, mm_mode_Bx,               &
!$acc&    ii_mode_By, jj_mode_By, kk_mode_By, mm_mode_By,               &
!$acc&    ii_mode_Bz, jj_mode_Bz, kk_mode_Bz, mm_mode_Bz )

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Initialize fluid variables on the CPU and then ship them over to the GPU.
! However, before we can do that, we have to do one pass of geometry-building
! on the CPU.

! ------------------------------
! Serial geometry building done here. That is why we have the
! "do_on_gpu = .FALSE." That causes all the subsequent subroutines
! to execute on the CPU.

      do_on_gpu = .FALSE.

! ------------------------------

      CALL STENCIL_INDEXING_TRICORNERS

     1       ( ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         ix_patch, iy_patch,

     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,

     1    imask_active_patch, istencil_patch)

! ------------------------------

      CALL COORDS_TO_WT_QUADPTS_FACES_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length)

! Build the quadrature points, weights and normals for the integration
! of magnetic field in the Multi-D Riemann solver.

      CALL COORDS_TO_WT_QUADPTS_BFLD_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

! ------------------------------

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

! ------------------------------

      CALL INIT_BASIS_FUNCTION_COEFS

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         imask_active_patch,

     1         v_quadpts_patch, v_quadwt_patch, v_volume_patch,
     1         v_centroid_patch, v_scale_length,

     1         Coef_Basis_Moments)

! ------------------------------
! If we get here, we do a one-time initialization of the MHD variables
! on the mesh. The initialization code has OpenMP and OpenACC, but
! is not worth initializing on the GPU. Safest to do this on the CPU.

! Now that we have built geometry on the CPU, we are ready to initialize
! the fluid variables on the CPU.

      a = cfl_coef
      cfl_coef = 0.05

!$acc update device ( cfl_coef)

      CALL INIT_PATCH_CUBE_SPHERE

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           iproblem, time, dtcur,
     1           rmin, rmax,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_quadpts_patch, v_quadwt_patch,

     1           U, V, entropygr, Bxf, Byf, Bzf,
     1           B0field,

     1           gamma_eos_gr,

     1           V_min_local, V_max_local)

      cfl_coef = a  ! Reset the CFL number to its original value.

! ------------------------------

! Since we have built everyting on CPU, we send
! "U, V, entropygr, Bxf, Byf, Bzf" to GPU.

!$acc update device ( U, V, entropygr, Bxf, Byf, Bzf)

! Since we have built everyting on CPU, we send "time, dtcur, cfl_coef" to GPU.

!$acc update device ( time, dtcur, cfl_coef)

! ------------------------------

      WRITE ( 6, *)"-------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Initialized time, dtcur = ", time, dtcur
      WRITE ( 6, *)" "
      WRITE ( 6, *)"-------------------------------------------------"

! ------------------------------
! Now we are ready to copy over variables to the GPU. These are
! "U, entropygr, Bxf, Byf, Bzf". They participate in the MHD update.

!$acc enter data copyin(                                                &
!$acc&    U, entropygr, V, Bxf, Byf, Bzf, B1_zc, B2_zc, B3_zc,          &
!$acc&    gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr)

! ------------------------------
! Now that we have set up the initial conditions and moved the MHD data
! over to the GPU, we put "do_on_gpu = .TRUE.". This allows all the
! Open ACC loops in all ACC-enabled subroutines to execute on the GPU.

! So far, we have been completely serial on CPU. Therefore, do one more
! update of all "(/RIEMANN_COM/)" variables on GPU.

      do_on_gpu = .TRUE.

!$acc update device  (/RIEMANN_COM/)
#if ( DEBUG_SILO == 0 )
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Main timestep loop begins here.

      DO istep = 1, ntstep

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      CALL CPU_TIME ( timing_info1)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Since "time, dtcur" etc. may have changed on CPU, synchronize them with
! GPU.

      do_on_gpu = .TRUE.

!$acc update device ( time, dtcur, cfl_coef)

!$acc update device (/RIEMANN_COM/)

! ----------------------------------------------------------------------
! Build the geometry.
! ----------------------------------------------------------------------

! At the tri-corners of a cube sphere mesh, we cannot index with "ixx + 1"
! etc. This subroutine works out the stencils for the zones that are one or two
! layers from a tri-corner. It maps the stencil back into the physical part of
! the mesh. Done on GPU.

      CALL STENCIL_INDEXING_TRICORNERS

     1       ( ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         ix_patch, iy_patch,

     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,

     1    imask_active_patch, istencil_patch)

! ----------------------------------------------------------------------

! Start with "vert_coords_patch" and build all the area-related quadrature points,
! weights and normals and areas for all the areas that bound each zone. Also build
! a length "v_scale_length" by which each zone can be scaled. Done on GPU.

      CALL COORDS_TO_WT_QUADPTS_FACES_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length)


! Build the quadrature points, weights and normals for the integration
! of magnetic field in the Multi-D Riemann solver.

      CALL COORDS_TO_WT_QUADPTS_BFLD_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

! ----------------------------------------------------------------------

! Start with "vert_coords_patch" and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone. Done on GPU.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

! ----------------------------------------------------------------------

! Now build the moments for that patch. These are useful for building the
! Taylor series basis functions into which the reconstruction problem
! is projected. Done on GPU.

      CALL INIT_BASIS_FUNCTION_COEFS

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         imask_active_patch,

     1         v_quadpts_patch, v_quadwt_patch, v_volume_patch,
     1         v_centroid_patch, v_scale_length,

     1         Coef_Basis_Moments)

! ----------------------------------------------------------------------
! End of building the geometry.
! ----------------------------------------------------------------------

! --------------------------------------------------
! Pad the boundaries in the radial direction if this is the innermost
! radial patch or the outermost radial patch.

! ------------------------------
! Most of the geomesh simulations require to implement the outflow  boundary
! conditions at the outer boundary. The PAD_BOUNDARY subroutine does that.

!      bcarr_r ( :, :, 1)  = 17   ! No inner boundary enforcement here.
      bcarr_r ( :, :, 1)  = 2    ! Enforce outflow conditions at Inner boundary.
      bcarr_r ( :, :, 2)  = 2    ! Enforce outflow conditions at Outer boundary.

!$acc update device ( bcarr_r)

      CALL PAD_BOUNDARY_R

     1    ( nr_zones, ntp_zones, NGhost,
     1      n_cc_components, n_cc_modes,
     1      n_fc_modes,

     1      vert_coords_patch, imask_active_patch,

     1      xf_area_patch, yf_area_patch, zf_area_patch,
     1      v_volume_patch,

     1      v_centroid_patch,

     1      U, entropygr, Bxf, Byf, Bzf,

     1      bcarr_r)

! --------------------

        ibdy_inner = 0   ! 1 == Impose Add_boundary at inner; 0 == No Add_Bdy
        ibdy_outer = 0   ! 1 == Impose Add_boundary at outer; 0 == No Add_Bdy

!$acc update device ( ibdy_inner, ibdy_outer )

      CALL ADD_BOUNDARY

     1         ( ibdy_inner, ibdy_outer,

     1           nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           iproblem, time, dtcur,
     1           rmin, rmax,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_centroid_patch,

     1           U, V, entropygr, Bxf, Byf, Bzf,

     1           gamma_eos_gr, bcarr_r,

     1           V_min_local, V_max_local)

! ----------------------------------------------------------------------
! Main MHD algorithm for a single stage of a Runge-Kutta timestep begins here.
! ----------------------------------------------------------------------

! Call limiter:-

! ------------------------------

      CALL LIMIT_ALL_SLOPES_P1

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           iz_patch, nr_patches,

     1         ii_mode, jj_mode, kk_mode,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch,

     1           v_scale_length,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc, Ex_zc, Ey_zc, Ez_zc, Wzc

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr, workgr,
     1         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,

     1         V_min_local, V_max_local)

! The way it has worked out is that the CPU has already updated
! "V_min_local, V_max_local" on the CPU. So no need to pull them over to CPU.

      WRITE ( 6, *)"-------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"After LIMIT_ALL_SLOPES_P1 "
      WRITE ( 6, *)" "
      WRITE ( 6, *)"V_min_local = ", V_min_local
      WRITE ( 6, *)" "
      WRITE ( 6, *)"V_max_local = ", V_max_local
      WRITE ( 6, *)" "
      WRITE ( 6, *)"-------------------------------------------------"
      CALL FLUSH ( 6)

! ----------------------------------------------------------------------

! Build all 1D Fluxes via calls to 1D Riemann solvers:-

! ------------------------------

      CALL MAKE_FACIAL_FLUXES

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         ibdy_outer,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           v_scale_length,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr,

     1         S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr,

     1         Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1         Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi, dU_dt,

     1         EntFlux_x, EntFlux_y, EntFlux_z,
     1         Efx, Efy, Efz)

! ----------------------------------------------------------------------

! Build all Electric fields via calls to 2D Riemann solvers:-
! The fluxes from the 2D Riemann solvers can also give more Multid aspect
! to the 1D fluid fluxes.

! ------------------------------

      CALL MAKE_EDGE_EFIELDS_P1

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1         ii_mode, jj_mode, kk_mode,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,
     1           xf_area_patch, xf_centroid_patch,

     1           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,
     1           yf_area_patch, yf_centroid_patch,

     1           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,
     1           zf_area_patch, zf_centroid_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch,

     1           v_scale_length,
     1           v_volume_patch, v_centroid_patch,

     1           Coef_Basis_Moments,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc, Ex_zc, Ey_zc, Ez_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr,
     1         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,

     1         S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr,
     1         Efx, Efy, Efz)


! ----------------------------------------------------------------------
! STEP 7) With geometry built, and with all fluxes in-hand, carry out
! the update.

      dtmove = dtcur

!$acc update device ( time, dtcur, dtmove )

      CALL UPDATE_RK_P1

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes, time, dtmove,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadnorml_patch, yf_quadnorml_patch,
     1           zf_quadnorml_patch,

     1           xf_area_patch, yf_area_patch, zf_area_patch,
     1           v_scale_length, v_volume_patch, v_centroid_patch,

     1         U, entropygr, Bxf, Byf, Bzf, Source,

     1         EntFlux_x, EntFlux_y, EntFlux_z,
     1         Efx, Efy, Efz)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Evaluate new timestep "dtcurp1" and move the timesteps. Increment "time"
! and make sure it is consistently available on the CPU and GPU.

      time = time + dtcur

! Send new time to GPU so that any user-supplied subroutine that keys off
! on "time" can use it.

!$acc update device ( time)

      CALL EVALUATE_TIMESTEP

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes, n_cc_modes_one_higher,
     1           n_fc_modes, n_cc_bfld_modes,

     1           time, dtcur, dtcurp1,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           v_centroid_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr)

      WRITE ( 6, *)"-------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"istep, time, dtcur, dtcurp1 = ",
     1              istep, time, dtcur, dtcurp1
      WRITE ( 6, *)" "
      WRITE ( 6, *)"-------------------------------------------------"
      CALL FLUSH ( 6)

      dtcur = dtcurp1  ! Swap in the new timestep. Done on CPU.

! Send new timestep "dtcur" to GPU.

!$acc update device ( dtcur, dtcurp1)

! ----------------------------------------------------------------------
! Simple serial imaging of one data slice done here. The idea is that we
! only bring needed data over to the CPU when imaging or full data dumps
! have to be done.
! ----------------------------------------------------------------------

#if ( GPU_EXECUTION == 0)

      IF ( ( time .GE. timeimage) .OR.
     1     ( time .GE. timestop) .OR.
     1     ( istep .EQ. ntstep) ) THEN

! ----------------------------------------------------------------------
! Image variables to be brought back from GPU to CPU.

!$acc update host ( V)

! ----------------------------------------------------------------------

! Images written out and file names incremented here. This is because the
! previous subroutine has just initialized all the primitive variables.
! For 3d data one may optionally choose to write out select slices of interest.

      timeimage = timeimage + dtimage

      WRITE ( 6, *)"Writing image files = ", rhofil," ", vlxfil," ",
     1            vlyfil," ", vlzfil," ", prsfil," ", bfxfil," ",
     1            bfyfil," ", bfzfil," ", speciefil
      WRITE ( 6, *)" "
      CALL FLUSH ( 6)

! -----

      OPEN ( unit = 99, file = rhofil)
      r4gr ( :, :, :) = V ( :, :, :, 1, 1)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ntp_zones),
     1                     iyy = 1, ntp_zones),
     1                     izz = nr_zones/2, nr_zones/2)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = vlxfil)
      r4gr ( :, :, :) = V ( :, :, :, 2, 1)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ntp_zones),
     1                     iyy = 1, ntp_zones),
     1                     izz = nr_zones/2, nr_zones/2)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = vlyfil)
      r4gr ( :, :, :) = V ( :, :, :, 3, 1)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ntp_zones),
     1                     iyy = 1, ntp_zones),
     1                     izz = nr_zones/2, nr_zones/2)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = vlzfil)
      r4gr ( :, :, :) = V ( :, :, :, 4, 1)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ntp_zones),
     1                     iyy = 1, ntp_zones),
     1                     izz = nr_zones/2, nr_zones/2)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = prsfil)
      r4gr ( :, :, :) = V ( :, :, :, 5, 1)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ntp_zones),
     1                     iyy = 1, ntp_zones),
     1                     izz = nr_zones/2, nr_zones/2)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = bfxfil)
      r4gr ( :, :, :) = V ( :, :, :, 6, 1)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ntp_zones),
     1                     iyy = 1, ntp_zones),
     1                     izz = nr_zones/2, nr_zones/2)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = bfyfil)
      r4gr ( :, :, :) = V ( :, :, :, 7, 1)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ntp_zones),
     1                     iyy = 1, ntp_zones),
     1                     izz = nr_zones/2, nr_zones/2)
      CLOSE ( 99)

! -----

      OPEN ( unit = 99, file = bfzfil)
      r4gr ( :, :, :) = V ( :, :, :, 8, 1)
      WRITE ( 99, *) ( ( ( r4gr ( ixx, iyy, izz), ixx = 1, ntp_zones),
     1                     iyy = 1, ntp_zones),
     1                     izz = nr_zones/2, nr_zones/2)
      CLOSE ( 99)

! ----------

      CALL INCNAME_8 ( rhofil, 1)
      CALL INCNAME_8 ( vlxfil, 1)
      CALL INCNAME_8 ( vlyfil, 1)
      CALL INCNAME_8 ( vlzfil, 1)
      CALL INCNAME_8 ( prsfil, 1)

      CALL INCNAME_8 ( bfxfil, 1)
      CALL INCNAME_8 ( bfyfil, 1)
      CALL INCNAME_8 ( bfzfil, 1)

! ----------------------------------------------------------------------

      END IF  ! End of "IF ( ( time .GE. timeimage) .OR."

#endif
! ----------------------------------------------------------------------
      CALL CPU_TIME ( timing_info2)
! ----------------------------------------------------------------------

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed this step = ", timing_info2 - timing_info1
      WRITE ( 6, *)" "

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      IF ( time .GE. timestop) THEN
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Stopping because time .GE. timestop", time, timestop
      WRITE ( 6, *)" "
      STOP
      END IF
      CALL FLUSH ( 6)

! ----------------------------------------------------------------------
! Main timestep loop ends here.
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO ! End of "istep = 1, ntstep"

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! #if ( DEBUG_SILO == 0)
#else

!$acc update host (                                                     &
!$acc&           vert_coords_patch, imask_active_patch,                 &
!$acc&                                                                  &
!$acc&           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch, &
!$acc&           xf_area_patch, xf_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch, &
!$acc&           yf_area_patch, yf_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch, &
!$acc&           zf_area_patch, zf_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           v_scale_length,                                        &
!$acc&                                                                  &
!$acc&           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,  &
!$acc&           v_volume_patch, v_centroid_patch,                      &
!$acc&                                                                  &
!$acc&           Coef_Basis_Moments, istencil_patch)


! ----------------------------------------------------------------------

      itri_patch ( 1) = itri_patch_LD
      itri_patch ( 2) = itri_patch_RD
      itri_patch ( 3) = itri_patch_RU
      itri_patch ( 4) = itri_patch_LU

      top_or_bottom ( 1: 2) = 1

      U0 ( :, :, :, :) = U  ( :, :, :, :, 1)
      Bxf0 ( :, :, :)  = Bxf ( :, :, :, 1)
      Byf0 ( :, :, :)  = Byf ( :, :, :, 1)
      Bzf0 ( :, :, :)  = Bzf ( :, :, :, 1)

      CALL WRITEBLOCK ( 1,                             ! block index
     1                  nr_zones, ntp_zones, NGhost,
     1                  top_or_bottom,
     1                  itri_patch,
     1                  vert_coords_patch,
     1                  n_cc_components,
     1                  U0, Bxf0, Byf0, Bzf0
     1                )

      U0 ( :, :, :, :) = 0.0
      Bxf0 ( :, :, :)  = 0.0
      Byf0 ( :, :, :)  = 0.0
      Bzf0 ( :, :, :)  = 0.0

      CALL READBLOCK  ( 1,                             ! block index
     1                  nr_zones, ntp_zones, NGhost,
     1                  top_or_bottom,
     1                  itri_patch,
     1                  vert_coords_patch,
     1                  n_cc_components,
     1                  U0, Bxf0, Byf0, Bzf0
     1                )

! ----------------------------------------------------------------------

      err_U ( :) = 0.0

      DO izz = 1,  nr_zones
      DO iyy = 1, ntp_zones
      DO ixx = 1, ntp_zones

        err_U = MAX ( err_U, ABS ( U  ( ixx, iyy, izz, :, 1)
     1                           - U0 ( ixx, iyy, izz, :)))

      END DO  ! End of "DO ixx = 1, ntp_zones"
      END DO  ! End of "DO iyy = 1, ntp_zones"
      END DO  ! End of "DO izz = 1,  nr_zones"

! ----------

      err_Bxf = 0.0
      err_Byf = 0.0
      err_Bzf = 0.0

      DO izz = 0,  nr_zones
      DO iyy = 0, ntp_zones
      DO ixx = 0, ntp_zones

        err_Bxf = MAX ( err_Bxf, ABS ( Bxf  ( ixx, iyy, izz, 1)
     1                               - Bxf0 ( ixx, iyy, izz)))

        err_Byf = MAX ( err_Byf, ABS ( Byf  ( ixx, iyy, izz, 1)
     1                               - Byf0 ( ixx, iyy, izz)))

        err_Bzf = MAX ( err_Bzf, ABS ( Bzf  ( ixx, iyy, izz, 1)
     1                               - Bzf0 ( ixx, iyy, izz)))

#if 0
      WRITE(6, *)"Bxf_write, Bxf_read = ",
     1            Bxf  ( ixx, iyy, izz, 1),
     1            Bxf0 ( ixx, iyy, izz)

      WRITE(6, *)"Byf_write - Byf_read = ",
     1            Byf  ( ixx, iyy, izz, 1) -
     1            Byf0 ( ixx, iyy, izz)

      WRITE(6, *)"Bzf_write - Bzf_read = ",
     1            Bzf  ( ixx, iyy, izz, 1) -
     1            Bzf0 ( ixx, iyy, izz)
#endif
      END DO  ! End of "DO ixx = 0, ntp_zones"
      END DO  ! End of "DO iyy = 0, ntp_zones"
      END DO  ! End of "DO izz = 0,  nr_zones"

! ------------------------------

      WRITE ( 6, *)"err_U   = ", err_U
      WRITE ( 6, *)"err_Bxf = ", err_Bxf
      WRITE ( 6, *)"err_Byf = ", err_Byf
      WRITE ( 6, *)"err_Bzf = ", err_Bzf

! ----------------------------------------------------------------------
#if 0
      CALL WRITEBLOCK ( 1,                             ! block index
     1                  nr_zones, ntp_zones, NGhost,
     1                  top_or_bottom,
     1                  itri_patch,
     1                  vert_coords_patch,
     1                  n_cc_components,
     1                  U, Bxf, Byf, Bzf
     1                )
#endif
! ----------------------------------------------------------------------
! #if ( DEBUG_SILO == 1)
#endif
! ----------------------------------------------------------------------
! We created a lot of temporary arrays that are used for processing a timestep
! on the GPU. Since we are at the end of the timestep, we can get rid of
! these arrays on the GPU.

!$acc end data
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END PROGRAM TEST_SINGLE_PATCH_TIMESTEPPING


