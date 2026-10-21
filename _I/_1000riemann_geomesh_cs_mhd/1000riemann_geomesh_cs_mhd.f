








      PROGRAM RIEMANN_GEOMESH_CS_MHD

! ----------------------------------------------------------------------
! Written by : Dinshaw Balsara and Sethupathy S.

      USE RIEMANN_3D_CS_AMR_SUPPORT

      IMPLICIT NONE

#include <cubesphere/000com/riemann.com>

#include <cubesphere/000com/cool.com>

      INTEGER nimages, ierr, ERRORCODE, my_rank, me, ilevel,
     1        ixx, iyy, izz, ii, jj, npatch_per_proc,
     1        ixxm1, ixxp1, iyym1, iyyp1, izzm1, izzp1, i, j,
     1        np, ipatch, isector, ix_patch, iy_patch, iz_patch,
     1        nr_geom_zones

      INTEGER nRK_stage, irk, iproblem, ndumpstep, nimgstep,
     1        itri_patch_LD, itri_patch_RD, itri_patch_RU,
     1        itri_patch_LU,

     1        istep, ntstep, filenum, ibdy_inner, ibdy_outer,
     1        nr_geom_patches

      REAL rmin, rmax, rint, rratio, time, dtcur, dtcurp1,
     1     a, b, dt_defcor, time_defcor, dtmove,
     1     timeimage, dtimage, timestop,
     1     timing_write1, timing_write2,
     1     timing_timestep1, timing_timestep2,

     1     geom_dr1

!$acc declare create ( time, dtcur, dtcurp1, ibdy_inner, ibdy_outer,    &
!$acc&                 dt_defcor, time_defcor, dtmove )

      REAL sector_vertices ( 4, 3, 6)

! ----------------------------------------------------------------------
! File names used for imaging.

      CHARACTER*13 silofil

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

      REAL timing_info1, timing_info2, alph_cak

! ----------------------------------------------------------------------
! Variables used in testing the MPI messaging.

      INTEGER iproc, nsilo_images

! ----------------------------------------------------------------------

      OPEN ( unit = 6, file = "riemann_geomesh.out")

! ----------------------------------------------------------------------
! Initialize MPI here.

      CALL MPI_INIT (ierr)


! Synchronize on entry.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Fix up the total number of processors "nimages" and my number in that
! space of processors, "me".

      CALL MPI_COMM_SIZE ( MPI_COMM_WORLD, nimages, ierr)

      CALL MPI_COMM_RANK ( MPI_COMM_WORLD, my_rank, ierr)
      me = my_rank + 1





! Create windows for message passing. Establishes handles for MPI-3
! one-sided get calls.
      CALL WINCREATE_MESSAGING ()

! ----------------------------------------------------------------------
! Create arrays that only need to reside on the GPU. These arrays will
! never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

! Also create fluid variable arrays that are needed by the algorithm and only need
! to reside on the GPU. These arrays are not directly used for the update
! and will never be needed again on the CPU. I.e. they are just work arrays that
! we need for processing the timestep on the GPU.

!$acc data create(                                                      &
!$acc&           xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch, &
!$acc&           xf_area_patch, xf_centroid_patch,                      &

!$acc&           yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch, &
!$acc&           yf_area_patch, yf_centroid_patch,                      &

!$acc&           zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch, &
!$acc&           zf_area_patch, zf_centroid_patch,                      &

!$acc&           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch, &
!$acc&           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch, &
!$acc&           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch, &

!$acc&           v_scale_length,                                        &

!$acc&           v_quadpts_patch, v_quadwt_patch,                       &
!$acc&           v_volume_patch, v_centroid_patch,                      &

!$acc&           isector, ipatch, ilevel, ix_patch, iy_patch, iz_patch, &

!$acc&           Coef_Basis_Moments,                                    &
!$acc&           itri_patch_LD, itri_patch_RD, itri_patch_RU,           &
!$acc&           itri_patch_LU, istencil_patch,                         &

!$acc&         U, entropygr, Bxf, Byf, Bzf,                             &
!$acc&         V, B1_zc, B2_zc, B3_zc, Wzc, Ex_zc, Ey_zc, Ez_zc,        &
!$acc&         Source, CGL_fbeta_gr, Un9, Pelectron_gr, Un10,           &
!$acc&         U_skinny, entropygr_skinny,                              &
!$acc&         Bxf_skinny, Byf_skinny, Bzf_skinny,                      &

!$acc&         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,       &
!$acc&         divvelgr, msonicgr, flattengr, workgr,                   &
!$acc&         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,                   &

!$acc&         V_min_local, V_max_local,                                &

!$acc&         S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr,          &
!$acc&         EntFlux_x, EntFlux_y, EntFlux_z,                         &
!$acc&         Efx, Efy, Efz,                                           &
!$acc&         bcarr_r,                                                 &
!$acc&         flag_entprs, flag_Tfloor,                                &

!$acc&         Flux_snd_r_patch, Flux_snd_l_patch, Flux_snd_u_patch,    &
!$acc&         Flux_snd_d_patch, Flux_snd_t_patch, Flux_snd_b_patch,    &

!$acc&         iproblem, rmin, rmax )





! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Initializing prs and rho floors, inner and outer radius and fluid variables.

      ilevel = 1  ! By default, we work on 1 level.

      rhofloor = 1.0e-12   ! Density floor.
      prsfloor = 1.0e-9    ! Pressure floor.
      smallnum = 1.0e-12   ! Effective small number used in code.

      pi = 4.0 * ATAN ( 1.0)

! Initialize "rmin, rmax" to reasonable values.
      rmin = 1.5
      rmax = 2.5

! If the ratio = 1.0 then there is no ratioing of the radial mesh, when
! the geometric mesh is choosen.
      radial_mesh_ratio = 2.0

      base_grid_rmin = rmin
      base_grid_rmax = rmax

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


! ----------
! Assign the unit vector for the inflow velocity at the outer boundary.
! It is used to correctly set the boundary condition at the outer boundary.
! By default it is assigned as zero.

      inflow_uvec ( 1:3) = 0.0

! ----------
! Initial time is always zero.
      time = 0.0

! ----------
! Decide which physical problem you want to initialize with "iproblem".
! iproblem = 0  for testing the MPI messaging or Regular testing
! iproblem = 2  MHD Blast;
! iproblem = 7  Parker wind test;
! iproblem = 9  SOD test;                iproblem = 10 Brio-Wu test
! iproblem = 13 VF slow solar wind;      iproblem = 16 CAK wind
! iproblem = 14 Heliosphere problem;

      iproblem = 14

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

























! These parameters control the checkpoint dumps and imaging.

! Start with initial time set to 0.Restart will refresh this from "riemann.in".
      time = 0.0

! Time at which we want to output the next image.Set this yourself for restart.
      timeimage = 0.0

! Number of timesteps you want to take. Code stops after these many steps.
      ntstep = 5000000

! Write out full restart dumps every "ndumpstep" steps. By default, all the
! variables are saved in 64 bit unformatted files called "dump0001" etc.
      ndumpstep = 500

! Write out images at every "nimgstep" timesteps.(not used)
      nimgstep = 10

! Time at which you want the simulation to stop.
      timestop = 0.20

! Choose the number of silo images that are expected to be written-out.
      nsilo_images = 5



















      ! TODO GOT TO HERE_________________________________







! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Populate iproblem specific values of magnetic field, angular velocity,
! sound speed etc., The subroutine which is called here, is present in
! the applications file, specific to the application.
! This part is, switched OFF in the main version of the code and
! Switched ON when running a specific application; if needed.

#if( CHOOSE_APPL_CAK == 1)

! For CAK winds problem:-

      CALL GET_CAK_PARAMETERS(
     1     iproblem, nsilo_images,
     1     rmin, rmax, gamma, gasconst, molewt, timestop,
     1     omega_passive, isospeed, dipole_B, alph_cak,
     1     Mstar, newtonsc, tempfloor, iflag_npoles, multipole_Bmag)

! Fill in the cooling table. different cooling functions can be read in by making
! small modifications to this subroutine.
      CALL FILL_CAK_COOLTABLE ( )

!$acc update device  (/COOL_COM/)

#elif( CHOOSE_APPL_HELIO == 1)

      CALL GET_HELIO_PARAMETERS (
     1     iproblem, nsilo_images,
     1     rmin, rmax, gamma, gasconst, molewt, timestop,
     1     omega_passive, isospeed, tempfloor,
     1     tau_phys_cgl, f_lowbd_cgl, tau_phys_sw, tau_phys_lism,
     1     neta_f_cgl, kappa_f_cgl, magB_small_cgl, eps_bdry_cgl,
     1     coef_Pelectron)

#endif


























! update "rmin, rmax".
      base_grid_rmin = rmin
      base_grid_rmax = rmax

! Write out images every "dtimage" amount of time. By default, each variable
! is saved in its own file name in 32 bit formatted form, like "rhoa0001" etc.

      dtimage = timestop / REAL ( nsilo_images)

!$acc update device( iproblem, rmin, rmax )



















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

      divgr_amr ( :, :, :, :) = 0.0

      flag_entprs_amr ( :, :, :, :) = 0.0
      flag_Tfloor_amr  ( :, :, :, :) = 0.0

! ----------------------------------------------------------------------
! By this time we have updated a lot of variables on the CPU. Let us
! migrate that data to the GPU. We still want to continue working
! on CPU, so we set "do_on_gpu = .FALSE.".

      do_on_gpu = .FALSE.

!$acc update device  (/RIEMANN_COM/)
























! ----------------------------------------------------------------------
! IMAGING/SILO  ;  IMAGING/SILO  ;  IMAGING/SILO
! ----------------------------------------------------------------------

! Set the variable names used in the silo file.

      silo_allvars ( 1) = "Density "
      silo_allvars ( 2) = "XMomentum "
      silo_allvars ( 3) = "YMomentum "
      silo_allvars ( 4) = "ZMomentum "
      silo_allvars ( 5) = "Energy "
      silo_allvars ( 6) = "XMagneticField "
      silo_allvars ( 7) = "YMagneticField "
      silo_allvars ( 8) = "ZMagneticField "

#if( IPDE_SYSTEM == 1 && NFLUID >= 1)
      silo_allvars ( 8 + 1) = "Tracer_1 "
#endif

#if( IPDE_SYSTEM == 1 && NFLUID >= 2)
      silo_allvars ( 8 + 2) = "Tracer_2 "
#endif

#if( IPDE_SYSTEM == 1 && NFLUID >= 3)
      silo_allvars ( 8 + 3) = "Tracer_3 "
#endif

#if( IPDE_SYSTEM == 1 && NFLUID >= 4)
      silo_allvars ( 8 + 4) = "Tracer_4 "
#endif

#if( IPDE_SYSTEM == 1 && NFLUID >= 5)
      silo_allvars ( 8 + 5) = "Tracer_5 "
#endif

#if( IPDE_SYSTEM == 1 && NFLUID >= 6)
      silo_allvars ( 8 + 6) = "Tracer_6 "
#endif

! ----------

#if( IPDE_SYSTEM == 2)

      silo_allvars ( 9) = "Pdiff_Ppar_Pperp "

#if ( NVAR >= 10)
      silo_allvars ( 10) = "Tracer_1 "
#endif

#if ( NVAR >= 11)
      silo_allvars ( 11) = "Tracer_2 "
#endif

#endif
! End of "IPDE_SYSTEM == 2".

! ----------

#if( IPDE_SYSTEM == 3)

      silo_allvars ( 9) = "Pdiff_Ppar_Pperp "

      silo_allvars ( 10) = "Prs_electron "

#if ( NVAR >= 11)
      silo_allvars ( 11) = "Tracer_1 "
#endif

#if ( NVAR >= 12)
      silo_allvars ( 12) = "Tracer_2 "
#endif

#if ( NVAR >= 13)
      silo_allvars ( 13) = "Tracer_3 "
#endif

#if ( NVAR >= 14)
      silo_allvars ( 14) = "Tracer_4 "
#endif

#if ( NVAR >= 15)
      silo_allvars ( 15) = "Tracer_5 "
#endif

#if ( NVAR >= 16)
      silo_allvars ( 16) = "Tracer_6 "
#endif


#endif
! End of "IPDE_SYSTEM == 3".

! ----------

#if( USE_B0B1_FORMULATION == 1)
      silo_allvars ( NVAR + NFLUID + 1) = "XTotalMagnField "
      silo_allvars ( NVAR + NFLUID + 2) = "YTotalMagnField "
      silo_allvars ( NVAR + NFLUID + 3) = "ZTotalMagnField "
#endif

! ----------------------------------------
! Visualizing the extra internal variables.

! --------------------
#if( USE_B0B1_FORMULATION == 1)
!-----------

#if ( SILO_EXTRA_DATA >= 1)
      silo_allvars ( NVAR + NFLUID + 3 + 1) = "xvar_data1 "
#endif
#if ( SILO_EXTRA_DATA >= 2)
      silo_allvars ( NVAR + NFLUID + 3 + 2) = "xvar_entropy "
#endif
#if ( SILO_EXTRA_DATA >= 3)
      silo_allvars ( NVAR + NFLUID + 3 + 3) = "xflag_entflag "
#endif
#if ( SILO_EXTRA_DATA >= 4)
      silo_allvars ( NVAR + NFLUID + 3 + 4) = "xflag_tempfloor "
#endif

! ----------
#else
! ----------

#if ( SILO_EXTRA_DATA >= 1)
      silo_allvars ( NVAR + NFLUID + 1) = "xvar_data1 "
#endif
#if ( SILO_EXTRA_DATA >= 2)
      silo_allvars ( NVAR + NFLUID + 2) = "xvar_entropy "
#endif
#if ( SILO_EXTRA_DATA >= 3)
      silo_allvars ( NVAR + NFLUID + 3) = "xflag_entprs "
#endif
#if ( SILO_EXTRA_DATA >= 4)
      silo_allvars ( NVAR + NFLUID + 4) = "xflag_tempfloor "
#endif

! ----------
#endif
! --------------------


#if ( SILO_EXTRA_DATA >= 5)
        WRITE( 6, *)" "
        WRITE( 6, *)"Stopping at - Set the silo variable names"
        WRITE( 6, *)"SILO_EXTRA_DATA is expected to be <= 4"
        WRITE( 6, *)"SILO_EXTRA_DATA = ", SILO_EXTRA_DATA
        WRITE( 6, *)" "
        CALL FLUSH( 6)
        ERROR STOP
#endif

! ----------------------------------------
! Required in "WRITE_SILO" routine.

      filenum = 1

! Set silo file name.

      silofil = "pltf0001.silo"

! ------------------------------
! The below filenames are deprecated relative to silo variables.

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

      IF ( me == 1 ) THEN

      WRITE ( 6, *)"--------------------------------------------------"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Starting problem with the following parameters"
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Chose physical problem given by iproblem = ",
     1              iproblem

      WRITE( 6, *)" "
      WRITE( 6, *)"iproblem, nsilo_images = ", iproblem, nsilo_images
      WRITE( 6, *)"rmin, rmax = ", rmin, rmax
      WRITE( 6, *)"timestop = ", timestop
      WRITE( 6, *)" "

      WRITE( 6, *)"omega_passive = ", omega_passive
      WRITE( 6, *)"dipole_B = ", dipole_B
      WRITE( 6, *)"isospeed, Mstar, newtonsc, tempfloor = ",
     1             isospeed, Mstar, newtonsc, tempfloor
      WRITE( 6, *)"iflag_npoles   = ", iflag_npoles
      WRITE( 6, *)"multipole_Bmag = ", multipole_Bmag
      WRITE( 6, *)" "

      CALL FLUSH ( 6)

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
      WRITE ( 6, *)"Number of processors used = ", nimages
      WRITE ( 6, *)" "
      WRITE ( 6, *)"--------------------------------------------------"
      CALL FLUSH ( 6)

      END IF ! End of "IF ( me == 1 ) THEN".






! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! The subroutine sequence ""BUILD_BASE_PATCHES, GENERATE_GLOBAL_RADIAL_LOG_MESH,
! BUILD_MODE_INDICES, Mesh_Sector_General, FILL_GHOST2D_GLOBAL" initializes
! the global mesh.

! These subroutines only depend on "rmin, rmax, smallnum,
! ntp_zones, nr_zones, NGhost, ntp_patches, nr_patches".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Here we build sibling info. for the base level mesh.

      CALL BUILD_BASE_PATCHES ()

      npatch_per_proc = (6 * Max_base_patches_per_sector) / nimages

! ----------------------------------------------------------------------
! Build a 1D mesh in the radial direction. This calls a generic on-processor
! serial code. Therefore, synchronize after exit.

! ------------------------------
#if ( CHOOSE_RADIAL_MESH == 0)
! ------------------------------
! Uniform mesh from "rmin" to "rmax". This can be trusted to work automatically.

      CALL GENERATE_GLOBAL_RADIAL_UNIFORM_MESH
     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global )

! ------------------------------
#elif ( CHOOSE_RADIAL_MESH == 1)
! ------------------------------
! Log scaled mesh from "rmin" to "rmax". This can be trusted to work automatically.

      CALL GENERATE_GLOBAL_RADIAL_LOG_MESH
     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global )

! ------------------------------
#elif ( CHOOSE_RADIAL_MESH == 2)
! ------------------------------
! Don't trust this one to work automatically.

      CALL GENERATE_GLOBAL_RADIAL_STRETCHED_MESH
     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global )

#elif ( CHOOSE_RADIAL_MESH == 3)
! Geometric meshing on the inside, log meshing on the outside.
! Don't trust this one to work automatically.

      radial_mesh_ratio = 1.02497   ! Geometric ratioing.
      geom_dr1          = 0.003     ! Size of first zone
      nr_geom_zones     = 90        ! number of radial patches with geometric mesh

      IF( ( nr_global .NE. 240) .OR. ( NINT( rmax) .NE. 15)) THEN
        WRITE( 6, *)" "
        WRITE( 6, *)"Stopping before GENERATE_GLOBAL_RADIAL_HYBRID_MESH"
        WRITE( 6, *)"Mismatch in nr_zones and nr_global or rmax"
        WRITE( 6, *)" "
        CALL FLUSH( 6)
        ERROR STOP
      END IF

      CALL GENERATE_GLOBAL_RADIAL_HYBRID_MESH

     1       ( nr_zones, nr_global, NGhost_global, nr_geom_zones,
     1         rmin, rmax, smallnum, radial_mesh_ratio, geom_dr1,
     1         rb_save_global, dr_save_global, rc_save_global)

! ------------------------------
#elif ( CHOOSE_RADIAL_MESH == 4)
! ------------------------------
! Set up a power law mesh. ! Don't trust this one to work automatically.

      CALL GENERATE_GLOBAL_RADIAL_POWERLAW_MESH

     1       ( nr_global, NGhost_global, rmin, rmax, smallnum,
     1         radial_mesh_ratio,
     1         rb_save_global, dr_save_global, rc_save_global)

! ------------------------------
#elif ( CHOOSE_RADIAL_MESH == 5)
! ------------------------------
! Set up mesh with uniform zones up to "rint" and log scaled zones on the outside.
! Will reset "rmax". Don't trust this one to work automatically.

        rint = 250  ! Intermediate radius; within it we have uniform zones.
        nr_geom_zones = 418  ! Total zones 688 = 270 inner + 418 outer

      CALL GENERATE_GLOBAL_RADIAL_UNIF_N_LOG_MESH
     1       ( nr_global, NGhost_global, nr_geom_zones,
     1         rmin, rmax, rint, smallnum,
     1         rb_save_global, dr_save_global, rc_save_global)

! ------------------------------
#elif ( CHOOSE_RADIAL_MESH >= 6 && CHOOSE_RADIAL_MESH <= 9)
! ------------------------------

! Build a 1D mesh in the radial direction. This calls a generic on-processor
! serial code. Therefore, synchronize after exit. This is Vladimir's favorite mesh.

! If you multiply this with "rmin" you get location of transition.
      radial_mesh_transition = 10.0

! Code tries to keep this as the ratio between largest to smallest radial zone.
      radial_mesh_ratio = 10.0

! Fraction of all points that lie in the uniform interior part.
      radial_mesh_fraction = 0.4

      CALL GENRADIALMESH
     1       (CHOOSE_RADIAL_MESH, nr_global, NGhost_global, rmin, rmax,
     1        rmin * radial_mesh_transition, radial_mesh_ratio,
     1        radial_mesh_fraction, rb_save_global, dr_save_global,
     1        rc_save_global)

! ------------------------------
#endif
! ------------------------------

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Build the powers for the polynomials so that we can zip through them
! automatically during stencil construction.

      CALL BUILD_MODE_INDICES

     1         ( n_cc_modes,

     1           ii_mode, jj_mode, kk_mode)

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ------------------------------
! Copy in gamma, molewt grid variables and then the modes of U and magnetic field.

!$acc enter data copyin(                                                &
!$acc&    vert_coords_patch, imask_active_patch,                        &

!$acc&    gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,            &

!$acc&    ii_mode, jj_mode, kk_mode)

! ----------------------------------------------------------------------

! ----------------------------------------------------------------------

! Build "sector_vertices" for all 6 sectors. Also build a global coordinate
! system on the UNIT sphere and store it in 2D mesh "coordinates_gr_global".

      coordinates_gr_global = 0.0

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL Mesh_Sector_General ( i, 1,
     1                   1, 1, ntp_global, NGhost_global,
     1sector_vertices ( :, :, i), coordinates_gr_global ( :, :, :, i) )

      END DO

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! Previous call to "Mesh_Sector_General" will only fill the interior of
! each sector. The ghost zone information will have to be imported from
! the neighboring sectors.

      DO i = 1, 6  ! Loop over the 6 sectors.

      CALL FILL_GHOST2D_GLOBAL

     1           ( i, ntp_global, NGhost_global,
     1             coordinates_gr_global)

      END DO

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Loop over patches, building their geometry and initializing MHD
! variables for all patches on this processor.
! ----------------------------------------------------------------------
! Since we are dealing with only one level mesh, each processor will work on
! "npatch_per_proc" number of patches. This variable has already been built.
! We loop over patches here.

! Note that "iactive_FromSector, isibling_number_t, isibling_sector_t, isibling_proc_t"
! and "iactive_FromProc" to be properly initialized on each processor. This
! was already done via calll to "BUILD_BASE_PATCHES".
! ----------------------------------------------------------------------








! START INITIALIZING GEOMETRY AND PROBLEM  ;  START INITIALIZING GEOMETRY AND PROBLEM









! The following subroutine sequence extracts patch coordinates from the
! global mesh, "EXTRACT_VERT_COORDS_PATCH, STENCIL_INDEXING_TRICORNERS".
! This extraction only needs to be done once and fills in
! "vert_coords_patch, imask_active_patch". Saved in "vert_coords_amr, imask_active_amr".

! From unpacking those two arrays, we can rebuild the entire geometry with the
! sequence, "COORDS_TO_WT_QUADPTS_FACES_ISO2, COORDS_TO_WT_QUADPTS_VOLUME_ISO2,
! INIT_BASIS_FUNCTION_COEFS".

! MHD Data is initialized after geometry with "INIT_PATCH_CUBE_SPHERE".

! ----------------------------------------------------------------------

      dtcurp1_reduce = 1.0e32  ! Set up to globally minimize the first timestep.

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      isector = iactive_FromProc ( ii, 1)
      ipatch = iactive_FromProc ( ii, 2)

! From "ipatch" (which is the logical patch number in "isector") we can
! extract "ix_patch, iy_patch, iz_patch".

      jj = ipatch
      iz_patch = 1 + (jj - 1) / (ntp_patches**2)

      jj = jj - (iz_patch - 1) * ntp_patches**2
      iy_patch = 1 + (jj - 1) / ntp_patches

      jj = jj - (iy_patch - 1) * ntp_patches
      ix_patch = jj

! ------------------------------

! Fill in "isector" and level information. Also fill in "ix_patch, iy_patch, iz_patch"
! information.
      isector_ilevel_amr ( 1, ii) = isector
      isector_ilevel_amr ( 2, ii) = 1  ! For now, only one level in the AMR hierarchy.

! Fill in "ix_patch, iy_patch, iz_patch" for this patch
      ixyz_patch_amr ( 1, ii) = ix_patch
      ixyz_patch_amr ( 2, ii) = iy_patch
      ixyz_patch_amr ( 3, ii) = iz_patch


! ----------------------------------------------------------------------
! Extract a patch and save it in "vert_coords_patch". Also make a mask for the
! active zones, "imask_active_patch", set to 1 for active zones. Make sure that
! "isector, ix_patch, iy_patch, iz_patch" and "ilevel, ilevel_max" have
! been set correctly.

      CALL EXTRACT_VERT_COORDS_PATCH

     1       ( ntp_global, nr_global, NGhost_global,
     1         ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         isector, ix_patch, iy_patch, iz_patch,
     1         ilevel, ilevel_max,

     1         coordinates_gr_global, rb_save_global,

     1         vert_coords_patch, imask_active_patch)

! ------------------------------
! Save "vert_coords_patch, imask_active_patch" for later reuse. These are
! also needed for imaging.

      vert_coords_amr ( :, :, :, :, ii) =
     1                              vert_coords_patch ( :, :, :, :)

      imask_active_amr ( :, :, :, ii) = imask_active_patch ( :, :, :)

! ------------------------------

! Understand that we are working on a single patch on each processor.
! We have just filled in "vert_coords_patch, imask_active_patch" which is
! the minimal amount of information needed for regenating geometry,
! quadrature points, weights etc.

! ----------------------------------------------------------------------
! At the tri-corners of a cube sphere mesh, we cannot index with "ixx + 1"
! etc. This subroutine works out the stencils for the zones that are one or two
! layers from a tri-corner. It maps the stencil back into the physical part of
! the mesh.

      CALL STENCIL_INDEXING_TRICORNERS

     1       ( ntp_zones, nr_zones, ntp_patches, nr_patches, NGhost,

     1         ix_patch, iy_patch,

     1    itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,

     1    imask_active_patch, istencil_patch)

! ------------------------------
! Save information on whether corners are tri-corners or quad-corners.
! Also, save stencil information because for cube-sphere this is only 2D
! and very lightweight.

      itri_patch_amr ( 1, ii) = itri_patch_LD
      itri_patch_amr ( 2, ii) = itri_patch_RD
      itri_patch_amr ( 3, ii) = itri_patch_RU
      itri_patch_amr ( 4, ii) = itri_patch_LU

      istencil_amr ( :, :, :, :, ii) = istencil_patch ( :, :, :, :)

! ----------------------------------------------------------------------
! All geometry data past this point can be built using "sector #, level #",
! "ix_patch, iy_patch, iz_patch", "itri_patch_LD, itri_patch_RD, itri_patch_RU,
! itri_patch_LU" and "vert_coords_patch, imask_active_patch, istencil_patch".













! ----------------------------------------------------------------------
! BEGIN BUILDING GEOMETRY  ;  BEGIN BUILDING GEOMETRY  ;  BEGIN BUILDING GEOMETRY
! ----------------------------------------------------------------------


















! Start with "vert_coords_patch, imask_active_patch"
! and build all the area-related quadrature points,
! weights and normals and areas for all the areas that bound each zone. Also build
! a length "v_scale_length" by which each zone can be scaled.

! This gives us interior weights and normals that will be used for flux evaluation.

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

! Build the quadrature points, weights and normals for the 1D edge integration
! of magnetic field in the Multi-D Riemann solver.

! This gives us weights and normals at the edges that will be used
! for electric field evaluation.

      CALL COORDS_TO_WT_QUADPTS_BFLD_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,
     1           yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,
     1           zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch)

! ----------------------------------------------------------------------

! Start with "vert_coords_patch, imask_active_patch"
! and build all the volume-related quadrature points,
! weights and tangents to the reference coordinates. Also build the volume
! of each zone and the centroid of each zone.

      CALL COORDS_TO_WT_QUADPTS_VOLUME_ISO2

     1         ( nr_zones, ntp_zones, NGhost,

     1           vert_coords_patch, imask_active_patch,

     1           v_quadpts_patch, v_quadwt_patch, v_quadtangent_patch,
     1           v_volume_patch, v_centroid_patch)

! ----------------------------------------------------------------------

! Now build the moments for that patch. These are useful for building the
! Taylor series basis functions into which the reconstruction problem is projected.

      CALL INIT_BASIS_FUNCTION_COEFS

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_modes,

     1         ii_mode, jj_mode, kk_mode,

     1         imask_active_patch,

     1         v_quadpts_patch, v_quadwt_patch, v_volume_patch,
     1         v_centroid_patch, v_scale_length,

     1         Coef_Basis_Moments)

! ----------------------------------------------------------------------
! Now we save the geometry information from "_patch" variables
! to "_amr" variables.

      xf_quadpts_amr ( :,:,:, :, :, ii)= xf_quadpts_patch ( :,:,:, :, :)
      yf_quadpts_amr ( :,:,:, :, :, ii)= yf_quadpts_patch ( :,:,:, :, :)
      zf_quadpts_amr ( :,:,:, :, :, ii)= zf_quadpts_patch ( :,:,:, :, :)

      xf_quadwt_amr ( :,:,:, :, ii) = xf_quadwt_patch ( :,:,:, :)
      yf_quadwt_amr ( :,:,:, :, ii) = yf_quadwt_patch ( :,:,:, :)
      zf_quadwt_amr ( :,:,:, :, ii) = zf_quadwt_patch ( :,:,:, :)

      xf_quadnorml_amr ( :,:,:, :, ii) = xf_quadnorml_patch ( :,:,:, :)
      yf_quadnorml_amr ( :,:,:, :, ii) = yf_quadnorml_patch ( :,:,:, :)
      zf_quadnorml_amr ( :,:,:, :, :, ii) =
     1                             zf_quadnorml_patch ( :,:,:, :, :)

      xf_area_amr ( :, :, :, ii) = xf_area_patch ( :, :, :)
      yf_area_amr ( :, :, :, ii) = yf_area_patch ( :, :, :)
      zf_area_amr ( :, :, :, ii) = zf_area_patch ( :, :, :)

      xf_centroid_amr ( :,:,:, :, ii) = xf_centroid_patch ( :,:,:, :)
      yf_centroid_amr ( :,:,:, :, ii) = yf_centroid_patch ( :,:,:, :)
      zf_centroid_amr ( :,:,:, :, ii) = zf_centroid_patch ( :,:,:, :)

! --------------------

      xb_quadpts_amr ( :,:,:, :, :, ii)= xb_quadpts_patch ( :,:,:, :, :)
      yb_quadpts_amr ( :,:,:, :, :, ii)= yb_quadpts_patch ( :,:,:, :, :)
      zb_quadpts_amr ( :,:,:, :, :, ii)= zb_quadpts_patch ( :,:,:, :, :)

      xb_quadwt_amr ( :,:,:, :, ii) = xb_quadwt_patch ( :,:,:, :)
      yb_quadwt_amr ( :,:,:, :, ii) = yb_quadwt_patch ( :,:,:, :)
      zb_quadwt_amr ( :,:,:, :, ii) = zb_quadwt_patch ( :,:,:, :)

      xb_quadnorml_amr ( :,:,:, :, ii) = xb_quadnorml_patch ( :,:,:, :)
      yb_quadnorml_amr ( :,:,:, :, ii) = yb_quadnorml_patch ( :,:,:, :)
      zb_quadnorml_amr ( :,:,:, :, :, ii) =
     1                             zb_quadnorml_patch ( :,:,:, :, :)

! --------------------

      v_scale_length_amr ( :, :, :, ii) = v_scale_length ( :, :, :)

      v_quadpts_amr ( :,:,:, :, :, ii) = v_quadpts_patch ( :,:,:, :, :)

      v_quadwt_amr ( :, :, :, :, ii) = v_quadwt_patch ( :, :, :, :)

      v_volume_amr ( :, :, :, ii) = v_volume_patch ( :, :, :)

      v_centroid_amr ( :, :, :, :, ii) = v_centroid_patch ( :, :, :, :)

! --------------------

      Coef_Basis_Moments_amr ( :, :, :, :, ii) =
     1                                  Coef_Basis_Moments ( :, :, :, :)

! ----------------------------------------------------------------------
! END BUILDING GEOMETRY  ;  END BUILDING GEOMETRY  ;  END BUILDING GEOMETRY
! ----------------------------------------------------------------------












! ----------------------------------------------------------------------
! START INITIALIZE ; START INITIALIZE ; START INITIALIZE ; START INITIALIZE
! ----------------------------------------------------------------------











! By this point, all geometry is built. Therefore, we are now in a position
! to initialize the MHD variables for this patch.

! ------------------------------
! The restart condition is removed, so as to fill the ghost zone information
! of the "B0" field.

!#if ( RESTART == 0)
! ------------------------------

      a = cfl_coef
      cfl_coef = AMIN1 ( 0.0001, a)

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

! Store local "dtcur" for each patch that resides on this processor.
      dtcurp1_amr ( ii) = dtcur

! For all the patches that reside on this processor, make a reduction of their
! "dtcur" value
      dtcurp1_reduce = AMIN1 ( dtcurp1_reduce, dtcur)

! ----------------------------------------------------------------------
! By this point, all "U, V, entropygr, Bxf, Byf, Bzf, gamma_eos_gr" for this
! "ipatch" are initialized in the interior of the mesh. But each processor
! can have multiple patches, which are all stored in
! "U_amr, entropygr_amr, Bxf_amr, Byf_amr, Bzf_amr".

      U_amr ( :, :, :, 1:n_cc_components, ii) =
     1      U ( :, :, :, 1:n_cc_components, 1)

      entropygr_amr ( :, :, :, ii) = entropygr ( :, :, :, 1)

      Bxf_amr ( :, :, :, ii) = Bxf ( :, :, :, 1)
      Byf_amr ( :, :, :, ii) = Byf ( :, :, :, 1)
      Bzf_amr ( :, :, :, ii) = Bzf ( :, :, :, 1)

      B0field_amr ( :, :, :, 1:3, ii) = B0field ( :, :, :, 1:3)

! ------------------------------
!#endif
! ------------------------------

! ----------------------------------------------------------------------
! END INITIALIZE ; END INITIALIZE ; END INITIALIZE ; END INITIALIZE ;
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Let all processors synchronize before going forward.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ------------------------------
#if ( RESTART == 0)
! ------------------------------

! Now make a reduction across all processors for "dtcurp1_reduce".

      CALL MPI_ALLREDUCE ( dtcurp1_reduce,   ! Input for reudction
     1                     dtcur,            ! Global minimum
     1                     1,                ! length of reduction variable
     1                     MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

      IF ( me .EQ. 1) THEN

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Initialized all patches, dtcur = ", dtcur
      WRITE ( 6, *)" "
      CALL FLUSH ( 6)

      END IF

! ------------------------------
#endif
! ------------------------------

! Let all processors synchronize before going forward.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! END INITIALIZING GEOMETRY AND PROBLEM  ;  END INITIALIZING GEOMETRY AND PROBLEM
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! End of looping over patches, building their geometry and initializing MHD
! variables for all patches on this processor.
! ----------------------------------------------------------------------
! Build, "is_bottom_or_top" information for each patch.
! Used in reading and writing of silo files.

      is_bottom_or_top ( :, :) = 0.0

      DO ii = 1, npatch_per_proc

        iz_patch = ixyz_patch_amr ( 3, ii)

        IF ( iz_patch == 1)           is_bottom_or_top ( 1, ii) = 1
        IF ( iz_patch == nr_patches)  is_bottom_or_top ( 2, ii) = 1

      END DO








! ----------------------------------------------------------------------
! START READING IN RESTART DUMP ; START READING IN RESTART DUMP ;
! ----------------------------------------------------------------------






! Read in a prior dump file if one has been written out. Increment dump
! file name.

! ------------------------------
#if ( RESTART == 1)
! ------------------------------

      IF (me == 1) THEN

      OPEN ( unit = 5, file = "riemann.in")

      READ ( 5, *) time, dtcur, jj
      READ ( 5, *) filenum
      CLOSE ( 5)

      a = time
      b = dtcur

      END IF

! Now broadcast "time, dtcur, filenum" to all processors.
      CALL MPI_Bcast ( time,    1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast ( dtcur,   1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast ( filenum, 1, MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

! Let all processors synchronize before going forward.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------
! Now read the restart dump.

      IF ( me == 1) THEN
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Reading dumpfile with filenum = ", filenum
      WRITE ( 6, *)" "
      END IF

! --------------------------------------------------

#if( USE_B0B1_FORMULATION == 1 )

      CALL READDATA ( npatch_per_proc,
     1                iactive_FromProc ( :, 2),
     1                nr_zones, ntp_zones, NGhost,
     1                is_bottom_or_top,
     1                itri_patch_amr,
     1                vert_coords_amr,
     1                n_cc_components + 3 + SILO_EXTRA_DATA,
     1                SILO_VARLEN, silo_allvars,
     1                U_silo_amr, Bxf_amr, Byf_amr, Bzf_amr,
     1                filenum, time
     1              )

      DO ii = 1, npatch_per_proc

            U_amr ( :, :, :, 1:n_cc_components, ii) =
     1 U_silo_amr ( :, :, :, 1:n_cc_components, ii)

#if ( SILO_EXTRA_DATA >= 2)
      entropygr_amr ( :, :, :, ii) =
     1   U_silo_amr ( :, :, :, n_cc_components + 3 + 2, ii)
#endif

      END DO

#else

      CALL READDATA ( npatch_per_proc,
     1                iactive_FromProc ( :, 2),
     1                nr_zones, ntp_zones, NGhost,
     1                is_bottom_or_top,
     1                itri_patch_amr,
     1                vert_coords_amr,
     1                n_cc_components + SILO_EXTRA_DATA,
     1                SILO_VARLEN, silo_allvars,
     1                U_silo_amr, Bxf_amr, Byf_amr, Bzf_amr,
     1                filenum, time
     1              )

      DO ii = 1, npatch_per_proc

            U_amr ( :, :, :, 1:n_cc_components, ii) =
     1 U_silo_amr ( :, :, :, 1:n_cc_components, ii)

#if ( SILO_EXTRA_DATA >= 2)
      entropygr_amr ( :, :, :, ii) =
     1   U_silo_amr ( :, :, :, n_cc_components + 2, ii)
#endif

      END DO

#endif

! --------------------------------------------------

      IF ( me == 1) THEN

        WRITE ( 6, *)"From riemann.in; time, dtcur, n_cc_components = ",
     1                a, b, jj
        WRITE ( 6, *)"From dump; time, dtcur, n_cc_components = ",
     1                time, dtcur, n_cc_components

        IF ( ABS ( time - a) >= smallnum) THEN
        WRITE ( 6, *)
     1  "time does not match between the restart dump and riemann.in",
     1  time, a, ABS ( time - a), smallnum
        CALL FLUSH ( 6)
        STOP
        END IF

        IF ( ABS ( dtcur - b) >= smallnum) THEN
        WRITE ( 6, *)
     1  "dtcur does not match between the restart dump and riemann.in",
     1  dtcur, b, ABS ( dtcur - b), smallnum
        CALL FLUSH ( 6)
        STOP
        END IF

        IF ( jj .NE. n_cc_components) THEN
        WRITE ( 6, *)"Stopping because : jj .NE. n_cc_components ",
     1                jj, n_cc_components
        CALL FLUSH ( 6)
        STOP
        END IF

      END IF

      IF ( me == 1) THEN
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Reading of restart file successful !!!"
      WRITE ( 6, *)"Entering time step loop with time, dtcur:",
     1              time, dtcur
      WRITE ( 6, *)" "
      END IF

      timeimage = time

      filenum = filenum + 1

! Let all processors synchronize before going forward.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ------------------------------
#endif
! ----------------------------------------------------------------------
! END READING IN RESTART DUMP ; END READING IN RESTART DUMP ;
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

















! Now that we have set up the initial conditions we put "do_on_gpu = .TRUE.".
! This allows all the Open ACC loops in all ACC-enabled subroutines to
! execute on the GPU.

! So far, we have been completely serial on CPU. Therefore, do one more
! update of all "(/RIEMANN_COM/)" variables on GPU.

      do_on_gpu = .TRUE.

!$acc update device (/RIEMANN_COM/)
!$acc update device ( time, dtcur, cfl_coef)

! ----------------------------------------------------------------------














! ----------------------------------------------------------------------
! START TIMESTEP LOOP  ;  START TIMESTEP LOOP  ;  START TIMESTEP LOOP
! ----------------------------------------------------------------------











      DO istep = 1, ntstep

      CALL CPU_TIME ( timing_timestep1)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

! Loop over Runge-Kutta stages.

      DO irk = 1, nRK_stage

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Synchronize all ghost zone data for "U_amr, entropygr_amr, Bxf_amr,Byf_amr,Bzf_amr".
! All these subroutines have "CALL MPI_Barrier" at begin and end.
! All these MPI calls are on CPU.

      CALL LOAD_LOCAL_BUFFERS_UXY ()
      CALL FETCH_BUFFERS_FILL_IN_UX ()

      CALL LOAD_LOCAL_BUFFERS_UXY ()
      CALL FETCH_BUFFERS_FILL_IN_UY ()

      CALL LOAD_LOCAL_BUFFERS_UXY ()
      CALL FETCH_BUFFERS_FILL_IN_UY_TOP ()

      CALL LOAD_LOCAL_BUFFERS_UZ ()
      CALL FETCH_BUFFERS_FILL_IN_UZ ()

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"LOAD_LOCAL_BUFFERS_UXY, istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif






! ----------------------------------------------------------------------
! START WRITING A DUMP FILE; 6 FILES FOR EACH OF THE SECTORS.
! ----------------------------------------------------------------------






! If we get here, we write out a dump. The ghost zones in the theta-phi
! directions are filled in. Not so in the r-direction.

      IF ( ( time .GE. timeimage) .OR. ( time .GE. timestop)) THEN

! ----------------------------------------------------------------------

      IF ( me == 1) THEN

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Writing image file with filenum = ", filenum
      WRITE ( 6, *)" "
      WRITE ( 6, *)"time, dtcur = ", time, dtcur
      WRITE ( 6, *)" "

      CALL FLUSH ( 6)


      OPEN ( unit = 5, file = "riemann.in")

      WRITE ( 5, *) time, dtcur, n_cc_components
      WRITE ( 5, *) filenum

      CALL FLUSH ( 5)
      CLOSE ( 5)

      END IF

! --------------------------------------------------

      CALL CPU_TIME ( timing_write1)

! ----------------------------------------
#if( USE_B0B1_FORMULATION == 1 )
! ----------------------------------------

      DO ii = 1, npatch_per_proc

        U_silo_amr ( :, :, :, 1:n_cc_components, ii) =
     1       U_amr ( :, :, :, 1:n_cc_components, ii)

        U_silo_amr ( :, :, :, n_cc_components+1:
     1                        n_cc_components+3, ii)=
     1 B0field_amr ( :, :, :, 1:3, ii) + U_silo_amr ( :, :, :, 6:8, ii)


#if ( SILO_EXTRA_DATA >= 1)
        U_silo_amr ( :, :, :, n_cc_components + 3 + 1, ii) =
     1                        divgr_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 2)
        U_silo_amr ( :, :, :, n_cc_components + 3 + 2, ii) =
     1                         entropygr_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 3)
        U_silo_amr ( :, :, :, n_cc_components + 3 + 3, ii) =
     1                     flag_entprs_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 4)
        U_silo_amr ( :, :, :, n_cc_components + 3 + 4, ii) =
     1                      flag_Tfloor_amr ( :, :, :, ii)
#endif

      END DO

      CALL WRITEDATA ( npatch_per_proc,
     1                 iactive_FromProc ( :, 2),
     1                 nr_zones, ntp_zones, NGhost,
     1                 is_bottom_or_top,
     1                 itri_patch_amr,
     1                 vert_coords_amr,
     1                 n_cc_components + 3 + SILO_EXTRA_DATA,
     1                 SILO_VARLEN, silo_allvars,
     1                 U_silo_amr, Bxf_amr, Byf_amr, Bzf_amr,
     1                 filenum, time
     1                )

! ----------------------------------------
#else
! ----------------------------------------

      DO ii = 1, npatch_per_proc

! ------------------------------
! This is only needed if we need properly floored data to be
! sent back for imaging.
#if ( CHOOSE_APPL_HELIO == 1)

        iz_patch = ixyz_patch_amr ( 3, ii)

        U ( :, :, :, 1:n_cc_components, 1) =
     1                 U_amr ( :, :, :, 1:n_cc_components, ii)

        entropygr ( :, :, :, 1) = entropygr_amr ( :, :, :, ii)


        CALL FLOOR_HELIO_PRS_1_LAYER

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,

     1           iz_patch,

     1         U, V, entropygr)


       U_amr ( :, :, :, 1:n_cc_components, ii) =
     1     U ( :, :, :, 1:n_cc_components, 1)

#endif
! ------------------------------

        U_silo_amr ( :, :, :, 1:n_cc_components, ii) =
     1       U_amr ( :, :, :, 1:n_cc_components, ii)

#if ( SILO_EXTRA_DATA >= 1)
        U_silo_amr ( :, :, :, n_cc_components + 1, ii) =
     1                        divgr_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 2)
        U_silo_amr ( :, :, :, n_cc_components + 2, ii) =
     1                         entropygr_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 3)
        U_silo_amr ( :, :, :, n_cc_components + 3, ii) =
     1                     flag_entprs_amr ( :, :, :, ii)
#endif
#if ( SILO_EXTRA_DATA >= 4)
        U_silo_amr ( :, :, :, n_cc_components + 4, ii) =
     1                      flag_Tfloor_amr ( :, :, :, ii)
#endif

! ------------------------------

      END DO

! ------------------------------

      CALL WRITEDATA ( npatch_per_proc,
     1                 iactive_FromProc ( :, 2),
     1                 nr_zones, ntp_zones, NGhost,
     1                 is_bottom_or_top,
     1                 itri_patch_amr,
     1                 vert_coords_amr,
     1                 n_cc_components + SILO_EXTRA_DATA,
     1                 SILO_VARLEN, silo_allvars,
     1                 U_silo_amr, Bxf_amr, Byf_amr, Bzf_amr,
     1                 filenum, time
     1                )

! ----------------------------------------
#endif

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"WRITEDATA, istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! ----------------------------------------

      CALL CPU_TIME ( timing_write2)

! --------------------------------------------------

      IF ( me == 1) THEN
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed writing image:", timing_write2 - timing_write1
      WRITE ( 6, *)" "
      CALL FLUSH ( 6)
      END IF

! ----------------------------------------------------------------------

      filenum = filenum + 1

      timeimage = timeimage + dtimage

! ----------------------------------------------------------------------

      END IF  ! End of "IF ( ( time .GE. timeimage) .OR. ( time .GE. timestop))"

! ----------------------------------------------------------------------
! This is the principal stopping point in this code.

      IF ( time .GE. timestop)  GO TO 12345

! ----------------------------------------------------------------------
! END OF WRITING A DUMP FILE; 6 FILES FOR EACH OF THE SECTORS.














! ----------------------------------------------------------------------
! BEGIN LOOP OVER PATCHES ON PROCESSOR   ;   BEGIN LOOP OVER PATCHES ON PROCESSOR
! ----------------------------------------------------------------------











! Now we loop over all "npatch_per_proc" patches on each processor and
! build fluxes and electric fields.

      dtcurp1_reduce = 1.0e32

      V_min_reduce ( :) = 1.0e32
      V_max_reduce ( :) = -1.0e32

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------










! ----------------------------------------------------------------------
! -STEP- 1) First reverse all that we wrote into the big "*_amr" arrays during
! initialization.












      isector = iactive_FromProc ( ii, 1)
      ipatch = iactive_FromProc ( ii, 2)

      ilevel = isector_ilevel_amr ( 2, ii)

      ix_patch = ixyz_patch_amr ( 1, ii)
      iy_patch = ixyz_patch_amr ( 2, ii)
      iz_patch = ixyz_patch_amr ( 3, ii)

! --------------------------------------------------

      vert_coords_patch ( :, :, :, :) =
     1                          vert_coords_amr ( :, :, :, :, ii)

      imask_active_patch ( :, :, :) =
     1                          imask_active_amr ( :, :, :, ii)

! --------------------------------------------------

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

      istencil_patch ( :, :, :, :) = istencil_amr ( :, :, :, :, ii)

! --------------------------------------------------
! Now we obtain the geometry information from "_amr" variables
! to "_patch" variables.

      xf_quadpts_patch ( :,:,:, :, :)= xf_quadpts_amr ( :,:,:, :, :, ii)
      yf_quadpts_patch ( :,:,:, :, :)= yf_quadpts_amr ( :,:,:, :, :, ii)
      zf_quadpts_patch ( :,:,:, :, :)= zf_quadpts_amr ( :,:,:, :, :, ii)

      xf_quadwt_patch ( :,:,:, :) = xf_quadwt_amr ( :,:,:, :, ii)
      yf_quadwt_patch ( :,:,:, :) = yf_quadwt_amr ( :,:,:, :, ii)
      zf_quadwt_patch ( :,:,:, :) = zf_quadwt_amr ( :,:,:, :, ii)

      xf_quadnorml_patch ( :,:,:, :) = xf_quadnorml_amr ( :,:,:, :, ii)
      yf_quadnorml_patch ( :,:,:, :) = yf_quadnorml_amr ( :,:,:, :, ii)
      zf_quadnorml_patch ( :,:,:, :, :) =
     1                              zf_quadnorml_amr ( :,:,:, :, :, ii)

      xf_area_patch ( :, :, :) = xf_area_amr ( :, :, :, ii)
      yf_area_patch ( :, :, :) = yf_area_amr ( :, :, :, ii)
      zf_area_patch ( :, :, :) = zf_area_amr ( :, :, :, ii)

      xf_centroid_patch ( :,:,:, :) = xf_centroid_amr ( :,:,:, :, ii)
      yf_centroid_patch ( :,:,:, :) = yf_centroid_amr ( :,:,:, :, ii)
      zf_centroid_patch ( :,:,:, :) = zf_centroid_amr ( :,:,:, :, ii)

! --------------------

      xb_quadpts_patch ( :,:,:, :, :)= xb_quadpts_amr ( :,:,:, :, :, ii)
      yb_quadpts_patch ( :,:,:, :, :)= yb_quadpts_amr ( :,:,:, :, :, ii)
      zb_quadpts_patch ( :,:,:, :, :)= zb_quadpts_amr ( :,:,:, :, :, ii)

      xb_quadwt_patch ( :,:,:, :) = xb_quadwt_amr ( :,:,:, :, ii)
      yb_quadwt_patch ( :,:,:, :) = yb_quadwt_amr ( :,:,:, :, ii)
      zb_quadwt_patch ( :,:,:, :) = zb_quadwt_amr ( :,:,:, :, ii)

      xb_quadnorml_patch ( :,:,:, :) = xb_quadnorml_amr ( :,:,:, :, ii)
      yb_quadnorml_patch ( :,:,:, :) = yb_quadnorml_amr ( :,:,:, :, ii)
      zb_quadnorml_patch ( :,:,:, :, :) =
     1                              zb_quadnorml_amr ( :,:,:, :, :, ii)

! --------------------

      v_scale_length ( :, :, :) = v_scale_length_amr ( :, :, :, ii)

      v_quadpts_patch ( :,:,:, :, :) = v_quadpts_amr ( :,:,:, :, :, ii)

      v_quadwt_patch ( :, :, :, :) = v_quadwt_amr ( :, :, :, :, ii)

      v_volume_patch ( :, :, :) = v_volume_amr ( :, :, :, ii)

      v_centroid_patch ( :, :, :, :) = v_centroid_amr ( :, :, :, :, ii)

! --------------------

      Coef_Basis_Moments ( :, :, :, :) =
     1                         Coef_Basis_Moments_amr ( :, :, :, :, ii)

! --------------------------------------------------
! To be efficient on GPUs, we need to store in variables like "U_skinny" and then
! move them to GPU.

      U_skinny ( :, :, :, 1:n_cc_components) =
     1     U_amr ( :, :, :, 1:n_cc_components, ii)

      entropygr_skinny ( :, :, :) =
     1     entropygr_amr ( :, :, :, ii)

      Bxf_skinny ( :, :, :) = Bxf_amr ( :, :, :, ii)
      Byf_skinny ( :, :, :) = Byf_amr ( :, :, :, ii)
      Bzf_skinny ( :, :, :) = Bzf_amr ( :, :, :, ii)

      flattengr ( :, :, :) = flatten_amr ( :, :, :, ii)

! ----------------------------------------------------------------------
! extract "mu, eps" for irk-step 2. If electron pressure is also evolved,
! we have a rate at which the electron pressure catches up with the ion pressure.

      IF (  irk .EQ. nRK_stage) THEN

      CGL_fbeta_gr ( :, :, :) = CGL_fbeta_gr_amr ( :, :, :, ii)

      Pelectron_gr ( :, :, :) = Pelectron_gr_amr ( :, :, :, ii)

!$acc update device( CGL_fbeta_gr, Pelectron_gr )

      END IF

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"Loaded on-node patches; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! --------------------------------------------------
! Now that all the "_patch" and "_skinny" type data has been updated on CPU,
! move it to GPU.

!$acc update device ( U_skinny, entropygr_skinny, Bxf_skinny,Byf_skinny,&
!$acc&       Bzf_skinny, flattengr,                                     &

!$acc&       vert_coords_patch, imask_active_patch, istencil_patch,     &

!$acc&       isector, ipatch, ilevel, ix_patch, iy_patch, iz_patch,     &
!$acc&       itri_patch_LD, itri_patch_RD, itri_patch_RU,               &
!$acc&       itri_patch_LU,                                             &

!$acc&       xb_quadpts_patch, xb_quadwt_patch, xb_quadnorml_patch,     &
!$acc&       xf_quadpts_patch, xf_quadwt_patch, xf_quadnorml_patch,     &
!$acc&       xf_area_patch, xf_centroid_patch,                          &

!$acc&       yb_quadpts_patch, yb_quadwt_patch, yb_quadnorml_patch,     &
!$acc&       yf_quadpts_patch, yf_quadwt_patch, yf_quadnorml_patch,     &
!$acc&       yf_area_patch, yf_centroid_patch,                          &

!$acc&       zb_quadpts_patch, zb_quadwt_patch, zb_quadnorml_patch,     &
!$acc&       zf_quadpts_patch, zf_quadwt_patch, zf_quadnorml_patch,     &
!$acc&       zf_area_patch, zf_centroid_patch,                          &

!$acc&       v_scale_length, v_quadpts_patch, v_quadwt_patch,           &
!$acc&       v_volume_patch, v_centroid_patch,                          &

!$acc&       Coef_Basis_Moments)

! --------------------------------------------------

! Now we copy all the "_skinny" variables back to the actual variables on the GPU

!$acc kernels
      U ( :, :, :, 1:n_cc_components, 1) =
     1       U_skinny ( :, :, :, 1:n_cc_components)
      entropygr ( :, :, :, 1) = entropygr_skinny ( :, :, :)

      Bxf ( :, :, :, 1) = Bxf_skinny ( :, :, :)
      Byf ( :, :, :, 1) = Byf_skinny ( :, :, :)
      Bzf ( :, :, :, 1) = Bzf_skinny ( :, :, :)
!$acc end kernels













! ----------------------------------------------------------------------
! -STEP- 2) With geometry built, now do limiting, predictor step, flux
! and field evaluation.
! --------------------------------------------------














! Pad the boundaries in the radial direction if this is the innermost
! radial patch or the outermost radial patch.

! ------------------------------
! Most of the geomesh simulations require to implement the outflow  boundary
! conditions at the outer boundary. The PAD_BOUNDARY subroutine does that.

! "2" - outflow boundary condition. "17" - no boundary condition imposed.

! "3" - inflow/outflow boundary condition based on "inflow_uvec". In other words,
! when the flow is incoming, the boundary condition behaves like inflow; when
! the flow is outgoing, the boundary condition behaves like outflow.

      bcarr_r ( :, :, 1) = 17   ! Inner boundary.
      bcarr_r ( :, :, 2) = 17   ! Outer boundary.

      IF ( iz_patch .EQ. 1         )  bcarr_r ( :, :, 1) = 2     ! Inner boundary.
      IF ( iz_patch .EQ. nr_patches)  bcarr_r ( :, :, 2) = 2     ! Outer boundary.


#if( CHOOSE_APPL_CAK == 1)
      IF ( iz_patch .EQ. 1         )  bcarr_r ( :, :, 1) = 17    ! Inner boundary.
      IF ( iz_patch .EQ. nr_patches)  bcarr_r ( :, :, 2) = 2     ! Outer boundary.
#elif( CHOOSE_APPL_HELIO == 1)
      IF ( iz_patch .EQ. 1         )  bcarr_r ( :, :, 1) = 17    ! Inner boundary.
      IF ( iz_patch .EQ. nr_patches)  bcarr_r ( :, :, 2) = 17     ! Outer boundary.
#endif

!$acc update device ( bcarr_r)

      IF ( ( iz_patch .EQ. 1) .OR. ( iz_patch .EQ. nr_patches)) THEN

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

      END IF

! ------------------------------
! Most of the geomesh simulations require the re-initialization of boundary
! conditions at the inner boundary. The ADD_BOUNDARY subroutine does that.

      ibdy_inner = 0
      ibdy_outer = 0

      IF ( iz_patch .EQ. 1         )  ibdy_inner = 0   ! 1 == Impose Add_boundary
      IF ( iz_patch .EQ. nr_patches)  ibdy_outer = 0   ! 1 == Impose Add_boundary


#if( CHOOSE_APPL_CAK == 1)
      IF ( iz_patch .EQ. 1         )  ibdy_inner = 1   ! 1 == Impose Add_boundary
      IF ( iz_patch .EQ. nr_patches)  ibdy_outer = 0   ! 1 == Impose Add_boundary
#elif( CHOOSE_APPL_HELIO == 1)
      IF ( iz_patch .EQ. 1         )  ibdy_inner = 1   ! 1 == Impose Add_boundary
      IF ( iz_patch .EQ. nr_patches)  ibdy_outer = 1   ! 1 == Impose Add_boundary
#endif

!$acc update device ( ibdy_inner, ibdy_outer )

      IF ( ( iz_patch .EQ. 1) .OR. ( iz_patch .EQ. nr_patches)) THEN

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

      END IF

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"ADD_BOUNDARY; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! ----------------------------------------------------------------------
! Part of STEP 2) also as STEP 5) for even time steps.
! Evolve the internal energy for the cooling function source term;
! using the operator-splitting. This is a stiff source term. It is evolved
! using the adaptive-stepsize RK scheme.

! For odd-steps, call the operator split-source after the MHD evolution.

#if ( USE_COOLING_FUNCTION == 1)

      IF ( ( irk .EQ. nRK_stage) .AND. ( MOD ( istep, 2) == 1)) THEN

      CALL OPERATOR_SPLIT_SOURCE

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes, time, dtcur,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length, v_volume_patch, v_centroid_patch,

     1           U, entropygr, Bxf, Byf, Bzf, flag_Tfloor)

#if ( SILO_EXTRA_DATA >= 4)
!$acc update host( flag_Tfloor )
      flag_Tfloor_amr ( :, :, :, ii) = flag_Tfloor ( :, :, :)
#endif

      END IF

#endif

! --------------------------------------------------















! Apply limiter to MHD variables.

#if (NUMORDER == 2)
      CALL LIMIT_ALL_SLOPES_P1
#elif (NUMORDER == 3)
      CALL LIMIT_ALL_SLOPES_P12
#endif

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
     1         V, B1_zc, B2_zc, B3_zc, Ex_zc, Ey_zc, Ez_zc, Wzc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr, workgr,
     1         Bxf_quadpts, Byf_quadpts, Bzf_quadpts,

     1         V_min_local, V_max_local)

! Store local min and max for each patch that resides on this processor.
      V_min_amr ( :, ii) = V_min_local ( :)
      V_max_amr ( :, ii) = V_max_local ( :)

! For all the patches that reside on this processor, make a reduction of their
! min and max values.
      V_min_reduce = AMIN1 ( V_min_reduce, V_min_local)
      V_max_reduce = AMAX1 ( V_max_reduce, V_max_local)

! ----------
! divergence, stored in "workgr"

#if ( SILO_EXTRA_DATA >= 1)

!$acc update host( workgr )
      divgr_amr ( :, :, :, ii) = workgr ( :, :, :)

#endif

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"LIMIT_ALL_SLOPES_P1; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! --------------------------------------------------
! Call routines that build the fluxes and electric fields. Mostly, these are just
! flux evaluations.
! We also return the extremal speeds "S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr"
! from the 1D Riemann solvers for direct use in the 2D Riemann solvers.

      ibdy_outer = 0  ! Default is don't do anything at outer boundary.

#if ( CHOOSE_APPL_HELIO == 1)

! Heliospheric problem has very unique outer radial boundary. Set "ibdy_outer" properly.
! ! 1 == Impose outer boundary as outflow with RS.

      IF ( iz_patch .EQ. nr_patches)  ibdy_outer = 1

#endif

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

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"MAKE_FACIAL_FLUXES; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! --------------------------------------------------

! Call routine that builds the electric fields via 2D Riemann solver.
! The resstricted 2D Riemann solver is called at "xb_quadpts_patch, yb_quadpts_patch,
! zb_quadpts_patch". It usess speeds "S_L_gr, S_R_gr, S_D_gr, S_U_gr, S_B_gr, S_T_gr"
! from the 1D Riemann solvers. The zone-centered, and time-centered, reconstructions of
! "Ex_zc, Ey_zc, Ez_zc" are also used to get the central part of the electric field.
! The facial fields "Bxf_quadpts, Byf_quadpts, Bzf_quadpts" are suitably updated in
! time and their differences provide the dissipation.

#if ( USE_FACIAL_EFIELDS != 1)

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

#endif

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"MAKE_EDGE_EFIELDS_P1; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! ----------------------------------------------------------------------
! Build non-stiff source terms such as, gravitational force,
! rotational forces, etc.,

      CALL BUILD_NON_STIFF_SOURCES

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,

     1           iproblem, time, dtmove,

     1           vert_coords_patch, imask_active_patch,

     1           v_scale_length,
     1           v_volume_patch, v_centroid_patch,

     1           U, V, Source, dU_dt)

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"BUILD_NON_STIFF; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif















! ----------------------------------------------------------------------
! -STEP- 3) With geometry built, and with all fluxes in-hand, carry out
! the update. Note that these fluxes (or electric fields) have not been
! synchronized across patch boundaries. So the update is not perfect and we
! will have to make it perfect with a defect correction step.
! ----------------------------------------------------------------------
















      IF ( ( nRK_stage .EQ. 2) .AND. ( irk .EQ. 1) ) THEN

! ----------------------------------------
! Here the hyperbolic terms push the solution forward in time by a time of "dtcur".
! "UPDATE_IMEX_RK_STEP_P1" keeps track of "irk", i.e. the stages in the RK
! timestepping.

        dtmove = dtcur  ! Predictor step of RK-2.

!$acc update device ( dtmove)

!$acc kernels
        U0_amr ( :, :, :, 1:n_cc_components, ii) =
     1       U ( :, :, :, 1:n_cc_components, 1)

        entropygr0_amr ( :, :, :, ii) = entropygr ( :, :, :, 1)

        Bxf0_amr ( :, :, :, ii) = Bxf ( :, :, :, 1)
        Byf0_amr ( :, :, :, ii) = Byf ( :, :, :, 1)
        Bzf0_amr ( :, :, :, ii) = Bzf ( :, :, :, 1)
!$acc end kernels

! ----------------------------------------

      END IF ! End of "IF ( ( nRK_stage .EQ. 2) .AND. ( irk .EQ. 1) ) THEN".

! --------------------------------------------------
! --------------------------------------------------

      IF ( ( nRK_stage .EQ. 2) .AND. ( irk .EQ. 2) ) THEN

! ----------------------------------------
! Here the hyperbolic terms push the solution forward in time by a time of "dtcur/2".
! The solution is the time-average of the original solution with which we started
! the timestep, and the solution from the predictor step.
! "UPDATE_IMEX_RK_STEP_P1" keeps track of "irk", i.e. the stages in the RK
! timestepping, so that it knows how to combine the source terms.

        dtmove =  0.5 * dtcur  ! Corrector step of RK-2

!$acc update device ( dtmove)

!$acc kernels

! ----------
#if( IPDE_SYSTEM == 1)
! ----------

        U ( :, :, :, 1:NVAR, 1) =
     1       ( U ( :, :, :, 1:NVAR, 1)
     1  + U0_amr ( :, :, :, 1:NVAR, ii)) * 0.5

#if ( NFLUID >= 1)
        U ( :, :, :, 9: NVAR + NFLUID, 1) =
     1       ( U ( :, :, :, 9: NVAR + NFLUID, 1)
     1  + U0_amr ( :, :, :, 9: NVAR + NFLUID, ii)) * 0.5
#endif

! ----------
#elif( IPDE_SYSTEM == 2)
! ----------

        U ( :, :, :, 1: 8, 1) =
     1       ( U ( :, :, :, 1: 8, 1)
     1  + U0_amr ( :, :, :, 1: 8, ii)) * 0.5

        Un9 ( :, :, :) = U0_amr ( :, :, :, 9, ii)   ! 9-th component - CGL

#if ( NVAR >= 10)
        U ( :, :, :, 10: NVAR, 1) =
     1       ( U ( :, :, :, 10: NVAR, 1)
     1  + U0_amr ( :, :, :, 10: NVAR, ii)) * 0.5
#endif

! ----------
#elif( IPDE_SYSTEM == 3)
! ----------

        U ( :, :, :, 1: 8, 1) =
     1       ( U ( :, :, :, 1: 8, 1)
     1  + U0_amr ( :, :, :, 1: 8, ii)) * 0.5

        Un9 ( :, :, :) = U0_amr ( :, :, :, 9, ii)   ! 9-th component - CGL

        Un10 ( :, :, :) = U0_amr ( :, :, :, 10, ii)  ! 10-th component - Electron pressure

#if ( NVAR >= 11)
        U ( :, :, :, 11: NVAR, 1) =
     1       ( U ( :, :, :, 11: NVAR, 1)
     1  + U0_amr ( :, :, :, 11: NVAR, ii)) * 0.5
#endif

! ----------
#endif
! ----------

        entropygr ( :, :, :, 1) = ( entropygr ( :, :, :, 1)
     1                       + entropygr0_amr ( :, :, :, ii)) * 0.5

        Bxf ( :, :, :, 1) = ( Bxf ( :, :, :, 1)
     1                 + Bxf0_amr ( :, :, :, ii)) * 0.5

        Byf ( :, :, :, 1) = ( Byf ( :, :, :, 1)
     1                 + Byf0_amr ( :, :, :, ii)) * 0.5

        Bzf ( :, :, :, 1) = ( Bzf ( :, :, :, 1)
     1                 + Bzf0_amr ( :, :, :, ii)) * 0.5

!$acc end kernels

! ----------------------------------------

      END IF ! End of "IF ( ( nRK_stage .EQ. 2) .AND. ( irk .EQ. 2) ) THEN".

! --------------------------------------------------

      dt_defcor = dtmove        ! backup the dtcur value for defect-corrector step.
      time_defcor = time        ! backup the time  value for defect-corrector step.

!$acc update device ( time, dtcur, dt_defcor, time_defcor )

      CALL UPDATE_IMEX_RK_STEP_P1

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_components, n_cc_modes,
     1         n_fc_modes, time, dtmove,

     1         irk,

     1         vert_coords_patch, imask_active_patch,

     1         xf_quadnorml_patch, yf_quadnorml_patch,
     1         zf_quadnorml_patch,
     1         xf_area_patch, yf_area_patch, zf_area_patch,
     1         v_scale_length, v_volume_patch, v_centroid_patch,

     1         U, entropygr, Bxf, Byf, Bzf, dU_dt,
     1         CGL_fbeta_gr, Un9, Pelectron_gr, Un10,

     1         EntFlux_x, EntFlux_y, EntFlux_z,
     1         Efx, Efy, Efz, flag_entprs)

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"UPDATE_IMEX_RK; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! --------------------------------------------------
! After the above update, we want to make sure that all variables have been
! brought within physical bounds. We do that here.

      CALL CLEANUP_AFTER_UPDATE

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           iz_patch, nr_patches,

     1           vert_coords_patch, imask_active_patch,

     1           xf_quadnorml_patch, xf_centroid_patch,
     1           yf_quadnorml_patch, yf_centroid_patch,
     1           zf_quadnorml_patch, zf_centroid_patch,

     1           v_scale_length, v_centroid_patch,

     1      itri_patch_LD, itri_patch_RD, itri_patch_RU, itri_patch_LU,
     1      istencil_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr)

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"CLEANUP_AFTER_UPDATE; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      IF (  irk .EQ. nRK_stage) THEN















! ----------------------------------------------------------------------
! -STEP- 4) Evaluate a new timestep for this patch. Store "dtcurp1" so that we can
! globally reduce it after this loop is done on ALL processors.


















      CALL EVALUATE_TIMESTEP

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           time, dtcur, dtcurp1,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length,

     1           v_centroid_patch,

     1         U, entropygr, Bxf, Byf, Bzf,
     1         V, B1_zc, B2_zc, B3_zc,

     1         gamma_eos_gr, gamma_soundspeed_gr, mean_mol_wt_gr,
     1         divvelgr, msonicgr, flattengr)

! Store local "dtcurp1" for each patch that resides on this processor.
! Inside the timestep subroutine, "dtcurp1" has been moved back from GPU to CPU.
      dtcurp1_amr ( ii) = dtcurp1

! For all the patches that reside on this processor, make a reduction of their
! "dtcurp1" value.
      dtcurp1_reduce = AMIN1 ( dtcurp1_reduce, dtcurp1)

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"EVALUATE_TIMESTEP; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif












! ----------------------------------------------------------------------
! -STEP- 5) Evolve the internal energy for the cooling function source term;
! using the operator-splitting. This is a stiff source term. It is evolved
! using the adaptive-stepsize RK scheme.











! For even-steps, call the operator split-source after the MHD evolution.

#if ( USE_COOLING_FUNCTION == 1)

      IF ( MOD ( istep, 2) == 0) THEN

      CALL OPERATOR_SPLIT_SOURCE

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes, time, dtcur,

     1           vert_coords_patch, imask_active_patch,
     1           v_scale_length, v_volume_patch, v_centroid_patch,

     1           U, entropygr, Bxf, Byf, Bzf, flag_Tfloor)

#if ( SILO_EXTRA_DATA >= 4)
!$acc update host( flag_Tfloor )
      flag_Tfloor_amr ( :, :, :, ii) = flag_Tfloor ( :, :, :)
#endif

      END IF

#endif

! ----------------------------------------------------------------------

      END IF ! End of "IF (  irk .EQ. nRK_stage) THEN".

! ----------------------------------------------------------------------













! ----------------------------------------------------------------------
! -STEP- 6) "U, entropygr, Bxf, Byf, Bzf" have been almost perfectly updated on GPU.
! So, copy them back to CPU. The only fluxes we need are the ones that fill in
! the faces for the defect correction. Therefore, we have an extra subroutine
! on the GPU that copies "Flux_snd_r, l, u, d, t, b_patch" back to CPU.














! First, we load all the "_skinny" variables; Second, use the flux and
! electric field arrays on GPU to load the "Flux_snd_r_patch" etc. on the GPU.
! Only then can we efficiently bring them back to CPU.

!$acc kernels
      U_skinny ( :, :, :, 1:n_cc_components) =
     1       U ( :, :, :, 1:n_cc_components, 1)

      entropygr_skinny ( :, :, :) = entropygr ( :, :, :, 1)

      Bxf_skinny ( :, :, :) = Bxf ( :, :, :, 1)
      Byf_skinny ( :, :, :) = Byf ( :, :, :, 1)
      Bzf_skinny ( :, :, :) = Bzf ( :, :, :, 1)
!$acc end kernels


      CALL FILL_FLUX_SND_PATCH

     1         ( nr_zones, ntp_zones, NGhost,
     1           n_cc_components, n_cc_modes,
     1           n_fc_modes,

     1           EntFlux_x, EntFlux_y, EntFlux_z,
     1           Efx, Efy, Efz,

     1           Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1           Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi,

     1           Flux_snd_r_patch, Flux_snd_l_patch, Flux_snd_u_patch,
     1           Flux_snd_d_patch, Flux_snd_t_patch, Flux_snd_b_patch)

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"FILL_FLUX_SND; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! ----------------------------------------------------------------------
!$acc update host ( U_skinny, entropygr_skinny,                         &
!$acc&              Bxf_skinny, Byf_skinny, Bzf_skinny,                 &
!$acc&      Flux_snd_r_patch, Flux_snd_l_patch, Flux_snd_u_patch,       &
!$acc&      Flux_snd_d_patch, Flux_snd_t_patch, Flux_snd_b_patch,       &
!$acc&      CGL_fbeta_gr, Pelectron_gr)

      U_amr ( :, :, :, 1:n_cc_components, ii) =
     1           U_skinny ( :, :, :, 1:n_cc_components)

      entropygr_amr ( :, :, :, ii) = entropygr_skinny ( :, :, :)

      Bxf_amr ( :, :, :, ii) = Bxf_skinny ( :, :, :)
      Byf_amr ( :, :, :, ii) = Byf_skinny ( :, :, :)
      Bzf_amr ( :, :, :, ii) = Bzf_skinny ( :, :, :)

      Flux_snd_r_amr ( :, :, :, ii) = Flux_snd_r_patch ( :, :, :)
      Flux_snd_l_amr ( :, :, :, ii) = Flux_snd_l_patch ( :, :, :)

      Flux_snd_u_amr ( :, :, :, ii) = Flux_snd_u_patch ( :, :, :)
      Flux_snd_d_amr ( :, :, :, ii) = Flux_snd_d_patch ( :, :, :)

      Flux_snd_t_amr ( :, :, :, ii) = Flux_snd_t_patch ( :, :, :)
      Flux_snd_b_amr ( :, :, :, ii) = Flux_snd_b_patch ( :, :, :)


! ----------------------------------------------------------------------
! Store "CGL_fbeta_gr" for irk-step 2, no need to be communicated across patches/procs.
! Do the same for "Pelectron_gr" which is only needed if we have CGL with electron pressure.

      IF (  irk .EQ. 1) THEN

!$acc update host ( CGL_fbeta_gr)
      CGL_fbeta_gr_amr ( :, :, :, ii) = CGL_fbeta_gr ( :, :, :)

!$acc update host ( Pelectron_gr)
      Pelectron_gr_amr ( :, :, :, ii) = Pelectron_gr ( :, :, :)

      END IF

! ----------------------------------------------------------------------
! We may never need this, but it helps to have the lines to fill in the fluxes.
! If bandwidth to main memory ever improves for GPUs, this would be an excellent step.

#if ( 0 == 1)

!$acc update host ( EntFlux_x, EntFlux_y,                               &
!$acc&              EntFlux_z, Efx, Efy, Efz)

      EntFlux_x_amr ( :, :, :, ii) = EntFlux_x ( :, :, :)
      EntFlux_y_amr ( :, :, :, ii) = EntFlux_y ( :, :, :)
      EntFlux_z_amr ( :, :, :, ii) = EntFlux_z ( :, :, :)

      Efx_amr ( :, :, :, ii) = Efx ( :, :, :)
      Efy_amr ( :, :, :, ii) = Efy ( :, :, :)
      Efz_amr ( :, :, :, ii) = Efz ( :, :, :)

#endif

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"Copy back; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! Let all patches synchronize before variable min, max reduction.
! That way, they have all updated their "V_min_reduce, V_max_reduce".

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! END LOOP OVER PATCHES ON PROCESSOR   ;   END LOOP OVER PATCHES ON PROCESSOR
! ----------------------------------------------------------------------









! -STEP- 7) Now make a reduction across all processors for "V_min_reduce, V_max_reduce".









      V_min_local = V_min_reduce
      V_max_local = V_max_reduce

      CALL MPI_ALLREDUCE ( V_min_local,          ! Input for reudction
     1                     V_min_reduce,         ! Global minimum
     1                     n_cc_components + 1,  ! length of reduction variable
     1                     MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

      CALL MPI_ALLREDUCE ( V_max_local,          ! Input for reudction
     1                     V_max_reduce,         ! Global maximum
     1                     n_cc_components + 1,  ! length of reduction variable
     1                     MPI_REAL8, MPI_MAX, MPI_COMM_WORLD, ierr)

! Write out the global min and max of all variables from processor 0.

      IF ( me .EQ. 1) THEN

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Step # = ", istep," time = ", time,
     1             " dtcur = ", dtcur
      WRITE ( 6, *)" "
      WRITE ( 6, *)"V_min = "
      WRITE ( 6, 27) V_min_reduce ( :)
      WRITE ( 6, *)" "
      WRITE ( 6, *)"V_max = "
      WRITE ( 6, 27) V_max_reduce ( :)
  27  FORMAT ( 5(1pe13.6,2x))
      WRITE ( 6, *)" "
      CALL FLUSH ( 6)

      END IF

! Let all patches synchronize.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------














! ----------------------------------------------------------------------
! -STEP- 8) Now recall that all the fluxes are slightly defective. To make
! them perfect, we have to synchronize flux and e-field data across patch boundaries.














! Synchronize all flux & e-field data for
! "EntFlux_x_amr, EntFlux_y_amr, EntFlux_z_amr,
! Efx_amr, Efy_amr, Efz_amr".
! All these subroutines have "CALL MPI_Barrier" at begin and end.

! Once we come here, all defects should start off as zeros, even if they are
! filled in at appropriate locations in the course of this step.

      DO ii = 1, npatch_per_proc

        Defect_EntFlux_x_amr ( :, :, :, ii) = 0.0
        Defect_EntFlux_y_amr ( :, :, :, ii) = 0.0
        Defect_EntFlux_z_amr ( :, :, :, ii) = 0.0

        Defect_Efx_amr ( :, :, :, ii) = 0.0
        Defect_Efy_amr ( :, :, :, ii) = 0.0
        Defect_Efz_amr ( :, :, :, ii) = 0.0

! Fluctuation terms.

        Defect_Fluct_x_Pl_amr ( :, :, :, :, ii) = 0.0
        Defect_Fluct_y_Pl_amr ( :, :, :, :, ii) = 0.0
        Defect_Fluct_z_Pl_amr ( :, :, :, :, ii) = 0.0

        Defect_Fluct_x_Mi_amr ( :, :, :, :, ii) = 0.0
        Defect_Fluct_y_Mi_amr ( :, :, :, :, ii) = 0.0
        Defect_Fluct_z_Mi_amr ( :, :, :, :, ii) = 0.0

      END DO

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! --------------------------------------------------
! -STEP- 8.a) Fill the boundary values of "Fluctuations_xyz_amr"

      CALL LOAD_AMR_FLUXES_AT_BDY ( )

! --------------------------------------------------
! -STEP- 8.b) Synchronize all flux & e-field data

      CALL LOAD_LOCAL_BUFFERS_Flux_XY ()

      CALL FETCH_BUFFERS_FILL_IN_Flux_X ()
      CALL FETCH_BUFFERS_FILL_IN_Flux_Y ()

      CALL LOAD_LOCAL_BUFFERS_Flux_Z ()
      CALL FETCH_BUFFERS_FILL_IN_Flux_Z ()

      CALL SYNCHRONIZE_EFZ_QUAD_CORNER ()

! --------------------------------------------------
! -STEP- 8.c) Unpack the defect fluxes and defect electric fields.

      CALL UNPACK_DEFECT_FLUXES ( )

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"UNPACK_DEFECT_FLUXES; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif










! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! -STEP- 9) The defects have been built in the previous step.
! This is the defect correction step. It operates just like a timestep.
! The defects have to be stored as if they are fluxes, which has already
! been done in the previous step. The defects only affect the boundary zones.
! This is such a low cost operation that it can be done entirely on CPU.















! ----------------------------------------------------------------------
! BEGIN LOOP OVER PATCHES ON PROCESSOR   ;   BEGIN LOOP OVER PATCHES ON PROCESSOR
! ----------------------------------------------------------------------
! Now we loop over all "npatch_per_proc" patches on each processor and
! just correct the defects.

      DO ii = 1, npatch_per_proc

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      isector = iactive_FromProc ( ii, 1)
      ipatch = iactive_FromProc ( ii, 2)

      ilevel = isector_ilevel_amr ( 2, ii)

      ix_patch = ixyz_patch_amr ( 1, ii)
      iy_patch = ixyz_patch_amr ( 2, ii)
      iz_patch = ixyz_patch_amr ( 3, ii)

! --------------------------------------------------

      vert_coords_patch ( :, :, :, :) =
     1                          vert_coords_amr ( :, :, :, :, ii)

      imask_active_patch ( :, :, :) =
     1                          imask_active_amr ( :, :, :, ii)

! --------------------------------------------------

      itri_patch_LD = itri_patch_amr ( 1, ii)
      itri_patch_RD = itri_patch_amr ( 2, ii)
      itri_patch_RU = itri_patch_amr ( 3, ii)
      itri_patch_LU = itri_patch_amr ( 4, ii)

      istencil_patch ( :, :, :, :) = istencil_amr ( :, :, :, :, ii)

! --------------------------------------------------
! This is entirely on CPU, so we can copy from "U_amr" to "U" and so on.
! Also copy in the minimal amount of geometry variables, because we will
! not be truly rebuilding geometry here.

      U ( :, :, :, 1:n_cc_components, 1) =
     1     U_amr ( :, :, :, 1:n_cc_components, ii)

      entropygr ( :, :, :, 1) =
     1     entropygr_amr ( :, :, :, ii)

      Bxf ( :, :, :, 1) = Bxf_amr ( :, :, :, ii)
      Byf ( :, :, :, 1) = Byf_amr ( :, :, :, ii)
      Bzf ( :, :, :, 1) = Bzf_amr ( :, :, :, ii)

      Source ( :, :, :, :) = 0.0     ! NO Source for defect-correction step.

      xf_quadnorml_patch ( :, :, :, :) =
     1                              xf_quadnorml_amr ( :, :, :, :, ii)

      yf_quadnorml_patch ( :, :, :, :) =
     1                              yf_quadnorml_amr ( :, :, :, :, ii)

      zf_quadnorml_patch ( :, :, :, :, :) =
     1                           zf_quadnorml_amr ( :, :, :, :, :, ii)


      xf_area_patch ( :, :, :) = xf_area_amr ( :, :, :, ii)
      yf_area_patch ( :, :, :) = yf_area_amr ( :, :, :, ii)
      zf_area_patch ( :, :, :) = zf_area_amr ( :, :, :, ii)


      v_scale_length ( :, :, :) = v_scale_length_amr ( :, :, :, ii)

      v_volume_patch ( :, :, :) = v_volume_amr ( :, :, :, ii)

      v_centroid_patch ( :, :, :, :) = v_centroid_amr ( :, :, :, :, ii)

! --------------------------------------------------
! The defects are used in the update step just like fluxes. Therefore,
! copy in "Defect_Fluct_x_Pl_amr" etc. into "Fluct_x_Pl" etc. Likewise, copy in
! "Defect_Efx_amr" etc. into "Efx" etc.

      EntFlux_x ( :, :, :) = Defect_EntFlux_x_amr ( :, :, :, ii)
      EntFlux_y ( :, :, :) = Defect_EntFlux_y_amr ( :, :, :, ii)
      EntFlux_z ( :, :, :) = Defect_EntFlux_z_amr ( :, :, :, ii)

      Efx ( :, :, :) = Defect_Efx_amr ( :, :, :, ii)
      Efy ( :, :, :) = Defect_Efy_amr ( :, :, :, ii)
      Efz ( :, :, :) = Defect_Efz_amr ( :, :, :, ii)

      Fluct_x_Pl ( :, :, :, :) = Defect_Fluct_x_Pl_amr ( :, :, :, :, ii)
      Fluct_y_Pl ( :, :, :, :) = Defect_Fluct_y_Pl_amr ( :, :, :, :, ii)
      Fluct_z_Pl ( :, :, :, :) = Defect_Fluct_z_Pl_amr ( :, :, :, :, ii)

      Fluct_x_Mi ( :, :, :, :) = Defect_Fluct_x_Mi_amr ( :, :, :, :, ii)
      Fluct_y_Mi ( :, :, :, :) = Defect_Fluct_y_Mi_amr ( :, :, :, :, ii)
      Fluct_z_Mi ( :, :, :, :) = Defect_Fluct_z_Mi_amr ( :, :, :, :, ii)

! --------------------------------------------------
! Defect correction step. Here, use the previously stored "time_defcor, dt_defcor".


      CALL UPDATE_DEFECT_IMEX_RK

     1       ( nr_zones, ntp_zones, NGhost,
     1         n_cc_components, n_cc_modes,
     1         n_fc_modes, time_defcor, dt_defcor,  ! note _defcor.

     1         irk,

     1         vert_coords_patch, imask_active_patch,

     1         xf_quadnorml_patch, yf_quadnorml_patch,
     1         zf_quadnorml_patch,
     1         xf_area_patch, yf_area_patch, zf_area_patch,
     1         v_scale_length, v_volume_patch, v_centroid_patch,

     1         U, entropygr, Bxf, Byf, Bzf, dU_dt,

     1         Fluct_x_Pl, Fluct_y_Pl, Fluct_z_Pl,
     1         Fluct_x_Mi, Fluct_y_Mi, Fluct_z_Mi,

     1         EntFlux_x, EntFlux_y, EntFlux_z,
     1         Efx, Efy, Efz, flag_entprs)


#if ( SILO_EXTRA_DATA >= 3)
      flag_entprs_amr ( :, :, :, ii) = flag_entprs ( :, :, :)
#endif

#if (DEBUG == 1)
      IF ( me == 1) THEN
      WRITE (6, *)"CALL UPDATE_DEFECT_IMEX; istep, irk, time, dtcur = ",
     1             istep, irk, time, dtcur
      CALL FLUSH ( 6)
      END IF
#endif

! --------------------------------------------------
! Now reverse the process; i.e. copy "U" into "U_amr" etc.

      U_amr ( :, :, :, 1:n_cc_components, ii) =
     1      U ( :, :, :, 1:n_cc_components, 1)

      entropygr_amr ( :, :, :, ii) = entropygr ( :, :, :, 1)

      Bxf_amr ( :, :, :, ii) = Bxf ( :, :, :, 1)
      Byf_amr ( :, :, :, ii) = Byf ( :, :, :, 1)
      Bzf_amr ( :, :, :, ii) = Bzf ( :, :, :, 1)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO ii = 1, npatch_per_proc".

! Let all patches synchronize before starting the next timestep.

      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO irk = 1, nRK_stage".

! ----------------------------------------------------------------------
! END LOOP OVER PATCHES ON PROCESSOR   ;   END LOOP OVER PATCHES ON PROCESSOR
! ----------------------------------------------------------------------
















! ----------------------------------------------------------------------
! -STEP- 10) Obtain globally minimized "dtcurp1_reduce" and broadcast it to
! all "dtcurp1". Also update "time".















! Let all patches synchronize before timestep reduction.
! That way, they have all updated their "dtcurp1_reduce".
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! Now make a reduction across all processors for "dtcurp1_reduce".

      CALL MPI_ALLREDUCE ( dtcurp1_reduce,   ! Input for reudction
     1                     dtcurp1,          ! Global minimum
     1                     1,                ! length of reduction variable
     1                     MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

      CALL CPU_TIME ( timing_timestep2)

      time = time + dtcur

      IF ( me .EQ. 1) THEN

      WRITE ( 6, *)" "
      WRITE ( 6, *)"Step # = ", istep," found dtcurp1 = ", dtcurp1
      WRITE ( 6, *)" "
      WRITE ( 6, *)"Timed this timestep = ", timing_timestep2
     1                                     - timing_timestep1
      WRITE ( 6, *)" "

      END IF

! ----------------------------------------

#if ( ACCURACY_ANALYSIS == 1)
      CALL MAKE_GLOBAL_ACCURACY_CS ( iproblem, rmin, rmax, time, dtcur)
#endif

! ----------------------------------------

      IF ( me .EQ. 1) THEN

      WRITE ( 6, *)"-------------------------------------------------"
      WRITE ( 6, *)" "

      END IF

! ----------------------------------------------------------------------
! Make the code to stop at the "time == timestop". This is only needed for
! accuracy analysis, so it is commented out in production code.
!      IF ( (time + dtcurp1) .GT. timestop)
!     1   dtcurp1 = MAX ( ( timestop - time), smallnum / 100.0)

      dtcur = dtcurp1

!$acc update device ( time, dtcur, dtcurp1)

! ----------------------------------------------------------------------

! Let all patches synchronize.
      CALL MPI_Barrier ( MPI_COMM_WORLD, ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END DO  ! End of "DO istep = 1, ntstep".

! ----------------------------------------------------------------------
! END TIMESTEP LOOP  ;  END TIMESTEP LOOP  ;  END TIMESTEP LOOP
! ----------------------------------------------------------------------

















! ----------------------------------------------------------------------
! We created a lot of temporary arrays that are used for processing a timestep
! on the GPU. Since we are at the end of the timestep, we can get rid of
! these arrays on the GPU.

12345  CONTINUE

!$acc end data
! ---------------------------------------------------------------------
! All the messaging is done; about to stop code. Now free up all the handles.

      CALL WINFREE_MESSAGING ()

      CALL MPI_FINALIZE (ierr)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      END PROGRAM RIEMANN_GEOMESH_CS_MHD




