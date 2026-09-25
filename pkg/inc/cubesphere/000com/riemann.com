! ----------------------------------------------------------------------
! ----------------------------------------------------------------------

      LOGICAL do_on_gpu

      INTEGER igeom, ientropyfix, iflag_monopole, i_dummy, 
     1        iflag_npoles ( 4)

      REAL pi, gamma, molewt, gasconst, rhofloor, prsfloor, smallnum,
     1     gamma_euler, molewt_euler, gasconst_euler,
     1     isospeed, barotropic_stiffness, lightspd, stefan,
     1     newtonsc, Mstar, omega_passive ( 3), omega_active ( 3),
     1     radial_mesh_ratio,
     1     radial_mesh_transition, radial_mesh_fraction,
     1     tempfloor,
     1     multipole_Bmag ( 4), inflow_uvec ( 3),
     1     tau_phys_cgl, f_lowbd_cgl, tau_phys_sw, tau_phys_lism,
     1     neta_f_cgl, kappa_f_cgl, magB_small_cgl, eps_bdry_cgl,
     1     source_tau, coef_Pelectron,

     1     cfl_coef,
     1     mc_coef, flatten_coef, prs_flatten_coef,
     1     fracdiff, divratio_einfeldt_rs,
     1     lapidus_threshold, lapidus_coef,
     1     blend_multid_flux, blend_mdhll_to_mdhllc,

     1     base_grid_rmin, base_grid_rmax,

     1     gamma_specie ( 0 : NFLUID), cp_specie ( 0 : NFLUID),
     1     cv_specie ( 0 : NFLUID), molewt_specie ( 0 : NFLUID),
     1     heat_of_formation ( 0 : NFLUID),

     1     gamma_specie_euler ( 0 : NFLUID_EULER),
     1     cp_specie_euler ( 0 : NFLUID_EULER),
     1     cv_specie_euler ( 0 : NFLUID_EULER),
     1     molewt_specie_euler ( 0 : NFLUID_EULER),
     1     heat_of_formation_euler ( 0 : NFLUID_EULER),

     1     monopole_u0, monopole_B ( 3), 
     1     dipole_B ( 3), qpole_B ( 3, 3)

      COMMON/RIEMANN_COM/ igeom, ientropyfix, iflag_monopole, i_dummy,
     1     iflag_npoles,
     1     pi, gamma, molewt, gasconst, rhofloor, prsfloor, smallnum,
     1     gamma_euler, molewt_euler, gasconst_euler,
     1     isospeed, barotropic_stiffness, lightspd, stefan,
     1     newtonsc, Mstar, omega_passive, omega_active,
     1     radial_mesh_ratio,
     1     radial_mesh_transition, radial_mesh_fraction,
     1     tempfloor,
     1     multipole_Bmag, inflow_uvec,
     1     tau_phys_cgl, f_lowbd_cgl, tau_phys_sw, tau_phys_lism,
     1     neta_f_cgl, kappa_f_cgl, magB_small_cgl, eps_bdry_cgl,
     1     source_tau, coef_Pelectron,

     1     cfl_coef,
     1     mc_coef, flatten_coef, prs_flatten_coef,
     1     fracdiff, divratio_einfeldt_rs,
     1     lapidus_threshold, lapidus_coef,
     1     blend_multid_flux, blend_mdhll_to_mdhllc,

     1     base_grid_rmin, base_grid_rmax,

     1     gamma_specie, cp_specie,
     1     cv_specie, molewt_specie,
     1     heat_of_formation,

     1     gamma_specie_euler, cp_specie_euler,
     1     cv_specie_euler, molewt_specie_euler,
     1     heat_of_formation_euler,

     1     monopole_u0, monopole_B, dipole_B, qpole_B,

     1     do_on_gpu

!$acc declare create (/RIEMANN_COM/)

! ----------------------------------------------------------------------
! ----------------------------------------------------------------------






