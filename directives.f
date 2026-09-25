
! ~~~main parameters~~~ !


! Choose the PDE System. (1) MHD; (2) CGL; (3) CGL with p_elec and 3 tracers
! "IPDE_SYSTEM == 1" for MHD
! "IPDE_SYSTEM == 2" for CGL System
! "IPDE_SYSTEM == 3" for CGL System with p_elec
#define IPDE_SYSTEM 2

! Define the order of accuracy; 2, 3 or 4.
#define NUMORDER 3



! ~~~BuildGeomesh_c3_src~~~ !


! Choice of Matmul loop (see Riemann_Integrate_CS)
#define DO_LOOP_FOR_MATMUL 0
! Whether to build the mesh in debug mode
#define DEBUG_GEOM 0



! ~~~application parameters~~~ !


! Choose application; so far only CAK wind problem and HELIO is added here; more application
! flags can be added, with the prefix "CHOOSE_APPL_".
! if the value(s) are zero, default iproblem runs from "Riemann_Module_Main_CS_src".
#define CHOOSE_APPL_CAK 0
#define CHOOSE_APPL_HELIO 1

! Only for heliosphere, we can choose monopole or split monopole. This also flips the
! Parker spiral in the hemispheres if split monopole is chosen.
#define SPLIT_MONOPOLE_HELIO 0

! Define the number of quadrature points at each face.
#define NAREA_QUAD_PTS 5

! Define the number of magnetic field integration points for the multid-rs
! At all orders we want 5 facial quadrature points.
#define NBFLD_QUAD_PTS 5

! Define the number of volumetric quadrature points that are needed for quick 4th
! order evaluation of source terms. (This is not as good as tensor product; but it is
! v.v. fast and sufficiently formally accurate.)
#define NQUICK_VOL_QUAD_PTS 6

! Choose whether logarithmic grid (1) or uniform grid (0) or geometrically
! ratioed mesh (2) or hybrid of geometric and log mesh (3) along the
! radial direction. There is also a mesh option (4) which can be anything between
! uniform and logarithmic. Suitable for Heliospheric problems. There is also
! a mesh option (5) which has uniform mesh inside and logarithmic mesh outside.
! (For hybrid, make further settings before call to "GENERATE_GLOBAL_RADIAL_HYBRID_MESH".)

! Florinski has his own favorite for heliospheric problems, which is (6,7,8,9)
! (6) for Uniform; (7) for Exponential; (8) for Power Law; (9) for Hybrid of uniform
! and exponential.

#define CHOOSE_RADIAL_MESH 9

! Choose whether we want adiabatic or isothermal versions of subroutines. They
! have slightly different subroutine interfaces. Isothermal case should be run
! without pressure positivity switch.
#define CHOOSE_ADIABATIC 1

! Just keep "NVAR 8" for basic MHD, "9" for CGL and "10+3" for CGL with p_e
#if( IPDE_SYSTEM == 1)
#define NVAR 8
#elif( IPDE_SYSTEM == 2)
#define NVAR 11
#elif( IPDE_SYSTEM == 3)
#define NVAR 12
#endif

! "NFLUID" is the number of MHD fluid species.
#define NFLUID 0
#define NFLUID_EULER 0

! We can optionally turn the 2D HLL RS into a 2D LLF RS by setting this ifdef
! to 1. This may increase the dissipation in the magnetic field, but it also
! stabilizes it. For anything but second order, please keep "USE_2D_LLF_RS_FOR_EZ == 0".
#define USE_2D_LLF_RS_FOR_EZ 1

! At second order, and only for second order, we can get by with only the zone-
! centered reconstruction of magnetic field variables. This is less expensive
! and works. Set to 0 for inexpensive version; 1 for expensive but better version.
#define FULL_BFLD_RECONSTRUCTION 0

! Sometimes, the magnetic field needs to be split into a strong constant part (B0)
! and a fluctuating part (B1) which can vary with time. The default is that this is not
! turned on so that the code defaults to (B0 = 0) all over.
! A function "GET_MAGNETIC_FIELD_B0B1" has to be provided with correct interface
! when "USE_B0B1_FORMULATION 1".
#define USE_B0B1_FORMULATION 0

! Usage of cooling function in the code. If it is defined "1", the cooling table
! is read from the data-file "cooltable.dat" and it is incorporated in the 
! operator split fashion.
#if( (CHOOSE_ADIABATIC == 1) && (CHOOSE_APPL_CAK == 1) )

#define USE_COOLING_FUNCTION 1
#define USE_INVERSE_COMPTON_COOLING 1 

#else

#define USE_COOLING_FUNCTION 0
#define USE_INVERSE_COMPTON_COOLING 0 

#endif

! "ACCURACY_ANALYSIS = 0 or 1" depending on whether we have an analytical function
! for the solution against which we are checking for accuracy.
#define ACCURACY_ANALYSIS 0

! For classical MHD we can use just the contact and Alfven waves with
! "COMPLETE_EIGENSYSTEM = 0". In the hydro limit, we replace Alfven with shear.
! To include Alfven and entropy waves, set "COMPLETE_EIGENSYSTEM = 1". Along
! with Alfven waves, the entropy waves carry shear in the flow.
! To include fast magnetosonic waves, set "COMPLETE_EIGENSYSTEM = 2".
! These carry the largest pressure jumps from shocks; but not needed in a Riemann solver.
#define COMPLETE_EIGENSYSTEM 0

! At second order, we can use facial fluxes to obtain edge-averaged electric fields.
! Best to keep zero and use the multid-RS for the edge-centered electric fields.
! Never set to 1, except for debugging.
#define USE_FACIAL_EFIELDS 0

! If RESTART is 0 we initialize the problem & then run. Otherwise, we
! assume a restart file has been specified in "riemann.in".
#define RESTART 0

! For problems with very strong shocks and large magnetic fields, one
! may have a loss of pressure positivity. The entropy formulation cures
! that and it is invoked by switching on the pressure positivity switch.
! It does cause a local loss of conservation in the zone where it is invoked.
! It is switched on by setting "PRESSURE_POSITIVITY" to 1 or 2 (0, default).
! When "PRESSURE_POSITIVITY == 1" we only check for positivity of density or
! pressure. When "PRESSURE_POSITIVITY == 2" we check all the variables.
#define PRESSURE_POSITIVITY 1

! If WRITE_DIAGNOSTICS is set to 1, the code writes out min and max of all
! variables at each timestep. Can be used with OpenMP, but not with MPI
! versions of code. Setting it to 0 turns off the writing of diagnostics. The
! diagnostics are evaluated and written out in subroutine "LIMIT_ALL_SLOPES".
#define WRITE_DIAGNOSTICS 1

! Sets the character length of the variable names in the silo file.
#define SILO_VARLEN 24

! Set GPU execution to "0" for CPU. when it is "1", it disables the error
! messages getting written out; inside the "acc routine seq" subroutines.  
! Without this flag the GPU code does not successfully compile.
#define GPU_EXECUTION 0

! Number of extra variables, written out in the silo file. This helps in
! monitoring the simulation. "1" write out divergence and other extra variables
! as needed (see below). "0" write out only the conserved variables.
#define SILO_EXTRA_DATA 1


! Write out more extra data as needed.
#if ( ( SILO_EXTRA_DATA >= 1) && ( PRESSURE_POSITIVITY >= 1))

#if ( USE_COOLING_FUNCTION == 1)

! Write out 1) divergence  2) entropy  3) entropy-use-flag  4) tempfloor-use-flag.
#undef SILO_EXTRA_DATA
#define SILO_EXTRA_DATA 4

#else

! Write out 1) divergence  2) entropy  3) entropy-use-flag
#undef SILO_EXTRA_DATA
#define SILO_EXTRA_DATA 3

#endif

#endif



! CAK wind parameters needed in Implicit Taylor Predictor and nowhere else.
#define USE_CAK_WIND 0


