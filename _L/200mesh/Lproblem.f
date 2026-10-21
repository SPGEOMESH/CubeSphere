

#include "directives.f"
!include <cubesphere/000com/090init_riemann.f>

#include <cubesphere/150util/160print_vect.f>

#include <cubesphere/660pde_3d/660pde_3d_force_pperp_ppar_to_hyperbolic_domain.f>
#include <cubesphere/660pde_3d/660pde_3d_prim_to_wavspd.f>
#include <cubesphere/660pde_3d/661pde_3d_prim_to_cons.f>
#include <cubesphere/660pde_3d/662pde_3d_cons_to_prim.f>
#include <cubesphere/700cak_helio/helio/Heliosphere_src.f>
#include <cubesphere/800applications_cs/801problem_generator_u_a.f>
#include <cubesphere/800applications_cs/820init_patch_cube_sphere.f>


